/*
Bullet Continuous Collision Detection and Physics Library
https://bulletphysics.org

This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from the use of this software.
Permission is granted to anyone to use this software for any purpose,
including commercial applications, and to alter it and redistribute it freely,
subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*/

#include "btVariationalProjectiveDynamics.h"
#include <stdio.h>

void btVariationalProjectiveDynamics::localStep(const btVariationalSystem& system, const btVec3d* x, btVec3d* rhs) const
{
	const int n = system.numNodes();
	const double invH2 = 1.0 / (system.m_dt * system.m_dt);
	// inertia: (M/h^2) xTilde
	for (int i = 0; i < n; i++) rhs[i] = system.m_xTilde[i] * (system.m_mass[i] * invH2);
	// tetrahedra: k V G^T p with p the projection of F
	for (int t = 0; t < system.m_tets.size(); t++)
	{
		const btVariationalTet& tet = system.m_tets[t];
		const btVariationalMaterial* material = system.m_bodies[tet.m_body].m_material;
		const btMat3d F = system.deformationGradient(tet, x);
		btMat3d target;
		double k = 0;
		if (material->projectiveDynamicsTarget(F, target, k))
			target *= k * tet.m_restVolume;
		else
			// not a projective energy: (k F - dPsi/dF) reproduces a step along the energy gradient
			target = (F * material->laplacianStiffness() - material->firstPiola(F)) * tet.m_restVolume;
		for (int j = 0; j < 4; j++) rhs[tet.m_n[j]] += target * tet.m_g[j];
	}
	// springs: the projection is the current direction at rest length
	for (int s = 0; s < system.m_springs.size(); s++)
	{
		const btVariationalSpring& sp = system.m_springs[s];
		const btVec3d d = x[sp.m_n[0]] - x[sp.m_n[1]];
		const double len = d.length();
		if (len < 1e-300) continue;
		const btVec3d p = d * (sp.m_stiffness * sp.m_restLength / len);
		rhs[sp.m_n[0]] += p;
		rhs[sp.m_n[1]] -= p;
	}
	for (int a = 0; a < system.m_attachments.size(); a++)
	{
		const btVariationalAttachment& at = system.m_attachments[a];
		const double k = at.m_hard ? at.m_penalty : at.m_stiffness;
		if (at.m_nodeB < 0)
			rhs[at.m_nodeA] += (at.m_target + at.m_c0) * k - at.m_lambda;
		else
		{
			const btVec3d p = at.m_c0 * k - at.m_lambda;
			rhs[at.m_nodeA] += p;
			rhs[at.m_nodeB] -= p;
		}
	}
}

void btVariationalProjectiveDynamics::projectColliders(const btVariationalSystem& system, btVec3d* x)
{
	const double margin = system.m_colliderMargin;
	for (int b = 0; b < system.m_bodies.size(); b++)
	{
		const btVariationalBody& body = system.m_bodies[b];
		for (int s = 0; s < body.m_colliders.size(); s++)
		{
			const btVariationalCollider& col = system.m_colliders[body.m_colliders[s]];
			const double mu = body.m_friction * col.m_friction;
			for (int i = body.m_firstNode; i < body.m_firstNode + body.m_numNodes; i++)
			{
				if (system.m_fixed[i]) continue;
				double d;
				btVec3d nrm;
				col.signedDistance(x[i], d, nrm);
				if (d >= margin) continue;
				const double depth = margin - d;
				x[i] += nrm * depth;
				if (!(mu > 0.0)) continue;
				// tangential motion relative to the collider during this step
				const btVec3d rel = (x[i] - system.m_xPrev[i]) - col.velocityAt(system.m_xPrev[i]) * system.m_dt;
				const btVec3d tangential = rel - nrm * rel.dot(nrm);
				const double len = tangential.length();
				if (len <= mu * depth)
					x[i] -= tangential;  // static friction holds the node
				else
					x[i] -= tangential * (mu * depth / len);
			}
		}
	}
}

void btVariationalProjectiveDynamics::solve(btVariationalSystem& system)
{
	const int n = system.numNodes();
	m_stats.reset();
	if (n == 0) return;
	m_rhs.resize(n);
	m_next.resize(n);
	btVec3d* x = &system.m_x[0];
	if (!m_laplacian.update(system))
	{
		// the matrix is not positive definite (should not happen); fall back to pure inertia
		for (int i = 0; i < n; i++) x[i] = system.m_xTilde[i];
		return;
	}
	// start from the inertial positions; fixed nodes sit at their targets
	for (int i = 0; i < n; i++) x[i] = system.m_xTilde[i];
	projectColliders(system, x);

	const bool keepColliders = system.m_evaluateColliders;
	system.m_evaluateColliders = false;  // rigid contact is handled by projection here
	for (int iter = 0; iter < m_iterations; iter++)
	{
		localStep(system, x, &m_rhs[0]);
		m_laplacian.subtractFixedCoupling(x, &m_rhs[0]);
		for (int i = 0; i < n; i++)
			if (system.m_fixed[i]) m_rhs[i] = system.m_fixedTarget[i];
		m_laplacian.solve(&m_rhs[0], &m_next[0]);
		projectColliders(system, &m_next[0]);
		double change = 0;
		for (int i = 0; i < n; i++)
		{
			const double c = (m_next[i] - x[i]).maxAbs();
			if (c > change) change = c;
			x[i] = m_next[i];
		}
		m_stats.m_iterations++;
		if (m_verbose) printf("  local/global %3d: max change %.3e  E %.9e\n", iter, change, system.energy(x));
		if (m_tolerance > 0.0 && change / system.m_dt < m_tolerance)
		{
			m_stats.m_converged = true;
			break;
		}
	}
	system.gradient(x, &m_rhs[0]);
	m_stats.m_finalGradient = system.gradientNorm(&m_rhs[0]);
	m_stats.m_energy = system.energy(x);
	system.m_evaluateColliders = keepColliders;
}
