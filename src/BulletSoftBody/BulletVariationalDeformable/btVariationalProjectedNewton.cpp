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

#include "btVariationalProjectedNewton.h"
#include <stdio.h>

//
// helpers shared by all integrators
//

double btVariationalIntegrator::advanceToInitialGuess(btVariationalSystem& system, bool freeNodesToo)
{
	const int n = system.numNodes();
	btVec3dArray p;
	p.resize(n);
	// kinematic nodes first: move them to their targets as far as contact allows
	double covered = 1.0;
	bool any = false;
	for (int i = 0; i < n; i++)
	{
		p[i] = system.m_fixed[i] ? system.m_fixedTarget[i] - system.m_x[i] : btVec3d(0, 0, 0);
		if (p[i].length2() > 0.0) any = true;
	}
	if (any)
	{
		covered = system.maxStepSize(&system.m_x[0], &p[0]);
		for (int i = 0; i < n; i++) system.m_x[i] += p[i] * covered;
	}
	if (!freeNodesToo) return covered;
	// free nodes: the inertial target is the better start in free flight, the previous positions are
	// better when stiff forces dominate (large time steps, strong kinematic motion). Take the one
	// with the lower energy.
	any = false;
	for (int i = 0; i < n; i++)
	{
		p[i] = system.m_fixed[i] ? btVec3d(0, 0, 0) : system.m_xTilde[i] - system.m_x[i];
		if (p[i].length2() > 0.0) any = true;
	}
	if (!any) return covered;
	const double alpha = system.maxStepSize(&system.m_x[0], &p[0]);
	const double energyStay = system.energy(&system.m_x[0]);
	for (int i = 0; i < n; i++) p[i] = system.m_x[i] + p[i] * alpha;
	if (system.energy(&p[0]) < energyStay)
		for (int i = 0; i < n; i++) system.m_x[i] = p[i];
	return covered;
}

double btVariationalIntegrator::lineSearch(const btVariationalSystem& system, btVec3d* x, const btVec3d* p, const btVec3d* g, double alphaMax,
										   double energy0, btVec3d* scratch, int& numSteps, double* energyOut)
{
	const int n = system.numNodes();
	double slope = 0;
	for (int i = 0; i < n; i++) slope += g[i].dot(p[i]);
	double alpha = alphaMax;
	for (int ls = 0; ls < 40; ls++)
	{
		for (int i = 0; i < n; i++) scratch[i] = x[i] + p[i] * alpha;
		const double E = system.energy(scratch);
		if (E <= energy0 + 1e-4 * alpha * slope)
		{
			for (int i = 0; i < n; i++) x[i] = scratch[i];
			if (energyOut) *energyOut = E;
			return alpha;
		}
		alpha *= 0.5;
		numSteps++;
	}
	if (energyOut) *energyOut = energy0;
	return 0.0;
}

void btVariationalBlockJacobi::build(const btVariationalHessian& H)
{
	const int n = H.m_nodeBlocks.size();
	m_inverse.resize(n);
	H.diagonalBlocks(&m_inverse[0]);
	for (int i = 0; i < n; i++)
	{
		const btMat3d& D = m_inverse[i];
		const double det = D.det();
		const double scale = fabs(D.m[0][0]) + fabs(D.m[1][1]) + fabs(D.m[2][2]);
		if (det > 1e-18 * scale * scale * scale && scale > 0.0)
			m_inverse[i] = D.inverse();
		else
			m_inverse[i] = btMat3d::identity() * (scale > 0.0 ? 3.0 / scale : 1.0);
	}
}

void btVariationalBlockJacobi::apply(const btVec3d* r, btVec3d* z) const
{
	for (int i = 0; i < m_inverse.size(); i++) z[i] = m_inverse[i] * r[i];
}

//
// Projected Newton
//

void btVariationalProjectedNewton::solve(btVariationalSystem& system)
{
	const int n = system.numNodes();
	m_stats.reset();
	m_numProjectedSolves = m_numExactSolves = 0;
	if (n == 0) return;
	m_gradient.resize(n);
	m_direction.resize(n);
	m_scratch.resize(n);
	btVec3d* x = &system.m_x[0];
	const char* fixed = &system.m_fixed[0];

	double covered = advanceToInitialGuess(system, true);
	double E = system.energy(x);
	for (int iter = 0; iter < m_maxIterations; iter++)
	{
		if (covered < 1.0)
		{
			// kinematic nodes that were blocked by contact keep trying to reach their targets
			covered = advanceToInitialGuess(system, false);
			E = system.energy(x);
		}
		system.gradient(x, &m_gradient[0]);
		const double gNorm = system.gradientNorm(&m_gradient[0]);
		if (iter == 0) m_stats.m_initialGradient = gNorm;
		m_stats.m_finalGradient = gNorm;

		// First try the exact Hessian (if enabled), fall back to the eigenvalue filtered one whenever it
		// does not produce a usable descent step.
		bool exact = m_projectOnDemand && system.m_hessianFilter != BT_VARIATIONAL_FILTER_NONE;
		bool done = false;
		for (;;)
		{
			system.hessian(x, m_hessian, exact ? BT_VARIATIONAL_FILTER_NONE : system.m_hessianFilter);
			for (int i = 0; i < n; i++) m_scratch[i] = -m_gradient[i];
			const btVariationalPreconditioner* pre = m_preconditioner;
			if (!pre)
			{
				m_blockJacobi.build(m_hessian);
				pre = &m_blockJacobi;
			}
			bool indefinite = false;
			m_stats.m_linearIterations += btVariationalConjugateGradient(m_hessian, &m_scratch[0], &m_direction[0], fixed, n, pre, m_linearTolerance, m_maxLinearIterations, &indefinite);
			double step = 0, slope = 0;
			for (int i = 0; i < n; i++)
			{
				if (fixed[i])
				{
					m_direction[i].setZero();
					continue;
				}
				slope += m_gradient[i].dot(m_direction[i]);
				const double a = m_direction[i].maxAbs();
				if (a > step) step = a;
			}
			if (exact && (indefinite || !(slope < 0.0)))
			{
				exact = false;
				continue;
			}
			if (step / system.m_dt < m_tolerance && covered >= 1.0)
			{
				m_stats.m_converged = true;
				done = true;
				break;
			}
			const double alphaMax = system.maxStepSize(x, &m_direction[0]);
			const double alpha = lineSearch(system, x, &m_direction[0], &m_gradient[0], alphaMax, E, &m_scratch[0], m_stats.m_lineSearchSteps, &E);
			if (alpha == 0.0 && exact)
			{
				exact = false;
				continue;
			}
			if (exact)
				m_numExactSolves++;
			else
				m_numProjectedSolves++;
			m_stats.m_iterations++;
			if (alpha < m_stats.m_minStepSize) m_stats.m_minStepSize = alpha;
			if (m_verbose) printf("  newton %3d (%s): |g| %.3e  step %.3e  alphaMax %.3e  alpha %.3e  E %.9e  cg %d\n", iter, exact ? "exact" : "filtered", gNorm, step, alphaMax, alpha, E, m_stats.m_linearIterations);
			if (alpha == 0.0) done = true;  // no further decrease possible at working precision
			break;
		}
		if (done) break;
	}
	m_stats.m_energy = E;
	if (m_stats.m_converged)
	{
		system.gradient(x, &m_gradient[0]);
		m_stats.m_finalGradient = system.gradientNorm(&m_gradient[0]);
	}
}
