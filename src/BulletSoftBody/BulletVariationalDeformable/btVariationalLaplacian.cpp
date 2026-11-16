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

#include "btVariationalLaplacian.h"

double btVariationalLaplacian::computeSignature(const btVariationalSystem& system) const
{
	// a cheap fingerprint of everything the matrix depends on
	double sig = system.m_dt * 1.0e3 + system.m_topologyVersion * 7.0 + system.numNodes() * 13.0;
	for (int b = 0; b < system.m_bodies.size(); b++)
	{
		const btVariationalBody& body = system.m_bodies[b];
		sig += (b + 1) * 1e-3 * (body.m_material->laplacianStiffness() + 3.0 * body.m_springStiffness + 5.0 * body.m_bendingStiffness);
	}
	for (int a = 0; a < system.m_attachments.size(); a++)
	{
		const btVariationalAttachment& at = system.m_attachments[a];
		sig += (a + 1) * (at.m_nodeA * 3.0 + at.m_nodeB * 5.0 + (at.m_hard ? at.m_penalty : at.m_stiffness) * 1e-3);
	}
	for (int i = 0; i < system.numNodes(); i++)
	{
		if (system.m_fixed[i]) sig += (i + 1) * 0.37;
		sig += system.m_mass[i] * 1e-2 * ((i % 97) + 1);
	}
	return sig;
}

bool btVariationalLaplacian::update(const btVariationalSystem& system)
{
	const double sig = computeSignature(system);
	if (m_valid && sig == m_signature) return true;
	const int n = system.numNodes();
	const double invH2 = 1.0 / (system.m_dt * system.m_dt);
	btAlignedObjectArray<int> rows, cols;
	btAlignedObjectArray<double> vals;
	rows.reserve(n + 16 * system.m_tets.size());
	cols.reserve(n + 16 * system.m_tets.size());
	vals.reserve(n + 16 * system.m_tets.size());
	for (int i = 0; i < n; i++)
	{
		rows.push_back(i);
		cols.push_back(i);
		vals.push_back(system.m_mass[i] * invH2);
	}
	for (int t = 0; t < system.m_tets.size(); t++)
	{
		const btVariationalTet& tet = system.m_tets[t];
		const double w = system.m_bodies[tet.m_body].m_material->laplacianStiffness() * tet.m_restVolume;
		for (int j = 0; j < 4; j++)
			for (int k = 0; k < 4; k++)
			{
				rows.push_back(tet.m_n[j]);
				cols.push_back(tet.m_n[k]);
				vals.push_back(w * tet.m_g[j].dot(tet.m_g[k]));
			}
	}
	for (int s = 0; s < system.m_springs.size(); s++)
	{
		const btVariationalSpring& sp = system.m_springs[s];
		for (int j = 0; j < 2; j++)
			for (int k = 0; k < 2; k++)
			{
				rows.push_back(sp.m_n[j]);
				cols.push_back(sp.m_n[k]);
				vals.push_back(j == k ? sp.m_stiffness : -sp.m_stiffness);
			}
	}
	for (int a = 0; a < system.m_attachments.size(); a++)
	{
		const btVariationalAttachment& at = system.m_attachments[a];
		const double k = at.m_hard ? at.m_penalty : at.m_stiffness;
		const int nodes[2] = {at.m_nodeA, at.m_nodeB};
		for (int j = 0; j < 2; j++)
			for (int l = 0; l < 2; l++)
			{
				if (nodes[j] < 0 || nodes[l] < 0) continue;
				rows.push_back(nodes[j]);
				cols.push_back(nodes[l]);
				vals.push_back(j == l ? k : -k);
			}
	}
	m_full.setFromTriplets(n, rows, cols, vals);

	// reduced matrix: fixed rows and columns become identity, their coupling is remembered
	m_fixed.resize(n);
	for (int i = 0; i < n; i++) m_fixed[i] = system.m_fixed[i];
	m_couplingFree.resize(0);
	m_couplingFixed.resize(0);
	m_couplingValue.resize(0);
	rows.resize(0);
	cols.resize(0);
	vals.resize(0);
	for (int i = 0; i < n; i++)
	{
		if (m_fixed[i])
		{
			rows.push_back(i);
			cols.push_back(i);
			vals.push_back(1.0);
			continue;
		}
		bool hasDiagonal = false;
		for (int p = m_full.m_rowStart[i]; p < m_full.m_rowStart[i + 1]; p++)
		{
			const int j = m_full.m_cols[p];
			if (m_fixed[j])
			{
				m_couplingFree.push_back(i);
				m_couplingFixed.push_back(j);
				m_couplingValue.push_back(m_full.m_vals[p]);
				continue;
			}
			rows.push_back(i);
			cols.push_back(j);
			vals.push_back(m_full.m_vals[p]);
			if (j == i && m_full.m_vals[p] > 0.0) hasDiagonal = true;
		}
		if (!hasDiagonal)
		{
			// a massless node that nothing is attached to; keep the matrix regular
			rows.push_back(i);
			cols.push_back(i);
			vals.push_back(1.0);
		}
	}
	btVariationalSparseMatrix reduced;
	reduced.setFromTriplets(n, rows, cols, vals);
	m_valid = m_cholesky.factorize(reduced);
	m_signature = sig;
	m_numFactorizations++;
	return m_valid;
}

void btVariationalLaplacian::subtractFixedCoupling(const btVec3d* x, btVec3d* rhs) const
{
	for (int c = 0; c < m_couplingFree.size(); c++) rhs[m_couplingFree[c]] -= x[m_couplingFixed[c]] * m_couplingValue[c];
}
