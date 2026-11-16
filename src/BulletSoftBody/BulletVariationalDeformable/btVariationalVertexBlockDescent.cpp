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

#include "btVariationalVertexBlockDescent.h"
#include "LinearMath/btThreads.h"
#include <stdio.h>

btVariationalVertexBlockDescent::btVariationalVertexBlockDescent()
	: m_iterations(20),
	  m_initialization(INIT_ADAPTIVE),
	  m_tolerance(0.0),
	  m_useChebyshev(false),
	  m_chebyshevRho(0.93),
	  m_rayleighDamping(0.0),
	  m_parallel(true),
	  m_coloringSignature(-1.0),
	  m_lastTopologyVersion(-1),
	  m_maxConstraintError(0.0)
{
}

void btVariationalVertexBlockDescent::buildColoring(const btVariationalSystem& system)
{
	// vertices that share an element (or a two-node attachment) must not have the same color
	const int n = system.numNodes();
	int numPairAttachments = 0;
	for (int a = 0; a < system.m_attachments.size(); a++)
		if (system.m_attachments[a].m_nodeB >= 0) numPairAttachments++;
	int numFixed = 0;
	for (int i = 0; i < n; i++) numFixed += system.m_fixed[i] ? 1 : 0;
	const double signature = system.m_topologyVersion * 1e6 + numPairAttachments * 1e3 + numFixed + n * 1e-3;
	if (signature == m_coloringSignature) return;
	m_coloringSignature = signature;

	btAlignedObjectArray<int> color, used;
	color.resize(n, -1);
	int numColors = 0;
	for (int i = 0; i < n; i++)
	{
		if (system.m_fixed[i]) continue;
#define BT_MARK_NEIGHBOR(node)                                    \
	{                                                             \
		const int c = color[node];                                \
		if (c >= 0) used[c] = i;                                  \
	}
		for (int p = system.m_nodeTetStart[i]; p < system.m_nodeTetStart[i + 1]; p++)
		{
			const btVariationalTet& tet = system.m_tets[system.m_nodeTets[p] >> 2];
			for (int j = 0; j < 4; j++) BT_MARK_NEIGHBOR(tet.m_n[j]);
		}
		for (int p = system.m_nodeSpringStart[i]; p < system.m_nodeSpringStart[i + 1]; p++)
		{
			const btVariationalSpring& sp = system.m_springs[system.m_nodeSprings[p] >> 1];
			BT_MARK_NEIGHBOR(sp.m_n[0]);
			BT_MARK_NEIGHBOR(sp.m_n[1]);
		}
		for (int p = system.m_nodeAttachStart[i]; p < system.m_nodeAttachStart[i + 1]; p++)
		{
			const btVariationalAttachment& at = system.m_attachments[system.m_nodeAttach[p] >> 1];
			BT_MARK_NEIGHBOR(at.m_nodeA);
			if (at.m_nodeB >= 0) BT_MARK_NEIGHBOR(at.m_nodeB);
		}
#undef BT_MARK_NEIGHBOR
		int c = 0;
		while (c < numColors && used[c] == i) c++;
		if (c == numColors)
		{
			numColors++;
			used.push_back(-1);
		}
		color[i] = c;
	}
	m_colorStart.resize(0);
	m_colorStart.resize(numColors + 1, 0);
	for (int i = 0; i < n; i++)
		if (color[i] >= 0) m_colorStart[color[i] + 1]++;
	for (int c = 0; c < numColors; c++) m_colorStart[c + 1] += m_colorStart[c];
	m_colorNodes.resize(m_colorStart[numColors]);
	btAlignedObjectArray<int> fill;
	fill.resize(numColors);
	for (int c = 0; c < numColors; c++) fill[c] = m_colorStart[c];
	for (int i = 0; i < n; i++)
		if (color[i] >= 0) m_colorNodes[fill[color[i]]++] = i;
}

double btVariationalVertexBlockDescent::updateVertex(const btVariationalSystem& system, btVec3d* x, int node) const
{
	btVec3d force;
	btMat3d hess, elastic;
	system.vertexForceAndHessian(node, x, force, hess, m_rayleighDamping > 0.0 ? &elastic : 0);
	if (m_rayleighDamping > 0.0)
	{
		// f -= gamma K v,  H += gamma / h K
		const btVec3d velocity = (x[node] - system.m_xPrev[node]) * (1.0 / system.m_dt);
		force -= (elastic * velocity) * m_rayleighDamping;
		hess += elastic * (m_rayleighDamping / system.m_dt);
	}
	btVec3d dx;
	if (!btSolveSpd3d(hess, force, dx))
	{
		// degenerate block (cannot happen with positive mass); fall back to a scaled gradient step
		const double scale = hess.trace();
		if (!(scale > 0.0)) return 0.0;
		dx = force * (3.0 / scale);
	}
	const double alpha = system.vertexMaxStepSize(node, x, dx, true);
	if (alpha < 1.0) dx *= alpha;
	x[node] += dx;
	return dx.length();
}

namespace
{
struct btVbdColorBody : public btIParallelSumBody
{
	const btVariationalVertexBlockDescent* m_vbd;
	const btVariationalSystem* m_system;
	btVec3d* m_x;
	const int* m_nodes;
	virtual btScalar sumLoop(int iBegin, int iEnd) const
	{
		double largest = 0;
		for (int k = iBegin; k < iEnd; k++)
		{
			const double d = m_vbd->updateVertex(*m_system, m_x, m_nodes[k]);
			if (d > largest) largest = d;
		}
		return btScalar(largest);
	}
};
}  // namespace

double btVariationalVertexBlockDescent::sweep(const btVariationalSystem& system, btVec3d* x)
{
	// vertices of a color share no element, so they do not read each other's positions
	bool parallel = m_parallel;
	for (int c = 0; c < system.m_colliders.size() && parallel; c++)
		if (!system.m_colliders[c].isConvexExact()) parallel = false;  // Bullet's sparse SDF cache is not thread safe
	double largest = 0;
	for (int c = 0; c + 1 < m_colorStart.size(); c++)
	{
		const int begin = m_colorStart[c], end = m_colorStart[c + 1];
		if (parallel && end - begin > 64)
		{
			// btParallelSum adds up the per-thread results; the sweep only needs a bound, so do the
			// reduction by hand on chunks
			btVbdColorBody body;
			body.m_vbd = this;
			body.m_system = &system;
			body.m_x = x;
			body.m_nodes = &m_colorNodes[0] + begin;
			const double sum = btParallelSum(0, end - begin, 64, body);
			if (sum > largest) largest = sum;  // an upper bound of the largest displacement
		}
		else
		{
			for (int k = begin; k < end; k++)
			{
				const double d = updateVertex(system, x, m_colorNodes[k]);
				if (d > largest) largest = d;
			}
		}
	}
	return largest;
}

void btVariationalVertexBlockDescent::solve(btVariationalSystem& system)
{
	const int n = system.numNodes();
	m_stats.reset();
	if (n == 0) return;
	btVec3d* x = &system.m_x[0];
	const bool firstStep = system.m_topologyVersion != m_lastTopologyVersion;
	m_lastTopologyVersion = system.m_topologyVersion;
	buildColoring(system);
	m_al.warmStart(system);

	// initial guess
	const bool barrier = system.m_contact != 0 || system.m_colliderModel == btVariationalSystem::COLLIDER_BARRIER;
	{
		advanceToInitialGuess(system, false);  // kinematic nodes
		const double h = system.m_dt;
		const double gNorm = system.m_gravity.length();
		const btVec3d gDir = gNorm > 0.0 ? system.m_gravity / gNorm : btVec3d(0, 0, 0);
		m_prevIterate.resize(n);
		for (int i = 0; i < n; i++)
		{
			btVec3d target = x[i];
			if (!system.m_fixed[i])
			{
				if (m_initialization == INIT_INERTIA)
					target = system.m_xTilde[i];
				else if (m_initialization == INIT_ADAPTIVE)
				{
					// acceleration of the previous step along gravity, clamped to [0, |g|]
					// (the very first step has no previous acceleration and uses full gravity)
					double a = firstStep ? gNorm : ((system.m_v[i] - system.m_vPrev[i]) * (1.0 / h)).dot(gDir);
					a = a < gNorm ? a : gNorm;
					a = a > 1e-5 ? a : 0.0;
					target = system.m_xPrev[i] + system.m_v[i] * h + gDir * (a * h * h);
				}
			}
			m_prevIterate[i] = target;
		}
		// with barrier contact the guess is approached node by node, as far as it stays feasible
		moveNodesTowards(system, &m_prevIterate[0]);
	}

	const bool chebyshev = m_useChebyshev && !barrier;
	if (chebyshev)
	{
		m_prevIterate.resize(n);
		m_prevPrevIterate.resize(n);
		for (int i = 0; i < n; i++) m_prevPrevIterate[i] = m_prevIterate[i] = x[i];
	}
	double omega = 1.0;
	// With deformable-deformable contact every sweep runs inside a trust region: the candidate pairs
	// are collected for a ball around each vertex, and no vertex may leave its ball during the sweep.
	// The candidates are rebuilt only when the vertices have used up half of their ball.
	// The ball starts small (the contact model enlarges it to a multiple of its activation distance)
	// and grows only if the vertices actually need the room, up to a limit: the number of candidates
	// rises steeply with the radius, and in the step of an impact an unlimited ball made single
	// steps twenty times more expensive than the rest. A vertex that needs to go further simply gets
	// there over several rebuilds.
	const double maxGrowth = 4.0;
	double radius = 0, travelled = 0, smallestRadius = 0;
	bool rebuild = true;
	for (int iter = 0; iter < m_iterations; iter++)
	{
		if (system.m_contact && rebuild)
		{
			radius = system.m_contact->prepareLocalSteps(x, radius);
			if (smallestRadius == 0.0) smallestRadius = radius;
			travelled = 0;
			rebuild = false;
		}
		const double change = sweep(system, x);
		travelled += change;
		if (travelled > 0.5 * radius)
		{
			radius = btMin(2.0 * travelled, maxGrowth * smallestRadius);
			rebuild = true;
		}
		m_stats.m_iterations++;
		if (chebyshev)
		{
			// x^(k+1) = omega_{k+1} (xbar^(k+1) - x^(k-1)) + x^(k-1)
			omega = iter == 0 ? 1.0 : (iter == 1 ? 2.0 / (2.0 - m_chebyshevRho * m_chebyshevRho) : 4.0 / (4.0 - m_chebyshevRho * m_chebyshevRho * omega));
			for (int i = 0; i < n; i++)
			{
				if (!system.m_fixed[i]) x[i] = (x[i] - m_prevPrevIterate[i]) * omega + m_prevPrevIterate[i];
				m_prevPrevIterate[i] = m_prevIterate[i];
				m_prevIterate[i] = x[i];
			}
		}
		m_maxConstraintError = m_al.updateDuals(system, x);
		if (m_verbose) printf("  vbd sweep %3d: max change %.3e  E %.9e  constraint error %.3e\n", iter, change, system.energy(x), m_maxConstraintError);
		if (m_tolerance > 0.0 && change / system.m_dt < m_tolerance)
		{
			m_stats.m_converged = true;
			break;
		}
	}
	btVariationalAugmentedLagrangian::finishStep(system);
	m_stats.m_energy = system.energy(x);
}
