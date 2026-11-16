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

#include "btVariationalIPC.h"
#include <stdio.h>

btVariationalIPC::btVariationalIPC()
	: m_dhat(0.0),
	  m_dhatRelative(1e-3),
	  m_stiffness(0.0),
	  m_adaptiveStiffness(true),
	  m_selfCollision(true),
	  m_friction(true),
	  m_frictionEpsilon(1e-2),
	  m_initialStateValid(true),
	  m_warned(false),
	  m_system(0),
	  m_dhatActive(0.0),
	  m_kappa(0.0),
	  m_kappaMax(0.0),
	  m_sceneDiagonal(1.0),
	  m_prevMinDistance2(0.0),
	  m_muted(false),
	  m_localRadius(0.0)
{
}

void btVariationalIPC::setCollisionGroup(const btSoftBody* psb, int group)
{
	for (int i = 0; i < m_groupBodies.size(); i++)
		if (m_groupBodies[i] == psb)
		{
			m_groupIds[i] = group;
			return;
		}
	m_groupBodies.push_back(psb);
	m_groupIds.push_back(group);
}

namespace
{
///Uniform grid with hashed cells over a set of axis aligned boxes
struct btHashGrid
{
	double m_invCell;
	int m_mask;
	btAlignedObjectArray<int> m_start, m_items;

	static int cellOf(double v, double invCell) { return (int)floor(v * invCell); }
	int bucket(int ix, int iy, int iz) const
	{
		return (int)(((unsigned int)ix * 73856093u) ^ ((unsigned int)iy * 19349663u) ^ ((unsigned int)iz * 83492791u)) & m_mask;
	}
	void build(int count, const btVec3d* boxMin, const btVec3d* boxMax, double cell)
	{
		m_invCell = 1.0 / cell;
		int size = 16;
		while (size < 2 * count) size *= 2;
		m_mask = size - 1;
		m_start.resize(0);
		m_start.resize(size + 1, 0);
		for (int pass = 0; pass < 2; pass++)
		{
			for (int i = 0; i < count; i++)
			{
				const int x0 = cellOf(boxMin[i][0], m_invCell), x1 = cellOf(boxMax[i][0], m_invCell);
				const int y0 = cellOf(boxMin[i][1], m_invCell), y1 = cellOf(boxMax[i][1], m_invCell);
				const int z0 = cellOf(boxMin[i][2], m_invCell), z1 = cellOf(boxMax[i][2], m_invCell);
				for (int ix = x0; ix <= x1; ix++)
					for (int iy = y0; iy <= y1; iy++)
						for (int iz = z0; iz <= z1; iz++)
						{
							const int b = bucket(ix, iy, iz);
							if (pass == 0)
								m_start[b + 1]++;
							else
								m_items[m_start[b]++] = i;
						}
			}
			if (pass == 0)
			{
				for (int b = 0; b < size; b++) m_start[b + 1] += m_start[b];
				m_items.resize(m_start[size]);
			}
			else
			{
				// the fill pass advanced every start to the end of its bucket; shift back
				for (int b = size; b > 0; b--) m_start[b] = m_start[b - 1];
				m_start[0] = 0;
			}
		}
	}
};

inline bool btBoxesOverlap(const btVec3d& aMin, const btVec3d& aMax, const btVec3d& bMin, const btVec3d& bMax)
{
	return aMin[0] <= bMax[0] && bMin[0] <= aMax[0] && aMin[1] <= bMax[1] && bMin[1] <= aMax[1] && aMin[2] <= bMax[2] && bMin[2] <= aMax[2];
}
}  // namespace

void btVariationalIPC::buildCandidates(const btVec3d* x, const btVec3d* p, double inflation) const
{
	const btVariationalSystem& sys = *m_system;
	const int n = sys.numNodes();
	const int numTris = sys.m_triangles.size(), numEdges = sys.m_edges.size();
	m_pairs.resize(0);
	m_boxMin.resize(n);
	m_boxMax.resize(n);
	for (int i = 0; i < n; i++)
	{
		btVec3d lo = x[i], hi = x[i];
		if (p)
			for (int c = 0; c < 3; c++)
			{
				const double e = x[i][c] + p[i][c];
				if (e < lo[c]) lo[c] = e;
				if (e > hi[c]) hi[c] = e;
			}
		m_boxMin[i] = lo - btVec3d(inflation, inflation, inflation);
		m_boxMax[i] = hi + btVec3d(inflation, inflation, inflation);
	}
	if (numTris > 0)
	{
		btVec3dArray triMin, triMax, edgeMin, edgeMax;
		triMin.resize(numTris);
		triMax.resize(numTris);
		edgeMin.resize(numEdges);
		edgeMax.resize(numEdges);
		double extent = 0;
		for (int t = 0; t < numTris; t++)
		{
			const int* tn = sys.m_triangles[t].m_n;
			btVec3d lo = m_boxMin[tn[0]], hi = m_boxMax[tn[0]];
			for (int k = 1; k < 3; k++)
				for (int c = 0; c < 3; c++)
				{
					if (m_boxMin[tn[k]][c] < lo[c]) lo[c] = m_boxMin[tn[k]][c];
					if (m_boxMax[tn[k]][c] > hi[c]) hi[c] = m_boxMax[tn[k]][c];
				}
			triMin[t] = lo;
			triMax[t] = hi;
			extent += (hi - lo).maxAbs();
		}
		for (int e = 0; e < numEdges; e++)
		{
			const int* en = sys.m_edges[e].m_n;
			for (int c = 0; c < 3; c++)
			{
				edgeMin[e][c] = m_boxMin[en[0]][c] < m_boxMin[en[1]][c] ? m_boxMin[en[0]][c] : m_boxMin[en[1]][c];
				edgeMax[e][c] = m_boxMax[en[0]][c] > m_boxMax[en[1]][c] ? m_boxMax[en[0]][c] : m_boxMax[en[1]][c];
			}
		}
		const double cell = extent > 0.0 ? extent / numTris : 1.0;
		btHashGrid triGrid, edgeGrid;
		triGrid.build(numTris, &triMin[0], &triMax[0], cell);
		edgeGrid.build(numEdges, &edgeMin[0], &edgeMax[0], cell);
		btAlignedObjectArray<int> stamp;

		// point-triangle candidates
		stamp.resize(numTris, -1);
		for (int i = 0; i < n; i++)
		{
			if (!sys.m_surfaceNode[i]) continue;
			const int x0 = btHashGrid::cellOf(m_boxMin[i][0], triGrid.m_invCell), x1 = btHashGrid::cellOf(m_boxMax[i][0], triGrid.m_invCell);
			const int y0 = btHashGrid::cellOf(m_boxMin[i][1], triGrid.m_invCell), y1 = btHashGrid::cellOf(m_boxMax[i][1], triGrid.m_invCell);
			const int z0 = btHashGrid::cellOf(m_boxMin[i][2], triGrid.m_invCell), z1 = btHashGrid::cellOf(m_boxMax[i][2], triGrid.m_invCell);
			for (int ix = x0; ix <= x1; ix++)
				for (int iy = y0; iy <= y1; iy++)
					for (int iz = z0; iz <= z1; iz++)
					{
						const int b = triGrid.bucket(ix, iy, iz);
						for (int q = triGrid.m_start[b]; q < triGrid.m_start[b + 1]; q++)
						{
							const int t = triGrid.m_items[q];
							if (stamp[t] == i) continue;
							stamp[t] = i;
							const btVariationalTriangle& tri = sys.m_triangles[t];
							if (tri.m_n[0] == i || tri.m_n[1] == i || tri.m_n[2] == i) continue;
							if (excluded(tri.m_body, sys.m_nodeBody[i])) continue;
							if (!btBoxesOverlap(m_boxMin[i], m_boxMax[i], triMin[t], triMax[t])) continue;
							if (sys.m_fixed[i] && sys.m_fixed[tri.m_n[0]] && sys.m_fixed[tri.m_n[1]] && sys.m_fixed[tri.m_n[2]]) continue;
							Pair pr;
							pr.m_n[0] = i;
							pr.m_n[1] = tri.m_n[0];
							pr.m_n[2] = tri.m_n[1];
							pr.m_n[3] = tri.m_n[2];
							pr.m_edgeEdge = false;
							pr.m_epsX = 0.0;
							pr.m_refDistance = -1.0;
							m_pairs.push_back(pr);
						}
					}
		}
		// edge-edge candidates
		stamp.resize(0);
		stamp.resize(numEdges, -1);
		for (int e = 0; e < numEdges; e++)
		{
			const btVariationalEdge& ea = sys.m_edges[e];
			const int x0 = btHashGrid::cellOf(edgeMin[e][0], edgeGrid.m_invCell), x1 = btHashGrid::cellOf(edgeMax[e][0], edgeGrid.m_invCell);
			const int y0 = btHashGrid::cellOf(edgeMin[e][1], edgeGrid.m_invCell), y1 = btHashGrid::cellOf(edgeMax[e][1], edgeGrid.m_invCell);
			const int z0 = btHashGrid::cellOf(edgeMin[e][2], edgeGrid.m_invCell), z1 = btHashGrid::cellOf(edgeMax[e][2], edgeGrid.m_invCell);
			for (int ix = x0; ix <= x1; ix++)
				for (int iy = y0; iy <= y1; iy++)
					for (int iz = z0; iz <= z1; iz++)
					{
						const int b = edgeGrid.bucket(ix, iy, iz);
						for (int q = edgeGrid.m_start[b]; q < edgeGrid.m_start[b + 1]; q++)
						{
							const int f = edgeGrid.m_items[q];
							if (f <= e || stamp[f] == e) continue;
							stamp[f] = e;
							const btVariationalEdge& eb = sys.m_edges[f];
							if (ea.m_n[0] == eb.m_n[0] || ea.m_n[0] == eb.m_n[1] || ea.m_n[1] == eb.m_n[0] || ea.m_n[1] == eb.m_n[1]) continue;
							if (excluded(ea.m_body, eb.m_body)) continue;
							if (!btBoxesOverlap(edgeMin[e], edgeMax[e], edgeMin[f], edgeMax[f])) continue;
							if (sys.m_fixed[ea.m_n[0]] && sys.m_fixed[ea.m_n[1]] && sys.m_fixed[eb.m_n[0]] && sys.m_fixed[eb.m_n[1]]) continue;
							Pair pr;
							pr.m_n[0] = ea.m_n[0];
							pr.m_n[1] = ea.m_n[1];
							pr.m_n[2] = eb.m_n[0];
							pr.m_n[3] = eb.m_n[1];
							pr.m_edgeEdge = true;
							pr.m_refDistance = -1.0;
							// mollifier threshold from the edge lengths at the start of the step
							pr.m_epsX = 1e-3 * (sys.m_xPrev[ea.m_n[0]] - sys.m_xPrev[ea.m_n[1]]).length2() * (sys.m_xPrev[eb.m_n[0]] - sys.m_xPrev[eb.m_n[1]]).length2();
							m_pairs.push_back(pr);
						}
					}
		}
	}
	// node -> pairs, encoded as 4 * pair + slot
	m_nodePairStart.resize(0);
	m_nodePairStart.resize(n + 1, 0);
	for (int q = 0; q < m_pairs.size(); q++)
		for (int j = 0; j < 4; j++) m_nodePairStart[m_pairs[q].m_n[j] + 1]++;
	for (int i = 0; i < n; i++) m_nodePairStart[i + 1] += m_nodePairStart[i];
	m_nodePairs.resize(4 * m_pairs.size());
	btAlignedObjectArray<int> fill;
	fill.resize(n);
	for (int i = 0; i < n; i++) fill[i] = m_nodePairStart[i];
	for (int q = 0; q < m_pairs.size(); q++)
		for (int j = 0; j < 4; j++) m_nodePairs[fill[m_pairs[q].m_n[j]]++] = 4 * q + j;
}

double btVariationalIPC::distanceLowerBound(const Pair& pair, const btVec3d* x) const
{
	if (pair.m_refDistance < 0.0 || m_reference.size() == 0) return -1.0;
	double off[4];
	for (int j = 0; j < 4; j++) off[j] = (x[pair.m_n[j]] - m_reference[pair.m_n[j]]).length();
	// the two primitives are (0 | 1 2 3) for point-triangle and (0 1 | 2 3) for edge-edge
	const double a = pair.m_edgeEdge ? (off[0] > off[1] ? off[0] : off[1]) : off[0];
	double b = off[2] > off[3] ? off[2] : off[3];
	if (!pair.m_edgeEdge && off[1] > b) b = off[1];
	return pair.m_refDistance - a - b;
}

bool btVariationalIPC::evaluatePair(const Pair& pair, const btVec3d* x, int order, PairValue& out, int onlySlot) const
{
	const btVec3d xs[4] = {x[pair.m_n[0]], x[pair.m_n[1]], x[pair.m_n[2]], x[pair.m_n[3]]};
	const btVariationalClosest closest = pair.m_edgeEdge ? btEdgeEdgeClosest(xs[0], xs[1], xs[2], xs[3]) : btPointTriangleClosest(xs[0], xs[1], xs[2], xs[3]);
	btVariationalPairDistance dist;
	btEvaluatePairDistance(xs, closest, 0, dist);
	const double s = dist.m_distanceSquared, sHat = m_dhatActive * m_dhatActive;
	if (s >= sHat) return false;
	if (!(s > 0.0))
	{
		out.m_energy = 1e300;  // touching: infeasible
		return true;
	}
	const double t = s - sHat, lg = log(s / sHat);
	const double B = -t * t * lg;
	// mollifier of nearly parallel edge-edge pairs
	double m = 1.0, dm = 0.0, d2m = 0.0, cross2 = 0.0;
	btVec3d u, v;
	if (pair.m_edgeEdge)
	{
		u = xs[1] - xs[0];
		v = xs[3] - xs[2];
		cross2 = u.cross(v).length2();
		if (cross2 < pair.m_epsX)
		{
			const double r = cross2 / pair.m_epsX;
			m = (2.0 - r) * r;
			dm = (2.0 - 2.0 * r) / pair.m_epsX;
			d2m = -2.0 / (pair.m_epsX * pair.m_epsX);
		}
	}
	out.m_energy = m_kappa * m * B;
	if (order < 1) return true;
	btEvaluatePairDistance(xs, closest, order, dist, onlySlot);
	const double dB = -2.0 * t * lg - t * t / s;
	const double d2B = -2.0 * lg - 4.0 * t / s + t * t / (s * s);
	const int count = closest.count();
	for (int j = 0; j < 4; j++) out.m_grad[j] = dist.m_grad[j] * (m_kappa * m * dB);
	if (order >= 2)
	{
		for (int j = 0; j < 4; j++)
			for (int k = 0; k < 4; k++) out.m_hess[j][k].setZero();
		for (int a = 0; a < count; a++)
			for (int b = 0; b < count; b++)
			{
				const int j = closest.m_slot[a], k = closest.m_slot[b];
				if (onlySlot >= 0 && (j != onlySlot || k != onlySlot)) continue;
				btMat3d M = dist.m_hess[j][k] * (m * dB);
				M.addScaledOuter(m * d2B, dist.m_grad[j], dist.m_grad[k]);
				out.m_hess[j][k] = M * m_kappa;
			}
	}
	if (dm != 0.0)
	{
		// c = |u x v|^2 and its derivatives with respect to (ea0, ea1, eb0, eb1)
		const double uv = u.dot(v), uu = u.length2(), vv = v.length2();
		const btVec3d gu = (u * vv - v * uv) * 2.0, gv = (v * uu - u * uv) * 2.0;
		const btVec3d gc[4] = {-gu, gu, -gv, gv};
		for (int j = 0; j < 4; j++) out.m_grad[j] += gc[j] * (m_kappa * dm * B);
		if (order >= 2)
		{
			btMat3d Huu = btMat3d::identity() * (2.0 * vv), Hvv = btMat3d::identity() * (2.0 * uu), Huv = btMat3d::identity() * (-2.0 * uv);
			Huu.addScaledOuter(-2.0, v, v);
			Hvv.addScaledOuter(-2.0, u, u);
			Huv.addScaledOuter(4.0, u, v);
			Huv.addScaledOuter(-2.0, v, u);
			const btMat3d Hvu = Huv.transpose();
			const double sign[4] = {-1, 1, -1, 1};
			for (int j = 0; j < 4; j++)
				for (int k = 0; k < 4; k++)
				{
					if (onlySlot >= 0 && (j != onlySlot || k != onlySlot)) continue;
					const btMat3d& Hc = j < 2 ? (k < 2 ? Huu : Huv) : (k < 2 ? Hvu : Hvv);
					btMat3d M = Hc * (sign[j] * sign[k] * dm * B);
					M.addScaledOuter(d2m * B, gc[j], gc[k]);
					M.addScaledOuter(dm * dB, gc[j], dist.m_grad[k]);
					M.addScaledOuter(dm * dB, dist.m_grad[j], gc[k]);
					out.m_hess[j][k] += M * m_kappa;
				}
		}
	}
	return true;
}

void btVariationalIPC::frictionKernel(const FrictionPair& fp, const btVec3d* x, double& energy, btVec3d* grad, btMat3d* hess) const
{
	// relative tangential displacement of the two contact points since the start of the step
	btVec3d u(0, 0, 0);
	for (int j = 0; j < 4; j++)
		if (fp.m_w[j] != 0.0) u += (x[fp.m_n[j]] - m_system->m_xPrev[fp.m_n[j]]) * fp.m_w[j];
	const btVec3d uT = u - fp.m_normal * u.dot(fp.m_normal);
	const double y = uT.length();
	const double eps = m_frictionEpsilon * m_system->m_dt;
	if (y >= eps)
	{
		energy += fp.m_force * y;
		if (grad) *grad = uT * (fp.m_force / y);
		if (hess)
		{
			hess->setZero();
			hess->addScaledIdentity(fp.m_force / y);
			hess->addScaledOuter(-fp.m_force / y, fp.m_normal, fp.m_normal);
			hess->addScaledOuter(-fp.m_force / (y * y * y), uT, uT);
		}
	}
	else
	{
		energy += fp.m_force * (-y * y * y / (3.0 * eps * eps) + y * y / eps + eps / 3.0);
		const double f1OverY = -y / (eps * eps) + 2.0 / eps;
		if (grad) *grad = uT * (fp.m_force * f1OverY);
		if (hess)
		{
			hess->setZero();
			hess->addScaledIdentity(fp.m_force * f1OverY);
			hess->addScaledOuter(-fp.m_force * f1OverY, fp.m_normal, fp.m_normal);
			if (y > 0.0) hess->addScaledOuter(-fp.m_force / (eps * eps * y), uT, uT);
		}
	}
}

void btVariationalIPC::buildFrictionPairs(const btVec3d* x)
{
	const btVariationalSystem& sys = *m_system;
	m_frictionPairs.resize(0);
	const double sHat = m_dhatActive * m_dhatActive;
	for (int q = 0; q < m_pairs.size(); q++)
	{
		const Pair& pair = m_pairs[q];
		const double mu = sqrt(sys.m_bodies[sys.m_nodeBody[pair.m_n[0]]].m_friction * sys.m_bodies[sys.m_nodeBody[pair.m_n[pair.m_edgeEdge ? 2 : 1]]].m_friction);
		if (!(mu > 0.0)) continue;
		const btVec3d xs[4] = {x[pair.m_n[0]], x[pair.m_n[1]], x[pair.m_n[2]], x[pair.m_n[3]]};
		const btVariationalClosest closest = pair.m_edgeEdge ? btEdgeEdgeClosest(xs[0], xs[1], xs[2], xs[3]) : btPointTriangleClosest(xs[0], xs[1], xs[2], xs[3]);
		btVariationalPairDistance dist;
		btEvaluatePairDistance(xs, closest, 0, dist);
		const double s = dist.m_distanceSquared;
		if (s >= sHat || !(s > 0.0)) continue;
		double m = 1.0;
		if (pair.m_edgeEdge)
		{
			const double cross2 = (xs[1] - xs[0]).cross(xs[3] - xs[2]).length2();
			if (cross2 < pair.m_epsX) m = (2.0 - cross2 / pair.m_epsX) * cross2 / pair.m_epsX;
		}
		const double t = s - sHat, dB = -2.0 * t * log(s / sHat) - t * t / s;
		const double d = sqrt(s);
		FrictionPair fp;
		btVec3d c(0, 0, 0);
		for (int j = 0; j < 4; j++)
		{
			fp.m_n[j] = pair.m_n[j];
			fp.m_w[j] = dist.m_weights[j];
			c += xs[j] * dist.m_weights[j];
		}
		fp.m_normal = c.normalized();
		fp.m_force = mu * (-m_kappa * m * dB * 2.0 * d);  // mu times |d(kappa b)/dd|
		if (fp.m_force > 0.0) m_frictionPairs.push_back(fp);
	}
	const int n = sys.numNodes();
	m_nodeFrictionStart.resize(0);
	m_nodeFrictionStart.resize(n + 1, 0);
	for (int q = 0; q < m_frictionPairs.size(); q++)
		for (int j = 0; j < 4; j++)
			if (m_frictionPairs[q].m_w[j] != 0.0) m_nodeFrictionStart[m_frictionPairs[q].m_n[j] + 1]++;
	for (int i = 0; i < n; i++) m_nodeFrictionStart[i + 1] += m_nodeFrictionStart[i];
	m_nodeFriction.resize(m_nodeFrictionStart[n]);
	btAlignedObjectArray<int> fill;
	fill.resize(n);
	for (int i = 0; i < n; i++) fill[i] = m_nodeFrictionStart[i];
	for (int q = 0; q < m_frictionPairs.size(); q++)
		for (int j = 0; j < 4; j++)
			if (m_frictionPairs[q].m_w[j] != 0.0) m_nodeFriction[fill[m_frictionPairs[q].m_n[j]]++] = 4 * q + j;
}

void btVariationalIPC::beginStep(const btVariationalSystem& system)
{
	m_system = &system;
	const int n = system.numNodes();
	m_frictionPairs.resize(0);
	m_nodeFrictionStart.resize(0);
	m_nodeFrictionStart.resize(n + 1, 0);
	m_reference.resize(0);
	if (n == 0) return;
	const btVec3d* x = &system.m_x[0];
	btVec3d lo = x[0], hi = x[0];
	double mass = 0;
	int numFree = 0;
	for (int i = 0; i < n; i++)
	{
		for (int c = 0; c < 3; c++)
		{
			if (x[i][c] < lo[c]) lo[c] = x[i][c];
			if (x[i][c] > hi[c]) hi[c] = x[i][c];
		}
		if (!system.m_fixed[i])
		{
			mass += system.m_mass[i];
			numFree++;
		}
	}
	m_sceneDiagonal = (hi - lo).length();
	if (!(m_sceneDiagonal > 0.0)) m_sceneDiagonal = 1.0;
	m_dhatActive = m_dhat > 0.0 ? m_dhat : m_dhatRelative * m_sceneDiagonal;
	m_bodyGroup.resize(0);
	m_bodyGroup.resize(system.m_bodies.size(), 0);
	for (int b = 0; b < system.m_bodies.size(); b++)
		for (int i = 0; i < m_groupBodies.size(); i++)
			if (m_groupBodies[i] == system.m_bodies[b].m_psb) m_bodyGroup[b] = m_groupIds[i];
	buildCandidates(x, 0, 0.5 * m_dhatActive);
	{
		double gap = 0;
		countActivePairs(x, &gap);
		m_initialStateValid = gap > 0.0;
		if (!m_initialStateValid && !m_warned)
		{
			printf("btVariationalIPC: the step starts from touching or intersecting surfaces; the barrier cannot recover from that. "
				   "Separate the bodies, or put jointed bodies into one collision group (setCollisionGroup).\n");
			m_warned = true;
		}
	}

	if (m_stiffness > 0.0)
	{
		m_kappa = m_stiffness;
		m_kappaMax = m_adaptiveStiffness ? 100.0 * m_stiffness : m_stiffness;
	}
	else
	{
		// Lower bound of the paper / IPC Toolkit: the barrier curvature at a distance of 1e-8 scene
		// diagonals must reach 1e11 times the average nodal mass. Their incremental potential is h^2
		// times ours, hence the division by h^2.
		const double sHat = m_dhatActive * m_dhatActive;
		double d0 = 1e-8 * m_sceneDiagonal;
		d0 *= d0;
		if (d0 >= sHat) d0 = 0.5 * sHat;
		const double t = d0 - sHat;
		const double d2B = -2.0 * log(d0 / sHat) - 4.0 * t / d0 + t * t / (d0 * d0);
		const double averageMass = numFree > 0 ? mass / numFree : 1.0;
		const double kappaMin = 1e11 * averageMass / (4.0 * d0 * d2B) / (system.m_dt * system.m_dt);
		m_kappaMax = 100.0 * kappaMin;
		m_kappa = 1.0;
		// balance the barrier gradient against the gradient of the rest of the potential
		btVec3dArray gB, gE;
		gB.resize(n);
		gE.resize(n);
		addGradient(x, &gB[0]);
		double bb = 0, be = 0;
		for (int i = 0; i < n; i++)
			if (!system.m_fixed[i]) bb += gB[i].length2();
		double kappa = kappaMin;
		if (bb > 0.0)
		{
			m_muted = true;
			system.gradient(x, &gE[0]);
			m_muted = false;
			for (int i = 0; i < n; i++)
				if (!system.m_fixed[i]) be += gB[i].dot(gE[i]);
			kappa = -be / bb;
		}
		m_kappa = kappa < kappaMin ? kappaMin : (kappa > m_kappaMax ? m_kappaMax : kappa);
	}
	countActivePairs(x, &m_prevMinDistance2);
	m_prevMinDistance2 *= m_prevMinDistance2;
	if (m_friction) buildFrictionPairs(x);
}

double btVariationalIPC::energy(const btVec3d* x) const
{
	if (m_muted) return 0.0;
	double E = 0;
	PairValue value;
	for (int q = 0; q < m_pairs.size(); q++)
		if (evaluatePair(m_pairs[q], x, 0, value)) E += value.m_energy;
	for (int q = 0; q < m_frictionPairs.size(); q++) frictionKernel(m_frictionPairs[q], x, E, 0, 0);
	return E;
}

void btVariationalIPC::addGradient(const btVec3d* x, btVec3d* g) const
{
	if (m_muted) return;
	PairValue value;
	for (int q = 0; q < m_pairs.size(); q++)
	{
		if (!evaluatePair(m_pairs[q], x, 1, value) || value.m_energy >= 1e300) continue;
		for (int j = 0; j < 4; j++) g[m_pairs[q].m_n[j]] += value.m_grad[j];
	}
	for (int q = 0; q < m_frictionPairs.size(); q++)
	{
		const FrictionPair& fp = m_frictionPairs[q];
		double unused = 0;
		btVec3d gvec;
		frictionKernel(fp, x, unused, &gvec, 0);
		for (int j = 0; j < 4; j++)
			if (fp.m_w[j] != 0.0) g[fp.m_n[j]] += gvec * fp.m_w[j];
	}
}

void btVariationalIPC::addHessian(const btVec3d* x, btVariationalHessian& H, btVariationalHessianFilter filter) const
{
	if (m_muted) return;
	PairValue value;
	for (int q = 0; q < m_pairs.size(); q++)
	{
		if (!evaluatePair(m_pairs[q], x, 2, value) || value.m_energy >= 1e300) continue;
		H.m_stencils.expand();
		btVariationalHessianStencil& st = H.m_stencils[H.m_stencils.size() - 1];
		st.m_H.setZero();
		for (int j = 0; j < 4; j++)
		{
			st.m_n[j] = m_pairs[q].m_n[j];
			for (int k = 0; k < 4; k++) st.m_H.addBlock3(j, k, value.m_hess[j][k]);
		}
		btFilterSymmetric(st.m_H, filter);
	}
	for (int q = 0; q < m_frictionPairs.size(); q++)
	{
		const FrictionPair& fp = m_frictionPairs[q];
		double unused = 0;
		btMat3d K;
		frictionKernel(fp, x, unused, 0, &K);
		H.m_stencils.expand();
		btVariationalHessianStencil& st = H.m_stencils[H.m_stencils.size() - 1];
		st.m_H.setZero();
		for (int j = 0; j < 4; j++)
		{
			st.m_n[j] = fp.m_n[j];
			for (int k = 0; k < 4; k++)
				if (fp.m_w[j] != 0.0 && fp.m_w[k] != 0.0) st.m_H.addBlock3(j, k, K * (fp.m_w[j] * fp.m_w[k]));
		}
	}
}

void btVariationalIPC::addVertexForceAndHessian(int node, const btVec3d* x, btVec3d& force, btMat3d& hess) const
{
	if (m_muted) return;
	PairValue value;
	for (int p = m_nodePairStart[node]; p < m_nodePairStart[node + 1]; p++)
	{
		const Pair& pair = m_pairs[m_nodePairs[p] >> 2];
		const int j = m_nodePairs[p] & 3;
		if (distanceLowerBound(pair, x) >= m_dhatActive) continue;  // cannot be active
		if (!evaluatePair(pair, x, 2, value, j) || value.m_energy >= 1e300) continue;
		force -= value.m_grad[j];
		// the diagonal block of a barrier can be indefinite (it is repulsive sideways); keep its positive part
		btMatNd<3> B;
		for (int r = 0; r < 3; r++)
			for (int c = 0; c < 3; c++) B(r, c) = 0.5 * (value.m_hess[j][j].m[r][c] + value.m_hess[j][j].m[c][r]);
		btFilterSymmetric(B, BT_VARIATIONAL_FILTER_CLAMP);
		for (int r = 0; r < 3; r++)
			for (int c = 0; c < 3; c++) hess.m[r][c] += B(r, c);
	}
	if (m_nodeFrictionStart.size() > node + 1)
	{
		for (int p = m_nodeFrictionStart[node]; p < m_nodeFrictionStart[node + 1]; p++)
		{
			const FrictionPair& fp = m_frictionPairs[m_nodeFriction[p] >> 2];
			const double w = fp.m_w[m_nodeFriction[p] & 3];
			double unused = 0;
			btVec3d gvec;
			btMat3d K;
			frictionKernel(fp, x, unused, &gvec, &K);
			force -= gvec * w;
			hess += K * (w * w);
		}
	}
}

double btVariationalIPC::maxStepSize(const btVec3d* x, const btVec3d* p) const
{
	m_reference.resize(0);
	buildCandidates(x, p, 0.5 * m_dhatActive);
	double alpha = 1.0;
	for (int q = 0; q < m_pairs.size(); q++)
	{
		const Pair& pair = m_pairs[q];
		const btVec3d xs[4] = {x[pair.m_n[0]], x[pair.m_n[1]], x[pair.m_n[2]], x[pair.m_n[3]]};
		const btVec3d ps[4] = {p[pair.m_n[0]], p[pair.m_n[1]], p[pair.m_n[2]], p[pair.m_n[3]]};
		alpha = btAdditiveCCD(xs, ps, pair.m_edgeEdge, alpha);
		if (alpha == 0.0) break;
	}
	return alpha;
}

double btVariationalIPC::prepareLocalSteps(const btVec3d* x, double radius, const btVec3d* direction)
{
	if (direction)
	{
		// every node stays on its own segment: swept boxes cover all positions, no trust region needed
		m_reference.resize(0);
		buildCandidates(x, direction, 0.5 * m_dhatActive);
		return radius;
	}
	if (radius < 2.0 * m_dhatActive) radius = 2.0 * m_dhatActive;
	m_localRadius = radius;
	const int n = m_system->numNodes();
	m_reference.resize(n);
	for (int i = 0; i < n; i++) m_reference[i] = x[i];
	buildCandidates(x, 0, 0.5 * m_dhatActive + radius);
	for (int q = 0; q < m_pairs.size(); q++)
	{
		Pair& pair = m_pairs[q];
		const btVec3d xs[4] = {x[pair.m_n[0]], x[pair.m_n[1]], x[pair.m_n[2]], x[pair.m_n[3]]};
		pair.m_refDistance = sqrt(pair.m_edgeEdge ? btEdgeEdgeDistanceSquared(xs) : btPointTriangleDistanceSquared(xs));
	}
	return radius;
}

double btVariationalIPC::vertexMaxStepSize(int node, const btVec3d* x, const btVec3d& dx, bool conservative) const
{
	double alpha = 1.0;
	const double len = dx.length();
	if (len == 0.0) return 1.0;
	if (m_reference.size() > node)
	{
		// stay within the radius the candidate set was built for
		const double room = m_localRadius - (x[node] - m_reference[node]).length();
		if (room <= 0.0) return 0.0;
		if (len > room) alpha = room / len;
	}
	for (int p = m_nodePairStart[node]; p < m_nodePairStart[node + 1]; p++)
	{
		const Pair& pair = m_pairs[m_nodePairs[p] >> 2];
		// far enough that even the whole step keeps the pair outside the zone the step control cares about
		if (distanceLowerBound(pair, x) - len * alpha >= 0.5 * m_dhatActive) continue;
		const btVec3d xs[4] = {x[pair.m_n[0]], x[pair.m_n[1]], x[pair.m_n[2]], x[pair.m_n[3]]};
		btVec3d ps[4];
		ps[m_nodePairs[p] & 3] = dx;
		alpha = conservative ? btAdditiveCCD(xs, ps, pair.m_edgeEdge, alpha, 0.9, 0.5 * m_dhatActive, 0.9) : btAdditiveCCD(xs, ps, pair.m_edgeEdge, alpha);
		if (alpha == 0.0) break;
	}
	return alpha;
}

bool btVariationalIPC::mayHaveContact(const btVec3d* x, const btVec3d* direction) const
{
	m_reference.resize(0);
	buildCandidates(x, direction, 0.5 * m_dhatActive);
	for (int q = 0; q < m_pairs.size(); q++)
	{
		const Pair& pair = m_pairs[q];
		const btVec3d xs[4] = {x[pair.m_n[0]], x[pair.m_n[1]], x[pair.m_n[2]], x[pair.m_n[3]]};
		// how much the two primitives can approach each other: displacements relative to their mean
		btVec3d mean(0, 0, 0);
		for (int j = 0; j < 4; j++) mean += direction[pair.m_n[j]] * 0.25;
		double len[4];
		for (int j = 0; j < 4; j++) len[j] = (direction[pair.m_n[j]] - mean).length();
		const double a = pair.m_edgeEdge ? (len[0] > len[1] ? len[0] : len[1]) : len[0];
		double b = len[2] > len[3] ? len[2] : len[3];
		if (!pair.m_edgeEdge && len[1] > b) b = len[1];
		const double d = sqrt(pair.m_edgeEdge ? btEdgeEdgeDistanceSquared(xs) : btPointTriangleDistanceSquared(xs));
		if (d - a - b < m_dhatActive) return true;
	}
	return false;
}

bool btVariationalIPC::updateStiffness(const btVec3d* x)
{
	if (!m_adaptiveStiffness) return false;
	double minDistance;
	countActivePairs(x, &minDistance);
	const double min2 = minDistance * minDistance;
	double eps = 1e-9 * m_sceneDiagonal;
	eps *= eps;
	bool changed = false;
	if (m_prevMinDistance2 < eps && min2 < eps && min2 < m_prevMinDistance2 && m_kappa < m_kappaMax)
	{
		m_kappa = 2.0 * m_kappa < m_kappaMax ? 2.0 * m_kappa : m_kappaMax;
		changed = true;
	}
	m_prevMinDistance2 = min2;
	return changed;
}

int btVariationalIPC::countActivePairs(const btVec3d* x, double* minDistance) const
{
	int active = 0;
	double min2 = 1e300;
	const double sHat = m_dhatActive * m_dhatActive;
	for (int q = 0; q < m_pairs.size(); q++)
	{
		const Pair& pair = m_pairs[q];
		const btVec3d xs[4] = {x[pair.m_n[0]], x[pair.m_n[1]], x[pair.m_n[2]], x[pair.m_n[3]]};
		const double s = pair.m_edgeEdge ? btEdgeEdgeDistanceSquared(xs) : btPointTriangleDistanceSquared(xs);
		if (s < sHat) active++;
		if (s < min2) min2 = s;
	}
	if (minDistance) *minDistance = min2 < 1e300 ? sqrt(min2) : 1e150;
	return active;
}

bool btVariationalIPC::hasIntersections(const btVec3d* x) const
{
	const btVariationalSystem& sys = *m_system;
	for (int t = 0; t < sys.m_triangles.size(); t++)
	{
		const int* tn = sys.m_triangles[t].m_n;
		const btVec3d a = x[tn[0]], e1 = x[tn[1]] - a, e2 = x[tn[2]] - a;
		for (int e = 0; e < sys.m_edges.size(); e++)
		{
			const int* en = sys.m_edges[e].m_n;
			if (en[0] == tn[0] || en[0] == tn[1] || en[0] == tn[2] || en[1] == tn[0] || en[1] == tn[1] || en[1] == tn[2]) continue;
			if (excluded(sys.m_edges[e].m_body, sys.m_triangles[t].m_body)) continue;  // not a pair this contact model handles
			// segment - triangle intersection (Moeller-Trumbore)
			const btVec3d o = x[en[0]], dir = x[en[1]] - o;
			const btVec3d h = dir.cross(e2);
			const double det = e1.dot(h);
			// an edge (numerically) parallel to the triangle plane does not cross it
			if (fabs(det) <= 1e-10 * e1.length() * e2.length() * dir.length()) continue;
			const btVec3d s = o - a;
			const double u = s.dot(h) / det;
			if (u < 0.0 || u > 1.0) continue;
			const btVec3d q = s.cross(e1);
			const double v = dir.dot(q) / det;
			if (v < 0.0 || u + v > 1.0) continue;
			const double tt = e2.dot(q) / det;
			if (tt >= 0.0 && tt <= 1.0) return true;
		}
	}
	return false;
}
