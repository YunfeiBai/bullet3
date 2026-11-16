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

#ifndef BT_VARIATIONAL_DISTANCE_H
#define BT_VARIATIONAL_DISTANCE_H

///Unsigned distances between the primitives of triangle meshes, with exact first and second
///derivatives, and additive continuous collision detection built on them.
///
///A contact pair is either point-triangle, with nodes (p, t0, t1, t2), or edge-edge, with nodes
///(ea0, ea1, eb0, eb1). Depending on where the closest points lie, its distance is the distance of a
///point-point, point-line, point-plane or line-line configuration of a subset of the four nodes
///(Section 4 of "Incremental Potential Contact", Li et al. 2020). The classification follows the IPC
///Toolkit; the derivatives are written out in closed form here.

#include "btVariationalMath.h"

enum btVariationalClosestKind
{
	BT_CLOSEST_POINT_POINT = 0,
	BT_CLOSEST_POINT_LINE,
	BT_CLOSEST_POINT_PLANE,
	BT_CLOSEST_LINE_LINE
};

///Which sub-primitives of a pair are closest. m_slot lists the pair nodes (0..3) in the order the
///reduced distance function takes them: point-point (a, b), point-line (p, e0, e1),
///point-plane (p, t0, t1, t2), line-line (ea0, ea1, eb0, eb1).
struct btVariationalClosest
{
	int m_kind;
	int m_slot[4];
	void set(int kind, int s0, int s1, int s2 = -1, int s3 = -1)
	{
		m_kind = kind;
		m_slot[0] = s0;
		m_slot[1] = s1;
		m_slot[2] = s2;
		m_slot[3] = s3;
	}
	int count() const { return m_kind == BT_CLOSEST_POINT_POINT ? 2 : (m_kind == BT_CLOSEST_POINT_LINE ? 3 : 4); }
};

///Squared distance of a pair and its derivatives with respect to the four pair nodes
struct btVariationalPairDistance
{
	double m_distanceSquared;
	btVec3d m_grad[4];
	btMat3d m_hess[4][4];
	///closest point difference as a combination of the nodes: c = sum_j m_weights[j] x_j, |c| = distance
	double m_weights[4];

	void clear(int order)
	{
		m_distanceSquared = 0;
		for (int j = 0; j < 4; j++)
		{
			m_weights[j] = 0;
			if (order >= 1) m_grad[j].setZero();
			if (order >= 2)
				for (int k = 0; k < 4; k++) m_hess[j][k].setZero();
		}
	}
};

///closest features of a point p and a triangle (t0, t1, t2); pair node order (p, t0, t1, t2)
inline btVariationalClosest btPointTriangleClosest(const btVec3d& p, const btVec3d& t0, const btVec3d& t1, const btVec3d& t2)
{
	btVariationalClosest c;
	const btVec3d normal = (t1 - t0).cross(t2 - t0);
	const btVec3d corner[3] = {t0, t1, t2};
	double along[3];
	for (int k = 0; k < 3; k++)
	{
		const btVec3d e = corner[(k + 1) % 3] - corner[k];
		const btVec3d d = p - corner[k];
		const btVec3d perp = e.cross(normal);  // in the triangle plane, pointing away from the triangle
		const double a = e.length2(), c2 = perp.length2();
		along[k] = a > 0.0 ? e.dot(d) / a : 0.0;
		const double off = c2 > 0.0 ? perp.dot(d) / c2 : 0.0;
		if (along[k] > 0.0 && along[k] < 1.0 && off >= 0.0)
		{
			c.set(BT_CLOSEST_POINT_LINE, 0, 1 + k, 1 + (k + 1) % 3);
			return c;
		}
	}
	if (along[0] <= 0.0 && along[2] >= 1.0)
		c.set(BT_CLOSEST_POINT_POINT, 0, 1);
	else if (along[1] <= 0.0 && along[0] >= 1.0)
		c.set(BT_CLOSEST_POINT_POINT, 0, 2);
	else if (along[2] <= 0.0 && along[1] >= 1.0)
		c.set(BT_CLOSEST_POINT_POINT, 0, 3);
	else
		c.set(BT_CLOSEST_POINT_PLANE, 0, 1, 2, 3);
	return c;
}

///closest features of two edges; pair node order (ea0, ea1, eb0, eb1)
inline btVariationalClosest btEdgeEdgeClosest(const btVec3d& ea0, const btVec3d& ea1, const btVec3d& eb0, const btVec3d& eb1)
{
	btVariationalClosest res;
	const btVec3d u = ea1 - ea0, v = eb1 - eb0, w = ea0 - eb0;
	const double a = u.length2(), b = u.dot(v), c = v.length2(), d = u.dot(w), e = v.dot(w);
	const double D = a * c - b * b;
	if (a == 0.0 && c == 0.0)
	{
		res.set(BT_CLOSEST_POINT_POINT, 0, 2);
		return res;
	}
	if (a == 0.0)
	{
		res.set(BT_CLOSEST_POINT_LINE, 0, 2, 3);
		return res;
	}
	if (c == 0.0)
	{
		res.set(BT_CLOSEST_POINT_LINE, 2, 0, 1);
		return res;
	}
	// which end (0, 1) or interior (2) of each edge is closest
	int ea = 2, eb = 2;
	const bool parallel = u.cross(v).length2() < 2.5e-16 * a * c;
	if (parallel)
	{
		// position of the ends of edge b along edge a
		const double alpha = (eb0 - ea0).dot(u) / a, beta = (eb1 - ea0).dot(u) / a;
		if (alpha < 0.0)
		{
			ea = (0.0 <= beta && beta <= 1.0) ? 2 : 0;
			eb = (beta <= alpha) ? 0 : (beta <= 1.0 ? 1 : 2);
		}
		else if (alpha > 1.0)
		{
			ea = (0.0 <= beta && beta <= 1.0) ? 2 : 1;
			eb = (beta >= alpha) ? 0 : (0.0 <= beta ? 1 : 2);
		}
		else
		{
			ea = 2;
			eb = 0;
		}
	}
	else
	{
		const double sN = b * e - c * d;
		double tN, tD;
		if (sN <= 0.0)
		{
			tN = e;
			tD = c;
			ea = 0;
		}
		else if (sN >= D)
		{
			tN = e + b;
			tD = c;
			ea = 1;
		}
		else
		{
			tN = a * e - b * d;
			tD = D;
		}
		if (tN <= 0.0)
		{
			eb = 0;
			ea = (-d <= 0.0) ? 0 : (-d >= a ? 1 : 2);
		}
		else if (tN >= tD)
		{
			eb = 1;
			ea = ((-d + b) <= 0.0) ? 0 : ((-d + b) >= a ? 1 : 2);
		}
	}
	if (ea < 2 && eb < 2)
		res.set(BT_CLOSEST_POINT_POINT, ea, 2 + eb);
	else if (ea < 2)
		res.set(BT_CLOSEST_POINT_LINE, ea, 2, 3);
	else if (eb < 2)
		res.set(BT_CLOSEST_POINT_LINE, 2 + eb, 0, 1);
	else
		res.set(BT_CLOSEST_LINE_LINE, 0, 1, 2, 3);
	return res;
}

///Evaluates the squared distance of the closest features. order = 0: value and weights only,
///1: plus gradient, 2: plus Hessian. x holds the four pair nodes. If onlySlot >= 0 the Hessian is
///only filled in for the diagonal block of that pair node (what a per-vertex solver needs).
inline void btEvaluatePairDistance(const btVec3d* x, const btVariationalClosest& c, int order, btVariationalPairDistance& out, int onlySlot = -1)
{
	out.clear(order);
	const int* s = c.m_slot;
	if (c.m_kind == BT_CLOSEST_POINT_POINT)
	{
		const btVec3d d = x[s[0]] - x[s[1]];
		out.m_distanceSquared = d.length2();
		out.m_weights[s[0]] = 1.0;
		out.m_weights[s[1]] = -1.0;
		if (order >= 1)
		{
			out.m_grad[s[0]] = d * 2.0;
			out.m_grad[s[1]] = d * -2.0;
		}
		if (order >= 2)
		{
			const btMat3d I2 = btMat3d::identity() * 2.0;
			out.m_hess[s[0]][s[0]] = I2;
			out.m_hess[s[1]][s[1]] = I2;
			out.m_hess[s[0]][s[1]] = I2 * -1.0;
			out.m_hess[s[1]][s[0]] = I2 * -1.0;
		}
		return;
	}
	if (c.m_kind == BT_CLOSEST_POINT_LINE)
	{
		// f(w, e) = |w|^2 - (w.e)^2 / |e|^2 with w = p - e0, e = e1 - e0
		const btVec3d w = x[s[0]] - x[s[1]], e = x[s[2]] - x[s[1]];
		const double a = w.dot(e), b = e.length2();
		const double t = a / b;
		out.m_distanceSquared = w.length2() - a * t;
		if (out.m_distanceSquared < 0.0) out.m_distanceSquared = 0.0;
		out.m_weights[s[0]] = 1.0;
		out.m_weights[s[1]] = -(1.0 - t);
		out.m_weights[s[2]] = -t;
		if (order >= 1)
		{
			const btVec3d gw = (w - e * t) * 2.0;
			const btVec3d ge = (e * t - w) * (2.0 * t);
			out.m_grad[s[0]] = gw;
			out.m_grad[s[1]] = -(gw + ge);
			out.m_grad[s[2]] = ge;
			if (order >= 2)
			{
				btMat3d Hww = btMat3d::identity() * 2.0;
				Hww.addScaledOuter(-2.0 / b, e, e);
				btMat3d Hwe = btMat3d::identity() * (-2.0 * t);
				Hwe.addScaledOuter(-2.0 / b, e, w);
				Hwe.addScaledOuter(4.0 * t / b, e, e);
				btMat3d Hee = btMat3d::identity() * (2.0 * t * t);
				Hee.addScaledOuter(-2.0 / b, w, w);
				Hee.addScaledOuter(4.0 * t / b, w, e);
				Hee.addScaledOuter(4.0 * t / b, e, w);
				Hee.addScaledOuter(-8.0 * t * t / b, e, e);
				const btMat3d Hew = Hwe.transpose();
				// p = w + e0, e0 = -(w + e) + ..., e1 = e + e0
				out.m_hess[s[0]][s[0]] = Hww;
				out.m_hess[s[0]][s[2]] = Hwe;
				out.m_hess[s[2]][s[0]] = Hew;
				out.m_hess[s[2]][s[2]] = Hee;
				out.m_hess[s[0]][s[1]] = (Hww + Hwe) * -1.0;
				out.m_hess[s[1]][s[0]] = (Hww + Hew) * -1.0;
				out.m_hess[s[2]][s[1]] = (Hew + Hee) * -1.0;
				out.m_hess[s[1]][s[2]] = (Hwe + Hee) * -1.0;
				out.m_hess[s[1]][s[1]] = Hww + Hwe + Hew + Hee;
			}
		}
		return;
	}
	// point-plane and line-line share f(w, e1, e2) = det[w e1 e2]^2 / |e1 x e2|^2
	btVec3d w, e1, e2;
	double cw[4] = {0, 0, 0, 0}, c1[4] = {0, 0, 0, 0}, c2[4] = {0, 0, 0, 0};  // w, e1, e2 as combinations of the reduced nodes
	if (c.m_kind == BT_CLOSEST_POINT_PLANE)
	{
		w = x[s[0]] - x[s[1]];
		e1 = x[s[2]] - x[s[1]];
		e2 = x[s[3]] - x[s[1]];
		cw[0] = 1;
		cw[1] = -1;
		c1[2] = 1;
		c1[1] = -1;
		c2[3] = 1;
		c2[1] = -1;
	}
	else
	{
		w = x[s[0]] - x[s[2]];
		e1 = x[s[1]] - x[s[0]];
		e2 = x[s[3]] - x[s[2]];
		cw[0] = 1;
		cw[2] = -1;
		c1[1] = 1;
		c1[0] = -1;
		c2[3] = 1;
		c2[2] = -1;
	}
	const btVec3d N = e1.cross(e2);
	const double q = N.length2(), sv = w.dot(N);
	out.m_distanceSquared = sv * sv / q;
	// closest point difference: the component of w along the normal, written with the node weights
	if (c.m_kind == BT_CLOSEST_POINT_PLANE)
	{
		// barycentric coordinates of the projection of p in the triangle
		const double a11 = e1.length2(), a12 = e1.dot(e2), a22 = e2.length2(), b1 = w.dot(e1), b2 = w.dot(e2);
		const double det = a11 * a22 - a12 * a12;
		const double beta1 = (a22 * b1 - a12 * b2) / det, beta2 = (a11 * b2 - a12 * b1) / det;
		out.m_weights[s[0]] = 1.0;
		out.m_weights[s[1]] = -(1.0 - beta1 - beta2);
		out.m_weights[s[2]] = -beta1;
		out.m_weights[s[3]] = -beta2;
	}
	else
	{
		// closest point parameters on the two lines: minimize |w + sa e1 - sb e2|
		const double a11 = e1.length2(), a12 = e1.dot(e2), a22 = e2.length2(), b1 = w.dot(e1), b2 = w.dot(e2);
		const double det = a11 * a22 - a12 * a12;
		const double sa = (a12 * b2 - a22 * b1) / det, sb = (a11 * b2 - a12 * b1) / det;
		out.m_weights[s[0]] = 1.0 - sa;
		out.m_weights[s[1]] = sa;
		out.m_weights[s[2]] = -(1.0 - sb);
		out.m_weights[s[3]] = -sb;
	}
	if (order < 1) return;
	const double e1e2 = e1.dot(e2), n1 = e1.length2(), n2 = e2.length2();
	// gradients of s = w.(e1 x e2) and q = |e1 x e2|^2 with respect to (w, e1, e2)
	const btVec3d ds[3] = {N, e2.cross(w), w.cross(e1)};
	const btVec3d dq[3] = {btVec3d(0, 0, 0), (e1 * n2 - e2 * e1e2) * 2.0, (e2 * n1 - e1 * e1e2) * 2.0};
	const double f1 = 2.0 * sv / q, f2 = sv * sv / (q * q);
	btVec3d g[3];
	for (int a = 0; a < 3; a++) g[a] = ds[a] * f1 - dq[a] * f2;
	const double* coef[3] = {cw, c1, c2};
	for (int j = 0; j < 4; j++)
		for (int a = 0; a < 3; a++)
			if (coef[a][j] != 0.0) out.m_grad[s[j]] += g[a] * coef[a][j];
	if (order < 2) return;
	// second derivatives
	btMat3d d2s[3][3], d2q[3][3];
	d2s[0][1] = btMat3d::skew(e2) * -1.0;
	d2s[1][0] = btMat3d::skew(e2);
	d2s[0][2] = btMat3d::skew(e1);
	d2s[2][0] = btMat3d::skew(e1) * -1.0;
	d2s[1][2] = btMat3d::skew(w) * -1.0;
	d2s[2][1] = btMat3d::skew(w);
	d2q[1][1] = btMat3d::identity() * (2.0 * n2);
	d2q[1][1].addScaledOuter(-2.0, e2, e2);
	d2q[2][2] = btMat3d::identity() * (2.0 * n1);
	d2q[2][2].addScaledOuter(-2.0, e1, e1);
	d2q[1][2] = btMat3d::identity() * (-2.0 * e1e2);
	d2q[1][2].addScaledOuter(4.0, e1, e2);
	d2q[1][2].addScaledOuter(-2.0, e2, e1);
	d2q[2][1] = d2q[1][2].transpose();
	const double h1 = 2.0 / q, h3 = 2.0 * sv / (q * q), h4 = 2.0 * sv * sv / (q * q * q);
	btMat3d H[3][3];
	int only = -1;  // position of the requested slot among the reduced nodes
	for (int j = 0; j < 4; j++)
		if (onlySlot >= 0 && s[j] == onlySlot) only = j;
	if (onlySlot >= 0 && only < 0) return;  // the requested node is not part of the closest features
	for (int a = 0; a < 3; a++)
		for (int b = 0; b < 3; b++)
		{
			if (only >= 0 && (coef[a][only] == 0.0 || coef[b][only] == 0.0)) continue;
			btMat3d M = d2s[a][b] * f1 - d2q[a][b] * f2;
			M.addScaledOuter(h1, ds[a], ds[b]);
			M.addScaledOuter(-h3, ds[a], dq[b]);
			M.addScaledOuter(-h3, dq[a], ds[b]);
			M.addScaledOuter(h4, dq[a], dq[b]);
			H[a][b] = M;
		}
	for (int j = 0; j < 4; j++)
		for (int k = 0; k < 4; k++)
		{
			if (onlySlot >= 0 && (s[j] != onlySlot || s[k] != onlySlot)) continue;
			btMat3d B;
			for (int a = 0; a < 3; a++)
				for (int b = 0; b < 3; b++)
				{
					const double cc = coef[a][j] * coef[b][k];
					if (cc != 0.0) B += H[a][b] * cc;
				}
			out.m_hess[s[j]][s[k]] = B;
		}
}

inline double btPointTriangleDistanceSquared(const btVec3d* x)
{
	btVariationalPairDistance d;
	btEvaluatePairDistance(x, btPointTriangleClosest(x[0], x[1], x[2], x[3]), 0, d);
	return d.m_distanceSquared;
}

inline double btEdgeEdgeDistanceSquared(const btVec3d* x)
{
	btVariationalPairDistance d;
	btEvaluatePairDistance(x, btEdgeEdgeClosest(x[0], x[1], x[2], x[3]), 0, d);
	return d.m_distanceSquared;
}

///Additive continuous collision detection (ACCD, "Codimensional Incremental Potential Contact",
///Li, Kaufman, Jiang 2021) for one pair moving linearly from x to x + p.
///Returns a step t in [0, tMax] that is free of contact. Every sub-step advances by at most the
///fraction "rescaling" of the current distance divided by the largest relative displacement, which
///cannot close the gap; the iteration stops once the distance has dropped below (1 - rescaling) of
///its initial value, so a gap always remains. Only distance evaluations are needed, which makes the
///test robust to degenerate configurations.
///
///If minDistanceCap > 0 the pair additionally keeps a separation of
///xi = min(minDistanceCap, keepFraction * initial distance): the step stops before the distance drops
///below xi. Solvers that move one vertex at a time use this to enter a contact barrier gradually.
inline double btAdditiveCCD(const btVec3d* x0, const btVec3d* p0, bool edgeEdge, double tMax, double rescaling = 0.9,
							double minDistanceCap = 0.0, double keepFraction = 0.0)
{
	btVec3d x[4], p[4];
	btVec3d mean(0, 0, 0);
	for (int j = 0; j < 4; j++) mean += p0[j] * 0.25;
	double len[4];
	for (int j = 0; j < 4; j++)
	{
		x[j] = x0[j];
		p[j] = p0[j] - mean;  // a common translation does not change distances
		len[j] = p[j].length();
	}
	const double maxDisp = edgeEdge ? (len[0] > len[1] ? len[0] : len[1]) + (len[2] > len[3] ? len[2] : len[3])
									: len[0] + (len[1] > len[2] ? (len[1] > len[3] ? len[1] : len[3]) : (len[2] > len[3] ? len[2] : len[3]));
	if (maxDisp == 0.0) return tMax;
	double d = sqrt(edgeEdge ? btEdgeEdgeDistanceSquared(x) : btPointTriangleDistanceSquared(x));
	if (!(d > 0.0)) return 0.0;
	double xi = 0.0;
	if (minDistanceCap > 0.0) xi = minDistanceCap < keepFraction * d ? minDistanceCap : keepFraction * d;
	const double gap = (1.0 - rescaling) * (d - xi);
	double toi = 0.0;
	for (int it = 0; it < 100000; it++)
	{
		const double advance = rescaling * (d - xi) / maxDisp;
		for (int j = 0; j < 4; j++) x[j] += p[j] * advance;
		d = sqrt(edgeEdge ? btEdgeEdgeDistanceSquared(x) : btPointTriangleDistanceSquared(x));
		if (toi > 0.0 && d - xi < gap) break;
		toi += advance;
		if (toi > tMax) return tMax;
	}
	return toi;
}

#endif  // BT_VARIATIONAL_DISTANCE_H
