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

#ifndef BT_VARIATIONAL_MATH_H
#define BT_VARIATIONAL_MATH_H

///Small double-precision linear algebra used by the variational deformable solvers.
///
///The optimization-based integrators in this module (projected Newton, barrier contact, ...) need
///more precision than a single precision btScalar provides, so their state is kept in double precision
///regardless of how Bullet itself is compiled. The types below are deliberately minimal.

#include "LinearMath/btVector3.h"
#include "LinearMath/btAlignedObjectArray.h"
#include <math.h>

struct btVec3d
{
	double m_v[3];

	btVec3d() { m_v[0] = m_v[1] = m_v[2] = 0.0; }
	btVec3d(double x, double y, double z)
	{
		m_v[0] = x;
		m_v[1] = y;
		m_v[2] = z;
	}
	explicit btVec3d(const btVector3& v)
	{
		m_v[0] = v.x();
		m_v[1] = v.y();
		m_v[2] = v.z();
	}
	btVector3 toBt() const { return btVector3(btScalar(m_v[0]), btScalar(m_v[1]), btScalar(m_v[2])); }

	double& operator[](int i) { return m_v[i]; }
	const double& operator[](int i) const { return m_v[i]; }

	btVec3d operator+(const btVec3d& o) const { return btVec3d(m_v[0] + o.m_v[0], m_v[1] + o.m_v[1], m_v[2] + o.m_v[2]); }
	btVec3d operator-(const btVec3d& o) const { return btVec3d(m_v[0] - o.m_v[0], m_v[1] - o.m_v[1], m_v[2] - o.m_v[2]); }
	btVec3d operator-() const { return btVec3d(-m_v[0], -m_v[1], -m_v[2]); }
	btVec3d operator*(double s) const { return btVec3d(m_v[0] * s, m_v[1] * s, m_v[2] * s); }
	btVec3d operator/(double s) const { return btVec3d(m_v[0] / s, m_v[1] / s, m_v[2] / s); }
	btVec3d& operator+=(const btVec3d& o)
	{
		m_v[0] += o.m_v[0];
		m_v[1] += o.m_v[1];
		m_v[2] += o.m_v[2];
		return *this;
	}
	btVec3d& operator-=(const btVec3d& o)
	{
		m_v[0] -= o.m_v[0];
		m_v[1] -= o.m_v[1];
		m_v[2] -= o.m_v[2];
		return *this;
	}
	btVec3d& operator*=(double s)
	{
		m_v[0] *= s;
		m_v[1] *= s;
		m_v[2] *= s;
		return *this;
	}
	double dot(const btVec3d& o) const { return m_v[0] * o.m_v[0] + m_v[1] * o.m_v[1] + m_v[2] * o.m_v[2]; }
	btVec3d cross(const btVec3d& o) const
	{
		return btVec3d(m_v[1] * o.m_v[2] - m_v[2] * o.m_v[1],
					   m_v[2] * o.m_v[0] - m_v[0] * o.m_v[2],
					   m_v[0] * o.m_v[1] - m_v[1] * o.m_v[0]);
	}
	double length2() const { return dot(*this); }
	double length() const { return sqrt(length2()); }
	///returns the zero vector when the length is (numerically) zero
	btVec3d normalized() const
	{
		double l = length();
		return (l > 1e-300) ? (*this) / l : btVec3d(0, 0, 0);
	}
	double maxAbs() const
	{
		double a = fabs(m_v[0]), b = fabs(m_v[1]), c = fabs(m_v[2]);
		return a > b ? (a > c ? a : c) : (b > c ? b : c);
	}
	void setZero() { m_v[0] = m_v[1] = m_v[2] = 0.0; }
};

inline btVec3d operator*(double s, const btVec3d& v) { return v * s; }

typedef btAlignedObjectArray<btVec3d> btVec3dArray;

///3x3 matrix, m[row][col]
struct btMat3d
{
	double m[3][3];

	btMat3d() { setZero(); }
	void setZero()
	{
		for (int r = 0; r < 3; r++)
			for (int c = 0; c < 3; c++) m[r][c] = 0.0;
	}
	static btMat3d identity()
	{
		btMat3d I;
		I.m[0][0] = I.m[1][1] = I.m[2][2] = 1.0;
		return I;
	}
	static btMat3d diag(double a, double b, double c)
	{
		btMat3d D;
		D.m[0][0] = a;
		D.m[1][1] = b;
		D.m[2][2] = c;
		return D;
	}
	static btMat3d fromColumns(const btVec3d& a, const btVec3d& b, const btVec3d& c)
	{
		btMat3d M;
		for (int r = 0; r < 3; r++)
		{
			M.m[r][0] = a[r];
			M.m[r][1] = b[r];
			M.m[r][2] = c[r];
		}
		return M;
	}
	///a * b^T
	static btMat3d outer(const btVec3d& a, const btVec3d& b)
	{
		btMat3d M;
		for (int r = 0; r < 3; r++)
			for (int c = 0; c < 3; c++) M.m[r][c] = a[r] * b[c];
		return M;
	}
	///cross product matrix: skew(v) * w == v x w
	static btMat3d skew(const btVec3d& v)
	{
		btMat3d M;
		M.m[0][1] = -v[2];
		M.m[0][2] = v[1];
		M.m[1][0] = v[2];
		M.m[1][2] = -v[0];
		M.m[2][0] = -v[1];
		M.m[2][1] = v[0];
		return M;
	}
	btVec3d col(int c) const { return btVec3d(m[0][c], m[1][c], m[2][c]); }
	btVec3d row(int r) const { return btVec3d(m[r][0], m[r][1], m[r][2]); }
	void setCol(int c, const btVec3d& v)
	{
		m[0][c] = v[0];
		m[1][c] = v[1];
		m[2][c] = v[2];
	}
	btMat3d operator+(const btMat3d& o) const
	{
		btMat3d R;
		for (int r = 0; r < 3; r++)
			for (int c = 0; c < 3; c++) R.m[r][c] = m[r][c] + o.m[r][c];
		return R;
	}
	btMat3d operator-(const btMat3d& o) const
	{
		btMat3d R;
		for (int r = 0; r < 3; r++)
			for (int c = 0; c < 3; c++) R.m[r][c] = m[r][c] - o.m[r][c];
		return R;
	}
	btMat3d operator*(double s) const
	{
		btMat3d R;
		for (int r = 0; r < 3; r++)
			for (int c = 0; c < 3; c++) R.m[r][c] = m[r][c] * s;
		return R;
	}
	btMat3d operator*(const btMat3d& o) const
	{
		btMat3d R;
		for (int r = 0; r < 3; r++)
			for (int c = 0; c < 3; c++)
				R.m[r][c] = m[r][0] * o.m[0][c] + m[r][1] * o.m[1][c] + m[r][2] * o.m[2][c];
		return R;
	}
	btVec3d operator*(const btVec3d& v) const
	{
		return btVec3d(m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
					   m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
					   m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2]);
	}
	///transpose(M) * v
	btVec3d mulT(const btVec3d& v) const
	{
		return btVec3d(m[0][0] * v[0] + m[1][0] * v[1] + m[2][0] * v[2],
					   m[0][1] * v[0] + m[1][1] * v[1] + m[2][1] * v[2],
					   m[0][2] * v[0] + m[1][2] * v[1] + m[2][2] * v[2]);
	}
	btMat3d& operator+=(const btMat3d& o)
	{
		for (int r = 0; r < 3; r++)
			for (int c = 0; c < 3; c++) m[r][c] += o.m[r][c];
		return *this;
	}
	btMat3d& operator-=(const btMat3d& o)
	{
		for (int r = 0; r < 3; r++)
			for (int c = 0; c < 3; c++) m[r][c] -= o.m[r][c];
		return *this;
	}
	btMat3d& operator*=(double s)
	{
		for (int r = 0; r < 3; r++)
			for (int c = 0; c < 3; c++) m[r][c] *= s;
		return *this;
	}
	///adds s * a * b^T
	void addScaledOuter(double s, const btVec3d& a, const btVec3d& b)
	{
		for (int r = 0; r < 3; r++)
		{
			double sa = s * a[r];
			m[r][0] += sa * b[0];
			m[r][1] += sa * b[1];
			m[r][2] += sa * b[2];
		}
	}
	void addScaledIdentity(double s)
	{
		m[0][0] += s;
		m[1][1] += s;
		m[2][2] += s;
	}
	btMat3d transpose() const
	{
		btMat3d R;
		for (int r = 0; r < 3; r++)
			for (int c = 0; c < 3; c++) R.m[r][c] = m[c][r];
		return R;
	}
	double det() const
	{
		return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) -
			   m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
			   m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
	}
	///cofactor matrix, i.e. the derivative of det(M) with respect to M
	btMat3d cofactor() const
	{
		btVec3d c0 = col(0), c1 = col(1), c2 = col(2);
		return fromColumns(c1.cross(c2), c2.cross(c0), c0.cross(c1));
	}
	///inverse via the adjugate; the caller is responsible for det() != 0
	btMat3d inverse() const
	{
		return cofactor().transpose() * (1.0 / det());
	}
	double trace() const { return m[0][0] + m[1][1] + m[2][2]; }
	///Frobenius inner product
	double ddot(const btMat3d& o) const
	{
		double s = 0;
		for (int r = 0; r < 3; r++)
			for (int c = 0; c < 3; c++) s += m[r][c] * o.m[r][c];
		return s;
	}
	double frobeniusNorm2() const { return ddot(*this); }
};

///Small dense fixed-size vector and (row-major) square matrix
template <int N>
struct btVecNd
{
	double v[N];
	btVecNd() { setZero(); }
	void setZero()
	{
		for (int i = 0; i < N; i++) v[i] = 0.0;
	}
	double& operator[](int i) { return v[i]; }
	const double& operator[](int i) const { return v[i]; }
	double dot(const btVecNd<N>& o) const
	{
		double s = 0;
		for (int i = 0; i < N; i++) s += v[i] * o.v[i];
		return s;
	}
};

template <int N>
struct btMatNd
{
	double m[N * N];
	btMatNd() { setZero(); }
	void setZero()
	{
		for (int i = 0; i < N * N; i++) m[i] = 0.0;
	}
	double& operator()(int r, int c) { return m[r * N + c]; }
	const double& operator()(int r, int c) const { return m[r * N + c]; }
	///adds s * a * b^T
	void addScaledOuter(double s, const btVecNd<N>& a, const btVecNd<N>& b)
	{
		for (int r = 0; r < N; r++)
		{
			double sa = s * a.v[r];
			if (sa == 0.0) continue;
			double* row = m + r * N;
			for (int c = 0; c < N; c++) row[c] += sa * b.v[c];
		}
	}
	btVecNd<N> operator*(const btVecNd<N>& x) const
	{
		btVecNd<N> y;
		for (int r = 0; r < N; r++)
		{
			const double* row = m + r * N;
			double s = 0;
			for (int c = 0; c < N; c++) s += row[c] * x.v[c];
			y.v[r] = s;
		}
		return y;
	}
	btMatNd<N>& operator+=(const btMatNd<N>& o)
	{
		for (int i = 0; i < N * N; i++) m[i] += o.m[i];
		return *this;
	}
	btMatNd<N>& operator*=(double s)
	{
		for (int i = 0; i < N * N; i++) m[i] *= s;
		return *this;
	}
	///3x3 sub-block access (block row br, block column bc)
	void addBlock3(int br, int bc, const btMat3d& B)
	{
		for (int r = 0; r < 3; r++)
			for (int c = 0; c < 3; c++) (*this)(3 * br + r, 3 * bc + c) += B.m[r][c];
	}
	btMat3d block3(int br, int bc) const
	{
		btMat3d B;
		for (int r = 0; r < 3; r++)
			for (int c = 0; c < 3; c++) B.m[r][c] = (*this)(3 * br + r, 3 * bc + c);
		return B;
	}
};

typedef btVecNd<9> btVec9d;
typedef btVecNd<12> btVec12d;
typedef btMatNd<9> btMat9d;
typedef btMatNd<12> btMat12d;

///Cyclic Jacobi eigenvalue iteration for a symmetric NxN matrix.
///On return evals holds the eigenvalues and the columns of V the matching orthonormal eigenvectors.
template <int N>
inline void btSymmetricEigenJacobi(const btMatNd<N>& A, double* evals, btMatNd<N>& V, int maxSweeps = 50)
{
	btMatNd<N> a = A;
	V.setZero();
	for (int i = 0; i < N; i++) V(i, i) = 1.0;
	for (int sweep = 0; sweep < maxSweeps; sweep++)
	{
		double off = 0, dia = 0;
		for (int p = 0; p < N; p++)
		{
			dia += a(p, p) * a(p, p);
			for (int q = p + 1; q < N; q++) off += a(p, q) * a(p, q);
		}
		if (off <= 1e-30 * (dia + off) || off == 0.0) break;
		for (int p = 0; p < N - 1; p++)
		{
			for (int q = p + 1; q < N; q++)
			{
				double apq = a(p, q);
				if (apq == 0.0) continue;
				double theta = (a(q, q) - a(p, p)) / (2.0 * apq);
				double t = (theta >= 0 ? 1.0 : -1.0) / (fabs(theta) + sqrt(theta * theta + 1.0));
				double c = 1.0 / sqrt(t * t + 1.0);
				double s = t * c;
				for (int k = 0; k < N; k++)
				{
					double akp = a(k, p), akq = a(k, q);
					a(k, p) = c * akp - s * akq;
					a(k, q) = s * akp + c * akq;
				}
				for (int k = 0; k < N; k++)
				{
					double apk = a(p, k), aqk = a(q, k);
					a(p, k) = c * apk - s * aqk;
					a(q, k) = s * apk + c * aqk;
				}
				for (int k = 0; k < N; k++)
				{
					double vkp = V(k, p), vkq = V(k, q);
					V(k, p) = c * vkp - s * vkq;
					V(k, q) = s * vkp + c * vkq;
				}
			}
		}
	}
	for (int i = 0; i < N; i++) evals[i] = a(i, i);
}

///How negative eigenvalues of a local Hessian are treated before it is used in a Newton-type solve
enum btVariationalHessianFilter
{
	BT_VARIATIONAL_FILTER_NONE = 0,  // keep the exact (possibly indefinite) Hessian
	BT_VARIATIONAL_FILTER_CLAMP,     // max(lambda, 0), the classic projected Newton [Teran et al. 2005]
	BT_VARIATIONAL_FILTER_ABS        // |lambda|, "Stabler Neo-Hookean Simulation" [Chen et al. 2024]
};

inline double btVariationalFilterEigenvalue(double lambda, btVariationalHessianFilter filter)
{
	if (filter == BT_VARIATIONAL_FILTER_CLAMP) return lambda > 0.0 ? lambda : 0.0;
	if (filter == BT_VARIATIONAL_FILTER_ABS) return fabs(lambda);
	return lambda;
}

///Replaces a symmetric matrix by its filtered version (see btVariationalHessianFilter).
///Only the leading n x n part is considered, so the same routine serves 6, 9 and 12 dof stencils.
template <int N>
inline void btFilterSymmetric(btMatNd<N>& H, btVariationalHessianFilter filter)
{
	if (filter == BT_VARIATIONAL_FILTER_NONE) return;
	double evals[N];
	btMatNd<N> V;
	btSymmetricEigenJacobi(H, evals, V);
	bool negative = false;
	for (int i = 0; i < N; i++)
		if (evals[i] < 0.0) negative = true;
	if (!negative) return;
	H.setZero();
	for (int i = 0; i < N; i++)
	{
		double l = btVariationalFilterEigenvalue(evals[i], filter);
		if (l == 0.0) continue;
		for (int r = 0; r < N; r++)
		{
			double lr = l * V(r, i);
			for (int c = 0; c < N; c++) H(r, c) += lr * V(c, i);
		}
	}
}

///Singular value decomposition F = U * diag(S) * V^T by one-sided Jacobi rotations.
///The result is "rotation variant": U and V are proper rotations (det = +1), the singular values are
///sorted by decreasing magnitude and S[2] carries the sign of det(F). This is the convention the
///eigen-analysis of isotropic energies relies on.
inline void btSvd3d(const btMat3d& F, btMat3d& U, btVec3d& S, btMat3d& V)
{
	btVec3d a[3] = {F.col(0), F.col(1), F.col(2)};
	btVec3d v[3] = {btVec3d(1, 0, 0), btVec3d(0, 1, 0), btVec3d(0, 0, 1)};
	for (int sweep = 0; sweep < 60; sweep++)
	{
		bool rotated = false;
		for (int p = 0; p < 2; p++)
		{
			for (int q = p + 1; q < 3; q++)
			{
				double alpha = a[p].length2(), beta = a[q].length2(), gamma = a[p].dot(a[q]);
				if (gamma == 0.0 || gamma * gamma <= 1e-32 * alpha * beta) continue;
				rotated = true;
				double zeta = (beta - alpha) / (2.0 * gamma);
				double t = (zeta >= 0 ? 1.0 : -1.0) / (fabs(zeta) + sqrt(1.0 + zeta * zeta));
				double c = 1.0 / sqrt(1.0 + t * t);
				double s = c * t;
				btVec3d ap = a[p], aq = a[q];
				a[p] = ap * c - aq * s;
				a[q] = ap * s + aq * c;
				btVec3d vp = v[p], vq = v[q];
				v[p] = vp * c - vq * s;
				v[q] = vp * s + vq * c;
			}
		}
		if (!rotated) break;
	}
	double sig[3] = {a[0].length(), a[1].length(), a[2].length()};
	// sort by decreasing singular value
	int order[3] = {0, 1, 2};
	for (int i = 0; i < 2; i++)
		for (int j = 0; j < 2 - i; j++)
			if (sig[order[j]] < sig[order[j + 1]])
			{
				int tmp = order[j];
				order[j] = order[j + 1];
				order[j + 1] = tmp;
			}
	btVec3d u[3], vs[3];
	for (int i = 0; i < 3; i++)
	{
		S[i] = sig[order[i]];
		vs[i] = v[order[i]];
		u[i] = a[order[i]];
	}
	const double tiny = 1e-12 * (S[0] > 1e-300 ? S[0] : 1.0);
	// build an orthonormal U, completing the basis where singular values vanish
	if (S[0] > tiny)
		u[0] = u[0] / S[0];
	else
		u[0] = btVec3d(1, 0, 0);
	if (S[1] > tiny)
	{
		u[1] = u[1] - u[0] * u[0].dot(u[1]);
		u[1] = u[1].normalized();
	}
	if (!(S[1] > tiny) || u[1].length2() < 0.5)
	{
		// any unit vector orthogonal to u[0]
		btVec3d e = fabs(u[0][0]) < 0.9 ? btVec3d(1, 0, 0) : btVec3d(0, 1, 0);
		u[1] = (e - u[0] * u[0].dot(e)).normalized();
	}
	btVec3d u2 = u[0].cross(u[1]);
	// orientation of the third singular vector as it came out of the iteration
	bool flipped = (S[2] > tiny) && (u[2].dot(u2) < 0.0);
	u[2] = u2;  // det(U) == +1 by construction
	// make V a rotation; flipping v2 flips the sign of the third singular triplet
	if (vs[0].cross(vs[1]).dot(vs[2]) < 0.0)
	{
		vs[2] = -vs[2];
		flipped = !flipped;
	}
	if (flipped) S[2] = -S[2];
	U = btMat3d::fromColumns(u[0], u[1], u[2]);
	V = btMat3d::fromColumns(vs[0], vs[1], vs[2]);
}

///Rotation of the polar decomposition F = R * S (closest rotation to F)
inline btMat3d btPolarRotation3d(const btMat3d& F)
{
	btMat3d U, V;
	btVec3d S;
	btSvd3d(F, U, S, V);
	return U * V.transpose();
}

///Solves the symmetric positive definite 3x3 system A x = b. Returns false if A is not usable.
inline bool btSolveSpd3d(const btMat3d& A, const btVec3d& b, btVec3d& x)
{
	double det = A.det();
	double scale = fabs(A.m[0][0]) + fabs(A.m[1][1]) + fabs(A.m[2][2]);
	if (!(det > 1e-18 * scale * scale * scale) || !(scale > 0.0)) return false;
	x = A.cofactor().mulT(b) * (1.0 / det);
	return true;
}

#endif  // BT_VARIATIONAL_MATH_H
