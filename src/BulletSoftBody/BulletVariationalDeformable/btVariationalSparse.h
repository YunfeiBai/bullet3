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

#ifndef BT_VARIATIONAL_SPARSE_H
#define BT_VARIATIONAL_SPARSE_H

#include "btVariationalMath.h"

///Symmetric sparse matrix in compressed row storage (the full pattern is stored, columns sorted).
class btVariationalSparseMatrix
{
public:
	int m_n;
	btAlignedObjectArray<int> m_rowStart;
	btAlignedObjectArray<int> m_cols;
	btAlignedObjectArray<double> m_vals;

	btVariationalSparseMatrix() : m_n(0) {}

	///Builds the matrix from (row, col, value) triplets; duplicate entries are summed.
	///The triplets must describe a symmetric matrix, i.e. contain both (i,j) and (j,i).
	void setFromTriplets(int n, const btAlignedObjectArray<int>& rows, const btAlignedObjectArray<int>& cols, const btAlignedObjectArray<double>& vals);

	///y = A * x, applied to each of the three components
	void multiply(const btVec3d* x, btVec3d* y) const;
	int nonZeros() const { return m_cols.size(); }
};

///Sparse Cholesky factorization A = L L^T of a symmetric positive definite matrix.
///
///The fill-reducing ordering is a nested dissection built from breadth-first level structures of the
///matrix graph, the numeric factorization is the classic up-looking algorithm driven by the
///elimination tree. Projective dynamics and the quasi-Newton solver factorize their constant system
///matrix once with this class and afterwards only pay for the triangular solves.
class btVariationalSparseCholesky
{
	int m_n;
	btAlignedObjectArray<int> m_perm;     // m_perm[new] = old
	btAlignedObjectArray<int> m_invPerm;  // m_invPerm[old] = new
	btAlignedObjectArray<int> m_colStart;
	btAlignedObjectArray<int> m_rowIdx;
	btAlignedObjectArray<double> m_Lx;
	mutable btVec3dArray m_work;

	void computeOrdering(const btVariationalSparseMatrix& A);

public:
	btVariationalSparseCholesky() : m_n(0) {}

	///returns false if the matrix is not (numerically) positive definite
	bool factorize(const btVariationalSparseMatrix& A);

	///solves A x = b for the three components at once; b and x may alias
	void solve(const btVec3d* b, btVec3d* x) const;

	int size() const { return m_n; }
	int nonZeros() const { return m_rowIdx.size(); }
};

#endif  // BT_VARIATIONAL_SPARSE_H
