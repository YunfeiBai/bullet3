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

#ifndef BT_VARIATIONAL_LAPLACIAN_H
#define BT_VARIATIONAL_LAPLACIAN_H

#include "btVariationalSystem.h"
#include "btVariationalSparse.h"
#include "btVariationalIntegrator.h"

///The constant, prefactorized system matrix at the heart of projective dynamics:
///
///   A = M / h^2 + sum_elements k_e V_e G_e^T G_e + springs + attachments
///
///Every element is replaced by its "as rigid as possible" quadratic (k/2)|F - R|^2, so the matrix does
///not depend on the deformation and has the same n x n pattern for the x, y and z coordinates. It is
///factorized once (sparse Cholesky) and reused until the time step, the mesh, the stiffnesses, the
///attachments or the set of fixed nodes change.
///
///"Fast Simulation of Mass-Spring Systems" (Liu et al. 2013) and "Projective Dynamics" (Bouaziz et
///al. 2014) solve with this matrix in every global step; "Quasi-Newton Methods for Real-Time
///Simulation of Hyperelastic Materials" (Liu, Bouaziz, Kavan 2017) uses it as the initial Hessian of
///L-BFGS. It can also precondition conjugate gradients.
class btVariationalLaplacian : public btVariationalPreconditioner
{
	btVariationalSparseMatrix m_full;  // including the rows of fixed nodes
	btVariationalSparseCholesky m_cholesky;
	btAlignedObjectArray<char> m_fixed;
	// entries A(i, j) with i free and j fixed: they move the fixed node positions to the right hand side
	btAlignedObjectArray<int> m_couplingFree, m_couplingFixed;
	btAlignedObjectArray<double> m_couplingValue;
	double m_signature;
	bool m_valid;
	int m_numFactorizations;

	double computeSignature(const btVariationalSystem& system) const;

public:
	btVariationalLaplacian() : m_signature(0), m_valid(false), m_numFactorizations(0) {}

	///(re)builds and factorizes the matrix if anything it depends on has changed; returns false if
	///the factorization failed
	bool update(const btVariationalSystem& system);

	///Solves A x = rhs on the free nodes. Fixed nodes are passed through (x_i = rhs_i), so for a
	///position solve the caller stores their target there and calls subtractFixedCoupling first; for
	///a direction solve the fixed entries of rhs are simply zero.
	void solve(const btVec3d* rhs, btVec3d* x) const { m_cholesky.solve(rhs, x); }

	///rhs_i -= sum_{j fixed} A(i,j) x_j for the free nodes i
	void subtractFixedCoupling(const btVec3d* x, btVec3d* rhs) const;

	///y = A x with the full matrix (rows of fixed nodes included)
	void multiply(const btVec3d* x, btVec3d* y) const { m_full.multiply(x, y); }

	virtual void apply(const btVec3d* r, btVec3d* z) const { m_cholesky.solve(r, z); }

	///diagonal entry of the full matrix
	double diagonal(int node) const
	{
		for (int p = m_full.m_rowStart[node]; p < m_full.m_rowStart[node + 1]; p++)
			if (m_full.m_cols[p] == node) return m_full.m_vals[p];
		return 0.0;
	}
	bool isFixed(int node) const { return m_fixed[node] != 0; }
	int numFactorizations() const { return m_numFactorizations; }
	int factorNonZeros() const { return m_cholesky.nonZeros(); }
};

///Two-level preconditioner for the Newton systems of the incremental potential:
///
///      M^-1 = D^-1 + P A^-1 P
///
///D holds the 3x3 diagonal blocks of the current Hessian (block Jacobi) and A is the prefactorized
///Laplacian of projective dynamics. The Laplacian solve removes the smooth error components that
///make block Jacobi slow on stiff or finely meshed bodies, at the price of one pair of triangular
///solves and without ever being refactorized. The mask P switches the global correction off for nodes
///whose actual stiffness is far above what the Laplacian knows about (active contact barriers, hard
///attachments); those are left to the block Jacobi part, which sees the true diagonal.
///
///Measured (see the module README): the number of conjugate gradient iterations roughly halves on
///contact-free solids, but with the block sparse matrix-vector product the triangular solves cost
///more than the products they save, so the solve as a whole came out slower than with block Jacobi
///alone on meshes of up to a few thousand nodes. It is an option, not a default.
class btVariationalLaplacianPreconditioner : public btVariationalPreconditioner
{
	const btVariationalLaplacian* m_laplacian;
	btVariationalBlockJacobi m_jacobi;
	btAlignedObjectArray<char> m_masked;
	mutable btVec3dArray m_tmp, m_sol;
	int m_numMasked;

public:
	double m_maskRatio;  // a node is masked when its Hessian diagonal exceeds this multiple of the Laplacian diagonal
	bool m_useLaplacian;  // false: plain block Jacobi

	btVariationalLaplacianPreconditioner() : m_laplacian(0), m_numMasked(0), m_maskRatio(4.0), m_useLaplacian(true) {}
	void setLaplacian(const btVariationalLaplacian* laplacian) { m_laplacian = laplacian; }
	int numMasked() const { return m_numMasked; }
	virtual void build(const btVariationalHessian& H);
	virtual void apply(const btVec3d* r, btVec3d* z) const;
};

#endif  // BT_VARIATIONAL_LAPLACIAN_H
