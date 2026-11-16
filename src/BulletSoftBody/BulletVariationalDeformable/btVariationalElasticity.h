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

#ifndef BT_VARIATIONAL_ELASTICITY_H
#define BT_VARIATIONAL_ELASTICITY_H

#include "btVariationalMath.h"

///Eigen-decomposition of the 9x9 energy Hessian d2Psi/dF2 of an isotropic material:
/// d2Psi/dF2 = sum_i m_lambda[i] * vec(m_Q[i]) * vec(m_Q[i])^T, with orthonormal eigen-matrices m_Q.
///Entries 0-2 are the twist modes, 3-5 the flip modes and 6-8 the scaling modes.
struct btVariationalEigenSystem
{
	double m_lambda[9];
	btMat3d m_Q[9];
};

///Hyperelastic material of a tetrahedral element, expressed in the deformation gradient F.
class btVariationalMaterial
{
public:
	virtual ~btVariationalMaterial() {}
	virtual const char* getName() const = 0;

	///strain energy density Psi(F)
	virtual double energyDensity(const btMat3d& F) const = 0;

	///first Piola-Kirchhoff stress dPsi/dF
	virtual btMat3d firstPiola(const btMat3d& F) const = 0;

	///analytic, unfiltered eigen-system of d2Psi/dF2
	virtual void eigenSystem(const btMat3d& F, btVariationalEigenSystem& es) const = 0;

	///Stiffness k of the constant Laplacian approximation (k/2)|F - R|^2 of this material. It weights
	///the global matrix of projective dynamics and the initial Hessian of the quasi-Newton solver.
	virtual double laplacianStiffness() const = 0;

	///Projective dynamics local step. If the energy has the form (k/2)|F - p(F)|^2 this returns true
	///and fills the projection p(F) and k; materials that are not of that form return false.
	virtual bool projectiveDynamicsTarget(const btMat3d& F, btMat3d& target, double& stiffness) const
	{
		(void)F;
		(void)target;
		(void)stiffness;
		return false;
	}

	///Contribution of one element to a single vertex whose shape gradient is g (F = sum_j x_j g_j^T):
	///Pg = dPsi/dF * g and the 3x3 block sum_{b,d} g_b g_d d2Psi/dF_ab dF_cd. Vertex block descent
	///only needs these blocks. The default implementation uses the clamped eigen-system.
	virtual void vertexBlock(const btMat3d& F, const btVec3d& g, btVec3d& Pg, btMat3d& block) const;

	///dense 9x9 Hessian in column-major vec(F) ordering (index = row + 3 * column)
	void hessian(const btMat3d& F, btMat9d& H, btVariationalHessianFilter filter) const;

	///builds the three twist and three flip eigen-matrices shared by all isotropic materials
	static void buildTwistAndFlip(const btMat3d& U, const btMat3d& V, btMat3d* Q);
};

///"Stable Neo-Hookean Flesh Simulation", Smith, de Goes, Kim, ACM Transactions on Graphics 37(2), 2018.
///
/// Psi = mu/2 (I_C - 3) + lambda/2 (J - alpha)^2 - mu/2 log(I_C + 1)        (Eq. 14 of the paper)
///
///with I_C = |F|^2 and J = det F. The energy is finite and rest-stable for inverted and degenerate
///elements and its Hessian has a closed-form eigen-system (Section 4 of the paper), which is what
///makes projected Newton and vertex block descent robust here. The logarithmic term can be switched
///off, which gives the simpler variant the authors use in production (see the HOBAK library).
class btVariationalStableNeoHookean : public btVariationalMaterial
{
	double m_youngsModulus, m_poissonRatio;
	double m_mu, m_lambda, m_alpha;  // reparameterized constants of the energy
	bool m_logTerm;

public:
	btVariationalStableNeoHookean(double youngsModulus, double poissonRatio, bool logTerm = true);
	void setParameters(double youngsModulus, double poissonRatio, bool logTerm);
	double getYoungsModulus() const { return m_youngsModulus; }
	double getPoissonRatio() const { return m_poissonRatio; }

	virtual const char* getName() const { return "StableNeoHookean"; }
	virtual double energyDensity(const btMat3d& F) const;
	virtual btMat3d firstPiola(const btMat3d& F) const;
	virtual void eigenSystem(const btMat3d& F, btVariationalEigenSystem& es) const;
	virtual double laplacianStiffness() const;
	///closed form: the determinant is linear in a single vertex, so the block is positive semi-definite as is
	virtual void vertexBlock(const btMat3d& F, const btVec3d& g, btVec3d& Pg, btMat3d& block) const;
};

///As-rigid-as-possible (co-rotational) material, Psi = mu |F - R|^2 with R the rotation closest to F.
///This is the tetrahedral strain constraint of "Projective Dynamics" (Bouaziz et al. 2014) with
///singular values projected to 1, so the local/global solver minimizes it exactly.
class btVariationalARAP : public btVariationalMaterial
{
	double m_youngsModulus, m_poissonRatio;
	double m_mu;

public:
	btVariationalARAP(double youngsModulus, double poissonRatio);
	void setParameters(double youngsModulus, double poissonRatio);

	virtual const char* getName() const { return "ARAP"; }
	virtual double energyDensity(const btMat3d& F) const;
	virtual btMat3d firstPiola(const btMat3d& F) const;
	virtual void eigenSystem(const btMat3d& F, btVariationalEigenSystem& es) const;
	virtual double laplacianStiffness() const { return 2.0 * m_mu; }
	virtual bool projectiveDynamicsTarget(const btMat3d& F, btMat3d& target, double& stiffness) const;
};

#endif  // BT_VARIATIONAL_ELASTICITY_H
