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

#ifndef BT_VARIATIONAL_PROJECTED_NEWTON_H
#define BT_VARIATIONAL_PROJECTED_NEWTON_H

#include "btVariationalIntegrator.h"

///Projected Newton with line search on the incremental potential.
///
///Every element Hessian is replaced by a positive semi-definite one before assembly, using the
///analytic eigen-system of the material ("Stable Neo-Hookean Flesh Simulation", Smith et al. 2018)
///and either clamping negative eigenvalues to zero or taking their absolute value ("Stabler
///Neo-Hookean Simulation", Chen et al. 2024; select with btVariationalSystem::m_hessianFilter).
///The search direction is therefore always a descent direction, and the backtracking line search,
///started at the largest intersection-free step when barrier contact is active, makes every
///iteration decrease the energy. This is also the solver "Incremental Potential Contact"
///(Li et al. 2020) is built on.
class btVariationalProjectedNewton : public btVariationalIntegrator
{
	btVariationalHessian m_hessian;
	btVariationalBlockJacobi m_blockJacobi;
	btVec3dArray m_gradient, m_direction, m_scratch;

public:
	int m_maxIterations;
	double m_tolerance;  // [m/s] stop when the Newton step, divided by the time step, is below this
	int m_maxLinearIterations;
	double m_linearTolerance;                             // relative residual of the conjugate gradient solve
	const btVariationalPreconditioner* m_preconditioner;  // optional override of the block Jacobi preconditioner
	///If set, an iteration first tries the exact Hessian and only falls back to the filtered one when
	///conjugate gradients meets negative curvature. Exact steps converge quadratically, while the
	///filtered Hessian only gives linear convergence wherever the filter is active.
	bool m_projectOnDemand;
	int m_numProjectedSolves, m_numExactSolves;  // statistics of the last solve

	btVariationalProjectedNewton()
		: m_maxIterations(50), m_tolerance(1e-3), m_maxLinearIterations(2000), m_linearTolerance(1e-2), m_preconditioner(0), m_projectOnDemand(true), m_numProjectedSolves(0), m_numExactSolves(0)
	{
	}

	virtual const char* getName() const { return "ProjectedNewton"; }
	virtual void solve(btVariationalSystem& system);
};

#endif  // BT_VARIATIONAL_PROJECTED_NEWTON_H
