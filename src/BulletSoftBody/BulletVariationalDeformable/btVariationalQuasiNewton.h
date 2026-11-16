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

#ifndef BT_VARIATIONAL_QUASI_NEWTON_H
#define BT_VARIATIONAL_QUASI_NEWTON_H

#include "btVariationalIntegrator.h"
#include "btVariationalLaplacian.h"

///"Quasi-Newton Methods for Real-Time Simulation of Hyperelastic Materials",
///Liu, Bouaziz, Kavan, ACM Transactions on Graphics 36(3), 2017.
///
///The paper observes that a projective dynamics iteration is a quasi-Newton step whose Hessian
///approximation is the constant matrix M/h^2 + L. That reading lifts the restriction to projective
///energies: for any hyperelastic material the same prefactorized matrix serves as the initial Hessian
///of L-BFGS, the curvature pairs of the last few iterations correct it, and a backtracking line search
///guarantees descent. The cost per iteration stays close to that of projective dynamics (one energy
///gradient and one pair of triangular solves) while the convergence is markedly faster, and
///general materials such as Stable Neo-Hookean are supported.
///
///Rigid collider contact enters through its penalty energy, as in the paper's implementation.
class btVariationalQuasiNewton : public btVariationalIntegrator
{
	btVariationalLaplacian m_laplacian;
	btVec3dArray m_gradient, m_newGradient, m_direction, m_scratch;
	btVec3dArray m_s, m_y;  // history, m_historySize blocks of numNodes entries (ring buffer)
	btAlignedObjectArray<double> m_rho, m_alpha;
	int m_historyCount, m_historyHead;

public:
	int m_historySize;    // number of curvature pairs kept (m in the paper, which uses 5)
	int m_maxIterations;
	double m_tolerance;   // [m/s] stop when the step, divided by the time step, is below this
	bool m_keepHistory;   // carry the curvature pairs over to the next time step instead of restarting

	btVariationalQuasiNewton()
		: m_historyCount(0), m_historyHead(0), m_historySize(5), m_maxIterations(20), m_tolerance(1e-3), m_keepHistory(false)
	{
	}

	virtual const char* getName() const { return "QuasiNewton(L-BFGS)"; }
	virtual void solve(btVariationalSystem& system);
	const btVariationalLaplacian& getLaplacian() const { return m_laplacian; }
};

#endif  // BT_VARIATIONAL_QUASI_NEWTON_H
