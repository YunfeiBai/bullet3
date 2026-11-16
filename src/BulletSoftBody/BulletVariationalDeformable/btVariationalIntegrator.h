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

#ifndef BT_VARIATIONAL_INTEGRATOR_H
#define BT_VARIATIONAL_INTEGRATOR_H

#include "btVariationalSystem.h"

///What a solver did in its last time step
struct btVariationalSolveStats
{
	int m_iterations;        // outer iterations (Newton steps, local/global iterations, sweeps, ...)
	int m_linearIterations;  // total conjugate gradient iterations, if any
	int m_lineSearchSteps;   // total backtracking steps
	double m_initialGradient;  // infinity norm of the gradient of the incremental potential at the initial guess
	double m_finalGradient;    // ... and at the returned iterate
	double m_energy;         // incremental potential at the returned iterate
	double m_minStepSize;    // smallest step length taken (1 if no line search was needed)
	bool m_converged;

	btVariationalSolveStats() { reset(); }
	void reset()
	{
		m_iterations = m_linearIterations = m_lineSearchSteps = 0;
		m_initialGradient = m_finalGradient = m_energy = 0.0;
		m_minStepSize = 1.0;
		m_converged = false;
	}
};

///Augmented Lagrangian bookkeeping for the hard attachments of a system, as introduced for vertex
///block descent by "Augmented Vertex Block Descent" (Giles, Diaz, Yuksel 2025). Every hard attachment
///carries a dual variable lambda and a penalty stiffness k; its energy is lambda.C + k/2 |C|^2.
///  warmStart():   lambda <- alpha gamma lambda, k <- gamma k (clamped), and the error present at the
///                 start of the step is only removed at the rate (1 - alpha) per step,
///  updateDuals(): lambda <- lambda + k C,  k <- k + beta |C|.
///Any solver can use this: the constraint error then converges to zero instead of settling at the
///compromise a plain penalty spring finds, also for extreme mass and stiffness ratios.
struct btVariationalAugmentedLagrangian
{
	bool m_enabled;
	double m_alpha;  // fraction of the start-of-step constraint error that is kept (0.99 in the paper)
	double m_beta;   // growth of the penalty stiffness per unit of constraint error and update
	double m_gamma;  // warm start decay of lambda and of the penalty stiffness
	double m_penaltyMin;
	double m_penaltyMax;
	///Optional multiplicative growth of the penalty stiffness for attachments whose error is still above
	///m_growthThreshold at an update. Vertex block descent keeps this at 1 (stiff penalties slow down
	///Gauss-Seidel); solvers with a global linear solve can afford a large penalty and need far fewer
	///dual updates with it.
	double m_growth;
	double m_growthThreshold;  // [m]

	btVariationalAugmentedLagrangian()
		: m_enabled(true), m_alpha(0.95), m_beta(1e5), m_gamma(0.99), m_penaltyMin(1.0), m_penaltyMax(1e9), m_growth(1.0), m_growthThreshold(1e-4)
	{
	}

	///call once at the start of a step (after btVariationalSystem::beginStep)
	void warmStart(btVariationalSystem& system) const;
	///dual update at the iterate x; returns the largest constraint error |C| of the hard attachments
	double updateDuals(btVariationalSystem& system, const btVec3d* x) const;
	///call once at the end of a step
	static void finishStep(btVariationalSystem& system);
	static bool hasHardAttachments(const btVariationalSystem& system);
};

///A time integrator minimizes the incremental potential of btVariationalSystem for one time step.
///On entry system.m_x holds the start-of-step positions x_n and system.m_xTilde the inertial target;
///on exit system.m_x holds x_{n+1}. Velocities are updated by the system afterwards.
class btVariationalIntegrator
{
public:
	btVariationalSolveStats m_stats;
	bool m_verbose;  // print per-iteration diagnostics

	btVariationalIntegrator() : m_verbose(false) {}
	virtual ~btVariationalIntegrator() {}
	virtual const char* getName() const = 0;
	virtual void solve(btVariationalSystem& system) = 0;

	///Initial guess of a step. Kinematic (fixed) nodes are moved to their targets as far as the
	///barrier terms allow; the fraction of the way they covered is returned. If freeNodesToo is set
	///the free nodes start from the inertial target or from x_n, whichever has the lower energy.
	static double advanceToInitialGuess(btVariationalSystem& system, bool freeNodesToo);

	///Moves every free node towards its target position as far as the barrier terms allow, one node
	///at a time. Nodes that are not near any contact take the whole step; without barrier terms this
	///is simply x = target. Unlike a single step length for all nodes, a body in resting contact does
	///not hold back a body in free flight.
	static void moveNodesTowards(btVariationalSystem& system, const btVec3d* target);

	///Backtracking (Armijo) line search along p starting at alphaMax. Returns the accepted step
	///length (0 if no decrease was found) and updates x in place.
	static double lineSearch(const btVariationalSystem& system, btVec3d* x, const btVec3d* p, const btVec3d* g, double alphaMax,
							 double energy0, btVec3d* scratch, int& numSteps, double* energyOut);
};

///Block Jacobi preconditioner: inverse of the 3x3 diagonal blocks
class btVariationalBlockJacobi : public btVariationalPreconditioner
{
	btAlignedObjectArray<btMat3d> m_inverse;

public:
	virtual void build(const btVariationalHessian& H);
	virtual void apply(const btVec3d* r, btVec3d* z) const;
	const btMat3d& inverseBlock(int node) const { return m_inverse[node]; }
};

#endif  // BT_VARIATIONAL_INTEGRATOR_H
