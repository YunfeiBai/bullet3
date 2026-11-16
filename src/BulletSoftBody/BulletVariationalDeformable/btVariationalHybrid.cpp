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

#include "btVariationalHybrid.h"

btVariationalHybrid::btVariationalHybrid()
	: m_tolerance(1e-2),
	  m_maxNewtonIterations(40),
	  m_maxQuasiNewtonIterations(150),
	  m_quasiNewtonToleranceScale(0.25),
	  m_useQuasiNewton(true),
	  m_laplacianMinNodes(600),
	  m_warmStartSweeps(0),
	  m_maxConstraintUpdates(3),
	  m_constraintTolerance(1e-4),
	  m_lastPath(PATH_NONE),
	  m_maxConstraintError(0.0)
{
	m_newton.m_linearTolerance = 0.05;
	// Let the penalties grow quickly while a joint is still open, but only up to a moderate value:
	// very stiff two-node penalties ruin the conditioning of the linear solves, and the dual
	// variables carry the constraint force anyway.
	m_al.m_growth = 4.0;
	m_al.m_penaltyMax = 1e7;
	m_vbd.m_al.m_enabled = false;  // the duals are updated here, not inside the warm start sweeps
}

const char* btVariationalHybrid::lastPathName() const
{
	switch (m_lastPath)
	{
		case PATH_QUASI_NEWTON:
			return "quasi-Newton";
		case PATH_NEWTON:
			return "Newton";
		case PATH_NEWTON_LAPLACIAN:
			return "Newton + Laplacian";
		case PATH_QUASI_NEWTON_THEN_NEWTON:
			return "quasi-Newton, then Newton";
		default:
			return "-";
	}
}

void btVariationalHybrid::solve(btVariationalSystem& system)
{
	const int n = system.numNodes();
	m_stats.reset();
	m_lastPath = PATH_NONE;
	m_maxConstraintError = 0.0;
	if (n == 0) return;
	btVec3d* x = &system.m_x[0];
	m_newton.m_verbose = m_quasiNewton.m_verbose = m_verbose;
	m_newton.m_tolerance = m_tolerance;
	m_newton.m_maxIterations = m_maxNewtonIterations;

	const bool hard = btVariationalAugmentedLagrangian::hasHardAttachments(system);
	m_al.m_growthThreshold = m_constraintTolerance;
	m_al.warmStart(system);

	// Contact-free and free of hard constraints: the Laplacian is a good model of the whole Hessian,
	// so L-BFGS on it gets there without assembling anything.
	const bool contact = system.mayHaveContact(x);
	bool continueFromIterate = false;
	if (m_useQuasiNewton && !contact && !hard)
	{
		m_quasiNewton.m_tolerance = m_tolerance * m_quasiNewtonToleranceScale;
		m_quasiNewton.m_maxIterations = m_maxQuasiNewtonIterations;
		m_quasiNewton.solve(system);
		m_stats = m_quasiNewton.m_stats;
		m_lastPath = PATH_QUASI_NEWTON;
		if (m_stats.m_converged)
		{
			btVariationalAugmentedLagrangian::finishStep(system);
			return;
		}
		continueFromIterate = true;  // out of budget (or contact appeared): Newton finishes the step
	}
	else if (m_warmStartSweeps > 0 && !contact)
	{
		m_vbd.m_iterations = m_warmStartSweeps;
		m_vbd.solve(system);
		continueFromIterate = true;
	}

	const bool useLaplacian = n >= m_laplacianMinNodes && m_laplacian.update(system);
	if (useLaplacian)
	{
		m_preconditioner.setLaplacian(&m_laplacian);
		m_newton.m_preconditioner = &m_preconditioner;
	}
	else
		m_newton.m_preconditioner = 0;

	const btVariationalSolveStats before = m_stats;
	for (int update = 0;; update++)
	{
		m_newton.m_skipInitialGuess = continueFromIterate;
		m_newton.solve(system);
		m_newton.m_skipInitialGuess = false;
		if (update == 0)
		{
			m_stats = m_newton.m_stats;
			m_stats.m_iterations += before.m_iterations;
			m_stats.m_lineSearchSteps += before.m_lineSearchSteps;
			if (m_lastPath == PATH_QUASI_NEWTON) m_stats.m_initialGradient = before.m_initialGradient;
		}
		else
		{
			m_stats.m_iterations += m_newton.m_stats.m_iterations;
			m_stats.m_linearIterations += m_newton.m_stats.m_linearIterations;
			m_stats.m_lineSearchSteps += m_newton.m_stats.m_lineSearchSteps;
			m_stats.m_finalGradient = m_newton.m_stats.m_finalGradient;
			m_stats.m_energy = m_newton.m_stats.m_energy;
			m_stats.m_converged = m_newton.m_stats.m_converged;
		}
		if (!hard) break;
		// augmented Lagrangian: move the constraint force into the dual variable and, if the joints
		// are still open, solve again from the current iterate (a short solve)
		m_maxConstraintError = m_al.updateDuals(system, x);
		if (m_maxConstraintError <= m_constraintTolerance || update + 1 >= m_maxConstraintUpdates) break;
		continueFromIterate = true;
	}
	m_lastPath = m_lastPath == PATH_QUASI_NEWTON ? PATH_QUASI_NEWTON_THEN_NEWTON : (useLaplacian ? PATH_NEWTON_LAPLACIAN : PATH_NEWTON);
	btVariationalAugmentedLagrangian::finishStep(system);
}
