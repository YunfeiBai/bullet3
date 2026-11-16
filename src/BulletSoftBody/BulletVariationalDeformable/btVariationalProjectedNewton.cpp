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

#include "btVariationalProjectedNewton.h"
#include <stdio.h>
#include "LinearMath/btQuickprof.h"

//
// helpers shared by all integrators
//

double btVariationalIntegrator::advanceToInitialGuess(btVariationalSystem& system, bool freeNodesToo)
{
	const int n = system.numNodes();
	btVec3dArray p;
	p.resize(n);
	// kinematic nodes first: move them to their targets as far as contact allows
	double covered = 1.0;
	bool any = false;
	for (int i = 0; i < n; i++)
	{
		p[i] = system.m_fixed[i] ? system.m_fixedTarget[i] - system.m_x[i] : btVec3d(0, 0, 0);
		if (p[i].length2() > 0.0) any = true;
	}
	if (any)
	{
		covered = system.maxStepSize(&system.m_x[0], &p[0]);
		for (int i = 0; i < n; i++) system.m_x[i] += p[i] * covered;
	}
	if (!freeNodesToo) return covered;
	// free nodes: the inertial target is the better start in free flight, the previous positions are
	// better when stiff forces dominate (large time steps, strong kinematic motion). Take the one
	// with the lower energy.
	const double energyStay = system.energy(&system.m_x[0]);
	for (int i = 0; i < n; i++) p[i] = system.m_x[i];
	moveNodesTowards(system, &system.m_xTilde[0]);
	if (!(system.energy(&system.m_x[0]) < energyStay))
		for (int i = 0; i < n; i++) system.m_x[i] = p[i];
	return covered;
}

void btVariationalIntegrator::moveNodesTowards(btVariationalSystem& system, const btVec3d* target)
{
	const int n = system.numNodes();
	const bool barrier = system.m_contact != 0 || (system.m_colliderModel == btVariationalSystem::COLLIDER_BARRIER && system.m_evaluateColliders);
	if (!barrier)
	{
		for (int i = 0; i < n; i++)
			if (!system.m_fixed[i]) system.m_x[i] = target[i];
		return;
	}
	btVec3dArray direction;
	direction.resize(n);
	double reach = 0;
	for (int i = 0; i < n; i++)
	{
		if (system.m_fixed[i]) continue;
		direction[i] = target[i] - system.m_x[i];
		const double d = direction[i].length();
		if (d > reach) reach = d;
	}
	if (reach == 0.0) return;
	if (system.m_contact) system.m_contact->prepareLocalSteps(&system.m_x[0], reach, &direction[0]);
	for (int i = 0; i < n; i++)
	{
		if (system.m_fixed[i]) continue;
		const btVec3d dx = direction[i];
		if (dx.length2() == 0.0) continue;
		const double alpha = system.vertexMaxStepSize(i, &system.m_x[0], dx, true);
		system.m_x[i] += dx * alpha;
	}
}

double btVariationalIntegrator::lineSearch(const btVariationalSystem& system, btVec3d* x, const btVec3d* p, const btVec3d* g, double alphaMax,
										   double energy0, btVec3d* scratch, int& numSteps, double* energyOut)
{
	const int n = system.numNodes();
	double slope = 0;
	for (int i = 0; i < n; i++) slope += g[i].dot(p[i]);
	double alpha = alphaMax;
	for (int ls = 0; ls < 40; ls++)
	{
		for (int i = 0; i < n; i++) scratch[i] = x[i] + p[i] * alpha;
		const double E = system.energy(scratch);
		if (E <= energy0 + 1e-4 * alpha * slope)
		{
			for (int i = 0; i < n; i++) x[i] = scratch[i];
			if (energyOut) *energyOut = E;
			return alpha;
		}
		alpha *= 0.5;
		numSteps++;
	}
	if (energyOut) *energyOut = energy0;
	return 0.0;
}

void btVariationalAugmentedLagrangian::warmStart(btVariationalSystem& system) const
{
	for (int a = 0; a < system.m_attachments.size(); a++)
	{
		btVariationalAttachment& at = system.m_attachments[a];
		if (!at.m_hard || !m_enabled)
		{
			at.m_c0.setZero();
			if (!at.m_hard) at.m_lambda.setZero();
			continue;
		}
		at.m_lambda *= m_alpha * m_gamma;
		at.m_penalty *= m_gamma;
		if (at.m_penalty < m_penaltyMin) at.m_penalty = m_penaltyMin;
		if (at.m_penalty > m_penaltyMax) at.m_penalty = m_penaltyMax;
		// Stabilization: C = C(x) - alpha C(x_n). The error is measured against the previous target, so
		// a moving target is followed exactly and only old error is bled off slowly.
		const btVec3d c0 = system.m_xPrev[at.m_nodeA] - (at.m_nodeB >= 0 ? system.m_xPrev[at.m_nodeB] : at.m_prevTarget);
		at.m_c0 = c0 * m_alpha;
	}
}

double btVariationalAugmentedLagrangian::updateDuals(btVariationalSystem& system, const btVec3d* x) const
{
	double worst = 0.0;
	for (int a = 0; a < system.m_attachments.size(); a++)
	{
		btVariationalAttachment& at = system.m_attachments[a];
		if (!at.m_hard) continue;
		const btVec3d raw = x[at.m_nodeA] - (at.m_nodeB >= 0 ? x[at.m_nodeB] : at.m_target);
		if (raw.length() > worst) worst = raw.length();
		if (!m_enabled) continue;
		const btVec3d C = raw - at.m_c0;
		at.m_lambda += C * at.m_penalty;
		at.m_penalty += m_beta * C.length();
		if (m_growth > 1.0 && C.length() > m_growthThreshold) at.m_penalty *= m_growth;
		if (at.m_penalty > m_penaltyMax) at.m_penalty = m_penaltyMax;
	}
	return worst;
}

void btVariationalAugmentedLagrangian::finishStep(btVariationalSystem& system)
{
	for (int a = 0; a < system.m_attachments.size(); a++) system.m_attachments[a].m_prevTarget = system.m_attachments[a].m_target;
}

bool btVariationalAugmentedLagrangian::hasHardAttachments(const btVariationalSystem& system)
{
	for (int a = 0; a < system.m_attachments.size(); a++)
		if (system.m_attachments[a].m_hard) return true;
	return false;
}

void btVariationalBlockJacobi::build(const btVariationalHessian& H)
{
	const int n = H.m_nodeBlocks.size();
	m_inverse.resize(n);
	H.diagonalBlocks(&m_inverse[0]);
	for (int i = 0; i < n; i++)
	{
		const btMat3d& D = m_inverse[i];
		const double det = D.det();
		const double scale = fabs(D.m[0][0]) + fabs(D.m[1][1]) + fabs(D.m[2][2]);
		if (det > 1e-18 * scale * scale * scale && scale > 0.0)
			m_inverse[i] = D.inverse();
		else
			m_inverse[i] = btMat3d::identity() * (scale > 0.0 ? 3.0 / scale : 1.0);
	}
}

void btVariationalBlockJacobi::apply(const btVec3d* r, btVec3d* z) const
{
	for (int i = 0; i < m_inverse.size(); i++) z[i] = m_inverse[i] * r[i];
}

//
// Projected Newton
//

void btVariationalProjectedNewton::solve(btVariationalSystem& system)
{
	const int n = system.numNodes();
	m_stats.reset();
	m_numProjectedSolves = m_numExactSolves = 0;
	if (n == 0) return;
	m_gradient.resize(n);
	m_direction.resize(n);
	m_scratch.resize(n);
	btVec3d* x = &system.m_x[0];
	const char* fixed = &system.m_fixed[0];

	btClock clock;
	double tGuess = 0, tGradient = 0, tHessian = 0, tLinear = 0, tCcd = 0, tSearch = 0;  // seconds, reported when verbose
#define BT_NEWTON_LAP(acc)                              \
	{                                                   \
		acc += clock.getTimeMicroseconds() * 1e-6;      \
		clock.reset();                                  \
	}
	double covered = m_skipInitialGuess ? advanceToInitialGuess(system, false) : advanceToInitialGuess(system, true);
	if (m_skipInitialGuess && system.m_contact)
	{
		// make the contact candidates valid for the current positions
		for (int i = 0; i < n; i++) m_direction[i].setZero();
		system.maxStepSize(x, &m_direction[0]);
	}
	double E = system.energy(x);
	BT_NEWTON_LAP(tGuess);
	int skipExact = 0;  // iterations for which the exact Hessian is not tried after it failed
	for (int iter = 0; iter < m_maxIterations; iter++)
	{
		if (covered < 1.0)
		{
			// kinematic nodes that were blocked by contact keep trying to reach their targets
			covered = advanceToInitialGuess(system, false);
			E = system.energy(x);
		}
		// adaptive barrier stiffness of IPC: stiffen when contacts keep closing in
		if (system.m_contact && system.m_contact->updateStiffness(x)) E = system.energy(x);
		system.gradient(x, &m_gradient[0]);
		const double gNorm = system.gradientNorm(&m_gradient[0]);
		if (iter == 0) m_stats.m_initialGradient = gNorm;
		m_stats.m_finalGradient = gNorm;

		// First try the exact Hessian (if enabled), fall back to the eigenvalue filtered one whenever it
		// does not produce a usable descent step.
		bool exact = m_projectOnDemand && system.m_hessianFilter != BT_VARIATIONAL_FILTER_NONE && skipExact == 0;
		if (skipExact > 0) skipExact--;
		bool done = false;
		BT_NEWTON_LAP(tGradient);
		for (;;)
		{
			system.hessian(x, m_hessian, exact ? BT_VARIATIONAL_FILTER_NONE : system.m_hessianFilter);
			BT_NEWTON_LAP(tHessian);
			for (int i = 0; i < n; i++) m_scratch[i] = -m_gradient[i];
			btVariationalPreconditioner* pre = m_preconditioner ? m_preconditioner : &m_blockJacobi;
			pre->build(m_hessian);
			bool indefinite = false;
			m_stats.m_linearIterations += btVariationalConjugateGradient(m_hessian, &m_scratch[0], &m_direction[0], fixed, n, pre, m_linearTolerance, m_maxLinearIterations, &indefinite);
			BT_NEWTON_LAP(tLinear);
			double step = 0, slope = 0;
			int stepNode = -1;
			for (int i = 0; i < n; i++)
			{
				if (fixed[i])
				{
					m_direction[i].setZero();
					continue;
				}
				slope += m_gradient[i].dot(m_direction[i]);
				const double a = m_direction[i].maxAbs();
				if (a > step)
				{
					step = a;
					stepNode = i;
				}
			}
			if (exact && (indefinite || !(slope < 0.0)))
			{
				// the exact Hessian is indefinite here; it usually stays so for a few iterations
				exact = false;
				skipExact = 3;
				continue;
			}
			if (step / system.m_dt < m_tolerance && covered >= 1.0)
			{
				m_stats.m_converged = true;
				done = true;
				break;
			}
			const double alphaMax = system.maxStepSize(x, &m_direction[0]);
			BT_NEWTON_LAP(tCcd);
			const double alpha = lineSearch(system, x, &m_direction[0], &m_gradient[0], alphaMax, E, &m_scratch[0], m_stats.m_lineSearchSteps, &E);
			BT_NEWTON_LAP(tSearch);
			if (alpha == 0.0 && exact)
			{
				exact = false;
				continue;
			}
			if (exact)
				m_numExactSolves++;
			else
				m_numProjectedSolves++;
			m_stats.m_iterations++;
			if (alpha < m_stats.m_minStepSize) m_stats.m_minStepSize = alpha;
			if (m_verbose)
				printf("  newton %3d (%s): |g| %.3e  step %.3e (node %d of body %d)  alphaMax %.3e  alpha %.3e  E %.9e  cg %d\n", iter, exact ? "exact" : "filtered", gNorm, step,
					   stepNode, stepNode >= 0 ? system.m_nodeBody[stepNode] : -1, alphaMax, alpha, E, m_stats.m_linearIterations);
			if (alpha == 0.0) done = true;  // no further decrease possible at working precision
			break;
		}
		if (done) break;
	}
#undef BT_NEWTON_LAP
	if (m_verbose)
		printf("  newton time [ms]: initial guess %.2f, gradient %.2f, hessian %.2f, linear solve %.2f, ccd %.2f, line search %.2f\n",
			   tGuess * 1e3, tGradient * 1e3, tHessian * 1e3, tLinear * 1e3, tCcd * 1e3, tSearch * 1e3);
	m_stats.m_energy = E;
	if (m_stats.m_converged)
	{
		system.gradient(x, &m_gradient[0]);
		m_stats.m_finalGradient = system.gradientNorm(&m_gradient[0]);
	}
}
