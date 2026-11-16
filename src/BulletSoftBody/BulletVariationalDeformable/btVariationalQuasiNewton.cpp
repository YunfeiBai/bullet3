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

#include "btVariationalQuasiNewton.h"
#include <stdio.h>

void btVariationalQuasiNewton::solve(btVariationalSystem& system)
{
	const int n = system.numNodes();
	m_stats.reset();
	if (n == 0) return;
	btVec3d* x = &system.m_x[0];
	const char* fixed = &system.m_fixed[0];
	m_gradient.resize(n);
	m_newGradient.resize(n);
	m_direction.resize(n);
	m_scratch.resize(n);
	const int m = m_historySize > 0 ? m_historySize : 0;
	if (m_s.size() != m * n || !m_keepHistory)
	{
		m_s.resize(m * n);
		m_y.resize(m * n);
		m_rho.resize(m);
		m_alpha.resize(m);
		m_historyCount = m_historyHead = 0;
	}
	if (!m_laplacian.update(system))
	{
		for (int i = 0; i < n; i++) x[i] = system.m_xTilde[i];
		return;
	}

	advanceToInitialGuess(system, true);
	double E = system.energy(x);
	system.gradient(x, &m_gradient[0]);
	m_stats.m_initialGradient = m_stats.m_finalGradient = system.gradientNorm(&m_gradient[0]);

	for (int iter = 0; iter < m_maxIterations; iter++)
	{
		// two-loop recursion; the initial Hessian is the prefactorized Laplacian
		btVec3d* q = &m_direction[0];
		for (int i = 0; i < n; i++) q[i] = -m_gradient[i];
		for (int h = 0; h < m_historyCount; h++)
		{
			const int slot = (m_historyHead - 1 - h + 2 * m) % m;  // newest first
			const btVec3d* s = &m_s[slot * n];
			const btVec3d* y = &m_y[slot * n];
			double sq = 0;
			for (int i = 0; i < n; i++) sq += s[i].dot(q[i]);
			m_alpha[slot] = m_rho[slot] * sq;
			for (int i = 0; i < n; i++) q[i] -= y[i] * m_alpha[slot];
		}
		m_laplacian.solve(q, &m_scratch[0]);
		btVec3d* p = &m_direction[0];
		for (int i = 0; i < n; i++) p[i] = fixed[i] ? btVec3d(0, 0, 0) : m_scratch[i];
		for (int h = m_historyCount - 1; h >= 0; h--)
		{
			const int slot = (m_historyHead - 1 - h + 2 * m) % m;  // oldest first
			const btVec3d* s = &m_s[slot * n];
			const btVec3d* y = &m_y[slot * n];
			double yr = 0;
			for (int i = 0; i < n; i++) yr += y[i].dot(p[i]);
			const double beta = m_rho[slot] * yr;
			for (int i = 0; i < n; i++) p[i] += s[i] * (m_alpha[slot] - beta);
		}

		double step = 0, slope = 0;
		for (int i = 0; i < n; i++)
		{
			slope += m_gradient[i].dot(p[i]);
			const double a = p[i].maxAbs();
			if (a > step) step = a;
		}
		if (!(slope < 0.0))
		{
			// should not happen with positive curvature pairs; restart from the plain Laplacian step
			m_historyCount = m_historyHead = 0;
			for (int i = 0; i < n; i++) m_scratch[i] = -m_gradient[i];
			m_laplacian.solve(&m_scratch[0], p);
			step = 0;
			slope = 0;
			for (int i = 0; i < n; i++)
			{
				if (fixed[i]) p[i].setZero();
				slope += m_gradient[i].dot(p[i]);
				const double a = p[i].maxAbs();
				if (a > step) step = a;
			}
		}
		if (step / system.m_dt < m_tolerance)
		{
			m_stats.m_converged = true;
			break;
		}
		const double alphaMax = system.maxStepSize(x, p);
		for (int i = 0; i < n; i++) m_scratch[i] = x[i];  // remember the iterate
		const double alpha = lineSearch(system, x, p, &m_gradient[0], alphaMax, E, &m_newGradient[0], m_stats.m_lineSearchSteps, &E);
		if (alpha < m_stats.m_minStepSize) m_stats.m_minStepSize = alpha;
		m_stats.m_iterations++;
		if (alpha == 0.0) break;
		system.gradient(x, &m_newGradient[0]);
		m_stats.m_finalGradient = system.gradientNorm(&m_newGradient[0]);
		if (m_verbose) printf("  l-bfgs %3d: |g| %.3e  step %.3e  alpha %.3e  E %.9e  pairs %d\n", iter, m_stats.m_finalGradient, step, alpha, E, m_historyCount);

		// curvature pair s = x_new - x_old, y = g_new - g_old; only kept if it has positive curvature
		if (m > 0)
		{
			btVec3d* s = &m_s[m_historyHead * n];
			btVec3d* y = &m_y[m_historyHead * n];
			double sy = 0, yy = 0;
			for (int i = 0; i < n; i++)
			{
				s[i] = x[i] - m_scratch[i];
				y[i] = m_newGradient[i] - m_gradient[i];
				sy += s[i].dot(y[i]);
				yy += y[i].length2();
			}
			if (sy > 1e-12 * yy && sy > 0.0)
			{
				m_rho[m_historyHead] = 1.0 / sy;
				m_historyHead = (m_historyHead + 1) % m;
				if (m_historyCount < m) m_historyCount++;
			}
		}
		for (int i = 0; i < n; i++) m_gradient[i] = m_newGradient[i];
	}
	m_stats.m_energy = E;
}
