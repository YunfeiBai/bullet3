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

#ifndef BT_VARIATIONAL_HYBRID_H
#define BT_VARIATIONAL_HYBRID_H

#include "btVariationalIntegrator.h"
#include "btVariationalProjectedNewton.h"
#include "btVariationalQuasiNewton.h"
#include "btVariationalVertexBlockDescent.h"
#include "btVariationalLaplacian.h"

///A solver that combines the methods of this module, each where it measurably helps.
///
///  - Backbone: projected Newton on the incremental potential with the eigenvalue filtered Hessians
///    of the Stable Neo-Hookean material (Smith et al. 2018; clamp or the absolute value filter of
///    Chen et al. 2024) and a line search. It is the most robust of the solvers and by far the
///    fastest one as soon as stiff contact is involved.
///  - Contact: the barrier, continuous collision detection and lagged friction of Incremental
///    Potential Contact (Li et al. 2020), if the system has a contact model. Every step is
///    intersection free.
///  - Contact-free steps: when no node can reach a collider or another surface during the step, the
///    quasi-Newton method of Liu, Bouaziz, Kavan 2017 (L-BFGS on the prefactorized matrix of
///    projective dynamics; Liu et al. 2013, Bouaziz et al. 2014) solves the step without assembling
///    any Hessian, about twice as fast as Newton. If it does not converge within its budget, Newton
///    takes over from where it stopped.
///  - Optionally the same prefactorized matrix, together with block Jacobi, preconditions the
///    conjugate gradient solves inside Newton (btVariationalLaplacianPreconditioner). This is off by
///    default: it halves the number of conjugate gradient iterations on a contact-free bar, but a
///    pair of sparse triangular solves costs more than the products it saves, and the step came out
///    15 to 25 percent slower at every size measured (272 to 4617 nodes).
///  - Hard attachments: the augmented Lagrangian of Augmented VBD (Giles et al. 2025) wraps the
///    solve, so joints close exactly instead of stretching like penalty springs.
///  - Optionally a few sweeps of vertex block descent (Chen et al. 2024) warm start a contact-free
///    Newton solve. This is off by default: it saved under ten percent where it helped and cost time
///    where contact was active.
///
///See the module README for the measurements behind these choices.
class btVariationalHybrid : public btVariationalIntegrator
{
public:
	enum Path
	{
		PATH_NONE = 0,
		PATH_QUASI_NEWTON,            // contact-free step solved by L-BFGS
		PATH_NEWTON,                  // Newton with block Jacobi preconditioned CG
		PATH_NEWTON_LAPLACIAN,        // Newton with the Laplacian + block Jacobi preconditioner
		PATH_QUASI_NEWTON_THEN_NEWTON  // L-BFGS ran out of budget and Newton finished the step
	};

	double m_tolerance;                 // [m/s] convergence tolerance of the Newton step
	int m_maxNewtonIterations;
	int m_maxQuasiNewtonIterations;     // budget of the contact-free fast path
	double m_quasiNewtonToleranceScale;  // the fast path uses m_tolerance times this (its step is measured in a stiffer metric)
	bool m_useQuasiNewton;              // enable the contact-free fast path
	bool m_useLaplacianPreconditioner;  // precondition Newton's linear solves with the projective dynamics matrix
	int m_warmStartSweeps;              // vertex block descent sweeps before a contact-free Newton solve
	int m_maxConstraintUpdates;         // augmented Lagrangian updates (re-solves) per step
	double m_constraintTolerance;       // [m] stop updating once the hard attachments are this accurate
	btVariationalAugmentedLagrangian m_al;

	btVariationalHybrid();

	virtual const char* getName() const { return "Hybrid"; }
	virtual void solve(btVariationalSystem& system);

	int lastPath() const { return m_lastPath; }
	const char* lastPathName() const;
	double maxConstraintError() const { return m_maxConstraintError; }
	btVariationalProjectedNewton& newton() { return m_newton; }
	btVariationalQuasiNewton& quasiNewton() { return m_quasiNewton; }

private:
	btVariationalProjectedNewton m_newton;
	btVariationalQuasiNewton m_quasiNewton;
	btVariationalVertexBlockDescent m_vbd;
	btVariationalLaplacian m_laplacian;
	btVariationalLaplacianPreconditioner m_preconditioner;
	int m_lastPath;
	double m_maxConstraintError;
};

#endif  // BT_VARIATIONAL_HYBRID_H
