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

#ifndef BT_VARIATIONAL_PROJECTIVE_DYNAMICS_H
#define BT_VARIATIONAL_PROJECTIVE_DYNAMICS_H

#include "btVariationalIntegrator.h"
#include "btVariationalLaplacian.h"

///Local/global solver of
///  "Fast Simulation of Mass-Spring Systems", Liu, Bargteil, O'Brien, Kavan, SIGGRAPH Asia 2013, and
///  "Projective Dynamics: Fusing Constraint Projections for Fast Simulation", Bouaziz, Martin, Liu,
///   Kavan, Pauly, SIGGRAPH 2014.
///
///Every elastic term is written as a distance to a constraint set, (k/2)|A x - p|^2:
///  - springs: p is the spring direction scaled to the rest length (Liu et al. 2013),
///  - tetrahedra: p is the rotation closest to the deformation gradient (the strain constraint of
///    Bouaziz et al. 2014, i.e. the btVariationalARAP material),
///  - attachments: p is the target point.
///One iteration alternates
///  local step:  with the positions fixed, project every element onto its constraint set (independent
///               per element), and
///  global step: with the projections fixed, minimize over the positions. This is a linear solve with
///               the constant matrix M/h^2 + sum k A^T A, which is factorized once and reused.
///Each iteration decreases the energy, without any line search, and costs one pair of sparse
///triangular solves.
///
///Materials that are not of projective form (e.g. Stable Neo-Hookean) are still accepted: their
///elements take a step with the same constant matrix along the true energy gradient. That step has no
///descent guarantee; use btVariationalQuasiNewton for such materials.
///
///Contact with rigid colliders is handled the way the projective dynamics paper suggests, as a
///projection after every global step: penetrating nodes are moved back to the surface and friction
///is applied to their tangential motion (Coulomb, position based).
class btVariationalProjectiveDynamics : public btVariationalIntegrator
{
	btVariationalLaplacian m_laplacian;
	btVec3dArray m_rhs, m_next;

	void localStep(const btVariationalSystem& system, const btVec3d* x, btVec3d* rhs) const;

public:
	int m_iterations;    // local/global iterations per time step
	double m_tolerance;  // [m/s] stop early when an iteration moves no node faster than this (0: never)

	btVariationalProjectiveDynamics() : m_iterations(10), m_tolerance(0.0) {}

	virtual const char* getName() const { return "ProjectiveDynamics"; }
	virtual void solve(btVariationalSystem& system);
	const btVariationalLaplacian& getLaplacian() const { return m_laplacian; }

	///moves nodes out of the rigid colliders and applies position-based Coulomb friction
	static void projectColliders(const btVariationalSystem& system, btVec3d* x);
};

#endif  // BT_VARIATIONAL_PROJECTIVE_DYNAMICS_H
