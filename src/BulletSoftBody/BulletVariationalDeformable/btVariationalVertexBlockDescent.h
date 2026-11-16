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

#ifndef BT_VARIATIONAL_VERTEX_BLOCK_DESCENT_H
#define BT_VARIATIONAL_VERTEX_BLOCK_DESCENT_H

#include "btVariationalIntegrator.h"

///"Vertex Block Descent", Chen, Liu, Yang, Yuksel, ACM Transactions on Graphics 43(4) (SIGGRAPH 2024),
///with the hard constraints of "Augmented Vertex Block Descent", Giles, Diaz, Yuksel, ACM Transactions
///on Graphics 44(4) (SIGGRAPH 2025).
///
///Block coordinate descent on the incremental potential with one block per vertex: all other vertices
///held fixed, the energy of vertex i is minimized by a single Newton step
///
///     x_i  <-  x_i + H_i^-1 f_i ,
///
///where f_i is the force on the vertex (inertia, elasticity, attachments, contact, friction) and H_i
///the 3x3 block of the Hessian. No global linear system is ever assembled or solved. Vertices are
///visited in Gauss-Seidel order; a graph coloring groups vertices that share no element, and the
///vertices of one color can be updated in parallel (btParallelFor). Each local step reduces the
///energy, which makes the method stable for any iteration budget.
///
///Implemented from the paper and its reference implementation (Gaia):
/// - adaptive initialization: x = x_n + h v_n + h^2 a~, a~ being the previous acceleration along
///   gravity clamped to [0, |g|] (resting bodies start from rest, falling ones from the inertial guess),
/// - optional Chebyshev acceleration of the iteration sequence,
/// - stiffness proportional (Rayleigh) damping, f -= gamma K v, H += gamma/h K with K the elastic block,
/// - smoothed friction evaluated per vertex.
///
///AVBD adds an augmented Lagrangian for hard constraints, here the hard attachments of the system:
///every constraint keeps a dual variable lambda and a penalty stiffness k. A vertex sees the force
///-(k C + lambda) and the Hessian k I; after every sweep lambda <- lambda + k C and k grows by
///beta |C|. Both are warm started with decay across time steps (lambda <- alpha gamma lambda,
///k <- gamma k) and the constraint error present at the start of a step is only removed at the rate
///(1 - alpha) per step, which avoids injecting energy. The constraints converge to zero error even
///when the attached masses and stiffnesses differ by orders of magnitude, which a plain penalty
///spring cannot do.
class btVariationalVertexBlockDescent : public btVariationalIntegrator
{
public:
	enum Initialization
	{
		INIT_PREVIOUS = 0,  // x = x_n
		INIT_INERTIA,       // x = x_n + h v_n + h^2 g
		INIT_ADAPTIVE       // the adaptive blend of the paper
	};

	int m_iterations;          // Gauss-Seidel sweeps per time step
	int m_initialization;
	double m_tolerance;        // [m/s] stop early when a sweep moves no vertex faster than this (0: never)
	bool m_useChebyshev;
	double m_chebyshevRho;     // estimate of the spectral radius of one sweep (0.93 in the reference)
	double m_rayleighDamping;  // gamma [s]
	bool m_parallel;           // update the vertices of a color with btParallelFor

	// augmented Lagrangian of AVBD
	bool m_augmentedLagrangian;
	double m_alpha;       // fraction of the start-of-step constraint error that is kept (stabilization)
	double m_beta;        // growth of the penalty stiffness per unit of constraint error and sweep
	double m_gamma;       // warm start decay of lambda and of the penalty stiffness
	double m_penaltyMin;
	double m_penaltyMax;

	btVariationalVertexBlockDescent();

	virtual const char* getName() const { return m_augmentedLagrangian ? "AugmentedVertexBlockDescent" : "VertexBlockDescent"; }
	virtual void solve(btVariationalSystem& system);

	int numColors() const { return m_colorStart.size() > 0 ? m_colorStart.size() - 1 : 0; }
	///largest constraint error |C| of the hard attachments after the last solve
	double maxConstraintError() const { return m_maxConstraintError; }

	///one Newton step on a single vertex; returns the length of the displacement
	double updateVertex(const btVariationalSystem& system, btVec3d* x, int node) const;
	///one Gauss-Seidel sweep over all colors; returns the largest vertex displacement
	double sweep(const btVariationalSystem& system, btVec3d* x);
	void buildColoring(const btVariationalSystem& system);

private:
	btAlignedObjectArray<int> m_colorStart, m_colorNodes;
	double m_coloringSignature;
	int m_lastTopologyVersion;
	btVec3dArray m_prevIterate, m_prevPrevIterate;
	double m_maxConstraintError;

	void warmStartConstraints(btVariationalSystem& system) const;
	void updateDuals(btVariationalSystem& system, const btVec3d* x);
};

#endif  // BT_VARIATIONAL_VERTEX_BLOCK_DESCENT_H
