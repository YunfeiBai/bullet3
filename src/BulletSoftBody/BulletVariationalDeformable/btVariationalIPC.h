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

#ifndef BT_VARIATIONAL_IPC_H
#define BT_VARIATIONAL_IPC_H

#include "btVariationalSystem.h"
#include "btVariationalDistance.h"

///"Incremental Potential Contact: Intersection- and Inversion-free, Large-Deformation Dynamics",
///Li, Ferguson, Schneider, Langlois, Zorin, Panozzo, Jiang, Kaufman, ACM Transactions on Graphics
///39(4) (SIGGRAPH 2020).
///
///Contact between the surfaces of the deformable bodies (including self contact) is a smooth barrier
///potential on the unsigned distance d of every point-triangle and edge-edge pair,
///
///      kappa * b(d^2) ,    b(s) = -(s - dhat^2)^2 ln(s / dhat^2)   for d < dhat, 0 otherwise.
///
///It vanishes smoothly at the activation distance dhat and grows without bound as the pair
///approaches contact, so a minimizer of the incremental potential that is reached through
///intersection-free steps is itself intersection free. The steps are kept intersection free by
///continuous collision detection: maxStepSize() returns the largest feasible step along a search
///direction and the solver's line search starts there. Pieces implemented here:
///  - closest-feature classification and distances with exact derivatives (btVariationalDistance.h),
///  - the mollifier for nearly parallel edge-edge pairs,
///  - barrier stiffness kappa chosen and adapted as in the paper / IPC Toolkit,
///  - positive semi-definite projection of every pair Hessian,
///  - additive CCD on all candidate pairs, with a spatial hash as broad phase,
///  - friction as a lagged, smoothed dissipative potential (Section 5 of the paper).
///
///Plug an instance into btVariationalSystem::m_contact. Projected Newton then is the solver of the
///paper; vertex block descent uses the same potential with per-vertex feasible steps. Rigid colliders
///are not part of this class: btVariationalSystem treats them through their signed distance, and its
///COLLIDER_BARRIER model applies the same log barrier to them.
class btVariationalIPC : public btVariationalContact
{
public:
	// ---- parameters
	double m_dhat;          // activation distance [m]; <= 0: m_dhatRelative times the scene diagonal
	double m_dhatRelative;  // (the paper's default is 1e-3)
	double m_stiffness;     // barrier stiffness kappa; <= 0: chosen automatically each step
	bool m_adaptiveStiffness;
	bool m_selfCollision;       // include pairs within one body
	bool m_friction;            // friction between deformables, coefficient sqrt(mu_a mu_b) of the two bodies
	double m_frictionEpsilon;   // [m/s] slip velocity below which friction behaves like static friction

	btVariationalIPC();

	// ---- btVariationalContact
	virtual void beginStep(const btVariationalSystem& system);
	virtual double energy(const btVec3d* x) const;
	virtual void addGradient(const btVec3d* x, btVec3d* g) const;
	virtual void addHessian(const btVec3d* x, btVariationalHessian& H, btVariationalHessianFilter filter) const;
	virtual void addVertexForceAndHessian(int node, const btVec3d* x, btVec3d& force, btMat3d& hess) const;
	virtual double maxStepSize(const btVec3d* x, const btVec3d* p) const;
	virtual double vertexMaxStepSize(int node, const btVec3d* x, const btVec3d& dx, bool conservative) const;
	virtual double prepareLocalSteps(const btVec3d* x, double radius, const btVec3d* direction = 0);
	virtual bool updateStiffness(const btVec3d* x);

	// ---- queries
	double getDhat() const { return m_dhatActive; }
	double getBarrierStiffness() const { return m_kappa; }
	int numCandidates() const { return m_pairs.size(); }
	///number of candidate pairs closer than dhat and the smallest distance among all candidates
	int countActivePairs(const btVec3d* x, double* minDistance) const;
	///true if any surface triangle is pierced by an edge of the collision mesh (O(n^2), for tests)
	bool hasIntersections(const btVec3d* x) const;

private:
	struct Pair
	{
		int m_n[4];
		bool m_edgeEdge;
		double m_epsX;         // mollifier threshold of edge-edge pairs
		double m_refDistance;  // distance at the reference positions of a local-step candidate set (< 0: unknown)
	};
	struct FrictionPair
	{
		int m_n[4];
		double m_w[4];
		btVec3d m_normal;
		double m_force;  // mu times the lagged normal force
	};
	struct PairValue
	{
		double m_energy;
		btVec3d m_grad[4];
		btMat3d m_hess[4][4];
	};

	const btVariationalSystem* m_system;
	double m_dhatActive, m_kappa, m_kappaMax, m_sceneDiagonal;
	double m_prevMinDistance2;
	mutable bool m_muted;
	// candidates (rebuilt for every search direction)
	mutable btAlignedObjectArray<Pair> m_pairs;
	mutable btAlignedObjectArray<int> m_nodePairStart, m_nodePairs;
	mutable btVec3dArray m_reference;  // positions the local-step candidates were built at
	mutable double m_localRadius;
	btAlignedObjectArray<FrictionPair> m_frictionPairs;
	btAlignedObjectArray<int> m_nodeFrictionStart, m_nodeFriction;
	// scratch of the broad phase
	mutable btVec3dArray m_boxMin, m_boxMax;

	void buildCandidates(const btVec3d* x, const btVec3d* p, double inflation) const;
	///barrier energy of one pair and, for order >= 1 / 2, its gradient / Hessian. Returns false if inactive.
	bool evaluatePair(const Pair& pair, const btVec3d* x, int order, PairValue& out, int onlySlot = -1) const;
	///lower bound of the current distance of a local-step candidate, from how far its nodes have moved
	double distanceLowerBound(const Pair& pair, const btVec3d* x) const;
	void frictionKernel(const FrictionPair& fp, const btVec3d* x, double& energy, btVec3d* grad, btMat3d* hess) const;
	void buildFrictionPairs(const btVec3d* x);
};

#endif  // BT_VARIATIONAL_IPC_H
