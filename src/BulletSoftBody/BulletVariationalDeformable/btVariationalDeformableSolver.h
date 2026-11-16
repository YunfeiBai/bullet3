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

#ifndef BT_VARIATIONAL_DEFORMABLE_SOLVER_H
#define BT_VARIATIONAL_DEFORMABLE_SOLVER_H

#include "BulletSoftBody/btDeformableBodySolver.h"
#include "btVariationalSystem.h"
#include "btVariationalIntegrator.h"
#include "btVariationalProjectedNewton.h"

///Deformable body solver for btDeformableMultiBodyDynamicsWorld that advances all soft bodies with an
///optimization time integrator: each step minimizes the incremental potential of implicit Euler
///(see btVariationalSystem) with the selected btVariationalIntegrator.
///
///Usage mirrors the stock deformable world, only the solver object differs:
///
///   btVariationalDeformableSolver* solver = new btVariationalDeformableSolver();
///   btDeformableMultiBodyConstraintSolver* sol = new btDeformableMultiBodyConstraintSolver();
///   sol->setDeformableSolver(solver);
///   world = new btDeformableMultiBodyDynamicsWorld(dispatcher, broadphase, sol, collisionConfig, solver);
///   solver->setIntegrator(&myIntegrator);          // projected Newton by default
///   solver->setMaterial(psb, &myMaterial);         // per soft body
///
///Soft bodies with tetrahedra are simulated as finite element solids, bodies with links only as
///mass-spring systems. Nodes with zero mass are kinematic: they follow btSoftBody::Node::m_x.
///Rigid bodies and other collision objects in the world act as colliders through their signed
///distance (planes, boxes, spheres and capsules exactly, other convex shapes through Bullet's sparse
///SDF; triangle meshes are not supported). With setTwoWayCoupling(true) dynamic rigid bodies
///receive the contact reaction as an impulse at the end of the step.
///
///The elastic forces (btDeformableLagrangianForce) and the contact pipeline of the stock solver are
///not used by this solver.
class btVariationalDeformableSolver : public btDeformableBodySolver
{
protected:
	struct BodySettings
	{
		btSoftBody* m_psb;
		const btVariationalMaterial* m_material;
		double m_springStiffness, m_bendingStiffness;
		int m_numNodes, m_numTets, m_numLinks, m_numFaces;
	};
	struct AttachmentSettings
	{
		btSoftBody* m_psbA;
		int m_nodeA;
		btSoftBody* m_psbB;
		int m_nodeB;
		btVector3 m_target;
		btScalar m_stiffness;
		bool m_hard;
	};
	struct PendingCollider
	{
		btSoftBody* m_psb;
		const btCollisionObject* m_object;
		const btCollisionShape* m_shape;
		btTransform m_localTransform;  // of the shape relative to the collision object
	};

	btVariationalSystem m_system;
	btVariationalIntegrator* m_integrator;
	btVariationalProjectedNewton m_defaultIntegrator;
	btVariationalStableNeoHookean m_defaultMaterial;
	btAlignedObjectArray<BodySettings> m_settings;
	btAlignedObjectArray<AttachmentSettings> m_attachmentSettings;
	btAlignedObjectArray<PendingCollider> m_pendingColliders;
	bool m_dirty;
	int m_firstPickingAttachment;
	double m_lastSolveSeconds;

	BodySettings& settingsOf(btSoftBody* psb);
	void rebuildIfNeeded();
	void syncPickingForces();

public:
	btVariationalDeformableSolver();
	virtual ~btVariationalDeformableSolver() {}

	// ---- configuration
	void setIntegrator(btVariationalIntegrator* integrator) { m_integrator = integrator ? integrator : &m_defaultIntegrator; }
	btVariationalIntegrator* getIntegrator() { return m_integrator; }
	///material of the tetrahedra of a soft body (not owned; must outlive the solver)
	void setMaterial(btSoftBody* psb, const btVariationalMaterial* material);
	///stiffness [N/m] of structural and bending links of a soft body without tetrahedra
	void setSpringStiffness(btSoftBody* psb, btScalar structural, btScalar bending);
	void setTwoWayCoupling(bool enable) { m_system.m_twoWayCoupling = enable; }
	btVariationalSystem& getSystem() { return m_system; }
	const btVariationalSystem& getSystem() const { return m_system; }
	///wall clock time of the last deformable solve in seconds
	double getLastSolveSeconds() const { return m_lastSolveSeconds; }

	///Ties a node to a world-space point. Hard attachments are enforced by the augmented Lagrangian of
	///AVBD when the vertex block descent or the hybrid integrator is used; the other integrators
	///treat them as springs of the given stiffness.
	int addAttachment(btSoftBody* psb, int node, const btVector3& target, btScalar stiffness, bool hard = false);
	///Ties two nodes (of the same or of different soft bodies) together.
	int addAttachment(btSoftBody* psbA, int nodeA, btSoftBody* psbB, int nodeB, btScalar stiffness, bool hard = false);
	void setAttachmentTarget(int attachment, const btVector3& target);
	void removeAllAttachments();

	// ---- btDeformableBodySolver interface
	virtual void reinitialize(const btAlignedObjectArray<btSoftBody*>& softBodies, btScalar dt);
	virtual void predictMotion(btScalar solverdt);
	virtual void processCollision(btSoftBody* softBody, const btCollisionObjectWrapper* collisionObjectWrap);
	virtual void processCollision(btSoftBody* softBody, btSoftBody* otherSoftBody)
	{
		// deformable-deformable contact is part of the variational system (btVariationalContact)
		(void)softBody;
		(void)otherSoftBody;
	}
	virtual void setConstraints(const btContactSolverInfo& infoGlobal) { (void)infoGlobal; }
	virtual btScalar solveContactConstraints(btCollisionObject** deformableBodies, int numDeformableBodies, const btContactSolverInfo& infoGlobal)
	{
		(void)deformableBodies;
		(void)numDeformableBodies;
		(void)infoGlobal;
		return 0;
	}
	virtual void setupDeformableSolve(bool implicit) { (void)implicit; }
	virtual void solveDeformableConstraints(btScalar solverdt);
	virtual void applyTransforms(btScalar timeStep);
	virtual void applyExplicitForce() {}
	virtual void setProjection() {}
	virtual void setLagrangeMultiplier() {}
	virtual void setGravity(const btVector3& gravity) { m_system.m_gravity = btVec3d(gravity); }
};

#endif  // BT_VARIATIONAL_DEFORMABLE_SOLVER_H
