/*
 Written by Xuchen Han <xuchenhan2015@u.northwestern.edu>
 
 Bullet Continuous Collision Detection and Physics Library
 Copyright (c) 2019 Google Inc. http://bulletphysics.org
 This software is provided 'as-is', without any express or implied warranty.
 In no event will the authors be held liable for any damages arising from the use of this software.
 Permission is granted to anyone to use this software for any purpose,
 including commercial applications, and to alter it and redistribute it freely,
 subject to the following restrictions:
 1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.
 2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.
 3. This notice may not be removed or altered from any source distribution.
 */

#ifndef BT_DEFORMABLE_MULTIBODY_CONSTRAINT_SOLVER_H
#define BT_DEFORMABLE_MULTIBODY_CONSTRAINT_SOLVER_H

#include "btDeformableBodySolver.h"
#include "BulletDynamics/Featherstone/btMultiBodyConstraintSolver.h"

class btDeformableBodySolver;

// btDeformableMultiBodyConstraintSolver extendsn btMultiBodyConstraintSolver to solve for the contact among rigid/multibody and deformable bodies. Notice that the following constraints
// 1. rigid/multibody against rigid/multibody
// 2. rigid/multibody against deforamble
// 3. deformable against deformable
// 4. deformable self collision
// 5. joint constraints
// are all coupled in this solve.
ATTRIBUTE_ALIGNED16(class)
btDeformableMultiBodyConstraintSolver : public btMultiBodyConstraintSolver
{
	btDeformableBodySolver* m_deformableSolver;

protected:
	// override the iterations method to include deformable/multibody contact
	//    virtual btScalar solveGroupCacheFriendlyIterations(btCollisionObject** bodies,int numBodies,btPersistentManifold** manifoldPtr, int numManifolds,btTypedConstraint** constraints,int numConstraints,const btContactSolverInfo& infoGlobal,btIDebugDraw* debugDrawer);

	// write the velocity of the the solver body to the underlying rigid body
	void solverBodyWriteBack(const btContactSolverInfo& infoGlobal);

	// write the velocity of the underlying rigid body to the the the solver body
	void writeToSolverBody(btCollisionObject * *bodies, int numBodies, const btContactSolverInfo& infoGlobal);

	// let each deformable body knows which solver body is in constact
	void pairDeformableAndSolverBody(btCollisionObject** bodies, int numBodies, int numDeformableBodies, const btContactSolverInfo& infoGlobal);

	// The multibody constraint rows (joint motors, joint limits, multibody contacts...) compute their
	// impulses from the velocity change accumulated in m_data.m_deltaVelocities, while the deformable
	// contacts read and write btMultiBody::m_deltaV. The multibody rows also apply their impulses to
	// btMultiBody::m_deltaV (DIRECTLY_UPDATE_VELOCITY_DURING_SOLVER_ITERATIONS is defined in
	// btMultiBodyConstraintSolver.h), so the deformable contacts see them, but the multibody rows would
	// not see the impulses of the deformable contacts: a joint motor that pushes a gripper finger
	// against a deformable object would keep applying the impulse it needs without the object, and
	// joint limits would not hold against deformable contacts. So after each iteration, the velocity
	// change caused by the deformable contacts is added to m_data.m_deltaVelocities, analogous to
	// solverBodyWriteBack() and writeToSolverBody() for rigid bodies.
	struct btMultiBodySyncData
	{
		btMultiBody* m_multiBody;
		int m_deltaVelIndex;  // offset of the multibody in m_data.m_deltaVelocities
		int m_offset;         // offset of the multibody in m_syncDeltaV
		int m_numDofs;        // number of velocity dofs, including the 6 base dofs
	};
	btAlignedObjectArray<btMultiBodySyncData> m_syncMultiBodies;
	btAlignedObjectArray<btScalar> m_syncDeltaV;  // btMultiBody::m_deltaV before the deformable contact solve

	// find the multibodies that are handled by the multibody constraint rows
	void setupMultiBodySync();
	// remember btMultiBody::m_deltaV before the deformable contact solve
	void backupMultiBodyDeltaV();
	// add the velocity change caused by the deformable contacts to m_data.m_deltaVelocities
	void writeDeformableImpulsesToMultiBodySolver();

	virtual void solveGroupCacheFriendlySplitImpulseIterations(btCollisionObject * *bodies, int numBodies, btCollisionObject** deformableBodies, int numDeformableBodies, btPersistentManifold** manifoldPtr, int numManifolds, btTypedConstraint** constraints, int numConstraints, const btContactSolverInfo& infoGlobal, btIDebugDraw* debugDrawer);

	virtual btScalar solveDeformableGroupIterations(btCollisionObject * *bodies, int numBodies, btCollisionObject** deformableBodies, int numDeformableBodies, btPersistentManifold** manifoldPtr, int numManifolds, btTypedConstraint** constraints, int numConstraints, const btContactSolverInfo& infoGlobal, btIDebugDraw* debugDrawer);

public:
	BT_DECLARE_ALIGNED_ALLOCATOR();

	void setDeformableSolver(btDeformableBodySolver * deformableSolver)
	{
		m_deformableSolver = deformableSolver;
	}

	virtual void solveDeformableBodyGroup(btCollisionObject * *bodies, int numBodies, btCollisionObject** deformableBodies, int numDeformableBodies, btPersistentManifold** manifold, int numManifolds, btTypedConstraint** constraints, int numConstraints, btMultiBodyConstraint** multiBodyConstraints, int numMultiBodyConstraints, const btContactSolverInfo& info, btIDebugDraw* debugDrawer, btDispatcher* dispatcher);
};

#endif /* BT_DEFORMABLE_MULTIBODY_CONSTRAINT_SOLVER_H */
