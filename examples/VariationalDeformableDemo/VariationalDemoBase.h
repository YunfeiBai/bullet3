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

#ifndef VARIATIONAL_DEMO_BASE_H
#define VARIATIONAL_DEMO_BASE_H

#include "btBulletDynamicsCommon.h"
#include "BulletSoftBody/btDeformableMultiBodyDynamicsWorld.h"
#include "BulletSoftBody/btSoftBody.h"
#include "BulletSoftBody/btSoftBodyHelpers.h"
#include "BulletSoftBody/btSoftBodyRigidBodyCollisionConfiguration.h"
#include "BulletSoftBody/BulletVariationalDeformable/btVariationalDeformableSolver.h"
#include "BulletSoftBody/BulletVariationalDeformable/btVariationalHelpers.h"
#include "../CommonInterfaces/CommonDeformableBodyBase.h"
#include "../CommonInterfaces/CommonParameterInterface.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

///Shared scaffolding of the variational deformable demos: a btDeformableMultiBodyDynamicsWorld driven by
///btVariationalDeformableSolver, helpers to add static geometry and soft bodies, and an on-screen line
///with the solver statistics.
class VariationalDemoBase : public CommonDeformableBodyBase
{
protected:
	btVariationalDeformableSolver* m_variationalSolver;
	btAlignedObjectArray<btVariationalMaterial*> m_materials;      // owned by the demo
	btAlignedObjectArray<btVariationalIntegrator*> m_integrators;  // owned by the demo
	btScalar m_fixedTimeStep;
	int m_stepCount;
	double m_averageSolveMs, m_averageIterations;
	char m_statusText[256];
	bool m_logToConsole;

	void createWorld(const btVector3& gravity)
	{
		m_guiHelper->setUpAxis(1);
		m_collisionConfiguration = new btSoftBodyRigidBodyCollisionConfiguration();
		m_dispatcher = new btCollisionDispatcher(m_collisionConfiguration);
		m_broadphase = new btDbvtBroadphase();
		m_variationalSolver = new btVariationalDeformableSolver();
		btDeformableMultiBodyConstraintSolver* sol = new btDeformableMultiBodyConstraintSolver();
		sol->setDeformableSolver(m_variationalSolver);
		m_solver = sol;
		m_dynamicsWorld = new btDeformableMultiBodyDynamicsWorld(m_dispatcher, m_broadphase, sol, m_collisionConfiguration, m_variationalSolver);
		m_dynamicsWorld->setGravity(gravity);
		getDeformableDynamicsWorld()->getWorldInfo().m_gravity = gravity;
		getDeformableDynamicsWorld()->getWorldInfo().m_sparsesdf.setDefaultVoxelsz(0.25);
		getDeformableDynamicsWorld()->getWorldInfo().m_sparsesdf.Reset();
		m_guiHelper->createPhysicsDebugDrawer(m_dynamicsWorld);
		// picking: spring stiffness in N/m of the mouse spring (three nodes are attached)
		m_pickingForceElasticStiffness = 2000;
		m_pickingForceDampingStiffness = 0;
		m_maxPickingForce = 1e10;
	}

	btRigidBody* addStaticBox(const btVector3& halfExtents, const btVector3& position, const btQuaternion& rotation = btQuaternion::getIdentity(), btScalar friction = 0.8)
	{
		btCollisionShape* shape = new btBoxShape(halfExtents);
		m_collisionShapes.push_back(shape);
		btTransform t;
		t.setIdentity();
		t.setOrigin(position);
		t.setRotation(rotation);
		btRigidBody* body = createRigidBody(0, t, shape);
		body->setFriction(friction);
		return body;
	}

	btRigidBody* addGround(btScalar friction = 0.8)
	{
		return addStaticBox(btVector3(50, 1, 50), btVector3(0, -1, 0), btQuaternion::getIdentity(), friction);
	}

	///registers a soft body with the world and assigns its material
	void addSoftBody(btSoftBody* psb, const btVariationalMaterial* material, btScalar friction = 0.5)
	{
		psb->getCollisionShape()->setMargin(0.01);
		psb->m_cfg.kDF = friction;
		getDeformableDynamicsWorld()->addSoftBody(psb);
		m_variationalSolver->setMaterial(psb, material);
	}

	btSoftBody* addTetGrid(const btVector3& center, const btVector3& halfExtents, int nx, int ny, int nz, btScalar density,
						   const btVariationalMaterial* material, btScalar friction = 0.5)
	{
		btSoftBody* psb = btVariationalHelpers::createTetGrid(getDeformableDynamicsWorld()->getWorldInfo(), center, halfExtents, nx, ny, nz, density);
		addSoftBody(psb, material, friction);
		return psb;
	}

	btVariationalSystem& system() { return m_variationalSolver->getSystem(); }

	///one-line description of the configured method, shown in front of the statistics
	virtual const char* methodName() const { return m_variationalSolver->getIntegrator()->getName(); }

	virtual void updateStatus()
	{
		const btVariationalSolveStats& stats = m_variationalSolver->getIntegrator()->m_stats;
		const double ms = m_variationalSolver->getLastSolveSeconds() * 1000.0;
		const double blend = m_stepCount < 10 ? 1.0 / (m_stepCount + 1) : 0.1;
		m_averageSolveMs += blend * (ms - m_averageSolveMs);
		m_averageIterations += blend * (stats.m_iterations - m_averageIterations);
		snprintf(m_statusText, sizeof(m_statusText), "%s | %d nodes, %d tets | %.1f iterations, %.1f ms per step",
				 methodName(), system().numNodes(), system().m_tets.size(), m_averageIterations, m_averageSolveMs);
		if (m_logToConsole && (m_stepCount % 30) == 0)
		{
			printf("[step %4d] %s | |grad| %.2e -> %.2e%s\n", m_stepCount, m_statusText, stats.m_initialGradient, stats.m_finalGradient, stats.m_converged ? "" : " (stopped at the iteration limit)");
			// per body: inverted elements, volume relative to rest, lowest node
			const btVariationalSystem& sys = system();
			for (int b = 0; b < sys.m_bodies.size(); b++)
			{
				const btVariationalBody& body = sys.m_bodies[b];
				int inverted = 0;
				double volume = 0, rest = 0, lowest = 1e30;
				for (int t = body.m_firstTet; t < body.m_firstTet + body.m_numTets; t++)
				{
					const double J = sys.deformationGradient(sys.m_tets[t], &sys.m_x[0]).det();
					if (J <= 0) inverted++;
					volume += J * sys.m_tets[t].m_restVolume;
					rest += sys.m_tets[t].m_restVolume;
				}
				for (int i = body.m_firstNode; i < body.m_firstNode + body.m_numNodes; i++)
					if (sys.m_x[i][1] < lowest) lowest = sys.m_x[i][1];
				printf("      body %d: %d inverted of %d tets, volume ratio %.4f, lowest y %.4f\n", b, inverted, body.m_numTets, rest > 0 ? volume / rest : 0.0, lowest);
			}
		}
	}

public:
	VariationalDemoBase(struct GUIHelperInterface* helper)
		: CommonDeformableBodyBase(helper),
		  m_variationalSolver(0),
		  m_fixedTimeStep(btScalar(1.) / btScalar(60.)),
		  m_stepCount(0),
		  m_averageSolveMs(0),
		  m_averageIterations(0),
		  m_logToConsole(getenv("BT_VARIATIONAL_DEMO_LOG") != 0)
	{
		m_statusText[0] = 0;
	}

	virtual ~VariationalDemoBase() {}

	virtual void stepSimulation(float deltaTime)
	{
		(void)deltaTime;
		// one implicit step per frame
		m_variationalSolver->getIntegrator()->m_verbose = m_logToConsole && getenv("BT_VARIATIONAL_DEMO_VERBOSE") != 0 && (m_stepCount % 20) == 19;
		m_dynamicsWorld->stepSimulation(m_fixedTimeStep, 0, m_fixedTimeStep);
		updateStatus();
		m_stepCount++;
	}

	virtual void renderScene()
	{
		CommonDeformableBodyBase::renderScene();
		if (m_guiHelper->getAppInterface() && m_statusText[0])
		{
			float color[4] = {0.1f, 0.1f, 0.1f, 1.0f};
			m_guiHelper->getAppInterface()->drawText(m_statusText, 260, 60, 0.45f, color);
		}
	}

	virtual void exitPhysics()
	{
		removePickingConstraint();
		for (int i = m_dynamicsWorld->getNumCollisionObjects() - 1; i >= 0; i--)
		{
			btCollisionObject* obj = m_dynamicsWorld->getCollisionObjectArray()[i];
			btRigidBody* body = btRigidBody::upcast(obj);
			if (body && body->getMotionState()) delete body->getMotionState();
			m_dynamicsWorld->removeCollisionObject(obj);
			delete obj;
		}
		for (int j = 0; j < m_forces.size(); j++) delete m_forces[j];
		m_forces.clear();
		for (int j = 0; j < m_collisionShapes.size(); j++) delete m_collisionShapes[j];
		m_collisionShapes.clear();
		delete m_dynamicsWorld;
		m_dynamicsWorld = 0;
		delete m_solver;
		m_solver = 0;
		delete m_variationalSolver;
		m_variationalSolver = 0;
		delete m_broadphase;
		m_broadphase = 0;
		delete m_dispatcher;
		m_dispatcher = 0;
		delete m_collisionConfiguration;
		m_collisionConfiguration = 0;
		for (int j = 0; j < m_integrators.size(); j++) delete m_integrators[j];
		m_integrators.clear();
		for (int j = 0; j < m_materials.size(); j++) delete m_materials[j];
		m_materials.clear();
	}
};

#endif  // VARIATIONAL_DEMO_BASE_H
