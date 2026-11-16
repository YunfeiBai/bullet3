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

#include "IncrementalPotentialContact.h"
#include "VariationalSceneHelpers.h"
#include "BulletSoftBody/BulletVariationalDeformable/btVariationalIPC.h"

///"Incremental Potential Contact" (Li et al. 2020): contact between the deformable bodies is a log
///barrier on the distance of every point-triangle and edge-edge pair, the step length of the Newton
///solver is limited by continuous collision detection, and friction is a lagged smooth potential.
///Bodies therefore never interpenetrate, at any time step size and for any stiffness.
///
///Soft cubes and slabs are dropped on top of each other; the ground is a barrier as well. The status
///line shows the number of active contact pairs and the smallest gap between any two surfaces, which
///stays strictly positive.
static btScalar gIpcFriction = 0.5;

class VariationalIncrementalPotentialContact : public VariationalDemoBase
{
	btVariationalIPC m_ipc;
	btVariationalProjectedNewton* m_newton;

public:
	VariationalIncrementalPotentialContact(struct GUIHelperInterface* helper) : VariationalDemoBase(helper), m_newton(0) {}

	void resetCamera()
	{
		float dist = 4.2;
		float pitch = -22;
		float yaw = 35;
		float targetPos[3] = {0, 0.7, 0};
		m_guiHelper->resetCamera(dist, yaw, pitch, targetPos[0], targetPos[1], targetPos[2]);
	}

	virtual const char* methodName() const { return "Incremental Potential Contact, projected Newton"; }

	virtual void updateStatus()
	{
		VariationalDemoBase::updateStatus();
		double gap = 0;
		const int active = m_ipc.countActivePairs(&system().m_x[0], &gap);
		const size_t len = strlen(m_statusText);
		snprintf(m_statusText + len, sizeof(m_statusText) - len, " | %d contact pairs, smallest gap %.2f mm", active, gap < 1e100 ? gap * 1000.0 : 0.0);
		if (m_logToConsole && (m_stepCount % 30) == 0)
			printf("      contact: %d active pairs of %d candidates, smallest gap %.3e, intersections: %s\n", active, m_ipc.numCandidates(), gap,
				   m_ipc.hasIntersections(&system().m_x[0]) ? "YES" : "none");
	}

	void initPhysics()
	{
		createWorld(btVector3(0, -9.8, 0));
		addGround();
		m_newton = new btVariationalProjectedNewton();
		m_newton->m_tolerance = 1e-2;
		m_newton->m_linearTolerance = 0.05;
		m_newton->m_maxIterations = 30;
		m_integrators.push_back(m_newton);
		m_variationalSolver->setIntegrator(m_newton);
		system().m_damping = 1.0;
		// the ground is a barrier too, so nothing ever sinks into it
		system().m_colliderModel = btVariationalSystem::COLLIDER_BARRIER;
		system().m_colliderMargin = 0.004;
		system().m_colliderStiffness = 100.0;
		m_ipc.m_dhat = 0.004;
		system().m_contact = &m_ipc;

		btVariationalStableNeoHookean* rubber = new btVariationalStableNeoHookean(4e5, 0.4);
		m_materials.push_back(rubber);
		VariationalDemoRandom rnd(7);
		// a leaning tower of cubes ...
		for (int i = 0; i < 5; i++)
		{
			const btScalar half = btScalar(0.22) - btScalar(0.02) * i;
			btSoftBody* cube = addTetGrid(btVector3(btScalar(0.06) * rnd.next(), btScalar(0.3) + btScalar(0.62) * i, btScalar(0.06) * rnd.next()),
										  btVector3(half, half, half), 2, 2, 2, 1000, rubber, gIpcFriction);
			variationalRotateAboutCentre(cube, btQuaternion(btVector3(0, 1, 0), btScalar(0.5) * rnd.next()));
		}
		// ... and two slabs that land across it
		for (int i = 0; i < 2; i++)
		{
			btSoftBody* slab = addTetGrid(btVector3(btScalar(0.5) - btScalar(1.0) * i, btScalar(3.6) + btScalar(0.5) * i, btScalar(0.1) * i),
										  btVector3(0.55, 0.06, 0.25), 6, 1, 3, 1000, rubber, gIpcFriction);
			variationalRotateAboutCentre(slab, btQuaternion(btVector3(0, 0, 1), btScalar(0.4) - btScalar(0.9) * i));
		}
		m_guiHelper->autogenerateGraphicsObjects(m_dynamicsWorld);

		if (m_guiHelper->getParameterInterface())
		{
			SliderParams friction("friction coefficient", &gIpcFriction);
			friction.m_minVal = 0;
			friction.m_maxVal = 1;
			m_guiHelper->getParameterInterface()->registerSliderFloatParameter(friction);
		}
	}

	virtual void stepSimulation(float deltaTime)
	{
		btSoftBodyArray& bodies = getDeformableDynamicsWorld()->getSoftBodyArray();
		for (int i = 0; i < bodies.size(); i++) bodies[i]->m_cfg.kDF = gIpcFriction;
		VariationalDemoBase::stepSimulation(deltaTime);
	}

	virtual void exitPhysics()
	{
		if (m_variationalSolver) system().m_contact = 0;
		VariationalDemoBase::exitPhysics();
	}
};

class CommonExampleInterface* VariationalIncrementalPotentialContactCreateFunc(struct CommonExampleOptions& options)
{
	return new VariationalIncrementalPotentialContact(options.m_guiHelper);
}
