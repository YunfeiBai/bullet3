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

#include "ProjectiveDynamics.h"
#include "VariationalDemoBase.h"
#include "BulletSoftBody/BulletVariationalDeformable/btVariationalProjectiveDynamics.h"

///"Fast Simulation of Mass-Spring Systems" (Liu et al. 2013) and "Projective Dynamics" (Bouaziz et al.
///2014): a mass-spring cloth and two tetrahedral solids are advanced by the same local/global solver.
///The global system matrix is factorized once; every iteration afterwards is a cheap projection per
///element followed by a pair of triangular solves, so the cost per step stays flat no matter how
///stiff the materials are. Lower the iteration count to see the typical trade: the result stays
///stable but looks softer.
static btScalar gPdIterations = 10;

class VariationalProjectiveDynamics : public VariationalDemoBase
{
	btVariationalProjectiveDynamics* m_pd;

public:
	VariationalProjectiveDynamics(struct GUIHelperInterface* helper) : VariationalDemoBase(helper), m_pd(0) {}

	void resetCamera()
	{
		float dist = 5.0;
		float pitch = -22;
		float yaw = 25;
		float targetPos[3] = {0.2, 1.0, 0};
		m_guiHelper->resetCamera(dist, yaw, pitch, targetPos[0], targetPos[1], targetPos[2]);
	}

	virtual const char* methodName() const { return "Projective Dynamics (local/global, prefactorized)"; }

	void initPhysics()
	{
		createWorld(btVector3(0, -9.8, 0));
		addGround();
		m_pd = new btVariationalProjectiveDynamics();
		m_integrators.push_back(m_pd);
		m_variationalSolver->setIntegrator(m_pd);
		system().m_damping = 0.5;

		// a sphere for the cloth to drape over
		{
			btCollisionShape* sphere = new btSphereShape(0.45);
			m_collisionShapes.push_back(sphere);
			btTransform t;
			t.setIdentity();
			t.setOrigin(btVector3(0, 0.9, 0.1));
			createRigidBody(0, t, sphere)->setFriction(0.6);
		}
		// mass-spring cloth pinned at two corners (Liu et al. 2013)
		{
			const btScalar s = 0.9, h = 1.9;
			btSoftBody* cloth = btSoftBodyHelpers::CreatePatch(getDeformableDynamicsWorld()->getWorldInfo(),
															   btVector3(-s, h, -s), btVector3(s, h, -s), btVector3(-s, h, s), btVector3(s, h, s),
															   32, 32, 1 + 2, true);
			cloth->generateBendingConstraints(2);
			cloth->setTotalMass(0.4);
			addSoftBody(cloth, 0, 0.4);
			m_variationalSolver->setSpringStiffness(cloth, 800.0, 4.0);
		}
		// tetrahedral solids with the strain constraint of projective dynamics (Bouaziz et al. 2014)
		btVariationalARAP* soft = new btVariationalARAP(2e5, 0.3);
		btVariationalARAP* stiff = new btVariationalARAP(2e6, 0.3);
		m_materials.push_back(soft);
		m_materials.push_back(stiff);
		{
			// cantilever
			btSoftBody* bar = addTetGrid(btVector3(-1.9, 1.3, 0.3), btVector3(0.6, 0.1, 0.15), 12, 2, 3, 1000, stiff);
			for (int i = 0; i < bar->m_nodes.size(); i++)
				if (bar->m_nodes[i].m_x.x() < btScalar(-2.499)) bar->setMass(i, 0);
		}
		{
			// block dropped on its corner
			btSoftBody* block = addTetGrid(btVector3(1.9, 1.4, 0.3), btVector3(0.3, 0.3, 0.3), 5, 5, 5, 1000, soft);
			block->rotate(btQuaternion(btVector3(1, 0, 1).normalized(), 0.7));
		}
		m_guiHelper->autogenerateGraphicsObjects(m_dynamicsWorld);

		if (m_guiHelper->getParameterInterface())
		{
			SliderParams iterations("local/global iterations", &gPdIterations);
			iterations.m_minVal = 1;
			iterations.m_maxVal = 50;
			iterations.m_clampToIntegers = true;
			m_guiHelper->getParameterInterface()->registerSliderFloatParameter(iterations);
		}
	}

	virtual void stepSimulation(float deltaTime)
	{
		m_pd->m_iterations = int(gPdIterations + btScalar(0.5));
		VariationalDemoBase::stepSimulation(deltaTime);
	}
};

class CommonExampleInterface* VariationalProjectiveDynamicsCreateFunc(struct CommonExampleOptions& options)
{
	return new VariationalProjectiveDynamics(options.m_guiHelper);
}
