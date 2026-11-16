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

#include "QuasiNewton.h"
#include "VariationalDemoBase.h"
#include "BulletSoftBody/BulletVariationalDeformable/btVariationalQuasiNewton.h"

///"Quasi-Newton Methods for Real-Time Simulation of Hyperelastic Materials" (Liu, Bouaziz, Kavan 2017):
///L-BFGS whose initial Hessian is the prefactorized projective dynamics matrix. Unlike projective
///dynamics this handles a general hyperelastic material; here it is Stable Neo-Hookean with a Poisson
///ratio of 0.47.
///
///A bar is clamped at both ends and the right clamp stretches it to almost twice its length and then
///compresses it, so the nearly incompressible material necks and bulges. A cube and a slab of the same
///material drop to the ground. The history size slider goes from 0 (the plain constant-matrix step,
///i.e. what projective dynamics would do) to 10 curvature pairs.
static btScalar gQnHistory = 5;
static btScalar gQnIterations = 12;

class VariationalQuasiNewton : public VariationalDemoBase
{
	btVariationalQuasiNewton* m_qn;
	btSoftBody* m_bar;
	btAlignedObjectArray<int> m_pullNodes;
	btAlignedObjectArray<btVector3> m_pullRest;
	btScalar m_time;

public:
	VariationalQuasiNewton(struct GUIHelperInterface* helper) : VariationalDemoBase(helper), m_qn(0), m_bar(0), m_time(0) {}

	void resetCamera()
	{
		float dist = 4.6;
		float pitch = -18;
		float yaw = 15;
		float targetPos[3] = {0.3, 1.0, 0};
		m_guiHelper->resetCamera(dist, yaw, pitch, targetPos[0], targetPos[1], targetPos[2]);
	}

	virtual const char* methodName() const { return "Quasi-Newton L-BFGS on the projective dynamics matrix, Stable Neo-Hookean"; }

	void initPhysics()
	{
		createWorld(btVector3(0, -9.8, 0));
		addGround();
		m_qn = new btVariationalQuasiNewton();
		m_qn->m_tolerance = 2e-3;
		m_integrators.push_back(m_qn);
		m_variationalSolver->setIntegrator(m_qn);
		system().m_damping = 1.0;
		// The Laplacian does not know about contact, so keep the contact penalty moderate
		// (the reference implementation of the paper uses soft collision springs as well).
		system().m_colliderStiffness = 30.0;
		system().m_colliderMargin = 0.02;

		btVariationalStableNeoHookean* rubber = new btVariationalStableNeoHookean(4e5, 0.47);
		m_materials.push_back(rubber);
		{
			m_bar = addTetGrid(btVector3(-0.4, 1.9, 0), btVector3(0.6, 0.15, 0.15), 12, 3, 3, 1000, rubber);
			for (int i = 0; i < m_bar->m_nodes.size(); i++)
			{
				const btVector3 x = m_bar->m_nodes[i].m_x;
				if (x.x() < btScalar(-0.999)) m_bar->setMass(i, 0);
				if (x.x() > btScalar(0.199))
				{
					m_bar->setMass(i, 0);
					m_pullNodes.push_back(i);
					m_pullRest.push_back(x);
				}
			}
		}
		addTetGrid(btVector3(-1.2, 0.9, 0.9), btVector3(0.3, 0.3, 0.3), 5, 5, 5, 1000, rubber);
		{
			btSoftBody* slab = addTetGrid(btVector3(1.3, 1.0, 0.9), btVector3(0.5, 0.1, 0.35), 8, 2, 6, 1000, rubber);
			slab->rotate(btQuaternion(btVector3(0, 0, 1), 0.9));
		}
		m_guiHelper->autogenerateGraphicsObjects(m_dynamicsWorld);

		if (m_guiHelper->getParameterInterface())
		{
			SliderParams history("L-BFGS history size", &gQnHistory);
			history.m_minVal = 0;
			history.m_maxVal = 10;
			history.m_clampToIntegers = true;
			m_guiHelper->getParameterInterface()->registerSliderFloatParameter(history);
			SliderParams iterations("max iterations", &gQnIterations);
			iterations.m_minVal = 1;
			iterations.m_maxVal = 60;
			iterations.m_clampToIntegers = true;
			m_guiHelper->getParameterInterface()->registerSliderFloatParameter(iterations);
		}
	}

	virtual void stepSimulation(float deltaTime)
	{
		m_qn->m_historySize = int(gQnHistory + btScalar(0.5));
		m_qn->m_maxIterations = int(gQnIterations + btScalar(0.5));
		// right clamp: x offset from +0.9 (stretched) to -0.3 (compressed) with a period of 8 seconds
		m_time += m_fixedTimeStep;
		const btScalar offset = btScalar(0.3) + btScalar(0.6) * btSin(SIMD_2_PI * m_time / btScalar(8));
		for (int i = 0; i < m_pullNodes.size(); i++)
		{
			btSoftBody::Node& n = m_bar->m_nodes[m_pullNodes[i]];
			n.m_x = m_pullRest[i] + btVector3(offset, 0, 0);
			n.m_q = n.m_x;
		}
		VariationalDemoBase::stepSimulation(deltaTime);
	}
};

class CommonExampleInterface* VariationalQuasiNewtonCreateFunc(struct CommonExampleOptions& options)
{
	return new VariationalQuasiNewton(options.m_guiHelper);
}
