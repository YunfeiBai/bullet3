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

#include "StableNeoHookean.h"
#include "VariationalDemoBase.h"

///"Stable Neo-Hookean Flesh Simulation" (Smith, de Goes, Kim, 2018) solved with projected Newton,
///optionally with the absolute-value eigenvalue filter of "Stabler Neo-Hookean Simulation"
///(Chen et al., 2024).
///
///The scene exercises what the material is designed for:
/// - a bar clamped at both ends whose right clamp is twisted back and forth by one and a half turns,
/// - a cube that starts squashed to zero thickness (every element degenerate),
/// - a cube that starts with all of its nodes at random positions (about half of the elements inverted).
///All three find their way back to a valid configuration and the last two drop to the ground.
static btScalar gSnhYoungsModulus = 5e5;
static btScalar gSnhPoissonRatio = 0.45;
static btScalar gSnhAbsFilter = 0;

class VariationalStableNeoHookean : public VariationalDemoBase
{
	btVariationalStableNeoHookean* m_material;
	btVariationalProjectedNewton* m_newton;
	btSoftBody* m_bar;
	btAlignedObjectArray<int> m_twistNodes;
	btAlignedObjectArray<btVector3> m_twistRest;
	btScalar m_time;

public:
	VariationalStableNeoHookean(struct GUIHelperInterface* helper)
		: VariationalDemoBase(helper), m_material(0), m_newton(0), m_bar(0), m_time(0)
	{
	}

	void resetCamera()
	{
		float dist = 4.2;
		float pitch = -18;
		float yaw = 20;
		float targetPos[3] = {0, 1.1, 0.3};
		m_guiHelper->resetCamera(dist, yaw, pitch, targetPos[0], targetPos[1], targetPos[2]);
	}

	virtual const char* methodName() const
	{
		return gSnhAbsFilter > btScalar(0.5) ? "Stable Neo-Hookean, projected Newton (abs filter)" : "Stable Neo-Hookean, projected Newton (clamp filter)";
	}

	void initPhysics()
	{
		createWorld(btVector3(0, -9.8, 0));
		addGround();

		m_material = new btVariationalStableNeoHookean(gSnhYoungsModulus, gSnhPoissonRatio);
		m_materials.push_back(m_material);
		m_newton = new btVariationalProjectedNewton();
		m_newton->m_tolerance = 5e-3;
		m_newton->m_linearTolerance = 0.05;
		m_newton->m_maxIterations = 30;
		m_integrators.push_back(m_newton);
		m_variationalSolver->setIntegrator(m_newton);
		system().m_damping = 1.0;

		// twisted bar, clamped at both ends
		{
			m_bar = addTetGrid(btVector3(0, 2.0, 0), btVector3(1.0, 0.15, 0.15), 16, 3, 3, 1000, m_material);
			for (int i = 0; i < m_bar->m_nodes.size(); i++)
			{
				const btVector3 x = m_bar->m_nodes[i].m_x;
				if (x.x() < btScalar(-0.999)) m_bar->setMass(i, 0);
				if (x.x() > btScalar(0.999))
				{
					m_bar->setMass(i, 0);
					m_twistNodes.push_back(i);
					m_twistRest.push_back(x);
				}
			}
		}
		// cube squashed flat: rest shape first, then collapse the positions
		{
			btSoftBody* psb = addTetGrid(btVector3(-1.6, 0.8, 1.2), btVector3(0.3, 0.3, 0.3), 4, 4, 4, 1000, m_material);
			for (int i = 0; i < psb->m_nodes.size(); i++)
			{
				psb->m_nodes[i].m_x.setY(0.8);
				psb->m_nodes[i].m_q = psb->m_nodes[i].m_x;
			}
		}
		// cube with scrambled nodes
		{
			btSoftBody* psb = addTetGrid(btVector3(1.6, 0.8, 1.2), btVector3(0.3, 0.3, 0.3), 4, 4, 4, 1000, m_material);
			unsigned int seed = 12345;
			for (int i = 0; i < psb->m_nodes.size(); i++)
			{
				btVector3 r;
				for (int c = 0; c < 3; c++)
				{
					seed = 1664525u * seed + 1013904223u;
					r[c] = btScalar((seed >> 8) & 0xffff) / btScalar(65535) - btScalar(0.5);
				}
				psb->m_nodes[i].m_x = btVector3(1.6, 0.8, 1.2) + r * btScalar(0.6);
				psb->m_nodes[i].m_q = psb->m_nodes[i].m_x;
			}
		}
		m_guiHelper->autogenerateGraphicsObjects(m_dynamicsWorld);

		if (m_guiHelper->getParameterInterface())
		{
			SliderParams youngs("Young's modulus", &gSnhYoungsModulus);
			youngs.m_minVal = 5e4;
			youngs.m_maxVal = 5e6;
			m_guiHelper->getParameterInterface()->registerSliderFloatParameter(youngs);
			SliderParams poisson("Poisson ratio", &gSnhPoissonRatio);
			poisson.m_minVal = 0.05;
			poisson.m_maxVal = 0.495;
			m_guiHelper->getParameterInterface()->registerSliderFloatParameter(poisson);
			SliderParams filter("abs eigenvalue filter (0/1)", &gSnhAbsFilter);
			filter.m_minVal = 0;
			filter.m_maxVal = 1;
			m_guiHelper->getParameterInterface()->registerSliderFloatParameter(filter);
		}
	}

	virtual void stepSimulation(float deltaTime)
	{
		m_material->setParameters(gSnhYoungsModulus, gSnhPoissonRatio, true);
		system().m_hessianFilter = gSnhAbsFilter > btScalar(0.5) ? BT_VARIATIONAL_FILTER_ABS : BT_VARIATIONAL_FILTER_CLAMP;
		// twist the right clamp about the bar axis: +-1.5 turns with a period of 12 seconds
		m_time += m_fixedTimeStep;
		const btScalar angle = btScalar(1.5) * SIMD_2_PI * btSin(SIMD_2_PI * m_time / btScalar(12));
		const btQuaternion q(btVector3(1, 0, 0), angle);
		const btVector3 axisPoint(0, 2.0, 0);
		for (int i = 0; i < m_twistNodes.size(); i++)
		{
			btSoftBody::Node& n = m_bar->m_nodes[m_twistNodes[i]];
			n.m_x = axisPoint + quatRotate(q, m_twistRest[i] - axisPoint);
			n.m_q = n.m_x;
		}
		VariationalDemoBase::stepSimulation(deltaTime);
	}
};

class CommonExampleInterface* VariationalStableNeoHookeanCreateFunc(struct CommonExampleOptions& options)
{
	return new VariationalStableNeoHookean(options.m_guiHelper);
}
