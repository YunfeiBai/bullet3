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

#include "VertexBlockDescent.h"
#include "VariationalSceneHelpers.h"
#include "BulletSoftBody/BulletVariationalDeformable/btVariationalVertexBlockDescent.h"
#include "BulletSoftBody/BulletVariationalDeformable/btVariationalIPC.h"

///"Vertex Block Descent" (Chen, Liu, Yang, Yuksel 2024): every vertex is moved by one Newton step on
///its own 3x3 block while the others are held fixed, color by color in Gauss-Seidel order. There is no
///global linear system at all, each sweep is cheap and reduces the energy, and the simulation stays
///stable for any number of sweeps.
///
///A batch of soft cubes rains onto the ground and onto each other. Their contact uses the barrier of
///the IPC demo with per-vertex feasible steps, so the pile is free of interpenetration as well. Use
///the slider to trade sweeps per step against accuracy: with few sweeps the bodies look softer and
///more damped, but nothing blows up.
static btScalar gVbdSweeps = 12;
static btScalar gVbdDamping = 0.0;

class VariationalVertexBlockDescent : public VariationalDemoBase
{
	btVariationalIPC m_ipc;
	btVariationalVertexBlockDescent* m_vbd;

public:
	VariationalVertexBlockDescent(struct GUIHelperInterface* helper) : VariationalDemoBase(helper), m_vbd(0) {}

	void resetCamera()
	{
		float dist = 4.8;
		float pitch = -28;
		float yaw = 30;
		float targetPos[3] = {0, 0.5, 0};
		m_guiHelper->resetCamera(dist, yaw, pitch, targetPos[0], targetPos[1], targetPos[2]);
	}

	virtual const char* methodName() const { return "Vertex Block Descent (Gauss-Seidel over colored vertices)"; }

	virtual void updateStatus()
	{
		VariationalDemoBase::updateStatus();
		double gap = 0;
		const int active = m_ipc.countActivePairs(&system().m_x[0], &gap);
		const size_t len = strlen(m_statusText);
		snprintf(m_statusText + len, sizeof(m_statusText) - len, " | %d colors, %d contact pairs", m_vbd->numColors(), active);
		if (m_logToConsole && (m_stepCount % 30) == 0)
			printf("      contact: %d active pairs of %d candidates, smallest gap %.3e, intersections: %s\n", active, m_ipc.numCandidates(), gap,
				   m_ipc.hasIntersections(&system().m_x[0]) ? "YES" : "none");
	}

	void initPhysics()
	{
		createWorld(btVector3(0, -9.8, 0));
		addGround();
		m_vbd = new btVariationalVertexBlockDescent();
		m_integrators.push_back(m_vbd);
		m_variationalSolver->setIntegrator(m_vbd);
		system().m_damping = 0.5;
		system().m_colliderModel = btVariationalSystem::COLLIDER_BARRIER;
		system().m_colliderMargin = 0.006;
		system().m_colliderStiffness = 100.0;
		m_ipc.m_dhat = 0.006;
		system().m_contact = &m_ipc;

		btVariationalStableNeoHookean* soft = new btVariationalStableNeoHookean(1.5e5, 0.4);
		btVariationalStableNeoHookean* firm = new btVariationalStableNeoHookean(8e5, 0.4);
		m_materials.push_back(soft);
		m_materials.push_back(firm);
		VariationalDemoRandom rnd(3);
		int count = 0;
		for (int layer = 0; layer < 3; layer++)
			for (int ix = 0; ix < 3; ix++)
				for (int iz = 0; iz < 2; iz++, count++)
				{
					// spaced by more than a cube diagonal, so the randomly rotated cubes start apart
					const btVector3 centre(btScalar(0.8) * (ix - 1) + btScalar(0.04) * rnd.next(), btScalar(0.45) + btScalar(0.85) * layer,
										   btScalar(0.8) * (iz - btScalar(0.5)) + btScalar(0.04) * rnd.next());
					btSoftBody* cube = addTetGrid(centre, btVector3(0.2, 0.2, 0.2), 2, 2, 2, 1000, (count % 2) ? firm : soft, 0.5);
					variationalRotateAboutCentre(cube, btQuaternion(btVector3(rnd.next(), rnd.next(), rnd.next()).normalized(), btScalar(0.6) * rnd.next()));
				}
		m_guiHelper->autogenerateGraphicsObjects(m_dynamicsWorld);

		if (m_guiHelper->getParameterInterface())
		{
			SliderParams sweeps("sweeps per step", &gVbdSweeps);
			sweeps.m_minVal = 1;
			sweeps.m_maxVal = 60;
			sweeps.m_clampToIntegers = true;
			m_guiHelper->getParameterInterface()->registerSliderFloatParameter(sweeps);
			SliderParams damping("stiffness damping (ms)", &gVbdDamping);
			damping.m_minVal = 0;
			damping.m_maxVal = 20;
			m_guiHelper->getParameterInterface()->registerSliderFloatParameter(damping);
		}
	}

	virtual void stepSimulation(float deltaTime)
	{
		m_vbd->m_iterations = int(gVbdSweeps + btScalar(0.5));
		m_vbd->m_rayleighDamping = gVbdDamping * 1e-3;
		VariationalDemoBase::stepSimulation(deltaTime);
	}

	virtual void exitPhysics()
	{
		if (m_variationalSolver) system().m_contact = 0;
		VariationalDemoBase::exitPhysics();
	}
};

class CommonExampleInterface* VariationalVertexBlockDescentCreateFunc(struct CommonExampleOptions& options)
{
	return new VariationalVertexBlockDescent(options.m_guiHelper);
}
