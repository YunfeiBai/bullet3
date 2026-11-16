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

#include "VariationalCombined.h"
#include "VariationalSceneHelpers.h"
#include "BulletSoftBody/BulletVariationalDeformable/btVariationalHybrid.h"
#include "BulletSoftBody/BulletVariationalDeformable/btVariationalIPC.h"

///One scene, one potential, the pieces of the variational module working together (btVariationalHybrid):
///
///  - a tower of Stable Neo-Hookean cubes of alternating stiffness (Smith et al. 2018) that rest on each
///    other through the barrier contact of Incremental Potential Contact (Li et al. 2020), with friction,
///  - a wrecking ball: a chain of light blocks carrying a block fifty times heavier, joined by hard
///    attachments that the augmented Lagrangian of Augmented VBD (Giles et al. 2025) keeps closed. Its
///    anchor swings and the ball knocks the tower over,
///  - a mass-spring cloth (the model of Liu et al. 2013) that falls over a soft cube and drapes; cloth,
///    cube and ground all interact through the same contact model,
///  - a clamped bar that is twisted back and forth by a full turn (inversion-safe material, filtered
///    Hessians; Smith et al. 2018, Chen et al. 2024).
///
///The hybrid integrator chooses how to solve each step: projected Newton with a continuous collision
///detection line search whenever contact is possible, and the quasi-Newton method of Liu et al. 2017
///on the prefactorized projective dynamics matrix (Bouaziz et al. 2014) for steps without contact or
///hard joints. In this scene something always touches, so every step takes the Newton path; the
///status line shows the path taken, the smallest gap between surfaces (nothing ever interpenetrates)
///and the joint error.
///
///The solver slider switches to the individual methods on the same scene for comparison:
///0 hybrid, 1 projected Newton alone (the joints are then mere springs and stretch), 2 vertex block
///descent with augmented Lagrangian joints (20 sweeps per step: the cheapest, not converged and
///therefore a little softer), 3 quasi-Newton alone (the barrier is poorly modelled by its constant
///matrix: many iterations).
static btScalar gCombinedSolver = 0;

class VariationalCombined : public VariationalDemoBase
{
	btVariationalIPC m_ipc;
	btVariationalHybrid* m_hybrid;
	btVariationalProjectedNewton* m_newton;
	btVariationalVertexBlockDescent* m_vbd;
	btVariationalQuasiNewton* m_quasiNewton;
	btAlignedObjectArray<int> m_anchorAttachments;
	btAlignedObjectArray<btVector3> m_anchorRest;
	btSoftBody* m_bar;
	btAlignedObjectArray<int> m_twistNodes;
	btAlignedObjectArray<btVector3> m_twistRest;
	btVector3 m_barCentre;
	btScalar m_time;
	int m_pathCount[5];

	static int corner(int n, int ix, int iy, int iz) { return (ix * n) + (n + 1) * ((iy * n) + (n + 1) * (iz * n)); }

public:
	VariationalCombined(struct GUIHelperInterface* helper)
		: VariationalDemoBase(helper), m_hybrid(0), m_newton(0), m_vbd(0), m_quasiNewton(0), m_bar(0), m_time(0)
	{
		for (int i = 0; i < 5; i++) m_pathCount[i] = 0;
	}

	void resetCamera()
	{
		float dist = 4.3;
		float pitch = -18;
		float yaw = 170;
		float targetPos[3] = {0.25, 0.85, 0};
		m_guiHelper->resetCamera(dist, yaw, pitch, targetPos[0], targetPos[1], targetPos[2]);
	}

	virtual const char* methodName() const
	{
		const int s = int(gCombinedSolver + btScalar(0.5));
		return s == 1 ? "Projected Newton alone" : (s == 2 ? "Vertex block descent alone" : (s == 3 ? "Quasi-Newton alone" : "Hybrid"));
	}

	virtual void updateStatus()
	{
		VariationalDemoBase::updateStatus();
		const btVariationalSystem& sys = system();
		double gap = 0, jointError = 0;
		const int active = m_ipc.countActivePairs(&sys.m_x[0], &gap);
		for (int a = 0; a < sys.m_attachments.size(); a++)
		{
			const btVariationalAttachment& at = sys.m_attachments[a];
			if (at.m_hard) jointError = btMax(jointError, (sys.m_x[at.m_nodeA] - (at.m_nodeB >= 0 ? sys.m_x[at.m_nodeB] : at.m_target)).length());
		}
		const bool hybrid = m_variationalSolver->getIntegrator() == m_hybrid;
		if (hybrid) m_pathCount[m_hybrid->lastPath()]++;
		const size_t len = strlen(m_statusText);
		snprintf(m_statusText + len, sizeof(m_statusText) - len, " | %d contacts, gap %.2f mm, joints %.2f mm%s%s", active, gap < 1e100 ? gap * 1000.0 : 0.0,
				 jointError * 1000.0, hybrid ? " | " : "", hybrid ? m_hybrid->lastPathName() : "");
		if (m_logToConsole && (m_stepCount % 30) == 0)
			printf("      %d active pairs, smallest gap %.3e, joint error %.3e, intersections: %s, hybrid steps by path: quasi-Newton %d, Newton %d, Newton+Laplacian %d, QN->Newton %d\n",
				   active, gap, jointError, m_ipc.hasIntersections(&sys.m_x[0]) ? "YES" : "none", m_pathCount[1], m_pathCount[2], m_pathCount[3], m_pathCount[4]);
	}

	void initPhysics()
	{
		createWorld(btVector3(0, -9.8, 0));
		addGround();
		m_hybrid = new btVariationalHybrid();
		m_newton = new btVariationalProjectedNewton();
		m_newton->m_tolerance = 1e-2;
		m_newton->m_linearTolerance = 0.05;
		m_vbd = new btVariationalVertexBlockDescent();
		m_vbd->m_iterations = 20;
		m_quasiNewton = new btVariationalQuasiNewton();
		m_quasiNewton->m_tolerance = 2.5e-3;
		m_quasiNewton->m_maxIterations = 100;
		m_integrators.push_back(m_hybrid);
		m_integrators.push_back(m_newton);
		m_integrators.push_back(m_vbd);
		m_integrators.push_back(m_quasiNewton);
		m_variationalSolver->setIntegrator(m_hybrid);

		btVariationalSystem& sys = system();
		sys.m_damping = 1.0;
		sys.m_colliderModel = btVariationalSystem::COLLIDER_BARRIER;
		sys.m_colliderMargin = 0.004;
		sys.m_colliderStiffness = 100.0;
		m_ipc.m_dhat = 0.004;
		sys.m_contact = &m_ipc;

		btVariationalStableNeoHookean* soft = new btVariationalStableNeoHookean(2e5, 0.42);
		btVariationalStableNeoHookean* firm = new btVariationalStableNeoHookean(1.5e6, 0.42);
		m_materials.push_back(soft);
		m_materials.push_back(firm);

		// tower of cubes with alternating stiffness
		for (int i = 0; i < 4; i++)
			addTetGrid(btVector3(0, btScalar(0.21) + btScalar(0.42) * i, 0), btVector3(0.2, 0.2, 0.2), 2, 2, 2, 1000, (i % 2) ? soft : firm, 0.5);

		// wrecking ball: three light blocks and a heavy one, hard joints corner to corner, swinging anchor
		{
			const int n = 2;
			const btScalar half = 0.12, stiffness = 5e4;
			const btVector3 top(-1.05, 2.1, 0);
			btSoftBody* previous = 0;
			for (int k = 0; k < 4; k++)
			{
				btSoftBody* block = addTetGrid(top - btVector3(0, half + btScalar(2) * half * k, 0), btVector3(half, half, half), n, n, n, k == 3 ? 15000 : 300, firm, 0.5);
				// the blocks of the chain touch at their joints, so they must not collide with each other
				m_ipc.setCollisionGroup(block, 1);
				for (int ix = 0; ix < 2; ix++)
					for (int iz = 0; iz < 2; iz++)
					{
						if (previous)
							m_variationalSolver->addAttachment(previous, corner(n, ix, 0, iz), block, corner(n, ix, 1, iz), stiffness, true);
						else
						{
							const btVector3 target = block->m_nodes[corner(n, ix, 1, iz)].m_x;
							m_anchorAttachments.push_back(m_variationalSolver->addAttachment(block, corner(n, ix, 1, iz), target, stiffness, true));
							m_anchorRest.push_back(target);
						}
					}
				previous = block;
			}
		}

		// cloth falling over a soft cube
		addTetGrid(btVector3(1.5, 0.21, 0.7), btVector3(0.2, 0.2, 0.2), 2, 2, 2, 1000, soft, 0.5);
		{
			const btScalar s = 0.5, h = 0.75;
			const btVector3 c(1.5, 0, 0.7);
			btSoftBody* cloth = btSoftBodyHelpers::CreatePatch(getDeformableDynamicsWorld()->getWorldInfo(), c + btVector3(-s, h, -s), c + btVector3(s, h, -s),
															   c + btVector3(-s, h, s), c + btVector3(s, h, s), 13, 13, 0, true);
			cloth->generateBendingConstraints(2);
			cloth->setTotalMass(3.0);  // a rubber mat rather than silk: a few kilograms per square metre
			addSoftBody(cloth, 0, 0.5);
			m_variationalSolver->setSpringStiffness(cloth, 4000.0, 20.0);
		}

		// twisted bar
		{
			m_barCentre = btVector3(0, 0.9, -1.3);
			m_bar = addTetGrid(m_barCentre, btVector3(0.7, 0.1, 0.1), 12, 2, 2, 1000, soft, 0.5);
			for (int i = 0; i < m_bar->m_nodes.size(); i++)
			{
				const btVector3 x = m_bar->m_nodes[i].m_x;
				if (x.x() < btScalar(-0.699)) m_bar->setMass(i, 0);
				if (x.x() > btScalar(0.699))
				{
					m_bar->setMass(i, 0);
					m_twistNodes.push_back(i);
					m_twistRest.push_back(x);
				}
			}
		}
		m_guiHelper->autogenerateGraphicsObjects(m_dynamicsWorld);

		if (m_guiHelper->getParameterInterface())
		{
			SliderParams solver("solver: 0 hybrid, 1 Newton, 2 VBD, 3 quasi-Newton", &gCombinedSolver);
			solver.m_minVal = 0;
			solver.m_maxVal = 3;
			solver.m_clampToIntegers = true;
			m_guiHelper->getParameterInterface()->registerSliderFloatParameter(solver);
		}
	}

	virtual void stepSimulation(float deltaTime)
	{
		const int s = int(gCombinedSolver + btScalar(0.5));
		m_variationalSolver->setIntegrator(s == 1 ? (btVariationalIntegrator*)m_newton : (s == 2 ? (btVariationalIntegrator*)m_vbd : (s == 3 ? (btVariationalIntegrator*)m_quasiNewton : (btVariationalIntegrator*)m_hybrid)));
		m_time += m_fixedTimeStep;
		// the anchor of the wrecking ball swings towards the tower after the scene has settled
		const btScalar swing = m_time > btScalar(1.5) ? btScalar(0.5) * (btScalar(1) - btCos(SIMD_2_PI * (m_time - btScalar(1.5)) / btScalar(3.5))) : btScalar(0);
		for (int i = 0; i < m_anchorAttachments.size(); i++)
			m_variationalSolver->setAttachmentTarget(m_anchorAttachments[i], m_anchorRest[i] + btVector3(swing, 0, 0));
		// twist the bar by up to one turn
		const btQuaternion q(btVector3(1, 0, 0), SIMD_2_PI * btSin(SIMD_2_PI * m_time / btScalar(10)));
		for (int i = 0; i < m_twistNodes.size(); i++)
		{
			btSoftBody::Node& n = m_bar->m_nodes[m_twistNodes[i]];
			n.m_x = m_barCentre + quatRotate(q, m_twistRest[i] - m_barCentre);
			n.m_q = n.m_x;
		}
		VariationalDemoBase::stepSimulation(deltaTime);
	}

	virtual void exitPhysics()
	{
		if (m_variationalSolver) system().m_contact = 0;
		VariationalDemoBase::exitPhysics();
	}
};

class CommonExampleInterface* VariationalCombinedCreateFunc(struct CommonExampleOptions& options)
{
	return new VariationalCombined(options.m_guiHelper);
}
