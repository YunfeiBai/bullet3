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

///Benchmark of the integrators of the variational deformable module on three scenes that each stress
///a different part of the problem. The numbers in the module README come from this program.
///
///   usage: Benchmark_btVariationalDeformable [resolution = 2] [tolerance in m/s = 1e-2]
///
///   twist: a clamped bar twisted by a full turn, no contact         (large deformation)
///   pile:  four cubes dropped onto each other and the ground        (barrier contact, friction)
///   chain: four light blocks and a 100 times heavier one on hard
///          attachments, swung sideways                              (constraints, mass ratio)
///
///   usage: Benchmark_btVariationalDeformable filters
///
///   compares the two eigenvalue filters of projected Newton (clamp negative eigenvalues to zero, or
///   take their absolute value) on a start with a large volume change and on a boundary-driven twist.

#include <stdio.h>
#include <stdlib.h>

#include "BulletCollision/CollisionDispatch/btCollisionObject.h"
#include "BulletCollision/CollisionShapes/btStaticPlaneShape.h"
#include "BulletSoftBody/btSoftBody.h"
#include "BulletSoftBody/BulletVariationalDeformable/btVariationalHelpers.h"
#include "BulletSoftBody/BulletVariationalDeformable/btVariationalHybrid.h"
#include "BulletSoftBody/BulletVariationalDeformable/btVariationalIPC.h"
#include "BulletSoftBody/BulletVariationalDeformable/btVariationalProjectedNewton.h"
#include "BulletSoftBody/BulletVariationalDeformable/btVariationalQuasiNewton.h"
#include "BulletSoftBody/BulletVariationalDeformable/btVariationalVertexBlockDescent.h"
#include "LinearMath/btQuickprof.h"

static btSoftBodyWorldInfo gWorldInfo;

enum SceneKind
{
	SCENE_TWIST = 0,
	SCENE_PILE,
	SCENE_CHAIN,
	NUM_SCENES
};

struct Scene
{
	btVariationalSystem m_sys;
	btVariationalIPC m_ipc;
	btVariationalStableNeoHookean m_material;
	btStaticPlaneShape m_plane;
	btCollisionObject m_ground;
	btAlignedObjectArray<btSoftBody*> m_bodies;
	btAlignedObjectArray<int> m_anchors;     // attachments of the chain to the world
	btAlignedObjectArray<btVec3d> m_anchorRest;
	btAlignedObjectArray<int> m_twistNodes;  // kinematic end of the bar
	btAlignedObjectArray<btVector3> m_twistRest;
	int m_kind;
	double m_time;

	Scene(int kind, int resolution) : m_material(1e6, 0.45), m_plane(btVector3(0, 1, 0), 0), m_kind(kind), m_time(0)
	{
		m_ground.setCollisionShape(&m_plane);
		m_ground.setFriction(1.0);
		m_sys.m_gravity = btVec3d(0, -9.8, 0);
		m_sys.m_damping = 1.0;
		if (kind == SCENE_TWIST)
		{
			btSoftBody* bar = btVariationalHelpers::createTetGrid(gWorldInfo, btVector3(0, 2, 0), btVector3(1.0, 0.15, 0.15), 8 * resolution, resolution + 1, resolution + 1, 1000);
			for (int i = 0; i < bar->m_nodes.size(); i++)
			{
				const btVector3 x = bar->m_nodes[i].m_x;
				if (x.x() < btScalar(-0.999)) bar->setMass(i, 0);
				if (x.x() > btScalar(0.999))
				{
					bar->setMass(i, 0);
					m_twistNodes.push_back(i);
					m_twistRest.push_back(x);
				}
			}
			m_bodies.push_back(bar);
		}
		else if (kind == SCENE_PILE)
		{
			for (int c = 0; c < 4; c++)
				m_bodies.push_back(btVariationalHelpers::createTetGrid(gWorldInfo, btVector3(btScalar(0.03) * c, btScalar(0.3) + btScalar(0.6) * c, btScalar(0.02) * c),
																	   btVector3(0.22, 0.22, 0.22), resolution + 1, resolution + 1, resolution + 1, 1000));
		}
		else
		{
			for (int k = 0; k < 5; k++)
				m_bodies.push_back(btVariationalHelpers::createTetGrid(gWorldInfo, btVector3(0, btScalar(3.05) - btScalar(0.3) * k, 0), btVector3(0.15, 0.15, 0.15), 2, 2, 2, k == 4 ? 30000 : 300));
		}
		for (int b = 0; b < m_bodies.size(); b++) m_sys.addBody(m_bodies[b], &m_material, 0, 0);
		m_sys.buildTopology();
		if (kind == SCENE_CHAIN)
		{
			// corner nodes of the 2x2x2 grids: index = ix + 3 * (iy + 3 * iz)
			for (int k = 0; k < 5; k++)
				for (int ix = 0; ix < 2; ix++)
					for (int iz = 0; iz < 2; iz++)
					{
						const int top = 2 * ix + 3 * (2 + 3 * 2 * iz), bottom = 2 * ix + 3 * (0 + 3 * 2 * iz);
						btVariationalAttachment at;
						at.m_hard = true;
						at.m_stiffness = at.m_penalty = 5e4;
						if (k == 0)
						{
							at.m_nodeA = m_sys.m_bodies[0].m_firstNode + top;
							at.m_target = at.m_prevTarget = m_sys.m_x[at.m_nodeA];
							m_anchors.push_back(m_sys.m_attachments.size());
							m_anchorRest.push_back(at.m_target);
						}
						else
						{
							at.m_nodeA = m_sys.m_bodies[k - 1].m_firstNode + bottom;
							at.m_nodeB = m_sys.m_bodies[k].m_firstNode + top;
						}
						m_sys.addAttachment(at);
					}
		}
		if (kind == SCENE_PILE)
		{
			m_sys.m_colliderModel = btVariationalSystem::COLLIDER_BARRIER;
			m_sys.m_colliderMargin = 0.004;
			m_sys.m_colliderStiffness = 100.0;
			btTransform identity;
			identity.setIdentity();
			const int c = m_sys.addCollider(&m_ground, &m_plane, identity, &gWorldInfo);
			for (int b = 0; b < m_sys.m_bodies.size(); b++)
			{
				m_sys.m_bodies[b].m_friction = 0.5;
				m_sys.m_bodies[b].m_colliders.push_back(c);
			}
			m_ipc.m_dhat = 0.004;
			m_sys.m_contact = &m_ipc;
		}
	}

	~Scene()
	{
		m_sys.m_contact = 0;
		for (int i = 0; i < m_bodies.size(); i++) delete m_bodies[i];
	}

	// moves the kinematic parts of the scene to the next time step
	void drive()
	{
		m_time += 1.0 / 60.0;
		if (m_kind == SCENE_TWIST)
		{
			const btVector3 centre(0, 2, 0);
			const btQuaternion q(btVector3(1, 0, 0), SIMD_2_PI * btSin(SIMD_2_PI * btScalar(m_time) / btScalar(8)));
			for (int k = 0; k < m_twistNodes.size(); k++) m_bodies[0]->m_nodes[m_twistNodes[k]].m_x = centre + quatRotate(q, m_twistRest[k] - centre);
			m_sys.readFromSoftBodies();
		}
		if (m_kind == SCENE_CHAIN)
		{
			const double offset = m_time > 1.0 ? 0.35 * sin(SIMD_2_PI * (m_time - 1.0) / 3.0) : 0.0;
			for (int k = 0; k < m_anchors.size(); k++) m_sys.m_attachments[m_anchors[k]].m_target = m_anchorRest[k] + btVec3d(offset, 0, 0);
		}
	}

	double jointError() const
	{
		double e = 0;
		for (int a = 0; a < m_sys.m_attachments.size(); a++)
		{
			const btVariationalAttachment& at = m_sys.m_attachments[a];
			e = btMax(e, (m_sys.m_x[at.m_nodeA] - (at.m_nodeB >= 0 ? m_sys.m_x[at.m_nodeB] : at.m_target)).length());
		}
		return e;
	}
};

// Newton iterations of the two eigenvalue filters, with and without trying the exact Hessian first
static void compareFilters()
{
	const char* names[4] = {"clamp", "clamp, exact Hessian first", "absolute value", "absolute value, exact Hessian first"};

	// 1. A free cube that starts as its rest shape scaled uniformly (plus 2 % noise), one quasi-static step.
	//    Iterations until the step is below 1e-6 m or the energy cannot be decreased any further.
	printf("free cube, start scaled uniformly, one quasi-static step (Newton iterations, at most 1000):\n");
	const double poissonRatios[2] = {0.3, 0.495};
	const double scales[4] = {0.2, 0.5, 2.0, 5.0};
	for (int p = 0; p < 2; p++)
		for (int k = 0; k < 4; k++)
		{
			printf("  Poisson's ratio %.3f, scale %.1f:", poissonRatios[p], scales[k]);
			for (int cfg = 0; cfg < 4; cfg++)
			{
				btSoftBody* psb = btVariationalHelpers::createTetGrid(gWorldInfo, btVector3(0, 0, 0), btVector3(0.5, 0.5, 0.5), 3, 3, 3, 1000);
				btVariationalStableNeoHookean material(1e6, poissonRatios[p]);
				btVariationalSystem sys;
				sys.addBody(psb, &material, 0, 0);
				sys.buildTopology();
				sys.m_gravity = btVec3d(0, 0, 0);
				sys.m_hessianFilter = cfg < 2 ? BT_VARIATIONAL_FILTER_CLAMP : BT_VARIATIONAL_FILTER_ABS;
				unsigned int state = 7;
				for (int i = 0; i < sys.numNodes(); i++)
				{
					btVec3d noise;
					for (int c = 0; c < 3; c++)
					{
						state = state * 1664525u + 1013904223u;
						noise[c] = (state >> 8) / 8388608.0 - 1.0;
					}
					sys.m_x[i] = sys.m_x[i] * scales[k] + noise * (0.02 * scales[k]);
				}
				const double dt = 1000.0;
				btVariationalProjectedNewton newton;
				newton.m_tolerance = 1e-6 / dt;
				newton.m_linearTolerance = 1e-6;
				newton.m_maxIterations = 1000;
				newton.m_projectOnDemand = (cfg % 2) == 1;
				sys.beginStep(dt);
				newton.solve(sys);
				sys.endStep();
				printf("  %s: %d (residual %.0e)%s", names[cfg], newton.m_stats.m_iterations, newton.m_stats.m_finalGradient, cfg < 3 ? "," : "\n");
				delete psb;
			}
		}

	// 2. A bar clamped at both ends whose one end is turned back and forth by one and a half turns, 60 steps per second.
	printf("bar twisted by one and a half turns, 360 steps of 1/60 s, tolerance 1e-3 m/s (Newton iterations per step):\n");
	for (int cfg = 0; cfg < 4; cfg++)
	{
		btSoftBody* psb = btVariationalHelpers::createTetGrid(gWorldInfo, btVector3(0, 2, 0), btVector3(1.0, 0.15, 0.15), 16, 3, 3, 1000);
		btAlignedObjectArray<int> twistNodes;
		btAlignedObjectArray<btVector3> twistRest;
		for (int i = 0; i < psb->m_nodes.size(); i++)
		{
			const btVector3 x = psb->m_nodes[i].m_x;
			if (x.x() < btScalar(-0.999)) psb->setMass(i, 0);
			if (x.x() > btScalar(0.999))
			{
				psb->setMass(i, 0);
				twistNodes.push_back(i);
				twistRest.push_back(x);
			}
		}
		btVariationalStableNeoHookean material(1e6, 0.45);
		btVariationalSystem sys;
		sys.addBody(psb, &material, 0, 0);
		sys.buildTopology();
		sys.m_gravity = btVec3d(0, -9.8, 0);
		sys.m_damping = 1.0;
		sys.m_hessianFilter = cfg < 2 ? BT_VARIATIONAL_FILTER_CLAMP : BT_VARIATIONAL_FILTER_ABS;
		btVariationalProjectedNewton newton;
		newton.m_tolerance = 1e-3;
		newton.m_maxIterations = 60;
		newton.m_projectOnDemand = (cfg % 2) == 1;
		int iterations = 0, unconverged = 0;
		double milliseconds = 0;
		const btVector3 centre(0, 2, 0);
		for (int s = 0; s < 360; s++)
		{
			const btQuaternion q(btVector3(1, 0, 0), btScalar(1.5) * SIMD_2_PI * btSin(SIMD_2_PI * btScalar(s + 1) / btScalar(720)));
			for (int k = 0; k < twistNodes.size(); k++) psb->m_nodes[twistNodes[k]].m_x = centre + quatRotate(q, twistRest[k] - centre);
			sys.readFromSoftBodies();
			sys.beginStep(1.0 / 60.0);
			btClock clock;
			newton.solve(sys);
			milliseconds += clock.getTimeMicroseconds() * 1e-3;
			sys.endStep();
			iterations += newton.m_stats.m_iterations;
			if (!newton.m_stats.m_converged) unconverged++;
		}
		printf("  %-36s: %5.1f iterations, %6.1f ms per step, %d steps stopped at the iteration limit\n", names[cfg], iterations / 360.0, milliseconds / 360.0, unconverged);
		delete psb;
	}
}

int main(int argc, char** argv)
{
	gWorldInfo.m_sparsesdf.Initialize();
	if (argc > 1 && argv[1][0] == 'f')
	{
		compareFilters();
		return 0;
	}
	const int resolution = argc > 1 ? atoi(argv[1]) : 2;
	const double tolerance = argc > 2 ? atof(argv[2]) : 1e-2;
	const char* sceneNames[NUM_SCENES] = {"twist", "pile ", "chain"};
	const int numSteps[NUM_SCENES] = {120, 150, 300};
	const char* solverNames[6] = {"Newton", "quasi-Newton (L-BFGS)", "VBD, 20 sweeps", "VBD, 100 sweeps", "Hybrid", "Hybrid, fast path off"};
	printf("resolution %d, Newton tolerance %g m/s, time step 1/60 s\n", resolution, tolerance);
	for (int kind = 0; kind < NUM_SCENES; kind++)
	{
		for (int solver = 0; solver < 6; solver++)
		{
			Scene scene(kind, resolution);
			btVariationalProjectedNewton newton;
			newton.m_tolerance = tolerance;
			newton.m_linearTolerance = 0.05;
			newton.m_maxIterations = 60;
			btVariationalQuasiNewton quasiNewton;
			quasiNewton.m_tolerance = tolerance * 0.25;  // its step is measured in a stiffer metric
			quasiNewton.m_maxIterations = 300;
			btVariationalVertexBlockDescent vbd;  // with the augmented Lagrangian of AVBD for hard attachments
			vbd.m_iterations = solver == 2 ? 20 : 100;
			btVariationalHybrid hybrid;
			hybrid.m_tolerance = tolerance;
			hybrid.m_useQuasiNewton = solver == 4;
			btVariationalIntegrator* integrators[6] = {&newton, &quasiNewton, &vbd, &vbd, &hybrid, &hybrid};
			btVariationalIntegrator* integrator = integrators[solver];

			int iterations = 0, unconverged = 0, intersectingSteps = 0, paths[5] = {0, 0, 0, 0, 0};
			double milliseconds = 0, residualSum = 0, jointMax = 0, jointEnd = 0;
			btVec3dArray g;
			g.resize(scene.m_sys.numNodes());
			for (int s = 0; s < numSteps[kind]; s++)
			{
				scene.drive();
				scene.m_sys.beginStep(1.0 / 60.0);
				btClock clock;
				integrator->solve(scene.m_sys);
				milliseconds += clock.getTimeMicroseconds() * 1e-3;
				// residual of the step: gradient of the incremental potential at the returned positions
				scene.m_sys.gradient(&scene.m_sys.m_x[0], &g[0]);
				residualSum += scene.m_sys.gradientNorm(&g[0]);
				if (kind == SCENE_PILE && scene.m_ipc.hasIntersections(&scene.m_sys.m_x[0])) intersectingSteps++;
				if (kind == SCENE_CHAIN)
				{
					jointEnd = scene.jointError();
					jointMax = btMax(jointMax, jointEnd);
				}
				scene.m_sys.endStep();
				iterations += integrator->m_stats.m_iterations;
				if (!integrator->m_stats.m_converged) unconverged++;
				if (solver >= 4) paths[hybrid.lastPath()]++;
			}
			printf("%s %5d nodes  %-22s: %7.1f ms/step, %6.1f iterations/step, mean residual %.2e", sceneNames[kind], scene.m_sys.numNodes(), solverNames[solver],
				   milliseconds / numSteps[kind], double(iterations) / numSteps[kind], residualSum / numSteps[kind]);
			if (kind == SCENE_PILE)
			{
				// height of the topmost cube: about 2.0 when the pile stands
				const btVariationalBody& last = scene.m_sys.m_bodies[scene.m_sys.m_bodies.size() - 1];
				double height = 0;
				for (int i = last.m_firstNode; i < last.m_firstNode + last.m_numNodes; i++) height += scene.m_sys.m_x[i][1] / last.m_numNodes;
				printf(", intersecting steps %d, top cube at y = %.2f", intersectingSteps, height);
			}
			if (kind == SCENE_CHAIN) printf(", joint error: max %.2e, final %.2e m", jointMax, jointEnd);
			if (solver >= 4) printf(" [steps by path: quasi-Newton %d, Newton %d, Newton + Laplacian %d, quasi-Newton then Newton %d]", paths[1], paths[2], paths[3], paths[4]);
			printf("\n");
		}
	}
	return 0;
}
