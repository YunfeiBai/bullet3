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

#include "VariationalShowcase.h"
#include "../VariationalSceneHelpers.h"
#include "VariationalMeshExport.h"
#include "BulletSoftBody/BulletVariationalDeformable/btVariationalHybrid.h"
#include "BulletSoftBody/BulletVariationalDeformable/btVariationalIPC.h"

///A larger scene for offline rendering (see render_showcase.py and the offline runner
///ShowcaseMain.cpp); with several thousand nodes it steps far below real time.
///
///  - A cloth (mass-spring, Liu et al. 2013) is tied by its corners between four posts and sags
///    into a hammock.
///  - Soft toys drop into it one after another: a teddy bear, a duck, a bunny, a banana, a boot, a
///    loaf, two tori and a ball, all Stable Neo-Hookean (Smith et al. 2018). The organic shapes
///    are coarse tetrahedral meshes from Bullet's data folder; for rendering, the detailed surface
///    models of the same folder are carried along by those tetrahedra.
///  - A wrecking ball hangs on a chain of thirteen links that are joined by hard attachments (the
///    augmented Lagrangian of Augmented VBD, Giles et al. 2025). It is held back at first, then let
///    go, and swings into the loaded hammock.
///  - Everything touches through the barrier contact of Incremental Potential Contact (Li et al.
///    2020) with friction: cloth against toys, toys against each other, the cloth against itself,
///    the chain against the cloth. Nothing interpenetrates at any time.
///
///The solver is btVariationalHybrid.
namespace
{
VariationalShowcaseOptions gShowcaseOptions;

const btScalar kShowcaseStep = btScalar(1.) / btScalar(120.);
}  // namespace

void VariationalShowcaseSetOptions(const VariationalShowcaseOptions& options)
{
	gShowcaseOptions = options;
}

class VariationalShowcase : public VariationalDemoBase
{
	struct Toy
	{
		btSoftBody* m_psb;
		const btVariationalMaterial* m_material;
		btScalar m_friction;
		btScalar m_spawnTime;
		int m_exportHandle;
		bool m_spawned;
	};
	struct Pin
	{
		int m_node;
		btVector3 m_start, m_end;
	};

	btVariationalIPC m_ipc;
	btVariationalHybrid* m_hybrid;
	VariationalMeshExporter m_exporter;
	btAlignedObjectArray<Toy> m_toys;
	btAlignedObjectArray<Pin> m_pins;
	btSoftBody* m_cloth;
	btSoftBody* m_ball;
	int m_holdAttachment;
	btScalar m_releaseTime;
	btScalar m_time;
	bool m_released;

	bool findDataFile(const char* name, char* path, size_t size) const
	{
		const char* prefixes[] = {gShowcaseOptions.m_dataPath, "data", "../data", "../../data", "../../../data", "../../../../data"};
		for (int k = 0; k < 6; k++)
		{
			if (!prefixes[k]) continue;
			snprintf(path, size, "%s/%s", prefixes[k], name);
			FILE* f = fopen(path, "r");
			if (f)
			{
				fclose(f);
				return true;
			}
		}
		return false;
	}

	// moves a freshly loaded mesh to its place: centred, scaled to the given largest dimension, rotated
	static void place(btSoftBody* psb, btScalar size, const btQuaternion& rotation, const btVector3& position)
	{
		btVector3 lo(SIMD_INFINITY, SIMD_INFINITY, SIMD_INFINITY), hi = -lo;
		for (int i = 0; i < psb->m_nodes.size(); i++)
		{
			lo.setMin(psb->m_nodes[i].m_x);
			hi.setMax(psb->m_nodes[i].m_x);
		}
		const btVector3 extent = hi - lo;
		const btScalar scale = size / btMax(extent.x(), btMax(extent.y(), extent.z()));
		psb->translate(-(lo + hi) * btScalar(0.5));
		psb->scale(btVector3(scale, scale, scale));
		psb->rotate(rotation);
		psb->translate(position);
		psb->initializeDmInverse();
		for (int i = 0; i < psb->m_nodes.size(); i++)
		{
			psb->m_nodes[i].m_v.setZero();
			psb->m_nodes[i].m_q = psb->m_nodes[i].m_x;
		}
	}

	// registers a toy that drops into the scene at the given time
	void addToy(btSoftBody* psb, const char* name, const char* look, btScalar youngsModulus, btScalar density, btScalar spawnTime,
				const btAlignedObjectArray<int>* quads = 0, int exportHandle = -1)
	{
		if (!psb) return;
		btVariationalStableNeoHookean* material = new btVariationalStableNeoHookean(youngsModulus, 0.42);
		m_materials.push_back(material);
		btVariationalHelpers::setLumpedMass(psb, density);
		Toy toy;
		toy.m_psb = psb;
		toy.m_material = material;
		toy.m_friction = 0.5;
		toy.m_spawnTime = spawnTime;
		toy.m_spawned = false;
		toy.m_exportHandle = exportHandle;
		m_ipc.excludePair(psb, psb);  // compact solids: no contact of a toy with itself
		if (m_exporter.isOpen() && exportHandle < 0)
			toy.m_exportHandle = quads ? m_exporter.addBody(psb, name, look, &(*quads)[0], quads->size() / 4, 4, 0, 2) : m_exporter.addBodyFromFaces(psb, name, look, 1);
		m_toys.push_back(toy);
	}

	// a toy from a tetrahedral mesh file; if there is a surface model for it, that is what gets rendered
	void addMeshToy(const char* vtkName, const char* objName, const char* texture, const char* name, const char* look, btScalar size,
					const btQuaternion& rotation, const btVector3& position, btScalar youngsModulus, btScalar density, btScalar spawnTime)
	{
		btSoftBody* psb = loadVtk(vtkName);
		if (!psb) return;
		int handle = -1;
		char path[1024];
		if (m_exporter.isOpen() && objName && findDataFile(objName, path, sizeof(path))) handle = m_exporter.addEmbeddedSurface(psb, name, look, path, texture, 1);
		place(psb, size, rotation, position);
		addToy(psb, name, look, youngsModulus, density, spawnTime, 0, handle);
	}

	btSoftBody* loadVtk(const char* name)
	{
		char path[1024];
		if (!findDataFile(name, path, sizeof(path)))
		{
			printf("VariationalShowcase: %s not found (set the data path); the scene goes without it\n", name);
			return 0;
		}
		return btSoftBodyHelpers::CreateFromVtkFile(getDeformableDynamicsWorld()->getWorldInfo(), path);
	}

	void addCapsule(const btVector3& a, const btVector3& b, btScalar radius, const char* look)
	{
		// (vertical capsules only: that is all the scene needs)
		btCapsuleShape* shape = new btCapsuleShape(radius, (b - a).length());
		m_collisionShapes.push_back(shape);
		btTransform t;
		t.setIdentity();
		t.setOrigin((a + b) * btScalar(0.5));
		createRigidBody(0, t, shape)->setFriction(0.6);
		m_exporter.addProp("capsule", look, a, b, radius);
	}

public:
	VariationalShowcase(struct GUIHelperInterface* helper)
		: VariationalDemoBase(helper), m_hybrid(0), m_cloth(0), m_ball(0), m_holdAttachment(-1), m_releaseTime(4.6), m_time(0), m_released(false)
	{
		m_fixedTimeStep = kShowcaseStep;
	}

	void resetCamera()
	{
		float dist = 5.2;
		float pitch = -22;
		float yaw = 200;
		float targetPos[3] = {-0.2, 1.1, 0};
		m_guiHelper->resetCamera(dist, yaw, pitch, targetPos[0], targetPos[1], targetPos[2]);
	}

	virtual const char* methodName() const { return "Showcase, hybrid"; }

	virtual void updateStatus()
	{
		const bool log = m_logToConsole;
		m_logToConsole = false;  // (the per body report of the base class is too long for this scene)
		VariationalDemoBase::updateStatus();
		m_logToConsole = log;
		const btVariationalSystem& sys = system();
		double jointError = 0;
		for (int a = 0; a < sys.m_attachments.size(); a++)
		{
			const btVariationalAttachment& at = sys.m_attachments[a];
			if (at.m_hard) jointError = btMax(jointError, (sys.m_x[at.m_nodeA] - (at.m_nodeB >= 0 ? sys.m_x[at.m_nodeB] : at.m_target)).length());
		}
		const size_t len = strlen(m_statusText);
		snprintf(m_statusText + len, sizeof(m_statusText) - len, " | t %.2f s, joints %.2f mm", double(m_time), jointError * 1000.0);
		if (log && (m_stepCount % 20) == 0)
		{
			const btVariationalSolveStats& stats = m_variationalSolver->getIntegrator()->m_stats;
			double gap = 0;
			const int active = m_ipc.countActivePairs(&sys.m_x[0], &gap);
			int inverted = 0;
			for (int t = 0; t < sys.m_tets.size(); t++)
				if (sys.deformationGradient(sys.m_tets[t], &sys.m_x[0]).det() <= 0) inverted++;
			printf("[step %4d, t %.3f] %d nodes | %d Newton iterations (%s), %d conjugate gradient iterations, %.0f ms | %d active pairs, gap %.2e | joints %.2e | %d inverted",
				   m_stepCount, double(m_time), sys.numNodes(), stats.m_iterations, stats.m_converged ? "converged" : "limit", stats.m_linearIterations,
				   m_variationalSolver->getLastSolveSeconds() * 1000.0, active, gap, jointError, inverted);
			if ((m_stepCount % 120) == 0) printf(" | intersections: %s", m_ipc.hasIntersections(&sys.m_x[0]) ? "YES" : "none");
			printf("\n");
			fflush(stdout);
		}
	}

	void initPhysics()
	{
		createWorld(btVector3(0, -9.8, 0));
		btSoftBodyWorldInfo& worldInfo = getDeformableDynamicsWorld()->getWorldInfo();
		addGround(1.0);
		m_hybrid = new btVariationalHybrid();
		// Newton stops when no node would move faster than this. The light nodes of the cloth react
		// to forces of a few hundredths of a newton, so a tighter tolerance mostly irons out wrinkles
		// nobody sees (Incremental Potential Contact uses one percent of the scene size per second).
		m_hybrid->m_tolerance = 0.04;
		m_integrators.push_back(m_hybrid);
		m_variationalSolver->setIntegrator(m_hybrid);

		btVariationalSystem& sys = system();
		sys.m_damping = 0.5;
		sys.m_colliderModel = btVariationalSystem::COLLIDER_BARRIER;
		sys.m_colliderMargin = 0.003;
		sys.m_colliderStiffness = 100.0;
		m_ipc.m_dhat = 0.003;
		sys.m_contact = &m_ipc;

		if (gShowcaseOptions.m_exportPath && !m_exporter.open(gShowcaseOptions.m_exportPath, 1.0 / (double(kShowcaseStep) * gShowcaseOptions.m_exportEvery)))
			printf("VariationalShowcase: cannot write to %s\n", gShowcaseOptions.m_exportPath);
		m_exporter.addProp("ground", "floor", btVector3(0, 0, 0), btVector3(0, 1, 0), 0);

		// ---- hammock: a cloth that starts flat and whose corners are pulled in towards the posts
		{
			const int res = gShowcaseOptions.m_clothResolution;
			const btScalar half = 0.70, height = 1.30, pin = 0.60, pinHeight = 1.25, post = 0.86, postHeight = 1.42;
			m_cloth = btSoftBodyHelpers::CreatePatch(worldInfo, btVector3(-half, height, -half), btVector3(half, height, -half),
													 btVector3(-half, height, half), btVector3(half, height, half), res, res, 0, true);
			m_cloth->generateBendingConstraints(2);
			m_cloth->setTotalMass(1.2);
			for (int cz = 0; cz < 2; cz++)
				for (int cx = 0; cx < 2; cx++)
				{
					const btVector3 sign(cx ? 1 : -1, 0, cz ? 1 : -1);
					const btVector3 shift(-sign.x() * (half - pin), pinHeight - height, -sign.z() * (half - pin));
					for (int j = 0; j < 2; j++)
						for (int i = 0; i < 2; i++)
						{
							Pin p;
							p.m_node = (cz ? res - 1 - j : j) * res + (cx ? res - 1 - i : i);
							p.m_start = m_cloth->m_nodes[p.m_node].m_x;
							p.m_end = p.m_start + shift;
							m_cloth->setMass(p.m_node, 0);
							m_pins.push_back(p);
						}
					const btVector3 postTop(sign.x() * post, postHeight, sign.z() * post);
					addCapsule(btVector3(postTop.x(), 0.03, postTop.z()), postTop - btVector3(0, 0.03, 0), 0.03, "wood");
					m_exporter.addProp("rope", "rope", postTop - btVector3(0, 0.05, 0), btVector3(sign.x() * pin, pinHeight, sign.z() * pin), 0.008);
				}
			addSoftBody(m_cloth, 0, 0.5);
			m_variationalSolver->setSpringStiffness(m_cloth, 3e4, 40.0);
			if (m_exporter.isOpen())
			{
				btAlignedObjectArray<int> quads;
				btAlignedObjectArray<btScalar> uv;
				for (int j = 0; j < res; j++)
					for (int i = 0; i < res; i++)
					{
						uv.push_back(btScalar(i) / (res - 1));
						uv.push_back(btScalar(j) / (res - 1));
					}
				for (int j = 0; j + 1 < res; j++)
					for (int i = 0; i + 1 < res; i++)
					{
						quads.push_back(j * res + i);
						quads.push_back((j + 1) * res + i);
						quads.push_back((j + 1) * res + i + 1);
						quads.push_back(j * res + i + 1);
					}
				m_exporter.setVisible(m_exporter.addBody(m_cloth, "cloth", "cloth", &quads[0], quads.size() / 4, 4, &uv[0], 1));
			}
		}

		// ---- wrecking ball on a chain, held back until its release
		{
			const btVector3 pivot(0.1, 3.25, 0);
			const btScalar angle = btRadians(55), linkLength = 0.16, linkRadius = 0.028, ballRadius = 0.2;
			const int numLinks = 13, n = 2, nb = 4;
			const btVector3 along(-btSin(angle), -btCos(angle), 0);
			const btQuaternion tilt(btVector3(0, 0, 1), -angle);  // turns +y towards the pivot
			btVariationalStableNeoHookean* linkMaterial = new btVariationalStableNeoHookean(2e7, 0.4);
			btVariationalStableNeoHookean* ballMaterial = new btVariationalStableNeoHookean(3e6, 0.45);
			m_materials.push_back(linkMaterial);
			m_materials.push_back(ballMaterial);
			const int top = n / 2 + (n + 1) * (n + (n + 1) * (n / 2)), bottom = n / 2 + (n + 1) * (0 + (n + 1) * (n / 2));
			btAlignedObjectArray<int> linkQuads, ballQuads;
			btVariationalHelpers::gridSurfaceQuads(n, n, n, linkQuads);
			btVariationalHelpers::gridSurfaceQuads(nb, nb, nb, ballQuads);
			btSoftBody* previous = 0;
			for (int k = 0; k < numLinks; k++)
			{
				btSoftBody* link = btVariationalHelpers::createTetEllipsoid(worldInfo, pivot + along * (linkLength * (k + btScalar(0.5))),
																			btVector3(linkRadius, linkLength / 2, linkRadius), n, n, n, 800);
				variationalRotateAboutCentre(link, tilt);
				link->initializeDmInverse();
				addSoftBody(link, linkMaterial, 0.4);
				m_ipc.excludePair(link, link);
				if (previous)
				{
					m_variationalSolver->addAttachment(previous, bottom, link, top, 2e5, true);
					m_ipc.excludePair(previous, link);  // they touch at the joint
				}
				else
					m_variationalSolver->addAttachment(link, top, pivot, 2e5, true);
				char name[32];
				snprintf(name, sizeof(name), "link%02d", k);
				if (m_exporter.isOpen()) m_exporter.setVisible(m_exporter.addBody(link, name, "chain", &linkQuads[0], linkQuads.size() / 4, 4, 0, 2));
				previous = link;
			}
			const btVector3 ballCentre = pivot + along * (linkLength * numLinks + ballRadius);
			m_ball = btVariationalHelpers::createTetEllipsoid(worldInfo, ballCentre, btVector3(ballRadius, ballRadius, ballRadius), nb, nb, nb, 1500);
			variationalRotateAboutCentre(m_ball, tilt);
			m_ball->initializeDmInverse();
			addSoftBody(m_ball, ballMaterial, 0.4);
			m_ipc.excludePair(m_ball, m_ball);
			const int ballTop = nb / 2 + (nb + 1) * (nb + (nb + 1) * (nb / 2)), ballCentreNode = nb / 2 + (nb + 1) * (nb / 2 + (nb + 1) * (nb / 2));
			m_variationalSolver->addAttachment(previous, bottom, m_ball, ballTop, 2e5, true);
			m_ipc.excludePair(previous, m_ball);
			m_holdAttachment = m_variationalSolver->addAttachment(m_ball, ballCentreNode, m_ball->m_nodes[ballCentreNode].m_x, 2e5, true);
			if (m_exporter.isOpen()) m_exporter.setVisible(m_exporter.addBody(m_ball, "ball", "steel", &ballQuads[0], ballQuads.size() / 4, 4, 0, 2));
			m_exporter.addProp("pivot", "steel", pivot, pivot + btVector3(0, 0.12, 0), 0.03);
		}

		// ---- toys, dropped in one by one
		{
			const btScalar h = 2.15, interval = 0.33;
			btScalar t = 1.3;
			btAlignedObjectArray<int> quads;
			addMeshToy("teddy.vtk", "teddy.obj", 0, "teddy", "toy_teddy", 0.44, btQuaternion(btVector3(1, 0, 0), 0.5), btVector3(0.0, h, 0.05), 5e4, 250, t);
			t += interval;

			btSoftBody* torus = btVariationalHelpers::createTetTorus(worldInfo, btVector3(0, 0, 0), 0.16, 0.055, 20, 3, 250);
			place(torus, 0.42, btQuaternion(btVector3(1, 0, 0), 0.5), btVector3(-0.2, h, -0.14));
			btVariationalHelpers::torusSurfaceQuads(20, 3, quads);
			addToy(torus, "torus_a", "toy_teal", 1.5e5, 250, t, &quads);
			t += interval;

			addMeshToy("duck.vtk", "duck.obj", "duckCM.png", "duck", "toy_duck", 0.38, btQuaternion(btVector3(0, 1, 0), 0.9), btVector3(0.2, h, 0.12), 1.2e5, 200, t);
			t += interval;

			btSoftBody* ball = btVariationalHelpers::createTetEllipsoid(worldInfo, btVector3(-0.12, h, 0.2), btVector3(0.11, 0.11, 0.11), 4, 4, 4, 250);
			btVariationalHelpers::gridSurfaceQuads(4, 4, 4, quads);
			addToy(ball, "ball_a", "toy_yellow", 1e5, 250, t, &quads);
			t += interval;

			addMeshToy("bunny.vtk", "bunny.obj", 0, "bunny", "toy_bunny", 0.38, btQuaternion(btVector3(1, 0, 0), -1.2), btVector3(0.05, h, -0.18), 6e4, 250, t);
			t += interval;

			addMeshToy("banana.vtk", "banana.obj", "banana.png", "banana", "toy_banana", 0.34, btQuaternion(btVector3(1, 0, 0), 1.2), btVector3(-0.16, h, 0.05), 2e5, 250, t);
			t += interval;

			torus = btVariationalHelpers::createTetTorus(worldInfo, btVector3(0, 0, 0), 0.13, 0.045, 20, 3, 250);
			place(torus, 0.34, btQuaternion(btVector3(0, 0, 1), 1.1), btVector3(0.16, h, -0.1));
			btVariationalHelpers::torusSurfaceQuads(20, 3, quads);
			addToy(torus, "torus_b", "toy_purple", 1.5e5, 250, t, &quads);
			t += interval;

			addMeshToy("boot.vtk", 0, 0, "boot", "toy_red", 0.34, btQuaternion(btVector3(0, 0, 1), 0.4), btVector3(-0.05, h, 0.16), 1.2e5, 250, t);
			t += interval;

			addMeshToy("bread.vtk", 0, 0, "bread", "toy_bread", 0.34, btQuaternion(btVector3(0, 1, 0), 0.8), btVector3(0.12, h, 0.0), 3e4, 150, t);
			t += interval;
			m_releaseTime = t + btScalar(0.4);
		}
		m_guiHelper->autogenerateGraphicsObjects(m_dynamicsWorld);
	}

	virtual void stepSimulation(float deltaTime)
	{
		// the corners of the cloth move in during the first 0.8 s
		const btScalar u = btMin(btScalar(1), (m_time + m_fixedTimeStep) / btScalar(0.8));
		const btScalar s = u * u * (btScalar(3) - btScalar(2) * u);
		for (int i = 0; i < m_pins.size(); i++)
		{
			btSoftBody::Node& n = m_cloth->m_nodes[m_pins[i].m_node];
			n.m_x = m_pins[i].m_start + (m_pins[i].m_end - m_pins[i].m_start) * s;
			n.m_q = n.m_x;
		}
		for (int i = 0; i < m_toys.size(); i++)
		{
			Toy& toy = m_toys[i];
			if (toy.m_spawned || m_time < toy.m_spawnTime) continue;
			addSoftBody(toy.m_psb, toy.m_material, toy.m_friction);
			toy.m_spawned = true;
			if (toy.m_exportHandle >= 0) m_exporter.setVisible(toy.m_exportHandle);
		}
		if (!m_released && m_time >= m_releaseTime)
		{
			m_variationalSolver->setAttachmentEnabled(m_holdAttachment, false);
			m_released = true;
		}
		VariationalDemoBase::stepSimulation(deltaTime);
		m_time += m_fixedTimeStep;
		if (m_exporter.isOpen() && (m_stepCount % gShowcaseOptions.m_exportEvery) == 0) m_exporter.writeFrame();
	}

	virtual void exitPhysics()
	{
		m_exporter.close();
		// toys that never entered the world belong to nobody else
		for (int i = 0; i < m_toys.size(); i++)
			if (!m_toys[i].m_spawned) delete m_toys[i].m_psb;
		m_toys.clear();
		if (m_variationalSolver) system().m_contact = 0;
		VariationalDemoBase::exitPhysics();
	}
};

class CommonExampleInterface* VariationalShowcaseCreateFunc(struct CommonExampleOptions& options)
{
	return new VariationalShowcase(options.m_guiHelper);
}
