// Contact between deformable bodies and rigid bodies / multibodies with joint motors, as in a
// robot gripper that grasps a deformable object.

#include <gtest/gtest.h>

#include "btBulletDynamicsCommon.h"
#include "BulletDynamics/Featherstone/btMultiBody.h"
#include "BulletDynamics/Featherstone/btMultiBodyJointMotor.h"
#include "BulletDynamics/Featherstone/btMultiBodyLinkCollider.h"
#include "BulletSoftBody/btDeformableBodySolver.h"
#include "BulletSoftBody/btDeformableGravityForce.h"
#include "BulletSoftBody/btDeformableMultiBodyConstraintSolver.h"
#include "BulletSoftBody/btDeformableMultiBodyDynamicsWorld.h"
#include "BulletSoftBody/btDeformableNeoHookeanForce.h"
#include "BulletSoftBody/btSoftBodyHelpers.h"
#include "BulletSoftBody/btSoftBodyRigidBodyCollisionConfiguration.h"

namespace
{
const btScalar kTimeStep = btScalar(1. / 240.);
const btScalar kBlockSize = btScalar(0.1);

struct DeformableScene
{
	btSoftBodyRigidBodyCollisionConfiguration m_collisionConfiguration;
	btCollisionDispatcher m_dispatcher;
	btDbvtBroadphase m_broadphase;
	btDeformableBodySolver m_deformableBodySolver;
	btDeformableMultiBodyConstraintSolver m_solver;
	btDeformableMultiBodyDynamicsWorld* m_world;
	btAlignedObjectArray<btCollisionShape*> m_shapes;
	btAlignedObjectArray<btDeformableLagrangianForce*> m_forces;
	btRigidBody* m_ground;
	btSoftBody* m_block;
	btMultiBody* m_finger;
	btMultiBodyJointMotor* m_motor;

	DeformableScene()
		: m_dispatcher(&m_collisionConfiguration), m_block(0), m_finger(0), m_motor(0)
	{
		m_solver.setDeformableSolver(&m_deformableBodySolver);
		m_world = new btDeformableMultiBodyDynamicsWorld(&m_dispatcher, &m_broadphase, &m_solver, &m_collisionConfiguration, &m_deformableBodySolver);
		const btVector3 gravity(0, -10, 0);
		m_world->setGravity(gravity);
		m_world->getWorldInfo().m_gravity = gravity;

		// static ground with its top face at y = 0
		btBoxShape* groundShape = new btBoxShape(btVector3(1, 0.1, 1));
		m_shapes.push_back(groundShape);
		btTransform groundTransform;
		groundTransform.setIdentity();
		groundTransform.setOrigin(btVector3(0, -0.1, 0));
		m_ground = new btRigidBody(btRigidBody::btRigidBodyConstructionInfo(0, 0, groundShape));
		m_ground->setWorldTransform(groundTransform);
		m_ground->setFriction(1);
		m_world->addRigidBody(m_ground);
	}

	~DeformableScene()
	{
		if (m_motor)
		{
			m_world->removeMultiBodyConstraint(m_motor);
			delete m_motor;
		}
		if (m_finger)
		{
			for (int i = 0; i < m_finger->getNumLinks(); ++i)
			{
				m_world->removeCollisionObject(m_finger->getLink(i).m_collider);
				delete m_finger->getLink(i).m_collider;
			}
			m_world->removeMultiBody(m_finger);
			delete m_finger;
		}
		if (m_block)
		{
			m_world->removeSoftBody(m_block);
			delete m_block;
		}
		for (int i = 0; i < m_forces.size(); ++i)
		{
			delete m_forces[i];
		}
		m_world->removeRigidBody(m_ground);
		delete m_ground;
		delete m_world;
		for (int i = 0; i < m_shapes.size(); ++i)
		{
			delete m_shapes[i];
		}
	}

	// a Neo-Hookean cube resting on the ground, 3x3x3 hexahedra split into 6 tetrahedra each
	void addBlock()
	{
		const int n = 3;
		btAlignedObjectArray<btVector3> positions;
		for (int k = 0; k <= n; ++k)
			for (int j = 0; j <= n; ++j)
				for (int i = 0; i <= n; ++i)
					positions.push_back(btVector3(i - 0.5 * n, j + 0.02, k - 0.5 * n) * (kBlockSize / n));
		m_block = new btSoftBody(&m_world->getWorldInfo(), positions.size(), &positions[0], 0);
		const int corners[6][4] = {{0, 1, 3, 7}, {0, 1, 5, 7}, {0, 2, 3, 7}, {0, 2, 6, 7}, {0, 4, 5, 7}, {0, 4, 6, 7}};
		for (int k = 0; k < n; ++k)
			for (int j = 0; j < n; ++j)
				for (int i = 0; i < n; ++i)
				{
					int cell[8];
					for (int c = 0; c < 8; ++c)
						cell[c] = (i + (c & 1)) + (n + 1) * ((j + ((c >> 1) & 1)) + (n + 1) * (k + ((c >> 2) & 1)));
					for (int t = 0; t < 6; ++t)
					{
						int v[4];
						for (int c = 0; c < 4; ++c)
							v[c] = cell[corners[t][c]];
						// the rest volume of a tetrahedron is signed, so use a positive orientation
						const btVector3& x0 = positions[v[0]];
						if ((positions[v[1]] - x0).cross(positions[v[2]] - x0).dot(positions[v[3]] - x0) < 0)
							btSwap(v[2], v[3]);
						m_block->appendTetra(v[0], v[1], v[2], v[3]);
						for (int a = 0; a < 4; ++a)
							for (int b = a + 1; b < 4; ++b)
								m_block->appendLink(v[a], v[b], 0, true);
					}
				}
		m_block->initializeDmInverse();
		m_block->m_tetraScratches.resize(m_block->m_tetras.size());
		m_block->m_tetraScratchesTn.resize(m_block->m_tetras.size());
		btSoftBodyHelpers::generateBoundaryFaces(m_block);

		// same settings as used by pybullet for volumetric deformable bodies
		m_block->getCollisionShape()->setMargin(0.002);
		m_block->setTotalMass(0.3);
		m_block->m_cfg.kKHR = 1;
		m_block->m_cfg.kCHR = 1;
		m_block->m_cfg.kDF = 1;
		m_block->m_cfg.collisions = btSoftBody::fCollision::SDF_RD | btSoftBody::fCollision::SDF_RDF | btSoftBody::fCollision::SDF_MDF;
		m_block->setCollisionFlags(0);
		m_block->initializeFaceTree();
		m_world->addSoftBody(m_block);

		m_forces.push_back(new btDeformableGravityForce(m_world->getGravity()));
		m_forces.push_back(new btDeformableNeoHookeanForce(3000, 4500, 0.02));
		for (int i = 0; i < m_forces.size(); ++i)
		{
			m_world->addForce(m_block, m_forces[i]);
		}
	}

	// a fixed-base multibody with a box shaped finger on a vertical prismatic joint (positive joint
	// position moves the finger down), driven by a velocity motor with limited force
	void addFinger(btScalar targetVelocity, btScalar maxForce)
	{
		const btVector3 halfExtents(0.03, 0.01, 0.03);
		const btScalar mass = 0.1;
		btVector3 inertia;
		btBoxShape* shape = new btBoxShape(halfExtents);
		m_shapes.push_back(shape);
		shape->calculateLocalInertia(mass, inertia);

		const btVector3 basePosition(0, 0.3, 0);
		const btScalar fingerCenterY = kBlockSize + 0.002 + halfExtents.y() + 0.005;
		m_finger = new btMultiBody(1, 0, btVector3(0, 0, 0), true, false);
		m_finger->setBasePos(basePosition);
		m_finger->setupPrismatic(0, mass, inertia, -1, btQuaternion::getIdentity(), btVector3(0, -1, 0),
								 btVector3(0, fingerCenterY - basePosition.y(), 0), btVector3(0, 0, 0), true);
		m_finger->finalizeMultiDof();
		m_finger->setLinearDamping(0);
		m_finger->setAngularDamping(0);
		m_world->addMultiBody(m_finger);

		btMultiBodyLinkCollider* collider = new btMultiBodyLinkCollider(m_finger, 0);
		collider->setCollisionShape(shape);
		btTransform tr;
		tr.setIdentity();
		tr.setOrigin(btVector3(0, fingerCenterY, 0));
		collider->setWorldTransform(tr);
		collider->setFriction(1);
		m_world->addCollisionObject(collider, btBroadphaseProxy::DefaultFilter, btBroadphaseProxy::AllFilter);
		m_finger->getLink(0).m_collider = collider;

		m_motor = new btMultiBodyJointMotor(m_finger, 0, 0, targetVelocity, maxForce * kTimeStep);
		m_world->addMultiBodyConstraint(m_motor);
		m_motor->finalizeMultiDof();
	}

	void step()
	{
		m_world->stepSimulation(kTimeStep, 1, kTimeStep);
	}

	btScalar minNodeHeight() const
	{
		btScalar y = SIMD_INFINITY;
		for (int i = 0; i < m_block->m_nodes.size(); ++i)
			y = btMin(y, m_block->m_nodes[i].m_x.y());
		return y;
	}

	btScalar maxNodeHeight() const
	{
		btScalar y = -SIMD_INFINITY;
		for (int i = 0; i < m_block->m_nodes.size(); ++i)
			y = btMax(y, m_block->m_nodes[i].m_x.y());
		return y;
	}

	btScalar maxNodeSpeed() const
	{
		btScalar v = 0;
		for (int i = 0; i < m_block->m_nodes.size(); ++i)
		{
			if (!(m_block->m_nodes[i].m_v.length() < SIMD_INFINITY))
				return SIMD_INFINITY;  // NaN
			v = btMax(v, m_block->m_nodes[i].m_v.length());
		}
		return v;
	}
};
}  // namespace

// A velocity controlled joint motor that pushes a finger onto a deformable block has to apply
// its maximum force once the finger rests on the block. This only works if the joint motor sees
// the impulses of the deformable contacts during the constraint solve.
TEST(DeformableMultiBodyContactTest, JointMotorPressesOnDeformable)
{
	DeformableScene scene;
	scene.addBlock();
	const btScalar maxForce = 4;
	scene.addFinger(0.05, maxForce);

	btScalar motorForce = 0;
	btScalar minFingerPos = SIMD_INFINITY;
	btScalar maxFingerPos = -SIMD_INFINITY;
	int numSamples = 0;
	for (int i = 0; i < 360; ++i)
	{
		scene.step();
		ASSERT_LT(scene.maxNodeSpeed(), 1.0) << "step " << i;
		if (i >= 240)
		{
			motorForce += scene.m_motor->getAppliedImpulse(0) / kTimeStep;
			minFingerPos = btMin(minFingerPos, scene.m_finger->getJointPos(0));
			maxFingerPos = btMax(maxFingerPos, scene.m_finger->getJointPos(0));
			++numSamples;
		}
	}
	motorForce /= numSamples;
	// the motor is saturated: the block stops the finger
	EXPECT_NEAR(motorForce, maxForce, 0.1 * maxForce);
	EXPECT_LT(maxFingerPos - minFingerPos, 0.001);
	// the finger squeezes the block, but does not sink into it
	const btScalar fingerBottom = scene.m_finger->getLink(0).m_collider->getWorldTransform().getOrigin().y() - 0.01;
	EXPECT_GT(fingerBottom, scene.maxNodeHeight() - 0.002);
	EXPECT_LT(scene.maxNodeHeight(), kBlockSize);
}

// A deformable block resting on the ground must neither sink into the ground nor jitter.
TEST(DeformableMultiBodyContactTest, DeformableRestsOnGround)
{
	DeformableScene scene;
	scene.addBlock();
	for (int i = 0; i < 480; ++i)
	{
		scene.step();
		ASSERT_LT(scene.maxNodeSpeed(), 1.0) << "step " << i;
		EXPECT_GT(scene.minNodeHeight(), -0.002) << "step " << i;
		if (i >= 240)
		{
			EXPECT_LT(scene.maxNodeSpeed(), 0.05) << "step " << i;
		}
	}
}

int main(int argc, char** argv)
{
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
