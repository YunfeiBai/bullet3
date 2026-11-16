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

#include "btVariationalDeformableSolver.h"
#include "BulletSoftBody/btDeformableMousePickingForce.h"
#include "BulletCollision/CollisionDispatch/btCollisionObjectWrapper.h"
#include "LinearMath/btQuickprof.h"

btVariationalDeformableSolver::btVariationalDeformableSolver()
	: m_integrator(&m_defaultIntegrator),
	  m_defaultMaterial(1e5, 0.4),
	  m_dirty(true),
	  m_firstPickingAttachment(0),
	  m_lastSolveSeconds(0)
{
}

btVariationalDeformableSolver::BodySettings& btVariationalDeformableSolver::settingsOf(btSoftBody* psb)
{
	for (int i = 0; i < m_settings.size(); i++)
		if (m_settings[i].m_psb == psb) return m_settings[i];
	BodySettings s;
	s.m_psb = psb;
	s.m_material = &m_defaultMaterial;
	s.m_springStiffness = 1000.0;
	s.m_bendingStiffness = 10.0;
	s.m_numNodes = s.m_numTets = s.m_numLinks = s.m_numFaces = -1;
	m_settings.push_back(s);
	return m_settings[m_settings.size() - 1];
}

void btVariationalDeformableSolver::setMaterial(btSoftBody* psb, const btVariationalMaterial* material)
{
	settingsOf(psb).m_material = material ? material : &m_defaultMaterial;
	const int b = m_system.findBody(psb);
	if (b >= 0)
	{
		m_system.m_bodies[b].m_material = settingsOf(psb).m_material;
		m_system.m_topologyVersion++;  // cached factorizations depend on the stiffness
	}
}

void btVariationalDeformableSolver::setSpringStiffness(btSoftBody* psb, btScalar structural, btScalar bending)
{
	BodySettings& s = settingsOf(psb);
	s.m_springStiffness = structural;
	s.m_bendingStiffness = bending;
	m_dirty = true;
}

int btVariationalDeformableSolver::addAttachment(btSoftBody* psb, int node, const btVector3& target, btScalar stiffness, bool hard)
{
	AttachmentSettings a;
	a.m_psbA = psb;
	a.m_nodeA = node;
	a.m_psbB = 0;
	a.m_nodeB = -1;
	a.m_target = target;
	a.m_stiffness = stiffness;
	a.m_hard = hard;
	m_attachmentSettings.push_back(a);
	m_dirty = true;
	return m_attachmentSettings.size() - 1;
}

int btVariationalDeformableSolver::addAttachment(btSoftBody* psbA, int nodeA, btSoftBody* psbB, int nodeB, btScalar stiffness, bool hard)
{
	AttachmentSettings a;
	a.m_psbA = psbA;
	a.m_nodeA = nodeA;
	a.m_psbB = psbB;
	a.m_nodeB = nodeB;
	a.m_target.setZero();
	a.m_stiffness = stiffness;
	a.m_hard = hard;
	m_attachmentSettings.push_back(a);
	m_dirty = true;
	return m_attachmentSettings.size() - 1;
}

void btVariationalDeformableSolver::setAttachmentTarget(int attachment, const btVector3& target)
{
	m_attachmentSettings[attachment].m_target = target;
	if (!m_dirty && attachment < m_system.m_attachments.size()) m_system.m_attachments[attachment].m_target = btVec3d(target);
}

void btVariationalDeformableSolver::removeAllAttachments()
{
	m_attachmentSettings.resize(0);
	m_dirty = true;
}

void btVariationalDeformableSolver::rebuildIfNeeded()
{
	bool rebuild = m_dirty || m_system.m_bodies.size() != m_softBodies.size();
	for (int i = 0; i < m_softBodies.size() && !rebuild; i++)
	{
		btSoftBody* psb = m_softBodies[i];
		const BodySettings& s = settingsOf(psb);
		if (m_system.m_bodies[i].m_psb != psb || s.m_numNodes != psb->m_nodes.size() || s.m_numTets != psb->m_tetras.size() ||
			s.m_numLinks != psb->m_links.size() || s.m_numFaces != psb->m_faces.size())
			rebuild = true;
	}
	if (!rebuild) return;
	m_system.clear();
	for (int i = 0; i < m_softBodies.size(); i++)
	{
		btSoftBody* psb = m_softBodies[i];
		BodySettings& s = settingsOf(psb);
		s.m_numNodes = psb->m_nodes.size();
		s.m_numTets = psb->m_tetras.size();
		s.m_numLinks = psb->m_links.size();
		s.m_numFaces = psb->m_faces.size();
		m_system.addBody(psb, s.m_material, s.m_springStiffness, s.m_bendingStiffness);
		// the solver does not track Bullet's sleeping heuristics
		psb->setActivationState(DISABLE_DEACTIVATION);
	}
	m_system.buildTopology();
	for (int a = 0; a < m_attachmentSettings.size(); a++)
	{
		const AttachmentSettings& src = m_attachmentSettings[a];
		const int bodyA = m_system.findBody(src.m_psbA);
		const int bodyB = src.m_psbB ? m_system.findBody(src.m_psbB) : -1;
		btVariationalAttachment at;
		at.m_nodeA = bodyA >= 0 ? m_system.m_bodies[bodyA].m_firstNode + src.m_nodeA : -1;
		at.m_nodeB = bodyB >= 0 ? m_system.m_bodies[bodyB].m_firstNode + src.m_nodeB : -1;
		at.m_target = btVec3d(src.m_target);
		at.m_stiffness = src.m_stiffness;
		at.m_hard = src.m_hard;
		at.m_penalty = src.m_stiffness;
		if (at.m_nodeA < 0 || (src.m_psbB && at.m_nodeB < 0))
		{
			// refers to a body that is not in the world (any more); keep the slot so indices stay valid
			at.m_nodeA = 0;
			at.m_nodeB = -1;
			at.m_stiffness = at.m_penalty = 0.0;
			at.m_hard = false;
		}
		m_system.m_attachments.push_back(at);
	}
	m_system.rebuildAttachmentLookup();
	m_firstPickingAttachment = m_system.m_attachments.size();
	m_dirty = false;
}

void btVariationalDeformableSolver::reinitialize(const btAlignedObjectArray<btSoftBody*>& softBodies, btScalar dt)
{
	// The world also calls non-virtual methods of the base class (backupVelocity, the island callback,
	// the mouse picking force bookkeeping), so its arrays and node indices must stay valid.
	btDeformableBodySolver::reinitialize(softBodies, dt);
	rebuildIfNeeded();
}

void btVariationalDeformableSolver::predictMotion(btScalar solverdt)
{
	m_pendingColliders.resize(0);
	const btScalar gravityReach = btScalar(m_system.m_gravity.length()) * solverdt * solverdt + btScalar(m_system.m_colliderMargin);
	for (int i = 0; i < m_softBodies.size(); i++)
	{
		btSoftBody* psb = m_softBodies[i];
		if (psb->m_bUpdateRtCst)
		{
			psb->m_bUpdateRtCst = false;
			psb->updateConstants();
			psb->m_fdbvt.clear();
		}
		psb->m_sst.sdt = solverdt * psb->m_cfg.timescale;
		psb->m_sst.isdt = 1 / psb->m_sst.sdt;
		psb->m_sst.velmrg = psb->m_sst.sdt * 3;
		psb->m_sst.radmrg = psb->getCollisionShape()->getMargin();
		psb->m_sst.updmrg = psb->m_sst.radmrg * (btScalar)0.25;
		psb->m_nodeRigidContacts.resize(0);
		psb->m_faceRigidContacts.resize(0);
		psb->m_faceNodeContacts.resize(0);
		psb->m_faceNodeContactsCCD.resize(0);

		// Enlarge the broadphase bounds by what the body can travel in this step, so that every rigid
		// object it may touch is reported to processCollision().
		btScalar maxSpeed = 0;
		for (int j = 0; j < psb->m_nodes.size(); j++)
		{
			btSoftBody::Node& n = psb->m_nodes[j];
			maxSpeed = btMax(maxSpeed, n.m_v.length());
			n.m_q = n.m_x;
			n.m_splitv.setZero();
			n.m_constrained = false;
		}
		psb->updateBounds();
		const btScalar reach = maxSpeed * solverdt + gravityReach;
		psb->m_bounds[0] -= btVector3(reach, reach, reach);
		psb->m_bounds[1] += btVector3(reach, reach, reach);
		if (psb->getBroadphaseHandle())
			psb->getWorldInfo()->m_broadphase->setAabb(psb->getBroadphaseHandle(), psb->m_bounds[0], psb->m_bounds[1], psb->getWorldInfo()->m_dispatcher);
		psb->updateNodeTree(true, true);
	}
}

void btVariationalDeformableSolver::processCollision(btSoftBody* softBody, const btCollisionObjectWrapper* wrap)
{
	const btCollisionObject* object = wrap->getCollisionObject();
	const btCollisionShape* root = object->getCollisionShape();
	PendingCollider pc;
	pc.m_psb = softBody;
	pc.m_object = object;
	if (root->getShapeType() == STATIC_PLANE_PROXYTYPE)
	{
		// planes arrive here as a stream of temporary triangles; use the plane itself
		pc.m_shape = root;
		pc.m_localTransform.setIdentity();
	}
	else if (root->isConcave() || !wrap->getCollisionShape()->isConvex())
	{
		return;  // triangle meshes have no signed distance field here
	}
	else
	{
		pc.m_shape = wrap->getCollisionShape();
		pc.m_localTransform = object->getWorldTransform().inverse() * wrap->getWorldTransform();
	}
	for (int i = 0; i < m_pendingColliders.size(); i++)
		if (m_pendingColliders[i].m_psb == softBody && m_pendingColliders[i].m_object == object && m_pendingColliders[i].m_shape == pc.m_shape) return;
	m_pendingColliders.push_back(pc);
}

void btVariationalDeformableSolver::syncPickingForces()
{
	// The example browser picks deformables with a btDeformableMousePickingForce registered at the
	// world. Mirror those as soft attachments of the three nodes of the picked face.
	bool changed = m_system.m_attachments.size() != m_firstPickingAttachment;
	m_system.m_attachments.resize(m_firstPickingAttachment);
	const btAlignedObjectArray<btDeformableLagrangianForce*>& forces = *getLagrangianForceArray();
	for (int f = 0; f < forces.size(); f++)
	{
		if (forces[f]->getForceType() != BT_MOUSE_PICKING_FORCE) continue;
		const btDeformableMousePickingForce* pick = static_cast<const btDeformableMousePickingForce*>(forces[f]);
		const btSoftBody::Face& face = pick->getFace();
		for (int b = 0; b < m_system.m_bodies.size(); b++)
		{
			const btVariationalBody& body = m_system.m_bodies[b];
			const btSoftBody::Node* base = &body.m_psb->m_nodes[0];
			if (face.m_n[0] < base || face.m_n[0] >= base + body.m_numNodes) continue;
			for (int k = 0; k < 3; k++)
			{
				btVariationalAttachment at;
				at.m_nodeA = body.m_firstNode + int(face.m_n[k] - base);
				at.m_target = btVec3d(pick->getMousePos());
				at.m_stiffness = pick->getElasticStiffness();
				m_system.m_attachments.push_back(at);
			}
			changed = true;
			break;
		}
	}
	if (changed) m_system.rebuildAttachmentLookup();
}

void btVariationalDeformableSolver::solveDeformableConstraints(btScalar solverdt)
{
	BT_PROFILE("btVariationalDeformableSolver::solveDeformableConstraints");
	if (m_system.numNodes() == 0) return;
	btClock clock;

	// rigid colliders found by the broadphase, at their end-of-step poses
	m_system.clearColliders();
	for (int i = 0; i < m_pendingColliders.size(); i++)
	{
		const PendingCollider& pc = m_pendingColliders[i];
		const int b = m_system.findBody(pc.m_psb);
		if (b < 0) continue;
		btTransform objectTransform = pc.m_object->getWorldTransform();
		const btRigidBody* rb = btRigidBody::upcast(pc.m_object);
		if (rb && !rb->isStaticOrKinematicObject())
			btTransformUtil::integrateTransform(pc.m_object->getWorldTransform(), rb->getLinearVelocity(), rb->getAngularVelocity(), solverdt, objectTransform);
		const int c = m_system.addCollider(pc.m_object, pc.m_shape, objectTransform * pc.m_localTransform, pc.m_psb->getWorldInfo());
		if (c < 0) continue;
		m_system.m_colliders[c].m_center = btVec3d(objectTransform.getOrigin());
		if (m_system.m_bodies[b].m_colliders.findLinearSearch(c) == m_system.m_bodies[b].m_colliders.size()) m_system.m_bodies[b].m_colliders.push_back(c);
	}
	for (int b = 0; b < m_system.m_bodies.size(); b++) m_system.m_bodies[b].m_friction = m_system.m_bodies[b].m_psb->m_cfg.kDF;
	for (int a = 0; a < m_attachmentSettings.size() && a < m_system.m_attachments.size(); a++)
		if (m_attachmentSettings[a].m_nodeB < 0) m_system.m_attachments[a].m_target = btVec3d(m_attachmentSettings[a].m_target);
	syncPickingForces();

	m_system.readFromSoftBodies();
	m_system.beginStep(solverdt);
	m_integrator->solve(m_system);
	m_system.endStep();

	if (m_system.m_twoWayCoupling)
	{
		for (int c = 0; c < m_system.m_colliders.size(); c++)
		{
			const btVariationalCollider& col = m_system.m_colliders[c];
			btRigidBody* rb = const_cast<btRigidBody*>(btRigidBody::upcast(col.m_object));
			if (!rb || rb->isStaticOrKinematicObject() || col.m_impulse.length2() == 0.0) continue;
			rb->activate();
			rb->applyCentralImpulse(col.m_impulse.toBt());
			rb->applyTorqueImpulse(col.m_torqueImpulse.toBt());
		}
	}
	m_lastSolveSeconds = clock.getTimeMicroseconds() * 1e-6;
}

void btVariationalDeformableSolver::applyTransforms(btScalar timeStep)
{
	(void)timeStep;
	m_system.writeToSoftBodies();
	for (int i = 0; i < m_softBodies.size(); i++)
	{
		btSoftBody* psb = m_softBodies[i];
		psb->updateNormals();
		psb->updateBounds();
		psb->interpolateRenderMesh();
	}
}
