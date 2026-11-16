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

#include "btVariationalSystem.h"
#include "BulletSoftBody/btSoftBody.h"
#include "BulletCollision/CollisionShapes/btStaticPlaneShape.h"
#include "BulletCollision/CollisionShapes/btBoxShape.h"
#include "BulletCollision/CollisionShapes/btSphereShape.h"
#include "BulletCollision/CollisionShapes/btCapsuleShape.h"
#include "BulletDynamics/Dynamics/btRigidBody.h"

//
// Hessian storage and conjugate gradients
//

void btVariationalHessian::reset(int numNodes)
{
	m_stencils.resize(0);
	m_nodeBlocks.resize(0);
	m_nodeBlocks.resize(numNodes);
}

void btVariationalHessian::multiply(const btVec3d* x, btVec3d* y, const char* fixed) const
{
	const int n = m_nodeBlocks.size();
	for (int i = 0; i < n; i++) y[i] = fixed[i] ? x[i] : m_nodeBlocks[i] * x[i];
	for (int s = 0; s < m_stencils.size(); s++)
	{
		const btVariationalHessianStencil& st = m_stencils[s];
		double xe[12], ye[12];
		int count = 0;
		for (int j = 0; j < 4; j++)
		{
			const int nj = st.m_n[j];
			if (nj < 0) break;
			count++;
			const bool active = !fixed[nj];
			xe[3 * j + 0] = active ? x[nj][0] : 0.0;
			xe[3 * j + 1] = active ? x[nj][1] : 0.0;
			xe[3 * j + 2] = active ? x[nj][2] : 0.0;
		}
		const int dim = 3 * count;
		for (int r = 0; r < dim; r++)
		{
			const double* row = st.m_H.m + 12 * r;
			double sum = 0;
			for (int c = 0; c < dim; c++) sum += row[c] * xe[c];
			ye[r] = sum;
		}
		for (int j = 0; j < count; j++)
		{
			const int nj = st.m_n[j];
			if (fixed[nj]) continue;
			y[nj][0] += ye[3 * j + 0];
			y[nj][1] += ye[3 * j + 1];
			y[nj][2] += ye[3 * j + 2];
		}
	}
}

void btVariationalHessian::diagonalBlocks(btMat3d* D) const
{
	const int n = m_nodeBlocks.size();
	for (int i = 0; i < n; i++) D[i] = m_nodeBlocks[i];
	for (int s = 0; s < m_stencils.size(); s++)
	{
		const btVariationalHessianStencil& st = m_stencils[s];
		for (int j = 0; j < 4; j++)
		{
			if (st.m_n[j] < 0) break;
			D[st.m_n[j]] += st.m_H.block3(j, j);
		}
	}
}

int btVariationalConjugateGradient(const btVariationalHessian& H, const btVec3d* b, btVec3d* x, const char* fixed, int numNodes,
								   const btVariationalPreconditioner* preconditioner, double relativeTolerance, int maxIterations,
								   bool* indefinite)
{
	if (indefinite) *indefinite = false;
	btVec3dArray r, z, p, Ap;
	r.resize(numNodes);
	z.resize(numNodes);
	p.resize(numNodes);
	Ap.resize(numNodes);
	double bNorm2 = 0;
	for (int i = 0; i < numNodes; i++)
	{
		x[i].setZero();
		r[i] = fixed[i] ? btVec3d(0, 0, 0) : b[i];
		bNorm2 += r[i].length2();
	}
	if (bNorm2 == 0.0) return 0;
	const double stop2 = relativeTolerance * relativeTolerance * bNorm2;
	if (preconditioner)
		preconditioner->apply(&r[0], &z[0]);
	else
		for (int i = 0; i < numNodes; i++) z[i] = r[i];
	double rz = 0;
	for (int i = 0; i < numNodes; i++)
	{
		if (fixed[i]) z[i].setZero();
		p[i] = z[i];
		rz += r[i].dot(z[i]);
	}
	int it = 0;
	for (; it < maxIterations; it++)
	{
		H.multiply(&p[0], &Ap[0], fixed);
		double pAp = 0;
		for (int i = 0; i < numNodes; i++)
			if (!fixed[i]) pAp += p[i].dot(Ap[i]);
		if (!(pAp > 0.0))
		{
			if (indefinite) *indefinite = true;
			break;  // not positive definite along p; keep what we have
		}
		const double alpha = rz / pAp;
		double rNorm2 = 0;
		for (int i = 0; i < numNodes; i++)
		{
			if (fixed[i]) continue;
			x[i] += p[i] * alpha;
			r[i] -= Ap[i] * alpha;
			rNorm2 += r[i].length2();
		}
		if (rNorm2 <= stop2)
		{
			it++;
			break;
		}
		if (preconditioner)
			preconditioner->apply(&r[0], &z[0]);
		else
			for (int i = 0; i < numNodes; i++) z[i] = r[i];
		double rzNew = 0;
		for (int i = 0; i < numNodes; i++)
		{
			if (fixed[i]) z[i].setZero();
			rzNew += r[i].dot(z[i]);
		}
		const double beta = rzNew / rz;
		rz = rzNew;
		for (int i = 0; i < numNodes; i++) p[i] = z[i] + p[i] * beta;
	}
	return it;
}

//
// Colliders
//

void btVariationalCollider::signedDistance(const btVec3d& x, double& dist, btVec3d& normal) const
{
	const btVec3d l = m_R.mulT(x - m_p);
	btVec3d nl(0, 1, 0);
	switch (m_type)
	{
		case PLANE:
		{
			dist = m_planeNormal.dot(l) - m_planeConstant;
			nl = m_planeNormal;
			break;
		}
		case BOX:
		{
			const btVec3d q(fabs(l[0]) - m_halfExtents[0], fabs(l[1]) - m_halfExtents[1], fabs(l[2]) - m_halfExtents[2]);
			if (q[0] <= 0.0 && q[1] <= 0.0 && q[2] <= 0.0)
			{
				int a = q[0] > q[1] ? (q[0] > q[2] ? 0 : 2) : (q[1] > q[2] ? 1 : 2);
				dist = q[a];
				nl = btVec3d(0, 0, 0);
				nl[a] = l[a] >= 0.0 ? 1.0 : -1.0;
			}
			else
			{
				btVec3d o(q[0] > 0.0 ? q[0] : 0.0, q[1] > 0.0 ? q[1] : 0.0, q[2] > 0.0 ? q[2] : 0.0);
				dist = o.length();
				nl = btVec3d(l[0] >= 0.0 ? o[0] : -o[0], l[1] >= 0.0 ? o[1] : -o[1], l[2] >= 0.0 ? o[2] : -o[2]) / dist;
			}
			break;
		}
		case SPHERE:
		{
			const double len = l.length();
			dist = len - m_radius;
			if (len > 1e-12) nl = l / len;
			break;
		}
		case CAPSULE:
		{
			btVec3d c(0, 0, 0);
			c[m_axis] = l[m_axis] < -m_halfHeight ? -m_halfHeight : (l[m_axis] > m_halfHeight ? m_halfHeight : l[m_axis]);
			const btVec3d d = l - c;
			const double len = d.length();
			dist = len - m_radius;
			if (len > 1e-12) nl = d / len;
			break;
		}
		default:
		{
			btVector3 nrm(0, 1, 0);
			dist = m_worldInfo->m_sparsesdf.Evaluate(l.toBt(), m_shape, nrm, btScalar(0));
			nl = btVec3d(nrm).normalized();
			break;
		}
	}
	normal = m_R * nl;
}

//
// System
//

btVariationalSystem::btVariationalSystem()
	: m_gravity(0, -9.8, 0),
	  m_dt(1.0 / 60.0),
	  m_hessianFilter(BT_VARIATIONAL_FILTER_CLAMP),
	  m_damping(0.0),
	  m_colliderModel(COLLIDER_PENALTY),
	  m_colliderStiffness(1000.0),
	  m_colliderMargin(0.01),
	  m_frictionEpsilon(1e-2),
	  m_twoWayCoupling(false),
	  m_contact(0),
	  m_topologyVersion(0)
{
}

void btVariationalSystem::clear()
{
	m_x.clear();
	m_xPrev.clear();
	m_v.clear();
	m_vPrev.clear();
	m_xTilde.clear();
	m_fixedTarget.clear();
	m_mass.clear();
	m_fixed.clear();
	m_nodeBody.clear();
	m_bodies.clear();
	m_tets.clear();
	m_springs.clear();
	m_triangles.clear();
	m_edges.clear();
	m_surfaceNode.clear();
	m_attachments.clear();
	m_colliders.clear();
	m_colliderLag.clear();
	m_colliderLagStart.clear();
	m_topologyVersion++;
}

int btVariationalSystem::addBody(btSoftBody* psb, const btVariationalMaterial* material, double springStiffness, double bendingStiffness)
{
	btVariationalBody body;
	body.m_psb = psb;
	body.m_material = material;
	body.m_springStiffness = springStiffness;
	body.m_bendingStiffness = bendingStiffness;
	body.m_friction = psb->m_cfg.kDF;
	body.m_firstNode = body.m_numNodes = 0;
	body.m_firstTet = body.m_numTets = 0;
	body.m_firstSpring = body.m_numSprings = 0;
	m_bodies.push_back(body);
	return m_bodies.size() - 1;
}

int btVariationalSystem::findBody(const btSoftBody* psb) const
{
	for (int b = 0; b < m_bodies.size(); b++)
		if (m_bodies[b].m_psb == psb) return b;
	return -1;
}

namespace
{
struct btSortedFace
{
	int m_key[3];   // sorted node indices
	int m_orig[3];  // original winding
	bool operator<(const btSortedFace& o) const
	{
		if (m_key[0] != o.m_key[0]) return m_key[0] < o.m_key[0];
		if (m_key[1] != o.m_key[1]) return m_key[1] < o.m_key[1];
		return m_key[2] < o.m_key[2];
	}
	bool sameAs(const btSortedFace& o) const { return m_key[0] == o.m_key[0] && m_key[1] == o.m_key[1] && m_key[2] == o.m_key[2]; }
};
struct btSortedEdge
{
	int m_a, m_b;
	bool operator<(const btSortedEdge& o) const { return m_a != o.m_a ? m_a < o.m_a : m_b < o.m_b; }
};
template <class T>
struct btLess
{
	bool operator()(const T& a, const T& b) const { return a < b; }
};

void btBuildNodeLookup(int numNodes, int numEntries, const int* nodeOfEntry, const int* valueOfEntry, btAlignedObjectArray<int>& start, btAlignedObjectArray<int>& list)
{
	start.resize(0);
	start.resize(numNodes + 1, 0);
	for (int e = 0; e < numEntries; e++) start[nodeOfEntry[e] + 1]++;
	for (int i = 0; i < numNodes; i++) start[i + 1] += start[i];
	list.resize(numEntries);
	btAlignedObjectArray<int> fill;
	fill.resize(numNodes);
	for (int i = 0; i < numNodes; i++) fill[i] = start[i];
	for (int e = 0; e < numEntries; e++) list[fill[nodeOfEntry[e]]++] = valueOfEntry[e];
}
}  // namespace

void btVariationalSystem::buildTopology()
{
	m_tets.resize(0);
	m_springs.resize(0);
	m_triangles.resize(0);
	m_edges.resize(0);
	int numNodes = 0;
	for (int b = 0; b < m_bodies.size(); b++)
	{
		m_bodies[b].m_firstNode = numNodes;
		m_bodies[b].m_numNodes = m_bodies[b].m_psb->m_nodes.size();
		numNodes += m_bodies[b].m_numNodes;
	}
	m_x.resize(numNodes);
	m_xPrev.resize(numNodes);
	m_v.resize(numNodes);
	m_vPrev.resize(numNodes);
	m_xTilde.resize(numNodes);
	m_fixedTarget.resize(numNodes);
	m_mass.resize(numNodes);
	m_fixed.resize(numNodes);
	m_nodeBody.resize(numNodes);
	m_surfaceNode.resize(0);
	m_surfaceNode.resize(numNodes, 0);

	for (int b = 0; b < m_bodies.size(); b++)
	{
		btVariationalBody& body = m_bodies[b];
		btSoftBody* psb = body.m_psb;
		const int first = body.m_firstNode;
		const btSoftBody::Node* base = &psb->m_nodes[0];
		for (int i = 0; i < body.m_numNodes; i++)
		{
			const btSoftBody::Node& n = psb->m_nodes[i];
			m_x[first + i] = btVec3d(n.m_x);
			m_xPrev[first + i] = m_x[first + i];
			m_fixedTarget[first + i] = m_x[first + i];
			m_v[first + i] = btVec3d(n.m_v);
			m_vPrev[first + i] = m_v[first + i];
			m_mass[first + i] = n.m_im > 0 ? 1.0 / n.m_im : 0.0;
			m_fixed[first + i] = n.m_im > 0 ? 0 : 1;
			m_nodeBody[first + i] = b;
		}

		// tetrahedra. Bullet stores the rest shape as Dm^-1 in btScalar precision. If the body is still
		// in its rest pose the rest shape is recomputed from the node positions in double precision.
		body.m_firstTet = m_tets.size();
		bool atRest = true;
		for (int t = 0; t < psb->m_tetras.size() && atRest; t++)
		{
			const btSoftBody::Tetra& tet = psb->m_tetras[t];
			const btVector3 c1 = tet.m_n[1]->m_x - tet.m_n[0]->m_x, c2 = tet.m_n[2]->m_x - tet.m_n[0]->m_x, c3 = tet.m_n[3]->m_x - tet.m_n[0]->m_x;
			const btMatrix3x3 Ds(c1.x(), c2.x(), c3.x(), c1.y(), c2.y(), c3.y(), c1.z(), c2.z(), c3.z());
			const btMatrix3x3 F = Ds * tet.m_Dm_inverse;
			for (int r = 0; r < 3; r++)
				for (int c = 0; c < 3; c++)
					if (!(btFabs(F[r][c] - (r == c ? btScalar(1) : btScalar(0))) < btScalar(1e-3))) atRest = false;
		}
		for (int t = 0; t < psb->m_tetras.size(); t++)
		{
			const btSoftBody::Tetra& src = psb->m_tetras[t];
			btVariationalTet tet;
			tet.m_body = b;
			for (int j = 0; j < 4; j++) tet.m_n[j] = first + int(src.m_n[j] - base);
			btMat3d DmInv;
			if (atRest)
			{
				const btMat3d Dm = btMat3d::fromColumns(m_x[tet.m_n[1]] - m_x[tet.m_n[0]], m_x[tet.m_n[2]] - m_x[tet.m_n[0]], m_x[tet.m_n[3]] - m_x[tet.m_n[0]]);
				const double det = Dm.det();
				if (fabs(det) < 1e-30) continue;  // degenerate element, skip
				DmInv = Dm.inverse();
				tet.m_restVolume = fabs(det) / 6.0;
			}
			else
			{
				for (int r = 0; r < 3; r++)
					for (int c = 0; c < 3; c++) DmInv.m[r][c] = src.m_Dm_inverse[r][c];
				tet.m_restVolume = fabs(double(src.m_element_measure));
				if (!(tet.m_restVolume > 0.0)) continue;
			}
			// F = Ds * DmInv = sum_{j=1..3} (x_j - x_0) * row_{j-1}(DmInv)
			tet.m_g[1] = DmInv.row(0);
			tet.m_g[2] = DmInv.row(1);
			tet.m_g[3] = DmInv.row(2);
			tet.m_g[0] = -(tet.m_g[1] + tet.m_g[2] + tet.m_g[3]);
			m_tets.push_back(tet);
		}
		body.m_numTets = m_tets.size() - body.m_firstTet;

		// bodies without tetrahedra (cloth, ropes) are mass-spring systems built from the links
		body.m_firstSpring = m_springs.size();
		if (body.m_numTets == 0)
		{
			for (int l = 0; l < psb->m_links.size(); l++)
			{
				const btSoftBody::Link& link = psb->m_links[l];
				btVariationalSpring s;
				s.m_body = b;
				s.m_n[0] = first + int(link.m_n[0] - base);
				s.m_n[1] = first + int(link.m_n[1] - base);
				s.m_restLength = link.m_rl;
				s.m_stiffness = link.m_bbending ? body.m_bendingStiffness : body.m_springStiffness;
				if (s.m_stiffness > 0.0 && s.m_restLength > 0.0) m_springs.push_back(s);
			}
		}
		body.m_numSprings = m_springs.size() - body.m_firstSpring;

		// collision surface: boundary faces of the tetrahedral mesh, or the faces of a surface mesh
		const int firstTriangle = m_triangles.size();
		if (body.m_numTets > 0)
		{
			btAlignedObjectArray<btSortedFace> faces;
			faces.reserve(4 * body.m_numTets);
			static const int tetFaces[4][3] = {{1, 2, 3}, {0, 3, 2}, {0, 1, 3}, {0, 2, 1}};
			for (int t = body.m_firstTet; t < body.m_firstTet + body.m_numTets; t++)
			{
				for (int f = 0; f < 4; f++)
				{
					btSortedFace sf;
					for (int k = 0; k < 3; k++) sf.m_orig[k] = sf.m_key[k] = m_tets[t].m_n[tetFaces[f][k]];
					for (int i = 0; i < 2; i++)
						for (int j = 0; j < 2 - i; j++)
							if (sf.m_key[j] > sf.m_key[j + 1])
							{
								int tmp = sf.m_key[j];
								sf.m_key[j] = sf.m_key[j + 1];
								sf.m_key[j + 1] = tmp;
							}
					faces.push_back(sf);
				}
			}
			faces.quickSort(btLess<btSortedFace>());
			for (int i = 0; i < faces.size();)
			{
				int j = i + 1;
				while (j < faces.size() && faces[j].sameAs(faces[i])) j++;
				if (j - i == 1)
				{
					btVariationalTriangle tri;
					tri.m_body = b;
					for (int k = 0; k < 3; k++) tri.m_n[k] = faces[i].m_orig[k];
					m_triangles.push_back(tri);
				}
				i = j;
			}
		}
		else
		{
			for (int f = 0; f < psb->m_faces.size(); f++)
			{
				btVariationalTriangle tri;
				tri.m_body = b;
				for (int k = 0; k < 3; k++) tri.m_n[k] = first + int(psb->m_faces[f].m_n[k] - base);
				m_triangles.push_back(tri);
			}
		}
		btAlignedObjectArray<btSortedEdge> edges;
		for (int t = firstTriangle; t < m_triangles.size(); t++)
		{
			for (int k = 0; k < 3; k++)
			{
				btSortedEdge e;
				e.m_a = m_triangles[t].m_n[k];
				e.m_b = m_triangles[t].m_n[(k + 1) % 3];
				if (e.m_a > e.m_b)
				{
					int tmp = e.m_a;
					e.m_a = e.m_b;
					e.m_b = tmp;
				}
				edges.push_back(e);
				m_surfaceNode[m_triangles[t].m_n[k]] = 1;
			}
		}
		edges.quickSort(btLess<btSortedEdge>());
		for (int i = 0; i < edges.size(); i++)
		{
			if (i > 0 && edges[i].m_a == edges[i - 1].m_a && edges[i].m_b == edges[i - 1].m_b) continue;
			btVariationalEdge e;
			e.m_body = b;
			e.m_n[0] = edges[i].m_a;
			e.m_n[1] = edges[i].m_b;
			m_edges.push_back(e);
		}
		// spring-only bodies without faces (ropes): every node can touch a collider
		if (body.m_numTets == 0 && psb->m_faces.size() == 0)
			for (int i = 0; i < body.m_numNodes; i++) m_surfaceNode[first + i] = 1;
	}

	// node -> element lookups
	{
		btAlignedObjectArray<int> nodeOf, valueOf;
		for (int t = 0; t < m_tets.size(); t++)
			for (int j = 0; j < 4; j++)
			{
				nodeOf.push_back(m_tets[t].m_n[j]);
				valueOf.push_back(4 * t + j);
			}
		btBuildNodeLookup(numNodes, nodeOf.size(), nodeOf.size() ? &nodeOf[0] : 0, valueOf.size() ? &valueOf[0] : 0, m_nodeTetStart, m_nodeTets);
		nodeOf.resize(0);
		valueOf.resize(0);
		for (int s = 0; s < m_springs.size(); s++)
			for (int j = 0; j < 2; j++)
			{
				nodeOf.push_back(m_springs[s].m_n[j]);
				valueOf.push_back(2 * s + j);
			}
		btBuildNodeLookup(numNodes, nodeOf.size(), nodeOf.size() ? &nodeOf[0] : 0, valueOf.size() ? &valueOf[0] : 0, m_nodeSpringStart, m_nodeSprings);
	}
	m_attachments.resize(0);
	rebuildAttachmentLookup();
	m_colliderLag.resize(0);
	m_topologyVersion++;
}

int btVariationalSystem::addAttachment(const btVariationalAttachment& attachment)
{
	m_attachments.push_back(attachment);
	rebuildAttachmentLookup();
	return m_attachments.size() - 1;
}

void btVariationalSystem::clearAttachments()
{
	m_attachments.resize(0);
	rebuildAttachmentLookup();
}

void btVariationalSystem::rebuildAttachmentLookup()
{
	btAlignedObjectArray<int> nodeOf, valueOf;
	for (int a = 0; a < m_attachments.size(); a++)
	{
		nodeOf.push_back(m_attachments[a].m_nodeA);
		valueOf.push_back(2 * a);
		if (m_attachments[a].m_nodeB >= 0)
		{
			nodeOf.push_back(m_attachments[a].m_nodeB);
			valueOf.push_back(2 * a + 1);
		}
	}
	btBuildNodeLookup(numNodes(), nodeOf.size(), nodeOf.size() ? &nodeOf[0] : 0, valueOf.size() ? &valueOf[0] : 0, m_nodeAttachStart, m_nodeAttach);
}

void btVariationalSystem::readFromSoftBodies()
{
	for (int b = 0; b < m_bodies.size(); b++)
	{
		const btVariationalBody& body = m_bodies[b];
		for (int i = 0; i < body.m_numNodes; i++)
		{
			const btSoftBody::Node& n = body.m_psb->m_nodes[i];
			const int g = body.m_firstNode + i;
			m_mass[g] = n.m_im > 0 ? 1.0 / n.m_im : 0.0;
			m_fixed[g] = n.m_im > 0 ? 0 : 1;
			m_fixedTarget[g] = btVec3d(n.m_x);
			if (m_fixed[g]) continue;
			// The solver state is authoritative (and more precise than btScalar). Only take over values
			// that were changed from the outside since the last write-back.
			const btVector3 px = m_x[g].toBt(), pv = m_v[g].toBt();
			if (px.x() != n.m_x.x() || px.y() != n.m_x.y() || px.z() != n.m_x.z()) m_x[g] = btVec3d(n.m_x);
			if (pv.x() != n.m_v.x() || pv.y() != n.m_v.y() || pv.z() != n.m_v.z()) m_v[g] = btVec3d(n.m_v);
		}
	}
}

void btVariationalSystem::writeToSoftBodies()
{
	for (int b = 0; b < m_bodies.size(); b++)
	{
		const btVariationalBody& body = m_bodies[b];
		for (int i = 0; i < body.m_numNodes; i++)
		{
			btSoftBody::Node& n = body.m_psb->m_nodes[i];
			const int g = body.m_firstNode + i;
			n.m_x = m_x[g].toBt();
			n.m_q = n.m_x;
			n.m_v = m_v[g].toBt();
			n.m_vn = n.m_v;
		}
	}
}

void btVariationalSystem::clearColliders()
{
	m_colliders.resize(0);
	for (int b = 0; b < m_bodies.size(); b++) m_bodies[b].m_colliders.resize(0);
}

int btVariationalSystem::addCollider(const btCollisionObject* object, const btCollisionShape* shape, const btTransform& endOfStepTransform, btSoftBodyWorldInfo* worldInfo)
{
	for (int c = 0; c < m_colliders.size(); c++)
		if (m_colliders[c].m_object == object && m_colliders[c].m_shape == shape) return c;
	btVariationalCollider col;
	col.m_object = object;
	col.m_shape = shape;
	col.m_worldInfo = worldInfo;
	col.m_radius = col.m_halfHeight = 0;
	col.m_axis = 1;
	col.m_planeConstant = 0;
	switch (shape->getShapeType())
	{
		case STATIC_PLANE_PROXYTYPE:
		{
			const btStaticPlaneShape* plane = static_cast<const btStaticPlaneShape*>(shape);
			col.m_type = btVariationalCollider::PLANE;
			col.m_planeNormal = btVec3d(plane->getPlaneNormal());
			col.m_planeConstant = plane->getPlaneConstant();
			break;
		}
		case BOX_SHAPE_PROXYTYPE:
		{
			col.m_type = btVariationalCollider::BOX;
			col.m_halfExtents = btVec3d(static_cast<const btBoxShape*>(shape)->getHalfExtentsWithMargin());
			break;
		}
		case SPHERE_SHAPE_PROXYTYPE:
		{
			col.m_type = btVariationalCollider::SPHERE;
			col.m_radius = static_cast<const btSphereShape*>(shape)->getRadius();
			break;
		}
		case CAPSULE_SHAPE_PROXYTYPE:
		{
			const btCapsuleShape* capsule = static_cast<const btCapsuleShape*>(shape);
			col.m_type = btVariationalCollider::CAPSULE;
			col.m_radius = capsule->getRadius();
			col.m_halfHeight = capsule->getHalfHeight();
			col.m_axis = capsule->getUpAxis();
			break;
		}
		default:
		{
			if (!shape->isConvex()) return -1;  // concave shapes have no signed distance here
			col.m_type = btVariationalCollider::GENERIC;
			break;
		}
	}
	const btMatrix3x3& R = endOfStepTransform.getBasis();
	for (int r = 0; r < 3; r++)
		for (int c = 0; c < 3; c++) col.m_R.m[r][c] = R[r][c];
	col.m_p = btVec3d(endOfStepTransform.getOrigin());
	col.m_center = col.m_p;
	const btRigidBody* rb = btRigidBody::upcast(object);
	if (rb)
	{
		col.m_linVel = btVec3d(rb->getLinearVelocity());
		col.m_angVel = btVec3d(rb->getAngularVelocity());
	}
	else
	{
		col.m_linVel = btVec3d(object->getInterpolationLinearVelocity());
		col.m_angVel = btVec3d(object->getInterpolationAngularVelocity());
	}
	col.m_friction = object->getFriction();
	m_colliders.push_back(col);
	return m_colliders.size() - 1;
}

void btVariationalSystem::beginStep(double dt)
{
	m_dt = dt;
	const int n = numNodes();
	for (int i = 0; i < n; i++)
	{
		m_xPrev[i] = m_x[i];
		if (m_fixed[i])
			m_xTilde[i] = m_fixedTarget[i];
		else
			m_xTilde[i] = m_x[i] + m_v[i] * dt + m_gravity * (dt * dt);
	}
	for (int c = 0; c < m_colliders.size(); c++)
	{
		m_colliders[c].m_impulse.setZero();
		m_colliders[c].m_torqueImpulse.setZero();
	}
	// (re)allocate the friction lag slots for the collider pairs of this step
	m_colliderLagStart.resize(m_bodies.size());
	int total = 0;
	for (int b = 0; b < m_bodies.size(); b++)
	{
		m_colliderLagStart[b] = total;
		total += m_bodies[b].m_colliders.size() * m_bodies[b].m_numNodes;
	}
	m_colliderLag.resize(0);
	m_colliderLag.resize(total);
	updateFrictionLag(&m_xPrev[0]);
	if (m_contact) m_contact->beginStep(*this);
}

void btVariationalSystem::endStep()
{
	const int n = numNodes();
	const double invDt = 1.0 / m_dt;
	for (int i = 0; i < n; i++)
	{
		m_vPrev[i] = m_v[i];
		m_v[i] = (m_x[i] - m_xPrev[i]) * invDt;
	}
	if (m_damping > 0.0)
	{
		const double keep = exp(-m_damping * m_dt);
		for (int b = 0; b < m_bodies.size(); b++)
		{
			const btVariationalBody& body = m_bodies[b];
			double mass = 0;
			btVec3d com(0, 0, 0), mom(0, 0, 0);
			for (int i = body.m_firstNode; i < body.m_firstNode + body.m_numNodes; i++)
			{
				mass += m_mass[i];
				com += m_x[i] * m_mass[i];
				mom += m_v[i] * m_mass[i];
			}
			if (!(mass > 0.0)) continue;
			com *= 1.0 / mass;
			const btVec3d vcom = mom * (1.0 / mass);
			btVec3d L(0, 0, 0);
			btMat3d I;
			for (int i = body.m_firstNode; i < body.m_firstNode + body.m_numNodes; i++)
			{
				const btVec3d r = m_x[i] - com;
				L += r.cross(m_v[i]) * m_mass[i];
				I.addScaledIdentity(m_mass[i] * r.length2());
				I.addScaledOuter(-m_mass[i], r, r);
			}
			btVec3d omega(0, 0, 0);
			if (!btSolveSpd3d(I, L, omega)) omega.setZero();
			for (int i = body.m_firstNode; i < body.m_firstNode + body.m_numNodes; i++)
			{
				if (m_fixed[i]) continue;
				const btVec3d rigid = vcom + omega.cross(m_x[i] - com);
				m_v[i] = rigid + (m_v[i] - rigid) * keep;
			}
		}
	}
	if (m_twoWayCoupling)
	{
		// reaction of the contact forces at the end-of-step state, as impulses on the colliders
		for (int b = 0; b < m_bodies.size(); b++)
		{
			const btVariationalBody& body = m_bodies[b];
			for (int s = 0; s < body.m_colliders.size(); s++)
			{
				btVariationalCollider& col = m_colliders[body.m_colliders[s]];
				for (int i = body.m_firstNode; i < body.m_firstNode + body.m_numNodes; i++)
				{
					if (m_fixed[i]) continue;
					double d;
					btVec3d nrm;
					col.signedDistance(m_x[i], d, nrm);
					double dE, d2E;
					colliderKernel(d, m_colliderStiffness * m_mass[i] / (m_dt * m_dt), dE, d2E);
					if (dE == 0.0) continue;
					const btVec3d impulse = nrm * (dE * m_dt);  // minus the impulse on the node
					col.m_impulse += impulse;
					col.m_torqueImpulse += (m_x[i] - col.m_center).cross(impulse);
				}
			}
		}
	}
}

double btVariationalSystem::colliderKernel(double d, double stiffness, double& dE, double& d2E) const
{
	dE = 0.0;
	d2E = 0.0;
	const double m = m_colliderMargin;
	if (d >= m) return 0.0;
	if (m_colliderModel == COLLIDER_PENALTY)
	{
		dE = -stiffness * (m - d);
		d2E = stiffness;
		return 0.5 * stiffness * (m - d) * (m - d);
	}
	// log barrier of IPC written in the distance
	if (d <= 0.0) return 1e300;
	const double t = d - m, lg = log(d / m);
	dE = -stiffness * (2.0 * t * lg + t * t / d);
	d2E = -stiffness * (2.0 * lg + 4.0 * t / d - t * t / (d * d));
	return -stiffness * t * t * lg;
}

void btVariationalSystem::addColliderNode(int b, int node, const btVec3d& xi, const btVec3d& xPrev, double& energy, btVec3d* grad, btMat3d* hess) const
{
	const btVariationalBody& body = m_bodies[b];
	const double stiffness = m_colliderStiffness * m_mass[node] / (m_dt * m_dt);
	for (int s = 0; s < body.m_colliders.size(); s++)
	{
		const btVariationalCollider& col = m_colliders[body.m_colliders[s]];
		double d;
		btVec3d nrm;
		col.signedDistance(xi, d, nrm);
		double dE, d2E;
		energy += colliderKernel(d, stiffness, dE, d2E);
		if (dE != 0.0)
		{
			if (grad) *grad += nrm * dE;
			if (hess) hess->addScaledOuter(d2E, nrm, nrm);
		}
		// friction with the lagged normal force (smoothed Coulomb model of IPC)
		const ColliderLag* lag = lagOf(b, s, node);
		const double mu = body.m_friction * col.m_friction;
		if (!lag || !(lag->m_normalForce > 0.0) || !(mu > 0.0)) continue;
		const btVec3d u = (xi - xPrev) - col.velocityAt(xPrev) * m_dt;
		const btVec3d uT = u - lag->m_normal * u.dot(lag->m_normal);
		const double y = uT.length();
		const double eps = m_frictionEpsilon * m_dt;
		const double scale = mu * lag->m_normalForce;
		if (y >= eps)
		{
			energy += scale * y;
			if (grad) *grad += uT * (scale / y);
			if (hess)
			{
				// (1/y) (P_T - t t^T), with P_T the tangent projector
				hess->addScaledIdentity(scale / y);
				hess->addScaledOuter(-scale / y, lag->m_normal, lag->m_normal);
				hess->addScaledOuter(-scale / (y * y * y), uT, uT);
			}
		}
		else
		{
			energy += scale * (-y * y * y / (3.0 * eps * eps) + y * y / eps + eps / 3.0);
			const double f1OverY = -y / (eps * eps) + 2.0 / eps;
			if (grad) *grad += uT * (scale * f1OverY);
			if (hess)
			{
				// f1/y (P_T - t t^T) + f1' t t^T  with f1' = f1/y - y/eps^2
				hess->addScaledIdentity(scale * f1OverY);
				hess->addScaledOuter(-scale * f1OverY, lag->m_normal, lag->m_normal);
				if (y > 0.0) hess->addScaledOuter(-scale / (eps * eps * y), uT, uT);
			}
		}
	}
}

void btVariationalSystem::updateFrictionLag(const btVec3d* x)
{
	if (m_colliderLag.size() == 0) return;
	for (int b = 0; b < m_bodies.size(); b++)
	{
		const btVariationalBody& body = m_bodies[b];
		for (int s = 0; s < body.m_colliders.size(); s++)
		{
			const btVariationalCollider& col = m_colliders[body.m_colliders[s]];
			ColliderLag* lag = &m_colliderLag[m_colliderLagStart[b] + s * body.m_numNodes];
			for (int i = 0; i < body.m_numNodes; i++)
			{
				const int g = body.m_firstNode + i;
				lag[i].m_normalForce = 0.0;
				if (m_fixed[g]) continue;
				double d, dE, d2E;
				col.signedDistance(x[g], d, lag[i].m_normal);
				colliderKernel(d, m_colliderStiffness * m_mass[g] / (m_dt * m_dt), dE, d2E);
				lag[i].m_normalForce = -dE;
			}
		}
	}
}

double btVariationalSystem::energy(const btVec3d* x) const
{
	const int n = numNodes();
	const double invH2 = 1.0 / (m_dt * m_dt);
	double E = 0;
	for (int i = 0; i < n; i++)
		if (!m_fixed[i]) E += 0.5 * invH2 * m_mass[i] * (x[i] - m_xTilde[i]).length2();
	E += elasticEnergy(x);
	for (int a = 0; a < m_attachments.size(); a++)
	{
		const btVariationalAttachment& at = m_attachments[a];
		const btVec3d C = x[at.m_nodeA] - (at.m_nodeB >= 0 ? x[at.m_nodeB] : at.m_target) - at.m_c0;
		const double k = at.m_hard ? at.m_penalty : at.m_stiffness;
		E += at.m_lambda.dot(C) + 0.5 * k * C.length2();
	}
	for (int b = 0; b < m_bodies.size(); b++)
	{
		const btVariationalBody& body = m_bodies[b];
		if (body.m_colliders.size() == 0) continue;
		for (int i = body.m_firstNode; i < body.m_firstNode + body.m_numNodes; i++)
			if (!m_fixed[i]) addColliderNode(b, i, x[i], m_xPrev[i], E, 0, 0);
	}
	if (m_contact) E += m_contact->energy(x);
	return E;
}

double btVariationalSystem::elasticEnergy(const btVec3d* x) const
{
	double E = 0;
	for (int t = 0; t < m_tets.size(); t++)
	{
		const btVariationalTet& tet = m_tets[t];
		E += tet.m_restVolume * m_bodies[tet.m_body].m_material->energyDensity(deformationGradient(tet, x));
	}
	for (int s = 0; s < m_springs.size(); s++)
	{
		const btVariationalSpring& sp = m_springs[s];
		const double e = (x[sp.m_n[0]] - x[sp.m_n[1]]).length() - sp.m_restLength;
		E += 0.5 * sp.m_stiffness * e * e;
	}
	return E;
}

void btVariationalSystem::gradient(const btVec3d* x, btVec3d* g) const
{
	const int n = numNodes();
	const double invH2 = 1.0 / (m_dt * m_dt);
	for (int i = 0; i < n; i++) g[i] = m_fixed[i] ? btVec3d(0, 0, 0) : (x[i] - m_xTilde[i]) * (invH2 * m_mass[i]);
	for (int t = 0; t < m_tets.size(); t++)
	{
		const btVariationalTet& tet = m_tets[t];
		const btMat3d P = m_bodies[tet.m_body].m_material->firstPiola(deformationGradient(tet, x));
		for (int j = 0; j < 4; j++) g[tet.m_n[j]] += (P * tet.m_g[j]) * tet.m_restVolume;
	}
	for (int s = 0; s < m_springs.size(); s++)
	{
		const btVariationalSpring& sp = m_springs[s];
		const btVec3d d = x[sp.m_n[0]] - x[sp.m_n[1]];
		const double len = d.length();
		if (len < 1e-300) continue;
		const btVec3d f = d * (sp.m_stiffness * (len - sp.m_restLength) / len);
		g[sp.m_n[0]] += f;
		g[sp.m_n[1]] -= f;
	}
	for (int a = 0; a < m_attachments.size(); a++)
	{
		const btVariationalAttachment& at = m_attachments[a];
		const btVec3d C = x[at.m_nodeA] - (at.m_nodeB >= 0 ? x[at.m_nodeB] : at.m_target) - at.m_c0;
		const btVec3d f = at.m_lambda + C * (at.m_hard ? at.m_penalty : at.m_stiffness);
		g[at.m_nodeA] += f;
		if (at.m_nodeB >= 0) g[at.m_nodeB] -= f;
	}
	for (int b = 0; b < m_bodies.size(); b++)
	{
		const btVariationalBody& body = m_bodies[b];
		if (body.m_colliders.size() == 0) continue;
		double unused = 0;
		for (int i = body.m_firstNode; i < body.m_firstNode + body.m_numNodes; i++)
			if (!m_fixed[i]) addColliderNode(b, i, x[i], m_xPrev[i], unused, &g[i], 0);
	}
	if (m_contact) m_contact->addGradient(x, g);
	for (int i = 0; i < n; i++)
		if (m_fixed[i]) g[i].setZero();
}

void btVariationalSystem::hessian(const btVec3d* x, btVariationalHessian& H, btVariationalHessianFilter filter) const
{
	const int n = numNodes();
	const double invH2 = 1.0 / (m_dt * m_dt);
	H.reset(n);
	for (int i = 0; i < n; i++) H.m_nodeBlocks[i].addScaledIdentity(m_fixed[i] ? 1.0 : invH2 * m_mass[i]);
	H.m_stencils.reserve(m_tets.size() + m_springs.size());
	for (int t = 0; t < m_tets.size(); t++)
	{
		const btVariationalTet& tet = m_tets[t];
		btVariationalEigenSystem es;
		m_bodies[tet.m_body].m_material->eigenSystem(deformationGradient(tet, x), es);
		H.m_stencils.expand();
		btVariationalHessianStencil& st = H.m_stencils[H.m_stencils.size() - 1];
		st.m_H.setZero();
		for (int j = 0; j < 4; j++) st.m_n[j] = tet.m_n[j];
		for (int i = 0; i < 9; i++)
		{
			const double l = tet.m_restVolume * btVariationalFilterEigenvalue(es.m_lambda[i], filter);
			if (l == 0.0) continue;
			btVec12d w;
			for (int j = 0; j < 4; j++)
			{
				const btVec3d wj = es.m_Q[i] * tet.m_g[j];
				w[3 * j] = wj[0];
				w[3 * j + 1] = wj[1];
				w[3 * j + 2] = wj[2];
			}
			st.m_H.addScaledOuter(l, w, w);
		}
	}
	for (int s = 0; s < m_springs.size(); s++)
	{
		const btVariationalSpring& sp = m_springs[s];
		const btVec3d d = x[sp.m_n[0]] - x[sp.m_n[1]];
		const double len = d.length();
		if (len < 1e-300) continue;
		const btVec3d dir = d / len;
		// k [ n n^T + (1 - L/l)(I - n n^T) ]; the second coefficient is the only eigenvalue that can be negative
		const double c = btVariationalFilterEigenvalue(1.0 - sp.m_restLength / len, filter);
		btMat3d B;
		B.addScaledIdentity(sp.m_stiffness * c);
		B.addScaledOuter(sp.m_stiffness * (1.0 - c), dir, dir);
		H.m_stencils.expand();
		btVariationalHessianStencil& st = H.m_stencils[H.m_stencils.size() - 1];
		st.m_H.setZero();
		st.m_n[0] = sp.m_n[0];
		st.m_n[1] = sp.m_n[1];
		st.m_n[2] = st.m_n[3] = -1;
		st.m_H.addBlock3(0, 0, B);
		st.m_H.addBlock3(1, 1, B);
		st.m_H.addBlock3(0, 1, B * -1.0);
		st.m_H.addBlock3(1, 0, B * -1.0);
	}
	for (int a = 0; a < m_attachments.size(); a++)
	{
		const btVariationalAttachment& at = m_attachments[a];
		const double k = at.m_hard ? at.m_penalty : at.m_stiffness;
		if (at.m_nodeB < 0)
		{
			H.m_nodeBlocks[at.m_nodeA].addScaledIdentity(k);
			continue;
		}
		H.m_stencils.expand();
		btVariationalHessianStencil& st = H.m_stencils[H.m_stencils.size() - 1];
		st.m_H.setZero();
		st.m_n[0] = at.m_nodeA;
		st.m_n[1] = at.m_nodeB;
		st.m_n[2] = st.m_n[3] = -1;
		btMat3d B = btMat3d::identity() * k;
		st.m_H.addBlock3(0, 0, B);
		st.m_H.addBlock3(1, 1, B);
		st.m_H.addBlock3(0, 1, B * -1.0);
		st.m_H.addBlock3(1, 0, B * -1.0);
	}
	for (int b = 0; b < m_bodies.size(); b++)
	{
		const btVariationalBody& body = m_bodies[b];
		if (body.m_colliders.size() == 0) continue;
		double unused = 0;
		for (int i = body.m_firstNode; i < body.m_firstNode + body.m_numNodes; i++)
			if (!m_fixed[i]) addColliderNode(b, i, x[i], m_xPrev[i], unused, 0, &H.m_nodeBlocks[i]);
	}
	if (m_contact) m_contact->addHessian(x, H, filter);
}

void btVariationalSystem::vertexForceAndHessian(int node, const btVec3d* x, btVec3d& force, btMat3d& hess) const
{
	const double invH2 = 1.0 / (m_dt * m_dt);
	force = (x[node] - m_xTilde[node]) * (-invH2 * m_mass[node]);
	hess.setZero();
	hess.addScaledIdentity(invH2 * m_mass[node]);
	for (int p = m_nodeTetStart[node]; p < m_nodeTetStart[node + 1]; p++)
	{
		const btVariationalTet& tet = m_tets[m_nodeTets[p] >> 2];
		const int j = m_nodeTets[p] & 3;
		btVec3d Pg;
		btMat3d block;
		m_bodies[tet.m_body].m_material->vertexBlock(deformationGradient(tet, x), tet.m_g[j], Pg, block);
		force -= Pg * tet.m_restVolume;
		hess += block * tet.m_restVolume;
	}
	for (int p = m_nodeSpringStart[node]; p < m_nodeSpringStart[node + 1]; p++)
	{
		const btVariationalSpring& sp = m_springs[m_nodeSprings[p] >> 1];
		const double sign = (m_nodeSprings[p] & 1) ? -1.0 : 1.0;
		const btVec3d d = x[sp.m_n[0]] - x[sp.m_n[1]];
		const double len = d.length();
		if (len < 1e-300) continue;
		const btVec3d dir = d / len;
		force -= dir * (sign * sp.m_stiffness * (len - sp.m_restLength));
		const double c = 1.0 - sp.m_restLength / len;
		const double cc = c > 0.0 ? c : 0.0;
		hess.addScaledIdentity(sp.m_stiffness * cc);
		hess.addScaledOuter(sp.m_stiffness * (1.0 - cc), dir, dir);
	}
	for (int p = m_nodeAttachStart[node]; p < m_nodeAttachStart[node + 1]; p++)
	{
		const btVariationalAttachment& at = m_attachments[m_nodeAttach[p] >> 1];
		const double sign = (m_nodeAttach[p] & 1) ? -1.0 : 1.0;
		const btVec3d C = x[at.m_nodeA] - (at.m_nodeB >= 0 ? x[at.m_nodeB] : at.m_target) - at.m_c0;
		const double k = at.m_hard ? at.m_penalty : at.m_stiffness;
		force -= (at.m_lambda + C * k) * sign;
		hess.addScaledIdentity(k);
	}
	const int b = m_nodeBody[node];
	if (m_bodies[b].m_colliders.size() > 0)
	{
		double unused = 0;
		btVec3d grad(0, 0, 0);
		addColliderNode(b, node, x[node], m_xPrev[node], unused, &grad, &hess);
		force -= grad;
	}
	if (m_contact) m_contact->addVertexForceAndHessian(node, x, force, hess);
}

double btVariationalSystem::maxStepSize(const btVec3d* x, const btVec3d* p) const
{
	double alpha = 1.0;
	if (m_colliderModel == COLLIDER_BARRIER)
	{
		// The distance to a convex obstacle is convex along a ray, so Newton steps on d(t) = s * d(0)
		// approach the time of impact from below. Other shapes use the Lipschitz bound |d'| <= |p|.
		const double keep = 0.2;  // fraction of the initial gap that is preserved
		for (int b = 0; b < m_bodies.size(); b++)
		{
			const btVariationalBody& body = m_bodies[b];
			for (int s = 0; s < body.m_colliders.size(); s++)
			{
				const btVariationalCollider& col = m_colliders[body.m_colliders[s]];
				for (int i = body.m_firstNode; i < body.m_firstNode + body.m_numNodes; i++)
				{
					const double pLen = p[i].length();
					if (pLen == 0.0) continue;
					double d0;
					btVec3d nrm;
					col.signedDistance(x[i], d0, nrm);
					if (d0 <= 0.0 || d0 > pLen * alpha + m_colliderMargin) continue;
					const double target = keep * d0;
					double t = 0.0, d = d0;
					for (int it = 0; it < 32 && t < alpha; it++)
					{
						const double rate = col.isConvexExact() ? -nrm.dot(p[i]) : pLen;
						if (rate <= 0.0)
						{
							t = alpha;  // moving away, and convexity keeps it that way
							break;
						}
						const double dt = (d - target) / rate;
						if (dt <= 1e-12 * (1.0 + t)) break;
						t += dt;
						if (t >= alpha) break;
						col.signedDistance(x[i] + p[i] * t, d, nrm);
						if (d <= target * (1.0 + 1e-9)) break;
					}
					if (t < alpha) alpha = t;
				}
			}
		}
	}
	if (m_contact)
	{
		const double a = m_contact->maxStepSize(x, p);
		if (a < alpha) alpha = a;
	}
	return alpha;
}

double btVariationalSystem::gradientNorm(const btVec3d* g) const
{
	double m = 0;
	for (int i = 0; i < numNodes(); i++)
		if (!m_fixed[i])
		{
			const double a = g[i].maxAbs();
			if (a > m) m = a;
		}
	return m;
}

double btVariationalSystem::kineticEnergy() const
{
	double E = 0;
	for (int i = 0; i < numNodes(); i++) E += 0.5 * m_mass[i] * m_v[i].length2();
	return E;
}

double btVariationalSystem::totalVolume(const btVec3d* x) const
{
	double V = 0;
	for (int t = 0; t < m_tets.size(); t++) V += m_tets[t].m_restVolume * deformationGradient(m_tets[t], x).det();
	return V;
}

double btVariationalSystem::colliderPenetration(const btVec3d* x) const
{
	double deepest = 0;
	for (int b = 0; b < m_bodies.size(); b++)
	{
		const btVariationalBody& body = m_bodies[b];
		for (int s = 0; s < body.m_colliders.size(); s++)
		{
			const btVariationalCollider& col = m_colliders[body.m_colliders[s]];
			for (int i = body.m_firstNode; i < body.m_firstNode + body.m_numNodes; i++)
			{
				double d;
				btVec3d nrm;
				col.signedDistance(x[i], d, nrm);
				if (-d > deepest) deepest = -d;
			}
		}
	}
	return deepest;
}
