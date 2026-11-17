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

#include "btVariationalHelpers.h"
#include "BulletSoftBody/btSoftBodyHelpers.h"

btSoftBody* btVariationalHelpers::createTetGrid(btSoftBodyWorldInfo& worldInfo, const btVector3& center, const btVector3& halfExtents,
												int nx, int ny, int nz, btScalar density)
{
	const int px = nx + 1, py = ny + 1, pz = nz + 1;
	btAlignedObjectArray<btVector3> pos;
	pos.resize(px * py * pz);
	for (int k = 0; k < pz; k++)
		for (int j = 0; j < py; j++)
			for (int i = 0; i < px; i++)
				pos[i + px * (j + py * k)] = center + btVector3(halfExtents.x() * (btScalar(2 * i) / nx - 1),
																  halfExtents.y() * (btScalar(2 * j) / ny - 1),
																  halfExtents.z() * (btScalar(2 * k) / nz - 1));
	btSoftBody* psb = new btSoftBody(&worldInfo, pos.size(), &pos[0], 0);
	// the six tetrahedra of a cell all contain its main diagonal; one per ordering of the three axes
	static const int perm[6][3] = {{0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}};
	const int stride[3] = {1, px, px * py};
	for (int k = 0; k < nz; k++)
		for (int j = 0; j < ny; j++)
			for (int i = 0; i < nx; i++)
			{
				const int v0 = i + px * (j + py * k);
				for (int p = 0; p < 6; p++)
				{
					int n0 = v0;
					int n1 = n0 + stride[perm[p][0]];
					int n2 = n1 + stride[perm[p][1]];
					int n3 = n2 + stride[perm[p][2]];
					// keep a positive orientation so that generated boundary faces point outwards
					const btVector3 a = pos[n1] - pos[n0], b = pos[n2] - pos[n0], c = pos[n3] - pos[n0];
					if (a.cross(b).dot(c) < 0)
					{
						int tmp = n1;
						n1 = n2;
						n2 = tmp;
					}
					psb->appendTetra(n0, n1, n2, n3);
				}
			}
	psb->initializeDmInverse();
	btSoftBodyHelpers::generateBoundaryFaces(psb);
	setLumpedMass(psb, density);
	return psb;
}

namespace
{
// maps a point of the cube [-1, 1]^3 to the unit ball: the cube shell through the point goes to the
// sphere of the same "radius" (the largest coordinate), using the area-friendly cube to sphere map
btVector3 btCubeToBall(const btVector3& p)
{
	const btScalar m = btMax(btFabs(p.x()), btMax(btFabs(p.y()), btFabs(p.z())));
	if (m < btScalar(1e-12)) return btVector3(0, 0, 0);
	const btVector3 q = p / m;  // on the surface of the unit cube
	const btScalar x2 = q.x() * q.x(), y2 = q.y() * q.y(), z2 = q.z() * q.z();
	const btVector3 onSphere(q.x() * btSqrt(btScalar(1) - y2 / 2 - z2 / 2 + y2 * z2 / 3),
							 q.y() * btSqrt(btScalar(1) - z2 / 2 - x2 / 2 + z2 * x2 / 3),
							 q.z() * btSqrt(btScalar(1) - x2 / 2 - y2 / 2 + x2 * y2 / 3));
	return onSphere * m;
}

// six tetrahedra around the main diagonal of a hexahedral cell given by its corner indices,
// corner[a + 2 * b + 4 * c] with a, b, c in {0, 1} the offsets along the three grid directions
void btAppendCellTetras(btSoftBody* psb, const btAlignedObjectArray<btVector3>& pos, const int* corner)
{
	static const int perm[6][3] = {{0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}};
	static const int bit[3] = {1, 2, 4};
	for (int p = 0; p < 6; p++)
	{
		int code = 0;
		int n[4];
		n[0] = corner[0];
		for (int step = 0; step < 3; step++)
		{
			code += bit[perm[p][step]];
			n[step + 1] = corner[code];
		}
		// keep a positive orientation so that generated boundary faces point outwards
		const btVector3 a = pos[n[1]] - pos[n[0]], b = pos[n[2]] - pos[n[0]], c = pos[n[3]] - pos[n[0]];
		if (a.cross(b).dot(c) < 0)
		{
			const int tmp = n[1];
			n[1] = n[2];
			n[2] = tmp;
		}
		psb->appendTetra(n[0], n[1], n[2], n[3]);
	}
}
}  // namespace

btSoftBody* btVariationalHelpers::createTetEllipsoid(btSoftBodyWorldInfo& worldInfo, const btVector3& center, const btVector3& radii,
													 int nx, int ny, int nz, btScalar density, btScalar roundness)
{
	const int px = nx + 1, py = ny + 1, pz = nz + 1;
	btAlignedObjectArray<btVector3> pos;
	pos.resize(px * py * pz);
	for (int k = 0; k < pz; k++)
		for (int j = 0; j < py; j++)
			for (int i = 0; i < px; i++)
			{
				const btVector3 inCube(btScalar(2 * i) / nx - 1, btScalar(2 * j) / ny - 1, btScalar(2 * k) / nz - 1);
				const btVector3 u = inCube * (btScalar(1) - roundness) + btCubeToBall(inCube) * roundness;
				pos[i + px * (j + py * k)] = center + btVector3(radii.x() * u.x(), radii.y() * u.y(), radii.z() * u.z());
			}
	btSoftBody* psb = new btSoftBody(&worldInfo, pos.size(), &pos[0], 0);
	for (int k = 0; k < nz; k++)
		for (int j = 0; j < ny; j++)
			for (int i = 0; i < nx; i++)
			{
				int corner[8];
				for (int c = 0; c < 8; c++) corner[c] = (i + (c & 1)) + px * ((j + ((c >> 1) & 1)) + py * (k + ((c >> 2) & 1)));
				btAppendCellTetras(psb, pos, corner);
			}
	psb->initializeDmInverse();
	btSoftBodyHelpers::generateBoundaryFaces(psb);
	setLumpedMass(psb, density);
	return psb;
}

btSoftBody* btVariationalHelpers::createTetTorus(btSoftBodyWorldInfo& worldInfo, const btVector3& center, btScalar majorRadius, btScalar minorRadius,
												 int segments, int n, btScalar density)
{
	const int p1 = n + 1, perSection = p1 * p1;
	btAlignedObjectArray<btVector3> pos;
	pos.resize(segments * perSection);
	for (int s = 0; s < segments; s++)
	{
		const btScalar phi = SIMD_2_PI * btScalar(s) / btScalar(segments);
		const btVector3 radial(btCos(phi), 0, btSin(phi));
		for (int j = 0; j < p1; j++)
			for (int i = 0; i < p1; i++)
			{
				// square to disc
				const btScalar a = btScalar(2 * i) / n - 1, b = btScalar(2 * j) / n - 1;
				const btScalar da = a * btSqrt(btScalar(1) - b * b / 2), db = b * btSqrt(btScalar(1) - a * a / 2);
				pos[s * perSection + i + p1 * j] = center + radial * (majorRadius + minorRadius * da) + btVector3(0, minorRadius * db, 0);
			}
	}
	btSoftBody* psb = new btSoftBody(&worldInfo, pos.size(), &pos[0], 0);
	for (int s = 0; s < segments; s++)
		for (int j = 0; j < n; j++)
			for (int i = 0; i < n; i++)
			{
				int corner[8];
				for (int c = 0; c < 8; c++) corner[c] = ((s + (c & 1)) % segments) * perSection + (i + ((c >> 1) & 1)) + p1 * (j + ((c >> 2) & 1));
				btAppendCellTetras(psb, pos, corner);
			}
	psb->initializeDmInverse();
	btSoftBodyHelpers::generateBoundaryFaces(psb);
	setLumpedMass(psb, density);
	return psb;
}

void btVariationalHelpers::gridSurfaceQuads(int nx, int ny, int nz, btAlignedObjectArray<int>& quads)
{
	quads.resize(0);
	const int px = nx + 1, py = ny + 1;
	const int counts[3] = {nx, ny, nz};
	// faces with constant coordinate `axis`; (u, v, axis) is a right handed triple of directions
	for (int axis = 0; axis < 3; axis++)
	{
		const int u = (axis + 1) % 3, v = (axis + 2) % 3;
		for (int side = 0; side < 2; side++)
			for (int b = 0; b < counts[v]; b++)
				for (int a = 0; a < counts[u]; a++)
				{
					int ijk[4][3];
					static const int du[4] = {0, 1, 1, 0}, dv[4] = {0, 0, 1, 1};
					for (int c = 0; c < 4; c++)
					{
						// outward normal is +axis on the far side (counter-clockwise in u, v) and -axis on the near side
						const int cc = side == 1 ? c : 3 - c;
						ijk[c][axis] = side == 1 ? counts[axis] : 0;
						ijk[c][u] = a + du[cc];
						ijk[c][v] = b + dv[cc];
					}
					for (int c = 0; c < 4; c++) quads.push_back(ijk[c][0] + px * (ijk[c][1] + py * ijk[c][2]));
				}
	}
}

void btVariationalHelpers::torusSurfaceQuads(int segments, int n, btAlignedObjectArray<int>& quads)
{
	quads.resize(0);
	const int p1 = n + 1, perSection = p1 * p1;
	// boundary of the cross-section grid as a closed loop of 4 n nodes, counter-clockwise in (i, j)
	btAlignedObjectArray<int> loop;
	for (int i = 0; i < n; i++) loop.push_back(i);
	for (int j = 0; j < n; j++) loop.push_back(n + p1 * j);
	for (int i = n; i > 0; i--) loop.push_back(i + p1 * n);
	for (int j = n; j > 0; j--) loop.push_back(p1 * j);
	for (int s = 0; s < segments; s++)
	{
		const int s1 = (s + 1) % segments;
		for (int k = 0; k < loop.size(); k++)
		{
			const int k1 = (k + 1) % loop.size();
			// i points outwards from the axis, j along the axis, the sweep goes from x towards z:
			// this order makes the normal point out of the tube
			quads.push_back(s * perSection + loop[k]);
			quads.push_back(s * perSection + loop[k1]);
			quads.push_back(s1 * perSection + loop[k1]);
			quads.push_back(s1 * perSection + loop[k]);
		}
	}
}

void btVariationalHelpers::setLumpedMass(btSoftBody* psb, btScalar density)
{
	btAlignedObjectArray<btScalar> mass;
	mass.resize(psb->m_nodes.size(), btScalar(0));
	const btSoftBody::Node* base = &psb->m_nodes[0];
	for (int t = 0; t < psb->m_tetras.size(); t++)
	{
		const btSoftBody::Tetra& tet = psb->m_tetras[t];
		const btVector3 a = tet.m_n[1]->m_x - tet.m_n[0]->m_x, b = tet.m_n[2]->m_x - tet.m_n[0]->m_x, c = tet.m_n[3]->m_x - tet.m_n[0]->m_x;
		const btScalar quarter = density * btFabs(a.cross(b).dot(c)) / btScalar(24);
		for (int j = 0; j < 4; j++) mass[int(tet.m_n[j] - base)] += quarter;
	}
	for (int i = 0; i < psb->m_nodes.size(); i++)
		if (psb->m_nodes[i].m_im > 0 && mass[i] > 0) psb->setMass(i, mass[i]);
}
