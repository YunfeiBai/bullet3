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
