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

#ifndef BT_VARIATIONAL_HELPERS_H
#define BT_VARIATIONAL_HELPERS_H

#include "BulletSoftBody/btSoftBody.h"

///Convenience functions for setting up soft bodies for the variational solvers
struct btVariationalHelpers
{
	///Axis aligned block of nx * ny * nz cells with six tetrahedra per cell (Kuhn subdivision, which
	///is conforming across cells). Boundary faces are generated for rendering and collision, and the
	///nodes receive lumped masses for the given density.
	static btSoftBody* createTetGrid(btSoftBodyWorldInfo& worldInfo, const btVector3& center, const btVector3& halfExtents,
									 int nx, int ny, int nz, btScalar density);

	///The same grid with its nodes moved so that the block becomes an ellipsoid with the given radii
	///(roundness 1), a box (roundness 0, identical to createTetGrid) or a rounded box in between.
	///Concentric cube shells of the grid map to concentric ellipsoid shells, which keeps the elements
	///well shaped. Node (i, j, k) has the index i + (nx + 1) * (j + (ny + 1) * k), as in createTetGrid.
	static btSoftBody* createTetEllipsoid(btSoftBodyWorldInfo& worldInfo, const btVector3& center, const btVector3& radii,
										  int nx, int ny, int nz, btScalar density, btScalar roundness = btScalar(1));

	///Solid torus around the y axis through center. The cross-section is a grid of n * n cells mapped
	///to a disc of the minor radius and swept around the axis in the given number of segments.
	///Node (i, j) of the cross-section of segment s has the index s * (n + 1)^2 + i + (n + 1) * j.
	static btSoftBody* createTetTorus(btSoftBodyWorldInfo& worldInfo, const btVector3& center, btScalar majorRadius, btScalar minorRadius,
									  int segments, int n, btScalar density);

	///Quadrilaterals of the surface of a body made by createTetGrid / createTetEllipsoid (four node
	///indices per quad, counter-clockwise seen from outside). The triangles Bullet generates for
	///collision split each of these quads in two; for rendering the quads are the better control
	///mesh, e.g. for subdivision surfaces.
	static void gridSurfaceQuads(int nx, int ny, int nz, btAlignedObjectArray<int>& quads);
	///... and of a body made by createTetTorus
	static void torusSurfaceQuads(int segments, int n, btAlignedObjectArray<int>& quads);

	///Lumped nodal masses: every tetrahedron gives a quarter of its mass to each of its nodes.
	///Nodes that are currently fixed (zero inverse mass) stay fixed.
	static void setLumpedMass(btSoftBody* psb, btScalar density);
};

#endif  // BT_VARIATIONAL_HELPERS_H
