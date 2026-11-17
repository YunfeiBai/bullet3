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

#ifndef VARIATIONAL_MESH_EXPORT_H
#define VARIATIONAL_MESH_EXPORT_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "BulletSoftBody/btSoftBody.h"
#include "BulletSoftBody/BulletVariationalDeformable/btVariationalMath.h"
#include "LinearMath/btAlignedObjectArray.h"

///Writes the surfaces of soft bodies frame by frame for offline rendering.
///
///  <directory>/scene.json   the bodies (surface polygons, material tag, texture coordinates, the frame
///                           from which a body exists) and static props such as posts and the ground
///  <directory>/frames.bin   for every frame the positions of all surface vertices of all bodies in
///                           the order of scene.json, as little endian 32 bit floats (x, y, z)
///
///The surface of a body is either made of its own nodes (addBody, addBodyFromFaces) or a separate,
///usually much finer triangle mesh that is carried along by the tetrahedra it lies in
///(addEmbeddedSurface): a coarse volume is simulated and the detailed surface is rendered.
///
///render_showcase.py reads this format.
class VariationalMeshExporter
{
	struct Body
	{
		const btSoftBody* m_psb;
		char m_name[48];
		char m_material[48];
		int m_faceSize;                           // 3 or 4
		btAlignedObjectArray<int> m_faces;        // indices into the exported vertices of the body
		btAlignedObjectArray<int> m_nodes;        // exported vertex -> node of the soft body
		btAlignedObjectArray<int> m_embedNodes;   // embedded surface: four nodes per exported vertex ...
		btAlignedObjectArray<float> m_embedWeights;  // ... and their weights (m_nodes is then just a counter)
		btAlignedObjectArray<float> m_uv;         // optional, two per exported vertex
		btAlignedObjectArray<float> m_cornerUv;   // optional, two per face corner (textures with seams)
		char m_texture[64];                       // optional image file name
		int m_firstFrame;                         // -1: not in the scene yet
		float m_smoothing;                        // hint for the renderer: subdivision levels
	};
	struct Prop
	{
		char m_type[16];
		char m_material[48];
		float m_a[3], m_b[3], m_radius;
	};
	btAlignedObjectArray<Body> m_bodies;
	btAlignedObjectArray<Prop> m_props;
	char m_directory[512];
	FILE* m_frames;
	int m_numFrames;
	bool m_sceneWritten;
	double m_framesPerSecond;

	static void copyName(char* destination, size_t size, const char* source)
	{
		strncpy(destination, source, size - 1);
		destination[size - 1] = 0;
	}

	void writeScene()
	{
		char path[640];
		snprintf(path, sizeof(path), "%s/scene.json", m_directory);
		FILE* f = fopen(path, "w");
		if (!f) return;
		int total = 0;
		for (int b = 0; b < m_bodies.size(); b++) total += m_bodies[b].m_nodes.size();
		fprintf(f, "{\n \"frames_per_second\": %g,\n \"num_frames\": %d,\n \"vertices_per_frame\": %d,\n \"bodies\": [\n", m_framesPerSecond, m_numFrames, total);
		for (int b = 0; b < m_bodies.size(); b++)
		{
			const Body& body = m_bodies[b];
			fprintf(f, "  {\"name\": \"%s\", \"material\": \"%s\", \"num_vertices\": %d, \"face_size\": %d, \"first_frame\": %d, \"smoothing\": %g,\n   \"faces\": [",
					body.m_name, body.m_material, body.m_nodes.size(), body.m_faceSize, body.m_firstFrame, body.m_smoothing);
			for (int i = 0; i < body.m_faces.size(); i++) fprintf(f, i ? ",%d" : "%d", body.m_faces[i]);
			fprintf(f, "]");
			if (body.m_uv.size() > 0)
			{
				fprintf(f, ",\n   \"uv\": [");
				for (int i = 0; i < body.m_uv.size(); i++) fprintf(f, i ? ",%.5f" : "%.5f", body.m_uv[i]);
				fprintf(f, "]");
			}
			if (body.m_cornerUv.size() > 0)
			{
				fprintf(f, ",\n   \"corner_uv\": [");
				for (int i = 0; i < body.m_cornerUv.size(); i++) fprintf(f, i ? ",%.5f" : "%.5f", body.m_cornerUv[i]);
				fprintf(f, "]");
			}
			if (body.m_texture[0]) fprintf(f, ",\n   \"texture\": \"%s\"", body.m_texture);
			fprintf(f, "}%s\n", b + 1 < m_bodies.size() ? "," : "");
		}
		fprintf(f, " ],\n \"props\": [\n");
		for (int p = 0; p < m_props.size(); p++)
		{
			const Prop& prop = m_props[p];
			fprintf(f, "  {\"type\": \"%s\", \"material\": \"%s\", \"a\": [%g, %g, %g], \"b\": [%g, %g, %g], \"radius\": %g}%s\n", prop.m_type, prop.m_material,
					prop.m_a[0], prop.m_a[1], prop.m_a[2], prop.m_b[0], prop.m_b[1], prop.m_b[2], prop.m_radius, p + 1 < m_props.size() ? "," : "");
		}
		fprintf(f, " ]\n}\n");
		fclose(f);
	}

public:
	VariationalMeshExporter() : m_frames(0), m_numFrames(0), m_sceneWritten(false), m_framesPerSecond(60)
	{
		m_directory[0] = 0;
	}
	~VariationalMeshExporter() { close(); }

	bool isOpen() const { return m_frames != 0; }
	int numFrames() const { return m_numFrames; }

	///the directory must exist
	bool open(const char* directory, double framesPerSecond)
	{
		close();
		copyName(m_directory, sizeof(m_directory), directory);
		m_framesPerSecond = framesPerSecond;
		char path[640];
		snprintf(path, sizeof(path), "%s/frames.bin", m_directory);
		m_frames = fopen(path, "wb");
		m_numFrames = 0;
		return m_frames != 0;
	}

	///Registers a body with its surface polygons (node indices of the soft body, faceSize per face).
	///Only the nodes that the polygons use are exported. uv, if given, has two entries per node of the
	///soft body. Returns a handle for setVisible().
	int addBody(const btSoftBody* psb, const char* name, const char* material, const int* faces, int numFaces, int faceSize, const btScalar* uv = 0, float smoothing = 1)
	{
		Body body;
		body.m_psb = psb;
		copyName(body.m_name, sizeof(body.m_name), name);
		copyName(body.m_material, sizeof(body.m_material), material);
		body.m_faceSize = faceSize;
		body.m_firstFrame = -1;
		body.m_smoothing = smoothing;
		body.m_texture[0] = 0;
		btAlignedObjectArray<int> exported;
		exported.resize(psb->m_nodes.size(), -1);
		for (int i = 0; i < numFaces * faceSize; i++)
		{
			const int node = faces[i];
			if (exported[node] < 0)
			{
				exported[node] = body.m_nodes.size();
				body.m_nodes.push_back(node);
				if (uv)
				{
					body.m_uv.push_back(float(uv[2 * node]));
					body.m_uv.push_back(float(uv[2 * node + 1]));
				}
			}
			body.m_faces.push_back(exported[node]);
		}
		m_bodies.push_back(body);
		return m_bodies.size() - 1;
	}

	///the triangles Bullet generated for a body (btSoftBodyHelpers::generateBoundaryFaces)
	int addBodyFromFaces(const btSoftBody* psb, const char* name, const char* material, float smoothing = 1)
	{
		btAlignedObjectArray<int> triangles;
		const btSoftBody::Node* base = &psb->m_nodes[0];
		for (int f = 0; f < psb->m_faces.size(); f++)
			for (int j = 0; j < 3; j++) triangles.push_back(int(psb->m_faces[f].m_n[j] - base));
		return addBody(psb, name, material, triangles.size() ? &triangles[0] : 0, triangles.size() / 3, 3, 0, smoothing);
	}

	///Uses the triangles of a Wavefront OBJ file as the surface of the body. Every vertex of the file
	///is tied to the tetrahedron of the body that contains it (or comes closest to containing it) and
	///follows that tetrahedron from then on. Call this while the body is still where its mesh file
	///put it, i.e. in the coordinates of the OBJ file. texture names an image for the texture
	///coordinates of the file, if it has any. Returns -1 if the file cannot be read.
	int addEmbeddedSurface(const btSoftBody* psb, const char* name, const char* material, const char* objPath, const char* texture = 0, float smoothing = 1)
	{
		FILE* f = fopen(objPath, "r");
		if (!f || psb->m_tetras.size() == 0)
		{
			if (f) fclose(f);
			return -1;
		}
		btAlignedObjectArray<btVec3d> positions;
		btAlignedObjectArray<float> uvs;
		Body body;
		body.m_psb = psb;
		copyName(body.m_name, sizeof(body.m_name), name);
		copyName(body.m_material, sizeof(body.m_material), material);
		copyName(body.m_texture, sizeof(body.m_texture), texture ? texture : "");
		body.m_faceSize = 3;
		body.m_firstFrame = -1;
		body.m_smoothing = smoothing;
		bool allCornersHaveUv = true;
		char line[1024];
		while (fgets(line, sizeof(line), f))
		{
			if (line[0] == 'v' && line[1] == ' ')
			{
				double x, y, z;
				if (sscanf(line + 2, "%lf %lf %lf", &x, &y, &z) == 3) positions.push_back(btVec3d(x, y, z));
			}
			else if (line[0] == 'v' && line[1] == 't')
			{
				double u = 0, v = 0;
				sscanf(line + 2, "%lf %lf", &u, &v);
				uvs.push_back(float(u));
				uvs.push_back(float(v));
			}
			else if (line[0] == 'f' && line[1] == ' ')
			{
				// corners "v", "v/vt", "v//vn" or "v/vt/vn"; polygons become triangle fans
				int vertex[64], coordinate[64], count = 0;
				for (char* token = strtok(line + 2, " \t\r\n"); token && count < 64; token = strtok(0, " \t\r\n"))
				{
					int v = atoi(token), t = 0;
					const char* slash = strchr(token, '/');
					if (slash && slash[1] != '/' && slash[1] != 0) t = atoi(slash + 1);
					vertex[count] = v > 0 ? v - 1 : positions.size() + v;
					coordinate[count] = t > 0 ? t - 1 : (t < 0 ? uvs.size() / 2 + t : -1);
					count++;
				}
				for (int k = 1; k + 1 < count; k++)
				{
					const int corner[3] = {0, k, k + 1};
					for (int c = 0; c < 3; c++)
					{
						body.m_faces.push_back(vertex[corner[c]]);
						const int t = coordinate[corner[c]];
						if (t < 0 || 2 * t + 1 >= uvs.size())
							allCornersHaveUv = false;
						else if (allCornersHaveUv)
						{
							body.m_cornerUv.push_back(uvs[2 * t]);
							body.m_cornerUv.push_back(uvs[2 * t + 1]);
						}
					}
				}
			}
		}
		fclose(f);
		if (!allCornersHaveUv) body.m_cornerUv.resize(0);
		if (body.m_cornerUv.size() == 0) body.m_texture[0] = 0;

		// barycentric coordinates in every tetrahedron; keep the one in which the smallest is largest
		const btSoftBody::Node* base = &psb->m_nodes[0];
		const int numTets = psb->m_tetras.size();
		btAlignedObjectArray<btMat3d> inverse;
		btAlignedObjectArray<btVec3d> origin;
		inverse.resize(numTets);
		origin.resize(numTets);
		for (int t = 0; t < numTets; t++)
		{
			const btSoftBody::Tetra& tet = psb->m_tetras[t];
			origin[t] = btVec3d(tet.m_n[0]->m_x);
			const btMat3d edges = btMat3d::fromColumns(btVec3d(tet.m_n[1]->m_x) - origin[t], btVec3d(tet.m_n[2]->m_x) - origin[t], btVec3d(tet.m_n[3]->m_x) - origin[t]);
			inverse[t] = fabs(edges.det()) > 1e-300 ? edges.inverse() : btMat3d();
		}
		for (int v = 0; v < positions.size(); v++)
		{
			int best = 0;
			double bestInside = -1e300;
			btVec3d bestWeights(0, 0, 0);
			for (int t = 0; t < numTets; t++)
			{
				const btVec3d w = inverse[t] * (positions[v] - origin[t]);
				const double w0 = 1.0 - w[0] - w[1] - w[2];
				double inside = w0 < w[0] ? w0 : w[0];
				if (w[1] < inside) inside = w[1];
				if (w[2] < inside) inside = w[2];
				if (inside > bestInside)
				{
					bestInside = inside;
					best = t;
					bestWeights = w;
				}
			}
			const btSoftBody::Tetra& tet = psb->m_tetras[best];
			for (int j = 0; j < 4; j++) body.m_embedNodes.push_back(int(tet.m_n[j] - base));
			body.m_embedWeights.push_back(float(1.0 - bestWeights[0] - bestWeights[1] - bestWeights[2]));
			for (int j = 0; j < 3; j++) body.m_embedWeights.push_back(float(bestWeights[j]));
			body.m_nodes.push_back(v);
		}
		m_bodies.push_back(body);
		return m_bodies.size() - 1;
	}

	///a body is written (and rendered) from the frame at which it is first made visible
	void setVisible(int body)
	{
		if (m_bodies[body].m_firstFrame < 0) m_bodies[body].m_firstFrame = m_numFrames;
	}

	///static scenery: "capsule" (segment a-b with a radius), "ground" (a: a point on it, b: its normal),
	///or anything else the render script knows how to draw from two points and a radius
	void addProp(const char* type, const char* material, const btVector3& a, const btVector3& b, btScalar radius)
	{
		Prop prop;
		copyName(prop.m_type, sizeof(prop.m_type), type);
		copyName(prop.m_material, sizeof(prop.m_material), material);
		for (int k = 0; k < 3; k++)
		{
			prop.m_a[k] = float(a[k]);
			prop.m_b[k] = float(b[k]);
		}
		prop.m_radius = float(radius);
		m_props.push_back(prop);
	}

	void writeFrame()
	{
		if (!m_frames) return;
		btAlignedObjectArray<float> buffer;
		for (int b = 0; b < m_bodies.size(); b++)
		{
			const Body& body = m_bodies[b];
			const bool embedded = body.m_embedNodes.size() > 0;
			for (int i = 0; i < body.m_nodes.size(); i++)
			{
				btVector3 x(0, 0, 0);
				if (embedded)
					for (int j = 0; j < 4; j++) x += body.m_psb->m_nodes[body.m_embedNodes[4 * i + j]].m_x * btScalar(body.m_embedWeights[4 * i + j]);
				else
					x = body.m_psb->m_nodes[body.m_nodes[i]].m_x;
				buffer.push_back(float(x.x()));
				buffer.push_back(float(x.y()));
				buffer.push_back(float(x.z()));
			}
		}
		if (buffer.size() > 0) fwrite(&buffer[0], sizeof(float), buffer.size(), m_frames);
		m_numFrames++;
		// keep the description current, so that a run that is still going (or was stopped) can be rendered
		if ((m_numFrames % 30) == 1)
		{
			fflush(m_frames);
			writeScene();
		}
	}

	void close()
	{
		if (!m_frames) return;
		fclose(m_frames);
		m_frames = 0;
		writeScene();
	}
};

#endif  // VARIATIONAL_MESH_EXPORT_H
