"""Tetrahedralizes a closed surface for the deformable examples.

    python make_tet_mesh.py surface.obj volume.vtk [edge_length_fac=0.16] [epsilon=0.006]

Reads the triangles of a Wavefront OBJ file, fills them with tetrahedra using fTetWild (through the
pytetwild package: pip install pytetwild pyvista) and writes a legacy ASCII VTK file of the kind that
btSoftBodyHelpers::CreateFromVtkFile reads. The output is deliberately coarse: its surface follows
the input only within epsilon times the diagonal of the bounding box. Simulate the coarse volume and
render the original surface, moved along with the tetrahedra it lies in (VariationalMeshExport.h
does that).

data/bunny.vtk, data/teddy.vtk and data/duck.vtk were made with this script and the default settings
from the OBJ files of the same names.
"""
import os
import sys
import numpy as np
import pytetwild


def read_obj_triangles(path):
    vertices, triangles = [], []
    for line in open(path):
        parts = line.split()
        if not parts:
            continue
        if parts[0] == "v":
            vertices.append([float(x) for x in parts[1:4]])
        elif parts[0] == "f":
            index = [int(p.split("/")[0]) for p in parts[1:]]
            index = [i - 1 if i > 0 else len(vertices) + i for i in index]
            for k in range(1, len(index) - 1):
                triangles.append([index[0], index[k], index[k + 1]])
    return np.array(vertices, dtype=np.float64), np.array(triangles, dtype=np.int32)


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        return 1
    source, target = sys.argv[1], sys.argv[2]
    edge_length_fac = float(sys.argv[3]) if len(sys.argv) > 3 else 0.16
    epsilon = float(sys.argv[4]) if len(sys.argv) > 4 else 0.006
    vertices, triangles = read_obj_triangles(source)
    points, tets = pytetwild.tetrahedralize(vertices, triangles, edge_length_fac=edge_length_fac, epsilon=epsilon,
                                           optimize=True, coarsen=True, quiet=True)
    if os.path.exists("__tracked_surface.stl"):  # fTetWild leaves this behind in the working directory
        os.remove("__tracked_surface.stl")
    tets = np.array(tets, dtype=np.int64)
    # positive orientation
    a, b, c = (points[tets[:, k]] - points[tets[:, 0]] for k in (1, 2, 3))
    volume = np.einsum("ij,ij->i", np.cross(a, b), c) / 6.0
    flip = volume < 0
    tets[flip, 1], tets[flip, 2] = tets[flip, 2].copy(), tets[flip, 1].copy()
    volume = np.abs(volume)
    edges = [a, b, c, b - a, c - a, c - b]
    mean_square = sum((e * e).sum(1) for e in edges) / 6.0
    quality = 6.0 * np.sqrt(2.0) * volume / mean_square ** 1.5  # 1 for a regular tetrahedron
    with open(target, "w") as f:
        f.write("# vtk DataFile Version 2.0\n%s, tetrahedralized with fTetWild\nASCII\nDATASET UNSTRUCTURED_GRID\n" % source.split("/")[-1])
        f.write("POINTS %d double\n" % len(points))
        for p in points:
            f.write("%.9g %.9g %.9g\n" % (p[0], p[1], p[2]))
        f.write("\nCELLS %d %d\n" % (len(tets), 5 * len(tets)))
        for t in tets:
            f.write("4 %d %d %d %d\n" % (t[0], t[1], t[2], t[3]))
        f.write("\nCELL_TYPES %d\n" % len(tets))
        f.write("10\n" * len(tets))
    print("%s: %d nodes, %d tetrahedra, volume %.6g, element quality mean %.2f, worst %.2f" %
          (target, len(points), len(tets), volume.sum(), quality.mean(), quality.min()))
    return 0


if __name__ == "__main__":
    sys.exit(main())
