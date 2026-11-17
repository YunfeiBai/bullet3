"""Renders the frames written by App_VariationalShowcase with Blender's Cycles.

    python render_showcase.py <export directory> <output directory> [options]
    blender -b -P render_showcase.py -- <export directory> <output directory> [options]

The first form needs Blender as a Python module (pip install bpy). Options:

    --data DIR        Bullet's data folder, for the textures of the duck and the banana
    --start S --end S time range in seconds of simulated time (default: 1.2 to the end)
    --fps N           frames per second of the video (default 30)
    --width W --height H   resolution (default 1280 x 720)
    --samples N       samples per pixel (default 48; the image is denoised)
    --frames A B      render only output frames A..B-1 (to split the work or to resume)
    --still T         render one frame at simulated time T and stop
    --camera NAME     "story" (default, the animated camera), or a fixed "wide" / "close" view

One PNG per frame is written to the output directory (frames that exist are skipped); turn them
into a video with e.g.

    ffmpeg -framerate 30 -i <output directory>/frame_%04d.png -c:v libx264 -crf 17 -pix_fmt yuv420p showcase.mp4

The simulation uses Bullet's convention of y pointing up; the scene is turned so that z is up.
"""
import json
import math
import os
import sys
import time

import bpy
import bmesh
import numpy as np
from mathutils import Vector


def parse_arguments():
    argv = sys.argv[1:]
    if "--" in argv:
        argv = argv[argv.index("--") + 1:]
    options = {"data": None, "start": 1.2, "end": None, "fps": 30, "width": 1280, "height": 720, "samples": 48,
               "frames": None, "still": None, "camera": "story"}
    positional = []
    i = 0
    while i < len(argv):
        a = argv[i]
        if a == "--frames":
            options["frames"] = (int(argv[i + 1]), int(argv[i + 2]))
            i += 3
        elif a.startswith("--") and a[2:] in options:
            key = a[2:]
            value = argv[i + 1]
            options[key] = value if key in ("data", "camera") else (float(value) if key in ("start", "end", "still", "fps") else int(value))
            i += 2
        else:
            positional.append(a)
            i += 1
    if len(positional) < 2:
        print(__doc__)
        sys.exit(1)
    options["export"], options["output"] = positional[0], positional[1]
    return options


def to_blender(points):
    """Bullet (x, y up, z) -> Blender (x, -z, z up)"""
    out = np.empty_like(points)
    out[..., 0] = points[..., 0]
    out[..., 1] = -points[..., 2]
    out[..., 2] = points[..., 1]
    return out


# ---------------------------------------------------------------- materials

def principled(name, color, roughness=0.5, metallic=0.0, sheen=0.0, coat=0.0, subsurface=0.0, specular=0.5):
    material = bpy.data.materials.new(name)
    material.use_nodes = True
    bsdf = material.node_tree.nodes["Principled BSDF"]

    def set_input(names, value):
        for n in names:
            if n in bsdf.inputs:
                bsdf.inputs[n].default_value = value
                return

    set_input(["Base Color"], (color[0], color[1], color[2], 1.0))
    set_input(["Roughness"], roughness)
    set_input(["Metallic"], metallic)
    set_input(["Sheen Weight", "Sheen"], sheen)
    set_input(["Coat Weight", "Clearcoat"], coat)
    set_input(["Specular IOR Level", "Specular"], specular)
    if subsurface > 0:
        set_input(["Subsurface Weight", "Subsurface"], subsurface)
        set_input(["Subsurface Radius"], (color[0] * 0.02 + 0.004, color[1] * 0.02 + 0.004, color[2] * 0.02 + 0.004))
        set_input(["Subsurface Scale"], 1.0)
    return material


def add_texture(material, path):
    tree = material.node_tree
    bsdf = tree.nodes["Principled BSDF"]
    image = tree.nodes.new("ShaderNodeTexImage")
    image.image = bpy.data.images.load(path)
    tree.links.new(image.outputs["Color"], bsdf.inputs["Base Color"])


def add_fabric_pattern(material, color_a, color_b, squares):
    """gingham: two sets of stripes that darken where they cross"""
    tree = material.node_tree
    bsdf = tree.nodes["Principled BSDF"]
    coordinates = tree.nodes.new("ShaderNodeTexCoord")
    checker = tree.nodes.new("ShaderNodeTexChecker")
    checker.inputs["Scale"].default_value = squares
    checker.inputs["Color1"].default_value = (color_a[0], color_a[1], color_a[2], 1)
    checker.inputs["Color2"].default_value = (color_b[0], color_b[1], color_b[2], 1)
    tree.links.new(coordinates.outputs["UV"], checker.inputs["Vector"])
    tree.links.new(checker.outputs["Color"], bsdf.inputs["Base Color"])
    # a fine weave in the normal
    weave = tree.nodes.new("ShaderNodeTexChecker")
    weave.inputs["Scale"].default_value = squares * 24
    bump = tree.nodes.new("ShaderNodeBump")
    bump.inputs["Strength"].default_value = 0.15
    tree.links.new(coordinates.outputs["UV"], weave.inputs["Vector"])
    tree.links.new(weave.outputs["Fac"], bump.inputs["Height"])
    tree.links.new(bump.outputs["Normal"], bsdf.inputs["Normal"])


def add_fuzz(material, scale=160.0, strength=0.5):
    """a little noise in the normal: plush instead of plastic"""
    tree = material.node_tree
    bsdf = tree.nodes["Principled BSDF"]
    noise = tree.nodes.new("ShaderNodeTexNoise")
    noise.inputs["Scale"].default_value = scale
    bump = tree.nodes.new("ShaderNodeBump")
    bump.inputs["Strength"].default_value = strength
    bump.inputs["Distance"].default_value = 0.002
    tree.links.new(noise.outputs["Fac"], bump.inputs["Height"])
    tree.links.new(bump.outputs["Normal"], bsdf.inputs["Normal"])


def make_materials(data_directory):
    m = {}
    m["floor"] = principled("floor", (0.40, 0.39, 0.38), roughness=0.5, specular=0.35)
    m["wood"] = principled("wood", (0.45, 0.27, 0.13), roughness=0.5)
    m["rope"] = principled("rope", (0.75, 0.66, 0.48), roughness=0.9)
    m["steel"] = principled("steel", (0.55, 0.57, 0.60), roughness=0.28, metallic=1.0)
    m["chain"] = principled("chain", (0.30, 0.31, 0.33), roughness=0.35, metallic=1.0)
    m["frame"] = principled("frame", (0.16, 0.17, 0.19), roughness=0.4, metallic=0.9)
    m["cloth"] = principled("cloth", (0.8, 0.8, 0.8), roughness=0.85, sheen=0.4, specular=0.2)
    add_fabric_pattern(m["cloth"], (0.80, 0.14, 0.12), (0.93, 0.90, 0.84), 14.0)
    m["toy_teddy"] = principled("toy_teddy", (0.42, 0.24, 0.11), roughness=0.95, sheen=0.9, specular=0.15)
    add_fuzz(m["toy_teddy"])
    m["toy_bunny"] = principled("toy_bunny", (0.93, 0.88, 0.86), roughness=0.9, sheen=0.8, specular=0.2)
    add_fuzz(m["toy_bunny"], strength=0.35)
    m["toy_duck"] = principled("toy_duck", (0.95, 0.75, 0.10), roughness=0.3, coat=0.4)
    m["toy_banana"] = principled("toy_banana", (0.92, 0.80, 0.25), roughness=0.5)
    m["toy_bread"] = principled("toy_bread", (0.78, 0.55, 0.30), roughness=0.8, subsurface=0.3)
    add_fuzz(m["toy_bread"], scale=60.0, strength=0.25)
    m["toy_red"] = principled("toy_red", (0.80, 0.08, 0.07), roughness=0.3, coat=0.5)
    m["toy_teal"] = principled("toy_teal", (0.03, 0.55, 0.52), roughness=0.3, coat=0.5, subsurface=0.2)
    m["toy_purple"] = principled("toy_purple", (0.42, 0.16, 0.70), roughness=0.3, coat=0.5, subsurface=0.2)
    m["toy_yellow"] = principled("toy_yellow", (0.95, 0.68, 0.05), roughness=0.35, coat=0.3, subsurface=0.2)
    m["toy_blue"] = principled("toy_blue", (0.08, 0.30, 0.85), roughness=0.35, coat=0.3, subsurface=0.2)
    m["default"] = principled("default", (0.6, 0.6, 0.6))
    return m


# ---------------------------------------------------------------- geometry

def new_object(name, mesh, material):
    obj = bpy.data.objects.new(name, mesh)
    bpy.context.scene.collection.objects.link(obj)
    mesh.materials.append(material)
    return obj


def build_body(body, positions, materials, data_directory):
    """positions: rest positions (n, 3) in Blender coordinates"""
    size = body["face_size"]
    faces = np.array(body["faces"], dtype=np.int64).reshape(-1, size)
    mesh = bpy.data.meshes.new(body["name"])
    mesh.from_pydata([tuple(p) for p in positions], [], [tuple(int(i) for i in f) for f in faces])
    mesh.update()
    # consistent outward normals (the bodies are closed surfaces; a cloth just needs consistency)
    bm = bmesh.new()
    bm.from_mesh(mesh)
    bmesh.ops.recalc_face_normals(bm, faces=bm.faces)
    bm.to_mesh(mesh)
    bm.free()
    if "uv" in body or "corner_uv" in body:
        layer = mesh.uv_layers.new(name="uv")
        if "corner_uv" in body:
            # (face orientation may have been flipped above: go through the loops by vertex)
            corner_uv = np.array(body["corner_uv"], dtype=np.float32).reshape(-1, size, 2)
            for polygon in mesh.polygons:
                original = faces[polygon.index]
                for loop_index in polygon.loop_indices:
                    vertex = mesh.loops[loop_index].vertex_index
                    k = int(np.nonzero(original == vertex)[0][0])
                    layer.data[loop_index].uv = (float(corner_uv[polygon.index, k, 0]), float(corner_uv[polygon.index, k, 1]))
        else:
            uv = np.array(body["uv"], dtype=np.float32).reshape(-1, 2)
            for loop in mesh.loops:
                layer.data[loop.index].uv = (float(uv[loop.vertex_index, 0]), float(uv[loop.vertex_index, 1]))
    mesh.polygons.foreach_set("use_smooth", [True] * len(mesh.polygons))
    material = materials.get(body["material"], materials["default"])
    texture = body.get("texture")
    if texture and data_directory and os.path.exists(os.path.join(data_directory, texture)):
        material = material.copy()
        add_texture(material, os.path.join(data_directory, texture))
    obj = new_object(body["name"], mesh, material)
    levels = int(body.get("smoothing", 1))
    if body["material"] == "cloth":
        solidify = obj.modifiers.new("thickness", "SOLIDIFY")
        solidify.thickness = 0.003
        solidify.offset = 0.0
    if levels > 0:
        subdivision = obj.modifiers.new("subdivision", "SUBSURF")
        subdivision.levels = levels
        subdivision.render_levels = levels
    return obj


def cylinder_between(name, a, b, radius, material, caps=True):
    a, b = Vector(a), Vector(b)
    axis = b - a
    mesh = bpy.data.meshes.new(name)
    bm = bmesh.new()
    bmesh.ops.create_cone(bm, cap_ends=True, segments=32, radius1=radius, radius2=radius, depth=axis.length)
    if caps:
        for end in (-0.5, 0.5):
            created = bmesh.ops.create_uvsphere(bm, u_segments=24, v_segments=12, radius=radius)
            bmesh.ops.translate(bm, verts=created["verts"], vec=(0, 0, end * axis.length))
    bm.to_mesh(mesh)
    bm.free()
    mesh.polygons.foreach_set("use_smooth", [True] * len(mesh.polygons))
    obj = new_object(name, mesh, material)
    obj.location = (a + b) * 0.5
    obj.rotation_mode = "QUATERNION"
    obj.rotation_quaternion = Vector((0, 0, 1)).rotation_difference(axis.normalized())
    return obj


def build_props(props, materials):
    for index, prop in enumerate(props):
        a = to_blender(np.array(prop["a"], dtype=np.float64))
        b = to_blender(np.array(prop["b"], dtype=np.float64))
        material = materials.get(prop["material"], materials["default"])
        kind = prop["type"]
        if kind == "ground":
            # a photo studio cove: the floor bends up into a wall all around, so there is no horizon
            flat, bend, wall, segments = 14.0, 5.0, 10.0, 96
            profile = [(0.0, 0.0), (flat, 0.0)]
            profile += [(flat + bend * math.sin(u), bend * (1.0 - math.cos(u))) for u in np.linspace(0.0, math.pi / 2, 13)[1:]]
            profile.append((flat + bend, bend + wall))
            vertices, faces = [(0.0, 0.0, 0.0)], []
            for radius, height in profile[1:]:
                for k in range(segments):
                    angle = 2.0 * math.pi * k / segments
                    vertices.append((radius * math.cos(angle), radius * math.sin(angle), height))
            for k in range(segments):
                faces.append((0, 1 + k, 1 + (k + 1) % segments))
            for ring in range(len(profile) - 2):
                for k in range(segments):
                    a0, a1 = 1 + ring * segments + k, 1 + ring * segments + (k + 1) % segments
                    faces.append((a0, a0 + segments, a1 + segments, a1))
            mesh = bpy.data.meshes.new("ground")
            mesh.from_pydata(vertices, [], faces)
            mesh.update()
            mesh.polygons.foreach_set("use_smooth", [True] * len(mesh.polygons))
            new_object("ground", mesh, material).location = (a[0], a[1], a[2])
        elif kind in ("capsule", "rope"):
            cylinder_between("%s%d" % (kind, index), a, b, prop["radius"], material)
        elif kind == "pivot":
            # the chain hangs from an eye bolt under a beam that stays above the picture
            cylinder_between("eye", a, (a[0], a[1], a[2] + 0.08), 0.02, materials["steel"])
            cylinder_between("beam", (a[0], a[1] - 6.0, a[2] + 0.14), (a[0], a[1] + 6.0, a[2] + 0.14), 0.06, materials["frame"])


# ---------------------------------------------------------------- light and camera

def build_lighting(scene):
    """studio lighting: a large warm key light, a weak cool fill from the other side, a rim light"""
    world = bpy.data.worlds.new("world")
    scene.world = world
    world.use_nodes = True
    background = world.node_tree.nodes["Background"]
    background.inputs["Color"].default_value = (0.80, 0.86, 1.0, 1)
    background.inputs["Strength"].default_value = 0.09
    for name, location, energy, size, color in (("key", (4.5, -5.0, 6.5), 950.0, 3.5, (1.0, 0.93, 0.84)),
                                                ("fill", (-5.5, -4.0, 3.0), 90.0, 6.0, (0.78, 0.86, 1.0)),
                                                ("rim", (-2.5, 6.0, 5.0), 430.0, 4.0, (1.0, 0.97, 0.92))):
        light = bpy.data.lights.new(name, "AREA")
        light.energy = energy
        light.size = size
        light.color = color
        obj = bpy.data.objects.new(name, light)
        scene.collection.objects.link(obj)
        obj.location = location
        obj.rotation_mode = "QUATERNION"
        obj.rotation_quaternion = (Vector((0, 0, 0.9)) - Vector(location)).to_track_quat("-Z", "Y")


# camera path: (simulated time, position, look-at point, focal length), Blender coordinates
CAMERA_KEYS = {
    "story": [
        (1.2, (1.7, -2.3, 2.9), (0.05, 0.0, 0.95), 42.0),   # looking down into the hammock as the toys drop in
        (4.0, (1.5, -2.6, 2.8), (0.05, 0.0, 0.95), 42.0),
        (5.0, (0.5, -4.4, 2.2), (-0.5, 0.0, 1.3), 34.0),    # back and to the left: the ball is let go
        (5.7, (0.9, -4.0, 1.75), (0.35, 0.0, 1.0), 34.0),   # the impact
        (6.4, (1.2, -4.6, 1.7), (0.9, 0.0, 0.8), 33.0),     # after the flying toys
        (7.3, (1.1, -4.5, 1.6), (0.8, 0.0, 0.6), 35.0),
        (8.9, (0.0, -4.7, 1.6), (0.65, 0.0, 0.65), 36.0),   # slowly around what is left
    ],
    "wide": [(0.0, (1.0, -7.5, 2.6), (-0.2, 0.0, 1.5), 30.0)],
    "close": [(0.0, (2.6, -3.6, 2.6), (0.0, 0.0, 0.95), 40.0)],
}


def smooth(u):
    return u * u * (3.0 - 2.0 * u)


def camera_pose(keys, t):
    if t <= keys[0][0]:
        return keys[0][1:]
    if t >= keys[-1][0]:
        return keys[-1][1:]
    for k in range(len(keys) - 1):
        t0, t1 = keys[k][0], keys[k + 1][0]
        if t0 <= t <= t1:
            u = smooth((t - t0) / (t1 - t0))
            blend = lambda a, b: tuple(x + (y - x) * u for x, y in zip(a, b))
            return blend(keys[k][1], keys[k + 1][1]), blend(keys[k][2], keys[k + 1][2]), keys[k][3] + (keys[k + 1][3] - keys[k][3]) * u
    return keys[-1][1:]


def main():
    options = parse_arguments()
    with open(os.path.join(options["export"], "scene.json")) as f:
        description = json.load(f)
    per_frame = description["vertices_per_frame"]
    raw = np.memmap(os.path.join(options["export"], "frames.bin"), dtype=np.float32, mode="r")
    num_frames = raw.size // (3 * per_frame)
    frames = raw[:num_frames * 3 * per_frame].reshape(num_frames, per_frame, 3)
    simulation_fps = description["frames_per_second"]
    print("%d simulated frames at %g per second, %d vertices each" % (num_frames, simulation_fps, per_frame))

    bpy.ops.wm.read_factory_settings(use_empty=True)
    scene = bpy.context.scene
    scene.render.engine = "CYCLES"
    scene.cycles.samples = options["samples"]
    scene.cycles.use_adaptive_sampling = True
    scene.cycles.adaptive_threshold = 0.02
    scene.cycles.use_denoising = True
    try:
        scene.cycles.denoiser = "OPENIMAGEDENOISE"
    except Exception:
        pass
    scene.cycles.max_bounces = 6
    scene.cycles.caustics_reflective = False
    scene.cycles.caustics_refractive = False
    scene.render.resolution_x = options["width"]
    scene.render.resolution_y = options["height"]
    scene.render.resolution_percentage = 100
    scene.render.image_settings.file_format = "PNG"
    scene.render.use_persistent_data = True
    try:
        scene.view_settings.view_transform = "Standard"  # keeps the toy colours saturated
    except Exception:
        pass

    materials = make_materials(options["data"])
    build_props(description["props"], materials)
    build_lighting(scene)
    # faint tiles on the flat part of the floor, to read the motion against
    floor_tree = materials["floor"].node_tree
    coordinates = floor_tree.nodes.new("ShaderNodeTexCoord")
    flatten = floor_tree.nodes.new("ShaderNodeMapping")
    flatten.inputs["Scale"].default_value = (1.0, 1.0, 0.0)
    tiles = floor_tree.nodes.new("ShaderNodeTexChecker")
    tiles.inputs["Scale"].default_value = 2.0
    tiles.inputs["Color1"].default_value = (0.41, 0.40, 0.39, 1)
    tiles.inputs["Color2"].default_value = (0.36, 0.355, 0.345, 1)
    floor_tree.links.new(coordinates.outputs["Object"], flatten.inputs["Vector"])
    floor_tree.links.new(flatten.outputs["Vector"], tiles.inputs["Vector"])
    # ... and none where the cove bends up into the wall
    height = floor_tree.nodes.new("ShaderNodeSeparateXYZ")
    floor_tree.links.new(coordinates.outputs["Object"], height.inputs["Vector"])
    on_wall = floor_tree.nodes.new("ShaderNodeMath")
    on_wall.operation = "GREATER_THAN"
    on_wall.inputs[1].default_value = 0.001
    floor_tree.links.new(height.outputs["Z"], on_wall.inputs[0])
    blend = floor_tree.nodes.new("ShaderNodeMixRGB")
    blend.inputs["Color2"].default_value = (0.38, 0.375, 0.365, 1)
    floor_tree.links.new(on_wall.outputs["Value"], blend.inputs["Fac"])
    floor_tree.links.new(tiles.outputs["Color"], blend.inputs["Color1"])
    floor_tree.links.new(blend.outputs["Color"], floor_tree.nodes["Principled BSDF"].inputs["Base Color"])

    objects = []
    offset = 0
    for body in description["bodies"]:
        count = body["num_vertices"]
        rest = to_blender(np.array(frames[0, offset:offset + count], dtype=np.float64))
        objects.append((build_body(body, rest, materials, options["data"]), offset, count, body["first_frame"]))
        offset += count

    camera_data = bpy.data.cameras.new("camera")
    camera_data.dof.use_dof = True
    camera_data.dof.aperture_fstop = 5.6
    camera = bpy.data.objects.new("camera", camera_data)
    scene.collection.objects.link(camera)
    scene.camera = camera
    camera.rotation_mode = "QUATERNION"
    keys = CAMERA_KEYS[options["camera"]]

    os.makedirs(options["output"], exist_ok=True)
    if options["still"] is not None:
        times = [options["still"]]
    else:
        end = options["end"] if options["end"] is not None else num_frames / simulation_fps
        times = [options["start"] + k / options["fps"] for k in range(int((end - options["start"]) * options["fps"]))]
    first, last = options["frames"] if options["frames"] else (0, len(times))
    started = time.time()
    rendered = 0
    for index in range(first, min(last, len(times))):
        path = os.path.join(options["output"], "frame_%04d.png" % index)
        if os.path.exists(path):
            continue
        t = times[index]
        frame = min(num_frames - 1, max(0, int(round(t * simulation_fps)) - 1))
        for obj, start, count, first_frame in objects:
            visible = 0 <= first_frame <= frame
            obj.hide_render = not visible
            if visible:
                obj.data.vertices.foreach_set("co", to_blender(np.array(frames[frame, start:start + count], dtype=np.float32)).ravel())
                obj.data.update()
        position, target, lens = camera_pose(keys, t)
        camera.location = position
        camera.rotation_quaternion = (Vector(target) - Vector(position)).to_track_quat("-Z", "Y")
        camera_data.lens = lens
        camera_data.dof.focus_distance = (Vector(target) - Vector(position)).length
        scene.render.filepath = path
        bpy.ops.render.render(write_still=True)
        rendered += 1
        elapsed = time.time() - started
        print("frame %d (t = %.3f s) done, %.1f s per frame, %d to go" % (index, t, elapsed / rendered, min(last, len(times)) - index - 1))
        sys.stdout.flush()


if __name__ == "__main__":
    main()
