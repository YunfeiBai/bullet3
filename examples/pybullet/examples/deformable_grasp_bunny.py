"""Pick and place of a deformable bunny with a parallel-jaw gripper.

A more demanding version of deformable_grasp.py: the gripper (x and y slides, z lift, a
wrist that turns around the vertical axis and two prismatic fingers, all driven by force
limited position control) picks up a volumetric Neo-Hookean bunny, lifts it, carries it
over to a tray while turning it by 90 degrees, shakes it, holds it still and puts it down
in the tray. The irregular shape gives contacts on curved, sloped surfaces, and the bunny
has to be held against gravity, inertia and the wrist rotation through friction alone.

Usage:
  python deformable_grasp_bunny.py                    # interactive GUI
  python deformable_grasp_bunny.py --direct --metrics bunny.json --video bunny.mp4

In --direct mode the script prints (and optionally saves) these metrics:
  intact                the bunny mesh did not blow up
  held                  the bunny was still in the gripper at the end of the hold phase
  placed                the bunny was carried over and ended up in the tray
  slip_mm               how far the bunny moved relative to the gripper between the start of
                        the transport and the end of the hold (0 for a perfect grasp)
  jitter_rms_mm_per_s   RMS speed of the bunny nodes relative to the gripper while the
                        gripper holds still (0 for a perfectly static grasp)
  penetration_max_mm    deepest penetration of a bunny node into a finger while it is held
  grip_force_*          mean / standard deviation of the finger motor forces while holding

bunny.vtk in pybullet_data is bunny.obj tetrahedralized with gmsh (1357 nodes); its elastic
forces are integrated explicitly, so stiffer materials need a smaller --timestep, e.g.
--mu 8000 --lam 12000 --timestep 0.002.
"""
import argparse
import json
import math
import os
import tempfile
import time

import numpy as np
import pybullet as p
import pybullet_data

from deformable_grasp import VideoWriter, box_signed_distance, smoothstep

BUNNY_LENGTH = 0.13  # [m] bunny.vtk has length 1, faces +x and its bottom is at z = 0
GRASP_X = -0.01  # [m] where the fingers grip the body (the centroid is at x = 0)
FINGER_CENTER_Z = 0.045  # [m] height of the finger centres when grasping
FINGER_THICKNESS = 0.012
FINGER_OPEN_GAP = 0.13  # distance between the finger pads when the gripper is open [m]
PALM_THICKNESS = 0.02
TRAY_POSITION = (0.32, 0.12, 0.0)
START_HEIGHT = 0.30  # palm height at the start [m]
CARRY_HEIGHT = 0.16  # how high the bunny is lifted [m]
CAMERA = (0.45, 150, -25)  # distance, yaw, pitch; the camera follows the gripper

SLIDE_X, SLIDE_Y, LIFT_Z, WRIST, LEFT, RIGHT = range(6)  # joint (and link) indices


def gripper_urdf(finger_width, finger_height):
  """A fixed-base gripper: x slide -> y slide -> z lift -> wrist (yaw) -> palm -> fingers."""
  open_y = 0.5 * (FINGER_OPEN_GAP + FINGER_THICKNESS)
  palm = (finger_width + 0.01, FINGER_OPEN_GAP + 2 * FINGER_THICKNESS + 0.02, PALM_THICKNESS)
  inertia = '<inertia ixx="{0}" ixy="0" ixz="0" iyy="{0}" iyz="0" izz="{0}"/>'
  box = '<geometry><box size="%g %g %g"/></geometry>'

  def carriage(name):
    return '<link name="%s"><inertial><mass value="0.5"/>%s</inertial></link>' % (
        name, inertia.format(1e-3))

  def prismatic(name, parent, child, axis):
    return ('<joint name="%s" type="prismatic"><parent link="%s"/><child link="%s"/>'
            '<axis xyz="%s"/><limit lower="-1" upper="1" effort="1000" velocity="2"/></joint>' %
            (name, parent, child, axis))

  finger_z = -0.5 * (PALM_THICKNESS + finger_height)
  finger_box = box % (finger_width, FINGER_THICKNESS, finger_height)
  fingers = ""
  for name, y, axis in (("left", -open_y, "0 1 0"), ("right", open_y, "0 -1 0")):
    fingers += """
  <link name="{name}_finger">
    <inertial><mass value="0.1"/>{inertia}</inertial>
    <visual>{box}<material name="dark"><color rgba="0.25 0.27 0.3 1"/></material></visual>
    <collision>{box}</collision>
  </link>
  <joint name="{name}_finger_joint" type="prismatic">
    <parent link="palm"/><child link="{name}_finger"/>
    <origin xyz="0 {y} {z}"/><axis xyz="{axis}"/>
    <limit lower="0" upper="{travel}" effort="100" velocity="1"/>
  </joint>""".format(name=name, inertia=inertia.format(5e-5), box=finger_box, y=y, z=finger_z,
                     axis=axis, travel=0.5 * FINGER_OPEN_GAP)
  return """<?xml version="1.0"?>
<robot name="parallel_jaw_gripper_with_wrist">
  <link name="base"/>
  {slide_x}{slide_y}{lift}
  <link name="palm">
    <inertial><mass value="0.3"/>{palm_inertia}</inertial>
    <visual>{palm}<material name="grey"><color rgba="0.55 0.57 0.6 1"/></material></visual>
    <collision>{palm}</collision>
  </link>
  {joint_x}{joint_y}{joint_z}
  <joint name="wrist" type="revolute">
    <parent link="lift"/><child link="palm"/><axis xyz="0 0 1"/>
    <limit lower="-3.2" upper="3.2" effort="100" velocity="10"/>
  </joint>{fingers}
</robot>
""".format(slide_x=carriage("slide_x"), slide_y=carriage("slide_y"), lift=carriage("lift"),
           palm_inertia=inertia.format(1e-3), palm=box % palm,
           joint_x=prismatic("slide_x", "base", "slide_x", "1 0 0"),
           joint_y=prismatic("slide_y", "slide_x", "slide_y", "0 1 0"),
           joint_z=prismatic("lift_z", "slide_y", "lift", "0 0 1"), fingers=fingers)


class Schedule(object):
  """Scripted pick and place."""
  phases = [  # name, start [s], end [s]
      ("settle", 0.0, 0.5),
      ("descend", 0.5, 1.3),
      ("close", 1.3, 2.1),
      ("lift", 2.1, 3.0),
      ("transport", 3.0, 4.6),  # move over the tray and turn the wrist by 90 degrees
      ("shake", 4.6, 5.4),  # two vertical oscillations of +-2 cm at 2.5 Hz
      ("hold", 5.4, 6.4),
      ("lower", 6.4, 7.2),
      ("release", 7.2, 7.7),
      ("retract", 7.7, 8.4),
  ]
  duration = phases[-1][2]

  def __init__(self, grasp_height, finger_close):
    self.grasp_height = grasp_height
    self.finger_close = finger_close

  def phase(self, t):
    for name, t0, t1 in self.phases:
      if t0 <= t < t1:
        return name
    return "done"

  def targets(self, t):
    """(x, y, z, wrist yaw, finger) joint position targets at time t."""
    transport = smoothstep(3.0, 4.6, t)
    x = GRASP_X + (TRAY_POSITION[0] - GRASP_X) * transport
    y = TRAY_POSITION[1] * transport
    z = START_HEIGHT + (self.grasp_height - START_HEIGHT) * smoothstep(0.5, 1.3, t)
    z += CARRY_HEIGHT * smoothstep(2.1, 3.0, t)
    if 4.6 <= t < 5.4:
      z += 0.02 * math.sin(2 * math.pi * 2.5 * (t - 4.6))
    z -= (CARRY_HEIGHT - 0.02) * smoothstep(6.4, 7.2, t)
    z += 0.15 * smoothstep(7.7, 8.4, t)
    yaw = 0.5 * math.pi * transport
    finger = self.finger_close * (smoothstep(1.3, 2.0, t) - smoothstep(7.2, 7.6, t))
    return x, y, z, yaw, finger


class Demo(object):

  def __init__(self, args):
    self.args = args
    self.finger_half = (0.5 * args.finger_width, 0.5 * FINGER_THICKNESS, 0.5 * args.finger_height)
    p.resetSimulation(p.RESET_USE_DEFORMABLE_WORLD)
    p.setGravity(0, 0, -9.81)
    p.setTimeStep(args.timestep)
    p.setPhysicsEngineParameter(numSolverIterations=args.iterations)
    p.setAdditionalSearchPath(pybullet_data.getDataPath())

    self.plane = p.loadURDF("plane.urdf")
    p.changeDynamics(self.plane, -1, lateralFriction=0.6)
    self.tray = p.loadURDF("tray/traybox.urdf", TRAY_POSITION, useFixedBase=True, globalScaling=0.5)

    urdf_path = os.path.join(tempfile.mkdtemp(prefix="deformable_grasp_bunny_"), "gripper.urdf")
    with open(urdf_path, "w") as f:
      f.write(gripper_urdf(args.finger_width, args.finger_height))
    self.gripper = p.loadURDF(urdf_path, [0, 0, 0], useFixedBase=True)
    for link in (LEFT, RIGHT):
      p.changeDynamics(self.gripper, link, lateralFriction=args.finger_friction)
    p.resetJointState(self.gripper, SLIDE_X, GRASP_X)
    p.resetJointState(self.gripper, LIFT_Z, START_HEIGHT)

    self.bunny = p.loadSoftBody("bunny.vtk",
                                basePosition=[0, 0, 0.003],
                                scale=BUNNY_LENGTH,
                                mass=args.bunny_mass,
                                useNeoHookean=1,
                                NeoHookeanMu=args.mu,
                                NeoHookeanLambda=args.lam,
                                NeoHookeanDamping=args.damping,
                                collisionMargin=0.002,
                                frictionCoeff=args.bunny_friction,
                                useSelfCollision=0,
                                useFaceContact=1)
    p.changeVisualShape(self.bunny, -1, rgbaColor=[0.93, 0.74, 0.55, 1])

    grasp_height = FINGER_CENTER_Z + 0.5 * (PALM_THICKNESS + args.finger_height)
    finger_close = 0.5 * FINGER_OPEN_GAP - args.grasp_half_width
    self.schedule = Schedule(grasp_height, finger_close)

  def command(self, t):
    x, y, z, yaw, finger = self.schedule.targets(t)
    for joint, target, force in ((SLIDE_X, x, 500), (SLIDE_Y, y, 500), (LIFT_Z, z, 500),
                                 (WRIST, yaw, 50)):
      p.setJointMotorControl2(self.gripper, joint, p.POSITION_CONTROL, targetPosition=target,
                              force=force, positionGain=0.2, velocityGain=1.0)
    for joint in (LEFT, RIGHT):
      p.setJointMotorControl2(self.gripper, joint, p.POSITION_CONTROL, targetPosition=finger,
                              force=self.args.grip_force, positionGain=0.2, velocityGain=1.0,
                              maxVelocity=0.2)

  def nodes(self):
    return np.asarray(p.getMeshData(self.bunny, -1, flags=p.MESH_DATA_SIMULATION_MESH)[1])

  def camera_target(self):
    palm = p.getLinkState(self.gripper, WRIST, computeForwardKinematics=1)[0]
    return [palm[0] + 0.04, palm[1] + 0.02, max(0.08, palm[2] - 0.11)]


def render(width, height, target):
  distance, yaw, pitch = CAMERA
  view = p.computeViewMatrixFromYawPitchRoll(cameraTargetPosition=target, distance=distance,
                                             yaw=yaw, pitch=pitch, roll=0, upAxisIndex=2)
  proj = p.computeProjectionMatrixFOV(fov=40, aspect=float(width) / height, nearVal=0.01,
                                      farVal=5.0)
  img = p.getCameraImage(width, height, view, proj, renderer=p.ER_TINY_RENDERER,
                         lightDirection=[0.6, -0.5, 1.2], shadow=1)
  return np.asarray(img[2], dtype=np.uint8).reshape(height, width, 4)


def simulate(args, demo, video=None):
  """Runs the scripted pick and place once and returns the per-step measurements."""
  dt = args.timestep
  frame_every = max(1, int(round(1.0 / (args.fps * dt))))
  prev_local = None
  samples = []
  for step in range(int(round(Schedule.duration / dt))):
    demo.command(step * dt)
    p.stepSimulation()
    t = (step + 1) * dt

    nodes = demo.nodes()
    palm_pos, palm_orn = p.getLinkState(demo.gripper, WRIST, computeForwardKinematics=1)[:2]
    rot = np.asarray(p.getMatrixFromQuaternion(palm_orn)).reshape(3, 3)
    local = (nodes - np.asarray(palm_pos)).dot(rot)  # bunny nodes in the palm frame
    fingers = [p.getLinkState(demo.gripper, j, computeForwardKinematics=1) for j in (LEFT, RIGHT)]
    penetration = max(
        float(-box_signed_distance(nodes, f[0], f[1], demo.finger_half).min()) for f in fingers)
    # RMS speed of the nodes relative to the gripper (in the palm frame, so that the wrist
    # rotation does not count as motion)
    rel_speed = float("nan")
    if prev_local is not None:
      rel_v = (local - prev_local) / dt
      rel_speed = float(np.sqrt((rel_v * rel_v).sum(axis=1).mean()))
    prev_local = local
    joints = [p.getJointState(demo.gripper, j) for j in (LEFT, RIGHT)]
    samples.append({
        "t": t,
        "phase": demo.schedule.phase(t - 1e-9),
        "com": nodes.mean(axis=0).tolist(),
        "com_in_gripper": local.mean(axis=0).tolist(),
        "extent": float(np.linalg.norm(nodes.max(axis=0) - nodes.min(axis=0))),
        "penetration": max(penetration, 0.0),
        "rel_speed": rel_speed,
        "grip_force": 0.5 * (abs(joints[0][3]) + abs(joints[1][3])),
    })

    if video is not None and step % frame_every == 0:
      video.add(render(args.width, args.height, demo.camera_target()))
    if not args.direct:
      time.sleep(dt)
  return samples


def compute_metrics(samples):
  by_phase = {}
  for s in samples:
    by_phase.setdefault(s["phase"], []).append(s)
  carried = by_phase["transport"] + by_phase["shake"] + by_phase["hold"]
  static = by_phase["hold"][len(by_phase["hold"]) // 2:]  # the gripper does not move
  start = np.asarray(by_phase["transport"][0]["com_in_gripper"])
  end = np.asarray(by_phase["hold"][-1]["com_in_gripper"])
  final_com = np.asarray(samples[-1]["com"])
  # the mesh did not blow up: it never grew to more than 1.5 times its size
  intact = max(s["extent"] for s in samples) < 1.5 * samples[0]["extent"]
  held = intact and np.linalg.norm(end) < 0.1
  in_tray = (np.abs(final_com[:2] - np.asarray(TRAY_POSITION[:2])).max() < 0.12 and
             final_com[2] < 0.1)
  return {
      "intact": bool(intact),
      "held": bool(held),
      "placed": bool(held and in_tray),
      "slip_mm": 1000.0 * float(np.linalg.norm(end - start)),
      "jitter_rms_mm_per_s": 1000.0 * float(np.sqrt(np.mean([s["rel_speed"]**2 for s in static]))),
      "jitter_max_mm_per_s": 1000.0 * max(s["rel_speed"] for s in static),
      "penetration_max_mm": 1000.0 * max(s["penetration"] for s in carried),
      "grip_force_mean_N": float(np.mean([s["grip_force"] for s in static])),
      "grip_force_std_N": float(np.std([s["grip_force"] for s in static])),
  }


def report(args, samples):
  metrics = compute_metrics(samples)
  print(json.dumps(metrics, indent=2))
  if args.metrics:
    with open(args.metrics, "w") as f:
      json.dump({"metrics": metrics, "config": vars(args),
                 "series": {k: [s[k] for s in samples] for k in samples[0]}}, f)
  return metrics


def run(args):
  video = VideoWriter(args.video, args.width, args.height, args.fps) if args.video else None
  if args.direct:
    p.connect(p.DIRECT)
    samples = simulate(args, Demo(args), video)
    if video is not None:
      video.close()
    p.disconnect()
    return report(args, samples)

  # GUI: replay the sequence until the window is closed
  p.connect(p.GUI)
  p.configureDebugVisualizer(p.COV_ENABLE_GUI, 0)
  metrics = None
  try:
    while p.isConnected():
      demo = Demo(args)
      p.resetDebugVisualizerCamera(CAMERA[0], CAMERA[1], CAMERA[2], demo.camera_target())
      samples = simulate(args, demo, video)
      if video is not None:
        video.close()
        video = None
      metrics = report(args, samples)
  except p.error:
    pass  # the window was closed
  return metrics


def parse_args(argv=None):
  parser = argparse.ArgumentParser(description=__doc__,
                                   formatter_class=argparse.RawDescriptionHelpFormatter)
  parser.add_argument("--direct", action="store_true", help="run headless (no GUI)")
  parser.add_argument("--timestep", type=float, default=1.0 / 240.0)
  parser.add_argument("--iterations", type=int, default=50, help="solver iterations")
  parser.add_argument("--grip-force", type=float, default=10.0, help="finger motor force [N]")
  parser.add_argument("--grasp-half-width", type=float, default=0.017,
                      help="finger target: half the distance between the pads when closed [m]")
  parser.add_argument("--finger-width", type=float, default=0.06, help="[m]")
  parser.add_argument("--finger-height", type=float, default=0.08, help="[m]")
  parser.add_argument("--bunny-mass", type=float, default=0.15, help="[kg]")
  parser.add_argument("--mu", type=float, default=5000.0, help="Neo-Hookean mu [Pa]")
  parser.add_argument("--lam", type=float, default=7500.0, help="Neo-Hookean lambda [Pa]")
  parser.add_argument("--damping", type=float, default=0.02, help="Neo-Hookean damping")
  parser.add_argument("--bunny-friction", type=float, default=1.0)
  parser.add_argument("--finger-friction", type=float, default=1.0)
  parser.add_argument("--metrics", help="write metrics and time series (json) to this file")
  parser.add_argument("--video", help="write an mp4 video to this file (needs ffmpeg)")
  parser.add_argument("--width", type=int, default=640)
  parser.add_argument("--height", type=int, default=480)
  parser.add_argument("--fps", type=int, default=30)
  return parser.parse_args(argv)


if __name__ == "__main__":
  run(parse_args())
