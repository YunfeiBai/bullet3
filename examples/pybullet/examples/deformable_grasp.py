"""Parallel-jaw gripper grasping, lifting, holding and transporting a deformable sphere.

The sphere is a volumetric (tetrahedral, Neo-Hookean FEM) deformable body. The gripper
is a fixed-base btMultiBody with four prismatic joints (x slide, z lift and two fingers)
driven by force-limited, position-controlled joint motors, as is common for robot
grippers. The script runs a scripted sequence: approach, close, lift, hold, transport
sideways, hold, lower and release.

Usage:
  python deformable_grasp.py                    # interactive GUI
  python deformable_grasp.py --direct --metrics grasp.json --video grasp.mp4

In --direct mode the script prints (and optionally saves) grasp-stability metrics:
  lifted                  the sphere was lifted off the ground and held
  slip_*_mm               how far the sphere slides down relative to the fingers
                          (0 for a perfect grasp)
  jitter_rms_mm_per_s     RMS speed of the sphere nodes relative to the gripper while
                          the gripper does not move (0 for a perfectly static grasp)
  penetration_max_mm      deepest penetration of a sphere node into a finger
  grip_force_*            mean / standard deviation of the finger motor forces

The elastic forces of the deformable body are integrated explicitly, so stiffer materials
(--mu, --lam) need a smaller --timestep, e.g. --mu 8000 --lam 12000 --timestep 0.001.
"""
import argparse
import json
import os
import subprocess
import tempfile
import time

import numpy as np
import pybullet as p
import pybullet_data

SPHERE_RADIUS = 0.05  # [m] ball.vtk is a unit sphere that gets scaled to this radius
FINGER_THICKNESS = 0.012  # [m]
FINGER_OPEN_GAP = 0.13  # distance between the finger pads when the gripper is open [m]
PALM_THICKNESS = 0.02  # [m]
LIFT_START = 0.25  # height of the palm at the start [m]
LIFT_HEIGHT = 0.20  # how high the sphere is lifted [m]
TRANSPORT_X = 0.20  # sideways transport distance [m]

SLIDE_X, LIFT_Z, LEFT, RIGHT = 0, 1, 2, 3  # joint (and link) indices of the gripper
CAMERA = (0.4, 85, -10)  # distance, yaw, pitch


def gripper_urdf(finger_width, finger_height):
  """A fixed-base gripper: x slide -> z lift -> palm -> two prismatic fingers."""
  fx, fy, fz = finger_width, FINGER_THICKNESS, finger_height
  open_y = 0.5 * FINGER_OPEN_GAP + 0.5 * FINGER_THICKNESS
  palm = (finger_width + 0.01, FINGER_OPEN_GAP + 2 * FINGER_THICKNESS + 0.02, PALM_THICKNESS)
  inertia = '<inertia ixx="{i}" ixy="0" ixz="0" iyy="{i}" iyz="0" izz="{i}"/>'
  box = '<geometry><box size="%g %g %g"/></geometry>'
  return """<?xml version="1.0"?>
<robot name="parallel_jaw_gripper">
  <link name="base"/>
  <link name="slide">
    <inertial><mass value="0.5"/>{i3}</inertial>
  </link>
  <link name="palm">
    <inertial><mass value="0.5"/>{i3}</inertial>
    <visual>{palm}<material name="grey"><color rgba="0.55 0.57 0.6 1"/></material></visual>
    <collision>{palm}</collision>
  </link>
  <link name="left_finger">
    <inertial><mass value="0.1"/>{i5}</inertial>
    <visual>{finger}<material name="dark"><color rgba="0.25 0.27 0.3 1"/></material></visual>
    <collision>{finger}</collision>
  </link>
  <link name="right_finger">
    <inertial><mass value="0.1"/>{i5}</inertial>
    <visual>{finger}<material name="dark"><color rgba="0.25 0.27 0.3 1"/></material></visual>
    <collision>{finger}</collision>
  </link>
  <joint name="slide_x" type="prismatic">
    <parent link="base"/><child link="slide"/><axis xyz="1 0 0"/>
    <limit lower="-1" upper="1" effort="1000" velocity="2"/>
  </joint>
  <joint name="lift_z" type="prismatic">
    <parent link="slide"/><child link="palm"/><axis xyz="0 0 1"/>
    <limit lower="-1" upper="1" effort="1000" velocity="2"/>
  </joint>
  <joint name="left_finger_joint" type="prismatic">
    <parent link="palm"/><child link="left_finger"/>
    <origin xyz="0 {ly} {fzc}"/><axis xyz="0 1 0"/>
    <limit lower="0" upper="{travel}" effort="100" velocity="1"/>
  </joint>
  <joint name="right_finger_joint" type="prismatic">
    <parent link="palm"/><child link="right_finger"/>
    <origin xyz="0 {ry} {fzc}"/><axis xyz="0 -1 0"/>
    <limit lower="0" upper="{travel}" effort="100" velocity="1"/>
  </joint>
</robot>
""".format(i3=inertia.format(i=1e-3), i5=inertia.format(i=5e-5), palm=box % palm,
           finger=box % (fx, fy, fz), ly=-open_y, ry=open_y,
           fzc=-0.5 * (PALM_THICKNESS + finger_height), travel=0.5 * FINGER_OPEN_GAP)


def smoothstep(t0, t1, t):
  """C1 ramp from 0 (t <= t0) to 1 (t >= t1)."""
  s = min(max((t - t0) / (t1 - t0), 0.0), 1.0)
  return s * s * (3.0 - 2.0 * s)


class Schedule(object):
  """Scripted gripper motion."""
  phases = [  # name, start [s], end [s]
      ("settle", 0.0, 0.5),
      ("descend", 0.5, 1.3),
      ("close", 1.3, 2.1),
      ("lift", 2.1, 3.1),
      ("hold", 3.1, 5.1),
      ("transport", 5.1, 6.1),
      ("hold2", 6.1, 7.1),
      ("lower", 7.1, 7.9),
      ("release", 7.9, 8.5),
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
    """(x, z, finger) joint position targets at time t."""
    z = LIFT_START + (self.grasp_height - LIFT_START) * smoothstep(0.5, 1.3, t)
    z += LIFT_HEIGHT * (smoothstep(2.1, 3.1, t) - smoothstep(7.1, 7.9, t))
    x = TRANSPORT_X * smoothstep(5.1, 6.1, t)
    finger = self.finger_close * (smoothstep(1.3, 2.0, t) - smoothstep(7.9, 8.3, t))
    return x, z, finger


def box_signed_distance(points, box_pos, box_orn, half_extents):
  """Exact signed distance from world points (N x 3) to an oriented box."""
  rot = np.asarray(p.getMatrixFromQuaternion(box_orn)).reshape(3, 3)
  local = (points - np.asarray(box_pos)).dot(rot)  # world -> box frame
  q = np.abs(local) - np.asarray(half_extents)
  return np.linalg.norm(np.maximum(q, 0.0), axis=1) + np.minimum(q.max(axis=1), 0.0)


class Demo(object):

  def __init__(self, args):
    self.args = args
    self.finger_half = (0.5 * args.finger_width, 0.5 * FINGER_THICKNESS, 0.5 * args.finger_height)
    p.resetSimulation(p.RESET_USE_DEFORMABLE_WORLD)
    p.setGravity(0, 0, -9.81)
    p.setTimeStep(args.timestep)
    p.setPhysicsEngineParameter(numSolverIterations=args.iterations,
                                useSplitImpulse=args.split_impulse)
    p.setAdditionalSearchPath(pybullet_data.getDataPath())

    self.plane = p.loadURDF("plane.urdf")
    p.changeDynamics(self.plane, -1, lateralFriction=0.6)

    urdf_path = os.path.join(tempfile.mkdtemp(prefix="deformable_grasp_"), "gripper.urdf")
    with open(urdf_path, "w") as f:
      f.write(gripper_urdf(args.finger_width, args.finger_height))
    self.gripper = p.loadURDF(urdf_path, [0, 0, 0], useFixedBase=True)
    for link in (LEFT, RIGHT):
      p.changeDynamics(self.gripper, link, lateralFriction=args.finger_friction)
    p.resetJointState(self.gripper, LIFT_Z, LIFT_START)

    self.sphere = p.loadSoftBody("ball.vtk",
                                 basePosition=[0, 0, SPHERE_RADIUS + 0.002],
                                 scale=SPHERE_RADIUS,
                                 mass=args.sphere_mass,
                                 useNeoHookean=1,
                                 NeoHookeanMu=args.mu,
                                 NeoHookeanLambda=args.lam,
                                 NeoHookeanDamping=args.damping,
                                 collisionMargin=0.002,
                                 frictionCoeff=args.sphere_friction,
                                 useSelfCollision=0,
                                 useFaceContact=1)
    p.changeVisualShape(self.sphere, -1, rgbaColor=[0.95, 0.45, 0.2, 1])

    # the finger centres are at the sphere centre while grasping
    finger_centre_below_palm = 0.5 * (PALM_THICKNESS + args.finger_height)
    grasp_height = SPHERE_RADIUS + finger_centre_below_palm
    # close the fingers until they would squeeze the sphere by `squeeze` per side; the force
    # limit of the finger motors decides how far they actually get
    finger_close = 0.5 * FINGER_OPEN_GAP - (SPHERE_RADIUS - args.squeeze)
    self.schedule = Schedule(grasp_height, finger_close)

  def command(self, t):
    x, z, finger = self.schedule.targets(t)
    for joint, target in ((SLIDE_X, x), (LIFT_Z, z)):
      p.setJointMotorControl2(self.gripper, joint, p.POSITION_CONTROL, targetPosition=target,
                              force=500, positionGain=0.2, velocityGain=1.0)
    for joint in (LEFT, RIGHT):
      p.setJointMotorControl2(self.gripper, joint, p.POSITION_CONTROL, targetPosition=finger,
                              force=self.args.grip_force, positionGain=0.2, velocityGain=1.0,
                              maxVelocity=0.2)

  def nodes(self):
    return np.asarray(p.getMeshData(self.sphere, -1, flags=p.MESH_DATA_SIMULATION_MESH)[1])

  def finger_states(self):
    return [p.getLinkState(self.gripper, link, computeLinkVelocity=1, computeForwardKinematics=1)
            for link in (LEFT, RIGHT)]

  def camera_target(self):
    palm = p.getLinkState(self.gripper, LIFT_Z, computeForwardKinematics=1)[0]
    return [palm[0], 0, max(palm[2] - 0.06, 0.06)]


class VideoWriter(object):
  """Pipes RGBA frames into ffmpeg."""

  def __init__(self, path, width, height, fps):
    self.proc = subprocess.Popen([
        "ffmpeg", "-loglevel", "error", "-y", "-f", "rawvideo", "-pix_fmt", "rgba", "-s",
        "%dx%d" % (width, height), "-r", str(fps), "-i", "-", "-pix_fmt", "yuv420p", "-vcodec",
        "libx264", "-crf", "18", path
    ], stdin=subprocess.PIPE)

  def add(self, rgba):
    self.proc.stdin.write(rgba.tobytes())

  def close(self):
    self.proc.stdin.close()
    self.proc.wait()


def render(width, height, target):
  view = p.computeViewMatrixFromYawPitchRoll(cameraTargetPosition=target, distance=CAMERA[0],
                                             yaw=CAMERA[1], pitch=CAMERA[2], roll=0, upAxisIndex=2)
  proj = p.computeProjectionMatrixFOV(fov=40, aspect=float(width) / height, nearVal=0.01,
                                      farVal=5.0)
  img = p.getCameraImage(width, height, view, proj, renderer=p.ER_TINY_RENDERER,
                         lightDirection=[1.0, -0.4, 1.2], shadow=1)
  return np.asarray(img[2], dtype=np.uint8).reshape(height, width, 4)


def simulate(args, demo, video=None):
  """Runs the scripted grasp once and returns the per-step measurements."""
  dt = args.timestep
  frame_every = max(1, int(round(1.0 / (args.fps * dt))))
  prev_nodes = None
  samples = []
  for step in range(int(round(Schedule.duration / dt))):
    demo.command(step * dt)
    p.stepSimulation()
    t = (step + 1) * dt

    nodes = demo.nodes()
    fingers = demo.finger_states()
    finger_mid = 0.5 * (np.asarray(fingers[0][0]) + np.asarray(fingers[1][0]))
    gripper_vel = 0.5 * (np.asarray(fingers[0][6]) + np.asarray(fingers[1][6]))
    # deepest penetration of a sphere node into a finger
    penetration = max(
        float(-box_signed_distance(nodes, f[0], f[1], demo.finger_half).min()) for f in fingers)
    # RMS speed of the sphere nodes relative to the gripper
    rel_speed = float("nan")
    if prev_nodes is not None:
      rel_v = (nodes - prev_nodes) / dt - gripper_vel
      rel_speed = float(np.sqrt((rel_v * rel_v).sum(axis=1).mean()))
    prev_nodes = nodes
    joints = [p.getJointState(demo.gripper, j) for j in (LEFT, RIGHT)]
    samples.append({
        "t": t,
        "phase": demo.schedule.phase(t - 1e-9),
        "com_z": float(nodes[:, 2].mean()),
        "rel_z": float(nodes[:, 2].mean() - finger_mid[2]),
        "penetration": max(penetration, 0.0),
        "rel_speed": rel_speed,
        "grip_force": 0.5 * (abs(joints[0][3]) + abs(joints[1][3])),
        "finger_q": 0.5 * (joints[0][0] + joints[1][0]),
    })

    if video is not None and step % frame_every == 0:
      video.add(render(args.width, args.height, demo.camera_target()))
    if not args.direct:
      time.sleep(dt)
  return samples


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


def compute_metrics(samples):
  by_phase = {}
  for s in samples:
    by_phase.setdefault(s["phase"], []).append(s)
  held = by_phase["lift"] + by_phase["hold"] + by_phase["transport"] + by_phase["hold2"]
  # the gripper does not move during the second half of the hold phases
  static = (by_phase["hold"][len(by_phase["hold"]) // 2:] +
            by_phase["hold2"][len(by_phase["hold2"]) // 2:])
  ref = by_phase["hold"][0]["rel_z"]  # the sphere is off the ground when the hold starts
  return {
      "lifted": bool(by_phase["hold2"][-1]["com_z"] > SPHERE_RADIUS + 0.5 * LIFT_HEIGHT),
      "slip_hold_and_transport_mm": 1000.0 * (ref - by_phase["hold2"][-1]["rel_z"]),
      "slip_max_mm": 1000.0 * max(ref - s["rel_z"] for s in by_phase["hold"] + by_phase["transport"] +
                                  by_phase["hold2"]),
      "jitter_rms_mm_per_s": 1000.0 * float(np.sqrt(np.mean([s["rel_speed"]**2 for s in static]))),
      "jitter_max_mm_per_s": 1000.0 * max(s["rel_speed"] for s in static),
      "penetration_max_mm": 1000.0 * max(s["penetration"] for s in held),
      "grip_force_mean_N": float(np.mean([s["grip_force"] for s in static])),
      "grip_force_std_N": float(np.std([s["grip_force"] for s in static])),
  }


def parse_args(argv=None):
  parser = argparse.ArgumentParser(description=__doc__,
                                   formatter_class=argparse.RawDescriptionHelpFormatter)
  parser.add_argument("--direct", action="store_true", help="run headless (no GUI)")
  parser.add_argument("--timestep", type=float, default=1.0 / 240.0)
  parser.add_argument("--iterations", type=int, default=50, help="solver iterations")
  parser.add_argument("--split-impulse", type=int, default=1, help="use split impulse (0/1)")
  parser.add_argument("--grip-force", type=float, default=10.0, help="finger motor force [N]")
  parser.add_argument("--squeeze", type=float, default=0.015,
                      help="finger travel past the undeformed sphere surface, per side [m]")
  parser.add_argument("--finger-width", type=float, default=0.03, help="[m]")
  parser.add_argument("--finger-height", type=float, default=0.07, help="[m]")
  parser.add_argument("--sphere-mass", type=float, default=0.1, help="[kg]")
  parser.add_argument("--mu", type=float, default=2000.0, help="Neo-Hookean mu [Pa]")
  parser.add_argument("--lam", type=float, default=3000.0, help="Neo-Hookean lambda [Pa]")
  parser.add_argument("--damping", type=float, default=0.02, help="Neo-Hookean damping")
  parser.add_argument("--sphere-friction", type=float, default=1.0)
  parser.add_argument("--finger-friction", type=float, default=1.0)
  parser.add_argument("--metrics", help="write metrics and time series (json) to this file")
  parser.add_argument("--video", help="write an mp4 video to this file (needs ffmpeg)")
  parser.add_argument("--width", type=int, default=640)
  parser.add_argument("--height", type=int, default=480)
  parser.add_argument("--fps", type=int, default=30)
  return parser.parse_args(argv)


if __name__ == "__main__":
  run(parse_args())
