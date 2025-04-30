import pybullet as p
import pybullet_data
import time

# Connect to PyBullet physics server with GUI
physicsClient = p.connect(p.GUI)

# Add the PyBullet data path for loading URDFs
p.setAdditionalSearchPath(pybullet_data.getDataPath())

# Set gravity
p.setGravity(0, 0, -9.81)

print("PyBullet connected and gravity set.")

# Load plane URDF
planeId = p.loadURDF("plane.urdf")

# Load sphere URDF
startPos = [0, 0, 1]  # Initial position (x, y, z)
startOrientation = p.getQuaternionFromEuler([0, 0, 0])
sphereId = p.loadURDF("sphere_small.urdf", startPos, startOrientation)
p.changeDynamics(sphereId, -1, restitution=0.9) # -1 refers to the base link

print("Plane and sphere loaded.")

# Simulation loop - run until the user closes the window
print("Simulation running. Close the PyBullet window to stop.")
while p.isConnected():
    p.stepSimulation()
    time.sleep(1./240.) # Simulate at 240Hz

# Disconnect from PyBullet
p.disconnect()

print("Simulation finished.")
