# Variational deformable bodies

Optimization-based time integration of deformable bodies for Bullet. Finite element solids (linear
tetrahedra) and mass-spring cloth are advanced by minimizing the incremental potential of implicit
Euler,

    E(x) = 1/(2 h^2) |x - x~|^2_M + elastic energy + attachments + contact + friction,
    x~ = x_n + h v_n + h^2 g,

and eight papers on how to pose and solve that problem are implemented behind one interface, so the
methods can be exchanged, compared and combined on the same scene. The module plugs into
`btDeformableMultiBodyDynamicsWorld` as a replacement for `btDeformableBodySolver`.

## Papers

| | Paper | Published | Code that was consulted | Implemented in |
|-|-------|-----------|-------------------------|----------------|
| 1 | *Fast Simulation of Mass-Spring Systems*, Liu, Bargteil, O'Brien, Kavan | SIGGRAPH Asia 2013 | [sam007961/FastMassSpring](https://github.com/sam007961/FastMassSpring) (third party) | `btVariationalProjectiveDynamics`, `btVariationalLaplacian` |
| 2 | *Projective Dynamics: Fusing Constraint Projections for Fast Simulation*, Bouaziz, Martin, Liu, Kavan, Pauly | SIGGRAPH 2014 | [EPFL-LGG/ShapeOp](https://github.com/EPFL-LGG/ShapeOp) | `btVariationalProjectiveDynamics`, `btVariationalARAP` |
| 3 | *Quasi-Newton Methods for Real-Time Simulation of Hyperelastic Materials*, Liu, Bouaziz, Kavan | ACM TOG 36(3), 2017 (presented at SIGGRAPH 2017) | [ltt1598/Quasi-Newton-Methods-for-Real-time-Simulation-of-Hyperelastic-Materials](https://github.com/ltt1598/Quasi-Newton-Methods-for-Real-time-Simulation-of-Hyperelastic-Materials) | `btVariationalQuasiNewton` |
| 4 | *Stable Neo-Hookean Flesh Simulation*, Smith, de Goes, Kim | ACM TOG 37(2), 2018 (presented at SIGGRAPH 2018) | the authors' supplemental code, [theodorekim/HOBAKv1](https://github.com/theodorekim/HOBAKv1) | `btVariationalStableNeoHookean`, `btVariationalProjectedNewton` |
| 5 | *Stabler Neo-Hookean Simulation: Absolute Eigenvalue Filtering for Projected Newton*, Chen, Liu, Levin, Zheng, Jacobson | SIGGRAPH 2024 | [honglin-c/abs-psd](https://github.com/honglin-c/abs-psd) | `BT_VARIATIONAL_FILTER_ABS` in `btVariationalMath.h`, selectable for every Hessian assembly |
| 6 | *Incremental Potential Contact: Intersection- and Inversion-free, Large-Deformation Dynamics*, Li, Ferguson, Schneider, Langlois, Zorin, Panozzo, Jiang, Kaufman | SIGGRAPH 2020 | [ipc-sim/IPC](https://github.com/ipc-sim/IPC), [ipc-sim/ipc-toolkit](https://github.com/ipc-sim/ipc-toolkit) | `btVariationalIPC`, `btVariationalDistance.h` |
| 7 | *Vertex Block Descent*, Chen, Liu, Yang, Yuksel | SIGGRAPH 2024 | [AnkaChan/Gaia](https://github.com/AnkaChan/Gaia) | `btVariationalVertexBlockDescent` |
| 8 | *Augmented Vertex Block Descent*, Giles, Diaz, Yuksel | SIGGRAPH 2025 | [savant117/avbd-demo2d](https://github.com/savant117/avbd-demo2d), [avbd-demo3d](https://github.com/savant117/avbd-demo3d) | `btVariationalAugmentedLagrangian`, used by `btVariationalVertexBlockDescent` and `btVariationalHybrid` |

The code here is written against the papers for Bullet's data structures. The repositories were read
for details the papers leave open (barrier stiffness adaptation, initialization and acceleration of
vertex block descent, the dual update of the augmented Lagrangian); no code was copied from them.

## Structure

| File | Contents |
|------|----------|
| `btVariationalDeformableSolver` | `btDeformableBodySolver` subclass: collects soft bodies, colliders, attachments and mouse picking from the Bullet world, runs one optimization step per time step, writes the result back |
| `btVariationalSystem` | The incremental potential: energy, gradient, Hessian (per-element stencils), per-vertex force and 3x3 block, step size limits. Rigid colliders as signed distance functions with penalty or barrier contact and lagged friction |
| `btVariationalElasticity` | Materials: Stable Neo-Hookean with its analytic eigen-system, and as-rigid-as-possible (the projective dynamics strain energy) |
| `btVariationalProjectedNewton` | Newton with eigenvalue-filtered Hessians (clamp or absolute value), preconditioned conjugate gradients and a line search that starts at the largest intersection-free step |
| `btVariationalProjectiveDynamics` | Local/global solver with a prefactorized system matrix |
| `btVariationalQuasiNewton` | L-BFGS with the prefactorized matrix as initial Hessian |
| `btVariationalVertexBlockDescent` | Colored Gauss-Seidel with one 3x3 Newton step per vertex, adaptive initialization, Chebyshev acceleration, augmented Lagrangian for hard attachments |
| `btVariationalIPC` | Barrier contact between deformable surfaces with continuous collision detection and friction |
| `btVariationalHybrid` | Chooses among the above per time step, see below |
| `btVariationalLaplacian`, `btVariationalSparse` | The constant matrix M/h^2 + L of projective dynamics and a sparse Cholesky factorization (nested dissection, up-looking) |
| `btVariationalDistance.h` | Point-triangle and edge-edge distances with exact derivatives, additive continuous collision detection |
| `btVariationalMath.h` | Double precision 3x3 / 9x9 / 12x12 algebra, rotation-variant SVD, symmetric eigen-solver, eigenvalue filters |

The solver state is kept in double precision whatever `btScalar` is: barrier contact works with
gaps that single precision cannot resolve.

## Use

    btVariationalDeformableSolver* solver = new btVariationalDeformableSolver();
    btDeformableMultiBodyConstraintSolver* sol = new btDeformableMultiBodyConstraintSolver();
    sol->setDeformableSolver(solver);
    world = new btDeformableMultiBodyDynamicsWorld(dispatcher, broadphase, sol, collisionConfiguration, solver);

    btVariationalStableNeoHookean rubber(1e6 /* Young's modulus */, 0.45 /* Poisson's ratio */);
    btSoftBody* psb = btVariationalHelpers::createTetGrid(world->getWorldInfo(), center, halfExtents, 8, 2, 2, 1000);
    world->addSoftBody(psb);
    solver->setMaterial(psb, &rubber);              // any btSoftBody with tetrahedra
    solver->setSpringStiffness(cloth, 4000, 20);    // bodies with links only are mass-spring systems

    btVariationalHybrid integrator;                 // or any other btVariationalIntegrator
    solver->setIntegrator(&integrator);
    btVariationalIPC contact;                       // optional: contact between deformables
    solver->getSystem().m_contact = &contact;
    solver->addAttachment(psb, node, target, 1e5, true);   // hard attachment to a world point

Nodes with zero mass are kinematic and follow `btSoftBody::Node::m_x`. Static, kinematic and dynamic
rigid bodies of the world act as colliders. Materials, integrators and the contact object are not
owned by the solver.

## Hybrid integrator

`btVariationalHybrid` is the result of measuring where each method wins (next section):

- Whenever contact is possible in a step, it takes projected Newton steps on the full potential
  (papers 4, 5, 6). Newton is the only method here that converges in a few iterations with active
  barriers.
- A step without any possible contact and without hard attachments goes to L-BFGS on the
  prefactorized matrix (papers 2, 3): no Hessian is assembled. If it runs out of its iteration
  budget, Newton continues from where it stopped.
- From 600 nodes on, the conjugate gradient solves inside Newton are preconditioned with the
  prefactorized projective dynamics matrix, combined with block Jacobi for the nodes whose stiffness
  that matrix does not know about (contact, joints).
- Hard attachments are enforced by the augmented Lagrangian of paper 8 around the solve: after each
  solve the constraint force moves into the dual variable and, while joints are still open, the
  step is solved again from the current iterate.
- A warm start with a few sweeps of vertex block descent (paper 7) is available but off by default:
  it saved under ten percent without contact and cost time with contact.

## Examples

ExampleBrowser, group "Variational Deformable Body" (`examples/VariationalDeformableDemo`):

| Example | Shows |
|---------|-------|
| Projective Dynamics | A cloth and two solids advanced by the local/global solver; slider for the iteration count |
| Quasi-Newton Hyperelastic | A nearly incompressible Stable Neo-Hookean bar stretched and compressed; slider for the L-BFGS history, where 0 is the plain constant-matrix step of projective dynamics |
| Stable Neo-Hookean | A twisted bar, a cube squashed flat and a cube with scrambled nodes recover; clamp and absolute-value filter |
| Incremental Potential Contact | Soft cubes and slabs pile up without interpenetration, with friction |
| Vertex Block Descent | A batch of soft cubes rains down, no global linear system |
| Augmented Vertex Block Descent | Two chains with a hundredfold heavier block at the end: hard joints stay closed, springs come apart |
| Combined: Hybrid Solver | A tower of cubes, a wrecking ball on hard joints, a cloth draped over a cube and a twisted bar in one scene and one potential; slider to switch to the individual methods |

Every example shows the solver, iteration count, time per step and residual on screen. Set the
environment variable `BT_VARIATIONAL_DEMO_LOG` to also print them, with per-body diagnostics, to the
console. An example can be started directly with
`App_ExampleBrowser --start_demo_name="Combined: Hybrid Solver"`.

## Tests

`test/BulletSoftBody/test_btVariationalDeformable.cpp` (built with `BUILD_UNIT_TESTS`, run with
`ctest -R Variational` or directly) checks, among others:

- stress and Hessian of the materials, gradient and Hessian of the whole potential, the distance
  functions and the contact barrier against finite differences,
- a hanging bar against linear elasticity (elongation within 3 %, lateral contraction gives the
  Poisson's ratio within 0.02), exact free fall and momentum, friction on an incline against the
  analytic sliding velocity, recovery from flattened and scrambled states,
- projective dynamics, L-BFGS and vertex block descent converge to the Newton solution of the same
  step, and the energy of projective dynamics and of vertex block descent never increases from one
  iteration to the next,
- stacked bodies never intersect (verified by an independent edge-triangle test) with Newton,
  vertex block descent and the hybrid integrator,
- hard attachments close to a fraction of a millimetre where springs of the same stiffness stretch
  by centimetres,
- a Bullet world with the solver steps, lands a body on a rigid box, keeps an attachment, and can be
  torn down.

## Measurements

`test/BulletSoftBody/benchmark_btVariationalDeformable.cpp` runs every integrator on three scenes
that each stress a different part of the problem (`Benchmark_btVariationalDeformable 2` and `4` for
the two mesh resolutions below). Time step 1/60 s, Young's modulus 1 MPa, Poisson's ratio 0.45,
Newton tolerance 0.01 m/s, one thread of a shared 4-core virtual machine (timings repeat to within
10 to 30 %). The residual is the largest force imbalance left at the end of a step.

**Twisted bar, no contact** (clamped at both ends, one end turned by up to a full turn)

| ms per step | 272 nodes | 1188 nodes | Mean residual at 1188 nodes |
|-------------|----------:|-----------:|----------------------------:|
| Newton | 18.7 | 222.6 | 13 N |
| Quasi-Newton (L-BFGS) | 6.0 | 41.4 | 3.8 N |
| Vertex block descent, 20 sweeps | 3.7 | 22.2 | 1300 N |
| Vertex block descent, 100 sweeps | 16.2 | 98.6 | 250 N |
| Hybrid (takes the quasi-Newton path) | 5.8 | 40.7 | 3.8 N |
| Hybrid without the fast path (Newton, Laplacian preconditioner from 600 nodes) | 20.2 | 166.2 | 10 N |

**Pile of four cubes dropped on each other, barrier contact with friction**

| ms per step | 256 nodes | 864 nodes | Outcome |
|-------------|----------:|----------:|---------|
| Newton | 10.9 | 64.1 | pile stands, 1.5 to 1.7 iterations per step |
| Quasi-Newton (L-BFGS) | 132.1 | 545.9 | pile stands, 87 to 114 iterations per step |
| Vertex block descent, 20 sweeps | 20.7 | 167.8 | pile slumps to the ground |
| Vertex block descent, 100 sweeps | 104.2 | 337.4 | stands at 256 nodes, slumps at 864 |
| Hybrid (takes the Newton path once contact is possible) | 11.4 | 60.5 | pile stands |

No step of any run had an intersection.

**Chain: four light blocks and a 100 times heavier one on hard attachments, swung sideways** (135 nodes)

| | ms per step | Largest joint error |
|-|------------:|--------------------:|
| Newton, attachments as springs | 2.5 | 103 mm |
| Quasi-Newton, attachments as springs | 0.2 | 103 mm |
| Augmented vertex block descent, 20 sweeps | 1.1 | 1.1 mm |
| Augmented vertex block descent, 100 sweeps | 5.5 | 0.08 mm |
| Hybrid (Newton inside the augmented Lagrangian) | 3.7 | 0.34 mm (0.12 mm at the end) |

What the numbers say:

- No method wins everywhere, which is the reason for the hybrid integrator.
- Without contact the constant matrix of projective dynamics is a good model of the Hessian:
  L-BFGS on it is 3 to 5 times faster than Newton at the same tolerance. With barrier contact the
  same method is about ten times slower than Newton, because the matrix knows nothing about the
  barrier.
- Vertex block descent is the cheapest per step and never blows up or intersects, but a fixed
  number of sweeps is far from converged on stiff materials, and further the finer the mesh: with
  20 sweeps the residual is two orders of magnitude above Newton's, and stacked bodies sag under
  their own weight until the pile is flat.
- Newton is the one method that takes barrier contact in its stride.
- For hard joints the augmented Lagrangian is what matters: with it the error drops from 10 cm
  to a fraction of a millimetre, whether the inner solver is vertex block descent or Newton. On
  this small scene augmented vertex block descent is the better deal: a third of the cost at 1 mm,
  or four times the accuracy at one and a half times the cost. The hybrid pays for full Newton
  solves.
- The Laplacian preconditioner cuts a quarter off the Newton time on the 1188 node bar. On the 864
  node pile the gain (6 %) is within the timing noise, and forced on for the 529 node combined
  example it made the step 10 % slower, hence the size threshold.

**Eigenvalue filters of projected Newton** (`Benchmark_btVariationalDeformable filters`)

Newton iterations for one quasi-static step of a free cube that starts as its rest shape scaled
uniformly, i.e. with a large volume change (stopped at 1000):

| Poisson's ratio | Initial scale | Clamp | Absolute value |
|----------------:|--------------:|------:|---------------:|
| 0.3 | 0.2 | 222 | 15 |
| 0.3 | 0.5 | 934 | 9 |
| 0.3 | 2 | over 1000 | 66 |
| 0.3 | 5 | 896 | 179 |
| 0.495 | 0.2 | over 1000 | over 1000 (16 when the exact Hessian is tried first) |
| 0.495 | 0.5 | over 1000 | 8 |
| 0.495 | 2 | over 1000 | 60 |
| 0.495 | 5 | over 1000 | over 1000 |

Newton iterations per step for a bar twisted by one and a half turns at 60 steps per second:

| | Filter always applied | Exact Hessian tried first |
|-|----------------------:|--------------------------:|
| Clamp | 14.0 | 9.4 |
| Absolute value | 20.3 | 11.6 |

In the regime paper 5 is about, a start far from equilibrium in volume, its absolute-value filter
converges one to two orders of magnitude faster than clamping in most of the cases. For motion that
is driven through boundary conditions at an ordinary time step, clamping needed fewer iterations
here, and trying the unfiltered Hessian first (falling back to the filter on negative curvature)
helped both. The defaults are therefore clamping with the exact Hessian tried first;
`btVariationalSystem::m_hessianFilter = BT_VARIATIONAL_FILTER_ABS` selects the other filter.

## Limitations

- Elements are linear tetrahedra and springs. There are no shell or rod elements; cloth bending is
  modelled by bending springs.
- Time integration is implicit Euler, which damps motion: a tumbling soft block loses a noticeable
  part of its spin within a second.
- Rigid colliders are planes, boxes, spheres and capsules (exact signed distances) and other convex
  shapes through Bullet's sparse signed distance field. Triangle meshes and other concave shapes are
  ignored. Dynamic rigid bodies are colliders that move; passing the contact reaction back to them
  (`setTwoWayCoupling`) is an explicit impulse at the end of the step and is untested.
- Barrier contact needs an intersection-free start with a positive distance. Bodies that touch
  where they are joined have to share a collision group (`btVariationalIPC::setCollisionGroup`).
  Friction is lagged once per time step.
- Projective dynamics treats rigid colliders by projection and has no contact between deformables.
- Vertex block descent with a fixed sweep count does not converge on stiff or finely meshed
  scenes; the result is stable and intersection-free but softer than the material (see the pile
  above).
- Everything is single threaded, apart from the optional `btParallelFor` over the vertices of a
  color in vertex block descent. The combined example (529 nodes, up to several hundred active
  contact pairs) takes roughly 50 to 150 ms per step on the machine the measurements were made on,
  so it does not run in real time there.
- The solver is reachable from C++ only; it is not wired into the physics server or pybullet, and
  it has no serialization.
- `btSoftBody` features outside nodes, links, tetrahedra and faces (anchors, clusters, pose
  matching, aerodynamics) are not used; Bullet's sleeping is disabled for bodies handled here.
