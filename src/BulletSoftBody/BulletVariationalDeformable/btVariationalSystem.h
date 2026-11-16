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

#ifndef BT_VARIATIONAL_SYSTEM_H
#define BT_VARIATIONAL_SYSTEM_H

#include "btVariationalMath.h"
#include "btVariationalElasticity.h"
#include "LinearMath/btTransform.h"

class btSoftBody;
class btCollisionObject;
class btCollisionShape;
class btVariationalContact;
struct btSoftBodyWorldInfo;

///Linear tetrahedron. The deformation gradient is F = sum_j x_j * m_g[j]^T.
struct btVariationalTet
{
	int m_n[4];
	btVec3d m_g[4];
	double m_restVolume;
	int m_body;
};

struct btVariationalSpring
{
	int m_n[2];
	double m_restLength;
	double m_stiffness;
	int m_body;
};

///Surface triangle and edge of the collision mesh (indices into the global node array)
struct btVariationalTriangle
{
	int m_n[3];
	int m_body;
};

struct btVariationalEdge
{
	int m_n[2];
	int m_body;
};

///Ties a node to a world-space target (m_nodeB < 0) or to another node.
///Soft attachments are springs. Hard attachments are enforced by the augmented Lagrangian of
///"Augmented Vertex Block Descent" (Giles, Diaz, Yuksel 2025): m_lambda is the dual variable and
///m_penalty the ramped penalty stiffness.
struct btVariationalAttachment
{
	int m_nodeA;
	int m_nodeB;
	btVec3d m_target;
	btVec3d m_prevTarget;  // target at the end of the previous step
	double m_stiffness;
	bool m_hard;
	btVec3d m_lambda;
	double m_penalty;
	btVec3d m_c0;  // constraint error at the start of the step (used by the AVBD stabilization)

	btVariationalAttachment() : m_nodeA(-1), m_nodeB(-1), m_stiffness(0), m_hard(false), m_penalty(0) {}
};

///Per soft body settings and bookkeeping
struct btVariationalBody
{
	btSoftBody* m_psb;
	const btVariationalMaterial* m_material;  // not owned
	double m_springStiffness;                 // stiffness of structural links (bodies without tetrahedra)
	double m_bendingStiffness;                // stiffness of bending links
	double m_friction;
	int m_firstNode, m_numNodes;
	int m_firstTet, m_numTets;
	int m_firstSpring, m_numSprings;
	btAlignedObjectArray<int> m_colliders;  // rigid colliders overlapping this body in the current step
};

///A rigid (static, kinematic or dynamic) Bullet collision object seen as an implicit surface
struct btVariationalCollider
{
	enum Type
	{
		PLANE,
		BOX,
		SPHERE,
		CAPSULE,
		GENERIC  // any other convex shape, evaluated through Bullet's sparse signed distance field
	};
	const btCollisionObject* m_object;
	const btCollisionShape* m_shape;
	btSoftBodyWorldInfo* m_worldInfo;
	int m_type;
	btMat3d m_R;  // rotation and origin of the shape at the end of the step
	btVec3d m_p;
	btVec3d m_center;  // origin of the rigid body the shape belongs to (pivot of m_angVel and of the reaction torque)
	btVec3d m_linVel, m_angVel;
	double m_friction;
	btVec3d m_halfExtents;
	double m_radius, m_halfHeight;
	int m_axis;
	btVec3d m_planeNormal;
	double m_planeConstant;
	btVec3d m_impulse, m_torqueImpulse;  // reaction accumulated for two-way coupling

	///signed distance (negative inside) and outward unit normal at a world-space point
	void signedDistance(const btVec3d& x, double& dist, btVec3d& normal) const;
	bool isConvexExact() const { return m_type != GENERIC; }
	btVec3d velocityAt(const btVec3d& x) const { return m_linVel + m_angVel.cross(x - m_center); }
};

///Element-by-element storage of the (filtered) Hessian of the incremental potential:
/// H = sum_i P_i^T m_blocks[i] P_i  +  blockdiag(m_nodeBlocks)
///Stencils have up to four nodes; unused slots are -1.
struct btVariationalHessianStencil
{
	int m_n[4];
	btMat12d m_H;
};

class btVariationalHessian
{
public:
	btAlignedObjectArray<btVariationalHessianStencil> m_stencils;
	btAlignedObjectArray<btMat3d> m_nodeBlocks;

	void reset(int numNodes);
	///y = H * x; rows and columns of fixed nodes are treated as identity
	void multiply(const btVec3d* x, btVec3d* y, const char* fixed) const;
	///3x3 diagonal blocks, used by the block Jacobi preconditioner
	void diagonalBlocks(btMat3d* D) const;
};

///Preconditioner interface for the conjugate gradient solver, z = M^-1 r
class btVariationalPreconditioner
{
public:
	virtual ~btVariationalPreconditioner() {}
	virtual void apply(const btVec3d* r, btVec3d* z) const = 0;
};

///Conjugate gradients on the element-by-element Hessian. Returns the number of iterations.
///If the matrix turns out not to be positive definite along a search direction the iteration stops
///with the iterate reached so far and *indefinite is set.
int btVariationalConjugateGradient(const btVariationalHessian& H, const btVec3d* b, btVec3d* x, const char* fixed, int numNodes,
								   const btVariationalPreconditioner* preconditioner, double relativeTolerance, int maxIterations,
								   bool* indefinite = 0);

class btVariationalSystem;

///Deformable-deformable (and self) contact potential that can be plugged into the system.
///See btVariationalIPC for the implementation of "Incremental Potential Contact".
class btVariationalContact
{
public:
	virtual ~btVariationalContact() {}
	///called once per time step before the solve, with x = start-of-step positions
	virtual void beginStep(const btVariationalSystem& system) = 0;
	virtual double energy(const btVec3d* x) const = 0;
	virtual void addGradient(const btVec3d* x, btVec3d* g) const = 0;
	virtual void addHessian(const btVec3d* x, btVariationalHessian& H, btVariationalHessianFilter filter) const = 0;
	virtual void addVertexForceAndHessian(int node, const btVec3d* x, btVec3d& force, btMat3d& hess) const = 0;
	///Largest alpha in [0,1] for which moving only this node by alpha * dx is free of intersections.
	///With conservative set the node also enters the contact zone only gradually (see btAdditiveCCD).
	virtual double vertexMaxStepSize(int node, const btVec3d* x, const btVec3d& dx, bool conservative) const = 0;
	///Called by solvers that move one node at a time, before a sweep in which no node will travel
	///further than radius from x (vertexMaxStepSize enforces that bound). If direction is given, node i
	///will instead only move along the segment from x_i to x_i + direction_i.
	///Returns the radius that is actually in effect (it may be enlarged).
	virtual double prepareLocalSteps(const btVec3d* x, double radius, const btVec3d* direction = 0)
	{
		(void)x;
		(void)direction;
		return radius;
	}
	///Gives the contact model the chance to adapt its stiffness to the current iterate. Returns true
	///if the energy changed as a result.
	virtual bool updateStiffness(const btVec3d* x)
	{
		(void)x;
		return false;
	}
	///largest alpha in [0,1] for which x + alpha * p is free of intersections
	virtual double maxStepSize(const btVec3d* x, const btVec3d* p) const = 0;
};

///The discrete problem all integrators of this module work on.
///
///One implicit Euler step is the minimization of the incremental potential
///
///   E(x) = 1/(2 h^2) |x - xTilde|^2_M  +  Psi_elastic(x)  +  E_springs(x)  +  E_attach(x)  +  E_contact(x)
///
///with xTilde = x_n + h v_n + h^2 g. This class owns the state, the mesh and the energy terms and
///offers the evaluations the different solvers need: the energy, its gradient, the eigenvalue
///filtered Hessian (element by element), per-vertex 3x3 blocks, and a feasible step bound.
class btVariationalSystem
{
public:
	enum ColliderModel
	{
		COLLIDER_PENALTY = 0,  // quadratic penalty on the signed distance (allows small penetration)
		COLLIDER_BARRIER       // IPC log barrier on the distance (never penetrates, needs a feasible start)
	};

	// ---- parameters
	btVec3d m_gravity;
	double m_dt;
	btVariationalHessianFilter m_hessianFilter;
	double m_damping;            // [1/s] damping of the deformation modes, rigid body motion is preserved
	int m_colliderModel;
	double m_colliderStiffness;  // penalty stiffness of a node, relative to its mass / h^2
	double m_colliderMargin;     // distance at which rigid contact starts to act
	double m_frictionEpsilon;    // [m/s] slip velocity below which friction behaves like static friction
	bool m_twoWayCoupling;       // accumulate reaction impulses for dynamic rigid bodies
	bool m_evaluateColliders;    // include the rigid collider terms in the evaluations (integrators that treat rigid contact by projection switch this off)

	// ---- state (one entry per node)
	btVec3dArray m_x, m_xPrev, m_v, m_vPrev, m_xTilde, m_fixedTarget;
	btAlignedObjectArray<double> m_mass;
	btAlignedObjectArray<char> m_fixed;
	btAlignedObjectArray<int> m_nodeBody;

	// ---- mesh
	btAlignedObjectArray<btVariationalBody> m_bodies;
	btAlignedObjectArray<btVariationalTet> m_tets;
	btAlignedObjectArray<btVariationalSpring> m_springs;
	btAlignedObjectArray<btVariationalTriangle> m_triangles;
	btAlignedObjectArray<btVariationalEdge> m_edges;
	btAlignedObjectArray<char> m_surfaceNode;
	btAlignedObjectArray<btVariationalAttachment> m_attachments;
	btAlignedObjectArray<btVariationalCollider> m_colliders;

	// node -> incident elements (compressed lists); tets are stored as 4 * tet + localIndex
	btAlignedObjectArray<int> m_nodeTetStart, m_nodeTets;
	btAlignedObjectArray<int> m_nodeSpringStart, m_nodeSprings;
	btAlignedObjectArray<int> m_nodeAttachStart, m_nodeAttach;

	btVariationalContact* m_contact;  // optional deformable-deformable contact (IPC), not owned
	int m_topologyVersion;            // incremented whenever the mesh is rebuilt

private:
	// lagged friction data of the node-collider pairs, indexed [body collider slot][node - firstNode]
	struct ColliderLag
	{
		btVec3d m_normal;
		double m_normalForce;
	};
	btAlignedObjectArray<ColliderLag> m_colliderLag;
	btAlignedObjectArray<int> m_colliderLagStart;  // per body: offset of its first slot

	const ColliderLag* lagOf(int body, int slot, int node) const
	{
		if (m_colliderLag.size() == 0) return 0;
		return &m_colliderLag[m_colliderLagStart[body] + slot * m_bodies[body].m_numNodes + (node - m_bodies[body].m_firstNode)];
	}
	double colliderKernel(double d, double stiffness, double& dE, double& d2E) const;
	double colliderNodeMaxStep(const btVariationalCollider& col, const btVec3d& xi, const btVec3d& pi, double alpha, bool conservative = false) const;
	void addColliderNode(int body, int node, const btVec3d& xi, const btVec3d& xPrev, double& energy, btVec3d* grad, btMat3d* hess) const;

public:
	btVariationalSystem();

	// ---- construction
	void clear();
	///registers a soft body; nodes, tetrahedra (or links) and faces are converted in buildTopology()
	int addBody(btSoftBody* psb, const btVariationalMaterial* material, double springStiffness, double bendingStiffness);
	int findBody(const btSoftBody* psb) const;
	void buildTopology();
	///adds an attachment and returns its index; rebuilds the node lookup lazily
	int addAttachment(const btVariationalAttachment& attachment);
	void clearAttachments();
	void rebuildAttachmentLookup();

	// ---- exchange with Bullet
	void readFromSoftBodies();  // picks up positions/velocities changed outside the solver and kinematic targets
	void writeToSoftBodies();
	void clearColliders();
	int addCollider(const btCollisionObject* object, const btCollisionShape* shape, const btTransform& endOfStepTransform, btSoftBodyWorldInfo* worldInfo);

	// ---- time stepping
	void beginStep(double dt);
	void endStep();
	///recomputes the lagged normal forces and contact normals the friction terms use
	void updateFrictionLag(const btVec3d* x);

	// ---- evaluations of the incremental potential
	int numNodes() const { return m_x.size(); }
	double energy(const btVec3d* x) const;
	void gradient(const btVec3d* x, btVec3d* g) const;
	///element-by-element Hessian with the given eigenvalue filter (m_hessianFilter by default)
	void hessian(const btVec3d* x, btVariationalHessian& H) const { hessian(x, H, m_hessianFilter); }
	void hessian(const btVec3d* x, btVariationalHessian& H, btVariationalHessianFilter filter) const;
	///descent force (-gradient) and 3x3 Hessian block of a single node, all other nodes held fixed.
	///elasticHess, if given, receives the part of the block that comes from the elements alone
	///(vertex block descent uses it for stiffness proportional damping).
	void vertexForceAndHessian(int node, const btVec3d* x, btVec3d& force, btMat3d& hess, btMat3d* elasticHess = 0) const;
	///Largest alpha in [0,1] such that moving only this node by alpha * dx keeps the barrier terms
	///feasible. With conservative set, a node outside a barrier zone stops at half the activation
	///distance and a node inside closes at most a tenth of its gap, so contact is entered gradually.
	double vertexMaxStepSize(int node, const btVec3d* x, const btVec3d& dx, bool conservative = false) const;
	///largest step alpha in [0,1] such that x + alpha * p stays feasible for the barrier terms
	double maxStepSize(const btVec3d* x, const btVec3d* p) const;
	///infinity norm of the gradient over the free nodes
	double gradientNorm(const btVec3d* g) const;

	// ---- helpers shared by the integrators
	btMat3d deformationGradient(const btVariationalTet& tet, const btVec3d* x) const
	{
		btMat3d F;
		for (int j = 0; j < 4; j++) F.addScaledOuter(1.0, x[tet.m_n[j]], tet.m_g[j]);
		return F;
	}
	double elasticEnergy(const btVec3d* x) const;
	double kineticEnergy() const;
	double totalVolume(const btVec3d* x) const;
	double colliderPenetration(const btVec3d* x) const;  // deepest penetration into any collider
};

#endif  // BT_VARIATIONAL_SYSTEM_H
