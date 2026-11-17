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

///Tests of the variational deformable module (src/BulletSoftBody/BulletVariationalDeformable):
///derivatives against finite differences, the solvers against analytic solutions and against each
///other, and the invariants the methods promise (no intersections, closed joints, monotone energy).

#include <gtest/gtest.h>
#include <math.h>

#include "btBulletDynamicsCommon.h"
#include "BulletSoftBody/btDeformableMultiBodyDynamicsWorld.h"
#include "BulletSoftBody/btDeformableMultiBodyConstraintSolver.h"
#include "BulletSoftBody/btSoftBody.h"
#include "BulletSoftBody/btSoftBodyHelpers.h"
#include "BulletSoftBody/btSoftBodyRigidBodyCollisionConfiguration.h"
#include "BulletSoftBody/BulletVariationalDeformable/btVariationalDeformableSolver.h"
#include "BulletSoftBody/BulletVariationalDeformable/btVariationalDistance.h"
#include "BulletSoftBody/BulletVariationalDeformable/btVariationalElasticity.h"
#include "BulletSoftBody/BulletVariationalDeformable/btVariationalHelpers.h"
#include "BulletSoftBody/BulletVariationalDeformable/btVariationalHybrid.h"
#include "BulletSoftBody/BulletVariationalDeformable/btVariationalIPC.h"
#include "BulletSoftBody/BulletVariationalDeformable/btVariationalProjectedNewton.h"
#include "BulletSoftBody/BulletVariationalDeformable/btVariationalProjectiveDynamics.h"
#include "BulletSoftBody/BulletVariationalDeformable/btVariationalQuasiNewton.h"
#include "BulletSoftBody/BulletVariationalDeformable/btVariationalSparse.h"
#include "BulletSoftBody/BulletVariationalDeformable/btVariationalVertexBlockDescent.h"

namespace
{
const double kDt = 1.0 / 60.0;

// small deterministic generator, so that the tests do not depend on the C library's rand()
struct Rng
{
	unsigned int m_state;
	explicit Rng(unsigned int seed) : m_state(seed) {}
	unsigned int bits()
	{
		m_state = m_state * 1664525u + 1013904223u;
		return m_state >> 8;
	}
	double uniform() { return bits() / 8388608.0 - 1.0; }  // in [-1, 1)
	int index(int n) { return int(bits() % (unsigned int)n); }
	btVec3d vec()
	{
		const double a = uniform(), b = uniform(), c = uniform();
		return btVec3d(a, b, c);
	}
	btMat3d mat()
	{
		btMat3d M;
		for (int r = 0; r < 3; r++)
			for (int c = 0; c < 3; c++) M.m[r][c] = uniform();
		return M;
	}
};

btSoftBodyWorldInfo& worldInfo()
{
	static btSoftBodyWorldInfo info;
	static bool initialized = false;
	if (!initialized)
	{
		info.m_sparsesdf.Initialize();
		initialized = true;
	}
	return info;
}

// a static ground plane y = 0 for systems that are stepped without a Bullet world
struct Ground
{
	btStaticPlaneShape m_shape;
	btCollisionObject m_object;
	Ground() : m_shape(btVector3(0, 1, 0), 0)
	{
		m_object.setCollisionShape(&m_shape);
		m_object.setFriction(1.0);
	}
	void addTo(btVariationalSystem& sys)
	{
		btTransform identity;
		identity.setIdentity();
		const int c = sys.addCollider(&m_object, &m_shape, identity, &worldInfo());
		for (int b = 0; b < sys.m_bodies.size(); b++) sys.m_bodies[b].m_colliders.push_back(c);
	}
};

void step(btVariationalSystem& sys, btVariationalIntegrator& integrator, double dt = kDt)
{
	sys.beginStep(dt);
	integrator.solve(sys);
	sys.endStep();
}

// restarts the system from its rest state and takes a few steps; returns the final positions
void runFromRest(btVariationalSystem& sys, btVariationalIntegrator& integrator, int steps, btVec3dArray& out)
{
	sys.readFromSoftBodies();
	for (int i = 0; i < sys.numNodes(); i++)
	{
		sys.m_x[i] = sys.m_fixedTarget[i];
		sys.m_v[i] = sys.m_vPrev[i] = btVec3d(0, 0, 0);
	}
	for (int s = 0; s < steps; s++) step(sys, integrator);
	out.resize(sys.numNodes());
	for (int i = 0; i < sys.numNodes(); i++) out[i] = sys.m_x[i];
}

double maxDifference(const btVec3dArray& a, const btVec3dArray& b)
{
	double m = 0;
	for (int i = 0; i < a.size(); i++) m = btMax(m, (a[i] - b[i]).maxAbs());
	return m;
}

double frobenius(const btMat3d& M) { return sqrt(M.frobeniusNorm2()); }

// bar along x, clamped at its left end
btSoftBody* createClampedBar(int nx = 12)
{
	btSoftBody* psb = btVariationalHelpers::createTetGrid(worldInfo(), btVector3(0, 1, 0), btVector3(0.6, 0.1, 0.1), nx, 2, 2, 1000);
	for (int i = 0; i < psb->m_nodes.size(); i++)
		if (psb->m_nodes[i].m_x.x() < btScalar(-0.599)) psb->setMass(i, 0);
	return psb;
}

btVec3d centreOfMass(const btVariationalSystem& s, int first, int count)
{
	btVec3d c(0, 0, 0);
	double m = 0;
	for (int i = first; i < first + count; i++)
	{
		c += s.m_x[i] * s.m_mass[i];
		m += s.m_mass[i];
	}
	return c / m;
}

btVec3d linearMomentum(const btVariationalSystem& s)
{
	btVec3d p(0, 0, 0);
	for (int i = 0; i < s.numNodes(); i++) p += s.m_v[i] * s.m_mass[i];
	return p;
}

btVec3d angularMomentum(const btVariationalSystem& s)
{
	const btVec3d c = centreOfMass(s, 0, s.numNodes());
	btVec3d L(0, 0, 0);
	for (int i = 0; i < s.numNodes(); i++) L += (s.m_x[i] - c).cross(s.m_v[i]) * s.m_mass[i];
	return L;
}

int countInvertedTets(const btVariationalSystem& sys)
{
	int inverted = 0;
	for (int t = 0; t < sys.m_tets.size(); t++)
		if (sys.deformationGradient(sys.m_tets[t], &sys.m_x[0]).det() <= 0) inverted++;
	return inverted;
}

double maxAttachmentError(const btVariationalSystem& sys)
{
	double err = 0;
	for (int a = 0; a < sys.m_attachments.size(); a++)
	{
		const btVariationalAttachment& at = sys.m_attachments[a];
		err = btMax(err, (sys.m_x[at.m_nodeA] - (at.m_nodeB >= 0 ? sys.m_x[at.m_nodeB] : at.m_target)).length());
	}
	return err;
}

// Compares gradient() with central differences of energy(), the unfiltered Hessian with central
// differences of gradient(), and the per-vertex force with the gradient.
struct DerivativeErrors
{
	double m_gradient, m_hessian, m_vertexForce;
	double m_compressed;  // block sparse form of the Hessian against its stencil form
	int m_indefiniteBlocks;
};

DerivativeErrors checkDerivatives(btVariationalSystem& sys, double perturbation, Rng& rng)
{
	DerivativeErrors e;
	const int n = sys.numNodes();
	btVec3dArray x, g, gp, gm, v, Hv, xp, xm;
	x.resize(n);
	g.resize(n);
	gp.resize(n);
	gm.resize(n);
	v.resize(n);
	Hv.resize(n);
	xp.resize(n);
	xm.resize(n);
	for (int i = 0; i < n; i++) x[i] = sys.m_x[i] + (sys.m_fixed[i] ? btVec3d(0, 0, 0) : rng.vec() * perturbation);
	sys.gradient(&x[0], &g[0]);
	const double h = 1e-6;
	double gscale = 1e-12;
	for (int i = 0; i < n; i++) gscale = btMax(gscale, g[i].maxAbs());
	e.m_gradient = 0;
	for (int s = 0; s < 60; s++)
	{
		const int i = rng.index(n), c = rng.index(3);
		if (sys.m_fixed[i]) continue;
		const double keep = x[i][c];
		x[i][c] = keep + h;
		const double Ep = sys.energy(&x[0]);
		x[i][c] = keep - h;
		const double Em = sys.energy(&x[0]);
		x[i][c] = keep;
		e.m_gradient = btMax(e.m_gradient, fabs((Ep - Em) / (2 * h) - g[i][c]) / gscale);
	}
	btVariationalHessian H;
	sys.hessian(&x[0], H, BT_VARIATIONAL_FILTER_NONE);
	for (int i = 0; i < n; i++) v[i] = sys.m_fixed[i] ? btVec3d(0, 0, 0) : rng.vec();
	H.multiply(&v[0], &Hv[0], &sys.m_fixed[0]);
	for (int i = 0; i < n; i++)
	{
		xp[i] = x[i] + v[i] * h;
		xm[i] = x[i] - v[i] * h;
	}
	sys.gradient(&xp[0], &gp[0]);
	sys.gradient(&xm[0], &gm[0]);
	double hscale = 1e-12;
	e.m_hessian = 0;
	for (int i = 0; i < n; i++)
		if (!sys.m_fixed[i]) hscale = btMax(hscale, Hv[i].maxAbs());
	for (int i = 0; i < n; i++)
		if (!sys.m_fixed[i]) e.m_hessian = btMax(e.m_hessian, ((gp[i] - gm[i]) * (1.0 / (2 * h)) - Hv[i]).maxAbs() / hscale);
	// the matrix the conjugate gradient solver works with: same product, fixed nodes included
	btVec3dArray w, Hw, Hc;
	w.resize(n);
	Hw.resize(n);
	Hc.resize(n);
	for (int i = 0; i < n; i++) w[i] = rng.vec();
	H.multiply(&w[0], &Hw[0], &sys.m_fixed[0]);
	H.compress(&sys.m_fixed[0]);
	H.multiplyCompressed(&w[0], &Hc[0]);
	e.m_compressed = 0;
	for (int i = 0; i < n; i++) e.m_compressed = btMax(e.m_compressed, (Hw[i] - Hc[i]).maxAbs() / hscale);
	e.m_vertexForce = 0;
	e.m_indefiniteBlocks = 0;
	for (int i = 0; i < n; i++)
	{
		if (sys.m_fixed[i]) continue;
		btVec3d f, tmp;
		btMat3d B;
		sys.vertexForceAndHessian(i, &x[0], f, B);
		e.m_vertexForce = btMax(e.m_vertexForce, (f + g[i]).maxAbs() / gscale);
		if (!btSolveSpd3d(B, f, tmp)) e.m_indefiniteBlocks++;
	}
	return e;
}

// independent reference distances (Ericson, Real-Time Collision Detection)
double referencePointTriangle(const btVec3d& p, const btVec3d& a, const btVec3d& b, const btVec3d& c)
{
	const btVec3d ab = b - a, ac = c - a, ap = p - a;
	const double d1 = ab.dot(ap), d2 = ac.dot(ap);
	if (d1 <= 0 && d2 <= 0) return (p - a).length2();
	const btVec3d bp = p - b;
	const double d3 = ab.dot(bp), d4 = ac.dot(bp);
	if (d3 >= 0 && d4 <= d3) return (p - b).length2();
	const double vc = d1 * d4 - d3 * d2;
	if (vc <= 0 && d1 >= 0 && d3 <= 0) return (p - (a + ab * (d1 / (d1 - d3)))).length2();
	const btVec3d cp = p - c;
	const double d5 = ab.dot(cp), d6 = ac.dot(cp);
	if (d6 >= 0 && d5 <= d6) return (p - c).length2();
	const double vb = d5 * d2 - d1 * d6;
	if (vb <= 0 && d2 >= 0 && d6 <= 0) return (p - (a + ac * (d2 / (d2 - d6)))).length2();
	const double va = d3 * d6 - d5 * d4;
	if (va <= 0 && (d4 - d3) >= 0 && (d5 - d6) >= 0) return (p - (b + (c - b) * ((d4 - d3) / ((d4 - d3) + (d5 - d6))))).length2();
	const double denom = 1.0 / (va + vb + vc);
	return (p - (a + ab * (vb * denom) + ac * (vc * denom))).length2();
}

double clamp01(double v) { return v < 0 ? 0 : (v > 1 ? 1 : v); }

double referenceSegmentSegment(const btVec3d& p1, const btVec3d& q1, const btVec3d& p2, const btVec3d& q2)
{
	const btVec3d d1 = q1 - p1, d2 = q2 - p2, r = p1 - p2;
	const double a = d1.dot(d1), e = d2.dot(d2), f = d2.dot(r), c = d1.dot(r), b = d1.dot(d2);
	const double denom = a * e - b * b;
	double s = denom > 1e-14 * a * e ? clamp01((b * f - c * e) / denom) : 0.0;
	double t = (b * s + f) / e;
	if (t < 0)
	{
		t = 0;
		s = clamp01(-c / a);
	}
	else if (t > 1)
	{
		t = 1;
		s = clamp01((b - c) / a);
	}
	return ((p1 + d1 * s) - (p2 + d2 * t)).length2();
}

// a light block pinned to the world at its top corners that carries a block 200 times heavier
// on four node-to-node attachments
struct HangingBlocks
{
	btSoftBody* m_light;
	btSoftBody* m_heavy;
	btVariationalStableNeoHookean m_material;
	btVariationalSystem m_sys;

	explicit HangingBlocks(bool hard) : m_material(1e6, 0.4)
	{
		m_light = btVariationalHelpers::createTetGrid(worldInfo(), btVector3(0, 2.0, 0), btVector3(0.1, 0.1, 0.1), 2, 2, 2, 100);
		m_heavy = btVariationalHelpers::createTetGrid(worldInfo(), btVector3(0, 1.7, 0), btVector3(0.2, 0.2, 0.2), 2, 2, 2, 2500);
		int top[4], bottom[4], heavyTop[4], numTop = 0, numBottom = 0, numHeavy = 0;
		for (int i = 0; i < m_light->m_nodes.size(); i++)
		{
			const btVector3 x = m_light->m_nodes[i].m_x;
			if (btFabs(btFabs(x.x()) - btScalar(0.1)) < btScalar(1e-4) && btFabs(btFabs(x.z()) - btScalar(0.1)) < btScalar(1e-4))
			{
				if (x.y() > btScalar(2.09)) top[numTop++] = i;
				if (x.y() < btScalar(1.91)) bottom[numBottom++] = i;
			}
		}
		for (int i = 0; i < m_heavy->m_nodes.size(); i++)
		{
			const btVector3 x = m_heavy->m_nodes[i].m_x;
			if (btFabs(btFabs(x.x()) - btScalar(0.2)) < btScalar(1e-4) && btFabs(btFabs(x.z()) - btScalar(0.2)) < btScalar(1e-4) && x.y() > btScalar(1.89)) heavyTop[numHeavy++] = i;
		}
		// move the heavy block's top corners onto the light block's bottom corners: the joints start closed
		int match[4];
		for (int c = 0; c < 4; c++)
		{
			int best = 0;
			btScalar bestDistance = btScalar(1e9);
			for (int d = 0; d < 4; d++)
			{
				const btScalar distance = (m_light->m_nodes[bottom[c]].m_x - m_heavy->m_nodes[heavyTop[d]].m_x).length();
				if (distance < bestDistance)
				{
					bestDistance = distance;
					best = d;
				}
			}
			match[c] = heavyTop[best];
			m_heavy->m_nodes[match[c]].m_x = m_light->m_nodes[bottom[c]].m_x;
		}
		m_sys.addBody(m_light, &m_material, 0, 0);
		m_sys.addBody(m_heavy, &m_material, 0, 0);
		m_sys.buildTopology();
		m_sys.m_gravity = btVec3d(0, -9.8, 0);
		m_sys.m_damping = 3.0;
		const double stiffness = 1e5;
		for (int c = 0; c < 4; c++)
		{
			btVariationalAttachment pin;
			pin.m_nodeA = top[c];
			pin.m_target = pin.m_prevTarget = m_sys.m_x[top[c]];
			pin.m_stiffness = pin.m_penalty = stiffness;
			pin.m_hard = hard;
			m_sys.addAttachment(pin);
			btVariationalAttachment joint;
			joint.m_nodeA = bottom[c];
			joint.m_nodeB = m_sys.m_bodies[1].m_firstNode + match[c];
			joint.m_stiffness = joint.m_penalty = stiffness;
			joint.m_hard = hard;
			m_sys.addAttachment(joint);
		}
	}
	~HangingBlocks()
	{
		delete m_light;
		delete m_heavy;
	}
};
}  // namespace

//
// Small dense linear algebra
//

GTEST_TEST(VariationalMath, SvdGivesProperRotations)
{
	Rng rng(1);
	const btMat3d I = btMat3d::identity();
	for (int it = 0; it < 3000; it++)
	{
		btMat3d F = rng.mat();
		if (it % 4 == 1) F.setCol(2, F.col(2) * 1e-7);         // nearly singular
		if (it % 4 == 2) F = I + F * 1e-3;                     // nearly repeated singular values
		if (it % 4 == 3) F.setCol(0, F.col(1) * 2.0);          // rank deficient
		btMat3d U, V;
		btVec3d S;
		btSvd3d(F, U, S, V);
		const double scale = 1.0 + frobenius(F);
		EXPECT_LT(frobenius(U * btMat3d::diag(S[0], S[1], S[2]) * V.transpose() - F), 1e-9 * scale);
		EXPECT_LT(frobenius(U.transpose() * U - I), 1e-9);
		EXPECT_LT(frobenius(V.transpose() * V - I), 1e-9);
		EXPECT_NEAR(U.det(), 1.0, 1e-9);
		EXPECT_NEAR(V.det(), 1.0, 1e-9);
		EXPECT_GE(S[0], S[1]);
		EXPECT_GE(S[1], fabs(S[2]));
		if (fabs(F.det()) > 1e-9) EXPECT_GT(S[2] * F.det(), 0.0);  // the sign of the determinant sits in S[2]
	}
}

GTEST_TEST(VariationalMath, PolarRotationIsClosestRotation)
{
	Rng rng(2);
	const btMat3d I = btMat3d::identity();
	for (int it = 0; it < 500; it++)
	{
		const btMat3d F = I + rng.mat() * 0.5;
		const btMat3d R = btPolarRotation3d(F);
		EXPECT_LT(frobenius(R.transpose() * R - I), 1e-9);
		EXPECT_NEAR(R.det(), 1.0, 1e-9);
		const btMat3d S = R.transpose() * F;  // the stretch must be symmetric
		EXPECT_LT(frobenius(S - S.transpose()), 1e-9);
	}
}

GTEST_TEST(VariationalMath, EigenvalueFiltersMakeMatricesPositiveSemiDefinite)
{
	Rng rng(3);
	for (int it = 0; it < 100; it++)
	{
		btMat9d A;
		for (int r = 0; r < 9; r++)
			for (int c = r; c < 9; c++)
			{
				const double value = rng.uniform();
				A(r, c) = value;
				A(c, r) = value;
			}
		double evals[9], filtered[9];
		btMat9d V, W;
		btSymmetricEigenJacobi(A, evals, V);
		// A V = V diag(evals)
		for (int k = 0; k < 9; k++)
			for (int r = 0; r < 9; r++)
			{
				double av = 0;
				for (int c = 0; c < 9; c++) av += A(r, c) * V(c, k);
				EXPECT_NEAR(av, evals[k] * V(r, k), 1e-9);
			}
		for (int f = 0; f < 2; f++)
		{
			const btVariationalHessianFilter filter = f == 0 ? BT_VARIATIONAL_FILTER_CLAMP : BT_VARIATIONAL_FILTER_ABS;
			btMat9d B = A;
			btFilterSymmetric(B, filter);
			btSymmetricEigenJacobi(B, filtered, W);
			double sumExpected = 0, sumFiltered = 0;
			for (int k = 0; k < 9; k++)
			{
				EXPECT_GT(filtered[k], -1e-9);
				sumExpected += f == 0 ? btMax(evals[k], 0.0) : fabs(evals[k]);
				sumFiltered += filtered[k];
			}
			EXPECT_NEAR(sumExpected, sumFiltered, 1e-8);
		}
	}
}

//
// Materials
//

namespace
{
void checkMaterial(const btVariationalMaterial& material, unsigned int seed)
{
	Rng rng(seed);
	const btMat3d I = btMat3d::identity();
	double worstStress = 0, worstHessian = 0, worstOrthogonality = 0, worstReconstruction = 0, mostNegative = 0;
	for (int it = 0; it < 400; it++)
	{
		btMat3d F = rng.mat();
		const int kind = it % 4;
		if (kind == 0) F = I + F * 0.3;  // moderate deformation
		if (kind == 1) F = F * 1.5;      // wild, often inverted
		if (kind == 2)
		{
			F = I + F * 0.3;  // inverted
			F.setCol(2, F.col(2) * -1.0);
		}
		if (kind == 3)
		{
			F = I + F * 0.2;  // squashed almost flat
			F.setCol(1, F.col(1) * 0.05);
		}
		const btMat3d P = material.firstPiola(F);
		const double h = 1e-6, pscale = frobenius(P) + 1.0;
		for (int r = 0; r < 3; r++)
			for (int c = 0; c < 3; c++)
			{
				btMat3d Fp = F, Fm = F;
				Fp.m[r][c] += h;
				Fm.m[r][c] -= h;
				worstStress = btMax(worstStress, fabs((material.energyDensity(Fp) - material.energyDensity(Fm)) / (2 * h) - P.m[r][c]) / pscale);
			}
		btMat9d H;
		material.hessian(F, H, BT_VARIATIONAL_FILTER_NONE);
		double hscale = 0;
		for (int i = 0; i < 81; i++) hscale += H.m[i] * H.m[i];
		hscale = sqrt(hscale) + 1.0;
		for (int c = 0; c < 3; c++)
			for (int r = 0; r < 3; r++)
			{
				btMat3d Fp = F, Fm = F;
				Fp.m[r][c] += h;
				Fm.m[r][c] -= h;
				const btMat3d dP = (material.firstPiola(Fp) - material.firstPiola(Fm)) * (1.0 / (2 * h));
				for (int c2 = 0; c2 < 3; c2++)
					for (int r2 = 0; r2 < 3; r2++) worstHessian = btMax(worstHessian, fabs(dP.m[r2][c2] - H(r2 + 3 * c2, r + 3 * c)) / hscale);
			}
		// the analytic eigen-system is orthonormal and reproduces the Hessian
		btVariationalEigenSystem es;
		material.eigenSystem(F, es);
		btMat9d R;
		R.setZero();
		for (int i = 0; i < 9; i++)
		{
			for (int j = 0; j < 9; j++) worstOrthogonality = btMax(worstOrthogonality, fabs(es.m_Q[i].ddot(es.m_Q[j]) - (i == j ? 1.0 : 0.0)));
			btVec9d q;
			for (int c = 0; c < 3; c++)
				for (int r = 0; r < 3; r++) q.v[r + 3 * c] = es.m_Q[i].m[r][c];
			R.addScaledOuter(es.m_lambda[i], q, q);
		}
		for (int i = 0; i < 81; i++) worstReconstruction = btMax(worstReconstruction, fabs(R.m[i] - H.m[i]) / hscale);
		// both filters give positive semi-definite Hessians
		for (int f = 0; f < 2; f++)
		{
			btMat9d Hf, V;
			double evals[9];
			material.hessian(F, Hf, f == 0 ? BT_VARIATIONAL_FILTER_CLAMP : BT_VARIATIONAL_FILTER_ABS);
			btSymmetricEigenJacobi(Hf, evals, V);
			for (int i = 0; i < 9; i++) mostNegative = btMin(mostNegative, evals[i] / hscale);
		}
	}
	EXPECT_LT(worstStress, 1e-6);
	EXPECT_LT(worstHessian, 1e-6);
	EXPECT_LT(worstOrthogonality, 1e-9);
	EXPECT_LT(worstReconstruction, 1e-9);
	EXPECT_GT(mostNegative, -1e-9);
	// the rest state is stress free
	EXPECT_LT(frobenius(material.firstPiola(I)), 1e-9 * material.laplacianStiffness());
}
}  // namespace

GTEST_TEST(VariationalElasticity, StableNeoHookeanDerivatives)
{
	checkMaterial(btVariationalStableNeoHookean(1e4, 0.45, true), 11);
	checkMaterial(btVariationalStableNeoHookean(1e4, 0.45, false), 12);
	checkMaterial(btVariationalStableNeoHookean(1e4, 0.495, true), 13);  // nearly incompressible
}

GTEST_TEST(VariationalElasticity, ArapDerivativesAndProjection)
{
	btVariationalARAP arap(1e4, 0.3);
	checkMaterial(arap, 14);
	// the projective dynamics form of the energy: P = stiffness (F - target)
	Rng rng(15);
	for (int it = 0; it < 100; it++)
	{
		const btMat3d F = btMat3d::identity() + rng.mat() * 0.4;
		btMat3d target;
		double stiffness;
		arap.projectiveDynamicsTarget(F, target, stiffness);
		EXPECT_LT(frobenius(arap.firstPiola(F) - (F - target) * stiffness), 1e-8 * stiffness);
	}
}

GTEST_TEST(VariationalElasticity, RestLinearizationIsLinearElasticity)
{
	// dP = 2 mu eps + lambda tr(eps) I for a small symmetric strain eps, with the Lame parameters of (E, nu)
	const double E = 1e4, nu = 0.45, mu = E / (2 * (1 + nu)), lambda = E * nu / ((1 + nu) * (1 - 2 * nu));
	Rng rng(16);
	btMat3d eps;
	for (int r = 0; r < 3; r++)
		for (int c = r; c < 3; c++)
		{
			const double value = rng.uniform();
			eps.m[r][c] = value;
			eps.m[c][r] = value;
		}
	const btMat3d I = btMat3d::identity();
	const btMat3d linear = eps * (2 * mu) + I * (lambda * eps.trace());
	for (int logTerm = 0; logTerm < 2; logTerm++)
	{
		btVariationalStableNeoHookean material(E, nu, logTerm == 1);
		const double t = 1e-6;
		const btMat3d dP = (material.firstPiola(I + eps * t) - material.firstPiola(I - eps * t)) * (1.0 / (2 * t));
		EXPECT_LT(frobenius(dP - linear) / frobenius(linear), 1e-6);
	}
}

//
// Mesh generators
//

namespace
{
double tetVolume(const btSoftBody* psb, double* smallest)
{
	double volume = 0;
	*smallest = 1e300;
	for (int t = 0; t < psb->m_tetras.size(); t++)
	{
		const btSoftBody::Tetra& tet = psb->m_tetras[t];
		const btVector3 a = tet.m_n[1]->m_x - tet.m_n[0]->m_x, b = tet.m_n[2]->m_x - tet.m_n[0]->m_x, c = tet.m_n[3]->m_x - tet.m_n[0]->m_x;
		const double v = double(a.cross(b).dot(c)) / 6.0;
		volume += v;
		*smallest = btMin(*smallest, v);
	}
	return volume;
}

// volume enclosed by the quads; positive if they are oriented outwards
double quadVolume(const btSoftBody* psb, const btAlignedObjectArray<int>& quads)
{
	double volume = 0;
	for (int q = 0; q + 3 < quads.size(); q += 4)
	{
		const btVector3 p0 = psb->m_nodes[quads[q]].m_x, p1 = psb->m_nodes[quads[q + 1]].m_x, p2 = psb->m_nodes[quads[q + 2]].m_x, p3 = psb->m_nodes[quads[q + 3]].m_x;
		volume += double(p0.dot(p1.cross(p2)) + p0.dot(p2.cross(p3))) / 6.0;
	}
	return volume;
}
}  // namespace

GTEST_TEST(VariationalHelpers, EllipsoidAndRoundedBox)
{
	btAlignedObjectArray<int> quads;
	const btVector3 radii(0.3, 0.1, 0.2);
	btSoftBody* ellipsoid = btVariationalHelpers::createTetEllipsoid(worldInfo(), btVector3(1, 2, 3), radii, 8, 4, 6, 1000);
	btVariationalHelpers::gridSurfaceQuads(8, 4, 6, quads);
	double smallest;
	const double volume = tetVolume(ellipsoid, &smallest), exact = 4.0 / 3.0 * SIMD_PI * 0.3 * 0.1 * 0.2;
	EXPECT_GT(smallest, 0.0);  // no inverted or flat element
	EXPECT_LT(volume, exact);  // the polyhedron is inscribed
	EXPECT_GT(volume, 0.93 * exact);
	EXPECT_EQ(2 * (8 * 4 + 4 * 6 + 6 * 8), quads.size() / 4);
	EXPECT_EQ(2 * (quads.size() / 4), ellipsoid->m_faces.size());  // Bullet's collision triangles split the same quads
	// closed and oriented outwards (the quads are not flat, and the tetrahedra may split one along
	// its other diagonal, hence the tolerance)
	EXPECT_NEAR(quadVolume(ellipsoid, quads), volume, 1e-3 * volume);
	double mass = 0;
	for (int i = 0; i < ellipsoid->m_nodes.size(); i++) mass += 1.0 / ellipsoid->m_nodes[i].m_im;
	EXPECT_NEAR(mass, 1000.0 * volume, 1e-4 * mass);
	// every surface node lies on the ellipsoid
	for (int q = 0; q < quads.size(); q++)
	{
		const btVector3 d = ellipsoid->m_nodes[quads[q]].m_x - btVector3(1, 2, 3);
		EXPECT_NEAR(btSqrt(d.x() * d.x() / (radii.x() * radii.x()) + d.y() * d.y() / (radii.y() * radii.y()) + d.z() * d.z() / (radii.z() * radii.z())), 1.0, 1e-4);
	}
	delete ellipsoid;

	// roundness 0 is the plain block
	btSoftBody* block = btVariationalHelpers::createTetEllipsoid(worldInfo(), btVector3(0, 0, 0), btVector3(0.2, 0.3, 0.1), 3, 4, 2, 1000, 0);
	btSoftBody* grid = btVariationalHelpers::createTetGrid(worldInfo(), btVector3(0, 0, 0), btVector3(0.2, 0.3, 0.1), 3, 4, 2, 1000);
	ASSERT_EQ(grid->m_nodes.size(), block->m_nodes.size());
	for (int i = 0; i < grid->m_nodes.size(); i++) EXPECT_LT((grid->m_nodes[i].m_x - block->m_nodes[i].m_x).length(), btScalar(1e-6));
	EXPECT_NEAR(tetVolume(block, &smallest), 0.4 * 0.6 * 0.2, 1e-6);
	delete block;
	delete grid;
}

GTEST_TEST(VariationalHelpers, Torus)
{
	btAlignedObjectArray<int> quads;
	btSoftBody* torus = btVariationalHelpers::createTetTorus(worldInfo(), btVector3(0, 1, 0), 0.15, 0.05, 24, 3, 1000);
	btVariationalHelpers::torusSurfaceQuads(24, 3, quads);
	double smallest;
	const double volume = tetVolume(torus, &smallest), exact = 2.0 * SIMD_PI * SIMD_PI * 0.15 * 0.05 * 0.05;
	EXPECT_GT(smallest, 0.0);
	EXPECT_LT(volume, exact);
	EXPECT_GT(volume, 0.92 * exact);
	EXPECT_EQ(24 * 4 * 3, quads.size() / 4);
	EXPECT_EQ(2 * (quads.size() / 4), torus->m_faces.size());
	EXPECT_NEAR(quadVolume(torus, quads), volume, 1e-3 * volume);
	// the mesh closes on itself around the ring: left alone it stays exactly as it is
	btVariationalStableNeoHookean material(1e5, 0.4);
	btVariationalSystem sys;
	sys.addBody(torus, &material, 0, 0);
	sys.buildTopology();
	sys.m_gravity = btVec3d(0, 0, 0);
	btVec3dArray rest;
	rest.resize(sys.numNodes());
	for (int i = 0; i < sys.numNodes(); i++) rest[i] = sys.m_x[i];
	btVariationalProjectedNewton newton;
	for (int s = 0; s < 3; s++) step(sys, newton);
	double moved = 0;
	for (int i = 0; i < sys.numNodes(); i++) moved = btMax(moved, (sys.m_x[i] - rest[i]).length());
	EXPECT_LT(moved, 1e-12);
	EXPECT_LT(sys.elasticEnergy(&sys.m_x[0]), 1e-20);
	delete torus;
}

//
// Sparse Cholesky
//

namespace
{
struct Triplets
{
	btAlignedObjectArray<int> m_rows, m_cols;
	btAlignedObjectArray<double> m_values;
	void add(int i, int j, double value)
	{
		m_rows.push_back(i);
		m_cols.push_back(j);
		m_values.push_back(value);
	}
};

// graph Laplacian of an nx * ny * nz grid with random positive weights, plus a diagonal shift
int addGridLaplacian(Triplets& t, int offset, int nx, int ny, int nz, double shift, Rng& rng)
{
	for (int k = 0; k < nz; k++)
		for (int j = 0; j < ny; j++)
			for (int i = 0; i < nx; i++)
			{
				const int a = offset + i + nx * (j + ny * k);
				t.add(a, a, shift);
				for (int d = 0; d < 3; d++)
				{
					const int i2 = i + (d == 0), j2 = j + (d == 1), k2 = k + (d == 2);
					if (i2 >= nx || j2 >= ny || k2 >= nz) continue;
					const int b = offset + i2 + nx * (j2 + ny * k2);
					const double w = 1.0 + 0.5 * rng.uniform();
					t.add(a, a, w);
					t.add(b, b, w);
					t.add(a, b, -w);
					t.add(b, a, -w);
				}
			}
	return nx * ny * nz;
}

double choleskyError(int n, const Triplets& t, Rng& rng)
{
	btVariationalSparseMatrix A;
	A.setFromTriplets(n, t.m_rows, t.m_cols, t.m_values);
	btVariationalSparseCholesky cholesky;
	if (!cholesky.factorize(A)) return 1e300;
	btVec3dArray reference, b, x;
	reference.resize(n);
	b.resize(n);
	x.resize(n);
	for (int i = 0; i < n; i++) reference[i] = rng.vec();
	A.multiply(&reference[0], &b[0]);
	cholesky.solve(&b[0], &x[0]);
	double error = 0, norm = 0;
	for (int i = 0; i < n; i++)
	{
		error += (x[i] - reference[i]).length2();
		norm += reference[i].length2();
	}
	return sqrt(error / norm);
}
}  // namespace

GTEST_TEST(VariationalSparse, CholeskySolvesLaplacians)
{
	Rng rng(21);
	{
		Triplets t;
		const int n = addGridLaplacian(t, 0, 12, 12, 12, 0.01, rng);
		EXPECT_LT(choleskyError(n, t, rng), 1e-8);
	}
	{
		Triplets t;  // long thin beam
		const int n = addGridLaplacian(t, 0, 150, 3, 3, 0.01, rng);
		EXPECT_LT(choleskyError(n, t, rng), 1e-8);
	}
	{
		Triplets t;  // disconnected pieces and rows that only have a diagonal entry
		int n = 0;
		for (int piece = 0; piece < 12; piece++) n += addGridLaplacian(t, n, 3 + piece % 4, 3 + piece % 4, 3 + piece % 4, 0.05, rng);
		for (int q = 0; q < 100; q++)
		{
			t.add(n, n, 2.0 + rng.uniform());
			n++;
		}
		EXPECT_LT(choleskyError(n, t, rng), 1e-8);
	}
}

GTEST_TEST(VariationalSparse, CholeskyRejectsIndefiniteMatrix)
{
	Triplets t;
	t.add(0, 0, 1);
	t.add(1, 1, -1);
	t.add(0, 1, 0.1);
	t.add(1, 0, 0.1);
	btVariationalSparseMatrix A;
	A.setFromTriplets(2, t.m_rows, t.m_cols, t.m_values);
	btVariationalSparseCholesky cholesky;
	EXPECT_FALSE(cholesky.factorize(A));
}

//
// The incremental potential
//

GTEST_TEST(VariationalSystem, DerivativesOfTheIncrementalPotential)
{
	btVariationalStableNeoHookean snh(2e5, 0.45), snhNoLog(2e5, 0.45, false);
	btVariationalARAP arap(2e5, 0.3);
	const btVariationalMaterial* materials[3] = {&snh, &snhNoLog, &arap};
	Rng rng(31);
	// variants 0-2: the three materials alone; 3: attachments, a ground plane and friction on top
	for (int variant = 0; variant < 4; variant++)
	{
		btSoftBody* psb = btVariationalHelpers::createTetGrid(worldInfo(), btVector3(0, 1, 0), btVector3(0.5, 0.15, 0.15), 5, 2, 2, 1000);
		for (int i = 0; i < psb->m_nodes.size(); i++)
			if (psb->m_nodes[i].m_x.x() < btScalar(-0.49)) psb->setMass(i, 0);
		btVariationalSystem sys;
		sys.addBody(psb, materials[variant % 3], 0, 0);
		sys.buildTopology();
		sys.m_gravity = btVec3d(0, -9.8, 0);
		btStaticPlaneShape plane(btVector3(0, 1, 0), 0);
		btCollisionObject object;
		object.setCollisionShape(&plane);
		object.setFriction(1.0);
		if (variant == 3)
		{
			sys.m_colliderMargin = 0.25;  // wide enough that many nodes are in contact
			sys.m_bodies[0].m_friction = 0.5;
			btTransform t;
			t.setIdentity();
			t.setOrigin(btVector3(0, 0.75, 0));
			t.setRotation(btQuaternion(btVector3(1, 0, 0), 0.1));
			sys.m_bodies[0].m_colliders.push_back(sys.addCollider(&object, &plane, t, &worldInfo()));
			btVariationalAttachment toWorld;
			toWorld.m_nodeA = 20;
			toWorld.m_target = btVec3d(0.2, 1.3, 0.1);
			toWorld.m_stiffness = 5e3;
			sys.addAttachment(toWorld);
			btVariationalAttachment nodeToNode;
			nodeToNode.m_nodeA = 30;
			nodeToNode.m_nodeB = 44;
			nodeToNode.m_stiffness = 2e3;
			nodeToNode.m_lambda = btVec3d(3, -2, 1);
			sys.addAttachment(nodeToNode);
		}
		for (int i = 0; i < sys.numNodes(); i++) sys.m_v[i] = rng.vec() * 0.5;
		sys.beginStep(kDt);
		// lagged friction data from a perturbed state, so that the friction terms are not trivial
		btVec3dArray lagged;
		lagged.resize(sys.numNodes());
		for (int i = 0; i < sys.numNodes(); i++) lagged[i] = sys.m_x[i] + rng.vec() * 0.01;
		sys.updateFrictionLag(&lagged[0]);
		for (int large = 0; large < 2; large++)
		{
			// the large perturbation inverts elements
			const DerivativeErrors e = checkDerivatives(sys, large ? 0.12 : 0.01, rng);
			EXPECT_LT(e.m_gradient, 1e-5) << "variant " << variant;
			EXPECT_LT(e.m_hessian, 1e-5) << "variant " << variant;
			EXPECT_LT(e.m_compressed, 1e-12) << "variant " << variant;
			EXPECT_LT(e.m_vertexForce, 1e-9) << "variant " << variant;
			EXPECT_EQ(0, e.m_indefiniteBlocks) << "variant " << variant;
		}
		delete psb;
	}
}

GTEST_TEST(VariationalSystem, CurvedColliderGradients)
{
	// Box, sphere and capsule colliders: the gradient is exact. (Their Hessian drops the curvature of
	// the collider surface on purpose, so it is not compared here.)
	btVariationalStableNeoHookean snh(2e5, 0.45);
	Rng rng(32);
	btSoftBody* psb = btVariationalHelpers::createTetGrid(worldInfo(), btVector3(0, 1, 0), btVector3(0.5, 0.15, 0.15), 5, 2, 2, 1000);
	btVariationalSystem sys;
	sys.addBody(psb, &snh, 0, 0);
	sys.buildTopology();
	sys.m_colliderMargin = 0.25;
	btBoxShape box(btVector3(0.3, 0.3, 0.3));
	btSphereShape sphere(0.25);
	btCapsuleShape capsule(0.2, 0.5);
	btCollisionObject boxObject, sphereObject, capsuleObject;
	boxObject.setCollisionShape(&box);
	sphereObject.setCollisionShape(&sphere);
	capsuleObject.setCollisionShape(&capsule);
	btTransform t;
	t.setIdentity();
	t.setOrigin(btVector3(0.3, 0.45, 0.1));
	t.setRotation(btQuaternion(btVector3(0, 0, 1), 0.3));
	sys.m_bodies[0].m_colliders.push_back(sys.addCollider(&boxObject, &box, t, &worldInfo()));
	t.setOrigin(btVector3(-0.2, 1.4, 0));
	sys.m_bodies[0].m_colliders.push_back(sys.addCollider(&sphereObject, &sphere, t, &worldInfo()));
	t.setOrigin(btVector3(0.1, 1.0, 0.45));
	sys.m_bodies[0].m_colliders.push_back(sys.addCollider(&capsuleObject, &capsule, t, &worldInfo()));
	sys.beginStep(kDt);
	const DerivativeErrors e = checkDerivatives(sys, 0.01, rng);
	EXPECT_LT(e.m_gradient, 1e-5);
	EXPECT_LT(e.m_compressed, 1e-12);
	EXPECT_LT(e.m_vertexForce, 1e-9);
	EXPECT_EQ(0, e.m_indefiniteBlocks);
	delete psb;
}

//
// Projected Newton and the Stable Neo-Hookean material
//

GTEST_TEST(VariationalProjectedNewton, HangingBarMatchesLinearElasticity)
{
	// A bar hanging under its own weight, solved quasi-statically (huge time step). Small strain theory:
	// elongation rho g L^2 / (2 E), lateral strain = -nu * axial strain.
	const double L = 1.0, E = 1e6, nu = 0.3, rho = 1000, g = 9.8;
	for (int logTerm = 0; logTerm < 2; logTerm++)
	{
		btSoftBody* psb = btVariationalHelpers::createTetGrid(worldInfo(), btVector3(0, -0.5, 0), btVector3(0.05, 0.5, 0.05), 2, 20, 2, rho);
		for (int i = 0; i < psb->m_nodes.size(); i++)
			if (psb->m_nodes[i].m_x.y() > btScalar(-1e-4)) psb->setMass(i, 0);
		btVariationalStableNeoHookean material(E, nu, logTerm == 1);
		btVariationalSystem sys;
		sys.addBody(psb, &material, 0, 0);
		sys.buildTopology();
		sys.m_gravity = btVec3d(0, -g, 0);
		btVariationalProjectedNewton newton;
		newton.m_tolerance = 1e-10;
		newton.m_linearTolerance = 1e-10;
		newton.m_maxIterations = 100;
		for (int s = 0; s < 3; s++) step(sys, newton, 100.0);
		EXPECT_TRUE(newton.m_stats.m_converged);
		double tip = 0, xmin = 1e9, xmax = -1e9;
		int count = 0;
		for (int i = 0; i < sys.numNodes(); i++)
		{
			const btScalar y0 = psb->m_nodes[i].m_q.y();  // rest position
			if (y0 < btScalar(-0.999))
			{
				tip += sys.m_x[i][1] + 1.0;
				count++;
			}
			if (btFabs(y0 + btScalar(0.5)) < btScalar(1e-3))
			{
				xmin = btMin(xmin, sys.m_x[i][0]);
				xmax = btMax(xmax, sys.m_x[i][0]);
			}
		}
		ASSERT_GT(count, 0);
		const double elongation = -tip / count, theory = rho * g * L * L / (2 * E);
		EXPECT_NEAR(elongation / theory, 1.0, 0.03);
		const double lateral = (xmax - xmin) / 0.1 - 1.0, axialAtMidHeight = rho * g * (L / 2) / E;
		EXPECT_NEAR(-lateral / axialAtMidHeight, nu, 0.02);
		delete psb;
	}
}

GTEST_TEST(VariationalProjectedNewton, FreeFlightConservesMomentum)
{
	btSoftBody* psb = btVariationalHelpers::createTetGrid(worldInfo(), btVector3(0, 5, 0), btVector3(0.3, 0.2, 0.1), 4, 3, 2, 1000);
	btVariationalStableNeoHookean material(1e5, 0.4);
	btVariationalSystem sys;
	sys.addBody(psb, &material, 0, 0);
	sys.buildTopology();
	sys.m_gravity = btVec3d(0, -9.8, 0);
	const btVec3d c0 = centreOfMass(sys, 0, sys.numNodes()), omega(1.0, 2.0, 3.0), v0(0.5, 1.0, -0.2);
	double mass = 0;
	for (int i = 0; i < sys.numNodes(); i++)
	{
		sys.m_v[i] = v0 + omega.cross(sys.m_x[i] - c0);
		mass += sys.m_mass[i];
	}
	const btVec3d L0 = angularMomentum(sys);
	btVariationalProjectedNewton newton;
	newton.m_tolerance = 1e-8;
	newton.m_linearTolerance = 1e-10;
	const int N = 60;
	for (int s = 0; s < N; s++) step(sys, newton);
	const btVec3d expectedMomentum = (v0 + sys.m_gravity * (N * kDt)) * mass;
	EXPECT_LT((linearMomentum(sys) - expectedMomentum).length() / expectedMomentum.length(), 1e-6);
	// implicit Euler: x_{k+1} = x_k + h v_{k+1}, so the centre of mass is c0 + h sum_k (v0 + g h k)
	const btVec3d expectedCentre = c0 + v0 * (N * kDt) + sys.m_gravity * (kDt * kDt * N * (N + 1) / 2.0);
	EXPECT_LT((centreOfMass(sys, 0, sys.numNodes()) - expectedCentre).length(), 1e-6);
	// implicit Euler dissipates: the spin may decay but must not grow, and the axis stays put
	const btVec3d L1 = angularMomentum(sys);
	EXPECT_LT(L1.length(), L0.length() * (1.0 + 1e-9));
	EXPECT_GT(L1.length(), L0.length() * 0.5);
	EXPECT_GT(L1.dot(L0) / (L1.length() * L0.length()), 0.999);
	EXPECT_NEAR(sys.totalVolume(&sys.m_x[0]) / (0.6 * 0.4 * 0.2), 1.0, 0.01);
	delete psb;
}

GTEST_TEST(VariationalProjectedNewton, RecoversFromDegenerateAndInvertedStates)
{
	// scenario 0: the cube is squashed into a plane, every element has zero volume.
	// scenario 1: every node starts at a random position, half of the elements are inverted.
	for (int scenario = 0; scenario < 2; scenario++)
		for (int f = 0; f < 2; f++)
		{
			Rng rng(41);
			btSoftBody* psb = btVariationalHelpers::createTetGrid(worldInfo(), btVector3(0, 0, 0), btVector3(0.5, 0.5, 0.5), 3, 3, 3, 1000);
			btVariationalStableNeoHookean material(1e5, 0.4);
			btVariationalSystem sys;
			sys.addBody(psb, &material, 0, 0);
			sys.buildTopology();
			sys.m_gravity = btVec3d(0, 0, 0);
			sys.m_damping = 1.0;
			sys.m_hessianFilter = f == 0 ? BT_VARIATIONAL_FILTER_CLAMP : BT_VARIATIONAL_FILTER_ABS;
			for (int i = 0; i < sys.numNodes(); i++)
			{
				if (scenario == 0)
					sys.m_x[i][1] = 0.0;
				else
					sys.m_x[i] = rng.vec() * 0.5;
			}
			EXPECT_GT(countInvertedTets(sys), sys.m_tets.size() / 4);
			const double initialEnergy = sys.elasticEnergy(&sys.m_x[0]);
			btVariationalProjectedNewton newton;
			newton.m_tolerance = 1e-4;
			for (int s = 0; s < 240; s++)
			{
				step(sys, newton);
				EXPECT_TRUE(newton.m_stats.m_converged);
			}
			EXPECT_EQ(0, countInvertedTets(sys)) << "scenario " << scenario << ", filter " << f;
			const double energyRatio = sys.elasticEnergy(&sys.m_x[0]) / initialEnergy;
			if (scenario == 0)
			{
				// the flat cube pops back into its rest shape
				EXPECT_LT(energyRatio, 1e-4) << "filter " << f;
				EXPECT_NEAR(sys.totalVolume(&sys.m_x[0]), 1.0, 0.01) << "filter " << f;
			}
			else
			{
				// The scrambled cube untangles. (It usually ends in the rest shape as well, but where a
				// chaotic start ends up is not something to pin down in a test: a coarse mesh has
				// other, non-inverted equilibria.)
				EXPECT_LT(energyRatio, 0.05) << "filter " << f;
			}
			delete psb;
		}
}

GTEST_TEST(VariationalProjectedNewton, ContactAndFrictionOnAnIncline)
{
	// gravity is tilted instead of the plane, so "downhill" is +x
	const double theta = 0.3;  // tan(theta) = 0.309
	for (int variant = 0; variant < 3; variant++)
	{
		const double slope = variant == 0 ? 0.0 : theta, mu = variant == 2 ? 0.1 : 0.6;
		btSoftBody* psb = btVariationalHelpers::createTetGrid(worldInfo(), btVector3(0, 0.25, 0), btVector3(0.2, 0.2, 0.2), 2, 2, 2, 1000);
		btVariationalStableNeoHookean material(2e5, 0.4);
		btVariationalSystem sys;
		sys.addBody(psb, &material, 0, 0);
		sys.buildTopology();
		sys.m_bodies[0].m_friction = mu;
		sys.m_gravity = btVec3d(9.8 * sin(slope), -9.8 * cos(slope), 0);
		sys.m_damping = 2.0;
		Ground ground;
		ground.addTo(sys);
		btVariationalProjectedNewton newton;
		newton.m_tolerance = 1e-5;
		const int steps = 240;
		int unconverged = 0;
		for (int s = 0; s < steps; s++)
		{
			step(sys, newton);
			if (!newton.m_stats.m_converged) unconverged++;
		}
		EXPECT_EQ(0, unconverged);
		double lowest = 1e9, mass = 0;
		for (int i = 0; i < sys.numNodes(); i++)
		{
			lowest = btMin(lowest, sys.m_x[i][1]);
			mass += sys.m_mass[i];
		}
		const double vx = linearMomentum(sys)[0] / mass;
		// the penalty contact lets the body sink into the margin, but not through the plane
		EXPECT_GT(lowest, 0.0);
		EXPECT_LT(lowest, sys.m_colliderMargin);
		if (variant == 0) EXPECT_LT(fabs(vx), 1e-6);  // flat ground: at rest
		if (variant == 1) EXPECT_LT(fabs(vx), 2e-2);  // mu > tan(theta): sticks (creeps below the smoothing velocity)
		if (variant == 2)
		{
			// mu < tan(theta): slides with a = g (sin(theta) - mu cos(theta)); the cube lands during the first steps
			const double expected = 9.8 * (sin(theta) - mu * cos(theta)) * steps * kDt;
			EXPECT_NEAR(vx, expected, 0.03 * expected);
		}
		delete psb;
	}
}

//
// Projective dynamics / fast mass-spring
//

GTEST_TEST(VariationalProjectiveDynamics, ConvergesToTheNewtonSolution)
{
	btSoftBody* psb = createClampedBar();
	btVariationalARAP arap(2e5, 0.3);
	btVariationalSystem sys;
	sys.addBody(psb, &arap, 0, 0);
	sys.buildTopology();
	sys.m_gravity = btVec3d(0, -9.8, 0);
	btVariationalProjectedNewton newton;
	newton.m_tolerance = 1e-9;
	newton.m_linearTolerance = 1e-10;
	newton.m_maxIterations = 200;
	btVec3dArray reference, solution;
	runFromRest(sys, newton, 10, reference);
	double previous = 1e300;
	const int iterations[4] = {5, 20, 80, 320};
	for (int k = 0; k < 4; k++)
	{
		btVariationalProjectiveDynamics pd;
		pd.m_iterations = iterations[k];
		pd.m_tolerance = 0;
		runFromRest(sys, pd, 10, solution);
		const double difference = maxDifference(reference, solution);
		EXPECT_LT(difference, previous);  // more iterations: closer to the exact step
		previous = difference;
		EXPECT_EQ(1, pd.getLaplacian().numFactorizations());  // the system matrix is factorized once
	}
	EXPECT_LT(previous, 1e-6);
	delete psb;
}

GTEST_TEST(VariationalProjectiveDynamics, EnergyDecreasesMonotonically)
{
	// a mass-spring cloth pinned at two corners (Liu et al. 2013) and an ARAP solid (Bouaziz et al. 2014)
	for (int scene = 0; scene < 2; scene++)
	{
		btSoftBody* psb;
		if (scene == 0)
		{
			psb = btSoftBodyHelpers::CreatePatch(worldInfo(), btVector3(-0.5, 1, -0.5), btVector3(0.5, 1, -0.5), btVector3(-0.5, 1, 0.5), btVector3(0.5, 1, 0.5), 9, 9, 1 + 2, true);
			psb->generateBendingConstraints(2);
			psb->setTotalMass(0.3);
		}
		else
			psb = createClampedBar();
		btVariationalARAP arap(2e5, 0.3);
		btVariationalSystem sys;
		sys.addBody(psb, &arap, 500.0, 5.0);
		sys.buildTopology();
		sys.m_gravity = btVec3d(0, -9.8, 0);
		btVariationalProjectiveDynamics warmUp;
		warmUp.m_iterations = 20;
		btVec3dArray unused;
		runFromRest(sys, warmUp, 5, unused);
		sys.beginStep(kDt);
		// solve() starts every time from the same initial guess, so the energy after k iterations must
		// not increase with k
		double previous = 1e300;
		for (int its = 1; its <= 30; its++)
		{
			btVariationalProjectiveDynamics pd;
			pd.m_iterations = its;
			pd.m_tolerance = 0;
			for (int i = 0; i < sys.numNodes(); i++) sys.m_x[i] = sys.m_xPrev[i];
			pd.solve(sys);
			EXPECT_LE(pd.m_stats.m_energy, previous * (1 + 1e-12) + 1e-12) << "scene " << scene << ", " << its << " iterations";
			previous = pd.m_stats.m_energy;
		}
		sys.endStep();
		delete psb;
	}
}

//
// Quasi-Newton
//

GTEST_TEST(VariationalQuasiNewton, ConvergesToTheNewtonSolution)
{
	btSoftBody* psb = createClampedBar();
	btVariationalStableNeoHookean snh(2e5, 0.45);
	btVariationalSystem sys;
	sys.addBody(psb, &snh, 0, 0);
	sys.buildTopology();
	sys.m_gravity = btVec3d(0, -9.8, 0);
	btVariationalProjectedNewton newton;
	newton.m_tolerance = 1e-9;
	newton.m_linearTolerance = 1e-10;
	newton.m_maxIterations = 200;
	btVec3dArray reference, solution;
	runFromRest(sys, newton, 10, reference);
	int iterations[2];
	for (int k = 0; k < 2; k++)
	{
		btVariationalQuasiNewton quasiNewton;
		quasiNewton.m_historySize = k == 0 ? 0 : 5;  // 0: plain preconditioned descent on the Laplacian
		quasiNewton.m_tolerance = 1e-7;
		quasiNewton.m_maxIterations = 400;
		runFromRest(sys, quasiNewton, 10, solution);
		EXPECT_TRUE(quasiNewton.m_stats.m_converged);
		EXPECT_LT(maxDifference(reference, solution), 1e-5);
		iterations[k] = quasiNewton.m_stats.m_iterations;
	}
	EXPECT_LT(iterations[1], iterations[0]);  // the L-BFGS history pays off
	delete psb;
}

//
// Vertex block descent
//

GTEST_TEST(VariationalVertexBlockDescent, ConvergesToTheNewtonSolution)
{
	btSoftBody* psb = createClampedBar();
	btVariationalStableNeoHookean snh(2e5, 0.45);
	btVariationalSystem sys;
	sys.addBody(psb, &snh, 0, 0);
	sys.buildTopology();
	sys.m_gravity = btVec3d(0, -9.8, 0);
	btVariationalProjectedNewton newton;
	newton.m_tolerance = 1e-9;
	newton.m_linearTolerance = 1e-10;
	newton.m_maxIterations = 200;
	btVec3dArray reference, solution;
	runFromRest(sys, newton, 10, reference);
	double plain[3];
	const int sweeps[3] = {10, 40, 640};
	for (int k = 0; k < 3; k++)
	{
		btVariationalVertexBlockDescent vbd;
		vbd.m_iterations = sweeps[k];
		runFromRest(sys, vbd, 10, solution);
		plain[k] = maxDifference(reference, solution);
	}
	EXPECT_LT(plain[1], plain[0]);
	EXPECT_LT(plain[2], 1e-6);
	// Chebyshev acceleration gets closer with the same number of sweeps
	btVariationalVertexBlockDescent accelerated;
	accelerated.m_iterations = 40;
	accelerated.m_useChebyshev = true;
	accelerated.m_chebyshevRho = 0.95;
	runFromRest(sys, accelerated, 10, solution);
	EXPECT_LT(maxDifference(reference, solution), plain[1]);
	delete psb;
}

GTEST_TEST(VariationalVertexBlockDescent, SweepsNeverIncreaseTheEnergy)
{
	btSoftBody* psb = createClampedBar();
	btVariationalStableNeoHookean snh(2e5, 0.45);
	btVariationalSystem sys;
	sys.addBody(psb, &snh, 0, 0);
	sys.buildTopology();
	sys.m_gravity = btVec3d(0, -9.8, 0);
	btVariationalVertexBlockDescent vbd;
	vbd.m_iterations = 30;
	btVec3dArray unused;
	runFromRest(sys, vbd, 5, unused);
	sys.beginStep(kDt);
	vbd.buildColoring(sys);
	EXPECT_GE(vbd.numColors(), 4);  // a tetrahedron needs four colors
	for (int i = 0; i < sys.numNodes(); i++)
		if (!sys.m_fixed[i]) sys.m_x[i] = sys.m_xTilde[i];
	double previous = sys.energy(&sys.m_x[0]);
	for (int k = 0; k < 40; k++)
	{
		vbd.sweep(sys, &sys.m_x[0]);
		const double energy = sys.energy(&sys.m_x[0]);
		EXPECT_LE(energy, previous * (1 + 1e-12) + 1e-12);
		previous = energy;
	}
	sys.endStep();
	delete psb;
}

GTEST_TEST(VariationalVertexBlockDescent, FreeFallAndLanding)
{
	btSoftBody* psb = btVariationalHelpers::createTetGrid(worldInfo(), btVector3(0, 1.0, 0), btVector3(0.2, 0.2, 0.2), 2, 2, 2, 1000);
	btVariationalStableNeoHookean material(2e5, 0.4);
	btVariationalSystem sys;
	sys.addBody(psb, &material, 0, 0);
	sys.buildTopology();
	sys.m_gravity = btVec3d(0, -9.8, 0);
	Ground ground;
	ground.addTo(sys);
	btVariationalVertexBlockDescent vbd;
	vbd.m_iterations = 10;
	const double y0 = centreOfMass(sys, 0, sys.numNodes())[1];
	const int N = 18;
	for (int s = 0; s < N; s++) step(sys, vbd);
	// the adaptive initialization reproduces free fall (implicit Euler: g h^2 N (N + 1) / 2)
	EXPECT_NEAR(y0 - centreOfMass(sys, 0, sys.numNodes())[1], 9.8 * kDt * kDt * N * (N + 1) / 2.0, 1e-3);
	for (int s = 0; s < 240; s++) step(sys, vbd);
	double lowest = 1e9, fastest = 0;
	for (int i = 0; i < sys.numNodes(); i++)
	{
		lowest = btMin(lowest, sys.m_x[i][1]);
		fastest = btMax(fastest, sys.m_v[i].length());
	}
	EXPECT_GT(lowest, 0.0);
	EXPECT_LT(lowest, sys.m_colliderMargin);
	EXPECT_LT(fastest, 0.05);  // came to rest on the ground
	delete psb;
}

GTEST_TEST(VariationalVertexBlockDescent, AugmentedLagrangianKeepsHardJointsClosed)
{
	// 160 kg hang from a 0.8 kg block. Springs of 1e5 N/m stretch by centimetres under that load; the
	// augmented Lagrangian of AVBD closes the joints with the same initial stiffness.
	double error[2];
	for (int hard = 0; hard < 2; hard++)
	{
		HangingBlocks scene(hard == 1);
		btVariationalVertexBlockDescent vbd;
		vbd.m_iterations = 20;
		vbd.m_al.m_enabled = hard == 1;
		for (int s = 0; s < 240; s++) step(scene.m_sys, vbd);
		error[hard] = maxAttachmentError(scene.m_sys);
	}
	EXPECT_GT(error[0], 1e-2);
	EXPECT_LT(error[1], 5e-4);
	EXPECT_LT(error[1], 0.05 * error[0]);
}

//
// Distances and continuous collision detection
//

GTEST_TEST(VariationalDistance, MatchesReferenceAndFiniteDifferences)
{
	Rng rng(51);
	for (int edgeEdge = 0; edgeEdge < 2; edgeEdge++)
	{
		double worstValue = 0, worstWeights = 0, worstGradient = 0, worstHessian = 0, worstSymmetry = 0;
		int kinds[4] = {0, 0, 0, 0};
		for (int it = 0; it < 20000; it++)
		{
			btVec3d x[4];
			for (int j = 0; j < 4; j++) x[j] = rng.vec();
			if (it % 5 == 1) x[0] = x[0] * 3.0;                                  // far away: vertex and edge regions
			if (it % 5 == 2) x[0] = (x[1] + x[2] + x[3]) / 3.0 + rng.vec() * 0.2;  // near the middle
			const btVariationalClosest closest = edgeEdge ? btEdgeEdgeClosest(x[0], x[1], x[2], x[3]) : btPointTriangleClosest(x[0], x[1], x[2], x[3]);
			kinds[closest.m_kind]++;
			btVariationalPairDistance d;
			btEvaluatePairDistance(x, closest, 2, d);
			const double reference = edgeEdge ? referenceSegmentSegment(x[0], x[1], x[2], x[3]) : referencePointTriangle(x[0], x[1], x[2], x[3]);
			worstValue = btMax(worstValue, fabs(d.m_distanceSquared - reference) / (1.0 + reference));
			btVec3d difference(0, 0, 0);  // the weights reproduce the vector between the closest points
			for (int j = 0; j < 4; j++) difference += x[j] * d.m_weights[j];
			worstWeights = btMax(worstWeights, fabs(difference.length2() - d.m_distanceSquared) / (1.0 + reference));
			if (it % 20 != 0) continue;
			// derivatives with the closest feature pair frozen (the distance is smooth within a region)
			const double h = 1e-6;
			double gscale = 1e-9, hscale = 1e-9;
			for (int j = 0; j < 4; j++)
			{
				gscale = btMax(gscale, d.m_grad[j].maxAbs());
				for (int k = 0; k < 4; k++)
					for (int r = 0; r < 3; r++)
						for (int c = 0; c < 3; c++) hscale = btMax(hscale, fabs(d.m_hess[j][k].m[r][c]));
			}
			for (int j = 0; j < 4; j++)
				for (int a = 0; a < 3; a++)
				{
					btVec3d xp[4], xm[4];
					for (int k = 0; k < 4; k++) xp[k] = xm[k] = x[k];
					xp[j][a] += h;
					xm[j][a] -= h;
					btVariationalPairDistance dp, dm;
					btEvaluatePairDistance(xp, closest, 1, dp);
					btEvaluatePairDistance(xm, closest, 1, dm);
					worstGradient = btMax(worstGradient, fabs((dp.m_distanceSquared - dm.m_distanceSquared) / (2 * h) - d.m_grad[j][a]) / gscale);
					for (int k = 0; k < 4; k++)
						for (int b = 0; b < 3; b++)
						{
							worstHessian = btMax(worstHessian, fabs((dp.m_grad[k][b] - dm.m_grad[k][b]) / (2 * h) - d.m_hess[k][j].m[b][a]) / hscale);
							worstSymmetry = btMax(worstSymmetry, fabs(d.m_hess[k][j].m[b][a] - d.m_hess[j][k].m[a][b]) / hscale);
						}
				}
		}
		EXPECT_LT(worstValue, 1e-9);
		EXPECT_LT(worstWeights, 1e-9);
		EXPECT_LT(worstGradient, 1e-5);
		EXPECT_LT(worstHessian, 1e-4);
		EXPECT_LT(worstSymmetry, 1e-12);
		// every kind of closest feature pair was exercised
		EXPECT_GT(kinds[BT_CLOSEST_POINT_POINT], 100);
		EXPECT_GT(kinds[BT_CLOSEST_POINT_LINE], 100);
	}
}

GTEST_TEST(VariationalDistance, AdditiveCcdNeverPassesThrough)
{
	Rng rng(52);
	int cut = 0;
	double smallestRatio = 1e9;
	for (int it = 0; it < 4000; it++)
	{
		const bool edgeEdge = (it % 2) == 1;
		btVec3d x[4], p[4];
		for (int j = 0; j < 4; j++)
		{
			x[j] = rng.vec();
			p[j] = rng.vec() * 2.0;
		}
		const double d0 = sqrt(edgeEdge ? btEdgeEdgeDistanceSquared(x) : btPointTriangleDistanceSquared(x));
		if (d0 < 1e-3) continue;
		const double t = btAdditiveCCD(x, p, edgeEdge, 1.0);
		ASSERT_GT(t, 0.0);
		ASSERT_LE(t, 1.0);
		if (t < 1.0) cut++;
		// the distance stays positive along the whole step [0, t]
		double dmin = 1e9;
		for (int k = 0; k <= 100; k++)
		{
			btVec3d y[4];
			for (int j = 0; j < 4; j++) y[j] = x[j] + p[j] * (t * k / 100.0);
			dmin = btMin(dmin, sqrt(edgeEdge ? btEdgeEdgeDistanceSquared(y) : btPointTriangleDistanceSquared(y)));
		}
		smallestRatio = btMin(smallestRatio, dmin / d0);
	}
	EXPECT_GT(cut, 100);              // the test did produce approaching pairs
	EXPECT_GT(smallestRatio, 0.05);  // the step stops at a fraction of the initial distance, never at zero
}

//
// Incremental potential contact
//

GTEST_TEST(VariationalIPC, BarrierDerivatives)
{
	Rng rng(61);
	btVariationalStableNeoHookean material(2e5, 0.4);
	// two cubes a small gap apart, rotated so that corners face faces and edges cross edges
	btSoftBody* a = btVariationalHelpers::createTetGrid(worldInfo(), btVector3(0, 0.25, 0), btVector3(0.25, 0.25, 0.25), 2, 2, 2, 1000);
	btSoftBody* b = btVariationalHelpers::createTetGrid(worldInfo(), btVector3(0.1, 0.78, 0.05), btVector3(0.25, 0.25, 0.25), 2, 2, 2, 1000);
	const btVector3 centre(0.1, 0.78, 0.05);
	const btQuaternion q = btQuaternion(btVector3(1, 0, 0), 0.08) * btQuaternion(btVector3(0, 1, 0), 0.5);
	for (int i = 0; i < b->m_nodes.size(); i++)
	{
		b->m_nodes[i].m_x = centre + quatRotate(q, b->m_nodes[i].m_x - centre);
		b->m_nodes[i].m_q = b->m_nodes[i].m_x;
	}
	btVariationalSystem sys;
	sys.addBody(a, &material, 0, 0);
	sys.addBody(b, &material, 0, 0);
	sys.buildTopology();
	sys.m_bodies[0].m_friction = sys.m_bodies[1].m_friction = 0.5;
	sys.m_gravity = btVec3d(0, -9.8, 0);
	btVariationalIPC ipc;
	ipc.m_dhat = 0.12;
	sys.m_contact = &ipc;
	for (int i = 0; i < sys.numNodes(); i++) sys.m_v[i] = rng.vec() * 0.3;
	sys.beginStep(kDt);
	EXPECT_TRUE(ipc.initialStateValid());
	const int n = sys.numNodes();
	btVec3dArray x, g, gp, gm, v, Hv, xp, xm;
	x.resize(n);
	g.resize(n);
	gp.resize(n);
	gm.resize(n);
	v.resize(n);
	Hv.resize(n);
	xp.resize(n);
	xm.resize(n);
	for (int i = 0; i < n; i++) x[i] = sys.m_x[i] + rng.vec() * 0.004;
	double smallestDistance;
	EXPECT_GT(ipc.countActivePairs(&x[0], &smallestDistance), 10);
	EXPECT_GT(smallestDistance, 0.0);
	EXPECT_FALSE(ipc.hasIntersections(&x[0]));
	sys.gradient(&x[0], &g[0]);
	const double h = 1e-7;
	double worstGradient = 0, gscale = 1e-12;
	for (int i = 0; i < n; i++) gscale = btMax(gscale, g[i].maxAbs());
	for (int s = 0; s < 150; s++)
	{
		const int i = rng.index(n), c = rng.index(3);
		const double keep = x[i][c];
		x[i][c] = keep + h;
		const double Ep = sys.energy(&x[0]);
		x[i][c] = keep - h;
		const double Em = sys.energy(&x[0]);
		x[i][c] = keep;
		worstGradient = btMax(worstGradient, fabs((Ep - Em) / (2 * h) - g[i][c]) / gscale);
	}
	EXPECT_LT(worstGradient, 1e-5);
	btVariationalHessian H;
	sys.hessian(&x[0], H, BT_VARIATIONAL_FILTER_NONE);
	for (int i = 0; i < n; i++) v[i] = rng.vec();
	H.multiply(&v[0], &Hv[0], &sys.m_fixed[0]);
	for (int i = 0; i < n; i++)
	{
		xp[i] = x[i] + v[i] * h;
		xm[i] = x[i] - v[i] * h;
	}
	sys.gradient(&xp[0], &gp[0]);
	sys.gradient(&xm[0], &gm[0]);
	double worstHessian = 0, hscale = 1e-12;
	for (int i = 0; i < n; i++) hscale = btMax(hscale, Hv[i].maxAbs());
	for (int i = 0; i < n; i++) worstHessian = btMax(worstHessian, ((gp[i] - gm[i]) * (1.0 / (2 * h)) - Hv[i]).maxAbs() / hscale);
	EXPECT_LT(worstHessian, 1e-5);
	// the filtered Hessian is positive definite along random directions, and its block sparse form
	// (many contact stencils overlap) gives the same products
	btVariationalHessian filtered;
	sys.hessian(&x[0], filtered, BT_VARIATIONAL_FILTER_CLAMP);
	filtered.compress(&sys.m_fixed[0]);
	EXPECT_LT(filtered.numCompressedBlocks(), 16 * filtered.m_stencils.size());
	for (int k = 0; k < 30; k++)
	{
		for (int i = 0; i < n; i++) v[i] = rng.vec();
		filtered.multiply(&v[0], &Hv[0], &sys.m_fixed[0]);
		filtered.multiplyCompressed(&v[0], &gp[0]);
		double quadratic = 0, scale = 1e-12, difference = 0;
		for (int i = 0; i < n; i++)
		{
			quadratic += v[i].dot(Hv[i]);
			scale = btMax(scale, Hv[i].maxAbs());
			difference = btMax(difference, (Hv[i] - gp[i]).maxAbs());
		}
		EXPECT_GT(quadratic, 0.0);
		EXPECT_LT(difference, 1e-12 * scale);
	}
	sys.m_contact = 0;
	delete a;
	delete b;
}

namespace
{
// a cube dropped onto a cube that rests on the ground; returns the number of steps with an intersection
struct StackResult
{
	int m_intersectingSteps, m_nodesBelowGround;
	double m_smallestDistance, m_upperCubeHeight, m_dhat;
};

StackResult dropCubeOnCube(btVariationalIntegrator& integrator, int steps)
{
	btVariationalStableNeoHookean material(2e5, 0.4);
	btSoftBody* a = btVariationalHelpers::createTetGrid(worldInfo(), btVector3(0, 0.26, 0), btVector3(0.25, 0.25, 0.25), 2, 2, 2, 1000);
	btSoftBody* b = btVariationalHelpers::createTetGrid(worldInfo(), btVector3(0.03, 0.9, 0.02), btVector3(0.2, 0.2, 0.2), 2, 2, 2, 1000);
	btVariationalSystem sys;
	sys.addBody(a, &material, 0, 0);
	sys.addBody(b, &material, 0, 0);
	sys.buildTopology();
	sys.m_bodies[0].m_friction = sys.m_bodies[1].m_friction = 0.5;
	sys.m_gravity = btVec3d(0, -9.8, 0);
	sys.m_damping = 1.0;
	sys.m_colliderModel = btVariationalSystem::COLLIDER_BARRIER;
	sys.m_colliderMargin = 0.005;
	sys.m_colliderStiffness = 100.0;
	Ground ground;
	ground.addTo(sys);
	btVariationalIPC ipc;
	ipc.m_dhat = 0.005;
	sys.m_contact = &ipc;
	StackResult result;
	result.m_intersectingSteps = result.m_nodesBelowGround = 0;
	result.m_smallestDistance = 1e9;
	for (int s = 0; s < steps; s++)
	{
		step(sys, integrator);
		if (ipc.hasIntersections(&sys.m_x[0])) result.m_intersectingSteps++;
		double smallest;
		ipc.countActivePairs(&sys.m_x[0], &smallest);
		result.m_smallestDistance = btMin(result.m_smallestDistance, smallest);
		for (int i = 0; i < sys.numNodes(); i++)
			if (sys.m_x[i][1] <= 0.0) result.m_nodesBelowGround++;
	}
	result.m_upperCubeHeight = centreOfMass(sys, sys.m_bodies[1].m_firstNode, sys.m_bodies[1].m_numNodes)[1];
	result.m_dhat = ipc.getDhat();
	sys.m_contact = 0;
	delete a;
	delete b;
	return result;
}
}  // namespace

GTEST_TEST(VariationalIPC, StackStaysIntersectionFreeWithNewton)
{
	btVariationalProjectedNewton newton;
	newton.m_tolerance = 1e-3;
	const StackResult r = dropCubeOnCube(newton, 180);
	EXPECT_EQ(0, r.m_intersectingSteps);
	EXPECT_EQ(0, r.m_nodesBelowGround);
	EXPECT_GT(r.m_smallestDistance, 0.0);
	// The lower cube (edge 0.5) hovers a few millimetres above the ground and is compressed by about
	// 12 mm under its own weight and the load; the upper cube (edge 0.4) rests a few millimetres above it.
	EXPECT_NEAR(r.m_upperCubeHeight, 0.693, 0.015);
}

GTEST_TEST(VariationalIPC, StackStaysIntersectionFreeWithVertexBlockDescent)
{
	btVariationalVertexBlockDescent vbd;
	vbd.m_iterations = 30;
	const StackResult r = dropCubeOnCube(vbd, 180);
	EXPECT_EQ(0, r.m_intersectingSteps);
	EXPECT_EQ(0, r.m_nodesBelowGround);
	EXPECT_GT(r.m_smallestDistance, 0.0);
	EXPECT_NEAR(r.m_upperCubeHeight, 0.693, 0.015);
}

GTEST_TEST(VariationalIPC, CollisionGroupsExcludeJoinedBodies)
{
	// two cubes that touch along a face: a barrier cannot start from zero distance, unless the pair is
	// taken out of the contact model through a common collision group
	btVariationalStableNeoHookean material(2e5, 0.4);
	for (int grouped = 0; grouped < 2; grouped++)
	{
		btSoftBody* a = btVariationalHelpers::createTetGrid(worldInfo(), btVector3(0, 1.0, 0), btVector3(0.1, 0.1, 0.1), 2, 2, 2, 1000);
		btSoftBody* b = btVariationalHelpers::createTetGrid(worldInfo(), btVector3(0, 0.8, 0), btVector3(0.1, 0.1, 0.1), 2, 2, 2, 1000);
		btVariationalSystem sys;
		sys.addBody(a, &material, 0, 0);
		sys.addBody(b, &material, 0, 0);
		sys.buildTopology();
		btVariationalIPC ipc;
		ipc.m_dhat = 0.005;
		if (grouped)
		{
			ipc.setCollisionGroup(a, 1);
			ipc.setCollisionGroup(b, 1);
		}
		sys.m_contact = &ipc;
		sys.beginStep(kDt);
		EXPECT_EQ(grouped == 1, ipc.initialStateValid());
		if (grouped)
		{
			double smallest;
			EXPECT_EQ(0, ipc.countActivePairs(&sys.m_x[0], &smallest));
		}
		sys.endStep();
		sys.m_contact = 0;
		delete a;
		delete b;
	}
}

GTEST_TEST(VariationalIPC, PairsAndSelfContactCanBeExcluded)
{
	btVariationalStableNeoHookean material(2e5, 0.4);
	// two cubes that touch along a face, as two links of a chain do at their joint
	for (int excluded = 0; excluded < 2; excluded++)
	{
		btSoftBody* a = btVariationalHelpers::createTetGrid(worldInfo(), btVector3(0, 1.0, 0), btVector3(0.1, 0.1, 0.1), 2, 2, 2, 1000);
		btSoftBody* b = btVariationalHelpers::createTetGrid(worldInfo(), btVector3(0, 0.8, 0), btVector3(0.1, 0.1, 0.1), 2, 2, 2, 1000);
		btSoftBody* c = btVariationalHelpers::createTetGrid(worldInfo(), btVector3(0, 0.598, 0), btVector3(0.1, 0.1, 0.1), 2, 2, 2, 1000);
		btVariationalSystem sys;
		sys.addBody(a, &material, 0, 0);
		sys.addBody(b, &material, 0, 0);
		sys.addBody(c, &material, 0, 0);
		sys.buildTopology();
		btVariationalIPC ipc;
		ipc.m_dhat = 0.005;
		if (excluded) ipc.excludePair(a, b);
		sys.m_contact = &ipc;
		sys.beginStep(kDt);
		EXPECT_EQ(excluded == 1, ipc.initialStateValid());
		if (excluded)
		{
			// the third cube, 2 mm below the second, is still in contact with it
			double smallest;
			EXPECT_GT(ipc.countActivePairs(&sys.m_x[0], &smallest), 0);
			EXPECT_NEAR(smallest, 0.002, 1e-6);
		}
		sys.endStep();
		sys.m_contact = 0;
		delete a;
		delete b;
		delete c;
	}
	// a ring whose hole is narrower than the activation distance touches itself in its rest shape
	for (int excluded = 0; excluded < 2; excluded++)
	{
		btSoftBody* ring = btVariationalHelpers::createTetTorus(worldInfo(), btVector3(0, 1, 0), 0.051, 0.05, 16, 2, 1000);
		btVariationalSystem sys;
		sys.addBody(ring, &material, 0, 0);
		sys.buildTopology();
		btVariationalIPC ipc;
		ipc.m_dhat = 0.005;
		if (excluded) ipc.excludePair(ring, ring);
		sys.m_contact = &ipc;
		sys.beginStep(kDt);
		double smallest;
		const int active = ipc.countActivePairs(&sys.m_x[0], &smallest);
		if (excluded)
			EXPECT_EQ(0, active);
		else
			EXPECT_GT(active, 0);
		sys.endStep();
		sys.m_contact = 0;
		delete ring;
	}
}

//
// Hybrid integrator
//

GTEST_TEST(VariationalHybrid, TakesTheQuasiNewtonPathWithoutContact)
{
	btSoftBody* psb = createClampedBar();
	btVariationalStableNeoHookean snh(2e5, 0.45);
	btVariationalSystem sys;
	sys.addBody(psb, &snh, 0, 0);
	sys.buildTopology();
	sys.m_gravity = btVec3d(0, -9.8, 0);
	btVariationalProjectedNewton newton;
	newton.m_tolerance = 1e-9;
	newton.m_linearTolerance = 1e-10;
	newton.m_maxIterations = 200;
	btVec3dArray reference, solution;
	runFromRest(sys, newton, 10, reference);
	btVariationalHybrid hybrid;
	hybrid.m_tolerance = 1e-4;
	runFromRest(sys, hybrid, 10, solution);
	EXPECT_EQ(btVariationalHybrid::PATH_QUASI_NEWTON, hybrid.lastPath());
	EXPECT_TRUE(hybrid.m_stats.m_converged);
	EXPECT_LT(maxDifference(reference, solution), 1e-4);
	// with the fast path switched off the same step is solved by Newton
	hybrid.m_useQuasiNewton = false;
	runFromRest(sys, hybrid, 10, solution);
	EXPECT_EQ(btVariationalHybrid::PATH_NEWTON, hybrid.lastPath());
	EXPECT_LT(maxDifference(reference, solution), 1e-4);
	const int plainIterations = hybrid.m_stats.m_linearIterations;
	// ... optionally with the projective dynamics matrix as preconditioner, which needs fewer
	// conjugate gradient iterations for the same result
	hybrid.m_useLaplacianPreconditioner = true;
	runFromRest(sys, hybrid, 10, solution);
	EXPECT_EQ(btVariationalHybrid::PATH_NEWTON_LAPLACIAN, hybrid.lastPath());
	EXPECT_LT(maxDifference(reference, solution), 1e-4);
	EXPECT_LT(hybrid.m_stats.m_linearIterations, plainIterations);
	delete psb;
}

GTEST_TEST(VariationalHybrid, ClosesHardJointsThatNewtonAloneLeavesOpen)
{
	double error[2];
	for (int k = 0; k < 2; k++)
	{
		HangingBlocks scene(true);
		btVariationalProjectedNewton newton;  // treats hard attachments as springs
		newton.m_tolerance = 1e-2;
		btVariationalHybrid hybrid;
		btVariationalIntegrator* integrator = k == 0 ? (btVariationalIntegrator*)&newton : (btVariationalIntegrator*)&hybrid;
		for (int s = 0; s < 120; s++) step(scene.m_sys, *integrator);
		error[k] = maxAttachmentError(scene.m_sys);
		if (k == 1)
		{
			EXPECT_EQ(btVariationalHybrid::PATH_NEWTON, hybrid.lastPath());
			EXPECT_NEAR(error[k], hybrid.maxConstraintError(), 1e-12);
		}
	}
	EXPECT_GT(error[0], 1e-2);
	EXPECT_LT(error[1], 5e-4);
}

GTEST_TEST(VariationalHybrid, TakesTheNewtonPathInContact)
{
	btVariationalHybrid hybrid;
	const StackResult r = dropCubeOnCube(hybrid, 120);
	EXPECT_EQ(btVariationalHybrid::PATH_NEWTON, hybrid.lastPath());
	EXPECT_EQ(0, r.m_intersectingSteps);
	EXPECT_EQ(0, r.m_nodesBelowGround);
	EXPECT_NEAR(r.m_upperCubeHeight, 0.693, 0.015);
}

//
// Integration into btDeformableMultiBodyDynamicsWorld
//

GTEST_TEST(VariationalDeformableSolver, StepsSoftBodiesInABulletWorld)
{
	btSoftBodyRigidBodyCollisionConfiguration* collisionConfiguration = new btSoftBodyRigidBodyCollisionConfiguration();
	btCollisionDispatcher* dispatcher = new btCollisionDispatcher(collisionConfiguration);
	btDbvtBroadphase* broadphase = new btDbvtBroadphase();
	btVariationalDeformableSolver* deformableSolver = new btVariationalDeformableSolver();
	btDeformableMultiBodyConstraintSolver* solver = new btDeformableMultiBodyConstraintSolver();
	solver->setDeformableSolver(deformableSolver);
	btDeformableMultiBodyDynamicsWorld* world = new btDeformableMultiBodyDynamicsWorld(dispatcher, broadphase, solver, collisionConfiguration, deformableSolver);
	world->setGravity(btVector3(0, -9.8, 0));
	world->getWorldInfo().m_gravity = btVector3(0, -9.8, 0);
	world->getWorldInfo().m_sparsesdf.setDefaultVoxelsz(0.25);
	world->getWorldInfo().m_sparsesdf.Reset();

	// static ground: a box whose top face is the plane y = 0
	btBoxShape* groundShape = new btBoxShape(btVector3(50, 1, 50));
	btTransform groundTransform;
	groundTransform.setIdentity();
	groundTransform.setOrigin(btVector3(0, -1, 0));
	btRigidBody* ground = new btRigidBody(btRigidBody::btRigidBodyConstructionInfo(0, 0, groundShape));
	ground->setWorldTransform(groundTransform);
	ground->setFriction(0.8);
	world->addRigidBody(ground);

	btVariationalStableNeoHookean material(2e5, 0.4);
	btSoftBody* falling = btVariationalHelpers::createTetGrid(world->getWorldInfo(), btVector3(0, 0.5, 0), btVector3(0.2, 0.2, 0.2), 2, 2, 2, 1000);
	btSoftBody* hanging = btVariationalHelpers::createTetGrid(world->getWorldInfo(), btVector3(2, 1.5, 0), btVector3(0.1, 0.1, 0.1), 2, 2, 2, 1000);
	btSoftBody* bodies[2] = {falling, hanging};
	for (int k = 0; k < 2; k++)
	{
		bodies[k]->getCollisionShape()->setMargin(0.01);
		bodies[k]->m_cfg.kDF = 0.5;
		world->addSoftBody(bodies[k]);
		deformableSolver->setMaterial(bodies[k], &material);
	}
	// the second cube hangs from one of its nodes
	const btVector3 anchor = hanging->m_nodes[0].m_x;
	deformableSolver->addAttachment(hanging, 0, anchor, 1e5, true);

	for (int s = 0; s < 150; s++) world->stepSimulation(btScalar(kDt), 0, btScalar(kDt));

	EXPECT_EQ(2, deformableSolver->getSystem().m_bodies.size());
	// the first cube has landed: its nodes sit just above the ground and it kept its shape
	btScalar lowest = btScalar(1e9), highest = btScalar(-1e9);
	for (int i = 0; i < falling->m_nodes.size(); i++)
	{
		lowest = btMin(lowest, falling->m_nodes[i].m_x.y());
		highest = btMax(highest, falling->m_nodes[i].m_x.y());
		EXPECT_LT(falling->m_nodes[i].m_v.length(), btScalar(0.05));
	}
	EXPECT_GT(lowest, btScalar(-1e-3));
	EXPECT_LT(lowest, btScalar(0.02));
	EXPECT_NEAR(highest - lowest, btScalar(0.4), btScalar(0.02));
	// the second cube still hangs from its anchor
	EXPECT_LT((hanging->m_nodes[0].m_x - anchor).length(), btScalar(2e-3));
	btScalar lowestHanging = btScalar(1e9);
	for (int i = 0; i < hanging->m_nodes.size(); i++) lowestHanging = btMin(lowestHanging, hanging->m_nodes[i].m_x.y());
	EXPECT_GT(lowestHanging, btScalar(0.9));

	// removing a body that an attachment refers to, and tearing the world down, must both be safe
	world->removeSoftBody(hanging);
	delete hanging;
	for (int s = 0; s < 5; s++) world->stepSimulation(btScalar(kDt), 0, btScalar(kDt));
	EXPECT_EQ(1, deformableSolver->getSystem().m_bodies.size());
	world->removeSoftBody(falling);
	delete falling;
	world->removeRigidBody(ground);
	delete ground;
	delete groundShape;
	delete world;
	delete solver;
	delete deformableSolver;
	delete broadphase;
	delete dispatcher;
	delete collisionConfiguration;
}

GTEST_TEST(VariationalDeformableSolver, JointsSurviveNewBodiesAndCanBeReleased)
{
	btSoftBodyRigidBodyCollisionConfiguration* collisionConfiguration = new btSoftBodyRigidBodyCollisionConfiguration();
	btCollisionDispatcher* dispatcher = new btCollisionDispatcher(collisionConfiguration);
	btDbvtBroadphase* broadphase = new btDbvtBroadphase();
	btVariationalDeformableSolver* deformableSolver = new btVariationalDeformableSolver();
	btDeformableMultiBodyConstraintSolver* solver = new btDeformableMultiBodyConstraintSolver();
	solver->setDeformableSolver(deformableSolver);
	btDeformableMultiBodyDynamicsWorld* world = new btDeformableMultiBodyDynamicsWorld(dispatcher, broadphase, solver, collisionConfiguration, deformableSolver);
	world->setGravity(btVector3(0, -9.8, 0));
	world->getWorldInfo().m_gravity = btVector3(0, -9.8, 0);
	world->getWorldInfo().m_sparsesdf.Initialize();
	btVariationalHybrid hybrid;
	deformableSolver->setIntegrator(&hybrid);

	// 8 kg hang from one node on a hard attachment that starts as a soft spring: only the dual
	// variable and the grown penalty of the augmented Lagrangian keep it from sagging by centimetres
	btVariationalStableNeoHookean material(1e6, 0.4);
	btSoftBody* hanging = btVariationalHelpers::createTetGrid(world->getWorldInfo(), btVector3(0, 2, 0), btVector3(0.1, 0.1, 0.1), 2, 2, 2, 1000);
	world->addSoftBody(hanging);
	deformableSolver->setMaterial(hanging, &material);
	const int node = 1 + 3 * (2 + 3 * 1);  // middle of the top face
	const btVector3 anchor = hanging->m_nodes[node].m_x;
	const int attachment = deformableSolver->addAttachment(hanging, node, anchor, 2e3, true);
	for (int s = 0; s < 90; s++) world->stepSimulation(btScalar(kDt), 0, btScalar(kDt));
	EXPECT_LT((hanging->m_nodes[node].m_x - anchor).length(), btScalar(3e-4));
	const double penalty = deformableSolver->getSystem().m_attachments[attachment].m_penalty;
	EXPECT_GT(penalty, 1e4);

	// a body that enters the world rebuilds the system; the joint must not notice
	btSoftBody* other = btVariationalHelpers::createTetGrid(world->getWorldInfo(), btVector3(3, 2, 0), btVector3(0.1, 0.1, 0.1), 1, 1, 1, 1000);
	world->addSoftBody(other);
	deformableSolver->setMaterial(other, &material);
	for (int s = 0; s < 3; s++)
	{
		world->stepSimulation(btScalar(kDt), 0, btScalar(kDt));
		EXPECT_LT((hanging->m_nodes[node].m_x - anchor).length(), btScalar(3e-4));
	}
	EXPECT_EQ(2, deformableSolver->getSystem().m_bodies.size());
	EXPECT_GE(deformableSolver->getSystem().m_attachments[attachment].m_penalty, 0.9 * penalty);

	// released, it falls
	deformableSolver->setAttachmentEnabled(attachment, false);
	for (int s = 0; s < 30; s++) world->stepSimulation(btScalar(kDt), 0, btScalar(kDt));
	EXPECT_LT(hanging->m_nodes[node].m_x.y(), anchor.y() - btScalar(1.0));

	world->removeSoftBody(hanging);
	world->removeSoftBody(other);
	delete hanging;
	delete other;
	delete world;
	delete solver;
	delete deformableSolver;
	delete broadphase;
	delete dispatcher;
	delete collisionConfiguration;
}

int main(int argc, char** argv)
{
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
