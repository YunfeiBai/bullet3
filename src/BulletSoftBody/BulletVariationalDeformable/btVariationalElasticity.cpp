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

#include "btVariationalElasticity.h"

void btVariationalMaterial::hessian(const btMat3d& F, btMat9d& H, btVariationalHessianFilter filter) const
{
	btVariationalEigenSystem es;
	eigenSystem(F, es);
	H.setZero();
	for (int i = 0; i < 9; i++)
	{
		double l = btVariationalFilterEigenvalue(es.m_lambda[i], filter);
		if (l == 0.0) continue;
		btVec9d q;
		for (int c = 0; c < 3; c++)
			for (int r = 0; r < 3; r++) q[r + 3 * c] = es.m_Q[i].m[r][c];
		H.addScaledOuter(l, q, q);
	}
}

void btVariationalMaterial::vertexBlock(const btMat3d& F, const btVec3d& g, btVec3d& Pg, btMat3d& block) const
{
	Pg = firstPiola(F) * g;
	btVariationalEigenSystem es;
	eigenSystem(F, es);
	block.setZero();
	for (int i = 0; i < 9; i++)
	{
		if (es.m_lambda[i] <= 0.0) continue;
		const btVec3d w = es.m_Q[i] * g;
		block.addScaledOuter(es.m_lambda[i], w, w);
	}
}

void btVariationalMaterial::buildTwistAndFlip(const btMat3d& U, const btMat3d& V, btMat3d* Q)
{
	const double s = 0.70710678118654752440;
	const btVec3d u0 = U.col(0), u1 = U.col(1), u2 = U.col(2);
	const btVec3d v0 = V.col(0), v1 = V.col(1), v2 = V.col(2);
	// twists: infinitesimal rotations about the three principal axes
	Q[0] = (btMat3d::outer(u2, v1) - btMat3d::outer(u1, v2)) * s;
	Q[1] = (btMat3d::outer(u0, v2) - btMat3d::outer(u2, v0)) * s;
	Q[2] = (btMat3d::outer(u0, v1) - btMat3d::outer(u1, v0)) * s;
	// flips: the symmetric counterparts
	Q[3] = (btMat3d::outer(u1, v2) + btMat3d::outer(u2, v1)) * s;
	Q[4] = (btMat3d::outer(u0, v2) + btMat3d::outer(u2, v0)) * s;
	Q[5] = (btMat3d::outer(u0, v1) + btMat3d::outer(u1, v0)) * s;
}

//
// Stable Neo-Hookean
//

btVariationalStableNeoHookean::btVariationalStableNeoHookean(double youngsModulus, double poissonRatio, bool logTerm)
{
	setParameters(youngsModulus, poissonRatio, logTerm);
}

void btVariationalStableNeoHookean::setParameters(double youngsModulus, double poissonRatio, bool logTerm)
{
	m_youngsModulus = youngsModulus;
	m_poissonRatio = poissonRatio;
	m_logTerm = logTerm;
	const double muLame = youngsModulus / (2.0 * (1.0 + poissonRatio));
	const double lambdaLame = youngsModulus * poissonRatio / ((1.0 + poissonRatio) * (1.0 - 2.0 * poissonRatio));
	// Section 3.4 of the paper: choose the constants so that the linearization at rest is
	// linear elasticity with the given Lame parameters.
	if (m_logTerm)
	{
		m_mu = (4.0 / 3.0) * muLame;
		m_lambda = lambdaLame + (5.0 / 6.0) * muLame;
		m_alpha = 1.0 + m_mu / m_lambda - m_mu / (4.0 * m_lambda);
	}
	else
	{
		m_mu = muLame;
		m_lambda = lambdaLame + muLame;
		m_alpha = 1.0 + m_mu / m_lambda;
	}
}

double btVariationalStableNeoHookean::energyDensity(const btMat3d& F) const
{
	const double Ic = F.frobeniusNorm2();
	const double Jm = F.det() - m_alpha;
	double psi = 0.5 * m_mu * (Ic - 3.0) + 0.5 * m_lambda * Jm * Jm;
	if (m_logTerm) psi -= 0.5 * m_mu * log(Ic + 1.0);
	// shift so that the rest state has zero energy (does not change forces)
	const double JmRest = 1.0 - m_alpha;
	double rest = 0.5 * m_lambda * JmRest * JmRest;
	if (m_logTerm) rest -= 0.5 * m_mu * log(4.0);
	return psi - rest;
}

btMat3d btVariationalStableNeoHookean::firstPiola(const btMat3d& F) const
{
	const double Ic = F.frobeniusNorm2();
	const double muEff = m_logTerm ? m_mu * (1.0 - 1.0 / (Ic + 1.0)) : m_mu;
	return F * muEff + F.cofactor() * (m_lambda * (F.det() - m_alpha));
}

void btVariationalStableNeoHookean::eigenSystem(const btMat3d& F, btVariationalEigenSystem& es) const
{
	btMat3d U, V;
	btVec3d S;
	btSvd3d(F, U, S, V);
	const double J = S[0] * S[1] * S[2];
	const double Ic = S.length2();
	const double muEff = m_logTerm ? m_mu * (1.0 - 1.0 / (Ic + 1.0)) : m_mu;
	const double cF = m_logTerm ? 2.0 * m_mu / ((Ic + 1.0) * (Ic + 1.0)) : 0.0;
	const double cJ = m_lambda * (J - m_alpha);

	// d2Psi/dF2 = muEff I + cF f f^T + lambda g g^T + cJ H_J with f = vec(F), g = vec(dJ/dF).
	// The twist/flip modes are eigenvectors of H_J (eigenvalues +sigma_k / -sigma_k) orthogonal to f and g.
	buildTwistAndFlip(U, V, es.m_Q);
	for (int k = 0; k < 3; k++)
	{
		es.m_lambda[k] = muEff + cJ * S[k];
		es.m_lambda[3 + k] = muEff - cJ * S[k];
	}
	// the remaining three modes live in the span of the scaling directions; reduce to a 3x3 problem
	const btVec3d c(S[1] * S[2], S[0] * S[2], S[0] * S[1]);
	btMatNd<3> A;
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			A(i, j) = cF * S[i] * S[j] + m_lambda * c[i] * c[j];
	A(0, 1) += cJ * S[2];
	A(1, 0) += cJ * S[2];
	A(0, 2) += cJ * S[1];
	A(2, 0) += cJ * S[1];
	A(1, 2) += cJ * S[0];
	A(2, 1) += cJ * S[0];
	for (int i = 0; i < 3; i++) A(i, i) += muEff;
	double ev[3];
	btMatNd<3> W;
	btSymmetricEigenJacobi(A, ev, W);
	const btVec3d u[3] = {U.col(0), U.col(1), U.col(2)};
	const btVec3d v[3] = {V.col(0), V.col(1), V.col(2)};
	for (int i = 0; i < 3; i++)
	{
		es.m_lambda[6 + i] = ev[i];
		btMat3d Q;
		for (int a = 0; a < 3; a++) Q.addScaledOuter(W(a, i), u[a], v[a]);
		es.m_Q[6 + i] = Q;
	}
}

void btVariationalStableNeoHookean::vertexBlock(const btMat3d& F, const btVec3d& g, btVec3d& Pg, btMat3d& block) const
{
	const double Ic = F.frobeniusNorm2();
	const double muEff = m_logTerm ? m_mu * (1.0 - 1.0 / (Ic + 1.0)) : m_mu;
	const double cF = m_logTerm ? 2.0 * m_mu / ((Ic + 1.0) * (Ic + 1.0)) : 0.0;
	const btVec3d Fg = F * g;
	const btVec3d Cg = F.cofactor() * g;
	Pg = Fg * muEff + Cg * (m_lambda * (F.det() - m_alpha));
	// det F is linear in the position of one vertex, hence the H_J term of the Hessian drops out
	block.setZero();
	block.addScaledIdentity(muEff * g.length2());
	block.addScaledOuter(m_lambda, Cg, Cg);
	if (cF != 0.0) block.addScaledOuter(cF, Fg, Fg);
}

double btVariationalStableNeoHookean::laplacianStiffness() const
{
	// The quasi-Newton paper weights the Laplacian with (approximately) 2 mu + lambda in Lame
	// parameters, i.e. the P-wave modulus of the linearized material.
	const double muLame = m_youngsModulus / (2.0 * (1.0 + m_poissonRatio));
	const double lambdaLame = m_youngsModulus * m_poissonRatio / ((1.0 + m_poissonRatio) * (1.0 - 2.0 * m_poissonRatio));
	return 2.0 * muLame + lambdaLame;
}

//
// ARAP
//

btVariationalARAP::btVariationalARAP(double youngsModulus, double poissonRatio)
{
	setParameters(youngsModulus, poissonRatio);
}

void btVariationalARAP::setParameters(double youngsModulus, double poissonRatio)
{
	m_youngsModulus = youngsModulus;
	m_poissonRatio = poissonRatio;
	m_mu = youngsModulus / (2.0 * (1.0 + poissonRatio));
}

double btVariationalARAP::energyDensity(const btMat3d& F) const
{
	return m_mu * (F - btPolarRotation3d(F)).frobeniusNorm2();
}

btMat3d btVariationalARAP::firstPiola(const btMat3d& F) const
{
	return (F - btPolarRotation3d(F)) * (2.0 * m_mu);
}

void btVariationalARAP::eigenSystem(const btMat3d& F, btVariationalEigenSystem& es) const
{
	btMat3d U, V;
	btVec3d S;
	btSvd3d(F, U, S, V);
	buildTwistAndFlip(U, V, es.m_Q);
	const btVec3d u[3] = {U.col(0), U.col(1), U.col(2)};
	const btVec3d v[3] = {V.col(0), V.col(1), V.col(2)};
	for (int i = 0; i < 3; i++)
	{
		es.m_Q[6 + i] = btMat3d::outer(u[i], v[i]);
		es.m_lambda[6 + i] = 2.0 * m_mu;
		es.m_lambda[3 + i] = 2.0 * m_mu;
	}
	// only the twist modes see the rotation gradient dR/dF, whose eigenvalue is 2 / (sigma_i + sigma_j)
	const double tiny = 1e-8;
	const double d0 = S[1] + S[2], d1 = S[0] + S[2], d2 = S[0] + S[1];
	es.m_lambda[0] = 2.0 * m_mu * (1.0 - 2.0 / (d0 > tiny ? d0 : tiny));
	es.m_lambda[1] = 2.0 * m_mu * (1.0 - 2.0 / (d1 > tiny ? d1 : tiny));
	es.m_lambda[2] = 2.0 * m_mu * (1.0 - 2.0 / (d2 > tiny ? d2 : tiny));
}

bool btVariationalARAP::projectiveDynamicsTarget(const btMat3d& F, btMat3d& target, double& stiffness) const
{
	target = btPolarRotation3d(F);
	stiffness = 2.0 * m_mu;
	return true;
}
