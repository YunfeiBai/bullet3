#ifndef VARIATIONAL_SCENE_HELPERS_H
#define VARIATIONAL_SCENE_HELPERS_H

#include "VariationalDemoBase.h"

///Small deterministic random number generator so that the demo scenes are reproducible
struct VariationalDemoRandom
{
	unsigned int m_seed;
	VariationalDemoRandom(unsigned int seed) : m_seed(seed) {}
	///uniform in [-1, 1]
	btScalar next()
	{
		m_seed = 1664525u * m_seed + 1013904223u;
		return btScalar((m_seed >> 8) & 0xffff) / btScalar(32767.5) - btScalar(1);
	}
};

///rotates a soft body about its own centre
inline void variationalRotateAboutCentre(btSoftBody* psb, const btQuaternion& q)
{
	btVector3 centre(0, 0, 0);
	for (int i = 0; i < psb->m_nodes.size(); i++) centre += psb->m_nodes[i].m_x;
	centre /= btScalar(psb->m_nodes.size());
	psb->translate(-centre);
	psb->rotate(q);
	psb->translate(centre);
}

#endif  // VARIATIONAL_SCENE_HELPERS_H
