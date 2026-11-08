#include <BulletSoftBody/btSoftBody.h>
#include <gtest/gtest.h>

static void expectMass(btScalar expected, btScalar actual)
{
	EXPECT_NEAR(expected, actual, btScalar(1e-5) * expected);
}

GTEST_TEST(BulletSoftBody, SetTotalMassFromTetrasLumpsRestVolume)
{
	// two tetrahedra sharing the face (0,1,2) with volumes 1/6 and 3/6, the second one inverted
	btVector3 x[5] = {btVector3(0, 0, 0), btVector3(1, 0, 0), btVector3(0, 1, 0), btVector3(0, 0, 1), btVector3(0, 0, -3)};
	btSoftBodyWorldInfo worldInfo;
	btSoftBody psb(&worldInfo, 5, x, 0);
	psb.appendTetra(0, 1, 2, 3);
	psb.appendTetra(0, 2, 1, 4);
	psb.setTotalMassFromTetras(8);

	// a quarter of each volume per node, in units of 1/24: nodes 0-2 get 1+3, node 3 gets 1, node 4 gets 3
	expectMass(2, psb.getMass(0));
	expectMass(2, psb.getMass(1));
	expectMass(2, psb.getMass(2));
	expectMass(0.5, psb.getMass(3));
	expectMass(1.5, psb.getMass(4));
	expectMass(8, psb.getTotalMass());
}

GTEST_TEST(BulletSoftBody, SetTotalMassFromTetrasKeepsUnusedNodesDynamic)
{
	btVector3 x[5] = {btVector3(0, 0, 0), btVector3(1, 0, 0), btVector3(0, 1, 0), btVector3(0, 0, 1), btVector3(5, 5, 5)};
	btSoftBodyWorldInfo worldInfo;
	btSoftBody psb(&worldInfo, 5, x, 0);
	psb.appendTetra(0, 1, 2, 3);
	psb.setTotalMassFromTetras(1);

	// node 4 is in no tetrahedron: it gets the smallest lumped mass instead of becoming static
	for (int i = 0; i < 5; ++i)
	{
		EXPECT_GT(psb.m_nodes[i].m_im, 0);
		expectMass(0.2, psb.getMass(i));
	}
}

GTEST_TEST(BulletSoftBody, SetTotalMassFromTetrasWithoutTetrasIsSetTotalMass)
{
	btVector3 x[3] = {btVector3(0, 0, 0), btVector3(1, 0, 0), btVector3(0, 1, 0)};
	btSoftBodyWorldInfo worldInfo;
	btSoftBody cloth(&worldInfo, 3, x, 0);
	btSoftBody reference(&worldInfo, 3, x, 0);
	cloth.appendFace(0, 1, 2);
	reference.appendFace(0, 1, 2);
	cloth.setMass(1, 3);
	reference.setMass(1, 3);
	cloth.setTotalMassFromTetras(10);
	reference.setTotalMass(10);

	for (int i = 0; i < 3; ++i)
	{
		EXPECT_EQ(reference.getMass(i), cloth.getMass(i));
	}
}

int main(int argc, char** argv)
{
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
