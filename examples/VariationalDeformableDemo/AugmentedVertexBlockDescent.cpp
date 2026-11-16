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

#include "AugmentedVertexBlockDescent.h"
#include "VariationalSceneHelpers.h"
#include "BulletSoftBody/BulletVariationalDeformable/btVariationalVertexBlockDescent.h"

///"Augmented Vertex Block Descent" (Giles, Diaz, Yuksel 2025): hard constraints inside vertex block
///descent through an augmented Lagrangian, with a dual variable and a growing penalty stiffness per
///constraint.
///
///Two identical chains hang from anchors that swing sideways: four light blocks and a block that is a
///hundred times heavier at the bottom, joined corner to corner. The left chain uses hard attachments
///(solved by AVBD), the right chain plain springs of the same stiffness. The springs stretch under the
///load and the blocks of the right chain come apart, while the hard joints of the left chain stay
///closed. The status line reports the largest joint error of both chains.
static btScalar gAvbdSweeps = 20;

class VariationalAugmentedVertexBlockDescent : public VariationalDemoBase
{
	btVariationalVertexBlockDescent* m_vbd;
	btAlignedObjectArray<int> m_anchorAttachments;
	btAlignedObjectArray<btVector3> m_anchorRest;
	btScalar m_time;
	int m_firstSoftAttachment;

	static int corner(int n, int ix, int iy, int iz) { return (ix * n) + (n + 1) * ((iy * n) + (n + 1) * (iz * n)); }

	void addChain(const btVector3& top, const btVariationalMaterial* material, bool hard)
	{
		const int n = 2;  // cells per block edge
		const btScalar half = 0.15, stiffness = 5e4;
		btSoftBody* previous = 0;
		for (int k = 0; k < 5; k++)
		{
			const bool heavy = k == 4;
			btSoftBody* block = addTetGrid(top - btVector3(0, half + btScalar(2) * half * k, 0), btVector3(half, half, half), n, n, n, heavy ? 30000 : 300, material);
			for (int ix = 0; ix < 2; ix++)
				for (int iz = 0; iz < 2; iz++)
				{
					if (previous)
						m_variationalSolver->addAttachment(previous, corner(n, ix, 0, iz), block, corner(n, ix, 1, iz), stiffness, hard);
					else
					{
						const btVector3 target = block->m_nodes[corner(n, ix, 1, iz)].m_x;
						m_anchorAttachments.push_back(m_variationalSolver->addAttachment(block, corner(n, ix, 1, iz), target, stiffness, hard));
						m_anchorRest.push_back(target);
					}
				}
			previous = block;
		}
	}

public:
	VariationalAugmentedVertexBlockDescent(struct GUIHelperInterface* helper)
		: VariationalDemoBase(helper), m_vbd(0), m_time(0), m_firstSoftAttachment(0)
	{
	}

	void resetCamera()
	{
		float dist = 3.4;
		float pitch = -8;
		float yaw = 0;
		float targetPos[3] = {0, 2.3, 0};
		m_guiHelper->resetCamera(dist, yaw, pitch, targetPos[0], targetPos[1], targetPos[2]);
	}

	virtual const char* methodName() const { return "Augmented VBD: hard joints (left) vs springs (right)"; }

	virtual void updateStatus()
	{
		VariationalDemoBase::updateStatus();
		double hardError = 0, softError = 0;
		const btVariationalSystem& sys = system();
		for (int a = 0; a < sys.m_attachments.size(); a++)
		{
			const btVariationalAttachment& at = sys.m_attachments[a];
			const double e = (sys.m_x[at.m_nodeA] - (at.m_nodeB >= 0 ? sys.m_x[at.m_nodeB] : at.m_target)).length();
			if (at.m_hard)
				hardError = btMax(hardError, e);
			else
				softError = btMax(softError, e);
		}
		const size_t len = strlen(m_statusText);
		snprintf(m_statusText + len, sizeof(m_statusText) - len, " | joint error: hard %.2f mm, springs %.1f mm", hardError * 1000.0, softError * 1000.0);
		if (m_logToConsole && (m_stepCount % 30) == 0) printf("      joint error: hard %.3e m, springs %.3e m\n", hardError, softError);
	}

	void initPhysics()
	{
		createWorld(btVector3(0, -9.8, 0));
		addGround();
		m_vbd = new btVariationalVertexBlockDescent();
		m_vbd->m_rayleighDamping = 2e-3;
		m_integrators.push_back(m_vbd);
		m_variationalSolver->setIntegrator(m_vbd);
		btVariationalStableNeoHookean* rubber = new btVariationalStableNeoHookean(2e6, 0.4);
		m_materials.push_back(rubber);
		// with the default camera +x is on the left of the screen
		addChain(btVector3(0.9, 3.2, 0), rubber, true);
		m_firstSoftAttachment = m_anchorAttachments.size();
		addChain(btVector3(-0.9, 3.2, 0), rubber, false);
		m_guiHelper->autogenerateGraphicsObjects(m_dynamicsWorld);

		if (m_guiHelper->getParameterInterface())
		{
			SliderParams sweeps("sweeps per step", &gAvbdSweeps);
			sweeps.m_minVal = 2;
			sweeps.m_maxVal = 60;
			sweeps.m_clampToIntegers = true;
			m_guiHelper->getParameterInterface()->registerSliderFloatParameter(sweeps);
		}
	}

	virtual void stepSimulation(float deltaTime)
	{
		m_vbd->m_iterations = int(gAvbdSweeps + btScalar(0.5));
		// the anchors swing sideways after a second of settling
		m_time += m_fixedTimeStep;
		const btScalar offset = m_time > btScalar(1) ? btScalar(0.35) * btSin(SIMD_2_PI * (m_time - btScalar(1)) / btScalar(3)) : btScalar(0);
		for (int i = 0; i < m_anchorAttachments.size(); i++)
			m_variationalSolver->setAttachmentTarget(m_anchorAttachments[i], m_anchorRest[i] + btVector3(offset, 0, 0));
		VariationalDemoBase::stepSimulation(deltaTime);
	}
};

class CommonExampleInterface* VariationalAugmentedVertexBlockDescentCreateFunc(struct CommonExampleOptions& options)
{
	return new VariationalAugmentedVertexBlockDescent(options.m_guiHelper);
}
