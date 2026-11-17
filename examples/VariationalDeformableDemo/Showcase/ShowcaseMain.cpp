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

///Runs the showcase scene (VariationalShowcase.cpp) without graphics and writes its surfaces frame by
///frame, for rendering with render_showcase.py.
///
///   usage: App_VariationalShowcase <output directory> [seconds = 10.5] [data directory] [steps per frame = 2] [cloth resolution = 44]
///
///The output directory must exist. The data directory is the "data" folder of Bullet, which holds
///the tetrahedral meshes of some of the toys; without it the scene runs with fewer toys. The time
///step is 1/120 s, so two steps per frame give 60 frames per second. A lower cloth resolution
///(e.g. 28) makes a quicker preview.

#include <stdio.h>
#include <stdlib.h>
#include "VariationalShowcase.h"
#include "../../CommonInterfaces/CommonExampleInterface.h"
#include "../../CommonInterfaces/CommonGUIHelperInterface.h"
#include "LinearMath/btQuickprof.h"

int main(int argc, char** argv)
{
	if (argc < 2)
	{
		printf("usage: %s <output directory> [seconds = 10.5] [data directory] [steps per frame = 2] [cloth resolution = 44]\n", argv[0]);
		return 1;
	}
	VariationalShowcaseOptions showcase;
	showcase.m_exportPath = argv[1];
	const double seconds = argc > 2 ? atof(argv[2]) : 10.5;
	showcase.m_dataPath = argc > 3 ? argv[3] : 0;
	showcase.m_exportEvery = argc > 4 ? atoi(argv[4]) : 2;
	if (argc > 5) showcase.m_clothResolution = atoi(argv[5]);
	VariationalShowcaseSetOptions(showcase);

	DummyGUIHelper noGraphics;
	CommonExampleOptions options(&noGraphics);
	CommonExampleInterface* example = VariationalShowcaseCreateFunc(options);
	example->initPhysics();
	const int steps = int(seconds * 120.0 + 0.5);
	btClock clock;
	for (int i = 0; i < steps; i++)
	{
		example->stepSimulation(1.f / 120.f);
		if ((i % 60) == 59)
		{
			const double elapsed = clock.getTimeMicroseconds() * 1e-6;
			printf("simulated %.2f of %.2f s in %.0f s (about %.0f s to go)\n", (i + 1) / 120.0, seconds, elapsed, elapsed * (steps - i - 1) / (i + 1));
			fflush(stdout);
		}
	}
	example->exitPhysics();
	delete example;
	return 0;
}
