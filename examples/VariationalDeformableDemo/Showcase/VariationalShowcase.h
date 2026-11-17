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
#ifndef VARIATIONAL_SHOWCASE_DEMO_H
#define VARIATIONAL_SHOWCASE_DEMO_H

///Settings of the showcase scene that are chosen from outside (the offline runner sets them from its
///command line; the defaults suit the example browser).
struct VariationalShowcaseOptions
{
	const char* m_dataPath;    // directory with Bullet's data files (boot.vtk, ...); 0: search around the working directory
	const char* m_exportPath;  // existing directory to write scene.json / frames.bin into; 0: no export
	int m_exportEvery;         // write every n-th step (the step is 1/120 s)
	int m_clothResolution;     // nodes along one side of the cloth (lower it for quick previews)
	VariationalShowcaseOptions() : m_dataPath(0), m_exportPath(0), m_exportEvery(2), m_clothResolution(44) {}
};
void VariationalShowcaseSetOptions(const VariationalShowcaseOptions& options);

class CommonExampleInterface* VariationalShowcaseCreateFunc(struct CommonExampleOptions& options);

#endif
