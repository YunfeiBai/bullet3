# Offline runner of the variational deformable showcase scene: it simulates without graphics and
# writes the surfaces frame by frame for rendering (see README.md in this folder).

INCLUDE_DIRECTORIES(
	${BULLET_PHYSICS_SOURCE_DIR}/src
	${BULLET_PHYSICS_SOURCE_DIR}/examples
)

LINK_LIBRARIES(
	BulletSoftBody BulletDynamics BulletCollision Bullet3Common LinearMath
)

IF (WIN32)
	ADD_EXECUTABLE(App_VariationalShowcase
		ShowcaseMain.cpp
		VariationalShowcase.cpp
		VariationalShowcase.h
		VariationalMeshExport.h
		${BULLET_PHYSICS_SOURCE_DIR}/build3/bullet.rc
	)
ELSE()
	ADD_EXECUTABLE(App_VariationalShowcase
		ShowcaseMain.cpp
		VariationalShowcase.cpp
		VariationalShowcase.h
		VariationalMeshExport.h
	)
ENDIF()

IF (INTERNAL_ADD_POSTFIX_EXECUTABLE_NAMES)
			SET_TARGET_PROPERTIES(App_VariationalShowcase PROPERTIES  DEBUG_POSTFIX "_Debug")
			SET_TARGET_PROPERTIES(App_VariationalShowcase PROPERTIES  MINSIZEREL_POSTFIX "_MinsizeRel")
			SET_TARGET_PROPERTIES(App_VariationalShowcase PROPERTIES  RELWITHDEBINFO_POSTFIX "_RelWithDebugInfo")
ENDIF(INTERNAL_ADD_POSTFIX_EXECUTABLE_NAMES)
