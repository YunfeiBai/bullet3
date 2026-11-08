INCLUDE_DIRECTORIES(
		"${PROJECT_SOURCE_DIR}/src"
		"${PROJECT_SOURCE_DIR}/test/gtest-1.7.0/include")

ADD_DEFINITIONS(-D_VARIADIC_MAX=10)

LINK_LIBRARIES(BulletSoftBody BulletDynamics BulletCollision LinearMath gtest)

IF (NOT WIN32)
	FIND_PACKAGE(Threads)
	LINK_LIBRARIES( ${CMAKE_THREAD_LIBS_INIT} )
ENDIF()

ADD_EXECUTABLE(Test_btDeformableMultiBodyContact test_btDeformableMultiBodyContact.cpp)

ADD_TEST(Test_btDeformableMultiBodyContact_PASS Test_btDeformableMultiBodyContact)

IF (INTERNAL_ADD_POSTFIX_EXECUTABLE_NAMES)
			SET_TARGET_PROPERTIES(Test_btDeformableMultiBodyContact PROPERTIES  DEBUG_POSTFIX "_Debug")
			SET_TARGET_PROPERTIES(Test_btDeformableMultiBodyContact PROPERTIES  MINSIZEREL_POSTFIX "_MinsizeRel")
			SET_TARGET_PROPERTIES(Test_btDeformableMultiBodyContact PROPERTIES  RELWITHDEBINFO_POSTFIX "_RelWithDebugInfo")
ENDIF(INTERNAL_ADD_POSTFIX_EXECUTABLE_NAMES)
