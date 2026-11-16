
INCLUDE_DIRECTORIES(
		"${PROJECT_SOURCE_DIR}/src"
		"${PROJECT_SOURCE_DIR}/test/gtest-1.7.0/include")

ADD_DEFINITIONS(-DUSE_GTEST)
ADD_DEFINITIONS(-D_VARIADIC_MAX=10)

LINK_LIBRARIES(BulletSoftBody BulletDynamics BulletCollision LinearMath gtest)

IF (NOT WIN32)
	FIND_PACKAGE(Threads)
	LINK_LIBRARIES( ${CMAKE_THREAD_LIBS_INIT} )
ENDIF()

ADD_EXECUTABLE(Test_btVariationalDeformable test_btVariationalDeformable.cpp)

ADD_TEST(Test_btVariationalDeformable_PASS Test_btVariationalDeformable)

# timings of the integrators on three reference scenes (not a test, see the module README)
ADD_EXECUTABLE(Benchmark_btVariationalDeformable benchmark_btVariationalDeformable.cpp)

IF (INTERNAL_ADD_POSTFIX_EXECUTABLE_NAMES)
			SET_TARGET_PROPERTIES(Test_btVariationalDeformable PROPERTIES  DEBUG_POSTFIX "_Debug")
			SET_TARGET_PROPERTIES(Test_btVariationalDeformable PROPERTIES  MINSIZEREL_POSTFIX "_MinsizeRel")
			SET_TARGET_PROPERTIES(Test_btVariationalDeformable PROPERTIES  RELWITHDEBINFO_POSTFIX "_RelWithDebugInfo")
			SET_TARGET_PROPERTIES(Benchmark_btVariationalDeformable PROPERTIES  DEBUG_POSTFIX "_Debug")
			SET_TARGET_PROPERTIES(Benchmark_btVariationalDeformable PROPERTIES  MINSIZEREL_POSTFIX "_MinsizeRel")
			SET_TARGET_PROPERTIES(Benchmark_btVariationalDeformable PROPERTIES  RELWITHDEBINFO_POSTFIX "_RelWithDebugInfo")
ENDIF(INTERNAL_ADD_POSTFIX_EXECUTABLE_NAMES)
