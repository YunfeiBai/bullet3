	project "Test_BulletSoftBody"

	kind "ConsoleApp"

	defines {"USE_GTEST"}

	includedirs
	{
		".",
		"../../src",
		"../gtest-1.7.0/include"
	}

	if os.is("Windows") then
		--see http://stackoverflow.com/questions/12558327/google-test-in-visual-studio-2012
		defines {"_VARIADIC_MAX=10"}
	end

	links {"BulletSoftBody", "BulletDynamics", "BulletCollision", "LinearMath", "gtest"}

	files {
		"test_btVariationalDeformable.cpp",
	}

	if os.is("Linux") then
		links {"pthread"}
	end
