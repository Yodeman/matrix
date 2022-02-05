add_executable(test mtest2.cpp)

target_include_directories(
	test PUBLIC
	"${PROJECT_BINARY_DIR}"
	)
