cmake_minimum_required(VERSION 3.5.0)
project(MyProject VERSION 1.0.0)

# Find the spdlog package
find_package(spdlog REQUIRED)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 14) # Adjust according to your requirements
set(default_build_type "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
    set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
        STRING "Choose the type of build." FORCE)
    # Set the possible values of build type for cmake-gui, ccmake
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
        "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# Add -O3 optimization flag
set(CMAKE_CXX_FLAGS_RELEASE "-O3" "-pthread" "-lgtest" "-lgtest_main")


# Include directories
# If you're using the single-header file directly, include the directory
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# Add source files
file(GLOB SOURCES src/*.cpp src/*.h)

# Add tests files
file(GLOB TESTS tests/*.cpp)

# Add executable run_code
add_executable(bmc ${SOURCES})
# Add executable run_tests
add_executable(runTests ${TESTS})

# Link spdlog to your executable
target_link_libraries(bmc PRIVATE spdlog::spdlog)

#add_executable(runTests tests/test_my_code.cpp src/my_code.cpp)
target_link_libraries(runTests ${GTEST_BOTH_LIBRARIES} pthread my_code)

# Ativar testes no CMake
enable_testing()
add_test(NAME runTests COMMAND runTests)

# Set the output directory for the executable
set_target_properties(bmc PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin"
)

# Set the output directory for the executable
set_target_properties(runTests PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin"
)