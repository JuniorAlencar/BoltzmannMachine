cmake_minimum_required(VERSION 3.5.0)
project(MyProject VERSION 1.0.0)

# Find the spdlog package
find_package(spdlog REQUIRED)

# Find the Google Test package
find_package(GTest REQUIRED)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 14)
set(default_build_type "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
    set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
        STRING "Choose the type of build." FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
        "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# Add -O3 optimization flag
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

# Include directories
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)

# Add source files (excluding main.cpp)
file(GLOB_RECURSE SOURCES "src/*.cpp" "src/*.h")
list(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp")

# Add executable for the main application
add_executable(bmc src/main.cpp ${SOURCES})
target_link_libraries(bmc PRIVATE spdlog::spdlog)

# Set the output directory for the main executable
set_target_properties(bmc PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin"
)

# Find all test files
file(GLOB TEST_FILES "tests/*.cpp")

# Create a separate executable for each test file
foreach(TEST_FILE ${TEST_FILES})
    # Get the filename without extension
    get_filename_component(TEST_NAME ${TEST_FILE} NAME_WE)
    
    # Add executable for each test
    add_executable(${TEST_NAME} ${TEST_FILE} ${SOURCES})

    # Link the necessary libraries
    target_link_libraries(${TEST_NAME} PRIVATE spdlog::spdlog ${GTEST_LIBRARIES} pthread)

    # Set the output directory for the test executables
    set_target_properties(${TEST_NAME} PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin"
    )

    # Register the test
    add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
endforeach()

# Enable testing
enable_testing()


