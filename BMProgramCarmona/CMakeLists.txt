cmake_minimum_required(VERSION 3.5.0)
project(Boltzmann_machine VERSION 1.0.0)

# Find the spdlog package
find_package(spdlog REQUIRED)

# Find the Boost package
find_package(Boost REQUIRED COMPONENTS filesystem system)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 11) # Adjust according to your requirements
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
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

# Include directories
# If you're using the single-header file directly, include the directory
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${Boost_INCLUDE_DIRS}) # Include Boost headers

# Add source files
file(GLOB SOURCES src/*.cpp src/*.h)

# Add executable
add_executable(bmc ${SOURCES})

# Link spdlog and Boost to your executable
target_link_libraries(bmc PRIVATE spdlog::spdlog ${Boost_LIBRARIES})

# Set the output directory for the executable
set_target_properties(bmc PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin"
)
