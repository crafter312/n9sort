cmake_minimum_required(VERSION 3.14)
project(EnergyLossCalculator LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Define variables for CMake build
set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(SRC_PARENT ${CMAKE_CURRENT_SOURCE_DIR}/../../src)

# List of source files
set(SOURCES main.cpp)
list(TRANSFORM SOURCES PREPEND ${SRC}/)

# Download and incorporate CATima automatically if not found locally
include(FetchContent)
FetchContent_Declare(
    catima
    GIT_REPOSITORY https://github.com/hrosiak/catima.git
    GIT_TAG        v1.90
)
FetchContent_MakeAvailable(catima)

# Define the Boost source
# We use the official Boost GitHub mirror to grab just the preprocessor module
FetchContent_Declare(
    boost_preprocessor
    GIT_REPOSITORY https://github.com/boostorg/preprocessor.git
    GIT_TAG        master # Or a specific tag like boost-1.84.0
)
FetchContent_MakeAvailable(boost_preprocessor)
set(BOOST_PP_INCLUDE "${boost_preprocessor_SOURCE_DIR}/include")

add_executable(eloss ${SOURCES})
target_link_libraries(eloss PRIVATE catima boost_preprocessor)
target_include_directories(eloss PRIVATE ${SRC_PARENT})
