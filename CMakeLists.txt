cmake_minimum_required(VERSION 3.16 FATAL_ERROR)

include(FetchContent)
include(ExternalProject)

project(n9sort)

# Enable C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

message("C++ compiler: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")

# Setup directories.
set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_INSTALL_LIBDIR ${PROJECT_BINARY_DIR}/lib)

# Uncomment the below compile time variable to enable debug outputs
#add_definitions(-DENABLE_DEBUG=true)

# Enable relativistic calculations
add_definitions(-Drel=true)

# Set project sources
set(SOURCES SortConfig.cpp Det.cpp Gobbi28.cpp wood.cpp histo.cpp HINP.cpp telescope.cpp elist.cpp solution.cpp pid.cpp ZApar.cpp einstein.cpp losses.cpp loss3.cpp correl2.cpp parType.cpp calibrate.cpp Input.cpp)
set(LIBHEADERS parType.h wood.h einstein.h solution.h)

list(TRANSFORM SOURCES PREPEND ${SRC}/)

# Define the Boost source
# We use the official Boost GitHub mirror to grab just the preprocessor module
#FetchContent_Declare(
#    boost_preprocessor
#    GIT_REPOSITORY https://github.com/boostorg/preprocessor.git
#    GIT_TAG        master # Or a specific tag like boost-1.84.0
#)
#FetchContent_MakeAvailable(boost_preprocessor)
#set(BOOST_PP_INCLUDE "${boost_preprocessor_SOURCE_DIR}/include")

# Locate the ROOT package and define a number of useful targets and variables
find_package(ROOT REQUIRED COMPONENTS RIO Tree Hist TreePlayer Core Imt Thread MultiProc)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(${ROOT_INCLUDE_DIRS})
message("ROOT include directory: ${ROOT_INCLUDE_DIRS}")
message("Available ROOT library components: ${ROOT_LIBRARIES}")

# Generate the dictionary using ROOT_GENERATE_DICTIONARY
#include_directories(${BOOST_PP_INCLUDE})
include_directories(${SRC})
ROOT_GENERATE_DICTIONARY(G__n9sort ${LIBHEADERS} MODULE n9sort LINKDEF ${SRC}/LinkDef.h)

# Create a custom target for the dictionary generation
add_custom_target(generate_dict ALL DEPENDS ${CMAKE_BINARY_DIR}/G__n9sort.cxx)

# Create a shared library to handle the generated sources (with sim class)
add_library(n9sort SHARED ${SOURCES} ${CMAKE_BINARY_DIR}/G__n9sort.cxx)
set_target_properties(n9sort PROPERTIES EXCLUDE_FROM_ALL TRUE)
#target_link_libraries(n9sort PRIVATE boost_preprocessor)
target_link_libraries(n9sort PUBLIC ROOT::RIO ROOT::Tree ROOT::Hist ROOT::TreePlayer ROOT::Core ROOT::Imt ROOT::Thread ROOT::MultiProc)
set_target_properties(n9sort PROPERTIES
	BUILD_RPATH "${CMAKE_BINARY_DIR}"
	INSTALL_RPATH "${CMAKE_BINARY_DIR}"
)

# Create sort executable
add_executable(sort ${SRC}/sort.cpp)
#target_link_libraries(sort PRIVATE boost_preprocessor)
target_link_libraries(sort PUBLIC n9sort ROOT::RIO ROOT::Tree ROOT::Hist ROOT::TreePlayer ROOT::Core ROOT::Imt ROOT::Thread ROOT::MultiProc)
set_target_properties(sort PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set_target_properties(sort PROPERTIES
	BUILD_RPATH "${CMAKE_BINARY_DIR}"
	INSTALL_RPATH "${CMAKE_BINARY_DIR}"
)



