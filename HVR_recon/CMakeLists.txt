# Minimum cmake version
cmake_minimum_required (VERSION 3.0.0)

# Name of project and that it is C++ only.
project (SeqanReconstruct CXX)

# Search for zlib as a dependency for SeqAn.
find_package (ZLIB)

# Load the SeqAn module and fail if not found.
find_package (SeqAn REQUIRED)

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# Add executable and link against SeqAn dependencies.
add_executable (SeqanReconstruct hypervariable_reconstruction.cpp)
target_link_libraries (SeqanReconstruct ${SEQAN_LIBRARIES})

#Call As Follows:
#cmake ../ -DCMAKE_PREFIX_PATH="$HOME/devel/seqan/util/cmake" -DSEQAN_INCLUDE_PATH="$HOME/devel/seqan/include"
