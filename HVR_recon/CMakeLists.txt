cmake_minimum_required(VERSION 3.13)

Project(ReconstructCDRs CXX)

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_PREFIX_PATH seqan/util/cmake)

set(SEQAN_INCLUDE_PATH seqan/include)


find_package(ZLIB)

find_package(SeqAn REQUIRED)

include_directories (${SEQAN_INCLUDE_DIRS})
add_definitions (${SEQAN_DEFINITIONS})

add_executable (ReconstructCDRs
				cdr_reconstruction.cpp 
				input_output_utilities.h 
				input_output_utilities.cpp
				alignment_utilities.h
				alignment_utilities.cpp
				declarations.h 
				parser.h)

target_link_libraries (ReconstructCDRs ${SEQAN_LIBRARIES})

