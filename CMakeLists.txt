cmake_minimum_required(VERSION 2.8)
cmake_policy(SET CMP0015 NEW)
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules")
include(AppendCompilerFlags)

project(SURF CXX C)

INCLUDE_DIRECTORIES(${CMAKE_HOME_DIRECTORY}/include
                    ${CMAKE_HOME_DIRECTORY}/external/fastpfor/headers/
                    ${CMAKE_BINARY_DIR}/external/sdsl-lite/include
                    ${CMAKE_BINARY_DIR}/external/sdsl-lite/external/libdivsufsort-2.0.1/include
                    )
LINK_DIRECTORIES(${CMAKE_BINARY_DIR}/external/sdsl-lite/lib)

append_cxx_compiler_flags("-msse4.2 -std=c++11 -Wall -DNDEBUG" "GCC" CMAKE_CXX_FLAGS)
append_cxx_compiler_flags("-O3 -ffast-math -funroll-loops" "GCC" CMAKE_CXX_FLAGS)
append_cxx_compiler_flags("-msse4.2 -std=c++11 -g -funroll-loops -DNDEBUG -stdlib=libc++" "CLANG" CMAKE_CXX_FLAGS)

ADD_SUBDIRECTORY(external/sdsl-lite)

ADD_LIBRARY(fastpfor_lib STATIC external/fastpfor/src/bitpacking.cpp
                                external/fastpfor/src/bitpackingaligned.cpp
                                external/fastpfor/src/bitpackingunaligned.cpp
                                external/fastpfor/src/simdbitpacking.cpp)

# # read the index configs
file(GLOB index_config_files RELATIVE ${CMAKE_HOME_DIRECTORY}/config/ "${CMAKE_HOME_DIRECTORY}/config/*.config")
foreach(f ${index_config_files})
	file(STRINGS ${CMAKE_HOME_DIRECTORY}/config/${f} config_contents)
	set(compile_defs "")
	foreach(keyvalue ${config_contents})
		string(REGEX REPLACE "^[ ]+" "" keyvalue ${keyvalue})
		string(REGEX MATCH "^[^=]+" key ${keyvalue})
		string(REPLACE "${key}=" "" value ${keyvalue})
		set(${key} "${value}")
		list(APPEND compile_defs ${key}=${value})
	endforeach(keyvalue)

	ADD_EXECUTABLE(surf_index-${NAME} src/surf_index.cpp)
	TARGET_LINK_LIBRARIES(surf_index-${NAME} sdsl divsufsort divsufsort64 pthread fastpfor_lib)
	set_property(TARGET surf_index-${NAME} PROPERTY COMPILE_DEFINITIONS IDXNAME="${NAME}" ${compile_defs})

	ADD_EXECUTABLE(surf_search-${NAME} src/surf_search.cpp)
	TARGET_LINK_LIBRARIES(surf_search-${NAME} sdsl divsufsort divsufsort64 pthread fastpfor_lib)
	set_property(TARGET surf_search-${NAME} PROPERTY COMPILE_DEFINITIONS IDXNAME="${NAME}" ${compile_defs})

endforeach(f)

ADD_EXECUTABLE(test src/test.cpp)
TARGET_LINK_LIBRARIES(test sdsl divsufsort divsufsort64 pthread)


