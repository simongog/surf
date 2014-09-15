cmake_minimum_required(VERSION 2.8)
cmake_policy(SET CMP0015 NEW)
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules")
include(AppendCompilerFlags)
include(ExternalProject)

project(SURF CXX C)

INCLUDE_DIRECTORIES(${CMAKE_HOME_DIRECTORY}/include
                    ${CMAKE_HOME_DIRECTORY}/external/fastpfor/headers/
                    ${CMAKE_HOME_DIRECTORY}/external/cppzmq/
                    ${CMAKE_HOME_DIRECTORY}/external/zeromq/include/
                    ${CMAKE_BINARY_DIR}/external/sdsl-lite/include
                    ${CMAKE_BINARY_DIR}/external/sdsl-lite/external/libdivsufsort-2.0.1/include
                    )

LINK_DIRECTORIES(${CMAKE_BINARY_DIR}/external/sdsl-lite/lib
				 ${CMAKE_BINARY_DIR}/external/zeromq/lib
					)

append_cxx_compiler_flags("-msse4.2 -std=c++11 -Wall -DNDEBUG" "GCC" CMAKE_CXX_FLAGS)
append_cxx_compiler_flags("-O3 -ffast-math -funroll-loops" "GCC" CMAKE_CXX_FLAGS)
append_cxx_compiler_flags("-msse4.2 -std=c++11 -g -funroll-loops -DNDEBUG -stdlib=libc++" "CLANG" CMAKE_CXX_FLAGS)


ADD_SUBDIRECTORY(external/zeromq)
SET_PROPERTY(DIRECTORY external/zeromq PROPERTY ZMQ_BUILD_TESTS FALSE)

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

    ADD_EXECUTABLE(surf_query-${NAME} src/surf_query.cpp)
    TARGET_LINK_LIBRARIES(surf_query-${NAME} sdsl divsufsort divsufsort64 pthread fastpfor_lib)
	set_property(TARGET surf_query-${NAME} PROPERTY COMPILE_DEFINITIONS IDXNAME="${NAME}" ${compile_defs})
endforeach(f)

ADD_EXECUTABLE(doc_lengths src/doc_lengths.cpp)
TARGET_LINK_LIBRARIES(doc_lengths sdsl)

ADD_EXECUTABLE(gen_patterns src/gen_patterns.cpp)
TARGET_LINK_LIBRARIES(gen_patterns sdsl divsufsort divsufsort64 pthread)

ADD_EXECUTABLE(gen_patterns_int src/gen_patterns_int.cpp)
TARGET_LINK_LIBRARIES(gen_patterns_int sdsl divsufsort divsufsort64 pthread)

ADD_EXECUTABLE(test src/test.cpp)
TARGET_LINK_LIBRARIES(test sdsl divsufsort divsufsort64 pthread)

ADD_EXECUTABLE(create_surf_collection tools/create_surf_collection.cpp)
TARGET_LINK_LIBRARIES(create_surf_collection sdsl divsufsort divsufsort64 pthread fastpfor_lib)

ADD_EXECUTABLE(convert_results_to_trec tools/convert_results_to_trec.cpp)
TARGET_LINK_LIBRARIES(convert_results_to_trec sdsl divsufsort divsufsort64 pthread fastpfor_lib)

ADD_EXECUTABLE(extract_documents tools/extract_documents.cpp)
TARGET_LINK_LIBRARIES(extract_documents sdsl divsufsort divsufsort64 pthread fastpfor_lib)

ADD_EXECUTABLE(extract_document tools/extract_document.cpp)
TARGET_LINK_LIBRARIES(extract_document sdsl divsufsort divsufsort64 pthread fastpfor_lib)

ADD_EXECUTABLE(surf_collection_info tools/surf_collection_info.cpp)
TARGET_LINK_LIBRARIES(surf_collection_info sdsl divsufsort divsufsort64 pthread fastpfor_lib)

ADD_EXECUTABLE(text2surf tools/text2surf.cpp)
TARGET_LINK_LIBRARIES(text2surf sdsl)

ADD_EXECUTABLE(int2surf tools/int2surf.cpp)
TARGET_LINK_LIBRARIES(int2surf sdsl)

ADD_EXECUTABLE(size_info tools/size_info_nn.cpp)
TARGET_LINK_LIBRARIES(size_info)
