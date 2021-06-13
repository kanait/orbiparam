if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "DO NOT BUILD in-tree.")
endif()

cmake_minimum_required ( VERSION 3.10 )
project( orbiparam )

SET(CMAKE_C_COMPILER_WORKS 1 CACHE INTERNAL "")
SET(CMAKE_CXX_COMPILER_WORKS 1 CACHE INTERNAL "")

add_executable( orbiparam
                MyMesh.hxx
                Orbifold.hxx
                orbiparam.cxx
                Param.hxx
                ShortestPathDijkstra.hxx
                ToMeshR.hxx
                VertexPQ.hxx
)
find_package(OpenGL REQUIRED)

target_include_directories( orbiparam
                            PUBLIC
                            ${PROJECT_SOURCE_DIR}/external/render_Eigen
                            ${PROJECT_SOURCE_DIR}/external/eigen
                            ${PROJECT_SOURCE_DIR}/external/stb
)

target_compile_options( orbiparam PUBLIC -O3 -Wall)
target_compile_features( orbiparam PUBLIC cxx_std_11 )

target_link_libraries( orbiparam
                       LINK_PUBLIC
)

if(UNIX)
  target_include_directories( orbiparam
                              PUBLIC
                              /usr/local/include
  )
  target_compile_options( orbiparam PUBLIC -O3 -Wall -fpermissive )
  link_directories( /usr/lib/x86_64-linux-gnu
                    /usr/local/lib
  )
  target_link_libraries( orbiparam
                         PUBLIC
			 -lOpenGL -lGLU
			 -lOpenMeshCore
			 -lglfw -lGLEW
  )
endif(UNIX)

# Set StartUp Project (Option)
set_property( DIRECTORY PROPERTY VS_STARTUP_PROJECT "orbiparam" )

