if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "DO NOT BUILD in-tree.")
endif()

cmake_minimum_required ( VERSION 3.11 )
project( orbiparam )

# Add module path
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

find_package(OpenGL REQUIRED)
find_package(GLEW REQUIRED)
find_package(GLFW REQUIRED)
find_package(OpenMesh REQUIRED)

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
                              ${OPENGL_INCLUDE_DIRS}
                              ${GLFW_INCLUDE_DIRS}
                              ${GLEW_INCLUDE_DIRS}
                              ${OPENMESH_INCLUDE_DIRS}
  )
  target_compile_options( orbiparam PUBLIC -O3 -Wall -fpermissive )
  link_directories( ${OPENGL_LIBRARY_DIR}
                    ${GLFW_LIBRARY_DIR}
                    ${GLEW_LIBRARY_DIR}
                    ${OPENMESH_LIBRARY_DIR}
  )
  target_link_libraries( orbiparam
                         PUBLIC
                         ${OPENMESH_LIBRARIES}
                         ${GLFW_LIBRARIES}
                         ${GLEW_LIBRARIES}
  )
endif(UNIX)

# Set StartUp Project (Option)
set_property( DIRECTORY PROPERTY VS_STARTUP_PROJECT "orbiparam" )

