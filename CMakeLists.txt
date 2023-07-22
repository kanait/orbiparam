if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "DO NOT BUILD in-tree.")
endif()

cmake_minimum_required ( VERSION 3.11 )
project( orbiparam )

# Add module path
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

find_package(OpenGL REQUIRED)
find_package(GLEW REQUIRED)
find_package(OpenMesh REQUIRED)
if(UNIX)
  find_package(GLFW REQUIRED)
endif(UNIX)
if(MSVC)
  find_package(glfw3 REQUIRED)
endif(MSVC)

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


if(MSVC)
  add_definitions(-DNOMINMAX)
  target_include_directories( orbiparam
                              PRIVATE
                              ${PROJECT_SOURCE_DIR}/external/render_Eigen
                              ${PROJECT_SOURCE_DIR}/external/eigen
                              ${PROJECT_SOURCE_DIR}/external/stb
                              ${GLEW_INCLUDE_DIR}
                              ${GLFW_INCLUDE_DIR}
                              ${OPENMESH_INCLUDE_DIRS}
                              )

  target_compile_options( orbiparam PRIVATE -O2 )
  target_compile_features( orbiparam PRIVATE cxx_std_14 )

  target_link_libraries( orbiparam
                         PRIVATE
                         GLEW::GLEW
                         glfw
                         ${OPENMESH_LIBRARIES}
                         )
endif(MSVC)

if(UNIX)
  target_include_directories( orbiparam
                              PRIVATE
                              ${PROJECT_SOURCE_DIR}/external/render_Eigen
                              ${PROJECT_SOURCE_DIR}/external/eigen
                              ${PROJECT_SOURCE_DIR}/external/stb
                              ${OPENGL_INCLUDE_DIRS}
                              ${GLFW_INCLUDE_DIRS}
                              ${GLEW_INCLUDE_DIRS}
                              ${OPENMESH_INCLUDE_DIRS}
  )

  target_compile_options( orbiparam PRIVATE -O3 -fpermissive )
  target_compile_features( ${PROJECT_NAME} PRIVATE cxx_std_14 )

  link_directories( ${OPENGL_LIBRARY_DIR}
                    ${GLFW_LIBRARY_DIR}
                    ${GLEW_LIBRARY_DIR}
                    ${OPENMESH_LIBRARY_DIR}
  )

  target_link_libraries( orbiparam
                         PRIVATE
                         ${OPENGL_LIBRARIES}
                         ${OPENMESH_LIBRARIES}
                         GLEW::GLEW
                         ${GLFW_LIBRARY}
  )
endif(UNIX)

# Set StartUp Project (Option)
set_property( DIRECTORY PROPERTY VS_STARTUP_PROJECT "orbiparam" )

