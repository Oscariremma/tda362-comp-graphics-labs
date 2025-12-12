cmake_minimum_required(VERSION 3.5)

project(pathtracer)

find_package(embree 2.12 REQUIRED)

# Workaround for NixOS embree2 package bug where paths are duplicated
# Extract the correct paths from the malformed EMBREE_INCLUDE_DIRS
string(REGEX REPLACE "^.*/nix/store/" "/nix/store/" EMBREE_INCLUDE_DIRS_FIXED "${EMBREE_INCLUDE_DIRS}")
string(REGEX REPLACE "^.*/nix/store/" "/nix/store/" EMBREE_LIBRARY_FIXED "${EMBREE_LIBRARY}")

include_directories(${EMBREE_INCLUDE_DIRS_FIXED})

find_package(OpenMP REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

# Find *all* shaders.
file(GLOB_RECURSE SHADERS
        "${CMAKE_CURRENT_SOURCE_DIR}/*.vert"
        "${CMAKE_CURRENT_SOURCE_DIR}/*.frag"
)
# Separate filter for shaders.
source_group("Shaders" FILES ${SHADERS})

# Build and link executable.
add_executable(${PROJECT_NAME}
        main.cpp
        Pathtracer.h
        Pathtracer.cpp
        sampling.h
        sampling.cpp
        HDRImage.h
        HDRImage.cpp
        embree.h
        embree.cpp
        material.h
        material.cpp
        ${SHADERS}
)

target_link_libraries(${PROJECT_NAME} labhelper ${EMBREE_LIBRARY_FIXED})
config_build_output()
