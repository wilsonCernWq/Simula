project (simula-fortran)
cmake_minimum_required (VERSION 2.6)

# use folders
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# global variables
set(SIMULA_BUILD_PATH ${PROJECT_BINARY_DIR})
set(FOX_PATH "utils/FoX-4.1.2")

# compile fox library
add_subdirectory(${FOX_PATH})

# compile main program
add_subdirectory(main)

