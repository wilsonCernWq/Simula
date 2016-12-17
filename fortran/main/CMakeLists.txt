# CMake project file
project (main)
cmake_minimum_required (VERSION 2.6)
enable_language (Fortran)

# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
    "Choose the type of build, options are: None Debug Release."
    FORCE)
endif (NOT CMAKE_BUILD_TYPE)

# default installation
get_filename_component (default_prefix "." ABSOLUTE)
set (CMAKE_INSTALL_PREFIX ${default_prefix} 
  CACHE 
  STRING "Choose the installation directory as"
  FORCE)
message(STATUS "installation firectory: ${default_prefix}")

# FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  # gfortran
  set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3 -cpp")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g -fbounds-check -cpp")
elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
  # ifort (untested)
  set (CMAKE_Fortran_FLAGS_RELEASE "-f77rtl -O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-f77rtl -Od")
elseif (Fortran_COMPILER_NAME MATCHES "g77")
  # g77
  set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3 -m32")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g -m32")
else (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
  message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
  message ("No optimized Fortran compiler flags are known, we just try -O2...")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
endif (Fortran_COMPILER_NAME MATCHES "gfortran.*")

# Add source files
file(GLOB_RECURSE SRCS "src/*.f90")
file(GLOB_RECURSE INDS "include/*")
foreach(f ${SRCS} ${INDS})
  message ("-- " ${f})
endforeach()

# include files
include_directories(
  "include"
  ${SIMULA_BUILD_PATH}/${FOX_PATH}/dom
  ${SIMULA_BUILD_PATH}/${FOX_PATH}/wxml
  ${SIMULA_BUILD_PATH}/${FOX_PATH}/wcml
  ${SIMULA_BUILD_PATH}/${FOX_PATH}/wkml
  ${SIMULA_BUILD_PATH}/${FOX_PATH}/common
  ${SIMULA_BUILD_PATH}/${FOX_PATH}/utils
  ${SIMULA_BUILD_PATH}/${FOX_PATH}/fsys
  ${SIMULA_BUILD_PATH}/${FOX_PATH}/sax
)

# executable
set(EXECUTABLE_OUTPUT_PATH ${SIMULA_BUILD_PATH})
add_executable ("simula.exe" ${SRCS} ${INDS})

# organize folder
macro(GroupSources curdir)
   file(GLOB children RELATIVE ${PROJECT_SOURCE_DIR}/${curdir} ${PROJECT_SOURCE_DIR}/${curdir}/*)
   foreach(child ${children})
          if(IS_DIRECTORY ${PROJECT_SOURCE_DIR}/${curdir}/${child})
                  GroupSources(${curdir}/${child})
          else()
                  string(REPLACE "/" "\\" groupname ${curdir})
                  string(REPLACE "src" "Source" groupname ${groupname})
            source_group(${groupname} FILES ${PROJECT_SOURCE_DIR}/${curdir}/${child})
          endif()
   endforeach()
endmacro()

GroupSources("include")
GroupSources("src")