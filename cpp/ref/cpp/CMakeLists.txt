cmake_minimum_required(VERSION 3.1.0)
project(Simula)

# check c++11 standard
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
  message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support")
endif()

# include path
include_directories("include")
include_directories("utils/rapidxml_1_13")
include_directories("utils/boost_1_60_0")

# source files
file(GLOB SOURCE_FILES "src/*.cpp")

# set executable
add_executable(simula ${SOURCE_FILES})

# copy additional files
file(COPY "convert.m" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})