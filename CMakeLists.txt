cmake_minimum_required(VERSION 3.12)
project(nekRK LANGUAGES CXX)

IF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    SET(CMAKE_INSTALL_PREFIX $ENV{HOME}/.local/nekRK CACHE PATH "Default install path" FORCE)
ENDIF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

#set(CMAKE_BUILD_TYPE Debug)

if(${CMAKE_BINARY_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR})
  MESSAGE(FATAL_ERROR "Error:  In-place builds are not supported. Please create a separate build directory")
endif(${CMAKE_BINARY_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR})

#Link:https://gitlab.kitware.com/cmake/community/wikis/doc/cmake/RPATH-handling
set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
#set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_RPATH};${CMAKE_INSTALL_PREFIX}/occa/lib")

add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/3rd_party/occa")

include_directories(mechanisms)

add_library(nekrk SHARED
 # List of dependencies :
 src/lib.cpp
)

find_package(MPI REQUIRED)
if(NOT MPI_FOUND)
  message(FATAL_ERROR "MPI is not found for one of the languages: CXX")
endif()

set (CMAKE_CXX_STANDARD 14)

target_include_directories (nekrk PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries (nekrk PUBLIC libocca)

install(TARGETS nekrk DESTINATION lib)
install(FILES README.md DESTINATION .)
install(FILES include/nekrk.h DESTINATION include)
install(DIRECTORY okl DESTINATION .)
install(DIRECTORY share/mechanisms DESTINATION share/)

add_executable(bk BK/bk.cpp)
target_link_libraries (bk  PRIVATE nekrk libocca)

install(TARGETS bk DESTINATION bin)
