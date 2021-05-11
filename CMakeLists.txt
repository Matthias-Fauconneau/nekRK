cmake_minimum_required(VERSION 3.11)
project(nekrk LANGUAGES CXX)
#find_package(OCCA REQUIRED)
include_directories(${OCCA_DIR}/include)
include_directories(${OCCA_DIR}/include)
link_directories(${OCCA_DIR}/lib)
find_library(NAMES occa PATHS ${OCCA_DIR}/lib)

include_directories(${CMAKE_SOURCE_DIR})
include_directories(${CMAKE_BINARY_DIR})
include_directories(mechanisms)

#set(CMAKE_BUILD_TYPE Debug)

add_library(nekrk SHARED
 # List of dependencies :
 src/lib.cpp
)

#find_package(OpenMP)
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

set_property(TARGET nekrk PROPERTY CXX_STANDARD 14)
target_include_directories (nekrk PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries (nekrk PUBLIC occa)
install(TARGETS nekrk DESTINATION lib)
install(FILES README.md DESTINATION .)
install(FILES include/nekrk.h DESTINATION include)
install(DIRECTORY okl DESTINATION .)
install(DIRECTORY share/mechanisms DESTINATION share/)

add_executable(bk1 BK1/bk1.cpp)
target_link_libraries (bk1 LINK_PUBLIC nekrk occa ${MPI_LIBRARIES})
set_target_properties(bk1 PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
set_target_properties(bk1 PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
install(PROGRAMS ${CMAKE_BINARY_DIR}/bk1 DESTINATION bin)

add_executable(bk2 BK2/bk2.cpp)
target_link_libraries (bk2 LINK_PUBLIC nekrk occa ${MPI_LIBRARIES})
set_target_properties(bk2 PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
set_target_properties(bk2 PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
install(PROGRAMS ${CMAKE_BINARY_DIR}/bk2 DESTINATION bin)
