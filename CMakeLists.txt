cmake_minimum_required(VERSION 3.11)
project(nekrk LANGUAGES CXX)

set(CMAKE_C_FLAGS "-g")
set(CMAKE_CXX_FLAGS "-g")
#set(CMAKE_C_FLAGS "-g -fstack-protector -fstack-protector-all -fsanitize=bounds -Werror -fsanitize=address")
#set(CMAKE_CXX_FLAGS "-g -fstack-protector -fstack-protector-all -fsanitize=bounds -Werror -fsanitize=address")

#find_package(OCCA REQUIRED)
include_directories(${OCCA_DIR}/include)
link_directories(${OCCA_DIR}/lib)
#message(${OCCA_DIR}/lib)
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

set_property(TARGET nekrk PROPERTY CXX_STANDARD 20) # Will automatically decay to the latest supported standard
target_include_directories (nekrk PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries (nekrk PUBLIC occa)
install(TARGETS nekrk DESTINATION lib)
install(FILES README.md DESTINATION share)
install(FILES include/nekrk.h DESTINATION include)
install(DIRECTORY okl DESTINATION .)
install(DIRECTORY mechanisms DESTINATION .)

add_executable(main main.cc)
target_link_libraries (main LINK_PUBLIC nekrk occa ${MPI_LIBRARIES})
set_target_properties(main PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
set_target_properties(main PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
install(PROGRAMS ${CMAKE_BINARY_DIR}/main DESTINATION bin)
