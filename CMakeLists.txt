#cmake -DCMAKE_PREFIX_PATH=~/libtorch ..
#cmake --build . --config Release

cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(MPF)

find_package(Torch REQUIRED)

find_package(MPI REQUIRED)

#target_link_libraries(${PROJECT_NAME} ${MPI_LIBRARIES})


include_directories(SYSTEM ${MPI_INCLUDE_PATH})
include_directories(./include)

file(GLOB SOURCES "./src/*.cpp")

add_executable(MPF ${SOURCES})

target_link_libraries(MPF ${TORCH_LIBRARIES} ${MPI_LIBRARIES})
set_property(TARGET MPF PROPERTY CXX_STANDARD 14)