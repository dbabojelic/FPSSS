cmake_minimum_required(VERSION 3.5)
project(FPSSS)

set(CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -march=native")


include_directories(${PROJECT_SOURCE_DIR}/minimizer)
add_subdirectory(minimizer)

include_directories(${PROJECT_SOURCE_DIR}/fasta)
add_subdirectory(fasta)

include_directories(${PROJECT_SOURCE_DIR}/filter)
add_subdirectory(filter)

include_directories(${PROJECT_SOURCE_DIR}/opal)
add_subdirectory(opal)


add_executable(FPSSS main.cpp)
target_link_libraries(${PROJECT_NAME} minimizer)
target_link_libraries(${PROJECT_NAME} fasta)
target_link_libraries(${PROJECT_NAME} filter)
target_link_libraries(${PROJECT_NAME} opal)
