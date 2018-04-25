cmake_minimum_required(VERSION 3.10)
project(FPSSS)

set(CMAKE_CXX_STANDARD 11)


include_directories(${PROJECT_SOURCE_DIR}/minimizer)
add_subdirectory(minimizer)

include_directories(${PROJECT_SOURCE_DIR}/fasta)
add_subdirectory(fasta)

include_directories(${PROJECT_SOURCE_DIR}/lis)
add_subdirectory(lis)


add_executable(FPSSS main.cpp)
target_link_libraries(${PROJECT_NAME} minimizer)
target_link_libraries(${PROJECT_NAME} fasta)
target_link_libraries(${PROJECT_NAME} lis)
