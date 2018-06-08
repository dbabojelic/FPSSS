# add libraries
add_library(opal opal.cpp opal.h ScoreMatrix.cpp ScoreMatrix.hpp)
# c++11 and march
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -march=native")
