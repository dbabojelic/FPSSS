
# add libraries
add_library(minimizer minimizer.h minimizer.cpp lshbox.h lshbox/lsh/rbslsh.h lshbox/basis.h)
# c++11
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")


