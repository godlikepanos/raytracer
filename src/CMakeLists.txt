add_definitions("-Wall -Wpedantic -Wsign-conversion -std=c++14")

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

subdirs(engine)
subdirs(programs/test)