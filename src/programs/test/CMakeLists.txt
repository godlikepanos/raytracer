file(GLOB_RECURSE SOURCES *.cpp)

add_executable(rt_test ${SOURCES})

target_link_libraries(rt_test rt)