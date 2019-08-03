file(GLOB_RECURSE SOURCES *.c)

add_executable(rt_test ${SOURCES})
include_directories(../..)

target_link_libraries(rt_test rt)
