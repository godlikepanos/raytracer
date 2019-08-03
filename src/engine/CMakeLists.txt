file(GLOB_RECURSE SOURCES *.c)
add_library(rt ${SOURCES})
target_link_libraries(rt m)
