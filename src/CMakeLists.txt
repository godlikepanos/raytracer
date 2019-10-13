add_definitions("-Wall -Wextra -Wpedantic -Wsign-conversion -std=c11 -mrdrnd -mrdseed")
include_directories(.)
subdirs(engine)
subdirs(programs/test)
