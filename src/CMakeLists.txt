add_definitions("-Wall -Wextra -Wpedantic -Wsign-conversion -std=c99")
include_directories(.)
subdirs(engine)
subdirs(programs/test)
