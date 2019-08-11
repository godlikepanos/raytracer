cmake_minimum_required(VERSION 2.8)
project(rt)

macro(install_executable exe)
	add_custom_command(TARGET ${exe} POST_BUILD
		COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/bin
		COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:${exe}> ${CMAKE_BINARY_DIR}/bin)
endmacro()

include_directories(.)
add_subdirectory(src)
