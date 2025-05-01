cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

project(RandomGRAPH)

add_executable(RandomGRAPH 
src/main.cpp 
src/GeantGenerator.cpp
src/Graph.cpp
src/GraphUtil.cpp
src/KGraph.cpp 
src/kGraphOperation.cpp
src/test.cpp


headers/base_def.h
headers/GeantGenerator.h
headers/Graph.h
headers/GraphUtil.h
headers/KGraph.h
headers/list.h
headers/test.h

)

file(GLOB TXT_FILES "kGraphFiles/triangle.txt")

source_group("kGraphFiles" FILES ${kGraphFiles})

add_custom_command(
TARGET RandomGRAPH POST_BUILD
COMMAND ${CMAKE_COMMAND} -E copy_directory
${CMAKE_CURRENT_SOURCE_DIR}/kGraphFiles
${CMAKE_CURRENT_BINARY_DIR}/kGraphFiles
)
