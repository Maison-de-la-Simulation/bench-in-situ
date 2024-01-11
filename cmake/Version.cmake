find_package(Git QUIET)
if(NOT GIT_FOUND)
  set(GIT_BUILD_STRING "N/A")
  set(GIT_BRANCH "N/A")
else(NOT GIT_FOUND)
  execute_process(COMMAND ${GIT_EXECUTABLE} describe --tags --always --dirty
    OUTPUT_VARIABLE GIT_BUILD_STRING
    OUTPUT_STRIP_TRAILING_WHITESPACE)
  execute_process(COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_BRANCH
    OUTPUT_STRIP_TRAILING_WHITESPACE)
endif(NOT GIT_FOUND)

execute_process(COMMAND date "+%d/%m/%y"
  OUTPUT_VARIABLE DATE_STRING
  OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND date "+%H:%M:%S"
  OUTPUT_VARIABLE TIME_STRING
  OUTPUT_STRIP_TRAILING_WHITESPACE)

configure_file(
  ${PROJECT_SOURCE_DIR}/src/HydroVersion.cpp.in
  ${PROJECT_BINARY_DIR}/src/HydroVersion.cpp
  @ONLY)

add_library(version ${PROJECT_BINARY_DIR}/src/HydroVersion.cpp)
target_include_directories(version PUBLIC ${PROJECT_SOURCE_DIR}/src)
target_compile_options(version PRIVATE -std=c++11)
