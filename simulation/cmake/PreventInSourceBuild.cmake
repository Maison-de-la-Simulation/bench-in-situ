string(COMPARE EQUAL ${CMAKE_SOURCE_DIR} ${CMAKE_BINARY_DIR} CMAKE_IN_SOURCE_BUILD)
if(CMAKE_IN_SOURCE_BUILD)
  message(FATAL_ERROR "In-source build not allowed. Please make a new directory (called a build directory) and run CMake from there. You may need to remove CMakeCache.txt and CMakeFiles.")
endif(CMAKE_IN_SOURCE_BUILD)
