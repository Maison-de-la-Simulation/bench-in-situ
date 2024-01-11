cmake_minimum_required(VERSION 3.9)

project(Euler VERSION 1.0 LANGUAGES C CXX)


list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
include(DefaultBuildType)
include(EulerCompileOptions)
include(PreventInSourceBuild)
include(Version)


# Float definitions
set(FLOAT_TYPES FLOAT32 FLOAT64 FLOAT64_EXTENDED)
set(${PROJECT_NAME}_FLOAT_TYPE FLOAT64 CACHE STRING "Float type. Options: ${FLOAT_TYPES}.")
set_property(CACHE ${PROJECT_NAME}_FLOAT_TYPE PROPERTY STRINGS ${FLOAT_TYPES})

# Integer definitions
set(INT_TYPES INT32 INT64)
set(${PROJECT_NAME}_INT_TYPE INT32 CACHE STRING "Int type. Options: ${INT_TYPES}.")
set_property(CACHE ${PROJECT_NAME}_INT_TYPE PROPERTY STRINGS ${INT_TYPES})

# Unit systems definitions
set(SYSTEMS_OF_UNITS CGS_UNITS SI_UNITS)
set(${PROJECT_NAME}_SYSTEM_OF_UNITS SI_UNITS CACHE STRING "System of units. Options: ${SYSTEMS_OF_UNITS}.")
set_property(CACHE ${PROJECT_NAME}_SYSTEM_OF_UNITS PROPERTY STRINGS ${SYSTEMS_OF_UNITS})

option(${PROJECT_NAME}_ENABLE_DOCUMENTATION "Enable ${PROJECT_NAME}'s documentation build." OFF)
option(${PROJECT_NAME}_ENABLE_EXPERIMENTAL "Enable experimental kernels build." OFF)
option(${PROJECT_NAME}_ENABLE_MPI_CUDA_AWARE "Enable MPI CUDA-aware build." OFF)
option(${PROJECT_NAME}_ENABLE_TESTS "Enable ${PROJECT_NAME}'s tests build." OFF)
option(${PROJECT_NAME}_ENABLE_TRILINOS "Enable Trilinos, use Kokkos from Trilinos." OFF)
option(${PROJECT_NAME}_ENABLE_PDI "Enable PDI, use PDI for output." OFF)

if(${PROJECT_NAME}_ENABLE_PDI)
  find_package(PDI REQUIRED)
endif(${PROJECT_NAME}_ENABLE_PDI)

if(${PROJECT_NAME}_ENABLE_TRILINOS)
  # Get Trilinos as one entity
  find_package(Trilinos REQUIRED)

  # Echo trilinos build info just for fun
  message(STATUS "Found Trilinos!  Here are the details: ")
  message(STATUS "   Trilinos_DIR = ${Trilinos_DIR}")
  message(STATUS "   Trilinos_VERSION = ${Trilinos_VERSION}")
  message(STATUS "   Trilinos_PACKAGE_LIST = ${Trilinos_PACKAGE_LIST}")
  message(STATUS "   Trilinos_LIBRARIES = ${Trilinos_LIBRARIES}")
  message(STATUS "   Trilinos_INCLUDE_DIRS = ${Trilinos_INCLUDE_DIRS}")
  message(STATUS "   Trilinos_LIBRARY_DIRS = ${Trilinos_LIBRARY_DIRS}")
  message(STATUS "   Trilinos_TPL_LIST = ${Trilinos_TPL_LIST}")
  message(STATUS "   Trilinos_TPL_INCLUDE_DIRS = ${Trilinos_TPL_INCLUDE_DIRS}")
  message(STATUS "   Trilinos_TPL_LIBRARIES = ${Trilinos_TPL_LIBRARIES}")
  message(STATUS "   Trilinos_TPL_LIBRARY_DIRS = ${Trilinos_TPL_LIBRARY_DIRS}")
  message(STATUS "   Trilinos_BUILD_SHARED_LIBS = ${Trilinos_BUILD_SHARED_LIBS}")
  message(STATUS "   Trilinos_CXX_COMPILER_FLAGS = ${Trilinos_CXX_COMPILER_FLAGS}")
  message(STATUS "   Trilinos_C_COMPILER_FLAGS = ${Trilinos_C_COMPILER_FLAGS}")
  message(STATUS "   Trilinos_Fortran_COMPILER_FLAGS = ${Trilinos_Fortran_COMPILER_FLAGS}")
  message(STATUS "End of Trilinos details\n")

  # Really dirty, this modifies CMAKE_CXX_FLAGS even for projects that don't need it
  # Hope in the near future, Trilinos will define targets that will use transitivity.
  set(CMAKE_CXX_FLAGS "${Trilinos_CXX_COMPILER_FLAGS} ${CMAKE_CXX_FLAGS}")
  set(CMAKE_C_FLAGS "${Trilinos_C_COMPILER_FLAGS} ${CMAKE_C_FLAGS}")
  set(CMAKE_Fortran_FLAGS "${Trilinos_Fortran_COMPILER_FLAGS} ${CMAKE_Fortran_FLAGS}")

  list(FIND Trilinos_PACKAGE_LIST Kokkos Kokkos_List_ID)
  if(Kokkos_List_ID EQUAL -1)
    message(FATAL_ERROR "Kokkos is not enabled inside Trilinos.")
  endif(Kokkos_List_ID EQUAL -1)

  add_library(kokkos INTERFACE)
  target_link_libraries(kokkos INTERFACE ${Trilinos_LIBRARIES})
  set(Kokkos_INCLUDE_DIRS_RET ${Trilinos_INCLUDE_DIRS})

  add_library(Trilinos INTERFACE)
  target_link_libraries(Trilinos INTERFACE ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES})
  target_include_directories(Trilinos SYSTEM INTERFACE ${Trilinos_INCLUDE_DIRS} ${Trilinos_TPL_INCLUDE_DIRS})
else(${PROJECT_NAME}_ENABLE_TRILINOS)
  add_subdirectory(${PROJECT_SOURCE_DIR}/lib/kokkos)
endif(${PROJECT_NAME}_ENABLE_TRILINOS)

add_subdirectory(${PROJECT_SOURCE_DIR}/lib/dms)
add_subdirectory(${PROJECT_SOURCE_DIR}/lib/inih)
add_subdirectory(${PROJECT_SOURCE_DIR}/src)

# Main file
add_executable(main src/main.cpp)
target_link_libraries(main ${PROJECT_NAME})
if(${PROJECT_NAME}_ENABLE_PDI)
  target_compile_definitions(main PUBLIC ${PROJECT_NAME}_ENABLE_PDI)
endif(${PROJECT_NAME}_ENABLE_PDI)
set_target_properties(main PROPERTIES ENABLE_EXPORTS TRUE)

if(${PROJECT_NAME}_ENABLE_DOCUMENTATION)
  find_package(Doxygen REQUIRED)

  doxygen_add_docs(doc
    ${PROJECT_SOURCE_DIR}/lib/dms/src
    ${PROJECT_SOURCE_DIR}/lib/inih
    ${PROJECT_SOURCE_DIR}/lib/kokkos/core/src
    ${PROJECT_SOURCE_DIR}/src)
endif(${PROJECT_NAME}_ENABLE_DOCUMENTATION)

if(${PROJECT_NAME}_ENABLE_TESTS)
  enable_testing()
  add_subdirectory(test)
endif(${PROJECT_NAME}_ENABLE_TESTS)
