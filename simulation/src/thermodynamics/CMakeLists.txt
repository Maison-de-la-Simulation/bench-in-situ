add_library(thermodynamics
  PerfectGas.cpp
  StiffenedGas.cpp
  ThermoParams.cpp)

target_include_directories(thermodynamics PUBLIC ${CMAKE_SOURCE_DIR}/src ${CMAKE_SOURCE_DIR}/src/utils ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_definitions(thermodynamics PUBLIC ${${PROJECT_NAME}_FLOAT_TYPE} ${${PROJECT_NAME}_INT_TYPE} ${${PROJECT_NAME}_SYSTEM_OF_UNITS})
target_link_libraries(thermodynamics PUBLIC inih kokkos)
euler_compile_options(thermodynamics)
