add_library(problems
  ImplodeProblem.cpp
  BlastProblem.cpp
  Blast_low_PB_Problem.cpp
  RiemannProblem.cpp
  ProblemFactory.cpp
  Riemann2dProblem.cpp
  RayleighTaylorProblem.cpp
  GreshoProblem.cpp
  OTangProblem.cpp
  SkewedShockProblem.cpp
  KevinHelmoltzMHDProblem.cpp
  ConvectionProblem.cpp)

target_include_directories(problems PUBLIC
  ${CMAKE_SOURCE_DIR}/src
  ${CMAKE_SOURCE_DIR}/src/utils
  ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_definitions(problems PUBLIC
  ${${PROJECT_NAME}_FLOAT_TYPE}
  ${${PROJECT_NAME}_INT_TYPE}
  ${${PROJECT_NAME}_SYSTEM_OF_UNITS})
target_link_libraries(problems PUBLIC inih dms kokkos)
euler_compile_options(problems)

if(${PROJECT_NAME}_ENABLE_MPI_CUDA_AWARE)
  target_compile_definitions(problems PUBLIC ${PROJECT_NAME}_ENABLE_MPI_CUDA_AWARE)
endif(${PROJECT_NAME}_ENABLE_MPI_CUDA_AWARE)
