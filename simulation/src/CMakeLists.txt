add_subdirectory(problems)
add_subdirectory(thermodynamics)

# Source files of euler project
add_library(${PROJECT_NAME}
  godunov_solver/FluxesAndUpdateKernelDispatch.cpp
  godunov_solver/Godunov.cpp
  godunov_solver/FluxesAndUpdateExecution.cpp
  godunov_solver/MusclReconstructionExecution.cpp
  Hydro.cpp
  HydroEngine.cpp
  HydroParams.cpp
  io/Reader.cpp
  io/Writer.cpp
  io/WriterVTK.cpp
  shared/ConservativeToPrimitiveExecution.cpp
  shared/PrimitiveToConservativeExecution.cpp
  shared/global_meanExecution.cpp
  shared/vp2Execution.cpp
  shared/TimeStepExecution.cpp
  shared/ConvectionSourceTermExecution.cpp
  shared/MagneticResistivityExecution.cpp
  SolverFactory.cpp)
target_include_directories(${PROJECT_NAME} PUBLIC
  ${CMAKE_SOURCE_DIR}/src
  ${CMAKE_SOURCE_DIR}/src/utils
  ${CMAKE_CURRENT_SOURCE_DIR}/shared)
target_link_libraries(${PROJECT_NAME} PUBLIC inih kokkos dms problems thermodynamics version)
target_compile_definitions(${PROJECT_NAME} PUBLIC
  ${${PROJECT_NAME}_FLOAT_TYPE}
  ${${PROJECT_NAME}_INT_TYPE}
  ${${PROJECT_NAME}_SYSTEM_OF_UNITS})
euler_compile_options(${PROJECT_NAME})
set_target_properties(${PROJECT_NAME} PROPERTIES ENABLE_EXPORTS TRUE)

# Temporary fix for Kokkos v3.4
if (Kokkos_ENABLE_AGGRESSIVE_VECTORIZATION)
  target_compile_definitions(${PROJECT_NAME} PUBLIC KOKKOS_ENABLE_AGGRESSIVE_VECTORIZATION)
endif()

if(${PROJECT_NAME}_ENABLE_MPI_CUDA_AWARE)
  target_compile_definitions(${PROJECT_NAME} PUBLIC ${PROJECT_NAME}_ENABLE_MPI_CUDA_AWARE)
endif(${PROJECT_NAME}_ENABLE_MPI_CUDA_AWARE)

if(${PROJECT_NAME}_ENABLE_PDI)
  target_sources(${PROJECT_NAME} PRIVATE io/ReaderPDI.cpp io/WriterPDI.cpp)
  target_link_libraries(${PROJECT_NAME} PUBLIC PDI::PDI_C)
  target_compile_definitions(${PROJECT_NAME} PUBLIC ${PROJECT_NAME}_ENABLE_PDI)
endif(${PROJECT_NAME}_ENABLE_PDI)
