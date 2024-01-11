function(euler_compile_options target)
  if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    target_compile_options(${target} PRIVATE
      -Wall
      -Wextra
      -pedantic
      -Wshadow
      -Wformat=2
      #-Wfloat-equal
      -Wwrite-strings)
  endif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")

  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_compile_options(${target} PRIVATE
      -Wall
      -Wextra
      -pedantic
      -Wshadow
      -Wformat=2
      #-Wfloat-equal
      -Wwrite-strings
      # Triggers too much warning because Kokkos headers are not "system"
      # -Weffc++
      -Wdouble-promotion
      # Triggers too much warning because Kokkos headers are not "system"
      # -Wconversion
      # -Wsign-conversion
      # -Wuseless-cast
      -Wlogical-op)

    if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 6)
      target_compile_options(${target} PRIVATE
        -Wduplicated-cond
        -Wnull-dereference)
    endif(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 6)

    if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 7)
      target_compile_options(${target} PRIVATE
        -Wduplicated-branches)
    endif(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 7)

    # Triggers too much warning because Kokkos headers are not "system"
    # if(NOT KOKKOS_ENABLE_CUDA)
    #   target_compile_options(${target} PRIVATE
    #     -Wold-style-cast)
    # endif(NOT KOKKOS_ENABLE_CUDA)
  endif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")

  if(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    target_compile_options(${target} PRIVATE
      -Wall
      -Wextra
      -pedantic
      -Wshadow
      -Wformat=2
      #-Wfloat-equal
      -Wwrite-strings)
  endif(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
endfunction(euler_compile_options)
