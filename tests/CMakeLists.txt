# CMakeLists.txt for test suite
cmake_minimum_required(VERSION 3.15)
project(NuclearPhysicsTests Fortran)

# Set Fortran standard
set(CMAKE_Fortran_FLAGS "-g -O0 -fcheck=all -fbacktrace -Wall -Wextra -fsanitize=address -fno-omit-frame-pointer")

# Set module directory
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
include_directories(${CMAKE_Fortran_MODULE_DIRECTORY})

# Find BLAS and LAPACK
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

# Core library
add_library(core
    ../core/kinds.f90
    ../core/constants.f90
    ../core/numeric_utils.f90
    ../core/rng.f90
)

# Linear algebra library (depends on BLAS/LAPACK)
add_library(linalg
    ../kernels/linear_algebra/dense_matrix.f90
    ../kernels/linear_algebra/eigen.f90
    ../kernels/linear_algebra/solve_linear.f90
    ../kernels/linear_algebra/sparse_matrix.f90
)
target_link_libraries(linalg core ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

# ODE solvers
add_library(ode
    ../kernels/ode/rk4.f90
    ../kernels/ode/dormand_prince.f90
    ../kernels/ode/backward_euler.f90
)
target_link_libraries(ode core linalg)

# Optimisation
add_library(optimisation
    ../kernels/optimisation/gradient_descent.f90
    ../kernels/optimisation/conjugate_newton.f90
    ../kernels/optimisation/quasi_newton.f90
    ../kernels/optimisation/constrained.f90
)
target_link_libraries(optimisation core)

# PDE solvers
add_library(pde
    ../kernels/pde/finite_difference.f90
    ../kernels/pde/finite_volume.f90
    ../kernels/pde/spectral.f90
)
target_link_libraries(pde core)

# Physics models
add_library(models
    ../models/heat_transfer.f90
    ../models/fluid_dynamics.f90
    ../models/nuclear_fission.f90
    ../models/pressure_dynamics.f90
)
target_link_libraries(models core linalg pde)

# Test executable
add_executable(test_suite test_suite.f90)
target_link_options(test_suite PRIVATE "-Wl,-z,noexecstack")
target_link_libraries(test_suite 
    core 
    linalg 
    ode 
    optimisation 
    pde 
    models
    ${BLAS_LIBRARIES} 
    ${LAPACK_LIBRARIES}
)

# Enable testing
enable_testing()
add_test(NAME FullTestSuite COMMAND test_suite)

# Custom target to run tests
add_custom_target(run_tests
    COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure --verbose
    DEPENDS test_suite
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

# Print configuration
message(STATUS "Fortran Compiler: ${CMAKE_Fortran_COMPILER}")
message(STATUS "Build Type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Module Directory: ${CMAKE_Fortran_MODULE_DIRECTORY}")