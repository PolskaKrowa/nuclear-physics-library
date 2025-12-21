set(KERNELS_SOURCES
    linear_algebra/dense_matrix.f90
    linear_algebra/eigen.f90
    linear_algebra/solve_linear.f90
    linear_algebra/sparse_matrix.f90
    ode/rk4.f90
    ode/dormand_prince.f90
    ode/backward_euler.f90
    optimisation/gradient_descent.f90
    optimisation/conjugate_newton.f90
    optimisation/quasi_newton.f90
    optimisation/constrained.f90
    pde/finite_difference.f90
    pde/finite_volume.f90
    pde/spectral.f90)

add_library(nuclear_physics_kernels ${KERNELS_SOURCES})

# Link against core library
target_link_libraries(nuclear_physics_kernels 
    PUBLIC nuclear_physics_core)

target_include_directories(nuclear_physics_kernels
    PUBLIC
    $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>
    $<INSTALL_INTERFACE:include>)

set_target_properties(nuclear_physics_kernels PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY})

# Export this target
install(TARGETS nuclear_physics_kernels
    EXPORT NuclearPhysicsTargets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib)
