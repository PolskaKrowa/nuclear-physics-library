set(CORE_SOURCES
    kinds.f90
    constants.f90
    numeric_utils.f90
    rng.f90)

add_library(nuclear_physics_core ${CORE_SOURCES})

target_include_directories(nuclear_physics_core
    PUBLIC
    $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>
    $<INSTALL_INTERFACE:include>)

set_target_properties(nuclear_physics_core PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY})

# Export this target
install(TARGETS nuclear_physics_core
    EXPORT NuclearPhysicsTargets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib)
