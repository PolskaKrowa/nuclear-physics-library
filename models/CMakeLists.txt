set(MODEL_SOURCES
    heat_transfer.f90
    two_phase_flow.f90
    multigroup_diffusion.f90
    cross_sections.f90
    burnup_depletion.f90)

add_library(nuclear_physics_models ${MODEL_SOURCES})

target_include_directories(nuclear_physics_models
    PUBLIC
    $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>
    $<INSTALL_INTERFACE:include>)

set_target_properties(nuclear_physics_models PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY})

# Export this target
install(TARGETS nuclear_physics_models
    EXPORT NuclearPhysicsTargets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib)