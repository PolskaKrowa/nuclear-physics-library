set(MODEL_SOURCES
    heat_transfer.F90
    two_phase_flow.F90
    multigroup_diffusion.F90
    cross_sections.F90
    burnup_depletion.F90)

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