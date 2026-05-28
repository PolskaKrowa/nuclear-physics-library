set(SUBSYSTEM_SOURCES
    vessel.f90
    fuel.f90
    crd.f90
    recirculation.f90
    main_steam.f90
    feedwater.f90
    rhr.f90
    condensate_loop.f90)

add_library(nuclear_physics_subsystems ${SUBSYSTEM_SOURCES})

target_include_directories(nuclear_physics_subsystems
    PUBLIC
    $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>
    $<INSTALL_INTERFACE:include>)

set_target_properties(nuclear_physics_subsystems PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY})

# Export this target
install(TARGETS nuclear_physics_subsystems
    EXPORT NuclearPhysicsTargets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib)