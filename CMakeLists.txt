cmake_minimum_required(VERSION 3.15)
project(NuclearPhysicsLibrary 
    VERSION 1.0.0
    DESCRIPTION "A library for nuclear reactor simulations"
    LANGUAGES Fortran C)

set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules")
set(CMAKE_FORTRAN_FLAGS "${CMAKE_FORTRAN_FLAGS} -cpp")

# Set default build type to Release
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Enable testing
enable_testing()

# Organize output directories
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

# Add subdirectories
add_subdirectory(core)
add_subdirectory(kernels)
add_subdirectory(c_interface)

# Installation configuration
install(DIRECTORY "${CMAKE_BINARY_DIR}/modules/"
    DESTINATION include
    FILES_MATCHING PATTERN "*.mod")

install(TARGETS nuclear_physics_core
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib)

install(TARGETS nuclear_physics_kernels
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib)

# Create pkg-config file
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/nuclear-physics-library.pc.in"
    "${CMAKE_CURRENT_BINARY_DIR}/nuclear-physics-library.pc"
    @ONLY)

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/nuclear-physics-library.pc"
    DESTINATION lib/pkgconfig)

# Export targets for use by other projects
install(EXPORT NuclearPhysicsTargets
    FILE NuclearPhysicsTargets.cmake
    NAMESPACE NuclearPhysics::
    DESTINATION lib/cmake/NuclearPhysics)

configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/NuclearPhysicsConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/NuclearPhysicsConfig.cmake"
    @ONLY)

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/NuclearPhysicsConfig.cmake"
    DESTINATION lib/cmake/NuclearPhysics)
