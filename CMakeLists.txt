cmake_minimum_required(VERSION 3.20)
project(NuclearPhysics VERSION 1.0.0 LANGUAGES Fortran C)

# Set Fortran standards and flags
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/mod)
include_directories(${CMAKE_Fortran_MODULE_DIRECTORY})

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -std=f2008 -fimplicit-none -Wall -Wextra")
    set(CMAKE_Fortran_FLAGS_DEBUG "-g -fcheck=all -fbacktrace")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -march=native")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -stand f08 -implicitnone -warn all")
    set(CMAKE_Fortran_FLAGS_DEBUG "-g -check all -traceback")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -xHost")
endif()

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -llapack -lopenblas")

add_compile_definitions(_POSIX_C_SOURCE=200809L)
add_compile_definitions(M_PI=3.14159265358979)

# Build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Find BLAS/LAPACK
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

# If CMake didn't find them, set defaults to OpenBLAS
if(NOT BLAS_FOUND)
    set(BLAS_LIBRARIES "/usr/lib/libopenblas.so")
    message(STATUS "Using OpenBLAS for BLAS: ${BLAS_LIBRARIES}")
endif()

if(NOT LAPACK_FOUND)
    set(LAPACK_LIBRARIES "/usr/lib/libopenblas.so")
    message(STATUS "Using OpenBLAS for LAPACK: ${LAPACK_LIBRARIES}")
endif()

# Core library
add_subdirectory(core)

# Kernel library (depends on core)
add_subdirectory(kernels)
target_link_libraries(nuclear_physics_kernels nuclear_physics_core ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})

add_subdirectory(tests)

# Models library (depends on kernels and core)
add_subdirectory(models)
target_link_libraries(nuclear_physics_models nuclear_physics_kernels nuclear_physics_core)

# C Interface library
add_library(nuclear_physics_c SHARED interface/c_interface.f90)
target_link_libraries(nuclear_physics_c 
    nuclear_physics_models 
    nuclear_physics_kernels 
    nuclear_physics_core
    ${LAPACK_LIBRARIES} 
    ${BLAS_LIBRARIES})

set_target_properties(nuclear_physics_c PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}
    VERSION ${PROJECT_VERSION}
    SOVERSION 1)

# Install C interface
install(TARGETS nuclear_physics_c
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib)

install(FILES include/nuclear_physics.h
    DESTINATION include)

# BWR Simulator example
add_executable(bwr_simulator examples/bwr_simulator.c)
target_link_libraries(bwr_simulator nuclear_physics_c m)
target_include_directories(bwr_simulator PRIVATE include)

# Install example
install(TARGETS bwr_simulator
    RUNTIME DESTINATION bin)

# Export configuration
install(EXPORT NuclearPhysicsTargets
    FILE NuclearPhysicsTargets.cmake
    NAMESPACE NuclearPhysics::
    DESTINATION lib/cmake/NuclearPhysics)

# Create package config
include(CMakePackageConfigHelpers)
configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/NuclearPhysicsConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/NuclearPhysicsConfig.cmake
    INSTALL_DESTINATION lib/cmake/NuclearPhysics)

write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/NuclearPhysicsConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion)

install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/NuclearPhysicsConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/NuclearPhysicsConfigVersion.cmake
    DESTINATION lib/cmake/NuclearPhysics)

# Installation summary
message(STATUS "")
message(STATUS "═══════════════════════════════════════════════════════")
message(STATUS "  Nuclear Physics Library Configuration")
message(STATUS "═══════════════════════════════════════════════════════")
message(STATUS "  Version: ${PROJECT_VERSION}")
message(STATUS "  Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "  Fortran compiler: ${CMAKE_Fortran_COMPILER_ID}")
message(STATUS "  C compiler: ${CMAKE_C_COMPILER_ID}")
message(STATUS "  Install prefix: ${CMAKE_INSTALL_PREFIX}")
message(STATUS "")
message(STATUS "  Libraries:")
message(STATUS "    - nuclear_physics_core")
message(STATUS "    - nuclear_physics_kernels")
message(STATUS "    - nuclear_physics_models")
message(STATUS "    - nuclear_physics_c (C interface)")
message(STATUS "")
message(STATUS "  Examples:")
message(STATUS "    - bwr_simulator (Boiling Water Reactor)")
message(STATUS "═══════════════════════════════════════════════════════")
message(STATUS "")