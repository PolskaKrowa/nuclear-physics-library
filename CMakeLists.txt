cmake_minimum_required(VERSION 3.20)
project(NuclearPhysics VERSION 0.2.17 LANGUAGES Fortran C)

# Set Fortran standards and flags
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/mod)
include_directories(${CMAKE_Fortran_MODULE_DIRECTORY})

set(CMAKE_CONFIGURATION_TYPES "Debug;Release;Profiling" CACHE STRING "" FORCE)

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -std=f2008 -fimplicit-none -Wall -Wextra -ffree-line-length-none")
    set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -fcheck=all -fbacktrace -fimplicit-none -fno-omit-frame-pointer -pedantic ")
    set(CMAKE_Fortran_FLAGS_RELEASE "-Ofast -march=native -mtune=native -flto -ffpe-summary=none -fopenmp")
    set(CMAKE_Fortran_FLAGS_PROFILING "-O3 -g -pg -fopenmp")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -stand f08 -implicitnone -warn all -ffree-line-length-none")
    set(CMAKE_Fortran_FLAGS_DEBUG "-g -check all -traceback")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -xHost")
endif()

add_compile_definitions(_POSIX_C_SOURCE=200809L)
add_compile_definitions(M_PI=3.14159265358979)

# Option to manually specify OpenBLAS path
set(OPENBLAS_PATH "/usr/lib/libopenblas.so" CACHE FILEPATH "Path to OpenBLAS library")

# Try to find BLAS/LAPACK using CMake's built-in finders
# First, try to explicitly use OpenBLAS
set(BLA_VENDOR OpenBLAS)
find_package(BLAS QUIET)
find_package(LAPACK QUIET)

# If not found, try generic search
if(NOT BLAS_FOUND OR NOT LAPACK_FOUND)
    set(BLA_VENDOR All)
    find_package(BLAS QUIET)
    find_package(LAPACK QUIET)
endif()

# If still not found, use manual fallback to OpenBLAS
if(NOT BLAS_FOUND)
    if(EXISTS "${OPENBLAS_PATH}")
        set(BLAS_LIBRARIES "${OPENBLAS_PATH}")
        set(BLAS_FOUND TRUE)
        message(STATUS "Manually set BLAS to OpenBLAS: ${BLAS_LIBRARIES}")
    else()
        message(FATAL_ERROR "BLAS not found! Please install OpenBLAS or set OPENBLAS_PATH")
    endif()
endif()

if(NOT LAPACK_FOUND)
    if(EXISTS "${OPENBLAS_PATH}")
        set(LAPACK_LIBRARIES "${OPENBLAS_PATH}")
        set(LAPACK_FOUND TRUE)
        message(STATUS "Manually set LAPACK to OpenBLAS: ${LAPACK_LIBRARIES}")
    else()
        message(FATAL_ERROR "LAPACK not found! Please install OpenBLAS or set OPENBLAS_PATH")
    endif()
endif()

# Diagnostic output
message(STATUS "BLAS libraries: ${BLAS_LIBRARIES}")
message(STATUS "LAPACK libraries: ${LAPACK_LIBRARIES}")

# Find OpenGL and GLUT (optional for visualisation)
option(BUILD_WITH_OPENGL "Build OpenGL visualisation examples" ON)

if(BUILD_WITH_OPENGL)
    find_package(OpenGL)
    find_package(GLUT)
    
    if(OpenGL_FOUND AND GLUT_FOUND)
        message(STATUS "OpenGL found: ${OPENGL_LIBRARIES}")
        message(STATUS "GLUT found: ${GLUT_LIBRARIES}")
        set(OPENGL_AVAILABLE TRUE)
    else()
        message(STATUS "OpenGL or GLUT not found - OpenGL examples will not be built")
        message(STATUS "  To enable, install: sudo apt-get install freeglut3-dev mesa-common-dev")
        set(OPENGL_AVAILABLE FALSE)
    endif()
else()
    set(OPENGL_AVAILABLE FALSE)
    message(STATUS "OpenGL support disabled by user")
endif()

# Core library
add_subdirectory(core)

# Kernel library (depends on core)
add_subdirectory(kernels)
target_link_libraries(nuclear_physics_kernels nuclear_physics_core ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})

add_subdirectory(tests)

# Physics-models library (depends on kernels and core).
# Lives under physics/ since the directory was renamed from models/ as part of
# the phase-1 reorg (see docs/design/bwr4_physics_spec/REORG_PLAN.md step 1).
# Library target name kept as nuclear_physics_models for C-ABI continuity.
add_subdirectory(physics)
target_link_libraries(nuclear_physics_models nuclear_physics_kernels nuclear_physics_core)

# Subsystems library — one module per spec chapter (vessel, fuel, crd,
# recirculation, main_steam, feedwater) implementing the 5-call contract
# from REORG_PLAN §3. Step 2 wires only the skeleton; subsystem-level
# behaviour migrates in steps 3-7.
add_subdirectory(subsystems)
# Step 11 (reorg): fuel.f90 imports two_phase_flow for the heat↔two-phase
# convection-coefficient helper, so subsystems must compile after models.
target_link_libraries(nuclear_physics_subsystems nuclear_physics_models nuclear_physics_core)

# C Interface library
add_subdirectory(interface)

# BWR Core Simulator with OpenGL (if available)
if(OPENGL_AVAILABLE)
    add_executable(bwr_core_simulator_opengl examples/bwr_core_simulator_opengl.f90)
    
    target_link_libraries(bwr_core_simulator_opengl 
        nuclear_physics_models 
        nuclear_physics_kernels 
        nuclear_physics_core
        ${LAPACK_LIBRARIES} 
        ${BLAS_LIBRARIES}
        ${OPENGL_LIBRARIES}
        ${GLUT_LIBRARIES})
    
    target_include_directories(bwr_core_simulator_opengl PRIVATE 
        include
        ${OPENGL_INCLUDE_DIR}
        ${GLUT_INCLUDE_DIR})
    
    # Add definitions for OpenGL
    target_compile_definitions(bwr_core_simulator_opengl PRIVATE USE_OPENGL)
    
    install(TARGETS bwr_core_simulator_opengl RUNTIME DESTINATION bin)
    
    message(STATUS "OpenGL visualisation enabled - bwr_core_simulator_opengl will be built")
endif()

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
message(STATUS "")
message(STATUS "  Nuclear Physics Library Configuration")
message(STATUS "═════════════════════════════════════════════════════")
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
message(STATUS "    - nuclear_physics_subsystems")
message(STATUS "    - nuclear_physics_c (C interface)")
message(STATUS "")
message(STATUS "  Examples:")
if(OPENGL_AVAILABLE)
    message(STATUS "    - bwr_core_simulator_opengl (3D visualisation)")
else()
    message(STATUS "    - bwr_core_simulator_opengl (not available - missing OpenGL/GLUT)")
endif()
message(STATUS "═════════════════════════════════════════════════════")
message(STATUS "")