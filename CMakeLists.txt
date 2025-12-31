cmake_minimum_required(VERSION 3.20)
project(NuclearPhysics VERSION 1.0.0 LANGUAGES Fortran C)

# Set Fortran standards and flags
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/mod)
include_directories(${CMAKE_Fortran_MODULE_DIRECTORY})

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -std=f2008 -fimplicit-none -Wall -Wextra")
    set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -fcheck=all -fbacktrace -fimplicit-none -fno-omit-frame-pointer -pedantic")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -march=native")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -stand f08 -implicitnone -warn all")
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
add_subdirectory(interface)

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