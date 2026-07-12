# cuda/CMakeLists.txt
#
# Builds the nuclear_physics_cuda library from .cuf files.
# Only included when COMPUTATION_MODE is GPU or HYBRID (see top-level CMakeLists.txt).

set(CUDA_SOURCES
    multigroup_diffusion_gpu.cuf
    reductions_gpu.cuf
    heat_transfer_gpu.cuf
    two_phase_flow_gpu.cuf
    burnup_depletion_gpu.cuf
    linear_algebra_gpu.cuf
)

# Tell CMake to treat .cuf as Fortran source so the nvfortran driver
# picks them up. This MUST be set in the same scope where the files
# are referenced (i.e. here, not in the parent CMakeLists.txt), because
# CMake tracks source file properties by the exact path string used.
# Setting it in the parent with absolute paths does not propagate to
# files referenced here by relative path.
set_source_files_properties(${CUDA_SOURCES} PROPERTIES LANGUAGE Fortran)

add_library(nuclear_physics_cuda ${CUDA_SOURCES})

# Explicitly set the linker language to Fortran. Without this, CMake
# cannot determine the linker language (because .cuf is not a recognised
# extension) and refuses to export the target with the error:
#   "Exporting the target ... is not allowed since its linker language
#    cannot be determined"
set_target_properties(nuclear_physics_cuda PROPERTIES
    LINKER_LANGUAGE Fortran
    Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}
    LINK_FLAGS "-cuda")

# Reuse the top-level module directory so .mod files coexist with host ones.
target_include_directories(nuclear_physics_cuda
    PUBLIC
    $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>
    $<INSTALL_INTERFACE:include>)

# Link against the host kinds module so the .cuf files see wp = real64.
# nuclear_physics_core provides the kinds module and compute_mode module.
target_link_libraries(nuclear_physics_cuda PUBLIC nuclear_physics_core)

# CUDA Fortran needs the CUDA runtime library. nvfortran links this
# automatically when -cuda is on the link line (set via LINK_FLAGS above),
# but we make it explicit so consumers of the library don't have to know.
if(NVFORTRAN_EXECUTABLE)
    target_link_options(nuclear_physics_cuda PRIVATE -cuda)
endif()

# Propagate the computation mode definition to the CUDA library.
if(COMPUTATION_MODE STREQUAL "GPU")
    target_compile_definitions(nuclear_physics_cuda PUBLIC COMPUTATION_MODE_GPU)
elseif(COMPUTATION_MODE STREQUAL "HYBRID")
    target_compile_definitions(nuclear_physics_cuda PUBLIC COMPUTATION_MODE_HYBRID)
endif()

install(TARGETS nuclear_physics_cuda
    EXPORT NuclearPhysicsTargets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib)