# cuda/CMakeLists.txt
#
# Builds the nuclear_physics_cuda library from .cuf files.
# Only included when BUILD_WITH_CUDA=ON (see top-level CMakeLists.txt).

set(CUDA_SOURCES
    multigroup_diffusion_gpu.cuf
    # Future phases:
    # heat_transfer_gpu.cuf
    # cross_sections_gpu.cuf
    # two_phase_flow_gpu.cuf
)

add_library(nuclear_physics_cuda ${CUDA_SOURCES})

# Reuse the top-level module directory so .mod files coexist with host ones.
target_include_directories(nuclear_physics_cuda
    PUBLIC
    $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>
    $<INSTALL_INTERFACE:include>)

set_target_properties(nuclear_physics_cuda PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}
    LINK_FLAGS "-cuda")

# Link against the host kinds module so the .cuf files see wp = real64.
# nuclear_physics_core provides the kinds module.
target_link_libraries(nuclear_physics_cuda PUBLIC nuclear_physics_core)

# CUDA Fortran needs the CUDA runtime library. nvfortran links this
# automatically when -cuda is on the link line (set via LINK_FLAGS above),
# but we make it explicit so consumers of the library don't have to know.
if(NVFORTRAN_EXECUTABLE)
    target_link_options(nuclear_physics_cuda PRIVATE -cuda)
endif()

install(TARGETS nuclear_physics_cuda
    EXPORT NuclearPhysicsTargets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib)
