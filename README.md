<div align="center">

# ⚛️ Nuclear Physics Library

### High-Performance Fortran Toolkit for Reactor Simulations & Multiphysics Studies

[![GitHub Actions Workflow Status](https://img.shields.io/github/actions/workflow/status/PolskaKrowa/nuclear-physics-library/.github%2Fworkflows%2Fmain.yml)](https://github.com/PolskaKrowa/nuclear-physics-library/actions)
[![License](https://img.shields.io/github/license/PolskaKrowa/nuclear-physics-library)](LICENSE)
[![Release](https://img.shields.io/github/v/release/PolskaKrowa/nuclear-physics-library)](https://github.com/PolskaKrowa/nuclear-physics-library/releases)
[![Fortran](https://img.shields.io/badge/Fortran-734F96?logo=fortran&lnvidiaogoColor=white)](https://fortran-lang.org/)
[![CUDA](https://img.shields.io/badge/NVIDIA_HPC_SDK-76B900?logo=nvidia&logoColor=white)](https://developer.nvidia.com/hpc-sdk)


[Features](#-features) • [Installation](#-installation) • [Documentation](#-module-organisation) • [Contributing](#-contributing) • [Roadmap](#-roadmap)

</div>

---

## 📚 Overview

Nuclear Physics Library is a comprehensive Fortran toolkit providing utilities, numerical solvers, and physics models to support the development of reactor simulations and multiphysics studies. Built for performance and accuracy, it offers a complete suite of mathematical kernels and physics models.

## ✨ Features

<details open>
<summary><b>🔧 Core Utilities</b></summary>

<br>

| Component | Description |
|-----------|-------------|
| **Precision Control** | Standardised kind parameters (single, double, quad precision) |
| **Physical Constants** | Mathematical and physical constants with high accuracy |
| **Random Number Generation** | Deterministic RNG with multiple distributions |
| **Numerical Utilities** | Array operations, interpolation, integration, and safety functions |

</details>

<details open>
<summary><b>🧮 Mathematical Kernels</b></summary>

<br>

### Linear Algebra
- ✅ Dense matrix operations using BLAS/LAPACK
- ✅ Eigenvalue solvers (symmetric, general, generalised)
- ✅ Linear system solvers (LU, Cholesky, iterative methods)
- ✅ Sparse matrix operations (CSR format)
- ✅ SVD, QR, and Schur decompositions

### ODE Solvers
| Method | Type | Use Case |
|--------|------|----------|
| **RK4** | Explicit | Classic fourth-order Runge-Kutta |
| **DOPRI5** | Adaptive | Dormand-Prince with adaptive step-size |
| **Backward Euler** | Implicit | Stiff equations |

### Optimisation Algorithms
- **Gradient Methods**: Gradient descent (with momentum, Adam)
- **Conjugate Gradient**: Fletcher-Reeves, Polak-Ribière
- **Quasi-Newton**: BFGS, L-BFGS
- **Constrained**: Projected gradient, penalty methods

### PDE Methods
- **Finite Difference**: Various accuracy orders and boundary conditions
- **Finite Volume**: Conservation law solvers with multiple flux schemes
- **Spectral Methods**: FFT-based derivatives and Poisson solvers

</details>

<details open>
<summary><b>⚡ Physics Models</b></summary>

<br>

| Model | Capabilities |
|-------|-------------|
| **Burnup Depletion** | Fuel burnup and isotopic depletion tracking |
| **Cross Sections** | Cross-section library with temperature and density feedback |
| **Heat Transfer** | Heat transfer between fuel and coolant |
| **Multigroup Diffsion** | Multi-group neutron diffusion solver |
| **Two-Phase Flow** | Two-Phase drift-flux model |

</details>

## 🚀 Installation

### Prerequisites

> [!IMPORTANT]
> Ensure you have the following installed:
> - **Fortran Compiler**: gfortran 9.0+ or Intel Fortran
> - **CMake**: Version 3.20 or higher
> - **BLAS/LAPACK**: For linear algebra operations
> - **MSYS2**: Required for Windows compilation

### Quick Start

```bash
# Clone the repository
git clone https://github.com/PolskaKrowa/nuclear-physics-library.git
cd nuclear-physics-library/

# Create and enter build directory
mkdir build && cd build

# Configure with CMake
cmake ..

# Build the library (using all available cores)
make -j$(nproc)

# Install to system (requires sudo)
sudo make install
```

### Custom Installation

```bash
cmake -DCMAKE_INSTALL_PREFIX=/path/to/install ..
make -j$(nproc)
make install
```

### Build Configurations

<table>
<tr>
<th>Configuration</th>
<th>Command</th>
<th>Use Case</th>
</tr>
<tr>
<td><b>Release</b></td>
<td><code>cmake -DCMAKE_BUILD_TYPE=Release ..</code></td>
<td>Production with optimisations</td>
</tr>
<tr>
<td><b>Debug</b></td>
<td><code>cmake -DCMAKE_BUILD_TYPE=Debug ..</code></td>
<td>Development with symbols</td>
</tr>
<tr>
<td><b>Intel Compiler</b></td>
<td><code>cmake -DCMAKE_Fortran_COMPILER=ifort ..</code></td>
<td>Intel-optimised builds</td>
</tr>
</table>

### Computation Modes (CPU / GPU / Hybrid)

> [!CAUTION]
> **nvfortran does not correctly generate intermediate .mod files to compile kinds.f90**
>
> This project's CMake based build system does not correctly generate the necessary module files from core/kinds.f90, Initially I thought it was a fault
> with nvfortran's handling of quad precision arithmetic, but that turned out to be false. nvfortran silently fails to compile the module and returns no error
> messages.
>
> ***Please do not create an issue if you fail to compile with GPU support. This is a known issue.***

This library supports three compilation modes that select where computation runs:

<table>
<tr>
<th>Mode</th>
<th>Command</th>
<th>Description</th>
</tr>
<tr>
<td><b>CPU</b> (default)</td>
<td><code>cmake -DCOMPUTATION_MODE=CPU ..</code></td>
<td>All computation on CPU (gfortran + OpenMP + OpenBLAS). No GPU required.</td>
</tr>
<tr>
<td><b>GPU</b></td>
<td><code>cmake -DCOMPUTATION_MODE=GPU ..</code></td>
<td>All GPU-eligible computation on NVIDIA GPU (nvfortran + CUDA Fortran). CPU handles only orchestration. Requires NVIDIA HPC SDK.</td>
</tr>
<tr>
<td><b>Hybrid</b></td>
<td><code>cmake -DCOMPUTATION_MODE=HYBRID ..</code></td>
<td>Work split between GPU and CPU based on problem size. Large parallel workloads go to GPU; small tasks stay on CPU. Requires NVIDIA HPC SDK.</td>
</tr>
</table>

> [!IMPORTANT]
> **GPU and Hybrid modes require:**
> - NVIDIA GPU (Compute Capability 6.0+)
> - NVIDIA HPC SDK (provides `nvfortran` with CUDA Fortran support)
> - Install from https://developer.nvidia.com/hpc-sdk

**What runs on the GPU:**

| Module | GPU Kernels |
|--------|-------------|
| Multigroup Diffusion | Jacobi source iteration sweep, fission source computation, flux reduction/normalization |
| Heat Transfer | Explicit heat equation step (7-point stencil + convection + cooling) |
| Two-Phase Flow | Per-cell void fraction, quality, slip ratio, CHF ratio |
| Burnup Depletion | Per-cell isotope depletion (Xe-135, I-135, Sm-149, Pm-149, U-235/238, Pu-239/241) |
| Linear Algebra | Matrix-vector, matrix-matrix multiply, AXPY, norm, dot product |
| Reductions | Sum, max, min, dot product (tree reduction with warp-level final step) |

**Hybrid mode tuning:**

The `compute_mode` module provides runtime tuning for hybrid mode:

```fortran
use compute_mode
call set_gpu_partition(0.7_wp)       ! 70% of work to GPU (default)
call set_gpu_min_workload(4096)      ! Skip GPU for problems < 4096 cells
call set_num_threads(4)              ! CPU thread count (0 = OMP default)
```

### Manual Library Specification

> [!NOTE]
> CMake may require manual specification of BLAS/LAPACK locations:

```bash
cmake -DBLAS_LIBRARIES=/path/to/libblas.so \
      -DLAPACK_LIBRARIES=/path/to/liblapack.so ..
```

> [!TIP]
> **OpenBLAS Support**: This programme supports OpenBLAS, which includes LAPACK in the same `.so` file. However, both libraries must still be passed to CMake for proper linking.

### Platform-Specific Dependencies

<details>
<summary><b>🐧 Ubuntu/Debian</b></summary>

```bash
sudo apt-get update
sudo apt-get install gfortran cmake libblas-dev liblapack-dev
```
</details>

<details>
<summary><b>🍎 macOS</b></summary>

```bash
brew install gcc cmake openblas lapack
```
</details>

<details>
<summary><b>🎩 Red Hat/CentOS</b></summary>

```bash
sudo yum install gcc-gfortran cmake blas-devel lapack-devel
```
</details>

## 📁 Module Organisation

```graphql
nuclear-physics-library/
├── 📦 core/                    # Foundational utilities
│   ├── kinds.f90              # Precision definitions
│   ├── constants.f90          # Physical/mathematical constants
│   ├── numeric_utils.f90      # Numerical utilities
│   ├── rng.f90                # Random number generation
│   └── compute_mode.f90       # CPU/GPU/Hybrid mode dispatch
│
├── 🧮 kernels/                 # Mathematical methods
│   ├── linear_algebra/        # Matrix operations
│   ├── ode/                   # ODE solvers
│   ├── optimisation/          # Optimisation algorithms
│   └── pde/                   # PDE methods
│
├── ⚡ models/                  # Physics models
│   ├── burnup_depletion.f90
│   ├── cross_sections.f90
│   ├── multigroup_diffusion.f90
│   ├── two_phase_flow.f90
│   └── heat_transfer.f90
│
└── 🎮 cuda/                    # CUDA Fortran GPU kernels
    ├── multigroup_diffusion_gpu.cuf
    ├── heat_transfer_gpu.cuf
    ├── two_phase_flow_gpu.cuf
    ├── burnup_depletion_gpu.cuf
    ├── reductions_gpu.cuf
    └── linear_algebra_gpu.cuf
```

## 📖 Dependencies

### Required
- ✅ Fortran compiler (gfortran 9.0+ or Intel Fortran)
- ✅ CMake 3.20+
- ✅ BLAS/LAPACK libraries

### Optional
- 🔄 FFTW3 (for optimised spectral methods)
- 🔄 MPI (for parallel simulations)
- 🔄 OpenMP (for shared-memory parallelisation)

## ⚡ Performance Considerations

### Compiler Optimisation Flags

The CMakeLists.txt file provides aggressive optimisation flags:

| Compiler | Flags |
|----------|-------|
| **GNU** | `-Ofast -march=native -mtune=native -flto` |
| **Intel** | `-O3 -xHost` |

### BLAS/LAPACK Performance

> [!TIP]
> For optimal performance, use these BLAS implementations:
> - **Intel MKL**: Best for Intel CPUs
> - **OpenBLAS**: Good general-purpose choice
> - **ATLAS**: Auto-tuned alternative

**Example with MKL:**
```bash
cmake -DBLA_VENDOR=Intel10_64lp ..
```

## 🤝 Contributing

Contributions are welcome! Please follow these steps:

1. 🍴 Fork the repository
2. 🌿 Create a feature branch (`git checkout -b feature/amazing-feature`)
3. 💾 Commit your changes (`git commit -m 'Add amazing feature'`)
4. 📤 Push to the branch (`git push origin feature/amazing-feature`)
5. 🔀 Open a Pull Request

### Code Style Guidelines

| Guideline | Requirement |
|-----------|-------------|
| Indentation | 4 spaces |
| Line Length | 100 characters maximum (excluding indentation) |
| Naming | Descriptive variable names |
| Documentation | Comments for complex algorithms |
| Structure | Follow existing module organisation |

## 📄 Licence

This project is licensed under the Apache V2.0 Licence - see the [LICENCE](LICENSE) file for details.

## 📝 Citation

If you use this library in your research, please cite:

```bibtex
@software{nuclear_physics_library,
  title = {Nuclear Physics Library: A Fortran Library for Reactor Simulation},
  author = {Stevenson Parker},
  year = {2025},
  url = {https://github.com/PolskaKrowa/nuclear-physics-library}
}
```

## 🙏 Acknowledgements

- LAPACK and BLAS developers for linear algebra routines
- The Fortran community for continued language development
- Contributors and users of this library

## 📬 Contact

- **Issues**: [GitHub Issues](https://github.com/PolskaKrowa/nuclear-physics-library/issues)
- **Discussions**: [GitHub Discussions](https://github.com/PolskaKrowa/nuclear-physics-library/discussions)

## 📋 Version History

<details open>
<summary><b>v0.4.0 (current)</b></summary>

- Added 3-mode compilation system: CPU, GPU, Hybrid GPU/CPU
- Added CUDA Fortran kernels for heat transfer, two-phase flow, burnup depletion, linear algebra, and reductions
- Enhanced multigroup diffusion GPU with fission source and flux reduction kernels
- Added `compute_mode` module for runtime GPU/CPU dispatch
- Performance: eliminated triple-copy BLAS pattern in dense_matrix.f90
- Performance: fixed sparse_matrix_mult dense-intermediate bug (Gustavson SparseGEMM)
- Performance: added OpenMP to heat_transfer, two_phase_flow, burnup_depletion, finite_difference
- Performance: refactored eigen inverse_iteration to use dgetrf+dgetrs (factor once)
- Performance: added LAPACK workspace queries for optimal lwork
- Performance: vectorised rng_normal_array (batched Box-Muller)
- Performance: removed dead identity matrix in BFGS, switched to BLAS dger
- Performance: optimised backward_euler FD Jacobian build (eliminated O(n²) copy)
</details>

<details>
<summary><b>v0.3.0</b></summary>

- Majorly reworked simulation (thanks to @fxllenfx)
- Added subsystems
- Fixed an issue with heat transfer using unstable solver
- Further improved performance
</details>

<details>
<summary><b>v0.2.0</b></summary>

- Enhanced neutronics models
- Improved fluid dynamics simulation
- Consistent testing suite
- Implement basic multithreading for neutron transport
</details>

<details>
<summary><b>v0.1.0</b></summary>

- Initial release
- Core utilities and mathematical kernels
- Basic physics models (fluid, heat, fission, pressure)
- Single-threaded implementation
</details>

---

<div align="center">

**Note**: This library is intended for educational and research purposes.

Founded by [Steve](https://github.com/PolskaKrowa)

</div>
