<div align="center">

# ⚛️ Nuclear Physics Library

### High-Performance Fortran Toolkit for Reactor Simulations & Multiphysics Studies

[![GitHub Actions Workflow Status](https://img.shields.io/github/actions/workflow/status/PolskaKrowa/nuclear-physics-library/.github%2Fworkflows%2Fmain.yml)](https://github.com/PolskaKrowa/nuclear-physics-library/actions)
[![License](https://img.shields.io/github/license/PolskaKrowa/nuclear-physics-library)](LICENSE)
[![Release](https://img.shields.io/github/v/release/PolskaKrowa/nuclear-physics-library)](https://github.com/PolskaKrowa/nuclear-physics-library/releases)
[![Fortran](https://img.shields.io/badge/Fortran-734F96?logo=fortran&logoColor=white)](https://fortran-lang.org/)

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
| **Fluid Dynamics** | Incompressible Navier-Stokes, projection method, natural convection, multi-region support |
| **Heat Transfer** | Diffusion-convection with sources, multi-material regions, coupled flow |
| **Nuclear Fission** | Point kinetics, multi-group diffusion, decay heat (ANS-5.1), 6-group precursors |
| **Pressure Dynamics** | Multiple EOS, acoustic waves, phase transitions, compressibility |

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
│   └── rng.f90                # Random number generation
│
├── 🧮 kernels/                 # Mathematical methods
│   ├── linear_algebra/        # Matrix operations
│   ├── ode/                   # ODE solvers
│   ├── optimisation/          # Optimisation algorithms
│   └── pde/                   # PDE methods
│
└── ⚡ models/                  # Physics models
    ├── burnup_depletion.f90
    ├── cross_sections.f90
    ├── multigroup_diffusion.f90
    ├── two_phase_flow.f90
    ├── fluid_dynamics.f90
    ├── heat_transfer.f90
    ├── nuclear_fission.f90
    └── pressure_dynamics.f90
```

## 📖 Dependencies

### Required
- ✅ Fortran compiler (gfortran 9.0+ or Intel Fortran)
- ✅ CMake 3.10+
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
| Line Length | 100 characters maximum |
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

## 🗺️ Roadmap

### 🔮 Planned Features
- [ ] MPI parallelisation for distributed computing
- [ ] Python bindings via f2py
- [ ] I/O and Visualisation (HDF5, VTK, etc.)
- [ ] Reactor control modules
- [ ] Material Properties database
- [ ] Detailed fuel modelling
- [ ] Enhanced testing suite (unit tests, benchmarks, verification tests)

### 🚧 In Progress
- [ ] Additional reactor geometries (hexagonal, cylindrical)
- [ ] Advanced turbulence models
- [x] OpenMP threading for shared-memory systems
- [x] Improved performance for real-time/faster-than-real-time simulation

## 📋 Version History

<details>
<summary><b>v0.2.0 (current)</b></summary>

- Enhanced neutronics models
- Improved fluid dynamics simulation
- Consistent testing suite
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

Made with ⚛️ by [Stevenson Parker](https://github.com/PolskaKrowa)

</div>
