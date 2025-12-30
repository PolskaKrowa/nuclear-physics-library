# Nuclear Physics Library
[![Build Status](https://github.com/PolskaKrowa/nuclear-physics-library/actions/workflows/build.yml/badge.svg)](https://github.com/PolskaKrowa/nuclear-physics-library/actions)

A somewhat high-performance Fortran library for nuclear reactor simulation and computational physics. This library provides a comprehensive suite of numerical methods, physics models, and mathematical utilities for developing multiphysics reactor simulations.

## Features

### Core Utilities
- **Precision Control**: Standardised kind parameters (single, double, quad precision)
- **Physical Constants**: Mathematical and physical constants with high accuracy
- **Random Number Generation**: Deterministic RNG with multiple distributions
- **Numerical Utilities**: Array operations, interpolation, integration, and safety functions

### Mathematical Kernels

#### Linear Algebra
- Dense matrix operations using BLAS/LAPACK
- Eigenvalue solvers (symmetric, general, generalised)
- Linear system solvers (LU, Cholesky, iterative methods)
- Sparse matrix operations (CSR format)
- SVD, QR, and Schur decompositions

#### ODE Solvers
- **RK4**: Classic fourth-order Runge-Kutta
- **Dormand-Prince (DOPRI5)**: Adaptive step-size method
- **Backward Euler**: Implicit method for stiff equations

#### Optimisation
- Gradient descent (with momentum, Adam)
- Conjugate gradient methods (Fletcher-Reeves, Polak-Ribière)
- Quasi-Newton methods (BFGS, L-BFGS)
- Constrained optimisation (projected gradient, penalty methods)

#### PDE Methods
- **Finite Difference**: Various accuracy orders and boundary conditions
- **Finite Volume**: Conservation law solvers with multiple flux schemes
- **Spectral Methods**: FFT-based derivatives and Poisson solvers

### Physics Models

#### Fluid Dynamics
- Incompressible Navier-Stokes equations
- Projection method for pressure-velocity coupling
- Natural convection (Boussinesq approximation)
- Multiple fluid regions with different properties

#### Heat Transfer
- Diffusion-convection equation with source terms
- Multiple material regions
- Various boundary conditions (Dirichlet, Neumann, convective)
- Coupled with fluid flow

#### Nuclear Fission
- Point kinetics equations with delayed neutrons
- Multi-group neutron diffusion
- Fission power density calculation
- Decay heat generation (ANS-5.1 standard)
- 6-group delayed neutron precursor model

#### Pressure Dynamics
- Multiple equations of state (ideal gas, incompressible, stiffened gas)
- Acoustic wave propagation
- Liquid/vapour phase transitions
- Compressibility effects

## Installation

### Prerequisites

- **Fortran Compiler**: gfortran 9.0+ or ifort
- **CMake**: Version 3.10 or higher
- **BLAS/LAPACK**: For linear algebra operations
- **MSYS2**: If compiling on Windows

### Building from Source

```bash
# Clone the repository
git clone https://github.com/PolskaKrowa/nuclear-physics-library.git
cd nuclear-physics-library/

# Create build directory
mkdir build && cd build

# Configure with CMake
cmake ..

# Build the library (using all available cores)
make -j$(nproc)

# Install to system (requires sudo)
sudo make install
```

### Custom Installation Path

```bash
cmake -DCMAKE_INSTALL_PREFIX=/path/to/install ..
make -j$(nproc)
make install
```

### Build Options

```bash
# Release build with optimisations
cmake -DCMAKE_BUILD_TYPE=Release ..

# Debug build with symbols
cmake -DCMAKE_BUILD_TYPE=Debug ..

# Specify compiler
cmake -DCMAKE_Fortran_COMPILER=ifort ..
```

## Module Organisation

```
nuclear-physics-library/
├── core/                    # Foundational utilities
│   ├── kinds.f90           # Precision definitions
│   ├── constants.f90       # Physical/mathematical constants
│   ├── numeric_utils.f90   # Numerical utilities
│   └── rng.f90            # Random number generation
│
├── kernels/                # Mathematical methods
│   ├── linear_algebra/    # Matrix operations
│   ├── ode/              # ODE solvers
│   ├── optimisation/     # Optimisation algorithms
│   └── pde/              # PDE methods
│
└── models/               # Physics models
    ├── fluid_dynamics.f90
    ├── heat_transfer.f90
    ├── nuclear_fission.f90
    └── pressure_dynamics.f90
```

## Dependencies

### Required
- Fortran compiler (gfortran 9.0+ or Intel Fortran)
- CMake 3.10+
- BLAS/LAPACK libraries

### Optional
- FFTW3 (for optimised spectral methods)
- MPI (for parallel simulations)
- OpenMP (for shared-memory parallelisation)

### Installing Dependencies (Ubuntu/Debian)

```bash
sudo apt-get update
sudo apt-get install gfortran cmake libblas-dev liblapack-dev
```

### Installing Dependencies (macOS)

```bash
brew install gcc cmake openblas lapack
```

### Installing Dependencies (Red Hat/CentOS)

```bash
sudo yum install gcc-gfortran cmake blas-devel lapack-devel
```

## Performance Considerations

### Compiler Optimisation Flags

For production builds, use aggressive optimisation:

```bash
cmake -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_Fortran_FLAGS="-O3 -march=native -ffast-math" ..
```

### BLAS/LAPACK Performance

For best performance, use optimised BLAS implementations:
- **Intel MKL**: Best for Intel CPUs
- **OpenBLAS**: Good general-purpose choice
- **ATLAS**: Auto-tuned alternative

Example with MKL:
```bash
cmake -DBLA_VENDOR=Intel10_64lp ..
```

## Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

### Code Style

- Use 4-space indentation
- Maximum line length: 100 characters
- Use descriptive variable names
- Add comments for complex algorithms
- Follow existing module structure

## Licence

This project is licensed under the Apache V2.0 Licence - see the LICENCE file for details.

## Citation

If you use this library in your research, please cite:

```bibtex
@software{nuclear_physics_library,
  title = {Nuclear Physics Library: A Fortran Library for Reactor Simulation},
  author = {Stevenson Parker},
  year = {2025},
  url = {https://github.com/PolskaKrowa/nuclear-physics-library}
}
```

## Acknowledgements

- LAPACK and BLAS developers for linear algebra routines
- The Fortran community for continued language development
- Contributors and users of this library

## Contact

- **Issues**: [GitHub Issues](https://github.com/PolskaKrowa/nuclear-physics-library/issues)
- **Discussions**: [GitHub Discussions](https://github.com/PolskaKrowa/nuclear-physics-library/discussions)

## Roadmap

### Planned Features
- [ ] MPI parallelisation for distributed computing
- [ ] OpenMP threading for shared-memory systems
- [ ] Python bindings via f2py
- [ ] Additional reactor geometries (hexagonal, cylindrical)
- [ ] Multi-group neutron diffusion
- [ ] Advanced turbulence models
- [ ] Xenon/samarium poisoning models
- [ ] Fuel depletion and burnup tracking

### In Progress
- [x] Comprehensive test suite
- [ ] User documentation and tutorials
- [ ] Example reactor simulations

## Version History

### v0.1.0 (Current)
- Initial release
- Core utilities and mathematical kernels
- Basic physics models (fluid, heat, fission, pressure)
- Single-threaded implementation

---

**Note**: This library is intended for educational and research purposes.
