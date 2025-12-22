# nuclear physics library
A library for the assistance of the development of nuclear reactor simulations

---
A lot of the code for this library was ~~stolen~~ borrowed from my Distributed Computing Project

## installation

Installing the library can be done using the following commands:

```bash
git clone https://github.com/PolskaKrowa/nuclear-physics-library.git
cd nuclear-physics-library/
mkdir build && cd build
cmake ..
make -j$(nproc)
sudo make install # install built library files to the system's gcc/gfortran compiler
```