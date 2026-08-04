# TempLat example programs

This directory contains example programs demonstrating the use of the TempLat library for lattice field theory simulations. Each example is self-contained and can be run independently.

To build the examples, create a build directory and run CMake:

```bash
mkdir build
cd build
cmake ..
make -j$(nproc)
```

You can then run the examples individually. For example, to run the Langevin dynamics example:

```bash
./langevin 128 100
```

To run with MPI, configure CMake with MPI support and use `mpirun`:

```bash
cmake -DMPI=ON ..
make -j$(nproc)
mpirun -np 4 ./langevin 128 100 
```

For more information on each example, please refer to the comments at the top of each source file.
The following examples are included:
- `average.cpp`: Computes the average of a field configuration with gaussian fluctuations over the lattice.
- `langevin.cpp`: Simulates Langevin dynamics for a scalar field with a quartic interaction.
- `classical_phi4.cpp`: Simulates classical dynamics of a scalar field with a quartic interaction.
- `classical_su2.cpp`: Simulates classical dynamics of an SU(2) gauge field on a lattice.
