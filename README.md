# Rapid Mineos

Uses [MINEOS](https://github.com/geodynamics/mineos) to calculate eigenfunctions but returns them in a buffer without writing files.

## Installation

Build system is configured by [cmake](https://cmake.org/).

```bash
mkdir build
cd build
cmake ..
make
```

This will produce the `rapid_mineos` executable in the `build` directory, which takes a path to a model file and output directory as arguments

## Usage

From the `build` directory, run

```bash
./rapid_mineos <path/to/earth/model> <output/directory>
```

The input earth model file should be of the same format as used by `minos_bran` in MINEOS.
