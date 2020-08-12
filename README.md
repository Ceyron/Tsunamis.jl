# Tsunamis.jl - A Finite Volume Shallow Water Equation Solver for Julia

Future Trends in HPC Seminar - Shallow Water Equation Solvers with Julia

This is a translation of elements of the Shallow Water Equation code from
https://github.com/TUM-I5/SWE.

## Installing Dependencies

In a Julia REPL environment run

    ] add NetCDF NCDatasets ArgParse

## Running the code

### Sequential Run

Make sure all dependencies are installed, then run

    julia single_node.jl

for the single-node sequential implementation.

For help add the "-h" flag.

### Parallel Run

Choose the number of julia workers or provide the julia executable with a machine file.

    julia -p <NUM_CORES> multi_node.jl

For help add the "-h" flag.

## Postprocessing

If the output flag is set (default) either one NetCDF (sequential) or multiple (parallel run) NetCDF files are created. They can be postprocessed e.g. with ParaView.