# INSE FIEM 2D #
Fast integral equation solver for the incompressible Navier-Stokes equations in 2D.

This code accompanies the paper "A Fast Integral Equation Method for the Two-Dimensional
Navier-Stokes Equations" by L. af Klinteberg, T. Askham and M. C. Kropinski.

Note that this software is not well-documented, and only partially user-friendly. It is mainly supplied as supplementary material for the abovementioned paper.

Also note that the fast direct solver for the modified Stokes equation is not included, due to licensing issues.

# Setup & Build

* Make sure you have the following installed:
    - Julia 0.6 (0.7 or higher will not work)
    - python-matplotlib
    - python-scipy
    - python-sympy
    - python-pip
    - Matlab
    - csh
    - hdf5-tools
* Run `julia julia/setup.jl` to install required Julia packages.
* To init submodules, run `git submodule update --init --recursive`
* Run `make` to build binaries.

# Running code

From the `julia` directory, run `julia <filename>`, e.g.

* `julia run_unit_tests.jl` to run unit test suite.
* `julia tests/test_XXXX.jl` to run a specific test.
* `julia examples/XXXXX.jl` to run one of the examples.
* `julia output_scripts/XXXXX.jl` to generate the plots shown in the paper.
