# NIRISS/WFSS Simulations

[![astropy](http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat)](http://www.astropy.org/)

Code to simulate NIRISS/WFSS data. As of 2018, all scripts are in a state of active development and results may vary as a result. Maintained by the NIRISS branch. 

Note
----

Very specific versions of the following required packages are needed in order for the simulation code to run:
* [grizli](https://github.com/jotaylor/grizli)
* the [JWST calibration pipeline](https://github.com/jotaylor/jwst)
* [webbpsf](https://github.com/mperrin/webbpsf)

To ensure package compatibility, one solution is to create a separate environment for simulation work:

    conda env create -f jwstsim_env.yaml -n jwstsim
    source activate jwstsim

