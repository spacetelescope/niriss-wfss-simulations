# NIRISS/WFSS Simulations

[![astropy](http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat)](http://www.astropy.org/)

Code to simulate NIRISS/WFSS data. As of 2018, all scripts are in a state of active development and results may vary as a result. Maintained by the NIRISS branch. 

Usage
-----

Very specific versions of the following required packages are needed in order for the simulation code to run:
* [grizli](https://github.com/jotaylor/grizli)
* the [JWST calibration pipeline](https://github.com/jotaylor/jwst)
* [webbpsf](https://github.com/mperrin/webbpsf)

To ensure package compatibility, one solution is to create a separate environment for simulation work:

    conda env create -f jwstsim_env.yaml -n jwstsim
    source activate jwstsim

There is a [known issue](https://github.com/spacetelescope/pysiaf/issues/69) with `pysiaf` that can occur and gives the following error:

    ImportError [...] Library not loaded: libxml2.2.dylib Reason: Incompatible library version: etree.[...] requires version 12.0.0 or later, but libxml2.2.dylib provides version 10.0.0

A stop-gap solution is to uninstall `pysiaf` and `lxml` and reinstall `pysiaf`:

    pip uninstall pysiaf
    pip uninstall lxml
    pip install pysiaf

