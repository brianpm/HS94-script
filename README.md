# HS94-script
A script to analyze Held-Suarez 1994 experimental data produced by the CAM-MPAS model

See https://www.cesm.ucar.edu/models/simpler-models/held-suarez.html for background on this test.
This script is an attempt to replace HS94 NCL script provided at that site with one based on Python and Xarrays.

The script is set up to run on:
a local system with a local history file path (hpath): here set to where I've stashed 5 history files on my laptop.

-OR-

on the casper server at NCAR, with hpath set to a specific directory on GLADE where a full 1200-day, HS94 dataset for CAM-MPAS running at 120 km (40962 columns and 32 levels resolution has been interpolated to a lat-lon grid.

The script presupposes that the geocat environment is available with Conda.

To install the geocat environment, follow these directions:
 >  conda create -n geocat -c conda-forge -c ncar geocat-comp geocat-datafiles
 >  conda activate geocat
 >  conda install -c conda-forge matplotlib cartopy jupyter
 >  conda install netcdf4

Running on your laptop would follow the pattern:

> conda activate geocat
> python interp.py
> conda deactivate

The basic run pattern in an HPC environment (for example on casper) is:

module load conda
conda activate geocat
python interp.py
conda deactivate


