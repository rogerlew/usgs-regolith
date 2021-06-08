README - Tests
==============

Description
-----------
This directory contains a python notebook *test_regolith.ipynb*, and input files for testing the accuracy of several intermediate and final outputs of the REGOLITH program.  The tests evaluate accuracy for synthetic terrain defined by intersecting sine waves.  Formulas to create the inputs and standards were derived symbolically using Mathematica 12.1.0.0 (Wolfram Research, Inc., 2020) and coded in python 3.6 to save the outputs as ASCII grids.  The tests compare outputs from REGOLITH against the standards to verify that the mean of the residuals agree within a specified tolerance.  The purpose of the tests is to allow users to verify that the code is running correctly on their system.   A second notebook, *SyntheticTerrain.ipynb*, helps the user visualize many of the outputs computed during the tests.

The synthetic terrain is defined by the formula

*z* = 50 + 25sin(&pi;*x*/50) + 25sin(&pi;*y*/50),

where, *z* is the ground-surface elevation in meters, *x* and *y* are horizontal coorinates in the east and north directions respectively. 

The python script tests the following quantities:
Aspect, Laplacian, Secant, Slope angle, LCSD, Gradient squared, NSD transport, NSD

Instructions for running the tests
----------------------------------
Compile REGOLITH either using the makefile provided in the 'src' directory, or the Fortran development environment of your choice.
Place either a copy of the REGOLITH executable file, or a symbolic link to it, in the test directory.
Launch Jupyter and Python on your system and open the notebook *test_regolith.ipynb*.
Run the notebook.
After a few seconds, test results will be printed below the last cell of the notebook indicating whether each of the eight functions passed or failed.

Notebook *SyntheticTerrain.ipynb*
---------------------------------------
After the tests have completed successfully, the notebook *SyntheticTerrain.ipynb* can be opened and run to generate contour plots of (1) the synthetic terrain, (2) various derivatives and (3) output from analytical versions of several soil depth models implemented in REGOLITH.  The notebook also computes histograms of differences beween outputs from the notebook's  analytical formulas and numerical values that were output previously by REGOLITH when you ran *test_regolith.ipynb*. 

References
-------------
Wolfram Research Inc., 2020, Mathematica, version 12.1: Champaign, IL.  
