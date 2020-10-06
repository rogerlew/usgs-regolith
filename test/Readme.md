README - Tests
==============

Description
-----------
This directory contains a python notebook and input files for testing the accuracy of several intermediate and final outputs of the REGOLITH program.  The tests evaluate accuracy for synthetic terrain defined by intersecting sine waves.  Formulas to create the inputs and standards were derived symbolically using Mathematica 12.1.0.0 and coded in python 3.6 to save the outputs as ASCII grids.  The tests compare outputs from REGOLITH against the standards to verify that the mean of the residuals agree within a specified tolerance.  The purpose of the tests is to allow users to verify that the code is running correctly on their system.  

The synthetic terrain is defined by the formula

*z* = 50 + 25sin(&pi;*x*/50) + 25sin(&pi;*y*/50),

where, *z* is the ground-surface elevation in meters, *x* and *y* are horizontal coorinates in the east and north directions respectively. 

The python script tests the following quantities:
Aspect, Laplacian, Secant, Slope angle, LCSD, Gradient squared, NSD transport, NSD

Instructions for running the tests
----------------------------------
Compile REGOLITH either using the makefile provided in the 'src' directory, or the Fortran development environment of your choice.
Place either a copy of the REGOLITH executable file, or a symbolic link to it, in the test directory.
Launch Jupyter and Python on your system and open the notebook test_regolith.ipynb.
Run the notebook.
After a few seconds, test results will be printed below the last cell of the notebook indicating whether each of the eight functions passed or failed.
