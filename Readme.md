README
======

Description
-----------
The computer program REGOLITH is a Fortran program designed for estimating soil thickness over a digital landscape.  The program implements various published empirical and computational models for soil thickness.  This command-line program is used in conjunction with Python scripts and shell scripts to prepare input files and visualize results.

Purpose and limitations
-----------------------
The program is designed to make rapid estimates of soil depth using several published models to to create input for spatially distributed slope stability models. The formulas are applied throughout the model domain to both upland hillslopes and upland valleys.  The program does not contain any code for statistically testing soil depth distributions against observed or measured soil depths.  

Compiling the program from source code
--------------------------------------
A simple makefile in the src directory can be used to compile the fortran source code.  The makefile has been tested on a linux system which had the GNU compiler collection (gcc) and gfortran installed.  The makefile should be usable after minor modifications with gfortran under Cygwin or OSX.  After downloading the repository, change directories to "src" and type `make` to compile regolith to create an executable file.

Using the program
-----------------
The program runs from the command line and can be compiled for any operating system that supports Fortran 95.  The user provides a small input file, *rg_in.txt*, that contains a list of model parameters and a list of path names for other input files.  The digital elevation model and related input files , slope, upslope contributing area, elevation index, plan-view curvature, and property zone are raster grids in ASCII grid format.  All input grids must have the same spatial reference, resolution, and footprint.  The program outputs a log file, an ASCII grid of computed soil depth, and optional grids of various intermediate values.  It will also output a grid of slope angle values if one does not exist or if the user specifies topographic smoothing.

Input parameter definitions
---------------------------
The following list explains the variable names, types, and permissible values of the input paramters in order of their appearance in the input file, *rg_in.txt*. 

Name,         Type,      Description
----------------------------------
- `title`        char,      Brief description of project, up to 224 characters
- `num_zones`    integer    Number of parameter zones to allow for variable soils across a landscape
- `max_zones`    integer    Highest zone number, to allow for using a subsample of a larger grid, where some zones of the larger grid are absent in the subsample.
- `h0`           float,     Characteristic soil depth, typically 0.3 - 0.5 m
- `sc`           float,     Angle of stability (in degrees)
- `dif_ratio`    float,     Diffusivity ratio, >0., (&rho;<sub>b</sub> x P<sub>0</sub>)/(&rho;<sub>s</sub> x D), where D is soil diffusivity, &rho;<sub>b</sub> and &rho;<sub>s</sub> are density of bedrock and soil, respectively, and P<sub>0</sub> is maximum bedrock lowering rate on a flat surface.
- `C1`           float,     Empirical constant, >0.
- `depth_max`    float,     Maximum soil depth, >0.
- `depth_min`    float,     Minimum soil depth, >=0.
- `C0`           float,     Calibration factor, >0.
- `C2`           float,     Empirical Constant, >0.
- `chan_thresh`  float,     Threshold (minumum) upslope contributing area at channel heads, >0.
- `chan_depth`   float,     Average steady-state depth of alluvium in steep channels (>=0.2*sc).
- `num_steps`    integer,   Number of (~1-cm) depth increments per unit depth (~100 - 500 for depth in meters)
- `trans_model`  char,      Alpha-numeric code to designate soil depth or transport model 
- Legal values of `trans_model`: DRS1, DRS2, DRS3, WNDX, LCSD, NSD, NSDA, NASD, or NDSD
- `hump_prod`    logical,   Use humped soil production function? (enabled if `.true.` or disabled if `.false.`)
- `power`        integer,   Exponent of DRS2 ploynomial or of uplsope area in NASD, NSDA, and WNDX models
- `elevfil`      char,      File path name of digital elevation grid (<=255 characters)
- `slopefil`     char,      File name of slope angle grid (in degrees) (<=255 characters)
- `flo_accfil`   char,      File name of flow-accumulation grid (<=255 characters)
- `ndxfil`       char,      File name of elevation index grid (<=255 characters)
- `pv_curvfil`   char,      File name of plan-view curvature grid (<=255 characters) 

- `zonfil`       char,      File name of parameter-zone grid (<=255 characters)
- `folder`       char,      Folder where output grid files will be stored (<=255 characters)
- `suffix`       char,      Identification code to be added to names of output files (<=16 char.)
- `lasc`         logical,   Specify grid file extension? (*.asc* if `.true.`, *.txt* if `.false`)
- `l_deriv`      logical,   Specify to output derivatives and related quantites
- `l_test`       logical,   Switch between test mode and production mode for the LCSD, NASD, NSD-A, and NSD models.  Test mode allows comparison with analytical solutions and production mode is for mapping soil depths using a modified form of the theoretical model. 
- `topoSmooth`   logical, Specify topographic smoothing (enabled if `.true.`,  diabled if `.false.`)
- `soilSmooth`   logical,  Specify smoothing of computed soil depth (enabled if `.true.`,  diabled if `.false.`)
- `n_points`     integer,  Window size, in number of grid cells, of a running average to approximate gaussian smoothing.

Suggested Input Parameter Values
---------------------------
Suggested Input Parameter Values for Empirical Models
-----------------------------------------------

| `trans_model` | `sc` (degrees)*| `C1` | `depth_max` (m)* | `depth_min` (m) | `C0` | `C2` | `power` | Reference |
| ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ |
| DRS1 | 30 - 60 | 0.124 | 2 - 5 | 0 - 0.1 | 5.0 | — | — | DeRose et al. (1991) |
| DRS2 | 30 - 60 | 0.019 - 0.022 | 2 - 5 | 0 - 0.1 | 1.57 - 1.68 | — | 3 | DeRose (1996) |
| DRS3 | 30 - 60 | 0.04 | 2 - 5 | 0 - 0.1 | 5.0 | 1.5 | — | Baum et al. (2011) |
| WNDX | 30 - 60 | — | 2 - 5 | 0 - 0.1 | 0.1- 0.3 | — | 1 | Ho et al. (2012) |

*Varies with terrain and/or climate

Suggested Input Parameter Values for Process-Based Models
-----------------------------------------------
`l_test` is available for verifying code output against analytical solutions. The production mode ranges for `dif_ratio` are based on tests on soils developed in volcanic and plutonic rocks in Puerto Rico while the test mode ranges were obtained from Pelletier & Rasmussen (2009) from tests on soils in Pima County, Arizona.

| `trans_model` | `h0` (m) | `sc` (degrees)* | `dif_ratio` `l_test` | `dif_ratio` | `depth_max` (m)* | `depth_min` (m) | `power` | Reference |
| ----- | ------ | ------ | ------ | ------ | ------ | ------- | ------ | ------ |
| LCSD | 0.3 - 0.5 | 30 - 60 | 1 - 5 | 0.005 - 0.1 | 2 - 5 | 0 - 0.1 | — | Pelletier & Rasmussen (2009) |
| NSD | 0.3 - 0.5 | 30 - 60 | 0.2 - 1 | 0.005 - 0.1 | 2 - 5 | 0 - 0.1 | — | Pelletier & Rasmussen (2009) |
| NSDA | 0.3 - 0.5 | 30 - 60 | 1 - 5 | 0.02 - 0.9 | 2 - 5 | 0 - 0.1 | 1 - 2 | Pelletier & Rasmussen (2009) |
| NASD | 0.3 - 0.5 | 30 - 60 | 1 - 5 | 0.15 - 1 | 2 - 5 | 0 - 0.1 | 1 - 2 | Pelletier & Rasmussen (2009) |
| NDSD | 0.3 - 0.5 | 30 - 60 | — | 0.1 - 3 | 2 - 5 | 0 - 0.1 | — | Pelletier & Rasmussen (2009) |

*Varies with terrain and/or climate


Suggested Input Parameter Values for all Models
-----------------------------------------------
| Parameter | Suggested Range |
| ------ | ------ |
| `chan_thresh` (m<sup>2</sup>) | 1500 - 2000 |
| `chan_depth` (m) | < 0.5 |


Inputs used by each soil depth model
------------------------------------

Model code,     Description and Reference(s)
-----------------------------------------------
- DRS1:         Empirical exponential slope dependence (DeRose et al. 1991)
- DRS2:         Empirical 3rd-degree polynomial slope dependence (DeRose et al. 1991)
- DRS3:         Empirical exponential slope and sign of plan-view curvature dependence (Baum et al. 2011)
- WNDX:         Modified soil wetness index, linear area dependent transport (Ho et al. 2012)
- LCSD:         Linear curvature and slope dependent transport (Pelletier & Rasmussen, 2009; Peletier et al. 2016)
- NSD:          Nonlinear slope dependent transport (Pelletier & Rasmussen, 2009)
- NSDA:         Nonlinear slope dependent transport (Pelletier & Rasmussen, 2009)
- NASD:         Nonlinear slope and area dependent transport (Pelletier & Rasmussen, 2009)
- NDSD:         Nonlinear slope and depth dependent transport (Pelletier & Rasmussen, 2009)

Soil Depth Model Formulas
------------------------------------
Empirical Model Formula Parameters and Regolith Input Parameter Names (if applicable)

- *d*<sub>*r*</sub>,     regolith depth (m)
- *C*<sub>*0*</sub>,     empirical constant (C0)
- *C*<sub>*1*</sub>,     calibration factor (C1)
- *C*<sub>*2*</sub>,     empirical constant (C2)
- &delta;,               mean slope angle (degrees)
- *p*,                   exponent of DRS2 polynomial or of upsloe area (power)
- &kappa;,               plan-view curvature of ground surface
- *A*,                   upslope contributing area normalized (m<sup>2</sup>)

Process Based Model Formula Parameters and Regolith Input Parameter Names (if applicable)
- *d*<sub>*r*</sub>,    regolith depth (m)
- *z*,                  ground surface (m)
- *D*,                  diffusivity ratio (dif_ratio)
- *h*<sub>*0*</sub>,    characteristic soil depth (h0, m)
- &theta;,              slope angle
- *S*<sub>*c*</sub>,    tangent of the angle of stability (sc, degrees)
- *A*,                  upslope contributing area (m<sup>2</sup>)
- *m*,                  exponent of upslope area (power)
- *h*<sub>*n*</sub>,      soil thickness normal to surface (m)


| `trans_model` | Regolith Depth Formula |
| ------ | ------ |
| DRS1 | ![equation](https://latex.codecogs.com/gif.latex?d_%7Br%7D%20%3D%20C_%7B0%7De%5E%7B-C_%7B1%7D%5Cdelta%7D) |
| DRS2 | ![equation](https://latex.codecogs.com/gif.latex?d_%7Br%7D%20%3D%20%28C_%7B0%7D-C_%7B1%7D%5Cdelta%29%5E%7Bp%7D) |
| DRS3 | ![equation](https://latex.codecogs.com/gif.latex?d_%7Br%7D%20%3D%20%28C_%7B0%7D-C_%7B2%7Dsgn%28%5Ckappa%20%29%29e%5E%7B-C_%7B1%7D%5Cdelta%20%7D) |
| WNDX | ![equation](https://latex.codecogs.com/gif.latex?d_%7Br%7D%20%3D%20C_%7B0%7Dln%28%5Cfrac%7BA%5E%7Bp%7D%7D%7Btan%28%5Cdelta%29%7D%29) |
| NSD | ![equation](https://latex.codecogs.com/gif.latex?d_%7Br%7D%20%3D%20%5Cfrac%7Bh_%7B0%7D%7D%7Bcos%28%5Ctheta%29%7Dln%28-%5Cfrac%7BD%7D%7Bcos%28%5Ctheta%29%5Ctriangledown%5Ccdot%5Cfrac%7B%5Ctriangledown%20z%7D%7B1-%28%5Cfrac%7B%7C%5Ctriangledown%20z%7C%7D%7BS_%7Bc%7D%7D%29%5E%7B2%7D%7D%7D%29) |
| NASD | ![equation](https://latex.codecogs.com/gif.latex?d_%7Br%7D%20%3D%20%5Cfrac%7Bh_%7B0%7D%7D%7Bcos%28%5Ctheta%29%7Dln%28-%5Cfrac%7BD%7D%7Bcos%28%5Ctheta%29%5Ctriangledown%5Ccdot%5Cfrac%7BA%5E%7Bm%7D%5Ctriangledown%20z%7D%7B1-%28%5Cfrac%7B%7C%5Ctriangledown%20z%7C%7D%7BS_%7Bc%7D%7D%29%5E%7B2%7D%7D%7D%29) |
| NSDA | ![equation](https://latex.codecogs.com/gif.latex?d_%7Br%7D%20%3D%20%5Cfrac%7Bh_%7B0%7D%7D%7Bcos%28%5Ctheta%29%7Dln%28-%5Cfrac%7BDA%5E%7B-m%7D%7D%7Bcos%28%5Ctheta%29%5Ctriangledown%5Ccdot%5Cfrac%7B%5Ctriangledown%20z%7D%7B1-%28%5Cfrac%7B%7C%5Ctriangledown%20z%7C%7D%7BS_%7Bc%7D%7D%29%5E%7B2%7D%7D%7D%29) |
| LCSD | ![equation](https://latex.codecogs.com/gif.latex?d_%7Br%7D%20%3D%20%5Cfrac%7Bh_%7B0%7D%7D%7Bcos%28%5Ctheta%29%7Dln%28-%5Cfrac%7BD%7D%7Bcos%28%5Ctheta%29%5Ctriangledown%5E%7B2%7Dz%7D%29) |
| NDSD | ![equation](https://latex.codecogs.com/gif.latex?d_%7Br%7D%20%3D%20%5Cfrac%7Bh_%7B0%7D%7D%7Bcos%28%5Ctheta%29%7Dln%28-%5Cfrac%7BD%7D%7Bcos%28%5Ctheta%29%5Ctriangledown%5Ccdot%5Cfrac%7Bh_%7Bn%7D%5Ctriangledown%20z%7D%7B1-%28%5Cfrac%7B%7C%5Ctriangledown%20z%7C%7D%7BS_%7Bc%7D%7D%29%5E%7B2%7D%7D%7D%29) |


Model code,    Required input parameters
-----------------------------------------------
- DRS1:         `sc, dif_ratio, depth_max, depth_min, chan_thresh, chan_depth`
- DRS2:         `sc, dif_ratio, depth_min, cal_fac, chan_thresh, chan_depth`
- DRS3:         `sc, dif_ratio, depth_max, depth_min, cal_fac, chan_thresh, chan_depth`
- WNDX:         `sc, depth_min, cal_fac, chan_thresh, chan_depth`
- LCSD:          `h0, sc, dif_ratio, depth_min, chan_thresh, chan_depth, hump_prod, l_test`
- NSD:          `h0, sc, dif_ratio, depth_min, chan_thresh, chan_depth, hump_prod, l_test`
- NSDA:         `h0, sc, dif_ratio, depth_min, chan_thresh, chan_depth, hump_prod, l_test`
- NASD:         `h0, sc, dif_ratio, depth_min, chan_thresh, chan_depth, hump_prod, l_test`
- NDSD:         `h0, sc, dif_ratio, depth_max, depth_min, chan_thresh, chan_depth, hump_prod, num_steps`

Model code,    Required input files
------------------------------------------
- DRS1:          `elevfil, slopefil, flo_accfil`
- DRS2:          `elevfil, slopefil, flo_accfil`
- DRS3:          `elevfil, slopefil, flo_accfil, pv_curvfil`
- WNDX:          `elevfil, slopefil, flo_accfil`
- LCSD:          `elevfil, slopefil, flo_accfil`
- NSD:           `elevfil, slopefil, flo_accfil`
- NSDA:          `elevfil, slopefil, flo_accfil`
- NASD:          `elevfil, slopefil, flo_accfil`
- NDSD:          `elevfil, slopefil, flo_accfil, ndxfil`

Input file variable name,   Contents
-------------------------------------
- `elevfil`        A gridded digital elevation model of type float 
- `slopefil`       Slope angle grid of type float
- `flo_accfil `    A grid of upslope contributing area of type float
- `ndxfil`         A grid of integer values ranking grid cells in the digital elevation model from highest (1) to lowest (maximum number of data cells in the grid).  This grid is used only by the NDSD model to ensure that soil depth is computed in order from peaks and ridge crests down the slope to successively lower points as explained by Pelletier and Rasmussen (2009).  
- ` pv_curvfil`    Plan-curvature grid of type float used only by the DRS3 model
- `zonfil`         A grid of integer values identifying zones within the model where similar values of input parameters are to be applied.

*Notes:* 
1. If `slopefil` does not exist or is not available, REGOLITH will generate a slope grid.
2. The word "none" can be used as a placeholder in *rg_in.txt* for names of files not needed for a particular soil model.

Optional property zones
-----------------------
REGOLITH has an option to read in a property zone grid for the area covered by the DEM.  This option will allow for application of a soil model over an area where bedrock type, microclimate, or other factors result in significant differences in soil development from one zone to another.  The option requires import of a zone grid, with integer values starting at 1 and ranging up to the number of zones present, to identify each area where unique soil model paramaters apply.  The user specifies a line of model input parameters for each zone.  A single soil model is applied to the entire grid, only model parameters are allowed to vary across the grid.

Optional smoothing
------------------
A simple low-pass filter routine has been added to smooth topographic data before computing soil depth and (or) to smooth computed soil depth output from any of the models.  Pelletier and Rasmussen (2009) stated that smoothing topographic data was necessary for their NASD and NSD models.  Filtering soil depth produced by the WNDX model can mitigate abrupt changes in soil depth near channel junctions where the upslope contributing area (flow accumulation) increases abruptly.

The filter algorithm computes the mean value of an N X N (where N is any odd, positive integer) patch of grid cells and replaces the original value at the center cell with the mean. In other words the filter is a running average across the grid applied successively in east-west and north south directions, respectively.  Along edges, at corners, and irregular boundaries the algorithm uses a subset of the N x N patch. Near irregular boundaries, no-data values are excluded from computation of the mean by reflecting interior values across the boundary.  The filter is applied four times to approximate the gaussian filter ﻿(Smith, 1997, see www.dspguide.com/ch24/3.htm).

A line at the end of initialization file, *rg_in.txt*, allows the user to specify whether to smooth the input elevation grid, or computed soil-depth grid by typing `.true.` or `.false.` for each option. The user also specifies the value of N at the end of this line.  Note that topographic smoothing is applied before computing soil depth or any of the intermediate arrays used to compute soil depth.  Soil smoothing is computed after soil depth is computed from the original (`topoSmooth=.false.`) or smoothed (`topoSmooth.true.`) DEM. 

Sample data
-----------

An example rg_in.txt file is provided in the main directory of this repository.  It is configured to apply test-mode parameters to the NSD soil depth model.  Sample data for synthetic terrain are provided in the directory data.  The synthetic terrain is modeled on orthogonal sine waves to represent topography having concave and convex features. Output of the test mode can be compared with analytically computed values of soil depth to confirm that the program is working correctly.

Additional sample data for a small basin in dissected topography of the Oregon Coast Range are available for free download (Baum and others, 2020).  These data are stored in GEOTIFF format, but can be readily converted to ASCII grid format supported by REGOLITH using commercial or open-source GIS software.

References cited
----------------
Baum, R.L., Godt, J.W., and Coe, J.A., 2011, Assessing susceptibility and timing of shallow landslide and debris flow initiation in the Oregon Coast Range, USA, In Genevois, R. Hamilton, D.L. and Prestininzi, A. (eds.) Proceedings of the Fifth International Conference on Debris Flow Hazards Mitigation—Mechanics, Prediction, and Assessment, Padua, Italy, June 7-11, 2011, p. 825-834. Rome: Casa Editrice Universitá La Sapienza (doi: 10.4408/IJEGE.2011-03.B-090).

Baum, R.L., Lewis, A.C., Coe, J.A., and Godt, J.W., 2020, Map and model input and output data for the north Charlotte Creek Basin, Douglas County, Oregon, for analysis of debris-flow initiation resulting from the storm of November 17 - 19, 1996: U.S. Geological Survey data release, https://doi.org/10.5066/P9QLVF5R.

DeRose, R.C.; Trustrum, N.A.; and Blaschke, P.M., 1991, Geomorphic change implied by regolith-slope relationships on steepland hillslopes, Taranaki, New Zealand. Catena, Vol. 18, pp. 489-514.

DeRose, R.C., 1996, Relationships between slope morphology, regolith depth, and the incidence of shallow landslides in eastern Taranaki hill country: Zeitschrift für Geomorphologie, Supplementband 105, pp. 49-60.

Ho, J.-Y.; Lee, K.T.; Chang, T.-C.; Wang, Z.-Y.; and Liao,
Y.-H., 2012, Influences of spatial distribution of soil thickness on shallow landslide prediction: Engineering Geology, Vol. 124, pp. 38–46.

Pelletier, J.D., and Rasmussen, C., 2009, Geomorphically based predictive mapping of soil thickness in upland watersheds: Water Resources Research, Vol. 45, W09417

Pelletier, J.D., Broxton, P.D., Hazenberg, P., Zeng, X., Troch,  P. A., Niu, G.-Y., Williams, Z., Brunke, M. A., and Gochis,  D., 2016, A gridded global data set of soil, immobile regolith, and sedimentary deposit thicknesses for regional and global land surface modeling, J. Adv. Model. Earth Syst., 8, 41–65, doi:10.1002/2015MS000526.

Smith, S.W., 1997, The Scientist and Engineer's Guide to Digital Signal Processing: www.DSPguide.com


