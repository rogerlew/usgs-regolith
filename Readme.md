README
======

Description
-----------
The computer program REGOLITH is a Fortran program designed for estimating soil thickness over a digital landscape.  The program implements various published empirical and computational models for soil thickness.  This command-line program is used in conjunction with Python scripts to prepare input files and visualize results.

Purpose and limitations
-----------------------
The program is designed to make rapid estimates of soil depth using several published models to to create input for spatially distributed slope stability models. The formulas are applied throughout the model domain to both hillslopes and valleys.  The program does not contain any code for statistically testing soil depth distributions against observed or measured soil depths.  

Compiling the program from source code
--------------------------------------
A simple makefile in the src directory can be used to compile the fortran source code.  The makefile has been tested on a linux system which had the GNU compiler collection (gcc) and gfortran installed.  The makefile should be usable after minor modifications with gfortran under Cygwin or OSX.  After downloading the repository, change directories to "src" and type `make` to compile REGOLITH to create an executable file.  If users do not have a compiler installed, an executable file (regolith.exe) is provided in the Releases section of the repository.

Using the program
-----------------
The program runs from the command line and can be compiled for any operating system that supports Fortran 95.  The user provides a small input file, *rg_in.txt*, that contains a list of model parameters and a list of path names for other input files.  The digital elevation model and related input files, slope, upslope contributing area, elevation index, plan-view curvature, and property zone are raster grids in ASCII grid format.  All input grids must have the same spatial reference, resolution, and footprint.  Input file requirements vary depending on the selected soil model; however, each model requires a user provided DEM. The NDSD model requires an elevation index grid, which can be created by the utility program TopoIndex (distributed with TRIGRS).  An updated version of TopoIndex is available at https://code.usgs.gov/usgs/landslides-trigrs.  For the DRS3 model, a user-supplied plan-view curvature grid, which can be created using GIS software, is also needed. The NASD, NSDA, and WNDX models require a flow-accumulation grid (count of upslope grid cells) which can additionally be created using GIS software.  While each model will utilize a grid of slope angle values, the program can compute this and optionally output the grid if one does not exist or if the user specifies topographic smoothing.  The outputted slope grid can be used in future iterations within the same study area to improve computation time.  In addition to a grid of slope angle values, the program outputs a log file, an ASCII grid of computed soil depth, and optional grids of various intermediate values.

Input parameter definitions
---------------------------
The following list explains the variable names, types, and permissible values of the input paramters in order of their appearance in the input file, *rg_in.txt*. 

Name,         Type,      Description
----------------------------------
- `title`        char,      Brief description of project, up to 224 characters
- `num_zones`    integer    Number of parameter zones to allow for variable soils across a landscape
- `max_zones`    integer    Highest zone number, to allow for using a subsample of a larger grid, where some zones of the larger grid are absent in the subsample.
- `h0`           float,     Characteristic soil depth, typically 0.3 - 0.5 m
- `theta_c`      float,     Angle of stability (in degrees) This parameter enters calculation of soil depth in the nonlinear process-based models.  A soil depth of zero or min_depth is assigned to slopes steeper than this angle for all models.
- `dif_ratio`    float,     Diffusivity ratio, >0., (&rho;<sub>b</sub> x P<sub>0</sub>)/(&rho;<sub>s</sub> x D), where D is soil diffusivity, &rho;<sub>b</sub> and &rho;<sub>s</sub> are density of bedrock and soil, respectively, and P<sub>0</sub> is maximum bedrock lowering rate on a flat surface.
- `depth_max`    float,     Maximum soil depth, >0.
- `depth_min`    float,     Minimum soil depth, >=0.
- `C0`           float,     Empirical constant, >0.
- `C1`           float,     Empirical constant, >0.
- `C2`           float,     Empirical constant, >0.
- `chan_thresh`  float,     Threshold (minumum) upslope contributing area at channel heads, >0.
- `chan_depth`   float,     Average steady-state depth of alluvium in steep channels (>=0.2*sc).
- `num_steps`    integer,   Number of (~1-cm) depth increments per unit depth (~100 - 500 for depth in meters)
- `trans_model`  char,      Alpha-numeric code to designate soil depth or transport model 
- Legal values of `trans_model`: DRS1, DRS2, DRS3, WNDX, LCSD, NSD, NSDA, NASD, or NDSD
- `hump_prod`    logical,   Use humped soil production function? (enabled if `.true.` or disabled if `.false.`)
- `power`        integer,   Exponent of DRS2 ploynomial or of uplsope area in NASD, NSDA, and WNDX models
- `elevfil`      char,      File path name of digital elevation grid (<=255 characters)
- `slopefil`     char,      File name of slope angle grid (in degrees) (<=255 characters)
- `flo_accfil`   char,      File name of flow-accumulation grid in units of number of cells, Arc convention (<=255 characters)
- `ndxfil`       char,      File name of elevation index grid (<=255 characters)
- `pv_curvfil`   char,      File name of plan-view curvature grid in units of curvature = curvature * -100, Arc convention (<=255 characters) 
- `zonfil`       char,      File name of parameter-zone grid (<=255 characters)
- `folder`       char,      Folder where output grid files will be stored (<=255 characters)
- `suffix`       char,      Identification code to be added to names of output files (<=16 char.)
- `lasc`         logical,   Specify grid file extension? (*.asc* if `.true.`, *.txt* if `.false`)
- `l_deriv`      logical,   Specify to output derivatives and related quantites
- `l_test`       logical,   Switch between test mode and production mode for the LCSD, NASD, NSDA, and NSD models.  Test mode should be used solely for comparison with analytical solutions.  Production mode is for mapping soil depths using a modified form of the theoretical model and is the preferred mode for applying the models to natural terrain.
- `topoSmooth`   logical, Specify topographic smoothing (enabled if `.true.`,  diabled if `.false.`)
- `soilSmooth`   logical,  Specify smoothing of computed soil depth (enabled if `.true.`,  diabled if `.false.`)
- `n_points`     integer,  Window size, in number of grid cells, of a running average to approximate gaussian smoothing.

Inputs used by each soil depth model
------------------------------------

Model code,     Description and reference(s)
-----------------------------------------------
- DRS1:         Empirical exponential slope dependence (DeRose et al. 1991)
- DRS2:         Empirical 3rd-degree polynomial slope dependence (DeRose et al. 1991)
- DRS3:         Empirical exponential slope and sign of plan-view curvature dependence (Baum et al. 2011)
- WNDX:         Modified soil wetness index, linear area dependent transport (Ho et al. 2012)
- LCSD:         Linear curvature and slope dependent transport (Pelletier & Rasmussen, 2009; Pelletier et al. 2016)
- NSD:          Nonlinear slope dependent transport (Pelletier & Rasmussen, 2009)
- NSDA:         Nonlinear slope dependent transport with linear area dependence (Pelletier & Rasmussen, 2009)
- NASD:         Nonlinear slope and area dependent transport (Pelletier & Rasmussen, 2009)
- NDSD:         Nonlinear slope and depth dependent transport (Pelletier & Rasmussen, 2009)

Soil depth model formulas
------------------------------------
Empirical model formula parameters and REGOLITH input parameter names (if applicable)

- *d*<sub>*r*</sub>,     regolith depth, m
- *C*<sub>*0*</sub>,     empirical constant (C0)
- *C*<sub>*1*</sub>,     empirical constant (C1)
- *C*<sub>*2*</sub>,     empirical constant (C2)
- &theta;,               mean slope angle, degrees
- *p*,                   exponent of DRS2 polynomial or of upslope area (power)
- &kappa;,               plan-view curvature of ground surface, Arc convention
- *A*,                   upslope contributing area, m<sup>2</sup>

Process mased model formula parameters and REGOLITH input parameter names (if applicable)
- *d*<sub>*r*</sub>,        regolith depth, m
- *z*,                      ground surface, m
- *D*,                      diffusivity ratio (dif_ratio)
- *h*<sub>*0*</sub>,        characteristic soil depth, m (h0)
- &theta;,                  mean slope angle, degrees
- *S*<sub>*c*</sub>,        critical slope, tangent of angle of stability, &theta;<sub>*c*</sub>, degrees (theta_c)
- *A*,                      upslope contributing area, m<sup>2</sup>
- *m*,                      exponent of upslope area (power)
- *h*<sub>*n*</sub>,        soil thickness normal to surface, m 

| `trans_model` | Regolith depth formula |
| ------ | ------ |
| DRS1 | ![equation](https://latex.codecogs.com/gif.latex?d_%7Br%7D%3DC_%7B0%7De%5E%7B-C_%7B1%7D%5Ctheta%7D) |
| DRS2 | ![equation](https://latex.codecogs.com/gif.latex?d_%7Br%7D%20%3D%20%28C_%7B0%7D-C_%7B1%7D%5Ctheta%29%5E%7Bp%7D) |
| DRS3 | ![equation](https://latex.codecogs.com/gif.latex?d_%7Br%7D%20%3D%20%28C_%7B0%7D-C_%7B2%7Dsgn%28%5Ckappa%20%29%29e%5E%7B-C_%7B1%7D%5Ctheta%20%7D) |
| WNDX | ![equation](https://latex.codecogs.com/gif.latex?d_%7Br%7D%20%3D%20C_%7B0%7Dln%28%5Cfrac%7BA%5E%7Bp%7D%7D%7Btan%28%5Ctheta%29%7D%29) |
| LCSD | ![equation](https://latex.codecogs.com/gif.latex?d_%7Br%7D%3D%5Cfrac%7Bh_%7B0%7D%7D%7Bcos%28%5Ctheta%29%7Dln%28-%5Cfrac%7BD_%7BLCSD%7D%7D%7Bcos%28%5Ctheta%29%5Ctriangledown%5E%7B2%7Dz%7D%29) |
| NSD | ![equation](https://latex.codecogs.com/gif.latex?d_%7Br%7D%3D%5Cfrac%7Bh_%7B0%7D%7D%7Bcos%28%5Ctheta%29%7Dln%28-%5Cfrac%7BD_%7BNSD%7D%7D%7Bcos%28%5Ctheta%29%5Ctriangledown%5Ccdot%20%5Cfrac%7B%5Ctriangledown%20z%7D%7B1-%28%5Cfrac%7B%7C%5Ctriangledown%20z%7C%7D%7BS_%7Bc%7D%7D%29%5E%7B2%7D%7D%7D%29) |
| NSDA | ![equation](https://latex.codecogs.com/gif.latex?d_%7Br%7D%3D%5Cfrac%7Bh_%7B0%7D%7D%7Bcos%28%5Ctheta%29%7Dln%28-%5Cfrac%7BD_%7BNSDA%7DA%5E%7B-m%7D%7D%7Bcos%28%5Ctheta%29%5Ctriangledown%5Ccdot%20%5Cfrac%7B%5Ctriangledown%20z%7D%7B1-%28%5Cfrac%7B%7C%5Ctriangledown%20z%7C%7D%7BS_%7Bc%7D%7D%29%5E%7B2%7D%7D%7D%29) |
| NASD | ![equation](https://latex.codecogs.com/gif.latex?d_%7Br%7D%3D%5Cfrac%7Bh_%7B0%7D%7D%7Bcos%28%5Ctheta%29%7Dln%28-%5Cfrac%7BD_%7BNASD%7D%7D%7Bcos%28%5Ctheta%29%5Ctriangledown%5Ccdot%20%5Cfrac%7BA%5E%7Bm%7D%5Ctriangledown%20z%7D%7B1-%28%5Cfrac%7B%7C%5Ctriangledown%20z%7C%7D%7BS_%7Bc%7D%7D%29%5E%7B2%7D%7D%7D%29) |
| NDSD | ![equation](https://latex.codecogs.com/gif.latex?d_%7Br%7D%3D%5Cfrac%7Bh_%7B0%7D%7D%7Bcos%28%5Ctheta%29%7Dln%28-%5Cfrac%7BD_%7BNDSD%7D%7D%7Bcos%28%5Ctheta%29%5Ctriangledown%5Ccdot%20%5Cfrac%7Bh_%7Bn%7D%5Ctriangledown%20z%7D%7B1-%28%5Cfrac%7B%7C%5Ctriangledown%20z%7C%7D%7BS_%7Bc%7D%7D%29%5E%7B2%7D%7D%7D%29) |
 

Model code,    Required input parameters
-----------------------------------------------
- DRS1:         `theta_c, depth_max, depth_min, C0, C1, chan_thresh, chan_depth`
- DRS2:         `theta_c, depth_min, depth_max, C0, C1, power, chan_thresh, chan_depth`
- DRS3:         `theta_c, depth_min, depth_max, depth_min, C0, C1, C2, chan_thresh, chan_depth`
- WNDX:         `theta_c, depth_min, depth_max, C0, power, chan_thresh, chan_depth`
- LCSD:          `h0, theta_c, dif_ratio, depth_max, depth_min, chan_thresh, chan_depth, hump_prod, l_test`
- NSD:          `h0, theta_c, dif_ratio, depth_max, depth_min, chan_thresh, chan_depth, hump_prod, l_test`
- NSDA:         `h0, theta_c, dif_ratio, depth_max, depth_min, chan_thresh, chan_depth, hump_prod, l_test`
- NASD:         `h0, theta_c, dif_ratio, depth_max, depth_min, chan_thresh, chan_depth, hump_prod, l_test`
- NDSD:         `h0, theta_c, dif_ratio, depth_max, depth_min, chan_thresh, chan_depth, hump_prod, num_steps`

Model code,    Required input files
------------------------------------------
- DRS1:          `elevfil, slopefil `
- DRS2:          `elevfil, slopefil `
- DRS3:          `elevfil, slopefil, pv_curvfil`
- WNDX:          `elevfil, slopefil, flo_accfil`
- LCSD:          `elevfil, slopefil `
- NSD:           `elevfil, slopefil `
- NSDA:          `elevfil, slopefil, flo_accfil`
- NASD:          `elevfil, slopefil, flo_accfil`
- NDSD:          `elevfil, slopefil, ndxfil`

Input file variable name,   Contents
-------------------------------------
- `elevfil`        A gridded digital elevation model of type float 
- `slopefil`       Slope angle grid of type float
- `flo_accfil `    A grid of upslope contributing area of type float in units of number of cells
- `ndxfil`         A grid of integer values ranking grid cells in the digital elevation model from highest (1) to lowest (maximum number of data cells in the grid).  This grid is used only by the NDSD model to ensure that soil depth is computed in order from peaks and ridge crests down the slope to successively lower points as explained by Pelletier and Rasmussen (2009).  
- ` pv_curvfil`    Plan-curvature grid of type float used only by the DRS3 model
- `zonfil`         A grid of integer values identifying zones within the model where similar values of input parameters are to be applied.

*Notes:* 
1. If `slopefil` does not exist or is not available, REGOLITH will generate a slope grid.
2. The word "none" can be used as a placeholder in *rg_in.txt* for names of files not needed for a particular soil model.

Output files
------------
The following table outlines the prefixes of the output files created within the program. Several intermediates, such as derivatives, can be output if logical `l_deriv` is set to .true. within the input file and outputs will vary based on the model. Additional suffixes may be added to output file names if certain logicals, such as `topoSmooth`, are set to .true..

| File prefix/suffix | Logical to trigger output | Description | Applicable models |
| ------ | ------ | ------ | ------ |
| `RG_Log_`, `RG_` | All model runs | Log file outlining program performance and soil depth output grid | All models |
| `RG_aspect_gs_`, `RG_slope_` | `l_deriv` .true., no slope file specified in input | Grids of angles of steepest descent and slope in degrees | All models |
| `RG_dzdxgs_`, `RG_dzdygs_` | `l_deriv` .true. | Grids of first derivatives in x or y direction | All models |
| `RG_d2zdx2gs_`, `RG_d2zdy2gs_`, `RG_del2gs_`  | `l_deriv` .true. | Grids of second derivatives in x or y direction and Laplacian | LCSD, NDSD |
| `RG_sec_delta_`, `RG_mag_del_z_sq_` | `l_deriv` .true. | Grids of secant delta (1/slope angle) and squared magnitude of z (elevation) | Process-based |
| `RG_nl_slope_fac_`, `RG_trans_x_`, `RG_trans_y_` | `l_deriv` .true. | Grids of slope factor and x or y component of transport factor | NSD, NSDA, NASD, NDSD |
| `RG_trans_(trans_model)_` | `l_deriv` .true. | Grid of transport factor | NSD, NSDA, NASD |
| `RG_Del_dotDelz_nlso_` | `l_deriv` .true. | Grid of derivative of dot product of derivative of z | NDSD |
| `_smo` | `topoSmooth`, `soilSmooth` .true. | Suffix added on smoothed grids | All models |
| `_hmp` | `hump_prod` .true. | Suffix added to output if humped production model was used | All models |
| `_test` | `l_test` .true. | Suffix added to output if model run in test mode | All models |


Suggested input parameter values for empirical models
-----------------------------------------------
The following table displays the full ranges of parameters used in running the program successfully with the empirical models based on results yielded from a variety of study areas with varying geological and climate conditions.  The subsequent tables display site-specific parameter ranges from within these study areas.  The ranges at continental glacial deposits in humid temperate settings, granitoid and gneiss in semi-arid and subalpine settings, and clastic sedimentary geology in humid temperate settings were determined at sites in Mukilteo, WA, Raymond, CO, and North Charlotte Creek, OR, respectively.  The parameters from submarine basalt and volcaniclastic in humid tropical settings and granitoid in human tropical settings were gathered from sites in Anasco, Lares, and Naranjito, Puerto Rico, and Utado, Puerto Rico, respectively, as provided in Tello (2020).  Additionally, parameters from a study area in Eastern Taranaki hill country (sandstone in humid temperate setting) in the North Island of New Zealand as gathered in DeRose et al. (1991) and DeRose (1996) and from a study area in the Tung-An watershed in southern Taiwan (sandstone and shale in a marine tropical to humid temperate setting) as gathered in Ho et al. (2012) are provided.  Note that depth_max and depth_min have been omitted from the subsequent tables as these will vary regardless of the climate and geology, however, a value still must be provided for these parameters when running the program with the empirical models.

| `trans_model` | `theta_c` (degrees) |`depth_max` (m) | `depth_min` (m) | `C0` | `C1` |  `C2` | `power` |
| ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ |
| DRS1 | 30 - 60 | 2 - 5 | 0 - 0.1 | 1 - 7 | 0.005 - 0.06 | — | — |
| DRS2 | 30 - 60 | 2 - 5 | 0 - 0.1 | 1 - 3 | 0.01 - 0.06 | — | 1 - 3 |
| DRS3 | 30 - 60 | 2 - 5 | 0 - 0.1 | 1 - 5 | 0.01 - 0.06 | 1.5 | — |
| WNDX | 30 - 60 | 2 - 5 | 0 - 0.1 | 0.1 - 1 | — | — | 1 - 2 |


Continental glacial deposits in humid temperate setting
| `trans_model` | `theta_c` (degrees) | `C0` | `C1` | `C2` | `power` |
| ------ | ------ | ------ | ------ | ------ | ------ |
| DRS1 | 40 - 60 | 1 - 3 | 0.01 - 0.03 | — | — | 
| DRS2 | 40 - 60 | 1 - 3 | 0.01 - 0.03 | — | 1 - 2 |
| DRS3 | 40 - 60 | 1 - 3 | 0.01 - 0.03 | 1 - 2 | — |
| WNDX | 40 - 60 | 0.2 - 0.5 | — | — | 2 |

Granitoid and gneiss in semi-arid subalpine setting
| `trans_model` | `theta_c` (degrees) | `C0` | `C1` | `C2` | `power` |
| ------ | ------ | ------ | ------ | ------ | ------ |
| DRS1 | 35 - 55 | 1 - 3 | 0.04 - 0.06 | — | — | 
| DRS2 | 35 - 55 | 1 - 3 | 0.05 - 0.07 | — | 3 |
| DRS3 | 35 - 55 | 1 - 3 | 0.04 - 0.06 | 1 | — |
| WNDX | 35 - 55 | 0.1 - 0.3 | — | — | 1 |

Clastic sedimentary in humid temperate setting
| `trans_model` | `theta_c` (degrees) | `C0` | `C1` | `C2` | `power` |
| ------ | ------ | ------ | ------ | ------ | ------ |
| DRS1 | 40 - 60 | 6 - 7 | 0.04 - 0.05 | — | — | 
| DRS2 | 40 - 60 | 1 - 3 | 0.01 - 0.03 | — | 3 |
| DRS3 | 40 - 60 | 1 - 3 | 0.01 - 0.03 | 1 - 2 | — |
| WNDX | 40 - 60 | 0.2 - 1 | — | — | 1 |

Submarine basalt and volcaniclastic in humid tropical setting, Tello (2020)
| `trans_model` | `theta_c` (degrees) | `C0` | `C1` | `C2` | `power` |
| ------ | ------ | ------ | ------ | ------ | ------ |
| DRS1 | 50 - 60 | 1 - 3 | 0.005 - 0.01 | — | — | 
| DRS2 | 50 - 60 | 2 - 3 | 0.02 - 0.04 | — | 3 |
| DRS3 | 50 - 60 | 1 - 3 | 0.01 | 1 - 4 | — |
| WNDX | 50 - 60 | 0.35 - 0.45 | — | — | 1 |

Granitoid in humid tropical setting, Tello (2020) 
| `trans_model` | `theta_c` (degrees) | `C0` | `C1` | `C2` | `power` |
| ------ | ------ | ------ | ------ | ------ | ------ |
| DRS1 | 50 - 60 | 1 - 3 | 0.02 - 0.04 | — | — | 
| DRS2 | 50 - 60 | 1 - 3 | 0.01 - 0.03 | — | 3 |
| WNDX | 50 - 60 | 0.1 - 0.3 | — | — | 1 |

Sandstone in humid temperate setting, DeRose (1996) and Baum et al. (2011)
| `trans_model` | `theta_c` (degrees) | `C0` | `C1` | `C2` | `power` |
| ------ | ------ | ------ | ------ | ------ | ------ |
| DRS1 | 30 - 60 | 5.0 | 0.124 | — | — | 
| DRS2 | 30 - 60 | 1.57 - 1.68 | 0.019 - 0.022 | — | 3 |
| DRS3 | 30 - 60 | 5.0 | 0.04 | 1.5 | — |

Clastic sedimentary in a marine tropical to humid temperate setting, Ho et al. (2012)
| `trans_model` | `theta_c` (degrees) | `C0` | `C1` | `C2` | `power` |
| ------ | ------ | ------ | ------ | ------ | ------ |
| WNDX | 30 - 60 | 0.1- 0.3 | — | — | 1 |
 
 
 
Suggested input parameter values for process-based models
-----------------------------------------------
The following table displays the full ranges of parameters used in running the program successfully with the process-based models based on results yielded from a variety of study areas with varying geological and climate conditions.  The subsequent tables display site-specific parameter ranges from within these study areas.  The ranges at continental glacial deposits in humid temperate settings, granitoid and gneiss in semi-arid and subalpine settings, and clastic sedimentary geology in humid temperate settings were conducted on sites in Mukilteo, WA, Raymond, CO, and North Charlotte Creek, OR, respectively.  The parameters from submarine basalt and volcaniclastic in humid tropical settings and granitoid in human tropical settings were gathered from sites in Anasco, Lares, and Naranjito, Puerto Rico, and Utado, Puerto Rico, respectively, as provided in Tello (2020). `l_test` is available for verifying code output against analytical solutions. The test mode ranges were obtained from Pelletier & Rasmussen (2009) from tests on soils in Pima County, Arizona.

| `trans_model` | `h0` (m) | `theta_c` (degrees) | `dif_ratio` | `depth_max` (m) | `depth_min` (m) | `power` |
| ----- | ------ | ------ | ------ | ------ | ------ | ------- |
| LCSD | 0.3 - 0.5 | 30 - 60 | 0.001 - 0.1 | 2 - 5 | 0 - 0.1 | — |
| NSD | 0.3 - 0.5 | 30 - 60 | 0.005 - 0.1 | 2 - 5 | 0 - 0.1 | — |
| NSDA | 0.3 - 0.5 | 30 - 60 | 0.02 - 5 | 2 - 5 | 0 - 0.1 | 1 - 3 |
| NASD | 0.3 - 0.5 | 30 - 60 | 0.15 - 1 | 2 - 5 | 0 - 0.1 | 1 - 3 |
| NDSD | 0.3 - 0.5 | 30 - 60 | 0.1 - 3 | 2 - 5 | 0 - 0.1 | — |


Continental glacial deposits in humid temperate setting
| `trans_model` |`theta_c` (degrees) | `dif_ratio` | `power` |
| ----- | ------ | ------ | ------ |
| LCSD | 40 - 60 | 0.001 - 0.02 | — |
| NSD | 40 - 60 | 0.005 - 0.03 | — |
| NSDA | 40 - 60 | 0.5 - 0.8 | 2 |
| NASD | 40 - 60 | 6 - 9 | 3 |
| NDSD | 40 - 60 | 0.08 - 1 | — |

Granitoid and gneiss in semi-arid subalpine setting
| `trans_model` | `theta_c` (degrees) | `dif_ratio` | `power` |
| ----- | ------ | ------ | ------ |
| LCSD | 35 - 55 | 0.01 - 0.04 | — |
| NSD | 35 - 55 | 0.04 - 0.07 | — |
| NSDA | 35 - 55 | 1 - 2 | 1 |
| NASD | 35 - 55 | 3 - 4 | 1 |
| NDSD | 35 - 55 | 0.03 - 0.05 | — |

Clastic sedimentary in humid temperate setting
| `trans_model` | `theta_c` (degrees) | `dif_ratio` | `power` |
| ----- | ------ | ------ | ------ |
| LCSD | 40 - 60 | 0.001 - 0.005 | — |
| NSD | 40 - 60 | 0.005 - 0.01 | — |
| NSDA | 40 - 60 | 4 - 5 | 3 |
| NASD | 40 - 60 | 4 - 5 | 2 |
| NDSD | 40 - 60 | 0.08 - 0.09 | — |

Submarine basalt and volcaniclastic in humid tropical setting, Tello (2020)
| `trans_model` | `theta_c` (degrees) | `dif_ratio` | `power` |
| ----- | ------ | ------ | ------ |
| NSD | 50 - 60 | 0.005 - 0.009 | — |
| NSDA | 50 - 60 | 0.02 - 0.15 | 3 |
| NASD | 50 - 60 | 0.15 - 0.25 | 3 |
| NDSD | 50 - 60 | 2 - 3 | — |

Granitoid in humid tropical setting, Tello (2020)
| `trans_model` | `theta_c` (degrees) | `dif_ratio` | `power` |
| ----- | ------ | ------ | ------ |
| NSD | 50 - 60 | 0.08 - 0.1 | — |
| NSDA | 50 - 60 | 0.9 | 3 |
| NASD | 50 - 60 | 1 - 2 | 3 |
| NDSD | 50 - 60 | 0.1 - 0.3 | — |

Test mode ranges for metamorphic and granite in temperate semi-arid Setting, Pelletier & Rasmussen (2009)
| `trans_model` | `theta_c` (degrees) | `dif_ratio` `l_test` | `power` |
| ----- | ------ | ------ | ------ |
| LCSD | 30 - 60 | 1 - 5 | — |
| NSD | 30 - 60 | 0.2 - 1 | — |
| NSDA | 30 - 60 | 1 - 5 | 1 - 2 |
| NASD | 30 - 60 | 1 - 5 | 1 - 2 |
 
Suggested input parameter values for all models
-----------------------------------------------
| Parameter | Suggested range |
| ------ | ------ |
| `chan_thresh` (m<sup>2</sup>) | 1500 - 2000 |
| `chan_depth` (m) | < 0.5 |

The `chan_thresh` parameter represents the minimum upslope contributing area at which channels initiate.  If the upslope contributing area as provided in the flow-accumulation grid is greater than this threshold and if the calculated soil depth is greater than the average depth of alluvium in channels steeper than 20% of the critical slope (`chan_depth`), the depths in these regions will be corrected to `chan_depth`.  To override (turn off) the adjustment of soil depth in steep channels, the user can enter a large value for the channel threshold. In most cases, `chan_thresh` > 10<sup>6</sup> will suffice.  Alternately, for models that are not area dependent, omitting the flow-accumulation grid will likewise override the depth adjustment in steep channels.
 
Optional property zones
-----------------------
REGOLITH has an option to read in a property zone grid for the area covered by the DEM.  This option will allow for application of a soil model over an area where bedrock type, microclimate, or other factors result in significant differences in soil development from one zone to another.  The option requires import of a zone grid, with integer values starting at 1 and ranging up to the number of zones present, to identify each area where unique soil model paramaters apply.  The user specifies a line of model input parameters for each zone.  A single soil model is applied to the entire grid, only model parameters are allowed to vary across the grid.

Optional smoothing
------------------
A simple low-pass filter routine has been added to smooth topographic data before computing soil depth and (or) to smooth computed soil depth output from any of the models.  Pelletier and Rasmussen (2009) stated that smoothing topographic data was necessary for their NASD and NSD models.  Filtering soil depth produced by the WNDX model can mitigate abrupt changes in soil depth near channel junctions where the upslope contributing area (flow accumulation) increases abruptly.

The filter algorithm computes the mean value of an N X N (where N is any odd, positive integer) patch of grid cells and replaces the original value at the center cell with the mean. In other words the filter is a running average across the grid applied successively in east-west and north south directions, respectively.  Along edges, at corners, and irregular boundaries the algorithm uses a subset of the N x N patch. Near irregular boundaries, no-data values are excluded from computation of the mean by reflecting interior values across the boundary.  The filter is applied four times to approximate the gaussian filter �(Smith, 1997, see www.dspguide.com/ch24/3.htm).

A line at the end of initialization file, *rg_in.txt*, allows the user to specify whether to smooth the input elevation grid, or computed soil-depth grid by typing `.true.` or `.false.` for each option. The user also specifies the value of N at the end of this line.  Note that topographic smoothing is applied before computing soil depth or any of the intermediate arrays used to compute soil depth.  Soil smoothing is computed after soil depth is computed from the original (`topoSmooth=.false.`) or smoothed (`topoSmooth.true.`) DEM. 

Sample data
-----------

An example rg_in.txt file is provided in the main directory of this repository.  It is configured to apply test-mode parameters to the NSD soil depth model.  Sample data for synthetic terrain are provided in the directory data.  The synthetic terrain is modeled on orthogonal sine waves to represent topography having concave and convex features. Output of the test mode can be compared with analytically computed values of soil depth to confirm that the program is working correctly.

Additional sample data for a small basin in dissected topography of the Oregon Coast Range are available for free download (Baum and others, 2020).  These data are stored in GEOTIFF format, but can be readily converted to ASCII grid format supported by REGOLITH using commercial or open-source GIS software.

Regolith iterations
-------------------
A Jupyter Notebook script Regolith_Iterations.ipynb is provided in the main directory of this repository and can be used in conjunction with the program to test varying ranges of parameters.  Running the script will iterate through these parameters by rewriting rg_in.txt with the varying parameters within the ranges specified by the user within the script.  REGOLITH will run for each set of parameters and output results with varying suffixes to denote the parameters used in each run.


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

Tello, M., 2020, Optimization of landslide susceptibility modeling: a Puerto Rico case study, 2020 - Mines Theses and Dissertations.


