README - Examples
==============

Description
-----------
This directory contains three input files and input data to display the differences in running the program with empirical versus process-based models as well as an example of running the program using data with multiple zones with varying parameters. The input data comes from a study area in North Charlotte Creek, Oregon. Parameters used for running the program with varying models were tested using the Regolith_iterations.ipynb notebook. The parameters yielding soil depth outputs most similar to the known characteristics of this study area were selected as the best fit parameters to this study area.

The input file "rg_inDRS3.txt" and "rg_inNASD.txt" display examples of input files for empirical and process-based models, respectively. The primary difference in these files is in line 12 of the input files. For empirical models, the parameters listed in this line are sc, depth_min, depth_max, C0, C1, and C2 while those for process-based models are sc, depth_min, depth_max, h0, dif_ratio, and hump_prod. Test mode, as specified in line 4 of the input file, can only be applied to process-based models, with the exception of NDSD. If test mode is specified as .true. when running the program with an empirical model or NDSD, this will have no effect on model performance. Additionally, the exponent of the DRS2 polynomial or of upslope contributing area in the NASD, NSDA, and WNDX models, as specified in line 15 of the input file, will not affect the performance of models that do not require this parameter when specified in the input file.

For study areas with zones of varying parameters, a model-parameter zone grid, with cells specifying the parameter zone, can be specified in the input file. An example of running the program with zones for the same study zone but using two parameter zones is displayed in "rg_inWNDXzones.txt." Several changes must be made to the input file when including parameter zones in order for it to be compatible with the program. First, a model-parameter grid ("zones.asc") must be included in the input file. The number of soil depth zones and the maximum number of zones (num_zones, max_zones) in line 8 of the input file must be changed to agree with the number of zones in the model-parameter grid. New parameter lines for each of the zones are then added until the number of parameter zones agrees with the specified number of zones, n, as follows:

zone, 1  
&#35;sc, depth_min, depth_max, C0, C1, C2 for empirical models or sc, depth_min, depth_max, h0, dif_ratio, hump_prod for process-based models  
55,	0.1,	5.0,	0.7,	1,   1  
zone, 2  
&#35;sc, depth_min, depth_max, C0, C1, C2 for empirical models or sc, depth_min, depth_max, h0, dif_ratio, hump_prod for process-based models  
55,	0.1,	5.0,	0.35,	1,   1  
...  
zone, n  
&#35;sc, depth_min, depth_max, C0, C1, C2 for empirical models or sc, depth_min, depth_max, h0, dif_ratio, hump_prod for process-based models  
55,	0.1,	5.0,	0.5,	1,   1  


