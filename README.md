## MOQ - Mountain Observation QC

#### Scripts (MATLAB) for running the quality control methods developed explicitly for mountain observations of irradiances. 
_Specifically targets:_

item 1 Burial of sensors by snow (MOQ_SOD_Detect.m, used in SW_OBS_QC_v2.m and SnowOnDome_StandAlone.m, stand alone version of function)

item 2 Shading (SW_Obs_QC_PlottingTool_SolGeo.m)

_Additionally:_

-other climate based checks (MOQ_LongShiLimits_Graph.m, based on Long and Shi 2008, An Automated Quality Assessment and Control Algorithm for Surface Radiation Measurements)

-master function (SW_Obs_QC_v2.m)

Notes on dependencies: 
item 1 Plotting transmissivity requires [colorbrewer](http://www.mathworks.com/matlabcentral/fileexchange/34087-cbrewer---colorbrewer-schemes-for-matlab)
functions

Item 2 Assumes times/dates are in [time_builder](github.com/klapo/time_tools) format
