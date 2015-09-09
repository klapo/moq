## MOQ - Mountain Observation QC

#### Quality control methods in MATLAB (SW_Obs_QC_v2.m) developed explicitly for mountain observations of irradiances. 

#####Specifically targets:

1. Burial of sensors by snow (MOQ_SOD_Detect.m) ![image of snow on dome @ Snoqalmie](https://github.com/klapo/moq/blob/master/TimeLapse.Snoqualmie.SnowOnDome.png)
2. Shading (SW_Obs_QC_PlottingTool_SolGeo.m)

#####Additionally:

-Includes climate based checks (MOQ_LongShiLimits_Graph.m) based on QCRad.
[Long and Shi 2008, An Automated Quality Assessment and Control Algorithm for Surface Radiation Measurements](http://www.arm.gov/publications/tech_reports/doe-sc-arm-tr-074.pdf)

#####Notes on dependencies:

1. Plotting transmissivity requires [colorbrewer](http://www.mathworks.com/matlabcentral/fileexchange/34087-cbrewer---colorbrewer-schemes-for-matlab)
functions

2. Assumes times/dates are in [time_builder](http://www.github.com/klapo/time_tools) format
