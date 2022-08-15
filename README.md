# Badas_Methane
Mathematical model (in Matlab) that simulates methane and carbon dioxide dynamics in a tropical peatland drainage canal. These codes are shared as part of the submission process of an academic journal article. Accompanying field data for these codes is available from Hydroshare: 

The inital commit contains three matlab scripts:

<b> Along_canal_nlimmultifit_DOM_scaled_Jan_Aug.m </b> 
This file fits the mathematical model of methane and carbon dioxide concentrations and stable carbon isotopes for field data.

<b> Error_Bands_Budget_Jan_Aug_new.m </b>  
This script take the best fit simulation and runs a monte-carlo error propagation routine and then makes budget plots that are in the manuscript. 

<b> GW_exp_fit.m </b>
This is a function called by the above codes that fits the depth of groundwater contribution to the canal based on field data.


