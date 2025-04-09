### INSTRUCTIONS
# Move your site boundary shapefile into the Data > Boundaries folder
# Input the arguments into the function at the bottom of the page:
# # shapefile: If only a single shapefile defines the site, enter the name of the shapefile. If multiple shapefiles, enter a list with all shapefile names
# site_ID: enter the four letter site code
# pcode: project code
# sample_size: number of desired base sample grids, if entered will override calculated 20%
# overdraw_size:  number of desired overdraw grids, if entered will override calculated 1/4 of sample size
# check the interactive map output (the first element of the list that is output from the function) 
# to see that everything looks normal and there is not too much clumping
# rerun the function if necessary
# NOTE: delete the entire output folder before rerunning to ensure everything is properly exported

source("Functions/GRTS_Function.R")
newSite <- grts_GBMP(shapefile = "StrattonPeake_Dorothy",
                  site_ID = "ASPD",
                  pcode = "GBM")

#view interavtive map to look at the sample
newSite[1]
