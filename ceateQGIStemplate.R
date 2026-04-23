#This script create the gpkg files needed to make site maps in QGIS. This is for older sites
#where GRTS samples were drawn from an ealier version of this script before we starting using 
#QGIS to make mapes

library(sf)

#First copy the gpkg file for the site into the project Data/GIS folder
#update the file path below accordingly
file.copy(from = "C:/Users/robinsonba/OneDrive - EC-EC/Documents/Projects/Avian Trend Monitoring/Grassland Bird Monitoring/Data/GIS files/GBMP Sites/ASMB/ASMB.gpkg",
          to = "Data/GIS/ASMB.gpkg")
#update site code
siteCode <-"ASMB"
#load the GPKG file
st_layers(dsn = "Data/GIS/ASMB.gpkg")
site = st_read(dsn = "Data/GIS/ASMB.gpkg",
               layer = paste0(siteCode, "_boundary"))
