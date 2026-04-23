#calculate centroid UTMs for section grids
library(sf)
library(dplyr)
st_layers("Data/GIS/SD.gpkg")
secSample <- st_read(dsn = "Data/GIS/SD.gpkg", layer = "SD_SecSample")
centroid <- st_centroid(secSample)
coords <- as.data.frame(st_coordinates(centroid))
coords$site <- paste0("NWA-", centroid$rank)
coords$rank <- centroid$rank
coords <- arrange(coords, rank)
write.csv(coords, "Output/SD_secCentroids.csv", row.names = F)
