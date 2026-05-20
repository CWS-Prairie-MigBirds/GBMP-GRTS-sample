#This script creates an export table for the Nature Counts point count app for older sites that 
#were created prior to app creation
library(sf)
library(dplyr)

#load ASMB gpkg file for the site with all survey locations
st_layers(dsn = "C:/Users/robinsonba/OneDrive - EC-EC/Documents/Projects/Avian Trend Monitoring/Grassland Bird Monitoring/Data/GIS files/GBMP Sites/ASMB/ASMB.gpkg")
pts <- st_read(dsn = "C:/Users/robinsonba/OneDrive - EC-EC/Documents/Projects/Avian Trend Monitoring/Grassland Bird Monitoring/Data/GIS files/GBMP Sites/ASMB/ASMB.gpkg",
               layer = "ASMB_PointSample")

#This Suffield MTA table is in a really old format so need to change a few things
pts <- pts |>
  arrange(rank) |>
  st_transform(crs = 4269)

#There are 50 sections, but it seems in 2022, we cut the sample down to 25, so query out the appropriate ones
pts <- pts |>
  filter(rank %in% c(1,28,24,27,4,9,5,23,46,15,49,12,8,16,13,25,14,26,45,3,2,17,6,11,10))

#plot to double check things match the map from 2022
pts_rank <- pts |>
  select(rank)
plot(pts_rank)

#extract coordinates
coords <- st_coordinates(pts)
pts <- pts |>
  mutate(PCODE = "GBM",
         SITE = "ASMB",
         grid = sprintf("%02d", rank),
         point = sprintf("%02d", pts),
         lat = coords[,2],
         lon = coords[,1]) |>
  mutate(label = paste(SITE, grid, point, sep = "-"))

nc_asmb <- pts %>% 
  dplyr::select(PCODE, SITE, label, lat, lon) %>% 
  rename(project_id = PCODE, 
         site_id = SITE, 
         point_id = label, 
         latitude = lat, 
         longitude = lon) %>%
  st_drop_geometry()


#Now repeat for Suffield NWA (SD)
st_layers(dsn = "C:/Users/robinsonba/OneDrive - EC-EC/Documents/Projects/Avian Trend Monitoring/Grassland Bird Monitoring/Data/GIS files/GBMP Sites/SD/SD.gpkg")
pts <- st_read(dsn = "C:/Users/robinsonba/OneDrive - EC-EC/Documents/Projects/Avian Trend Monitoring/Grassland Bird Monitoring/Data/GIS files/GBMP Sites/SD/SD.gpkg",
               layer = "SD_SecSample_pts") |>
  st_transform(crs = 4269)

plot(select(pts, rank))

coords <- st_coordinates(pts)
nc_sd <- pts |>
  mutate(project_id = "GBM",
         site_id = "SD",
         latitude = coords[,2],
         longitude = coords[,1]) |>
  mutate(point_id = paste(site_id, STN, sep = "-")) |>
  select(colnames(nc_asmb)) |>
  st_drop_geometry()

nc_all <- rbind(nc_asmb, nc_sd)
#export individual tables for each site to go in the site folders, then a combined table for Catherine
write.csv(nc_asmb, paste0("Output/", "ASMB", "_NatureCounts_Upload.csv"), row.names = F)
write.csv(nc_sd, paste0("Output/", "SD", "_NatureCounts_Upload.csv"), row.names = F)
write.csv(nc_all, paste0("Output/", "GBM_2026", "_NatureCounts_Upload.csv"), row.names = F)

















