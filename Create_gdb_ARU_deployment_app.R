#This script creates a gdb file for the ARU deployment app for any GBMP sites that are being sampled with ARUs
#This example is fo SGOV

#load existing GDB or gpkg file created when GRTS sample was drawn

#Govenlock will have some odd locations removed. We've switched to doing only 8 ARUs per section, including 1, 3, 6, 8, 9, 11, 14, 16.
#But there are a few odd exceptions for Govenlock due to past sampling and site boundary
sf <- st_read(dsn = "C:/GBM2026.gdb",
              layer = "GBM2026") |>
  mutate(stn = as.numeric(sub("^[^-]+-[^-]+-", "", station)),
         grid = substr(station, 1, 7)) |>
  arrange(station) |>
  #Remove northern 8 points in section 3 because they are outside of ranch boundary
  filter(!(grid == "SGOV-03" & stn %in% c(1:8))) |>
  #remove sourthern 8 points from grid 9 because they have never been sampled in the past (possibly access issues)
  filter(!(grid == "SGOV-09" & stn %in% c(9:16))) |>
  #Remove eastern locations in grid 10
  filter(!(grid == "SGOV-10" & stn %in% c(1,2,5,6,9,10,13,14))) |>
  #for remaining sites, remove the normal numbers 
  filter(grid %in% c("SGOV-03", "SGOV-09", "SGOV-10") | (!(grid %in% c("SGOV-03", "SGOV-09")) & stn %in% c(1,3,6,8,9,11,14,16)))

#plot to make sure it worked
plot(select(sf, grid))

#export
st_write(sf, "C:/Users/robinsonba/OneDrive - EC-EC/Documents/Projects/Avian Trend Monitoring/Grassland Bird Monitoring/Data/GIS files/GBMP Sites/SNAS_SBCR_SGOV/SGOV_ARUapp.gdb")
plot(sf)
