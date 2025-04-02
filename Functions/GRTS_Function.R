
          ### CWS GRASSLAND BIRD MONITORING PROGRAM ###
  ### Generalized Random Tessellation Stratified (GRTS) sample ###
      # Barry Robinson, Janine McManus, Josiah Van Egmond #
                  # last updated: March 2025 #

### INFO =======================================================================
#*This is a custom function to generate a ranked random sample oflegal land sections or quarter sections in Alberta, Saskatchewan or
#*Manitoba and establish point count locations in 4x4 grids spaced 400 m apart.

### Functionality to be added:
    # make pretty map
    # generate sample from only legal land description


#temp input values to test function
shapefile = "HandHillER"
shapefile = "AndersonNorthLease"
shapefile = "dundurn"
shapefile = "hole_in_the_wall"
shapefile = list("AndersonNorthLease", "AndersonSouthLease")
shapefile = "manitou"
site_ID = "AHHL"
grts_GBMP <- function(shapefile,
                      site_ID, 
                      pcode, 
                      sample_size = NULL, 
                      overdraw_size = NULL)  { 
   
  
#Load packages; required packages will be automatically installed and loaded
  source("Functions/InstallLoadPackages.R")
  using(
    "tidyverse",
    "parallel", 
    "sf", 
    #"raster",
    "gfcanalysis", 
    "spsurvey", 
    "pgirmess", 
    "dplyr", 
    "ggplot2", 
    "matrixStats",
    "tmap",
    "tmaptools",
    #"OpenStreetMap",
    "basemaps",
    "geodata"
  )
 
##### FILE SETUP ===============================================================
  
  #create output directory for project
  dir.create(paste0("./Output/", site_ID),showWarnings = F)
    
  #switch off spherical geometry
  sf_use_s2(FALSE)
  
  #load in site shapefile, if a list of shapefiles is provided, combine
  if(is.list(shapefile)) {
    site <- lapply(shapefile, function(x) {return(st_read(dsn = "Data/Boundaries", layer = x))})
    site <- do.call("bind_rows", site)
  } else {
    site <- st_read(dsn = "Data/Boundaries", layer = shapefile)
  }
  
  #some sites have boundary files that have multiple polygons that need to be grouped together
  #This step also simplifies the attribute table to have only 1 field: site name
  site$name <- site_ID
  site <- site %>% group_by(name) %>% dplyr::summarise()
  
  
##### Transform to UTM NAD83 projection to allow calculation of coordinates for point count locations========================================
  
  #if projection is missing, set to WGS84
  crs_info <- st_crs(site)
  if(is.na(crs_info)) {
    st_crs(site) <- 4326
  }
  
  #Determine UTM zone
  #if projection is lat/long
  if(grepl("+proj=longlat", crs_info$proj4string)) {
    centroid <- st_centroid(site)
    centroid <- st_coordinates(centroid)                  #calculate coordinates of centroid
    zone <- utm_zone(centroid[1], centroid[2])            #calculate UTM zone using centroid coordinates
    utm_zone_num <- as.numeric(sub("[C-Z]$", "", zone))   #drop letter from zone, since sites are only in northern hemisphere
  }
  
  #if projection is utm
  if(grepl("+proj=utm", crs_info$proj4string)) {
    utm_zone_num <- as.numeric(str_extract(crs_info$proj4string, "(\\d+)"))
  }
  
  #define projection with NAD83 datum and corresponding UTM Zone, transform if needed
  epsg <- 26900 + utm_zone_num
  crs_srid <- as.numeric(str_extract(crs_info$srid, "(\\d+)"))
  if(crs_srid != epsg) {
    site <- st_transform(site, crs = epsg)
  }
  
  ##BRobinson April 24: Moved this to the end to ensure everything is exported in same projection
  
  # #export site boundary gpkg for later map use
  # 
  # st_write(obj = site, 
  #          dsn = paste0("./output/", site_ID, "/", site_ID, ".gpkg"), 
  #          layer = paste0(site_ID, "_boundary"), 
  #          driver = "GPKG",
  #          append = TRUE)  
  

##### GET PROVINCE GRID FOR SITE ===============================================
  
  # load in Canada province data, filter only prairie region provinces
  provinces <- c("Alberta", "Manitoba", "Saskatchewan")
  
  #check if Canada shapefile exists before downloading
  if(!file.exists("Data/gadm/gadm41_CAN_1_pk.rds")) {
    can <- gadm(country="CAN", level=1, path="Data")
  } else {
    can <- readRDS("Data/gadm/gadm41_CAN_1_pk.rds")
  }
  
  #filter out prairie provinces
  can <- st_as_sf(can)  
  can <- filter(can, NAME_1 %in% provinces)

  #transform Canada province data to site crs and intersect
  can <- st_transform(can, crs = st_crs(site))
  can <- st_intersection(can, site)

  # read in province grid data and change legal land description column name to LLD
  # empty grid object so province grids can be bound to it if site crosses two provinces
  
  grid <- NULL
  
  if ("Alberta" %in% can$NAME_1) {
    grid <- st_read(dsn = "Data/Grid/section_grid_BCR11.gpkg", layer = "AB") %>%
      rename(LLD = DESCRIPTOR)
    grid <- st_transform(grid, crs = st_crs(site))
  }
  
  if ("Saskatchewan" %in% can$NAME_1) {
    grid_tmp <- st_read(dsn = "Data/Grid/section_grid_BCR11.gpkg", layer = "SK")
    grid_tmp <- st_transform(grid_tmp, crs = st_crs(site))
    
    # Merge into the main grid object
    if (is.null(grid)) {
      grid <- grid_tmp
    } else {
      grid <- bind_rows(grid, grid_tmp)
    }
  }
  
  if ("Manitoba" %in% can$NAME_1) {
    grid_tmp <- st_read(dsn = "Data/Grid/section_grid_BCR11.gpkg", layer = "MB") %>%
      rename(LLD = LEG_DESC)
    grid_tmp <- st_transform(grid_tmp, crs = st_crs(site))
    
    # Merge into the main grid object
    if (is.null(grid)) {
      grid <- grid_tmp
    } else {
      grid <- bind_rows(grid, grid_tmp)
    }
  }
  
  #remove temp grid to clear up RAM
  rm(grid_tmp)
  gc()

  
##### PROCESS PROVINCE SECTION GRID ============================================
  
  #Query out all land sections that intersect with each site polygon
  gridToKeep <- st_intersection(x = site, y = grid)
  gridToKeep$area = st_area(gridToKeep)
  
  
  #Remove sections with <half of area within the site
  # 1 mile^2 = 259000 m2, divide by 2 for area of half a section
  #using 2.05 instead of 2 to have some tolerance for areas that are slightly less than half a section due to imperfect polygons and rounding error
  gridToKeep = gridToKeep %>% dplyr::filter(as.numeric(area) >= 2590000/2.05)
  
  #GBMP monitoring design is different for large (>5 sections) and small (<= 5 sections) sites
  #Determine site size
  
  if(nrow(gridToKeep) <= 5) {
    site_type = "small" 
  } else {
    site_type = "large"
  }

#### If site is small, proceed with quarter section grid 
  if (site_type == "small") {
   
   grid <- NULL

   if ("Alberta" %in% can$NAME_1) {
     grid <- st_read(dsn = "Data/Grid/quarter_grid_BCR11.gpkg", layer = "AB") %>%
       rename(LLD = DESCRIPTIO)
     grid <- st_transform(grid, crs = st_crs(site))
   }

   if ("Saskatchewan" %in% can$NAME_1) {
     grid_tmp <- st_read(dsn = "Data/Grid/quarter_grid_BCR11.gpkg", layer = "SK")
     grid_tmp <- grid_tmp %>% unite(col = "LLD", c("QSECT", "PSECT", "PTWP", "PRGE", "PMER"), sep = "-")
     grid_tmp <- st_transform(grid_tmp, crs = st_crs(site))

     #Merge into the main grid object
     if (is.null(grid)) {
       grid <- grid_tmp
     } else {
       grid <- bind_rows(grid, grid_tmp)
     }
   }

   if ("Manitoba" %in% can$NAME_1) {
     grid_tmp <- st_read(dsn = "Data/Grid/quarter_grid_BCR11.gpkg", layer = "MB") %>%
       rename(LLD = LEGAL_DESC) %>% 
       rename(zone_1 = ZONE)
     grid_tmp <- st_transform(grid_tmp, crs = st_crs(site))

     # Merge into the main grid object
     if (is.null(grid)) {
       grid <- grid_tmp
     } else {
       grid <- bind_rows(grid, grid_tmp)
     }
   }
   
   #remove temp grid to clear up RAM
   rm(grid_tmp)
   gc()
   
   #intersect quarter grid with site boundary
   gridToKeep <- st_intersection(x = site, y = grid)
   
   #calculate areas to filter out partial grids
   gridToKeep$area = st_area(gridToKeep)
   
   #Remove quater sections with <half of area within the site
   gridToKeep = gridToKeep %>% dplyr::filter(as.numeric(area) >= 647497/2.05)
  }
  
  ###filter out section/quarter section grids by legal land description
  fullgrid <- grid %>% dplyr::filter(LLD %in% c(gridToKeep$LLD))
  
  centroid <- st_centroid(fullgrid)
  
  fullgrid$CentX <- st_coordinates(centroid, byid = T, id = centroid$LLD)[,1]
  fullgrid$CentY <- st_coordinates(centroid, byid = T, id = centroid$LLD)[,2]
  
  centroid$CentX <- st_coordinates(centroid, byid = T, id = centroid$LLD)[,1]
  centroid$CentY <- st_coordinates(centroid, byid = T, id = centroid$LLD)[,2]

  # #write section grid file  
  # st_write(obj = fullgrid, 
  #          dsn = paste0("./output/", site_ID, "/", site_ID, ".gpkg"),
  #          layer = paste0(site_ID, "_SecGrid"), 
  #          driver = "GPKG",
  #          append = TRUE)
  
  
  rm(grid) #remove grid to save RAM
  gc() #garbage collection to free up RAM
  
  
##### DRAW RANDOM SAMPLE =======================================================
  
  #calculate sample size
  #Large sites: either 20% of available grid cells or 5 cells, whichever is greater
  #Small sites: include all grid cells from fullgrid so that all quarter sections are ranked.
  #or custom input
  #overdraw size is 1/4 of sample size, may be too large for big sites
  #can manually enter overdraw sample size
  
  #BROBINSON: These ifelse statements need to be reviewed relative to the sample size rules we've established. As they are, they are more complex that the simple rule written above
  N <- nrow(gridToKeep)
  if(is.null(sample_size)) {
    if(site_type == "large") {
      if(0.2*N > 5) {
        n <- round(0.2*N)
      } else { #if(0.2*n < 5 & n >5){ #BRobinson: This isn't needed anymore because site size is now defined above
        n <- 5
      }
    } else {
      n = N
    }
  } else {
    n <- sample_size
  }
  
  if(is.null(overdraw_size)) {
    if(site_type == "small") {
      n.over = 0
    } else {
      n.over <- ceiling(0.25*N)
    }
  } else {
    n.over <- overdraw_size
  }
  
  #draw GRTS sample
  sample = grts(sframe = centroid, n_base = n, seltype = "equal", n_over = n.over)
  
  #get samples out of list
  sample_b <- sample$sites_base
  sample_o <- sample$sites_over
  sample <- rbind(sample_b, sample_o)
  
  #Create a rank column
  sample$rank <- 1:nrow(sample)
  
  #create section sample grids
  sec_sample <- dplyr::filter(fullgrid, LLD %in% c(sample$LLD)) 
  
  #add rank column to sec_sample
  rank <- st_drop_geometry(sample[c("LLD", "rank")])
  sec_sample <- left_join(sec_sample, rank, by = "LLD")

  
##### CALCULATE POINT COUNT LOCATIONS ==========================================
  
  # calculate XY coordinates for point count locations within the each grid drawn for the sample
  
  if(site_type =="small") { #small sites have 4 points per quarter section grid
    grid.exp <- st_drop_geometry(sample)
    grid.exp <- grid.exp[rep(seq_len(nrow(grid.exp)), 4),]
    #order by rank
    grid.exp <- grid.exp[order(grid.exp$rank),]
    #add field to indicate point number within each section (1-4)
    grid.exp$pts <- as.vector(rep(1:4, nrow(sample)))
    #create fields for X and Y coordinates of points
    grid.exp$Xcl <- 0
    grid.exp$Ycl <- 0
    
    #Calculate coordinates for each point based on it's distance from the centroid
    grid.exp$Xcl <- ifelse(grid.exp$pts==2 | grid.exp$pts==4,
                           grid.exp$CentX +200, grid.exp$CentX-200)
    grid.exp$Ycl <- ifelse(grid.exp$pts==2 | grid.exp$pts==1,
                           grid.exp$CentY +200, grid.exp$CentY-200)
    
 
    
  } else { #large sites have 16 points per section grid
    grid.exp <- st_drop_geometry(sample)
    grid.exp <- grid.exp[rep(seq_len(nrow(grid.exp)), 16),]
    #order by rank
    grid.exp <- grid.exp[order(grid.exp$rank),]
    #add field to indicate point number within each section (1-16)
    grid.exp$pts <- as.vector(rep(1:16, nrow(sample)))
    #create fields for X and Y coordinates of points
    grid.exp$Xcl <- 0
    grid.exp$Ycl <- 0
    
    #Calculate coordinates for each point based on it's distance from the centroid
    grid.exp$Xcl <- ifelse(grid.exp$pts==1 | grid.exp$pts==5 | grid.exp$pts==9 | grid.exp$pts==13, grid.exp$CentX+600, grid.exp$Xcl)
    grid.exp$Xcl <- ifelse(grid.exp$pts==2 | grid.exp$pts==6 | grid.exp$pts==10 | grid.exp$pts==14, grid.exp$CentX+200, grid.exp$Xcl)
    grid.exp$Xcl <- ifelse(grid.exp$pts==3 | grid.exp$pts==7 | grid.exp$pts==11 | grid.exp$pts==15, grid.exp$CentX-200, grid.exp$Xcl)
    grid.exp$Xcl <- ifelse(grid.exp$pts==4 | grid.exp$pts==8 | grid.exp$pts==12 | grid.exp$pts==16, grid.exp$CentX-600, grid.exp$Xcl)
    grid.exp$Ycl <- ifelse(grid.exp$pts==1 | grid.exp$pts==2 | grid.exp$pts==3 | grid.exp$pts==4, grid.exp$CentY+600, grid.exp$Ycl)
    grid.exp$Ycl <- ifelse(grid.exp$pts==5 | grid.exp$pts==6 | grid.exp$pts==7 | grid.exp$pts==8, grid.exp$CentY+200, grid.exp$Ycl)
    grid.exp$Ycl <- ifelse(grid.exp$pts==9 | grid.exp$pts==10 | grid.exp$pts==11 | grid.exp$pts==12, grid.exp$CentY-200, grid.exp$Ycl)
    grid.exp$Ycl <- ifelse(grid.exp$pts==13 | grid.exp$pts==14 | grid.exp$pts==15 | grid.exp$pts==16, grid.exp$CentY-600, grid.exp$Ycl)
    }
  
  
  
  #Modify table so that it can be brought into GBM Access database and ArcMap. 
  # ***MODIFY SITE, PCODE, Datum, and Zone***
  
  grid.exp$SITE <- site_ID 
  grid.exp$PCODE <- pcode
  # grid.exp$Datum <- "NAD83"   
  # grid.exp$Zone <- utm_zone_num
  grid.exp$STN <- paste(sprintf("%02d", grid.exp$rank),sprintf("%02d",grid.exp$pts),sep="-")
  grid.exp$label <- paste(grid.exp$SITE,grid.exp$STN,sep="-")
  
  #convert to sf object
  pts <- st_as_sf(grid.exp, coords = c("Xcl", "Ycl"))
  
  #update point crs to match site
  st_crs(pts) = crs(site)
  
  #clip out any points that fall outside of original study site
  pts <- st_intersection(pts, site)
  pts <- dplyr::select(pts, label, everything())
 
  
##### EXPORT ALL GEOSPATIAL FILES TO THE GPKG SITE FILE ========================
  #Convert everything back to WGS84 projection
  gpkg_export <- lapply(list(site = site, fullgrid = fullgrid, sec_sample = sec_sample, pts = pts), FUN = st_transform, crs = 4326)
  
  #export site boundary file
  st_write(obj = site, 
           dsn = paste0("./output/", site_ID, "/", site_ID, ".gpkg"), 
           layer = paste0(site_ID, "_boundary"), 
           driver = "GPKG",
           append = TRUE)
  
  #export section grid file  
  st_write(obj = fullgrid, 
           dsn = paste0("./output/", site_ID, "/", site_ID, ".gpkg"),
           layer = paste0(site_ID, "_SecGrid"), 
           driver = "GPKG",
           append = TRUE)
  
  #export section sample file
  st_write(obj = sec_sample, 
           dsn = paste0("./output/", site_ID, "/", site_ID, ".gpkg"), 
           layer = paste0(site_ID, "_SecSample"), 
           driver = "GPKG", 
           append = TRUE
  )
  
  #export the point location file
  st_write(obj = pts, 
           dsn = paste0("./output/", site_ID, "/", site_ID, ".gpkg"), 
           layer = paste0(site_ID, "_PointSample"), 
           driver = "GPKG",
           append = TRUE
  )
  
##### WRITE GPX FILE ===========================================================
  
  #transform points to WGS83 for export to gpx
  latlong <- st_transform(pts, crs = 4326)
  
  #get coordinates from file
  latlong$lon <- st_coordinates(latlong)[,1]
  latlong$lat <- st_coordinates(latlong)[,2]
   
  #filter dataframe for export to gpx
  latlong2 <- latlong %>% 
    dplyr::select(label, lon, lat) %>% 
    rename(name = label) 

   
  ### write_gpx function ###
  # custom function to write GPX file
  #adapted by Josiah from https://rdrr.io/github/s-u/snippets/src/R/gpx.R 
  
   write_gpx <- function(lat, lon, time = NULL, name = NULL, out_file) {
     o <- c('<gpx version="1.1" creator="R">')
     
     if (is.null(time)) {
       for (i in seq_along(lat)) {
         if (is.null(name)) {
           o <- c(o, paste('<wpt lat="', lat[i], '" lon="', lon[i], '" />', sep=''))
         } else {
           o <- c(o, paste('<wpt lat="', lat[i], '" lon="', lon[i], '"><name>', name[i], '</name></wpt>', sep=''))
         }
       }
     } else {
       for (i in seq_along(lat)) {
         if (is.null(name)) {
           o <- c(o, paste('<wpt lat="', lat[i], '" lon="', lon[i], '"><time>', 
                           paste(gsub(' ','T', as.character(time[i])), 'Z', sep=''), 
                           '</time></wpt>', sep=''))
         } else {
           o <- c(o, paste('<wpt lat="', lat[i], '" lon="', lon[i], '"><time>', 
                           paste(gsub(' ','T', as.character(time[i])), 'Z', sep=''), '</time><name>', 
                           name[i], '</name></wpt>', sep=''))
         }
       }
     }
     
     o <- c(o, '</gpx>')
     
     if (is.character(out_file) || inherits(out_file, "connection")) 
       cat(o, file=out_file, sep='\n')
   }
   
  
  #run function
  write_gpx(lat = latlong2$lat, 
            latlong2$lon, 
            name = latlong2$name, 
            out_file = paste0("./output/", site_ID, "/", site_ID, "_GPX.gpx")
            )
   
  
  
  
##### EXPORT STATION TABLE FOR GBM ACCESS DATABASE ===============================
  #(append to tblStations)
  ##BRobinson Nov 6, 2023: no need to recalculate coordinates in UTMs since they were originally in UTMs
  #Xcl and Ycl witnin pts should still have correct coordinates in UTMs
  utm2 <- pts
  utm2$Xcl <- st_coordinates(pts)[,1]
  utm2$Ycl <- st_coordinates(pts)[,2]
  
  #BRObinson Nov 6, 2023: CentX and CentY are not needed in this table. However, X and Y, which are the 
  #coordinates in lat/long (NAD83), needs to be included in this export
  utm2 = subset(utm2, select= c(PCODE, SITE, STN, label, CentX, CentY, Xcl, Ycl, Datum, Zone))
  

  write.csv(utm2, 
            paste0("./output/", site_ID, "/tblStations_", site_ID, ".csv"), 
            row.names=F
            )
  
  #import this into GBM Access database (append to tblSections)

  
  
##### EXPORT TABLE FOR NATURE COUNTS DATABASE ==================================
  
  nc <- latlong %>% dplyr::select(PCODE, SITE, label, lat, lon) %>% 
    rename(project_id = PCODE, 
           site_id = SITE, 
           point_id = label, 
           latitude = lat, 
           longitude = lon) 
  
  
  write_csv(nc, 
            paste0("./output/", site_ID, "/", site_ID, "_NatureCounts_Upload.csv"))
  
  

### CREATE QGIS EXPORTS FOR MAKING MAP =========================================
  
  #read QGIS style layer
  style <- st_read("./data/GIS/layer_styles.gpkg", layer = "layer_styles")
  
  
  #### MAKE MORE EFFICIENT
  bbox <- st_buffer(site, 15000)
  bbox <- st_transform(bbox, crs = st_crs(site))
  
    road <- NULL
    
  if ("Alberta" %in% can$NAME_1) {
    AB_road <- st_read(dsn = "./data/GIS/road_network.gpkg", layer = "AB_road")
    AB_road <- st_transform(AB_road, crs = st_crs(site)) 
    AB_road <- st_crop(AB_road, st_bbox(bbox))
  

    # Merge into the main road object
    if (is.null(road)) {
      road <- AB_road
    } else {
      road <- bind_rows(road, AB_road)
    }
  }

  if ("Saskatchewan" %in% can$NAME_1) {
    SK_road <- st_read(dsn = "./data/GIS/road_network.gpkg", layer = "SK_road")
    SK_road <- st_transform(SK_road, crs = st_crs(site))
    SK_road <- st_crop(SK_road, st_bbox(bbox))
  
    # Merge into the main road object
    if (is.null(road)) {
      road <- SK_road
    } else {
      road <- bind_rows(road, SK_road)
    }
  }
  
    if ("Manitoba" %in% can$NAME_1) {
      MB_road <- st_read(dsn = "./data/GIS/road_network.gpkg", layer = "MB_road")
      MB_road <- st_transform(MB_road, crs = st_crs(site)) 
      MB_road <- st_crop(MB_road, st_bbox(bbox))
    
    
    # Merge into the main road object
    if (is.null(road)) {
      road <- MB_road
    } else {
      road <- bind_rows(road, MB_road)
    }
  }
    
    
  hwy <- NULL
    
    if ("Alberta" %in% can$NAME_1) {
      AB_hwy <- st_read(dsn = "./data/GIS/road_network.gpkg", layer = "AB_hwy")
      AB_hwy <- st_transform(AB_hwy, crs = st_crs(site)) 
      AB_hwy <- st_crop(AB_hwy, st_bbox(bbox))
    
    
    # Merge into the main road object
    if (is.null(hwy)) {
      hwy <- AB_hwy
    } else {
      hwy <- bind_rows(hwy, AB_hwy)
    }
  }
    
    
    if ("Saskatchewan" %in% can$NAME_1) {
      SK_hwy <- st_read(dsn = "./data/GIS/road_network.gpkg", layer = "SK_hwy")
      SK_hwy <- st_transform(SK_hwy, crs = st_crs(site))
      SK_hwy <- st_crop(SK_hwy, st_bbox(bbox))
    
    
    # Merge into the main hwy object
    if (is.null(hwy)) {
      hwy <- SK_hwy
    } else {
      hwy <- bind_rows(hwy, SK_hwy)
    }
  }
    
    if ("Manitoba" %in% can$NAME_1) {
      MB_hwy <- st_read(dsn = "./data/GIS/road_network.gpkg", layer = "MB_hwy")
      MB_hwy <- st_transform(MB_hwy, crs = st_crs(site)) 
      MB_hwy <- st_crop(MB_hwy, st_bbox(bbox))
    
    
    # Merge into the main hwy object
    if (is.null(hwy)) {
      hwy <- MB_hwy
    } else {
      hwy <- bind_rows(hwy, MB_hwy)
    }
  }
    gc()
  
    
  st_write(obj = site, 
           dsn = paste0("./output/", site_ID, "/", site_ID, "_QGIS.gpkg"), 
           layer = "boundary", 
           driver = "GPKG",
           append = TRUE)
  
  st_write(obj = fullgrid, 
           dsn = paste0("./output/", site_ID, "/", site_ID, "_QGIS.gpkg"),
           layer = "sec_grid", 
           driver = "GPKG",
           append = TRUE)
  
  st_write(obj = sec_sample, 
           dsn = paste0("./output/", site_ID, "/", site_ID, "_QGIS.gpkg"), 
           layer = "grid_labels", 
           driver = "GPKG", 
           append = TRUE
  )
  
  st_write(obj = road, 
           dsn = paste0("./output/", site_ID, "/", site_ID, "_QGIS.gpkg"), 
           layer = "road", 
           driver = "GPKG", 
           append = TRUE
  )
  
  st_write(obj = hwy, 
           dsn = paste0("./output/", site_ID, "/", site_ID, "_QGIS.gpkg"), 
           layer = "hwy", 
           driver = "GPKG", 
           append = TRUE
  )
  
  st_write(obj = pts, 
           dsn = paste0("./output/", site_ID, "/", site_ID, "_QGIS.gpkg"), 
           layer = "pts", 
           driver = "GPKG",
           append = TRUE
  )
  
  st_write(obj = style,
           dsn = paste0("./output/", site_ID, "/", site_ID, "_QGIS.gpkg"), 
           layer = "layer_styles", 
           driver = "GPKG",
           append = TRUE
  )
  

##### CREATE INTERACTIVE MAP ===================================================
  
  
  
  tmap_mode("view")
  
  tmap_options(check.and.fix = TRUE)
  
  map_interactive <- tm_shape(site) + tm_borders(col = "black",lwd = 3) + 
    tm_basemap("Esri.WorldImagery") +
    tm_shape(secToKeep) + tm_borders(col = "#393939") + 
    tm_shape(road) + tm_lines(col = "orange") + 
    tm_shape(hwy) + tm_lines(col = "yellow", lwd = 3) +
    tm_shape(pts) + tm_dots(col = "#e8e8e8", border.col = "#393939") + 
    tm_shape(sec_sample) + tm_borders(alpha = 0) + tm_text(text = "rank")
  
                                                         
                                                           #   col = "black",
                                                           # fontface = 2,
                                                           # shadow = TRUE)
  
  
  # return number of base sample grids and overdraw grids
  grids <- c(print(paste0("Base Sample Size: ", print(n))), print(paste0("Overdraw Sample Size: ", print(n.over))))
  
  #o output of number of grids and render map
  return(list(map_interactive, grids))
  
}






# shapefile: enter the name of the shapefile, not including the suffix .shp
# site_name_2: if the boundary includes another shapefile
# site_ID: enter the four letter site code
# pcode: project code
# sample_size: number of desired base sample grids, if entered will override calculated 20%
# overdraw_size:  number of desired overdraw grids, if entered will override calculated 1/4 of sample size

grts_GBMP(site_name = "masefield",
          site_name_2 = NULL,
          site_ID = "SMAS", 
          pcode = "GBMP",
          sample_size = 10,
          overdraw_size = 5)

