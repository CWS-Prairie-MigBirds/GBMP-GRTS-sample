
          ### CWS GRASSLAND BIRD MONITORING PROGRAM ###
  ### Generalized Random Tessellation Stratified (GRTS) sample ###

      # Barry Robinson, Janine McManus, Josiah Van Egmond #
                  # last updated: July 2023 #

### INFO =======================================================================
#    This is a self-contained project to generate a ranked random sample of...
# ...legal land sections or quarter sections in Alberta, Saskatchewan or...
#....Manitoba with point count locations 400 m apart.



### INSTRUCTIONS
# Move your site boundary shapefile into the data > boundaries folder
# Run the function line
# Input the arguments into the function at the bottom of the page:
    # site_name: the name of the shapefile in the boundaries folder, without suffix .shp
    # site_ID: the 4-letter site code, starting with the province letter
    # pcode: project code
    # sample_size: if you want to manually input number of grids, otherwise the sample size will be 20% of available grids
    # overdraw_size: to manually input the number of overdraw grids, otherwise will be 25% of sample size
# check the interactive map output to see that everything looks normal and there is not too much clumping
# rerun the function if necessary
# NOTE: delete the entire output folder before rerunning to ensure everything is properly exported



### Functionality to be added:
    # make pretty map
    # generate sample from only legal land description



### FUNCTION ===================================================================


### run this line: 
grts_GBMP <- function(site_name, 
                      site_name_2 = NULL, 
                      site_ID, 
                      pcode, 
                      sample_size = NULL, 
                      overdraw_size = NULL)  { 
   
  
##### LOAD PACKAGES ============================================================
  
  # custom function to check for and automatically load packages
  
  using<-function(...) {
    libs<-unlist(list(...))
    req<-unlist(lapply(libs,require,character.only=TRUE))
    need<-libs[req==FALSE]
    if(length(need)>0){ 
      install.packages(need)
      lapply(need,require,character.only=TRUE)
    }
  }
  
  # required packages will be automatically installed and loaded
  
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
  dir.create(paste0("./output/", site_ID))
    
  
  sf_use_s2(FALSE)
  
  #load in site shapefile
  site <- st_read(dsn = "./data/boundaries", layer = site_name)
  # site <- st_read("./data/boundaries/test.shp")

    ##some sites have boundary files that have multiple polygons that need to be grouped together
  site$name <- "site_name"
  site <- site %>% group_by(name) %>% dplyr::summarise()
  
  
  
##### SELECT CRS FOR SITE FROM LOCATION ========================================
  
  # check for missing crs, set to WGS84 to get coordinates for utm_zone function - JVE
  
  if (is.na(st_crs(site))) {
    st_crs(site) <- 4326
  }
  

 #find coords for center of site
  centroid <- st_centroid(site) %>% 
    st_transform(centroid, crs = 4326) %>%  
    st_coordinates(centroid)
  
  #find utm zone
  zone <- utm_zone(centroid[1], centroid[2])            #calculate UTM zone using centroid coordinates
  utm_zone_num <- as.numeric(sub("[C-Z]$", "", zone))   #drop letter from zone, since sites are only in northern hemisphere
  epsg <- 26900 + utm_zone_num
  
  #transform site to correct crs. crs of other objects will be based off site crs
  site <- st_transform(site, crs = epsg)
  
  #if there is a second shapefile
  
  if(is.null(site_name_2) == FALSE) {
    site_2 <- st_read(dsn = "./data/boundaries", layer = site_name_2)
    site_2 <- st_transform(site_2, crs = st_crs(site))
    site <- bind_rows(site, site_2)
    
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
  
  can <- gadm(country="CAN", level=1, path="./Data")
  can <- st_as_sf(can)  
  can <- filter(can, NAME_1 %in% provinces)

  #transform Canada province data to site crs and intersect
  can <- st_transform(can, crs = st_crs(site))
  can <- st_intersection(can, site)

  # read in province grid data and change legal land description column name to LLD

  # empty grid object so province grids can be bound to it if site crosses two provinces
  
  grid <- NULL
  
  if ("Alberta" %in% can$NAME_1) {
    grid_A <- st_read(dsn = "./data/grid/section_grid_BCR11.gpkg", layer = "AB") %>%
      rename(LLD = DESCRIPTOR)
    grid_A <- st_transform(grid_A, crs = st_crs(site))
    
    # Merge into the main grid object
    if (is.null(grid)) {
      grid <- grid_A
    } else {
      grid <- bind_rows(grid, grid_A)
    }
  }
  
  if ("Saskatchewan" %in% can$NAME_1) {
    grid_S <- st_read(dsn = "./data/grid/section_grid_BCR11.gpkg", layer = "SK")
    grid_S <- st_transform(grid_S, crs = st_crs(site))
    
    # Merge into the main grid object
    if (is.null(grid)) {
      grid <- grid_S
    } else {
      grid <- bind_rows(grid, grid_S)
    }
  }
  
  if ("Manitoba" %in% can$NAME_1) {
    grid_M <- st_read(dsn = "./data/grid/section_grid_BCR11.gpkg", layer = "MB") %>%
      rename(LLD = LEG_DESC)
    grid_M <- st_transform(grid_M, crs = st_crs(site))
    
    # Merge into the main grid object
    if (is.null(grid)) {
      grid <- grid_M
    } else {
      grid <- bind_rows(grid, grid_M)
    }
  }
  
  

##### PROCESS PROVINCE SECTION GRID ============================================
  
  #Query out all land sections that intersect with each site polygon
  SecGrid <- st_intersection(x = site, y = grid)
  SecGrid$area = st_area(SecGrid)
  
  
  #Need to remove slivers of sections! 
  #But then we also need to filter the full grid list by Legal description
  # That way we have full sections for creating the centroids and grids.
  
  # 1 sq mi = 259000 m2, divide by 2 for half grid...
  #...using 2.05 instead of 2 to have some tolerance for areas that are slightly... 
  #...less than half a section due to imperfect polygons and rounding error
  
  secToKeep = SecGrid %>% dplyr::filter(as.numeric(area) >= 2590000/2.05)
 
  
  # large or small site? 
  
     if (nrow(secToKeep) <= 5) {
      site_type = "small" 
      } 
    else {
      site_type = "large"
      }
  
  

  
##### SMALL SITES GRID ===========================================================  
 
  #   check if it's a small site, if so load in appropriate province grids,...
  #...transform them and bind them together
  
  if (site_type == "small") {
   
   suppressWarnings(rm(grid, grid_A, grid_S, grid_M)) #remove grid to save RAM
   gc() #garbage collection to free up RAM

   grid <- NULL

   if ("Alberta" %in% can$NAME_1) {
     grid_A <- st_read(dsn = "./data/grid/quarter_grid_BCR11.gpkg", layer = "AB") %>%
       rename(LLD = DESCRIPTIO)
     grid_A <- st_transform(grid_A, crs = st_crs(site))

     # Merge into the main grid object
     if (is.null(grid)) {
       grid <- grid_A
     } else {
       grid <- bind_rows(grid, grid_A)
     }
   }

   if ("Saskatchewan" %in% can$NAME_1) {
     grid_S <- st_read(dsn = "./data/grid/quarter_grid_BCR11.gpkg", layer = "SK")
     grid_S <- grid_S %>% unite(col = "LLD", c("QSECT", "PSECT", "PTWP", "PRGE", "PMER"), sep = "-")
     grid_S <- st_transform(grid_S, crs = st_crs(site))

     # Merge into the main grid object
     if (is.null(grid)) {
       grid <- grid_S
     } else {
       grid <- bind_rows(grid, grid_S)
     }
   }


   if ("Manitoba" %in% can$NAME_1) {
     grid_M <- st_read(dsn = "./data/grid/quarter_grid_BCR11.gpkg", layer = "MB") %>%
       rename(LLD = LEGAL_DESC) %>% 
       rename(zone_1 = ZONE)
     grid_M <- st_transform(grid_M, crs = st_crs(site))

     # Merge into the main grid object
     if (is.null(grid)) {
       grid <- grid_M
     } else {
       grid <- bind_rows(grid, grid_M)
     }
   }

   


    # intersect quarter grid with site boundary
    grid <- st_intersection(x = site, y = grid)
    
    # calculate areas to filter out partial grids
    grid$area = st_area(grid)
    
    # need to remove slivers of sections!
    secToKeep = grid %>% dplyr::filter(as.numeric(area) >= 647497/2.05)
    
    
  }
  

    
  ### filter by legal land description
  fullgrid <- grid %>% dplyr::filter(LLD %in% c(secToKeep$LLD))
  
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
  
  
  suppressWarnings(rm(grid, grid_A, grid_S, grid_M)) #remove grid to save RAM
  gc() #garbage collection to free up RAM
  
  
##### DRAW RANDOM SAMPLE =======================================================
  
  #calculate sample size, either 20% of available grid cells or 5 cells, whichever is greater
  #or custom input
  #overdraw size is 1/4 of sample size, may be too large for big sites
  #can manually enter overdraw sample size
  
  #BROBINSON: These ifelse statements need to be reviewed relative to the sample size rules we've established. As they are, they are more complex that the simple rule written above
  if(is.null(sample_size)) {
    if(site_type == "small") {
      n <- nrow(secToKeep)
      if(0.2*n >= 10) {
        n <- round(0.2*n)
      } else { #if(0.2*n < 5 & n > 10){ #BRobinson: This isn't needed anymore because site size is now defined above
        n <- 10
      }
    } else {
      n <- nrow(secToKeep)
      if(0.2*n >= 5) {
        n <- round(0.2*n)
      } else { #if(0.2*n < 5 & n >5){ #BRobinson: This isn't needed anymore because site size is now defined above
        n <- 5
      }
    }
  } else {
    n <- sample_size
  }
  
  
  
  if(is.null(overdraw_size)) {
    n.over <- ceiling(0.25*n)
  } else {
    n.over <- overdraw_size
  }
  
  
  ##now we can run the GRTS sample, 
  
  sample = grts(sframe = centroid, n_base = n, seltype = "equal", n_over = n.over)
  
  ##get samples out of list
  
  sample_b <- sample$sites_base
  sample_o <- sample$sites_over
  
  sample2 <- rbind(sample_b, sample_o)
  
  #Check to see if the draw is funky before proceeding, if heavily clumped, run grts() again.

  #Create a rank column
  
  sample2$rank <- 1:nrow(sample2)
  

  # create section sample grids
  
  sec_sample <- fullgrid %>% dplyr::filter(LLD %in% c(sample2$LLD)) 
  
  #add rank column to sec_sample
  
  rank_col <- sample2[c("LLD", "rank")] %>% 
    st_drop_geometry()
  
  sec_sample <- left_join(as.data.frame(sec_sample), rank_col, by = "LLD")
  sec_sample <- st_as_sf(sec_sample)
  
  # #write section sample shapefile
  # st_write(obj = sec_sample, 
  #          dsn = paste0("./output/", site_ID, "/", site_ID, ".gpkg"), 
  #          layer = paste0(site_ID, "_SecSample"), 
  #          driver = "GPKG", 
  #          append = TRUE
  #          )
  
  ###Create csv table of the section sample to be imported into Access (Append to tblSections)
  
  stable <- st_drop_geometry(sample2)
  
  # UPDATE site code
  stable$SITE <- site_ID #change to correct code for site
  

  #export table for GBM Access database
  write.csv(subset(stable, select = -c(replsite, caty)), 
            file=paste0("./output/", site_ID, "/tblSections_", site_ID, ".csv"), 
            row.names=F,
            ) 
   
  #import this into GBM Access database (append to tblSections)
  
  

  
  
  
##### CALCULATE POINT COUNT LOCATIONS ==========================================
  
  # calculate XY coordinates for point count locations...
  #...within the each section drawn for the sample
  
  
  
  if(site_type =="small") { ##### small site point count locations, 4 per grid
    
    grid.exp <- stable[rep(seq_len(nrow(stable)), 4),]
    #order by rank
    grid.exp <- grid.exp[order(grid.exp$rank),]
    #add field to indicate point number within each section (1-4)
    grid.exp$pts <- as.vector(rep(1:4, nrow(stable)))
    #create fields for X and Y coordinates of points
    grid.exp$Xcl <- 0
    grid.exp$Ycl <- 0
    
    #Calculate coordinates for each point based on it's distance from the centroid
    
    grid.exp$Xcl <- ifelse(grid.exp$pts==2 | grid.exp$pts==4,
                           grid.exp$CentX +200, grid.exp$CentX-200)
    grid.exp$Ycl <- ifelse(grid.exp$pts==2 | grid.exp$pts==1,
                           grid.exp$CentY +200, grid.exp$CentY-200)
    
 
    
  } else { ##### large site point count locations, 16 per grid
  
  
  grid.exp <- stable[rep(seq_len(nrow(stable)), 16),]
  #order by rank
  grid.exp <- grid.exp[order(grid.exp$rank),]
  #add field to indicate point number within each section (1-16)
  grid.exp$pts <- as.vector(rep(1:16, nrow(stable)))
  #create fields for X and Y coordinates of points
  grid.exp$Xcl <- 0
  grid.exp$Ycl <- 0
  
  #Calculate coordinates for each point based on it's distance from the centroid
  grid.exp$Xcl <- ifelse(grid.exp$pts==1 | grid.exp$pts==5 | grid.exp$pts==9 | grid.exp$pts==13,
                         grid.exp$CentX+600, grid.exp$Xcl)
  grid.exp$Xcl <- ifelse(grid.exp$pts==2 | grid.exp$pts==6 | grid.exp$pts==10 | grid.exp$pts==14,
                         grid.exp$CentX+200, grid.exp$Xcl)
  grid.exp$Xcl <- ifelse(grid.exp$pts==3 | grid.exp$pts==7 | grid.exp$pts==11 | grid.exp$pts==15,
                         grid.exp$CentX-200, grid.exp$Xcl)
  grid.exp$Xcl <- ifelse(grid.exp$pts==4 | grid.exp$pts==8 | grid.exp$pts==12 | grid.exp$pts==16,
                         grid.exp$CentX-600, grid.exp$Xcl)
  grid.exp$Ycl <- ifelse(grid.exp$pts==1 | grid.exp$pts==2 | grid.exp$pts==3 | grid.exp$pts==4,
                         grid.exp$CentY+600, grid.exp$Ycl)
  grid.exp$Ycl <- ifelse(grid.exp$pts==5 | grid.exp$pts==6 | grid.exp$pts==7 | grid.exp$pts==8,
                         grid.exp$CentY+200, grid.exp$Ycl)
  grid.exp$Ycl <- ifelse(grid.exp$pts==9 | grid.exp$pts==10 | grid.exp$pts==11 | grid.exp$pts==12,
                         grid.exp$CentY-200, grid.exp$Ycl)
  grid.exp$Ycl <- ifelse(grid.exp$pts==13 | grid.exp$pts==14 | grid.exp$pts==15 | grid.exp$pts==16,
                         grid.exp$CentY-600, grid.exp$Ycl)
  
  }
  
  
  
  #Modify table so that it can be brought into GBM Access database and ArcMap. 
  # ***MODIFY SITE, PCODE, Datum, and Zone***
  
  grid.exp$SITE <- site_ID 
  grid.exp$PCODE <- pcode
  grid.exp$Datum <- "NAD83"   
  grid.exp$Zone <- zone
  grid.exp$STN <- paste(sprintf("%02d", grid.exp$rank),sprintf("%02d",grid.exp$pts),sep="-")
  grid.exp$label <- paste(grid.exp$SITE,grid.exp$STN,sep="-")
  
  #convert to sf object
  pts <- st_as_sf(grid.exp, coords = c("Xcl", "Ycl"))
  
  #update point crs to match site
  st_crs(pts) = crs(site)
  
  #clip out any points that fall outside of original study site
  pts <- st_intersection(pts, site)
  
  pts <- pts %>% dplyr::select(label, everything())
 
###BRobinson_Nov 6, 2023: this is exporting the point shapefile in the projection of the section grid,
  #which could be different from province to province. Save exports until the end after all 
  #shapefiles have been converted back to WGS83
  
  
  #  #write the point sample shapefile
  # 
  # st_write(obj = pts, 
  #          dsn = paste0("./output/", site_ID, "/", site_ID, ".gpkg"), 
  #          layer = paste0(site_ID, "_PointSample"), 
  #          driver = "GPKG",
  #          append = TRUE
  #          )
  
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






# site_name: enter the name of the shapefile, not including the suffix .shp
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

