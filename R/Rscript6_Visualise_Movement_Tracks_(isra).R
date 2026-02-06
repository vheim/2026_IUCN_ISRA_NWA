# TOP OF CODE ----
### ====================================================================================================
### Project:    Large bodied hammerhead complex in the western North Atlantic
### Analysis:   Animate movement tracks of tagged animals
### Script:     ~~2026_IUCN_ISRA_NWA/R/Rscript6_Visualise_Movement_Tracks.R
### Author:     Vital Heim
### Version:    1.0
### ====================================================================================================

### ....................................................................................................
### Content: this script contains the code to visualise time series location data
### such as animal tracking data.
### 
### The code uses a ggplot2 visualisation approach, but base R options will be implemented
### in the future
###
### Data is visualised using a combination of tracking data, shapefiles, and rasters.
###            
### Please see TODO list at end of script for open issues.
### ....................................................................................................

### ....................................................................................................
### [A] Ready environment, load packages ----
### ....................................................................................................

# A1: clear memory ----

rm(list = ls())

# A2: isntall and load necessary packages ----

## if first time
# install.packages("tidyverse")
# install.packages("magrittr")
# install.packages("rnaturalearth") # if first time
# install.packages("rnaturalearthdata") # if first time
# install.packages("ncdf4")
# install.packages("stars")
# install.packages("raster")
# install.packages("terra")
# install.packages("tidyterra")
# install.packages("scales")
# install.packages("ggnewscale")
# install.packages("sf")
# remotes::install_github("SimonDedman/gbm.auto")
# install.packages("gbm.auto")
# install.packages("ggnewscale")
# install.packages("beepr")

## load packages
library(tidyverse)
library(magrittr)
library(rnaturalearth)
library(rnaturalearthdata)
library(ncdf4)
library(stars)
library(raster)
library(terra)
library(tidyterra)
library(scales)
library(ggnewscale)
library(sf)
# library(gbm.auto)
library(beepr)

# A3: Specify needed functions ----

## FUNCTION 1: to create custom colour ranges for plotting multiple individuals----
## Please note: if you have more than 30 individuals in your dataframe, this function will reuse
## colours
my_color_range <- function(colrange_start, colrange_end, n_colors, 
                           colrange_mid1 = NULL, colrange_mid2 = NULL,
                           alpha = 1, reverse = FALSE, max_colors = 200) {
  
  # Create the palette based on how many colors are specified
  if (is.null(colrange_mid1) && is.null(colrange_mid2)) {
    # Two-color gradient
    palette_func <- colorRampPalette(c(colrange_start, colrange_end), alpha = alpha)
  } else if (!is.null(colrange_mid1) && is.null(colrange_mid2)) {
    # Three-color gradient  
    palette_func <- colorRampPalette(c(colrange_start, colrange_mid1, colrange_end), alpha = alpha)
  } else if (!is.null(colrange_mid1) && !is.null(colrange_mid2)) {
    # Four-color gradient
    palette_func <- colorRampPalette(c(colrange_start, colrange_mid1, colrange_mid2, colrange_end), alpha = alpha)
  } else {
    # Invalid combination (mid2 without mid1)
    stop("If using colrange_mid2, you must also specify colrange_mid1")
  }
  
  # Generate colors up to the maximum or requested amount, whichever is smaller
  colors_to_generate <- min(n_colors, max_colors)
  base_colors <- palette_func(colors_to_generate)
  
  # If we need more colors than the maximum, cycle through the palette
  if (n_colors > max_colors) {
    # Calculate how many full cycles we need and remainder
    full_cycles <- n_colors %/% max_colors
    remainder <- n_colors %% max_colors
    
    # Create the full color vector by repeating and adding remainder
    if (remainder > 0) {
      colors <- c(rep(base_colors, full_cycles), base_colors[1:remainder])
    } else {
      colors <- rep(base_colors, full_cycles)
    }
  } else {
    colors <- base_colors
  }
  
  # Reverse if requested
  if (reverse) {
    colors <- rev(colors)
  }
  
  return(colors)
}

## FUNCTION 2: select fitted or predicted track files ----
select_and_read_file <- function(folder_path, pattern = "\\.rds$") {
  # Get list of files in the folder
  file_list <- list.files(
    path = folder_path,
    pattern = pattern,
    full.names = FALSE
  )
  
  # Check if any files were found
  if (length(file_list) == 0) {
    stop("No files matching the pattern were found in the specified folder.")
  }
  
  # Display files with numbers in yellow
  cat("\033[33mAvailable files:\033[0m\n")
  for (i in seq_along(file_list)) {
    cat("\033[33m", i, ": ", file_list[i], "\033[0m\n", sep = "")
  }
  
  # Prompt user to choose
  file_choice <- as.integer(readline(prompt = "Enter the number of the file you want to import: "))
  
  # Validate choice
  if (is.na(file_choice) || file_choice < 1 || file_choice > length(file_list)) {
    stop("Invalid selection. Please run the function again and choose a valid number.")
  }
  
  # Get selected filename
  selected_file <- file_list[file_choice]
  selected_path <- file.path(folder_path, selected_file)
  
  # Determine tracktype based on filename
  if (grepl("fitted", selected_file, ignore.case = TRUE)) {
    tracktype <<- "fitted"
  } else if (grepl("segmented", selected_file, ignore.case = TRUE)) {
    tracktype <<- "predicted"
  } else {
    tracktype <<- NA
  }
  
  cat("\033[33mReading:", selected_path, "\033[0m\n")
  cat("\033[33mtracktype set to:", tracktype, "\033[0m\n")
  
  # Read and return the selected file
  return(readRDS(selected_path))
}

## FUNCTION3: choose the species you want to visualise if multiple ----
select_species <- function(df, species_col = "species") {
  species_choices <- unique(df[[species_col]])
  
  choice <- menu(
    species_choices,
    title = "Select a species by number:"
  )
  
  if (choice == 0) {
    stop("No species selected.")
  }
  
  sp.f <- species_choices[choice]
  return(sp.f)
}


# A4: Specify data and saveloc ----

YOUR_IP <- "NA" # add your IP address, server name or similar if you connect via a shared drive

## Project folder
projloc <- file.path("/",YOUR_IP, "Science","Projects_current", "2026_IUCN_ISRA_NWA")

## Input data
dataloc <- file.path("/",YOUR_IP, "Science","Projects_current", "2026_IUCN_ISRA_NWA", "Data_input")
bathyloc <- file.path("/",YOUR_IP, "Science", "Data_raw", "Bathymetry_maps_GMRT", "GMRTv4_3_0_20250120topo_wider_NWA.grd")
misc_shapefileloc <- file.path("/",YOUR_IP, "Science", "Data_raw", "Shapefiles")

## Output data
saveloc <- file.path("/",YOUR_IP, "Science","Projects_current","2026_IUCN_ISRA_NWA","Data_output", "Maps_and_tracks")
### check if saveloc already present
if (!dir.exists(saveloc)) {
  dir.create(saveloc, recursive = TRUE)
  cat("New saveloc was created at:", saveloc, "\n")
} else {
  cat("Saveloc is already present","\n")
}

# A5: Define universal variables (e.g. for plotting) ----

# Plotting variables

## bathymetry colors
shallow <- "#D3E5E8"
deep <- "#2B628B"
depth <- c(deep, shallow)

## speciescolours
slewcol <- "#EDA904"
smokcol <- "#70AB27"
szygcol <- "#E26306"

## if you need individual track colours and want to use a range of colours
# colrange_start <- "#FFFFCC"
# colrange_mid1 <- "#FEB24C"
# colrange_mid2 <- "#E31A1C"
# colrange_end <- "#800026"

colrange_start <- "#FFFFCC"
colrange_mid1 <- "#8FD900"
colrange_mid2 <- "#0F9D00"
colrange_end <- "#006B29"

# A6: Define universal options, variables, etc. (e.g. for plotting) ----

options(timeout = 3000) # manually increase time out threshold (needed when downloading basemap)
options(scipen=999) # so that R doesn't act up for pit numbers  
options(warn=1) #set this to two if you have warnings that need to be addressed as errors via traceback()
options(error = function() beep(9))  # give warning noise if it fails

### ....................................................................................................
### [B] Data import ----
### ....................................................................................................

# B1: Movement data ----

## movement data
mov_tracks <- select_and_read_file(
  folder_path = file.path(dataloc, "Maps_and_tracks")
) |>
  dplyr::rename(
    shark = id
  ) %>%
  dplyr::mutate(
    shark = as.character(shark)
  )

### check if there are underscores in shark id column, i.e. there should be if you are using segmented tracks
if (exists("tracktype") && !is.na(tracktype) && tracktype == "predicted") {
  # Check if any underscores are present in the shark column
  if (any(grepl("_", mov_tracks$shark))) {
    # Rename the original shark column to shark_raw
    mov_tracks$shark_raw <- mov_tracks$shark
    
    # Create new shark column without underscore and following numbers
    mov_tracks$shark <- sub("_.*$", "", mov_tracks$shark)
    
    cat("\033[33mtracktype is 'predicted': Underscores removed. Original values stored in 'shark_raw' column.\033[0m\n")
  } else {
    cat("\033[33mtracktype is 'predicted': No underscores detected. Dataframe unchanged.\033[0m\n")
  }
} else {
  cat("\033[33mtracktype is 'fitted' or undefined: Dataframe unchanged.\033[0m\n")
}

## metadata for length etc.
meta <- read_csv(file.path(dataloc, "Telemetry", "Data_Sphyrna_SPOT_tags_metadata.csv"))|>
  dplyr::filter( # keep species of movement data
    species %in% c("S.mokarran", "S.lewini","S.zygaena"),
    group %in% c("Florida Keys", "Jupiter", "Jupiter, FL", "Marquesas", "South Carolina", "Tampa")
  ) |>
  dplyr::select( # keep needed columns only
    ptt_id,
    datetime_deployment_local,
    stl,
    sex,
    species
  ) |>
  dplyr::mutate(
    shark = as.character(ptt_id)
  )

## combine movement data with metadata
mov_tracks %<>% left_join(meta)  # , by = join_by(shark == id) # doesn't work naming columns, has gotten worse.

## VERY IMPORTANT: The movement data needs to have time stamps in ascending order,
## double check that this is true by ordering the dataframe
mov_tracks %<>%
  dplyr::arrange(
    shark,
    date
  )

# B2: Basemap data and shapefiles ----

# *B2.1: Landmasses basemap ----
## Basemap data for landmasses- HQ
## Option 1 - required downloaded shapefile
### define resolution: 1=c, 2=l,3=i,4=h,5=f, 1:5 increasing quality
# res <- "f"
### read in worldmap
# world <- sf::st_read(dsn = paste0(shapefileloc,"/",res,"/GSHHS_", res, "_L1.shp"), layer = paste0("GSHHS_", res, "_L1"), quiet = TRUE) # read in worldmap

## World shapefile - LQ
world <- rnaturalearth::ne_countries(scale = 10, returnclass = "sf")

# *B2.2: Bathymetry maps and other rasters ----

bathyR <- raster::raster(bathyloc)
## todo: use rast() rather than raster(), rast() is from the terra package, and terra has more options and will replace raster

#SpatExtent : -84.3068857607634, -74.5949697079866, 22.996942201359, 35.7020417853476 (xmin, xmax, ymin, ymax)
proj4string(bathyR) # from raster package
#> proj4string(bathy)
#[1] "+proj=longlat +datum=WGS84 +no_defs"

## find min and max values
raster::setMinMax(bathyR)
# minValue(bathyR); maxValue(bathyR)
### we are working with ocean depth data, i.e. anyhting above 0 meter elevation we dont really need and can set to NA
bathyR <- clamp(bathyR, upper = 10, useValues = F) # this is really just an aesthetic thing, needs to be adjusted if you print a legend that goes between 0 to -8000 meters

## Extent of the raster
# bbox(bathyR)
# min       max
# s1 -99.10547 -65.51367
# s2  22.56938  42.90172

## Dimension
# dim(bathyR)
#> dim(bathy)
#[1] 2777 3822    1

## Resolution
# res(bathyR)
#[1] [1] 0.008789062 0.007321693

## Visualise
plot(bathyR)

## make a df needed for later plotting
raster.df <- as.data.frame(bathyR, xy = T)

## prep for TERRA PACKAGE.....
bathyterra <- terra::rast(bathyR)
plot(bathyterra)

# *B2.2: Other shapefiles, e.g. EEZ or closure boundaries ----

### Bahamian EEZ
# bah_eez <- read_sf("//Sharktank/Science/Data_raw/Shapefiles/Bahamas/Bahamas_EEZ_boundary.shp")
# sf::st_crs(bah_eez)
# bah_eez <- sf::st_transform(bah_eez, st_crs = proj4string(bathyR)) # bring to same CRS as bathymetry raster

### US State boundaries
usstates <- sf::st_read(dsn = file.path(misc_shapefileloc, "USA", "US_State_Boundaries", "US_State_Boundaries.shp"), quiet = TRUE) # read in worldmap

### Federal vs. state waters USA
fed_state_boundary <- sf::st_read(dsn = file.path(misc_shapefileloc, "USA", "US_Federal_State_Waters", "Federal_Waters.shp"), quiet = TRUE) # read in worldmap

### US EEZ
eez_usa <- sf::st_read(dsn = file.path(misc_shapefileloc, "USA", "US_EEZ","eez.shp"), quiet = TRUE) # read in worldmap


### ....................................................................................................
### [C] Data housekeeping and preparation ----
### ....................................................................................................

# C1: housekeeping detection dataframe ----

## filter the needed columns
mov_tracks_f <- mov_tracks %>%
  # dplyr::mutate( # define any additional variables that might be handy for plotting
  #   month = format(dateEST, format = "%b", tz = "US/Eastern"), # get rid of year for later season variable definitions
  #   season_bahamas = with(.,case_when(# summer/winter based on temperature by van Zinnicq Bergmann et al. 2022 ; summer =  1st Jun to 30th Nov, winter = 1st Dec to 31st May
  #     month %in% c("May", "Jun", "Jul", "Aug", "Sep", "Oct") ~ "wet",
  #     month %in% c("Nov", "Dec", "Jan", "Feb", "Mar", "Apr") ~ "dry",
  #     TRUE ~ "seasonnonexistent")
  #   )
  # ) %>%
  dplyr::select(
    shark,
    date,
    lon,
    lat,
    species
  )

# C2: filter your dataframe by species if needed ----

## Define species and filter data
sp.f <- select_species(mov_tracks_f)

## filter
df_plot <- mov_tracks_f %>%
  dplyr::filter(
    if (sp.f == "NA") TRUE else species %in% sp.f
  )

# C3: prepare parameters for static map ----

## plot extent
min(df_plot$lon);max(df_plot$lon)
min(df_plot$lat);max(df_plot$lat)

xmin <- ceiling(min(df_plot$lon)-1);xmin
xmax <- ceiling(max(df_plot$lon)+1);xmax
ymin <- floor(min(df_plot$lat)-1);ymin
ymax <- ceiling(max(df_plot$lat)+1);ymax

xlabs = seq(xmin, xmax, 5)
ylabs = seq(ymin, ymax, 5)

### ....................................................................................................
### [D] Plot your movement tracks ----
### ....................................................................................................

## !! NEEDS UDPATEfor BASE R SOLUTION!!##

# Pre-step: create clean species label
sp_clean <- gsub("\\.", "", sp.f)

# D1: static map - group/species ----

trackmap <- ggplot() +
  
  # bathymetry raster
  ## stars package
  # stars::geom_stars(data = stars::st_downsample(bathystar, 50) |> sf::st_transform(proj4string(bathyR)), inherit.aes = FALSE) + # choose your epsg code accordinlgy (here EPSG:3857 is for WGS 84 / Pseudo-Mercator -- Spherical Mercator, Google Maps, OpenStreetMap, Bing, ArcGIS, ESRI) - personally don't like it takes ages even with downsample
  ## raster package
  # ggplot2::geom_raster(data = raster.df , aes(x = x, y = y, fill = layer)) +
  ## terra & tidyterra package - by far performs best
  tidyterra::geom_spatraster(data = bathyterra) +
  ## add a legend/gradient to your bathymetry raster
  ggplot2::scale_fill_gradientn(colors = depth, guide = "none") +   # guide = "none" prevents legend
  # # labs(x=NULL, y=NULL,
  #      fill = 'Depth [m]',
  #      color = 'Depth [m]')+
  # theme(panel.grid = element_blank(), legend.title = element_text(face = "bold")) +
  
  #coord_quickmap() +
  
  # start a new scale
  new_scale_fill() +
  
  # add shapefiles that need to go behind tracks (e.g. eez boundaries)
  ## us eez
  ggplot2::geom_sf(data = eez_usa, colour = "black", fill = NA, size = .25) +
  ## us state vs federal waters
  # ggplot2::geom_sf(data = fed_state_boundary, colour = "black", fill = NA, size = .25) +
  
  # lines and points
  ## add track positions
  geom_point(data = df_plot,
             aes(x=lon,y=lat, group = shark),
             # group = seq_along(Index), # somehow needed if you want to keep data points, i.e. if you want keep_last = T in gganimate::animate()
             fill="black",
             shape = 21,
             alpha = 1, size = 1) + #original point size 1.15, new is .5
  ## add movement paths
  geom_path(data = df_plot,
            aes(x=lon,y=lat, group = shark, color = shark),
            alpha = 1, linewidth = .4) + #original linewidth ist .5,
  ## colour individual tracks
  scale_color_manual(name = "Shark-ID", # sets legend name
                     values = my_color_range(colrange_start, colrange_end,
                                             n_colors = length(sort(unique(df_plot$shark))),
                                             colrange_mid1, colrange_mid2),
                     # guide = "none",
                     drop = F) +
  
  # add a symbol rather than a simple data point
  # ggimage::geom_image(aes(image = image), size = 0.1)+ # NOT READY YET
  
  # basemap
  ggplot2::geom_sf(data = world, fill = "gray80", color = "black", size = .25, inherit.aes = F) +
  
  # raster
  #geom_polygon(data = r, aes(x = long, y = lat)) +
  
  # additional shapefiles that go on top of terrestrial or marine data
  ggplot2::geom_sf(data = usstates, colour = "black", fill = NA, size = .25) +
  
  # define plot limits
  ggplot2::coord_sf(xlim = c(xmin, xmax),
                    ylim = c(ymin+1, ymax),
                    expand = T) +
  
  # plot axis labels
  scale_x_continuous(labels = function(x) paste0(x, '\u00B0', "W")) +
  scale_y_continuous(labels = function(x) paste0(x, '\u00B0', "N")) +
  
  # formatting
  # scale_fill_viridis_c(option = "inferno")+
  # scale_color_viridis_c(option = "inferno")+
  
  ## continous variables
  # scale_fill_gradientn(colors = seasoncol) +
  # scale_color_gradientn(colors = seasoncol) +
  # scale_size_continuous(range = c(0.1,10))+
  
  ## categorical variables
  ### shape and legend for seasons
  # scale_shape_manual(values = seasonsym,
  #                    name= "Season", # sets legend name
  #                    drop = F) +
  ### fill colour for shark symbols
  # scale_fill_manual(name = "Shark-ID",
  #                   values = all_cols,
  #                   drop = F,
  #                   # guide = "none" # removes legend
  #                   guide = guide_legend(override.aes = list(fill=all_cols, shape = 22, color = "black")) # use this if you use a wake effect in the animation or if you want a dot legend for individuals
  # ) +
  
  # Add scale bar
  # ggspatial::annotation_scale(location = "bl",  # bottom-left
  #                             width_hint = 0.3,   # width as fraction of plot
  #                             style = "bar",      # or "ticks"
  #                             unit_category = "metric") +  # uses km automatically
  # 
  # # Add north arrow
  # ggspatial::annotation_north_arrow(location = "tr",  # top-right
  #                                   which_north = "true", 
  #                                   # style = "north_arrow_minimal",
  #                                   height = unit(1, "cm"),
  #                                   width = unit(1, "cm")) +

  # Define your theme aesthethics
  theme(panel.grid = element_blank(), #legend.title = element_text(face = "bold"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = .75), # Add a black border around the plot
        # panel.background = element_rect(fill = "#C1D4E1"), # as alternative as GMRT background takes ages to animate
        # panel.background = element_rect(fill = "white"), # as alternative as GMRT background takes ages to animate
        text = element_text(family = "serif", face = "plain"), # all text to Times New Roman look-a-like
        plot.background = element_rect(fill = "white", color = "white"),    # Set overall plot background to white
        axis.text.x = element_text(size = 15), # change the font size of x.axis text
        axis.text.y = element_text(size = 15), # change the font size of y.axis text
        axis.title.x = element_blank(), # removes x axis title (i.e. "lon")
        axis.title.y = element_blank(), # removes y axis title (i.e. "lat")
        legend.title = element_text(size = 12, face = "bold"), # change the font size of the legend titles
        legend.text = element_text(size = 10),# change the font size of legend text
        legend.background = element_blank(),
        legend.key = element_blank(),
        plot.title = element_text(face = "bold", size = 14)
  ); trackmap

## save your map
ggsave(file.path(saveloc, paste0("Movement_tracks_CTCRW_", tracktype,"_",sp_clean,".tiff")), width = 21, height = 21, units = "cm", dpi = 300)

# D2: static map - individual ----

## Get all unique shark IDs
ind <- unique(df_plot$shark)

## Create subdirectory for individual tracks
individual_maps_dir <- file.path(saveloc, paste0("Individual_tracks_", sp_clean, "_",tracktype))
dir.create(individual_maps_dir, showWarnings = FALSE, recursive = TRUE)

## Loop through each animal and save individually
for (shark_id in ind) {
  
  # Filter data for current shark
  df_current <- df_plot[df_plot$shark == shark_id, ]
  
  # Create plot for current shark
  ind_trackmap <- ggplot() +
    
    # bathymetry raster
    tidyterra::geom_spatraster(data = bathyterra) +
    ggplot2::scale_fill_gradientn(colors = depth, guide = "none") +
    
    # start a new scale
    new_scale_fill() +
    
    # add shapefiles that need to go behind tracks (e.g. eez boundaries)
    ggplot2::geom_sf(data = eez_usa, colour = "black", fill = NA, size = .25) +
    
    # lines and points
    ## add track positions
    geom_point(data = df_current,
               aes(x=lon, y=lat),
               fill="black",
               shape = 21,
               alpha = 1, size = 1.15) +
    ## add movement paths
    geom_path(data = df_current,
              aes(x=lon, y=lat),
              alpha = 1, linewidth = .5, color = "#800026") +
    ## colour individual tracks
    # scale_color_manual(name = "Shark-ID",
    #                    values = my_color_range(colrange_start, colrange_end,
    #                                            n_colors = length(sort(unique(df_plot$shark))),
    #                                            colrange_mid1, colrange_mid2),
    #                    guide = "none",
    #                    drop = F) +
    
    # basemap
    ggplot2::geom_sf(data = world, fill = "gray80", color = "black", size = .25, inherit.aes = F) +
    
    # additional shapefiles
    ggplot2::geom_sf(data = usstates, colour = "black", fill = NA, size = .25) +
    
    # define plot limits
    ggplot2::coord_sf(xlim = c(xmin, xmax),
                      ylim = c(ymin+1, ymax),
                      expand = T) +
    
    # plot axis labels
    scale_x_continuous(labels = function(x) paste0(x, '\u00B0', "W")) +
    scale_y_continuous(labels = function(x) paste0(x, '\u00B0', "N")) +
    
    # Define theme aesthetics
    theme(panel.grid = element_blank(),
          panel.border = element_rect(color = "black", fill = NA, linewidth = .75),
          text = element_text(family = "serif", face = "plain"),
          plot.background = element_rect(fill = "white", color = "white"),
          axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.title = element_text(size = 12, face = "bold"),
          legend.text = element_text(size = 10),
          legend.background = element_blank(),
          legend.key = element_blank(),
          plot.title = element_text(face = "bold", size = 14)
    ) +
    
    # Add title with shark ID
    ggtitle(paste("Shark ID:", shark_id))
  
  # Save with shark ID in filename
  ggsave(file.path(individual_maps_dir, paste0("Movement_tracks_CTCRW_", tracktype, "_", sp_clean, "_", shark_id, ".tiff")), 
         plot = ind_trackmap,
         width = 21, height = 21, units = "cm", dpi = 300)
  
  # Optional: print progress
  print(paste("Saved map for shark:", shark_id))
}

# END OF CODE ----
# FOR NOW - TO BE CONTINUED

# TODO LIST ----
# TODO 1: figure out north arrow and scalebar properly
# TODO 2: figure out federal/state water boundary for map
# TODO 3: colour individual tracks individually

