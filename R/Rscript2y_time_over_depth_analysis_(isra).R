# TOP OF CODE ----
### ====================================================================================================
### Project:    Large bodied hammerhead complex in the western North Atlantic
### Analysis:   Assocaite horizontal movement data with bathymetry data
### Script:     ~2026_IUCN_ISRA_NWA/R/Rscript2y_time_over_depth_analysis_(isra).R
### Author:     Vital Heim
### Version:    1.0
### ====================================================================================================

### ....................................................................................................
### Content: this code allows to append bathymetry information from the ETOPO 2022 dataset hosted on
### the NOAA server to coordiantes of satellite-derived movement data
### ....................................................................................................

### ....................................................................................................
### [A] Ready environment, load packages ----
### ....................................................................................................

# A1: clear memory ----

## clear memory
rm(list = ls())

## if you are working on a network server
YOUR_IP <- "192.168.1.144" # add your IP address, server name or similar if you connect via a shared drive

## remove potential tmp.tif file from previous session
if (file.exists("tmp.tif")) file.remove("tmp.tif")

# A2: install and load necessary packages ----

## if first time or needed updates
# install.packages("tidyverse")
# install.packages("magrittr")
# install.packages("lubridate")
# install.packages("sp")
# install.packages("sf")
# install.packages("dplyr")
# install.packages("beepr")
# install.packages("raster")
# install.packages("terra")
# install.packages("reshape2")
# install.packages("maptools")
# install.packages("data.table")
# install.packages("ggplot2")
# install.packages("shape")
# install.packages("fields")
# install.packages("marmap")

## load packages and source needed functions
library(tidyverse)
library(magrittr)
library(lubridate)
library(sp)
library(sf)
library(dplyr)
library(beepr)
library(raster)
library(terra)
library(reshape2)
library(data.table)
library(ggplot2)
library(shape)
library(fields)
library(marmap)

# A3: Specify needed functions ----

## *A3.1: extract lat/lon extents from your movement data
get_bbox <- function(data, lat_col = "lat", lon_col = "lon", buffer = 1) {
  list(
    lon1 = min(data[[lon_col]], na.rm = TRUE) - buffer,
    lon2 = max(data[[lon_col]], na.rm = TRUE) + buffer,
    lat1 = min(data[[lat_col]], na.rm = TRUE) - buffer,
    lat2 = max(data[[lat_col]], na.rm = TRUE) + buffer
  )
}

# A4: Specify data and saveloc ----

## Project folder
projloc <- file.path("/",YOUR_IP, "Science","Projects_current", "2026_IUCN_ISRA_NWA")

## Input data
dataloc <- file.path("/",YOUR_IP, "Science","Projects_current", "2026_IUCN_ISRA_NWA","Data_input")

## Output data
saveloc <- file.path("/", YOUR_IP,"Science","Projects_current", "2026_IUCN_ISRA_NWA","Data_output", "Time_over_depth") # Adjust this
### check if saveloc already present
if (!dir.exists(saveloc)) {
  dir.create(saveloc, recursive = TRUE)
  cat("New saveloc was created at:", saveloc, "\n")
} else {
  cat("Saveloc is already present","\n")
}

# A5: Define universal variables (e.g. for plotting) ----

## speciescolours
slewcol <- "#EDA904"
smokcol <- "#70AB27"
szygcol <- "#E26306"

# A6: Define universal options, variables, etc. (e.g. for plotting) ----

options(timeout = 3000) # manually increase time out threshold (needed when downloading basemap)
options(scipen=999) # so that R doesn't act up for pit numbers  
options(warn=2) #set this to two if you have warnings that need to be addressed as errors via traceback()
options(error = function() beep(9))  # give warning noise if it fails

### ....................................................................................................
### [B] Data import ----
### ....................................................................................................

# B1: Data import ----

## movement data - filtered and standardised
hammers <- readRDS(file = file.path(dataloc, "Time_over_depth", "Data_aniMotum_CRW_output_segmented_rerouted_proj_WGS84_converted_with_coord_CIs_with_Argosfilter_data.rds")) |>
  mutate(shark = as.numeric(str_sub(id, start = 1, end = str_locate(id, "\\_")[,1] - 1)))

## metadata for species
meta <- read_csv(file.path(dataloc,"Telemetry", "Data_Sphyrna_SPOT_tags_metadata.csv"))|>
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
  dplyr::rename(
    shark = ptt_id)

### deal with local time stamp for deployment
# meta$datetime_deployment <- as.POSIXct(meta$datetime_deployment_local,format="%Y/%m/%d %H:%M:%S",tz="US/Eastern") # Cabo follows the mountain standard time (UTC-7), that does not have a time change. You can specify as "MST" or as 
# attr(meta$datetime_deployment, "tzone") <- "UTC"

## combine movement data with metadata
hammers %<>% left_join(meta)  # , by = join_by(shark == id) # doesn't work naming columns, has gotten worse.

## VERY IMPORTANT: The movement data needs to have time stamps in ascending order,
## double check that this is true by ordering the dataframe
hammers %<>%
  dplyr::arrange(
    shark,
    date
  ) %>%
  dplyr::select( # and since we are at it, the old id column can be removed
    shark,
    date,
    lon,
    lat,
    stl,
    sex,
    species,
    -id
  )

### ....................................................................................................
### [C] Extract bathymetry data  ----
### ....................................................................................................

# C1: define spatial extent of your movement data ----

bbox <- get_bbox(hammers, buffer = 1)  # buffer adds padding in degrees

# C2: obtain bathymetry data from NOAA server ----
setwd(tempdir())  # switch to temp directory where R always has permissions - othwerwise code below acting up
bathy <- marmap::getNOAA.bathy(lon1 = bbox$lon1, lon2 = bbox$lon2,
                               lat1 = bbox$lat1, lat2 = bbox$lat2,
                               resolution = 1) # resolution of the grid, in minutes (default is 4)
setwd("C:/Users/vital/SharktankDrive/Project_code/2026_IUCN_ISRA_NWA") #setwd back

# C3: explore bathymetry data ----

## summary
# summary.bathy(bathy)

## visualisation
# plot.bathy(bathy, image = T, land = T)
# scaleBathy(bathy, deg = 2, x = "bottomleft", inset = 5)

### ..........................................................................................
### [D] Combine bathymetry data with movement data ----
### ..........................................................................................

## the get_depth() function can be used to retrieve depth information by either
## clicking on the map or by providing a set of longitude/latitude pairs. This is
## helpfullto get depth information along a GPS track record for instance. If the
## argument distance is set to TRUE, the haversine distance (in km) from the first
## data point on will also be computed

# D1: overview of tracking data on bathymetry raster

plot(bathy, image = TRUE, land = TRUE, n=1)
## add location estimates
points(hammers$lon, hammers$lat, pch = 21, col = "black",
       bg = "yellow", cex = 1.3)

# D2: associate depth info with tracking points

hammers_depth <- hammers %>%
  dplyr::bind_cols(marmap::get.depth(bathy, hammers[, c("lon", "lat")],
                                     locator  = FALSE,
                                     distance = TRUE) %>% dplyr::select(depth, dist.km)
  ) %>%
  dplyr::mutate( # due to resolution of bathymetry we have some locations that are >0, which obvsly makes no sense, set them to NA
    depth = dplyr::if_else(depth >= 0, NA, depth)
  )

## visualise data points with associated depth on bahtymetry contur map
### create a contour plot for the bathymetry and add a scale
# tiff(file.path(saveloc, "Tracks_w_depth_info_coloured_all.tiff"), width = 24, height = 18, units = "cm", res = 300)
# par(mai=c(1,1,1,1.5))
# plot(bathy, lwd = c(0.3, 1), lty = c(1, 1),
#      deep = c(-4500, 0), shallow = c(-50, 0), step = c(500, 0),
#      col = c("grey", "black"), drawlabels = c(FALSE, FALSE))
# ### add scale bar
# marmap::scaleBathy(bathy, deg = 2, x = "bottomleft", inset = 5)
# ### set color palette for depth
# mx <- ceiling(abs(min(hammers_depth$depth, na.rm = TRUE)))
# col.points <- femmecol(mx)
# col.points[is.na(hammers_depth$depth)] <- "black"
# # plot points and color depth scale
# points(hammers_depth[,3:4], col = "black", bg = col.points[abs(hammers_depth$depth)],
#        pch = 21, cex = 1.5)
# colorlegend(zlim = c(mx, 0), col = rev(col.points),
#             main = "Depth [m]", posx = c(0.85, 0.88))
# dev.off()

## above code is a pain in the neck with colour legend, try to rescale depth to a smaller palette.
tiff(file.path(saveloc, "Tracks_w_depth_info_coloured_all_palette_corrected.tiff"), width = 24, height = 18, units = "cm", res = 300)
par(mai=c(1,1,1,1.5))
plot(bathy, lwd = c(0.3, 1), lty = c(1, 1),
     deep = c(-4500, 0), shallow = c(-50, 0), step = c(500, 0),
     col = c("grey", "black"), drawlabels = c(FALSE, FALSE))
### add scale bar
marmap::scaleBathy(bathy, deg = 2, x = "bottomleft", inset = 5)
### set color palette for depth
mx <- ceiling(abs(min(hammers_depth$depth, na.rm = TRUE)))
n_cols <- 500  # number of color steps
col.points <- femmecol(n_cols)
## Map depths to color indices within 1:n_cols
depth_vals <- abs(hammers_depth$depth)
col_indices <- round((depth_vals / mx) * (n_cols - 1)) + 1
col_indices <- pmin(pmax(col_indices, 1), n_cols)  # clamp to valid range

bg_cols <- col.points[col_indices]
bg_cols[is.na(hammers_depth$depth)] <- "black"
# plot points and color depth scale
points(hammers_depth[,3:4], col = "black", bg = bg_cols,
       pch = 21, cex = 1.5)
colorlegend(zlim = c(mx, 0), col = rev(col.points),
            main = expression(bold("Depth [m]")), posx = c(0.85, 0.88))
dev.off()

### ....................................................................................................
### [E] Summarise time over depth data by species ----
### ....................................................................................................

# E1: mean depth across individuals, sex and population

## individual-level
mean_depth_ind <- hammers_depth %>%
  dplyr::group_by(
    shark
  ) %>%
  dplyr::summarise(
    mean_depth = mean(depth, na.rm = T),
    sd_depth = sd(depth, na.rm = T),
    max_depth = min(depth, na.rm = T),
    min_depth = max(depth, na.rm = T)
  )

write.csv(mean_depth_ind, file.path(saveloc, "Summary_mean_depth_by_individual.csv"), row.names = F)

## sex-level
# mean_depth_sex <- hammers_depth %>%
#   dplyr::group_by(
#     sex
#   ) %>%
#   dplyr::summarise(
#     mean_depth = mean(depth, na.rm = T),
#     sd_depth = sd(depth, na.rm = T),
#     max_depth = min(depth, na.rm = T),
#     min_depth = max(depth, na.rm = T)
#   )
# 
# write.csv(mean_depth_sex, file.path(saveloc, "Summary_mean_depth_by_sex.csv"), row.names = F)

## overall
mean_depth_species <- hammers_depth %>%
  dplyr::group_by(
    species
  ) %>%
  dplyr::summarise(
    mean_depth = mean(depth, na.rm = T),
    sd_depth = sd(depth, na.rm = T),
    max_depth = min(depth, na.rm = T),
    min_depth = max(depth, na.rm = T)
  )

write.csv(mean_depth_species, file.path(saveloc, "Summary_mean_depth_by_species.csv"), row.names = F)

# E2: ToD - % time over depth approach ----

## remove NA values in your depth values
# sum(is.na(hammers_depth$depth)) #233
hammers_depth_nona <- hammers_depth %>% filter(depth >= min(hammers_depth$depth, na.rm = T) & depth <= 0)

## Define depth bin breaks and clean labels
# breaks <- seq(0, -3500, by = -100)
# labels <- paste0(abs(breaks[-length(breaks)]), "-", abs(breaks[-1]), "m")
breaks <- c(seq(0, -1000, by = -100), seq(-1500, -3500, by = -500))
labels <- paste0(abs(breaks[-length(breaks)]), "-", abs(breaks[-1]), "m"); labels

## Create depth bin column in df
hammers_depth_nona <- hammers_depth_nona %>%
  mutate(depth_bin = cut(abs(depth),
                         breaks = abs(breaks),
                         labels = labels,
                         include.lowest = T,
                         right = F)) # for depth bins, left-closed [0, 100) makes more sense because: depth of exactly 100m logically belongs in the "100-200m" bin

## Calcualte % in each bin and extend all ID x bin combinations, fill missing with 0
## we can either calculate the species mean directly, but this treats each detection as equal, 
##regardless of individual. This is fine if you want to know "of all detections of species X, 
## what % occurred in each depth bin?"
## However, Calculating individual means first (two-step approach) is more statistically appropriate if 
## individuals vary a lot in sampling effort — e.g., one shark was detected 1000 times and another only 10 times. 
## In that case, the heavily detected individual dominates the species-level pattern. The two-step would be:

tod_ind <- hammers_depth_nona %>%
  dplyr::group_by(
    shark,
    species,
    depth_bin
  ) %>%
  dplyr::summarise(
    n = n(),
    .groups = "drop"
  ) %>%
  tidyr::complete(
    nesting(
      shark, 
      species),
    depth_bin,
    fill = list(n=0)
  ) %>%
  dplyr::group_by(
    shark
  ) %>%
  dplyr::mutate(
    tod_pct = n/sum(n)*100
    ) %>%
  dplyr::ungroup()
  
write.csv(tod_ind, file.path(saveloc, "Summary_time_over_depth_by_individual.csv"), row.names =F)

tod_hammers <- tod_ind %>%
  dplyr::group_by(
    species, 
    depth_bin
    ) %>%
  dplyr::summarise(
    tod_mean = mean(tod_pct),
    tod_median = median(tod_pct),
    tod_sd = sd(tod_pct),
    tod_min = min(tod_pct),
    tod_max = max(tod_pct),
    n_individuals = n_distinct(shark),
    .groups = "drop"
  )

write.csv(tod_hammers, file.path(saveloc, "Summary_time_over_depth_by_species.csv"), row.names =F)

### ....................................................................................................
### [F] Visualise your animals' time over depth patterns ----
### ....................................................................................................

# F1: plot timeseries plot of location over depth by individual ----

df_plot <- hammers_depth_nona

## define colour palette anew
color.points <- colorRampPalette(c("#00008F", "#0000FF", "#00FFFF",
                                   "#FFFF00", "#FF0000", "#800000"))(3500)

## define your unique ids
ids <- sort(unique(df_plot$shark))

for (i in seq_along(ids)) {
  
  ## create new directory for plots
  dir.create(file.path(saveloc, "individual_timeseries_plots"), showWarnings = FALSE)
  
  shark_id <- ids[i]
  
  # Progress indicator
  cat(sprintf("[%d/%d] Processing shark %s... ", i, length(ids), shark_id))
  
  ## filter individual
  df_ind <- df_plot %>% filter(shark == ids[i])
  
  ## create species variable for filename
  spp_ind <- gsub("\\.", "", unique(df_ind$species))
  
  ## Open plot window 
  tiff(file.path(saveloc, "individual_timeseries_plots", paste0("Timeseries_plot_ToD_",shark_id, "_", spp_ind, "_", Sys.Date(),".tiff")), width =12, height = 6, units = "cm", res = 300)
  
  ## We only need 1 column per individual as no sub-groupings (e.g. season)
  par(mfrow = c(1, 1))
  
  ## define plot margins
  par(mar = c(3,3.5,.5,1.15), # A numeric vector of length 4, which sets the margin sizes in the following order: bottom, left, top, and right. The default is c(5.1, 4.1, 4.1, 2.1).
      mgp = c(2, 0.3,0), # A numeric vector of length 3, which sets the axis label locations relative to the edge of the inner plot window. The first value represents the location the labels (i.e. xlab and ylab in plot), the second the tick-mark labels, and third the tick marks. The default is c(3, 1, 0).
      las = 1) # numeric value indicating the orientation of the tick mark labels and any other text added to a plot after its initialization. The options are as follows: always parallel to the axis (the default, 0), always horizontal (1), always perpendicular to the axis (2), and always vertical (3).
  
  ## plot your data
  plot(df_ind$date, abs(df_ind$depth),
       ylim = c(3500, 0),
       type = "n",            # empty plot first
       xlab = "",
       ylab = "",
       cex.axis = .65,
       # main = paste("Mako", ids[i]),
       axes = FALSE)
  
  ## reference depth lines
  abline(h = c(0, 100, 200, 300, 400, 500, 1000, 1500, 2000, 2500, 3000, 3500),
         col = "grey80", lwd = 0.4, lty = 2)
  
  ## timeseries location over depth  - line
  lines(df_ind$date, abs(df_ind$depth), col = "grey20", lwd = 0.5)
  
  ## timeseries location over depth - points
  points(df_ind$date, abs(df_ind$depth),
         pch = 19, cex = 0.5,
         col = color.points[round(abs(df_ind$depth)) + 1])
  
  ## label vertical axis
  axis(2, at = seq(0, 3500, by = 500), labels = seq(0, 3500, by = 500),
       las = 1, cex.axis = 0.55, tck = -0.035)
  mtext(expression(bold("Depth [m]")), side = 2, line = 2.6, cex = 0.65, las = 0)
  
  ## label horizontal axis
  date_ticks <- pretty(df_ind$date, n = 4)
  date_ticks <- date_ticks[date_ticks >= min(df_ind$date) & date_ticks <= max(df_ind$date)]
  date_labels <- format(date_ticks, "%b-%y")
  date_ticks <- date_ticks[!duplicated(date_labels)]
  date_labels <- unique(date_labels)
  axis.POSIXct(1, at = date_ticks, labels = date_labels, las = 1, cex.axis = 0.55, tck = -0.035)
  axis.POSIXct(1, at = c(min(df_ind$date), max(df_ind$date)),
               labels = c("Start track", "End track"),
               las = 1, cex.axis = 0.55, tck = -0.035)
  mtext(expression(bold("Tracking period")), side = 1, line = 1.6, cex = 0.65)
  
  ## add legend - doesnt work yet
  # image.plot(legend.only = TRUE,
  #            zlim = c(0, 4200),
  #            col = color.points,
  #            legend.shrink = 0.5,
  #            legend.width = 0.8,
  #            # smallplot = c(0.05, 0.08, 0.1, 0.4),  # bottom left: c(x1, x2, y1, y2) in npc
  #            axis.args = list(cex.axis = 0.7, at = seq(0, 4200, by = 1000)),
  #            legend.args = list(text = "Depth (m)", cex = 0.8, side = 2, line = 1.5)
  # )
  
  box()
  
  dev.off()
}

# F2: plot time over depth histograms by species ----
## (individual-level plotting tbc)

# D1: plot TAT horizontally for each species

## plot for each species individuall
for (spp_f in unique(tod_hammers$species)) {
  
  ## prepare new directory for plots
  dir.create(file.path(saveloc, "tod_species_level"), showWarnings = FALSE)
  
  ## Create plot df
  tod_plot <- tod_hammers %>%
  dplyr::filter(
    species %in% spp_f
  )
  ### need to reverse the df as we want descending depths not ascending depth
  tod_plot <- tod_plot[nrow(tod_plot):1, ]

  ## define depth bin labels
  plot_labels <- gsub("m", "", rev(levels(tod_plot$depth_bin)))
  
  ## add a spp label for filename without "."
  spp <- gsub("\\.", "", spp_f) # createa a character string to name files when saving, "." might cause issues downstream

  ## Open plot window
  tiff(file.path(saveloc, "tod_species_level", paste0("TOD_Histo_", spp, "_", Sys.Date(),".tiff")), width = 14, height = 10, units = "cm", res = 300)

  ## Define plot margins
  ## Margins should be minimal on the left, and top, close to default bottom and right
  par(mar = c(3,3.5,1.5,.75), # A numeric vector of length 4, which sets the margin sizes in the following order: bottom, left, top, and right. The default is c(5.1, 4.1, 4.1, 2.1).
      mgp = c(2, 0.75,0), # A numeric vector of length 3, which sets the axis label locations relative to the edge of the inner plot window. The first value represents the location the labels (i.e. xlab and ylab in plot), the second the tick-mark labels, and third the tick marks. The default is c(3, 1, 0).
      las = 1) # numeric value indicating the orientation of the tick mark labels and any other text added to a plot after its initialization. The options are as follows: always parallel to the axis (the default, 0), always horizontal (1), always perpendicular to the axis (2), and always vertical (3).

  ## Barplot
  bp_tod <- barplot(height = tod_plot$tod_mean,
                  , horiz = T
                  , xlim = c(0,100)
                  , las = 1
                  , xlab = "", ylab = ""
                  #, xaxt = "n"
                  , yaxt = "n"
                  , cex.axis = 0.65
                  , border = "black"
                  , col = ifelse(tod_plot$species == "S.lewini", slewcol,
                                 ifelse(tod_plot$species == "S.mokarran", smokcol,
                                        ifelse(tod_plot$species == "S.zygaena", szygcol, "black"))))

  ## add error bars
  # segments(tod_plot$tod_mean - tod_plot$tod_sd, bp_tod, tod_plot$tod_mean + tod_plot$tod_sd , bp_tod,
  #          lwd = 1)
  arrows(tod_plot$tod_mean - tod_plot$tod_sd, bp_tod, 
         tod_plot$tod_mean + tod_plot$tod_sd , bp_tod,
         lwd = 1, angle = 90,
         code = 3, length = 0.05)

  ## label vertical axis
  axis(2, at = bp_tod, labels = plot_labels, cex.axis = 0.65, tck = -0.035, las = 1) # tickmarks and labels

  ## label horizontal axis
  mtext(expression(bold("[%] of time")), side = 1, line = 1.6, cex = 0.65)

  ## Add unit of vertical axis above and centered
  mtext(expression(bold("Depth [m]")), side = 3, line = -0.5, at = -5.5, cex = 0.65)

  ## save it
  dev.off()

  cat(sprintf("Saved plot for %s\n", spp_f))
}

# F3: plot movement tracks coloured by location over depth ----

## prepare depth bin colours for more aesthetically pleasing plotting
depth_breaks_map <- c(0, 100, 200, 300, 400, 500, Inf)
depth_labels_map <- c("0-100", "100-200", "200-300", "300-400", "400-500", "500+")

## Use shades of blue for depth bin colours
depth_cols <- c("#E0F3FF", "#9ECAE1", "#4292C6", "#2171B5", "#08519C", "#08306B")

## Assign locations with depth bin value with corresponding colour
hammers_depth_nona <- hammers_depth_nona %>%
  mutate(
    depth_bin_map = cut(abs(depth),
                    breaks = depth_breaks_map,
                    labels = depth_labels_map,
                    include.lowest = T, right = F), # for depth bins, left-closed [0, 100) makes more sense because: depth of exactly 100m logically belongs in the "100-200m" bin
    pt_col = depth_cols[as.integer(depth_bin_map)]
  )

## *F3.1.: individual-level location over depth maps ----
for (shark_id in unique(hammers_depth_nona$shark)) {
  
  ## prepare new directory for plots
  dir.create(file.path(saveloc, "Map_loc_over_depth_individual"), showWarnings = FALSE)
  
  ## filter individual data
  df_map_ind <- hammers_depth_nona %>% filter(shark == shark_id)
  
  ## prepare species variable for file name
  spp    <- gsub("\\.", "", unique(df_map_ind$species))
  
  ## open device
  tiff(file.path(saveloc, "Map_loc_over_depth_individual", paste0("Map_location_over_depth_", shark_id, "_", spp, "_", Sys.Date(), ".tiff")),
       width = 12, height = 8, units = "cm", res = 300)
  
  ## define plot margins
  par(mar = c(3,3,1.5,.75), # A numeric vector of length 4, which sets the margin sizes in the following order: bottom, left, top, and right. The default is c(5.1, 4.1, 4.1, 2.1).
      mgp = c(2, 0.75,0), # A numeric vector of length 3, which sets the axis label locations relative to the edge of the inner plot window. The first value represents the location the labels (i.e. xlab and ylab in plot), the second the tick-mark labels, and third the tick marks. The default is c(3, 1, 0).
      las = 1)
  
  ## plot bathymetry map with depth contours
  plot(bathy, lwd = c(0.3, 1), lty = c(1, 1),
       deep = c(-1000, 0), shallow = c(-200, 0), step = c(100, 0), # plot isobaths from 200m to 1000m in 200m steps
       col = c("grey", "black"), drawlabels = c(F, F),
       las = 1, cex.axis = .65,
       xlab = "", ylab = "")
  
  ## add a scale bar
  marmap::scaleBathy(bathy, deg = 3, x = "bottomleft", inset = 3, cex = .65)
  
  ## add location estimates of your animal
  points(df_map_ind[, 3:4], col = "black", bg = df_map_ind$pt_col,
         pch = 21, cex = .75)
  
  ## add a depth legend
  legend("bottomright",
         legend = depth_labels_map,
         pt.bg = depth_cols,
         pch = 21, col = "black",
         title = expression(bold("Depth [m]")),
         cex = 0.7, bty = "n", pt.cex = 1.3)
  
  ## label the x axis
  mtext(expression(bold("Longitude [°W]")), side = 1, line = 1.6, cex = 0.65)
  
  ## label the y axis
  mtext(expression(bold("Latitude [°N]")), side = 2, line = 1.6, cex = 0.65, las = 0)
  
  ## add a title
  title(main = bquote(bolditalic(.(paste(unique(df_map_ind$species), "-", shark_id)))), 
        cex.main = 0.9, adj = 0)
  
  ## save your map
  dev.off()
  
  ## print progress
  cat(sprintf("Saved individual map: %s\n", shark_id))
}

## *F3.2.: species-level location over depth maps ----
for (spp.f in unique(hammers_depth_nona$species)) {
  
  ## prepare new directory for plots
  dir.create(file.path(saveloc, "Map_loc_over_depth_species"), showWarnings = FALSE)
  
  ## filter species data
  df_map_spp <- hammers_depth_nona %>% filter(species == spp.f)
  
  ## prepare species variable for filename
  spp    <- gsub("\\.", "", spp.f)
  
  tiff(file.path(saveloc, "Map_loc_over_depth_species", paste0("Map_location_over_depth_", spp, "_", Sys.Date(), ".tiff")),
       width = 12, height = 8, units = "cm", res = 300)
  
  # define plot margins
  par(mar = c(3,3,1.5,.75), # A numeric vector of length 4, which sets the margin sizes in the following order: bottom, left, top, and right. The default is c(5.1, 4.1, 4.1, 2.1).
      mgp = c(2, 0.75,0), # A numeric vector of length 3, which sets the axis label locations relative to the edge of the inner plot window. The first value represents the location the labels (i.e. xlab and ylab in plot), the second the tick-mark labels, and third the tick marks. The default is c(3, 1, 0).
      las = 1)
  
  ## plot bathymetry map with depth contours
  plot(bathy, lwd = c(0.3, 1), lty = c(1, 1),
       deep = c(-1000, 0), shallow = c(-200, 0), step = c(100, 0), # plot isobaths from 200m to 1000m in 200m steps
       col = c("grey", "black"), drawlabels = c(FALSE, FALSE),
       cex.axis = .65, ylab ="", xlab = "")
  
  ## add scalebar
  marmap::scaleBathy(bathy, deg = 3, x = "bottomleft", inset = 3, cex = .65)
  
  ## add location estimates of your species
  points(df_map_spp[, 3:4], col = "black", bg = df_map_spp$pt_col,
         pch = 21, cex = .75)
  
  ## add a depth legend
  legend("bottomright",
         legend = depth_labels_map,
         pt.bg = depth_cols,
         pch = 21, col = "black",
         title = expression(bold("Depth [m]")),
         cex = 0.7, bty = "n", pt.cex = 1.3)
  
  ## label the x axis
  mtext(expression(bold("Longitude [°W]")), side = 1, line = 1.6, cex = 0.65)
  
  ## label the y axis
  mtext(expression(bold("Latitude [°N]")), side = 2, line = 1.6, cex = 0.65, las = 0)
  
  ## add a title
  title(main = bquote(bolditalic(.(paste(unique(df_map_spp$species))))), 
        cex.main = 0.9, adj = 0)
  
  ## save your map
  dev.off()
  
  ## print your progress
  cat(sprintf("Saved species map: %s\n", spp.f))
}

# END OF CODE ----

## TODO: individual TOD histogram plots