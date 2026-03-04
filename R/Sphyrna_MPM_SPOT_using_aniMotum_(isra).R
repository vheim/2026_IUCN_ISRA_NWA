# TOP OF CODE ----
### ====================================================================================================
### Project:    Large bodied hammerhead complex in the western North Atlantic
### Analysis:   Processing and cleaning satellite telemetry data of fin-mounted SPOT tags for further steps
### Script:     ~2026_IUCN_ISRA_NWA/R/Rscript2_CTCRW_SSM_SPOT_using_aniMotum_(isra).R
### Author:     Vital Heim
### Version:    1.0
### ====================================================================================================

### ....................................................................................................
### Content: this R script contains the code to filter, clean and process raw data collected by fin-
###          mounted Smart Position and Temperature (SPOT) transmitters.
###          This script contains the code fit a continuous time correlated random walk model (CTCRW
###          in a state space model (SSM) framework using sda-prefiltered data and to fit a 
###          move persistence model to obtain estimates of a continuous-valued behavioural index along 
###          individual tracks (Auger-Méthé et al. 2017; Jonsen et al. 2019).
###
###          NOTE: this script is under construction and specifically tailored to MPM fitting and 
###          extraction. This means that if you run this script without running Rscript2 in this
###          project first, additional steps (e.g. validadation, residual plots, etc.) still need to be added.
###        
###          The code is based on the functions provided within the animotum (Jonsen et al. 2023) package:
###          > citation("aniMotum")
###          Ian Jonsen, W James Grecian, Lachlan Phillips, Gemma Carroll, Clive R. McMahon,
###          Robert G. Harcourt, Mark A. Hindell, and Toby A. Patterson (2023) aniMotum, an R
###          package for animal movement data: Rapid quality control, behavioural estimation
###          and simulation.  Methods in Ecology and Evolution DOI: 10.1111/2041-210X.14060
###          The script contains additional code for pre-fitting processing and cleaning of the data.
### ....................................................................................................

### ....................................................................................................
### [A] Setwd, paths and parameters ----
### ....................................................................................................

# A1: clear memory ----

rm(list = ls())

## Detach all packages except base ones
# lapply(paste0('package:', names(sessionInfo()$otherPkgs)),
# detach, character.only = TRUE, unload = TRUE, force = TRUE)

# A2: load necessary packages ----

# library("remotes")
## install if first time
# install.packages("tidyverse")
# install.packages("magrittr")
# install.packages("TMB", type = "source") # if package version inconsistency detected during loading of aniMotum
# install.packages("aniMotum",
# repos = c("https://cloud.r-project.org",
# "https://ianjonsen.r-universe.dev"),
# dependencies = TRUE) #foiegras was removed from CRAN and replaced with Animotum on 12-12-2022
# install.packages("patchwork")
# install.packages("sf")
# install.packages("sp")
# install.packages("rnaturalearth")
# install.packages("rnaturalearthdata")
# remotes::install_github("ropensci/rnaturalearthhires")
# install.packages("ggplot2")
# install.packages("xts")
# install.packages("trip")
# install.packages("ggspatial")
# install.packages("virids")

## load
library(tidyverse)
library(magrittr)
library(aniMotum)
library(patchwork)
library(sf)
library(sp)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)
library(ggplot2)
# library(xts)
# library(trip)
library(ggspatial)
library(viridis)

# A3: Specify needed functions

## Data import: import speedfiltered data and validate that only the most current data is used.
#' @param dataloc Character string specifying the directory path
#' @param speedfilter Character string for the file pattern prefix, i.e. the sda-filter used in previous scripts
#'
#' @return A data frame containing the imported data if successful, NULL otherwise

import_sda_data <- function(dataloc, speedfilter) {
  
  # Check number of data files in directory
  data_files <- list.files(path = dataloc,
                           pattern = paste0(speedfilter, "_.+.R"),
                           full.names = FALSE)
  
  n_files <- length(data_files)
  
  if (n_files == 0) {
    stop("No data files found in directory: ", dataloc)
    
  } else if (n_files == 1) {
    # Print success message in green
    cat("\033[32m✓ Single file available and imported:", data_files[1], "\033[0m\n")
    
    # Run the import code
    mydets_f <- list.files(path = dataloc,
                           pattern = paste0(speedfilter, "_.+.R"),
                           full.names = TRUE) %>%
      purrr::map_dfr(readRDS)
    
    return(mydets_f)
    
  } else {
    # Print warning in red
    cat("\033[31m⚠ WARNING: Multiple files found in the import folder!\033[0m\n")
    cat("\033[31mFiles detected:\033[0m\n")
    
    # Get full file info including modification times
    full_paths <- list.files(path = dataloc,
                             pattern = paste0(speedfilter, "_.+.R"),
                             full.names = TRUE)
    
    file_info <- file.info(full_paths)
    file_info$filename <- basename(rownames(file_info))
    file_info <- file_info[order(file_info$mtime, decreasing = TRUE), ]
    
    # List all files with timestamps
    for (i in 1:nrow(file_info)) {
      timestamp <- format(file_info$mtime[i], "%Y-%m-%d %H:%M:%S")
      if (i == 1) {
        cat("\033[31m  -", file_info$filename[i], "(", timestamp, ") [MOST RECENT]\033[0m\n")
      } else {
        cat("\033[31m  -", file_info$filename[i], "(", timestamp, ")\033[0m\n")
      }
    }
    
    # Get most recent file
    most_recent_file <- rownames(file_info)[1]
    most_recent_name <- basename(most_recent_file)
    most_recent_time <- format(file_info$mtime[1], "%Y-%m-%d %H:%M:%S")
    
    cat("\033[33m→ The most recent file is:", most_recent_name, "(", most_recent_time, ")\033[0m\n")
    
    # Ask for confirmation
    response <- readline(prompt = "Do you want to import this file? (y/n): ")
    
    if (tolower(trimws(response)) == "y") {
      cat("\033[32m✓ Importing:", most_recent_name, "\033[0m\n")
      
      # Import only the most recent file
      mydets_f <- readRDS(most_recent_file)
      
      return(mydets_f)
    } else {
      cat("\033[31m✗ Import cancelled by user.\033[0m\n")
      return(NULL)
    }
  }
}

# A4: Specify data and saveloc ----

YOUR_IP <- "NA" # add your IP address, server name or similar if you connect via a shared drive

## Project folder
projloc <- file.path("/",YOUR_IP, "Science","Projects_current", "2026_IUCN_ISRA_NWA")

## Input data
dataloc <- file.path("/",YOUR_IP, "Science","Projects_current", "2026_IUCN_ISRA_NWA","Data_input")

## Output data
saveloc <- file.path("/", YOUR_IP,"Science","Projects_current", "2026_IUCN_ISRA_NWA","Data_output", "MPM") # Adjust this
### check if saveloc already present
if (!dir.exists(saveloc)) {
  dir.create(saveloc, recursive = TRUE)
  cat("New saveloc was created at:", saveloc, "\n")
} else {
  cat("Saveloc is already present","\n")
}

# A5: Define universal variables (e.g. for plotting) ----

# NA

# A6: Define universal options, variables, etc. (e.g. for plotting) ----

options(scipen=999) # so that R doesn't act up for pit numbers  
options(rosm.cachedir = NULL) # prevents large cache folders to be created when plotting map tiles

### ....................................................................................................
### [B] SSM with prefiltered data ----
### ....................................................................................................

# Remember: if you use prefiltered data, this means raw argos data was filtered to remove
# spurious detections based on  a speed-distance-angle algorithm. Based on the script/file either
# the sda filter from the package "argosfilter" or "trip" was used.
# In either case, the files were already adjusted for near duplicated observations so this does not
# need to be done again here in section B.

# Define first which sda filter you used

speedfilter <- "Argosfilter" # choices: "Argosfilter" or "TripSDA"

# B1: Import the filtered data ----

# mydets_f <- list.files(path = dataloc,
#                        pattern = paste0(speedfilter,"_.+.R"),
#                        full.names = TRUE ) %>%
#   purrr::map_dfr(readRDS) # make sure that there is only the most up to date speedfilter file in the directory

## Important: there should only be a single rds file from the speedfilter script. Check this and make sure you only import the
## most current data.

mydets_f <- import_sda_data(dataloc = file.path(dataloc, "CTCRW"), speedfilter = speedfilter)

## check the data
# print(n = 1000, mydets_f %>% group_by(id) %>% dplyr::summarise(n = n()))
## A tibble: 69 × 2
# id         n
# <chr>  <int>
#   1 174515   129
# 2 175427   114
# 3 179468   798
# 4 179471   250
# 5 179472    38
# 6 180910   444
# 7 180912   876
# 8 180914   748
# 9 182746    36
# 10 183619  1048
# 11 183620   493
# 12 183621   537
# 13 183622   137
# 14 183624    72
# 15 198201   736
# 16 198202   481
# 17 198203   218
# 18 198204   217
# 19 198205   238
# 20 198206   299
# 21 210602   311
# 22 210603   789
# 23 210604    18
# 24 210605   264
# 25 210606   739
# 26 222134   258
# 27 222136    51
# 28 222137   877
# 29 222138   263
# 30 233536   724
# 31 233537   982
# 32 233538  1120
# 33 235278   466
# 34 235281   782
# 35 244597   687
# 36 244598   692
# 37 244600   616
# 38 244601  1197
# 39 244602   208
# 40 244603   103
# 41 244605  1451
# 42 244606  1991
# 43 244609   361
# 44 244610   788
# 45 244611   337
# 46 261302  2107
# 47 261303  1495
# 48 261304  1151
# 49 261305   619
# 50 261306   318
# 51 261309  1014
# 52 261742    82
# 53 264015  1609
# 54 264019   147
# 55 264021   356
# 56 264022   277
# 57 264023   262
# 58 264024   243
# 59 264025   480
# 60 264026   771
# 61 286887   275
# 62 286888   957
# 63 286889  1284
# 64 286890   119
# 65 286895   287
# 66 286896   409
# 67 286897   837
# 68 286898   146
# 69 286899   309

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
  dplyr::mutate(
    ptt_id = as.character(ptt_id)
  )|>
  dplyr::rename(
    id = ptt_id,
    FishLengthCm = stl)

## combine movement data with metadata
mydets_f %<>% left_join(meta)  # , by = join_by(shark == id) # doesn't work naming columns, has gotten worse.

## VERY IMPORTANT: The movement data needs to have time stamps in ascending order,
## double check that this is true by ordering the dataframe
mydets_f %<>%
  dplyr::arrange(
    id,
    date
  )

# B2: housekeeping, tidy data and define column classes ----

## the fit_ssm and functions in the animotum package expect a data.frame or tibble in
## the following format if (!) error elipse information is available. If the data was collected
## using other methods than the CLS Argos' Kalman filter model other formats might be needed.
## Check: https://ianjonsen.github.io/aniMotum/articles/Overview.html

## We need:
## 'id'
## 'date'
## 'loc' -location class
## 'lon'
## 'lat'
## 'smaj' - error semi-major axis
## 'smin' - error Semi-minor axis
## 'eor' - error ellipse orientation

## Basic housekeeping
det_f <-mydets_f %>%
  dplyr::select( # select relevant columns, here: id, date, location class (lc), lon, lat,
    id,
    date,
    lc,
    lon,
    lat,
    smaj,
    smin,
    eor,
    species
  )

## check if there are Z locations left (should not!)
print(paste0("There are ", length(which(det_f$lc == "Z")), " Z-locations left."))

# B4: Known/tagging locations ----

## we already added tagging locations in a previous step. If you have no FastLoc data
## and no GPS location classes, then the tagging locations can be filtered
## by going for lc == G. If there is GPS data this section needs to be updated using the
## type == "user" argument.

## If this is not the case and the smaj, smin, eor data is already present from the previous
## script this can be commented out

# det_f <- det_f %>%
#   dplyr::mutate(
#     smaj = ifelse(lc == 'G', 50,
#                   smaj),
#     smin = ifelse(lc == 'G', 50,
#                   smin),
#     eor = ifelse(lc == 'G', 0,
#                  eor)
#   )

# B5: Deal with parametrically logical but biologically unreasonable locations ----

## Only execute this section if above speedfilters do not exclude all spurious
## detections satisfyingly.

## Plot filtered tracks individually. Plot the tracks of the desired speedfilter.

### Basemap
# esri_ocean <- paste0('https://services.arcgisonline.com/arcgis/rest/services/',
#                      'Ocean/World_Ocean_Base/MapServer/tile/${z}/${y}/${x}.jpeg')
# ### Colors
# nb.cols <- length(unique(det_f$id))
# mycolors <- viridis(nb.cols)
# 
# ### Tracks
# sf_locs <- sf::st_as_sf(det_f, coords = c("lon","lat")) %>%
#   sf::st_set_crs(4326)
# 
# sf_lines <- sf_locs %>%
#   dplyr::arrange(id, date) %>%
#   sf::st_geometry() %>%
#   sf::st_cast("MULTIPOINT",ids = as.integer(as.factor(sf_locs$id))) %>%
#   sf::st_cast("MULTILINESTRING") %>%
#   sf::st_sf(id = as.factor(unique(det_f$id)))
# 
# sf_points <- sf_locs %>%
#   dplyr::arrange(id, date) %>%
#   sf::st_geometry() %>%
#   sf::st_cast("MULTIPOINT",ids = as.integer(as.factor(sf_locs$id))) %>%
#   sf::st_sf(id = as.factor(unique(det_f$id)))
# 
# for (i in 1:length(unique(det_f$id))){
# 
#   # plot the filtered locations by individual
#   ggplot() +
#     annotation_map_tile(type = esri_ocean,zoomin = 1,progress = "none") +
#     layer_spatial(sf_points[i,], size = 0.5) +
#     layer_spatial(sf_lines[i,], size = 0.75,aes(color = id)) +
#     #scale_x_continuous(expand = expansion(mult = c(.6, .6))) +
#     scale_colour_manual(values = mycolors[i]) +
#     theme() +
#     ggtitle("Speedfiltered (argosfilter::sdafilter()) Argos Location Paths",
#             subtitle = paste("Sphyrna spp. from the U.S. Atlantic (n = ", length(unique(det_f$id)), ")"))
# 
#   # save it
#   ggsave(file.path(saveloc,paste0("Speedfiltered_Argos_data_individual_", i, ".tiff")),
#          width = 21, height = 15, units = "cm", device ="tiff", dpi=300)
# }

## Notes code: if the code throws some errors/warning at you. However, there are not impacting
## the output, especially since the output only serves finding locations that are biologically
## unreasonable.

## Notes from comparing individual plots.
## Most individuals are fine. Some have locations on land in Andros, but given that
## the island is small in relation to available instrument accuray, these can be left in for now.
## Locations on land will later be dealt with by re-routing the tracks around land barriers

## 175427 has 1 spurious location at around -82.5°W and 28°N
## 2023-07-22 00:57:32 at lon. -36.06180 and lat: 34.09160 ###

## 179472 has two locations at around -85W
## one at 2019-11-16 00:03:57
## one at 2019-11-22 11:51:16

## 182746 has 1 location at -55 W
## at 2021-09-27 02:12:30

## 244611 has 1 location at -69.5 W

## 264024 has 1 location at -81.75W

## Filter out these segments manually
det_f <- det_f %>%
  filter(!(  
    #175427
    (id == "175427" & date == as.POSIXct("2025-02-28 15:18:10", tz = "UTC")) |
    #179472
    (id == "179472" & date == as.POSIXct("2019-11-16 00:03:57", tz = "UTC")) |
    (id == "179472" & date == as.POSIXct("2019-11-22 11:51:16", tz = "UTC")) |
    #182746
    (id == "182746" & date == as.POSIXct("2021-09-27 02:12:30", tz = "UTC")) |
    #235278
    (id == "235278" & date == as.POSIXct("2023-07-03 02:00:25", tz = "UTC")) |
    (id == "235278" & date == as.POSIXct("2023-07-02 01:11:49", tz = "UTC")) |
    (id == "235278" & date == as.POSIXct("2023-07-02 00:47:37", tz = "UTC")) |  
    #244611
    (id == "244611" & date == as.POSIXct("2024-02-17 15:55:20", tz = "UTC")) |
    #264024
    (id == "264024" & date == as.POSIXct("2025-02-14 22:40:38", tz = "UTC")) 
  ))

## Visualise manually improved tracks again:
#' sf_locs_clean <- sf::st_as_sf(det_f, coords = c("lon","lat")) %>%
#'   sf::st_set_crs(4326)
#' 
#' sf_lines_clean <- sf_locs_clean %>%
#'   dplyr::arrange(id, date) %>%
#'   sf::st_geometry() %>%
#'   sf::st_cast("MULTIPOINT",ids = as.integer(as.factor(det_f$id))) %>%
#'   sf::st_cast("MULTILINESTRING") %>%
#'   sf::st_sf(id = as.factor(unique(det_f$id)))
#' 
#' sf_points_clean <- sf_locs_clean %>%
#'   dplyr::arrange(id, date) %>%
#'   sf::st_geometry() %>%
#'   sf::st_cast("MULTIPOINT",ids = as.integer(as.factor(sf_locs_clean$id))) %>%
#'   sf::st_sf(id = as.factor(unique(det_f$id)))
#' 
#' #esri_ocean <- paste0('https://services.arcgisonline.com/arcgis/rest/services/',
#' #'Ocean/World_Ocean_Base/MapServer/tile/${z}/${y}/${x}.jpeg')
#' 
#' ### ALL
#' ggplot() +
#'   annotation_map_tile(type = esri_ocean,zoomin = 1,progress = "none") +
#'   layer_spatial(sf_points_clean, size = 0.5) +
#'   layer_spatial(sf_lines_clean, size = 0.75,aes(color = id)) +
#'   scale_x_continuous(expand = expansion(mult = c(.6, .6))) +
#'   scale_fill_manual(values = mycolors) +
#'   theme() +
#'   ggtitle("Speedfiltered (argosfilter::sdafilter()) and manually checked Argos Location Paths",
#'           subtitle = paste("Sphyrna spp. from the U.S. Atlantic (n = ", length(unique(det_f$id)), ")"))
#' ggsave(file.path(saveloc,"Speedfiltered_Argos_data_all_individuals_clean.tiff"),
#'        width = 21, height = 15, units = "cm", device ="tiff", dpi=150)
#' 
#' ### IND
#' for (i in 1:length(unique(det_f$id))){
#' 
#'   # plot the filtered locations by individual
#'   ggplot() +
#'     annotation_map_tile(type = esri_ocean,zoomin = 1,progress = "none") +
#'     layer_spatial(sf_points_clean[i,], size = 0.5) +
#'     layer_spatial(sf_lines_clean[i,], size = 0.75,aes(color = id)) +
#'     #scale_x_continuous(expand = expansion(mult = c(.6, .6))) +
#'     scale_colour_manual(values = mycolors[i]) +
#'     theme() +
#'     ggtitle("Speedfiltered (argosfilter::sdafilter()) and manually checked Argos Location Paths",
#'             subtitle = paste("Sphyrna spp. from the U.S. Atlantic (n = ", length(unique(det_f$id)), ")"))
#' 
#'   # save it
#'   ggsave(file.path(saveloc,paste0("Speedfiltered_Argos_data_individual_", i, "_clean.tiff")),
#'          width = 21, height = 15, units = "cm", device ="tiff", dpi=150)
#' # 


# B6: Calculate the time difference between detections in days and segment tracks ----

## Depending if all sharks survived and/or all tasg were deployed there might be NA values in the dataset
## check
check <- det_f[is.na(det_f$date),] # if these are NAs coming from post-release mortality sharks or undeployed tags, delete them

nrow(check) # if this is 0 you are good to continue

## Calculate the time difference between detections in days
det_f$tdiff.days <- unlist(tapply(det_f$date, INDEX = det_f$id,
                                  FUN = function(x) c(0, `units<-`(diff(x), "days"))))

## Find the maximum values
max(det_f$tdiff.days)
# [1] 149.5232

## Assign different segments if the time difference is larger than the cutoff
## we choose a cut off of 12 days
## Make sure to retain original ptt info and arrange df by timestamp
det_seg <- det_f %>%
  dplyr::group_by(id) %>%
  dplyr::mutate(id = paste0(id,"_", 1+cumsum(tdiff.days >= 12))) %>%
  dplyr::mutate(shark = as.numeric(str_sub(id, start = 1, end = str_locate(id, "\\_")[,1] - 1))) %>%
  dplyr::arrange( # arrange by timestamp by individual so df can be used for fit_() functions
    shark,
    date
  )

## Let's check if there are sharks that have the tagging data as single segment
print(n = 1000,
      det_seg %>%
        dplyr::group_by(shark)%>%
        dplyr::slice(1:2)
)
## currently two sharks: 179472, 210603 have that issue
## we need to add the tagging event manually later on
## to do so we extract the tagging location of these sharks
# extra_tagging <- det_seg %>%
#   # dplyr::filter(
#   #   shark %in% c(235283, 261743)
#   # ) %>%
#   dplyr::group_by(
#     shark
#   ) %>%
#   dplyr::slice(1) # only get first row, i.e. tagging date

## Based on Logan et al. 2020, segments of <20 days should be removed
## Count number of rows per id
tracklengths <- det_seg %>%
  group_by(id) %>%
  summarise(
    num_locs = n(),
    start_date = min(as.Date(date)),
    end_date = max(as.Date(date)),
  ) %>%
  mutate(tracklength_in_days = as.numeric(end_date - start_date) + 1)

print(tracklengths, n = 999)

#View(tracklengths)

### Filter out short track segments
del_obs <- dplyr::filter(tracklengths, num_locs < 12 | tracklength_in_days < 20) # find suitable parameters
# del_obs <- dplyr::filter(tracklengths, tracklength_in_days < 15) # find suitable parameters

print(n = 999, del_obs); nrow(del_obs); sum(del_obs$num_locs)
# [1] 86 segments
# [1] 675 detections

det_seg_tl <- det_seg %>%
  inner_join(., tracklengths) %>%
  dplyr::filter( # remove tracks shorter than X days and/or less than Y total observations
    !(tracklength_in_days < 20 | num_locs < 12)
    # !(tracklength_in_days < 12)
  ) %>%
  dplyr::select( #remove unnecessary columns
    -start_date,
    -end_date
    #-tracklength_in_days
  )

## Find the distribution of time gaps between detections
#ddet_f_cleaner$tdiff.days[ddet_f_cleaner$tdiff.days < 0] <- 0
hist <- hist(det_seg_tl$tdiff.days
             , breaks = c(seq(from = 0, to = max(ceiling(det_seg_tl$tdiff.days)), by = .5))
             , plot = F)

hist$density <- hist$counts / sum(hist$counts)*100
hist$density

plot(hist, freq = F
     , ylim = c(0, 100)
     , col = "skyblue") # approximately 84 % of detections are less or equal to 12 hours apart. Use this info for track predictions.

## nearly 85% of detections are 12h or less apart, for SSM use 12 hour intervals

# B7: fit MPM using aniMotum ----

### make a map of your raw data before model fitting
world <- ne_countries(scale = 10, returnclass = "sf")

ggplot(data = world) +
  geom_sf() +
  geom_point(data = det_seg_tl, aes(x = lon, y = lat, colour = id), size = 2, shape = 20) +
  coord_sf(xlim = c(min(det_seg_tl$lon - 0.5), max(det_seg_tl$lon + 0.5)),
           ylim = c(min(det_seg_tl$lat - 0.5), max(det_seg_tl$lat + 0.5)), expand = F)

## Define parameters

## The fit_ssm() expects following parameters
# d = a data frame of observations including Argos KF error ellipse info (when present)
# vmax = max travel rate (m/s) passed to sda to identify outlier locations
# ang = angles (deg) of outlier location "spikes"
# distlim =lengths (m) of outlier location "spikes"
# spdf = (logical) turn trip::sda on (default; TRUE) or off
# min.dt = minimum allowable time difference between observations; dt <= min.dt will be ignored by the SSM
# pf = just pre-filter the data, do not fit the SSM (default is FALSE)
# model =  fit either a simple random walk ("rw") or correlated random walk ("crw") as a continuous-time process model
# time.step = options: 1) the regular time interval, in hours, to predict to; 2) a vector of prediction times, possibly not regular, must be specified as a data.frame with id and
# POSIXt dates; 3) NA - turns off prediction and locations are only estimated at observation times.
# scale = scale location data for more efficient optimization. This should rarely be needed (default = FALSE)
# emf = optionally supplied data.frame of error multiplication factors for Argos location quality classes. Default behaviour is to use the factors supplied in foieGras::emf()
# map =  a named list of parameters as factors that are to be fixed during estimation, e.g., list(psi = factor(NA))
# parameters = a list of initial values for all model parameters and unobserved states, default is to let sfilter specify these. Only play with this if you know what you are doing
# fit.to.subset = fit the SSM to the data subset determined by prefilter (default is TRUE)
# control = list of control settings for the outer optimizer (see ssm_control for details)
# inner.control = list of control settings for the inner optimizer (see MakeADFun for additional details)
# verbose = [Deprecated] use ssm_control(verbose = 1) instead, see ssm_control for details
# optim = [Deprecated] use ssm_control(optim = "optim") instead, see ssm_control for details
# optMeth = [Deprecated] use ssm_control(method = "L-BFGS-B") instead, see ssm_control for details
# lpsi [Deprecated] use ssm_control(lower = list(lpsi = -Inf)) instead, see ssm_control for details

## Data
## fit_ssm expects 'd' to be a dataframe or tiblle or sf-tibble (with projection info) with 5,7, or 8 columns
## The data should have  5 columns in the following order: "id", "date", "lc", "lon", "lat".
## Where "date" can be a POSIX object or text string in YYYY-MM-DD HH:MM:SS format
## Argos Kalman Filter (or Kalman Smoother) data should have 8 columns, including the above 5 plus
## "smaj", "smin", "eor" that contain Argos error ellipse variables (in m for "smaj", "smin" and deg for "eor").

## Do you want to calculate fitted observations across entire track?
entire_track <- "No" # change between "Yes" and "No", "No" means you are predicting observations at a given time interval using the segmented tracks

if (entire_track == "Yes"){ ## DO NOT change this
  dpf <- det_f
} else (dpf <- det_seg_tl)

## Remember: According to Lea et al. 2013

## "because each raw position has a different error field according to its Argos location class, we needed to
## decide the most probable location for each point within its error field. We achieved this by using a Bayesian
## state-space model (SSM) that adjusted the filtered tracks by producing regular positions based on the Argos location class,
## mean turning angle, and autocorrelation in speed and direction, producing the most probable track through the error fields"

## "Argos tracks only have locations for when the sharks were at the surface; consequently there is high
## variability in the number of locations in a given area, as a result of the shark’s varied surfacing behaviour
## rather than because of its actual location. This would introduce a bias into the analysis of time spent in different
## areas. To correct this bias, linear interpolation was used to normalise the transmission fre- quency by generating
## points at 12 hour intervals along track gaps of <20 days."

## Speedfilter
## For the prefiltering of the data we need a speedfilter. While other Sphyrna papers use the same speed
## as Vaudo et al. 2017 for mako sharks, this possibly is too high.
## Based on Ryan et al. 2015: https://link.springer.com/article/10.1007/s00227-015-2670-4
## we could use speed as function of FL, i.e. speed = 1xFL [m] * s^-1
## Or we could use modelled cruising speed from Payne et al. 2017, i.e. 2.1 m/s

## Prefilter
#vmax = 2.1 # based on Payne et al. 2016
#ang = c(15,25) # Vaudo et al. 2017, values are internal spikes, i.e. for values you need to 165 = 180 - 15, 155 = 180 - 25
#distlim = c(5000, 8000) # Vaudo et al. 2017
spdf = F
pf = F


## Model
model = "mp" # choose between rw, crw, mp available in fit_ssm()

## Fitted vs. predicted, i.e. choice of normalisation time-step
if(entire_track == "Yes"){
  time.step = NA # NA turns the time step off and estimates locations at observation times only
} else (time.step = 12)


## Optimizers
optim = "optim"
maxit = 2000
verbose = 2

## Prepare directories for mpm plots
dir.create(file.path(saveloc,"plots"), showWarnings = T)

## Get unique species within your df
species_list <- unique(dpf$species); species_list

## Initialize storage for all results
all_species_results <- list()
master_predicted_list <- list()
master_rerouted_list <- list()

## Calculate MPM for each individual across all species
for (sp in species_list) {
  
  cat("Processing Species:", sp, "\n")
  
  # Filter data for this species
  species_data <- dpf %>%
    filter(species == sp)
  
  # Get unique IDs for this species
  shark_ids <- unique(species_data$id)
  cat("Number of individuals or track segments:", length(shark_ids), "\n")
  
  # Initialize lists to store results for this species
  species_results <- list()
  species_predicted_list <- list()
  species_rerouted_list <- list()
  
  # Loop through each individual in this species
  for (shark_id in shark_ids) {
    
    cat("\n  Processing ID:", shark_id, "\n")
    
    # Filter data for this individual
    individual_data <- species_data %>%
      filter(id == shark_id)
    
    # Fit SSM with error handling
    tryCatch({
      mod.mp_pf <- aniMotum::fit_ssm(
        as.data.frame(individual_data),
        spdf = spdf,
        min.dt = 0,
        pf = pf,
        model = model,
        time.step = time.step,
        control = ssm_control(
          optim = optim,
          maxit = maxit,
          verbose = verbose),
        map = list(psi = factor(NA))
      )
      
      # Create timeseries plot
      plot_timeseries <- file.path(saveloc, "plots", paste0(sp, "_", shark_id, "_MPM_predicted_timeseries.png"))
      png(plot_timeseries, width = 25, height = 20, units = "cm", res = 300)
      print(plot(mod.mp_pf, what = "predicted", type = 3, normalise = TRUE)) # add print() to force plot rendering as loop ploitting is acting up
      dev.off()
      
      # create 2d plot
      plot_2d <- file.path(saveloc, "plots", paste0(sp, "_", shark_id, "_MPM_predicted_2D.png"))
      png(plot_2d, width = 25, height = 20, units = "cm", res = 300)
      print(map(mod.mp_pf, what = "p", normalise = TRUE, silent = TRUE)) # add print() to force plot rendering as loop ploitting is acting up
      dev.off()
      
      
      # Re-route locations that fall on land
      mod.mp_pf_rr <- route_path(
        mod.mp_pf,
        what = if(entire_track == "Yes") {"fitted"} else {"predicted"},
        map_scale = 50,
        dist = 5000,
        buffer = 0.5,
        centroids = TRUE,
        append = TRUE
      )
      
      # create 2d plot for rerouted locations
      plot_2d_rr <- file.path(saveloc, "plots", paste0(sp, "_", shark_id, "_MPM_rerouted_2D.png"))
      png(plot_2d_rr, width = 25, height = 20, units = "cm", res = 300)
      print(map(mod.mp_pf_rr, what = "rerouted", normalise = TRUE, silent = TRUE)) # add print() to force plot rendering as loop ploitting is acting up
      dev.off()
      
      # Extract predicted locations
      loc_predicted <- grab(mod.mp_pf, what = "predicted", normalise = T, group = F)
      # Extract rerouted locations
      loc_rerouted <- grab(mod.mp_pf_rr, what = "rerouted", normalise = T, group = F)
      
      # Store results
      species_results[[shark_id]] <- list(
        id = shark_id,
        species = sp,
        model = mod.mp_pf,
        rerouted_model = mod.mp_pf_rr
      )
      
      # Store dataframes for this species
      species_predicted_list[[shark_id]] <- loc_predicted
      species_rerouted_list[[shark_id]] <- loc_rerouted
      
      # Store dataframes for master (all species)
      master_predicted_list[[paste0(sp, "_", shark_id)]] <- loc_predicted
      master_rerouted_list[[paste0(sp, "_", shark_id)]] <- loc_rerouted
      
      cat("    Successfully processed\n")
      
      cat("    ✗ Error:", e$message, "\n")
    })
  }
  
  # Combine all predicted locations for this species
  species_predicted_df <- species_predicted_list %>%
    bind_rows()
  
  # Combine all rerouted locations for this species
  species_rerouted_df <- species_rerouted_list %>%
    bind_rows()
  
  # Save species-level outputs
  if (nrow(species_predicted_df) > 0) {
    saveRDS(species_predicted_df, file.path(saveloc, paste0(sp, "_predicted_all_individuals.rds")))
    write_csv(species_predicted_df, file.path(saveloc, paste0(sp, "_predicted_all_individuals.csv")))
    cat("Saved predicted locations for species:", sp, "\n")
  }
  
  if (nrow(species_rerouted_df) > 0) {
    saveRDS(species_rerouted_df, file.path(saveloc, paste0(sp, "_rerouted_all_individuals.rds")))
    write_csv(species_rerouted_df, file.path(saveloc, paste0(sp, "_rerouted_all_individuals.csv")))
    cat("Saved rerouted locations for species:", sp, "\n")
  }
  
  cat("Total successful individuals:", length(species_results), "\n")
  
  # Store this species results
  all_species_results[[sp]] <- list(
    species = sp,
    results = species_results,
    predicted_data = species_predicted_df,
    rerouted_data = species_rerouted_df
  )
}

# B8: save your master files with all species ----

## Combine master dataframes (all species, all individuals)
### predicted
master_predicted_df <- master_predicted_list %>%
  bind_rows()
### rerouted
master_rerouted_df <- master_rerouted_list %>%
  bind_rows()

## Save master outputs
### predicted
saveRDS(master_predicted_df, file.path(saveloc, "MPM_master_predicted_all_species_all_individuals.rds"))
write_csv(master_predicted_df, file.path(saveloc, "MPM_master_predicted_all_species_all_individuals.csv"))
### rerouted
saveRDS(master_rerouted_df, file.path(saveloc, "MPM_master_rerouted_all_species_all_individuals.rds"))
write_csv(master_rerouted_df, file.path(saveloc, "MPM_master_rerouted_all_species_all_individuals.csv"))

## Save complete results object with models
saveRDS(all_species_results, file.path(saveloc, "All_species_MPM_model_results.rds"));cat("ALL PROCESSING COMPLETE!\n");cat("Species processed:", length(species_list), "\n")

## NOTE: this script is tailored for MPM fitting and extraction of move persistence estimates.
## Validation of model fitting (i.e. ssm-predicted location) is not yet added, as Rscript2 already contains the
## validation steps. If this Rscript is used exclusively (without Rscript2) the validation has to be added.

# END OF CODE ----


 
# # Residual plots are important for validating models, but classical Pearson residuals,
# # for example, are not appropriate for state-space models. Instead, a one-step-ahead prediction
# # residual, provides a useful if computationally demanding alternative.
# # In `aniMotum`, prediction residuals from state-space model fits are calculated using the
# # `osar` function and can be visualized as time-series plots, Q-Q plots, or autocorrelation
# # functions:
# 
# ## calculate & plot residuals
# 
# for (i in 1:length(IDs)){
# 
#   # subset individual tracks
#   res.s <- osar(mod.crw_pf_rr[i,]) ## change accordingly
# 
#   # open plot window
#   png(file = file.path(saveloc,"Validation", paste0("Prediction_", model,"_fitted","_rerouted_residuals_for_validation_", IDs[i],"_fitted_with_", optim, "_",maxit,"iterations_", speedfilter, "_filter.png")), res = 150, height = 20, width = 30, units = "cm")
# 
#   #plot
#   print((plot(res.s, type = "ts") | plot(res.s, type = "qq")) / # when combining plots with patchwork inside a non-interactive context (like saving to a file), you need to explicitly print() the result.
#           (plot(res.s, type = "acf") | plot_spacer()))
#   
#   # save it
#   dev.off()
# }
# 
