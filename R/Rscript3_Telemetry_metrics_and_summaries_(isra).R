# TOP OF CODE ----
### ====================================================================================================
### Project:    Large bodied hammerhead complex in the western North Atlantic
### Analysis:   Processing and cleaning satellite telemetry data of fin-mounted SPOT tags for further steps
### Script:     ~2026_IUCN_ISRA_NWA/R/Rscript3_Telemetry_metrics_and_summaries_(isra).R
### Author:     Vital Heim
### Version:    1.0
### ====================================================================================================

### ....................................................................................................
### Content: this script contains the code to quickly extract some movement metrics such as days
###          at liberty, etc. from your satellite telemetry data
### ....................................................................................................

### ....................................................................................................
### [A] Ready environment, load packages ----
### ....................................................................................................

# A1: clear memory ----

rm(list = ls())

# A2: load necessary packages ----

## if for the first time
# install.packages("remotes")
# library(remotes)
# install.packages("tidyverse")
# install.packages("magrittr")
# install.packages("lubridate")
# install.packages("data.table")
# install.packages("sf")
# install.packages("geosphere")
# install.packages("amt")

## load
library(tidyverse)
library(magrittr)
library(lubridate)
library(data.table)
library(sf)
library(geosphere)
library(amt)

# A3: Specify needed functions

## NA

# A4: Specify data and saveloc ----

YOUR_IP <- "NA" # add your IP address, server name or similar if you connect via a shared drive

## Project folder
projloc <- file.path("/",YOUR_IP, "Science","Projects_current", "2026_IUCN_ISRA_NWA")

## Input data
dataloc <- file.path("/",YOUR_IP, "Science","Projects_current", "2026_IUCN_ISRA_NWA","Data_input", "Telemetry")

## Output data
saveloc <- file.path("/", YOUR_IP,"Science","Projects_current", "2026_IUCN_ISRA_NWA","Data_output", "Summaries", "Summaries") # Adjust this
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

### ...................................................................................................
### [A] Data preparation
### ...................................................................................................

# A1 step: import data ----

## Tagging metadata
TAG <- read.csv(file = file.path(dataloc, "Data_Sphyrna_SPOT_tags_metadata.csv"),sep=",",dec=".",header=T,na.strings=c(""," ",NA))
TAG$datetime_deployment_local <- as.POSIXct(TAG$datetime_deployment_local,format="%Y/%m/%d %H:%M:%S",tz="US/Eastern")
attr(TAG$datetime_deployment, "tzone") <- "UTC"
TAG$last_loc <- as.POSIXct(TAG$last_loc,format="%Y/%m/%d %H:%M:%S",tz="UTC")
TAG <- dplyr::select(TAG, ptt_id, group, species, sex, datetime_deployment, deployment_lat, deployment_lon, last_loc, status, pcl, fl, stl)
colnames(TAG) <- c("id", "group", "species", "sex", "date", "lat", "lon", "last_loc", "status", "pcl", "fl", "stl")
TAG$id <- as.character(TAG$id)

## Detections
DET <- readRDS(file.path("/", YOUR_IP,"Science","Projects_current", "2026_IUCN_ISRA_NWA","Data_output","CTCRW", "fitted", paste0("Data_aniMotum_CTCRW_output_","fitted","_non-projected_with_","Argosfilter", "_data.rds")))

# A2: basic housekeeping

## Tagging
TAG_f <- TAG
TAG_f$lat <- round(as.numeric(TAG_f$lat),3)
TAG_f$lon <- round(as.numeric(TAG_f$lon),3)
TAG_f %<>%
  dplyr::filter( # if needed
    species %in% c("S.mokarran", "S.lewini","S.zygaena"),
    group %in% c("Florida Keys", "Jupiter", "Jupiter, FL", "Marquesas", "South Carolina", "Tampa")
  ) %>%
  dplyr::mutate( #combine lat/lon
    latlong = paste0(lat, ", ", lon)
  ) %>%
  dplyr::select(
    id,
    species,
    sex,
    # pcl,
    fl,
    stl,
    group,
    date,
    lat,
    lon,
    latlong,
    last_loc,
    status
  )

## Detections
DET_f <- DET %>%
  dplyr::left_join( # add species and group info to the movement data for later filtering
    TAG_f %>% dplyr::select(id, species, group), by = "id"
  ) %>%
  dplyr::arrange(
    id,
    date
  ) %>%
  dplyr::mutate(
    # ptt_id = as.numeric(str_sub(id, start = 1, end = str_locate(id, "\\_")[,1] - 1))
  ) %>%
  dplyr::select(
    id,
    # ptt_id,
    date,
    lon,
    lat,
    species
  )

### ...................................................................................................
### [C] Calculate track duraion, i.e. days at liberty ----
### ...................................................................................................

# C1.1: get first and last detections per animal using the locations/detections data

# first_obs <- DET_f %>%
#   group_by(id) %>%
#   slice_head(n = 1) %>%
#   ungroup() %>%
#   rename_with(~ paste0("first_", .x), -id)  # Add "first_" prefix to all columns except id
# 
# last_obs <- DET_f %>%
#   group_by(id) %>%
#   slice_tail(n = 1) %>%
#   ungroup() %>%
#   rename_with(~ paste0("last_", .x), -id)   # Add "last_" prefix to all columns except id
# 
# ## join to tagging df
# firstlast_obs <- TAG_f %>%
#   left_join(
#     first_obs, by = "id"
#     ) %>%
#   left_join(
#     last_obs, by = "id"
#   ) %>% 
#   dplyr::arrange(
#     id
#   )
# 
# ## calculate trackduration
# firstlast_obs$first_date <- as.Date(TAG_f$date, format = "%Y-%m-%d")
# firstlast_obs$last_date <- as.Date(TAG_f$last_loc, format = "%Y-%m-%d")
# firstlast_obs$trackduration <- as.numeric(difftime(firstlast_obs$last_date,firstlast_obs$first_date, units = c("days")))+1 # add 1 to include last location day too

# C1.2: calculate number days between tagging and last loc based on metadata

## prepare columns
# TAG_f$start_tag <- as.Date(TAG_f$date, format = "%Y-%m-%d")
# TAG_f$end_tag <- as.Date(TAG_f$last_loc, format = "%Y-%m-%d")

## calculate differences
# TAG_f$trackduration <- as.numeric(difftime(TAG_f$end,TAG_f$start, units = c("days")))+1 # add 1 to include last location day too

## ATTENTION: be aware, that there is a risk that sharks were fished and some last_loc info is from tags on land
## In the next section, also calculate the days at liberty based on the filtered location data.

### ...................................................................................................
### [D] Caluclate track lenghts ----
### ...................................................................................................

# D1: tracklengths in km using the geosphere package ----

DET_sum <- DET_f %>%
  group_by(id) %>%
  mutate(
    distance_m = distGeo(cbind(lon, lat),
                         cbind(lag(lon), lag(lat))),
    start_obs = as.Date(min(date)),
    end_obs = as.Date(max(date)),
    total_days = as.numeric(difftime(as.Date(max(date)), 
                                     as.Date(min(date)), 
                                     units = "days")) + 1
  ) %>%
  summarise(
    start_obs = first(start_obs),
    end_obs = first(end_obs),
    trackduration_days = first(total_days),
    tracklength_km = sum(distance_m, na.rm = TRUE) / 1000
    
  )

# D2: using the amt package ----

# fishlist <- unique(DET$ptt_id)
#
# for (i in fishlist) { #i <- fishlist[9]
#   DETi <- DET[which(DET$ptt_id == i),] #subset to each shark
#   print(paste0(which(fishlist == i), " of ", length(fishlist), "; adding steplength data to ", i))
#   #   track <- amt::mk_track(DETi,
#                          .x = lon,
#                          .y = lat,
#                          .t = date, # was DateTimeUTCmin5
#                          crs = 4326) %>%  #crs = sp::CRS("+init=epsg:4326")) %>%
#     transform_coords(3395)# only needed here as we need projected crs in meters for SL but not for TA
#   # transform_coords(6931)# only needed here as we need projected crs in meters for SL but not for TA
#   stps <- steps(track, lonlat = F)
#   # "b2.3: step  length ----
#   ## (sl;  in  CRSunits), 6931=m ## CHANGE THIS BASED ON YOUR transform_coords() epsg
#   DETi$StepLengthKm <- c(NA, (stps$sl_/1000)) #add NA first since results start from second row
#   # *b2.6: save data to df ----
#   setDT(DET) # convert to data.table without copy
#   setDT(DETi) # convert to data.table without copy
#   # join and update "hammers" by reference, i.e. without copy
#   # https://stackoverflow.com/questions/44930149/replace-a-subset-of-a-data-frame-with-dplyr-join-operations
#   # If you want to return only df_nonai that have a matching hammers (i.e. rows where the key is in both tables), set the nomatch argument of data.table to 0.
#   # nomatch isn't relevant together with :=, ignoring nomatch
#   DET[DETi, on = c("date", "ptt_id"), StepLengthKm := i.StepLengthKm]
#
#   # kinda flat, low depth range. See weekly pdfs.
#   # dives per Y hours < Z  (depends on definition of ‘dive’) = 0
#   # distance covered per day > A? Do histogram of distances, then per
#   # behaviour type to see if distance covered is sig higher when transitioning, e.g.
# }

# D3: calculate the days at liberty using filtered movement data 

### ...................................................................................................
### [E] Combine dataframes and summarise ----
### ...................................................................................................

# E1: make final df ----

final_all <- TAG_f %>%
  dplyr::left_join(
    DET_sum, by = "id"
  ) %>%
  dplyr::arrange(
    id
  )

write.csv(final_all, file.path(saveloc, "Data_summaries_duration_and_tracklengths_all.csv"), row.names = F)

final_survive <- TAG_f %>%
  dplyr::inner_join(
    DET_sum, by = "id"
  ) %>%
  dplyr::arrange(
    id
  )

write.csv(final_survive, file.path(saveloc, "Data_summaries_duration_and_tracklengths_survivors_only.csv"), row.names = F)

# E2: make summaries

## if we calculate summaries across all sharks, not only active tags, we need to take into account pprms
final_all <- final_all %>%
  dplyr::mutate( # replace pprms trackduration (<= 1 days) or lengths in km (NA) with 0
    across(c(trackduration_days, tracklength_km), ~ ifelse(is.na(.x) | .x <= 1, 1, .x))
  )

## tracklengths
avg_TL <- final_all %>%
  group_by(
    species,
    # sex
  ) %>%
  dplyr::summarise(
    n = n(),
    meanTL = mean(tracklength_km),
    sdTL = sd(tracklength_km),
    minTL = min(tracklength_km),
    maxTL = max(tracklength_km)
  );avg_TL

write.csv(avg_TL, file.path(saveloc, "Data_summaries_mean_tracklengths_by_species_incl_pprm.csv"), row.names = F)

## trackdurations
avg_TDur <- final_all %>%
  group_by(
    species
  ) %>%
  dplyr::summarise(
    n = n(),
    meanTDur = mean(trackduration_days),
    sdTDur = sd(trackduration_days),
    minTDur = min(trackduration_days),
    maxTDur = max(trackduration_days)
  ); avg_TDur

write.csv(avg_TDur, file.path(saveloc, "Data_summaries_mean_trackdurations_by_species_incl_pprm.csv"), row.names = F)

## sizes - take 1
avg_size <- final_all %>%
  group_by(
    species,
    sex
  ) %>%
  dplyr::summarise(
    n = n(),
    mean_stl = mean(stl, na.rm = T),
    sd_stl = sd(stl, na.rm = T),
    min_stl = min(stl, na.rm = T),
    max_stl = max(stl, na.rm = T),
    #mean_fl = mean(fl, na.rm = T), # too many FL measurements missing
    #sd_fll = sd(fl, na.rm = T),
    #min_fl = min(fl, na.rm = T),
    #max_fl = max(fl, na.rm = T)
  ); avg_size

write.csv(avg_size, file.path(saveloc, "Data_summaries_mean_STLs_by_species_sex.csv"), row.names = F)

## sizes - take 2
avg_size_spp <- final_all %>%
  group_by(
    species
  ) %>%
  dplyr::summarise(
    n = n(),
    mean_stl = mean(stl, na.rm = T),
    sd_stl = sd(stl, na.rm = T),
    min_stl = min(stl, na.rm = T),
    max_stl = max(stl, na.rm = T),
    #mean_fl = mean(fl, na.rm = T), # too many FL measurements missing
    #sd_fll = sd(fl, na.rm = T),
    #min_fl = min(fl, na.rm = T),
    #max_fl = max(fl, na.rm = T)
  ); avg_size_spp

write.csv(avg_size_spp, file.path(saveloc, "Data_summaries_mean_STLs_by_species.csv"), row.names = F)

## numbers
tag_nrs <- final_all %>%
  group_by(
    species,
    sex
  ) %>%
  dplyr::summarise(
    n = n()
  ); tag_nrs

write.csv(tag_nrs, file.path(saveloc, "Data_summaries_tagged_sharks_numbers_species_sex.csv"), row.names = F)

# END OF CODE ----