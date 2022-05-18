library(dplyr)
library(lubridate)
library(ggplot2)
library(gridExtra)
library(grid)
library(sp)
library(sf)
library(tmap)

#### Load Data ####
slopefiles = list.files(path = "Data/WQ_Slopes", pattern = "V2")
slopefiles = paste0("Data/WQ_Slopes/", slopefiles)
slopes = lapply(slopefiles, read.csv)

apr = rbind(slopes[[1]], slopes[[2]], slopes[[3]])
feb = rbind(slopes[[4]], slopes[[5]], slopes[[6]], slopes[[7]])
jul = rbind(slopes[[8]], slopes[[9]], slopes[[10]], slopes[[11]])

# Store the slope/difference thresholds
tempslope = 0.5
tempdiff = 0.1
ecslope = 3
ecdiff = 0.6

# Make the dfs into a list to loop through
slopedfs = list(feb, jul, apr)

#### ID Anomalies ####

# Use mutate to create anomaly columns (0 : background, 1: anomaly, NA : no slope
for (i in 1:length(slopedfs)){
    slopedfs[[i]]$TIMESTAMP <- as_datetime(slopedfs[[i]]$TIMESTAMP, tz = "America/Los_Angeles")
    slopedfs[[i]] <- slopedfs[[i]] %>% mutate(anom_dtemp = as.factor(case_when(is.na(temp_deepslope) ~ "NA",
                                                          abs(temp_deepslope) > tempslope & abs(temp_deepdiff) > tempdiff ~ "1",
                                                          TRUE ~ "0")))
    
    slopedfs[[i]] <- slopedfs[[i]] %>% mutate(anom_stemp = as.factor(case_when(is.na(temp_surfslope) ~ "NA",
                                                          abs(temp_surfslope) > tempslope & abs(temp_surfdiff) > tempdiff ~ "1",
                                                          TRUE ~ "0")))
    
    slopedfs[[i]] <- slopedfs[[i]] %>% mutate(anom_dec = as.factor(case_when(is.na(ec_deepslope) ~ "NA",
                                                          abs(ec_deepslope) > ecslope & abs(ec_deepdiff) > ecdiff ~ "1",
                                                          TRUE ~ "0")))
    
    slopedfs[[i]] <- slopedfs[[i]] %>% mutate(anom_sec = as.factor(case_when(is.na(ec_surfslope) ~ "NA",
                                                          abs(ec_surfslope) > ecslope & abs(ec_surfdiff) > ecdiff ~ "1",
                                                          TRUE ~ "0")))
}

# Put data back into individual dataframes
feb = slopedfs[[1]]
jul = slopedfs[[2]]
apr = slopedfs[[3]]

#### Plot anomalies to check the thresholds ####
anoms_timeseries <- function(df, month){
  for (i in unique(day(df$TIMESTAMP))){
    temp = df[day(df$TIMESTAMP) == i,]
    plotted <- grid.arrange(
      (ggplot(temp)+geom_point(aes(x = TIMESTAMP, y = surfcond, color = anom_sec))),
      (ggplot(temp)+geom_point(aes(x = TIMESTAMP, y = deepcond, color = anom_dec))),
      (ggplot(temp)+geom_point(aes(x = TIMESTAMP, y = surftemp, color = anom_stemp))),
      (ggplot(temp)+geom_point(aes(x = TIMESTAMP, y = deeptemp, color = anom_dtemp))),
      top=textGrob(paste0(month, i)),
      ncol = 1
    )
    ggsave(paste0(month, '_', i,'.png'), plot = plotted, width = 10.5, height = 8, units = 'in')
  }
}

anoms_timeseries(feb, 'February')
anoms_timeseries(jul, 'July')
anoms_timeseries(apr, 'April')

#### Average Uncorrected Values for Neil ####
# Calculate average uncorrected background values for FloaTEM processing
mean(jul$surf_Cuncorr_Avg[jul$anom_sec == 0])
mean(jul$deep_Cuncorr_Avg[jul$anom_dec == 0])
mean(apr$surf_Cuncorr_Avg[apr$anom_sec == 0])
mean(apr$deep_Cuncorr_Avg[apr$anom_dec == 0])

#### Make Data Spatial ####
coordinates(feb) <- c('Lon', 'Lat')
feb_sf <- st_as_sf(feb) %>% st_set_crs(4326) %>% st_transform(3857)

coordinates(jul) <- c('Lon', 'Lat')
jul_sf <- st_as_sf(jul) %>% st_set_crs(4326) %>% st_transform(3857)

coordinates(apr) <- c('Lon', 'Lat')
apr_sf <- st_as_sf(apr) %>% st_set_crs(4326) %>% st_transform(3857)

#### Check Locations for Points Off River ####
crpoly <- st_read('~/columbiariver/public-data/02-raw-data/poly_hanfordreach.shp')
ggplot(crpoly)+geom_sf()+geom_sf(data = feb_sf, pch = '.', col = 'red')
ggplot(crpoly)+geom_sf()+geom_sf(data = jul_sf, pch = '.', col = 'red')
ggplot(crpoly)+geom_sf()+geom_sf(data = apr_sf, pch = '.', col = 'red')
# Only April has points where GPS was off

#### Clean April Data ####
# Use a 400 m buffer to see which points are outside the river
buff <- crpoly %>% st_union() %>% st_buffer(400)
onriver <- lengths(st_intersects(apr_sf, buff)) > 0

# Take the points not on the river out of the dataframe
apr_sf <- apr_sf[onriver, ]
offriver_apr <- apr_sf[!onriver, ]

#### Left/Right Profiles ####
# Define mapping CRS
epsg3857 <- CRS(SRS_string = "EPSG:3857")

# Load profile shapefiles
leftprof <- st_read('~/columbiariver/field-data/03-processed-data/cr_leftshore.shp') %>% st_transform(3857) %>% as_Spatial()
slot(leftprof, "proj4string") <- epsg3857
rightprof <- st_read('~/columbiariver/field-data/03-processed-data/cr_rightshore.shp') %>% st_transform(3857) %>% as_Spatial()
slot(rightprof, "proj4string") <- epsg3857
profs <- rbind(st_as_sf(leftprof), st_as_sf(rightprof))

# Split February
# Add column for which profile each point is closer to
feb_sf$side <-  st_nearest_feature(feb_sf, profs)

feb_left <- filter(feb_sf, side == 1)
feb_left <- as_Spatial(feb_left)
slot(feb_left, "proj4string") <- epsg3857

feb_right <- filter(feb_sf, side == 2)
feb_right <- as_Spatial(feb_right)
slot(feb_right, "proj4string") <- epsg3857

# Split July
jul_left <- filter(jul_sf, day == 20 | day == 22)
jul_left <- as_Spatial(jul_left)
slot(jul_left, "proj4string") <- epsg3857

jul_right <- filter(jul_sf, day == 21 | day == 23)
jul_right <- as_Spatial(jul_right)
slot(jul_right, "proj4string") <- epsg3857

# Split April
apr_left <- rbind(apr_sf[day(apr_sf$TIMESTAMP) == 19,], apr_sf[day(apr_sf$TIMESTAMP) == 21 & apr_sf$TIMESTAMP < '2022-04-21 14:27:00',])
apr_left <- as_Spatial(apr_left)
slot(apr_left, "proj4string") <- epsg3857

apr_right <- rbind(apr_sf[day(apr_sf$TIMESTAMP) == 20,], apr_sf[day(apr_sf$TIMESTAMP) == 21 & apr_sf$TIMESTAMP > '2022-04-21 14:27:00',])
apr_right <- as_Spatial(apr_right)
slot(apr_right, "proj4string") <- epsg3857

