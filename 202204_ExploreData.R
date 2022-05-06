library(here)
library(dplyr)
library(sp)
library(sf)
library(lubridate)
library(ggplot2)
library(xts)
library(data.table)
library(basemaps)
library(rgeos)

#### Load WQ Data ####
aprwq <- read.csv('~/columbiariver/field-data/03-processed-data/apr2022_cleandata.csv')
aprwq <- rename(aprwq, LAT = Lat, LON = Lon)
aprwq$TIMESTAMP <- as_datetime(aprwq$TIMESTAMP, tz = "America/Los_Angeles")

crpoly <- st_read('~/columbiariver/public-data/02-raw-data/poly_hanfordreach.shp')

#### Make data spatial ####
coordinates(aprwq) <- c('LON', 'LAT')
apr_sf <- st_as_sf(aprwq) %>% st_set_crs(4326) %>% st_transform(3857) 


# Some of the points appear to be on the river (but not in the polygon) and some points are way outside
# Use a 400 m buffer to see which points are outside the river
buff <- crpoly %>% st_union() %>% st_buffer(400)
onriver <- lengths(st_intersects(apr_sf, buff)) > 0

# Take the points not on the river out of the dataframe
apr_sf <- apr_sf[onriver, ]
offriver_apr <- apr_sf[!onriver, ]

#### Snap data to profiles ####
# Create crs object
epsg3857 <- CRS(SRS_string = "EPSG:3857")

# Load profile shapefiles
leftprof <- st_read('~/columbiariver/field-data/03-processed-data/cr_leftshore.shp') %>% st_transform(3857) %>% as_Spatial()
slot(leftprof, "proj4string") <- epsg3857
rightprof <- st_read('~/columbiariver/field-data/03-processed-data/cr_rightshore.shp') %>% st_transform(3857) %>% as_Spatial()
slot(rightprof, "proj4string") <- epsg3857


# Split WQ data into left and right bank
apr_left <- rbind(apr_sf[day(apr_sf$TIMESTAMP) == 19,], apr_sf[day(apr_sf$TIMESTAMP) == 21 & apr_sf$TIMESTAMP < '2022-04-21 14:27:00',])
apr_left <- as_Spatial(apr_left)
slot(apr_left, "proj4string") <- epsg3857

apr_right <- rbind(apr_sf[day(apr_sf$TIMESTAMP) == 20,], apr_sf[day(apr_sf$TIMESTAMP) == 21 & apr_sf$TIMESTAMP > '2022-04-21 14:27:00',])
apr_right <- as_Spatial(apr_right)
slot(apr_right, "proj4string") <- epsg3857

# Snap WQ data to the bank profiles
# Create columns of coordinates (need to store the original point location)
apr_left$x.orig <- coordinates(apr_left)[,1]
apr_left$y.orig <- coordinates(apr_left)[,2]
apr_right$x.orig <- coordinates(apr_right)[,1]
apr_right$y.orig <- coordinates(apr_right)[,2]

# Calculate distance of points along the profile
apr_left$profdist <- gProject(leftprof, apr_left)
apr_right$profdist <- gProject(rightprof, apr_right)
  
# Create spatial point dataframes of the snapped points 
apr_left_snap <- gInterpolate(leftprof, apr_left$profdist)  
apr_right_snap <- gInterpolate(rightprof, apr_right$profdist)
  
# Join the data back to the spatial point dataframes of the snapped points and make them into sf objects (for plotting)
apr_left_snapped <- cbind(apr_left@data, apr_left_snap)
coordinates(apr_left_snapped) <- ~x+y
apr_left_snapped <- apr_left_snapped %>% st_as_sf() %>% st_set_crs(epsg3857)

apr_right_snapped <- cbind(apr_right@data, apr_right_snap)
coordinates(apr_right_snapped) <- ~x+y
apr_right_snapped <- apr_right_snapped %>% st_as_sf() %>% st_set_crs(epsg3857)

#### Temp/Cond Profiles ####

# Plot temperature and conductivity vs distance along profile
for (df in list(apr_left_snapped, apr_right_snapped)){
  cond <- ggplot(df) + geom_point(aes(x = TIMESTAMP, y = surfcond_corr, color = "Surface")) +
    geom_point(aes(x = TIMESTAMP, y = deepcond_corr, color = "Deep")) + 
    labs(title = paste0("April ", d, " 2022"), x="", y ="Conductivity (us/cm)", color = "")+ scale_color_manual(values = c("Surface" = "grey", "Deep" = "black"))
  
  temp <- ggplot(df) + geom_point(aes(x = TIMESTAMP, y = surf_WTemp109_Avg, color = "Surface")) +
    geom_point(aes(x = TIMESTAMP, y = deep_WTemp109_Avg, color = "Deep")) + 
    labs(title = paste0("April ", d, " 2022"), x="", y ="Temperature (C)", color = "")+ scale_color_manual(values = c("Surface" = "coral", "Deep" = "red"))
  
  ggsave(paste0('202204', df,'_cond.png'), plot = cond, width = 8, height = 2.5, units = 'in')
  ggsave(paste0('202204', df,'_temp.png'), plot = temp, width = 8, height = 2.5, units = 'in')}

# Plot temperature and conductivity vs time

for (d in 19:21){
  day <- aprdata[day(aprdata$TIMESTAMP) == d,]
  
  cond <- ggplot(day) + geom_point(aes(x = TIMESTAMP, y = surfcond_corr, color = "Surface")) +
    geom_point(aes(x = TIMESTAMP, y = deepcond_corr, color = "Deep")) + 
    labs(title = paste0("April ", d, " 2022"), x="", y ="Conductivity (us/cm)", color = "")+ scale_color_manual(values = c("Surface" = "grey", "Deep" = "black"))
  
  temp <- ggplot(day) + geom_point(aes(x = TIMESTAMP, y = surf_WTemp109_Avg, color = "Surface")) +
    geom_point(aes(x = TIMESTAMP, y = deep_WTemp109_Avg, color = "Deep")) + 
    labs(title = paste0("April ", d, " 2022"), x="", y ="Temperature (C)", color = "")+ scale_color_manual(values = c("Surface" = "coral", "Deep" = "red"))
  
  ggsave(paste0('202204', d,'_cond.png'), plot = cond, width = 8, height = 2.5, units = 'in')
  ggsave(paste0('202204', d,'_temp.png'), plot = temp, width = 8, height = 2.5, units = 'in')
}

#### Map Radon Samples ####

# Check that all of the radon data survived (start w/ CR105 and end w/ CR130 - the other samples were for deep holes)
c(apr_left_snapped[!is.na(apr_left_snapped$Corrected.Mean..pCi.L.), ]$Tape.Color..., apr_right_snapped[!is.na(apr_right_snapped$Corrected.Mean..pCi.L.), ]$Tape.Color...)


#### Average Data to 10s ####
# Average data in 10 second chunks for faster plotting
apr_sf$X <- st_coordinates(apr_sf)[,1]
apr_sf$Y <- st_coordinates(apr_sf)[,2]
aprdata <- st_drop_geometry(apr_sf)

dat <- as.data.table(aprdata[c('TIMESTAMP', 'surfcond_corr', 'deepcond_corr', 'surf_WTemp109_Avg', 'deep_WTemp109_Avg', 'X', 'Y')])
aprdata_10s <- as.data.frame(dat[,lapply(.SD, mean),.(TIMESTAMP = round_date(TIMESTAMP, "10 seconds"))])
