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
library(tmap)
library(leaflet)

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
Left <- cbind(apr_left@data, apr_left_snap)
coordinates(Left) <- ~x+y
Left <- Left %>% st_as_sf() %>% st_set_crs(epsg3857)

Right <- cbind(apr_right@data, apr_right_snap)
coordinates(Right) <- ~x+y
Right <- Right %>% st_as_sf() %>% st_set_crs(epsg3857)

# Create marker points to help match profile with map
leftprof_lab <- gInterpolate(leftprof, seq(0,105000, by = 15000)) %>% st_as_sf() %>% 
  cbind(seq(0,105000, by = 15000)) %>% rename(marks = seq.0..105000..by...15000.)

rightprof_lab <- gInterpolate(rightprof, seq(0,105000, by = 15000)) %>% st_as_sf() %>% 
  cbind(seq(0,105000, by = 15000)) %>% rename(marks = seq.0..105000..by...15000.)

leftprof_lab$marks <- paste0(leftprof_lab$marks, ' m')
rightprof_lab$marks <- paste0(rightprof_lab$marks, ' m')

#### Temp/Cond Profiles ####

# Plot temperature and conductivity vs distance along profile
profiles = list(Left, Right)
profnames = c('Left', 'Right')

for (i in 1:2){
  cond <- ggplot(profiles[[i]]) + geom_point(aes(x = profdist, y = surfcond_corr, color = "Surface"), size = 0.2) + 
    geom_point(aes(x = profdist, y = deepcond_corr, color = "Deep"), size = 0.2) + 
    labs(title = paste0("April ", profnames[i]), x="", y ="Conductivity (us/cm)", color = "") + 
    scale_color_manual(values = c("Surface" = "grey", "Deep" = "black")) +
    lims(y = c(130, 157)) + guides(colour = guide_legend(override.aes = list(size=5)))
  
  temp <- ggplot(profiles[[i]]) + geom_point(aes(x = profdist, y = surf_WTemp109_Avg, color = "Surface"), size = 0.2) + 
    geom_point(aes(x = profdist, y = deep_WTemp109_Avg, color = "Deep"), size = 0.2) + 
    labs(title = paste0("April ", profnames[i]), x="", y ="Temperature (C)", color = "") +
    scale_color_manual(values = c("Surface" = "coral", "Deep" = "red")) + guides(colour = guide_legend(override.aes = list(size=5)))
  
  ggsave(paste0('202204_', profnames[i],'_cond.png'), plot = cond, width = 8, height = 2.5, units = 'in')
  ggsave(paste0('202204_', profnames[i],'_temp.png'), plot = temp, width = 8, height = 2.5, units = 'in')}

# Plot temperature and conductivity vs time
for (d in 19:21){
  day <- aprdata[day(aprdata$TIMESTAMP) == d,]
  
  cond <- ggplot(day) + geom_point(aes(x = TIMESTAMP, y = surfcond_corr, color = "Surface"), size = 0.2) +
    geom_point(aes(x = TIMESTAMP, y = deepcond_corr, color = "Deep"), size = 0.2) + 
    labs(title = paste0("April ", d, " 2022"), x="", y ="Conductivity (us/cm)", color = "")+ scale_color_manual(values = c("Surface" = "grey", "Deep" = "black"))
  
  temp <- ggplot(day) + geom_point(aes(x = TIMESTAMP, y = surf_WTemp109_Avg, color = "Surface"), size = 0.2) +
    geom_point(aes(x = TIMESTAMP, y = deep_WTemp109_Avg, color = "Deep"), size = 0.2) + 
    labs(title = paste0("April ", d, " 2022"), x="", y ="Temperature (C)", color = "")+ scale_color_manual(values = c("Surface" = "coral", "Deep" = "red"))
  
  ggsave(paste0('202204', d,'_cond.png'), plot = cond, width = 8, height = 2.5, units = 'in')
  ggsave(paste0('202204', d,'_temp.png'), plot = temp, width = 8, height = 2.5, units = 'in')
}

#### Map Radon Samples ####

# Check that all of the radon data survived (start w/ CR105 and end w/ CR130 - the other samples were for deep holes)
c(apr_left_snapped[!is.na(apr_left_snapped$Corrected.Mean..pCi.L.), ]$Tape.Color..., apr_right_snapped[!is.na(apr_right_snapped$Corrected.Mean..pCi.L.), ]$Tape.Color...)

# Subset the radon samples
apr_rn <- filter(apr_sf, !is.na(apr_sf$Corrected.Mean..pCi.L.))
apr_rn$Corrected.Mean..pCi.L. <- round(apr_rn$Corrected.Mean..pCi.L., digits = 2)
# Set up tmap
tmap_mode("view")
radonmap <- tm_shape(apr_rn) + tm_dots(col = "Corrected.Mean..pCi.L.", id = "Corrected.Mean..pCi.L.", title = "Rn (pCi/L)", size = 0.05)
tmap_save(radonmap, filename = "202204_Rn.html")

#### Average Data to 10s ####
# Average data in 10 second chunks for faster plotting
Left$X <- st_coordinates(Left)[,1]
Left$Y <- st_coordinates(Left)[,2]
Right$X <- st_coordinates(Right)[,1]
Right$Y <- st_coordinates(Right)[,2]
Left_data <- st_drop_geometry(Left)
Right_data <- st_drop_geometry(Right)

dat_L <- as.data.table(Left_data[c('TIMESTAMP', 'surfcond_corr', 'deepcond_corr', 'surf_WTemp109_Avg', 'deep_WTemp109_Avg', 'X', 'Y')])
dat_R <- as.data.table(Right_data[c('TIMESTAMP', 'surfcond_corr', 'deepcond_corr', 'surf_WTemp109_Avg', 'deep_WTemp109_Avg', 'X', 'Y')])

Left_10s <- as.data.frame(dat_L[,lapply(.SD, mean),.(TIMESTAMP = round_date(TIMESTAMP, "10 seconds"))])
Right_10s <- as.data.frame(dat_R[,lapply(.SD, mean),.(TIMESTAMP = round_date(TIMESTAMP, "10 seconds"))])

# Make the averaged data spatial
coordinates(Left_10s) <- ~X+Y
coordinates(Right_10s) <- ~X+Y
Left_10s <- Left_10s %>% st_as_sf() %>% st_set_crs(3857)
Right_10s <- Right_10s %>% st_as_sf() %>% st_set_crs(3857)

#### Map the averaged conductivity and temperature data ####
tm_shape(Left_10s) + tm_dots(col="deepcond_corr", border.lwd = 0, style = "cont", palette = "viridis") 
