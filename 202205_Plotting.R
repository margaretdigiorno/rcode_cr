library(dplyr)
library(lubridate)
library(ggplot2)
library(gridExtra)
library(grid)
library(sp)
library(sf)
library(tmap)
library(rgeos)
library(basemaps)
library(reshape2)
library(facetscales)
library(ggpubr)
library(scales)
library(latex2exp)
library(patchwork)
library(geosphere)
library(ggrepel)

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
                                                          abs(temp_deepslope) > tempslope & abs(temp_deepdiff) > tempdiff ~ "Deep Anomaly",
                                                          TRUE ~ "Deep Background")))
    
    slopedfs[[i]] <- slopedfs[[i]] %>% mutate(anom_stemp = as.factor(case_when(is.na(temp_surfslope) ~ "NA",
                                                          abs(temp_surfslope) > tempslope & abs(temp_surfdiff) > tempdiff ~ "Surface Anomaly",
                                                          TRUE ~ "Surface Background")))
    
    slopedfs[[i]] <- slopedfs[[i]] %>% mutate(anom_dec = as.factor(case_when(is.na(ec_deepslope) ~ "NA",
                                                          abs(ec_deepslope) > ecslope & abs(ec_deepdiff) > ecdiff ~ "Deep Anomaly",
                                                          TRUE ~ "Deep Background")))
    
    slopedfs[[i]] <- slopedfs[[i]] %>% mutate(anom_sec = as.factor(case_when(is.na(ec_surfslope) ~ "NA",
                                                          abs(ec_surfslope) > ecslope & abs(ec_surfdiff) > ecdiff ~ "Surface Anomaly",
                                                          TRUE ~ "Surface Background")))
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
      top=textGrob(paste0(month,' ', i)),
      ncol = 1
    )
    #ggsave(paste0(month, '_', i,'.png'), plot = plotted, width = 10.5, height = 8, units = 'in')
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
feb_left$side <- "left"
feb_left <- as_Spatial(feb_left)
slot(feb_left, "proj4string") <- epsg3857

feb_right <- filter(feb_sf, side == 2)
feb_right$side <- "right"
feb_right <- as_Spatial(feb_right)
slot(feb_right, "proj4string") <- epsg3857

# Split July
jul_left <- filter(jul_sf, day(jul_sf$TIMESTAMP) == 20 | day(jul_sf$TIMESTAMP) == 22)
jul_left$side <- "left"
jul_left <- as_Spatial(jul_left)
slot(jul_left, "proj4string") <- epsg3857

jul_right <- filter(jul_sf, day(jul_sf$TIMESTAMP) == 21 | day(jul_sf$TIMESTAMP) == 23)
jul_right$side <- "right"
jul_right <- as_Spatial(jul_right)
slot(jul_right, "proj4string") <- epsg3857

# Split April
apr_left <- rbind(apr_sf[day(apr_sf$TIMESTAMP) == 19,], apr_sf[day(apr_sf$TIMESTAMP) == 21 & apr_sf$TIMESTAMP < '2022-04-21 14:27:00',])
apr_left$side <- "left"
apr_left <- as_Spatial(apr_left)
slot(apr_left, "proj4string") <- epsg3857

apr_right <- rbind(apr_sf[day(apr_sf$TIMESTAMP) == 20,], apr_sf[day(apr_sf$TIMESTAMP) == 21 & apr_sf$TIMESTAMP > '2022-04-21 14:27:00',])
apr_right$side <- "right"
apr_right <- as_Spatial(apr_right)
slot(apr_right, "proj4string") <- epsg3857

#### Snap Data to the Profiles ####
alldata <- list(feb_left, feb_right, jul_left, jul_right, apr_left, apr_right)

for (i in 1:length(alldata)){
  # save the original coordinates
  alldata[[i]]$x.orig <- coordinates(alldata[[i]])[,1]
  alldata[[i]]$y.orig <- coordinates(alldata[[i]])[,2]
  
  # snap data to the profiles (if i is even right else left)
  if (i %% 2 == 0){
    alldata[[i]]$profdist <- gProject(rightprof, alldata[[i]])
    alldata[[i]] <- cbind(alldata[[i]]@data, gInterpolate(rightprof, alldata[[i]]$profdist))
  } else {
    alldata[[i]]$profdist <- gProject(leftprof, alldata[[i]])
    alldata[[i]] <- cbind(alldata[[i]]@data, gInterpolate(leftprof, alldata[[i]]$profdist))
  }
  
  alldata[[i]]$profdist <- alldata[[i]]$profdist/1000
  
  # make the data frame sf 
  coordinates(alldata[[i]]) <- ~x+y
  alldata[[i]] <- alldata[[i]] %>% st_as_sf() %>% st_set_crs(epsg3857)
  }

# Add a column of distance between original and snapped wq locations
for (i in 1:6){
  orig <- alldata[[i]][,c('x.orig', 'y.orig')] %>% st_drop_geometry()
  coordinates(orig) <- ~x.orig+y.orig
  orig <- orig %>% st_as_sf() %>% st_set_crs(3857)
  alldata[[i]]$dist <- st_distance(alldata[[i]], orig, by_element = T) %>% as.numeric()
}
#### Map Marker Points ####

# Create marker points to help match profile with map
leftprof_lab <- gInterpolate(leftprof, seq(0,105000, by = 15000)) %>% st_as_sf() %>% 
  cbind(seq(0,105000, by = 15000)) %>% rename(marks = seq.0..105000..by...15000.)

rightprof_lab <- gInterpolate(rightprof, seq(0,105000, by = 15000)) %>% st_as_sf() %>% 
  cbind(seq(0,105000, by = 15000)) %>% rename(marks = seq.0..105000..by...15000.)

leftprof_lab$marks <- leftprof_lab$marks/1000
rightprof_lab$marks <- rightprof_lab$marks/1000
leftprof_lab$marks <- paste0(leftprof_lab$marks, ' km')
rightprof_lab$marks <- paste0(rightprof_lab$marks, ' km')

# Create map of marker points
areaext <- read_sf("~/columbiariver/public-data/03-processed-data/areabbox.shp") %>% st_transform(3857)
set_defaults(map_service = "esri", map_type = "world_imagery", ext = areaext)
basemap_ggplot(areaext)+geom_sf(data = st_crop(crpoly, areaext), fill = 'gray', alpha = 0.2)+ geom_sf(data = leftprof_lab, aes(col = 'Left')) + 
  geom_sf(data = rightprof_lab, aes(color = 'Right'))+
  geom_sf_label(data = leftprof_lab[1:3,], aes(label = marks), nudge_y = 3000, nudge_x = -7000)+
  geom_sf_label(data = leftprof_lab[4:8,], aes(label = marks), nudge_x = 8000)+
  ggtitle("Profile Distances")+labs(color="", x="",y="")+ guides(colour = guide_legend(override.aes = list(size=5)))
#ggsave('markermap.png', width = 6, height = 4, units = 'in')

#### Average deep hole data to single point for each ####

# Load data and convert timestamp class
holedata <- read.csv("Data/aprholes_inwater_cleandata.csv")
holedata$TIMESTAMP <- as_datetime(holedata$TIMESTAMP, tz = "America/Los_Angeles")
coordinates(holedata) <- c('Lon', 'Lat')
holedata <- st_as_sf(holedata) %>% st_set_crs(4326) %>% st_transform(3857)

# Add a segment column so we can group by hole
holedata$s <- cumsum(c(TRUE,diff(holedata$TIMESTAMP)>=300))

# Look at the deep data (can change column to see other attributes)
# We might want to think about selecting times right around the sample time since the sensors were in the water
# longer than we were over each hole (but they look relatively consistent so leaving it for now)
tm_shape(holedata) +
  tm_dots("deeptemp", style = "cont", palette = "viridis") +
  tm_facets(by = "s")

# List the columns that we want to average
tomean <- c('TIMESTAMP', 'surftemp', 'deeptemp', 'deepcond', 'surfcond', 'Corrected.Mean..pCi.L.', 'Quality')

# Group by hole number and summarise using mean/removing NAs
holeavg <- holedata %>% group_by(s) %>% summarise_at(tomean, mean, na.rm = TRUE)

# Find the centroid of all the points for each hole
holeavg$geometry <- st_centroid(holeavg$geometry)

# Sort the deep hole data by whether it's closer to the left or right side
holeavg$side <- st_nearest_feature(holeavg, profs)

holeavg_l <- filter(holeavg, side == 1)
holeavg_l$side <- "left"
holeavg_l <- as_Spatial(holeavg_l)
slot(holeavg_l, "proj4string") <- epsg3857

holeavg_r <- filter(holeavg, side == 2)
holeavg_r$side <- "right"
holeavg_r <- as_Spatial(holeavg_r)
slot(holeavg_r, "proj4string") <- epsg3857

# Save the original coordinates
holeavg_l$x.orig <- coordinates(holeavg_l)[,1]
holeavg_l$y.orig <- coordinates(holeavg_l)[,2]
holeavg_r$x.orig <- coordinates(holeavg_r)[,1]
holeavg_r$y.orig <- coordinates(holeavg_r)[,2]

# Snap hole data to profiles
holeavg_l$profdist <- gProject(leftprof, holeavg_l)
holeavg_l <- cbind(holeavg_l@data, gInterpolate(leftprof, holeavg_l$profdist))

holeavg_r$profdist <- gProject(rightprof, holeavg_r)
holeavg_r <- cbind(holeavg_r@data, gInterpolate(rightprof, holeavg_r$profdist))

holeavg_r$profdist <- holeavg_r$profdist/1000
holeavg_l$profdist <- holeavg_l$profdist/1000

# Make the hole dfs SF
coordinates(holeavg_l) <- ~x+y
coordinates(holeavg_r) <- ~x+y
holeavg_l <- holeavg_l %>% st_as_sf() %>% st_set_crs(epsg3857)
holeavg_r <- holeavg_r %>% st_as_sf() %>% st_set_crs(epsg3857)

#### Prepare data for faceting ####

# All dfs need the same number of columns so that we can use facets
names <- colnames(alldata[[3]])
alldata <- lapply(alldata, select, all_of(names))

# Make an empty dataframe and add all the data to it
combined <- data.frame()

for (i in 1:length(alldata)){
  combined <- rbind(combined, alldata[[i]])
}

# Make month a factor (for nice facet labels)
combined$month <- combined$TIMESTAMP %>% month() %>% factor(levels = c(2,7,4), ordered = FALSE)
levels(combined$month) <- c('2021: February', '2021: July', '2022: April')

# Keep just the columns that we need for plotting
ss <- select(combined, TIMESTAMP, deepcond, surfcond, deeptemp, surftemp, anom_dtemp, anom_stemp, 
             anom_dec, anom_sec, Corrected.Mean..pCi.L., Quality, profdist, month, side, dist, geometry)

# Make separate left/right data frames
ss_l <- ss[ss$side == "left",]
ss_r <- ss[ss$side == "right",]
ptsize <- 0.5

#### Locations for vlines ####
# We want to add vertical lines to the plots so we can see which data was collected on a different day
leftmax <- st_drop_geometry(ss_l) %>% group_by(month, day(TIMESTAMP)) %>% summarise(max = max(TIMESTAMP))
leftmax$side <- "Left"
leftmax$profdist <- NA
for (i in 1:nrow(leftmax)){
  leftmax[i, 'profdist'] <- ss_l[ss_l$TIMESTAMP == leftmax$max[i],]$profdist
}
rightmax <- st_drop_geometry(ss_r) %>% group_by(month, day(TIMESTAMP)) %>% summarise(max = max(TIMESTAMP))
rightmax$side <- "Right"
rightmax$profdist <- 0.00
for (i in 1:nrow(rightmax)){
  rightmax[i, 'profdist'] <- ss_r[ss_r$TIMESTAMP == rightmax$max[i],]$profdist
}

allmax <- rbind(leftmax, rightmax)

# Add the locations of known surface water inflows
# Load the times of known inflows
inflowtimes <- readxl::read_excel('Data/surfaceinflowtimes.xlsx', col_names = c('time', 'label'))
inflowtimes$time <- inflowtimes$time %>% as_datetime() %>% force_tz(tzone = "America/Los_Angeles")

# Find the profile distance of the closest WQ point in time (using df 'ss') for each inflow
library(data.table)
setDT(ss)[,join_date := TIMESTAMP]
setDT(inflowtimes)[,join_date := time]
joined <- ss[inflowtimes, on = .(join_date), roll = "nearest"]

# Load the .csv of wasteway locations
wwloc <- st_read('Data/wasteways.csv') %>% st_set_crs(4326) %>% st_transform(3857) %>% as_Spatial()

# Snap wasteway locations to the left profile and calculate prof dist
wwloc$profdist <- gProject(leftprof, wwloc)
wwloc <- cbind(wwloc@data, gInterpolate(leftprof, wwloc$profdist))
wwloc$profdist <- wwloc$profdist/1000
coordinates(wwloc) <- ~x+y

# Look at the two dataframes (not that many so fine to manually combine)
wwloc$side <- 'left'
wwloc <- select(wwloc@data, name, profdist, side) %>% rename(label = name)
sw_comb <- select(joined, label, profdist, side) %>% rbind(wwloc) %>% arrange(profdist)

# Use observed locations of WWs w/USGS names and remove duplicates
sw_comb[2, 'label'] <- 'WB10 WW'
sw_comb[6, 'label'] <- 'PE 16.4 WW'
sw_comb[8, 'label'] <- 'SW Inflow'
sw_comb[14, 'label'] <- 'Pasco WW'
sw_comb[15, 'label'] <- 'Esquatzel WW'
sw_comb[17, 'label'] <- 'SW Inflow'
sw_comb <- sw_comb[-c(3, 5, 7, 9, 13, 16, 18),]

#### EC/Temp vs. Profile Distance (All Months) ####

# Plot and save spc v profile distance for left and right
cond_l <- ggplot(ss_l[!is.na(ss_l$anom_dec),]) + geom_point(aes(x=profdist, y = surfcond, color = anom_sec), size = ptsize)+ 
  geom_point(aes(x=profdist, y = deepcond, color = anom_dec), size = ptsize)+ ylim(125,155)+
  scale_color_manual(values = c("Deep Background" = "black", "Deep Anomaly" = "red", 
                                "Surface Background" = "gray", "Surface Anomaly" = "coral"))+
  labs(x =  "Distance Along Profile (km)", y = "Specific Conductance (us/cm)", color = "", title = "Left Shore") + 
  facet_grid(month ~.) + guides(colour = guide_legend(override.aes = list(size=3)))
ggsave('cond_l.png', plot = cond_l, width = 14, height = 8, units = 'in')

cond_r <- ggplot(ss_r[!is.na(ss_r$anom_dec),]) + geom_point(aes(x=profdist, y = surfcond, color = anom_sec), size = ptsize)+
  geom_point(aes(x=profdist, y = deepcond, color = anom_dec), size = ptsize) + ylim(125,155) +
  scale_color_manual(values = c("Deep Background" = "black", "Deep Anomaly" = "red", 
                                "Surface Background" = "gray", "Surface Anomaly" = "coral"))+
  labs(x =  "Distance Along Profile (km)", y = "Specific Conductance (us/cm)", color = "", title = "Right Shore") + 
  facet_grid(month ~.)+ guides(colour = guide_legend(override.aes = list(size=3)))
#ggsave('cond_r.png', plot = cond_r, width = 14, height = 8, units = 'in')

# Plot and save temp vs profile distance for left and right
  # Individually specify scales to give each facet the same magnitude
scales_y <- list(
  '2021: February' = scale_y_continuous(limits = c(1, 4)),
  '2021: July' = scale_y_continuous(limits = c(18.75, 21.75)),
  '2022: April' = scale_y_continuous(limits = c(6, 9))
)

temp_l <- ggplot(ss_l[!is.na(ss_l$anom_dtemp),]) + geom_point(aes(x=profdist, y = surftemp, color = anom_stemp), size = ptsize)+
  geom_point(aes(x=profdist, y = deeptemp, color = anom_dtemp), size = ptsize)+
  scale_color_manual(values = c("Deep Background" = "black", "Deep Anomaly" = "red", 
                                "Surface Background" = "gray", "Surface Anomaly" = "coral"))+
  labs(x =  "Distance Along Profile (km)", y = "Temperature (C)", color = "", title = "Left Shore") + 
  facet_grid_sc(month ~., scales = list(y = scales_y))+ guides(colour = guide_legend(override.aes = list(size=3)))
ggsave('temp_l.png', plot = temp_l, width = 14, height = 8, units = 'in')

temp_r <- ggplot(ss_r[!is.na(ss_r$anom_dtemp),]) + geom_point(aes(x=profdist, y = surftemp, color = anom_stemp), size = ptsize)+
  geom_point(aes(x=profdist, y = deeptemp, color = anom_dtemp), size = ptsize)+
  scale_color_manual(values = c("Deep Background" = "black", "Deep Anomaly" = "red", 
                                "Surface Background" = "gray", "Surface Anomaly" = "coral"))+
  labs(x =  "Distance Along Profile (km)", y = "Temperature (C)", color = "", title = "Right Shore") + 
  facet_grid_sc(month ~., scales = list(y = scales_y)) + guides(colour = guide_legend(override.aes = list(size=3)))
ggsave('temp_r.png', plot = temp_r, width = 14, height = 8, units = 'in')

#### Plot SpC/Temp vs Prof Dist with Deep Holes (April) ####

# Make the plots for each side/analyte
cond_lhole <- ggplot(ss_l[(!is.na(ss_l$anom_dec)) & (ss_l$month == "2022: April"),]) + geom_point(aes(x=profdist, y = surfcond, color = anom_sec), size = ptsize)+ 
  geom_point(aes(x=profdist, y = deepcond, color = anom_dec), size = ptsize)+ ylim(125,155)+
  labs(x =  "Distance Along Profile (km)", y = "Specific Conductance (us/cm)", title = "Left Shore Holes: April 2022", color = "", fill = "") + 
  guides(colour = guide_legend(override.aes = list(size=3))) +
  geom_point(data = holeavg_l, aes(x = profdist, y = deepcond, fill = "Hole: Deep"), size = 3, pch = 25, color = "white") + 
  geom_point(data = holeavg_l, aes(x = profdist, y = surfcond, fill = 'Hole: Surface'), size = 3, pch = 25, color = "white")+
  geom_vline(data = allmax[(allmax$month == "2022: April") & (allmax$side == "Left"),], aes(xintercept = profdist))+
  scale_color_manual(values = c("Deep Background" = "black", "Deep Anomaly" = "red", 
                                "Surface Background" = "gray", "Surface Anomaly" = "coral")) + 
  scale_fill_manual(values = c("Hole: Deep" = "black", "Hole: Surface" = "gray"))

cond_rhole <- ggplot(ss_r[(!is.na(ss_r$anom_dec)) & (ss_r$month == "2022: April"),]) + geom_point(aes(x=profdist, y = surfcond, color = anom_sec), size = ptsize)+ 
  geom_point(aes(x=profdist, y = deepcond, color = anom_dec), size = ptsize)+ ylim(125,155)+
  labs(x =  "Distance Along Profile (km)", y = "Specific Conductance (us/cm)", title = "Right Shore Holes: April 2022", color = "", fill = "") + 
  guides(colour = guide_legend(override.aes = list(size=3))) +
  geom_vline(data = allmax[(allmax$month == "2022: April") & (allmax$side == "Right"),], aes(xintercept = profdist))+
  geom_point(data = holeavg_r, aes(x = profdist, y = deepcond, fill = "Hole: Deep"), size = 3, pch = 25, color = "white") + 
  geom_point(data = holeavg_r, aes(x = profdist, y = surfcond, fill = 'Hole: Surface'), size = 3, pch = 25, color = "white")+
  scale_color_manual(values = c("Deep Background" = "black", "Deep Anomaly" = "red", 
                                "Surface Background" = "gray", "Surface Anomaly" = "coral")) + 
  scale_fill_manual(values = c("Hole: Deep" = "black", "Hole: Surface" = "gray"))

temp_lhole <- ggplot(ss_l[(!is.na(ss_l$anom_dec)) & (ss_l$month == "2022: April"),])+ geom_point(aes(x=profdist, y = surftemp, color = anom_stemp), size = ptsize)+
  geom_point(aes(x=profdist, y = deeptemp, color = anom_dtemp), size = ptsize) +
  labs(x =  "Distance Along Profile (km)", y = "Temperature (C)", color = "", fill = "") +
  guides(colour = guide_legend(override.aes = list(size=3))) + ylim(6,8)+
  geom_vline(data = allmax[(allmax$month == "2022: April") & (allmax$side == "Left"),], aes(xintercept = profdist))+
  geom_point(data = holeavg_l, aes(x = profdist, y = deeptemp, fill = "Hole: Deep"), size = 3, pch = 25, color = "white") + 
  geom_point(data = holeavg_l, aes(x = profdist, y = surftemp, fill = 'Hole: Surface'), size = 3, pch = 25, color = "white")+
  scale_color_manual(values = c("Deep Background" = "black", "Deep Anomaly" = "red", 
                                "Surface Background" = "gray", "Surface Anomaly" = "coral")) + 
  scale_fill_manual(values = c("Hole: Deep" = "black", "Hole: Surface" = "gray"))

temp_rhole <- ggplot(ss_r[(!is.na(ss_r$anom_dec)) & (ss_r$month == "2022: April"),])+ geom_point(aes(x=profdist, y = surftemp, color = anom_stemp), size = ptsize)+
  geom_point(aes(x=profdist, y = deeptemp, color = anom_dtemp), size = ptsize) +
  labs(x =  "Distance Along Profile (km)", y = "Temperature (C)", color = "", fill = "") +
  guides(colour = guide_legend(override.aes = list(size=3))) + ylim(6,8)+
  geom_vline(data = allmax[(allmax$month == "2022: April") & (allmax$side == "Right"),], aes(xintercept = profdist))+
  geom_point(data = holeavg_r, aes(x = profdist, y = deeptemp, fill = "Hole: Deep"), size = 3, pch = 25, color = "white") + 
  geom_point(data = holeavg_r, aes(x = profdist, y = surftemp, fill = 'Hole: Surface'), size = 3, pch = 25, color = "white")+
  scale_color_manual(values = c("Deep Background" = "black", "Deep Anomaly" = "red", 
                                "Surface Background" = "gray", "Surface Anomaly" = "coral")) + 
  scale_fill_manual(values = c("Hole: Deep" = "black", "Hole: Surface" = "gray"))

# Arrange and save the plots
ggsave('leftholes.png', plot = ggarrange(cond_lhole, temp_lhole, common.legend = T, ncol = 1, legend = "right"), width = 14, height = 6, units = 'in')
ggsave('rightholes.png', plot = ggarrange(cond_rhole, temp_rhole, common.legend = T, ncol = 1, legend = "right"), width = 14, height = 6, units = 'in')

#### Make plots that include Rn ####

# Bind the hole radons so they can get plotted in the facets (add month for faceting)
  # Left profile
holeavg_l$month <- "2022: April"
holeavg_l$anom_dec <- "Deep Anomaly"
ss_l[setdiff(names(holeavg_l), names(ss_l))] <- NA
holeavg_l[setdiff(names(ss_l), names(holeavg_l))] <- NA
ss_l_rn <- rbind(ss_l, holeavg_l)
ss_l_rn <- ss_l_rn[(!is.na(ss_l_rn$Corrected.Mean..pCi.L.)) & (ss_l_rn$Quality == 1),]


  # Right profile
holeavg_r$month <- "2022: April"
holeavg_r$anom_dec <- "Deep Anomaly"
ss_r[setdiff(names(holeavg_r), names(ss_r))] <- NA
holeavg_r[setdiff(names(ss_r), names(holeavg_r))] <- NA
ss_r_rn <- rbind(ss_r, holeavg_r)
ss_r_rn <- ss_r_rn[(!is.na(ss_r_rn$Corrected.Mean..pCi.L.)) & (ss_r_rn$Quality == 1),]

rnplot <- function(m, side, side_l){
  rn <- ggplot(side[side$month == m,], aes(x = profdist)) + geom_point(aes(y = Corrected.Mean..pCi.L.)) +
          geom_vline(data = allmax[(allmax$month == m) & (allmax$side == side_l),], aes(xintercept = profdist, linetype = "End of Day"))+
          geom_vline(data = sw_comb[sw_comb$side == tolower(side_l)], aes(xintercept = profdist, linetype = "SW Inflow")) +
          ylim(0,35) + xlim(0,117) + scale_linetype_manual(values = c(1,2)) +
          labs(x = "", y = "Rn 222 (pCi/L)", title = paste0(m, ' (', side_l,')'), linetype = "") +
          theme(legend.position="none")
  if (side_l == "Left"){
    return(rn + ggtext::geom_richtext(data = filter(sw_comb[sw_comb$side == tolower(side_l)], grepl('WW', label)), aes(x = profdist, y = +Inf, label = "WW"), angle = 90, hjust = 1, size = 3))
  } else{
    return(rn)
  }
}

rnplot_l <- lapply(unique(ss_l_rn$month), rnplot, ss_l_rn, "Left")
rnplot_r <- lapply(unique(ss_r_rn$month), rnplot, ss_r_rn, "Right")

condplot <- function(m, side, side_l){
  ggplot(side[(!is.na(side$anom_dec)) & (side$month == m) & (side$dist < 500),]) + geom_point(aes(x=profdist, y = surfcond, color = anom_sec), size = ptsize)+ 
    geom_point(aes(x=profdist, y = deepcond, color = anom_dec), size = ptsize)+ ylim(125,155)+
    geom_vline(data = allmax[(allmax$month == m) & (allmax$side == side_l),], aes(xintercept = profdist, linetype = "End of Day"))+
    geom_vline(data = sw_comb[sw_comb$side == tolower(side_l)], aes(xintercept = profdist, linetype = "SW Inflow")) +
    guides(colour = guide_legend(override.aes = list(size=3))) +
    scale_color_manual(values = c("Deep Background" = "black", "Deep Anomaly" = "red", 
                                  "Surface Background" = "gray", "Surface Anomaly" = "coral")) + 
    scale_linetype_manual(values = c(1,2)) +
    labs(x =  "Profile Distance (km)", y = "SpC (us/cm)", color = "", linetype = "") + 
    xlim(0,117) #
  
}

# Make the specific conductance plots
cond_l_plot <- lapply(unique(ss_l$month), condplot, ss_l, "Left")
cond_r_plot <- lapply(unique(ss_r$month), condplot, ss_r, "Right")

# Function for temperature plots
tempplot <- function(m, side, side_l){
  ggplot(side[(!is.na(side$anom_dtemp)) & (side$month == m) & (side$dist < 500),]) +
    geom_point(aes(x=profdist, y = surftemp, color = anom_stemp), size = ptsize) +
    geom_point(aes(x=profdist, y = deeptemp, color = anom_dtemp), size = ptsize) +
    geom_vline(data = allmax[(allmax$month == m) & (allmax$side == side_l),], aes(xintercept = profdist, linetype = "End of Day"))+
    geom_vline(data = sw_comb[sw_comb$side == tolower(side_l)], aes(xintercept = profdist, linetype = "SW Inflow")) +
    labs(x = "Profile Distance (km)", y = "Temp (C)", color = "") +
    guides(colour = guide_legend(override.aes = list(size=3))) +
    scale_color_manual(values = c("Deep Background" = "black", "Deep Anomaly" = "red", 
                                  "Surface Background" = "gray", "Surface Anomaly" = "coral")) + 
    scale_linetype_manual(values = c(1,2)) + xlim(0,117) + scales_y[m] + theme(legend.position="none")
}

# Make the temperature plots
temp_l_plot <- lapply(unique(ss_l$month), tempplot, ss_l, "Left")
temp_r_plot <- lapply(unique(ss_r$month), tempplot, ss_r, "Right")

# Make plots showing radon and specific conductance (or temperature) for all months (to compare)
rn_cond_l <- ggarrange(rnplot_l[[1]], cond_l_plot[[1]], rnplot_l[[2]], cond_l_plot[[2]], rnplot_l[[3]],  cond_l_plot[[3]], ncol = 1, common.legend = T, legend = "right", heights = c(2,3,2,3,2,3))
rn_cond_r <- ggarrange(rnplot_r[[1]], cond_r_plot[[1]], rnplot_r[[2]], cond_r_plot[[2]], rnplot_r[[3]],  cond_r_plot[[3]], ncol = 1, common.legend = T, legend = "right", heights = c(2,3,2,3,2,3))

ggsave('rn_cond_l.png', plot = rn_cond_l, width = 8, height = 12, units = 'in')
ggsave('rn_cond_r.png', plot = rn_cond_r, width = 8, height = 12, units = "in")

rn_temp_l <- ggarrange(rnplot_l[[1]], temp_l_plot[[1]], rnplot_l[[2]], temp_l_plot[[2]], rnplot_l[[3]],  temp_l_plot[[3]], ncol = 1, common.legend = T, legend = "right", heights = c(2,3,2,3,2,3))
rn_temp_r <- ggarrange(rnplot_r[[1]], temp_r_plot[[1]], rnplot_r[[2]], temp_r_plot[[2]], rnplot_r[[3]],  temp_r_plot[[3]], ncol = 1, common.legend = T, legend = "right", heights = c(2,3,2,3,2,3))

ggsave('rn_temp_l.png', plot = rn_temp_l, width = 8, height = 12, units = 'in')
ggsave('rn_temp_r.png', plot = rn_temp_r, width = 8, height = 12, units = 'in')

# Make plots showing radon and specific conductance and temperature (one month at a time)

rn_feb_l <- ggarrange(rnplot_l[[1]], temp_l_plot[[1]], cond_l_plot[[1]], ncol = 1, common.legend = T, legend = "right", heights = c(2,3,3))
rn_jul_l <- ggarrange(rnplot_l[[2]], temp_l_plot[[2]], cond_l_plot[[2]], ncol = 1, common.legend = T, legend = "right", heights = c(2,3,3))
rn_apr_l <- ggarrange(rnplot_l[[3]], temp_l_plot[[3]], cond_l_plot[[3]], ncol = 1, common.legend = T, legend = "right", heights = c(2,3,3))

ggsave('rn_feb_l.png', plot = rn_feb_l, width = 12, height = 6, units = 'in')
ggsave('rn_jul_l.png', plot = rn_jul_l, width = 12, height = 6, units = 'in')
ggsave('rn_apr_l.png', plot = rn_apr_l, width = 12, height = 6, units = 'in')

rn_feb_r <- ggarrange(rnplot_r[[1]], temp_r_plot[[1]], cond_r_plot[[1]], ncol = 1, common.legend = T, legend = "right", heights = c(2,3,3))
rn_jul_r <- ggarrange(rnplot_r[[2]], temp_r_plot[[2]], cond_r_plot[[2]], ncol = 1, common.legend = T, legend = "right", heights = c(2,3,3))
rn_apr_r <- ggarrange(rnplot_r[[3]], temp_r_plot[[3]], cond_r_plot[[3]], ncol = 1, common.legend = T, legend = "right", heights = c(2,3,3))

ggsave('rn_feb_r.png', plot = rn_feb_r, width = 12, height = 6, units = 'in')
ggsave('rn_jul_r.png', plot = rn_jul_r, width = 12, height = 6, units = 'in')
ggsave('rn_apr_r.png', plot = rn_apr_r, width = 12, height = 6, units = 'in')


#### Geophysical Data ####

# Load data
geophys = list.files(path = "Data/floatem_inversions/")
geophys = paste0('Data/floatem_inversions/', geophys)
gp <- lapply(geophys, read.table, header = T)

# Transform data and snap to profiles

for (i in 1:length(gp)){
  # Make spatial, transform, make spatial again
  coordinates(gp[[i]]) <- ~UTMX+UTMY
  gp[[i]] <- gp[[i]] %>% st_as_sf() %>% st_set_crs(32611) %>% st_transform(3857) %>% as_Spatial()
  slot(gp[[i]], "proj4string") <- epsg3857
  
  # save the original coordinates
  gp[[i]]$x.orig <- coordinates(gp[[i]])[,1]
  gp[[i]]$y.orig <- coordinates(gp[[i]])[,2]
  
  # snap data to the profiles (if i is even right else left)
  if (i %% 2 == 0){
    gp[[i]]$profdist <- gProject(rightprof, gp[[i]])
    gp[[i]] <- cbind(gp[[i]]@data, gInterpolate(rightprof, gp[[i]]$profdist))
  } else {
    gp[[i]]$profdist <- gProject(leftprof, gp[[i]])
    gp[[i]] <- cbind(gp[[i]]@data, gInterpolate(leftprof, gp[[i]]$profdist))
  }
  
  gp[[i]]$profdist <- gp[[i]]$profdist/1000
  
  gp[[i]] <- gp[[i]] %>% rowwise(RECORD) %>%  mutate(mean25 = mean(c(RHO_I_2, RHO_I_3, RHO_I_4, RHO_I_5)))
  # make the data frame sf 
  coordinates(gp[[i]]) <- ~x+y
  gp[[i]] <- gp[[i]] %>% st_as_sf() %>% st_set_crs(epsg3857)
}

# Add a column of distance between original and snapped geophys locations
for (i in 1:4){
  orig <- gp[[i]][,c('x.orig', 'y.orig')] %>% st_drop_geometry()
  coordinates(orig) <- ~x.orig+y.orig
  orig <- orig %>% st_as_sf() %>% st_set_crs(3857)
  gp[[i]]$dist <- st_distance(gp[[i]], orig, by_element = T) %>% as.numeric()
}

#### Line Plots Comparing Bed/5 m Resistivity ####
lineplots <- function(df){
  ggplot(df, aes(x = profdist)) + geom_line(aes(y = RHO_I_2, color = "Layer 2"), lwd = 1) + geom_line(aes(y = mean25, color = "First ~ 5m of Sediment"), lwd =1) +
    labs(x = "Profile Distance (km)", y = TeX('$\\log_{10}(Resistivity \\; \\Omega \\cdot m)$'), color = "") + 
    scale_y_continuous(trans = 'log10', limits = c(5,300)) + xlim(0,115)+
    annotation_logticks(sides = "l") 

}

lines_gp <- lapply(gp, lineplots)
names <- c('April Left', 'April Right', 'July Left', 'July Right')

for (i in 1:length(lines_gp)){
  plot <- lines_gp[[i]] + ggtitle(names[i])
  ggsave(paste0('reslayers_', gsub(" ", "", names[i]),'.png'), plot = plot, width = 14, height = 3, units = 'in')
}

#### Plot the resistivity data using geom_rect ####

# Set up function to add all of the rectangles (returns the ggplot)
plotres <- function(df){
  data <- df[(df$profdist > 0),]
  data$xmin <- append(0, data$profdist[1:length(data$profdist)-1])
  scalemin <- quantile((1/as.matrix(st_drop_geometry(select(data, RHO_I_1 : RHO_I_20))))*10000, 0.05)
  scalemax <- quantile((1/as.matrix(st_drop_geometry(select(data, RHO_I_1 : RHO_I_20))))*10000, 0.95)
  plot <- filter(data, dist < 500)
  noplot <- filter(data, !(dist < 500))
  res <- ggplot(data)
  for (i in 1:20){
    above <- filter(plot, !!sym(paste0('DEP_BOT_',i)) < DOI_STANDARD)
    below <- filter(plot, !!sym(paste0('DEP_BOT_',i)) > DOI_STANDARD)
    res <- res + geom_rect(data = above, aes(xmin = xmin, xmax = profdist, 
                                             ymin = !!sym(paste0('DEP_TOP_',i)),
                                             ymax = !!sym(paste0('DEP_BOT_',i)),
                                             fill = 1/(!!sym(paste0('RHO_I_',i)))*10000)) +
      
                geom_rect(data = below, aes(xmin = xmin, xmax = profdist, 
                                            ymin = !!sym(paste0('DEP_TOP_',i)),
                                            ymax = !!sym(paste0('DEP_BOT_',i)),
                                            fill = 1/(!!sym(paste0('RHO_I_',i)))*10000), alpha = 0.2)
  }
  res <- res + geom_rect(data = noplot, aes(xmin = xmin, xmax = profdist, ymin = DEP_TOP_1, ymax = DEP_BOT_20), fill = 'gray')+
    scale_y_reverse() + labs(y = "Depth (m)", x = "Profile Distance (km)") +
    scale_fill_viridis_c('EC (us/cm)', option = 'turbo', trans = "log10", limits = c(scalemin, scalemax), oob = squish)+
    xlim(0,117)
  
  return(res)
  
  }

# Apply the function to all of the geophysical data
ecplots <- lapply(gp, plotres)
ec_l <- ecplots[c(3,1)]
ec_r <- ecplots[c(4,2)]

# Combine the resistivity plots with the radon and conductivity data for July/April
for (i in 2:3){
  plotl <- rnplot_l[[i]]/temp_l_plot[[i]]/cond_l_plot[[i]]/ec_l[[i-1]] + plot_layout(guides = "collect", heights = c(1,1,1,2))
  plotr <- rnplot_r[[i]]/temp_r_plot[[i]]/cond_r_plot[[i]]/ec_r[[i-1]] + plot_layout(guides = "collect", heights = c(1,1,1,2))
  ggsave(paste0('wgeophys_l_', i, '.png'), plot = plotl, width = 12, height = 6, units = 'in')
  ggsave(paste0('wgeophys_r_', i, '.png'), plot = plotr, width = 12, height = 6, units = 'in')
}

# Plot February data with July resistivity then with April resistivity
for (i in 1:2){
  plotl <- rnplot_l[[1]]/temp_l_plot[[1]]/cond_l_plot[[1]]/ec_l[[i]] + plot_layout(guides = "collect", heights = c(1,1,1,2))
  plotr <- rnplot_r[[1]]/temp_r_plot[[1]]/cond_r_plot[[1]]/ec_r[[i]] + plot_layout(guides = "collect", heights = c(1,1,1,2))
  ggsave(paste0('feb_l_', i, 'res.png'), plot = plotl, width = 12, height = 6, units = 'in')
  ggsave(paste0('feb_r_', i, 'res.png'), plot = plotr, width = 12, height = 6, units = 'in')
}

# Extract start/end times
alltimes <- rbind(select(st_drop_geometry(holedata), TIMESTAMP), select(ss, TIMESTAMP))
weeks <- alltimes %>% group_by(month(TIMESTAMP)) %>% summarise(max = max(TIMESTAMP), min = min(TIMESTAMP))
write.csv(weeks, file = 'ends_byweek.csv')
days <- alltimes %>% group_by(month(TIMESTAMP), day(TIMESTAMP)) %>% summarise(min = min(TIMESTAMP))
max <- alltimes %>% group_by(month(TIMESTAMP), day(TIMESTAMP)) %>% summarise(max = max(TIMESTAMP))
daysminmax <- cbind(days, max[,'max'])
write.csv(daysminmax, file = 'ends_byday.csv')
