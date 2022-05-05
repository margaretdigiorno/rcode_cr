library(here)
library(dplyr)
library(sp)
library(sf)
library(lubridate)
library(ggplot2)

#### Load Data ####
aprwq <- read.csv('~/columbiariver/field-data/03-processed-data/apr2022_cleandata.csv')
aprwq <- rename(aprwq, LAT = Lat, LON = Lon)
aprwq$TIMESTAMP <- as_datetime(aprwq$TIMESTAMP, tz = "America/Los_Angeles")

crpoly <- st_read('~/columbiariver/public-data/02-raw-data/poly_hanfordreach.shp')



#### Make data spatial ####
coordinates(aprwq) <- c('LON', 'LAT')
apr_sf <- st_as_sf(aprwq) %>% st_set_crs(4326) %>% st_transform(3857) 

# Map the data
plot(crpoly$geometry)
plot(apr_sf$geometry, add = TRUE, pch = '.', col = 'red')

# Some of the points appear to be on the river (but not in the polygon) and some points are way outside
# Use a 200 m buffer to see which points are outside the river
buff <- crpoly %>% st_union() %>% st_buffer(400)
onriver <- lengths(st_intersects(apr_sf, buff)) > 0

# Take the points not on the river out of the dataframe
apr_sf <- apr_sf[onriver, ]
offriver_apr <- apr_sf[!onriver, ]

aprdata <- st_drop_geometry(apr_sf)

# Plot temperature and conductivity vs time

for (d in 19:21){
  day <- aprdata[day(aprdata$TIMESTAMP) == d,]
  
  cond <- ggplot(day) + geom_point(aes(x = TIMESTAMP, y = surfcond_corr, color = "Surface")) +
    geom_point(aes(x = TIMESTAMP, y = deepcond_corr, color = "Deep")) + 
    labs(title = paste0("April ", d, " 2022"), x="", y ="Conductivity (us/cm)", color = "")+ scale_color_manual(values = c("Surface" = "grey", "Deep" = "black"))
  
  temp <- ggplot(day) + geom_point(aes(x = TIMESTAMP, y = surf_WTemp109_Avg, color = "Surface")) +
    geom_point(aes(x = TIMESTAMP, y = deep_WTemp109_Avg, color = "Deep")) + 
    labs(title = paste0("April ", d, " 2022"), x="", y ="Temperature (C)", color = "")+ scale_color_manual(values = c("Surface" = "coral", "Deep" = "red"))
  
  ggsave(paste0('202204', d,'_cond.png'), plot = cond, width = 10, height = 3, units = 'in')
  ggsave(paste0('202204', d,'_temp.png'), plot = temp, width = 10, height = 3, units = 'in')
}



