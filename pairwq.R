library(sp)
library(sf)
library(ggplot2)
library(basemaps)
library(here)

# load area extent shapefile and transform  to 3857 for mapping
areaext <- read_sf("~/columbiariver/public-data/03-processed-data/areabbox.shp") %>% st_transform(3857)

# load the water quality data
febwq <- read.csv('~/columbiariver/field-data/03-processed-data/feb2021_cleandata.csv')
febwq <- rename(febwq, LAT = Lat, LON = Lon)
julwq <- read.csv('~/columbiariver/field-data/03-processed-data/jul2021_cleandata.csv')
julwq <- rename(julwq, LAT = Lat, LON = Lon)
aprwq <- read.csv('~/columbiariver/field-data/03-processed-data/apr2022_cleandata.csv')
aprwq <- rename(aprwq, LAT = Lat, LON = Lon)