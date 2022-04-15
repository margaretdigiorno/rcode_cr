library(sp)
library(sf)
library(ggplot2)
library(basemaps)
library(here)


leftbank <- read.table(paste0(here('Data'), '/LEFTBANK_MINFILTER_pro_I02_MOD_inv.xyz'), header = TRUE)
rightbank <- read.table(paste0(here('Data'), '/RIGHTBANK_MINFILTER_pro_I02_MOD_inv.xyz'), header = TRUE)

coordinates(leftbank) <- c('UTMX', 'UTMY')
coordinates(rightbank) <- c('UTMX', 'UTMY')
lb_sf <- st_as_sf(leftbank)
rb_sf <- st_as_sf(rightbank)
lb_sf <- st_set_crs(lb_sf, 32611)
rb_sf <- st_set_crs(rb_sf, 32611)
lb_3857 <- st_transform(lb_sf, "EPSG:3857")
rb_3857 <- st_transform(rb_sf, "EPSG:3857")
top <- st_transform(top, "EPSG:3857")
bottom <- st_transform(bottom, "EPSG:3857")
set_defaults(map_service = "esri", map_type = "world_imagery")

topmap <- basemap_ggplot(top) + geom_sf(data = st_crop(lb_3857, top), size = 0.05, color = "white") + geom_sf(data = st_crop(rb_3857, top), size = 0.05, color = "white")
ggsave('floatemtop.pdf', plot = topmap, device = "pdf", width = 10.5, height = 8, units = "in")

botmap <- basemap_ggplot(bottom) + geom_sf(data = st_crop(lb_3857, bottom), size = 0.05, color = "white") + geom_sf(data = st_crop(rb_3857, bottom), size = 0.05, color = "white")
ggsave('floatembot.pdf', plot = botmap, device = "pdf", width = 8, height = 10.5, units = "in")
