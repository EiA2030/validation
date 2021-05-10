
library(raster)
library(osmdata)
library(sf)
# Oyo State
oyo <- opq(bbox = 'Oyo NG', timeout = 25*10) %>%
  add_osm_feature(key = 'name', value = "Oyo", value_exact = T) %>%
  add_osm_feature(key = 'admin_level', value = "4", value_exact = T) %>%
  osmdata_sf()
## UNION both and create a new AOI
pol <- suppressWarnings(st_union(oyo$osm_multipolygons, crs = 'EPSG:4326'))
oyo <- suppressWarnings(vect(st_as_sf(pol), st_crs(pol)))

# # NGA geom
# nga <- raster::getData('GADM' , country="NGA", level=0, path = '/media/data/work/iita/EiA2030/validation/data')
# WroldClim data
clim <- rast(raster::getData('worldclim', var='bio', res=0.5, lon=7, lat=2, path = '/media/data/work/iita/EiA2030/validation/data'))
# plot(mask(crop(clim$bio1_26,oyo),oyo), main="Annual Mean Temperature")
# plot(oyo, add = T)
# SRTM data
srtm0 <- raster::getData('SRTM', lon=7, lat=2, path = '/media/data/work/iita/EiA2030/validation/data/srtm/')
srtm1 <- raster::getData('SRTM', lon=2, lat=10, path = '/media/data/work/iita/EiA2030/validation/data/srtm/')
srtm2 <- raster::getData('SRTM', lon=6, lat=10, path = '/media/data/work/iita/EiA2030/validation/data/srtm/')
srtm3 <- raster::getData('SRTM', lon=10, lat=10, path = '/media/data/work/iita/EiA2030/validation/data/srtm/')
srtm4 <- raster::getData('SRTM', lon=2, lat=14, path = '/media/data/work/iita/EiA2030/validation/data/srtm/')
srtm5 <- raster::getData('SRTM', lon=6, lat=14, path = '/media/data/work/iita/EiA2030/validation/data/srtm/')
srtm6 <- raster::getData('SRTM', lon=10, lat=14, path = '/media/data/work/iita/EiA2030/validation/data/srtm/')
srtm <- rast(terra::mosaic(srtm0, srtm1, srtm2, srtm3, srtm4, srtm5, srtm6, fun='mean'))
# plot(srtm)
# plot(nga, add = T)
# SoilGrids data
source('soilgrids.R')
soilSOC.nga <- soilgrids250_data(par = "soc", depth = "0-5", xmin = 2, ymin = 4, xmax = 12, ymax = 12, path = '/media/data/work/iita/EiA2030/validation/data/soil/soc/')
plot(soilSOC.nga)
soilN.nga <- soilgrids250_data(par = "nitrogen", depth = "0-5", xmin = 2, ymin = 4, xmax = 12, ymax = 12, path = '/media/data/work/iita/EiA2030/validation/data/soil/nitro/')
soilPH.nga <- soilgrids250_data(par = "phh2o", depth = "0-5", xmin = 2, ymin = 4, xmax = 12, ymax = 12, path = '/media/data/work/iita/EiA2030/validation/data/soil/ph/')
soilCEC.nga <- soilgrids250_data(par = "cec", depth = "0-5", xmin = 2, ymin = 4, xmax = 12, ymax = 12, path = '/media/data/work/iita/EiA2030/validation/data/soil/cec/')
soilSand.nga <- soilgrids250_data(par = "sand", depth = "0-5", xmin = 2, ymin = 4, xmax = 12, ymax = 12, path = '/media/data/work/iita/EiA2030/validation/data/soil/sand/')



r.scale <- function(x){(x-raster::cellStats(x, "min"))/(raster::cellStats(x, "max")-raster::cellStats(x, "min"))}
kluster <- raster::stack(terra::crop(r.scale(raster::raster(clim$bio1_26)),as(oyo, 'Spatial')),
                         terra::crop(r.scale(raster::raster(clim$bio12_26)),as(oyo, 'Spatial')),
                         terra::crop(r.scale(raster::raster(clim$bio15_26)),as(oyo, 'Spatial')),
                         terra::crop(r.scale(raster::raster(clim$bio16_26)),as(oyo, 'Spatial')),
                         terra::crop(r.scale(terra::resample(raster::raster(srtm),raster::raster(clim$bio1_26),'ngb')),as(oyo, 'Spatial')),
                         terra::crop(r.scale(terra::resample(raster::raster(soilSOC.nga),raster::raster(clim$bio1_26),'ngb')),as(oyo, 'Spatial')),
                         terra::crop(r.scale(terra::resample(raster::raster(soilN.nga),raster::raster(clim$bio1_26),'ngb')),as(oyo, 'Spatial')),
                         terra::crop(r.scale(terra::resample(raster::raster(soilPH.nga),raster::raster(clim$bio1_26),'ngb')),as(oyo, 'Spatial')),
                         terra::crop(r.scale(terra::resample(raster::raster(soilCEC.nga),raster::raster(clim$bio1_26),'ngb')),as(oyo, 'Spatial')),
                         terra::crop(r.scale(terra::resample(raster::raster(soilSand.nga),raster::raster(clim$bio1_26),'ngb')),as(oyo, 'Spatial')))
kluster <- stackApply(kluster, nlayers(kluster), fun = sum)
nr <- getValues(kluster)
set.seed(99)
kmncluster <- kmeans(na.omit(nr), centers = 5, iter.max = 500, nstart = 5, algorithm="Lloyd")
knr <- setValues(kluster, kmncluster$cluster)
mycolor <- c("#fef65b","#ff0000", "#daa520","#0000ff","#0000ff","#00ff00","#cbbeb5",
             "#c3ff5b", "#ff7373", "#00ff00", "#808080")
plot(mask(knr,as(oyo, 'Spatial')), main = 'Agro-Environmental Clusters (K-means)', col = mycolor )
plot(oyo,add=T)

