library(raster)
library(osmdata)
library(sf)
library(terra)
# Oyo State
ebonyi <- opq(bbox = 'Ebonyi NG', timeout = 25*10) %>%
  add_osm_feature(key = 'name', value = "Ebonyi", value_exact = T) %>%
  add_osm_feature(key = 'admin_level', value = "4", value_exact = T) %>%
  osmdata_sf()
# Osun State
osun <- opq(bbox = 'Osun NG', timeout = 25*10) %>%
  add_osm_feature(key = 'name', value = "Osun", value_exact = T) %>%
  add_osm_feature(key = 'admin_level', value = "4", value_exact = T) %>%
  osmdata_sf()
enugu <- opq(bbox = 'Enugu NG', timeout = 25*10) %>%
  add_osm_feature(key = 'name', value = "Enugu", value_exact = T) %>%
  add_osm_feature(key = 'admin_level', value = "4", value_exact = T) %>%
  osmdata_sf()
## UNION both and create a new AOI
pol <- suppressWarnings(st_union(st_union(oyo$osm_multipolygons[1], osun$osm_multipolygons[1], crs = 'EPSG:4326'), enugu$osm_multipolygons[1], crs = 'EPSG:4326'))
pol <- suppressWarnings(vect(pol))
plot(pol)

# # NGA geom
nga <- raster::getData('GADM' , country="NGA", level=0, path = 'data')
# WroldClim data
clim <- rast(raster::getData('worldclim', var='bio', res=0.5, lon=7, lat=2, path = 'data'))
plot(pol)
plot(terra::mask(terra::crop(clim$bio1_26, pol, snap = 'in'), pol), main="Annual Mean Temperature", add = T)
plot(pol, add = T)

# SRTM data
srtm0 <- terra::rast('data/srtm/srtm_37_10.tif')
srtm1 <- terra::rast('data/srtm/srtm_37_11.tif')
srtm2 <- terra::rast('data/srtm/srtm_38_10.tif')
srtm3 <- terra::rast('data/srtm/srtm_38_11.tif')
srtm4 <- terra::rast('data/srtm/srtm_39_10.tif')
srtm5 <- terra::rast('data/srtm/srtm_39_11.tif')
srtm <- terra::mosaic(srtm0, srtm1, srtm2, srtm3, srtm4, srtm5, fun='mean')
plot(pol)
plot(terra::mask(terra::crop(srtm, pol, snap = 'in'), pol), main="Elevation (m)", add = T)
plot(pol, add = T)

# SoilGrids data
source('soilgrids.R')
# soilSOC.nga <- soilgrids250_data(par = "soc", depth = "0-5", xmin = 2, ymin = 4, xmax = 8, ymax = 10, path = 'data/soil/soc/')
plot(soilSOC.nga)
soilN.nga <- soilgrids250_data(par = "nitrogen", depth = "0-5", xmin = 2, ymin = 4, xmax = 8, ymax = 10, path = 'data/soil/nitro/')
soilPH.nga <- soilgrids250_data(par = "phh2o", depth = "0-5", xmin = 2, ymin = 4, xmax = 8, ymax = 10, path = 'data/soil/ph/')
soilCEC.nga <- soilgrids250_data(par = "cec", depth = "0-5", xmin = 2, ymin = 4, xmax = 8, ymax = 10, path = 'data/soil/cec/')
# soilSand.nga <- soilgrids250_data(par = "sand", depth = "0-5", xmin = 2, ymin = 4, xmax = 8, ymax = 10, path = 'data/soil/sand/')



r.scale <- function(x){(x-raster::cellStats(x, "min"))/(raster::cellStats(x, "max")-raster::cellStats(x, "min"))}
kluster <- raster::stack(terra::crop(r.scale(raster::raster(terra::crop(clim$bio1_26, as(pol, 'Spatial')))), as(pol, 'Spatial')),
                         terra::crop(r.scale(raster::raster(terra::crop(clim$bio12_26, as(pol, 'Spatial')))), as(pol, 'Spatial')),
                         terra::crop(r.scale(raster::raster(terra::crop(clim$bio15_26, as(pol, 'Spatial')))), as(pol, 'Spatial')),
                         terra::crop(r.scale(raster::raster(terra::crop(clim$bio16_26, as(pol, 'Spatial')))), as(pol, 'Spatial')),
                         terra::crop(r.scale(raster::raster(terra::resample(terra::crop(srtm, as(pol, 'Spatial')), clim$bio1_26, 'bilinear'))), as(pol, 'Spatial')),
                         terra::crop(r.scale(raster::raster(terra::resample(terra::crop(soilSOC.nga, as(pol, 'Spatial')), clim$bio1_26, 'near'))), as(pol, 'Spatial')),
                         terra::crop(r.scale(raster::raster(terra::resample(terra::crop(soilN.nga, as(pol, 'Spatial')), clim$bio1_26, 'near'))), as(pol, 'Spatial')),
                         terra::crop(r.scale(raster::raster(terra::resample(terra::crop(soilPH.nga, as(pol, 'Spatial')), clim$bio1_26, 'near'))), as(pol, 'Spatial')),
                         terra::crop(r.scale(raster::raster(terra::resample(terra::crop(soilCEC.nga, as(pol, 'Spatial')), clim$bio1_26, 'near'))), as(pol, 'Spatial')))
                         
                         # terra::crop(r.scale(raster::raster(clim$bio12_26)),as(pol, 'Spatial')),
                         # terra::crop(r.scale(raster::raster(clim$bio15_26)),as(pol, 'Spatial')),
                         # terra::crop(r.scale(raster::raster(clim$bio16_26)),as(pol, 'Spatial')),
                         # terra::crop(r.scale(terra::resample(raster::raster(srtm),raster::raster(clim$bio1_26),'ngb')),as(pol, 'Spatial')),
                         # terra::crop(r.scale(terra::resample(raster::raster(soilSOC.nga),raster::raster(clim$bio1_26),'ngb')),as(pol, 'Spatial')),
                         # terra::crop(r.scale(terra::resample(raster::raster(soilN.nga),raster::raster(clim$bio1_26),'ngb')),as(pol, 'Spatial')),
                         # terra::crop(r.scale(terra::resample(raster::raster(soilPH.nga),raster::raster(clim$bio1_26),'ngb')),as(pol, 'Spatial')),
                         # terra::crop(r.scale(terra::resample(raster::raster(soilCEC.nga),raster::raster(clim$bio1_26),'ngb')),as(pol, 'Spatial')))
# terra::crop(r.scale(terra::resample(raster::raster(soilSand.nga),raster::raster(clim$bio1_26),'ngb')),as(pol, 'Spatial')))
kluster.k <- stackApply(kluster, nlayers(kluster), fun = sum)
nr <- getValues(kluster.k)
set.seed(99)
kmncluster <- kmeans(na.omit(nr), centers = 10, iter.max = 500, nstart = 5, algorithm="Lloyd")
knr <- setValues(kluster, kmncluster$cluster)
mycolor <- c("#fef65b","#ff0000", "#daa520","#0000ff","#0000ff","#00ff00","#cbbeb5",
             "#c3ff5b", "#ff7373", "#00ff00")
plot(mask(knr$layer.1,as(pol, 'Spatial')), main = 'Agro-Environmental Clusters (K-means)', col = mycolor)
plot(pol,add=T)
