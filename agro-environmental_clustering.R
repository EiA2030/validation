envKlusters <- function(nk, aoi){
  pol <- terra::vect(aoi)
  pol <- terra::project(pol,'EPSG:4326')
  xmin <- terra::ext(pol)[1]
  ymin <- terra::ext(pol)[3]
  xmax <- terra::ext(pol)[2]
  ymax <- terra::ext(pol)[4]
  tmp <- tempdir()
  CLIM <- terra::rast(raster::getData('worldclim', var='bio', res=0.5, lon = ((xmin+xmax)/2), lat = ((ymin+ymax)/2), path = tmp))
  unlink(tmp)
  SRTM <- terra::rast(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, res = c(0.0008333333, 0.0008333333))
  for (x in seq(xmin, xmax, by = 5)) {
    for (y in seq(ymin, ymax, by = 5)) {
      x2 = x + 5
      y2 = y + 5
      tmp <- tempfile()
      s <- terra::rast(raster::getData('SRTM', lon = ((x+x2)/2), lat = ((y+y2)/2)), path = tmp)
      SRTM <- terra::mosaic(s,SRTM,fun="mean")
      unlink(tmp)
    }
  }
  devtools::source_url("https://raw.githubusercontent.com/EiA2030/source_data/main/R/soilgrids250_download.R")
  tmp <- tempdir()
  SOC <- soilgrids250_data(par = "soc", depth = "0-5", xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax, path = paste(tmp, "/", sep = ""))
  unlink(tmp)
  tmp <- tempdir()
  N <- soilgrids250_data(par = "nitrogen", depth = "0-5", xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax, path = paste(tmp, "/", sep = ""))
  unlink(tmp)
  tmp <- tempdir()
  PH <- soilgrids250_data(par = "phh2o", depth = "0-5", xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax, path = paste(tmp, "/", sep = ""))
  unlink(tmp)
  tmp <- tempdir()
  CEC <- soilgrids250_data(par = "cec", depth = "0-5", xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax, path = paste(tmp, "/", sep = ""))
  unlink(tmp)
  r.scale <- function(x){(x-raster::cellStats(x, "min"))/(raster::cellStats(x, "max")-raster::cellStats(x, "min"))}
  kluster <- raster::stack(terra::crop(r.scale(raster::raster(CLIM$bio1_26)),as(pol, 'Spatial')),
                           terra::crop(r.scale(raster::raster(CLIM$bio12_26)),as(pol, 'Spatial')),
                           terra::crop(r.scale(raster::raster(CLIM$bio15_26)),as(pol, 'Spatial')),
                           terra::crop(r.scale(raster::raster(CLIM$bio16_26)),as(pol, 'Spatial')),
                           terra::crop(r.scale(terra::resample(raster::raster(SRTM),raster::raster(CLIM$bio1_26),'ngb')),as(pol, 'Spatial')),
                           terra::crop(r.scale(terra::resample(raster::raster(SOC),raster::raster(CLIM$bio1_26),'ngb')),as(pol, 'Spatial')),
                           terra::crop(r.scale(terra::resample(raster::raster(N),raster::raster(CLIM$bio1_26),'ngb')),as(pol, 'Spatial')),
                           terra::crop(r.scale(terra::resample(raster::raster(PH),raster::raster(CLIM$bio1_26),'ngb')),as(pol, 'Spatial')),
                           terra::crop(r.scale(terra::resample(raster::raster(CEC),raster::raster(CLIM$bio1_26),'ngb')),as(pol, 'Spatial')))
  kluster.k <- raster::stackApply(kluster, raster::nlayers(kluster), fun = sum)
  nr <- raster::getValues(kluster.k)
  set.seed(99)
  kmncluster <- kmeans(na.omit(nr), centers = nk, iter.max = 500, nstart = 5, algorithm="Lloyd")
  knr <- raster::setValues(kluster, kmncluster$cluster)
  clusters <- knr
  return(terra::rast(clusters))
}