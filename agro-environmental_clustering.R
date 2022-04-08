envKlusters <- function(aoi, user, gcs_bucket){
  pol <- terra::vect(aoi)
  pol <- terra::project(pol,'EPSG:4326')
  xmin <- terra::ext(pol)[1]
  ymin <- terra::ext(pol)[3]
  xmax <- terra::ext(pol)[2]
  ymax <- terra::ext(pol)[4]
  tmp <- tempdir()
  CLIM <- terra::rast(raster::getData('worldclim', var='bio', res=0.5, lon = ((xmin+xmax)/2), lat = ((ymin+ymax)/2), path = tmp))
  unlink(tmp)
  tmp <- tempdir()
  SRTM <- download.file(paste('https://portal.opentopography.org/API/globaldem?demtype=SRTMGL1',
                              '&south=', ymin,
                              '&north=', ymax,
                              '&west=', xmin,
                              '&east=', xmax,
                              '&outputFormat=GTiff', sep = ''),
                        paste(tmp, 'srtm.tif', sep = '/'))
  SRTM <- terra::rast(paste(tmp, 'srtm.tif', sep = '/'))
  unlink(tmp)
  devtools::source_url("https://raw.githubusercontent.com/EiA2030/source_data/main/R/soilgrids250_download.R")
  tmp <- tempdir()
  SOC <- soilgrids250_data(par = "soc", depth = "5-15", xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax)
  unlink(tmp)
  tmp <- tempdir()
  N <- soilgrids250_data(par = "nitrogen", depth = "5-15", xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax)
  unlink(tmp)
  tmp <- tempdir()
  PH <- soilgrids250_data(par = "phh2o", depth = "5-15", xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax)
  unlink(tmp)
  tmp <- tempdir()
  CEC <- soilgrids250_data(par = "cec", depth = "5-15", xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax)
  unlink(tmp)
  tmp <- tempdir()
  clay <- soilgrids250_data(par = "clay", depth = "5-15", xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax)
  unlink(tmp)
  tmp <- tempdir()
  sand <- soilgrids250_data(par = "sand", depth = "5-15", xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax)
  unlink(tmp)
  # Load Vegetation (MODIS)
  require(rgee)
  ee_Initialize(email = user, gcs = T)
  bb <- ee$Geometry$BBox(xmin,ymin,xmax,ymax)
  clipbb <- function(image) {
    return(image$clip(bb))
  }
  means <- function(col, b, sdate = '2020-01-01', edate = '2021-01-01', bb){
    img <- ee$ImageCollection(col)$
      filterDate(sdate, edate)$
      filterBounds(bb)$
      select(b)$
      map(clipbb)$
      mean()
    img <- img$expression('BAND*0.0001', list(BAND = img$select(b)))
    task_img <- ee_image_to_gcs(img, bucket = 'iita_transform_bucket', scale = 1000, fileFormat = 'GEO_TIFF', region = bb, fileNamePrefix = paste(b,'_', sep = ''))
    task_img$start()
    ee_monitoring(task_img)
    tmp <- tempdir()
    ras <- ee_gcs_to_local(task = task_img, dsn = paste(tmp, b, sep = '/'), public = TRUE, overwrite = TRUE)
    unlink(tmp)
    ras <- terra::rast(ras)
    return(ras)
  }
  evi <- means(col = "MODIS/006/MOD13A2", b = "EVI", bb = bb)
  ndvi <- means(col = "MODIS/006/MOD13A2", b = "NDVI", bb = bb)
  gpp <- means(col = "MODIS/006/MOD17A2H", b = "Gpp", bb = bb)
  npp <- means(col = "MODIS/006/MOD17A3HGF", b = "Npp", bb = bb)
  # Mask: PAs, Urban, water
  bb <- ee$Geometry$BBox(xmin,ymin,xmax,ymax)
  valU <- ee$List(list(2))
  valW <- ee$List(list(1))
  urban <- ee$Image("JRC/GHSL/P2016/BUILT_LDSMT_GLOBE_V1")$select('built')$rename(list('mask'))$clip(bb)
  water <- ee$Image("JRC/GSW1_3/GlobalSurfaceWater")$select('occurrence')$rename(list('mask'))$clip(bb)
  pas <- ee$FeatureCollection('WCMC/WDPA/current/polygons')
  urban <- urban$updateMask(urban$neq(ee$Image$constant(valU))$reduce(ee$Reducer$anyNonZero()))
  water <- water$updateMask(water$gt(ee$Image$constant(valW))$reduce(ee$Reducer$anyNonZero()))
  pas <- pas$filter(ee$Filter$notNull(list("WDPAID")))$reduceToImage(properties = list("WDPAID"), reducer = ee$Reducer$first())$rename(list('mask'))$clip(bb)
  urban <- urban$divide(urban)$ceil()$reproject(crs = "EPSG:4326", scale = 1000)
  water <- water$divide(water)$ceil()$reproject(crs = "EPSG:4326", scale = 1000)
  pas <- pas$divide(pas)$ceil()$reproject(crs = "EPSG:4326", scale = 1000)
  mask <- ee$ImageCollection$fromImages(list(urban, water, pas))$mosaic()
  task_img <- ee_image_to_gcs(mask, bucket = gcs_bucket, scale = 1000, fileFormat = 'GEO_TIFF', region = bb, fileNamePrefix = 'MASK')
  task_img$start()
  ee_monitoring(task_img)
  tmp <- tempdir()
  mask.gee <- ee_gcs_to_local(task = task_img, dsn = paste(tmp, 'MASK.tif', sep = '/'), public = TRUE, overwrite = TRUE)
  unlink(tmp)
  mask <- terra::rast(mask.gee)
  mask <- terra::crop(terra::resample(mask, CLIM[[1]],'ngb'), as(pol, 'Spatial'))
  # Raster clusters with K-means
  r.scale <- function(x){(x-raster::cellStats(x, "min"))/(raster::cellStats(x, "max")-raster::cellStats(x, "min"))}
  kluster <- raster::stack(r.scale(terra::mask(terra::crop(raster::raster(CLIM[[1]]),as(pol, 'Spatial')), raster::raster(mask), inverse = T)),
                           r.scale(terra::mask(terra::crop(raster::raster(CLIM[[12]]),as(pol, 'Spatial')), raster::raster(mask), inverse = T)),
                           r.scale(terra::mask(terra::crop(raster::raster(CLIM[[15]]),as(pol, 'Spatial')), raster::raster(mask), inverse = T)),
                           r.scale(terra::mask(terra::crop(raster::raster(CLIM[[16]]),as(pol, 'Spatial')), raster::raster(mask), inverse = T)),
                           r.scale(terra::mask(terra::crop(terra::resample(raster::raster(SRTM),raster::raster(CLIM[[1]]),'ngb'),as(pol, 'Spatial')), raster::raster(mask), inverse = T)),
                           r.scale(terra::mask(terra::crop(terra::resample(raster::raster(SOC),raster::raster(CLIM[[1]]),'ngb'),as(pol, 'Spatial')), raster::raster(mask), inverse = T)),
                           r.scale(terra::mask(terra::crop(terra::resample(raster::raster(N),raster::raster(CLIM[[1]]),'ngb'),as(pol, 'Spatial')), raster::raster(mask), inverse = T)),
                           r.scale(terra::mask(terra::crop(terra::resample(raster::raster(PH),raster::raster(CLIM[[1]]),'ngb'),as(pol, 'Spatial')), raster::raster(mask), inverse = T)),
                           r.scale(terra::mask(terra::crop(terra::resample(raster::raster(CEC),raster::raster(CLIM[[1]]),'ngb'),as(pol, 'Spatial')), raster::raster(mask), inverse = T)),
                           r.scale(terra::mask(terra::crop(terra::resample(raster::raster(clay),raster::raster(CLIM[[1]]),'ngb'),as(pol, 'Spatial')), raster::raster(mask), inverse = T)),
                           r.scale(terra::mask(terra::crop(terra::resample(raster::raster(sand),raster::raster(CLIM[[1]]),'ngb'),as(pol, 'Spatial')), raster::raster(mask), inverse = T)),
                           r.scale(terra::mask(terra::crop(terra::resample(raster::raster(evi),raster::raster(CLIM[[1]]),'ngb'),as(pol, 'Spatial')), raster::raster(mask), inverse = T)),
                           r.scale(terra::mask(terra::crop(terra::resample(raster::raster(ndvi),raster::raster(CLIM[[1]]),'ngb'),as(pol, 'Spatial')), raster::raster(mask), inverse = T)),
                           r.scale(terra::mask(terra::crop(terra::resample(raster::raster(gpp),raster::raster(CLIM[[1]]),'ngb'),as(pol, 'Spatial')), raster::raster(mask), inverse = T)),
                           r.scale(terra::mask(terra::crop(terra::resample(raster::raster(npp),raster::raster(CLIM[[1]]),'ngb'),as(pol, 'Spatial')), raster::raster(mask), inverse = T)))
  # Run PCA
  set.seed(25)
  rpc <- RStoolbox::rasterPCA(kluster)
  rpc$model$loadings
  loadings <- rpc$model$loadings
  class(loadings) <- "matrix"
  # Apply weighting for each co-variable
  kluster.w <- raster::stack()
  for (i in 1:raster::nlayers(kluster)) {
    if(loadings[i,1] > 0.1 | loadings[i,1] < -0.1){kluster.w <- raster::stack(kluster.w, loadings[i,1]*kluster[[i]])}
  }
  kluster.kw <- raster::stackApply(kluster.w, raster::nlayers(kluster.w), fun = sum)
  nr.w <- raster::getValues(kluster.kw)
  optim.wss <- factoextra::fviz_nbclust(as.data.frame(sample(na.omit(nr.w), size = length(na.omit(nr.w))/10)), kmeans, method = 'wss')
  #? Which variation is smaller than variation with the previous element?
  nc <- which(dplyr::lag(optim.wss$data$y, n = 1L)-optim.wss$data$y > dplyr::lag(dplyr::lag(optim.wss$data$y, n = 1L)-optim.wss$data$y, n=1L))[1]
  set.seed(99)
  kmncluster.w <- kmeans(na.omit(nr.w), centers = nc, iter.max = 500, nstart = 2, algorithm="Lloyd")
  knr.w <- terra::rast(raster::setValues(kluster.w, kmncluster.w$cluster))
  return(list(aoi = pol, mask = mask, clusters = knr.w))
}