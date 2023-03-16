envKlusters <- function(aoi, sample.size = 100){
  pol <- terra::vect(aoi)
  pol <- terra::project(pol,'EPSG:4326')
  xmin <- terra::ext(pol)[1]
  ymin <- terra::ext(pol)[3]
  xmax <- terra::ext(pol)[2]
  ymax <- terra::ext(pol)[4]
  tmp <- dirname(aoi)
  system(paste('python3 ./agro-environmental_clustering.py', xmin, ymin, xmax, ymax, tmp, sep = " "))
  bio01 <- terra::rast(paste0(tmp, "/bio01.tif"))
  bio12 <- terra::rast(paste0(tmp, "/bio12.tif"))
  srtm <- terra::rast(paste0(tmp, "/srtm.tif"))
  soc <- terra::rast(paste0(tmp, "/soc.tif"))
  ph <- terra::rast(paste0(tmp, "/ph.tif"))
  cec <- terra::rast(paste0(tmp, "/cec.tif"))
  clay <- terra::rast(paste0(tmp, "/clay.tif"))
  sand <- terra::rast(paste0(tmp, "/sand.tif"))
  ndvi <- terra::rast(paste0(tmp, "/ndvi.tif"))
  # Raster clusters with K-means
  r.scale <- function(x){(x-raster::cellStats(x, "min"))/(raster::cellStats(x, "max")-raster::cellStats(x, "min"))}
  kluster <- raster::stack(r.scale(raster::raster(bio01)),
                           r.scale(raster::raster(bio12)),
                           r.scale(raster::raster(srtm)),
                           r.scale(raster::raster(soc)),
                           r.scale(raster::raster(ph)),
                           r.scale(raster::raster(cec)),
                           r.scale(raster::raster(clay)),
                           r.scale(raster::raster(sand)),
                           r.scale(raster::raster(ndvi)))
  # Run PCA
  set.seed(25)
  rpc <- RStoolbox::rasterPCA(kluster)
  # rpc$model$loadings
  loadings <- rpc$model$loadings
  class(loadings) <- "matrix"
  # Apply weighting for each co-variable
  kluster.w <- raster::stack()
  for (i in 1:raster::nlayers(kluster)) {
    if(loadings[i,1] > 0.1 | loadings[i,1] < -0.1){kluster.w <- raster::stack(kluster.w, loadings[i,1]*kluster[[i]])}
  }
  kluster.kw <- raster::stackApply(kluster.w, raster::nlayers(kluster.w), fun = sum)
  nr.w <- raster::getValues(kluster.kw)
  # The optimal number of cluster is highly dependent on the sample size. ^Size -> ^Klusters and viceversa
  optim.wss <- factoextra::fviz_nbclust(as.data.frame(sample(na.omit(nr.w), size = length(na.omit(nr.w)) / sample.size)), kmeans, method = 'wss', k.max = 50)
  #? Which variation is smaller than variation with the previous element?
  nc <- which(dplyr::lag(optim.wss$data$y, n = 1L)-optim.wss$data$y > dplyr::lag(dplyr::lag(optim.wss$data$y, n = 1L)-optim.wss$data$y, n=1L))[1]
  set.seed(99)
  kmncluster.w <- kmeans(na.omit(nr.w), centers = nc, iter.max = 500, nstart = 2, algorithm="Lloyd")
  knr.w <- terra::rast(raster::setValues(kluster.w, kmncluster.w$cluster))
  terra::writeRaster(knr.w[[1]], paste0(tmp, '/agro-env_clusters.tif'), overwrite=TRUE)
  # return(list(aoi = pol, clusters = knr.w))
}
