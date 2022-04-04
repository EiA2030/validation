zoning <- function(shape = NULL, gee.user = NULL, gcs.bucket = NULL){
  # require(terra)
  require(factoextra)
  # require(dplyr)
  require(raster)
  
  # Create a temporal output directory
  path <- tempdir()
  # Download input data to temporal directory
  download.vars(aoi = shape, user = gee.user, gcs.bucket = gcs.bucket, dest.path = path)
  
  # Format input AOI
  pol <- terra::vect(shape)
  pol <- terra::project(pol,'EPSG:4326')
  # Read mask
  mask <- terra::rast(paste(path, "MASK.tif", sep = "/"))
  # Read input data
  inputs <- list.files(path = path)[!grepl("MASK.tif", list.files(path = path))]
  r.scale <- function(x){(x-raster::cellStats(x, "min"))/(raster::cellStats(x, "max")-raster::cellStats(x, "min"))}
  inputs <- lapply(inputs, function(x) r.scale(raster::raster(terra::mask(terra::resample(terra::rast(paste(path,x,sep = "/")), mask, method = 'near'), mask, inverse = TRUE))))
  
  # Create stack
  cube <- raster::stack()
  for (layer in inputs) {cube <- raster::addLayer(cube, layer)}

  set.seed(29)
  klust_samp <- raster::sampleRandom(cube, raster::ncell(cube)/5) # Extract 20% of values in each variable
  
  # PCA
  pca <- stats::prcomp(klust_samp, scale. = TRUE)
  pc.coms <- seq_len(length(pca$sdev)) # Number of PCs
  eigen.vals <- pca$sdev # Eigen Values
  comp.loading <- summary(pca)[[6]][c(seq(2,(length(pc.coms)*3)-1, 3))] # Loadings of each PC
  cumm.loading <- summary(pca)[[6]][c(seq(3,(length(pc.coms)*3), 3))] # Cummulative Loadings
  # Assign weights to PCs with eigenvalue >= 1
  n <- 1
  weight.ass <- c()
  while (n <= length(pca$sdev[pca$sdev>=1])) {
    weight.ass <- c(weight.ass, (comp.loading[[n]]/cumm.loading[[length(pca$sdev[pca$sdev>=1])]])*100)
    n <- n + 1
  }
  
  # Construct stack with the selected PCs
  cube_pca <- r.scale(raster::predict(cube, pca, index = 1:length(pca$sdev[pca$sdev>=1])))
  # Assign weights to the values in each selected PC
  n <- 1
  m <- 2
  c <- raster::stack()
  while (n <= length(pca$sdev[pca$sdev>=1])) {
    c.c <- (((comp.loading[[m]]/comp.loading[[3*length(pca$sdev[pca$sdev>=1])]])*100)*cube_pca[[n]])
    c <- raster::addLayer(c, c.c)
    n <- n + 1
    m <- m + 3
  }
  # Add PCs to create a single layer
  c <- r.scale(sum(c))
  
  # Extract values from the PCA processed output
  nr.w <- raster::getValues(c)
  # Create a WSS elbow method graph to find the optimal number of clusters
  optim.wss <- factoextra::fviz_nbclust(as.data.frame(sample(na.omit(nr.w), size = length(na.omit(nr.w))/5)), kmeans, method = 'wss', k.max = 50)
  #? Which variation is smaller than variation with the previous element?
  nc <- which(dplyr::lag(optim.wss$data$y, n = 1L) - optim.wss$data$y < 0.5)[1] # Suing a slope of 0.5. This is an assumption, and can be changed...
  
  # Perform final K-means clustering with the optimal number of clusters
  set.seed(99)
  nr.w[is.na(nr.w)] <- -9999
  kmncluster.w <- kmeans(nr.w, centers = nc, iter.max = 500, nstart = 2, algorithm="Lloyd")
  knr.w <- raster::setValues(raster::raster(cube[[1]]), kmncluster.w$cluster)
  
  # Write final output
  output <- terra::mask(terra::rast(knr.w), mask, inverse = TRUE)
  terra::writeRaster(output, "output.tif", overwrite = TRUE)
}

# Example
# zoning(shape = "path/to/aoi.shp", gee.user = "validated_gee_user@gmail.com", gcs.bucket = "GCS_bucket")
