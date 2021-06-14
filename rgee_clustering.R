envKlusters <- function(aoi, user){
  pol <- terra::vect(aoi)
  pol <- terra::project(pol,'EPSG:4326')
  xmin <- terra::ext(pol)[1]
  ymin <- terra::ext(pol)[3]
  xmax <- terra::ext(pol)[2]
  ymax <- terra::ext(pol)[4]
  
  # Initialize Google Earth Engine
  require(rgee)
  ee_Initialize(email = user, gcs = T, drive = T)
  
  # Create AOI (geometry)
  geom <- sf::st_read(aoi)$geometry %>% sf_as_ee()
  clipbbox <- function(image) {
    return(image$clip(geom))
  }
  
  # Load data
  # Climate (WorldClim)
  amt <- ee$Image('WORLDCLIM/V1/BIO')$select('bio01')$clip(geom)
  prc <- ee$Image('WORLDCLIM/V1/BIO')$select('bio12')$clip(geom)
  cvr <- ee$Image('WORLDCLIM/V1/BIO')$select('bio15')$clip(geom)
  wet <- ee$Image('WORLDCLIM/V1/BIO')$select('bio16')$clip(geom)
  # SRTM
  srtm <- ee$Image('CGIAR/SRTM90_V4')$select('elevation')$clip(geom)
  slp <- ee$Terrain$slope(srtm)
  # Soil data
  SOC <- ee$Image("projects/soilgrids-isric/ocd_mean")$select('ocd_0-5cm_mean')$clip(geom)
  N <- ee$Image("projects/soilgrids-isric/nitrogen_mean")$select('nitrogen_0-5cm_mean')$clip(geom)
  PH <- ee$Image("projects/soilgrids-isric/phh2o_mean")$select('phh2o_0-5cm_mean')$clip(geom)
  CEC <- ee$Image("projects/soilgrids-isric/cec_mean")$select('cec_0-5cm_mean')$clip(geom)
  sand <- ee$Image("projects/soilgrids-isric/sand_mean")$select('sand_0-5cm_mean')$clip(geom)
  clay <- ee$Image("projects/soilgrids-isric/clay_mean")$select('clay_0-5cm_mean')$clip(geom)
  # Vegetation (EVI, NDVI, NPP, GPP)
  veg.means <- function(col, b, sdate = '2020-01-01', edate = '2021-01-01', bbox){
    img <- ee$ImageCollection(col)$
      filterDate(sdate, edate)$
      filterBounds(bbox)$
      select(b)$
      map(clipbbox)$
      mean()
    img <- img$expression('BAND*0.0001', list(BAND = img$select(b)))
    return(img)
  }
  evi <- veg.means(col = "MODIS/006/MOD13A2", b = "EVI", bbox = geom)
  ndvi <- veg.means(col = "MODIS/006/MOD13A2", b = "NDVI", bbox = geom)
  gpp <- veg.means(col = "MODIS/006/MOD17A2H", b = "Gpp", bbox = geom)
  npp <- veg.means(col = "MODIS/006/MOD17A3HGF", b = "Npp", bbox = geom)
 
  # Rescale function
  rescale = function(img) {
    mm <- img$reduceRegion(reducer = ee$Reducer$minMax(), bestEffort = TRUE, geometry = geom)
    mM <- mm$rename(mm$keys(), list('max', 'min'))
    return(mM)
  }
  
  # Raster clusters with K-means
  stack <- ee$ImageCollection$fromImages(list(
    amt$unitScale(ee$Number(rescale(amt$reproject(crs = 'EPSG:4326', scale = 1000))$values()$get(1)),
                            rescale(amt$reproject(crs = 'EPSG:4326', scale = 1000))$values()$get(0))$toFloat()$rename('amt'),
    prc$unitScale(ee$Number(rescale(prc$reproject(crs = 'EPSG:4326', scale = 1000))$values()$get(1)),
                            rescale(prc$reproject(crs = 'EPSG:4326', scale = 1000))$values()$get(0))$toFloat()$rename('prc'),
    cvr$unitScale(ee$Number(rescale(cvr$reproject(crs = 'EPSG:4326', scale = 1000))$values()$get(1)),
                            rescale(cvr$reproject(crs = 'EPSG:4326', scale = 1000))$values()$get(0))$toFloat()$rename('cvr'),
    wet$unitScale(ee$Number(rescale(wet$reproject(crs = 'EPSG:4326', scale = 1000))$values()$get(1)),
                            rescale(wet$reproject(crs = 'EPSG:4326', scale = 1000))$values()$get(0))$toFloat()$rename('wet'),
    srtm$unitScale(ee$Number(rescale(srtm$reproject(crs = 'EPSG:4326', scale = 1000))$values()$get(1)),
                             rescale(srtm$reproject(crs = 'EPSG:4326', scale = 1000))$values()$get(0))$toFloat()$rename('srtm'),
    slp$unitScale(ee$Number(rescale(slp$reproject(crs = 'EPSG:4326', scale = 1000))$values()$get(1)),
                            rescale(slp$reproject(crs = 'EPSG:4326', scale = 1000))$values()$get(0))$toFloat()$rename('slp'),
    SOC$unitScale(ee$Number(rescale(SOC$reproject(crs = 'EPSG:4326', scale = 1000))$values()$get(1)),
                            rescale(SOC$reproject(crs = 'EPSG:4326', scale = 1000))$values()$get(0))$toFloat()$rename('SOC'),
    N$unitScale(ee$Number(rescale(N$reproject(crs = 'EPSG:4326', scale = 1000))$values()$get(1)),
                          rescale(N$reproject(crs = 'EPSG:4326', scale = 1000))$values()$get(0))$toFloat()$rename('N'),
    PH$unitScale(ee$Number(rescale(PH$reproject(crs = 'EPSG:4326', scale = 1000))$values()$get(1)),
                           rescale(PH$reproject(crs = 'EPSG:4326', scale = 1000))$values()$get(0))$toFloat()$rename('PH'),
    CEC$unitScale(ee$Number(rescale(CEC$reproject(crs = 'EPSG:4326', scale = 1000))$values()$get(1)),
                            rescale(CEC$reproject(crs = 'EPSG:4326', scale = 1000))$values()$get(0))$toFloat()$rename('CEC'),
    sand$unitScale(ee$Number(rescale(sand$reproject(crs = 'EPSG:4326', scale = 1000))$values()$get(1)),
                             rescale(sand$reproject(crs = 'EPSG:4326', scale = 1000))$values()$get(0))$toFloat()$rename('sand'),
    clay$unitScale(ee$Number(rescale(clay$reproject(crs = 'EPSG:4326', scale = 1000))$values()$get(1)),
                             rescale(clay$reproject(crs = 'EPSG:4326', scale = 1000))$values()$get(0))$toFloat()$rename('clay'),
    evi$unitScale(ee$Number(rescale(evi$reproject(crs = 'EPSG:4326', scale = 1000))$values()$get(1)),
                            rescale(evi$reproject(crs = 'EPSG:4326', scale = 1000))$values()$get(0))$toFloat()$rename('evi'),
    ndvi$unitScale(ee$Number(rescale(ndvi$reproject(crs = 'EPSG:4326', scale = 1000))$values()$get(1)),
                             rescale(ndvi$reproject(crs = 'EPSG:4326', scale = 1000))$values()$get(0))$toFloat()$rename('ndvi'),
    gpp$unitScale(ee$Number(rescale(gpp$reproject(crs = 'EPSG:4326', scale = 1000))$values()$get(1)),
                            rescale(gpp$reproject(crs = 'EPSG:4326', scale = 1000))$values()$get(0))$toFloat()$rename('gpp'),
    npp$unitScale(ee$Number(rescale(npp$reproject(crs = 'EPSG:4326', scale = 1000))$values()$get(1)),
                            rescale(npp$reproject(crs = 'EPSG:4326', scale = 1000))$values()$get(0))$toFloat()$rename('npp')))

  # Multi-band image with bands as variables
  arrayImage <- stack$toBands()$select(stack$toBands()$bandNames())$toArray()
  # Calculate covariance for PCA (https://developers.google.com/earth-engine/guides/arrays_eigen_analysis)
  covar <- arrayImage$reduceRegion(reducer = ee$Reducer$covariance(), bestEffort = TRUE, geometry = geom)
  covarArray <- ee$Array(covar$get('array'))
  eigens <- covarArray$eigen()
  eigenVectors <- eigens$slice(1, 1)
  # Dictionary with weights of each variable
  principalComponents <- ee$Image(eigenVectors)$matrixMultiply(arrayImage$toArray(1))
  pcImage <- principalComponents$
    arrayProject(list(0))$
    arrayFlatten(list(list('pc_amt', 'pc_prc', 'pc_cvr', 'pc_wet',
                           'pc_srtm', 'pc_slp',
                           'pc_SOC', 'pc_N', 'pc_PH', 'pc_CEC', 'pc_sand', 'pc_clay',
                           'pc_evi', 'pc_ndvi', 'pc_gpp', 'pc_npp')))
  weights <- pcImage$reduceRegion(reducer = ee$Reducer$mean(), bestEffort = TRUE, geometry = geom)
  # Rename dict keys
  oldKeys = weights$keys()
  newKeys = oldKeys$map(
    ee_utils_pyfunc(
      function(x) ee$String(x)$replace('pc_','')
    )
  )
  newWeights <- weights$rename(from = oldKeys, to = newKeys)
  # Apply weight to each variable
  times <- function(img) {
    img <- img$multiply(ee$Number(newWeights$get(img$bandNames()$get(0))))
    return(img)
  }
  klusterW <- stack$map(times)
  klusterKW <- klusterW$toBands()
  
  # Start clustering Weighting
  # Training dataset.
  training <- klusterKW$sample(region = geom, scale = 1000, numPixels = 5000)
  # Instantiate the clusterer and train it.
  clusterer <- ee$Clusterer$
    wekaXMeans(minClusters = 1, maxClusters = 100, maxIterations = 500, seed = 99)$
    train(training)
  # Cluster the input using the trained clusterer.
  result <- klusterKW$cluster(clusterer)
  
  # Export results
  task_img <- ee_image_to_gcs(result, bucket = 'iita_transform_bucket', fileFormat = 'GEO_TIFF', region = geom, crs = 'EPSG:4326', scale = 1000, fileNamePrefix = paste('test_deleteme_', sep = ''))
  task_img$start()
  ee_monitoring(task_img)
  tmp <- tempdir()
  ras <- ee_gcs_to_local(task = task_img, dsn = paste(tmp, 'cluster', sep = '/'), public = TRUE, overwrite = TRUE)
  unlink(tmp)
  ras <- terra::rast(ras)
  return(ras)
}
