# Example running environmental clustering
``example <- envKlusters(aoi = 'AOI.shp', user = 'google_earth_engine_activated_account@gmail.com', gcs_bucket = 'name_of_GoogleCloudStorage_bucket')``

# Basic plot
``plot(terra::mask(example$clusters$layer.1, example$aoi), main = 'Agro-Environmental Clusters (K-means)', col = rainbow(terra::minmax(example$clusters$layer.1)[2]))``
``plot(example$mask, col = 'white', legend = F, add = T)``
``plot(example$pol,add=T)``

# Interactive plot (Leaflet)
``library(leaflet)``
``leaflet() %>%``
``  addTiles() %>% # Add default OpenStreetMap map tiles``
``  addRasterImage(raster::mask(raster::raster(example$clusters$layer.1), sf::st_as_sf(as(example$aoi,'Spatial'))), colors = rainbow(terra::minmax(example$clusters$layer.1)[2]), opacity = 0.7, project = FALSE) %>%``
``  addPolygons(data = sf::st_as_sf(as(example$aoi,'Spatial')), weight = 2, fillColor = 'transparent')``