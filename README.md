# Agro-Environmental clustering tool
This tools allows users and use cases to define the initial set of agro-environmental clusters in the targeting areas (or areas of study). The tool works as follows:

1. Define the target area: The user should provide the tool with a vector (shapefile or other) of the geographic region of interest.
2. Data gathering: Automatically, the tool collects all necessary agronomic and environmental data (soil, climate, vegetation) to run the clustering analysis.
3. Data aggregation: All data is re-scaled and aggregated into a single dataset as the sum of all variables.
4. PCA: A principal component analysis (PCA) is performed to identify variables explaining most of the spatial variation. Variables explaining more than 0.1 of the variation (in both directions) are kept all other variables are dropped.
5. Weighting: Values of each variable are determined based on the PCA weighting. A new aggregated dataset is created based on it.
6. Cluster optimisation: The tool searches for the optimal number of clusters based on the 'Within-Sum-Squares' variation (elbow method).
7. K-means clustering: The final map with the agro-environmental clusters is produced.

The tool finally returns a list with the AOI, the masked areas in the analisys and the final clusters (map)

## Example running environmental clustering
Here is an example of how the tool can be run. *The tools is under development and not optimized. Large sampling areas result in failing. please run only for Administrative Level 4 or lower scales*

``example <- envKlusters(aoi = 'AOI.shp', user = 'google_earth_engine_activated_account@gmail.com', gcs_bucket = 'name_of_GoogleCloudStorage_bucket')``

Simple plot

``plot(terra::mask(example$clusters$layer.1, example$aoi), main = 'Agro-Environmental Clusters (K-means)', col = rainbow(terra::minmax(example$clusters$layer.1)[2]))
plot(example$mask, col = 'white', legend = F, add = T)
plot(example$pol,add=T)``

Interactive plot (Leaflet)

``library(leaflet)
leaflet() %>%
  addTiles() %>% # Add default OpenStreetMap map tiles
  addRasterImage(raster::mask(raster::raster(example$clusters$layer.1), sf::st_as_sf(as(example$aoi,'Spatial'))), colors = rainbow(terra::minmax(example$clusters$layer.1)[2]), opacity = 0.7, project = FALSE) %>%
  addPolygons(data = sf::st_as_sf(as(example$aoi,'Spatial')), weight = 2, fillColor = 'transparent')``