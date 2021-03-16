
# Intro a Grass, video no. 3 ----
# Parte reutilizable del script ----
# Cargar paquetes
library(rgrass7)
library(sp)
use_sp()
library(sf)
library(mapview)

gisdbase <- 'grass-data-test' #Base de datos de GRASS GIS
wd <- getwd() #Directorio de trabajo
wd
loc <- initGRASS(gisBase = "/usr/lib/grass78/",
                 home = wd,
                 gisDbase = paste(wd, gisdbase, sep = '/'),
                 location = 'guayubin',
                 mapset = "PERMANENT",
                 override = TRUE)

# Imprimir fuentes en la region
execGRASS(
  'g.list',
  flags = 't',
  parameters = list(
    type = c('raster', 'vector')
  )
)

# Limpiar archivo de bloqueo del conjunto de mapas de GRASS
unlink_.gislock()




# Fin de la parte reutilizable


# Video no. 4, Definir proyección de la región de GRASS GIS, importar fuente y utilizarla para definir extensión y resolución. Cómo ver la ayuda de las funciones ----
# Muestra la definición de la región
gmeta()

# Definir ruta del DEM
dem <- 'datos-fuente/srtm_dem_cuenca_guayubin.tif'
execGRASS(
  cmd = 'g.proj',
  flags = c('t','c'),
  georef = dem)

# Muestra la definición de la región modificada
gmeta()

#r.in.gdal importa la fuente a GRASS
execGRASS(
  cmd = 'r.in.gdal',
  flags=c('overwrite','quiet'),
  parameters=list(
    input=dem,
    output='dem'
  )
)

# Actualizar la extensión de la región al DEM, sólo por precaución
execGRASS(
  cmd = 'g.region',
  parameters=list(
    raster = 'dem',
    align = 'dem'
  )
)

# Muestra la región de Grass again
gmeta()

# Importar vectorial a la región de Grass
demext <- 'datos-fuente/srtm_dem_cuenca_guayubin.geojson'
execGRASS(
  cmd = 'v.in.ogr',
  flags=c('overwrite','quiet'),
  parameters=list(
    input = demext,
    output = 'dem_extent'
  )
)

# Imprimir lista de mapas ráster y vectoriales dentro en la región/localización activa
execGRASS(
  'g.list',
  flags = 't',
  parameters = list(
    type = c('raster', 'vector')
  )
)

# Ver los addons disponibles en el repositorio oficial de GRASS GIS, incluyendo descripción
execGRASS(
  cmd = 'g.extension',
  flags = 'c'
)

# Consultar la ayuda de una función
parseGRASS("r.in.gdal")

# Consultar la ayuda de una función. Segunda alternativa
system('r.in.gdal --help')

# Limpiar archivo de bloqueo del conjunto de mapas de GRASS
unlink_.gislock()


# Video no. 5, Explorar datos espaciales básicos entre GRASS y R ----
# Imprimir lista de mapas ráster y vectoriales dentro en la región/localización activa
execGRASS(
  'g.list',
  flags = 't',
  parameters = list(
    type = c('raster', 'vector')
  )
)

# Cargar en R el DEM (mapa ráster)
library(sp)
use_sp()
dem_sp <- readRAST('dem')
op <- par()
plot(dem_sp)

# Cargar a R el mapa vectorial de una cuenca que se encuentra alojado fuera de GRASS, hacer el plot y representar la cuenca del rio Guayubin superpuesta
library(sf)
rutaguayubin <- 'datos-fuente/cuenca_guayubin.geojson'
guayubin <- st_read(rutaguayubin)
plot(dem_sp)
plot(guayubin, add=T, col='transparent', border='black', lwd=5);par(op[c('mfrow','mar')])

# Analizar el DEM dentro de la cuenca de guayubin
library(raster)
dem_r0 <- raster(dem_sp)
dem_r1 <- crop(dem_r0, guayubin)
dem_guayu <- mask(dem_r1, guayubin)
plot(dem_guayu)
summary(dem_guayu)
hist(dem_guayu)

# Obtener variables de terreno básicas con el paquete raster dentro de R
pend_guayu <- terrain(x = dem_guayu, opt = 'slope', unit = 'degrees')
plot(pend_guayu)
summary(pend_guayu)
hist(pend_guayu)

# Obtener la misma variable de terreno con GRASS GIS
writeVECT(as_Spatial(guayubin), 'guayubin', v.in.ogr_flags='quiet')
execGRASS(
  "g.region",
  parameters=list(
    vector = "guayubin"
  )
)

execGRASS(
  "r.mask",
  flags = c('verbose','overwrite','quiet'),
  parameters = list(
    vector = 'guayubin'
  )
)

execGRASS(
  cmd = 'r.slope.aspect',
  flags = c('overwrite','quiet'),
  parameters = list(
    elevation='dem',
    slope='slope',
    aspect='aspect',
    pcurvature='pcurv',
    tcurvature='tcurv')
)

pend_guayu_g <- readRAST('slope')
plot(pend_guayu_g);par(op[c('mfrow','mar')])
summary(pend_guayu_g)
summary(pend_guayu)
gmeta()

execGRASS(
  "g.region",
  parameters=list(
    raster = "dem"
  )
)

execGRASS(
  "r.mask",
  flags = c('r','quiet')
)

gmeta()


# Video 6. Calcular parámetros hidrográficos con r.watershed. Visualizar con leaflet ----

# Imprimir lista de mapas ráster y vectoriales dentro en la región/localización activa 
# Está en el archivo reusable como (# Imprimir fuentes en la region)
execGRASS(
  'g.list',
  flags = 't',
  parameters = list(
    type = c('raster', 'vector')
  )
)
# Calcular parámetros hidrográficos de interés usando r.watershed
execGRASS(
  "r.watershed",
  flags = c('overwrite','quiet'),
  parameters = list(
    elevation = "dem",
    accumulation = "accum-de-rwshed",
    stream = "stream-de-rwshed",
    drainage = "drainage-dir-de-rwshed",
    basin = 'basins',
    half_basin = 'half-basins',
    threshold = 80
  )
)

## Traer capas a R

# Usar Spatial*
library(sp)
use_sp()
# Paquete manejo de los raster
library(raster)
# DEM
dem <- raster(readRAST('dem'))
# Basins
basins <- raster(readRAST('basins'))
# Stream network
stream <- raster(readRAST('stream-de-rwshed'))
stream3857 <- projectRaster(stream, crs = CRS("+init=epsg:3857"), method = 'ngb')
# Generar un vectorial de extensión de capa en EPSG:4326
e <- extent(stream)
e <- as(e, 'SpatialPolygons')
proj4string(e) <- CRS("+init=epsg:32619")
e <- spTransform(e, CRSobj = CRS("+init=epsg:4326"))

# Visualizar capas con leaflet
library(leaflet)
library(leafem)
leaflet() %>%
  addProviderTiles(providers$Stamen.Terrain, group = 'terrain') %>%
  addRasterImage(dem, group='DEM', opacity = 0.5) %>%
  addRasterImage(
    ratify(basins),
    group='basins', opacity = 0.7,
    colors = sample(rep(RColorBrewer::brewer.pal(12, 'Set3'),1000))) %>% 
  addRasterImage(stream3857, project = F, group='str', opacity = 0.7, method = 'ngb', colors = 'blue') %>% 
  addLayersControl(
    overlayGroups = c('terrain','DEM','basins','str'),
    options = layersControlOptions(collapsed=FALSE)) %>% 
  addHomeButton(extent(e), 'Ver todo')


# Video 7, Extraer una cuenca con r.water.outlet. Visualizar con mapview y leaflet ----
# Imprimir lista de mapas ráster y vectoriales dentro en la región/localización activa (está en el reproducible)
execGRASS(
  'g.list',
  flags = 't',
  parameters = list(
    type = c('raster', 'vector')
  )
)

# Obtener las coordenadas de la desembocadura de la cuenca de interés
library(mapview)
mapview(
  stream3857, method='ngb', col.regions = 'blue',
  legend = FALSE, label = FALSE, maxpixels =  910425
)

# Convertir las coordenadas lat/lon a EPSG:32619
my_trans <- function(coords = NULL) {
  require(sp)
  pt <- SpatialPoints(matrix(coords, ncol = 2), CRS("+init=epsg:4326"))
  foo <- spTransform(pt, CRSobj = CRS("+init=epsg:32619"))
  bar <- as.vector(coordinates(foo))
  return(bar)
}
guayu_out <- my_trans(coords = c(-71.40246,19.67306))
guayu_out

# Extraer la cuenca de interés
execGRASS(
  "r.water.outlet",
  flags = c('overwrite','quiet'),
  parameters = list(
    input = 'drainage-dir-de-rwshed',
    output = 'guayubin-basin',
    coordinates = guayu_out
  )
)

# Convertir la cuenca a vectorial en GRASS
execGRASS(
  "r.to.vect",
  flags = c('overwrite','quiet'),
  parameters = list(
    input = 'guayubin-basin',
    output = 'guayubin_basin',
    type = 'area'
  )
)

# Mostrar lista nuevamente
execGRASS(
  'g.list',
  flags = 't',
  parameters = list(
    type = c('raster', 'vector')
  )
)

# Traer a R la cuenca del rio guayubin
guayu_bas <- readVECT('guayubin_basin')
guayu_bas
plot(guayu_bas)
guayu_bas4326 <- spTransform(guayu_bas, CRSobj = CRS("+init=epsg:4326"))
leaflet() %>% 
  addProviderTiles(providers$Stamen.Terrain) %>%
  addRasterImage(stream, opacity = 0.7, method = 'ngb', colors = 'blue') %>% 
  addPolygons(data = guayu_bas4326) %>% 
  leafem::addHomeButton(extent(guayu_bas4326), 'Ver cuenca')

# Video 8, Extraer una red drenaje con r.stream.extract. Visualizar con leaflet ----
# Imprimir lista de mapas ráster y vectoriales dentro en la región/localización activa
execGRASS(
  'g.list',
  flags = 't',
  parameters = list(
    type = c('raster', 'vector')
  )
)

# Usar la cuenca del rio guayubin como máscara
execGRASS(
  "r.mask",
  flags = c('verbose','overwrite','quiet'),
  parameters = list(
    vector = 'guayubin_basin'
  )
)

# Extraer la red de drenaje de la cuenca de interés
execGRASS(
  "r.stream.extract",
  flags = c('overwrite','quiet'),
  parameters = list(
    elevation = 'dem',
    threshold = 80,
    stream_raster = 'guayubin-stream-de-rstr',
    stream_vector = 'guayubin_stream_de_rstr'
  )
)

# Mostrar lista nuevamente
execGRASS(
  'g.list',
  flags = 't',
  parameters = list(
    type = c('raster', 'vector')
  )
)

# Traer a R la red de drenaje del rio guayubin
guayu_net <- readVECT('guayubin_stream_de_rstr', ignore.stderr = T)
guayu_net
plot(guayu_net)
guayu_net4326 <- spTransform(guayu_net, CRSobj = CRS("+init=epsg:4326"))
guayu_net4326
guayu_centroid <- coordinates(rgeos::gCentroid(guayu_bas4326))
guayu_centroid
guayu_net_r <- raster(readRAST('guayubin-stream-de-rstr'))
guayu_net_r
guayu_net_r3857 <- projectRaster(guayu_net_r, crs = CRS("+init=epsg:3857"), method = 'ngb')
guayu_net_r3857
leaflet() %>% 
  setView(lng = guayu_centroid[1], lat = guayu_centroid[2], zoom = 11) %>%
  addProviderTiles(providers$Stamen.Terrain, group = 'terrain') %>%
  addRasterImage(guayu_net_r3857, opacity = 0.7, method = 'ngb', colors = 'grey20', group = 'str_raster') %>% 
  addPolylines(data = guayu_net4326, weight = 3, opacity = 0.7, group = 'str_vect') %>% 
  leafem::addHomeButton(extent(guayu_net4326), 'Ver todo') %>% 
  addLayersControl(
    overlayGroups = c('terrain','str_vect','str_raster'),
    options = layersControlOptions(collapsed=FALSE)) 

# Video 9, Orden de red y razón de bifurcación explicados ----
