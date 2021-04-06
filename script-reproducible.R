
# Intro a Grass, video no. 3 ----
## Parte reutilizable del script ----
## Cargar paquetes
library(rgrass7)
library(sp)
use_sp()
library(sf)
library(raster)
library(leaflet)
library(leafem)
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

## Imprimir fuentes en la region
execGRASS(
  'g.list',
  flags = 't',
  parameters = list(
    type = c('raster', 'vector')
  )
)

## Limpiar archivo de bloqueo del conjunto de mapas de GRASS
unlink_.gislock()

## Fin de la parte reutilizable

# Video no. 4, Definir proyección de la región de GRASS GIS, importar fuente y utilizarla para definir extensión y resolución. Cómo ver la ayuda de las funciones ----
## Muestra la definición de la región
gmeta()

## Definir ruta del DEM
dem <- 'datos-fuente/srtm_dem_cuenca_guayubin.tif'
execGRASS(
  cmd = 'g.proj',
  flags = c('t','c'),
  georef = dem)

## Muestra la definición de la región modificada
gmeta()

## r.in.gdal importa la fuente a GRASS
execGRASS(
  cmd = 'r.in.gdal',
  flags=c('overwrite','quiet'),
  parameters=list(
    input=dem,
    output='dem'
  )
)

## Actualizar la extensión de la región al DEM, sólo por precaución
execGRASS(
  cmd = 'g.region',
  parameters=list(
    raster = 'dem',
    align = 'dem'
  )
)

## Muestra la región de Grass again
gmeta()

## Importar vectorial a la región de Grass
demext <- 'datos-fuente/srtm_dem_cuenca_guayubin.geojson'
execGRASS(
  cmd = 'v.in.ogr',
  flags=c('overwrite','quiet'),
  parameters=list(
    input = demext,
    output = 'dem_extent'
  )
)

## Imprimir lista de mapas ráster y vectoriales dentro en la región/localización activa
execGRASS(
  'g.list',
  flags = 't',
  parameters = list(
    type = c('raster', 'vector')
  )
)

## Ver los addons disponibles en el repositorio oficial de GRASS GIS, incluyendo descripción
execGRASS(
  cmd = 'g.extension',
  flags = 'c'
)

## Consultar la ayuda de una función
parseGRASS("r.in.gdal")

## Consultar la ayuda de una función. Segunda alternativa
system('r.in.gdal --help')

## Limpiar archivo de bloqueo del conjunto de mapas de GRASS
unlink_.gislock()

# Video no. 5, Explorar datos espaciales básicos entre GRASS y R ----
## Imprimir lista de mapas ráster y vectoriales dentro en la región/localización activa
execGRASS(
  'g.list',
  flags = 't',
  parameters = list(
    type = c('raster', 'vector')
  )
)

## Cargar en R el DEM (mapa ráster)
library(sp)
use_sp()
dem_sp <- readRAST('dem')
op <- par()
plot(dem_sp)

## Cargar a R el mapa vectorial de una cuenca que se encuentra alojado fuera de GRASS, hacer el plot y representar la cuenca del rio Guayubin superpuesta
library(sf)
rutaguayubin <- 'datos-fuente/cuenca_guayubin.geojson'
guayubin <- st_read(rutaguayubin)
plot(dem_sp)
plot(guayubin, add=T, col='transparent', border='black', lwd=5);par(op[c('mfrow','mar')])

## Analizar el DEM dentro de la cuenca de guayubin
library(raster)
dem_r0 <- raster(dem_sp)
dem_r1 <- crop(dem_r0, guayubin)
dem_guayu <- mask(dem_r1, guayubin)
plot(dem_guayu)
summary(dem_guayu)
hist(dem_guayu)

## Obtener variables de terreno básicas con el paquete raster dentro de R
pend_guayu <- terrain(x = dem_guayu, opt = 'slope', unit = 'degrees')
plot(pend_guayu)
summary(pend_guayu)
hist(pend_guayu)

## Obtener la misma variable de terreno con GRASS GIS
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

## Imprimir lista de mapas ráster y vectoriales dentro en la región/localización activa 
## Está en el archivo reusable como (# Imprimir fuentes en la region)
execGRASS(
  'g.list',
  flags = 't',
  parameters = list(
    type = c('raster', 'vector')
  )
)
## Calcular parámetros hidrográficos de interés usando r.watershed
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

## Usar Spatial*
library(sp)
use_sp()
## Paquete manejo de los raster
library(raster)
## DEM
dem <- raster(readRAST('dem'))
## Basins
basins <- raster(readRAST('basins'))
## Stream network
stream <- raster(readRAST('stream-de-rwshed'))
stream3857 <- projectRaster(stream, crs = CRS("+init=epsg:3857"), method = 'ngb')
## Generar un vectorial de extensión de capa en EPSG:4326
e <- extent(stream)
e <- as(e, 'SpatialPolygons')
proj4string(e) <- CRS("+init=epsg:32619")
e <- spTransform(e, CRSobj = CRS("+init=epsg:4326"))

## Visualizar capas con leaflet
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
## Imprimir lista de mapas ráster y vectoriales dentro en la región/localización activa (está en el reproducible)
execGRASS(
  'g.list',
  flags = 't',
  parameters = list(
    type = c('raster', 'vector')
  )
)

## Obtener las coordenadas de la desembocadura de la cuenca de interés
library(mapview)
mapview(
  stream3857, method='ngb', col.regions = 'blue',
  legend = FALSE, label = FALSE, maxpixels =  910425
)

## Convertir las coordenadas lat/lon a EPSG:32619
my_trans <- function(coords = NULL) {
  require(sp)
  pt <- SpatialPoints(matrix(coords, ncol = 2), CRS("+init=epsg:4326"))
  foo <- spTransform(pt, CRSobj = CRS("+init=epsg:32619"))
  bar <- as.vector(coordinates(foo))
  return(bar)
}
guayu_out <- my_trans(coords = c(-71.40021,19.66387))
guayu_out

## Extraer la cuenca de interés
execGRASS(
  "r.water.outlet",
  flags = c('overwrite','quiet'),
  parameters = list(
    input = 'drainage-dir-de-rwshed',
    output = 'guayubin-basin',
    coordinates = guayu_out
  )
)

## Convertir la cuenca a vectorial en GRASS
execGRASS(
  "r.to.vect",
  flags = c('overwrite','quiet'),
  parameters = list(
    input = 'guayubin-basin',
    output = 'guayubin_basin',
    type = 'area'
  )
)

## Mostrar lista nuevamente
execGRASS(
  'g.list',
  flags = 't',
  parameters = list(
    type = c('raster', 'vector')
  )
)

## Traer a R la cuenca del rio guayubin
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
## Imprimir lista de mapas ráster y vectoriales dentro en la región/localización activa
execGRASS(
  'g.list',
  flags = 't',
  parameters = list(
    type = c('raster', 'vector')
  )
)

## Usar la cuenca del rio guayubin como máscara
execGRASS(
  "r.mask",
  flags = c('verbose','overwrite','quiet'),
  parameters = list(
    vector = 'guayubin_basin'
  )
)

## Extraer la red de drenaje de la cuenca de interés
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

## Mostrar lista nuevamente
execGRASS(
  'g.list',
  flags = 't',
  parameters = list(
    type = c('raster', 'vector')
  )
)

## Traer a R la red de drenaje del rio guayubin
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


# Video 10,  Orden de red y análisis hortoniano usando r.stream*. Visualizar con leaflet ----
## Imprimir lista de mapas ráster y vectoriales dentro en la región/localización activa (está en el reproducible)
execGRASS(
  'g.list',
  flags = 't',
  parameters = list(
    type = c('raster', 'vector')
  )
)

## Crear mapa de dirección de flujo a partir de r.stream
execGRASS(
  "r.stream.extract",
  flags = c('overwrite','quiet'),
  parameters = list(
    elevation = 'dem',
    threshold = 80,
    direction = 'drainage-dir-de-rstr'
  )
)

## Crear mapas de órdenes de red
execGRASS(
  "r.stream.order",
  flags = c('overwrite','quiet'),
  parameters = list(
    stream_rast = 'guayubin-stream-de-rstr',
    direction = 'drainage-dir-de-rstr',
    elevation = 'dem',
    accumulation = 'accum-de-rwshed',
    stream_vect = 'order_all',
    strahler = 'order-strahler',
    horton = 'order-horton',
    shreve = 'order-shreve',
    hack = 'order-hack-gravelius',
    topo = 'order-topology'
  )
)

## Visualizar la red con leaflet
## Simbología única
order <- readVECT('order_all')
order4326 <- spTransform(order, CRSobj = CRS("+init=epsg:4326"))
leaflet() %>% 
  addProviderTiles(providers$Stamen.Terrain, group = 'terrain') %>%
  addPolylines(
    data = order4326, weight = 3, opacity = 0.7, group = 'order',
    label = ~as.character(strahler),
    highlightOptions = highlightOptions(color = "white",
                                        weight = 5, bringToFront = F, opacity = 1),
    labelOptions = labelOptions(noHide = T,
                                style = list(
                                  "font-size" = "8px",
                                  "background" = "rgba(255, 255, 255, 0.5)",
                                  "background-clip" = "padding-box",
                                  "padding" = "1px"))) %>% 
  leafem::addHomeButton(extent(order4326), 'Ver todo') %>% 
  addLayersControl(
    overlayGroups = c('terrain','order'),
    options = layersControlOptions(collapsed=FALSE))

## Simbología aplicando grosor según orden de red
leaflet() %>% 
  addProviderTiles(providers$Stamen.Terrain, group = 'terrain') %>%
  addPolylines(
    data = order4326, weight = order4326$strahler*1.5, opacity = 0.7, group = 'order',
    label = ~as.character(strahler),
    highlightOptions = highlightOptions(color = "white",
                                        weight = 5, bringToFront = F, opacity = 1),
    labelOptions = labelOptions(noHide = F)) %>% 
  leafem::addHomeButton(extent(order4326), 'Ver todo') %>% 
  addLayersControl(
    overlayGroups = c('terrain','order'),
    options = layersControlOptions(collapsed=FALSE))

## Delimitar cuencas según orden de red de Strahler

## Obtener órdenes de red mínimo y máximo
## Estadísticas para obtener los valores mínimo y máximo del orden de red de Strahler
rinfo.ordstra <- execGRASS(
  'r.info',
  flags = 'r',
  parameters = list(
    map = 'order-strahler'
  )
)

## Órdenes de red mínimo y máximo
minmaxord <- as.numeric(
  stringr::str_extract_all(
    attributes(rinfo.ordstra)$resOut,
    "[0-9]+"
  )
)
minmaxord

## Delimitar cuencas, convertirlas de ráster a vectorial
sapply(
  min(minmaxord):max(minmaxord),
  function(x){
    execGRASS(
      "r.stream.basins",
      flags = c('overwrite','c','quiet'),
      parameters = list(
        direction = 'drainage-dir-de-rstr',
        stream_rast = 'order-strahler',
        cats = as.character(x),
        basins = paste0('r-stream-basins-',x)
      )
    )
    execGRASS(
      "r.to.vect",
      flags=c('overwrite','quiet'),
      parameters = list(
        input = paste0('r-stream-basins-',x),
        output = paste0('r_stream_basins_',x),
        type = 'area'
      )
    )
  }

)
## Representar las cuencas con leaflet
sapply(
  min(minmaxord):max(minmaxord),
  function(x){
    assign(
      paste0('orden', x),
      spTransform(readVECT(paste0('r_stream_basins_',x)), CRSobj = CRS("+init=epsg:4326")),
      envir = .GlobalEnv)
  }
)

paleta <- RColorBrewer::brewer.pal(12, 'Set3')
leaflet() %>% 
  addProviderTiles(providers$Stamen.Terrain, group = 'terrain') %>%
  addPolygons(data = orden5, stroke = T, weight = 2,
              color = ~paleta, fillOpacity = 0.4, group = 'O5') %>% 
   addPolygons(data = orden4, stroke = T, weight = 2,
              color = ~paleta, fillOpacity = 0.4, group = 'O4') %>% 
  addPolygons(data = orden3, stroke = T, weight = 2,
              color = ~paleta, fillOpacity = 0.4, group = 'O3') %>%
  addPolygons(data = orden2, stroke = T, weight = 2,
              color = ~paleta, fillOpacity = 0.4, group = 'O2') %>%
  addPolygons(data = orden1, stroke = T, weight = 2,
              color = ~paleta, fillOpacity = 0.4, group = 'O1') %>%
  addPolylines(
    data = order4326, weight = order4326$strahler*1.5,
    opacity = 0.7, group = 'str_order') %>%
  leafem::addHomeButton(extent(order4326), 'Ver todo') %>% 
  addLayersControl(
    overlayGroups = c('terrain','O1','O2','O3','O4','O5','str_order'),
    options = layersControlOptions(collapsed=FALSE))

## Estadísticas de red resumidas por orden de red
execGRASS(
  "r.stream.stats",
  flags = c('overwrite','quiet','o'),
  parameters = list(
    stream_rast = 'order-strahler',
    direction = 'drainage-dir-de-rstr',
    elevation = 'dem',
    output = 'guayu_stats.txt'
  )
)
file.show('guayu_stats.txt')
d <- read.csv("guayu_stats.txt", skip=1, header=TRUE)
plot(num_of_streams~order, data=d, log="y")
mod <- lm(log10(num_of_streams)~order, data=d)
abline(mod)
text(2, 20, 'logN=2.064-0.544u')
rb <- 1/10^mod$coefficients[[2]]
rb

## Estadísticas de red ampliadas
execGRASS(
  "r.stream.stats",
  flags = c('overwrite','quiet'),
  parameters = list(
    stream_rast = 'order-strahler',
    direction = 'drainage-dir-de-rstr',
    elevation = 'dem',
    output = 'guayu_stats_expanded.txt'
  )
)
file.show('guayu_stats_expanded.txt')


# Video 11, Calcular índices de concavidad y perfiles longitudinales de cursos fluviales----
## Imprimir lista de mapas ráster y vectoriales dentro en la región/localización activa
execGRASS(
  'g.list',
  flags = 't',
  parameters = list(
    type = c('raster', 'vector')
  )
)

## Obtener coordenada
mapview(order, col.regions = 'blue', legend = FALSE)

## Obtener cursos más largos (cargar función propia)
devtools::source_url('https://raw.githubusercontent.com/geofis/rgrass/master/lfp_network.R') #Cargada como función "LfpNetwork"
LfpNetwork(
  xycoords = my_trans(c(-71.40047,19.662755)),
  suffix = 'Gua',
  stream_vect = 'order_all',
  direction = 'drainage-dir-de-rstr'
)

##Imprimir lista de mapas ráster y vectoriales
execGRASS(
  'g.list',
  flags = 't',
  parameters = list(
    type = c('raster', 'vector')
  )
)

## Representar con leaflet
lfp <- readVECT('LfpNetwork_lfp_all_final_Gua')
lfp4326 <- spTransform(lfp, CRSobj = CRS("+init=epsg:4326"))
leaflet() %>%
  addProviderTiles(providers$Stamen.Terrain, group = 'terrain') %>%
  addPolylines(
    data = lfp4326, weight = 3, opacity = 0.7, group = 'order',
    label = ~as.character(cat),
    highlightOptions = highlightOptions(color = "white",
                                        weight = 5, bringToFront = F, opacity = 1),
    labelOptions = labelOptions(noHide = T,
                                style = list(
                                  "font-size" = "8px",
                                  "background" = "rgba(255, 255, 255, 0.5)",
                                  "background-clip" = "padding-box",
                                  "padding" = "1px"))) %>% 
  leafem::addHomeButton(extent(lfp4326), 'Ver todo')

## Exportar a KML
execGRASS(
  'v.out.ogr',
  flags = c('overwrite','quiet'),
  parameters = list(
    input = 'LfpNetwork_lfp_all_final_Gua',
    output = 'lfp_kml.kml',
    format = 'KML',
    dsco = 'NameField=cat'
  )
)

## Obtención de perfiles longitudinales e índices de concavidad
source('lfp_profiles_concavity.R') #Cargado como función "LfpProfilesConcavity"
guaybin_conv_prof <- LfpProfilesConcavity(
  xycoords = my_trans(c(-71.40047,19.662755)),
  network = 'LfpNetwork_lfp_all_final_Gua',
  prefix = 'Gyb',
  dem = 'dem',
  direction = 'drainage-dir-de-rstr',
  crs = '+init=epsg:32619',
  smns = 0.5,
  nrow = 3)

##Mostrar resultados
guayubin_conv_prof$profiles
guayubin_conv_prof$concavityindex
guayubin_conv_prof$dimensionlessprofiles

## Tabla dx/dy, tanto en metros como adimensional. Útiles para construir perfiles por cuenta propia
guayubin_conv_prof$lengthzdata %>% tibble::as.tibble()
guayubin_conv_prof$lengthzdatadmnls %>% tibble::as.tibble()