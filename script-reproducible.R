
# Intro a Grass, video no. 3 ----
# Parte reutilizable del script ----
# Cargar paquetes
library(rgrass7)
library(sp)
library(sf)

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

# Retomar región de Grass Gis creada en pasos previos ----

{r; include=FALSE}
source(
  knitr::purl(
    'proyection-importar-fuente-extension.Rmd',
    output=tempfile()
  )
)


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
