# Mask Functions
library( 'sp' )
library( 'raster' )
library( 'geosphere' )

# ----------------------------------------------------------
# flip_a_matrix 
# Brief: generate a fliped matrix by a given matrix
flip_a_matrix <- function( x ) {
  apply( x, 2, rev )
}

# ----------------------------------------------------------
# rasterize_GlobalExt 
# Brief: generate a fliped matrix by a given matrix
rasterize_GlobalExt <- function( x, proj = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" ) {
  rasterized <- raster( x, xmn = -180, xmx = 180, ymn = -90, ymx = 90, crs = proj)
  return( rasterized )   
}

# ----------------------------------------------------------
