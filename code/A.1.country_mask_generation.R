library( 'sp' )
library( 'raster' )
library( 'rgeos' )
library( 'rgdal' )
library( 'geosphere' )
library( 'maptools' )
library( 'plyr' )

source( './code/common/mask_functions.R' )

# ----------------------------------------------------
# 0. basic settings

grid_resolution <- 0.5
mask_out_dir_name <- paste0( paste0( 'mask_out_', gsub( '.', '', grid_resolution, fixed = T ) ) )
mask_out <- paste0( './intermediate-out/', mask_out_dir_name )
if( !dir.exists( mask_out ) ) { 
dir.create( mask_out, showWarnings = TRUE, recursive = FALSE, mode = "0777" )
}

supplements_out <- './intermediate-out/supplements_out' 

# ----------------------------------------------------
# 1. load in the shapefile 
start <- proc.time() 
worldshp <- shapefile( './input/shp/CEDS_country_scheme/country_boundary.shp' )
time_passed <- proc.time( ) - start
print( paste0( 'time used for loading the shp ', time_passed[ 3 ], ' s' ) )


# ----------------------------------------------------
# 2. generate location index table 
locationIndexGeneration <- function( shp, grid_res ) { 
  
  # extract bbox for each feature in the given shapefile
  total_feature_number <- nrow( shp@data )
  bbox_res_list <- lapply( 1 : total_feature_number, function( i ) { 
    bbox_mat <- shp[ i, ]@bbox
    xmin <- round( bbox_mat[ 'x', 'min' ], 6 ) 
    ymin <- round( bbox_mat[ 'y', 'min' ], 6 ) 
    xmax <- round( bbox_mat[ 'x', 'max' ], 6 ) 
    ymax <- round( bbox_mat[ 'y', 'max' ], 6 ) 
    out_df <- data.frame( xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax, stringsAsFactors = F )
    } )
  bbox_df <- do.call( 'rbind', bbox_res_list )
  
  # compute bounding coordinates for xmin, ymin, xmax, ymax 
  boundingMin <- function( min_coor, grid_res ) {
    factor <- 1 / grid_res
    min_bounding <- floor(  min_coor * factor ) / factor
    return( min_bounding )
  }
  boundingMax <- function( max_coor, grid_res ){
    factor <- 1 / grid_res
    max_bounding <- ceiling( max_coor * factor ) / factor
    return( max_bounding ) 
  }
  
  xmin_boundings <- unlist( lapply( bbox_df[ , 'xmin' ], boundingMin, grid_res ) )
  ymin_boundings <- unlist( lapply( bbox_df[ , 'ymin' ], boundingMin, grid_res ) )
  xmax_boundings <- unlist( lapply( bbox_df[ , 'xmax' ], boundingMax, grid_res ) )
  ymax_boundings <- unlist( lapply( bbox_df[ , 'ymax' ], boundingMax, grid_res ) )
  
  xmin_bounding_lut <- seq( -180, ( 180 - grid_res ), grid_res )
  xmax_bounding_lut <- seq( ( -180 + grid_res ), 180, grid_res )  
  ymin_bounding_lut <- seq( ( 90 - grid_res ) , -90, -grid_res )
  ymax_bounding_lut <- seq( 90, ( -90 + grid_res ), -grid_res )
  
  start_row <- unlist( lapply( ymax_boundings, function( ymax_bounding ) { 
    which( ymax_bounding_lut %in% ymax_bounding )
    } ) )
  end_row <- unlist( lapply( ymin_boundings, function( ymin_bounding ) { 
    which( ymin_bounding_lut %in% ymin_bounding)
    } ) )
  start_col <- unlist( lapply( xmin_boundings, function( xmin_bounding ) { 
    which( xmin_bounding_lut %in% xmin_bounding )
    } ) )
  end_col <- unlist( lapply( xmax_boundings, function( xmax_bounding ) { 
    which( xmax_bounding_lut %in% xmax_bounding )
    } ) )
  
  location_index_df <- data.frame( iso = shp@data$iso, 
                                   start_row = start_row, 
                                   end_row = end_row,
                                   start_col = start_col, 
                                   end_col = end_col,
                                   xmin_bounding = xmin_boundings,
                                   xmax_bounding = xmax_boundings, 
                                   ymin_bounding = ymin_boundings, 
                                   ymax_bounding = ymax_boundings,
                                   stringsAsFactors = F )
  return( location_index_df )
  }

location_index <- locationIndexGeneration( worldshp, grid_resolution ) 

write.csv( location_index, paste0( './intermediate-out/supplements_out/', 'country_location_index_', gsub( '.', '', grid_resolution, fixed = T ), '.csv' ), row.names = F )

# ----------------------------------------------------
# 3. generate the country masks
region_list <- location_index$iso 

# define the temp mask generation function 
tempMaskGene <- function( region, location_index, grid_resolution ) { 
  
  indexline <- location_index[ location_index$iso == region, ]
  xmin_bounding <- indexline$xmin_bounding 
  ymin_bounding <- indexline$ymin_bounding 
  xmax_bounding <- indexline$xmax_bounding 
  ymax_bounding <- indexline$ymax_bounding 
  
  # generate TopoGrid
  # first define three parameters for TopoGrid
  topo_offset <- c( ( xmin_bounding + grid_resolution / 2 ), 
                    ( ymin_bounding + grid_resolution / 2 ) )
  topo_res <- c( grid_resolution, grid_resolution )
  col_row_num <- c( ( xmax_bounding - xmin_bounding ) / grid_resolution,
                    ( ymax_bounding - ymin_bounding ) / grid_resolution )
  col_row_num <- round( col_row_num )
  # then generate the TopoGrid
  grid <- GridTopology( topo_offset, topo_res, col_row_num )
  
  # convert the topogrid into SpatialPolygon
  poly <- as( grid, 'SpatialPolygons' )  
  
  # tempid in polydf: fake a id for later on join
  polydf <- SpatialPolygonsDataFrame( poly, 
                                      data.frame( 
                                        tempid = 1:( col_row_num[ 1 ] * 
                                                       col_row_num[ 2 ] ) ),
                                      match.ID = F)
  # select the region shp 
  regionshp <- worldshp[ worldshp$iso == region, ]
  regionshp <- gBuffer( regionshp, byid = T, width = 0 )
  unique( gIsValid( regionshp, byid = T, reason = T ) )

  # set projection
  projection( polydf ) <- projection( regionshp )
  
  # intersection
  intpoly <- intersect( polydf, regionshp[  regionshp$iso == region , ] )
  # calculate the area
  intareas <- areaPolygon( intpoly ) 
  # add the calculated area to the intpoly data frame 
  intpoly_df <- intpoly@data
  intpoly_df$area <- intareas
  
  # areas for the grd
  polyareas <-  areaPolygon( polydf )
  polydf_df <- polydf@data
  polydf_df$area <- polyareas
  
  # merge two attribute tables together by tempid
  polyatt <- merge( polydf_df, intpoly_df, by = 'tempid', all.x = T )
  polyatt <- polyatt[ order( polyatt$tempid ), ]
  
  calc_df <- polyatt[ , c( 'area.x', 'area.y' ) ]
  calc_df$area.y[ is.na( calc_df$area.y ) ] <- 0 
  calc_df$area_percent <- calc_df$area.y / calc_df$area.x
  
  percent_mat <- matrix( calc_df$area_percent, col_row_num[ 2 ], col_row_num[ 1 ], byrow = T )
  
  tempname <- paste0( region, '_mask' )
  assign( tempname, percent_mat )
  # save the masks as R objects
  save( list = tempname,
        file = paste0( mask_out, '/', tempname ) )
  
  print( region )
}

# lapply the function to region_list 
invisible( lapply( region_list, tempMaskGene, location_index, grid_resolution ) )


