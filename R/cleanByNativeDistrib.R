#' Functions to filter occurrences based on a system of geographical units developped by The International Working Group on Taxonomic Databases for Plant Sciences (TDWG) at approximately "country" level and upwards.
#'
#' Clean the species occurrence locations not occurring within its native range.
#'
#' @param point_data A two-column data.frame of occurrence records coordinates.
#' @param species_name A two character string (Genus species) specifying the name of the species to check.
#' @param species_id A character string specifying the id of the species in ipni or Kew world checklist database.
#' @param initial_level A numeric integer between 1 and 4. The initial unit level of the TWDG at which the data must be checked against. Default initial level is 2.
#' @param force.output A logical. Should the initial point dataset be returned in case where all points have been removed by the cleaning ? Default is FALSE.
#' @param sp A logical. Should the results be returned as a SpatialPointsDataFrame ? Default is FALSE.
#' @param ... Additional parameters to be passed to \code{TDWGinfo()}
#' @return A two-column data.frame or a SpatialPointsDataFrame object with point coordinates cleaned by the species native geographic range defined by the TDWG and according to specified criteria.
#' @author <i.ondo@kew.org>
#' @export
cleanByNativeDistrib <- function(point_data,
                                 species_name,
                                 species_id = NULL,
                                 initial_level = 2,
                                 force.output = FALSE,
                                 sp = FALSE,
                                 verbose = TRUE,...){

  call.fun 	<- match.call(expand.dots=TRUE)
  tmp.args    <- c("",'point_data','species_name','species_id','initial_level','force.output','sp','verbose','use_name_matching','full_data')
  call.tmp    <- call.fun[match(tmp.args, names(call.fun),nomatch=0)]
  call.tmp$status <- 'native'
  call.tmp$which_skip <- c(1,2,3)
  call.tmp$sp <- sp
  call.tmp$verbose <- FALSE
  call.tmp[[1]] 	<- as.name("TDWGinfo")
  NativeDistrib	<- eval(call.tmp, parent.frame())
  #-----------------------------
  #= Non available distribution
  #-----------------------------
  if(!inherits(NativeDistrib,"list")){
    if(verbose)
      warning(paste0("Distribution data for ",species_name," are not available in the TDWG database."));flush.console()
    if(force.output){
      if(verbose)
        warning("Returning initial point dataset.");flush.console()
      return(point_data)
    }
    return(NA)
  }
  #-------------------------
  #= Non native distribution
  #-------------------------
  # if(!NativeDistrib$hasNativeDistrib){
  # if(verbose)
  # warning(paste0("Distribution data for ",species_name," are not distributed over its native geographic range according to the TDWG database and the given criteria."));flush.console()
  # if(force.output){
  # if(verbose)
  # warning("Returning initial point dataset.");flush.console()
  # return(point_data)
  # }
  # if(verbose)
  # warning("Returning NA.");flush.console()
  # return(NA)
  # }
  #-------------------------
  #= Native distribution
  #-------------------------
  # if(sp){
  # return(
  # SpatialPointsDataFrame(	coords = NativeDistrib$cleaned_point_data,
  # data = data.frame(species=species_name,
  # hasNativeDistrib=NativeDistrib$hasNativeDistrib
  # ),
  # proj4string=NativeDistrib$crs
  # )
  # )
  # }
  return(NativeDistrib$cleaned_point_data)
}

