#' Functions to filter occurrences based on a system of geographical units developped by The International Working Group on Taxonomic Databases for Plant Sciences (TDWG) at approximately "country" level and upwards.
#'
#' Clean the species occurrence locations not occurring within a given range.
#'
#' @param point_data A two-column data.frame of occurrence records coordinates.
#' @param species_name A two character string (Genus species) specifying the name of the species to check.
#' @param status A character string specifying the geographic range where to search for the species. Can be : 'native', 'introduced' or 'both'. Default is 'both'.
#' @param initial_level A numeric integer between 1 and 4. The initial unit level of the TWDG at which the data must be checked against. Default initial level is 2.
#' @param force.output A logical. Should the initial point dataset be returned in case where all points have been removed by the cleaning ? Default is FALSE.
#' @param sp A logical. Should the results be returned as a SpatialPointsDataFrame ? Default is TRUE. Ignored if \code{cleanOcc=FALSE}
#' @param ... Additional parameters to be passed to \code{TDWGinfo()}
#' @return A two-column data.frame or a SpatialPointsDataFrame object with point coordinates cleaned by the species geographic range defined by the TDWG and according to specified criteria.
#' @author <i.ondo@kew.org>
#' @export
coordinatesRangeCleaner <- function(point_data,
                                    species_name,
                                    species_id = NULL,
                                    status = 'both',
                                    initial_level = 2,
                                    force.output = FALSE,
                                    sp = FALSE,
                                    verbose = TRUE,...){

  call.fun 	<- match.call(expand.dots=TRUE)
  tmp.args    <- c("",'point_data','species_name','species_id','status','initial_level','force.output','sp','verbose','use_name_matching','full_data')
  call.tmp    <- call.fun[match(tmp.args, names(call.fun),nomatch=0)]
  call.tmp$verbose <- FALSE

  if(!status%in%c('native','introduced','both'))
    stop("Argument 'status' must be one of 'native', 'introduced' or 'both'.")

  cleaner 		<- switch(status, 'native'= "cleanByNativeDistrib", 'introduced'="cleanByIntroducedDistrib", 'both' ="cleanByKnownDistrib")
  call.tmp[[1]] 	<- as.name(cleaner)

  return(eval(call.tmp, parent.frame()))
}
