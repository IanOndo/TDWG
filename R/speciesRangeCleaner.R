#' Functions to filter species based on a system of geographical units developped by The International Working Group on Taxonomic Databases for Plant Sciences (TDWG) at approximately "country" level and upwards.
#'
#' Filter the species whose the occurrence distribution data don't match expected distribution data defined by the user.
#'
#' @param point_data A two-column data.frame of occurrence records coordinates.
#' @param species_name A two character string (Genus species) specifying the name of the species to check.
#' @param species_id A character string specifying the id of the species in ipni or Kew world checklist database.
#' @param status A character string specifying the geographic range where to search for the species. Can be : 'native', 'introduced' or 'both'
#' @param point_fraction A numeric between 0 and 1. The minimum fraction of occurrence records that must be included within the geographic range defined in \code{status}.
#' @param unit_fraction A numeric between 0 and 1. The minimum fraction of geographical units that must include occurrence records.
#' @param range_filling A numeric between 0 and 100. The minimum percent of the geographic area covered by the occurrence records distribution.
#' @param grid_resol A numeric string specifying the resolution of the raster grid to be used when calculation the percent coverage of geographical units by point distribution.
#'        Must be given in m. Default is 10000 i.e. 10 km
#' @param initial_level A numeric integer between 1 and 4. The initial unit level of the TWDG at which the data must be checked against. Default initial level is 2.
#' @param by_id A logical. Should the criteria be applied to each geographical units ? Default is FALSE. Not implemented yet.
#' @param recursive A logical. Should the criteria be also applied at the subsequent lower units level ? Default is FALSE.
#' @param sf A logical. Should the results be returned as a sf object ? Default is TRUE. Ignored if \code{cleanOcc=FALSE}
#' @param ... Additional parameters to be passed to \code{TDWGinfo()}
#' @return By default the function returns a logical:
#'         TRUE if the distribution of species occurrences match the \code{status} range required, given the criteria.
#'         FALSE otherwise.
#' @author <i.ondo@kew.org>
#' @export
speciesRangeCleaner <- function(point_data,
                                species_name,
                                species_id = NULL,
                                status = 'both',
                                point_fraction = 0.50,
                                unit_fraction = NULL,
                                range_filling = NULL,
                                grid_resol = 10000,
                                initial_level = 2,
                                backbone='wcvp',
                                by_id = FALSE,
                                recursive = FALSE,
                                sf = FALSE,
                                verbose = TRUE,...){

  call.fun 	<- match.call(expand.dots=TRUE)
  tmp.args    <- c("",'point_data','species_name','species_id','status','point_fraction', 'unit_fraction','range_filling','grid_resol','initial_level','backbone','by_id','recursive','verbose','use_name_matching','full_data')
  call.tmp    <- call.fun[match(tmp.args, names(call.fun),nomatch=0)]
  call.tmp$cleanOcc <- FALSE
  call.tmp$verbose <- FALSE

  if(!status%in%c('native','introduced','both'))
    stop("Argument 'status' must be one of 'native', 'introduced' or 'both'.")

  fun 			<- switch(status, 'native'= "hasNativeDistrib", 'introduced'="hasIntroducedDistrib", 'both' ="hasKnownDistrib")
  call.tmp[[1]] 	<- as.name(fun)
  return(eval(call.tmp, parent.frame()))
}
