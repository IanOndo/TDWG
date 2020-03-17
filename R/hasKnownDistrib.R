#' Functions to filter/select species/occurrences based on a system of geographical units
#' developped by The International Working Group on Taxonomic Databases for Plant Sciences (TDWG) at approximately "country" level and upwards.
#'
#' Tests whether or not the species occurrence locations are distributed over its known range.
#'
#' @param point_data A two-column data.frame of occurrence records coordinates.
#' @param species_name A two character string (Genus species) specifying the name of the species to check.
#' @param species_id A character string specifying the id of the species in ipni or Kew world checklist database.
#' @param point_fraction A numeric between 0 and 1. The minimum fraction of occurrence records that must be included within the geographic range defined in \code{status}.
#' @param unit_fraction A numeric between 0 and 1. The minimum fraction of geographical units that must include occurrence records.
#' @param range_filling A numeric between 0 and 100. The minimum percent of the geographic area covered by the occurrence records distribution.
#' @param grid_resol A numeric string specifying the resolution of the raster grid to be used when calculation the percent coverage of geographical units by point distribution.
#'        Must be given in m. Default is 10000 i.e. 10 km
#' @param initial_level A numeric integer between 1 and 4. The initial unit level of the TWDG at which the data must be checked against. Default is level 2.
#' @param by_id A logical. Should the criteria be applied to each geographical units ? Default is FALSE. Not implemented yet.
#' @param recursive A logical. Should the criteria be also applied at the subsequent lower units level ? Default is FALSE.
#' @param cleanOcc A logical. Should the dataset of cleaned occurrence records be returned also ? Default is FALSE.
#' @param sp A logical. Should the results be returned as a SpatialPointsDataFrame ? Default is TRUE. Ignored if \code{cleanOcc=FALSE}
#' @param ... Additional parameters to be passed to \code{TDWGinfo()}
#' @return By default the function returns a logical:
#'         \code{TRUE} if the distribution of species occurrences match the \code{status} range required, given the criteria;
#'         \code{FALSE} otherwise.
#'         If \code{cleanOcc=TRUE}, a list of two elements is returned: \code{$hasIntroducedDistrib} returns the default logical as described just above; and \code{$cleaned_point_data} returns
#'         the point dataset cleaned by the species introduced geographic range defined by the TDWG and according to specified criteria.
#' @author <i.ondo@kew.org>
#' @export
hasKnownDistrib <- function(point_data,
                            species_name,
                            species_id = NULL,
                            point_fraction = 0.50,
                            unit_fraction = NULL,
                            range_filling = NULL,
                            grid_resol = 10000,
                            initial_level = 2,
                            by_id = FALSE,
                            recursive = FALSE,
                            cleanOcc = FALSE,
                            sp = FALSE,
                            verbose = TRUE,...){

  # if(!hasDistrib(species_name, verbose=FALSE)){
  # if(verbose)
  # warning(paste0("Distribution data for ",species_name," are not available in the TDWG database."))
  # return(NA)
  # }

  call.fun 	<- match.call(expand.dots=TRUE)
  tmp.args    <- c("",'point_data','species_name','species_id','grid_resol','initial_level','by_id','recursive','sp','verbose','use_name_matching','full_data')
  call.tmp    <- call.fun[match(tmp.args, names(call.fun),nomatch=0)]
  call.tmp$status <- 'both'
  if(!is.null(range_filling))
    call.tmp$which_skip <- NULL
  call.tmp$sp		<- FALSE
  call.tmp$verbose<- FALSE
  call.tmp[[1]] 	<- as.name("TDWGinfo")
  tdwg_info		<- eval(call.tmp, parent.frame())

  #----------------
  #=0 Check inputs
  #----------------
  if(!inherits(point_fraction,"numeric")||(point_fraction<0 | point_fraction>1))
    stop("Argument 'point_fraction' must be a fraction number between 0 and 1.")
  if(!is.null(unit_fraction)){
    if(!inherits(unit_fraction,"numeric")||(unit_fraction<0 | unit_fraction>1))
      stop("Argument 'unit_fraction' must be a fraction number between 0 and 1.")
  }
  if(!is.null(range_filling)){
    if(!inherits(range_filling,"numeric")||(range_filling<0 | range_filling>100))
      stop("Argument 'range_filling' must be a number between 0 and 100.")
  }
  # test if the fraction of occurrence condition is satisfied
  if(recursive)
    test_point_fraction = any(tdwg_info[[4]][,'frac_occ']>=point_fraction)
  else
    test_point_fraction = tdwg_info[[3]][,'frac_occ']>=point_fraction
  test_known = test_point_fraction

  # test if the proportion of tdwg units condition is satisfied
  if(!is.null(unit_fraction)){
    if(recursive)
      test_unit_fraction = any(tdwg_info[[4]][,'frac_units']>=unit_fraction)
    else
      test_unit_fraction = tdwg_info[[3]][,'frac_units']>=unit_fraction
    test_known = test_known & test_unit_fraction
  }
  # test if the coverage of tdwg units condition is satisfied
  if(!is.null(range_filling)){
    if(recursive)
      test_range_filling = any(sapply(tdwg_info[[4]][,'perc_cover'], function(x) any(x>=range_filling)))
    else
      test_range_filling = tdwg_info[[3]][,'perc_cover']>=range_filling
    test_known = test_known & test_range_filling
  }

  if(cleanOcc && !is.null(tdwg_info$cleaned_point_data)){
    if(sp){

      id_x_lon 	<- grep(pattern = "[Ll][Oo][Nn]",x = names(tdwg_info$cleaned_point_data))[1]
      id_y_lat 	<- grep(pattern = "[Ll][Aa][Tt]",x = names(tdwg_info$cleaned_point_data))[1]
      coordHeaders <- c(id_x_lon, id_y_lat)

      return(
        SpatialPointsDataFrame(	coords=tdwg_info$cleaned_point_data[, coordHeaders],
                                data=data.frame(species=if(recursive) rep(species_name, times= nrow(tdwg_info[[4]])) else species_name,
                                                status=if(recursive) rep(status, times=nrow(tdwg_info[[4]])) else status,
                                                hasKnownDistrib=if(recursive) rep(test_known,times=nrow(tdwg_info[[4]])) else test_known,
                                                levels=if(recursive) rep(initial_level, times=nrow(tdwg_info[[4]])) else min(initial_level+1, 4),
                                                point_fraction=if(recursive) tdwg_info[[4]][,'frac_occ'] else tdwg_info[[3]][,'frac_occ'],
                                                unit_fraction=if(recursive) tdwg_info[[4]][,'frac_units'] else tdwg_info[[3]][,'frac_units'],
                                                range_filling=if(recursive) tdwg_info[[4]][,'perc_cover'] else tdwg_info[[3]][,'perc_cover']
                                ),
                                proj4string=tdwg_info$crs
        )
      )
    }
    return(list(hasKnownDistrib=test_known, crs=tdwg_info$crs, cleaned_point_data=tdwg_info$cleaned_point_data))
  }
  return(test_known)
}

