#' Functions developed to filter/select species/occurrences based on a system of geographical units
#' developped by The International Working Group on Taxonomic Databases for Plant Sciences (TDWG) at approximately "country" level and upwards.
#'
#' Generates a list of parameters that can be used to filter species and/or occurrence locations based on geographic units recording the native range of species
#'
#' @param point_data A two-column data.frame of occurrence records coordinates.
#' @param species_name A two character string (Genus species) specifying the name of the species to check.
#' @param species_id A character string specifying the id of the species in ipni or Kew world checklist.
#' @param grid_resol A numeric string specifying the resolution of the raster grid to be used when calculating the percent coverage of geographical units by point distribution.
#'        Must be given in m. Default is 10000 i.e. 10 km. Ignored if number 3 is specified in \code{which_skip} argument.
#' @param use_name_matching A logical. Should the species name be matched with the Kew checklist (i.e. look for direct correspondence) ? Default is FALSE.
#' @param initial_level A numeric integer between 1 and 4. The initial unit level of the TWDG at which the data must be checked against. Default is level 2.
#' @param by_id A logical. Should the criteria be applied to each geographical units ? Default is FALSE.
#' @param which_skip A numeric value or a vector of length 2 specifying which point distribution parameter(s) should be skipped from computations.
#'        Set to \code{2} to skip the computation of the proportion of tdwg units sampled; \code{3} to skip the computation of the percent coverage of the tdwg units sampled or\code{NULL} to compute all parameters.
#' @param recursive A logical. Should the criteria be also applied at the subsequent lower units level ? Default is FALSE.
#' @param sf A logical. Should the results be returned as a sf object ? Default is TRUE.
#' @param proj_utm A character string specifying the projection of the raster grid to be used when calculating the percent coverage of geographical units by point distribution.
#' @param A logical. Should the info be re-computed ? Default is FALSE. 
#' @param full_data A logical. Should the initial cleaned dataset be returned ? or the cleaned records only ? Default is FALSE i.e. only cleaned records are returned.
#' @return A sf object if \code{sf=TRUE}, a list object otherwise.
#' @author <i.ondo@kew.org>
#' @export
TDWGinfo <- function(point_data,
                     species_name,
                     species_id = NULL,
                     grid_resol = 10000,
                     status	=	'both',
                     use_name_matching = FALSE,
                     initial_level = 2,
                     backbone='wcvp',
                     by_id = FALSE,
                     which_skip = 3,
                     recursive = FALSE,
                     sf = TRUE,
                     proj_utm = "+proj=eqearth",
                     full_data = FALSE,
                     force = FALSE,
                     verbose = FALSE,...){

  if(verbose){
    message('#=================')
    message('#=0 Check inputs')
    message('#=================')
  }
  if(missing(point_data))
    stop("Argument 'point_data' is missing")
  if(!inherits(point_data,c("data.frame","data.table")))
    stop("Argument 'point_data' must be a data.frame or a data.table object")
  if(nrow(point_data)<1 || ncol(point_data)<2)
    stop(paste("Invalid format for the occurrence data.frame: point_data has",nrow(point_data),"rows and",ncol(point_data),"columns."))
  if(missing(species_name))
    stop("Argument 'species_name' is missing")
  if(!inherits(species_name,"character"))
    stop("Argument 'species_name' must be a character string.")
  if(nchar(species_name)==0L)
    stop("Invalid species name")

  call.fun 	<- match.call(expand.dots=FALSE)
  tmp.args    <- c("",'point_data','species_name', 'status','use_name_matching','grid_resol','initial_level','backbone','by_id','which_skip','recursive','sf','proj_utm','full_data')
  call.tmp    <- call.fun[match(tmp.args, names(call.fun),nomatch=0)]

  tdwg_info_fn <- file.path(tempdir(),gsub(" ","_",paste0("TDWGinfo_",species_name,".rds")))
  if(file.exists(tdwg_info_fn) & !force){
    tdwg_info  	<- readRDS(tdwg_info_fn)
    if(identical(call.tmp,tdwg_info[[1]]))
      return(tdwg_info[[2]])
  }

  if(verbose){
    message('#=============================')
    message('#= 1. Get species identifiers')
    message('#=============================')
  }
  #------------------------------------------
  #= 1.a get database id from the checklist
  #------------------------------------------
  if(backbone=='wcvp'){
    kew_checklist <- kew_wcvp
    kew_lookup_table <- kew_wcvp_distrib
    if(use_name_matching && data.table::key(kew_checklist)!="taxon_name") data.table::setkey(kew_checklist,"taxon_name")
    if(!use_name_matching && data.table::key(kew_checklist)!="taxon_name") data.table::setkey(kew_checklist,"taxon_name")
  }else{
    kew_checklist <- kew_powo
    kew_lookup_table <- kew_powo_distrib
    if(use_name_matching && data.table::key(kew_checklist)!="full_name_without_family") data.table::setkey(kew_checklist,"full_name_without_family")
    if(!use_name_matching && data.table::key(kew_checklist)!="acc_full_name_without_family") data.table::setkey(kew_checklist,"acc_full_name_without_family")
  }

  # if(is.na(hasDistrib(species_name,verbose=FALSE))){
  # if(verbose)
  # warning(paste("Unable to find a database id for species",species_name))
  # return(list(species=species_name, status=status, levels=c(initial_level, if(recursive & initial_level<4) initial_level+1), point_fraction=NULL, unit_fraction=NULL, range_filling=NULL, crs=NULL, cleaned_point_data = NULL))
  # }
  # if(!hasDistrib(species_name,verbose=FALSE)){
  # if(verbose)
  # warning(paste0("Distribution data for ",species_name," are not available in the TDWG database."))
  # return(list(species=species_name, status=status, levels=c(initial_level, if(recursive & initial_level<4) initial_level+1), point_fraction=NULL, unit_fraction=NULL, range_filling=NULL, crs=NULL, cleaned_point_data = NULL))
  # }

  # db_id <- try(kew_checklist[species_name][["accepted_db_id"]], silent=TRUE)
  # if(inherits(db_id,"try-error")||is.na(db_id)){
  #   cat('failed\n')
  #   if(verbose)
  #     warning(paste("Unable to find a database id for species",species_name))
  #   return(list(species=species_name, status=status, levels=c(initial_level, if(recursive & initial_level<4) initial_level+1), point_fraction=NULL, unit_fraction=NULL, range_filling=NULL, crs=NULL, cleaned_point_data = NULL))
  # }
  # if(verbose) cat("done\n")

  #----------------------------------------------
  #= 1.a get geographic code
  #----------------------------------------------
  if(verbose) message("> 1.a getting geographic code ...")
  # set geographic level and get polygons
  if(!inherits(initial_level,"numeric") || (initial_level<1 | initial_level>4)){
    message('failed')
    stop("Argument 'initial level' must be an integer between 1 and 4")
  }
  dots.args = list(...)

  switch(as.integer(initial_level),
         {
           geographic_level 	<- "continent_code_l1"
           level_code 			<- "LEVEL1_COD"
           geographic_polygon <- tdwg_level1
         },
         {
           geographic_level 	<- "region_code_l2"
           level_code 			<- "LEVEL2_COD"
           geographic_polygon <- tdwg_level2
         },
         {
           geographic_level 	<- "area_code_l3"
           level_code 			<- "LEVEL3_COD"
           geographic_polygon <- tdwg_level3
         },
         {
           geographic_level 	<- "area_code_l3" # not available yet
           level_code 			<- "level4_cod"
           geographic_polygon <- tdwg_level4
         }
  )
  do.recursive = recursive & initial_level<4
  if(do.recursive){
    sub_level = initial_level+1
    switch(as.integer(max(sub_level-1, 1)),
           {
             geographic_sub_level 	<- "region_code_l2"
             sub_level_code 			<- "LEVEL2_COD"
             sub_geographic_polygon <- tdwg_level2
           },
           {
             geographic_sub_level 	<- "area_code_l3"
             sub_level_code 			<- "LEVEL3_COD"
             sub_geographic_polygon <- tdwg_level3
           },
           {
             geographic_sub_level 	<- "area_code_l3" # not available yet
             sub_level_code 			<- "level4_cod"
             sub_geographic_polygon <- tdwg_level4
           }
    )
  }

  # subset or not according to native/introduced status
  if(!status%in%c('native','introduced','both'))
    stop("Argument 'status' must be one of 'native', 'introduced' or 'both'.")

  geographic_code <- get_geocode(species_name, species_id=species_id, status=status, level=initial_level, backbone=backbone, use_name_matching=use_name_matching)

  if(do.recursive)
    sub_geographic_code <- get_geocode(species_name, species_id=species_id, status=status, level=sub_level, backbone=backbone, use_name_matching=use_name_matching)

  # if(status =='both'){
  #   geographic_code <- try(kew_lookup_table[db_id][[geographic_level]],silent=TRUE)
  #   if(do.recursive)
  #     sub_geographic_code <- try(kew_lookup_table[db_id][[geographic_sub_level]],silent=TRUE)
  # }else if(status =='native'){
  #   geographic_code <- try(kew_lookup_table[db_id][introduced==0][[geographic_level]],silent=TRUE)
  #   if(do.recursive)
  #     sub_geographic_code <- try(kew_lookup_table[db_id][introduced==0][[geographic_sub_level]],silent=TRUE)
  # }else if(status =='introduced'){
  #   geographic_code <- try(kew_lookup_table[db_id][introduced==1][[geographic_level]],silent=TRUE)
  #   if(do.recursive)
  #     sub_geographic_code <- try(kew_lookup_table[db_id][introduced==1][[geographic_sub_level]],silent=TRUE)
  # }

  if(is.null(geographic_code)){
    if(verbose){
      message('failed')
      warning(paste0("Unable to find geo-data for species '",species_name,"' with the status '",status,"'"))
    }
    return(list(species=species_name, status=status, levels=c(initial_level, if(do.recursive) initial_level+1), point_fraction=NULL, unit_fraction=NULL, range_filling=NULL, crs=NULL, cleaned_point_data = NULL))
  }
  geographic_code = unique(geographic_code)
  if(do.recursive){
    sub_geographic_code	= unique(sub_geographic_code)
  }
  if(verbose) message("done\n")
  if(verbose){
    message('#=========================================================')
    message('#= 2. Get species polygons\n')
    message('#=========================================================')
  }
  if(verbose) message("> 2. getting species polygons ...")
  #tdwg_polygons <- eval(parse(text=sprintf("subset(geographic_polygon, %s %%in%% geographic_code)",level_code)))
  tdwg_polygons <- geographic_polygon %>%
    dplyr::filter(!!as.name(level_code) %in% geographic_code)
  stopifnot(inherits(tdwg_polygons,"sf"))
  if(do.recursive){
    tdwg_sub_level_polygons <- sub_geographic_polygon %>%
      dplyr::filter(!!as.name(sub_level_code) %in% sub_geographic_code)
  }
    #tdwg_sub_level_polygons <- eval(parse(text=sprintf("subset(sub_geographic_polygon, %s %%in%% sub_geographic_code)",sub_level_code)))
  
  if(verbose) message("done\n")
  if(!is.null(which_skip)){
    if(!all(inherits(which_skip,'numeric')) & !all(which_skip %in% c(2,3)) & length(which_skip)<=2)
      stop("Argument 'which.skip' must be an integer value or vector with value(s) 1, 2 or 3.")
  }
  if(verbose){
    if(!is.null(which_skip) & all(c(1,2,3) %in% which_skip) ){
      message("> skipping computation of points distribution parameters...\n")
    }else{
      message('#====================================================================')
      message('#= 3. Compute points distribution parameters')
      message('#====================================================================')
    }
  }

  #-----------------------------------
  #=3.a fraction of occurrence points
  #-----------------------------------

  point_data		<- setup_point_data(point_data, full_data=full_data)

  if(nrow(point_data)==0){
    if(verbose)
      warning(paste("No points available after quick setup for species:",species_name))
    initial_level_data <- as.data.frame(array(NA, dim=c(if(by_id & length(sf::st_geometry(tdwg_polygons))>1) length(sf::st_geometry(tdwg_polygons)) else 1, 4)))
    colnames(initial_level_data) <- c("frac_occ","frac_units","perc_cover","region")
    out = list(species=species_name, status=status, initial_level_data, crs= sf::st_crs(tdwg_polygons)$proj4string, cleaned_point_data = NULL)
    names(out)[3] <- paste0("level",initial_level)

    if(do.recursive){
      sub_level_data <- as.data.frame(array(0, dim=c(length(sf::st_geometry(tdwg_sub_level_polygons)), 4)))
      colnames(sub_level_data) <- c("frac_occ","frac_units","perc_cover","region")
      out = append(out, list(sub_level_data), after=3)
      names(out)[4] <- paste0("level",sub_level)
    }
    return(out)
  }
  nb_pts_total 	<- nrow(point_data)
  nb_tdwg_units <- length(sf::st_geometry(tdwg_polygons))
  id_x_lon 	<- grep(pattern = "[Ll][Oo][Nn]|^[Xx]$",x = names(point_data))[1]
  id_y_lat 	<- grep(pattern = "[Ll][Aa][Tt]|^[Yy]$",x = names(point_data))[1]
  coordHeaders <- c(id_x_lon, id_y_lat)
  sp.pts			  <- sf::st_as_sf(point_data[,coordHeaders], coords=c(1,2), crs=sf::st_crs(tdwg_polygons))
  is_in_range  	<- sf::st_intersects(sp.pts, sf::st_geometry(tdwg_polygons %>% dplyr::select(as.name(level_code))), sparse=TRUE)%>%lengths()>0L

  # create empty array to store points distribution parameters at initial level
  initial_level_data <- as.data.frame(array(0, dim=c(if(by_id & nb_tdwg_units>1) nb_tdwg_units else 1, 4)))
  colnames(initial_level_data) <- c("frac_occ","frac_units","perc_cover","region")

  do.frac_occ = !(!is.null(which_skip) & 1 %in% which_skip)
  if(verbose){
    if(!do.frac_occ)
      message("> skipping computation of the fraction of occurrences lying within TDWG units...\n")
  }
  if(do.frac_occ){
    if(verbose) message("> computing the fraction of occurrences lying within TDWG units...\n")

    point_fraction <- sum(is_in_range)/nb_pts_total
    if(point_fraction==0){
      if(verbose) message("No points has been found within the known range. Exiting...\n")
      out = list(species=species_name,
                 status=status,
                 eval(parse(text=sprintf("level%g=initial_level_data",initial_level))),
                 crs= sf::st_crs(tdwg_polygons)$proj4string,
                 cleaned_point_data = NULL)
      names(out)[3] <- paste0("level",initial_level)

      if(do.recursive){
        sub_level_data <- as.data.frame(array(0, dim=c(length(sf::st_geometry(tdwg_sub_level_polygons)), 4)))
        colnames(sub_level_data) <- c("frac_occ","frac_units","perc_cover","region")
        out = append(out, list(sub_level_data), after=3)
        names(out)[4] <- paste0("level",sub_level)
      }
      return(out)
    }

    if(!by_id){
      initial_level_data[1,1]	<- point_fraction
      initial_level_data[1,4] <- "Total"
    }
    else if(by_id){
      if(nb_tdwg_units>1){
        for(k in tdwg_polygons %>% dplyr::pull(as.name(level_code))){
          nb_pts_in_sub_range <- sum(lengths(sf::st_intersects(sp.pts, dplyr::filter(tdwg_polygons %>% dplyr::select(as.name(level_code)), !!as.name(level_code)==k), sparse=TRUE))>0L)
          initial_level_data[match(k,tdwg_polygons %>% dplyr::pull(as.name(level_code))),c(1,4)] <- c(nb_pts_in_sub_range/nb_pts_total,k)
        }
      }
      else{
        initial_level_data[1,1]	<- point_fraction
        initial_level_data[1,4] <- "Total"
      }
    }

    if(do.recursive){
      nb_tdwg_sub_units 	<- length(sf::st_geometry(tdwg_sub_level_polygons))
      # create empty array to store points distribution parameters at sub level
      sub_level_data <- as.data.frame(array(0, dim=c(nb_tdwg_sub_units, 4)))
      colnames(sub_level_data) <- c("frac_occ","frac_units","perc_cover","region")

      is_in_sub_units    <- vector(length=nb_tdwg_sub_units)
      for(k in tdwg_sub_level_polygons %>% dplyr::pull(as.name(sub_level_code))){
        nb_pts_in_sub_units <- sum(lengths(sf::st_intersects(sp.pts, dplyr::filter(tdwg_sub_level_polygons %>% dplyr::select(as.name(sub_level_code)), !!as.name(sub_level_code)==k), sparse=TRUE))>0L)
        sub_level_data[match(k,tdwg_sub_level_polygons %>% dplyr::pull(as.name(sub_level_code))),c(1,4)] <- c(nb_pts_in_sub_units/nb_pts_total, k)
      }
    }
    if(verbose) message("done")
  }

  #-----------------------------------
  #=3.b proportion of tdwg units
  #-----------------------------------
  do.frac_units = !(!is.null(which_skip) & 2 %in% which_skip)
  if(verbose){
    if(!do.frac_units)
      message("> skipping computation of the proportion of the tdwg units sampled ...\n")
  }
  if(do.frac_units){
    if(verbose) message("> computing the proportion of TDWG units in which the species is found...")

    is_in_units <- vector(length=nb_tdwg_units)
    for(k in tdwg_polygons %>% dplyr::pull(as.name(level_code)))
      is_in_units[match(k,tdwg_polygons %>% dplyr::pull(as.name(level_code)))] <- any(lengths(sf::st_intersects(sp.pts, dplyr::filter(tdwg_polygons %>% dplyr::select(as.name(level_code)), !!as.name(level_code)==k), sparse=TRUE))>0L)

    if(by_id)
      initial_level_data[,2] <- rep(sum(is_in_units)/nb_tdwg_units, times=nb_tdwg_units)
    else
      initial_level_data[1,2] <- sum(is_in_units)/nb_tdwg_units

    if(do.recursive){
      nb_tdwg_sub_units <- length(tdwg_sub_level_polygons)
      is_in_sub_units   <- vector(length=nb_tdwg_sub_units)
      for(k in tdwg_sub_level_polygons %>% dplyr::pull(as.name(sub_level_code))){
        is_in_sub_units[match(k,tdwg_sub_level_polygons %>% dplyr::pull(as.name(sub_level_code)))] <- any(lengths(sf::st_intersects(sp.pts, dplyr::filter(tdwg_sub_level_polygons %>% dplyr::select(as.name(sub_level_code)), !!as.name(sub_level_code)==k), sparse=TRUE))>0L)
      }
      sub_level_data[,2] <- sum(is_in_sub_units)/nb_tdwg_sub_units
    }

    if(verbose) message("done")
  }

  #-----------------------------------
  #=3.c percent coverage of tdwg units
  #-----------------------------------
  do.perc_cover = !(!is.null(which_skip) & 3 %in% which_skip)
  if(verbose){
    if(!do.perc_cover)
      message("> skipping computation of the percent coverage of the tdwg units ...\n")
  }
  if(do.perc_cover){
    if(verbose) message("> computing the percent coverage of the tdwg units by occurrence points...")

    # keep point inside the tdwg
    sp.pts_in_range	<- sp.pts[is_in_range,] # should have >0 pts since the fraction of point inside the tdwg unit >0
    # transform in utm
    sp.pts_utm_in_range <- sf::st_transform(sp.pts_in_range, sf::st_crs(proj_utm))
    # create raster template
    r_init	<- if(!exists('r_init')) { r_init <<- terra::init(terra::rast(),fun=0);r_init} else r_init
    r_temp_fn 		<- file.path(tempdir(),paste0(grid_resol,"_",gsub(" ","_",proj_utm),".rds"))
    if(file.exists(r_temp_fn))
      r_temp <- readRDS(r_temp_fn)
    else{
      r_temp  <- terra::project(r_init, y=proj_utm, res=grid_resol)
      saveRDS(r_temp,file=r_temp_fn )
    }
    # compute the percent coverage of the tdwg area occupied by occurrence points
    tdwg_utm_polygons_fn <- file.path(tempdir(),paste0('tdwg_utm_polygons_',grid_resol,"_",gsub(" ","_",proj_utm),".rds"))
    if(file.exists(tdwg_utm_polygons_fn))
      tdwg_utm_polygons <- readRDS(tdwg_utm_polygons_fn)
    else{
      tdwg_utm_polygons	<- sf::st_transform(tdwg_polygons, sf::st_crs(proj_utm))
      saveRDS(tdwg_utm_polygons,file=tdwg_utm_polygons_fn )
    }
    # TODO: use terra::expanse() in the future
    identity <- function(x) return(x)
    one <- function(x) return(1)
    g <- terra::rasterize(terra::vect(tdwg_utm_polygons),r_temp,field=level_code)
    tdwg_units_area <- terra::expanse(r_temp, unit="km", byValue=TRUE, zones=g)[,"area"]
    gg <- terra::rasterize(terra::vect(sp.pts_utm_in_range),r_temp, fun=one)
    sp.pts_area_cell 	<- terra::expanse(r_temp, unit="km", byValue=TRUE, zones=gg)[,"area"]
    #cell_area			<- (terra::xres(r_temp)*terra::yres(r_temp))
    #tdwg_units_area		<- if(nb_tdwg_units>1) lengths(terra::extract(r_temp,terra::vect(tdwg_utm_polygons), cells=TRUE)[,'cells'])*cell_area else length(terra::extract(r_temp,terra::vect(tdwg_utm_polygons),cells=TRUE)[,'cells'])*cell_area
    #sp.pts_area_cell 	<- rep(cell_area, times=length(unique(terra::extract(r_temp, terra::vect(sp.pts_utm_in_range), cells=TRUE)[,'cells'])))

    if(!by_id){
      initial_level_data[1,3]	<- (sum(sp.pts_area_cell,na.rm=TRUE)/sum(tdwg_units_area))*100
    }
    else if(by_id){
      if(nb_tdwg_units>1){
        range_filling <- vector(length=nb_tdwg_units)
        for(k in tdwg_utm_polygons %>% dplyr::pull(as.name(level_code))){
          is_in_id_units   <- sf::st_intersects(sp.pts_utm_in_range, 
                                                dplyr::filter(tdwg_utm_polygons %>% dplyr::select(as.name(level_code)), !!as.name(level_code)==k), sparse=TRUE) %>% lengths() >0L
          if(any(is_in_id_units))
            initial_level_data[match(k,tdwg_utm_polygons %>% dplyr::pull(as.name(level_code))),3] <- (sum(sp.pts_area_cell[is_in_id_units],na.rm=TRUE)/tdwg_units_area[match(k,tdwg_utm_polygons %>% dplyr::pull(as.name(level_code)))])*100
        }
      }
      else
        initial_level_data[1,3]	<- (sum(sp.pts_area_cell,na.rm=TRUE)/sum(tdwg_units_area))*100
    }

    if(do.recursive){
      tdwg_sub_level_utm_polygons_fn <- file.path(dirname(rasterTmpFile()),paste0('tdwg_sub_level_utm_polygons_',grid_resol,"_",gsub(" ","_",proj_utm),".rds"))
      if(file.exists(tdwg_sub_level_utm_polygons_fn)){
        tdwg_sub_level_utm_polygons <- readRDS(tdwg_sub_level_utm_polygons_fn)
      }else{
        tdwg_sub_level_utm_polygons	<- sf::st_transform(tdwg_sub_level_polygons, sf::st_crs(proj_utm))
        saveRDS(tdwg_sub_level_utm_polygons,file=tdwg_sub_level_utm_polygons_fn)
      }
      g <- terra::rasterize(terra::vect(tdwg_sub_level_utm_polygons),r_temp,field=level_code)
      tdwg_sub_units_area <- terra::expanse(r_temp, unit="km", byValue=TRUE, zones=g)[,"area"]
      sub_range_filling 	<- vector(length=nb_tdwg_sub_units)
      for(k in tdwg_sub_level_utm_polygons %>% dplyr::pull(as.name(sub_level_code))){
        is_in_id_sub_units   <- sf::st_intersects(sp.pts_utm_in_range, 
                                                  dplyr::filter(tdwg_sub_level_utm_polygons %>% dplyr::select(as.name(sub_level_code)), !!as.name(level_code)==k), sparse=TRUE)%>% lengths()>0L
        if(any(is_in_id_sub_units))
          sub_level_data[match(k,tdwg_sub_level_utm_polygons %>% dplyr::pull(as.name(sub_level_code))),3] <- (sum(sp.pts_area_cell[is_in_id_sub_units],na.rm=TRUE)/tdwg_sub_units_area[match(k,tdwg_sub_level_utm_polygons %>% dplyr::pull(as.name(sub_level_code)))])*100
      }
    }

    if(verbose) message("done\n")
  }

  out =list(species=species_name, status=status, initial_level_data, crs=sf::st_crs(sp.pts)$proj4string, cleaned_point_data=if(sf){if(full_data) sf::st_sf(data=point_data[is_in_range,], geom=sf::st_geometry(sp.pts[is_in_range,])) else sp.pts[is_in_range,]} else if(full_data) point_data[is_in_range,] else sf::st_coordinates(sp.pts[is_in_range,]))
  names(out)[3] <- paste0("level",initial_level)

  if(do.recursive){
    out = append(out, list(sub_level_data), after=3)
    names(out)[4] <- paste0("level",sub_level)
  }

  # save info
  saveRDS(list(call.tmp,out), file=tdwg_info_fn)

  return(out)
}

