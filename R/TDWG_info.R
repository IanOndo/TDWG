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
#' @param sp A logical. Should the results be returned as a SpatialPointsDataFrame ? Default is TRUE.
#' @param proj_utm A character string specifying the projection of the raster grid to be used when calculating the percent coverage of geographical units by point distribution.
#' @param full_data A logical. Should the initial cleaned dataset be returned ? or the cleaned records only ? Default is FALSE i.e. only cleaned records are returned.
#' @return A SpatialPointsDataFrame object if \code{sp=TRUE}, a list object otherwise.
#' @author <i.ondo@kew.org>
#' @export
TDWGinfo <- function(point_data,
                     species_name,
                     species_id = NULL,
                     grid_resol = 10000,
                     status	=	'both',
                     use_name_matching = FALSE,
                     initial_level = 2,
                     by_id = FALSE,
                     which_skip = 3,
                     recursive = FALSE,
                     sp = TRUE,
                     proj_utm = "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs",
                     full_data = FALSE,
                     verbose = FALSE,...){

  if(verbose){
    cat('#=================\n')
    cat('#=0 Check inputs\n')
    cat('#=================\n')
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
  tmp.args    <- c("",'point_data','species_name', 'status','use_name_matching','grid_resol','initial_level','by_id','which_skip','recursive','sp','proj_utm','full_data')
  call.tmp    <- call.fun[match(tmp.args, names(call.fun),nomatch=0)]

  tdwg_info_fn <- file.path(dirname(raster::rasterTmpFile()),gsub(" ","_",paste0("TDWGinfo_",species_name,".rds")))
  if(file.exists(tdwg_info_fn)){
    tdwg_info  	<- readRDS(tdwg_info_fn)
    if(identical(call.tmp,tdwg_info[[1]]))
      return(tdwg_info[[2]])
  }

  if(verbose){
    cat('#=============================\n')
    cat('#= 1. Get species identifiers\n')
    cat('#=============================\n')
  }
  #------------------------------------------
  #= 1.a get database id from the checklist
  #------------------------------------------
  # if(verbose) cat("> 1.a getting database id from the checklist ...")
  # if(!exists('kew_checklist',envir=.GlobalEnv)){
  #   if(is.null(checklist)){
  #     dots.args = list(...)
  #     if(!"checklist" %in% names(dots.args))
  #       stop(paste0("Unable to locate the species checklist. Please provide a valid path to Kew checklist file."))
  #     else if(file.exists(dot.args[['checklist']]))
  #       kew_checklist <<- data.table::fread(dot.args[['checklist']], key=if(use_name_matching) "full_name_without_family" else "acc_full_name_without_family",showProgress=FALSE)
  #   }else if(!file.exists(checklist))
  #     stop(paste0("Unable to locate file: ",checklist,". Please provide a valid path to Kew checklist file."))
  #   else kew_checklist <<- data.table::fread(checklist, key=if(use_name_matching) "full_name_without_family" else "acc_full_name_without_family",showProgress=FALSE)
  # }else{
    if(use_name_matching && data.table::key(kew_checklist)!="full_name_without_family") data.table::setkey(kew_checklist,"full_name_without_family")
    if(!use_name_matching && data.table::key(kew_checklist)!="acc_full_name_without_family") data.table::setkey(kew_checklist,"acc_full_name_without_family")
  # }

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
  if(verbose) cat("> 1.a getting geographic code ...")
  # set geographic level and get polygons
  if(!inherits(initial_level,"numeric") || (initial_level<1 | initial_level>4)){
    cat('failed\n')
    stop("Argument 'initial level' must be an integer between 1 and 4")
  }
  dots.args = list(...)
  # if(length(grep("tdwg_level[1-4]",names(dots.args)))==0L){
  #   cat('failed\n')
  #   stop("The path(s) to the tdwg ESRI shapefiles is (are) missing. Please provide a valid path(s)")
  # }
  # if(recursive & initial_level<4){
  #   if(length(grep("tdwg_level[1-4]",names(dots.args)))<2){
  #     cat('failed\n')
  #     stop("Please provide 2 paths to the tdwg Esri shapefiles when argument recursive is set to TRUE.")
  #   }
  # }
  # provided_levels = as.numeric(stringr::str_extract(grep("tdwg_level[1-4]",names(dots.args),value=TRUE),"[1-4]"))
  # if(!initial_level %in% provided_levels){
  #   cat('failed\n')
  #   stop(paste("Please provide a path to the tdwg shapefile level",initial_level))
  # }

  switch(as.integer(initial_level),
         {
           geographic_level 	<- "continent_code_l1"
           level_code 			<- "LEVEL1_COD"
           # if(!exists('tdwg_level1')){
           #   if(length(dots.args)>0L){
           #     if("tdwg_level1" %in% names(dots.args)){
           #       if(file.exists(dots.args[["tdwg_level1"]])){
           #         tdwg_level1	<<-  rgdal::readOGR(dsn=dirname(dots.args[["tdwg_level1"]]), layer=gsub(".shp","",basename(dots.args[["tdwg_level1"]])),verbose=FALSE)
           #         geographic_polygon <- tdwg_level1
           #       }
           #     }
           #   }else{
           #     fileURL	<- "https://raw.githubusercontent.com/tdwg/wgsrpd/master/geojson/level1.geojson"
           #     tdwg_level1 <<- geojsonio::geojson_read(fileURL, what="sp")
           #     geographic_polygon <- tdwg_level1
           #   }
           # }else
             geographic_polygon <- tdwg_level1
         },
         {
           geographic_level 	<- "region_code_l2"
           level_code 			<- "LEVEL2_COD"
           # if(!exists('tdwg_level2')){
           #   if(length(dots.args)>0L){
           #     if("tdwg_level2" %in% names(dots.args)){
           #       if(file.exists(dots.args[["tdwg_level2"]])){
           #         tdwg_level2	<<-  rgdal::readOGR(dsn=dirname(dots.args[["tdwg_level2"]]), layer=gsub(".shp","",basename(dots.args[["tdwg_level2"]])),verbose=FALSE)
           #         geographic_polygon <- tdwg_level2
           #       }
           #     }
           #   }else{
           #     fileURL	<- "https://raw.githubusercontent.com/tdwg/wgsrpd/master/geojson/level2.geojson"
           #     tdwg_level2 <<- geojsonio::geojson_read(fileURL, what="sp")
           #     geographic_polygon <- tdwg_level2
           #   }
           # }else
             geographic_polygon <- tdwg_level2
         },
         {
           geographic_level 	<- "area_code_l3"
           level_code 			<- "LEVEL3_COD"
           # if(!exists('tdwg_level3')){
           #   if(length(dots.args)>0L){
           #     if("tdwg_level3" %in% names(dots.args)){
           #       if(file.exists(dots.args[["tdwg_level3"]])){
           #         tdwg_level3	<<-  rgdal::readOGR(dsn=dirname(dots.args[["tdwg_level3"]]), layer=gsub(".shp","",basename(dots.args[["tdwg_level3"]])),verbose=FALSE)
           #         geographic_polygon <- tdwg_level3
           #       }
           #     }
           #   }else{
           #     fileURL	<- "https://raw.githubusercontent.com/tdwg/wgsrpd/master/geojson/level3.geojson"
           #     tdwg_level3 <<- geojsonio::geojson_read(fileURL, what="sp")
           #     geographic_polygon <- tdwg_level3
           #   }
           # }else
             geographic_polygon <- tdwg_level3
         },
         {
           geographic_level 	<- "area_code_l3" # not available yet
           level_code 			<- "level4_cod"
           # if(!exists('tdwg_level4')){
           #   if(length(dots.args)>0L){
           #     if("tdwg_level4" %in% names(dots.args)){
           #       if(file.exists(dots.args[["tdwg_level4"]])){
           #         tdwg_level4	<<-  rgdal::readOGR(dsn=dirname(dots.args[["tdwg_level4"]]), layer=gsub(".shp","",basename(dots.args[["tdwg_level4"]])),verbose=FALSE)
           #         geographic_polygon <- tdwg_level4
           #       }
           #     }
           #   }else{
           #     fileURL	<- "https://raw.githubusercontent.com/tdwg/wgsrpd/master/geojson/level4.geojson"
           #     tdwg_level4 <<- geojsonio::geojson_read(fileURL, what="sp")
           #     geographic_polygon <- tdwg_level4
           #   }
           # }else
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
             # if(!exists('tdwg_level2')){
             #   if(length(dots.args)>0L){
             #     if("tdwg_level2" %in% names(dots.args)){
             #       if(file.exists(dots.args[["tdwg_level2"]])){
             #         tdwg_level2	<<-  rgdal::readOGR(dsn=dirname(dots.args[["tdwg_level2"]]), layer=gsub(".shp","",basename(dots.args[["tdwg_level2"]])),verbose=FALSE)
             #         sub_geographic_polygon <- tdwg_level2
             #       }
             #     }
             #   }else{
             #     fileURL	<- "https://raw.githubusercontent.com/tdwg/wgsrpd/master/geojson/level2.geojson"
             #     tdwg_level2 <<- geojsonio::geojson_read(fileURL, what="sp")
             #     sub_geographic_polygon <- tdwg_level2
             #   }
             # }else
               sub_geographic_polygon <- tdwg_level2
           },
           {
             geographic_sub_level 	<- "area_code_l3"
             sub_level_code 			<- "LEVEL3_COD"
             # if(!exists('tdwg_level3')){
             #   if(length(dots.args)>0L){
             #     if("tdwg_level3" %in% names(dots.args)){
             #       if(file.exists(dots.args[["tdwg_level3"]])){
             #         tdwg_level3	<<-  rgdal::readOGR(dsn=dirname(dots.args[["tdwg_level3"]]), layer=gsub(".shp","",basename(dots.args[["tdwg_level3"]])),verbose=FALSE)
             #         sub_geographic_polygon <- tdwg_level3
             #       }
             #     }
             #   }else{
             #     fileURL	<- "https://raw.githubusercontent.com/tdwg/wgsrpd/master/geojson/level3.geojson"
             #     tdwg_level3 <<- geojsonio::geojson_read(fileURL, what="sp")
             #     sub_geographic_polygon <- tdwg_level3
             #   }
             # }else
               sub_geographic_polygon <- tdwg_level3
           },
           {
             geographic_sub_level 	<- "area_code_l3" # not available yet
             sub_level_code 			<- "level4_cod"
             # if(!exists('tdwg_level4')){
             #   if(length(dots.args)>0L){
             #     if("tdwg_level4" %in% names(dots.args)){
             #       if(file.exists(dots.args[["tdwg_level4"]])){
             #         tdwg_level4	<<-  rgdal::readOGR(dsn=dirname(dots.args[["tdwg_level4"]]), layer=gsub(".shp","",basename(dots.args[["tdwg_level4"]])),verbose=FALSE)
             #         sub_geographic_polygon <- tdwg_level4
             #       }
             #     }
             #   }else{
             #     fileURL	<- "https://raw.githubusercontent.com/tdwg/wgsrpd/master/geojson/level4.geojson"
             #     tdwg_level4 <<- geojsonio::geojson_read(fileURL, what="sp")
             #     sub_geographic_polygon <- tdwg_level4
             #   }
             # }else
               sub_geographic_polygon <- tdwg_level4
           }
    )
  }
  # if(!exists('kew_lookup_table',envir=.GlobalEnv)){
  #   if(is.null(lookup_table)){
  #     dots.args = list(...)
  #     if(!"lookup_table" %in% names(dots.args)){
  #       if(verbose) cat('failed\n')
  #       stop(paste0("Unable to locate the species lookup table. Please provide a valid path to Kew lookup table file."))
  #     }else if(file.exists(dot.args[['lookup_table']]))
  #       kew_lookup_table <<- data.table::fread(lookup_table,sep="\t", key=c("db_id"),showProgress=FALSE)
  #   }else if(!file.exists(lookup_table)){
  #     if(verbose) cat('failed\n')
  #     stop(paste0("Unable to locate file: ",lookup_table,". Please provide a valid path to Kew lookup_table file."))
  #   }else kew_lookup_table <<- data.table::fread(lookup_table,sep="\t", key=c("db_id"),showProgress=FALSE)
  # }
  # subset or not according to native/introduced status
  if(!status%in%c('native','introduced','both'))
    stop("Argument 'status' must be one of 'native', 'introduced' or 'both'.")

  geographic_code <- get_geocode(species_name, species_id=species_id, status=status, level=initial_level, use_name_matching=use_name_matching)

  if(do.recursive)
    sub_geographic_code <- get_geocode(species_name, species_id=species_id, status=status, level=sub_level, use_name_matching=use_name_matching)

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
      cat('failed\n')
      warning(paste0("Unable to find geo-data for species '",species_name,"' with the status '",status,"'"))
    }
    return(list(species=species_name, status=status, levels=c(initial_level, if(do.recursive) initial_level+1), point_fraction=NULL, unit_fraction=NULL, range_filling=NULL, crs=NULL, cleaned_point_data = NULL))
  }
  geographic_code = unique(geographic_code)
  if(do.recursive){
    sub_geographic_code	= unique(sub_geographic_code)
  }
  if(verbose) cat("done\n\n")
  if(verbose){
    cat('#=========================================================\n')
    cat('#= 2. Get species polygons\n')
    cat('#=========================================================\n')
  }
  if(verbose) cat("> 2. getting species polygons ...")
  tdwg_polygons <- eval(parse(text=sprintf("subset(geographic_polygon, %s %%in%% geographic_code)",level_code)))
  stopifnot(inherits(tdwg_polygons,"SpatialPolygonsDataFrame"))
  if(do.recursive)
    tdwg_sub_level_polygons <- eval(parse(text=sprintf("subset(sub_geographic_polygon, %s %%in%% sub_geographic_code)",sub_level_code)))
  if(verbose) cat("done\n\n")
  if(!is.null(which_skip)){
    if(!all(inherits(which_skip,'numeric')) & !all(which_skip %in% c(2,3)) & length(which_skip)<=2)
      stop("Argument 'which.skip' must be an integer value or vector with value(s) 1, 2 or 3.")
  }
  if(verbose){
    if(!is.null(which_skip) & all(c(1,2,3) %in% which_skip) ){
      cat("> skipping computation of points distribution parameters...\n")
    }else{
      cat('#====================================================================\n')
      cat('#= 3. Compute points distribution parameters\n')
      cat('#====================================================================\n')
    }
  }

  #-----------------------------------
  #=3.a fraction of occurrence points
  #-----------------------------------

  point_data		<- setup_point_data(point_data, full_data=full_data)

  if(nrow(point_data)==0){
    if(verbose)
      warning(paste("No points available after quick setup for species:",species_name))
    initial_level_data <- as.data.frame(array(NA, dim=c(if(by_id & length(tdwg_polygons)>1) length(tdwg_polygons) else 1, 4)))
    colnames(initial_level_data) <- c("frac_occ","frac_units","perc_cover","region")
    out = list(species=species_name, status=status, initial_level_data, crs= crs(tdwg_polygons), cleaned_point_data = NULL)
    names(out)[3] <- paste0("level",initial_level)

    if(do.recursive){
      sub_level_data <- as.data.frame(array(0, dim=c(length(tdwg_sub_level_polygons), 4)))
      colnames(sub_level_data) <- c("frac_occ","frac_units","perc_cover","region")
      out = append(out, list(sub_level_data), after=3)
      names(out)[4] <- paste0("level",sub_level)
    }
    return(out)
  }
  nb_pts_total 	<- nrow(point_data)
  nb_tdwg_units <- length(tdwg_polygons)
  id_x_lon 	<- grep(pattern = "[Ll][Oo][Nn]",x = names(point_data))[1]
  id_y_lat 	<- grep(pattern = "[Ll][Aa][Tt]",x = names(point_data))[1]
  coordHeaders <- c(id_x_lon, id_y_lat)
  sp.pts			  <- sp::SpatialPoints(point_data[,coordHeaders], crs(tdwg_polygons))
  is_in_range  	<- complete.cases(sp::over(sp.pts, tdwg_polygons[level_code]))

  # create empty array to store points distribution parameters at initial level
  initial_level_data <- as.data.frame(array(0, dim=c(if(by_id & nb_tdwg_units>1) nb_tdwg_units else 1, 4)))
  colnames(initial_level_data) <- c("frac_occ","frac_units","perc_cover","region")

  do.frac_occ = !(!is.null(which_skip) & 1 %in% which_skip)
  if(verbose){
    if(!do.frac_occ)
      cat("> skipping computation of the fraction of occurrences lying within TDWG units...\n")
  }
  if(do.frac_occ){
    if(verbose) cat("> computing the fraction of occurrences lying within TDWG units...\n")

    point_fraction <- sum(is_in_range)/nb_pts_total
    if(point_fraction==0){
      if(verbose) cat("No points has been found within the known range. Exiting...\n")
      out = list(species=species_name, status=status, eval(parse(text=sprintf("level%g=initial_level_data",initial_level))), crs= crs(tdwg_polygons), cleaned_point_data = NULL)
      names(out)[3] <- paste0("level",initial_level)

      if(do.recursive){
        sub_level_data <- as.data.frame(array(0, dim=c(length(tdwg_sub_level_polygons), 4)))
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
        for(k in tdwg_polygons[level_code][[1]]){
          nb_pts_in_sub_range <- sum(complete.cases(sp::over(sp.pts, subset(tdwg_polygons[level_code],eval(parse(text=paste0(level_code,"==",k)))))))
          initial_level_data[match(k,tdwg_polygons[level_code][[1]]),c(1,4)] <- c(nb_pts_in_sub_range/nb_pts_total,k)
        }
      }
      else{
        initial_level_data[1,1]	<- point_fraction
        initial_level_data[1,4] <- "Total"
      }
    }

    if(do.recursive){
      nb_tdwg_sub_units 	<- length(tdwg_sub_level_polygons)
      # create empty array to store points distribution parameters at sub level
      sub_level_data <- as.data.frame(array(0, dim=c(nb_tdwg_sub_units, 4)))
      colnames(sub_level_data) <- c("frac_occ","frac_units","perc_cover","region")

      is_in_sub_units    <- vector(length=nb_tdwg_sub_units)
      for(k in tdwg_sub_level_polygons[sub_level_code][[1]]){
        nb_pts_in_sub_units <- sum(complete.cases(sp::over(sp.pts, subset(tdwg_sub_level_polygons[sub_level_code],eval(parse(text=paste0(sub_level_code,"=='",k,"'")))))))
        sub_level_data[match(k,tdwg_sub_level_polygons[sub_level_code][[1]]),c(1,4)] <- c(nb_pts_in_sub_units/nb_pts_total, k)
      }
    }
    if(verbose) cat("done\n")
  }

  #-----------------------------------
  #=3.b proportion of tdwg units
  #-----------------------------------
  do.frac_units = !(!is.null(which_skip) & 2 %in% which_skip)
  if(verbose){
    if(!do.frac_units)
      cat("> skipping computation of the proportion of the tdwg units sampled ...\n")
  }
  if(do.frac_units){
    if(verbose) cat("> computing the proportion of TDWG units in which the species is found...")

    is_in_units <- vector(length=nb_tdwg_units)
    for(k in paste(tdwg_polygons[level_code][[1]]))
      is_in_units[match(k,tdwg_polygons[level_code][[1]])] <- any(complete.cases(sp::over(sp.pts, subset(tdwg_polygons[level_code],eval(parse(text=paste0(level_code,"==",k)))))))

    if(by_id)
      initial_level_data[,2] <- rep(sum(is_in_units)/nb_tdwg_units, times=nb_tdwg_units)
    else
      initial_level_data[1,2] <- sum(is_in_units)/nb_tdwg_units

    if(do.recursive){
      nb_tdwg_sub_units <- length(tdwg_sub_level_polygons)
      is_in_sub_units   <- vector(length=nb_tdwg_sub_units)
      for(k in tdwg_sub_level_polygons[sub_level_code][[1]]){
        is_in_sub_units[match(k,tdwg_sub_level_polygons[sub_level_code][[1]])] <- any(complete.cases(sp::over(sp.pts, subset(tdwg_sub_level_polygons[sub_level_code],eval(parse(text=paste0(sub_level_code,"=='",k,"'")))))))
      }
      sub_level_data[,2] <- sum(is_in_sub_units)/nb_tdwg_sub_units
    }

    if(verbose) cat("done\n")
  }

  #-----------------------------------
  #=3.c percent coverage of tdwg units
  #-----------------------------------
  do.perc_cover = !(!is.null(which_skip) & 3 %in% which_skip)
  if(verbose){
    if(!do.perc_cover)
      cat("> skipping computation of the percent coverage of the tdwg units ...\n")
  }
  if(do.perc_cover){
    if(verbose) cat("> computing the percent coverage of the tdwg units by occurrence points...")

    # keep point inside the tdwg
    sp.pts_in_range	<- sp.pts[is_in_range] # should have >0 pts since the fraction of point inside the tdwg unit >0
    # transform in utm
    sp.pts_utm_in_range <- sp::spTransform(sp.pts_in_range, CRS(proj_utm))
    #require(geosphere)
    # create raster template
    r_init	<- if(!exists('r_init')) { r_init <<- init(raster());r_init} else r_init
    r_temp_fn 		<- file.path(dirname(rasterTmpFile()),paste0(grid_resol,"_",gsub(" ","_",proj_utm),".rds"))
    if(file.exists(r_temp_fn))
      r_temp <- readRDS(r_temp_fn)
    else{
      r_temp  <- projectRaster(r_init, res=grid_resol, crs=CRS(proj_utm))
      saveRDS(r_temp,file=r_temp_fn )
    }
    #sp_patch <- rasterToPolygons(rasterize(sp.pts_utm_in_range, r_temp, fun='count'), fun=function(x) x>0)
    # # create raster template
    # if(length(sp.pts_in_range)>1){
    # bottomLeft.coo = raster::extent(sp.pts_in_range)[c(1,3)];bottomRight.coo = raster::extent(sp.pts_in_range)[c(2,3)];topLeft.coo = raster::extent(sp.pts_in_range)[c(1,4)]# ; topRight = extent(sp.pts_in_range)[c(2,4)] # coordinates of the bounding box corners
    # xdist = geosphere::distGeo(p1=bottomLeft.coo, p2=bottomRight.coo)# longitude distance in m
    # ydist = geosphere::distGeo(p1=bottomLeft.coo, p2=topLeft.coo)# latitude distance in m
    # }
    # else{
    # xdist=6378137
    # ydist=6378137
    # }
    # if(!inherits(grid_resol,'numeric') || grid_resol<=0)
    # stop("Argument 'grid_resol' must be an integer >0")
    # nc = ceiling(xdist/grid_resol) # expected number of columns
    # nr = ceiling(ydist/grid_resol) # expected number of rows
    # rtmp = raster::raster(nrows=nr, ncols=nc, crs=crs(sp.pts_in_range))
    # r <- setExtent(rtmp, ext=if(length(sp.pts_in_range)>1) extend(extent(sp.pts_in_range),2) else extent(-180, 180, -90, 90))
    # raster::res(r) = min(raster::res(r))
    # print(res(r))
    # compute the approximate surface area of cells
    #r_area <- raster::area(r_temp)
    # extract the surface area of cells of occurrence points (in km2)
    # sp.pts_area_cell <- r_area[unique(extract(r_area, sp.pts_utm_in_range, cellnumber=TRUE)[,'cells'])]

    # compute the percent coverage of the tdwg area occupied by occurrence points
    tdwg_utm_polygons_fn <- file.path(dirname(rasterTmpFile()),paste0('tdwg_utm_polygons_',grid_resol,"_",gsub(" ","_",proj_utm),".rds"))
    if(file.exists(tdwg_utm_polygons_fn))
      tdwg_utm_polygons <- readRDS(tdwg_utm_polygons_fn)
    else{
      tdwg_utm_polygons	<- sp::spTransform(tdwg_polygons, CRS(proj_utm))
      saveRDS(tdwg_utm_polygons,file=tdwg_utm_polygons_fn )
    }
    cell_area			<- (xres(r_temp)*yres(r_temp))
    tdwg_units_area		<- if(nb_tdwg_units>1) lengths(raster::cellFromPolygon(r_temp,tdwg_utm_polygons))*cell_area else length(raster::cellFromPolygon(r_temp,tdwg_utm_polygons))*cell_area
    sp.pts_area_cell 	<- rep(cell_area, times=length(unique(raster::extract(r_temp, sp.pts_utm_in_range, cellnumber=TRUE)[,'cells'])))

    if(!by_id){
      initial_level_data[1,3]	<- (sum(sp.pts_area_cell,na.rm=TRUE)/sum(tdwg_units_area))*100
    }
    else if(by_id){
      if(nb_tdwg_units>1){
        range_filling <- vector(length=nb_tdwg_units)
        for(k in tdwg_utm_polygons[level_code][[1]]){
          is_in_id_units   <- complete.cases(sp::over(sp.pts_utm_in_range, subset(tdwg_utm_polygons[level_code],eval(parse(text=paste0(level_code,"==",k))))))
          if(any(is_in_id_units))
            initial_level_data[match(k,tdwg_utm_polygons[level_code][[1]]),3] <- (sum(sp.pts_area_cell[is_in_id_units],na.rm=TRUE)/tdwg_units_area[match(k,tdwg_utm_polygons[level_code][[1]])])*100
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
        tdwg_sub_level_utm_polygons	<- sp::spTransform(tdwg_sub_level_polygons, CRS(proj_utm))
        saveRDS(tdwg_sub_level_utm_polygons,file=tdwg_sub_level_utm_polygons_fn)
      }
      tdwg_sub_units_area <- if(nb_tdwg_sub_units > 1) lengths(raster::cellFromPolygon(r_temp,tdwg_sub_level_utm_polygons))*cell_area else length(raster::cellFromPolygon(r_temp,tdwg_sub_level_utm_polygons))*cell_area #raster::area(tdwg_sub_level_utm_polygons)/1e06
      sub_range_filling 	<- vector(length=nb_tdwg_sub_units)
      for(k in tdwg_sub_level_utm_polygons[sub_level_code][[1]]){
        is_in_id_sub_units   <- complete.cases(sp::over(sp.pts_utm_in_range, subset(tdwg_sub_level_utm_polygons[sub_level_code],eval(parse(text=paste0(sub_level_code,"=='",k,"'"))))))
        if(any(is_in_id_sub_units))
          sub_level_data[match(k,tdwg_sub_level_utm_polygons[sub_level_code][[1]]),3] <- (sum(sp.pts_area_cell[is_in_id_sub_units],na.rm=TRUE)/tdwg_sub_units_area[match(k,tdwg_sub_level_utm_polygons[sub_level_code][[1]])])*100
      }
    }

    if(verbose) cat("done\n\n")
  }
  # returns a SpatialPointsDataFrame
  # if(sp){
  # return(
  # SpatialPointsDataFrame(	coords=if(point_fraction[1]==0) data.frame(lon=0,lat=0) else sp.pts[is_in_range]@coords,
  # data=data.frame(species=species_name,
  # status=status,
  # initial_level=initial_level,
  # sub_level =if(do.recursive) sub_level else NA,
  # point_fraction=point_fraction,
  # unit_fraction=unit_fraction,
  # range_filling=if(!by_id) range_filling else NA
  # ),
  # proj4string=crs(sp.pts)
  # )
  # )
  # }
  # return list
  out =list(species=species_name, status=status, initial_level_data, crs=crs(sp.pts), cleaned_point_data=if(sp) sp.pts[is_in_range] else if(full_data) point_data[is_in_range,] else sp.pts[is_in_range]@coords)
  names(out)[3] <- paste0("level",initial_level)

  if(do.recursive){
    out = append(out, list(sub_level_data), after=3)
    names(out)[4] <- paste0("level",sub_level)
  }

  # save info
  saveRDS(list(call.tmp,out), file=tdwg_info_fn)

  return(out)
}

