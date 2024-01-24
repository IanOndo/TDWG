#' Functions to filter species occurrences or select species based on a system of geographical units developped by The International Working Group on Taxonomic Databases for Plant Sciences (TDWG) at approximately "country" level and upwards.
#'
#' Filter or select species occurrences or species not occurring within a given range.
#'
#' @param point_data A two-column data.frame of occurrence records coordinates, a path to a directory or a file containing the occurrence records coordinates of the species.
#' @param what A character string specifying whether the filtering must be applied on species or their occurrences. Can be either \code{'species'} or \code{'occurrences'}
#' @param working_dir A character string specifying the path to the output directory.
#' @param species_name A two character string (Genus species) specifying the name of the species to check. Ignored if point_data is a directory or if \code{fileAsNames=TRUE}.
#' @param species_id A character string specifying the id of the species in ipni or Kew world checklist database.
#' @param status A character string specifying the geographic range where to search for the species. Can be : 'native', 'introduced' or 'both'
#' @param initial_level A numeric integer between 1 and 4. The initial unit level of the TWDG at which the data must be checked against. Default initial level is 2.
#' @param by_id A logical. Should the criteria be applied to each geographical units ? Default is FALSE. Ignored if \code{what='occurrences'}
#' @param recursive A logical. Should the criteria be also applied at the subsequent lower units level ? Default is FALSE. Ignored if \code{what='occurrences'}
#' @param force.output A logical. Should the initial point dataset be returned in case where all points have been removed by the cleaning ? Default is FALSE.
#' @param sf A logical. Should the results be returned as a sf object ? Default is FALSE. Ignored if \code{cleanOcc=FALSE}
#' @param fileAsNames A logical. Should species name be retrieved from file names ? Default is FALSE. Ignored if point_data is a path to a directory.
#' @param do.parallel A logical. Should computations be run in parallel. Default is FALSE.
#' @param ncores A numeric integer specifying the number of cores to use in parallel processing. Ignored if \code{do.parallel=FALSE}.
#' @param save.outputs A logical. Should the results be saved on the disk. Default is TRUE. The results wil be saved in the directory given in \code{working_dir} argument.
#'                                Ignored if \code{what='species'}.
#' @param ... Additional parameters.
#'            If \code{what='occurrences'} ...is ignored.
#'            If \code{what='species'}, ... should include parameters for \code{speciesRangeCleaner()}, such as:
#'            the point_fraction which is the minimum fraction of occurrence records that must be included within the geographic range defined in \code{status};
#'            the unit_fraction which is the minimum fraction of geographical units that must include occurrence records;
#'            the range_filling which is the minimum percent of the geographic area covered by the occurrence records distribution;
#'            or the grid_resol which specifies the resolution of the raster grid to be used when calculating the percentage of geographic range covered by occurrence records.
#'            Or other additional parameters to be passed to \code{TDWGinfo()}
#'            See \code{TDWGinfo()} documentation or examples below for more details.
#' @return If \code{what='species'}, the function returns a \code{TRUE} if the distribution of species occurrences follows the range specified in \code{status}, and the given criteria. \code{FALSE} otherwise.
#'         If \code{what='occurrences'}, the function returns a list whose elements are the species occurrence records cleaned by the range specified in \code{status}.
#' @author <i.ondo@kew.org>
#' @seealso \code{TDWGinfo()}, \code{speciesRangeCleaner()}
#' @examples wk_dir = <path to a csv species files directory >
#'			out = rangeCleaner(point_data=wk_dir,
#'				working_dir=wk_dir,
#'				what='occurrences',
#'				status='native',
#'				do.parallel=FALSE)
#' @export
rangeCleaner <- function(point_data,
                         what = "species",
                         working_dir = NULL,
                         species_name = NULL,
                         species_id = NULL,
                         col_species = NULL,
                         status = 'both',
                         initial_level = 2,
                         backbone='wcvp',
                         force.output = FALSE,
                         sf = FALSE,
                         fileAsNames = FALSE,
                         do.parallel = FALSE,
                         ncores = NULL,
                         save.outputs = TRUE,
                         verbose = TRUE,...) {

  if(verbose){
    cat('#====================#\n')
    cat('#=1. Set user input  #\n')
    cat('#====================#\n')
  }
  work_dir_flag = if(is.null(working_dir)) FALSE else tryCatch(dir.exists(working_dir), error=function(err) FALSE)
  if(missing(point_data))
    stop("Argument 'point_data' is missing")
  dir_flag = tryCatch(dir.exists(point_data), error=function(err) FALSE)
  file_flag= tryCatch(file.exists(point_data), error=function(err) FALSE) && !tryCatch(dir.exists(point_data), error=function(err) FALSE)
  data_flag = inherits(point_data,c("data.frame","data.table"))

  if(verbose){
    cat('#---------------------------------\n')
    cat('#=1.a Make a list of all species\n')
    cat('#---------------------------------\n')
  }
  if(!work_dir_flag){
    if(verbose)
      warning(paste0("Unable to find working directory: '",working_dir,"'. Setting working directory to current one: ",getwd()));flush.console()
    working_dir = getwd()
  }
  if(dir_flag){

    dn <- dir(working_dir, pattern="^Cleaned_data_by_", full.names=TRUE)
    if(length(dn)>0L)
      unlink(dn, recursive=TRUE)


    list.species <- list.files(point_data, pattern="\\.csv$", full.names=TRUE, recursive=TRUE) # select species
    if(length(list.species)==0L)
      stop(paste0("Unable to find csv files in directory:",point_data,". Please provide a directory with csv files."))
    empty.string <- nchar(list.species)==0L
    if(any(empty.string)){
      if(verbose)
        warning(paste("Removing",sum(empty.string),"file(s) with no species names."))
      list.species <- list.species[!empty.string]
    }
    if(!is.null(species_name)){
      if(length(species_name)!=length(list.species))
        stop(paste("If provided, argument 'species_name' must have the same length as the number of files in",point_data))
      diff.list.species = setdiff(list.species,species_name)
      if(length(diff.list.species)>0L)
        if(verbose) warning(paste(length(diff.list.species),"species from 'species_name' was/were not found in",point_data))
      list.species = list.species[list.species%in%species_name]
    }
  }else if(file_flag){
    if(!grepl("\\.csv$", point_data))
      stop(paste0(basename(point_data),"is not a csv file. Please provide a csv file."))
    if(fileAsNames & !is.null(species_name))
      stop("Arguments 'species_name' and 'fileAsNames' are mutually exclusive. You must either set 'fileAsNames' to FALSE or 'species_name' to NULL")
    if(!is.null(species_name) & length(species_name)!=length(point_data))
      stop("Arguments 'point_data' and 'species_name' must have the same length.")
    if(is.null(species_name) & !fileAsNames){
      occ_data 		<- data.table::fread(point_data, showProgress=FALSE)
      id_sp_name		<- grep(pattern = "^[Ss][Pp]|[Bb][Ii][Nn][Oo][Mm]",x = names(occ_data))[1]
      if(length(id_sp_name)==0L)
        stop("Unable to find a species name for the data. Please specified a name for the species with 'species_name' or using file name (i.e. fileAsNames set to TRUE)")
      sp.colname <- colnames(occ_data)[id_sp_name]
      setkeyv(occ_data, sp.colname)
      list.species 	<- c(as.character(unique(occ_data[[sp.colname]]))) # select species
    }
    else
      list.species <- gsub("\\.csv","",basename(point_data))
  }
  else if(data_flag){
    id_sp_name	<- ifelse(!is.null(col_species) && (col_species > 0 & col_species <= ncol(point_data)), col_species, grep(pattern = "^[Ss][Pp]|[Bb][Ii][Nn][Oo][Mm]",x = names(point_data))[1])
    if(length(id_sp_name)==0L)
      stop("Unable to find a species name for the species name column. Please add a valid species name column to your data.")
    sp.colname <- colnames(point_data)[id_sp_name]
    if(!inherits(point_data,"data.table")) data.table::setDT(point_data)
    if(!data.table::haskey(point_data)) data.table::setkeyv(point_data, sp.colname)
    list.species 	<- c(as.character(unique(point_data[[sp.colname]])))
    if(!is.null(species_name) && length(species_name)!=length(list.species))
      stop("Arguments 'species_name' must have the same length as the number of species detected in argument 'point_data'.")
    if(!is.null(species_name) && species_name != list.species){
      if(verbose)
        warning(paste0("Argument 'species_name' and the species name detected in point_data are different: ",species_name," vs ", list.species,". Overriding detected name..."))
    }
  }
  else{
    stop(paste0("Unable to find directory or file '",point_data,"'. Please provide a valid path"))
  }

  use_sp_id = FALSE
  if(!is.null(species_id)){
    if(length(list.species) != length(species_id))
      stop("Argument 'species id' must have the same length as the number of species to be cleaned.")
    if(length(names(species_id))==0L)
      stop("Argument 'species_id' must be a named vector")
    use_sp_id = TRUE
  }

  if(verbose){
    cat("#---------------------------------\n")
    cat("#=1.b Set up environment\n")
    cat("#---------------------------------\n")
  }
  call.fun 	<- match.call(expand.dots=TRUE)
  tmp.args    <- c('species_name','species_id','status','point_fraction', 'unit_fraction','range_filling','grid_resol','initial_level','backbone','by_id','which_skip','recursive','force.output','verbose','use_name_matching','proj_utm','full_data')
  call.tmp    <- call.fun[match(tmp.args, names(call.fun),nomatch=0)]
  cleaner		<- switch(what, 'species' = "speciesRangeCleaner", 'occurrences'="coordinatesRangeCleaner")
  cleaner.args<- as.list(call.tmp)

  if(dir_flag | file_flag){
    arg_data = if(file_flag) dirname(point_data) else point_data
    if(arg_data==working_dir)
      working_dir <- file.path(working_dir, paste0("Cleaned_data_by_",if(cleaner.args$status=="both") "known" else cleaner.args$status,"_range"))
  }
  # if(!exists('kew_checklist',envir=.GlobalEnv)){
  #   if(is.null(checklist) || !file.exists(checklist))
  #     stop(paste0("Unable to locate file: ",checklist,". Please provide a valid path to Kew checklist file."))
  #   kew_checklist <<- data.table::fread(checklist, key=c("acc_full_name_without_family"), showProgress=FALSE)
  # }
  #
  # if(!exists('kew_lookup_table',envir=.GlobalEnv)){
  #   if(is.null(lookup_table) || !file.exists(lookup_table))
  #     stop(paste0("Unable to locate file: ",lookup_table,". Please provide a valid path to Kew lookup table file."))
  #   kew_lookup_table <<- data.table::fread(lookup_table,sep="\t",key=c("db_id"), showProgress=FALSE)
  # }

  if(verbose){
    cat("#---------------------------------\n")
    cat("#=1.c Set up computing\n")
    cat("#---------------------------------\n")
  }
  if(do.parallel & length(list.species)>1){
    # Parallel processing
    toExport <- c("cleaner", "cleaner.args","working_dir","use_sp_id")# send objects to cluster nodes
    # toExport <- c("cleaner", "cleaner.args", "setup_point_data","hasDistrib","TDWGinfo","working_dir")# send objects to cluster nodes
    # switch(status,
    #        'native'= {
    #          toExport <- append(toExport,"hasNativeDistrib")
    #        },
    #        'introduced'= {
    #          toExport <- append(toExport,"hasIntroducedDistrib")
    #        },
    #        'both' = {
    #          toExport <- append(toExport,"hasKnownDistrib")
    #        }
    # )
    # switch(cleaner,
    #        'speciesRangeCleaner' = {
    #          toExport <- append(toExport,"speciesRangeCleaner")
    #        },
    #
    #        'coordinatesRangeCleaner'= {
    #          toExport <- append(toExport,"coordinatesRangeCleaner")
    #          toExport <- append(toExport,switch(status, 'native'='cleanByNativeDistrib', 'introduced'='cleanByIntroducedDistrib', 'both' ='cleanByKnownDistrib'))
    #        }
    # )
    if(file_flag)
      toExport <- append(toExport,c("occ_data"))
    if(data_flag)
      toExport <- append(toExport,c("point_data"))

    ncores = if(is.null(ncores)) parallel::detectCores()-1 else min(ncores,parallel::detectCores()) # number of cores to use
    doParallel::registerDoParallel(ncores)
  }else{
    toExport <- NULL
    # Sequential processing
    foreach::registerDoSEQ()
  }

  if(verbose){
    cat("#========================#\n")
    cat("#=2. Process\n")
    cat("#========================#\n")
    cat(paste0('#',paste(rep('-',times=100),collapse="")))
    cat('\n')
    cat(paste0('> Start Run on: ',Sys.time()),'\n\n')
  }

  tryCatch({

    require(doParallel)

    started.at <- Sys.time()

      out <- foreach::foreach(k = list.species, .packages = c("TDWG","geosphere", "rgdal", "geojsonio","stringr", "raster", "parallel", "foreach","data.table"), .export=toExport, .errorhandling = 'pass') %dopar% {

        # Get data for just that species
        if(file.exists(k)){
          species.data <- tryCatch(data.table::fread(k, showProgress=FALSE),error=function(err) return(NULL))
          k = stringr::str_replace_all(basename(k),c(".csv"="","_"=" "))
        }else if(exists('occ_data')){
          species.data 	<- tryCatch({
            occ_data[k]
          },
          error = function(err) {
            tryCatch(occ_data[J(k)],error = function(err) return(NULL))
          })
        }else{
          species.data 	<- tryCatch({
            point_data[k]
          },
          error = function(err) {
            tryCatch(point_data[J(k)],error = function(err) return(NULL))
          })
        }
        # if an error occurred during the reading returns an NULL
        if (is.null(species.data) || ncol(species.data) < 2){
          if(verbose) warning(paste0("Cannot read file from species ",k))
          if(save.outputs & cleaner == "coordinatesRangeCleaner"){
            dn <- file.path(working_dir,"problems")
            if(!dir.exists(dn))
              dir.create(dn, recursive =TRUE)
            fn <-file.path(dn,'Cannot_read_file.csv')
            if(!file.exists(fn))
              write.table(data.frame(Species=k), file=fn, row.names=FALSE, col.names = !file.exists(fn), sep=",",  append=TRUE)
            else if(!k %in% read.csv(fn)$Species)
              write.table(data.frame(Species=k), file=fn, row.names=FALSE, col.names = !file.exists(fn), sep=",",  append=TRUE)
          }
          return(NULL)
        }
        # if the species has no points
        if(nrow(species.data)< 1){
          if(verbose) warning(paste0("Species '",k,"' has no points"))
          if(save.outputs & cleaner == "coordinatesRangeCleaner"){
            dn <- file.path(working_dir,"problems")
            if(!dir.exists(dn))
              dir.create(dn, recursive =TRUE)
            fn <-file.path(dn,'Has_no_points.csv')
            if(!file.exists(fn))
              write.table(data.frame(Species=k), file=fn, row.names=FALSE, col.names = !file.exists(fn), sep=",",  append=TRUE)
            else if(!k %in% read.csv(fn)$Species)
              write.table(data.frame(Species=k), file=fn, row.names=FALSE, col.names = !file.exists(fn), sep=",",  append=TRUE)
          }
          return(NULL)
        }
        # Print species name
        if(verbose) cat(sprintf("Species name: %s\n\n", k));flush.console()
        cleaner.args[['point_data']] 	<- species.data
        cleaner.args[['species_name']]	<- k
        if(use_sp_id){
          species_id <- eval(cleaner.args$species_id,parent.frame())
          if(length(species_id)>1L)
            species_id <- species_id[k]
          cleaner.args[['species_id']]	<- species_id
        }
        #if(is.null(cleaner.args[['species_name']]))  cleaner.args[['species_name']]	<- k
        res = do.call(cleaner, cleaner.args)
        if(cleaner=="coordinatesRangeCleaner"){
          if(is.null(res) || nrow(res)==0L) return(NULL)
          if(save.outputs){
            dn <- working_dir
            if(!dir.exists(dn))
              dir.create(dn, recursive =TRUE)
            write.csv(res, file.path(dn,paste0(k,'.csv')), row.names=FALSE)
          }
        }
        return(res)
      }
      finished.at <- Sys.time()
      time.elapsed <- finished.at - started.at

  }, finally={
    if(do.parallel & length(list.species)>1)
      doParallel::stopImplicitCluster()
  })

  if(verbose){
    cat('...End of computations...\n\n')
    cat(paste0('> End of Run on: ',Sys.time()),'\n\n')
    cat(paste0('> Running time: ',as.numeric(time.elapsed),' ', attr(time.elapsed,'units'),'\n'))
    cat(paste0('#',paste(rep('-',times=100),collapse="")))
    cat('\n')
  }

  if(cleaner=="speciesRangeCleaner"){
    out <- do.call('c', out)
    out <- setNames(out,gsub("\\.csv","",basename(list.species)))
    return(out)
  }
  out <- setNames(out,gsub("\\.csv","",basename(list.species)))
  return(out)
}
