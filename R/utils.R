#' @export
setup_point_data <- function(point_data, full_data=FALSE){
  #--------------------------------------
  #= 1. Retrieve coordinates
  #--------------------------------------
  id_x_lon 	<- grep(pattern = "[Ll][Oo][Nn]",x = names(point_data))[1]
  id_y_lat 	<- grep(pattern = "[Ll][Aa][Tt]",x = names(point_data))[1]
  coordHeaders <- c(id_x_lon, id_y_lat)

  #--------------------------------------
  #= 2. Make sure coordinates are numeric
  #--------------------------------------
  if(inherits(point_data,"data.table")) setDF(point_data)
  point_data[, coordHeaders][, !sapply(point_data[, coordHeaders], is.numeric)] <-sapply(point_data[, coordHeaders][, !sapply(point_data[, coordHeaders], is.numeric)], function(f)  as.numeric(levels(f))[f])

  #--------------------------------------
  #= 3. Clean coordinates
  #--------------------------------------
  NA_to_remove <- !complete.cases(point_data[, coordHeaders]) # find NA's
  Dup_to_remove <- duplicated(point_data[, coordHeaders]) # find duplicates
  Rows_to_remove <- NA_to_remove | Dup_to_remove
  point_data.cleaned	<- point_data[!Rows_to_remove, ] # remove Na's and duplicates

  if(full_data)
    return(point_data.cleaned)

  return(point_data.cleaned[, coordHeaders])
}

#' @export
hasDistrib <- function(species_name, checklist=NULL, use_name_matching=FALSE, verbose=TRUE) sapply(species_name, .hasdistrib, checklist=checklist, use_name_matching=use_name_matching, verbose=verbose, USE.NAMES=TRUE)
.hasdistrib <- function(species_name, checklist=NULL, use_name_matching=FALSE, verbose=TRUE){

  if(missing(species_name))
    stop("Argument 'species_name' is missing")
  if(!inherits(species_name,"character"))
    stop("Argument 'species_name' must be a character string.")
  if(nchar(species_name)==0L)
    stop("Invalid species name")

  #----------------------------------------------
  #= Get distribution info from the checklist
  #----------------------------------------------
  if(!exists('kew_checklist',envir=.GlobalEnv)){
    if(is.null(checklist) || !file.exists(checklist))
      stop(paste0("Unable to locate file: ",checklist,". Please provide a valid path to Kew checklist file."))
    kew_checklist <<- data.table::fread(checklist, key=if(use_name_matching) "full_name_without_family" else "acc_full_name_without_family",showProgress=FALSE)
  }else{
    if(use_name_matching && key(kew_checklist)!="full_name_without_family") setkey(kew_checklist,"full_name_without_family")
    if(!use_name_matching && key(kew_checklist)!="acc_full_name_without_family") setkey(kew_checklist,"acc_full_name_without_family")
  }

  distrib <- try(kew_checklist[species_name][['hasdistribution']])
  if(inherits(distrib,"error")){
    if(verbose)
      warning(paste("Unable to find species",species_name,"in the TDWG database. Returning NA"));flush.console()
    return(NA)
  }
  return(any(distrib=="Y"))
}
