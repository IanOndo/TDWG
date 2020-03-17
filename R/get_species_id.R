#' Functions to get the species IPNI ID from the Kew checklist.
#'
#' Get the IPNI id of a species that have geographic distribution reported.
#'
#' @param species_name A two character string (Genus species) specifying the name of the species to search for.
#' @param status Optional. A character string specifying the type of geographic region where to search for the species. Can be : 'native', 'introduced' or 'both'
#' @param exclude Optional. A character string or a vector of names specifying the geographic region(s) to exclude from the search.
#'                          Can be : "NORTHERN AMERICA","SOUTHERN AMERICA", "ANTARCTICA", "EUROPE", "AFRICA", "AUSTRALASIA","ASIA-TEMPERATE", "ASIA-TROPICAL" or "PACIFIC"
#' @export
get_species_id <- function(species_name, status='native', level=2, exclude=NULL){

  if(key(kew_checklist)!="full_name_without_family") setkey(kew_checklist,'full_name_without_family')
  best_guess 	<- unique(kew_checklist[species_name][["acc_full_name_without_family"]])
  geo_level 	<- switch(level, "continent_code_l1", "region_code_l2", "area_code_l3", "area_code_l3")
  fun 		<- function(x) all(!x$continent %in% exclude)

  if(is.na(best_guess)){	# try last updates from api ipni
    results <- queryDistributionId(species_name,status=ifelse(status=='native','natives',status))
    if(is.null(results)) return(NULL)
    if(key(kew_checklist)!='ipni_id') data.table::setkey(kew_checklist,'ipni_id')
    db_id 		<- kew_checklist[results][['db_id']]
    geo_data 	<- kew_lookup_table[, fun(.SD), by = db_id][V1==TRUE]
    if(length(geo_data[['db_id']])==0L || all(is.na(geo_data[['db_id']]))) return(NULL)
    if(!haskey(geo_data)) data.table::setkey(geo_data, 'db_id')
    ipni_id 	<- results[keyvalExists(geo_data, db_id)]
    if(length(ipni_id)==0L) return(NULL)
    return(ipni_id)
  }
  if(length(best_guess)>1){
    db_ids <- ifelse(best_guess=='Unplaced Unplaced', kew_checklist[species_name][["db_id"]], kew_checklist[species_name][["accepted_db_id"]])
    has_distrib	 <- which(keyvalExists(kew_lookup_table, db_ids, how=NULL))
    if(length(has_distrib)==0L){	# try last updates from api ipni
      results <- queryDistributionId(species_name,status=ifelse(status=='native','natives',status))
      if(is.null(results)) return(NULL)
      if(key(kew_checklist)!='ipni_id') data.table::setkey(kew_checklist,'ipni_id')
      db_id 		<- kew_checklist[results][['db_id']]
      geo_data 	<- kew_lookup_table[, fun(.SD), by = db_id][V1==TRUE]
      if(length(geo_data[['db_id']])==0L || all(is.na(geo_data[['db_id']]))) return(NULL)
      if(!haskey(geo_data)) data.table::setkey(geo_data, 'db_id')
      ipni_id 	<- results[keyvalExists(geo_data, db_id)]
      if(length(ipni_id)==0L) return(NULL)
      return(ipni_id)
    }else if(length(has_distrib)>0L){
      db_ids_with_distrib <- db_ids[has_distrib]
      geo_ids <- switch(status, 'introduced' ={
        tmp_data <- kew_lookup_table[introduced==1][, fun(.SD), by = db_id][V1==TRUE]
        try(db_ids_with_distrib[keyvalExists(tmp_data, db_ids_with_distrib)], silent=TRUE)
      },
      'native' = {
        tmp_data <- kew_lookup_table[introduced==0][, fun(.SD), by = db_id][V1==TRUE]
        try(db_ids_with_distrib[keyvalExists(tmp_data, db_ids_with_distrib)], silent=TRUE)

      },
      'both' = {
        tmp_data <- kew_lookup_table[, fun(.SD), by = db_id][V1==TRUE]
        try(db_ids_with_distrib[keyvalExists(tmp_data, db_ids_with_distrib)], silent=TRUE)

      })
      if(inherits(geo_ids,'try-error') || length(geo_ids)==0L || all(is.na(geo_ids))){ # try last updates from the api
        results <- queryDistributionId(species_name,status=ifelse(status=='native','natives',status))
        if(is.null(results)) return(NULL)
        if(key(kew_checklist)!='ipni_id') data.table::setkey(kew_checklist,'ipni_id')
        db_id 		<- kew_checklist[results][['db_id']]
        geo_data 	<- kew_lookup_table[, fun(.SD), by = db_id][V1==TRUE]
        if(length(geo_data[['db_id']])==0L || all(is.na(geo_data[['db_id']]))) return(NULL)
        if(!haskey(geo_data)) data.table::setkey(geo_data, 'db_id')
        ipni_id 	<- results[keyvalExists(geo_data, db_id)]
        if(length(ipni_id)==0L) return(NULL)
        return(ipni_id)
      }
      return(unique(na.omit(geo_ids)))
    }else{
      results <- queryDistributionId(species_name,status=ifelse(status=='native','natives',status))
      if(is.null(results)) return(NULL)
      if(key(kew_checklist)!='ipni_id') data.table::setkey(kew_checklist,'ipni_id')
      db_id 		<- kew_checklist[results][['db_id']]
      geo_data 	<- kew_lookup_table[, fun(.SD), by = db_id][V1==TRUE]
      if(length(geo_data[['db_id']])==0L || all(is.na(geo_data[['db_id']]))) return(NULL)
      if(!haskey(geo_data)) data.table::setkey(geo_data, 'db_id')
      ipni_id 	<- results[keyvalExists(geo_data, db_id)]
      if(length(ipni_id)==0L) return(NULL)
      return(ipni_id)
    }
  }else{
    db_id 			<- ifelse(best_guess=='Unplaced Unplaced', kew_checklist[species_name][["db_id"]], kew_checklist[species_name][["accepted_db_id"]])
    if(is.na(db_id)){
      results <- queryDistributionId(species_name,status=ifelse(status=='native','natives',status))
      if(is.null(results)) return(NULL)
      if(key(kew_checklist)!='ipni_id') data.table::setkey(kew_checklist,'ipni_id')
      db_id 		<- kew_checklist[results][['db_id']]
      geo_data 	<- kew_lookup_table[, fun(.SD), by = db_id][V1==TRUE]
      if(length(geo_data[['db_id']])==0L || all(is.na(geo_data[['db_id']]))) return(NULL)
      if(!haskey(geo_data)) data.table::setkey(geo_data, 'db_id')
      ipni_id 	<- results[keyvalExists(geo_data, db_id)]
      return(ipni_id)
    }
    has_distrib	 	<- which(keyvalExists(kew_lookup_table, db_id, how=NULL))
    if(length(has_distrib)>0L){
      db_id_with_distrib <- db_id[has_distrib]
      geo_ids <- switch(status, 'introduced' ={
        tmp_data <- kew_lookup_table[introduced==1][, fun(.SD), by = db_id][V1==TRUE]
        try(db_id_with_distrib[keyvalExists(tmp_data, db_id_with_distrib)], silent=TRUE)
      },
      'native' = {
        tmp_data <- kew_lookup_table[introduced==0][, fun(.SD), by = db_id][V1==TRUE]
        try(db_id_with_distrib[keyvalExists(tmp_data, db_id_with_distrib)], silent=TRUE)
        },
      'both' = {
        tmp_data <- kew_lookup_table[, fun(.SD), by = db_id][V1==TRUE]
        try(db_ids_with_distrib[keyvalExists(tmp_data, db_ids_with_distrib)], silent=TRUE)
      })
      if(inherits(geo_ids,'try-error') || length(geo_ids)==0L || all(is.na(geo_ids))){ # try last updates from the api
        results <- queryDistributionId(species_name,status=ifelse(status=='native','natives',status))
        if(is.null(results)) return(NULL)
        if(key(kew_checklist)!='ipni_id') data.table::setkey(kew_checklist,'ipni_id')
        db_id 		<- kew_checklist[results][['db_id']]
        geo_data 	<- kew_lookup_table[, fun(.SD), by = db_id][V1==TRUE]
        if(length(geo_data[['db_id']])==0L || all(is.na(geo_data[['db_id']]))) return(NULL)
        if(!haskey(geo_data)) data.table::setkey(geo_data, 'db_id')
        ipni_id 	<- results[keyvalExists(geo_data, db_id)]
        return(ipni_id)
      }
      return(unique(na.omit(geo_ids)))
    }else{
      results <- queryDistributionId(species_name,status=ifelse(status=='native','natives',status))
      if(is.null(results)) return(NULL)
      if(key(kew_checklist)!='ipni_id') data.table::setkey(kew_checklist,'ipni_id')
      db_id 		<- kew_checklist[results][['db_id']]
      geo_data 	<- kew_lookup_table[, fun(.SD), by = db_id][V1==TRUE]
      if(length(geo_data[['db_id']])==0L || all(is.na(geo_data[['db_id']]))) return(NULL)
      if(!haskey(geo_data)) data.table::setkey(geo_data, 'db_id')
      ipni_id 	<- results[keyvalExists(geo_data, db_id)]
      return(ipni_id)
    }
  }
}

