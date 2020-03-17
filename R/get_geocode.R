#' Functions to get the geographical unit codes of the TDWG regions the species belong to.
#'
#' Get the TDWG code of the region for a given species and geographic level.
#'
#' @param species_name A two character string (Genus species) specifying the name of the species to search for.
#' @param species_id A character string specifying the id of the species in ipni or Kew world checklist.
#' @param status Optional. A character string specifying the type of geographic region where to search for the species. Can be : 'native', 'introduced' or 'both'
#' @param level A numeric integer between 1 and 4. The geographic level of the TWDG at which the code must be extracted from. Default is level 2.
#' @param use_name_matching A logical. Should the species names be matched with the Kew checklist (i.e. try to guess actual accepted name) ? Default is FALSE.
#' @param exclude Optional. A character string or a vector of names specifying the geographic region(s) to exclude from the search.
#'                          Can be : "NORTHERN AMERICA","SOUTHERN AMERICA", "ANTARCTICA", "EUROPE", "AFRICA", "AUSTRALASIA","ASIA-TEMPERATE", "ASIA-TROPICAL" or "PACIFIC"
#' @export
get_geocode <- function(species_name=NULL, species_id=NULL, status='native', level=2, use_name_matching =FALSE, exclude=NULL, time_delay=3){

  use_sp_id = FALSE
  geo_level 	<- switch(level, "continent_code_l1", "region_code_l2", "area_code_l3", "area_code_l3")

  if(use_name_matching && data.table::key(kew_checklist)!="full_name_without_family") data.table::setkey(kew_checklist,"full_name_without_family")
  if(!use_name_matching && data.table::key(kew_checklist)!="acc_full_name_without_family") data.table::setkey(kew_checklist,"acc_full_name_without_family")

  if(is.null(species_id)){
    if(is.null(species_name))
        stop("Please provide at least a species name or a species id")
    best_guess 	<- unique(kew_checklist[species_name][["acc_full_name_without_family"]])
  }else{
    best_guess 	<- unique(subset(kew_checklist, ipni_id==species_id | db_id==species_id)[["acc_full_name_without_family"]])
    use_sp_id =TRUE
  }

  if(is.na(best_guess)){	# try last updates from api ipni
    results <- queryDistribution(species_name,status=ifelse(status=='native','natives',status),use_acc_name=use_name_matching)
    if(is.null(results)) return(NULL)
    geo_codes <- switch(status, 'introduced' = try(kew_lookup_table[area_code_l3 %in% results$tdwgCode][introduced==1][!continent %in% exclude][[geo_level]],silent=TRUE),
                        'native' = try(kew_lookup_table[area_code_l3 %in% results$tdwgCode][introduced==0][!continent %in% exclude][[geo_level]],silent=TRUE),
                        'both' = try(kew_lookup_table[area_code_l3 %in% results$tdwgCode][!continent %in% exclude][[geo_level]],silent=TRUE))
    if(inherits(geo_codes,'try-error') || length(geo_codes)==0L || all(is.na(geo_codes))) return(NULL)
    return(unique(na.omit(geo_codes)))
  }

  if(length(best_guess)>1){
    db_ids <-  ifelse(best_guess=='Unplaced Unplaced', kew_checklist[species_name][["db_id"]], if(use_sp_id) subset(kew_checklist, ipni_id==species_id | db_id==species_id)[["accepted_db_id"]] else kew_checklist[species_name][["accepted_db_id"]])
    has_distrib	 <- which(keyvalExists(kew_lookup_table, db_ids, how=NULL))
    if(length(has_distrib)==0L){	# try last updates from api ipni
      Sys.sleep(time_delay) # sleep before trying to access api server
      results <- queryDistribution(species_name,status=ifelse(status=='native','natives',status))
      if(is.null(results)) return(NULL)
      geo_codes <- switch(status, 'introduced' = try(kew_lookup_table[area_code_l3 %in% results$tdwgCode][introduced==1][!continent %in% exclude][[geo_level]],silent=TRUE),
                          'native' = try(kew_lookup_table[area_code_l3 %in% results$tdwgCode][introduced==0][!continent %in% exclude][[geo_level]],silent=TRUE),
                          'both' = try(kew_lookup_table[area_code_l3 %in% results$tdwgCode][!continent %in% exclude][[geo_level]],silent=TRUE))
      if(inherits(geo_codes,'try-error') || length(geo_codes)==0L || all(is.na(geo_codes))) return(NULL)
      return(unique(na.omit(geo_codes)))
    }else if(length(has_distrib)>0L){
      db_ids_with_distrib <- db_ids[has_distrib]
      geo_codes <- switch(status, 'introduced' = try(kew_lookup_table[db_ids_with_distrib][introduced==1][!continent %in% exclude][[geo_level]],silent=TRUE),
                          'native' = try(kew_lookup_table[db_ids_with_distrib][introduced==0][!continent %in% exclude][[geo_level]],silent=TRUE),
                          'both' = try(kew_lookup_table[db_ids_with_distrib][!continent %in% exclude][[geo_level]],silent=TRUE))
      if(inherits(geo_codes,'try-error') || length(geo_codes)==0L || all(is.na(geo_codes))){ # try last updates from the api
        Sys.sleep(time_delay) # sleep before trying to access api server
        results <- queryDistribution(species_name,status=ifelse(status=='native','natives',status))
        if(is.null(results)) return(NULL)
        geo_codes <- switch(status, 'introduced' = try(kew_lookup_table[area_code_l3 %in% results$tdwgCode][introduced==1][!continent %in% exclude][[geo_level]],silent=TRUE),
                            'native' = try(kew_lookup_table[area_code_l3 %in% results$tdwgCode][introduced==0][!continent %in% exclude][[geo_level]],silent=TRUE),
                            'both' = try(kew_lookup_table[area_code_l3 %in% results$tdwgCode][!continent %in% exclude][[geo_level]],silent=TRUE))
        if(inherits(geo_codes,'try-error') || length(geo_codes)==0L || all(is.na(geo_codes))) return(NULL)
        return(unique(na.omit(geo_codes)))
      }
      return(unique(na.omit(geo_codes)))
    }else{
      Sys.sleep(time_delay) # sleep before trying to access api server
      results <- queryDistribution(species_name,status=ifelse(status=='native','natives',status))
      if(is.null(results)) return(NULL)
      geo_codes <- switch(status, 'introduced' = try(kew_lookup_table[area_code_l3 %in% results$tdwgCode][introduced==1][!continent %in% exclude][[geo_level]],silent=TRUE),
                          'native' = try(kew_lookup_table[area_code_l3 %in% results$tdwgCode][introduced==0][!continent %in% exclude][[geo_level]],silent=TRUE),
                          'both' = try(kew_lookup_table[area_code_l3 %in% results$tdwgCode][!continent %in% exclude][[geo_level]],silent=TRUE))
      if(inherits(geo_codes,'try-error') || length(geo_codes)==0L || all(is.na(geo_codes))) return(NULL)
      return(unique(na.omit(geo_codes)))
    }
  }else{
    db_id 			<- if(best_guess=='Unplaced Unplaced') kew_checklist[species_name][["db_id"]] else if(use_sp_id) subset(kew_checklist, ipni_id==species_id | db_id==species_id)[["accepted_db_id"]] else kew_checklist[species_name][["accepted_db_id"]]
    if(all(is.na(db_id))){
      Sys.sleep(time_delay) # sleep before trying to access api server
      results <- queryDistribution(species_name,status=ifelse(status=='native','natives',status),use_acc_name=use_name_matching)
      if(is.null(results)) return(NULL)
      geo_codes <- switch(status, 'introduced' = try(kew_lookup_table[area_code_l3 %in% results$tdwgCode][introduced==1][!continent %in% exclude][[geo_level]],silent=TRUE),
                          'native' = try(kew_lookup_table[area_code_l3 %in% results$tdwgCode][introduced==0][!continent %in% exclude][[geo_level]],silent=TRUE),
                          'both' = try(kew_lookup_table[area_code_l3 %in% results$tdwgCode][!continent %in% exclude][[geo_level]],silent=TRUE))
      if(inherits(geo_codes,'try-error') || length(geo_codes)==0L || all(is.na(geo_codes))) return(NULL)
      return(unique(na.omit(geo_codes)))
    }
    has_distrib	 	<- which(keyvalExists(kew_lookup_table, db_id, how=NULL))
    if(length(has_distrib)>0L){
      geo_codes <- switch(status, 'introduced' = try(kew_lookup_table[db_id][introduced==1][!continent_code_l1 %in% exclude][[geo_level]],silent=TRUE),
                          'native' = try(kew_lookup_table[db_id][introduced==0][!continent %in% exclude][[geo_level]],silent=TRUE),
                          'both' = try(kew_lookup_table[db_id][!continent %in% exclude][[geo_level]],silent=TRUE))
      if(inherits(geo_codes,'try-error') || length(geo_codes)==0L || all(is.na(geo_codes))){ # try last updates from the api
        Sys.sleep(time_delay) # sleep before trying to access api server
        results <- queryDistribution(species_name,status=ifelse(status=='native','natives',status),use_acc_name=use_name_matching)
        if(is.null(results)) return(NULL)
        geo_codes <- switch(status, 'introduced' = try(kew_lookup_table[area_code_l3 %in% results$tdwgCode][introduced==1][!continent %in% exclude][[geo_level]],silent=TRUE),
                            'native' = try(kew_lookup_table[area_code_l3 %in% results$tdwgCode][introduced==0][!continent %in% exclude][[geo_level]],silent=TRUE),
                            'both' = try(kew_lookup_table[area_code_l3 %in% results$tdwgCode][!continent %in% exclude][[geo_level]],silent=TRUE))
        if(inherits(geo_codes,'try-error') || length(geo_codes)==0L || all(is.na(geo_codes))) return(NULL)
        return(unique(na.omit(geo_codes)))
      }
      return(unique(na.omit(geo_codes)))
    }else{
      Sys.sleep(time_delay) # sleep before trying to access api server
      results <- queryDistribution(species_name,status=ifelse(status=='native','natives',status),use_acc_name=use_name_matching)
      if(is.null(results)) return(NULL)
      geo_codes <- switch(status, 'introduced' = try(kew_lookup_table[area_code_l3 %in% results$tdwgCode][introduced==1][!continent %in% exclude][[geo_level]],silent=TRUE),
                          'native' = try(kew_lookup_table[area_code_l3 %in% results$tdwgCode][introduced==0][!continent %in% exclude][[geo_level]],silent=TRUE),
                          'both' = try(kew_lookup_table[area_code_l3 %in% results$tdwgCode][!continent %in% exclude][[geo_level]],silent=TRUE))
      if(inherits(geo_codes,'try-error') || length(geo_codes)==0L || all(is.na(geo_codes))) return(NULL)
      return(unique(na.omit(geo_codes)))
    }
  }
}

