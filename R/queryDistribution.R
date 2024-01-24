#' Functions to get the geographic distribution information of a plant species using the Plant Of the World Online (POWO) API.
#'
#' Get information about the geographic distribution of a plant species.
#'
#' @param binomial A two character string (Genus species) specifying the name of the species to search for.
#' @param ipni_id A character string specifying the IPNI id of the species
#' @param status Optional. A character string specifying the geographic range where to search for the species. Can be : 'natives' or 'introduced'
#' @param use_acc_name A logical. Should the species name be considered 'accepted' ? Default is TRUE.
#' @export
queryDistribution <- function(binomial, ipni_id = NULL, status='natives', verbose=TRUE, use_acc_name=TRUE, ...){

  # Set the API link

  POWO_API = 'http://www.plantsoftheworldonline.org/api/2'

  if(is.null(ipni_id)){

    # Look up the IPNI id

    powo.url <- paste0(POWO_API,'/','search?perPage=500&cursor=%2A&q=', RCurl::curlEscape(binomial)) #

    # Try and catch the POWO results

    powo.json <- try(RCurl::getURLAsynchronous(powo.url),silent = TRUE)

    if(inherits(powo.json,'try-error')) return(NULL)

    # Parse them

    results <- try(jsonlite::parse_json(powo.json, simplifyVector = TRUE)$results, silent = TRUE)

    acc_idx = if(use_acc_name) results$accepted else rep(TRUE, times=length(results$accepted))

    if(inherits(results,'try-error') || !('fqId' %in% names(results)) || !any(acc_idx)){return(NULL)} else {

      id <- results$fqId[acc_idx] # Get all ids

    }

  } else {

    id <- paste0('urn:lsid:ipni.org:names:',ipni_id)

  }

  # Query API

  powo.url <- paste0(POWO_API,'/','taxon/',id,'?fields=distribution') #

  # Try and catch the POWO results

  powo.json <- try(RCurl::getURLAsynchronous(powo.url),silent = TRUE)

  if(inherits(powo.json,'try-error')) return(NULL)

  # Parse url to check which species id has distribution information that fits
  if(status=='both')
    hasdistrib <- sapply(powo.json, function(x) length(jsonlite::parse_json(x, simplifyVector = TRUE)[['distribution']])>0L, USE.NAMES=FALSE)
  else
    hasdistrib <- sapply(powo.json, function(x) length(jsonlite::parse_json(x, simplifyVector = TRUE)[['distribution']][[status]])>0L, USE.NAMES=FALSE)

  if(!any(hasdistrib))
    return(NULL)

  # Get the first
  powo.json <- powo.json[hasdistrib][1]

  # Parse to get the distribution information
  results <- try({
    if(status=='both'){
      res <- jsonlite::parse_json(powo.json,simplifyVector = TRUE)$distribution
      if(!is.null(res) && is.list(res))
        do.call(rbind, c(jsonlite::parse_json(powo.json,simplifyVector = TRUE)$distribution, list(make.row.names=FALSE)))
      else
        res
    }
    else
     jsonlite::parse_json(powo.json,simplifyVector = TRUE)$distribution[[status]]
    }, silent = TRUE)

  # Return results

  if(is.null(results)) { return(NULL) } else {return(results)}

}


#' Functions to get the IPNI id of a plant species using the Plant Of the World Online (POWO) API.
#'
#' Get the IPNI id of a species with geographic distribution information.
#'
#' @param binomial A two character string (Genus species) specifying the name of the species to search for.
#' @param status Optional. A character string specifying the geographic range where to search for the species. Can be : 'natives' or 'introduced'
#' @param use_acc_name A logical. Should the species name be considered 'accepted' ? Default is TRUE.
#' @export
queryDistributionId <- function(binomial, status='natives', use_acc_name=TRUE,...){

  # Set the API link

  POWO_API = 'http://www.plantsoftheworldonline.org/api/2'

  # Look up the IPNI id

  powo.url <- paste0(POWO_API,'/','search?perPage=500&cursor=%2A&q=', RCurl::curlEscape(binomial) ) #

  # Try and catch the POWO results

  powo.json <- try(RCurl::getURLAsynchronous(powo.url),silent = TRUE)

  if(inherits(powo.json,'try-error')) return(NULL)

  # Parse them

  results <- try(jsonlite::parse_json(powo.json,simplifyVector = TRUE)$results, silent = TRUE)
  
  if(inherits(results,'try-error')) return(NULL)

  acc_idx = if(use_acc_name) results$accepted else rep(TRUE, times=length(results$accepted))

  if(inherits(results,'try-error') || !('fqId' %in% names(results)) || !any(acc_idx)){return(NULL)} else {

    id <- results$fqId[acc_idx] # Get all ids

  }

  # Query API

  powo.url <- paste0(POWO_API,'/','taxon/',id,'?fields=distribution') #

  # Try and catch the POWO results

  powo.json <- try(RCurl::getURLAsynchronous(powo.url),silent = TRUE)

  if(inherits(powo.json,'try-error')) return(NULL)

  # Parse them to check which one has distribution information that fits
  if(status=='both')
    hasdistrib <- sapply(powo.json, function(x) length(jsonlite::parse_json(x, simplifyVector = TRUE)[['distribution']])>0L, USE.NAMES=FALSE)
  else
    hasdistrib <- sapply(powo.json, function(x) length(jsonlite::parse_json(x, simplifyVector = TRUE)[['distribution']][[status]])>0L, USE.NAMES=FALSE)

  if(!any(hasdistrib))
    return(NULL)

  ipni_ids <- stringr::str_extract(id[hasdistrib],"(?<=\\:)[[:digit:]].+")

  return(setNames(ipni_ids,rep("ipni_id",length(results))))

}
