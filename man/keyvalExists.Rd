% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/keyvalExists.R
\name{keyvalExists}
\alias{keyvalExists}
\title{Check for a key value whether or not a key value exist}
\usage{
keyvalExists(DT, keyval, how = NULL)
}
\arguments{
\item{DT}{A data.table object. DT must have a key set.}

\item{keyval}{A character/numeric vector of key values.}

\item{how}{A character string specifying the type of logical output to return comparison to perform. Either \code{'all'} or \code{'any'}}
}
\value{
A logical value if \code{how} is set to 'all' or 'any'. If \code{how} is set to NULL (default), a logical vector is returned. TRUE if the key value exists in DT, otherwise FALSE.
}
\description{
Check whether or not a key value exist in a data.table
}
