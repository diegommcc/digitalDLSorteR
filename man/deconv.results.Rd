% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R
\docType{methods}
\name{deconv.results}
\alias{deconv.results}
\alias{deconv.results,DigitalDLSorter-method}
\alias{deconv.results<-}
\alias{deconv.results<-,DigitalDLSorter-method}
\title{Get and set \code{deconv.results} slot in a
\code{\linkS4class{DigitalDLSorter}} object}
\usage{
deconv.results(object, name.data = NULL)

deconv.results(object, name.data = NULL) <- value
}
\arguments{
\item{object}{\code{\linkS4class{DigitalDLSorter}} object.}

\item{name.data}{Name of the data. If \code{NULL} (by default), all
results contained in the \code{deconv.results} slot are returned.}

\item{value}{List whose names are the reference of the stored results.}
}
\description{
Get and set \code{deconv.results} slot in a
\code{\linkS4class{DigitalDLSorter}} object
}
