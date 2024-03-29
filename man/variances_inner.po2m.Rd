% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PO2PLS_functions.R
\name{variances_inner.po2m}
\alias{variances_inner.po2m}
\title{Calculate the variance covariance matrix of the estimated PO2PLS parameters}
\usage{
variances_inner.po2m(fit, X, Y)
}
\arguments{
\item{fit}{A PO2PLS fit of class po2m}

\item{X}{Numeric matrix. Other types will be coerced to matrix with \code{as.matrix} (if this is possible)}

\item{Y}{Numeric matrix. Other types will be coerced to matrix with \code{as.matrix} (if this is possible)}

\item{type_var}{String. Type of covariance matrix sought}
}
\value{
A covariance matrix and standard errors

A vector with the standard errors for B per component
}
\description{
Calculate the variance covariance matrix of the estimated PO2PLS parameters
}
\keyword{B}
\keyword{Calculate}
\keyword{coefficient}
\keyword{errors}
\keyword{for}
\keyword{inner}
\keyword{internal}
\keyword{relation}
\keyword{standard}
\keyword{the}
