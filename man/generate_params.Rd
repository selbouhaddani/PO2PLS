% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PO2PLS_functions.R
\name{generate_params}
\alias{generate_params}
\title{Generate parameter values of a PO2PLS model}
\usage{
generate_params(
  X,
  Y,
  r,
  rx,
  ry,
  alpha = 0.1,
  type = c("random", "o2m", "unit")
)
}
\arguments{
\item{X}{Numerical data matrix or positive integer. This parameter should either be a dataset \eqn{X} or the number of desired \eqn{X} variables.}

\item{Y}{Numerical data matrix or positive integer. This parameter should either be a dataset \eqn{Y} or the number of desired \eqn{Y} variables.}

\item{r}{Positive integer. Number of joint PLS components. Must be positive!}

\item{rx}{Non-negative integer. Number of orthogonal components in \eqn{X}. Can be 0}

\item{ry}{Non-negative integer. Number of orthogonal components in \eqn{Y}. Can be 0}

\item{alpha}{Numeric vector. The length should be either one or three, with each entry between 0 and 1. It represents the proportion of noise relative to the variation of \eqn{X}, \eqn{Y}, and \eqn{U}, respectively. If only one number is given, it is used for all three parts.}

\item{type}{Character. Should be one of "random", "o2m" or "unit". Specifies which kind of parameters should be generated. If "o2m" is chosen, \code{X} and \code{Y} should be data matrices.}
}
\value{
A list with
\describe{
\item{W}{\eqn{X} joint loadings}
\item{Wo}{\eqn{X} specific loadings}
\item{C}{\eqn{Y} joint loadings}
\item{Co}{\eqn{Y} specific loadings}
\item{B}{Regression matrix of \eqn{U} on \eqn{T}}
\item{SigT}{Covariance matrix of \eqn{T}}
\item{SigTo}{Covariance matrix of \eqn{To}}
\item{SigUo}{Covariance matrix of \eqn{Uo}}
\item{SigH}{Covariance matrix of \eqn{H}}
\item{sig2E}{Variance of \eqn{E}}
\item{sig2F}{Variance of \eqn{F}}
}
}
\description{
Generate parameter values of a PO2PLS model
}
\details{
A list of PO2PLS parameters are generated based on the value of \code{type}:
\describe{
\item{\code{type="random"}}{Variance parameters are randomly sampled from a uniform distribution on 1 and 3 (1 and 4 for \eqn{B}).}
\item{\code{type="o2m"}}{O2PLS is fitted to \code{X} and \code{Y} first using \code{\link{o2m}} from the OmicsPLS package, and the corresponding PO2PLS parameters are derived from the result.}
\item{\code{type="unit"}}{The diagonal of each covariance matrix is a decreasing sequence from the number of components to one.}
}
}
