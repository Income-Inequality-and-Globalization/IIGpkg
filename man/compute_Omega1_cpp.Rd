% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{compute_Omega1_cpp}
\alias{compute_Omega1_cpp}
\title{Sum in the inverse of the posterior variance for loadings/partial effects}
\usage{
compute_Omega1_cpp(invOmega0, availableObs, selectR, fPost, w_regs, Viarray)
}
\arguments{
\item{invOmega0}{inverse prior matrix}

\item{availableObs}{available observations}

\item{selectR}{selector matrix for missing y}

\item{fPost}{matrix of factors (backward sampled states via FFBS)}

\item{w_regs}{regressor matrix}

\item{Viarray}{VCOV array (npara x npara x TT) of cross-sectional unit i}
}
\value{
summation of kronecker products of appropriate dimension
}
\description{
Sum in the inverse of the posterior variance for loadings/partial effects
}
