% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HG Functions.R
\name{FIpmfHG.Tx.perfect}
\alias{FIpmfHG.Tx.perfect}
\title{Fisher Information for \eqn{T_x} (Perfect Case)}
\usage{
FIpmfHG.Tx.perfect(N, barN, Tx, b, method = c("AD", "ND", "PMF-HI", "PMF-SM"))
}
\arguments{
\item{N}{Integer. Total population size.}

\item{barN}{Numeric. Average group size.}

\item{Tx}{Integer. Number of target individuals in the population.}

\item{b}{Integer. Total number of groups.}

\item{method}{Character. Method to calculate the derivative. Options are:
\itemize{
\item "AD" (Analytic Derivative)
\item "ND" (Numerical Derivative)
\item "PMF" (PMF-based approximation)
}}
}
\value{
Numeric. Fisher Information for the given \eqn{T_x}.
}
\description{
Calculates Fisher Information for the perfect case using different methods.
}
