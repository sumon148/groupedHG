% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BN Functions.R
\name{ExpVarBN.imperfect}
\alias{ExpVarBN.imperfect}
\title{Expected Value and Variance for Grouped Binomial (Imperfect Case at either group or item level)}
\usage{
ExpVarBN.imperfect(N, Tx, barN, b, delta, lambda, type = c("group", "item"))
}
\arguments{
\item{N}{Numeric. The total population size.}

\item{Tx}{Numeric. The number of target individuals in the population.}

\item{barN}{Numeric. The average number of individuals per group.}

\item{b}{Numeric. The total number of groups.}

\item{delta}{Numeric. The group or item-level false negative rate (probability that a positive case is classified as negative).}

\item{lambda}{Numeric. The group or item-level false positive rate (probability that a negative case is classified as positive).}

\item{type}{Character. The type of calculation: either "group" for group-level or "item" for item-level calculation. Default is "group".}
}
\value{
A list containing:
\itemize{
\item \code{expectation}: The expected number of detected target individuals based on the chosen type.
\item \code{variance}: The variance of the number of detected target individuals based on the chosen type.
}
The returned list contains the expectation and variance for the selected model type (group-level or item-level).
}
\description{
This function calculates the expected value and variance for a binomial distribution in the imperfect detection case,
with the option to compute the results at either the group level or item level. The function routes the calculation
to the appropriate sub-function depending on the specified type ("group" or "item").
}
\details{
This unified function selects the appropriate method for calculating expectation and variance based on the \code{type} parameter.
If \code{type} is "group", the function calls \code{ExpVarBN.imperfect.group}, which calculates the expectation and variance at the group level.
If \code{type} is "item", the function calls \code{ExpVarBN.imperfect.item}, which calculates the expectation and variance at the item level.
The function then returns the results for the selected type of model.
}
\examples{
# Example: Expectation and variance calculation for group-level model
ExpVarBN.imperfect(N = 1000, Tx = 10, barN = 50, b = 5, delta = 0.1, lambda = 0.05, type = "group")

# Example: Expectation and variance calculation for item-level model
ExpVarBN.imperfect(N = 1000, Tx = 10, barN = 50, b = 5, delta = 0.1, lambda = 0.05, type = "item")

}
