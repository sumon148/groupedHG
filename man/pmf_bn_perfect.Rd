% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BN Functions.R
\name{pmf_bn_perfect}
\alias{pmf_bn_perfect}
\title{Probability Mass Function for Grouped Binomial Distribution under Perfect Test}
\usage{
pmf_bn_perfect(ty, Tx, b, barN, N)
}
\arguments{
\item{ty}{Integer. Number of groups where the target outcome is observed.}

\item{Tx}{Integer. Number of target individuals in the population.}

\item{b}{Integer. Total number of groups.}

\item{barN}{Integer. Average number of individuals per group.}

\item{N}{Integer. Total population size.}
}
\value{
Numeric. The probability mass function (PMF) value for the given parameters.
}
\description{
Computes the probability mass function (PMF) for a grouped Binomial distribution
under the assumption of a perfect test (no false positives or false negatives).
}
\details{
The function calculates the probability of observing \code{ty} positive groups in \code{b} total groups,
using binomial probabilities for detecting target individuals in a group.
}
\examples{
pmf_bn_perfect(ty = 3, Tx = 10, b = 5, barN = 2, N = 20)

}
