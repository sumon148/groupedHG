% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HG Functions.R
\name{ExpVarHG.imperfect.group}
\alias{ExpVarHG.imperfect.group}
\title{Calculate Expectation and Variance with Group-Level Sensitivity and Specificity}
\usage{
ExpVarHG.imperfect.group(N, Tx, barN, b, delta, lambda)
}
\arguments{
\item{N}{Total number of items.}

\item{Tx}{Number of contaminated items.}

\item{barN}{Number of items per group.}

\item{b}{Number of groups.}

\item{delta}{Sensitivity of the test.}

\item{lambda}{Specificity of the test.}
}
\value{
A list with the following elements:
\describe{
\item{V}{Probability that a single group tests positive.}
\item{W}{Covariance between two groups.}
\item{expectation}{Expected number of positive groups.}
\item{variance}{Variance of the number of positive groups.}
}
}
\description{
Computes the expectation and variance of positive test groups using
hypergeometric probabilities at the group level.
}
