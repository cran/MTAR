% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{rawdata}
\alias{rawdata}
\title{Example1: individual-level data}
\format{A list with 3 sublists:
\describe{
  \item{geno.dat}{a list of genotype data, each sublist is \eqn{n \times m} matrix for each study}
  \item{traits.dat}{a list of trait data, each sublist is \eqn{n \times K} matrix for each study}
  \item{cov.dat}{a list of covariates data, each sublist is \eqn{n \times D} matrix for each study}}
There are 3 studies, 3 continuous traits and 10 rare variants. Specifically, there are 1500 subjects in study1, but each subject only has one trait measurement. In study2 and study3, the sample size is 500 and each subject has two or three traits measurements.}
\usage{
data(rawdata)
}
\description{
An example of calculating U and V for MTAR given the individual-level data set.
}
\keyword{datasets}
