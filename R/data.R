#' PNPLA2 Gene
#'
#' PNPLA2 Gene Inofmation required by MTAR function.
#' @docType data
#' @usage data(MTAR.example)
#' @format A list with 6 sublists:
#' \describe{
#'   \item{annotation}{some annotation information of gene PNPLA2}
#'   \item{U}{a list containing the summary statistics for each trait}
#'   \item{V}{a list containing the covariance matrix of summary statistics for each trait}
#'   \item{MAF}{the minor allele frequency of all the rare variants in this gene}
#'   \item{genetic_cor.trait}{the genetic correlation among traits}
#'   \item{zeta}{the overlapping estimation matrix, which is approximated by the sample correlation matrix of the scaled summary statistics over a large number of independent null common SNPs}}
"MTAR.example"


#' Example1: individual-level data
#'
#' An example of calculating U and V for MTAR given the individual-level data set.
#' @docType data
#' @usage data(rawdata)
#' @format A list with 3 sublists:
#' \describe{
#'   \item{geno.dat}{a list of genotype data, each sublist is \eqn{n \times m} matrix for each study}
#'   \item{traits.dat}{a list of trait data, each sublist is \eqn{n \times K} matrix for each study}
#'   \item{cov.dat}{a list of covariates data, each sublist is \eqn{n \times D} matrix for each study}}
#' There are 3 studies, 3 continuous traits and 10 rare variants. Specifically, there are 1500 subjects in study1, but each subject only has one trait measurement. In study2 and study3, the sample size is 500 and each subject has two or three traits measurements.
"rawdata"

#' Example2: the summary statistics and their variance
#'
#' An example of calculating U and V for MTAR given the summary statistics and its variance.
#' @docType data
#' @usage data(varU.example)
#' @format A list with 3 sublists:
#' \describe{
#'   \item{U}{a list of summary statistics, each sublist is \eqn{m \times K} matrix for each study}
#'   \item{varU}{a list of variance of summary statistics, each sublist is \eqn{m \times K} matrix for each study}
#'   \item{R}{a SNP correlation matrix for the union of SNPs appearing in all the studies}}
"varU.example"

#' Example3: the genetic effect estimates and their standard errors
#'
#' An example of calculating U and V for MTAR given the genetic effect estimates and their standard errors.
#' @docType data
#' @usage data(beta.example)
#' @format A list with 3 sublists:
#' \describe{
#'   \item{Beta}{a list of genetic effect estimates, each sublist is \eqn{m \times K} matrix for each study}
#'   \item{Beta.se}{a list of standard errors of genetic effect estimates, each sublist is \eqn{m \times K} matrix for each study}
#'   \item{R}{a SNP correlation matrix for the union of SNPs appearing in all the studies}}
"beta.example"

#' Example4: the summary statistics of 737 common and null SNPs
#'
#' An example of estimating matrix \eqn{\zeta} given the summary statistics information of null and common SNPs.
#' @docType data
#' @usage data(zeta.example)
#' @format A list with 1 sublist:
#' \describe{
#'   \item{Zscore}{a list, each sublist contains the Z-scores of 737 null and common SNPs for each trait}
#'   }
"zeta.example"
