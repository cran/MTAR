## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----load, echo=FALSE----------------------------------------------------
library(MTAR)

## ---- out.width = "650px", echo=FALSE, fig.cap="Fig 1: An overview of MTAR workflow. Light blue rectangle represents necessary input. Dark blue rectangle denotes the final output of MTAR function. Gray rectangle denotes the intermediate parameters."----
knitr::include_graphics("workflow.png")

## ---- out.width = "500px", echo=FALSE, fig.cap="Fig 2: Venn  diagram of sample overlap among traits in three studies."----
knitr::include_graphics("sampleoverlap.png")

## ----sumstats, eval=TRUE, message=FALSE----------------------------------
data("rawdata")
attach(rawdata)
head(traits.dat$Study1)
head(cov.dat$Study1)
head(geno.dat$Study1)
obs.stat <- Get_UV_from_data(traits = traits.dat,
                         covariates = cov.dat,
                         genotype = geno.dat,
                         covariance = TRUE)
obs.stat$U
detach(rawdata)

## ----cMTAR_Vse,eval = TRUE, message=FALSE--------------------------------
data("varU.example")
attach(varU.example)
obs.stat <- Get_UV_from_varU(U = U, varU = varU, R= R)
obs.stat$U
detach(varU.example)

## ----Beta, eval=TRUE-----------------------------------------------------
data("beta.example")
attach(beta.example)
obs.stat <- Get_UV_from_beta(Beta = Beta, Beta.se = Beta.se, R = R)
detach(beta.example)

## ----zeta, eval=FALSE----------------------------------------------------
#  data("zeta.example")
#  attach(zeta.example)
#  # Downloading independent common SNPs from 1000Genome data set.
#  githubURL <- "https://github.com/lan/MTAR/blob/master/indp_snps.1KG.rda?raw=true"
#  utils::download.file(githubURL,"1kgfile")
#  load("1kgfile")
#  zeta1 <- Get_zeta(Zscore = Zscore, Indp_common_snp = indp_snps.1KG)
#  detach(zeta.example)

## ----MTAR----------------------------------------------------------------
data("MTAR.example")
names(MTAR.example)

## ----cMTAR, warning=FALSE, message=FALSE, eval = TRUE--------------------
attach(MTAR.example)
pval <-  MTAR(U = U, V = V, MAF = MAF, genetic_cor.trait = genetic_cor.trait, zeta = zeta)
pval
detach(MTAR.example)

