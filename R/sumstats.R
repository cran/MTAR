#' Compute the summary statistics given the individual-level data
#'
#' This function allows you to calculate the score summary statistics \eqn{U} and their covariance matrix \eqn{V} for MTAR, given the traits, covariates and genotype data sets. If one trait only takes no more than two values, this function will treat this trait as a binary trait.
#' @param traits a numeric list, each sublist containing trait information for each study. In each study, a numeric \eqn{n \times K} matrix with each row as an independent individual and each column as a separate trait. If subject \eqn{i} is without trait \eqn{k}, the corresponding value is set as NA. The number of traits in each study can be different but the names of traits are required.
#' @param covariates a numeric list, each sublist containing covariates information for each study. In each study, a numeric \eqn{n \times D} matrix with each row as an independent individual and each column as a covariate.
#' @param genotype a numeric list, each sublist containing genotype information for each study. In each study, a numeric \eqn{n \times m} matrix with each row as an independent individual and each column as a SNP Each genotype should be coded as 0, 1, 2, and NA for AA, Aa, aa, and missing, where A and a represents a major and minor allele, respectively. The number of SNPs in each study can be different but the names of SNPs are required. Also, the number of studies must be the same in genotype, covariates and traits lists. The order of subject ID must be the same among traits, covariates, and genotype within each study.
#' @param covariance a logical value indicating whether to calculate the covariance matrix of score summary statistics \eqn{U}. The default value is TRUE. If covariance is set as FALSE, then only the diagonal values of the covairance matrix are calculated, which is faster. In estimating the zeta matrix correcting for overlap samples, we recommend set covariance as FALSE in calculating the summary statistics for common variants. Since the number of common variants may be large and time-consuming.
#' @return A list containing summary statistics for each trait. If covariance is TRUE, the score summary statistics \eqn{U} and its covariance matrix \eqn{V} are returned. Otherwise, only \eqn{U} and the diagonal elements of covariance matrix  are returned.
#' @author Lan Luo
#' @references Hu, Y.J., Berndt, S.I., Gustafsson, S., Ganna, A., Mägi, R., Wheeler, E., Feitosa, M.F., Justice, A.E., Monda, K.L., Croteau-Chonka, D.C. and Day, F.R., 2013. Meta-analysis of gene-level associations for rare variants based on single-variant statistics. The American Journal of Human Genetics, 93(2), pp.236-248.
#' @references Tang, Z.Z. and Lin, D.Y., 2015. Meta-analysis for discovering rare-variant associations: statistical methods and software programs. The American Journal of Human Genetics, 97(1), pp.35-53.
#' @export
#' @examples
#' data(rawdata)
#' attach(rawdata)
#' obs.stat <- Get_UV_from_data(traits = traits.dat,
#'                              covariates = cov.dat,
#'                              genotype = geno.dat,
#'                              covariance = TRUE)
#' obs.stat
#' detach(rawdata)
Get_UV_from_data <- function(traits, covariates, genotype, covariance = TRUE){
  # study.list <- names(traits)

  if(length(covariates) != length(traits) | length(genotype) != length(traits)) {
    stop("The number of studies is not equal among traits, covariates or genotype")
  }

  obs.stat.bystudy <- list()
  trait.list <- NULL
  snp.list <- NULL
  for(study_id in 1:length(traits)) {
    trait.bystudy <- traits[[study_id]]
    covariate.bystudy <- covariates[[study_id]]
    genotype.bystudy <- genotype[[study_id]]
    if(is.null(colnames(trait.bystudy))){
      stop(paste0("No input of names of traits in study ", study_id, ", calculation of summary statistics cannot continue."))
    }
    trait.list <- c(trait.list, colnames(trait.bystudy))
    if(is.null(colnames(genotype.bystudy))){
      stop(paste0("No input of names of SNPs in study ", study_id, ", calculation of summary statistics cannot continue."))
    }
    snp.list <- c(snp.list, colnames(genotype.bystudy))
    obs.stat.bystudy[[study_id]] <- cal.sumstats(traits = trait.bystudy, covariates = covariate.bystudy, genotype = genotype.bystudy, covariance = covariance)
  }
  # names(obs.stat.bystudy) <- study.list
  trait.list <- unique(trait.list)
  snp.list <- unique(snp.list)

  if(length(traits) == 1) {
    message(paste0("There is 1 study in the data set and ", length(trait.list), " traits in total."))
  } else{
    message(paste0("There are ", length(traits), " studies in the data set and ", length(trait.list), " traits in total."))
  }

  obs.stat <- list()
  m <- length(snp.list)
  for(trait_id in 1:length(trait.list)) {
    # if one trait shows in multiple studies, add them up
    U.bytrait <- numeric(m)
    names(U.bytrait) <- snp.list
    if(covariance) {
      V.bytrait <- matrix(0, m, m)
      colnames(V.bytrait) <- rownames(V.bytrait) <- snp.list
    }else{
      V.bytrait <- numeric(m)
      names(V.bytrait) <- snp.list
    }

    for(study_id in 1:length(traits)) {
      sublist_id <- which(names(obs.stat.bystudy[[study_id]]) == trait.list[trait_id])
      if(length(sublist_id) != 0) {
        sublist <- obs.stat.bystudy[[study_id]][[sublist_id]]
        sublist_snps <- names(sublist$U)
        sublist_order <- order(match(sublist_snps, snp.list))
        list_order <- which(snp.list%in% sublist_snps)
        U.bytrait[list_order] <- U.bytrait[list_order] + sublist$U[sublist_order]
        if(covariance) {
          V.bytrait[list_order, list_order] <- V.bytrait[list_order, list_order] + sublist$V[sublist_order, sublist_order]
        }else{
          V.bytrait[list_order] <- V.bytrait[list_order] + sublist$V[sublist_order]
        }
      }
    }
    obs.stat[[trait_id]] <- list(U = U.bytrait, V = V.bytrait)
  }
  names(obs.stat) <- trait.list
  ret <- convert_toUV(obs.stat)
  return(ret)
}

cal.sumstats <- function(traits, covariates, genotype, covariance = TRUE){
  K <- ncol(traits)
  trait.list <- colnames(traits)

  obs.stat <- list()
  n <- nrow(traits)
  snp.list <- colnames(genotype)
  m <- length(snp.list)
  if(covariance) {
    for(k in 1:K){
      Y <- traits[, k]
      valid.id <- which(!is.na(Y)) # extract valid value
      trait.dat <- Y[valid.id]

      if(length(unique(trait.dat)) <= 2) {
        if(length(unique(trait.dat)) < 2) {
          warning(paste0("There is only one value in trait ", trait.list[k]))
        }
        type <- "binary"
        trait.dat <- as.factor(trait.dat)
        levels(trait.dat) <- c(0, 1)[1:length(unique(trait.dat))]
        trait.dat <- as.numeric(as.character(trait.dat))
      }else{
        type <- "continuous"
      }

      geno.dat <- genotype[valid.id, ]
      cov.dat <- model.matrix(trait.dat ~ covariates[valid.id, ])

      if(type == "continuous"){
        message(paste0("Calculating summary statistics for continuous trait", k))
        # continuous traits
        W <- t(geno.dat) %*% geno.dat-
          (t(geno.dat) %*% cov.dat) %*%
          solve(t(cov.dat) %*% cov.dat) %*%
          (t(cov.dat) %*% geno.dat)

        mod <- summary(lm(trait.dat ~ cov.dat))
        s2 <- mod$sigma^2
        resid <- mod$residuals
        U <- as.vector(t(geno.dat) %*% resid) / s2
        V <- W / s2
        names(U) <- snp.list
        if(length(snp.list) > 2) {
          colnames(V) <- rownames(V) <- snp.list
        }else{
          names(V) <- snp.list
        }
        obs.stat[[k]] <- list(U = U, V = V)
      }else{
        message(paste0("Calculating summary statistics for binary trait ", k))
        mod <- glm(trait.dat ~ cov.dat, family = "binomial")

        sum.U <- numeric(m)
        sum1 <- matrix(0, m, m)
        sum2 <- matrix(0, m, ncol(cov.dat))
        sum3 <- matrix(0, ncol(cov.dat), ncol(cov.dat))
        sum4 <- matrix(0, ncol(cov.dat), m)
        for(i in 1:length(trait.dat)){
          lambdaX <- as.numeric(mod$coefficients[-2] %*% cov.dat[i, ])
          b1 <- exp(lambdaX) / (1 + exp(lambdaX))
          b2 <- exp(lambdaX) / ((1 + exp(lambdaX)))^2

          U.part1 <- (trait.dat[i] - b1) * geno.dat[i, ]
          U.part1[is.na(U.part1)] <- 0
          V.part1 <- b2 * geno.dat[i, ] %*% t(geno.dat[i, ])
          V.part1[is.na(V.part1)] <- 0
          V.part2 <- b2 * geno.dat[i, ] %*% t(cov.dat[i, ])
          V.part2[is.na(V.part2)] <- 0
          V.part3 <- b2 * cov.dat[i, ] %*% t(cov.dat[i, ])
          V.part3[is.na(V.part3)] <- 0
          V.part4 <- b2 * cov.dat[i, ] %*% t(geno.dat[i, ])
          V.part4[is.na(V.part4)] <- 0

          sum.U <- sum.U + U.part1
          sum1 <- sum1 + V.part1
          sum2 <- sum2 + V.part2
          sum3 <- sum3 + V.part3
          sum4 <- sum4 + V.part4
        }
        U <- sum.U
        V <- sum1 - sum2 %*% solve(sum3) %*% sum4
        names(U) <- snp.list
        if(length(snp.list) > 2) {
          colnames(V) <- rownames(V) <- snp.list
        }else{
          names(V) <- snp.list
        }
        obs.stat[[k]] <- list(U = U, V = V)
      }
    }
    names(obs.stat) <- trait.list
  } else{
    for(k in 1:K) {
      Y <- traits[, k]
      valid.id <- which(!is.na(Y)) # extract valid value
      trait.dat <- Y[valid.id]

      if(length(unique(trait.dat)) <= 2) {
        if(length(unique(trait.dat)) < 2) {
          warning(paste0("There is only one value in trait ", trait.list[k]))
        }
        type <- "binary"
        trait.dat <- as.factor(trait.dat)
        levels(trait.dat) <- c(0, 1)[1:length(unique(trait.dat))]
        trait.dat <- as.numeric(as.character(trait.dat))
      }else{
        type <- "continuous"
      }

      geno.dat <- genotype[valid.id, ]
      cov.dat <- model.matrix(trait.dat ~ covariates[valid.id, ])

      if(type == "continuous"){
        message(paste0("Calculating summary statistics for continuous trait", k))
        # continuous traits
        mod <- summary(lm(trait.dat ~ cov.dat))
        s2 <- mod$sigma^2
        resid <- mod$residuals
        U <- as.vector((t(geno.dat) %*% resid)/s2)
        V <- numeric(m)
        for(j in 1:m){
          G <- as.matrix(geno.dat[, j])
          V[j] <- as.numeric(t(G) %*% G - (t(G) %*% cov.dat) %*%
                               solve(t(cov.dat) %*% cov.dat) %*%
                               (t(cov.dat) %*% G))/s2
        }
        names(U) <- names(V) <- snp.list
        obs.stat[[k]] <- list(U = U, V = V)
      } else {
        message(paste0("Calculating summary statistics for binary trait", k))
        # binary traits
        mod <- glm(trait.dat ~ cov.dat, family = "binomial")

        U <- numeric(m)
        V <- numeric(m)
        for(j in 1:m){
          G <- as.matrix(geno.dat[, j])
          sum1 <- matrix(0, 1, 1)
          sum2 <- matrix(0, 1, ncol(cov.dat))
          sum3 <- matrix(0, ncol(cov.dat), ncol(cov.dat))
          sum4 <- matrix(0, ncol(cov.dat), 1)
          for(i in 1:n){
            lambdaX <- as.numeric(mod$coefficients[-2] %*% cov.dat[i, ])
            b1 <- exp(lambdaX) / (1 + exp(lambdaX))
            b2 <- exp(lambdaX) / ((1 + exp(lambdaX)))^2

            U.part1 <- (trait.dat[i] - b1) * G[i]
            U.part1[is.na(U.part1)] <- 0
            V.part1 <- b2 * G[i] %*% t(G[i])
            V.part1[is.na(V.part1)] <- 0
            V.part2 <- b2 * G[i] %*% t(cov.dat[i, ])
            V.part2[is.na(V.part2)] <- 0
            V.part3 <- b2 * cov.dat[i, ] %*% t(cov.dat[i, ])
            V.part3[is.na(V.part3)] <- 0
            V.part4 <- b2 * cov.dat[i, ] %*% t(G[i])
            V.part4[is.na(V.part4)] <- 0

            U[j] <- U[j] + U.part1
            sum1 <- sum1 + V.part1
            sum2 <- sum2 + V.part2
            sum3 <- sum3 + V.part3
            sum4 <- sum4 + V.part4
          }
          V[j] <- sum1 - sum2 %*% solve(sum3) %*% sum4
        }
        names(U) <- names(V) <- snp.list
        obs.stat[[k]] <- list(U = U, V = V)
      }
    }
    names(obs.stat) <- trait.list
  }
  return(obs.stat)
}

#' Compute the summary statistics given the score statistics and their variance.
#'
#' This function allows you to calculate the score summary statistics \eqn{U} and their covariance matrix \eqn{V} for MTAR, given the score summary statistics and their variance.
#' @param U a numeric list, each sublist containing score summary statistics \eqn{U} for each study. In each study, a numeric \eqn{m \times K} matrix with each row as a SNP and each column as a separate trait. The number of traits and the number of SNPs in each study can be different but their names are required.
#' @param varU a numeric list, each sublist containing the variance of score summary statistics information for each study. In each study, a numeric \eqn{m \times K} matrix with each row as a SNP and each column as a separate trait.
#' @param R a SNP correlation matrix, which should contain the correlation of all the SNPs in these studies.
#' @return A list containing summary statistics for each traits, the score summary statistics \eqn{U} and their covariance matrix \eqn{V}.
#' @author Lan Luo
#' @export
#' @examples
#' data("varU.example")
#' attach(varU.example)
#' obs.stat <- Get_UV_from_varU(U = U, varU = varU, R= R)
#' obs.stat
#' detach(varU.example)

Get_UV_from_varU <- function(U, varU, R) {

  if(length(U) != length(varU)) {
    stop("The number of studies is not equal in U and varU")
  }
  trait.list <- NULL
  snp.list <- NULL
  for(study_id in 1:length(U)) {
    if(is.null(colnames(U[[study_id]]))){
      stop(paste0("No input of names of traits in study ", study_id, ", calculation of summary statistics cannot continue."))
    }
    trait.list <- c(trait.list, colnames(U[[study_id]]))
    if(is.null(rownames(U[[study_id]]))){
      stop(paste0("No input of names of SNPs in study ", study_id, ", calculation of summary statistics cannot continue."))
    }
    snp.list <- c(snp.list, rownames(U[[study_id]]))
  }

  trait.list <- unique(trait.list)
  snp.list <- unique(snp.list)

  obs.stat <- list()
  m <- length(snp.list)
  for(trait_id in 1:length(trait.list)) {
    # if one trait shows in multiple studies, add them up
    U.bytrait <- numeric(m)
    names(U.bytrait) <- snp.list
    V.bytrait <- matrix(0, m, m)
    colnames(V.bytrait) <- rownames(V.bytrait) <- snp.list

    for(study_id in 1:length(U)) {
      sublist_id <- which(colnames(U[[study_id]]) == trait.list[trait_id])
      if(length(sublist_id) != 0) {
        sublist <- U[[study_id]][, sublist_id]
        varU.sublist <- varU[[study_id]][, sublist_id]
        sublist_snps <- names(sublist)
        sublist_order <- order(match(sublist_snps, snp.list))
        list_order <- which(snp.list %in% sublist_snps)
        U.bytrait[list_order] <- U.bytrait[list_order] + sublist[sublist_order]
        R_order <- order(match(colnames(R)[colnames(R) %in% sublist_snps], snp.list))
        V.sublist <- diag(sqrt(varU.sublist[sublist_order])) %*% R[R_order, R_order] %*% diag(sqrt(varU.sublist[sublist_order]))
        colnames(V.sublist) <- rownames(V.sublist) <- names(varU.sublist)
        V.bytrait[list_order, list_order] <- V.bytrait[list_order, list_order] + V.sublist[sublist_order, sublist_order]
      }
    }
    obs.stat[[trait_id]] <- list(U = U.bytrait, V = V.bytrait)
  }
  names(obs.stat) <- trait.list
  ret <- convert_toUV(obs.stat)
  return(ret)
}

#' Compute the summary statistics given the genetic effect estimates and their standard errors
#'
#' This function allows you to calculate the score summary statistics \eqn{U} and their covariance matrix \eqn{V} for MTAR, given the genetic effect estimates and their standard errors
#' @param Beta a numeric list, each sublist containing estimation information of genetic effect estimates \eqn{\beta} for each study. In each study, a numeric \eqn{m \times K} matrix with each row as a SNP and each column as a separate trait. The number of traits and the number of SNPs in each study can be different but their names are required.
#' @param Beta.se a numeric list, each sublist containing the standard error of estimators information for each study. In each study, a numeric \eqn{m \times K} matrix with each row as a SNP and each column as a separate trait.
#' @param R a SNP correlation matrix, which should contain the correlation of all the SNPs in these studies.
#' @return A list containing summary statistics for each traits, including the score summary statistics \eqn{U} and their covariance matrix \eqn{V}.
#' @author Lan Luo
#' @export
#' @examples
#' data("beta.example")
#' attach(beta.example)
#' obs.stat <- Get_UV_from_beta(Beta = Beta, Beta.se = Beta.se, R = R)
#' detach(beta.example)


Get_UV_from_beta <- function(Beta, Beta.se, R) {

  if(length(Beta) != length(Beta.se)) {
    stop("The number of studies is not equal in Beta and Beta.se")
  }
  trait.list <- NULL
  snp.list <- NULL
  for(study_id in 1:length(Beta)) {
    if(is.null(colnames(Beta[[study_id]]))){
      stop(paste0("No input of names of traits in study ", study_id, ", calculation of summary statistics cannot continue."))
    }
    trait.list <- c(trait.list, colnames(Beta[[study_id]]))
    if(is.null(rownames(Beta[[study_id]]))){
      stop(paste0("No input of names of SNPs in study ", study_id, ", calculation of summary statistics cannot continue."))
    }
    snp.list <- c(snp.list, rownames(Beta[[study_id]]))
  }

  trait.list <- unique(trait.list)
  snp.list <- unique(snp.list)

  obs.stat <- list()
  m <- length(snp.list)
  for(trait_id in 1:length(trait.list)) {
    # if one trait shows in multiple studies, add them up
    U.bytrait <- numeric(m)
    names(U.bytrait) <- snp.list
    V.bytrait <- matrix(0, m, m)
    colnames(V.bytrait) <- rownames(V.bytrait) <- snp.list

    for(study_id in 1:length(Beta)) {
      sublist_id <- which(colnames(Beta[[study_id]]) == trait.list[trait_id])
      if(length(sublist_id) != 0) {
        sublist <- Beta[[study_id]][, sublist_id]
        se.sublist <- Beta.se[[study_id]][, sublist_id]
        sublist_snps <- names(sublist)
        sublist_order <- order(match(sublist_snps, snp.list))
        list_order <- which(snp.list %in% sublist_snps)

        U.var <- 1/(se.sublist)^2
        U <- sublist * U.var
        R_order <- order(match(colnames(R)[colnames(R) %in% sublist_snps], snp.list))
        V.sublist <- diag(sqrt(U.var[sublist_order])) %*% R[R_order, R_order] %*% diag(sqrt(U.var[sublist_order]))
        colnames(V.sublist) <- rownames(V.sublist) <- names(U.var)

        U.bytrait[list_order] <- U.bytrait[list_order] + U[sublist_order]
        V.bytrait[list_order, list_order] <- V.bytrait[list_order, list_order] + V.sublist[sublist_order, sublist_order]
      }
    }
    obs.stat[[trait_id]] <- list(U = U.bytrait, V = V.bytrait)
  }
  names(obs.stat) <- trait.list
  ret <- convert_toUV(obs.stat)
  return(ret)
}

#' Calculate Covariances of Z-scores between Traits from Overlapping Samples
#'
#' This function allows you to estimate the matrix \eqn{\zeta} to adjust for the potential sample overlap in the data set. Here we applied LD pruning (\eqn{r^2 < 0.1} in 500kb region) on 1000 Genome genotype dataset (hg19) as a list of reference independent SNPs. The SNP ID is chr:pos.
#' @param Zscore a numeric list, each sublist containing a vector of Z scores of SNPs with minor allele frequency (MAF) larger than 0.05. The chr:pos for each SNP is required.
#' @param pval_cutoff a numeric value indicating the cutoff threshold of p-values. The default value is 0.05. Variants with p-value less than or equal to this threshold will be automatically removed.
#' @param Indp_common_snp a numeric list of independent common SNPs
#' @return A \eqn{K \times K} matrix \eqn{\zeta}, where \eqn{K} is the number of traits.
#' @author Lan Luo
#' @export
#' @examples
#' \donttest{data(zeta.example)
#' attach(zeta.example)
#' # Downloading independent common SNPs from 1000Genome data set.
#' githubURL <- "https://github.com/lan/MTAR/blob/master/indp_snps.1KG.rda?raw=true"
#' utils::download.file(githubURL,"1kgfile")
#' load("1kgfile")
#' zeta1 <- Get_zeta(Zscore = Zscore, Indp_common_snp = indp_snps.1KG)
#' zeta1
#' detach(zeta.example)
#' }

Get_zeta <- function(Zscore, pval_cutoff = 0.05, Indp_common_snp){

  K <- length(Zscore)
  zscore <- list()
  for(k in 1:K) {
    pvalue <- 2 * pnorm(-abs(Zscore[[k]]))
    zscore[[k]] <- Zscore[[k]][names(Zscore[[k]]) %in% Indp_common_snp & pvalue > pval_cutoff]
  }

  subset_snps <- list()
  for(iter in 1:length(zscore)) {
    subset_snps[[iter]] <- names(zscore[[iter]])
  }
  Indp_null_common_snp <- Reduce(intersect, subset_snps)

  zscore.mat <- NULL
  for(k in 1:K) {
    zscore.mat <- cbind(zscore.mat, zscore[[k]][which(names(zscore[[k]]) %in% Indp_null_common_snp)])
  }
  if(nrow(zscore.mat) < 500){
    message(paste0("There are only ", nrow(zscore.mat), " common null independent SNPs remaining in calculating zeta, the result may be less accurate."))
  }else{
    message(paste0("There are ", nrow(zscore.mat), " common null independent SNPs remaining in calculating zeta."))
  }
  zeta <- cor(zscore.mat)
  colnames(zeta) <- rownames(zeta) <- names(Zscore)
  return(zeta)
}
convert_toUV <- function(obs.stat){
  K <- length(obs.stat)
  U <- list()
  V <- list()
  for(k in 1:K) {
    U[[k]] <- obs.stat[[k]]$U
    V[[k]] <- obs.stat[[k]]$V
  }
  names(U) <- names(V) <- names(obs.stat)
  ret <- list(U = U, V = V)
  return(ret)
}
