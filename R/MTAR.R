#' Multiple-Traits Analysis of Rare-Variant Association Test
#'
#' Test for association between a set of rare SNPs and multiple traits with input of summary statistics, possibly from overlap samples. The input number of SNPs in each trait can be different, MTAR function will analyze the union of SNPs that show up in at least one trait, and automatically handle the non-polymorphic SNPs.
#' @param U a numeric list, each sublist containing summary statistics U for each traits. The SNP IDs must be provided.
#' @param V a numeric list, each sublist containing the corresponding covariance matrix of summary statistics. If your original summary statistics are other format, please use Get_UV_from_data, Get_UV_from_varU or Get_UV_from_beta to generate the summary statistics \eqn{U} and \eqn{V} for MTAR.
#' @param MAF a numeric vector containing minor allele frequency for the SNPs show up in at least one trait. The SNP IDs must be provided.
#' @param MAF_UB a numeric value indicating the cutoff threshold of minor allele frequency for SNPs The default value is 0.05.
#' @param zeta a numeric matrix containing the sample correlation of Z-scores over a large number of independent null common SNPs across genome. The default value is NULL, where MTAR assumes there are no overlap samples. However, if there is overlapping in subjects, zeta must be provided. zeta can be estimated using MTAR::Get_zeta.
#' @param genetic_cor.trait a numeric matrix containing the genetic correlation among traits. The default value of genetic_cor.trait is NULL, where an exchangeable correlation structure with the correlation coefficient denoted by rho.trait (\eqn{\rho_2}) is assumed. In this case, there is no difference between cMTAR and iMTAR.
#' @param rho.SNP a numeric vector containing all the possible values of \eqn{\rho_1}. The default value is c(0, 0.5, 1).
#' @param rho.trait a numeric vector containing all the possible values of \eqn{\rho_2}. The default value is c(0, 0.5, 1).
#' @param weight.SNP a numeric vector containing the parameters in Beta density function to calculate the weight among SNPs. The default value is c(1, 25).
#' @details MTAR assumes that the genetic effect estimates \eqn{\beta} has covariance matrix \eqn{B}, which is a Kronecker product of two pre-specified matrices: among-variant effect covariance \eqn{B_1} and among-trait effect covariance \eqn{B_2}. An exchaneable correlation structure with the correlation coefficient denoted by rho.SNPs (\eqn{\rho_1}) for \eqn{B_1} is assumed. The default MTAR requires the input of genetic correlation matrix genetoc_cor.trait, if missing, then an exchaneable correlation structure for rho.trait (\eqn{\rho_2}) is assumed. The default weight of \eqn{B_1} is \eqn{dBeta(MAF, 1, 25)}, which can be changed freely by users.
#' @return a list of p-values of MTAR-O, cMTAR, iMTAR and cctP as well as ancillary information. Here cctP is the Cauchy-combined p-value of SKAT and burden tests with default weight \eqn{dBeta(MAF, 1, 25)}.
#' @references Liu, Y., Chen, S., Li, Z., Morrison, A.C., Boerwinkle, E. and Lin, X., 2019. ACAT: A fast and powerful p value combination method for rare-variant analysis in sequencing studies. The American Journal of Human Genetics, 104(3), pp.410-421.
#' @references Liu, Y. and Xie, J., 2019. Cauchy combination test: a powerful test with analytic p-value calculation under arbitrary dependency structures. Journal of the American Statistical Association, pp.1-18.
#' @references Luo, L., Shen, J., Zhang, H., Chhibber, A., Mehrotra, D. V., Tang, Z., 2019. Multi-trait analysis of rare-variant association summary statistics using MTAR.
#' @author Lan Luo
#' @export
#' @examples
#' data(MTAR.example)
#' attach(MTAR.example)
#' pval <- MTAR(U = U, V = V, MAF = MAF, genetic_cor.trait = genetic_cor.trait,
#'              zeta = zeta)
#' pval
#' detach(MTAR.example)
#'
MTAR <- function(U, V, MAF, MAF_UB = 0.05, zeta = NULL,
                 genetic_cor.trait = NULL,
                 rho.SNP = c(0, 0.5, 1), rho.trait = c(0, 0.5, 1),
                 weight.SNP = c(1, 25)) {
    KA <- genetic_cor.trait
    K <- length(U)

    # MAF filtering
    if(any(MAF <= MAF_UB)) {
      MAF <- MAF[MAF <= MAF_UB]
      U1 <- list()
      V1 <- list()
      for(k in 1:K) {
        U1[[k]] <- U[[k]][names(U[[k]]) %in% names(MAF)]
        V1[[k]] <- V[[k]][rownames(V[[k]]) %in% names(MAF),
                          colnames(V[[k]]) %in% names(MAF)]
      }
      names(U1) <- names(V1) <- names(U)
      U <- U1
      V <- V1
      message(paste0("There are ", length(MAF), " rare variants with MAF less or equal to ", MAF_UB))
    }

    if(is.null(names(unlist(U))) | is.null(names(MAF))) {
      stop("SNPs IDs are required for U, V, and MAF.")
    }
    snp.list <- NULL
    for(k in 1:K) {
      snp.list <- c(snp.list, names(U[[k]]))
    }
    snp.list <- unique(snp.list)

    if(!all(snp.list %in% names(MAF))) {
      stop("The minor allele frequencies of some SNPs are missing, please provide the MAF for all SNPs.")
    }
    SKAT.p <- numeric(K)
    Burden.p <- numeric(K)
    for (k in 1:K) {
      order1 <- order(match(names(U[[k]])[names(U[[k]]) %in% names(MAF)], names(MAF)))
      order2 <- which(names(MAF) %in% names(U[[k]]))
      SKAT.p[k] <- RVAS(U = U[[k]][order1], V = V[[k]][order1,order1], MAF = MAF[order2], weight.beta = c(1, 25), type = "variance")
      Burden.p[k] <- RVAS(U = U[[k]][order1], V = V[[k]][order1,order1], MAF = MAF[order2], weight.beta = c(1, 25), type = "mean")
    }
    cMTAR <- MTAR.main(U = U, V = V, MAF = MAF, snp.list = snp.list, KA = KA, zeta = zeta, iMTAR = FALSE, rho.SNP = rho.SNP, rho.trait = rho.trait, weight.SNP = weight.SNP)
    iMTAR <- MTAR.main(U = U, V = V, MAF = MAF, snp.list = snp.list, KA = KA, zeta = zeta, iMTAR = TRUE, rho.SNP = rho.SNP, rho.trait = rho.trait, weight.SNP = weight.SNP)
    cctP <- ACAT(c(SKAT.p, Burden.p))
    p.list <- c(cMTAR$p, iMTAR$p, cctP)
    MTARO <- ACAT(p.list[!is.na(p.list)])
    pval <- list(MTARO = MTARO,
                 cMTAR = cMTAR, iMTAR = iMTAR,
                 cctP = cctP)
  return(pval)
}

MTAR.main <- function(U, V, MAF, KA = NULL, snp.list, zeta = NULL, iMTAR = FALSE, rho.SNP, rho.trait, weight.SNP){

  lambda <- expand.grid(lambdaA = 1 - rho.trait,
                        lambdaC = rho.SNP)

  K <- length(U)
  m <- length(MAF)

  if(iMTAR) {
    JA <- diag(K)
  } else{
    JA <- matrix(1, nrow = K, ncol = K)
  }

  if (m < 2) {
    stop("MTAR requires more than one SNP as input.")
  }
  if (is.null(KA)) {
    message("Without the input of genetic correlation information, an exchangeable correlation structure among traits is assumed in MTAR.")
    KA <- diag(K)
    JA <- matrix(1, K, K)
  }
  if (is.null(zeta)) {
    message("Without the input of zeta, MTAR assumes there are no overlap samples in the dataset.")
  }

  R.bytraits <- list()
  for(k in 1:K) {
    R.bytraits[[k]] <- cov2cor(V[[k]])
  }

  R.C <- matrix(0, length(snp.list), length(snp.list))
  diag(R.C) <- 1
  colnames(R.C) <- rownames(R.C) <- snp.list
  row_col.id <- which(lower.tri(R.C), arr.ind = TRUE)

  for(iter in 1:nrow(row_col.id)) {
    R.val <- 0
    row_col.name <- c(rownames(R.C)[row_col.id[iter, 1]], colnames(R.C)[row_col.id[iter, 2]])
    for(k in 1:K) {
      row.id <- which(rownames(R.bytraits[[k]]) == row_col.name[1])
      col.id <- which(colnames(R.bytraits[[k]]) == row_col.name[2])
      if(length(row.id) != 0 & length(col.id) != 0) {
        R.C[row_col.id[iter, 1], row_col.id[iter, 2]] <- R.bytraits[[k]][row.id, col.id]
        break
      }
    }
  }
  R.C <- Matrix::forceSymmetric(R.C, uplo="L")

  obs.stat <- list()
  U.complete <- list()
  for (k in 1:K) {
    U.bytrait <- numeric(length(snp.list))
    V.bytrait <- matrix(0, length(snp.list), length(snp.list))
    names(U.bytrait) <- snp.list
    colnames(V.bytrait) <- rownames(V.bytrait) <- snp.list
    order1 <- order(match(names(U[[k]])[names(U[[k]]) %in% snp.list], snp.list))
    order2 <- which(snp.list %in% names(U[[k]]))
    U.bytrait[order2] <- U[[k]][order1]
    V.bytrait[order2, order2] <- V[[k]][order1, order1]
    obs.stat[[k]] <- list(U = U.bytrait, V = V.bytrait)
    U.complete[[k]] <- U.bytrait
  }
  nonpolymorphic <- FALSE
  if (any(unlist(U.complete) == 0)) {
    update.obs.stat <- list()

    for (k in 1:K) {
      ind <- which(obs.stat[[k]]$U == 0)
      if (length(ind) != 0 ) {
        update.obs.stat[[k]] <- list(U = obs.stat[[k]]$U[-ind],
                                     V = obs.stat[[k]]$V[-ind, -ind, drop = FALSE])
      }else{
        update.obs.stat[[k]] <- list(U = obs.stat[[k]]$U,
                                     V = obs.stat[[k]]$V)
      }
    }

    update.R.C <- nonpolymorphic.fn(R.C, obs.stat)
    nonpolymorphic <- TRUE
  }

  if(nonpolymorphic) {
    obs.stat1 <- update.obs.stat
  } else{
    obs.stat1 <- obs.stat
  }

  if(!is.null(zeta)) {
    ## without calculating the invser of each V
    V.sqrt <- list()
    for(k in 1:K) {
      V.sqrt[[k]] <- diag(length(obs.stat1[[k]]$U))
      diag(V.sqrt[[k]]) <- sqrt(diag(obs.stat1[[k]]$V))
    }
    ind <- expand.grid(row.ind = 1:K, column.ind = 1:K)
    Ucov <- list()
    for(iter in 1:nrow(ind)) {
      k1 <- ind[iter, 1]
      k2 <- ind[iter, 2]
      if(nonpolymorphic) {
        Ucov[[iter]] <- V.sqrt[[k1]] %*% update.R.C$update.mat[[iter]] %*% V.sqrt[[k2]]
      }else{
        Ucov[[iter]] <- V.sqrt[[k1]] %*% R.C %*% V.sqrt[[k2]]
      }
    }

    U.cov <- NULL
    mlist <- list()
    for(col_id in 1:K) {
      diag.ind <- which(ind$column.ind == col_id & ind$row.ind == col_id)
      U.cov.col <- NULL
      index <- which(ind$column.ind == col_id)
      row_id <- ind[index, 1]
      for(iter in 1:length(index)) {
        # print(c(row_id[iter], col_id))
        U.cov.col <- rbind(U.cov.col, as.matrix(zeta[row_id[iter], col_id] * Ucov[[index[iter]]]))
      }
      U.cov <- cbind(U.cov, U.cov.col)
      mlist[[col_id]] <- Ucov[[diag.ind]]
    }
    U.cov <- as.matrix(U.cov)
    U.inv <- try(MASS::ginv(U.cov), silent = T)
    if(inherits(U.inv, "try-error")) {
      warning("The covariance matrix of U is exactly singular, MASS:ginv() function doesn't work here.")
      MTAR.cct.p <- list(p = NA,rho1.min = NA, rho2.min = NA)
      return(MTAR.cct.p)
    }
    V.diag <- as.matrix(Matrix::bdiag(mlist))
    Sigma.inv <- V.diag %*% MASS::ginv(U.cov) %*% V.diag
  } else {
    mlist <- list()
    for(k in 1:K){
      mlist[[k]] <- obs.stat1[[k]]$V
    }
    U.cov <- V.diag <- Sigma.inv <- as.matrix(Matrix::bdiag(mlist))
    U.inv <- try(MASS::ginv(U.cov), silent = T)
    if(inherits(U.inv, "try-error")) {
      warning("The covariance matrix of U is exactly singular, MASS:ginv() function doesn't work here.")
      MTAR.cct.p <- list(p = NA,rho1.min = NA, rho2.min = NA)
      return(MTAR.cct.p)
    }
  }
  U.alltraits <- NULL
  for(k in 1:K) {
    U.alltraits <- c(U.alltraits, obs.stat1[[k]]$U)
  }
  U.alltraits <- as.matrix(U.alltraits)

  KC <- diag(m)
  MAF <- MAF[order(match(names(MAF), snp.list))]
  diag(KC) <- Beta.Weights(MAF, weights.beta = weight.SNP)
  JC <- matrix(1, nrow = m, ncol = m)

  p.obs <- NULL
  for(ii in 1:nrow(lambda)) {
    lambdaA <- lambda[ii, 1]
    A <- (1 - lambdaA) * KA + lambdaA * JA

    lambdaC <- lambda[ii, 2]
    C <- KC %*% ((1 - lambdaC) * diag(m) + lambdaC * JC) %*% KC
    if(nonpolymorphic) {
      update.C <- nonpolymorphic.fn(C, obs.stat)$update.mat
      ind <- expand.grid(row.ind = 1:K, column.ind = 1:K)
      B <- NULL
      for(col_id in 1:K) {
        B.col <- NULL
        index <- which(ind$column.ind == col_id)
        row_id <- ind[index, 1]
        for(iter in 1:length(index)) {
          # print(c(row_id[iter], col_id))
          B.col <- rbind(B.col, A[row_id[iter], col_id] * update.C[[index[iter]]])
        }
        B <- cbind(B, B.col)
      }
    } else{
      B <- kronecker(A, C)
    }

    R <- chol(B, pivot = TRUE)
    r <- attr(R, 'rank')
    if (r < nrow(B)) R[(r+1):nrow(B), (r+1):nrow(B)] <- 0
    oo <- order(attr(R, 'pivot'))
    unpivQ <- R[, oo]
    rho <- try(Get_Lambda(unpivQ %*% Sigma.inv %*% t(unpivQ)), silent = TRUE)

    if(!inherits(rho, "try-error")) {
      Q <- as.numeric(t(U.alltraits) %*% MASS::ginv(U.cov) %*% V.diag %*% B %*% V.diag %*% MASS::ginv(U.cov) %*% U.alltraits)
      pvalues <- pchisqsum2(Q = Q, lambda=rho, method = c("integration"), acc=1e-20)$p
      if(!(pvalues < 0.99 & pvalues > 0)){
        pvalues <- Get_PValue.Lambda(rho, Q)$p.value
      }
      p.obs <- c(p.obs, pvalues)
    }
  }
  p.obs <- ifelse(p.obs == 1, 0.999, p.obs)

  # MTAR.cct.p <- ACAT(p.obs)
  MTAR.cct.p <- list(p = ACAT(p.obs),
                     rho1.min = lambda[which.min(p.obs),2], rho2.min = 1 - lambda[which.min(p.obs),1])
  return(MTAR.cct.p)

}

Beta.Weights<-function(MAF, weights.beta = c(1, 25), Cutoff=1, Is.MAF=TRUE){
  n<-length(MAF)
  weights<-rep(0,n)
  Sign<-rep(1,n)

  IDX1<-which(MAF > 0.5)
  if(length(IDX1) > 0){
    Sign[IDX1]<--1
    MAF[IDX1]<-1-MAF[IDX1]
  }

  IDX_0<-union(which(MAF == 0), which(MAF > Cutoff))
  if(length(IDX_0) == n){
    #stop("No polymorphic SNPs")
    weights<-rep(0,n)
  } else if( length(IDX_0) == 0){
    weights<-dbeta(MAF,weights.beta[1],weights.beta[2])
  } else {
    weights[-IDX_0]<-dbeta(MAF[-IDX_0],weights.beta[1],weights.beta[2])
  }

  weights = weights * Sign
  return(weights)

}

nonpolymorphic.fn <- function(mat, obs.stat) {
  K <- length(obs.stat)
  mat1 <- list()
  ind <- expand.grid(row.ind = 1:K, column.ind = 1:K)
  for(iter in 1:nrow(ind)) {
    k1 <- ind[iter, 1]
    k2 <- ind[iter, 2]
    ind1 <- which(obs.stat[[k1]]$U == 0)
    ind2 <- which(obs.stat[[k2]]$U == 0)
    if(length(ind1) != 0 & length(ind2) != 0) {
      mat1[[iter]] <- mat[-ind1, -ind2, drop = FALSE]
    } else if(length(ind1) != 0 & length(ind2) == 0) {
      mat1[[iter]] <- mat[-ind1, , drop = FALSE]
    } else if(length(ind1) == 0 & length(ind2) != 0) {
      mat1[[iter]] <- mat[, -ind2, drop = FALSE]
    } else {
      mat1[[iter]] <- mat
    }
  }
  update.ret <- list(update.mat = mat1, row_col.info = ind)
  return(update.ret)
}


RVAS <- function(U, V, MAF, weight.beta = NULL, type = "mean", grid = NULL){
  if(is.null(weight.beta)) {
    weight <- rep(1/length(MAF), length(MAF))
  } else {
    weight <- Beta.Weights(MAF, weights.beta = weight.beta)
  }
  Score <- U * weight
  SMat.Summary <- t(t(V * weight) * weight)
  if(type == "mean") {
    # run Burden
    ret <- Met_SKAT_Get_Pvalue(Score, SMat.Summary, r.corr = 1)$p.value
  }else if(type == "variance") {
    # run SKAT
    ret <- Met_SKAT_Get_Pvalue(Score, SMat.Summary, r.corr = 0)$p.value
  }else if(type == "optimal"){
    if(is.null(grid)) {
      grid <- seq(0, 1, by = 0.1)
    }
    # run SKATO
    ret <- Met_SKAT_Get_Pvalue(Score, SMat.Summary, r.corr = grid)$p.value
  }
  ret <- ifelse(ret == 1, 0.999, ret)
  return(ret)
}

