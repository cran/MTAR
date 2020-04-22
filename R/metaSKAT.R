# non-exported functions from SKAT package
# https://github.com/leeshawn/SKAT
# date: 01/23/2020

Get_Lambda<-function (K)
{
  out.s <- eigen(K, symmetric = TRUE, only.values = TRUE)
  lambda1 <- out.s$values
  IDX1 <- which(lambda1 >= 0)
  IDX2 <- which(lambda1 > mean(lambda1[IDX1])/1e+05)
  if (length(IDX2) == 0) {
    stop("No Eigenvalue is bigger than 0!!")
  }
  lambda <- lambda1[IDX2]
  return(lambda)
}
Get_Liu_PVal.MOD.Lambda <- function (Q.all, lambda, log.p = FALSE)
{
  param <- Get_Liu_Params_Mod_Lambda(lambda)
  Q.Norm <- (Q.all - param$muQ)/param$sigmaQ
  Q.Norm1 <- Q.Norm * param$sigmaX + param$muX
  p.value <- pchisq(Q.Norm1, df = param$l, ncp = param$d,
                    lower.tail = FALSE, log.p = log.p)
  return(p.value)
}

Get_Liu_Params_Mod_Lambda <- function (lambda)
{
  c1 <- rep(0, 4)
  for (i in 1:4) {
    c1[i] <- sum(lambda^i)
  }
  muQ <- c1[1]
  sigmaQ <- sqrt(2 * c1[2])
  s1 = c1[3]/c1[2]^(3/2)
  s2 = c1[4]/c1[2]^2
  beta1 <- sqrt(8) * s1
  beta2 <- 12 * s2
  type1 <- 0
  if (s1^2 > s2) {
    a = 1/(s1 - sqrt(s1^2 - s2))
    d = s1 * a^3 - a^2
    l = a^2 - 2 * d
  }
  else {
    type1 <- 1
    l = 1/s2
    a = sqrt(l)
    d = 0
  }
  muX <- l + d
  sigmaX <- sqrt(2) * a
  re <- list(l = l, d = d, muQ = muQ, muX = muX, sigmaQ = sigmaQ,
             sigmaX = sigmaX)
  return(re)
}

Get_Liu_PVal.MOD.Lambda.Zero <- function (Q, muQ, muX, sigmaQ, sigmaX, l, d)
{
  Q.Norm <- (Q - muQ)/sigmaQ
  Q.Norm1 <- Q.Norm * sigmaX + muX
  temp <- c(0.05, 10^-10, 10^-20, 10^-30, 10^-40, 10^-50,
            10^-60, 10^-70, 10^-80, 10^-90, 10^-100)
  out <- qchisq(temp, df = l, ncp = d, lower.tail = FALSE)
  IDX <- max(which(out < Q.Norm1))
  pval.msg <- sprintf("Pvalue < %e", temp[IDX])
  return(pval.msg)
}

Get_Liu_Params_Mod <- function (c1)
{
  muQ <- c1[1]
  sigmaQ <- sqrt(2 * c1[2])
  s1 = c1[3]/c1[2]^(3/2)
  s2 = c1[4]/c1[2]^2
  beta1 <- sqrt(8) * s1
  beta2 <- 12 * s2
  type1 <- 0
  if (s1^2 > s2) {
    a = 1/(s1 - sqrt(s1^2 - s2))
    d = s1 * a^3 - a^2
    l = a^2 - 2 * d
  }
  else {
    type1 <- 1
    l = 1/s2
    a = sqrt(l)
    d = 0
  }
  muX <- l + d
  sigmaX <- sqrt(2) * a
  re <- list(l = l, d = d, muQ = muQ, muX = muX, sigmaQ = sigmaQ,
             sigmaX = sigmaX)
  return(re)
}

SKAT_Optimal_Integrate_Func_Davies <- function (x, pmin.q, param.m, r.all)
{
  n.r <- length(r.all)
  n.x <- length(x)
  temp1 <- param.m$tau %x% t(x)
  temp <- (pmin.q - temp1)/(1 - r.all)
  temp.min <- apply(temp, 2, min)
  re <- rep(0, length(x))
  for (i in 1:length(x)) {
    min1 <- temp.min[i]
    if (min1 > sum(param.m$lambda) * 10^4) {
      temp <- 0
    }
    else {
      min1.temp <- min1 - param.m$MuQ
      sd1 <- sqrt(param.m$VarQ - param.m$VarRemain)/sqrt(param.m$VarQ)
      min1.st <- min1.temp * sd1 + param.m$MuQ
      dav.re <- CompQuadForm::davies(min1.st, param.m$lambda, acc = 10^(-6))
      temp <- dav.re$Qq
      if (dav.re$ifault != 0) {
        stop("dav.re$ifault is not 0")
      }
    }
    if (temp > 1) {
      temp = 1
    }
    re[i] <- (1 - temp) * dchisq(x[i], df = 1)
  }
  return(re)
}
SKAT_Optimal_PValue_Liu <- function (pmin.q, param.m, r.all, pmin = NULL)
{
  re <- integrate(SKAT_Optimal_Integrate_Func_Liu, lower = 0,
                  upper = 40, subdivisions = 2000, pmin.q = pmin.q, param.m = param.m,
                  r.all = r.all, abs.tol = 10^-25)
  pvalue <- 1 - re[[1]]
  if (!is.null(pmin)) {
    if (pmin * length(r.all) < pvalue) {
      pvalue = pmin * length(r.all)
    }
  }
  return(pvalue)
}
SKAT_Optimal_Integrate_Func_Liu <- function (x, pmin.q, param.m, r.all)
{
  n.r <- length(r.all)
  n.x <- length(x)
  temp1 <- param.m$tau %x% t(x)
  temp <- (pmin.q - temp1)/(1 - r.all)
  temp.min <- apply(temp, 2, min)
  temp.q <- (temp.min - param.m$MuQ)/sqrt(param.m$VarQ) *
    sqrt(2 * param.m$Df) + param.m$Df
  re <- pchisq(temp.q, df = param.m$Df) * dchisq(x, df = 1)
  return(re)
}

Get_Liu_Params <- function (c1)
{
  muQ <- c1[1]
  sigmaQ <- sqrt(2 * c1[2])
  s1 = c1[3]/c1[2]^(3/2)
  s2 = c1[4]/c1[2]^2
  beta1 <- sqrt(8) * s1
  beta2 <- 12 * s2
  type1 <- 0
  if (s1^2 > s2) {
    a = 1/(s1 - sqrt(s1^2 - s2))
    d = s1 * a^3 - a^2
    l = a^2 - 2 * d
  }
  else {
    type1 <- 1
    a = 1/s1
    d = 0
    l = 1/s1^2
  }
  muX <- l + d
  sigmaX <- sqrt(2) * a
  re <- list(l = l, d = d, muQ = muQ, muX = muX, sigmaQ = sigmaQ,
             sigmaX = sigmaX)
  return(re)
}

Get_Liu_PVal <- function (Q, W, Q.resampling = NULL)
{
  Q.all <- c(Q, Q.resampling)
  A1 <- W/2
  A2 <- A1 %*% A1
  c1 <- rep(0, 4)
  c1[1] <- sum(diag(A1))
  c1[2] <- sum(diag(A2))
  c1[3] <- sum(A1 * t(A2))
  c1[4] <- sum(A2 * t(A2))
  param <- Get_Liu_Params(c1)
  Q.Norm <- (Q.all - param$muQ)/param$sigmaQ
  Q.Norm1 <- Q.Norm * param$sigmaX + param$muX
  p.value <- pchisq(Q.Norm1, df = param$l, ncp = param$d,
                    lower.tail = FALSE)
  p.value.resampling = NULL
  if (length(Q.resampling) > 0) {
    p.value.resampling <- p.value[-1]
  }
  re <- list(p.value = p.value[1], param = param, p.value.resampling = p.value.resampling)
  return(re)
}

Get_PValue.Lambda <- function (lambda, Q)
{
  n1 <- length(Q)
  p.val <- rep(0, n1)
  p.val.liu <- rep(0, n1)
  is_converge <- rep(0, n1)
  p.val.liu <- Get_Liu_PVal.MOD.Lambda(Q, lambda)
  for (i in 1:n1) {
    out <- CompQuadForm::davies(Q[i], lambda, acc = 10^(-6))
    p.val[i] <- out$Qq
    is_converge[i] <- 1
    if (length(lambda) == 1) {
      p.val[i] <- p.val.liu[i]
    }
    else if (out$ifault != 0) {
      is_converge[i] <- 0
    }
    if (p.val[i] > 1 || p.val[i] <= 0) {
      is_converge[i] <- 0
      p.val[i] <- p.val.liu[i]
    }
  }
  p.val.msg = NULL
  p.val.log = NULL
  if (p.val[1] == 0) {
    param <- Get_Liu_Params_Mod_Lambda(lambda)
    p.val.msg <- Get_Liu_PVal.MOD.Lambda.Zero(Q[1], param$muQ,
                                              param$muX, param$sigmaQ, param$sigmaX, param$l,
                                              param$d)
    p.val.log <- Get_Liu_PVal.MOD.Lambda(Q[1], lambda, log.p = TRUE)[1]
  }
  return(list(p.value = p.val, p.val.liu = p.val.liu, is_converge = is_converge,
              p.val.log = p.val.log, pval.zero.msg = p.val.msg))
}

SKAT_Optimal_Each_Q <- function (param.m, Q.all, r.all, lambda.all, method = NULL)
{
  n.r <- length(r.all)
  c1 <- rep(0, 4)
  n.q <- dim(Q.all)[1]
  pval <- matrix(rep(0, n.r * n.q), ncol = n.r)
  pmin.q <- matrix(rep(0, n.r * n.q), ncol = n.r)
  param.mat <- NULL
  for (i in 1:n.r) {
    Q <- Q.all[, i]
    r.corr <- r.all[i]
    lambda.temp <- lambda.all[[i]]
    c1[1] <- sum(lambda.temp)
    c1[2] <- sum(lambda.temp^2)
    c1[3] <- sum(lambda.temp^3)
    c1[4] <- sum(lambda.temp^4)
    param.temp <- Get_Liu_Params_Mod(c1)
    muQ <- param.temp$muQ
    varQ <- param.temp$sigmaQ^2
    df <- param.temp$l
    Q.Norm <- (Q - muQ)/sqrt(varQ) * sqrt(2 * df) + df
    pval[, i] <- pchisq(Q.Norm, df = df, lower.tail = FALSE)
    if (!is.null(method)) {
      if (method == "optimal.mod" || method == "optimal.adj" ||
          method == "optimal.moment.adj") {
        pval[, i] <- Get_PValue.Lambda(lambda.temp,
                                       Q)$p.value
      }
    }
    param.mat <- rbind(param.mat, c(muQ, varQ, df))
  }
  pmin <- apply(pval, 1, min)
  for (i in 1:n.r) {
    muQ <- param.mat[i, 1]
    varQ <- param.mat[i, 2]
    df <- param.mat[i, 3]
    q.org <- qchisq(1 - pmin, df = df)
    q.q <- (q.org - df)/sqrt(2 * df) * sqrt(varQ) + muQ
    pmin.q[, i] <- q.q
  }
  out <- list(pmin = pmin, pval = pval, pmin.q = pmin.q)
  return(out)
}

SKAT_Optimal_PValue_Davies <- function (pmin.q, param.m, r.all, pmin = NULL)
{
  re <- try(integrate(SKAT_Optimal_Integrate_Func_Davies,
                      lower = 0, upper = 40, subdivisions = 1000, pmin.q = pmin.q,
                      param.m = param.m, r.all = r.all, abs.tol = 10^-25),
            silent = TRUE)
  if (class(re) == "try-error") {
    re <- SKAT_Optimal_PValue_Liu(pmin.q, param.m, r.all,
                                  pmin)
    return(re)
  }
  pvalue <- 1 - re[[1]]
  if (!is.null(pmin)) {
    if (pmin * length(r.all) < pvalue) {
      pvalue = pmin * length(r.all)
    }
  }
  return(pvalue)
}

SKAT_META_Optimal_Get_Q<-function(Score, r.all){
  # no change
  n.r<-length(r.all)
  Q.r<-rep(0,n.r)

  for(i in 1:n.r){
    r.corr<-r.all[i]
    Q.r[i]<-(1-r.corr) * sum(Score^2) + r.corr * sum(Score)^2
  }
  Q.r = Q.r /2

  re<-list(Q.r=Q.r)
  return(re)
}
SKAT_META_Optimal_Param<-function(Phi,r.all){
  # no change
  p.m<-dim(Phi)[2]
  r.n<-length(r.all)

  # ZMZ
  Z.item1.1<- Phi %*% rep(1,p.m)
  ZZ<-Phi
  ZMZ<- Z.item1.1 %*% t(Z.item1.1) / sum(ZZ)

  # W3.2 Term : mixture chisq
  W3.2.t<-ZZ - ZMZ
  lambda<-Get_Lambda(W3.2.t)

  # W3.3 Term : variance of remaining ...
  W3.3.item<-sum(ZMZ *(ZZ-ZMZ)) * 4

  # tau term
  z_mean_2<- sum(ZZ)/p.m^2
  tau1<- sum(ZZ %*% ZZ) / p.m^2 / z_mean_2

  # Mixture Parameters
  MuQ<-sum(lambda)
  VarQ<-sum(lambda^2) *2 + W3.3.item
  KerQ<-sum(lambda^4)/(sum(lambda^2))^2 * 12
  Df<-12/KerQ

  # W3.1 Term : tau1 * chisq_1
  tau<-rep(0,r.n)
  for(i in 1:r.n){
    r.corr<-r.all[i]
    term1<-p.m^2*r.corr * z_mean_2 + tau1 * (1-r.corr)
    tau[i]<-term1
  }

  out<-list(MuQ=MuQ,VarQ=VarQ,KerQ=KerQ,lambda=lambda,VarRemain=W3.3.item,Df=Df,tau=tau,
    z_mean_2=z_mean_2, p.m=p.m,
    tau.1 = tau1,
    tau.2= p.m*z_mean_2 )

  #param2<<-out
  return(out)
}
SKAT_META_Optimal_Get_Pvalue<-function(Q.all, Phi, r.all, method){
  # no change
  n.r<-length(r.all)
  n.q<-dim(Q.all)[1]
  p.m<-dim(Phi)[2]

  lambda.all<-list()
  for(i in 1:n.r){
    r.corr<-r.all[i]
    R.M<-diag(rep(1-r.corr,p.m)) + matrix(rep(r.corr,p.m*p.m),ncol=p.m)
    L<-chol(R.M,pivot=TRUE)
    Phi_rho<- L %*% (Phi %*% t(L))
    lambda.all[[i]]<-Get_Lambda(Phi_rho)
  }

  # Get Mixture param
  param.m <- SKAT_META_Optimal_Param(Phi,r.all)
  Each_Info <- SKAT_Optimal_Each_Q(param.m, Q.all, r.all,
    lambda.all, method=method)
  pmin.q<-Each_Info$pmin.q
  pval <- rep(0,n.q)

  # added
  pmin<-Each_Info$pmin

  for(i in 1:n.q){
    pval[i]<-SKAT_Optimal_PValue_Davies(pmin.q[i,],param.m,r.all, pmin[i])
  }

  # Check the pval
  # Since SKAT-O is between burden and SKAT, SKAT-O p-value should be <= min(p-values) * 2
  # To correct conservatively, we use min(p-values) * 3

  multi<-3
  if(length(r.all) < 3){
    multi<-2
  }

  for(i in 1:n.q){
    pval.each<-Each_Info$pval[i,]
    IDX<-which(pval.each > 0)

    pval1<-min(pval.each) * multi
    if(pval[i] <= 0 || length(IDX) < length(r.all)){
      pval[i]<-pval1
    }
    # if pval==0, use nonzero min each.pval as p-value
    if(pval[i] == 0){
      if(length(IDX) > 0){
        pval[i] = min(pval.each[IDX])
      }
    }

  }

  return(list(p.value=pval,p.val.each=Each_Info$pval))
}
SKAT_META_Optimal <- function(Score, Phi, r.all, method="davies"){
  # no change
  # if r.all >=0.999 ,then r.all = 0.999
  IDX<-which(r.all >= 0.999)
  if(length(IDX) > 0){
    r.all[IDX]<-0.999
  }

  p.m<-dim(Phi)[2]
  n.r<-length(r.all)

  ###########################################
  # Compute Q.r and Q.r.res
  ##########################################
  out.Q <- SKAT_META_Optimal_Get_Q(Score, r.all)
  Q.res=NULL
  Q.all<-rbind(out.Q$Q.r, Q.res)

  ##################################################
  # Compute P-values
  #################################################

  out<-SKAT_META_Optimal_Get_Pvalue(Q.all, Phi/2, r.all, method)

  param<-list(p.val.each=NULL, q.val.each=NULL)
  param$p.val.each<-out$p.val.each[1,]
  param$q.val.each<-Q.all[1,]
  param$rho<-r.all
  param$minp<-min(param$p.val.each)

  id_temp<-which(param$p.val.each == min(param$p.val.each))
  id_temp1<-which(param$rho >= 0.999) # treat rho > 0.999 as 1
  if(length(id_temp1) > 0){
    param$rho[id_temp1] = 1
  }

  param$rho_est<-param$rho[id_temp]
  p.value<-out$p.value[1]
  re<-list(p.value = p.value, param=param)
  return(re)
}
Met_SKAT_Get_Pvalue<-function(Score, Phi, r.corr = 0){
  method <- "davies"
  # change SKAT
  Q.res = NULL
  p.m<-nrow(Phi)
  # if Phi==0
  if(sum(abs(Phi)) == 0){
    warning("No polymorphic SNPs!",call.=FALSE)
    return(list(p.value=1, p.value.resampling= NULL, pval.zero.msg=NULL))
  }

  if(length(Phi) <=1){
    r.corr=0
  } else{
    if(ncol(Phi) <=10){
      if(qr(Phi)$rank <= 1){
        r.corr=0
      }
    }
  }

  if(length(r.corr) > 1){
    re = SKAT_META_Optimal(Score, Phi, r.corr, method=method)
    return(re)
  }

  if (r.corr == 0){
    Q<-sum(Score^2)/2
  } else if (r.corr==1){
    Q <- SKAT_META_Optimal_Get_Q(Score, r.corr)$Q.r
    a<- as.matrix(sum(Phi))
    re <- Get_Liu_PVal(Q, a, Q.res)
    return(re)
  } else {
    # like r.corr = 0.1 or 0.2
    Q <- SKAT_META_Optimal_Get_Q(Score, r.corr)$Q.r
    R.M<-diag(rep(1-r.corr,p.m)) + matrix(rep(r.corr,p.m*p.m),ncol=p.m)
    L<-chol(R.M,pivot=TRUE)
    Phi<- L %*% (Phi %*% t(L))
  }
  lambda <- Get_Lambda(Phi/2)
  re1 <- pchisqsum2(Q, lambda, method = "integration", acc=1e-20)$p
  if(re1 <= 0 | re1 >= 1){
    re1 <- Get_PValue.Lambda(lambda, Q)$p.value
  }
  re <- list(p.value = re1)
  # re<-SKAT:::Get_Davies_PVal(Q, Phi)$p.value #result is same as SKAT:::Get_PValue.Lambda
  return(re)
}

# non-exported functions from SeqMeta package
# https://github.com/cran/seqMeta
# date: 04/22/2020

pchisqsum2 <- function (Q, lambda, delta = rep(0, length(lambda)), method = c("saddlepoint",
                                                                              "integration", "liu"), acc = 1e-07)
{
  method <- match.arg(method)
  delta <- delta[lambda > 0]
  lambda <- lambda[lambda > 0]
  if (method == "saddlepoint") {
    p = saddle(Q, lambda, delta)
    if (is.na(p)) {
      method <- "integration"
    }
    else {
      return(list(p = p, errflag = 0))
    }
  }
  if (method == "integration") {
    tmp <- CompQuadForm::davies(q = Q, lambda = lambda,
                                delta = delta, acc = acc)
    if (tmp$ifault > 0) {
      lambda <- zapsmall(lambda, digits = 2)
      delta <- delta[lambda > 0]
      lambda <- lambda[lambda > 0]
      tmp <- CompQuadForm::farebrother(q = Q, lambda = lambda,
                                       delta = delta)
    }
    Qq <- if ("Qq" %in% names(tmp))
      tmp$Qq
    else tmp$res
    return(list(p = Qq, errflag = 0))
  }
  if (method == "liu") {
    tmp <- CompQuadForm::liu(Q, lambda = lambda, delta = delta)
    return(list(p = tmp, errflag = 0))
  }
}

