allocate.lspoisregmem <- function(data, beta0, n, p, q, testtype) {
  tt <- match(testtype, c("lrt", "score", "robustscore", "wald", "robustwald"))
  if (is.na(tt) == TRUE)
    return (NA)
  if (tt == 1) {
    return(allocate.lspoisregmem.lrt(data = data,
                                    beta0 = beta0,
                                    n = n,
                                    p = p,
                                    q = q))
    
  }
  if (tt == 2) {
    return (allocate.lspoisregmem.score(data = data,
                                       beta0 = beta0,
                                       n = n,
                                       p = p,
                                       q = q))
  }
  if (tt == 3) {
    return (allocate.lspoisregmem.robustscore(data = data,
                                             beta0 = beta0,
                                             n = n,
                                             p = p,
                                             q = q))
  }
  if (tt == 4) {
    return (allocate.lspoisregmem.wald(data = data,
                                      beta0 = beta0,
                                      n = n,
                                      p = p,
                                      q = q))
  }
  if (tt == 5) {
    return (allocate.lspoisregmem.robustwald(data = data,
                                            beta0 = beta0,
                                            n = n,
                                            p = p,
                                            q = q))
  }
  return(NA)
}

##############################################################################
####################     Fit                ##################################
##############################################################################

allocate.lspoisregmem.fit <- function(data, beta0, n, p, q) {
  y <- data[,1]
  xl <- as.matrix(data[,1:p])
  xl[,1] <- 1.  
  ql0 <- matrix(0, n, p)
  rtl0 <- matrix(0, p, p)
  rtl0inv <- matrix(0, p, p)
  xtx0 <- matrix(0, p + q, p + q)
  xtx0inv <- matrix(0, p, p)
  yp0 <- numeric(n)
  w0 <- numeric(n)
  w0inv <- numeric(n)
  xlw0 <- matrix(0, n, p)
  k0 <- numeric(p)
  abx <- numeric(n)
  expabx <- numeric(n)
  beta <- numeric(p + q)
  bt <- numeric(p)
  bb <- numeric(q)
  betat <- numeric(p)
  betab <- numeric(q)
  h <- matrix(0, p, q)
  k <- numeric(p)
  qr <- matrix(0, n, q)
  rtr <- matrix(0, p, q)
  rbr <- matrix(0, q, q)
  scoret <- numeric(p)
  scoreb <- numeric(q)
  t <- matrix(0, n, q)
  xrw <- matrix(0, n, q)
  yp <- numeric(n)
  zt <- numeric(p)
  zb <- numeric(q)
  initlspoisregfit(y = y,
                   xl = xl,
                   beta0 = beta0,
                   ql0 = ql0,
                   rtl0 = rtl0,
                   rtl0inv = rtl0inv,
                   xtx0 = xtx0,
                   xtx0inv = xtx0inv,
                   yp0 = yp0,
                   w0 = w0,
                   w0inv = w0inv,
                   xlw0 = xlw0,
                   abx = abx,
                   expabx = expabx)
  retval <- list(family = "poisson",
                 y = y,
                 xl = xl,
                 beta0 = beta0,
                 ql0 = ql0,
                 rtl0 = rtl0,
                 yp0 = yp0,
                 k0 = k0,
                 w0 = w0,
                 w0inv = w0inv,
                 beta = beta,
                 bt = bt,
                 bb = bb,
                 betat = betat,
                 betab = betab,
                 abx = abx,
                 expabx = expabx,
                 h = h,
                 k = k,
                 qr = qr,
                 rtr = rtr,
                 rbr = rbr,
                 scoret = scoret,
                 scoreb = scoreb,
                 t = t,
                 xrw = xrw,
                 yp = yp,
                 zt = zt,
                 zb = zb)
  class(retval) <- "fitdata"
  return(retval)
}

##############################################################################
####################     Likelihood Ratio   ##################################
##############################################################################

allocate.lspoisregmem.lrt <- function(data, beta0, n, p, q) {
  fitdata <- allocate.lspoisregmem.fit(data = data,
                                      beta0 = beta0,
                                      n = n,
                                      p = p,
                                      q = q)
  lnfacty <- sum(lgamma(fitdata$y + 1))
  loglike <- numeric(2)
  testvalue <- numeric(1)
  loglike[1] <- lspoisreglikelihood(y = fitdata$y,
                                    abx = fitdata$abx,
                                    expabx = fitdata$expabx,
                                    lnfacty = lnfacty)
  retval <- list(fitdata = fitdata,
                 family = "poisson",
                 testtype = "lrt",
                 lnfacty = lnfacty,
                 loglike = loglike,
                 testvalue = testvalue)
  class(retval) <- "lsregmem"
  return (retval)
}

##############################################################################
####################     Score              ##################################
##############################################################################

allocate.lspoisregmem.score <- function(data, beta0, n, p, q) {
  y <- data[,1]
  xl <- as.matrix(data[,1:p])
  xl[,1] <- 1.
  ql0 <- matrix(0, n, p)
  rtl0 <- matrix(0, p, p)
  rtl0inv <- matrix(0, p, p)
  xtx0 <- matrix(0, p + q, p + q)
  xtx0inv <- matrix(0, p, p)
  yp0 <- numeric(n)
  w0 <- numeric(n)
  abx <- numeric(n)
  expabx <- numeric(n)
  xlw0 <- matrix(0, n, p)
  xrw0 <- matrix(0, n, q)
  info0 <- matrix(0, q, q)
  score0 <- numeric(q)
  testvalue <- matrix(0, 1, 1)
  initlspoisregscore(y = y,
                     xl = xl,
                     beta = beta0,
                     ql0 = ql0,
                     rtl0 = rtl0,
                     rtl0inv = rtl0inv,
                     xtx0 = xtx0,
                     xtx0inv = xtx0inv,
                     yp0 = yp0,
                     w0 = w0,
                     abx = abx,
                     expabx = expabx,
                     xlw0 = xlw0)
  retval <- list(fitdata = NA,
                 family = "poisson",
                 testtype = "score",
                 w0 = w0,
                 xlw0 = xlw0,
                 xtx0inv = xtx0inv,
                 resids0 = yp0,
                 xtx0 = xtx0,
                 xrw0 = xrw0,
                 info0 = info0,
                 score0 = score0,
                 testvalue = testvalue)
  class(retval) <- c("lsregmem")
  return (retval)
}

##############################################################################
####################     Robust Score       ##################################
##############################################################################

allocate.lspoisregmem.robustscore <- function(data, beta0, n, p, q)  {
  y <- data[,1]
  xl <- as.matrix(data[,1:p])
  xl[,1] <- 1.
  ql0 <- matrix(0, n, p)
  rtl0 <- matrix(0, p, p)
  rtl0inv <- matrix(0, p, p)
  xtx0 <- matrix(0, p + q, p + q)
  xtx0inv <- matrix(0, p, p)
  yp0 <- numeric(n)
  w0 <- numeric(n)
  abx <- numeric(n)
  expabx <- numeric(n)
  expabxp1 <- numeric(n)
  expitabx <- numeric(n)
  xlw0 <- matrix(0, n, p)
  xrw0 <- matrix(0, n, q)
  xlr0 <- matrix(0, n, p)
  xrw0 <- matrix(0, n, q)
  xrr0 <- matrix(0, n, q)
  uut <- matrix(0, p + q, p +q)
  cmat <- matrix(0, q, p + q)
  robustinfo0 <- matrix(0, q, q)
  robustscore0 <- numeric(q)
  testvalue <- matrix(0, 1, 1)
  initlspoisregrobustscore(y = y,
                           xl = xl,
                           beta = beta0,
                           ql0 = ql0,
                           rtl0 = rtl0,
                           rtl0inv = rtl0inv,
                           xtx0 = xtx0,
                           xtx0inv = xtx0inv,
                           yp0 = yp0,
                           w0 = w0,
                           abx = abx,
                           expabx = expabx,
                           xlw0 = xlw0,
                           xlr0 = xlr0,
                           uut = uut)
  retval <- list(fitdata = NA,
                 family = "poisson",
                 testtype = "robustscore",
                 resids0 = yp0,
                 w0 = w0,
                 xlw0 = xlw0,
                 xlr0 = xlr0,
                 xtx0 = xtx0,
                 xtx0inv = xtx0inv,
                 xrw0 = xrw0,
                 xrr0 = xrr0,
                 uut = uut,
                 cmat = cmat,
                 robustinfo0 = robustinfo0,
                 robustscore0 = robustscore0,
                 testvalue = testvalue)
  class(retval) <- "lsregmem"
  return (retval)
}  

##############################################################################
####################     Wald               ##################################
##############################################################################

allocate.lspoisregmem.wald <- function(data, beta0, n, p, q) {
  fitdata <- allocate.lspoisregmem.fit(data = data,
                                      beta0 = beta0,
                                      n = n,
                                      p = p,
                                      q = q)
  xw <- matrix(0, n, p + q)
  w <- numeric(n)
  xtx <- matrix(0, p + q, p + q)
  xtxinv <- matrix(0, p + q, p + q)
  qrq <- matrix(0, n, p + q)
  qrr <- matrix(0, p + q, p + q)
  rinv <- matrix(0, p + q, p + q)
  testvalue <- matrix(0, 1, 1)
  retval <- list(family = "poisson",
                 testtype = "wald",
                 fitdata = fitdata,
                 xw = xw,
                 w = w,
                 q = qrq,
                 r = qrr,
                 rinv = rinv,
                 xtx = xtx,
                 xtxinv = xtxinv,
                 testvalue = testvalue)
  class(retval) <- "lsregmem"
  return(retval)
}

##############################################################################
####################     Robust Wald        ##################################
##############################################################################

allocate.lspoisregmem.robustwald <- function(data, beta0, n, p, q) {
  fitdata <- allocate.lspoisregmem.fit(data = data,
                                      beta0 = beta0,
                                      n = n,
                                      p = p,
                                      q = q)
  xw <- matrix(0, n, p + q)
  w <- numeric(n)
  xtx <- matrix(0, p + q, p + q)
  xtxinv <- matrix(0, p + q, p + q)
  xrr <- matrix(0, n, p + q)
  hws2 <- matrix(0, p + q, p + q)
  qrq <- matrix(0, n, p + q)
  qrr <- matrix(0, p + q, p + q)
  rinv <- matrix(0, p + q, p + q)
  testvalue <- matrix(0, 1, 1)
  retval <- list(family = "poisson",
                 testtype = "robustwald",
                 fitdata = fitdata,
                 xw = xw,
                 w = w,
                 q = qrq,
                 r = qrr,
                 rinv = rinv,
                 xtx = xtx,
                 xtxinv = xtxinv,
                 xrr = xrr,
                 hws2 = hws2,
                 testvalue = testvalue)
  class(retval) <- "lsregmem"
  return(retval)
}