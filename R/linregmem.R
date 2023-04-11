allocate.lslinregmem <- function(data, beta0, n, p, q, testtype) {
  tt <- match(testtype, c("lrt", "score", "robustscore", "wald", "robustwald"))
  if (is.na(tt) == TRUE)
    return (NA)
  if (tt == 1) {
    return(allocate.lslinregmem.lrt(data = data,
                                    beta0 = beta0,
                                    n = n,
                                    p = p,
                                    q = q))
    
  }
  if (tt == 2) {
    return (allocate.lslinregmem.score(data = data,
                                       beta0 = beta0,
                                       n = n,
                                       p = p,
                                       q = q))
  }
  if (tt == 3) {
    return (allocate.lslinregmem.robustscore(data = data,
                                             beta0 = beta0,
                                             n = n,
                                             p = p,
                                             q = q))
  }
  if (tt == 4) {
    return (allocate.lslinregmem.wald(data = data,
                                      beta0 = beta0,
                                      n = n,
                                      p = p,
                                      q = q))
  }
  if (tt == 5) {
    return (allocate.lslinregmem.robustwald(data = data,
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

allocate.lslinregmem.fit <- function(data, beta0, n, p, q) {
  y <- data[,1]
  xl <- as.matrix(data[,1:p])
  xl[,1] <- 1.  
  ql <- matrix(0, n, p)
  rtl <- matrix(0, p, p)
  bt <- numeric(p)
  bb <- numeric(q)
  betat <- numeric(p)
  betab <- numeric(q)
  qr <- matrix(0, n, q)
  rtr <- matrix(0, p, q)
  rbr <- matrix(0, q, q)
  xtx0 <- matrix(0, p + q, p + q)
  xtx0inv <- matrix(0, p, p)
  h <- matrix(0, p, q)
  t <- matrix(0, n, q)
  zb <- numeric(q)
  initlslinregfit(xl = xl,
                  ql = ql,
                  rtl = rtl)
  retval <- list(family = "gaussian",
                 y = y,
                 xl = xl,
                 beta0 = beta0,
                 ql = ql,
                 rtl = rtl,
                 bt = bt,
                 bb = bb,
                 betat = betat,
                 betab = betab,
                 qr = qr,
                 rtr = rtr,
                 rbr = rbr,
                 h = h,
                 t = t,
                 zb = zb)
  class(retval) <- "fitdata"
  return(retval)
}

allocate.lslinregmem2.fit <- function(data, beta0, n, p, q) {
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
  initlslinregfit(xl = xl,
                  ql = ql0,
                  rtl = rtl0)
  retval <- list(family = "gaussian",
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

allocate.lslinregmem.lrt <- function(data, beta0, n, p, q) {
  fitdata <- allocate.lslinregmem.fit(data = data,
                                      beta0 = beta0,
                                      n = n,
                                      p = p,
                                      q = q)
  resids <- numeric(n)
  loglike <- numeric(2)
  testvalue <- numeric(1)
  initlslinreglrt(y = fitdata$y,
                  xl = fitdata$xl,
                  beta = fitdata$beta0,
                  resids = resids,
                  loglike = loglike)
  retval <- list(fitdata = fitdata,
                 family = "gaussian",
                 testtype = "lrt",
                 resids = resids,
                 loglike = loglike,
                 testvalue = testvalue)
  class(retval) <- "lsregmem"
  return (retval)
}

##############################################################################
####################     Score              ##################################
##############################################################################

allocate.lslinregmem.score <- function(data, beta0, n, p, q) {
  y <- data[,1]
  xl <- as.matrix(data[,1:p])
  xl[,1] <- 1.
  xtx <- matrix(0, p + q, p + q)
  xtx0inv <- matrix(0, p, p)
  resids0 <- numeric(n)
  s2 <- numeric(1)
  info <- matrix(0, q, q)
  score <- numeric(q)
  testvalue <- matrix(0, 1, 1)
  initlslinregscore(y = y,
                    xl = xl,
                    beta = beta0,
                    resids0 = resids0,
                    s2 = s2,
                    xtx = xtx,
                    xtx0inv = xtx0inv)
  retval <- list(fitdata = NA,
                 family = "gaussian",
                 testtype = "score",
                 xl = xl,
                 xtx0inv = xtx0inv,
                 resids0 = resids0,
                 s2 = s2,
                 xtx = xtx,
                 info0 = info,
                 score0 = score,
                 testvalue = testvalue)
  class(retval) <- c("lsregmem")
  return (retval)
}

##############################################################################
####################     Robust Score       ##################################
##############################################################################

allocate.lslinregmem.robustscore <- function(data, beta0, n, p, q)  {
  y <- data[,1]
  xl <- as.matrix(data[,1:p])
  xl[,1] <- 1.
  xtx <- matrix(0, p + q, p + q)
  xtx0inv <- matrix(0, p, p)
  resids0 <- numeric(n)
  s2 <- numeric(1)
  xlr0 <- matrix(0, n, p)
  xrr0 <- matrix(0, n, q)
  uut <- matrix(0, p + q, p + q)
  cmat <- matrix(0, q, p + q)
  info <- matrix(0, q, q)
  score <- numeric(q)
  testvalue <- matrix(0, 1, 1)
  initlslinregrobustscore(y = y,
                          xl = xl,
                          beta = beta0,
                          resids0 = resids0,
                          xlr0 = xlr0,
                          s2 = s2,
                          xtx = xtx,
                          xtx0inv = xtx0inv,
                          uut = uut,
                          cmat = cmat)
  retval <- list(fitdata = NA,
                 family = "gaussian",
                 testtype = "robustscore",
                 xl = xl,
                 xlr0 = xlr0,
                 xtx0inv = xtx0inv,
                 resids0 = resids0,
                 xrr0 = xrr0,
                 uut = uut,
                 cmat = cmat,
                 info = info,
                 score = score,
                 testvalue = testvalue)
  class(retval) <- "lsregmem"
  return (retval)
}  

##############################################################################
####################     Wald               ##################################
##############################################################################

allocate.lslinregmem.wald <- function(data, beta0, n, p, q) {
  fitdata <- allocate.lslinregmem.fit(data = data,
                                      beta0 = beta0,
                                      n = n,
                                      p = p,
                                      q = q)
  xtxinv <- matrix(0, p + q, p + q)
  rtlinv <- matrix(0, p, p)
  rtrinv <- matrix(0, p , q)
  rbrinv <- matrix(0, q, q)
  resids <- numeric(n)
  s2 <- numeric(1)
  testvalue <- matrix(0, 1, 1)
  initlslinregwald(rtl = fitdata$rtl,
                   rtlinv = rtlinv)
  retval <- list(family = "gaussian",
                 testtype = "wald",
                 fitdata = fitdata,
                 rtlinv = rtlinv,
                 resids = resids,
                 xtxinv = xtxinv,
                 rtrinv = rtrinv,
                 rbrinv = rbrinv,
                 s2 = s2,
                 testvalue = testvalue)
  class(retval) <- "lsregmem"
  return(retval)
}

##############################################################################
####################     Robust Wald        ##################################
##############################################################################

allocate.lslinregmem.robustwald <- function(data, beta0, n, p, q) {
  fitdata <- allocate.lslinregmem.fit(data = data,
                                      beta0 = beta0,
                                      n = n,
                                      p = p,
                                      q = q)
  xtxinv <- matrix(0, p + q, p + q)
  rtlinv <- matrix(0, p, p)
  rtrinv <- matrix(0, p , q)
  rbrinv <- matrix(0, q, q)
  xlr <- matrix(0, n, p)
  xrr <- matrix(0, n, q)
  hwinfo <- matrix(0, p + q, p + q)
  resids <- numeric(n)
  s2 <- numeric(1)
  testvalue <- matrix(0, 1, 1)
  initlslinregwald(rtl = fitdata$rtl,
                   rtlinv = rtlinv)
  retval <- list(family = "gaussian",
                 testtype = "robustwald",
                 fitdata = fitdata,
                 rtlinv = rtlinv,
                 resids = resids,
                 xtxinv = xtxinv,
                 rtrinv = rtrinv,
                 rbrinv = rbrinv,
                 xlr = xlr,
                 xrr = xrr,
                 s2 = s2,
                 hwinfo = hwinfo,
                 testvalue = testvalue)
  class(retval) <- "lsregmem"
  return(retval)
}
