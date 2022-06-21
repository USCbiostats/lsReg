#' Title
#'
#' @param lsregmem Memory allocated for large scale regression
#' @param xr Additional variables to be tested
#'
#' @return Nonzero value indicates an error
#' @export
#'
#' @examples
runtest <- function(lsregmem, xr) {
  if (class(lsregmem) != "lsregmem")
    return(1)
  if (is.list(lsregmem$fitdata) == TRUE) {
#    lsreg.fit(fitdata = lsregmem$fitdata,
#              xr = xr)
    if (lsregmem$fitdata$family == "gaussian") {
      lslinregfit(y = lsregmem$fitdata$y,
                  xl = lsregmem$fitdata$xl,
                  xr = xr,
                  beta0 = lsregmem$fitdata$beta0,
                  ql = lsregmem$fitdata$ql,
                  rtl = lsregmem$fitdata$rtl,
                  bt = lsregmem$fitdata$bt,
                  bb = lsregmem$fitdata$bb,
                  qr = lsregmem$fitdata$qr,
                  rtr = lsregmem$fitdata$rtr,
                  rbr = lsregmem$fitdata$rbr,
                  h = lsregmem$fitdata$h,
                  t = lsregmem$fitdata$t,
                  zb = lsregmem$fitdata$zb)
    } else if (lsregmem$fitdata$family == "binomial") {
      lslogregfit(y = lsregmem$fitdata$y,
                  xl = lsregmem$fitdata$xl,
                  xr = xr,
                  beta0 = lsregmem$fitdata$beta0,
                  yp0 = lsregmem$fitdata$yp0,
                  ql = lsregmem$fitdata$ql0,
                  rtl = lsregmem$fitdata$rtl0,
                  k0 = lsregmem$fitdata$k0,
                  w0 = lsregmem$fitdata$w0,
                  winv0 = lsregmem$fitdata$w0inv,
                  beta = lsregmem$fitdata$beta,
                  bt = lsregmem$fitdata$bt,
                  bb = lsregmem$fitdata$bb,
                  betat = lsregmem$fitdata$betat,
                  betab = lsregmem$fitdata$betab,
                  abx = lsregmem$fitdata$abx,
                  expabx = lsregmem$fitdata$expabx,
                  expabxp1 = lsregmem$fitdata$expabxp1,
                  expitabx = lsregmem$fitdata$expitabx,
                  h = lsregmem$fitdata$h,
                  k = lsregmem$fitdata$k,
                  qr = lsregmem$fitdata$qr,
                  rtr = lsregmem$fitdata$rtr,
                  rbr = lsregmem$fitdata$rbr,
                  scoret = lsregmem$fitdata$scoret,
                  scoreb = lsregmem$fitdata$scoreb,
                  t = lsregmem$fitdata$t,
                  xrw = lsregmem$fitdata$xrw,
                  yp = lsregmem$fitdata$yp,
                  zt = lsregmem$fitdata$zt,
                  zb = lsregmem$fitdata$zb)
    } else if (lsregmem$fitdata$family == "poisson") {
      lspoisregfit(y = lsregmem$fitdata$y,
                   xl = lsregmem$fitdata$xl,
                   xr = xr,
                   beta0 = lsregmem$fitdata$beta0,
                   yp0 = lsregmem$fitdata$yp0,
                   ql = lsregmem$fitdata$ql0,
                   rtl = lsregmem$fitdata$rtl0,
                   k0 = lsregmem$fitdata$k0,
                   w0 = lsregmem$fitdata$w0,
                   winv0 = lsregmem$fitdata$w0inv,
                   beta = lsregmem$fitdata$beta,
                   bt = lsregmem$fitdata$bt,
                   bb = lsregmem$fitdata$bb,
                   betat = lsregmem$fitdata$betat,
                   betab = lsregmem$fitdata$betab,
                   abx = lsregmem$fitdata$abx,
                   expabx = lsregmem$fitdata$expabx,
                   h = lsregmem$fitdata$h,
                   k = lsregmem$fitdata$k,
                   qr = lsregmem$fitdata$qr,
                   rtr = lsregmem$fitdata$rtr,
                   rbr = lsregmem$fitdata$rbr,
                   scoret = lsregmem$fitdata$scoret,
                   scoreb = lsregmem$fitdata$scoreb,
                   t = lsregmem$fitdata$t,
                   xrw = lsregmem$fitdata$xrw,
                   yp = lsregmem$fitdata$yp,
                   zt = lsregmem$fitdata$zt,
                   zb = lsregmem$fitdata$zb)
    }
  }
  if (lsregmem$testtype == "lrt") {
    if (lsregmem$family == "gaussian") {
      return(lslinreglrt(y = lsregmem$fitdata$y,
                         xl = lsregmem$fitdata$xl,
                         xr = xr,
                         bt = lsregmem$fitdata$bt,
                         bb = lsregmem$fitdata$bb,
                         resids = lsregmem$resids,
                         loglike = lsregmem$loglike,
                         testvalue = lsregmem$testvalue))
    } else if (lsregmem$family == "binomial") {
      return(lslogreglrt(y = lsregmem$fitdata$y,
                         abx = lsregmem$fitdata$abx,
                         expabxp1 = lsregmem$fitdata$expabxp1,
                         loglike = lsregmem$loglike,
                         testvalue = lsregmem$testvalue))
    } else if (lsregmem$family == "poisson") {
      return(lspoisreglrt(y = lsregmem$fitdata$y,
                          abx = lsregmem$fitdata$abx,
                          expabx = lsregmem$fitdata$expabx,
                          lnfacty = lsregmem$lnfacty,
                          loglike = lsregmem$loglike,
                          testvalue = lsregmem$testvalue))
    }
    return (1)
  } else if (lsregmem$testtype == "score") {
    if (lsregmem$family == "binomial" || lsregmem$family == "poisson") {
      return (score_test(xr = xr,
                         w0 = lsregmem$w0,
                         xlw0 = lsregmem$xlw0,
                         xtxinv0 = lsregmem$xtx0inv,
                         resids0 = lsregmem$resids0,
                         xtx0 = lsregmem$xtx0,
                         xrw0 = lsregmem$xrw0,
                         info0 = lsregmem$info0,
                         score0 = lsregmem$score0,
                         testvalue = lsregmem$testvalue))
    } else if (lsregmem$family == "gaussian") {
      return (lslinregscore(xl = lsregmem$xl,
                            xr = xr,
                            xtx0inv = lsregmem$xtx0inv,
                            resids0 = lsregmem$resids0,
                            s2 = lsregmem$s2,
                            xtx = lsregmem$xtx,
                            info = lsregmem$info,
                            score = lsregmem$score,
                            testvalue = lsregmem$testvalue))
    }
    return (1)
  } else if (lsregmem$testtype == "robustscore") {
    if (lsregmem$family == "binomial" || lsregmem$family == "poisson") {
      return (robust_score_test(xr = xr,
                                resids0 = lsregmem$resids0,
                                w0 = lsregmem$w0,
                                xlw0 = lsregmem$xlw0,
                                xlr0 = lsregmem$xlr0,
                                xtx0 = lsregmem$xtx0,
                                xtx0inv = lsregmem$xtx0inv,
                                xrw0 = lsregmem$xrw0,
                                xrr0 = lsregmem$xrr0,
                                uut = lsregmem$uut,
                                cmat = lsregmem$cmat,
                                robscore0 = lsregmem$robustscore0,
                                robinfo0 = lsregmem$robustinfo0,
                                testvalue = lsregmem$testvalue))
    } else if (lsregmem$family == "gaussian") {
      return (lslinregrobustscore(xl = lsregmem$xl,
                                  xlr0 = lsregmem$xlr0,
                                  xr = xr,
                                  xtx0inv = lsregmem$xtx0inv,
                                  resids0 = lsregmem$resids0,
                                  xrr0 = lsregmem$xrr0,
                                  uut = lsregmem$uut,
                                  cmat = lsregmem$cmat,
                                  info = lsregmem$info,
                                  score = lsregmem$score,
                                  testvalue = lsregmem$testvalue))
    }
    return (1)
  } else if (lsregmem$testtype == "wald") {
    if (lsregmem$family == "binomial") {
      return(log_wald_test(xl = lsregmem$fitdata$xl,
                           xr = xr,
                           betab = lsregmem$fitdata$betab,
                           expitabx = lsregmem$fitdata$expitabx,
                           w = lsregmem$w,
                           xw = lsregmem$xw,
                           q = lsregmem$q,
                           r = lsregmem$r,
                           rinv = lsregmem$rinv,
                           xtx = lsregmem$xtx,
                           xtxinv = lsregmem$xtxinv,
                           testvalue = lsregmem$testvalue))
    } else if (lsregmem$family == "poisson") {
      return(pois_wald_test(xl = lsregmem$fitdata$xl,
                            xr = xr,
                            betab = lsregmem$fitdata$betab,
                            expabx = lsregmem$fitdata$expabx,
                            w = lsregmem$w,
                            xw = lsregmem$xw,
                            q = lsregmem$q,
                            r = lsregmem$r,
                            rinv = lsregmem$rinv,
                            xtx = lsregmem$xtx,
                            xtxinv = lsregmem$xtxinv,
                            testvalue = lsregmem$testvalue))
    } else if (lsregmem$family == "gaussian") {
      return (lslinregwaldtest(y = lsregmem$fitdata$y,
                               xl = lsregmem$fitdata$xl,
                               xr = xr,
                               bt = lsregmem$fitdata$bt,
                               bb = lsregmem$fitdata$bb,
                               rtlinv = lsregmem$rtlinv,
                               rtr = lsregmem$fitdata$rtr,
                               rbr = lsregmem$fitdata$rbr,
                               resids = lsregmem$resids,
                               s2 = lsregmem$s2,
                               rtrinv = lsregmem$rtrinv,
                               rbrinv = lsregmem$rbrinv,
                               xtxinv = lsregmem$xtxinv,
                               testvalue = lsregmem$testvalue))
    }
    return (1)
  } else if (lsregmem$testtype == "robustwald") {
    if (lsregmem$family == "binomial") {
      return(log_robustwald_test(xl = lsregmem$fitdata$xl,
                                 xr = xr,
                                 betab = lsregmem$fitdata$betab,
                                 expitabx = lsregmem$fitdata$expitabx,
                                 resids = lsregmem$fitdata$yp,
                                 w = lsregmem$w,
                                 xw = lsregmem$xw,
                                 xrr = lsregmem$xrr,
                                 q = lsregmem$q,
                                 r = lsregmem$r,
                                 rinv = lsregmem$rinv,
                                 xtx = lsregmem$xtx,
                                 xtxinv = lsregmem$xtxinv,
                                 hws2 = lsregmem$hws2,
                                 testvalue = lsregmem$testvalue))
    } else if (lsregmem$family == "poisson") {
      return(pois_robustwald_test(xl = lsregmem$fitdata$xl,
                                  xr = xr,
                                  betab = lsregmem$fitdata$betab,
                                  expabx = lsregmem$fitdata$expabx,
                                  resids = lsregmem$fitdata$yp,
                                  w = lsregmem$w,
                                  xw = lsregmem$xw,
                                  xrr = lsregmem$xrr,
                                  q = lsregmem$q,
                                  r = lsregmem$r,
                                  rinv = lsregmem$rinv,
                                  xtx = lsregmem$xtx,
                                  xtxinv = lsregmem$xtxinv,
                                  hws2 = lsregmem$hws2,
                                  testvalue = lsregmem$testvalue))
    } else if (lsregmem$family == "gaussian") {
      return (lslinregrobustwaldtest(y = lsregmem$fitdata$y,
                                     xl = lsregmem$fitdata$xl,
                                     xr = xr,
                                     bt = lsregmem$fitdata$bt,
                                     bb = lsregmem$fitdata$bb,
                                     rtlinv = lsregmem$rtlinv,
                                     rtr = lsregmem$fitdata$rtr,
                                     rbr = lsregmem$fitdata$rbr,
                                     xlr = lsregmem$xlr,
                                     xrr = lsregmem$xrr,
                                     resids = lsregmem$resids,
                                     s2 = lsregmem$s2,
                                     rtrinv = lsregmem$rtrinv,
                                     rbrinv = lsregmem$rbrinv,
                                     xtxinv = lsregmem$xtxinv,
                                     hwinfo = lsregmem$hwinfo,
                                     testvalue = lsregmem$testvalue))
    }
    return (1)
  }
  return(1)
}

likelihoodratio.test <- function(lsregmem, xr) {
  if (lsregmem$family == "gaussian") {
    return(lslinreglrt(y = lsregmem$fitdata$y,
                       xl = lsregmem$fitdata$xl,
                       xr = xr,
                       bt = lsregmem$fitdata$bt,
                       bb = lsregmem$fitdata$bb,
                       resids = lsregmem$resids,
                       loglike = lsregmem$loglike,
                       testvalue = lsregmem$testvalue))
  }
  if (lsregmem$family == "binomial") {
    return(lslogreglrt(y = lsregmem$fitdata$y,
                       abx = lsregmem$fitdata$abx,
                       expabxp1 = lsregmem$fitdata$expabxp1,
                       loglike = lsregmem$loglike,
                       testvalue = lsregmem$testvalue))
  }
  if (lsregmem$family == "poisson") {
    return(lspoisreglrt(y = lsregmem$fitdata$y,
                        abx = lsregmem$fitdata$abx,
                        expabx = lsregmem$fitdata$expabx,
                        lnfacty = lsregmem$lnfacty,
                        loglike = lsregmem$loglike,
                        testvalue = lsregmem$testvalue))
  }
  return(1)
}

score.test <- function(lsregmem, xr) {
  if (lsregmem$family == "binomial" || lsregmem$family == "poisson")
    return (score_test(xr = xr,
                       w0 = lsregmem$w0,
                       xlw0 = lsregmem$xlw0,
                       xtxinv0 = lsregmem$xtx0inv,
                       resids0 = lsregmem$resids0,
                       xtx0 = lsregmem$xtx0,
                       xrw0 = lsregmem$xrw0,
                       info0 = lsregmem$info0,
                       score0 = lsregmem$score0,
                       testvalue = lsregmem$testvalue))
  if (lsregmem$family == "gaussian")
    return (lslinregscore(xl = lsregmem$xl,
                          xr = xr,
                          xtx0inv = lsregmem$xtx0inv,
                          resids0 = lsregmem$resids0,
                          s2 = lsregmem$s2,
                          xtx = lsregmem$xtx,
                          info = lsregmem$info,
                          score = lsregmem$score,
                          testvalue = lsregmem$testvalue))
  return (1)
}

robustscore.test <- function(lsregmem, xr) {
  if (lsregmem$family == "binomial" || lsregmem$family == "poisson")
    return (robust_score_test(xr = xr,
                              resids0 = lsregmem$resids0,
                              w0 = lsregmem$w0,
                              xlw0 = lsregmem$xlw0,
                              xlr0 = lsregmem$xlr0,
                              xtx0 = lsregmem$xtx0,
                              xtx0inv = lsregmem$xtx0inv,
                              xrw0 = lsregmem$xrw0,
                              xrr0 = lsregmem$xrr0,
                              uut = lsregmem$uut,
                              cmat = lsregmem$cmat,
                              robscore0 = lsregmem$robustscore0,
                              robinfo0 = lsregmem$robustinfo0,
                              testvalue = lsregmem$testvalue))
  if (lsregmem$family == "gaussian")
    return (lslinregrobustscore(xl = lsregmem$xl,
                                xlr0 = lsregmem$xlr0,
                                xr = xr,
                                xtx0inv = lsregmem$xtx0inv,
                                resids0 = lsregmem$resids0,
                                xrr0 = lsregmem$xrr0,
                                uut = lsregmem$uut,
                                cmat = lsregmem$cmat,
                                info = lsregmem$info,
                                score = lsregmem$score,
                                testvalue = lsregmem$testvalue))
  return(1)
}

wald.test <- function(lsregmem, xr) {
  if (lsregmem$family == "binomial")
    return(log_wald_test(xl = lsregmem$fitdata$xl,
                         xr = xr,
                         betab = lsregmem$fitdata$betab,
                         expitabx = lsregmem$fitdata$expitabx,
                         w = lsregmem$w,
                         xw = lsregmem$xw,
                         q = lsregmem$q,
                         r = lsregmem$r,
                         rinv = lsregmem$rinv,
                         xtx = lsregmem$xtx,
                         xtxinv = lsregmem$xtxinv,
                         testvalue = lsregmem$testvalue))
  if (lsregmem$family == "poisson")
    return(pois_wald_test(xl = lsregmem$fitdata$xl,
                          xr = xr,
                          betab = lsregmem$fitdata$betab,
                          expabx = lsregmem$fitdata$expabx,
                          w = lsregmem$w,
                          xw = lsregmem$xw,
                          q = lsregmem$q,
                          r = lsregmem$r,
                          rinv = lsregmem$rinv,
                          xtx = lsregmem$xtx,
                          xtxinv = lsregmem$xtxinv,
                          testvalue = lsregmem$testvalue))
  if (lsregmem$family == "gaussian")
    return (lslinregwaldtest(y = lsregmem$fitdata$y,
                             xl = lsregmem$fitdata$xl,
                             xr = xr,
                             bt = lsregmem$fitdata$bt,
                             bb = lsregmem$fitdata$bb,
                             rtlinv = lsregmem$rtlinv,
                             rtr = lsregmem$fitdata$rtr,
                             rbr = lsregmem$fitdata$rbr,
                             resids = lsregmem$resids,
                             s2 = lsregmem$s2,
                             rtrinv = lsregmem$rtrinv,
                             rbrinv = lsregmem$rbrinv,
                             xtxinv = lsregmem$xtxinv,
                             testvalue = lsregmem$testvalue))
  return(1)
}

robustwald.test <- function(lsregmem, xr) {
  if (lsregmem$family == "binomial")
    return(log_robustwald_test(xl = lsregmem$fitdata$xl,
                               xr = xr,
                               betab = lsregmem$fitdata$betab,
                               expitabx = lsregmem$fitdata$expitabx,
                               resids = lsregmem$fitdata$yp,
                               w = lsregmem$w,
                               xw = lsregmem$xw,
                               xrr = lsregmem$xrr,
                               q = lsregmem$q,
                               r = lsregmem$r,
                               rinv = lsregmem$rinv,
                               xtx = lsregmem$xtx,
                               xtxinv = lsregmem$xtxinv,
                               hws2 = lsregmem$hws2,
                               testvalue = lsregmem$testvalue))
  if (lsregmem$family == "poisson")
    return(pois_robustwald_test(xl = lsregmem$fitdata$xl,
                                xr = xr,
                                betab = lsregmem$fitdata$betab,
                                expabx = lsregmem$fitdata$expabx,
                                resids = lsregmem$fitdata$yp,
                                w = lsregmem$w,
                                xw = lsregmem$xw,
                                xrr = lsregmem$xrr,
                                q = lsregmem$q,
                                r = lsregmem$r,
                                rinv = lsregmem$rinv,
                                xtx = lsregmem$xtx,
                                xtxinv = lsregmem$xtxinv,
                                hws2 = lsregmem$hws2,
                                testvalue = lsregmem$testvalue))
  if (lsregmem$family == "gaussian")
    return (lslinregrobustwaldtest(y = lsregmem$fitdata$y,
                                   xl = lsregmem$fitdata$xl,
                                   xr = xr,
                                   bt = lsregmem$fitdata$bt,
                                   bb = lsregmem$fitdata$bb,
                                   rtlinv = lsregmem$rtlinv,
                                   rtr = lsregmem$fitdata$rtr,
                                   rbr = lsregmem$fitdata$rbr,
                                   xlr = lsregmem$xlr,
                                   xrr = lsregmem$xrr,
                                   resids = lsregmem$resids,
                                   s2 = lsregmem$s2,
                                   rtrinv = lsregmem$rtrinv,
                                   rbrinv = lsregmem$rbrinv,
                                   xtxinv = lsregmem$xtxinv,
                                   hwinfo = lsregmem$hwinfo,
                                   testvalue = lsregmem$testvalue))
  
  return(1)
}