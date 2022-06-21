// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
int initlspoisregfit(const arma::vec &y,
                     const arma::mat &xl,
                     const arma::vec &beta0,
                     arma::mat &ql0,
                     arma::mat &rtl0,
                     arma::mat &rtl0inv,
                     arma::mat &xtx0,
                     arma::mat &xtx0inv,
                     arma::vec &yp0,
                     arma::vec &w0,
                     arma::vec &w0inv,
                     arma::mat &xlw0,
                     arma::vec &abx,
                     arma::vec &expabx) {
  abx = xl * beta0;
  expabx = exp(abx);
  w0 = sqrt(expabx);
  w0inv = 1. / w0;
  xlw0 = arma::diagmat(w0) * xl;
  xtx0.submat(0,0,xl.n_cols -1, xl.n_cols - 1) = xlw0.t() * xlw0;
  qr_econ(ql0, rtl0, xlw0);
  rtl0inv = inv(trimatu(rtl0));
  xtx0inv = rtl0inv * rtl0inv.t();
  yp0 = y - expabx;

  return 0;
}

// [[Rcpp::export]]
int lspoisregfit(const arma::vec &y,
                 const arma::mat &xl,
                 const arma::mat &xr,
                 const arma::vec &beta0,
                 const arma::vec &yp0,
                 const arma::mat &ql,
                 const arma::mat &rtl,
                 const arma::vec &k0,
                 const arma::vec &w0,
                 const arma::vec &winv0,
                 arma::vec &beta,
                 arma::vec &bt,
                 arma::vec &bb,
                 arma::vec &betat,
                 arma::vec &betab,
                 arma::vec &abx,
                 arma::vec &expabx,
                 arma::mat &h,
                 arma::vec &k,
                 arma::mat &qr,
                 arma::mat &rtr,
                 arma::mat &rbr,
                 arma::vec &scoret,
                 arma::vec &scoreb,
                 arma::mat &t,
                 arma::mat &xrw,
                 arma::vec &yp,
                 arma::vec &zt,
                 arma::vec &zb) {
  int p;
  int q;
  int ncov;
  
  p = xl.n_cols;
  q = xr.n_cols;
  ncov = xl.n_cols + q;
  
  beta.subvec(0, ncov - 1) = arma::zeros(ncov);
  beta.subvec(0, p - 1) = beta0;
  
  xrw.cols(0, xrw.n_cols - 1) = xr.each_col() % w0;
  rtr = ql.t() * xrw;
  t = xrw - ql * rtr;
  qr_econ(qr, rbr, t);
  if (solve(h, rtl, rtr, arma::solve_opts::no_approx) == false) {
    beta[0] = 3.;
    return 3;
  }
  
  betat = beta0;
  betab = arma::zeros(q);
  scoreb = xr.t() * yp0;
  if (max(abs(scoreb)) < 1e-6) {
    abx = xl * betat + xr * betab;
    expabx = exp(abx);
    return 0;
  }
  
  yp = yp0 % winv0;
  k = k0;
  for (int i = 0; i < 10; ++i) {
    zb = qr.t() * yp;
    if (solve(bb, rbr, zb, arma::solve_opts::no_approx) == false) {
      beta[0] = 4;
      return 4;
    }
    
    bt = k - h * bb;
    betat += bt;
    betab += bb;
    abx = xl * betat + xr * betab;
    expabx = exp(abx);
    yp = y - expabx;
    scoret = xl.t() * yp;
    scoreb = xr.t() * yp;
    if (max(abs(scoret)) < 1e-6 && max(abs(scoreb)) < 1e-6) {
      beta.subvec(0, p - 1) = betat;
      beta.subvec(p, ncov - 1) = betab;
      return 0;
    }
    
    yp %= winv0;
    zt = ql.t() * yp;
    if (solve(k, rtl, zt, arma::solve_opts::no_approx) == false) {
      beta[0] = 5;
      return 5;
    }
  }
  
  beta[0] = 6.;
  
  return 6;
}

/////////////////////////////////////////////////////////////////////////////
//                   Likelihood Ratio                                      //
/////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
double lspoisreglikelihood(const arma::vec &y,
                           const arma::vec &abx,
                           const arma::vec &expabx,
                           const arma::vec &lnfacty) {
  double loglike;
  loglike = sum(abx % y) - sum(expabx) - lnfacty[0];
  return loglike;
}

// [[Rcpp::export]]
int lspoisreglrt(const arma::vec &y,
                 const arma::vec &abx,
                 const arma::vec &expabx,
                 const arma::vec &lnfacty,
                 arma::vec &loglike,
                 arma::vec &testvalue) {
  loglike[1] = lspoisreglikelihood(y, abx, expabx, lnfacty);
  testvalue[0] = 2*(loglike[1] - loglike[0]);
  return 0;
}

/////////////////////////////////////////////////////////////////////////////
//                   Score                                                 //
/////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
int initlspoisregscore(const arma::vec &y,
                       const arma::mat &xl,
                       const arma::vec &beta,
                       arma::mat &ql0,
                       arma::mat &rtl0,
                       arma::mat &rtl0inv,
                       arma::mat &xtx0,
                       arma::mat &xtx0inv,
                       arma::vec &yp0,
                       arma::vec &w0,
                       arma::vec &abx,
                       arma::vec &expabx,
                       arma::mat &xlw0) {
  abx = xl * beta;
  expabx = exp(abx);
  w0 = sqrt(expabx);
  xlw0 = arma::diagmat(w0) * xl;
  xtx0.submat(0,0,xl.n_cols -1, xl.n_cols - 1) = xlw0.t() * xlw0;
  qr_econ(ql0, rtl0, xlw0);
  rtl0inv = inv(trimatu(rtl0));
  xtx0inv = rtl0inv * rtl0inv.t();
  yp0 = y - expabx;
  
  return 0;
}

/////////////////////////////////////////////////////////////////////////////
//                   Robust Score                                          //
/////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
int initlspoisregrobustscore(const arma::vec &y,
                             const arma::mat &xl,
                             const arma::vec &beta,
                             arma::mat &ql0,
                             arma::mat &rtl0,
                             arma::mat &rtl0inv,
                             arma::mat &xtx0,
                             arma::mat &xtx0inv,
                             arma::vec &yp0,
                             arma::vec &w0,
                             arma::vec &abx,
                             arma::vec &expabx,
                             arma::mat &xlw0,
                             arma::mat &xlr0,
                             arma::mat &uut) {
  initlspoisregscore(y, xl, beta,
                     ql0, rtl0, rtl0inv,
                     xtx0, xtx0inv, yp0, w0,
                     abx, expabx, xlw0);
  xlr0 = arma::diagmat(yp0) * xl;
  uut.submat(0,0,xl.n_cols -1, xl.n_cols - 1) = xlr0.t() * xlr0;
  return 0;
}

/////////////////////////////////////////////////////////////////////////////
//                   Wald                                                  //
/////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
int pois_wald_test(const arma::mat &xl,
                   const arma::mat &xr,
                   const arma::vec &betab,
                   const arma::vec &expabx,
                   arma::vec &w,
                   arma::mat &xw,
                   arma::mat &q,
                   arma::mat &r,
                   arma::mat &rinv,
                   arma::mat &xtx,
                   arma::mat &xtxinv,
                   arma::mat &testvalue) {
  int p, n;
  
  p = xl.n_cols;
  n = xtx.n_cols;
  w = sqrt(expabx);
  xw.cols(0,p-1) = xl.each_col() % w;
  xw.cols(p,n-1) = xr.each_col() % w;
  xtx = xw.t() * xw;
  qr_econ(q, r, xw);
  rinv = inv(trimatu(r));
  xtxinv = rinv * rinv.t();
  if (xr.n_cols == 1)
    testvalue(0,0) = betab(0) / sqrt(xtxinv(n-1,n-1));
  else
    testvalue = betab.t() * inv(xtxinv.submat(p,p,n-1,n-1)) * betab;
  return 0;
}

/////////////////////////////////////////////////////////////////////////////
//                   Robust Wald                                           //
/////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
int pois_robustwald_test(const arma::mat &xl,
                         const arma::mat &xr,
                         const arma::vec &betab,
                         const arma::vec &expabx,
                         const arma::vec &resids,
                         arma::vec &w,
                         arma::mat &xw,
                         arma::mat &xrr,
                         arma::mat &q,
                         arma::mat &r,
                         arma::mat &rinv,
                         arma::mat &xtx,
                         arma::mat &xtxinv,
                         arma::mat &hws2,
                         arma::mat &testvalue) {
  int p, n;
  
  p = xl.n_cols;
  n = xtx.n_cols;
  
  w = sqrt(expabx);
  xw.cols(0,p-1) = xl.each_col() % w;
  xw.cols(p,n-1) = xr.each_col() % w;
  xtx = xw.t() * xw;
  qr_econ(q, r, xw);
  rinv = inv(trimatu(r));
  xtxinv = rinv * rinv.t();
  
  xrr.cols(0,p-1) = xl.each_col() % resids;
  xrr.cols(p,n-1) = xr.each_col() % resids;
  hws2 = xtxinv * xrr.t() * xrr * xtxinv;
  
  if (xr.n_cols == 1)
    testvalue(0,0) = betab(0) / sqrt(hws2(n-1,n-1));
  else
    testvalue = betab.t() * inv(hws2.submat(p,p,n-1,n-1)) * betab;
  return 0;
}