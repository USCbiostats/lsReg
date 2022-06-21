// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

const double LOG2PIP1 = log(M_PI + M_PI) + 1.;

// [[Rcpp::export]]
int initlslinreg2(const arma::vec &y,
                  const arma::mat &xl,
                  arma::mat &ql0,
                  arma::mat &rtl0,
                  arma::mat &rtlinv0,
                  arma::vec &resids0,
                  arma::vec &w0,
                  arma::mat &wxl0,
                  arma::mat &xtx0,
                  arma::mat &xtxinv0,
                  arma::vec &k,
                  arma::vec &zt,
                  arma::vec &sigma2,
                  arma::vec &s2,
                  arma::vec &loglike) {
  qr_econ(ql0, rtl0, xl);
  zt = ql0.t() * y;
  solve(k, rtl0, zt);
  resids0 = y - xl * k;
  sigma2[0] = arma::dot(resids0, resids0) / xl.n_rows;
  s2[0] = sigma2[0] * xl.n_rows / (xl.n_rows - xl.n_cols);
  w0.fill(sqrt(s2[0]));
  wxl0 = sqrt(s2[0]) * xl;
  xtx0.submat(0, 0, xl.n_cols - 1, xl.n_cols - 1) = trans(wxl0) * wxl0;
  qr_econ(ql0, rtl0, wxl0);
  rtlinv0 = arma::inv(trimatu(rtl0));
  xtxinv0 = rtlinv0 * rtlinv0.t();
  loglike[0] = -0.5 * xl.n_rows * (log(sigma2[0]) + LOG2PIP1);
  
  return 0;
}

// [[Rcpp::export]]
int lslinregscore2(const arma::mat &xr,
                   const arma::vec &w0,
                   const arma::mat &wxl0,
                   const arma::mat &xtxinv0,
                   const arma::vec &resids0,
                   arma::mat &wxr0,
                   arma::mat &xltxr0,
                   arma::mat &info0,
                   arma::vec &score0,
                   arma::mat &scoretest) {
  wxr0 = arma::diagmat(w0) * xr;
  xltxr0 = wxl0.t() * wxr0;
  info0 = wxr0.t()*wxr0 - xltxr0.t() * xtxinv0 * xltxr0;
  score0 = resids0.t() * xr;
  if (xr.n_cols == 1)
    scoretest(0,0) = score0(0) / sqrt(info0(0,0));
  else
    scoretest = score0.t() * inv(info0) * score0;
  return 0;
}

