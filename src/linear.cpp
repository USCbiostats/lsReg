// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

const double LOG2PIP1 = log(M_PI + M_PI) + 1.;

/////////////////////////////////////////////////////////////////////////////
//                   Fit                                                   //
/////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
int initlslinregfit(const arma::mat &xl,
                    arma::mat &ql,
                    arma::mat &rtl) {
  qr_econ(ql, rtl, xl);

  return 0;
}

// [[Rcpp::export]]
int lslinregfit(const arma::vec &y,
                const arma::mat &xl,
                const arma::mat &xr,
                const arma::vec &beta0,
                const arma::mat &ql,
                const arma::mat &rtl,
                arma::vec &bt,
                arma::vec &bb,
                arma::mat &qr,
                arma::mat &rtr,
                arma::mat &rbr,
                arma::mat &h,
                arma::mat &t,
                arma::vec &zb) {
  rtr = trans(ql) * xr;
  t = xr - ql * rtr;
  qr_econ(qr, rbr, t);
  solve(h, rtl, rtr);
  zb = qr.t() * y;
  arma::mat bbt = bb.t();
  solve(bbt, rbr, zb);
  bt = beta0 - h*bbt;
  bb = bbt.t();
  
  return 0;
}

/////////////////////////////////////////////////////////////////////////////
//                   Likelihood Ratio                                      //
/////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
int initlslinreglrt(const arma::vec &y,
                    const arma::mat &xl,
                    const arma::vec &beta,
                    arma::vec &resids,
                    arma::vec &loglike) {
  double sigma2;
  
  resids = y - xl*beta;
  sigma2 = arma::dot(resids, resids)/ xl.n_rows;
  loglike[0] = -0.5 * xl.n_rows * (log(sigma2) + LOG2PIP1);
  
  return 0;
}

// [[Rcpp::export]]
int lslinreglrt(const arma::vec &y,
                const arma::mat &xl,
                const arma::mat &xr,
                const arma::vec &bt,
                const arma::vec &bb,
                arma::vec &resids,
                arma::vec &loglike,
                arma::vec &testvalue) {
  double sigma2;
  
  resids = y - xl*bt - xr*bb;
  sigma2 = arma::dot(resids, resids)/ xl.n_rows;
  loglike[1] = -0.5 * xl.n_rows * (log(sigma2) + LOG2PIP1);
  testvalue[0] = 2*(loglike[1] - loglike[0]);
  
  return 0;
}

/////////////////////////////////////////////////////////////////////////////
//                   Score                                                 //
/////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
int initlslinregscore(const arma::vec &y,
                      const arma::mat &xl,
                      const arma::vec &beta,
                      arma::vec &resids0,
                      arma::vec &s2,
                      arma::mat &xtx,
                      arma::mat &xtx0inv) {
  int p;
  arma::mat ql, rtl, rtlinv;
  
  p = xl.n_cols;
  resids0 = y - xl*beta;
  s2[0] = arma::dot(resids0, resids0)/ (xl.n_rows - xl.n_cols);
  xtx.submat(0,0,p-1,p-1) = xl.t() * xl;
  qr_econ(ql, rtl, xl);
  rtlinv = inv(trimatu(rtl));
  xtx0inv = rtlinv * rtlinv.t();

  return 0;
}

// [[Rcpp::export]]
int lslinregscore(const arma::mat &xl,
                  const arma::vec &xr,
                  const arma::mat &xtx0inv,
                  const arma::vec &resids0,
                  const arma::vec &s2,
                  arma::mat &xtx,
                  arma::mat &info,
                  arma::vec &score,
                  arma::mat &testvalue) {
  int p, n;

  p = xl.n_cols;
  n = xl.n_cols + xr.n_cols;

  xtx.submat(p,0,n-1,p-1) = xr.t() * xl;
  xtx.submat(0,p,p-1,n-1) = xtx.submat(p,0,n-1,p-1).t();
  xtx.submat(p,p,n-1,n-1) = xr.t() * xr;
  
  info = (xtx.submat(p,p,n-1,n-1) - xtx.submat(p,0,n-1,p-1)*xtx0inv*xtx.submat(0,p,p-1,n-1))*s2[0];
  score = arma::sum(arma::diagmat(resids0)*xr, 0);
  if (xr.n_cols == 1)
    testvalue(0,0) = score[0] / sqrt(info[0]);
  else
    testvalue = score.t() * inv(info) * score;
  
  return 0;
}

/////////////////////////////////////////////////////////////////////////////
//                   Robust Score                                          //
/////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
int initlslinregrobustscore(const arma::vec &y,
                            const arma::mat &xl,
                            const arma::vec &beta,
                            arma::vec &resids0,
                            arma::mat &xlr0,
                            arma::vec &s2,
                            arma::mat &xtx,
                            arma::mat &xtx0inv,
                            arma::mat &uut,
                            arma::mat &cmat) {
  int p;

  p = xl.n_cols;
  initlslinregscore(y, xl, beta, resids0, s2, xtx, xtx0inv);
  xlr0 = arma::diagmat(resids0) * xl;
  uut.submat(0,0,p-1,p-1) = xlr0.t() * xlr0;
  cmat.submat(0,p,cmat.n_rows-1,cmat.n_cols-1).eye();  
  return 0;
}

// [[Rcpp::export]]
int lslinregrobustscore(const arma::mat &xl,
                        const arma::mat &xlr0,
                        const arma::mat &xr,
                        const arma::mat &xtx0inv,
                        const arma::vec &resids0,
                        arma::mat &xrr0,
                        arma::mat &uut,
                        arma::mat &cmat,
                        arma::vec &score,
                        arma::mat &info,
                        arma::vec &testvalue) {
  int n, p;
  
  p = xlr0.n_cols;
  n = xlr0.n_cols + xrr0.n_cols;

  xrr0 = diagmat(resids0) * xr;  
  uut.submat(0,p,p-1,n-1) = xlr0.t() * xrr0;
  uut.submat(p,0,n-1,p-1) = uut.submat(0,p,p-1,n-1).t();
  uut.submat(p,p,n-1,n-1) = xrr0.t() * xrr0;
  cmat.submat(0,0,n-p-1,p-1) = -xr.t() * xl * xtx0inv;
  score = sum(xr.each_col() % resids0);
  info = cmat * uut * cmat.t();
  
  if (xr.n_cols == 1)
    testvalue(0,0) = score(0) / sqrt(info(0,0));
  else
    testvalue = score.t() * pinv(info) * score;
  
  return 0;  
}

/////////////////////////////////////////////////////////////////////////////
//                   Wald                                                  //
/////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
int initlslinregwald(const arma::mat &rtl,
                     arma::mat &rtlinv) {
  rtlinv = inv(trimatu(rtl));
  return 0;
}

// [[Rcpp::export]]
int lslinregwaldtest(const arma::vec &y,
                     const arma::mat &xl,
                     const arma::mat &xr,
                     const arma::vec &bt,
                     const arma::vec &bb,
                     const arma::mat &rtlinv,
                     const arma::mat &rtr,
                     const arma::mat &rbr,
                     arma::vec &resids,
                     arma::vec &s2,
                     arma::mat &rtrinv,
                     arma::mat &rbrinv,
                     arma::mat &xtxinv,
                     arma::mat &testvalue) {
  int p, n;
  
  p = xl.n_cols;
  n = p + xr.n_cols;
  resids = y - xl*bt - xr*bb;
  s2 = arma::dot(resids,resids)/(xl.n_rows - xl.n_cols - xr.n_cols);
  rbrinv = inv(trimatu(rbr));
  rtrinv = -rtlinv * rtr * rbrinv;
  xtxinv.submat(0,0,p-1,p-1) = rtlinv * rtlinv.t() + rtrinv * rtrinv.t();
  xtxinv.submat(0,p,p-1,n-1) = rtrinv * rbrinv.t();
  xtxinv.submat(p,0,n-1,p-1) = xtxinv.submat(0,p,p-1,n-1).t();
  xtxinv.submat(p,p,n-1,n-1) = rbrinv * rbrinv.t();
  if (xr.n_cols == 1)
    testvalue(0,0) = bb(0) / sqrt(xtxinv(n-1,n-1) * s2[0]);
  else
    testvalue = (bb.t() * pinv(xtxinv.submat(p,p,n-1,n-1)) * bb) / s2[1];

  return 0;
}
                     
/////////////////////////////////////////////////////////////////////////////
//                   Robust Wald                                           //
/////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
int lslinregrobustwaldtest(const arma::vec &y,
                           const arma::mat &xl,
                           const arma::mat &xr,
                           const arma::vec &bt,
                           const arma::vec &bb,
                           const arma::mat &rtlinv,
                           const arma::mat &rtr,
                           const arma::mat &rbr,
                           arma::mat &xlr,
                           arma::mat &xrr,
                           arma::vec &resids,
                           arma::vec &s2,
                           arma::mat &rtrinv,
                           arma::mat &rbrinv,
                           arma::mat &xtxinv,
                           arma::mat &hwinfo,
                           arma::mat &testvalue) {
  int p, n;
  
  p = xl.n_cols;
  n = p + xr.n_cols;
  resids = y - xl*bt - xr*bb;
  s2 = arma::dot(resids,resids)/(xl.n_rows - xl.n_cols - xr.n_cols);
  rbrinv = inv(trimatu(rbr));
  rtrinv = -rtlinv * rtr * rbrinv;
  xtxinv.submat(0,0,p-1,p-1) = rtlinv * rtlinv.t() + rtrinv * rtrinv.t();
  xtxinv.submat(0,p,p-1,n-1) = rtrinv * rbrinv.t();
  xtxinv.submat(p,0,n-1,p-1) = xtxinv.submat(0,p,p-1,n-1).t();
  xtxinv.submat(p,p,n-1,n-1) = rbrinv * rbrinv.t();
  xlr = diagmat(resids) * xl;
  xrr = diagmat(resids) * xr;
  hwinfo.submat(0, 0, p-1, p-1) = xlr.t() * xlr;
  hwinfo.submat(0, p, p-1, n-1) = xlr.t() * xrr;
  hwinfo.submat(p, 0, n-1, p-1) = hwinfo.submat(0, p, p-1, n-1).t();
  hwinfo.submat(p, p, n-1, n-1) = xrr.t() * xrr;
  hwinfo.submat(0, 0, n-1, n-1) = xtxinv * hwinfo * xtxinv;
  if (xr.n_cols == 1)
    testvalue(0,0) = bb(0) / sqrt(hwinfo(n-1,n-1));
  else
    testvalue = bb.t() * pinv(hwinfo.submat(p,p,n-1,n-1)) * bb;
  
  return 0;
}
