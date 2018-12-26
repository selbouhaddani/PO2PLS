#include <iostream>
#include <RcppArmadillo.h>
//#include <boost/timer/timer.hpp>
//#include <cmath>
// [[Rcpp::depends(RcppArmadillo)]]
//// [[Rcpp::depends(BH)]]

using namespace std;
using namespace std;
using namespace Rcpp;
using namespace arma;


// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
List E_step_testC(arma::mat dataXY, int p, int q, int r, int rx, int ry, int N,
                arma::mat SigmaZ, arma::mat GammaEF, arma::mat invZtilde,
                arma::mat Gamma, arma::mat GGef, arma::mat EZc, arma::mat Szz,
                arma::mat W, arma::mat C, arma::mat Wo, arma::mat Co, arma::mat SigT,
                arma::mat SigH, double sig2E, double sig2F){

  mat X = dataXY.cols(0, p-1);
  mat Y = dataXY.cols(p, p+q-1);

  mat mu_EF = dataXY;
  mu_EF = mu_EF - dataXY * (GammaEF * invZtilde) * Gamma.t();

  mat tempp = join_vert(join_horiz(W, join_horiz(zeros<mat>(p,r), join_horiz(Wo, zeros<mat>(p,ry)))),
                        join_horiz(zeros<mat>(q,r), join_horiz(C*0, join_horiz(zeros<mat>(q,rx), Co*0))));
  double Cee = trace(tempp.t() * tempp * invZtilde) / p + pow(norm(mu_EF.cols(0, p-1),"fro"),2) / N / p ;

  tempp = join_vert(join_horiz(W*0, join_horiz(zeros<mat>(p,r), join_horiz(Wo*0, zeros<mat>(p,ry)))),
                    join_horiz(zeros<mat>(q,r), join_horiz(C, join_horiz(zeros<mat>(q,rx), Co))));
  double Cff = trace(tempp.t() * tempp * invZtilde) / q + pow(norm(mu_EF.cols(p, p+q-1),"fro"),2) / N / q ;

  mat covH = join_vert(W*0, C * SigH) ;
  mat covHEF = join_vert(W*0, C * SigH / sig2F) ;

  mat mu_H = dataXY  *  covHEF ;
  mu_H = mu_H - (dataXY  *  (GammaEF  *  invZtilde))  *  (trans(Gamma)  *  covHEF) ;

  mat Chh = SigH;
  Chh = Chh - trans(covH)  *  covHEF ;
  Chh = Chh + (trans(covH)  *  GammaEF  *  invZtilde)  *  (trans(Gamma)  *  covHEF) ;
  Chh = Chh + trans(mu_H)*mu_H / N ;

  double logdet = log(det(eye<mat>(2*r+rx+ry,2*r+rx+ry) + GGef * SigmaZ)) + p * log(sig2E) + q * log(sig2F) ;
  double XYinvS = pow(norm(join_horiz(X / sqrt(sig2E), Y / sqrt(sig2F)), "fro"),2) ;
  XYinvS = XYinvS - trace(trans(dataXY  *  GammaEF) * (dataXY  *  GammaEF)  *  invZtilde) ;

  double loglik = N*(p+q)*log(2*datum::pi) + N * logdet + XYinvS ;
  loglik = - loglik/2 ;

  double comp_log = -N / 2 * (p+q) * log(2 * datum::pi) ;
  comp_log = comp_log - N/2*(p*log(sig2E)+q*log(sig2F)) ;
  comp_log = comp_log - N/2*pow(norm(join_horiz(X/sqrt(sig2E), Y/sqrt(sig2F)), "fro"),2) ;
  comp_log = comp_log + N*trace(trans(EZc) * dataXY * GammaEF) ;
  comp_log = comp_log - N/2*trace(GGef * Szz) ;

  List ret;
  ret["EZc"] = EZc;
  ret["Szz"] = Szz;
  ret["mu_T"] = EZc.cols(0,r-1);
  ret["mu_U"] = EZc.cols(r,2*r-1);
  ret["mu_To"] = EZc.cols(2*r,2*r+rx-1);
  ret["mu_Uo"] = EZc.cols(2*r+rx,2*r+rx+ry-1);
  ret["Stt"] = Szz(span(0,r-1),span(0,r-1));
  ret["Suu"] = Szz(span(r,2*r-1),span(r,2*r-1));
  ret["Stoto"] = Szz(span(2*r,2*r+rx-1),span(2*r,2*r+rx-1));
  ret["Suouo"] = Szz(span(2*r+rx,2*r+rx+ry-1),span(2*r+rx,2*r+rx+ry-1));
  ret["Sut"] = Szz(span(r,2*r-1), span(0,r-1));
  ret["Stto"] = Szz(span(0,r-1),span(2*r, 2*r+rx-1));
  ret["Suuo"] = Szz(span(r,2*r-1),span(2*r+rx,2*r+rx+ry-1));
  ret["See"] = Cee;
  ret["Sff"] = Cff;
  ret["Shh"] = Chh;
  ret["loglik"] = loglik;
  ret["comp_log"] = comp_log;

  return ret;
}


