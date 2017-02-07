#ifndef OBJCORE_H_
#define OBJCORE_H_

#include <iostream>
#include <string>
#include <RcppArmadillo.h>

using namespace std;

// [[Rcpp::depends(RcppArmadillo)]]

Rcpp::List snw_cnc_core_obj_core(arma::colvec beta, arma::colvec q, arma::mat ax, arma::mat axExpand, arma::colvec y, arma::vec gCounts, double K);
Rcpp::List bn_logistic_glrt_obj_core(arma::colvec beta, double f0, arma::colvec y, arma::colvec g, arma::colvec a, arma::mat ax, arma::mat ax_expand, arma::colvec g_expand, arma::colvec a_expand, double K, double Fst, unsigned N, unsigned n0, unsigned n1, unsigned c0);
Rcpp::List double_logistic_glrt_obj_core(arma::colvec beta, arma::colvec theta, arma::colvec y, arma::colvec g, arma::mat ax, arma::mat ax_g, arma::mat ax_expand, arma::mat ax_g_expand, arma::colvec g_expand, double w0, double w1);

#endif /* OBJCORE_H_ */