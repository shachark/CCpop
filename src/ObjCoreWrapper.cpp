#include <iostream>
#include "ObjCoreWrapper.h"

using namespace std; 

// [[Rcpp::export]]
Rcpp::List snw_cnc_core_obj_core(arma::vec beta, arma::vec q, arma::mat ax, arma::mat axExpand, 
  arma::vec y, arma::vec gCounts, double K)
{
    arma::vec eta = ax * beta;  
    arma::vec exp_eta = exp(eta);
    arma::vec P1 = exp_eta / (1 + exp_eta);
    arma::vec P0 = 1 - P1;
    arma::vec p1 = exp_eta / ((1 + exp_eta) % (1 + exp_eta));
    arma::vec p0 = -p1;
    
    unsigned n = ax.n_rows, nCoeffsAltern = ax.n_cols;
    arma::vec P = arma::zeros<arma::vec>(n)+1;
    arma::vec p = arma::zeros<arma::vec>(n);
    
    arma::uvec y0Idx = find(y == 0);
    arma::uvec y1Idx = find(y == 1);
    arma::uvec ymIdx = find_nonfinite(y);
    
    for (unsigned i = 0; i<y0Idx.n_rows; i++) {
      P(y0Idx(i)) = P0(y0Idx(i));
      p(y0Idx(i)) = p0(y0Idx(i));
    }
    for (unsigned i = 0; i<y1Idx.n_rows; i++) {
      P(y1Idx(i)) = P1(y1Idx(i));
      p(y1Idx(i)) = p1(y1Idx(i));
    }
    
    // The CC-only part    
    double ll_cc = sum(log(P));
    unsigned nr_g_levels = gCounts.n_elem;
    arma::mat temp = trans(p / P) * ax;
    arma::mat temp2 = arma::zeros<arma::rowvec>(nr_g_levels);
    arma::vec dll_cc = trans(join_rows(trans(p / P) * ax, arma::zeros<arma::rowvec>(nr_g_levels)));
    
    // The CCpop G part
    double ll_ccpopG = sum(gCounts % log(q));
    arma::mat dll_ccpopG = join_cols(arma::zeros<arma::vec>(nCoeffsAltern), gCounts / q);
    
    // The CCpop E part
    eta = axExpand * beta;
    exp_eta = exp(eta);
    P1 = exp_eta / (1 + exp_eta);
    P0 = 1 - P1;
    p1 = exp_eta / ((1 + exp_eta) % (1 + exp_eta));
    p0 = -p1;
    unsigned n0 = y0Idx.n_elem;
    unsigned n1 = y1Idx.n_elem;
    double w0 = n0 / (1 - K);
    double w1 = n1 / K;
    arma::vec Pw = w0 * P0 + w1 * P1;
    arma::vec pw = w0 * p0 + w1 * p1;
    
    arma::vec Q = trans(vectorise(repmat(q, 1, n), 1));
    
    arma::vec c0mQPw = ymIdx.n_elem + sum(reshape(Q % Pw, n, nr_g_levels),1);
    arma::mat Qpw = arma::zeros<arma::mat>(axExpand.n_rows,axExpand.n_cols);
    for (unsigned i = 0; i<axExpand.n_cols; i++) {
      Qpw.col(i) = (Q % pw) % axExpand.col(i);
    }
    
    arma::mat qPw;
    if (nr_g_levels == 2) {
      qPw = join_rows(Pw.rows(0, n-1), Pw.rows(n, 2*n-1));
      Qpw = Qpw.rows(0, n-1) + Qpw.rows(n, 2*n-1);
    } else {
      qPw = join_rows(join_rows(Pw.rows(0, n-1), Pw.rows(n, 2*n-1)), Pw.rows(2*n, 3*n-1));
      Qpw = Qpw.rows(0, n-1) + Qpw.rows(n, 2*n-1) + Qpw.rows(2*n, 3*n-1);
    }
    
    double ll_ccpopE = sum(log(c0mQPw));
    arma::mat QpwqPw = join_rows(Qpw, qPw);
    arma::mat QpwqPwc0mQPw = arma::zeros<arma::mat>(n, QpwqPw.n_cols);
    for (unsigned i = 0; i < QpwqPw.n_cols; i++) {
      QpwqPwc0mQPw.col(i) = QpwqPw.col(i) / c0mQPw;
    }
    arma::vec dll_ccpop_e = trans(sum(QpwqPwc0mQPw, 0));
    
    // Sum up
    double obj = ll_ccpopE - ll_cc - ll_ccpopG;
    arma::vec grad = dll_ccpop_e - dll_cc - dll_ccpopG;

    Rcpp::List returnList;
    returnList["objective"] = obj;
    returnList["gradient"] = Rcpp::wrap(grad);
    return (returnList);
}

// [[Rcpp::export]]
Rcpp::List bn_logistic_glrt_obj_core(arma::colvec beta, double f0, arma::colvec y, arma::colvec g, 
  arma::colvec a, arma::mat ax, arma::mat ax_expand, arma::colvec g_expand, arma::colvec a_expand, 
  double K, double Fst, unsigned N, unsigned n0, unsigned n1, unsigned c0)
{
  unsigned n = ax.n_rows, N_expand = ax_expand.n_rows, n_coeffs_altern = beta.n_rows;

  arma::vec eta = ax * beta;  
  arma::vec exp_eta = exp(eta);
  arma::vec P = exp_eta / (1 + exp_eta);
  arma::vec p = exp_eta / ((1 + exp_eta) % (1 + exp_eta));

  for (unsigned i = 0; i < n; i++) {
    if (y(i) == 0) {
      P(i) = 1 - P(i);
      p(i) = -p(i);
    }
  }

  // The CC-only part    
  double ll_cc = sum(log(P));
  arma::vec dll_cc = join_cols(trans(ax) * (p / P), arma::zeros<arma::vec>(1));

  // The CCpop G|E part
  double tmp_V = Fst * f0 * (1 - f0);
  double d_tmp_V = Fst * (1 - 2 * f0);
  arma::vec q = arma::zeros<arma::vec>(N), Q = arma::zeros<arma::vec>(N);

  for (unsigned i = 0; i < N; i++) {
    switch ((int)g(i)) {
    case 0:
      Q(i) = 1 + f0 * (f0 - 2) + tmp_V + 2 * a(i) * tmp_V * (a(i) - 1);
      q(i) = 2 * f0 - 2 + (1 - 2 * a(i) + 2 * a(i) * a(i)) * d_tmp_V;
      break;
    case 1:
      Q(i) = 2 * (f0 - tmp_V - f0 * f0 + 2 * a(i) * tmp_V * (1 - a(i)));
      q(i) = 2 * (1 - 2 * f0 + (-1 + 2 * a(i) - 2 * a(i) * a(i)) * d_tmp_V);
      break;
    case 2:
      Q(i) = tmp_V + f0 * f0 - 2 * a(i) * tmp_V * (1 - a(i));
      q(i) = 2 * f0 + (1 - 2 * a(i) + 2 * a(i) * a(i)) * d_tmp_V;
      break;
    }
  }

  double ll_ccpop_ge = sum(log(Q));
  arma::vec tmpsc = arma::zeros<arma::vec>(1); tmpsc(0) = sum(q / Q);
  arma::vec dll_ccpop_ge = join_cols(arma::zeros<arma::vec>(n_coeffs_altern), tmpsc);

  // The CCpop E part
  eta = ax_expand * beta;
  exp_eta = exp(eta);
  arma::vec P1 = exp_eta / (1 + exp_eta);
  arma::vec P0 = 1 - P1;
  arma::vec p1 = exp_eta / ((1 + exp_eta) % (1 + exp_eta));
  arma::vec p0 = -p1;
  double w0 = n0 / (1 - K);
  double w1 = n1 / K;
  arma::vec Pw = w0 * P0 + w1 * P1;
  arma::vec pw = w0 * p0 + w1 * p1;

  q = arma::zeros<arma::vec>(N_expand);
  Q = arma::zeros<arma::vec>(N_expand);

  for (unsigned i = 0; i < N_expand; i++) {
    switch ((int)(g_expand(i))) {
    case 0:
      Q(i) = 1 + f0 * (f0 - 2) + tmp_V + 2 * a_expand(i) * tmp_V * (a_expand(i) - 1);
      q(i) = 2 * f0 - 2 + (1 - 2 * a_expand(i) + 2 * a_expand(i) * a_expand(i)) * d_tmp_V;
      break;
    case 1:
      Q(i) = 2 * (f0 - tmp_V - f0 * f0 + 2 * a_expand(i) * tmp_V * (1 - a_expand(i)));
      q(i) = 2 * (1 - 2 * f0 + (-1 + 2 * a_expand(i) - 2 * a_expand(i) * a_expand(i)) * d_tmp_V);
      break;
    case 2:
      Q(i) = tmp_V + f0 * f0 - 2 * a_expand(i) * tmp_V * (1 - a_expand(i));
      q(i) = 2 * f0 + (1 - 2 * a_expand(i) + 2 * a_expand(i) * a_expand(i)) * d_tmp_V;
      break;
    }
  }

  unsigned nr_g_levels = N_expand / N;
  arma::vec QPw = sum(reshape(Q % Pw, N, nr_g_levels), 1);
  arma::vec qPw = sum(reshape(q % Pw, N, nr_g_levels), 1);

  arma::mat Qpw = arma::zeros<arma::mat>(ax_expand.n_rows, ax_expand.n_cols);
  for (unsigned i = 0; i < ax_expand.n_cols; i++) {
    Qpw.col(i) = (Q % pw) % ax_expand.col(i);
  }
  if (nr_g_levels == 2) {
    Qpw = Qpw.rows(0, N-1) + Qpw.rows(N, 2*N-1);
  } else {
    Qpw = Qpw.rows(0, N-1) + Qpw.rows(N, 2*N-1) + Qpw.rows(2*N, 3*N-1);
  }
  arma::vec c0mQPw = c0 + QPw;

  double ll_ccpop_e = sum(log(c0mQPw));
  arma::mat tmp = join_rows(Qpw, qPw);
  tmp.each_col() /= c0mQPw;
  arma::vec dll_ccpop_e = trans(sum(tmp, 0));

  // Sum up
  double obj = ll_ccpop_e - ll_cc - ll_ccpop_ge;
  arma::vec grad = dll_ccpop_e - dll_cc - dll_ccpop_ge;
  
  Rcpp::List returnList;
  returnList["objective"] = obj;
  returnList["gradient"] = Rcpp::wrap(grad);
  return (returnList);
}

// [[Rcpp::export]]
Rcpp::List double_logistic_glrt_obj_core(arma::colvec beta, arma::colvec theta, arma::colvec y, 
  arma::colvec g, arma::mat ax, arma::mat ax_g, arma::mat ax_expand, arma::mat ax_g_expand, 
  arma::colvec g_expand, double w0, double w1)
{
  unsigned n = ax.n_rows, n_coeffs_g = theta.n_rows, N = g.n_rows, N_expand = g_expand.n_rows;
  unsigned nr_g_levels = N_expand / N;
  unsigned nr_g_alleles = (nr_g_levels == 3) ? 2 : 1;
  unsigned n_coeffs_altern = beta.n_rows;
    
  arma::vec eta = ax * beta;  
  arma::vec exp_eta = exp(eta);
  arma::vec P = exp_eta / (1 + exp_eta);
  arma::vec p = exp_eta / ((1 + exp_eta) % (1 + exp_eta));
  
  for (unsigned i = 0; i < n; i++) {
    if (y(i) == 0) {
      P(i) = 1 - P(i);
      p(i) = -p(i);
    }
  }
  
  // The CC-only part
  double ll_cc = sum(log(P));
  arma::vec dll_cc = join_cols(trans(ax) * (p / P), arma::zeros<arma::colvec>(n_coeffs_g));

  // The CCpop G|E part
  arma::vec eta_g = ax_g * theta;  
  arma::vec exp_eta_g = exp(eta_g);
  arma::vec maf_g = exp_eta_g / (1 + exp_eta_g);

  arma::vec q, Q = arma::zeros<arma::vec>(N);
  
  if (nr_g_alleles == 2) {
    q = 2 * exp_eta_g / ((1 + exp_eta_g) % (1 + exp_eta_g));
    
    for (unsigned i = 0; i < N; i++) {
      switch ((int)g(i)) {
      case 0:
        Q(i) = (1 - maf_g(i)) * (1 - maf_g(i));
        q(i) *= (maf_g(i) - 1);
        break;
      case 1:
        Q(i) = 2 * maf_g(i) * (1 - maf_g(i));
        q(i) *= (1 - 2 * maf_g(i));
        break;
      case 2:
        Q(i) = maf_g(i) * maf_g(i);
        q(i) *= maf_g(i);
        break;
      }
    }
  } else {
    q = exp_eta_g / ((1 + exp_eta_g) % (1 + exp_eta_g));
    
    for (unsigned i = 0; i < N; i++) {
      if (g(i) == 0) {
        Q(i) = (1 - maf_g(i));
        q(i) = -q(i);
      } else {
        Q(i) = maf_g(i);
      }
    }
  }

  double ll_ccpop_ge = sum(log(Q));
  arma::vec dll_ccpop_ge = join_cols(arma::zeros<arma::vec>(n_coeffs_altern), trans(ax_g) * (q / Q));

  // The CCpop E part
  eta = ax_expand * beta;
  exp_eta = exp(eta);
  arma::vec P1 = exp_eta / (1 + exp_eta);
  arma::vec P0 = 1 - P1;
  arma::vec p1 = exp_eta / ((1 + exp_eta) % (1 + exp_eta));
  arma::vec p0 = -p1;
  arma::vec Pw = w0 * P0 + w1 * P1;
  arma::vec pw = w0 * p0 + w1 * p1;
  unsigned c0 = N - n;

  eta_g = ax_g_expand * theta;  
  exp_eta_g = exp(eta_g);
  maf_g = exp_eta_g / (1 + exp_eta_g);

  Q = arma::zeros<arma::vec>(N_expand);
  
  if (nr_g_alleles == 2) {
    q = 2 * exp_eta_g / ((1 + exp_eta_g) % (1 + exp_eta_g));
    
    for (unsigned i = 0; i < N_expand; i++) {
      switch ((int)g_expand(i)) {
      case 0:
        Q(i) = (1 - maf_g(i)) * (1 - maf_g(i));
        q(i) *= (maf_g(i) - 1);
        break;
      case 1:
        Q(i) = 2 * maf_g(i) * (1 - maf_g(i));
        q(i) *= (1 - 2 * maf_g(i));
        break;
      case 2:
        Q(i) = maf_g(i) * maf_g(i);
        q(i) *= maf_g(i);
        break;
      }
    }
  } else {
    q = exp_eta_g / ((1 + exp_eta_g) % (1 + exp_eta_g));
    
    for (unsigned i = 0; i < N_expand; i++) {
      if (g_expand(i) == 0) {
        Q(i) = (1 - maf_g(i));
        q(i) = -q(i);
      } else {
        Q(i) = maf_g(i);
      }
    }
  }

  arma::vec QPw = sum(reshape(Q % Pw, N, nr_g_levels), 1);
  arma::mat tmp1 = ax_expand; tmp1.each_col() %= (Q % pw);
  arma::mat tmp2 = ax_g_expand; tmp2.each_col() %= (q % Pw);
  arma::mat Qpw, qPw;
  if (nr_g_levels == 2) {
    Qpw = tmp1.rows(0, N-1) + tmp1.rows(N, 2*N-1);
    qPw = tmp2.rows(0, N-1) + tmp2.rows(N, 2*N-1);
  } else {
    Qpw = tmp1.rows(0, N-1) + tmp1.rows(N, 2*N-1) + tmp1.rows(2*N, 3*N-1);
    qPw = tmp2.rows(0, N-1) + tmp2.rows(N, 2*N-1) + tmp2.rows(2*N, 3*N-1);
  }
  arma::vec c0mQPw = c0 + QPw;
  arma::mat tmp3 = join_rows(Qpw, qPw);
  tmp3.each_col() /= c0mQPw;
  
  double ll_ccpop_e = sum(log(c0mQPw));
  arma::vec dll_ccpop_e = trans(sum(tmp3, 0));

  // Sum up
  double obj = ll_ccpop_e - ll_cc - ll_ccpop_ge;
  arma::vec grad = dll_ccpop_e - dll_cc - dll_ccpop_ge;
  
  Rcpp::List returnList;
  returnList["objective"] = obj;
  returnList["gradient"] = Rcpp::wrap(grad);
  return (returnList);
}
