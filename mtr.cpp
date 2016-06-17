#include <Rcpp.h>
#include <cmath>
#include <Rmath.h>
#include <string>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
const double llik_nb(const NumericVector par, const NumericVector x) {
  // log-likelihood of negative Binomial (NB) model
  double llik, r = par[0], q = par[1];
  if (r <= 0.0 || q <= 0.0 || q > 1.0) {
    llik = -INFINITY;
  }
  else {
    llik = sum(dnbinom(x, r, q, true));
  }
  return llik;
}

// [[Rcpp::export]]
const double llik_nb_rt(const NumericVector par, const NumericVector x, const NumericVector t) {
  // log-likelihood of NB model w/ trend in r (size)
  if (x.size() != t.size()) {
    stop("Lengths of x and t have to be equal.");
  }
  const double r0 = par[0], q = par[1], theta = par[2];
  const NumericVector rt = r0 + theta * t;
  double llik;
  if (is_true(any(rt <= 0.0)) || q <= 0.0 || q > 1.0) {
    llik = -INFINITY;
  }
  else {
    int n = x.size();
    llik = 0.0;
    for (int i = 0; i < n; i++) {
      NumericVector xi = NumericVector::create(x[i]);
      llik += dnbinom(xi, rt[i], q, true)[0];
    }
  }
  return llik;
}

// [[Rcpp::export]]
const double llik_nb_qt(const NumericVector par, const NumericVector x, const NumericVector t) {
  // log-likelihood of NB model w/ trend in q (prob)
  if (x.size() != t.size()) {
    stop("Lengths of x and t have to be equal.");
  }
  const double r = par[0], q0 = par[1], theta = par[2];
  const NumericVector qt = q0 + theta * t;
  double llik;
  if (is_true(any(qt <= 0.0)) || is_true(any(qt > 1.0)) || r <= 0.0) {
    llik = -INFINITY;
  }
  else {
    int n = x.size();
    llik = 0.0;
    for (int i = 0; i < n; i++) {
      NumericVector xi = NumericVector::create(x[i]);
      llik += dnbinom(xi, r, qt[i], true)[0];
    }
  }
  return llik;
}

// [[Rcpp::export]]
const double llik_nb_rk(const double r1, const double r2, const double q, const int k, const NumericVector x) {
  // log-likelihood of NB model w/ chgpt in r (size)
  const int n = x.size();
  double llik;
  if (r1 <= 0.0 || r2 <= 0.0 || q <= 0.0 || q > 1.0 || k < 1 || k > n) {
    llik = -INFINITY;
  }
  else if (k == n) {
    // no chgpt
    llik = sum(dnbinom(x, r1, q, true));
  }
  else {
    // a chgpt less than n
    IntegerVector ind1 = seq_len(k) - 1, 
      ind2 = seq_len(n) - 1,
      ind3 = setdiff(ind2, ind1);
    NumericVector x1 = x[ind1], x2 = x[ind3];
    llik = sum(dnbinom(x1, r1, q, true)) + sum(dnbinom(x2, r2, q, true));
  }
  return llik;
}

// [[Rcpp::export]]
const double llik_nb_rk_fix_k(const NumericVector par, const int k, const NumericVector x) {
  // wrapper of llik_nb_rk() for optim() w/ fixed k
  if (par.size() != 3) {
    stop("par has to be of length 3.");
  }
  const double llik = llik_nb_rk(par[0], par[1], par[2], k, x);
  return llik;
}

// [[Rcpp::export]]
const double llik_nb_qk(const double r, const double q1, const double q2, const int k, const NumericVector x) {
  // log-likelihood of NB model w/ chgpt in q (prob)
  const int n = x.size();
  double llik;
  if (r <= 0.0 || q1 <= 0.0 || q1 > 1.0 || q2 <= 0.0 || q2 > 1.0 || k < 1 || k > n) {
    llik = -INFINITY;
  }
  else if (k == n) {
    // no chgpt
    llik = sum(dnbinom(x, r, q1, true));
  }
  else {
    // a chgpt less than n
    IntegerVector ind1 = seq_len(k) - 1, 
      ind2 = seq_len(n) - 1,
      ind3 = setdiff(ind2, ind1);
    NumericVector x1 = x[ind1], x2 = x[ind3];
    llik = sum(dnbinom(x1, r, q1, true)) + sum(dnbinom(x2, r, q2, true));
  }
  return llik;
}

// [[Rcpp::export]]
const double llik_nb_qk_fix_k(const NumericVector par, const int k, const NumericVector x) {
  // wrapper of llik_nb_qk() for optim() w/ fixed k
  if (par.size() != 3) {
    stop("par has to be of length 3.");
  }
  const double llik = llik_nb_qk(par[0], par[1], par[2], k, x);
  return llik;
}

// [[Rcpp::export]]
const double llik_nb_rqk(const double r1, const double r2, const double q1, const double q2, const int k, const NumericVector x) {
  // log-likelihood of NB model w/ chgpt in r & q simultaneously
  const int n = x.size();
  double llik;
  if (r1 <= 0.0 || r2 <= 0.0 || q1 <= 0.0 || q1 > 1.0 || q2 <= 0.0 || q2 > 1.0 || k < 1 || k > n) {
    llik = -INFINITY;
  }
  else if (k == n) {
    // no chgpt
    llik = sum(dnbinom(x, r1, q1, true));
  }
  else {
    // a chgpt less than n
    IntegerVector ind1 = seq_len(k) - 1, 
      ind2 = seq_len(n) - 1,
      ind3 = setdiff(ind2, ind1);
    NumericVector x1 = x[ind1], x2 = x[ind3];
    llik = sum(dnbinom(x1, r1, q1, true)) + sum(dnbinom(x2, r2, q2, true));
  }
  return llik;
}

// [[Rcpp::export]]
const double llik_nb_rqk_fix_k(const NumericVector par, const int k, const NumericVector x) {
  // wrapper of llik_nb_rqk() for optim() w/ fixed k
  if (par.size() != 4) {
    stop("par has to be of length 4.");
  }
  const double llik = llik_nb_rqk(par[0], par[1], par[2], par[3], k, x);
  return llik;
}

// [[Rcpp::export]]
const double llik_nb_rkqk(const double r1, const double r2, const double q1, const double q2, const double kr, const double kq, const NumericVector x) {
  // log-likelihood of NB model w/ chgpt in r & q separately
  const int n = x.size();
  double llik;
  if (r1 <= 0.0 || r2 <= 0.0 || q1 <= 0.0 || q1 > 1.0 || q2 <= 0.0 || q2 > 1.0 || kr < 1 || kr > n || kq < 1 || kq > n) {
    llik = -INFINITY;
  }
  else if (kr == n && kq == n) { 
    // no chgpt for both r & q
    llik = sum(dbinom(x, r1, q1, true));
  }
  else if (kr == n) {
    // no chgpt for r
    IntegerVector ind1 = seq_len(kq) - 1,
      ind2 = seq_len(n) - 1,
      ind3 = setdiff(ind2, ind1);
    NumericVector x1 = x[ind1], x2 = x[ind3];
    llik = sum(dnbinom(x1, r1, q1, true)) + sum(dnbinom(x2, r1, q2, true));
  }
  else if (kq == n) {
    // no chgpt for q
    IntegerVector ind1 = seq_len(kr) - 1,
      ind2 = seq_len(n) - 1,
      ind3 = setdiff(ind2, ind1);
    NumericVector x1 = x[ind1], x2 = x[ind3];
    llik = sum(dnbinom(x1, r1, q1, true)) + sum(dnbinom(x2, r2, q1, true));
  }
  else if (kr == kq) {
    // both chgpts coincide and are less than n
    IntegerVector ind1 = seq_len(kr) - 1,
      ind2 = seq_len(n) - 1,
      ind3 = setdiff(ind2, ind1);
    NumericVector x1 = x[ind1], x2 = x[ind3];
    llik = sum(dnbinom(x1, r1, q1, true)) + sum(dnbinom(x2, r2, q2, true));
  }
  else if (kr < kq) {
    // chgpt for r before chgpt for q, both chgpts smaller than n
    IntegerVector ind1 = seq_len(kr) - 1,
      ind2 = seq_len(kq) - 1,
      ind3 = seq_len(n) - 1,
      ind4 = setdiff(ind2, ind1),
      ind5 = setdiff(ind3, ind2);
    NumericVector x1 = x[ind1], x2 = x[ind4], x3 = x[ind5];
    llik = sum(dnbinom(x1, r1, q1, true)) + sum(dnbinom(x2, r2, q1, true)) + sum(dnbinom(x3, r2, q2, true));
  }
  else if (kr > kq) {
    // chgpt for q before chgpt for r, both chgpts smaller than n
    IntegerVector ind1 = seq_len(kq) - 1,
      ind2 = seq_len(kr) - 1,
      ind3 = seq_len(n) - 1,
      ind4 = setdiff(ind2, ind1),
      ind5 = setdiff(ind3, ind2);
    NumericVector x1 = x[ind1], x2 = x[ind4], x3 = x[ind5];
    llik = sum(dnbinom(x1, r1, q1, true)) + sum(dnbinom(x2, r1, q2, true)) + sum(dnbinom(x3, r2, q2, true));
  }
  return llik;
}

// [[Rcpp::export]]
const double llik_nb_rkqk_fix_k(const NumericVector par, const int kr, const int kq, const NumericVector x) {
  // wrapper of llik_nb_rkqk() for optim() w/ fixed kr & kq
  if (par.size() != 4) {
    stop("par has to be of length 4.");
  }
  const double llik = llik_nb_rkqk(par[0], par[1], par[2], par[3], kr, kq, x);
  return llik;
}
