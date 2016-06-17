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
    // a chgpt between 1 and n exclusive
    IntegerVector ind1 = seq_len(k) - 1, ind2 = seq_len(n) - 1;
    IntegerVector ind3 = setdiff(ind2, ind1);
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

