// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>
#include <cmath>
#include <Rmath.h>
#include <string>
using namespace Rcpp;
using namespace std;
using namespace arma;

/*
  For all changepoint models, we assume the following priors:
  r ~ Gamma(a, b), r1 ~ Gamma(a1, b1), r2 ~ Gamma(a2, b2) (all default to 0.01)
  q, q1, q2 ~ Unif(0, 1)
  k, kr, kq ~ discrete Uniform[1, n]
*/ 

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
const double llik_nb_rqk(const double r1, const double r2, const double q1, const double q2, const int k, const NumericVector x) {
  // log-likelihood of NB model w/ chgpt in r & q simultaneously
  const int n = x.size();
  double llik;
  if (r1 <= 0.0 || r2 <= 0.0 || q1 <= 0.0 || q1 > 1.0 || q2 <= 0.0 || q2 > 1.0 || k < 1 || k >= n) {
    llik = -INFINITY;
  }
  else {
    // a chgpt less than n
    IntegerVector indn = seq_len(n) - 1,
      ind1 = seq_len(k) - 1,
      ind2 = setdiff(indn, ind1);
    NumericVector x1 = x[ind1], x2 = x[ind2];
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
const double lpost_nb_rqk(const double r1, const double r2, const double q1, const double q2, const int k, const NumericVector x, const double a1 = 0.01, const double b1 = 0.01, const double a2 = 0.01, const double b2 = 0.01) {
  // log-posterior of NB model w/ chgpt in r & q simultaneously
  double lpost = llik_nb_rqk(r1, r2, q1, q2, k, x) + dgamma(NumericVector::create(r1), a1, 1.0 / b1, true)[0] + dgamma(NumericVector::create(r2), a2, 1.0 / b2, true)[0];
  if (lpost != lpost) {
    lpost = -INFINITY;
  }
  return lpost;
}

// [[Rcpp::export]]
const double llik_r(const double r, const double beta, const double lam, const double a, const double b) {
  double llik;
  if (r <= 0.0 || beta <= 0.0 || lam <= 0.0 || a <= 0.0 || b <= 0.0) {
    llik = -INFINITY;
  }
  else{
    llik = r * log(beta) + (r - 1.0) * log(lam) + (a - 1.0) * log(r) - b * r - lgamma(r);
  }
  return llik;
}

// [[Rcpp::export]]
const NumericMatrix mwg_nb_lamk(const double r1, const double r2, const double beta1, const double beta2, const int k, const NumericVector x, const double s_r1, const double s_r2, const int N = 1e+6, const int thin = 1, const int burnin = 0, const double a1 = 0.01, const double b1 = 0.01, const double c1 = 0.01, const double d1 = 0.01, const double a2 = 0.01, const double b2 = 0.01, const double c2 = 0.01, const double d2 = 0.01) {
  // Metropolis-within-Gibbs for NB-as-Poisson-Gamma-mixture model
  const int n = x.size();
  double lam1_curr, lam2_curr, r1_curr = r1, r1_prop, r2_curr = r2, r2_prop, beta1_curr = beta1, beta2_curr = beta2, log_alpha;
  int k_curr = k;
  NumericMatrix par_mat(N, 7);
  const IntegerVector indn = seq_len(n) - 1, seq_k = seq_len(n - 1);
  IntegerVector ind1, ind2;
  NumericVector x1, x2, cx = cumsum(x), exponent(n - 1), seq_unscaled(n - 1), seq_scaled(n - 1);
  cx = cx[seq_k - 1];
  int i, j;
  RNGScope scope;
  for (i  = 0; i < N * thin + burnin; i++) {
    ind1 = seq_len(k_curr) - 1;
    ind2 = setdiff(indn, ind1);
    x1 = x[ind1];
    x2 = x[ind2];
    lam1_curr = rgamma(1, r1_curr + sum(x1), 1.0 / (beta1_curr + (double) k_curr))[0]; // 1
    lam2_curr = rgamma(1, r2_curr + sum(x2), 1.0 / (beta2_curr + (double) (n - k_curr)))[0]; // 2
    exponent = (lam2_curr - lam1_curr) * (NumericVector) seq_k + (log(lam1_curr) - log(lam2_curr)) * cx;
    seq_unscaled = exp(exponent - max(exponent));
    seq_unscaled = ifelse(seq_unscaled != seq_unscaled, 0.0, seq_unscaled); // underflow gives nan
    seq_scaled = seq_unscaled / sum(seq_unscaled); // the probabilities
    k_curr = Rcpp::RcppArmadillo::sample(seq_k, 1, false, seq_scaled)[0]; // 3
    beta1_curr = rgamma(1, r1_curr + c1, 1.0 / (lam1_curr + d1))[0]; // 4
    beta2_curr = rgamma(1, r2_curr + c2, 1.0 / (lam2_curr + d2))[0]; // 5
    r1_prop = exp(rnorm(1, log(r1_curr), s_r1))[0]; // 6
    log_alpha = llik_r(r1_prop, beta1_curr, lam1_curr, a1, b1) - llik_r(r1_curr, beta1_curr, lam1_curr, a1, b1) + log(r1_prop) - log(r1_curr);
    if (log(runif(1))[0] < log_alpha) {
      r1_curr = r1_prop;
    }
    r2_prop = exp(rnorm(1, log(r2_curr), s_r2))[0]; // 7
    log_alpha = llik_r(r2_prop, beta2_curr, lam2_curr, a2, b2) - llik_r(r2_curr, beta2_curr, lam2_curr, a2, b2) + log(r2_prop) - log(r2_curr);
    if (log(runif(1))[0] < log_alpha) {
      r2_curr = r2_prop;
    }
    if (i >= burnin && (i - burnin + 1) % thin == 0) {
      j = (i - burnin + 1) / thin - 1;
      par_mat(j, _) = NumericVector::create(lam1_curr, lam2_curr, beta1_curr, beta2_curr, r1_curr, r2_curr, (double) k_curr);
      if ((j + 1) % 100 == 0) {
        Rcout << j + 1 << endl;
      }
    }
  }
  return par_mat;
}

// [[Rcpp::export]]
const double llik_r_beta(const double r, const double beta, const double lam, const double a, const double b, const double c, const double d) {
  double llik;
  if (r <= 0.0 || beta <= 0.0 || lam <= 0.0 || a <= 0.0 || b <= 0.0 || c <= 0.0 || d <= 0.0) {
    llik = -INFINITY;
  }
  else{
    llik = r * log(beta) + (r - 1.0) * log(lam) + (a - 1.0) * log(r) - b * r - lgamma(r) + (c - 1.0) * log(beta) - beta * d;
  }
  return llik;
}

// [[Rcpp::export]]
const NumericMatrix mwg_nb_lamk_block(const double r1, const double r2, const double beta1, const double beta2, const int k, const NumericVector x, const double s_r1, const double s_r2, const double s_beta1, const double s_beta2, const double rho1, const double rho2, const int N = 1e+6, const int thin = 1, const int burnin = 0, const double a1 = 0.01, const double b1 = 0.01, const double c1 = 0.01, const double d1 = 0.01, const double a2 = 0.01, const double b2 = 0.01, const double c2 = 0.01, const double d2 = 0.01) {
  // Metropolis-within-Gibbs for NB-as-Poisson-Gamma-mixture model, with block updating for (beta1, r1) & (beta2, r2)
  const int n = x.size();
  double lam1_curr, lam2_curr, r1_curr = r1, r1_prop, r2_curr = r2, r2_prop, beta1_curr = beta1, beta1_prop, beta2_curr = beta2, beta2_prop, log_alpha;
  int k_curr = k;
  NumericMatrix par_mat(N, 7);
  const IntegerVector indn = seq_len(n) - 1, seq_k = seq_len(n - 1);
  IntegerVector ind1, ind2;
  NumericVector x1, x2, cx = cumsum(x), exponent(n - 1), seq_unscaled(n - 1), seq_scaled(n - 1);
  cx = cx[seq_k - 1];
  int i, j;
  RNGScope scope;
  for (i  = 0; i < N * thin + burnin; i++) {
    ind1 = seq_len(k_curr) - 1;
    ind2 = setdiff(indn, ind1);
    x1 = x[ind1];
    x2 = x[ind2];
    lam1_curr = rgamma(1, r1_curr + sum(x1), 1.0 / (beta1_curr + (double) k_curr))[0]; // 1
    lam2_curr = rgamma(1, r2_curr + sum(x2), 1.0 / (beta2_curr + (double) (n - k_curr)))[0]; // 2
    exponent = (lam2_curr - lam1_curr) * (NumericVector) seq_k + (log(lam1_curr) - log(lam2_curr)) * cx;
    seq_unscaled = exp(exponent - max(exponent));
    seq_unscaled = ifelse(seq_unscaled != seq_unscaled, 0.0, seq_unscaled); // underflow gives nan
    seq_scaled = seq_unscaled / sum(seq_unscaled); // the probabilities
    k_curr = Rcpp::RcppArmadillo::sample(seq_k, 1, false, seq_scaled)[0]; // 3
    beta1_prop = rnorm(1, beta1_curr, s_beta1)[0]; // 4
    r1_prop = rnorm(1, (r1_curr + s_r1 / s_beta1 * rho1 * (beta1_prop - beta1_curr)), sqrt(1.0 - rho1 * rho1) * s_r1)[0]; // 5
    log_alpha = llik_r_beta(r1_prop, beta1_prop, lam1_curr, a1, b1, c1, d1) - llik_r_beta(r1_curr, beta1_curr, lam1_curr, a1, b1, c1, d1);
    if (log(runif(1))[0] < log_alpha) {
      r1_curr = r1_prop;
      beta1_prop = beta1_prop;
    }
    beta2_prop = rnorm(1, beta2_curr, s_beta2)[0]; // 6
    r2_prop = rnorm(1, (r2_curr + s_r2 / s_beta2 * rho2 * (beta2_prop - beta2_curr)), sqrt(1.0 - rho2 * rho2) * s_r2)[0]; // 7
    log_alpha = llik_r_beta(r2_prop, beta2_prop, lam2_curr, a2, b2, c2, d2) - llik_r_beta(r2_curr, beta2_curr, lam2_curr, a2, b2, c2, d2);
    if (log(runif(1))[0] < log_alpha) {
      r2_curr = r2_prop;
      beta2_prop = beta2_prop;
    }
    if (i >= burnin && (i - burnin + 1) % thin == 0) {
      j = (i - burnin + 1) / thin - 1;
      par_mat(j, _) = NumericVector::create(lam1_curr, lam2_curr, beta1_curr, beta2_curr, r1_curr, r2_curr, (double) k_curr);
      if ((j + 1) % 100 == 0) {
        Rcout << j + 1 << endl;
      }
    }
  }
  return par_mat;
}



