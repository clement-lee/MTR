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
  if (par.size() != 2) {
    stop("par has to be of length 2.");
  }
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
const double lpost_nb(const double r, const double q, const NumericVector x, const double a = 0.01, const double b = 0.01) {
  // log-posterior density of NB model
  const NumericVector par = NumericVector::create(r, q);
  double llik = llik_nb(par, x) + dgamma(NumericVector::create(r), a, 1.0 / b, true)[0];
  if (llik != llik) { // NaN
    llik = -INFINITY;
  }
  return llik;
}

// [[Rcpp::export]]
List mwg_nb(const double r, const double q, const NumericVector x, const double s_r, const int N = 1e+6, const int thin = 1, const int burnin = 0, const double a = 0.01, const double b = 0.01) {
  // Metropolis-within-Gibbs for NB model
  const int n = x.size();
  if (r <= 0.0 || q <= 0.0) {
    stop("All parameters have to be positive.");
  }
  double r_curr = r, r_prop, q_curr = q, log_alpha;
  NumericVector r_vec(N), q_vec(N), llik_vec(N), par(2);
  int i, j;
  RNGScope scope;
  for (i = 0; i < N * thin + burnin; i++) {
    // update r
    r_prop = exp(rnorm(1, log(r_curr), s_r))[0];
    log_alpha = lpost_nb(r_prop, q_curr, x, a, b) - lpost_nb(r_curr, q_curr, x, a, b) + log(r_prop) - log(r_curr);
    if (log(runif(1))[0] < log_alpha) {
      r_curr = r_prop;
    }
    // update q
    q_curr = rbeta(1, (double) n * r_curr + 1.0, sum(x) + 1.0)[0];
    if (i >= burnin && (i - burnin + 1) % thin == 0) {
      j = (i - burnin + 1) / thin - 1;
      r_vec[j] = r_curr;
      q_vec[j] = q_curr;
      llik_vec[j] = llik_nb(NumericVector::create(r_curr, q_curr), x);
      if ((j + 1) % 100 == 0) {
        Rcout << j + 1 << endl;
      }
    }
  }
  DataFrame output = DataFrame::create(Named("r") = r_vec, Named("q") = q_vec, Named("llik") = llik_vec);
  return output;
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
const double llik_nb_rqk(const NumericVector par, const NumericVector x) {
  // log-likelihood of NB model w/ chgpt in r & q simultaneously
  if (par.size() != 5) {
    stop("par has to be length 5.");
  }
  const int n = x.size(), k = (int) par[4];
  double llik;
  NumericVector par1 = NumericVector::create(par[0], par[1]), par2 = NumericVector::create(par[2], par[3]);
  if (k < 1 || k > n) {
    llik = -INFINITY;
  }
  else if (k == n) {
    llik = llik_nb(par1, x);
  }
  else {
    NumericVector x1 = head(x, k), x2 = tail(x, n - k);
    llik = llik_nb(par1, x1) + llik_nb(par2, x2);
  }
  return llik;
}

// [[Rcpp::export]]
const double llik_nb_rqk_fix_k(const NumericVector par, const int k, const NumericVector x) {
  // wrapper of llik_nb_rqk() for optim() w/ fixed k
  if (par.size() != 4) {
    stop("par has to be of length 4.");
  }
  const NumericVector par0 = NumericVector::create(par[0], par[1], par[2], par[3], k);
  const double llik = llik_nb_rqk(par0, x);
  return llik;
}

// [[Rcpp::export]]
List mwg_nb_rqk(const double r1, const double r2, const double q1, const double q2, const int k, const NumericVector x, const double s_r1, const double s_r2, const int N = 1e+6, const int thin = 1, const int burnin = 0, const double a1 = 0.01, const double b1 = 0.01, const double a2 = 0.01, const double b2 = 0.01) {
  // Metropolis-within-Gibbs for NB model w/ chgpt in r & q simultaneously
  const int n = x.size();
  if (k < 1 || k >= n) {
    stop("k has to be a positive integer smaller than n.");
  }
  if (r1 <= 0.0 || r2 <= 0.0 || q1 <= 0.0 || q2 <= 0.0) {
    stop("All parameters have to be positive.");
  }
  double r1_curr = r1, r1_prop, r2_curr = r2, r2_prop, q1_curr = q1, q2_curr = q2, log_alpha;
  int k_curr = k;
  NumericVector r1_vec(N), r2_vec(N), q1_vec(N), q2_vec(N), llik_vec(N), par1(2), par2(2);
  IntegerVector k_vec(N);
  const IntegerVector seq__n = seq_len(n - 1); // 1 to (n - 1); *__n means excluding the nth element
  NumericVector x1, x2, x__n = head(x, n - 1), lgamma1(n - 1), lgamma2(n - 1), seq_unscaled(n - 1), seq_scaled(n - 1);
  int i, j;
  RNGScope scope;
  for (i = 0; i < N * thin + burnin; i++) {
    // update k
    lgamma1 = lgamma(x__n + r1_curr) - lgamma(r1_curr) + r1_curr * log(q1_curr) + x__n * log(1.0 - q1_curr);
    lgamma2 = lgamma(x__n + r2_curr) - lgamma(r2_curr) + r2_curr * log(q2_curr) + x__n * log(1.0 - q2_curr);
    NumericVector exponent = cumsum(lgamma1 - lgamma2); // have to define everytime
    seq_unscaled = exp(exponent - max(exponent));
    seq_unscaled = ifelse(seq_unscaled != seq_unscaled, 0.0, seq_unscaled); // underflow gives nan
    seq_scaled = seq_unscaled / sum(seq_unscaled); // the probabilities
    k_curr = Rcpp::RcppArmadillo::sample(seq__n, 1, false, seq_scaled)[0];
    // update r1 & r2
    x1 = head(x, k_curr);
    r1_prop = exp(rnorm(1, log(r1_curr), s_r1))[0];
    log_alpha = lpost_nb(r1_prop, q1_curr, x1, a1, b1) - lpost_nb(r1_curr, q1_curr, x1, a1, b1) + log(r1_prop) - log(r1_curr);
    if (log(runif(1))[0] < log_alpha) {
      r1_curr = r1_prop;
    }
    x2 = tail(x, n - k_curr);
    r2_prop = exp(rnorm(1, log(r2_curr), s_r2))[0];
    log_alpha = lpost_nb(r2_prop, q2_curr, x2, a2, b2) - lpost_nb(r2_curr, q2_curr, x2, a2, b2) + log(r2_prop) - log(r2_curr);
    if (log(runif(1))[0] < log_alpha) {
      r2_curr = r2_prop;
    }
    // update q1 & q2
    q1_curr = rbeta(1, (double) k_curr * r1_curr + 1.0, sum(x1) + 1.0)[0];
    q2_curr = rbeta(1, (double) (n - k_curr) * r2_curr + 1.0, sum(x2) + 1.0)[0];
    if (i >= burnin && (i - burnin + 1) % thin == 0) {
      j = (i - burnin + 1) / thin - 1;
      r1_vec[j] = r1_curr;
      r2_vec[j] = r2_curr;
      q1_vec[j] = q1_curr;
      q2_vec[j] = q2_curr;
      k_vec[j] = k_curr;
      par1 = NumericVector::create(r1_curr, q1_curr);
      par2 = NumericVector::create(r2_curr, q2_curr);
      llik_vec[j] = llik_nb(par1, x1) + llik_nb(par2, x2);
      if ((j + 1) % 100 == 0) {
        Rcout << j + 1 << endl;
      }
    }
  }
  DataFrame output = DataFrame::create(Named("r1") = r1_vec, Named("r2") = r2_vec, Named("q1") = q1_vec, Named("q2") = q2_vec, Named("k") = k_vec, Named("llik") = llik_vec);
  return output;
}
