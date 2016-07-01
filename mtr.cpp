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
    IntegerVector indn = seq_len(n) - 1,
      ind1 = seq_len(k) - 1,
      ind2 = setdiff(indn, ind1);
    NumericVector x1 = x[ind1], x2 = x[ind2];
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
const double lpost_nb_rk(const double r1, const double r2, const double q, const int k, const NumericVector x, const double a1 = 0.01, const double b1 = 0.01, const double a2 = 0.01, const double b2 = 0.01) {
  // log-posterior of NB model w/ chgpt in r (size)
  double lpost = llik_nb_rk(r1, r2, q, k, x) + dgamma(NumericVector::create(r1), a1, 1.0 / b1, true)[0] + dgamma(NumericVector::create(r2), a2, 1.0 / b2, true)[0];
  // uniform priors for q & k won't matter
  if (lpost != lpost) { // NaN
    // this happens when r is outside support
    lpost = -INFINITY;
  }
  return lpost;
}

// [[Rcpp::export]]
List rwm_nb_rk(const double r1, const double r2, const double q, const int k, const NumericVector x, const double s_r1, const double s_r2, const double s_k, const int N = 1e+6, const int thin = 1, const int burnin = 0, const double a1 = 0.01, const double b1 = 0.01, const double a2 = 0.01, const double b2 = 0.01) {
  // componentwise RWM of NB model w/ chgpt in r (size)
  double r1_curr = r1, r1_prop, r2_curr = r2, r2_prop, q_curr = q, k_curr = k, k_prop;
  double lpost_curr = lpost_nb_rk(r1_curr, r2_curr, q_curr, k_curr, x, a1, b1, a2, b2), lpost_prop, log_alpha;
  NumericMatrix par_mat(N, 4);
  NumericVector lpost_vec(N), log_u(3);
  const int n = x.size();
  int i, j;
  RNGScope scope;
  for (i = 0; i < N * thin + burnin; i++) {
    log_u = log(runif(3)); // for acceptance-rejection
    // update r1
    r1_prop = exp(rnorm(1, log(r1_curr), s_r1)[0]);
    lpost_prop = lpost_nb_rk(r1_prop, r2_curr, q_curr, k_curr, x, a1, b1, a2, b2);
    log_alpha = lpost_prop - lpost_curr + log(r1_prop) - log(r1_curr);
    if (log_u[0] < log_alpha) {
      r1_curr = r1_prop;
      lpost_curr = lpost_prop;
    }
    // update r2
    r2_prop = exp(rnorm(1, log(r2_curr), s_r2)[0]);
    lpost_prop = lpost_nb_rk(r1_curr, r2_prop, q_curr, k_curr, x, a1, b1, a2, b2);
    log_alpha = lpost_prop - lpost_curr + log(r2_prop) - log(r2_curr);
    if (log_u[1] < log_alpha) {
      r2_curr = r2_prop;
      lpost_curr = lpost_prop;
    }
    // update q
    q_curr = rbeta(1, k_curr * r1_curr + (n - k_curr) * r2_curr + 1.0, sum(x) + 1.0)[0];
    lpost_curr = lpost_nb_rk(r1_curr, r2_curr, q_curr, k_curr, x, a1, b1, a2, b2);
    // update k
    k_prop = (int) round(rnorm(1, k_curr + 0.0, s_k)[0]);
    lpost_prop = lpost_nb_rk(r1_curr, r2_curr, q_curr, k_prop, x, a1, b1, a2, b2);
    log_alpha = lpost_prop - lpost_curr;
    if (log_u[2] < log_alpha) {
      k_curr = k_prop;
      lpost_curr = lpost_prop;
    }
    // save for output
    if (i >= burnin && (i - burnin + 1) % thin == 0) {
      j = (i - burnin + 1) / thin - 1;
      par_mat(j, _) = NumericVector::create(r1_curr, r2_curr, q_curr, k_curr);
      lpost_vec[j] = lpost_curr;
    }
  }
  // return
  List L;
  L["par"] = par_mat;
  L["lpost"] = lpost_vec;
  return L;
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
    IntegerVector indn = seq_len(n) - 1,
      ind1 = seq_len(k) - 1,
      ind2 = setdiff(indn, ind1);
    NumericVector x1 = x[ind1], x2 = x[ind2];
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
const double lpost_nb_qk(const double r, const double q1, const double q2, const int k, const NumericVector x, const double a = 0.01, const double b = 0.01) {
  // log-posterior of NB model w/ chgpt in q (prob)
  double lpost = llik_nb_qk(r, q1, q2, k, x) + dgamma(NumericVector::create(r), a, 1.0 / b, true)[0];
  if (lpost != lpost) {
    lpost = -INFINITY;
  }
  return lpost;
}

// [[Rcpp::export]]
List rwm_nb_qk(const double r, const double q1, const double q2, const int k, const NumericVector x, const double s_r, const double s_k, const int N = 1e+6, const int thin = 1, const int burnin = 0, const double a = 0.01, const double b = 0.01) {
  // componentwise RWM of NB model w/ chgpt in q (prob)
  double r_curr = r, r_prop, q1_curr = q1, q2_curr = q2, k_curr = k, k_prop;
  double lpost_curr = lpost_nb_qk(r_curr, q1_curr, q2_curr, k_curr, x, a, b), lpost_prop, log_alpha;
  NumericMatrix par_mat(N, 4);
  NumericVector lpost_vec(N), log_u(2), x1, x2;
  const int n = x.size();
  const IntegerVector indn = seq_len(n) - 1;
  IntegerVector ind1, ind2;
  int i, j;
  RNGScope scope;
  for (i = 0; i < N * thin + burnin; i++) {
    log_u = log(runif(2));
    // update r
    r_prop = exp(rnorm(1, log(r_curr), s_r)[0]);
    lpost_prop = lpost_nb_qk(r_prop, q1_curr, q2_curr, k_curr, x, a, b);
    log_alpha = lpost_prop - lpost_curr + log(r_prop) - log(r_curr);
    if (log_u[0] < log_alpha) {
      r_curr = r_prop;
      lpost_curr = lpost_prop;
    }
    // update q1
    ind1 = seq_len(k_curr) - 1;
    x1 = x[ind1];
    q1_curr = rbeta(1, k_curr * r_curr + 1.0, sum(x1) + 1.0)[0];
    lpost_curr = lpost_nb_qk(r_curr, q1_curr, q2_curr, k_curr, x, a, b);
    // update q2
    ind2 = setdiff(indn, ind1);
    x2 = x[ind2];
    q2_curr = rbeta(1, (n - k_curr) * r_curr + 1.0, sum(x2) + 1.0)[0];
    lpost_curr = lpost_nb_qk(r_curr, q1_curr, q2_curr, k_curr, x, a, b);
    // update k
    k_prop = (int) round(rnorm(1, k_curr + 0.0, s_k)[0]);
    lpost_prop = lpost_nb_qk(r_curr, q1_curr, q2_curr, k_prop, x, a, b);
    log_alpha = lpost_prop - lpost_curr;
    if (log_u[1] < log_alpha) {
      k_curr = k_prop;
      lpost_curr = lpost_prop;
    }
    // save for output
    if (i >= burnin && (i - burnin + 1) % thin == 0) {
      j = (i - burnin + 1) / thin - 1;
      par_mat(j, _) = NumericVector::create(r_curr, q1_curr, q2_curr, k_curr);
      lpost_vec[j] = lpost_curr;
    }
  }
  // return
  List L;
  L["par"] = par_mat;
  L["lpost"] = lpost_vec;
  return L;
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
List rwm_nb_rqk(const double r1, const double r2, const double q1, const double q2, const int k, const NumericVector x, const double s_r1, const double s_r2, const double s_k, const int N = 1e+6, const int thin = 1, const int burnin = 0, const double a1 = 0.01, const double b1 = 0.01, const double a2 = 0.01, const double b2 = 0.01) {
  // componentwise RWM of NB model w/ chgpt in r & q simultaneously
  double r1_curr = r1, r1_prop, r2_curr = r2, r2_prop, q1_curr = q1, q2_curr = q2, k_curr = k, k_prop;
  double lpost_curr = lpost_nb_rqk(r1_curr, r2_curr, q1_curr, q2_curr, k_curr, x, a1, b1, a2, b2), lpost_prop, log_alpha;
  NumericMatrix par_mat(N, 5);
  NumericVector lpost_vec(N), log_u(3), x1, x2;
  const int n = x.size();
  const IntegerVector indn = seq_len(n) - 1;
  IntegerVector ind1, ind2;
  int i, j;
  RNGScope scope;
  for (i = 0; i < N * thin + burnin; i++) {
    log_u = log(runif(3));
    // update r1
    r1_prop = exp(rnorm(1, log(r1_curr), s_r1)[0]);
    lpost_prop = lpost_nb_rqk(r1_prop, r2_curr, q1_curr, q2_curr, k_curr, x, a1, b1, a2, b2);
    log_alpha = lpost_prop - lpost_curr + log(r1_prop) - log(r1_curr);
    if (log_u[0] < log_alpha) {
      r1_curr = r1_prop;
      lpost_curr = lpost_prop;
    }
    // update r2
    r2_prop = exp(rnorm(1, log(r2_curr), s_r2)[0]);
    lpost_prop = lpost_nb_rqk(r1_curr, r2_prop, q1_curr, q2_curr, k_curr, x, a1, b1, a2, b2);
    log_alpha = lpost_prop - lpost_curr + log(r2_prop) - log(r2_curr);
    if (log_u[1] < log_alpha) {
      r2_curr = r2_prop;
      lpost_curr = lpost_prop;
    }
    // update q1
    ind1 = seq_len(k_curr) - 1;
    x1 = x[ind1];
    q1_curr = rbeta(1, k_curr * r1_curr + 1.0, sum(x1) + 1.0)[0];
    lpost_curr = lpost_nb_rqk(r1_curr, r2_curr, q1_curr, q2_curr, k_curr, x, a1, b1, a2, b2);
    // update q2
    ind2 = setdiff(indn, ind1);
    x2 = x[ind2];
    q2_curr = rbeta(1, (n - k_curr) * r2_curr + 1.0, sum(x2) + 1.0)[0];
    lpost_curr = lpost_nb_rqk(r1_curr, r2_curr, q1_curr, q2_curr, k_curr, x, a1, b1, a2, b2);
    // update k
    k_prop = (int) round(rnorm(1, k_curr + 0.0, s_k)[0]);
    lpost_prop = lpost_nb_rqk(r1_curr, r2_curr, q1_curr, q2_curr, k_prop, x, a1, b1, a2, b2);
    log_alpha = lpost_prop - lpost_curr;
    if (log_u[2] < log_alpha) {
      k_curr = k_prop;
      lpost_curr = lpost_prop;
    }
    // save for output
    if (i >= burnin && (i - burnin + 1) % thin == 0) {
      j = (i - burnin + 1) / thin - 1;
      par_mat(j, _) = NumericVector::create(r1_curr, r2_curr, q1_curr, q2_curr, k_curr);
      lpost_vec[j] = lpost_curr;
    }
  }
  // return
  List L;
  L["par"] = par_mat;
  L["lpost"] = lpost_vec;
  return L;
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

// [[Rcpp::export]]
const double lpost_nb_rkqk(const double r1, const double r2, const double q1, const double q2, const int kr, const int kq, const NumericVector x, const double a1 = 0.01, const double b1 = 0.01, const double a2 = 0.01, const double b2 = 0.01) {
  // log-likelihood of NB model w/ chgpt in r & q separately
  double lpost = llik_nb_rkqk(r1, r2, q1, q2, kr, kq, x) + dgamma(NumericVector::create(r1), a1, 1.0 / b1, true)[0] + dgamma(NumericVector::create(r2), a2, 1.0 / b2, true)[0];
  if (lpost != lpost) {
    lpost = -INFINITY;
  }
  return lpost;
}

// [[Rcpp::export]]
List rwm_nb_rkqk(const double r1, const double r2, const double q1, const double q2, const int kr, const int kq, const NumericVector x, const double s_r1, const double s_r2, const double s_kr, const double s_kq, const int N = 1e+6, const int thin = 1, const int burnin = 0, const double a1 = 0.01, const double b1 = 0.01, const double a2 = 0.01, const double b2 = 0.01) {
  // componentwise RWM of NB model w/ chgpt in r & q separately
  double r1_curr = r1, r1_prop, r2_curr = r2, r2_prop, q1_curr = q1, q2_curr = q2, kr_curr = kr, kr_prop, kq_curr = kq, kq_prop;
  double lpost_curr = lpost_nb_rkqk(r1, r2, q1, q2, kr, kq, x, a1, b1, a2, b2), lpost_prop, log_alpha, aq1, aq2;
  NumericMatrix par_mat(N, 6);
  NumericVector lpost_vec(N), log_u(4), x1, x2, x3;
  const int n = x.size();
  const IntegerVector indn = seq_len(n) - 1;
  IntegerVector ind1, ind2, ind3;
  int i, j;
  RNGScope scope;
  for (i = 0; i < N * thin + burnin; i++) {
    log_u = log(runif(4));
    // update r1
    r1_prop = exp(rnorm(1, log(r1_curr), s_r1)[0]);
    lpost_prop = lpost_nb_rkqk(r1_prop, r2_curr, q1_curr, q2_curr, kr_curr, kq_curr, x, a1, b1, a2, b2);
    log_alpha = lpost_prop - lpost_curr + log(r1_prop) - log(r1_curr);
    if (log_u[0] < log_alpha) {
      r1_curr = r1_prop;
      lpost_curr = lpost_prop;
    }
    // update r2
    r2_prop = exp(rnorm(1, log(r2_curr), s_r2)[0]);
    lpost_prop = lpost_nb_rkqk(r1_curr, r2_prop, q1_curr, q2_curr, kr_curr, kq_curr, x, a1, b1, a2, b2);
    log_alpha = lpost_prop - lpost_curr + log(r2_prop) - log(r2_curr);
    if (log_u[1] < log_alpha) {
      r2_curr = r2_prop;
      lpost_curr = lpost_prop;
    }
    // update q1
    ind1 = seq_len(kq_curr) - 1;
    x1 = x[ind1];
    if (kr_curr <= kq_curr) {
      aq1 = kr_curr * r1_curr + (kq_curr - kr_curr) * r2_curr + 1.0;
    }
    else {
      aq1 = kq_curr * r1_curr + 1.0;
    }
    q1_curr = rbeta(1, aq1, sum(x1) + 1.0)[0];
    lpost_curr = lpost_nb_rkqk(r1_curr, r2_curr, q1_curr, q2_curr, kr_curr, kq_curr, x, a1, b1, a2, b2);
    // update q2
    ind2 = setdiff(indn, ind1);
    x2 = x[ind2];
    if (kr_curr <= kq_curr) {
      aq2 = (n - kq_curr) * r2_curr + 1.0;
    }
    else {
      aq2 = (kr_curr - kq_curr) * r1_curr + (n - kr_curr) * r2_curr + 1.0;
    }
    q2_curr = rbeta(1, aq2, sum(x2) + 1.0)[0];
    lpost_curr = lpost_nb_rkqk(r1_curr, r2_curr, q1_curr, q2_curr, kr_curr, kq_curr, x, a1, b1, a2, b2);
    // update kr
    kr_prop = (int) round(rnorm(1, kr_curr + 0.0, s_kr)[0]);
    lpost_prop = lpost_nb_rkqk(r1_curr, r2_curr, q1_curr, q2_curr, kr_prop, kq_curr, x, a1, b1, a2, b2);
    log_alpha = lpost_prop - lpost_curr;
    if (log_u[2] < log_alpha) {
      kr_curr = kr_prop;
      lpost_curr = lpost_prop;
    }
    // update kq
    kq_prop = (int) round(rnorm(1, kq_curr + 0.0, s_kq)[0]);
    lpost_prop = lpost_nb_rkqk(r1_curr, r2_curr, q1_curr, q2_curr, kr_curr, kq_prop, x, a1, b1, a2, b2);
    log_alpha = lpost_prop - lpost_curr;
    if (log_u[3] < log_alpha) {
      kq_curr = kq_prop;
      lpost_curr = lpost_prop;
    }
    // save for output
    if (i >= burnin && (i - burnin + 1) % thin == 0) {
      j = (i - burnin + 1) / thin - 1;
      Rcout << j + 1 << endl;
      par_mat(j, _) = NumericVector::create(r1_curr, r2_curr, q1_curr, q2_curr, kr_curr, kq_curr);
      lpost_vec[j] = lpost_curr;
    }
  }
  // return
  List L;
  L["par"] = par_mat;
  L["lpost"] = lpost_vec;
  return L;
}

// [[Rcpp::export]]
const double llik_nb_rqk12(const double r1, const double r2, const double r3, const double q1, const double q2, const double q3, const int k1, const int k2, const NumericVector x) {
  // log-likelihood of NB model w/ 2 chgpts, both in r & q simultaneously
  const int n = x.size();
  double llik;
  NumericVector r = NumericVector::create(r1, r2, r3), q = NumericVector::create(q1, q2, q3);
  IntegerVector k = IntegerVector::create(k1, k2);
  if (is_true(any(r <= 0.0)) || is_true(any(q <= 0.0)) || is_true(any(q > 1.0)) || is_true(any(k < 1)) || is_true(any(k > n)) || k1 >= k2) {
    llik = -INFINITY;
  }
  else if (k2 == n) {
    // 1 chgpt
    llik = llik_nb_rqk(r1, r2, q1, q2, k1, x);
  }
  else {
    // k2 less than n, k1 less than k2
    IntegerVector ind1 = seq_len(k1) - 1,
      ind2 = seq_len(k2) - 1,
      ind3 = seq_len(n) - 1,
      ind4 = setdiff(ind2, ind1),
      ind5 = setdiff(ind3, ind2);
    NumericVector x1 = x[ind1], x2 = x[ind4], x3 = x[ind5];
    llik = sum(dnbinom(x1, r1, q1, true)) + sum(dnbinom(x2, r2, q2, true)) + sum(dnbinom(x3, r3, q3, true));
  }
  return llik;
}

// [[Rcpp::export]]
const double lpost_nb_rqk12(const double r1, const double r2, const double r3, const double q1, const double q2, const double q3, const int k1, const int k2, const NumericVector x, const double a1 = 0.01, const double b1 = 0.01, const double a2 = 0.01, const double b2 = 0.01, const double a3 = 0.01, const double b3 = 0.01) {
  // log-posterior of NB model w/ 2 chgpts, both in r & q simultaneously
  double lpost = llik_nb_rqk12(r1, r2, r3, q1, q2, q3, k1, k2, x) + dgamma(NumericVector::create(r1), a1, 1.0 / b1, true)[0] + dgamma(NumericVector::create(r2), a2, 1.0 / b2, true)[0] + dgamma(NumericVector::create(r3), a3, 1.0 / b3, true)[0];
  // add joint prior of k1 & k2?
  if (lpost != lpost) {
    lpost = -INFINITY;
  }
  return lpost;
}

// [[Rcpp::export]]
List rwm_nb_rqk12(const double r1, const double r2, const double r3, const double q1, const double q2, const double q3, const int k1, const int k2, const NumericVector x, const double s_r1, const double s_r2, const double s_r3, const int N = 1e+6, const int thin = 1, const int burnin = 0, const double a1 = 0.01, const double b1 = 0.01, const double a2 = 0.01, const double b2 = 0.01, const double a3 = 0.01, const double b3 = 0.01) {
  // componentwise RWM of NB model w/ 2 chgpts, both in r & q simultaneously
  double r1_curr = r1, r1_prop, r2_curr = r2, r2_prop, r3_curr = r3, r3_prop, q1_curr = q1, q1_prop, q2_curr = q2, q2_prop, q3_curr = q3, q3_prop, k1_curr = k1, k1_prop, k2_curr = k2, k2_prop;
  double lpost_curr = lpost_nb_rqk12(r1_curr, r2_curr, r3_curr, q1_curr, q2_curr, q3_curr, k1_curr, k2_curr, x, a1, b1, a2, b2, a3, b3), lpost_prop, log_alpha, aq1, aq2, aq3;
  NumericMatrix par_mat(N, 8);
  NumericVector lpost_vec(N), log_u(6), x1, x2, x3;
  List L;
  L["asdf"] = 1.0;
  return L;
/*  // componentwise RWM of NB model w/ chgpt in r & q separately
  double r1_curr = r1, r1_prop, r2_curr = r2, r2_prop, q1_curr = q1, q2_curr = q2, kr_curr = kr, kr_prop, kq_curr = kq, kq_prop;
  double lpost_curr = lpost_nb_rkqk(r1, r2, q1, q2, kr, kq, x, a1, b1, a2, b2), lpost_prop, log_alpha, aq1, aq2;
  NumericMatrix par_mat(N, 6);
  NumericVector lpost_vec(N), log_u(4), x1, x2, x3;
  const int n = x.size();
  const IntegerVector indn = seq_len(n) - 1;
  IntegerVector ind1, ind2, ind3;
  int i, j;
  RNGScope scope;
  for (i = 0; i < N * thin + burnin; i++) {
    log_u = log(runif(4));
    // update r1
    r1_prop = exp(rnorm(1, log(r1_curr), s_r1)[0]);
    lpost_prop = lpost_nb_rkqk(r1_prop, r2_curr, q1_curr, q2_curr, kr_curr, kq_curr, x, a1, b1, a2, b2);
    log_alpha = lpost_prop - lpost_curr + log(r1_prop) - log(r1_curr);
    if (log_u[0] < log_alpha) {
      r1_curr = r1_prop;
      lpost_curr = lpost_prop;
    }
    // update r2
    r2_prop = exp(rnorm(1, log(r2_curr), s_r2)[0]);
    lpost_prop = lpost_nb_rkqk(r1_curr, r2_prop, q1_curr, q2_curr, kr_curr, kq_curr, x, a1, b1, a2, b2);
    log_alpha = lpost_prop - lpost_curr + log(r2_prop) - log(r2_curr);
    if (log_u[1] < log_alpha) {
      r2_curr = r2_prop;
      lpost_curr = lpost_prop;
    }
    // update q1
    ind1 = seq_len(kq_curr) - 1;
    x1 = x[ind1];
    if (kr_curr <= kq_curr) {
      aq1 = kr_curr * r1_curr + (kq_curr - kr_curr) * r2_curr + 1.0;
    }
    else {
      aq1 = kq_curr * r1_curr + 1.0;
    }
    q1_curr = rbeta(1, aq1, sum(x1) + 1.0)[0];
    lpost_curr = lpost_nb_rkqk(r1_curr, r2_curr, q1_curr, q2_curr, kr_curr, kq_curr, x, a1, b1, a2, b2);
    // update q2
    ind2 = setdiff(indn, ind1);
    x2 = x[ind2];
    if (kr_curr <= kq_curr) {
      aq2 = (n - kq_curr) * r2_curr + 1.0;
    }
    else {
      aq2 = (kr_curr - kq_curr) * r1_curr + (n - kr_curr) * r2_curr + 1.0;
    }
    q2_curr = rbeta(1, aq2, sum(x2) + 1.0)[0];
    lpost_curr = lpost_nb_rkqk(r1_curr, r2_curr, q1_curr, q2_curr, kr_curr, kq_curr, x, a1, b1, a2, b2);
    // update kr
    kr_prop = (int) round(rnorm(1, kr_curr + 0.0, s_kr)[0]);
    lpost_prop = lpost_nb_rkqk(r1_curr, r2_curr, q1_curr, q2_curr, kr_prop, kq_curr, x, a1, b1, a2, b2);
    log_alpha = lpost_prop - lpost_curr;
    if (log_u[2] < log_alpha) {
      kr_curr = kr_prop;
      lpost_curr = lpost_prop;
    }
    // update kq
    kq_prop = (int) round(rnorm(1, kq_curr + 0.0, s_kq)[0]);
    lpost_prop = lpost_nb_rkqk(r1_curr, r2_curr, q1_curr, q2_curr, kr_curr, kq_prop, x, a1, b1, a2, b2);
    log_alpha = lpost_prop - lpost_curr;
    if (log_u[3] < log_alpha) {
      kq_curr = kq_prop;
      lpost_curr = lpost_prop;
    }
    // save for output
    if (i >= burnin && (i - burnin + 1) % thin == 0) {
      j = (i - burnin + 1) / thin - 1;
      Rcout << j + 1 << endl;
      par_mat(j, _) = NumericVector::create(r1_curr, r2_curr, q1_curr, q2_curr, kr_curr, kq_curr);
      lpost_vec[j] = lpost_curr;
    }
  }
  // return
  List L;
  L["par"] = par_mat;
  L["lpost"] = lpost_vec;
  return L; */
}

// [[Rcpp::export]]
List gibbs_p_lamk(const double lam1, const double lam2, const int k, const NumericVector x, const int N = 1e+6, const int thin = 1, const int burnin = 0, const double a1 = 0.01, const double b1 = 0.01, const double a2 = 0.01, const double b2 = 0.01) {
  const int n = x.size();
  double lam1_curr = lam1, lam2_curr = lam2;
  int k_curr = k;
  NumericMatrix par_mat(N, 3);
  const IntegerVector indn = seq_len(n) - 1, seq_k = seq_len(n);
  IntegerVector ind1, ind2;
  NumericVector x1, x2, cx = cumsum(x), seq_unscaled(n), seq_scaled(n), exponent(n);
  int i, j;
  RNGScope scope;
  for (i = 0; i < N * thin + burnin; i++) {
    ind1 = seq_len(k_curr) - 1;
    ind2 = setdiff(indn, ind1);
    x1 = x[ind1];
    x2 = x[ind2];
    lam1_curr = rgamma(1, a1 + sum(x1), 1.0 / (b1 + (double) k_curr))[0];
    lam2_curr = rgamma(1, a2 + sum(x2), 1.0 / (b2 + (double) (n - k_curr)))[0];
    exponent = (lam2_curr - lam1_curr) * (NumericVector) seq_k + (log(lam1_curr) - log(lam2_curr)) * cx;
    seq_unscaled = exp(exponent - max(exponent));
    seq_unscaled = ifelse(seq_unscaled != seq_unscaled, 0.0, seq_unscaled); // underflow gives nan
    seq_scaled = seq_unscaled / sum(seq_unscaled); // the probabilities
    k_curr = Rcpp::RcppArmadillo::sample(seq_k, 1, false, seq_scaled)[0];
    if (i >= burnin && (i - burnin + 1) % thin == 0) {
      j = (i - burnin + 1) / thin - 1;
      par_mat(j, _) = NumericVector::create(lam1_curr, lam2_curr, (double) k_curr);
      //      if ((j + 1) % 100 == 0) {
      //        Rcout << j + 1 << endl;
      //      }
    }
  }
  List L;
  L["par"] = par_mat;
  return L;
}

