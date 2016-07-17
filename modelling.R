### Models:
### 0) Simple Poisson, par: lambda
### 1) Simple negative Binomial (NB), par: (r, q)
### 2) NB w/ linear trend in r, par: (r_i = r_0 + theta * t_i, q)
### 2.5) NB w/ linear trend in q, par: (r, q_i = q_0 + theta * t_i)
### 3) NB w/ 1 chgpt in r & q simul., par: (r1, r2, q1, q2, k)

remove(list = ls())
source("extraction.R") # run the whole script
library(Rcpp)
library(RcppArmadillo)
sourceCpp("mtr.cpp") # load required llik() functions etc.
optim.ctrl <- list(fnscale = -1, reltol = 1e-10, maxit = 5000)
x <- df0.days$all
t0 <- seq_along(x) # time indices i.e. support of chgpt k



### LIKELIHOOD: FIT models 0, 1, 2, 2.5 & 3
mle.p <- x %>% mean # sample mean = MLE of lambda, no numerical computing required
obj0.nb <- optim(c(3/7, 0.5), llik_nb, x = x, control = optim.ctrl)
mle.nb <- obj0.nb$par
obj0.rt <- optim(c(mle.nb, 1e-5), llik_nb_rt, x = x, t = t0, control = optim.ctrl)
mle.rt <- obj0.rt$par
obj0.qt <- optim(c(mle.nb, 1e-5), llik_nb_qt, x = x, t = t0, control = optim.ctrl)
mle.qt <- obj0.qt$par
system.time({
    l0.rqk <- sapply(t0, function(k) {
        print(k)
        par0.rqk <- c(mle.nb[1], mle.nb[1], mle.nb[2], mle.nb[2])
        obj0.rqk <- optim(par0.rqk, llik_nb_rqk_fix_k, x = x, k = k, control = optim.ctrl)
        return(obj0.rqk)
    }, simplify = F)
}) # ~780s on Fujitsu
ind.rqk <- l0.rqk %>% sapply(extract2, "value") %>% which.max
obj0.rqk <- l0.rqk[[ind.rqk]]
mle.rqk <- c(l0.rqk[[ind.rqk]]$par, t0[ind.rqk])



### LIKELIHOOD: COMPARE estimates w/ data
## counts of count i.e. no. of days w/ a particular # of incidents
df0.counts <- seq(0, x %>% max) %>%
    data_frame(count = .) %>%
    ## i) all
    left_join(df0.days %>% count(all), by = c("count" = "all")) %>% 
    mutate(all = ifelse(is.na(n), 0L, n)) %>% 
    select(-n) %>% 
    ## ii) signal
    left_join(df0.days %>% count(signal), by = c("count" = "signal")) %>% 
    mutate(signal = ifelse(is.na(n), 0L, n)) %>%
    select(-n) %>% 
    ## iii) train
    left_join(df0.days %>% count(train), by = c("count" = "train")) %>%
    mutate(train = ifelse(is.na(n), 0L, n)) %>% 
    select(-n) %>%
    ## iv) all, estimated by Poisson & negative Binomial
    mutate(all.est.p = count %>% dpois(mle.p) * n0.days,
           all.est.nb = count %>% dnbinom(mle.nb[1], mle.nb[2]) * n0.days)



### LIKELIHOOD: DIAGNOSTICS & SELECTION
## Pearson chi-squared goodness-of-fit test
df0.counts %>% 
    summarise(Poisson = sum((all - all.est.p) ^ 2 / all.est.p), 
              Negative.Binomial = sum((all - all.est.nb) ^ 2 / all.est.nb),
              Chi.squared = qchisq(0.95, df0.counts %>% nrow - 1)) %>% 
    print
## likelihood ratio tests & AIC
data_frame(cv.rt = qchisq(0.95, mle.rt %>% length - mle.nb %>% length),
           stat.rt = 2 * (obj0.rt$value - obj0.nb$value),
           AIC.rt = 2 * (mle.rt %>% length - obj0.rt$value),
           cv.qt = qchisq(0.95, mle.qt %>% length - mle.nb %>% length),
           stat.qt = 2 * (obj0.qt$value - obj0.nb$value),
           AIC.qt = 2 * (mle.qt %>% length - obj0.qt$value),
           cv.rqk = qchisq(0.95, mle.rqk %>% length - mle.nb %>% length),
           stat.rqk = 2 * (obj0.rqk$value - obj0.nb$value),
           AIC.rqk = 2 * (mle.rqk %>% length - obj0.rqk$value)) %>%
    print
### interestingly, all models are "significant"



### LIKELIHOOD: PREPARE new data_frame for vis
df1.days <- df0.days %>% 
    mutate(r.rt = (mle.rt[1] + mle.rt[3] * t0),
           q.rt = mle.rt[2],
           mean.rt = r.rt * (1.0 - q.rt) / q.rt, 
           r.qt = mle.qt[1],
           q.qt = (mle.qt[2] + mle.qt[3] * t0),
           mean.qt = r.qt * (1.0 - q.qt) / q.qt,
           r.rqk = ifelse(t0 <= mle.rqk[5], mle.rqk[1], mle.rqk[2]),
           q.rqk = ifelse(t0 <= mle.rqk[5], mle.rqk[3], mle.rqk[4]),
           mean.rqk = r.rqk * (1.0 - q.rqk) / q.rqk)
write_csv(df1.days, "ts_days.csv")



### BAYESIAN: FIT models 1, 2 & 3
## 1)
set.seed(1000)
system.time(mwg0.nb <- mwg_nb(2.5, 0.2, x, 0.2, 1e+4, 100, 5e+4))
## ~0.00118676s / iteration on avignon; ~0.0014014s / iteration on Fujitsu
write_csv(mwg0.nb, "mwg_nb.csv")
## 2)
set.seed(2000)
system.time(mwg0.rt <- mwg_nb_rt(2.5, 0.2, 0.0, x, t0, 0.25, 3e-4, 1e+4, 100, 5e+4))
## ~0.00270456s / iteration on avignon
write_csv(mwg0.rt, "mwg_rt.csv")
## 3)
set.seed(12345)
system.time(mwg0.rqk <- mwg_nb_rqk(2.5, 0.9, 0.2, 0.2, 100, x, 0.4, 0.35, 1e+4, 100, 5e+4))
## ~0.0013112s / iteration (no thinning) on XPS; ~0.00176656s / iteration on avignon
write_csv(mwg0.rqk, "mwg_rqk.csv")



### BAYESIAN: SELECTION via Bayes factor
offset <- mwg0.rqk$llik %>% max
## compare model 3 w/ model 1
K31 <- exp(mwg0.rqk$llik - offset) %>% mean / exp(mwg0.nb$llik - offset) %>% mean # ~100
print(c(K31, 2 * log(K31))) # strength of evidence: strong
## compare model 3 w/ model 2
K32 <- exp(mwg0.rqk$llik - offset) %>% mean / exp(mwg0.rt$llik - offset) %>% mean # ~6
print(c(K32, 2 * log(K32))) # strength of evidence: positive



### BAYESIAN: COMPARE estimates w/ data
## counts of count i.e. no of days w/ a particular # of incidents
mat0.rt <- mwg0.rt$r0 %>% # rt = r0 + theta * t
    matrix(., length(.), length(t0)) + 
    outer(mwg0.rt$theta, t0, FUN = "*")
mat0.q <- mwg0.rt$q %>% 
    matrix(., length(.), length(t0))
df1.counts <- seq(0, x %>% max) %>%
    data_frame(count = .) %>%
    ## i) data
    left_join(df0.days %>% count(all), by = c("count" = "all")) %>% 
    mutate(all = ifelse(is.na(n), 0L, n)) %>% 
    select(-n) %>%
    ## ii) estimated by model 1
    mutate(all.est.nb =
           count %>%
           sapply(dnbinom, mwg0.nb$r, mwg0.nb$q, simplify = F) %>%
           sapply(mean, simplify = F) %>%
           sapply("*", n0.days)
           ) %>% 
    ## iii) estimated by model 2
    mutate(all.est.rt = 
           count %>% 
           sapply(dnbinom, mat0.rt, mat0.q, simplify = F) %>% 
           sapply(mean, simplify = F) %>% 
           sapply("*", n0.days)
           ) %>% 
    ## iv) estimated by model 3
    mutate(all.est.rqk =
           count %>%
           sapply(dnbinom_rqk, mwg0.rqk$r1, mwg0.rqk$r2, mwg0.rqk$q1, mwg0.rqk$q2, mwg0.rqk$k, n0.days, simplify = F) %>%
           sapply(mean, simplify = F) %>%
           sapply("*", n0.days)
           ) %>% 
    print
write_csv(df1.counts, "counts_days.csv")

df1.freq <- df1.counts %>% 
    transmute(count, 
              all.est.nb = all.est.nb / n0.days,
              all.est.rt = all.est.rt / n0.days,
              all.est.rqk = all.est.rqk / n0.days,
              all.est.rqk.after = 
              count %>% 
              sapply(dnbinom, mwg0.rqk$r2, mwg0.rqk$q2, simplify = F) %>% 
              sapply(mean)
              ) %>% 
    print


