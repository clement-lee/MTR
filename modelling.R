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



### LIKELIHOOD: FIT models 0, 1, 2 & 2.5
mle.p <- x %>% mean # sample mean = MLE of lambda, no numerical computing required
obj0.nb <- optim(c(3/7, 0.5), llik_nb, x = x, control = optim.ctrl)
mle.nb <- obj0.nb$par
t0 <- seq_along(x)
obj0.rt <- optim(c(mle.nb, 1e-5), llik_nb_rt, x = x, t = t0, control = optim.ctrl)
mle.rt <- obj0.rt$par
obj0.qt <- optim(c(mle.nb, 1e-5), llik_nb_qt, x = x, t = t0, control = optim.ctrl)
mle.qt <- obj0.qt$par



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
           all.est.nb = count %>% dnbinom(mle.nb[1], mle.nb[2]) * n0.days,
           summand.p = (all - all.est.p) ^ 2 / all.est.p, # for pearson chi-sq test
           summand.nb = (all - all.est.nb) ^ 2 / all.est.nb) # same as summand.p
## we put it here as it is NOT needed for vis via ggplot2



### LIKELIHOOD: TESTS
## Pearson chi-squared goodness-of-fit test
(df0.counts %>% 
    summarise(Poisson = sum(summand.p), 
              Negative.Binomial = sum(summand.nb),
              Chi.squared = qchisq(0.95, df0.counts %>% nrow - 1)))
## likelihood ratio tests
data_frame(critical.value = qchisq(0.95, 1),
           test.stat.rt = 2 * (obj0.rt$value - obj0.nb$value),
           test.stat.qt = 2 * (obj0.qt$value - obj0.nb$value)) %>%
    print
### interestingly, theta is "significant" in both models



### LIKELIHOOD: PREPARE new data_frame for vis (OBSOLETE?)
df1.days <- df0.days %>% 
    mutate(t = t0, 
           rt = (mle.rt[1] + mle.rt[3] * t0),
           q = mle.rt[2],
           mean.rt = rt * (1.0 - q)/ q, 
           r = mle.qt[1],
           qt = (mle.qt[2] + mle.qt[3] * t0),
           mean.qt = r * (1.0 - qt) / qt)



### LIKELIHOOD: FIT model 3
l0.rqk <- list()
system.time({
    for (i in seq_along(t0)) {
        print(i)
        par0.rqk <- c(mle.nb[1], mle.nb[1], mle.nb[2], mle.nb[2])
        obj0.rqk <- optim(par0.rqk, llik_nb_rqk_fix_k, x = x, k = i, control = optim.ctrl)
        l0.rqk[[i]] <- data_frame(
            llik = obj0.rqk$value,
            r1 = obj0.rqk$par[1],
            r2 = obj0.rqk$par[2],
            q1 = obj0.rqk$par[3],
            q2 = obj0.rqk$par[4]
        )
    }
})
df0.rqk <- l0.rqk %>% 
    bind_rows %>% 
    bind_cols(df0.days, .) %>% 
    mutate(
        k = which.max(llik),
        r1 = r1[k],
        r2 = r2[k],
        q1 = q1[k],
        q2 = q2[k],
        r = ifelse(op_date <= op_date[k], r1, r2),
        q = ifelse(op_date <= op_date[k], q1, q2),
        mean = r * (1.0 - q) / q
    )
## maybe save df0.rqk as csv file for easy access?



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
print(c(K31, 2 * log(K31)))
## compare model 3 w/ model 2
K32 <- exp(mwg0.rqk$llik - offset) %>% mean / exp(mwg0.rt$llik - offset) %>% mean # ~6
print(c(K32, 2 * log(K32)))



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
           sapply(mean) %>%
           sapply("*", n0.days)
           ) %>% 
    print

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


