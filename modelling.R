remove(list = ls())
source("extraction.R") # run the whole script
library(Rcpp)
library(RcppArmadillo)
sourceCpp("mtr.cpp") # load required llik() functions etc.
optim.ctrl <- list(fnscale = -1, reltol = 1e-10, maxit = 5000)



### fitting Poisson & negative Binomial model
x <- df0.days$all
mle.p <- x %>% mean # sample mean = MLE of lambda, no numerical computing required
obj0.nb <- optim(c(3/7, 0.5), llik_nb, x = x, control = optim.ctrl)
mle.nb <- obj0.nb$par



### counts of count i.e. no. of days w/ a particular number of incidents
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



### Pearson chi-squared goodness-of-fit test
(df0.counts %>% 
    summarise(Poisson = sum(summand.p), 
              Negative.Binomial = sum(summand.nb),
              Chi.squared = qchisq(0.95, df0.counts %>% nrow - 1)))



### linear trend model
t0 <- seq_along(x)
obj0.nb.rt <- optim(c(mle.nb, 1e-5), llik_nb_rt, x = x, t = t0, control = optim.ctrl)
mle.nb.rt <- obj0.nb.rt$par
obj0.nb.qt <- optim(c(mle.nb, 1e-5), llik_nb_qt, x = x, t = t0, control = optim.ctrl)
mle.nb.qt <- obj0.nb.qt$par



### likelihood ratio test
(data_frame(critical.value = qchisq(0.95, 1),
            test.stat.nb.rt = 2 * (obj0.nb.rt$value - obj0.nb$value),
            test.stat.nb.qt = 2 * (obj0.nb.qt$value - obj0.nb$value)))
### interestingly, theta is "significant" in both models



### new data_frame for vis
df1.days <- df0.days %>% 
    mutate(t = t0, 
           rt = (mle.nb.rt[1] + mle.nb.rt[3] * t0),
           q = mle.nb.rt[2],
           mean.rt = rt * (1.0 - q)/ q, 
           r = mle.nb.qt[1],
           qt = (mle.nb.qt[2] + mle.nb.qt[3] * t0),
           mean.qt = r * (1.0 - qt) / qt)



### chgpt in r, likelihood method
l0.rk <- list()
for (i in seq_along(t0)) {
    print(i)
    par0.rk <- c(mle.nb[1], mle.nb[1], mle.nb[2])
    obj0.rk <- optim(par0.rk, llik_nb_rk_fix_k, x = x, k = i, control = optim.ctrl)
    l0.rk[[i]] <- data_frame(
        llik = obj0.rk$value,
        r1 = obj0.rk$par[1],
        r2 = obj0.rk$par[2],
        q = obj0.rk$par[3]
    )
}
df0.rk <- l0.rk %>% 
    bind_rows %>% 
    bind_cols(df0.days, .) %>% 
    mutate(
        k = which.max(llik),
        r1 = r1[k],
        r2 = r2[k],
        q = q[k],
        r = ifelse(op_date <= op_date[k], r1, r2),
        mean = r * (1.0 - q) / q
    )



### chgpt in r, Bayesian method
system.time(rwm0.rk <- rwm_nb_rk(2.5, 0.9, 0.2, 100, x, 0.375, 0.375, 200, 2e+4, 10, 5e+4)) # ~0.01s / 3 iterations (no thinning)



### chgpt in q, likelihood method
l0.qk <- list()
for (i in seq_along(t0)) {
    print(i)
    par0.qk <- c(mle.nb[1], mle.nb[2], mle.nb[2])
    obj0.qk <- optim(par0.qk, llik_nb_qk_fix_k, x = x, k = i, control = optim.ctrl)
    l0.qk[[i]] <- data_frame(
        llik = obj0.qk$value,
        r = obj0.qk$par[1],
        q1 = obj0.qk$par[2],
        q2 = obj0.qk$par[3]
    )
}
df0.qk <- l0.qk %>% 
    bind_rows %>% 
    bind_cols(df0.days, .) %>% 
    mutate(
        k = which.max(llik),
        r = r[k],
        q1 = q1[k],
        q2 = q2[k],
        q = ifelse(op_date <= op_date[k], q1, q2),
        mean = r * (1.0 - q) / q
    )



### chgpt in q, Bayesian method
system.time(rwm0.qk <- rwm_nb_qk(0.9, 0.2, 0.2, 100, x, 0.375, 200, 2e+4, 10, 5e+4)) # ~0.01s / 3 iterations (no thinning)



### chgpt in r & q simultaneously, likelihood method
l0.rqk <- list()
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



### chgpt in r & q simultaneously, Bayesian method
system.time(rwm0.rqk <- rwm_nb_rqk(2.5, 0.9, 0.2, 0.2, 100, x, 0.375, 0.375, 400, 2e+4, 10, 5e+4)) # ~0.004s / iteration (no thinning)

system.time(mwg0.lamk <- mwg_nb_lamk(2.5, 0.9, 0.15, 0.15, 100, x, 0.375, 0.375, 2e+4, 10, 5e+4)) # ~0.0004s / iteration (no thinning)


### pred. dist. for daily count
l0.pred.rqk <- sapply(0:6, dnbinom_rqk, rwm0.rqk$par[,1], rwm0.rqk$par[,2], rwm0.rqk$par[,3], rwm0.rqk$par[,4], rwm0.rqk$par[,5], n0.days, simplify = F) # daily
v0.pred.rqk <- sapply(l0.pred, mean) # combine with existing data_frames?
## have to think about how to work out weekly numbers



### chgpt in r & q separately, likelihood method
if (FALSE) {
    r1.rkqk <- r2.rkqk <- 
    q1.rkqk <- q2.rkqk <- 
    l0.rkqk <- matrix(NA, length(t0), length(t0))
    for (i in seq_along(t0)) {
        for (j in seq_along(t0)) {
            print(c(i, j))
            par0.rkqk <- c(mle.nb[1], mle.nb[1], mle.nb[2], mle.nb[2])
            obj0.rkqk <- optim(par0.rkqk, llik_nb_rkqk_fix_k, x = x, kr = i, kq = j, control = optim.ctrl)
            r1.rkqk[i, j] <- obj0.rkqk$par[1]
            r2.rkqk[i, j] <- obj0.rkqk$par[2]
            q1.rkqk[i, j] <- obj0.rkqk$par[3]
            q2.rkqk[i, j] <- obj0.rkqk$par[4]
            l0.rkqk[i, j] <- obj0.rkqk$value
        }
    }
}



### chgpt in r & q separately, Bayesian method
system.time(rwm0.rkqk <- rwm_nb_rkqk(2.5, 0.9, 0.2, 0.2, 500, 500, x, 0.375, 0.375, 300, 300, 2e+4, 100, 1e+5)) # ~0.06s per 11 iterations (no thinning)



### chgpt in lambda, Bayesian method
set.seed(123)
system.time(m0.lam <- gibbs_p_lamk(0.5, 0.5, 500, x, 1e+4, 100, 2e+4))



### 2 chgpts in r & q simultaneously, Bayesian method
set.seed(234)
system.time(rwm0.rqk12 <- rwm_nb_rqk12(1.46, 0.31, 3.32, 0.80, 0.67, 0.91, 690, 740, x, 0.5, 0.75, 0.5, 60, 60, 0.5, 5e+5, 1, 2e+4))


