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
system.time(mwg0.lamk <- mwg_nb_lamk(2.5, 0.9, 0.2, 0.2, 100, x, 4, 4, 5e+4, 100, 5e+4)) # ~0.0004s / iteration (no thinning)
# ~1s / 808 iterations on Fujitsu (no thinning)

system.time(mwg1.lamk <- mwg_nb_lamk_block(2.5, 0.9, 0.2, 0.2, 100, x, 0.1, 0.1, 0.1, 0.1, 0.8, 0.8, 2e+5, 1, 0)) # ~0.00123s / iteration (no thinning) on Fujitsu



### pred. dist. for daily count
l0.pred.rqk <- sapply(0:6, dnbinom_rqk, rwm0.rqk$par[,1], rwm0.rqk$par[,2], rwm0.rqk$par[,3], rwm0.rqk$par[,4], rwm0.rqk$par[,5], n0.days, simplify = F) # daily
v0.pred.rqk <- sapply(l0.pred, mean) # combine with existing data_frames?
## have to think about how to work out weekly numbers




