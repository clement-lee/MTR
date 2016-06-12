remove(list = ls())
source("extraction.R") # run the whole script
optim.ctrl <- list(fnscale = -1, reltol = 1e-10, maxit = 5000)



### fitting Poisson & negative Binomial model
mle.p <- df0.weeks$all %>% mean # sample mean = MLE of lambda, no numerical computing required
obj0.nb <- optim(c(3, 0.5), llik_nb, x = df0.weeks$all, control = optim.ctrl)
mle.nb <- obj0.nb$par



### counts of count i.e. no. of weeks w/ a particular number of incidents
df0.counts <- seq(0, df0.weeks$all %>% max) %>%
    data_frame(count = .) %>%
    ## i) all
    left_join(df0.weeks %>% count(all), by = c("count" = "all")) %>% 
    mutate(all = ifelse(is.na(n), 0L, n)) %>% 
    select(-n) %>% 
    ## ii) signal
    left_join(df0.weeks %>% count(signal), by = c("count" = "signal")) %>% 
    mutate(signal = ifelse(is.na(n), 0L, n)) %>%
    select(-n) %>% 
    ## iii) train
    left_join(df0.weeks %>% count(train), by = c("count" = "train")) %>%
    mutate(train = ifelse(is.na(n), 0L, n)) %>% 
    select(-n) %>%
    ## iv) all, estimated by Poisson & negative Binomial
    mutate(all.est.p = count %>% dpois(mle.p) * n0.weeks,
           all.est.nb = count %>% dnbinom(mle.nb[1], mle.nb[2]) * n0.weeks,
           summand.p = (all - all.est.p) ^ 2 / all.est.p, # for pearson chi-sq test
           summand.nb = (all - all.est.nb) ^ 2 / all.est.nb) # same as summand.p
## we put it here as it is NOT needed for vis via ggplot2



### Pearson chi-squared goodness-of-fit test
(df0.counts %>% 
    summarise(Poisson = sum(summand.p), 
              Negative.Binomial = sum(summand.nb),
              Chi.squared = qchisq(0.95, df0.counts %>% nrow - 1)))



### linear trend model
t0 <- seq_along(df0.weeks$all)
obj0.nb.rt <- optim(c(mle.nb, 1e-4), llik_nb_rt, x = df0.weeks$all, t = t0, control = optim.ctrl)
mle.nb.rt <- obj0.nb.rt$par
obj0.nb.qt <- optim(c(mle.nb, 1e-4), llik_nb_qt, x = df0.weeks$all, t = t0, control = optim.ctrl)
mle.nb.qt <- obj0.nb.qt$par



### likelihood ratio test
(data_frame(critical.value = qchisq(0.95, 1),
            test.stat.nb.rt = 2 * (obj0.nb.rt$value - obj0.nb$value),
            test.stat.nb.qt = 2 * (obj0.nb.qt$value - obj0.nb$value)))
### interestingly, theta is "significant" in both models



### new data_frame for vis
df1.weeks <- df0.weeks %>% 
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
    obj0.rk <- optim(par0.rk, llik_nb_rk_fix_k, x = df0.weeks$all, t = t0, k = i, control = optim.ctrl)
    l0.rk[[i]] <- data_frame(
        llik = obj0.rk$value,
        r1 = obj0.rk$par[1],
        r2 = obj0.rk$par[2],
        q = obj0.rk$par[3]
    )
}
df0.rk <- l0.rk %>% 
    bind_rows %>% 
    bind_cols(df0.weeks, .) %>% 
    mutate(
        k = which.max(llik),
        llik = llik[k],
        r1 = r1[k],
        r2 = r2[k],
        q = q[k],
        r = ifelse(week <= week[k], r1, r2),
        mean = r * (1.0 - q) / q
    )
    


