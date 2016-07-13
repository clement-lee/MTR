remove(list = ls())
source("extraction.R") # run the whole script
library(Rcpp)
library(RcppArmadillo)
sourceCpp("mtr.cpp") # load required llik() functions etc.
optim.ctrl <- list(fnscale = -1, reltol = 1e-10, maxit = 5000)



### LIKELIHOOD: fitting Poisson & negative Binomial model
x <- df0.days$all
mle.p <- x %>% mean # sample mean = MLE of lambda, no numerical computing required
obj0.nb <- optim(c(3/7, 0.5), llik_nb, x = x, control = optim.ctrl)
mle.nb <- obj0.nb$par



### LIKELIHOOD: counts of count i.e. no. of days w/ a particular # of incidents
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



### LIKELIHOOD: Pearson chi-squared goodness-of-fit test
(df0.counts %>% 
    summarise(Poisson = sum(summand.p), 
              Negative.Binomial = sum(summand.nb),
              Chi.squared = qchisq(0.95, df0.counts %>% nrow - 1)))



### LIKELIHOOD: linear trend models & likelihood ratio test
t0 <- seq_along(x)
obj0.nb.rt <- optim(c(mle.nb, 1e-5), llik_nb_rt, x = x, t = t0, control = optim.ctrl)
mle.nb.rt <- obj0.nb.rt$par
obj0.nb.qt <- optim(c(mle.nb, 1e-5), llik_nb_qt, x = x, t = t0, control = optim.ctrl)
mle.nb.qt <- obj0.nb.qt$par
data_frame(critical.value = qchisq(0.95, 1),
           test.stat.nb.rt = 2 * (obj0.nb.rt$value - obj0.nb$value),
           test.stat.nb.qt = 2 * (obj0.nb.qt$value - obj0.nb$value)) %>%
    print
### interestingly, theta is "significant" in both models



### LIKELIHOOD: new data_frame for vis
df1.days <- df0.days %>% 
    mutate(t = t0, 
           rt = (mle.nb.rt[1] + mle.nb.rt[3] * t0),
           q = mle.nb.rt[2],
           mean.rt = rt * (1.0 - q)/ q, 
           r = mle.nb.qt[1],
           qt = (mle.nb.qt[2] + mle.nb.qt[3] * t0),
           mean.qt = r * (1.0 - qt) / qt)



### LIKELIHOOD: chgpt in r & q simultaneously
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



### BAYESIAN: simple NB model (model 1)
set.seed(1000)
system.time(mwg0.nb <- mwg_nb(2.5, 0.2, x, 0.2, 1e+4, 10, 5e+4))
## ~0.00118676s / iteration on avignon; ~0.0014014s / iteration on Fujitsu
write_csv(mwg0.nb, "mwg_nb.csv")



### BAYESIAN: chgpt in r & q simultaneously (model 2)
set.seed(12345)
system.time(mwg0.rqk <- mwg_nb_rqk(2.5, 0.9, 0.2, 0.2, 100, x, 0.4, 0.35, 1e+4, 10, 5e+4))
## ~0.0013112s / iteration (no thinning) on XPS; ~0.00176656s / iteration on avignon
write_csv(mwg0.rqk, "mwg_rqk.csv")



### BAYESIAN: Bayes factor for model 2 compared to model 1
offset <- mwg0.rqk$llik %>% max
K12 <- exp(mwg0.rqk$llik - offset) %>% mean / exp(mwg0.nb$llik - offset) %>% mean
2 * log(K12)



### BAYESIAN: counts of count i.e. no of days w/ a particular # of incidents
df1.counts <- seq(0, x %>% max) %>%
    data_frame(count = .) %>%
    ## i) all
    left_join(df0.days %>% count(all), by = c("count" = "all")) %>% 
    mutate(all = ifelse(is.na(n), 0L, n)) %>% 
    select(-n) %>%
    ## ii) all, estimated by simple NB model
    mutate(all.est.nb =
           count %>%
           sapply(dnbinom, mwg0.nb[,1], mwg0.nb[,2], simplify = F) %>%
           sapply(mean) %>%
           sapply("*", n0.days)
           ) %>% 
    ## iii) all, estimated by 1-chgpt simul. NB model
    mutate(all.est.nb.rqk =
           count %>%
           sapply(dnbinom_rqk, mwg0.rqk[,1], mwg0.rqk[,2], mwg0.rqk[,3], mwg0.rqk[,4], mwg0.rqk[,5], n0.days, simplify = F) %>%
           sapply(mean) %>%
           sapply("*", n0.days)
           ) 


