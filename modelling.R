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




