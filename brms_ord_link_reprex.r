# reprex
rm(list=ls())
library(dplyr)
library(brms)

generate.data <- function(seed=1, n=200, p=0.5, alpha=0, beta=c(1, -0.5), sigma=1){
  set.seed(seed)
  z1 <- sample(c(0,1), size=n, replace=TRUE, prob=c(1-p, p))
  z2 <- rnorm(n, 0, 1)
  y <- rnorm(n, alpha+beta[1]*z1 + beta[2]*z2, sigma)
  data <- data.frame(y=y, z1=z1, z2=z2)
  return(data)
}

dat10levs <- generate.data(seed=45232,n=100) %>% mutate(y=cut_number(y,10))

fit_logit <- brm(as.numeric(y) ~ z1+z2, data=dat10levs, family=cumulative("logit"),
                 chains=1, iter=1000, seed=5, refresh=0) 

fit_probit <- brm(as.numeric(y) ~ z1+z2, data=dat10levs, family=cumulative("probit"),
                  chains=1, iter=1000, seed=5, refresh=0, init_r=0.001,
                  prior=set_prior("normal(0,5)",class="b"))

prior_summary(fit_probit)

# don't need probit_approx if probit is working
#fit_probit_app <- brm(as.numeric(y) ~ z1+z2, data=dat10levs, family=cumulative("probit_approx"),
#                      chains=1, iter=1000, seed=5, refresh=0)


#try manually specifying intercept prior to be more informative
fit_cloglog <- brm(as.numeric(y) ~ z1+z2, data=dat10levs, family=cumulative("cloglog"),
                   chains=1, iter=1000, seed=5, refresh=0, init_r=0.000001,
                   prior=set_prior("normal(0,5)",class="b"))

fit_cloglog_sratio <- brm(as.numeric(y) ~ z1+z2, data=dat10levs, family=sratio("cloglog"),
                   chains=1, iter=1000, seed=5, refresh=0)

fitted(fit_cloglog_sratio)


get_prior(as.numeric(y) ~ z1+z2, data=dat10levs, family=cumulative("cloglog"))
