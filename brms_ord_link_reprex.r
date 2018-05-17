# reprex
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
                  chains=1, iter=1000, seed=5, refresh=0, inits="0")

fit_probit_app <- brm(as.numeric(y) ~ z1+z2, data=dat10levs, family=cumulative("probit_approx"),
                      chains=1, iter=1000, seed=5, refresh=0, inits="0")

fit_cloglog <- brm(as.numeric(y) ~ z1+z2, data=dat10levs, family=cumulative("cloglog"),
                   chains=1, iter=1000, seed=5, refresh=0, inits="0")
