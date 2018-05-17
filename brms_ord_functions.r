rm(list=ls())

library(brms)
library(rms)
library(dplyr)

# modified from Liu et al. sim code
generate.data.2 <- function(seed=1, n=200, p=0.5, alpha=0, beta=c(1, -0.5), sigma=1){
  set.seed(seed)
  z1 <- sample(c(0,1), size=n, replace=TRUE, prob=c(1-p, p))
  y <- rnorm(n, alpha+beta[1]*z1, sigma)
  data <- data.frame(y=y, z1=z1)
  return(data)
}

dat100levs <- generate.data.2(seed=2458,n=100,beta=c(2,0))
dat100levs$y <- factor(dat100levs$y)

#with(subset(dat100levs,z1==1),hist( as.numeric(levels(dat100levs$y)) ))
#with(subset(dat100levs,z1==0),hist( as.numeric(levels(dat100levs$y)) ))
#hist(dat100levs$y)

# call this once to distribute chains across cpu cores:
options(mc.cores=parallel::detectCores())

#possible links: "logit", "probit", "probit_approx", "cloglog", "cauchit"
brm_fit_logit<-brm(as.numeric(y) ~ z1, data=dat100levs, family=cumulative(link="logit"),
         control=list(adapt_delta=0.95))

if (0){
#! get error using probit, probit_approx, cloglog
#Rejecting initial value:
#  Log probability evaluates to log(0), i.e. negative infinity.
#  Stan can't start sampling from this initial value.
# --> check link_disc; link for discrimination in ordinal models
# --> check inits; probit, probit_approx, cloglog don't work with inits="0"

brm_fit_probit<-brm(as.numeric(y) ~ z1, data=dat100levs, family=cumulative(link="probit"),
         control=list(adapt_delta=0.95))

make_stancode(as.numeric(y) ~ z1, data=dat100levs, family=cumulative(link="probit"))
make_standata(as.numeric(y) ~ z1, data=dat100levs, family=cumulative(link="probit"))

brm_fit_probit_app<-brm(as.numeric(y) ~ z1, data=dat100levs, family=cumulative(link="probit_approx"),
         control=list(adapt_delta=0.95))

brm_fit_cloglog<-brm(as.numeric(y) ~ z1, data=dat100levs, family=cumulative(link="cloglog"),
         control=list(adapt_delta=0.95))

#! cauchit takes an extremely long time (> 1 hr), but works
# link_disc = "identity"
if(0){
fit<-brm(as.numeric(y) ~ z1, data=dat100levs, family=cumulative(link="cauchit"),
         control=list(adapt_delta=0.95))
fit
}

}

#true dist. for z1=0
norm_cdf_z10 <- curve(pnorm(x,0,1),-2.3,3.6)
n_cdf_0 <- data.frame(yval=norm_cdf_z10$x, cdf=norm_cdf_z10$y)

#true dist. for z1=1
norm_cdf_z11 <- curve(pnorm(x,2,1),-2.3,3.6)
n_cdf_1 <- data.frame(yval=norm_cdf_z11$x, cdf=norm_cdf_z11$y)

#inverse logit
#expit <- function(y) exp(y)/(1+exp(y))

# get inverse link function
fit$formula$family$linkinv()

## -- estimate conditional CDF -- ##

# check docs - how to get smooth SE

getCDF<-function(brmfit,newdata,...){
  require(dplyr)
  
  #check that cumulative model used
  #check for brmfit and newdata
  # other checks?
  
  truey <- as.numeric( levels(brmfit$data$y) )
  
  fitvals <- fitted(brmfit,newdata=newdata,summary=FALSE,...)
  
  # Var1 is MCMC draw, 
  # Var2 is same for all (can be dropped),
  # Var3 is intercepts (alpha_1, alpha_2, etc.)
  fv_tab<-as.data.frame.table(fitvals) %>% select(-Var2)
  
  fv_tab %>% group_by(Var1) %>% mutate(cdf=cumsum(Freq)) %>% 
    ungroup() %>%
    group_by(Var3) %>% 
    dplyr::summarize(mn_cdf=mean(cdf),
              med_cdf=median(cdf),
              cdf_q2.5=quantile(cdf,probs=0.025),
              cdf_q97.5=quantile(cdf,probs=0.975)) %>% 
    ungroup() %>% mutate(yval=truey) %>% 
    select(-Var3)
  
}

cdf_0<-getCDF(brm_fit_logit,newdata=data.frame(z1=0)) %>% mutate(z1=0)
cdf_1<-getCDF(brm_fit_logit,newdata=data.frame(z1=1)) %>% mutate(z1=1)

cdf <- rbind(cdf_0,cdf_1)

#! see ggalt stat=stepribbon or RcmdrPlugin.KMggplot2 geom_stepribbon for stepped ribbon
#! google: geom_ribbon geom_step

png("cdf1.png")
cdf %>% ggplot(aes(group=z1))+
  geom_ribbon(aes(x=yval, ymin=cdf_q2.5,ymax=cdf_q97.5),fill="grey30", alpha=0.4)+
  geom_step(aes(x=yval,y=mn_cdf))+
  geom_line(data=n_cdf_0,aes(x=yval,y=cdf),color="blue",inherit.aes = FALSE) + 
  geom_line(data=n_cdf_1,aes(x=yval,y=cdf),color="blue",inherit.aes = FALSE) +
  xlab("") + ylab("Conditional CDF")
dev.off()

# how are %tile derived? what is relationship to Est.Error

# could get SE of F_y = G^{-1}(alpha_i-beta X) using delta method


## -- estimate conditional mean -- ##

getMean<-function(brmfit,newdata,...){
  require(dplyr)
  
  truey <- as.numeric( levels(brmfit$data$y) )
  
  fitvals <- fitted(brmfit,newdata=newdata,summary=FALSE,...)
  
  # Var1 is MCMC draw, 
  # Var2 is same for all (can be dropped),
  # Var3 is intercepts (alpha_1, alpha_2, etc.)
  fv_tab <- as.data.frame.table(fitvals) %>% select(-Var2)
  
  n_samp<-nsamples(brmfit)
  
  fv_tab %>% arrange(Var1) %>% 
    mutate(y=rep(truey,n_samp),fy_Py=Freq*y) %>% group_by(Var1) %>% 
    dplyr::summarize(mn=sum(fy_Py))  %>% ungroup() %>%
    dplyr::summarize(mean_mn=mean(mn),
              med_mn=median(mn),
              sd_mn=sd(mn),
              mn_q2.5=quantile(mn,probs=0.025),
              mn_q97.5=quantile(mn,probs=0.975)) 
  
}

getMean(brm_fit_logit,newdata=data.frame(z1=0))
getMean(brm_fit_logit,newdata=data.frame(z1=1))


## -- estimate conditional quantiles -- ##


getQuantile<-function(brmfit,newdata,q,...){
  
  
}