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

if(0){ #scratch for conditional CDF
truey <- as.numeric( levels(fit$data$y) )
n_ints <- length(unique(truey))-1
  
# mean of draws
coef_z0 <- fixef(fit,summary=TRUE)[1:n_ints,1]
cdf_z0 <- expit(coef_z0)

coef_z1 <- fixef(fit,summary=TRUE)["z1",1]
cdf_z1 <- expit(coef_z0-coef_z1)

cdfdat <- data.frame(cbind(yval=truey,cdf_0=c(cdf_z0,1),cdf_1=c(cdf_z1,1) ))
  
#median of draws
coef_z0_m<-fixef(fit,summary=TRUE,robust=TRUE)[1:n_ints,1]
cdf_z0_m <- expit(coef_z0_m)
coef_z1_m<-fixef(fit,summary=TRUE,robust=TRUE)["z1",1]
cdf_z1_m<-expit(coef_z0_m-coef_z1_m)
cdfdat_m<-data.frame(cbind(yval=truey,cdf_0=c(cdf_z0_m,1),cdf_1=c(cdf_z1_m,1) ))

#plot
ggplot(cdfdat)+geom_step(aes(x=yval,y=cdf_0)) + 
  geom_step(aes(x=yval,y=cdf_1))+
  geom_line(data=n_cdf_0,aes(x=yval,y=cdf),color="blue") + 
  geom_line(data=n_cdf_1,aes(x=yval,y=cdf),color="blue")

# - add estimate of SE - #

coef_z0_all <- t(fixef(fit,summary=FALSE))[1:n_ints,]
cdf_z0_all <- expit(coef_z0_all)

#note mean from summary=TRUE same as rowMeans over all draws
cbind(coef_z0, rowMeans(coef_z0_all))

#median over all draws
cbind(coef_z0_m,apply(coef_z0_all,1,median))

#coef_z0_all[1:10,1:5]
#cdf_z0_all[1:10,1:5]

cdf_z0_bnds<-apply(cdf_z0_all,1,quantile,probs=c(0.025,0.975))

#fix this to use matrix operations
coef_z1_all<-t(fixef(fit,summary=FALSE))["z1",]

coef_z0_all[1:3,1:5]
matrix(coef_z1_all[1:5],nrow=3,ncol=5,byrow=TRUE)

predict(fit)

#dim(coef_z0_all)
#coef_z0_all[,1:5]-matrix(coef_z1_all[1:5],nrow=99,ncol=5,byrow=TRUE)

cdf_z1_all<-expit(coef_z0_all - matrix(coef_z1_all,nrow=99,ncol=4000,byrow=TRUE))
cdf_z1_bnds<-apply(cdf_z1_all,1,quantile,probs=c(0.025,0.975))

cdfdat_all<-data.frame(cbind(yval=truey,cdf_0=c(cdf_z0,1),
                             cdf_0_2.5=c(cdf_z0_bnds[1,],1),
                             cdf_0_97.5=c(cdf_z0_bnds[2,],1),
                             cdf_1=c(cdf_z1,1),
                             cdf_1_2.5=c(cdf_z1_bnds[1,],1),
                             cdf_1_97.5=c(cdf_z1_bnds[2,],1)))

ggplot(cdfdat_all)+
  geom_ribbon(aes(x=yval, ymin=cdf_0_2.5,ymax=cdf_0_97.5),fill="grey30",alpha=0.4)+
    geom_step(aes(x=yval,y=cdf_0)) +
  geom_ribbon(aes(x=yval, ymin=cdf_1_2.5,ymax=cdf_1_97.5),fill="grey30",alpha=0.4)+
    geom_step(aes(x=yval,y=cdf_1))+ 
  geom_line(data=n_cdf_0,aes(x=yval,y=cdf),color="blue") + 
  geom_line(data=n_cdf_1,aes(x=yval,y=cdf),color="blue")


fy_z0_brm_0<-fitted(fit,newdata=data.frame(z1=0))

library(tidyr)
library(dplyr)
fy_z0_brm_1<-as.data.frame.table(fy_z0_brm_0) %>% spread(Var2,Freq) %>% select(-Var1)

cumsum(fy_z0_brm_1$Estimate)

ff<-fitted(fit,newdata=data.frame(z1=0),summary=FALSE)

dim(ff)

#Var1 is draws, Var2 is same for all can be dropped, Var3 is intercepts (alpha_1, alpha_2, etc.)
fff<-as.data.frame.table(ff) %>% select(-Var2)

length(table(fff$Var1))

fff_dat<-fff %>% group_by(Var1) %>% 
  mutate(cdf=cumsum(Freq)) %>% ungroup() %>%
  group_by(Var3) %>% 
  summarize(mn=mean(cdf),
            qnt2.5=quantile(cdf,probs=0.025),
            qnt97.5=quantile(cdf,probs=0.975)) %>% 
  ungroup() %>% mutate(yval=truey)


fff %>% filter(Var3=="A") %>% summary(mn=mean(Freq))
fff %>% group_by(Var3) %>% summarize(mn=mean(Freq),sd=sd(Freq))

fy_z0_brm<-predict(fit,newdata=data.frame(z1=0))
fy_z1_brm<-predict(fit,newdata=data.frame(z1=1))

Fy_z0<-cumsum(fy_z0_brm)
Fy_z1<-cumsum(fy_z1_brm)

cdfdat_brms<-data.frame(cbind(yval=truey,cdf_0=Fy_z0,cdf_1=Fy_z1))

ggplot(cdfdat_brms)+
  geom_step(aes(x=yval,y=cdf_0))+
  geom_step(aes(x=yval,y=cdf_1))+ 
  geom_line(data=n_cdf_0,aes(x=yval,y=cdf),color="blue") + 
  geom_line(data=n_cdf_1,aes(x=yval,y=cdf),color="blue")




}

#use brms functions

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

cdf %>% ggplot(aes(group=z1))+
  geom_ribbon(aes(x=yval, ymin=cdf_q2.5,ymax=cdf_q97.5),fill="grey30", alpha=0.4)+
  geom_step(aes(x=yval,y=mn_cdf))+
  geom_line(data=n_cdf_0,aes(x=yval,y=cdf),color="blue",inherit.aes = FALSE) + 
  geom_line(data=n_cdf_1,aes(x=yval,y=cdf),color="blue",inherit.aes = FALSE)


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


if (0){
#f(y_1)=F(y_1) and f(y_i)=F(y_i)-F(y_i-1) for i>1
fy_z0<-c(cdfdat_all$cdf_0[1],diff(cdfdat_all$cdf_0))

#mean

# Z1=0
sum(truey*fy_z0)

getmn<-function(i){
  cdfi<-c(cdf_z0_all[,i],1)
  fy<-c(cdfi[1],diff(cdfi))
  mn<-sum(truey*fy)
return(mn)
}

means_z0<-sapply(1:4000,getmn)
sd(means_z0) # sigma/sqrt(n)
1/sqrt(100)
qplot(means_z0)


# Z1=1
fy_z1<-c(cdfdat_all$cdf_1[1],diff(cdfdat_all$cdf_1))
sum(truey*fy_z1)


getmn2<-function(i){
  cdfi<-c(cdf_z1_all[,i],1)
  fy<-c(cdfi[1],diff(cdfi))
  mn<-sum(truey*fy)
  return(mn)
}

means_z1<-sapply(1:4000,getmn2)
sd(means_z1) # sigma/sqrt(n)
1/sqrt(100)
qplot(means_z1)


# with brms
fy_z0_brm<-predict(fit,newdata=data.frame(z1=0))
fy_z1_brm<-predict(fit,newdata=data.frame(z1=1))


sum(truey*fy_z0_brm)
sum(truey*fy_z1_brm)

}

## -- estimate conditional quantiles -- ##


getQuantile<-function(brmfit,newdata,q,...){
  
  
}