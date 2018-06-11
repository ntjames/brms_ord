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

# --> try different ordinal families
  
brm_fit_probit<-brm(as.numeric(y) ~ z1, data=dat100levs, family=cumulative(link="probit"),
         control=list(adapt_delta=0.95), init_r=0.00001,
         prior=c(prior(student_t(3,0,1),class="Intercept"),
                 prior(normal(0,1),class="b")) )

get_prior(as.numeric(y) ~ z1, data=dat100levs, family=cumulative(link="probit"))
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

# to get inverse link function from brms fit object
# fit$formula$family$linkinv()

## -- estimate conditional CDF -- ##

# check docs - how to get smooth SE

# helper function to create labels for newdata rows
lfunc<-function(x){
  # add err if x<=0
  if (x<=26){
    out<-LETTERS[1:x]
  } else {
    addtl<-x %/% 26
    av<-1:addtl
    moreLETTERS<-paste0(LETTERS,sort(rep(av,26)))
    nLETTERS<-c(LETTERS,moreLETTERS)
    out<-nLETTERS[1:x]
  }
  return(factor(out))
}

# --> just get posterior distribution rather than summary??
# how are %tile derived? what is relationship to Est.Error

getCDF<-function(brmfit,newdata,...){
  require(dplyr)
  
  #check that cumulative model used
  if(class(brmfit)!="brmsfit") stop("brmsfit object must be used")
  if(brm_fit_logit$family$family!="cumulative") stop("fit needs to be from cumulative() family")
  #check for brmfit and newdata
  # other checks?
  
  truey <- as.numeric( levels(brmfit$data$y) )
  # array (num MCMC draws x num rows in newdata x num params (int+vars))
  fitvals <- fitted(brmfit,newdata=newdata,summary=FALSE,...)
  
  #add group label to rows of newdata
  nd<-newdata %>% mutate(Var2=lfunc(nrow(newdata)))
  
  # Var1 is MCMC draw 
  # Var2 is row of newdata
  # Var3 is intercepts (alpha_1, alpha_2, etc.)
  fv_tab<-as.data.frame.table(fitvals)
  
  #within each covar group (Var2) and MCMC sample (Var1), get cdf
  fv_tab %>% group_by(Var2, Var1) %>% mutate(cdf=cumsum(Freq)) %>% 
    ungroup() %>%
    group_by(Var2, Var3) %>% 
    dplyr::summarize(mn_cdf=mean(cdf),
              med_cdf=median(cdf),
              cdf_q2.5=quantile(cdf,probs=0.025),
              cdf_q97.5=quantile(cdf,probs=0.975)) %>% 
    ungroup() %>% mutate(yval=rep(truey,nrow(newdata))) %>%
    full_join(., nd, by="Var2")
  
}

cond_cdf<-getCDF(brm_fit_logit,newdata=data.frame(z1=c(0,1)))

#! see ggalt stat=stepribbon or RcmdrPlugin.KMggplot2 geom_stepribbon for stepped ribbon
#! google: geom_ribbon geom_step
cond_cdf %>% ggplot(aes(group=z1))+
  geom_ribbon(aes(x=yval, ymin=cdf_q2.5,ymax=cdf_q97.5),fill="grey30", alpha=0.4)+
  geom_step(aes(x=yval,y=mn_cdf))+
  geom_line(data=n_cdf_0,aes(x=yval,y=cdf),color="blue",inherit.aes = FALSE) + 
  geom_line(data=n_cdf_1,aes(x=yval,y=cdf),color="blue",inherit.aes = FALSE) +
  xlab("") + ylab("Conditional CDF")

## -- estimate conditional mean -- ##

# just get posterior dist and post-process to get intervals, etc?
getMean<-function(brmfit,newdata,...){
  require(dplyr)
  
  truey <- as.numeric( levels(brmfit$data$y) )
  fitvals <- fitted(brmfit,newdata=newdata,summary=FALSE,...)
  
  #add group label to rows of newdata
  nd <- newdata %>% mutate(Var2=lfunc(nrow(newdata)))
  
  # Var1 is MCMC draw 
  # Var2 is row of newdata
  # Var3 is intercepts (alpha_1, alpha_2, etc.)
  fv_tab <- as.data.frame.table(fitvals)
  
  n_samp<-nsamples(brmfit)
  
  fv_tab %>% arrange(Var2, Var1) %>% 
    mutate(y=rep(truey,n_samp*nrow(nd)),fy_Py=Freq*y) %>% 
    group_by(Var2, Var1) %>% 
    dplyr::summarize(mn=sum(fy_Py)) %>%
    dplyr::summarize(mean_mn=mean(mn),
                     med_mn=median(mn),
                     sd_mn=sd(mn),
                     mn_q2.5=quantile(mn,probs=0.025),
                     mn_q97.5=quantile(mn,probs=0.975)) %>%
    full_join(., nd, by="Var2")

}

cond_mean<-getMean(brm_fit_logit,newdata=data.frame(z1=c(0,1)) )
cond_mean

## -- estimate conditional quantiles -- ##

# head(posterior_samples(brm_fit_logit))

# --> add multiple probs?? - no, just get function for 1 and vectorize??
# add fix for low quantiles (e.g. 0.0001)

getQuantile<-function(brmfit,newdata,q,...){
  require(dplyr)

  truey <- as.numeric( levels(brmfit$data$y) )
  fitvals <- fitted(brmfit,newdata=newdata,summary=FALSE,...)
  
  #add group label to rows of newdata
  nd<-newdata %>% mutate(Var2=lfunc(nrow(newdata)))
  
  # Var1 is MCMC draw 
  # Var2 is row of newdata
  # Var3 is intercepts (alpha_1, alpha_2, etc.)
  fv_tab<-as.data.frame.table(fitvals) 
  
  out<-fv_tab %>% group_by(Var2, Var1) %>% 
    mutate(cdf=cumsum(Freq),
           idx.1 = max(which(cdf<=q)),
           idx.2 = min(which(cdf>=q)),
           cdf.1 = cdf[idx.1],
           cdf.2 = cdf[idx.2]) %>%
    filter(Var3=="A") %>%
    mutate(idx.y1.cdf=ifelse(idx.1==-Inf,0,cdf.1),
           idx.y2.cdf=ifelse(idx.2==-Inf,0,cdf.2), 
           idx.y1=ifelse(idx.1==-Inf,-Inf,truey[idx.1]),
           idx.y2=ifelse(idx.2==-Inf,-Inf,truey[idx.2]),
           qtile=ifelse(idx.1==idx.2,idx.y1,
                        (idx.y2-idx.y1)/(idx.y2.cdf - idx.y1.cdf)*(q-idx.y1.cdf) + idx.y1))  %>%
    group_by(Var2) %>% 
    dplyr::summarize(mean_qtile=mean(qtile),
                     med_qtile=median(qtile),
                     sd_qtile=sd(qtile),
                     qtile_q2.5=quantile(qtile,probs=0.025),
                     qtile_q97.5=quantile(qtile,probs=0.975)) %>%
    full_join(., nd, by="Var2")
  
  return(out)
}

cond_med<-getQuantile(brm_fit_logit,newdata=data.frame(z1=c(0,1)),q=0.5)
cond_med


if (0){
group_by(Var2) %>% 
  dplyr::summarize(mn=sum(fy_Py)) %>%
  dplyr::summarize(mean_mn=mean(mn),
                   med_mn=median(mn),
                   sd_mn=sd(mn),
                   mn_q2.5=quantile(mn,probs=0.025),
                   mn_q97.5=quantile(mn,probs=0.975)) %>%
  full_join(., nd, by="Var2")

brmfit<-brm_fit_logit
newdata0=data.frame(z1=0)
newdata=data.frame(z1=c(0,1))
truey <- as.numeric( levels(brmfit$data$y) )
#truey <- c(-Inf, truey)

fitvals <- fitted(brmfit,newdata=newdata,summary=FALSE)
prob<-0.5

# Var1 is MCMC draw, 
# Var2 is same for all (can be dropped),
# Var3 is intercepts (alpha_1, alpha_2, etc.)
fv_tab<-as.data.frame.table(fitvals) %>% select(-Var2)

foo<-fv_tab %>% group_by(Var1) %>% mutate(cdf=cumsum(Freq),intercept=1:n())

# code for two cdf draws
#cdf_foo<-foo %>% filter(Var1=="A"|Var1=="B") 
#cdf_foo<-foo %>% ungroup() %>% filter(Var1=="A"|Var1=="B") 
#cdf_foo<-c(0,cdf_foo) #??

# idx.1 will be -inf if less than smallest observed cdf 
# is this case idx.y1 should be -Inf as well
# not sure ifelse is needed for idx.y2

# --> expand for multiple rows of newdata
tmp<-foo %>% mutate(idx.1=max(which(cdf<=prob)),
                        idx.2=min(which(cdf>=prob)),
                        idx.y1=ifelse(idx.1==-Inf,-Inf,truey[idx.1]),
                        idx.y2=ifelse(idx.2==-Inf,-Inf,truey[idx.2]),
                        idx.y1.cdf=ifelse(idx.1==-Inf,0,cdf[idx.1]),
                        idx.y2.cdf=ifelse(idx.2==-Inf,0,cdf[idx.2]),
                        quantile=ifelse(idx.1==idx.2,idx.y1,
                                        (idx.y2-idx.y1)/(idx.y2.cdf - idx.y1.cdf)*(prob-idx.y1.cdf) + idx.y1)) %>% 
  filter(intercept==1) %>% select(Var1, quantile)

tmp %>% ggplot(aes(x=quantile)) + geom_histogram()

#idx.1<-max(which(cdf_foo<=probs[1]))
#idx.2<-min(which(cdf_foo>=probs[1]))

#not sure why ifelse needed. 
# when is index>length(truey)?? maybe with duplicate ys??
idx.y1 <- ifelse(idx.1>length(truey), Inf, truey[idx.1])
idx.y2 <- ifelse(idx.2>length(truey), Inf, truey[idx.2])

#when is idx=0, err when prob is less than minimum?
# make sure these vals are correct
idx.y1.cdf <- cdf_foo[idx.1]
idx.y2.cdf <- cdf_foo[idx.2]


# old
#idx.y1.cdf <- ifelse(idx.1==0, 0, cdf_foo[cbind(1:dim(newdata)[1], idx.1)])
#idx.y2.cdf <- ifelse(idx.2>length(truey), 1, cdf_foo[cbind(1:dim(newdata)[1], idx.2)])

  
#index.y1.cdf <- ifelse(index.1==0, 0, m.cdf[cbind(1:dim(new.data)[1], index.1)])
#index.y2.cdf <- ifelse(index.2>length(order.y), 1, m.cdf[cbind(1:dim(new.data)[1], index.2)])

quantil <- ifelse(idx.1==idx.2, idx.y1, 
                       (idx.y2-idx.y1)/(idx.y2.cdf - idx.y1.cdf)*(probs[1]-idx.y1.cdf) + idx.y1) 
quantil <- ifelse(is.infinite(quantil), max(truey), quantil)
}
