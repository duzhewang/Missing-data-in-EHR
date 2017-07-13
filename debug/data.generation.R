#----------------------------------
# This file is used to 
# 1.generate the complete simulated data
# 2.set the prior
# Last updated date: 7/11/2017
#----------------------------------
rm(list=ls())
set.seed(0311)
library(MASS) ##use function mvrnorm
library(LaplacesDemon) ##use function rcat
library(MCMCpack) ##use function rinvgamma, default rate=1
library(mvtnorm)  ##use function dmvnorm, rmvt
library(pscl)  ##use function rigamma 
#library(metRology)  ##use function rt.scaled 
library(BayesLogit) ##use rpg
library(coda)
library(LearnBayes)

#---------------------------------
#      Global set-up       
#---------------------------------
n=100 ##number of subject
K=4 ##number of latent classes

##assign each subject the number of tracked quarters
T=sample(5:40, size=n, replace=TRUE) 

#--------------------------------------
#   Complete simulated data generation
#--------------------------------------
##set up true values of parameters
eta_sim=c(0, 0.7, 1.6, 0.9)
beta_sim=c(-0.4, 0.5)
M_sim=c(0, -0.6, 0.6, 1.2) 
v_sim=c(0.5, -0.3)
sgmr2_sim=1                ##true value of variance of b_{i}
sgm2_sim=1                 ##true value of variance of epsilon_{it}
E_sim=1                    ##true value of variance of e_{i}

##generate time-invariate covariate V_{i}
V_sim=rnorm(n=n, mean=0, sd=1)

##simulate latent class c[i] for each subject 
Pi_sim=matrix(0, n, K )   ##Pi[i,l]=P(c_{i}=l|V_{i})
for (i in 1:n){
  for (l in 1:K){
    Pi_sim[i,l]=exp(V_sim[i]*eta_sim[l])/sum(exp(V_sim[i]*eta_sim))
  }
}
c_sim=apply(Pi_sim, 1, function(x) rcat(n=1, p=x)) ##assign latent class to each subject

##generate random effect b_{i}
b_sim=rnorm(n,mean=0, sd=sqrt(sgmr2_sim))  

##set D_{i}, D^{*}_{i} and D^{**}_{i}
D=matrix(0, max(T), 2*n)              ##each two columns for one subject
for (i in 1:n){
  D[1:T[i], 2*i-1]=rep(1, T[i])
  D[1:T[i], 2*i]=rep(V_sim[i], T[i])
}

D_star=matrix(0, max(T), n)   
for (i in 1:n){
  D_star[1:T[i],i]=1:T[i]
}

D_dstar=matrix(0, max(T), n)
for (i in 1:n){
  D_dstar[1:T[i],i]=rep(1, T[i])
}

##sum_D is used in updating beta
sum_D=matrix(0, 2, 2)
for (i in 1:n){
  sum_D=sum_D+t(D[1:T[i],c(2*i-1, 2*i)])%*%D[1:T[i], c(2*i-1, 2*i)] 
}

##simulate X_{it}
X_sim=matrix(0, nrow=max(T), ncol = n)  ##each column is one subject
for (i in 1:n){
  X_sim[1:T[i],i]=D[1:T[i],c(2*i-1,2*i)]%*%beta_sim+
    D_star[1:T[i],i]*M_sim[c_sim[i]]+
    D_dstar[1:T[i],i]*b_sim[i]+          
    rnorm(n=T[i], mean = 0, sd=sqrt(sgm2_sim))
}


##generate random effect e_{i}
e_sim=rnorm(n, mean=0, sd=sqrt(E_sim))

##generate y_{it}
Y_sim=matrix(0,nrow=max(T), ncol=n)     ##each column is one subject
for (i in 1:n){
  top= exp(cbind(D[1:T[i],2*i], X_sim[1:T[i],i])%*%v_sim+rep(e_sim[i],T[i]))
  bot=rep(1, T[i])+top
  Y_sim[1:T[i],i]=rbinom(n=T[i], size=1, p=top/bot)   
}

#---------------------
#     PRIORS    
#---------------------
beta_pri=10^{4} 
M_pri=10^{4}
sgm2_pri=0.001
sgmr2_pri=0.001
v_pri=10^{4}
E_pri=0.001
eta_pri=10^{4}
