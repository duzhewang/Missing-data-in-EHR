#--------------------------------------------------------------------------------
# This script is used to generate the simulated data with 20% MCAR missing values
# Last updated date: 7/27/2017
#--------------------------------------------------------------------------------
rm(list=ls())
set.seed(0311)
library(MASS) ##use function mvrnorm
library(LaplacesDemon) ##use function rcat
library(MCMCpack) ##use function rinvgamma, default rate=1
library(mvtnorm)  ##use function dmvnorm, rmvt
library(pscl)  ##use function rigamma 
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
eta_sim=c(0, 0.5, 1.5, 1)
beta_sim=c(-0.4, 0.5)
M_sim=c(0, -0.6, 0.6, 1.2) 
v_sim=c(0.5, -0.3)
sgmr2_sim=1                ##true value of variance of b_{i}
sgm2_sim=1                 ##true value of variance of epsilon_{it}
E_sim=1                    ##true value of variance of e_{i}

##generate time-invariate covariate V_{i}
V_sim=rnorm(n=n, mean=0, sd=1)

##simulate latent class c[i] for each subject 
Pi_sim=exp(V_sim%o%eta_sim)
c_sim=apply(Pi_sim, 1, function(x) rcat(n=1, p=x/sum(x)))

##generate random effect b_{i}
b_sim=rnorm(n,mean=0, sd=sqrt(sgmr2_sim))  

##set D_{i}, D^{*}_{i} and D^{**}_{i}
D=matrix(0, sum(T), 2)
D[,1]=rep(1,sum(T))
D[,2]=rep(V_sim,T)

D_star=vector(length=sum(T))
for(i in 1:n){
  D_star[(sum(T[1:i])-T[i]+1):sum(T[1:i])]=1:T[i]
}

D_dstar=rep(1, sum(T))

##simulate X_{it}
BigM=rep(M_sim[c_sim],T)
Bigb=rep(b_sim, T)
X_sim=D%*%beta_sim+D_star*BigM+D_dstar*Bigb+rnorm(n=sum(T), mean=0, sd=sqrt(sgm2_sim))

##generate random effect e_{i}
e_sim=rnorm(n, mean=0, sd=sqrt(E_sim))

##generate y_{it}
BigVX=matrix(0, sum(T), 2)
BigVX[, 1]=rep(V_sim,T)
BigVX[, 2]=X_sim
top=exp(BigVX%*%v_sim+rep(e_sim,T))
bot=1+top
Y_sim=rbinom(n=sum(T), size=1, p=top/bot)

#--------------------------------
# SET MISSING VALUES: 20% MCAR    
#--------------------------------
mr=0.2 ##missing rate
R_sim=rbinom(sum(T), 1, 1-mr)  ##R=0 represents missing
#create a data frame
SVXYR=data.frame(Subject=rep(1:n, T), V=rep(V_sim,T), X=X_sim, Y=Y_sim, R=R_sim) 





