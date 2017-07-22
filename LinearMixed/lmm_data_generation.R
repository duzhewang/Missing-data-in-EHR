#------------------------------------------------------------------------
# lmm_data_generation.R is used to 
# 1.generate the complete simulated data by linear mixed effects model
# 2.set the prior
# Last updated date: 7/18/2017
#-----------------------------------------------------------------------
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
#library(doParallel)

#---------------------------------
#      Global set-up       
#---------------------------------
n=100 ##number of subject

##assign each subject the number of tracked quarters
T=sample(5:40, size=n, replace=TRUE) 

#--------------------------------------
#   Complete simulated data generation
#--------------------------------------
##set up true values of parameters
beta_sim=c(-0.4, 0.5)
sgmr2_sim=sgmr2sim               ##true value of variance of b_{i}
sgm2_sim=sgm2sim                ##true value of variance of epsilon_{it}

##generate time-invariate covariate V_{i}
V_sim=rnorm(n=n, mean=0, sd=1)

##generate random effect b_{i}
b_sim=rnorm(n,mean=0, sd=sqrt(sgmr2_sim))  

##set D_{i} and D^{**}_{i}
D=matrix(0, sum(T), 2)
D[,1]=rep(1,sum(T))
D[,2]=rep(V_sim,T)

D_dstar=rep(1, sum(T))

##simulate X_{it}
Bigb=rep(b_sim, T)
X_sim=D%*%beta_sim+D_dstar*Bigb+rnorm(n=sum(T), mean=0, sd=sqrt(sgm2_sim))

#------------
# PRIORS    
#------------
beta_pri=10^{4} 
sgm2_pri=0.001
sgmr2_pri=0.001



