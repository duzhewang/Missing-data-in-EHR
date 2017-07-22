#----------------------------------------------------
# lmm_mcmc_run.R is used to do diagnostics
# Last updated date: 7/18/2017
#----------------------------------------------------
setwd("/Users/peterwang/Desktop/Research/missingdata/Project/code/LinearMixed")

##run n_iter iterations
system.time(source("lmm_mcmc_update.R")) 

burnin=5000

##posterior mean
(posterior.mean.beta=apply(beta_keep[-(1:burnin),],2, mean))
(posterior.mean.sgmr2=mean(sgmr2_keep[-(1:burnin)]))
(posterior.mean.sgm2=mean(sgm2_keep[-(1:burnin)]))

traceplot(x=as.mcmc(beta_keep[-(1:burnin),1]), ylab="beta_1")
traceplot(x=as.mcmc(beta_keep[-(1:burnin),2]), ylab="beta_2")
traceplot(x=as.mcmc(sgmr2_keep[-(1:burnin)]), ylab="sgmr2")
traceplot(x=as.mcmc(sgm2_keep[-(1:burnin)]), ylab="sgm2")







