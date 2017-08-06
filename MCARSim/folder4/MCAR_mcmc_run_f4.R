#------------------------------------
# This script is used to run mcmc
# Last updated date: 7/27/2017
#------------------------------------
setwd("/Users/peterwang/Desktop/Research/missingdata/Project/code/MCARSim/folder4")

##running n_iter iterations
system.time(source("MCAR_mcmc_update_f4.R")) 

#evaluate the method of imputing initial value of missing X_{it}
maxdiff.ini.X
mindiff.ini.X

burnin=5000
##posterior mean
(posterior.mean.eta=apply(eta_keep[-(1:burnin),],2, mean))
(posterior.mean.M=apply(M_keep[-(1:burnin),],2, mean))
(posterior.mean.v=apply(v_keep[-(1:burnin),],2, mean))
(posterior.mean.beta=apply(beta_keep[-(1:burnin),],2, mean))
(posterior.mean.sgmr2=mean(sgmr2_keep[-(1:burnin)]))
(posterior.mean.sgm2=mean(sgm2_keep[-(1:burnin)]))
(posterior.mean.E=mean(E_keep[-(1:burnin)]))

##mean of imputed X
MI.mean.X=apply(X_keep[-(1:burnin),], 2, mean)

##difference with the true X
(diff=MI.mean.X-(SVXYR$X)[R_sim==0])
max(diff)
min(diff)

##traceplots after burn-in
traceplot(x=as.mcmc(eta_keep[-(1:burnin),2]), ylab="eta_2")
traceplot(x=as.mcmc(eta_keep[-(1:burnin),3]), ylab="eta_3")
traceplot(x=as.mcmc(eta_keep[-(1:burnin),4]), ylab="eta_4")
traceplot(x=as.mcmc(v_keep[-(1:burnin),1]), ylab="v_1")
traceplot(x=as.mcmc(v_keep[-(1:burnin),2]), ylab="v_2")
traceplot(x=as.mcmc(beta_keep[-(1:burnin),1]), ylab="beta_1")
traceplot(x=as.mcmc(beta_keep[-(1:burnin),2]), ylab="beta_2")
traceplot(x=as.mcmc(M_keep[-(1:burnin),2]), ylab="M_2")
traceplot(x=as.mcmc(M_keep[-(1:burnin),3]), ylab="M_3")
traceplot(x=as.mcmc(M_keep[-(1:burnin),4]), ylab="M_4")
traceplot(x=as.mcmc(sgmr2_keep[-(1:burnin)]), ylab="sgmr2")
traceplot(x=as.mcmc(sgm2_keep[-(1:burnin)]), ylab="sgm2")
traceplot(x=as.mcmc(E_keep[-(1:burnin)]), ylab="E")











