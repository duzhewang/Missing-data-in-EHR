#-----------------------------------------------------
# this script is used to run mcmc and do diagnostics
# Last updated date: 7/18/2017
#----------------------------------------------------
setwd("/Users/peterwang/Desktop/Research/missingdata/Project/code/CompleteDataSim/method2")

##run n_iter iterations
system.time(source("cpl_mcmc_update_m2.R")) 

burnin=5000

##posterior mean
(posterior.mean.eta=apply(eta_keep[-(1:burnin),],2, mean))
(posterior.mean.M=apply(M_keep[-(1:burnin),],2, mean))
(posterior.mean.v=apply(v_keep[-(1:burnin),],2, mean))
(posterior.mean.beta=apply(beta_keep[-(1:burnin),],2, mean))
(posterior.mean.sgmr2=mean(sgmr2_keep[-(1:burnin)]))
(posterior.mean.sgm2=mean(sgm2_keep[-(1:burnin)]))
(posterior.mean.E=mean(E_keep[-(1:burnin)]))

##traceplots after burn-in
path="/Users/peterwang/Desktop/Research/missingdata/Project/code/CompleteDataSim/method2/results"

jpeg(file=paste(path, 'eta2.png', sep="/"))
traceplot(x=as.mcmc(eta_keep[-(1:burnin),2]), ylab="eta_2")
dev.off()

jpeg(file=paste(path, 'eta3.png', sep="/"))
traceplot(x=as.mcmc(eta_keep[-(1:burnin),3]), ylab="eta_3")
dev.off()

jpeg(file=paste(path, 'eta4.png', sep="/"))
traceplot(x=as.mcmc(eta_keep[-(1:burnin),4]), ylab="eta_4")
dev.off()

jpeg(file=paste(path, 'v1.png', sep="/"))
traceplot(x=as.mcmc(v_keep[-(1:burnin),1]), ylab="v_1")
dev.off()

jpeg(file=paste(path, 'v2.png', sep="/"))
traceplot(x=as.mcmc(v_keep[-(1:burnin),2]), ylab="v_2")
dev.off()

jpeg(file=paste(path, 'beta1.png', sep="/"))
traceplot(x=as.mcmc(beta_keep[-(1:burnin),1]), ylab="beta_1")
dev.off()

jpeg(file=paste(path, 'beta2.png', sep="/"))
traceplot(x=as.mcmc(beta_keep[-(1:burnin),2]), ylab="beta_2")
dev.off()

jpeg(file=paste(path, 'M2.png', sep="/"))
traceplot(x=as.mcmc(M_keep[-(1:burnin),2]), ylab="M_2")
dev.off()

jpeg(file=paste(path, 'M3.png', sep="/"))
traceplot(x=as.mcmc(M_keep[-(1:burnin),3]), ylab="M_3")
dev.off()

jpeg(file=paste(path, 'M4.png', sep="/"))
traceplot(x=as.mcmc(M_keep[-(1:burnin),4]), ylab="M_4")
dev.off()

jpeg(file=paste(path, 'sgmr2.png', sep="/"))
traceplot(x=as.mcmc(sgmr2_keep[-(1:burnin)]), ylab="sgmr2")
dev.off()

jpeg(file=paste(path, 'sgm2.png', sep="/"))
traceplot(x=as.mcmc(sgm2_keep[-(1:burnin)]), ylab="sgm2")
dev.off()

jpeg(file=paste(path, 'E.png', sep="/"))
traceplot(x=as.mcmc(E_keep[-(1:burnin)]), ylab="E")
dev.off()







