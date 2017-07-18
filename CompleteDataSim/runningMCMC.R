#-------------------------------------------
# This file is used to run mcmc
# Last updated date: 7/12/2017
#-------------------------------------------
setwd("/Users/peterwang/Desktop/Research/missingdata/Project/code/CompleteDataSim")

##running n_iter iterations
system.time(source("mcmc.R")) 

##traceplots after burn-in
burnin=5000

jpeg(file="/Users/peterwang/Desktop/Research/missingdata/Project/code/traceplots/eta2.png")
traceplot(x=as.mcmc(eta_keep[-(1:burnin),2]), ylab="eta_2")
dev.off()

jpeg(file="/Users/peterwang/Desktop/Research/missingdata/Project/code/traceplots/eta3.png")
traceplot(x=as.mcmc(eta_keep[-(1:burnin),3]), ylab="eta_3")
dev.off()

jpeg(file="/Users/peterwang/Desktop/Research/missingdata/Project/code/traceplots/eta4.png")
traceplot(x=as.mcmc(eta_keep[-(1:burnin),4]), ylab="eta_4")
dev.off()

jpeg(file="/Users/peterwang/Desktop/Research/missingdata/Project/code/traceplots/v1.png")
traceplot(x=as.mcmc(v_keep[-(1:burnin),1]), ylab="v_1")
dev.off()

jpeg(file="/Users/peterwang/Desktop/Research/missingdata/Project/code/traceplots/v2.png")
traceplot(x=as.mcmc(v_keep[-(1:burnin),2]), ylab="v_2")
dev.off()

jpeg(file="/Users/peterwang/Desktop/Research/missingdata/Project/code/traceplots/beta1.png")
traceplot(x=as.mcmc(beta_keep[-(1:burnin),1]), ylab="beta_1")
dev.off()

jpeg(file="/Users/peterwang/Desktop/Research/missingdata/Project/code/traceplots/beta2.png")
traceplot(x=as.mcmc(beta_keep[-(1:burnin),2]), ylab="beta_2")
dev.off()

jpeg(file="/Users/peterwang/Desktop/Research/missingdata/Project/code/traceplots/M2.png")
traceplot(x=as.mcmc(M_keep[-(1:burnin),2]), ylab="M_2")
dev.off()

jpeg(file="/Users/peterwang/Desktop/Research/missingdata/Project/code/traceplots/M3.png")
traceplot(x=as.mcmc(M_keep[-(1:burnin),3]), ylab="M_3")
dev.off()

jpeg(file="/Users/peterwang/Desktop/Research/missingdata/Project/code/traceplots/M4.png")
traceplot(x=as.mcmc(M_keep[-(1:burnin),4]), ylab="M_4")
dev.off()

jpeg(file="/Users/peterwang/Desktop/Research/missingdata/Project/code/traceplots/sgmr2.png")
traceplot(x=as.mcmc(sgmr2_keep[-(1:burnin)]), ylab="sgmr2")
dev.off()

jpeg(file="/Users/peterwang/Desktop/Research/missingdata/Project/code/traceplots/sgm2.png")
traceplot(x=as.mcmc(sgm2_keep[-(1:burnin)]), ylab="sgm2")
dev.off()

jpeg(file="/Users/peterwang/Desktop/Research/missingdata/Project/code/traceplots/E.png")
traceplot(x=as.mcmc(E_keep[-(1:burnin)]), ylab="E")
dev.off()








