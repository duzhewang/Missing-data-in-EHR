#------------------------------------
# This script is used to run mcmc
# Last updated date: 7/15/2017
#------------------------------------
setwd("/Users/peterwang/Desktop/Research/missingdata/Project/code/MCARSim")

##running n_iter iterations
system.time(source("mcmc_update.R")) 

##traceplots after burn-in
burnin=5000

jpeg(file="/Users/peterwang/Desktop/Research/missingdata/Project/code/MCARSim/traceplots/eta2.png")
traceplot(x=as.mcmc(eta_keep[-(1:burnin),2]), ylab="eta_2")
dev.off()

jpeg(file="/Users/peterwang/Desktop/Research/missingdata/Project/code/MCARSim/traceplots/eta3.png")
traceplot(x=as.mcmc(eta_keep[-(1:burnin),3]), ylab="eta_3")
dev.off()

jpeg(file="/Users/peterwang/Desktop/Research/missingdata/Project/code/MCARSim/traceplots/eta4.png")
traceplot(x=as.mcmc(eta_keep[-(1:burnin),4]), ylab="eta_4")
dev.off()

jpeg(file="/Users/peterwang/Desktop/Research/missingdata/Project/code/MCARSim/traceplots/v1.png")
traceplot(x=as.mcmc(v_keep[-(1:burnin),1]), ylab="v_1")
dev.off()

jpeg(file="/Users/peterwang/Desktop/Research/missingdata/Project/code/MCARSim/traceplots/v2.png")
traceplot(x=as.mcmc(v_keep[-(1:burnin),2]), ylab="v_2")
dev.off()

jpeg(file="/Users/peterwang/Desktop/Research/missingdata/Project/code/MCARSim/traceplots/beta1.png")
traceplot(x=as.mcmc(beta_keep[-(1:burnin),1]), ylab="beta_1")
dev.off()

jpeg(file="/Users/peterwang/Desktop/Research/missingdata/Project/code/MCARSim/traceplots/beta2.png")
traceplot(x=as.mcmc(beta_keep[-(1:burnin),2]), ylab="beta_2")
dev.off()

jpeg(file="/Users/peterwang/Desktop/Research/missingdata/Project/code/MCARSim/traceplots/M2.png")
traceplot(x=as.mcmc(M_keep[-(1:burnin),2]), ylab="M_2")
dev.off()

jpeg(file="/Users/peterwang/Desktop/Research/missingdata/Project/code/MCARSim/traceplots/M3.png")
traceplot(x=as.mcmc(M_keep[-(1:burnin),3]), ylab="M_3")
dev.off()

jpeg(file="/Users/peterwang/Desktop/Research/missingdata/Project/code/MCARSim/traceplots/M4.png")
traceplot(x=as.mcmc(M_keep[-(1:burnin),4]), ylab="M_4")
dev.off()

jpeg(file="/Users/peterwang/Desktop/Research/missingdata/Project/code/MCARSim/traceplots/sgmr2.png")
traceplot(x=as.mcmc(sgmr2_keep[-(1:burnin)]), ylab="sgmr2")
dev.off()

jpeg(file="/Users/peterwang/Desktop/Research/missingdata/Project/code/MCARSim/traceplots/sgm2.png")
traceplot(x=as.mcmc(sgm2_keep[-(1:burnin)]), ylab="sgm2")
dev.off()

jpeg(file="/Users/peterwang/Desktop/Research/missingdata/Project/code/MCARSim/traceplots/E.png")
traceplot(x=as.mcmc(E_keep[-(1:burnin)]), ylab="E")
dev.off()








