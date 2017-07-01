#-------------------------------------------
# This file is for testing MCMC running time 
# and MCMC diagnostics
# Last updated date: 7/1/2017
#-------------------------------------------
setwd("/Users/peterwang/Desktop/Research/missingdata/Project/code")
library(coda)

##running n_iter iterations
system.time(source("PolyaGamma.R")) 

##export the results
write.csv(eta_keep, file="/Users/peterwang/Desktop/Research/missingdata/Project/code/results/eta.csv")
write.csv(v_keep, file="/Users/peterwang/Desktop/Research/missingdata/Project/code/results/v.csv")
write.csv(beta_keep, file="/Users/peterwang/Desktop/Research/missingdata/Project/code/results/beta.csv")
write.csv(M_keep, file="/Users/peterwang/Desktop/Research/missingdata/Project/code/results/M.csv")
write.csv(sgmr2_keep, file="/Users/peterwang/Desktop/Research/missingdata/Project/code/results/sgmr2.csv")
write.csv(sgm2_keep, file="/Users/peterwang/Desktop/Research/missingdata/Project/code/results/sgm2.csv")
write.csv(E_keep, file="/Users/peterwang/Desktop/Research/missingdata/Project/code/results/E.csv")
write.csv(c_keep, file="/Users/peterwang/Desktop/Research/missingdata/Project/code/results/c.csv")
write.csv(b_keep, file="/Users/peterwang/Desktop/Research/missingdata/Project/code/results/b.csv")
write.csv(e_keep, file="/Users/peterwang/Desktop/Research/missingdata/Project/code/results/e.csv")


##traceplot after burn-in
burnin=5000
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














