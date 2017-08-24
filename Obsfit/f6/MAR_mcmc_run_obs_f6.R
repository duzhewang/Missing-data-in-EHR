#-----------------------------------------
# This script is used to run mcmc
# Last updated date: 8/20/2017
#-----------------------------------------
rm(list=ls())
source("/Users/peterwang/Dropbox/missingdata/Project/code/functions.R")

##running n_iter iterations
system.time(source("/Users/peterwang/Dropbox/missingdata/Project/code/Obsfit/f6/MAR_mcmc_update_obs_f6.R")) 

##save results 
write.csv(c_keep, file="/Users/peterwang/Dropbox/missingdata/Project/code/Obsfit/f6/c_keep.csv")
write.csv(b_keep, file="/Users/peterwang/Dropbox/missingdata/Project/code/Obsfit/f6/b_keep.csv")
write.csv(e_keep, file="/Users/peterwang/Dropbox/missingdata/Project/code/Obsfit/f6/e_keep.csv")
write.csv(eta_keep, file="/Users/peterwang/Dropbox/missingdata/Project/code/Obsfit/f6/eta_keep.csv")
write.csv(v_keep, file="/Users/peterwang/Dropbox/missingdata/Project/code/Obsfit/f6/v_keep.csv")
write.csv(beta_keep, file="/Users/peterwang/Dropbox/missingdata/Project/code/Obsfit/f6/beta_keep.csv")
write.csv(M_keep, file="/Users/peterwang/Dropbox/missingdata/Project/code/Obsfit/f6/M_keep.csv")
write.csv(sgmr2_keep, file="/Users/peterwang/Dropbox/missingdata/Project/code/Obsfit/f6/sgmr2_keep.csv")
write.csv(sgm2_keep, file="/Users/peterwang/Dropbox/missingdata/Project/code/Obsfit/f6/sgm2_keep.csv")
write.csv(E_keep, file="/Users/peterwang/Dropbox/missingdata/Project/code/Obsfit/f6/E_keep.csv")

##calculate posterior mean after burn-in
burnin=3000
pm(burnin)

##traceplots
path="/Users/peterwang/Dropbox/missingdata/Project/code/Obsfit/f6/traceplots"
source("/Users/peterwang/Dropbox/missingdata/Project/code/diagnostics/traceplot.R")

##autocorrelation plots after burn-in
path="/Users/peterwang/Dropbox/missingdata/Project/code/Obsfit/f6/autocorplots"
source("/Users/peterwang/Dropbox/missingdata/Project/code/diagnostics/autocorplot.R")

##mean plots
path="/Users/peterwang/Dropbox/missingdata/Project/code/Obsfit/f6/meanplots"
source("/Users/peterwang/Dropbox/missingdata/Project/code/diagnostics/meanplot.R")

##effective sample size
source("/Users/peterwang/Dropbox/missingdata/Project/code/diagnostics/ess.R", echo = TRUE)

##classification probability
classification_prob_table=matrix(0, K, n)
for(l in 1:K){
  classification_prob_table[l,]=apply(c_keep[-(1:burnin),], 2, function(x) length(which(x==l))/(n_iter-burnin))
}
write.csv(classification_prob_table, file="/Users/peterwang/Dropbox/missingdata/Project/code/Obsfit/f6/classtable.csv")









