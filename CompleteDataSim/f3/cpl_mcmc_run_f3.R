#----------------------------------
# this script is used to run mcmc 
# Last updated date: 8/16/2017
#----------------------------------
rm(list=ls())
source("/Users/peterwang/Dropbox/missingdata/Project/code/functions.R")

##run n_iter iterations
system.time(source("/Users/peterwang/Dropbox/missingdata/Project/code/CompleteDataSim/f3/cpl_mcmc_update_f3.R")) 

##save results 
write.csv(c_keep, file="/Users/peterwang/Dropbox/missingdata/Project/code/CompleteDataSim/f3/c_keep.csv")
write.csv(b_keep, file="/Users/peterwang/Dropbox/missingdata/Project/code/CompleteDataSim/f3/b_keep.csv")
write.csv(e_keep, file="/Users/peterwang/Dropbox/missingdata/Project/code/CompleteDataSim/f3/e_keep.csv")
write.csv(eta_keep, file="/Users/peterwang/Dropbox/missingdata/Project/code/CompleteDataSim/f3/eta_keep.csv")
write.csv(v_keep, file="/Users/peterwang/Dropbox/missingdata/Project/code/CompleteDataSim/f3/v_keep.csv")
write.csv(beta_keep, file="/Users/peterwang/Dropbox/missingdata/Project/code/CompleteDataSim/f3/beta_keep.csv")
write.csv(M_keep, file="/Users/peterwang/Dropbox/missingdata/Project/code/CompleteDataSim/f3/M_keep.csv")
write.csv(sgmr2_keep, file="/Users/peterwang/Dropbox/missingdata/Project/code/CompleteDataSim/f3/sgmr2_keep.csv")
write.csv(sgm2_keep, file="/Users/peterwang/Dropbox/missingdata/Project/code/CompleteDataSim/f3/sgm2_keep.csv")
write.csv(E_keep, file="/Users/peterwang/Dropbox/missingdata/Project/code/CompleteDataSim/f3/E_keep.csv")

##calculate posterior mean after burn-in
burnin=3000
pm(burnin)

##traceplots
path="/Users/peterwang/Dropbox/missingdata/Project/code/CompleteDataSim/f3/traceplots"
source("/Users/peterwang/Dropbox/missingdata/Project/code/diagnostics/traceplot.R")

##autocorrelation plots after burn-in
path="/Users/peterwang/Dropbox/missingdata/Project/code/CompleteDataSim/f3/autocorplots"
source("/Users/peterwang/Dropbox/missingdata/Project/code/diagnostics/autocorplot.R")

##mean plots
path="/Users/peterwang/Dropbox/missingdata/Project/code/CompleteDataSim/f3/meanplots"
source("/Users/peterwang/Dropbox/missingdata/Project/code/diagnostics/meanplot.R")

##effective sample size
source("/Users/peterwang/Dropbox/missingdata/Project/code/diagnostics/ess.R", echo = TRUE)




