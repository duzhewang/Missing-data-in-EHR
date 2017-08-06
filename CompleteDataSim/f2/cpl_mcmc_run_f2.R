#----------------------------------
# this script is used to run mcmc 
# Last updated date: 8/2/2017
#----------------------------------
rm(list=ls())
source("/Users/peterwang/Dropbox/missingdata/Project/code/functions.R")

##run n_iter iterations
system.time(source("/Users/peterwang/Dropbox/missingdata/Project/code/CompleteDataSim/f2/cpl_mcmc_update_f2.R")) 

##calculate posterior mean after burn-in
burnin=3000
pm(burnin)

##traceplots
path="/Users/peterwang/Dropbox/missingdata/Project/code/CompleteDataSim/f2/traceplots"
source("/Users/peterwang/Dropbox/missingdata/Project/code/diagnostics/traceplot.R")

##autocorrelation plots after burn-in
path="/Users/peterwang/Dropbox/missingdata/Project/code/CompleteDataSim/f2/autocorplots"
source("/Users/peterwang/Dropbox/missingdata/Project/code/diagnostics/autocorplot.R")

##mean plots
path="/Users/peterwang/Dropbox/missingdata/Project/code/CompleteDataSim/f2/meanplots"
source("/Users/peterwang/Dropbox/missingdata/Project/code/diagnostics/meanplot.R")

##effective sample size
source("/Users/peterwang/Dropbox/missingdata/Project/code/diagnostics/ess.R", echo = TRUE)
