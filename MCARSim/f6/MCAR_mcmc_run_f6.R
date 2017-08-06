#------------------------------------
# This script is used to run mcmc
# Last updated date: 8/2/2017
#------------------------------------
rm(list=ls())
source("/Users/peterwang/Dropbox/missingdata/Project/code/functions.R")

##running n_iter iterations
system.time(source("/Users/peterwang/Dropbox/missingdata/Project/code/MCARSim/f6/MCAR_mcmc_update_f6.R")) 

##calculate posterior mean after burn-in
burnin=3000
pm(burnin)

##traceplots
path="/Users/peterwang/Dropbox/missingdata/Project/code/MCARSim/f6/traceplots"
source("/Users/peterwang/Dropbox/missingdata/Project/code/diagnostics/traceplot.R")

##autocorrelation plots after burn-in
path="/Users/peterwang/Dropbox/missingdata/Project/code/MCARSim/f6/autocorplots"
source("/Users/peterwang/Dropbox/missingdata/Project/code/diagnostics/autocorplot.R")

##mean plots
path="/Users/peterwang/Dropbox/missingdata/Project/code/MCARSim/f6/meanplots"
source("/Users/peterwang/Dropbox/missingdata/Project/code/diagnostics/meanplot.R")

##effective sample size
source("/Users/peterwang/Dropbox/missingdata/Project/code/diagnostics/ess.R", echo = TRUE)

##mean of imputed X
MI.mean.X=apply(X_keep[-(1:burnin),], 2, mean)
##absolute relative difference with the true missing X
abs.rela.diff=abs((MI.mean.X-(datfrm$X)[R_sim==0])/((datfrm$X)[R_sim==0]))
(num_ge_1=length(which(abs.rela.diff>1)))
max(abs.rela.diff)
hist(abs.rela.diff[abs.rela.diff<1], breaks=10, xlab = "absolute relative difference")

##classification probability
classification_prob_table=matrix(0, K, n)
for(l in 1:K){
  classification_prob_table[l,]=apply(c_keep[-(1:burnin),], 2, function(x) length(which(x==l))/(n_iter-burnin))
}
write.csv(classification_prob_table, file="/Users/peterwang/Dropbox/missingdata/Project/code/MCARSim/f6/classtable.csv")









