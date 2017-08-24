#-----------------------------------------
# This script is used to run mcmc in MAR
# Last updated date: 8/18/2017
#-----------------------------------------
rm(list=ls())
source("/Users/dwang282/Desktop/f4/functions.R")

##running n_iter iterations
system.time(source("/Users/dwang282/Desktop/f4/MAR_mcmc_update_f4.R")) 

##save results 
write.csv(c_keep, file="/Users/dwang282/Desktop/f4/c_keep.csv")
write.csv(b_keep, file="/Users/dwang282/Desktop/f4/b_keep.csv")
write.csv(e_keep, file="/Users/dwang282/Desktop/f4/e_keep.csv")
write.csv(eta_keep, file="/Users/dwang282/Desktop/f4/eta_keep.csv")
write.csv(v_keep, file="/Users/dwang282/Desktop/f4/v_keep.csv")
write.csv(beta_keep, file="/Users/dwang282/Desktop/f4/beta_keep.csv")
write.csv(M_keep, file="/Users/dwang282/Desktop/f4/M_keep.csv")
write.csv(sgmr2_keep, file="/Users/dwang282/Desktop/f4/sgmr2_keep.csv")
write.csv(sgm2_keep, file="/Users/dwang282/Desktop/f4/sgm2_keep.csv")
write.csv(E_keep, file="/Users/dwang282/Desktop/f4/E_keep.csv")
write.csv(X_keep, file="/Users/dwang282/Desktop/f4/X_keep.csv")

##calculate posterior mean after burn-in
burnin=3000
pm(burnin)

##traceplots
path="/Users/dwang282/Desktop/f4/traceplots"
source("/Users/dwang282/Desktop/f4/diagnostics/traceplot.R")



##autocorrelation plots after burn-in
path="/Users/dwang282/Desktop/f4/autocorplots"
source("/Users/dwang282/Desktop/f4/diagnostics/autocorplot.R")

##mean plots
path="/Users/dwang282/Desktop/f4/meanplots"
source("/Users/dwang282/Desktop/f4/diagnostics/meanplot.R")

##effective sample size
source("/Users/dwang282/Desktop/f4/diagnostics/ess.R", echo = TRUE)

##mean of imputed X
MI.mean.X=apply(X_keep[-(1:burnin),], 2, mean)
##absolute relative difference with the true missing X
abs.rela.diff=abs((MI.mean.X-(datfrm$X)[R_sim==0])/((datfrm$X)[R_sim==0]))
(num_missing=length(R_sim)-sum(R_sim))
(overall_mr=num_missing/length(R_sim))                ## overall missing rate
(num_ge_1=length(which(abs.rela.diff>1)))
(num_le_1=length(which(abs.rela.diff<1)))
max(abs.rela.diff)
hist(abs.rela.diff[abs.rela.diff<1], breaks=10, xlab = "absolute relative difference")

##classification probability
classification_prob_table=matrix(0, K, n)
for(l in 1:K){
  classification_prob_table[l,]=apply(c_keep[-(1:burnin),], 2, function(x) length(which(x==l))/(n_iter-burnin))
}
write.csv(classification_prob_table, file="/Users/dwang282/Desktop/f4/classtable.csv")









