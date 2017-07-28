#------------------------------------
#  Run Sub-model 3 with missing value
#  Last updated date: 7/25/2017
#------------------------------------
setwd("/Users/peterwang/Desktop/Research/missingdata/Project/code/submodel3")

##run n_iter iterations
system.time(source("submodel3_mcmc_update_case2.R")) 

#evaluate the method of imputing initial value of missing X_{it}
maxdiff.ini.X
mindiff.ini.X

burnin=5000
##posterior mean
(posterior.mean.v=apply(v_keep[-(1:burnin),],2, mean))
(posterior.mean.E=mean(E_keep[-(1:burnin)]))


##mean of imputed X
MI.mean.X=apply(X_keep[-(1:burnin),], 2, mean)
##difference with the true X
(diff=MI.mean.X-(SVXYR$X)[R_sim==0])
max(diff)
min(diff)


##traceplots after burn-in
traceplot(x=as.mcmc(v_keep[-(1:burnin),1]), ylab="v_1")
traceplot(x=as.mcmc(v_keep[-(1:burnin),2]), ylab="v_2")
traceplot(x=as.mcmc(E_keep[-(1:burnin)]), ylab="E")


