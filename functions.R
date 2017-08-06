#------------------------------
# A collection of functions
# last updated date: 8/2/2017
#------------------------------
library(coda)

#1. given any vector x with non-negative elements, sample a value from the categorical distribution 
#which is propotional to each element
assign_c=function(x){
  l=length(x)
  if (all(x==0) ){
    latentclass=rcat(n=1, p=rep(1/l, l))
  } else {
    latentclass=rcat(n=1, p=x/sum(x))
  }
  return(latentclass)
}

#2. normalize a matrix to make the sum of each row be 1
prop_mat=function(x){
  n=nrow(x)
  m=ncol(x)
  new_mat=matrix(0, n, m)
  for(i in 1:n){
    new_mat[i, ]=prop.table(x[i, ])
  }
  return(new_mat)
}

#3. running mean plots
runmeanplot=function(x){
  meanseq=cumsum(x)/seq_along(x)
  traceplot(as.mcmc(meanseq))
}

#4. calculate posterior mean of parameters after burn-in
pm=function(burnin){
  pm.eta=apply(eta_keep[-(1:burnin),],2, mean)
  pm.M=apply(M_keep[-(1:burnin),],2, mean)
  pm.v=apply(v_keep[-(1:burnin),],2, mean)
  pm.beta=apply(beta_keep[-(1:burnin),],2, mean)
  pm.sgmr2=mean(sgmr2_keep[-(1:burnin)])
  pm.sgm2=mean(sgm2_keep[-(1:burnin)])
  pm.E=mean(E_keep[-(1:burnin)])
  return(list=c(pm.eta=pm.eta, pm.M=pm.M, pm.v=pm.v, pm.beta=pm.beta, pm.sgmr2=pm.sgmr2, pm.sgm2=pm.sgm2, pm.E=pm.E))
}











