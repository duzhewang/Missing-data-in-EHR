#------------------------------------------------------
# lmm_mcmc_update.R is used to update mcmc iteration
# Last updated date: 7/18/2017
#------------------------------------------------------
setwd("/Users/peterwang/Desktop/Research/missingdata/Project/code/LinearMixed")
source("lmm_data_generation.R")

#------------------------
# SET INITIAL VALUES
#------------------------
inits=list(beta=c(-0.5, 0.8),
           sgmr2=sgmr2ini,
           sgm2=sgm2ini, 
           b=rnorm(n, mean=0, sd=sqrt(sgmr2ini))     ##initial value of random effect b_{i}
)

#-----------------------
#   SET-UP OF ITERATION                       
#-----------------------
#number of iterations
n_iter=10000

##variable names in the iteration
beta=inits$beta
sgmr2=inits$sgmr2
sgm2=inits$sgm2
b=inits$b
X=X_sim

##recording structure, each row is one iteration
beta_keep=matrix(0, nrow=n_iter, ncol=2)  
sgmr2_keep=rep(0, n_iter)  
sgm2_keep=rep(0, n_iter)   
b_keep=matrix(0, nrow=n_iter, ncol=n) 

crossD=crossprod(D)         ## t(D)%*%D for updating beta

#------------------
# RUN ITERATIONS
#------------------
for (m in 1:n_iter){  
  
  ##sample beta
  Bigb=rep(b, T)
  sum_beta=crossprod(D, X-D_dstar*Bigb)
  var_beta=solve((1/beta_pri)*diag(2)+(1/sgm2)*crossD)
  mean_beta=(1/sgm2)*(var_beta%*%sum_beta)
  beta=mvrnorm(n=1, mu=mean_beta, Sigma = var_beta)
  
  ##sample sigma^{2}
  sum_sgm2=sum((X-D%*%beta-D_dstar*Bigb)^2)
  shape_sgm2=(1/2)*sum(T)+sgm2_pri
  scale_sgm2=(1/2)*sum_sgm2+sgm2_pri
  sgm2=rigamma(n=1, a=shape_sgm2, b=scale_sgm2)
  
  ##sample sigma_{r}^{2}
  sgmr2=rigamma(n=1, a=(n/2)+sgmr2_pri, b=(1/2)*sum(b^{2})+sgmr2_pri) 
  
  ##sample b[i]
  for (i in 1:n){
    b_index=c(rep(0,i-1),1,rep(0, n-i))
    var_b=((1/sgmr2)+(1/sgm2)*T[i])^{-1}
    mean_b=(1/sgm2)*var_b*sum((X-D%*%beta)*rep(b_index,T))
    b[i]=rnorm(1, mean=mean_b, sd=sqrt(var_b))
  }  
  

  ##record parameters
  b_keep[m, ]=b
  beta_keep[m, ]=beta
  sgmr2_keep[m]=sgmr2
  sgm2_keep[m]=sgm2

} ##iteration ends


