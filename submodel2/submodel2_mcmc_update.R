#----------------------------------------------------------
# Sub-model 2: only update beta, M, c[i], b[i], sgmr2, sgm2 
# Last updated date: 7/18/2017
#----------------------------------------------------------
setwd("/Users/peterwang/Desktop/Research/missingdata/Project/code/submodel2")
source("submodel2_data_generation.R")

#------------
# PRIORS    
#------------
beta_pri=10^{4} 
M_pri=10^{4}
sgm2_pri=0.001
sgmr2_pri=0.001
eta_pri=10^{4}

#------------------------
# SET INITIAL VALUES
#------------------------
inits=list(beta=c(-0.5, 0.8),
           M=c(0, 0.5, 1.4, -0.4),          
           sgmr2=1.5^2,
           sgm2=1.5^2, 
           b=rnorm(n, mean=0, sd=1.5)               ##initial value of random effect b_{i}
)

#-----------------------
#   SET-UP OF ITERATION                       
#-----------------------
#number of iterations
n_iter=10000

##variable names in the iteration
eta=eta_sim
beta=inits$beta
M=inits$M
sgmr2=inits$sgmr2
sgm2=inits$sgm2
c=vector(length=n)
b=inits$b
X=X_sim

##recording structure, each row is one iteration
beta_keep=matrix(0, nrow=n_iter, ncol=2)  
M_keep=matrix(0, nrow=n_iter, ncol=K)    
sgmr2_keep=rep(0, n_iter)  
sgm2_keep=rep(0, n_iter)   
c_keep=matrix(0, nrow=n_iter, ncol=n) 
b_keep=matrix(0, nrow=n_iter, ncol=n) 

## for updating c 
Pi_c=exp(V_sim%o%eta)
for(i in 1:n){
  Pi_c[i, ]=prop.table(Pi_c[i,]) 
}

f_c=matrix(0,n, K)          ## for updating c
crossD=crossprod(D)         ## t(D)%*%D for updating beta

#------------------
# RUN ITERATIONS
#------------------
for (m in 1:n_iter){  
  
  ##sample c_{i} for all i
  for(i in 1:n){
    for (l in 1:K){
      pos=(sum(T[1:i])-T[i]+1):sum(T[1:i])
      f_c[i,l]=dmvnorm(X[pos], mean=D[pos,]%*%beta+D_star[pos]*M[l]+D_dstar[pos]*b[i],
                       sigma=sgm2*diag(T[i]) )
    }
  }
  
  for(i in 1:n){
    prob=(Pi_c*f_c)[i, ]
    if (all(prob==0) ){
      c[i]=rcat(n=1, p=rep(1/K, K))
    } else {
      c[i]=rcat(n=1, p=prob/sum(prob) )
    }
  }
  
  ##sample M2, M3 and M4
  Bigb=rep(b, T)
  for(l in 2:K){
    c_index=rep(0, n)
    index=which(c==l)
    c_index[index]=1
    var_M=((1/M_pri)+(1/sgm2)*sum((D_star*rep(c_index,T))^2))^{-1}
    mean_M=(1/sgm2)*var_M*sum((D_star*(X-D%*%beta-D_dstar*Bigb))*rep(c_index,T))
    M[l]=rnorm(1, mean=mean_M, sd=sqrt(var_M))
  }
  
  ##sample beta
  BigM=rep(M[c],T)
  sum_beta=crossprod(D, X-D_star*BigM-D_dstar*Bigb)
  var_beta=solve((1/beta_pri)*diag(2)+(1/sgm2)*crossD)
  mean_beta=(1/sgm2)*(var_beta%*%sum_beta)
  beta=mvrnorm(n=1, mu=mean_beta, Sigma = var_beta)
  
  ##sample b[i]
  for (i in 1:n){
    b_index=c(rep(0,i-1),1,rep(0, n-i))
    var_b=((1/sgmr2)+(1/sgm2)*T[i])^{-1}
    mean_b=(1/sgm2)*var_b*sum((X-D%*%beta-D_star*BigM)*rep(b_index,T))
    b[i]=rnorm(1, mean=mean_b, sd=sqrt(var_b))
  }  
  
  ##sample sigma_{r}^{2}
  sgmr2=rigamma(n=1, a=(n/2)+sgmr2_pri, b=(1/2)*sum(b^{2})+sgmr2_pri)
  
  ##sample sigma^{2}
  Bigb=rep(b, T)
  sum_sgm2=sum((X-D%*%beta-D_star*BigM-D_dstar*Bigb)^2)
  shape_sgm2=(1/2)*sum(T)+sgm2_pri
  scale_sgm2=(1/2)*sum_sgm2+sgm2_pri
  sgm2=rigamma(n=1, a=shape_sgm2, b=scale_sgm2)
  
  ##record parameters
  c_keep[m, ]=c
  b_keep[m, ]=b
  beta_keep[m, ]=beta
  M_keep[m, ]=M
  sgmr2_keep[m]=sgmr2
  sgm2_keep[m]=sgm2

} ##iteration ends


