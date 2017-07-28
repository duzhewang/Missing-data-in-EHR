#--------------------------------------------
# Sub-model 3: only update v, e_{i}, and E 
# Last updated date: 7/25/2017
#--------------------------------------------
setwd("/Users/peterwang/Desktop/Research/missingdata/Project/code/submodel3")
source("submodel3_data_generation.R")

#------------
# PRIORS    
#------------
v_pri=10^{4}
E_pri=0.001

#------------------------
# SET INITIAL VALUES
#------------------------
inits=list(v=c(0.3, -1),
           E=1.5^2, 
           e=rnorm(n, mean=0, sd=1.5) ##initial value of random effect e_{i}
)

#-----------------------
#   SET-UP OF ITERATION                       
#-----------------------
#number of iterations
n_iter=10000

##variable names in the iteration
v=inits$v
E=inits$E
e=inits$e
X=X_sim

##recording structure, each row is one iteration
v_keep=matrix(0, nrow=n_iter, ncol=2)    
E_keep=rep(0, n_iter)     
e_keep=matrix(0, nrow=n_iter, ncol=n) 

omega_v=rep(0, sum(T))      ##for updating v

#------------------
# RUN ITERATIONS
#------------------
for (m in 1:n_iter){  
  
  ##sample v
  tilting_v=BigVX%*%v+rep(e,T)
  for (i in 1:sum(T)){
    omega_v[i]=rpg(num=1, h=1, z=tilting_v[i])  ## sample w_{it}^{*}
  }
  k_v=Y_sim-omega_v*rep(e, T)-1/2
  S_v=solve((1/v_pri)*diag(2)+crossprod(BigVX,diag(omega_v))%*%BigVX )
  m_v=S_v%*%crossprod(BigVX,k_v)
  v=mvrnorm(n=1, mu=m_v, Sigma = S_v)
  
  ##sample e_{i}
  for(i in 1:n){
    e_index=c(rep(0,i-1),1,rep(0, n-i))
    S_e=((1/E)+sum(omega_v*rep(e_index, T)))^{-1}
    m_e=S_e*sum((Y_sim-omega_v*(BigVX%*%v)-1/2)*rep(e_index, T))
    e[i]=rnorm(n=1, mean=m_e, sd=sqrt(S_e))
  }
  
  ##sample E
  E=rigamma(n=1, a=n/2+E_pri, b=(1/2)*sum(e^{2})+E_pri)
  
  ##record parameters
  e_keep[m, ]=e
  v_keep[m, ]=v
  E_keep[m]=E
  
} ##iteration ends

