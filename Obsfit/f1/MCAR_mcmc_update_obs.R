#------------------------------------------------------------------------------
# This script is used to update mcmc using observed data only and NOT impute x_{it}
# Last updated date: 8/7/2017
#------------------------------------------------------------------------------
source("/Users/peterwang/Dropbox/missingdata/Project/code/Obsfit/f1/MCAR_data_generation_obs.R")

#number of iterations
n_iter=5000
#----------------------------
# PRIORS AND INITIAL VALUES    
#----------------------------
beta_pri=10^{4} 
M_pri=10^{4}
sgm2_pri=0.001
sgmr2_pri=0.001
v_pri=10^{4}
E_pri=0.001
eta_pri=10^{4}

inits=list(eta=c(0, 0.2, 1.3, 0.7), 
           beta=c(-0.5, 0.8),
           M=c(0, 0.5, 1.4, -0.4),          
           v=c(0.3, -1),
           sgmr2=1.5^2,
           sgm2=1.5^2, 
           E=1.5^2, 
           b=rnorm(n, mean=0, sd=1.5),             ##initial value of random effect b_{i}
           e=rnorm(n, mean=0, sd=1.5)              ##initial value of random effect e_{i}
)

#-----------------------
#   SET-UP OF ITERATION                       
#-----------------------
##variable names in the iteration
c=vector(length=n)
eta=inits$eta
M=inits$M
beta=inits$beta
b=inits$b
sgmr2=inits$sgmr2
sgm2=inits$sgm2
v=inits$v
e=inits$e
E=inits$E
X=X_sim

##recording structure, each row is one iteration
eta_keep=matrix(0, nrow=n_iter, ncol=K)  
v_keep=matrix(0, nrow=n_iter, ncol=2)    
beta_keep=matrix(0, nrow=n_iter, ncol=2)  
M_keep=matrix(0, nrow=n_iter, ncol=K)    
sgmr2_keep=rep(0, n_iter)  
sgm2_keep=rep(0, n_iter)   
E_keep=rep(0, n_iter)     
c_keep=matrix(0, nrow=n_iter, ncol=n) 
b_keep=matrix(0, nrow=n_iter, ncol=n) 
e_keep=matrix(0, nrow=n_iter, ncol=n) 

Pi_c=matrix(0, n, K)        ## for updating c 
f_c=matrix(0,n, K)          ## for updating c
LD=matrix(0, sum(T), 2)     ## for updating beta
La=vector(length=sum(T))    ## for updating beta
var_b=vector(length=n)      ## for updating b
mean_b=vector(length=n)     ## for updating b
omega_v=rep(0, sum(T))      ## for updating v
BigVX=matrix(0, sum(T), 2)  ## for updating v
BigVX[, 1]=rep(V_sim,T)
S_e=vector(length=n)        ## for updating e
m_e=vector(length=n)        ## for updating e

#------------------
# RUN ITERATIONS
#------------------
for (m in 1:n_iter){  
  
  ##sample c_{i} for all i
  Pi_c=prop_mat(exp(V_sim%o%eta))
  
  for(i in 1:n){
    for (l in 1:K){
      pos=(sum(T[1:i])-T[i]+1):sum(T[1:i])
      f_c[i,l]=dmvnorm(X[pos], mean=D[pos,]%*%beta+D_star[pos]*M[l]+rep(b[i],T[i]),
                       sigma=sgm2*diag(T[i]) )
    }
  }
  
  for(i in 1:n){
    prob=(Pi_c*f_c)[i, ]
    c[i]=assign_c(prob)
  }
  
  ##sample eta2, eta3 and eta4
  for(l in 2:K){
    logsum=log(apply(exp(V_sim%o%eta[-l]),1,sum))
    r=V_sim*eta[l]-logsum
    w=rpg(num=n, h=1, z=r)
    k_eta=as.numeric(c==l)-1/2+w*logsum
    S_eta=(1/eta_pri+t(V_sim)%*%diag(w)%*%V_sim)^{-1}
    m_eta=S_eta*crossprod(V_sim, k_eta)    
    eta[l]=rnorm(n=1,mean=m_eta, sd=sqrt(S_eta))
  }
  
  ##sample M2, M3 and M4
  Bigb=rep(b, T)
  for(l in 2:K){
    c_index=as.numeric(c==l)  ##select i s.t. c[i]=l
    pos=which(rep(c_index, T)==1)
    var_M=((1/M_pri)+(1/sgm2)*sum((D_star[pos])^2))^{-1}
    mean_M=(1/sgm2)*var_M*crossprod(D_star[pos],(X-D%*%beta-D_dstar*Bigb)[pos])
    M[l]=rnorm(1, mean=mean_M, sd=sqrt(var_M))
  }
  
  ##sample beta: use Cholesky decomposition
  for (i in 1:n){
    posD=which(datfrm$Subject==i)
    InvW=sgmr2*matrix(1,T[i],T[i])+sgm2*diag(T[i])
    InvCholW=solve(t(chol(InvW)))
    LD[posD,]=InvCholW%*%D[posD, ]
    La[posD]=InvCholW%*%(X[posD]-D_star[posD]*M[c[i]])
  }
  var_beta=solve((1/beta_pri)*diag(2)+ crossprod(LD))
  mean_beta=var_beta%*%crossprod(LD, La)
  beta=as.vector(rmvnorm(n=1, mean=mean_beta, sigma=var_beta))

  ##sample b[i]
  BigM=rep(M[c],T)
  for (i in 1:n){
    pos=which(datfrm$Subject==i)
    var_b[i]=((1/sgmr2)+(1/sgm2)*T[i])^{-1}
    mean_b[i]=(1/sgm2)*var_b[i]*sum((X-D%*%beta-D_star*BigM)[pos])
  }  
  b=rnorm(n)*sqrt(var_b)+mean_b
  
  ##sample sigma_{r}^{2}
  sgmr2=rigamma(n=1, a=(n/2)+sgmr2_pri, b=(1/2)*crossprod(b)+sgmr2_pri)
  
  ##sample sigma^{2}
  Bigb=rep(b, T)
  sum_sgm2=crossprod(X-D%*%beta-D_star*BigM-D_dstar*Bigb)
  shape_sgm2=(1/2)*sum(T)+sgm2_pri
  scale_sgm2=(1/2)*sum_sgm2+sgm2_pri
  sgm2=rigamma(n=1, a=shape_sgm2, b=scale_sgm2)
  
  ##sample v
  BigVX[,2]=X
  tilting_v=BigVX%*%v+rep(e,T)
  omega_v=rpg(num=sum(T), h=1, z=tilting_v)      ## sample w_{it}^{*}
  k_v=Y_sim-omega_v*rep(e, T)-1/2
  S_v=solve((1/v_pri)*diag(2)+crossprod(BigVX,diag(omega_v))%*%BigVX)
  m_v=S_v%*%crossprod(BigVX,k_v)
  v=as.vector(rmvnorm(n=1, mean=m_v, sigma=S_v))

  ##sample e_{i}
  for(i in 1:n){
    pos=which(datfrm$Subject==i)
    S_e[i]=((1/E)+sum(omega_v[pos]))^{-1}
    m_e[i]=S_e[i]*sum(Y_sim[pos]-1/2-omega_v[pos]*(BigVX[pos, ]%*%v))
  }
  e=rnorm(n)*sqrt(S_e)+m_e
  
  ##sample E
  E=rigamma(n=1, a=n/2+E_pri, b=(1/2)*crossprod(e)+E_pri)

  ##record parameters
  c_keep[m, ]=c
  b_keep[m, ]=b
  e_keep[m, ]=e
  eta_keep[m,]=eta
  v_keep[m, ]=v
  beta_keep[m, ]=beta
  M_keep[m, ]=M
  sgmr2_keep[m]=sgmr2
  sgm2_keep[m]=sgm2
  E_keep[m]=E
 
} ##iteration ends

