#----------------------------------
# BAYESIAN METHOD SIMULATIOIN 
# Last updated date: 7/1/2017
#----------------------------------
rm(list=ls())
set.seed(0311)
library(MASS) ##use function mvrnorm
library(LaplacesDemon) ##use function rcat
library(MCMCpack) ##use function rinvgamma, default rate=1
library(mvtnorm)  ##use function dmvnorm, rmvt
library(pscl)  ##use function rigamma 
#library(metRology)  ##use function rt.scaled 
library(BayesLogit) ##use rpg

#---------------------------------
#      Global set-up       
#---------------------------------
n=100 ##number of subject
K=4 ##number of latent classes

##assign each subject the number of tracked quarters
T=sample(5:40, size=n, replace=TRUE) 

#--------------------------------------
#   Complete simulated data generation
#--------------------------------------
##set up true values of parameters
eta_sim=c(0, 0.7, 1.6, 0.9)
beta_sim=c(-0.4, 0.5)
M_sim=c(0, -0.6, 0.6, 1.2) 
v_sim=c(0.5, -0.3)
sgmr2_sim=1                ##true value of variance of b_{i}
sgm2_sim=1                 ##true value of variance of epsilon_{it}
E_sim=1                    ##true value of variance of e_{i}

##generate time-invariate covariate V_{i}
V_sim=rnorm(n=n, mean=0, sd=1)

##simulate latent class c[i] for each subject 
Pi_sim=matrix(0, n, K )   ##Pi[i,l]=P(c_{i}=l|V_{i})
for (i in 1:n){
  for (l in 1:K){
    Pi_sim[i,l]=exp(V_sim[i]*eta_sim[l])/sum(exp(V_sim[i]*eta_sim))
  }
}
c_sim=apply(Pi_sim, 1, function(x) rcat(n=1, p=x)) ##assign latent class to each subject

##generate random effect b_{i}
b_sim=rnorm(n,mean=0, sd=sqrt(sgmr2_sim))  

##set D_{i}, D^{*}_{i} and D^{**}_{i}
D=matrix(0, max(T), 2*n)              ##each two columns for one subject
for (i in 1:n){
  D[1:T[i], 2*i-1]=rep(1, T[i])
  D[1:T[i], 2*i]=rep(V_sim[i], T[i])
}

D_star=matrix(0, max(T), n)   
for (i in 1:n){
  D_star[1:T[i],i]=1:T[i]
}

D_dstar=matrix(0, max(T), n)
for (i in 1:n){
  D_dstar[1:T[i],i]=rep(1, T[i])
}

##sum_D is used in updating beta
sum_D=matrix(0, 2, 2)
for (i in 1:n){
  sum_D=sum_D+t(D[1:T[i],c(2*i-1, 2*i)])%*%D[1:T[i], c(2*i-1, 2*i)] 
}

##simulate X_{it}
X_sim=matrix(0, nrow=max(T), ncol = n)  ##each column is one subject
for (i in 1:n){
  X_sim[1:T[i],i]=D[1:T[i],c(2*i-1,2*i)]%*%beta_sim+
                  D_star[1:T[i],i]*M_sim[c_sim[i]]+
                  D_dstar[1:T[i],i]*b_sim[i]+          
                  rnorm(n=T[i], mean = 0, sd=sqrt(sgm2_sim))
}

##generate random effect e_{i}
e_sim=rnorm(n, mean=0, sd=sqrt(E_sim))

##generate y_{it}
Y_sim=matrix(0,nrow=max(T), ncol=n)     ##each column is one subject
for (i in 1:n){
  top= exp(cbind(D[1:T[i],2*i], X_sim[1:T[i],i])%*%v_sim+rep(e_sim[i],T[i]))
  bot=rep(1, T[i])+top
  Y_sim[1:T[i],i]=rbinom(n=T[i], size=1, p=top/bot)   
}


#-------------------------------------------
#          PRIORS AND INITIAL VALUES    
#-------------------------------------------
beta_pri=10^{4} 
M_pri=10^{4}
sgm2_pri=0.001
sgmr2_pri=0.001
v_pri=10^{4}
E_pri=0.001
eta_pri=10^{4}

##initial values
inits=list(eta1=0, eta2=0.5, eta3=1.3, eta4=1,
           beta=c(-0.5, 0.8),
           M1=0, M2=-1, M3=1.2, M4=0.6, 
           v=c(0.3, -1),
           sgmr2=1.5^2,
           sgm2=1.5^2, 
           E=1.5^2, 
           b=rnorm(n, mean=0, sd=1.5),             ##initial value of random effect b_{i}
           e=rnorm(n, mean=0, sd=1.5)              ##initial value of random effect e_{i}
)

#set initial latent class for each subject
eta_ini=c(inits$eta1,inits$eta2,inits$eta3,inits$eta4)
Pi_ini=matrix(0, n, K )  
for (i in 1:n){
  for (l in 1:K){
    Pi_ini[i,l]=exp(V_sim[i]*eta_ini[l])/sum(exp(V_sim[i]*eta_ini))
  }
}
c_ini=apply(Pi_ini, 1, function(x) rcat(n=1, p=x)) 

#--------------------------------------------------
#        GIBBS SAMPLING-ITERATIONS                       
#--------------------------------------------------
n_iter=10000         ##number of iterations

##variable names in the iteration
eta=c(inits$eta1,inits$eta2,inits$eta3,inits$eta4)
v=inits$v
beta=inits$beta
M=c(inits$M1,inits$M2,inits$M3,inits$M4)
sgmr2=inits$sgmr2
sgm2=inits$sgm2
E=inits$E
c=c_ini
b=inits$b
e=inits$e
X=X_sim    ##rename the simulated compelete X_{it}

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

Pi_c=matrix(0, n, K)   ## for updating c 
f_c=matrix(0,n, K)     ## for updating c

PG_eta=matrix(0, n, K-1)    ## for w_{il} in updating eta_{l}
k_eta=matrix(0, n, K-1)     ## for k_{l} in updating eta_{l}
PG_v=matrix(0,n, max(T))    ## for w_{it} in updating v and e_{i}
k_v_mat=matrix(0,n,max(T))  ## for k_{v} in updating v
B_v=matrix(0,sum(T), 2)     ## for updating v 
for(i in 1:n){
  B_v[(sum(T[1:i])-T[i]+1):sum(T[1:i]),1]=rep(V_sim[i],T[i])
  B_v[(sum(T[1:i])-T[i]+1):sum(T[1:i]),2]=X[1:T[i],i]
}

for (m in 1:n_iter){  ##iteration starts
  
  ##sample c_{i} for all i
  for (i in 1:n){
    for (l in 1:K){
      Pi_c[i,l]=exp(V_sim[i]*eta[l])/sum(exp(V_sim[i]*eta))
    }
  }
  
  for(i in 1:n){
    for (l in 1:K){
      f_c[i,l]=dmvnorm(X[1:T[i],i], mean=D[1:T[i],c(2*i-1,2*i)]%*%beta+
                                    D_star[1:T[i],i]*M[l]+
                                    D_dstar[1:T[i],i]*b[i],
                                    sigma=sgm2*diag(T[i])
      )
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
  

  ##sample beta
  var_beta=solve((1/beta_pri)*diag(2)+(1/sgm2)*sum_D)
  sum_beta=rep(0, 2)
  for (i in 1:n){
    sum_beta=sum_beta+t(D[1:T[i],c(2*i-1, 2*i)])%*%(X[1:T[i], i]-D_star[1:T[i],i]*M[c[i]]-D_dstar[1:T[i],i]*b[i])
  }
  mean_beta=(1/sgm2)*var_beta%*%sum_beta
  beta=mvrnorm(n=1, mu=mean_beta, Sigma = var_beta)
  
  ##sample M2, M3 and M4
  M_mat1=matrix(0, n, K)
  M_mat2=matrix(0, n, K)
  for(l in 2:K){
    index=which(c==l)
    num_index=length(index)
    for(i in 1:num_index){
      M_mat1[i,l]=t(D_star[1:T[index[i]], index[i]])%*%D_star[1:T[index[i]], index[i]]
      M_mat2[i,l]=t(D_star[1:T[index[i]], index[i]])%*%(X[1:T[index[i]],index[i]]-
                                                          D[1:T[index[i]], c(2*index[i]-1, 2*index[i])]%*%beta-
                                                          b[index[i]]*D_dstar[1:T[index[i]], index[i]])
    }
    var_M=((1/M_pri)+(1/sgm2)*colSums(M_mat1)[l])^{-1}
    mean_M=(1/sgm2)*var_M*colSums(M_mat2)[l]
    M[l]=rnorm(1, mean=mean_M, sd=sqrt(var_M))
  }
  
  ##sample sigma^{2}
  sum_sgm2=0
  for (i in 1:n){
    sum_sgm2=sum_sgm2+t(X[1:T[i], i]-D[1:T[i], c(2*i-1, 2*i)]%*%beta-D_star[1:T[i],i]*M[c[i]]-D_dstar[1:T[i],i]*b[i])%*%
      (X[1:T[i], i]-D[1:T[i], c(2*i-1, 2*i)]%*%beta-D_star[1:T[i],i]*M[c[i]]-D_dstar[1:T[i],i]*b[i])
  }
  rate_sgm2=(1/2)*sum_sgm2+sgm2_pri
  sgm2=rigamma(n=1, alpha=(1/2)*sum(T)+sgm2_pri, beta=rate_sgm2)
  
  ##sample sigma_{r}^{2}
  sgmr2=rigamma(n=1, alpha=(n/2)+sgmr2_pri, beta=(1/2)*sum(b^{2})+sgmr2_pri) 
  
  ##sample b[i]
  for (i in 1:n){
    var_b=((1/sgmr2)+(1/sgm2)*t(D_dstar[1:T[i],i])%*%D_dstar[1:T[i],i])^{-1}
    mean_b=(1/sgm2)*var_b*t(D_dstar[1:T[i],i])%*%(X[1:T[i],i]-D[1:T[i],c(2*i-1, 2*i)]%*%beta-D_star[1:T[i],i]*M[c[i]])
    b[i]=rnorm(1, mean=mean_b, sd=sqrt(var_b))
  }  

  ##sample E
  E=rigamma(n=1, alpha=n/2+E_pri, beta=(1/2)*sum(e^{2})+E_pri)
  
  ##sample eta2, eta3 and eta4
  for(l in 2:K){
    for(i in 1:n){
      tilting_eta=V_sim[i]*eta[l]-log(sum(exp(V_sim[i]*eta[-l])) )
      PG_eta[i,l-1]=rpg(num=1, h=1, z=tilting_eta)
      index_eta=ifelse(c[i]==l,1,0)
      k_eta[i,l-1]=index_eta-1/2+PG_eta[i,l-1]*log(sum(exp(V_sim[i]*eta[-l])))
    }
    S_eta=(1/eta_pri+t(V_sim)%*%diag(PG_eta[,l-1])%*%V_sim)^{-1}
    m_eta=S_eta*t(V_sim)%*%k_eta[,l-1]
    eta[l]=rnorm(n=1,mean=m_eta, sd=sqrt(S_eta))
  }
  
  ##sample v
  for(i in 1:n){
    for (t in 1:T[i]){
      PG_v[i, t]=rpg(num=1,h=1,z=V_sim[i]*v[1]+X[t,i]*v[2]+e[i])   ##sample w_{it}^{*}
      k_v_mat[i,t]=Y_sim[t,i]-1/2-PG_v[i,t]*e[i]
    }
  }
  omega_v=NULL
  for(i in 1:n){
    omega_v=c(omega_v,PG_v[i,1:T[i]])
  }
  k_v=NULL
  for(i in 1:n){
    k_v=c(k_v,k_v_mat[i,1:T[i]])
  }
  S_v=solve((1/v_pri)*diag(2)+t(B_v)%*%diag(omega_v)%*%B_v)
  m_v=S_v%*%t(B_v)%*%k_v
  v=mvrnorm(n=1, mu=m_v, Sigma = S_v)
  
  ##sample e_{i}
  for(i in 1:n){
    S_e=((1/E)+ sum(PG_v[i,1:T[i]]) )^{-1}
    B_k_e=sum(Y_sim[1:T[i],i])-T[i]*(1/2)-sum(PG_v[i,1:T[i]]*(B_v[(sum(T[1:i])-T[i]+1):sum(T[1:i]), ]%*%v))  ##B_{i}^{\top}k_{e_{i}}
    m_e=S_e*B_k_e
    e[i]=rnorm(n=1, mean=m_e, sd=sqrt(S_e))
  }
  
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




