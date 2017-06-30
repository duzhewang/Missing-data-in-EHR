#----------------------------------
# BAYESIAN METHOD IN MISSING DATA  
# Last updated date: 6/26/2017
#----------------------------------
rm(list=ls())
library(MASS) ##use function mvrnorm
library(LaplacesDemon) ##use function rcat
library(MCMCpack) ##use function rinvgamma, default rate=1
library(mvtnorm)  ##use function dmvnorm, rmvt
library(pscl)  ##use function rigamma 
library(metRology)  ##use function rt.scaled 
library(Rmpfr)   ##use mpfr


#---------------------------------
#      Global set-up       
#---------------------------------
n=100 ##number of subject
K=4 ##number of latent classes

##assign each subject the number of tracked quarters
T=sample(5:40, size=100, replace=TRUE) 

#--------------------------------
#   Complete data generation
#--------------------------------
##set up true values of eta, beta, M, v, sgm2 and sgmr2
eta_sim=c(0, 0.7, 1.6, 0.9)
beta_sim=c(-0.4, 0.5)
M_sim=c(-1.6, -0.6, 0.6, 1.2) 
v_sim=c(0.5, -0.3)
sgmr2_sim=1                ##true value of variance of b_{i}
sgm2_sim=1                 ##true value of variance of epsilon_{it}
E_sim=1                    ##true value of variance of e_{i}

##generate time-invariate covariate V_{i}
V_sim=rnorm(n=100, mean=0, sd=1)

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
D=matrix(0, max(T), 2*n)
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
X_sim=matrix(0, nrow=max(T), ncol = n) ##each column is one subject
for (i in 1:n){
  X_sim[1:T[i],i]=D[1:T[i],c(2*i-1,2*i)]%*%beta_sim+
                  D_star[1:T[i],i]*M_sim[c_sim[i]]+
                  D_dstar[1:T[i],i]*b_sim[i]+          ##each column is one subject
                  rnorm(n=T[i], mean = 0, sd=sqrt(sgm2_sim))
}

##generate random effect e_{i}
e_sim=rnorm(n, mean=0, sd=sqrt(E_sim))

##generate y_{it}
Y_sim=matrix(0,nrow=max(T), ncol=n)
for (i in 1:n){
  top= exp(cbind(D[1:T[i],2*i], X_sim[1:T[i],i])%*%v_sim+rep(e_sim[i],T[i]))
  bot=rep(1, T[i])+top
  Y_sim[1:T[i],i]=rbinom(n=T[i], size=1, p=top/bot)   ##each column is one subject
}

#-------------------------------------------
#          PRIORS AND INITIAL VALUES    
#-------------------------------------------
beta_pri=10^{4}   ##variance of normal prior
M_pri=10^{4}
sgm2_pri=0.001
sgmr2_pri=0.001
v_pri=10^{4}
E_pri=0.001
eta_pri=10^{4}


inits=list(eta1=0, eta2=0.5, eta3=1.3, eta4=1,
           beta=c(-0.5, 0.8),
           M1=-2, M2=-1, M3=1.2, M4=0.6, 
           v=c(0.3, -1),
           sgmr2=1.5^2,
           sgm2=1.5^2, 
           E=1.5^2, 
           b=rnorm(n, mean=0, sd=1.5),             ##initial value of random effect b_{i}
           e=rnorm(n, mean=0, sd=1.5)              ##initial value of random effect e_{i}
)

#set initial value of latent class for each subject
eta_ini=c(inits$eta1,inits$eta2,inits$eta3,inits$eta4)
Pi_ini=matrix(0, n, K )  
for (i in 1:n){
  for (l in 1:K){
    Pi_ini[i,l]=exp(V_sim[i]*eta_ini[l])/sum(exp(V_sim[i]*eta_ini))
  }
}
c_ini=apply(Pi_ini, 1, function(x) rcat(n=1, p=x)) 

#----------------------------------------------------------
#            GIBBS SAMPLING-ITERATIONS                       
#      order of iteration: c[i],eta,b[i],e[i],
#                          v,beta,M,sgmr2, sgm2, E  
#----------------------------------------------------------
n_iter=3           ##number of iterations

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
X=X_sim    ##rename the simulated compelete X_{it} data as X

##recording structure, each row is one iteration
eta_keep=matrix(0, nrow=n_iter, ncol=3)  
v_keep=matrix(0, nrow=n_iter, ncol=2)    
beta_keep=matrix(0, nrow=n_iter, ncol=2)  
M_keep=matrix(0, nrow=n_iter, ncol=4)    
sgmr2_keep=rep(0, n_iter)  
sgm2_keep=rep(0, n_iter)   
E_keep=rep(0, n_iter)     
c_keep=matrix(0, nrow=n_iter, ncol=n) 
b_keep=matrix(0, nrow=n_iter, ncol=n) 
e_keep=matrix(0, nrow=n_iter, ncol=n) 

Pi_c=matrix(0, n, K)   ## for updating c 
f_c=matrix(0,n, K)     ## for updating c


for (m in 1:n_iter){  ##iteration 

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
  
 # c=apply(Pi_c*f_c, 1, function(x) rcat(n=1, p=x/sum(x))) ##assign latent class to each subject

  ##sample eta: Metropolis-within-Gibbs
  for (l in 2:K){
      target=function(x){     ##define a target function
      prodtarget=1
      for (i in 1:n){
        prodtarget=ifelse(c[i]==l,prodtarget*exp(V_sim[i]*x)/(sum(exp(V_sim[i]*eta[-l]))+exp(V_sim[i]*x)),
                          prodtarget*exp(V_sim[i]*eta[c[i]])/(sum(exp(V_sim[i]*eta[-l]))+exp(V_sim[i]*x)) 
                         )
      }
      prodtarget=prodtarget*dnorm(x, mean = 0, sd=sqrt(eta_pri), log=FALSE)
      return(prodtarget)
    }
    eta_prop=rt.scaled(1, df=3, mean=eta[l], sd=1)
    ratio=(target(eta_prop)/target(eta[l]))*(dt.scaled(x=eta[l], df=3, mean=eta_prop, sd=1 )/dt.scaled(x=eta_prop, df=3, mean=eta[l],sd=1))
    eta[l]=ifelse(runif(1)<ratio, eta_prop, eta[l] )
  }
  
  ##sample b[i]
  for (i in 1:n){
    var_b=((1/sgmr2)+(1/sgm2)*t(D_dstar[1:T[i],i])%*%D_dstar[1:T[i],i])^{-1}
    mean_b=(1/sgm2)*var_b*t(D_dstar[1:T[i],i])%*%(X[1:T[i],i]-D[1:T[i],c(2*i-1, 2*i)]%*%beta-D_star[1:T[i],i]*M[c[i]])
    b[i]=rnorm(1, mean=mean_b, sd=sqrt(var_b))
  }  
  
  ##sample e[i]: Metropolis-within-Gibbs
  for(i in 1:n){
    e_target=function(x){
      e_prodtarget=1
      for (t in 1:T[i]){
        e_prodtarget=ifelse(Y_sim[t,i]==1,e_prodtarget*exp(V_sim[i]*v[1]+X[t,i]*v[2]+x)/(1+exp(V_sim[i]*v[1]+X[t,i]*v[2]+x)),
                            e_prodtarget*1/(1+exp(V_sim[i]*v[1]+X[t,i]*v[2]+x)) 
                           )
      }      
      e_prodtarget=e_prodtarget*dnorm(x, mean = 0, sd=sqrt(E), log=FALSE)
      return(e_prodtarget)
    }
    e_prop=rt.scaled(1, df=3, mean=e[i], sd=1)
    e_ratio=(e_target(e_prop)/e_target(e[i]))*(dt.scaled(x=e[i], df=3, mean=e_prop,sd=1)/dt.scaled(x=e_prop, df=3, mean=e[i],sd=1)) 
    e[i]=ifelse(runif(1)<e_ratio, e_prop, e[i] ) ## test accepting or not
  }
  
  ##sample v: Metropolis-within-Gibbs
  v_target=function(x){
    v_prodtarget=matrix(1, max(T), n)
    for(i in 1:n){
      for (t in 1:T[i]){
        v_prodtarget[t,i]=ifelse(Y_sim[t,i]==1, exp(V_sim[i]*x[1]+X[t,i]*x[2]+e[i])/(1+exp(V_sim[i]*x[1]+X[t,i]*x[2]+e[i])), 
                                 1/(1+exp(V_sim[i]*x[1]+X[t,i]*x[2]+e[i])) 
                                )
      }
    }
    v_prodtarget=prod(v_prodtarget)*dmvnorm(x, mean=rep(0,2), sigma=v_pri*diag(2))
    return(v_prodtarget)
  }
  v_prop=rmvt(1, sigma=1/3*diag(2), df=3, delta=v, type="shifted") ##propose a new value by multivariate t distribution
  v_ratio=(v_target(v_prop)/v_target(v))*(dmvt(v, delta=v_prop, sigma=1/3*diag(2), df=3, log=FALSE, type="shifted")/dmvt(v_prop, delta=v, sigma=1/3*diag(2), df=3, log=FALSE, type="shifted")) 
  if(runif(1)<v_ratio){
    v=v_prop
  } else{
    v=v
  }
  
  ##sample beta
  var_beta=solve((1/beta_pri)*diag(2)+(1/sgm2)*sum_D)
  sum_beta=rep(0, 2)
  for (i in 1:n){
    sum_beta=sum_beta+t(D[1:T[i],c(2*i-1, 2*i)])%*%(X[1:T[i], i]-D_star[1:T[i],i]*M[c[i]]-D_dstar[1:T[i],i]*b[i])
  }
  mean_beta=(1/sgm2)*var_beta%*%sum_beta
  beta=mvrnorm(n=1, mu=mean_beta, Sigma = var_beta)
  
  ##sample M
  M_mat1=matrix(0, n, K)
  M_mat2=matrix(0, n, K)
  for(l in 1:K){
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
  
  ##sample sigma_{r}^{2}
  sgmr2=rigamma(n=1, alpha=(n/2)+sgmr2_pri, beta=1/2*sum(b^{2})+sgmr2_pri) 
  
  ##sample sigma^{2}
  sum_sgm2=0
  for (i in 1:n){
    sum_sgm2=sum_sgm2+t(X[1:T[i], i]-D[1:T[i], c(2*i-1, 2*i)]%*%beta-D_star[1:T[i],i]*M[c[i]]-D_dstar[1:T[i],i]*b[i])%*%
                       (X[1:T[i], i]-D[1:T[i], c(2*i-1, 2*i)]%*%beta-D_star[1:T[i],i]*M[c[i]]-D_dstar[1:T[i],i]*b[i])
  }
  rate_sgm2=1/2*sum_sgm2+sgm2_pri
  sgm2=rigamma(n=1, alpha=1/2*sum(T)+sgm2_pri, beta=rate_sgm2)
  
  ##sample E
  E=rigamma(n=1, alpha=n/2+E_pri, beta=1/2*sum(e^{2})+E_pri)
  
  ##record parameters
  c_keep[m, ]=c
  eta_keep[m,]=eta[2:K]
  b_keep[m, ]=b
  e_keep[m, ]=e
  v_keep[m, ]=v
  beta_keep[m, ]=beta
  M_keep[m, ]=M
  sgmr2_keep[m]=sgmr2
  sgm2_keep[m]=sgm2
  E_keep[m]=E
  
} ##for m-th iteration




