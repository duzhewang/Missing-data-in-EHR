##########################################
# BAYESIAN METHOD IN MISSING DATA  
# CURRENT VERSION: 6/14/2017
##########################################
rm(list=ls())
library(MASS) ##use function mvrnorm
library(LaplacesDemon) ##use function rcat
library(MCMCpack) ##use function rinvgamma, default rate=1
library(mvtnorm)  ##use function dmvnorm, rmvt
library(pscl)  ##use function rigamma 
library(metRology)  ##use function rt.scaled 

###########################################
##        Part 1: global set-up          ##
###########################################
n=100 ##number of subject
T=5 ##number of tracked quarters
K=4 ##number of latent classes


#############################################
##Part 2: complete original data generation##
#############################################

##generate time-invariate covariate V_{i} 
V_sim=rep(0, n)
for (i in 1:100){
  V_sim[i]=rnorm(1, mean=0,sd=1)
}

##generate true vaules of eta for sub-model 1 
eta_sim=c(0, rnorm(K-1, mean=0, sd=3/2))


##generate true values of beta for sub-model 2 
beta_sim=mvrnorm(n=1, mu=rep(0, 2), Sigma=diag(2))

##generate M1, M2, M3, M4 for sub-model 2 
M_sim=c(0, rnorm(K-1, mean=0, sd=1))

##generate v for sub-model 3 
v_sim=mvrnorm(n=1, mu=rep(0, 2), Sigma=diag(2))

##generate latent class c[i] for each subject 
Pi_sim=matrix(0, n, K )   ##Pi[i,l]=P(c_{i}=l)
for (i in 1:n){
  for (l in 1:K){
    Pi_sim[i,l]=exp(V_sim[i]*eta_sim[l])/sum(exp(V_sim[i]*eta_sim))
  }
}
c_sim=apply(Pi_sim,1, function(x) rcat(n=1, p=x)) ##assign latent class to each subject

##generate x[it] 
b_sim=rnorm(n,mean=0, sd=1) ##generate random effect b[i]
X_sim=matrix(0, nrow=n, ncol = T)
for (i in 1:n){
  for (t in 1:T){
    epsilon=rnorm(n=1, mean=0, sd=1)
    X_sim[i,t]=V_sim[i]*beta_sim[1]+t*beta_sim[2]+M_sim[c_sim[i]]+b_sim[i]+epsilon
  }
}

##generate binary random variable y[it] 
Prob_sim=matrix(0,n,T)
Y_sim=matrix(0,n,T)
e_sim=rnorm(n,mean=0, sd=1) ##random efect e[i] in sub-model 3
for(i in 1:n){
  for (t in 1:T){
    Prob_sim[i,t]=exp(V_sim[i]*v_sim[1]+X_sim[i,t]*v_sim[2]+e_sim[i])/(1+exp(V_sim[i]*v_sim[1]+X_sim[i,t]*v_sim[2]+e_sim[i]))
    Y_sim[i,t]=rbinom(1,1,Prob_sim[i,t])
  }
}


######################################################################
##    visualization of the complete original data           
#df=data.frame(V_sim, X_sim, Y_sim)
#colnames(df)=c("V","X1","X2","X3","X4","X5","Y1","Y2","Y3","Y4","Y5")
#head(df)
######################################################################


#############################################
##  Part 3: Missing data generation-MCAR   ##
#############################################

##CASE1: missing rate=10% at all time points
mr=0.1  #missing rate
R_sim=matrix(rbinom(n*T,1,1-mr), n, T)  ##R[i,t]=0 represents missing
XM_sim=X_sim
XM_sim[R_sim==0]=NA ##missing data represented by NA


########################################################
##   Part 4: GIBBS SAMPLING-SET INITIAL VALUES        ##
########################################################

XM=XM_sim   #copy the simulated missing data set x[it] 
MI=data.matrix(which(R_sim==0, arr.ind = T),rownames.force = NA) ##missing index, which row and which column

##set up initial values of parameters from their prior distributions
inits=list(eta1=0, eta2=rnorm(1, mean=0, sd=3/2), eta3=rnorm(1, mean=0, sd=3/2), eta4=rnorm(1, mean=0, sd=3/2),
           v=c(rnorm(1, mean=0, sd=1),rnorm(1, mean=0, sd=1)), beta=c(rnorm(1, mean=0, sd=1),rnorm(1, mean=0, sd=1)),
           M1=0, M2=rnorm(1, mean=0, sd=1), M3=rnorm(1, mean=0, sd=1), M4=rnorm(1, mean=0, sd=1), 
           sgmr2=rinvgamma(1,shape = 1), sgm2=rinvgamma(1,shape = 1), E=rinvgamma(1,shape = 1)
) 

##set up initial values of missing data x[it]: first generate c[i] using initial values of eta by model (1), and then
##generate initial values of missing x[it] by model (2). 
eta_ini=c(inits$eta1, inits$eta2, inits$eta3, inits$eta4)  
beta_1_ini=inits$beta[1]
beta_2_ini=inits$beta[2]
M_ini=c(inits$M1, inits$M2, inits$M3, inits$M4)
sgmr_ini=sqrt(inits$sgmr2) 
sgm_ini=sqrt(inits$sgm2)

Pi_ini=matrix(0, n, K )   
for (i in 1:n){
  for (l in 1:K){
    Pi_ini[i,l]=exp(V_sim[i]*eta_ini[l])/sum(exp(V_sim[i]*eta_ini))
  }
}
c_ini=apply(Pi_ini,1, function(x) rcat(n=1, p=x)) ##assign initial latent class to each subject

b_ini=rnorm(n,mean=0, sd=sgmr_ini) ##generate initial random effect b[i] using initial sigma_{r}

##generate the initial missing x[i,t]
for (i in 1:dim(MI)[1]){
  XM[MI[i,1], MI[i,2]]=V_sim[MI[i,1]]*beta_1_ini+MI[i,2]*beta_2_ini+M_ini[c_ini[MI[i,1]]]+b_ini[MI[i,1]]
  +rnorm(1, mean=0, sd=sgm_ini)
}


##generate the initial e[i] using the initial E
e_ini=rnorm(n, mean=0, sd=sqrt(inits$E))

## create D_{i}, D_{i}^{*}, D_{i}^{**} in equation (17)
D=vector("list", n)
for ( i in 1:n){
  D[[i]]=cbind(rep(V_sim[i], T), 1:T)
  
}

D_star=vector("list",n)
for ( i in 1:n){
  D_star[[i]]=rep(1, T)
}

D_dstar=vector("list",n)
for ( i in 1:n){
  D_dstar[[i]]=rep(1, T)
}

sum_D=matrix(0, 2, 2)  
for (i in 1:n){
<<<<<<< HEAD
  sum_D=sum_D+ t(D[[i]])%*%D[[i]]   
=======
  sum_D=sum_D+ t(D[[i]])%*%D[[i]]   ##sum of D[i] to be used for updating beta 
>>>>>>> update case 1
}

####################################################################
####       Part 5: GIBBS SAMPLING-ITERATIONS                       
##order of iteration: c[i],eta,b[i],e[i],v,beta,M,sigma_{r}^{2}
##                    sigma^{2}, E, missing x[it]
####################################################################
n_iter=3 ##number of iterations
eta=c(inits$eta1,inits$eta2,inits$eta3,inits$eta4)
v=inits$v
beta=inits$beta
M=c(inits$M1,inits$M2,inits$M3,inits$M4)
sgmr2=inits$sgmr2
sgm2=inits$sgm2
E=inits$E
b=b_ini
c=c_ini
e=e_ini

##recording structure
eta_keep=matrix(0, nrow=n_iter, ncol=3)  ##record eta2, eta3, eta4
v_keep=matrix(0, nrow=n_iter, ncol=2)    ##record v of size 2
beta_keep=matrix(0, nrow=n_iter, ncol=2)  ##record beta of size 2
M_keep=matrix(0, nrow=n_iter, ncol=3)    ##record M2, M3, M4
sgmr2_keep=rep(0, n_iter)  ##record sigma_{r}^{2}
sgm2_keep=rep(0, n_iter)   ##record sigma^{2}
E_keep=rep(0, n_iter)     ##record E
c_keep=matrix(0, nrow=n_iter, ncol=n) ##record latent class for each subject, one row one iteration
b_keep=matrix(0, nrow=n_iter, ncol=n) ##record b[i], one row  one iteration
e_keep=matrix(0, nrow=n_iter, ncol=n) ##record e[i], one row one iteration
Pi1_keep=matrix(0, n, K )    ##equation (1)
Pi_c_keep=matrix(0, n, K)    ##used to update c


for (m in 1:n_iter){
  
  ##sample c[i] for all i
  for (i in 1:n){
    for (l in 1:K){
      Pi1_keep[i,l]=exp(V_sim[i]*eta[l])/sum(exp(V_sim[i]*eta)) ##equation (1)
    }
  }
  for (i in 1:n){
    mid=sapply(2:K, function(l) dmvnorm(XM[i,], mean=D[[i]]%*%beta+D_star[[i]]*M[l]+D_dstar[[i]]*b[i],sigma=sgm2*diag(T)))
    bot=Pi1_keep[i,1]*dmvnorm(XM[i,],mean=D[[i]]%*%beta+D_dstar[[i]]*b[i],sigma=sgm2*diag(T))+sum(Pi1_keep[i,-1]*mid) 
    for (l in 1:K){
      top=ifelse(l==1,Pi1_keep[i,l]*dmvnorm(XM[i,], mean=D[[i]]%*%beta+D_dstar[[i]]*b[i],sigma=sgm2*diag(T)),
                 Pi1_keep[i,l]*dmvnorm(XM[i,], mean=D[[i]]%*%beta+D_star[[i]]*M[l]+D_dstar[[i]]*b[i],sigma=sgm2*diag(T)))
      Pi_c_keep[i,l]=top/bot            
    }
  }
  c=apply(Pi_c_keep, 1, function(x) rcat(n=1, p=x)) ##assign latent class to each subject
  
  ##sample eta: Metropolis-within-Gibbs
  for (l in 2:K){
    ##define a target function
    target=function(x){
      prodtarget=1
      for (i in 1:n){
        prodtarget=ifelse(c[i]==l,prodtarget*exp(V_sim[i]*x)/( sum(exp(V_sim[i]*eta[-l]))+exp(V_sim[i]*x) ),
                          prodtarget*exp(V_sim[i]*eta[c[i]])/(sum(exp(V_sim[i]*eta[-l]))+exp(V_sim[i]*x) ) )
      }
      prodtarget=prodtarget*dnorm(x, mean = 0, sd=3/2, log=FALSE)
      return(prodtarget)
    }
    
    eta_prop=rt.scaled(1, df=3, mean=eta[l], sd=1)
    ratio=(target(eta_prop)/target(eta[l]))*(dt.scaled(x=eta[l], df=3, mean=eta_prop, sd=1 )/dt.scaled(x=eta_prop, df=3, mean=eta[l],sd=1))
    eta[l]=ifelse(runif(1)<ratio, eta_prop, eta[l] )
  }
  
  ##sample b[i]
  for (i in 1:n){
    var_b=((1/sgmr2)+(1/sgm2)*t(D_dstar[[i]])%*%D_dstar[[i]])^{-1}
    mean_b=(1/sgm2)*var_b*t(D_dstar[[i]])%*%(XM[i,]-D[[i]]%*%beta-D_star[[i]]*M[c[i]])
    b[i]=rnorm(1, mean=mean_b, sd=sqrt(var_b))
  }  
  
  ##sample e[i]: Metropolis-within-Gibbs
  for(i in 1:n){
    e_target=function(x){
      e_prodtarget=1
      for (t in 1:T){
        e_prodtarget=ifelse(Y_sim[i,t]==1,e_prodtarget*exp(V_sim[i]*v[1]+XM[i,t]*v[2]+x )/(1+exp(V_sim[i]*v[1]+XM[i,t]*v[2]+x )),
                            e_prodtarget* 1/(1+exp(V_sim[i]*v[1]+XM[i,t]*v[2]+x)) )
      }
      e_prodtarget=e_prodtarget*dnorm(x, mean = 0, sd=sqrt(E), log=FALSE)
      return(e_prodtarget)
    }
    
    e_prop=rt.scaled(1, df=3, mean=e[i], sd=1) ##propose a new value
    e_ratio=(e_target(e_prop)/e_target(e[i]))*(dt.scaled(x=e[i], df=3, mean=e_prop,sd=1)/dt.scaled(x=e_prop, df=3, mean=e[i],sd=1)) 
    e[i]=ifelse(runif(1)<e_ratio, e_prop, e[i] ) ## test accepting or not
  }
  
  ##sample v: Metropolis-within-Gibbs
  v_target=function(x){
    v_prodtarget=matrix(1, n, T)
    for(i in 1:n){
      for (t in 1:T){
        v_prodtarget[i,t]=ifelse(Y_sim[i,t]==1, exp(V_sim[i]*x[1]+XM[i,t]*x[2]+e[i])/(1+exp(V_sim[i]*x[1]+XM[i,t]*x[2]+e[i])), 
                                 1/(1+exp(V_sim[i]*x[1]+XM[i,t]*x[2]+e[i])) )
      }
    }
    v_prodtarget=prod(v_prodtarget)*dmvnorm(x, mean=rep(0,2), sigma=diag(2))
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
  var_beta=solve(diag(2)+(1/sgm2)*sum_D)
  sum_beta=rep(0, 2)
  for (i in 1:n){
    sum_beta=sum_beta+t(D[[i]])%*%(XM[i]-D_star[[i]]*M[c[i]]-D_dstar[[i]]*b[i])
  }
  mean_beta=(1/sgm2)*var_beta%*%sum_beta
  beta=mvrnorm(n=1, mu=mean_beta, Sigma = var_beta)
  
  ##sample M
  sum_M=rep(0,3)
  for(l in 2:K){
    for(i in 1:length(which(c==l))){
      sum_M[l-1]=sum_M[l-1]+sum(XM[which(c==l)[i]]-D[[which(c==l)[i]]]%*%beta-D_dstar[[which(c==l)[i]]]*b[which(c==l)[i]] )  
    }
  }
  for (l in 2:K){
    var_M=(1+(1/sgm2)*sum(c==l)*5)^{-1}
    mean_M=(1/sgm2)*var_M*(sum_M[l-1])
    M[l]=rnorm(1, mean=mean_M, sd=sqrt(var_M)) 
  }
  
  ##sample sigma_{r}^{2}
  sgmr2=rigamma(n=1, alpha=(n/2)+1, beta=1/2*sum(b^{2})+1) 
  
  ##sample sigma^{2}
  sum_rate_sgm2=0
  for (i in 1:n){
    sum_rate_sgm2=sum_rate_sgm2+sum((XM[i]-D[[i]]%*%beta-D_star[[i]]*M[c[i]]-D_dstar[[i]]*b[i])^{2})  
  }
  rate_sgm2=1/2*sum_rate_sgm2+1
  sgm2=rigamma(n=1, alpha=n*T/2+1, beta = rate_sgm2)
  
  ##sample E
  E=rigamma(n=1, alpha=n/2+1, beta=1/2*sum(e^{2})+1)
  
  ##impute missing x[it]
  for (i in 1:dim(MI)[1]){
    XM[MI[i,1], MI[i,2]]=V_sim[MI[i,1]]*beta[1]+MI[i,2]*beta[2]+M[c[MI[i,1]]]+b[MI[i,1]]
    +rnorm(1, mean=0, sd=sqrt(sgm2))
  }
  
  ##record parameters
  c_keep[m, ]=c
  eta_keep[m,]=eta[2:K]
  b_keep[m, ]=b
  e_keep[m, ]=e
  v_keep[m, ]=v
  beta_keep[m, ]=beta
  M_keep[m, ]=M[2:K]
  sgmr2_keep[m]=sgmr2
  sgm2_keep[m]=sgm2
  E_keep[m]=E
  
} ##for m-th iteration








