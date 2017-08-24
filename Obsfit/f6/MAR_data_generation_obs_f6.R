#-----------------------------------------------------
# This script is used to generate the simulated data 
# Last updated date: 8/20/2017
#-----------------------------------------------------
source("/Users/peterwang/Dropbox/missingdata/Project/code/cpl_data_generation.R")

#--------------------------------
# SET MISSING VALUES: MAR    
#--------------------------------
age_sim=sample(20:90,size = n, replace = TRUE)
gender_sim=rbinom(n, size=1, p=0.5) ##gender=1 represents female
age_coef=c(-0.2, -0.1, 0, 0.2, 0.3, 0.4)
gender_coef=-0.1
alpha=0.1 ##the intercept is set to 0.1 in order to have around 50% missing probability for each subject. 

age_ind=matrix(0,n,6)
for (i in 1:n){
  if (age_sim[i]<=39){
    age_ind[i, ]=c(1,0,0,0,0,0)
  } else if ( age_sim[i]>39 & age_sim[i]<=49){
    age_ind[i, ]=c(0,1,0,0,0,0)
  } else if ( age_sim[i]>49 & age_sim[i]<=59){
    age_ind[i, ]=c(0,0,1,0,0,0)
  } else if (age_sim[i]>59 & age_sim[i]<=69){
    age_ind[i, ]=c(0,0,0,1,0,0)
  } else if (age_sim[i]>69 & age_sim[i]<=79){
    age_ind[i, ]=c(0,0,0,0,1,0)
  } else
    age_ind[i, ]=c(0,0,0,0,0,1)
}

##probability of not missing 
prob_sim=1/(1+exp(alpha+age_ind%*%age_coef+gender_sim*gender_coef))
R_sim=rbinom(sum(T), 1, prob=rep(prob_sim, T))  ##R=0 represents missing
#create a data frame
datfrm=data.frame(Subject=rep(1:n, T), V=rep(V_sim,T), 
                  X=X_sim, Y=Y_sim, D=D, D_star=D_star, D_dstar=D_dstar, R=R_sim) 

# check if observations of X_{i} are all missing
#checkNA=vector(length=n)
#for (i in 1:n){
#  pos=which(datfrm$Subject==i)
#  checkNA[i]=all(R_sim[pos]==0)
#}
#which(checkNA==TRUE)


#--------------------------
# USE OBSERVED DATA ONLY
#--------------------------
#new data frame consists of observed data
datfrm=datfrm[R_sim==1,]
#observed X and corresponding Y
X_sim=datfrm$X
Y_sim=datfrm$Y
#new D, D_star and D_dstar only using observed data
D=as.matrix(cbind(datfrm$D.1, datfrm$D.2))
D_star=datfrm$D_star
D_dstar=datfrm$D_dstar
#new T
for(i in 1:n){
  T[i]=length(which(datfrm$Subject==i)) 
}



