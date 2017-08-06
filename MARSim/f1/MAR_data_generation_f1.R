#----------------------------------------------------------------------
# This script is used to generate the simulated data in  MAR setting
# and number of subjects is 500
# Last updated date: 8/3/2017
#----------------------------------------------------------------------
source("/Users/peterwang/Dropbox/missingdata/Project/code/cpl_data_generation.R")

#--------------------------------
# SET MISSING VALUES: MAR    
#--------------------------------
alpha_sim=c(-0.5,1.5)
BigVY_sim=matrix(0, sum(T), 2)
BigVY_sim[,1]=rep(V_sim,T)
BigVY_sim[,2]=Y_sim
prob_sim=exp(BigVY_sim%*%alpha_sim)/(1+exp(BigVY_sim%*%alpha_sim))
R_sim=rbinom(sum(T), 1, prob=prob_sim)  ##R=0 represents missing
#create a data frame
datfrm=data.frame(Subject=rep(1:n, T), V=rep(V_sim,T), X=X_sim, Y=Y_sim, R=R_sim) 
X_sim[R_sim==0]=NA

# check if observations of X_{i} are all missing
#checkNA=vector(length=n)
#for (i in 1:n){
#  pos=which(datfrm$Subject==i)
#  checkNA[i]=all(R_sim[pos]==0)
#}
#which(checkNA==TRUE)

