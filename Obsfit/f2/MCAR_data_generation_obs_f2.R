#---------------------------------------------------------------------------------
# This script is used to generate the simulated data with 50% MCAR missing values
# and number of subjects is 500
# Last updated date: 8/7/2017
#---------------------------------------------------------------------------------
source("/Users/peterwang/Dropbox/missingdata/Project/code/cpl_data_generation.R")

#--------------------------------
# SET MISSING VALUES: 50% MCAR    
#--------------------------------
mr=0.5 ##missing rate
R_sim=rbinom(sum(T), 1, 1-mr)  ##R=0 represents missing
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
# USE ONLY OBSERVED DATA
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

