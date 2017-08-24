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
datfrm=data.frame(Subject=rep(1:n, T), V=rep(V_sim,T), X=X_sim, Y=Y_sim, R=R_sim) 
X_sim[R_sim==0]=NA

# check if observations of X_{i} are all missing
#checkNA=vector(length=n)
#for (i in 1:n){
#  pos=which(datfrm$Subject==i)
#  checkNA[i]=all(R_sim[pos]==0)
#}
#which(checkNA==TRUE)








