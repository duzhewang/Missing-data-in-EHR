#----------------------------------------------------------------------
# This script is used to generate the simulated data in MAR setting
# Last updated date: 8/20/2017
#----------------------------------------------------------------------
source("/Users/dwang282/Desktop/f4/cpl_data_generation.R")

#--------------------------
# SET MISSING VALUES: MAR    
#--------------------------
# we assume older patients may be less likely to have A1C measured 
# while female patients may be more likely to have A1C measured. 
# age is a categorical variable with 6 groups: <=39, 40-49, 50-59, 60-69, 70-79 and >=80. 

age_sim=sample(20:90,size = n, replace = TRUE)
gender_sim=rbinom(n, size=1, p=0.5) ##gender=1 represents female
age_coef=c(-0.2, -0.1, 0, 0.2, 0.3, 0.4)
gender_coef=-0.1
alpha=0.1 ##the intercept is set to 0.1 in order to have around 50% missing values. 

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
datfrm=data.frame(Subject=rep(1:n, T), V=rep(V_sim,T), X=X_sim, Y=Y_sim, R=R_sim) 
X_sim[R_sim==0]=NA

# check if observations of X_{i} are all missing
#checkNA=vector(length=n)
#for (i in 1:n){
# pos=which(datfrm$Subject==i)
#  checkNA[i]=all(R_sim[pos]==0)
#}
#which(checkNA==TRUE)

