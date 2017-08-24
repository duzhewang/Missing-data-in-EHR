# Missing-data-in-EHR

## About R files

| file names    | descriptions  |
| ------------- |:-------------:|
| cpl_data_generation.R | generate simulated complete data set | 
| functions.R     | a collection of functions which are used in other files     |  
| diagnostics folder | R scripts used to do MCMC diagnostics, including traceplots, autocorrelation plots, running mean plots and effective sample size     |  
| CompleteDataSim folder   | Complete data simulation includes Study 1,2 and 3  |
| MCARSim folder     | MCAR simulation includes Study 4, 5 and 6    |  
| MARSim folder    | MAR simulation includes Study 7, 8, 9 and 10   |  
| Obsfit folder     | simulation only using observed data includes Study 11-16|   



In each simulation study folder of CompleteDataSim, there are two R scripts, one is used to update MCMC and the other is used to run MCMC, save results and do diagnostics. See [here](https://github.com/dzwang91/Missing-data-in-EHR/tree/master/CompleteDataSim/f1) as an example. 

For other kinds of simulation studies(MCAR,MAR and Obsfit), there are three R scripts in each study folder. See [here](https://github.com/dzwang91/Missing-data-in-EHR/tree/master/MCARSim/f6) as an example.  **'data_generation_R'** is used to generate data for specific setting based on **'cpl_data_generation.R'**. **'mcmc_update.R'** is used to set up the initial values and sample posterior in mcmc, **'mcmc_run.R'** is used to run mcmc and do diagonistics. 

In **'cpl_data_generation.R'** file, true values can be changed here.
```
#---------------------------------
#      Global set-up.     
#---------------------------------
n=500 ##number of subject
K=4 ##number of latent classes
minQ=5   ##minimum number of tracked quarters
maxQ=45  ##maximum number of tracker quarters
eta_sim=c(0, 0.5, 1.5, 1)
beta_sim=c(-0.4, 0.5)
M_sim=c(0, -0.6, 0.6, 1.2) 
v_sim=c(0.5, -0.3)
sgmr2_sim=1                ##true value of variance of b_{i}
sgm2_sim=1                 ##true value of variance of epsilon_{it}
E_sim=1                    ##true value of variance of e_{i}
```

In **'mcmc_update.R'** files, initial values can be changed here.
```
#---------------------------
# PRIORS AND INITIAL VALUES   
#--------------------------
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

```


## About simulation studies

See [here](https://github.com/dzwang91/Missing-data-in-EHR/blob/master/simulations.pdf) for all simulation studies. 









