# Missing-data-in-EHR

For privacy, The detailed design of this simulation is described, updated weekly and uploaded to Box.com. 

-**Some important variable names**:
1. v_sim: synthetic time-invariant covariate
2. eta_sim: true value of eta
3. beta_sim: true value of beta
4. M_sim: true value of M
5. v_sim: true value of v
6. c_sim:true latent class for each subject
7. b_sim: random effect used to generate the original data
8. X_sim: synthetic data 
9. Y_sim: synthetic data
10. XM_sim: synthetic data with missing values

-**Variable names used in the Gibbs sampling iteration**:
    XM, eta, v, beta, M, sgmr2, sgm2, E, b, c, e

-**Updates**:

1. As of 6/14/2017: testing Case 1 listed in Table 1 of Report 3.
2. 6/15/2017: test issue: codes are unstable, sometimes work, sometimes doesn't. 



