#------------------------
# effective sample size
#------------------------
effectiveSize(as.mcmc(eta_keep[-(1:burnin),]))
effectiveSize(as.mcmc(v_keep[-(1:burnin),]))
effectiveSize(as.mcmc(beta_keep[-(1:burnin),]))
effectiveSize(as.mcmc(M_keep[-(1:burnin),]))
effectiveSize(as.mcmc(sgmr2_keep[-(1:burnin)]))
effectiveSize(as.mcmc(sgm2_keep[-(1:burnin)]))
effectiveSize(as.mcmc(E_keep[-(1:burnin)]))
