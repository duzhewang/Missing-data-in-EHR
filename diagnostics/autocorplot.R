#--------------------------------------
# autocorrelation plots after burn-in
#--------------------------------------
jpeg(file=paste(path, 'eta2.png', sep="/"))
autocorr.plot(x=as.mcmc(eta_keep[-(1:burnin),2]), main="eta_2")
dev.off()

jpeg(file=paste(path, 'eta3.png', sep="/"))
autocorr.plot(x=as.mcmc(eta_keep[-(1:burnin),3]), main="eta_3")
dev.off()

jpeg(file=paste(path, 'eta4.png', sep="/"))
autocorr.plot(x=as.mcmc(eta_keep[-(1:burnin),4]), main="eta_4")
dev.off()

jpeg(file=paste(path, 'v1.png', sep="/"))
autocorr.plot(x=as.mcmc(v_keep[-(1:burnin),1]), main="v_1")
dev.off()

jpeg(file=paste(path, 'v2.png', sep="/"))
autocorr.plot(x=as.mcmc(v_keep[-(1:burnin),2]), main="v_2")
dev.off()

jpeg(file=paste(path, 'beta1.png', sep="/"))
autocorr.plot(x=as.mcmc(beta_keep[-(1:burnin),1]), main="beta_1")
dev.off()

jpeg(file=paste(path, 'beta2.png', sep="/"))
autocorr.plot(x=as.mcmc(beta_keep[-(1:burnin),2]), main="beta_2")
dev.off()

jpeg(file=paste(path, 'M2.png', sep="/"))
autocorr.plot(x=as.mcmc(M_keep[-(1:burnin),2]), main="M_2")
dev.off()

jpeg(file=paste(path, 'M3.png', sep="/"))
autocorr.plot(x=as.mcmc(M_keep[-(1:burnin),3]), main="M_3")
dev.off()

jpeg(file=paste(path, 'M4.png', sep="/"))
autocorr.plot(x=as.mcmc(M_keep[-(1:burnin),4]), main="M_4")
dev.off()

jpeg(file=paste(path, 'sgmr2.png', sep="/"))
autocorr.plot(x=as.mcmc(sgmr2_keep[-(1:burnin)]), main="sgmr2")
dev.off()

jpeg(file=paste(path, 'sgm2.png', sep="/"))
autocorr.plot(x=as.mcmc(sgm2_keep[-(1:burnin)]), main="sgm2")
dev.off()

jpeg(file=paste(path, 'E.png', sep="/"))
autocorr.plot(x=as.mcmc(E_keep[-(1:burnin)]), main="E")
dev.off()







