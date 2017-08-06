#--------------------------
# traceplots, no burn-in
#--------------------------
jpeg(file=paste(path, 'eta2.png', sep="/"))
traceplot(x=as.mcmc(eta_keep[,2]), ylab="eta_2")
dev.off()

jpeg(file=paste(path, 'eta3.png', sep="/"))
traceplot(x=as.mcmc(eta_keep[,3]), ylab="eta_3")
dev.off()

jpeg(file=paste(path, 'eta4.png', sep="/"))
traceplot(x=as.mcmc(eta_keep[,4]), ylab="eta_4")
dev.off()

jpeg(file=paste(path, 'v1.png', sep="/"))
traceplot(x=as.mcmc(v_keep[,1]), ylab="v_1")
dev.off()

jpeg(file=paste(path, 'v2.png', sep="/"))
traceplot(x=as.mcmc(v_keep[,2]), ylab="v_2")
dev.off()

jpeg(file=paste(path, 'beta1.png', sep="/"))
traceplot(x=as.mcmc(beta_keep[,1]), ylab="beta_1")
dev.off()

jpeg(file=paste(path, 'beta2.png', sep="/"))
traceplot(x=as.mcmc(beta_keep[,2]), ylab="beta_2")
dev.off()

jpeg(file=paste(path, 'M2.png', sep="/"))
traceplot(x=as.mcmc(M_keep[,2]), ylab="M_2")
dev.off()

jpeg(file=paste(path, 'M3.png', sep="/"))
traceplot(x=as.mcmc(M_keep[,3]), ylab="M_3")
dev.off()

jpeg(file=paste(path, 'M4.png', sep="/"))
traceplot(x=as.mcmc(M_keep[,4]), ylab="M_4")
dev.off()

jpeg(file=paste(path, 'sgmr2.png', sep="/"))
traceplot(x=as.mcmc(sgmr2_keep), ylab="sgmr2")
dev.off()

jpeg(file=paste(path, 'sgm2.png', sep="/"))
traceplot(x=as.mcmc(sgm2_keep), ylab="sgm2")
dev.off()

jpeg(file=paste(path, 'E.png', sep="/"))
traceplot(x=as.mcmc(E_keep), ylab="E")
dev.off()

