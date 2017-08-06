#-----------------------------------
# mean plots after burn-in
#-----------------------------------
source("/Users/peterwang/Dropbox/missingdata/Project/code/functions.R")
jpeg(file=paste(path, 'eta2.png', sep="/"))
runmeanplot(eta_keep[-(1:burnin),2])
dev.off()

jpeg(file=paste(path, 'eta3.png', sep="/"))
runmeanplot(eta_keep[-(1:burnin),3])
dev.off()

jpeg(file=paste(path, 'eta4.png', sep="/"))
runmeanplot(eta_keep[-(1:burnin),4])
dev.off()

jpeg(file=paste(path, 'v1.png', sep="/"))
runmeanplot(v_keep[-(1:burnin),1])
dev.off()

jpeg(file=paste(path, 'v2.png', sep="/"))
runmeanplot(v_keep[-(1:burnin),2])
dev.off()

jpeg(file=paste(path, 'beta1.png', sep="/"))
runmeanplot(beta_keep[-(1:burnin),1])
dev.off()

jpeg(file=paste(path, 'beta2.png', sep="/"))
runmeanplot(beta_keep[-(1:burnin),2])
dev.off()

jpeg(file=paste(path, 'M2.png', sep="/"))
runmeanplot(M_keep[-(1:burnin),2])
dev.off()

jpeg(file=paste(path, 'M3.png', sep="/"))
runmeanplot(M_keep[-(1:burnin),3])
dev.off()

jpeg(file=paste(path, 'M4.png', sep="/"))
runmeanplot(M_keep[-(1:burnin),4])
dev.off()

jpeg(file=paste(path, 'sgmr2.png', sep="/"))
runmeanplot(sgmr2_keep[-(1:burnin)])
dev.off()

jpeg(file=paste(path, 'sgm2.png', sep="/"))
runmeanplot(sgm2_keep[-(1:burnin)])
dev.off()

jpeg(file=paste(path, 'E.png', sep="/"))
runmeanplot(E_keep[-(1:burnin)])
dev.off()

