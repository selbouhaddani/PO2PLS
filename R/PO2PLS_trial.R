f <- function(none){
#set.seed(13254)
#library(OmicsPLS)
#library(parallel)

#simu_test<-parSapply(cl, 1:4, function(niks){
p = 7
q = 6
r = 2
rx = 3
ry = 1
parms = generate_params(matrix(0,1,p),matrix(0,1,q),r,rx,ry, type = 'r')
# parms$Wo = 0*parms$Wo
# parms$Co = 0*parms$Co
# parms$SigTo = 0*parms$SigTo
# parms$SigUo = 0*parms$SigUo
Dat = generate_data(1e3, parms)
X = scale(Dat[,1:p], scale = F)
Y = scale(Dat[,-(1:p)], scale = F)
#parms

#fitppls <- (PPLS::PPLS_simult(X,Y,r,1e4,1e-6))
fit_o2m = o2m(X,Y,r,rx,ry,stripped=T)
fit <- PO2PLS(X, Y, r, rx, ry, 5e4, tol = 1e-6, 'o2m')
#refit_i <- 0
#while(min(diff(fit$logl)) < -1e-5){
#  refit_i <- refit_i + 1
#  fit <- PO2PLS(X, Y, r, rx, ry, 1e3, tol = 1e-6, 'random')
#  print(refit_i)
#}
#plot(diff(fit$logl), type='l'); abline(h=0)

cmax <- function(A,B){
  apply(Reduce(crossprod, list(A,B)), 2, function(e) max(abs(e)))
}
# max_po2m <- cmax(fit$par$W, parms$W)
# max_o2m <- cmax(fit_o2m$W., parms$W)
Wmax <- cmax(fit$par$W, parms$W) - cmax(fit_o2m$W., parms$W)
Cmax <- cmax(fit$par$C, parms$C) - cmax(fit_o2m$C., parms$C)
Womax <- cmax(fit$par$Wo, parms$Wo) - cmax(orth(fit_o2m$P_Y), parms$Wo)
Comax <- cmax(fit$par$Co, parms$Co) - cmax(orth(fit_o2m$P_X), parms$Co)
#cat("Is PO2PLS better than O2PLS? \n")
list(W = Wmax, C = Cmax, Wo=Womax, Co=Comax, Negs = any(diff(fit$logl)<0))
#return(max_po2m - max_o2m)
#})
#print(rowMeans(simu_test))
}

library(parallel)
library(tidyverse)
library(OmicsPLS)
source('R/PO2PLS_functions.R')
system.time(outp <- parallelsugar::mclapply(mc.cores=4, 1:200, f))
outp2 <- cbind(data.frame(DifComp=t(sapply(outp, function(e) e$W))), data.frame(NegDif = sapply(outp, function(e) e$Negs)))
outp2 %>% group_by(NegDif) %>% summarise(m1=mean(DifComp.2))
