f <- function(none){
#set.seed(13254)
library(OmicsPLS)

p = 30
q = 20
r = 5
rx = 5
ry = 0
parms = generate_params(matrix(0,1,p),matrix(0,1,q),r,rx,ry, type = 'r')
Dat = generate_data(1e3, parms)
X = scale(Dat[,1:p], scale = F)
Y = scale(Dat[,-(1:p)], scale = F)

time_o2m <- system.time(print(pryr::mem_change(fit_o2m <- o2m(X,Y,r,rx,ry,stripped=T))))
time_po2m <- system.time(print(pryr::mem_change(fit <- PO2PLS(X, Y, r, rx, ry, 5e3, tol = 1e-6, 'o2m'))))

cmax <- function(A,B){
  apply(Reduce(crossprod, list(A,B)), 2, function(e) max(abs(e)))
}

Wmax <- cmax(fit$par$W, parms$W) - cmax(fit_o2m$W., parms$W)
Cmax <- cmax(fit$par$C, parms$C) - cmax(fit_o2m$C., parms$C)
Womax <- cmax(fit$par$Wo, parms$Wo) - cmax(orth(fit_o2m$P_Y), parms$Wo)
Comax <- cmax(fit$par$Co, parms$Co) - cmax(orth(fit_o2m$P_X), parms$Co)

list(W = Wmax, C = Cmax, Wo=Womax, Co=Comax, Negs = any(diff(fit$logl)<0), time = c(o2m = time_o2m, po2m <- time_po2m))

}

library(parallel)
library(tidyverse)
library(OmicsPLS)
source('R/PO2PLS_functions.R')
system.time(outp <- parallelsugar::mclapply(mc.cores=4, 1:200, f))
outp2 <- cbind(data.frame(DifComp=t(sapply(outp, function(e) e$W))), data.frame(NegDif = sapply(outp, function(e) e$Negs)))
outp2 %>% group_by(NegDif) %>% summarise(m1=mean(DifComp.2))
