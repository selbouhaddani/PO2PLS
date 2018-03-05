f <- function(niets){
#set.seed(13254)
library(OmicsPLS)
distr_list <- list(norm=rnorm,
                   t=function(n) rt(n, df=2),
                   pois=function(n) rpois(n, lambda=1),
                   binom=function(n) rbinom(n, size = 2, prob = .25))
N = 10
p = 1e4
q = 5
r = 1
rx = 1
ry = 1
noise_alpha = 0.1
distr_name <- names(distr_list)[1]
distr = distr_list[[distr_name]]
#load(paste0('parms_N500_p',p,'_a',noise_alpha*100,'_s',ifelse(noise_alpha == 0.5, 87548, 29867),'L.RData'))
parms <- generate_params(p,q,r,rx,ry, type = 'r',alpha=noise_alpha)
Dat = generate_data(N, parms, distr)
X = scale(Dat[,1:p], scale = F)
Y = scale(Dat[,-(1:p)], scale = F)

time_o2m <- system.time(print(pryr::mem_change(fit_o2m <- o2m(X,Y,r,rx,ry,stripped=T))))
time_po2m <- system.time(print(pryr::mem_change(fit <- PO2PLS(X = X, Y = Y, r = r, rx = rx, ry = ry,
                                                              steps = 5e3, tol = 1e-6,
                                                              init_param = 'o2m'))))

cmax <- function(A,B){
  apply(Reduce(crossprod, list(A,B)), 2, function(e) max(abs(e)))
}

Wmax <- cmax(fit$par$W, parms$W) - cmax(fit_o2m$W., parms$W)
Cmax <- cmax(fit$par$C, parms$C) - cmax(fit_o2m$C., parms$C)
Womax <- cmax(fit$par$Wo, parms$Wo) - cmax(orth(fit_o2m$P_Y), parms$Wo)
Comax <- cmax(fit$par$Co, parms$Co) - cmax(orth(fit_o2m$P_X), parms$Co)

list(fits = list(PO2PLS = fit, O2PLS = fit_o2m),
     W = Wmax, C = Cmax, Wo=Womax, Co=Comax, Negs = any(diff(fit$logl)<0),
     time = c(o2m = time_o2m[3], po2m = time_po2m[3]))

}

# library(parallel)
# library(tidyverse)
# library(OmicsPLS)
# source('R/PO2PLS_functions.R')
# system.time(outp <- parallelsugar::mclapply(mc.cores=4, 1:10, f))
# outp2 <- cbind(data.frame(DifComp=t(sapply(outp, function(e) e$W))), data.frame(NegDif = sapply(outp, function(e) e$Negs)))
# outp2 %>% group_by(NegDif) %>% summarise(m1=median(DifComp.1), s1 = mad(DifComp.1))
