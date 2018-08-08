coverW <- function(W0, W){
  W0 <- W0[,apply(crossprod(W0,W)^2,2,which.max)]
  W0 <- W0 %*% sign(crossprod(W0,W)*diag(1,ncol(W)))
}

seBoot <- function(K, X, Y, fit){
  require(parallel)
  require(OmicsPLS)
  require(magrittr)
  require(fBasics)
  tic <- proc.time()
  r <- ncol(fit$par$W)
  rx <- ncol(fit$par$Wo)*sign(ssq(fit$par$Wo))
  ry <- ncol(fit$par$Co)*sign(ssq(fit$par$Co))
  bootsam <- matrix(sample(1:nrow(X), size=K*nrow(X),replace=TRUE),nrow=K)
  reps <- parallelsugar::mclapply(mc.cores=detectCores(), 1:nrow(bootsam), function(i) PO2PLS(X[bootsam[i,],], Y[bootsam[i,],], r, rx, ry, 1e3, 1e-4))
  repsW <- lapply(reps, function(e) coverW(e$par$W,fit$par$W))
  repsC <- lapply(reps, function(e) coverW(e$par$C,fit$par$C))
  repsWo <- lapply(reps, function(e) coverW(e$par$C,fit$par$Wo))
  repsCo <- lapply(reps, function(e) coverW(e$par$C,fit$par$Co))
  seW <- sapply(1:r, function(i.comp) matrix(lapply(1:ncol(X), function(iii) lapply(1:K, function(ii) repsW[[ii]][iii,i.comp])) %>% unlist,K) %>% colSds)
  seC <- sapply(1:r, function(i.comp) matrix(lapply(1:ncol(Y), function(iii) lapply(1:K, function(ii) repsC[[ii]][iii,i.comp])) %>% unlist,K) %>% colSds)
  seWo <- sapply(1:rx, function(i.comp) matrix(lapply(1:ncol(X), function(iii) lapply(1:K, function(ii) repsWo[[ii]][iii,i.comp])) %>% unlist,K) %>% colSds)
  seCo <- sapply(1:ry, function(i.comp) matrix(lapply(1:ncol(Y), function(iii) lapply(1:K, function(ii) repsCo[[ii]][iii,i.comp])) %>% unlist,K) %>% colSds)
  list(repsW=repsW, repsC=repsC, repsWo=repsWo, repsCo=repsCo, seW = seW, seC = seC,seWo = seWo, seCo = seCo, time = (proc.time() - tic)[3])
}
seEmpir <- function(K, X, Y, fit){
  require(parallel)
  require(OmicsPLS)
  require(magrittr)
  require(fBasics)
  tic <- proc.time()
  N <- nrow(X)
  r <- ncol(fit$par$W)
  rx <- ncol(fit$par$Wo)
  ry <- ncol(fit$par$Co)
  reps <- parallelsugar::mclapply(mc.cores=detectCores(), 1:K, function(i) {
    Dat = generate_data(N, parms, distr)
    X = scale(Dat[,1:p], scale = F)
    Y = scale(Dat[,-(1:p)], scale = F)
    PO2PLS(X, Y, r, rx, ry, 1e3, 1e-4)
  })
  repsW <- lapply(reps, function(e) coverW(e$par$W,fit$par$W))
  repsC <- lapply(reps, function(e) coverW(e$par$C,fit$par$C))
  repsWo <- lapply(reps, function(e) coverW(e$par$C,fit$par$Wo))
  repsCo <- lapply(reps, function(e) coverW(e$par$C,fit$par$Co))
  seW <- sapply(1:r, function(i.comp) matrix(lapply(1:ncol(X), function(iii) lapply(1:K, function(ii) repsW[[ii]][iii,i.comp])) %>% unlist,K) %>% colSds)
  seC <- sapply(1:r, function(i.comp) matrix(lapply(1:ncol(Y), function(iii) lapply(1:K, function(ii) repsC[[ii]][iii,i.comp])) %>% unlist,K) %>% colSds)
  seWo <- sapply(1:rx, function(i.comp) matrix(lapply(1:ncol(X), function(iii) lapply(1:K, function(ii) repsWo[[ii]][iii,i.comp])) %>% unlist,K) %>% colSds)
  seCo <- sapply(1:ry, function(i.comp) matrix(lapply(1:ncol(Y), function(iii) lapply(1:K, function(ii) repsCo[[ii]][iii,i.comp])) %>% unlist,K) %>% colSds)
  list(repsW=repsW, repsC=repsC, repsWo=repsWo, repsCo=repsCo, seW = seW, seC = seC,seWo = seWo, seCo = seCo, time = (proc.time() - tic)[3])
}

f <- function(niets){
#set.seed(13254)
library(OmicsPLS)
distr_list <- list(norm=rnorm,
                   t=function(n) rt(n, df=2),
                   pois=function(n) rpois(n, lambda=1),
                   binom=function(n) rbinom(n, size = 2, prob = .25))
N = 5000
p = 25
q = 15
r = 5
rx = 5
ry = 0
noise_alpha = 0.5
distr_name <- names(distr_list)[1]
distr = distr_list[[distr_name]]
#load(paste0('parms_N500_p',p,'_a',noise_alpha*100,'_s',ifelse(noise_alpha == 0.5, 87548, 29867),'L.RData'))
parms <- generate_params(p,q,r,rx,ry, type = 'r',alpha=noise_alpha)
#parms$SigH = diag(rep(mean(diag(parms$SigH)),r))
#parms$B <- parms$B + diag(10,r)
#tmp_factor = parms$B^2 + parms$SigH
#parms$B = parms$B %*% solve(sqrt(tmp_factor))
#parms$SigH = parms$SigH %*% solve(tmp_factor)
Dat = generate_data(N, parms, distr)
X = scale(Dat[,1:p], scale = F)
Y = scale(Dat[,-(1:p)], scale = F)

time_o2m <- system.time(print(pryr::mem_change(fit_o2m <- o2m(X,Y,r,rx,ry,stripped=T))))
time_po2m <- system.time(print(pryr::mem_change(fit <- PO2PLS(X = X, Y = Y, r = r, rx = rx, ry = ry,
#                                                              steps = 1e3, tol = 1e-4,
                                                              init_param = 'o2m'))))

cmax <- function(A,B){
  apply(Reduce(crossprod, list(A,B)), 2, function(e) max(abs(e)))
}

(Wmax <- cmax(fit$par$W, parms$W) - cmax(fit_o2m$W., parms$W))
(Cmax <- cmax(fit$par$C, parms$C) - cmax(fit_o2m$C., parms$C))
(Womax <- cmax(fit$par$Wo, parms$Wo) - cmax(orth(fit_o2m$P_Y), parms$Wo))
(Comax <- cmax(fit$par$Co, parms$Co) - cmax(orth(fit_o2m$P_X), parms$Co))

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
library(tidyverse)
N
K = 20
SEs10 = unlist(seBoot(K,X[1:10,],Y[1:10,],fit)[5:8])
SEs10 <- SEs10[SEs10>0]
aSEs10 <- sqrt(diag(as.matrix(Matrix::nearPD(-solve(variances.PO2PLS(fit,Dat[1:10,])$Iobs))$mat)))
SEs50 = unlist(seBoot(K,X[1:50,],Y[1:50,],fit)[5:8])
SEs50 <- SEs50[SEs50>0]
aSEs50 <- sqrt(diag(as.matrix(Matrix::nearPD(-solve(variances.PO2PLS(fit,Dat[1:50,])$Iobs))$mat)))
SEs500 = unlist(seBoot(K,X[1:500,],Y[1:500,],fit)[5:8])
SEs500 <- SEs500[SEs500>0]
aSEs500 <- sqrt(diag(as.matrix(Matrix::nearPD(-solve(variances.PO2PLS(fit,Dat[1:500,])$Iobs))$mat)))
SEs5000 = unlist(seBoot(K,X,Y,fit)[5:8])
SEs5000 <- SEs5000[SEs5000>0]
aSEs5000 <- sqrt(diag(as.matrix(Matrix::nearPD(-solve(variances.PO2PLS(fit,Dat)$Iobs))$mat)))

seDat <- as.matrix(rbind(t(c(SEs10,xxType="Bootstrap",xxScenario="N=10")),t(c(aSEs10,"Asymptotic","N=10"))))
seDat %<>% rbind(rbind(t(c(SEs50,xxType="Bootstrap",xxScenario="N=50")),t(c(aSEs50,"Asymptotic","N=50"))))
seDat %<>% rbind(rbind(t(c(SEs500,xxType="Bootstrap",xxScenario="N=500")),t(c(aSEs500,"Asymptotic","N=500"))))
seDat %<>% rbind(rbind(t(c(SEs5000,xxType="Bootstrap",xxScenario="N=5000")),t(c(aSEs5000,"Asymptotic","N=5000"))))
colnames(seDat) %<>% str_sub(start = 3)
seDat %<>% as_data_frame
seDat %<>% mutate_at(vars(-Type,-Scenario), as.numeric)
seDat %<>% mutate_at(vars(Type,Scenario), as.factor)
seDat <- reshape2::melt(seDat, id.var=c("Type","Scenario"))
ggplot(data = seDat, aes(x=variable, y=value)) +
  geom_point(aes(col=Type, shape=Type), size=2) +
  facet_grid(Scenario ~ ., scales = "free_y") +
  theme_bw()



library(magrittr)
outp <- replicate(50,{
Dat = generate_data(N, parms, distr)
X = scale(Dat[,1:p], scale = F)
Y = scale(Dat[,-(1:p)], scale = F)
fit_NULL <- PO2PLS(X,Y,3,0,0,tol = 1e-2)
fit_ALT <- PO2PLS(X,Y,4,0,0,tol = 1e-2)
return(2*(tail(fit_ALT$logl,1) - tail(fit_NULL$logl,1)))
}) %>% invisible
hist(outp, freq = F)
