#' PO2PLS: Probabilistic Two-way Orthogonal Partial Least Squares
#'
#' This package implements the O2PLS method in a probabilistic framework.
#' @author
#' Said el Bouhaddani (\email{s.el_bouhaddani@@lumc.nl}),
#' Jeanine Houwing-Duistermaat (\email{J.J.Houwing@@lumc.nl}),
#' Geurt Jongbloed (\email{G.Jongbloed@@tudelft.nl}),
#' Szymon Kielbasa (\email{S.M.Kielbasa@@lumc.nl}),
#' Hae-Won Uh (\email{H.Uh@@lumc.nl}).
#'
#' Maintainer: Said el Bouhaddani (\email{s.el_bouhaddani@@lumc.nl}).
#'
#' @docType package
#' @name PO2PLS-package
#' @keywords Probabilistic-O2PLS
#' @import OmicsPLS
NULL

#' @export
blockm<-function(A,B,C)
  #input: Matrices A,B,C
  #output: the block matrix
  # A    B
  #t(B)  C
{
  M = rbind(cbind(A,B),cbind(t(B),C))
  return(M)
}

#' @export
generate_params <- function(X, Y, r, rx, ry, alpha = 0.1, type=c('o2m','random')){
  type=match.arg(type)
  p = ifelse(is.matrix(X) | type != "random", ncol(X), X)
  q = ifelse(is.matrix(Y) | type != "random", ncol(Y), Y)
  if(type=="o2m"){
    return(with(OmicsPLS::o2m(X, Y, r, rx, ry, stripped=TRUE),{
      list(
        W = W.,
        Wo = suppressWarnings(OmicsPLS::orth(P_Yosc.)),
        C = C.,
        Co = suppressWarnings(OmicsPLS::orth(P_Xosc.)),
        B = abs(cov(Tt,U)%*%solve(cov(Tt)))*diag(1,r),
        SigT = cov(Tt)*diag(1,r),
        SigTo = sign(rx)*cov(T_Yosc.)*diag(1,max(1,rx)),
        SigUo = sign(ry)*cov(U_Xosc.)*diag(1,max(1,ry)),
        SigH = cov(H_UT)*diag(1,r),
        sig2E = (ssq(X)-ssq(Tt)-ssq(T_Yosc.))/prod(dim(X)) + 0.01,
        sig2F = (ssq(Y)-ssq(U)-ssq(U_Xosc.))/prod(dim(Y)) + 0.01
      )}))
  }
  if(type=="random"){
    outp <- list(
      W = OmicsPLS::orth(matrix(rnorm(p*r), p, r)+1),
      Wo = suppressWarnings(sign(rx)*OmicsPLS::orth(matrix(rnorm(p*max(1,rx)), p, max(1,rx))+seq(-p/2,p/2,length.out = p))),
      C = OmicsPLS::orth(matrix(rnorm(q*r), q, r)+1),
      Co = suppressWarnings(sign(ry)*OmicsPLS::orth(matrix(rnorm(q*max(1,rx)), q, max(1,ry))+seq(-q/2,q/2,length.out = q))),
      B = diag(sort(runif(r,1,4),decreasing = TRUE),r),
      SigT = diag(sort(runif(r,1,3),decreasing = TRUE),r),
      SigTo = sign(rx)*diag(sort(runif(max(1,rx),1,3),decreasing = TRUE),max(1,rx)),
      SigUo = sign(ry)*diag(sort(runif(max(1,ry),1,3),decreasing = TRUE),max(1,ry))
    )
    outp$SigH = diag(alpha/(1-alpha)*(mean(diag(outp$SigT%*%outp$B))),r) #cov(H_UT)*diag(1,r),
    with(outp, {
      c(outp,
        sig2E = alpha/(1-alpha)*(mean(diag(SigT)) + mean(diag(SigTo)))/p,
        sig2F = alpha/(1-alpha)*(mean(diag(SigT%*%B^2 + SigH)) + mean(diag(SigUo)))/q)
    })
  }
}

#' @export
generate_data <- function(N, params, distr = rnorm){
  W = params$W
  C = params$C
  Wo = params$Wo
  Co = params$Co
  B = params$B
  SigT = params$SigT
  SigTo = params$SigTo
  SigUo = params$SigUo
  SigH = params$SigH
  sig2E = params$sig2E
  sig2F = params$sig2F
  SigU = SigT%*%B^2 + SigH

  p = nrow(W)
  q = nrow(C)
  r = ncol(W)
  rx = ncol(Wo)
  ry = ncol(Co)

  Gamma = rbind(cbind(W, matrix(0,p,r), Wo, matrix(0,p,ry)),
                cbind(matrix(0,q,r), C, matrix(0,q,rx), Co))
  VarZ = blockm(
    blockm(
      blockm(SigT, SigT%*%B, SigU),
      matrix(0,2*r,rx), SigTo),
    matrix(0,2*r+rx,ry), SigUo)

  # MASS::mvrnorm(n = N,
  #               mu = rep(0,p+q),
  #               Sigma = Gamma %*% VarZ %*% t(Gamma) +
  #                 diag(rep(c(sig2E,sig2F),c(p,q))))

  Z <- scale(matrix(distr(N*(2*r+rx+ry)), N))
  Z <- Z %*% chol(VarZ)
  Z[,2*r+1:rx] <- sign(ssq(Wo))*Z[,2*r+1:rx]
  Z[,2*r+rx+1:ry] <- sign(ssq(Co))*Z[,2*r+rx+1:ry]

  EF <- cbind(scale(matrix(distr(N*p), N))*sqrt(sig2E), scale(matrix(distr(N*q), N))*sqrt(sig2F))

  Z %*% t(Gamma) + EF

}

#' @export
Lemma <- function(X, SigmaZ, invZtilde, Gamma, sig2E, sig2F, p, q, r, rx, ry){
  GammaEF <- Gamma
  GammaEF[1:p,c(1:r,2*r+1:rx)] <- 1/sig2E* GammaEF[1:p,c(1:r,2*r+1:rx)]
  GammaEF[-(1:p),c(r+1:r,2*r+rx+1:ry)] <- 1/sig2F* GammaEF[-(1:p),c(r+1:r,2*r+rx+1:ry)]

  #invSEF <- diag(1/diag(SigmaEF))
  #invS <- invSEF - invSEF %*% Gamma %*% solve(solve(SigmaZ) + t(Gamma)%*%invSEF%*%Gamma) %*% t(Gamma) %*% invSEF
  GGef <- t(Gamma) %*% GammaEF
  VarZc <- SigmaZ - (t(Gamma %*% SigmaZ) %*% GammaEF) %*% SigmaZ +
    (t(Gamma %*% SigmaZ) %*% GammaEF) %*% solve(solve(SigmaZ) + GGef) %*% GGef %*% SigmaZ

  EZc <- X %*% (GammaEF %*% SigmaZ)
  EZc <- EZc - X %*% ((GammaEF %*% invZtilde)  %*% (GGef %*% SigmaZ))
  # solve(t(0))
  return(list(EZc = EZc, VarZc = VarZc))
}



#' @export
E_step <- function(X, Y, params){
  ## retrieve parameters
  W = params$W
  C = params$C
  Wo = params$Wo
  Co = params$Co
  B = params$B
  SigT = params$SigT
  SigTo = (ssq(Wo)>0)*params$SigTo + 1e-10*(ssq(Wo)==0)
  SigUo = (ssq(Co)>0)*params$SigUo + 1e-10*(ssq(Co)==0)
  SigH = params$SigH
  sig2E = params$sig2E
  sig2F = params$sig2F
  SigU = SigT%*%B^2 + SigH

  ## define dimensions
  N = nrow(X)
  p = nrow(W)
  q = nrow(C)
  r = ncol(W)
  rx = ncol(Wo)
  ry = ncol(Co)

  ## concatenate data
  dataXY <- cbind(X,Y)

  ## Gamma is the generalized loading matrix, with PO2PLS structure
  Gamma = rbind(cbind(W, matrix(0,p,r), Wo, matrix(0,p,ry)),
                cbind(matrix(0,q,r), C, matrix(0,q,rx), Co))
  ## Gamma multiplied by inverse SigmaEF
  GammaEF <- Gamma
  GammaEF[1:p,c(1:r,2*r+1:rx)] <- 1/sig2E* GammaEF[1:p,c(1:r,2*r+1:rx)]
  GammaEF[-(1:p),c(r+1:r,2*r+rx+1:ry)] <- 1/sig2F* GammaEF[-(1:p),c(r+1:r,2*r+rx+1:ry)]
  GGef <- t(Gamma) %*% GammaEF

  ## diagonal cov matrix of (E,F), hopefully NOT NEEDED
  # SigmaEF = diag(rep(c(sig2E,sig2F),c(p,q)))
  ## ALMOST diagonal cov matrix of (T,U,To,Uo)
  SigmaZ = blockm(
    blockm(
      blockm(SigT, SigT%*%B, SigU),
      matrix(0,2*r,rx), SigTo),
    matrix(0,2*r+rx,ry), SigUo)

  ## inverse middle term lemma
  invZtilde <- solve(solve(SigmaZ) + GGef)

  ## Calculate conditional expectations with efficient lemma
  # print(all.equal(invS,invS_old))
  tmp <- Lemma(dataXY, SigmaZ, invZtilde, Gamma, sig2E, sig2F,p,q,r,rx,ry)
  # print(all.equal(Lemma(cbind(X,Y), SigmaZ, invS, Gamma, sig2E, sig2F,p,q,r,rx,ry), Lemma_old(cbind(X,Y), SigmaZ, invS, Gamma)))

  ## Define Szz as expected crossprod of Z
  VarZc = tmp$VarZc
  EZc = tmp$EZc
  Szz = VarZc + crossprod(EZc)/N

  ## For compatibility
  # invEF_Gamma <- rbind(cbind(W, matrix(0,p,r), Wo, matrix(0,p,ry))/sig2E,
  #                      cbind(matrix(0,q,r), C, matrix(0,q,rx), Co)/sig2F)
  # inv2EF_Gamma <- rbind(cbind(W, matrix(0,p,r), Wo, matrix(0,p,ry))/(sig2E^2),
  #                       cbind(matrix(0,q,r), C, matrix(0,q,rx), Co)/(sig2F^2))
  # invZtilde <- solve(solve(SigmaZ) +
  #                      t(Gamma) %*% rbind(cbind(W, matrix(0,p,r), Wo, matrix(0,p,ry))/sig2E,
  #                                         cbind(matrix(0,q,r), C, matrix(0,q,rx), Co)/sig2F))
  #
  # invS_Gamma <- invEF_Gamma - invEF_Gamma %*% invZtilde %*% crossprod(invEF_Gamma,Gamma)

  ## inverse in middle term in lemma
  # invZtilde <- solve(solve(SigmaZ) +
  #                      t(Gamma) %*% rbind(cbind(W, matrix(0,p,r), Wo, matrix(0,p,ry))/sig2E,
  #                                         cbind(matrix(0,q,r), C, matrix(0,q,rx), Co)/sig2F))

  ## Calculate cond mean of E,F
  # invS_covEF <- diag(1,p+q) - invEF_Gamma %*% invZtilde %*% t(Gamma)
  # covEF = rbind(diag(sig2E,p), diag(0,q,p))
  # mu_EF_old = dataXY - dataXY %*% invEF_Gamma %*% invZtilde %*% t(Gamma)
  mu_EF = dataXY
  mu_EF <- mu_EF - (dataXY %*% (GammaEF %*% invZtilde)) %*% t(Gamma)

  ## Calculate immediately expected crossprod of E,F
  # Ceeff_old = SigmaEF - t(SigmaEF) %*% invS_covEF + crossprod(mu_EF) / N
  # Ceeff = Gamma %*% solve(solve(SigmaZ) + t(Gamma)%*%invSEF%*%Gamma) %*% t(Gamma) +
  #   crossprod(mu_EF) / N
  # print(all.equal(mu_EF_old, mu_EF))
  # print(all.equal(Ceeff_old, Ceeff))

  ## Take trace of the matrix
  # Cee_old <- sum(diag(Ceeff_old[1:p,1:p]))/p
  Cee <- sum(diag(
    crossprod(rbind(cbind(W, matrix(0,p,r), Wo, matrix(0,p,ry)),
                    cbind(matrix(0,q,r), 0*C, matrix(0,q,rx), 0*Co)))%*%invZtilde
  ))/p + OmicsPLS::ssq(mu_EF[,1:p])/N/p
  # Cff_old <- sum(diag(Ceeff_old[-(1:p),-(1:p)]))/q
  Cff <- sum(diag(
    crossprod(rbind(cbind(0*W, matrix(0,p,r), 0*Wo, matrix(0,p,ry)),
                    cbind(matrix(0,q,r), C, matrix(0,q,rx), Co)))%*%invZtilde
  ))/q + OmicsPLS::ssq(mu_EF[,-(1:p)])/N/q
  # cat('Cee\n'); print(all.equal(Cee_old,Cee));
  # print(all.equal(Cff_old,Cff))

  # Cee <- sum(diag(
  #   crossprod(rbind(cbind(W, matrix(0,p,r), Wo, matrix(0,p,ry)),
  #                   cbind(matrix(0,q,r), 0*C, matrix(0,q,rx), 0*Co)))%*%invZtilde
  # ))/p + OmicsPLS::ssq(mu_EF[,1:p])/N/p
  # Cff <- sum(diag(
  #   crossprod(rbind(cbind(0*W, matrix(0,p,r), 0*Wo, matrix(0,p,ry)),
  #                   cbind(matrix(0,q,r), C, matrix(0,q,rx), Co)))%*%invZtilde
  # ))/q + OmicsPLS::ssq(mu_EF[,-(1:p)])/N/q


  #covE = rbind(diag(sig2E,p), diag(0,q,p))
  #mu_E = cbind(X,Y) %*% invS %*% covE
  #Cee = sum(diag(diag(sig2E,p) - t(covE) %*% invS %*% covE + crossprod(mu_E) / N))/p

  #covF = rbind(diag(0,p,q), diag(sig2F,q))
  #mu_F = cbind(X,Y) %*% invS %*% covF
  #Cff = sum(diag(diag(sig2F,q) - t(covF) %*% invS %*% covF + crossprod(mu_F) / N))/q

  covH = rbind(0*W, C%*%SigH)
  covHEF = rbind(0*W, C%*%SigH/sig2F)
  # invS_covH <- (covH/sig2F - invEF_Gamma %*% invZtilde %*% crossprod(invEF_Gamma,covH))
  # mu_H_old = dataXY %*% invS_covH
  # mu_H_old <- dataXY %*% invS %*% covH
  mu_H <- dataXY %*% covHEF
  mu_H <- mu_H - (dataXY %*% (GammaEF %*% invZtilde)) %*% (t(Gamma) %*% covHEF)
  # Chh_old = SigH - t(covH) %*% invS_covH + crossprod(mu_H) / N
  # Chh_old <- SigH - t(covH) %*% invS %*% covH + crossprod(mu_H) / N
  Chh <- SigH
  Chh <- Chh - t(covH) %*% covHEF
  Chh <- Chh + (t(covH) %*% GammaEF %*% invZtilde) %*% (t(Gamma) %*% covHEF)
  Chh <- Chh + crossprod(mu_H) / N
  # print(all.equal(mu_H_old, mu_H))
  # print(all.equal(Chh_old, Chh))

  ## diagonal cov matrix of (X,Y), hopefully NOT NEEDED
  # SigmaXY = Gamma %*% SigmaZ %*% t(Gamma) + SigmaEF
  ## INVERSE diagonal cov matrix of (E,F), hopefully NOT NEEDED
  # invSEF <- diag(1/diag(SigmaEF))
  ## INVERSE diagonal cov matrix of (X,Y), hopefully NOT NEEDED
  # invS <- invSEF - invSEF %*% Gamma %*% solve(solve(SigmaZ) + t(Gamma)%*%invSEF%*%Gamma) %*% t(Gamma) %*% invSEF
  # if(use_lemma == TRUE){solve(t(0))}
  ## log of det SigmaXY
  logdet <- log(det(diag(2*r+rx+ry) + GGef%*%SigmaZ))+p*log(sig2E)+q*log(sig2F)
  ## representation of SigmaXY %*% invS
  XYinvS <- ssq(cbind(X/sqrt(sig2E), Y/sqrt(sig2F)))
  XYinvS <- XYinvS - sum(diag(crossprod(dataXY %*% GammaEF) %*% invZtilde))
  ## Log likelihood
  loglik = N*(p+q)*log(2*pi) + N * logdet + XYinvS
  loglik = - loglik/2
  # print(all.equal(XYinvS, sum(diag(dataXY %*% invS %*% t(dataXY)))))
  # solve(t(0))

  # comp_log <- - N/2*(p+q)*log(2*pi)
  # comp_log <- comp_log - N/2*(p*log(sig2E)+q*log(sig2F))
  # comp_log <- comp_log - N/2*ssq(cbind(X/sqrt(sig2E), Y/sqrt(sig2F)))
  # comp_log <- comp_log + N*sum(diag(crossprod(EZc,dataXY)%*%GammaEF))
  # comp_log <- comp_log - N/2*sum(diag(GGef%*%Szz))
  list(
    EZc = EZc,
    Szz = Szz,
    mu_T = matrix(EZc[,1:r],N,r),
    mu_U = matrix(EZc[,r+1:r],N,r),
    mu_To = matrix(EZc[,2*r+1:rx],N,rx),
    mu_Uo = matrix(EZc[,2*r+rx+1:ry],N,ry),
    Stt = matrix(Szz[1:r, 1:r],r,r),
    Suu = matrix(Szz[r+1:r, r+1:r],r,r),
    Stoto = matrix(Szz[2*r+1:rx, 2*r+1:rx],rx,rx),
    Suouo = matrix(Szz[2*r+rx+1:ry, 2*r+rx+1:ry],ry,ry),
    Sut = matrix(Szz[r+1:r, 1:r],r,r),
    Stto = matrix(Szz[1:r, 2*r+1:rx],r,rx),
    Suuo = matrix(Szz[r+1:r, 2*r+rx+1:ry],r,ry),
    See = Cee,
    Sff = Cff,
    Shh = Chh,
    loglik = loglik,
    comp_log = 0#comp_log
  )
}

#' @export
M_step <- function(E_fit, params, X, Y, orth_type = c("SVD","QR"),
                   multiCM = FALSE, steps_multiCM = 10, verbose = FALSE,
                   debug = FALSE){
  orth_x = OmicsPLS::ssq(params$Wo) > 0
  orth_y = OmicsPLS::ssq(params$Co) > 0
  #print(E_fit[-(1:4)])
  orth_type = match.arg(orth_type)
  # params_old <- params
  with(E_fit,{
    tmp <- E_fit$loglik

    # N = nrow(X)
    r = ncol(mu_T)
    rx = ncol(mu_To)
    ry = ncol(mu_Uo)
    params_old <- params
    #Q_old <- E_step(X, Y, params_old)$com
    params$B = Sut %*% solve(Stt) * diag(1,r)
    params$SigT = Stt*diag(1,r)
    params$SigTo = Stoto*diag(1,rx)
    params$SigUo = Suouo*diag(1,ry)
    params$SigH = Shh*diag(1,r)#abs(Suu - 2*Sut%*%params_old$B + Stt%*%params_old$B^2)
    params$sig2E = See#abs(mean(diag(crossprod(X) - 2*crossprod(X,mu_T)%*%t(params_old$W) -
    #  2*crossprod(X,mu_To)%*%t(params_old$Wo) + 2*params_old$W%*%Stto%*%t(params_old$Wo) +
    #  params_old$W%*%Stt%*%t(params_old$W) + params_old$Wo%*%Stoto%*%t(params_old$Wo)))/N)
    params$sig2F = Sff#abs(mean(diag(crossprod(Y) - 2*crossprod(Y,mu_U)%*%t(params_old$C) -
    #  2*crossprod(Y,mu_Uo)%*%t(params_old$Co) + 2*params_old$C%*%Suuo%*%t(params_old$Co) +
    #  params_old$C%*%Suu%*%t(params_old$C) + params_old$Co%*%Suouo%*%t(params_old$Co)))/N)
    # #    solve(t(0))
    if(multiCM) tmp2 <- E_step(X, Y, params)$loglik

    params$W = OmicsPLS::orth(t(X) %*% mu_T - params$Wo%*%t(Stto),type = orth_type)
    params$C = OmicsPLS::orth(t(Y) %*% mu_U - params$Co%*%t(Suuo),type = orth_type)
    # params$W = OmicsPLS::orth(t(X - mu_To %*% t(params$Wo)) %*% mu_T, type = orth_type)
    # params$C = OmicsPLS::orth(t(Y - mu_Uo %*% t(params$Co)) %*% mu_U, type = orth_type)

    if(multiCM) tmp3 <- E_step(X, Y, params)$loglik

    params$Wo = suppressWarnings(orth_x*OmicsPLS::orth(t(X) %*% mu_To - params$W%*%Stto,type = orth_type))
    params$Co = suppressWarnings(orth_y*OmicsPLS::orth(t(Y) %*% mu_Uo - params$C%*%Suuo,type = orth_type))
    # params$Wo = suppressWarnings(orth_x*OmicsPLS::orth(t(X - mu_T %*% t(params$W)) %*% mu_To, type = orth_type))
    # params$Co = suppressWarnings(orth_y*OmicsPLS::orth(t(Y - mu_U %*% t(params$C)) %*% mu_Uo, type = orth_type))

    # objfun2 <- function(Warg){list(value = E_step(X,Y,within(params,{W=as.matrix(Warg$Qt[1:p,1:r])}))$com)}
    # params$W <- orth(GrassmannOptim(objfun2, W = list(Qt = orth(cbind(params$W,matrix(rnorm(p*(p-r)),p))), dim=c(r,p)), max_iter = 1e2)$Qt[,1:r])
    # objfun2 <- function(Warg){list(value = E_step(X,Y,within(params,{C=as.matrix(Warg$Qt[1:q,1:r])}))$com)}
    # params$C <- orth(GrassmannOptim(objfun2, W = list(Qt = orth(cbind(params$C,matrix(rnorm(q*(q-r)),q))), dim=c(r,q)), max_iter = 1e2)$Qt[,1:r])
    # objfun2 <- function(Warg){list(value = E_step(X,Y,within(params,{Wo=as.matrix(Warg$Qt[1:p,1:rx])}))$com)}
    # params$Wo <- orth(GrassmannOptim(objfun2, W = list(Qt = orth(cbind(params$Wo,matrix(rnorm(p*(p-rx)),p))), dim=c(rx,p)), max_iter = 1e2)$Qt[,1:rx])
    # objfun2 <- function(Warg){list(value = E_step(X,Y,within(params,{Co=as.matrix(Warg$Qt[1:q,1:ry])}))$com)}
    # params$Co <- orth(GrassmannOptim(objfun2, W = list(Qt = orth(cbind(params$Co,matrix(rnorm(q*(q-ry)),q))), dim=c(ry,q)), max_iter = 1e2)$Qt[,1:ry])

    # Gamma <- (crossprod(cbind(X,Y), EZc) %*% solve(Szz))
    # params$W <- orth(Gamma[1:p,1:r])
    # params$Wo <- orth(Gamma[1:p,2*r+1:rx])
    # params$C <- orth(Gamma[-(1:p),r+1:r])
    # params$Co <- orth(Gamma[-(1:p),2*r+rx+1:ry])

    if(multiCM) tmp4 <- E_step(X, Y, params)$loglik

    if(debug) solve(t(0))

    if(multiCM && any(c(tmp4 < tmp, tmp4 < tmp))){
      if(verbose){
        cat("\n Old likelihood after var's: \n")
        cat(tmp2)

        cat(" \n Intermediate likelihood: \n")
        cat(tmp3)

        cat("\n New likelihood: \n")
        cat(tmp4)

        cat("\n \n Conclusion: \n Failed to converge in step ")
        cat(which(c(tmp3 < tmp2, tmp4 < tmp3)))


        #if(tmp2 <= tmp) params[c("W","C")] <- params_old[c("W","C")]
        #if(tmp3 <= tmp2) params[c("Wo","Co")] <- params_old[c("Wo","Co")]

        cat("\n --------------------- \n ")
      }
      if(steps_multiCM > 0){
        for(i in 1:steps_multiCM){
          params$W = OmicsPLS::orth(t(X) %*% mu_T - params$Wo%*%t(Stto),type = orth_type)
          params$C = OmicsPLS::orth(t(Y) %*% mu_U - params$Co%*%t(Suuo),type = orth_type)
          # params$W = OmicsPLS::orth(t(X - mu_To %*% t(params$Wo)) %*% mu_T, type = orth_type)
          # params$C = OmicsPLS::orth(t(Y - mu_Uo %*% t(params$Co)) %*% mu_U, type = orth_type)
          params$Wo = suppressWarnings(orth_x*OmicsPLS::orth(t(X) %*% mu_To - params$W%*%Stto,type = orth_type))
          params$Co = suppressWarnings(orth_y*OmicsPLS::orth(t(Y) %*% mu_Uo - params$C%*%Suuo,type = orth_type))
          # params$Wo = suppressWarnings(orth_x*OmicsPLS::orth(t(X - mu_T %*% t(params$W)) %*% mu_To, type = orth_type))
          # params$Co = suppressWarnings(orth_y*OmicsPLS::orth(t(Y - mu_U %*% t(params$C)) %*% mu_Uo, type = orth_type))
        }
      }
    }

    params
  })
}

#' @export
jitter_params <- function(params, amount = NULL){
  suppressWarnings(params[1:4] <- lapply(params[1:4], function(e) sign(ssq(e))*orth(jitter(e,amount = 1))))
  params
}

#' @export
diagnostics.PO2PLS <- function(th, th0){
  c(
    W = max(abs(crossprod(th$W,th0$W))),
    C = max(abs(crossprod(th$C,th0$C))),
    Wo = max(abs(crossprod(th$Wo,th0$Wo))),
    Co = max(abs(crossprod(th$Co,th0$Co))),
    varTo_T = sum(diag(th$SigTo))/sum(diag(th$SigT))/ncol(th$Wo)*ncol(th$W),
    varUo_U = sum(diag(th$SigUo))/sum(diag(th$SigT%*%th$B+th$SigH))/ncol(th$Co)*ncol(th$C),
    varU_T = sum(diag(th$SigT%*%th$B+th$SigH))/sum(diag(th$SigT))
  )
}

#' @export
PO2PLS <- function(X, Y, r, rx, ry, steps = 1e2, tol = 1e-6, init_param='o2m',
                   orth_type = "SVD", multiCM = FALSE, steps_multiCM = 10,
                   random_restart = FALSE, verbose = FALSE){
  if(all(c("W","Wo","C","Co","B","SigT","SigTo","SigUo","SigH","sig2E","sig2F") %in% names(init_param))) {cat('using old fit \n'); params <- init_param}
  else {params <- generate_params(X, Y, r, rx, ry, type = init_param)}
  #params <- parms2
  #params$Wo <- params$Wo
  #params$Co <- params$Co
  #diags <- matrix(NA,0,7)
  logl = 0*0:steps
  tic <- proc.time()
  print(paste('started',date()))

  i_rr <- 0
  random_restart_original <- random_restart
  random_restart <- TRUE
  while(random_restart){

    if(i_rr > 0) {
      message("Log-likelihood: ", logl[i+1])
      message(paste("random restart no",i_rr))
    }
    params_max <- params
    for(i in 1:steps){
      E_next = E_step(X, Y, params)
      params_next = M_step(E_next, params, X, Y, orth_type = orth_type,
                           multiCM = multiCM, steps_multiCM = steps_multiCM,
                           verbose = verbose)
      # if( i == 300) M_step(E_next, params, X, Y, orth_type = orth_type,
      #                      multiCM = multiCM, steps_multiCM = steps_multiCM,
      #                      verbose = verbose, debug = TRUE)

      params_next$B <- abs(params_next$B)
      #params_next$C <- params_next$C %*% diag(signB,r)
      #ordSB <- order(diag(params_next$SigT %*% params_next$B), decreasing = TRUE)
      #params_next$W <- params_next$W[,ordSB]
      #params_next$C <- params_next$C[,ordSB]

      # diags <- rbind(diags,diagnostics.PO2PLS(params_next, parms))

      if(i == 1) logl[1] = E_next$logl
      logl[i+1] = E_next$logl# - err[i]
       #sum(mapply(function(e,f) OmicsPLS::mse(e, f), e=parms, f = parms_next))
      if(i > 1 && abs(logl[i+1]-logl[i]) < tol) break
      if(i %in% c(1e2, 1e3, 5e3, 1e4, 4e4)) {
        print(data.frame(row.names = 1, steps = i, time = unname(proc.time()-tic)[3], diff = logl[i+1]-logl[i], logl = logl[i+1]))
      }
      if(logl[i+1] > max(logl[1:i])) params_max <- params_next
      params = params_next
    }
    if(!any(diff(logl[-1]) < -1e-10) | !random_restart_original) {
      random_restart = FALSE
      break
    }
    i_rr <- i_rr + 1
    params <- jitter_params(params)
    params[-(1:4)] <- generate_params(X, Y, r, rx, ry, type = 'r')[-(1:4)]
  }
  # params <- params_max
  signB <- sign(diag(params$B))
  params$B <- params$B %*% diag(signB,r)
  params$C <- params$C %*% diag(signB,r)
  ordSB <- order(diag(params$SigT %*% params$B), decreasing = TRUE)
  params$W <- params$W[,ordSB]
  params$C <- params$C[,ordSB]
  #message("Nr steps was ", i, "; error was ", signif(err[i+1],4))
  message("Nr steps was ", i)
  message("Negative increments: ", any(diff(logl[-1]) < -1e-10),
          "; Last increment: ", signif(logl[i+1]-logl[i],4))
  message("Log-likelihood: ", logl[i+1])
  outputt <- list(params = params_next, logl = logl[0:i+1][-1])
  class(outputt) <- "PO2PLS"
  return(outputt)
  #list(params = params_next, err = err[1:i])
}

plot_accur.PO2PLS <- function(fit){
  library(ggplot2)
  library(gridExtra)
  fit_o2m <- o2m(X,Y,ncol(parms$W),ncol(parms$Wo),ncol(parms$Co))
  g1 <- ggplot(reshape2::melt(fit$diags[,1:4]), aes(x=Var1,y=value)) + geom_line(aes(col=Var2,linetype=grepl('o',Var2)))
  g2 <- ggplot(reshape2::melt(fit$diags[,5:6]), aes(x=Var1,y=value)) + geom_line(aes(col=Var2))
  g3 <- ggplot(reshape2::melt(fit$diags[,7]), aes(x=1:nrow(fit$diags),y=value)) + geom_line()
  g4 <- qplot(x=1:length(fit$logl), y=fit$logl, geom='line')
  grid.arrange(g1,g2,g3,g4)
  print("### MAX ABS CROSSPROD WITH TRUE LOADINGS")
  print(apply(fit$diags,2,function(e) c(min=which.min(e), max=which.max(e))))
  print("### MAX VALUES FOR O2PLS")
  print(c(
    W = max(abs(crossprod(fit_o2m$W.,parms$W))),
    C = max(abs(crossprod(fit_o2m$C.,parms$C))),
    Wo = max(abs(crossprod(orth(fit_o2m$P_Y),parms$Wo))),
    Co = max(abs(crossprod(orth(fit_o2m$P_X),parms$Co)))
  ))
  print("### MAX CROSSPROD JOINT AND ORTHOGONAL SPACE")
  print(c(W=max(abs(crossprod(parms$W,parms$Wo))), C=max(abs(crossprod(parms$C,parms$Co)))))
}

cov.PO2PLS <- function(fit){
  with(fit$par,
       blockm(W%*%SigT%*%t(W)+Wo%*%SigTo%*%t(Wo) ,
              W%*%SigT%*%B%*%t(C) ,
              C%*%(SigT%*%B^2+SigH)%*%t(C)+Co%*%SigUo%*%t(Co)))
}
