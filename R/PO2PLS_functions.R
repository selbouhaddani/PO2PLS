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
generate_params <- function(X, Y, r, rx, ry, type=c('o2m','random')){
  p = ncol(X)
  q = ncol(Y)
  type=match.arg(type)
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
        sig2E = 0.05,
        sig2F = 0.05
      )}))
  }
  if(type=="random"){
    list(
      W = OmicsPLS::orth(matrix(runif(p*r), p, r)),
      Wo = suppressWarnings(sign(rx)*OmicsPLS::orth(matrix(runif(p*max(1,rx)), p, max(1,rx)))),
      C = OmicsPLS::orth(matrix(runif(q*r), q, r)),
      Co = suppressWarnings(sign(ry)*OmicsPLS::orth(matrix(runif(q*max(1,rx)), q, max(1,ry)))),
      B = diag(sort(runif(r,1,4),decreasing = TRUE),r),
      SigT = diag(sort(runif(r,1,3),decreasing = TRUE),r),
      SigTo = sign(rx)*diag(sort(runif(max(1,rx),1,3),decreasing = TRUE),max(1,rx)),
      SigUo = sign(ry)*diag(sort(runif(max(1,ry),1,3),decreasing = TRUE),max(1,ry)),
      SigH = diag(0.1,r), #cov(H_UT)*diag(1,r),
      sig2E = 0.05,
      sig2F = 0.05
    )
  }
}

#' @export
generate_data <- function(N, params){
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

  MASS::mvrnorm(n = N,
                mu = rep(0,p+q),
                Sigma = Gamma %*% VarZ %*% t(Gamma) +
                  diag(rep(c(sig2E,sig2F),c(p,q))))
}

#' @export
Lemma <- function(X, SigmaZ, invS, Gamma){
  Gamma <- Gamma %*% SigmaZ
  VarZc <- SigmaZ - t(Gamma) %*% invS %*% Gamma
  EZc <- X %*% invS %*% Gamma
  return(list(EZc = EZc, VarZc = VarZc))
}

#' @export
E_step <- function(X, Y, params, use_lemma = FALSE){
  W = params$W
  C = params$C
  Wo = params$Wo
  Co = params$Co
  B = params$B
  SigT = params$SigT
  SigTo = (ssq(Wo)>0)*params$SigTo + 0.0001*(ssq(Wo)==0)
  SigUo = (ssq(Co)>0)*params$SigUo + 0.0001*(ssq(Co)==0)
  SigH = params$SigH
  sig2E = params$sig2E
  sig2F = params$sig2F
  SigU = SigT%*%B^2 + SigH

  N = nrow(X)
  p = nrow(W)
  q = nrow(C)
  r = ncol(W)
  rx = ncol(Wo)
  ry = ncol(Co)

  dataXY <- cbind(X,Y)

  Gamma = rbind(cbind(W, matrix(0,p,r), Wo, matrix(0,p,ry)),
                cbind(matrix(0,q,r), C, matrix(0,q,rx), Co))
  SigmaEF = diag(rep(c(sig2E,sig2F),c(p,q)))
  SigmaZ = blockm(
    blockm(
      blockm(SigT, SigT%*%B, SigU),
      matrix(0,2*r,rx), SigTo),
    matrix(0,2*r+rx,ry), SigUo)

  SigmaXY = Gamma %*% SigmaZ %*% t(Gamma) + SigmaEF
  invS <- solve(SigmaXY)
  if(use_lemma) tmp <- Lemma(cbind(X,Y), SigmaZ, invS, Gamma)
  invEF_Gamma <- rbind(cbind(W, matrix(0,p,r), Wo, matrix(0,p,ry))/sig2E,
                       cbind(matrix(0,q,r), C, matrix(0,q,rx), Co)/sig2F)
  inv2EF_Gamma <- rbind(cbind(W, matrix(0,p,r), Wo, matrix(0,p,ry))/(sig2E^2),
                        cbind(matrix(0,q,r), C, matrix(0,q,rx), Co)/(sig2F^2))
  invZtilde <- solve(solve(SigmaZ) +
                       t(Gamma) %*% rbind(cbind(W, matrix(0,p,r), Wo, matrix(0,p,ry))/sig2E,
                                          cbind(matrix(0,q,r), C, matrix(0,q,rx), Co)/sig2F))

  invS_Gamma <- invEF_Gamma - invEF_Gamma %*% invZtilde %*% crossprod(invEF_Gamma,Gamma)

  VarZc <- SigmaZ - t(Gamma %*% SigmaZ) %*% invS_Gamma %*% SigmaZ
  EZc <- dataXY %*% invS_Gamma %*% SigmaZ
  if(use_lemma) VarZc = tmp$VarZc
  if(use_lemma) EZc = tmp$EZc
  Szz = VarZc + crossprod(EZc)/N

  #invS_covEF <- diag(1,p+q) - invEF_Gamma %*% invZtilde %*% t(Gamma)
  #covEF = rbind(diag(sig2E,p), diag(0,q,p))
  mu_EF = dataXY - dataXY %*% invEF_Gamma %*% invZtilde %*% t(Gamma)
  #Ceeff = SigmaEF - t(SigmaEF) %*% invS_covEF + crossprod(mu_EF) / N

  #Cee <- sum(diag(Ceeff[1:p,1:p]))/p
  #Cff <- sum(diag(Ceeff[-(1:p),-(1:p)]))/q

  Cee <- sum(diag(
    crossprod(rbind(cbind(W, matrix(0,p,r), Wo, matrix(0,p,ry)),
                    cbind(matrix(0,q,r), 0*C, matrix(0,q,rx), 0*Co)))%*%invZtilde
  ))/p + OmicsPLS::ssq(mu_EF[,1:p])/N/p
  Cff <- sum(diag(
    crossprod(rbind(cbind(0*W, matrix(0,p,r), 0*Wo, matrix(0,p,ry)),
                    cbind(matrix(0,q,r), C, matrix(0,q,rx), Co)))%*%invZtilde
  ))/q + OmicsPLS::ssq(mu_EF[,-(1:p)])/N/q
  #covE = rbind(diag(sig2E,p), diag(0,q,p))
  #mu_E = cbind(X,Y) %*% invS %*% covE
  #Cee = sum(diag(diag(sig2E,p) - t(covE) %*% invS %*% covE + crossprod(mu_E) / N))/p

  #covF = rbind(diag(0,p,q), diag(sig2F,q))
  #mu_F = cbind(X,Y) %*% invS %*% covF
  #Cff = sum(diag(diag(sig2F,q) - t(covF) %*% invS %*% covF + crossprod(mu_F) / N))/q

  covH = rbind(0*W, C%*%SigH)
  invS_covH <- (covH/sig2F - invEF_Gamma %*% invZtilde %*% crossprod(invEF_Gamma,covH))
  mu_H = dataXY %*% invS_covH
  Chh = SigH - t(covH) %*% invS_covH + crossprod(mu_H) / N

  # covH = rbind(0*W, C%*%SigH)
  # mu_H = cbind(X,Y) %*% invS %*% covH
  # Chh = SigH - t(covH) %*% invS %*% covH + crossprod(mu_H) / N

  #solve(t(0))
  loglik = N*(p+q)*log(2*pi) +
   N * c(determinant(SigmaXY)$mod) +
   sum(diag((crossprod(cbind(X,Y)))%*%invS))
  loglik = - loglik/2

  comp_log <- -sum(diag(crossprod(cbind(X,Y)) - 2*crossprod(cbind(X,Y),EZc)%*%t(Gamma) + Gamma %*% Szz %*% t(Gamma)))
  list(
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
    comp_log = comp_log
  )
}

#' @export
M_step <- function(E_fit, params, X, Y){
  orth_x = OmicsPLS::ssq(params$Wo) > 0
  orth_y = OmicsPLS::ssq(params$Co) > 0
  #print(E_fit[-(1:4)])
  with(E_fit,{
    N = nrow(X)
    r = ncol(mu_T)
    rx = ncol(mu_To)
    ry = ncol(mu_Uo)
    params_old <- params
    #Q_old <- E_step(X, Y, params_old)$com
    params$W = OmicsPLS::orth(t(X) %*% mu_T - params$Wo%*%t(Stto),type = 'SVD')
    #L_W <- E_step(X,Y,params)$log; cat("Only W ", L_W, "  --  "); cat(E_step(X,Y,params)$comp-Q_old, "  -*-  ")
    params$Wo = suppressWarnings(orth_x*OmicsPLS::orth(t(X) %*% mu_To - params_old$W%*%Stto,type = 'SVD'))
    #L_Wo <- E_step(X,Y,params)$log; cat("Only Wo ", L_Wo, "  --  "); cat(E_step(X,Y,within(params,{W=params_old$W}))$comp-Q_old)
    #params$Wo = suppressWarnings(orth_x*OmicsPLS::orth(t(X) %*% mu_To - params$W%*%Stto,type = 'SVD'))
    #L_W_Wo <- E_step(X,Y,params)$log; cat("Wo after W ", L_W_Wo, "  --  "); cat(E_step(X,Y,params)$comp)
    #cat("\n")
    params$C = OmicsPLS::orth(t(Y) %*% mu_U - params$Co%*%t(Suuo),type = 'SVD')
    params$Co = suppressWarnings(orth_y*OmicsPLS::orth(t(Y) %*% mu_Uo - params$C%*%Suuo,type = 'SVD'))
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
    params
  })
}

#' @export
PO2PLS <- function(X, Y, r, rx, ry, steps = 1e2, tol = 1e-6, init_param='o2m', use_lemma = FALSE){
  params <- generate_params(X, Y, r, rx, ry, type = init_param)
  #params <- parms2
  params$Wo <- params$Wo
  params$Co <- params$Co
  err = logl = 0*0:steps
  for(i in 1:steps){
    E_next = E_step(X, Y, params, use_lemma = use_lemma)
    params_next = M_step(E_next, params, X, Y)
    #parms_next[-1] = params[-1]
    # if(i == 1) err[1] = mse(params_next[[1]],parms[[1]])
    # err[i+1] = mse(
    #   params_next[[1]],
    #   parms[[1]]%*%sign(abs(crossprod(params_next[[1]],parms[[1]]))>0.5))
    if(i == 1) logl[1] = E_next$logl
    logl[i+1] = E_next$logl# - err[i]
     #sum(mapply(function(e,f) OmicsPLS::mse(e, f), e=parms, f = parms_next))
    if(i > 1 && (logl[i+1]-logl[i]) < tol) break
    #if( (err[i]<-sum(mapply(mse, params, params_next))) < tol ) {params = params_next; break}
    params = params_next
  }
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
  list(params = params_next, err = err[1:i], logl = logl[0:i+1][-1])
  #list(params = params_next, err = err[1:i])
}

Lemma_old <- function(X, SigmaZ, SigmaE, Gamma){
  Gamma <- Gamma %*% SigmaZ
  VarZc <- solve(solve(SigmaZ) + t(Gamma) %*% solve(SigmaE) %*% Gamma)
  EZc <- X %*% solve(SigmaE) %*% Gamma %*% VarZc
  return(list(EZc = EZc, VarZc = VarZc))
}

E_stepold <- function(X, Y, params){
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

  N = nrow(X)
  p = nrow(W)
  q = nrow(C)
  r = ncol(W)
  rx = ncol(Wo)
  ry = ncol(Co)

  Gamma = rbind(cbind(W, matrix(0,p,r), Wo, matrix(0,p,ry)),
                cbind(matrix(0,q,r), C, matrix(0,q,rx), Co))
  Gamma_noise = rbind(1/sig2E*cbind(W, matrix(0,p,r), Wo, matrix(0,p,ry)),
                      1/sig2F*cbind(matrix(0,q,r), C, matrix(0,q,rx), Co))
  VarZ = blockm(
    blockm(
      blockm(SigT, SigT%*%B, SigU),
      matrix(0,2*r,rx), SigTo),
    matrix(0,2*r+rx,ry), SigUo)
  SigmaXY = Gamma %*% VarZ %*% t(Gamma) + diag(rep(c(sig2E,sig2F),c(p,q)))
  VarZc = solve(solve(VarZ) + t(Gamma) %*% Gamma_noise)
  #VarZc = VarZ - VarZ %*% t(Gamma) %*% solve(SigmaXY) %*% Gamma %*% VarZ
  EZc = cbind(1/sig2E*X,1/sig2F*Y) %*% Gamma %*% VarZc
  #EZc = cbind(X,Y) %*% solve(SigmaXY) %*% Gamma %*% VarZ
  Szz = VarZc + crossprod(EZc)/N

  loglik = N*(p+q)*log(2*pi) +
    c(determinant(SigmaXY)$mod) +
    sum(diag(crossprod(cbind(X,Y))%*%SigmaXY))
  loglik = - loglik/2
  list(
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
    loglik = loglik
  )
}
