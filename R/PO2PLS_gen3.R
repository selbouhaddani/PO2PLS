library(OmicsPLS)
library(PO2PLS)

# p <- 10
# q <- 11
# p3 <- 12
# r <- 1
# rx <- 1
# ry <- 1
# prm <- generate_params(p, q, r, rx, ry)
# prm3 <- generate_params(p, p3, r, rx, ry)
# 
# prm
# prm3


gen_par3 <- function(p1, p2, p3, r, rx1, rx2, rx3, alpha = 0.1){
  prm <- generate_params(p1, p2, r, rx1, rx2, alpha)
  prm3 <- generate_params(p3, p2, r, rx3, rx2, alpha)
  list(W1 = prm$W, 
       W2 = prm$C,
       W3 = prm3$W,
       Wo1 = prm$Wo,
       Wo2 = prm$Co,
       Wo3 = prm3$Wo,
       SigT = prm$SigT,
       SigTo1 = prm$SigTo,
       SigTo2 = prm$SigUo,
       SigTo3 = prm3$SigTo,
       sig2E1 = prm$sig2E,
       sig2E2 = prm$sig2F,
       sig2E3 = prm3$sig2E,
       beta_T = sort(runif(r, -1, 1))
       )
}

gen_dat3 <- function(N, params, alpha_out = 0.1, sparse = FALSE){
  W1 <- params$W1
  W2 <- params$W2
  W3 <- params$W3
  Wo1 <- params$Wo1
  Wo2 <- params$Wo2
  Wo3 <- params$Wo3
  
  r <- ncol(W1)
  rx1 <- ncol(Wo1)
  rx2 <- ncol(Wo2)
  rx3 <- ncol(Wo3)
  p1 <- nrow(W1)
  p2 <- nrow(W2)
  p3 <- nrow(W3)
  
  if(sparse){
    W1[-(1:(p1/4)),] <- 0
    W2[-(1:(p2/4)),] <- 0
    W3[-(1:(p3/4)),] <- 0
  }
  # print(W1)
  
  SigT = params$SigT
  SigTo1 = params$SigTo1 + 1e-06 * SigT[1] * (params$SigTo1[1] == 0)
  SigTo2 = params$SigTo2 + 1e-06 * SigT[1] * (params$SigTo2[1] == 0)
  SigTo3 = params$SigTo3 + 1e-06 * SigT[1] * (params$SigTo3[1] == 0)
  
  N <- 3/2 * N
  
  Tt <- matrix(rnorm(N * r), N, r) %*% chol(SigT)
  To1 <- matrix(rnorm(N * rx1), N, rx1) %*% chol(SigTo1)
  To2 <- matrix(rnorm(N * rx2), N, rx2) %*% chol(SigTo2)
  To3 <- matrix(rnorm(N * rx3), N, rx3) %*% chol(SigTo3)
  
  E1 <- matrix(rnorm(N * p1), N, p1) * sqrt(params$sig2E1)
  E2 <- matrix(rnorm(N * p2), N, p2) * sqrt(params$sig2E2)
  E3 <- matrix(rnorm(N * p3), N, p3) * sqrt(params$sig2E3)
  
  X1 <- Tt %*% t(W1) + To1 %*% t(Wo1) + E1
  X2 <- Tt %*% t(W2) + To2 %*% t(Wo2) + E2
  X3 <- Tt %*% t(W3) + To3 %*% t(Wo3) + E3
  
  outc <- Tt %*% params$beta_T
  outc <- outc + rnorm(N, 0, sd = sqrt( alpha_out/(1-alpha_out)*var(outc) ))
  
  indx1 <- order(outc)[1:(N/3)]
  indx2 <- order(-outc)[1:(N/3)]
  indxN <- c(indx1, indx2)
  
  outc_bin <- 1*(outc > 0)
  outc_bin <- outc_bin[indxN]
  
  #print((indxN))
  #print(str(X1))
  
  X1 <- X1[indxN,]
  X2 <- X2[indxN,]
  X3 <- X3[indxN,]
  
  return(list(X1 = X1, X2 = X2, X3 = X3, outc = outc_bin))
  
}


# prm <- gen_par3(10, 11, 12, 2, 1, 1, 1, alpha = 0.1)
# dat <- gen_dat3(100, prm, sparse = TRUE)



