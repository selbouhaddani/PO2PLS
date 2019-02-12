#' @keyword internal
#' @export
tr <- function(X){
  sum(diag(X))
}

# output:
# LVs
# parameters
# meta_data: time, call, type, total var
#' @export
PO2PLS_to_po2m <- function(fit, X, Y){
  p <- ncol(X)
  q <- ncol(Y)
  SigU <- with(fit$par, SigT %*% B^2 + SigH)
  LVs <- E_step(X, Y, fit$par)
  R2s <- c(
    R2X = with(fit$par, tr(SigT + SigTo)/tr(SigT+SigTo+p*sig2E)),
    R2Y = with(fit$par, tr(SigU + SigUo)/tr(SigU+SigUo+q*sig2F)),
    R2Xj = with(fit$par, tr(SigT)/tr(SigT+SigTo+p*sig2E)),
    R2Yj = with(fit$par, tr(SigU)/tr(SigU+SigUo+q*sig2F)),
    R2Yhat = with(fit$par, tr(SigT %*% B)/tr(SigU))
  )
  outp <- list(parameters = fit$par,
    latent_vars = LVs[3:15],
    meta_data = list(
      loglikelihood = fit$log,
      explained_vars = R2s,
      comps = with(fit$par, c(r=ncol(W), rx=ncol(Wo),ry=ncol(Co))),
      time = fit$flags$time,
      call = fit$flags$call,
      convergence = fit$flags$converg
    ))
  class(outp) <- "po2m"
  return(outp)
}

# print.po2m
# type fit
# number of components
# time
# loglikelihood sequence + last nr + last increment
print.po2m <- function(x, ...){
  cat('\n')
  cat(paste0("Call: ", x$meta_data$call))
  cat('\n')

  cat('\n')
  cat(paste0("This is a Probabilistic ", ifelse(sum(x$meta_data$comps[2:3])==0, "", "O2"), "PLS fit"))
  cat('\n')

#  cat('\n')
  cat(paste0("With r=",x$meta$comps[1],", rx=",x$meta$comps[2]," and ry=", x$meta$comps[3]," components"))
  cat('\n')

#  cat('\n')
  cat(paste0("Elapsed time: ", signif(x$meta$time,3), " sec"))
  cat('\n')

#  cat('\n')
  cat(paste0("The log-likelihood is ", round(tail(x$meta$log)[1],2), ", with last increment ", signif(tail(diff(x$meta$log))[1],5)))
  cat('\n')
  cat(ifelse(any(diff(x$meta$log)<0), "**WARNING**: negative increments in log-likelihood", "Success: no negative increments in log-likelihood"))

  cat('\n')
  cat(paste0("EM algorithm ", ifelse(x$meta$converg,"","*NOT*")," converged after ", length(x$meta$log), " steps"))
  cat('\n')

}

#
summary.po2m <- function(object, ...){
  outp <- "\n Not implemented yet \n"
  class(outp) <- "summary.po2m"
}

#
print.summary.po2m <- function(x, ...){
  cat(x)
}

