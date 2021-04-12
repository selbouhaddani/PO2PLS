#' Calculate trace of a matrix
#'
#' @param X a square matrix
#'
#' @keywords internal
#' @export
tr <- function(X){
  sum(diag(X))
}

#' Convert old PO2PLS fit to new-style PO2PLS fit
#'
#' Internal use only
#'
#' @param fit A PO2PLS of old style
#' @inheritParams PO2PLS
#'
#' @keywords internal
#' @export
PO2PLS_to_po2m <- function(fit, X, Y){
  p <- ncol(X)
  q <- ncol(Y)
  comps <- with(fit$par, c(r=ncol(W), rx=sign(ssq(Wo))*ncol(Wo),ry=sign(ssq(Co))*ncol(Co)))
  ssqX <- ssq(X)
  ssqY <- ssq(Y)
  SigU <- with(fit$par, SigT %*% B^2 + SigH)
  LVs <- E_step(X, Y, fit$par)
  row.names(LVs$mu_T) <- row.names(LVs$mu_U) <- row.names(LVs$mu_To) <- row.names(LVs$mu_Uo) <- row.names(X)
  R2s <- c(
    R2X = with(fit$par, (tr(SigT) + tr(SigTo))/(tr(SigT)+tr(SigTo)+p*sig2E)),
    R2Y = with(fit$par, (tr(SigU) + tr(SigUo))/(tr(SigU)+tr(SigUo)+q*sig2F)),
    R2Xj = with(fit$par, tr(SigT)/(tr(SigT)+tr(SigTo)+p*sig2E)),
    R2Yj = with(fit$par, tr(SigU)/(tr(SigU)+tr(SigUo)+q*sig2F)),
    R2Xhat = NA,
    R2Yhat = with(fit$par, tr(SigT %*% B^2)/tr(SigU))
  )
  loglk <- list(all_vals = fit$log, last_val = tail(fit$log,1),
                df = unname((p+1)*(comps[1]+comps[2])+(q+1)*(comps[1]+comps[3])+2*(comps[1]+1)))
  class(loglk) <- "loglik.po2m"
  outp <- list(parameters = fit$par,
    latent_vars = LVs[3:6] %>% as_tibble,
    meta_data = list(
      loglikelihood = loglk,
      explained_vars = R2s,
      comps = comps,
      time = ifelse(is.null(fit$flags$time),NA,fit$flags$time),
      call = fit$flags$call,
      convergence = fit$flags$converg,
      ssqs = c(X=ssqX, Y=ssqY)
    ))
  class(outp) <- "po2m"
  return(outp)
}

#' Prints an overview of a PO2PLS fit
#'
#' @param x A PO2PLS fit of class po2m
#' @param digits Number of decimals to output
#' @param ... For consistency
#'
#' @export
print.po2m <- function(x, digits = 3, ...){
#  cat('\n')
  cat(paste0("Call: ", x$meta_data$call %>% deparse(500L)))
  cat('\n')

  cat('\n')
  cat(paste0("This is a Probabilistic ", ifelse(sum(x$meta_data$comps[2:3])==0, "", "O2"), "PLS fit"))
  cat('\n')

#  cat('\n')
  cat(paste0("  with r=",x$meta$comps[1],", rx=",x$meta$comps[2]," and ry=", x$meta$comps[3]," components"))
  cat('\n')

#  cat('\n')
  cat(paste0("NB Data dimensions: N=",nrow(x$latent_vars[[1]]),", p=",nrow(x$param$W)," and q=", nrow(x$par$C)))
  cat('\n')

#  cat('\n')
  cat(paste0("Elapsed time: ", signif(x$meta$time,3), " sec"))
  cat('\n')

#  cat('\n')
  cat(paste0("The log-likelihood is ", round(x$meta$log$last,2), ", with last increment ", signif(tail(diff(x$meta$log$all),1),5)))
  cat('\n')
  cat(ifelse(any(diff(x$meta$log$all)<0), "**WARNING**: negative increments in log-likelihood", "NB No negative increments in log-likelihood"))

  cat('\n')
  cat(paste0("EM algorithm ", ifelse(x$meta$converg,"","*NOT*")," converged after ", length(x$meta$log$all), " steps"))
  cat('\n')

}

#' Calculates a summary of a PO2PLS fit
#'
#' @param object A PO2PLS fit of class po2m
#' @inheritParams print.po2m
#'
#' @export
summary.po2m <- function(object, digits = 3, ...){
#  outp <- "Not implemented yet \n"
  outp <- with(object, {
    list(
      Comps = meta_data$comps,
      Vars = with(meta_data, c(explained_vars, R2Ypred = explained_vars[6]/explained_vars[4])),
      Loglik = meta_data$loglikelihood,
      Time = meta_data$time,
      Digits = digits,
      Call = meta_data$call,
      Covcoeff = with(parameters, diag(SigT %*% B))
    )
  })
  class(outp) <- "summary.po2m"
  return(outp)
}

#' Prints a summary of a PO2PLS fit
#'
#' @inheritParams print.po2m
#'
#' @export
print.summary.po2m <- function(x, digits = 3, ...){
  digits <- x$Digits
  cat("==== *** Summary of fit *** ========\n\n")
  cat(paste0("This is a ** Probabilistic ", ifelse(sum(x$meta_data$comps[2:3])==0, "", "O2"), "PLS ** fit.\n"))
  cat(paste0("  Call: ", x$Call %>% deparse),'\n')
  cat("This fit took ", round(x$Time, digits), " seconds.\n\n")

  cat(" === Explained Variances in each part:\n")
  R2dat <- data.frame(X = c(Joint = x$Vars["R2Xj"],
                        Specific = x$Vars["R2X"]-x$Vars["R2Xj"],
                        Noise = 1 - x$Vars["R2X"],
                        Predictive = NA),
                  Y = c(Joint = x$Vars["R2Yj"],
                        Specific = x$Vars["R2Y"]-x$Vars["R2Yj"],
                        Noise = 1 - x$Vars["R2Y"],
                        Predictive = x$Vars["R2Yhat"]))
  row.names(R2dat) <- c("Joint", "Specific", "Noise", "Predictive")
  print(round(R2dat, digits))

  cat("\nCovariance coefficients SigmaT*B:", round(x$Covcoeff, digits), "\n\n")
  cat(" === Log-likelihood: \n")
  print.loglik.po2m(x$Loglik)

  cat()
}

#' Prints log-likelihood statistics of a PO2PLS fit
#'
#' @param x A log-likelihood object of class loglik.po2m, produced within a PO2PLS fit
#' @inheritParams print.po2m
#'
#' @export
print.loglik.po2m <- function(x, digits=3, ...){
  cat('Log-likelihood value is:', x$last,'\n')
  cat('Last increment was', tail(diff(x$all),1),'\n')
  cat("Degrees of freedom is", x$df, "\n")
  cat("  Note: Run LRT(fit, fit0) to compare two nested PO2PLS models \n")
}

#' Likelihood ratio test for two nested PO2PLS fits
#'
#' Perform a likelihood ratio test to compare two nested fits from a PO2PLS model.
#'
#' @param fit A PO2PLS fit of class po2m
#' @param fit0 A second PO2PLS fit of class po2m
#' @param digits Number of decimals to output
#' @param ... For consistency
#'
#' @return The chi-square statistics (2 times the difference in log-likelihood), the degrees of freedom and the p-value of the comparison.
#'
#' @rdname LRT
#' @export
LRT <- function(fit, fit0, digits=3, ...) UseMethod("LRT")

#' @inherit LRT
#'
#' @rdname LRT
#' @export
LRT.po2m <- function(fit, fit0, digits=3, ...){
  x <- fit$meta$log
  y <- fit0$meta$log
  if(inherits(x, "loglik.po2m") & inherits(y, "loglik.po2m")){
    cat("Difference in log-likelihood:", (2*abs(x$last - y$last)) %>% round(digits),'\n')
    cat("Difference in degrees of freedom:", abs(x$df - y$df) %>% round(digits), "\n")
    cat("Comparing to a Chi-square distribution... \n")
    cat("P-value of H0 <Smaller model equals larger model> is: \n",
        (1-pchisq(2*abs(x$last - y$last), abs(x$df - y$df))) %>% round(digits), "\n")
    return(invisible(c(statistic = 2*abs(x$last - y$last),
                       df = abs(x$df - y$df),
                       pvalue=1-pchisq(2*abs(x$last - y$last), abs(x$df - y$df)))))
  } else stop("Something went wrong, maybe fit$meta_data$loglikelihood is NULL?")
}

