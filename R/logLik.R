#' Extract the Number of Observations from a Fit of the msel Function.
#' @description Extract the number of observations from a model fit
#' of the \code{\link[switchSelection]{msel}} function.
#' @param object object of class "msel"
#' @param ... further arguments (currently ignored)
#' @details Unobservable values of continuous equations are included into
#' the number of observations.
#' @return A single positive integer number.
nobs.msel <- function(object, ...)
{
  if (length(list(...)) > 0)
  {
    warning("Additional arguments passed through ... are ignored.")   
  }
  
  return(object$control_lnL$n_obs)
}

#' Extract Log-Likelihood from a Fit of the msel Function. 
#' @description Extract Log-Likelihood from a model fit
#' of the \code{\link[switchSelection]{msel}} function.
#' @param object object of class "msel"
#' @param ... further arguments (currently ignored)
#' @details If \code{estimator == "2step"} in 
#' \code{\link[switchSelection]{msel}} then function may return
#' \code{NA} value since two-step estimator of covariance matrix may be
#' not positively defined.
#' @return Returns an object of class 'logLik'.
logLik.msel <- function (object, ...)
{
  if (length(list(...)) > 0)
  {
    warning("Additional arguments passed through ... are ignored.")   
  }
  
  lnL <- object$logLik
  attr(lnL, "class") <- "logLik"
  attr(lnL, "df") <- length(as.vector(object$par))
  
  return(lnL)
}