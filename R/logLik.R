#' Extract the Number of Observations from a Fit of the mvoprobit Function.
#' @description Extract the number of observations from a model fit
#' of the \code{\link[switchSelection]{mvoprobit}} function.
#' @param object object of class "mvoprobit"
#' @param ... further arguments (currently ignored)
#' @details Unobservable values of continuous equations are included into
#' the number of observations.
#' @return A single positive integer number.
nobs.mvoprobit <- function(object, ...)
{
  if (length(list(...)) > 0)
  {
    warning("Additional arguments passed through ... are ignored.")   
  }
  
  return(object$control_lnL$n_obs)
}

#' Extract Log-Likelihood from a Fit of the mvoprobit Function. 
#' @description Extract Log-Likelihood from a model fit
#' of the \code{\link[switchSelection]{mvoprobit}} function.
#' @param object object of class "mvoprobit"
#' @param ... further arguments (currently ignored)
#' @details If \code{estimator == "2step"} in 
#' \code{\link[switchSelection]{mvoprobit}} then function may return
#' \code{NA} value since two-step estimator of covariance matrix may be
#' not positively defined.
#' @return Returns an object of class 'logLik'.
logLik.mvoprobit <- function (object, ...)
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

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

#' Extract the Number of Observations from a Fit of the mnprobit Function.
#' @description Extract the number of observations from a model fit
#' of the \code{\link[switchSelection]{mnprobit}} function.
#' @param object object of class "mnprobit"
#' @param ... further arguments (currently ignored)
#' @details Unobservable values of continuous equations are included into
#' the number of observations.
#' @return A single positive integer number.
nobs.mnprobit <- function(object, ...)
{
  if (length(list(...)) > 0)
  {
    warning("Additional arguments passed through ... are ignored.")   
  }
  
  return(object$control_lnL$n_obs)
}

#' Extract Log-Likelihood from a Fit of the mnprobit Function. 
#' @description Extract Log-Likelihood from a model fit
#' of the \code{\link[switchSelection]{mnprobit}} function.
#' @param object object of class "mnprobit"
#' @param ... further arguments (currently ignored)
#' @details If \code{estimator == "2step"} in 
#' \code{\link[switchSelection]{mnprobit}} then function may return
#' \code{NA} value since two-step estimator of covariance matrix may be
#' not positively defined.
#' @return Returns an object of class 'logLik'.
logLik.mnprobit <- function(object, ...)
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