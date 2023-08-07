#' Extract Residual Standard Deviation 'Sigma'
#' @description Extract standard deviations of random errors of continuous
#' equations of \code{\link[switchSelection]{mvoprobit}} function.
#' @param object object of class "mvoprobit".
#' @param use.fallback logical, passed to \code{nobs} (currently ignored).
#' @param ... further arguments (currently ignored).
#' @param eq2 index of continuous equation
#' @param regime regime of continuous equation
#' @details Available only if \code{estimator = "ml"} or all \code{degrees}
#' values are equal to \code{1}.
#' @returns Returns estimates of the standard deviations 
#' of \eqn{\varepsilon_{i}}. 
#' If \code{eq2 = k} then estimates only for \eqn{k}-th continuous equation are 
#' returned. If in addition \code{regime = r} then estimate of 
#' \eqn{\sqrt{Var(\varepsilon_{ri})}} is returned. 
#' Herewith if \code{regime} is not \code{NULL} and \code{eq2} is \code{NULL} 
#' it is assumed that \code{eq2 = 1}.
sigma.mvoprobit <- function (object, use.fallback = TRUE, ..., 
                             regime = NULL, eq2 = NULL)
{
  # Validation
  if (object$estimator != "ml")
  {
    if (object$cov_type != "parametric")
    {
      stop (paste0("Availabe only for maximum-likelihood estimator or if ",
                   "parametric estimator of the asymptotic covariance matrix ",
                   "is used."))
    }
  }
  
  # Estimate variance
  val <- coef(object, type = "var", eq2 = eq2, regime = regime)

  # Get standard deviation
  if (is.list(val))
  {
    val <- lapply(val, function(x){sqrt(x)})
  }
  else
  {
    val <- sqrt(val)
  }
  
  # Return the result
  return(val)
}

#' Extract Residual Standard Deviation 'Sigma'
#' @description Extract standard deviations of random errors of continuous
#' equation of \code{\link[switchSelection]{mnprobit}} function.
#' @param object object of class "mnprobit".
#' @param use.fallback logical, passed to \code{nobs} (currently ignored).
#' @param ... further arguments (currently ignored).
#' @param regime regime of continuous equation
#' @details Available only if \code{estimator = "ml"}.
#' @returns Returns the estimates of standard deviations 
#' of \eqn{\varepsilon_{i}}.
#' If \code{regime = r} then estimate of 
#' \eqn{\sqrt{Var(\varepsilon_{ri})}} is returned.
sigma.mnprobit <- function (object, use.fallback = TRUE, ..., regime = NULL)
{
  # Validation
  if (object$estimator != "ml")
  {
    stop (paste0("Availabe only for maximum-likelihood estimator."))
  }
  
  # Estimate variance
  val <- sqrt(coef(object, type = "var", regime = regime))
  
  # Get standard deviation
  if (is.list(val))
  {
    val <- lapply(val, function(x){sqrt(x)})
  }
  else
  {
    val <- sqrt(val)
  }
  
  # Return the result
  return(val)
}