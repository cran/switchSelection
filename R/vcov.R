# Estimate asymptotic covariance matrix of the maximum-likelihood estimator
vcov_ml <- function(object, type, n_cores, n_sim)
{
  # List of output values
  out <- list(vcov = NULL)
  
  # Get some variables
  par <- object$par
  cov_type <- type
  estimator <- object$estimator
  control_lnL <- object$control_lnL
  n_par <- control_lnL$n_par

  # Estimate asymptotic covariance matrix
  H <- NULL
  H_inv <- NULL
  J <- NULL
  cov <- diag(rep(1, n_par))
  if (!is.matrix(cov_type))
  {
    if (cov_type != "no")
    {
      if ((cov_type == "sandwich") | (cov_type == "hessian"))
      {
        tryCatch(
          {
            cov_type_old <- cov_type
            cov_type <- "gop"
            H <- gena::gena.hessian(gr = ifelse(is(object = object, 
                                                   class2 = "mvoprobit"),
                                                lnL_mvoprobit,
                                                lnL_mnprobit),
                                    par = par,
                                    gr.args = list(n_sim = n_sim, 
                                                   n_cores = n_cores,
                                                   control_lnL = control_lnL, 
                                                   out_type = "grad"))
            out$H <- H
            H_inv <- qr.solve(H, tol = 1e-16)
            cov_type <- cov_type_old
          },
          error = function(e) {
            warning(paste0("Problems with numeric hessian calculation. ",
                           "Therefore 'cov_type' has been changed to 'gop'."))
          }
        )
      }
      if ((cov_type == "sandwich") | (cov_type == "gop"))
      {
        if (is(object = object, class2 = "mvoprobit"))
        {
          J <- lnL_mvoprobit(par = par,
                             n_sim = n_sim, n_cores = n_cores,
                             control_lnL = control_lnL, 
                             out_type = "jac")
        }
        if (is(object = object, class2 = "mnprobit"))
        {
          J <- lnL_mnprobit(par = par,
                            n_sim = n_sim, n_cores = n_cores,
                            control_lnL = control_lnL, 
                            out_type = "jac")
        }
        out$J <- J
      }
      if (cov_type == "sandwich")
      {
        cov <- H_inv %*% t(J) %*% J %*% H_inv
      }
      if (cov_type == "hessian")
      {
        cov <- -H_inv
      }
      
      if (any(is.na(cov)))
      {
        warning(paste0("Can't calculate covariance matrix of type '", 
                       cov_type, "'. ", 
                       "Therefore gop covariance matrix will be used instead."))
        cov_type <- "gop"
      }
      
      if (cov_type == "gop")
      {
        tryCatch(
        {
          cov <- qr.solve(t(J) %*% J, tol = 1e-16)
        },
        error = function(e) {
          warning(paste0("Problems with numeric hessian calculation. ",
                         "Therefore 'cov_type' has been changed to 'NO'."))
        }
       )
      }
    }
  } 
  
  # Aggregate some output
  out$vcov <- cov
  
  return(out)
}

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

#' Calculate Variance-Covariance Matrix for a mvoprobit Object.
#' @description Return the variance-covariance matrix of the parameters of
#' mvoprobit model.
#' @param object an object of class \code{mvoprobit}.
#' @param ... further arguments (currently ignored).
#' @param type character representing the type of the asymptotic covariance 
#' matrix estimator. If 
#' \code{estimator} argument of \code{\link[switchSelection]{mvoprobit}} is 
#' \code{"ml"} then \code{type} may be changed to any value available for input 
#' argument \code{cov_type} of \code{\link[switchSelection]{mvoprobit}}. 
#' Otherwise \code{type} will coincide with \code{cov_type} output value 
#' of \code{\link[switchSelection]{mvoprobit}}.
#' @param regime non-negative integer representing the regime of the two-step
#' procedure for which covariance matrix should be returned.
#' If \code{estimator = "2step"} and \code{regime = NULL} or 
#' \code{is.na(regime)} then covariance
#' matrix of the first step parameters' estimator will be returned.
#' Otherwise the function estimates covariance matrix for the second step
#' parameters associated with corresponding regime.
#' @param n_sim integer representing the number of GHK draws when there are
#' more than 3 ordered equations. Otherwise alternative (much more efficient) 
#' algorithms will be used to calculate multivariate normal probabilities.
#' @param n_cores positive integer representing the number of CPU cores used for 
#' parallel computing. If possible it is highly recommend to set it equal to
#' the number of available physical cores especially when the system of
#' ordered equations has 2 or 3 equations.
#' @details Argument \code{type} is closely related to the argument 
#' \code{cov_type} of \code{\link[switchSelection]{mvoprobit}} function. 
#' See 'Details' and 'Usage' sections of 
#' \code{\link[switchSelection]{mvoprobit}} 
#' for more information on \code{cov_type} argument.
#' @return Returns numeric matrix which represents estimate of the asymptotic 
#' covariance matrix of model's parameters.
vcov.mvoprobit <- function(object, ..., 
                           type = object$cov_type, 
                           regime = NULL,
                           n_cores = object$other$n_cores,
                           n_sim = object$other$n_sim)
{
  # Validate dots
  if (length(list(...)) > 0)
  {
    warning("Additional arguments passed through ... are ignored.")   
  }
  
  if (object$estimator == "ml")
  {
    if (type == object$cov_type)
    {
      return (object$cov)
    }
    return (vcov_ml(object = object, type = type, 
                    n_cores = n_cores, n_sim = n_sim)$vcov)
  }
  
  if (object$estimator == "2step")
  {
    if (type != object$cov_type)
    {
      stop(paste0("If estimator is '2step' then argument type should be the ",
                  "same as 'object$cov_type'."))
    }
    if (is.null(regime) | is.na(regime))
    {
      return (object$model1$cov)
    }
    else
    {
      if ((regime >= 0) & (regime < object$control_lnL$n_regimes[1]))
      {
        return (object$cov_2step[[regime + 1]])
      }
      else
      {
        stop("Incorrect value of the 'regime' input argument.")
      }
    }
  }
  
  stop("Incorrect 'object$estimator' value.")
}

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

#' Calculate Variance-Covariance Matrix for a mnprobit Object.
#' @description Return the variance-covariance matrix of the parameters of
#' mnprobit model.
#' @param object an object of class \code{mnprobit}.
#' @param ... further arguments (currently ignored).
#' @param type character representing the type of the asymptotic covariance 
#' matrix estimator. If 
#' @param regime non-negative integer representing the regime of the two-step
#' procedure for which covariance matrix should be returned.
#' If \code{estimator = "2step"} and \code{regime = NULL} then covariance
#' matrix of the first step parameters' estimator will be returned.
#' Otherwise the function estimates covariance matrix for the second step
#' parameters associated with corresponding regime.
#' @param n_sim integer representing the number of GHK draws when there are
#' more then 3 ordered equations. Otherwise alternative (much more efficient) 
#' algorithms will be used to calculate multivariate normal probabilities.
#' @param n_cores positive integer representing the number of CPU cores used for 
#' parallel computing. If possible it is highly recommend to set it equal to
#' the number of available physical cores especially when the system of
#' ordered equations has 2 or 3 equations.
#' \code{estimator} argument of \code{\link[switchSelection]{mnprobit}} is 
#' \code{"ml"} then \code{type} may be changed to any value available for input 
#' argument \code{cov_type} of \code{\link[switchSelection]{mnprobit}}. 
#' Otherwise \code{type} will coincide with \code{cov_type} output value 
#' of \code{\link[switchSelection]{mnprobit}}.
#' @details Argument \code{type} is closely related to the argument 
#' \code{cov_type} of \code{\link[switchSelection]{mnprobit}} function. 
#' See 'Details' and 'Usage' sections of 
#' \code{\link[switchSelection]{mnprobit}} 
#' for more information on \code{cov_type} argument.
#' @return Returns numeric matrix which represents estimate of the asymptotic 
#' covariance matrix of model's parameters.
vcov.mnprobit <- function(object, ..., 
                          type = object$cov_type,
                          regime = NULL,
                          n_cores = object$other$n_cores,
                          n_sim = object$other$n_sim)
{
  return(vcov.mvoprobit(object = object, type = type, regime = regime,
                        n_cores = n_cores, n_sim = n_sim))
}