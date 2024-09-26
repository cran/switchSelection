#' Bootstrap for msel function
#' @description Function \code{\link[switchSelection]{bootstrap_msel}} 
#' provides bootstrap estimates of the parameters of the model estimated via
#' the \code{\link[switchSelection]{msel}} function.
#' Function \code{\link[switchSelection]{bootstrap_combine_msel}} allows to
#' combine several objects of class \code{'bootstrap_msel'}.
#' @name bootstrap
#' @template bootstrap_param_Template
#' @template bootstrap_details_Template
#' @template bootstrap_return_Template
#' @template bootstrap_examples_Template
bootstrap_msel <- function(object, 
                           iter     = 100,
                           opt_type = "optim", 
                           opt_args = NULL,
                           is_ind   = FALSE,
                           n_sim    = 1000, 
                           n_cores  = 1)
{
  # The number of observations
  n <- nrow(object$data)
  
  # List to store the results
  out <- list(par = matrix(NA, nrow = iter, ncol = object$other$n_par))
  if (is_ind)
  {
    out$ind <- matrix(NA, nrow = n, ncol = iter)
  }
  out$iter <- iter
  
  # Main bootstrap routine
  for (i in 1:iter)
  {
    # Select the indexes of the observations to include into the sample
    ind <- sample(x = 1:n, size = n, replace = TRUE)
    if (is_ind)
    {
      object$ind[, i] <- ind
    }
      
    # Estimate the model
    model <- msel(formula        = object$formula, 
                  formula2       = object$formula2,
                  formula3       = object$formula3,
                  data           = object$data[ind, ], 
                  groups         = object$groups, 
                  groups2        = object$groups2,
                  groups3        = object$groups3,
                  marginal       = object$marginal,
                  opt_type       = opt_type, 
                  opt_args       = opt_args,
                  start          = object$par,
                  estimator      = object$estimator,
                  cov_type       = "no",
                  degrees        = object$degrees,
                  degrees3       = object$degrees3,
                  n_sim          = n_sim,
                  n_cores        = n_cores,
                  regularization = object$other$regularization,
                  type3          = object$type3)
      
      # Store the results
      out$par[i, ] <- model$par
  }
    
  # Assign the class
  class(out) <- "bootstrap_msel"
  
  # Covariance matrix
  out$cov  <- cov(out$par)
  
  # Return the results
  return(out)
}

#' Update msel object with the new estimates
#' @description This function updates parameters of the model estimated via 
#' \code{\link[switchSelection]{msel}} function.
#' @param object an object of class \code{'msel'}.
#' @param par a vector of parameters which substitutes \code{object$par} and
#' used to update the estimates i.e., \code{object$coef}, \code{object$cuts} and
#' others.
#' @details It may be useful to apply this function to the bootstrap
#' estimates of \code{\link[switchSelection]{bootstrap_msel}}.
#' @return This function returns an object \code{object} of class \code{'msel'}
#' in which \code{object$par} is substituted with \code{par}. Also, \code{par} 
#' is used to update the estimates i.e., \code{object$coef}, \code{object$cuts} 
#' and others.
update_msel <- function(object, par)
{
  object <- par_msel(object = object, par = par, type = "object")
  return(object)
}

#' @name bootstrap
bootstrap_combine_msel <- function(...)
{
  # Collect the models into the list
  b <- list(...)
  
  # Retrieve the first set of bootstrap results
  out <- b[[1]]
  if (length(b) == 1)
  {
    return (out)
  }
  
  # Combine the first set with other sets
  for (i in length(b))
  {
    for (j in names(b[[i]]))
    {
      if (is.list(b[[i]][[j]]))
      {
        out[[j]] <- c(out[[j]], b[[i]][[j]])
      }
    }
    out$par  <- rbind(out$par, b[[i]]$par)
    out$iter <- out$iter + b[[i]]$iter
  }
  out$cov <- cov(out$par)
  
  # Return the results
  return(out)
}