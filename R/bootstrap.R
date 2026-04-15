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
  # Deal with several objects
  n_objects <- NULL
  if (!is(object = object, class2 = "msel"))
  {
    # Get the number of objects
    n_objects <- length(object)
    
    # Routine for several objects
    if (n_objects > 1)
    {
      # Check the class of the object
      for (i in 1:n_objects)
      {
        if (!is(object = object[[i]], class2 = "msel"))
        {
          stop(paste0("Argument object[[", i, "]] is not of class 'msel'."))
        }
      }
    }
  }
  else
  {
    n_objects <- 1
    object <- list(object)
  }
  
  # Combine the data
  data <- object[[1]]$data
  if (n_objects > 1)
  {
    for (i in 2:n_objects)
    {
      data <- suppressMessages(dplyr::full_join(x = data, y = object[[i]]$data))
    }
  }
  
  # The number of observations
  n <- nrow(data)
  
  # Lists to store the results
  out <- vector(mode = "list", length = n_objects)
  for (i in 1:n_objects)
  {
    # Prepare the output
    out[[i]] <- list(par = matrix(NA, 
                                  nrow = iter, 
                                  ncol = object[[i]]$other$n_par))
    if (is_ind)
    {
      out[[i]]$ind <- matrix(NA, nrow = n, ncol = iter)
    }
    
    # Save the number of iterations
    out[[i]]$iter <- iter
    
    # Assign the class
    class(out[[i]]) <- "bootstrap_msel"
    class(out)      <- "crossmodel_bootstrap_msel"
  }
  
  # Main bootstrap routine
  for (i in 1:iter)
  {
    # Select the indexes of the observations to include into the sample
    ind <- sample(x = 1:n, size = n, replace = TRUE)
    
    # For each model
    for (j in 1:n_objects)
    {
      # Save the indexes if needed
      if (is_ind)
      {
        out[[j]]$ind[, i] <- ind
      }
        
      # Estimate the model
      model <- msel(formula        = object[[j]]$formula, 
                    formula2       = object[[j]]$formula2,
                    formula3       = object[[j]]$formula3,
                    data           = data[ind, ], 
                    groups         = object[[j]]$groups, 
                    groups2        = object[[j]]$groups2,
                    groups3        = object[[j]]$groups3,
                    marginal       = object[[j]]$marginal,
                    opt_type       = opt_type, 
                    opt_args       = opt_args,
                    start          = object[[j]]$par,
                    estimator      = object[[j]]$estimator,
                    cov_type       = "no",
                    degrees        = object[[j]]$degrees,
                    degrees3       = object[[j]]$degrees3,
                    n_sim          = n_sim,
                    n_cores        = n_cores,
                    regularization = object[[j]]$other$regularization,
                    type3          = object[[j]]$type3,
                    cluster        = NA)
        
        # Store the results
        out[[j]]$par[i, ] <- model$par
    }
  }
  
  # Covariance matrix
  for (i in 1:n_objects)
  {
    out[[i]]$cov  <- cov(out[[i]]$par)
  }
  
  # Return the results
  if (n_objects == 1)
  {
    out <- out[[1]]
  }
  return(out)
}

#' Update msel object with the new estimates
#' @description This function updates parameters of the model estimated via 
#' \code{\link[switchSelection]{msel}} function.
#' @param object an object of class \code{'msel'}.
#' @param par a vector of parameters which substitutes \code{object$par} and
#' is used to update the estimates i.e., \code{object$coef}, \code{object$cuts} 
#' and others.
#' @details It may be useful to apply this function to the bootstrap
#' estimates of \code{\link[switchSelection]{bootstrap_msel}}.
#' @return This function returns an object of class \code{'msel'}
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
  for (i in 2:length(b))
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