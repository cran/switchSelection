regularization_validate <- function(regularization = NULL, 
                                    n_par, estimator)
{
  # -----------------------------------------------------------
  # Note that indexes are adjusted inside opt_switchSelection
  # -----------------------------------------------------------
  
  # Return if NULL
  if (is.null(regularization) | (length(regularization) == 0))
  {
    return (NULL)
  }
  
  # Stop if incorrect estimator is used
  if (estimator != "ml")
  {
    stop ("Regularization available only for maximum-likelihood estimator.")
  }
  
  # Check that regularization is a list
  if (!is.list(regularization))
  {
    stop ("Argument 'regularization' should be a list.")
  }
  
  # Validate names
  names(regularization)[which(names(regularization) == 
                        "ridge_loc")] <- "ridge_location"
  names(regularization)[which(names(regularization) == 
                        "lasso_loc")] <- "lasso_location"
  names(regularization)[which(names(regularization) == 
                        "ridge_indexes")] <- "ridge_ind"
  names(regularization)[which(names(regularization) == 
                        "lasso_indexes")] <- "lasso_ind"
    
  # Check types of regularization
  is_ridge <- hasName(x = regularization, name = "ridge_ind")
  is_lasso <- hasName(x = regularization, name = "lasso_ind")

  # Get the number of regularized parameters
  n_ridge <- 0
  n_lasso <- 0
  if (is_ridge)
  {
    n_ridge <- length(regularization$ridge_ind)
  }
  if (is_lasso)
  {
    n_lasso <- length(regularization$lasso_ind)
  }
    
  # Validate indexes
  if ((n_ridge > n_par) | (n_ridge < 0))
  {
    stop ("Invalid indexes in 'ridge_ind'.")
  }
  if ((n_lasso > n_par) | (n_lasso < 0))
  {
    stop ("Invalid indexes in 'lasso_ind'.")
  }
    
  # Validate scales
  n_ridge_scale <- 0
  n_lasso_scale <- 0
    # ridge
  if (is_ridge)
  {
    if (!hasName(x = regularization, name = "ridge_scale"))
    {
      stop ("Argument 'regularization' should have 'ridge_scale' element.")
    }
    n_ridge_scale <- length(regularization$ridge_scale)
    if (n_ridge_scale != n_ridge)
    {
      if (n_ridge_scale == 1)
      {
        regularization$ridge_scale <- rep(regularization$ridge_scale, n_ridge)
      }
      else
      {
        stop ("Invalid length of 'ridge_scale'.")
      }
    }
  }
    # lasso
  if (is_lasso)
  {
    if (!hasName(x = regularization, name = "lasso_scale"))
    {
      stop ("Argument 'regularization' should have 'lasso_scale' element.")
    }
    n_lasso_scale <- length(regularization$lasso_scale)
    if (n_lasso_scale != n_lasso)
    {
      if (n_lasso_scale == 1)
      {
        regularization$lasso_scale <- rep(regularization$lasso_scale, n_lasso)
      }
      else
      {
        stop ("Invalid length of 'lasso_scale'.")
      }
    }
  }

  # Validate locations
  n_ridge_location <- 0
  n_lasso_location <- 0
    # ridge
  if (is_ridge)
  {
    if (!hasName(x = regularization, name = "ridge_location"))
    {
      n_ridge_location <- n_ridge
      regularization$ridge_location <- rep(0, n_ridge)
    }
    else
    {
      n_ridge_location <- length(regularization$ridge_location)
      if (n_ridge_location != n_ridge)
      {
        if (n_ridge_location == 1)
        {
          regularization$ridge_location <- rep(regularization$ridge_location, 
                                               n_ridge)
        }
        else
        {
          stop ("Invalid length of 'ridge_location'.")
        }
      }
    }
  }
    # lasso
  if (is_lasso)
  {
    if (!hasName(x = regularization, name = "lasso_location"))
    {
      n_lasso_location <- n_lasso
      regularization$lasso_location <- rep(0, n_lasso)
    }
    else
    {
      n_lasso_location <- length(regularization$lasso_location)
      if (n_lasso_location != n_lasso)
      {
        if (n_lasso_location == 1)
        {
          regularization$lasso_location <- rep(regularization$lasso_location, 
                                               n_lasso)
        }
        else
        {
          stop ("Invalid length of 'lasso_location'.")
        }
      }
    }
  }

  return(regularization)
}