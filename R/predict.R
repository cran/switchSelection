#' Predict method for mvoprobit function
#' @description Predicted values based on object of class 'mvoprobit'.
#' @template predict.mvoprobit_param_Template
#' @template predict.mvoprobit_details_Template
#' @template predict.mvoprobit_return_Template
predict.mvoprobit <- function(object, ..., 
                              newdata = NULL, 
                              given_ind = numeric(),
                              group = NULL,
                              group2 = NULL,
                              type = ifelse(is.null(group2), "prob", "val"),
                              me = NULL,
                              eps = NULL,
                              control = list(),
                              se = FALSE)
{
  # Special routine for standard errors
  if (se)
  {
    # Get arguments of the function and change se
    # argument to false to prevent eternal cycle
    predict_args <- c(as.list(environment()), list(...))
    predict_args$se <- FALSE
    
    # Apply delta method
    predict_fn <- function(object)
    {
      predict_args$object <- object
      do.call(what = predict.mvoprobit, args = predict_args)
    }
    out <- delta_method(object = object, fn = predict_fn)
    
    return(out)
  }
  
  # Some values
  n_par <- length(object$par)
  n_eq <- object$control_lnL$n_eq
  n_eq2 <- object$control_lnL$n_eq2
  n_groups <- object$control_lnL$n_groups
  
  # Validation
  
  # group
  if (!is.null(group))
  {
    if (any((group %% 1) != 0))
    {
      stop(paste0("Invalid 'group' argument. Please, insure that 'group' is ",
                  "a vector of integers."))
    }
    if (length(group) != n_eq)
    {
      stop(paste0("Invalid 'group' argument. Please, insure that length ",
                  "of 'group' equals to the number of ordered equations i.e. ",
                  "'length(group) == length(formula)'."))
    }
    if (any(group < -1))
    {
      stop(paste0("Invalid 'group' argument. Please, insure that it ", 
                  "does not contain any negative values other than -1."))
    }
    if (any(group >= n_groups))
    {
      stop(paste0("Invalid 'group' argument. Please, insure that it ", 
                  "does not contain any values greater than the number of ", 
                  "the maximum category of corresponding equation."))
    }
  }
  
  # Get information on marginal distributions
  marginal <- object$marginal
  is_marginal <- length(object$marginal) > 0
  
  # Transform marginal argument into mnorm format
  marginal_mnorm <- object$marginal
  if (is_marginal)
  {
    for (i in 1:n_eq)
    {
      if (length(object$marginal_par[[i]]) > 0)
      {
        marginal_mnorm[[i]] <- object$marginal_par[[i]]
      }
    }
  }
  
  # Special routine to differentiate lambda respect to all arguments
  if (is.null(group) & is.null(group2) & 
      (type == "dlambda") & is.null(newdata))
  {
    delta <- sqrt(.Machine$double.eps) * 10
    lambda <- predict.mvoprobit(object = object, type = "lambda")
    dlambda <- array(dim = c(nrow(lambda), ncol(lambda), n_par))
    # respect to coefficients
    coef_ind <- object$control_lnL$coef_ind
    for (i in 1:n_eq)
    {
      for (j in 1:length(object$coef[[i]]))
      {
        par_old_tmp <- object$coef[[i]][j] 
        delta_tmp <- delta * abs(par_old_tmp)
        object$coef[[i]][j] <- par_old_tmp + delta_tmp
        lambda_tmp <- predict.mvoprobit(object = object, type = "lambda")
        dlambda[, , coef_ind[[i]][j] + 1] <- (lambda_tmp - lambda) / 
          delta_tmp
        object$coef[[i]][j] <- par_old_tmp
      }
    }
    # respect to cuts
    cuts_ind <- object$control_lnL$cuts_ind
    for (i in 1:n_eq)
    {
      for (j in 1:length(object$cuts[[i]]))
      {
        par_old_tmp <- object$cuts[[i]][j]
        delta_tmp <- delta * abs(par_old_tmp)
        object$cuts[[i]][j] <- par_old_tmp + delta_tmp
        lambda_tmp <- predict.mvoprobit(object = object, type = "lambda")
        dlambda[, , cuts_ind[[i]][j] + 1] <- (lambda_tmp - lambda) /
          delta_tmp
        object$cuts[[i]][j] <- par_old_tmp
      }
    }
    # respect to coefficients of variance equation
    coef_var_ind <- object$control_lnL$coef_var_ind
    is_het <- object$other$is_het
    for (i in 1:n_eq)
    {
      if (is_het[i])
      {
        for (j in 1:length(coef_var_ind[[i]]))
        {
          par_old_tmp <- object$coef_var[[i]][j]
          delta_tmp <- delta * abs(par_old_tmp)
          object$coef_var[[i]][j] <- par_old_tmp + delta_tmp
          lambda_tmp <- predict.mvoprobit(object = object, type = "lambda")
          dlambda[, , coef_var_ind[[i]][j] + 1] <- (lambda_tmp - lambda) / 
            delta_tmp
          object$coef_var[[i]][j] <- par_old_tmp
        }
      }
    }
    # respect to sigma
    if (n_eq > 1)
    {
      sigma_ind_mat <- object$control_lnL$sigma_ind_mat
      for (i in 1:(n_eq - 1))
      {
        for (j in (i + 1):n_eq)
        {
          par_old_tmp <- object$sigma[i, j]
          delta_tmp <- delta * abs(par_old_tmp)
          object$sigma[i, j] <- par_old_tmp + delta_tmp
          object$sigma[j, i] <- object$sigma[i, j]
          lambda_tmp <- predict.mvoprobit(object = object, type = "lambda")
          dlambda[, , sigma_ind_mat[i, j] + 1] <- (lambda_tmp - lambda) / 
            delta_tmp
          object$sigma[i, j] <- par_old_tmp
          object$sigma[j, i] <- par_old_tmp
        }
      }
    }
    # respect to parameters of marginal distribution
    if (is_marginal)
    {
      marginal_par_n <- object$control_lnL$marginal_par_n
      marginal_par_ind <- object$control_lnL$marginal_par_ind
      for (i in 1:n_eq)
      {
        if (marginal_par_n[i] > 0)
        {
          for (j in 1:marginal_par_n[i])
          {
            par_old_tmp <- object$marginal_par[[i]][j]
            delta_tmp <- delta * abs(par_old_tmp)
            object$marginal_par[[i]][j] <- par_old_tmp + delta_tmp
            lambda_tmp <- predict.mvoprobit(object = object, type = "lambda")
            dlambda[, , marginal_par_ind[[i]][j] + 1] <- (lambda_tmp - lambda) / 
              delta_tmp
            object$marginal_par[[i]][j] <- par_old_tmp
          }
        }
      }
    }
    return(dlambda)
  }
  
  # Special routine for calculation of lambda for each observation
  if (is.null(group) & is.null(group2) & 
      (type == "lambda") & is.null(newdata))
  {
    data <- object$data
    n <- nrow(data)
    lambda <- NULL
    if (!hasName(control, "adj2"))
    {
      lambda <- matrix(0, nrow = n, ncol = object$control_lnL$n_eq)
    }
    else
    {
      lambda <- array(dim = c(object$control_lnL$n_eq, 
                              object$control_lnL$n_eq, n))
      lambda[, ,] <- 0
    }
    groups <- object$groups
    ind_g <- lapply(object$control_lnL$ind_g, function(x){x + 1})
    ind_eq <- lapply(object$control_lnL$ind_eq, function(x){x + 1})
    for (i in 1:length(ind_g))
    {
      if (!hasName(control, "adj2"))
      {
        lambda[ind_g[[i]], ind_eq[[i]]] <- predict(object, 
                                                   newdata = data[ind_g[[i]], ],
                                                   group = groups[i, ],
                                                   type = "lambda",
                                                   control = control)
      }
      else
      {
        lambda[ind_eq[[i]], 
               ind_eq[[i]], 
               ind_g[[i]]] <- predict(object, 
                                      newdata = data[ind_g[[i]], ],
                                      group = groups[i, ],
                                      type = "lambda",
                                      control = control)
      }
    }
    return(lambda)
  }
  
  # For continuous variables
  if (!is.null(group2))
  {
    # Get data
    if (is.null(newdata))
    {
      W_mean <- as.matrix(object$W_mean)
      W_var <- as.matrix(object$W_var)
      X <- object$X
    }
    else
    {
      if (is.data.frame(newdata))
      {
        newdata_tmp <- switchSelection::mvoprobit(formula = object$formula, 
                                                  formula2 = object$formula2,
                                                  data = newdata, 
                                                  control = list(out_type = "data"))
      }
      W_mean <- as.matrix(newdata_tmp$W_mean)
      W_var <- as.matrix(newdata_tmp$W_var)
      X <- newdata_tmp$X
    }
    
    # Prepare matrix to store predictions
    y_pred <- matrix(NA, nrow = nrow(W_mean[[1]]), ncol = length(group2))
  }
  
  # Calculation of marginal effects
  if (!is.null(me))
  {
    # If several marginal effects should be calculated
    n_me <- length(me)
    if (n_me > 1)
    {
      list_return <- list()
      for (i in 1:n_me)
      {
        list_return[[me[i]]] <- predict(object, ..., newdata = newdata,
                                        given_ind = given_ind, group = group,
                                        type = type, me = me[i], eps = eps,
                                        control = control)
      }
      return(list_return)
    }
    # Prepare data
    if (is.null(newdata))
    {
      newdata <- object$data
    }
    else
    {
      newdata <- switchSelection::mvoprobit(formula = object$formula, 
                                            formula2 = object$formula2,
                                            data = newdata, 
                                            control = list(out_type = "data"))
      newdata <- newdata$data
    }
    
    # Determine the type of marginal effect
    is_discrete <- length(eps) > 1
    
    # Adjust epsilon if need
    if (is.null(eps))
    {
      eps <- newdata[[me]] * sqrt(.Machine$double.eps)
    }
    
    # Calculate values before and after
    if (is_discrete)
    {
      newdata[me] <- eps[1]
    }
    p0 <- predict.mvoprobit(object, 
                            newdata = newdata,
                            given_ind = given_ind,
                            group = group,
                            group2 = group2,
                            type = type,
                            me = NULL,
                            control = control)
    if (is_discrete)
    {
      newdata[me] <- eps[2]
    }
    else
    {
      newdata[me] <- newdata[me] + eps
    }
    p1 <- predict.mvoprobit(object, 
                            newdata = newdata,
                            given_ind = given_ind,
                            group = group,
                            group2 = group2,
                            type = type,
                            me = NULL,
                            control = control)
    
    # Estimate marginal effect
    val <- p1 - p0
    if (!is_discrete)
    {
      val <- val / eps
    }
    return(val)
  }
  
  # Prediction for continuous equation
  if (!is.null(group2))
  {
    # Get indexes of equations that are not omitted
    if (!is.null(group))
    {
      ind_eq <- which(group != -1)
    }
    
    # Part for conditional prediction if need
    if (!is.null(group))
    {
      lambda <- predict(object, ..., newdata = newdata,
                        given_ind = given_ind, group = group,
                        type = "lambda",
                        control = control)
    }
    
    # Unconditional predictions
    for (i in 1:length(group2))
    {
      if (group2[i] != -1)
      {
        X_i <- as.matrix(X[[i]])
        y_pred[, i] <- X_i %*% matrix(object$coef2[[i]][group2[i] + 1, ], 
                                      ncol = 1)
        # conditional predictions (on ordered equations)
        if (!is.null(group))
        {
          if (object$estimator == "ml")
          {
            y_pred[, i] <- y_pred[, i] + lambda %*% 
              matrix(object$cov2[[i]][group2[i] + 1, ind_eq], 
                     ncol = 1)
          }
          else
          {
            degrees <- object$degrees[group2 + 1, ]
            for (j in 1:n_eq)
            {
              lambda_pow <- 1
              if (degrees[j] >= 1)
              {
                for (j1 in 1:degrees[j])
                {
                  lambda_pow <- lambda_pow * lambda[, j]
                  y_pred[, i] <- y_pred[, i] + 
                    lambda_pow * 
                    object$coef_lambda[[group2 + 1]][[j]][j1]
                }
              }
            }
          }
        }
      }
    }
    
    # Return prediction for continuous equation
    if (type == "val")
    {
      return(y_pred[, group2 != -1])
    }
  }
  
  # Assign group value for some types of return value
  if (((type == "li") | (type == "sd")) & is.null(group))
  {
    group <- rep(1, object$control_lnL$n_eq)
  }
  
  # Store the data
  W_mean <- NULL
  W_var <- NULL
  if (is.null(newdata))
  {
    W_mean <- as.matrix(object$W_mean)
    W_var <- as.matrix(object$W_var)
  }
  else
  {
    if (is.data.frame(newdata))
    {
      newdata <- switchSelection::mvoprobit(formula = object$formula, 
                                            formula2 = object$formula2,
                                            data = newdata, 
                                            control = list(out_type = "data"))
    }
    W_mean <- as.matrix(newdata$W_mean)
    W_var <- as.matrix(newdata$W_var)
  }
  n_obs <- nrow(W_mean[[1]])
  
  # Return marginal probabilities for each equation
  # if there is no particular group
  if (is.null(group) & (type == "prob"))
  {
    n_obs <- object$control_lnL$n_obs
    n_eq <- object$control_lnL$n_eq
    prob <- matrix(NA, nrow = n_obs, ncol = n_eq)
    colnames(prob) <- paste0("P(y", 1:n_eq, "=1)")
    group <- rep(-1, n_eq)
    
    for (i in 1:n_eq)
    {
      group[i] <- 1
      prob[, i] <- predict(object, 
                           newdata = newdata, 
                           type = "prob",
                           group = group,
                           control = control)
      group[i] <- -1
    }
    
    warning(paste0("Since 'group' argument is not specified, marginal ",
                   "probabilities of '1' will be estimated for each ",
                   "ordered equation and 'type' will be set to 'prob'."))
    
    return(prob)
  }
  
  # Some additional variables
  is_eq <- group != -1
  ind_eq <- which(is_eq)
  ind_eq_omit <- which(!is_eq)
  n_eq_g <- length(ind_eq)
  
  # Adjust for marginal distributions of observable equations only
  marginal_mnorm_g <- marginal_mnorm
  if (is_marginal)
  {
    marginal <- marginal[ind_eq]
    marginal_mnorm_g <- marginal_mnorm[ind_eq]
  }
  
  # Validation
  if (length(given_ind) >= n_eq_g)
  {
    stop(paste("At least one component should be unconditioned.",
               "Please, insure that 'length(given_ind)' is smaller than",
               "the number of observable equations."))
  }
  
  # Calculate conditional probabilities if need
  if (length(given_ind) > 0)
  {
    if (type != "prob")
    {
      warning(paste0("Since 'given_ind' has been provided then 'type'",
                     "will be coerced to 'prob'."))
    }
    
    group_given <- group
    group_given[-given_ind] <- -1
    
    p_intersection <- predict(object, 
                              newdata = newdata, 
                              type = "prob",
                              group = group,
                              control = control)
    p_given <- predict(object, 
                       newdata = newdata, 
                       type = "prob",
                       group = group_given,
                       control = control)
    
    return(p_intersection / p_given)
  }
  
  # Get some variables from the model
  coef <- object$coef
  coef_var <- object$coef_var
  sigma <- object$sigma
  cuts <- object$cuts
  control_lnL <- object$control_lnL
  n_eq <- control_lnL$n_eq
  n_cuts_eq <- control_lnL$n_cuts_eq
  is_het <- control_lnL$is_het
  
  # Calculate linear index
  li_lower <- matrix(NA, nrow = n_obs, ncol = n_eq_g)
  li_upper <- matrix(NA, nrow = n_obs, ncol = n_eq_g)
  li_mean <- matrix(NA, nrow = n_obs, ncol = n_eq_g)
  li_var <- matrix(1, nrow = n_obs, ncol = n_eq_g)
  for (i in 1:n_eq_g)
  {
    # Calculate linear indexes for mean and variance parts
    li_mean[, i] <- W_mean[[ind_eq[i]]] %*% coef[[ind_eq[i]]]
    if (is_het[ind_eq[i]])
    {
      li_var[, i] <- exp(W_var[[ind_eq[i]]] %*% coef_var[[ind_eq[i]]])
    }
    # Adjust lower limits for cuts and heteroscedasticity
    if (group[ind_eq[i]] != 0)
    {
      li_lower[, i] <- cuts[[ind_eq[i]]][group[ind_eq[i]]] - 
        li_mean[, i]
      if (is_het[ind_eq[i]])
      {
        li_lower[, i] <- li_lower[, i, drop = FALSE] / 
          li_var[, i, drop = FALSE]
      }
    }
    else
    {
      li_lower[, i] <- -Inf
    }
    # Adjust upper limits for cuts and heteroscedasticity
    if (group[ind_eq[i]] != n_cuts_eq[[ind_eq[i]]])
    {
      li_upper[, i] <- cuts[[ind_eq[i]]][group[ind_eq[i]] + 1] - 
        li_mean[, i, drop = FALSE]
      if (is_het[ind_eq[i]])
      {
        li_upper[, i] <- li_upper[, i, drop = FALSE] / 
          li_var[, i, drop = FALSE]
      }
    }
    else
    {
      li_upper[, i] <- Inf
    }
  }
  
  # Return linear indexes if need
  if (type == "li")
  {
    return(li_mean)
  }
  
  # Return standard deviations if need
  if (type == "sd")
  {
    return(li_var)
  }
  
  # Get covariance matrix
  sigma_g <- sigma[ind_eq, ind_eq, drop = FALSE]
  if (type == "prob")
  {
    # Calculate probabilities
    prob <- mnorm::pmnorm(lower = li_lower, upper = li_upper, 
                          mean = rep(0, n_eq_g), sigma = sigma_g,
                          marginal = marginal_mnorm_g)$prob
    return(prob)
  }
  
  if (type == "lambda")
  {
    # Calculate lambdas for conditional expectations
    grads <- mnorm::pmnorm(lower = li_lower, upper = li_upper, 
                           mean = rep(0, n_eq_g), sigma = sigma_g,
                           log = TRUE, 
                           grad_upper = TRUE, 
                           grad_lower = TRUE,
                           grad_sigma = hasName(control, name = "adj2"),
                           marginal = marginal_mnorm_g,
                           grad_marginal = is_marginal,
                           grad_marginal_prob = is_marginal)
    lambda_lower <- NULL
    lambda_upper <- NULL
    if (!is_marginal)
    {
      lambda_lower <- -grads$grad_lower
      lambda_upper <- grads$grad_upper
    }
    else
    {
      lambda_lower <- -grads$grad_lower_marginal
      lambda_upper <- grads$grad_upper_marginal
    }
    if (hasName(control, name = "adj"))
    {
      if (control$adj)
      {
        li_lower_adj <- li_lower
        li_lower_adj[is.infinite(li_lower)] <- 0
        li_upper_adj <- li_upper
        li_upper_adj[is.infinite(li_upper)] <- 0
        lambda <- lambda_lower * li_lower_adj - lambda_upper * li_upper_adj
      }
      else
      {
        stop("Parameter 'adj' in 'control' should be 'TRUE'.")
      }
    }
    else
    {
      lambda <- lambda_lower - lambda_upper
    }
    if (hasName(control, name = "adj2"))
    {
      lambda2 <- grads$grad_sigma
      return(lambda2)
    }
    
    return(lambda)
  }
  
  stop("Wrong 'type' argument.")
}

#' Predict method for mnprobit function
#' @template predict.mnprobit_param_Template
#' @template predict.mnprobit_details_Template
#' @template predict.mnprobit_return_Template
predict.mnprobit <- function(object, ..., 
                             newdata = NULL,
                             alt = 1, 
                             regime = -1,
                             type = ifelse(is.null(regime) | (regime == -1), 
                                           "prob", "val"),
                             alt_obs = "all",
                             me = NULL, eps = NULL,
                             control = list(),
                             se = FALSE)
{   
  # Special routine for standard errors
  if (se)
  {
    # Get arguments of the function and change se
    # argument to false to prevent eternal cycle
    predict_args <- c(as.list(environment()), list(...))
    predict_args$se <- FALSE
    
    # Apply delta method
    predict_fn <- function(object)
    {
      predict_args$object <- object
      do.call(what = predict.mnprobit, args = predict_args)
    }
    out <- delta_method(object = object, fn = predict_fn)
    
    return(out)
  }
  
  # Get some variables
  coef <- object$coef
  sigma <- object$sigma
  n_alt <- object$control_lnL$n_alt
  n_obs <- object$control_lnL$n_obs
  n_par <- length(object$par)
  n_coef <- object$control_lnL$n_coef
  n_regimes <- object$n_regimes
  is2 <- n_regimes > 0
  
  if ((type == "prob") & is.null(alt))
  {
    probs <- matrix(NA, ncol = n_alt, nrow = n_obs)
    for (i in 1:n_alt)
    {
      probs[, i] <- predict(object, type = "prob", 
                            alt = i, newdata = newdata)
    }
    colnames(probs) <- paste0("P(z=", 1:n_alt, ")")
    rownames(probs) <- 1:n_obs
    return(probs)
  }
  
  # Special routine to differentiate lambda respect to all arguments
  if (is.null(alt) & (regime == -1) & 
      (type == "dlambda") & is.null(newdata))
  {
    delta <- sqrt(.Machine$double.eps) * 10
    lambda <- predict.mnprobit(object = object, type = "lambda", alt = NULL)
    dlambda <- array(dim = c(nrow(lambda), ncol(lambda), n_par))
    # respect to coefficients
    coef_ind_alt <- object$control_lnL$coef_ind_alt
    for (i in 1:(n_alt - 1))
    {
      for (j in 1:n_coef)
      {
        par_old_tmp <- object$coef[j, i] 
        delta_tmp <- delta * abs(par_old_tmp)
        object$coef[j, i]  <- par_old_tmp + delta_tmp
        lambda_tmp <- predict.mnprobit(object = object, 
                                       type = "lambda", alt = NULL)
        dlambda[, , coef_ind_alt[j, i] + 1] <- (lambda_tmp - lambda) / 
          delta_tmp
        object$coef[j, i]  <- par_old_tmp
      }
    }
    # respect to sigma
    if (n_alt > 2)
    {
      sigma_ind <- object$control_lnL$sigma_ind
      counter <- 1
      for (i in 1:(n_alt - 1))
      {
        for (j in 1:i)
        {
          if (!((i == 1) & (j == 1)))
          {
            par_old_tmp <- object$sigma[i, j]
            delta_tmp <- delta * abs(par_old_tmp)
            object$sigma[i, j] <- par_old_tmp + delta_tmp
            object$sigma[j, i] <- object$sigma[i, j]
            lambda_tmp <- predict.mnprobit(object = object, 
                                           type = "lambda", alt = NULL)
            dlambda[, , sigma_ind[counter] + 1] <- (lambda_tmp - lambda) /
              delta_tmp
            counter <- counter + 1
            object$sigma[i, j] <- par_old_tmp
            object$sigma[j, i] <- par_old_tmp
          }
        }
      }
    }
    return(dlambda)
  }
  
  # Calculate all alternative specific lambdas
  if (is.null(alt) & (type == "lambda"))
  {
    lambda <- matrix(NA, nrow = n_obs, ncol = n_alt - 1)
    for (i in 1:n_alt)
    {
      ind_alt <- object$control_lnL$ind_alt[[i]] + 1
      lambda[ind_alt, ] <- predict(object, alt = i, 
                                   type = "lambda", alt_obs = "alt")
    }
    return(lambda)
  }
  
  # Calculation of marginal effects
  if (!is.null(me))
  {
    # If several marginal effects should be calculated
    n_me <- length(me)
    if (n_me > 1)
    {
      list_return <- list()
      for (i in 1:n_me)
      {
        list_return[[me[i]]] <- predict(object, newdata = newdata,
                                        alt = alt, regime = regime,
                                        alt_obs = alt_obs,
                                        type = type, me = me[i], eps = eps)
      }
      return(list_return)
    }
    # Prepare data
    if (is.null(newdata))
    {
      newdata <- object$data
    }
    else
    {
      newdata <- switchSelection::mnprobit(formula = object$formula, 
                                           formula2 = object$formula2,
                                           data = newdata, 
                                           control = list(out_type = "data"))
      newdata <- newdata$data
    }
    
    # Determine the type of marginal effect
    is_discrete <- length(eps) > 1
    
    # Adjust epsilon if need
    if (is.null(eps))
    {
      eps <- newdata[[me]] * sqrt(.Machine$double.eps)
    }
    
    # Calculate values before and after
    if (is_discrete)
    {
      newdata[me] <- eps[1]
    }
    p0 <- predict.mnprobit(object, newdata = newdata,
                           alt = alt, regime = regime,
                           alt_obs = alt_obs,
                           type = type, me = NULL, eps = NULL)
    if (is_discrete)
    {
      newdata[me] <- eps[2]
    }
    else
    {
      newdata[me] <- newdata[me] + eps
    }
    p1 <- predict.mnprobit(object, ..., newdata = newdata,
                           alt = alt, regime = regime,
                           alt_obs = alt_obs,
                           type = type, me = NULL, eps = NULL)
    
    # Estimate marginal effect
    val <- p1 - p0
    if (!is_discrete)
    {
      val <- val / eps
    }
    return(val)
  }
  
  # Store the data
  if (!is.null(newdata))
  {
    df <- model.frame(object$formula, data = newdata, na.action = na.omit)
    W <- as.matrix(cbind(1, df[, -1]))
    z <- as.vector(df[, 1])
    df2 <- NULL
    X <- NULL
    y <- NULL
    if (is2)
    {
      df2 <- model.frame(object$formula2, data = newdata, na.action = na.omit)
      X <- as.matrix(cbind(1, df2[, -1]))
      y <- as.vector(df2[, 1])
    }
    n_obs <- nrow(W)
  }
  else
  {
    W <- as.matrix(object$W)
    z <- as.vector(object$z)
    X <- NULL
    y <- NULL
    if (is2)
    {
      X <- as.matrix(object$X)
      y <- as.vector(object$y)
    }
  }
  
  # # Slightly transform sigma_alt if need
  # if (exists("cor_bounds", where = control))
  # {
  #   sigma_cor <- cov2cor(sigma)
  #   for (i in 2:(n_alt - 1))
  #   {
  #     for (j in 1:(i - 1))
  #     {
  #       if (sigma[i, j] > control$cor_bounds[2])
  #       {
  #         sigma_cor[i, j] <- control$cor_bounds[2]
  #         sigma_cor[j, i] <- sigma_cor[i, j]
  #       }
  #       if (sigma_cor[i, j] < control$cor_bounds[1])
  #       {
  #         sigma_cor[i, j] <- control$cor_bounds[1]
  #         sigma_cor[j, i] <- sigma_cor[i, j]
  #       }
  #     }
  #   }
  #   sigma_sd <- diag(sqrt(diag(sigma)))
  #   sigma <- sigma_sd %*% sigma_cor %*% sigma_sd
  # }
  
  # Select appropriate observations
  if (alt_obs != "all")
  {
    if (alt_obs == "alt")
    {
      alt_obs <- alt
    }
    alt_ind = which(z %in% alt_obs)
    n_obs <- length(alt_ind)
    W <- W[alt_ind, ]
    z <- z[alt_ind]
  }
  
  # Unconditional predictions of continuous variable
  coef2 <- NULL
  y_pred <- NULL
  coef_lambda <- NULL
  if (is2 & (type == "val"))
  {
    coef_lambda <- object$coef_lambda
    coef2 <- object$coef2
    if (regime != -1)
    {
      y_pred <- X %*% matrix(coef2[, regime + 1], ncol = 1)
    }
    if (is.null(alt))
    {
      return(y_pred)
    }
  }
  
  # Calculate linear index
  li <- matrix(NA, nrow = n_obs, ncol = n_alt - 1)
  for (i in 1:(n_alt - 1))
  {
    li[, i] <- W %*% coef[, i]
  }
  
  # Return linear indexes if enough
  if (type == "li")
  {
    return(li[, alt])
  }
  
  # Create transformation matrix
  transform_mat <- matrix(0, nrow = n_alt, ncol = n_alt)
  diag(transform_mat) <- 1
  if (alt != n_alt)
  {
    transform_mat[, alt] <- -1
  }
  transform_mat <- transform_mat[-alt, , drop = FALSE]
  transform_mat <- transform_mat[, -n_alt, drop = FALSE]
  
  # Construct covariance matrix for alternative
  sigma_alt <- transform_mat %*% sigma %*% t(transform_mat)
  
  # Calculate the differences between linear indexes
  li_diff <- matrix(NA, n_obs, n_alt - 1);
  if (alt == n_alt)
  {
    li_diff <- -li;
  }
  else
  {
    li_diff[, n_alt - 1] <- li[, alt];
    li_diff[, -(n_alt - 1)] <- apply(li[, -alt, drop = FALSE], 2, 
                                     function(x)
                                     {
                                       return(li[, alt] - x)
                                     }
    ) 
    # counter <- 0
    # for (j in 1:(n_alt - 1))
    # {
    #   if (alt != j)
    #   {
    #     li_diff[, counter] = li[, alt] - li[, j]
    #     counter <- counter + 1
    #   }
    # }
  }
  
  # Calculate the probabilities
  lower_neg_inf <- matrix(-Inf, nrow = n_obs, ncol = n_alt - 1)
  prob = mnorm::pmnorm(lower = lower_neg_inf, upper = li_diff,
                       mean = rep(0, n_alt - 1), sigma = sigma_alt,
                       is_validation = FALSE)$prob
  
  if (type == "prob")
  {
    return(prob)
  }
  
  # Calculate lambdas
  lambda <- mnorm::pmnorm(lower = lower_neg_inf, upper = li_diff,
                          mean = rep(0, n_alt - 1), sigma = sigma_alt,
                          grad_upper = TRUE,
                          log = TRUE,
                          is_validation = FALSE)$grad_upper
  A <- matrix(0, nrow = n_alt - 1, ncol = n_alt - 1)
  for (i in 1:(n_alt - 1))
  {
    for (j in 1:(n_alt - 1))
    {
      if (i == alt)
      {
        A[i, j] <- 1
      }
      if (((i < alt) & (i == j)) | 
          ((i > alt) & (i == (j + 1))))
      {
        A[i, j] <- -1
      }
    }
  }
  lambda <- lambda %*% t(A)
  
  # return the results
  if (type == "lambda")
  {
    return(lambda)
  }
  
  # Conditional predictions
  if (is2 & (regime != -1) & (type == "val"))
  {
    # conditional predictions (on continuous equations)
    if (!is.null(alt))
    {
      degrees <- object$degrees[regime + 1, ]
      for (j in 1:(n_alt - 1))
      {
        lambda_pow <- 1
        if (degrees[j] >= 1)
        {
          for (j1 in 1:degrees[j])
          {
            lambda_pow <- lambda_pow * lambda[, j]
            y_pred <- y_pred + lambda_pow * coef_lambda[[regime + 1]][[j]][j1]
          }
        }
      }
    }
  }
  
  # Return prediction for continuous equation
  if (type == "val")
  {
    return(y_pred)
  }
  
  # If everything should be returned
  out <- list(li = li,
              li_diff = li_diff,
              prob = prob,
              lambda = lambda)
  
  return("No return value has been specified")
}