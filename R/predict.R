#' Predict method for msel function
#' @description Predicted values based on the object of class 'msel'.
#' @template predict.msel_param_Template
#' @template predict.msel_details_Template
#' @template predict.msel_return_Template
predict.msel <- function(object, ..., 
                         newdata   = NULL, 
                         given_ind = numeric(),
                         group     = NA,
                         group2    = NA,
                         group3    = NA,
                         type      = ifelse(any(is.na(group2)), 
                                            "prob", "val"),
                         me        = NULL,
                         eps       = NULL,
                         control   = list(),
                         test      = FALSE,
                         exogenous = NULL)
{
  # -------------------------------------------------------
  # Deal with the data and the variables
  # -------------------------------------------------------

  # Deal with the groups
  is_na_group  <- c(any(is.na(group)), any(is.na(group2)), any(is.na(group3)))
  
  # Get some values
  groups    <- object$groups
  groups2   <- object$groups2
  groups3   <- object$groups3
  formula   <- object$formula
  formula2  <- object$formula2
  formula3  <- object$formula3
  ind_g     <- object$other$ind_g
  ind_eq    <- object$other$ind_eq
  n_par     <- length(object$par)
  n_eq      <- object$other$n_eq
  n_eq2     <- object$other$n_eq2
  n_eq3     <- object$other$n_eq3
  is1       <- object$other$is1
  is2       <- object$other$is2
  is3       <- object$other$is3
  n_groups  <- object$other$n_groups
  estimator <- object$estimator
  type3     <- object$type3
  
  # Provide 'data' as 'newdata' if need
  is_newdata <- TRUE
  if (is.null(newdata))
  {
    is_newdata <- FALSE
    newdata    <- object$data
  }
  
  # Add the intercept to the newdata
  if (is2 | is3)
  {
    newdata[, '(Intercept)'] <- 1
  }

  # Deal with the 'exogenous' argument
  n_exogenous  <- length(exogenous)
  is_exogenous <- n_exogenous > 0
  if (is_exogenous)
  {
    newdata    <- exogenous_fn(exogenous = exogenous, newdata = newdata)
    is_newdata <- TRUE
  }

  # Recalculate selectivity terms for the two-step estimator
  if ((type == "val") & (estimator == "2step") & 
      ((is1 & !is_na_group[1]) | (is3 & !is_na_group[3])))
  {
    # Set newdata
    is_newdata <- TRUE
    
    # Remove the data associated with the selectivity terms
    for (v in 1:n_eq2)
    {
      formula2_terms <- attr(terms(formula2[[v]]), "term.labels")
      formula2_terms <- formula2_terms[formula2_terms %in% colnames(newdata)]
      newdata[, formula2_terms[grepl("lambda", 
                                     formula2_terms, 
                                     fixed = TRUE)]] <- NULL
      if (is1)
      {
        for (j in seq_len(n_eq))
        {
          newdata[[paste0("lambda", j)]] <- 0
        }
      }
      if (is3)
      {
        for (j in seq_len(n_eq3 - (type3 == "probit")))
        {
          newdata[[paste0("lambda", j, "_mn")]] <- 0
        }
      }
    }
    
    # Calculate the selectivity terms of the ordinal equations
    lambda <- NULL
    if (is1)
    {
      lambda <- predict(object, group = group, group3 = -1, 
                        type = "lambda", newdata = newdata)
      for (j in seq_len(n_eq))
      {
        if (group[j] != -1)
        {
          newdata[[paste0("lambda", j)]] <- lambda[, j]
        }
      }
    }
    
    # Calculate the selectivity terms of the multinomial equation
    lambda_mn <- NULL
    if (is3)
    {
      lambda_mn <- predict(object, group = group, group3 = group3, 
                           type = "lambda_mn", newdata = newdata)
      for (j in seq_len(n_eq3 - (type3 == "probit")))
      {
        if (group3 != -1)
        {
          newdata[[paste0("lambda", j, "_mn")]] <- lambda_mn[, j]
        }
      }
    }
  }

  # Adjust the new data
  if (is_newdata)
  {
    # Deal with the omitted variables
    all_names                <- object$other$all_names
    lambda_names             <- object$other$lambda_names
    omitted_names            <- all_names[!(all_names %in% colnames(newdata))]
    if (estimator == "2step")
    {
      omitted_names <- omitted_names[!(omitted_names %in% lambda_names)]
    }
    newdata[, omitted_names] <- 0

    # Remove the missing values
    newdata <- complete_msel(object = object, data = newdata)
  }

  # Split the newdata according to the formulas
  if (!is_newdata)
  {
    W_mean <- object$W_mean
    W_var  <- object$W_var
    W_mn   <- object$W_mn
    X      <- object$X
    y      <- object$y
  }
  else
  {
    newdata_list <- data_msel(object = object, data = newdata)
    W_mean       <- newdata_list$W_mean
    W_var        <- newdata_list$W_var
    X            <- newdata_list$X
    W_mn         <- newdata_list$W_mn
    z            <- newdata_list$z
    y            <- newdata_list$y
    z_mn         <- newdata_list$z_mn
    newdata_list <- NULL
  }
  
  # Get the number of observations
  n_obs <- nrow(newdata)
  
  # Get the information on the marginal distributions
  marginal    <- object$marginal
  is_marginal <- length(object$marginal) > 0
  
  # Get estimates of the parameters and some related variables
  coef        <- object$coef
  coef_var    <- object$coef_var
  coef3       <- object$coef3
  sigma       <- object$sigma
  sigma3      <- object$sigma3
  cuts        <- object$cuts
  sigma_omit  <- object$other$sigma_omit
  n_eq        <- object$other$n_eq
  n_cuts_eq   <- object$other$n_cuts_eq
  is_het      <- object$other$is_het
  
  # Estimate the value
  if (all(is_na_group))
  {
    # Get some variables
    groups   <- object$groups
    groups2  <- object$groups2
    groups3  <- object$groups3
    n_groups <- object$other$n_groups
    ind_g    <- object$other$ind_g
    
    # Change the groups related data if need
    if (is_newdata)
    {
      groups_list <- groups_msel(object  = object, data    = newdata, 
                                 groups  = groups, groups2 = groups2,
                                 groups3 = groups3)
      groups      <- groups_list$groups
      groups2     <- groups_list$groups2
      groups3     <- groups_list$groups3
      n_groups    <- groups_list$n_groups
      ind_g       <- groups_list$ind_g
    }
    
    # Calculate the value for all groups and
    # aggregate the result
    val <- NULL
    if (n_groups > 1)
    {
      for (i in 1:n_groups)
      {
        groups_i  <- NA
        groups2_i <- NA
        groups3_i <- NA
        if (is1)
        {
          groups_i  <- groups[i, ]
        }
        if (is2)
        {
          groups2_i <- groups2[i, ]
        }
        if (is3)
        {
          groups3_i <- groups3[i]
        }
        val_tmp <- predict(object, 
                           newdata = newdata[ind_g[[i]], ],
                           group   = groups_i,
                           group2  = groups2_i,
                           group3  = groups3_i,
                           type    = type)
        if (i == 1)
        {
          val           <- matrix(0, nrow = n_obs, ncol = ncol(val_tmp))
          colnames(val) <- colnames(val_tmp)
        }
        val[ind_g[[i]], ] <- val_tmp
      }
    }
    
    # Return the results
    return(val)
  }
  
  # -------------------------------------------------------
  # Statistical test
  # -------------------------------------------------------
  
  # Deal with the input
  test_fn <- NULL
  if (!is.logical(test))
  {
    test_fn <- test
    test    <- TRUE
  }
  
  # Perform the test
  if (test)
  {
    # Get arguments of the function and change 'test'
    # argument to false to prevent the eternal cycle
    predict_args         <- c(as.list(environment()), list(...))
    predict_args$newdata <- newdata
    predict_args$test    <- FALSE

    # Apply the test
    predict_fn <- function(object)
    {
      predict_args$object <- object
      val <- do.call(what = predict.msel, args = predict_args)
      if (!is.null(test_fn))
      {
        if (is.function(test_fn))
        {
          val <- test_fn(val)
        }
        else
        {
          val <- val[, test_fn]
        }
      }
      return(val)
    }
    out <- test_msel(object = object, fn = predict_fn)
    
    return(out)
  }
  
  # -------------------------------------------------------
  # Validation
  # -------------------------------------------------------
  
  # Validate type
  if (type == "lp")
  {
    type <- "li"
  }
  if (type == "lp_mn")
  {
    type <- "li_mn"
  }
  
  # Validate group
  if (!is_na_group[1])
  {
    if (any((group %% 1) != 0))
    {
      stop(paste0("Invalid 'group' argument. Please, insure that 'group' is ",
                  "a vector of integers.\n"))
    }
    if (length(group) != n_eq)
    {
      stop(paste0("Invalid 'group' argument. Please, insure that length ",
                  "of 'group' equals to the number of ordered equations i.e. ",
                  "'length(group) == length(formula)'.\n"))
    }
    if (any(group < -1))
    {
      stop(paste0("Invalid 'group' argument. Please, insure that it ", 
                  "does not contain any negative values other than -1.\n"))
    }
    if (any(group >= n_groups))
    {
      stop(paste0("Invalid 'group' argument. Please, insure that it ", 
                  "does not contain any values greater than the number of ", 
                  "the maximum category of corresponding equation.\n"))
    }
  }
  
  # Validate group2
  if (is2 & (type == "val") & !is_na_group[1] & is_na_group[2])
  {
    group2 <- rep(0, n_eq2)
    warning(paste0("It is assumed that 'group2' is a vector of zeros.\n"))
  }
  if (is2 & !is_na_group[2] & (length(group2) != n_eq2))
  {
    stop(paste0("Invalid 'group2' argument. ",
                "It should be a vector of length ", n_eq2, ".\n"))
  }
  
  # Validate group3
  if (!is_na_group[3])
  {
    if (((group3 %% 1) != 0) & length(group3) > 1)
    {
      stop(paste0("Invalid 'group3' argument. Please, insure that 'group3' is ",
                  "an integer.\n"))
    }
    if (any(group3 < -1))
    {
      stop(paste0("Invalid 'group3' argument. Please, insure that it ", 
                  "is not a negative value different from -1.\n"))
    }
    if (group3 > (n_eq3 - 1))
    {
      stop(paste0("Invalid 'group3' argument. Please, insure that its ", 
                  "value is not larget than the number of the alternatives ",
                  "in the multinomial equation.\n"))
    }
  }
  
  # Additional validation for groups and groups3
  if (is1 & is3)
  {
    if (is_na_group[1])
    {
      group          <- rep(-1, n_eq)
      is_na_group[1] <- FALSE
    }
    if (is_na_group[3])
    {
      group3         <- -1
      is_na_group[3] <- FALSE
    }
  }
  
  # Validate control
  is_scores <- FALSE
  if (hasName(control, "is_scores"))
  {
    is_scores <- control$is_scores
  }
  
  # Assign 'group' for some 'type'
  if (((type == "li") | (type == "sd")) & is_na_group[1])
  {
    group          <- rep(0, n_eq)
    is_na_group[1] <- FALSE
  }
  if ((type == "li_mn") & is_na_group[3])
  {
    group3         <- 0
    is_na_group[3] <- FALSE
  }
  
  # Some additional variables
  is_eq       <- group != -1
  ind_eq      <- which(is_eq)
  ind_eq_omit <- which(!is_eq)
  n_eq_g      <- length(ind_eq)
  y_names     <- object$other$y_names
  
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
  
  # Adjust for the marginal distributions of the observable equations only
  marginal_mnorm_g <- marginal_mnorm
  if (is_marginal)
  {
    marginal         <- marginal[ind_eq]
    marginal_mnorm_g <- marginal_mnorm[ind_eq]
  }
  
  # Covariance matrix for the multinomial probit
  sigma3_alt    <- NULL
  transform_mat <- NULL
  if (is3 & (type3 == "probit"))
  {
    transform_mat <- matrix(0, nrow = n_eq3, ncol = n_eq3)
    diag(transform_mat) <- 1
    if (group3 != (n_eq3 - 1))
    {
      transform_mat[, group3 + 1] <- -1
    }
    transform_mat <- transform_mat[-(group3 + 1), , drop = FALSE]
    transform_mat <- transform_mat[, -n_eq3, drop = FALSE]
    sigma3_alt    <- transform_mat %*% sigma3 %*% t(transform_mat)
  }
  
  # -------------------------------------------------------
  # Conditional probabilities
  # -------------------------------------------------------
  
  # Calculate the conditional probabilities
  if (length(given_ind) > 0)
  {
    if (length(given_ind) >= n_eq_g)
    {
      stop(paste("At least one component should be unconditioned.",
                 "Please, insure that 'length(given_ind)' is smaller than",
                 "the number of observable equations."))
    }
    
    if (type != "prob")
    {
      warning(paste0("Since 'given_ind' has been provided then 'type'",
                     "will be coerced to 'prob'."))
    }
    
    group_given             <- group
    group_given[-given_ind] <- -1
    
    # Probability of the intersection
    p_intersection <- predict(object, 
                              newdata = newdata, 
                              type    = "prob",
                              group   = group,
                              group3  = group3,
                              control = control)
    
    # Probability of the condition
    p_given <- predict(object, 
                       newdata = newdata, 
                       type    = "prob",
                       group   = group_given,
                       group3  = group3,
                       control = control)
    
    # Conditional probability
    p_cond <- p_intersection / p_given
    
    # Provide the name for the conditional probability
    is_eq_ng            <- is_eq
    is_eq_ng[given_ind] <- FALSE
    colnames(p_cond)[1] <- paste0("P(",
                                  paste0(object$other$z_names[is_eq_ng], "=",
                                         group[is_eq_ng], collapse = ", "),
                                  "|",
                                  paste0(object$other$z_names[given_ind], "=",
                                         group[given_ind], collapse = ", "),
                                  ")")
    rownames(p_cond)   <- 1:length(p_cond)
    
    # Return the results
    return(p_cond)
  }
  
  # -------------------------------------------------------
  # Marginal effects
  # -------------------------------------------------------
  
  # Calculate marginal effects
  if (!is.null(me))
  {
    # If several marginal effects should be calculated
    n_me <- length(me)
    if (n_me > 1)
    {
      list_return <- list()
      for (i in 1:n_me)
      {
        list_return[[me[i]]] <- predict(object, ..., 
                                        newdata   = newdata,
                                        given_ind = given_ind, 
                                        group     = group,
                                        group3    = group3,
                                        type      = type, 
                                        me        = me[i], 
                                        eps       = eps,
                                        control   = control)
      }
      return(list_return)
    }
    
    # Determine the type of the marginal effect
    is_discrete <- length(eps) > 1
    
    # Adjust epsilon if need
    if (is.null(eps))
    {
      eps <- newdata[[me]] * sqrt(.Machine$double.eps)
    }
    
    # Calculate the values before and after the increment
    if (is_discrete)
    {
      newdata[me] <- eps[1]
    }
    p0 <- predict.msel(object, 
                            newdata   = newdata,
                            given_ind = given_ind,
                            group     = group,
                            group2    = group2,
                            group3    = group3,
                            type      = type,
                            me        = NULL,
                            control   = control)
    if (is_discrete)
    {
      newdata[me] <- eps[2]
    }
    else
    {
      newdata[me] <- newdata[me] + eps
    }
    p1 <- predict.msel(object, 
                       newdata   = newdata,
                       given_ind = given_ind,
                       group     = group,
                       group2    = group2,
                       group3    = group3,
                       type      = type,
                       me        = NULL,
                       control   = control)
    
    # Estimate marginal effect
    val <- p1 - p0
    if (!is_discrete)
    {
      val <- val / eps
    }
    
    # Provide a name for the output
    if (is.matrix(val))
    { 
      colnames(val) <- paste0("d", colnames(val), "/", "d", me)
      if (is_discrete)
      {
        colnames(val) <- paste0(colnames(val), " from ", 
                                eps[1], " to ", eps[2])
      }
    }
    
    # Return marginal effect
    return(val)
  }

  # -------------------------------------------------------
  # Linear predictors (indexes)
  # -------------------------------------------------------
  
  # Calculate the linear predictors (indexes) of the ordinal equations
  li_lower <- NULL
  li_upper <- NULL
  li_mean  <- NULL
  li_var   <- NULL
  if (is1)
  {
    if (n_eq_g > 0)
    {
      li_lower <- matrix(NA, nrow = n_obs, ncol = n_eq_g)
      li_upper <- matrix(NA, nrow = n_obs, ncol = n_eq_g)
      li_mean  <- matrix(NA, nrow = n_obs, ncol = n_eq_g)
      li_var   <- matrix(1,  nrow = n_obs, ncol = n_eq_g)
      for (i in 1:n_eq_g)
      {
        # Calculate the linear predictors (indexes) for the 
        # mean and variance parts
        li_mean[, i] <- W_mean[[ind_eq[i]]] %*% coef[[ind_eq[i]]]
        if (is_het[ind_eq[i]])
        {
          li_var[, i] <- exp(W_var[[ind_eq[i]]] %*% coef_var[[ind_eq[i]]])
        }
        # Adjust lower limits for the cuts and heteroscedasticity
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
        # Adjust upper limits for the cuts and heteroscedasticity
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
    }
    
    # Return the linear predictors (indexes)
    if (type == "li")
    {
      colnames(li_mean) <- object$other$z_names[is_eq]
      rownames(li_mean) <- 1:nrow(li_mean)
      return(li_mean)
    }
    
    # Return the standard deviations
    if (type == "sd")
    {
      colnames(li_var) <- object$other$z_names[is_eq]
      rownames(li_var) <- 1:nrow(li_var)
      return(li_var)
    }
  }
  
  # Calculate the linear predictors (indexes) of the multinomial equation
  li_mn   <- NULL
  li_diff <- NULL
  if (is3 & (group3 != -1))
  {
    # Linear predictors (indexes)
    li_mn <- matrix(0, nrow = n_obs, ncol = n_eq3)
    if (type %in% c("li_mn", "prob_mn", "lambda_mn"))
    {
      for (i in seq_len(n_eq3 - 1))
      {
        li_mn[, i] <- W_mn %*% t(coef3[i, , drop = FALSE])
      }
    }
    # Differences between the linear predictors (indexes)
    if (type3 == "probit")
    {
      li_mn <- li_mn[, -n_eq3, drop = FALSE]
      li_diff <- matrix(NA, n_obs, n_eq3 - 1)
      if ((group3 + 1) == n_eq3)
      {
        li_diff <- -li_mn
      }
      else
      {
        li_diff[, n_eq3 - 1]    <- li_mn[, group3 + 1, drop = FALSE]
        li_diff[, -(n_eq3 - 1)] <- apply(li_mn[, -(group3 + 1), 
                                                 drop = FALSE], 2, 
                                         function(x)
                                         {
                                           return(li_mn[, (group3 + 1)] - x)
                                         })
      }
    }
    
    # Return the linear predictors (indexes)
    if (type == "li_mn")
    {
      colnames(li_mn) <- paste0("alt", 0:(n_eq3 - 1))
      rownames(li_mn) <- 1:nrow(li_mn)
      return(li_mn)
    }
    
    # Return the differences
    if (type == "li_diff")
    {
      return(li_diff)
    }
  }
  
  # -------------------------------------------------------
  # Probabilities
  # -------------------------------------------------------
  
  # Covariances between random errors of the observable
  # ordinal equations
  sigma_g <- sigma
  if (is1)
  {
    sigma_g <- sigma[ind_eq, ind_eq, drop = FALSE]
    if (any(sigma_omit[ind_eq, ind_eq] == 1))
    {
      sigma_g[sigma_omit[ind_eq, ind_eq]] <- 0
      warning("Unidentified covariances are set to 0.")
    }
  }
  
  # Calculate the probabilities of the ordered equations
  prob <- 0
  if (is1 & (type == "prob") & (n_eq_g > 0))
  {
    # Estimate the probability
    prob <- mnorm::pmnorm(lower    = li_lower, 
                          upper    = li_upper, 
                          mean     = rep(0, n_eq_g), 
                          sigma    = sigma_g,
                          marginal = marginal_mnorm_g)$prob
    
    # Provide the name for the probability
    colnames(prob)[1] <- paste0("P(",
                                paste0(object$other$z_names[is_eq], "=",
                                       group[is_eq], collapse = ", "),
                                ")")
    rownames(prob)    <- 1:nrow(prob)
    
    # Return the result
    return(prob)
  }
  
  # Calculate the probabilities of the multinomial equation
  prob_mn <- 0
  if (is3 & (group3 != -1) &
      ((type == "prob_mn") | (type == "lambda_mn")))
  {
    # Multinomial logit
    if (type3 == "logit")
    {
      li_exp    <- exp(li_mn)
      prob_mn <- sweep(li_exp, 1, rowSums(li_exp), FUN = '/')
    }
    # Multinomial probit
    if (type3 == "probit")
    {
      lower_neg_inf <- matrix(-Inf, nrow = n_obs, ncol = n_eq3 - 1)
      prob_mn       <- mnorm::pmnorm(lower         = lower_neg_inf, 
                                     upper         = li_diff,
                                     mean          = rep(0, n_eq3 - 1), 
                                     sigma         = sigma3_alt,
                                     is_validation = FALSE)$prob
    }
  }
  
  # Return the result
  if  (type == "prob_mn")
  {
    if (type3 == "logit")
    {
      prob_mn <- prob_mn[, (group3 + 1), drop = FALSE]
    }
    colnames(prob_mn) <- paste0("P(", object$other$z_mn_names, 
                                " = ", group3, ")")
    rownames(prob_mn) <- 1:nrow(prob_mn)
    return(prob_mn)
  }
  
  # -------------------------------------------------------
  # Lambda
  # -------------------------------------------------------

  # Calculate lambda of the ordered equations
  lambda <- NULL
  if (is1)
  {
    lambda <- matrix(0, nrow = n_obs, ncol = n_eq_g)
    if (n_eq_g > 0 & ((type == "lambda") | 
                      ((type == "val") & (estimator == "ml"))))
    {
      # Prepare the matrix to store the final results
      lambda_mat <- matrix(0, nrow = n_obs, ncol = n_eq)
      colnames(lambda_mat) <- paste0("lambda", 1:n_eq)
      if (n_eq_g > 0)
      {
        # Calculate lambdas for the conditional expectations
        grads <- mnorm::pmnorm(lower              = li_lower, 
                               upper              = li_upper, 
                               mean               = rep(0, n_eq_g), 
                               sigma              = sigma_g,
                               log                = TRUE, 
                               grad_upper         = TRUE, 
                               grad_lower         = TRUE,
                               grad_sigma         = FALSE,
                               marginal           = marginal_mnorm_g,
                               grad_marginal      = is_marginal,
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
            li_lower_adj                        <- li_lower
            li_lower_adj[is.infinite(li_lower)] <- 0
            li_upper_adj                        <- li_upper
            li_upper_adj[is.infinite(li_upper)] <- 0
            lambda                              <- lambda_lower * li_lower_adj - 
                                                   lambda_upper * li_upper_adj
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
      }
      
      # Return lambda
      if (type == "lambda")
      {
        lambda_mat[, ind_eq] <- lambda
        return(lambda_mat)
      }
    }
    else
    {
      if (type == "lambda")
      {
        return(matrix(NA, nrow = n_obs, ncol = n_eq))
      }
    }
  }
  
  # Calculate lambda for the multinomial equation
  lambda_mn <- NULL
  if (is3 & (group3 == -1) & (type == "lambda_mn"))
  {
    if (type3 == "logit")
    {
      lambda_mn <- matrix(0, nrow = n_obs, ncol = n_eq3)
    }
    if (type3 == "probit")
    {
      lambda_mn <- matrix(0, nrow = n_obs, ncol = n_eq3 - 1)
    }
    return(lambda_mn)
  }
  if (is3 & (group3 != -1) & (type == "lambda_mn"))
  {
    lambda_mn <- matrix(0, nrow = n_obs, ncol = n_eq3)
    # Multinomial logit
    if (type3 == "logit")
    {
      lambda_mn                  <- matrix(NA, nrow = n_obs, ncol = n_eq3)
      lambda_mn[, group3 + 1]    <- -log(prob_mn[, group3 + 1])
      lambda_mn[, -(group3 + 1)] <- prob_mn[, -(group3 + 1)] * 
                                      log(prob_mn[, -(group3 + 1)]) / 
                                      (1 - prob_mn[, -(group3 + 1)])
    }
    # Multinomial probit
    if (type3 == "probit")
    {
      lower_neg_inf <- matrix(-Inf, nrow = n_obs, ncol = n_eq3 - 1)
      lambda_mn     <- mnorm::pmnorm(lower         = lower_neg_inf, 
                                     upper         = li_diff,
                                     mean          = rep(0, n_eq3 - 1), 
                                     sigma         = sigma3_alt,
                                     grad_upper    = TRUE,
                                     log           = TRUE, 
                                     is_validation = FALSE)$grad_upper
      A <- matrix(0, nrow = n_eq3 - 1, ncol = n_eq3 - 1)
      for (i in 1:(n_eq3 - 1))
      {
        for (j in 1:(n_eq3 - 1))
        {
          if (i == (group3 + 1))
          {
            A[i, j] <- 1
          }
          if (((i < (group3 + 1)) & (i == j)) | 
              ((i > (group3 + 1)) & (i == (j + 1))))
          {
            A[i, j] <- -1
          }
        }
      }
      lambda_mn <- lambda_mn %*% t(A)
    }
    
    # Return the result
    if (type == "lambda_mn")
    {
      return(lambda_mn)
    }
  }
  
  # -------------------------------------------------------
  # Predictions for the continuous equations
  # -------------------------------------------------------
  
  # Prediction for the continuous equations
  y_pred    <- NULL
  scores    <- vector(mode = "list", length = n_eq2)
  if (!is_na_group[2])
  {
    # Matrix to store the predictions
    y_pred           <- matrix(NA, nrow = n_obs, ncol = n_eq2)
    colnames(y_pred) <- object$other$y_names
    
    # Predictions
    for (i in 1:n_eq2)
    {
      X_i <- as.matrix(X[[i]])
      if (group2[i] != -1)
      {
        # substitute data with the predictions
        if (i >= 2)
        {
          # predictions of the previous equation
          for (v0 in 1:(i - 1))
          {
            if (group2[v0] != -1)
            {
              if (y_names[v0] %in% colnames(X_i))
              {
                X_i[, y_names[v0]] <- y_pred[, v0]
              }
            }
          }
        }
        # predictions
        y_pred[, i] <- X_i %*% matrix(object$coef2[[i]][group2[i] + 1, ], 
                                      ncol = 1)
        if (!is_na_group[1])
        {
          if (object$estimator == "ml")
          {
            y_pred[, i] <- y_pred[, i] + lambda %*% 
                           matrix(object$cov2[[i]][group2[i] + 1, ind_eq], 
                                  ncol = 1)
          }
        }
      }
      # provide the names
      is1_cond <- FALSE
      if (!is_na_group[1])
      {
        if (any(is_eq))
        {
          is1_cond <- TRUE
        }
      }
      is3_cond <- FALSE
      if (!is_na_group[3])
      {
        if (group3 != -1)
        {
          is3_cond <- TRUE
        }
      }
      E_y_names  <- paste0("E(", object$other$y_names)
      if (!is_na_group[2])
      {
        group2_tmp                   <- group2
        group2_tmp[group2_tmp == -1] <- ""
        E_y_names                    <- paste0(E_y_names, group2_tmp)
      }
      cond_names <- ""
      if (is1_cond | is3_cond)
      {
        E_y_names <- paste0(E_y_names, "|")
        if (is1_cond)
        {
          cond_names <- paste0(cond_names, 
                               paste0(object$other$z_names[is_eq], "=", 
                                      group[is_eq], collapse = ", "))
          if (is3_cond)
          {
            cond_names <- paste0(cond_names, ",")
          }
        }
        if (is3_cond)
        {
          cond_names <- paste0(object$other$z_mn_names, "=", 
                               group3, collapse = ", ")
        }
      }
      E_y_names           <- paste0(paste0(E_y_names, cond_names), ")")
      colnames(y_pred)    <- E_y_names
      rownames(y_pred)    <- 1:nrow(y_pred)
      # calculate the scores
      if (is_scores)
      {
        scores[[i]] <- as.vector(y[, i] - y_pred[, i]) * X_i
        scores[[i]][is.infinite(scores[[i]]) | is.na(scores[[i]])] <- 0
      }
    }
  }
    
  # If the scores should be returned
  if (is_scores)
  {
    names(scores) <- object$other$y_names
    return(scores)
  }
    
  # Return prediction for the continuous equation
  if (type == "val")
  {
    return(y_pred)
  }
  
  # -------------------------------------------------------
  # Normally user should not pass to this section
  # -------------------------------------------------------

  stop("Wrong 'type' argument.")
}