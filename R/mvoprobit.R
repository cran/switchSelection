#' Multivariate ordered probit model with heteroscedasticity 
#' and (non-random) sample selection.
#' @description This function allows to estimate parameters of multivariate
#' ordered probit model and it's extenstions. 
#' It is possible to account for heteroscedastic
#' variances, non-normal marginal distributions of random errors 
#' (under Gaussian copula) and (non-random) sample selection i.e. when some 
#' categories of 
#' particular dependent variables are observable only under some specific values 
#' of other dependent variables. Also it is possible to include continuous
#' equations to get multivariate generalization of endogenous switching model.
#' In this case both maximum-likelihood and two-step 
#' (similar to Heckman's method) estimation procedures are implemented.
#' @template mvoprobit_param_Template
#' @template mvoprobit_details_Template
#' @template mvoprobit_return_Template
#' @template mvoprobit_references_Template
#' @template mvoprobit_examples1_Template
#' @template mvoprobit_examples2_Template
#' @template mvoprobit_examples3_Template
#' @template mvoprobit_examples4_Template
#' @template mvoprobit_examples5_Template
#' @template mvoprobit_examples6_Template
#' @template mvoprobit_examples7_Template
#'
mvoprobit <- function(formula,
                      formula2 = NULL,
                      data = NULL,
                      groups = NULL,
                      groups2 = NULL,
                      marginal = list(),
                      opt_type = "optim",
                      opt_args = NULL,
                      start = NULL,
                      estimator = "ml",
                      cov_type = ifelse(estimator == "ml", 
                                        "sandwich", "parametric"),
                      degrees = NULL,
                      n_sim = 1000, 
                      n_cores = 1,
                      control = NULL)
{
  # Get first step model if it has been provided
  model1 <- NULL
  if (as.character(class(formula)) == "mvoprobit")
  {
    model1 <- formula
    formula <- model1$formula
  }
  
  # List to store the output
  out <- list()
  class(out) <- "mvoprobit"
  
  # Validation
  
    # data
  if (!is.data.frame(data))
  {
    stop(paste0("Argument 'data' is missing or wrong. ",
                "Please, insure that 'data' is dataframe ",
                "containing variables described in 'formula'",
                ifelse(is.null(formula2), "", " and 'formula2'"),
                ".\n"))
  }
  
      # estimator
  estimator <- tolower(estimator)
  if (estimator %in% c("two-step", "twostep","2-step", "2st", "2"))
  {
    warning(paste0("It is assumed that estimator '", estimator,
                   "' is '2step'."))
    estimator <- "2step"
  }
  if (estimator %in% c("mle", "likelihood","maximum-likelihood", "1step"))
  {
    warning(paste0("It is assumed that estimator '", estimator,
                   "' is 'ml'."))
    estimator <- "ml"
  }
  if (!(estimator %in% c("2step", "ml")))
  {
    stop(paste0("Wrong 'estimator' argument. ",
                "Please, insure that it is either '2step' or 'ml'."))
  }
  if ((estimator == "2step") & (is.list(formula2) & (length(formula2) != 1)))
  {
    stop(paste0("Currently '2step' estimator available only ",
                "when there is a single continuous equation. ",
                "Please, insure that 'length(formula2) == 1'."))
  }
  
    # groups
  if (!is.null(groups))
  {
    if (!is.matrix(groups))
    {
      groups <- as.matrix(groups)
    }
    if (ncol(groups) != length(formula))
    {
      stop(paste0("Argument 'groups' is wrong. ",
                  "Please, insure that 'groups' is a matrix ",
                  "which number of columns equals to the length ",
                  "of 'formula'.",
                  "\n"))
    }
  }
  
    # opt_type
  opt_type_vec <- c("optim", "gena", "pso")
  if (!(opt_type %in% opt_type_vec))
  {
    stop(paste0("Argument 'opt_type' is wrong. ",
                "Please, insure that it is one of: ",
                paste0(opt_type_vec, collapse = ", "),
                "\n"))
  }
  
    # cov_type
  cov_type <- tolower(cov_type)
  cov_type_vec <- NULL
  if (estimator == "ml")
  {
    cov_type_vec <- c("sandwich", "hessian", "gop", "no")
  }
  else
  {
    cov_type_vec <- c("parametric", "nonparametric")
  }
  if (!is.matrix(cov_type))
  {
    if (!(cov_type %in% cov_type_vec))
    {
      stop(paste0("Argument 'cov_type' is wrong. ",
                  "Please, insure that it is one of: ",
                  paste0(cov_type_vec, collapse = ", "),
                  "\n"))
    }
  }
  
    # formula
  if (!is.list(formula))
  {
    tryCatch(formula <- as.formula(formula),
             error = function(e) {stop("Invalid 'formula' argument.")})
    formula <- list(formula)
  }
  
    # formula2
  is2 <- FALSE
  if (!is.null(formula2))
  {
    is2 <- TRUE
    if (!is.list(formula2))
    {
      tryCatch(formula2 <- as.formula(formula2),
               error = function(e) {stop("Invalid 'formula2' argument.")})
      formula2 <- list(formula2)
    }
  }
  
    # marginal
  n_marginal <- length(marginal)
  if (n_marginal > 0)
  {
    # transform vector of characters into the list of NULL values
    if (!is.list(marginal))
    {
      if (is.character(marginal))
      {
        marginal_list <- vector(mode = "list", length = n_marginal)
        names(marginal_list) <- marginal
        marginal <- marginal_list
      }
      else
      {
        stop("Argument 'marginal' should be either character vector or a list.")
      }
    }
    if (n_marginal != length(formula))
    {
      stop(paste0("The number of elements in 'marginal' should be the same ",
                  "as the number of ordered equations. ",
                  "Please, insure that 'length(marginal) == length(formula)' ",
                  "or 'marginal' is an empty character vector."))
    }
  }
  
    # degrees
  if (is.null(degrees))
  {
    degrees <- rep(1, length(formula))
  }
  if (any((degrees %% 1) != 0) | any(degrees < 0))
  {
    stop(paste0("Invalid 'degrees' argument. Please, insure that ", 
                "it is a vector of non-negative integers."))
  }
  if (!is.null(formula2) & !is.null(degrees))
  {
    if ((length(degrees) != length(formula)) & !is.matrix(degrees))
    {
      stop(paste0("Number of elements in 'degrees' should be equal to the ",
                  "number of ordered equations. ",
                  "Please, insure that 'length(degrees) == length(formula)'"))
    }
    if (!is.numeric(degrees) & !is.null(degrees))
    {
      degrees <- as.numeric(degrees)
    }
  }
  if ((cov_type =="parametric") & any(degrees != 1))
  {
    cov_type <- "nonparametric"
    warning(paste0("Since some 'degrees' are not equal to one 'cov_type' ",
                   "will be set to 'nonparametric'."))
  }
  
  #---

  # Get the number of equations
  n_eq <- length(formula)
  n_eq2 <- length(formula2)
  
  # Estimate the first step if need
  if (estimator == "2step")
  {
    if (is.null(model1))
    {
      model1 <- mvoprobit(formula = formula, 
                          groups = groups,
                          data = data, 
                          estimator = "ml",
                          marginal = marginal, 
                          opt_type = opt_type,
                          opt_args = opt_args,
                          n_sim = n_sim, 
                          n_cores = n_cores,
                          control = list(cov_type = cov_type))
    }
    out$model1 <- model1
    data <- model1$data
    lambda1 <- model1$lambda
    lambda1_adj <- predict(model1, type = "lambda", 
                           control = list(adj = TRUE))
    if (n_eq > 1)
    {
      lambda1_two <- predict(model1, type = "lambda", 
                             control = list(adj2 = TRUE))
    }
    for (i in 1:model1$control_lnL$n_eq)
    {
      data[[paste0("lambda", i)]] <- lambda1[, i]
      data[[paste0("lambda_adj", i)]] <- lambda1_adj[, i]
      if (n_eq > 1)
      {
        for (j in (1:model1$control_lnL$n_eq)[-i])
        {
          data[[paste0("lambda_two", i, "_", j)]] <- lambda1_two[i, j, ]
        }
      }
    }
    dlambda <- predict(model1, type = "dlambda")
    for (i in 1:model1$control_lnL$n_eq)
    {
      for (j in 1:length(model1$par))
      {
        data[[paste0("dlambda", i, "_", j)]] <- dlambda[, i, j]
      }
    }
  }

  # Deal with control variables
  out_type <- "default"
  if (!is.null(control))
  {
    if (exists("out_type", where = control))
    {
      out_type <- control$out_type
    }
  }
  
  # Coerce formula to type 'formula' and check whether it has "|" symbol
  # which indicates the presence of heteroscedasticity.
  # If there is heteroscedasticity then store separate formulas for
  # mean and variance equations
  is_het <- rep(FALSE, n_eq)
  formula_mean <- vector(mode = "list", length = n_eq)
  formula_var <- vector(mode = "list", length = n_eq)
  for (i in 1:n_eq)
  {
    formula[[i]] <- as.formula(formula[[i]])
    frm_tmp <- paste(formula[[i]][2], formula[[i]][3], sep = '~')
    is_het[i] <- grepl(x = frm_tmp, pattern = "|", fixed = TRUE)
    if (is_het[i])
    {
      frm_tmp <- formula_split(formula[[i]])
      formula_mean[[i]] <- frm_tmp[[1]]
      formula_var[[i]] <- frm_tmp[[2]]
    }
    else
    {
      formula_mean[[i]] <- formula[[i]]
    }
  }
  
  # Get separate dataframe for each equation
  # insuring that all dataframes include the 
  # same observations
    # prepare dataframes
  df <- vector(mode = "list", length = n_eq)
  df_mean <- vector(mode = "list", length = n_eq)
  df_var <- vector(mode = "list", length = n_eq)
  df2 <- NULL
    # indexes of observations to be included
  complete_is <- matrix(NA, 
                        nrow = nrow(data), 
                        ncol = n_eq)
  complete_is2 <- NULL
    # reserve memory for some variables if need
  if (is2)
  {
    df2 <- vector(mode = "list", length = n_eq2)
    complete_is2 <- matrix(NA, 
                           nrow = nrow(data), 
                           ncol = n_eq2)
  }
    # ordered
  frm_tmp <- NULL
  for (i in 1:n_eq)
  {
    if (is_het[i])
    {
      frm_tmp <- switchSelection::formula_merge(
        switchSelection::formula_split(formula[[i]]), 
        type = "var-terms")
    }
    else
    {
      frm_tmp <- formula[[i]]
    }
    df[[i]] <- model.frame(frm_tmp, data, na.action = na.pass)
    z_i_tmp <- as.vector(df[[i]][, 1])
    complete_is[, i] <- complete.cases(df[[i]]) | ((z_i_tmp == -1) & !is.na(z_i_tmp))
  }
    # continuous
  if (is2)
  {
    for (i in 1:n_eq2)
    {
      df2[[i]] <- model.frame(formula2[[i]], data, na.action = na.pass)
      y_i_tmp <- as.vector(df2[[i]][, 1])
      complete_is2[, i] <- complete.cases(df2[[i]]) | is.infinite(y_i_tmp)
    }
  }
    # merge
  complete_ind_vec <- rowSums(complete_is) == n_eq
  complete_ind_vec2 <- rep(TRUE, nrow(data))
  complete_ind_vec3 <- rep(TRUE, nrow(data))
  if (is2)
  {
    complete_ind_vec2 <- rowSums(complete_is2) == n_eq2
  }
  complete_ind <- which(complete_ind_vec & 
                        complete_ind_vec2)
    # preserve only complete observations
  data <- data[complete_ind, ]
    # seperate dataframes for mean and variance equations
      # of ordered variables
  for (i in 1:n_eq)
  {
    df_mean[[i]] <- model.frame(formula_mean[[i]], data, na.action = na.pass)
    if (is_het[i])
    {
      df_var[[i]] <- model.frame(formula_var[[i]], data, na.action = na.pass)
    }
  }
    # dataframe for continuous equations
  if (is2)
  {
    for (i in 1:n_eq2)
    {
      df2[[i]] <- model.frame(formula2[[i]], data, na.action = na.pass)
    }
  }

  # Calculate number of observations
  n_obs <- length(complete_ind)

  # Extract dependent and exogeneous variables
    # ordered variables
  z <- matrix(NA, nrow = n_obs, ncol = n_eq)
  W_mean <- vector("mode" = "list", length = n_eq)
  W_var <- vector("mode" = "list", length = n_eq)
  for (i in 1:n_eq)
  {
    z[, i] <- as.vector(df_mean[[i]][, 1])
    W_mean[[i]] <- as.matrix(df_mean[[i]][, -1])
    if (is_het[i])
    {
      W_var[[i]] <- as.matrix(df_var[[i]][, -1, drop = FALSE])
    }
    else
    {
      W_var[[i]] <- matrix()
    }
  }
    # continuous variables
  y <- matrix(NA)
  X <- list(matrix())
  if (is2)
  {
    y <- matrix(NA, nrow = n_obs, ncol = n_eq2)
    X <- vector("mode" = "list", length = n_eq2)
    for (i in 1:n_eq2)
    {
      y[, i] <- as.vector(df2[[i]][, 1])
      X[[i]] <- cbind(1, as.matrix(df2[[i]][, -1, drop = FALSE]))
      colnames(X[[i]])[1] <- "(Intercept)"
    }
  }

  # Get names of the dependent variables
  z_names <- rep(NA, n_eq)
  y_names <- NULL
  for (i in 1:n_eq)
  {
    z_names[i] <- as.character(formula[[i]][[2]])
  }
  if (n_eq2 > 0)
  {
    y_names <- rep(NA, n_eq2)
    for (i in 1:n_eq2)
    {
      y_names[i] <- all.vars(formula2[[i]])[1]
    }
  }

  # Return data if need
  if (out_type == "data")
  {
    out <- list(n_obs = n_obs,
                dependent = z,
                W_mean = W_mean,
                W_var = W_var,
                X = X,
                complete_ind = complete_ind,
                data = data)

    return(out)
  }
  
  # Calculate the number of cuts for each equation
  n_cuts_eq <- rep(NA, n_eq)
  for (i in 1:n_eq)
  {
    n_cuts_eq[i] <- max(z[, i])
  }
  
  # Set combinations by default if need
  if (is.null(groups))
  {
    groups <- as.matrix(t(hpa::polynomialIndex(n_cuts_eq + 1) - 1))
  }
  n_groups_unique <- nrow(groups)
  
  # Set groups for continuous equations if need
  if (is.null(groups2))
  {
    if (is2)
    {
      groups_tmp <- groups
      groups2_tmp <- as.matrix(t(hpa::polynomialIndex(rep(1, n_eq2)) - 1))
      for (i in 1:nrow(groups2_tmp))
      {
        groups2 <- rbind(groups2,
                         matrix(rep(groups2_tmp[i, ], 
                                    each = n_groups_unique), 
                                nrow = n_groups_unique))
        if (i > 1)
        {
          groups <- rbind(groups, groups_tmp)
        }
      }
    }
    else
    {
      groups2 <- matrix()
    }
  }

  # Get indexes of observations for each group
  ind_g <- NULL
  if (is2)
  {
    y_tmp <- y
    y_tmp[is.infinite(y)] <- -1
    y_tmp[!is.infinite(y)] <- 0
    groups2_tmp <- groups2
    groups2_tmp[groups2_tmp != -1] <- 0
    ind_g <- findGroup(cbind(z, y_tmp), cbind(groups, groups2_tmp))
  }
  else
  {
    ind_g <- findGroup(z, groups)
  }

  # Calculate the number of groups
  n_groups <- length(ind_g)
  
  # Calculate number of observations in each group
  n_obs_g <- rep(NA, n_groups)
  for (i in 1:n_groups)
  {
    n_obs_g[i] <- length(ind_g[[i]])
  }
  n_obs_g_0 <- n_obs_g == 0
  
  # Find groups with zero observable equations
  ind_o_groups <- rep(FALSE, n_groups)
  for (i in 1:n_groups)
  {
    if (any(groups[i, ] != -1))
    {
      ind_o_groups[i] <- TRUE
    }
  }
  
  # Remove groups with zero observations
  # and without observable equations
  ind_g_include <- !n_obs_g_0 & ind_o_groups
  groups <- groups[ind_g_include, , drop = FALSE]
  if (is2)
  {
    groups2 <- groups2[ind_g_include, , drop = FALSE]
  }
  n_obs_g <- n_obs_g[ind_g_include]
  ind_g <- ind_g[ind_g_include]
  n_groups <- nrow(groups[, , drop = FALSE])

  # Find unobservable groups for continuous variables
  if (is2)
  {
    for (i in n_groups)
    {
      for (j in 1:n_eq2)
      {
        if (all(is.infinite(y[ind_g[[i]], j]) | 
                is.na(y[ind_g[[i]], j])))
        {
          groups2[i, j] <- -1
        }
      }
    }
  }

  # Calculate the number of observed equations per group
  n_eq_g <- vector(mode = "numeric", length = n_groups)
  n_eq2_g <- rep(0, n_groups)
  n_eq_all_g <- n_eq_g
  for (i in 1:n_groups)
  {
    n_eq_g[i] <- sum(groups[i, , drop = FALSE] != -1)
  }
  if (is2)
  {
    n_eq2_g <- vector(mode = "numeric", length = n_groups)
    for (i in 1:n_groups)
    {
      n_eq2_g[i] <- sum(groups2[i, , drop = FALSE] != -1)
      n_eq_all_g[i] <- n_eq_g[i] + n_eq2_g[i]
    }
  }

  # Get indexes of observed equations for each group
  ind_eq <- vector(mode = "list", length = n_groups)
  ind_eq2 <- list(NA)
  ind_eq_all <- ind_eq
  for (i in 1:n_groups)
  {
    ind_eq[[i]] <- which(groups[i, , drop = FALSE] != -1)
  }
  if (is2)
  {
    ind_eq2 <- vector(mode = "list", length = n_groups)
    for (i in 1:n_groups)
    {
      ind_eq2[[i]] <- which(groups2[i, , drop = FALSE] != -1)
      ind_eq_all[[i]] <- c(ind_eq[[i]], ind_eq2[[i]] + n_eq)
    }
  }
  else
  {
    ind_eq_all <- ind_eq
  }
  # Get the number of coefficients (regressors)
  # for each equation and determine their indexes
  # in parameters vector
  n_coef <- vector("mode" = "numeric", length = n_eq)
  coef_ind <- vector("mode" = "list", length = n_eq)
  for(i in 1:n_eq)
  {
    n_coef[i] <- ncol(W_mean[[i]])
    if (i != 1)
    {
      coef_ind[[i]] <- (coef_ind[[i - 1]][n_coef[i - 1]] + 1):
                       ((coef_ind[[i - 1]][n_coef[i - 1]] + n_coef[i]))
    }
    else
    {
      coef_ind[[i]] <- 1:n_coef[i]
    }
  }

  # Calculate total number of estimated coefficients
  n_coef_total <- sum(n_coef)
  
  # Get the number of the covariance matrix 
  # elements to be estimated
  n_sigma <- (n_eq ^ 2 - n_eq) / 2
  
  # Get the number of cut points to be considered
  n_cuts <- sum(n_cuts_eq)
  
  # Get indexes of sigma elements in parameters vector
  sigma_ind <- (n_coef_total + 1):(n_coef_total + n_sigma)

  # Store indexes of sigma elements in a matrix form
  sigma_ind_mat <- matrix(NA, nrow = n_eq, ncol = n_eq)
  if (n_eq >= 2)
  {
    counter <- sigma_ind[1]
    for (i in 2:n_eq)
    {
      for (j in 1:(i - 1))
      {
        sigma_ind_mat[i, j] <- counter
        sigma_ind_mat[j, i] <- sigma_ind_mat[i, j]
        counter <- counter + 1
      }
    }
  }

  # Store indexes of cuts in a list form
  cuts_ind <- vector("mode" = "list", length = n_eq)
  n_par <- n_coef_total + n_sigma
  for(i in 1:n_eq)
  {
    for (j in 1:n_cuts_eq[i])
    {
      n_par <- n_par + 1
      cuts_ind[[i]][j] <- n_par
    }
  }
  
  # Deal with coefficients of variance equation
  n_coef_var <- vector("mode" = "numeric", length = n_eq)
  coef_var_ind <- vector("mode" = "list", length = n_eq)
  for(i in 1:n_eq)
  {
    if(is_het[i])
    {
      n_coef_var[i] <- ncol(W_var[[i]])
      coef_var_ind[[i]] <- (n_par + 1):(n_par + n_coef_var[i])
      n_par <- n_par + n_coef_var[i]
    }
    else
    {
      coef_var_ind[[i]] <- numeric()
    }
  }
  
  # Calculate the number of regimes for each continuous variable
  n_regimes <- vector(mode = "numeric", length = 0)
  if (is2)
  {
    n_regimes <- vector(mode = "numeric", length = n_eq2)
    for (i in 1:n_eq2)
    {
      n_regimes[i] <- max(groups2[, i] + 1)
    }
  }
  
  # Convert degrees to matrix
  if (estimator == "2step")
  {
    if (is.matrix(degrees))
    {
      if (nrow(degrees) != n_regimes)
      {
        stop(paste0("If 'degrees' is a matrix then it's number of rows should ",
                    "be equal to the number of regimes that is ", n_regimes[1]))
      }
    }
    else
    {
      degrees <- matrix(rep(degrees, n_regimes[1]), 
                        nrow = n_regimes[1], byrow = TRUE)
    }
  }

  # Deal with coefficients of continuous equations
  n_coef2 <- vector(mode = "numeric", length = 0)
  coef2_ind <- list(matrix())
  if (is2)
  {
    n_coef2 <- vector(mode = "numeric", length = n_eq2)
    coef2_ind <- vector(mode = "list", length = n_eq2)
    for(i in 1:n_eq2)
    {
      n_coef2[i] <- ncol(X[[i]])
      # coefficients for each regime are located in a different rows
      coef2_ind[[i]] <- matrix((n_par + 1):(n_par + n_coef2[i] * n_regimes[i]), 
                               nrow = n_regimes[i], ncol = n_coef2[i], 
                               byrow = TRUE)
      n_par <- n_par + length(coef2_ind[[i]])
    }
  }
  
  # Determine structure of covariances between continuous
  # equations for different regimes
  regimes <- as.matrix(t(hpa::polynomialIndex(n_regimes - 1)))
  n_regimes_total <- nrow(regimes)
  
  # If need add covariances related to continuous equations
  var_y_ind <- list(NA)
  cov_y_ind <- list(matrix(0))
  sigma2_ind <- list(matrix(0))
  sigma2_ind_mat <- list(matrix(0))
  if (is2)
  {
    # variances
    for (i in 1:n_eq2)
    {
      var_y_ind[[i]] <- (n_par + 1):(n_par + n_regimes[i])
      n_par <- n_par + n_regimes[i]
    }
    # covariances with continuous equations
    for (i in 1:n_eq2)
    {
      cov_y_ind[[i]] <- matrix(ncol = n_eq, nrow = n_regimes[i])
      for (j in 1:n_regimes[i])
      {
        cov_y_ind[[i]][j, ] <- (n_par + 1):(n_par + n_eq)
        n_par <- n_par + n_eq
      }
    }
    # covariances between continuous equations
    if (n_eq2 >= 2)
    {
      n_sigma2 <- (n_eq2 ^ 2 - n_eq2) / 2
      sigma2_ind <- vector(mode = "list", length = n_sigma2)
      regimes_pair <- vector(mode = "list", length = n_sigma2)
      sigma2_ind_mat <- vector(mode = "list", length = n_groups)
      for (g in 1:n_groups)
      {
        sigma2_ind_mat[[g]] <- matrix(NA, nrow = n_eq2, ncol = n_eq2)
      }
      # indexes for each covariance depending on regime
      counter <- 1
      for (i in 2:n_eq2)
      {
        for (j in 1:(i - 1))
        {
          sigma2_ind[[counter]] <- numeric()
          regimes_pair[[counter]] <- numeric()
          # note that i > j so in regimes_pair equation with
          # the greater index goes first
          pairs <- unique(groups2[, c(i, j)])
          pairs <- pairs[(pairs[, 1] != -1) & (pairs[, 2] != -1), , 
                         drop = FALSE]
          regimes_pair_n <- ifelse(is.matrix(pairs), nrow(pairs), 0)
          if (regimes_pair_n > 0)
          {
            regimes_pair[[counter]] <- pairs
            sigma2_ind[[counter]] <- (n_par + 1):(n_par + regimes_pair_n)
            n_par <- n_par + regimes_pair_n
            # indexes of covariances depending on group
            for (g in 1:n_groups)
            {
              cond <- which((regimes_pair[[counter]][, 1] == groups2[g, i]) &
                            (regimes_pair[[counter]][, 2] == groups2[g, j]))
              if (length(cond) > 0)
              {
                sigma2_ind_mat[[g]][i, j] <- sigma2_ind[[counter]][cond]
                sigma2_ind_mat[[g]][j, i] <- sigma2_ind[[counter]][cond]
              }
            }
          }
          counter <- counter + 1
        }
      }
    }
  }

  # Parameters associated with marginal distributions
  is_marginal <- length(marginal) > 0
  marginal_par_ind <- list(NA)
  marginal_par_n <- rep(0, n_eq)
  marginal_names <- character()
  if (is_marginal)
  {
    marginal_names <- names(marginal)
    for (i in 1:n_eq)
    {
      marginal_par_n[i] <- ifelse((length(marginal[[i]]) != 0) &
                                  !is.null(marginal[[i]]), 
                                  as.numeric(marginal[[i]]), 0)
      if (marginal_par_n[i] > 0)
      {
        marginal_par_ind[[i]] <- (n_par + 1):(n_par + marginal_par_n[i])
        n_par <- n_par + length(marginal_par_ind[[i]])
      }
    }
  }
  
  # Round the number of parameters to prevent bugs
  marginal_par_n <- round(marginal_par_n)

  # Store the information into the list
  control_lnL <- list(n_par = n_par,
                      n_obs = n_obs,
                      n_eq = n_eq,
                      n_eq2 = n_eq2,
                      n_coef = n_coef,
                      n_coef2 = n_coef2,
                      n_regimes = n_regimes,
                      coef_ind = lapply(coef_ind, function(x){x - 1}),
                      coef_var_ind = lapply(coef_var_ind, function(x){x - 1}),
                      coef2_ind = lapply(coef2_ind, function(x){x - 1}),
                      sigma_ind = sigma_ind - 1,
                      sigma2_ind = lapply(sigma2_ind, function(x){x - 1}),
                      sigma_ind_mat = sigma_ind_mat - 1,
                      sigma2_ind_mat = lapply(sigma2_ind_mat, 
                                              function(x){x - 1}),
                      var_y_ind = lapply(var_y_ind, function(x){x - 1}),
                      cov_y_ind = lapply(cov_y_ind, function(x){x - 1}),
                      cuts_ind = lapply(cuts_ind, function(x){x - 1}),
                      marginal_par_ind = lapply(marginal_par_ind, 
                                                function(x){x - 1}),
                      n_groups = n_groups,
                      ind_g = lapply(ind_g, function(x){x - 1}),
                      ind_eq = lapply(ind_eq, function(x){x - 1}),
                      ind_eq2 = lapply(ind_eq2, function(x){x - 1}),
                      ind_eq_all = lapply(ind_eq_all, function(x){x - 1}),
                      n_cuts_eq = n_cuts_eq,
                      groups = groups,
                      groups2 = groups2,
                      n_obs_g = n_obs_g,
                      n_eq_g = n_eq_g,
                      n_eq2_g = n_eq2_g,
                      n_eq_all_g = n_eq_all_g,
                      is_het = is_het,
                      is2 = is2,
                      y = y,
                      X = X,
                      W = W_mean,
                      W_var = W_var,
                      marginal_names = marginal_names,
                      marginal_par_n = marginal_par_n)
  
  # Second step procedure if need
  par <- NULL
  cov_2step <- NULL
  model2_list <- NULL
  if (estimator == "2step")
  {
    coef_lambda <- vector(mode = "list", length = max(groups2[, 1]) + 1)
    par <- rep(0, n_par)
    cov_2step <- vector(mode = "list", length = n_regimes[1])
    # parameters of the first step
    for (i in 1:n_eq)
    {
      par[cuts_ind[[i]]] <- model1$cuts[[i]]
      par[coef_ind[[i]]] <- model1$coef[[i]]
      par[coef_var_ind[[i]]] <- model1$coef_var[[i]]
      if (is_marginal)
      {
        if (model1$control_lnL$marginal_par_n[i] > 0)
        {
          par[marginal_par_ind[[i]]] <- model1$par[
            model1$control_lnL$marginal_par_ind[[i]] + 1]
        }
      }
    }
    par[sigma_ind] <- model1$par[model1$control_lnL$sigma_ind + 1]
    # parameters of the second step
    model2_list <- vector(mode = "list", length = n_regimes[1])
    # Separate two-step procedure for each regime
    for (i in 1:n_regimes[1])
    {
      ind_regime <- NULL
      for (j in 1:n_groups)
      {
        if (groups2[j] == (i - 1))
        {
          ind_regime <- c(ind_regime, ind_g[[j]])
        }
      }
      # Set formula for particular regime if need
      formula2_tmp <- formula2[[1]]
      if (!is.null(degrees))
      {
        for (j in 1:n_eq)
        {
          if (degrees[i, j, drop = FALSE] != 0)
          {
            for (j1 in 1:degrees[i, j, drop = FALSE])
            {
              formula2_tmp <- update(formula2_tmp, paste0("~. + I(lambda", 
                                                            j, "^", j1, ")"))
            }
          }
        }
      }
      # Estimate least squares regression
      data_regime <- data[ind_regime, ]
      model2_list[[i]] <- lm(formula2_tmp, data = data_regime, 
                             na.action = na.exclude)
      # Save coefficients
      coef2_tmp <- coef(model2_list[[i]])
      ind2_tmp <- coef2_ind[[1]][i, ]
      coef_lambda_tmp <- vector(mode = "list", length = n_eq)
      counter_tmp <- length(ind2_tmp)
      for (t in 1:n_eq)
      {
        if (degrees[i, t] > 0)
        {
          coef_lambda_tmp[[t]] <- coef2_tmp[(counter_tmp + 1): 
                                            (counter_tmp + degrees[i, t])]
        }
        counter_tmp <- counter_tmp + degrees[i, t]
      }
      out$coef_lambda[[i]] <- coef_lambda_tmp
      par[ind2_tmp] <- coef2_tmp[1:length(ind2_tmp)]
      # Calculate unconditional variance
        # calculate squared residuals
      resid <- resid(model2_list[[i]])
      resid2 <- resid ^ 2
        # number of observations
      n_obs_regime <- length(na.omit(resid))
        # get data for these observations
        # get matrices of exogeneous variables
      X_2step <- as.matrix(model.frame(data = data_regime, 
                                       formula = formula2_tmp))
      X_2step[, 1] <- 1
      colnames(X_2step)[1] <- "(Intercept)"
        # estimate common part of the covariance matrix
      X_2step_mat <- qr.solve(t(X_2step) %*% X_2step, tol = 1e-16) %*% t(X_2step)
        # prepare matrix to store estimate
      cov_2step[[i]] <- matrix(1, 
                               nrow = ncol(X_2step),
                               ncol = ncol(X_2step))
      rownames(cov_2step[[i]]) <- names(coef2_tmp)
      colnames(cov_2step[[i]]) <- names(coef2_tmp)
        # derivatives of lambdas
      dlambda <- matrix(0, 
                        nrow = length(ind_regime), 
                        ncol = length(model1$par))
      for (t in 1:length(model1$par))
      {
        for (t1 in 1:n_eq)
        {
          if (degrees[i, t1] != 0)
          {
            for (t2 in 1:degrees[i, t1])
            {
              dlambda[, t] <- dlambda[, t] + 
                              t2 * (data_regime[[paste0("lambda", t1)]] ^ (t2 - 1)) *
                              data_regime[[paste0("dlambda", t1, "_", t)]] *
                              coef_lambda_tmp[[t1]][t2]
            }
          }
        }
      }
        # parametric or nonparametric estimation
      sigma2_cond_triangle <- NULL
      if (cov_type == "parametric")
      {
        par[cov_y_ind[[1]][i, ]] <- coef2_tmp[(length(ind2_tmp) + 1):
                                              (length(ind2_tmp) + n_eq)]
        sigma_adj_p1 <- 0
        sigma_adj_p2 <- 0
        sigma_adj_p3 <- 0
        for (j in 1:n_eq)
        {
          sigma_adj_p1 <- sigma_adj_p1 + 
                          data_regime[[paste0("lambda_adj", j)]] * 
                          (par[cov_y_ind[[1]][i, j]] ^ 2)
          sigma_adj_p2 <- sigma_adj_p2 + 
                          data_regime[[paste0("lambda", j)]] * 
                          par[cov_y_ind[[1]][i, j]]
          for (j1 in (1:n_eq)[-j])
          {
            sigma_adj_p3 <- sigma_adj_p3 + 
                            par[cov_y_ind[[1]][i, j]] * 
                            (par[cov_y_ind[[1]][i, j1]] - 
                             par[sigma_ind_mat[j, j1]] * 
                             par[cov_y_ind[[1]][i, j]]) *
                            data_regime[[paste0("lambda_two", j, "_", j1)]]
          }
        }
        sigma_adj_p2_square <- sigma_adj_p2 ^ 2
        sigma2_uncond <- mean(resid2) - mean(sigma_adj_p1) + 
                         mean(sigma_adj_p2_square) - mean(sigma_adj_p3)
        par[var_y_ind[[1]][i]] <- sigma2_uncond
        sigma2_cond_obs <- sigma2_uncond + sigma_adj_p1 - 
                           sigma_adj_p2_square + sigma_adj_p3
        # Adjust covariance matrix for heteroscedasticity
        sigma2_cond_triangle <- diag(sigma2_cond_obs)
      }
      else 
      {
        # Adjust covariance matrix for heteroscedasticity
        sigma2_cond_triangle <- diag(resid2)
      }
      # Adjust covariance matrix for estimates from the first step
      sigma2_cond_triangle <- sigma2_cond_triangle + (dlambda %*% 
                                                      model1$cov %*% 
                                                      t(dlambda))
      cov_2step[[i]] <- X_2step_mat %*% 
                        sigma2_cond_triangle %*% 
                        t(X_2step_mat)
      }
    }
  # Create a starting point for the numeric optimization
  if (is.null(start) & (estimator == "ml"))
  {
    start <- rep(0, n_par)
    for (i in 1:n_eq)
    {
      # cuts
      for (j in 1:n_cuts_eq[i])
      {
        start[cuts_ind[[i]]][j] <- qnorm(mean(z[, i] <= (j - 1)))
      }
      # coefficients
      z_tmp <- z[, i]
      z_tmp[z_tmp == -1] <- NA
      z_tmp_table <- table(z[, i])
      z_tmp_ind <- which.min(abs((cumsum(z_tmp_table) / 
                                  sum(z_tmp_table)) - 0.5)) - 1
      z_tmp_cond <- z_tmp <= z_tmp_ind
      z_tmp[z_tmp_cond] <- 0
      z_tmp[!z_tmp_cond] <- 1
      data_tmp <- data[complete_ind, ]
      data_tmp[z_names[i]] <- z_tmp
      model_tmp <- glm(formula_mean[[i]], data = data_tmp,
                       family = binomial(link = "logit"))
      start[coef_ind[[i]]] <- coef(model_tmp)[-1] * (sqrt(3) / pi)
    }
    if (n_eq2 > 0)
    {
      for (i in 1:n_eq2)
      {
        data_tmp <- data.frame(y = y[, i], X[[i]])
        for (j in 1:n_regimes[i])
        {
          ind_g_include_tmp <- which(groups2[, i] == (j - 1))
          ing_g_tmp <- NULL
          for (t in ind_g_include_tmp)
          {
            ing_g_tmp <- c(ing_g_tmp, ind_g[[t]])
          }
          data_tmp_g <- data_tmp[ing_g_tmp, ]
          model_tmp <- lm(y ~ . + 0, data = data_tmp_g[!is.na(data_tmp_g$y) & 
                                                       !is.infinite(data_tmp_g$y), ])
          start[coef2_ind[[i]][j, ]] <- coef(model_tmp)
          start[var_y_ind[[i]][j]] <- sigma(model_tmp) ^ 2
        }
      }
    }
    if (is_marginal)
    {
      for (i in 1:n_eq)
      {
        if ((marginal_names[i] == "PGN") | (marginal_names[i] == "hpa"))
        {
          start[marginal_par_ind[[i]]] <- 1e-8
        }
        if (marginal_names[i] == "student")
        {
          start[marginal_par_ind[[i]]] <- 30
        }
      }
    }
  }

  # Code to test derivatives
  # f0 <- lnL_mvoprobit(par = start,
  #                     n_sim = n_sim, n_cores = n_cores,
  #                     control_lnL = control_lnL)
  # grad.a <- lnL_mvoprobit(par = start,
  #                         n_sim = n_sim, n_cores = n_cores,
  #                         control_lnL = control_lnL,
  #                         out_type = "grad")
  # grad.n <- rep(NA, n_par)
  # for (i in 1:n_par)
  # {
  #   delta <- 1e-6
  #   start.delta <- start
  #   start.delta[i] <- start[i] + delta
  #   f1 <- lnL_mvoprobit(par = start.delta,
  #                      n_sim = n_sim, n_cores = n_cores,
  #                      control_lnL = control_lnL)
  #   grad.n[i] <- (f1 - f0) / delta
  # }
  # return(cbind(as.vector(grad.a), as.vector(grad.n)))
  
  # Optimization for one step procedure
  if (estimator == "ml")
  {
    # Perform the optimization routine
    opt <- optim(par = start,
                 method = "BFGS",
                 fn = lnL_mvoprobit,
                 hessian = TRUE,
                 gr = grad_mvoprobit,
                 control = list(maxit = 10000000,
                                fnscale = -1,
                                reltol = 1e-10,
                                abstol = 1e-10),
                 n_sim = n_sim, n_cores = n_cores,
                 control_lnL = control_lnL)
  
  
    # Add genetic optimization if need
    if (opt_type == "gena")
    {
      opt <- gena::gena(fn = lnL_mvoprobit,
                        gr = grad_mvoprobit,
                        pop.initial = opt$par,
                        mutation.method = "percent", 
                        maxiter = 100, info = TRUE,
                        lower = -2 * abs(opt$par), 
                        upper = 2 * abs(opt$par),
                        n_sim = n_sim, n_cores = n_cores,
                        hybrid.prob = 0.2,
                        control_lnL = control_lnL)
    }
    
    # Add genetic optimization if need
    if (opt_type == "pso")
    {
      opt <- gena::pso(fn = lnL_mvoprobit,
                       gr = grad_mvoprobit,
                       pop.initial = opt$par,
                       maxiter = 100, info = TRUE,
                       lower = -2 * abs(opt$par), 
                       upper =  2 * abs(opt$par),
                       n_sim = n_sim, n_cores = n_cores,
                       hybrid.prob = 0,
                       control_lnL = control_lnL)
    }
    
    # Store parameters estimates
    par <- opt$par
  }
  
  # Store parameters
  coef <- vector(mode = "list", length = n_eq)
  coef_var <- vector(mode = "list", length = n_eq)
  cuts <- vector(mode = "list", length = n_eq)
  for (i in 1:n_eq)
  {
    coef[[i]] <- par[coef_ind[[i]]]
    coef_var[[i]] <- par[coef_var_ind[[i]]]
    cuts[[i]] <- par[cuts_ind[[i]]]
  }
  coef2 <- NULL
  var_y <- NULL
  cov_y <- NULL
  if (n_eq2 > 0)
  {
    coef2 <- vector(mode = "list", length = n_eq2)
    var_y <- vector(mode = "list", length = n_eq2)
    cov_y <- vector(mode = "list", length = n_eq2)
    for (i in 1:n_eq2)
    {
      coef2[[i]] <- matrix(nrow = n_regimes[i], ncol = n_coef2[i])
      var_y[[i]] <- par[var_y_ind[[i]]]
      cov_y[[i]] <- matrix(nrow = n_regimes[i], ncol = n_eq)
      for (j in 1:n_regimes[i])
      {
        coef2[[i]][j ,] <- par[coef2_ind[[i]][j ,]]
        cov_y[[i]][j ,] <- par[cov_y_ind[[i]][j ,]]
      }
      rownames(coef2[[i]]) <- paste0("regime ", 0:(n_regimes[i] - 1))
      colnames(coef2[[i]]) <- colnames(X[[i]])
      rownames(cov_y[[i]]) <- paste0("regime ", 0:(n_regimes[i] - 1))
      colnames(cov_y[[i]]) <- 1:n_eq
    }
    names(coef2) <- y_names
  }

  # Store covariance matrix between ordered equations
  sigma <- matrix(NA, nrow = n_eq, ncol = n_eq)
  diag(sigma) <- rep(1, n_eq)
  if (n_eq >= 2)
  {
    sigma_vec <- par[sigma_ind]
    sigma_vec_ind <- matrix(NA, nrow = n_sigma, ncol = 2)
    counter <- 1
    for (i in 2:n_eq)
    {
      for (j in 1:(i - 1))
      {
        sigma_vec_ind[counter, ] <- c(i, j)
        sigma[i, j] <- sigma_vec[counter]
        sigma[j, i] <- sigma[i, j]
        counter <- counter + 1
      }
    }
  }

  # Store covariance matrix between continuous equations
  sigma2 <- NULL
  counter <- 1
  if (n_eq2 >= 2)
  {
    sigma2 <- vector(mode = "list", length = n_sigma2)
    for (i in 2:n_eq2)
    {
      for (j in 1:(i - 1))
      {
        sigma2[[counter]] <- par[sigma2_ind[[counter]]]
        counter <- counter + 1
      }
    }
  }
  
  # Store marginal distribution parameters
  marginal_par <- vector(mode = "list", length = n_eq)
  if (is_marginal)
  {
    for (i in 1:n_eq)
    {
      if (marginal_par_n[i] > 0)
      {
        marginal_par[[i]] <- par[marginal_par_ind[[i]]]
      }
    }
  }
  
  # Estimate asymptotic covariance matrix
  H <- NULL
  H_inv <- NULL
  J <- NULL
  cov <- diag(rep(1, n_par))
  if (!is.matrix(cov_type) & (estimator != "2step"))
  {
    if (cov_type != "no")
    {
      if ((cov_type == "sandwich") | (cov_type == "hessian"))
      {
        H <- gena::gena.hessian(gr = lnL_mvoprobit,
                                par = par,
                                gr.args = list(n_sim = n_sim, n_cores = n_cores,
                                               control_lnL = control_lnL, 
                                               out_type = "grad"))
        out$H <- H
        tryCatch(
        {
          H_inv <- qr.solve(H, tol = 1e-16)
        },
          error = function(e) {
          warning(paste0("Problems with hessian calculation. ",
                         "Standard errors may ",
                         "be estimated inaccurately. Please, try to change ",
                         "'cov_type' argument to 'gop' ",
                         "or change 'opt_type' to 'gena' or 'pso'."))
          }
        )
        if (is.null(H_inv))
        {
          H_inv <- solve(H)
        }
      }
      if ((cov_type == "sandwich") | (cov_type == "gop"))
      {
        J <- lnL_mvoprobit(par = par,
                           n_sim = n_sim, n_cores = n_cores,
                           control_lnL = control_lnL, 
                           out_type = "jac")
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
        cov <- qr.solve(t(J) %*% J, tol = 1e-16)
      }
    }
  } 
  else
  {
    # Manually provided covariance
    if (estimator != "2step")
    {
      cov <- cov_type
    }
    else
    {
      # Covariance of the two-step procedure
      # parameters of the first step
      for (i in 1:n_eq)
      {
        cov[cuts_ind[[i]], 
            cuts_ind[[i]]] <- model1$cov[model1$control_lnL$cuts_ind[[i]] + 1, 
                                         model1$control_lnL$cuts_ind[[i]] + 1]
        cov[coef_ind[[i]], 
            coef_ind[[i]]] <- model1$cov[model1$control_lnL$coef_ind[[i]] + 1, 
                                         model1$control_lnL$coef_ind[[i]] + 1]
        cov[coef_var_ind[[i]], 
            coef_var_ind[[i]]] <- model1$cov[model1$control_lnL$coef_var_ind[[i]] + 1, 
                                             model1$control_lnL$coef_var_ind[[i]] + 1]
        if (is_marginal)
        {
          if (model1$control_lnL$marginal_par_n[i] > 0)
          {
            cov[marginal_par_ind[[i]], 
                marginal_par_ind[[i]]] <- model1$cov[
                  model1$control_lnL$marginal_par_ind[[i]] + 1, 
                  model1$control_lnL$marginal_par_ind[[i]] + 1]
          }
        }
      }
      cov[sigma_ind, sigma_ind] <- model1$cov[model1$control_lnL$sigma_ind + 1, 
                                              model1$control_lnL$sigma_ind + 1]
      # Covariances of the second step
      for (i in 1:n_regimes[1])
      {
        ind2_tmp <- coef2_ind[[1]][i, ]
        cov[ind2_tmp, ind2_tmp] <- cov_2step[[i]][1:length(ind2_tmp),
                                                  1:length(ind2_tmp)]
        if (cov_type == "parametric")
        {
          cov[cov_y_ind[[1]][i, ], 
              cov_y_ind[[1]][i, ]] <- cov_2step[[i]][-(1:length(ind2_tmp)),
                                                     -(1:length(ind2_tmp))]
        }
      }
    }
  }
  se <- sqrt(diag(cov))

  # Calculate p-values
  z_value <- par / se
  p_value <- rep(NA, n_par)
  for (i in 1:n_par)
  {
    p_value[i] <- 2 * min(pnorm(z_value[i]), 1 - pnorm(z_value[i]))
  }
  
  # Construct output table
  tbl_coef <- vector(mode = "list", length = n_eq)
  tbl_coef_var <- vector(mode = "list", length = n_eq)
  tbl_cuts <- vector(mode = "list", length = n_eq)
  tbl_sigma <- vector(mode = "list", length = 1)
  tbl_coef2 <- NULL
  tbl_var_y <- NULL
  tbl_cov_y <- NULL
  tbl_sigma2 <- NULL
  tbl_marginal_par <- NULL
  if (n_eq2 > 0)
  {
    tbl_coef2 <- vector(mode = "list", length = n_eq2)
    tbl_var_y <- vector(mode = "list", length = n_eq2)
    tbl_cov_y <- vector(mode = "list", length = n_eq2)
    tbl_sigma2 <- vector(mode = "list", length = n_eq2 * (n_eq2 - 1) / 2)
  }
  for (i in 1:n_eq)
  {
    # coefficients of mean equation
    tbl_coef[[i]] <- as.matrix(cbind(Estimate = coef[[i]],
                                     Std_Error = se[coef_ind[[i]]],
                                     z_value = z_value[coef_ind[[i]]],
                                     p_value = p_value[coef_ind[[i]]]))
    rownames(tbl_coef[[i]]) <- colnames(W_mean[[i]])
    names(coef[[i]]) <- colnames(W_mean[[i]])
    # coefficients of variance equation
    if (is_het[i])
    {
      tbl_coef_var[[i]] <- as.matrix(cbind(Estimate = coef_var[[i]],
                                        Std_Error = se[coef_var_ind[[i]]],
                                        z_value = z_value[coef_var_ind[[i]]],
                                        p_value = p_value[coef_var_ind[[i]]]))
      rownames(tbl_coef_var[[i]]) <- colnames(W_var[[i]])
      names(coef_var[[i]]) <- colnames(W_var[[i]])
    }
    # cuts
    tbl_cuts[[i]] <- as.matrix(cbind(Estimate = cuts[[i]],
                                     Std_Error = se[cuts_ind[[i]]],
                                     z_value = z_value[cuts_ind[[i]]],
                                     p_value = p_value[cuts_ind[[i]]]))
    rownames(tbl_cuts[[i]]) <- paste0("cut", 1:n_cuts_eq[i])
    names(cuts[[i]]) <- rownames(tbl_cuts[[i]])
  }
  # covariance matrix of ordered equations
  tbl_sigma <- NULL
  if (n_eq >= 2)
  {
    sigma_vec_names <- paste0("sigma", sigma_vec_ind[, 1, drop = FALSE], 
                                       sigma_vec_ind[, 2, drop = FALSE])
    tbl_sigma <- as.matrix(cbind(Estimate = par[sigma_ind],
                                 Std_Error = se[sigma_ind],
                                 z_value = z_value[sigma_ind],
                                 p_value = p_value[sigma_ind]))
    rownames(tbl_sigma) <- sigma_vec_names
  }
  # coef2
  if (n_eq2 > 0)
  {
    for (i in 1:n_eq2)
    {
      tbl_coef2[[i]] <- vector(mode = "list", length = n_regimes[i])
      tbl_cov_y[[i]] <- vector(mode = "list", length = n_regimes[i])
      for (j in 1:n_regimes[i])
      {
        tbl_coef2[[i]][[j]] <- as.matrix(cbind(Estimate = coef2[[i]][j, ],
                                               Std_Error = se[coef2_ind[[i]][j, ]],
                                               z_value = z_value[coef2_ind[[i]][j, ]],
                                               p_value = p_value[coef2_ind[[i]][j, ]]))
        rownames(tbl_coef2[[i]][[j]]) <- colnames(X[[i]])
        names(tbl_coef2[[i]]) <- paste0("regime ", 0:(n_regimes[i] - 1))
        
        tbl_cov_y[[i]][[j]] <- as.matrix(cbind(Estimate = cov_y[[i]][j, ],
                                               Std_Error = se[cov_y_ind[[i]][j, ]],
                                               z_value = z_value[cov_y_ind[[i]][j, ]],
                                               p_value = p_value[cov_y_ind[[i]][j, ]]))
        rownames(tbl_cov_y[[i]][[j]]) <- z_names
        names(tbl_coef2[[i]]) <- paste0("regime ", 0:(n_regimes[i] - 1))
      }
      tbl_var_y[[i]] <- as.matrix(cbind(Estimate = var_y[[i]],
                                        Std_Error = se[var_y_ind[[i]]],
                                        z_value = z_value[var_y_ind[[i]]],
                                        p_value = p_value[var_y_ind[[i]]]))
      rownames(tbl_var_y[[i]]) <- rep("var", n_regimes[i])
    }
    names(tbl_coef2) <- paste0("Equation ", 1:n_eq2)
  }
  # Covariances between continuous equations
  if (n_eq2 >= 2)
  {
    counter <- 1
    for (i in 2:n_eq2)
    {
      for (j in 1:(i - 1))
      {
        tbl_sigma2[[counter]] <- as.matrix(cbind(Estimate = sigma2[[counter]],
                                                 Std_Error = se[sigma2_ind[[counter]]],
                                                 z_value = z_value[sigma2_ind[[counter]]],
                                                 p_value = p_value[sigma2_ind[[counter]]]))
        names(tbl_sigma2)[counter] <- paste0(y_names[i], " and ", y_names[j])
        rownames(tbl_sigma2[[counter]]) <- rep("", nrow(tbl_sigma2[[counter]]))
        for (t in 1:nrow(regimes_pair[[counter]]))
        {
          rownames(tbl_sigma2[[counter]])[t] <- paste0("(",
                                                       regimes_pair[[counter]][t, 1], 
                                                       ",",
                                                       regimes_pair[[counter]][t, 2],
                                                       ")")
        }
        counter <- counter + 1
      }
    }
  }
  # marginal_par
  if (is_marginal)
  {
    tbl_marginal_par <- vector(mode = "list", length = n_eq)
    for (i in 1:n_eq)
    {
      if (marginal_par_n[i] > 0)
      {
        tbl_marginal_par[[i]] <- as.matrix(
          cbind(Estimate = marginal_par[[i]],
          Std_Error = se[marginal_par_ind[[i]] ],
          z_value = z_value[marginal_par_ind[[i]]],
          p_value = p_value[marginal_par_ind[[i]]]))
        rownames(tbl_marginal_par[[i]]) <- paste0("par", 1:marginal_par_n[i])
        names(marginal_par[[i]]) <- paste0("par", 1:marginal_par_n[i])
      }
    }
  }
  # Aggregate the results into the table
  tbl = list(coef = tbl_coef,
             coef_var = tbl_coef_var,
             cuts = tbl_cuts,
             sigma = tbl_sigma,
             coef2 = tbl_coef2,
             var_y = tbl_var_y,
             cov_y = tbl_cov_y,
             sigma2 = tbl_sigma2,
             marginal_par = tbl_marginal_par)
  # standard errors for lambdas
  if (estimator == "2step")
  {
    tbl_lambda <- vector(mode = "list", length = n_regimes[1])
    for (i in 1:n_regimes[1])
    {
      n_coef2_tmp <- length(coef2_ind[[1]][i, ])
      estimate_lambda <- coef(model2_list[[i]])[-(1:n_coef2_tmp)]
      se_lambda <- sqrt(diag(cov_2step[[i]][-(1:n_coef2_tmp), 
                                            -(1:n_coef2_tmp), 
                                            drop = FALSE]))
      z_lambda <- estimate_lambda / se_lambda
      n_lambda <- length(estimate_lambda)
      p_value_lambda <- rep(NA, n_lambda)
      for (j in 1:n_lambda)
      {
        p_value_lambda[j] <- 2 * min(pnorm(z_value[j]), 1 - pnorm(z_lambda[j]))
      }
      tbl_lambda[[i]] <- as.matrix(cbind(Estimate = estimate_lambda,
                                         Std_Error = se_lambda,
                                         z_value = z_lambda,
                                         p_value = p_value_lambda))
      rownames(tbl_lambda[[i]]) <- names(estimate_lambda)
    }
    names(tbl_lambda) <- paste0("regime ", 0:(n_regimes[1] - 1))
    tbl$lambda <- tbl_lambda
  }
  
  # Calculate likelihood value
  loglLik_val <- NA
  if (estimator == "ml")
  {
    loglLik_val <- opt$value
  }
  else
  {
    if (is.null(degrees))
    {
      tryCatch(
        {
          loglLik_val <- lnL_mvoprobit(par = par, 
                                   control_lnL = control_lnL)
          if (loglLik_val <= (-1e+100))
          {
            loglLik_val <- NA
            stop()
          }
        },
        error = function(e) {
          warning(paste0("Second step estimates of covariances and variances ",
                         "can't be combined into estimate of the covariance ",
                         "matrix of random errors. Therefore 'NA' will ",
                         "be returned for Log-likelihood value. ",
                         "It is a common situtation. Estimates ",
                         " and p-values are ",
                         "still interpretable in a usual way."))
        }
      )
    }
  }

  # Store the output
  out$par <- par
  out$coef <- coef
  out$coef_var <- coef_var
  out$coef2 <- coef2
  out$sigma <- sigma
  out$sigma2 <- sigma2
  out$cuts <- cuts
  out$marginal_par <- marginal_par
  out$ind = list(coef = coef_ind,
                 coef_var = coef_var_ind,
                 cuts = cuts_ind,
                 sigma = sigma_ind,
                 g = ind_g,
                 eq = ind_eq,
                 var_y = var_y_ind,
                 sigma2 = sigma2_ind,
                 coef2 = coef2_ind,
                 cov_y = cov_y_ind,
                 sigma2 = sigma2_ind,
                 marginal_par = marginal_par_ind)
  out$logLik <- loglLik_val
  out$W_mean <- W_mean
  out$W_var <- W_var
  out$X <- X
  out$dependent <- z
  out$control_lnL <- control_lnL
  out$formula <- formula
  out$formula2 <- formula2
  out$data_list <- df
  out$data <- data
  out$cov <- cov
  out$cov_2step <- cov_2step
  out$se <- se
  out$p_value <- p_value
  out$tbl <- tbl
  out$groups <- groups
  out$groups2 <- groups2
  out$var_y <- var_y
  out$cov_y <- cov_y
  out$start <- start
  out$marginal <- marginal
  out$twostep <- model2_list
  out$estimator <- estimator
  out$degrees <- as.matrix(degrees)
  out$cov_type <- cov_type
  
  # Estimate lambda for sample selection models
  out$lambda <- predict(out, type = "lambda")

  return(out)
}

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
#' @param object object of class "mnprobit"
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

#' Summary for an Object of Class mvoprobit
#' @description Provides summary for an object of class 'mvoprobit'.
#' @param object object of class "mvoprobit"
#' @param ... further arguments (currently ignored)
#' @details This function just changes the class of the 'mvoprobit'
#' object to 'summary.mvoprobit'.
#' @return Returns an object of class 'summary.mvoprobit'.
summary.mvoprobit <- function(object, ...) 
{
  if (length(list(...)) > 0)
  {
    warning("Additional arguments passed through ... are ignored.")   
  }
  
  class(object) <- "summary.mvoprobit";
  
  return(object)
}

#' Print summary for an Object of Class mvoprobit
#' @description Prints summary for an object of class 'mvoprobit'.
#' @param x object of class "mvoprobit"
#' @param ... further arguments (currently ignored)
#' @return The function returns input argument \code{x} changing
#' it's class to \code{lrtest}.
print.summary.mvoprobit <- function(x, ...)
{
  class(x) <- "mvoprobit"
  tbl_coef <- x$tbl$coef
  tbl_coef_var <- x$tbl$coef_var
  tbl_sigma <- x$tbl$sigma
  tbl_cuts <- x$tbl$cuts
  tbl_coef2 <- x$tbl$coef2
  tbl_var_y <- x$tbl$var_y
  tbl_cov_y <- x$tbl$cov_y
  tbl_sigma2 <- x$tbl$sigma2
  tbl_marginal_par = x$tbl$marginal_par
  tbl_lambda <- NULL
  if (hasName(x$tbl, "lambda"))
  {
    tbl_lambda <- x$tbl$lambda
  }
  n_eq <- x$control_lnL$n_eq
  n_eq2 <- x$control_lnL$n_eq2
  coef_ind <- x$ind$coef
  coef_var_ind <- x$ind$coef_var
  cuts_ind <- x$ind$cuts
  coef2_ind <- x$ind$coef2
  sigma_ind <- x$ind$sigma
  sigma2_ind <- x$ind$sigma2
  var_y_ind <- x$ind$var_y
  cov_y_ind <- x$ind$cov_y
  marginal_par_ind <- x$ind$marginal_par
  n_obs <- x$control_lnL$n_obs
  is_het <- x$control_lnL$is_het
  n_regimes <- x$control_lnL$n_regimes
  
  marginal <- x$marginal
  marginal_par <- x$marginal_par
  marginal_names <- names(marginal)
  is_marginal <- length(marginal) > 0

  cat(paste0(ifelse(n_eq2 == 0, 
                    ifelse(n_eq > 1, "Multivariate ", "Univariate "), 
                    paste0("Selection mechanism is ", 
                           ifelse(n_eq > 1, "multivariate ", "univariate "))),
             ifelse(any(is_het), "heteroscedastic ", ""),
             "ordered ", ifelse(is_marginal, "choice", "probit"), 
             " model\n"))
  if ((n_eq > 1) | (n_eq2 > 0))
  {
    cat(paste0("There are ", n_eq, " ordered ", 
               ifelse(n_eq2 > 0, paste0("and ", n_eq2, " continuous "), ""),
               "equations\n"))
  }
  cat("---\n")
  if (!is.na(x$logLik))
  {
    cat(paste0("Log-likelihood = ", round(x$logLik, 4), "\n"))
  }
  if (!is.na(AIC(x)))
  {
    cat(paste0("AIC = ", round(AIC(x), 4), "\n"))
  }
  cat(paste0("Observations = ", n_obs, "\n"))

  stars <- switchSelection::starsVector(x$p_value)

  for (i in 1:n_eq)
  {
    cat("---\n")
    if ((n_eq >= 2) | (n_eq2 >= 1))
    {
      # cat(paste0("Ordered equation ", i, "\n"))
      cat(paste0(x$formula[[i]][[2]], " equation\n"))
      cat("-\n");
    }
    if (is_marginal)
    {
      cat(paste0("Distribution: ", marginal_names[i], "\n"))
      cat("-\n")
    }
    cat("Coefficients: \n")
    print(as.table(cbind(round(tbl_coef[[i]], 4), 
                         stars[coef_ind[[i]]])))
    if (is_het[i])
    {
      cat("-\n")
      cat("Coefficients of variance equation: \n")
      print(as.table(cbind(round(tbl_coef_var[[i]], 4), 
                           stars[coef_var_ind[[i]]])))
    }
    cat("-\n")
    cat("Cuts: \n")
    print(as.table(cbind(round(tbl_cuts[[i]], 4), 
                         stars[cuts_ind[[i]]])))
    if (length(marginal_par[[i]]) > 0)
    {
      cat("-\n")
      cat(paste0("Parameters of ", marginal_names[i], " distribution: \n"))
      print(as.table(cbind(round(tbl_marginal_par[[i]], 4), 
                           stars[marginal_par_ind[[i]]])))
    }
  }
  
  if (n_eq2 > 0)
  {
    for (i in 1:n_eq2)
    {
      cat("---\n")
      # cat(paste0("Continuous equation ", i, "\n"))
      cat(paste0(all.vars(x$formula2[[i]])[1], " equation\n"))
      tbl_var_y[[i]] <- cbind(round(tbl_var_y[[i]], 4), stars[var_y_ind[[i]]])
      for (j in 1:n_regimes[i])
      {
        cat("--\n")
        if (n_regimes[i] > 1)
        {
          cat(paste0("Regime ", j - 1, ":\n"))
          cat("-\n")
        }
        if (x$estimator == "2step")
        {
          r_squared <- round(summary(x$twostep[[j]])$r.squared, 4)
          r_adj_squared <- round(summary(x$twostep[[j]])$adj.r.squared, 4)
          loocv_val <- round(loocv(x$twostep[[j]]), 4)
          cat(paste0("R-squared: ", r_squared, "\n"))
          cat(paste0("Adjusted R-squared: ", r_adj_squared, "\n"))
          cat(paste0("Leave-one-out cross-validation RMSE: ", loocv_val, "\n"))
          cat("-\n")
        }
        cat("Coefficients:\n")
        print(as.table(cbind(round(tbl_coef2[[i]][[j]], 4), 
                             stars[coef2_ind[[i]][j, ]])))
        cat("-\n")
        if (x$estimator == "2step")
        {
          if (x$cov_type == "parametric")
          {
            cat(paste0("Variance stimate ", tbl_var_y[[i]][j, 1], "\n"))
            cat("standard errors and p-values are unavailable \n")
            cat("for variance when two-step estimator is used\n")
            cat("-\n")
            cat(paste0("Covariances with ordered equations:\n"))
            print(as.table(cbind(round(tbl_cov_y[[i]][[j]], 4), 
                                 stars[cov_y_ind[[i]][j, ]])))
          }
          else
          {
            if (length(tbl_lambda[[j]]) > 1)
            {
              stars_lambda <- switchSelection::starsVector(
                tbl_lambda[[j]][, "p_value"])
              cat(paste0("Selectivity correction terms:\n"))
              print(as.table(cbind(round(tbl_lambda[[j]], 4), 
                                   as.vector(stars_lambda))))
            }
          }
        }
        else
        {
          cat(paste0("Variance:\n"))
          print(as.table(tbl_var_y[[i]][j, , drop = FALSE]))
          cat("-\n")
          cat(paste0("Covariances with ordered equations:\n"))
          print(as.table(cbind(round(tbl_cov_y[[i]][[j]], 4), 
                               stars[cov_y_ind[[i]][j, ]])))
        }
      }
    }
  }
  
  if (n_eq >= 2)
  {
    cat("---\n")
    cat("Correlations between ordered equations: \n")
    print(as.table(cbind(round(tbl_sigma, 4), 
                         stars[sigma_ind])))
  }
  
  if (n_eq2 >= 2)
  {
    cat("---\n")
    cat("Covariances between continuous equations: \n")
    counter <- 1
    for (i in 2:n_eq2)
    {
      for (j in 1:(i - 1))
      {
        cat("-\n")
        cat(paste0("Between ", x$formula2[[i]][[2]], 
                   " and ", x$formula2[[j]][[2]], "\n"))
        print(as.table(cbind(round(tbl_sigma2[[counter]], 4), 
                             stars[sigma2_ind[[counter]]])))
        counter <- counter + 1
      }
    }
  }
  
  cat("---\n")
  cat("Signif. codes:  0 '***' 0.01 '**' 0.05 '*' 0.1 \n")
  
  return(x)
}

#' Print for an Object of Class mvoprobit
#' @description Prints summary for an object of class 'mvoprobit'.
#' @param x object of class "mvoprobit"
#' @param ... further arguments (currently ignored)
#' @return The function returns input argument \code{x}.
print.mvoprobit <- function(x, ...)
{
  print.summary.mvoprobit(x)
  return(x)
}

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
                              control = list())
{
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
    for (i in 1:n_eq)
    {
      if (length(coef_var_ind[[i]]) != 0)
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
      lambda <- matrix(NA, nrow = n, ncol = object$control_lnL$n_eq)
    }
    else
    {
      lambda <- array(dim = c(object$control_lnL$n_eq, 
                              object$control_lnL$n_eq, n))
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
                           matrix(object$cov_y[[i]][group2[i] + 1, ind_eq], 
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
      return(y_pred)
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
  if (is.null(group))
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
  if (is_marginal)
  {
    marginal <- marginal[ind_eq]
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
                          marginal = marginal)$prob
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
                           marginal = marginal,
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