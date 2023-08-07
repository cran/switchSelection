#' Multivariate ordered probit model with heteroscedasticity 
#' and (non-random) sample selection.
#' @description This function allows to estimate parameters of multivariate
#' ordered probit model and its extensions. 
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
#' @template mvoprobit_examples_cps_Template
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
                      control = NULL,
                      regularization = NULL)
{
  # Get first step model if it has been provided
  model1 <- NULL
  if (is(object = formula, class2 = "mvoprobit"))
  {
    model1 <- formula
    formula <- model1$formula
  }
  
  # List to store the output
  out <- list()
  class(out) <- "mvoprobit"
  out$other <- list(n_cores = n_cores, n_sim = n_sim)
  
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
    cov_type_vec <- c("parametric", "nonparametric", "no")
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
  out$other$is2 <- is2
  
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
      for (j in seq_len(model1$control_lnL$n_eq)[-i])
      {
        data[[paste0("lambda_two", i, "_", j)]] <- lambda1_two[i, j, ]
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
    complete_is[, i] <- complete.cases(df[[i]]) | ((z_i_tmp == -1) & 
                                       !is.na(z_i_tmp))
  }
    # continuous
  if (is2)
  {
    for (i in 1:n_eq2)
    {
      df2[[i]] <- model.frame(formula2[[i]], data, 
                              na.action = na.pass)
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
    df_mean[[i]] <- model.frame(formula_mean[[i]], data, 
                                na.action = na.pass)
    if (is_het[i])
    {
      df_var[[i]] <- model.frame(formula_var[[i]], data, 
                                 na.action = na.pass)
    }
  }
    # dataframe for continuous equations
  if (is2)
  {
    for (i in 1:n_eq2)
    {
      df2[[i]] <- model.frame(formula2[[i]], data, 
                              na.action = na.pass)
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
  else
  {
    if (is2)
    {
      if (is.vector(groups2))
      {
        groups2 <- matrix(groups2, ncol = 1)
      }
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
  for (i in 1:n_groups)
  {
    n_eq_g[i] <- sum(groups[i, , drop = FALSE] != -1)
  }
  n_eq_all_g <- n_eq_g
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
  out$other$n_regimes <- n_regimes
  
  # Convert degrees to matrix
  if (estimator == "2step")
  {
    if (is.matrix(degrees))
    {
      if (nrow(degrees) != n_regimes)
      {
        stop(paste0("If 'degrees' is a matrix then its number of rows should ",
                    "be equal to the number of regimes that is ", n_regimes[1]))
      }
    }
    else
    {
      degrees <- matrix(rep(degrees, n_regimes[1]), 
                        nrow = n_regimes[1], byrow = TRUE)
    }
  }
  
  # Check that there is no non-zero degrees corresponding to constant lambda
  if ((estimator == "2step") & !is.null(degrees))
  {
    for (i in 1:n_regimes[1])
    {
      for (j in 1:n_eq)
      {
        if (all(groups[groups2[, 1] == (i - 1), j] == -1))
        {
          degrees[i, j] <- 0
          cov_type <- "nonparametric"
          warning(paste0("Since in regime ", i - 1, " value of the ", j,
                         " selection equation is always unobservable ",
                         " it is assumed",
                         " that 'degrees[", i - 1, ", ", j, "] = 0' and ",
                         "'cov_type = 'nonparametric''."))
        }
      }
    }
  }

  # Deal with coefficients of continuous equations
  n_coef2 <- vector(mode = "numeric", length = 0)
  coef2_ind <- list(matrix(1))
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
  var2_ind <- list(NA)
  cov2_ind <- list(matrix(0))
  sigma2_ind <- list(matrix(0))
  sigma2_ind_mat <- list(matrix(0))
  if (is2)
  {
    # variances
    for (i in 1:n_eq2)
    {
      var2_ind[[i]] <- (n_par + 1):(n_par + n_regimes[i])
      n_par <- n_par + n_regimes[i]
    }
    # covariances with continuous equations
    for (i in 1:n_eq2)
    {
      cov2_ind[[i]] <- matrix(ncol = n_eq, nrow = n_regimes[i])
      for (j in 1:n_regimes[i])
      {
        cov2_ind[[i]][j, ] <- (n_par + 1):(n_par + n_eq)
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
  out$other$is_marginal <- is_marginal
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
                      var2_ind = lapply(var2_ind, function(x){x - 1}),
                      cov2_ind = lapply(cov2_ind, function(x){x - 1}),
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
  out$control_lnL <- control_lnL
  
  # Save indexes
  out$ind = list(coef = coef_ind,
                 coef_var = coef_var_ind,
                 cuts = cuts_ind,
                 sigma = sigma_ind,
                 g = ind_g,
                 eq = ind_eq,
                 var2 = var2_ind,
                 sigma2 = sigma2_ind,
                 coef2 = coef2_ind,
                 cov2 = cov2_ind,
                 sigma2 = sigma2_ind,
                 marginal_par = marginal_par_ind)
  
  # Validate regularization
  regularization <- regularization_validate(regularization = regularization, 
                                            n_par = n_par, 
                                            estimator = estimator)
  
  # Second step procedure if need
  par <- NULL
  cov_2step <- NULL
  model2_list <- NULL
  coef_lambda <- NULL
  coef_lambda_row <- NULL
  if (estimator == "2step")
  {
    coef_lambda <- vector(mode = "list", length = max(groups2[, 1]) + 1)
    coef_lambda_row <- vector(mode = "list", length = max(groups2[, 1]) + 1)
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
      # indexes of observations corresponding
      # to the regime
      ind_regime <- NULL
      for (j in 1:n_groups)
      {
        if (groups2[j] == (i - 1))
        {
          ind_regime <- c(ind_regime, ind_g[[j]])
        }
      }
      # regime specific data
      data_regime <- data[ind_regime, ]
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
      model2_list[[i]] <- lm(formula2_tmp, data = data_regime, 
                             na.action = na.exclude)
      # Save coefficients
      coef2_tmp <- coef(model2_list[[i]])
      n_coef2_tmp <- length(coef2_ind[[1]][i, ])
      coef_lambda_row[[i]] <- coef(model2_list[[i]])[-(1:n_coef2_tmp)]
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
        names(coef_lambda_tmp)[t] <- z_names[t]
        counter_tmp <- counter_tmp + degrees[i, t]
      }
      coef_lambda[[i]] <- coef_lambda_tmp
      names(coef_lambda)[[i]] <- paste("regime ", i - 1)
      par[ind2_tmp] <- coef2_tmp[1:length(ind2_tmp)]
      # Calculate unconditional variance
        # calculate squared residuals
      resid <- resid(model2_list[[i]])
      resid2 <- resid ^ 2
        # number of observations
      n_obs_regime <- length(na.omit(resid))
        # get data for these observations
        # get matrices of exogenous variables
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
        par[cov2_ind[[1]][i, ]] <- coef2_tmp[(length(ind2_tmp) + 1):
                                              (length(ind2_tmp) + n_eq)]
        sigma_adj_p1 <- 0
        sigma_adj_p2 <- 0
        sigma_adj_p3 <- 0
        for (j in 1:n_eq)
        {
          sigma_adj_p1 <- sigma_adj_p1 + 
                          data_regime[[paste0("lambda_adj", j)]] * 
                          (par[cov2_ind[[1]][i, j]] ^ 2)
          sigma_adj_p2 <- sigma_adj_p2 + 
                          data_regime[[paste0("lambda", j)]] * 
                          par[cov2_ind[[1]][i, j]]
          for (j1 in (1:n_eq)[-j])
          {
            sigma_adj_p3 <- sigma_adj_p3 + 
                            par[cov2_ind[[1]][i, j]] * 
                            (par[cov2_ind[[1]][i, j1]] - 
                             par[sigma_ind_mat[j, j1]] * 
                             par[cov2_ind[[1]][i, j]]) *
                            data_regime[[paste0("lambda_two", j, "_", j1)]]
          }
        }
        sigma_adj_p2_square <- sigma_adj_p2 ^ 2
        sigma2_uncond <- mean(resid2) - mean(sigma_adj_p1) + 
                         mean(sigma_adj_p2_square) - mean(sigma_adj_p3)
        par[var2_ind[[1]][i]] <- sigma2_uncond
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
      #return(list(dlambda = dlambda, X_2step_mat = X_2step_mat, 
      #            sigma2_cond_triangle = sigma2_cond_triangle))
      cov_2step[[i]] <- X_2step_mat %*% 
                        sigma2_cond_triangle %*% 
                        t(X_2step_mat)
      }
  }
  out$cov_2step <- cov_2step
  out$coef_lambda <- coef_lambda
  
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
    for (i in seq_len(n_eq2))
    {
      data_tmp <- data.frame(y = y[, i], X[[i]])
      # control for infinite endogenous regressors
      data_tmp <- data_tmp[rowSums(is.infinite(X[[i]])) == 0, ]
      for (j in 1:n_regimes[i])
      {
        ind_g_include_tmp <- which(groups2[, i] == (j - 1))
        ing_g_tmp <- NULL
        for (t in ind_g_include_tmp)
        {
          ing_g_tmp <- c(ing_g_tmp, ind_g[[t]])
        }
        data_tmp_g <- data_tmp[ing_g_tmp, ]
        model_tmp <- lm(y ~ . + 0, 
                        data = data_tmp_g[!is.na(data_tmp_g$y) & 
                                          !is.infinite(data_tmp_g$y), ])
        start[coef2_ind[[i]][j, ]] <- coef(model_tmp)
        start[var2_ind[[i]][j]] <- sigma(model_tmp) ^ 2
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
  opt <- NULL
  if (estimator == "ml")
  {
    opt <- opt_switchSelection(opt_args = opt_args, control_lnL = control_lnL,
                               n_sim = n_sim, n_cores = n_cores,
                               type = "mvoprobit", start = start, 
                               opt_type = opt_type, 
                               regularization = regularization)
    par <- opt$par
  }
  
  # Store parameters into variables
  par_list <- par_mvoprobit(par = par, 
                            is_het = is_het, is_marginal = is_marginal,
                            n_coef2 = n_coef2,
                            n_sigma = n_sigma, n_sigma2 = n_sigma2,
                            coef_ind = coef_ind, coef_var_ind = coef_var_ind, 
                            cuts_ind = cuts_ind, coef2_ind = coef2_ind, 
                            sigma_ind = sigma_ind, var2_ind = var2_ind, 
                            cov2_ind = cov2_ind, sigma2_ind = sigma2_ind, 
                            marginal_par_ind = marginal_par_ind,
                            marginal_par_n = marginal_par_n,
                            y_names = y_names, z_names = z_names, 
                            n_eq = n_eq, n_eq2 = n_eq2, 
                            n_regimes = n_regimes, regimes_pair = regimes_pair,
                            colnames_X = lapply(X, colnames), 
                            colnames_W_mean = lapply(W_mean, colnames), 
                            colnames_W_var = lapply(W_var, colnames))
  par <- par_list$par
  coef <- par_list$coef
  coef_var <- par_list$coef_var
  cuts <- par_list$cuts
  coef2 <- par_list$coef2
  sigma <- par_list$sigma
  var2 <- par_list$var2
  cov2 <- par_list$cov2
  sigma2 <- par_list$sigma2
  marginal_par <- par_list$marginal_par
  sigma_vec_ind <- par_list$sigma_vec_ind
  
  # Store variables in the output list
  out$par <- par
  out$coef <- coef
  out$coef_var <- coef_var
  out$cuts <- cuts
  out$coef2 <- coef2
  out$sigma <- sigma
  out$var2 <- var2
  out$cov2 <- cov2
  out$sigma2 <- sigma2
  out$marginal_par <- marginal_par
  out$other$sigma_vec_ind <- sigma_vec_ind
  
  # Estimate asymptotic covariance matrix
    # Prepare some values
  H <- NULL
  H_inv <- NULL
  J <- NULL
  cov <- diag(rep(1, n_par))
    # Manual covariance matrix
  if (is.matrix(cov_type))
  {
    cov <- cov_type
  }
    # Maximum-likelihood covariance matrix
  if (!is.matrix(cov_type) & (estimator == "ml"))
  {
    vcov_object <- vcov_ml(object = out, type = cov_type, 
                           n_cores = n_cores, n_sim = n_sim)
    cov <- vcov_object$vcov
    if (hasName(vcov_object, "H"))
    {
      H <- vcov_object$H
      out$H <- H
    }
    if (hasName(vcov_object, "J"))
    {
      J <- vcov_object$J
      out$J <- J
    }
  }
    # Two-step covariance matrix
  if (!is.matrix(cov_type) & (estimator == "2step"))
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
        cov[cov2_ind[[1]][i, ], 
            cov2_ind[[1]][i, ]] <- cov_2step[[i]][-(1:length(ind2_tmp)),
                                                   -(1:length(ind2_tmp))]
      }
    }
  }
  
  # Create tables to store the results in a format used by summary function
  tbl_list <- tbl_mvoprobit(par = par, cov = cov, 
                            n_par = n_par, n_eq = n_eq, n_eq2 = n_eq2,
                            n_cuts_eq = n_cuts_eq,
                            n_regimes = n_regimes, regimes_pair = regimes_pair,
                            coef = coef, coef_ind = coef_ind,
                            coef_var = coef_var, coef_var_ind = coef_var_ind,
                            cuts = cuts, cuts_ind = cuts_ind,
                            sigma = sigma, sigma_ind = sigma_ind, 
                            sigma_vec_ind = sigma_vec_ind,
                            coef2 = coef2, coef2_ind = coef2_ind,
                            var2 = var2, var2_ind = var2_ind,
                            cov2 = cov2, cov2_ind = cov2_ind,
                            sigma2 = sigma2, sigma2_ind = sigma2_ind,
                            marginal_par = marginal_par, 
                            marginal_par_ind = marginal_par_ind, 
                            marginal_par_n = marginal_par_n,
                            is_marginal = is_marginal, is_het = is_het, 
                            z_names = z_names, y_names = y_names,
                            estimator = estimator, 
                            coef_lambda_row = coef_lambda_row,
                            cov_2step = cov_2step)
  out$tbl <- tbl_list$tbl
  out$se <- tbl_list$se
  out$p_value <- tbl_list$p_value
  
  # Calculate likelihood value
  logLik_val <- NA
  if (estimator == "ml")
  {
    if (is.null(regularization))
    {
      logLik_val <- opt$value
    }
    else
    {
      # because of regularization
      logLik_val <- lnL_mvoprobit(par = par, control_lnL = control_lnL)
    }
  }
  else
  {
    if (is.null(degrees))
    {
      tryCatch(
        {
          logLik_val <- lnL_mvoprobit(par = par, control_lnL = control_lnL)
          if (logLik_val <= (-1e+100))
          {
            logLik_val <- NA
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
  out$logLik <- logLik_val
  out$W_mean <- W_mean
  out$W_var <- W_var
  out$X <- X
  out$dependent <- z
  out$formula <- formula
  out$formula2 <- formula2
  out$data_list <- df
  out$data <- data
  out$cov <- cov
  out$cov_2step <- cov_2step
  colnames(groups) <- z_names
  out$groups <- groups
  colnames(groups2) <- y_names
  out$groups2 <- groups2
  out$start <- start
  out$marginal <- marginal
  out$twostep <- model2_list
  out$estimator <- estimator
  out$degrees <- as.matrix(degrees)
  out$cov_type <- cov_type
  out$other$z_names <- z_names
  out$other$y_names <- y_names
  out$other$coef_lambda_row <- coef_lambda_row
  if (n_eq2 >= 2)
  {
    out$other$regimes_pair <- regimes_pair
  }
  out$other$is_het <- is_het
  
  # Estimate lambda for sample selection models
  out$lambda <- predict(out, type = "lambda")

  return(out)
}