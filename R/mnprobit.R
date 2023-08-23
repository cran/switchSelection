#' Multinomial probit model
#' @description This function estimates parameters of multinomial probit
#' model and sample selection model with continuous outcome and multinomial
#' probit selection mechanism.
#' @template mnprobit_details_Template
#' @template mnprobit_param_Template
#' @template mnprobit_examples_cps_Template
#' @template mnprobit_examples1_Template
#' @template mnprobit_examples2_Template
#' @template mvoprobit_references_Template
#' @template mnprobit_return_Template
mnprobit <- function(formula, 
                     formula2 = NULL,
                     data,
                     regimes = NULL,
                     opt_type = "optim",
                     opt_args = NULL,
                     start = NULL,
                     estimator = "ml",
                     cov_type = "sandwich",
                     degrees = NULL,
                     n_sim = 1000, 
                     n_cores = 1,
                     control = NULL,
                     regularization = NULL)
{   
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
  
  #---
  
  
  # Output list
  out <- list()
  class(out) <- "mnprobit"
  out$other <- list(n_cores = n_cores, n_sim = n_sim)
  
  # Deal with control variables
  out_type <- "default"
  if (!is.null(control))
  {
    if (exists("out_type", where = control))
    {
      out_type <- control$out_type
    }
  }
  
  # Get first step model if it has been provided
  model1 <- NULL
  if (is(object = formula, class2 = "mnprobit"))
  {
    model1 <- formula
    formula <- model1$formula
  }
  
  # Estimate the first step if need
  if (estimator == "2step")
  {
    if (is.null(model1))
    {
      model1 <- mnprobit(formula = formula, 
                         data = data, 
                         estimator = "ml",
                         opt_type = opt_type,
                         opt_args = opt_args,
                         n_sim = n_sim, 
                         cov_type = cov_type,
                         n_cores = n_cores)
    }
    n_alt <- model1$n_alt
    out$model1 <- model1
    data <- model1$data
    lambda1 <- model1$lambda
    for (i in 1:(model1$control_lnL$n_alt - 1))
    {
      data[[paste0("lambda", i)]] <- lambda1[, i]
    }
    dlambda <- predict(model1, type = "dlambda", alt = NULL)
    for (i in 1:(model1$control_lnL$n_alt - 1))
    {
      for (j in 1:length(model1$par))
      {
        data[[paste0("dlambda", i, "_", j)]] <- dlambda[, i, j]
      }
    }
  }

  # Coerce formula to type 'formula'
  # and check whether there are any
  # continuous equations
  tryCatch(formula <- as.formula(formula),
           error = function(e) {stop("Invalid 'formula' argument.")})
  is2 <- FALSE
  if (!is.null(formula2))
  {
    tryCatch(formula2 <- as.formula(formula2),
             error = function(e) {stop("Invalid 'formula2' argument.")})
    is2 <- TRUE
  }
  
  # Get dataframes according to the formula
  df <- model.frame(formula, data, na.action = na.pass)
  df2 <- NULL 
  if (is2)
  {
    df2 <- model.frame(formula2, data, na.action = na.pass)
  }
  
  # Deal with regimes
  if (!is2)
  {
    regimes <- rep(-1, max(df[, 1], na.rm = TRUE))
  }
  else
  {
    if ((length(regimes) == 1) | is.null(regimes))
    {
      regimes <- rep(0, max(df[, 1], na.rm = TRUE))
    }
  }
  n_regimes <- ifelse(is2, max(regimes) + 1, 0)
  regimes_names <- NULL
  if (n_regimes > 0)
  {
    regimes_names <-  paste0("regime ", 0:(n_regimes - 1))
  }
  
  # Get indexes of complete observations
  z_tmp <- as.vector(df[, 1])
  complete_is <- complete.cases(df) | ((z_tmp == -1) & !is.na(z_tmp))
  y_tmp <- NULL
  complete_is2 <- rep(TRUE, length(complete_is))
  if (is2)
  {
    y_tmp <- as.vector(df[, 2])
    complete_is2 <- complete.cases(df2) | is.infinite(y_tmp)
  }
  complete_ind <- which(complete_is & complete_is2)
  
  # Remove incomplete observations from data
  data <- data[complete_ind, ]
  df <- model.frame(formula, data, na.action = na.pass)
  if (is2)
  {
    df2 <- model.frame(formula2, data, na.action = na.pass)
  }
  
  # Extract dependent variables and regressors
    # multinomial equation
  z <- df[, 1]                                             # dependent variable
  W <- as.matrix(cbind(1, df[, -1, drop = FALSE]))         # regressors
  colnames(W)[1] <- "(Intercept)"
    # continuous equation
  y <- numeric()
  X <- matrix()
  if (is2)
  {
    y <- df2[, 1]
    X <- as.matrix(cbind(1, df2[, -1, drop = FALSE]))
  }
  colnames(X)[1] <- "(Intercept)"
  
  # Get the number of coefficients (regressors)
  n_coef <- ncol(W)
  n_coef2 <- 0
  if (is2)
  {
    n_coef2 <- ncol(X)
  }
  
  # Get all the alternatives
  alt <- sort(unique(z))
  
  # Calculate the number of alternatives
  n_alt <- length(alt)
  alt_names <- paste0("Alternative ", 1:n_alt)
  if (length(regimes) != n_alt)
  {
    stop(paste0("Incorrect 'regimes' argument. It should be a vector of the ",
                "same length as the number of alternatives."))
  }
  
  # Validate degrees
  if (is.null(degrees))
  {
    degrees <- rep(1, n_alt - 1)
  }
  if (any((degrees %% 1) != 0) | any(degrees < 0))
  {
    stop(paste0("Invalid 'degrees' argument. Please, insure that ", 
                "it is a vector of non-negative integers."))
  }
  
  # Estimate the number of observations
  n_obs <- length(z)
  
  # Return data if need
  if (out_type == "data")
  {
    out <- list(n_obs = n_obs,
                W = W,
                z = z,
                X = X,
                y = y,
                complete_ind = complete_ind,
                data = data)
    
    return(out)
  }
  
  # Get indexes of observations associated with  
  # different alternatives and regimes
  ind_alt <- vector(mode = "list", length = n_alt)
  ind_regime <- vector(mode = "list", length = n_regimes)
  for (i in 1:n_alt)
  {
    ind_alt[[i]] <- which(z == i)
    if (is2)
    {
      if (is.infinite(y[ind_alt[[i]][1]]))
      {
        regimes[i] <- -1
      }
      if (regimes[i] != -1)
      {
        ind_regime[[regimes[i] + 1]] <- c(ind_regime[[regimes[i] + 1]], 
                                          ind_alt[[i]])
      }
    }
  }

  # Get the number of the covariance matrix
  # elements to be estimated
  n_sigma <- ((n_alt - 1) ^ 2 - (n_alt - 1)) / 2 + (n_alt - 2)
  
  # Get the number of regression coefficients
  n_coef_total <- n_coef * (n_alt - 1)
  n_coef2_total <- n_coef2 * n_regimes
  
  # The number of covariance matrix parameters associated with
  # the presence of continuous equation
  n_cov2 <- ifelse(is2, n_alt * n_regimes, 0)
  
  # Calculate the total number of estimated parameters
  # without continuous equations
  n_par <- n_coef_total + n_sigma
  
  # Get indexes of coefficients
  coef_ind <- 1:n_coef_total
  
  # Store the coefficients indexes for each alternative
  coef_ind_alt <- matrix(1, nrow = n_coef, ncol = n_alt - 1)
  for (i in 1:(n_alt - 1)) 
  {
    coef_ind_alt[, i] <- ((i - 1) * n_coef + 1):(i * n_coef)
  }
  
  # Get indexes of the covariances between alternatives
  sigma_ind <- 1
  sigma_ind_mat <- matrix(1)
  if (n_sigma > 0)
  {
    sigma_ind <- (n_coef_total + 1):n_par
    sigma_ind_mat <- matrix(1, nrow = n_sigma, ncol = 2)
  }
  if (n_alt > 2)
  {
    counter <- 1
    for (i in 1:(n_alt - 1))
    {
      for (j in 1:i)
      {
        if (!((i == 1) & (j == 1)))
        {
          sigma_ind_mat[counter, 1] = i;
          sigma_ind_mat[counter, 2] = j;
          counter <- counter + 1
        }
      }
    }
  }
  
  # Parameters associated with continuous equations
  coef2_ind_regime <- matrix(1)
  cov2_ind_regime <- matrix(1)
  var2_ind_regime <- 1
  if (is2)
  {
    coef2_ind_regime <- matrix(1, nrow = n_coef2, ncol = n_regimes)
    cov2_ind_regime <- matrix(1, nrow = n_alt - 1, ncol = n_regimes)
    var2_ind_regime <- vector(mode = "numeric", length = n_regimes)
    for (i in 1:n_regimes) 
    {
      coef2_ind_regime[, i] <- (n_par + 1):(n_par + n_coef2)
      n_par <- n_par + n_coef2
    }
    for (i in 1:n_regimes) 
    {
      cov2_ind_regime[, i] <- (n_par + 1):(n_par + n_alt - 1)
      n_par <- n_par + n_alt - 1
    }
    var2_ind_regime <- (n_par + 1):(n_par + n_regimes)
    n_par <- n_par + n_regimes
  }
  
  # Store the information into the list
  control_lnL <- list(n_par = n_par,
                      n_alt = n_alt,
                      n_regimes = n_regimes,
                      n_obs = n_obs,
                      n_coef = n_coef,
                      n_coef2 = n_coef2,
                      n_sigma = n_sigma,
                      coef_ind = coef_ind - 1,
                      sigma_ind = sigma_ind - 1,
                      sigma_ind_mat = sigma_ind_mat - 1,
                      coef_ind_alt = coef_ind_alt - 1,
                      coef2_ind_regime = coef2_ind_regime - 1,
                      cov2_ind_regime = cov2_ind_regime - 1,
                      var2_ind_regime = var2_ind_regime - 1,
                      ind_alt = lapply(ind_alt, function(x){x - 1}),
                      ind_regime = lapply(ind_regime, function(x){x - 1}),
                      regimes = regimes,
                      z = as.vector(z),
                      W = as.matrix(W),
                      y = y,
                      X = X)
  out$control_lnL <- control_lnL
  
  # Convert degrees to matrix
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
    degrees <- matrix(rep(degrees, n_regimes), 
                      nrow = n_regimes, byrow = TRUE)
  }
  
  # Create a starting point for the numeric 
  # optimization routine
  if (is.null(start) & (estimator != "2step"))
  {
    start <- rep(0, n_par)
    cumsum_i <- 0
    for (i in 2:(n_alt - 1))                               # starting values
    {                                                      # for variances
      cumsum_i <- cumsum_i + i                             # should be equal
      start[n_coef_total + cumsum_i] <- 1                  # to one
    }
    if (is2)
    {
      for (i in 1:n_regimes)
      {
        model_tmp <- lm(formula2, data = data[ind_regime[[i]], ])
        start[coef2_ind_regime[, i]] <- coef(model_tmp)
        start[var2_ind_regime[i]] <- sigma(model_tmp) ^ 2
      }
    }
  }

  # Code to test derivatives
  if (n_sim == -1)
  {
    f0 <- lnL_mnprobit(par = start,
                       n_sim = n_sim, n_cores = n_cores,
                       control_lnL = control_lnL)
    grad.a <- lnL_mnprobit(par = start,
                           n_sim = n_sim, n_cores = n_cores,
                           control_lnL = control_lnL,
                           out_type = "grad")
    grad.n <- rep(NA, n_par)
    for (i in 1:n_par)
    {
      delta <- 1e-6
      start.delta <- start
      start.delta[i] <- start[i] + delta
      f1 <- lnL_mnprobit(par = start.delta,
                         n_sim = n_sim, n_cores = n_cores,
                         control_lnL = control_lnL)
      grad.n[i] <- (f1 - f0) / delta
    }
    return(cbind(analytical = as.vector(grad.a),
                 numeric = as.vector(grad.n)))
  }
  
  # Validate regularization
  regularization <- regularization_validate(regularization = regularization, 
                                            n_par = n_par, 
                                            estimator = estimator)

  # Two-step estimation procedure
  par <- NULL
  cov_2step <- NULL
  model2_list <- NULL
  coef_lambda <- vector(mode = "list", length = n_regimes)
  coef_lambda_row <- vector(mode = "list", length = n_regimes)
  names(coef_lambda) <- regimes_names
  if (estimator == "2step")
  {
    par <- rep(0, n_par)
    cov_2step <- vector(mode = "list", length = n_regimes)
    # parameters of the first step
    for (i in 1:(n_alt - 1))
    {
      par[coef_ind_alt[, i]] <- model1$coef[, i]
    }
    par[sigma_ind] <- model1$par[sigma_ind]
    # parameters of the second step
    model2_list <- vector(mode = "list", length = n_regimes)
    # separate two-step procedure for each regime
    formula2_tmp <- formula2
    for (i in 1:n_regimes)
    {
      # set formula for particular regime if need
      if (!is.null(degrees))
      {
        formula2_tmp <- formula2
        for (j in 1:(n_alt - 1))
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
      # estimate least squares regression
      data_regime <- data[ind_regime[[i]], ]
      model2_list[[i]] <- lm(formula2_tmp, data = data_regime, 
                             na.action = na.exclude)
      # save coefficients
      coef2_tmp <- coef(model2_list[[i]])
      n_coef2_tmp <- length(coef2_ind_regime[, i])
      coef_lambda_row[[i]] <- coef(model2_list[[i]])[-(1:n_coef2_tmp)]
      ind2_tmp <- coef2_ind_regime[, i]
      coef_lambda_tmp <- vector(mode = "list", length = n_alt - 1)
      counter_tmp <- length(ind2_tmp)
      for (t in 1:(n_alt - 1))
      {
        if (degrees[i, t] > 0)
        {
          coef_lambda_tmp[[t]] <- coef2_tmp[(counter_tmp + 1): 
                                            (counter_tmp + degrees[i, t])]
        }
        counter_tmp <- counter_tmp + degrees[i, t]
      }
      coef_lambda[[i]] <- coef_lambda_tmp
      names(coef_lambda[[i]]) <- alt_names[-n_alt]
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
                      nrow = length(ind_regime[[i]]), 
                      ncol = length(model1$par))
    for (t in 1:length(model1$par))
    {
      for (t1 in 1:(n_alt - 1))
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
      # adjust covariance matrix for heteroscedasticity
    sigma2_cond_triangle <- diag(resid2)
      # adjust covariance matrix for estimates from the first step
    sigma2_cond_triangle <- sigma2_cond_triangle + (dlambda %*% 
                                                    model1$cov %*% 
                                                    t(dlambda))
    cov_2step[[i]] <- X_2step_mat %*% 
                      sigma2_cond_triangle %*% 
                      t(X_2step_mat)
    }
  }

  # One step estimation procedure
  opt <- NULL
  if (estimator == "ml")
  {
    opt <- opt_switchSelection(opt_args = opt_args, control_lnL = control_lnL,
                               n_sim = n_sim, n_cores = n_cores, 
                               type = "mnprobit", start = start, 
                               opt_type = opt_type,
                               regularization = regularization)
    par <- opt$par
  }
  
  # Store parameters into variables
  par_list <- par_mnprobit(par = par, n_alt = n_alt, n_coef = n_coef, 
                           n_coef2 = n_coef2, n_regimes = n_regimes,
                           coef_ind_alt = coef_ind_alt, sigma_ind = sigma_ind, 
                           coef2_ind_regime = coef2_ind_regime, 
                           var2_ind_regime = var2_ind_regime, 
                           cov2_ind_regime = cov2_ind_regime,
                           coef_lambda = coef_lambda,
                           is2 = is2, degrees = degrees, 
                           estimator = estimator,
                           alt_names = alt_names, regimes_names = regimes_names,
                           colnames_W = colnames(W), colnames_X = colnames(X))
  par <- par_list$par
  coef <- par_list$coef
  coef2 <- par_list$coef2
  var2 <- par_list$var2
  cov2 <- par_list$cov2
  coef_lambda <- par_list$coef_lambda
  sigma <- par_list$sigma
  sigma_vec <- par_list$sigma_vec
  
  # Store variables in the output list
  out$par <- par
  out$coef <- coef
  if (!is.null(degrees))
  {
    out$coef2 <- coef2
    out$var2 <- var2
    out$cov2 <- cov2
  }
  out$coef_lambda <- coef_lambda
  out$sigma <- sigma
  out$other$sigma_vec <- sigma_vec

  # Estimate asymptotic covariance matrix
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
    # Covariances of the first step
    for (i in 1:(n_alt - 1))
    {
      cov[coef_ind_alt[, i],
          coef_ind_alt[, i]] <- model1$cov[model1$control_lnL$coef_ind_alt[, i] + 1,
                                           model1$control_lnL$coef_ind_alt[, i] + 1]
    }
    cov[sigma_ind, sigma_ind] <- model1$cov[model1$control_lnL$sigma_ind + 1,
                                            model1$control_lnL$sigma_ind + 1]
    # Covariances of the second step
    for (i in 1:n_regimes)
    {
      ind2_tmp <- coef2_ind_regime[, i]
      cov[ind2_tmp, ind2_tmp] <- cov_2step[[i]][1:length(ind2_tmp),
                                                1:length(ind2_tmp)]
      # if (cov_type == "parametric")
      # {
      #   cov[cov_y_ind[[1]][i, ],
      #       cov_y_ind[[1]][i, ]] <- cov_2step[[i]][-(1:length(ind2_tmp)),
      #                                              -(1:length(ind2_tmp))]
      # }
    }
  }
  
  # Create tables to store the results in a format used by summary function
  tbl_list <- tbl_mnprobit(par = par, n_par = n_par, cov = cov, 
                           n_alt = n_alt, n_coef = n_coef, 
                           n_coef2 = n_coef2, n_regimes = n_regimes,
                           coef_ind_alt = coef_ind_alt, 
                           sigma_ind = sigma_ind, sigma_ind_mat = sigma_ind_mat,
                           coef2_ind_regime = coef2_ind_regime, 
                           var2_ind_regime = var2_ind_regime,
                           cov2_ind_regime = cov2_ind_regime,
                           alt_names = alt_names, coef_lambda = coef_lambda,
                           is2 = is2, degrees = degrees, estimator = estimator,
                           cov_2step = cov_2step,
                           coef = coef, sigma_vec = sigma_vec, 
                           coef2 = coef2, var2 = var2, cov2 = cov2,
                           coef_lambda_row = coef_lambda_row,
                           regimes_names = regimes_names,
                           colnames_W = colnames(W), colnames_X = colnames(X))
  out$tbl <- tbl_list$tbl
  out$se <- tbl_list$se
  out$p_value <- tbl_list$p_value

  # Output list
  out$cov <- cov
  out$cov_type <- cov_type
  out$coef <- coef
  out$sigma = sigma
  out$logLik <- NA
  if (estimator != "2step")
  {
    if (is.null(regularization))
    {
      out$logLik <- opt$value
    }
    else
    {
      # because of regularization
      out$logLik <- lnL_mnprobit(par = par, control_lnL = control_lnL)
    }
  }
  out$W <- W
  out$z <- z
  out$n_obs <- n_obs
  out$formula <- formula
  out$formula2 <- formula2
  out$regimes <- regimes
  out$n_regimes <- n_regimes
  out$n_alt <- n_alt
  out$data <- data
  out$cov_2step <- cov_2step
  out$twostep <- model2_list
  out$estimator <- estimator
  out$degrees <- degrees
  out$coef_lambda <- coef_lambda
  if (!is.null(degrees))
  {
    out$degrees <- as.matrix(degrees)
    out$X <- X
    out$y <- y
  }
  out$alt_names <- alt_names
  out$regimes_names <- regimes_names
  out$other$is2 <- is2
  out$other$coef_lambda_row <- coef_lambda_row
  out$other$z_names <- colnames(W)
  out$other$y_names <- colnames(X)
  
  # Predict lambdas
  out$lambda <- predict(out, type = "lambda", alt = NULL)
  
  # return
  return(out)
}