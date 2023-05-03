#' Multinomial probit model
#' @description This function estimates parameters of multinomial probit
#' model and sample selection model with continuous outcome and multinomial
#' probit selection mechanism.
#' @template mnprobit_details_Template
#' @template mnprobit_param_Template
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
                     control = NULL)
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
  if (as.character(class(formula)) == "mvoprobit")
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
    regimes <- regimes <- rep(-1, max(df[, 1]))
  }
  else
  {
    if ((length(regimes) == 1) | is.null(regimes))
    {
      regimes <- rep(0, max(df[, 1]))
    }
  }
  n_regimes <- ifelse(is2, max(regimes) + 1, 0)
  
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
  W <- as.matrix(cbind(1, df[, -1]))                       # regressors
  colnames(W)[1] <- "(Intercept)"
    # continuous equation
  y <- numeric()
  X <- matrix()
  if (is2)
  {
    y <- df2[, 1]
    X <- as.matrix(cbind(1, df2[, -1]))
  }
  
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
  coef_ind_alt <- matrix(NA, nrow = n_coef, ncol = n_alt - 1)
  for (i in 1:(n_alt - 1)) 
  {
    coef_ind_alt[, i] <- ((i - 1) * n_coef + 1):(i * n_coef)
  }
  
  # Get indexes of covariances between alternatives
  sigma_ind <- (n_coef_total + 1):n_par
  sigma_ind_mat <- matrix(NA, nrow = n_sigma, ncol = 2)
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
  coef2_ind_regime <- matrix()
  cov2_ind_regime <- matrix()
  var2_ind_regime <- numeric()
  if (is2)
  {
    coef2_ind_regime <- matrix(NA, nrow = n_coef2, ncol = n_regimes)
    cov2_ind_regime <- matrix(NA, nrow = n_alt - 1, ncol = n_regimes)
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
                      n_coef_total = n_coef_total,
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
  
  # Convert degrees to matrix
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

  # Two-step estimation procedure
  par <- NULL
  cov_2step <- NULL
  model2_list <- NULL
  coef_lambda <- vector(mode = "list", length = n_regimes)
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
    opt <- optim(par = start,
                 method = "BFGS",
                 fn = lnL_mnprobit,
                 gr = grad_mnprobit,
                 control = list(maxit = 1000000,
                                fnscale = -1,
                                reltol = 1e-10),
                 n_sim = n_sim, n_cores = n_cores,
                 control_lnL = control_lnL)
    
    # Add genetic optimization if need
    if (opt_type == "gena")
    {
      opt <- gena::gena(fn = lnL_mnprobit,
                        gr = grad_mnprobit,
                        pop.initial = opt$par,
                        mutation.method = "percent", 
                        maxiter = 100, info = TRUE,
                        lower = -2 * abs(opt$par), 2 * abs(opt$par),
                        n_sim = n_sim, n_cores = n_cores,
                        control_lnL = control_lnL)
    }
    
    # Store parameters estimates
    par <- opt$par
  }
  
  # Store the coefficients
  coef <- matrix(NA, nrow = n_coef, ncol = n_alt - 1)
  for (i in 1:(n_alt - 1))
  {
    coef[, i] <- par[coef_ind_alt[, i]]
  }
  
  # Store the covariance matrix
  sigma <- matrix(NA, nrow = n_alt - 1, ncol = n_alt - 1)
  sigma[1, 1] = 1;
  sigma_vec <- NULL
  if (n_alt > 2)
  {
    counter <- 1;
    sigma_vec <- par[sigma_ind]
    for (i in 1:(n_alt - 1))
    {
      for (j in 1:i)
      {
        if (!((i == 1) & (j == 1)))
        {
          sigma[i, j] <- sigma_vec[counter]
          sigma[j, i] <- sigma[i, j];
          counter <- counter + 1
        }
      }
    }
  }
  
  # Store parameters of continuous equation
  coef2 <- NULL
  var2 <- NULL
  cov2 <- NULL
  if (is2)
  {
    coef2 <- matrix(NA, ncol = n_regimes, nrow = n_coef2)
    var2 <- vector(mode = "numeric", length = n_regimes)
    cov2 <- matrix(NA, nrow = n_alt - 1, ncol = n_regimes)
    for (i in 1:n_regimes)
    {
      coef2[, i] <- par[coef2_ind_regime[, i]]
      var2[i] <- par[var2_ind_regime[i]]
      cov2[, i] <- par[cov2_ind_regime[, i]]
      if (estimator == "ml")
      {
        coef_lambda[[i]] <- vector(mode = "list", length = n_alt - 1)
        for (j in 1:(n_alt - 1))
        {
          coef_lambda[[i]][[j]] <- cov2[j, i]
        }
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
        H <- gena::gena.hessian(gr = lnL_mnprobit,
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
        J <- lnL_mnprobit(par = par,
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
  tbl_coef <- vector(mode = "list", length = n_alt - 1)
  tbl_sigma <- matrix()
  tbl_coef2 <- vector(mode = "list", length = n_regimes)
  tbl_var2 <- vector(mode = "list", length = n_regimes)
  tbl_cov2 <- vector(mode = "list", length = n_regimes)
  # coefficients
  for (i in 1:(n_alt - 1))
  {
    tbl_coef[[i]] <- as.matrix(cbind(Estimate = coef[, i],
                                     Std_Error = se[coef_ind_alt[, i]],
                                     z_value = z_value[coef_ind_alt[, i]],
                                     p_value = p_value[coef_ind_alt[, i]]))
    rownames(tbl_coef[[i]]) <- colnames(W)
    names(coef[, i]) <- colnames(W)
  }
  # covariances
  if (n_alt > 2)
  {
    sigma_vec_names <- paste0("cov(", sigma_ind_mat[, 1, drop = FALSE], 
                              ",",
                                       sigma_ind_mat[, 2, drop = FALSE],
                              ")")
    tbl_sigma <- as.matrix(cbind(Estimate = sigma_vec,
                                 Std_Error = se[sigma_ind],
                                 z_value = z_value[sigma_ind],
                                 p_value = p_value[sigma_ind]))
    rownames(tbl_sigma) <- sigma_vec_names
  }
  # Continuous equation
  if (is2)
  {
    for (i in 1:n_regimes)
    {
      #coefficients
      tbl_coef2[[i]] <- as.matrix(cbind(Estimate = coef2[, i],
                                        Std_Error = se[coef2_ind_regime[, i]],
                                        z_value = z_value[coef2_ind_regime[, i]],
                                        p_value = p_value[coef2_ind_regime[, i]]))
      rownames(tbl_coef2[[i]]) <- colnames(X)
      names(coef2[, i]) <- colnames(X)
      # variances
      tbl_var2[[i]] <- as.matrix(cbind(Estimate = var2[i],
                                       Std_Error = se[var2_ind_regime[i]],
                                       z_value = z_value[var2_ind_regime[i]],
                                       p_value = p_value[var2_ind_regime[i]]))
      rownames(tbl_var2[[i]]) <- "var"
      # covariances
      tbl_cov2[[i]] <- as.matrix(cbind(Estimate = cov2[, i],
                                       Std_Error = se[cov2_ind_regime[, i]],
                                       z_value = z_value[cov2_ind_regime[, i]],
                                       p_value = p_value[cov2_ind_regime[, i]]))
      rownames(tbl_cov2[[i]]) <- paste0("cov(", 1:(n_alt - 1), ")")
    }
  }
    # standard errors for lambdas
  if (estimator == "2step")
  {
    tbl_lambda <- vector(mode = "list", length = n_regimes[1])
    for (i in 1:n_regimes[1])
    {
      n_coef2_tmp <- length(coef2_ind_regime[, i])
      estimate_lambda <- coef(model2_list[[i]])[-(1:n_coef2_tmp)]
      se_lambda <- sqrt(diag(cov_2step[[i]][-(1:n_coef2_tmp), 
                                            -(1:n_coef2_tmp), 
                                            drop = FALSE]))
      z_lambda <- estimate_lambda / se_lambda
      n_lambda <- length(estimate_lambda)
      p_value_lambda <- rep(NA, n_lambda)
      for (j in 1:n_lambda)
      {
        p_value_lambda[j] <- 2 * min(pnorm(z_lambda[j]), 1 - pnorm(z_lambda[j]))
      }
      tbl_lambda[[i]] <- as.matrix(cbind(Estimate = estimate_lambda,
                                         Std_Error = se_lambda,
                                         z_value = z_lambda,
                                         p_value = p_value_lambda))
      rownames(tbl_lambda[[i]]) <- names(estimate_lambda)
    }
    names(tbl_lambda) <- paste0("regime ", 0:(n_regimes- 1))
    out$tbl$lambda <- tbl_lambda
  }
  # Save the tables
  out$tbl <- list(coef = tbl_coef,
                  sigma = tbl_sigma,
                  coef2 = tbl_coef2,
                  var2 = tbl_var2,
                  cov2 = tbl_cov2)
  
  # Standard errors for lambdas
  if (estimator == "2step")
  {
    tbl_lambda <- vector(mode = "list", length = n_regimes)
    for (i in 1:n_regimes)
    {
      n_coef2_tmp <- length(coef2_ind_regime[, i])
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
    names(tbl_lambda) <- paste0("regime ", 0:(n_regimes - 1))
    out$tbl$lambda <- tbl_lambda
  }

  # Output list
  out$par <- par
  out$cov <- cov
  out$coef <- coef
  out$sigma = sigma
  out$logLik <- NA
  if (estimator != "2step")
  {
    out$logLik <- opt$value
  }
  out$W <- W
  out$z <- z
  out$control_lnL <- control_lnL
  out$n_obs <- n_obs
  out$formula <- formula
  out$formula2 <- formula2
  out$regimes <- regimes
  out$n_regimes <- n_regimes
  out$n_alt <- n_alt
  out$p_value <- p_value
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
    out$coef2 <- coef2
    out$var2 <- var2
    out$cov2 <- cov2
  }
  class(out) <- "mnprobit"
  
  # Predict lambdas
  out$lambda <- predict(out, type = "lambda", alt = NULL)
  
  # return
  return(out)
}

#' Multinomial probit model
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
                             control = list())
{   
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

#' Extract the Number of Observations from a Fit of the mnprobit Function.
#' @description Extract the number of observations from a model fit
#' of the \code{\link[switchSelection]{mnprobit}} function.
#' @param object object of class "mnprobit"
#' @param ... further arguments (currently ignored)
#' @details Unobservable values of continuous equations are included into
#' the number of observations.
#' @return A single positive integer number.
nobs.mnprobit <- function(object, ...)
{
  if (length(list(...)) > 0)
  {
    warning("Additional arguments passed through ... are ignored.")   
  }
  
  return(object$control_lnL$n_obs)
}

#' Extract Log-Likelihood from a Fit of the mnprobit Function. 
#' @description Extract Log-Likelihood from a model fit
#' of the \code{\link[switchSelection]{mnprobit}} function.
#' @param object object of class "mnprobit"
#' @param ... further arguments (currently ignored)
#' @details If \code{estimator == "2step"} in 
#' \code{\link[switchSelection]{mnprobit}} then function may return
#' \code{NA} value since two-step estimator of covariance matrix may be
#' not positively defined.
#' @return Returns an object of class 'logLik'.
logLik.mnprobit <- function(object, ...)
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

#' Summary for an Object of Class mnprobit
#' @description Provides summary for an object of class 'mnprobit'.
#' @param object object of class "mnprobit"
#' @param ... further arguments (currently ignored)
#' @details This function just changes the class of the 'mnprobit'
#' object to 'summary.mnprobit'.
#' @return Returns an object of class 'summary.mnprobit'.
summary.mnprobit <- function(object, ...) 
{
  if (length(list(...)) > 0)
  {
    warning("Additional arguments passed through ... are ignored.")   
  }
  
  class(object) <- "summary.mnprobit";
  
  return(object)
}

#' Print summary for an Object of Class mnprobit
#' @description Prints summary for an object of class 'mnprobit'.
#' @param x object of class "mnprobit"
#' @param ... further arguments (currently ignored)
#' @return The function returns input argument \code{x} changing
#' it's class to \code{lrtest}.
print.summary.mnprobit <- function(x, ...)
{
  class(x) <- "mnprobit"
  estimator <- x$estimator
  
  tbl_coef <- x$tbl$coef
  tbl_sigma <- x$tbl$sigma
  tbl_coef2 <- x$tbl$coef2
  tbl_var2 <- x$tbl$var2
  tbl_cov2 <- x$tbl$cov2
  tbl_lambda <- NULL
  if (estimator == "2step")
  {
    tbl_lambda <- x$tbl$lambda
  }
  
  n_regimes <- x$n_regimes
  n_alt <- x$n_alt
  is2 <- n_regimes > 0
  
  coef_ind_alt <- x$control_lnL$coef_ind_alt + 1
  sigma_ind <- x$control_lnL$sigma_ind + 1
  coef2_ind_regime <- x$control_lnL$coef2_ind_regime + 1
  var2_ind_regime <- x$control_lnL$var2_ind_regime + 1
  cov2_ind_regime <- x$control_lnL$cov2_ind_regime + 1
  
  if (is2)
  {
    cat("Model with multinomial probit selection mechanism\n")
  }
  else
  {
    cat("Multinomial probit model\n")
  }
  cat(paste0("There are ", n_alt, " alternatives",
             ifelse(is2, paste0(" and ", n_regimes, " regimes\n"), ".")))
  cat("---\n")
  if (estimator == "ml")
  {
    if (!is.na(x$logLik))
    {
      cat(paste0("Log-likelihood = ", round(x$logLik, 4), "\n"))
    }
    if (!is.na(AIC(x)))
    {
      cat(paste0("AIC = ", round(AIC(x), 4), "\n"))
    }
  }
  cat(paste0("Observations = ", nobs(x), "\n"))
  
  # stars
  stars <- switchSelection::starsVector(x$p_value)
  
  # coefficients
  for (i in 1:(n_alt - 1))
  {
    cat("---\n")
    cat(paste0("Alternative ", i, "\n"))
    cat("-\n");
    cat("Coefficients: \n")
    print(as.table(cbind(round(tbl_coef[[i]], 4), 
                         stars[coef_ind_alt[, i]])))
  }
  
  # covariances
  cat("---\n")
  cat("Covarianes between alternatives: \n")
  print(as.table(cbind(round(tbl_sigma, 4), 
                        stars[sigma_ind])))
  
  if (is2)
  {
    for (i in 1:n_regimes)
    {
      cat("---\n")
      cat(paste0("Regime ", i - 1, "\n"))
      cat("-\n");
      cat("Coefficients: \n")
      print(as.table(cbind(round(tbl_coef2[[i]], 4), 
                           stars[coef2_ind_regime[, i]])))
      if (estimator == "ml")
      {
        cat("-\n");
        cat("Variance: \n")
        print(as.table(cbind(round(tbl_var2[[i]], 4), 
                             stars[var2_ind_regime[i]])))
        cat("-\n");
        cat("Covariances with alternatives: \n")
        print(as.table(cbind(round(tbl_cov2[[i]], 4), 
                             stars[cov2_ind_regime[, i]])))
      }
      else
      {
        stars_lambda <- switchSelection::starsVector( 
          tbl_lambda[[i]][, "p_value"])
        cat("-\n");
        cat(paste0("Selectivity correction terms:\n"))
        print(as.table(cbind(round(tbl_lambda[[i]], 4), 
                             as.vector(stars_lambda))))
      }
    }
  }
  
  cat("---\n")
  cat("Signif. codes:  0 '***' 0.01 '**' 0.05 '*' 0.1 \n")
  
  return(x)
}

#' Print for an Object of Class mnprobit
#' @description Prints summary for an object of class 'mnprobit'.
#' @param x object of class "mnprobit"
#' @param ... further arguments (currently ignored)
#' @return The function returns input argument \code{x}.
print.mnprobit <- function(x, ...)
{
  print.summary.mnprobit(x)
  return(x)
}