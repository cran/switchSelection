# Differentiate function respect to the parameters
# estimated via msel function
deriv_msel <- function(object, 
                       fn, 
                       fn_args = list(), 
                       eps     = max(1e-4, sqrt(.Machine$double.eps) * 10),
                       type    = "default",
                       n_sim   = 10000,
                       n_cores = 1)
{
  # Set some properties of the object
  object$n_sim   <- n_sim
  object$n_cores <- n_cores
  
  # Get some variables
  n_par      <- object$other$n_par
  n_eq       <- object$other$n_eq
  n_eq2      <- object$other$n_eq2
  n_eq3      <- object$other$n_eq3
  is1        <- object$other$is1
  is2        <- object$other$is2
  is3        <- object$other$is3
  type3      <- object$type3
  sigma_omit <- object$other$sigma_omit
  
  # Deal with eps
  if (length(eps) == 1)
  {
    eps <- rep(eps, n_par)
  }
  
  # Estimate value of the function at initial point
  fn_args$object <- object
  fn_val         <- do.call(what = fn, args = fn_args)
  if (is.matrix(fn_val))
  {
    if (ncol(fn_val) > 1)
    {
      stop ("Function 'fn' should return a vector or a single column matrix.")
    }
  }
  n_val <- length(fn_val)
  
  # Prepare the output matrix
  out           <- matrix(0, nrow = n_val, ncol = n_par)
  colnames(out) <- 1:ncol(out)
  
  # Differentiate respect to the coefficients
  coef_ind <- object$ind$coef
  if (is1)
  {
    for (i in 1:n_eq)
    {
      for (j in 1:length(object$coef[[i]]))
      {
        if (eps[coef_ind[[i]][j]] != 0)
        {
          # save initial value
          par_old_tmp                     <- object$coef[[i]][j] 
          # prepare increment
          eps_tmp                         <- eps[coef_ind[[i]][j]] * 
                                             abs(par_old_tmp)
          # plus
          fn_args$object$coef[[i]][j]     <- par_old_tmp + eps_tmp
          fn_plus                         <- do.call(what = fn, args = fn_args)
          # minus
          fn_args$object$coef[[i]][j]     <- par_old_tmp - eps_tmp
          fn_minus                        <- do.call(what = fn, args = fn_args)
          # derivative
          out[, coef_ind[[i]][j]]         <- (fn_plus - fn_minus) / 
                                             (2 * eps_tmp)
          # names
          colnames(out)[coef_ind[[i]][j]] <- paste0("coef",  j, " of ",
                                                    object$other$z_names[i])
          # set initial value
          fn_args$object$coef[[i]][j]     <- par_old_tmp
        }
      }
    }
  }
  
  # Differentiate respect to the cuts
  cuts_ind <- object$ind$cuts
  if (is1)
  {
    for (i in 1:n_eq)
    {
      for (j in 1:length(object$cuts[[i]]))
      {
        if (eps[cuts_ind[[i]][j]] != 0)
        {
          # save initial value
          par_old_tmp <- object$cuts[[i]][j]
          # prepare increment
          eps_tmp <- eps[cuts_ind[[i]][j]] * abs(par_old_tmp)
          # plus
          fn_args$object$cuts[[i]][j] <- par_old_tmp + eps_tmp
          fn_plus <- do.call(what = fn, args = fn_args)
          # minus
          fn_args$object$cuts[[i]][j] <- par_old_tmp - eps_tmp
          fn_minus <- do.call(what = fn, args = fn_args)
          # derivative
          out[, cuts_ind[[i]][j]] <- (fn_plus - fn_minus) / (2 * eps_tmp)
          # names
          colnames(out)[cuts_ind[[i]][j]] <- paste0("cut",  j, " of ",
                                                    object$other$z_names[i])
          # set initial value
          fn_args$object$cuts[[i]][j] <- par_old_tmp
        }
      }
    }
  }
  
  # Differentiate respect to the coefficients of the variance equation
  coef_var_ind <- object$ind$coef_var
  if (is1)
  {
    for (i in 1:n_eq)
    {
      if (object$other$is_het[i])
      {
        for (j in 1:length(coef_var_ind[[i]]))
        {
          if (eps[coef_var_ind[[i]][j]] != 0)
          {
            # save initial value
            par_old_tmp <- object$coef_var[[i]][j]
            # prepare increment
            eps_tmp <- eps[coef_var_ind[[i]][j]] * abs(par_old_tmp)
            # plus
            fn_args$object$coef_var[[i]][j] <- par_old_tmp + eps_tmp
            fn_plus <- do.call(what = fn, args = fn_args)
            # minus
            fn_args$object$coef_var[[i]][j] <- par_old_tmp - eps_tmp
            fn_minus <- do.call(what = fn, args = fn_args)
            # derivative
            out[, coef_var_ind[[i]][j]] <- (fn_plus - fn_minus) / (2 * eps_tmp)
            # names
            colnames(out)[coef_var_ind[[i]][j]] <- paste0("coef_var", j, " of ",
                                                          object$other$z_names[i])
            # set initial value
            fn_args$object$coef_var[[i]][j] <- par_old_tmp
          }
        }
      }
    }
  }
  
  # Differentiate respect to the elements of the covariance matrix
  # of the ordinal equations
  if (n_eq > 1)
  {
    sigma_ind_mat <- object$ind$sigma_mat
    for (i in 1:(n_eq - 1))
    {
      for (j in (i + 1):n_eq)
      {
        if (eps[sigma_ind_mat[i, j]] != 0 & (sigma_omit[i, j] != 1))
        {
          # save initial value
          par_old_tmp <- object$sigma[i, j]
          # prepare increment
          eps_tmp <- eps[sigma_ind_mat[i, j]] * abs(par_old_tmp)
          # plus
          fn_args$object$sigma[i, j] <- par_old_tmp + eps_tmp
          fn_args$object$sigma[j, i] <- fn_args$object$sigma[i, j]
          fn_plus <- do.call(what = fn, args = fn_args)
          # minus
          fn_args$object$sigma[i, j] <- par_old_tmp - eps_tmp
          fn_args$object$sigma[j, i] <- fn_args$object$sigma[i, j]
          fn_minus <- do.call(what = fn, args = fn_args)
          # derivative
          out[, sigma_ind_mat[i, j]] <- (fn_plus - fn_minus) / (2 * eps_tmp)
          # names
          colnames(out)[sigma_ind_mat[i, j]] <- paste0("cov(",
                                                        object$other$z_names[i],
                                                        ",",
                                                        object$other$z_names[j],
                                                        ")")
          # set initial value
          fn_args$object$sigma[i, j] <- par_old_tmp
          fn_args$object$sigma[j, i] <- par_old_tmp
        }
      }
    }
  }
  
  # Differentiate respect to the parameters of the marginal distribution
  if (object$other$is_marginal)
  {
    marginal_par_n   <- object$other$marginal_par_n
    marginal_par_ind <- object$ind$marginal_par
    for (i in 1:n_eq)
    {
      if (marginal_par_n[i] > 0)
      {
        for (j in 1:marginal_par_n[i])
        {
          if (eps[marginal_par_ind[[i]][j]] != 0)
          {
            # save initial value
            par_old_tmp <- object$marginal_par[[i]][j]
            # prepare increment
            eps_tmp <- eps[marginal_par_ind[[i]][j]] * abs(par_old_tmp)
            # plus
            fn_args$object$marginal_par[[i]][j] <- par_old_tmp + eps_tmp
            fn_plus <- do.call(what = fn, args = fn_args)
            # minus
            fn_args$object$marginal_par[[i]][j] <- par_old_tmp - eps_tmp
            fn_minus <- do.call(what = fn, args = fn_args)
            # derivative
            out[, marginal_par_ind[[i]][j]] <- (fn_plus - fn_minus) / 
                                               (2 * eps_tmp)
            # set initial value
            fn_args$object$marginal_par[[i]][j] <- par_old_tmp
          }
        }
      }
    }
  }
  
  # Derivatives respect to the parameters of the continuous equations
  if (is2)
  {
    n_eq2     <- object$other$n_eq2 
    n_regimes <- object$other$n_regimes
    
    # Differentiate respect to the coefficients of the continuous equations
    coef2_ind <- object$ind$coef2
    for (i in 1:n_eq2)
    {
      for (j in 1:n_regimes[i])
      {
        for (t in 1:length(object$coef2[[i]][j, ]))
        {
          if (eps[coef2_ind[[i]][j, t]] != 0)
          {
            # save initial value
            par_old_tmp <- object$coef2[[i]][j, t] 
            # prepare increment
            eps_tmp <- eps[coef2_ind[[i]][j, t]] * abs(par_old_tmp)
            # plus
            fn_args$object$coef2[[i]][j, t] <- par_old_tmp + eps_tmp
            fn_plus <- do.call(what = fn, args = fn_args)
            # minus
            fn_args$object$coef2[[i]][j, t] <- par_old_tmp - eps_tmp
            fn_minus <- do.call(what = fn, args = fn_args)
            # derivative
            out[, coef2_ind[[i]][j, t]] <- (fn_plus - fn_minus) / (2 * eps_tmp)
            # names
            colnames(out)[coef2_ind[[i]][j, t]] <- paste0("coef2(", 
                                                          j - 1, ",", 
                                                          t, ")", " of ",
                                                          object$other$y_names[i])
            # set initial value
            fn_args$object$coef2[[i]][j, t] <- par_old_tmp
          }
        }
      }
    }
    
    # Differentiate respect to the variances of the continuous equations
    if (object$estimator == "ml")
    {
      var2_ind <- object$ind$var2
      for (i in 1:n_eq2)
      {
        for (j in 1:n_regimes[i])
        {
          if (eps[var2_ind[[i]][j]] != 0)
          {
            # save initial value
            par_old_tmp <- object$var2[[i]][j] 
            # prepare increment
            eps_tmp <- eps[var2_ind[[i]][j]] * abs(par_old_tmp)
            # plus
            fn_args$object$var2[[i]][j] <- par_old_tmp + eps_tmp
            fn_plus <- do.call(what = fn, args = fn_args)
            # plus
            fn_args$object$var2[[i]][j] <- par_old_tmp - eps_tmp
            fn_minus <- do.call(what = fn, args = fn_args)
            # derivative
            out[, var2_ind[[i]][j]] <- (fn_plus - fn_minus) / (2 * eps_tmp)
            # names
            colnames(out)[var2_ind[[i]][j]] <- paste0("var(", 
                                                       object$other$y_names[i],
                                                       ")", j - 1)
            # set initial value
            fn_args$object$var2[[i]][j] <- par_old_tmp
          }
        }
      }
    }
    
    # Differentiate respect to the covariances between the continuous 
    # and ordinal equations
    if ((object$estimator == "ml") & is1)
    {
      cov2_ind <- object$ind$cov2
      for (i in 1:n_eq2)
      {
        for (j in 1:n_regimes[i])
        {
          for (t in 1:n_eq)
          {
            if (eps[cov2_ind[[i]][j, t]] != 0)
            {
              # save initial value
              par_old_tmp <- object$cov2[[i]][j, t] 
              # prepare increment
              eps_tmp <- eps[cov2_ind[[i]][j, t]] * abs(par_old_tmp)
              # plus
              fn_args$object$cov2[[i]][j, t] <- par_old_tmp + eps_tmp
              fn_plus <- do.call(what = fn, args = fn_args)
              # minus
              fn_args$object$cov2[[i]][j, t] <- par_old_tmp - eps_tmp
              fn_minus <- do.call(what = fn, args = fn_args)
              # derivative
              out[, cov2_ind[[i]][j, t]] <- (fn_plus - fn_minus) / (2 * eps_tmp)
              # names
              colnames(out)[cov2_ind[[i]][j, t]] <- paste0("cov(", 
                                                            object$other$y_names[i],
                                                            "(", j, ")", ",", 
                                                            object$other$z_names[t],
                                                            ")")
              # set initial value
              fn_args$object$cov2[[i]][j, t] <- par_old_tmp
            }
          }
        }
      }
    }
    
    # Differentiate respect to the covariances between regimes of the 
    # continuous equations
    if (n_eq2 > 1 & (object$estimator == "ml"))
    {
      sigma2_ind <- lapply(object$control_lnL$sigma2_ind, function(x){x + 1})
      counter <- 1
      regimes_pair <- object$other$regimes_pair
      for (i in 2:n_eq2)
      {
        for (j in 1:(i - 1))
        {
          for (t in 1:nrow(regimes_pair[[counter]]))
          {
            if (eps[sigma2_ind[[counter]][t]] != 0)
            {
              # save initial value
              par_old_tmp <- object$sigma2[[counter]][t]
              # prepare increment
              eps_tmp <- eps[sigma2_ind[[counter]][t]] * abs(par_old_tmp)
              # plus
              fn_args$object$sigma2[[counter]][t] <- par_old_tmp + eps_tmp
              fn_plus <- do.call(what = fn, args = fn_args)
              # minus
              fn_args$object$sigma2[[counter]][t] <- par_old_tmp - eps_tmp
              fn_minus <- do.call(what = fn, args = fn_args)
              # derivative
              out[, sigma2_ind[[counter]][t]] <- (fn_plus - fn_minus) / 
                                                 (2 * eps_tmp)
              # names
              colnames(out)[sigma2_ind[[counter]]] <- paste0(
                  "cov(",
                   object$other$y_names[i],
                   "(",
                   regimes_pair[[counter]][t, 1],
                   ")", ",",
                   "(",
                   regimes_pair[[counter]][t, 2],
                   ")",
                   object$other$y_names[j],
                   ")")
              # set initial value
              fn_args$object$sigma2[[counter]][t] <- par_old_tmp
            }
            counter <- counter + 1
          }
        }
      }
    }
  }
  
  # Differentiate respect to the coefficients of the multinomial equations
  if (is3)
  {
    coef3_ind <- object$ind$coef3
    for (i in 1:(n_eq3 - 1))
    {
      for (j in 1:length(object$coef3[i, ]))
      {
        if (eps[coef3_ind[i, j]] != 0)
        {
          # save initial value
          par_old_tmp                     <- object$coef3[i, j]
          # prepare increment
          eps_tmp                         <- eps[coef3_ind[i, j]] * 
                                             abs(par_old_tmp)
          # plus
          fn_args$object$coef3[i, j]      <- par_old_tmp + eps_tmp
          fn_plus                         <- do.call(what = fn, args = fn_args)
          # minus
          fn_args$object$coef3[i, j]      <- par_old_tmp - eps_tmp
          fn_minus                        <- do.call(what = fn, args = fn_args)
          # derivative
          out[, coef3_ind[i, j]]          <- (fn_plus - fn_minus) / 
                                             (2 * eps_tmp)
          # names
          colnames(out)[coef3_ind[i, j]] <- paste0("coef",  j, 
                                                   " of alt", i)
          # set initial value
          fn_args$object$coef3[i, j]     <- par_old_tmp
        }
      }
    }
  }
  
  # Differentiate respect to the covariances between the alternatives of 
  # the multinomial equations
  if (is3 & (n_eq3 > 2) & (type3 == "probit"))
  {
    sigma3_ind <- object$ind$sigma3
    counter    <- 1
    for (i in seq_len(n_eq3 - 1))
    {
      for (j in 1:i)
      {
        if (!((i == 1) & (j == 1)))
        {
          if (eps[sigma3_ind[counter]] != 0)
          {
            # save initial value
            par_old_tmp <- object$sigma3[i, j]
            # prepare increment
            eps_tmp                            <- eps[sigma3_ind[counter]] * 
                                                  abs(par_old_tmp)
            # plus
            fn_args$object$sigma3[i, j]        <- par_old_tmp + eps_tmp
            fn_args$object$sigma3[j, i]        <- fn_args$object$sigma3[j, i]
            fn_plus                            <- do.call(what = fn, 
                                                          args = fn_args)
            # minus
            fn_args$object$sigma3[i, j]        <- par_old_tmp - eps_tmp
            fn_args$object$sigma3[j, i]        <- fn_args$object$sigma3[j, i]
            fn_minus                           <- do.call(what = fn, 
                                                          args = fn_args)
            # derivative
            out[, sigma3_ind[counter]]         <- (fn_plus - fn_minus) / 
                                                  (2 * eps_tmp)
            # names
            colnames(out)[sigma3_ind[counter]] <- paste0("cov(alt", i, 
                                                         ",alt", j, ")")
            # set initial value
            fn_args$object$sigma3[i, j]        <- par_old_tmp
            fn_args$object$sigma3[j, i]        <- par_old_tmp
          }
          counter <- counter + 1
        }
      }
    }
  }
  
  # Short version of the output
  if (type %in% c("grad", "derivative", "jac", "gradient", "jacobian"))
  {
    return(out)
  }

  out <- list(grad = out, val = fn_val)

  return(out)
}

#' Tests and confidence intervals for the parameters estimated by 
#' the msel function
#' @description This function conducts various statistical tests and calculates
#' confidence intervals for the parameters of the model estimated via the 
#' \code{\link[switchSelection]{msel}} function.
#' @template test_msel_param_Template
#' @template test_msel_details_Template
#' @template test_msel_return_Template
#' @template test_msel_references_Template
#' @template test_msel_examples_Template
test_msel <- function(object, 
                      fn, 
                      fn_args   = list(), 
                      test      = "t",
                      method    = "classic",
                      ci        = "classic",
                      cl        = 0.95,
                      se_type   = "dm",
                      trim      = 0,
                      vcov      = object$cov,
                      iter      = 100,
                      generator = rnorm,
                      bootstrap = NULL,
                      par_ind   = 1:object$control_lnL$n_par,
                      eps       = max(1e-4, sqrt(.Machine$double.eps) * 10),
                      n_sim     = 1000,
                      n_cores   = 1)
{
  # -------------------------------------------------------
  # Likelihood ratio test
  # -------------------------------------------------------
  
  if (length(object) == 2)
  {
    out <- lrtest_msel(model1 = object[[1]], model2 = object[[2]])
    return(out)
  }
  
  # -------------------------------------------------------
  # Setup
  # -------------------------------------------------------
  
  # List to aggregate the output
  out <- list()
  
  # Get the number of parameters of the model
  n_par <- object$control_lnL$n_par
  
  # Values for the output
  val     <- NULL
  stat    <- NULL
  se      <- NULL
  p_value <- NULL
  lwr     <- NULL
  upr     <- NULL
  
  # Provide the parameters
  object$n_sim   <- n_sim
  object$n_cores <- n_cores
  fn_args$object <- object
  
  # -------------------------------------------------------
  # Validation
  # -------------------------------------------------------
  
  # Validate test
  test     <- tolower(test)
  test_vec <- c("t", "wald")
  if (test %in% c("waldtest", "wald test", "aggregate", "wald-test"))
  {
    test <- "wald"
    warning("It is assumed that 'test' is 'wald' so Wald test is used.")
  }
  if (test %in% c("delta method", "deltamethod", "individual", "ind",
                  "t-test", "t_test"))
  {
    test <- "t"
    warning("It is assumed that 'test' is 't' so t-test is used.")
  }
  if (!(test %in% test_vec))
  {
    stop(paste0("Argument 'test' is wrong. ",
                "Please, insure that it is one of: ",
                paste0(test_vec, collapse = ", "),
                ".\n"))
  }
  
  # Validate method
  method     <- tolower(method)
  method_vec <- c("classic", "bootstrap", "score")
  if (method %in% c("asymptotic", "parametric", "normal"))
  {
    method <- "classic"
    warning("It is assumed that 'method' is 'classic'.")
  }
  if (method %in% c("boot", "nonparametric", "non-parametric",
                    "percentile"))
  {
    method <- "bootstrap"
    warning("It is assumed that 'method' is 'bootstrap'.")
  }
  if (method %in% c("score bootstrap", "bootstrap score"))
  {
    method <- "score"
    warning("It is assumed that 'method' is 'score'.")
  }
  if (!(method %in% method_vec))
  {
    stop(paste0("Argument 'method' is wrong. ",
                "Please, insure that it is one of: ",
                paste0(method_vec, collapse = ", "),
                ".\n"))
  }
  
  # Validate ci
  ci     <- tolower(ci)
  ci_vec <- c("classic", "percentile", "bc")
  if (ci %in% c("asymptotic", "normal"))
  {
    ci <- "classic"
    warning("It is assumed that 'ci' is 'classic'.")
  }
  if (ci %in% c("bootstrap", "bootstrap percentile", "percentile bootstrap"))
  {
    ci <- "percentile"
    warning("It is assumed that 'ci' is 'percentile'.")
  }
  if (ci %in% c("bias corrected", "bias-corrected", "bias_corrected",
                "Efron", "Efron1982", "Efron 1982"))
  {
    ci <- "bc"
    warning("It is assumed that 'ci' is 'bc'.")
  }
  if (!(ci %in% ci_vec))
  {
    stop(paste0("Argument 'ci' is wrong. ",
                "Please, insure that it is one of: ",
                paste0(ci_vec, collapse = ", "),
                ".\n"))
  }
  
  # Validate cl
  if (!is.numeric(cl) | (length(cl) > 1))
  {
    stop("Invalid 'cl' value. It should be a numeric value of length 1.")
  }
  if ((cl <= 0) | (cl >= 1))
  {
    stop("Invalid 'cl' value. It should be a numeric value between 0 and 1.")
  }
  
  # Validate se_type
  se_type     <- tolower(se_type)
  se_type_vec <- c("dm", "bootstrap")
  if (se_type %in% c("delta-method", "delta method", "delta_method"))
  {
    se_type <- "dm"
    warning("It is assumed that 'se_type' is 'dm'.")
  }
  if (se_type %in% c("boot"))
  {
    se_type <- "bootstrap"
    warning("It is assumed that 'se_type' is 'bootstrap'.")
  }
  if ((se_type == "bootstrap") & is.null(bootstrap))
  {
    stop("Invalide 'se_type' since 'bootstrap' is 'NULL'.")
  }
  if (!(se_type %in% se_type_vec))
  {
    stop(paste0("Argument 'se_type' is wrong. ",
                "Please, insure that it is one of: ",
                paste0(se_type_vec, collapse = ", "),
                ".\n"))
  }
  
  # Validate trim
  if ((trim < 0) | (trim >= 1) | !is.numeric(trim) | (length(trim) > 1))
  {
    stop("Parameter 'trim' should be a numeric value between 0 and 1")
  }
  
  # Validate vcov
  if (!is.matrix(vcov) | any(dim(vcov) != dim(object$cov)))
  {
    stop(paste0("Invalid 'vcov' value. ",
                "It should be a square matrix of appropriate size."))
  }
  
  # Validate iter
  iter <- round(iter)
  if (iter <= 0)
  {
    stop(paste0("Invalid 'iter' value. ",
                "It should be a positive integer."))
  }
  
  # Validate generator
  if (!is.function(generator))
  {
    stop(paste0("Invalid 'generator' value. ",
                "It should be a function."))
  }
  if (!("n" %in% names(formals(generator))))
  {
    stop(paste0("Function 'generator' should have 'n' input argument ",
                "responsible for the number of generated random variables."))
  }
  
  # Validate par_ind
  par_ind_unique <- unique(par_ind)
  if (length(par_ind) != length(par_ind_unique))
  {
    warning("Duplicates of some indexes have been removed from 'par_ind'.")
    par_ind <- par_ind_unique
  }
  out_of_range <- (par_ind < 1) | (par_ind > n_par)
  if (any(out_of_range))
  {
    warning(paste0("Indexes ", paste(par_ind[out_of_range], collapse = ", "), 
                   " are out of range so they have been removed",
                   " from 'par_ind'."))
    par_ind <- par_ind[!out_of_range]
  }
  if (length(par_ind) == 0)
  {
    stop("Input argument 'par_ind' should not be empty.")
  }
  
  # Validate eps
  if (length(eps) == 1)
  {
    eps <- rep(eps, length(par_ind))
  }
  if (length(eps) != length(par_ind))
  {
    stop(paste0("If 'eps' has more than 1 element then ",
                "length of 'eps' and 'par_ind' should be the same. "))
  }
  eps_tmp          <- rep(0, n_par)
  eps_tmp[par_ind] <- eps
  eps              <- eps_tmp
  eps[is.na(eps)]  <- 0
  
  # Validate object
  if (!is(object = object, class2 = "msel"))
  {
    stop("Invalid 'object' argument. It should be an object of class 'msel'.")
  }
  
  # Validate bootstrap
  is_bootstrap     <- !is.null(bootstrap)
  if (is_bootstrap)
  {
    if (!is(object = bootstrap, class2 = "bootstrap_msel"))
    {
      stop(paste0("Invalid 'bootstrap' argument. ",
                  "It should be an object of class 'bootstrap_msel'."))
    }
  }
  if (((method == "bootstrap")         | 
       (ci %in% c("percentile", "bc")) | 
       (se_type == "bootstrap"))        &
      !is_bootstrap)
  {
    stop("Argument 'bootstrap' should be provided.")
  }
  out$is_bootstrap <- is_bootstrap
  
  # Additional validation
  if ((test == "t") & (method == "score"))
  {
    stop("Score bootstrap is not available for the t-test.\n")
  }
  
  # -------------------------------------------------------
  # Values of the function and the derivatives
  # -------------------------------------------------------
  
  # Get the values of the function and its derivatives
  d        <- NULL  # derivatives
  val      <- NULL  # values of the function
  n_val    <- NULL  # dimensions of the function
  if ((method %in% c("classic", "score")) | (ci %in% c("classic")) |
      (se_type == "dm"))
  {
    d_list <- deriv_msel(object  = object,  fn  = fn, 
                         fn_args = fn_args, eps = eps)
    d      <- d_list$grad
    val    <- d_list$val
    n_val  <- nrow(d)
  }
  else
  {
    val   <- do.call(what = fn, args = fn_args)
    n_val <- length(val)
  }
  if (is.null(names(val)))
  {
    names(val) <- 1:n_val
  }
  
  # -------------------------------------------------------
  # Bootstrap values
  # -------------------------------------------------------
  
  val_bootstrap <- NULL   # values of the function for each bootstrap sample
  n_bootstrap   <- NULL   # the number of bootstrap iterations
  if (!is.null(bootstrap))
  {
    n_bootstrap   <- bootstrap$iter
    val_bootstrap <- matrix(NA, nrow = n_bootstrap, ncol = n_val)
    for (i in 1:n_bootstrap)
    {
      fn_args$object     <- update_msel(object, par = bootstrap$par[i, ])
      val_bootstrap[i, ] <- do.call(what = fn, args = fn_args)
    }
    fn_args$object <- object
  }
  out$n_bootstrap <- n_bootstrap
  
  # -------------------------------------------------------
  # Standard errors
  # -------------------------------------------------------
  
  # Calculate standard error via the delta method
  se     <- vector(mode = "numeric", length = n_val)
  fn_cov <- NULL
  if (se_type == "dm")
  {
    fn_cov <- d %*% vcov %*% t(d)
    se     <- sqrt(diag(fn_cov))
  }
  
  # Estimate standard error via the bootstrap
  if (se_type == "bootstrap")
  {
    # Apply the trimming if need
    val_bootsrap_adj <- val_bootstrap
    if(trim != 0)
    {
      trimmed                    <- rep(FALSE, n_bootstrap)
      val_bootstrap_adj          <- sweep(x      = val_bootstrap,           
                                          MARGIN = 2, 
                                          STATS  = colMeans(val_bootstrap), 
                                          FUN    = "-")
      trimmed_val                 <- rowSums(val_bootstrap_adj ^ 2)
      trimmed_crit                <- quantile(trimmed_val, probs = 1 - trim, 
                                              type = 1)
      trimmed                     <- trimmed_val > trimmed_crit
      val_bootsrap_adj[trimmed, ] <- 0 
    }
    
    # Calculate the standard errors
    fn_cov <- cov(val_bootsrap_adj)
    se     <- sqrt(diag(fn_cov))
  }
  
  # -------------------------------------------------------
  # Classic and bootstrap t-test
  # -------------------------------------------------------
  
  # Conduct t-test
  if (test == "t")
  {
    # Estimate test statistics
    stat <- val / se
    
    # Calculate p-values
    p_value <- rep(NA, n_val)
    if (method == "classic")
    {
      for (i in 1:n_val)
      {
        p_value[i] <- 2 * min(pnorm(stat[i]), 1 - pnorm(stat[i]))
      }
    }
    if (method == "bootstrap")
    {
      for (i in 1:n_val)
      {
        stat_bootstrap <- val_bootstrap[, i] / se[i]
        p_value[i]     <- mean(abs(stat_bootstrap - stat[i]) > abs(stat[i]))
      }
    }
  }
  
  # -------------------------------------------------------
  # Classic and bootstrap Wald test
  # -------------------------------------------------------
  
  # Conduct Wald test
  if ((test == "wald") & (method %in% c("classic", "bootstrap")))
  {
    # Prepare the values
    p_value <- NA
    stat    <- NA
    
    # Estimate test statistics
    val         <- matrix(val, ncol = 1)
    fn_cov_inv  <- qr.solve(fn_cov, tol = 1e-16)
    stat        <- as.numeric(t(val) %*% fn_cov_inv %*% val)
    
    # Calculate p-values
    if (method == "classic")
    {
      p_value <- 1 - pchisq(stat, df = n_val)
    }
    if (method == "bootstrap")
    {
      stat_bootstrap <- vector(mode = "numeric", length = n_bootstrap)
      for (i in 1:n_bootstrap)
      {
        val_diff          <- matrix(val_bootstrap[i, ] - val, ncol = 1)
        stat_bootstrap[i] <- as.numeric(t(val_diff) %*% 
                                        fn_cov_inv  %*% 
                                        val_diff)
      }
      p_value <- mean(stat_bootstrap > stat)
    }
  }
  
  # -------------------------------------------------------
  # Score bootstrap Wald test
  # -------------------------------------------------------

  # Conduct score bootstrap Wald test
  if ((test == "wald") & (method == "score"))
  {
    # Calculate Jacobian and Hessian if need
    if (!hasName(object, name = "J") | !hasName(object, name = "H"))
    {
      vcov_ml_list  <- vcov_ml(object, type = "sandwich", 
                               n_cores = 1, n_sim = 1000)
      object$J      <- vcov_ml_list$J
      object$H      <- vcov_ml_list$H
    }
    # initial statistic
    scores <- object$J
    H      <- object$H
    n_obs  <- nrow(scores)
    val    <- matrix(val, ncol = 1) / sqrt(n_obs)
    H_inv  <- qr.solve(H, tol = 1e-16)
    A      <- d %*% H_inv
    stat   <- t(val) %*% 
              qr.solve(A %*% cov(scores) %*% t(A), tol = 1e-16) %*% 
              val
    stat    <- as.numeric(stat)
    # bootstrapped statistics 
    stat_bootstrap <- rep(NA, iter)
    for (i in 1:iter)
    {
      weights            <- generator(n = nrow(scores))
      scores_adj         <- weights * scores
      S                  <- A %*% colSums(scores_adj) / sqrt(n_obs)
      stat_bootstrap[i]  <- t(S) %*% qr.solve(A %*% cov(scores_adj) %*% t(A), 
                                              tol = 1e-16) %*% S
    }
    p_value <- mean(stat_bootstrap >= stat)
  }
  
  # -------------------------------------------------------
  # Confidence intervals
  # -------------------------------------------------------
  
  is_ci     <- FALSE
  if (test == "t")
  {
    is_ci      <- TRUE
    cl_p_lower <- (1 - cl) / 2
    cl_p_upper <- 1 - cl_p_lower
    
    # Classic asymptotic confidence interval
    if (ci == "classic")
    {
      z_q <- qnorm(cl + (1 - cl) / 2) * se
      lwr <- val - z_q
      upr <- val + z_q
    }
    
    # Percentile bootstrap confidence interval
    if (ci == "percentile")
    {
      lwr      <- vector(mode = "numeric", length = n_val)
      upr      <- lwr
      for (i in 1:n_val)
      {
        lwr[i]   <- quantile(val_bootstrap[, i], probs = cl_p_lower, type = 1)
        upr[i]   <- quantile(val_bootstrap[, i], probs = cl_p_upper, type = 1)
      }
    }
    
    # Bias-corrected percentile bootstrap confidence interval
    if (ci == "bc")
    {
      lwr <- vector(mode = "numeric", length = n_val)
      upr <- lwr
      for (i in 1:n_val)
      {
        adj_tmp     <- 2 * qnorm(mean(val_bootstrap[, i] <= val[i]))
        cl_bc_lower <- pnorm(qnorm(cl_p_lower) + adj_tmp)
        cl_bc_upper <- pnorm(qnorm(cl_p_upper) + adj_tmp)
        lwr[i]      <- quantile(val_bootstrap[, i], 
                                probs = cl_bc_lower, type = 1)
        upr[i]      <- quantile(val_bootstrap[, i], 
                                probs = cl_bc_upper, type = 1)
      }
    }
  }
  out$is_ci  <- is_ci
  
  # -------------------------------------------------------
  # Output
  # -------------------------------------------------------
  
  # Aggregate the output
  if ((test == "wald") | (n_val == 1))
  {
    out$tbl <- data.frame(row.names = "test")
  }
  else
  {
    out$tbl <- data.frame(row.names = 1:n_val)
  }
  out$tbl$stat <- stat

  if (test == "t")
  {
    out$tbl$val       <- val
    out$tbl$se        <- se
    rownames(out$tbl) <- 1:n_val
  }

  out$test    <- test
  out$method  <- method
  out$se_type <- se_type
  out$ci      <- ci
  out$cl      <- cl
  out$n_val   <- n_val
  if (is_ci)
  {
    out$tbl$lwr <- lwr
    out$tbl$upr <- upr
  }
  out$tbl$p_value <- p_value

  # Some other output
  if (method == "score")
  {
    out$iter <- iter
  }
  
  # Return the results
  class(out)  <- "test_msel"
  
  return(out)
}

#' Summary for an Object of Class delta_method
#' @description Provides summary for an object of class 'delta_method'.
#' @param object object of class 'delta_method'
#' @param ... further arguments (currently ignored)
#' @param is_legend a logical; if \code{TRUE} then additional information
#' is shown.
#' @return Returns an object of class 'summary.delta_method'.
summary.test_msel <- function(object, ..., is_legend = TRUE)
{
  if (length(list(...)) > 0)
  {
    warning("Additional arguments passed through ... are ignored.")   
  }
  
  class(object)                         <- "summary.test_msel"
  attr(x = object, which = "is_legend") <- is_legend
  
  return(object)
}

#' Print summary for an Object of Class test_msel
#' @description Prints summary for an object of class 'test_msel'.
#' @param x object of class 'test_msel'
#' @param ... further arguments (currently ignored)
#' @param is_legend a logical; if \code{TRUE} then additional information
#' is shown.
#' @return The function returns input argument \code{x}.
print.summary.test_msel <- function(x, ..., is_legend = TRUE)
{
  if (length(list(...)) > 0)
  {
    warning("Additional arguments passed through ... are ignored.")   
  }
  
  # Assign the class
  class(x) <- "test_msel"
  
  # Assign the is_legend value
  if (!is.null(attr(x = x, "is_legend")))
  {
    is_legend <- attr(x = x, which = "is_legend")
  }
  
  # Header
  if (x$test == "t")
  {
    cat("The results of the t-test\n")
  }
  if (x$test == "wald")
  {
    if (x$method == "score")
    {
      cat("The results of the score bootstrap Wald test\n")
    }
    else
    {
      cat("The results of the Wald test\n")
    }
  }
  cat("---\n")
  
  # Results
  printCoefmat(x$tbl, signif.legend = TRUE, has.Pvalue = TRUE)
  
  if (!is_legend)
  {
    return(x)
  }
  
  # Footnote
  cat("---\n")
  if (x$test == "t")
  {
    cat("H0: val = 0\n")
  }
  if (x$test == "wald")
  {
    cat("H0: val = (0,...,0)\n")
  }
  
  # method
  if (x$method == "classic")
  {
    if (x$test == "t")
    {
      cat(paste0("p-values are calculated under the assumption of the \n",
                 "asymptotic normality of the test statistic.\n"))
    }
    if (x$test == "wald")
    {
      cat(paste0("p-values are calculated under the assumption that \n",
                 "asymptotic distribution of the test statistic is \n",
                 "chi-squared with ", x$n_val, " degrees of freedom.\n"))
    }
  }
  if (x$method == "bootstrap")
  {
    cat("p-values are calculated via the bootstrap.\n")
  }
  
  # se_type
  if (x$se_type == "dm")
  {
    if (x$test == "t")
    {
      cat("Standard errors are estimated via the delta method.\n")
    }
    if (x$test == "wald")
    {
      cat(paste0("Covariance matrix of the outputs of the 'fn'\n",
                 "has been estimated via the delta method.\n"))
    }
  }
  if (x$se_type == "bootstrap")
  {
    if (x$test == "t")
    {
      cat("Standard errors are estimated via the bootstrap.\n")
    }
    if (x$test == "wald")
    {
      cat(paste0("Covariance matrix of the outputs of the 'fn'\n",
                 "has been estimated via the bootstrap.\n"))
    }
  }
  
  # ci
  if (x$is_ci)
  {
    if (x$ci == "classic")
    {
      cat(paste0("Classic ", x$cl * 100, "% confidence interval is used.\n"))
    }
    if (x$ci == "percentile")
    {
      cat(paste0("Percentile bootstrap ", x$cl * 100, 
                 "% confidence interval is used.\n"))
    }
    if (x$ci == "bc")
    {
      cat(paste0("Bias corrected percentile bootstrap ", x$cl * 100, 
                 "% confidence interval is used.\n"))
    }
  }
  
  # bootstrap
  if (x$is_bootstrap)
  {
    cat(paste0("The number of the bootstrap iterations is ", 
               x$n_bootstrap, ".\n"))
  }
  if (x$method == "score")
  {
    cat(paste0("The number of the score bootstrap iterations is ", 
               x$iter, ".\n"))
  }
  
  return(x)
}