# Differentiate function respect to parameters
# estimated via mvoprobit function
deriv_mvoprobit <- function(object, fn, fn_args = list(), 
                            eps = max(1e-4, sqrt(.Machine$double.eps) * 10),
                            type = "default")
{
  # Get some variables
  n_par <- object$control_lnL$n_par
  n_eq <- object$control_lnL$n_eq
  
  # Estimate value of the function at initial point
  fn_args$object <- object
  fn_val <- do.call(what = fn, args = fn_args)
  if (is.matrix(fn_val))
  {
    if (ncol(fn_val) > 1)
    {
      stop ("Function 'fn' should return a vector or a single column matrix.")
    }
  }
  n_val <- length(fn_val)
  
  # Prepare output matrix
  out <- matrix(NA, nrow = n_val, ncol = n_par)
  colnames(out) <- 1:ncol(out)
  
  # Differentiate respect to coefficients
  coef_ind <- lapply(object$control_lnL$coef_ind, function(x){x + 1})
  for (i in 1:n_eq)
  {
    for (j in 1:length(object$coef[[i]]))
    {
      # save initial value
      par_old_tmp <- object$coef[[i]][j] 
      # prepare increment
      eps_tmp <- eps * abs(par_old_tmp)
      # plus
      fn_args$object$coef[[i]][j] <- par_old_tmp + eps_tmp
      fn_plus <- do.call(what = fn, args = fn_args)
      # minus
      fn_args$object$coef[[i]][j] <- par_old_tmp - eps_tmp
      fn_minus <- do.call(what = fn, args = fn_args)
      # derivative
      out[, coef_ind[[i]][j]] <- (fn_plus - fn_minus) / (2 * eps_tmp)
      # names
      colnames(out)[coef_ind[[i]][j]] <- paste0("coef",  j, " of ",
                                                object$other$z_names[i])
      # set initial value
      fn_args$object$coef[[i]][j] <- par_old_tmp
    }
  }
  
  # Differentiate respect to cuts
  cuts_ind <- lapply(object$control_lnL$cuts_ind, function(x){x + 1})
  for (i in 1:n_eq)
  {
    for (j in 1:length(object$cuts[[i]]))
    {
      # save initial value
      par_old_tmp <- object$cuts[[i]][j]
      # prepare increment
      eps_tmp <- eps * abs(par_old_tmp)
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
  
  # Differentiate respect to coefficients of variance equation
  coef_var_ind <- lapply(object$control_lnL$coef_var_ind, function(x){x + 1})
  for (i in 1:n_eq)
  {
    if (object$control_lnL$is_het[i])
    {
      for (j in 1:length(coef_var_ind[[i]]))
      {
        # save initial value
        par_old_tmp <- object$coef_var[[i]][j]
        # prepare increment
        eps_tmp <- eps * abs(par_old_tmp)
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
  
  # Differentiate respect to elements of covariance matrix
  # of ordered equations
  if (n_eq > 1)
  {
    sigma_ind_mat <- object$control_lnL$sigma_ind_mat + 1
    for (i in 1:(n_eq - 1))
    {
      for (j in (i + 1):n_eq)
      {
        # save initial value
        par_old_tmp <- object$sigma[i, j]
        # prepare increment
        eps_tmp <- eps * abs(par_old_tmp)
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
  
  # Differentiate respect to parameters of marginal distribution
  if (object$other$is_marginal)
  {
    marginal_par_n <- object$control_lnL$marginal_par_n
    marginal_par_ind <- lapply(object$control_lnL$marginal_par_ind, 
                               function(x){x + 1})
    for (i in 1:n_eq)
    {
      if (marginal_par_n[i] > 0)
      {
        for (j in 1:marginal_par_n[i])
        {
          # save initial value
          par_old_tmp <- object$marginal_par[[i]][j]
          # prepare increment
          eps_tmp <- eps * abs(par_old_tmp)
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
  
  # Stuff for continuous equations
  if (object$other$is2)
  {
    n_eq2 <- object$control_lnL$n_eq2 
    n_regimes <- object$control_lnL$n_regimes
    
    # Differentiate respect to coefficients of continuous equations
    coef2_ind <- lapply(object$control_lnL$coef2_ind, function(x){x + 1})
    for (i in 1:n_eq2)
    {
      for (j in 1:n_regimes[i])
      {
        for (t in 1:length(object$coef2[[i]][j, ]))
        {
          # save initial value
          par_old_tmp <- object$coef2[[i]][j, t] 
          # prepare increment
          eps_tmp <- eps * abs(par_old_tmp)
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
    
    # Differentiate respect to variances of continuous equations
    var2_ind <- lapply(object$control_lnL$var2_ind, function(x){x + 1})
    for (i in 1:n_eq2)
    {
      for (j in 1:n_regimes[i])
      {
        # save initial value
        par_old_tmp <- object$var2[[i]][j] 
        # prepare increment
        eps_tmp <- eps * abs(par_old_tmp)
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
    
    # Differentiate respect to covariances
    # between continuous and ordered equations
    cov2_ind <- lapply(object$control_lnL$cov2_ind, function(x){x + 1})
    for (i in 1:n_eq2)
    {
      for (j in 1:n_regimes[i])
      {
        for (t in 1:n_eq)
        {
          # save initial value
          par_old_tmp <- object$cov2[[i]][j, t] 
          # prepare increment
          eps_tmp <- eps * abs(par_old_tmp)
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
    
    # Differentiate respect to covariances
    # between regime of continuous equations
    if (n_eq2 > 1)
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
            # save initial value
            par_old_tmp <- object$sigma2[[counter]][t]
            # prepare increment
            eps_tmp <- eps * abs(par_old_tmp)
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
  
  # Short version of the output
  if (type %in% c("grad", "derivative", "jac", "gradient", "jacobian"))
  {
    return(out)
  }

  out <- list(grad = out, val = fn_val)

  return(out)
}

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

# Differentiate function respect to parameters
# estimated via mnprobit function
deriv_mnprobit <- function(object, fn, fn_args = list(), 
                           eps = max(1e-4, sqrt(.Machine$double.eps) * 10),
                           type = "default")
{
  # Get some variables
  n_par <- object$control_lnL$n_par
  n_alt <- object$control_lnL$n_alt
  n_coef <- object$control_lnL$n_coef
  
  # Estimate value of the function at initial point
  fn_args$object <- object
  fn_val <- do.call(what = fn, args = fn_args)
  if (is.matrix(fn_val))
  {
    if (ncol(fn_val) > 1)
    {
      stop ("Function 'fn' should return a vector or a single column matrix.")
    }
  }
  n_val <- length(fn_val)
  
  # Prepare output matrix
  out <- matrix(NA, nrow = n_val, ncol = n_par)
  colnames(out) <- 1:ncol(out)
  
  # Differentiate respect to coefficients
  coef_ind_alt <- object$control_lnL$coef_ind_alt + 1
  for (i in 1:(n_alt - 1))
  {
    for (j in 1:n_coef)
    {
      # save initial value
      par_old_tmp <- object$coef[j, i] 
      # prepare increment
      eps_tmp <- eps * abs(par_old_tmp)
      # plus
      fn_args$object$coef[j, i] <- par_old_tmp + eps_tmp
      fn_plus <- do.call(what = fn, args = fn_args)
      # plus
      fn_args$object$coef[j, i] <- par_old_tmp - eps_tmp
      fn_minus <- do.call(what = fn, args = fn_args)
      # derivative
      out[, coef_ind_alt[j, i]] <- (fn_plus - fn_minus) / (2 * eps_tmp)
      # names
      colnames(out)[coef_ind_alt[j, i]] <- paste0("coef(", i, ",", j, ")")
      # set initial value
      fn_args$object$coef[j, i] <- par_old_tmp
    }
  }
  
  # Differentiate respect to covariances between the
  # alternatives of the multinomial equations
  if (n_alt > 2)
  {
    sigma_ind <- object$control_lnL$sigma_ind + 1
    counter <- 1
    for (i in 1:(n_alt - 1))
    {
      for (j in 1:i)
      {
        if (!((i == 1) & (j == 1)))
        {
          # save initial value
          par_old_tmp <- object$sigma[i, j]
          # prepare increment
          eps_tmp <- eps * abs(par_old_tmp)
          # plus
          fn_args$object$sigma[i, j] <- par_old_tmp + eps_tmp
          fn_args$object$sigma[j, i] <- fn_args$object$sigma[j, i]
          fn_plus <- do.call(what = fn, args = fn_args)
          # minus
          fn_args$object$sigma[i, j] <- par_old_tmp - eps_tmp
          fn_args$object$sigma[j, i] <- fn_args$object$sigma[j, i]
          fn_minus <- do.call(what = fn, args = fn_args)
          # derivative
          out[, sigma_ind[counter]] <- (fn_plus - fn_minus) / (2 * eps_tmp)
          # names
          colnames(out)[sigma_ind[counter]] <- paste0("cov(", i, ",", j, ")")
          # set initial value
          fn_args$object$sigma[i, j] <- par_old_tmp
          fn_args$object$sigma[j, i] <- par_old_tmp
          # counter
          counter <- counter + 1
        }
      }
    }
  }
  
  # Stuff for continuous equations
  is2 <- object$other$is2
  if (is2)
  {
    # Get some variables
    n_regimes <- object$n_regimes
      
    # Differentiate respect to coefficients
    coef2_ind_regimes <- object$control_lnL$coef2_ind_regime + 1
    for (i in 1:n_regimes)
    {
      for (j in 1:n_coef)
      {
        # save initial value
        par_old_tmp <- object$coef2[j, i] 
        # prepare increment
        eps_tmp <- eps * abs(par_old_tmp)
        # plus
        fn_args$object$coef2[j, i] <- par_old_tmp + eps_tmp
        fn_plus <- do.call(what = fn, args = fn_args)
        # minus
        fn_args$object$coef2[j, i] <- par_old_tmp - eps_tmp
        fn_minus <- do.call(what = fn, args = fn_args)
        # derivative
        out[, coef2_ind_regimes[j, i]] <- (fn_plus - fn_minus) / (2 * eps_tmp)
        # names
        colnames(out)[coef2_ind_regimes[j, i]] <- paste0("coef2(", i, 
                                                         ",", j, ")")
        # set initial value
        fn_args$object$coef2[j, i] <- par_old_tmp
      }
    }
    # Differentiate respect to variances
    var2_ind_regime <- object$control_lnL$var2_ind_regime + 1
    for (i in 1:n_regimes)
    {
      # save initial value
      par_old_tmp <- object$var2[i] 
      # prepare increment
      eps_tmp <- eps * abs(par_old_tmp)
      # plus
      fn_args$object$var2[i] <- par_old_tmp + eps_tmp
      fn_plus <- do.call(what = fn, args = fn_args)
      # minus
      fn_args$object$var2[i] <- par_old_tmp - eps_tmp
      fn_minus <- do.call(what = fn, args = fn_args)
      # derivative
      out[, var2_ind_regime[i]] <- (fn_plus - fn_minus) / (2 * eps_tmp)
      # names
      colnames(out)[var2_ind_regime[i]] <- paste0("var(", i - 1, ")")
      # set initial value
      fn_args$object$var2[i] <- par_old_tmp
    }

    # Differentiate respect to covariances between
    # multinomial and continuous equations
    cov2_ind_regime <- object$control_lnL$cov2_ind_regime + 1
    for (i in 1:n_regimes)
    {
      for (j in 1:(n_alt - 1))
      {
        # save initial value
        par_old_tmp <- object$cov2[j, i] 
        # prepare increment
        eps_tmp <- eps * abs(par_old_tmp)
        # plus
        fn_args$object$cov2[j, i] <- par_old_tmp + eps_tmp
        fn_args$object$coef_lambda[[i]][[j]] <- fn_args$object$cov2[j, i]
        fn_plus <- do.call(what = fn, args = fn_args)
        # minus
        fn_args$object$cov2[j, i] <- par_old_tmp - eps_tmp
        fn_args$object$coef_lambda[[i]][[j]] <- fn_args$object$cov2[j, i]
        fn_minus <- do.call(what = fn, args = fn_args)
        # derivative
        out[, cov2_ind_regime[j, i]] <- (fn_plus - fn_minus) / (2 * eps_tmp)
        # names
        colnames(out)[cov2_ind_regime[j, i]] <- paste0("cov2(", i - 1, 
                                                        ",", j, ")")
        # set initial value
        fn_args$object$cov2[j, i] <- par_old_tmp
        fn_args$object$coef_lambda[[i]][[j]] <- par_old_tmp
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

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

#' Delta method for mvoprobit and mnprobit functions
#' @description This function uses delta method to estimate standard errors
#' of functions of the estimator of the parameters of 
#' \code{mnprobit} and \code{mvoprobit}
#' functions if maximum-likelihood estimator has been used.
#' @param object an object of class 'mvoprobit' or 'mnprobit'.
#' @param fn function which returns a numeric vector and should depend on the 
#' elements of \code{object}. This elements should be accessed via 
#' \code{\link[switchSelection]{coef.mvoprobit}} and 
#' \code{\link[switchSelection]{coef.mnprobit}} functions. 
#' Also it is possible to use \code{\link[switchSelection]{predict.mvoprobit}} 
#' and \code{\link[switchSelection]{predict.mnprobit}} functions.
#' The first argument of \code{fn} should be \code{object}.
#' Therefore \code{coef} and \code{predict} functions in \code{fn} should also
#' depend on \code{object}.
#' @param fn_args list of additional arguments of \code{fn}.
#' @param eps positive numeric value representing the increment used for
#' numeric differentiation of \code{fn}.
#' @param cl numeric value between \code{0} and \code{1} representing
#' a confidence level of the confidence interval.
#' @details Numeric differentiation is used to estimate derivatives of
#' \code{fn} respect to various parameters of 
#' \code{\link[switchSelection]{mvoprobit}}  and 
#' \code{\link[switchSelection]{mnprobit}} functions.
#' 
#' This function may be used only if \code{object$estimator = "ml"}.
#' @return This function returns an object of class \code{delta_method}
#' that is a matrix which columns are as follows:
#' \itemize{
#' \item \code{val} - output of the \code{fn} function.
#' \item \code{se} - numeric vector such that \code{se[i]} represents standard
#' error associated with \code{val[i]}.
#' \item \code{p_value} - numeric vector such that \code{p_value[i]} represents
#' p-value of the two-sided significance test associated with \code{val[i]}.
#' \item \code{lwr} - realization of the lower (left) bound of the 
#' confidence interval.
#' \item \code{upr} - realization of the upper (right) bound of the 
#' confidence interval.
#' }
#' An object of class \code{delta_method} has implementation of 
#' \code{summary} method 
#' \code{\link[switchSelection]{summary.delta_method}}.
#' @template delta_method_examples_Template
delta_method <- function(object, fn, fn_args = list(), 
                        eps = max(1e-4, sqrt(.Machine$double.eps) * 10), 
                        cl = 0.95)
{
  if (object$estimator != "ml")
  {
    stop ("Available only for maximum-likelihood estimator.")
  }
  
  # Validate confidence level
  if ((cl <= 0) | (cl >= 1))
  {
    stop("Invalid 'cl' value. It should be between 0 and 1.")
  }
  
  # Get covariance matrix
  cov <- vcov(object)
  
  # Get values of derivatives
  d_list <- NULL
  if (is(object = object, class2 = "mvoprobit"))
  {
    d_list <- deriv_mvoprobit(object = object, fn = fn, 
                              fn_args = fn_args, eps = eps)
  }
  if (is(object = object, class2 = "mnprobit"))
  {
    d_list <- deriv_mnprobit(object = object, fn = fn, 
                             fn_args = fn_args, eps = eps)
  }
  d <- d_list$grad
  val <- d_list$val
  n <- nrow(d)
  if (is.null(names(val)))
  {
    names(val) <- 1:n
  }
  
  # Use delta method
  se <- vector(mode = "numeric", length = n)
  for (i in 1:n)
  {
    se[i] <- sqrt(d[i, , drop = FALSE] %*% cov %*% t(d[i, , drop = FALSE]))
  }
  
  # Get p-values
  z_value <- val / se
  p_value <- rep(NA, n)
  for (i in 1:n)
  {
    p_value[i] <- 2 * min(pnorm(z_value[i]), 1 - pnorm(z_value[i]))
  }
  
  # Get bounds of the confidence interval
  b <- qnorm(cl + (1 - cl) / 2) * se
  lwr <- val - b
  upr <- val + b
  
  # Aggregate the output
  names(se) <- names(val)
  names(p_value) <- names(val)
  out <- cbind(val, se, lwr, upr, p_value)
  colnames(out) <- c("val", "se", "lwr", "upr", "p_value")
  class(out) <- "delta_method"
  
  return(out)
}

#' Summary for an Object of Class delta_method
#' @description Provides summary for an object of class 'delta_method'.
#' @param object object of class 'delta_method'
#' @param ... further arguments (currently ignored)
#' @return Returns an object of class 'summary.delta_method'.
summary.delta_method <- function(object, ...)
{
  if (length(list(...)) > 0)
  {
    warning("Additional arguments passed through ... are ignored.")   
  }
  
  class(object) <- "summary.delta_method"
  
  return(object)
}

#' Print summary for an Object of Class delta_method
#' @description Prints summary for an object of class 'delta_method'.
#' @param x object of class 'delta_method'
#' @param ... further arguments (currently ignored)
#' @return The function returns input argument \code{x}.
print.summary.delta_method <- function(x, ...)
{
  if (length(list(...)) > 0)
  {
    warning("Additional arguments passed through ... are ignored.")   
  }
  
  class(x) <- "delta_method"
  
  printCoefmat(x, signif.legend = TRUE, has.Pvalue = TRUE)
  
  return(x)
}