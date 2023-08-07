#' Summary for an Object of Class mvoprobit
#' @description Provides summary for an object of class 'mvoprobit'.
#' @param object object of class 'mvoprobit'
#' @param ... further arguments (currently ignored)
#' @param vcov positively defined numeric matrix representing
#' asymptotic variance-covariance matrix of the estimator to be
#' used for calculation of standard errors and p-values. It may also be a 
#' character. Then \code{\link[switchSelection]{vcov.mvoprobit}} function
#' will be used which input argument \code{type} will be set to \code{vcov}.
#' If \code{estimator = "2step"} then \code{vcov} should be an estimate of the 
#' asymptotic covariance matrix of the first step estimator.
#' @param show_ind logical; if \code{TRUE} then indexes of parameters will be
#' shown. Particularly, these indexes may be used in \code{ind} element of
#' \code{regularization} parameter of \code{\link[switchSelection]{mvoprobit}}.
#' @details If \code{vcov} is \code{NULL} then this function just changes the 
#' class of the 'mvoprobit' object to 'summary.mvoprobit'. Otherwise it 
#' additionally changes \code{object$cov} to \code{vcov} and use it to
#' recalculate \code{object$se}, \code{object$p_value} and \code{object$tbl} 
#' values. It also adds the value of \code{ind} argument to the object.
#' @return Returns an object of class 'summary.mvoprobit'.
summary.mvoprobit <- function(object, ..., vcov = NULL, show_ind = FALSE) 
{
  if (length(list(...)) > 0)
  {
    warning("Additional arguments passed through ... are ignored.")   
  }
  
  # Add information on the indexes
  object$show_ind <- show_ind
  
  if (object$estimator == "2step")
  {
    if (!is.null(vcov))
    {
      warning("Argument 'vcov' is ignored for two-step estimator.")
    }
    vcov <- NULL
  }
  
  if (!is.null(vcov))
  {
    if (is.character(vcov))
    {
      vcov <- vcov.mvoprobit(object, type = vcov)
    }
    if (any(dim(vcov) != dim(object$cov)))
    {
      stop("Incorrect size of 'vcov' matrix.")
    }
    object$cov <- vcov
    tbl_list <- tbl_mvoprobit(par = object$par, cov = object$cov, 
                              n_par = object$control_lnL$n_par, 
                              n_eq = object$control_lnL$n_eq, 
                              n_eq2 = object$control_lnL$n_eq2,
                              n_cuts_eq = object$control_lnL$n_cuts_eq,
                              n_regimes = object$control_lnL$n_regimes,
                              regimes_pair = object$other$regimes_pair,
                              coef = object$coef, coef_ind = object$ind$coef,
                              coef_var = object$coef_var, 
                              coef_var_ind = object$ind$coef_var,
                              cuts = object$cuts, cuts_ind = object$ind$cuts,
                              sigma = object$sigma, 
                              sigma_ind = object$ind$sigma,
                              sigma_vec_ind = object$other$sigma_vec,
                              coef2 = object$coef2, 
                              coef2_ind = object$ind$coef2,
                              var2 = object$var2, var2_ind = object$ind$var2,
                              cov2 = object$cov2, cov2_ind = object$ind$cov2,
                              sigma2 = object$sigma2, 
                              sigma2_ind = object$ind$sigma2,
                              marginal_par = object$marginal_par, 
                              marginal_par_ind = object$ind$marginal_par, 
                              marginal_par_n = object$control_lnL$marginal_par_n,
                              is_marginal = object$other$is_marginal, 
                              is_het = object$other$is_het, 
                              z_names = object$other$z_names, 
                              y_names = object$other$y_names,
                              estimator = object$estimator, 
                              coef_lambda_row = object$other$coef_lambda_row, 
                              cov_2step = object$cov_2step)
    object$tbl <- tbl_list$tbl
    object$se <- tbl_list$se
    object$p_value <- tbl_list$p_value
  }
  
  class(object) <- "summary.mvoprobit"
  
  return(object)
}

#' Print summary for an Object of Class mvoprobit
#' @description Prints summary for an object of class 'mvoprobit'.
#' @param x object of class 'mvoprobit'
#' @param ... further arguments (currently ignored)
#' @return The function returns \code{x}.
print.summary.mvoprobit <- function(x, ...)
{
  show_ind <- x$show_ind
  
  class(x) <- "mvoprobit"
  tbl_coef <- x$tbl$coef
  tbl_coef_var <- x$tbl$coef_var
  tbl_sigma <- x$tbl$sigma
  tbl_cuts <- x$tbl$cuts
  tbl_coef2 <- x$tbl$coef2
  tbl_var2 <- x$tbl$var2
  tbl_cov2 <- x$tbl$cov2
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
  var2_ind <- x$ind$var2
  cov2_ind <- x$ind$cov2
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
             " model \n"))
  if ((n_eq > 1) | (n_eq2 > 0))
  {
    cat(paste0("There are ", n_eq, " ordered ", 
               ifelse(n_eq2 > 0, paste0("and ", n_eq2, " continuous "), ""),
               "equations \n"))
  }
  cat("--- \n")
  if (!is.na(x$logLik))
  {
    cat(paste0("Log-likelihood = ", round(x$logLik, 4), " \n"))
  }
  if (!is.na(AIC(x)))
  {
    cat(paste0("AIC = ", round(AIC(x), 4), "\n"))
  }
  cat(paste0("Observations = ", n_obs, "\n"))
  
  for (i in 1:n_eq)
  {
    cat("--- \n")
    if ((n_eq >= 2) | (n_eq2 >= 1))
    {
      cat(paste0(x$formula[[i]][[2]], " equation \n"))
      cat("- \n");
    }
    if (is_marginal)
    {
      cat(paste0("Distribution: ", marginal_names[i], "\n"))
      cat("- \n")
    }
    cat("Coefficients: \n")
    printCoefmat(remove_column(tbl_coef[[i]], "ind", !show_ind), 
                 has.Pvalue = TRUE, signif.legend = FALSE, P.values = TRUE)
    if (is_het[i])
    {
      cat("- \n")
      cat("Coefficients of variance equation: \n")
      printCoefmat(remove_column(tbl_coef_var[[i]], "ind", !show_ind),
                   has.Pvalue = TRUE, signif.legend = FALSE, P.values = TRUE)
    }
    cat("- \n")
    cat("Cuts: \n")
    printCoefmat(remove_column(tbl_cuts[[i]], "ind", !show_ind),
                 has.Pvalue = TRUE, signif.legend = FALSE, P.values = TRUE)
    if (length(marginal_par[[i]]) > 0)
    {
      cat("- \n")
      cat(paste0("Parameters of ", marginal_names[i], " distribution: \n"))
      printCoefmat(remove_column(tbl_marginal_par[[i]], "ind", !show_ind),
                   has.Pvalue = TRUE, signif.legend = FALSE, P.values = TRUE)
    }
  }
  
  if (n_eq2 > 0)
  {
    for (i in 1:n_eq2)
    {
      cat("--- \n")
      cat(paste0(all.vars(x$formula2[[i]])[1], " equation \n"))
      for (j in 1:n_regimes[i])
      {
        cat("-- \n")
        if (n_regimes[i] > 1)
        {
          cat(paste0("Regime ", j - 1, ": \n"))
          cat("- \n")
        }
        if (x$estimator == "2step")
        {
          r_squared <- round(summary(x$twostep[[j]])$r.squared, 4)
          r_adj_squared <- round(summary(x$twostep[[j]])$adj.r.squared, 4)
          loocv_val <- round(loocv(x$twostep[[j]]), 4)
          cat(paste0("R-squared: ", r_squared, "\n"))
          cat(paste0("Adjusted R-squared: ", r_adj_squared, "\n"))
          cat(paste0("Leave-one-out cross-validation RMSE: ", loocv_val, "\n"))
          cat("- \n")
        }
        cat("Coefficients: \n")
        printCoefmat(remove_column(tbl_coef2[[i]][[j]], "ind", !show_ind),
                     has.Pvalue = TRUE, signif.legend = FALSE, P.values = TRUE)
        cat("- \n")
        if (x$estimator == "2step")
        {
          if (x$cov_type == "parametric")
          {
            cat(paste0("Variance stimate: ", tbl_var2[[i]][j, 1], "\n"))
            cat("standard errors and p-values are unavailable \n")
            cat("for variance when two-step estimator is used \n")
            cat("- \n")
            cat(paste0("Covariances with ordered equations: \n"))
            printCoefmat(remove_column(tbl_cov2[[i]][[j]], "ind", !show_ind),
                         has.Pvalue = TRUE, signif.legend = FALSE, 
                         P.values = TRUE)
          }
          else
          {
            if (length(tbl_lambda[[j]]) > 1)
            {
              cat(paste0("Selectivity correction terms: \n"))
              printCoefmat(tbl_lambda[[j]], has.Pvalue = TRUE, 
                           tst.ind = 3, signif.legend = FALSE, P.values = TRUE)
            }
          }
        }
        else
        {
          cat(paste0("Variance: \n"))
          printCoefmat(remove_column(tbl_var2[[i]][j, , drop = FALSE], "ind", 
                                     !show_ind),
                       has.Pvalue = TRUE, signif.legend = FALSE, 
                       P.values = TRUE)
          cat("- \n")
          cat(paste0("Covariances with ordered equations: \n"))
          printCoefmat(remove_column(tbl_cov2[[i]][[j]], "ind", !show_ind),
                       has.Pvalue = TRUE, signif.legend = FALSE, 
                       P.values = TRUE)
        }
      }
    }
  }
  
  if (n_eq >= 2)
  {
    cat("--- \n")
    cat("Correlations between ordered equations: \n")
    printCoefmat(remove_column(tbl_sigma, "ind", !show_ind), 
                 has.Pvalue = TRUE, signif.legend = FALSE, P.values = TRUE)
  }
  
  if (n_eq2 >= 2)
  {
    cat("--- \n")
    cat("Covariances between continuous equations: \n")
    counter <- 1
    for (i in 2:n_eq2)
    {
      for (j in 1:(i - 1))
      {
        cat("- \n")
        cat(paste0("Between ", x$formula2[[i]][[2]], 
                   " and ", x$formula2[[j]][[2]], "\n"))
        printCoefmat(remove_column(tbl_sigma2[[counter]], "ind", !show_ind),
                     has.Pvalue = TRUE, signif.legend = FALSE, P.values = TRUE)
        counter <- counter + 1
      }
    }
  }
  
  if (getOption("show.signif.stars"))
  {
    cat("--- \n")
    cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 \n")
  }
  
  return(x)
}

#' Print for an Object of Class mvoprobit
#' @description Prints information on the object of class 'mvoprobit'.
#' @param x object of class 'mvoprobit'
#' @param ... further arguments (currently ignored)
#' @return The function returns \code{NULL}.
print.mvoprobit <- function(x, ...)
{
  if (length(list(...)) > 0)
  {
    warning("Additional arguments passed through ... are ignored.")   
  }
  
  cat("This object of class 'mvoprobit' contains the following elements: \n")
  cat("--- \n")
  print(summary.default(x))
  cat("--- \n")
  cat(paste0("Structure of the model i.e. correspondence between possible \n",
             "values of ordered equations and regimes of continuous ",
             "equations: \n"))
  cat("--- \n")
  groups_all <- x$groups
  colnames
  if (x$other$is2)
  {
    groups_all <- cbind(groups_all, x$groups2)
  }
  groups_all <- cbind(groups_all, obs = x$control_lnL$n_obs_g)
  rownames(groups_all) <- rep("", nrow(groups_all))
  print(groups_all)
  cat("obs - the number of observations with corresponding values. \n")
  cat("--- \n")
  cat("Use 'summary()' to get more detailed information on the output. \n")
  return(NULL)
}

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

#' Summary for an Object of Class mnprobit
#' @description Provides summary for an object of class 'mnprobit'.
#' @param object object of class 'mnprobit'
#' @param ... further arguments (currently ignored)
#' @param vcov positively defined numeric matrix representing
#' asymptotic variance-covariance matrix of the estimator to be
#' used for calculation of standard errors and p-values. It may also be a 
#' character. Then \code{\link[switchSelection]{vcov.mnprobit}} function
#' will be used which input argument \code{type} will be set to \code{vcov}.
#' If \code{estimator = "2step"} then \code{vcov} should be an estimate of the 
#' asymptotic covariance matrix of the first step estimator.
#' @param show_ind logical; if \code{TRUE} then indexes of parameters will be
#' shown. Particularly, these indexes may be used in \code{ind} element of
#' \code{regularization} parameter of \code{\link[switchSelection]{mvoprobit}}.
#' @details If \code{vcov} is \code{NULL} then this function just changes the 
#' class of the 'mnprobit' object to 'summary.mnprobit'. Otherwise it 
#' additionally changes \code{object$cov} to \code{vcov} and use it to
#' recalculate \code{object$se}, \code{object$p_value} and \code{object$tbl} 
#' values. It also adds the value of \code{ind} argument to the object.
#' @return Returns an object of class 'summary.mnprobit'.
summary.mnprobit <- function(object, ..., vcov = NULL, show_ind = FALSE) 
{
  if (length(list(...)) > 0)
  {
    warning("Additional arguments passed through ... are ignored.")   
  }
  
  # Add information on the indexes
  object$show_ind <- show_ind
  
  if (object$estimator == "2step")
  {
    if (!is.null(vcov))
    {
      warning("Argument 'vcov' is ignored for two-step estimator.")
    }
    vcov <- NULL
  }
  
  if (!is.null(vcov))
  {
    if (is.character(vcov))
    {
      vcov <- vcov.mnprobit(object, type = vcov)
    }
    if (any(dim(vcov) != dim(object$cov)))
    {
      stop("Incorrect size of 'vcov' matrix.")
    }
    object$cov <- vcov
    tbl_list <- tbl_mnprobit(par = object$par, cov = object$cov, 
                             n_par = object$control_lnL$n_par,
                             n_alt = object$n_alt,
                             n_coef = object$control_lnL$n_coef, 
                             n_coef2 = object$control_lnL$n_coef2,
                             n_regimes = object$control_lnL$n_regimes,
                             coef_ind_alt = object$control_lnL$coef_ind_alt + 1, 
                             sigma_ind = object$control_lnL$sigma_ind + 1, 
                             sigma_ind_mat = object$control_lnL$sigma_ind_mat + 1,
                             coef2_ind_regime = object$control_lnL$coef2_ind_regime + 1, 
                             var2_ind_regime = object$control_lnL$var2_ind_regime + 1, 
                             cov2_ind_regime = object$control_lnL$cov2_ind_regime + 1,
                             alt_names = object$alt_names, 
                             coef_lambda = object$coef_lambda, 
                             is2 = object$other$is2, 
                             degrees = object$degrees,
                             estimator = object$estimator,
                             cov_2step = object$cov_2step,
                             coef = object$coef, 
                             sigma_vec = object$other$sigma_vec,
                             coef2 = object$coef2,
                             var2 = object$var2,
                             cov2 = object$cov2,
                             coef_lambda_row = object$other$coef_lambda_row,
                             regimes_names = object$regimes_names,
                             colnames_W = object$other$z_names, 
                             colnames_X = object$other$y_names)
    object$tbl <- tbl_list$tbl
    object$se <- tbl_list$se
    object$p_value <- tbl_list$p_value
  }
  
  class(object) <- "summary.mnprobit";
  
  return(object)
}

#' Print summary for an Object of Class mnprobit
#' @description Prints summary for an object of class 'mnprobit'.
#' @param x object of class 'mnprobit'
#' @param ... further arguments (currently ignored)
#' @return The function returns \code{x}.
print.summary.mnprobit <- function(x, ...)
{
  show_ind <- x$show_ind
  
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
  
  alt_names <- x$other$alt_names
  regimes_names <- x$other$regimes_names
  
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
    cat("Model with multinomial probit selection mechanism \n")
  }
  else
  {
    cat("Multinomial probit model \n")
  }
  cat(paste0("There are ", n_alt, " alternatives",
             ifelse(is2, paste0(" and ", n_regimes, " regimes"), ". \n")))
  cat("--- \n")
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
  
  # coefficients
  for (i in 1:(n_alt - 1))
  {
    cat("--- \n")
    cat(paste0("Alternative ", i, "\n"))
    cat("- \n");
    cat("Coefficients: \n")
    printCoefmat(remove_column(tbl_coef[[i]], "ind", !show_ind),
                 has.Pvalue = TRUE, signif.legend = FALSE, P.values = TRUE)
    cat("- \n");
    cat(paste0("Observations = ", length(x$control_lnL$ind_alt[[i]]), "\n"))
  }
  
  # covariances
  cat("--- \n")
  cat("Covariances between alternatives: \n")
  printCoefmat(remove_column(tbl_sigma, "ind", !show_ind), 
               has.Pvalue = TRUE, signif.legend = FALSE, P.values = TRUE)
  
  if (is2)
  {
    for (i in 1:n_regimes)
    {
      cat("--- \n")
      if (n_regimes > 1)
      {
        cat(paste0("Regime ", i - 1, "\n"))
      }
      else
      {
        cat(paste0(all.vars(x$formula2)[1], "\n"))
      }
      cat("- \n");
      if (x$estimator == "2step")
      {
        r_squared <- round(summary(x$twostep[[i]])$r.squared, 4)
        r_adj_squared <- round(summary(x$twostep[[i]])$adj.r.squared, 4)
        loocv_val <- round(loocv(x$twostep[[i]]), 4)
        cat(paste0("R-squared: ", r_squared, "\n"))
        cat(paste0("Adjusted R-squared: ", r_adj_squared, "\n"))
        cat(paste0("Leave-one-out cross-validation RMSE: ", loocv_val, "\n"))
        cat("- \n")
      }
      cat("Coefficients: \n")
      printCoefmat(remove_column(tbl_coef2[[i]], "ind", !show_ind), 
                   has.Pvalue = TRUE, signif.legend = FALSE, P.values = TRUE)
      if (estimator == "ml")
      {
        cat("- \n");
        cat("Variance: \n")
        printCoefmat(remove_column(tbl_var2[[i]], "ind", !show_ind), 
                     has.Pvalue = TRUE, signif.legend = FALSE, P.values = TRUE)
        cat("- \n");
        cat("Covariances with alternatives: \n")
        printCoefmat(remove_column(tbl_cov2[[i]], "ind", !show_ind), 
                     has.Pvalue = TRUE, signif.legend = FALSE, P.values = TRUE)
      }
      else
      {
        cat("- \n");
        cat(paste0("Selectivity correction terms: \n"))
        printCoefmat(tbl_lambda[[i]], has.Pvalue = TRUE, 
                     tst.ind = 3, signif.legend = FALSE, P.values = TRUE)
      }
      cat("- \n");
      cat(paste0("Observations = ", length(x$control_lnL$ind_regime[[i]]), 
                 "\n"))
    }
  }
  
  if (getOption("show.signif.stars"))
  {
    cat("--- \n")
    cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 \n")
  }
  
  return(x)
}

#' Print for an Object of Class mnprobit
#' @description Prints information on the object of class 'mnprobit'.
#' @param x object of class 'mnprobit'
#' @param ... further arguments (currently ignored)
#' @return The function returns input argument \code{NULL}.
print.mnprobit <- function(x, ...)
{
  if (length(list(...)) > 0)
  {
    warning("Additional arguments passed through ... are ignored.")   
  }
  
  cat("This object of class 'mnprobit' contains the following elements: \n")
  cat("--- \n")
  print(summary.default(x))
  cat("--- \n")
  cat("Use 'summary()' to get more detailed information on the output.")
  return(NULL)
}