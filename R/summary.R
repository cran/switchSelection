#' Summary for an Object of Class msel
#' @description Provides summary for an object of class 'msel'.
#' @param object object of class 'msel'
#' @param ... further arguments (currently ignored)
#' @param vcov positively defined numeric matrix representing
#' asymptotic variance-covariance matrix of the estimator to be
#' used for calculation of standard errors and p-values. It may also be a 
#' character. Then \code{\link[switchSelection]{vcov.msel}} function
#' will be used which input argument \code{type} will be set to \code{vcov}.
#' If \code{estimator = "2step"} then \code{vcov} should be an estimate of the 
#' asymptotic covariance matrix of the first step estimator.
#' @param show_ind logical; if \code{TRUE} then indexes of parameters will be
#' shown. Particularly, these indexes may be used in \code{ind} element of
#' \code{regularization} parameter of \code{\link[switchSelection]{msel}}.
#' @details If \code{vcov} is \code{NULL} then this function just changes the 
#' class of the 'msel' object to 'summary.msel'. Otherwise it 
#' additionally changes \code{object$cov} to \code{vcov} and use it to
#' recalculate \code{object$se}, \code{object$p_value} and \code{object$tbl} 
#' values. It also adds the value of \code{ind} argument to the object.
#' @return Returns an object of class 'summary.msel'.
summary.msel <- function(object, ..., vcov = NULL, show_ind = FALSE) 
{
  if (length(list(...)) > 0)
  {
    warning("Additional arguments passed through ... are ignored.")   
  }
  
  # Add information on the indexes
  object$show_ind <- show_ind
  
  if (!is.null(vcov))
  {
    if (is.character(vcov))
    {
      vcov <- vcov.msel(object, type = vcov)
    }
    if (any(dim(vcov) != dim(object$cov)))
    {
      stop("Incorrect size of 'vcov' matrix.")
    }
    tbl_list       <- tbl_msel(object = object, vcov = vcov)
    object$tbl     <- tbl_list$tbl
    object$se      <- tbl_list$se
    object$p_value <- tbl_list$p_value
  }
  
  class(object) <- "summary.msel"
  
  return(object)
}

#' Print summary for an Object of Class msel
#' @description Prints summary for an object of class 'msel'.
#' @param x object of class 'msel'
#' @param ... further arguments (currently ignored)
#' @return The function returns \code{x}.
print.summary.msel <- function(x, ...)
{
  # Assign the class
  class(x) <- "msel"
  
  # Get some data
  show_ind         <- x$show_ind
  tbl_coef         <- x$tbl$coef
  tbl_coef_var     <- x$tbl$coef_var
  tbl_sigma        <- x$tbl$sigma
  tbl_cuts         <- x$tbl$cuts
  tbl_coef2        <- x$tbl$coef2
  tbl_lambda       <- x$tbl$lambda
  tbl_var2         <- x$tbl$var2
  tbl_cov2         <- x$tbl$cov2
  tbl_sigma2       <- x$tbl$sigma2
  tbl_sigma3       <- x$tbl$sigma3
  tbl_marginal_par <- x$tbl$marginal_par
  tbl_coef3        <- x$tbl$coef3
  n_eq             <- x$other$n_eq
  n_eq2            <- x$other$n_eq2
  n_eq3            <- x$other$n_eq3
  coef_ind         <- x$ind$coef
  coef_var_ind     <- x$ind$coef_var
  cuts_ind         <- x$ind$cuts
  coef2_ind        <- x$ind$coef2
  sigma_ind        <- x$ind$sigma
  sigma2_ind       <- x$ind$sigma2
  sigma3_ind       <- x$ind$sigma3
  var2_ind         <- x$ind$var2
  cov2_ind         <- x$ind$cov2
  marginal_par_ind <- x$ind$marginal_par
  n_obs            <- x$other$n_obs
  is_het           <- x$other$is_het
  n_regimes        <- x$other$n_regimes
  marginal         <- x$marginal
  marginal_par     <- x$marginal_par
  marginal_names   <- names(marginal)
  is_marginal      <- length(marginal) > 0
  is1              <- x$other$is1
  is2              <- x$other$is2
  is3              <- x$other$is3
  type3            <- x$type3
  estimator        <- x$estimator
  coef_lambda_ind  <- x$other$coef_lambda_ind
  z_mn_names       <- x$other$z_mn_names
  
  # Provide the name to the model
  cat(paste0(ifelse(estimator == "ml", "Maximum-likelihood ", "Two-step "),
             "estimator is used\n"))
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
  
  if (is1)
  {
    cat("--- \n")
    cat("Ordinal equations \n")
    for (i in 1:n_eq)
    {
      cat("-- \n")
      cat(paste0(x$formula[[i]][[2]], " equation \n"))
      cat("- \n");
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
        cat("Coefficients of the variance equation: \n")
        printCoefmat(remove_column(tbl_coef_var[[i]], "ind", !show_ind),
                     has.Pvalue = TRUE, signif.legend = FALSE, P.values = TRUE)
      }
      cat("- \n")
      cat("Cuts (thresholds): \n")
      printCoefmat(remove_column(tbl_cuts[[i]], "ind", !show_ind),
                   has.Pvalue = TRUE, signif.legend = FALSE, P.values = TRUE)
      if (length(marginal_par[[i]]) > 0)
      {
        cat("- \n")
        cat(paste0("Parameters of the ", marginal_names[i], 
                   " distribution: \n"))
        printCoefmat(remove_column(tbl_marginal_par[[i]], "ind", !show_ind),
                     has.Pvalue = TRUE, signif.legend = FALSE, P.values = TRUE)
      }
    }
    if (n_eq >= 2)
    {
      cat("--- \n")
      cat("Correlations between the ordinal equations: \n")
      printCoefmat(remove_column(tbl_sigma, "ind", !show_ind), 
                   has.Pvalue = TRUE, signif.legend = FALSE, P.values = TRUE)
    }
  }
  
  if (is3)
  {
    cat("--- \n")
    cat(paste0(z_mn_names[1], " multinomial ", type3, " equation \n"))
    
    # Coefficients
    for (i in 1:(n_eq3 - 1))
    {
      cat("-- \n")
      cat(paste0("Alternative ", i - 1, "\n"))
      cat("- \n")
      cat(paste0("Coefficients: \n"))
      printCoefmat(remove_column(tbl_coef3[[i]], "ind", !show_ind),
                   has.Pvalue = TRUE, signif.legend = FALSE, 
                   P.values = TRUE)
    }
    
    # Covariances between the alternatives
    if (type3 == "probit")
    {
      cat("- \n")
      cat("Covariances between the alternatives: \n")
      printCoefmat(remove_column(tbl_sigma3, "ind", !show_ind), 
                   has.Pvalue = TRUE, signif.legend = FALSE, P.values = TRUE)
    }
  }
  
  if (is2)
  {
    cat("--- \n")
    cat("Continuous equations \n")
    for (i in 1:n_eq2)
    {
      cat("-- \n")
      cat(paste0(all.vars(x$formula2[[i]])[1], " equation \n"))
      for (j in 1:n_regimes[i])
      {
        cat("-- \n")
        if (n_regimes[i] > 1)
        {
          cat(paste0("Regime ", j - 1, ": \n"))
          cat("- \n")
        }
        if (estimator == "2step")
        {
          r_squared <- round(summary(x$twostep[[i]][[j]])$r.squared, 4)
          r_adj_squared <- round(summary(x$twostep[[i]][[j]])$adj.r.squared, 4)
          loocv_val <- round(loocv(x$twostep[[i]][[j]]), 4)
          cat(paste0("R-squared: ", r_squared, "\n"))
          cat(paste0("Adjusted R-squared: ", r_adj_squared, "\n"))
          cat(paste0("Leave-one-out cross-validation RMSE: ", loocv_val, "\n"))
          cat("- \n")
        }
        cat("Coefficients: \n")
        printCoefmat(remove_column(tbl_coef2[[i]][[j]], "ind", !show_ind),
                     has.Pvalue = TRUE, signif.legend = FALSE, P.values = TRUE)
        if ((estimator == "2step"))
        {
          if (length(coef_lambda_ind[[i]]) > 0)
          {
            cat("- \n")
            cat("Selectivity terms: \n")
            printCoefmat(remove_column(tbl_lambda[[i]][[j]], "ind", !show_ind),
                         has.Pvalue = TRUE, signif.legend = FALSE, 
                         P.values = TRUE)
          }
        }
        if (estimator != "2step")
        {
          cat("- \n")
          cat(paste0("Variance: \n"))
          printCoefmat(remove_column(tbl_var2[[i]][j, , drop = FALSE], "ind", 
                                     !show_ind),
                       has.Pvalue = TRUE, signif.legend = FALSE, 
                       P.values = TRUE)
          if (is1)
          {
            cat("- \n")
            cat(paste0("Covariances with the ordinal equations: \n"))
            printCoefmat(remove_column(tbl_cov2[[i]][[j]], "ind", !show_ind),
                         has.Pvalue = TRUE, signif.legend = FALSE, 
                         P.values = TRUE)
          }
        }
      }
    }
  }
  
  if ((n_eq2 >= 2) & (estimator == "ml"))
  {
    cat("--- \n")
    cat("Covariances between the continuous equations: \n")
    counter <- 1
    for (i in 2:n_eq2)
    {
      for (j in 1:(i - 1))
      {
        cat("- \n")
        cat(paste0("Between ", x$formula2[[i]][[2]], 
                   " and ", x$formula2[[j]][[2]], "\n"))
        if (length(x$other$regimes_pair[[counter]]) > 0)
        {
          printCoefmat(remove_column(tbl_sigma2[[counter]], "ind", !show_ind),
                       has.Pvalue = TRUE, signif.legend = FALSE, 
                       P.values = TRUE)
        }
        else
        {
          cat("Unidentified\n")
        }
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

#' Print for an Object of Class msel
#' @description Prints information on the object of class 'msel'.
#' @param x object of class 'msel'
#' @param ... further arguments (currently ignored)
#' @return The function returns \code{NULL}.
print.msel <- function(x, ...)
{
  if (length(list(...)) > 0)
  {
    warning("Additional arguments passed through ... are ignored.")   
  }
  
  cat("This object of class 'msel' contains the following elements: \n")
  cat("--- \n")
  print(summary.default(x))
  cat("--- \n")
  print(struct_msel(x))
  cat("--- \n")
  cat("Use 'summary()' to get more detailed information on the output. \n")
  return(NULL)
}

#' Structure of the Object of Class msel
#' @description Prints information on the structure of the model.
#' @param x object of class 'msel'
#' @return The function returns a numeric matrix which columns are 
#' \code{groups}, \code{groups2}, \code{groups3} correspondingly. It also has
#' additional (last) column with the number of observations associated with the
#' corresponding combinations of the groups.
struct_msel <- function(x)
{
  # Prepare the variable to store the output
  out <- NULL
  
  # Combine the groups
  if (x$other$is1)
  {
    out <- cbind(out, x$groups)
  }
  if (x$other$is2)
  {
    out <- cbind(out, x$groups2)
  }
  if (x$other$is3)
  {
    out                             <- cbind(out, x$groups3)
    colnames(out)[ncol(out)] <- "alt"
  }
  
  # Add the number of observations
  out <- cbind(out, obs = x$control_lnL$n_obs_g)
  rownames(out) <- rep("", nrow(out))
  
  # Assign the class
  class(out) <- "struct_msel"
  
  # Return the result
  return(out)
}

#' Print for an Object of Class struct_msel
#' @description Prints information on the object of class 'struct_msel'.
#' @param x object of class 'struct_msel'
#' @param ... further arguments (currently ignored)
#' @return The function returns \code{NULL}.
print.struct_msel <- function(x, ...)
{
  if (length(list(...)) > 0)
  {
    warning("Additional arguments passed through ... are ignored.")   
  }
  
  cat(paste0("Structure of the model i.e. correspondence between the possible ",
             "values\n",
             "of the ordinal equations and the regimes of the ",
             "continuous equations: \n"))
  cat("--- \n")
  attr(x,"class") <- NULL
  print.default(x)
  cat("obs - the number of observations with corresponding values. \n")
  
  return(NULL)
}