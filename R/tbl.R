# Create output tables msel
tbl_msel <- function(object, par = NULL, vcov = NULL)
{
  # Get some variables
  if (!is.null(par))
  {
    object <- par_msel(object = object, par = par, type = "object")
  }
  if (!is.null(vcov))
  {
    object$cov <- vcov
  }
  par              <- object$par
  cov              <- object$cov
  n_par            <- object$other$n_par
  estimator        <- object$estimator
  is1              <- object$other$is1
  is2              <- object$other$is2
  is3              <- object$other$is3
  is_het           <- object$other$is_het
  n_eq             <- object$other$n_eq
  n_eq2            <- object$other$n_eq2
  n_eq3            <- object$other$n_eq3
  n_cuts_eq        <- object$other$n_cuts_eq
  n_regimes        <- object$other$n_regimes
  regimes_pair     <- object$other$regimes_pair
  coef             <- object$coef
  coef_var         <- object$coef_var
  cuts             <- object$cuts
  sigma            <- object$sigma
  coef2            <- object$coef2
  var2             <- object$var2
  cov2             <- object$cov2
  sigma2           <- object$sigma2
  sigma3           <- object$sigma3
  marginal_par     <- object$marginal_par
  coef3            <- object$coef3
  coef_ind         <- object$ind$coef
  coef_var_ind     <- object$ind$coef_var
  cuts_ind         <- object$ind$cuts 
  sigma_ind        <- object$ind$sigma
  sigma_vec_ind    <- object$other$sigma_vec_ind
  coef2_ind        <- object$ind$coef2
  var2_ind         <- object$ind$var2
  cov2_ind         <- object$ind$cov2
  sigma2_ind       <- object$ind$sigma2
  sigma3_ind       <- object$ind$sigma3
  sigma3_ind_mat   <- object$ind$sigma3_mat
  marginal_par_ind <- object$ind$marginal_par
  marginal_par_n   <- object$other$marginal_par_n
  coef3_ind        <- object$ind$coef3
  is_marginal      <- object$other$is_marginal
  z_names          <- object$other$z_names
  y_names          <- object$other$y_names
  type3            <- object$type3
  coef_lambda_ind  <- object$other$coef_lambda_ind
  
  # Calculate standard errors
  se <- sqrt(diag(cov))
  
  # Calculate p-values
  z_value <- par / se
  p_value <- rep(NA, n_par)
  for (i in 1:n_par)
  {
    p_value[i] <- 2 * min(pnorm(z_value[i]), 1 - pnorm(z_value[i]))
  }
  
  # Construct output table
  tbl_coef         <- vector(mode = "list", length = n_eq)
  tbl_coef_var     <- vector(mode = "list", length = n_eq)
  tbl_cuts         <- vector(mode = "list", length = n_eq)
  tbl_sigma        <- vector(mode = "list", length = 1)
  tbl_coef2        <- vector(mode = "list", length = n_eq2)
  tbl_lambda       <- vector(mode = "list", length = n_eq2)
  tbl_var2         <- vector(mode = "list", length = n_eq2)
  tbl_cov2         <- vector(mode = "list", length = n_eq2)
  tbl_sigma2       <- vector(mode = "list", length = n_eq2 * (n_eq2 - 1) / 2)
  tbl_marginal_par <- NULL
  tbl_coef3        <- NULL
  tbl_sigma3       <- NULL
  
  # Ordered equations
  if (is1)
  {
    for (i in 1:n_eq)
    {
      # coefficients of the mean equation
      tbl_coef[[i]] <- as.matrix(cbind(Estimate = coef[[i]],
                                       Std_Error = se[coef_ind[[i]]],
                                       z_value = z_value[coef_ind[[i]]],
                                       ind = coef_ind[[i]],
                                       p_value = p_value[coef_ind[[i]]]))
      rownames(tbl_coef[[i]]) <- names(coef[[i]])
      names(tbl_coef) <- z_names
      
      # coefficients of the variance equation
      if (is_het[i])
      {
        tbl_coef_var[[i]] <- as.matrix(cbind(Estimate = coef_var[[i]],
                                             Std_Error = se[coef_var_ind[[i]]],
                                             z_value = z_value[coef_var_ind[[i]]],
                                             ind = coef_var_ind[[i]],
                                             p_value = p_value[coef_var_ind[[i]]]))
        rownames(tbl_coef_var[[i]]) <- names(coef_var[[i]])
      }
      names(tbl_coef_var) <- z_names
      
      # cuts
      tbl_cuts[[i]] <- as.matrix(cbind(Estimate = cuts[[i]],
                                       Std_Error = se[cuts_ind[[i]]],
                                       z_value = z_value[cuts_ind[[i]]],
                                       ind = cuts_ind[[i]],
                                       p_value = p_value[cuts_ind[[i]]]))
      rownames(tbl_cuts[[i]]) <- paste0("cut", 1:n_cuts_eq[i])
    }
  }
  
  # covariance matrix of ordered equations
  tbl_sigma <- NULL
  if (is1)
  {
    if (n_eq >= 2)
    {
      sigma_vec_names <- paste0("cov(", 
                                z_names[sigma_vec_ind[, 1, drop = FALSE]], ", ",
                                z_names[sigma_vec_ind[, 2, drop = FALSE]],
                                ")")
      tbl_sigma <- as.matrix(cbind(Estimate = par[sigma_ind],
                                   Std_Error = se[sigma_ind],
                                   z_value = z_value[sigma_ind],
                                   ind = sigma_ind,
                                   p_value = p_value[sigma_ind]))
      rownames(tbl_sigma) <- sigma_vec_names
    }
  }
  
  # Continuous equations
  if (is2)
  {
    for (i in 1:n_eq2)
    {
      tbl_coef2[[i]]  <- vector(mode = "list", length = n_regimes[i])
      tbl_lambda[[i]] <- vector(mode = "list", length = n_regimes[i])
      tbl_cov2[[i]]   <- vector(mode = "list", length = n_regimes[i])
      for (j in 1:n_regimes[i])
      {
        # coefficients
        tbl_coef2[[i]][[j]] <- as.matrix(cbind(Estimate = coef2[[i]][j, ],
                                               Std_Error = se[coef2_ind[[i]][j, ]],
                                               z_value = z_value[coef2_ind[[i]][j, ]],
                                               ind = coef2_ind[[i]][j, ],
                                               p_value = p_value[coef2_ind[[i]][j, ]]))
        rownames(tbl_coef2[[i]][[j]]) <- names(coef2[[i]][j, ])
        names(tbl_coef2[[i]]) <- paste0("regime ", 0:(n_regimes[i] - 1))
        if ((estimator == "2step"))
        {
          if (length(coef_lambda_ind[[i]]) > 0)
          {
            tbl_lambda[[i]][[j]] <- tbl_coef2[[i]][[j]][coef_lambda_ind[[i]], ,
                                                        drop = FALSE]
            tbl_coef2[[i]][[j]]  <- tbl_coef2[[i]][[j]][-coef_lambda_ind[[i]], ,
                                                        drop = FALSE]
          }
        }
        
        # covariances with ordered equations
        if (is1 & (estimator == "ml"))
        {
          tbl_cov2[[i]][[j]] <- as.matrix(cbind(Estimate = cov2[[i]][j, ],
                                                 Std_Error = se[cov2_ind[[i]][j, ]],
                                                 z_value = z_value[cov2_ind[[i]][j, ]],
                                                 ind = cov2_ind[[i]][j, ],
                                                 p_value = p_value[cov2_ind[[i]][j, ]]))
          rownames(tbl_cov2[[i]][[j]]) <- z_names
          names(tbl_coef2[[i]]) <- paste0("regime ", 0:(n_regimes[i] - 1))
        }
      }
      if (estimator != "2step")
      {
        tbl_var2[[i]] <- as.matrix(cbind(Estimate = var2[[i]],
                                          Std_Error = se[var2_ind[[i]]],
                                          z_value = z_value[var2_ind[[i]]],
                                          ind = var2_ind[[i]],
                                          p_value = p_value[var2_ind[[i]]]))
        rownames(tbl_var2[[i]]) <- rep("var", n_regimes[i])
      }
    }
    names(tbl_coef2) <- paste0("Equation ", 1:n_eq2)
  }
  
  # Covariances between continuous equations
  if ((n_eq2 >= 2) & (estimator == "ml"))
  {
    counter <- 1
    for (i in 2:n_eq2)
    {
      for (j in 1:(i - 1))
      {
        tbl_sigma2[[counter]] <- as.matrix(cbind(Estimate = sigma2[[counter]],
                                                 Std_Error = se[sigma2_ind[[counter]]],
                                                 z_value = z_value[sigma2_ind[[counter]]],
                                                 ind = sigma2_ind[[counter]],
                                                 p_value = p_value[sigma2_ind[[counter]]]))
        names(tbl_sigma2)[counter] <- paste0(y_names[i], " and ", y_names[j])
        rownames(tbl_sigma2[[counter]]) <- rep("", nrow(tbl_sigma2[[counter]]))
        if (length(regimes_pair[[counter]]) > 0)
        {
          for (t in 1:nrow(regimes_pair[[counter]]))
          {
            rownames(tbl_sigma2[[counter]])[t] <- paste0("(",
                                                         regimes_pair[[counter]][t, 1], 
                                                         ",",
                                                         regimes_pair[[counter]][t, 2],
                                                         ")")
          }
        }
        counter <- counter + 1
      }
    }
  }
  
  # Parameters of the marginal distributions
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
                ind = marginal_par_ind[[i]],
                p_value = p_value[marginal_par_ind[[i]]]))
        rownames(tbl_marginal_par[[i]]) <- paste0("par", 1:marginal_par_n[i])
        names(marginal_par[[i]]) <- paste0("par", 1:marginal_par_n[i])
      }
    }
  }
  
  # Coefficients of the multinomial equations
  if (is3)
  {
    tbl_coef3        <- vector(mode = "list", length = n_eq3 - 1)
    names(tbl_coef3) <- paste0("Alternative ", 1:(n_eq3 - 1))
    for (i in 1:(n_eq3 - 1))
    {
      tbl_coef3[[i]] <- as.matrix(cbind(Estimate = coef3[i, ],
                                        Std_Error = se[coef3_ind[i, ]],
                                        z_value = z_value[coef3_ind[i, ]],
                                        ind = coef3_ind[i, ],
                                        p_value = p_value[coef3_ind[i, ]]))
    }
  }
  
  # Covariances of the multinomial equations
  if (is3 & (n_eq3 > 2) & (type3 == "probit"))
  {
    sigma3_vec_names     <- paste0("cov(", sigma3_ind_mat[, 1, drop = FALSE], 
                                   ",",
                                   sigma3_ind_mat[, 2, drop = FALSE],
                                   ")")
    tbl_sigma3           <- as.matrix(cbind(Estimate  = par[sigma3_ind],
                                            Std_Error = se[sigma3_ind],
                                            z_value   = z_value[sigma3_ind],
                                            ind       = sigma3_ind,
                                            p_value   = p_value[sigma3_ind]))
    rownames(tbl_sigma3) <- sigma3_vec_names
  }
  
  # Aggregate the results into the table
  tbl <- list(coef         = tbl_coef,
              coef_var     = tbl_coef_var,
              cuts         = tbl_cuts,
              sigma        = tbl_sigma,
              coef2        = tbl_coef2,
              lambda       = tbl_lambda,
              var2         = tbl_var2,
              cov2         = tbl_cov2,
              sigma2       = tbl_sigma2,
              marginal_par = tbl_marginal_par,
              coef3        = tbl_coef3,
              sigma3       = tbl_sigma3)
  
  # Store the output
  out <- list(tbl = tbl, se = se, p_value = p_value)
  
  return (out)
}