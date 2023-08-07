# Create output tables mvoprobit
tbl_mvoprobit <- function(par, cov, 
                          n_par, n_eq, n_eq2, n_cuts_eq,
                          n_regimes, regimes_pair,
                          coef, coef_ind,
                          coef_var, coef_var_ind,
                          cuts, cuts_ind,
                          sigma, sigma_ind, sigma_vec_ind,
                          coef2, coef2_ind,
                          var2, var2_ind,
                          cov2, cov2_ind,
                          sigma2, sigma2_ind,
                          marginal_par, marginal_par_ind, marginal_par_n,
                          is_marginal, is_het, z_names, y_names,
                          estimator, coef_lambda_row, cov_2step)
{
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
  tbl_coef <- vector(mode = "list", length = n_eq)
  tbl_coef_var <- vector(mode = "list", length = n_eq)
  tbl_cuts <- vector(mode = "list", length = n_eq)
  tbl_sigma <- vector(mode = "list", length = 1)
  tbl_coef2 <- NULL
  tbl_var2 <- NULL
  tbl_cov2 <- NULL
  tbl_sigma2 <- NULL
  tbl_marginal_par <- NULL
  tbl_lambda <- NULL
  if (n_eq2 > 0)
  {
    tbl_coef2 <- vector(mode = "list", length = n_eq2)
    tbl_var2 <- vector(mode = "list", length = n_eq2)
    tbl_cov2 <- vector(mode = "list", length = n_eq2)
    tbl_sigma2 <- vector(mode = "list", length = n_eq2 * (n_eq2 - 1) / 2)
  }
  
  # Ordered equations
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
  
  # covariance matrix of ordered equations
  tbl_sigma <- NULL
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
  
  # Continuous equations
  if (n_eq2 > 0)
  {
    for (i in 1:n_eq2)
    {
      tbl_coef2[[i]] <- vector(mode = "list", length = n_regimes[i])
      tbl_cov2[[i]] <- vector(mode = "list", length = n_regimes[i])
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
        
        # covariances with ordered equations
        tbl_cov2[[i]][[j]] <- as.matrix(cbind(Estimate = cov2[[i]][j, ],
                                               Std_Error = se[cov2_ind[[i]][j, ]],
                                               z_value = z_value[cov2_ind[[i]][j, ]],
                                               ind = cov2_ind[[i]][j, ],
                                               p_value = p_value[cov2_ind[[i]][j, ]]))
        rownames(tbl_cov2[[i]][[j]]) <- z_names
        names(tbl_coef2[[i]]) <- paste0("regime ", 0:(n_regimes[i] - 1))
      }
      tbl_var2[[i]] <- as.matrix(cbind(Estimate = var2[[i]],
                                        Std_Error = se[var2_ind[[i]]],
                                        z_value = z_value[var2_ind[[i]]],
                                        ind = var2_ind[[i]],
                                        p_value = p_value[var2_ind[[i]]]))
      rownames(tbl_var2[[i]]) <- rep("var", n_regimes[i])
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
                                                 ind = sigma2_ind[[counter]],
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
  
  # Parameters of marginal distribution
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
  
  # Standard errors of coefficients of lambdas
  if (estimator == "2step")
  {
    tbl_lambda <- vector(mode = "list", length = n_regimes[1])
    for (i in 1:n_regimes[1])
    {
      n_coef2_tmp <- length(coef2_ind[[1]][i, ])
      estimate_lambda <- coef_lambda_row[[i]]
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
    names(tbl_lambda) <- paste0("regime ", 0:(n_regimes[1] - 1))
  }
  
  # Aggregate the results into the table
  tbl <- list(coef = tbl_coef,
              coef_var = tbl_coef_var,
              cuts = tbl_cuts,
              sigma = tbl_sigma,
              coef2 = tbl_coef2,
              var2 = tbl_var2,
              cov2 = tbl_cov2,
              sigma2 = tbl_sigma2,
              marginal_par = tbl_marginal_par,
              lambda = tbl_lambda)
  
  # Store the output
  out <- list(tbl = tbl, se = se, p_value = p_value)
  
  return (out)
}

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

tbl_mnprobit <- function(par, cov, n_par = n_par,
                         n_alt, n_coef, n_coef2, n_regimes,
                         coef_ind_alt, sigma_ind, sigma_ind_mat,
                         coef2_ind_regime, var2_ind_regime, cov2_ind_regime,
                         alt_names, coef_lambda, is2, degrees, estimator,
                         cov_2step,
                         coef, sigma_vec, coef2, var2, cov2,
                         coef_lambda_row, regimes_names,
                         colnames_W, colnames_X)
{
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
  tbl_coef <- vector(mode = "list", length = n_alt - 1)
  tbl_sigma <- matrix()
  tbl_coef2 <- vector(mode = "list", length = n_regimes)
  tbl_var2 <- vector(mode = "list", length = n_regimes)
  tbl_cov2 <- vector(mode = "list", length = n_regimes)
  tbl_lambda <- NULL
  
  # Coefficients of the multinomial equation
  for (i in 1:(n_alt - 1))
  {
    tbl_coef[[i]] <- as.matrix(cbind(Estimate = coef[, i],
                                     Std_Error = se[coef_ind_alt[, i]],
                                     z_value = z_value[coef_ind_alt[, i]],
                                     ind = coef_ind_alt[, i],
                                     p_value = p_value[coef_ind_alt[, i]]))
    rownames(tbl_coef[[i]]) <- colnames_W
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
                                 ind = sigma_ind,
                                 p_value = p_value[sigma_ind]))
    rownames(tbl_sigma) <- sigma_vec_names
  }
  
  # Coefficients of the continuous equations
  if (is2)
  {
    for (i in 1:n_regimes)
    {
      #coefficients
      tbl_coef2[[i]] <- as.matrix(cbind(Estimate = coef2[, i],
                                        Std_Error = se[coef2_ind_regime[, i]],
                                        z_value = z_value[coef2_ind_regime[, i]],
                                        ind = coef2_ind_regime[, i],
                                        p_value = p_value[coef2_ind_regime[, i]]))
      rownames(tbl_coef2[[i]]) <- colnames_X
      # variances
      tbl_var2[[i]] <- as.matrix(cbind(Estimate = var2[i],
                                       Std_Error = se[var2_ind_regime[i]],
                                       z_value = z_value[var2_ind_regime[i]],
                                       ind = var2_ind_regime[i],
                                       p_value = p_value[var2_ind_regime[i]]))
      rownames(tbl_var2[[i]]) <- "var"
      # covariances
      tbl_cov2[[i]] <- as.matrix(cbind(Estimate = cov2[, i],
                                       Std_Error = se[cov2_ind_regime[, i]],
                                       z_value = z_value[cov2_ind_regime[, i]],
                                       ind = cov2_ind_regime[, i],
                                       p_value = p_value[cov2_ind_regime[, i]]))
      rownames(tbl_cov2[[i]]) <- paste0("cov(", 1:(n_alt - 1), ")")
    }
  }
  
  # Coefficients for lambdas
  if (estimator == "2step")
  {
    tbl_lambda <- vector(mode = "list", length = n_regimes)
    for (i in 1:n_regimes)
    {
      n_coef2_tmp <- length(coef2_ind_regime[, i])
      estimate_lambda <- coef_lambda_row[[i]]
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
    names(tbl_lambda) <- regimes_names
  }
  
  # Save the tables
  tbl <- list(coef = tbl_coef,
              sigma = tbl_sigma,
              coef2 = tbl_coef2,
              var2 = tbl_var2,
              cov2 = tbl_cov2,
              lambda = tbl_lambda)
  
  # Store the output
  out <- list(tbl = tbl, se = se, p_value = p_value)
  
  return (out)
}