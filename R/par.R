# Store parameters for mvoprobit
par_mvoprobit <- function(par, 
                          is_marginal, is_het,
                          n_coef2, n_sigma, n_sigma2,
                          coef_ind, coef_var_ind, cuts_ind, coef2_ind, 
                          sigma_ind, var2_ind, cov2_ind, sigma2_ind, 
                          marginal_par_ind, marginal_par_n,
                          y_names, z_names, 
                          n_eq, n_eq2, n_regimes, regimes_pair,
                          colnames_X, colnames_W_mean, colnames_W_var)
                                
{
  # Coefficients of ordered equations
  coef <- vector(mode = "list", length = n_eq)
  coef_var <- vector(mode = "list", length = n_eq)
  cuts <- vector(mode = "list", length = n_eq)
  for (i in 1:n_eq)
  {
    coef[[i]] <- par[coef_ind[[i]]]
    names(coef[[i]]) <- colnames_W_mean[[i]]
    if (is_het[i])
    {
      coef_var[[i]] <- par[coef_var_ind[[i]]]
      names(coef_var[[i]]) <- colnames_W_var[[i]]
    }
    cuts[[i]] <- par[cuts_ind[[i]]]
    names(cuts[[i]]) <- paste0("cut ", 1:length(cuts[[i]]))
  }
  names(coef) <- z_names
  names(coef_var) <- z_names
  names(cuts) <- z_names
  
  # Coefficients of continuous equations and covariances between
  # continuous and ordered equations
  coef2 <- NULL
  var2 <- NULL
  cov2 <- NULL
  if (n_eq2 > 0)
  {
    coef2 <- vector(mode = "list", length = n_eq2)
    var2 <- vector(mode = "list", length = n_eq2)
    cov2 <- vector(mode = "list", length = n_eq2)
    for (i in 1:n_eq2)
    {
      coef2[[i]] <- matrix(nrow = n_regimes[i], ncol = n_coef2[i])
      var2[[i]] <- par[var2_ind[[i]]]
      names(var2[[i]]) <- paste0("regime ", 0:(n_regimes[i] - 1))
      cov2[[i]] <- matrix(nrow = n_regimes[i], ncol = n_eq)
      for (j in 1:n_regimes[i])
      {
        coef2[[i]][j ,] <- par[coef2_ind[[i]][j ,]]
        cov2[[i]][j ,] <- par[cov2_ind[[i]][j ,]]
      }
      rownames(coef2[[i]]) <- paste0("regime ", 0:(n_regimes[i] - 1))
      colnames(coef2[[i]]) <- colnames_X[[i]]
      rownames(cov2[[i]]) <- paste0("regime ", 0:(n_regimes[i] - 1))
      colnames(cov2[[i]]) <- z_names
    }
    names(coef2) <- y_names
    names(cov2) <- y_names
    names(var2) <- y_names
  }
  
  # Covariance matrix between ordered equations
  sigma <- matrix(NA, nrow = n_eq, ncol = n_eq)
  sigma_vec_ind <- NULL
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
  colnames(sigma) <- z_names
  rownames(sigma) <- z_names
  
  # Covariance matrix between continuous equations
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
        names(sigma2)[counter] <- paste0(y_names[i], ", ", y_names[j])
        for (t in 1:nrow(regimes_pair[[counter]]))
        {
          names(sigma2[[counter]])[t] <- paste0("(",
                                                regimes_pair[[counter]][t, 1], 
                                                ",",
                                                regimes_pair[[counter]][t, 2],
                                                ")")
        }
        counter <- counter + 1
      }
    }
  }
  
  # Marginal distribution parameters
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
  
  # Store the output
  out <- list(par = par,
              coef = coef, coef_var = coef_var,
              cuts = cuts,
              coef2 = coef2, var2 = var2, cov2 = cov2,
              sigma = sigma, sigma2 = sigma2,
              marginal_par = marginal_par,
              sigma_vec_ind = sigma_vec_ind)
}

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

par_mnprobit <- function(par, n_alt, n_coef, n_coef2, n_regimes,
                         coef_ind_alt, sigma_ind, 
                         coef2_ind_regime, var2_ind_regime, cov2_ind_regime,
                         coef_lambda, is2, degrees, estimator,
                         alt_names, regimes_names, colnames_W, colnames_X)
{
  # Store the coefficients
  coef <- matrix(NA, nrow = n_coef, ncol = n_alt - 1)
  for (i in 1:(n_alt - 1))
  {
    coef[, i] <- par[coef_ind_alt[, i]]
  }
  rownames(coef) <- colnames_W
  colnames(coef) <- alt_names[-n_alt]
  
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
  rownames(sigma) <- alt_names[-n_alt]
  colnames(sigma) <- alt_names[-n_alt]
  
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
    rownames(coef2) <- colnames_X
    colnames(coef2) <- regimes_names
    names(var2) <- regimes_names
    colnames(cov2) <- regimes_names
    rownames(cov2) <- alt_names[-n_alt]
  }
  
  # Store the output
  out <- list(par = par,
              coef = coef, coef_lambda = coef_lambda,
              sigma = sigma, sigma_vec = sigma_vec)
  if (!is.null(degrees))
  {
    out$coef2 <- coef2
    out$var2 <- var2
    out$cov2 <- cov2
  }
  
  return (out)
}