# Store parameters for msel
par_msel <- function(object, par, type = "par")
{
  # Get some variables
  is_het           <- object$other$is_het
  is_marginal      <- object$other$is_marginal
  n_coef2          <- object$control_lnL$n_coef2
  n_coef3          <- object$control_lnL$n_coef3
  n_sigma          <- object$other$n_sigma
  n_sigma2         <- object$other$n_sigma2
  coef_ind         <- object$ind$coef
  coef_var_ind     <- object$ind$coef_var
  cuts_ind         <- object$ind$cuts
  coef2_ind        <- object$ind$coef2
  sigma_ind        <- object$ind$sigma
  var2_ind         <- object$ind$var2
  cov2_ind         <- object$ind$cov2
  sigma2_ind       <- object$ind$sigma2
  sigma3_ind       <- object$ind$sigma3
  marginal_par_ind <- object$ind$marginal_par
  marginal_par_n   <- object$control_lnL$marginal_par_n
  sigma_omit       <- object$control_lnL$sigma_omit
  coef3_ind        <- object$ind$coef3
  y_names          <- object$other$y_names
  z_names          <- object$other$z_names
  n_eq             <- object$control_lnL$n_eq
  n_eq2            <- object$control_lnL$n_eq2
  n_eq3            <- object$other$n_eq3
  n_regimes        <- object$control_lnL$n_regimes
  regimes_pair     <- object$other$regimes_pair
  colnames_X       <- object$other$colnames_X
  colnames_W_mean  <- object$other$colnames_W_mean
  colnames_W_var   <- object$other$colnames_W_var
  colnames_W_mn    <- object$other$colnames_W_mn
  estimator        <- object$estimator
  degrees          <- object$degrees
  is1              <- object$other$is1
  is2              <- object$other$is2
  is3              <- object$other$is3
  type3            <- object$type3
  
  # Coefficients of the ordered equations
  coef     <- vector(mode = "list", length = n_eq)
  coef_var <- vector(mode = "list", length = n_eq)
  cuts     <- vector(mode = "list", length = n_eq)
  if (is1)
  {
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
    names(coef)     <- z_names
    names(coef_var) <- z_names
    names(cuts)     <- z_names
  }
  
  # Coefficients of the continuous equations
  coef2 <- list()
  if (is2)
  {
    coef2 <- vector(mode = "list", length = n_eq2)
    for (i in 1:n_eq2)
    {
      coef2[[i]] <- matrix(nrow = n_regimes[i], ncol = n_coef2[i])
      for (j in 1:n_regimes[i])
      {
        coef2[[i]][j ,] <- par[coef2_ind[[i]][j ,]]
      }
      rownames(coef2[[i]]) <- paste0("regime ", 0:(n_regimes[i] - 1))
      colnames(coef2[[i]]) <- colnames_X[[i]]
    }
    names(coef2) <- y_names
  }
  
  # Variances of the continuous equations
  var2 <- list()
  if (is2 & (estimator == "ml"))
  {
    var2 <- vector(mode = "list", length = n_eq2)
    for (i in 1:n_eq2)
    {
      var2[[i]] <- par[var2_ind[[i]]]
      names(var2[[i]]) <- paste0("regime ", 0:(n_regimes[i] - 1))
    }
    names(var2)  <- y_names
  }
  
  # Covariances between the continuous and the ordered equations
  cov2 <- list()
  if (is1 & is2 & (estimator == "ml"))
  {
    cov2 <- vector(mode = "list", length = n_eq2)
    for (i in 1:n_eq2)
    {
      cov2[[i]] <- matrix(nrow = n_regimes[i], ncol = n_eq)
      for (j in 1:n_regimes[i])
      {
        cov2[[i]][j ,] <- par[cov2_ind[[i]][j ,]]
      }
      rownames(cov2[[i]]) <- paste0("regime ", 0:(n_regimes[i] - 1))
      colnames(cov2[[i]]) <- z_names
    }
    names(cov2)  <- y_names
  }
  
  # Covariances between the ordered equations
  sigma <- matrix(NA, nrow = n_eq, ncol = n_eq)
  sigma_vec_ind <- numeric()
  if (is1)
  {
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
          if (sigma_omit[i, j] == 0)
          {
            sigma_vec_ind[counter, ] <- c(i, j)
            sigma[i, j] <- sigma_vec[counter]
            sigma[j, i] <- sigma[i, j]
            counter <- counter + 1
          }
        }
      }
    }
    colnames(sigma) <- z_names
    rownames(sigma) <- z_names
  }
  
  # Covariances between the continuous equations
  sigma2  <- list()
  counter <- 1
  if ((n_eq2 >= 2) & (estimator == "ml"))
  {
    sigma2 <- vector(mode = "list", length = n_sigma2)
    for (i in 2:n_eq2)
    {
      for (j in 1:(i - 1))
      {
        sigma2[[counter]] <- par[sigma2_ind[[counter]]]
        names(sigma2)[counter] <- paste0(y_names[i], ", ", y_names[j])
        if (length(regimes_pair[[counter]]) > 0)
        {
          for (t in 1:nrow(regimes_pair[[counter]]))
          {
            names(sigma2[[counter]])[t] <- paste0("(",
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
  
  # Marginal distribution parameters
  marginal_par <- vector(mode = "list", length = n_eq)
  if (is_marginal)
  {
    for (i in 1:n_eq)
    {
      if (marginal_par_n[i] > 0)
      {
        marginal_par[[i]]        <- par[marginal_par_ind[[i]]]
        names(marginal_par[[i]]) <- paste0("par", 1:marginal_par_n[i])
      }
    }
    names(marginal_par) <- z_names
  }
  
  
  # Coefficients of the multinomial equations
  coef3 <- numeric()
  if (is3)
  {
    coef3 <- matrix(NA, nrow = n_eq3 - 1, ncol = n_coef3)
    rownames(coef3) <- 1:(n_eq3 - 1)
    for (i in 1:(n_eq3 - 1))
    {
      coef3[i, ] <- par[coef3_ind[i, ]]
    }
    colnames(coef3) <- colnames_W_mn
    rownames(coef3) <- paste0("alt", 0:(n_eq3 - 2))
  }
  
  # Covariances of the multinomial equation
  sigma3 <- matrix(1)
  if(is3 & (type3 == "probit"))
  {
    sigma3       <- matrix(NA, nrow = n_eq3 - 1, ncol = n_eq3 - 1)
    sigma3[1, 1] <- 1
    sigma3_vec   <- NULL
    if (n_eq3 > 2)
    {
      counter <- 1;
      sigma3_vec <- par[sigma3_ind]
      for (i in 1:(n_eq3 - 1))
      {
        for (j in 1:i)
        {
          if (!((i == 1) & (j == 1)))
          {
            sigma3[i, j] <- sigma3_vec[counter]
            sigma3[j, i] <- sigma3[i, j];
            counter <- counter + 1
          }
        }
      }
    }
    colnames(sigma3) <- paste0("alt", 0:(n_eq3 - 2))
    rownames(sigma3) <- colnames(sigma3)
  }
  
  # If the model should be updated
  if (type == "object")
  {
    object$par                 <- par
    object$coef                <- coef
    object$coef_var            <- coef_var
    object$cuts                <- cuts
    object$coef2               <- coef2
    object$var2                <- var2
    object$cov2                <- cov2
    object$sigma               <- sigma
    object$sigma2              <- sigma2
    object$marginal_par        <- marginal_par
    object$other$sigma_vec_ind <- sigma_vec_ind
    object$coef3               <- coef3
    object$sigma3              <- sigma3
    
    return(object)
  }
  
  # Store the output
  out <- list(par           = par,           coef         = coef, 
              coef_var      = coef_var,      cuts         = cuts,                 
              coef2         = coef2,         var2         = var2,                 
              cov2          = cov2,          sigma        = sigma,               
              sigma2        = sigma2,        marginal_par = marginal_par, 
              sigma_vec_ind = sigma_vec_ind, coef3        = coef3,
              sigma3        = sigma3)
  
  # Return the result
  return(out)
}