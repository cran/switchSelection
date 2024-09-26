# Split the data according to the equations
data_msel <- function(object, data = NULL)
{
  # Get some variables
  is1          <- object$other$is1
  is2          <- object$other$is2
  is3          <- object$other$is3
  is_het       <- object$other$is_het
  estimator    <- object$estimator
  n_eq         <- object$other$n_eq
  n_eq2        <- object$other$n_eq2
  n_eq3        <- object$other$n_eq3
  formula      <- object$formula
  formula2     <- object$formula2
  formula3     <- object$formula3
  formula_mean <- object$other$formula_mean
  formula_var  <- object$other$formula_var
  if (is.null(data))
  {
    object <- object$data
  }
  n_obs <- nrow(data)
  
  # Seperate the dataframes for the mean and variance equations
  # of the ordered equations
  df_mean <- NULL
  df_var  <- NULL
  if (is1)
  {
    df_mean <- vector(mode = "list", length = n_eq)
    df_var  <- vector(mode = "list", length = n_eq)
    for (i in seq_len(n_eq))
    {
      df_mean[[i]] <- model.frame(formula_mean[[i]], data, 
                                  na.action = na.pass)
      if (is_het[i])
      {
        df_var[[i]] <- model.frame(formula_var[[i]], data, 
                                   na.action = na.pass)
      }
    }
  }
  # dataframe for the continuous equations
  df2 <- NULL
  if (is2)
  {
    df2 <- vector(mode = "list", length = n_eq2)
    for (i in seq_len(n_eq2))
    {
      df2[[i]] <- model.frame(formula2[[i]], data, na.action = na.pass)
    }
  }
  # dataframe for the multinomial equation
  df3 <- NULL
  if (is3)
  {
    df3 <- model.frame(formula3, data, na.action = na.pass)
  }
  
  # Extract the variables of the ordinal equations
  z      <- matrix(NA)
  W_mean <- list(matrix())
  W_var  <- list(matrix())
  if (is1)
  {
    z      <- matrix(NA, nrow = n_obs, ncol = n_eq)
    W_mean <- vector("mode" = "list", length = n_eq)
    W_var  <- vector("mode" = "list", length = n_eq)
    for (i in seq_len(n_eq))
    {
      z[, i]      <- as.vector(df_mean[[i]][, 1])
      W_mean[[i]] <- as.matrix(df_mean[[i]][, -1, drop = FALSE])
      if (is_het[i])
      {
        W_var[[i]] <- as.matrix(df_var[[i]][, -1, drop = FALSE])
      }
      else
      {
        W_var[[i]] <- matrix()
      }
    }
    z[z == -1] <- NA
  }
  
  # Extract the variables of the continuous equations
  y <- matrix(NA)
  X <- list(matrix())
  if (is2)
  {
    y <- matrix(NA, nrow = n_obs, ncol = n_eq2)
    X <- vector("mode" = "list", length = n_eq2)
    for (i in seq_len(n_eq2))
    {
      y[, i]              <- as.vector(df2[[i]][, 1])
      X[[i]]              <- cbind(1, as.matrix(df2[[i]][, -1, drop = FALSE]))
      colnames(X[[i]])[1] <- "(Intercept)"
    }
  }
  
  # Extract the variables of the multinomial equation
  z_mn              <- vector(mode = "numeric")
  W_mn              <- matrix(NA)
  if (is3)
  {
    z_mn              <- df3[, 1]
    z_mn[z_mn == -1]  <- NA
    W_mn              <- as.matrix(cbind(1, df3[, -1, drop = FALSE]))
    colnames(W_mn)[1] <- "(Intercept)"
  }
  
  # Store the data to the output
  out <- list(W_mean = W_mean, W_var = W_var, X    = X,   W_mn = W_mn, 
              z      = z,      y     = y,     z_mn = z_mn)
  return(out)
}