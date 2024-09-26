# Get the names associated with the data
names_msel <- function(object, data = NULL)
{
  # Get some variables
  is1       <- object$other$is1
  is2       <- object$other$is2
  is3       <- object$other$is3
  is_het    <- object$other$is_het
  estimator <- object$estimator
  n_eq      <- object$other$n_eq
  n_eq2     <- object$other$n_eq2
  n_eq3     <- object$other$n_eq3
  formula   <- object$formula
  formula2  <- object$formula2
  formula3  <- object$formula3
  y         <- object$y
  z         <- object$z  
  z_mn      <- object$z_mn
  X         <- object$X
  W_mean    <- object$W_mean
  W_var     <- object$W_var
  W_mn      <- object$W_mn
  
  if (is.null(data))
  {
    object <- object$data
  }
  n_obs <- nrow(data)
  
  # Get names of the dependent variables
  z_names      <- NA
  y_names      <- NA
  z_mn_names   <- NA
  if (is1)
  {
    z_names <- rep(NA, n_eq)
    for (i in seq_len(n_eq))
    {
      z_names[i] <- as.character(formula[[i]][[2]])
    }
  }
  if (is2)
  {
    y_names <- rep(NA, n_eq2)
    for (i in seq_len(n_eq2))
    {
      y_names[i] <- all.vars(formula2[[i]])[1]
    }
  }
  if (is3)
  {
    z_mn_names <- all.vars(formula3)[1]
  }
  
  # Get names of the independent variables
  colnames_X                <- lapply(X,      colnames)
  colnames_W_mean           <- lapply(W_mean, colnames)
  colnames_W_var            <- lapply(W_var,  colnames)
  colnames_W_mn             <- colnames(W_mn)
  
  # Return the result
  out <- list(z_names         = z_names,    y_names    = y_names, 
              z_mn_names      = z_mn_names, colnames_X = colnames_X,
              colnames_W_mean = colnames_W_mean, 
              colnames_W_var  = colnames_W_var, 
              colnames_W_mn   = colnames_W_mn)
}