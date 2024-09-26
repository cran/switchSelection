# Get the groups related variables
groups_msel <- function(object, data, groups, groups2, groups3)
{
  # Get some variables
  is1         <- object$other$is1
  is2         <- object$other$is2
  is3         <- object$other$is3
  is_het      <- object$other$is_het
  estimator   <- object$estimator
  n_eq        <- object$other$n_eq
  n_eq2       <- object$other$n_eq2
  formula     <- object$formula
  formula2    <- object$formula2
  formula3    <- object$formula3
  is_na_group <- object$other$is_na_group
  if (is.null(data))
  {
    object <- object$data
  }
  n_obs <- nrow(data)
  
  # The number of observations
  n_obs <- nrow(data)
  
  # Get the dependent variables of the ordinal equations
  z <- matrix(NA)
  if (is1)
  {
    z <- matrix(NA, nrow = n_obs, ncol = n_eq)
    for (i in 1:n_eq)
    {
      z[, i] <- data[, all.vars(formula[[i]])[1]]
    }
  }
  
  # Get the dependent variables of the continuous equations
  y <- matrix(NA)
  if (is2)
  {
    y <- matrix(NA, nrow = n_obs, ncol = n_eq2)
    for (i in 1:n_eq2)
    {
      y[, i] <- data[, all.vars(formula2[[i]])[1]]
    }
  }
  
  # Get the dependent variable of the muiltinomial equation
  z_mn <- vector(mode = "numeric")
  if (is3)
  {
    z_mn <- data[, all.vars(formula3)[1]]
  }
  
  # Calculate the number of cuts for each equation
  n_cuts_eq <- rep(0, n_eq)
  if (is1)
  {
    for (i in seq_len(n_eq))
    {
      n_cuts_eq[i] <- max(na.omit(z[, i]))
    }
  }

  # Calculate the number of the alternatives in the 
  # multinomial equation
  alt   <- 0
  n_eq3 <- 0
  if (is3)
  {
    alt <- sort(unique(na.omit(z_mn)))
    n_eq3 <- length(alt[alt != -1])
  }
  
  # Set the groups by default
  groups_all <- NULL
  if ((is_na_group[1] & is1) | 
      (is_na_group[2] & is2) |
      (is_na_group[3] & is3))
  {
    # Generate all groups
    groups_pows <- NULL
    if (is1)
    {
      groups_pows <- n_cuts_eq + 1
    }
    if (is2)
    {
      groups_pows <- c(groups_pows, rep(1, n_eq2))
    }
    if (is3)
    {
      groups_pows <- c(groups_pows, n_eq3)
    }
    groups_all <- as.matrix(t(hpa::polynomialIndex(groups_pows) - 1))

    # Split the groups
    counter <- 1
    if (is1)
    {
      groups <- groups_all[, counter:n_eq, drop = FALSE]
      counter <- n_eq + 1
    }
    if (is2)
    {
      groups2 <- groups_all[, counter:(counter + n_eq2 - 1), drop = FALSE]
      counter <- counter + n_eq2
    }
    if (is3)
    {
      groups3 <- groups_all[, counter, drop = TRUE]
      counter <- counter + 1
    }
  }
  
  # Get the indexes of the observations for each group
  ind_g <- NULL
  
  # Get the indexes of the observations corresponding to each group
  groups123 <- NULL
  data123   <- NULL
  if (is1)
  {
    z_tmp           <- z
    z_tmp[is.na(z)] <- -1
    groups123       <- cbind(groups123, groups)
    data123         <- cbind(data123, z_tmp)
  }
  if (is2)
  {
    y_tmp                        <- y
    y_tmp[is.na(y)]              <- -1
    y_tmp[!is.na(y)]             <- 0
    groups2_tmp                  <- groups2
    groups2_tmp[groups2_tmp > 0] <- 0
    groups123                    <- cbind(groups123, groups2_tmp)
    data123                      <- cbind(data123, y_tmp)
  }
  if (is3)
  {
    z_mn_tmp              <- z_mn
    z_mn_tmp[is.na(z_mn)] <- -1
    groups123             <- cbind(groups123, groups3)
    data123               <- cbind(data123, z_mn_tmp)
  }

  # Get the indexes
  ind_g <- findGroup(data123, groups123)
  
  # Calculate the number of groups
  n_groups <- length(ind_g)
  
  # Calculate the number of observations in each group
  n_obs_g <- rep(0, n_groups)
  for (i in seq_len(n_groups))
  {
    n_obs_g[i] <- length(ind_g[[i]])
  }
  n_obs_g_0 <- n_obs_g == 0
  
  # Get the groups without the observable categories
  g_unobs_1 <- rep(TRUE, n_groups)
  g_unobs_2 <- rep(TRUE, n_groups)
  g_unobs_3 <- rep(TRUE, n_groups)
  for (i in seq_len(n_groups))
  {
    if (is1)
    {
      g_unobs_1[i] <- all(groups[i, , drop = FALSE] == -1)
    }
    if (is2)
    {
      g_unobs_2[i] <- all(groups2[i, , drop = FALSE] == -1)
    }
  }
  if (is3)
  {
    g_unobs_3 <- groups3 == -1
  }
  g_unobs <- g_unobs_1 & g_unobs_2 & g_unobs_3
  
  # Remove the groups with zero observations
  ind_g_include <- !(n_obs_g_0 | g_unobs)
  if (is1)
  {
    groups <- groups[ind_g_include, , drop = FALSE]
  }
  if (is2)
  {
    groups2 <- groups2[ind_g_include, , drop = FALSE]
  }
  if (is3)
  {
    groups3 <- groups3[ind_g_include]
  }
  n_obs_g            <- n_obs_g[ind_g_include]
  ind_g              <- ind_g[ind_g_include]
  n_groups           <- length(n_obs_g)
  
  # Calculate the number of the observed equations per group
  n_eq_g     <- rep(0, n_groups)
  n_eq_all_g <- n_eq_g
  n_eq2_g    <- rep(0, n_groups)
  # ordered equations
  if (is1)
  {
    for (i in seq_len(n_groups))
    {
      n_eq_g[i] <- sum(groups[i, , drop = FALSE] != -1)
    }
    n_eq_all_g <- n_eq_g
  }
  # continuos equations
  if (is2)
  {
    n_eq2_g <- vector(mode = "numeric", length = n_groups)
    for (i in seq_len(n_groups))
    {
      n_eq2_g[i]    <- sum(groups2[i, , drop = FALSE] != -1)
      n_eq_all_g[i] <- n_eq_g[i] + n_eq2_g[i]
    }
  }
  
  # Get the indexes of the observed equations for each group
  ind_eq     <- vector(mode = "list", length = n_groups)
  ind_eq2    <- vector(mode = "list", length = n_groups)
  ind_eq3    <- vector(mode = "list", length = n_groups)
  ind_eq_all <- vector(mode = "list", length = n_groups)
  # ordered equations
  if (is1)
  {
    for (i in seq_len(n_groups))
    {
      ind_eq[[i]] <- which(groups[i, , drop = FALSE] != -1)
    }
    ind_eq_all <- ind_eq
  }
  # continuous equations
  if (is2)
  {
    for (i in seq_len(n_groups))
    {
      ind_eq2[[i]]    <- which(groups2[i, , drop = FALSE] != -1)
      ind_eq_all[[i]] <- c(ind_eq[[i]], ind_eq2[[i]] + n_eq)
    }
  }
  # multinomial equations
  # does not include ind_eq_all for now but
  # may be added later
  if (is3)
  {
    ind_eq3 <- which(groups3 != -1)
  }
  
  # Get the vector of all indexes
  ind_g_all <- NULL
  for (i in 1:n_groups)
  {
    ind_g_all <- c(ind_g_all, ind_g[[i]])
  }
  ind_g_all <- sort(ind_g_all)
  # Important: the vector of all indexes of the second step may be
  #            used to select the scores from the first step
  
  # Logical values related to the observations 
  # associated with some group
  ind_g_logical            <- rep(FALSE, n_obs)
  ind_g_logical[ind_g_all] <- TRUE
  ind_g_logical_cumsum     <- cumsum(ind_g_logical)
  
  # Preserve only the observations associated with some group
  data <- data[ind_g_all, , drop = FALSE]
  for (i in 1:n_groups)
  {
    ind_g[[i]] <- ind_g_logical_cumsum[ind_g[[i]]]
  }
  n_obs <- nrow(data)
  
  # Return the results
  out <- list(groups     = groups,     groups2    = groups2,
              groups3    = groups3,    ind_g      = ind_g,
              ind_eq     = ind_eq,     ind_eq2    = ind_eq2,
              ind_eq3    = ind_eq3,    ind_eq_all = ind_eq_all,
              n_cuts_eq  = n_cuts_eq,  n_eq3      = n_eq3,
              n_groups   = n_groups,   data       = data,
              n_obs      = n_obs,      n_obs_g    = n_obs_g,
              n_eq_g     = n_eq_g,     n_eq2_g    = n_eq2_g,
              n_eq_all_g = n_eq_all_g, ind_g_all  = ind_g_all)
}