# Estimate the asymptotic covariance matrix of the maximum-likelihood estimator
vcov_ml <- function(object, type, n_cores, n_sim)
{
  # Validate a type of the estimator
  type <- tolower(type)
  
  # List of output values
  out <- list(vcov = NULL)
  
  # Get some variables
  par            <- object$par
  cov_type       <- type
  estimator      <- object$estimator
  control_lnL    <- object$control_lnL
  n_par          <- object$other$n_par
  regularization <- object$other$regularization
  n_obs          <- object$other$n_obs
  is_cluster     <- object$other$is_cluster
  n_cluster      <- object$other$n_cluster
  ind_cluster    <- object$other$ind_cluster
  rm(object)

  # Estimate asymptotic covariance matrix
  H         <- NULL
  J         <- NULL
  J_cluster <- NULL
  cov       <- matrix(NA, nrow = n_par, ncol = n_par)
  if (!is.matrix(cov_type))
  {
    if (cov_type != "no")
    {
      if (cov_type %in% c("sandwich", "hessian", "mm"))
      {
        tryCatch(
          {
            cov_type_old <- cov_type
            cov_type     <- "gop"
            gr.args      <- list(n_sim       = n_sim, 
                                 n_cores     = n_cores,
                                 control_lnL = control_lnL, 
                                 out_type    = "grad")
            H            <- NULL
            if (cov_type_old == "mm")
            {
              gr.args$regularization <- regularization
            }
            H        <- gena::gena.hessian(gr      = lnL_msel,
                                           par     = par,
                                           gr.args = gr.args,
                                           is.gc   = TRUE)
            out$H    <- H
            cov_type <- cov_type_old
          },
          error = function(e) {
            warning(paste0("Problems with numeric hessian calculation. ",
                           "Therefore 'cov_type' has been changed to 'gop'."))
          }
        )
        gc()
      }
      if (cov_type %in% c("sandwich", "gop", "mm"))
      {
        J <- NULL
        if (cov_type == "mm")
        {
          J <- lnL_msel(par            = par,
                        n_sim          = n_sim, 
                        n_cores        = n_cores,
                        control_lnL    = control_lnL, 
                        out_type       = "jac",
                        regularization = regularization)
        }
        else
        {
          J <- lnL_msel(par            = par,
                        n_sim          = n_sim, 
                        n_cores        = n_cores,
                        control_lnL    = control_lnL, 
                        out_type       = "jac")
        }
        out$J <- J
        gc()
        
        # Estimate clustered scores
        J_cluster <- matrix(0, nrow = n_cluster, ncol = n_par)
        if (is_cluster)
        {
          for (i in 1:n_cluster)
          {
            J_cluster[i, ] <- colSums(J[ind_cluster == i, , drop = FALSE])
          }
        }
        out$J_cluster <- J_cluster
      }
      if ((cov_type == "sandwich") | (cov_type == "mm"))
      {
        if (is_cluster)
        {
          cov <- tcrossprod(qr.solve(H, t(J_cluster), tol = 1e-16))
        }
        else
        {
          cov <- tcrossprod(qr.solve(H, t(J), tol = 1e-16))
        }
      }
      if (cov_type == "hessian")
      {
        cov <- -qr.solve(H, tol = 1e-16)
      }
      
      if (any(is.na(cov)) & (cov_type != "gop"))
      {
        warning(paste0("Can't calculate the covariance matrix of type '", 
                       cov_type, "'. ", 
                       "Therefore gop covariance matrix will be used instead."))
        cov_type <- "gop"
      }
      
      if (cov_type == "gop")
      {
        tryCatch(
        {
          H     <- crossprod(J)
          J_tmp <- qr.R(qr(J, LAPACK = TRUE))
          cov   <- chol2inv(diag(sign(diag(J_tmp))) %*% J_tmp)
          out$H <- H
        },
        error = function(e) {
          warning(paste0("Problems with numeric hessian calculation. ",
                         "Therefore 'cov_type' has been changed to 'no'."))
        }
       )
      }
    }
  } 
  
  # Aggregate some output
  out$vcov <- cov
  gc()
  
  return(out)
}

# Estimate the asymptotic covariance matrix of the two-step estimator
vcov_2step <- function(object, type = "default", n_cores = 1, n_sim = 10000,
                       is_cov_simple = FALSE)
{
  # Some variables
  n_eq        <- object$other$n_eq
  n_eq2       <- object$other$n_eq2
  n_eq3       <- object$other$n_eq3
  n_par       <- object$other$n_par
  n_obs       <- object$other$n_obs
  is1         <- object$other$is1
  is3         <- object$other$is3
  type3       <- object$type3
  n_sigma     <- object$other$n_sigma
  is_cluster  <- object$other$is_cluster
  n_cluster   <- object$other$n_cluster
  ind_cluster <- object$other$ind_cluster
  
  # Calculate the scores of the second step
  scores <- scores_2step(object, n_cores = n_cores, n_sim = n_sim)
  
  # Calculate the derivatives with respect to the scores of the second step
  scores_list <- NULL
  if (is_cov_simple)
  {
    # Preliminary values
    n_groups   <- object$other$n_groups
    coef2_ind  <- object$ind$coef2
    groups2    <- object$groups2
    y_names    <- object$other$y_names
    colnames_x <- object$other$colnames_X
    par        <- object$par
    y_pred     <- object$y_pred
    y_names    <- object$other$y_names
    y          <- object$y
    
    # Manual calculations with some simplifications
    H <- matrix(0, nrow = n_par, ncol = n_par)
    for (v in 1:n_eq2)
    {
      X           <- object$X[[v]]
      colnames(X) <- colnames_x[[v]]
      X_Y_names   <- y_names[which(y_names %in% colnames_x[[v]])]
      if (length(X_Y_names) > 0)
      {
        X[, X_Y_names] <- y_pred[, X_Y_names]
      }
      for (i in 1:n_groups)
      {
        # Scores with respect to their own coefficients
        if (groups2[i, v] >= 0)
        {
          X0          <- X[object$other$ind_g[[i]], , drop = FALSE]
          ind         <- coef2_ind[[v]][groups2[i, v] + 1, ]
          H[ind, ind] <- H[ind, ind] - crossprod(X0)
        }
        # Scores with respect to other coefficients
        if (v < n_eq2)
        {
          for (v1 in (v + 1):n_eq2)
          {
            if (groups2[i, v1] != -1)
            {
              indy         <- which(colnames_x[[v1]] == y_names[v])
              X1           <- object$X[[v1]][object$other$ind_g[[i]], , drop = FALSE]
              colnames(X1) <- colnames_x[[v1]]
              X1_Y_names   <- y_names[which(y_names %in% colnames_x[[v1]])]
              if (length(X1_Y_names) > 0)
              {
                X1[, X1_Y_names] <- y_pred[object$other$ind_g[[i]], 
                                           X1_Y_names, drop = FALSE]
              }
              if (length(indy) > 0)
              {
                # (y-yh)' * x
                ind1         <- coef2_ind[[v1]][groups2[i, v1] + 1, ]
                H_new        <- (t(X0) %*% X1) * par[ind1[indy]]
                H[ind1, ind] <- H[ind1, ind]  - t(H_new)
                # (y - yh) * x'
                H_new <- colSums(sweep(X0, MARGIN = 1,
                                 STATS = y[object$other$ind_g[[i]], v1] -
                                         y_pred[object$other$ind_g[[i]], v1],
                                 FUN = "*"))
                H[ind1[indy], ind] <- H[ind1[indy], ind] + H_new
              }
            }
          }
        }
        gc()
      }
    }
    scores_list <- list(val  = matrix(colSums(scores), ncol = 1),
                        grad = H)
  }
  else
  {
    # Full calculations which may take much time in large samples
    scores_list <- deriv_msel(object  = object, 
                              fn      = scores_2step,
                              fn_args = list(type     = "aggregate",
                                             n_cores  = n_cores, 
                                             n_sim    = n_sim))
  }
  scores_sum           <- scores_list$val
  scores_jac           <- scores_list$grad
  rownames(scores_jac) <- colnames(scores_jac)
  names(scores_sum)    <- colnames(scores_jac)

  # Collect the estimates of the first step
  if (is1 | is3)
  {
    # Jacobian and average Hessian
    H1 <- object$model1$H
    J1 <- object$model1$J[object$other$ind_g_all, ]

    # Indexes of the first step in the final model
    coef_ind          <- object$ind$coef
    cuts_ind          <- object$ind$cuts
    coef_var_ind      <- object$ind$coef_var
    marginal_par_ind  <- object$ind$marginal_par
    sigma_ind         <- object$ind$sigma
    coef3_ind         <- object$ind$coef3
    sigma3_ind        <- object$ind$sigma3
    
    # Indexes of the first step in first step model
    coef_ind1         <- object$model1$ind$coef
    cuts_ind1         <- object$model1$ind$cuts
    coef_var_ind1     <- object$model1$ind$coef_var
    marginal_par_ind1 <- object$model1$ind$marginal_par
    sigma_ind1        <- object$model1$ind$sigma
    coef3_ind1        <- object$model1$ind$coef3
    sigma3_ind1       <- object$model1$ind$sigma3
    
    # Assign the scores from the first step
      # ordered equations
    if (is1)
    {
      l  <- list()
      l1 <- list()
      if (n_sigma > 0)
      {
        l[[1]]  <- sigma_ind
        l1[[1]] <- sigma_ind1
      }
      for (i in 1:n_eq)
      {
        l[[length(l) + 1]]   <- coef_ind[[i]]
        l1[[length(l1) + 1]] <- coef_ind1[[i]]
        l[[length(l) + 1]]   <- cuts_ind[[i]]
        l1[[length(l1) + 1]] <- cuts_ind1[[i]]
        if (object$other$is_het[i])
        {
          l[[length(l) + 1]]   <- coef_var_ind[[i]]
          l1[[length(l1) + 1]] <- coef_var_ind1[[i]]
        }
        if (object$other$marginal_par_n[i] > 0)
        {
          l[[length(l) + 1]]   <- marginal_par_ind[[i]]
          l1[[length(l1) + 1]] <- marginal_par_ind1[[i]]
        }
        for(j in seq_len(length(l)))
        {
          scores[, l[[j]]] <- J1[, l1[[j]]]
          for(t in 1:length(l))
          {
            scores_jac[l[[j]], l[[t]]] <- H1[l1[[j]], l1[[t]]]
          }
        }
      }
    }
      # multinomial equations
    if (is3)
    {
      if (type3 == "probit")
      {
        scores[, sigma3_ind] <- J1[, sigma3_ind1]
      }
      for (i in seq_len(n_eq3 - 1))
      {
        scores[, coef3_ind[i, ]] <- J1[, coef3_ind1[i, ]]
        for(j in 1:(n_eq3 - 1))
        {
          scores_jac[coef3_ind[i, ], coef3_ind[j, ]] <- H1[coef3_ind1[i, ], 
                                                           coef3_ind1[j, ]]
          if (type3 == "probit")
          {
            scores_jac[coef3_ind[j, ], sigma3_ind] <- H1[coef3_ind1[j, ], 
                                                         sigma3_ind1]
            scores_jac[sigma3_ind, coef3_ind[j, ]] <- H1[sigma3_ind1, 
                                                         coef3_ind1[j, ]]
            scores_jac[sigma3_ind, sigma3_ind]     <- H1[sigma3_ind1, 
                                                         sigma3_ind1]
          }
        }
      }
    }
    gc()
  }

  # Estimate clustered scores
  scores_cluster <- matrix(0, nrow = n_cluster, ncol = n_par)
  if (is_cluster)
  {
    for (i in 1:n_cluster)
    {
      scores_cluster[i, ] <- colSums(scores[ind_cluster == i, , drop = FALSE])
    }
  }

  # Estimate asymptotic covariance matrix
  vcov <- NULL
  if (is_cluster)
  {
    vcov <- tcrossprod(qr.solve(scores_jac, t(scores_cluster), tol = 1e-16))
  }
  else
  {
    vcov <- tcrossprod(qr.solve(scores_jac, t(scores), tol = 1e-16))
  }
  colnames(vcov) <- NULL
  rownames(vcov) <- NULL

  # Return the results
  return(list(scores_jac     = scores_jac,
              scores         = scores,
              scores_sum     = scores_sum,
              vcov           = vcov,
              scores_cluster = scores_cluster))
}

# Estimate scores associated with the second step
scores_2step <- function(object, type = "obs", n_cores = 1, n_sim = 1000)
{
  # Get some variables
  n_par           <- object$other$n_par
  n_eq            <- object$other$n_eq
  n_eq2           <- object$other$n_eq2
  n               <- object$other$n_obs
  groups          <- object$groups
  groups2         <- object$groups2
  groups3         <- object$groups3
  n_groups        <- object$other$n_groups
  ind_g           <- object$ind$g
  coef2_ind       <- object$ind$coef2
  degrees         <- object$degrees
  is1             <- object$other$is1
  is2             <- object$other$is2
  is3             <- object$other$is3

  # Matrix to store the scores
  scores <- matrix(0, nrow = n, ncol = n_par)
  
  # Get conditional predictions with a new model
  for (i in 1:n_groups)
  {
    groups_i  <- NA
    groups3_i <- NA
    if (is1)
    {
      groups_i <- groups[i, ]
    }
    if (is3)
    {
      groups3_i <- groups3[i]
    }
    scores_group <- predict(object, 
                            type    = "val", 
                            newdata = object$data[ind_g[[i]], , drop = FALSE],
                            group   = groups_i, 
                            group2  = groups2[i, ],
                            group3  = groups3_i,
                            control = list(is_scores = TRUE))
    for (v in 1:n_eq2)
    {
      coef2_ind_regime                     <- coef2_ind[[v]][groups2[i, v] + 1, ]
      scores[ind_g[[i]], coef2_ind_regime] <- scores_group[[v]]
    }
  }
  
  if (type == "aggregate")
  {
    scores_sum <- as.matrix(colSums(scores), ncol = 1)
    return(scores_sum)
  }
  
  return(scores)
}

#' Calculate Variance-Covariance Matrix for a msel Object.
#' @description Return the variance-covariance matrix of the parameters of
#' msel model.
#' @param object an object of class \code{msel}.
#' @param ... further arguments (currently ignored).
#' @param type character representing the type of the asymptotic covariance 
#' matrix estimator. It takes the same values as \code{cov_type} parameter of
#' the \code{\link[switchSelection]{msel}} function.
#' @param n_sim integer representing the number of GHK draws when there are
#' more than 3 ordered equations. Otherwise alternative (much more efficient) 
#' algorithms will be used to calculate multivariate normal probabilities.
#' @param n_cores positive integer representing the number of CPU cores used for 
#' parallel computing. If possible it is highly recommended to set it equal to
#' the number of available physical cores especially when the system of
#' ordered equations has 2 or 3 equations.
#' @param recalculate logical; if \code{TRUE} then covariance matrix will be
#' recalculated even if 'type' is the same as 'cov_type' input argument
#' of the model.
#' @param cluster an object which takes the same values as the \code{cluster} 
#' argument of the \code{\link[switchSelection]{msel}} function. 
#' The only exception is \code{cluster = NULL} which implies 
#' \code{cluster = object$cluster}.
#' @details Argument \code{type} is closely related to the argument 
#' \code{cov_type} of \code{\link[switchSelection]{msel}} function. 
#' See 'Details' and 'Usage' sections of \code{\link[switchSelection]{msel}} 
#' for more information on \code{cov_type} argument.
#' 
#' The user may also estimate asymptotic covariance matrix of the parameters
#' of several models. For more information, see paragraph 'Cross-model tests'
#' in 'Details' section of \code{\link[switchSelection]{test_msel}}.
#' 
#' @return Returns numeric matrix which represents estimate of the asymptotic 
#' covariance matrix of model's parameters. If \code{object} is a list of
#' models then rows and columns of this matrix will have names \code{"mipj"}
#' indicating that corresponding element is associated with \code{j}-th 
#' parameter of the \code{i}-th model.
vcov.msel <- function(object, ..., 
                      type        = object$cov_type, 
                      n_cores     = object$other$n_cores,
                      n_sim       = object$other$n_sim,
                      recalculate = FALSE,
                      cluster     = NULL)
{
  # Validate dots
  if (length(list(...)) > 0)
  {
    warning("Additional arguments passed through ... are ignored.")   
  }
  
  # Deal with several objects
  if (!is(object = object, class2 = "msel"))
  {
    return (vcov_combine(object))
  }
  
  # Clusters
  if (!is.null(cluster))
  {
    object$cluster           <- cluster
    cluster_list             <- cluster_msel(object = object, cluster = cluster)
    object$other$ind_cluster <- cluster_list$ind_cluster
    object$other$n_cluster   <- cluster_list$n_cluster
    object$other$is_cluster  <- cluster_list$is_cluster
    recalculate              <- TRUE
  }
  
  # If no recalculations are required
  if ((type == object$cov_type) & !recalculate)
  {
    return (object$cov)
  }
  
  # Maximum-likelihood
  if (object$estimator == "ml")
  {
    return (vcov_ml(object  = object,  type  = type, 
                    n_cores = n_cores, n_sim = n_sim)$vcov)
  }
  
  # Two-step
  if (object$estimator == "2step")
  {
    return(vcov_2step(object = object,n_cores = n_cores, n_sim = n_sim)$vcov)
  }
  
  # Error
  stop("Incorrect 'object$estimator' value.")
}

#' Combine the asymptotic covariance matrices
#' @param object a list of objects of class \code{msel}.
#' @param type the same as the \code{type} argument of the 
#' \code{\link[switchSelection]{vcov.msel}} function.
vcov_combine <- function(object, type = "cov")
{
  # Get the number of objects
  n_objects <- length(object)
  
  # Total number of parameters and indexes of the parameters
  n_par     <- 0
  ind_par   <- vector(mode = "list", length = n_objects)
  for (i in 1:n_objects)
  {
    # Validation of the class
    if (!is(object = object[[i]], class2 = "msel"))
    {
      stop(paste0("Argument object[[", i, "]] is not of class 'msel'."))
    }
    
    # Validation of the covariance matrix type
    if (!(object[[i]]$cov_type %in% c("sandwich", "mm")))
    {
      stop(paste0("Covariance matrix type of object[[", i, "]] should be ",
                  "either 'sandwich' or 'mm', but '", object[[i]]$cov_type,
                  "' has been supplied."))
    }
    
    # Validation of the crossind variable
    if (!hasName(x = object[[i]]$data, "crossind"))
    {
      stop (paste0("There is no 'crossind' variable in data of object[[", 
                   i, "]]."))
    }
    
    # Indexes
    ind_start    <- n_par + 1
    n_par        <- n_par + object[[i]]$other$n_par
    ind_end      <- n_par
    ind_par[[i]] <- ind_start:ind_end
  }
  
  # Create indexes
  crossind <- numeric()
  for (i in 1:n_objects)
  {
    crossind <- c(crossind, object[[i]]$data$crossind)
  }
  crossind <- unique(crossind)
  
  # Total number of observations
  n <- max(crossind)
  
  # Store the Jacobians and Hessians
  J <- matrix(0, nrow = n,     ncol = n_par)
  H <- matrix(0, nrow = n_par, ncol = n_par)
  for (i in 1:n_objects)
  {
    J[object[[i]]$data$crossind, ind_par[[i]]] <- object[[i]]$J
    H[ind_par[[i]], ind_par[[i]]]              <- object[[i]]$H
  }
  
  # Remove bad values
  J_remove  <- apply(J, 1, function(x) any(is.na(x) | is.infinite(x)))
  if (any(J_remove))
  {
    J <- J[!J_remove, ]
    n <- nrow(J)
    warning("NA, NaN or Inf values in J")
  }

  # Check whether clustering is needed
  is_cluster <- FALSE
  for (i in 1:n_objects)
  {
    if ("crosscluster" %in% colnames(object[[i]]$data))
    {
      is_cluster <- TRUE
    }
    else
    {
      if (is_cluster)
      {
        stop (paste0("There is no 'crosscluster' variable in data of object[[",
                     i, "]]. If clustering of standard errors is not needed, ",
                     "then remove this variable from data of each model. ",
                     "Otherwise, add this variable to data of each model."))
      }
    }
  }
  
  # Stuff for clustered data
  if (is_cluster)
  {
    # Combine indexes of the clusters
    cluster_ind <- vector("numeric", length = n)
    for (i in 1:n_objects)
    {
      cluster_ind[object[[i]]$data$crossind] <- object[[i]]$data$crosscluster
    }
    cluster_unique <- unique(cluster_ind)
    n_cluster      <- length(cluster_unique)
    
    # Estimate clustered scores
    J_cluster  <- matrix(0, nrow = n_cluster, ncol = ncol(J))
    for (i in 1:n_cluster)
    {
      J_cluster[i, ] <- colSums(J[cluster_ind == cluster_unique[i], , 
                                  drop = FALSE])
    }
    
    # Rewrite initial scores with clustered scores
    J <- J_cluster
    n <- nrow(J)
  }
    
  # Remove rows of zeros
  J_zero_row <- apply(X = J, MARGIN = 1, FUN = function(x) all(x == 0))
  J          <- J[!J_zero_row, ]
  n          <- nrow(J)
  
  # Sandwich estimator
  cov <- tcrossprod(qr.solve(H, t(J), tol = 1e-16))
  
  # Names
  colnames(cov) <- 1:n_par
  rownames(cov) <- 1:n_par
  for (i in 1:n_objects)
  {
    cov_names                   <- paste0("m", i, "p", 
                                          1:object[[i]]$other$n_par)
    colnames(cov)[ind_par[[i]]] <- cov_names
    rownames(cov)[ind_par[[i]]] <- cov_names
  }
  
  # Extended return
  if (type == "all")
  {
    out <- list(cov = cov, J = J, H = H, n = n, n_par = n_par)
    return (out)
  }
  
  # Return the results
  return(cov)
}