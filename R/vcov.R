# Estimate asymptotic covariance matrix of the maximum-likelihood estimator
vcov_ml <- function(object, type, n_cores, n_sim)
{
  # List of output values
  out <- list(vcov = NULL)
  
  # Get some variables
  par            <- object$par
  cov_type       <- type
  estimator      <- object$estimator
  control_lnL    <- object$control_lnL
  n_par          <- object$other$n_par
  regularization <- object$other$regularization

  # Estimate asymptotic covariance matrix
  H     <- NULL
  H_inv <- NULL
  J     <- NULL
  cov   <- diag(rep(1, n_par))
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
                                           gr.args = gr.args)
            out$H    <- H
            H_inv    <- qr.solve(H, tol = 1e-16)
            cov_type <- cov_type_old
          },
          error = function(e) {
            warning(paste0("Problems with numeric hessian calculation. ",
                           "Therefore 'cov_type' has been changed to 'gop'."))
          }
        )
      }
      if (cov_type %in% c("sandwich", "gop", "mm"))
      {
        J <- NULL
        if (cov_type == "mm")
        {
          J <- lnL_msel(par         = par,
                        n_sim       = n_sim, 
                        n_cores     = n_cores,
                        control_lnL = control_lnL, 
                        out_type    = "jac",
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
      }
      if ((cov_type == "sandwich") | (cov_type == "mm"))
      {
        cov <- H_inv %*% t(J) %*% J %*% H_inv
      }
      if (cov_type == "hessian")
      {
        cov <- -H_inv
      }
      
      if (any(is.na(cov)))
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
          cov <- qr.solve(t(J) %*% J, tol = 1e-16)
        },
        error = function(e) {
          warning(paste0("Problems with numeric hessian calculation. ",
                         "Therefore 'cov_type' has been changed to 'NO'."))
        }
       )
      }
    }
  } 
  
  # Aggregate some output
  out$vcov <- cov
  
  return(out)
}

# Estimate the asymptotic covariance matrix of the two-step estimator
vcov_2step <- function(object, type = "default", n_cores = 1, n_sim = 10000)
{
  # Some variables
  n_eq    <- object$other$n_eq
  n_eq2   <- object$other$n_eq2
  n_eq3   <- object$other$n_eq3
  n_obs   <- object$other$n_obs
  is1     <- object$other$is1
  is3     <- object$other$is3
  type3   <- object$type3
  n_sigma <- object$other$n_sigma
  
  # Calculate the scores of the second step
  scores <- scores_2step(object, n_cores = n_cores, n_sim = n_sim)
  
  # Calculate the derivatives respect to the scores of the second step
  scores_list <- deriv_msel(object  = object, 
                            fn      = scores_2step, 
                            fn_args = list(type    = "aggregate",
                                           n_cores = n_cores, 
                                           n_sim   = n_sim))
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
    marginal_par_ind  <- object$ind$marginal_par_ind
    sigma_ind         <- object$ind$sigma
    coef3_ind         <- object$ind$coef3
    sigma3_ind        <- object$ind$sigma3
    
    # Indexes of the first step in first step model
    coef_ind1         <- object$model1$ind$coef
    cuts_ind1         <- object$model1$ind$cuts
    coef_var_ind1     <- object$model1$ind$coef_var
    marginal_par_ind1 <- object$model1$ind$marginal_par_ind
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
  }
  
  # Estimate asymptotic covariance matrix
  scores_cov      <- cov(scores)
  vcov            <- qr.solve(t(scores_jac) %*% 
                              qr.solve(scores_cov, tol = 1e-16) %*% 
                              scores_jac, tol = 1e-16) * n_obs
  colnames(vcov)  <- NULL
  rownames(vcov)  <- NULL
  
  # Return the results
  return(list(scores_cov  = scores_cov, 
              scores_jac  = scores_jac,
              scores      = scores,
              scores_sum  = scores_sum,
              vcov        = vcov))
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
#' parallel computing. If possible it is highly recommend to set it equal to
#' the number of available physical cores especially when the system of
#' ordered equations has 2 or 3 equations.
#' @param recalculate logical; if \code{TRUE} then covariance matrix will be
#' recalculated even if 'type' is the same as 'cov_type' input argument
#' of the model.
#' @details Argument \code{type} is closely related to the argument 
#' \code{cov_type} of \code{\link[switchSelection]{msel}} function. 
#' See 'Details' and 'Usage' sections of \code{\link[switchSelection]{msel}} 
#' for more information on \code{cov_type} argument.
#' @return Returns numeric matrix which represents estimate of the asymptotic 
#' covariance matrix of model's parameters.
vcov.msel <- function(object, ..., 
                      type        = object$cov_type, 
                      n_cores     = object$other$n_cores,
                      n_sim       = object$other$n_sim,
                      recalculate = FALSE)
{
  # Validate dots
  if (length(list(...)) > 0)
  {
    warning("Additional arguments passed through ... are ignored.")   
  }
  
  if (object$estimator == "ml")
  {
    if ((type == object$cov_type) & !recalculate)
    {
      return (object$cov)
    }
    return (vcov_ml(object = object, type = type, 
                    n_cores = n_cores, n_sim = n_sim)$vcov)
  }
  
  if (object$estimator == "2step")
  {
    if (type != object$cov_type)
    {
      stop(paste0("If estimator is '2step' then argument type should be the ",
                  "same as 'object$cov_type'."))
    }
    return(object$cov)
  }
  
  stop("Incorrect 'object$estimator' value.")
}