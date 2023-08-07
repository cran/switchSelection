#' Coefficients extraction method for mvoprobit.
#' @description Extract coefficients and other estimates from mvoprobit object.
#' @param object object of class "mvoprobit".
#' @param ... further arguments (currently ignored).
#' @param eq integer representing an index of the ordered equation.
#' @param eq2 integer representing an index of the continuous equation.
#' @param regime integer representing a regime of the continuous equation.
#' @param type character representing a type of the output. Possible options
#' are \code{"coef"}, \code{"coef2"}, \code{"cov"}, \code{"cov1"}, \code{"var"},
#' \code{"cov2"}, \code{"cov3"}, \code{coef_lambda} and \code{marginal}.
#' See 'Details' for additional information.
#' @template coef_mvoprobit_details_Template
#' @returns See 'Details' section.
coef.mvoprobit <- function (object, ..., eq = NULL, eq2 = NULL, 
                            regime = NULL, type = "coef")
{
  # Validate dots
  if (length(list(...)) > 0)
  {
    warning("Additional arguments passed through ... are ignored.")   
  }
  
  # Some variables
  n_eq <- length(object$coef)
  n_eq2 <- length(object$coef2)
  is2 <- hasName(object, "formula2")
  n_regimes <- object$other$n_regimes
  z_names = object$other$z_names
  y_names = object$other$y_names
  
  # Validation of type
  type <- tolower(type)
  if ((type == "coefficients") | (type == "coeff"))
  {
    type <- "coef"
    warning("It is assumed that 'type' is 'coef'.")
  }
  if ((type == "coefficients2") | (type == "coeff2"))
  {
    type <- "coef2"
    warning("It is assumed that 'type' is 'coef2'.")
  }
  if ((type == "coef") & !is.null(eq2))
  {
    type <- "coef2"
    warning("Since 'eq2' is not 'NULL' it is assumed that 'type' is 'coef2'.")
  }
  if ((type == "lambda") | (type == "coefficients_lambda") | 
      (type == "coeff_lambda") | (type == "lambdas"))
  {
    type <- "coef_lambda"
    warning("It is assumed that 'type' is 'coef_lambda'.")
  }
  if ((type == "coefficients_var") | (type == "coeff_var"))
  {
    type <- "coef_var"
    warning("It is assumed that 'type' is 'coef_var'.")
  }
  if (type == "marginal_par")
  {
    type <- "marginal"
    warning("It is assumed that 'type' is 'marginal'.")
  }
  
  # Convert eq into numeric if need
  if (!is.null(eq))
  {
    if (is.character(eq))
    {
      eq <- which(z_names %in% eq)
    }
    if (any(eq < 1) | any(eq > n_eq))
    {
      stop("Incorrect 'eq' value.")
    }
  }
  
  # Convert eq2 into numeric if need
  if (is2)
  {
    if (!is.null(eq2))
    {
      if (!is.numeric(eq2))
      {
        eq2 <- which(y_names %in% eq2)
      }
      if (any(eq2 < 0) | any(eq2 > n_eq2))
      {
        stop("Incorrect 'eq2' value.")
      }
    }
  }
  
  # Validate regime
  if (is2)
  {
    if (!is.null(regime) & !is.null(eq2))
    {
      if (any(regime < 0) | any(regime >= n_regimes[eq2]))
      {
        stop("Incorrect 'regime' value.")
      }
    }
  }
  
  # Coefficients of the ordered equations
  if (type == "coef")
  {
    if (is.null(eq))
    {
      return (object$coef)
    }
    return (object$coef[[eq]])
  }
  
  # Coefficients of continuous equations
  if (type == "coef2")
  {
    if (!is2)
    {
      stop("Available only when there is at least one continuous equation.")
    }
    if (is.null(eq2) & !is.null(regime))
    {
      eq2 <- 1
    }
    if (is.null(eq2))
    {
      return (object$coef2)
    }
    if (is.null(regime))
    {
      return (object$coef2[[eq2]])
    }
    return (object$coef2[[eq2]][regime + 1, ])
  }
  
  # Coefficients of the variance equation
  if (type == "coef_var")
  {
    if (is.null(eq))
    {
      return (object$coef_var)
    }
    return (object$coef_var[[eq]])
  }
  
  # Coefficients of selectivity correction terms (lambda)
  if (type == "coef_lambda")
  {
    if (object$estimator != "2step")
    {
      stop("Available only for two-step estimator.")
    }
    if ((n_eq2 >= 2) | !is2)
    {
      stop("Available only when there is a single continuous equation.")
    }
    if (is.null(regime))
    {
      return (object$coef_lambda)
    }
    if (is.null(eq))
    {
      return (object$coef_lambda[[regime + 1]])
    }
    return (object$coef_lambda[[regime + 1]][[eq]])
  }
  
  # Covariances
  # asymptotic covariance matrix of estimates
  if ((type == "vcov") | (type == "cov"))
  {
    return (object$cov)
  }
  # covariances between ordered equations
  if (type == "cov1")
  {
    if (length(eq) > 1)
    {
      return(object$sigma[eq[1], eq[2]])
    }
    return(object$sigma)
  }
  # covariances between continuous and ordered equations
  if (type == "cov12")
  {
    if (any(object$degrees != 1))
    {
      stop("Not available if some 'degrees' are different from 1.")
    }
    if (!is.null(eq) & !is.null(regime) & is.null(eq2))
    {
      eq2 <- 1
    }
    if (is.null(eq) | is.null(eq2))
    {
      return (object$cov2)
    }
    if (is.null(regime))
    {
      return (object$cov2[[eq2]])
    }
    return (object$cov2[[eq2]][regime + 1, eq])
  }
  # variances of continuous equations
  if ((type == "var") | (type == "cov2") | (type == "variance"))
  {
    if (!is2)
    {
      stop("Available only when there is a continous equation.")
    }
    if (any(object$degrees != 1))
    {
      stop("Not available if some 'degrees' are different from 1.")
    }
    if (is.null(eq2) & !is.null(regime))
    {
      eq2 <- 1
    }
    if (is.null(eq2))
    {
      if (n_eq2 == 1)
      {
        return (object$var2[[1]])
      }
      return (object$var2)
    }
    if (is.null(regime))
    {
      return (object$var2[[eq2]])
    }
    if ((regime < 0) | (regime >= n_regimes[eq2]))
    {
      stop("Incorrect 'regime' value.")
    }
    return(object$var2[[eq2]][regime[1] + 1])
  }
  # covariances between continuous equations
  if (type == "cov3")
  {
    if (n_eq2 < 2)
    {
      stop("Available only when there are more than one continuous equations.")
    }
    if (is.null(eq2))
    {
      return (object$sigma2)
    }
    if ((length(regime) != 2) | (length(eq2) != 2))
    {
      stop("Arguments 'regime' and 'eq2' should be of length '2'.")
    }
    if (eq2[1] == eq2[2])
    {
      stop("Indexes of continuous equations 'eq2' should be distinct.")
    }
    counter <- 1
    for (i in 2:n_eq2)
    {
      for (j in 1:(i - 1))
      {
        if (((i == eq2[1]) & (j == eq2[2])) | ((i == eq2[2]) & (j == eq2[1])))
        {
          pairs <- object$other$regimes_pair[[counter]]
          ind <- which(((pairs[, 1] == regime[1]) & (pairs[, 2] == regime[2])) |
                         ((pairs[, 1] == regime[2]) & (pairs[, 2] == regime[1])))
          if (length(ind) > 0)
          {
            return (object$sigma2[[counter]][ind])
          }
          return (NA)
        }
        counter <- counter + 1
      }
    }
  }
  
  # Parameters of marginal distribution
  if (type == "marginal")
  {
    if (is.null(eq))
    {
      return(object$marginal_par)
    }
    return(object$marginal_par[[eq]])
  }
  
  stop("Incorrect 'type' argument.")
}

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

#' Coefficients extraction method for mnprobit.
#' @description Extract coefficients and other estimates from mnprobit object.
#' @param object object of class "mnprobit"
#' @param ... further arguments (currently ignored)
#' @param alt integer representing index of the alternative
#' @param regime integer representing regime of the continuous equation
#' @param type character representing the type of the output. Possible options
#' are \code{"coef"}, \code{"coef2"}, \code{"cov1"}, \code{"var"},
#' \code{"cov2"}, \code{coef_lambda}.
#' See 'Details' for additional information.
#' @template coef_mnprobit_details_Template
#' @returns See 'Details' section.
coef.mnprobit <- function (object, ..., alt = NULL, 
                           regime = NULL, type = "coef")
{
  # Validate dots
  if (length(list(...)) > 0)
  {
    warning("Additional arguments passed through ... are ignored.")   
  }
  
  # Some variables
  n_alt <- object$n_alt
  n_regimes <- object$control_lnL$n_regimes
  is2 <- hasName(object, "formula2")
  
  # Validation of type
  type <- tolower(type)
  if ((type == "coefficients") | (type == "coeff"))
  {
    type <- "coef"
    warning("It is assumed that 'type' is 'coef'.")
  }
  if ((type == "coefficients2") | (type == "coeff2"))
  {
    type <- "coef2"
    warning("It is assumed that 'type' is 'coef2'.")
  }
  if ((type == "coef") & !is.null(regime))
  {
    type <- "coef2"
    warning("Since 'regime' is not 'NULL' it is assumed that 'type' is 'coef2'.")
  }
  if ((type == "lambda") | (type == "coefficients_lambda") | 
      (type == "coeff_lambda") | (type == "lambdas"))
  {
    type <- "coef_lambda"
    warning("It is assumed that 'type' is 'coef_lambda'.")
  }
  
  # Validate regime
  if (is2)
  {
    if (!is.null(regime))
    {
      if ((regime < 0) | (regime >= n_regimes))
      {
        stop("Incorrect 'regime' value.")
      }
    }
  }
  
  # Validate alternative
  if (!is.null(alt))
  {
    if (any(alt == n_alt))
    {
      stop(paste0("For the alternative ", n_alt, " all coefficients are ",
                  "constrained to 0 for identification purposes."))
    }
    if (any(alt < 1) | any(alt >= n_alt))
    {
      stop("Incorrect 'alt' value.")
    }
  }
  
  # Coefficients of the multinomial equations
  if (type == "coef")
  {
    if (is.null(alt))
    {
      return (object$coef)
    }
    return (object$coef[, alt])
  }
  
  # Coefficients of the continuous equations
  if (type == "coef2")
  {
    if (is.null(regime))
    {
      return (object$coef2)
    }
    return (object$coef2[, regime + 1])
  }
  
  # Coefficients of selectivity correction terms (lambda)
  if (type == "coef_lambda")
  {
    if (object$estimator != "2step")
    {
      stop("Available only for two-step estimator.")
    }
    if (is.null(regime))
    {
      return (object$coef_lambda)
    }
    if (is.null(alt))
    {
      return (object$coef_lambda[[regime + 1]])
    }
    return (object$coef_lambda[[regime + 1]][[alt]])
  }
  
  # Estimate of the covariance matrix of
  # random errors of multinomial equations
  if (type == "cov1")
  {
    if (length(alt) == 1)
    {
      return (object$sigma)
    }
    return (object$sigma[alt[1], alt[2]])
  }
  
  # Estimates of the variances of
  # random errors of continuous equation
  # in different regimes
  if ((type == "cov2") | (type == "var"))
  {
    if (!is2)
    {
      stop("Available only when there is a continous equation.")
    }
    if (object$estimator != "ml")
    {
      stop("Available only for maximum-likelihood estimator.")
    }
    if (is.null(regime))
    {
      return (object$var2)
    }
    return (object$var2[regime + 1])
  }
  
  # Covariances between random errors of
  # multinomial and continuous equations
  if (type == "cov12")
  {
    if (!is2)
    {
      stop("Available only when there is a continous equation.")
    }
    if (object$estimator != "ml")
    {
      stop("Available only for maximum-likelihood estimator.")
    }
    if (is.null(regime) & is.null(alt))
    {
      return (object$cov2)
    }
    if (is.null(regime) & !is.null(alt))
    {
      return (object$cov2[alt, ])
    }
    if (!is.null(regime) & is.null(alt))
    {
      return (object$cov2[, regime + 1])
    }
    return (object$cov2[alt, regime + 1])
  }
  
  # Asymptotic covariance matrix of the estimator
  if ((type == "cov") | (type == "vcov"))
  {
    return(object$cov)
  }
  
  stop("Incorrect 'type' argument.")
}