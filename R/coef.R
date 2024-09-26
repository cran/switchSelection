#' Coefficients extraction method for msel.
#' @description Extract coefficients and other estimates from msel object.
#' @param object an object of class "msel".
#' @param ... further arguments (currently ignored).
#' @param eq an integer representing the index of the ordered equation.
#' @param eq2 an integer representing the index of the continuous equation.
#' @param eq3 an integer representing the index of the alternative of the
#' multinomial equation.
#' @param regime an integer representing a regime of the continuous equation.
#' @param type a character representing a type of the output. Possible options
#' are \code{"coef"}, \code{"coef2"}, \code{coef_lambda}, 
#' \code{"coef_var"}, \code{"coef3"}, \code{"cuts"}, \code{"cov"}, 
#' \code{"cov1"}, \code{"var"}, \code{"cov2"}, \code{"cov3"}, 
#' and \code{marginal}.
#' See 'Details' for additional information.
#' @template coef_msel_details_Template
#' @returns See 'Details' section.
coef.msel <- function (object, ..., 
                       eq = NULL, eq2 = NULL, eq3 = NULL,
                       regime = NULL, type = "coef")
{
  # Validate dots
  if (length(list(...)) > 0)
  {
    warning("Additional arguments passed through ... are ignored.")   
  }
  
  # Some variables
  n_eq            <- object$other$n_eq
  n_eq2           <- object$other$n_eq2
  n_eq3           <- object$other$n_eq3
  is1             <- object$other$is1
  is2             <- object$other$is2
  is3             <- object$other$is1
  n_regimes       <- object$other$n_regimes
  z_names         <- object$other$z_names
  y_names         <- object$other$y_names
  estimator       <- object$estimator
  coef_lambda_ind <- object$other$coef_lambda_ind
  
  # Validation of the type
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
  if ((type == "lambda") | (type == "coef2_lambda"))
  {
    type <- "coef_lambda"
    warning("It is assumed that 'type' is 'coef_lambda'.")
  }
  if ((type == "coefficients3") | (type == "coeff3") |
      (type == "coef_mn")       | (type == "coef_mult"))
  {
    type <- "coef3"
    warning("It is assumed that 'type' is 'coef3'.")
  }
  if ((type == "coef") & !is.null(eq2))
  {
    type <- "coef2"
    warning("Since 'eq2' is not 'NULL' it is assumed that 'type' is 'coef2'.")
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
  if (type %in% c("cut", "threshold", "thresholds"))
  {
    type <- "cuts"
    warning("It is assumed that 'type' is 'cuts'.")
  }
  if (type %in% c("sigma", "sigma1"))
  {
    type <- "cov1"
    warning("It is assumed that 'type' is 'cov1'.")
  }
  if (type %in% c("sigma2"))
  {
    type <- "cov2"
    warning("It is assumed that 'type' is 'cov2'.")
  }
  if (type %in% c("sigma3"))
  {
    type <- "cov3"
    warning("It is assumed that 'type' is 'cov3'.")
  }
  
  # Convert eq into numeric if need
  if (is1)
  {
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
  
  # Convert eq3 into numeric if need
  if (is3)
  {
    if (!is.null(eq3))
    {
      if (any(eq3 < 0) | any(eq3 >= n_eq3))
      {
        stop("Incorrect 'eq3' value.")
      }
    }
  }
  
  # Validate the regime
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
  
  # Coefficients of the ordinal equations
  if (type == "coef")
  {
    if (is.null(eq))
    {
      return (object$coef)
    }
    return (object$coef[[eq]])
  }
  
  # Coefficients of the continuous equations
  if (type %in% c("coef2", "coef_lambda", "coef2_all"))
  {
    if (!is2)
    {
      stop("Available only when there is at least one continuous equation.")
    }
    if (estimator == "2step")
    {
      if (type != "coef2_all")
      {
        for (v in 1:n_eq2)
        {
          if (length(coef_lambda_ind[[v]]) > 0)
          {
            s                 <- ifelse(type == "coef_lambda", 1, -1)
            object$coef2[[v]] <- object$coef2[[v]][, s * coef_lambda_ind[[v]],
                                                     drop = FALSE]
          }
        }
      }
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
  
  # Cuts
  if (type == "cuts")
  {
    if (is.null(eq))
    {
      return (object$cuts)
    }
    return (object$cuts[[eq]])
  }
  
  # Asymptotic covariance matrix of the estimates
  if ((type == "vcov") | (type == "cov"))
  {
    return (object$cov)
  }
  
  # Covariances between the ordinal equations
  if (type == "cov1")
  {
    if (length(eq) > 1)
    {
      return(object$sigma[eq[1], eq[2], drop = FALSE])
    }
    return(object$sigma)
  }
  
  # Covariances between the continuous and ordinal equations
  if (type == "cov12")
  {
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
    return (object$cov2[[eq2]][regime + 1, eq, drop = FALSE])
  }
  
  # Variances of the continuous equations
  if ((type == "var") | (type == "variance"))
  {
    if (!is2)
    {
      stop("Available only when there is a continous equation.\n")
    }
    if (estimator != "ml")
    {
      stop("Available only when 'estimator = ml'.\n")
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
  
  # Covariances between the continuous equations
  if (type == "cov2")
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
      stop("Indexes of the continuous equations 'eq2' should be distinct.")
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
  
  # Parameters of the marginal distribution
  if (type == "marginal")
  {
    if (is.null(eq))
    {
      return(object$marginal_par)
    }
    return(object$marginal_par[[eq]])
  }
  
  # Coefficients of the multinomial equation
  if (type == "coef3")
  {
    if (is.null(eq3))
    {
      return (object$coef3)
    }
    return (object$coef3[eq3 + 1, ])
  }
  
  # Estimate of the covariance matrix of the 
  # random errors of the multinomial equations
  if (type == "cov3")
  {
    if (length(eq3) == 1)
    {
      return (object$sigma3)
    }
    return (object$sigma3[eq3[1] + 1, eq3[2] + 1, drop = FALSE])
  }
  
  stop("Incorrect 'type' argument.")
}