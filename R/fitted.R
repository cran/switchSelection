#' Extract Model Fitted Values
#' @description Extracts fitted values from 'msel' object
#' @param object object of class 'msel'.
#' @param ... further arguments (currently ignored).
#' @param newdata an optional data frame in which to look for variables 
#' with which to predict. If omitted, the original data frame used. 
#' This data frame should contain values of dependent variables even if they 
#' are not actually needed for prediction (simply assign them with 0 values).
#' @returns Returns a numeric matrix. Its columns which names coincide with the
#' names of the ordinal and multinomial equations provide the index of the most 
#' probable category for each observation.
#' Columns which names coincide with the names of the continuous equations 
#' provide conditional expectations of the dependent variables in observable 
#' regimes for each observation.
fitted.msel <- function (object, ..., newdata = NULL)
{
  # Validate dots
  if (length(list(...)) > 0)
  {
    warning("Additional arguments passed through ... are ignored.")   
  }
  
  # Some variables
  n_eq       <- object$other$n_eq
  n_eq2      <- object$other$n_eq2
  n_eq3      <- object$other$n_eq3
  is1        <- object$other$is1
  is2        <- object$other$is2
  is3        <- object$other$is3
  n_regimes  <- object$other$n_regimes
  z_names    <- object$other$z_names
  y_names    <- object$other$y_names
  z_mn_names <- object$other$z_mn_names
  if (is.null(newdata))
  {
    newdata <- object$data
  }
  n          <- nrow(newdata)
  n_cuts_eq  <- object$other$n_cuts_eq
  
  # Output variable
  out           <- matrix(NA, nrow = n, ncol = n_eq + n_eq2 + is3)
  colnames(out) <- rep("val", ncol(out))
  rownames(out) <- 1:n
  
  # Most probable categories of the ordinal equations
  if (is1)
  {
    group <- rep(-1, n_eq)
    for (i in 1:n_eq)
    {
      probs <- matrix(NA, nrow = n, ncol = n_cuts_eq[i] + 1)
      for (j in 0:n_cuts_eq[i])
      {
        group[i]         <- j
        probs[, j + 1]   <- predict(object, type = "prob", newdata = newdata, 
                                    group = group)
        group[i]         <- -1
      }
      probs_obs               <- complete.cases(probs)
      probs_tmp               <- probs
      probs_tmp[!probs_obs, ] <- 0
      out[, i]                <- apply(X = probs_tmp, MARGIN = 1, 
                                       FUN = which.max) - 1
      out[!probs_obs, i]      <- NA
    }
    colnames(out)[1:n_eq] <- z_names
  }
  
  # Most probable categories of the multinomial equation
  if (is3)
  {
    probs <- matrix(NA, nrow = n, ncol = n_eq3)
    for (j in 1:n_eq3)
    {
      probs[, j] <- predict(object, type = "prob_mn", 
                            newdata = newdata, group3 = j - 1)
    }
    out[, ncol(out)]         <- apply(X = probs, MARGIN = 1, 
                                      FUN = which.max) - 1
    colnames(out)[ncol(out)] <- z_mn_names
  }
  
  # Fitted values of the continuous equations
  if (is2)
  {
    ind                <- (n_eq + 1):(n_eq + n_eq2)
    out[, ind]         <- predict(object, type = "val", newdata = newdata)
    colnames(out)[ind] <- y_names
  }
  
  # Return the result
  return(out)
}