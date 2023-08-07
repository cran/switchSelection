#' Extract Model Fitted Values
#' @description Extracts fitted values from 'mvoprobit' object
#' @param object object of class 'mvoprobit'.
#' @param ... further arguments (currently ignored).
#' @param newdata an optional data frame in which to look for variables 
#' with which to predict. If omitted, the original data frame used. 
#' This data frame should contain values of dependent variables even if they 
#' are not actually needed for prediction (simply assign them with 0 values).
#' @returns Returns a data frame. Its columns which names coincide with the
#' names of the ordered equations provide an index of the most probable
#' category. Columns which names coincide with the
#' names of the continuous equations provide uncinditional expectations of
#' the dependent variables in available regimes.
fitted.mvoprobit <- function (object, ..., newdata = NULL)
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
  if (is.null(newdata))
  {
    newdata <- object$data
  }
  n <- nrow(newdata)
  
  # Output variable
  out <- data.frame(matrix(NA, nrow = n, ncol = n_eq + sum(n_regimes)))
  
  # Fitted values of the ordered equations
  for (i in 1:n_eq)
  {
    max_i <- max(object$groups[, i])
    probs_i <- matrix(NA, nrow = n, ncol = max_i + 1)
    for (j in 0:max_i)
    {
      group <- rep(-1, n_eq)
      group[i] <- j
      probs_i[, j + 1] <- predict(object, group = group, 
                                  type = "prob", newdata = newdata)
    }
    out[, i] <- apply(X = probs_i, MARGIN = 1, FUN = which.max) - 1
  }
  colnames(out)[1:n_eq] <- z_names
  
  # Fitted values of the continuous equations
  counter <- n_eq + 1
  for (i in seq_len(n_eq2))
  {
    for (j in 1:n_regimes[i])
    {
      group2 <- rep(-1, n_eq2)
      group2[i] <- j - 1
      out[, counter] <- predict(object, type = "val", 
                                group2 = group2, newdata = newdata)
      colnames(out)[counter] <- paste0(y_names[i], j - 1)
      counter <- counter + 1
    }
  }
  
  # Return the result
  return(out)
}

#' Extract Model Fitted Values
#' @description Extracts fitted values from 'mnprobit' object
#' @param object object of class 'mnprobit'.
#' @param ... further arguments (currently ignored).
#' @param newdata an optional data frame in which to look for variables 
#' with which to predict. If omitted, the original data frame used. 
#' This data frame should contain values of dependent variables even if they 
#' are not actually needed for prediction (simply assign them with 0 values).
#' @returns Returns a data frame. Its first column provides an index of the
#' most probable alternative. Columns which names coincide with the
#' names of the continuous equation provide unconditional expectation of
#' the dependent variable in available regimes.
fitted.mnprobit <- function (object, ..., newdata = NULL)
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
  if (is.null(newdata))
  {
    newdata <- object$data
  }
  n <- nrow(newdata)
  
  # Output variable
  out <- data.frame(matrix(NA, nrow = n, ncol = 1 + n_regimes))
  
  # Most probable alternative
  probs <- matrix(NA, nrow = n, ncol = n_alt)
  for (i in 1:n_alt)
  {
    probs[, i] <- predict(object, alt = i, type = "prob", newdata = newdata)
  }
  out[, 1] <- apply(X = probs, MARGIN = 1, FUN = which.max)
  colnames(out)[1] <- "Alternative"
  
  # Fitted values of the continuous equation
  for (i in seq_len(n_regimes))
  {
    out[, i + 1] <- predict(object, alt = NULL, regime = i - 1, 
                            type = "val", newdata = newdata)
    colnames(out)[i + 1] <- paste0("Regime ", i - 1)
  }
  
  # Return the result
  return(out)
}