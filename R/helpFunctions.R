#' Stars for p-values
#' @description This function assigns stars 
#' (associated with different significance levels) to p-values.
#' @param p_value vector of values between 0 and 1 representing p-values.
#' @details Three stars are assigned to p-values not greater than 0.01.
#' Two stars are assigned to p-values greater than 0.01 and not greater 
#' than 0.05.
#' One star is assigned to p-values greater than 0.05 and not greater than 0.1.
#' @return The function returns a string vector of stars assigned according
#' to the rules described in 'Details' section.
#' @examples 
#' p_value <- c(0.002, 0.2, 0.03, 0.08)
#' starsVector(p_value)
starsVector <- function(p_value)
{
  n <- length(p_value)
  stars <- rep(NA, n)
  
  for (i in 1:n)
  {
    if (p_value[i] > 0.1)
    {
      stars[i] <- ""
      next
    }
    
    if ((p_value[i] <= 0.1) & (p_value[i] > 0.05))
    {
      stars[i] <- "*"
      next
    }
    
    if ((p_value[i] <= 0.05) & (p_value[i] > 0.01))
    {
      stars[i] <- "**"
      next
    }
    if (p_value[i] <= 0.01)
    {
      stars[i] <- "***"
      next
    }
  }
  
  return(stars)
}

remove_column <- function(mat, name, action = TRUE)
{
  if (!action)
  {
    return (mat)
  }
  mat <- mat[, -which(colnames(mat) == name), drop = FALSE]
}

#' Modify exogenous variables in data frame
#' @description Change some values of the exogenous variables in a data frame.
#' @param exogenous list such that \code{exogenous[[i]]} represents the value 
#' (or a vector of values of the same size as \code{nrow(newdata)}) which will 
#' be exogenously assigned to the variable \code{names(exogenous)[[i]]} in 
#' \code{newdata} i.e., \code{newdata[, names(exogenous)[i]] <- exogenous[[i]]}.
#' If \code{newdata} is \code{NULL} and \code{exogenous} is not \code{NULL} then
#' \code{newdata} is set to \code{object$data}.
#' This argument is especially useful for causal inference when some endogenous 
#' (dependent) variables should be exogenously assigned with some values i.e.,
#' in the right hand side of the \code{formula} and \code{formula2}.
#' The purpose of the \code{exogeneous} argument is just a \code{convenience} so
#' equivalently it is possible to exogenously provide the values to variables
#' via \code{newdata} argument.
#' @param newdata data frame.
#' @details This function changes \code{exogenous} variables in \code{newdata}.
#' @return The function returns data frame which is similar to \code{newdata} 
#' but some values of this data frame are set according to \code{exogenous}.
exogenous_fn <- function(exogenous, newdata)
{
  n_exogenous <- length(exogenous)
  exogenous_names <- names(exogenous)
  if (!is.list(exogenous))
  {
    stop ("Argument 'exogenous' should be a list.")
  }
  n_obs <- nrow(newdata)
  for (i in 1:n_exogenous)
  {
    n_exogenous_i <- length(exogenous[[i]])
    if (!is.numeric(exogenous[[i]]) | ((n_exogenous_i > 1) & 
                                       (n_exogenous_i != n_obs)))
    {
      stop (paste0("Elements of 'exogenous' argument should be ",
                   "numeric values of appropriate size but ",
                   "exogenous[[", i, "]] is not."))
    }
    if (!(exogenous_names[i] %in% colnames(newdata)))
    {
      stop (paste0("Incorrect value of exogenous[[", i, "]] since there ",
                   "is no variable ", exogenous_names[i], " in the data."))
    }
    newdata[, exogenous_names[i]] <- exogenous[[i]]
  }
  
  return (as.data.frame(newdata))
}