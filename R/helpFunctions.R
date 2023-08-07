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