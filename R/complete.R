# Get the data without the missing observations
complete_msel <- function(object, data = NULL)
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
  if (is.null(data))
  {
    object <- object$data
  }
  
  # Convert the dependent variable to unobservable if
  # there is some unobservable independent variables
    # ordered equations
  if (is1)
  {
    frm_tmp <- NULL
    for (i in seq_len(n_eq))
    {
      if (is_het[i])
      {
        frm_tmp <- switchSelection::formula_merge(
          switchSelection::formula_split(formula[[i]]), 
          type = "var-terms")
      }
      else
      {
        frm_tmp <- formula[[i]]
      }
      df <- model.frame(frm_tmp, data, na.action = na.pass)
      data[!complete.cases(df), names(df)[1]] <- NA
    }
  }
    # continuous equations
  if (is2)
  {
    for (i in seq_len(n_eq2))
    {
      df2 <- model.frame(formula2[[i]], data, na.action = na.pass)
      data[!complete.cases(df2), names(df2)[1]] <- NA
    }
  }
    # multinomial equations
  if (is3)
  {
    df3 <- model.frame(formula3, data, na.action = na.pass)
    data[!complete.cases(df3), names(df3)[1]] <- NA
  }
  
  # Return the data
  return(data)
}