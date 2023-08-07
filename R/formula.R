#' Formulas of mnprobit model.
#' @description Provides formulas associated with the 
#' object of class 'mnprobit'.
#' @param x object of class 'mnprobit'.
#' @param ... further arguments (currently ignored).
#' @param type character; 
#' if \code{type = "formula"} or \code{type = 1} then function returns 
#' a formulas of multinomial equation. 
#' If \code{type = "formula2"} or \code{type = 2} then function returns 
#' a formula of continuous equation.
#' @return Returns a formula.
formula.mnprobit <- function(x, ..., type = "formula") 
{
  if (length(list(...)) > 0)
  {
    warning("Additional arguments passed through ... are ignored.")   
  }
  
  # Validation
  if (type == 1)
  {
    type <- "formula"
  }
  if (type == 2)
  {
    type <- "formula2"
  }
  if (!(type %in% c("formula", "formula2")))
  {
    stop("Incorrect 'type' argument.")
  }
  
  # Formula of multinomial equation
  if (type == "formula")
  {
    return (x$formula)
  }
  
  # Formula of continuous equation
  if (type == "formula2")
  {
    if (!x$other$is2)
    {
      return (NULL)
    }
    return (x$formula2)
  }
  
  return (NULL)
}

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

#' Formulas of mvoprobit model.
#' @description Provides formulas associated with the 
#' object of class 'mvoprobit'.
#' @param x object of class 'mvoprobit'.
#' @param ... further arguments (currently ignored).
#' @param type character; 
#' if \code{type = "formula"} or \code{type = 1} then function returns formulas 
#' of ordered equations. 
#' If \code{type = "formula2"} or \code{type = 2} then function returns formulas 
#' of continuous equations.
#' @param eq positive integer representing the index of the equation which
#' formula should be returned. If \code{NULL} (default) then formulas for each
#' equation will be returned as a list which \code{i}-th element associated 
#' with \code{i}-th equation.
#' @return Returns a formula or a list of formulas depending on \code{eq} value.
formula.mvoprobit <- function(x, ..., type = "formula", eq = NULL) 
{
  if (length(list(...)) > 0)
  {
    warning("Additional arguments passed through ... are ignored.")   
  }
  
  # Validation
  if (type == 1)
  {
    type <- "formula"
  }
  if (type == 2)
  {
    type <- "formula2"
  }
  if (type == "formula")
  {
    if (!is.null(eq))
    {
      if ((eq <= 0) | (eq > x$control_lnL$n_eq))
      {
        stop("incorrect 'eq' value.")
      }
    }
  }
  if (type == "formula2")
  {
    if (!x$other$is2)
    {
      return (NULL)
    }
    if (!is.null(eq))
    {
      if ((eq <= 0) | (eq > x$control_lnL$n_eq2))
      {
        stop("incorrect 'eq' value.")
      }
    }
  }
  if (!(type %in% c("formula", "formula2")))
  {
    stop("Incorrect 'type' argument.")
  }
  
  # Formula of ordered equation
  if (type == "formula")
  {
    if (x$control_lnL$n_eq == 1)
    {
      eq <- 1
    }
    if (!is.null(eq))
    {
      return(x$formula[[eq]])
    }
    return(x$formula)
  }
  
  # Formula of continuous equation
  if (type == "formula2")
  {
    if (x$control_lnL$n_eq2 == 1)
    {
      eq <- 1
    }
    if (!is.null(eq))
    {
      return(x$formula2[[eq]])
    }
    return(x$formula2)
  }
  
  return (NULL)
}

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

#' Merge formulas
#' @description This function merges all variables of several formulas
#' into a single formula.
#' @param ... formulas to be merged such that there is a single element
#' on the left hand side and various elements on the right hand side.
#' @param type string representing the type of merge to be used.
#' If \code{type = "all"} then both right hand side and left hand side
#' elements of the formulas will be merged on the right hand side.
#' If \code{type = "terms"} then only right hand side elements of the
#' formulas will be merged on the right hand side.
#' If \code{type = "var-terms"} then the result is the same as in case
#' when \code{type = "terms"} but there will be left hand side element
#' of the first formula on the left hand side of the merged formula.
#' @details Merged formulas should have a single element on the left hand
#' side and voluntary number of elements on the right hand side.
#' @return This function returns a formula which form depends on 
#' \code{type} input argument value. See 'Details' for additional information.
#' @examples 
#' # Consider three formulas
#' f1 <- as.formula("y1 ~ x1 + x2")
#' f2 <- as.formula("y2 ~ x2 + x3")
#' f3 <- as.formula("y3 ~ y2 + x6")
#' # Merge these formulas in a various ways
#' formula_merge(f1, f2, f3, type = "all")
#' formula_merge(f1, f2, f3, type = "terms")
#' formula_merge(f1, f2, f3, type = "var-terms")
#' 
formula_merge <- function(..., type = "all")
{
  formulas <- list(...)
  if (is.list(formulas[[1]]))
  {
    formulas <- formulas[[1]]
  }
  f.y <- NULL
  f.x <- NULL
  for (i in 1:length(formulas))
  {
    f.x <- c(f.x, labels(terms(formulas[[i]])))
    f.y <- c(f.y, all.vars(formulas[[i]])[1])
  }
  
  if (type == "all")
  {
    f.c <- unique(c(f.x, f.y))
  }
  if (type %in% c("terms", "var-terms"))
  {
    f.c <- unique(f.x)
  }
  
  frm <- do.call(paste, c(as.list(f.c), sep = " + "))
  frm <- paste("~ ", frm)
  if(type == "var-terms")
  {
    frm <- paste(f.y[1], frm)
  }
  frm <- as.formula(frm)
  
  return(frm)
}

#' Split formula by symbol
#' @description This function splits one formula into two formulas
#' by symbol.
#' @param formula an object of class \code{formula}.
#' @param symbol a string that is used to split \code{formula} into
#' two formulas.
#' @details The \code{symbol} should be on the right hand side of
#' the formula. 
#' @return This function returns a list of two formulas.
#' @examples 
#' formula_split("y ~ x1 + x2 | x2 + x3")
#' formula_split("y ~ x1 + x2 : x2 + x3", symbol = ":")
#' 
formula_split <- function(formula, symbol = "|")
{
  formula <- as.formula(formula)
  fs <- paste(formula[2], formula[3], sep = '~')
  formula <- strsplit(x = fs, split = symbol, fixed = TRUE)
  
  f1y <- strsplit(x = formula[[1]][1], split = "~", fixed = TRUE) 
  y <- f1y[[1]][1]
  f1 <- f1y[[1]][2]
  f2 <- formula[[1]][2]
  
  f1 <- as.formula(paste(y, "~", f1))
  f2 <- as.formula(paste(y, "~", f2))
  
  f_out <- list(f1 = f1,
                f2 = f2)
  return(f_out)
}

