#' Formulas of msel model.
#' @description Provides formulas associated with the 
#' object of class 'msel'.
#' @param x object of class 'msel'.
#' @param ... further arguments (currently ignored).
#' @param type character; 
#' if \code{type = "formula"} or \code{type = 1} then function returns formulas 
#' of the ordinal equations. 
#' If \code{type = "formula2"} or \code{type = 2} then function returns formulas 
#' of the continuous equations.
#' If \code{type = "formula3"} or \code{type = 3} then function returns formula 
#' of the multinomial equation.
#' @param eq positive integer representing the index of the equation which
#' formula should be returned. If \code{NULL} (default) then formulas for each
#' equation will be returned as a list which \code{i}-th element associated 
#' with \code{i}-th equation.
#' @return Returns a formula or a list of formulas depending on \code{eq} value.
formula.msel <- function(x, ..., type = "formula", eq = NULL) 
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
  if (type == 3)
  {
    type <- "formula3"
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
  if (type == "formula3")
  {
    if (!x$other$is3)
    {
      return (NULL)
    }
  }
  if (!(type %in% c("formula", "formula2", "formula3")))
  {
    stop("Incorrect 'type' argument.")
  }
  
  # Formula of the ordinal equation
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
  
  # Formula of the continuous equation
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
  
  # Formula of the multinomial equation
  if (type == "formula3")
  {
    return(x$formula3)
  }
  
  return (NULL)
}

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

# Create formulas
formula_msel <- function(object, formula, formula2, formula3, 
                         degrees = NULL, degrees3 = NULL)
{
  # Get some variables
  is1        <- object$other$is1
  is2        <- object$other$is2
  is3        <- object$other$is3
  estimator  <- object$estimator
  n_eq       <- object$other$n_eq
  n_eq2      <- object$other$n_eq2
  n_eq3      <- object$other$n_eq3
  
  # Include lambda into formula2 according to degrees
  if ((estimator == "2step") & is1 & (!is.null(degrees)))
  {
    for (v in seq_len(n_eq2))
    {
      for (j in seq_len(n_eq))
      {
        if (degrees[v, j, drop = FALSE] != 0)
        {
          for (j1 in seq_len(degrees[v, j]))
          {
            # non-interaction terms
            if (j1 == 1)
            {
              formula2[[v]] <- update(formula2[[v]], paste0("~. + lambda", j))
            }
            else
            {
              formula2[[v]] <- update(formula2[[v]], 
                                      paste0("~. + I(lambda", j, "^", j1, ")"))
            }
          }
        }
      }
    }
  }
  
  # Include lambda_mn into formula3 according to degrees3
  if ((estimator == "2step") & is3 & (!is.null(degrees3)))
  {
    for (v in seq_len(n_eq2))
    {
      for (j in seq_len(ncol(degrees3)))
      {
        if (degrees3[v, j, drop = FALSE] != 0)
        {
          for (j1 in seq_len(degrees3[v, j]))
          {
            # non-interaction terms
            if (j1 == 1)
            {
              formula2[[v]] <- update(formula2[[v]], 
                                      paste0("~. + lambda", j, "_mn"))
            }
            else
            {
              formula2[[v]] <- update(formula2[[v]], 
                                      paste0("~. + I(lambda", j, "_mn", 
                                             "^", j1, ")"))
            }
          }
        }
      }
    }
  }
  
  # Coerce formula to type 'formula' and check whether it has "|" symbol
  # which indicates the presence of heteroscedasticity.
  # If there is heteroscedasticity then store separate formulas for
  # mean and variance equations.
  is_het       <- rep(FALSE, n_eq)
  formula_mean <- vector(mode = "list", length = n_eq)
  formula_var  <- vector(mode = "list", length = n_eq)
  if (is1)
  {
    for (i in seq_len(n_eq))
    {
      formula[[i]] <- as.formula(formula[[i]])
      frm_tmp      <- paste(formula[[i]][2], formula[[i]][3], sep = '~')
      is_het[i]    <- grepl(x = frm_tmp, pattern = "|", fixed = TRUE)
      if (is_het[i])
      {
        frm_tmp           <- formula_split(formula[[i]])
        formula_mean[[i]] <- frm_tmp[[1]]
        formula_var[[i]]  <- frm_tmp[[2]]
      }
      else
      {
        formula_mean[[i]] <- formula[[i]]
      }
    }
  }
  
  # Find selectivity terms
  coef_lambda_ind <- list()
  if (estimator == "2step")
  {
    coef_lambda_ind <- vector(mode = "list", length = n_eq2)
    for (v in 1:n_eq2)
    {
      formula2_terms <- attr(terms(formula2[[v]]), "term.labels")
      # add 1 for the intercept
      coef_lambda_ind[[v]] <- which(grepl("lambda", formula2_terms, 
                                          fixed = TRUE)) + 1
    }
  }
  
  # Return the results
  out <- list(is_het       = is_het,       formula     = formula,
              formula2     = formula2,     formula3    = formula3,
              formula_mean = formula_mean, formula_var = formula_var,
              coef_lambda_ind = coef_lambda_ind)
}