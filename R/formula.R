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

