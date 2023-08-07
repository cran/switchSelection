#' Likelihood ratio test 
#' @description This function performs chi-squared test for nested models.
#' @param model1 the first model.
#' @param model2 the second model.
#' @details Arguments \code{model1} and \code{model2} should be objects
#' of class that has implementations of 
#' \code{\link[stats]{logLik}} and 
#' \code{\link[stats]{nobs}} methods. It is assumed that either \code{model1}
#' is nested into \code{model2} or vice versa. More precisely it is assumed
#' that the model with smaller log-likelihood value is nested into the model
#' with greater log-likelihood value.
#' @return The function returns an object of class \code{'lrtest'} that is
#' a list with the following elements:
#' \itemize{
#' \item \code{n1} - the number of observations in the first model.
#' \item \code{n2} - the number of observations in the second model.
#' \item \code{ll1} - log-likelihood value of the first model.
#' \item \code{ll2} - log-likelihood value of the second model.
#' \item \code{df1} - the number of parameters in the first model.
#' \item \code{df2} - the number of parameters in the second model.
#' \item \code{restrictions} - the number of restrictions in the nested model.
#' \item \code{value} - chi-squared test statistic value.
#' \item \code{p_value} - p-value of the chi-squared test.
#' }
#' @examples 
#' set.seed(123)
#' # Generate data according to linear regression
#' n <- 100
#' eps <- rnorm(n)
#' x1 <- runif(n)
#' x2 <- runif(n)
#' y <- x1 + 0.2 * x2 + eps
#' # Estimate full model
#' model1 <- lm(y ~ x1 + x2)
#' # Estimate restricted (nested) model
#' model2 <- lm(y ~ x1)
#' # Likelihood ratio test results
#' lrtest(model1, model2)
lrtest <- function(model1, model2)
{
  out <- list()
  
  out$n1 <- nobs(model1)
  out$n2 <- nobs(model2)
  
  if (out$n1 != out$n2)
  {
    stop("Models should have the same number of observations.")
  }
  
  out$ll1 <- logLik(model1)
  out$ll2 <- logLik(model2)
  
  out$df1 <- attr(out$ll1, which = "df")
  out$df2 <- attr(out$ll2, which = "df")
  
  if ((which.max(c(out$ll1, out$ll2)) != 
       which.max(c(out$df1, out$df2))) |
      (out$df1 == out$df2))
  {
    stop(paste0("Models should be nested."))
  }
  
  out$value <- as.numeric(2 * abs(out$ll1 - out$ll2))
  out$restrictions <- abs(out$df1 - out$df2)
  out$p_value <- 1 - pchisq(out$value, df = out$restrictions)
  class(out) <- "lrtest";
  
  return(out)
}

#' Print Method for Likelihood Ratio Test
#' @description Prints summary for an object of class 'lrtest'.
#' @param x object of class "lrtest".
#' @param ... further arguments (currently ignored).
#' @return The function returns input argument \code{x}.
print.lrtest <- function(x, ...)
{
  cat("LR-test results\n")
  cat("---\n")
  cat(paste0("It is assumed that model",
             which.min(c(x$df1, x$df2)),
             " is nested into model",
             which.max(c(x$df1, x$df2)), "\n"))
  cat(paste0("H0: model", which.min(c(x$df1, x$df2)), " is correct\n"))
  cat("--\n")
  cat(paste0("Restrictions = ", x$restrictions, "\n"))
  cat(paste0("Test statistic = ", round(x$value, 3), "\n"))
  cat(paste0("p-value = ", round(x$p_value, 3), 
             starsVector(x$p_value),
             "\n"))
  if (getOption("show.signif.stars"))
  {
    cat("--\n")
    cat("Signif. codes:  0 '***' 0.01 '**' 0.05 '*' 0.1 \n")
  }
  return(x)
}

#' Summary Method for Likelihood Ratio Test
#' @description Provides summary for an object of class 'lrtest'.
#' @param object object of class "lrtest"
#' @param ... further arguments (currently ignored)
#' @details This function just changes the class of the 'lrtest'
#' object to 'summary.lrtest'.
#' @return Returns an object of class 'summary.lrtest'.
summary.lrtest <- function(object, ...) 
{
  if (length(list(...)) > 0)
  {
    warning("Additional arguments passed through ... are ignored.")   
  }
  
  class(object) <- "summary.lrtest";
  
  return(object)
}

#' Print Summary Method for Likelihood Ratio Test
#' @description Prints summary for an object of class 'lrtest'.
#' @param x object of class "lrtest"
#' @param ... further arguments (currently ignored)
#' @return The function returns input argument \code{x} changing
#' it's class to \code{lrtest}.
print.summary.lrtest <- function(x, ...)
{
  class(x) <- "lrtest"
  print(x)
  return(x)
}
