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
#' 
#' Arguments \code{model1} and \code{model2} may be the lists of models.
#' If \code{model1} is a list of models then it is assumed that the number
#' of degrees of freedom and log-likelihood of the first model are just a sum
#' of degrees of freedom and log-likelihoods of the models in this list.
#' Similarly for \code{model2}. 
#' 
#' If \code{model1} or \code{model2} is a list then the number of observations 
#' of the associated models are calculated as the sum of the numbers of 
#' observations of the models in corresponding lists.
#' However sometimes it may be misleading. For example, when bivariate probit
#' model (full) is tested against two independent probit models (restricted). 
#' Then it will be assumed that the number of observations in the restricted
#' model is twice the number of observations in the full model that is not the
#' case. 
#' Fortunately it will not affect the results of the likelihood ratio test.
#' @return The function returns an object of class \code{'lrtest_msel'} that is
#' a list with the following elements:
#' \itemize{
#' \item \code{n1} - the number of observations in the first model.
#' \item \code{n2} - the number of observations in the second model.
#' \item \code{ll1} - log-likelihood value of the first model.
#' \item \code{ll2} - log-likelihood value of the second model.
#' \item \code{df1} - the number of parameters in the first model.
#' \item \code{df2} - the number of parameters in the second model.
#' \item \code{restrictions} - the number of restrictions in the nested model.
#' \item \code{value} - chi-squared (likelihood ratio) test statistic value.
#' \item \code{p_value} - p-value of the chi-squared (likelihood ratio) test.
#' }
#' @examples 
#' set.seed(123)
#' # Generate data according to linear regression
#' n   <- 100
#' eps <- rnorm(n)
#' x1  <- runif(n)
#' x2  <- runif(n)
#' y   <- x1 + 0.2 * x2 + eps
#' # Estimate full model
#' model1 <- lm(y ~ x1 + x2)
#' # Estimate restricted (nested) model
#' model2 <- lm(y ~ x1)
#' # Likelihood ratio test results
#' lrtest_msel(model1, model2)
lrtest_msel <- function(model1, model2)
{
  # List to save the results
  out <- list()
  
  # Control for multiple models input
  multiple1 <- !(paste0("nobs.", class(model1)) %in% 
                 methods(class = class(model1)))
  multiple2 <- !(paste0("nobs.", class(model2)) %in% 
                 methods(class = class(model2)))
  
  # First model
  if (!multiple1)
  {
    out$n1 <- nobs(model1)                                 
    out$ll1 <- logLik(model1)                            
    out$df1 <- attr(out$ll1, which = "df")                 
  }
  else
  {
    out$n1 <- 0
    out$ll1 <- 0
    out$df1 <- 0
    for (i in 1:length(model1))
    {
      out$n1 <- out$n1 + as.numeric(nobs(model1[[i]])) 
      ll1 <- logLik(model1[[i]])
      out$ll1 <- out$ll1 + as.numeric(ll1)                  
      out$df1 <- out$df1 + as.numeric(attr(ll1, which = "df"))
    }
  }
  
  # Second model
  if (!multiple2)
  {
    out$n2 <- nobs(model2)                                 
    out$ll2 <- logLik(model2)                             
    out$df2 <- attr(out$ll2, which = "df")                
  }
  else
  {
    out$n2 <- 0
    out$ll2 <- 0
    out$df2 <- 0
    for (i in 1:length(model2))
    {
      out$n2 <- out$n2 + as.numeric(nobs(model2[[i]]))    
      ll2 <- logLik(model2[[i]])
      out$ll2 <- out$ll2 + as.numeric(ll2) 
      out$df2 <- out$df2 + as.numeric(attr(ll2, which = "df"))  
    }
  }
  
  # Validation
  if (!(multiple1 | multiple2))
  {
    if (out$n1 != out$n2)
    {
      stop("Models should have the same number of observations.")
    }
    
    if ((which.max(c(out$ll1, out$ll2)) != 
         which.max(c(out$df1, out$df2))) |
        (out$df1 == out$df2))
    {
      stop(paste0("Models should be nested."))
    }
  }
  
  # Calculation of the test statistics and p-value
  out$value <- as.numeric(2 * abs(out$ll1 - out$ll2))
  out$restrictions <- abs(out$df1 - out$df2)
  out$p_value <- 1 - pchisq(out$value, df = out$restrictions)
  class(out) <- "lrtest_msel";
  
  return(out)
}

#' Print Method for Likelihood Ratio Test
#' @description Prints summary for an object of class 'lrtest_msel'.
#' @param x object of class "lrtest_msel".
#' @param ... further arguments (currently ignored).
#' @return The function returns the input argument \code{x}.
print.lrtest_msel <- function(x, ...)
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
  if (x$p_value >= 0.00001)
  {
    cat(paste0("p-value = ", round(x$p_value, 5), starsVector(x$p_value), "\n"))
  }
  else
  {
    cat(paste0("p-value < 0.00001 ", starsVector(x$p_value), "\n"))
  }
  if (getOption("show.signif.stars"))
  {
    cat("--\n")
    cat("Signif. codes:  0 '***' 0.01 '**' 0.05 '*' 0.1 \n")
  }
  return(x)
}

#' Summary Method for Likelihood Ratio Test
#' @description Provides summary for an object of class 'lrtest_msel'.
#' @param object object of class "lrtest_msel"
#' @param ... further arguments (currently ignored)
#' @details This function just changes the class of the 'lrtest_msel'
#' object to 'summary.lrtest_msel'.
#' @return Returns an object of class 'summary.lrtest_msel'.
summary.lrtest_msel <- function(object, ...) 
{
  if (length(list(...)) > 0)
  {
    warning("Additional arguments passed through ... are ignored.")   
  }
  
  class(object) <- "summary.lrtest_msel";
  
  return(object)
}

#' Print Summary Method for Likelihood Ratio Test
#' @description Prints summary for an object of class 'lrtest_msel'.
#' @param x object of class "lrtest_msel"
#' @param ... further arguments (currently ignored)
#' @return The function returns input argument \code{x} changing
#' it's class to \code{lrtest_msel}.
print.summary.lrtest_msel <- function(x, ...)
{
  class(x) <- "lrtest_msel"
  print(x)
  return(x)
}
