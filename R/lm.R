#' Bootstrap covariance matrix for least squares estimates of linear regression
#' @description This function calculates bootstrapped covariance matrix
#' for least squares estimates of linear regression. The estimates should be
#' obtained via \code{lm} function.
#' @param model object of class \code{lm}.
#' @param iter positive integer representing the number of bootstrap iterations.
#' @details Calculations may take long time for high \code{iter} value.
#' @return This function returns a bootstrapped covariance matrix of the
#' least squares estimator.
#' @examples 
#' set.seed(123)
#' # Generate data according to linear regression
#' n <- 20
#' eps <- rnorm(n)
#' x <- runif(n)
#' y <- x + eps
#' # Estimate the model
#' model <- lm(y ~ x)
#' # Calculate bootstrap covariance matrix
#' boot(model, iter = 50)
boot <- function(model, iter = 100)
{
  if (class(model)[1] != "lm")
  {
    stop(paste0("Argument 'model' should be an object of class 'lm', ", 
                "not '", class(model)[1], "'.\n"))
  }
  if (((iter %% 1) != 0) | (iter <= 0))
  {
    stop(paste0("Wrong 'iter' argument. ", 
                "Please, insture that it is a positive integer.\n"))
  }
  data <- model.frame(model)
  coef.names <- names(data)
  n <- nrow(data)
  model.boot <- vector(mode = "list", length = iter)
  coef.boot <- matrix(NA, nrow = iter, ncol = ncol(data))
  for (i in 1:iter)
  {
    ind.boot <- sample(1:n, size = n, replace = TRUE)
    data.boot <- data[ind.boot, ]
    model.boot[[i]] <- lm(paste0(coef.names[1], " ~ ."), 
                          data = data.boot)
    coef.boot[i, ] <- coef(model.boot[[i]])
  }
  coef.names[1] <- "(Intercept)"
  cov.boot <- cov(coef.boot)
  colnames(cov.boot) <- coef.names
  rownames(cov.boot) <- coef.names
  return(cov.boot)
}

#' Leave-one-out cross-validation
#' @description This function calculates root mean squared error (RMSE) for
#' leave-one-out cross-validation of linear regression estimated via
#' least squares method.
#' @param fit object of class \code{lm}.
#' @details Fast analytical formula is used.
#' @return This function returns a numeric value representing root mean squared 
#' error (RMSE) of leave-one-out cross-validation (LOOCV).
#' @examples 
#' set.seed(123)
#' # Generate data according to linear regression
#' n <- 100
#' eps <- rnorm(n)
#' x <- runif(n)
#' y <- x + eps
#' # Estimate the model
#' model <- lm(y ~ x)
#' # Perform cross-validation
#' loocv(model)
loocv <- function(fit)
{
  if (class(fit)[1] != "lm")
  {
    stop(paste0("Argument 'fit' should be an object of class 'lm', ", 
                "not '", class(fit)[1], "'.\n"))
  }
  return(sqrt(mean((residuals(fit) / (1 - lm.influence(fit)$h)) ^ 2)))
}