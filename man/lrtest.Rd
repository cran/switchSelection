% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lrtest.R
\name{lrtest}
\alias{lrtest}
\title{Likelihood ratio test}
\usage{
lrtest(model1, model2)
}
\arguments{
\item{model1}{the first model.}

\item{model2}{the second model.}
}
\value{
The function returns an object of class \code{'lrtest'} that is
a list with the following elements:
\itemize{
\item \code{n1} - the number of observations in the first model.
\item \code{n2} - the number of observations in the second model.
\item \code{ll1} - log-likelihood value of the first model.
\item \code{ll2} - log-likelihood value of the second model.
\item \code{df1} - the number of parameters in the first model.
\item \code{df2} - the number of parameters in the second model.
\item \code{restrictions} - the number of restrictions in the nested model.
\item \code{value} - chi-squared test statistic value.
\item \code{p_value} - p-value of the chi-squared test.
}
}
\description{
This function performs chi-squared test for nested models.
}
\details{
Arguments \code{model1} and \code{model2} should be objects
of class that has implementations of 
\code{\link[stats]{logLik}} and 
\code{\link[stats]{nobs}} methods. It is assumed that either \code{model1}
is nested into \code{model2} or vice versa. More precisely it is assumed
that the model with smaller log-likelihood value is nested into the model
with greater log-likelihood value.
}
\examples{
set.seed(123)
# Generate data according to linear regression
n <- 100
eps <- rnorm(n)
x1 <- runif(n)
x2 <- runif(n)
y <- x1 + 0.2 * x2 + eps
# Estimate full model
model1 <- lm(y ~ x1 + x2)
# Estimate restricted (nested) model
model2 <- lm(y ~ x1)
# Likelihood ratio test results
lrtest(model1, model2)
}
