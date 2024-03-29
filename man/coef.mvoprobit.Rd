% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coef.R
\name{coef.mvoprobit}
\alias{coef.mvoprobit}
\title{Coefficients extraction method for mvoprobit.}
\usage{
\method{coef}{mvoprobit}(object, ..., eq = NULL, eq2 = NULL, regime = NULL, type = "coef")
}
\arguments{
\item{object}{object of class "mvoprobit".}

\item{...}{further arguments (currently ignored).}

\item{eq}{integer representing an index of the ordered equation.}

\item{eq2}{integer representing an index of the continuous equation.}

\item{regime}{integer representing a regime of the continuous equation.}

\item{type}{character representing a type of the output. Possible options
are \code{"coef"}, \code{"coef2"}, \code{"cov"}, \code{"cov1"}, \code{"var"},
\code{"cov2"}, \code{"cov3"}, \code{coef_lambda} and \code{marginal}.
See 'Details' for additional information.}
}
\value{
See 'Details' section.
}
\description{
Extract coefficients and other estimates from mvoprobit object.
}
\details{
Consider notations from the 'Details' section of 
\code{\link[switchSelection]{mvoprobit}}.

Suppose that \code{type = "coef"}. Then estimates of \eqn{\gamma_{j}} 
coefficients are returned for each \eqn{j\in\{1,...,J\}}. 
If \code{eq = j} then only estimates of \eqn{\gamma_{j}} coefficients 
are returned.

Suppose that \code{type = "coef_var"}. Then estimates of \eqn{\gamma_{j}^{*}} 
coefficients are returned for each \eqn{j\in\{1,...,J\}}. 
If \code{eq = j} then only estimates of \eqn{\gamma_{j}^{*}} coefficients 
are returned.

Suppose that \code{type = "coef2"}. Then estimates of \eqn{\beta_{r}} 
coefficients are returned for each \eqn{r\in\{0,...,R - 1\}}. 
If \code{eq2 = k} then only estimates for the \eqn{k}-th continuous equation 
are returned. If \code{regime = r} then estimates of \eqn{\beta_{r}}
coefficients are returned for the \code{eq2}-th continuous equation. 
Herewith if \code{regime} is not \code{NULL} and \code{eq2} is \code{NULL} 
it is assumed that \code{eq2 = 1}.

Suppose that \code{type = "cov"}. Then estimate of the asymptotic covariance
matrix of the estimator is returned. Note that this estimate depends
on the \code{cov_type} argument of \code{\link[switchSelection]{mvoprobit}}.

Suppose that \code{type = "cov1"}. Then estimate of the covariance matrix of
\eqn{u_{i}} is returned. If \code{eq = c(a, b)} then the function returns
\eqn{(a, b)}-th element of this matrix i.e. an element from 
\code{a}-th row and \code{b}-th column.

Suppose that \code{type = "cov12"}. Then estimates of covariances between
\eqn{u_{i}} and \eqn{\varepsilon_{i}} are returned. If \code{eq2 = k} then
covariances with random errors of the \code{k}-th continuous equation are
returned. If in addition \code{eq = j} and \code{regime = r} then the 
function returns estimate of \eqn{Cov(u_{ji}, \varepsilon_{ri})} for the
\code{k}-th equation. If \code{eq2 = NULL} it is assumed that 
\code{eq2 = 1}.

Suppose that \code{type = "var"} or \code{type = "cov2"}. Then estimates of
the variances of \eqn{\varepsilon_{i}} are returned. If \code{eq2 = k}
then estimates only for \eqn{k}-th continuous equation are returned.
If in addition \code{regime = r} then estimate of \eqn{Var(\varepsilon_{ri})}
is returned. Herewith if \code{regime} is not \code{NULL} and
\code{eq2} is \code{NULL} it is assumed that \code{eq2 = 1}.

Suppose that \code{type = "cov3"}. Then estimates of the covariances between
random errors of different equations in different regimes are returned.
If \code{eq2 = c(a, b)} and \code{regime = c(c, d)} then function returns
an estimate of the covariance of random errors of the 
\code{a}-th and \code{b}-th 
continuous equations in regimes \code{c} and \code{d} correspondingly. 
If this covariance is not identifiable then \code{NA} value is returned.

Suppose that \code{type = "coef_lambda"}. Then estimates of the coefficients
for \eqn{\hat{\lambda}^{t}_{ji}} are returned i.e. 
estimates of \eqn{\tau_{jt}} for each regime. 
If \code{regime = r} then estimates are returned for the \eqn{r}-th 
regime. If in addition \code{eq = j} then only estimates for this \eqn{j} 
are returned.
}
