% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fitted.R
\name{fitted.mnprobit}
\alias{fitted.mnprobit}
\title{Extract Model Fitted Values}
\usage{
\method{fitted}{mnprobit}(object, ..., newdata = NULL)
}
\arguments{
\item{object}{object of class 'mnprobit'.}

\item{...}{further arguments (currently ignored).}

\item{newdata}{an optional data frame in which to look for variables 
with which to predict. If omitted, the original data frame used. 
This data frame should contain values of dependent variables even if they 
are not actually needed for prediction (simply assign them with 0 values).}
}
\value{
Returns a data frame. Its first column provides an index of the
most probable alternative. Columns which names coincide with the
names of the continuous equation provide unconditional expectation of
the dependent variable in available regimes.
}
\description{
Extracts fitted values from 'mnprobit' object
}
