#' A subset of data from Current Population Survey (CPS).
#'
#' Labor market data on 18,253 middle age (25-54 years) married women
#' in the year 2022.
#'
#' @docType data
#'
#' @format
#' A data frame with 18,253 rows and 13 columns. It contains information on
#' wages and some socio-demographic characteristics of middle age
#' (25-54 years) married women:
#' \describe{
#'   \item{age}{age of individual measured in years.}
#'   \item{lwage}{logarithm of hourly wage.}
#'   \item{slwage}{logarithm of hourly wage of a spouse.}
#'   \item{work}{binary variable for employment status 
#'   (0 - unemployed, 1 - employed).}
#'   \item{swork}{binary variable for employment status 
#'   of a spouse (0 - unemployed, 1 - employed).}
#'   \item{nchild}{the number of children under age 5.}
#'   \item{health}{subjective health status (1 - poor, 2 - fair, 3 - good, 
#'   4 - very good, 5 - excellent).}
#'   \item{basic}{binary variable which equals 1 for those who have graduated 
#'   from high school or has at least some college or has associated degree and 
#'   does not have any higher level of education, 0 - otherwise.}
#'   \item{bachelor}{binary variable which equals 1 for those whose highest 
#'   education level is a bachelor degree.}
#'   \item{master}{binary variable which equals 1 for those whose highest 
#'   education level is a master degree.}
#'   \item{sbasic}{the same as basic but for a spouse.}
#'   \item{sbachelor}{the same as bachelor but for a spouse.}
#'   \item{smaster}{the same as master but for a spouse.} 
#'   ...
#' }
#' @usage data(cps)
#' @references Flood S, King M, Rodgers R, Ruggles S, Warren R, 
#' Westberry M (2022). Integrated Public Use Microdata Series, 
#' Current Population Survey: Version 10.0 [dataset].
#' doi: 10.18128/D030.V10.0.
#' @source <https://www.census.gov/programs-surveys/cps.html>
#' @examples
#' \donttest{
#' data(cps)
#' model <- mvoprobit(work ~ age + bachelor + master, data = cps)
#' summary(model)
#' } 
"cps"