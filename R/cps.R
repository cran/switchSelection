#' A subset of data from Current Population Survey (CPS).
#'
#' Labor market data on 18,253 middle age (25-54 years) married women
#' in the year 2022.
#'
#' @docType data
#'
#' @format
#' A data frame with 18,253 rows and 23 columns. It contains information on
#' wages and some socio-demographic characteristics of the middle age
#' (25-54 years) married women:
#' \describe{
#'   \item{age}{the age measured in years.}
#'   \item{sage}{the same as age but for a spouse.}
#'   \item{work}{a binary variable for the employment status 
#'   (0 - unemployed, 1 - employed).}
#'   \item{swork}{the same as work but for a spouse.}
#'   \item{nchild}{the number of children under age 5.}
#'   \item{snchild}{the same as nchild but for a spouse.}
#'   \item{health}{subjective health status 
#'   (1 - poor, 2 - fair, 3 - good, 4 - very good, 5 - excellent).}
#'   \item{shealth}{the same as health but for a spouse.}
#'   \item{basic}{a binary variable which equals 1 for those who have graduated 
#'   from high school or has at least some college or has associated degree and 
#'   does not have any higher level of education, 0 - otherwise.}
#'   \item{bachelor}{a binary variable which equals 1 for those whose highest 
#'   education level is a bachelor degree.}
#'   \item{master}{a binary variable which equals 1 for those whose highest 
#'   education level is a master degree.}
#'   \item{sbasic}{the same as basic but for a spouse.}
#'   \item{sbachelor}{the same as bachelor but for a spouse.}
#'   \item{smaster}{the same as master but for a spouse.} 
#'   \item{educ}{a categorical variable for the level of education such that 
#'   educ = 0 if basic = 1, 
#'   educ = 1 if bachelor = 1 and 
#'   educ = 2 if master = 1.}
#'   \item{seduc}{the same as educ but for a spouse.}
#'   \item{weeks}{a total number of weeks worked durning the year.}
#'   \item{sweeks}{the same as weeks but for a spouse.}
#'   \item{hours}{a usual number of working hours per week.}
#'   \item{shours}{the same as hours but for a spouse.}
#'   \item{wage}{the wage of the individual.}
#'   \item{swage}{the same as wage but for a spouse.}
#'   \item{lwage}{an inverse hyperbolic sine transformation of the hourly wage.}
#'   \item{slwage}{the same as lwage but for a spouse.}
#'   \item{state}{a state of residence.}
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
#' model <- msel(work ~ age + bachelor + master, data = cps)
#' summary(model)
#' } 
"cps"