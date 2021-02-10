#' Real-time data set for Sweden.
#'
#' A dataset containing real-time data for mixed and quarterly frequencies.
#'
#' @format A mixed-frequency data set of five Swedish macroeconomic variables.
#' \describe{
#'   \item{unemp}{harmonized unemployment rate (source: OECD)}
#'   \item{infl}{inflation rate (source: OECD)}
#'   \item{ip}{industrial production (source: OECD)}
#'   \item{eti}{economic tendency indicator (source: National Institute of Economic Research)}
#'   \item{gdp}{GDP growth (source: Statistics Sweden)}
#' }
#' @references
#' OECD (2016) MEI Archive: Revisions Analysis Dataset.\cr
#' Billstam, M., Fr\''{a}nd\'{e}n, J., Samuelsson, J., \"{O}sterholm, P. (2016) Quasi-Real-Time Data of the Economic Tendency Survey. Working Paper No. 143, National Institute of Economic Research.
#' Statistics Sweden (2016) Revisions, expenditure approach and hours worked at each release.
#'
"mf_sweden"

#' US Macroeconomic Data Set
#'
#' A dataset containing mixed-frequency data from FRED for three US macroeconomic variables.
#'
#' @format A list with components:
#' \describe{
#'   \item{CPIAUCSL}{inflation rate}
#'   \item{UNRATE}{unemployment rate}
#'   \item{GDPC1}{GDP growth rate}
#' }
"mf_usa"
