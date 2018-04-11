#' Atmospheric pressure and wind curl data
#'
#' Data on derived winds and wind curl from Pacific Fisheries Environmental
#' Laboratory, pfeg.noaa.gov, collected from 15 stations on the west coast of North America.
#'
#' #' @format A data frame with 12885 rows and 6 variables:
#' \describe{
#'   \item{year}{year, 1946 to 2017}
#'   \item{month}{month, 1 to 12}
#'   \item{lat}{Latitude, in degrees N. First three characters are degrees}
#'   \item{lon}{Longitude, in degrees W. First three characters are degrees}
#'   \item{atmos_pres}{representing atmospheric pressure in millibars ((value +10,000)/10)}
#'   \item{wind_curl}{the wind stress curl or value X 1.0E-10 = curl index in dynes/centimeter**2/centimeter}
#' }
#'
#' @docType data
#'
#' @usage data(wind)
#'
#' @keywords datasets
#'
#' @references
#' (\href{https://www.pfeg.noaa.gov/products/PFEL/modeled/indices/transports/transports.html})
#'
#' @source \href{http://qtlarchive.org/db/q?pg=projdetails&proj=moore_2013b}{QTL Archive}
#'
#' @examples
#' data(wind)
"wind"
