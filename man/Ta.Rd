% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LST.R
\name{Ta}
\alias{Ta}
\title{Mean atmospheric temperature}
\usage{
Ta(To = To, mod = mod)
}
\arguments{
\item{To}{Near-surface air temperature (°C) of the date when Landsat passed over the study area}

\item{mod}{A string specifying which model to use. It can be anyone of "USA 1976 Standard" or
"Tropical Region" or "Mid-latitude Summer Region" or "Mid-latitude Winter Region"}
}
\value{
Mean atmospheric temperature (K)
}
\description{
This function calculates mean atmospheric temperature (Ta) using near-surface air
temperature (To)
}
\examples{
Ta(To = 26, mod = "Mid-latitude Winter Region")
}
