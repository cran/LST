#' At-Sensor Temperature or brightness temperature
#'
#' This function calculates at-Sensor Temperature or brightness temperature
#' @importFrom raster brick as.data.frame rasterFromXYZ
#' @param Landsat_10 Raster* object, Landsat band 10
#' @param Landsat_11 Raster* object, Landsat band 11
#' @return A list containing brightness temperature corresponding to Landsat band 10 and Landsat band 11
#' @export
#' @examples
#' a <- raster::raster(ncol=100, nrow=100)
#' set.seed(2)
#' raster::values(a) = runif(10000, min=27791, max=30878)
#'
#' b <- raster::raster(ncol=100, nrow=100)
#' set.seed(2)
#' raster::values(b) = runif(10000, min=25686, max=28069)
#'
#' BT(Landsat_10 = a, Landsat_11 = b)

BT <- function(Landsat_10 = Landsat_10, Landsat_11 = Landsat_10){
  #Read the TIR band
  tir_10 <- raster::brick(Landsat_10)
  tir_11 <- raster::brick(Landsat_11)

  #Top of Atmospheric Spectral Radiance
  l_lambda_10 <- 3.3420*10^-4*tir_10 + 0.1
  l_lambda_11 <- 3.3420*10^-4*tir_11 + 0.1

  #l_lambda = ML * Qcal + AL
  #l_lambda = TOA Spectral Reflectance (watts/m2*Srad µm),
  #ML= The band-specific multiplicative rescaling factor,
  #Qcal is the Band 10 image,
  #AL is the band-specific additive rescaling factor
  #The value for ML and AL can be found in MTL.txt file

  #Conversion of Radiance to At-Sensor Temperature
  #K1_CONSTANT_BAND_10 = 774.8853
  #K2_CONSTANT_BAND_10 = 1321.0789
  BT_10 <- (1321.0789/(log(774.8853/l_lambda_10 + 1)))

  #K1_CONSTANT_BAND_11 <- 480.8883
  #K2_CONSTANT_BAND_11 <- 1201.1442
  BT_11 <- (1201.1442/(log(480.8883/l_lambda_11 + 1)))

  return(list(BT_10, BT_11))
}

#' NDVI
#'
#' Function for NDVI calculation
#' @importFrom raster brick
#' @param Red Raster* object, red band of remote sensing imagery
#' @param NIR Raster* object, NIR band of remote sensing imagery
#' @return RasterLayer
#' @export
#' @examples
#' red <- raster::raster(ncol=100, nrow=100)
#' set.seed(2)
#' raster::values(red) = runif(10000, min=0.1, max=0.4)
#'
#' NIR <- raster::raster(ncol=100, nrow=100)
#' set.seed(2)
#' raster::values(NIR) = runif(10000, min=0.1, max=0.6)
#'
#' NDVI(Red = red, NIR = NIR)

NDVI <- function(Red, NIR){
  #Read the Red and NIR band
  red <- raster::brick(Red)
  nir <- raster::brick(NIR)

  #NDVI calculation
  ndvi <- (NIR - Red)/(NIR + Red)
  return(ndvi)
}

#' Proportion of vegetation or fractional vegetation cover
#'
#' Calculation of the proportion of vegetation or fractional vegetation cover from NDVI
#' @param NDVI Raster* object, NDVI calculated from remote sensing imagery
#' @param minNDVI = 0.2 (Ref. Sobrino et al. 2004)
#' @param maxNDVI = 0.5 (Ref. Sobrino et al. 2004)
#' @return RasterLayer
#' @export
#' @examples
#' NDVI <- raster::raster(ncol=100, nrow=100)
#' set.seed(2)
#' raster::values(NDVI) = runif(10000, min=0.02, max=0.8)
#' pv(NDVI = NDVI, minNDVI = 0.2, maxNDVI = 0.5)

pv <- function(NDVI, minNDVI, maxNDVI){
  #Read the Red and NIR band
  ndvi <- raster::brick(NDVI)

  #pv calculation
  pv <- ((ndvi - minNDVI)/(maxNDVI - minNDVI))^2
  return(pv)
}
#' Land Surface Emissivity according to Van de Griend and Owe 1993
#'
#' This fumction calculates Land Surface Emissivity according to Van de Griend and Owe 1993
#' @param NDVI Raster* object, NDVI calculated from remote sensing imagery
#' @return RasterLayer
#' @export
#' @examples
#' NDVI <- raster::raster(ncol=100, nrow=100)
#' set.seed(2)
#' raster::values(NDVI) = runif(10000, min=0.02, max=0.8)
#' E_VandeGriend(NDVI)

E_VandeGriend <- function(NDVI){
  #Read the NDVI band
  ndvi <- raster::brick(NDVI)

  #Emissivity calculation
  E_VandeGriend <- 1.094 + 0.047*log(ndvi)
  #E_VandeGriend[E_VandeGriend>1] <- 1
  return(E_VandeGriend)
}

#' Land Surface Emissivity according to Valor and Caselles 1996
#'
#' This fumction calculates Land Surface Emissivity according to Valor and Caselles 1996
#' @param NDVI Raster* object, NDVI calculated from remote sensing imagery
#' @return RasterLayer
#' @export
#' @examples
#' NDVI <- raster::raster(ncol=100, nrow=100)
#' set.seed(2)
#' raster::values(NDVI) = runif(10000, min=0.02, max=0.8)
#' E_Valor(NDVI)

E_Valor <- function(NDVI){
  #Read the Red and NIR band
  ndvi <- raster::brick(NDVI)

  #pv calculation
  pv <- ((ndvi - 0.2)/(0.5 - 0.2))^2
  pv_df <- raster::as.data.frame(pv, xy = T)
  #Emissivity calculation
  E_Valor <- 0.985*pv_df[3] + 0.960*(1 - pv_df[3]) + 0.06*pv_df[3]*(1 - pv_df[3])
  #E_Valor[E_Valor>1] <- 1
  E_Valor_xyz <- cbind.data.frame(pv_df[1:2], E_Valor)
  E_Valor_raster <- raster::rasterFromXYZ(E_Valor_xyz)
  return(E_Valor_raster)
}

#' Atmospheric transmittance calculation
#'
#' This fumction calculates Atmospheric transmittance from near-surface air temperature (To, °C) and relative humidity (RH, %) of the date when Landsat passed over the study area
#' @param To Near-surface air temperature (°C) of the date when Landsat passed over the study area
#' @param RH relative humidity (%) of the date when Landsat passed over the study area
#' @return Atmospheric transmittance
#' @export
#' @examples
#' tau(To = 26, RH = 42)

tau <- function(To = To, RH = To){
  To_K <- To + 273.15
  RH_fraction <- RH/100
  W <- 0.0981*10*0.6108*exp((17.27*(To_K-273.15))/
                              (237.3+(To_K-273.15)))*RH_fraction+0.1697
  tau <- -0.0164*W^2-0.04203*W+0.9715
  return(tau)
}

#' Land Surface Emissivity according to Sobrino 2007
#'
#' This fumction calculates Land Surface Emissivity according to Sobrino 2007
#' @param red Raster* object, red band of remote sensing imagery
#' @param NDVI Raster* object, NDVI calculated from remote sensing imagery
#' @return RasterLayer
#' @export
#' @examples
#' red <- raster::raster(ncol=100, nrow=100)
#' set.seed(2)
#' raster::values(red) = runif(10000, min=0.1, max=0.4)
#' NDVI <- raster::raster(ncol=100, nrow=100)
#' set.seed(2)
#' raster::values(NDVI) = runif(10000, min=0.02, max=0.8)
#' E_Sobrino(red = red, NDVI = NDVI)

E_Sobrino <- function(red = red, NDVI = NDVI){
  #Read the Red and NIR band
  ndvi <- raster::brick(NDVI)

  #pv calculation
  pv <- ((ndvi - 0.2)/(0.5 - 0.2))^2
  E1_Sobrino <- 0.004*pv + 0.986
  E_Sobrino <- raster::raster(ndvi)
  red <- (0.979 - 0.035*red)
  #where NDVI < 0.2, take the value from red band otherwise use NA:
  E_Sobrino[] = ifelse(ndvi[]<0.2, red[], NA)

  #where NDVI is over 0.5 set E to 0.99, otherwise use whatever was in E:
  E_Sobrino[] = ifelse(ndvi[]>0.5, 0.99, E_Sobrino[])

  #if NDVI is between 0.2 and 0.5 take the value from b otherwise use whatever was in E:
  E_Sobrino[] = ifelse(ndvi[]>=0.2 & ndvi[]<=0.5, E1_Sobrino[], E_Sobrino[])

  return(E_Sobrino)
}

#' Mono window algorithm
#'
#' This fumction calculates Land Surface Temperature using mono window algorithm
#' @param BT Raster* object, brightness temperature
#' @param tau Atmospheric transmittance
#' @param E Raster* object, NDVI calculated from remote sensing imagery
#' @param To Near-surface air temperature (°C) of the date when Landsat passed over the study area
#' @return RasterLayer
#' @export
#' @examples
#' BTemp <- raster::raster(ncol=100, nrow=100)
#' set.seed(2)
#' raster::values(BTemp) = runif(10000, min=298, max=305)
#' E <- raster::raster(ncol=100, nrow=100)
#' set.seed(2)
#' raster::values(E) = runif(10000, min=0.96, max=0.99)
#' MWA(BT = BTemp, tau = 0.86, E = E, To = 26)

MWA <- function(BT = BT, tau = tau, E = E, To = To){
  BT <- raster::brick(BT)
  E <- raster::brick(E)
  To_K <- To + 273.15
  #Mean atmospheric temperature
  Ta <- 16.011 + 0.9262*To_K
  C <- E*tau
  D <- (1-tau)*(1+(1-E)*tau)
  T_MWA <- (-67.355351*(1-C-D)+(0.458606*(1-C-D)+C+D)*BT-D*Ta)/C
  return(T_MWA)
}

