#' Allometric equations
#'
#' Functions to estimate biomass from diameter, height and wood density. If specifying a new function please include all these parameters, even if they are not used (as other functions will pass all these arguments to the function)
#'
#' @param d numeric vector giving tree diameter (DBH) in cm
#' @param h numeric vector giving tree height (m)
#' @param wd numeric vector giving wood density in g/m^3


#' @describeIn AGBChv14 Chave et al. 2014 equation with height (see Chave J, Rejou-Mechain M, Burquez A et al. 2014. Improved allometric models to estimate the aboveground biomass of tropical trees. Global Change Biology 20: 3177-3190. doi: 10.1111/gcb.12629)
AGBChv14<-function(d,h,wd){
	0.0673 *(wd * (d/10)^2* h)^0.976/1000
}

#' @describeIn AGBChv14 Basal area
BA<-function(d,h,wd){
	pi*(d/2000)^2
}

#' @describeIn AGBChv14 Basal area weighted by wood density
wdBA<-function(d,h,wd){
	((pi*(d/20)^2)*wd)/1000
}

#' @describeIn AGBChv14 Moist forest equation (without height) from Chave et al. 2005 (see Chave, J., Andalo, C., Brown, S. et al. Tree allometry and improved estimation of carbon stocks and balance in tropical forests. Oecologia 145, 87–99 (2005). https://doi.org/10.1007/s00442-005-0100-x)
AGBChv05M<-function(d,h,wd){
	wd * exp (-1.499 + (2.148*log(d/10))+ (0.207*(log(d/10))^2)- (0.0281*(log(d/10))^3))/1000
}

#' @describeIn AGBChv14 Height-based tree fern equation from Tiepolo et al.2002 (see Tiepolo, G., M. Calmon, et al. (2002). Measuring  and monitoring  carbon  stocks  at  the  guaraqueçaba  climate  action project, Paraná, Brazil. International Symposium on Forest Carbon Sequestration and Monitoring, Extension Serie Taiwan Forestry Research Institute: 98-115.)
TiepoloTreeFern<-function(d,h,wd){
  agb<--4266348/(1-(2792284*exp(-0.313677*h)))
  agb<-agb/1000
  agb
}

#' @describeIn AGBChv14 Diameter only palm equation from Goodman et al.2013 (see Goodman, R. et al. 2013. Amazon palm biomass and allometry, Forest Ecology and Management 310: 994-1004 https://doi.org/10.1016/j.foreco.2013.09.045)
GoodmanPalm<-function(d,h,wd){
  ln.AGB<--3.3488+(2.7483*log(d/10))
  return((exp(ln.AGB))/1000)
}



