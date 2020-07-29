#' Estimate tree height from a height-diameter model
#'
#' @description Function to estimate height from height-diameter model parameters
#' @param dbh Numeric vector of diameters (DBH, mm).
#' @param data Dataframe containing columns reference in arguments below (e.g. model parameters).
#' @param model Character vector giving the type of height-diameter to use. Should be one of Weibull, Loglog, Env, Obs, MicMent, Hoerl, Gompertz, Polynomial or User. Obs uses measured heights.
#' User uses user defined function supplied to user.func.
#' @param a Numeric. Height-diameter model parameter.
#' @param b Numeric. Height-diameter model parameter.
#' @param c Numeric. Height-diameter model parameter.
#' @param d Numeric. Height-diameter model parameter.
#' @param height Numeric. Tree height (used when model is "Obs").
#' @param ChaveE Environmental stress parameter for height from Chave et al. 2014 Global Change Biology.
#' @param user.func User defined height-diameter function. Must use on or more of the parameters defined above.
#' @return Numeric vector with estimated heights.
#' @author Martin Sullivan

#' @export

height.mod<-function(dbh,data,model="Model",a="a_par",b="b_par",c="c_par",d="d_par",height="Height",ChaveE="ChaveE",user.func=NULL){
if(is.numeric(dbh)){
	dbh<-dbh
}else{
	dbh<-data[,dbh]
}

dbh<-dbh/10

model<-data[,model]
if(length(setdiff(model,c("Weibull","Loglog","Env","Obs","MicMent","Hoerl","Gompertz","Polynomial","User")))>0){
	stop("Model is not one of Weibull, Loglog, Env, Obs, MicMent, Hoerl, Gompertz, Polynomial or User. Please check your height.data object and see if any models don't match, e.g. because of a spelling mistake")
}

if(a %in% names(data)){
	a<-data[,a]
}

if(b %in% names(data)){
	b<-data[,b]
}

if(c %in% names(data)){
	c<-data[,c]
}

if(d %in% names(data)){
	d<-data[,d]
}

if(height %in% names(data)){
	height<-data[,height]
}

if(ChaveE %in% names(data)){
	ChaveE<-data[,ChaveE]
}



res<-rep(NA,length(dbh))
res[model=="Weibull"]<-a[model=="Weibull"]*(1-exp(-b[model=="Weibull"]*(dbh[model=="Weibull"])^c[model=="Weibull"]))
res[model=="Loglog"]<-exp(a[model=="Loglog"]+b[model=="Loglog"]*log(dbh[model=="Loglog"]))
if(ChaveE %in% names(data)){
res[model=="Env"]<-exp(0.893-ChaveE[model=="Env"]+(0.76*log(dbh[model=="Env"]))-(0.034*log(dbh[model=="Env"])^2))
}
res[model=="Obs"]<-height[model=="Obs"]
res[model=="MicMent"]<-(a[model=="MicMent"]*dbh[model=="MicMent"])/(b[model=="MicMent"]+dbh[model=="MicMent"])
res[model=="Hoerl"]<-a[model=="Hoerl"]*(b[model=="Hoerl"]^dbh[model=="Hoerl"])*(dbh[model=="Hoerl"]^c[model=="Hoerl"])
res[model=="Gompertz"]<-a[model=="Gompertz"]*exp(-exp(b[model=="Gompertz"]-(c[model=="Gompertz"]*dbh[model=="Gompertz"])))
res[model=="Polynomial"]<-a[model=="Polynomial"]+(b[model=="Polynomial"]*dbh[model=="Polynomial"])+(c[model=="Polynomial"]*dbh[model=="Polynomial"]^2)
if(!is.null(user.func)){
res[model=="User"]<-user.func(dbh[model=="User"],a[model=="User"],b[model=="User"],c[model=="User"],d[model=="User"])
}
res<-as.numeric(res)
res
}

