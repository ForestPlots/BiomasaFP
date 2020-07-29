#' @title hd.fit
#' @description Fit height-diameter models for each level in a dataset. This is a wrapper for \code{\link{fit.weib.models}}.
#' @param htdata Dataframe containing height measurements.
#' @param level Level to fit height-diameter model. One of "Plot", "Cluster", "BioR" and "Continent".
#' @param weib.only Logical. Should the search for the best model be constrained to Weibull models. If TRUE, log-log models are only selected when neither Weibull models converge. Default is FALSE.
#' @param thresh thresh argument to be passed to /code{fit.weib.models}
#' @param forest.type Logical. If true (default), level is constrained to be within a particular forest type/
#' @return Dataframe with estimated parameters, number of height measurments used to fit the models and the error at estimating AGB compared to AGB estimates using observed heights (in units of Mg per tree). Parameters of the 'best' model (i.e. lowest AGB prediction error) are also given.
#' @author Martin Sullivan

#' @export

hd.fit<-function(htdata,level,weib.only=FALSE,thresh=20,forest.type=TRUE){
htdata$level<-level
htdata$FT<-paste(htdata$ForestMoistureID,htdata$ForestEdaphicHeightID,htdata$ForestElevationHeightID,sep="")
if(forest.type==TRUE){
htdata$ID<-ifelse(htdata$level=="Plot",htdata$PlotID,ifelse(htdata$level=="Cluster",paste(htdata$ClusterID,htdata$FT,sep="_"),ifelse(htdata$level=="BioR",paste(htdata$BiogeographicalRegionID,htdata$FT,sep="_"),paste(htdata$Continent,htdata$FT,sep="_"))))
}else{
htdata$ID<-ifelse(htdata$level=="Plot",htdata$PlotID,ifelse(htdata$level=="Cluster",htdata$ClusterID,ifelse(htdata$level=="BioR",htdata$BiogeographicalRegionID,htdata$Continent)))
}
htdata2<-split(htdata,f=htdata$ID)
a<-lapply(htdata2,fit.weib.models,weib.only=weib.only,thresh=thresh)
b<-do.call(rbind,a)
id.vec<-rownames(b)
b<-data.frame(b)
b$ID<-id.vec
return(b)
}
