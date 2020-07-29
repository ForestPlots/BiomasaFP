#' @title height.constrain
#' @description Modify selection of fitted height-diameter models to avoid models that overpredict over a threshold maximum height
#' @param fit Fitted height-diameter models object returned by \code{\link{hd.fit}}.
#' @param data Object returned by \code{\link{mergefp}}.
#' @param level Level to fit height-diameter model. One of "Plot", "Cluster", "BioR" and "Continent". Must match level used whe constructing fit object.
#' @param dbh Name of column containing diameter measurements. Default is "D4".
#' @param max.size Maximum tree height. If all models predict the height of the largest tree in the plot to be below this, then the model with the lowest AGB prediction error will be used. Otherwise models predicting heights above this will be excluded from selection. If all models predict heights above this, then the model giving the lowest maximum height will be selected.
#' @param na.too.big Logical. Should parameters be set to NA if all models exceed height treshold. Default is FALSE.
#' @param forest.type Logical. Should ID include forest type.
#' @return Dataframe with components returned by \code{\link{hd.fit}}, but with addtional columns giving the maximum diameter of the largest tree in the plot, the predicted maximum height using each model, and a column indicating whether the choice of best model was modified (0 = no change, 1 = previous best model excluded due to exceeding threshold, 2 = all models exceeded height threshold).
#' @author Martin Sullivan

#' @export

height.constrain<-function(fit,data,level,dbh="D4",max.size=80,na.too.big=FALSE,forest.type=TRUE){
data$level<-level
data$FT<-paste(data$ForestMoistureID,data$ForestEdaphicHeightID,data$ForestElevationHeightID,sep="")
if(forest.type==TRUE){
data$ID<-ifelse(data$level=="Plot",data$PlotID,ifelse(data$level=="Cluster",paste(data$ClusterID,data$FT,sep="_"),ifelse(data$level=="BioR",paste(data$BiogeographicalRegionID,data$FT,sep="_"),paste(data$Continent,data$FT,sep="_"))))
}else{
data$ID<-ifelse(data$level=="Plot",data$PlotID,ifelse(data$level=="Cluster",data$ClusterID,ifelse(data$level=="BioR",data$BiogeographicalRegionID,data$Continent)))
}
max.diam<-tapply(data[,dbh],data$ID,max,na.rm=T)
fit$Max.diameter<-as.numeric(max.diam[match(fit$ID,names(max.diam))])
fit$Max.best<-with(fit,height.mod(Max.diameter,Best_a,Best_b,Best_c))
fit$Max.Weib<-with(fit,height.mod(Max.diameter,a_Weib,b_Weib,c_Weib))
fit$Max.WeightWeib<-with(fit,height.mod(Max.diameter,a_WeightWeib,b_WeightWeib,c_WeightWeib))
fit$Max.Loglog<-with(fit,height.mod(Max.diameter,a_Loglog,b_Loglog,c_Loglog))
fit$Bad.best<-ifelse(fit$Max.best>max.size,1,0)
fit$Bad.Weib<-ifelse(fit$Max.Weib>max.size,1000,1)
fit$Bad.WeightWeib<-ifelse(fit$Max.WeightWeib>max.size,1000,1)
fit$Bad.Loglog<-ifelse(fit$Max.Loglog>max.size,1000,1)
fit$Weib.error2<-abs(fit$Weibull_error*fit$Bad.Weib)
fit$WeightWeib.error2<-abs(fit$Weibull.weight_error*fit$Bad.WeightWeib)
fit$Loglog.error2<-abs(fit$Loglog_error*fit$Bad.Loglog)
fit$BestMod2<-as.numeric(apply(fit[,c("Weib.error2","WeightWeib.error2","Loglog.error2")],1,function(x)grep(min(x,na.rm=T),x)))
fit$All.bad<-as.numeric(apply(fit[,c("Bad.Weib","Bad.WeightWeib","Bad.Loglog")],1,function(x)as.numeric(sum(x,na.rm=T)/length(x[!is.na(x)])>=1000)))
fit$Min.height<-as.numeric(apply(fit[,c("Max.Weib","Max.WeightWeib","Max.Loglog")],1,function(x)grep(min(x,na.rm=T),x)))
fit$BestMod3<-as.numeric(ifelse(fit$All.bad==0,fit$BestMod2,fit$Min.height))
fit$Modification<-ifelse(fit$BestMod==fit$BestMod3,0,ifelse(fit$BestMod3==as.numeric(fit$BestMod2),1,2))
fit$BestMod<-fit$BestMod3
fit$Best_a<-ifelse(fit$BestMod==1,fit[,1],ifelse(fit$BestMod==2,fit[,4],fit[,7]))
fit$Best_b<-ifelse(fit$BestMod==1,fit[,2],ifelse(fit$BestMod==2,fit[,5],fit[,8]))
fit$Best_c<-ifelse(fit$BestMod==1,fit[,3],ifelse(fit$BestMod==2,fit[,6],fit[,9]))
fit$Max.best<-with(fit,height.mod(Max.diameter,Best_a,Best_b,Best_c))
if(na.too.big==TRUE){
fit$Best_a<-ifelse(fit$All.bad==1,NA,fit$Best_a)
fit$Best_b<-ifelse(fit$All.bad==1,NA,fit$Best_b)
fit$Best_c<-ifelse(fit$All.bad==1,NA,fit$Best_c)
}
fit<-fit[,c(1:23,35)]
return(fit)
}
