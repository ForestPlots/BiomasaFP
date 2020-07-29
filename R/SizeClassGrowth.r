#' @title SizeClassGrowth
#' @description Estimate mean or median growth rate for size classes 100-200mm DBH, 200-400mm DBH and >400mm DBH.
#' @param xdataset Object returned by \code{\link{mergefp}}
#' @param method Either "median" (default) or "mean". If any value other than "median" is entered the mean will be calculated.
#' @param dbh Name of column with diameter measurement to use. Default to "D4".
#' @return A data frame with PlotViewID, medain/mean growth rates for each size class and mean wood density in each plot
#' @author Martin Sullivan

#' @export


SizeClassGrowth<-function(xdataset, method="median",dbh ="D4"){
	AGBData <- xdataset
	AGBData$DBH.prev<-AGBData[match(paste(AGBData$TreeID,AGBData$PlotViewID,AGBData$Census.No-1),paste(AGBData$TreeID,AGBData$PlotViewID,AGBData$Census.No)),dbh]
	AGBData$Delta.DBH<-AGBData[,dbh]-AGBData$DBH.prev
	AGBData$Census.prev<-AGBData[match(paste(AGBData$TreeID,AGBData$PlotViewID,AGBData$Census.No-1),paste(AGBData$TreeID,AGBData$PlotViewID,AGBData$Census.No)),"Census.Mean.Date"]
	AGBData$Delta.time<-AGBData$Census.Mean.Date-AGBData$Census.prev
   	Class1<-AGBData[AGBData[,dbh]>=50 & AGBData[,dbh]<200,]
	Class2<-AGBData[AGBData[,dbh]>=200 & AGBData[,dbh]<400,]
	Class3<-AGBData[AGBData[,dbh]>=400,]
	if(method=="median"){
		G1<-aggregate(Delta.DBH/Delta.time ~ PlotViewID, data=Class1,FUN=function(x)median(x,na.rm=T))
		G2<-aggregate(Delta.DBH/Delta.time ~ PlotViewID, data=Class2,FUN=function(x)median(x,na.rm=T))
		G3<-aggregate(Delta.DBH/Delta.time ~ PlotViewID, data=Class3,FUN=function(x)median(x,na.rm=T))
	}else{
		G1<-aggregate(Delta.DBH/Delta.time ~ PlotViewID, data=Class1,FUN=function(x)mean(x,na.rm=T))
		G2<-aggregate(Delta.DBH/Delta.time ~ PlotViewID, data=Class2,FUN=function(x)mean(x,na.rm=T))
		G3<-aggregate(Delta.DBH/Delta.time ~ PlotViewID, data=Class3,FUN=function(x)mean(x,na.rm=T))
	}
		WD<-aggregate(WD~PlotViewID,data=AGBData,FUN=function(x)mean(x,na.rm=T))
	GA<-merge(G1,G2,by="PlotViewID",all.x=T)
	GB<-merge(GA,G3,by="PlotViewID",all.x=T)
	GC<-merge(GB,WD,by="PlotViewID",all.x=T)
	names(GC)<-c("PlotViewID","Class1","Class2","Class3","MeanWD")
	GC
}
