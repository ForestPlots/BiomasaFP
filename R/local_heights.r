#' @title local.heights
#' @description Wrapper function for estimating heights using local allometric equations.
#' @param xdataset Object returned by \code{\link{mergefp}}.
#' @param dbh Character giving the name of the column contining the diameter measurements. Defaults to "D4".
#' @param use.height.constrain Logical. Should height constrain function be used when selecting best model. Defaul is FALSE.
#' @param max.size max.size argument to \code{\link{height.constrain}}.
#' @param thresh thresh argument to \code{\link{fit.weib.models}}.
#' @param na.too.big na.too.big argumet to \code{\link{height.constrain}}.
#' @param laser.cor Logical. Should laser heights be adjusted to correct for bias. Default is FALSE.
#' @param no.plot Logial. Should the best parameters exclude plot and cluster level. Default is TRUE.
#' @return Dataframe with best parameters for height-diameter models for each row in xdataset.
#' @author Martin Sullivan

#' @export

local.heights<-function(xdataset,dbh="D4",use.height.constrain=FALSE,max.size=80,thresh=20,na.too.big=FALSE,laser.cor=FALSE,no.plot=TRUE){
data<-xdataset
data$FT<-paste(data$ForestMoistureID,data$ForestEdaphicHeightID,data$ForestElevationHeightID,sep="")
data$ID.clust<-paste(data$ClusterID,data$FT,sep="_")
data$ID.clust[is.na(data$ClusterID)]<-NA
data$ID.bio<-paste(data$BiogeographicalRegionID,data$FT,sep="_")
data$ID.continent<-paste(data$Continent,data$FT,sep="_")
#Filter to select height data
data$Alive <- as.numeric (ifelse(data$F2==1, 1 ,0))
TreesHt <- data[data$Height>0 & data$Alive==1 & data$D1>90 & data$D1<5000    & !is.na(data$Height) & data$Height<90
                   & data$Monocot==0 &  !is.na(data$F5), ]
TreesHt$Method<- ifelse(TreesHt$F5==1 | TreesHt$F5==2,1,
                        ifelse(TreesHt$F5==6,6,
                               ifelse( TreesHt$F5==3, 3,4
                               )
                        )
	)

# Exclude F3-3
TreesHt <- TreesHt[ grepl('3',TreesHt$F3)==FALSE,  ]
# Exclude F4 not like '60' or '0' or '2'
TreesHt <- TreesHt[ grepl('0',TreesHt$F4)==TRUE|grepl('06',TreesHt$F4)==TRUE |grepl('2',TreesHt$F4)==TRUE ,  ]
# Exclude f1 flag1= b, c, d, j, k
htdata <- TreesHt[ grepl('b',TreesHt$F1)==FALSE & grepl('c',TreesHt$F1)==FALSE &
                            grepl('d',TreesHt$F1)==FALSE & grepl('j',TreesHt$F1)==FALSE&
                            grepl('k',TreesHt$F1)==FALSE
                    ,  ]

#Exclude treeswith method 1
htdata <- htdata[!(htdata$Method==1),]

#Correct laser measurements for bias
if(laser.cor==TRUE){
htdata$Height<-ifelse(htdata$Method==4,2.24*htdata$Height^0.8,htdata$Height)
}

#Make column with selected diameter
htdata$Diameter<-htdata[,dbh]

#Fit h-d models
plot.fit<-hd.fit(htdata,"Plot",thresh=thresh)
clust.fit<-hd.fit(htdata,"Cluster",thresh=thresh)
bio.fit<-hd.fit(htdata,"BioR",weib.only=TRUE,thresh=thresh)
cont.fit<-hd.fit(htdata,"Continent",weib.only=TRUE,thresh=thresh)
clust.all.fit<-hd.fit(htdata,"Cluster",thresh=thresh,forest.type=FALSE)
bio.all.fit<-hd.fit(htdata,"BioR",weib.only=TRUE,thresh=thresh,forest.type=FALSE)
cont.all.fit<-hd.fit(htdata,"Continent",weib.only=TRUE,thresh=thresh,forest.type=FALSE)

#Constrain by heights
if(use.height.constrain==TRUE){
plot.fit<-height.constrain(plot.fit,data,level="Plot",dbh=dbh,max.size=max.size,na.too.big=na.too.big)
clust.fit<-height.constrain(clust.fit,data,level="Cluster",dbh=dbh,max.size=max.size,na.too.big=na.too.big)
bio.fit<-height.constrain(bio.fit,data,level="BioR",dbh=dbh,max.size=max.size,na.too.big=na.too.big, weib.only=TRUE)
cont.fit<-height.constrain(cont.fit,data,level="Continent",dbh=dbh,max.size=max.size,na.too.big=na.too.big, weib.only=TRUE)
clust.all.fit<-height.constrain(clust.all.fit,data,level="Cluster",dbh=dbh,max.size=max.size,na.too.big=na.too.big,forest.type=FALSE)
bio.all.fit<-height.constrain(bio.all.fit,data,level="BioR",dbh=dbh,max.size=max.size,na.too.big=na.too.big,forest.type=FALSE, weib.only=TRUE)
cont.all.fit<-height.constrain(cont.all.fit,data,level="Continent",dbh=dbh,max.size=max.size,na.too.big=na.too.big,forest.type=FALSE, weib.only=TRUE)
}

#Match in best parameters
data$a_Plot<-plot.fit[match(data$PlotID,plot.fit$ID),"Best_a"]
data$b_Plot<-plot.fit[match(data$PlotID,plot.fit$ID),"Best_b"]
data$c_Plot<-plot.fit[match(data$PlotID,plot.fit$ID),"Best_c"]
data$Plot.BestMod<-plot.fit[match(data$PlotID,plot.fit$ID),"BestMod"]

data$a_ClustF<-clust.fit[match(data$ID.clust,clust.fit$ID),"Best_a"]
data$b_ClustF<-clust.fit[match(data$ID.clust,clust.fit$ID),"Best_b"]
data$c_ClustF<-clust.fit[match(data$ID.clust,clust.fit$ID),"Best_c"]
data$ClustF.BestMod<-clust.fit[match(data$ID.clust,clust.fit$ID),"BestMod"]


data$a_BioRF<-bio.fit[match(data$ID.bio,bio.fit$ID),"Best_a"]
data$b_BioRF<-bio.fit[match(data$ID.bio,bio.fit$ID),"Best_b"]
data$c_BioRF<-bio.fit[match(data$ID.bio,bio.fit$ID),"Best_c"]
data$BioRF.BestMod<-bio.fit[match(data$ID.bio,bio.fit$ID),"BestMod"]

data$a_ContF<-cont.fit[match(data$ID.cont,cont.fit$ID),"Best_a"]
data$b_ContF<-cont.fit[match(data$ID.cont,cont.fit$ID),"Best_b"]
data$c_ContF<-cont.fit[match(data$ID.cont,cont.fit$ID),"Best_c"]
data$ContF.BestMod<-cont.fit[match(data$ID.cont,cont.fit$ID),"BestMod"]


data$a_Clust<-clust.all.fit[match(data$ClusterID,clust.all.fit$ID),"Best_a"]
data$b_Clust<-clust.all.fit[match(data$ClusterID,clust.all.fit$ID),"Best_b"]
data$c_Clust<-clust.all.fit[match(data$ClusterID,clust.all.fit$ID),"Best_c"]
data$Clust.BestMod<-clust.all.fit[match(data$ClusterID,clust.all.fit$ID),"BestMod"]

data$a_BioR<-bio.all.fit[match(data$BiogeographicalRegionID,bio.all.fit$ID),"Best_a"]
data$b_BioR<-bio.all.fit[match(data$BiogeographicalRegionID,bio.all.fit$ID),"Best_b"]
data$c_BioR<-bio.all.fit[match(data$BiogeographicalRegionID,bio.all.fit$ID),"Best_c"]
data$BioR.BestMod<-bio.all.fit[match(data$BiogeographicalRegionID,bio.all.fit$ID),"BestMod"]


data$a_Cont<-cont.all.fit[match(data$Continent,cont.all.fit$ID),"Best_a"]
data$b_Cont<-cont.all.fit[match(data$Continent,cont.all.fit$ID),"Best_b"]
data$c_Cont<-cont.all.fit[match(data$Continent,cont.all.fit$ID),"Best_c"]
data$Cont.BestMod<-cont.all.fit[match(data$Continent,cont.all.fit$ID),"BestMod"]

#Find finest geographic scale matched. Note using a parameter to determine parameter presence, as c parameters are NA for log-log models
if(no.plot==FALSE){
data$a_Best<-ifelse(!is.na(data$a_Plot),data$a_Plot,
	ifelse(!is.na(data$a_ClustF),data$a_ClustF,
		ifelse(!is.na(data$a_BioRF),data$a_BioRF,
			ifelse(!is.na(data$a_ContF),data$a_ContF,
				ifelse(!is.na(data$a_Clust),data$a_Clust,
					ifelse(!is.na(data$a_BioR),data$a_BioR,data$a_Cont))))))

data$b_Best<-ifelse(!is.na(data$a_Plot),data$b_Plot,
	ifelse(!is.na(data$a_ClustF),data$b_ClustF,
		ifelse(!is.na(data$a_BioRF),data$b_BioRF,
			ifelse(!is.na(data$a_ContF),data$b_ContF,
				ifelse(!is.na(data$a_Clust),data$b_Clust,
					ifelse(!is.na(data$a_BioR),data$b_BioR,
				data$b_Cont))))))
data$c_Best<-ifelse(!is.na(data$a_Plot),data$c_Plot,
	ifelse(!is.na(data$a_ClustF),data$c_ClustF,
		ifelse(!is.na(data$a_BioRF),data$c_BioRF,
			ifelse(!is.na(data$a_ContF),data$c_ContF,
				ifelse(!is.na(data$a_Clust),data$c_Clust,
					ifelse(!is.na(data$a_BioR),data$c_BioR,
			data$c_Cont))))))

}else{

data$a_Best<-ifelse(!is.na(data$a_BioRF),data$a_BioRF,
			ifelse(!is.na(data$a_ContF),data$a_ContF,
					ifelse(!is.na(data$a_BioR),data$a_BioR,data$a_Cont)))

data$b_Best<-ifelse(!is.na(data$a_BioRF),data$b_BioRF,
			ifelse(!is.na(data$a_ContF),data$b_ContF,
					ifelse(!is.na(data$a_BioR),data$b_BioR,
				data$b_Cont)))
data$c_Best<-ifelse(!is.na(data$a_BioRF),data$c_BioRF,
			ifelse(!is.na(data$a_ContF),data$c_ContF,
					ifelse(!is.na(data$a_BioR),data$c_BioR,
			data$c_Cont)))
}

return(list(data,plot.fit,clust.fit,bio.fit,cont.fit,cont.all.fit))
}

