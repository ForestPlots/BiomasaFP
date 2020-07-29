#'	Biomass gains and losses incorporating Kohyama et al. 2019 census interval correction
#'
#'	Function to estimate AGWP and losses following equations 18 and 19 of Kohyama et al. 2019 Forest Ecology and Management. Acts as a wrapper for \code{turnover}.
#' @param xdataset Object returned by \code{mergefp}.
#' @param AGBEquation Allometric equation function to use when estimating AGB, for example \code{\link{AGBChv14}}.
#' @param dbh Name of column containing diameter data. Default is "D4".
#' @param height.data Dataframe containing model type and parameters for height diameter model. These are matched on to the xdataset object at tree level. If NULL (default), then regional height-diameter equations (from Felpaush et al. 2012) are used.
#' @param palm.eq Logical. If TRUE the family level diameter based equation from Goodman et al 2013 is used to estimate AGB of monocots. Otherwise AGB of monocots is estimated using the allometric equation supplied to AGBEquation.
#' @param treefern Logical. If TRUE the height based equation from Tiepolo et al. 2002 is used to estimate the biomass of treeferns. If FALSE biomass is esimated using the allometric equation supplied to AGBEquation.
#' @return Dataframe with AGWP and biomass losses for each census interval. NA for census 1. Both AGWP and losses incorporate unobserved growth, and are in units of Mg/ha/year.
#' @author Martin Sullivan. Turnover function by T.S. Kohyama, T.I. Kohyama and D. Sheil.
#' @references Kohyama et al. 2019 Estimating net biomass production and loss from repeated measurements of trees in forests and woodlands: Formulae, biases and recomendations. Forest Ecology and Management 433: 729-740.


SummaryAGWP_Kohyama<-function(xdataset, AGBEquation, dbh ="D4",height.data=NULL,palm.eq=TRUE,treefern=TRUE){
if(all.equal(AGBEquation,BA)==TRUE & palm.eq==TRUE){
stop("Palm equation set to TRUE when calculating basal area. Please set palm.eq to FALSE")
}

AGBData<-CalcAGB (xdataset, dbh,height.data=height.data,AGBFun=AGBEquation)
 	if(palm.eq==TRUE){
 		AGBData[AGBData$Monocot==1,]$AGBind<-GoodmanPalm(AGBData[AGBData$Monocot==1,dbh])
 		AGBData[AGBData$Monocot==1,]$AGBDead<-GoodmanPalm(AGBData[AGBData$Monocot==1,paste(dbh,"_D",sep="")])
 	}
	if(treefern==TRUE){
	  AGBData[AGBData$Family=="Cyatheaceae",]$AGBind<-TiepoloTreeFern(h=AGBData[AGBData$Family=="Cyatheaceae","HtF"])
	  AGBData[AGBData$Family=="Cyatheaceae",]$AGBDead<-TiepoloTreeFern(h=AGBData[AGBData$Family=="Cyatheaceae","Htd"])
	}

AGBData$Census.prev<-AGBData[match(paste(AGBData$TreeID,AGBData$Census.No-1),paste(AGBData$TreeID,AGBData$Census.No)),"Census.Mean.Date"]
AGBData$Delta.time<-AGBData$Census.Mean.Date-AGBData$Census.prev
AGBData$AGB.prev<-AGBData[match(paste(AGBData$TreeID,AGBData$Census.No-1),paste(AGBData$TreeID,AGBData$Census.No)),"AGBind"]

res<-unique(data.frame("PlotViewID"=AGBData$PlotViewID,"Census.No"=AGBData$Census.No))
res$AGWP.ha.yr<-NA
res$AGBMort.ha.yr<-NA

for(i in 1:nrow(res)){
pv<-res$PlotViewID[i]
cns<-res$Census.No[i]

tmp<-AGBData[AGBData$PlotViewID==pv & AGBData$Census.No==cns,]

Intvl <- mean(tmp$Delta.time,na.rm=T)
if(is.na(Intvl))next
Area<-tmp$PlotArea[1]
tmp$AGBind[is.na(tmp$AGBind)]<-0
tmp$AGB.prev[is.na(tmp$AGB.prev)]<-0
dsurv<-tmp[tmp$AGBind>0 & tmp$AGB.prev>0,]

# Select species more than 6 survivors
# ISSUE - HOW DO WE DEAL WITH INDETS - currently treated as seperate species
n = table(dsurv$Species)
sp_list = names(subset(n, n[] > 6))

# aggregate species not in sp_list
tmp$Species = ifelse(tmp$Species %in% sp_list, as.character(tmp$Species), 'Others')
tmp$Species = factor(tmp$Species, levels=c(sp_list, 'Others'))

#Initial biomass of survivors
tmp$xs <- with(tmp,ifelse(AGBind>0 & AGB.prev>0, AGB.prev,0))

koy = turnover(tmp$AGB.prev, tmp$AGBind, tmp$xs, intvl = Intvl, area = Area, subpop = tmp$Species)
res$AGWP.ha.yr[i]<-koy[[1]]
res$AGBMort.ha.yr[i]<-koy[[2]]
}
res$PlotCode<-xdataset[match(res$PlotViewID,xdataset$PlotViewID),"PlotCode"]
res
}



