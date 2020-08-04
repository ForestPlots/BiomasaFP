#'	Biomass gains and losses incorporating Kohyama et al. 2020 census interval correction
#'
#'	Function to estimate AGWP, biomass loss and stem turnover at stand level by sumarising species level estimates of vital rates. Species-level method from Kohyama et al. 2020 Journal of Ecology. Acts as a wrapper for \code{productivity}.
#' @param xdataset Object returned by \code{mergefp}.
#' @param AGBEquation Allometric equation function to use when estimating AGB, for example \code{\link{AGBChv14}}.
#' @param dbh Name of column containing diameter data. Default is "D4".
#' @param height.data Dataframe containing model type and parameters for height diameter model. These are matched on to the xdataset object at tree level. If NULL (default), then regional height-diameter equations (from Felpaush et al. 2012) are used.
#' @param palm.eq Logical. If TRUE the family level diameter based equation from Goodman et al 2013 is used to estimate AGB of monocots. Otherwise AGB of monocots is estimated using the allometric equation supplied to AGBEquation.
#' @param treefern Logical. If TRUE the height based equation from Tiepolo et al. 2002 is used to estimate the biomass of treeferns. If FALSE biomass is esimated using the allometric equation supplied to AGBEquation.
#' @return Dataframe with AGWP and biomass losses for each census interval. NA for census 1. Both AGWP and losses incorporate unobserved growth, and are in units of Mg/ha/year.
#' @author Martin Sullivan. Turnover function by T.S. Kohyama et al.
#' @references Kohyama et al. 2019 Estimating net biomass production and loss from repeated measurements of trees in forests and woodlands: Formulae, biases and recomendations. Forest Ecology and Management 433: 729-740.


SummaryAGWP_Kohyama2<-function(xdataset, AGBEquation, dbh ="D4",height.data=NULL,palm.eq=TRUE,treefern=TRUE){
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
  AGBData$D.prev<-AGBData[match(paste(AGBData$TreeID,AGBData$Census.No-1),paste(AGBData$TreeID,AGBData$Census.No)),dbh]
  AGBData$D.cur<-AGBData[,dbh]

  res<-unique(data.frame("PlotViewID"=AGBData$PlotViewID,"Census.No"=AGBData$Census.No))
  res$AGWP.ha.yr<-NA
  res$AGBMort.ha.yr<-NA
  res$StemRec<-NA
  res$StemMort<-NA
  res$AGB.ha<-NA
  res$Stems.ha<-NA

  for(i in 1:nrow(res)){
    pv<-res$PlotViewID[i]
    cns<-res$Census.No[i]

    tmp<-AGBData[AGBData$PlotViewID==pv & AGBData$Census.No==cns,]

    Intvl <- mean(tmp$Delta.time,na.rm=T)
    if(is.na(Intvl))next
    Area<-tmp$PlotArea[1]
    tmp$AGBind[is.na(tmp$AGBind)]<-0
    tmp$AGB.prev[is.na(tmp$AGB.prev)]<-0

    # Select species with > 10 stems
    n = table(tmp$Species)
    sp_list = names(subset(n, n[] > 10))
    tmp$Species = ifelse(tmp$Species %in% sp_list, as.character(tmp$Species), 'Others')
    tmp$Species = factor(tmp$Species, levels=c(sp_list, 'Others'))

    tmp2<-split(tmp,f=tmp$Species)
    sbpop<-lapply(tmp2,function(x)productivity(x$D.prev,x$D.cur,x$AGB.prev,x$AGBind,Intvl,Area,x$Alive,x$Recruit))
    sbpop.sum<-as.data.frame(do.call(rbind,sbpop))
    loss<-sum(sbpop.sum$l*sbpop.sum$B,na.rm=T)
    prod<-sum(sbpop.sum$p*sbpop.sum$B,na.rm=T)
    stem.r<-weighted.mean(sbpop.sum$r,sbpop.sum$N,na.rm=T)
    stem.m<-weighted.mean(sbpop.sum$m,sbpop.sum$N,na.rm=T)
    biomass<-sum(sbpop.sum$B,na.rm=T)
    stems<-sum(sbpop.sum$N,na.rm=T)
    res$AGWP.ha.yr[i]<-prod
    res$AGBMort.ha.yr[i]<-loss
    res$StemRec[i]<-stem.r
    res$StemMort[i]<-stem.m
    res$AGB.ha[i]<-biomass
    res$Stems.ha[i]<-stems
  }
  res$PlotCode<-xdataset[match(res$PlotViewID,xdataset$PlotViewID),"PlotCode"]
  res
}



