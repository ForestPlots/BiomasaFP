#' @title SummaryAGWP
#' @description Function to estimate AGWP by plotview and census interval, including estimating unobserved recruitment and growth of trees that died between census periods using the Talbot et al. (2014) census interval correction.
#' @param xdataset Object returned by \code{\link{mergefp}}.
#' @param AGBEquation Allometric equation function to use when estimating AGB, for example \code{\link{AGBChv14}}.
#' @param dbh Name of column containing diameter data. Default is "D4".
#' @param rec.meth Method used to estimate AGWP of recruits. If 0 (default), estimates growth from starting diameter of 0mm. If another value is provided, then growth is estimated from a starting diameter of 100mm.
#' @param height.data Dataframe containing model type and parameters for height diameter model. These are matched on to the xdataset object at tree level. If NULL (default), then regional height-diameter equations (from Felpaush et al. 2012) are used.
#' @param palm.eq Logical. If TRUE the family level diameter based equation from Goodman et al 2013 is used to estimate AGB of monocots. Otherwise AGB of monocots is estimated using the allometric equation supplied to AGBEquation.
#' @param treefern Logical. If TRUE the height based equation from Tiepolo et al. 2002 is used to estimate the biomass of treeferns. If FALSE biomass is esimated using the allometric equation supplied to AGBEquation.
#' @return A data frame with PlotViewID, CensusNo, and observed and unobserved elements of AGWP, stem dynamics and AGB mortality.
#' @author Martin Sullivan, Gabriela Lopez-Gonzalez
#' @references Talbot et al. 2014. Methods to estimate aboveground wood productivity from long-term forest inventory plots. Forest Ecology and Management 320: 30-38.

#' @export







SummaryAGWP <- function (xdataset, AGBEquation, dbh ="D4",rec.meth=0,height.data=NULL,palm.eq=TRUE,treefern=TRUE){
if(all.equal(AGBEquation,BA)==TRUE & palm.eq==TRUE){
stop("Palm equation set to TRUE when calculating basal area. Please set palm.eq to FALSE")
}
 	#Use AGBEquation to set Chv14 argument
	AGBData<-CalcAGB (xdataset, dbh,height.data=height.data,AGBFun=AGBEquation)

 	if(palm.eq==TRUE){
 		AGBData[AGBData$Monocot==1,]$AGBind<-GoodmanPalm(AGBData[AGBData$Monocot==1,dbh])
 		AGBData[AGBData$Monocot==1,]$AGBDead<-GoodmanPalm(AGBData[AGBData$Monocot==1,paste(dbh,"_D",sep="")])
 	}
	if(treefern==TRUE){
	  AGBData[AGBData$Family=="Cyatheaceae",]$AGBind<-TiepoloTreeFern(h=AGBData[AGBData$Family=="Cyatheaceae","HtF"])
	  AGBData[AGBData$Family=="Cyatheaceae",]$AGBDead<-TiepoloTreeFern(h=AGBData[AGBData$Family=="Cyatheaceae","Htd"])
	}
 	  IndAL <- aggregate (Alive/PlotArea ~ PlotViewID + Census.No,  data = AGBData, FUN=sum )
         AGBAlive <-aggregate (AGBind/PlotArea ~ PlotViewID + Census.No, data = AGBData, FUN=sum )
 	  AGBData$Census.prev<-AGBData[match(paste(AGBData$TreeID,AGBData$Census.No-1),paste(AGBData$TreeID,AGBData$Census.No)),"Census.Mean.Date"]
 	  AGBData$Delta.time<-AGBData$Census.Mean.Date-AGBData$Census.prev




 	  #Estimate AGB of recruits at 100mm
 	  Height100<-height.mod(dbh=100,data=AGBData)
	AGBData$AGB.100<-AGBEquation(d=100,h=Height100,wd=AGBData$WD)
 	if(palm.eq==TRUE){
 			if(length(AGBData[AGBData$Monocot==1,]$AGB.100)>0){
 				AGBData[AGBData$Monocot==1,]$AGB.100<-exp(-3.3488+(2.7483*log(10)))/1000
 			}
 	}
 	  AGBData$AGWPRec<-AGBData$AGBind-AGBData$AGB.100


 	  #Calculate growth of surviving trees
 	  AGBData$AGB.Prev<-AGBData[match(paste(AGBData$TreeID,AGBData$PlotViewID,AGBData$Census.No-1),paste(AGBData$TreeID,AGBData$PlotViewID,AGBData$Census.No)),"AGBind"]
 	  AGBData$Delta.AGB<-AGBData$AGBind-AGBData$AGB.Prev
 	  AGBData$Survive2<-ifelse(AGBData$Census.No>1 & AGBData$Recruit==0 & AGBData$Alive==1 & AGBData$Snapped==0,1,0)
 	  #Includes snapped trees - for stem dynamics
	  AGBData$Survive<-ifelse(AGBData$Census.No>1 & AGBData$Recruit==0 & AGBData$Alive==1 ,1,0)
 	  AGB.surv<-aggregate(Delta.AGB/PlotArea ~ PlotViewID + Census.No,  data = AGBData[AGBData$Survive2==1,], FUN=sum )
  	  IndSurv<-aggregate(Survive/PlotArea ~ PlotViewID + Census.No,  data = AGBData, FUN=sum )
 	  #Calculate mean WD per census
 	  WD<-aggregate(WD ~ PlotViewID + Census.No,  data = AGBData[AGBData$Alive==1,], FUN=function(x)mean(x,na.rm=T))


 	  # Find Recruits
         Recruits<- AGBData[AGBData$Recruit==1,]
         if(rec.meth==100){
 	AGBRec <- aggregate (AGWPRec/PlotArea ~ PlotViewID + Census.No,  data = Recruits, FUN=sum ) 	# NOTE: divides by PlotArea, so you get value per ha !!
 	}else{
 	AGBRec <- aggregate (AGBind/PlotArea ~ PlotViewID + Census.No,  data = Recruits, FUN=sum )
 	}
         colnames(AGBRec) <- c('PlotViewID','Census.No','AGBRec.PlotArea')
         IndRec <- aggregate (Recruit/PlotArea ~ PlotViewID + Census.No,  data = Recruits, FUN=sum )

         #merge recruit information
         Recs<- merge(AGBRec, IndRec, by =  c('PlotViewID','Census.No'), all.x=TRUE)

         # Dead stems  get only stems that are dead and have an AGB_D, to count only dead trees ones
         DeadTrees <-AGBData[AGBData$Dead==1 & !is.na(AGBData$D1_D),]
         #Match in size class growth rate
 	  IndDead <- aggregate (Dead/PlotArea ~ PlotViewID + Census.No,  data = DeadTrees, FUN=sum )
         AGBDe <-aggregate (AGBDead/PlotArea ~  PlotViewID + Census.No, data = AGBData, FUN=sum )
         #merge Dead trees
         Deads <- merge(AGBDe, IndDead, by = c('PlotViewID','Census.No'),all.x=TRUE)

 	  # Unobserved growth of dead trees
 	  growth.rate<-SizeClassGrowth(xdataset,dbh=dbh)
 	  dead2<-merge(DeadTrees,growth.rate,by="PlotViewID",all.x=T)
         dead2$DBH.death<-ifelse(dead2[,paste(dbh,"_D",sep="")]<200,(dead2$Delta.time/2)*dead2$Class1,
 		ifelse(dead2[,paste(dbh,"_D",sep="")]<400,(dead2$Delta.time/2)*dead2$Class2,
 			(dead2$Delta.time/2)*dead2$Class3))
 	  dead2$DBH.death<-dead2$DBH.death+dead2[,paste(dbh,"_D",sep="")]
 	  dead2<-dead2[!is.na(dead2$TreeID),]
 	  #Height at death
 	  dead2$Height.dead<-height.mod(dead2$DBH.death,data=dead2)
 	  #AGB at death
	  dead2$AGB.death2<-AGBEquation(d=dead2$DBH.death,h=dead2$Height.dead,wd=dead2$WD)

 	dead2<-dead2[!is.na(dead2$Monocot),]
 	 if(palm.eq==TRUE){
 		dead2$AGB.death2[dead2$Monocot==1]<-GoodmanPalm(dead2[dead2$Monocot==1,paste(dbh,"_D",sep="")])
 	 }
 	if(treefern==TRUE){
 	    dead2[dead2$Family=="Cyatheaceae",]$AGB.death2<-TiepoloTreeFern(h=dead2[dead2$Family=="Cyatheaceae","Height.dead"])
 	  }
 	 dead2$Unobs.dead<-dead2$AGB.death2-dead2$AGBDead
 	 unobs.dead<-aggregate(Unobs.dead/PlotArea~PlotViewID + Census.No, data = dead2, FUN=sum )


 	# merge all summaries
         mergeAGBAlive <- merge (AGBAlive, IndAL, by = c('PlotViewID','Census.No'),all.x=TRUE)
         mergeB <-  merge(mergeAGBAlive, Recs, by = c('PlotViewID','Census.No'),all.x=TRUE)
         mergeC <- merge (mergeB, Deads, by=  c('PlotViewID','Census.No'),all.x=TRUE)
 	  mergeD<-merge (mergeC,AGB.surv,  by = c('PlotViewID','Census.No'),all.x=TRUE)
         mergeE<-merge (mergeD,IndSurv,  by = c('PlotViewID','Census.No'),all.x=TRUE)
       	mergeF<-merge(mergeE,unobs.dead, by = c('PlotViewID','Census.No'),all.x=TRUE)
 	mergeF<-merge(mergeF,WD,by=c('PlotViewID','Census.No'),all.x=TRUE)
 	#Correct names of merge F to be legal
 	names(mergeF)<-sub("/",".",names(mergeF))	 						# NOTE: .PlotArea replaces /PlotArea, so actually values per ha (plotarea weighted)



 	#Estimate growth of unobserved recruits - 1st merge growth rate data with allometric parameter data
 	col.keep<-intersect(c("PlotViewID","Census.No","a_par","b_par","c_par","d_par","ChaveE","Height","Delta.time","Model"),names(AGBData))
 	unobs.dat<-unique(AGBData[!is.na(AGBData$Delta.time),col.keep])
 	unobs.dat<-unobs.dat[!duplicated(paste(unobs.dat$PlotViewID,unobs.dat$Census.No)),]
 	unobs.dat<-merge(unobs.dat,growth.rate,by="PlotViewID",all.x=T)
 	#Then merge with existing data to get number of recruits and deaths
 	tmp<-merge(mergeF,unobs.dat,by = c('PlotViewID','Census.No'),all.x=TRUE)
 	#Get number of stems in previous census
 	tmp$PrevStems<-tmp[match(paste(tmp$PlotViewID,tmp$Census.No-1),paste(tmp$PlotViewID,tmp$Census.No)),"Alive.PlotArea"]
 	#Estimate recruitment and mortality rates
 	stems<-tmp$PrevStems
 	death<-tmp$Dead.PlotArea
 	recruits<-tmp$Recruit.PlotArea
 	interval<-tmp$Delta.time
 	tmp$rec.rate<-recruits/stems/interval
 	tmp$mort.rate<-death/stems/interval
	#tmp$Model<-"Weibull" # Need to update
 	tmp$Model[is.na(tmp$Model)]<-"Weibull"
 	#tmp$rec.rate <- replace(tmp$rec.rate,is.na(tmp$rec.rate),0)													# CHANGE: added, otherwise function doesn't work as crashes on plotviews with 0 recruitment
 	#tmp$mort.rate <- replace(tmp$mort.rate,is.na(tmp$mort.rate),0)
 	#Get time weighted mean of rec and mort rates for each plot
 	a<-split(tmp,f=tmp$PlotViewID)
 	rec1<-unlist(lapply(a,function(x)weighted.mean(x$rec.rate,x$Delta.time,na.rm=T)))
 	mort1<-unlist(lapply(a,function(x)weighted.mean(x$mort.rate,x$Delta.time,na.rm=T)))
 	rec.means<-data.frame("PlotViewID"=names(rec1),"Rec.mean"=as.numeric(rec1),stringsAsFactors=F)
 	mort.means<-data.frame("PlotViewID"=names(mort1),"Mort.mean"=as.numeric(mort1),stringsAsFactors=F)
	plot.means<-merge(rec.means,mort.means,by="PlotViewID",all.x=T)
 	tmp<-merge(tmp,plot.means,by="PlotViewID",all.x=T)
 	#Growth of unobserved recruits
 	unobs.D<-tmp$Class1*tmp$Delta.time*(1/3)
 	#D at death
 	D.death<-100+unobs.D
 	#Get height at death
 	height.death<-height.mod(D.death,data=tmp)
 	#Height at point of recruitment
 	height.rec<-height.mod(100,data=tmp)
 	#WD_mean<-tmp$MeanWD # Commented out - mean across all censuses
 	WD_mean<-tmp$WD # Use mean in plot and census
 	#Estimate AGB at death
 	AGB.death<-AGBEquation(d=D.death,h=height.death,wd=WD_mean)
 	#Estimate AGB at point of recruitment
 	AGB.rec<-AGBEquation(d=100,h=height.rec,wd=WD_mean)
 	#Get number of unobserved recruits
 	unobs.rec<-tmp$Delta.time*tmp$Rec.mean*tmp$Mort.mean*tmp$PrevStems
 	#Get growth of unobs recruits
 	#Option for recruits starting at 100
 	if(rec.meth==100){
 		unobs.growth<-(unobs.rec*(AGB.death-AGB.rec))
 	}else{
 		unobs.growth<-(unobs.rec*AGB.death)
	}
 	tmp$UnobsGrowth<-unobs.growth
 	tmp$UnobsRec<-unobs.rec
 	#Simplyfy tmp and merge with MergeF
 	tmp<-tmp[,c("PlotViewID","Census.No","rec.rate","mort.rate","UnobsGrowth","UnobsRec","Delta.time")]
 	mergeG<-merge (mergeF,tmp,  by = c('PlotViewID','Census.No'),all.x=TRUE)


 	#Create output file
         SummaryB <- mergeG


         SummaryB <- SummaryB[order(SummaryB$PlotViewID, SummaryB$Census.No, decreasing=FALSE), ]
 	#Add AGWP per ha
 	SummaryB$AGWP.PlotArea<-with(SummaryB,ifelse(is.na(Delta.AGB.PlotArea),0,Delta.AGB.PlotArea)+ifelse(is.na(AGBRec.PlotArea),0,AGBRec.PlotArea)+ifelse(is.na(Unobs.dead.PlotArea),0,Unobs.dead.PlotArea)+ifelse(is.na(UnobsGrowth),0,UnobsGrowth)) 						# NOTE: .PlotArea means plotarea weighted so per ha

 	SummaryB$AGWP.PlotArea.Year<-SummaryB$AGWP.PlotArea/SummaryB$Delta.time

 	#Neaten up column names
 	names(SummaryB)<-c("PlotViewID","Census.No","AGB.ha","Stems.ha",
		"AGWPrec.ha","Recruit.ha","AGBmort.ha","Mortality.ha","AGWPsurv.ha","SurvivingStems.ha","UnobsAGWPmort.ha","Mean.WD",			# NOTE: from mergeF
		"Recruitment.stem.year","Mortality.stem.year","UnobsAGWPrec.ha","UnobsRecruits.ha","CensusInterval",						# NOTE: from tmp, all per ha as PlotArea weighted !
		"AGWP.ha","AGWP.ha.year")
	SummaryB$PlotCode<-xdataset[match(SummaryB$PlotViewID,xdataset$PlotViewID),"PlotCode"]

 	SummaryB
 }



