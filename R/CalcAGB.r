#' Estimate aboveground biomass of trees
#'
#' \code{CalcAGB} estimates the aboveground biomass of each tree in the merged ForestPlots output dataset using the selected allometric equation.
#' Height is estimated using height-diameter allometries using either Feldpausch et al. (2011) regional Weibull equation parameters,
#' or parameters and allometric equation methods supplied using the height.data argument. This includes the option to use measured
#' heights. See \code{\link{height.mod}} for details of implemented height-diameter equations and how to set up user defined ones.
#' The function adds columns  to the dataset with the biomass information for all alive trees.
#' This function needs a dataset with the following information: PlotViewID, PlotID, TreeID, CensusNo, Diameter (DBH1-DBH4), Wood density (WD) and Allometric RegionID.The function assumes that the diameter used is DBH4, unless other DBH is selected.
#' See ForestPlots.net documentation for more information.
#' @references
#' Feldpausch TR, Banin L, Phillips OL, Baker TR, Lewis SL et al. 2011. Height-diameter allometry of tropical forest trees. Biogeosciences 8 (5):1081-1106. doi:10.5194/bg-8-1081-2011.
#'
#' Chave J, Coomes DA, Jansen S, Lewis SL, Swenson NG, Zanne AE. 2009. Towards a worldwide wood economics spectrum. Ecology Letters 12(4): 351-366. http://dx.doi.org/10.1111/j.1461-0248.2009.01285.x
#'
#' Zanne AE, Lopez-Gonzalez G, Coomes DA et al. 2009. Data from: Towards a worldwide wood economics spectrum. Dryad Digital Repository. http://dx.doi.org/10.5061/dryad.234

#'
#' @param xdataset a dataset for estimating biomass
#' @param dbh Name of column containing tree diameter (DBH, in mm).
#' @param height.data Dataframe containing model type and parameters for height diameter model. These are matched on to the \code{xdataset} object at tree level. If NULL (default), then regional height-diameter equations (from Felpaush et al. 2012) are used.
#' @param AGBFun Allometric equation function to use when estimating AGB, for example \code{\link{AGBChv14}}.
#'
#' @export
#' @author Gabriela Lopez-Gonzalez, Martin Sullivan


CalcAGB <- function (xdataset, dbh = "D4",height.data=NULL,AGBFun=AGBChv14){
         cdf <- xdataset
         ## Clean file
         cdf <- CleaningCensusInfo(xdataset)
         # Get Weibull Parameters
     	if(is.null(height.data)){
     		WHP <- WeibullHeightParameters
    		cdf <- merge(cdf, WHP, by = "AllometricRegionID", all.x = TRUE)
		cdf$Model<-"Weibull"
    	 }else{
 		cdf<-merge(cdf,height.data,by="TreeID",all.x=TRUE)
     	}
         #Estimate height
         cdf$HtF <- ifelse(cdf$D1 > 0 | cdf$Alive == 1, height.mod(cdf[,dbh],data=cdf), NA)
         #Add dead and recruits when codes are improved
         dbh_d <- paste(dbh,"_D", sep="")
         cdf$Htd <- ifelse(cdf$CensusStemDied==cdf$Census.No, height.mod(cdf[,dbh_d],data=cdf), NA)

         # Calculate AGB by stem Alive type
         cdf$AGBind <- ifelse(cdf$D1>0 & cdf$Alive == 1 & (cdf$CensusStemDied>cdf$Census.No | is.na(cdf$IsSnapped)),
			AGBFun(d=cdf[,dbh],h=cdf$HtF,wd=cdf$WD),
                              NA)
         #cdf$AGBAl <-  ifelse(cdf$Alive == 1, cdf$AGBind, NA)
         #cdf$AGBRec <- ifelse(cdf$NewRecruit == 1, cdf$AGBind, NA)
         cdf$AGBDead <-ifelse(cdf$CensusStemDied==cdf$Census.No,
			AGBFun(d=cdf[,dbh_d],h=cdf$Htd,wd=cdf$WD), NA)

         cdf

 }

