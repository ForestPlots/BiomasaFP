
#' Read and merge data downloaded from ForestPlots.net
#'
#' \code{mergefp} merges the 3 datasets needed for estimating biomass: Census Data (trees),
#' Metadata (md), wood density from each individual tree (wd). The three datasets have to be donwloaded from ForestPlots.net.
#' This function can use data frames and as well as file paths as input objects.
#'
#' @param trees an individual tree data file downloaded from the "Advanced Search". Can be a dataframe or a character indicating the file path.
#' @param md a metadata file downloaded from the "Query Library". Can be a dataframe or a character indicating the file path.
#' @param wd a wood density file, with the wood density for each individual tree by PlotViewID downloaded from the "Query Library". Can be a dataframe or a character indicating the file path.
#'
#' @author Gabriela Lopez-Gonzalez, Tim Baker, Martin Sullivan
#' @references Lopez-Gonzalez G, Lewis SL, Burkitt M. and Phillips OL. 2011. ForestPlots.net: a web application and research tool to manage and analyse tropical forest plot data. Journal of Vegetation Science 22: 610–613. doi: 10.1111/j.1654-1103.2011.01312.x
#'
#' Chave J, Coomes DA, Jansen S, Lewis SL, Swenson NG, Zanne AE. 2009. Towards a worldwide wood economics spectrum. Ecology Letters 12(4): 351-366. http://dx.doi.org/10.1111/j.1461-0248.2009.01285.x
#'
#' Zanne AE, Lopez-Gonzalez G, Coomes DA et al. 2009. Data from: Towards a worldwide wood economics spectrum. Dryad Digital Repository. http://dx.doi.org/10.5061/dryad.234
#'
#'
#' @export




mergefp <- function (trees,md,wd) {

        #1 read datasets (Tim can add comments)

        if(is.data.frame(trees)){
                CD<-trees
        }else{
                CD <- read.csv(trees, header = TRUE, stringsAsFactors = FALSE)
        }
        if(is.data.frame(md)){
                MT<-md
        }else{
                MT <- read.csv(md, header = TRUE, stringsAsFactors = FALSE)
        }
        if(is.data.frame(wd)){
                W<-wd
        }else{
                W <- read.csv(wd, header = TRUE, stringsAsFactors = FALSE)
        }


        #2 clean census data file

        #find lines that are in wrong format: blank or with text
        CD$PlotID<- as.numeric(CD$PlotID)
        RNAS<-CD[!is.na(CD$PlotID),]


        #Set as numeric values that need to be numeric (Tree Census information)
        RNAS$Census.Mean.Date <- as.numeric(RNAS$Census.Date)# Mark changed the header
        RNAS$Census.No <- as.numeric(RNAS$Census.No)
        RNAS$PlotID <- as.numeric(RNAS$PlotID)
        RNAS$PlotViewID <- as.numeric(RNAS$PlotViewID)
        RNAS$TreeID <- as.numeric(RNAS$TreeID)
        RNAS$FamilyAPGID <- as.numeric(RNAS$FamilyAPGID)
        RNAS$GenusID <- as.numeric(RNAS$GenusID)
        RNAS$SpeciesID <- as.numeric(RNAS$SpeciesID)
        RNAS$D1<- as.numeric(gsub("=","",RNAS$D1))
        RNAS$DPOMtMinus1<- as.numeric(gsub("=","",RNAS$DPOMtMinus1))
        RNAS$D2<- as.numeric(gsub("=","",RNAS$D2))
        RNAS$D3<- as.numeric(gsub("=","",RNAS$D3))
        RNAS$D4<- as.numeric(gsub("=","",RNAS$D4))
        RNAS$POM<- as.numeric(gsub("=","",RNAS$POM))
        RNAS$F1<- (gsub("=","",RNAS$F1))
        RNAS$F2<- (gsub("=","",RNAS$F2))
        RNAS$F3<- (gsub("=","",RNAS$F3))
        RNAS$F4<- (gsub("=","",RNAS$F4))
        RNAS$Height<-as.numeric(gsub("=","",RNAS$Height))
        RNAS$F5<- as.numeric(gsub("=","",RNAS$F5))
        RNAS$CI<-ifelse(RNAS$CI=="",NA, (gsub("=","",RNAS$CI )))
        RNAS$LI<-ifelse(RNAS$LI=="",NA, (gsub("=","",RNAS$LI )))



        #Remove Duplicated treeids due to vouchers; Select columns that will be used in the package
        RNASu <- RNAS[,c('PlotID','Plot.Code','Country','Census.Mean.Date',
                         'Census.No','PlotViewID','TreeID','FamilyAPGID',
                         'Family','GenusID','Genus','SpeciesID','Species',
                         'D1','DPOMtMinus1','D2','D3','D4','POM','F1','F2',
                         'F3','F4','Height','F5', 'CI','LI'
        )]
        #extract unique elements, Unique dataset for census data
        Udfile <-unique(RNASu)
        #Set as numeric values that need to be numeric (WD values)
        W$PlotID <- as.numeric(W$PlotID)
        W$PlotViewID <- as.numeric (W$PlotViewID)
        W$TreeID<- as.numeric (W$TreeID)
        W$WD <- as.numeric(gsub("=","", W$WD))
        #head(W)
        #Set as numeric values that need to be numeric (Metadata values)
        MT$PlotID <- as.numeric(gsub("=","",MT$PlotID))
        MT$PlotViewID <- as.numeric (MT$PlotViewID)
        MT$Altitude <- as.numeric (MT$Altitude)
        MT$LatitudeDecimal <- as.numeric (MT$LatitudeDecimal)
        MT$LongitudeDecimal <- as.numeric (MT$LongitudeDecimal)
        MT$PlotArea<- as.numeric(gsub("=","",MT$PlotArea))
        MT$AllometricRegionID <-as.numeric(gsub("=","",MT$AllometricRegionID))
        MT$AllometricRegionID <-as.numeric(gsub("NULL","NA",MT$AllometricRegionID))

        MT$ClusterID <- as.numeric(gsub("=","",MT$ClusterID))
        MT$ClusterID <- as.numeric(gsub("NULL","NA",MT$ClusterID))

        MT$ForestMoistureID <- as.numeric(gsub("=","",MT$ForestMoistureID))
        MT$ForestMoistureID <- as.numeric(gsub("NULL","NA",MT$ForestMoistureID))

        MT$ForestEdaphicID <- as.numeric(gsub("=","",MT$ForestEdaphicID))
        MT$ForestEdaphicID <- as.numeric(gsub("NULL","NA",MT$ForestEdaphicID))

        MT$ForestEdaphicHeightID <- as.numeric(gsub("=","",MT$ForestEdaphicHeightID))
        MT$ForestEdaphicHeightID <- as.numeric(gsub("NULL","NA",MT$ForestEdaphicHeightID))

        MT$ForestElevationID <- as.numeric(gsub("=","",MT$ForestElevationID))
        MT$ForestElevationID <- as.numeric(gsub("NULL","NA",MT$ForestElevationID))

        MT$ForestElevationHeightID <- as.numeric(gsub("=","",MT$ForestElevationHeightID))
        MT$ForestElevationHeightID <- as.numeric(gsub("NULL","NA",MT$ForestElevationHeightID))

        MT$BiogeographicalRegionID <- as.numeric(gsub("=","",MT$BiogeographicalRegionID))
        MT$BiogeographicalRegionID <- as.numeric(gsub("NULL","NA",MT$BiogeographicalRegionID))

        # Select MT columns to merge
        MT <- MT[, c('PlotID','PlotViewID','PlotCode','Altitude','LatitudeDecimal','LongitudeDecimal',
                     'PlotArea', 'ClusterName', 'ClusterID', 'ForestMoistureName',
                     'ForestMoistureID', 'ForestEdaphicName', 'ForestEdaphicID',
                     'ForestEdaphicHeightName',	'ForestEdaphicHeightID','ForestElevationName',
                     'ForestElevationID','ForestElevationHeightName', 'ForestElevationHeightID',
                     'BiogeographicalRegionName','BiogeographicalRegionID','RegionName','Continent','AllometricRegionID')]

        #head(MT)


        #merge information
        dataset <- merge(Udfile, MT, by = c('PlotViewID','PlotID'), x.all=TRUE)
        #head(dataset)
        datasetb <- merge(dataset, W, by = c('PlotViewID','TreeID','PlotID','PlotCode'), x.all=TRUE)
        #head(datasetb)
        #Select columns for mergefp output
        datasetc <- datasetb[,c('Continent', 'Country', 'PlotID', 'PlotCode'
                                ,'PlotViewID','LatitudeDecimal','LongitudeDecimal','Altitude','PlotArea',
                                'TreeID','FamilyAPGID','Family','GenusID','Genus','SpeciesID','Species', 'Census.No', 'Census.Mean.Date',
                                'D1','DPOMtMinus1','D2','D3','D4','POM','F1','F2','F3','F4','Height','F5','LI','CI','WD',
                                'AllometricRegionID',
                                'ClusterID',
                                'ForestMoistureID','ForestEdaphicID',
                                'ForestEdaphicHeightID',
                                'ForestElevationID','ForestElevationHeightID',
                                'BiogeographicalRegionID'
        )]
        #head(datasetc)

        # add Palm/monocot category  for records from the following families: 'Arecaceae','Strelitziaceae','Poaceae','Cyatheaceae'
        datasetc$Monocot <-  ifelse ( grepl('Arecaceae',datasetc$Family)==TRUE |
                                           grepl('Strelitziaceae',datasetc$Family)==TRUE |
                                           grepl('Poaceae',datasetc$Family)==TRUE |
                                           grepl('Cyatheaceae',datasetc$Family)==TRUE
                                   ,1,0)


        datasetc$PomChange <- ifelse (grepl('6',datasetc$F4)==TRUE,1,0)
        datasetc$Recruit   <- ifelse (grepl('n',datasetc$F1)==TRUE,1,0)

        datasetc

        #head(datasetc)
}

#' @rdname mergefp
#' @param dfmerged Merged dataframe returned by \code{mergefp}.


CleaningCensusInfo <- function (dfmerged) {

        dfm <- dfmerged
        #test using test merged file
        #dfm <- (datasetc)
        #head(dfm)

        #DeadorSnapped: First appearance of a tree that is dead or snapped
        # First gets all trees records that indicate a tree has died or is snapped
        #The following code can be removed once Alive and Snapped are added to the download


        dfm$Alive <- as.numeric (ifelse(dfm$F2==1, 1 ,0))

        dfm$Snapped <- as.numeric(ifelse( grepl("k",dfm$F1)== TRUE, 1, 0))
        #head(dfm) to test if Alive and Snapped have been aded correctly

        #Change Census.No to CensusNo once new format is implemented
        DeadorSnapped <-data.frame (dfm[(dfm$Alive == 0 | dfm$Snapped ==1) , c('PlotViewID', 'Census.No', 'TreeID',  'D1', 'Alive', 'Snapped')])
        #head(dfm)

        # if dead or snapped is empty
        if (nrow(DeadorSnapped)==0) {

                Clean <-dfm

                #head (CleanA)'CensusStemDied','Dead','D1_D','D2_D','D3_D','D4_D'
                Clean$CensusStemDied <- NA
                Clean$Dead <- NA
                Clean$D1_D <- NA
                Clean$D2_D <- NA
                Clean$D3_D <- NA
                Clean$D4_D <- NA
                Clean$IsSnapped <-NA## just added GL
                Clean
        }


        else

                # First time a dead or census tree  appears
                #nrow(DeadorSnapped)

        {
                DeadFirstCensusA <- data.frame( aggregate (Census.No ~ PlotViewID + TreeID, data = DeadorSnapped,  min ))



                # merge with #censusstemdied to retrieve if it is dead or if it is k

                DeadFirstCensusB<- merge(DeadFirstCensusA, DeadorSnapped,by= c('TreeID','PlotViewID','Census.No'), all.x = TRUE)
                DeadFirstCensusC<- DeadFirstCensusB[,c('PlotViewID','Census.No', 'TreeID','Snapped')]
                colnames (DeadFirstCensusC) <-c('PlotViewID','CensusNoDead', 'TreeID', 'IsSnapped')
                DeadFirstCensus <- DeadFirstCensusC

                #write.csv (DeadFirstCensus, file = 'dead.csv')
                #head(DeadFirstCensus)

                # Merged with PlotData for selecting censuses before tree died or snapped

                PlotDataDeadA <- merge(dfm, DeadFirstCensus, by= c('TreeID','PlotViewID'), all = TRUE)
                #head (PlotDataDeadA)
                #Merged with Plot Data and select previous Census of when trees died/snapped
                #Change Census.No to CensusNo once new format is implemented
                PlotDataDeadB <- PlotDataDeadA[PlotDataDeadA$Census.No == PlotDataDeadA$CensusNoDead-1 & !is.na(PlotDataDeadA$CensusNoDead), c('PlotViewID','CensusNoDead', 'TreeID','D1', 'D2','D3','D4','IsSnapped')]
                ## Added a Flat to indicate if tree was dead or k, if =1 then k

                PlotDataDeadC <- PlotDataDeadB[,c('PlotViewID','CensusNoDead', 'TreeID','D1', 'D2','D3','D4','IsSnapped')]

                colnames(PlotDataDeadC)<- c('PlotViewID', 'CensusNoDead', 'TreeID', 'D1_D','D2_D','D3_D','D4_D','IsSnapped')
                #head (PlotDataDeadB)
                #Clean Dataset
                #Change Census.No to CensusNo once new format is implemented

                CleanA <- merge (dfm, PlotDataDeadC, by= c('TreeID', 'PlotViewID'), all= TRUE)
                #head(CleanA)
                #      CleanA$Dead <- ifelse(CleanA$CensusNoDead==CleanA$Census.No,1, NA) Censuswhen tree/stem contributes to dead biomass (either as dead or as snapped)
                # Changed this statement to display the census when the tree died in each line, as then it can be used for ifelse estatement of AGB of individuals
                CleanA$CensusStemDied <- CleanA$CensusNoDead
		#Change no data value for CensusStemDied to 9999 - this means that trees that have not yet died with have CensusStemDied>CensusNo
                #GL: Removed this statment as somehow it was breaking the code
		#CleanA$CensusStemDied[is.na(CleanA$CensusStemDied)]<-9999
		#CleanA$CensusNoDead[is.na(CleanA$CensusNoDead)]<-9999
		#Trees that were never snapped are NA for IsSnapped - correct by changing to 0. GL: reverted to old statement as all AGB functions use is.NA snapped
		#CleanA$IsSnapped[is.na(CleanA$IsSnapped)]<-0
                CleanA$Dead <- ifelse(CleanA$F1==0 & CleanA$CensusNoDead==CleanA$Census.No,1,
                                      ifelse (CleanA$CensusNoDead>CleanA$Census.No,0, NA))
                CleanA$D1_D <- ifelse(CleanA$CensusNoDead==CleanA$Census.No,CleanA$D1_D, NA)
                CleanA$D2_D <- ifelse(CleanA$CensusNoDead==CleanA$Census.No,CleanA$D2_D, NA)
                CleanA$D3_D <- ifelse(CleanA$CensusNoDead==CleanA$Census.No,CleanA$D3_D, NA)
                CleanA$D4_D <- ifelse(CleanA$CensusNoDead==CleanA$Census.No,CleanA$D4_D, NA)

                # Alive status is corrected for trees that are back to live

                CleanA$Alive <- ifelse(CleanA$F2==1 & is.na(CleanA$IsSnapped),1,
                                        ifelse(CleanA$CensusNoDead>CleanA$Census.No  & CleanA$IsSnapped==0,1,
                                                ifelse(CleanA$CensusNoDead>CleanA$Census.No & CleanA$IsSnapped==1,1,
                                                        ifelse(CleanA$CensusNoDead==CleanA$Census.No & CleanA$IsSnapped==0,0,
                                                                ifelse(CleanA$CensusNoDead==CleanA$Census.No & CleanA$IsSnapped==1,1,
                                                                        ifelse(CleanA$CensusNoDead<CleanA$Census.No & CleanA$IsSnapped==1 & CleanA$F2==1, 1,
                                                                                ifelse(CleanA$CensusNoDead<CleanA$Census.No & CleanA$IsSnapped==1 & CleanA$F1==0,0,NA
                                                                                )))))))



                #head(CleanA)
                Clean <- CleanA[ , c( 'Continent', 'Country', 'PlotCode', 'PlotID','PlotViewID', 'LatitudeDecimal',
                                      'LongitudeDecimal', 'Altitude',	'PlotArea',
                                      'TreeID',	'FamilyAPGID',	'Family',
                                      'GenusID', 'Genus', 'SpeciesID',	'Species',
                                      'Census.No', 'Census.Mean.Date',
                                      'D1', 'DPOMtMinus1', 'D2', 'D3',	'D4', 'POM', 'F1', 'F2',
                                      'F3', 'F4', 'Height', 'F5', 'LI',	'CI',
                                      'WD', 'AllometricRegionID', 'ClusterID',
                                      'ForestMoistureID', 'ForestEdaphicID',
                                      'ForestEdaphicHeightID',	'ForestElevationID',
                                      'ForestElevationHeightID', 'BiogeographicalRegionID',
                                      'Monocot', 'PomChange',  'Alive','Snapped' ,'Recruit','CensusStemDied','Dead','D1_D','D2_D','D3_D','D4_D','IsSnapped'
                )]
                # Discuss if this version should be implemented with all the subplot and t1 information
                #Clean <- CleanA[ , c('TreeID','PlotViewID','ContinentName','CountryID','CountryName','AllometricRegionID','PlotID','PlotCode','PlotViewPlotCensusID','CensusNo','MeanDecimalDate', 'Subplot_Standard','x_standard','y_standard','SubPlotT1','SubPlotT2','x','y','FamilyAPGID','FamilyAPGName','GenusID','GenusName','SpeciesID','FullSpeciesName',        'TagNumber','DBH1','DPOMtMinus1',DBH2','DBH3','DBH4','POM','Flag1','Flag2','Flag3','Flag4','CI','LI','Alive','NewRecruit','POMChange','AliveNormal',  'MultipleStem',        'Snapped','WD','CensusStemDied','Dead','DBH1_D','DBH2_D','DBH3_D','DBH4_D','Altitude','LatitudeDecimal', 'LongitudeDecimal','PlotArea')]

                Clean
        }

}





