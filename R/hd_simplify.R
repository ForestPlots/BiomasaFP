#' Select height-diameter allometry parameters
#'
#' Selects height-diameter allometry parameters at the desired level from the object returned by \code{\link{local.heights}}. Height-diameter allometric models are fitted at different levels of the dataset (e.g. cluster, biogeographic region) and can be either selected as the finest-scale available (\code{param.type}="Best", default) or at a specified level
#'
#' @param height.data Object returned by \code{\link{local.heights}}
#' @param param.type Level at which to return height-diameter parameters. One of "Best" (defualt), "ClusterF","BioRF" or "ContinentF

hd.simplify<-function(height.data,param.type="Best"){
		p1<-paste("a",param.type,sep="_")
 		p2<-paste("b",param.type,sep="_")
 		p3<-paste("c",param.type,sep="_")
 		hd<-height.data[,c("TreeID",p1,p2,p3)]
 		hd<-unique(hd)
 		names(hd)<-c("TreeID","a_par","b_par","c_par")
		hd$Model<-ifelse(is.na(hd$c_par),"Loglog","Weibull")
		hd
}


