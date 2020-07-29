#' Calculation of instantanious mortality and recruitment rates - used in SummaryAGWP2
#'
#' Accessory functions used by SummaryAGWP2.
#' \code{life.and.death} is called within \code{unobs.recs}.
#' @param time Census interval
#' @param mort.pc Per capita mortaility rate
#' @param rec.area Area-based recruitment rate


life.and.death<-function(time,mort.pc,rec.area){
if(is.na(time)){
unobs.rec<-NA
mn.life<-NA
Life.span.dead<-NA
}else{
dead<-c()
lspan<-c()
for(i in 1:(round(time)-1)){
dead2<-c()
for(z in 1:i){
dead2[z]<-rec.area-(((1-mort.pc)^z)*rec.area)
}
life<-data.frame("LS"=seq(1,length(dead2)),"Prop.dead"=c(dead2[1],diff(dead2)))
dead[i]<-rec.area-(((1-mort.pc)^i)*rec.area)
lspan[i]<-weighted.mean(life$LS,life$Prop.dead)
}
mn.life<-weighted.mean(lspan,dead)
unobs.rec<-sum(dead)
Life.span.dead=lspan[length(lspan)]
}
return(c("N.unobs.rec"=unobs.rec,"Life.span.unobs"=mn.life,"Life.span.dead"=Life.span.dead))
}

#' @rdname life.and.death
#' @param dat Dataframe returned by \code{CalcAGB}.

unobs.recs<-function(dat){
alive<-tapply(dat$Alive,list(dat$PlotViewID,dat$Census.No),sum)
alive2<-alive[,-1]
alive1<-alive[,-ncol(alive)]
recruit<-tapply(dat$Recruit,list(dat$PlotViewID,dat$Census.No),sum)
recruit<-recruit[,-1]
survive<-alive2-recruit
date<-tapply(dat$Census.Mean.Date,list(dat$PlotViewID,dat$Census.No),function(x)x[1])
date1<-date[,-ncol(date)]
date2<-date[,-1]
time<-date2-date1
if(is.matrix(time)==FALSE){
time<-as.matrix(time)
}
area<-tapply(dat$PlotArea,list(dat$PlotViewID,dat$Census.No),function(x)x[1])
area<-area[,-ncol(area)]

Ma<-(alive1/area)*(1-(survive/alive1)^(1/time))
Mc<-1-(survive/alive1)^(1/time)
Ra<-(Mc/area)*(alive2-survive)/(1-(1-Mc)^time)

unobs.rec<-matrix(nrow=nrow(time),ncol=ncol(time))
LS.unobs.rec<-matrix(nrow=nrow(time),ncol=ncol(time))
LS.dead<-matrix(nrow=nrow(time),ncol=ncol(time))
dimnames(unobs.rec)<-list(dimnames(alive)[[1]],seq(1:ncol(time)))
dimnames(LS.unobs.rec)<-list(dimnames(alive)[[1]],seq(1:ncol(time)))
dimnames(LS.dead)<-list(dimnames(alive)[[1]],seq(1:ncol(time)))

#With census weighted mean for mortality and recruitment
for(i in 1:nrow(time)){
for(z in 1:ncol(time)){
M1<-weighted.mean(Mc[i,],time[i,],na.rm=T)
R1<-weighted.mean(Ra[i,],time[i,],na.rm=T)
tmp<-life.and.death(time[i,z],M1,R1)
unobs.rec[i,z]<-tmp[1]
LS.unobs.rec[i,z]<-tmp[2]
LS.dead[i,z]<-tmp[3]
}
}
return(list(unobs.rec,LS.unobs.rec,LS.dead,time))
}
