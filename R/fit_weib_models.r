#' @title fit.weib.models
#' @description Function to fit Weibull and log-log height-diameter models.
#' @param x Dataframe with columns names Height (consisting of measured heights), Diameter (measured diameters) and WD (wood density). These will be present if an object return by \code{\link{mergefp}} is used.
#' @param thresh Minimum number of height observations needed to attempt to fit models
#' @param weib.only Logical. Should the search for the best model be constrained to Weibull models. If TRUE, log-log models are only selected when neither Weibull models converge. Default is FALSE.
#' @return Dataframe with estimated parameters, number of height measurments used to fit the models and the error at estimating AGB compared to AGB estimates using observed heights (in units of Mg per tree). Parameters of the 'best' model (i.e. lowest AGB prediction error) are also given.
#' @author Martin Sullivan

#' @export

fit.weib.models<-function(x,thresh=20,weib.only=FALSE){
x<-x[,c("Height","Diameter","WD")]
ntrees<-nrow(x)
if(ntrees<thresh){
wf.param<-c("NA","NA","NA")
wfw.param<-c("NA","NA","NA")
lf.param<-c("NA","NA","NA")
result<-rep(NA,17)
result[14]<-ntrees
names(result)<-c("a_Weib","b_Weib","c_Weib","a_WeightWeib","b_WeightWeib","c_WeightWeib","a_Loglog","b_Loglog","c_Loglog","BestMod","Weibull_error","Weibull.weight_error","Loglog_error","Ntrees","Best_a","Best_b","Best_c")
}else{

wf<-try(nls(Height ~ a*(1-exp(-b*(Diameter/10)^c)),data=x,
                  na.action=na.omit,
                  start = c(a=25, b= 0.05, c= 0.7)),silent=TRUE)
if(class(wf)=="try-error"){
x$wf.height<-NA
wf.param<-c("NA","NA","NA")
}else{
x$wf.height<-predict(wf)
wf.param<-coef(wf)
}

weight.vec<-(pi*(x$Diameter/2000)^2)+1
wfw<-try(nls (Height ~ a*(1-exp(-b*(Diameter/10)^c)),
                  data=x,weights=weight.vec,
                  na.action=na.omit,
                  start = c(a=25, b= 0.05, c= 0.7)),silent=TRUE)
if(class(wfw)=="try-error"){
x$wfw.height<-NA
wfw.param<-c("NA","NA","NA")
}else{
x$wfw.height<-predict(wfw)
wfw.param<-coef(wfw)
}

lf<-lm(log(Height)~log(Diameter/10),data=x)
x$lf.height<-exp(predict(lf))
lf.param<-c(coef(lf),NA)

wf.agb<-sum(0.0673*(x$WD*(x$Diameter/10)^2*x$wf.height)^0.976/1000,na.rm=T)
wfw.agb<-sum(0.0673*(x$WD*(x$Diameter/10)^2*x$wfw.height)^0.976/1000,na.rm=T)
lf.agb<-sum(0.0673*(x$WD*(x$Diameter/10)^2*x$lf.height)^0.976/1000,na.rm=T)
obs.agb<-sum(0.0673*(x$WD*(x$Diameter/10)^2*x$Height)^0.976/1000,na.rm=T)

#Sum of all NAs=0
if(wf.agb==0)wf.agb<-NA
if(wfw.agb==0)wfw.agb<-NA

wf.error<-wf.agb-obs.agb
wfw.error<-wfw.agb-obs.agb
lf.error<-lf.agb-obs.agb

error<-c(wf.error,wfw.error,lf.error)
names(error)<-c("Weibull","Weibull.weight","Loglog")
abs.error<-abs(error)
best.mod<-as.numeric(grep(min(abs.error,na.rm=T),abs.error))
if(weib.only==TRUE){
if(!is.na(abs.error[1])|!is.na(abs.error[2])){
best.mod<-as.numeric(grep(min(abs.error[1:2],na.rm=T),abs.error[1:2]))
}
}
result<-as.numeric(c(wf.param,wfw.param,lf.param,best.mod,error/ntrees,ntrees))
names(result)<-c("a_Weib","b_Weib","c_Weib","a_WeightWeib","b_WeightWeib","c_WeightWeib","a_Loglog","b_Loglog","c_Loglog","BestMod",paste(names(error),"error",sep="_"),"Ntrees")
best_a<-ifelse(best.mod==1,result[1],ifelse(best.mod==2,result[4],result[7]))
best_b<-ifelse(best.mod==1,result[2],ifelse(best.mod==2,result[5],result[8]))
best_c<-ifelse(best.mod==1,result[3],ifelse(best.mod==2,result[6],result[9]))
best<-c(best_a,best_b,best_c)
names(best)<-c("Best_a","Best_b","Best_c")
result<-c(result,best)
}
return(result)
}
