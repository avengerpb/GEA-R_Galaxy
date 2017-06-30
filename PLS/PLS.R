
PLSR<-function(datos,file_nameCOV,ExpDes=ExpDes,BiplotFormat=BiplotFormat,traits=traits,selected=selected, typecov=typecov, tmp=tmp){

rm(list=ls(all=TRUE))
library(lsmeans)
library(reshape)
library(nlme)
library(plsgenomics)
##############################################################################################################
#InputRespFile(crear una subset llamado raw con estas variables: ENV,GEN,REP,BLOCK,variable)
#InputCovFile(crear un subset con todas las covariables a usar)
#BiplotFormat (PDF,PNG,WINMETA)
orgdir=getwd()

for (b in 1:length(traits)){ 
  raw=datos
  index <- which(as.character(raw$Loc)%in%selected)
  if(length(index)>0)	raw <- raw[index,]
#Modify so that it read file directly from uploaded input  
  if (typecov=="env"){
  Covariate=file_nameCOV
  Covariate=Covariate[order(Covariate$ENV),]
  indexCOV <- which(as.character(Covariate$ENV)%in%selected)
  if(length(indexCOV)>0)	Covariate <- Covariate[indexCOV,]
  }
  
  if (typecov=="gen" || typecov=="markers"){
  Covariate=file_nameCOV
  Covariate=Covariate[order(Covariate$GEN),]
  }
#Modify Loc=ENV,Entry=GEN,Rep=REP,Block=BLOCK
  if (ExpDes=="LSMeans"){
    raw=as.data.frame(cbind(ENV=raw$ENV,GEN=raw$GEN,YLD=as.numeric(as.character(raw[,traits[b]]))))
    raw$GEN=as.factor(raw$GEN)
    raw$ENV=as.factor(raw$ENV)
    raw$YLD=as.numeric(as.character(raw$YLD))
    varname=traits[b]
    colnames(raw)[3]=c("YLD")
  }
  
  if (ExpDes=="RCB"){
    raw=as.data.frame(cbind(ENV=raw$ENV,GEN=raw$GEN,REP=raw$REP,YLD=as.numeric(as.character(raw[,traits[b]]))))
    raw$GEN=as.factor(raw$GEN)
    raw$ENV=as.factor(raw$ENV)
    raw$REP=as.factor(raw$REP)
    raw$YLD=as.numeric(as.character(raw$YLD))
    varname=traits[b]
    colnames(raw)[4]=c("YLD")
  }
  
  if (ExpDes=="Lattice"){
    raw=as.data.frame(cbind(ENV=raw$ENV,GEN=raw$GEN,REP=raw$REP,BLOCK=raw$BLOCK,YLD=as.numeric(as.character(raw[,traits[b]]))))
    raw$GEN=as.factor(raw$GEN)
    raw$ENV=as.factor(raw$ENV)
    raw$REP=as.factor(raw$REP)
    raw$BLOCK=as.factor(raw$BLOCK)
    raw$YLD=as.numeric(as.character(raw$YLD))
    varname=traits[b]
    colnames(raw)[5]=c("YLD")
  }
  
  if (ExpDes=="RCB"){
    
    #Analysis of Varianza Using All Data
    ammi=lm(YLD ~ REP+ENV+GEN,data=raw,na.action=na.omit )
    res=ammi$residuals
    res=as.data.frame(cbind(raw$GEN,raw$ENV,raw$REP,res))
    colnames(res)=c("GEN","ENV","REP","RESIDUAL")
	if (typecov=="env") {outres=cast(res,ENV~GEN,mean,value="RESIDUAL")       }
	if (typecov=="gen" || typecov=="markers") {outres=cast(res,GEN~ENV,mean,value="RESIDUAL")       }
  }
  
  if (ExpDes=="Lattice"){
    
	if (typecov=="env"){
    raw$ENV=as.factor(as.numeric(raw$ENV))
    covparms=list()
    lsgen=list()
  
  for (i in 1:length(levels(raw$ENV))){
    raw2<<-raw[raw$ENV==levels(raw$ENV)[i],]
    ver=lme(YLD ~ REP+GEN,data=raw2,random=~1|BLOCK/REP,na.action=na.omit)  
    lsgen[[i]]=cbind(ENV=rep(levels(raw$ENV)[i],length(levels(raw$GEN))),as.matrix(summary(lsmeans(ver,"GEN")))[,1:2])
    
  }
  
  lsgeng=lsgen[[1]]
  for (i in 2:length(levels(raw$ENV))){
    lsgeng=rbind(lsgeng,lsgen[[i]])
  }
  
  #Adjusted Means Using Proc Mixed for AMMI
  meansa=as.data.frame(cbind(as.numeric(lsgeng[,1]),as.numeric(lsgeng[,2]),as.numeric(lsgeng[,3])))
  colnames(meansa)=c("ENV","GEN","YLD")
  meansa$ENV=as.factor(meansa$ENV)
  meansa$GEN=as.factor(meansa$GEN)
  
  envgeno=lm(YLD ~ ENV+GEN,data=meansa,na.action=na.omit)
  res=envgeno$residuals
  
  res=as.data.frame(cbind(as.numeric(lsgeng[,2]),as.numeric(lsgeng[,1]),round(res,5)))
  colnames(res)=c("GEN","ENV","RESIDUAL")
  res$GEN=as.factor(res$GEN)
  res$ENV=as.factor(res$ENV)
  res$RESIDUAL=as.numeric(as.character(res$RESIDUAL))
  outres=cast(res,ENV~GEN,mean,value="RESIDUAL")
  }
  
  if(typecov=="gen" || typecov=="markers"){
  raw$GEN=as.factor(as.numeric(raw$GEN))
  covparms=list()
  lsgen=list()

for (i in 1:length(levels(raw$GEN))){
  raw2<<-raw[raw$GEN==levels(raw$GEN)[i],]
  ver=lme(YLD ~ REP+ENV,data=raw2,random=~1|BLOCK/REP,na.action=na.omit)  
  lsgen[[i]]=cbind(GEN=rep(levels(raw$GEN)[i],length(levels(raw$ENV))),as.matrix(summary(lsmeans(ver,"ENV")))[,1:2])
  
}

lsgeng=lsgen[[1]]
for (i in 2:length(levels(raw$GEN))){
  lsgeng=rbind(lsgeng,lsgen[[i]])
}

#Adjusted Means Using Proc Mixed for AMMI
meansa=as.data.frame(cbind(as.numeric(lsgeng[,1]),as.numeric(lsgeng[,2]),as.numeric(lsgeng[,3])))
colnames(meansa)=c("GEN","ENV","YLD")
meansa$ENV=as.factor(meansa$ENV)
meansa$GEN=as.factor(meansa$GEN)

envgeno=lm(YLD ~ ENV+GEN,data=meansa,na.action=na.omit)
res=envgeno$residuals

res=as.data.frame(cbind(as.numeric(lsgeng[,2]),as.numeric(lsgeng[,1]),round(res,5)))
colnames(res)=c("ENV","GEN","RESIDUAL")
res$GEN=as.factor(res$GEN)
res$ENV=as.factor(res$ENV)
res$RESIDUAL=as.numeric(as.character(res$RESIDUAL))
outres=cast(res,GEN~ENV,mean,value="RESIDUAL")
  }
  }
  
  if (ExpDes=="LSMeans"){
    ammi=lm(YLD ~ ENV+GEN,data=raw,na.action=na.omit )
    res=ammi$residuals
    res=as.data.frame(cbind(raw$GEN,raw$ENV,res))
    colnames(res)=c("GEN","ENV","RESIDUAL")
    if (typecov=="env") {outres=cast(res,ENV~GEN,mean,value="RESIDUAL")       }
	if (typecov=="gen" || typecov=="markers") {outres=cast(res,GEN~ENV,mean,value="RESIDUAL")}       
  }
  
  ##For ALL  
  nGen=length(levels(raw$GEN))
  nEnv=length(levels(raw$ENV))
  nCov=dim(Covariate)[2]-1
  
  if (typecov=="env"){
  outres2=as.matrix(scale(outres,center=T,scale=T)[1:nEnv,1:nGen])
  cov2=as.matrix(scale(Covariate[-1],center=T,scale=T)[1:nEnv,1:nCov])
  }  
  if (typecov=="gen" || typecov=="markers"){
  outres2=as.matrix(scale(outres,center=T,scale=T)[1:nGen,1:nEnv])
  cov2=as.matrix(Covariate[-1])[1:nGen,1:nCov]
  }
  
  anapls=pls.regression(cov2,outres2,Xtest=cov2, unit.weights=T)
  anaplsP=anapls$P[,1:2]
  colnames(anaplsP)=c("X-loadingsf1","X-loadingsf2")
  rownames(anaplsP)=colnames(cov2)
  anaplsP=anaplsP/max(abs(anaplsP))
  anaplsP[,2]=(-1)*anaplsP[,2]
  
  anaplsQ=anapls$Q[,1:2]
  colnames(anaplsQ)=c("Y-loadingsf1","Y-loadingsf2")
  rownames(anaplsQ)=colnames(outres2)
  anaplsQ=anaplsQ/max(abs(anaplsQ))
  anaplsQ[,2]=(-1)*anaplsQ[,2]
  
  anaplsT=anapls$T[,1:2]
  colnames(anaplsT)=c("X-Scoresf1","X-Scoresf2")
  rownames(anaplsT)=as.character(Covariate[,1])
  anaplsT=anaplsT/max(abs(anaplsT))
  anaplsT[,2]=(-1)*anaplsT[,2]
  
  cbind.fill<-function(...){
    nm <- list(...) 
    nm<-lapply(nm, as.matrix)
    n <- max(sapply(nm, nrow)) 
    do.call(cbind, lapply(nm, function (x) 
      rbind(x, matrix(, n-nrow(x), ncol(x))))) 
  }
  
  tot=cbind.fill(colnames(cov2),round(anaplsP,5),colnames(outres2),round(anaplsQ,5),as.character(Covariate[,1]),round(anaplsT,5))
  rownames(tot)=NULL
  if (typecov=="env"){
  colnames(tot)[1]="Covariates"
  colnames(tot)[4]="GEN"
  colnames(tot)[7]="ENV"
  }
  if (typecov=="gen" || typecov=="markers"){
  colnames(tot)[1]="Covariates"
  colnames(tot)[4]="ENV"
  colnames(tot)[7]="GEN"
  }
  
  exis=tcrossprod(anapls$T,anapls$P)
  f1=(t(anapls$T[,1])%*%anapls$T[,1]%*%(t(anapls$P[,1])%*%anapls$P[,1]))/sum(diag(t(exis)%*%exis))*100
  f1=format(f1,digits=4)
  f2=(t(anapls$T[,2])%*%anapls$T[,2]%*%(t(anapls$P[,2])%*%anapls$P[,2]))/sum(diag(t(exis)%*%exis))*100
  f2=format(f2,digits=4)
  
  nop=which((anaplsP[,1]<=0.5 & anaplsP[,1]>=(-0.5)) & (anaplsP[,2]<=0.5 & anaplsP[,2]>=(-0.5)))
  
  newfold=paste("OutputPLS&",paste(tmp,paste(ExpDes,paste("&",traits[b],sep=""),sep=""),sep=""),sep="")
  sal=paste(orgdir,paste("/",newfold,sep=""),sep="")
  dir.create(sal)
  setwd(sal)
  
  write.csv(tot,paste("Final Scores PLS",paste(traits[b],paste(ExpDes,".csv"))))
  
if (BiplotFormat=="pdf"){pdf(paste("Biplot for PLS",paste(traits[b],paste(ExpDes,".pdf"))))}
if (BiplotFormat=="png"){png(paste("Biplot for PLS",paste(traits[b],paste(ExpDes,".png"))))}
if (BiplotFormat=="wmf"){win.metafile(paste("Biplot for PLS",paste(traits[b],paste(ExpDes,".wmf"))))}
 
if (typecov=="env"){ 
biplot(anaplsQ,anaplsT,cex = 0.8,xlab=paste("factor1 (",paste(f1,"%)",sep=""),sep=""),
       ylab=paste("factor2 (",paste(f2,"%)",sep=""),sep=""),main="Biplot for PLS using mixed models",
       xlim=c(-1,1),ylim=c(-1,1),col.main="blue",col.lab="blue")	 
if(length(nop)!=0){	   
text(anaplsP[-nop,1],anaplsP[-nop,2],rownames(anaplsP)[-nop],col="green",cex=0.6)
}else{text(anaplsP[,1],anaplsP[,2],rownames(anaplsP),col="green",cex=0.6)}
abline(h=0)
abline(v=0)}
#yaxt="n",xaxt="n",  esta instruccion la quite para que salieran los valores de los ejes
if (typecov=="gen" || typecov=="markers"){ 
biplot(anaplsT,anaplsQ,cex = 0.8,xlab=paste("factor1 (",paste(f1,"%)",sep=""),sep=""),
       ylab=paste("factor2 (",paste(f2,"%)",sep=""),sep=""),main="Biplot for PLS using mixed models",
       xlim=c(-1,1),ylim=c(-1,1),col.main="blue",col.lab="blue")
if(length(nop)!=0){
text(anaplsP[-nop,1],anaplsP[-nop,2],rownames(anaplsP)[-nop],col="green",cex=0.6)
}else{text(anaplsP[,1],anaplsP[,2],rownames(anaplsP),col="green",cex=0.6)}
abline(h=0)
abline(v=0)}


dev.off()


  
 } 
} 
