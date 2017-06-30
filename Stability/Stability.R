
#Modify datos=datos1,selected=selected1,ExpDes=ExpDes1,BiplotFormat=BiplotFormat1,traits=traits1,LSMeans_par=LSMeans_par1,sta.parms=sta.parms1
STABILITY<-function(tmp=tmp1,datos=datos1,selected=selected1,ExpDes=ExpDes1,BiplotFormat=BiplotFormat1,traits=traits1,LSMeans_par=LSMeans_par1,sta.parms=sta.parms1){
library(reshape)
library(lsmeans)
library(ggplot2)

typedesign=ExpDes
  orgdir=getwd()

index <- which(as.character(datos$Loc)%in%selected)
if(length(index)>0)	datos <- datos[index,]

#Modify Loc = ENV,Rep=REP,Block= BLOCK and Loc = ENV,Entry = GEN
for (b in 1:length(traits)){  

  if (typedesign=="LSMeans"){
    datos$Entry=as.factor(datos$GEN)
    datos$Loc=as.factor(datos$ENV)
    raw=datos[order(datos$Loc,datos$Entry),c("Loc","Entry",traits[b])]
	varname=traits[b]
    colnames(raw)=c("ENV","GEN","YLD")
	raw$YLD=as.numeric(as.character(raw$YLD))
  repi=LSMeans_par[3]
  cme=LSMeans_par[1]
  }
  
  if (typedesign=="RCB"){
    datos$Entry=as.factor(datos$GEN)
    datos$Loc=as.factor(datos$ENV)
    datos$Rep=as.factor(datos$REP)
	repi<-max(as.numeric(datos$Rep))
	raw=datos[order(datos$Loc,datos$Entry,datos$Rep),c("Loc","Entry","Rep",traits[b])]
    varname=traits[b]
    colnames(raw)=c("ENV","GEN","REP","YLD")
	raw$YLD=as.numeric(as.character(raw$YLD))
    general=lm(YLD~ 1+ ENV + REP%in%ENV + GEN + GEN:ENV, data=raw)
	ver1=anova(general)	
	cme=ver1[5,3]	
  }
  
  if (typedesign=="Lattice"){
    datos$Entry=as.factor(datos$GEN)
    datos$Loc=as.factor(datos$ENV)
    datos$Rep=as.factor(datos$REP)
	datos$Block=as.factor(datos$BLOCK)
	repi=max(as.numeric(datos$Rep))
    raw=datos[order(datos$Loc,datos$Entry,datos$Rep,datos$Block),c("Loc","Entry","Rep","Block",traits[b])]
    varname=traits[b]
    colnames(raw)=c("ENV","GEN","REP","BLOCK","YLD")
	raw$YLD=as.numeric(as.character(raw$YLD))
	general=lm(YLD~ 1+ ENV + REP%in%ENV + BLOCK%in%(REP%in%ENV)+ GEN + GEN:ENV, data=raw)
	ver1=anova(general)
	cme=ver1[6,3]	
  }
  

datastability=cast(raw,GEN~ENV,mean,value="YLD",na.rm=T)
geno=as.character(datastability[,1])

datastability=datastability[,-1]

ind <- apply(datastability, 1, function(x) all(is.na(x)))
datastability <- datastability[ !ind, ]
geno<-geno[!ind]
ind <- apply(datastability, 2, function(x) all(is.na(x)))
datastability <- datastability[ ,!ind]


sdg=apply(datastability,1,sd,na.rm=T)
avg=apply(datastability,1,mean,na.rm=T)
a=dim(datastability)[2]
g=dim(datastability)[1]
yi.=as.vector(apply(datastability,1,mean,na.rm=T))
y.j=as.vector(apply(datastability,2,mean,na.rm=T))
y..=sum(datastability,na.rm=T)/(a*g)
Ij=y.j-y..

#CV Francis (CVFran)CVFran=(sdg/avg)*100

CVFran=(sdg/avg)*100
#############################################
##Eberhart y Russell (S2di)

zj=y.j-y..

sum7=0
for (i in 1:a){
  sum7=(as.vector(datastability[,i]))*zj[i]+sum7
}
biFinWil=sum7/sum(zj*zj)

sum4=0
for (i in 1:a){
  sum4=(datastability[,i]-(yi.))**2+sum4
}

S2di=(((sum4)-(biFinWil**2)*(sum(Ij*Ij)))/(a-2))-(cme/repi)

dataest=datastability
for ( i in 1:g){
dataest[i,]=yi.[i]+biFinWil[i]*Ij
}

##S2di= (((apply(dataest,1,sum)-(yi.**2/a))-(sum7**2/sum(zj*zj)))/(a-2))-(cme/repi)

R2=yi.
for (i in 1:g){
R2[i]=sum((dataest[i,]-sum(dataest[i,])/a)**2)/sum((datastability[i,]-sum(dataest[i,])/a)**2)
}

  newfold=paste("OutputStability&",paste(tmp1,paste(typedesign,paste("&",traits[b],sep=""),sep=""),sep=""),sep="")
  sal=paste(orgdir,paste("/",newfold,sep=""),sep="")
  dir.create(sal)
  setwd(sal)

Result=round(cbind(yi.,sdg,CVFran,biFinWil,S2di,R2),4)
R1=rbind(c("Mean","Sd","CV(%)","bi","S2di","R2"),Result)
colnames(R1)=c("*","*","Francis","Eberhart&Russell","*","*")
rownames(R1)=c("GEN",geno)
######################################
if(sta.parms[1]==TRUE){
#Shukla Stability (ri2Shukla & SiShukla)
sum1=0
for (i in 1:a){
sum1=(as.vector(datastability[,i]-yi.)-rep(y.j[i],dim(datastability)[1])+y..)**2+sum1
}

sum2=0
for (i in 1:a){
for (j in 1:g){
sum2=(as.vector(datastability[j,i]-yi.[j])-rep(y.j[i],dim(datastability)[1])+y..)**2+sum2
}
}
ri2Shukla=(g*(g-1)*sum1-sum2)/((a-1)*(g-1)*(g-2))

zj=y.j-y..

sum5=0
for (i in 1:a){
  sum5=(as.vector(datastability[,i])-rep(y.j[i],dim(datastability)[1]))*zj[i]+sum5
}
bi=sum5/sum(zj**2)

sum6=0
for (i in 1:a){
  sum6=(as.vector(datastability[,i]-yi.)-rep(y.j[i],dim(datastability)[1])+y..-(bi*zj[i]))**2+sum6
}

SiShukla=sum6

R1=cbind(R1,c("ri2",round(ri2Shukla,4)))
colnames(R1)[dim(R1)[2]]=c("Shuckla")
}
###########################################
if(sta.parms[2]==TRUE){
##Perkings and Jinks (Bi & DJi)
sum9=0
for (i in 1:a){
  sum9=(datastability[,i]-y.j[i]-yi.+y..)**2+sum9
}

Wi=sum9

sum3=0
for (i in 1:a){
  sum3=(datastability[,i]-y.j[i]-yi.+y..)*Ij[i]+sum3
}

Bi=sum3/sum(Ij*Ij)

DJi=(1/(a-2))*(Wi-(Bi**2)*(sum(Ij*Ij)))

R1=cbind(R1,cbind(c("Bi",round(Bi,4)),c("DJi",round(DJi,4))))
colnames(R1)[(dim(R1)[2]-1)]=c("Perkins&Jinks")
colnames(R1)[dim(R1)[2]]=c("*")

}
###########################################
if(sta.parms[3]==TRUE){
#Wrickes ecovalence (Wi)
sum9=0
for (i in 1:a){
  sum9=(datastability[,i]-y.j[i]-yi.+y..)**2+sum9
}

Wi=sum9

R1=cbind(R1,c("Wi",round(Wi,4)))
colnames(R1)[dim(R1)[2]]=c("Wricke's Ecovalence")
}
#############################################
if (sta.parms[4]==TRUE){
#Superiority measure Linn and Binns (Pi)
maxa=as.vector(apply(datastability,2,max))

sum8=0
for (i in 1:a){
  sum8=(datastability[,i]-maxa[i])**2+sum8
}
Pi=sum8/(2*a)

R1=cbind(R1,c("Pi",round(Pi,4)))
colnames(R1)[dim(R1)[2]]=c("Superiority Measure")

}
###########################################
if(sta.parms[5]==TRUE){
#Non parametric Hunh (rank2)
rank1=datastability
for (i in 1:a){
rank1[,i]=rank(datastability[,i])
}

ranki.=round(as.vector(apply(rank1,1,mean)),0)

a1=a-1
sum10=0
for (i in 1:a){
  sum10=(rank1[,i]-ranki.)**2+sum10
}

sum12=0
for (i in 1:a){
  sum12=abs(rank1[,i]-ranki.)+sum12
}
#Si2=round(sum10/sum12,2)
Si2=round(sum10/(a-1),2)


sum11=0
for (i in 1:a1){
  j=i+1
  sum11=abs(rank1[,i]-rank1[,j])+sum11
}
Si1=round((2*sum11)/(a*a1),2)

R1=cbind(R1,cbind(c("Si(1)",round(Si1,4)),c("Si2",round(Si2,4))))
colnames(R1)[(dim(R1)[2]-1)]=c("Non parametric Nassar&Huehn")
colnames(R1)[dim(R1)[2]]=c("*")
}
###########################################
write.csv(R1,paste("StabilityCoefficients",paste(traits[b],".csv")))

if (BiplotFormat=="pdf"){pdf(paste("PlotCV",paste(traits[b],".pdf")))}
if (BiplotFormat=="png"){png(paste("PlotCV",paste(traits[b],".png")))}
if (BiplotFormat=="wmf"){win.metafile(paste("PlotCV",paste(traits[b],".wmf")))}

#Plot CV%
rankingCV=data.frame(cbind(yi.,CVFran))
indexrankCV=which(rankingCV[,1]>y.. & rankingCV[,2]<mean(CVFran))
#good=rep("NS",length(CVFran))
#good[indexrankCV]="Good performance and stable"
#rankingCV$good=good
#rankingCV$geno=geno
mcv=mean(CVFran)
plot(yi.,CVFran,xlab="Mean",ylab="CV(%)",type="n")
if (length(indexrankCV)!=0){
text(yi.[-indexrankCV],CVFran[-indexrankCV],col="black",cex=0.8,labels=geno[-indexrankCV])
text(yi.[indexrankCV],CVFran[indexrankCV],col="red",cex=0.8,labels=geno[indexrankCV])
}
if (length(indexrankCV)==0){
text(yi.,CVFran,col="black",cex=0.8,labels=geno)
}
legend("bottomright", c("Good performance and stable"), col = c("red"),text.col = c("red"),pch = c(NA),bty="n",cex=0.6)
abline(v=y..)
abline(h=mcv)

#p=ggplot(rankingCV, aes(x=yi.,y=CVFran,label=geno,colour=good))+geom_text(size=3)+xlab("Mean")+ylab("CV(%)")+theme(legend.position = "bottom")+scale_colour_manual(values=c("red","blue"))+theme(legend.title = element_blank())
#p=p+geom_vline(aes(xintercept = y..))+geom_hline(aes(yintercept = mcv))
#print(p)
dev.off()

if (length(which(S2di%in%c(Inf,-Inf)==T))!=0){
print(paste("In",traits[b],"can not make the PlotEberhart&RussellCoeff"))
print("because not able to calculate the variance")
}
if (length(which(biFinWil%in%c(Inf,-Inf)==T))!=0){
print(paste("In",traits[b],"can not make the PlotEberhart&RussellCoeff"))
print("because not able to calculate the regression coefficients")
}
if (length(which(S2di%in%c(Inf,-Inf)==T))==0 && length(which(biFinWil%in%c(Inf,-Inf)==T))==0){

if (BiplotFormat=="pdf"){pdf(paste("PlotEberhart&RussellCoeff",paste(traits[b],".pdf")))}
if (BiplotFormat=="png"){png(paste("PlotEberhart&RussellCoeff",paste(traits[b],".png")))}
if (BiplotFormat=="wmf"){win.metafile(paste("PlotEberhart&RussellCoeff",paste(traits[b],".wmf")))}

#Plot Eberhart&Russell Coefficients
tcal=(biFinWil-1)/(sqrt(abs(S2di)/(sum(Ij*Ij))))
#tcal=(biFinWil-1)/(sqrt(abs(S2di)))
ttable=abs(qt((0.05/2),((repi*a)-1)))
adap=which(tcal>=ttable)
Fcal=(((sum4)-(biFinWil**2)*(sum(Ij*Ij)))/(a-2))/(cme/repi)
#Fcal=(a-1)*(sqrt(abs(S2di))/0.0001)**2
Ftable=qf((0.05/2),(a-1),(a*(g-1)))
#Ftable=qchisq((0.05/2),((repi*a)-1))
estab=which(Fcal<Ftable)

plot(S2di,biFinWil,xlab="Variability (S2di)",ylab="Coeff Reg (bi)",type="n")
legend("topright", c("no significant", "adaptable", "stable","adapt&stable"), col = c("black", "red", "blue","green"),
       text.col = c("black", "red", "blue","green"),pch = c(NA, NA, NA,NA),bty="n",cex=0.6)       
abline(h=1)
abline(v=0)
#print(length(adap))
#print(length(estab))
#print(tcal)
#print(ttable)
if (length(adap)==0 & length(estab)==0){
text(S2di,biFinWil,col="black",cex=0.8,labels=geno)
adapr=0
estabr=0
bothr=0
nosig=geno
}

if (length(adap)!=0 & length(estab)!=0){
both=intersect(adap,estab)
   if (length(both)!=0){
   text(S2di[both],biFinWil[both],col="green",cex=0.8,labels=geno[both])
   adap1=adap[-which(adap%in%both)]
   estab1=estab[-which(estab%in%both)]   
      if (length(adap1)==0 & length(estab1)==0){
	      if (length(both)!=length(geno)) {text(S2di[-both],biFinWil[-both],col="black",cex=0.8,labels=geno[-both])}}
      if (length(adap1)!=0 & length(estab1)!=0){
	      text(S2di[adap1],biFinWil[adap1],col="red",cex=0.8,labels=geno[adap1])
          text(S2di[estab1],biFinWil[estab1],col="blue",cex=0.8,labels=geno[estab1])
          amb=unique(c(adap1,estab1,both))
          if (length(amb)!=length(geno)) {text(S2di[-amb],biFinWil[-amb],col="black",cex=0.8,labels=geno[-amb])}}
      if (length(adap1)!=0 & length(estab1)==0){
          text(S2di[adap1],biFinWil[adap1],col="red",cex=0.8,labels=geno[adap1])
          amb=unique(c(adap1,both))
          if (length(amb)!=length(geno)) {text(S2di[-amb],biFinWil[-amb],col="black",cex=0.8,labels=geno[-amb])}}
	  if (length(adap1)==0 & length(estab1)!=0){
          text(S2di[estab1],biFinWil[estab1],col="blue",cex=0.8,labels=geno[estab1])
          amb=unique(c(estab1,both))	 
          if (length(amb)!=length(geno)) {text(S2di[-amb],biFinWil[-amb],col="black",cex=0.8,labels=geno[-amb])}} 
adapr=geno[adap]
estabr=geno[estab]
bothr=geno[both]
nosig=geno[-union(adap,estab)]		  
   }
   if (length(both)==0){
   text(S2di[adap],biFinWil[adap],col="red",cex=0.8,labels=geno[adap])
   text(S2di[estab],biFinWil[estab],col="blue",cex=0.8,labels=geno[estab])
   amb=unique(c(adap,estab))
   if (length(amb)!=length(geno)) {text(S2di[-amb],biFinWil[-amb],col="black",cex=0.8,labels=geno[-amb])}
adapr=geno[adap]
estabr=geno[estab]
bothr=0
nosig=geno[-union(adap,estab)]		  
   }
   


}

if (length(adap)!=0 & length(estab)==0){
text(S2di[adap],biFinWil[adap],col="red",cex=0.8,labels=geno[adap])
if (length(geno)!=length(adap)) text(S2di[-adap],biFinWil[-adap],col="black",cex=0.8,labels=geno[-adap])

adapr=geno[adap]
estabr=0
bothr=0
if (length(geno)!=length(adap)) nosig=geno[-adap]
if (length(geno)==length(adap)) nosig=0

}

if (length(adap)==0 & length(estab)!=0){
text(S2di[estab],biFinWil[estab],col="blue",cex=0.8,labels=geno[estab])
if (length(geno)!=length(estab)) text(S2di[-estab],biFinWil[-estab],col="black",cex=0.8,labels=geno[-estab])

adapr=0
estabr=geno[estab]
bothr=0
if (length(geno)!=length(estab)) nosig=geno[-estab]
if (length(geno)==length(estab)) nosig=0

}
dev.off()


cbind.fill<-function(...){
    nm <- list(...) 
    nm<-lapply(nm, as.matrix)
    n <- max(sapply(nm, nrow)) 
    do.call(cbind, lapply(nm, function (x) 
      rbind(x, matrix(, n-nrow(x), ncol(x))))) 
  }
  
  adapsta=cbind.fill(adapr,estabr,bothr,nosig)
  colnames(adapsta)=c("GenAdaptable","GenStable","GenAdap&Stable","GenNoSignificant")
write.csv(adapsta,paste("Adapt&StableGEN",paste(traits[b],".csv")))
}

}

}
