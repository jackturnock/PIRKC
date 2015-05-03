BasicModelFits<-function(REP,PAR,DAT,BIO,BootCIs=NA)
{
#==read in REP and DAT files
temp<-grep("start/end",DAT)
starYear<-as.numeric(unlist(strsplit(DAT[temp+1],split="\t")))[1]
temp<-grep("survey",DAT)[1]
yearsDat<-as.numeric(unlist(strsplit(DAT[temp+1],split="\t")))[1]
temp<-grep("survey",DAT)[2]
yearsData<-as.numeric(unlist(strsplit(DAT[temp+1],split="\t")))
BegOffset<-yearsData[1]-starYear
temp<-grep("LengthBins",DAT)
LenBins<-as.numeric(unlist(strsplit(DAT[temp+1],split="\t")))
LenBins<-LenBins[!is.na(LenBins)]
lenbinN<-length(LenBins)
Years<-seq(starYear+BegOffset,starYear+BegOffset+yearsDat-1)

temp<-grep("pred female surv",REP)
predSurF<-as.numeric(unlist(strsplit(REP[temp+1],split=" ")))
temp<-grep("pred male surv",REP)
predSurM<-as.numeric(unlist(strsplit(REP[temp+1],split=" ")))

temp<-grep("obs female surv",REP)
obsSurF<-as.numeric(unlist(strsplit(REP[temp+1],split=" ")))
temp<-grep("obs male surv",REP)
obsSurM<-as.numeric(unlist(strsplit(REP[temp+1],split=" ")))

temp<-grep("CV",DAT)[1]
survMcv<-as.numeric(unlist(strsplit(DAT[temp+1],split="\t")))
temp<-grep("CV",DAT)[2]
survFcv<-as.numeric(unlist(strsplit(DAT[temp+1],split="\t")))

temp<-grep("pred catch N",REP)
predCatN<-as.numeric(unlist(strsplit(REP[temp+1],split=" ")))
temp<-grep("obs catch N",REP)
obsCatN<-as.numeric(unlist(strsplit(REP[temp+1],split=" ")))

TrawlYears<-seq(1991,yearsData[(length(yearsData)-1)])
temp<-grep("pred trawl bycatch fem",REP)
predTrawlF<-as.numeric(unlist(strsplit(REP[temp+1],split=" ")))[2:(yearsDat+BegOffset+1)]
temp<-grep("pred trawl bycatch male",REP)
predTrawlM<-as.numeric(unlist(strsplit(REP[temp+1],split=" ")))[2:(yearsDat+BegOffset+1)]
predTrawl<-predTrawlF+predTrawlM

temp<-grep("obs trawl bycatch",REP)
obsTrawl<-as.numeric(unlist(strsplit(REP[temp+1],split=" ")))[2:(length(TrawlYears)+1)]

temp<-grep("sample",DAT)[1]
MsurvNsamp<-as.numeric(unlist(strsplit(DAT[temp+1],split="\t")))
temp<-grep("sample",DAT)[2]
FsurvNsamp<-as.numeric(unlist(strsplit(DAT[temp+1],split="\t")))

temp<-grep("obs surv M at len",REP)
obsSurMlen<-matrix(as.numeric(unlist(strsplit(REP[(temp+1):(temp+yearsDat)],split=" "))),ncol=lenbinN+1,byrow=T)
temp<-grep("obs surv F at len",REP)
obsSurFlen<-matrix(as.numeric(unlist(strsplit(REP[(temp+1):(temp+yearsDat)],split=" "))),ncol=lenbinN+1,byrow=T)

temp<-grep("pred surv F at len",REP)
predSurFlen<-matrix(as.numeric(unlist(strsplit(REP[(temp+1):(temp+yearsDat+BegOffset)],split=" "))),ncol=lenbinN+1,byrow=T)
temp<-grep("pred surv M at len",REP)
predSurMlen<-matrix(as.numeric(unlist(strsplit(REP[(temp+1):(temp+yearsDat+BegOffset)],split=" "))),ncol=lenbinN+1,byrow=T)

temp<-grep("B35",REP)
B35<-as.numeric(unlist(strsplit(REP[temp+1],split=" ")))
temp<-grep("male spawning",REP)
mspbio<-as.numeric(unlist(strsplit(REP[temp+1],split=" ")))
temp<-grep("female spawning",REP)
fspbio<-as.numeric(unlist(strsplit(REP[temp+1],split=" ")))

#==Par file quantities
cnt<-grep("mean_log_rec",PAR)
meanREC<-as.numeric(PAR[cnt+1])
cnt<-grep("rec_dev",PAR)
recDevs<-as.numeric(PAR[(cnt+1+BegOffset):(cnt+yearsDat+BegOffset)])

fmortYrs<-seq(1993,1998)
cnt<-grep("log_avg_fmort_dir",PAR)
logFdir<-as.numeric(PAR[cnt+1])
cnt<-grep("fmort_dir_dev",PAR)
fmort_dir_dev<-as.numeric(PAR[(cnt+1):(cnt+length(fmortYrs))])

cnt<-grep("log_avg_fmort_trawl",PAR)
logFtrawl<-as.numeric(PAR[cnt+1])
cnt<-grep("fmort_trawl_dev",PAR)
fmort_trawl_dev<-as.numeric(PAR[(cnt+1):(cnt+length(TrawlYears))])

#==estimated selectivities
cnt<-grep("srv_q",PAR)
srvQ<-as.numeric(PAR[cnt+1])
cnt<-grep("srv_sel50",PAR)
sel50surv<-as.numeric(PAR[cnt+1])
cnt<-grep("srv_sel95",PAR)
sel95surv<-as.numeric(PAR[cnt+1])

cnt<-grep("fish_sel50",PAR)
sel50fish<-as.numeric(PAR[cnt+1])
cnt<-grep("fish_sel95",PAR)
sel95fish<-as.numeric(PAR[cnt+1])

SurvSel<-srvQ/( 1 + exp( -1*log(19)*(LenBins-sel50surv)/(sel95surv-sel50surv)))
FishSel<-1/( 1 + exp( -1*log(19)*(LenBins-sel50fish)/(sel95fish-sel50fish)))

cnt<-grep("af",PAR)
af<-as.numeric(PAR[cnt+1])
cnt<-grep("am",PAR)[2]
am<-as.numeric(PAR[cnt+1])
cnt<-grep("bf",PAR)
bf<-as.numeric(PAR[cnt+1])
cnt<-grep("bm",PAR)
bm<-as.numeric(PAR[cnt+1])

#=========================
# CIs calculated from CVs from Bob
#=================================
#==variances of survey estimates (Bob assumes poisson)
varM<-(log(obsSurM[2:(yearsDat+1)])*survMcv)^2
varF<-(log(obsSurF[2:(yearsDat+1)])*survFcv)^2
sdM<-log((survMcv^2)+1)
sdF<-log((survFcv^2)+1)

mCIup<-exp(log(obsSurM[2:(yearsDat+1)])+2*sdM)
mCIdn<-exp(log(obsSurM[2:(yearsDat+1)])-2*sdM)

fCIup<-exp(log(obsSurF[2:(yearsDat+1)])+2*sdF)
fCIdn<-exp(log(obsSurF[2:(yearsDat+1)])-2*sdF)

dev.new()
par(xpd=F)
par(mfrow=c(4,1),mar=c(.1,5,1,1),oma=c(2,0,0,0))
plot(obsSurF[2:(yearsDat+1)]/10000~Years,ylim=c(0,max(obsSurF[2:(yearsDat+1)])/10000),las=1,
ylab="Survey F numbers
 (10000s)",xaxt='n')

#==three year Avg
#==use the CV to weight by inverse variance
useF<-obsSurF[2:(yearsDat+1)]
runAvgF<-rep(0,length(useF))
for(x in 2:(length(useF)-1))
 runAvgF[x]<-weighted.mean(useF[(x-1):(x+1)],w=1/varF[(x-1):(x+1)])

lines(runAvgF/10000~Years,lty=2,col=2)
lines(predSurF[(2+BegOffset):((yearsDat+1)+BegOffset)]/10000~Years)
for(j in 1:length(fCIdn))
 segments(x0=Years[j],x1=Years[j],y0=fCIup[j]/10000,y1=fCIdn[j]/10000)

#==plot observed and predicted males
plot(obsSurM[2:(yearsDat+1)]/10000~Years,ylim=c(0,max(obsSurM[2:(yearsDat+1)])/10000),xaxt='n',
las=1,ylab="Survey M numbers
 (10000s)")
lines(predSurM[(2+BegOffset):((yearsDat+1)+BegOffset)]/10000~Years)
for(j in 1:length(mCIup))
 segments(x0=Years[j],x1=Years[j],y0=mCIup[j]/10000,y1=mCIdn[j]/10000)
cbind(predSurM[(2+BegOffset):((yearsDat+1)+BegOffset)]/10000,Years)

#==3 year running average inverse variance weighted
useM<-obsSurM[2:(yearsDat+1)]
runAvgM<-rep(0,length(useM))
for(x in 2:(length(useM)-1))
 runAvgM[x]<-weighted.mean(useM[(x-1):(x+1)],w=1/varM[(x-1):(x+1)])

lines(runAvgM/10000~Years,lty=2,col=2)

plot(obsCatN[2:(yearsDat+1)]/1000~Years,ylab="Catch (1000s)",las=1,xaxt='n')
lines(predCatN[2:(yearsDat+1)]/1000~Years)
temp<-grep("obs surv M at len",REP)
obsSurMlen<-matrix(as.numeric(unlist(strsplit(REP[(temp+1):(temp+yearsDat)],split=" "))),nrow=yearsDat,byrow=T)
apply(obsSurMlen,1,sum,na.rm=T)

plot(obsTrawl~TrawlYears,ylab="Trawl (t)",las=1,xlim=c(min(Years),max(Years)))
lines(predTrawl[(1+BegOffset):(yearsDat+BegOffset)]~Years)
cbind(predSurM[(2+BegOffset):((yearsDat+1)+BegOffset)]/10000,Years)

#==========================================
#==Find the tier 4 OFL using the running average
#==========================================
fspTier4<-BIO$MT_MATURE_FEMALE
fspTier4cv<-BIO$CV_MT_MATURE_FEMALE
mspTier4<-BIO$MT_MATURE_MALE
mspTier4cv<-BIO$CV_MT_MATURE_MALE

varM<-(log(mspTier4+0.0001)*mspTier4cv)^2
varF<-(log(fspTier4+0.0001)*fspTier4cv)^2
sdM<-log((mspTier4cv)+1)
sdF<-log((fspTier4cv)+1)

useM<-mspTier4
runAvgMb<-rep(0,length(useM))
for(x in 2:(length(useM)-1))
 runAvgMb[x]<-weighted.mean(useM[(x-1):(x+1)],w=1/varM[(x-1):(x+1)])

x<-length(useM)
 runAvgMb[x]<-weighted.mean(useM[(x-1):(x)],w=1/varM[(x-1):(x)])

Tier4BMSY<-mean(runAvgMb[17:length(runAvgM)])
Tier4FOFL<-0.18
CurB<-runAvgMb[length(runAvgMb)]*exp(-4*.18/12)
testB<-CurB/Tier4BMSY
use4fOFL<-0.18
beta<-0.25
alpha<-0.1

if(testB<1 & testB>=beta)
 use4fOFL<-Tier4FOFL*((testB-alpha)/(1-alpha))
if(testB<beta)
 use4fOFL<-0

OFL4<-(1-exp(-use4fOFL))*CurB*exp(-4*.18/12)

BootMMB<-rnorm(100000,log(CurB),sqrt(sdM[yearsDat]^2+0.3^2))
OFL4dist<-(1-exp(-use4fOFL))*exp(BootMMB)*exp(-4*.18/12)
temp<-sort(OFL4dist)
temp[0.49*length(temp)]
temp[0.05*length(temp)]
temp[0.95*length(temp)]

dev.new()
par(mar=c(5,5,1,1))
hist(OFL4dist,main="",las=1,xlab="Tier 4 OFL distribution (t)",ylab="",
col='grey',xlim=c(0,12500),cex=.8)
mtext(side=2,"Frequency",line=3.4)

#================================
# bootstrapped CIs (if available)
#================================
if(!is.na(BootCIs))
{
dev.new()
par(xpd=F)
par(mfrow=c(4,1),mar=c(.1,5,1,1),oma=c(2,0,0,0))
plot(obsSurF[2:(yearsDat+1)]/10000~Years,ylim=c(0,max(obsSurF[2:(yearsDat+1)])/10000),las=1,
ylab="Survey F numbers
 (10000s)",xaxt='n')

#==three year Avg
#==use the CV to weight by inverse variance
useF<-obsSurF[2:(yearsDat+1)]
runAvgF<-rep(0,length(useF))
for(x in 2:(length(useF)-1))
 runAvgF[x]<-weighted.mean(useF[(x-1):(x+1)],w=1/varF[(x-1):(x+1)])
 x<-length(useF)
 runAvgF[x]<-weighted.mean(useF[(x-1):(x)],w=1/varF[(x-1):(x)])
lines(runAvgF/10000~Years,lty=2,col=2)
lines(predSurF[(2+BegOffset):((yearsDat+1)+BegOffset)]/10000~Years)
for(j in 1:nrow(CIsurvey))
 segments(x0=Years[j],x1=Years[j],y0=CIsurvey[j,2]/10000,y1=CIsurvey[j,3]/10000)

#==plot observed and predicted males
plot(obsSurM[2:(yearsDat+1)]/10000~Years,ylim=c(0,max(obsSurM[2:(yearsDat+1)])/10000),xaxt='n',
#las=1,ylab="Survey M numbers (10000s)",xaxt='n')
las=1,ylab="Survey M numbers
 (10000s)")
lines(predSurM[(2+BegOffset):((yearsDat+1)+BegOffset)]/10000~Years)
for(j in 1:length(mCIup))
 segments(x0=Years[j],x1=Years[j],y0=CIsurvey[j,4]/10000,y1=CIsurvey[j,5]/10000)

#==3 year running average inverse variance weighted
useM<-obsSurM[2:(yearsDat+1)]
runAvgM<-rep(0,length(useM))
for(x in 2:(length(useM)-1))
 runAvgM[x]<-weighted.mean(useM[(x-1):(x+1)],w=1/varM[(x-1):(x+1)])
 x<-length(useF)
 runAvgM[x]<-weighted.mean(useM[(x-1):(x)],w=1/varM[(x-1):(x)])
lines(runAvgM/10000~Years,lty=2,col=2)

temp<-grep("pred catch N",REP)
predCatN<-as.numeric(unlist(strsplit(REP[temp+1],split=" ")))
temp<-grep("obs catch N",REP)
obsCatN<-as.numeric(unlist(strsplit(REP[temp+1],split=" ")))

plot(obsCatN[2:(yearsDat+1)]/1000~Years,ylab="Catch (1000s)",las=1,xaxt='n')
lines(predCatN[2:(yearsDat+1)]/1000~Years)
temp<-grep("obs surv M at len",REP)
obsSurMlen<-matrix(as.numeric(unlist(strsplit(REP[(temp+1):(temp+yearsDat)],split=" "))),nrow=yearsDat,byrow=T)
apply(obsSurMlen,1,sum,na.rm=T)

plot(obsTrawl~TrawlYears,ylab="Trawl (t)",las=1,xlim=c(min(Years),max(Years)))
lines(predTrawl[(1+BegOffset):(yearsDat+BegOffset)]~Years)
}

##====LEN FREQ PLOTS==========

dev.new()
mat<-matrix(seq(1,yearsDat-13),ncol=3)
layout(mat)
par(mar=c(0,0,0,0),oma=c(0,0,2,0))
for(x in 14:yearsDat)
 {
  plot(obsSurMlen[x,],axes=F,ylim=c(0,.2),type="l")
  lines(predSurMlen[x+BegOffset,],lty=2,col=2)
  legend("topleft",legend=Years[x],bty='n')
  legend("topright",legend=round(MsurvNsamp[x]),bty='n')  
 }
mtext(side=3,outer=T,"Males")
dev.new()
mat<-matrix(seq(1,yearsDat-13),ncol=3)
layout(mat)
par(mar=c(0,0,0,0),oma=c(0,0,2,0))
for(x in 14:yearsDat)
 {
  plot(obsSurFlen[x,],axes=F,ylim=c(0,.2),type="l")
  lines(predSurFlen[x+BegOffset,],lty=2,col=2)
  legend("topleft",legend=Years[x],bty='n')
  legend("topright",legend=round(MsurvNsamp[x]),bty='n')  
 }
mtext(side=3,outer=T,"Females")
#==================================================
#   Par file quantities
#================================================
dev.new()
par(mfrow=c(4,1),mar=c(1,4,1,1),oma=c(3,1,1,1))
TotRec<-exp(meanREC+recDevs)
TotFdir<-exp(logFdir+fmort_dir_dev)
TotFdir<-c(rep(0,18),TotFdir,rep(0,(yearsDat-24)))
TotFtrawl<-exp(logFtrawl+fmort_trawl_dev)

plot(TotRec[1:yearsDat]/10000~Years[1:yearsDat],type="l",las=1,ylab="Recruitment (10000s)")
plot(TotFdir~Years,type="l",las=1,ylab="Directed F",xlim=c(min(Years),max(Years)))
plot(TotFtrawl~TrawlYears,type="l",las=1,ylab="Trawl F",xlim=c(min(Years),max(Years)))

plot(SurvSel~LenBins,type="l",xlab="Length",ylab="Selectivity")
mtext(side=1,"Length",line=2)
lines(FishSel~LenBins,lty=2,col=2)
legend("topleft",bty='n',col=c(1,2),lty=c(1,2),legend=c("Survey","Fishery"))


#==plot estimated growth
dev.new()
par(mar=c(4,4,1,1))
groInc<-(am+bm*LenBins)-LenBins
plot(groInc~LenBins,type="l",ylim=c(0,40),ylab="Growth increment (mm)",
xlab="Length")

groInc<-(af+bf*LenBins)-LenBins
lines(groInc~seq(37.5,207.5,5),type="l",ylim=c(0,40),ylab="Growth increment (mm)",
xlab="Length",lty=2,col=2)
legend("topleft",bty='n',lty=c(1,2),col=c(1,2),c("Males","Females"))


#============================
# male spawning biomass and reference points
#=========================
dev.new()
plot(mspbio[6:(length(Years)+5)]/1000000~Years,type="l",las=2,
ylab="Male spawning biomass (1000kg)",xaxt='n')
axis(1)
abline(h=B35/1000000,lty=2,col=2)
legend("topleft",paste(Years[length(Years)] ,"MMB = ",mspbio[(length(Years)+1)]/1000000," (t)"),bty='n')
legend("bottomright",lty=2,col=2,"BMSY",bty='n')


#============================
# size transition matrix, fraction recruiting, molting, selectivity
#=========================
temp<-grep("#size transition matrix",REP)[1]
SizeTransMat<-matrix(as.numeric(unlist(strsplit(REP[(temp+1):(temp+(lenbinN))],split=" "))),ncol=lenbinN+1,byrow=T)
temp<-grep("fraction recruiting",REP)
FracRec<-as.numeric(unlist(strsplit(REP[temp+1],split=" ")))
FracRec<-FracRec[!is.na(FracRec)]

dev.new()
par(mfrow=c(2,2),mar=c(.1,.1,4,.1),oma=c(5,5,0,5))
plot(0,ylim=c(0,nrow(SizeTransMat)+1),xaxt='n',xlim=c(min(LenBins),
max(LenBins)),yaxt='n',ylab="Length",xlab="Length",main="Size transition matrix")
counter<-0
for(x in 1:nrow(SizeTransMat))
{
 lines((SizeTransMat[x,2:36]+counter)~LenBins)
 counter<-counter+1
}
axis(side=1,at=LenBins,labels=LenBins,las=1)
axis(side=2,at=seq(0,nrow(SizeTransMat)-1,4),labels=LenBins[seq(1,length(LenBins),4)],las=1)
mtext(side=2,"Pre-molt length",line=3.3)

barplot(FracRec,names.arg=LenBins,main="Fraction recruiting",xlim=c(0,10),yaxt='n')
axis(side=4,las=2)

#=TAKE THESE FROM .CTL and .DAT FILE===============
MaleMolt<-1 - 1./(1.+exp(-1.*.101*(LenBins-140.19)))
TrawlSel<- 1./(1.+exp(-1.*.053*(LenBins-164.84)))
maleMat<-1/(1 + (5.842 * 10^14) * exp((LenBins+2.5) * -0.288))
femMat<- 1/(1 + (1.416 * 10^13) * exp((LenBins+2.5) * -0.297))

plot(MaleMolt~LenBins,type="l",ylab="Probability of molting",xlab="Length",las=1)
mtext(side=2,line=3,"Probability",outer=T,adj=.15)
legend("topleft",bty='n',c("Male molting probability","Female maturity","Male maturity"),
col=c(1,2,3),lty=c(1,2,3),cex=.8)
lines(femMat~LenBins,col=2,lty=2)
lines(maleMat~LenBins,col=3,lty=3)

plot(TrawlSel~LenBins,type="l",ylab="Trawl selectivity",xlab="Length",las=1,ylim=c(0,1),yaxt='n')
axis(side=4,las=2)
legend("bottomleft",bty='n',c("Trawl selectivity (bycatch)","Fishery selectivity"),
col=c(1,2),lty=c(1,2))
mtext(side=1,"Length",outer=T,line=2)
lines(FishSel~LenBins,lty=2,col=2)


#==========================
# make a table of estimated parameters
#==========================
cnt<-grep("#",PAR)
parNames<-(PAR[cnt+1])
parNames<-parNames[-1]
cnt<-grep("#",PAR)
parVals<-(PAR[cnt+2])
parVals<-parVals[-1]

estPars<-round(matrix(rbind(as.numeric(parVals)),ncol=18,nrow=1),2)
colnames(estPars)<-parNames
outs<-t(estPars)[-c(7,9,11),]
write.csv(outs,paste(curDir,"/pars.csv",sep=""))


#============================
# table of estimated abundance, abundance, biomass, 
# running average, recruitment, catch/biomass
#=========================
#==plot model fits
fspTier4<-BIO$MT_MATURE_FEMALE
fspTier4cv<-BIO$CV_MT_MATURE_FEMALE
mspTier4<-BIO$MT_MATURE_MALE
mspTier4cv<-BIO$CV_MT_MATURE_MALE

varM<-(log(mspTier4+0.001)*mspTier4cv)^2
varF<-(log(fspTier4+0.001)*fspTier4cv)^2
varF[11]<-50
sdM<-log((mspTier4cv)+1)
sdF<-log((fspTier4cv)+1)
mCIup<-exp(log(mspTier4+0.001)+2*sdM)
mCIdn<-exp(log(mspTier4+0.001)-2*sdM)
fCIup<-exp(log(fspTier4+0.001)+2*sdF)
fCIdn<-exp(log(fspTier4+0.001)-2*sdF)

#==3 year running average inverse variance weighted
useM<-mspTier4
runAvgMb<-rep(0,length(useM))
for(x in 2:(length(useM)-1))
 runAvgMb[x]<-weighted.mean(useM[(x-1):(x+1)],w=1/varM[(x-1):(x+1)])
x<-length(useM)
 runAvgMb[x]<-weighted.mean(useM[(x-1):(x)],w=1/varM[(x-1):(x)])

useF<-fspTier4
runAvgFb<-rep(0,length(useF))
for(x in 2:(length(useF)-1))
 runAvgFb[x]<-weighted.mean(useF[(x-1):(x+1)],w=1/varF[(x-1):(x+1)])
x<-length(useM)
 runAvgFb[x]<-weighted.mean(useF[(x-1):(x)],w=1/varF[(x-1):(x)])

#==output tables
PredsTable<-cbind(round(TotRec),round(fspbio[(2+BegOffset):(length(Years)+1+BegOffset)]/1000000),
round(mspbio[(2+BegOffset):(length(Years)+1+BegOffset)]/1000000),
round(predSurF[(2+BegOffset):((yearsDat+1)+BegOffset)]/1000,1),
round(predSurM[(2+BegOffset):((yearsDat+1)+BegOffset)]/1000,1))
colnames(PredsTable)<-c("Recruitment","MMB","Survey N (females)","Survey N (males)")
write.csv(PredsTable,paste(curDir,"table7.csv",sep=""))

PredsTable2<-cbind(Years,round(runAvgF/1000),round(runAvgM/1000),
round(runAvgFb),round(runAvgMb))
colnames(PredsTable2)<-c("Year","Females","Males","FMB","MMB")
write.csv(PredsTable2,paste(curDir,"table8.csv",sep=""))

#==========================================
#==Make summary table for the front end
#==========================================
# OFL,BMSY,CurrentMMB,B/BMSY,

firstMat<-matrix(ncol=5,nrow=3)
temp<-grep("OFL",REP)
firstMat[2,1]<-round(as.numeric(unlist(strsplit(REP[temp+1],split="\t")))[2]/1000000)
temp<-grep("B35",REP)
firstMat[2,2]<-round(as.numeric(unlist(strsplit(REP[temp+1],split="\t")))[1]/1000000)
temp<-grep("FOFL",REP)
FOFL<-round(as.numeric(unlist(strsplit(REP[temp+1],split="\t")))[1],2)
temp<-grep("Bmsy_tier4",REP)
firstMat[3,2]<-round(as.numeric(unlist(strsplit(REP[temp+1],split="\t")))[1]/1000000)
Bmsy_tier4<-firstMat[3,2]
temp<-grep("OFL_tier4",REP)
firstMat[3,1]<-round(as.numeric(unlist(strsplit(REP[temp+1],split="\t")))[2]/1000000)
currMMB<-mspbio[(length(Years)+1+BegOffset)]/1000000
currMMB4<-runAvgMb[length(runAvgMb)]*exp(-3*.18/12)

colnames(firstMat)<-c("OFL","BMSY","MMB","Ratio","gamma")
rownames(firstMat)<-c("RunAvg/Tier4","Assess/Tier3","Assess/Tier4")
firstMat[3,4]<-round(currMMB/(Bmsy_tier4),2)
firstMat[2,4]<-round(currMMB/(B35/1000000),2)
firstMat[1,4]<-round(currMMB4/Tier4BMSY,2)

firstMat[2,3]<-round(currMMB)
firstMat[1,3]<-round(currMMB4)
firstMat[3,3]<-round(currMMB)

firstMat[2,5]<-1
firstMat[1,5]<-1
firstMat[3,5]<-1

firstMat[1,2]<-round(Tier4BMSY)
firstMat[1,1]<-round(OFL4)
write.csv(firstMat,paste(curDir,"IntroTable.csv",sep=""))
}

