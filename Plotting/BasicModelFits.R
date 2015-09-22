BasicModelFits<-function(REP,PAR,DAT,BIO,BootCIs=NA,mcEval,STDfile,PSV,CAT)
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

temp<-grep("sel_fish_discard",REP)
SelFishDiscard<-as.numeric(unlist(strsplit(REP[temp+1],split=" ")))

#==read in Catch data for Tier 4 running average
temp<-grep("trawl bycatch",CAT)[2]
TrawlBycatch<-as.numeric(unlist(strsplit(CAT[temp+1],split="\t")))
TrawlBycatch<-TrawlBycatch[!is.na(TrawlBycatch)]

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
TotFtrawl<-exp(logFtrawl+fmort_trawl_dev)

#==estimated selectivities
cnt<-grep("srv_q",PAR)
srvQ<-as.numeric(PAR[cnt+1:(BegOffset+yearsDat)])
cnt<-grep("srv_sel50",PAR)
sel50surv<-as.numeric(PAR[cnt+1])
cnt<-grep("srv_sel95",PAR)
sel95surv<-as.numeric(PAR[cnt+1])

cnt<-grep("fish_sel50",PAR)
sel50fish<-as.numeric(PAR[cnt+1])
cnt<-grep("fish_sel95",PAR)
sel95fish<-as.numeric(PAR[cnt+1])

SurvSel<-1/( 1 + exp( -1*log(19)*(LenBins-sel50surv)/(sel95surv-sel50surv)))
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
#==variances of survey estimates 
sdM<-sqrt(log((survMcv^2)+1))
sdF<-sqrt(log((survFcv^2)+1))

varM<-sdM^2
varF<-sdF^2
varF[varF=='NaN']<-max(varF,na.rm=T)
varM[varM==0]<-max(varM,na.rm=T)
varF[varM==0]<-max(varF,na.rm=T)

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
x<-length(useF)
 runAvgF[x]<-weighted.mean(useF[(x-1):(x)],w=1/useF[(x-1):(x)])

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
x<-length(useM)
 runAvgM[x]<-weighted.mean(useM[(x-1):(x)],w=1/varM[(x-1):(x)])

lines(runAvgM/10000~Years,lty=2,col=2)

plot(obsCatN[2:(yearsDat+1)]/1000~Years,ylab="Catch (1000s)",las=1,xaxt='n')
lines(predCatN[2:(yearsDat+1)]/1000~Years)
temp<-grep("obs surv M at len",REP)
obsSurMlen<-matrix(as.numeric(unlist(strsplit(REP[(temp+1):(temp+yearsDat)],split=" "))),nrow=yearsDat,byrow=T)
apply(obsSurMlen,1,sum,na.rm=T)

plot(obsTrawl~TrawlYears,ylab="Trawl (t)",las=1,xlim=c(min(Years),max(Years)))
lines(predTrawl[(1+BegOffset):(yearsDat+BegOffset)]~Years)
#cbind(predSurM[(2+BegOffset):((yearsDat+1)+BegOffset)]/10000,Years)

#==========================================
#==Find the tier 4 OFL using the running average
#==========================================
fspTier4<-BIO$BIOMASS_FEMALE_MATURE
fspTier4cv<-BIO$CV_BIOMASS_FEMALE_MATURE
mspTier4<-BIO$BIOMASS_MALE_GE120
mspTier4cv<-BIO$CV_BIOMASS_MALE_GE120

#==variances of survey estimates 
sdM<-sqrt(log((mspTier4cv^2)+1))
sdF<-sqrt(log((fspTier4cv^2)+1))

varM<-sdM^2
varF<-sdF^2
varF[varF=='NaN']<-max(varF,na.rm=T)
varM[varM==0]<-max(varM,na.rm=T)
varF[varM==0]<-max(varF,na.rm=T)

mCIup<-exp(log(mspTier4+0.001)+2*sdM)
mCIdn<-exp(log(mspTier4+0.001)-2*sdM)
fCIup<-exp(log(fspTier4+0.001)+2*sdF)
fCIdn<-exp(log(fspTier4+0.001)-2*sdF)

useM<-mspTier4
runAvgMb<-rep(0,length(useM))
for(x in 2:(length(useM)-1))
 runAvgMb[x]<-weighted.mean(useM[(x-1):(x+1)],w=1/varM[(x-1):(x+1)])

x<-length(useM)
 runAvgMb[x]<-weighted.mean(useM[(x-1):(x)],w=1/varM[(x-1):(x)])

Tier4BMSY<-mean(useM[17:length(useM)]*exp(-.18*(8/12)))
Tier4FOFL<-0.18
CurB<-runAvgMb[length(runAvgMb)]
ProjTrawlF<-mean(TotFtrawl[(length(TotFtrawl)-2):length(TotFtrawl)])
CurB<-CurB*exp(-4*.18/12)
CurB<-CurB*exp(-ProjTrawlF)
testB<-CurB/Tier4BMSY
use4fOFL<-0.18
beta<-0.25
alpha<-0.1

if(testB<1 & testB>=beta)
 use4fOFL<-Tier4FOFL*((testB-alpha)/(1-alpha))
if(testB<beta)
 use4fOFL<-0

OFL4tier4<-(1-exp(-use4fOFL))*CurB

BootMMB<-rnorm(1000000,log(CurB),sqrt(sdM[yearsDat]^2+0.3^2))
OFL4dist<-(1-exp(-use4fOFL))*exp(BootMMB)*exp(-4*.18/12)
temp<-sort(OFL4dist)
OFL405<-temp[.95*length(temp)]
OFL495<-temp[.05*length(temp)]
ABCtier4avg<-temp[.49*length(temp)]

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
 x<-length(useM)
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
plotSquares<-seq(1,3*ceiling((yearsDat-13)/3))
mat<-matrix(plotSquares,ncol=3)
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
plotSquares<-seq(1,3*ceiling((yearsDat-13)/3))
mat<-matrix(plotSquares,ncol=3)
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
par(mfrow=c(5,1),mar=c(.1,5,.1,1),oma=c(3,1,1,1))
TotRec<-exp(meanREC+recDevs)
TotFdir<-exp(logFdir+fmort_dir_dev)
TotFdir<-c(rep(0,18),TotFdir,rep(0,(yearsDat-24)))
TotFtrawl<-exp(logFtrawl+fmort_trawl_dev)

plot(TotRec[1:yearsDat]/10000~Years[1:yearsDat],type="l",las=1,ylab="Recruitment
 (10000s)",xaxt='n')
plot(TotFdir~Years,type="l",las=1,ylab="Directed F",xlim=c(min(Years),max(Years)),xaxt='n')
plot(TotFtrawl~TrawlYears,type="l",las=1,ylab="Trawl F",xlim=c(min(Years),max(Years)),xaxt='n')
plot(srvQ[6:(length(yearsData)+5)]~yearsData,type="l",las=1,ylab="Survey 
 catchability",xlim=c(min(Years),max(Years)),ylim=c(0,1))
par(mar=c(.1,5,2,1))
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
legend("topleft",bty='n',c("Trawl selectivity (bycatch)","Fishery selectivity","Discard selectivity"),
col=c(1,2,4),lty=c(1,2,4))
mtext(side=1,"Length",outer=T,line=2)
lines(FishSel~LenBins,lty=2,col=2)
lines(SelFishDiscard[2:length(SelFishDiscard)]~LenBins,lty=4,col=4)

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

#=====================================================
# TRANSPLANTED MCMC OUTPUT CODE HERE SO THAT ABC CAN BE OUTPUT IN SAME TABLES
#=====================================================

#====================
# pull out quantities from MCMC
#====================
temp<-grep("start/end",DAT)
starYear<-as.numeric(unlist(strsplit(DAT[temp+1],split="\t")))[1]
temp<-grep("start/end",DAT)
endYear<-as.numeric(unlist(strsplit(DAT[temp+1],split="\t")))[2]
temp<-grep("survey",DAT)[1]
yearsDat<-as.numeric(unlist(strsplit(DAT[temp+1],split="\t")))[1]
temp<-grep("survey",DAT)[2]
yearsData<-as.numeric(unlist(strsplit(DAT[temp+1],split="\t")))
BegOffset<-yearsData[1]-starYear
lenbinN<-35
LenBins<-seq(37.5,207.5,5)
Years<-seq(starYear+BegOffset,starYear+BegOffset+yearsDat-1)

TrawlYears<-seq(1991,yearsData[(length(yearsData)-1)])
NtrawlYr<-length(TrawlYears)

allYears<-seq(starYear,endYear)
lenBins<-seq(37.5,207.5,5)

#==must pull from the .STD file for this...CHANGE WHEN YOU CAN
tmp<-grep("Spbio",STDfile$name)
cnt<-length(allYears)
spbioDat<-cbind(rep("Spbio",cnt*2),
		    c(STDfile[tmp[1]:(tmp[1]+cnt-1),3],
                STDfile[(tmp[1]+length(tmp)/2):(tmp[1]+length(tmp)/2+cnt-1),3]),
		    c(STDfile[tmp[1]:(tmp[1]+cnt-1),4],
                STDfile[(tmp[1]+length(tmp)/2):(tmp[1]+length(tmp)/2+cnt-1),4]),
		    c(rep(1,cnt),rep(2,cnt)))
PSV <- file(paste(curDir,"/PIRKCv2.psv",sep=""), "rb")
nopar1 <- readBin(PSV, what=integer(), n=1)
mcmc1 <- readBin(PSV, what=numeric(), n=nopar1 * 10000)
mcmc1 <- matrix(mcmc1, byrow=TRUE, ncol=nopar1)

objfun<-mcEval[,ncol(mcEval)]
burnin<-0.1
thin<-1
library(coda)
dev.new()
par(mfrow=c(1,2))
plot(objfun)
hist(objfun)
thinOb<-objfun[(burnin*length(objfun)):length(objfun)]
thinOb1<-thinOb[seq(1,length(thinOb),thin)]
geweke.diag(thinOb1)


dev.new()
par(mfrow=c(1,2))
plot(thinOb1,type="b")
acf(thinOb1)

mcEval<-mcEval[(burnin*nrow(mcEval)):nrow(mcEval),]
mcEval<-mcEval[seq(1,nrow(mcEval),thin),]

Bmsy<-mcEval[,1]
F35<-mcEval[,2]
FOFL<-mcEval[,3]
OFL<-mcEval[,4]
median(OFL)
median(FOFL)

#==makes this so it can be reused each year because a 
#==new value is added for recruitment and trawl F
#==there is one for McEval and one for mcmc1, don't get them confused
NtrawlYr+yearsDat+15
untilTrawl<-13
untilRec<-14+NtrawlYr 
afterRec<-untilRec+yearsDat+5

betarec<-mcEval[,afterRec+9]
alpharec<-mcEval[,afterRec+8]
growthbetaM<-mcEval[,afterRec+6]
growthbetaF<-mcEval[,afterRec+7]
bm<-mcEval[,afterRec+5]
bf<-mcEval[,afterRec+4]
am<-mcEval[,afterRec+3]
af<-mcEval[,afterRec+2]

BMSY4<-mcEval[,afterRec+10]
FOFL4<-mcEval[,afterRec+11]
OFL4<-mcEval[,afterRec+12]


dev.new()
par(mfrow=c(1,2),mar=c(4,1,.2,.2),oma=c(0,3,0,0))
hist(BMSY4/1000000,main='',col='grey',las=1,xlab="Tier 4 Bmsy")
temp<-sort(BMSY4/1000000)
BMSY405<-temp[.05*length(temp)]
BMSY495<-temp[.95*length(temp)]

hist(OFL4/1000000,main='',col='grey',las=1,xlab="Tier 4 OFL")
temp<-sort(OFL4/1000000)
OFL405int<-temp[.05*length(temp)]
OFL495int<-temp[.95*length(temp)]
ABCtier4<-temp[.49*length(temp)]


dev.new()
par(mfrow=c(1,4),mar=c(2,2,.2,.2),oma=c(2,6,1,1))
hist(betarec,main='',col='grey')
hist(alpharec,main='',col='grey')
hist(growthbetaM,main='',col='grey')
hist(growthbetaF,main='',col='grey')

dev.new()
par(mfrow=c(1,4),mar=c(2,2,.2,.2),oma=c(5,2,1,1))
hist(bm,main='',col='grey')
mtext(side=2,"Frequency",line=2)
mtext(side=1,expression(b[male]),line=3)
hist(bf,main='',col='grey')
mtext(side=1,expression(b[female]),line=3)
hist(am,main='',col='grey')
mtext(side=1,expression(a[male]),line=3)
hist(af,main='',col='grey')
mtext(side=1,expression(a[female]),line=3)
mtext(outer=T,side=1,expression(Growth(Length)==a[sex]+(Length)*b[sex]),line=3)


RecInd<-yearsDat+15+NtrawlYr
mcmc1[1,]
MCMCavgRec<-mcmc1[,RecInd+1]
MCMCrecDevs<-mcmc1[,(RecInd+2):(RecInd+1+yearsDat+5)]
MCMCrecMat<-matrix(nrow=nrow(mcmc1),ncol=ncol(MCMCrecDevs))
for(i in 1:nrow(MCMCrecMat))
 MCMCrecMat[i,]<-exp(MCMCavgRec[i]+MCMCrecDevs[i,])

FdirInd<-yearsDat+5+2
MCMCavgF<-mcmc1[,(FdirInd+1)]
MCMCfDevs<-mcmc1[,(FdirInd+2):(FdirInd+7)]
MCMCfMat<-matrix(nrow=nrow(mcmc1),ncol=ncol(MCMCfDevs))
for(i in 1:nrow(MCMCrecMat))
 MCMCfMat[i,]<-exp(MCMCavgF[i]+MCMCfDevs[i,])


FtrawlInd<-yearsDat+5+2+7
MCMCavgFt<-mcmc1[,(FtrawlInd+1)]
MCMCfDevst<-mcmc1[,(FtrawlInd+2):(FtrawlInd+1+NtrawlYr)]
MCMCftMat<-matrix(nrow=nrow(mcmc1),ncol=ncol(MCMCfDevst))
for(i in 1:nrow(MCMCrecMat))
 MCMCftMat[i,]<-exp(MCMCavgFt[i]+MCMCfDevst[i,])

SelInd<-yearsDat+5
MCMCsel50<-mcmc1[,SelInd+1]
MCMCsel95<-mcmc1[,SelInd+2]
MCMCselMat<-matrix(nrow=nrow(mcmc1),ncol=length(lenBins))
for(i in 1:nrow(MCMCselMat))
 MCMCselMat[i,]<-1/(1+exp(-log(19)*(lenBins-MCMCsel50[i])/(MCMCsel95[i]-MCMCsel50[i])))

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

dev.new()

par(mfrow=c(2,2),mar=c(.1,1,1,1),oma=c(2,4,0,0))
multiplier<-1.5
fems<-spbioDat[,1]=="Spbio" & spbioDat[,4]=="1"
FspBio<-as.numeric(spbioDat[fems,2])/1000000
FspBioSD<-as.numeric(spbioDat[fems,3])/1000000

Fspbio95<-FspBio+2*FspBioSD
Fspbio05<-FspBio-2*FspBioSD
Fspbio75<-FspBio+.67*FspBioSD
Fspbio25<-FspBio-.67*FspBioSD

plot(FspBio~allYears,ylim=c(0,multiplier*max(fspTier4,na.rm=T)),las=1,xaxt='n',type="l")
polygon(x=c(allYears,rev(allYears)),y=c(Fspbio95,rev(Fspbio05)),col="grey85",border=F)
polygon(x=c(allYears,rev(allYears)),y=c(Fspbio75,rev(Fspbio25)),col="grey45",border=F)
legend("topleft",bty='n',"Females")

plot(fspTier4~Years,yaxt='n',xaxt='n',pch=16,ylim=c(0,multiplier*max(fspTier4,na.rm=T)))
for(j in 1:length(Years))
 segments(x0=Years[j],x1=Years[j],y0=fCIup[j],y1=fCIdn[j])
lines(runAvgFb~Years,yaxt='n',xaxt='n',lty=2,col=2)
legend("topleft",bty='n',pch=c(16,NA),lty=c(NA,2),col=c(1,2),c("Observed","3 year average"))

male<-spbioDat[,1]=="Spbio"& spbioDat[,4]=="2"
MspBio<-as.numeric(spbioDat[male,2])/1000000
MspBioSD<-as.numeric(spbioDat[male,3])/1000000

Mspbio95<-MspBio+2*MspBioSD
Mspbio05<-MspBio-2*MspBioSD
Mspbio75<-MspBio+.67*MspBioSD
Mspbio25<-MspBio-.67*MspBioSD

plot(MspBio~allYears,ylim=c(0,multiplier*max(mspTier4,na.rm=T)),las=1,type="l")
polygon(x=c(allYears,rev(allYears)),y=c(Mspbio95,rev(Mspbio05)),col="grey85",border=F)
polygon(x=c(allYears,rev(allYears)),y=c(Mspbio75,rev(Mspbio25)),col="grey45",border=F)
mtext(side=2,outer=T,"Mature biomass (t)",line=2.25)
legend("topleft",bty='n',"Males")
legend("topright",col=3,lty=3,bty='n',expression(B[paste(35,"%")]))
text(x=1980,y=1800000000,expression(B[paste(35,"%")]))

abline(h=max(Bmsy)/1000000,lty=3,col=3)
abline(h=min(Bmsy)/1000000,lty=3,col=3)
cbind(MspBio,allYears)

plot(mspTier4~Years,yaxt='n',pch=16,ylim=c(0,multiplier*max(mspTier4,na.rm=T)))
for(j in 1:length(Years))
 segments(x0=Years[j],x1=Years[j],y0=mCIup[j],y1=mCIdn[j])
lines(runAvgMb~Years,yaxt='n',xaxt='n',lty=2,col=2)

cbind(runAvgMb,Years)
cbind(runAvgFb,Years)
cbind(MspBio,allYears)
cbind(FspBio,allYears)
#==========================================
# histograms of management quantities
temp<-sort(Bmsy/1000000)
BMSY05<-temp[.05*length(temp)]
BMSY95<-temp[.95*length(temp)]

dev.new()
par(mfrow=c(3,1),mar=c(1,1,1,1),oma=c(2,4,0,0))
temp<-sort((MspBio)[45]/(Bmsy/1000000))
BRAT05<-temp[.05*length(temp)]
BRAT95<-temp[.95*length(temp)]
Brat<-hist(temp,plot=F)
temp<-Brat$counts/sum(Brat$counts)

barplot(temp,names.arg=Brat$mids)
text(x=15,y=.15,expression(B[current]),cex=1.5)
text(x=15.,y=.12,expression(B[paste(35,"%")]),cex=1.5)
lines(x=c(14,16),y=c(.135,.135))

temp<-sort(F35)
F3505<-temp[.05*length(temp)]
F3595<-temp[.95*length(temp)]
Brat<-hist(F35,plot=F)
temp<-Brat$counts/sum(Brat$counts)
barplot(temp,names.arg=Brat$mids)
text(x=8,y=.2,expression(F[paste(35,"%")]),cex=1.5)

temp<-sort(OFL/1000000)
OFL05<-temp[.05*length(temp)]
OFL95<-temp[.95*length(temp)]
ABCtier3<-temp[.49*length(temp)]

Brat<-hist(OFL/1000000,plot=F)
temp<-Brat$counts/sum(Brat$counts)
barplot(temp,names.arg=Brat$mids)
text(x=13.,y=.2,"OFL",cex=1.5)
median(OFL)
mtext(outer=T,side=2,"Posterior probability",line=2)


CatchLenyrs<-seq(1975,2013)
lenbinN<-35
LenBins<-seq(37.5,207.5,5)

dev.new()
par(mfrow=c(4,1),mar=c(1,4,1,1),oma=c(3,1,1,1))

sortCImat<-apply(MCMCrecMat,2,sort)/1000
TimeStep<-Years
MedVals<-apply(MCMCrecMat,2,median)/1000
sortCImat<-sortCImat[,6:ncol(sortCImat)]
MedVals<-MedVals[6:length(MedVals)]

quant95<-sortCImat[nrow(sortCImat)*.95,]
quant5<-sortCImat[nrow(sortCImat)*.05,]
quant25<-sortCImat[nrow(sortCImat)*.25,]
quant75<-sortCImat[nrow(sortCImat)*.75,]

plot(MedVals~TimeStep,ylim=c(min(sortCImat),max(quant95)),type="l",las=1,ylab="")
mtext(side=2,line=3.2,"Recruitment (10000s)",cex=.7)
polygon(x=c(TimeStep,rev(TimeStep)),y=c(quant95,rev(quant5)),col="grey85",border=F)
polygon(x=c(TimeStep,rev(TimeStep)),y=c(quant25,rev(quant75)),col="grey55",border="grey55")
lines(MedVals~TimeStep)

sortCImat<-apply(MCMCfMat,2,sort)
TimeStep<-Years
quant95<-rep(0,yearsDat)
quant95[19:24]<-sortCImat[nrow(sortCImat)*.95,]
quant5<-rep(0,yearsDat)
quant5[19:24]<-sortCImat[nrow(sortCImat)*.05,]
quant25<-rep(0,yearsDat)
quant25[19:24]<-sortCImat[nrow(sortCImat)*.25,]
quant75<-rep(0,yearsDat)
quant75[19:24]<-sortCImat[nrow(sortCImat)*.75,]

plot(0,ylim=c(0,max(quant95)),xlim=c(min(Years),max(Years)),type="l",las=1,ylab="")
mtext(side=2,line=3.2,"Fishing mortality",cex=.7)
polygon(x=c(TimeStep,rev(TimeStep)),y=c(quant95,rev(quant5)),col="grey85",border="grey85")
polygon(x=c(TimeStep,rev(TimeStep)),y=c(quant25,rev(quant75)),col="grey55",border="grey55")
lines(MedVals~TimeStep)
#plot(TotFdir~Years,type="l",las=1,ylab="Directed F",ylim=c(0,max(TotFdir)))

sortCImat<-apply(MCMCftMat,2,sort)
TimeStep<-Years
quant95<-rep(0,yearsDat)
quant95[(length(quant95)-NtrawlYr+1):yearsDat]<-sortCImat[nrow(sortCImat)*.95,]
quant5<-rep(0,yearsDat)
quant5[(length(quant95)-NtrawlYr+1):yearsDat]<-sortCImat[nrow(sortCImat)*.05,]
quant25<-rep(0,yearsDat)
quant25[(length(quant95)-NtrawlYr+1):yearsDat]<-sortCImat[nrow(sortCImat)*.25,]
quant75<-rep(0,yearsDat)
quant75[(length(quant95)-NtrawlYr+1):yearsDat]<-sortCImat[nrow(sortCImat)*.75,]

plot(0,ylim=c(0,max(MCMCftMat)),xlim=c(min(Years),max(Years)),type="l",las=1,ylab="")
mtext(side=2,line=3.2,"Bycatch mortality",cex=.7)
polygon(x=c(TimeStep,rev(TimeStep)),y=c(quant95,rev(quant5)),col="grey85",border="grey85")
polygon(x=c(TimeStep,rev(TimeStep)),y=c(quant25,rev(quant75)),col="grey55",border="grey55")

sortCImat<-apply(MCMCselMat,2,sort)
MedVals<-apply(MCMCselMat,2,median)
quant95<-sortCImat[nrow(sortCImat)*.95,]
quant5<-sortCImat[nrow(sortCImat)*.05,]
quant25<-sortCImat[nrow(sortCImat)*.25,]
quant75<-sortCImat[nrow(sortCImat)*.75,]

plot(MedVals~LenBins,ylim=c(min(sortCImat),max(quant95)),type="l",las=1,ylab="")
mtext(side=2,line=3.2,"Survey selectivity",cex=.7)
polygon(x=c(LenBins,rev(LenBins)),y=c(quant95,rev(quant5)),col="grey85",border=F)
polygon(x=c(LenBins,rev(LenBins)),y=c(quant25,rev(quant75)),col="grey55",border=F)

#============================
# size transition matrix and fraction recruiting
#=========================
temp<-grep("#size transition matrix",REP)
SizeTransMat<-matrix(as.numeric(unlist(strsplit(REP[(temp+1):(temp+(lenbinN))],split=" "))),ncol=lenbinN+1,byrow=T)
temp<-grep("fraction recruiting",REP)
FracRec<-as.numeric(unlist(strsplit(REP[temp+1],split=" ")))
FracRec<-FracRec[!is.na(FracRec)]
#snowFrac<-FracRec[!is.na(FracRec)]
dev.new()
barplot(FracRec,names.arg=LenBins,ylab="Fraction recruiting",xlim=c(0,10))
dev.new()
plot(0,ylim=c(0,nrow(SizeTransMat)+1),xlim=c(min(LenBins),max(LenBins)),yaxt='n',
xlab="Length",ylab="Length")
counter<-0
for(x in 1:nrow(SizeTransMat))
{
 lines((SizeTransMat[x,2:36]+counter)~LenBins)
 counter<-counter+1
}
axis(side=2,at=seq(0,nrow(SizeTransMat)-1,2),labels=LenBins[seq(1,length(LenBins),2)],las=1)



#=================================
#==output tables
#=================================

PredsTable<-cbind(round(TotRec),round(fspbio[(2+BegOffset):(length(Years)+1+BegOffset)]/1000000),
round(mspbio[(2+BegOffset):(length(Years)+1+BegOffset)]/1000000),
round(predSurF[(2+BegOffset):((yearsDat+1)+BegOffset)]/1000,1),
round(predSurM[(2+BegOffset):((yearsDat+1)+BegOffset)]/1000,1))
colnames(PredsTable)<-c("Recruitment","FMB","MMB","Survey N (females)","Survey N (males)")
write.csv(PredsTable,paste(curDir,"/table7.csv",sep=""))

PredsTable2<-cbind(Years,round(runAvgF/1000),round(runAvgM/1000),
round(runAvgFb),round(runAvgMb))
colnames(PredsTable2)<-c("Year","Females","Males","FMB","MMB")
write.csv(PredsTable2,paste(curDir,"/table8.csv",sep=""))

#==========================================
#==Make summary table for the front end
#==========================================
# OFL,BMSY,CurrentMMB,B/BMSY,

firstMat<-matrix(ncol=6,nrow=3)
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
currMMB4<-runAvgMb[length(runAvgMb)]

colnames(firstMat)<-c("OFL","BMSY","MMB","Ratio","gamma","ABC")
rownames(firstMat)<-c("RunAvg/Tier4","Assess/Tier3","Assess/Tier4")
firstMat[3,4]<-round(currMMB/(Bmsy_tier4),2)
firstMat[2,4]<-round(currMMB/(B35/1000000),2)
firstMat[1,4]<-round(currMMB4/Tier4BMSY,2)
firstMat[,6]<-c(ABCtier4avg,ABCtier3,ABCtier4)*0.75

firstMat[2,3]<-round(currMMB)
firstMat[1,3]<-round(currMMB4)
firstMat[3,3]<-round(currMMB)

firstMat[2,5]<-1
firstMat[1,5]<-1
firstMat[3,5]<-1

firstMat[1,2]<-round(Tier4BMSY)
firstMat[1,1]<-round(OFL4tier4)
write.csv(firstMat,paste(curDir,"/IntroTable.csv",sep=""))


#==write the variances
varMat<-matrix(nrow=3,ncol=8)
colnames(varMat)<-c("OFL","OFL","BMSY","BMSY","Ratio","Ratio","F35","F35")
rownames(varMat)<-c("Tier4/RunAvg","Tier4/Int","Tier3/Int")
varMat[1,]<-c(OFL405,OFL495,NA,NA,NA,NA,NA,NA)
varMat[2,]<-c(OFL405int,OFL495int,BMSY405,BMSY495,NA,NA,NA,NA)
varMat[3,]<-c(OFL05,OFL95,BMSY05,BMSY95,BRAT05,BRAT95,F3505,F3595)
write.csv(varMat,paste(curDir,"/VarTable.csv",sep=""))

}




