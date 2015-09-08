
MCMCoutputs<-function(REP,PAR,DAT,BIO,mcEval,STDfile,PSV)
{

#====================
# pull out quantities
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
burnin<-0.25
thin<-1

dev.new()
par(mfrow=c(1,2))
plot(objfun)
hist(objfun)
thinOb<-objfun[(burnin*length(objfun)):length(objfun)]
thinOb1<-thinOb[seq(1,length(thinOb),thin)]

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
temp[.05*length(temp)]
temp[.95*length(temp)]

hist(OFL4/1000000,main='',col='grey',las=1,xlab="Tier 4 OFL")
temp<-sort(OFL4/1000000)
temp[.05*length(temp)]
temp[.95*length(temp)]
temp[.49*length(temp)]


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

#==plot model fits

fspTier4<-BIO$BIOMASS_FEMALE_MATURE
fspTier4cv<-BIO$CV_BIOMASS_FEMALE_MATURE
mspTier4<-BIO$BIOMASS_MALE_POSTRECRUIT
mspTier4cv<-BIO$CV_BIOMASS_MALE_POSTRECRUIT

varM<-(log(mspTier4+0.001)*mspTier4cv)^2
varF<-(log(fspTier4+0.001)*fspTier4cv)^2
varF[varF=='NaN']<-max(varF,na.rm=T)
varM[varM==0]<-max(varM,na.rm=T)
varF[varM==0]<-max(varF,na.rm=T)

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
temp[.05*length(temp)]
temp[.95*length(temp)]


dev.new()
par(mfrow=c(3,1),mar=c(1,1,1,1),oma=c(2,4,0,0))
temp<-sort((MspBio)[45]/(Bmsy/1000000))
temp[.05*length(temp)]
temp[.95*length(temp)]
Brat<-hist(temp,plot=F)
temp<-Brat$counts/sum(Brat$counts)

barplot(temp,names.arg=Brat$mids)
text(x=15,y=.15,expression(B[current]),cex=1.5)
text(x=15.,y=.12,expression(B[paste(35,"%")]),cex=1.5)
lines(x=c(14,16),y=c(.135,.135))

temp<-sort(F35)
temp[.05*length(temp)]
temp[.95*length(temp)]
Brat<-hist(F35,plot=F)
temp<-Brat$counts/sum(Brat$counts)
barplot(temp,names.arg=Brat$mids)
text(x=8,y=.2,expression(F[paste(35,"%")]),cex=1.5)

temp<-sort(OFL/1000000)
temp[.05*length(temp)]
temp[.95*length(temp)]
temp[.49*length(temp)]
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
}
