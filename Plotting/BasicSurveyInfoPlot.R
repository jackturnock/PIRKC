
#Figure 6 Total number of crab and stations with crab
##Figure 7 CIs from bootstrapped and given Bob's data
#Figure 8 and 9 density of males and females in mot recent year
#Figure 10-13 numbers at length and density

library(plotrix)
survDAT<-read.csv("C:/Users/Cody/Desktop/RKC/Red King Crab/PIRKC_survey.csv")
survDAT<-read.csv("C:/Users/Cody/Desktop/Latest RKC/PIRKC_DATA/CRABHAUL_PRIBRKC_TimeSeries.csv")
drvYear<-as.numeric(substr(survDAT$CRUISE,1,4))
SurvYR<-unique(drvYear)
AllStation<-unique(survDAT$GIS_STATION)

#==plot GIS stations
AllStnLoc<-matrix(ncol=2,nrow=length(AllStation))
for(w in 1:length(AllStation))
 {
  temp<-survDAT[survDAT$GIS_STATION==AllStation[w],]
  AllStnLoc[w,1]<-temp$MID_LATITUDE[1]
  AllStnLoc[w,2]<-temp$MID_LONGITUDE[1]
 }
dev.new(width=15,height=10)
par(mar=c(.1,.1,.1,.1))
plot(0,ylim=c(min(AllStnLoc[,1]),max(AllStnLoc[,1])),xlim=c(min(AllStnLoc[,2]),max(AllStnLoc[,2])))
for(w in 1:length(AllStation))
 text(AllStation[w],x=AllStnLoc[w,2],y=AllStnLoc[w,1],cex=.5)

#========================
#==find survey densities (total)
#========================
nmiSurv<-140350
DensityM<-matrix(nrow=length(SurvYR),ncol=length(AllStation))
DensityF<-matrix(nrow=length(SurvYR),ncol=length(AllStation))
StationYr<-matrix(nrow=length(SurvYR),ncol=length(AllStation))
for(y in 1:length(SurvYR))
{
 yrDAT<-survDAT[drvYear==SurvYR[y],]
 fileyr<-SurvYR[y]
 stationsUNQ<-(unique(yrDAT$GIS_STATION))
#==density at station
 for(j in 1:length(stationsUNQ))
  {
   stationALL<-yrDAT[yrDAT$GIS_STATION==stationsUNQ[j],]
   StationYr[y,j]<-as.character(stationsUNQ[j])
   Hauls<-(unique(stationALL$HAUL))
   Males<-stationALL[stationALL$SEX==1,]
   Females<-stationALL[stationALL$SEX==2,]
 #==densities across hauls in crabs per km^2
  tempDensM<-NULL
  tempDensF<-NULL
  for(k in 1:length(Hauls))
   {
    SampFactM<-Males$SAMPLING_FACTOR[which(Males$HAUL==Hauls[k])[1]] 
    AreaSweptM<-Males$AREA_SWEPT_VARIABLE[which(Males$HAUL==Hauls[k])[1]]
    tempDensM<-length(Males$HAUL==Hauls[k])*SampFactM/AreaSweptM

    SampFactF<-Females$SAMPLING_FACTOR[which(Females$HAUL==Hauls[k])[1]] 
    AreaSweptF<-Females$AREA_SWEPT_VARIABLE[which(Females$HAUL==Hauls[k])[1]]
    tempDensF<-length(Females$HAUL==Hauls[k])*SampFactF/AreaSweptF
    }
  DensityM[y,j]<-mean(tempDensM)
  DensityF[y,j]<-mean(tempDensF)
   }
}


#===================================
#  Densities only Pribolof district
#===================================
#  S of 58.39, W of 168
#===================================

DensityMpi<-matrix(nrow=length(SurvYR),ncol=length(AllStation))
DensityFpi<-matrix(nrow=length(SurvYR),ncol=length(AllStation))
StationYrPi<-matrix(nrow=length(SurvYR),ncol=length(AllStation))
TotalObsCrabF<-matrix(nrow=length(SurvYR),ncol=length(AllStation))
TotalObsCrabM<-matrix(nrow=length(SurvYR),ncol=length(AllStation))
HaulsRec<-matrix(nrow=2,ncol=length(SurvYR))
CaughtCrabM<-rep(0,length(SurvYR))
CaughtCrabF<-rep(0,length(SurvYR))

BinSize<-5
LengthBin<-seq(0,210,BinSize)
NatLengthF<-matrix(ncol=(length(LengthBin)-1),nrow=length(SurvYR))
NatLengthM<-matrix(ncol=(length(LengthBin)-1),nrow=length(SurvYR))
LenFreqF<-matrix(ncol=length(LengthBin)-1,nrow=length(SurvYR))
LenFreqM<-matrix(ncol=length(LengthBin)-1,nrow=length(SurvYR))

#===========================================================
# calculate number of crab and stations with crab
#=============================================================

for(y in 1:length(SurvYR))
{
 yrDAT<-survDAT[drvYear==SurvYR[y],]
 yrDAT<-yrDAT[yrDAT$MID_LONGITUDE< -168 & yrDAT$MID_LATITUDE < 58.65,]
 fileyr<-SurvYR[y]
 stationsUNQ<-(unique(yrDAT$GIS_STATION))
 yrDATm<-yrDAT[!is.na(yrDAT$SEX==1)&(yrDAT$SEX==1)=="TRUE",]
 yrDATf<-yrDAT[!is.na(yrDAT$SEX==2)&(yrDAT$SEX==2)=="TRUE",]

#==length frequencies
 tempM<-weighted.hist(yrDATm$LENGTH,breaks=LengthBin,w = yrDATm$SAMPLING_FACTOR, plot=F)
 tempF<-weighted.hist(yrDATf$LENGTH,breaks=LengthBin,w = yrDATf$SAMPLING_FACTOR, plot=F)

 NatLengthM[y,]<-tempM$counts
 NatLengthF[y,]<-tempF$counts
 LenFreqF[y,]<-tempF$density
 LenFreqM[y,]<-tempM$density

 CaughtCrabM[y]<-sum(tempM$counts)
 CaughtCrabF[y]<-sum(tempF$counts)

 HaulsRec[1,y]<-length(unique(yrDATm$GIS_STATION))
 HaulsRec[2,y]<-length(unique(yrDATf$GIS_STATION))

#==density at station
 for(j in 1:length(stationsUNQ))
  {
   stationALL<-yrDAT[yrDAT$GIS_STATION==stationsUNQ[j],]
   StationYrPi[y,j]<-as.character(stationsUNQ[j])
   Hauls<-(unique(stationALL$HAUL))
   #HaulsRec[y]<-HaulsRec[y]+sum(!is.na(Hauls))
   Males<-stationALL[stationALL$SEX==1,]
   Females<-stationALL[stationALL$SEX==2,]

 #==densities across hauls in crabs per km^2
  tempDensM<-NULL
  tempDensF<-NULL
  tempTotF<-NULL
  tempTotM<-NULL
  for(k in 1:length(Hauls))
   {
    SampFactM<-Males$SAMPLING_FACTOR[which(Males$HAUL==Hauls[k])[1]] 
    AreaSweptM<-Males$AREA_SWEPT_VARIABLE[which(Males$HAUL==Hauls[k])[1]]
    tempDensM<-append(tempDensM,length(Males$HAUL==Hauls[k])*SampFactM/AreaSweptM)
    tempTotM<-append(tempDensM,length(Males$HAUL==Hauls[k])*SampFactM)
    tempLenM<-Males[Males$HAUL==Hauls[k],]$LENGTH

    SampFactF<-Females$SAMPLING_FACTOR[which(Females$HAUL==Hauls[k])[1]] 
    AreaSweptF<-Females$AREA_SWEPT_VARIABLE[which(Females$HAUL==Hauls[k])[1]]
    tempDensF<-append(tempDensF,length(Females$HAUL==Hauls[k])*SampFactF/AreaSweptF)
    tempTotF<-append(tempDensF,length(Females$HAUL==Hauls[k])*SampFactF)
    tempLenF<-Females[Females$HAUL==Hauls[k],]$LENGTH

    }
  DensityMpi[y,j]<-mean(tempDensM)
  DensityFpi[y,j]<-mean(tempDensF)
  HaulsRec[y]<-sum(!is.na(DensityMpi[y,])|!is.na(DensityFpi[y,]))
  TotalObsCrabM[y,j]<-sum(tempTotM)
  TotalObsCrabF[y,j]<-sum(tempTotF)
   }
}

#==how many observations in each year?
#==how many hauls in each year?
DensMpi<-apply(DensityMpi,1,mean,na.rm=T)
DensFpi<-apply(DensityFpi,1,mean,na.rm=T)
TotF<-apply(TotalObsCrabF,1,sum,na.rm=T)
TotM<-apply(TotalObsCrabM,1,sum,na.rm=T)

par(mfrow=c(3,1),mar=c(.1,3,.1,.1))
plot(DensMpi~SurvYR,type="l")
lines(DensFpi~SurvYR,lty=2,col=2)
plot(TotM,type="l")
lines(TotF,lty=2,col=2)
plot(HaulsRec,type="l")

#==What is the area of the year in which the largest area was sampled?

par(mfrow=c(2,1),mar=c(.1,4,.1,.1),oma=c(3,0,0,0))
plot(CaughtCrabM~SurvYR,type="l",las=1,ylab="Total crabs observed",xaxt='n')
lines(CaughtCrabF~SurvYR,lty=2,col=2)
plot(HaulsRec[1,]~SurvYR,type="l",las=1,ylab="Stations with crab")
lines(HaulsRec[2,]~SurvYR,lty=2,col=2)




plot(TotM~TotF)
abline(a=0,b=1)
#===========================
#  plot densities by year
#============================

library("maps")
library("lattice")
library(PBSmapping)
library(mapdata)    #some additional hires data
library(maptools)   #useful tools such as reading shapefiles
library(mapproj)

PdfOut<-0
if(PdfOut==1)
 pdf(file="C:/Users/Cody/Desktop/PIRKC_survey_density.pdf")
if(PdfOut==2)
 pdf(file="C:/Users/Cody/Desktop/PIRKC_survey_density_zoom.pdf")

for(y in 1:length(SurvYR))
{
 library("maps")
 library("lattice")
 state.map <- map('worldHires', xlim=c(-175, -155.9), ylim=c(50, 65.5),
                  plot = FALSE, fill = TRUE,col='grey85')
#  S of 58.39, W of 168
 state.map <- map('worldHires', xlim=c(-173, -167.5), ylim=c(52, 58),
                  plot = FALSE, fill = TRUE,col='grey85')

state.info <- data.frame(density = DensityM[y,],
 long = AllStnLoc[match(StationYr[y,],AllStation,),2],
 lat = AllStnLoc[match(StationYr[y,],AllStation,),1])

state.info <- data.frame(density = DensityF[y,],
 long = AllStnLoc[match(StationYr[y,],AllStation,),2],
 lat = AllStnLoc[match(StationYr[y,],AllStation,),1])

state.info<-state.info[rowSums(is.na(state.info))!=3, ]
state.info[is.na(state.info),1]<-0

 panel.3dmap <- function(..., rot.mat, distance, xlim,
     ylim, zlim, xlim.scaled, ylim.scaled, zlim.scaled) {
     scaled.val <- function(x, original, scaled) {
         scaled[1] + (x - original[1]) * diff(scaled)/diff(original)
     }
     m <- ltransform3dto3d(rbind(scaled.val(state.map$x,
         xlim, xlim.scaled), scaled.val(state.map$y, ylim,
         ylim.scaled), zlim.scaled[1]), rot.mat, distance)
     panel.lines(m[1, ], m[2, ], col = "grey76")
 }

 pl <-cloud(density ~ long + lat, state.info, panel.3d.cloud = function(...) {
     panel.3dmap(...)
    panel.3dscatter(...)
 }, type = "h", scales = list(draw = TRUE), zoom = 1.1,
     xlim = state.map$range[1:2], ylim = state.map$range[3:4],
     xlab = NULL, ylab = NULL, zlab = NULL, aspect = c(diff(state.map$range[3:4])/diff(state.map$range[1:2]),
         0.3), panel.aspect = 0.75, lwd = 2, screen = list(z = 30,
         x = -60), par.settings = list(axis.line = list(col = "transparent"),
         box.3d = list(col = "transparent", alpha = 0)))

 print(pl)
 par(new=TRUE)
 plot.new()
 legend("topleft",legend=SurvYR[y],bty='n')
}
if(PdfOut>0)
 dev.off()

#=================================================
#====  bootstrap 95 ci for density================
#================================================
library(boot)
bM95CI<-matrix(nrow=4,ncol=length(SurvYR))
bM05CI<-matrix(nrow=4,ncol=length(SurvYR))
bootMeanMbig<-matrix(nrow=4,ncol=length(SurvYR))

for(d in 1:length(SurvYR))
 for(f in 1:4)
  {
 dataM<-DnsAtStnMbig[f,,d,1]
 dataM<-(dataM[-which(is.na(dataM))] )

# Statistic (MLE)
 samp.mean = function(dat) return(mean(dat))

# Bootstrap
 boots.outM <- boot(data=dataM, statistic=function(d, ind){samp.mean(d[ind])}, R = 10000)
 bootMeanMbig[f,d]<-boots.outM$t0

# 4 types of Bootstrap confidence intervals
# tempMci<-boot.ci(boots.outM, conf = 0.95, type = "norm")$normal

 tempMci<-boot.ci(boots.outM, conf = 0.95, type = "basic")$basic

 bM95CI[f,d]<-tempMci[5]
 bM05CI[f,d]<-tempMci[4]
  }


#==total numbers in an area by year
SurveyBigMdens<-apply(DnsAtStnMbig,c(1,3),mean,na.rm=T)
totSurvMbig<-matrix(nrow=4,ncol=length(SurvYR))
BootTotSurvMbig<-matrix(nrow=4,ncol=length(SurvYR))
for(j in 1:length(SurvYR))
 for(k in 1:nrow(totSurvMbig))
 {
  totSurvMbig[k,j]<-sum(propStFin[k]*nmiSurv*SurveyBigMdens[k,j])
  BootTotSurvMbig[k,j]<-sum(propStFin[k]*nmiSurv*bootMeanMbig[k,j])
 }
 
dev.new()
par(mfrow=c(1,2),mar=c(2,2,2,2))
plot(totSurvMbig[1,]~SurvYR,type="l",ylim=c(0,max(totSurvMbig)))
for(s in 2:4)
 lines(totSurvMbig[s,]~SurvYR,lty=s,col=s)
plot(BootTotSurvMbig[1,]~SurvYR,type="l",ylim=c(0,max(totSurvMbig)))
for(s in 2:4)
 lines(BootTotSurvMbig[s,]~SurvYR,lty=s,col=s)

#==proportion of crab in an area by year
PropInArea<-matrix(nrow=4,ncol=length(SurvYR))
   totCrab<-apply(totSurvMbig,2,sum)
for(j in 1:length(SurvYR))
 for(k in 1:nrow(BootTotSurvMbig))
   PropInArea[k,j]<-BootTotSurvMbig[k,j]/totCrab[j]


#==find the centroid of the fishery==
SurveyQuadCent<-rep(0,ncol(SurveyBigMaleCentroid))
for(i in 1:length(SurveyQuadCent))
 {
  if(SurveyBigMaleCentroid[2,i]<=58.5 & SurveyBigMaleCentroid[3,i] >=101){SurveyQuadCent[i]<-3}
  if(SurveyBigMaleCentroid[2,i]<=58.5 & SurveyBigMaleCentroid[3,i] <101){SurveyQuadCent[i]<-4}
  if(SurveyBigMaleCentroid[2,i]>58.5 & SurveyBigMaleCentroid[3,i] >=101){SurveyQuadCent[i]<-1}
  if(SurveyBigMaleCentroid[2,i]>58.5 & SurveyBigMaleCentroid[3,i] <101){SurveyQuadCent[i]<-2}
}
#==matches a lat and long with a depth
LatLongDep<-matrix(0,nrow=2,ncol=3)
counter<-1
for(y in 1:length(SurvYR))
{

yrDAT<-survDAT[drvYear==SurvYR[y],]

 for(x in 1:nrow(yrDAT))
 {
  if(yrDAT$MID_LATITUDE[x] != LatLongDep[counter,1] & yrDAT$MID_LONGITUDE[x] != LatLongDep[counter,2])
   {
    LatLongDep<-rbind(LatLongDep,c(yrDAT$MID_LATITUDE[x],yrDAT$MID_LONGITUDE[x],yrDAT$BOTTOM_DEPTH[x]))
    counter<-counter+1
   }
 }
}

# find catch centroids and depth
CatYR<-seq(1985,2011)
CatchBigMaleCentroid<-matrix(nrow=3,ncol=length(CatYR))

for(y in 1:length(CatYR))
{
fileyr<-CatYR[y]
file<-paste("C:/Users/Cody Szuwalski/Desktop/Spatial MSE/SnowCatchCSV/",fileyr,"snowcatch.csv",sep="")
catDAT<-read.csv(file)

long<-as.numeric(substring(catDAT$Stat.Area,1,2))
lat<-as.numeric(substring(catDAT$Stat.Area,3,6))
catDAT$Round.Pounds[is.na(catDAT$Round.Pounds)]<-0
CatchBigMaleCentroid[1,y]<--1*(weighted.mean(long,catDAT$Round.Pounds,na.rm=T)+100)
CatchBigMaleCentroid[2,y]<-weighted.mean(lat,catDAT$Round.Pounds,na.rm=T)/100

sortedline<-LatLongDep[order(LatLongDep[,1]),]
sortedline[,2]<-as.numeric(substring(as.character(sortedline[,2]),3))

depthLAT<-rep(0,nrow(catDAT))
for(k in 1:length(depthLAT))
{
#long<-long[-which(is.na(long))]
#lat<-lat[-which(is.na(lat))]
 testlat<-as.numeric(substring(lat[k],1,2))
 testmin<-as.numeric(substring(lat[k],3,4))
 maxlat<-testlat+1
 minlat<-testlat
 if(!is.na(testmin) & testmin>0) {minlat<-testlat+0.5}
 if(!is.na(testmin) & testmin==0){maxlat<-testlat+0.5}

 minlong<-as.numeric(long[k])
 maxlong<-as.numeric(long[k])+1

whichlong<- which(sortedline[,2]<maxlong & sortedline[,2]>minlong)
whichlat<- which(sortedline[,1]<maxlat& sortedline[,1]>minlat)
depthLAT[k]<-mean(sortedline[intersect(whichlong,whichlat),3])
if(is.na(testmin)){depthLAT[k]<-NA}
}

CatchBigMaleCentroid[3,y]<-weighted.mean(depthLAT,catDAT$Round.Pounds,na.rm=T)
}

CatchQuadCent<-rep(0,ncol(CatchBigMaleCentroid))
for(i in 1:length(CatchQuadCent))
 {
  if(CatchBigMaleCentroid[2,i]<=58.5 & CatchBigMaleCentroid[3,i] >=101){CatchQuadCent[i]<-3}
  if(CatchBigMaleCentroid[2,i]<=58.5 & CatchBigMaleCentroid[3,i] <101){CatchQuadCent[i]<-4}
  if(CatchBigMaleCentroid[2,i]>58.5 & CatchBigMaleCentroid[3,i] >=101){CatchQuadCent[i]<-1}
  if(CatchBigMaleCentroid[2,i]>58.5 & CatchBigMaleCentroid[3,i] <101){CatchQuadCent[i]<-2}
}
M100<-LatLongDep[LatLongDep[,3]<101 & LatLongDep[,3]>99,1:2]
M100s<-M100[order(M100[,1],M100[,2]),]
ernst<-read.csv("C:/Users/Cody Szuwalski/Desktop/Spatial MSE/ErnstMovement.csv")
ernst<-rbind(ernst,rep(0,3),rep(0,3))
MovePatterns<-c(ernst[1:7,2],SurveyQuadCent)

#correlation of catch centroid with ice edge latitude
IceData<-read.csv("C:/Users/Cody Szuwalski/Desktop/Spatial MSE/Ice Cover Index.csv")
cor(CatchBigMaleCentroid[2,1:25],IceData[,2])
cor(CatchBigMaleCentroid[2,1:25],IceData[,3])

plot(CatchBigMaleCentroid[2,1:25]~IceData[,2])

mod<-lm(CatchBigMaleCentroid[2,1:25] ~ IceData[,2]*IceData[,3])
summary(mod)
anova(mod)

mod<-lm(CatchBigMaleCentroid[2,1:25] ~ IceData[,2]+IceData[,3])
summary(mod)
anova(mod)

#north south distance traveled
mean((SurveyBigMaleCentroid[2,8:34]-CatchBigMaleCentroid[2,1:27])*40/(4*30))

# the ice extent is significantly correlated with the centroid of catch
#=========================================
# find the new distribution assuming equal dispersion and movement to the fishery
# not used as of now
#=========================================
ObsPropMove<-matrix(nrow=length(seq(8,length(SurvYR))),ncol=16)
for(y in 8:length(SurvYR))
{
yrDAT<-survDAT[drvYear==SurvYR[y],]
fileyr<-SurvYR[y]
BigMales<-yrDAT[yrDAT$SEX==1  & !is.na(yrDAT$WIDTH)& yrDAT$WIDTH<160 & yrDAT$WIDTH>101,]

SWshelf<-BigMales[BigMales$BOTTOM_DEPTH>=101 & BigMales$MID_LATITUDE<=58.5,]
NWshelf<-BigMales[BigMales$BOTTOM_DEPTH>=101 & BigMales$MID_LATITUDE>58.5,]
SEshelf<-BigMales[BigMales$BOTTOM_DEPTH<101 & BigMales$MID_LATITUDE<=58.5,]
NEshelf<-BigMales[BigMales$BOTTOM_DEPTH<101 & BigMales$MID_LATITUDE>58.5,]

SWshlLat<-SWshelf$MID_LATITUDE - (SurveyBigMaleCentroid[2,y] - CatchBigMaleCentroid[2,(y-7)])
SWshlLon<-SWshelf$MID_LONGITUDE - (SurveyBigMaleCentroid[1,y] - CatchBigMaleCentroid[1,(y-7)])
NWshlLat<-NWshelf$MID_LATITUDE - (SurveyBigMaleCentroid[2,y] - CatchBigMaleCentroid[2,(y-7)])
NWshlLon<-NWshelf$MID_LONGITUDE - (SurveyBigMaleCentroid[1,y] - CatchBigMaleCentroid[1,(y-7)])
SEshlLat<-SEshelf$MID_LATITUDE - (SurveyBigMaleCentroid[2,y] - CatchBigMaleCentroid[2,(y-7)])
SEshlLon<-SEshelf$MID_LONGITUDE - (SurveyBigMaleCentroid[1,y] - CatchBigMaleCentroid[1,(y-7)])
NEshlLat<-NEshelf$MID_LATITUDE - (SurveyBigMaleCentroid[2,y] - CatchBigMaleCentroid[2,(y-7)])
NEshlLon<-NEshelf$MID_LONGITUDE - (SurveyBigMaleCentroid[1,y] - CatchBigMaleCentroid[1,(y-7)])

tempName<-c("NW","NE","SW","SE")

for(z in 1:length(tempName))
{
  shelfLt<-get(paste(tempName[z],"shlLat",sep=""))
  shelfLn<-get(paste(tempName[z],"shlLon",sep=""))
  tempDepth<-rep(0,length(shelfLt))
   origShelf_Samp<-get(paste(tempName[z],"shelf",sep=""))$SAMPLING_FACTOR
  for(x in 1:length(shelfLt))
   {
   dist<-sqrt((shelfLt[x]-sortedline[,1])^2 + (shelfLn[x]-(-1*(100+sortedline[,2])))^2)
   tempDepth[x]<-sortedline[which(dist==min(dist,na.rm=T)),3]
   }
 tempFull<-cbind(shelfLt,shelfLn,tempDepth)

ToSE<-sum((tempFull[,1]<=58.5 & tempFull[,3]<=100)*origShelf_Samp,na.rm=T)/sum(origShelf_Samp)
ToSW<-sum((tempFull[,1]<=58.5 & tempFull[,3]>100)*origShelf_Samp,na.rm=T)/sum(origShelf_Samp)
ToNW<-sum((tempFull[,1]>58.5 & tempFull[,3]>100)*origShelf_Samp,na.rm=T)/sum(origShelf_Samp)
ToNE<-sum((tempFull[,1]>58.5 & tempFull[,3]<=100)*origShelf_Samp,na.rm=T)/sum(origShelf_Samp)
ObsPropMove[(y-7),(4*(z-1))+1]<-ToNW
ObsPropMove[(y-7),(4*(z-1))+2]<-ToNE
ObsPropMove[(y-7),(4*(z-1))+3]<-ToSW
ObsPropMove[(y-7),(4*(z-1))+4]<-ToSE
}
}

output<-round(ObsPropMove,2)
colnames(output)<-c("NW:NW","NW:NE","NW:SW","NW:SE","NE:NW","NE:NE","NE:SW","NE:SE",
"SW:NW","SW:NE","SW:SW","SW:SE","SE:NW","SE:NE","SE:SW","SE:SE")
rownames(output)<-seq(1985,2012)

write.csv(output,"C:/Users/Cody Szuwalski/Desktop/MovementProb.csv")

#=========================================
# find the new distribution assuming equal dispersion and movement to the fishery
# MOVEMENT BACK TO SURVEY GROUNDS
#=========================================
ObsPropMoveAfter<-matrix(nrow=length(seq(8,length(SurvYR))),ncol=16)
tempName<-c("NW","NE","SW","SE")

for(y in 8:length(SurvYR))
{
yrDAT<-survDAT[drvYear==SurvYR[y],]
fileyr<-SurvYR[y]
BigMales<-yrDAT[yrDAT$SEX==1  & !is.na(yrDAT$WIDTH)& yrDAT$WIDTH<160 & yrDAT$WIDTH>101,]

SWshelf<-BigMales[BigMales$BOTTOM_DEPTH>=101 & BigMales$MID_LATITUDE<=58.5,]
NWshelf<-BigMales[BigMales$BOTTOM_DEPTH>=101 & BigMales$MID_LATITUDE>58.5,]
SEshelf<-BigMales[BigMales$BOTTOM_DEPTH<101 & BigMales$MID_LATITUDE<=58.5,]
NEshelf<-BigMales[BigMales$BOTTOM_DEPTH<101 & BigMales$MID_LATITUDE>58.5,]

#TTHIS IS WRONG???
SWshlLat<-SWshelf$MID_LATITUDE + (SurveyBigMaleCentroid[2,y+1] - CatchBigMaleCentroid[2,(y-7)])
SWshlLon<-SWshelf$MID_LONGITUDE + (SurveyBigMaleCentroid[1,y+1] - CatchBigMaleCentroid[1,(y-7)])
NWshlLat<-NWshelf$MID_LATITUDE + (SurveyBigMaleCentroid[2,y+1] - CatchBigMaleCentroid[2,(y-7)])
NWshlLon<-NWshelf$MID_LONGITUDE + (SurveyBigMaleCentroid[1,y+1] - CatchBigMaleCentroid[1,(y-7)])
SEshlLat<-SEshelf$MID_LATITUDE + (SurveyBigMaleCentroid[2,y+1] - CatchBigMaleCentroid[2,(y-7)])
SEshlLon<-SEshelf$MID_LONGITUDE + (SurveyBigMaleCentroid[1,y+1] - CatchBigMaleCentroid[1,(y-7)])
NEshlLat<-NEshelf$MID_LATITUDE + (SurveyBigMaleCentroid[2,y+1] - CatchBigMaleCentroid[2,(y-7)])
NEshlLon<-NEshelf$MID_LONGITUDE + (SurveyBigMaleCentroid[1,y+1] - CatchBigMaleCentroid[1,(y-7)])

for(z in 1:length(tempName))
 {
   shelfLt<-get(paste(tempName[z],"shlLat",sep=""))
   shelfLn<-get(paste(tempName[z],"shlLon",sep=""))
   origShelf_Samp<-get(paste(tempName[z],"shelf",sep=""))$SAMPLING_FACTOR
   tempDepth<-rep(0,length(shelfLt))
  for(x in 1:length(shelfLt))
   {
   dist<-sqrt((shelfLt[x]-sortedline[,1])^2 + (shelfLn[x]-(-1*(100+sortedline[,2])))^2)
   tempDepth[x]<-sortedline[which(dist==min(dist,na.rm=T)),3]
   }
  tempFull<-cbind(shelfLt,shelfLn,tempDepth)

  ToSE<-sum((tempFull[,1]<=58.5 & tempFull[,3]<=100)*origShelf_Samp,na.rm=T)/sum(origShelf_Samp)
  ToSW<-sum((tempFull[,1]<=58.5 & tempFull[,3]>100)*origShelf_Samp,na.rm=T)/sum(origShelf_Samp)
  ToNW<-sum((tempFull[,1]>58.5 & tempFull[,3]>100)*origShelf_Samp,na.rm=T)/sum(origShelf_Samp)
  ToNE<-sum((tempFull[,1]>58.5 & tempFull[,3]<=100)*origShelf_Samp,na.rm=T)/sum(origShelf_Samp)
  ObsPropMoveAfter[(y-7),(4*(z-1))+1]<-ToNW
  ObsPropMoveAfter[(y-7),(4*(z-1))+2]<-ToNE
  ObsPropMoveAfter[(y-7),(4*(z-1))+3]<-ToSW
  ObsPropMoveAfter[(y-7),(4*(z-1))+4]<-ToSE
}
}

outputAft<-round(ObsPropMoveAfter,2)
colnames(outputAft)<-c("NW:NW","NW:NE","NW:SW","NW:SE","NE:NW","NE:NE","NE:SW","NE:SE",
"SW:NW","SW:NE","SW:SW","SW:SE","SE:NW","SE:NE","SE:SW","SE:SE")
rownames(outputAft)<-seq(1985,2012)
library(rgl)
write.csv(outputAft,"C:/Users/Cody Szuwalski/Desktop/MovementProbAfter.csv")



yrng<-seq(57.9,79.1,0.1)
xrng<-seq(54.5,62.7,0.1)
zmat<-matrix(ncol=length(xrng),nrow=length(yrng))
for(d in 1:length(yrng))
 for(f in 1:length(xrng))
  {
   tmpS<-sortedline[sortedline[,1]>xrng[f] & sortedline[,1]<xrng[f+1] & sortedline[,2]>yrng[d] & sortedline[,2]<yrng[d+1],]
   if(length(tmpS)>0)
   {
    try(zmat[d,f]<-(tmpS[3]),TRUE) 
    try(zmat[d,f]<-tmpS[,3],TRUE)
   }
  }

filled.contour(x=xrng,y=yrng,z=zmat)

library(maps)
library(PBSmapping)
library(maptools)
library(maps)       #basic mapping functions and some data
library(mapdata)    #some additional hires data
library(maptools)   #useful tools such as reading shapefiles
library(mapproj)

#====================
#== MOVEMENT FROM SURVEY TO FISHERY
#======================
par(mfrow=c(1,2),mar=c(2,2,2,2),oma=c(3,3,2,1))
#map('worldHires', xlim=c(-175, -160), ylim=c(56, 61))
plot(0,xlim=c(-174,-166.5),ylim=c(56,59.5),ylab="Latitude",xlab="Longitude")
for(w in 1:ncol(CatchBigMaleCentroid))
 points(CatchBigMaleCentroid[2,w]~CatchBigMaleCentroid[1,w],pch=21,bg='red')
for(w in 1:ncol(SurveyBigMaleCentroid))
 points(SurveyBigMaleCentroid[2,w]~SurveyBigMaleCentroid[1,w],pch=21,bg='blue')
for(w in 1:ncol(CatchBigMaleCentroid))
 points(CatchBigMaleCentroid[2,w]~CatchBigMaleCentroid[1,w],pch=21,bg='red')

for(w in 1:ncol(CatchBigMaleCentroid))
 lines(y=c(CatchBigMaleCentroid[2,w],SurveyBigMaleCentroid[2,w+7]),x=c(CatchBigMaleCentroid[1,w],SurveyBigMaleCentroid[1,w+7]))
thinObs<-seq(1,nrow(M100s),1)
lines(loess(M100s[thinObs,1]~M100s[thinObs,2]),lty=2,col=2)
#points(M100[,1]~M100[,2])
#lines(M100s[,1]~M100s[,2],lty=2,col=2)
library(calibrate)
textxy(CatchBigMaleCentroid[1,],CatchBigMaleCentroid[2,],CatYR,cx=.61)
textxy(SurveyBigMaleCentroid[1,],SurveyBigMaleCentroid[2,],SurvYR,cx=.61)
abline(h=58.5,lty=2,col=2)
#legend("topright",col=c('blue','red','black','red'),pch=c(16,16,NA,NA),
#lty=c(NA,NA,1,2), bty='n',legend=c("Survey","Fishery","Movement","Area Division"),
#bg=c(2,4,NA,NA))

mtext("Survey to Fishery",side=3,cex=1.5,line=.5)

#====================
#== MOVEMENT FROM FISHERY TO SURVEY
#======================
#map('worldHires', xlim=c(-175, -160), ylim=c(56, 61))
plot(0,xlim=c(-174,-166.5),ylim=c(56,59.5),ylab="Latitude",xlab="Longitude")
for(w in 1:ncol(CatchBigMaleCentroid))
 points(CatchBigMaleCentroid[2,w]~CatchBigMaleCentroid[1,w],pch=21,bg='red')
for(w in 1:ncol(SurveyBigMaleCentroid))
 points(SurveyBigMaleCentroid[2,w]~SurveyBigMaleCentroid[1,w],pch=21,bg='blue')
for(w in 1:ncol(CatchBigMaleCentroid))
 points(CatchBigMaleCentroid[2,w]~CatchBigMaleCentroid[1,w],pch=21,bg='red')

for(w in 1:ncol(CatchBigMaleCentroid))
 lines(y=c(CatchBigMaleCentroid[2,w],SurveyBigMaleCentroid[2,w+8]),x=c(CatchBigMaleCentroid[1,w],SurveyBigMaleCentroid[1,w+8]))
thinObs<-seq(1,nrow(M100s),1)
lines(loess(M100s[thinObs,1]~M100s[thinObs,2]),lty=2,col=2)
#points(M100[,1]~M100[,2])
#lines(M100s[,1]~M100s[,2],lty=2,col=2)
library(calibrate)
textxy(CatchBigMaleCentroid[1,],CatchBigMaleCentroid[2,],CatYR,cx=.61)
textxy(SurveyBigMaleCentroid[1,],SurveyBigMaleCentroid[2,],SurvYR,cx=.61)
abline(h=58.5,lty=2,col=2)
legend("topright",col=c('blue','red','black','red'),pch=c(16,16,NA,NA),
lty=c(NA,NA,1,2), bty='n',legend=c("Survey","Fishery","Movement","Area Division"),
bg=c(2,4,NA,NA))

mtext("Fishery to Survey",side=3,cex=1.5,line=.5)
mtext("Longitude",side=1,outer=T,cex=1.5,line=.5)
mtext("Latitude",side=2,outer=T,cex=1.5,line=.5)
#============================================
#=== Figure for demonstrating assumptions of movement
#===========================================
23, 27
for(y in 8:length(SurvYR))
{
dev.new(width=9,height=5)
temp<-DnsAtStnMbig[,,y,]
fileyr<-CatYR[y-7]
file<-paste("C:/Users/Cody Szuwalski/Desktop/Spatial MSE/SnowCatchCSV/",fileyr,"snowcatch.csv",sep="")
catDAT<-read.csv(file)
#apply(DnsAtStnMbig,c(1,3),mean,na.rm=T)
#temp[,,1]


# aggreagates duplicate stat ares
long<-as.numeric(substring(catDAT$Stat.Area,1,2))
lat<-as.numeric(substring(catDAT$Stat.Area,3,6))
modlong<--1*(100+long)
modlat<-lat/100
tempCatch<-cbind(modlong,modlat,catDAT$Round.Pounds)

plotLoc<-unique(paste(modlong,modlat))
unqLoc<-unique(cbind(modlong,modlat))
rawLoc<-paste(modlong,modlat)
plotCatch<-rep(0,length(plotLoc))
for(w in 1:length(plotLoc))
{
 selCat<-match(plotLoc[w],rawLoc)
 for(r in 1:length(selCat))
  plotCatch[w]<-plotCatch[w]+tempCatch[selCat[r],3]
}


meanDensYr<-matrix(ncol=200,nrow=4)
for(i in 1:4)
 for(j in 1:200)
   meanDensYr[i,j]<-mean(temp[i,j,],na.rm=T)

#pairs mean density with a latitude, longitude and station   
area1<-cbind(meanDensYr[1,],latsLonsSta[1,,y,])
area2<-cbind(meanDensYr[2,],latsLonsSta[2,,y,])
area3<-cbind(meanDensYr[3,],latsLonsSta[3,,y,])
area4<-cbind(meanDensYr[4,],latsLonsSta[4,,y,])
area1<-area1[!is.na(area1[,1]),]
area2<-area2[!is.na(area2[,1]),]
area3<-area3[!is.na(area3[,1]),]
area4<-area4[!is.na(area4[,1]),]

#aggregates and organizes for plotting
totAreas<-rbind(area1,area2,area3,area4)
latGrid<-seq(min(totAreas[,2]),max(totAreas[,2])+0.01,0.5)
lonGrid<-seq(min(totAreas[,3]),max(totAreas[,3])+0.01,0.5)
PlotZmat<-matrix(nrow=length(latGrid),ncol=length(lonGrid))
for(x in 1:length(latGrid))
 for(z in 1:length(lonGrid))
   PlotZmat[x,z]<-log(sum(totAreas[totAreas[,2]>=latGrid[x] & totAreas[,2]<latGrid[x+1] & totAreas[,3]>=lonGrid[z] & totAreas[,3]<lonGrid[z+1],1]) +0.01)

#filled contour maps
plot.new()
par(new="TRUE",plt=c(0.1,0.5,.15,.85))
filled.contour3(x=lonGrid,y=latGrid,z=t((PlotZmat)),ylim=c(55,60),xlim=c(-177.5,-162.5),axes=F)
axis(1)
axis(2)
lines(loess(M100s[thinObs,1]~M100s[thinObs,2]),lty=2,col=1,lwd=3)
abline(h=58.5,lty=2,col=1,lwd=3)
mtext(side=3,"Survey",cex=2,line=.25)
mtext(side=2,"Latitude",cex=2,line=2)
legend("topright",legend=SurvYR[y],bty='n')
#get a vector for the differences
SurvCatchDiffLat<-SurveyBigMaleCentroid[2,y] - CatchBigMaleCentroid[2,(y-7)]
SurvCatchDiffLon<-SurveyBigMaleCentroid[2,y] - CatchBigMaleCentroid[2,(y-7)]

NewLatGrid<-latGrid-SurvCatchDiffLat
NewLonGrid<-lonGrid-SurvCatchDiffLon

par(new=TRUE)
symbols(x=unqLoc[,1],y=unqLoc[,2],circles=log(plotCatch),inches=1/10,ylim=c(55,60),xlim=c(-177.5,-162.5),ylab="",xlab="",yaxt='n',xaxt='n')


par(new="TRUE",plt=c(0.55,0.95,.15,.85))
filled.contour3(x=NewLonGrid,y=NewLatGrid,z=t((PlotZmat)),ylim=c(55,60),xlim=c(-177.5,-162.5),axes=F)
axis(1)
lines(loess(M100s[thinObs,1]~M100s[thinObs,2]),lty=2,col=1,lwd=3)
abline(h=58.5,lty=2,col=1,lwd=3)
mtext(side=3,"Fishery",cex=2,line=.25)
mtext(side=1,outer=T,"Longitude",cex=2,line=-1.5,adj=.54)

par(new=TRUE)
#symbols(x=modlong,y=modlat,circles=log(catDAT$Round.Pounds),inches=1/5,ylim=c(55,60),xlim=c(-177.5,-162.5),axes=F,ylab="",xlab="")
symbols(x=unqLoc[,1],y=unqLoc[,2],circles=log(plotCatch),inches=1/8,ylim=c(55,60),xlim=c(-177.5,-162.5),ylab="",xlab="",yaxt='n',xaxt='n')

}
#=============================================
# FIGURE 3IN MANUSCRIPT
#=============================

source("http://wiki.cbr.washington.edu/qerm/sites/qerm/images/2/25/Filled.legend.R")
dev.new(width=9,height=9)
plot.new()
for(y in c(23,27))
{

temp<-DnsAtStnMbig[,,y,]
fileyr<-CatYR[y-7]
file<-paste("C:/Users/Cody Szuwalski/Desktop/Spatial MSE/SnowCatchCSV/",fileyr,"snowcatch.csv",sep="")
catDAT<-read.csv(file)
#apply(DnsAtStnMbig,c(1,3),mean,na.rm=T)
#temp[,,1]


# aggreagates duplicate stat ares
long<-as.numeric(substring(catDAT$Stat.Area,1,2))
lat<-as.numeric(substring(catDAT$Stat.Area,3,6))
modlong<--1*(100+long)
modlat<-lat/100
tempCatch<-cbind(modlong,modlat,catDAT$Round.Pounds)

plotLoc<-unique(paste(modlong,modlat))
unqLoc<-unique(cbind(modlong,modlat))
rawLoc<-paste(modlong,modlat)
plotCatch<-rep(0,length(plotLoc))
for(w in 1:length(plotLoc))
{
 selCat<-match(plotLoc[w],rawLoc)
 for(r in 1:length(selCat))
  plotCatch[w]<-plotCatch[w]+tempCatch[selCat[r],3]
}


meanDensYr<-matrix(ncol=200,nrow=4)
for(i in 1:4)
 for(j in 1:200)
   meanDensYr[i,j]<-mean(temp[i,j,],na.rm=T)

#pairs mean density with a latitude, longitude and station   
area1<-cbind(meanDensYr[1,],latsLonsSta[1,,y,])
area2<-cbind(meanDensYr[2,],latsLonsSta[2,,y,])
area3<-cbind(meanDensYr[3,],latsLonsSta[3,,y,])
area4<-cbind(meanDensYr[4,],latsLonsSta[4,,y,])
area1<-area1[!is.na(area1[,1]),]
area2<-area2[!is.na(area2[,1]),]
area3<-area3[!is.na(area3[,1]),]
area4<-area4[!is.na(area4[,1]),]

#aggregates and organizes for plotting
totAreas<-rbind(area1,area2,area3,area4)
latGrid<-seq(min(totAreas[,2]),max(totAreas[,2])+0.01,0.5)
lonGrid<-seq(min(totAreas[,3]),max(totAreas[,3])+0.01,0.5)
PlotZmat<-matrix(nrow=length(latGrid),ncol=length(lonGrid))
for(x in 1:length(latGrid))
 for(z in 1:length(lonGrid))
   PlotZmat[x,z]<-log(sum(totAreas[totAreas[,2]>=latGrid[x] & totAreas[,2]<latGrid[x+1] & totAreas[,3]>=lonGrid[z] & totAreas[,3]<lonGrid[z+1],1]) +0.01)

#filled contour maps

if(y==27)par(new="TRUE",plt=c(0.075,0.475,.1,.5))
if(y==23)par(new="TRUE",plt=c(0.075,0.475,.55,.95))

filled.contour3(x=lonGrid,y=latGrid,z=t((PlotZmat)),ylim=c(55,59.35),xlim=c(-177.5,-164),axes=F)
if(y==27)axis(1)
axis(2)
lines(loess(M100s[thinObs,1]~M100s[thinObs,2]),lty=2,col=1,lwd=3)
abline(h=58.5,lty=2,col=1,lwd=3)
if(y==23)mtext(side=3,"Survey",cex=1.75,line=.25)
mtext(side=2,"Latitude",cex=1.75,line=-1.1,outer=T)

#get a vector for the differences
SurvCatchDiffLat<-SurveyBigMaleCentroid[2,y] - CatchBigMaleCentroid[2,(y-7)]
SurvCatchDiffLon<-SurveyBigMaleCentroid[2,y] - CatchBigMaleCentroid[2,(y-7)]

NewLatGrid<-latGrid-SurvCatchDiffLat
NewLonGrid<-lonGrid-SurvCatchDiffLon

par(new=TRUE)
symbols(x=unqLoc[,1],y=unqLoc[,2],circles=log(plotCatch),inches=1/8,ylim=c(55,59.35),xlim=c(-177.5,-164),ylab="",xlab="",yaxt='n',xaxt='n')


if(y==27)par(new="TRUE",plt=c(0.525,0.925,.1,.5))
if(y==23)par(new="TRUE",plt=c(0.525,0.925,.55,.95))
filled.contour3(x=NewLonGrid,y=NewLatGrid,z=t((PlotZmat)),ylim=c(55,59.35),xlim=c(-177.5,-164),axes=F)
if(y==27)axis(1)
lines(loess(M100s[thinObs,1]~M100s[thinObs,2]),lty=2,col=1,lwd=3)
abline(h=58.5,lty=2,col=1,lwd=3)
if(y==23)mtext(side=3,"Fishery",cex=1.75,line=.25)
mtext(side=1,outer=T,"Longitude",cex=1.75,line=-1.5,adj=.54)
legend("topright",legend=SurvYR[y],bty='n')
par(new=TRUE)
symbols(x=unqLoc[,1],y=unqLoc[,2],circles=log(plotCatch),inches=1/8,ylim=c(55,59.35),xlim=c(-177.5,-164),ylab="",xlab="",yaxt='n',xaxt='n')

}

par(new="TRUE",plt=c(.93,.95,0.25,0.75))
filled.legend(x=NewLonGrid,y=NewLatGrid,z=t((PlotZmat)),zlim=c(-4.6,9.8))
mtext(side=3,"log(Density)",line=1)




#given the beginning area of the centroid, the temperatures (maybe these should be
#based on the whole area instead of weighted by where there are crab?)
# in the adjacent areas, what direction does the centroid move? what about
# after the fishery?

 library(rgl)
 hist3d(rnorm(2500),rnorm(2500),alpha=0.4,nclass=7,scale=30)
 hist3d(BigMales$MID_LATITUDE,BigMales$MID_LONGITUDE,alpha=0.4,nclass=7,scale=30)
#this will give three? different movement patterns in each year, for each movement
#

#fit to the centroid of abundances of large males in the objective function?

#is the slope of the line between the centroids correlated to the cold pool?
# take the centroid and spread of large males at survey
#move it to the centroid at the fishery
# find the proportion of the population ineach quadrant, regardless 
#of the spread at the time of fishery--this accounts for the fact 
#that the ice edge limits the fishery?

#take the survey data, add/subtract the centroid from it, recalculate the proport

filled.contour3 <-
  function (x = seq(0, 1, length.out = nrow(z)),
            y = seq(0, 1, length.out = ncol(z)), z, xlim = range(x, finite = TRUE), 
            ylim = range(y, finite = TRUE), zlim = range(z, finite = TRUE), 
            levels = pretty(zlim, nlevels), nlevels = 20, color.palette = cm.colors, 
            col = color.palette(length(levels) - 1), plot.title, plot.axes, 
            key.title, key.axes, asp = NA, xaxs = "i", yaxs = "i", las = 1, 
            axes = TRUE, frame.plot = axes,mar, ...) 
{
  # modification by Ian Taylor of the filled.contour function
  # to remove the key and facilitate overplotting with contour()
  # further modified by Carey McGilliard and Bridget Ferris
  # to allow multiple plots on one page

  if (missing(z)) {
    if (!missing(x)) {
      if (is.list(x)) {
        z <- x$z
        y <- x$y
        x <- x$x
      }
      else {
        z <- x
        x <- seq.int(0, 1, length.out = nrow(z))
      }
    }
    else stop("no 'z' matrix specified")
  }
  else if (is.list(x)) {
    y <- x$y
    x <- x$x
  }
  if (any(diff(x) <= 0) || any(diff(y) <= 0)) 
    stop("increasing 'x' and 'y' values expected")
 # mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
 # on.exit(par(par.orig))
 # w <- (3 + mar.orig[2]) * par("csi") * 2.54
 # par(las = las)
 # mar <- mar.orig
 plot.new()
 # par(mar=mar)
  plot.window(xlim, ylim, "", xaxs = xaxs, yaxs = yaxs, asp = asp)
  if (!is.matrix(z) || nrow(z) <= 1 || ncol(z) <= 1) 
    stop("no proper 'z' matrix specified")
  if (!is.double(z)) 
    storage.mode(z) <- "double"
  .Internal(filledcontour(as.double(x), as.double(y), z, as.double(levels), 
                          col = col))
  if (missing(plot.axes)) {
    if (axes) {
      title(main = "", xlab = "", ylab = "")
      Axis(x, side = 1)
      Axis(y, side = 2)
    }
  }
  else plot.axes
  if (frame.plot) 
    box()
  if (missing(plot.title)) 
    title(...)
  else plot.title
  invisible()
}