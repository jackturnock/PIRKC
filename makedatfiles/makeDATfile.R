MakeDATfile<-function(survDAT,numDAT,BinSize,firstbin,lastbin,USEBOB=1,bootCIout=0)
{
require(plotrix)
drvYear<-as.numeric(substr(survDAT$CRUISE,1,4))
SurvYR<-unique(drvYear)
yearRange<-SurvYR
LenBin1<-firstbin
FinBin<-lastbin

#==find the area of the PI surveyed grounds
RKCSta<-unique(survDAT[survDAT$MID_LONGITUDE< -168 & survDAT$MID_LATITUDE < 58.65,]$GIS_STATION)
#==just stations that have had observations of crab (22)
PIarea<-22*400
AllStation<-RKCSta

#===================================
#  Densities only Pribolof district (S of 58.39, W of 168)
#===================================
DensityMpi<-matrix(nrow=length(SurvYR),ncol=length(AllStation))
DensityFpi<-matrix(nrow=length(SurvYR),ncol=length(AllStation))
StationYrPi<-matrix(nrow=length(SurvYR),ncol=length(AllStation))
TotalObsCrabF<-matrix(nrow=length(SurvYR),ncol=length(AllStation))
TotalObsCrabM<-matrix(nrow=length(SurvYR),ncol=length(AllStation))
HaulsRec<-rep(0,length(SurvYR))

LengthBin<-seq(7.5,210,BinSize)
histBreaks<-LengthBin-BinSize/2
NatLengthF<-matrix(0,ncol=(length(LengthBin)-1),nrow=length(SurvYR))
NatLengthM<-matrix(0,ncol=(length(LengthBin)-1),nrow=length(SurvYR))
LenFreqF<-matrix(ncol=length(LengthBin)-1,nrow=length(SurvYR))
LenFreqM<-matrix(ncol=length(LengthBin)-1,nrow=length(SurvYR))

HasCatch<-NULL
HaulNums<-NULL

for(y in 1:length(SurvYR))
{
 yrDAT<-survDAT[drvYear==SurvYR[y],]
 yrDAT<-yrDAT[yrDAT$MID_LONGITUDE< -168 & yrDAT$MID_LATITUDE < 58.65,]
 fileyr<-SurvYR[y]
 stationsUNQ<-RKCSta
 yrDATm<-yrDAT[!is.na(yrDAT$SEX==1)&(yrDAT$SEX==1)=="TRUE" & yrDAT$LENGTH>35,]
 yrDATf<-yrDAT[!is.na(yrDAT$SEX==2)&(yrDAT$SEX==2)=="TRUE" & yrDAT$LENGTH>35,]

#==must deal with repeated hauls
#==loop over GIS STATIONS to create a length vector
#==year 36, b==3 for example.
 tempGIS<-unique(yrDATm$GIS_STATION)
 for(b in 1:length(tempGIS))
  {
  subDatm<-yrDATm[yrDATm$GIS_STATION==tempGIS[b],]
  subDatf<-yrDATf[yrDATf$GIS_STATION==tempGIS[b],]
  tempHAULm<-unique(subDatm$HAUL)
  tempHAULf<-unique(subDatm$HAUL)

#==if there are no repeated hauls, proceed as normal 
  if(length(tempHAULm)==1)
  {
   tempMhist<-weighted.hist(subDatm$LENGTH,breaks=histBreaks,w = subDatm$SAMPLING_FACTOR, plot=F)
   NatLengthM[y,]<-NatLengthM[y,]+tempMhist$counts
  }
  if(length(tempHAULf)==1)
  {
   tempFhist<-weighted.hist(subDatf$LENGTH,breaks=histBreaks,w = subDatf$SAMPLING_FACTOR, plot=F)
   NatLengthF[y,]<-NatLengthF[y,]+tempFhist$counts
  }
 #==males====
  if(length(tempHAULm)>1)
  {
   HaulSampNm<-rep(0,length(tempHAULm))
  #==find the avg number of samples==
  #==the number of observations in a haul * samplingfactor==
   for(q in 1:length(HaulSampNm))
    HaulSampNm[q]<-nrow(subDatm[subDatm$HAUL==tempHAULm[q],])*subDatm[subDatm$HAUL==tempHAULm[q],]$SAMPLING_FACTOR[1]
  #==find the length frequency ==
    tempMhist<-weighted.hist(subDatm$LENGTH,breaks=histBreaks,w = subDatm$SAMPLING_FACTOR, plot=F)
  #==mult avg sampN by length freq, add to total N at len (will give fractions)==
   if(any(!is.na(HaulSampNm)))
    NatLengthM[y,]<-NatLengthM[y,]+(tempMhist$density * mean(HaulSampNm,na.rm=T))
   }

  #==females==
  if(length(tempHAULf)>1)
  {
   HaulSampNf<-rep(0,length(tempHAULf))
   for(q in 1:length(HaulSampNf))
    HaulSampNf[q]<-nrow(subDatf[subDatf$HAUL==tempHAULf[q],])*subDatf[subDatf$HAUL==tempHAULf[q],]$SAMPLING_FACTOR[1]
    temphFist<-weighted.hist(subDatf$LENGTH,breaks=histBreaks,w = subDatf$SAMPLING_FACTOR, plot=F)
   if(any(!is.na(HaulSampNf)))
    NatLengthF[y,]<-NatLengthF[y,]+(tempFhist$density * mean(HaulSampNf,na.rm=T))
   }
  }
 DensityMpi[y,1:length(stationsUNQ)]<-0
 DensityFpi[y,1:length(stationsUNQ)]<-0
 TotalObsCrabM[y,1:length(stationsUNQ)]<-0
 TotalObsCrabF[y,1:length(stationsUNQ)]<-0

#==density at station==
 for(j in 1:length(stationsUNQ))
  {
   stationALL<-yrDAT[!is.na(match(yrDAT$GIS_STATION,stationsUNQ[j])),]
   StationYrPi[y,j]<-as.character(stationsUNQ[j])
   Hauls<-(unique(stationALL$HAUL))
   HaulsRec[y]<-HaulsRec[y]+sum(!is.na(Hauls))
   Males<-stationALL[stationALL$SEX==1,]
   Females<-stationALL[stationALL$SEX==2,]

 #==densities across hauls in crabs per km^2==
  tempDensM<-NULL
  tempDensF<-NULL
  tempTotF<-NULL
  tempTotM<-NULL

  for(k in 1:length(Hauls))
   {
    HaulNums<-append(HaulNums,length(Hauls))
    SampFactM<-Males$SAMPLING_FACTOR[which(Males$HAUL==Hauls[k])[1]] 
    AreaSweptM<-Males$AREA_SWEPT_VARIABLE[which(Males$HAUL==Hauls[k])[1]]
    tempDensM<-append(tempDensM,length(Males$HAUL==Hauls[k])*SampFactM/AreaSweptM)
    tempTotM<-append(tempTotM,length(Males$HAUL==Hauls[k])*SampFactM)
    tempLenM<-Males[Males$HAUL==Hauls[k],]$LENGTH

    SampFactF<-Females$SAMPLING_FACTOR[which(Females$HAUL==Hauls[k])[1]] 
    AreaSweptF<-Females$AREA_SWEPT_VARIABLE[which(Females$HAUL==Hauls[k])[1]]
    tempDensF<-append(tempDensF,length(Females$HAUL==Hauls[k])*SampFactF/AreaSweptF)
    tempTotF<-append(tempTotF,length(Females$HAUL==Hauls[k])*SampFactF)
    tempLenF<-Females[Females$HAUL==Hauls[k],]$LENGTH
    }

  if(!is.na(tempDensM)) 
  {
   DensityMpi[y,j]<-mean(tempDensM)
   TotalObsCrabM[y,j]<-sum(tempTotM)
   HasCatch<-append(HasCatch,as.character(stationsUNQ[j]))
  }

#==if there are multiple hauls at a station, use the average density==
  if(!is.na(tempDensF))
   {
   DensityFpi[y,j]<-mean(tempDensF)
   TotalObsCrabF[y,j]<-sum(tempTotF)
   HasCatch<-append(HasCatch,as.character(stationsUNQ[j]))
   }
  }
}

#==Calculate length frequencies==
LenFreqFtemp<-matrix(ncol=length(LengthBin)-1,nrow=length(SurvYR))
LenFreqMtemp<-matrix(ncol=length(LengthBin)-1,nrow=length(SurvYR))
for(x in 1:nrow(LenFreqM))
 {
  LenFreqMtemp[x,]<-NatLengthM[x,]/max(NatLengthM[x,],na.rm=T)
  LenFreqFtemp[x,]<-NatLengthF[x,]/max(NatLengthF[x,],na.rm=T)
 }

LenFreqF<-matrix(0,ncol=length(LengthBin)-1,nrow=length(SurvYR))
LenFreqM<-matrix(0,ncol=length(LengthBin)-1,nrow=length(SurvYR))
for(x in 1:nrow(LenFreqFtemp))
 {
  # bisection method to find a constant to multiply the natlen by
  # to make it a distribution
   a1<-0.00001
   a2<-2.1
   tempVec<-LenFreqMtemp[x,]
  if(!is.na(tempVec))
  {
   for(y in 1:30)
    {
     val1<-(a2+a1)/2
     chkIt<-sum(tempVec*val1)
     if(chkIt>1)
      a2<-val1
     if(chkIt<1)
      a1<-val1
     print(chkIt)
    }
    LenFreqM[x,]<-LenFreqMtemp[x,]*val1
   }
 #==for females now
   a1<-0.00001
   a2<-2.1
   tempVec<-LenFreqFtemp[x,]
   if(!is.na(tempVec))
   {
   for(y in 1:30)
    {
     val1<-(a2+a1)/2
     chkIt<-sum(tempVec*val1)
     if(chkIt>1)
      a2<-val1
     if(chkIt<1)
      a1<-val1
     print(chkIt)
    }
   LenFreqF[x,]<-LenFreqFtemp[x,]*val1
   }
}

#==total crab in year, total stations with crab
nSampMlenf<-apply(NatLengthM,1,sum,na.rm=T)
nSampFlenf<-apply(NatLengthF,1,sum,na.rm=T)
nSampStM<-apply(TotalObsCrabM>0,1,sum,na.rm=T)
nSampStF<-apply(TotalObsCrabF>0,1,sum,na.rm=T)

nSampMlenfUSE<-rep(200,length(nSampMlenf))
nSampFlenfUSE<-rep(200,length(nSampFlenf))
for(g in 1:length(nSampMlenfUSE))
{
 if(nSampMlenf[g]<200) nSampMlenfUSE[g]<-nSampMlenf[g]
 if(nSampFlenf[g]<200) nSampFlenfUSE[g]<-nSampFlenf[g]
}

#=================================================
#====  bootstrap 95 ci for density================
#================================================
if(bootCIout==1)
{
require(boot)
bM95CI<-rep(0,length(SurvYR))
bM05CI<-rep(0,length(SurvYR))
bootMeanMbig<-rep(0,length(SurvYR))

bF95CI<-rep(0,length(SurvYR))
bF05CI<-rep(0,length(SurvYR))
bootMeanFbig<-rep(0,length(SurvYR))
samp.mean = function(dat) return(mean(dat))

for(d in 1:length(SurvYR))
  {
 dataM<-0
 if(sum(DensityMpi[d,],na.rm=T)>0) 
 {
  dataM<-DensityMpi[d,]
  dataM<-(dataM[which(!is.na(dataM))] )

  boots.outM <- boot(data=dataM, statistic=function(d, ind){samp.mean(d[ind])}, R = 10000)
  bootMeanMbig[d]<-boots.outM$t0
  tempMci<-boot.ci(boots.outM, conf = 0.95, type = "all")$percent
  if(length(tempMci)>0)
  {
   bM95CI[d]<-tempMci[5]
   bM05CI[d]<-tempMci[4]
  }
 }

 dataF<-0
 if(sum(DensityFpi[d,],na.rm=T)>0) 
  {
   dataF<-DensityFpi[d,]
   dataF<-(dataF[which(!is.na(dataF))] )

  boots.outF <- boot(data=dataF, statistic=function(d, ind){samp.mean(d[ind])}, R = 10000)
  bootMeanFbig[d]<-boots.outF$t0
  tempFci<-boot.ci(boots.outF, conf = 0.95, type = "all")$percent
  if(length(tempFci)>0)
  {
  bF95CI[d]<-tempFci[5]
  bF05CI[d]<-tempFci[4]
  }
 }
}

write.csv(cbind(bF05CI,bF95CI,bM05CI,bM95CI)*PIarea,"CVs_SurveyN.csv")
}


#==select survey numbers and CV to use
if(USEBOB==1)
{
useSurvF<-numDat$ABUNDANCE_FEMALE_TOTAL
useSurvM<-numDat$ABUNDANCE_MALE_TOTAL
useFSurveyCV<-numDat$CV_ABUNDANCE_FEMALE_TOTAL
useMSurveyCV<-numDat$CV_ABUNDANCE_MALE_TOTAL
}
if(USEBOB==0)
{
YrDensM<-apply(DensityMpi,1,mean,na.rm=T)
YrDensF<-apply(DensityFpi,1,mean,na.rm=T)
useSurvF<-PIarea*YrDensF
useSurvM<-PIarea*YrDensM

#==bring the code over to calculate different CVs
useFSurveyCV<-numDat$CV_ABUNDANCE_FEMALE_TOTAL
useMSurveyCV<-numDat$CV_ABUNDANCE_MALE_TOTAL
}

#======================
#  make .DAT file
#=====================
#==natural mortality
M<-0.18
useFSurveyCV[11]<-1

#==Recruit to fishery: 135-149mm new shell
tmp<-which(apply(LenFreqM,2,sum,na.rm=T)==0)
cutLenBins<-tmp[length(tmp)]+1

#== .DAT file==
file.create("PIRKC.DAT")
cat("#Pribilof Island Red King Crab assessment","\n",file="PIRKC.DAT")
cat("#","\n",file="PIRKC.DAT",append=TRUE)
cat("#start/end year","\n",file="PIRKC.DAT",append=TRUE)
cat(c(min(yearRange)-5,max(yearRange)),"\n",file="PIRKC.DAT",append=TRUE)
cat("#number of length bins","\n",file="PIRKC.DAT",append=TRUE)
cat(ncol(LenFreqM[,cutLenBins:ncol(LenFreqM)]),"\n",file="PIRKC.DAT",append=TRUE)

cat("#survey years","\n",file="PIRKC.DAT",append=TRUE)
cat(length(yearRange),"\n",file="PIRKC.DAT",append=TRUE)
cat("#years for survey","\n",file="PIRKC.DAT",append=TRUE)
cat(yearRange,"\n",file="PIRKC.DAT",append=TRUE)
cat("#length freq from survey years","\n",file="PIRKC.DAT",append=TRUE)
cat(length(yearRange),"\n",file="PIRKC.DAT",append=TRUE)
cat("#years for length freq survey","\n",file="PIRKC.DAT",append=TRUE)
cat(yearRange,"\n",file="PIRKC.DAT",append=TRUE)
cat("#Male survey sample sizes","\n",file="PIRKC.DAT",append=TRUE)
cat(nSampMlenfUSE,"\n",file="PIRKC.DAT",append=TRUE)
cat("#female survey sample sizes","\n",file="PIRKC.DAT",append=TRUE)
cat(nSampFlenfUSE,"\n",file="PIRKC.DAT",append=TRUE)

cat("#total survey M numbers","\n",file="PIRKC.DAT",append=TRUE)
cat(useSurvM,"\n",file="PIRKC.DAT",append=TRUE)
cat("#male survey CV","\n",file="PIRKC.DAT",append=TRUE)
cat(useMSurveyCV,"\n",file="PIRKC.DAT",append=TRUE)
cat("#total survey F numbers","\n",file="PIRKC.DAT",append=TRUE)
cat(useSurvF,"\n",file="PIRKC.DAT",append=TRUE)
cat("#female survey CV","\n",file="PIRKC.DAT",append=TRUE)
cat(useFSurveyCV,"\n",file="PIRKC.DAT",append=TRUE)

cat("#male survey length frequency","\n",file="PIRKC.DAT",append=TRUE)
for(i in 1:nrow(LenFreqM))
 cat(LenFreqM[i,cutLenBins:ncol(LenFreqM)],"\n",file="PIRKC.DAT",append=TRUE)
cat("#female survey length frequency","\n",file="PIRKC.DAT",append=TRUE)
for(i in 1:nrow(LenFreqM))
 cat(LenFreqF[i,cutLenBins:ncol(LenFreqM)],"\n",file="PIRKC.DAT",append=TRUE)

cat("#LengthBins","\n",file="PIRKC.DAT",append=TRUE)
cat(tempMhist$mids[cutLenBins:length(tempMhist$mids)],"\n",file="PIRKC.DAT",append=TRUE)

cat("#Maturity at length females","\n",file="PIRKC.DAT",append=TRUE)
PropMatF<-1/(1+(1.416*10^13)*exp(tempMhist$mids[cutLenBins:length(tempMhist$mids)]*-0.297))
cat(PropMatF,"\n",file="PIRKC.DAT",append=TRUE)
cat("#Maturity at length males","\n",file="PIRKC.DAT",append=TRUE)
PropMatM<-1/(1+(5.842*10^14)*exp(tempMhist$mids[cutLenBins:length(tempMhist$mids)]*-.288))
cat(PropMatM,"\n",file="PIRKC.DAT",append=TRUE)

}
