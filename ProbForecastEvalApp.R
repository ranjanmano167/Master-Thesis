source("C:/Users/ci1mmu/QUICS/7_Reference/7_Misc/Thesis/4_Performance Measures/R_Scripts/Gen_util.R")
source("C:/Users/ci1mmu/QUICS/7_Reference/7_Misc/Thesis/5_Error Model/Scripts/ErrorModelv3.R")
source("C:/Users/ci1mmu/QUICS/7_Reference/7_Misc/Thesis/5_Error Model/Scripts/ErrorModelOri.R")
source("C:/Users/ci1mmu/QUICS/7_Reference/7_Misc/Thesis/5_Error Model/Scripts/ProbForecastEval.R")
directoryCal="C:/Users/ci1mmu/QUICS/7_Reference/7_Misc/Thesis/4_Performance Measures/Data/Calibration/"
directoryVal="C:/Users/ci1mmu/QUICS/7_Reference/7_Misc/Thesis/4_Performance Measures/Data/Validation/"
directoryObs="C:/Users/ci1mmu/QUICS/7_Reference/7_Misc/Thesis/4_Performance Measures/Data/FDC/"
stationFileDir="C:/Users/ci1mmu/QUICS/7_Reference/7_Misc/Thesis/4_Performance Measures/Station_Data.csv"
Sys.setenv(TZ='GMT')

stationsData = read.table(stationFileDir,header=FALSE, sep=",")
stationsID=as.character(stationsData[[1]])

stList=c(1:5)
ltList=c(1,3,6,12,24)

#Comparison of preliminary verification
taus=c(0.05,0.25, 0.5, 0.75,0.95)

#Case 1
NameCaseOne="QR-ORI"
EMTypeCaseOne=2  #1-NQT 2-Original
SimWeightsCaseOne=0 #0-No weigthts, 1-Rank based weights, 2-Value based weights 
extrapMethodCaseOne=1 #1-linear regression, 2-log regression

#Case 2
NameCaseTwo="QR-NQT"
EMTypeCaseTwo=1  #1-NQT 2-Original
SimWeightsCaseTwo=0 #0-No weigthts, 1-Rank based weights, 2-Value based weights 
extrapMethodCaseTwo=1 #1-linear regression, 2-log regression

#Case 3
NameCaseThree="QR-WT"
EMTypeCaseThree=2  #1-NQT 2-Original
SimWeightsCaseThree=1 #0-No weigthts, 1-Rank based weights, 2-Value based weights 
extrapMethodCaseThree=1 #1-linear regression, 2-log regression




meanDiffORIAll=NULL
meanDiff2All=NULL
meanDiff3All=NULL

for (j in stList) {
  
  #Calibration
  Data=preprocdata(directoryCal,stationsID[j],na.rm=TRUE)
  obs=Data[[1]]
  sim=Data[[2]]
  
  obsZoo=converttoZoo(obs)
  simZoo=converttoZoo(sim)
  
  
  obsData=obs[,2:25]
  simData=sim[,2:25]
  err=obsData-simData
  
  errorModelORI=quantileCalibOri(obs=obsData, sim=simData,taus=taus, weights=SimWeightsCaseOne)
  
  if (EMTypeCaseTwo==1) {
    errorModel2=quantileCalib(obsData, simData,taus=taus, NQTweights=SimWeightsCaseTwo)
  } else if (EMTypeCaseTwo==2){
    errorModel2=quantileCalibOri(obs=obsData, sim=simData,taus=taus, weights=SimWeightsCaseTwo)
  }
  
  if (EMTypeCaseThree==1) {
    errorModel3=quantileCalib(obsData, simData,taus=taus, NQTweights=SimWeightsCaseThree)
  } else if (EMTypeCaseThree==2){
    errorModel3=quantileCalibOri(obs=obsData, sim=simData,taus=taus, weights=SimWeightsCaseThree)
  }
  
  #forecast
  dataFore=preprocdata(directoryVal,stationsID[j],startDate="2009-12-22 01:00",na.rm=TRUE)
  obsFore=dataFore[[1]]
  simFore=dataFore[[2]]
  
  obsForeData=obsFore[,2:25]
  simForeData=simFore[,2:25]
  
  PreliVerifyORI=PICPPredict(obsData=obsForeData,upThresQuantile=0,simData=simForeData,errorModel=errorModelORI,errModType=EMTypeCaseOne,extrapMethod=extrapMethodCaseOne)
  PreliVerify2=PICPPredict(obsData=obsForeData,upThresQuantile=0,simData=simForeData,errorModel=errorModel2,errModType=EMTypeCaseTwo,extrapMethod=extrapMethodCaseTwo)
  PreliVerify3=PICPPredict(obsData=obsForeData,upThresQuantile=0,simData=simForeData,errorModel=errorModel3,errModType=EMTypeCaseThree,extrapMethod=extrapMethodCaseThree)
  
  
  absDiffFifetyORI=abs(PreliVerifyORI[[1]][,1]-50)
  absDiffNinetyORI=abs(PreliVerifyORI[[1]][,2]-90)
  meanDiffFifetyORI=mean(absDiffFifetyORI)
  meanDiffNinetyORI=mean(absDiffNinetyORI)
  meanDiffORI=cbind(meanDiffFifetyORI,meanDiffNinetyORI)
  meanDiffORIAll=rbind(meanDiffORIAll,meanDiffORI)
  
  absDiffFifety2=abs(PreliVerify2[[1]][,1]-50)
  absDiffNinety2=abs(PreliVerify2[[1]][,2]-90)
  meanDiffFifety2=mean(absDiffFifety2)
  meanDiffNinety2=mean(absDiffNinety2)
  meanDiff2=cbind(meanDiffFifety2,meanDiffNinety2)
  meanDiff2All=rbind(meanDiff2All,meanDiff2)
  
  absDiffFifety3=abs(PreliVerify3[[1]][,1]-50)
  absDiffNinety3=abs(PreliVerify3[[1]][,2]-90)
  meanDiffFifety3=mean(absDiffFifety3)
  meanDiffNinety3=mean(absDiffNinety3)
  meanDiff3=cbind(meanDiffFifety3,meanDiffNinety3)
  meanDiff3All=rbind(meanDiff3All,meanDiff3)
  
  #colList1=c("lightblue1","lightblue","lightblue3","royalblue","mediumblue","navy")
  colList1=c("lightblue1","lightblue","steelblue3","royalblue","navy")
  #colList2=c("seagreen1","seagreen3","seagreen4","darkgreen")
  #colList3=c("hotpink","red1","red3","red4")
  
  #colListOri=colList#rainbow(6)#c("grey75","grey60","grey45","grey30","grey15","grey0")
  
  png(paste("C:/Users/ci1mmu/QUICS/7_Reference/7_Misc/Thesis/5_Error Model/Verification/Preliminary/PICP&MPI_",j, ".png"),width=480*1.25,height=480*1.25)
  plot(PreliVerifyORI[[1]][ltList,1],PreliVerifyORI[[2]][ltList,1],
       xlim=c(0,100),ylim=c(min(PreliVerifyORI[[2]],PreliVerify2[[2]],PreliVerify3[[2]]),max(PreliVerifyORI[[2]],PreliVerify2[[2]],PreliVerify3[[2]])),
       main=paste("PICP-MPI graph , St:",j ),xlab="PICP [%]", ylab="MPI [m3/s]",xaxt = "n",
       col=colList1,pch=19,cex=2)
  axis(1,seq(0,100,10))
  grid(nx=FALSE,ny=NULL)
  points(PreliVerifyORI[[1]][ltList,2],PreliVerifyORI[[2]][ltList,2],
         xlim=c(0,100),col=colList1,pch=19,cex=2)
  
  #case2
  points(PreliVerify2[[1]][ltList,1],PreliVerify2[[2]][ltList,1],
         xlim=c(0,100),col=colList1,pch=15,cex=2)
  points(PreliVerify2[[1]][ltList,2],PreliVerify2[[2]][ltList,2],
         xlim=c(0,100),col=colList1,pch=15,cex=2)
  
  #case3
  points(PreliVerify3[[1]][ltList,1],PreliVerify3[[2]][ltList,1],
         xlim=c(0,100),col=colList1,pch=17,cex=2)
  points(PreliVerify3[[1]][ltList,2],PreliVerify3[[2]][ltList,2],
         xlim=c(0,100),col=colList1,pch=17,cex=2)
  
  abline(v=c(50,90),lty=c(2,2))
  #legend("topleft",legend=rep(paste0("LT:",ltList,"h"),3),col=rep(colList1,3),
  #       pch=rep(c(19,15,17),1,each=5),pt.cex=2,ncol=3,title=paste(NameCaseOne,"-",NameCaseTwo,"-",NameCaseThree),title.adj=0.2)
  dev.off()
}


#Probability forecast verification

########CRPSS#########
CRPSSStations=NULL
for (j in stList) {
  
  #Calibration
  Data=preprocdata(directoryCal,stationsID[j],na.rm=TRUE)
  obs=Data[[1]]
  sim=Data[[2]]
  
  obsZoo=converttoZoo(obs)
  simZoo=converttoZoo(sim)
  
  obsData=obs[,2:25]
  simData=sim[,2:25]
  err=obsData-simData
  
  errorModel=errModelRPS(obsData, simData, taus=c(0.025,0.5,0.975),errModType=2, type=1, extrap=1, weights=0)
  
  #forecast
  dataFore=preprocdata(directoryVal,stationsID[j],startDate="2009-12-22 01:00",na.rm=TRUE)
  obsFore=dataFore[[1]]
  simFore=dataFore[[2]]
  
  obsNew=obsFore[,2:25]
  simNew=simFore[,2:25]
  
  simClimChar=getQuantile(directory=directoryObs,stationsID=stationsID,st=j, quantile=c(0.14,0.5,0.86))
  simClimMean=simClimChar[[2]]
  simClimSD=((simClimChar[[2]]-simClimChar[[1]])+(simClimChar[[3]]-simClimChar[[2]]))/2
  
  m=NROW(simNew)
  quantile=c(0,0.05,seq(0.1,0.9,0.1),0.95)
  
  CRPSSLeadTime=NULL
  for (i in ltList) {
    CRPSSQuantile=NULL
    for (k in 1:length(quantile)) {
      CRPSInput=CRPSData(obsNew[,i],simNew[,i],errorModel,errModType=2, type=1,extrap=1,quantile=quantile[k],leadTime=i)
      crps=crps(CRPSInput[,1],CRPSInput[,2:3])
      crpsRef=crps(CRPSInput[,1],c(simClimMean,simClimSD))
      crpss=1-crps[[2]]/crpsRef[[2]]
      CRPSSQuantile=cbind(CRPSSQuantile,crpss)
    }
    CRPSSLeadTime=rbind(CRPSSLeadTime,CRPSSQuantile)
  }
  CRPSSStations=round(abind(CRPSSStations,CRPSSLeadTime,along=3),3)
}

CRPSS_NQT=CRPSSStations
CRPSS_ORI=CRPSSStations
CRPSS_NQTW=CRPSSStations

#Plots

for (j in stList) {
  
  
  png(paste("C:/Users/ci1mmu/QUICS/7_Reference/7_Misc/Thesis/5_Error Model/Verification/Extensive/CRPSS_",j, ".png"),width=720,height=480)
  par(mfrow=c(2,2))
  for (i in 1:NROW(CRPSS_NQT)) {
    plot(c(0.05,seq(0.1,0.9,0.1),0.95), CRPSS_NQT[i,,j],xlim=c(0,1),ylim=c(0.6,1),lwd=2,
         type="b", pch=1,xaxt = "n", main=paste("Continuous Ranked Probability Skill Score (CRPSS), St:",j,", LT (hrs):",ltList[i]),
         xlab="prob. of non-exceedance [-]",ylab="RPSS [-]",col="blue") 
    axis(1,seq(0,1,0.1))
    grid(nx=0.1,ny=NULL)
    points(c(0.05,seq(0.1,0.9,0.1),0.95),CRPSS_ORI[i,,j],type="b",pch=0,lwd=2,col="red")
    points(c(0.05,seq(0.1,0.9,0.1),0.95),CRPSS_NQTW[i,,j],type="b",pch=3,lwd=2,col="limegreen")
    #points(c(0.05,seq(0.1,0.9,0.1),0.95),CRPSSLog[i,,j],type="b",pch=8,lwd=2,col="seagreen4")
    
  }
  legend("topleft",legend=c("QR-NQT","QR-Original","QR-Weights"),pch=c(1,0,3),lty=c(1,1,1),
         col=c("blue","red","limegreen"))
  dev.off()
}

#finalplot

for (i in 1:5) {
  png(paste("C:/Users/ci1mmu/QUICS/7_Reference/7_Misc/Thesis/5_Error Model/Verification/Extensive/CRPSS_",i, ".png"),width=720,height=480)
  plot(NULL, NULL,xlim=c(0,1),ylim=c(0,1),lwd=2,
       type="b", pch=1,xaxt = "n", main=paste("Continuous Ranking Probability Skill Score (CRPSS), LT (hrs):",ltList[i]),
       xlab="Prob. of non-exceedance [-]",ylab="CRPSS [-]",cex.axis=1.5,cex.lab = 1.5,cex.main=1.5) 
  axis(1,seq(0,1,0.1),cex.axis=1.5,cex.lab = 1.5)
  abline(h=seq(0.1,0.9,0.1), lty=3,col="gray")
  abline(v=seq(0,1,0.1), lty=3,col="gray")
  
  colList=c("purple","red","darkgreen","royalblue","limegreen",
            "orange","violetred","turquoise","tomato","dimgrey") 
  for(j in stList) {
    points(c(0,0.05,seq(0.1,0.9,0.1),0.95),CRPSS_ORI[i,,j],type="b",pch=1,lwd=2,col=colList[j])
  }
  legend("bottomright",legend=c("st-1","st-2","st-3","st-4","st-5"),pch=c(1),lty=1,lwd=2,cex=1.5,
         col=colList,bty="n")
  dev.off()
}





##########RPSS#############

RPSSStations=NULL

for (j in stList) {
  
  #Calibration
  Data=preprocdata(directoryCal,stationsID[j],na.rm=TRUE)
  obs=Data[[1]]
  sim=Data[[2]]
  
  obsZoo=converttoZoo(obs)
  simZoo=converttoZoo(sim)
  
  
  obsData=obs[,2:25]
  simData=sim[,2:25]
  err=obsData-simData
  
  errorModel=errModelRPS(obsData, simData, taus=seq(0.04,0.96,0.04),errModType=2, type=1, extrap=1, weights=0)
  
  #forecast
  dataFore=preprocdata(directoryVal,stationsID[j],startDate="2009-12-22 01:00",na.rm=TRUE)
  obsFore=dataFore[[1]]
  simFore=dataFore[[2]]
  
  obsNew=obsFore[,2:25]
  simNew=simFore[,2:25]
  
  simBins=getQuantile(directory=directoryObs,stationsID=stationsID,st=j, quantile=c(seq(0.1,0.9,0.1),0.95))
  
  quantile=c(0,0.05,seq(0.1,0.9,0.1),0.95)
  
  RPSSLeadTime=NULL
  for (i in ltList) {
    RPSSQuantile=NULL
    for (k in 1:length(quantile)) {
      RPSInput=RPSData(obsNew[,i],simNew[,i],simBins=simBins, errorModel,errModType=2, type=1,extrap=1,quantile=quantile[k],leadTime=i)
      rps=rps(obs=RPSInput[[1]],pred=t(RPSInput[[2]]))
      rpss=rps[[2]]
      RPSSQuantile=cbind(RPSSQuantile,rpss)
    }
    RPSSLeadTime=rbind(RPSSLeadTime,RPSSQuantile)
  }
  RPSSStations=round(abind(RPSSStations,RPSSLeadTime,along=3),3)
}

RPSS_NQT=RPSSStations
RPSS_ORI=RPSSStations
RPSS_NQTW=RPSSStations

#Plot:RPSS

for (j in 1:5) {
  
  
  png(paste("C:/Users/ci1mmu/QUICS/7_Reference/7_Misc/Thesis/5_Error Model/Verification/Extensive/RPSS_",j, ".png"),width=720*1.5,height=480*1.5)
  par(mfrow=c(3,2))
  for (i in 1:NROW(RPSSBasic)) {
    plot(c(0.05,seq(0.1,0.9,0.1),0.95), RPSS_NQT[i,,j],xlim=c(0,1),ylim=c(0,1),lwd=2,
         type="b", pch=1,xaxt = "n", main=paste("Ranked Probability Skill Score (RPSS), St:",j,", LT (hrs):",ltList[i]),
         xlab="prob. of non-exceedance [-]",ylab="RPSS [-]",col="blue") 
    axis(1,seq(0,1,0.1))
    grid(nx=0.1,ny=NULL)
    points(c(0.05,seq(0.1,0.9,0.1),0.95),RPSS_ORI[i,,j],type="b",pch=0,lwd=2,col="red")
    points(c(0.05,seq(0.1,0.9,0.1),0.95),RPSS_NQTW[i,,j],type="b",pch=3,lwd=2,col="limegreen")
    #points(c(0.05,seq(0.1,0.9,0.1),0.95),RPSSLog[i,,j],type="b",pch=8,lwd=2,col="limegreen")
    
  }
  legend("topleft",legend=c("QR-NQT","QR-Original","QR-Weights"),pch=c(1,0,3),lty=c(1,1,1),
         col=c("blue","red","limegreen"))
  dev.off()
}

#finalplot

for (i in 1:5) {
  png(paste("C:/Users/ci1mmu/QUICS/7_Reference/7_Misc/Thesis/5_Error Model/Verification/Extensive/RPSS_",i, ".png"),width=720,height=480)
  plot(NULL, NULL,xlim=c(0,1),ylim=c(0,1),lwd=2,
       type="b", pch=1,xaxt = "n", main=paste("Ranking Probability Skill Score (RPSS), LT (hrs):",ltList[i]),
       xlab="Prob. of non-exceedance [-]",ylab="RPSS [-]",cex.axis=1.5,cex.lab = 1.5,cex.main=1.5) 
  axis(1,seq(0,1,0.1),cex.axis=1.5,cex.lab = 1.5)
  abline(h=seq(0.1,0.9,0.1), lty=3,col="gray")
  abline(v=seq(0,1,0.1), lty=3,col="gray")
  
  colList=c("purple","red","darkgreen","royalblue","limegreen",
            "orange","violetred","turquoise","tomato","dimgrey") 
  for(j in stList) {
    points(c(0,0.05,seq(0.1,0.9,0.1),0.95),RPSS_ORI[i,,j],type="b",pch=1,lwd=2,col=colList[j])
  }
  legend("bottomright",legend=c("st-1","st-2","st-3","st-4","st-5"),pch=c(1),lty=1,lwd=2,cex=1.5,
         col=colList,bty="n")
  dev.off()
}

######Brier skillscore######

brierSSStations=NULL

for (j in stList) {
  
  #Calibration
  Data=preprocdata(directoryCal,stationsID[j],na.rm=TRUE)
  obs=Data[[1]]
  sim=Data[[2]]
  
  obsZoo=converttoZoo(obs)
  simZoo=converttoZoo(sim)
  
  obsData=obs[,2:25]
  simData=sim[,2:25]
  err=obsData-simData
  
  errorModel=errModelRPS(obsData, simData, taus=seq(0.05,0.95,0.05),errModType=1, type=1, extrap=1, weights=0)
  #   errorModelNQTW=errModelRPS(obsData, simData, taus=seq(0.1,0.9,0.1),errModType=1, type=1, extrap=1, weights=1)
  
  #forecast
  dataFore=preprocdata(directoryVal,stationsID[j],startDate="2009-12-22 01:00",na.rm=TRUE)
  obsFore=dataFore[[1]]
  simFore=dataFore[[2]]
  
  obsNew=obsFore[,2:25]
  simNew=simFore[,2:25]
  
  quantiles=getQuantile(directory=directoryObs,stationsID=stationsID,st=j, quantile=c(0.05,seq(0.1,0.9,0.1),0.95))
  #png(paste("C:/Users/ci1mmu/QUICS/7_Reference/7_Misc/Thesis/5_Error Model/Verification/Extensive/Reliability/Rel_",j, ".png"),width=480*1.5,height=480*1.5)
  #par(mfrow=c(2,1))
  
  brierSSLT=NULL
  for (i in ltList) {
    obslt=obsNew[,i]
    simlt=simNew[,i]
    
    brierSSQuantile=NULL
    for (k in 1:length(quantiles)){
      brierData=brierSSData(obs=obslt,sim=simlt,eventthres=quantiles[k],errorModel,errModType=1,type=1, extrap=1,leadTime=i)
      
      brierSS=brier(obs=brierData[,2], pred=brierData[,1])[[4]]
      
      brierSSQuantile=cbind(brierSSQuantile,brierSS)
    }
    brierSSLT=rbind(brierSSLT,brierSSQuantile)
    
}
brierSSStations=round(abind(brierSSStations,brierSSLT,along=3),3)
}

brier_NQT=brierSSStations
brier_ORI=brierSSStations
brier_NQTW=brierSSStations

#BSS Plots

#comparison
for (j in stList) {
  
  
  png(paste("C:/Users/ci1mmu/QUICS/7_Reference/7_Misc/Thesis/5_Error Model/Verification/Extensive/brierSS_",j, ".png"),width=720,height=480)
  par(mfrow=c(3,2))
  for (i in 1:NROW(brier_NQT)) {
    plot(c(0.05,seq(0.1,0.9,0.1),0.95), brier_NQT[i,,j],xlim=c(0,1),ylim=c(0,1),lwd=2,
         type="b", pch=1,xaxt = "n", main=paste("Brier Skill Score (BSS), St:",j,", LT (hrs):",ltList[i]),
         xlab="prob. of non-exceedance [-]",ylab="BSS [-]",col="blue") 
    axis(1,seq(0,1,0.1))
    grid(nx=0.1,ny=NULL)
    points(c(0.05,seq(0.1,0.9,0.1),0.95),brier_ORI[i,,j],type="b",pch=0,lwd=2,col="red")
    points(c(0.05,seq(0.1,0.9,0.1),0.95),brier_NQTW[i,,j],type="b",pch=3,lwd=2,col="limegreen")
    #points(c(0.05,seq(0.1,0.9,0.1),0.95),CRPSSLog[i,,j],type="b",pch=8,lwd=2,col="seagreen4")
    
  }
  legend("topleft",legend=c("QR-NQT","QR-ORI","QR-NQTW"),pch=c(1,0,3),lty=c(1,1,1),
         col=c("blue","red","limegreen"))
  dev.off()
}

#finalplot

for (i in 1:4) {
  png(paste("C:/Users/ci1mmu/QUICS/7_Reference/7_Misc/Thesis/5_Error Model/Verification/Extensive/brierSS_",i, ".png"),width=720,height=480)
  plot(NULL, NULL,xlim=c(0,1),ylim=c(-1.5,1),lwd=2,
       type="b", pch=1,xaxt = "n", main=paste("Brier Skill Score (BSS), LT (hrs):",ltList[i]),
       xlab="Threshold Quantiles [-]",ylab="BSS [-]",cex.axis=1.5,cex.lab = 1.5,cex.main=1.5) 
  axis(1,seq(0,1,0.1),cex.axis=1.5,cex.lab = 1.5)
  grid(nx=0.1,ny=NULL)
  abline(a=0,b=0)
  abline(h=seq(0.1,0.9,0.1), lty=3,col="gray")
  abline(v=seq(0,1,0.1), lty=3,col="gray")
  
  colList=c("purple","red","darkgreen","royalblue","limegreen",
            "orange","violetred","turquoise","tomato","dimgrey") 
  for(j in stList) {
    points(c(0.05,seq(0.1,0.9,0.1),0.95),brier_ORI[i,,j],type="b",pch=1,lwd=2,col=colList[j])
  }
  legend("bottomright",legend=c("st-1","st-2","st-3","st-4","st-5"),pch=c(1),lty=1,lwd=2,cex=1.5,
         col=colList,bty="n")
  dev.off()
  }




#reliability diagram
for (j in stList) {
  
  #Calibration
  Data=preprocdata(directoryCal,stationsID[j],na.rm=TRUE)
  obs=Data[[1]]
  sim=Data[[2]]
  
  obsZoo=converttoZoo(obs)
  simZoo=converttoZoo(sim)
  
  obsData=obs[,2:25]
  simData=sim[,2:25]
  err=obsData-simData
  
  errorModel=errModelRPS(obsData, simData, taus=seq(0.04,0.96,0.04),errModType=2, type=1, extrap=1, weights=1)
  #   errorModelNQTW=errModelRPS(obsData, simData, taus=seq(0.1,0.9,0.1),errModType=1, type=1, extrap=1, weights=1)
  
  #forecast
  dataFore=preprocdata(directoryVal,stationsID[j],startDate="2009-12-22 01:00",na.rm=TRUE)
  obsFore=dataFore[[1]]
  simFore=dataFore[[2]]
  
  obsNew=obsFore[,2:25]
  simNew=simFore[,2:25]
  
  quantiles=getQuantile(directory=directoryObs,stationsID=stationsID,st=j, quantile=seq(0.1,0.9,0.2))
  
  for (i in ltList) {
    
    obslt=obsNew[,i]
    simlt=simNew[,i]
    
    png(paste("C:/Users/ci1mmu/QUICS/7_Reference/7_Misc/Thesis/5_Error Model/Verification/Extensive/Reliability/Rel_",j,i, ".png"),width=480*1.5,height=480*2)
    par(mfrow=c(2,1))
    
    par(mar=c(0.2,4.1,4.1,2.1))
    plot(NULL,NULL,xlim=c(0,1),ylim=c(0,1),,cex.axis=1.5,cex.lab = 1.5,cex.main=1.5,
         main=paste("Reliability Diagram, St:",j,", LT (hrs):",i),
         xlab="Forecast probability",ylab="Observed Frequency",xaxt="n")
    
    grid(nx=NULL,ny=NULL)    
    colours=c("purple","red","darkgreen","royalblue","limegreen",
              "orange","violetred","turquoise","tomato","dimgrey") 
    
    for (k in 1:length(quantiles)){
      
      reliabData=relDiaData(obs=obslt,sim=simlt,errModType=2,eventthres=quantiles[k],probBins=seq(0,1,0.1),errorModel,leadTime=i)
      points(reliabData[,1],reliabData[,4],type='b',col=colours[k],pch=20,lwd=2)
      
    }
    abline(a=0,b=1,col="black",lwd=2)
    
    
    par(mar=c(4.1,4.1,0.2,2.1))
    plot(NULL,NULL,xlim=c(0,1),ylim=c(0,1500),cex.axis=1.5,cex.lab = 1.5,cex.main=1.5,
         xlab="Forecast probability",ylab="Sample Size",)
    grid(nx=NULL,ny=NA) 
    for (k in 1:length(quantiles)){
      
      reliabData=relDiaData(obs=obslt,sim=simlt,errModType=2,eventthres=quantiles[k],probBins=seq(0,1,0.1),errorModel,leadTime=i) 
      points(reliabData[,1],reliabData[,5],type='b',col=colours[k],pch=20,lwd=2)
    }
    legend("top", legend=seq(0.1,0.9,0.2),title="Threshold quantiles", bty="n",col=colours[1:5],lty=1,lwd=2,xjust=0.5,pch=20,cex=1.5,ncol=5)
    dev.off()
  }
}


#PICP and MPI for Probability of non-exceedence of 0.75
PICPNQTAll75=NULL
PICPNQTWAll75=NULL
PICPORIAll75=NULL
PICPORIWAll75=NULL

MPINQTAll75=NULL
MPINQTWAll75=NULL
MPIORIAll75=NULL
MPIORIWAll75=NULL

for (j in stList) {
  quantiles=c(0.025,seq(0.05,0.45,0.05),0.475,0.525,seq(0.55,0.95,0.05),0.975)
  
  #Calibration
  Data=preprocdata(directoryCal,stationsID[j],na.rm=TRUE)
  obs=Data[[1]]
  sim=Data[[2]]
  
  obsZoo=converttoZoo(obs)
  simZoo=converttoZoo(sim)
  
  obsData=obs[,2:25]
  simData=sim[,2:25]
  err=obsData-simData
  
  errorModelNQT=errModelRPS(obsData, simData, taus=quantiles,errModType=1, type=1, extrap=1, weights=0)
  errorModelNQTW=errModelRPS(obsData, simData, taus=quantiles,errModType=1, type=1, extrap=1, weights=1)
  errorModelORI=errModelRPS(obsData, simData, taus=quantiles,errModType=2, type=1, extrap=1, weights=0)
  errorModelORIW=errModelRPS(obsData, simData, taus=quantiles,errModType=2, type=1, extrap=1, weights=1)
  
  #forecast
  dataFore=preprocdata(directoryVal,stationsID[j],startDate="2009-12-22 01:00",na.rm=TRUE)
  obsFore=dataFore[[1]]
  simFore=dataFore[[2]]
  
  obsNew=obsFore[,2:25]
  simNew=simFore[,2:25]
  
  PICPnMPINQT=PICPPredictAll(obsNew,simNew,upThresQuantile=0.75,confBands=quantiles,errorModel=errorModelNQT,errModType=1,extrapMethod=1) 
  PICPnMPINQTW=PICPPredictAll(obsNew,simNew,upThresQuantile=0.75,confBands=quantiles,errorModel=errorModelNQTW,errModType=1,extrapMethod=1) 
  PICPnMPIORI=PICPPredictAll(obsNew,simNew,upThresQuantile=0.75,confBands=quantiles,errorModel=errorModelORI,errModType=2,extrapMethod=1)
  PICPnMPIORIW=PICPPredictAll(obsNew,simNew,upThresQuantile=0.75,confBands=quantiles,errorModel=errorModelORIW,errModType=2,extrapMethod=1) 
  
  PICPNQT=PICPnMPINQT[[2]]
  PICPNQTW=PICPnMPINQTW[[2]]
  PICPORI=PICPnMPIORI[[2]]
  PICPORIW=PICPnMPIORIW[[2]]
  
  MPINQT=PICPnMPINQT[[3]]
  MPINQTW=PICPnMPINQTW[[3]]
  MPIORI=PICPnMPIORI[[3]]
  MPIORIW=PICPnMPIORIW[[3]]
  
  PICPNQTAll75=abind(PICPNQTAll75,PICPNQT,along=3)
  PICPNQTWAll75=abind(PICPNQTWAll75,PICPNQTW,along=3)
  PICPORIAll75=abind(PICPORIAll75,PICPORI,along=3)
  PICPORIWAll75=abind(PICPORIWAll75,PICPORIW,along=3)
  
  MPINQTAll75=abind(MPINQTAll75,MPINQT,along=3)
  MPINQTWAll75=abind(MPINQTWAll75,MPINQTW,along=3)
  MPIORIAll75=abind(MPIORIAll75,MPIORI,along=3)
  MPIORIWAll75=abind(MPIORIWAll75,MPIORIW,along=3)
}





####PICP plots

#PICP all
stList=c(1,4)
ltList=c(1,3,6,12,24)

for (j in stList) {
  for (i in ltList){
  plot(NULL,NULL,type="l",col="red",xlim=c(0,100),ylim=c(0,100),
       main=paste("PICP",j,i),lwd=2)
  abline(a=0,b=1,lty=2)
  points(c(0.05,seq(0.1,0.9,0.1),0.95)*100,rev(PICPNQTAll[,i,j]),type="l",col="red",lwd=2)
  points(c(0.05,seq(0.1,0.9,0.1),0.95)*100,rev(PICPNQTWAll[,i,j]),type="l",col="seagreen",lwd=2)
  points(c(0.05,seq(0.1,0.9,0.1),0.95)*100,rev(PICPORIAll[,i,j]),type="l",col="blue",lwd=2)
  points(c(0.05,seq(0.1,0.9,0.1),0.95)*100,rev(PICPORIWAll[,i,j]),type="l",col="purple",lwd=2)
  }
}

#PICP mean
for (j in stList) {
  png(paste("C:/Users/ci1mmu/QUICS/7_Reference/7_Misc/Thesis/5_Error Model/Verification/Extensive/Reliability/Rel_",j, ".png"),width=480*1.25,height=480*1.25)
  par(mfrow=c(1,1))
  plot(NULL,NULL,type="l",col="red",xlim=c(0,100),ylim=c(0,100),
       main=paste("PICP, St:",j ),lwd=2, xlab="Confidence interval [%]",
       ylab="PICP [%]")
  abline(a=0,b=1,lty=2)
  grid(nx=NULL,ny=NULL)
  points(c(0.05,seq(0.1,0.9,0.1),0.95)*100,rev(rowMeans(PICPNQTAll75[,,j])),type="b",col="red",lwd=2)
  points(c(0.05,seq(0.1,0.9,0.1),0.95)*100,rev(rowMeans(PICPNQTWAll75[,,j])),type="b",col="orange",lwd=2)
  points(c(0.05,seq(0.1,0.9,0.1),0.95)*100,rev(rowMeans(PICPORIAll75[,,j])),type="b",col="blue",lwd=2)
  points(c(0.05,seq(0.1,0.9,0.1),0.95)*100,rev(rowMeans(PICPORIWAll75[,,j])),type="b",col="seagreen4",lwd=2)
  legend("bottomright",legend=c("QR-ORI","QR-NQT","QR-WT"),col=c("blue","red","seagreen4"),lty=1,pch=1,lwd=2,bty="n")
  dev.off()
}

#MPI mean
for (j in stList) {
  png(paste("C:/Users/ci1mmu/QUICS/7_Reference/7_Misc/Thesis/5_Error Model/Verification/Extensive/Reliability/Rel_",j, ".png"),width=480*1.25,height=480*1.25)
  par(mfrow=c(1,1))
  plot(NULL,NULL,type="l",col="red",xlim=c(0,100),ylim=c(0,100),
       main=paste("PICP, St:",j ),lwd=2, xlab="Confidence interval [%]",
       ylab="PICP [%]")
  #abline(a=0,b=1,lty=2)
  grid(nx=NULL,ny=NULL)
  points(c(0.05,seq(0.1,0.9,0.1),0.95)*100,rev(rowMeans(MPINQTAll75[,,j])),type="b",col="red",lwd=2)
  points(c(0.05,seq(0.1,0.9,0.1),0.95)*100,rev(rowMeans(MPINQTWAll75[,,j])),type="b",col="orange",lwd=2)
  points(c(0.05,seq(0.1,0.9,0.1),0.95)*100,rev(rowMeans(MPIORIAll75[,,j])),type="b",col="blue",lwd=2)
  points(c(0.05,seq(0.1,0.9,0.1),0.95)*100,rev(rowMeans(MPIORIWAll75[,,j])),type="b",col="seagreen4",lwd=2)
  legend("bottomright",legend=c("QR-ORI","QR-NQT","QR-WT"),col=c("blue","red","seagreen4"),lty=1,pch=1,lwd=2,bty="n")
  dev.off()
}



##Bin
# ##both mean MPI and mean PICP in same plot
# for (j in stList) {
#   png(paste("C:/Users/ci1mmu/QUICS/7_Reference/7_Misc/Thesis/5_Error Model/Verification/Extensive/Reliability/Rel_",j, ".png"),width=480*1.25,height=480*1.25)
#   PICPMeanCom=cbind(rev(rowMeans(PICPORIAll75[,,j])),rev(rowMeans(PICPNQTAll75[,,j])),rev(rowMeans(PICPORIWAll75[,,j])))
#   MPIMeanCom=cbind(rev(rowMeans(MPIORIAll75[,,j])),rev(rowMeans(MPINQTAll75[,,j])),rev(rowMeans(MPIORIWAll75[,,j])))
#   conBounds=c(0.05,seq(0.1,0.9,0.1),0.95)*100
#   
#   
#   twoord.stackplot(lx=conBounds,ldata=PICPMeanCom,rx=conBounds,rdata=MPIMeanCom,
#                    ltype="b",rtype="b",lcol=c("blue","red","seagreen4"),rcol=c("blue","red","seagreen4"),
#                    lwd=2,lty=2,
#                    main=paste("Mean MPI and Mean PICP, St:",j ))
#   grid(nx=NULL,ny=NULL)
#   abline(a=0,b=1)
#   dev.off()
# }
# 
# par(mar = c(5,5,2,5))
# plot(ts(x), ylab = "x")
# par(new = T)
# plot(ts(y), col = "red", axes = F, xlab = NA, ylab = NA)
# axis(side = 4)
# mtext(side = 4, line = 3, "y")
# 
# 
# for (j in stList) {
#   png(paste("C:/Users/ci1mmu/QUICS/7_Reference/7_Misc/Thesis/5_Error Model/Verification/Extensive/Reliability/Rel_",j, ".png"),width=480*1.25,height=480*1.25)
#   PICPMeanCom=cbind(rev(rowMeans(PICPORIAll70[,,j])),rev(rowMeans(PICPNQTAll70[,,j])),rev(rowMeans(PICPORIWAll70[,,j])))
#   MPIMeanCom=cbind(rev(rowMeans(MPIORIAll70[,,j])),rev(rowMeans(MPINQTAll70[,,j])),rev(rowMeans(MPIORIWAll70[,,j])))
#   conBounds=c(0.05,seq(0.1,0.9,0.1),0.95)*100
#   
#   
#   par(mar = c(5,5,2,5))
#   plot(conBounds,PICPMeanCom[,1], ylab = "x",yaxt="n")
#   axis(side=2,at=seq(0,100,20),xpd=TRUE)
#   par(new = T)
#   plot(conBounds,MPIMeanCom[,1], col = "red", axes = F, xlab = NA, ylab = NA,ylim=c(0,100))
#   axis(side = 4)
#   mtext(side = 4, line = 3, "y")
#   
#   dev.off()
# }
# 
# 
# Reliability diagrams
# for (j in stList) {
#   
#   #Calibration
#   Data=preprocdata(directoryCal,stationsID[j],na.rm=TRUE)
#   obs=Data[[1]]
#   sim=Data[[2]]
#   
#   obsZoo=converttoZoo(obs)
#   simZoo=converttoZoo(sim)
#   
#   obsData=obs[,2:25]
#   simData=sim[,2:25]
#   err=obsData-simData
#   
#   errorModelNQT=errModelRPS(obsData, simData, taus=seq(0.1,0.9,0.1),errModType=1, type=1, extrap=1, weights=0)
#   errorModelNQTW=errModelRPS(obsData, simData, taus=seq(0.1,0.9,0.1),errModType=1, type=1, extrap=1, weights=1)
#   
#   #forecast
#   dataFore=preprocdata(directoryVal,stationsID[j],startDate="2009-12-22 01:00",na.rm=TRUE)
#   obsFore=dataFore[[1]]
#   simFore=dataFore[[2]]
#   
#   obsNew=obsFore[,2:25]
#   simNew=simFore[,2:25]
#   
#   png(paste("C:/Users/ci1mmu/QUICS/7_Reference/7_Misc/Thesis/5_Error Model/Verification/Extensive/Reliability/Rel_",j, ".png"),width=480*1.5,height=480*1.5)
#   par(mfrow=c(2,2))
#   for (i in ltList) {
#     dataRelPlotNQT=dataReliabPlot(obsNew,simNew,errorModelNQT,errModType=1,extrapMethod=1,leadTime=i)
#     dataRelPlotNQTW=dataReliabPlot(obsNew,simNew,errorModelNQTW,errModType=1,extrapMethod=1,leadTime=i)
#     plot(c(dataRelPlotNQT[[1]],1),dataRelPlotNQT[[2]],type='b',xlim=c(0,1),ylim=c(0,1),
#          xlab="forecast probability")
#     points(c(dataRelPlotNQTW[[1]],1),dataRelPlotNQTW[[2]],type='b',col="green")
#     abline(a=0,b=1,col="blue")
#   }
#   dev.off()
# }

