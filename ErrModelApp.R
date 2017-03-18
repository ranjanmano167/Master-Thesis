source("C:/Users/ci1mmu/QUICS/7_Reference/7_Misc/Thesis/4_Performance Measures/R_Scripts/Gen_util.R")
source("C:/Users/ci1mmu/QUICS/7_Reference/7_Misc/Thesis/5_Error Model/Scripts/ErrorModelv3.R")
source("C:/Users/ci1mmu/QUICS/7_Reference/7_Misc/Thesis/5_Error Model/Scripts/ErrorModelOri.R")
directoryCal="C:/Users/ci1mmu/QUICS/7_Reference/7_Misc/Thesis/4_Performance Measures/Data/Calibration/"
directoryVal="C:/Users/ci1mmu/QUICS/7_Reference/7_Misc/Thesis/4_Performance Measures/Data/Validation/"
stationFileDir="C:/Users/ci1mmu/QUICS/7_Reference/7_Misc/Thesis/4_Performance Measures/Station_Data.csv"
Sys.setenv(TZ='GMT')

stationsData = read.table(stationFileDir,header=FALSE, sep=",")
stationsID=as.character(stationsData[[1]])

taus=c(0,0.25,0.5,0.75,0.95) 
#taus=seq(0.04,0.96,0.04)
stList=c(4) #station list (1-Q3420,2-Q3530,3-Q3570,4-Q3725,5-Q3900)
#ltList=6 #lead time list
ltList=c(1,6,12,24)
eventDate="2010-09-19 01:00:00" # Which event to be considered for CI vs leadtime plot
simList=c(seq(0,2300,100),3000,4000,5000) # for Err Vs Sim original domain


#NQT plots
for (j in stList) {
  
  Data=preprocdata(directoryCal,stationsID[j],na.rm=TRUE)
  obs=Data[[1]]
  sim=Data[[2]]
  
  obsZoo=converttoZoo(obs)
  simZoo=converttoZoo(sim)
  
  obsData=obs[,2:25]
  simData=sim[,2:25]
  err=obsData-simData
  
  errorModel=quantileCalib(obsData, simData,taus=taus, NQTweights=0)
  
  
  png(paste("C:/Users/ci1mmu/QUICS/7_Reference/7_Misc/Thesis/5_Error Model/NQTPlots/NQT_",j, ".png"),width=480*2,height=480*2)
  par(mfrow=c(2,2))
  
  for (i in ltList) { 
    plotNorSimVsErr(errormodel,st=j,lt=i)
    
  }
  legend("bottomleft",bty="n", 
         legend=c("Median","50% CI","90% CI"),
         col ="blue", lty=c(1,2,3),lw=2)
  
  dev.off()
}

#Err Vs Sim original domain,
#One case
for (j in stList) {
  
  Data=preprocdata(directoryCal,stationsID[j],na.rm=TRUE)
  obs=Data[[1]]
  sim=Data[[2]]
  
  obsZoo=converttoZoo(obs)
  simZoo=converttoZoo(sim)
  
  obsData=obs[,2:25]
  simData=sim[,2:25]
  err=obsData-simData
  
  
  errorModel=quantileCalib(obsData, simData,taus=taus)
  png(paste("C:/Users/ci1mmu/QUICS/7_Reference/7_Misc/Thesis/5_Error Model/SimVsErrPlots/Ris_",j, ".png"),width=720,height=480)
  par(mfrow=c(1,1))
  for (i in ltList) {
    ylimit=max(abs(min(err)),abs(max(err)))
    plot(simData[,i],err[,i],ylim=c(-ylimit,ylimit),xlim=c(0,max(simData)),
         col='gray60', main=paste("Residuals, St:",j,", LT (hrs):",i),
         xlab="Predicted discharge [m3/s]",ylab="Residuals [m3/s]")
    plotErrRisOri(errorModel,simList,lt=i)
  }
  legend("bottomleft",bty="n", 
         legend=c("Median","50% CI","90% CI"),
         col ="blue", lty=c(1,2,3),cex=1.5,lw=2)
  dev.off()
}

#Err Vs Sim original domain,
#multi-case
for (j in stList) {
  
  Data=preprocdata(directoryCal,stationsID[j],na.rm=TRUE)
  obs=Data[[1]]
  sim=Data[[2]]
  
  obsZoo=converttoZoo(obs)
  simZoo=converttoZoo(sim)
  
  obsData=obs[,2:25]
  simData=sim[,2:25]
  err=obsData-simData
  
  errorModelORI=quantileCalibOri(obs=obsData, sim=simData,taus=taus,weights=0)
  errorModelNQT=quantileCalib(obsData, simData,taus=taus)
  errorModelNQTW=quantileCalib(obsData, simData,taus=taus,NQTweights=1)
  
  png(paste("C:/Users/ci1mmu/QUICS/7_Reference/7_Misc/Thesis/5_Error Model/SimVsErrPlots/Ris_",j, ".png"),width=720*1.5,height=480*1.5)
  par(mfrow=c(2,2))
  for (i in ltList) {
   
    ylimit=max(abs(min(err[,ltList])),abs(max(err[,ltList])))
    
    plot(simData[,i],err[,i],
         main=(paste("Quantile Regression, St:",j,", LT (hrs):",i)),
         xlim=c(0,max(simData)),ylim=c(-ylimit,ylimit),
         #xlim=c(0,20),ylim=c(-50,50),
         xlab="Predicted discharge [m3/s]",ylab="Residuals [m3/s]",
         col="gray60")
    plotSimVsErrOri(errorModel=errorModelORI,sim=simData,err=err,st=j,lt=i,colour="blue")
    plotErrRisOri(errorModel=errorModelNQT,simList,lt=i,colour="red")
    plotErrRisOri(errorModel=errorModelNQTW,simList,lt=i,colour="seagreen4")
    if (i==1) {
      legend("topleft",bty="n", 
             legend=c("QR-ORI:Median","QR-ORI:50% CI",
                      "QR-NQT:Median","QR-NQT:50% CI",
                      "QR-NQTW:Median","QR-NQTW:50% CI"),
#              "QR-NQT:Median","QR-NQT:50% CI","QR-NQT:90% CI",
#              "QR-WT:Median","QR-WT:50% CI","QR-WT:90% CI"),
             col =c(rep("blue",2),rep("red",2),rep("seagreen4",2)), lty=rep(1:2,3),cex=1,lw=2,ncol=3)
    }
  }
 
  dev.off()
}


#Confidence interval vs events - simgle case

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
  
  errorModel=quantileCalib(obsData, simData,taus=taus)
  
  #forecast
  dataCheck=preprocdata(directoryVal,stationsID[j],startDate="2010-09-01 01:00",na.rm=TRUE)
  obsCheck=dataCheck[[1]]
  simCheck=dataCheck[[2]]
  
  obsCheckData=obsCheck
  simCheckData=simCheck
  
  png(paste("C:/Users/ci1mmu/QUICS/7_Reference/7_Misc/Thesis/5_Error Model/FloodEventsPlots/FE_",j, ".png"),width=720*1.5,height=480*1.5)
  par(mfrow=c(2,2))
  for (i in ltList){
    plotCIlt(errorModel,simAll=simCheckData,obsAll=obsCheckData, st=j,lt=i)
  }
  legend("topleft",bty="n", 
         legend=c("Deterministic forecast","50% CI","90% CI","Observation"),
         col = c("navy", "green3", "green1","palevioletred4"),
         #col = c("navy", "steelblue3", "lightblue1","palevioletred4"),
         pch=c(NA,22,22,16),lty=c(1,NA,NA,NA),
         #pt.bg=c(NA, "steelblue3", "lightblue1",NA),
         pt.bg=c(NA, "green3", "green1",NA),
         pt.cex=c(1,3,3,1))
  dev.off()
}

#Confidence interval vs events - Multi case

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
  
  errorModelORI=quantileCalibOri(obs=obsData,  sim=simData,taus=taus,weights=0)
  errorModelNQT=quantileCalib(obsData, simData,taus=taus,NQTweights=0)
  errorModelNQTW=quantileCalib(obsData, simData,taus=taus,NQTweights=1)
  
  #forecast
  dataCheck=preprocdata(directoryVal,stationsID[j],startDate="2010-10-01 01:00",na.rm=TRUE)
  obsCheck=dataCheck[[1]]
  simCheck=dataCheck[[2]]
  
  obsCheckData=obsCheck
  simCheckData=simCheck
  
  png(paste("C:/Users/ci1mmu/QUICS/7_Reference/7_Misc/Thesis/5_Error Model/FloodEventsPlots/FE_",j, ".png"),width=720*1.5,height=480*1.5)
  par(mfrow=c(2,2))
  for (i in ltList){
    plotCIltOri(errorModelORI,simAll=simCheckData,obsAll=obsCheckData, st=j,lt=i)
    plotCIltLines(errorModelNQT,simAll=simCheckData,obsAll=obsCheckData, st=j,lt=i)
    plotCIltLines(errorModelNQTW,simAll=simCheckData,obsAll=obsCheckData, st=j,lt=i)
  }
  legend("topleft",bty="n", 
         legend=c("Deterministic forecast","QR-ORI:90% CI","QR-NQT:90% CI","QR-NQTW:90% CI","Observation"),
         col = c("navy", "lightblue1", "red","seagreen4","palevioletred4"),
         pch=c(NA,22,NA,NA,16),lty=c(1,NA,1,1,NA),pt.bg=c(NA, "lightblue1", NA,NA,NA),
         pt.cex=c(NA,3,NA,NA,1),lwd=c(1,NA,2,2,NA))
  dev.off()
}

#Confidence Interval vs lead time
for (j in stList) {
  
  Data=preprocdata(directoryCal,stationsID[j],na.rm=TRUE)
  obs=Data[[1]]
  sim=Data[[2]]
  
  obsZoo=converttoZoo(obs)
  simZoo=converttoZoo(sim)
  
  obsData=obs[,2:25]
  simData=sim[,2:25]
  err=obsData-simData
  
  errorModel=quantileCalib(obsData, simData,taus=taus)
  
  #forecast
  dataCheck=preprocdata(directoryVal,stationsID[j],startDate="2010-09-01 01:00",na.rm=TRUE)
  obsCheck=dataCheck[[1]]
  simCheck=dataCheck[[2]]
  
  obsCheckZoo=converttoZoo(obsCheck)
  simCheckZoo=converttoZoo(simCheck)
  
  
  simNew=as.numeric(simCheckZoo[eventDate])
  obsNew=as.numeric(obsCheckZoo[eventDate])
  
  #simNew=simData[500,]
  #obsNew=obsData[500,]
  
  
  png(paste("C:/Users/ci1mmu/QUICS/7_Reference/7_Misc/Thesis/5_Error Model/ConIntPlots/CI_",j, ".png"),width=720,height=480)
  plotCI(errorModel,simNew,obsNew=obsNew, st=j)
  legend("bottomright",bty="n", 
         legend=c("Observation","Deterministic forecast","50% CI","90% CI"),
         col = c("red","navy", "lightblue", "lightblue1"),
         pch=c(1,NA,22,22),lty=c(NA,1,NA,NA), pt.bg=c(NA,NA, "lightblue", "lightblue1"),
         pt.cex=c(1,1,3,3),lwd=c(NA,2,NA,NA),pt.lwd=c(2,NA,NA,NA))
  dev.off()
}

#Quantile regression plots (slope and intersection)

for (j in stList) {
  
  Data=preprocdata(directoryCal,stationsID[j],na.rm=TRUE)
  obs=Data[[1]]
  sim=Data[[2]]
  
  obsZoo=converttoZoo(obs)
  simZoo=converttoZoo(sim)
  
  obsData=obs[,2:25]
  simData=sim[,2:25]
  err=obsData-simData
  
  errorModel=quantileCalib(obsData, simData,taus=taus)
  
  png(paste("C:/Users/ci1mmu/QUICS/7_Reference/7_Misc/Thesis/5_Error Model/QuanRegFitPlots/QRF_",j, ".png"),width=720,height=480*1.5)
  par(mfrow=c(4,2))
  
  for (i in ltList) {
    plotQuanRegFit(errorModel,taus=taus,lt=i)
  }
  
  dev.off()
}

#Extrapolation plots

for (j in stList) {
  
  Data=preprocdata(directoryCal,stationsID[j],na.rm=TRUE)
  obs=Data[[1]]
  sim=Data[[2]]
  
  obsZoo=converttoZoo(obs)
  simZoo=converttoZoo(sim)
  
  obsData=obs[,2:25]
  simData=sim[,2:25]
  err=obsData-simData
  
  errorModel=quantileCalib(obsData, simData,taus=taus)
  
  png(paste("C:/Users/ci1mmu/QUICS/7_Reference/7_Misc/Thesis/5_Error Model/ExtrapPlots/ExPlots_",j, ".png"),width=540,height=480*2)
  par(mfrow=c(2,1))
  
  for (i in ltList){
    plotExtAllData(obsData,simData,errorModel,st=j,lt=i)
  }
  legend("bottomright",legend=c("Linear extrapolation","Log extrapolation"),lty=1,col=c("blue","red"),bty="n")
  dev.off()
}