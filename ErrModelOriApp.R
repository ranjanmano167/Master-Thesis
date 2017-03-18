source("C:/Data/Thesis/4_Performance Measures/R_Scripts/Gen_util.R")
source("C:/Data/Thesis/5_Error Model/Scripts/ErrorModelOri.R")
directoryCal="C:/Data/Thesis/4_Performance Measures/Data/Calibration/"
directoryVal="C:/Data/Thesis/4_Performance Measures/Data/Validation/"
stationFileDir="C:/Data/Thesis/4_Performance Measures/Station_Data.csv"
Sys.setenv(TZ='GMT')

stationsData = read.table(stationFileDir,header=FALSE, sep=",")
stationsID=as.character(stationsData[[1]])

taus=c(0.05,0.25, 0.5, 0.75,0.95)
#taus=seq(0.04,0.96,0.04)
#taus=c(0.025,seq(0.05,0.45,0.05),0.475,0.525,seq(0.55,0.95,0.05),0.975)
stList=c(1)
ltList=c(6)
eventDate="2010-09-19 13:00:00"
simList=c(seq(0,2300,100),3000,4000,5000)


#Sim Vs Err
for (j in stList) {
  
  Data=preprocdata(directoryCal,stationsID[j],na.rm=TRUE)
  obs=Data[[1]]
  sim=Data[[2]]
  
  obsZoo=converttoZoo(obs)
  simZoo=converttoZoo(sim)
  
  obsData=obs[,2:25]
  simData=sim[,2:25]
  err=obsData-simData
  
  errorModel=quantileCalibOri(obs=obsData, sim=simData,taus=taus,weights=0)
  
  png(paste("C:/Data/Thesis/5_Error Model/SimVsErrPlots/Ris_",j, ".png"),width=720,height=480)
  par(mfrow=c(1,1))
  
  for (i in ltList) { 
    ylimit=max(abs(min(err)),abs(max(err)))
    
    plot(simData[,i],err[,i],
         main=(paste("Quantile Regression, St:",j,", LT (hrs):",i)),
         xlim=c(0,max(simData)),ylim=c(-ylimit,ylimit),
         #xlim=c(600,1000),ylim=c(-100,100),
         xlab="Predicted discharge [m3/s]",ylab="Residuals [m3/s]",
         col="gray60")
    plotSimVsErrOri(errorModel,sim=simData,err=err,st=j,lt=i)
    
  }
  legend("bottomleft",bty="n", 
         legend=c("Median","50% CI","90% CI"),
         col ="blue", lty=c(1,2,3),lw=2)
  
  dev.off()
}



#Confidence interval vs events

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
  
  errorModel=quantileCalibOri(obsData, simData,taus=taus)
  
  #forecast
  dataCheck=preprocdata(directoryVal,stationsID[j],startDate="2010-09-01 01:00",na.rm=TRUE)
  obsCheck=dataCheck[[1]]
  simCheck=dataCheck[[2]]
  
  obsCheckData=obsCheck
  simCheckData=simCheck
  
  png(paste("C:/Data/Thesis/5_Error Model/FloodEventsPlots/FE_",j, ".png"),width=720*1.5,height=480*1.5)
  par(mfrow=c(2,2))
  for (i in ltList){
    plotCIltOri(errorModel,simAll=simCheckData,obsAll=obsCheckData, st=j,lt=i)
  }
  legend("topleft",bty="n", 
         legend=c("Deterministic forecast","50% CI","90% CI","Observation"),
         col = c("navy", "red", "hotpink","palevioletred4"),
         pch=c(NA,22,22,16),lty=c(1,NA,NA,NA),pt.bg=c(NA, "red", "hotpink",NA),
         pt.cex=c(1,3,3,1))
  dev.off()
}

#Confidence Interval plots
for (j in stList) {
  
  Data=preprocdata(directoryCal,stationsID[j],na.rm=TRUE)
  obs=Data[[1]]
  sim=Data[[2]]
  
  obsZoo=converttoZoo(obs)
  simZoo=converttoZoo(sim)
  
  obsData=obs[,2:25]
  simData=sim[,2:25]
  err=obsData-simData
  
  errorModel=quantileCalibOri(obsData, simData,taus=taus)
  
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
  
  
  png(paste("C:/Data/Thesis/5_Error Model/ConIntPlots/CI_",j, ".png"),width=720,height=480)
  plotCIOri(errorModel,simNew=simNew,obsNew=obsNew, st=j)
  legend("bottomright",bty="n", 
         legend=c("Deterministic forecast","50% CI","90% CI"),
         col = c("navy", "lightblue", "lightblue1"),
         pch=c(NA,22,22),lty=c(1,NA,NA), pt.bg=c(NA, "lightblue", "lightblue1"),
         pt.cex=c(1,3,3),lwd=c(2,NA,NA))
  dev.off()
}

#Quantile regression plots

for (j in stList) {
  
  Data=preprocdata(directoryCal,stationsID[j],na.rm=TRUE)
  obs=Data[[1]]
  sim=Data[[2]]
  
  obsZoo=converttoZoo(obs)
  simZoo=converttoZoo(sim)
  
  obsData=obs[,2:25]
  simData=sim[,2:25]
  err=obsData-simData
  
  errorModel=quantileCalibOri(obsData, simData,taus=taus)
  
  png(paste("C:/Data/Thesis/5_Error Model/QuanRegFitPlots/QRF_",j, ".png"),width=720,height=480*1.5)
  par(mfrow=c(4,2))
  
  for (i in ltList) {
    plotQuanRegFitOri(errorModel,taus=taus,lt=i)
  }
  
  dev.off()
}

