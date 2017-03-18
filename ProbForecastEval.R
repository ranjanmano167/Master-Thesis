library(xts)
library(abind)
library(quantreg)
library(MASS)
library(verification)



#quantiles from observation

getQuantile=function(directory=directoryObs,stationsID=stationsID,st, quantile=seq(0.1,0.9,0.1)) {
  fileName<-paste0(directory,'obs_', stationsID[st])
  fileName<-strsplit(fileName,'.dfs0')
  fileName<-paste0(fileName,'.csv')
  obsData<-read.csv(fileName, header=FALSE)
  obsData[,1]<-as.POSIXct(obsData[,1], format="%d/%m/%Y %H:%M")
  obsData[obsData<0]<-NA
  obsData=na.omit(obsData)
  quantiles=quantile(obsData[,2], quantile)
  return(quantiles)
}


#Data filtering based on Quantile

dataFilt=function(obs,sim,quantile) {
  lowLimit=quantile(obs,prob=quantile)
  dataComb=data.frame(obs,sim)
  dataCombSort=dataComb[order(dataComb[,1]),]
  dataCombFilt=subset(dataCombSort,dataCombSort[,1]>lowLimit)
  return(dataCombFilt)
}

#Data filtering based on lower threshold
dataFilt2=function(obs,sim,thresLow) {
  lowLimit=thresLow
  dataComb=data.frame(obs,sim)
  dataCombSort=dataComb[order(dataComb[,1]),]
  dataCombFilt=subset(dataCombSort,dataCombSort[,1]>lowLimit)
  return(dataCombFilt)
}

#Data filtering based on upper threshold
dataFilt3=function(obs,sim,thresUp) {
  upLimit=thresUp
  dataComb=data.frame(obs,sim)
  dataCombSort=dataComb[order(dataComb[,1]),]
  dataCombFilt=subset(dataCombSort,dataCombSort[,1]<upLimit)
  return(dataCombFilt)
}


#errorModel for RPS and CRPS calculation
errModelRPS=function(obsData,simData,taus=seq(0.04,0.96,0.04),errModType=1, type=1,extrap=1, weights=0) {
  if(errModType==1) {
    errorModel=quantileCalib(obsData, simData, taus=taus, extrap=extrap, NQTweights=weights)
  } else if (errModType==2) {
    errorModel=quantileCalibOri(obsData, simData, taus=taus, type=type, weights=weights)
  }
  return(errorModel)
}

##PICP and MPI calculation

#errModType=1: ErrorModel in Gaussian domain, 2: ErrrorModel in original domain
#extrapMethod: 1-Linear, 2-Log
#upThresQuantile-upper threshold quantile

PICPPredict=function(obsData,simData,upThresQuantile=0,errorModel,errModType=1,extrapMethod=1) {
  
  n=NCOL(simData)
  
  PICP=NULL
  MPI=NULL
  
  for (j in 1:n) {
    
    obs=obsData[,j]
    sim=simData[,j]
    
    
    filtData=dataFilt(obs=obs,sim=sim,quantile=upThresQuantile)

    obsFilt=filtData[,1]
    simFilt=filtData[,2]
    m=length(simFilt)
        
    if (errModType==1) {
      conBound=sapply(simFilt,function(x) CIEst(errorModel=errorModel,simVal=x,extrap=extrapMethod,lt=j)+x)
    } else if (errModType==2) {
      conBound=sapply(simFilt,function(x) CIEstOri(errorModel=errorModel,simVal=x,lt=j)+x)
    }
    
    
    conBoundFifty=conBound[c(2,4),]
    conBoundNinety=conBound[c(1,5),]
    
    conBoundDiffFifty=sum(apply(conBoundFifty,2,function(x) x[2]-x[1]))
    conBoundDiffNinety=sum(apply(conBoundNinety,2,function(x) x[2]-x[1]))
    
    conBoundfiftyObs=rbind(obsFilt,conBoundFifty)
    conBoundNinetyObs=rbind(obsFilt,conBoundNinety)
    
    rankFifty=apply(conBoundfiftyObs,2,function(x) rank(x,ties.method="first")[1])
    rankNinety=apply(conBoundNinetyObs,2,function(x) rank(x,ties.method="first")[1])
    
    countFifty=sum(rankFifty==2)
    countNinety=sum(rankNinety==2)
    
    PICPFifty=countFifty/m*100
    PICPNinety=countNinety/m*100
    PICPComb=c(PICPFifty,PICPNinety)
    
    MPIFifty=conBoundDiffFifty/m
    MPINinety=conBoundDiffNinety/m
    MPIComb=c(MPIFifty,MPINinety)
    
    PICP=rbind(PICP,PICPComb)  
    MPI=rbind(MPI,MPIComb)
    
    colnames(PICP)=c("PICP-50% CI","PICP-90% CI")
    colnames(MPI)=c("MPI-50% CI","MPI-90% CI")
    
    rownames(PICP)=NULL
    rownames(MPI)=NULL
  }
  return(list(PICP,MPI))
}



#PICP and MPI All
#confBands=taus to predict confidence band from

PICPPredictAll=function(obsData,simData,upThresQuantile=0.7,confBands=c(0.05,0.1,0.2,0.3,0.4,0.6,0.7,0.8,0.9,0.95),errorModel,errModType=1,extrapMethod=1) {
  
  n=NCOL(simData)
  p=length(confBands)
  
  PICP3=NULL
  MPI3=NULL
  
  
  for (j in 1:n) {
    
    obs=obsData[,j]
    sim=simData[,j]
    
    
    filtData=dataFilt(obs=obs,sim=sim,quantile=upThresQuantile)
    
    obsFilt=filtData[,1]
    simFilt=filtData[,2]
    m=length(simFilt)
    
    if (errModType==1) {
      conBound=sapply(simFilt,function(x) CIEst(errorModel,simVal=x,extrap=extrapMethod,lt=j)+x)
    } else if (errModType==2) {
      conBound=sapply(simFilt,function(x) CIEstOri(errorModel,simVal=x,lt=j)+x)
    }
    
    PICP2=NULL
    MPI2=NULL
    
    for (k in 1:(p/2)) {
      
      conBound1=conBound[c(k,(p-k+1)),]
      conBoundDiff=sum(apply(conBound1,2,function(x) x[2]-x[1]))
      conBoundObs=rbind(obsFilt,conBound1)
      rankCB=apply(conBoundObs,2,function(x) rank(x,ties.method="first")[1])
      countCB=sum(rankCB==2)
      
      PICP1=countCB/m*100
      MPI1=conBoundDiff/m
      
      PICP2=rbind(PICP2,PICP1)  
      MPI2=rbind(MPI2,MPI1)
      
      rownames(PICP2)=NULL
      rownames(MPI2)=NULL
    }
    PICP3=cbind(PICP3,PICP2)
    MPI3=cbind(MPI3,MPI2)
    
  }
  
  return(list(confBands,PICP3,MPI3))
}



#CRPS data
CRPSData=function(obs,sim,errorModel,errModType=1,type=1, extrap=1,quantile,leadTime) {
  dataNew=dataFilt(obs,sim,quantile)
  obsNew=dataNew[,1]
  simNew=dataNew[,2]
  
  if (errModType==1) {
    conBound=sapply(simNew,function(x) CIEst(errorModel,simVal=x,extrap=extrap,lt=leadTime)+x)
  } else if (errModType==2) {
    conBound=sapply(simNew,function(x) CIEstOri(errorModel,simVal=x, type=type, lt=leadTime)+x)
  }
  meanSimNorm=conBound[2,]
  SDSimNorm=(((conBound[2,]-conBound[1,])/2)+((conBound[3,]-conBound[2,])/2))/2
  return(cbind(obsNew,meanSimNorm,SDSimNorm))
  
}

#RPSS data
RPSData=function(obs,sim,simBins,errorModel,errModType=1,type=1, extrap=1,quantile,leadTime) {
  dataNew=dataFilt(obs,sim,quantile)
  obsNew=dataNew[,1]
  simNew=dataNew[,2]
  
  if (errModType==1) {
    conBound=sapply(simNew,function(x) CIEst(errorModel,simVal=x,extrap=extrap,lt=leadTime)+x)
  } else if (errModType==2) {
    conBound=sapply(simNew,function(x) CIEstOri(errorModel,simVal=x, type=type, lt=leadTime)+x)
  }
  
  simBinsMatrix=matrix(rep(simBins,length(obsNew)),nrow=length(simBins))
  simBinsObs=rbind(obsNew,simBinsMatrix)
  obsCat=apply(simBinsObs,2,function(x) rank(x,ties.method="max")[1])
  probMatrix=matrix(0,ncol=length(obsNew),nrow=(length(simBins)+1))
  
  for (i in 1:NCOL(conBound)) {
    for (j in 1:NROW(conBound)){
      CONBOUND=conBound[j,i]
      simBins2=c(CONBOUND,simBins)
      CONBOUNDRank=rank(simBins2,ties.method="max")[1]
      probMatrix[CONBOUNDRank,i]=probMatrix[CONBOUNDRank,i]+1
    }
  }
  
  PROBMATRIX=round(probMatrix/nrow(conBound),2)
  
  return(list(obsCat,PROBMATRIX))
  
}


#reliability diagram data
relDiaData=function(obs,sim,eventthres=100,probBins=seq(0,1,0.1),errorModel,errModType=1,type=1, extrap=1,leadTime) {
  
  if (errModType==1) {
    conBound=sapply(sim,function(x) CIEst(errorModel,simVal=x,extrap=extrap,lt=leadTime)+x)
  } else if (errModType==2) {
    conBound=sapply(sim,function(x) CIEstOri(errorModel,simVal=x, type=type, lt=leadTime)+x)
  }
  simYes=apply(conBound,2,function(x) sum(x>eventthres)/nrow(conBound))
  simNo=1-simYes
  
  obsYes=(obs>eventthres)*1
  obsYesFac=factor(obsYes,levels=c(1,0))
  simYesBin=cut(simYes,breaks=probBins,include.lowest = TRUE)
  obsSimYes=data.frame(simYesBin,obsYesFac)
  sumTable=table(obsSimYes)
  sumTableMatrix=as.data.frame.matrix(sumTable)
  sampleSize=sumTableMatrix[,1]+sumTableMatrix[,2]
  samplesizePer=sampleSize/sum(sampleSize)
  obsFreq=round(sumTableMatrix[,1]/sampleSize,3)
  midPointProb=(tail(probBins,-1)+head(probBins,-1))/2
  
  ReliabData=cbind(midPointProb,sumTableMatrix[,1:2],obsFreq,sampleSize,samplesizePer)
  colnames(ReliabData)=c("Fore Prob","Yes","No","Obs Freq","Sample Size","Sample Size %")
  
  return(ReliabData)
  
}

#Brier skillscore data
brierSSData=function(obs,sim,eventthres=100,errorModel,errModType=1,type=1, extrap=1,leadTime) {
  
  if (errModType==1) {
    conBound=sapply(sim,function(x) CIEst(errorModel,simVal=x,extrap=extrap,lt=leadTime)+x)
  } else if (errModType==2) {
    conBound=sapply(sim,function(x) CIEstOri(errorModel,simVal=x, type=type, lt=leadTime)+x)
  }
  
  simYes=apply(conBound,2,function(x) sum(x>eventthres)/nrow(conBound))
  obsYes=(obs>eventthres)*1
  
  brierSSData=cbind(simYes,obsYes)
  
  return(brierSSData)
  
}
#Bin

# #Inverse ecdf to interpolate for given x value
# invEcdf <- function(f){
#   x <- environment(f)$x
#   y <- environment(f)$y
#   approxfun(y, x,rule=2)
# }
# 
# #Random forecasts
# foreRandom=function(obs){
#   funEcdf=ecdf(obs)
#   foreRan=invEcdf(funEcdf)(runif(length(obs)))
#   return(foreRan)
# }

# #Reliability diagram
# dataReliabPlot=function(obsData,simData,errorModel,errModType=1,extrapMethod=1,leadTime=1) {
#   
#   obs=obsData[,leadTime]
#   sim=simData[,leadTime]
#   m=NROW(sim)
#   
#   if (errModType==1) {
#     nonExcBound=sapply(sim,function(x) CIEst(errorModel,simVal=x,extrap=extrapMethod,lt=j)+x)
#   } else if (errModType==2) {
#     nonExcBound=sapply(sim,function(x) CIEstOri(errorModel,simVal=x,lt=j)+x)
#   }
#   
#   nEBoundObs=rbind(obs,nonExcBound)
#   
#   rankObs=apply(nEBoundObs,2,function(x) rank(x,ties.method="first")[1])
#   
#   rankFactor=factor(rankObs,levels=1:(nrow(nonExcBound)+1))
#   
#   rankTable=table(rankFactor)
#   
#   probTable=rankTable/m
#   
#   obsFreq=cumsum(probTable)
#   
#   return(list(errorModel[[1]],obsFreq))
# }