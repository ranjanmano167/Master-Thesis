library(xts)
library(abind)
library(quantreg)
library(MASS)
#---------------Quantile regression in Gaussian domain-------------------
#v1-->v2: 
#1.Solved crossing problem by making NQT risiduals constant just before any crossing may occur

#v2-->v3:
#1.Introduction of log extrapolation


#--------calibration of error model------

##taus: Required quantiles
##NQTWeights:weights for the NQT simulated values, 0-No weigthts, 1-Rank based weights, 2-Value based weights 
##calibNQTLowerLim: Lower limit below which the simulated values will not be considered for QR 
##extrap:extrapolation method, 1-linear regression, 2-log regression

extrapMethod=1

quantileCalib = function (obsData, simData, taus,interpNumber=1000,calibNQTLowerLim=NA, extrap=extrapMethod, NQTweights=0, saveFilePath=NULL) {
  
  # default quantiles if NULL is specified
  
  if (length(taus)==0){ 
    taus=c(.05, .25, 0.5, .75, .95)
  }
  
  err=obsData-simData
  
  #Q-Q plot points
  simNorm=apply(simData, 2, function(x) qqnorm(sort(x), plot.it = F))
  errNorm=apply(err, 2, function(x) qqnorm(sort(x), plot.it = F))
  
  #interpolation
  simAppNorm<-mapply(function(x,y) approx(x$y,y$x,method = "linear",n=interpNumber),simNorm,simNorm)
  errAppNorm<-mapply(function(x,y) approx(x$y,y$x,method = "linear",n=interpNumber),errNorm,errNorm)
  
  
  #3d array
  simNQTMap=NULL
  errNQTMap=NULL
  
  
  for (i in 1:length(simAppNorm[1,])){
    simNQTMap=abind(simNQTMap, matrix(unlist(simAppNorm[,i]), ncol=2), along=3)
    errNQTMap=abind(errNQTMap, matrix(unlist(errAppNorm[,i]), ncol=2), along=3)
  } 
  
  
  # Linear extrapolation of mapping tables
  extrapSimUp=NULL  # Upper limit of simulated NQT map
  extrapSimLow=NULL # Lower end of simulated NQT map
  extrapErrUp=NULL  # Upper end of error NQT map
  extrapErrLow=NULL # lower end of error NQT map
  
  
  tailPointsSim=10 # Points used in linear regression of tail points in simulated data
  tailPointsErr=10 # Points used in linear regression of tail points in error
  maxiter=50
  
  
  for (i in 1:length(simNorm)){
    extrapSimLow=rbind(extrapSimLow,rlm(head(simNorm[[i]][[1]],tailPointsSim)~head(simNorm[[i]][[2]],tailPointsSim),maxit=maxiter)[[1]])
    extrapErrLow=rbind(extrapErrLow, rlm(head(errNorm[[i]][[2]],tailPointsErr)~head(errNorm[[i]][[1]],tailPointsErr),maxit=maxiter)[[1]])  
    
    if (extrap==1) {
      extrapSimUp=rbind(extrapSimUp, rlm(tail(simNorm[[i]][[1]],tailPointsSim)~tail(simNorm[[i]][[2]],tailPointsSim),maxit=maxiter)[[1]])
      extrapErrUp=rbind(extrapErrUp, rlm(tail(errNorm[[i]][[2]],tailPointsErr)~tail(errNorm[[i]][[1]],tailPointsErr),maxit=maxiter)[[1]])
      
    } else if (extrap==2) {
      extrapSimUp=rbind(extrapSimUp, lm(tail(simNorm[[i]][[1]],tailPointsSim)~log(tail(simNorm[[i]][[2]],tailPointsSim)))[[1]])
      extrapErrUp=rbind(extrapErrUp, lm(tail(errNorm[[i]][[2]],tailPointsErr)~exp(tail(errNorm[[i]][[1]],tailPointsErr)))[[1]])  
    }
    
  }
  
  
  # joining to one matrix
  extrapErr=matrix(c(extrapErrLow,extrapErrUp), ncol=4)
  extrapSim=matrix(c(extrapSimLow,extrapSimUp), ncol=4)
  
  colnames(extrapSim)=c("isect - Low","slope - Low", "isect - up","slope - up")
  colnames(extrapErr)=c("isect - Low","slope - Low", "isect - up","slope - up")
  
  qrCoeffs=NULL
  
  for (i in 1:length(errAppNorm[1,])){
    # mapping to normal space
    
    simNormData=approx(simAppNorm[[1,i]], simAppNorm[[2,i]], xout=simData[,i] ,method= "linear")
    errNormData=approx(errAppNorm[[1,i]], errAppNorm[[2,i]], xout=err[,i] ,method= "linear")
    
    
    #linear quantile regresiion 
    
    #setting a threshold value to fit QR only above that value
    
    simNormDatafilt=simNormData$y
    simNormDatafilt[simNormDatafilt<calibNQTLowerLim]=NA
    errNormDatafilt=errNormData$y
    errNormDatafilt[is.na(simNormDatafilt)]=NA
    
    #Introducing weights 
    
    simWeights=simNormDatafilt
    simWeights[simWeights<(0)]=0
    simWeights[simWeights>(0)]=1
    
#     if (NQTweights==2){
#       simWeights=(simNormDatafilt+abs(min(simNormDatafilt)))/(max(simNormDatafilt)-min(simNormDatafilt))
#     } else if (NQTweights==1) {
#       simWeights=rank(simNormDatafilt)/length(simNormDatafilt)
#     } else if (NQTweights==0){
#       simWeights=NULL
#     }
    
    if (i>1){
      k=coef(rq(errNormDatafilt ~  simNormDatafilt, weights=simWeights, tau = taus))
      rownames(k) = c("Intersept", "Slope")
      qrCoeffs=abind(qrCoeffs,k, along=3)
    }
    
    else {
      qrCoeffs=coef(rq(errNormDatafilt ~ simNormDatafilt, weights=simWeights, tau = taus))
      rownames(qrCoeffs) = c("Intersept", "Slope")
    }
  }
  
  #result
  errModel=list(taus, qrCoeffs,simNQTMap,  errNQTMap, extrapSim, extrapErr)
  
  if (is.null(saveFilePath)==FALSE){
    try(save(errModel, file=saveFilePath))
  }
  return(errModel)
}

#-------application/validation of error model------
#lt-leadtime
#sdLimit=Smoothing boundaries at lower and upper levels before extrapolation
#noncrossing: If crossing need to be avoided or not (TRUE/FALSE)

CIEst=function(errorModel,simVal,lt,extrap=extrapMethod,sdLimit=3,nonCrossing=FALSE) {
  
  simNQTMap=errorModel[[3]][,,lt]
  errNQTMap=errorModel[[4]][,,lt]
  rqFits=errorModel[[2]][,,lt]
  extrapSim=errorModel[[5]][lt,]
  extrapErr=errorModel[[6]][lt,]
  
  if (nonCrossing){
    m=NCOL(rqFits)
    NQTcutAll=NULL
    for (i in 1:(m-1)){
      NQTcut=(rqFits[1,i]-rqFits[1,(i+1)])/(rqFits[2,(i+1)]-rqFits[2,i])
      NQTcutAll=cbind(NQTcutAll,NQTcut)
    }
    NQTLowerLim=(max(NQTcutAll[NQTcutAll<0])+0.1)
  } else {
    NQTLowerLim=-10000
  }
  
  
  #Transformation of simulated values to Gaussian
  simNor=approx(simNQTMap[,1],simNQTMap[,2],xout=simVal,method="linear")
  
  #Extrapolation of outside values to gaussian domain
  
  if(is.na(simNor$y)) {
    
    
    if(simNor$x>max(simNQTMap[,1])){
      
      if(extrap==1){
        simNor$y=simNor$x*extrapSim[[4]]+extrapSim[[3]]
      } else if (extrap==2){
        simNor$y=log(simNor$x)*extrapSim[[4]]+extrapSim[[3]]
      }
      
    } else if(simNor$x<min(simNQTMap[,1])){
      simNor$y=simNor$x*extrapSim[[2]]+extrapSim[[1]]
    }
  }
  
  #Error in Gaussian domain 
  errAll=NULL
  
  for (j in 1:NCOL(rqFits)){
    if (simNor$y>NQTLowerLim){ 
      errNorm=(simNor$y*rqFits[2,j]+rqFits[1,j])
    }
    else {
      errNorm=(NQTLowerLim*rqFits[2,j]+rqFits[1,j])
    }
    
    #Inverse transformation
    errNew=approx(errNQTMap[,2],errNQTMap[,1],xout=errNorm,method='linear',rule=1)
    
    #Extrapolation to original domain
    if(is.na(errNew$y)) {
      if(errNew$x>max(errNQTMap[,2])){
        
        if (extrap==1) {
          errNew$y=errNorm*extrapErr[[4]]+extrapErr[[3]]
        } else if (extrap==2) {
          errNew$y=exp(errNorm)*extrapErr[[4]]+extrapErr[[3]]
        }
        
      } else {
        errNew$y=errNorm*extrapErr[[2]]+extrapErr[[1]]
      }
    }
    
    #Smoothing at upper and lower levels
    if(errNew$x>sdLimit) {
      extrapWeight=(errNew$x-sdLimit)/(max(errNQTMap[,2])-sdLimit)
      if (extrap==1) {
        errNew$y=errNew$y*(1-extrapWeight)+(extrapErr[[4]]*errNew$x+extrapErr[[3]])*extrapWeight
      } else if (extrap==2) {
        errNew$y=errNew$y*(1-extrapWeight)+(extrapErr[[4]]*exp(errNew$x)+extrapErr[[3]])*extrapWeight
      }
    } else if(errNew$x<(-sdLimit)) {
      extrapWeight=(errNew$x+sdLimit)/(min(errNQTMap[,2])+sdLimit)
      errNew$y=errNew$y*(1-extrapWeight)+(extrapErr[[2]]*errNew$x+extrapErr[[1]])*extrapWeight
    }
    
    errAll=cbind(errAll,errNew$y)
  }
  
  return(errAll)
  
}

#plots

#1. Normalised Sim Vs Normalised Err with quantile regression lines

plotNorSimVsErr=function(errormodel,st=NULL, lt,nonCrossing=FALSE)  {
  
  rqFits=errorModel[[2]][,,lt]
  simNQTMap=errorModel[[3]][,,lt]
  errNQTMap=errorModel[[4]][,,lt]
  extrapSim=errorModel[[5]][lt,]
  extrapErr=errorModel[[6]][lt,]
  m=NCOL(rqFits)
  ltyp=c(4,3,2,1,2,3,4)
  n=(7-m)/2 #used to choose line type
  
  simNormData=approx(simNQTMap[,1], simNQTMap[,2], xout=simData[,lt], method= "linear")
  errNormData=approx(errNQTMap[,1], errNQTMap[,2], xout=err[,lt], method= "linear")
  
  
  plot(simNormData$y,errNormData$y,
       main=(paste("Quantile calib., St:",st,", LT (hrs):",lt)),
       xlab="NQT:Predicted discharge [-]",ylab="NQT:Residuals [-]",
       col="gray60")
  
  for (j in 1:m) {
    a=rqFits[1,j]
    b=rqFits[2,j]
    
    if (nonCrossing){
      
      NQTcutAll=NULL
      
      for (i in 1:(m-1)){
        NQTcut=(rqFits[1,i]-rqFits[1,(i+1)])/(rqFits[2,(i+1)]-rqFits[2,i])
        NQTcutAll=cbind(NQTcutAll,NQTcut)
      }
      
      NQTLowerLim=max((max(NQTcutAll[NQTcutAll<0])+0.1),min(simNormData$y))
      
      xx=c(NQTLowerLim,max(simNormData$y))
      yy=c((b*NQTLowerLim+a),(b*max(simNormData$y)+a))
      lines(xx,yy,col='blue',lty=ltyp[j+n],lw=2)
      #lines(xx,yy,col='blue',lty=1,lw=2)
      
      xxx=c(min(simNormData$y),NQTLowerLim)
      yyy=c((b*NQTLowerLim+a),(b*NQTLowerLim+a))
      lines(xxx,yyy,col='blue',lty=ltyp[j+n],lw=2)
      #lines(xxx,yyy,col='blue',lty=1,lw=2)
      
    } else {
      xx=c(min(simNormData$y),max(simNormData$y))
      yy=c((b*min(simNormData$y)+a),(b*max(simNormData$y)+a))
      lines(xx,yy,col='blue',lty=ltyp[j+n],lw=2)
      #lines(xx,yy,col='blue',lty=1,lw=2)
    }
  }
}


#2. confidence Interval Plots

#2.1 CI vs leadtime
plotCI=function(errorModel,simNew,obsNew=NULL, st=NULL) {
  
  m=length(simNew)
  n=NCOL(errorModel[[2]])
  
  CI=NULL
  for (i in 1:m){
    CILt=CIEst(errorModel,simVal=simNew[i],lt=i)
    CI=rbind(CI,CILt)
  }
  
  a=(CI[,(n-1)/2]+as.numeric(simNew))
  b=(CI[,(n+3)/2]+as.numeric(simNew))
  c=(CI[,(n-3)/2]+as.numeric(simNew))
  d=(CI[,(n+5)/2]+as.numeric(simNew))
  e=as.numeric(obsNew)
  ymin=0 #min(a,b,c,d,e)
  ymax=max(1.5*max(simNew),max(obsNew))#max(a,b,c,d,e)
  
  plot(1:m,simNew,col='red',ylim=c(ymin,ymax),type='n', main=paste("Forecast with confidence intervals, St:",st), 
       ylab="Discharge[m3/s]", xlab="Lead time[hours]")
  polygon(c(1:m,rev(1:m)),c(a,rev(b)),col='lightblue',border = NA)
  polygon(c(1:m,rev(1:m)),c(a,rev(c)),col='lightblue1',border = NA)
  polygon(c(1:m,rev(1:m)),c(b,rev(d)),col='lightblue1',border = NA)
  
  points(1:m,simNew,type='l',col='navy',lw=2)
  
  #for (i in 1:n) {
  #  x=(CI[,i]+as.numeric(simNew))
  #  points(1:m,x,type='l')
  #}
  
  if(is.null(obsNew)=='FALSE'){
    points(1:m,e,type='p',col='red',lw=2)
  }
  
}

#2.2 CI vs events
plotCIlt=function(errorModel,simAll,obsAll=NULL, st=NULL,lt) {
  
  simNew=simAll[,2:25]
  obsNew=obsAll[,2:25]
  
  m=NROW(simNew)
  n=NCOL(errorModel[[2]])
  
  CI=NULL
  for (i in 1:m){
    CIAll=CIEst(errorModel,simVal=simNew[i,lt],lt=lt)
    CI=rbind(CI,CIAll)
  }
  
  a=(CI[,(n-1)/2]+as.numeric(simNew[,lt]))
  b=(CI[,(n+3)/2]+as.numeric(simNew[,lt]))
  c=(CI[,(n-3)/2]+as.numeric(simNew[,lt]))
  d=(CI[,(n+5)/2]+as.numeric(simNew[,lt]))
  e=as.numeric(obsNew[,lt])
  ymin=0 #min(a,b,c,d,e)
  ymax=max(1.25*max(simNew),max(obsNew))#max(a,b,c,d,e)
  xVal=simAll[,1]+(lt*3600)
  
  plot(xVal,simNew[,lt],col='red',ylim=c(ymin,ymax),type='n', 
       main=paste("Forecast with confidence intervals, St:",st,", LT (hrs):",lt), 
       ylab="Discharge[m3/s]", xaxt = "n",xlab="Date")
  axis.POSIXct(1,xVal,format = "%d/%m/%y")
  grid(nx=NULL,ny=NULL)
  polygon(c(xVal,rev(xVal)),c(a,rev(b)),col='lightblue1',border = NA)
  polygon(c(xVal,rev(xVal)),c(a,rev(c)),col='lightblue1',border = NA)
  polygon(c(xVal,rev(xVal)),c(b,rev(d)),col='lightblue1',border = NA)
 
#   polygon(c(xVal,rev(xVal)),c(a,rev(b)),col='green3',border = NA)
#   polygon(c(xVal,rev(xVal)),c(a,rev(c)),col='green1',border = NA)
#   polygon(c(xVal,rev(xVal)),c(b,rev(d)),col='green1',border = NA)
  
  points(xVal,simNew[,lt],type='l',col='navy')
  
  #for (i in 1:n) {
  #  x=(CI[,i]+as.numeric(simNew))
  #  points(1:m,x,type='l')
  #}
  
  if(is.null(obsNew)=='FALSE'){
    points(xVal,e,type='p',col='palevioletred4',pch=16)
  }
  
}

#2.3 CI vs events lines only
plotCIltLines=function(errorModel,simAll,obsAll=NULL, st=NULL,lt) {
  
  simNew=simAll[,2:25]
  obsNew=obsAll[,2:25]
  
  m=NROW(simNew)
  n=NCOL(errorModel[[2]])
  
  CI=NULL
  for (i in 1:m){
    CIAll=CIEst(errorModel,simVal=simNew[i,lt],lt=lt)
    CI=rbind(CI,CIAll)
  }
  
  a=(CI[,(n-1)/2]+as.numeric(simNew[,lt]))
  b=(CI[,(n+3)/2]+as.numeric(simNew[,lt]))
  c=(CI[,(n-3)/2]+as.numeric(simNew[,lt]))
  d=(CI[,(n+5)/2]+as.numeric(simNew[,lt]))
  e=as.numeric(obsNew[,lt])
  ymin=0 #min(a,b,c,d,e)
  ymax=max(1.25*max(simNew),max(obsNew))#max(a,b,c,d,e)
  xVal=simAll[,1]+(lt*3600)
  
#   plot(xVal,simNew[,lt],col='red',ylim=c(ymin,ymax),type='n', 
#        main=paste("Forecast with confidence intervals, St:",st,", LT (hrs):",lt), 
#        ylab="Discharge[m3/s]", xaxt = "n",xlab="Date")
#   axis.POSIXct(1,xVal,format = "%d/%m/%y")
#   points(xVal,a,col='red',type="l",lwd=2)
#   points(xVal,b,col='red',type="l",lwd=2)
  points(xVal,c,col='red',type="l",lwd=2)
  points(xVal,d,col='red',type="l",lwd=2)
#   points(xVal,simNew[,lt],col='green4',type='l',lwd=2)
  
  #for (i in 1:n) {
  #  x=(CI[,i]+as.numeric(simNew))
  #  points(1:m,x,type='l')
  #}
  
#   if(is.null(obsNew)=='FALSE'){
#     points(xVal,e,type='p',col='red',pch=16)
#   }
  
}

#3. Plots of  quantile slope and gradient

plotQuanRegFit=function(errorModel,taus, lt) {
  rqFits=errorModel[[2]][,,lt]
  plot(taus,rqFits[2,],type='b',xlim=c(0,1),ylim=c(-1,1),xlab='Quantile',ylab='Slope',main=paste('Quantile regression fit, LT (hrs):',lt))
  plot(taus,rqFits[1,],type='b',xlim=c(0,1),ylim=c(-2,2),xlab='Quantile',ylab='Intersection')
}

#4.Plots of error risiduals in original domain
plotErrRisOri=function(errorModel,simList,lt,colour="blue") {
  rqFits=errorModel[[2]][,,lt]
  risSim=NULL
  for(i in 1:length(simList)) {
    risSim=rbind(risSim,CIEst(errorModel,simList[i],lt))
  }
  for (j in 1:ncol(risSim)) {
    
    #colo=c("lightblue1","lightblue","mediumblue","navy","mediumblue","lightblue","lightblue1")
    ltyp=c(4,3,2,1,2,3,4)
    n=(7-NCOL(rqFits))/2
    
    points(simList,risSim[,j],type="l",col=colour,lty=ltyp[j+n],lw=2)
  }
}


#5.Extrapolation plots
plotExt=function(obsData,simData,errorModel,extrap=extrapMethod, st=NULL,lt) {
  
  err=obsData-simData
  
  #Q-Q plot points
  simNorm=apply(simData, 2, function(x) qqnorm(sort(x), plot.it = F))
  errNorm=apply(err, 2, function(x) qqnorm(sort(x), plot.it = F))
  
  #plot-discharge
  plot(tail(simNorm[[lt]][[2]],10),tail(simNorm[[lt]][[1]],10),
       xlab="Predicted discharge [m3/s]",ylab="NQT:Predicted discharge  [-]",
       main=paste("Sim. Extrapolation , St:",st,", LT (hrs):",lt))
  if (extrap==1){
    abline(a=errorModel[[5]][lt,3],b=errorModel[[5]][lt,4],col="blue")
  } else if (extrap==2) {
    points(tail(simNorm[[lt]][[2]],10),(errorModel[[5]][lt,3]+(log(tail(simNorm[[lt]][[2]],10))*errorModel[[5]][lt,4])),
           type="l",col="red")
  }
  
  
  #plot-error
  plot(tail(errNorm[[lt]][[1]],10),tail(errNorm[[lt]][[2]],10),
       xlab="NQT-Risiduals [-]",ylab="Risiduals  [m3/s]",
       main=paste("Error Extrapolation , St:",st,", LT (hrs):",lt))
  if (extrap==1){
    abline(a=errorModel[[6]][lt,3],b=errorModel[[6]][lt,4],col="blue")
  } else if (extrap==2) {
    points(tail(errNorm[[lt]][[1]],10),(errorModel[[6]][lt,3]+(exp(tail(errNorm[[lt]][[1]],10))*errorModel[[6]][lt,4])),
           type="l",col="red")
  }
}

#6.Extrapolation plots -two types in one plot
plotExtNew=function(obsData,simData,errorModel,st=NULL,lt) {
  
  err=obsData-simData
  
  #Q-Q plot points
  simNorm=apply(simData, 2, function(x) qqnorm(sort(x), plot.it = F))
  errNorm=apply(err, 2, function(x) qqnorm(sort(x), plot.it = F))
  
  #simNorm=apply(simData, 2, function(x) qqnorm(c(sort(x),1200), plot.it = F))
  #errNorm=apply(err, 2, function(x) qqnorm(c(sort(x),200), plot.it = F))
  
  #plot-discharge
  plot(tail(simNorm[[lt]][[2]],10),tail(simNorm[[lt]][[1]],10),
       xlab="Predicted discharge [m3/s]",ylab="NQT:Predicted discharge  [-]",
       main=paste("Sim. Extrapolation , St:",st,", LT (hrs):",lt))
  
  #abline(a=errorModel[[5]][lt,3],b=errorModel[[5]][lt,4],col="blue")
  abline(a=errorModel[[5]][lt,3],b=errorModel[[5]][lt,4],col="blue")
  logFit=predict(lm( tail(simNorm[[lt]][[1]],10)~log(tail(simNorm[[lt]][[2]],10)) ))
  lines(tail(simNorm[[lt]][[2]],10),logFit,col='red')
  
  #plot-error
  plot(tail(errNorm[[lt]][[1]],10),tail(errNorm[[lt]][[2]],10),
       xlab="NQT-Risiduals [-]",ylab="Risiduals  [m3/s]",
       main=paste("NQT-Error , St:",st,", LT (hrs):",lt))
  
  abline(a=errorModel[[6]][lt,3],b=errorModel[[6]][lt,4],col="blue")
  expFit=predict(lm( tail(errNorm[[lt]][[2]],10)~exp(tail(errNorm[[lt]][[1]],10))))
  lines(tail(errNorm[[lt]][[1]],10),expFit,col='red')
}


#7.Extrapolation plots all data
plotExtAllData=function(obsData,simData,errorModel,extrap=extrapMethod, st=NULL,lt) {
  
  err=obsData-simData
  
  #Q-Q plot points
  simNorm=apply(simData, 2, function(x) qqnorm(sort(x), plot.it = F))
  errNorm=apply(err, 2, function(x) qqnorm(sort(x), plot.it = F))
  
  #plot-discharge
  plot(simNorm[[lt]][[2]],simNorm[[lt]][[1]],
       xlab="Predicted discharge [m3/s]",ylab="NQT:Predicted discharge  [-]",
       main=paste("NQT-Discharge , St:",st,", LT (hrs):",lt))
  grid(nx=NULL,ny=NULL)
#   if (extrap==1){
#     abline(a=errorModel[[5]][lt,3],b=errorModel[[5]][lt,4],col="blue")
#   } else if (extrap==2) {
#     points(tail(simNorm[[lt]][[2]],10),(errorModel[[5]][lt,3]+(log(tail(simNorm[[lt]][[2]],10))*errorModel[[5]][lt,4])),
#            type="l",col="red")
#   }
  
  
  #plot-error
  plot(errNorm[[lt]][[1]],errNorm[[lt]][[2]],
       xlab="NQT-Risiduals [-]",ylab="Risiduals  [m3/s]",
       main=paste("NQT-Error , St:",st,", LT (hrs):",lt))

  grid(nx=NULL,ny=NULL)
#   if (extrap==1){
#     abline(a=errorModel[[6]][lt,3],b=errorModel[[6]][lt,4],col="blue")
#   } else if (extrap==2) {
#     points(tail(errNorm[[lt]][[1]],10),(errorModel[[6]][lt,3]+(exp(tail(errNorm[[lt]][[1]],10))*errorModel[[6]][lt,4])),
#            type="l",col="red")
#   }
}





##Bin

# #5.Extrapolation plots
# plotExt=function(obsData,simData,errorModel,st=NULL,lt) {
#   
#   err=obsData-simData
#   simExtrap=tail(sort(simData),5)
#   errExtrap=tail(sort(err),5)
#   
#   #Q-Q plot points
#   simNorm=apply(simData, 2, function(x) qqnorm(sort(x), plot.it = F))
#   simNormEx=apply(simExtrap, 2, function(x) qqnorm(sort(x), plot.it = F))
#   errNorm=apply(err, 2, function(x) qqnorm(sort(x), plot.it = F))
#   errNormEx=apply(errExtrap, 2, function(x) qqnorm(sort(x), plot.it = F))
#   
#   #plot-discharge
#   plot(tail(simNorm[[lt]][[2]],20),tail(simNorm[[lt]][[1]],20),
#        xlim=c(min(tail(simNorm[[lt]][[2]],20)),max(simNormEx[[lt]][[2]])),
#        ylim=c(min(tail(simNorm[[lt]][[1]],20)),max(simNormEx[[lt]][[1]])),
#        xlab="Predicted discharge [m3/s]",ylab="NQT:Predicted discharge  [-]",
#        main=paste("Sim. Extrapolation , St:",st,", LT (hrs):",lt))
#   
#   abline(a=errorModel[[5]][lt,3],b=errorModel[[5]][lt,4],col="blue")
#   points(simNormEx[[lt]][[2]],simNormEx[[lt]][[1]],col="red")
#   
#   #plot-error
#   plot(tail(errNorm[[lt]][[1]],20),tail(errNorm[[lt]][[2]],20),
#        xlim=c(min(tail(errNorm[[lt]][[1]],20)),max(errNormEx[[lt]][[1]])),
#        ylim=c(min(tail(errNorm[[lt]][[2]],20)),max(errNormEx[[lt]][[2]])),
#        xlab="NQT-Risiduals [-]",ylab="Risiduals  [m3/s]",
#        main=paste("Error Extrapolation , St:",st,", LT (hrs):",lt))
#   
#   abline(a=errorModel[[6]][lt,3],b=errorModel[[6]][lt,4],col="blue")
#   points(errNormEx[[lt]][[1]],errNormEx[[lt]][[2]],col="red")
#   
# }

