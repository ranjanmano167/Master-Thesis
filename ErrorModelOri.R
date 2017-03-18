library(xts)
library(abind)
library(quantreg)
library(MASS)

#--------Quantile regression in original domain-------------

#Devoloping/calibrating error model

##qrType=1 ---linear QR
##qrType=2 ---non-linear QR (piecewise linear)
##Weights: 0-No weigthts, 1-Rank based weights, 2-Value based weights when type=1 is selected
##calibLowerLim: Lower limit below which the simulated values will not be considered for QR when type=1 is selected
##breakPoints: Breaking points ¨for piecewise QR as Quantiles when type=2 is selected

qrtype=1


quantileCalibOri=function (obs, sim, taus,type=qrtype, calibLowerLim=NA, weights=0,breakPoints=c(.20,.80),saveFilePath=NULL) {
  
  if (length(taus)==0){ 
    taus=c(.05, .25, 0.5, .75, .95)
  }
  
  err=obs-sim
  m=NCOL(obs)
  
  extrapUp=NULL  # Upper extrapolation linear regression coefficients
  extrapLow=NULL # Lower extrapolation linear regression coefficients
  
  
  headPoints=10 # Points used in linear regression of tail points in simulated data
  tailPoints=10 # Points used in linear regression of tail points in error
  maxiter=50
  
  breakValuesAll=NULL
  qrCoeffs=NULL
  qrCoeffsLow=NULL
  qrCoeffsHigh=NULL
  qrCoeffsMed=NULL
  
  for (i in 1:m){
    simAndErr=data.frame(sim[,i],err[,i])
    simAndErrSort=simAndErr[order(simAndErr[,1]),]
    extrapUp=rbind(extrapUp, rlm(tail(simAndErrSort[,2],tailPoints)~tail(simAndErrSort[,1],tailPoints),maxit=maxiter)[[1]])
    extrapLow=rbind(extrapLow,rlm(head(simAndErrSort[,2],headPoints)~head(simAndErrSort[,1],headPoints),maxit=maxiter)[[1]])
    
    #linear quantile regresiion 
    
    if (type==1) {
      
      #setting a threshold value to fit QR only above that value
      simfilt=sim[,i]
      simfilt[simfilt<calibLowerLim]=NA
      errfilt=err[,i]
      errfilt[is.na(simfilt)]=NA
      
      #Introducing weights 
      if (weights==2){
        #simWeights=(simfilt+abs(min(simfilt)))/(max(simfilt)-min(simfilt))
        simWeights=(simfilt+abs(min(simfilt)))/(max(simfilt)-min(simfilt))
      } else if (weights==1) {
        simWeights=rank(simfilt)/length(simfilt)
      } else if (weights==0){
        simWeights=NULL
      }
      
      
      k=coef(rq(errfilt ~  simfilt, weights=simWeights, tau = taus))
      rownames(k) = c("Intersept", "Slope")
      qrCoeffs=abind(qrCoeffs,k, along=3)
      
    } else if (type==2) {
      
      #Introducing piecewise QR
      breakValues=quantile(simAndErrSort[,1],probs=breakPoints)
      breakValuesAll=rbind(breakValuesAll,breakValues)
      simAndErrLow=simAndErrSort[simAndErrSort[,1]<breakValues[[1]],]
      simAndErrHigh=simAndErrSort[simAndErrSort[,1]>breakValues[[2]],]
      simAndErrMed=simAndErrSort[simAndErrSort[,1]>breakValues[[1]] & simAndErrSort[,1]<breakValues[[2]],]
      
      kLow=coef(rq(simAndErrLow[,2] ~  simAndErrLow[,1], tau = taus))
      rownames(kLow) = c("Intersept", "Slope")
      qrCoeffsLow=abind(qrCoeffsLow,kLow, along=3)
      
      kHigh=coef(rq(simAndErrHigh[,2] ~  simAndErrHigh[,1], tau = taus))
      rownames(kHigh) = c("Intersept", "Slope")
      qrCoeffsHigh=abind(qrCoeffsHigh,kHigh, along=3)
      
      kMed=coef(rq(simAndErrMed[,2] ~  simAndErrMed[,1], tau = taus))
      rownames(kMed) = c("Intersept", "Slope")
      qrCoeffsMed=abind(qrCoeffsMed,kMed, along=3)
    }
    
  }
  
  if (type==1) {
    results=list(taus, qrCoeffs, extrapUp, extrapLow)
  } else if (type==2) {
    results=list(taus, breakValuesAll, qrCoeffsLow,qrCoeffsHigh,qrCoeffsMed, extrapUp, extrapLow)
  }
  return(results)
}

#-------application/validation of error model------

#noncrossing: If crossing need to be avoided or not (TRUE/FALSE)
#lt-leadtime

CIEstOri=function(errorModel,type=qrtype, simVal,lt,nonCrossing=TRUE) {
  
  if (type==1) {
    
    rqFits=errorModel[[2]][,,lt]
    
    if (nonCrossing){
      m=NCOL(rqFits)
      lineCutAll=NULL
      for (i in 1:(m-1)){
        lineCut=(rqFits[1,i]-rqFits[1,(i+1)])/(rqFits[2,(i+1)]-rqFits[2,i])
        lineCutAll=cbind(lineCutAll,lineCut)
      }
      LowerLim=(max(lineCutAll[lineCutAll<200])+5)
      
      errAll=NULL
      
      for (j in 1:NCOL(rqFits)){
        if (simVal>LowerLim){ 
          errNew=(simVal*rqFits[2,j]+rqFits[1,j])
        }
        else {
          errNew=(LowerLim*rqFits[2,j]+rqFits[1,j])
        }
        
        errAll=cbind(errAll,errNew)
      }
    } else {
      
      errAll=NULL
      
      for (j in 1:NCOL(rqFits)){
        errNew=(simVal*rqFits[2,j]+rqFits[1,j])
        errAll=cbind(errAll,errNew)
      }
      
    }
    
  } else if (type==2) {
    
    breakValues=errorModel[[2]][lt,]
    rqFitsLow=errorModel[[3]][,,lt]
    rqFitsHigh=errorModel[[4]][,,lt]
    rqFitsMed=errorModel[[5]][,,lt]
    
    
    errAll=NULL
    
    for (j in 1:NCOL(rqFitsLow)){
      if (simVal<breakValues[[1]]){ 
        errNew=(simVal*rqFitsLow[2,j]+rqFitsLow[1,j])
      } else if (simVal>breakValues[[1]] & simVal<breakValues[[2]]){
        errNew=(simVal*rqFitsMed[2,j]+rqFitsMed[1,j])
      } else if (simVal>breakValues[[2]]) {
        errNew=(simVal*rqFitsHigh[2,j]+rqFitsHigh[1,j])
      }
      
      errAll=cbind(errAll,errNew)
    }
    
  }
  return(errAll) 
}

#plots

#1 Sim vs Err
plotSimVsErrOri=function(sim, err, errorModel,type=qrtype, st=NULL, lt,nonCrossing=TRUE,colour="blue")  {
  
  if(type==1){
#     ylimit=max(abs(min(err)),abs(max(err)))
#     
#     plot(sim[,lt],err[,lt],
#          main=(paste("Quantile Regression, St:",st,", LT (hrs):",lt)),
#          xlim=c(0,max(sim)),ylim=c(-ylimit,ylimit),
#          #xlim=c(600,1000),ylim=c(-100,100),
#          xlab="Predicted discharge [m3/s]",ylab="Residuals [m3/s]",
#          col="gray60")
    
    rqFits=errorModel[[2]][,,lt]
    
    if (nonCrossing){
      m=NCOL(rqFits)
      lineCutAll=NULL
      for (i in 1:(m-1)){
        lineCut=(rqFits[1,i]-rqFits[1,(i+1)])/(rqFits[2,(i+1)]-rqFits[2,i])
        lineCutAll=cbind(lineCutAll,lineCut)
      }
      LowerLim=max((max(lineCutAll[lineCutAll<200])+5),min(sim[,lt]))
      
      
      for (i in 1:m) {
        a=rqFits[1,i]
        b=rqFits[2,i]
        
        #colo=c("lightblue1","lightblue","mediumblue","navy","mediumblue","lightblue","lightblue1")
        ltyp=c(4,3,2,1,2,3,4)
        n=(7-NCOL(rqFits))/2
        
        xx=c(LowerLim,4*max(sim[,lt]))
        yy=c((b*LowerLim+a),(b*4*max(sim[,lt])+a))
        lines(xx,yy,col=colour,lty=ltyp[i+n],lw=2)
        #lines(xx,yy,col=colour,lty=1,lw=2)
        
        xxx=c(min(sim[,lt]),LowerLim)
        yyy=c((b*LowerLim+a),(b*LowerLim+a))
        lines(xxx,yyy,col=colour,lty=ltyp[i+n],lw=2)
        #lines(xxx,yyy,col=colour,lty=1,lw=2)
      }
      
    } else {
      m=NCOL(rqFits)
      
      for (i in 1:m) {
        a=rqFits[1,i]
        b=rqFits[2,i]
        
        #colo=c("lightblue1","lightblue","mediumblue","navy","mediumblue","lightblue","lightblue1")
        ltyp=c(4,3,2,1,2,3,4)
        n=(7-NCOL(rqFits))/2
        
        xx=c(min(sim[,lt]),4*max(sim[,lt]))
        yy=c((b*min(sim[,lt])+a),(b*4*max(sim[,lt])+a))
        lines(xx,yy,col=colour,lty=ltyp[i+n],lw=2)
        #lines(xx,yy,col=colour,lty=2,lw=2)
      }
    }
    
  } else if (type==2){
    
    breakValues=errorModel[[2]][lt,]
    rqFitsLow=errorModel[[3]][,,lt]
    rqFitsHigh=errorModel[[4]][,,lt]
    rqFitsMed=errorModel[[5]][,,lt]
    
    ylim=c(min(err),max(err))
    plot(sim[,lt],err[,lt],
         main=(paste("Scatter plot, St:",st,", LT (hrs):",lt)),
         #xlim=c(0,max(sim)),ylim=c(min(err),max(err)),
         #xlim=c(0,200),ylim=c(-200,200),
         xlim=c(30,80),ylim=c(-25,25),
         xlab="Predicted discharge [m3/s]",ylab="Residuals [m3/s]",
         col="gray60")
    
    for (i in 1:NCOL(rqFitsLow)) {
      aLow=rqFitsLow[1,i]
      bLow=rqFitsLow[2,i]
      
      aHigh=rqFitsHigh[1,i]
      bHigh=rqFitsHigh[2,i]
      
      aMed=rqFitsMed[1,i]
      bMed=rqFitsMed[2,i]
      
      #colo=c("lightblue1","lightblue","mediumblue","navy","mediumblue","lightblue","lightblue1")
      ltyp=c(4,3,2,1,2,3,4)
      n=(7-NCOL(rqFitsLow))/2
      
      xx=c(min(sim[,lt]), breakValues[[1]])
      yy=c((bLow*min(sim[,lt])+aLow),(bLow*breakValues[[1]]+aLow))
      lines(xx,yy,col=colour,lty=ltyp[i+n],lw=2)
      #lines(xx,yy,col=colour,lty=1,lw=2)
      
      xxx=c(breakValues[[1]],breakValues[[2]])
      yyy=c((bMed*breakValues[[1]]+aMed),(bMed*breakValues[[2]]+aMed))
      lines(xxx,yyy,col=colour,lty=ltyp[i+n],lw=2)
      #lines(xxx,yyy,col=colour,lty=1,lw=2)
      
      xxxx=c(breakValues[[2]],max(4*sim[,lt]))
      yyyy=c((bHigh*breakValues[[2]]+aHigh),(bHigh*4*max(sim[,lt])+aHigh))
      lines(xxxx,yyyy,col=colour,lty=ltyp[i+n],lw=2)
      #lines(xxxx,yyyy,col=colour,lty=1,lw=2)
      
    }
    
  }
  
}


#2.1 CI vs Lead time
plotCIOri=function(errorModel,type=qrtype,simNew,obsNew=NULL, st=NULL) {
  
  m=length(simNew)
  
  if (type==1){
    n=NCOL(errorModel[[2]])
  } else if (type==2){
    n=NCOL(errorModel[[3]])
  }
  
  
  CI=NULL
  for (i in 1:m){
    CILt=CIEstOri(errorModel,type=qrtype,simVal=simNew[i], lt=i)
    CI=rbind(CI,CILt)
  }
  
  a=(CI[,(n-1)/2]+as.numeric(simNew))
  b=(CI[,(n+3)/2]+as.numeric(simNew))
  c=(CI[,(n-3)/2]+as.numeric(simNew))
  d=(CI[,(n+5)/2]+as.numeric(simNew))
  e=as.numeric(obsNew)
  ymin=0 #min(a,b,c,d,e)
  ymax=max(2*max(simNew),max(obsNew))#max(a,b,c,d,e)
  
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
plotCIltOri=function(errorModel,type=qrtype,simAll,obsAll=NULL, st=NULL,lt) {
  
  simNew=simAll[,2:25]
  obsNew=obsAll[,2:25]
  
  m=NROW(simNew)
  
  if (type==1){
    n=NCOL(errorModel[[2]])
  } else if (type==2){
    n=NCOL(errorModel[[3]])
  }
  
  CI=NULL
  for (i in 1:m){
    CIAll=CIEstOri(errorModel,type=qrtype,simVal=simNew[i,lt],lt=lt)
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
  grid(nx=NA,ny=NULL)
#   polygon(c(xVal,rev(xVal)),c(a,rev(b)),col='red',border = NA)
#   polygon(c(xVal,rev(xVal)),c(a,rev(c)),col='hotpink',border = NA)
#   polygon(c(xVal,rev(xVal)),c(b,rev(d)),col='hotpink',border = NA)

  polygon(c(xVal,rev(xVal)),c(a,rev(b)),col='lightblue',border = NA)
  polygon(c(xVal,rev(xVal)),c(a,rev(c)),col='lightblue',border = NA)
  polygon(c(xVal,rev(xVal)),c(b,rev(d)),col='lightblue',border = NA)


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
plotCIltLinesOri=function(errorModel,type=qrtype,simAll,obsAll=NULL, st=NULL,lt) {
  
  simNew=simAll[,2:25]
  obsNew=obsAll[,2:25]
  
  m=NROW(simNew)
  
  if (type==1){
    n=NCOL(errorModel[[2]])
  } else if (type==2){
    n=NCOL(errorModel[[3]])
  }
  
  CI=NULL
  for (i in 1:m){
    CIAll=CIEstOri(errorModel,type=qrtype,simVal=simNew[i,lt],lt=lt)
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
#   points(xVal,a,col='seagreen4',type="l",lwd=2)
#   points(xVal,b,col='seagreen4',type="l",lwd=2)
  points(xVal,c,col='seagreen4',type="l",lwd=2)
  points(xVal,d,col='seagreen4',type="l",lwd=2)
#   points(xVal,simNew[,lt],col='red3',type='l')
#   
  #for (i in 1:n) {
  #  x=(CI[,i]+as.numeric(simNew))
  #  points(1:m,x,type='l')
  #}
  
#   if(is.null(obsNew)=='FALSE'){
#     points(xVal,e,type='p',col='red',pch=16)
#   }
#   
}


#3. Plots for  quantile slope and gradient

plotQuanRegFitOri=function(errorModel,type=qrtype,taus, lt) {
  if (type==1) {
    rqFits=errorModel[[2]][,,lt]
  } else if (type==2) {
    rqFits=errorModel[[4]][,,lt]
  }
  plot(taus,rqFits[2,],type='b',xlim=c(0,1),ylim=c(-1,1),xlab='Quantile',ylab='Slope',main=paste('Quantile regression fit, LT (hrs):',lt))
  plot(taus,rqFits[1,],type='b',xlim=c(0,1),ylim=c(min(rqFits[1,]),max(rqFits[1,])),xlab='Quantile',ylab='Intersection')
}






