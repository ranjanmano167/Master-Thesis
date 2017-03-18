#To calculate gof measures for a time series with different lead times
#obs=observed time series with lead time increasing along rows
#sim=simulated time series with lead time increasing along rows
#lthres=lower threshold to be considered, observation (and corrosponding simulation) lower 
#         than this value will be removed (default=0)
#uthres=upper threshold to be considered, observation (and corrosponding simulation) higher 
#         than this value will be removed (default=max(obs,sim))

library(hydroGOF)

gofMeasures<-function (obs=obs,sim=sim,m=NCOL(obs),lthres=0,uthres=max(obs,sim,na.rm=TRUE)){
  
  gofLeadtime<-matrix(nrow=20,ncol=m)
  
  for(i in 1:m) {
    
    obsf<-obs[,i]     
    simf<-sim[,i]
    
    # to set threshold values
    
    obsf[obsf>uthres]<-NA
    obsf[obsf<lthres]<-NA
    
    
    #to skip columns with all N/A values
    
    sum<-obsf+simf
    if (mean(sum,na.rm=TRUE)=='NaN'){
      gofLeadtime[,i]<-NA
    } else if ((sum(sum[,1],na.rm=TRUE))==(mean(sum[,1],na.rm=TRUE))){
      gofLeadtime[,i]<-NA
    } else {  
      gofLeadtime[,i]<-gof(simf,obsf,na.rm=TRUE)
    }
    
  }
    return(gofLeadtime)  
}


