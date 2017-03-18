source("C:/Data/Thesis/4_Performance Measures/R_Scripts/Gen_util.R")
source("C:/Data/Thesis/4_Performance Measures/R_Scripts/Final/Mainfunctions_gofandCont.R")
library(hydroGOF)
library(xts)
directory="C:/Data/Thesis/4_Performance Measures/Data/"
stationFileDir="C:/Data/Thesis/4_Performance Measures/Station_Data.csv"
GOFMeasuresDir="C:/Data/Thesis/4_Performance Measures/GOF.csv"
PlotsDir="C:/Data/Thesis/4_Performance Measures/GOFplots"

GOFMeasures=read.table(GOFMeasuresDir, header=FALSE, sep=",")
stationsData = read.table(stationFileDir,header=FALSE, sep=",")
stationsID=as.character(stationsData[[1]])
#gofallst<-array(dim=c(20,))
gofall<-array(dim=c(20,24,5))
gofthr<-array(dim=c(dim(gofall),5))

A<-c(1:5)

for(j in A){
        
        #warning levels
        
        low= stationsData[j,4]
        medium= stationsData[j,5]
        yellow= stationsData[j,6]
        orange= stationsData[j,7]
        red= stationsData[j,8]
      
        Data=preprocdata(directory,stationsID[j])
        obs=Data[[1]]
        sim=Data[[2]]
        
        obs<-converttoZoo(obs)
        sim<-converttoZoo(sim)
        
        gofall[,,j]<-gofMeasures(obs,sim)
        
}

gofNoDA<-gofall
gofFOE<-gofall
gofSOE<-gofall


#plot 1- GOF of Different stations


gof_seq<-c(1,6)

for (i in gof_seq) {
        png(paste("C:/Data/Thesis/4_Performance Measures/Results/GOFplots/Stations/GOFst_",i, ".png"),width=720,height=480)
        ymin<-min(gofall[i,,],na.rm=TRUE)
        ymax<-ceiling(max(gofall[i,,],na.rm=TRUE))
        par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
        plot(1:24,gofall[i,,1], ylim=c(ymin,ymax), xlab="Lead time (hrs)",ylab=as.character(GOFMeasures[i,1]))
        points(1:24,gofall[i,,2],pch=2)
        points(1:24,gofall[i,,3],pch=3)
        points(1:24,gofall[i,,4],pch=4)
        points(1:24,gofall[i,,5],pch=5)
        labels<-c('station-1','station-2','station-3','station-4','station-5')
        legend("topright", inset=c(0,0), labels, bty='n',  pch = c(1:5) )
        dev.off()
        
}


#plot 2- GOF of Different thresholds


gofthr[,,1]<-gofall # No threshold
gofall[]<-NA

gofthr[,,2]<-gofall # medium
gofall[]<-NA

gofthr[,,3]<-gofall # yellow
gofall[]<-NA

gofthr[,,4]<-gofall # orange
gofall[]<-NA

gofthr[,,5]<-gofall # red
gofall[]<-NA


A<-c(1:5)

for(j in A){
 
  gof_seq<-c(2,4,5,6,9,10,17,18)

  for (i in gof_seq) {
        png(paste("C:/Data/Thesis/4_Performance Measures/Results/GOFplots/Threshold/GOF_",i, ".png"))
        ymin<-min(gofthr[i,,3:5], na.rm=TRUE)
        ymax<-ceiling(max(gofthr[i,,3:5], na.rm=TRUE))
        par(mar=c(5.1, 4.1, 4.1, 10.1), xpd=TRUE)
        plot(1:5,gofthr[i,,3], ylim=c(ymin,ymax), xlab="Stations", ylab=as.character(GOFMeasures[i,1]))
        points(1:5,gofthr[i,,4],pch=2)
        points(1:5,gofthr[i,,5],pch=3)
        #points(1:12,gofthr[i,1:12,j,4],pch=4)
        #points(1:12,gofthr[i,1:12,j,5],pch=5)
        labels<-c('Threshold-yellow','Threshold-orange','Threshold-red' )
        legend("topright", inset=c(-0.5,0), labels,  pch = c(1:5))
        dev.off()
        
}
}


#different setups

gof_seq<-c(2,4,5,6,9,10,17,18)

for (i in gof_seq) {
  for (j in 1:5) {
  png(paste("C:/Data/Thesis/4_Performance Measures/Results/GOFplots/ModelSetup/GOF_",j,i, ".png"))
  ymin<-min(gofNoDA[i,,j],gofFOE[i,,j],gofSOE[i,,j],gofSOEn[i,,j],na.rm=TRUE)
  ymax<-ceiling(max(gofNoDA[i,,j],gofFOE[i,,j],gofSOE[i,,j],gofSOEn[i,,j],na.rm=TRUE))
  par(mar=c(5.1, 4.1, 4.1, 10.1), xpd=TRUE)
  plot(1:12,gofNoDA[i,,j], ylim=c(ymin,ymax), xlab="Lead time (hrs)",ylab=as.character(GOFMeasures[i,1]))
  points(1:12,gofFOE[i,,j],pch=2)
  points(1:12,gofSOE[i,,j],pch=3)
  points(1:12,gofSOEn[i,,j],pch=4)
  labels<-c('NoDA','First order','Second order','Second order-New')
  legend("topright", inset=c(-0.475,0), labels,  pch = c(1:4) )
  dev.off()
  
}
}

#different setups 2



#par(mfrow=c(1,2))

gof_seq<-c(2)

for (i in gof_seq) {
  st=c(1:5)
  
  for (j in st) {
    png(paste("C:/Data/Thesis/4_Performance Measures/Results/GOFplots/ModelSetup/GOF_",j,i, ".png"), width=720,height=480)
    ymin<-min(gofNoDA[i,,j],gofFOE[i,,j],gofSOE[i,,j],na.rm=TRUE)
    ymax<-ceiling(max(gofNoDA[i,,j],gofFOE[i,,j],gofSOE[i,,j],na.rm=TRUE))
    par(mar=c(5.1, 4.1, 4.1, 2.1), xpd=TRUE)
    plot(1:24,gofNoDA[i,,j], ylim=c(0,ymax), xlim=c(1,24), xlab="Lead time (hrs)",ylab=as.character(GOFMeasures[i,1]), main=paste0('Station :',j))
    points(1:24,gofFOE[i,,j],pch=2)
    points(1:24,gofSOE[i,,j],pch=3)
    labels<-c('NoDA','First order','Second order')
    legend("bottomright", inset=c(0,0), labels,bty='n', pch = c(1:3) )
    dev.off()

  }
  
}




dev.off()

#mean values
gof_seq<-c(18)

for (i in gof_seq) {
    png(paste("C:/Data/Thesis/4_Performance Measures/Results/GOFplots/ModelSetup/GOF_",i, ".png"))
    ymin<-min(colMeans(gofNoDA[i,,]),colMeans(gofFOE[i,,]),colMeans(gofSOE[i,,]),colMeans(gofSOEn[i,,]),na.rm=TRUE)
    ymax<-ceiling(max(colMeans(gofNoDA[i,,]),colMeans(gofFOE[i,,]),colMeans(gofSOE[i,,]),colMeans(gofSOEn[i,,]), na.rm=TRUE))
    #par(mar=c(5.1, 4.1, 4.1, 10.1), xpd=TRUE)
    plot(1:5,colMeans(gofNoDA[i,,]), ylim=c(ymin,ymax), xlab="stations",ylab=as.character(GOFMeasures[i,1]))
    points(1:5,colMeans(gofFOE[i,,]),pch=2)
    points(1:5,colMeans(gofSOE[i,,]),pch=3)
    labels<-c('NoDA','First order','Second order','Second order(New)' )
    legend("topright", inset=c(0,0), labels,  pch = c(1:3) )
    dev.off()
    
  }

