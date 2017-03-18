#contingency table
source("C:/Data/Thesis/4_Performance Measures/R_Scripts/Gen_util.R")
library(xts)
directory="C:/Data/Thesis/4_Performance Measures/Data/"
stationFileDir="C:/Data/Thesis/4_Performance Measures/Station_Data.csv"
GOFMeasuresDir="C:/Data/Thesis/4_Performance Measures/GOF.csv"
PlotsDir="C:/Data/Thesis/4_Performance Measures/GOFplots"

Data=preprocdata(directory,stationsID[j],na.rm=TRUE)
obs=Data[[1]]
sim=Data[[2]]

obs=obs[,2:13]
sim=sim[,2:13]

a<-c(0,100,200,290)
obsmax<-apply(obs,1,max)
simmax<-apply(sim,1,max)

for (i in 1:NROW(obsmax)) {
  obsmax[i]<-max(a[which((obsmax[i]-a)>=0)])
}

for (i in 1:NROW(simmax)) {
  simmax[i]<-max(a[which((simmax[i]-a)>=0)])
}


for (i in 1:NROW(obsmax)) {
  obsmax[i]<-which((obsmax[1]-obs[1,])<=0)[1])
}





for (i in 1:NROW(obsmax)) {
  obsmax[i]<-which((obsmax[1]-obs[1,])<=0)[1]
}



b<-250

threslevel=function(Warninglevels=a,b) {
  threslevel<-max(a[which((b-a)>0)])
  return(threslevel)
}
  
  max(a[which((obsmax-a)>0)])
apply(obsmax,1,threslevel(a,))