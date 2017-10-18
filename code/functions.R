#functions to calculate sample dry mass and density

LoadHarvestFiles<-function(){
  
  s1<-read.csv("data/samp1data_201402.csv")
  s2<-read.csv("data/samp2data_201408.csv")
  s3<-read.csv("data/samp3data_201508.csv")
  
  #add an extra column to make all parallel
  s1$notes<-NA
  
  #add column for time
  s1$time<-6
  s2$time<-12
  s3$time<-24
  
  #bind everything together
  s.data<-rbind(s1,s2,s3)
  
  return(s.data)
  
}


CalcTotalDryMass<-function(data){
  data$dryMassCalc <- NA
  x <- which(data$drill == 'no'); data$dryMassCalc[x] <- data$dryMass[x]
  x <- which(data$drill == 'yes'); data$dryMassCalc[x] <- (data$wetWeight[x] * data$dryMass[x]) / data$drilledWeight[x]
  return(data)
}

CalcDensity<-function(data){
  data$densityCalc[x] <- data$drilledWright[x] / data$volMass[x]
  return(data)
}

ReorgDataFrame<-function(data){
  
  require(tidyr)
  
  #add a species column
  data1<-separate(data, unique, into=c("species","extraCode"), 4, remove=FALSE)
  
  #add a size column
  data1$size<-NA
  data1[tolower(data1$species) == data1$species,"size"]<-"small"
  data1[tolower(data1$species) != data1$species,"size"]<-"large"
  
  #organize columns
  data2<-data1[,c("unique","species","size","time","totalSampleDryMass","density","fruitingBodies","insectDamage","drill")]
  
  return(data2)
  
}