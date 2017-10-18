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
  data$totalSampleDryMass <- NA
  x <- which(data$drill == 'no'); data$dtotalSampleDryMass[x] <- data$dryMass[x]
  x <- which(data$drill == 'yes'); data$totalSampleDryMass[x] <- (data$wetWeight[x] * data$dryMass[x]) / data$drilledWeight[x]
  return(data)
}

CalcDensity<-function(data){
  data$density <- NA
  data$density <- data$drilledWeight / as.numeric(data$volMass)
  return(data)
}

ReorgDataFrame<-function(data){
  
  require(tidyr)
  
  #add a species column
  data1<-separate(data, unique, into=c("Species","extraCode"), 4, remove=FALSE)
  
  #add a size column
  data1$size<-NA
  data1[tolower(data1$species) == data1$species,"size"]<-"small"
  data1[tolower(data1$species) != data1$species,"size"]<-"large"
  
  #rename
  data2<-rename(data1, "fruiting"="fruitingBodies","insects"="insectDamage")
  
  #organize columns
  data3<-data2[,c("unique","Species","size","time","totalSampleDryMass","density","fruiting","insects","drill")]
  
  return(data3)
  
}