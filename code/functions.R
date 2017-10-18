#functions to calculate sample dry mass and density

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