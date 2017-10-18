#functions to calculate sample dry mass and density

LoadHarvestFiles<-function(){
  
  # functions for data cleaning
  read.samp1 <- function(){
    samp1 <- read.csv('data/samp1data_201402.csv', stringsAsFactors=F)
    samp1$notes <- ''
    samp1$typesInsects <- ''
    samp1$weightForVol <- samp1$dryMass
    samp1$wetWeightForMass <- apply(samp1[, c('wetWeight', 'wetWeightExcess')], 1, sum, na.rm=T)
    samp1[is.na(samp1$wetWeight), 'wetWeightForMass'] <- NA
    samp1$time <- 7
    return(samp1)
  }
  
  read.samp2 <- function(){
    samp2 <- read.csv('data/samp2data_201408.csv', stringsAsFactors=F)
    samp2$typesInsects <- ''
    samp2$weightForVol <- samp2$dryMass
    samp2$wetWeightForMass <- apply(samp2[, c('wetWeight', 'wetWeightExcess')], 1, sum, na.rm=T)
    samp2[is.na(samp2$wetWeight), 'wetWeightForMass'] <- NA
    samp2$time <- 13
    return(samp2)
  }
  
  read.samp3 <- function(){
    samp3 <- read.csv('data/samp3data_201508.csv', stringsAsFactors=F)
    samp3$typesInsects <- ''
    temp <- strsplit(samp3$wetWeightExcess, '+', fixed=T)
    samp3$wetWeightExcess <- sapply(temp, function(x){if(length(x) == 2) sum(as.numeric(x)) else if(length(x) == 1) as.numeric(x) else NA})
    samp3$weightForVol <- samp3$dryMass
    samp3$wetWeightForMass <- apply(samp3[, c('wetWeight', 'wetWeightExcess')], 1, sum, na.rm=T)
    samp3[is.na(samp3$wetWeight), 'wetWeightForMass'] <- NA
    samp3$time <- 25
    samp3$volMass <- as.numeric(samp3$volMass) # this needs to be removed after fixing this column for time 3
    return(samp3)
  }
  
  read.samp4 <- function(){
    samp4 <- read.csv('data/samp4data_201608_quantitative.csv', stringsAsFactors=F)
    names(samp4) <- gsub('wetWeightExcess..g.', 'wetWeightExcess', names(samp4))
    temp <- strsplit(gsub('^\\(', '', samp4$drilledWeight), ' total) ', fixed=T)
    samp4$weightForVol <- sapply(temp, function(x){if(length(x) == 2) as.numeric(x)[2] else as.numeric(x)[1]})
    samp4[samp4$weightForVol == 0, 'weightForVol'] <- samp4[samp4$weightForVol == 0, 'wetWeight']
    samp4$wetWeightForMass <- apply(samp4[, c('wetWeight', 'wetWeightExcess')], 1, sum, na.rm=T)
    samp4[is.na(samp4$wetWeight), 'wetWeightForMass'] <- NA
    samp4[is.na(samp4$total.dry), 'total.dry'] <- samp4[is.na(samp4$total.dry), 'dryMass.piece.used.to.do.vol.mass.']
    samp4$dryMass <- apply(samp4[, c('dry.WWE', 'total.dry')], 1, sum)
    samp4.1 <- read.csv('data/samp4data_201608_qualitative.csv', stringsAsFactor=F)
    names(samp4.1) <- gsub('notes', 'notes1', names(samp4.1))
    samp4 <- merge(samp4, samp4.1)
    samp4$notes <- with(samp4, paste(notes, notes1, sep=' -- '))
    samp4$notes1 <- samp4$dry.WWE <- samp4$dryMass.piece.used.to.do.vol.mass. <- NULL
    samp4$time <- 37
    return(samp4[, names(s3)])
  }
  
  # load data
  s1 <- read.samp1()
  s2 <- read.samp2()
  s3 <- read.samp3()
  s4 <- read.samp4()
  
  #bind everything together
  s.data<-rbind(s1,s2,s3,s4)
  
  return(s.data)
  
}


CalcTotalDryMass<-function(data){
  data$totalSampleDryMass <- NA
  x <- which(data$drill == 'no'); data$totalSampleDryMass[x] <- data$dryMass[x]
  x <- which(data$drill == 'yes'); data$totalSampleDryMass[x] <- (data$weightForVol[x] * data$dryMass[x]) / data$wetWeightForMass[x]
  return(data)
}

CalcDensity<-function(data){
  data$density <- NA
  data$density <- data$weightForVol / as.numeric(data$volMass)
  return(data)
}

ReorgDataFrame<-function(data){
  
  require(tidyr)
  
  #add a species column
  data1<-separate(data, unique, into=c("Species","extraCode"), 4, remove=FALSE)
  
  #add a size column
  data1$size<-NA
  data1[tolower(data1$Species) == data1$Species,"size"]<-"small"
  data1[tolower(data1$Species) != data1$Species,"size"]<-"large"
  
  #rename
  data2<-rename(data1, "fruiting"="fruitingBodies","insects"="insectDamage")
  
  #organize columns
  data3<-data2[,c("unique","Species","size","time","totalSampleDryMass","density","fruiting","insects","drill")]
  
  return(data3)
  
}
