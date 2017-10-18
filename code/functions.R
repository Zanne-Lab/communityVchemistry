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