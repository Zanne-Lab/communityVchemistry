
# load mass loss data and calculate percent mass remaining (pmr)

process_initial_large<-function(df,size){
  
  if ("Dry mass total (g)" %in% names(df)) {
    df <- rename(df,`Dry mass (g)`=`Dry mass total (g)`)
  }
  
  df %>%
    mutate(dry_mass_content=`Dry mass (g)`/`Fresh mass (g)`) %>%
    filter(!is.na(dry_mass_content)) %>%
    group_by(Species) %>%
    summarize(dry_mass_prop=mean(dry_mass_content,na.rm=T),n()) -> moisture
  
  df %>%
    left_join(moisture) %>%
    mutate(totalSampleDryMass=`Fresh mass (g)`*dry_mass_prop,size=size,density=NA,time=0,fruiting=NA,insects=NA,drill=NA) %>%
    select(unique, Species, size,time,totalSampleDryMass,density,fruiting,insects,drill) %>%
    rename("species"="Species") -> df_out
  
  
  #TODO ADD A REAL CALCULATION OF WOOD DENSITY ONCE WE UNDERSTAND HOW TO DO THAT
  
  return(df_out)
}

process_initial_small<-function(df,size){
  
  if ("Dry mass total (g)" %in% names(df)) {
    df <- rename(df,`Dry mass (g)`=`Dry mass total (g)`)
  }
  
  length <- 10 # LENGTH OF STEMS IN CM, SHOULD PROBABLY BE HIGHER UP SOMEWHERE, IF NEEDED ELSEWHERE
  
  # CALCULATE INITIAL PROPORTION DRY MASS
  df %>%
    mutate(dry_mass_content=ifelse(
      !is.na(`Dry mass (g)`),
      `Dry mass (g)`/`Fresh mass (g)`,
      (`Dry mass wood (g)`+`Dry mass bark (g)`)/`Fresh mass (g)`
    )) %>%
    filter(!is.na(dry_mass_content)) %>%
    group_by(Species) %>%
    summarize(dry_mass_prop=mean(dry_mass_content,na.rm=T),n()) -> moisture
  
  # WOOD DENSITY, XYLEM DENSITY, AND BARK DENSITY
  df %>%
    mutate(xylem_density=`Dry mass wood (g)`/`Volume (g)`) %>%
    mutate(total_volume=(pi*(`Diameter.wbark (mm)`/20)^2*length)) %>% # converting to cm here
    mutate(total_density=ifelse(
      !is.na(`Dry mass (g)`),
      `Dry mass (g)`/total_volume,
      (`Dry mass wood (g)`+`Dry mass bark (g)`)/total_volume
    )) %>%
    mutate(bark_volume=(pi*(`Diameter.wbark (mm)`/20)^2*length)-(pi*(`Diameter.nobark (mm)`/20)^2*length)) %>% # convert to cm
    mutate(bark_density=`Dry mass bark (g)`/bark_volume)->df
  
  df %>%
    left_join(moisture) %>%
    mutate(totalSampleDryMass=`Fresh mass (g)`*dry_mass_prop,size=size,time=0,fruiting=NA,insects=NA,drill=NA) %>%
    select(unique, Species, size,time,totalSampleDryMass,bark_density,xylem_density,total_density,fruiting,insects,drill) %>%
    rename("species"="Species") -> df_out
  
  return(df_out)
}

read_in_initial_mass <- function(){
  require(readr)
  require(dplyr)
  big <- read_csv("data/covariates_bigStems.csv")
  small <- read_csv("data/covariates_smallStems.csv")
  
  #look for missing data in olst small
  filter(small, Species == 'olst') -> tmp
  #View(tmp) #no samples with `Dry mass total (g)`
  
  big_out <- process_initial_large(big,"large")
  small_out <- process_initial_small(small,"small")
  
  df_out<-bind_rows(big_out,small_out)
  return(df_out)
} 

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
  # wet weight excess measured separately on multiple pieces
  temp <- strsplit(samp3$wetWeightExcess, '+', fixed=T)
  samp3$wetWeightExcess <- sapply(temp, function(x){if(length(x) == 2) sum(as.numeric(x)) else if(length(x) == 1) as.numeric(x) else NA})
  rm(temp)
  # when recording wood volume for small stems broken into two, could not measure on whole piece so additional wet weights recorded for relevant fragment
  # following code is for sorting this out
  temp <- strsplit(samp3$volMass, 'gpiece=', fixed=T)
  samp3$volMass <- sapply(temp, function(x) {if(length(x) == 2) x[2] else x[1]})
  samp3$volMass <- as.numeric(samp3$volMass)
  x <- is.na(samp3$weightForVol)
  samp3$weightForVol[x] <- samp3$dryMass[x]
  rm(temp, x)
  # include excess wood in bag (fragments falling off during transit) for wet mass of whole harvested piece
  # dry mass weighed for excess and nonexcess all at once
  samp3$wetWeightForMass <- apply(samp3[, c('wetWeight', 'wetWeightExcess')], 1, sum, na.rm=T)
  samp3[is.na(samp3$wetWeight), 'wetWeightForMass'] <- NA
  samp3$time <- 25
  return(samp3)
}

read.samp4 <- function(){
  samp4 <- read.csv('data/samp4data_201608_quantitative.csv', stringsAsFactors=F)
  # ensure consistency in column names
  names(samp4) <- gsub('wetWeightExcess..g.', 'wetWeightExcess', names(samp4))
  # when recording wood volume, could not measure on whole piece so additional wet weights recorded for relevant fragment
  # following code is for sorting this out
  temp <- strsplit(gsub('^\\(', '', samp4$drilledWeight), ' total) ', fixed=T)
  # samp4$weightForVol <- sapply(temp, function(x){if(length(x) == 2) as.numeric(x)[2] else as.numeric(x)[1]})
  # x <- samp4$weightForVol %in% 0
  # samp4[x, 'weightForVol'] <- samp4[x, 'wetWeight']
  samp4$weightForVol <- samp4$dryMass.piece.used.to.do.vol.mass.
  samp4$drilledWeight <- sapply(temp, function(x){as.numeric(x)[1]})
  samp4[samp4$drill == 'no', 'drilledWeight'] <- NA
  samp4[samp4$drill == 'yes' & samp4$wetWeight == 0, 'drilledWeight'] <- NA
  rm(temp)
  # include excess wood in bag (fragments falling off during transit) for wet and dry mass of whole harvested piece
  samp4$wetWeightForMass <- apply(samp4[, c('wetWeight', 'wetWeightExcess')], 1, sum, na.rm=T)
  samp4[is.na(samp4$wetWeight), 'wetWeightForMass'] <- NA  # two samples missing from plot
  samp4[is.na(samp4$total.dry), 'total.dry'] <- samp4[is.na(samp4$total.dry), 'dryMass.piece.used.to.do.vol.mass.']
  samp4$dryMass <- apply(samp4[, c('dry.WWE', 'total.dry')], 1, sum, na.rm=T)
  samp4[with(samp4, which(is.na(dry.WWE) & is.na(total.dry))), 'dryMass'] <- NA
  # include damage scoring data
  samp4.1 <- read.csv('data/samp4data_201608_qualitative.csv', stringsAsFactor=F)
  names(samp4.1) <- gsub('notes', 'notes1', names(samp4.1))
  samp4 <- merge(samp4, samp4.1)
  samp4$notes <- with(samp4, paste(notes, notes1, sep=' -- '))
  samp4$notes1 <- samp4$dry.WWE <- samp4$dryMass.piece.used.to.do.vol.mass. <- NULL
  samp4$time <- 37
  #select columns to keep
  keepCols<-c("order","unique","drill","wetWeight","fruitingBodies","wetWeightExcess",
              "drilledWeight","volMass","volMassRetained","insectDamage","weightForVol","dryMass",
              "notes","typesInsects","wetWeightForMass","time")
  samp4<-samp4[, keepCols]
  return(samp4)
}

read.samp5 <- function(){
  samp5 <- read.csv('data/samp5data_201806_quantitative.csv', stringsAsFactors=F)
  # ensure consistency in column names
  names(samp5) <- gsub('dryMassVolPiece', 'weightForVol', names(samp5))
  names(samp5) <- gsub('totalDryMass', 'total.dry', names(samp5))
  # 'totalDryMass' only included value if 'dryMassVolPiece' measured on fragment, sum of fragment and other wood pieces
  # include mass of fragment where this column is NA
  samp5[is.na(samp5$total.dry), 'total.dry'] <- samp5[is.na(samp5$total.dry), 'weightForVol']
  # include excess wood in bag (fragments falling off during transit) for wet and dry mass of whole harvested piece
  samp5$wetWeightForMass <- apply(samp5[, c('wetWeight', 'wetWeightExcess')], 1, sum, na.rm=T)
  samp5[is.na(samp5$wetWeight), 'wetWeightForMass'] <- NA
  samp5$dryMass <- apply(samp5[, c('dry.WWE', 'total.dry')], 1, sum, na.rm=T)
  samp5[with(samp5, which(is.na(dry.WWE) & is.na(total.dry))), 'dryMass'] <- NA
  # include damage scoring data
  samp5.1 <- read.csv('data/samp5data_201806_qualitative.csv', stringsAsFactor=F)
  names(samp5.1) <- gsub('notes', 'notes1', names(samp5.1))
  names(samp5.1) <- gsub('typeInsects', 'typesInsects', names(samp5.1))
  samp5 <- merge(samp5, samp5.1)
  samp5$notes <- with(samp5, paste(notes, notes1, sep=' -- '))
  samp5$notes1 <- samp5$dry.WWE <- samp5$dryMass.piece.used.to.do.vol.mass. <- NULL
  samp5$time <- 59
  #select columns to keep
  keepCols<-c("order","unique","drill","wetWeight","fruitingBodies","wetWeightExcess",
              "drilledWeight","volMass","volMassRetained","insectDamage","weightForVol","dryMass",
              "notes","typesInsects","wetWeightForMass","time")
  samp5<-samp5[, keepCols]
  samp5$wetWeight <- samp5$wetWeight / 1000
  samp5$wetWeightExcess <- samp5$wetWeightExcess / 1000
  samp5$drilledWeight <- samp5$drilledWeight / 1000
  samp5$volMass <- samp5$volMass / 1000
  samp5$volMassRetained <- samp5$volMassRetained / 1000
  samp5$weightForVol <- samp5$weightForVol / 1000
  samp5$dryMass <- samp5$dryMass / 1000
  samp5$wetWeightForMass <- samp5$wetWeightForMass / 1000
  return(samp5)
}

harvest_CalcTotalDryMass<-function(data){
  data$totalSampleDryMass <- NA
  x <- which(data$drill == 'no' | data$dryMass == 0 | data$weightForVol == 0); data$totalSampleDryMass[x] <- data$dryMass[x]
  x <- which(data$drill == 'yes' & !is.na(data$drilledWeight) & data$wetWeightForMass > 0); data$totalSampleDryMass[x] <- (data$drilledWeight[x] * data$dryMass[x]) / data$wetWeightForMass[x]
  x <- which(data$drill == 'yes' & is.na(data$drilledWeight) & data$wetWeightForMass > 0); data$totalSampleDryMass[x] <- data$dryMass[x]
  data$totalSampleDryMass <- round(data$totalSampleDryMass, 2)
  return(data)
}

harvest_CalcDensity<-function(data){
  data$total_density <- NA
  data$total_density <- round(data$weightForVol / data$volMass, 2)
  return(data)
}

harvest_ReorgDataFrame<-function(data){
  require(tidyr)
  
  #add a species column
  data<-separate(data, unique, into=c("species","extraCode"), 4, remove=FALSE)
  
  #add a size column
  data$size<-"large"
  data[tolower(data$species) == data$species,"size"]<-"small"
  
  #rename and select columns
  data %>%
    rename("fruiting"="fruitingBodies","insects"="insectDamage") %>%
    select(order, unique, species, size, time, totalSampleDryMass, total_density, fruiting, insects, drill, notes) -> data
  
  return(data)
  
}

LoadHarvestFiles<-function(){
  
  # load data
  s1 <- read.samp1()
  s2 <- read.samp2()
  s3 <- read.samp3()
  s4 <- read.samp4()
  s5 <- read.samp5()
  
  #bind everything together
  s.data<-rbind(s1,s2,s3,s4,s5)
  
  #calculate total dry mass
  s.data<-harvest_CalcTotalDryMass(s.data)
  
  #calculate density
  s.data<-harvest_CalcDensity(s.data)
  
  #reorganize data frame
  s.data<-harvest_ReorgDataFrame(s.data)
  
  #fix issues with problematic samples
  #duplicates -- s.data[duplicated(s.data$unique, fromLast=F), ]; s.data[duplicated(s.data$unique, fromLast=T), ]
  s.data[s.data$unique == 'ripi3k' & s.data$time == 37, c('order', 'unique', 'notes')] <- c(NA, 'ripi3b', '')  # ripi3b taken instead of ripi3k at 37
  s.data <- filter(s.data, !(unique == 'ALLI311' & time == 37))  # not found at 37, harvested later
  s.data <- filter(s.data, !(unique == 'baae1a' & time == 37))  # not found at 37, harvested later
  #very high density values -- filter(s.data, total_density > 1)
  
  
  #check for missing data
  #filter(s.data, is.na(totalSampleDryMass))
  filter(s.data, notes =="all wwe -- all wet weight excess")  # sample not intact so no drilling and density not measured
  filter(s.data, is.na(totalSampleDryMass), notes !="all wwe -- all wet weight excess") #for missing samples, we don't know if they rotted away completely or were moved/missing, so entered as NA
  #some looked to be completely rotted, entered totalSampleDryMass as 0
  x <- c(grep('rotted', s.data$notes, value=T), grep('crumbled', s.data$notes, value=T))
  s.data[s.data$notes %in% x & is.na(s.data$totalSampleDryMass), 'totalSampleDryMass'] <- 0
  rm(x)
  
  return(s.data)
  
}

Calc_massRemaining<-function(initial_mass, harvest_mass){
  
  mass.data<-bind_rows(initial_mass, harvest_mass)
  
  mass.data %>%
    filter(time==0) %>%
    rename(timeZeroDensity=density) %>%
    rename(timeZeroMass=totalSampleDryMass) %>%
    select(unique,timeZeroMass,timeZeroDensity)->time_zero
  
  mass.data %>%
    left_join(time_zero,by="unique") %>%
    mutate(pmr=totalSampleDryMass/timeZeroMass) %>%
    mutate(SpeciesCode=tolower(species)) -> plotting_df
  
  plotting_df %>% 
    filter(!is.na(pmr)) %>%
    select(unique, species, SpeciesCode, size, time, pmr) %>%
    rename(code=species) %>%
    rename(species=SpeciesCode) -> tmp
  
  #add codeStem based on deployment lookup table that links codeStem to unique
  deployment <- read_csv("data/deployment.csv")
  deployment %>%
    rename("code"="species") %>%
    mutate(codeStem=paste0(code, Stem)) %>%
    select(codeStem, unique) -> deploy.indx
  tmp %>%
    left_join(deploy.indx) -> pause
  
  # can the Stem values be double digits?
  # deployment %>%
  #   rename("code"="species") %>%
  #   mutate(codeStem=paste0(code, Stem)) %>%
  #   filter(Stem > 9)
  # yes, but only leer have Stem == 10 in this dataset
  
  #some of these unique ids are not in the deployment lookup table because they are time0 samples
  pause %>%
    filter(is.na(codeStem)) %>%
    separate(unique, into = c("codePart","Stem","leftover"), sep = c(4,5), remove = F) %>%
    mutate(codeStem = paste0(code, Stem)) -> test
  #sum(test$codePart != test$code) #all match
  #look for potential stem '10's that got split
  #test %>% filter(Stem == "1") -> test.sub
  #unique(test.sub$leftover) #there aren't any leading 0s, so it looks like there weren't any hidden 10's that got split
  t0.indx <- test[,c("codeStem","unique")]
  UNIQ <- unique(test$unique)
  for(i in 1:length(UNIQ)){
    fill <- t0.indx[t0.indx$unique == UNIQ[i],"codeStem"]
    pause[pause$unique == UNIQ[i],"codeStem"] <- fill$codeStem
  }
  sum(is.na(pause$codeStem)) # if 0, then all samples have an assigned codeStem
  return_df <- pause
  
  return(return_df)
}


# average pmr for each stem and timepoint

AvgPMR_byStem<-function(plotting_df){
  
  #average pmr by codeStem and time
  plotting_df %>%
    select(codeStem, time, pmr) %>%
    group_by(codeStem, time) %>%
    summarize(mean.pmr=mean(pmr, na.rm=TRUE),
              sd.pmr=sd(pmr, na.rm=TRUE),
              n.pmr=length(pmr)) -> pmr.byStem.df
  
  #check out the number of stemSamples per timestep here before I remove that information again...
  
  #make wide amd add back code, species, size
  indx <- unique(plotting_df[,c("codeStem","code","species","size")])
  pmr.byStem.df %>%
    select(-c(sd.pmr, n.pmr)) %>%
    spread(key=time, value=mean.pmr) %>%
    left_join(indx) %>%
    rename('time0'=`0`,
           'time7'=`7`,
           'time13'=`13`,
           'time25'=`25`,
           'time37'=`37`, 
           'time59'=`59`) -> pmr.byStem.df.w
  
  return(pmr.byStem.df.w)
  
}

n.PMR_byStem<-function(plotting_df){
  
  indx <- unique(plotting_df[,c("codeStem","code","species","size")])
  
  #number of pmr samples by codeStem and time
  plotting_df %>%
    select(codeStem, time, pmr) %>%
    group_by(codeStem, time) %>%
    summarize(n.pmr=length(pmr)) %>%
    spread(key=time, value=n.pmr) %>%
    left_join(indx) %>%
    rename('time0'=`0`,
           'time7'=`7`,
           'time13'=`13`,
           'time25'=`25`,
           'time37'=`37`) -> n.pmr.byStem.df.w
  
  return(n.pmr.byStem.df.w)
  
}


# calculate decay trajectory fits for each species+size

Calc_R2<-function(ne_fits_df){
  
  pred<-ne_fits_df[['predicted']]
  mass<-ne_fits_df[['mass']]
  
  sstot<-sum((mass - mean(mass))^2)
  ssres<-sum((pred - mass)^2)
  r2<-1-(ssres/sstot)
  
  return(r2)
}

fit_all_curves<-function(df_in, stemSamples){
  
  #negative expon fit
  ne_fits <- lapply(split(df_in, factor(df_in$code)),function(x){
    fit_litter(time = x$time/12, 
               mass.remaining = x$pmr, model = c("neg.exp"), iters = 500)
  })
  
  #weibull fit
  w.fits <- lapply(split(df_in, factor(df_in$code)), function(x){
    fit_litter(time = x$time/12, 
               mass.remaining = x$pmr, model = c("weibull"), iters = 500)
  })
  
  #create bootstrap distrib of k and get the 95% CI
  # bootk<-lapply(ne_fits, function(x){
  #   
  #   bootmat<-bootstrap_parameters(x, nboot=500) #1st column are iteration of the param, right, but what is in the second col?
  #   meanboot<-mean(bootmat[,1])
  #   sdboot<-sd(bootmat[,1])
  #   seboot<-sdboot/sqrt(length(bootmat[,1]))
  #   upperboot<-meanboot+(1.96*seboot)
  #   lowerboot<-meanboot-(1.96*seboot)
  #   
  #   result<-list(upper=upperboot, lower=lowerboot)
  #   return(result)
  # })
  # tmp<-lapply(bootk, function(x) data.frame(upper=x[[1]], lower=x[[2]]))
  # k.upper<-unlist(lapply(tmp, function(x) x[['upper']]))
  # k.lower<-unlist(lapply(tmp, function(x) x[['lower']]))
  
  #create a 95% CI bootstrap distribution of t70
  #not sure how to do this...
  
  #put everything together in a df
  k<-unlist(lapply(ne_fits, function(x) x$optimFit$par))
  t70<-unlist(lapply(ne_fits, function(x) time_to_prop_mass_remaining(x,threshold.mass=0.70)))
  neg.exp.aic<-unlist(lapply(ne_fits, function(x) x$fitAICc))
  ne.r2<-unlist(lapply(ne_fits, function(x) Calc_R2(x)))
  w.t70<-unlist(lapply(w.fits, function(x) time_to_prop_mass_remaining(x,threshold.mass=0.70)))
  w.aic<-unlist(lapply(w.fits, function(x) x$fitAICc))
  alpha<-unlist(lapply(w.fits, function(x) x$optimFit$par[[2]]))
  w.r2<-unlist(lapply(w.fits, function(x) Calc_R2(x)))
  spdf<-data.frame(k=k,
                   #k.upper=k.upper,
                   #k.lower=k.lower,
                   t70=t70,
                   w.t70 = w.t70,
                   ne.r2 = ne.r2,
                   ne.aic = neg.exp.aic,
                   w.aic = w.aic,
                   alpha = alpha)
  spdf$code<-rownames(spdf)
  
  #annotate df with species, size
  stemSamples %>%
    select(code, species, size) -> indx
  indx <- unique(indx)
  spdf %>%
    left_join(indx) -> spdf
  
  return(spdf)
}


# compare estimated t70 based on the negative exponential vs weibull model

comparePlot_ne_weibull <- function(decayfits){
  
  p <- ggplot(decayfits, aes(x=t70,y=w.t70, col=size))+
   geom_point()+
   labs(x="Time to 30% mass loss (negative exponential)",
        y="Time to 30% mass loss (Weibull)")+
   geom_abline(slope=1,intercept=0,linetype="dashed") + theme_bw()

  return(p)

}


