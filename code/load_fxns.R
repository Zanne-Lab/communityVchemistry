#load data and clean initial dataframes

# Harvest data and initial mass

read_in_initial_mass <- function(){
  library(readr)
  library(dplyr)
  big <- read_csv("data/covariates_bigStems.csv")
  small <- read_csv("data/covariates_smallStems.csv")
  
  big_out <- process_initial_file(big,"large")
  small_out <- process_initial_file(small,"small")
  
  df_out<-bind_rows(big_out,small_out)
  return(df_out)
} 

process_initial_file<-function(df,size){
  
  if ("Dry mass total (g)" %in% names(df)) {
    df <- rename(df,`Dry mass (g)`=`Dry mass total (g)`)
  }
  
  df %>%
    mutate(dry_mass_content=`Dry mass (g)`/`Fresh mass (g)`) %>%
    filter(!is.na(dry_mass_content)) %>%
    group_by(Species) %>%
    summarize(dry_mass_prop=mean(dry_mass_content,na.rm=T),n()) -> moisture
  
  #TODO ADD A REAL CALCULATION OF WOOD DENSITY ONCE WE UNDERSTAND HOW TO DO THAT
  df %>%
    left_join(moisture) %>%
    mutate(totalSampleDryMass=`Fresh mass (g)`*dry_mass_prop,size=size,density=NA,time=0,fruiting=NA,insects=NA,drill=NA) %>%
    select(unique, Species, size,time,totalSampleDryMass,density,fruiting,insects,drill) -> df_out
  
  return(df_out)
}

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


# Initial chemistry data

load_waterPercent<-function(){
  require(dplyr)
  
  # read in initial covariate data
  covar.big <-read.csv('data/covariates_bigStems.csv', stringsAsFactor = F)
  covar.small <-read.csv('data/covariates_smallStems.csv', stringsAsFactor = F)
  
  # calculate water content for each species, size class
  water.percent <- c(with(covar.big, (Fresh.mass..g. - Dry.mass..g.) / Dry.mass..g.),
                     with(covar.small, (Fresh.mass..g. - Dry.mass.total..g.) / Dry.mass.total..g.))
  water.percent <- data.frame(code=c(covar.big$Species, covar.small$Species),
                              StemSize=factor(c(rep('large', nrow(covar.big)), rep('small', nrow(covar.small)))),
                              water.percent, stringsAsFactors=F)
  
  #why is this so sparse?
  #View(water.percent)
  
  ## aggregate by code
  group_by(water.percent, code) %>%
    summarize(meanWaterPerc = mean(water.percent, na.rm=TRUE),
              sdWaterPerc = sd(water.percent, na.rm=TRUE)) -> water.percent.agg
              
  return(water.percent.agg)
  
}

load_densityNbarkthick<-function(){
  require(dplyr)
  
  # read in initial covariate data
  covar.small <-read.csv('data/covariates_smallStems.csv', stringsAsFactor = F)
  
  # calculate wood density and bark thickness for each species
  #only have these values measured on small stems
  
  covar.small$density.gpercm3 <- with(covar.small, Dry.mass.wood..g. / Volume..g.)
  covar.small$barkthickness.mm <- with(covar.small, Diameter.wbark..mm. - Diameter.nobark..mm.)
  
  ## aggregate by species (all small size)
  group_by(covar.small, Species) %>%
    summarize(meanDensity = mean(density.gpercm3, na.rm=TRUE),
              sdDensity = sd(density.gpercm3, na.rm=TRUE),
              meanBarkthick = mean(barkthickness.mm, na.rm=TRUE),
              sdBarkthick = sd(barkthickness.mm, na.rm=TRUE)) -> physTraits.agg
  colnames(physTraits.agg)[1]<-"code"
  
  return(physTraits.agg)
  
}

load_XRF<-function(){
  require(dplyr)
  
  data <- read.csv('data/sequencing_T0/NextGenSeqencing_2016_Sample_MasterSpreadsheet.csv', stringsAsFactors=F)
  data <- data[!data$SampleCode == 'blank', ] # extract 'blank' and 'mock' samples from 'mat.otu', delete from 'data'
  
  # create dataframe containing metadata for initial sequencing and XRF data
  meta <- data[, c('SampleCode', 'StemSize', 'mgSample', 'NucleicAcidConc', 'ExtractionDate')]
  meta$code <- substr(meta$SampleCode, 1, 4)
  # add column to 'meta' indicating whether data obtained from independent/composite sample
  meta$compositeSample <- T
  meta$compositeSample[grep('[0-9A-Z]$', meta$SampleCode)] <- F #if there is a number on the back of the code, then it was composited
  
  #isolate XRF cols
  df.xrf<-data.frame(SampleCode=data$SampleCode, data[, 26:ncol(data)])
  indx<-meta[,c("SampleCode","StemSize","code","compositeSample")]
  df.xrf1<-left_join(indx,df.xrf)
  
  ## identify samples that are not composited and do the aggregation by code
  temp<-df.xrf1[df.xrf1$compositeSample==FALSE, ]
  group_by(temp, code) %>%
    summarize(meanP = mean(P, na.rm=TRUE),
              sdP = sd(P, na.rm=TRUE),
              
              meanK = mean(K, na.rm=TRUE),
              sdK = sd(K, na.rm=TRUE),
              
              meanCa = mean(Ca, na.rm=TRUE),
              sdCa = sd(Ca, na.rm=TRUE),
              
              meanMn = mean(Mn, na.rm=TRUE),
              sdMn = sd(Mn, na.rm=TRUE),
              
              meanFe = mean(Fe, na.rm=TRUE),
              sdFe = sd(Fe, na.rm=TRUE),
              
              meanZn = mean(Zn, na.rm=TRUE),
              sdZn = sd(Zn, na.rm=TRUE)
              
              ) -> temp.agg

  aggSamps<-temp.agg[,c("code","meanP","meanK","meanCa","meanMn","meanFe","meanZn")]
  colnames(aggSamps)<-c("code","P","K","Ca","Mn","Fe","Zn")
  
  #combine composite and aggregated samples
  compSamps<-df.xrf1[df.xrf1$compositeSample==TRUE,c("code","P","K","Ca","Mn","Fe","Zn")]
  xrfSamps<-data.frame(rbind(aggSamps,compSamps))
  
  return(xrfSamps)
}

load_CN<-function(){
  require(dplyr)
  
  # read in CN data
  cndata <- read.csv('data/CN/JEFF_POWELL_CN_DATA_DIVERSITY_ROT_AUG_2013.csv', stringsAsFactors=F)
  colnames(cndata)<-c("sampleID","comments","mass","n.perc","c.perc")
  
  #create meta from XRF meta data
  data <- read.csv('data/sequencing_T0/NextGenSeqencing_2016_Sample_MasterSpreadsheet.csv', stringsAsFactors=F)
  data <- data[!data$SampleCode == 'blank', ] # extract 'blank' and 'mock' samples from 'mat.otu', delete from 'data'
  meta <- data[, c('SampleCode', 'StemSize', 'mgSample', 'NucleicAcidConc', 'ExtractionDate')]
  meta$code <- substr(meta$SampleCode, 1, 4)
  # add column to 'meta' indicating whether data obtained from independent/composite sample
  meta$compositeSample <- T
  meta$compositeSample[grep('[0-9A-Z]$', meta$SampleCode)] <- F #if there is a number on the back of the code, then it was composited
  
  ## fix sample names for composited samples
  indx<-meta[,c("SampleCode","StemSize","code","compositeSample")]
  temp<-merge(cndata, indx, by.x="sampleID",by.y="SampleCode", all.x=TRUE)
  
  x<-temp[,"sampleID"]
  x
  xx<-unlist(strsplit(x,"pooled"))
  xx
  xx[xx=="cali "]<-"cali"
  xx[xx=="acelextra"]<-"acel"
  temp[,"code"]<-substr(xx, 1,4)
  temp[,"compositeSample"]<-TRUE
  
  
  
  
  ## look up code in meta to fill in species and speciesStemSize
  temp<-temp[,c("sampleID","n.perc","c.perc","StemSize","code","compositeSample")]
  indx<-unique(meta[,c("code","species","speciesStemSize")])
  temp<-merge(temp,indx, by="code", all.x=TRUE)
  ## identify the composite samples
  pooledCodes<-codes$code %in% temp[temp$compositeSample==TRUE,"code"] 
  codes$code[pooledCodes]
  ## identify samples that are not composited, but for which there is a pooled sample
  tmp<-temp[temp$compositeSample==FALSE & temp$code %in% codes$code[pooledCodes], ]
  ddply(tmp, ~code, summarize, numReps=sum(!is.na(sampleID))) #number of samples per code.. there are not 3 for all of them, so I am just going to use the composite values
  ## identify samples that are not composited and there is no pooled sample
  tmp<-temp[temp$compositeSample==FALSE & !temp$code %in% codes$code[pooledCodes], ]
  tmp$speciesStemSize<-paste(tmp$species, tmp$StemSize, sep="__")
  # aggregate values by code
  tmp.agg <- summaryBy(n.perc + c.perc ~ speciesStemSize, data=tmp, FUN=c(mean, sd), na.rm=T)
  indx<-unique(tmp[,c("speciesStemSize","code")])
  tmp2<-left_join(tmp.agg, indx)
  aggSamps<-tmp2[,c("speciesStemSize","n.perc.mean","c.perc.mean","code")]
  colnames(aggSamps)<-c("speciesStemSize","n.perc","c.perc","code")
  #combine composite and aggregated samples
  compSamps<-temp[temp$compositeSample==TRUE,c("speciesStemSize","n.perc","c.perc","code")]
  cnSamps<-data.frame(rbind(aggSamps,compSamps))
  
}







####### NEED TO MODIFY




  
 
  
  # integrate species-level estimates of wood chemistry data into traits table ([NO EUSC DATA FOR SMALL SIZE CLASS] maybe only use from 'small' size class? - no interaction observed in adonis analysis [in 'analysis_T0.R'], only for K in univariate analysis)


 
  
  
  #merge together xrfSamps and cnSamps
#   xrfcnSamps<-left_join(xrfSamps,cnSamps)
#   
#   # merge together water percent, traits, and xrf.mean
#   tmp<-left_join(water.percent.agg, xrfcnSamps)
#   tmp<-separate(tmp, col=speciesStemSize, into=c("species","StemSize"), sep="__", remove=FALSE)
#   traits<-left_join(tmp, physTraits.agg)
#   traits.mean<-traits[,c("species","StemSize",
#                          "density.gpercm3.mean","barkthickness.mm.mean","water.percent.mean",
#                          "P","K","Ca","Mn","Fe","Zn","n.perc","c.perc")]
#   colnames(traits.mean)<-c("species","size",
#                            "density","barkthick","waterperc",
#                            "P","K","Ca","Mn","Fe","Zn","N","C")
#   dim(traits.mean)
#   tmp<-gather(traits.mean, key=trait, value=value, 3:13)
#   ddply(tmp, ~trait, summarize,
#         max=range(value, na.rm=TRUE)[1],
#         min=range(value, na.rm=TRUE)[2])
#   
#   
# 
# 
# 
# 
# 
# 
# 
# 
# #based on woodDecay_dataPrep_sequencing_T0_refactor.R
# load_matotu<-function(){
#   
#   # read in taxonomy
#   tax <- read.delim('data/sequencing_T0/DP16_tax.txt', stringsAsFactors=F)
#   
#   # read in OTU table (uclust output) and convert to matrix (rows=samples, columns=OTUs)
#   data.otu <- read.csv('data/sequencing_T0/DP16_OTUtable.csv')
#   mat.otu <- as.matrix(data.otu[, 2:ncol(data.otu)]); rownames(mat.otu) <- data.otu[, 1]
#   
#   # read in dataframe that contains sample information (also used to create meta and xrf)
#   data <- read.csv('data/sequencing_T0/NextGenSeqencing_2016_Sample_MasterSpreadsheet.csv', stringsAsFactors=F)
#   # extract 'blank' and 'mock' samples from 'mat.otu', delete from 'data'
#   data <- data[!data$SampleCode == 'blank', ]
#   
#   # re-label rownames in 'mat.otu' with sample codes
#   rownames(mat.otu) <- gsub('.', '-', rownames(mat.otu), fixed=T)
#   rownames(mat.otu)[match(data$NextGenID, rownames(mat.otu))] <- data$SampleCode
#   
#   # extract 'blank' and 'mock' samples from 'mat.otu', delete from 'data'
#   blank <- mat.otu[grep('blank', rownames(mat.otu)), ]
#   mock <- mat.otu['mock', ]
#   mat.otu <- mat.otu[-c(grep('blank', rownames(mat.otu)), grep('mock', rownames(mat.otu))), ]
#   
#   # otus, taxa in mock (select cut-off of >=9 reads in a sample)
#   mock <- data.frame(reads=sort(mock[mock > 0]))
#   mock <- cbind(mock, tax[match(rownames(mock), tax$OTU), 'species'])
#   # mock
#   mat.otu[mat.otu < 9] <- 0
#   
#   # re-order rows in 'mat.otu' to match rows in 'data'
#   mat.otu <- mat.otu[match(data$SampleCode, rownames(mat.otu)), ]
#   all(rownames(mat.otu) == data$SampleCode)  # is TRUE
#   
#   return(mat.otu)
#   
# }
# 
# #based on woodDecay_dataPrep_traits_T0.R
# load_traits_mean <- function(harvest, meta) {
#   
#   # get "codes" from harvest
#   spl <- split(harvest[, c('code', 'species', 'family', 'site','size')], harvest$code) 
#   codes <- do.call(rbind, lapply(spl, function(x) {x[1, ]})) 
#   codes$code <- as.character(codes$code) 
#   
#   # read in initial covariate data
#   covar.big <-read.csv('data/covariates_bigStems.csv', stringsAsFactor = F)
#   covar.small <-read.csv('data/covariates_smallStems.csv', stringsAsFactor = F)
#   
#   # calculate water content for each species, size class (SIGNIFICANT INTERACTION TERM, HOW TO COMBINE; IS THIS EVEN INFORMATIVE?)
#   water.percent <- c(with(covar.big, (Fresh.mass..g. - Dry.mass..g.) / Dry.mass..g.),
#                      with(covar.small, (Fresh.mass..g. - Dry.mass.total..g.) / Dry.mass.total..g.))
#   
#   water.percent <- data.frame(code=c(covar.big$Species, covar.small$Species),
#                               StemSize=factor(c(rep('large', nrow(covar.big)), rep('small', nrow(covar.small)))),
#                               water.percent, stringsAsFactors=F)
#   #m1 <- lm(water.percent ~ species * StemSize, data=water.percent)
#   #Anova(m1)
#   #visreg(m1, xvar='species', by='StemSize', overlay=T)
#   indx<-codes[,c("code","size","species")]
#   water.percent.tmp <- left_join(indx, water.percent)
#   water.percent.tmp$speciesStemSize<-paste(water.percent.tmp$species, water.percent.tmp$StemSize, sep="__") #why is this data so sparse?
#   ## by species x stemSize
#   water.percent.agg <-
#     summaryBy(
#       water.percent ~ speciesStemSize,
#       data = water.percent.tmp,
#       FUN = c(mean, sd),
#       na.rm = T
#     )
#   dim(water.percent.agg)
#   
#   # calculate wood density and bark thickness for each species
#   #only have these values measured on small stems
#   covar.small$density.gpercm3 <-
#     with(covar.small, Dry.mass.wood..g. / Volume..g.)
#   #hist(covar.small$density.gpercm3)
#   covar.small$barkthickness.mm <-
#     with(covar.small, Diameter.wbark..mm. - Diameter.nobark..mm.)
#   #hist(covar.small$barkthickness.mm)
#   covar.small <-
#     left_join(covar.small, codes, by = c('Species' = 'code'))
#   # summary(lm(density.gpercm3 ~ species, data=covar.small))
#   # summary(lm(barkthickness.mm ~ species, data=covar.small))
#   physTraits.agg <-
#     summaryBy(
#       density.gpercm3 + barkthickness.mm ~ species,
#       data = covar.small,
#       FUN = c(mean, sd),
#       na.rm = T
#     )
#   dim(physTraits.agg)
#   
#   # integrate species-level estimates of wood chemistry data into traits table ([NO EUSC DATA FOR SMALL SIZE CLASS] maybe only use from 'small' size class? - no interaction observed in adonis analysis [in 'analysis_T0.R'], only for K in univariate analysis)
#   # read in dataframe that contains sample information and xrf data (also used to create meta)
#   data <- read.csv('data/sequencing_T0/NextGenSeqencing_2016_Sample_MasterSpreadsheet.csv', stringsAsFactors=F)
#   df.xrf<-data.frame(SampleCode=data$SampleCode, data[, 26:ncol(data)])
#   df.xrf<-df.xrf[df.xrf$SampleCode!="blank",]
#   indx<-meta[,c("SampleCode","speciesStemSize","compositeSample")]
#   df.xrf1<-left_join(indx,df.xrf)
#   # deal with composite samples in xrf
#   ## add code to df.xrf1
#   indx<-meta[,c("SampleCode","StemSize","code","species")]
#   df.xrf2<-left_join(df.xrf1, indx)
#   ## identify the composite samples
#   pooledCodes<-codes$code %in% df.xrf2[df.xrf2$compositeSample==TRUE,"SampleCode"] 
#   codes$code[pooledCodes]
#   ## identify samples that are not composited, but for which there is a pooled sample
#   temp<-df.xrf2[df.xrf2$compositeSample==FALSE & df.xrf2$code %in% codes$code[pooledCodes], ]
#   ddply(temp, ~code, summarize, numReps=sum(!is.na(SampleCode))) #number of samples per code.. there are not 3 for all of them, so I am just going to use the composite values
#   ## identify samples that are not composited and there is no pooled sample
#   temp<-df.xrf2[df.xrf2$compositeSample==FALSE & !df.xrf2$code %in% codes$code[pooledCodes], ]
#   # aggregate values by code
#   temp.agg <- summaryBy(P + K + Ca + Mn + Fe + Zn ~ speciesStemSize, data=temp, FUN=c(mean, sd), na.rm=T)
#   indx<-unique(temp[,c("speciesStemSize","code")])
#   temp2<-left_join(temp.agg, indx)
#   aggSamps<-temp2[,c("speciesStemSize","P.mean","K.mean","Ca.mean","Mn.mean","Fe.mean","Zn.mean","code")]
#   colnames(aggSamps)<-c("speciesStemSize","P","K","Ca","Mn","Fe","Zn","code")
#   #combine composite and aggregated samples
#   compSamps<-df.xrf2[df.xrf2$compositeSample==TRUE,c("speciesStemSize","P","K","Ca","Mn","Fe","Zn","code")]
#   xrfSamps<-data.frame(rbind(aggSamps,compSamps))
#   
#   # read in CN data
#   cndata <- read.csv('data/CN/JEFF_POWELL_CN_DATA_DIVERSITY_ROT_AUG_2013.csv', stringsAsFactors=F)
#   colnames(cndata)<-c("sampleID","comments","mass","n.perc","c.perc")
#   # deal with composite samples in xrf
#   ## fix sample names for composited/pooled samples
#   indx<-meta[,c("SampleCode","StemSize","code","species","speciesStemSize","compositeSample")]
#   temp<-merge(cndata, indx, by.x="sampleID",by.y="SampleCode", all.x=TRUE)
#   x<-temp[is.na(temp$species),"sampleID"]
#   xx<-unlist(strsplit(x,"pooled"))
#   xx<-tolower(xx)
#   xx[xx=="cali "]<-"cali"
#   xx[xx=="acelextra"]<-"acel"
#   temp[is.na(temp$species),"code"]<-xx
#   temp[is.na(temp$species),"StemSize"]<-"small"
#   temp[is.na(temp$species),"compositeSample"]<-TRUE
#   ## look up code in meta to fill in species and speciesStemSize
#   temp<-temp[,c("sampleID","n.perc","c.perc","StemSize","code","compositeSample")]
#   indx<-unique(meta[,c("code","species","speciesStemSize")])
#   temp<-merge(temp,indx, by="code", all.x=TRUE)
#   ## identify the composite samples
#   pooledCodes<-codes$code %in% temp[temp$compositeSample==TRUE,"code"] 
#   codes$code[pooledCodes]
#   ## identify samples that are not composited, but for which there is a pooled sample
#   tmp<-temp[temp$compositeSample==FALSE & temp$code %in% codes$code[pooledCodes], ]
#   ddply(tmp, ~code, summarize, numReps=sum(!is.na(sampleID))) #number of samples per code.. there are not 3 for all of them, so I am just going to use the composite values
#   ## identify samples that are not composited and there is no pooled sample
#   tmp<-temp[temp$compositeSample==FALSE & !temp$code %in% codes$code[pooledCodes], ]
#   tmp$speciesStemSize<-paste(tmp$species, tmp$StemSize, sep="__")
#   # aggregate values by code
#   tmp.agg <- summaryBy(n.perc + c.perc ~ speciesStemSize, data=tmp, FUN=c(mean, sd), na.rm=T)
#   indx<-unique(tmp[,c("speciesStemSize","code")])
#   tmp2<-left_join(tmp.agg, indx)
#   aggSamps<-tmp2[,c("speciesStemSize","n.perc.mean","c.perc.mean","code")]
#   colnames(aggSamps)<-c("speciesStemSize","n.perc","c.perc","code")
#   #combine composite and aggregated samples
#   compSamps<-temp[temp$compositeSample==TRUE,c("speciesStemSize","n.perc","c.perc","code")]
#   cnSamps<-data.frame(rbind(aggSamps,compSamps))
#   
#   #merge together xrfSamps and cnSamps
#   xrfcnSamps<-left_join(xrfSamps,cnSamps)
#   
#   # merge together water percent, traits, and xrf.mean
#   tmp<-left_join(water.percent.agg, xrfcnSamps)
#   tmp<-separate(tmp, col=speciesStemSize, into=c("species","StemSize"), sep="__", remove=FALSE)
#   traits<-left_join(tmp, physTraits.agg)
#   traits.mean<-traits[,c("species","StemSize",
#                          "density.gpercm3.mean","barkthickness.mm.mean","water.percent.mean",
#                          "P","K","Ca","Mn","Fe","Zn","n.perc","c.perc")]
#   colnames(traits.mean)<-c("species","size",
#                            "density","barkthick","waterperc",
#                            "P","K","Ca","Mn","Fe","Zn","N","C")
#   dim(traits.mean)
#   tmp<-gather(traits.mean, key=trait, value=value, 3:13)
#   ddply(tmp, ~trait, summarize,
#         max=range(value, na.rm=TRUE)[1],
#         min=range(value, na.rm=TRUE)[2])
#   
#   
#   return(traits.mean)
# }
