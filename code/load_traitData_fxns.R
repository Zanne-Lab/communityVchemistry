
# species+size-level trait data

load_waterPercent.perGwetmass<-function(){
  require(dplyr)
  
  # read in initial covariate data
  covar.big <-read.csv('data/covariates_bigStems.csv', stringsAsFactor = F)
  covar.small <-read.csv('data/covariates_smallStems.csv', stringsAsFactor = F)
  
  # (1) calculate water content for each species, size class
  # g water per g dry mass
  water.percent <- c(with(covar.big, ((Fresh.mass..g. - Dry.mass..g.) / Fresh.mass..g.) * 100),
                     with(covar.small, ((Fresh.mass..g. - Dry.mass.total..g.) / Fresh.mass..g.)*100))
  water.percent <- data.frame(code=c(covar.big$Species, covar.small$Species),
                              StemSize=factor(c(rep('large', nrow(covar.big)), rep('small', nrow(covar.small)))),
                              water.percent, stringsAsFactors=F)
  
  ## aggregate by code
  group_by(water.percent, code) %>%
    summarize(meanWaterPerc = mean(water.percent, na.rm=TRUE),
              sdWaterPerc = sd(water.percent, na.rm=TRUE)) -> waterPercent
  
  #remove sd cols
  waterPercent<-waterPercent[,c("code","meanWaterPerc")]
  colnames(waterPercent)<-c("code","waterperc")
  length(unique(waterPercent$code)) #34
  
  return(waterPercent)
  
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
              sdBarkthick = sd(barkthickness.mm, na.rm=TRUE)) -> densityNbarkthick
  colnames(densityNbarkthick)[1]<-"code"
  
  #remove sd cols
  densityNbarkthick<-densityNbarkthick[,c("code","meanDensity","meanBarkthick")]
  colnames(densityNbarkthick)<-c("code","density","barkthick")
  length(unique(densityNbarkthick$code)) #22
  
  return(densityNbarkthick)
  
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
  
  group_by(df.xrf1, code) %>%
    summarize(compos=paste(unique(compositeSample), collapse="_")) -> summ
  #why are there some species+size that are composited == TRUE and FALSE?
  
  ## aggregate by code
  group_by(df.xrf1, code) %>%
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
              
    ) -> xrf
  
  xrf<-xrf[,c("code","meanP","meanK","meanCa","meanMn","meanFe","meanZn")]
  colnames(xrf)<-c("code","P","K","Ca","Mn","Fe","Zn")
  
  return(xrf)
}

load_CN<-function(){
  require(dplyr)
  
  # read in CN data
  cndata <- read.csv('data/CN/JEFF_POWELL_CN_DATA_DIVERSITY_ROT_AUG_2013.csv', stringsAsFactors=F)
  colnames(cndata)<-c("sampleID","comments","mass","n.perc","c.perc")
  
  #create meta from XRF meta data
  data <- read.csv('data/sequencing_T0/NextGenSeqencing_2016_Sample_MasterSpreadsheet.csv', stringsAsFactors=F)
  data <- data[!data$SampleCode == 'blank', ] # extract 'blank' and 'mock' samples from 'mat.otu', delete from 'data'
  meta <- data[, c('SampleCode', 'StemSize', 'StemCode', 'mgSample', 'NucleicAcidConc', 'ExtractionDate')]
  meta$code <- substr(meta$SampleCode, 1, 4)
  meta$Stem<-as.numeric(substr(meta$SampleCode, 5, 6)) # these numbers correspond to Stem, right? NAs are because character missing or character was a capital letter (not number)
  
  # add column to 'meta' indicating whether data obtained from independent/composite sample
  meta$compositeSample <- T
  meta$compositeSample[grep('[0-9A-Z]$', meta$SampleCode)] <- F #if there is a number on the back of the code, then it was composited
  
  ## fix code for composited samples
  indx<-meta[,c("SampleCode","StemSize","Stem","code","compositeSample")]
  temp<-merge(cndata, indx, by.x="sampleID",by.y="SampleCode", all.x=TRUE)
  temp[is.na(temp$compositeSample),"compositeSample"]<-TRUE
  x<-temp[temp$compositeSample==TRUE,"sampleID"]
  xx<-unlist(strsplit(x,"pooled"))
  xx<-tolower(xx)
  xx[xx=="cali "]<-"cali"
  xx[xx=="acelextra"]<-"acel"
  temp[temp$compositeSample==TRUE,"code"]<-substr(xx, 1,4)
  temp[temp$compositeSample==TRUE,"StemSize"]<-"small"
  
  ## make codeStem col
  temp$codeStem<-paste(temp$code, temp$Stem, sep="")
  
  # aggregate values by code
  group_by(temp, code) %>%
    summarize(meanN = mean(n.perc, na.rm=TRUE),
              sdN = sd(n.perc, na.rm=TRUE),
              
              meanC = mean(c.perc, na.rm=TRUE),
              sdC = sd(c.perc, na.rm=TRUE)
              
    ) -> temp.agg
  aggSamps<-temp.agg[,c("code","meanN","meanC")]
  colnames(aggSamps)<-c("code","n.perc","c.perc")
  cn<-data.frame(aggSamps)
  length(unique(cn$code)) #33
  
  return(cn)
  
}

mergeTraitData<-function(){
  
  #load code-aggregated trait data
  waterPercent<-load_waterPercent.perGwetmass() ##### this is in units of g water per g of wet mass x 100
  densityNbarkthick<-load_densityNbarkthick()
  xrf<-load_XRF()
  cn<-load_CN()
  
  # merge together water percent, traits, and xrf.mean by 'code'
  tmp<-left_join(waterPercent, densityNbarkthick)
  tmp<-left_join(tmp, xrf)
  traits<-left_join(tmp, cn)
  
  # rename columns
  traits<-rename(traits, "N"="n.perc", "C"="c.perc")
  
  # use species-level small-stem estimates of density and barkthick for large-stem samples
  traits$species<-tolower(traits$code)
  traits$size<-"small"
  traits[tolower(traits$code)!=traits$code,"size"]<-"large"
  largeCodes<-as.data.frame(traits[traits$size=="large","code"])[,1]
  for(i in 1:length(largeCodes)){
    curr.code<-largeCodes[i]
    curr.species<-tolower(curr.code)
    filter(traits, species==curr.species) %>%
      filter(size=="small") -> curr.row
    curr.data<-data.frame(curr.row[,c("density","barkthick")])
    traits[traits$code == curr.code, c("density","barkthick")]<-curr.data
  }
  
  traits<-traits[,c("code", "species","size",
                    "waterperc","density","barkthick","P" ,"K" ,"Ca" ,"Mn","Fe","Zn" ,"N","C")]
  
  # #summarize trait ranges
  # traits.long<-as.data.frame(gather(traits, key=trait, value=value, -(1:3)))
  # group_by(traits.long, trait) %>%
  #   summarize(max=range(value, na.rm=TRUE)[1],
  #             min=range(value, na.rm=TRUE)[2])
  
  return(traits)
  
}


# stem-level trait data

load_waterPercent.perGwetmass_byStem<-function(){
  require(dplyr)
  
  # read in initial covariate data
  covar.big <-read.csv('data/covariates_bigStems.csv', stringsAsFactor = F)
  covar.small <-read.csv('data/covariates_smallStems.csv', stringsAsFactor = F)
  
  # (1) calculate water content for each species, size class
  # g water per g dry mass
  water.percent <- c(with(covar.big, ((Fresh.mass..g. - Dry.mass..g.) / Fresh.mass..g.) * 100),
                     with(covar.small, ((Fresh.mass..g. - Dry.mass.total..g.) / Fresh.mass..g.)*100))
  water.percent <- data.frame(code=c(covar.big$Species, covar.small$Species),
                              Stem=c(covar.big$Stem, covar.small$Stem),
                              StemSize=factor(c(rep('large', nrow(covar.big)), rep('small', nrow(covar.small)))),
                              water.percent, stringsAsFactors=F)
  water.percent$codeStem<-paste(water.percent$code, water.percent$Stem, sep="")
  
  ## aggregate by code+stem
  group_by(water.percent, codeStem) %>%
    summarize(meanWaterPerc = mean(water.percent, na.rm=TRUE),
              sdWaterPerc = sd(water.percent, na.rm=TRUE)) -> waterPercent
  
  #remove sd cols
  waterPercent<-waterPercent[,c("codeStem","meanWaterPerc")]
  colnames(waterPercent)<-c("codeStem","waterperc")
  length(unique(waterPercent$codeStem)) #117
  
  return(waterPercent)
  
}

load_densityNbarkthick_byStem<-function(){
  require(dplyr)
  
  # read in initial covariate data
  covar.small <-read.csv('data/covariates_smallStems.csv', stringsAsFactor = F)
  
  # calculate wood density and bark thickness for each species
  #only have these values measured on small stems
  
  covar.small$density.gpercm3 <- with(covar.small, Dry.mass.wood..g. / Volume..g.)
  covar.small$barkthickness.mm <- with(covar.small, Diameter.wbark..mm. - Diameter.nobark..mm.)
  
  ## aggregate by code+Stem (all small size)
  covar.small$codeStem<-paste(as.character(covar.small$Species), covar.small$Stem, sep="")
  group_by(covar.small, codeStem) %>%
    summarize(meanDensity = mean(density.gpercm3, na.rm=TRUE),
              sdDensity = sd(density.gpercm3, na.rm=TRUE),
              meanBarkthick = mean(barkthickness.mm, na.rm=TRUE),
              sdBarkthick = sd(barkthickness.mm, na.rm=TRUE)) -> densityNbarkthick
  
  #remove sd cols
  densityNbarkthick<-densityNbarkthick[,c("codeStem","meanDensity","meanBarkthick")]
  colnames(densityNbarkthick)<-c("codeStem","density","barkthick")
  length(unique(densityNbarkthick$codeStem)) #80
  
  return(densityNbarkthick)
  
}

load_XRF_byStem<-function(){
  require(dplyr)
  
  data <- read.csv('data/sequencing_T0/NextGenSeqencing_2016_Sample_MasterSpreadsheet.csv', stringsAsFactors=F)
  data <- data[!data$SampleCode == 'blank', ] # extract 'blank' and 'mock' samples from 'mat.otu', delete from 'data'
  
  # create dataframe containing metadata for initial sequencing and XRF data
  meta <- data[, c('SampleCode', 'StemSize', 'mgSample', 'NucleicAcidConc', 'ExtractionDate')]
  meta$code <- substr(meta$SampleCode, 1, 4)
  meta$Stem<-as.numeric(substr(meta$SampleCode, 5, 6)) # these numbers correspond to Stem, right? NAs are because character missing or character was a capital letter (not number)
  
  # add column to 'meta' indicating whether data obtained from independent/composite sample
  meta$compositeSample <- T
  meta$compositeSample[grep('[0-9A-Z]$', meta$SampleCode)] <- F #if there is a number on the back of the code, then it was composited
  
  #isolate XRF cols
  df.xrf<-data.frame(SampleCode=data$SampleCode, data[, 26:ncol(data)])
  indx<-meta[,c("SampleCode","StemSize","code","Stem","compositeSample")]
  df.xrf1<-left_join(indx,df.xrf)
  
  #make a codeStem column
  df.xrf1$codeStem<-paste(df.xrf1$code, df.xrf1$Stem, sep="")
  group_by(df.xrf1, codeStem) %>%
    summarize(compos=paste(unique(compositeSample), collapse="_")) -> summ
  #why are there some species+size that are composited == TRUE and FALSE?
  
  ## aggregate by code
  group_by(df.xrf1, codeStem) %>%
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
              
    ) -> xrf
  
  xrf<-xrf[,c("codeStem","meanP","meanK","meanCa","meanMn","meanFe","meanZn")]
  colnames(xrf)<-c("codeStem","P","K","Ca","Mn","Fe","Zn")
  
  return(xrf)
}

load_CN_byStem<-function(){
  require(dplyr)
  
  # read in CN data
  cndata <- read.csv('data/CN/JEFF_POWELL_CN_DATA_DIVERSITY_ROT_AUG_2013.csv', stringsAsFactors=F)
  colnames(cndata)<-c("sampleID","comments","mass","n.perc","c.perc")
  
  #create meta from XRF meta data
  data <- read.csv('data/sequencing_T0/NextGenSeqencing_2016_Sample_MasterSpreadsheet.csv', stringsAsFactors=F)
  data <- data[!data$SampleCode == 'blank', ] # extract 'blank' and 'mock' samples from 'mat.otu', delete from 'data'
  meta <- data[, c('SampleCode', 'StemSize', 'mgSample', 'NucleicAcidConc', 'ExtractionDate')]
  meta$code <- substr(meta$SampleCode, 1, 4)
  meta$Stem<-as.numeric(substr(meta$SampleCode, 5, 6)) # these numbers correspond to Stem, right? NAs are because character missing or character was a capital letter (not number)
  
  # add column to 'meta' indicating whether data obtained from independent/composite sample
  meta$compositeSample <- T
  meta$compositeSample[grep('[0-9A-Z]$', meta$SampleCode)] <- F #if there is a number on the back of the code, then it was composited
  
  ## fix code for composited samples
  indx<-meta[,c("SampleCode","StemSize","Stem","code","compositeSample")]
  temp<-merge(cndata, indx, by.x="sampleID",by.y="SampleCode", all.x=TRUE)
  temp[is.na(temp$compositeSample),"compositeSample"]<-TRUE
  x<-temp[temp$compositeSample==TRUE,"sampleID"]
  xx<-unlist(strsplit(x,"pooled"))
  xx<-tolower(xx)
  xx[xx=="cali "]<-"cali"
  xx[xx=="acelextra"]<-"acel"
  temp[temp$compositeSample==TRUE,"code"]<-substr(xx, 1,4)
  temp[temp$compositeSample==TRUE,"StemSize"]<-"small"
  
  ## make codeStem col
  temp$codeStem<-paste(temp$code, temp$Stem, sep="")
  
  # aggregate values by code
  group_by(temp, codeStem) %>%
    summarize(meanN = mean(n.perc, na.rm=TRUE),
              sdN = sd(n.perc, na.rm=TRUE),
              
              meanC = mean(c.perc, na.rm=TRUE),
              sdC = sd(c.perc, na.rm=TRUE)
              
    ) -> temp.agg
  aggSamps<-temp.agg[,c("codeStem","meanN","meanC")]
  colnames(aggSamps)<-c("codeStem","n.perc","c.perc")
  cn<-data.frame(aggSamps)
  length(unique(cn$codeStem)) #80
  
  return(cn)
  
}

mergeTraitData_byStem<-function(){
  
  #load code-aggregated trait data
  waterPercent<-load_waterPercent.perGwetmass_byStem() ##### this is in units of g water per g of wet mass x 100
  densityNbarkthick<-load_densityNbarkthick_byStem()
  xrf<-load_XRF_byStem()
  cn<-load_CN_byStem()
  
  # merge together water percent, traits, and xrf.mean by 'codeStem'
  tmp<-full_join(waterPercent, densityNbarkthick)
  tmp<-full_join(tmp, xrf)
  traits<-full_join(tmp, cn)
  traits[is.nan(traits$waterperc),"waterperc"]<-NA
  traits[is.nan(traits$density),"density"]<-NA
  traits[is.nan(traits$barkthick),"barkthick"]<-NA
  
  # rename columns
  traits<-rename(traits, "N"="n.perc", "C"="c.perc")
  
  #most complete codeStem trait dataset
  traits.codeStem<-traits[!grepl("NA", traits$codeStem),]
  
  #add code and size
  traits.codeStem$code <- substr(traits.codeStem$codeStem, 1, 4)
  traits.codeStem$size<-"small"
  traits.codeStem[traits.codeStem$code!=tolower(traits.codeStem$code),"size"]<-"large"
  
  traits.codeStem<-traits.codeStem[,c("codeStem", "code","size",
                                      "waterperc","density","barkthick","P" ,"K" ,"Ca" ,"Mn","Fe","Zn" ,"N","C")]
  
  return(traits.codeStem)
  
}







### not used
# load_waterPercent.perGdrymass<-function(){
#   require(dplyr)
#   
#   # read in initial covariate data
#   covar.big <-read.csv('data/covariates_bigStems.csv', stringsAsFactor = F)
#   covar.small <-read.csv('data/covariates_smallStems.csv', stringsAsFactor = F)
#   
#   # (1) calculate water content for each species, size class
#   # g water per g dry mass
#   water.percent <- c(with(covar.big, ((Fresh.mass..g. - Dry.mass..g.) / Dry.mass..g.) * 100),
#                      with(covar.small, ((Fresh.mass..g. - Dry.mass.total..g.) / Dry.mass.total..g.) * 100))
#   water.percent <- data.frame(code=c(covar.big$Species, covar.small$Species),
#                               StemSize=factor(c(rep('large', nrow(covar.big)), rep('small', nrow(covar.small)))),
#                               water.percent, stringsAsFactors=F)
#   
#   ## aggregate by code
#   group_by(water.percent, code) %>%
#     summarize(meanWaterPerc = mean(water.percent, na.rm=TRUE),
#               sdWaterPerc = sd(water.percent, na.rm=TRUE)) -> waterPercent
#   
#   #remove sd cols
#   waterPercent<-waterPercent[,c("code","meanWaterPerc")]
#   colnames(waterPercent)<-c("code","waterperc")
#   length(unique(waterPercent$code)) #34
#   
#   return(waterPercent)
#   
# }

