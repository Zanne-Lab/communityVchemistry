#load data and clean initial dataframes

#############
# Sample codes

load_stemSamples<-function(){
  
  deployment <- read_csv("data/deployment.csv")
  deployment<-rename(deployment, "code"="species") #Code column
  
  #summarize by code
  deployment %>%
    group_by(code) %>%
    summarize(nStems=length(unique(unique))) %>%
    mutate(species=tolower(code)) -> deploy.new
  
  #add size
  deploy.new$size<-"large"
  deploy.new[tolower(deploy.new$code) == deploy.new$code, "size"]<-"small"
  
  #add species info
  species <- read_csv("data/species.csv")
  species %>%
    mutate(species=tolower(Code)) %>%
    select(-Code)-> species.new
  stemSamples<-left_join(deploy.new, species.new)
  
  return(stemSamples)
}


#############
# Sample mass & volume data

read_in_initial_mass <- function(){
  require(readr)
  require(dplyr)
  big <- read_csv("data/covariates_bigStems.csv")
  small <- read_csv("data/covariates_smallStems.csv")
  
  #look for missing data in olst small
  filter(small, Species == 'olst') -> tmp
  #View(tmp) #no samples with `Dry mass total (g)`
  
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
    select(unique, Species, size,time,totalSampleDryMass,density,fruiting,insects,drill) %>%
    rename("species"="Species") -> df_out
  
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
  samp4$weightForVol <- sapply(temp, function(x){if(length(x) == 2) as.numeric(x)[2] else as.numeric(x)[1]})
  x <- samp4$weightForVol %in% 0
  samp4[x, 'weightForVol'] <- samp4[x, 'wetWeight']
  rm(temp, x)
  # include excess wood in bag (fragments falling off during transit) for wet and dry mass of whole harvested piece
  samp4$wetWeightForMass <- apply(samp4[, c('wetWeight', 'wetWeightExcess')], 1, sum, na.rm=T)
  samp4[is.na(samp4$wetWeight), 'wetWeightForMass'] <- NA
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

CalcTotalDryMass<-function(data){
  data$totalSampleDryMass <- NA
  x <- which(data$drill == 'no' | data$dryMass == 0 | data$weightForVol == 0); data$totalSampleDryMass[x] <- data$dryMass[x]
  x <- which(data$drill == 'yes' & data$dryMass > 0 & data$weightForVol > 0); data$totalSampleDryMass[x] <- (data$weightForVol[x] * data$dryMass[x]) / data$wetWeightForMass[x]
  data$totalSampleDryMass <- round(data$totalSampleDryMass, 2)
  return(data)
}

CalcDensity<-function(data){
  data$density <- NA
  data$density <- round(data$weightForVol / data$volMass, 2)
  return(data)
}

ReorgDataFrame<-function(data){
  require(tidyr)
  
  #add a species column
  data<-separate(data, unique, into=c("species","extraCode"), 4, remove=FALSE)
  
  #add a size column
  data$size<-"large"
  data[tolower(data$species) == data$species,"size"]<-"small"
  
  #rename and select columns
  data %>%
    rename("fruiting"="fruitingBodies","insects"="insectDamage") %>%
    select(unique, species, size, time, totalSampleDryMass, density, fruiting, insects, drill, notes) -> data
  
  return(data)
  
}

LoadHarvestFiles<-function(){
  
  # load data
  s1 <- read.samp1()
  s2 <- read.samp2()
  s3 <- read.samp3()
  s4 <- read.samp4()
  
  #bind everything together
  s.data<-rbind(s1,s2,s3,s4)

  #calculate total dry mass
  s.data<-CalcTotalDryMass(s.data)

  #calculate density
  s.data<-CalcDensity(s.data)
  
  #reorganize data frame
  s.data<-ReorgDataFrame(s.data)
  
  #check for missing data
  #filter(s.data, is.na(totalSampleDryMass))
  filter(s.data, is.na(totalSampleDryMass), notes =="all wwe -- all wet weight excess") #what does this note mean?
  filter(s.data, is.na(totalSampleDryMass), notes !="all wwe -- all wet weight excess") #for missing samples, we don't know if they rotted away completely or were moved/missing, so entered as NA

  return(s.data)
  
}




#############
# Initial wood trait data

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
  
  #why is this so sparse? other stems were deployed into field
  #View(water.percent)
  
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

#issue --> why are there some species+size that are composited == TRUE and FALSE?
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
  meta <- data[, c('SampleCode', 'StemSize', 'mgSample', 'NucleicAcidConc', 'ExtractionDate')]
  meta$code <- substr(meta$SampleCode, 1, 4)
  # add column to 'meta' indicating whether data obtained from independent/composite sample
  meta$compositeSample <- T
  meta$compositeSample[grep('[0-9A-Z]$', meta$SampleCode)] <- F #if there is a number on the back of the code, then it was composited
  
  ## fix code for composited samples
  indx<-meta[,c("SampleCode","StemSize","code","compositeSample")]
  temp<-merge(cndata, indx, by.x="sampleID",by.y="SampleCode", all.x=TRUE)
  temp[is.na(temp$compositeSample),"compositeSample"]<-TRUE
  x<-temp[temp$compositeSample==TRUE,"sampleID"]
  xx<-unlist(strsplit(x,"pooled"))
  xx<-tolower(xx)
  xx[xx=="cali "]<-"cali"
  xx[xx=="acelextra"]<-"acel"
  temp[temp$compositeSample==TRUE,"code"]<-substr(xx, 1,4)
  temp[temp$compositeSample==TRUE,"StemSize"]<-"small"
  
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
  waterPercent<-load_waterPercent()
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
  traits$species_lower<-tolower(traits$code)
  traits$size<-"small"
  traits[tolower(traits$code)!=traits$code,"size"]<-"large"
  largeCodes<-as.data.frame(traits[traits$size=="large","code"])[,1]
  for(i in 1:length(largeCodes)){
    curr.code<-largeCodes[i]
    curr.species_lower<-tolower(curr.code)
    filter(traits, species_lower==curr.species_lower) %>%
      filter(size=="small") -> curr.row
    curr.data<-data.frame(curr.row[,c("density","barkthick")])
    traits[traits$code == curr.code, c("density","barkthick")]<-curr.data
  }
  
  traits<-traits[,c("code", "species_lower","size",
                    "waterperc","density","barkthick","P" ,"K" ,"Ca" ,"Mn","Fe","Zn" ,"N","C")]
  
  #summarize trait ranges
  traits.long<-as.data.frame(gather(traits, key=trait, value=value, -(1:3)))
  group_by(traits.long, trait) %>%
    summarize(max=range(value, na.rm=TRUE)[1],
              min=range(value, na.rm=TRUE)[2])
  
  return(traits)
  
}


#############
# Initial microbial community data

load_matotu<-function(){
  
  # read in OTU table (uclust output) and convert to matrix (rows=samples, columns=OTUs)
  data.otu <- read.csv('data/sequencing_T0/DP16_OTUtable.csv', stringsAsFactors = FALSE)
  mat.otu <- as.matrix(data.otu[, 2:ncol(data.otu)]); rownames(mat.otu) <- data.otu[, 1]
  
  sum(colSums(mat.otu)==0) # if 0, then there are no empty columns
  
  # read in dataframe that contains sample information (also used to create meta and xrf)
  data <- read.csv('data/sequencing_T0/NextGenSeqencing_2016_Sample_MasterSpreadsheet.csv', stringsAsFactors=F)
  # extract 'blank' and 'mock' samples from 'mat.otu', delete from 'data'
  data <- data[!data$SampleCode == 'blank', ]
  
  # re-label rownames in 'mat.otu' with sample codes
  rownames(mat.otu) <- gsub('.', '-', rownames(mat.otu), fixed=T)
  rownames(mat.otu)[match(data$NextGenID, rownames(mat.otu))] <- data$SampleCode
  
  # extract 'blank' and 'mock' samples from 'mat.otu', delete from 'data'
  blank <- mat.otu[grep('blank', rownames(mat.otu)), ]
  mock <- mat.otu['mock', ]
  mat.otu <- mat.otu[-c(grep('blank', rownames(mat.otu)), grep('mock', rownames(mat.otu))), ]
  
  # otus, taxa in mock (select cut-off of >=9 reads in a sample)
  tax <-read.delim('data/sequencing_T0/DP16_tax.txt', stringsAsFactors = F)
  mock <- data.frame(reads=sort(mock[mock > 0]))
  mock <- cbind(mock, tax[match(rownames(mock), tax$OTU), 'species'])
  # mock
  mat.otu[mat.otu < 9] <- 0
  
  # re-order rows in 'mat.otu' to match rows in 'data'
  mat.otu <- mat.otu[match(data$SampleCode, rownames(mat.otu)), ]
  all(rownames(mat.otu) == data$SampleCode)  # is TRUE
  
  return(mat.otu)
  
}

load_seqSamples<-function(mat.otu, stemSamples){
  stemSamples %>% select(code, species, size) -> codeindx
  seq_sampName<-row.names(mat.otu)
  seq_indx<-data.frame(seq_sampName=seq_sampName, code=substr(seq_sampName, 1, 4))
  seqSamples<-left_join(seq_indx, codeindx)
  return(seqSamples)
}

add_oomycetes<-function(fung.otu){
  
  # read in OTU table (uclust output) and convert to matrix (rows=samples, columns=OTUs)
  data.otu <- read.csv('data/sequencing_T0/OTUtable_oomycetes_20171020.csv', row.names=1)
  data.df<-data.frame(seqSamp=row.names(data.otu), data.otu)
  
  #make mat.otu of fungal taxa a dataframe
  fung.df<-data.frame(seqSamp=row.names(fung.otu), fung.otu)
  
  #merge by seqSamp
  comm.df<-left_join(fung.df, data.df)
  
  #make NAs into 0s
  comm.df[is.na(comm.df)]<-0
  
  #make dataframe into a matrix again
  row.names(comm.df)<-comm.df$seqSamp
  comm.mat<-comm.df[,-1]
  comm.mat<-as.matrix(comm.mat)
  
  return(comm.mat)
}

#need to add oomycete taxa
load_TaxAndFunguild <- function() {
  require(dplyr)
  
  funguild <-read.delim('data/sequencing_T0/DP16_funguild.txt', stringsAsFactors = F)
  tax <-read.delim('data/sequencing_T0/DP16_tax.txt', stringsAsFactors = F)
  
  # merge the dataframes by OTUId
  colnames(tax)[1] <- "OTUId" #make this match the column name in funguild
  taxAndFunguild <- full_join(tax, funguild)
  
  return(taxAndFunguild)
}


# calc_matotu_summStats<-function(mat.otu, meta){
#   
#   #total number of OTUs
#   totalOTUs<-dim(mat.otu)[2]
#   
#   #mean sample richness
#   richness<-apply(mat.otu,MARGIN = 1,function(x) sum(x>0))
#   meanRichness<-mean(richness)
#   seRichness<-sd(richness)/sqrt(length(richness))
#   
#   #mean sample evenness
#   H <- diversity(mat.otu) 
#   S <- specnumber(mat.otu) ## rowSums(BCI > 0) does the same...
#   J <- H/log(S) #Pielou's evenness (J)
#   meanJ<-mean(J)
#   seJ<-sd(J)/sqrt(length(J))
#   
#   #mean number of reads per sample
#   meanReads<-mean(rowSums(mat.otu))
#   seReads<-sd(rowSums(mat.otu))/sqrt(length(rowSums(mat.otu)))
#   seReads
#   
#   summaryStats<-data.frame(label=c("totalOTUs","meanRichness","seRichness","meanEvenness","seEvenness"),
#                            value=c(totalOTUs, meanRichness, seRichness, meanJ, seJ))
#   round(summaryStats$value, digits=4)
#   write.csv(summaryStats, file="output/matotu_summary.csv")
# }

plot_sampleEffortCurves<-function(mat.otu){
  
  pdf(file="output/sampleEffortCurve.pdf", width=5, height=5)
  
  rarecurve(mat.otu, step=100,
            xlab="Number of reads per sample", 
            ylab="Cumulative number of OTUs", label=FALSE)
  dev.off()
  
}


##############

#ggplot theme
mytheme <- theme_bw(base_size = 10, base_family = "Helvetica") +
  theme(panel.border = element_rect(colour = "black"),      #put a black box around the plotting area
        axis.line = element_line(colour = "black"),                 #axis lines are in black
        panel.grid.major = element_blank(),                         #turn off the gridlines
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(face='bold.italic', hjust=0.05),         #turn off the x axis facet labels
        strip.text.y = element_text(face='bold.italic', hjust=0.05),
        strip.background = element_rect(fill = 'white', colour='black'),    #make y axis facet labels be italic and top justified
        legend.key = element_blank(),                               #turn off box around legend
        plot.title=element_text(hjust=0, vjust=0.5, face='bold'), #style and position of the panel label
        plot.margin = unit(c(0.05,0.05,0.05,0.05),"in")
  )
        
        