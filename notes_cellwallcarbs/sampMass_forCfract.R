#

data<-read.csv("data/approxSampleMass_forCfract.csv")

#identify sample meta data
library(tidyr)
library(dplyr)
stemSamples<-read.csv("derived_data/stemSamples.csv")
seqSamples<-read.csv("derived_data/seqSamples.csv") 

#first of all, remember that seqSamples and trait samples (xrf) have 2 unique codeStems
seqSamples %>%
  filter(!codeStem %in% stemSamples$codeStem) %>%
  separate(seq_sampName, into=c("drop","seq.stem"), 4, remove=FALSE) %>%
  filter(grepl("[1-9]", seq.stem)) -> extra.codeStems
all.codeStems <- c(as.character(stemSamples$codeStem), as.character(extra.codeStems$seq_sampName)) # list of all unique codeStems

#add sample meta data and reorganize dataframe
code.indx<-unique(stemSamples[,c("code","species","size","Family","Binomial","site")])
data %>%
  separate(sampleID, into=c("code","other"), sep=4, remove=FALSE) %>%
  left_join(code.indx) %>%
  mutate(codeStem=ifelse(sampleID %in% all.codeStems, as.character(sampleID), NA)) %>%
  mutate(wt_in_grams0=ifelse(wt_in_grams < 0, 0, wt_in_grams)) %>% 
  mutate(compositeSamp = ifelse(is.na(codeStem), TRUE, FALSE)) %>%
  select(sampleID, wt_in_grams, wt_in_grams0, code, codeStem, compositeSamp, species, size, Family, Binomial, site) %>%
  arrange(sampleID) -> samp.wts

#add trait data



#reps from all species and sizes
samp.wts %>%
  select(sampleID, code, compositeSamp) %>%
  group_by(code, compositeSamp) %>%
  summarize(numSamps = length(sampleID)) -> summ.all
dim(summ.all)[1] == length(unique(summ.all$code))
summ.all %>%
  filter(compositeSamp == TRUE) %>%
  select(code) -> comp.codes
# 33 unique codes
# 10 are composited codes -- all from the small size class
sum(summ.all$numSamps) # 82 samples in total

#non-composited samples only
samp.wts %>%
  filter(compositeSamp == FALSE) %>%
  select(sampleID, code, compositeSamp) %>%
  group_by(code) %>%
  summarize(numSamps = length(sampleID)) %>%
  arrange(code) -> summ.noncomposited
dim(summ.noncomposited)[1] # 23 unique codes
sum(summ.noncomposited$numSamps) # 70 samples in total

# small size only
samp.wts %>%
  filter(size == "small") %>%
  select(sampleID, code, compositeSamp) %>%
  group_by(code, compositeSamp) %>%
  summarize(numSamps = length(sampleID)) -> summ.small
dim(summ.small)[1] == length(unique(summ.small$code))
dim(summ.small)[1] # 21 unique codes
summ.small %>%
  filter(compositeSamp == TRUE) %>%
  select(code) -> comp.codes
# 10 are composited codes -- all from the small size class
sum(summ.small$numSamps) # 45 samples in total

# large size only
samp.wts %>%
  filter(size == "large") %>%
  select(sampleID, code, compositeSamp) %>%
  group_by(code, compositeSamp) %>%
  summarize(numSamps = length(sampleID)) -> summ.large
dim(summ.large)[1] == length(unique(summ.large$code))
dim(summ.large)[1] # 12 unique codes
sum(summ.large$numSamps) # 37 samples in total






#---------------------------------#
  
#summarize study design
samp.wts %>%
  select(sampleID, Binomial, size, codeStem) %>%
  arrange(Binomial, size, codeStem) %>%
  group_by(Binomial, size) %>%
  summarize(numSamps = length(sampleID)) -> tmp
write.table(tmp, file = "data/sawDust_sampleSummary.txt")

#total number of samples
dim(samp.wts)[1]

#average sample weight
range(samp.wts$wt_in_grams)
mean(samp.wts$wt_in_grams0)

#sum sample weight by code
samp.wts %>%
  group_by(code) %>%
  summarize(pooled.wt=sum(wt_in_grams0)) %>%
  left_join(code.indx) %>%
  arrange(pooled.wt) -> samp.wts.byCode

#average pooled sample weight
hist(samp.wts.byCode$pooled.wt)
mean(samp.wts.byCode$pooled.wt)

#number of samples if pooled by code
dim(samp.wts.byCode)[1]

#samples less than 0.20 g
samp.wts %>%
  filter(wt_in_grams < 0.15)



