#

data<-read.csv("data/approxSampleMass_forCfract.csv")

#identify sample meta data
library(tidyr)
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
  select(sampleID, bag, wt_in_grams, wt_in_grams0, code, codeStem, species, size, Family, Binomial, site) %>%
  arrange(sampleID) -> samp.wts

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



