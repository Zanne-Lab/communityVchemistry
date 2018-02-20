#stem sample meta data
load_stemSamples<-function(){
  
  deployment <- read_csv("data/deployment.csv")
  deployment<-rename(deployment, "code"="species") #Code column
  
  #summarize by stem
  deployment %>%
    group_by(code, Stem) %>%
    summarize(num.unique=length(unique(unique))) %>%
    mutate(codeStem=paste(code, Stem, sep="")) %>%
    mutate(species=tolower(code)) %>%
    mutate(size=ifelse(code == tolower(code), 'small','large')) -> deploy.new
  
  #add species info
  species <- read_csv("data/species.csv")
  species[species$Binomial == "Ricinocarpus pinifolius", "Binomial"] <- "Ricinocarpos pinifolius" # fix a misspelled name
  species %>%
    mutate(species=tolower(Code)) %>%
    rename("site"="Collection site") %>%
    select(species, Family, Binomial, site) -> species.new
  deploy.new %>%
    left_join(species.new) -> stemSamples
  
  return(stemSamples)  

}

#OTU table
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
  #mock
  mat.otu[mat.otu < 9] <- 0
  
  # otus, taxa in blank
  blank <- data.frame(reads=sort(blank[blank > 0]))
  blank <- cbind(blank, tax[match(rownames(blank), tax$OTU), 'species'])
  #mat.otu[,'ITSall_OTUd_3713'] # the most abundant OTU in the blank does not show up in any of the samples
  
  # re-order rows in 'mat.otu' to match rows in 'data'
  mat.otu <- mat.otu[match(data$SampleCode, rownames(mat.otu)), ]
  all(rownames(mat.otu) == data$SampleCode)  # is TRUE
  
  #how many Glomeromycota are there?
  GlomOTUs<-tax[tax$phylum=="Glomeromycota","OTU"]
  glom.mat.otu<-mat.otu[,colnames(mat.otu) %in% GlomOTUs]
  rowSums(glom.mat.otu)!=0
  tmp<-glom.mat.otu[,colSums(glom.mat.otu)!=0] #get rid of cols and rows without reads
  glom.mat.otu.sel<-tmp[rowSums(tmp)!=0,]
  
  return(mat.otu)
  
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

#create sequence sample meta data table
load_seqSamples<-function(mat.otu, stemSamples){
  
  #identify sequence sampleIDs
  seq_indx<-data.frame(seq_sampName=row.names(mat.otu))
  
  #merge by codeStem
  stem.indx<-stemSamples[,c("codeStem","code","Stem")]
  left_join(seq_indx, stem.indx, by=c("seq_sampName"="codeStem")) %>%
    mutate(codeStem = ifelse(!is.na(Stem), paste(code, Stem, sep=""), NA)) %>%
    select(seq_sampName, codeStem) -> seqSamples.tmp
  #unique(seqSamples.tmp$codeStem) no stem-specific ripi
  seqSamples.tmp %>%
    separate(seq_sampName, into=c("code","extra"), 4, remove=FALSE) %>%
    select(-extra) -> seqSamples.tmp
  
  #add back in code-level information for seq_samples that have been pooled by code
  code.indx <- unique(stemSamples[,c("code","species","size")])
  code.indx %>%
    left_join(seqSamples.tmp) -> seqSamples
  
  return(seqSamples)
}

#taxon lookup info
load_TaxAndFunguild <- function(comm.otu) {
  require(dplyr)
  
  # load fungal OTU info
  funguild <-read.delim('data/sequencing_T0/DP16_funguild.txt', stringsAsFactors = F)
  tax <-read.delim('data/sequencing_T0/DP16_tax.txt', stringsAsFactors = F)
  
  # merge the dataframes by OTUId
  colnames(tax)[1] <- "OTUId" #make this match the column name in funguild
  taxAndFunguild <- left_join(tax, funguild)
  
  # create a kingdom column
  taxAndFunguild$kingdom<-NA
  taxAndFunguild[grepl("Fungi", taxAndFunguild$taxonomy),"kingdom"]<-"Fungi"
  
  # add oomycete OTUs as rows
  ooOTUs<-colnames(comm.otu)[!colnames(comm.otu) %in% taxAndFunguild$OTUId]
  oo.df<-data.frame(OTUId=ooOTUs, kingdom="Protist")
  taxAndFunguild<-bind_rows(taxAndFunguild, oo.df)
  
  # reorder taxAndFunguild to make OTU table
  o<-match(colnames(comm.otu), taxAndFunguild$OTUId)
  o.taxAndFunguild<-taxAndFunguild[o,]
  sum(o.taxAndFunguild$OTUId != colnames(comm.otu)) #this need to be 0
  
  # only use FUNGuild info with confidence ranking of Probable or Highly Probable
  o.taxAndFunguild[!o.taxAndFunguild$Confidence.Ranking %in% c("Probable","Highly Probable"),c("Trophic.Mode","Guild")]<-"unclassified"
  
  # select cols 
  o.taxAndFunguild %>%
    select(OTUId, taxonomy, kingdom, phylum, genus, species, 
           Trophic.Mode, Guild) -> o.taxAndFunguild
  
  #clean Trophic.Mode
  #unique(o.taxAndFunguild$Trophic.Mode)
  
  #clean Guild
  #unique(o.taxAndFunguild$Guild)
  o.taxAndFunguild[o.taxAndFunguild$Guild=="NULL","Guild"]<-"unclassified"
  
  #clean oomycetes
  o.taxAndFunguild[o.taxAndFunguild$kingdom=="Protist", c("taxonomy","phylum","genus","species")]<-"unclassified"
  
  return(o.taxAndFunguild)
}

#OTU sample-effort curves
plot_sampleEffortCurves<-function(mat.otu){
  
  pdf(file="output/sampleEffortCurve.pdf", width=5, height=5)
  
  rarecurve(mat.otu, step=100,
            xlab="Number of reads per sample", 
            ylab="Cumulative number of OTUs", label=FALSE)
  dev.off()
  
}
