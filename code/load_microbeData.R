#stem sample meta data
load_stemSamples <- function(){
  
  deployment <- read_csv("data/deployment.csv")
  deployment <- rename(deployment, "code"="species") #Code column
  
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

add_oomycetes <- function(fung.otu){

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

  #dim(fung.otu)
  #dim(comm.mat) #check that there are more columns for oomycete OTUs
  
  return(comm.mat)
}

#copied from fungal_wood_endophtyes repo
load_matotu <- function(){
  
  # read in OTU table (uclust output) and convert to matrix (rows=samples, columns=OTUs)
  data.otu <- read.csv('data/sequencing_T0/DP16_OTUtable.csv', stringsAsFactors = FALSE)
  mat.otu <- as.matrix(data.otu[, 2:ncol(data.otu)]); rownames(mat.otu) <- data.otu[, 1]
  # remove extra blank sample
  mat.otu <- mat.otu[-which(rownames(mat.otu) == 'X..blank'), ]
  mat.otu <- mat.otu[, colSums(mat.otu) > 0]
  
  sum(colSums(mat.otu)==0) # if 0, then there are no empty columns
  
  # read in dataframe that contains sample information (also used to create meta and xrf)
  data <- read.csv('data/sequencing_T0/NextGenSeqencing_2016_Sample_MasterSpreadsheet.csv', stringsAsFactors=F)
  
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
  mock <- cbind(mock, tax[match(rownames(mock), tax$qseqid), 'species'])
  #mock
  mat.otu[mat.otu < 9] <- 0
  
  # otus, taxa in blank
  blank <- data.frame(reads=sort(blank[blank > 0]))
  blank <- cbind(blank, tax[match(rownames(blank), tax$qseqid), 'species'])
  #mat.otu[,'ITSall_OTUd_3713'] # the most abundant OTU in the blank does not show up in any of the samples
  
  # re-order rows in 'mat.otu' to match rows in 'data' after deleting 'blank' from 'data'
  data <- data[!data$SampleCode == 'blank', ]
  mat.otu <- mat.otu[match(data$SampleCode, rownames(mat.otu)), ]
  all(rownames(mat.otu) == data$SampleCode)  # is TRUE
  
  # add oomycetes
  mat.otu <- add_oomycetes(fung.otu = mat.otu)
  sum(grepl("ITSoo", colnames(mat.otu)) == T) # verify that 13 oomycete OTUs have been added
  
  #trim OTUs that do not show up in ANY of the samples
  x <- apply(mat.otu > 0, 2, sum)
  mat.otu <- mat.otu[, x >= 1]
  
  return(mat.otu)
  
}

#copied from fungal_wood_endophtyes repo
load_TaxAndFunguild <- function(comm.otu.tmp){
  
  comm.otu.tmp <- comm.otu
  
  # load fungal OTU info
  funguild <-read.delim('data/sequencing_T0/DP16_funguild.txt', stringsAsFactors = F)
  tax <-read.delim('data/sequencing_T0/DP16_tax.txt', stringsAsFactors = F)
  
  # merge the dataframes by OTUId
  colnames(tax)[1] <- "OTUId" #make this match the column name in funguild
  tax %>%
    left_join(funguild) -> taxAndFunguild
  
  # apply trophic assignment for 1 genus to all
  GENUS <- unique(taxAndFunguild$genus)
  #i<-62
  for(i in 1:length(GENUS)){
    
    # isolate everything with the same genus
    criteria <- taxAndFunguild$genus == GENUS[i] & !is.na(taxAndFunguild$genus)
    select.rows <- taxAndFunguild[criteria,]
    
    # if it is more than 2 OTUs...
    dim(select.rows)[1] > 1
    select.rows
    # if there is something in Trophic.Mode
    criteria <- sum(!is.na(select.rows$Trophic.Mode))
    criteria > 0
    
    if(dim(select.rows)[1] > 1 & criteria > 0){
      
      troph.options <- unique(select.rows$Trophic.Mode)
      troph.options <- troph.options[!is.na(troph.options)]
      
      length(troph.options)
      if(length(troph.options) == 1){
        
        troph.apply <- troph.options
        taxAndFunguild[criteria, "Trophic.Mode"] <- troph.apply
        
      }else{
        print(paste("WARNING", i))
        print(troph.options)
      }
      
    }
  }
  
  # currently there is not oomycete information in the taxAndFunguild table
  sum(grepl("ITSoo", taxAndFunguild$OTUId))
  
  # add oomycete OTUs as rows
  ooOTUs <- colnames(comm.otu.tmp)[grepl("ITSoo", colnames(comm.otu.tmp))]
  oo.df<-data.frame(OTUId=ooOTUs, kingdom="Protist")
  taxAndFunguild<-bind_rows(taxAndFunguild, oo.df)
  sum(grepl("ITSoo", taxAndFunguild$OTUId))
  
  # identify OTUs from comm.otu not found in taxAndFunguild (probably plant DNA)
  taxAndFunguild %>%
    filter(is.na(coverage)) # the oomycetes don't have coverage data; everything else should
  taxAndFunguild %>%
    filter(kingdom != "Protist") %>%
    filter(!OTUId %in% colnames(comm.otu.tmp)) %>%
    select(OTUId, kingdom) -> nonfungalOTUs.df
  taxAndFunguild %>%
    filter(!OTUId %in% nonfungalOTUs.df$OTUId) -> tmp
  # also delete OTUs with suspect coverage (probably nonfungal)
  tmp %>%
    filter(coverage > 0.9 | kingdom == "Protist") -> tmp2 # this is the step that accidently removed the oomycetes!!!!!
  taxAndFunguild <- tmp2
  sum(grepl("ITSoo", taxAndFunguild$OTUId))

  # delete OTUs from taxAndFunguild if not found in comm.otu (only found in blanks, mock, or very infrequently)
  tmp <- comm.otu.tmp[, colnames(comm.otu.tmp) %in% taxAndFunguild$OTUId]
  dim(comm.otu.tmp); dim(tmp) # gets rid of more than 700 OTUs
  comm.otu.tmp <- tmp
  sum(grepl("ITSoo", colnames(comm.otu.tmp)))
  
  # reorder taxAndFunguild to make OTU table
  o<-match(colnames(comm.otu.tmp), taxAndFunguild$OTUId)
  o.taxAndFunguild<-taxAndFunguild[o,]
  sum(o.taxAndFunguild$OTUId != colnames(comm.otu.tmp)) #this need to be 0
  
  # only use FUNGuild info with confidence ranking of Probable or Highly Probable
  criteria <- !o.taxAndFunguild$Confidence.Ranking %in% c("Probable","Highly Probable")
  o.taxAndFunguild[criteria, c("Trophic.Mode","Guild")] <-"unclassified"
  
  # select cols 
  o.taxAndFunguild %>%
    select(OTUId, taxonomy, kingdom, phylum, family, genus, species, 
           Trophic.Mode, Guild) -> o.taxAndFunguild
  
  # o.taxAndFunguild %>%
  #   filter(genus == "unclassified") %>%
  #   filter(Trophic.Mode != "unclassified")
  
  #clean Trophic.Mode
  #unique(o.taxAndFunguild$Trophic.Mode)
  
  #clean Guild
  unique(o.taxAndFunguild$Guild)
  o.taxAndFunguild[o.taxAndFunguild$Guild=="NULL","Guild"]<-"unclassified"
  
  #clean oomycetes
  o.taxAndFunguild[o.taxAndFunguild$kingdom=="Protist", c("taxonomy","phylum","family","genus")]<-"unclassified"
  o.taxAndFunguild[o.taxAndFunguild$kingdom=="Protist", c("species")]<-"unclassified_Protist"
  
  #clean taxa
  o.taxAndFunguild %>%
    filter(is.na(genus)) %>%
    filter(Trophic.Mode != "unclassified") -> tmp
  dim(tmp) # why were these matched in FUNGuild?  I'm going to remove the FUNGuild designation
  o.taxAndFunguild[is.na(o.taxAndFunguild$phylum), 'phylum'] <- 'unclassified'
  o.taxAndFunguild[is.na(o.taxAndFunguild$family), 'family'] <- 'unclassified'
  o.taxAndFunguild[is.na(o.taxAndFunguild$genus), 'genus'] <- 'unclassified'
  o.taxAndFunguild[is.na(o.taxAndFunguild$species), 'species'] <- 'unclassified'
  criteria <- o.taxAndFunguild$genus == "unclassified"
  o.taxAndFunguild[criteria, c("Trophic.Mode","Guild")] <- "unclassified"
  
  #clean species
  #species column should have [genus]_sp if the genus is known
  criteria <- o.taxAndFunguild$genus != "unclassified" & o.taxAndFunguild$species == "unclassified"
  o.taxAndFunguild %>%
    mutate(species_fake = paste(genus, "sp", sep = "_")) %>%
    mutate(species_new = ifelse(genus != "unclassified" & species == "unclassified", 
                                species_fake, species)) %>%
    select(-species_fake) %>%
    select(-species) %>%
    rename('species'='species_new') -> o.taxAndFunguild
  
  #fix weird characters in 'Montagnula_aloës'
  o.taxAndFunguild[grep('Montagnula', o.taxAndFunguild$species), "species"] <- "Montagnula_aloes"
  
  #add numbers to repeated names in species to indicate that they are different OTUs
  o.taxAndFunguild %>%
    filter(grepl("_sp", species)) %>%
    group_by(species) %>%
    summarize(n = length(species)) %>%
    filter(n > 1) -> spfake_indx
  spfake_indx
  
  for(i in 1:dim(spfake_indx)[1]){
    curr.sp <- as.character(spfake_indx[i,"species"])
    curr.sp
    sp.vec <- o.taxAndFunguild[o.taxAndFunguild$species == curr.sp, "species"]
    sp.vec
    new.sp.vec <- paste(sp.vec, 1:length(sp.vec), sep="")
    new.sp.vec
    o.taxAndFunguild[o.taxAndFunguild$species == curr.sp, "species"] <- new.sp.vec
  }
  #simplify the OTUId and create an annotated OTUId column using species
  o.taxAndFunguild %>%
    separate(OTUId, into=c("drop","drop1","OTUId_num"), remove = FALSE) %>%
    select(-drop) %>% select(-drop1) %>%
    mutate(OTUId_simp = paste("OTU", OTUId_num, sep="_")) %>%
    select(-OTUId_num) %>%
    mutate(OTUId_ann = ifelse(species == "unclassified",
                              OTUId_simp, species)) %>%
    select(-OTUId_simp) -> o.taxAndFunguild
  
  sum(grepl("ITSoo", o.taxAndFunguild$OTUId))
  
  return(o.taxAndFunguild)
}

#copied from fungal_wood_endophtyes repo
clean_comm<-function(comm.otu.tmp, taxAndFunguild){
  
  # delete OTUs from comm.otu not found in taxAndFunguild
  comm.otu <- comm.otu.tmp[, colnames(comm.otu.tmp) %in% taxAndFunguild$OTUId]
  
  return(comm.otu)
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

load_MicrobeCollection <- function(stemSamples){
  
  # OTU table
  comm.otu <- load_matotu()
  dim(comm.otu) # 3795 OTUs, includes oomycetes
  sum(grepl("ITSoo", colnames(comm.otu)))
  
  # taxon table
  taxAndFunguild <- load_TaxAndFunguild(comm.otu)
  dim(taxAndFunguild) # 3021 OTUs because removed OTUs with suspect coverage, non-fungal?
  
  sum(grepl("ITSoo", taxAndFunguild$OTUId))
  taxAndFunguild %>%
    filter(kingdom == "Protist")
  
  # cleaning
  comm.otu <- clean_comm(comm.otu, taxAndFunguild) # need to also remove these OTUs from the OTU table
  dim(comm.otu) # good. also 3021 OTUs
  
  # sequencing sample look up table
  seqSamples <- load_seqSamples(comm.otu, stemSamples) #create sequence sample meta data table
  dim(seqSamples) # 105 samples (1 failed sequencing)
  seqSamples[!seqSamples$seq_sampName %in% row.names(comm.otu),] # failed sample
  
  microb.data <- list(taxAndFunguild = taxAndFunguild,
       comm.otu = comm.otu,
       seqSamples = seqSamples)
  
  return(microb.data)
  
}

removeRareOTUs <- function(comm.otu){
  
  minSamps<-floor(dim(comm.otu)[1] * .2) #how many samples is 20% of them?
  dat1 <- apply(comm.otu>0, 2, sum) #how many samples does each OTU show up in?
  comm.otu.trimmed <- comm.otu[,dat1>minSamps]
  
  print(paste("Keep", dim(comm.otu.trimmed)[2], "of", dim(comm.otu)[2], "OTUs"))
  return(comm.otu.trimmed)
  
}
