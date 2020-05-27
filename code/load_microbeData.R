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

# used inside of load_matotu()
add_oomycetes <- function(fung.otu){

  # read in OTU table (uclust output) and convert to matrix (rows=samples, columns=OTUs)
  # data.otu <- read.csv('data/sequencing_T0/OTUtable_oomycetes_20171020.csv', row.names=1)
  data.otu <- read.csv('data/sequencing_T0/DP16_T0_ITSoo_20200527_OTUtable.csv', row.names=1)
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

#modified from fungal_wood_endophtyes repo
load_TaxAndFunguild <- function(){
  
  # load fungal OTU info
  funguild <-read.delim('data/sequencing_T0/DP16_funguild.txt', stringsAsFactors = F)
  tax <-read.delim('data/sequencing_T0/DP16_tax.txt', stringsAsFactors = F)
  
  # merge the dataframes by OTUId
  colnames(tax)[1] <- "OTUId" #make this match the column name in funguild
  tax %>%
    left_join(funguild) -> taxAndFunguild
  
  return(taxAndFunguild)
}
  
clean_TaxTab <- function(taxAndFunguild){
  
  # select cols
  taxAndFunguild %>%
    select(OTUId, taxonomy, kingdom, phylum, class, order, family, genus, species, 
           Trophic.Mode, Guild, Confidence.Ranking) -> tmp
  
  # only use FUNGuild info with confidence ranking of Probable or Highly Probable
  tmp %>%
    mutate(Trophic.Mode = ifelse(Trophic.Mode %in% c("-","NULL"), "unclassified", Trophic.Mode)) %>%
    mutate(Guild = ifelse(Guild %in% c("-","NULL"), "unclassified", Guild)) %>%
    mutate(Trophic.Mode = ifelse(Confidence.Ranking %in% c("Probable","Highly Probable"), Trophic.Mode, "unclassified")) %>%
    mutate(Guild = ifelse(Confidence.Ranking %in% c("Probable","Highly Probable"), Guild, "unclassified")) -> tmp
  #unique(tmp$Trophic.Mode)
  #unique(tmp$Guild)
  
  # is there an issue where a FUNGuild assignment hasn't been applied across everything in that genus?
  tmp %>%
    group_by(genus) %>%
    summarize(uniq.TM = paste(unique(Trophic.Mode), collapse = "_"),
              n.uniq.TM = length(unique(Trophic.Mode))) %>%
    filter(n.uniq.TM > 1) # huh. there are a lot of these for some reason.
  # apply trophic assignment for 1 genus to all
  taxAndFunguild <- tmp
  GENUS <- unique(taxAndFunguild$genus)
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
  
  # check that there are no OTUs that are unclassified to genus but mysteriously have a FUNGuild assignment
  taxAndFunguild %>%
    filter(is.na(genus)) %>%
    filter(Trophic.Mode != "unclassified")
  # why were these matched in FUNGuild?  I'm going to remove the FUNGuild designation
  taxAndFunguild[is.na(taxAndFunguild$genus), c("Trophic.Mode","Guild")] <- "unclassified"
  
  # clean taxonomy data (replace NAs with "unclassified")
  taxAndFunguild %>%
    mutate(phylum = ifelse(is.na(phylum), "unclassified", phylum)) %>%
    mutate(family = ifelse(is.na(class), "unclassified", class)) %>%
    mutate(family = ifelse(is.na(order), "unclassified", order)) %>%
    mutate(family = ifelse(is.na(family), "unclassified", family)) %>%
    mutate(genus = ifelse(is.na(genus), "unclassified", genus)) %>%
    mutate(species = ifelse(is.na(species), "unclassified", species)) -> taxAndFunguild
  
  #clean oomycetes
  taxAndFunguild %>%
    filter(kingdom == "Protist")
  taxAndFunguild[taxAndFunguild$kingdom == "Protist", c("taxonomy","species")] <- "unclassified_Protist"
  
  #clean species
  #species column should have [genus]_sp if the genus is known
  taxAndFunguild %>%
    mutate(species_fake = paste(genus, "sp", sep = "_")) %>%
    mutate(species_new = ifelse(genus != "unclassified" & species == "unclassified", 
                                species_fake, species)) %>%
    select(-species_fake) %>%
    select(-species) %>%
    rename('species'='species_new') -> taxAndFunguild
  
  #fix weird characters in 'Montagnula_aloës'
  taxAndFunguild[grep('Montagnula', taxAndFunguild$species), "species"] <- "Montagnula_aloes"
  
  #add numbers to repeated names in species to indicate that they are different OTUs
  taxAndFunguild %>%
    filter(grepl("_sp", species)) %>%
    group_by(species) %>%
    summarize(n = length(species)) %>%
    filter(n > 1) -> spfake_indx
  spfake_indx
  for(i in 1:dim(spfake_indx)[1]){
    curr.sp <- as.character(spfake_indx[i,"species"])
    sp.vec <- taxAndFunguild[taxAndFunguild$species == curr.sp, "species"]
    new.sp.vec <- paste(sp.vec, 1:length(sp.vec), sep="")
    taxAndFunguild[taxAndFunguild$species == curr.sp, "species"] <- new.sp.vec
  }
  #simplify the OTUId and create an annotated OTUId column using species
  taxAndFunguild %>%
    separate(OTUId, into=c("drop","drop1","OTUId_num"), remove = FALSE) %>%
    select(-drop) %>% select(-drop1) %>%
    mutate(OTUId_simp = paste("OTU", OTUId_num, sep="_")) %>%
    select(-OTUId_num) %>%
    mutate(OTUId_ann = ifelse(species == "unclassified",
                              OTUId_simp, species)) %>%
    select(-OTUId_simp) -> taxAndFunguild
  
  sum(grepl("ITSoo", taxAndFunguild$OTUId)) # oomycetes are still here
  dim(taxAndFunguild) # 3301 OTUs
  
  return(taxAndFunguild)
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

#now works! and the oomycetes are in there!
load_MicrobeCollection <- function(stemSamples){
  
  # OTU table
  comm.otu <- load_matotu()
  dim(comm.otu) # 3795 OTUs, includes oomycetes
  sum(grepl("ITSoo", colnames(comm.otu))) # 13 oomycete OTUs
  
  # taxon table
  taxAndFunguild <- load_TaxAndFunguild()
  dim(taxAndFunguild) # info for 5284 OTUs
  
  # add oomycete info to the taxon table
  sum(grepl("ITSoo", taxAndFunguild$OTUId)) # currently there is not oomycete information in the taxAndFunguild table
  # add oomycete OTUs as rows
  ooOTUs <- colnames(comm.otu)[grepl("ITSoo", colnames(comm.otu))]
  oo.df <- data.frame(OTUId=ooOTUs, kingdom="Protist")
  taxAndFunguild <- bind_rows(taxAndFunguild, oo.df) #warning is ok
  sum(grepl("ITSoo", taxAndFunguild$OTUId))
  
  # identify OTUs in community matrix that are not found in the taxon table -- remove these.. likely plant OTUs
  sum(colnames(comm.otu) %in% taxAndFunguild$OTUId) #3301 OTUs shared
  sum(!colnames(comm.otu) %in% taxAndFunguild$OTUId) #494 OTUs in community matrix that are not in taxon matrix
  rm.OTUs <- colnames(comm.otu)[!colnames(comm.otu) %in% taxAndFunguild$OTUId]
  rm.OTUs # remove these 494 OTUs since they are probably non-fungal and non-oomycete
  comm.otu <- comm.otu[,!colnames(comm.otu) %in% rm.OTUs]
  dim(comm.otu) # 3301 OTUs, includes oomycetes
  sum(grepl("ITSoo", colnames(comm.otu))) # 13 oomycete OTUs
  
  # remove OTUs from the taxon table that are not relevant (i.e. only found in blanks, mock, or very infrequently)
  sum(!taxAndFunguild$OTUId %in% colnames(comm.otu)) # 1996 OTUs in taxon matrix that are not in the community matrix
  rm.OTUs <- taxAndFunguild[!taxAndFunguild$OTUId %in% colnames(comm.otu),"OTUId"]
  rm.OTUs # remove these 1996 OTUs from the taxon table since they are no longer in the cleaned community OTU matrix
  taxAndFunguild <- taxAndFunguild[!taxAndFunguild$OTUId %in% rm.OTUs,]
  dim(taxAndFunguild) # 3301 OTUs, includes oomycetes
  sum(grepl("ITSoo", taxAndFunguild$OTUId)) # 13 oomycete OTUs
  
  # clean the taxon table
  taxAndFunguild <- clean_TaxTab(taxAndFunguild = taxAndFunguild)
  dim(taxAndFunguild)
  
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
