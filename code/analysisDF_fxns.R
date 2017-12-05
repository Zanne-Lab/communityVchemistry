
CreateTraitPMRpair<-function(respVar, traits.stem, traits.code, pmr_byStem){
  
  #make a dataframe using the current time point's pmr and remove NAs
  pmr_byStem %>%
    select_("codeStem", respVar) %>%
    rename_("curr.pmr" = respVar) %>%
    filter(!is.na(curr.pmr)) -> pmr.noNAs
  
  #subset the trait matrix using these unique codeStems
  traits.stem %>%
    filter(codeStem %in% pmr.noNAs$codeStem) -> curr.traits
  
  #stem-level trait set
  curr.traits %>%
    filter(!is.na(waterperc) & !is.na(P) & !is.na(K) & !is.na(Ca) & !is.na(Mn) & !is.na(Fe) & !is.na(Zn) & !is.na(N) & !is.na(C)) -> curr.traits
  
  #add species-level traits
  traits.code %>%
    select(code, barkthick, density) %>%
    rename('barkthick_smspp'='barkthick',
           'density_smspp'='density')-> select.traits.code
  curr.traits %>%
    left_join(select.traits.code) -> curr.traits
  
  #get rid of pmr rows for which there is missing trait data
  pmr.noNAs %>%
    filter(codeStem %in% curr.traits$codeStem) -> curr.pmr

  #merge the dataframes
  curr.df<-left_join(curr.pmr, curr.traits) 
  
  #add code and species and size
  curr.df<-separate(curr.df, col=codeStem, into=c("code","Stem"), sep=4, remove=FALSE)
  curr.df$species<-tolower(curr.df$code)
  curr.df$size<-"large"
  curr.df[curr.df$code == tolower(curr.df$code),"size"]<-"small"
  
  return(curr.df)
  
}

AverageOTUabund_byCode<-function(comm.otu, seqSamples){
  
  #loop through each code and average the OTU abundances
  CODE<-unique(seqSamples$code)
  meanOTUabund.list<-list()
  sdOTUabund.list<-list()
  for(i in 1:length(CODE)){
    #identify the seq_sampNames to aggregate
    curr.seqSamps<-seqSamples[seqSamples$code==CODE[i],]$seq_sampName 
    #filter 
    curr.comm.otu<-comm.otu[row.names(comm.otu) %in% curr.seqSamps,]
    #summarize
    meanOTUabund.list[[i]]<-apply(curr.comm.otu, 2, mean, na.rm=TRUE)
    sdOTUabund.list[[i]]<-apply(curr.comm.otu, 2, sd, na.rm=TRUE)
  }
  names(meanOTUabund.list)<-CODE
  names(sdOTUabund.list)<-CODE
  meanOTUabund.list %>% combine() %>% as.data.frame %>% t() -> meanOTUabund
  sdOTUabund.list %>% combine() %>% as.data.frame %>% t() -> sdOTUabund
  
  #get rid of empty OTU cols
  #sum(colSums(meanOTUabund)!=0)
  meanOTUabund.trim<-meanOTUabund[,colSums(meanOTUabund)!=0]
  
  return(meanOTUabund.trim)
  
}

fitNcrossval_WAPLS<-function(curr.comm, curr.respVar){
  fit<-WAPLS(y=curr.comm, x=curr.respVar)
  fit.cv<-crossval(fit, cv.method="loo")
  return(fit.cv)
}

reformatMatrix<-function(commmat){
  newmat<-matrix(as.numeric(as.matrix(commmat)), ncol=dim(commmat)[2], nrow=dim(commmat)[1])
  return(newmat)
}

CreateCommPMRpair<-function(respVar, comm.mat, pmr_byStem){
  
  #make a dataframe using the current time point's pmr and remove NAs
  pmr_byStem %>%
    select_("codeStem", respVar) %>%
    rename_("curr.time"= respVar) %>%
    filter(!is.na(curr.time)) -> pmr.noNAs
  
  #subset the community matrix using these unique codeStems
  curr.comm<-comm.mat[row.names(comm.mat) %in% pmr.noNAs$codeStem,]
  curr.comm<-curr.comm[,colSums(curr.comm) != 0] #get rid of any OTU columns with all 0s
  
  #get rid of pmr rows for which there is missing community data
  pmr.noNAs %>%
    filter(codeStem %in% row.names(curr.comm)) -> curr.pmr
  
  #make sure the row orders match
  ord<-match(row.names(curr.comm), curr.pmr$codeStem)
  curr.pmr<-curr.pmr[ord,]
  
  #reformat the community matrix
  curr.comm.reform<-reformatMatrix(curr.comm)
  row.names(curr.comm.reform)<-row.names(curr.comm)
  colnames(curr.comm.reform)<-colnames(curr.comm)
  
  modelDat.list<-list(pmr=curr.pmr, comm=curr.comm.reform)
  
  return(modelDat.list)
  
}

CreateCommTraitResidpair<-function(respVar, comm.mat, traitResiduals, sampleName){
  
  #select the current respVar's traitResid
  traitResiduals %>%
    filter(resp==respVar) %>%
    filter(!is.na(resid)) %>%
    rename_('sampleName'=sampleName) -> traitResid.noNAs
  
  #subset the community matrix using these unique codeStems
  curr.comm<-comm.mat[row.names(comm.mat) %in% traitResid.noNAs$sampleName,]
  curr.comm<-curr.comm[,colSums(curr.comm) != 0] #get rid of any OTU columns with all 0s
  
  #get rid of traitResid rows for which there is missing community data
  traitResid.noNAs %>%
    filter(sampleName %in% row.names(curr.comm)) -> curr.traitResid
  
  #make sure the row orders match
  ord<-match(row.names(curr.comm), curr.traitResid$sampleName)
  curr.traitResid<-curr.traitResid[ord,]
  
  #reformat the community matrix
  curr.comm.reform<-reformatMatrix(curr.comm)
  row.names(curr.comm.reform)<-row.names(curr.comm)
  colnames(curr.comm.reform)<-colnames(curr.comm)
  
  modelDat.list<-list(traitresid=curr.traitResid, comm=curr.comm.reform)
  
  return(modelDat.list)
  
}

Calc_richOTU<-function(taxAndFunguild, comm.otu){
  
  # turn the OTU table into presence/absence
  comm.otu.pa<-comm.otu
  comm.otu.pa[comm.otu.pa > 0] <- 1 # convert to presence/absence
  
  # calculate the richness of OTUs per sample and save it in a df
  comm.rich<-rowSums(comm.otu.pa)
  comm.rich.df<-data.frame(seq_sampName=names(comm.rich), sub_rich=comm.rich)
  
  # use the seqSample to add code, size, species
  comm.rich.df<-left_join(comm.rich.df, seqSamples)
  
  return(comm.rich.df)
}

Calc_H.OTU<-function(taxAndFunguild, comm.otu){
  
  #calculate the shannon diversity
  require(vegan)
  comm.H<-diversity(comm.otu, index="shannon")
  comm.H.df<-data.frame(seq_sampName=names(comm.H), sub_rich=comm.H)
  
  # use the seqSample to add code, size, species
  comm.H.df<-left_join(comm.H.df, seqSamples)
  
  return(comm.H.df)
}

Calc_richOTUtype<-function(colNam, grepTerm, taxAndFunguild, comm.otu){
  
  # identify OTUs that match that type
  tf.col<-as.matrix(taxAndFunguild[,colNam])[,1]
  sub.otus<-taxAndFunguild[grepl(grepTerm,tf.col),"OTUId"]
  
  # subset the OTU table by the select otus
  sub.otu<-comm.otu[,colnames(comm.otu) %in% sub.otus$OTUId]
  
  # turn the OTU table into presence/absence
  sub.otu.pa<-sub.otu
  sub.otu.pa[sub.otu.pa > 0] <- 1 # convert to presence/absence
  
  # calculate the richness of select OTUs per sample and save it in a df
  sub.rich<-rowSums(sub.otu.pa)
  sub.rich.df<-data.frame(seq_sampName=names(sub.rich), sub_rich=sub.rich)
  
  # use the seqSample to add code, size, species
  sub.rich.df<-left_join(sub.rich.df, seqSamples)
  
  return(sub.rich.df)
}

CreateRichDecayfits_df<-function(otutype.df, decayfits){
  
  # summarize by Code
  otutype.df %>%
    group_by(code) %>%
    summarize(mean=mean(sub_rich),
              se=sd(sub_rich)/sqrt(length(sub_rich)),
              upper=mean+se,
              lower=mean-se) -> summ.otutype.df
  
  # join with decay trajectory data
  rich.decayfits<-left_join(decayfits, summ.otutype.df)
  
  return(rich.decayfits)
  
}