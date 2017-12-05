HowManyOfThese<-function(otuIDs, taxAndFunguild, comm.otu){
  
  #subset the OTU table
  curr.mat.otu<-comm.otu[,colnames(comm.otu) %in% otuIDs$OTUId]
  
  #get rid of cols and rows without reads
  tmp<-curr.mat.otu[,colSums(curr.mat.otu)!=0] 
  curr.mat.otu.sel<-tmp[rowSums(tmp)!=0,]
  
  #save the mini OTU table and taxon look up table
  results<-list(otu=curr.mat.otu.sel,
                tax=taxAndFunguild[taxAndFunguild$OTUId %in% colnames(curr.mat.otu.sel),c("OTUId","phylum","genus","species")])
  
  return(results)
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

Create_rich_spdf_DF<-function(otutype.df, spdf){
  
  # summarize by Code
  otutype.df %>%
    group_by(code) %>%
    summarize(mean=mean(sub_rich),
              se=sd(sub_rich)/sqrt(length(sub_rich)),
              upper=mean+se,
              lower=mean-se) -> summ.otutype.df
  
  # join with decay trajectory data
  rich.spdf<-left_join(spdf, summ.otutype.df)
  
  return(rich.spdf)
  
}

Create_rich_tr_DF<-function(trait.residuals.list, rich.df){
  
  rich.tr.list<-list()
  for(i in 1:length(trait.residuals.list)){
    curr.trait.residuals<-trait.residuals.list[[i]]
    rich.tr.list[[i]]<-curr.trait.residuals %>% left_join(rich.df)
  }
  names(rich.tr.list)<-names(trait.residuals.list)
  
  return(rich.tr.list)
  
}

Plot_richOTUtype<-function(rich.spdf, valueCol_vec, otutypeNam){
  
  # create xlab
  xlab<-paste(otutypeNam,"richness", sep=" ")
  
  # create plot
  pList<-list()
  for(i in 1:length(valueCol_vec)){
    curr.tmp<-rich.spdf
    colnames(curr.tmp)[colnames(curr.tmp) == valueCol_vec[i]]<-"yCol"
    pList[[i]]<-
      ggplot(curr.tmp, aes(x=mean, y=yCol, color=species, shape=size)) + 
      geom_point() + 
      geom_errorbarh(aes(xmin=lower, xmax=upper)) +
      xlab(xlab) + ylab(valueCol_vec[i])
  }
  names(pList)<-valueCol_vec
  
  return(pList)
}

FitResid_diversity<-function(rich.spdf, trait.residuals){
  
  #merge the trait residuals and diversity dataframes
  rich.spdf %>%
    left_join(trait.residuals) %>%
    filter(!is.na(mean)) -> rich.spdf.tr
  
  #fit models
  mod.r2<-lm(ne.r2~size+mean, data=rich.spdf.tr)
  mod.k<-lm(k~size+mean, data=rich.spdf.tr)
  mod.t70<-lm(t70~size+mean, data=rich.spdf.tr)
  #mod.alpha<-lm(alpha~size+mean, data=rich.spdf.tr)
  
  mod.list<-list(r2=mod.r2, k=mod.k, t70=mod.t70)
  
  return(mod.list)
  
}

FitResid_diversity_byStem<-function(rich.tr.list){
  
  mod.list<-list()
  for(i in 1:length(rich.tr.list)){
    curr.df<-rich.tr.list[[i]]
    mod.list[[i]]<-lm(trait.resid ~ size + sub_rich, data=curr.df)
  }
  names(mod.list)<-names(rich.tr.list)
  
  return(mod.list)
}

PlotResid_diversity<-function(rich.spdf, trait.residuals, xlab){
  
  #merge the trait residuals and diversity dataframes
  rich.spdf %>%
    left_join(trait.residuals) %>%
    filter(!is.na(mean)) -> rich.spdf.tr
  
  #fit models
  p.r2<-ggplot(rich.spdf.tr, aes(x=mean, y=ne.r2, color=species, shape=size)) + geom_point() + xlab(xlab)
  p.k<-ggplot(rich.spdf.tr, aes(x=mean, y=k, color=species, shape=size)) + geom_point() + xlab(xlab)
  p.t70<-ggplot(rich.spdf.tr, aes(x=mean, y=t70, color=species, shape=size)) + geom_point() + xlab(xlab)
  p.alpha<-ggplot(rich.spdf.tr, aes(x=mean, y=alpha, color=species, shape=size)) + geom_point() + xlab(xlab)
  
  p.list<-list(r2=p.r2, k=p.k, t70=p.t70, alpha=p.alpha)
  
  return(p.list)
  
}

CalcNumPairsPresent<-function(presentOTUs, corPairs.df){
  
  # define the number of pairs in the correlation list
  numCorPairs<-dim(corPairs.df)[1]
  pairPresent.list<-list()
  
  # identify whether each set of OTUs in a pair are present in the current community
  for(i in 1:numCorPairs){
    curr.pair<-corPairs.df[i, c("otu1","otu2")]
    curr.result<-curr.pair$otu1 %in% presentOTUs & curr.pair$otu2 %in% presentOTUs
    pairPresent.list[[i]]<-data.frame(curr.pair, pairPresent=curr.result)
  }
  pairPresent.list %>% 
    lapply(as.data.frame) %>% 
    bind_rows() -> pairPresent.df
  
  #summarize the number of pairs present in the current community
  numPairsPresent<-sum(pairPresent.df$pairPresent)
  
  return(numPairsPresent)
}

IDcorrelatedOTUsinEachSampl<-function(taxaAndFunguild, comm.otu){
  
  #1. load residual signif correlated OTUs from boral Env+LV model
  residCor<-load_boralResidCors(taxAndFunguild)
  
  #2. identify signif positively & negatively correlated OTU pairs
  pos.residCor<- residCor %>% filter(sig.residCor > 0)
  neg.residCor<- residCor %>% filter(sig.residCor < 0)
  
  #3. convert OTU table to presence/absence
  comm.otu.pa<-comm.otu
  comm.otu.pa[comm.otu.pa > 0] <- 1 # convert to presence/absence
  
  #4. for each community, identify the OTUs present
  numSeqSamps<-dim(comm.otu)[1]
  presentOTUs.list<-list()
  for(i in 1:numSeqSamps){
    curr.comm<-comm.otu.pa[i,]
    presentOTUs.list[[i]]<-names(curr.comm[curr.comm == 1])
  }
  names(presentOTUs.list)<-row.names(comm.otu)
  #presentOTUs.list
  
  #5. for each community, look up how many pairs are present
  #the lapply lines are super slow!!!!!!!
  
  #total number of OTU pairs possible in a sample
  totalPossiblePairs.list<-list()
  for(i in 1:length(presentOTUs.list)){
    curr.samp<-presentOTUs.list[[i]]
    tmp<-t(combn(curr.samp,2))
    totalPossiblePairs.list[[i]]<-dim(tmp)[1]
  }
  totalPossiblePairs.list %>%
    lapply(as.data.frame) %>%
    bind_rows() %>%
    bind_cols(data.frame(names(presentOTUs.list))) -> totalPossiblePairs.df
  colnames(totalPossiblePairs.df)<-c("numPairs_allPossible","seq_sampName")
  
  #positively correlated pairs
  pairsPresent.pos.list<-lapply(presentOTUs.list, CalcNumPairsPresent, corPairs.df=pos.residCor)
  pairsPresent.pos.list %>%
    lapply(as.data.frame) %>%
    bind_rows() %>%
    bind_cols(data.frame(names(presentOTUs.list))) -> pairsPresent.pos.df
  colnames(pairsPresent.pos.df)<-c("numPairs_posCorr","seq_sampName")
  
  #negatively correlated pairs
  pairsPresent.neg.list<-lapply(presentOTUs.list, CalcNumPairsPresent,corPairs.df=neg.residCor)
  pairsPresent.neg.list %>%
    lapply(as.data.frame) %>%
    bind_rows() %>%
    bind_cols(data.frame(names(presentOTUs.list)))-> pairsPresent.neg.df
  colnames(pairsPresent.neg.df)<-c("numPairs_negCorr","seq_sampName")
  
  #6. merge the total, positive and negative pairs data into 1 df
  totalPossiblePairs.df %>%
    left_join(pairsPresent.pos.df) %>%
    left_join(pairsPresent.neg.df) -> pairsPresent.df
  
  #7. calculate the freq of positive and negative pairs in each sample
  pairsPresent.df %>%
    mutate(percPairs_posCorr = (numPairs_posCorr / numPairs_allPossible) *100) %>%
    mutate(percPairs_negCorr = (numPairs_negCorr / numPairs_allPossible) *100) -> pairsPresent.df
  
  # select and reorder columns
  pairsPresent.df %>%
    select(seq_sampName, percPairs_posCorr, percPairs_negCorr) -> pairsPresent.df
  
  write.csv(pairsPresent.df, file="derived_data/pairsPresent.csv")
  
}

Plot_signifCor<-function(sign, valueCol_vec, pairsPresent.df, spdf, seqSamples){
  
  # 1. select the correlation sign you'd like to look at
  if(sign=='pos'){
    pairsPresent.df %>%
      select(seq_sampName, percPairs_posCorr) %>%
      rename('percPairs'='percPairs_posCorr') -> pairsPresent.sel.df
  }
  if(sign=='neg'){
    pairsPresent.df %>%
      select(seq_sampName, percPairs_negCorr) %>%
      rename('percPairs'='percPairs_negCorr') -> pairsPresent.sel.df
  }
  
  # 2. summarize by Code
  pairsPresent.sel.df %>%
    left_join(seqSamples) %>%
    group_by(code) %>%
    summarize(mean=mean(percPairs),
              se=sd(percPairs)/sqrt(length(percPairs)),
              upper=mean+se,
              lower=mean-se) -> summ.pairsPresent.sel.df
  
  # 3. join with decay trajectory data
  tmp<-left_join(spdf, summ.pairsPresent.sel.df)
  
  # 4. create xlab
  if(sign=='pos'){
    xlab<-"Positively correlated OTU pairs (%)"
  }
  if(sign=='neg'){
    xlab<-"Negatively correlated OTU pairs (%)"
  }
  
  # 5. create plots
  pList<-list()
  for(i in 1:length(valueCol_vec)){
    curr.tmp<-tmp
    colnames(curr.tmp)[colnames(curr.tmp) == valueCol_vec[i]]<-"yCol"
    pList[[i]]<-
      ggplot(curr.tmp, aes(x=mean, y=yCol, color=species, shape=size)) + 
      geom_point() + 
      geom_errorbarh(aes(xmin=lower, xmax=upper)) +
      xlab(xlab) + ylab(valueCol_vec[i])
  }
  names(pList)<-valueCol_vec
  
  return(pList)
}

reformatMatrix<-function(commmat){
  newmat<-matrix(as.numeric(as.matrix(commmat)), ncol=dim(commmat)[2], nrow=dim(commmat)[1])
  return(newmat)
}

CreateCommPMRpair<-function(timePoint, comm.mat, pmr.byStem.df.w){
  
  #make a dataframe using the current time point's pmr and remove NAs
  pmr.byStem.df.w %>%
    select_("codeStem",timePoint) %>%
    rename_("curr.time"=timePoint) %>%
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

CreateCommTraitResidpair<-function(timePoint, comm.mat, traitResid){
  
  #there should be no NAs in the current time point's traitResid
  traitResid %>%
    filter(!is.na(trait.resid)) -> traitResid.noNAs
  
  #subset the community matrix using these unique codeStems
  curr.comm<-comm.mat[row.names(comm.mat) %in% traitResid.noNAs$codeStem,]
  curr.comm<-curr.comm[,colSums(curr.comm) != 0] #get rid of any OTU columns with all 0s
  
  #get rid of traitResid rows for which there is missing community data
  traitResid.noNAs %>%
    filter(codeStem %in% row.names(curr.comm)) -> curr.traitResid
  
  #make sure the row orders match
  ord<-match(row.names(curr.comm), curr.traitResid$codeStem)
  curr.traitResid<-curr.traitResid[ord,]
  
  #reformat the community matrix
  curr.comm.reform<-reformatMatrix(curr.comm)
  row.names(curr.comm.reform)<-row.names(curr.comm)
  colnames(curr.comm.reform)<-colnames(curr.comm)
  
  modelDat.list<-list(traitresid=curr.traitResid, comm=curr.comm.reform)
  
  return(modelDat.list)
  
}

CreateCommTraitResidpair_code<-function(comm.mat, trait.residuals){
  
  trait.residuals
  comm.mat<-meanOTUabund.trim
  
  #code code as a character
  trait.residuals$code<-as.character(trait.residuals$code)
  
  #look for missing codes in comm.mat that exist in tr
  trait.residuals[!trait.residuals$code %in% row.names(comm.mat),"code"] #none missing
  
  #look for missing codes in tr that exist in comm.mat
  row.names(comm.mat)[!row.names(comm.mat) %in% trait.residuals$code] #missing olst from trait.residuals
  
  #trim tr
  trait.residuals.trim<-trait.residuals[trait.residuals$code %in% row.names(comm.mat),] #trim
  
  #trim comm.mat
  comm.mat.trim<-comm.mat[row.names(comm.mat) %in% trait.residuals.trim$code,] #trim
  
  #align rows
  ord<-match(row.names(comm.mat.trim), trait.residuals.trim$code)
  trait.residuals.trim.o<-trait.residuals.trim[ord,]
  rowOrderWarn<-sum(trait.residuals.trim.o$code != row.names(comm.mat.trim)) #this needs to by 0
  rowOrderWarn
  
  #make sure there are no empty OTU cols
  comm.mat.trim<-comm.mat.trim[,colSums(comm.mat.trim)!=0]
  
  if(rowOrderWarn !=0){
    modelDat.list<-"rows are out of order"
  }else{
    modelDat.list<-list(traitresid=trait.residuals.trim.o, comm=comm.mat.trim)
  }

  return(modelDat.list)
  
}

CreateTraitPMRpair<-function(timePoint, traits.codeStem, traits.mean, pmr.byStem.df.w){
  
  #make a dataframe using the current time point's pmr and remove NAs
  pmr.byStem.df.w %>%
    select_("codeStem",timePoint) %>%
    rename_("curr.time"=timePoint) %>%
    filter(!is.na(curr.time)) -> pmr.noNAs
  
  #subset the trait matrix using these unique codeStems
  traits.codeStem %>%
    filter(codeStem %in% pmr.noNAs$codeStem) -> curr.traits
  
  #stem-level trait set
  curr.traits %>%
    filter(!is.na(waterperc) & !is.na(P) & !is.na(K) & !is.na(Ca) & !is.na(Mn) & !is.na(Fe) & !is.na(Zn) & !is.na(N) & !is.na(C)) -> curr.traits
  
  #add species-level traits
  traits.mean %>%
    select(code, barkthick, density) %>%
    rename('barkthick_smspp'='barkthick',
           'density_smspp'='density')-> select.traits.mean
  curr.traits %>%
    left_join(select.traits.mean) -> curr.traits
  
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
