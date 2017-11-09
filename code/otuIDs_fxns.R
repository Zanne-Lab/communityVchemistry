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
  mod.alpha<-lm(alpha~size+mean, data=rich.spdf.tr)
  
  mod.list<-list(r2=mod.r2, k=mod.k, t70=mod.t70, alpha=mod.alpha)
  
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
