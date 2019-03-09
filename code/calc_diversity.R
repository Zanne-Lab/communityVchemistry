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
    temp<-apply(curr.comm.otu, 2, mean, na.rm=TRUE)
    temp[is.nan(temp)] <- NA
    temp2<-t(data.frame(temp))
    temp3<-data.frame(code=CODE[i], temp2)
    meanOTUabund.list[[i]]<-temp3
    temp<-apply(curr.comm.otu, 2, sd, na.rm=TRUE)
    temp[is.nan(temp)] <- NA
    temp2<-t(data.frame(temp))
    temp3<-data.frame(code=CODE[i], temp2)
    sdOTUabund.list[[i]]<-temp3
  }
  meanOTUabund <- bind_rows(meanOTUabund.list) 
  #warnings() #these are about binding characters and factors
  
  #make the code col into rownames
  tmp<-data.frame(meanOTUabund[,-1])
  row.names(tmp)<-meanOTUabund$code
  
  #get rid of empty OTU cols
  meanOTUabund.trim<-tmp[,colSums(tmp, na.rm = T)!=0]
  
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
  sub.otu<-comm.otu[,colnames(comm.otu) %in% sub.otus]
  
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
