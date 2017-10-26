
Calc_richOTUtype<-function(colNam, grepTerm, taxAndFunguild, comm.otu){
  
  # identify OTUs that match that type
  sub.otus<-taxAndFunguild[grepl(grepTerm,taxAndFunguild[,colNam]),"OTUId"]
  
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

Plot_richOTUtype<-function(otutype.df, valueCol_vec, otutypeNam, spdf){
  
  # summarize by Code
  otutype.df %>%
    group_by(code) %>%
    summarize(mean=mean(sub_rich),
              se=sd(sub_rich)/sqrt(length(sub_rich)),
              upper=mean+se,
              lower=mean-se) -> summ.otutype.df
  
  # join with decay trajectory data
  tmp<-left_join(spdf, summ.otutype.df)
  
  # create xlab
  xlab<-paste(otutypeNam,"richness", sep=" ")
  
  # create plot
  pList<-list()
  for(i in 1:length(valueCol_vec)){
    colnames(tmp)[colnames(tmp) == valueCol_vec[i]]<-"yCol"
    pList[[i]]<-
      ggplot(tmp, aes(x=mean, y=yCol, color=species, shape=size)) + 
      geom_point() + 
      geom_errorbarh(aes(xmin=lower, xmax=upper)) +
      xlab(xlab) + ylab(valueCol_vec[i])
  }
  names(pList)<-valueCol_vec
  
  return(pList)
}

