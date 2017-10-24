
#create a sample Tab for the sequence samples --- probs want to move this to load_fxns and not use mass.data
CreateSeqSampTab<-function(mass.data){
  sampTab<-unique(mass.data[,c("Species","size")])
  sampTab<-rename(sampTab, "code"="Species")
  seq_sampName<-row.names(comm.otu)
  indx<-data.frame(seq_sampName=seq_sampName, code=substr(seq_sampName, 1, 4))
  sampTab<-left_join(indx, sampTab)
  
  return(sampTab)
}

extract_uniquePairDists<-function(dist.mat){
  x<-dist.mat
  rowCol <- expand.grid(rownames(x), colnames(x))
  labs <- rowCol[as.vector(upper.tri(x,diag=F)),]
  df <- cbind(labs, x[upper.tri(x,diag=F)])
  colnames(df) <- c("sp1","sp2","dist")
  
  return(df)
}

Calc_commDists<-function(sampTab, comm.otu){
  
  #1. create sample table for the community data
  #sampTab
  
  #2. calc dist, make it long format
  dist <- vegdist(as.matrix(comm.otu), method = "bray", upper=TRUE, diag=TRUE)  #calculate jaccard index
  dist.mat<-as.matrix(dist)
  dist.df<-extract_uniquePairDists(dist.mat) #make the distance matrix long and add metadata
  dist.df<-rename(dist.df, "sampID1"="sp1", "sampID2"="sp2")
  
  #3. annotate distances with sample info
  indx<-rename(sampTab, "sampID"="seq_sampName")
  left_join(dist.df, indx, by=c("sampID1" = "sampID")) %>%
    rename("code1"="code",
           "size1"="size") -> dist.df1 # use the index to add info for 1st sample
  left_join(dist.df1, indx, by=c("sampID2" = "sampID")) %>%
    rename("code2"="code",
           "size2"="size") -> dist.df2 # use the index to add info for the 2nd sample
  
  #4. remove distances between size classes
  comm.dist<-subset(dist.df2, size1 == size2)
  
  return(comm.dist)
  
}

#--- probs want to move this to load_fxns and not use mass.data
AddCodeID<-function(sampTab){
  
  sampTab<-mutate(sampTab, species=tolower(code))
  indx<-unique(sampTab[,c("code","species","size")])
  spdf<-left_join(spdf, indx) #put code into the dataframe
  
  #fix an NA code 
  spdf[is.na(spdf$code),"code"]<-"eusc"
  
  return(spdf)
  
}

Calc_decayParamDiffs<-function(valueCol, spdf, sampTab){
  
  # do difference calculations
  spdf<-rename_(spdf, "select.decayParam"=valueCol)
  v <- spdf$select.decayParam
  z <- outer(v,v,'-') # sp1 - sp2 = dist
  colnames(z)<-spdf$code
  row.names(z)<-spdf$code
  dist.df<-extract_uniquePairDists(z) #make the distance matrix long and add metadata
  dist.df<-rename(dist.df, "code1"="sp1", "code2"="sp2")
  
  # add back species + size identifiers
  indx<-unique(sampTab[,c("code","size")])
  left_join(dist.df, indx, by=c("code1" = "code")) %>%
    rename("size1"="size") -> dist.df1 # use the index to add info for 1st sample
  left_join(dist.df1, indx, by=c("code2" = "code")) %>%
    rename("size2"="size") -> dist.df2 # use the index to add info for the 2nd sample
  
  #4. remove distances between size classes
  decayparam.dist<-subset(dist.df2, size1 == size2)
  
  return(decayparam.dist)
  
}

SummarizeCommDist_byCodePair<-function(comm.dist){
  
  comm.dist$codePair<-paste(comm.dist$code1, comm.dist$code2, sep="_")
  comm.df<-comm.dist[,c("codePair","dist")]
  comm.df<-rename(comm.df, "comm_dist"="dist")
  summ.comm_dist<-group_by(comm.df, codePair) %>%
    summarize(mean=mean(comm_dist),
              se=sd(comm_dist)/sqrt(length(comm_dist)),
              lower=mean-se,
              upper=mean+se)
  summ.comm_dist<-rename(summ.comm_dist, "mean_comm_dist"="mean",
                         "lower_comm_dist"="lower",
                         "upper_comm_dist"="upper")
  
  return(summ.comm_dist)
}

MergeCommNDecaydists_byCodePair<-function(decayparam.dist, summ.comm_dist){
  
  # elongate aic.dist by adding the forward and reverse combinations of codePair
  decayparam.dist$codePair<-paste(decayparam.dist$code1, decayparam.dist$code2, sep="_")
  decayparam.dist$codePair_rev<-paste(decayparam.dist$code2, decayparam.dist$code1, sep="_")
  decayparam.df<-decayparam.dist[,c("codePair","codePair_rev","dist")]
  decayparam.df<-rename(decayparam.df, "decayparam_dist"="dist")
  decayparam.df_forward<-decayparam.df[,c("codePair","decayparam_dist")]
  decayparam.df_reverse<-decayparam.df[,c("codePair_rev","decayparam_dist")]
  uniqCodePairs<-unique(c(decayparam.df$codePair, decayparam.df$codePair_rev))
  
  #why are all these code pairs missing from aic.df?
  missingCodePairs<-summ.comm_dist$codePair[!unique(summ.comm_dist$codePair) %in% uniqCodePairs]
  # df<-data.frame(codePair=missingCodePairs)
  #  df %>%
  #   separate(codePair, into=c("code1","code2"), remove=FALSE) %>%
  #   filter(code1 != code2) #1. aic does not include code pairs within the same species+size
  
  summ.comm_dist %>%
    left_join(decayparam.df_forward) %>% #join with the forward versions of decayparam.df$codePair
    rename("decayparam_dist_forward"="decayparam_dist") %>%
    left_join(decayparam.df_reverse, by=c("codePair"="codePair_rev")) %>%
    rename("decayparam_dist_reverse"="decayparam_dist") -> join.dist
  
  filter(join.dist, is.na(decayparam_dist_forward) & !is.na(decayparam_dist_reverse)) %>%
    mutate(decayparam_dist=decayparam_dist_reverse) -> reverse.rows
  filter(join.dist, !is.na(decayparam_dist_forward) & is.na(decayparam_dist_reverse)) %>%
    mutate(decayparam_dist=decayparam_dist_forward) -> forward.rows
  join.dist<-bind_rows(reverse.rows, forward.rows)
  join.dist<-join.dist[,c("codePair","mean_comm_dist","lower_comm_dist","upper_comm_dist","decayparam_dist")]

  #add size back
  join.dist$size<-"large"
  join.dist[tolower(join.dist$codePair) == join.dist$codePair,"size"]<-"small"
  
  return(join.dist)
  
}

