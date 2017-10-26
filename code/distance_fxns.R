
extract_uniquePairDists<-function(dist.mat){
  x<-dist.mat
  rowCol <- expand.grid(rownames(x), colnames(x))
  labs <- rowCol[as.vector(upper.tri(x,diag=F)),]
  df <- cbind(labs, x[upper.tri(x,diag=F)])
  colnames(df) <- c("sp1","sp2","dist")
  return(df)
}

Calc_commDists<-function(seqSamples, comm.otu, distType){
  
  # calc dist
  if(distType=="bray"){
    dist <- vegdist(as.matrix(comm.otu), method = "bray", upper=TRUE, diag=TRUE)
    dist.mat<-as.matrix(dist)
  }
  
  if(distType=="jaccard"){
    comm.otu.pa<-comm.otu
    comm.otu.pa[comm.otu.pa > 0] <- 1 # convert to presence/absence
    dist <- vegdist(as.matrix(comm.otu.pa), method = "jaccard", binary=TRUE, upper=TRUE, diag=TRUE)
    dist.mat<-as.matrix(dist)
  }
  
  if(distType=="richness"){
    comm.otu.pa<-comm.otu
    comm.otu.pa[comm.otu.pa > 0] <- 1 # convert to presence/absence
    richness<-rowSums(comm.otu.pa)
    z <- outer(richness,richness,'-') # sp1 - sp2 = dist
    dist.mat<-abs(z)
  }
  
  # make distances into long format
  dist.df<-extract_uniquePairDists(dist.mat) #make the distance matrix long and add metadata
  dist.df<-rename(dist.df, "sampID1"="sp1", "sampID2"="sp2")
  
  # annotate distances with sample info
  seqSamples %>%
    select(seq_sampName, code, size) %>%
    rename("sampID"="seq_sampName") -> indx
  left_join(dist.df, indx, by=c("sampID1" = "sampID")) %>%
    rename("code1"="code",
           "size1"="size") -> dist.df1 # use the index to add info for 1st sample
  left_join(dist.df1, indx, by=c("sampID2" = "sampID")) %>%
    rename("code2"="code",
           "size2"="size") -> dist.df2 # use the index to add info for the 2nd sample
  
  # remove distances between size classes
  comm.dist<-subset(dist.df2, size1 == size2)
  
  #keep just one of the size columns, since they are now the same
  select(comm.dist, -size2) %>%
    rename("size"="size1") -> comm.dist
  
  return(comm.dist)
  
}

Calc_decayParamDiffs<-function(valueCol, spdf, stemSamples){
  
  # do difference calculations
  spdf<-rename_(spdf, "select.decayParam"=valueCol)
  v <- spdf$select.decayParam
  z <- outer(v,v,'-') # sp1 - sp2 = dist
  colnames(z)<-spdf$code
  row.names(z)<-spdf$code
  z<-abs(z)
  
  #make it long formate
  dist.df<-extract_uniquePairDists(z) #make the distance matrix long and add metadata
  dist.df<-rename(dist.df, "code1"="sp1", "code2"="sp2")
  
  # add back species + size identifiers
  stemSamples %>%
    select(code, size) -> indx
  left_join(dist.df, indx, by=c("code1" = "code")) %>%
    rename("size1"="size") %>%
    left_join(indx, by=c("code2" = "code")) %>%
    rename("size2"="size") -> dist.df # use the index to add info for the 2nd sample

  # remove distances between size classes, keep just one of the size columns, since they are now the same
  dist.df %>%
    filter(size1 == size2) %>%
    select(-size2) %>%
    rename("size"="size1") -> decayparam.dist
  
  return(decayparam.dist)
  
}

Calc_woodTraitDist <- function(traits.mean){
  
  # identify rows with no missing values
  x <- complete.cases(traits.mean[,-(1:2)]) 
  traits.mean1<-traits.mean[x,-(1:3)]
  
  # log-transform and scale, do PCA and select the first 3 axis
  traits.scaled <- apply(log(traits.mean1+10), 2, scale)
  pc <- princomp(traits.scaled)  # stats package
  # pc=rda(traits.scaled) # vegan package
  pc.scores <- pc$scores[, 1:3]  # the first 3 axes
  row.names(pc.scores)<-traits.mean[x,]$code # make a unique identifier for each row
  
  # calculate the euclidean distance in trait space
  pc.dist.mat <- dist(pc.scores, method = "euclidean", diag=TRUE, upper=TRUE)
  mat.traitDist<-as.matrix(pc.dist.mat)
  
  # make it long format
  traitDist.l <- extract_uniquePairDists(dist.mat=mat.traitDist) #make it long
  colnames(traitDist.l)<-c("code1","code2","woodTraitDist")
  
  #add size info
  traitDist.l$size1<-"large"
  traitDist.l[tolower(traitDist.l$code1) == traitDist.l$code1,"size1"]<-"small"
  traitDist.l$size2<-"large"
  traitDist.l[tolower(traitDist.l$code2) == traitDist.l$code2,"size2"]<-"small"
  
  #get rid of distances between size classes
  traitDist.l<-filter(traitDist.l, size1 == size2)
  
  #keep just one of the size columns, since they are now the same
  select(traitDist.l, -size2) %>%
    rename("size"="size1") -> traitDist.l
  
  return(traitDist.l)
}

SummarizeCommDist_byCodePair<-function(comm.dist){
  
  comm.dist %>%
    mutate(codePair=paste(code1, code2, sep="_")) %>%
    select(codePair, dist) %>%
    group_by(codePair) %>%
    summarize(mean=mean(dist),
              se=sd(dist)/sqrt(length(dist)),
              lower=mean-se,
              upper=mean+se) -> summ.comm_dist
  
  #add back size column
  summ.comm_dist$size<-"large"
  summ.comm_dist[tolower(summ.comm_dist$codePair) == summ.comm_dist$codePair,"size"]<-"small"
  
  return(summ.comm_dist)
}

PrepDecayDistForMerge<-function(decayparam.dist){
  
  # elongate decayparm.dist by adding the forward and reverse combinations of codePair
  decayparam.dist %>%
    mutate(codePair=paste(code1, code2, sep="_")) %>%
    mutate(codePair_rev=paste(code2, code1, sep="_")) %>%
    select(codePair, codePair_rev, dist) %>%
    rename("decayparam_dist"="dist") -> decayparam.dist.el
  
  # make forward and reverse dataframes
  df_forward<-select(decayparam.dist.el, codePair, decayparam_dist)
  df_reverse<-select(decayparam.dist.el, codePair_rev, decayparam_dist) %>%
    rename("codePair"="codePair_rev")
  
  # add codePairs where the codes match and have a decayparm distance of 0
  same_codePair<-paste(stemSamples$code, stemSamples$code, sep="_")
  df_same<-data.frame(codePair=same_codePair, decayparam_dist=0, stringsAsFactors = FALSE)
  
  dfs<-list(df_forward=df_forward, df_reverse=df_reverse, df_same=df_same)
  return(dfs)
}

MergeCommNDecaydists_byCodePair<-function(decayparam.dist, summ.comm_dist){
  
  #prep decay dist df for merging by repeating forward and reverse codePairs and adding distance = 0 rows
  dfs<-PrepDecayDistForMerge(decayparam.dist)
  
  #codePairs missing?
  uniqCodePairs<-unique(c(dfs[['df_forward']]$codePair, dfs[['df_reverse']]$codePair, dfs[['df_same']]$codePair))
  missingCodePairs<-summ.comm_dist$codePair[!unique(summ.comm_dist$codePair) %in% uniqCodePairs]
  missingCodePairs
  
  #join community distances with forward and reverse decay distance dataframes
  summ.comm_dist %>%
    left_join(dfs[['df_forward']]) %>% 
    rename("decayparam_dist_forward"="decayparam_dist") %>%
    left_join(dfs[['df_reverse']]) %>%
    rename("decayparam_dist_reverse"="decayparam_dist") %>%
    left_join(dfs[['df_same']]) -> join.dist
  
  #find the values that matched the forward, reverse, and same codePairs and put it all together
  filter(join.dist, is.na(decayparam_dist_forward) & !is.na(decayparam_dist_reverse)) %>%
    mutate(decayparam_dist=decayparam_dist_reverse) -> reverse.rows
  filter(join.dist, !is.na(decayparam_dist_forward) & is.na(decayparam_dist_reverse)) %>%
    mutate(decayparam_dist=decayparam_dist_forward) -> forward.rows
  filter(join.dist, !is.na(decayparam_dist)) -> same.rows
  join.dist<-bind_rows(reverse.rows, forward.rows, same.rows)
  
  #select columns and rename
  join.dist %>%
    select(codePair, size, mean, lower, upper, decayparam_dist) %>%
    rename("mean_comm_dist"="mean",
           "lower_comm_dist"="lower",
           "upper_comm_dist"="upper") -> join.dist
  
  return(join.dist)
  
}

MergeWoodNDecaydists_byCodePair<-function(decayparam.dist, traits.dist){
  
  #make a codePairs column for traits.dist
  traits.dist %>%
    mutate(codePair=paste(code1, code2, sep="_")) -> traits.dist
  
  #prep decay dist df for merging by repeating forward and reverse codePairs and adding distance = 0 rows
  dfs<-PrepDecayDistForMerge(decayparam.dist)
  
  #codePairs missing?
  uniqCodePairs<-unique(c(dfs[['df_forward']]$codePair, dfs[['df_reverse']]$codePair, dfs[['df_same']]$codePair))
  missingCodePairs<-traits.dist$codePair[!unique(traits.dist$codePair) %in% uniqCodePairs]
  missingCodePairs
  
  #join community distances with forward and reverse decay distance dataframes
  traits.dist %>%
    left_join(dfs[['df_forward']]) %>% 
    rename("decayparam_dist_forward"="decayparam_dist") %>%
    left_join(dfs[['df_reverse']]) %>%
    rename("decayparam_dist_reverse"="decayparam_dist") %>%
    left_join(dfs[['df_same']]) -> join.dist
  
  #find the values that matched the forward and reverse codePairs and put it all together
  filter(join.dist, is.na(decayparam_dist_forward) & !is.na(decayparam_dist_reverse)) %>%
    mutate(decayparam_dist=decayparam_dist_reverse) -> reverse.rows
  filter(join.dist, !is.na(decayparam_dist_forward) & is.na(decayparam_dist_reverse)) %>%
    mutate(decayparam_dist=decayparam_dist_forward) -> forward.rows
  filter(join.dist, !is.na(decayparam_dist)) -> same.rows
  join.dist<-bind_rows(reverse.rows, forward.rows, same.rows)
  
  #select key columns
  join.dist<-select(join.dist, codePair, size, woodTraitDist, decayparam_dist)
  
  #select columns
  join.dist %>%
    select(codePair, size, woodTraitDist, decayparam_dist) -> join.dist

  return(join.dist)
  
}

MergeCommNWoodTraitdists_byCodePair<-function(traits.dist, summ.comm_dist){
  
  #prep trait dist df for merging by repeating forward and reverse codePairs and adding distance = 0 rows
  # elongate decayparm.dist by adding the forward and reverse combinations of codePair
  traits.dist %>%
    mutate(codePair=paste(code1, code2, sep="_")) %>%
    mutate(codePair_rev=paste(code2, code1, sep="_")) %>%
    select(codePair, codePair_rev, woodTraitDist) -> traits.dist.el
  
  # make forward and reverse dataframes
  df_forward<-select(traits.dist.el, codePair, woodTraitDist)
  df_reverse<-select(traits.dist.el, codePair_rev, woodTraitDist) %>%
    rename("codePair"="codePair_rev")
  
  # add codePairs where the codes match and have a decayparm distance of 0
  same_codePair<-paste(stemSamples$code, stemSamples$code, sep="_")
  df_same<-data.frame(codePair=same_codePair, woodTraitDist=0, stringsAsFactors = FALSE)

  dfs<-list(df_forward=df_forward, df_reverse=df_reverse, df_same=df_same)
  
  #join community distances with forward and reverse decay distance dataframes
  summ.comm_dist %>%
    left_join(dfs[['df_forward']]) %>% 
    rename("woodTraitDist_forward"="woodTraitDist") %>%
    left_join(dfs[['df_reverse']]) %>%
    rename("woodTraitDist_reverse"="woodTraitDist") %>%
    left_join(dfs[['df_same']]) -> join.dist
  
  #find the values that matched the forward, reverse, and same codePairs and put it all together
  filter(join.dist, is.na(woodTraitDist_forward) & !is.na(woodTraitDist_reverse)) %>%
    mutate(woodTraitDist=woodTraitDist_reverse) -> reverse.rows
  filter(join.dist, !is.na(woodTraitDist_forward) & is.na(woodTraitDist_reverse)) %>%
    mutate(woodTraitDist=woodTraitDist_forward) -> forward.rows
  filter(join.dist, !is.na(woodTraitDist)) -> same.rows
  join.dist<-bind_rows(reverse.rows, forward.rows, same.rows)
  
  #select columns and rename
  join.dist %>%
    select(codePair, size, mean, lower, upper, woodTraitDist) %>%
    rename("mean_comm_dist"="mean",
           "lower_comm_dist"="lower",
           "upper_comm_dist"="upper") -> join.dist
  
  return(join.dist)
  
}


Plot_DistvDist<-function(comm_decay.distList, distType, xCol){
  
  #figure out the x data and label
  if(xCol=="mean_comm_dist"){
    xlab<-paste("Microb comm distance (", distType,")", sep="")
    comm_decay.distList.newCol<-lapply(comm_decay.distList, rename, "Xdist"="mean_comm_dist")
  }
  if(xCol=="woodTraitDist"){
    xlab<-"Wood trait distance"
    comm_decay.distList.newCol<-lapply(comm_decay.distList, rename, "Xdist"="woodTraitDist")
  }
  
  #initalize the ggplots
  pList<-lapply(comm_decay.distList.newCol, ggplot, aes(x=Xdist, y=decayparam_dist, color=size))
  ylabs<-paste("delta",names(pList))
  
  #add plot characteristics
  if(xCol=="mean_comm_dist"){
    for(i in 1:length(pList)){
      pList[[i]]<-pList[[i]] +
        geom_point() + 
        geom_errorbarh(aes(xmin=lower_comm_dist, xmax=upper_comm_dist)) +
        xlab(xlab) + ylab(ylabs[i])
    }
  }
  if(xCol=="woodTraitDist"){
    for(i in 1:length(pList)){
      pList[[i]]<-pList[[i]] +
        geom_point() + 
        xlab(xlab) + ylab(ylabs[i])
    }
  }
  return(pList)
}

Make_commDistvDist_Fig<-function(distType, valueCol_vec, seqSamples, stemSamples, comm.otu, spdf){
  
  #calculate pairwise community distances
  comm.dist<-Calc_commDists(seqSamples, comm.otu, distType=distType) #2. calc the distances
  summ.comm_dist<-SummarizeCommDist_byCodePair(comm.dist) #3. summarize the distances by codePair
  
  #calculate pairwise decay parameter distances
  decayDistList<-list()
  for(i in 1:length(valueCol_vec)){
    valueCol<-valueCol_vec[i]
    decayDistList[[i]]<-Calc_decayParamDiffs(valueCol, spdf, stemSamples)
  }
  names(decayDistList)<-valueCol_vec
  
  #merge community and decay param distances into the same dataframe
  comm_decay.distList<-lapply(X=decayDistList, FUN=MergeCommNDecaydists_byCodePair, summ.comm_dist=summ.comm_dist)
  names(comm_decay.distList)<-names(decayDistList)
  
  #make plots
  pList<-Plot_DistvDist(comm_decay.distList, distType, xCol="mean_comm_dist")
  
  return(pList)
}

Make_woodTraitDistvDist_Fig<-function(valueCol_vec, seqSamples, stemSamples, traits.mean, spdf){
  
  #calculate pairwise wood trait distances
  traits.dist<-Calc_woodTraitDist(traits.mean)
  
  #calculate pairwise decay parameter distances
  decayDistList<-list()
  for(i in 1:length(valueCol_vec)){
    valueCol<-valueCol_vec[i]
    decayDistList[[i]]<-Calc_decayParamDiffs(valueCol, spdf, stemSamples)
  }
  names(decayDistList)<-valueCol_vec
  
  #merge wood trait and decay param distances into the same dataframe
  trait_decay.distList<-lapply(X=decayDistList, FUN=MergeWoodNDecaydists_byCodePair, traits.dist=traits.dist)
  names(trait_decay.distList)<-names(decayDistList)
  
  #make plots
  pList<-Plot_DistvDist(trait_decay.distList, distType="euclid", xCol="woodTraitDist")
  
  return(pList)
}

