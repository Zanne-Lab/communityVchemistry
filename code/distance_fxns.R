
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
  
  #keep just one of the size columns, since they are now the same
  select(comm.dist, -size2) %>%
    rename("size"="size1") -> comm.dist
  
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
  z<-abs(z)
  dist.df<-extract_uniquePairDists(z) #make the distance matrix long and add metadata
  dist.df<-rename(dist.df, "code1"="sp1", "code2"="sp2")
  
  # add back species + size identifiers
  indx<-unique(sampTab[,c("code","size")])
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

SummarizeCommDist_byCodePair<-function(comm.dist){
  
  comm.dist %>%
    mutate(codePair=paste(code1, code2, sep="_")) %>%
    select(codePair, dist) %>%
    group_by(codePair) %>%
    summarize(mean=mean(dist),
              se=sd(dist)/sqrt(length(dist)),
              lower=mean-se,
              upper=mean+se) %>%
    rename("mean_comm_dist"="mean",
           "lower_comm_dist"="lower",
           "upper_comm_dist"="upper") -> summ.comm_dist
  
  #add back size column
  summ.comm_dist$size<-"large"
  summ.comm_dist[tolower(summ.comm_dist$codePair) == summ.comm_dist$codePair,"size"]<-"small"
  
  return(summ.comm_dist)
}

MergeCommNDecaydists_byCodePair<-function(decayparam.dist, summ.comm_dist){
  
  # elongate aic.dist by adding the forward and reverse combinations of codePair
  decayparam.dist %>%
    mutate(codePair=paste(code1, code2, sep="_")) %>%
    mutate(codePair_rev=paste(code2, code1, sep="_")) %>%
    select(codePair, codePair_rev, dist) %>%
    rename("decayparam_dist"="dist") -> decayparam.dist.el
  
  #make forward and reverse dataframes
  df_forward<-select(decayparam.dist.el, codePair,decayparam_dist)
  df_reverse<-select(decayparam.dist.el, codePair_rev,decayparam_dist)
  
  #code pairs missing?
  uniqCodePairs<-unique(c(decayparam.dist.el$codePair, decayparam.dist.el$codePair_rev))
  missingCodePairs<-summ.comm_dist$codePair[!unique(summ.comm_dist$codePair) %in% uniqCodePairs]
  missingCodePairs
  
  #join community distances with forward and reverse dataframes
  summ.comm_dist %>%
    left_join(df_forward) %>% 
    rename("decayparam_dist_forward"="decayparam_dist") %>%
    left_join(df_reverse, by=c("codePair"="codePair_rev")) %>%
    rename("decayparam_dist_reverse"="decayparam_dist") -> join.dist
  
  #find the values that matched the forward and reverse codePairs and put it all together
  filter(join.dist, is.na(decayparam_dist_forward) & !is.na(decayparam_dist_reverse)) %>%
    mutate(decayparam_dist=decayparam_dist_reverse) -> reverse.rows
  filter(join.dist, !is.na(decayparam_dist_forward) & is.na(decayparam_dist_reverse)) %>%
    mutate(decayparam_dist=decayparam_dist_forward) -> forward.rows
  join.dist<-bind_rows(reverse.rows, forward.rows)
  
  #select key columns
  join.dist<-select(join.dist, codePair, size, 
                    mean_comm_dist, lower_comm_dist, upper_comm_dist,
                    decayparam_dist)
  
  return(join.dist)
  
}

Calc_woodTraitDist <- function(traits.mean){
  
  # calculate wood functional trait distance in multivariate space 
  
  # identify rows with no missing values
  x <- complete.cases(traits.mean[,-(1:2)]) 
  traits.mean1<-traits.mean[x,-(1:3)]
  
  # did you lose any species doing that?
  length(unique(traits.mean$species_lower)); length(unique(traits.mean[x,]$species_lower)) #yes, lost 1 species
  #unique(meta$species)[!unique(meta$species) %in% unique(traits.mean[x,]$species)] #this species is missing
  # Olax stricta is missing because it doesn't have a waterperc value and it is only represented in small stem samples 
  
  #log-transform and scale, do PCA and take the first 3 axis, measures euc
  traits.scaled <- apply(log(traits.mean1+10), 2, scale)
  pc <- princomp(traits.scaled)  # stats package
  # pc=rda(traits.scaled) # vegan package
  pc.scores <- pc$scores[, 1:3]  # the first 3 axes
  
  # make a unique identifier for each row
  row.names(pc.scores)<-traits.mean[x,]$code
  pc.dist.mat <- dist(pc.scores, method = "euclidean", diag=TRUE, upper=TRUE)
  mat.traitDist<-as.matrix(pc.dist.mat)
  
  #make it long
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

MergeWoodNDecaydists_byCodePair<-function(decayparam.dist, traits.dist){
  
  #make a codePairs column for traits.dist
  traits.dist %>%
    mutate(codePair=paste(code1, code2, sep="_")) -> traits.dist
  
  # elongate aic.dist by adding the forward and reverse combinations of codePair
  decayparam.dist %>%
    mutate(codePair=paste(code1, code2, sep="_")) %>%
    mutate(codePair_rev=paste(code2, code1, sep="_")) %>%
    select(codePair, codePair_rev, dist) %>%
    rename("decayparam_dist"="dist") -> decayparam.dist.el
  
  #make forward and reverse dataframes
  df_forward<-select(decayparam.dist.el, codePair,decayparam_dist)
  df_reverse<-select(decayparam.dist.el, codePair_rev,decayparam_dist)
  
  #code pairs missing?
  uniqCodePairs<-unique(c(decayparam.dist.el$codePair, decayparam.dist.el$codePair_rev))
  missingCodePairs<-traits.dist$codePair[!unique(traits.dist$codePair) %in% uniqCodePairs]
  missingCodePairs
  
  #join trait distances with forward and reverse dataframes
  traits.dist %>%
    left_join(df_forward) %>% #join with the forward versions of decayparam.df$codePair
    rename("decayparam_dist_forward"="decayparam_dist") %>%
    left_join(df_reverse, by=c("codePair"="codePair_rev")) %>%
    rename("decayparam_dist_reverse"="decayparam_dist") -> join.dist
  
  #find the values that matched the forward and reverse codePairs and put it all together
  filter(join.dist, is.na(decayparam_dist_forward) & !is.na(decayparam_dist_reverse)) %>%
    mutate(decayparam_dist=decayparam_dist_reverse) -> reverse.rows
  filter(join.dist, !is.na(decayparam_dist_forward) & is.na(decayparam_dist_reverse)) %>%
    mutate(decayparam_dist=decayparam_dist_forward) -> forward.rows
  join.dist<-bind_rows(reverse.rows, forward.rows)
  
  #select key columns
  join.dist<-select(join.dist, codePair, size, woodTraitDist, decayparam_dist)
  
  return(join.dist)
  
}


