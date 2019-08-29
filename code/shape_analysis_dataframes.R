
#########
# Create response + predictor dataframes
CreateTraitPMRpair<-function(respVar, traits.stem, traits.code, pmr_byStem, traitVars.stem){
  
  #make a dataframe using the current time point's pmr and remove NAs
  pmr_byStem %>%
    select("codeStem", respVar) %>%
    rename("curr.pmr" = respVar) %>%
    filter(!is.na(curr.pmr)) -> pmr.noNAs
  #subset the trait matrix using these unique codeStems
  traits.stem %>%
    filter(codeStem %in% pmr.noNAs$codeStem) -> curr.traits
  #make sure there are no NAs in waterperc or chemistry data
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
  
  #isolate trait data and scale it
  curr.df %>%
    select(c("codeStem","code","species","curr.pmr","size", traitVars.stem)) -> select.traits
  select.traits <- select.traits[complete.cases(select.traits),]
  matonly <- select.traits[,!colnames(select.traits) %in% c("codeStem","code", "species", "curr.pmr","size")]
  matonly.s <- scale(as.matrix(matonly))
  result <- data.frame(select.traits[,c("codeStem","code", "species", "curr.pmr","size")], matonly.s)
  
  return(result)
  
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

CreateRichTraitResid_df<-function(otutype.df, trait.residuals){
  
  #summarize by code
  otutype.df %>%
    group_by(code) %>%
    summarize(mean=mean(sub_rich),
              se=sd(sub_rich)/sqrt(length(sub_rich)),
              upper=mean+se,
              lower=mean-se) -> summ.otutype.df
  
  # join with trait residual data
  rich.traitresid<-left_join(trait.residuals, summ.otutype.df)
  
  #add size col
  colnames(rich.traitresid)
  rich.traitresid$size <- "large"
  rich.traitresid[rich.traitresid$code == tolower(rich.traitresid$code), "size"] <- "small"
  
  return(rich.traitresid)
}

CreateCommTraitpair <- function(comm.otu, traits, sampleName){
  
  #identify the column with the unique row info (sampleName)
  excludeCols<-c("codeStem","code","species")
  traits %>%
    rename_("sampleName"=sampleName) -> traits
  traits<-traits[,!colnames(traits) %in% excludeCols]
  
  # match the no-nas trait data with community data
  ord<-match(row.names(comm.otu), traits$sampleName)
  traits.noNAs.o<-traits[ord,]
  comm.otu.trim <- comm.otu[!is.na(traits.noNAs.o$sampleName),]
  traits.noNAs.o.trim<-traits.noNAs.o[!is.na(traits.noNAs.o$sampleName),]
  
  # make the trait df a matrix with only predictor values
  traits.noNAs.o.trim %>%
    select(-sampleName) %>%
    data.frame() -> traits.mat
  result<-list(comm =  comm.otu.trim,
               traits = traits.mat)
  
  test <- sum(row.names(comm.otu.trim) != traits.noNAs.o.trim$sampleName) # this needs to be 0
  if(test == 0){
    return(result)
  }else{
    return("non-matching")
  }
}


#########
# Model formulation

ModelFit_manyYs<-function(y, rhs, curr.data){
  
  #create model formula
  string<-paste(y, " ~ ", rhs)
  fmla<-as.formula(string)
  
  #fit full model
  mod.full<-lm(formula=fmla, data=curr.data)
  
  #return a list with the best model for each response var
  return(mod.full)
}

fitNcrossval_WAPLS<-function(curr.comm, curr.respVar){
  fit<-WAPLS(y=curr.comm, x=curr.respVar)
  fit.cv<-crossval(fit, cv.method="loo")
  return(fit.cv)
}

ordistep_wrapper<-function(datasets){
  
  cap.env <- capscale(datasets[['comm']] ~ ., data = datasets[['traits']], distance='bray') # full model
  mod0.env <- capscale(datasets[['comm']] ~ 1, data = datasets[['traits']], distance='bray') # set up the null cases with no predictors
  step.env <- ordistep(mod0.env, scope=formula(cap.env)) # model selection -- this take a while
  
  return(step.env)
  
}






