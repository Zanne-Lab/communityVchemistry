
ModelFit_manyYs<-function(y, rhs, curr.data){
  
  #create model formula
  string<-paste(y, " ~ ", rhs)
  fmla<-as.formula(string)
  
  #fit full model
  mod.full<-lm(formula=fmla, data=curr.data)
  
  #return a list with the best model for each response var
  return(mod.full)
}

PullLmCoefs<-function(sum.list, respvars){
  
  coefs.list<-lapply(sum.list, function(x){
    df<-data.frame(term=row.names(x$coefficients), x$coefficients)
    colnames(df)<-c("term","est","se","t.value","pval")
    return(df)
  })
  names(coefs.list)<-respvars
  coefs.df<-bind_rows(coefs.list, .id="respvar")
  #add pval stars
  coefs.df$stars<-""
  coefs.df[coefs.df$pval < 0.05, "stars"]<-"*"
  coefs.df[coefs.df$pval < 0.01, "stars"]<-"**"
  coefs.df[coefs.df$pval < 0.001, "stars"]<-"***"
  #create print vector
  coefs.df$printvec<-paste(round(coefs.df$est, digits=3), "+/-", 
                           round(coefs.df$se, digits=3),
                           coefs.df$stars)
  
  return(coefs.df)
}

PullLmFitStats<-function(sum.list, respvars){
  
  fitstats.list<-lapply(sum.list, function(x){
    df<-data.frame(Fstat = round(x$fstatistic['value'], digits=2),
                   numdf = x$fstatistic['numdf'],
                   dendf = x$fstatistic['dendf'],
                   r.squared = round(x$r.squared, digits=2))
    return(df)
  })
  names(fitstats.list)<-respvars
  fitstats.df<-bind_rows(fitstats.list)
  row.names(fitstats.df)<-respvars
  fitstat.df<-data.frame(t(fitstats.df))
  fitstat.df$term<-row.names(fitstat.df)
  
  return(fitstat.df)
}

MakeLmSummaryTable<-function(respvars, mod.list, termorder){
  
  #summarize
  sum.list<-lapply(mod.list, summary)
  
  #pull term coefs
  coefs.df<-PullLmCoefs(sum.list, respvars)
  
  #pull model fit stats
  fitstat.df<-PullLmFitStats(sum.list, respvars)
  
  #reformat
  coefs.df %>%
    select(respvar, term, printvec) %>%
    spread(respvar, printvec) -> coefs.df
  if(termorder=='stemTraits'){
    row.order<-c("(Intercept)","sizesmall","waterperc","density_smspp","barkthick_smspp","C","N","P","Mn")
  }
  if(termorder=='codeTraits'){
    row.order<-c("(Intercept)","sizesmall","waterperc","barkthick","C","N","Ca","Zn")
  }
  if(termorder=='none'){
    row.order<-coefs.df$term
  }
  coefs.df %>%
    slice(match(row.order, term)) -> coefs.df
  colorder<-c("term", respvars)
  coefs.df<-coefs.df[,colorder]
  prettyTab<-rbind(coefs.df, rep(NA, length(colnames(coefs.df))), fitstat.df)
  
  return(prettyTab)
  
}

ExtractResids<-function(mod.list, dataset.list, sampleName){
  
  resid.dfs<-list()
  for(i in 1:length(mod.list)){
    mod <- mod.list[[i]]
    dataset <- dataset.list[[i]]
    resid.dfs[[i]]<-data.frame(resid=mod$resid, sampleName=dataset[,sampleName])
  }
  names(resid.dfs)<-names(mod.list)
  resid.df<-bind_rows(resid.dfs, .id="resp")
  
  return(resid.df)
}

Do_coxTests<-function(mod.stem.list, mod.code.list, respvars){
  cox.sum.list<-list()
  for(i in 1:length(mod.stem.list)){
    curr.mod.stem<-mod.stem.list[[i]]
    curr.mod.code<-mod.code.list[[i]]
    cox.sum.list[[i]]<-coxtest(curr.mod.stem, curr.mod.code)
  }
  names(cox.sum.list)<-respvars
  return(cox.sum.list)
}

Do_compareR2<-function(mod.stem.list, mod.code.list, respvars){
  stem.r2s<-unlist(lapply(mod.stem.list, function(x) round(summary(x)$r.squared, digits=2)))
  code.r2s<-unlist(lapply(mod.code.list, function(x) round(summary(x)$r.squared, digits=2)))
  summ.r2s<-data.frame(respvars, code.r2s, stem.r2s)
  return(summ.r2s)
}

MakeSummaryTable_comcomp<-function(wapls.out, respvars){
  
  tmp<-lapply(wapls.out, function(x) x["Comp01",])
  tmp<-do.call(cbind,lapply(tmp,data.frame))
  colnames(tmp)<-respvars
  prettyTab<-data.frame(stat=row.names(tmp), round(tmp, digits=2))
  return(prettyTab)
}

MakeSummaryTable_diversity<-function(respvars, mod.list){
  
  #summarize
  sum.list<-lapply(mod.list, summary)
  
  #pull term coefs
  coefs.df<-PullLmCoefs(sum.list, respvars)
  
  #pull model fit stats
  fitstat.df<-PullLmFitStats(sum.list, respvars)
  
  #reformat
  coefs.df %>%
    select(respvar, term, printvec) %>%
    spread(respvar, printvec) -> coefs.df
  prettyTab<-rbind(coefs.df, rep(NA, length(colnames(coefs.df))), fitstat.df)
  
  return(prettyTab)
  
}
