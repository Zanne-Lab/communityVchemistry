


#########
# LM summary stats

PullLmCoefs<-function(sum.list, respvars){

  # extract coefs from summary table
  coefs.list<-lapply(sum.list, function(x){
    df<-data.frame(term=row.names(x$coefficients), x$coefficients)
    colnames(df)<-c("term","est","se","t.value","pval")
    return(df)
  })
  names(coefs.list) <- respvars
  coefs.df <- list_to_df(coefs.list)
  coefs.df %>%
    dplyr::rename('respvar'='source') -> coefs.df
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
  
  # extract model fit stats
  fitstats.list<-lapply(sum.list, function(x){
    df<-data.frame(Fstat = round(x$fstatistic['value'], digits=2),
                   numdf = x$fstatistic['numdf'],
                   dendf = x$fstatistic['dendf'],
                   r.squared = round(x$r.squared, digits=2))
    return(df)
  })
  fitstats.df <- list_to_df(fitstats.list)
  fitstat.df<-data.frame(t(fitstats.df))
  fitstat.df$term<-row.names(fitstat.df)
  fitstat.df %>%
    filter(term != 'source') -> fitstat.df
  
  return(fitstat.df)
}

PullAnova<-function(mod.list, respvars){
  
  # extract anova table (type II) and eta.sq for each predictor variable
  #require(sjstats) #this should already be loaded
  #Eta^2 (or Eta-sqr) is the proportion of variance variance associated with one of more main effects
  #Eta^2 = SSeffect / SStotal
  anova.list <- lapply(mod.list, function(x){
    df <- anova_stats(car::Anova(x, type = 2))
    df %>%
      select(term, sumsq, df, statistic, p.value, etasq)-> df
    return(df)
  })
  
  # turn into a df and filter for selected respvars
  anova.df <- list_to_df(anova.list)
  anova.df %>%
    rename('respvar'='source') %>%
    filter(respvar %in% unlist(respvars)) -> anova.df
  
  return(anova.df)
}

MakeLm_plottingDF <- function(mod.list, respvars){
  
  # extract model R2
  fitstat.df<-PullLmFitStats(mod.list, unlist(respvars))
  fitstat.df %>%
    gather(key = "respvar", value = "value", -c(term)) %>%
    filter(term == "r.squared") %>%
    select(-term) %>%
    rename('r.squared'=value) -> r2.df
  
  # extract model coefs
  coefs.df<-PullLmCoefs(mod.list, unlist(respvars))
  coefs.df %>%
    left_join(r2.df) %>%
    select(respvar, term, est, se) %>%
    filter(term != "(Intercept)")-> coefs.df
  
  # extract anova table (type II) and eta.sq for each predictor variable
  anova.df<-PullAnova(mod.list, unlist(respvars))
  anova.df %>%
    select(respvar, term, p.value, etasq) %>%
    filter(term != "Residuals") %>%
    mutate(term = ifelse(term == "size", "sizesmall", term)) -> anova.df
  
  # combine model coefs and anova
  coefs.df %>%
    full_join(anova.df) -> plotting.df

  result <- list(plotting.df = plotting.df, 
                 r2.df = r2.df)
  
  return(result)
  
}

MakeLmSummaryTable<-function(mod.list, respvars){
  
  #summarize
  sum.list<-lapply(mod.list, summary)
  
  #pull term coefs
  coefs.df<-PullLmCoefs(sum.list, respvars)
  
  #pull model fit stats
  fitstat.df<-PullLmFitStats(sum.list, respvars)
  
  #reformat
  coefs.df %>%
    mutate(print.est = paste(signif(est, 2), stars)) %>%
    select(respvar, term, print.est) %>%
    spread(respvar, print.est) -> coefs.df
  
  prettyTab<-rbind(coefs.df, rep(NA, length(colnames(coefs.df))), fitstat.df)
  
  return(prettyTab)
  
}


#########
# LM residuals (for use in other analyses)

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


#########
# LM fit comparisons

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


#########
# WAPLS summary stats

MakeSummaryTable_comcomp<-function(wapls.out, respvars){
  
  tmp<-lapply(wapls.out, function(x) x["Comp01",])
  tmp<-do.call(cbind,lapply(tmp,data.frame))
  colnames(tmp)<-respvars
  prettyTab<-data.frame(stat=row.names(tmp), round(tmp, digits=2))
  return(prettyTab)
}


#########
# db-RDA summary stats

extract_constrainedInertia_proport<-function(dbrda.obj){
  constrained.eig<-dbrda.obj$CCA$tot.chi
  unconstrained.eig<-dbrda.obj$CA$tot.chi
  constrained.inert.prop<-constrained.eig/(constrained.eig+unconstrained.eig)
  
  return(constrained.inert.prop)
}

anova.margin.table<-function(dbrda.obj){
  
  # anova table by margin
  an.obj<-anova(dbrda.obj, by="margin")
  an.df<-data.frame(term=row.names(an.obj), 
                    df=an.obj$Df,
                    Fval=round(an.obj$F, digits=2),
                    pval=an.obj$`Pr(>F)`)
  
  # inertia
  constr.inertia<-dbrda.obj$CCA$tot.chi
  unconstr.inertia<-dbrda.obj$CA$tot.chi
  inertia.df<-data.frame(term=c("Constrained inertia","Unconstrained inertia"),
                         df=c(round(constr.inertia, digits=2),
                              round(unconstr.inertia, digits=2)),
                         Fval=c(NA, NA),
                         pval=c(NA, NA))
  
  prettyTab<-rbind(an.df, rep(NA, length(colnames(an.df))), inertia.df)
  prettyTab
  
  return(prettyTab)
}

