
doAnalysis_traits_explain_decayParams <- function(decayfits, traits.code, code.respVars, use.cache){
  
  #merge traits and response variables into 1 df
  vars <- c("code","species","size", unlist(code.respVars))
  decayfits %>%
    select(!!vars) %>%
    left_join(traits.code) %>% 
    filter(!is.na(P)) %>%
    filter(!is.na(waterperc)) -> decayfits.traits
  
  #set up full models
  rhs <- "size + waterperc + density + barkthick + P + K + Ca + Mn + Fe + Zn + N + C"
  mod.full.list<-lapply(code.respVars, ModelFit_manyYs, rhs, curr.data=decayfits.traits) # summaryTable_fxns.R
  
  #do stepwise model selection
  if(use.cache == F){
    mod.select.list<-lapply(mod.full.list, function(x) {
      x.updated<-update(x, . ~ ., data = model.frame(x))
      mod.select<-step(x.updated, direction="backward")
      return(mod.select)
      })
    names(mod.select.list) <- respVars
    saveRDS(mod.select.list, file = "derived_data/modSelect.RData")
  }else{
    mod.select.list <- readRDS(file = "derived_data/modSelect.RData")
  }

  #make table
  #prettyTab <- MakeLmSummaryTable(respvars=unlist(respVars), mod.list=mod.select.list) # summaryTable_fxns.R
  #write.csv(prettyTab, "output/traits_explain_decayparams.csv")
  
  #extract and save the residuals
  dataset.list <- rep(list(decayfits.traits), length(code.respVars))
  names(mod.select.list) <- code.respVars
  traitResiduals.code<-ExtractResids(mod.list=mod.select.list, dataset.list=dataset.list, sampleName = "code") # summaryTable_fxns.R
  
  result.traitsDecay <- list(
       decayfits.traits = decayfits.traits,
       respVars = code.respVars,
       mod.select.list = mod.select.list,
       traitResiduals.code = traitResiduals.code)
  
  return(result.traitsDecay)
  
}

makefig__traits_explain_decayParams <- function(decayfits, traits.code, code.respVars, use.cache){
  
  #do analysis
  result.traitsDecay <- doAnalysis_traits_explain_decayParams(decayfits, traits.code, code.respVars, use.cache)
  mod.select.list <- result.traitsDecay$mod.select.list
  respVars <- result.traitsDecay$respVars
  
  #plot
  sum.list<-lapply(mod.select.list, summary)
  coefs.df<-PullLmCoefs(sum.list, unlist(respVars))
  fitstat.df<-PullLmFitStats(sum.list, unlist(respVars))
  fitstat.df %>%
    gather(key = "respvar", value = "value", -c(term)) %>%
    filter(term == "r.squared") %>%
    select(-term) %>%
    rename('r.squared'=value) -> r2.df
  coefs.df %>%
    left_join(r2.df) %>%
    filter(respvar %in% c("alpha","beta","w.t70","w.r2")) -> coefs.df
  
  ggplot(coefs.df, aes(x = respvar, y = term, fill = r.squared)) +
    geom_tile() + geom_text(aes(label = round(est, digits = 4))) +
    xlab("Response variable") + ylab("Wood trait") +
    theme_classic()
  
  ggsave(filename = "output/figures/maintext/traits_explain_decayparams.pdf", width = 5, height = 6)
  
  
}

makeDF_variation_densityNbarkthick <- function(traits.stem, traits.code, pmr_byStem){
  
  # prep a dataset with pmr, stem-level and code-level density and barkthick
  traits.stem %>%
    filter(size=="small") %>%
    select(codeStem, code, density, barkthick) %>%
    filter(!is.na(density) & !is.na(barkthick)) %>%
    rename("density_stem"="density",
           "barkthick_stem"="barkthick") -> db.stem
  traits.code %>%
    filter(size=="small") %>%
    select(code, density, barkthick) %>%
    filter(!is.na(density) & !is.na(barkthick)) %>%
    rename("density_code"="density",
           "barkthick_code"="barkthick") -> db.code
  db.stem %>%
    left_join(db.code) %>%
    left_join(pmr_byStem) -> db.df
  
  return(db.df)
  
}

makefig__variation_densityNbarkthick <- function(traits.stem, traits.code, pmr_byStem){
  
  db.df <- makeDF_variation_densityNbarkthick(traits.stem, traits.code, pmr_byStem)
    
  # look at the within-species variation in density and barkthick
  p.denVar<-ggplot(db.df, 
                   aes(x=reorder(code, density_code), y=density_stem)) + 
    geom_point() + theme_bw() + coord_flip() +
    xlab("Wood species") + ylab("Density (g/cm3)")
  p.barVar<-ggplot(db.df, 
                   aes(x=reorder(code, barkthick_code), y=barkthick_stem)) + 
    geom_point() + theme_bw() + coord_flip() +
    xlab("Wood species") + ylab("Bark thickness (mm)") +
    theme(axis.text.y = element_blank(), 
          axis.ticks.y = element_blank(), 
          axis.title.y = element_blank(),
          plot.margin = unit(c(0,1,0,0), "lines"),
          plot.background = element_blank())
  
  pdf("output/figures/supplementary/variation_densityNbarkthick.pdf", width=5, height=5)
  grid.newpage()
  grid.draw(cbind(ggplotGrob(p.denVar), ggplotGrob(p.barVar), size = "last"))
  dev.off()
  
}

extr_cox_results <-function(x){
  df <- data.frame(mod = row.names(x), pval = x$`Pr(>|z|)`) #annotation <- c("M1 = stem, M2 = code")
  return(df)
}

doAnalysis_variation_densityNbarkthick <- function(traits.stem, traits.code, pmr_byStem, stem.respVars){
  
  db.df <- makeDF_variation_densityNbarkthick(traits.stem, traits.code, pmr_byStem)
  
  #set up models
  mod.stem_density.list <- lapply(stem.respVars, ModelFit_manyYs, rhs="density_stem", curr.data=db.df) # summaryTable_fxns.R
  mod.code_density.list <- lapply(stem.respVars, ModelFit_manyYs, rhs="density_code", curr.data=db.df) # summaryTable_fxns.R
  mod.stem_barkthick.list <- lapply(stem.respVars, ModelFit_manyYs, rhs="barkthick_stem", curr.data=db.df) # summaryTable_fxns.R
  mod.code_barkthick.list <- lapply(stem.respVars, ModelFit_manyYs, rhs="barkthick_code", curr.data=db.df) # summaryTable_fxns.R
  
  #cox test to compare non-nested models
  # if the first model (M1) contains the correct set of regressors, then a fit of the regressors from the second model (M2) to the fitted values from M1 should have no further explanatory value. But if it has, it can be concluded that M1 does not contain the correct set of regressors
  cox.density <- Do_coxTests(mod.stem_density.list, 
                             mod.code_density.list, 
                             respvars=unlist(stem.respVars))
  # stem level density information improves estimates for time7 and time59
  cox.barkthick <- Do_coxTests(mod.stem_barkthick.list, 
                               mod.code_barkthick.list, 
                               respvars=unlist(stem.respVars))
  
  cox.density.df <- list_to_df(lapply(cox.density, extr_cox_results))
  cox.density.df$meas <- "density"
  cox.barkthick.df <- list_to_df(lapply(cox.barkthick, extr_cox_results))
  cox.barkthick.df$meas <- "barkthick"
  cox.df <- rbind(cox.density.df, cox.barkthick.df)
  
  return(cox.df)
  
}

doAnalysis_traits_explain_pmr <- function(pmr_byStem, traits.code, traits.stem, stem.respVars, use.cache){
  
  # create dataframes
  
  # stem-level waterperc, xrf, and CN and (small) species-level barkthickness and density
  datasets<-lapply(stem.respVars, function(x) {
    result<-CreateTraitPMRpair(x, traits.stem, traits.code, pmr_byStem) #analysisDF_fxns.R
    return(result)
  })
  names(datasets)<-stem.respVars
  
  #set up full models
  rhs <- "size + waterperc + density_smspp + barkthick_smspp + P + K + Ca + Mn + Fe + Zn + N + C"
  lhs <- "curr.pmr"
  mod.full.list<-lapply(datasets, function(x) {
    result<-ModelFit_manyYs(y=lhs, rhs=rhs, curr.data=x)
    return(result)
  })
  
  #do stepwise model selection
  if(use.cache == F){
    mod.select.list<-lapply(mod.full.list, function(x) {
      x.updated<-update(x, . ~ ., data = model.frame(x))
      mod.select<-step(x.updated, direction="backward")
      return(mod.select)
    })
    saveRDS(mod.select.list, file = "derived_data/modSelect_stem.RData")
  }else{
    mod.select.list <- readRDS(file = "derived_data/modSelect_stem.RData")
  }
  
  #make table
  #prettyTab <- MakeLmSummaryTable(respvars=unlist(stem.respVars), mod.list=mod.select.list) # summaryTable_fxns.R
  #write.csv(prettyTab, "output/traitsummary_stem.csv")
  
  #extract and save the residuals
  traitResiduals.stem<-ExtractResids(mod.list=mod.select.list, dataset.list=datasets, sampleName = "codeStem") # summaryTable_fxns.R
  
  result.traitsPMR <- list(
    datasets = datasets,
    respVars = stem.respVars,
    mod.select.list = mod.select.list,
    traitResiduals.stem = traitResiduals.stem)
  
  return(result.traitsPMR)
  
}

makefig__traits_explain_pmr <- function(traits.stem, traits.code, pmr_byStem, stem.respVars, use.cache){
  
  #do analysis
  result.traitsPMR <- doAnalysis_traits_explain_pmr(pmr_byStem, traits.code, traits.stem, stem.respVars, use.cache)
  mod.select.list <- result.traitsPMR$mod.select.list
  respVars <- result.traitsPMR$respVars
  
  #plot
  sum.list<-lapply(mod.select.list, summary)
  coefs.df<-PullLmCoefs(sum.list, unlist(respVars))
  fitstat.df<-PullLmFitStats(sum.list, unlist(respVars))
  fitstat.df %>%
    gather(key = "respvar", value = "value", -c(term)) %>%
    filter(term == "r.squared") %>%
    select(-term) %>%
    rename('r.squared'=value) -> r2.df
  coefs.df %>%
    left_join(r2.df) -> coefs.df
  
  coefs.df$respvar<- factor(coefs.df$respvar, levels = c("time7","time13","time25","time37","time59"))
  
  ggplot(coefs.df, aes(x = respvar, y = term, fill = r.squared)) +
    geom_tile() + geom_text(aes(label = round(est, digits = 4))) +
    xlab("Response variable") + ylab("Wood trait") +
    theme_classic()
  
  ggsave(filename = "output/figures/maintext/traits_explain_pmr.pdf", width = 5, height = 6)
  
}
