#--------------------------------------#
# code-level traits

doAnalysis_traits_explain_decayParams <- function(decayfits, traits.code, code.respVars, traitVars, 
                                                  use.cache, save.cache){
  
  #isolate trait data and scale it
  traits.code %>%
    select(c("code","species","size",traitVars)) -> select.traits.code
  select.traits.code <- select.traits.code[complete.cases(select.traits.code),]
  matonly <- select.traits.code[,!colnames(select.traits.code) %in% c("code", "species", "size")]
  matonly.s <- scale(as.matrix(matonly))
  select.traits.code<- data.frame(select.traits.code[,c("code","species","size")], matonly.s)
  
  #isolate decay data
  vars <- c("code","species","size", unlist(code.respVars))
  decayfits %>%
    select(vars) -> select.decayfits
  
  #merge traits and response variables into 1 df
  select.decayfits %>%
    left_join(select.traits.code) %>% 
    filter(!is.na(P)) %>%
    filter(!is.na(waterperc)) -> decayfits.traits
  
  #set up full models
  rhsVars <- paste(c(traitVars))
  rhs <- paste(rhsVars, collapse = " + ")
  mod.full.list<-lapply(code.respVars, ModelFit_manyYs, rhs, curr.data=decayfits.traits) # summaryTable_fxns.R
  
  #do stepwise model selection
  if(use.cache == F){
    mod.select.list<-lapply(mod.full.list, function(x) {
      x.updated<-update(x, . ~ ., data = model.frame(x))
      mod.select<-step(x.updated, direction="backward")
      return(mod.select)
      })
    names(mod.select.list) <- code.respVars
    if(save.cache == T){
      saveRDS(mod.select.list, file = "derived_data/modSelect.RData")
    }
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

makefig__traits_explain_decayParams <- function(decayfits, traits.code, code.respVars, traitVars, 
                                                use.cache, save.cache, cfract){
  
  #do analysis
  if(cfract == F){
    result.traitsDecay <- doAnalysis_traits_explain_decayParams(decayfits, traits.code, code.respVars, traitVars, 
                                                                use.cache, save.cache)
  }else{
    result.traitsDecay <- doAnalysis_cfract_explain_decayParams(decayfits, traits.code, code.respVars, traitVars, 
                                                                use.cache, save.cache)
  }
  
  mod.select.list <- result.traitsDecay$mod.select.list
  
  #set up plotting dataframe
  plot.obj <- MakeLm_plottingDF(mod.list = mod.select.list, respvars = result.traitsDecay$respVars)
  #plot.obj$plotting.df
  #plot.obj$r2.df
  
  #define wood trait levels
  trait.levels <- c("sizesmall", traitVars)
  plot.obj$plotting.df$term <- factor(plot.obj$plotting.df$term, levels = rev(trait.levels))
  
  #define response variable levels
  respvar.levels <- c("w.t50","alpha","beta","w.r2","t50","k","ne.r2")
  plot.obj$plotting.df %>%
    filter(respvar %in% respvar.levels) %>%
    mutate(respvar = factor(respvar, levels = respvar.levels)) -> plot.obj$plotting.df
  
  #plot
  p <- ggplot(plot.obj$plotting.df, aes(x = respvar, y = term, fill = etasq)) +
    geom_tile(color = "black") + 
    geom_text(aes(label = round(est, digits = 2))) +
    xlab("Response variable") + ylab("Wood trait") +
    theme_classic() +
    scale_fill_distiller(direction = 1)
  p
  
  result <- list(p = p, p.r2 = plot.obj$r2.df)
  return(result)
  
}

#--------------------------------------#
# stem-level traits

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

doAnalysis_traits_explain_pmr <- function(datasets, stem.respVars, traitVars.stem, 
                                          use.cache, save.cache){
  
  #set up full models
  rhsVars <- paste(c("size",traitVars.stem))
  rhs <- paste(rhsVars, collapse = " + ")
  lhs <- "curr.pmr"
  mod.full.list<-lapply(datasets, function(x){
    result <- ModelFit_manyYs(y = lhs, rhs=rhs, curr.data=x)
    return(result)
  })
  
  #do stepwise model selection
  if(use.cache == F){
    mod.select.list<-lapply(mod.full.list, function(x) {
      x.updated<-update(x, . ~ ., data = model.frame(x))
      mod.select<-step(x.updated, direction="backward")
      return(mod.select)
    })
    if(save.cache == T){
      saveRDS(mod.select.list, file = "derived_data/modSelect_stem.RData")
    }
  }else{
    mod.select.list <- readRDS(file = "derived_data/modSelect_stem.RData")
  }
  
  #make table
  #prettyTab <- MakeLmSummaryTable(respvars=unlist(stem.respVars), mod.list=mod.select.list) # summaryTable_fxns.R
  #write.csv(prettyTab, "output/traitsummary_stem.csv")
  
  #extract and save the residuals
  traitResiduals.stem <- ExtractResids(mod.list=mod.select.list, dataset.list=datasets, sampleName = "codeStem") # summaryTable_fxns.R
  
  result.traitsPMR <- list(
    datasets = datasets,
    respVars = stem.respVars,
    mod.select.list = mod.select.list,
    traitResiduals.stem = traitResiduals.stem)
  
  return(result.traitsPMR)
  
}

makefig__traits_explain_pmr <- function(datasets, stem.respVars, traitVars.stem, 
                                        use.cache, save.cache){
  
  #do analysis
  result.traitsPMR <- doAnalysis_traits_explain_pmr(datasets, stem.respVars, traitVars.stem, 
                                                    use.cache, save.cache)
  mod.select.list <- result.traitsPMR$mod.select.list
  
  #set up plotting dataframe
  plot.obj <- MakeLm_plottingDF(mod.list = mod.select.list, respvars = result.traitsPMR$respVars)
  #plot.obj$plotting.df
  #plot.obj$r2.df
  
  #define wood trait levels
  trait.levels <- c("sizesmall", traitVars.stem)
  plot.obj$plotting.df$term <- factor(plot.obj$plotting.df$term, levels = rev(trait.levels))
  
  #define response variable levels
  respvar.levels <- unlist(result.traitsPMR$respVars)
  plot.obj$plotting.df %>%
    filter(respvar %in% respvar.levels) %>%
    mutate(respvar = factor(respvar, levels = respvar.levels)) -> plot.obj$plotting.df
  
  #plot
  p <- ggplot(plot.obj$plotting.df, aes(x = respvar, y = term, fill = etasq)) +
    geom_tile(color = "black") + 
    geom_text(aes(label = round(est, digits = 2))) +
    xlab("Response variable") + ylab("Wood trait") +
    theme_classic() +
    scale_fill_distiller(direction = 1)
  p
  #ggsave(filename = "output/figures/maintext/traits_explain_pmr.pdf", width = 5, height = 6)
  
  
  result <- list(p = p, p.r2 = plot.obj$r2.df)
  return(result)
}

#--------------------------------------#
# code-level Cfractions

doAnalysis_cfract_explain_decayParams <- function(decayfits, traits.code, code.respVars, traitVars, 
                                                  use.cache, save.cache){
  
  #isolate trait data and scale it
  traits.code %>%
    select(c("code","species","size",traitVars)) -> select.traits.code
  select.traits.code <- select.traits.code[complete.cases(select.traits.code),]
  matonly <- select.traits.code[,!colnames(select.traits.code) %in% c("code", "species", "size")]
  matonly.s <- scale(as.matrix(matonly))
  select.traits.code<- data.frame(select.traits.code[,c("code","species","size")], matonly.s)
  
  #isolate decay data
  vars <- c("code","species","size", unlist(code.respVars))
  decayfits %>%
    select(vars) -> select.decayfits
  
  #merge traits and response variables into 1 df
  select.decayfits %>%
    left_join(select.traits.code) -> decayfits.traits
  decayfits.traits <- decayfits.traits[complete.cases(decayfits.traits),]
  
  #set up full models
  rhsVars <- paste(c(traitVars))
  rhs <- paste(rhsVars, collapse = " + ")
  mod.full.list<-lapply(code.respVars, ModelFit_manyYs, rhs, curr.data=decayfits.traits) # summaryTable_fxns.R
  
  #do stepwise model selection
  if(use.cache == F){
    mod.select.list<-lapply(mod.full.list, function(x) {
      x.updated<-update(x, . ~ ., data = model.frame(x))
      mod.select<-step(x.updated, direction="backward")
      return(mod.select)
    })
    names(mod.select.list) <- code.respVars
    if(save.cache == T){
      saveRDS(mod.select.list, file = "derived_data/modSelect_cfract.RData")
    }
  }else{
    mod.select.list <- readRDS(file = "derived_data/modSelect_cfract.RData")
  }
  
  #extract and save the residuals
  dataset.list <- rep(list(decayfits.traits), length(code.respVars))
  names(mod.select.list) <- code.respVars
  traitResiduals.code <- ExtractResids(mod.list=mod.select.list, dataset.list=dataset.list, sampleName = "code") # summaryTable_fxns.R
  
  result.traitsDecay <- list(
    decayfits.traits = decayfits.traits,
    respVars = code.respVars,
    mod.select.list = mod.select.list,
    traitResiduals.code = traitResiduals.code)
  
  return(result.traitsDecay)
  
}

#--------------------------------------#
# stem-level Cfractions

CreateCfractPMRpair<-function(respVar, traits.stem.cfract, pmr_byStem, traitVars.cfract){
  
  #make a dataframe using the current time point's pmr and remove NAs
  pmr_byStem %>%
    select("codeStem", respVar) %>%
    rename("curr.pmr" = respVar) %>%
    filter(!is.na(curr.pmr)) -> pmr.noNAs
  #subset the trait matrix using these unique codeStems
  traits.stem.cfract %>%
    filter(codeStem %in% pmr.noNAs$codeStem) -> curr.traits
  #make sure there are no NAs 
  curr.traits <- curr.traits[complete.cases(curr.traits),]
  
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
    select(c("codeStem","code","species","curr.pmr","size", traitVars.cfract)) -> select.traits
  select.traits <- select.traits[complete.cases(select.traits),]
  matonly <- select.traits[,!colnames(select.traits) %in% c("codeStem","code", "species", "curr.pmr","size")]
  matonly.s <- scale(as.matrix(matonly))
  result <- data.frame(select.traits[,c("codeStem","code", "species", "curr.pmr","size")], matonly.s)
  
  return(result)
  
}

doAnalysis_cfract_explain_pmr <- function(datasets, stem.respVars, traitVars.cfract, 
                                          use.cache, save.cache){
    
  #set up full models
  rhsVars <- paste(c(traitVars.cfract))
  rhs <- paste(rhsVars, collapse = " + ")
  lhs <- "curr.pmr"
  mod.full.list<-lapply(datasets, function(x){
    result <- ModelFit_manyYs(y = lhs, rhs=rhs, curr.data=x)
    return(result)
  })
  
  #do stepwise model selection
  if(use.cache == F){
    mod.select.list<-lapply(mod.full.list, function(x) {
      x.updated<-update(x, . ~ ., data = model.frame(x))
      mod.select<-step(x.updated, direction="backward")
      return(mod.select)
    })
    if(save.cache == T){
      saveRDS(mod.select.list, file = "derived_data/modSelect_stem_cfract.RData")
    }
  }else{
    mod.select.list <- readRDS(file = "derived_data/modSelect_stem_cfract.RData")
  }
  
  #make table
  #prettyTab <- MakeLmSummaryTable(respvars=unlist(stem.respVars), mod.list=mod.select.list) # summaryTable_fxns.R
  #write.csv(prettyTab, "output/traitsummary_stem.csv")
  
  #extract and save the residuals
  traitResiduals.stem <- ExtractResids(mod.list=mod.select.list, dataset.list=datasets, sampleName = "codeStem") # make_summaryTables.R
  
  result.traitsPMR <- list(
    datasets = datasets,
    respVars = stem.respVars,
    mod.select.list = mod.select.list,
    traitResiduals.stem = traitResiduals.stem)
  
  return(result.traitsPMR)
  
}

makefig__cfract_explain_pmr <- function(datasets, stem.respVars, traitVars.cfract, 
                                        use.cache, save.cache){
  
  #do analysis
  result.traitsPMR <- doAnalysis_cfract_explain_pmr(datasets, stem.respVars, traitVars.cfract, 
                                                    use.cache, save.cache)
  mod.select.list <- result.traitsPMR$mod.select.list
  
  #set up plotting dataframe
  plot.obj <- MakeLm_plottingDF(mod.list = mod.select.list, respvars = result.traitsPMR$respVars)
  #plot.obj$plotting.df
  #plot.obj$r2.df
  
  #define wood trait levels
  trait.levels <- c(traitVars.cfract)
  plot.obj$plotting.df$term <- factor(plot.obj$plotting.df$term, levels = rev(trait.levels))
  
  #define response variable levels
  respvar.levels <- unlist(result.traitsPMR$respVars)
  plot.obj$plotting.df %>%
    filter(respvar %in% respvar.levels) %>%
    mutate(respvar = factor(respvar, levels = respvar.levels)) -> plot.obj$plotting.df
  
  #plot
  p <- ggplot(plot.obj$plotting.df, aes(x = respvar, y = term, fill = etasq)) +
    geom_tile(color = "black") + 
    geom_text(aes(label = round(est, digits = 2))) +
    xlab("Response variable") + ylab("Wood trait") +
    theme_classic() +
    scale_fill_distiller(direction = 1)
  p
  #ggsave(filename = "output/figures/maintext/traits_explain_pmr.pdf", width = 5, height = 6)
  
  
  result <- list(p = p, p.r2 = plot.obj$r2.df)
  return(result)
}


