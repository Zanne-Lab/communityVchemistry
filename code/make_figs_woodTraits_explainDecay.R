#--------------------------------------#
# code-level traits

doAnalysis_traits_explain_decayParams <- function(decayfits, traits.code, code.respVars, traitVars, cfract){
  
  #isolate trait data and scale it
  traits.code %>%
    ungroup(code) %>%
    select(c("code","species","size",traitVars)) %>%
    na.omit %>%
    mutate_at(vars(traitVars), list(s = scale2)) -> select.traits.code  # see helper_fxns.R
  
  #isolate decay data
  decayfits %>%
    select(c("code","species","size", unlist(code.respVars))) -> select.decayfits
  
  #merge traits and response variables into 1 df
  select.decayfits %>%
    left_join(select.traits.code) %>%
    na.omit -> decayfits.traits
  #remove the na.action attribute so that it does not interfere with prediction downstream
  attr(decayfits.traits, "na.action") <- NULL
  
  #set up full models
  scaled.traitvars <- paste(traitVars,"s",sep = "_")
  if(cfract == FALSE){
    rhsVars <- paste(c("size",scaled.traitvars))
  }else{
    rhsVars <- paste(c(scaled.traitvars))
  }
  rhs <- paste(rhsVars, collapse = " + ")
  print(rhs)
  
  # fit models
  code.respVars %>%
    purrr::set_names() %>%
    map(~fit.lm(y = .x, rhs, data = decayfits.traits)) -> mod.full.list  # see helper_fxns.R
  
  #do stepwise model selection
  mod.full.list %>%
    map(~backward_selection(.x, data = decayfits.traits)) -> mod.select.list  # see helper_fxns.R
  
  #extract and save the residuals
  mod.select.list %>%
    map(~residuals(.x)) %>%
    map(~data.frame(resid=.x, code = decayfits.traits$code)) %>%
    bind_rows(.id = "resp") -> residuals
  
  #create prediction dataframes for each trait and model
  scaled.traitvars %>%
    purrr::set_names() %>%
    map(~preddat_fun_bysize_allmodels(models = mod.select.list,
                                      data = decayfits.traits, # see helper_fxns.R
                                      curr.traitVar = .x,
                                      data.list = F)) -> preddat
  
  #save model fit stats
  mod.select.list %>%
    map(~summary(.x)) %>%
    map(~data.frame(Fstat = round(.x$fstatistic['value'], digits=2),
                    numdf = .x$fstatistic['numdf'],
                    dendf = .x$fstatistic['dendf'],
                    r.squared = round(.x$r.squared, digits=2))) %>%
    list_to_df() %>%
    select(-source) %>%
    t() -> fitstats
  

  result.traitsDecay <- list(
    data = decayfits.traits,
    respVars = code.respVars,
    models = mod.select.list,
    residuals = residuals,
    preddat = preddat,
    fitstats = fitstats)
  
  return(result.traitsDecay)
  
}

makefig__traits_explain_decayParams <- function(result.list, traitVars, cfract){
  
  #set up plotting dataframe
  plot.obj <- MakeLm_plottingDF(mod.list = result.list$models, 
                                respvars = result.list$respVars)
  #plot.obj$plotting.df
  #plot.obj$r2.df
  
  #define wood trait levels
  if(cfract == F){
    trait.levels <- c("sizesmall", paste(traitVars, "s", sep = "_"))
    trait.labels <- c("sizesmall", traitVars)
  }else{
    trait.levels <- c(paste(traitVars, "s", sep = "_"))
    trait.labels <- c(traitVars)
  }
  plot.obj$plotting.df$term <- factor(plot.obj$plotting.df$term, 
                                      levels = rev(trait.levels),
                                      labels = rev(trait.labels))
  
  #define response variable levels
  respvar.levels <- c("w.t50","beta","alpha", "w.r2","t50","k","ne.r2")
  respvar.labels <- c("Years to 50%\nmass loss","Scale param.","Shape param.","Weibull R2",
                      "t50","k","ne.r2")
  plot.obj$plotting.df %>%
    filter(respvar %in% respvar.levels) %>%
    mutate(respvar = factor(respvar, levels = respvar.levels, 
                            labels = respvar.labels)) -> plot.obj$plotting.df
  
  #plot
  p <- ggplot(plot.obj$plotting.df, aes(x = respvar, y = term, fill = etasq)) +
    geom_tile(color = "black") + 
    geom_text(aes(label = round(est, digits = 2)), size = 3) +
    xlab("Response variable") + ylab("Wood trait") +
    theme_classic() +
    scale_fill_distiller(direction = 1) +
    scale_y_discrete(limits = rev(trait.labels)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
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
  
  result.list <- list(p.denVar = p.denVar, p.barVar = p.barVar)
  return(result.list)
  
}

extr_cox_results <-function(x){
  df <- data.frame(mod = row.names(x), pval = x$`Pr(>|z|)`) #annotation <- c("M1 = stem, M2 = code")
  return(df)
}

doAnalysis_variation_densityNbarkthick <- function(traits.stem, traits.code, pmr_byStem, stem.respVars){
  
  db.df <- makeDF_variation_densityNbarkthick(traits.stem, traits.code, pmr_byStem)
  
  # set up models
  stem.respVars %>%
    purrr::set_names() %>%
    map(~fit.lm(y = .x, rhs = "density_stem", data = db.df)) -> mod.stem_density.list  # see helper_fxns.R
  
  stem.respVars %>%
    purrr::set_names() %>%
    map(~fit.lm(y = .x, rhs = "density_code", data = db.df)) -> mod.code_density.list  # see helper_fxns.R
  
  stem.respVars %>%
    purrr::set_names() %>%
    map(~fit.lm(y = .x, rhs = "barkthick_stem", data = db.df)) -> mod.stem_barkthick.list  # see helper_fxns.R
  
  stem.respVars %>%
    purrr::set_names() %>%
    map(~fit.lm(y = .x, rhs = "barkthick_code", data = db.df)) -> mod.code_barkthick.list  # see helper_fxns.R
  
  
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
  
  return(curr.df)
  
}

doAnalysis_traits_explain_pmr <- function(datasets, stem.respVars, traitVars.stem, cfract){
  
  #isolate trait data and scale it
  datasets %>%
    map(~ungroup(.x, codeStem)) %>%
    map(~mutate_at(.x, vars(traitVars.stem), list(s = scale2))) -> datasets.s # see helper_fxns.R
  
  #set up full models
  scaled.traitvars <- paste(traitVars.stem ,"s",sep = "_")
  if(cfract == F){
    rhsVars <- paste(c("size",scaled.traitvars))
  }else{
    rhsVars <- paste(c(scaled.traitvars))
  }
  rhs <- paste(rhsVars, collapse = " + ")
  print(rhs)
  lhs <- "curr.pmr"
  
  # fit models
  datasets.s %>%
    map(~fit.lm(y = lhs, rhs, data = .x)) -> mod.full.list  # see helper_fxns.R
  
  #do stepwise model selection
  map2(.x = mod.full.list, 
       .y = datasets.s, ~backward_selection(.x, data = .y)) -> mod.select.list  # see helper_fxns.R
  
  #extract and save the residuals
  mod.select.list %>%
    map(~residuals(.x)) %>%
    map2(.y = datasets.s, 
         ~data.frame(resid=.x, codeStem = .y$codeStem)) %>%
    bind_rows(.id = "resp") -> residuals
  
  #create prediction dataframes for each trait and model
  scaled.traitvars %>%
    purrr::set_names() %>%
    map(~preddat_fun_bysize_allmodels(models = mod.select.list,
                                      data = datasets.s, # see helper_fxns.R
                                      curr.traitVar = .x,
                                      data.list = T)) -> preddat
  
  #save model fit stats
  mod.select.list %>%
    map(~summary(.x)) %>%
    map(~data.frame(Fstat = round(.x$fstatistic['value'], digits=2),
                    numdf = .x$fstatistic['numdf'],
                    dendf = .x$fstatistic['dendf'],
                    r.squared = round(.x$r.squared, digits=2))) %>%
    list_to_df() %>%
    select(-source) %>%
    t() -> fitstats

  result.traitsPMR <- list(
    data = datasets,
    respVars = stem.respVars,
    models = mod.select.list,
    residuals = residuals,
    preddat = preddat,
    fitstats = fitstats)
  
  return(result.traitsPMR)
  
}

makefig__traits_explain_pmr <- function(result.list, traitVars.stem, cfract){
  
  #set up plotting dataframe
  plot.obj <- MakeLm_plottingDF(mod.list = result.list$models, respvars = result.list$respVars)
  #plot.obj$plotting.df
  #plot.obj$r2.df
  
  #define wood trait levels
  traitVars.stem.s <- c(paste(traitVars.stem, "s", sep = "_"))
  if(cfract == F){
    trait.levels <- c("sizesmall", traitVars.stem.s)
    trait.labels <- c("sizesmall", traitVars.stem)
  }else{
    trait.levels <- c(traitVars.stem.s)
    trait.labels <- c(traitVars.stem)
  }
  plot.obj$plotting.df$term <- factor(plot.obj$plotting.df$term, 
                                      levels = rev(trait.levels),
                                      labels = rev(trait.labels))
  
  #define response variable levels
  respvar.levels <- unlist(result.list$respVars)
  plot.obj$plotting.df %>%
    filter(respvar %in% respvar.levels) %>%
    mutate(respvar = factor(respvar, levels = respvar.levels)) -> plot.obj$plotting.df
  
  #plot
  plot.obj$plotting.df
  p <- ggplot(plot.obj$plotting.df, aes(x = respvar, y = term, fill = etasq)) +
    geom_tile(color = "black") + 
    geom_text(aes(label = round(est, digits = 2)), size = 3) +
    xlab("Response variable") + ylab("Wood trait") +
    theme_classic() +
    scale_fill_distiller(direction = 1) +
    scale_y_discrete(limits = rev(trait.labels))
  p
  
  result <- list(p = p, p.r2 = plot.obj$r2.df)
  return(result)
}

#--------------------------------------#
# code-level Cfractions



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
  
  return(curr.df)
  
}

