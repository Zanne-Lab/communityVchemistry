#--------------------------------------#
# use code-level traits to explain decay

#updated 7.28.20
doAnalysis_traits_explain_decayParams <- function(decayfits, traits.code, 
                                                  code.respVars, traitVars, slice){
  
  # 1- set up df
  #isolate trait data and scale it
  traits.code %>%
    select(c("code","species","size",traitVars)) -> select.traits.code
  if(slice == 2){
    select.traits.code %>%
      filter(size == "large") -> select.traits.code
  }
  select.traits.code %>%
    mutate_at(vars(traitVars), list(s = scale2)) -> select.traits.code
  select.traits.code
  #isolate decay data
  decayfits %>%
    select(c("code","species","size", unlist(code.respVars))) -> select.decayfits
  if(slice == 2){
    select.decayfits %>%
      filter(size == "large") -> select.decayfits
  }
  #merge traits and response variables into 1 df
  select.decayfits %>%
    left_join(select.traits.code) -> decayfits.traits
  
  #set up full models
  scaled.traitvars <- paste(traitVars,"s",sep = "_")
  if(slice == 1){ # slice1 = no C fraction data, include samll
    rhsVars <- paste(c("size",scaled.traitvars))
  }else{
    rhsVars <- paste(c(scaled.traitvars))
  }
  rhs <- paste(rhsVars, collapse = " + ")
  print(rhs)
  
  # 2 - fit cv glmnet
  require(glmnet)
  var.list <- list()
  cvfit.list <- list()
  code.respVars
  i<-1
  for(i in 1:length(code.respVars)){
    y <- decayfits.traits[,code.respVars[[i]]]
    x <- data.frame(decayfits.traits[,rhsVars], stringsAsFactors = F)
    if(slice == 1){ # slice1 = no C fraction data, include samll
      x %>%
        mutate(size = ifelse(size == "small", 0, 1)) -> x
    }
    x <- as.matrix(x)
    y <- as.matrix(y)
    #fit <- glmnet(x = x, y = y, family="gaussian", standardize = F)
    #plot(fit, xvar = "lambda", label = T)
    cvfit <- cv.glmnet(x = x, y = y, family = "gaussian", standardize = F,
                       nfolds = nrow(y), grouped = FALSE)
    cvfit
    vars <- extract.lambda_uni(cvfit, s = "lambda.min")
    if(is.null(vars)){
      var.list[[i]] <- "None"
    }else{
      var.list[[i]] <- vars
      cvfit.list[[i]] <- cvfit
    }
  }
  var.list
  cvfit.list
  
  # 3 - fit cv glmnet final model
  # if there were any empty models, add back all the variables
  # var.list
  var.list.tmp <- var.list
  selection <- var.list %in% "None"
  length(selection)
  i<-1
  for(i in 1:length(selection)){
    if(selection[i] == TRUE){
      var.list.tmp[[i]] <- paste0(traitVars,"_s")
    }
  }
  var.list.tmp
  
  mod.list <- list()
  i<-1
  for(i in 1:length(code.respVars)){
    y <- decayfits.traits[,code.respVars[[i]]]
    x <- data.frame(decayfits.traits[,var.list.tmp[[i]]], stringsAsFactors = F)
    if("size" %in% colnames(x)){
      x %>%
        mutate(size = ifelse(size == "small", 0, 1)) -> x
    }
    x <- as.matrix(x)
    y <- as.matrix(y)
    fml<- as.formula(paste("y~",paste(colnames(x),collapse = "+")))
    fml
    df <- data.frame(y, x, stringsAsFactors = F)
    mod <- lm(fml, data = df)
    mod.list[[i]] <- mod
  }
  names(mod.list) <- code.respVars
  
  # 4 - extract and save the residuals from final model
  mod.list %>%
    map(~residuals(.x)) %>%
    map(~data.frame(resid=.x, code = decayfits.traits$code)) %>%
    list_to_df() %>%
    dplyr::rename('resp'='source') -> residuals
  
  # 5 - save model fit stats
  mod.list %>%
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
    models = mod.list,
    vars = var.list,
    residuals = residuals,
    fitstats = fitstats)
  
  return(result.traitsDecay)
  
}

#--------------------------------------#
# use stem-level traits to expain decay

#updated 7.28.20
doAnalysis_traits_explain_pmr <- function(stem.respVars, traits.stem, pmr_byStem, traitVars, slice){
  
  # 1- set up df
  #isolate trait data and scale it
  traits.stem %>%
    select(c("codeStem","code","species","size", traitVars)) -> select.traits.stem
  if(slice == 2){
    select.traits.stem %>%
      filter(size == "large") -> select.traits.stem
  }
  select.traits.stem %>%
    mutate_at(vars(traitVars), list(s = scale2)) -> select.traits.stem
  #isolate decay data
  pmr_byStem %>%
    select(c("codeStem", "code","species","size", unlist(stem.respVars))) -> select.pmr
  if(slice == 2){
    select.pmr %>%
      filter(size == "large") -> select.pmr
  }
  #merge traits and response variables into 1 df
  select.pmr %>%
    left_join(select.traits.stem) -> pmr.traits
  
  #set up full models
  scaled.traitvars <- paste(traitVars,"s",sep = "_")
  if(slice == 1){ # slice1 = no C fraction data, include samll
    rhsVars <- paste(c("size",scaled.traitvars))
  }else{
    rhsVars <- paste(c(scaled.traitvars))
  }
  rhs <- paste(rhsVars, collapse = " + ")
  print(rhs)
  
  # 2 - fit cv glmnet
  require(glmnet)
  var.list <- list()
  cvfit.list <- list()
  stem.respVars
  i<-1
  for(i in 1:length(stem.respVars)){
    y <- pmr.traits[,stem.respVars[[i]]]
    x <- data.frame(pmr.traits[,rhsVars], stringsAsFactors = F)
    if(slice == 1){ # slice1 = no C fraction data, include samll
      x %>%
        mutate(size = ifelse(size == "small", 0, 1)) -> x
    }
    x <- as.matrix(x)
    y <- as.matrix(y)
    data <- cbind(x, y)
    data.nonas <- data[complete.cases(data),]
    x.nonas <- data.nonas[,colnames(data.nonas) %in% colnames(x)]
    y.nonas <- data.nonas[,!colnames(data.nonas) %in% colnames(x)]
    #fit <- glmnet(x = x.nonas, y = y.nonas, family="gaussian")
    #plot(fit, xvar = "lambda", label = T)
    cvfit <- cv.glmnet(x = x.nonas, y = y.nonas, family = "gaussian", standardize = F,
                       nfolds = nrow(y), grouped = FALSE)
    #plot(cvfit)
    vars <- extract.lambda_uni(cvfit, s = "lambda.min")
    if(is.null(vars)){
      var.list[[i]] <- "None"
    }else{
      var.list[[i]] <- vars
      cvfit.list[[i]] <- cvfit
    }
  }
  var.list
  cvfit.list
  
  # 3 - fit cv glmnet final model
  # if there were any empty models, add back all the variables
  # var.list
  var.list.tmp <- var.list
  selection <- var.list %in% "None"
  selection
  i<-1
  for(i in 1:length(selection)){
    if(selection[i] == TRUE){
      var.list.tmp[[i]] <- paste0(traitVars,"_s")
    }
  }
  var.list.tmp
  
  mod.list <- list()
  samples.list <- list()
  i<-1
  for(i in 1:length(stem.respVars)){
    y <- pmr.traits[,stem.respVars[[i]]]
    x <- data.frame(pmr.traits[,var.list.tmp[[i]]], stringsAsFactors = F)
    if("size" %in% colnames(x)){
      x %>%
        mutate(size = ifelse(size == "small", 0, 1)) -> x
    }
    x <- as.matrix(x)
    y <- as.matrix(y)
    data <- cbind(x, y)
    data.nonas <- data[complete.cases(data),]
    fml <- paste(stem.respVars[[i]], "~", paste0(var.list.tmp[[i]], collapse = "+"))
    df <- data.frame(data.nonas, stringsAsFactors = F)
    mod <- lm(fml, data = df)
    mod.list[[i]] <- mod
    samples.list[[i]] <- pmr.traits[complete.cases(data),"codeStem"]
    
  }
  names(mod.list) <- stem.respVars
  names(samples.list) <- stem.respVars
  
  # 4 - extract and save the residuals from final model
  mod.list %>%
    map(~residuals(.x)) %>%
    map2(.y = samples.list, .f = ~data.frame(resid=.x, codeStem =.y)) %>%
    list_to_df() %>%
    dplyr::rename('resp'='source') -> residuals
  
  # 5 - save model fit stats
  mod.list %>%
    map(~summary(.x)) -> summary.list
  summary.list %>%
    map(~is.null(.x$fstatistic)) -> empty
  not.empty <- summary.list[!unlist(empty)]
  not.empty %>%
    map(~data.frame(Fstat = round(.x$fstatistic['value'], digits=2),
                    numdf = .x$fstatistic['numdf'],
                    dendf = .x$fstatistic['dendf'],
                    r.squared = round(.x$r.squared, digits=2))) %>%
    list_to_df() %>%
    select(-source) %>%
    t() -> fitstats
  fitstats
  if(sum(unlist(empty)) != 0){
    empty.mods <- summary.list[unlist(empty)]
    empty.mods %>%
      map(~data.frame(Fstat = NA,
                      numdf = .x$df[1],
                      dendf = .x$df[2],
                      r.squared = round(.x$r.squared, digits=2))) %>%
      list_to_df() %>%
      select(-source) %>%
      t() -> fitstats.empty
    fitstats <- cbind(fitstats, fitstats.empty)
  }
  
  result.traitsDecay <- list(
    data = pmr.traits,
    respVars = stem.respVars,
    models = mod.list,
    vars = var.list,
    residuals = residuals,
    fitstats = fitstats)
  
  
  return(result.traitsDecay)
  
}



#--------------------------------------#
# rationalize aggregation of trait data

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




