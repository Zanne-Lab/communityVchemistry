# source a bunch of files from the same folder; code taken from R help -> source()
sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}

# turn a list of dataframes into a long dataframe
list_to_df <- function(mylist){
  
  # make a vector of row ids that correspond to the list names
  rowid.indx <- lapply(mylist, function(x) dim(x)[1])
  sourceVec.list <- list()
  for(i in 1:length(rowid.indx)){
    sourceName <- names(rowid.indx)[i]
    numRows <- rowid.indx[[i]]
    sourceVec.list[[i]] <- rep(sourceName, numRows)
  }
  rowVec <- unlist(sourceVec.list)
  
  # combine into df
  df <- data.frame(do.call(rbind, mylist))
  df$source <- rowVec
  
  return(df)
}

# scale and center
scale2 <- function(x, na.rm = FALSE){
  x.s <- (x - mean(x, na.rm = na.rm)) / sd(x, na.rm)
}

# fit lm using strings
fit.lm <- function(y, rhs, data){
  string<-paste(y, " ~ ", rhs)
  fmla <- as.formula(string)
  mod.full <- lm(formula = fmla, data = data)
  return(mod.full)
}

# backward selection on lm
backward_selection <- function(x, data){
  
  x.updated <- update(x, . ~ ., data = data)
  mod.select <- stats::step(x.updated, direction="backward", trace = 0)
  final.mod <- lm(formula(mod.select), data = data)
  
  return(final.mod)
}

# create prediction datasets for lm added variable plots using a dataset, add the unscale curr.var to the df
preddat_fun_bysize <- function(data, allvars, curr.var){
  
  #identify the response and predictor variables, and categorical predictors
  respvar <- allvars[1] # this is the response variable
  xvars <- allvars[-1]
  cont.vars <- xvars[!xvars %in% "size"]
  curr.var.unscaled <- strsplit(curr.var,"_")[[1]][1]
  
  #find the mean of each col other than the response var (var)
  data %>%
    group_by(size) %>%
    summarise_at(vars(cont.vars), mean, na.rm = T) %>%
    select(-curr.var) -> sums
  
  # select the response var and add colums that just have predictor means
  data %>%
    select(respvar, size, curr.var, curr.var.unscaled) %>%
    left_join(sums) -> result
  
  return(result)
}

# create prediction datasets and get new fitted values using all models in list
preddat_fun_bysize_allmodels <- function(models, data, curr.traitVar, data.list){
  
  # create prediction data
  if(data.list == T){
    models %>%
      map(~ all.vars(formula(.x))) %>%
      keep(~ curr.traitVar %in% .x) -> curr.models.vars
    
    names(data) %in% names(curr.models.vars)
    curr.data <- data[names(data) %in% names(curr.models.vars)]
    curr.data
    map2(.x = curr.models.vars, .y = curr.data, 
         ~ preddat_fun_bysize(data = .y, allvars = .x, 
                              curr.var = curr.traitVar)) -> preddat.list
  }else{
    models %>%
      map(~ all.vars(formula(.x))) %>%
      keep(~ curr.traitVar %in% .x) %>%
      map(~ preddat_fun_bysize(data = data, allvars = .x, 
                               curr.var = curr.traitVar)) -> preddat.list
  }
  if(length(preddat.list) != 0){ # if the current trait shows up in any of the final models...
    
    #identify the unscaled current traitVar
    curr.traitVar.unscaled <- strsplit(curr.traitVar,"_")[[1]][1]
    
    # use the models and new data to predict
    kept.models <- models[names(models) %in% names(preddat.list)]
    map2(.x = kept.models, 
         .y = preddat.list, ~augment(.x, newdata = .y)) %>%
      map(~mutate(.x, 
                  lower = .fitted - 2*.se.fit,
                  upper = .fitted + 2*.se.fit)) %>%
      map(~select(.x, 1, size, curr.traitVar, curr.traitVar.unscaled, .fitted, lower, upper)) %>%
      map(~rename(.x, 'respval' = 1)) %>%
      list_to_df() %>%
      rename('respvar'='source') -> preddat.fitted
    
    # add data for response var if the curr.traitVar was not ultimately included in the final model
    rm.models <- models[!names(models) %in% names(preddat.list)]
    missing.respvars <- names(rm.models)
    missing.respvars
    # if there are these models, add data
    if(length(missing.respvars) != 0){
      if(data.list == T){
        missing.data <- data[names(data) %in% missing.respvars]
        missing.data %>%
          map(~select(.x, size, curr.pmr, curr.traitVar, curr.traitVar.unscaled)) %>%
          list_to_df() %>%
          rename('respvar'='source') %>%
          rename('respval'='curr.pmr') %>%
          mutate(.fitted = NA) %>%
          mutate(lower = NA) %>%
          mutate(upper = NA) %>%
          select(respval, size, curr.traitVar, curr.traitVar.unscaled, .fitted, lower, upper, respvar) -> missing.dat
      }else{
        data %>%
          select(size, missing.respvars, curr.traitVar, curr.traitVar.unscaled) %>%
          gather(key = "respvar",value = "respval", -c(size, curr.traitVar, curr.traitVar.unscaled)) %>%
          mutate(.fitted = NA) %>%
          mutate(lower = NA) %>%
          mutate(upper = NA) %>%
          select(respval, size, curr.traitVar, curr.traitVar.unscaled, .fitted, lower, upper, respvar) -> missing.dat
      }
      preddat.fitted <- rbind(preddat.fitted, missing.dat)
    }
    
  }else{
    preddat.fitted <- data.frame(NA)
  }
  
  
  
  return(preddat.fitted)
}




