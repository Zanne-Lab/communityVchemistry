Calc_R2<-function(ne_fits_df){
  
  pred<-ne_fits_df[['predicted']]
  mass<-ne_fits_df[['mass']]
  
  sstot<-sum((mass - mean(mass))^2)
  ssres<-sum((pred - mass)^2)
  r2<-1-(ssres/sstot)
  
  return(r2)
}

fit_all_curves<-function(df_in){
  
  df_in %>%
    unite(SpeciesCode,col=sp_size,size,sep="_") -> df
  
  #negative expon fit
  ne_fits <- lapply(split(df,factor(df$sp_size)),function(x){
    fit_litter(time = x$time/12, 
               mass.remaining = x$pmr, model = c("neg.exp"), iters = 500)
  })
  
  #weibold fit
  w.fits <- lapply(split(df,factor(df$sp_size)),function(x){
    fit_litter(time = x$time/12, 
               mass.remaining = x$pmr, model = c("weibull"), iters = 500)
  })
  
  #create bootstrap distrib of k and get the 95% CI
  bootk<-lapply(ne_fits, function(x){
    
    bootmat<-bootstrap_parameters(x, nboot=500) #1st column are iteration of the param, right, but what is in the second col?
    meanboot<-mean(bootmat[,1])
    sdboot<-sd(bootmat[,1])
    seboot<-sdboot/sqrt(length(bootmat[,1]))
    upperboot<-meanboot+(1.96*seboot)
    lowerboot<-meanboot-(1.96*seboot)
    
    result<-list(upper=upperboot, lower=lowerboot)
    return(result)
  })
  tmp<-lapply(bootk, function(x) data.frame(upper=x[[1]], lower=x[[2]]))
  k.upper<-unlist(lapply(tmp, function(x) x[['upper']]))
  k.lower<-unlist(lapply(tmp, function(x) x[['lower']]))
  
  #create a 95% CI bootstrap distribution of t70
  #not sure how to do this...
  
  k<-unlist(lapply(ne_fits, function(x) x$optimFit$par))
  t70<-unlist(lapply(ne_fits, function(x) time_to_prop_mass_remaining(x,threshold.mass=0.70)))
  neg.exp.aic<-unlist(lapply(ne_fits, function(x) x$fitAICc))
  ne.r2<-unlist(lapply(ne_fits, function(x) Calc_R2(x)))
  w.t70<-unlist(lapply(w.fits, function(x) time_to_prop_mass_remaining(x,threshold.mass=0.70)))
  w.aic<-unlist(lapply(w.fits, function(x) x$fitAICc))
  alpha<-unlist(lapply(w.fits, function(x) x$optimFit$par[[2]]))
  w.r2<-unlist(lapply(w.fits, function(x) Calc_R2(x)))
  
  spdf<-data.frame(k=k,
                   k.upper=k.upper,
                   k.lower=k.lower,
                   t70=t70,
                   ne.r2 = ne.r2,
                   alpha = alpha)
  spdf$sp_size<-rownames(spdf)
  
  spdf %>%
    separate(sp_size,sep="_",into = c("species","size"))->spdf
  
  return(spdf)
}
