Calc_R2<-function(ne_fits_df){
  
  pred<-ne_fits_df[['predicted']]
  mass<-ne_fits_df[['mass']]
  
  sstot<-sum((mass - mean(mass))^2)
  ssres<-sum((pred - mass)^2)
  r2<-1-(ssres/sstot)
  
  return(r2)
}

fit_all_curves<-function(df_in){
  
  df_in<-plotting_df
  
  df_in %>%
    unite(SpeciesCode,col=sp_size,size,sep="_") -> df
  
  ne_fits <- lapply(split(df,factor(df$sp_size)),function(x){
    fit_litter(time = x$time/12, 
               mass.remaining = x$pmr, model = c("neg.exp"), iters = 500)
  })
  
  w.fits <- lapply(split(df,factor(df$sp_size)),function(x){
    fit_litter(time = x$time/12, 
               mass.remaining = x$pmr, model = c("weibull"), iters = 500)
  })
  
  spdf<-data.frame(k=unlist(lapply(ne_fits,
                                   function(x)x$optimFit$par)),
                   t70=unlist(lapply(ne_fits,function(x) time_to_prop_mass_remaining(x,threshold.mass=0.70))),
                   neg.exp.aic=unlist(lapply(ne_fits,
                                             function(x)x$fitAICc)),
                   ne.r2 = unlist(lapply(ne_fits,function(x)Calc_R2(x))),
                   w.t70=unlist(lapply(w.fits,function(x)             time_to_prop_mass_remaining(x,threshold.mass=0.70))),
                   w.aic=unlist(lapply(w.fits,
                                       function(x)x$fitAICc)),
                   alpha=unlist(lapply(w.fits,
                                       function(x)x$optimFit$par[[2]])),
                   w.r2 = unlist(lapply(w.fits,function(x)Calc_R2(x)))#check parameter order
                   
  )
  
  spdf$sp_size<-rownames(spdf)
  
  
  spdf %>%
    separate(sp_size,sep="_",into = c("species","size"))->spdf
  return(spdf)
}
