fit_all_curves<-function(df_in){
  df_in %>%
    unite(SpeciesCode,col=sp_size,size,sep="_") -> df
  
  ne_fits <- lapply(split(df,factor(df$sp_size)),function(x){
    fit_litter(time = x$time/12, 
               mass.remaining = x$pmr, model = c("neg.exp"), iters = 1000)
  })
  
  spdf<-data.frame(k=unlist(lapply(ne_fits,function(x)x$optimFit$par)))
  spdf$sp_size<-row.names(spdf)
  
  spdf %>%
    separate(sp_size, c("species", "size"))%>%
    ggplot(aes(x=size,y=k,col=species))+geom_point()
  
  spdf$sp_size<-row.names(spdf)
  
  
  w.fits <- lapply(split(df,factor(df$sp_size)),function(x){
    fit_litter(time = x$time/12, 
               mass.remaining = x$pmr, model = c("weibull"), iters = 1000)
  })
  
  spdf<-data.frame(k=unlist(lapply(ne_fits,
                                   function(x)x$optimFit$par)),
                   t70=unlist(lapply(ne_fits,function(x)             time_to_prop_mass_remaining(x,threshold.mass=0.70))),
                   neg.exp.aic=unlist(lapply(ne_fits,
                                             function(x)x$fitAICc)),
                   w.t70=unlist(lapply(w.fits,function(x)             time_to_prop_mass_remaining(x,threshold.mass=0.70))),
                   w.aic=unlist(lapply(w.fits,
                                       function(x)x$fitAICc)),
                   alpha=unlist(lapply(w.fits,
                                       function(x)x$optimFit$par[[2]]))#check parameter order
                   
  )
  
  spdf$sp_size<-rownames(spdf)
  
  
  spdf %>%
    separate(sp_size,sep="_",into = c("species","size"))->spdf
  return(spdf)
}