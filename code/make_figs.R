# make figures

## decayPatterns

makefig__compare_t60 <- function(decayfits){
  
  decayfits %>%
    select(Binomial, size, t60, w.t60) %>%
    gather(key = "variable", value ="value", -c(Binomial,size)) %>%
    mutate(value.round = round(value, digits = 1)) -> plot.df
  plot.df$variable <- recode(plot.df$variable, "t60" = "neg.exp", "w.t60" = "weibull")
  
  ggplot(plot.df, aes(x = variable, y = reorder(Binomial, value), 
                      fill = value, label = value.round)) +
    geom_tile() +
    geom_text() +
    facet_grid(~size) +
    scale_fill_distiller(direction = 1) +
    theme_bw() +
    xlab("Years to 40% mass loss") + ylab("") +
    guides(fill = F)
  
  ggsave(filename = "output/figures/supplementary/compare_t60.pdf", width = 5, height = 6)
  
}

makefig__compare_aic <- function(decayfits){
  
  decayfits %>%
    select(code, Binomial, size, ne.aic, w.aic) %>%
    mutate(diff.aic = ne.aic - w.aic) -> tmp
  
  ggplot(tmp, aes(y = reorder(Binomial, diff.aic), x = diff.aic, color = size)) +
    geom_point() +
    geom_vline(xintercept = 0, linetype = 2) +
    ylab("Code") + xlab("Delta AIC (>0 favors Weibull)") +
    theme_bw()
  
  ggsave(filename = "output/figures/supplementary/compare_aic.pdf", width = 5, height = 6)
  
}

makefig__decayfits <- function(decayfits){
  
  #plot the weibull fits
  decayfits %>%
    select(Binomial, size, code, 
           alpha, beta,
           mean1, lower1, upper1, mean2, lower2, upper2, 
           w.t60, w.aic, w.r2) -> decayfits.w
  
  #order Binomial by w.t60
  decayfits.w %>%
    select(Binomial, w.t60) %>%
    group_by(Binomial) %>%
    summarize(mean = mean(w.t60)) %>%
    arrange(mean)-> binom.order
  decayfits.w$Binomial <- factor(decayfits.w$Binomial, levels = as.character(binom.order$Binomial))
  
  # plot
  p.t60 <- ggplot(decayfits.w, aes(y = Binomial, x = w.t60, 
                                   color = size)) +
    geom_point() +
    scale_color_manual(values = c("black","darkgray")) +
    ylab("") + xlab("Years to 40% mass loss") +
    theme_bw() + guides(color = F)
  
  p.scale <- ggplot(decayfits.w, aes(y = Binomial, x = beta, color = size)) +
    geom_point() +
    geom_errorbarh(aes(xmin = lower1, xmax = upper1)) +
    #geom_point(aes(x = mean1), shape = 2, alpha = .5, size = 3) +
    scale_color_manual(values = c("black","darkgray")) +
    ylab("") + xlab("Scale parameter") +
    theme_bw() + guides(color = F)
  
  p.shape <- ggplot(decayfits.w, aes(y = Binomial, x = alpha, color = size)) +
    geom_point() +
    geom_errorbarh(aes(xmin = lower2, xmax = upper2)) +
    #geom_point(aes(x = mean2), shape = 2, alpha = .5, size = 3) +
    scale_color_manual(values = c("black","darkgray")) +
    ylab("") + xlab("Shape parameter") +
    theme_bw() + guides(color = F)

  pdf(file = "output/figures/maintext/decayfits.pdf", width = 8, height = 6)
  grid.arrange(p.t60,
               p.scale +
                 theme(axis.line=element_blank(),
                       axis.text.y=element_blank(),
                       axis.title.y=element_blank()), 
               p.shape + 
                 theme(axis.line=element_blank(),
                       axis.text.y=element_blank(),
                       axis.title.y=element_blank()),
               ncol=3, widths=c(2,1,1)
  )
  dev.off()
  
}

makefig__pmr_byStem <- function(pmr_byStem){
  
  #make a long version of pmr.byStem.df.w
  pmr_byStem %>%
    gather(key="time",value="mean.pmr", 2:7, na.rm=TRUE) %>%
    separate(time, into=c("drop","timeNum"), sep="time") %>%
    mutate(timeNum=as.numeric(timeNum)) %>%
    mutate(species=tolower(code)) %>%
    mutate(size = case_when(code == species ~ "small", 
                            code != species ~ "large")) %>%
    select(-drop) -> pmr.byStem.df
  
  #order the wood species by t70 to match species order in previous figure
  unique(pmr.byStem.df$timeNum)
  woodSp.o <- levels(reorder(decayfits$species, -decayfits$w.t60))
  pmr.byStem.df$species <- factor(pmr.byStem.df$species, levels=rev(woodSp.o))
  
  #add underlying pmr data points
  pmr %>%
    select(codeStem, code, species, size, time, pmr) %>%
    rename('timeNum' = 'time') -> pmr.tmp
  
  p.pmrBystem <- ggplot(pmr.byStem.df, 
                        aes(x=timeNum/12, y=mean.pmr*100, 
                            group=codeStem,
                            linetype=size)) + 
    geom_line() + 
    geom_point(inherit.aes = F, data = pmr.tmp, 
               aes(x= timeNum/12, y = pmr*100), 
               alpha = .2, pch = 16) +
    facet_wrap(~species) +
    xlab("Time (years)") + ylab("Mass remaining (%)") + 
    theme_bw() +
    scale_linetype_manual(values=c(2,1)) + 
    guides(linetype=FALSE)
  
  pdf("output/figures/maintext/pmr_byStem.pdf", width=6, height=6)
  p.pmrBystem
  dev.off()
  
}

## woodTraits_explainDecay.Rmd

makefig__traits_explain_decayParams <- function(decayfits, traits.code){
  
  #merge traits and response variables into 1 df
  decayfits %>%
    select(code, species, size, 
           ne.r2, k, t60, 
           w.r2, alpha, beta, w.t60) %>%
    left_join(traits.code) %>% 
    filter(!is.na(P)) %>%
    filter(!is.na(waterperc)) -> decayfits.traits
  
  #set up full models
  respVars <- list("k","t60","ne.r2", "alpha","beta","w.t60","w.r2")
  rhs <- "size + waterperc + density + barkthick + P + K + Ca + Mn + Fe + Zn + N + C"
  mod.full.list<-lapply(respVars, ModelFit_manyYs, rhs, curr.data=decayfits.traits) # summaryTable_fxns.R
  
  #do stepwise model selection #uncomment if data changes
  # mod.select.list<-lapply(mod.full.list, function(x) {
  #   x.updated<-update(x, . ~ ., data = model.frame(x))
  #   mod.select<-step(x.updated, direction="backward")
  #   return(mod.select)
  #   })
  # names(mod.select.list) <- respVars
  # #saveRDS(mod.select.list, file = "derived_data/modSelect.RData")
  mod.select.list <- readRDS(file = "derived_data/modSelect.RData")
  
  #make table
  #prettyTab <- MakeLmSummaryTable(respvars=unlist(respVars), mod.list=mod.select.list) # summaryTable_fxns.R
  #write.csv(prettyTab, "output/traits_explain_decayparams.csv")
  
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

makefig__variation_densityNbarkthick <- function(traits.stem, traits.code, pmr_byStem){
  
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

makefig__traits_explain_pmr <- function(traits.stem, traits.code, pmr_byStem){
  
  respVars <- list("time7", "time13", "time25", "time37","time59")
  
  # create dataframes
  # stem-level waterperc, xrf, and CN and (small) species-level barkthickness and density
  datasets<-lapply(respVars, function(x) {
    result<-CreateTraitPMRpair(x, traits.stem, traits.code, pmr_byStem) #analysisDF_fxns.R
    return(result)
  })
  names(datasets)<-respVars
  
  #set up full models
  rhs <- "size + waterperc + density_smspp + barkthick_smspp + P + K + Ca + Mn + Fe + Zn + N + C"
  lhs <- "curr.pmr"
  mod.full.list<-lapply(datasets, function(x) {
    result<-ModelFit_manyYs(y=lhs, rhs=rhs, curr.data=x)
    return(result)
  })
  
  # #do stepwise model selection
  # mod.select.list<-lapply(mod.full.list, function(x) {
  #   x.updated<-update(x, . ~ ., data = model.frame(x)) 
  #   mod.select<-step(x.updated, direction="backward")
  #   return(mod.select)
  #   })
  # saveRDS(mod.select.list, file = "derived_data/modSelect_stem.RData")
  mod.select.list <- readRDS(file = "derived_data/modSelect_stem.RData")
  
  #summarize the best models
  #prettyTab <- MakeLmSummaryTable(respvars=unlist(respVars), mod.list=mod.select.list) # summaryTable_fxns.R
  #write.csv(prettyTab, "output/traitsummary_stem.csv")
  
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