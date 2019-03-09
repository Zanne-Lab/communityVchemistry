
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
  
  ggsave(filename = "output/figures/maintext/pmr_byStem.pdf", plot = p.pmrBystem, width=6, height=6)
  
}