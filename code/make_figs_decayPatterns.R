
makefig__compare_t50_aic <- function(decayfits){
  
  #(1) compare t50
  decayfits %>%
    select(species, size, t50, w.t50) %>%
    gather(key = "variable", value ="value", -c(species,size)) %>%
    mutate(value.round = round(value, digits = 1)) -> plot.df
  plot.df$variable <- recode(plot.df$variable, "t50" = "neg.exp", "w.t50" = "weibull")
  
  plot.df %>%
    group_by(species) %>%
    summarize(val = mean(value)) %>%
    arrange(-val) -> species.order
  species.order$species
  
  p1 <- ggplot(plot.df, aes(x = variable, y = reorder(species, value), 
                      fill = value, label = value.round)) +
    geom_tile() +
    geom_text() +
    facet_grid(~size) +
    scale_fill_distiller(direction = 1) +
    theme_bw() +
    xlab("Years to 50% mass loss by model type") + ylab("") +
    guides(fill = F)
  p1
  #ggsave(filename = "output/figures/supplementary/compare_t50.pdf", width = 5, height = 6)
  
  #(2) compare aic
  decayfits %>%
    select(code, species, size, ne.aic, w.aic) %>%
    mutate(diff.aic = ne.aic - w.aic) %>%
    mutate(species = factor(species, levels = rev(species.order$species))) -> tmp
  
  p2 <- ggplot(tmp, aes(y = species, x = diff.aic, color = size)) +
    geom_point() +
    geom_vline(xintercept = 0, linetype = 2) +
    ylab(" ") + xlab("Delta AIC (>0 favors Weibull)") +
    theme_bw() + guides(color = F) +
    scale_color_manual(values = c("darkgray","black"))
  p2
  
  #ggsave(filename = "output/figures/supplementary/compare_aic.pdf", width = 5, height = 6)
  
  pdf("output/figures/supplementary/compare_aic_t50.pdf", width=8, height=4)
  grid.arrange(p1, p2, ncol = 2)
  dev.off()
  
}

makefig__decayfits <- function(decayfits){
  
  #plot the weibull fits
  decayfits %>%
    select(species, size, code, 
           alpha, beta,
           mean1, lower1, upper1, mean2, lower2, upper2, 
           w.t50, w.aic, w.r2) -> decayfits.w
  
  stemSamples <- load_stemSamples()
  stemSamples %>%
    select(code, Binomial) %>%
    unique() -> tmp
  decayfits.w %>%
    left_join(tmp) -> decayfits.w
  
  #order species by w.t60
  decayfits.w %>%
    select(Binomial, w.t50) %>%
    group_by(Binomial) %>%
    summarize(mean = mean(w.t50,na.rm=TRUE)) %>%
    arrange(mean)-> binom.order
  decayfits.w$Binomial <- factor(decayfits.w$Binomial, levels = as.character(binom.order$Binomial))
  
  # plot
  p.t50 <- ggplot(decayfits.w, aes(y = Binomial, x = w.t50, 
                                   color = size)) +
    geom_point() +
    scale_color_manual(values = c("darkgray","black")) +
    ylab("") + xlab("Years to 50% mass loss") +
    theme_bw() + guides(color = F)
  
  p.scale <- ggplot(decayfits.w, aes(y = species, x = beta, color = size)) +
    geom_point() +
    geom_errorbarh(aes(xmin = lower1, xmax = upper1)) +
    #geom_point(aes(x = mean1), shape = 2, alpha = .5, size = 3) +
    scale_color_manual(values = c("darkgray","black")) +
    ylab("") + xlab("Scale parameter") +
    theme_bw() + guides(color = F)
  
  p.shape <- ggplot(decayfits.w, aes(y = species, x = alpha, color = size)) +
    geom_point() +
    geom_errorbarh(aes(xmin = lower2, xmax = upper2)) +
    #geom_point(aes(x = mean2), shape = 2, alpha = .5, size = 3) +
    scale_color_manual(values = c("darkgray","black")) +
    ylab("") + xlab("Shape parameter") +
    theme_bw() + guides(color = F)
  
  pdf(file = "output/figures/maintext/decayfits.pdf", width = 8, height = 6)
  grid.arrange(p.t50,
               p.scale +
                 theme(axis.line=element_blank(),
                       axis.text.y=element_blank(),
                       axis.title.y=element_blank()), 
               p.shape + 
                 theme(axis.line=element_blank(),
                       axis.text.y=element_blank(),
                       axis.title.y=element_blank()),
               ncol=3, widths=c(1.75,1,1)
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
  woodSp.o <- levels(reorder(decayfits$species, -decayfits$w.t50))
  pmr.byStem.df$species <- factor(pmr.byStem.df$species, levels=rev(woodSp.o))
  
  #add underlying pmr data points
  pmr %>%
    select(codeStem, code, species, size, time, pmr) %>%
    rename('timeNum' = 'time') -> pmr.tmp
  
  p.pmrBystem <- ggplot(pmr.byStem.df, 
                        aes(x=timeNum/12, y=mean.pmr*100, 
                            group=codeStem,
                            color = size)) + 
    geom_line() + 
    geom_point(inherit.aes = F, data = pmr.tmp, 
               aes(x= timeNum/12, y = pmr*100, color = size), pch = 16, size = 1) +
    facet_wrap(~species) +
    xlab("Time (years)") + ylab("Mass remaining (%)") + 
    theme_classic() +
    scale_color_manual(name = "Size", values=c("darkgray","black")) +
    guides(color = F, linetype = F)
  p.pmrBystem
  
  ggsave(filename = "output/figures/maintext/pmr_byStem.pdf", plot = p.pmrBystem, width=6, height=6)
  
}
