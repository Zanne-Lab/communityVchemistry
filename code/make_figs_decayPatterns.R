
# plot heatmap table of mean t50 for each species+size given neg. exp and weibull model (Fig S2a)
makefig__compare_t50_tab <- function(decayfits){
  
  #(1) compare t50
  decayfits %>%
    select(Binomial, species, size, t50, w.t50) %>%
    gather(key = "variable", value ="value", -c(Binomial, species,size)) %>%
    mutate(value.round = round(value, digits = 1)) -> plot.df
  plot.df$variable <- recode(plot.df$variable, "t50" = "neg.exp", "w.t50" = "weibull")
  
  # order Binomial by mean t50
  decayfits %>%
    group_by(Binomial) %>%
    summarize(mean = mean(w.t50)) %>%
    arrange(mean) -> order.binom
  plot.df$Binomial <- factor(plot.df$Binomial, levels = order.binom$Binomial)
  
  # order size
  plot.df$size <- factor(plot.df$size, levels = c("small","large"))
  
  p1 <- ggplot(plot.df, aes(x = variable, y = Binomial, 
                            fill = value, label = value.round)) +
    geom_tile() +
    geom_text() +
    facet_grid(~size) +
    scale_fill_distiller(direction = 1) +
    theme_bw() +
    xlab("Years to 50% mass loss by model type") + ylab("") +
    guides(fill = F)
  p1
  
  return(p1)
}

# plot delta AIC for weibull vs neg exp model for each species+size (Fig S2b)
makefig__compare_t50_aic <- function(decayfits){
  
  #(2) compare aic
  decayfits %>%
    select(Binomial, code, species, size, ne.aic, w.aic) %>%
    mutate(diff.aic = ne.aic - w.aic) -> plot.df
  
  # order Binomial by mean t50
  decayfits %>%
    group_by(Binomial) %>%
    summarize(mean = mean(w.t50)) %>%
    arrange(mean) -> order.binom
  plot.df$Binomial <- factor(plot.df$Binomial, levels = order.binom$Binomial)
  
  
  p2 <- ggplot(plot.df, aes(y = Binomial, x = diff.aic, color = size)) +
    geom_point() +
    geom_vline(xintercept = 0, linetype = 2) +
    ylab(" ") + xlab("Delta AIC (>0 favors Weibull)") +
    theme_bw() + guides(color = F) +
    scale_color_manual(values = c("darkgray","black"))
  p2
  
  return(p2)
}

# plot summary of decay fits for each species+size (Fig 1)
makefig__decayfits <- function(decayfits){
  
  #plot the weibull fits
  decayfits %>%
    select(Binomial, species, size, code, 
           alpha, beta,
           lower1.beta, upper1.beta, lower2.alpha, upper2.alpha, 
           w.t50, w.aic, w.r2) -> decayfits.w
  
  # order Binomial by mean t50
  decayfits %>%
    group_by(Binomial) %>%
    summarize(mean = mean(w.t50)) %>%
    arrange(mean) -> order.binom
  decayfits.w$Binomial <- factor(decayfits.w$Binomial, levels = order.binom$Binomial)
  
  # plot
  p.t50 <- ggplot(decayfits.w, aes(y = Binomial, x = w.t50, color = size)) +
    geom_point() +
    scale_color_manual(values = c("darkgray","black")) +
    ylab("") + xlab("Years to 50% mass loss") +
    theme_bw() + guides(color = F)
  p.t50
  
  p.1beta <- ggplot(decayfits.w, aes(y = Binomial, x = beta, color = size)) +
    geom_point() +
    geom_errorbarh(aes(xmin = lower1.beta, xmax = upper1.beta)) +
    scale_color_manual(values = c("darkgray","black")) +
    ylab("") + xlab("Scale parameter") +
    theme_bw() + guides(color = F)
  p.1beta
  
  p.2alpha <- ggplot(decayfits.w, aes(y = Binomial, x = alpha, color = size)) +
    geom_point() +
    geom_errorbarh(aes(xmin = lower2.alpha, xmax = upper2.alpha)) +
    scale_color_manual(values = c("darkgray","black")) +
    ylab("") + xlab("Shape parameter") +
    theme_bw() + guides(color = F)
  p.2alpha
  
  p.fit <- ggplot(decayfits.w, aes(y = Binomial, x = w.r2, color = size)) +
    geom_point() +
    scale_color_manual(values = c("darkgray","black")) +
    ylab("") + xlab("Weibull model fit (R2)") +
    theme_bw() + guides(color = F)
  p.fit

  plot.list <- list(p.t50 = p.t50, 
                    p.1beta = p.1beta, 
                    p.2alpha = p.2alpha, p.fit = p.fit)
  return(plot.list)
  
}

# plot percent mass remaining over time, overlay model fits (Fig S3)
makefig__pmr_byStem <- function(pmr, decayfits, param.labels){
  
  # pmr
  pmr %>%
    select(codeStem, code, species, size, time, pmr) -> pmr.tmp
  
  # add species binomial to pmr dataframes
  decayfits %>%
    select(code, species, size, Binomial) %>%
    unique() -> indx
  pmr.tmp %>%
    left_join(indx) -> pmr.tmp
  
  # add parameter estimates
  decayfits %>%
    select(code, mean1.beta, mean2.alpha) -> params
  
  # y = exp(-(x/beta)^alpha) # weibull function
  # make predictions using fitted parameters
  max.time <- max(pmr$time)/12
  time.vec <- seq(0, max.time, 1e-01)
  weibull.fxn <- function(x, beta, alpha){
    y = exp(-(x/beta)^alpha)
    return(y)
  }
  
  n.fits <- dim(params)[1]
  y.pred.list <- list()
  for(i in 1:n.fits){
    df <- data.frame(x = time.vec)
    df$beta <- params[i,"mean1.beta"]$mean1.beta
    df$alpha <- params[i,"mean2.alpha"]$mean2.alpha
    df$y.pred <- weibull.fxn(x = df$x, 
                          beta = df$beta, 
                          alpha = df$alpha)
    y.pred.list[[i]] <- df
  }
  names(y.pred.list) <- params$code
  y.pred.df <- list_to_df(y.pred.list)
  y.pred.df %>%
    dplyr::rename('code'='source') %>%
    left_join(indx) -> y.pred.df

  #use decayfits to order the wood species by t50 to match species order in previous figure
  decayfits %>%
    group_by(Binomial) %>%
    summarize(mean = mean(w.t50)) %>%
    arrange(mean) -> order.binom
  pmr.tmp$Binomial <- factor(pmr.tmp$Binomial, levels=rev(order.binom$Binomial))
  y.pred.df$Binomial <- factor(y.pred.df$Binomial, levels = rev(order.binom$Binomial))
  
  # plot
  y.pred.df %>%
    select(Binomial, size, beta, alpha) %>%
    group_by(Binomial) %>%
    summarize(beta.label = paste(round(unique(beta), digits = 1), collapse =", "),
              alpha.label = paste(round(unique(alpha), digits = 1), collapse = ", ")) %>%
    mutate(beta.label = paste0("Beta = ", beta.label)) %>%
    mutate(alpha.label = paste0("Alpha = ", alpha.label)) -> labels
  
  if(param.labels == TRUE){
    p <- ggplot(pmr.tmp, aes(x=time/12, y=pmr*100, group=code, color = size)) + 
      geom_point(pch = 16, size = 1) +
      geom_line(data = y.pred.df, aes(x = x, y = y.pred*100, color = size), inherit.aes = F) +
      geom_text(inherit.aes = F, x = 1, y = 30, data = labels, 
                aes(label = beta.label), 
                color = 2, size = 3) +
      geom_text(inherit.aes = F, x = 1, y = 15, data = labels, 
                aes(label = alpha.label), 
                color = 2, size = 3) +
      facet_wrap(~Binomial) +
      xlab("Time (years)") + ylab("Mass remaining (%)") + 
      theme_classic() +
      scale_color_manual(name = "Size", values=c("darkgray","black")) +
      guides(color = F, linetype = F)
  }else{
    p <- ggplot(pmr.tmp, aes(x=time/12, y=pmr*100, group=code, color = size)) + 
      geom_point(pch = 16, size = 1) +
      geom_line(data = y.pred.df, aes(x = x, y = y.pred*100, color = size), inherit.aes = F) +
      # geom_text(inherit.aes = F, x = 1, y = 30, data = labels, 
      #           aes(label = beta.label), 
      #           color = 2, size = 3) +
      # geom_text(inherit.aes = F, x = 1, y = 15, data = labels, 
      #           aes(label = alpha.label), 
      #           color = 2, size = 3) +
      facet_wrap(~Binomial) +
      xlab("Time (years)") + ylab("Mass remaining (%)") + 
      theme_classic() +
      scale_color_manual(name = "Size", values=c("darkgray","black")) +
      guides(color = F, linetype = F)
  }
  

  return(p)
}

# plot scale vs shapefor each species+size
makefig__scaleVshape <- function(decayfits){
  
  #plot the weibull fits
  decayfits %>%
    select(Binomial, species, size, code, 
           alpha, beta,
           lower1.beta, upper1.beta, lower2.alpha, upper2.alpha, 
           w.t50, w.aic, w.r2) -> decayfits.w
  
  # order Binomial by mean t50
  decayfits %>%
    group_by(Binomial) %>%
    summarize(mean = mean(w.t50)) %>%
    arrange(mean) -> order.binom
  decayfits.w$Binomial <- factor(decayfits.w$Binomial, levels = order.binom$Binomial)
  
  # plot
  p <- ggplot(decayfits.w, aes(y = alpha, x = beta, color = w.r2, label = code)) +
    geom_point(aes(size = w.t50)) +
    geom_text(nudge_x = .1, nudge_y = .05) +
    #geom_errorbarh(aes(xmin = lower1.beta, xmax = upper1.beta)) +
    #geom_errorbar(aes(ymin = lower2.alpha, ymax = upper2.alpha)) +
    ylab("Shape parameter") + xlab("Scale parameter") +
    theme_bw()
  p
  #scale_color_manual(values = c("darkgray","black")) + guides(color = F)
  
  return(p)
}

# weibull behavior
plot_weibull <- function(pmr, decayfits){
  
  weibull.fxn <- function(x, beta, alpha){
    y = exp(-(x/beta)^alpha)
    return(y)
  }
  
  # time
  max.time <- max(pmr$time)/12
  time.vec <- seq(0, max.time, 1e-01)
  # beta
  #min.beta <- round(min(decayfits$beta), digits = 1)
  min.beta <- 2
  #max.beta <- round(max(decayfits$beta), digits = 1)
  max.beta <- 10
  #mean.beta <- round(mean(decayfits$beta), digits = 1)
  mean.beta <- 5
  
  # alpha
  #min.alpha <- round(min(decayfits$alpha), digits = 1)
  min.alpha <- 0.5
  #max.alpha <- round(max(decayfits$alpha), digits = 1)
  max.alpha <- 6
  #mean.alpha <- round(mean(decayfits$alpha), digits = 1)
  mean.alpha = 3

  #predictions
  beta.l.alpha.mean <- weibull.fxn(x = time.vec, beta = min.beta, alpha = mean.alpha)
  beta.h.alpha.mean<- weibull.fxn(x = time.vec, beta = max.beta, alpha = mean.alpha)
  beta.mean.alpha.l<- weibull.fxn(x = time.vec, beta = mean.beta, alpha = min.alpha)
  beta.mean.alpha.h<- weibull.fxn(x = time.vec, beta = mean.beta, alpha = max.alpha)
  beta.mean.alpha.mean <- weibull.fxn(x = time.vec, beta = mean.beta, alpha = mean.alpha)
  beta.l.alpha.h <- weibull.fxn(x = time.vec, beta = min.beta, alpha = max.alpha)
  beta.h.alpha.l <- weibull.fxn(x = time.vec, beta = max.beta, alpha = min.alpha)
  beta.l.alpha.l <- weibull.fxn(x = time.vec, beta = min.beta, alpha = min.alpha)
  beta.h.alpha.h <- weibull.fxn(x = time.vec, beta = max.beta, alpha = max.alpha)
  
  df <- data.frame(x = time.vec, 
                   beta.l.alpha.mean, beta.h.alpha.mean,
                   beta.mean.alpha.l, beta.mean.alpha.h,
                   beta.mean.alpha.mean,
                   beta.l.alpha.h, beta.h.alpha.l,
                   beta.l.alpha.l, beta.h.alpha.h)
  df %>%
    gather(key = "fit.type", value = "y", -c(x)) %>%
    separate(fit.type, into = c("drop","beta.param","drop2","alpha.param")) %>%
    select(-c(drop, drop2)) %>%
    mutate(beta.param = ifelse(beta.param == "l",  paste0("beta = ",min.beta),
                               ifelse(beta.param == "h", paste0("beta = ",max.beta),
                                      paste0("beta = ", mean.beta)))) %>%
    mutate(alpha.param = ifelse(alpha.param == "l",  paste0("alpha = ",min.alpha),
                               ifelse(alpha.param == "h", paste0("alpha = ",max.alpha),
                                      paste0("alpha = ", mean.alpha)))) %>%
    mutate(fit.label = paste0("beta = ",beta.param, ", alpha = ", alpha.param)) -> plot.df

  plot.df$beta.param <- factor(plot.df$beta.param, levels = c("beta = 2", "beta = 5", "beta = 10"))
  plot.df$alpha.param <- factor(plot.df$alpha.param, levels = c("alpha = 0.5", "alpha = 3", "alpha = 6"))
  
  
  p <- ggplot(plot.df, aes(x = x, y = y, label = fit.label)) +
    geom_line(size =  3) +
    facet_grid(alpha.param~ beta.param) +
    xlab("Time (yrs)") + ylab("Mass remaining (%)") +
    theme_bw()
  
  return(p)
  
}







