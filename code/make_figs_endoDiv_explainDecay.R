
summarize_endoDiv <- function(taxAndFunguild, comm.otu){
  
  # summarize the diversity in each sample
  rich.df<-Calc_richOTU(taxAndFunguild, comm.otu)
  H.df<-Calc_H.OTU(taxAndFunguild, comm.otu)
  sapro.df<-Calc_richOTUtype(colNam="Trophic.Mode", grepTerm="Sapro", taxAndFunguild, comm.otu=comm.otu)  # analysisDF_fxns.R
  basidio.df<-Calc_richOTUtype(colNam="phylum", grepTerm="Basid", taxAndFunguild, comm.otu=comm.otu)
  path.df<-Calc_richOTUtype(colNam="Trophic.Mode", grepTerm="Patho", taxAndFunguild, comm.otu=comm.otu)
  oomy.df<-Calc_richOTUtype(colNam="kingdom", grepTerm="Protist", taxAndFunguild, comm.otu=comm.otu) # have to do this with the full community matrix because there are so few of these guys
  richList<-list(rich.df, H.df, sapro.df, basidio.df, path.df, oomy.df)
  richNames<-c("Richness","ShannonsH","Sapro.rich","Basidio.rich","Patho.rich","Oomy.rich")
  names(richList)<-richNames
  
  return(richList)
}

doAnalysis_endoDiv_explainDecay <- function(taxAndFunguild, comm.otu, decayfits, code.respVars){
  
  # taxAndFunguild
  # comm.otu
  # decayfits
  # code.respVars
  
  richList <- summarize_endoDiv(taxAndFunguild, comm.otu)
  
  #create a set of richness and decay dataframes
  rich.decayfits<-lapply(richList, function(x) {CreateRichDecayfits_df(otutype.df=x, decayfits=decayfits)})
  
  # fit models
  rhs <- "size * mean"
  # size = small/large
  # mean = mean OTU richness
  lhs <- code.respVars
  
  richType.list<-list()
  for(i in 1:length(rich.decayfits)){
    resp.list<-list()
    for(k in 1:length(lhs)){
      resp.list[[k]]<-ModelFit_manyYs(y = lhs[[k]], rhs = rhs, curr.data = rich.decayfits[[i]]) 
    }
    names(resp.list)<-lhs
    richType.list[[i]]<-resp.list
  }
  names(richType.list) <- names(richList)
  
  #summarize
  richType.list
  pretty.list<-lapply(richType.list, MakeLmSummaryTable, respvars=lhs)
  pretty.df <- list_to_df(pretty.list)
  
  # #identify significant relationships
  # pretty.df %>%
  #   filter(term %in% c('mean','sizesmall','sizesmall:mean')) %>%
  #   gather(key = "respvar", value = "coef", -c(source, term)) %>%
  #   mutate(signif = ifelse(grepl("*", coef, fixed = T), T, F)) %>%
  #   filter(signif == T)
  
  result.list <- list(rich.decayfits = rich.decayfits,
                      pretty.df = pretty.df)
  
  return(result.list)
}

doAnalysis_endoDiv_explainPMR <- function(taxAndFunguild, comm.otu, pmr_byStem, stem.respVars){
  
  richList <- summarize_endoDiv(taxAndFunguild, comm.otu)
  
  #create a set of richness and decay dataframes
  rich.decayfits<-lapply(richList, function(x) {
    left_join(pmr_byStem, x) %>% 
      filter(!is.na(seq_sampName)) -> result
    return(result)
  })
  
  # fit models
  rhs <- "size * sub_rich"
  lhs <- stem.respVars
  
  richType.list<-list()
  for(i in 1:length(rich.decayfits)){
    resp.list<-list()
    for(k in 1:length(lhs)){
      resp.list[[k]]<-ModelFit_manyYs(y = lhs[[k]], rhs = rhs, curr.data = rich.decayfits[[i]]) 
    }
    names(resp.list)<-lhs
    richType.list[[i]]<-resp.list
  }
  names(richType.list)<-names(richList)
  
  #summarize
  pretty.list<-lapply(richType.list, MakeLmSummaryTable, respvars=lhs)
  pretty.df <- list_to_df(pretty.list)
  
  # #identify significant relationships
  # pretty.df %>%
  #   filter(term %in% c('sizesmall','sizesmall:sub_rich','sub_rich')) %>%
  #   gather(key = "respvar", value = "coef", -c(source, term)) %>%
  #   mutate(signif = ifelse(grepl("*", coef, fixed = T), T, F)) %>%
  #   filter(signif == T)
  
  result.list <- list(rich.decayfits = rich.decayfits,
       pretty.df = pretty.df)
  
  return(result.list)
  
}

makefig__RichShannonPMR.stem <- function(taxAndFunguild, comm.otu, pmr_byStem, stem.respVars){
  
  result.list.stem <- doAnalysis_endoDiv_explainPMR(taxAndFunguild, comm.otu, pmr_byStem, stem.respVars)
  names(result.list.stem$rich.decayfits)
  
  tmp <- result.list.stem$rich.decayfits[['Richness']]
  tmp %>%
    select(codeStem, time25, time37, sub_rich, size) %>%
    gather(key = "time", value = "pmr", -c(codeStem, sub_rich, size)) -> tmp1
  p.rich<-ggplot(tmp1, aes(x = sub_rich, y = pmr)) +
    geom_point() +
    geom_smooth(method = "lm") +
    facet_grid(time~size)+
    xlab("OTU richness") + ylab("Mass remaining after X months (%)") +
    theme_bw()
p.rich

tmp <- result.list.stem$rich.decayfits[['ShannonsH']]
tmp %>%
  select(codeStem, time25, time37, sub_rich, size) %>%
  gather(key = "time", value = "pmr", -c(codeStem, sub_rich, size)) -> tmp1
p.h<-ggplot(tmp1, aes(x = sub_rich, y = pmr)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_grid(time~size)+
  xlab("Shannon's H") + ylab("Mass remaining after X months (%)") +
  theme_bw()
p.h

pdf(file = "output/figures/supplementary/richness_pmr.pdf", width = 10, height = 5)
grid.arrange(p.rich, p.h, ncol = 2)  
dev.off()
 
}

makefig__pathoRichPMR.stem <- function(taxAndFunguild, comm.otu, pmr_byStem, stem.respVars){
  
  result.list.stem <- doAnalysis_endoDiv_explainPMR(taxAndFunguild, comm.otu, pmr_byStem, stem.respVars)
  
  tmp <- result.list.stem$rich.decayfits[['Patho.rich']]
  p <- ggplot(tmp, aes(x = sub_rich, y = time37)) +
    geom_point() +
    geom_smooth(method = "lm") +
    facet_grid(~size)+
    xlab("Pathotroph OTU richness") + ylab("Mass remaining after 37 months (%)") +
    theme_bw()
  ggsave(filename = "output/figures/supplementary/richness_patho_pmr.pdf", plot = p, width = 5, height = 3.5)
  
}

doAnalysis_endoDiv_explainDecayResids <- function(taxAndFunguild, comm.otu, traitResiduals.code, code.respVars){
  
  richList <- summarize_endoDiv(taxAndFunguild, comm.otu)
  
  #merge the diversity and trait residuals data
  rich.traitresid.list<-lapply(richList, function(x) {CreateRichTraitResid_df(otutype.df = x, trait.residuals = traitResiduals.code)})
  
  #fit models
  rhs <- "size * mean"
  lhs <- code.respVars
  
  richType.list<-list()
  for(i in 1:length(rich.traitresid.list)){
    resp.list<-list()
    for(k in 1:length(lhs)){
      df.sub<-subset(rich.traitresid.list[[i]], resp == lhs[[k]]) #subset the dataset based on the current resp var
      resp.list[[k]]<-ModelFit_manyYs(y = "resid", rhs = rhs, curr.data = df.sub) #run models
    }
    names(resp.list)<-lhs
    richType.list[[i]]<-resp.list
  }
  names(richType.list)<-names(richList)
  
  #summarize
  pretty.list<-lapply(richType.list, MakeLmSummaryTable, respvars=lhs)
  pretty.df <- list_to_df(pretty.list)
  
  # #identify significant relationships
  # pretty.df %>%
  #   filter(term %in% c('mean','sizesmall','sizesmall:mean')) %>%
  #   gather(key = "respvar", value = "coef", -c(source, term)) %>%
  #   mutate(signif = ifelse(grepl("*", coef, fixed = T), T, F)) %>%
  #   filter(signif == T)
  
  result.list <- list(rich.traitresid.list = rich.traitresid.list,
                      pretty.df = pretty.df)
  
  return(result.list)
  
}

makefig__saproRichDecayResids <- function(taxAndFunguild, comm.otu, traitResiduals.code, code.respVars){
  
  result.list <- doAnalysis_endoDiv_explainDecayResids(taxAndFunguild, comm.otu, traitResiduals.code, code.respVars)
  
  tmp <- result.list$rich.traitresid.list[['Sapro.rich']]
  tmp %>%
    filter(resp %in% c("beta","t60")) -> tmp
  p <- ggplot(tmp, aes(x = mean, y = resid, color = size)) +
    geom_point() +
    geom_smooth(method = "lm") +
    xlab("Saprotroph richness") + ylab("Residuals") +
    facet_wrap(~resp, scales = "free") +
    theme_bw()
  p
  ggsave(filename = "output/figures/supplementary/richness_sapro_decayResids.pdf", plot = p, width = 5, height = 3.5)
  
}

doAnalysis_endoDiv_explainPMRResids <- function(taxAndFunguild, comm.otu, traitResiduals.stem, stem.respVars){
  
  richList <- summarize_endoDiv(taxAndFunguild, comm.otu)
  
  #merge the diversity and trait residuals data
  rich.traitresid.list<-lapply(richList, function(x) {traitResiduals.stem %>% left_join(x)})
  
  #fit models
  rhs <- "size * sub_rich"
  lhs <- stem.respVars
  
  richType.list<-list()
  for(i in 1:length(rich.traitresid.list)){
    resp.list<-list()
    for(k in 1:length(lhs)){
      df.sub<-subset(rich.traitresid.list[[i]], resp == lhs[[k]]) #subset the dataset based on the current resp var
      resp.list[[k]]<-ModelFit_manyYs(y = "resid", rhs = rhs, curr.data = df.sub) #run models
    }
    names(resp.list)<-lhs
    richType.list[[i]]<-resp.list
  }
  names(richType.list)<-names(richList)
  
  #summarize
  pretty.list<-lapply(richType.list, MakeLmSummaryTable, respvars=lhs)
  pretty.df <- list_to_df(pretty.list)
  
  result.list.stem <- list(rich.traitresid.list = rich.traitresid.list,
                           pretty.df = pretty.df)
  
  return(result.list.stem)
  
}

makefig__richPMRResids <- function(taxAndFunguild, comm.otu, traitResiduals.stem, stem.respVars){
  
  result.list.stem <- doAnalysis_endoDiv_explainPMRResids(taxAndFunguild, comm.otu, traitResiduals.stem, stem.respVars)
  
  tmp.rich <- result.list.stem$rich.traitresid.list[['Richness']]
  tmp.basid <- result.list.stem$rich.traitresid.list[['Basidio.rich']]
  tmp.rich$sub_rich_type <- "All"
  tmp.basid$sub_rich_type <- "Basidiomycota"
  tmp <- rbind(tmp.rich, tmp.basid)
  tmp %>%
    select(codeStem, resp, sub_rich, sub_rich_type, size, resid) %>%
    filter(resp == "time25") -> tmp1
  
  p <- ggplot(tmp1, aes(x = sub_rich, y = resid, color = size)) +
    geom_point() +
    geom_smooth(method = "lm") +
    facet_wrap(~sub_rich_type, scales = "free") +
    xlab("OTU richness") + ylab("Residuals -  pmr time25")
  p
  ggsave(filename = "output/figures/supplementary/richness_PMRResids.pdf", plot = p, width = 5, height = 3.5)
  
}

makefig__path_richPMRResids <- function(taxAndFunguild, comm.otu, traitResiduals.stem, stem.respVars){
  
  result.list.stem <- doAnalysis_endoDiv_explainPMRResids(taxAndFunguild, comm.otu, traitResiduals.stem, stem.respVars)
  
  tmp <- result.list.stem$rich.traitresid.list[['Patho.rich']]
  tmp %>%
    select(codeStem, resp, sub_rich, size, resid) %>%
    filter(resp == "time37") -> tmp1
  
  p <- ggplot(tmp1, aes(x = sub_rich, y = resid, color = size)) +
    geom_point() +
    geom_smooth(method = "lm") +
    xlab("Pathotroph OTU richness") + ylab("Residuals - pmr time37") +
    theme_bw()
  p
  ggsave(filename = "output/figures/supplementary/richness_patho_time37Resids.pdf", plot = p, width = 5, height = 3.5)
  
}
