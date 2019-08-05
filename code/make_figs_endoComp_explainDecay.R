
plot_sampleEffortCurves<-function(mat.otu){
  
  require(vegan)
  pdf(file="output/figures/supplementary/sampleEffortCurve.pdf", width=5, height=5)
  
  rarecurve(mat.otu, step=100,
            xlab="Number of reads per sample", 
            ylab="Cumulative number of OTUs", label=FALSE)
  dev.off()
  
}

calc_codeOTUabund <- function(comm.otu, seqSamples, use.cache){
  
  comm.otu.trimmed <- removeRareOTUs(comm.otu)
  
  #average OTU abundances by code
  if(use.cache == F){
    codeOTUabund <- AverageOTUabund_byCode(comm.otu=comm.otu, seqSamples=seqSamples) #analysisDF_fxns.R... this take a while
    write.csv(codeOTUabund, file="derived_data/codeOTUabund.csv")
    codeOTUabund.trim <- AverageOTUabund_byCode(comm.otu=comm.otu.trimmed, seqSamples=seqSamples) #analysisDF_fxns.R
    write.csv(codeOTUabund.trim, file="derived_data/codeOTUabund_trim.csv")
  }else{
    codeOTUabund<-read.csv("derived_data/codeOTUabund.csv", row.names=1)
    codeOTUabund.trim<-read.csv("derived_data/codeOTUabund_trim.csv", row.names=1)
  }
  
  #remove eusc because there is no community data for this code
  #codeOTUabund[row.names(codeOTUabund) == "eusc",1:10]
  codeOTUabund <- codeOTUabund[row.names(codeOTUabund) != "eusc",]
  #codeOTUabund.trim[row.names(codeOTUabund.trim) == "eusc",1:10]
  codeOTUabund.trim <- codeOTUabund.trim[row.names(codeOTUabund.trim) != "eusc",]
  
  codeOTU.list <- list(codeOTUabund = codeOTUabund,
       codeOTUabund.trim = codeOTUabund.trim)
  return(codeOTU.list)
  
}

doAnalysis_endoComp_explainDecay <- function(comm.otu, seqSamples, decayfits, code.respVars, use.cache){
  
  require(rioja)
  
  # calculate code-level OTU abundances w/ and w/o trimming the OTU table
  codeOTU.list <- calc_codeOTUabund(comm.otu, seqSamples, use.cache)
  codeOTUabund <- codeOTU.list$codeOTUabund
  codeOTUabund.trim <- codeOTU.list$codeOTUabund.trim
  
  #order the rows of decayfits.trim so that they line up with the community mats
  sum(row.names(codeOTUabund) != row.names(codeOTUabund.trim)) #this needs to be 0
  ord<-match(row.names(codeOTUabund), decayfits$code)
  decayfits.o<-decayfits[ord,]
  sum(row.names(codeOTUabund) != decayfits.o$code) #this needs to be 0
  
  #fit models
  cvfit.notrim<-lapply(code.respVars, function(x) {fitNcrossval_WAPLS(curr.comm = codeOTUabund, curr.respVar = decayfits.o[[x]])}) #analysisDF_fxns.R
  cvfit.trim<-lapply(code.respVars, function(x) {fitNcrossval_WAPLS(curr.comm = codeOTUabund.trim, curr.respVar = decayfits.o[[x]])}) #analysisDF_fxns.R
  
  #perform randomization t-test to test the significance of a cross-validated model
  wapls.out.notrim<-lapply(cvfit.notrim, rand.t.test)
  wapls.out.trim<-lapply(cvfit.trim, rand.t.test)
  
  #make summary table
  prettyTab.notrim<-MakeSummaryTable_comcomp(wapls.out=wapls.out.notrim, respvars=unlist(code.respVars))
  prettyTab.trim<-MakeSummaryTable_comcomp(wapls.out=wapls.out.trim, respvars=unlist(code.respVars))
  prettyTab.notrim$trim<-"no"
  prettyTab.trim$trim<-"yes"
  prettyTabs.code<-rbind(prettyTab.trim, prettyTab.notrim)
  
  cvfit.results.code <- list(cvfit.notrim = cvfit.notrim,
                            cvfit.trim = cvfit.trim,
                            prettyTabs = prettyTabs.code)
  
  return(cvfit.results.code)
  
}

doAnalysis_endoComp_explainPMR <- function(comm.otu, pmr_byStem, stem.respVars){
  
  require(rioja)
  
  comm.otu.trimmed <- removeRareOTUs(comm.otu)
  
  # create dataframes
  datasets.notrim<-lapply(stem.respVars, function(x) {CreateCommPMRpair(x, comm.mat=comm.otu, pmr_byStem)} ) #analysisDF_fxns.R
  names(datasets.notrim)<-unlist(stem.respVars)
  datasets.trim<-lapply(stem.respVars, function(x) {CreateCommPMRpair(x, comm.mat=comm.otu.trimmed, pmr_byStem)} ) #analysisDF_fxns.R
  names(datasets.trim)<-unlist(stem.respVars)
  
  #fit models
  cvfit.notrim<-lapply(datasets.notrim, function(x) {fitNcrossval_WAPLS(curr.comm = x[['comm']], curr.respVar = x[['pmr']][['curr.time']])})
  cvfit.trim<-lapply(datasets.trim, function(x) {fitNcrossval_WAPLS(curr.comm = x[['comm']], curr.respVar = x[['pmr']][['curr.time']])})
  
  #perform randomization t-test to test the significance of a cross-validated model....
  wapls.out.notrim<-lapply(cvfit.notrim, rand.t.test)
  wapls.out.trim<-lapply(cvfit.trim, rand.t.test)
  
  #make summary table
  prettyTab.notrim<-MakeSummaryTable_comcomp(wapls.out=wapls.out.notrim, respvars=unlist(stem.respVars))
  prettyTab.trim<-MakeSummaryTable_comcomp(wapls.out=wapls.out.trim, respvars=unlist(stem.respVars))
  prettyTab.notrim$trim<-"no"
  prettyTab.trim$trim<-"yes"
  prettyTabs.stem<-rbind(prettyTab.trim, prettyTab.notrim)
  
  cvfit.results.stem<- list(cvfit.notrim = cvfit.notrim,
                            cvfit.trim = cvfit.trim,
                            prettyTabs = prettyTabs.stem)
  
  return(cvfit.results.stem)
  
}

maketab_endoComp_explain <- function(cvfit.results.code, cvfit.results.stem, code.respVars, stem.respVars){
  # merge code and stem-level summary tables
  tmp <- left_join(cvfit.results.code$prettyTabs, cvfit.results.stem$prettyTabs)
  commTab <- tmp[,c("trim","stat", unlist(code.respVars), unlist(stem.respVars))]
  
  return(commTab)
}

extractcoefs_wapls_score_time37 <- function(cvfit.results.stem, taxAndFunguild){
  
  #extract the significant component coefficients
  coef.comp <- cvfit.results.stem$cvfit.notrim$time37$coefficients[,'Comp01']
  coef.comp.df <- data.frame(OTUId=names(coef.comp), coefComp=coef.comp)
  
  #create df of OTU taxon and guild info matched with coef value
  coef.comp.df %>% 
    left_join(taxAndFunguild) -> coef.comp.ann
  
  return(coef.comp.ann)
  
}

makefig__wapls_score_time37 <- function(comm.otu, pmr_byStem, stem.respVars, taxAndFunguild){
  
  cvfit.results.stem <- doAnalysis_endoComp_explainPMR(comm.otu, pmr_byStem, stem.respVars)
  
  coef.comp.ann <- extractcoefs_wapls_score_time37(cvfit.results.stem, taxAndFunguild)
  
  #(1) plot distribution of WA-PLS scores for each OTU
  quant <- quantile(coef.comp.ann$coefComp, c(.01, .99))
  p1 <- ggplot(coef.comp.ann, aes(x = reorder(OTUId, coefComp), y = coefComp)) +
    geom_point(alpha = .5, pch = 16) +
    theme_bw() +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    xlab("OTU identity") + ylab("WA-PLS score (pmr at time37)") +
    geom_hline(yintercept = quant[1], linetype = 2) +
    geom_hline(yintercept = quant[2], linetype = 2) 
  p1
  #ggsave(filename = "output/figures/supplementary/wapls_score_time37.pdf", plot = p, width = 4, height = 4)
  
  
  # #(2) plot by phylo and functional groupings
  # #shorten x axis labels
  # orgNames<-levels(factor(coef.comp.ann$Trophic.Mode))
  # newNames<-gsub('troph', '', orgNames)
  # coef.comp.ann %>% mutate(Trophic.Mode.short = factor(Trophic.Mode, labels = newNames)) -> coef.comp.ann
  # 
  # p.troph<-ggplot(coef.comp.ann, aes(x=Trophic.Mode.short, y=coefComp)) + 
  #   geom_jitter(alpha=.5, pch=16) +
  #   xlab("Trophic mode") + ylab("WA-PLS score (pmr at time37)") +
  #   theme_classic() +
  #   geom_hline(yintercept=0, linetype=2) +
  #   theme(axis.text.x = element_text(angle = 70, hjust = 1))
  # 
  # p.king<-ggplot(coef.comp.ann, aes(x=kingdom, y=coefComp)) + 
  #   geom_jitter(alpha=.5) + 
  #   xlab("Kingdom") + ylab("Component coefficient estimate") +
  #   theme_classic() +
  #   geom_hline(yintercept=0, linetype=2) +
  #   theme(axis.text.y = element_blank(), 
  #         axis.ticks.y = element_blank(), 
  #         axis.title.y = element_blank(),
  #         plot.margin = unit(c(0,1,0,0), "lines"),
  #         plot.background = element_blank(),
  #         axis.text.x = element_text(angle = 70, hjust = 1))
  # 
  # p.phylum<-ggplot(coef.comp.ann, aes(x=phylum, y=coefComp)) + 
  #   geom_jitter(alpha=.5) + 
  #   xlab("Phylum") + ylab("Component coefficient estimate") +
  #   theme_classic() +
  #   geom_hline(yintercept=0, linetype=2) +
  #   theme(axis.text.y = element_blank(), 
  #         axis.ticks.y = element_blank(), 
  #         axis.title.y = element_blank(),
  #         plot.margin = unit(c(0,1,0,0), "lines"),
  #         plot.background = element_blank(),
  #         axis.text.x = element_text(angle = 70, hjust = 1))
  # 
  # #pdf("output/figures/supplementary/wapls_score_time37_otuCats.pdf", width=10, height=6)
  # grid.newpage()
  # grid.draw(cbind(ggplotGrob(p.troph), ggplotGrob(p.king), ggplotGrob(p.phylum), size = "last"))
  # #dev.off()
  
  #(3) plot by OTU's boral scores
  # size was included as a roweffect in the boral model, so there are no OTU-specific estimates
  xVar.df<-read.csv("data/xVar_OTUcoefs.csv", row.names=1) #import wood trait coefs estimated by boral in the wooddecay repo
  xVar.df %>%
    left_join(coef.comp.ann) %>% # left join is going to drop the oomycetes that weren't included in the boral analysis
    filter(runNum=="runNum1") -> tmp.df #just look at the first model run estimates
  
  tmp.df %>%
    select(OTUId, phylum, Trophic.Mode, coefComp, waterperc, P) %>%
    gather(key="trait", value="coefEst", -(1:4)) %>%
    arrange(-coefComp) -> tmp.df1
  tmp.df1%>%
    filter(trait == 'waterperc') -> tmp
  
  p1 <- ggplot(tmp, aes(x = as.numeric(reorder(OTUId, coefComp)), y = coefComp)) +
    geom_point() +
    xlab("OTU identity") + ylab("WA-PLS score (pmr at time37)") +
    theme(axis.text.x=element_blank(),
                       axis.ticks.x=element_blank())
  p1 
  
  p2<-ggplot(tmp, aes(y=coefComp, x=coefEst)) + 
    geom_point() + 
    ylab("") + 
    xlab("Assoc. w/ water content (boral coef)")
  p2
  #ggsave(filename = "output/figures/supplementary/wapls_score_time37_boral.pdf", plot = p, width = 6, height = 4)
  
  #waterperc
  tmp.df1 %>%
    filter(trait == "waterperc") -> tmp
  mod <- lm(coefEst ~ coefComp, data = tmp)
  summary(mod)
  
  pdf("output/figures/supplementary/wapls_score_time37.pdf", width=8, height=4)
  grid.newpage()
  grid.draw(cbind(ggplotGrob(p1), ggplotGrob(p2), size = "last"))
  dev.off()

  
  
}

doAnalysis_endoComp_explainDecayResids <- function(comm.otu, seqSamples, decayfits, code.respVars, traitResiduals.code, use.cache){
  
  require(rioja)
  
  # calculate code-level OTU abundances w/ and w/o trimming the OTU table
  codeOTU.list <- calc_codeOTUabund(comm.otu, seqSamples, use.cache)
  codeOTUabund <- codeOTU.list$codeOTUabund
  codeOTUabund.trim <- codeOTU.list$codeOTUabund.trim
  
  # create dataframes
  datasets.notrim<-lapply(code.respVars, function(x) {
    CreateCommTraitResidpair(respVar=x, 
                             comm.mat=codeOTUabund, 
                             traitResiduals = traitResiduals.code,
                             sampleName = "code")} ) #analysisDF_fxns.R
  names(datasets.notrim) <- unlist(code.respVars)
  datasets.trim<-lapply(code.respVars, function(x) {
    CreateCommTraitResidpair(respVar=x, 
                             comm.mat=codeOTUabund.trim, 
                             traitResiduals = traitResiduals.code,
                             sampleName = "code")} ) #analysisDF_fxns.R
  names(datasets.trim)<-unlist(code.respVars)
  
  #fit models
  cvfit.notrim<-lapply(datasets.notrim, function(x) {
    fitNcrossval_WAPLS(curr.comm = x[['comm']], 
                       curr.respVar = x[['traitresid']][['resid']])
  })
  cvfit.trim<-lapply(datasets.trim, function(x) {
    fitNcrossval_WAPLS(curr.comm = x[['comm']], 
                       curr.respVar = x[['traitresid']][['resid']])
  })
  
  #perform randomization t-test to test the significance of a cross-validated model....
  wapls.out.notrim<-lapply(cvfit.notrim, rand.t.test)
  wapls.out.trim<-lapply(cvfit.trim, rand.t.test)
  
  #make summary table
  prettyTab.notrim<-MakeSummaryTable_comcomp(wapls.out=wapls.out.notrim, respvars=unlist(code.respVars))
  prettyTab.trim<-MakeSummaryTable_comcomp(wapls.out=wapls.out.trim, respvars=unlist(code.respVars))
  prettyTab.notrim$trim<-"no"
  prettyTab.trim$trim<-"yes"
  prettyTabs.code<-rbind(prettyTab.trim, prettyTab.notrim)
  
  cvfit.results.code <- list(cvfit.notrim = cvfit.notrim,
                             cvfit.trim = cvfit.trim,
                             prettyTabs = prettyTabs.code)
  
  return(cvfit.results.code)

}

doAnalysis_endoComp_explainPMRResids <- function(comm.otu, pmr_byStem, stem.respVars, traitResiduals.stem){
  
  require(rioja)
  
  # create dataframes
  datasets.notrim<-lapply(stem.respVars, function(x) {
    CreateCommTraitResidpair(respVar=x, 
                             comm.mat=comm.otu, 
                             traitResiduals = traitResiduals.stem,
                             sampleName = "codeStem")} ) #analysisDF_fxns.R
  names(datasets.notrim)<-unlist(stem.respVars)
  
  comm.otu.trimmed <- removeRareOTUs(comm.otu)
  datasets.trim<-lapply(stem.respVars, function(x) {
    CreateCommTraitResidpair(respVar=x, 
                             comm.mat=comm.otu.trimmed, 
                             traitResiduals = traitResiduals.stem,
                             sampleName = "codeStem")} ) #analysisDF_fxns.R
  names(datasets.trim)<-unlist(stem.respVars)
  
  #fit models
  cvfit.notrim<-lapply(datasets.notrim, function(x) {
    fitNcrossval_WAPLS(curr.comm = x[['comm']], 
                       curr.respVar = x[['traitresid']][['resid']])
  })
  cvfit.trim<-lapply(datasets.trim, function(x) {
    fitNcrossval_WAPLS(curr.comm = x[['comm']], 
                       curr.respVar = x[['traitresid']][['resid']])
  })
  
  #perform randomization t-test to test the significance of a cross-validated model....
  wapls.out.notrim<-lapply(cvfit.notrim, rand.t.test)
  wapls.out.trim<-lapply(cvfit.trim, rand.t.test)
  
  #make summary table
  prettyTab.notrim<-MakeSummaryTable_comcomp(wapls.out=wapls.out.notrim, respvars=unlist(stem.respVars))
  prettyTab.trim<-MakeSummaryTable_comcomp(wapls.out=wapls.out.trim, respvars=unlist(stem.respVars))
  prettyTab.notrim$trim<-"no"
  prettyTab.trim$trim<-"yes"
  prettyTabs.stem<-rbind(prettyTab.trim, prettyTab.notrim)
  
  cvfit.results.stem<- list(cvfit.notrim = cvfit.notrim,
                            cvfit.trim = cvfit.trim,
                            prettyTabs = prettyTabs.stem)
  
  return(cvfit.results.stem)
  
}

doAnalysis_endoComp_woodTraits <- function(comm.otu, seqSamples, traits.code, use.cache){
  
  # calculate code-level OTU abundances w/ and w/o trimming the OTU table
  codeOTU.list <- calc_codeOTUabund(comm.otu, seqSamples, use.cache)
  codeOTUabund <- codeOTU.list$codeOTUabund
  codeOTUabund.trim <- codeOTU.list$codeOTUabund.trim
  
  # remove trait NAs
  traits.code %>% 
    filter(!is.na(waterperc) & !is.na(P)) -> traits.noNAs
  
  # make pair of matching datasets with the community mat and the trait mat
  datasets.notrim.code<-CreateCommTraitpair(comm.otu = codeOTUabund, traits=traits.noNAs, sampleName = "code")
  datasets.trim.code<-CreateCommTraitpair(comm.otu = codeOTUabund.trim, traits=traits.noNAs, sampleName = "code")
  
  # choose a model by permutation tests in a constrained ordination
  if(use.cache == F){
    mod.nt.code<-ordistep_wrapper(datasets=datasets.notrim.code) #this can take a while
    saveRDS(mod.nt.code, file = "derived_data/modSelect_nt_code.RData")
    mod.t.code<-ordistep_wrapper(datasets=datasets.trim.code)
    saveRDS(mod.t.code, file = "derived_data/modSelect_t_code.RData")
  }else{
    mod.nt.code <- readRDS(file = "derived_data/modSelect_nt_code.RData")
    mod.t.code <- readRDS(file = "derived_data/modSelect_t_code.RData")
  }
  
  #variance inflation factor (VIF) which is 1 for completely independent variables, and values above 10 or 20 (depending on your taste) are regarded as highly multicollinear (dependent on others)
  #vif.cca(mod.nt.code) # none are highly multicollinear
  #vif.cca(mod.t.code)
  
  mod.list <- list(mod.nt.code = mod.nt.code,
                   mod.t.code = mod.t.code)
  
  return(mod.list)
  
}

makefigs__endoComp_woodTraits <- function(comm.otu, seqSamples, traits.code, traits.stem, use.cache){
  
  #(1) by code
  mod.list<- doAnalysis_endoComp_woodTraits(comm.otu, seqSamples, traits.code, use.cache)
  mod.nt.code <- mod.list$mod.nt.code
  mod.t.code <- mod.list$mod.t.code
  # proportion of constrained variance (inertia)
  prop.constr.nt.code<-paste("prop. constr. =", round(extract_constrainedInertia_proport(mod.nt.code), digits=2))
  prop.constr.t.code<-paste("prop. constr. =", round(extract_constrainedInertia_proport(mod.t.code), digits=2))
  
  
  #(2) by stem
  mod.list <- doAnalysis_endoComp_woodTraits.stem(comm.otu, traits.code, traits.stem, use.cache)
  mod.nt.stem <- mod.list$mod.nt.stem
  mod.t.stem <- mod.list$mod.t.stem
  # proportion of constrained variance (inertia)
  prop.constr.nt.stem<-paste("prop. constr. =", round(extract_constrainedInertia_proport(mod.nt.stem), digits=2))
  prop.constr.t.stem<-paste("prop. constr. =", round(extract_constrainedInertia_proport(mod.t.stem), digits=2))
  
  pdf("output/figures/supplementary/dbRDA_traits.pdf", width=12, height=6)
  par(mfrow=c(1,2))
  
  # plot(mod.nt.code, display = c("wa","bp")) # constrained with best model
  # mtext(prop.constr.nt.code, side=3, adj=0.9, line=-1.5, col=4)
  # title('Full community')
  
  plot(mod.t.code, display = c("wa","bp")) # constrained with best model
  mtext(prop.constr.t.code, side=3, adj=0.9, line=-1.5, col=4)
  title('By code')
  
  # plot(mod.nt.stem, display = c("wa","bp")) # constrained with best model
  # mtext(prop.constr.nt.stem, side=3, adj=0.9, line=-1.5, col=4)
  # title('Full community')
  
  plot(mod.t.stem, display = c("wa","bp")) # constrained with best model
  mtext(prop.constr.t.stem, side=3, adj=0.9, line=-1.5, col=4)
  title('By stem')
  
  dev.off()
  
  
}

maketabs__endoComp_woodTraits <- function(comm.otu, seqSamples, traits.code, use.cache){
  
  mod.list<- doAnalysis_endoComp_woodTraits(comm.otu, seqSamples, traits.code, use.cache)
  mod.nt.code <- mod.list$mod.nt.code
  mod.t.code <- mod.list$mod.t.code
  
  #save summary tables
  an.nt.code <- anova.margin.table(dbrda.obj=mod.nt.code)
  an.t.code<-anova.margin.table(dbrda.obj=mod.t.code)
  tab.list <- list(an.nt.code = an.nt.code,
                   an.t.code = an.t.code)
  
  return(tab.list)
}

doAnalysis_endoComp_woodTraits.stem <- function(comm.otu, traits.code, traits.stem, use.cache){
  
  comm.otu.trimmed <- removeRareOTUs(comm.otu)
  
  # put together best stem-level trait matrix, remove NAs
  traits.code <- data.frame(traits.code)
  traits.code %>%
    select(code, barkthick, density) %>%
    rename('barkthick_smspp'='barkthick',
           'density_smspp'='density')-> select.traits.code
  traits.stem %>%
    left_join(select.traits.code) %>%
    select(-density) %>%
    select(-barkthick) %>%
    filter(!is.na(waterperc) & !is.na(P)) -> traits.stem.plus
  
  # make pair of matching datasets with the community mat and the trait mat
  datasets.notrim.stem<-CreateCommTraitpair(comm.otu = comm.otu, traits = traits.stem.plus, sampleName = "codeStem")
  datasets.trim.stem<-CreateCommTraitpair(comm.otu = comm.otu.trimmed, traits = traits.stem.plus, sampleName = "codeStem")
  
  # choose a model by permutation tests in a constrained ordination
  if(use.cache == F){
    mod.nt.stem<-ordistep_wrapper(datasets=datasets.notrim.stem) #this can take a while
    saveRDS(mod.nt.stem, file = "derived_data/modSelect_nt_stem.RData")
    mod.t.stem<-ordistep_wrapper(datasets=datasets.trim.stem)
    saveRDS(mod.t.stem, file = "derived_data/modSelect_t_stem.RData")
  }else{
    mod.nt.stem <- readRDS(file = "derived_data/modSelect_nt_stem.RData")
    mod.t.stem <- readRDS(file = "derived_data/modSelect_t_stem.RData")
  }
  
  #variance inflation factor (VIF) which is 1 for completely independent variables, and values above 10 or 20 (depending on your taste) are regarded as highly multicollinear (dependent on others)
  #vif.cca(mod.nt.stem) # none are highly multicollinear
  #vif.cca(mod.t.stem)
  
  mod.list <- list(mod.nt.stem = mod.nt.stem,
                   mod.t.stem = mod.t.stem)
  
  return(mod.list)
  
}

maketabs__endoComp_woodTraits.stem <- function(comm.otu, traits.code, traits.stem, use.cache){
  
  mod.list <- doAnalysis_endoComp_woodTraits.stem(comm.otu, traits.code, traits.stem, use.cache)
  mod.nt.stem <- mod.list$mod.nt.stem
  mod.t.stem <- mod.list$mod.t.stem
  
  #save summary tables
  an.nt.stem <- anova.margin.table(dbrda.obj=mod.nt.stem)
  an.t.stem<-anova.margin.table(dbrda.obj=mod.t.stem)
  tab.list <- list(an.nt.stem = an.nt.stem,
                   an.t.stem = an.t.stem)
  
  return(tab.list)
}

