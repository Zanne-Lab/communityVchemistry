
# load wood phylogentic data

fix_problem_species<-function(tree, prob_species, dontchoose = wood_names){
  
  require(taxonlookup)
  
  for(i in 1:length(prob_species)){
    genus<-taxonlookup:::split_genus(prob_species[i])
    replace<-sample(tree$tip.label[grepl(genus,tree$tip.label)],1)
    while (replace%in%dontchoose) replace<-sample(tree$tip.label[grepl(genus,tree$tip.label)],1)
    tree$tip.label[tree$tip.label==replace]<-prob_species[i]
  }
  return(tree)
}

load_zanne_tree<-function(){
  
  require(phytools)
  require(diversitree)
  
  zae <- read.tree("data/zanne1.1.tre")
  stemSamples <- load_stemSamples()
  wood_names<-sub(" ", "_", unique(stemSamples$Binomial))
  set.seed(42)
  zae_mod<-fix_problem_species(zae, 
                               prob_species = c("Melaleuca_decora","Petrophile_pulchella","Persoonia_nutans","Callistemon_linearis","Eucalyptus_sclerophylla","Isopogon_anemonifolius","Hakea_sericea","Olax_stricta","Jacksonia_scoparia"), 
                               dontchoose = wood_names)
  zanneTree<-diversitree:::drop.tip.fixed(phy = zae_mod,zae_mod$tip.label[!zae_mod$tip.label%in%wood_names])
  
  return(zanneTree)
  
}

# raw trait datasets

load_waterPercent.perGwetmass<-function(){
  
  require(dplyr)
  
  # read in initial covariate data
  covar.big <-read.csv('data/covariates_bigStems.csv', stringsAsFactor = F)
  covar.small <-read.csv('data/covariates_smallStems.csv', stringsAsFactor = F)
  
  # clean the column names
  covar.big %>%
    rename('Diameter.cm'='Diameter..cm.',
           'Fresh.mass.g'='Fresh.mass..g.',
           'Dry.mass.g'='Dry.mass..g.') -> covar.big
  
  covar.small %>%
    rename('Diameter.wbark.mm'='Diameter.wbark..mm.',
           'Fresh.mass.g'='Fresh.mass..g.',
           'Diameter.nobark.mm'='Diameter.nobark..mm.',
           'Volume.g'='Volume..g.',
           'Dry.mass.wood.g'='Dry.mass.wood..g.',
           'Dry.mass.bark.g'='Dry.mass.bark..g.',
           'Dry.mass.total.g'='Dry.mass.total..g.') -> covar.small
  
  # (1) calculate water content for each species, size class
  
  covar.big %>%
    mutate(waterperc = ((Fresh.mass.g - Dry.mass.g) / Fresh.mass.g) * 100) %>%
    filter(!is.na(waterperc)) %>%
    select(Species, Stem, unique, waterperc) %>%
    mutate(size = "large") -> big.waterperc
  
  covar.small %>%
    mutate(Dry.mass.total.g = Dry.mass.wood.g + Dry.mass.bark.g) %>% # add Dry.mass.wood.g and Dry.mass.bark.g to get Dry.mass.total.g
    mutate(waterperc = ((Fresh.mass.g - Dry.mass.total.g) / Fresh.mass.g) * 100) %>%
    filter(!is.na(waterperc)) %>%
    select(Species, Stem, unique, waterperc) %>%
    mutate(size = "small") -> small.waterperc
  
  water.perc <- rbind(big.waterperc, small.waterperc)
  water.perc %>%
    rename('code'='Species') %>%
    select(unique, code, size, Stem, waterperc) -> water.perc
  
  #unique(water.perc$Stem) #all numeric
  #length(unique(water.perc$code)) #34
  
  return(water.perc)
}

load_densityNbarkthick<-function(){
  
  require(dplyr)
  
  # read in initial covariate data
  covar.small <-read.csv('data/covariates_smallStems.csv', stringsAsFactor = F)
  # only measured on small diameter samples
  
  # clean the column names
  covar.small %>%
    rename('Diameter.wbark.mm'='Diameter.wbark..mm.',
           'Fresh.mass.g'='Fresh.mass..g.',
           'Diameter.nobark.mm'='Diameter.nobark..mm.',
           'Volume.g'='Volume..g.',
           'Dry.mass.wood.g'='Dry.mass.wood..g.',
           'Dry.mass.bark.g'='Dry.mass.bark..g.',
           'Dry.mass.total.g'='Dry.mass.total..g.') -> covar.small
  
  # calculate wood density and bark thickness for each species
  covar.small %>%
    mutate(density = Dry.mass.wood.g / Volume.g) %>%
    mutate(barkthick = Diameter.wbark.mm - Diameter.nobark.mm) %>%
    filter(!is.na(density) | !is.na(barkthick))%>%
    rename('code'='Species') %>%
    mutate(size = ifelse(code == tolower(code), "small","large")) %>%
    select(unique, code, size, Stem, density, barkthick) -> densityBarkthick
  
  #unique(densityBarkthick$Stem) # all numeric
  #length(unique(densityBarkthick$code)) #22
  
  return(densityBarkthick)
  
}

load_XRF_meta <- function(){
  
  data <- read.csv('data/sequencing_T0/NextGenSeqencing_2016_Sample_MasterSpreadsheet.csv', stringsAsFactors=F)
  data <- data[!data$SampleCode == 'blank', ] # extract 'blank' and 'mock' samples from 'mat.otu', delete from 'data'
  
  # create dataframe containing metadata for initial sequencing and XRF data
  #if there is a number or letter on the back of the SampleCode, then it was composited
  data %>%
    select(SampleCode, StemSize, mgSample, NucleicAcidConc, ExtractionDate) %>%
    separate(SampleCode, into = c("code","Stem_maybe"), sep = 4, remove = FALSE) %>%
    mutate(compositeSample = ifelse(grepl('[0-9A-Z]$', SampleCode), FALSE, TRUE)) -> meta #stem assignment are only numbers
  meta %>%
    select(SampleCode, code, StemSize, Stem_maybe, compositeSample) -> meta.indx
  
  return(meta.indx)
}

load_XRF <- function(){
  
  require(dplyr)
  
  data <- read.csv('data/sequencing_T0/NextGenSeqencing_2016_Sample_MasterSpreadsheet.csv', stringsAsFactors=F)
  data <- data[!data$SampleCode == 'blank', ] # extract 'blank' and 'mock' samples from 'mat.otu', delete from 'data'
  
  # create dataframe containing metadata for initial sequencing and XRF data
  #if there is a number or letter on the back of the SampleCode, then it was composited
  meta.indx <- load_XRF_meta()
  
  #isolate XRF cols
  data %>%
    select(SampleCode, StemSize, P, K, Ca, Mn, Fe, Zn) %>%
    left_join(meta.indx) %>%
    rename('unique'='SampleCode',
           'size'='StemSize') %>%
    mutate(Stem = ifelse(grepl('[0-9]$', Stem_maybe), Stem_maybe, NA)) %>%
    filter(!is.na(P)) %>%
    select(unique, code, size, Stem, compositeSample, P, K, Ca, Mn, Fe, Zn) -> xrf
  
  # length(unique(xrf$code)) #33
  
  return(xrf)
}

load_CN<-function(){
  
  require(dplyr)
  
  # read in CN data
  cndata <- read.csv('data/CN/JEFF_POWELL_CN_DATA_DIVERSITY_ROT_AUG_2013.csv', stringsAsFactors=F)
  colnames(cndata)<-c("SampleID","comments","mass","N","C")
  
  #create meta from XRF meta data
  meta.indx <- load_XRF_meta()
  
  # add meta to cndata
  cndata %>%
    separate(SampleID, into = c("code","Stem_maybe"), sep = 4, remove = F) %>%
    mutate(Stem_maybe = ifelse(grepl("pooled", Stem_maybe), NA, Stem_maybe)) %>%
    mutate(code = ifelse(code == "Cali", "cali", code)) %>% #this species is only represented by small diams
    mutate(code = ifelse(code == "Olst", "olst", code)) %>% #this species is only represented by small diams
    mutate(code = ifelse(code == "Basp", "basp", code)) %>% #this species is only represented by small diams
    mutate(SampleCode = ifelse(is.na(Stem_maybe), code, paste0(code, Stem_maybe))) %>%
    select(SampleID, code, SampleCode, N, C) %>%
    left_join(meta.indx) %>%
    rename('size'='StemSize',
           'unique'='SampleCode') %>%
    mutate(Stem = ifelse(grepl('[0-9]$', Stem_maybe), Stem_maybe, NA)) %>%
    select(unique, code, size, Stem, N, C, compositeSample) -> cndata.clean
  
  #length(unique(cndata.clean$code)) #33
  #all.codes <-unique(stemSamples$code)
  #all.codes[!all.codes %in% unique(cndata.clean$code)] #missing eusc
  
  return(cndata.clean)
  
}

load_Cfract<-function(){
  
  cfract <- read_csv(file = "data/Marissa_Chemistry_Data_export.csv")
  cfract <- cfract[complete.cases(cfract),] # get rid of empty row from the excel file
  
  #I think there is a typo and JASC2 was entered twice instead of JASC3. My original datasheet that I sent included JASC3. 
  #It would be good to double check this with Shawn's lab.
  
  # fix typo
  cfract %>%
    mutate(codeStem = ifelse(codeStem == "JASC2" &  labID %in% c(47,48), 
                             "JASC3", 
                             codeStem)) -> cfract
  
  #Look for analytical reps that are way different
  num.out <- length(unique(cfract$labID))/2
  cfract %>%
    mutate(rep = paste0("rep",rep(c(1,2), times = num.out))) %>%
    select(codeStem, rep, perc.Total.lignin, Arabinose, Rhamnose, Galactose, Glucose, Xylose, Mannose, Total) %>%
    gather(key = "type", value = "value", -c(codeStem, rep)) %>%
    spread(key = rep, value = value) %>%
    mutate(diff = abs(rep1 - rep2)) -> tmp
  #plot difference between analytical rep1 and rep2
  ggplot(tmp, aes(x = codeStem, y = diff)) +
    geom_point() +
    facet_wrap(~type, scales = "free_x") +
    theme_bw() +
    coord_flip()
  
  #Go back and investigate the JASC/Arabinose and ACPA/Xylose outliers. 
  #For now, just average the analytical reps
  
  # average the 2 analytical reps
  num.out <- length(unique(cfract$labID))/2
  cfract %>%
    mutate(rep = paste0("rep",rep(c(1,2), times = num.out))) %>%
    select(codeStem, rep, perc.Total.lignin, Arabinose, Rhamnose, Galactose, Glucose, Xylose, Mannose, Total) %>%
    gather(key = "type", value = "value", -c(codeStem, rep)) %>%
    spread(key = rep, value = value) %>%
    mutate(diff = abs(rep1 - rep2)) %>%
    mutate(avg = (rep1 + rep2) / 2) %>%
    select(codeStem, type, avg) %>%
    spread(key = type, value = avg) -> cfract.w
  
  # add unique, code, size, Stem
  cfract.w %>%
    separate(codeStem, into = c("code","Stem"), sep = 4, remove = F) %>%
    mutate(size = ifelse(code == tolower(code), "small","large")) %>%
    rename('unique'=codeStem) -> cfract.wt
  
  return(cfract.wt)
}

# adds C fraction data
mergeTraitData<-function(){
  
  #load raw trait datasets
  waterperc <- load_waterPercent.perGwetmass() ##### this is in units of g water per g of wet mass x 100
  densityNbarkthick <- load_densityNbarkthick()
  xrf <- load_XRF()
  cn <- load_CN()
  cfract <-load_Cfract()
  
  #make long
  waterperc %>%
    mutate(trait = "waterperc") %>%
    rename('trait.val'='waterperc') %>%
    mutate(compositeSample = FALSE) %>%
    select(unique, code, size, Stem, compositeSample, trait, trait.val) -> waterperc.l
  densityNbarkthick %>%
    gather(key = "trait", value = "trait.val", c(density, barkthick)) %>%
    mutate(compositeSample = FALSE) %>%
    select(unique, code, size, Stem, compositeSample, trait, trait.val) -> densityNbarkthick.l
  xrf %>%
    gather(key = "trait", value = "trait.val", c(P, K, Ca, Mn, Fe, Zn)) %>%
    select(unique, code, size, Stem, compositeSample, trait, trait.val) -> xrf.l
  cn %>%
    gather(key = "trait", value = "trait.val", c(C, N)) %>%
    select(unique, code, size, Stem, compositeSample, trait, trait.val) -> cn.l
  cfract %>%
    select(-Total) %>%
    rename('Lignin'='perc.Total.lignin') %>%
    gather(key = "trait", value = "trait.val", -c(unique, code, Stem, size)) %>%
    mutate(compositeSample = FALSE) %>%
    select(unique, code, size, Stem, compositeSample, trait, trait.val) -> cfract.l
  
  trait.data.l <- rbind(waterperc.l, densityNbarkthick.l, xrf.l, cn.l, cfract.l)
  
  return(trait.data.l)
  
}

traitcol.order <- function(){
  id.cols <- c("code","species","size")
  phys.cols <- c("waterperc","density","barkthick")
  nutr.cols <- c("C","N","P","K","Ca","Fe","Mn","Zn")
  cfract.cols <- c("Arabinose","Galactose","Glucose", "Mannose", "Rhamnose","Xylose","Lignin")
  stem.trait.cols <- c("waterperc","density_smspp","barkthick_smspp", nutr.cols)
  trait.order <- list(id.cols = id.cols, 
                      phys.cols = phys.cols, 
                      nutr.cols = nutr.cols, 
                      cfract.cols = cfract.cols,
                      stem.trait.cols = stem.trait.cols)
  
  return(trait.order)
}

# code-level traits

#if fill.densitybark == TRUE, then use small stem estimates to approximate large stem density and bark thickness
trait.means_byCode <- function(stemSamples, fill.densitybark){
  
  #load trait data
  trait.data.l <- mergeTraitData()
  
  #summarize
  trait.data.l %>%
    group_by(code, trait) %>%
    summarize(val = mean(trait.val, na.rm=T)) %>%
    spread(key = trait, value = val) -> traitmeans.code
  
  #add species and size
  samp.indx <- unique(stemSamples[,c("code","species","size")])
  traitmeans.code %>%
    left_join(samp.indx) -> traitmeans.code
  
  #reorder the columns
  order <- traitcol.order()
  trait.order <- c(order$id.cols, order$phys.cols, order$nutr.cols, order$cfract.cols)
  traitmeans.code %>%
    select(trait.order) -> traitmeans.code
  
  if(fill.densitybark == TRUE){
    # use species-level small-stem estimates of density and barkthick for large-stem samples
    
    #pull out the small-stem estimates
    traitmeans.code %>%
      filter(size == "small") %>%
      select(code, density, barkthick) %>%
      mutate(species = code) %>%
      mutate(size = "large") -> filler.data
    tmp <- data.frame(filler.data)
    tmp %>%
      select(-code) %>%
      mutate(code = toupper(species)) -> filler.data
    
    #identify which species overlap in the large size class
    fillcodes <- unique(filler.data$code)
    allcodes <- unique(traitmeans.code$code)
    
    #loop through the large samples that can be filled in
    CODE <- allcodes[allcodes %in% fillcodes]
    for(i in 1:length(CODE)){
      fill.vals <- filler.data[filler.data$code == CODE[i], c("density","barkthick")]
      traitmeans.code[traitmeans.code$code == CODE[i], c("density","barkthick")] <- fill.vals
    }
    
  }
  return(traitmeans.code)
  
}

trait.sds_byCode <- function(stemSamples){
  
  #load trait data
  trait.data.l <- mergeTraitData()
  
  #summarize
  trait.data.l %>%
    group_by(code, trait) %>%
    summarize(val = sd(trait.val, na.rm=T)) %>%
    spread(key = trait, value = val) -> traitsds.code
  
  #add species and size
  samp.indx <- unique(stemSamples[,c("code","species","size")])
  traitsds.code %>%
    left_join(samp.indx) -> traitsds.code

  #reorder the columns
  order <- traitcol.order()
  trait.order <- c(order$id.cols, order$phys.cols, order$nutr.cols, order$cfract.cols)
  traitsds.code %>%
    select(trait.order) -> traitsds.code
  
  return(traitsds.code)
  
}

trait.n_byCode <- function(stemSamples){
  
  #load trait data
  trait.data.l <- mergeTraitData()
  
  #summarize
  trait.data.l %>%
    group_by(code, trait) %>%
    summarize(val = sum(!is.na(trait.val))) %>%
    spread(key = trait, value = val) -> traitn.code
  
  #add species and size
  samp.indx <- unique(stemSamples[,c("code","species","size")])
  traitn.code %>%
    left_join(samp.indx) -> traitn.code
  
  #reorder the columns
  order <- traitcol.order()
  trait.order <- c(order$id.cols, order$phys.cols, order$nutr.cols, order$cfract.cols)
  traitn.code %>%
    select(trait.order) -> traitn.code
  
  return(traitn.code)
  
}

# stem-level traits

trait.means_byStem <- function(stemSamples){
  
  #load trait data
  trait.data.l <- mergeTraitData()
  
  #summarize
  trait.data.l %>%
    filter(!is.na(Stem)) %>%
    mutate(codeStem = paste0(code, Stem)) %>%
    group_by(codeStem, trait) %>%
    summarize(val = mean(trait.val, na.rm=T)) %>%
    spread(key = trait, value = val) -> traitmeans.stem
  
  #add species and size
  samp.indx <- unique(stemSamples[,c("codeStem","species","size","code")])
  traitmeans.stem %>%
    left_join(samp.indx) -> traitmeans.stem
  
  #reorder the columns
  order <- traitcol.order()
  trait.order <- c(order$id.cols, order$phys.cols, order$nutr.cols, order$cfract.cols)
  if(sum(!trait.order %in% colnames(traitmeans.stem)) == 0){
    traitmeans.stem %>%
      select(c("codeStem", trait.order)) -> traitmeans.stem
  }else{
    print("warning: trait columns are out of order")
  }

  #note
  #there are 2 stem-level samples that are in the traits df and endophytes df but we're in the deployment df
  # acpa2 and lepa4
  # manually fill in the species and size column for these samples
  tmp <- traitmeans.stem
  tmp[tmp$codeStem == "acpa2", c("code","species","size")] <- c("acpa","acpa","small")
  tmp[tmp$codeStem == "lepa4", c("code","species","size")] <- c("lepa","lepa","small")
  traitmeans.stem <- tmp
  
  return(traitmeans.stem)
  
}

trait.sds_byStem <- function(stemSamples){
  
  #load trait data
  trait.data.l <- mergeTraitData()
  
  #summarize
  trait.data.l %>%
    filter(!is.na(Stem)) %>%
    mutate(codeStem = paste0(code, Stem)) %>%
    group_by(codeStem, trait) %>%
    summarize(val = sd(trait.val, na.rm=T)) %>%
    spread(key = trait, value = val) -> traitsds.stem
  
  #add species and size
  samp.indx <- unique(stemSamples[,c("codeStem","species","size","code")])
  traitsds.stem %>%
    left_join(samp.indx) -> traitsds.stem
  
  #reorder the columns
  order <- traitcol.order()
  trait.order <- c(order$id.cols, order$phys.cols, order$nutr.cols, order$cfract.cols)
  if(sum(!trait.order %in% colnames(traitsds.stem)) == 0){
    traitsds.stem %>%
      select(c("codeStem", trait.order)) -> traitsds.stem
  }else{
    print("warning: trait columns are out of order")
  }
  
  #note
  #there are 2 stem-level samples that are in the traits df and endophytes df but we're in the deployment df
  # acpa2 and lepa4
  # manually fill in the species and size column for these samples
  tmp <- traitsds.stem
  tmp[tmp$codeStem == "acpa2", c("code","species","size")] <- c("acpa","acpa","small")
  tmp[tmp$codeStem == "lepa4", c("code","species","size")] <- c("lepa","lepa","small")
  traitsds.stem <- tmp
  
  return(traitsds.stem)
  
}

trait.n_byStem <- function(stemSamples){
  
  #load trait data
  trait.data.l <- mergeTraitData()
  
  #summarize
  trait.data.l %>%
    filter(!is.na(Stem)) %>%
    mutate(codeStem = paste0(code, Stem)) %>%
    group_by(codeStem, trait) %>%
    summarize(val = sum(!is.na(trait.val))) %>%
    spread(key = trait, value = val) -> traitn.stem
  
  #add species and size
  samp.indx <- unique(stemSamples[,c("codeStem","species","size","code")])
  traitn.stem %>%
    left_join(samp.indx) -> traitn.stem
  
  #reorder the columns
  order <- traitcol.order()
  trait.order <- c(order$id.cols, order$phys.cols, order$nutr.cols, order$cfract.cols)
  if(sum(!trait.order %in% colnames(traitn.stem)) == 0){
    traitn.stem %>%
      select(c("codeStem", trait.order)) -> traitn.stem
  }else{
    print("warning: trait columns are out of order")
  }
  
  #note
  #there are 2 stem-level samples that are in the traits df and endophytes df but we're in the deployment df
  # acpa2 and lepa4
  # manually fill in the species and size column for these samples
  tmp <- traitn.stem
  tmp[tmp$codeStem == "acpa2", c("code","species","size")] <- c("acpa","acpa","small")
  tmp[tmp$codeStem == "lepa4", c("code","species","size")] <- c("lepa","lepa","small")
  traitn.stem <- tmp
  
  return(traitn.stem)
  
}






### not used
# load_waterPercent.perGdrymass<-function(){
#   require(dplyr)
#   
#   # read in initial covariate data
#   covar.big <-read.csv('data/covariates_bigStems.csv', stringsAsFactor = F)
#   covar.small <-read.csv('data/covariates_smallStems.csv', stringsAsFactor = F)
#   
#   # (1) calculate water content for each species, size class
#   # g water per g dry mass
#   water.percent <- c(with(covar.big, ((Fresh.mass..g. - Dry.mass..g.) / Dry.mass..g.) * 100),
#                      with(covar.small, ((Fresh.mass..g. - Dry.mass.total..g.) / Dry.mass.total..g.) * 100))
#   water.percent <- data.frame(code=c(covar.big$Species, covar.small$Species),
#                               StemSize=factor(c(rep('large', nrow(covar.big)), rep('small', nrow(covar.small)))),
#                               water.percent, stringsAsFactors=F)
#   
#   ## aggregate by code
#   group_by(water.percent, code) %>%
#     summarize(meanWaterPerc = mean(water.percent, na.rm=TRUE),
#               sdWaterPerc = sd(water.percent, na.rm=TRUE)) -> waterPercent
#   
#   #remove sd cols
#   waterPercent<-waterPercent[,c("code","meanWaterPerc")]
#   colnames(waterPercent)<-c("code","waterperc")
#   length(unique(waterPercent$code)) #34
#   
#   return(waterPercent)
#   
# }

