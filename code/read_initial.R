


read_in_initial_mass <- function(){
  library(readr)
  library(dplyr)
  big <- read_csv("data/covariates_bigStems.csv")
  small <- read_csv("data/covariates_smallStems.csv")
 
  big_out <- process_initial_file(big,"big")
  small_out <- process_initial_file(small,"small")
  
  df_out<-bind_rows(big_out,small_out)
  return(df_out)
} 
  
  
  
process_initial_file<-function(df,size){
  
  if ("Dry mass total (g)" %in% names(df)) {
    df <- rename(df,`Dry mass (g)`=`Dry mass total (g)`)
  }
  
  df %>%
     mutate(dry_mass_content=`Dry mass (g)`/`Fresh mass (g)`) %>%
     filter(!is.na(dry_mass_content)) %>%
     group_by(Species) %>%
     summarize(dry_mass_prop=mean(dry_mass_content,na.rm=T),n()) -> moisture
  
  df %>%
    left_join(moisture) %>%
    mutate(totalSampleDryMass=`Fresh mass (g)`*dry_mass_prop,size=size,density=NA,time=0,fruiting=NA,insects=NA,drill=NA) %>%
    select(unique, Species, size,time,totalSampleDryMass,density,fruiting,insects,drill) -> df_out
  
  return(df_out)
}
