#check for missing data

#issue 31 -- There are 2 unique codeStem ids that are found in the trait data (xrf sample names) and the sequence data, but not in the stemSamples data (deployment sample names).  
#These codeStem ids are not found in the percent mass loss data.  
#Because the main goal is to analyze decay responses, I'm going to continue to leave these codeStems out of the stemSamples dataframe. 
#Is it possible that stem id numbers got switched? Something to follow-up on.

checkMissing_stems <- function(traits.stem, seqSamples, pmr_byStem){
  
  traits.stem %>%
    filter(!codeStem %in% stemSamples$codeStem) #there are 2 stems that are in the traits df but are missing from the stem lookup table...
  
  seqSamples %>%
    filter(!codeStem %in% stemSamples$codeStem) %>%
    separate(seq_sampName, into=c("drop","seq.stem"), 4, remove=FALSE) %>%
    filter(grepl("[1-9]", seq.stem)) -> result #these same 2 stems are in the sequence samples but are missing from the stem lookup table...
  result
  
  pmr_byStem %>%
    filter(!codeStem %in% stemSamples$codeStem) # no stem samples are missing from the decay data
  
  return(result)
}

