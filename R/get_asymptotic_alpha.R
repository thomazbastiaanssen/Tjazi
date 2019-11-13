#Functional wrapper to get the output of the wonderful iNext library in the format I use to pipe into ggplot2. 

library(iNEXT)
get_asymptotic_alpha = function(species, groups = c() ){

  #create output df
  alpha_diversity <- data.frame(chao1              = rep(NA, ncol(species)), 
                                asymptotic_simps   = rep(NA, ncol(species)), 
                                asymptotic_shannon = rep(NA, ncol(species)))
  row.names(alpha_diversity) = colnames(species)
  
  
  for(samp in 1:ncol(species)){
    
    raw_observation                            = species[,samp]
    observations                               = raw_observation[raw_observation!= 0]
    print(paste("Processing sample ", samp, "of ", ncol(species)))
    alpha_diversity[samp,"chao1"]              = ChaoSpecies(x = observations,datatype = "abundance")$Estimator
    alpha_diversity[samp,"asymptotic_simps"]   =  EstSimpson(x = observations,datatype = "abundance")$Estimator
    alpha_diversity[samp,"asymptotic_shannon"] = ChaoEntropy(x = observations,datatype = "abundance")$Estimator
    print(paste("finished with sample ", samp, " out of ", ncol(species)))
  }
  if(length(metadata) == ncol(species)){
    alpha_diversity <-(cbind(alpha_diversity, group = groups))
  }
  return(alpha_diversity)
}

