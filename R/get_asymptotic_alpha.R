#' Functional wrapper to get the output of the wonderful iNext library in the format I use to pipe into ggplot2.
#' @export

get_asymptotic_alpha = function(species, groups = c(), verbose = TRUE ){

  #create output df
  alpha_diversity <- data.frame(chao1              = rep(NA, ncol(species)),
                                asymptotic_simps   = rep(NA, ncol(species)),
                                asymptotic_shannon = rep(NA, ncol(species)))
  row.names(alpha_diversity) = colnames(species)


  for(samp in 1:ncol(species)){

    raw_observation                            = species[,samp]
    observations                               = raw_observation[raw_observation!= 0]
    if(verbose){
      print(paste("Processing sample ", samp, "of ", ncol(species)))
    }
    alpha_diversity[samp,"chao1"]              = iNEXT::ChaoSpecies(x = observations,datatype = "abundance")$Estimator
    alpha_diversity[samp,"asymptotic_simps"]   = iNEXT::EstSimpson( x = observations,datatype = "abundance")$Estimator
    alpha_diversity[samp,"asymptotic_shannon"] = iNEXT::ChaoEntropy(x = observations,datatype = "abundance")$Estimator
    if(verbose){
      print(paste("finished with sample ", samp, " out of ", ncol(species)))
    }
  }
  if(length(groups) == ncol(species)){
    alpha_diversity <-(cbind(alpha_diversity, group = groups))
  }
  return(alpha_diversity)
}

