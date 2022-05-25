#' Functional wrapper to get the output of the wonderful iNext library in the format I use to pipe into ggplot2.
#' @export
get_asymptotic_alpha = function(species, verbose = TRUE){

  if(any(colSums(species) < 1)){
    if(verbose){print("reads seem to be relative abundance. This is not ideal. Applying transformation to address this.")}
    species <- apply(species, 2, un_cpm)
  }
  #create output df, use progress bar if available

  alpha_diversity <- if(requireNamespace("pbapply", quietly = TRUE)) {
    pbapply::pbapply(X = species,
                     MARGIN =  2,
                     FUN = div_applier) } else {
               apply(X = species,
                     MARGIN =  2,
                     FUN = div_applier)
                     }
  alpha_diversity = data.frame(alpha_diversity)

  colnames(alpha_diversity)  <- colnames(species)

  alpha_diversity = cbind("Index" = c("Chao1", "Simpson Index", "Shannon Entropy"), alpha_diversity)

  return(alpha_diversity)
}

#' Functional wrapper to get the output of the wonderful iNext library in the format I use to pipe into ggplot2.
#'
div_applier <- function(x){
  out_div <- vector(mode = "numeric", length = 3)

  out_div[1] <- iNEXT::ChaoRichness(x)$Estimator
  out_div[2] <- iNEXT::ChaoSimpson( x)$Estimator
  out_div[3] <- iNEXT::ChaoEntropy( x)$Estimator

  return(out_div)
}


#' Undo CPM transformation
#'
un_cpm <- function(x) 1E6 * x/sum(x)

