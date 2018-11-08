guess_counts <- function(vect, lowestCount = 1){
  #   Try to restore a compositional vector to its former count glory. 
  #   *WARNING*
  #   In the rare case that the lowest count wasn't zero, this will fail. 
  if(min(vect != 0)){
    output <- round(vect*lowestCount/min(vect))
    cat("There don't seem to be any zeroes in this batch, multiplied all values by ",
        lowestCount/min(vect),
        ".\n",
        sep = "")
  }
  else{
    lowest <- unique(vect[order(vect)])[2]
    output <- round(vect*lowestCount/lowest)
    cat("There were zeroes in this batch, multiplied all values by ",
        lowestCount/lowest,
        ".\n", 
        sep = "")
  }
  return(output)
}
