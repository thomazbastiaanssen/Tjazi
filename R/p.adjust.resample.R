#'Adjust p-values using resampling
#'@description Wrapper around `p.adjust`, for simple use will do the same.
#'@param p numeric vector of p-values
#'@param method correction method, a character string. Can be abbreviated. see `p.adjust.methods`
#'@param n.samples Integer. How many times should p-values be resampled
#'@return A numeric vector of corrected p-values (of the same length as p, with names copied from p).
#'@export
#'
p.adjust.resamp <- function(p, method = "BH", n.samples = 100){
  stopifnot("'method' should be one of “holm”, “hochberg”, “hommel”, “bonferroni”, “BH”, “BY”, “fdr”, “none”. " = method %in% p.adjust.methods)
  qs <- if(requireNamespace("pbapply", quietly = TRUE)) {
    pbapply::pbsapply(X = p, FUN = function(x){p_resamp_q(x = x, p = p, n = n.samples, method = method)}) }
  else {
    sapply(X = p, FUN = function(x){p_resamp_q(x = x, p = p, n = n.samples, method = method)}) }

  return(qs)
}

#' Compute a single q-value through resampling.
#'@param x the p-value to be corrected
#'@param p numeric vector of the other p-values to sample from.
#'@param method correction method, a character string. Can be abbreviated. see `p.adjust.methods`
#'@param n.samples Integer. How many times should p-values be resampled
#'@return A corrected p-value.
#'
p_resamp_q <- function(x, p, n, method){
  qs <- replicate(n = n, expr =
                    p.adjust(p = c(x, sample(p, size = (length(p) -1), replace = T)),
                             method = method)[1],
                  simplify = T)
  return(median(qs))
}
