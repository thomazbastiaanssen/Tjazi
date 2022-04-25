#' Impute zeroes and perform a centered log-ratio (CLR) transformation
#' @description Microbiome data is compositional. When compositional data is examined using non-compositional methods, many problems arise.
#' Performing a centered log-ratio transformation is a reasonable way to address these problems reasonably well.
#' \cr \cr
#' A major problem with this approach is that microbiome data typically contains lots of zeroes and the logarithm of zero is undefined.
#' Here, we implemented a few methods discussed by Lubbe \emph{et al.} 2021 to replace zeroes with non-zero values in such a way that the structure of the data remains reasonably well preserved.
#' \cr \cr
#' Some of these methods (namely 'logunif' and 'runif') involve imputing small values between 0 and the lowest non-zero value in the dataset.
#' For these methods, we have implemented a resampling approach in order to stabilize the inter-run variability. See \code{method} for more information.
#' @param counts A compositional count table.
#' @param samples_are Either "cols" or "rows". Default is "cols". Denotes whether the columns or rows depict individual samples.
#' @param method The method for zero imputation. One of \code{"logunif"}, \code{"unif"} or \code{"const"}.
#' \code{'logunif'} samples small numbers from a log-uniform distribution, whereas \code{'unif'} samples from a uniform one. On the other hand, \code{"const"} simply replaces zeroes with \code{0.65 * [the lowest value]}.
#' @param replicates An integer. For the two random sampling methods, if this is larger than 1, every zero will be imputed that many times. The median of the CLR of all those replicates will be returned. If \code{method} is set to \code{"const"}, replicates will be automatically set to 1 as no random numbers are generated.
#' @return A CLR-transformed count table.
#' @references Sugnet Lubbe, Peter Filzmoser, Matthias Templ (2021)
#' \emph{Comparison of zero replacement strategies for compositional data with large numbers of zeros.}
#' doi:https://doi.org/10.1016/j.chemolab.2021.104248
#' @export
clr_lite = function(counts, samples_are = "cols", method = "logunif", replicates = 1000)
{
  temp_counts = counts

  if(! method %in% c("logunif", "unif", "const"))
  {stop("`method` must be exactly `logunif`, `unif` or `const`")}

  if(method == "const"){replicates = 1}

  if(samples_are == "rows"){
    temp_counts = data.frame(t(temp_counts))
  }

  temp_counts = apply(X          = temp_counts,
                      MARGIN     = 2,
                      FUN        = clr_imputed,
                      method     = method,
                      replicates = replicates)

  if(samples_are == "rows"){
    temp_counts = data.frame(t(temp_counts))
  }

  clr_counts = data.frame(temp_counts)
  rownames(clr_counts) = rownames(counts)
  colnames(clr_counts) = colnames(counts)
  return(clr_counts)
}

#' compute CLR using Aitchison's method
#' @description See \code{clr_lite}.
#' @seealso \code{\link{clr_lite}}
#' @param x A vector of compositional data without zeroes.
#' @return A vector of CLR-transformed data
#'
anansi_compute_clr = function(x){
  #compute CLR using Aitchison's method
  return(log(x/exp(mean(log(x)))))
}

#' Replace zeroes with non-zero values in order to perform a CLR-transformation
#' @description See \code{\link{clr_lite}}.
#' @seealso \code{\link{clr_lite}}
#' @param vec A vector of compositional data that may contain zeroes.
#' @param method The method for zero imputation. One of "logunif", "unif" or "const".
#' @return A vector with all the zeroes replaced with non-zero values.
#' @importFrom stats runif
#'
impute_zeroes = function(vec, method = "logunif"){
  if(! method %in% c("logunif", "unif", "const")){stop("`method` must be exactly `logunif`, `unif` or `const`")}

  #Find detection limit
  DL = min(vec[vec != 0])
  if(method == "logunif"){
    vec[vec == 0] = DL/(10^(runif(n = sum(vec == 0), min =  0, max = 1)))
  }
  else if(method == "unif"){
    vec[vec == 0] = runif(n = sum(vec == 0), min =  0.1*DL, max = DL)
  }
  else if(method == "const"){
    vec[vec == 0] = 0.65 * DL
  }
  return(vec)
}

#' Resample random values, perform CLR over each iteration and return the median result.
#' @description See \code{clr_lite}.
#' @seealso \code{\link{clr_lite}}
#' @param vec A vector of compositional data that may contain zeroes.
#' @param method The method for zero imputation. One of "logunif", "unif" or "const".
#' @param replicates A positive integer. Default is 1000. Controls how many replicates the median should be taken over.
#' @return a vector of CLR-transformed data
#' @importFrom stats median
#'
clr_imputed = function(vec, method = "logunif", replicates = 1000){
  if(! method %in% c("logunif", "unif", "const")){stop("`method` must be exactly `logunif`, `unif` or `const`")}
  return(apply(replicate(replicates, anansi_compute_clr(impute_zeroes(vec = vec, method = method))), 1, median))
}

#' @rdname clr_lite
#' @section Functions:
#' \code{clr_c:}
#'  A wrapper for \code{clr_lite(counts, method = "const", replicates = 1)}.
#' @export
#'
clr_c <- function(counts, samples_are = "cols"){
  clr_lite(counts, samples_are = samples_are, method = "const", replicates = 1)
}

#' @rdname clr_lite
#' @section Functions:
#' \code{clr_unif:}
#'  A wrapper for \code{clr_lite(counts, method = "unif")}.
#' @export
#'
clr_unif <- function(counts, samples_are = "cols", replicates = 1000){
  clr_lite(counts, samples_are = samples_are, method = "unif", replicates = replicates)
}

#' @rdname clr_lite
#' @section Functions:
#' \code{clr_logunif:}
#'  A wrapper for \code{clr_lite(counts, method = "logunif")}.
#' @export
#'
clr_logunif <- function(counts, samples_are = "cols", replicates = 1000){
  clr_lite(counts, samples_are = samples_are, method = "logunif", replicates = replicates)
}
