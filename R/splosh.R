########################################################################################################
# This R code implements the SPLOSH method for the analysis of microarray data
# described by Pounds and Cheng (Improving False Discovery Rate Estimation - Bioinformatics 2004).
#
# Last Update: July 14, 2004
#
# Function name: splosh
# Arguments: pvals - a vector of p-values obtained by applying a statistical hypothesis
#                    test to the expression values for each gene individually
#            LC - (optional) a list containing tuning parameters for the loess algorithm,
#                 see help(loess.control) for more details, defaults to S-plus defaults
#            dplaces - (optional) round p-values to this many decimal places to avoid
#                      numerical difficulties that can arise by division by small numbers, default = 3
#            notes - (optional) may pass any string describing notes regarding the analysis, default = NULL
#
# Outputs: a list with the following elements
#         cfdr - a vector of the cfdr estimates at the values specified by the user in pvals
#         q - a vector of corresponding monotone fdr estimates, similar to Storey's q-value
#         cdf - vector of the SPLOSH cdf estimate evaluated at pvals
#         fp - vector of estimated no. of false positives incurred by setting threshold
#              at corresponding value of pvals
#         fn - vector of estimated no. of false negatives incurred by setting threshold
#              at corresponding value of pvals
#         toterr - vector giving the sum of fp and fn for each entry
#         ebp - the empirical Bayes posterior that the null hypothesis is true, see Efron et al. (JASA 2001)
#         pi - estimate of the mixing parameter giving the probability of the null hypothesis
#         ord - a vector giving index numbers to order results according to ascending p-values
#         dplaces - echoes value of dplaces used in algorithm
#         LC - echoes value of LC used in algorithm
#         fdr.method - gives bibliographic reference for the method
#         notes - returns notes regarding the analysis
###########################################################################################################
splosh<-function(pvals,LC=loess.control(),dplaces=3,notes=NULL)

{
  notes<-c(notes,"SPLOSH used to estimate FDR.")
  x<-round(pvals,dplaces)                                      # Prevent numerical difficulties
  n<-length(pvals)                                             # Obtain no. of points
  r<-rank(x)                                                   # rank observations
  ux<-sort(unique(x))                                          # ordered unique p-values
  ur<-(sort(unique(r)))-.5                                     # ordered unique ranks
  #plot(ux,ur)
  if (max(ux)<1)                                               # edge adjustments
  {
    ux<-c(ux,1)
    ur<-c(ur,n)
  }
  if (min(ux)>0)
  {
    ux<-c(0,ux)
    ur<-c(0,ur)
  }
  nux<-length(ux)                                              # No. of unique points
  ur<-ur/n                                                     # Divide ranks by n
  dx<-diff(ux)                                                 # spacings of unique p-values
  dr<-diff(ur)                                                 # spacings of unique ranks
  dFdx<-exp(log(dr)-log(dx))                                   # spacing-interval slopes
  mr<-(ur[1:(nux-1)]+ur[2:nux])/(2)                            # Unitized Medians of Unique Rank Intervals
  mx<-(ux[1:(nux-1)]+ux[2:nux])/2                              # Unitized Medians of Unique p-value Intervals
  tr<-asin(2*(mr-.5))                                          # Arc-Sine Transform Rank Midpoints
  fit<-loess(log(dFdx)~tr,loess.control=LC)                    # Lowess of log-transformed EDQF and transformed rank-midpoints
  tr2<-asin(2*(ur-.5))                                         # Arc-Sine Transform Unique
  yhat<-(predict(fit,tr2))                                     # Obtain Estimated Derivative of CDF (up to constant multiplier)
  yhat[c(1,nux)]<-approx(tr2[2:(nux-1)],
                         yhat[2:(nux-1)],
                         xout=tr2[c(1,nux)],
                         rule=3)$y                             # Extrapolate to get endpoint estimates
  pdf1<-exp(yhat)                                              # Obtain PDF up to a constant
  trap<-0.5*(pdf1[1:(nux-1)]+pdf1[2:nux])*dx                   # Trapezoid rule terms
  const<-sum(trap)                                             # Constant via trapezoid rule
  pdf2<-pdf1/const                                             # Adjust pdf by constant
  pdf<-approx(ux,pdf2,xout=pvals,rule=3)$y                     # Obtain pdf estimate for p-values
  cdf<-approx(ux[2:nux],
              cumsum(trap),
              xout=pvals,
              rule=2)$y/const                                  # Obtain CDF
  pi<-min(pdf)                                                 # Take pi as min of PDF
  p.pi<-max(pvals[pdf==pi])                                    # largest p-value at which pi occurs
  cdf.pi<-max(cdf[pvals==p.pi])                                # CDF evaluated at p.pi
  cfdr<-pi*pvals/cdf                                           # Estimate cFDR for p-values
  ebp<-pi/pdf                                                  # Estimate EBP of null for p-values
  cfdr[pvals==0]<-ebp[pvals==0]                                # L'Hospital's Rule
  cfdr[cfdr>1]<-1                                              # Set max(fdr)=1
  qval<-cfdr                                                   # Initialize q-value vector
  ordering<-order(pvals)                                       # Report Results
  qval[rev(ordering)]<-cummin(qval[rev(ordering)])             # Compute q-value
  fp<-n*cfdr*cdf                                               # Compute Estimated No. False Positives
  fn<-n*(cdf.pi-cdf-pi*(p.pi-pvals))*(pvals<=p.pi)             # Compute Estimated No. False Negatives
  toterr<-fp+fn                                                # Compute Estimated Total No. Errors

  return(list(cfdr=matrix(cfdr,n,1),                 # cfdr estimate
              q=matrix(qval,n,1),                    # q-value estimate
              cdf=matrix(cdf,n,1),                   # cdf estimate
              fp=matrix(fp,n,1),                     # estimate of no. of false postives
              fn=matrix(fn,n,1),                     # estimate of no. of false negatives
              toterr=matrix(toterr,n,1),             # estimated total no. of errors
              pdf=matrix(pdf,n,1),                   # estimate of pdf
              ebp=matrix(ebp,n,1),                   # estimate of ebp
              pi=pi,                                 # estimate of null proportion
              p.pi=p.pi,                             # p-value at which estimated pdf achieves minimum
              cdf.pi=cdf.pi,                         # cdf evaluated at p.pi
              ord=ordering,                          # ordering vector
              dplaces=dplaces,                       # rounded to this no. of decimal places
              LC=LC,                                 # loess.control options used
              fdr.method="Pounds and Cheng (2004)",  # Bibliographic Information
              notes=notes))                          # Notes regarding analysis
}
