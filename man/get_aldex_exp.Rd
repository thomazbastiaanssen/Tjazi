% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_aldex_exp.R
\name{get_aldex_exp}
\alias{get_aldex_exp}
\title{Extract the CLR-transformed values from the ALDEX2}
\usage{
get_aldex_exp(clr, useMC = TRUE, verbose = TRUE)
}
\description{
Microbiome data is compositional. When compositional data is examined using non-compositional methods, many problems arise.
Performing a centered log-ratio transformation is a reasonable way to address these problems reasonably well. This particular method takes advantage of the excellent ALDEX2 package.
}
