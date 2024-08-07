% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zscore_and_pvalue.R
\name{zscore_and_pvalue}
\alias{zscore_and_pvalue}
\title{zscore_and_pvalue}
\usage{
zscore_and_pvalue(sce, ppr, switching_genes, cpu = 1)
}
\arguments{
\item{sce}{A Single Cell Experiment object,
containing a matrix of your samples binary gene expression,
which has been filtered to only include switching genes,
using PathPinpointR::subset_switching_genes().}

\item{ppr}{An object of class ppr, must be the same sample.}

\item{switching_genes}{Genes which switch through the trajectory,}

\item{cpu}{Number of cores to use for parallel processing,
(default is one less than available cores, 4 works best).}
}
\value{
the zscore_and_pvalue of the samples raw ppr score.
}
\description{
Calculate the zscore and the p-value of the raw ppr score,
by comparing to 2000 random samples.
}
