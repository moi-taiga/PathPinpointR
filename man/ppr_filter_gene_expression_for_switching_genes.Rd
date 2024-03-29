% Generated by roxygen2: do not edit by hand
% Please edit documentation in
%   R/ppr_filter_gene_expression_for_switching_genes.R
\name{ppr_filter_gene_expression_for_switching_genes}
\alias{ppr_filter_gene_expression_for_switching_genes}
\title{ppr_filter Gene Expression for Switching Genes}
\usage{
ppr_filter_gene_expression_for_switching_genes(
  binary_counts_matrix,
  reference.sg
)
}
\arguments{
\item{binary_counts_matrix}{a binary expression matrix from your sample.}

\item{reference.sg}{Genes which switch through the trajectory as identified by GeneSwitches.}
}
\value{
a reduced binary expression matrix filtered to only include selected switching genes
}
\description{
Create a reduced binary expression matrix for only the selected switching genes,
binary_counts_matrix is from the sample DATA and gs_scorer_genes is from Atlas Data.
}
