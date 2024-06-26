% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{getGeneSymbols}
\alias{getGeneSymbols}
\title{Get gene symbols}
\usage{
getGeneSymbols(dataf, organism = "")
}
\arguments{
\item{dataf}{data frame obtained by the network with igraph}

\item{organism}{hsapiens or mmusculus}
}
\value{
data frame
}
\description{
This function is used to get the human gene symbols from the gene Ensembl
IDs in the input data frame.
}
\examples{

# Load the example data frame
data(network_df)

# Get the gene symbols corresponding to the gene IDs in the data frame
gene_symbols <- getGeneSymbols(network_df, organism = "hsapiens")

}
