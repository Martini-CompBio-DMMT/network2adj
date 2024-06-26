% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{writeAdjFile}
\alias{writeAdjFile}
\title{Get the adjacency matrix file as requested by VIPER}
\usage{
writeAdjFile(df_list, file_name = "")
}
\arguments{
\item{df_list}{list of genes with their interactors}

\item{file_name}{how to nominate the file}
}
\value{
.adj file
}
\description{
This function is used to obtain the adjacency matrix file in the format
requested by VIPER to obtain the regulon with the "aracne2regulon" function.
}
\examples{

# Load the list of regulators
data(regulators_list)

# Load the example data frame
data(network_df)

# Get the gene symbols corresponding to the gene IDs in the data frame
gene_symbols <- getGeneSymbols(network_df, organism = "hsapiens")

# Get the converted data frame
converted_df <- getGeneSymbolDataframe(network_df, gene_symbols,
organism = "hsapiens")

# Get the list of the regulators with their interactors
genes <- c("ADAR", "ADARB1")
interactors <- getInteractorsList(converted_df, regulators_list, genes)

writeAdjFile(interactors, file_name = "interactors")

}
