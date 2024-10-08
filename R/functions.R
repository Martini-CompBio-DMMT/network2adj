#' Get the data frame from the network
#'
#' This function is used to get a data frame from the file obtained by the
#' network reconstruction with ARACNe.
#'
#' @usage getDf(net)
#' @param net network obtained from ARACNe
#' @return data frame
#' @examples
#'
#' # Load the example network
#' data(network)
#'
#' # Get the data frame from the network obtained by ARACNe
#' net_df <- getDf(network)
#'
#' @importFrom igraph graph_from_adjacency_matrix as_data_frame
#' @export

getDf <- function(net) {

  g  <- graph_from_adjacency_matrix(net, weighted=TRUE, mode="undirected")
  dataf <- as_data_frame(g)

  dataf

}

#' Get gene symbols
#'
#' This function is used to get the human gene symbols from the gene Ensembl
#' IDs in the input data frame.
#'
#' @usage getGeneSymbols(dataf, organism = "")
#' @param dataf data frame obtained by the network with igraph
#' @param organism hsapiens or mmusculus
#' @return data frame
#' @examples
#'
#' # Load the example data frame
#' data(network_df)
#'
#' # Get the gene symbols corresponding to the gene IDs in the data frame
#' gene_symbols <- getGeneSymbols(network_df, organism = "hsapiens")
#'
#' @importFrom gprofiler2 gconvert gorth
#' @export

getGeneSymbols <- function(dataf, organism = "") {

  all_genes <- unique(c(dataf$from, dataf$to))
  if (organism == "hsapiens") {
    gene_symbols <- gconvert(all_genes, organism = "hsapiens",
                             target="ENTREZGENE", filter_na = F)
  }

  if (organism == "mmusculus") {
    gene_symbols <- gorth(query = all_genes, source_organism = "mmusculus",
                          target_organism = "hsapiens", filter_na = F)
  }

  gene_symbols

}

#' Convert gene Ensembl IDs in the corresponding gene symbols
#'
#' This function is used to convert the gene IDs of the data frame in the
#' corresponding gene symbols.
#'
#' @usage getGeneSymbolDataframe(dataf, symbols, organism = "")
#' @param dataf data frame obtained by the network
#' @param symbols data frame with the gene IDs and the corresponding symbols
#' @param organism hsapiens or mmusculus
#' @return data frame
#' @examples
#'
#' # Load the example data frame
#' data(network_df)
#'
#' # Get the gene symbols corresponding to the gene IDs in the data frame
#' gene_symbols <- getGeneSymbols(network_df, organism = "hsapiens")
#'
#' # Get the converted data frame
#' converted_df <- getGeneSymbolDataframe(network_df, gene_symbols,
#' organism = "hsapiens")
#'
#' @export

getGeneSymbolDataframe <- function(dataf, symbols, organism = "") {

  if (organism == "hsapiens") {
    input_indices_1 <- match(dataf[, 1], symbols$input)
    input_indices_2 <- match(dataf[, 2], symbols$input)
    target_genes_1 <- list(symbols$target[input_indices_1])
    target_genes_2 <- list(symbols$target[input_indices_2])
    dataf[, 1] <- target_genes_1
    dataf[, 2] <- target_genes_2
    dataf <- na.omit(dataf)
  }

  if (organism == "mmusculus") {
    input_indices_1 <- match(dataf[, 1], symbols$input)
    input_indices_2 <- match(dataf[, 2], symbols$input)
    target_genes_1 <- list(symbols$ortholog_name[input_indices_1])
    target_genes_2 <- list(symbols$ortholog_name[input_indices_2])
    dataf[, 1] <- target_genes_1
    dataf[, 2] <- target_genes_2

    dataf <- na.omit(dataf)
    dataf <- dataf[-grep("N/A", dataf$from), ]
    dataf <- dataf[-grep("N/A", dataf$to), ]
  }

  dataf

}

#' Convert gene Ensembl IDs in the corresponding gene symbols
#'
#' This function is used to convert the gene IDs of the expression matrix in
#' the corresponding gene symbols.
#'
#' @usage getGeneSymbolMatrix(vsd, symbols, organism = "")
#' @param vsd matrix containing the normalized expression values
#' @param symbols data frame with the gene IDs and the corresponding symbols
#' @param organism hsapiens or mmusculus
#' @return matrix to be used with VIPER
#' @examples
#'
#' # Load the example data frame
#' data(network_df)
#'
#' # Load the example matrix
#' data(vsd_matrix)
#'
#' # Get the gene symbols corresponding to the gene IDs in the data frame
#' gene_symbols <- getGeneSymbols(network_df, organism = "hsapiens")
#'
#' # Get the converted matrix
#' converted_matrix <- getGeneSymbolMatrix(vsd_matrix, gene_symbols,
#' organism = "hsapiens")
#'
#' @export

getGeneSymbolMatrix <- function(vsd, symbols, organism = "") {

  if (organism == "hsapiens") {
    for (i in 1:nrow(vsd)) {
      if (row.names(vsd)[i] %in% symbols$input) {
        newname <- symbols$target[which(symbols$input == row.names(vsd)[i])][1]
        if (is.na(newname) == FALSE) {
          if (!(newname %in% row.names(vsd))) {
            row.names(vsd)[i] <- newname
          }
          else {
            vsd <- (vsd[-i,])
          }
        }
        else {
          vsd <- (vsd[-i,])
        }
      }
    }

    vsd <- vsd[-grep("^ENSG", row.names(vsd)), ]
  }

  if (organism == "mmusculus") {
    for (i in 1:nrow(vsd)) {
      if (row.names(vsd)[i] %in% symbols$input) {
        newname <- symbols$ortholog_name[which(symbols$input
                                               == row.names(vsd)[i])][1]
        if (is.na(newname) == FALSE) {
          if (!(newname %in% row.names(vsd))) {
            row.names(vsd)[i] <- newname
          }
          else {
            vsd <- (vsd[-i,])
          }
        }
        else {
          vsd <- (vsd[-i,])
        }
      }
    }
    vsd <- vsd[-grep("N/A", row.names(vsd)), ]
  }

  mat <- as.matrix(vsd)
  mat

}

#' Get a list of genes with their interactors
#'
#' This function extracts gets the list of regulators, and it adds the genes
#' for which we want to compute the activity with VIPER.
#' Then it filters the data frame obtained by the network to consider only the
#' regulators of the gene of interest and, for each regulator, it gets a list of
#' their interactors and the corresponding scores.
#'
#' @usage getInteractorsList(dataf, reg_list, target_genes)
#' @param dataf data frame obtained by the network
#' @param reg_list list with the regulator genes
#' @param target_genes list of the genes of interest
#' @return list
#' @examples
#'
#' # Load the list of regulators
#' data(regulators_list)
#'
#' # Load the example data frame
#' data(network_df)
#'
#' # Get the gene symbols corresponding to the gene IDs in the data frame
#' gene_symbols <- getGeneSymbols(network_df, organism = "hsapiens")
#'
#' # Get the converted data frame
#' converted_df <- getGeneSymbolDataframe(network_df, gene_symbols,
#' organism = "hsapiens")
#'
#' # Get the list of the regulators with their interactors
#' genes <- c("ADAR", "ADARB1")
#' interactors <- getInteractorsList(converted_df, regulators_list, genes)
#'
#' @importFrom dplyr filter
#' @export

getInteractorsList <- function(dataf, reg_list, target_genes) {

  genes_selected <- c(reg_list, target_genes)

  df_new <- data.frame()
  for (g in genes_selected) {
    df1 <- filter(dataf, dataf[, 1] == g)
    df2 <- filter(dataf, dataf[, 2] == g)
    df2 <- df2[, c(2, 1, 3)]
    colnames(df2) <- c("from", "to", "weight")
    df3 <- rbind(df1, df2)
    df_new <- rbind(df_new, df3)
  }

  df_list <- split(df_new, df_new$from)
  df_list <- lapply(df_list, function(dataf) {dataf[, -1]})

  lengths <- c()
  for (tg in target_genes) {
    lengths <- c(lengths, dim(df_list[[tg]])[1])
  }

  min_length <- min(lengths)

  df_list <- Filter(function(dataf) nrow(dataf) >= min_length, df_list)
  df_list

}

#' Get the adjacency matrix file as requested by VIPER
#'
#' This function is used to obtain the adjacency matrix file in the format
#' requested by VIPER to obtain the regulon with the "aracne2regulon" function.
#'
#' @usage writeAdjFile(df_list, file_name = "")
#' @param df_list list of genes with their interactors
#' @param file_name how to nominate the file
#' @return .adj file
#' @examples
#'
#' # Load the list of regulators
#' data(regulators_list)
#'
#' # Load the example data frame
#' data(network_df)
#'
#' # Get the gene symbols corresponding to the gene IDs in the data frame
#' gene_symbols <- getGeneSymbols(network_df, organism = "hsapiens")
#'
#' # Get the converted data frame
#' converted_df <- getGeneSymbolDataframe(network_df, gene_symbols,
#' organism = "hsapiens")
#'
#' # Get the list of the regulators with their interactors
#' genes <- c("ADAR", "ADARB1")
#' interactors <- getInteractorsList(converted_df, regulators_list, genes)
#'
#' writeAdjFile(interactors, file_name = "interactors")
#'
#' @import utils
#' @export

writeAdjFile <- function(df_list, file_name = "") {

  file_txt <- file.path(getwd(), paste0(file_name, ".txt"))
  file_adj <- file.path(getwd(), paste0(file_name, ".adj"))

  if (file.exists(file_txt)) {
    file.remove(file_txt)
  }

  if (file.exists(file_adj)) {
    file.remove(file_adj)
  }

  i = 1
  file(file_txt, "w")
  while (i < length(df_list)) {
    a <- c(names(df_list[i]), paste(df_list[[i]]$to, df_list[[i]]$weight))
    b <- write.table(x = rbind(a), file = file_txt, row.names = F,
                     col.names = F, quote = F, append = T)
    i = i+1
  }

  txt_file <- readLines(file_txt)
  adj_file <- file_adj
  modified_file <- gsub(" ", "\t", txt_file)
  writeLines(modified_file, adj_file)

}
