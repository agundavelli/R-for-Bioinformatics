library('tidyverse')
library('RColorBrewer')

#' Read the expression data "csv" file as a dataframe, not tibble
#'
#' @param filename (str): the path of the file to read
#' @param delimiter (str): generalize the function so it can read in data with
#'   your choice of delimiter
#'
#' @return A dataframe containing the example intensity data with rows as probes
#'   and columns as samples
#' @export
#'
#' @examples
read_data <- function(intensity_data, delimiter) {
  # use read.csv from base R because it returns a dataframe
  # read_csv returns a tibble
  intensity_df <- read.csv(intensity_data, sep = delimiter)
  return(intensity_df)
}

#' Define a function to calculate the proportion of variance explained by each PC
#'
#' @param pca_results (obj): the results returned by `prcomp()`
#'
#' @return A vector containing the values of the variance explained by each PC
#' @export
#'
#' @examples
calculate_variance_explained <- function(pca_results) {
  variance = pca_results$sdev**2
  explained_variance = variance/sum(variance)
  return(explained_variance)
}

#' Define a function that takes in the variance values and the PCA results to
#' make a tibble with PC names, variance explained by each PC, and the
#' cumulative sum of variance explained. These columns should be named 
#' "principal_components", "variance_explained", and "cumulative", respectively.
#' 
#'
#'
#' @param pca_ve (vector): the vector generated in the previous function with
#'   the variance explained values
#' @param pca_results (object): the results returned by `prcomp()`
#' @return A tibble that contains the names of the PCs, the individual variance
#'   explained, and the cumulative variance explained with names described above
#' @export
#' @examples 
make_variance_tibble <- function(pca_ve, pca_results) {
  pca_tib <- tibble(
    variance_explained = pca_ve,
    principal_components = factor(str_c("PC", 1:length(pca_ve)), str_c("PC", 1:length(pca_ve))),
    cumulative = cumsum(pca_ve)
  )
    return(pca_tib)
}



#' Define a function to create a biplot of PC1 vs. PC2 labeled by
#' SixSubTypesClassification
#'
#' @param metadata (str): The path to the proj_metadata.csv file
#' @param pca_results (obj): The results returned by `prcomp()`
#'
#' @return A ggplot consisting of a scatter plot of PC1 vs PC2 labeled by
#'   SixSubTypesClassification found in the metadata
#' @export
#'
#' @examples
make_biplot <- function(metadata, pca_results) {
  metadata_tib <- read_csv(metadata)
  metadata_tib <- metadata_tib %>% 
    select(geo_accession, SixSubtypesClassification)
  pca_results_tib <- as_tibble(pca_results$x, rownames = "geo_accession") %>% 
    left_join(metadata_tib, by = "geo_accession")
  
  pc_biplot <- ggplot(pca_results_tib, aes(x=PC1, y=PC2, color = SixSubtypesClassification, main = "PC1 vs PC2 Biplot")) +
    geom_point() 
  
  return(pc_biplot)
}

#' Define a function to return a list of probeids filtered by signifiance
#'
#' @param diff_exp_tibble (tibble): A tibble containing the differential expression results
#' @param fdr_threshold (float): an appropriate FDR threshold, we will use a
#'   value of .01. This is the column "padj" in the tibble.
#'
#' @return A list with the names of the probeids passing the fdr_threshold
#' @export
#'
#' @examples
list_significant_probes <- function(diff_exp_tibble, fdr_threshold) {
  filtered_diff_exp_tibble <- diff_exp_tibble %>% filter(padj < fdr_threshold)
  return(filtered_diff_exp_tibble$probeid)
}

#' Define a function that uses the list of significant probeids to return a
#' matrix with the intensity values for only those probeids.
#' @param intensity (dataframe): The dataframe of intensity data generated in
#'   part 1
#' @param sig_ids_list (list/vector): The list of differentially expressed
#'   probes generated in part 6
#'
#' @return A `matrix()` of the probe intensities for probes in the list of
#'   significant probes by FDR determined in the previous function.
#'
#' @export
#'
#' @examples
return_de_intensity <- function(intensity, sig_ids_list) {
  # dataframe of the significant probes by FDR in the probe intensities data
  filtered_intensity_matrix <- intensity %>%
    filter(rownames(intensity) %in% sig_ids_list)
  filtered_intensity_matrix <- as.matrix(filtered_intensity_matrix)
  return(filtered_intensity_matrix)
}

#' Define a function that takes the intensity values for significant probes and
#' creates a color-blind friendly heatmap
#'
#' @param de_intensity (matrix): The matrix of intensity values for significant
#'   differentially expressed probes returned in part 7
#' @param num_colors (int): The number of colors in a specificed RColorBrewer
#'   palette
#' @param palette (str): The name of the chosen RColorBrewer palette
#'
#' @return A heatmap displaying the intensity values for the differentially
#'   expressed probes
#' @export
#'
#' @examples
plot_heatmap <- function(de_intensity, num_colors, palette) {
  intensity_heatmap <- heatmap(de_intensity, col = brewer.pal(num_colors, palette), main = "Heatmap of Intensity Values")
    return(intensity_heatmap)
}

