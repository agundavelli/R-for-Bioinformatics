#Imports
library(tidyverse)
library(DESeq2)

#' Load a tsv located at specific location `filename` into a tibble
#'
#'
#' @param filename (str): the path to a specific file (ie 'file/path/to/file.tsv')
#'
#' @return tibble: a (g x 1+m) tibble with a 'gene' column followed by
#' sample names as column names.
#'
#' @note Column 'gene' should be first and the only column to contain strings.
#' Data in sample_name columns CANNOT be strings
#'
#' @example `verse_counts <- read_data('verse_counts.tsv')`

read_data <- function(filename){
  verse_counts_tibble <- read_tsv(filename)
  return (verse_counts_tibble) 
}


#' Filter out genes with zero variance
#'
#'
#' @param verse_counts tibble: a (g x 1+m) tibble with a 'gene' column followed
#' by m raw counts columns with sample names as column names.
#'
#' @return tibble: a (n x 1+m) tibble with a 'gene' column followed by m columns
#' of raw counts with genes that have zero variance across samples removed
#'
#' @note (g >= n)
#'
#' @example `filtered_counts <- filter_zero_var_genes(verse_counts)`

filter_zero_var_genes <- function(verse_counts) {
  # calculate variance for each gene across samples, store in a variance column, 
  # filter out the genes with zero variance, and remove variance column
  filtered_verse_counts <- verse_counts %>%
    rowwise() %>%
    mutate(variance = var(c_across(-gene))) %>%
    ungroup() %>% 
    filter(variance > 0) %>%
    select(-variance)
  # Filter out genes with zero variance
  return(filtered_verse_counts)
}


#' Extract time point information from sample name
#'
#'
#' @param str string: sample name from count data.
#'
#' @return string: string character representing sample time point
#'
#' @example `timepoint_from_sample("vAd_1")`
#' output:`"Ad"`

timepoint_from_sample <- function(x) {
  timepoint <- substr(x, 2, 3)
  return(timepoint)
}


#' Grab sample replicate number from sample name
#'
#'
#' @param str  string: sample name from count data.
#'
#' @return string: string character represent sample replicate number
#'
#' @example `sample_replicate("vAd_1")`
#' output: `"1"`

sample_replicate <- function(x) {
  replicate_number = substr(x, 5, 5)
  return(replicate_number)
}


#' Generate sample-level metadata from sample names.
#'
#' Will include columns named "sample", "timepoint", and "replicate" that store
#' sample names, sample time points, and sample replicate, respectively.
#'
#'
#' @param sample_names vector: character vector of length (_S_) consisting of sample
#' names from count data.
#'
#' @return tibble: a (_S_ x 3) tibble with column names "sample",
#' "timepoint", and "replicate". "sample"holds sample_names; "timepoint"
#' stores sample time points; and "replicate" stores sample replicate
#'
#' @note _S_ < m
#'
#' @example `meta <- meta_info_from_labels(colnames(count_data)[colnames(count_data)!='gene'])`

meta_info_from_labels <- function(sample_names) {
  metadata <- tibble(sample = sample_names, 
                     timepoint = sapply(sample_names, timepoint_from_sample), 
                     replicate = sapply(sample_names, sample_replicate))
  return(metadata)
}


#' Calculate total read counts for each sample in a count data.
#'
#'
#' @param count_data tibble: a (n x 1+m) tibble with a 'gene' column followed
#' by m raw counts columns of read counts
#'
#' @return tibble or named vector of read totals from each sample. Vectors must
#' be length `_S_ `, a tibble can be `(1 x _S_)` with sample names as columns
#' names OR `(_S_ x 2)` with columns ("sample", "value")
#'
#' @examples `get_library_size(count_data)`

get_library_size <- function(count_data) {
  read_totals <- count_data %>% 
    select(-gene) %>%
    summarize(across(everything(), ~ sum(.)))
  return(read_totals)
}


#' Normalize raw count data to counts per million WITH pseudocounts using the
#' following formula:
#'     count / (sample_library_size/10^6)
#'
#' @param count_data tibble: a (n x 1+m) tibble with a 'gene' column followed
#' by m raw counts columns of read counts
#'
#' @param count_data tibble: a (n x 1+m) tibble with a 'gene' column followed
#' by m columns of cpm normalized read counts
#'
#' @examples
#' `normalize_by_cpm(count_data)`

normalize_by_cpm <- function(count_data) {
  library_size <- colSums(select(count_data, -gene))
  normalized_data <- count_data %>%
    mutate(across(-gene, ~ . / ((library_size / 1e6)[cur_column()])))
  return(normalized_data)
}

#' Normalize raw count matrix using DESeq2
#'
#' @param count_data tibble: a (n x 1+m) tibble with a 'gene' column followed
#' by m raw counts columns of read counts

#' @param meta_data tibble: sample-level information tibble corresponding to the
#' count matrix columns
#'
#' @return tibble: DESeq2 normalized count matrix
#' @export
#'
#' @examples
#' `deseq_normalize(count_data, meta_data)`
deseq_normalize <- function(count_data, meta_data) {
  # select all columns except the gene column
  count_mat <- as.matrix(count_data[, -1])
  row.names(count_mat) <- count_data$gene
  
  dds <- DESeqDataSetFromMatrix(
    countData = count_mat,
    colData = meta_data,
    design = ~1
  )
  
  dds <- DESeq(dds)
  
  normalized_counts <- counts(dds, normalized = TRUE)
  
  # convert the matrix to a tibble with the 'gene' column
  normalized_counts <- as_tibble(normalized_counts, rownames = "gene")
  
  return(normalized_counts)
}


#' Perform and plot PCA using processed data.
#'
#' PCA is performed over genes, and samples should be colored by time point.
#' Both `y` and `x` axis should have percent of explained variance included.
#'
#'
#' @param data tibble: a (n x _S_) data set
#' @param meta tibble: sample-level meta information (_S_ x 3)
#' @param title string: title for plot
#'
#' @return ggplot: scatter plot showing each sample in the first two PCs.
#'
#' @examples
#' `plot_pca(data, meta, "Raw Count PCA")`

plot_pca <- function(data, meta, title="") {
  pca <- data.frame(prcomp(t(data))$x)
  pca <- pca %>% data.frame(cbind(sample = c(rownames(pca))))
  time1 <- factor(meta$timepoint)
  plot <- ggplot(pca, aes(PC1, PC2)) +
    geom_point(aes(color = time1)) +
    theme(legend.position = "top") +
    ggtitle(title)
  return(plot)
}


#' Plot gene count distributions for each sample using boxplots.
#'
#'
#' @param data tibble: a (n x _S_) data set
#' @param scale_y_axis boolean: whether to scale the `y` axis to log10 values.
#' Default is FALSE, and y-axis will not be transformed.
#' @param title string: title to give the chart.
#'
#' @return ggplot: boxplot show gene count distributions for each sample
#'
#' @example `plot_sample_distributions(data, scale_y_axis=TRUE, title='Raw Count Distributions')`

plot_sample_distributions <- function(data, scale_y_axis=FALSE, title="") {
  data_long <- data %>%
    pivot_longer(cols = everything(), names_to = "sample", values_to = "read_totals")
  p <- ggplot(data_long, aes(x = sample, y = read_totals, fill = sample)) +
    geom_boxplot() +
    labs(
      title = title,
      x = "sample",
      y = "counts"
    ) +
    theme_minimal()
  
  # optionally scale the y-axis to log10
  if (scale_y_axis == TRUE) {
    p <- p + scale_y_continuous(trans = "log10")
  }
  return(p)
}


#' Plot relationship between mean read counts and variability over all genes.
#'
#'
#' @param data tibble: a (n x _S_) data set
#' @param scale_y_axis boolean: whether to scale to y-axis to log10 values. Default
#' is false, and the y-axis will not be transformed.
#' @param title string: title to give the chart.
#'
#' @return ggplot: A scatter plot where the x-axis is the rank of gene ordered by mean
#' count over all samples, and the y-axis is the observed variance of the
#' given gene. Each dot should have their transparency increased. The scatter
#' plot should also be accompanied by a line representing the average mean and
#' variance values.
#'
#' @example `plot_variance_vs_mean(data, scale_y_axis=TRUE, title='variance vs mean (raw counts)')`

plot_variance_vs_mean <- function(data, scale_y_axis=FALSE, title="") {
  # calculate the means and variances of all the numerical columns
  mean_counts <- apply(data[,-1], 1, mean)
  variance_counts <- apply(data[,-1], 1, var)
  df <- data.frame(mean_counts, variance_counts) %>%
    mutate(rank = rank(mean_counts, ties.method = 'average'))

  plot <- ggplot(df, aes(x=rank, y=variance_counts)) +
    geom_point() +
    geom_smooth(method = "gam", color = "red") +
    labs(title = title,
    x = 'Rank(Mean)',
    y = 'Variance') +
    scale_y_continuous(trans = "log10")
  return(plot)
}

