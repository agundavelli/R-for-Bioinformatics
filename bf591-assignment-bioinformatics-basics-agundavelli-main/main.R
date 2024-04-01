#!/usr/bin/Rscript
## Author: Taylor Falk
## tfalk@bu.edu
## BU BF591
## Assignment Week 2

#### Bioconductor ####
# it is standard among R packages to define libraries and packages at the 
# beginning of a script. Also note that a package should NOT be installed every 
# time a script runs.
# The bioconductor repository has installation instructions for biomaRt: 
# https://bioconductor.org/install/

if (!require("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}
if (!require("biomaRt", quietly = TRUE)){
  BiocManager::install("biomaRt")
}
library(biomaRt)
library(tidyverse)

#### Loading and processing data ####
#' Load Expression Data
#'
#' @param filepath A text string of the full filepath to the file to load.
#'
#' @return A tibble containing the data loaded from the CSV in `filepath`. 
#' 
#' @examples 
#' `data <- load_expression('data/example_intensity_data_subset.csv')`
load_expression <- function(filepath) {
  intensity_data <- readr::read_csv(filepath)
  return(intensity_data)
}

#' Filter 15% of the gene expression values.
#'
#' @param tibble A tibble of expression values, rows by probe and columns by sample.
#'
#' @return A tibble of affymetrix probe names from the input expression data tibble. 
#' These names match the rows with 15% or more of the expression values about log2(15).
#' 
#' @details This is similar to the filters being implemented in BF528's project 1. 
#' We do not necessarily want to capture all parts of the assay in our analysis, so 
#' filters like this serve to reduce the noise and amount of data to examine.
#'
#' @examples `samples <- filter_15(data_tib)`
#' 
filter_15 <- function(tibble){
  # create new column to store number of entries in each row that are greater 
  # than log2(15)
  # select each row and count number of entries greater than log2(15) and store
  # in column
  # 
  tibble <- mutate(tibble, Number_of_Entries = rowSums(dplyr::select(tibble, -1) > log2(15)))
  tibble <- mutate(tibble, Percent_of_Entries = Number_of_Entries / (ncol(tibble)-2))
  filtered_tibble <- tibble %>% filter(Percent_of_Entries > 0.15)
  filtered_tibble <- filtered_tibble %>% 
    select(-all_of(c("Number_of_Entries", "Percent_of_Entries")))
  return(filtered_tibble)
}

#### Gene name conversion ####

#' Convert affymetrix array names into hgnc_symbol IDs using biomaRt. Inputs and 
#' outputs will likely not be the same size.
#'
#' @param affy_tib A single column tibble of strings containing array names.
#'
#' @return A 2 column tibble that contains affy IDs in the first column,
#' and their corresponding HGNC gene ID in the second column. Note that not all affy IDs 
#' will necessarily correspond with a gene ID, and one gene may have multiple affy IDs.
#' 
#' @details Connecting to ensembl via biomaRt can be...hit or miss...so you may 
#' want to check if data was correctly returned (or if it was just empty). The 
#' `getBM()` function may not accept a tibble, so you might need to convert your 
#' input into a flat vector.
#'
affy_to_hgnc <- function(affy_vector) {
  ensembl <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
  tibble_2 <- as_tibble(getBM(attributes = c('affy_hg_u133_plus_2', 'hgnc_symbol'),
    filters = 'affy_hg_u133_plus_2',
    values = pull(affy_vector),
    mart = ensembl))
  return(tibble_2)
}

#### ggplot ####

#' Reduce a tibble of expression data to only the rows in good_genes or bad_genes.
#'
#' @param expr_tibble A tibble holding the expression data, each row corresponding to
#' one affymetrix probe ID and each column to a sample.
#' @param names_ids A two column tibble that associates affy IDs with HGNC gene IDs. 
#' Generated `with affy_to_hgnc()`.
#' @param good_genes A list of gene names stored as a vector of strings.
#' @param bad_genes A list of gene names stored as a vector of strings.
#'
#' @return A tibble with two additional columns added:
#' 1. HGNC gene IDs 
#' 2. Does the gene is this row fall into "good" or "bad" genes?
#' This tibble should be reduced to only rows present in good or bad genes. All
#' other rows can be discarded.
#' 
#' @details In order to plot only our genes of interest, we need to rearrange our 
#' data to include only the elements we want to see. We also want to add to columns, 
#' one that associates the probeids with the HGNC gene name, and one that says if 
#' that gene is in the good or bad sets of genes.
#'
reduce_data <- function(expr_tibble, names_ids, good_genes, bad_genes){
  names_ids <- names_ids %>%
    mutate(gene_set = NA_character_) %>%
    mutate(gene_set = case_when(hgnc_symbol %in% good_genes ~ "good", hgnc_symbol %in% bad_genes ~ "bad", TRUE ~ gene_set))
  
  expr_tibble <- left_join(expr_tibble, names_ids, by = c("probe" = "affy_hg_u133_plus_2")) %>%
    select(probe, hgnc_symbol, gene_set, everything())
  
  expr_tibble <- filter(expr_tibble, !is.na(gene_set))
  
  return(expr_tibble)
}


#' Convert the tibble from wide to long format.
#'
#' @param tibble A tibble of expression data in wide format with information about
#' good and bad genes, gene names, sample names, and expression values.
#'
#' @return A tibble in long format containing the same information.
#'
#' @details This function's primary objective is to reformat the tibble from a 
#' wide format, where there are separate columns for sample name and expression value, 
#' to a long format, where sample names are stored in a column named 'sample' 
#' and their values stored in another column named 'value'.
#'
convert_to_long <- function(tibble) {
  # Convert wide format to long format
  long_tibble <- tibble %>% tidyr::pivot_longer(
    c(GSM1, GSM2),
    names_to = "sample",
    values_to = "value"
  )
  return(long_tibble)
}


