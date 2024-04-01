  library('tidyverse')
  library('SummarizedExperiment')
  library('DESeq2')
  library('biomaRt')
  library('testthat')
  #BiocManager::install("fgsea")
  library("fgsea")

#' Function to generate a SummarizedExperiment object with counts and coldata
#' to use in DESeq2
#'
#' @param csv_path (str): path to the file verse_counts.tsv
#' @param metafile (str): path to the metadata sample_metadata.csv
#' @param selected_times (list): list of sample timepoints to use
#' 
#'   
#' @return SummarizedExperiment object with subsetted counts matrix
#'   and sample data. Ensure that the timepoints column used as input 
#'   to the model design has 'vP0' set as the reference factor level. Your 
#'   colData dataframe should have columns named samplename and timepoint.
#' @export
#'
#' @examples se <- make_se('verse_counts.tsv', 'sample_metadata.csv', c('vP0', 'vAd'))
make_se <- function(counts_csv, metafile_csv, selected_times) {
  # read in counts_csv file
  data <- as.data.frame(readr::read_tsv(counts_csv))
  counts <- data[, names(data) != "gene"]
  row.names(counts) <- data$gene
  
  # read in metadata
  metadata <- as.data.frame(readr::read_csv(metafile_csv))
  
  # filter metadata to only have the selected times
  selected_timepoints <- metadata[metadata$timepoint %in% selected_times, ]
  colData <- DataFrame(samplename = selected_timepoints$samplename,
                       timepoint = selected_timepoints$timepoint,
                       row.names = selected_timepoints$samplename)
  
  # set reference level to vP0
  colData$timepoint <- relevel(factor(coldata$timepoint), ref = "vP0")

  # convert counts into matrix
  counts <- as.matrix(counts[, selected_timepoints$samplename])
  
  # create SummarizedExperiment object with counts stored as an assay
  se <- SummarizedExperiment(
  assays = list(counts = counts),
  colData = colData)
  
  return(se)
}

#' Function that runs DESeq2 and returns a named list containing the DESeq2
#' results as a dataframe and the dds object returned by DESeq2
#'
#' @param se (obj): SummarizedExperiment object containing counts matrix and
#' coldata
#' @param design: the design formula to be used in DESeq2
#'
#' @return list with DESeqDataSet object after running DESeq2 and results from
#'   DESeq2 as a dataframe
#' @export
#'
#' @examples results <- return_deseq_res(se, ~ timepoint)
return_deseq_res <- function(se, design) {
  dds <- DESeqDataSet(se, design = design)
  dds <- DESeq(dds)
  deseq_results <- results(dds)
  deseq_list <- list(results = as.data.frame(deseq_results), dds = dds)
  return(deseq_list)
}

#' Function that takes the DESeq2 results dataframe, converts it to a tibble and
#' adds a column to denote plotting status in volcano plot. Column should denote
#' whether gene is either 1. Significant at padj < .10 and has a positive log
#' fold change, 2. Significant at padj < .10 and has a negative log fold change,
#' 3. Not significant at padj < .10. Have the values for these labels be UP,
#' DOWN, NS, respectively. The column should be named `volc_plot_status`.
#'
#' @param deseq2_res (df): results from DESeq2 
#' @param padj_threshold (float): threshold for considering significance (padj)
#'
#' @return Tibble with all columns from DESeq2 results and one additional column
#'   labeling genes by significant and up-regulated, significant and
#'   downregulated, and not significant at padj < .10.
#'   
#' @export
#'
#' @examples labeled_results <- label_res(res, .10)
label_res <- function(deseq2_res, padj_threshold) {
  deseq2_res_tib <- as_tibble(deseq2_res)
  deseq2_res_tib <- deseq2_res_tib %>% 
    mutate(genes = rownames(deseq2_res)) %>%
    mutate(volc_plot_status = case_when(
      padj < padj_threshold & log2FoldChange > 0 ~ "UP",
      padj < padj_threshold & log2FoldChange < 0 ~ "DOWN",
      TRUE ~ "NS"
    ))       
    return(deseq2_res_tib)
}

#' Function to plot the unadjusted p-values as a histogram
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one additional
#' column denoting status in volcano plot
#'
#' @return ggplot: a histogram of the raw p-values from the DESeq2 results
#' @export
#'
#' @examples pval_plot <- plot_pvals(labeled_results)
plot_pvals <- function(labeled_results) {
  pval_plot <- ggplot(data = labeled_results, aes(x = pvalue)) +
    geom_histogram(binwidth = 0.02, color = "black", fill = "lightblue") +
    labs(title = "Histogram of raw p-values obtained from DE analysis (vP0 vs vAd)", 
         x = "p-value", 
         y = "count") +
    theme_minimal()
  
  return(pval_plot)
}

#' Function to plot the log2foldchange from DESeq2 results in a histogram
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one additional
#' column denoting status in volcano plot
#' @param padj_threshold (float): threshold for considering significance (padj)
#'
#' @return ggplot: a histogram of log2FC values from genes significant at padj 
#' threshold of 0.1
#' @export
#'
#' @examples log2fc_plot <- plot_log2fc(labeled_results, .10)
plot_log2fc <- function(labeled_results, padj_threshold) {
  significant_genes <- labeled_results %>% filter(padj < padj_threshold)
  log2fc_plot <- ggplot(significant_genes, aes(x = log2FoldChange)) +
    geom_histogram(binwidth = 0.15, fill = "lightblue", color = "black") +
    labs(title = "Histogram of Log2FoldChanges for DE Genes (vP0 vs. vAd)",
         x = "log2FoldChange",
         y = "count") +
    theme_minimal()
  
  return(log2fc_plot)
}

#' Function to make scatter plot of normalized counts for top ten genes ranked
#' by ascending padj
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one
#'   additional column denoting status in volcano plot
#' @param dds_obj (obj): The object returned by running DESeq (dds) containing
#' the updated DESeqDataSet object with test results
#' @param num_genes (int): Number of genes to plot
#'
#' @return ggplot: a scatter plot with the normalized counts for each sample for
#' each of the top ten genes ranked by ascending padj
#' @export
#'
#' @examples norm_counts_plot <- scatter_norm_counts(labeled_results, dds, 10)
scatter_norm_counts <- function(labeled_results, dds_obj, num_genes){
  top_genes <- arrange(labeled_results, padj) 
  top_genes <- top_genes %>%  head(num_genes)
  
  top_gene_names <- top_genes$genes
  
  subset_dds <- dds_obj[2]
  dds <- estimateSizeFactors(dds_obj[[2]])
  
  counts_df <- as.data.frame(counts(dds, normalized = TRUE))
  
  counts_df <- counts_df %>%
    rownames_to_column(var = "Sample") %>%
    pivot_longer(-Sample, names_to = "Gene", values_to = "Normalized_Counts")
  scatter_plot <- ggplot(counts_df, aes(x = Sample, y = Normalized_Counts, color = Gene)) +
    geom_point() +
    labs(
      title = "Plot of Log10(normalized counts) for top 10 DE genes",
      x = sample,
      y = "log10(norm_counts)"
    ) +
    theme_minimal()
  
  return(scatter_plot)
}

#' Function to generate volcano plot from DESeq2 results
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one
#'   additional column denoting status in volcano plot
#'
#' @return ggplot: a scatterplot (volcano plot) that displays log2foldchange vs
#'   -log10(padj) and labeled by status
#' @export
#'
#' @examples volcano_plot <- plot_volcano(labeled_results)
#' 
plot_volcano <- function(labeled_results) {
  volcano_plot <- ggplot(data = labeled_results, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(color = volc_plot_status)) +
    labs(
      title = "Volcano plot of DESeq2 differential expression results (vP0 vs. vAd)",
      x = "log2FoldChange", 
      y = "-log10(padj)") +
    theme_minimal()
  
  return(volcano_plot)
}

#' Function to generate a named vector ranked by log2FC descending
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one
#'   additional column denoting status in volcano plot
#' @param id2gene_path (str): Path to the file containing the mapping of
#' ensembl IDs to MGI symbols
#'
#' @return Named vector with gene symbols as names, and log2FoldChange as values
#' ranked in descending order
#' @export
#'
#' @examples rnk_list <- make_ranked_log2fc(labeled_results, 'data/id2gene.txt')

make_ranked_log2fc <- function(labeled_results, id2gene_path) {
  # Read id2gene.txt file
  id2gene <- read.table(id2gene_path, header = FALSE, sep = "\t")
  # Set column names in the table to genes and IDs
  colnames(id2gene) <- c("genes", "ID")
  
  # Use left join to merge the id column of the id2gene table with the labeled results tibble
  merged_results <- labeled_results %>%
    left_join(id2gene, by = c("genes" = "genes")) %>%
    arrange(desc(log2FoldChange))
  
  # Create a named vector with log2FoldChange values and corresponding IDs
  ranked_log2fc <- merged_results$log2FoldChange
  names(ranked_log2fc) <- merged_results$ID
  ranked_log2fc <- na.omit(ranked_log2fc)
  
  return(ranked_log2fc)
}

#' Function to run fgsea with arguments for min and max gene set size
#'
#' @param gmt_file_path (str): Path to the gene sets of interest in GMT format
#' @param rnk_list (named vector): Named vector generated previously with gene 
#' symbols and log2Fold Change values in descending order
#' @param min_size (int): Minimum number of genes in gene sets to be allowed
#' @param max_size (int): Maximum number of genes in gene sets to be allowed
#'
#' @return Tibble of results from running fgsea
#' @export
#'
#' @examples fgsea_results <- run_fgsea('data/m2.cp.v2023.1.Mm.symbols.gmt', rnk_list, 15, 500)
run_fgsea <- function(gmt_file_path, rnk_list, min_size, max_size) {
  gene_sets_of_interest <- fgsea::gmtPathways(gmt_file_path)
  fgsea_results <- fgsea(gene_sets_of_interest, rnk_list, minSize = min_size, maxSize = max_size)
  fgsea_results <- as_tibble(fgsea_results)
  return(fgsea_results)
}

#' Function to plot top ten positive NES and top ten negative NES pathways
#' in a barchart
#'
#' @param fgsea_results (tibble): the fgsea results in tibble format returned by
#'   the previous function
#' @param num_paths (int): the number of pathways for each direction (top or
#'   down) to include in the plot. Set this at 10.
#'
#' @return ggplot with a barchart showing the top twenty pathways ranked by positive
#' and negative NES
#' @export
#'
#' @examples fgsea_plot <- top_pathways(fgsea_results, 10)
top_pathways <- function(fgsea_results, num_paths){
  # select the top positive and negative pathways
  top_positive <- fgsea_results %>% 
    filter(NES > 0) %>%
    arrange(desc(NES)) %>% 
    mutate(pos_or_neg = "positive") %>% 
    head(num_paths)
  
  top_negative <- fgsea_results %>% 
    filter(NES < 0) %>%
    arrange(NES) %>% 
    mutate(pos_or_neg = "negative") %>% 
    head(num_paths)
  
  # combine the top positive and top negative pathways
  top_pathways_df <- dplyr::bind_rows(top_positive, top_negative)
  
  # create a ggplot barchart
  pathway_plot <- ggplot(top_pathways_df, mapping = aes(x = NES, y = reorder(pathway, NES), fill = pos_or_neg)) +
    geom_bar(stat = "identity") +
    labs(
      title = "GSEA results for Hallmark MSIgDB gene set",
      x = "Normalized Enrichment Score (NES)", ) +
    theme_minimal()
  
  return(pathway_plot)
}

