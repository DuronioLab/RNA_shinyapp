## List of all CRAN repos packages that are needed (installed and loaded)
CRANpackages = c("BiocManager", "ggplot2", "gridExtra", "RColorBrewer", "pheatmap",
                 "dplyr", "readr", "colorspace", "gplots", "DT", "heatmaply", "ggrepel", "Rcpp",
                 "stringr", "ggside", "stringr", "tidyr", "tibble", "shiny", "ggforce",
                 "ggVennDiagram")

## List of all Bioconductor repos packages that are needed (installed and loaded)
BioCpackages = c("DESeq2", "Gviz","GenomicRanges", "TxDb.Dmelanogaster.UCSC.dm6.ensGene",
                 "Rsamtools", "clusterProfiler", "org.Dm.eg.db")

CRANpackage.check <- lapply(
  CRANpackages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)
    }
  }
)

BioCpackage.check <- lapply(
  BioCpackages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      BiocManager::install(x)
      library(x, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)
    }
  }
)

## specify source file containing all function calls
source(file = "./extra_functions.R")

## Run the DESeq2 pipeline
source(file = "./diff_exp.R")

print("Libraries loaded and files sourced.")

#Make the individual log2FC tables
log2fc_samples <-comparisons
for (i in log2fc_samples) {
  i_data <- paste(i, "_resOrdered", sep='')
  i_data <- subset(eval(as.name(i_data)), select = c(log2FoldChange))
  i_data <- tibble::rownames_to_column(as.data.frame(i_data))
  i <- paste(i, "Log2FC", sep="_")
  colnames(i_data)[2] <- paste(i)
  assign(i, i_data)
}

## create table containing all logFC values
first_name <- paste(log2fc_samples[1], "_Log2FC", sep="")
all_logfc <-eval(as.symbol(first_name))
for (i in 2:length(log2fc_samples)) {
  i_name <- paste(log2fc_samples[i], "_Log2FC", sep="")
  #if(!exists("i_name")){
  i_data <- eval(as.symbol(i_name))
  all_logfc <- dplyr::left_join(all_logfc, i_data, by = "rowname")
  #}
}
all_logfc <- tibble::column_to_rownames(all_logfc, var = "rowname")
pre_all_logfc <- all_logfc
print("Prepared Log2 Fold Change tables")


## Possibly delete later (check if logf2fc heatmaps still works after deleting)



#if(!exists("geneSelector")){
  temp_selected <- dplyr::select(gtf_genes, gene_symbol, gene_id)
  all_logfc <- dplyr::left_join(tibble::rownames_to_column(all_logfc), temp_selected, by = c("rowname" = "gene_symbol"))
  colnames(all_logfc)[1] <- "gene_symbol"

  
  for(i in 1:nrow(all_logfc)){
    if(is.na(all_logfc$gene_id[i])){
      fixed <- get_geneID(all_logfc$gene_symbol[i], gtf, current_fb_ids)
      if(length(fixed) == 1){
        all_logfc$gene_id[i] <- fixed
      }
    }else{
      
    }
  }
  geneSelector <- TRUE
  rm(temp_selected)
#}

## Make the norm_counts
if(!exists("norm_counts")){
  vsd <- vst(dds, blind=TRUE)
  norm_counts <- assay(vsd)
}

## replace tmp with modified current_fb_ids???
# tmp2 <- data.frame(gene_symbol = current_fb_ids$current_symbol,
#                    gene_id = current_fb_ids$primary_FBid)
  
## Make annotated fbgn_counts vectors
if(!exists("fbgn_counts")){
  fbgn_counts <- dplyr::left_join(tibble::rownames_to_column(as.data.frame(norm_counts)), tmp, by = c("rowname" = "gene_symbol"))
  
  for(i in 1:nrow(fbgn_counts)){
    if(is.na(fbgn_counts$gene_id[i])){
      fixed <- get_geneID(fbgn_counts$rowname[i], gtf, current_fb_ids)
      if(length(fixed) == 1){
        fbgn_counts$gene_id[i] <- fixed
      }
    }else{
    }
  }
  rownames(fbgn_counts) <- fbgn_counts[,1]
  fbgn_counts <- fbgn_counts[,-1]
}

repdata2 <- data.frame()

for (i in 1:nrow(sample_table)){
  repdata2[i,1] <- sample_table$sample[i]
  repdata2[i,2] <- sample_table$rep[i]
}
repList2 <- c(unique(repdata2[,1]))

## Intron annotations for the gene browser

## Fix to work with any txdb for genome browser

# if(!exists("introns")){
#   txdb <- makeTxDbFromGFF(paste("./input_data/",annotation_file, sep = ""))
#   introns <- intronsByTranscript(txdb)
#   introns <- unlist(introns)
#   introns <- unique(introns)
# }
print("loaded txdb and intron functions")

# Ensuring that some outputs are saved
raw_res_table_out <- data.frame()
cluster_table_out <- data.frame()
fc_image_out <- list()
gene_exp_image_out <- list()
cluster_box_image_out <- list()
volcano_image_out <- list()
hm_data_image_out <- matrix()
hm_cluster_image_out <- data.frame()
hm_info_image_out <- data.frame()


print("loaded empty data containers for output")

print("finished with all global data")


