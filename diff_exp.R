## set the working directory to current folder
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## import the sample sheet output by snakemake pipeline
if(!exists("sample_table")){
  sample_table <- read.table("./input_data/sampleSheet.tsv", sep= "\t", header=TRUE)
}
print("Sample sheet imported.")

## add in the change to the featurecounts column name (.Aligned.sortedBy.....) Also change
# the ordering to match sample sheet and feature counts


## import count table and remove the first 5 columns that contain annotation information, leaving only the read counts
if(!exists("count_table")){
  count_table <- readr::read_tsv("./input_data/featureCounts.txt", comment="#")
  count_table <- count_table %>% tibble::column_to_rownames("Geneid")
  gene_counts <- as.matrix(count_table[,6:ncol(count_table)])
  colnames(gene_counts) <- sub("_dm6.*", "", colnames(gene_counts))
  
  ## generate a data frame with information about each sample
  coldata <- as.matrix(sample_table$sample)
  coldata <- data.frame(sample_table$sample)
  rownames(coldata) <- sample_table$bam
  colnames(coldata) <- c("sample")
}
print("Count data imported.")

##Warning: might not work! It should but I've been having issues with it dropping info
coldata <- coldata[ order(row.names(coldata)),, drop = FALSE]
gene_counts <- gene_counts[ , order(colnames(gene_counts)), drop = FALSE]


## Perform DESeq2 analysis and save if not already created
if(!file.exists("./output_data/dds.rds")){
  dds_data <- DESeqDataSetFromMatrix(gene_counts, coldata, design = ~ sample)
  
  ## filter out genes with less than 50 reads across all samples
  keep <- rowSums(counts(dds_data)) >= 50
  dds_data_filtered <- dds_data[keep,]
  
  ## test for differential gene expression
  dds2 <- DESeq(dds_data_filtered)
  dds <- DESeq(dds_data)
  
  ## write dds2 object to Rdata file
  saveRDS(dds2, file = "output_data/dds2.rds")
  saveRDS(dds, file = "output_data/dds.rds")
}else{
  dds2 <- readRDS("output_data/dds2.rds")
  dds <- readRDS("output_data/dds.rds")
}


print("DESeq2 data imported or generated")

##Get the sample names
sample_names <- as.matrix(unique(coldata))
sample_names <- sample_names[1:length(sample_names)]

if(!exists("choices")){
  ##Make pair-wise comparisons for all samples with a pvalue threshold of 0.05, also save a list of the comparisons
  comparisons <- c()
  compList <- c()
  
  i<-1
  while(i <= length(sample_names)){
    j<-1
    while(j <= length(sample_names)){ 
      #if(sample_names[i]!=sample_names[j] & !exists(paste(sample_names[j],"_vs_",sample_names[i],"_resOrdered",sep = ""))){
      
      ##Important: Checks for and only does one of two comparisons, unless uncommented out and the above is commented out.    
      if(sample_names[i]!=sample_names[j]){
        
        comp <-paste(sample_names[i], "_vs_",sample_names[j], sep = "")
        comparisons <- c(comparisons, comp)
        
        compName <- paste(sample_names[i], "vs", sample_names[j], sep = " ")
        compList <- c(compList, compName)
        
        temp_name <- paste(i, "_vs_",j,"_res", sep = "")
        temp_res <- results(dds2, contrast = c("sample", sample_names[i], sample_names[j]), alpha = 0.05)
        temp_name_ord <- paste(sample_names[i], "_vs_",sample_names[j],"_resOrdered", sep = "")
        temp_ordered <- temp_res[order(temp_res$pvalue),]
        assign(temp_name_ord, temp_ordered)
        
        fn <- paste("./output_data/", temp_name_ord, ".txt", sep = "")
        write.table(as.data.frame(temp_ordered) %>% tibble::rownames_to_column("gene_symbol"), file=fn, sep="\t", row.names = FALSE)
      }
      j <- j+1
    }
    i<-i+1
  }
  
  #Constructs the "choices" named vector list for the ui.R script
  choices <- comparisons
  names(choices) <-compList
}
print("Pair-wise comparisons generated.")


#### look at overlap of differentially expressed genes between different comparisons

all_diff_fbgns <- list()
all_diff_fbgns_index <- list()
j <- 1
for (i in comparisons) {
  ## for each comparison, create a new table containing only significantly differentially expressed genes
  i_name <- paste(i, "_resOrdered", sep = '')
  i_data <- eval(as.symbol(i_name))
  sig <- which_diff(i_data, p_threshold = 0.05, log2fc_threshold = 1)
  new_name <- paste(i_name, "_sig", sep = "")
  assign(new_name, sig)
  
  ## create a list of the gene IDs from the table above
  sig_names <- rownames(sig)
  listname <- paste(i_name, "_sig_fbgn", sep = "")
  assign(listname, sig_names)
  all_diff_fbgns_index[[j]] <- listname
  #all_diff_fbgns <- append(all_diff_fbgns, listname)
  all_diff_fbgns[[j]] <- sig_names
  j <- j+1
}

## create table summarizing overlap of DE genes between comparisons

all_diff_fbgns_names <- c()
for (i in 1:length(all_diff_fbgns[])){
  all_diff_fbgns_names <- c(all_diff_fbgns_names, all_diff_fbgns[[i]])
  #print(length(all_diff_fbgns_names))
}
all_diff_fbgns_names <- unique(all_diff_fbgns_names)
is_fbgn_diff <- data.frame(all_diff_fbgns_names, stringsAsFactors = FALSE)
colnames(is_fbgn_diff)[1] <- "gene_symbol"
fbgn_diff <- data.frame(matrix(, nrow=length(is_fbgn_diff[,]), ncol=length(comparisons)))
rownames(fbgn_diff) <- is_fbgn_diff[,]


for (i in 1:length(comparisons)) {
  i_name <- paste(comparisons[i], "_resOrdered_sig_fbgn", sep = '')
  i_values <- as.data.frame(eval(as.symbol(i_name)))
  logical_column <- is_fbgn_diff[,] %in% i_values[,]
  fbgn_diff[,i] <- logical_column
  colnames(fbgn_diff)[i] <- comparisons[i]
}
fbgn_diff <- tibble::rownames_to_column(fbgn_diff, var = "gene_symbol")

fbgn_diff_genes <- fbgn_diff
for(i in 1:nrow(fbgn_diff_genes)){
  for(j in 2:ncol(fbgn_diff_genes)){
    if(fbgn_diff_genes[i,j] == TRUE){
      fbgn_diff_genes[i,j] <- fbgn_diff_genes[i,1]
    }else if(fbgn_diff_genes[i,j] == FALSE){
      fbgn_diff_genes[i,j] <- NA
    }
  }
}

print("DE gene tables summarized.")


#### write a results table of significant genes with additional information about each gene
## import dm6 GTF file and add columns for gene_symbol and gene_ID

if(!exists("gtf_genes")){
  gtf_fn <- "./input_data/dmel-all-r6.20.ucsc.gtf"
  gtf <- readGTF(gtf_fn)
  gtf$gene_symbol <- GTF.attributes(gtf, "gene_symbol")
  gtf$gene_id <- GTF.attributes(gtf, "gene_id")
  gtf_genes <- dplyr::filter(gtf, feature=="gene")
}
print("GTF object loaded")

## import flybase file containing current and previous gene IDs and gene symbols. This file will be used to update outdated gene names from published data
## Check to see if the fb_synonym_fb_2021_03.tsv file is up to date and/or named correctly!!

if(!exists("current_fb_ids")){
  #current_fb_ids <- read.table("./input_data/fb_synonym_fb_2021_03.tsv", sep="\t", as.is = TRUE,header = TRUE, quote = "", comment.char = "", skip = 5)
  suppressWarnings(current_fb_ids <- readr::read_tsv("./input_data/fb_synonym_fb_2021_03.tsv", skip = 5, quote = "", col_names = TRUE))
  current_fb_ids <- current_fb_ids %>% rename(primary_FBid = "##primary_FBid")
  current_fb_ids <- current_fb_ids %>% rename(fullname_synonym.s. = 'fullname_synonym(s)')
  current_fb_ids <- current_fb_ids %>% rename(symbol_synonym.s. = 'symbol_synonym(s)')
  
  current_fb_ids <- dplyr::filter(current_fb_ids, organism_abbreviation == "Dmel")
  current_fb_ids <- dplyr::filter(current_fb_ids, grepl('FBgn', current_fb_ids$primary_FBid))
}
print("Flybase Gene IDs loaded")

for (i in comparisons) {
  ## for each comparison, create a new table containing only significantly differentially expressed genes
  i_name <- paste(i, "_resOrdered_sig", sep = '')
  
  #Only make the new file if it doesn't already exist, also runs once if tmp is missing
  checkFile <- paste("./output_data/", i_name, "_final.txt",sep = "")
  if(!file.exists(checkFile) | !exists("tmp")){
    i_data <- eval(as.symbol(i_name))
    i_data <- as.data.frame(i_data)
    i_data <- tibble::rownames_to_column(i_data, var = "gene_symbol")
    
    ## add new column to results containing the FBgn gene id
    tmp <- dplyr::select(gtf_genes, gene_symbol, gene_id)
    new_table <- dplyr::left_join(i_data, tmp, by = "gene_symbol")
    
    ## add columns to results table indicating whether each gene was also differentially expressed in another comparison
    new_table <- dplyr::left_join(new_table, fbgn_diff, by = "gene_symbol")
    
    ## create a final output table, save it to an object, and write it to a file
    new_name <- paste(i_name, "_final.txt", sep = "")
    assign(new_name, new_table)
    fn <- paste("./output_data/", new_name, sep = "")
    write.table(new_table, file=fn, sep="\t", row.names = FALSE)
    rm(new_table)
  }
}


print("Significance tables saved to ./output_data/ directory.")

## Extra functions to grab data for the server.R file

bamFileList <- c()
for (i in sample_table$bam){
  bamName <- strsplit(i, "/")[[1]][2]
  bamName <- strsplit(bamName, "\\.")
  bamName <- paste("./input_data/", bamName[[1]][1], "_sorted.bam", sep = "")
  bamFileList <- c(bamFileList, bamName)
}

for (i in sample_table$baseName){
  nameIndex <- grep(i, bamFileList)
  temp_name <- paste(i, "_bam", sep = "")
  assign(temp_name, bamFileList[nameIndex])
}
