## this file contains all the function definitions for the RNAseq analysis

## Perform heirarchical clustering
perform_hc <- function(what = "hc", cluster_n = 3, input_counts = NA){
  
  # If no specific counts have been input, just take the norm_counts
  suppressWarnings(
    if(is.na(input_counts)){
      norm_counts_clust <- norm_counts
      norm_counts_clust <- as.data.frame(norm_counts_clust)
      cluster_data <- norm_counts_clust
    }else{
      norm_counts_clust <- as.data.frame(input_counts)
      cluster_data <- norm_counts_clust
    }
  )
  ## compute distance matrix
  sampleDists <- dist(cluster_data)
  
  ## perform hierarchecal clustering
  hc <- hclust(sampleDists)
  
  ## extract cluster assignments for each gene
  hc.cutree <- cutree(hc, cluster_n)
  clusters <- data.frame(clusters = as.factor(hc.cutree))
  clusters <<- tibble::rownames_to_column(clusters, var = "gene_id")
  cluster_data$clusters <- clusters$clusters
  
  if(what == "hc"){
    return(hc)
  }else{
    return(cluster_data)
  }
}


## Perform k-means clustering
perform_kc <- function(what = "kc", cluster_n = 3, input_counts = NA){
  
  # If no specific counts have been input, just take the norm_counts
  suppressWarnings(
    if(is.na(input_counts)){
      norm_counts_clust <- norm_counts
      norm_counts_clust <- as.data.frame(norm_counts_clust)
      cluster_data <- norm_counts_clust
    }else{
      norm_counts_clust <- as.data.frame(input_counts)
      cluster_data <- norm_counts_clust
    }
  )
  
  set.seed(20)
  kc <- kmeans(cluster_data, centers=cluster_n, nstart = 1000, iter.max = 20)
  kClusters <- as.factor(x = kc$cluster)
  clusters <- data.frame(clusters = kClusters)
  clusters <- tibble::rownames_to_column(clusters, var = "gene_id")
  cluster_data$clusters <- clusters$clusters
  if(what == "kc"){
    return(kc)
  }else{
    return(cluster_data)
  }
  
}



## define function to get list of genes differentially expressed in any comparison. x is a DEseq results table
which_diff <- function(x, p_threshold = 0.05, log2fc_threshold = 1) {
  meets_pvalue <- subset(x, padj < p_threshold)
  sig_genes <- subset(meets_pvalue, abs(log2FoldChange) > log2fc_threshold)
  return(sig_genes)
}

## define function to read GTF file into a table
readGTF <- function(gtf_file) {
  gtf_annotation <- read.table(gtf_file, sep = "\t", header = FALSE, as.is = TRUE)
  colnames(gtf_annotation) <- c("seqname","source","feature","start","end","score","strand","frame","attribute")
  return(gtf_annotation)
}

## define function to extract a particular attribute for all genes. x is a GTF file, attr is the attribute you want to extract (e.g. "gene_symbol")
GTF.attributes <- function(x, attr, attrsep = ';') {
  atts <- strsplit(x$attribute, split = attrsep, fixed = TRUE)
  get_attr <- function(x, attr) {
    attr_pair <- x[grep(attr, x)]
    if (length(attr_pair) > 0){
      attr_value <- unlist(strsplit(attr_pair, split = '\\s'))
      attr_value <- attr_value[attr_value != ""]
      return(attr_value[2])
    } else {
      return('NA')
    }
  }
  all_attrs <- sapply(atts, get_attr, attr)
  return(all_attrs)
  
  
}

## define function to take a data frame containing flybase gene names and add a column with the most up-to-date flybase gene ID (FBgn)
## this function takes as input a synonym table which can be downloaded from flybase. Only the "Primary_FBid", and "symbol_synonym" columns of the synonym table are required
## function arguments: 
##x is a data frame, gene_id_column
## gtf is a gtf file used to retrieve IDs for genes where the symbol is not out of date
## gene_column specifies the name of the column of x that contains the gene name
## synonym table is a table containing synonymous gene names/IDs. this can be retrieved from flybase

update_geneID <- function(x,gene_column, gtf, synonym_table) {
  ## print the name of the object the function is updating
  nm<- deparse(substitute(x))
  print(paste("update_geneID input object:", nm))
  
  ## get a gtf table to add gene_ids for gene names that haven't changed from dm3 to dm6
  gtf$gene_symbol <- GTF.attributes(gtf, "gene_symbol")
  gtf$gene_id <- GTF.attributes(gtf, "gene_id")
  gtf_genes <- dplyr::filter(gtf, feature=="gene")
  
  ## first check if the gene_symbols are in the gtf file, if they are, add a column containing the flybase gene_id
  tmp <- select(gtf_genes, gene_symbol, gene_id)
  name_cols <- c("gene_symbol")
  names(name_cols) <- gene_column
  x <- left_join(x, tmp, by=name_cols)
  
  ## for genes that were not found in the gtf table, check the flybase synonym table. only add gene_ids for genes with a single perfect match in the fb synonym table
  for (i in 1:nrow(x)) {
    if (is.na(x[i,"gene_id"])) {
      old_gene_name <- x[i,gene_column]
      synonym_index <- grep(old_gene_name, synonym_table$symbol_synonym.s., fixed = TRUE)
      exact_hits <- 0    
      
      for (j in synonym_index) {
        synonyms <- strsplit(synonym_table[j,"symbol_synonym.s."], ",")[[1]]
        if (old_gene_name %in% synonyms) {
          current_id <- synonym_table[j,"primary_FBid"]
          exact_hits <- exact_hits + 1
        }
      }
      if (exact_hits == 1) {
        x[i, "gene_id"] <- current_id
      } else if (exact_hits > 1) {
        x[i, "gene_id"] <- "ambiguous"
      } else if (exact_hits == 0) {
        x[i, "gene_id"] <- NA
      }
      
    }
  }
  print(paste("number of gene_ids not assigned by update_geneID:",length(which(is.na(x$gene_id))) + length(which(x$gene_id == "ambiguous"))))
  return(x)
  
}


## function to calculate R-squared values
rsq <- function(x,y) {
  rsq <- cor(x,y) ^ 2
  rounded <- round(rsq, 4)
  return(rounded)
}


## function to remove all columns from a table that have the same gene_id. x is a data frame and gene_column is the name of the column containing the gene_id
rm_dup_genes <- function(x, gene_column) {
  dup_rows <- vector("list", nrow(x))
  for (i in 1:nrow(x)) {
    dups <- which(x[,gene_column] %in% x[i,gene_column])
    if (length(dups) > 1) {
      dup_rows[[i]] <- dups
    }
  }
  dup_rows <- unlist(dup_rows)
  deduped_x <- x[-dup_rows,]
  in_rows <- nrow(x)
  out_rows <- nrow(deduped_x)
  dropped_rows <- in_rows - out_rows
  print(paste(dropped_rows,"duplicate rows were removed from the table, leaving", out_rows))
  return(deduped_x)
}


## function to compute the coefficient of variation of a vector or matrix
cv <- function(x) {
  sd <- sd(x)
  mean <- mean(x)
  cv <- sd/mean * 100
  return(cv)
}

## function to perform fisher's exact test for enrichment with a group of differentially expressed genes
## arguments are as follows:
# x: a results table produced by diff_exp
# test_group: the set of genes to use for enrichment testing. Has three possible values: up, down, or diff
# is_diff_column: a logical column of x specifying if the gene is differentially expressed
# fc_column: a column of x containg a foldchange value
# test_column: a column of x indicating if genes are a member of the group being tested for enrichment
# test_criteria: the values in test_column indicating group membership
# log2fc_threshold: the logfc threshold that was used to determine differential expression. default is 1
fisher_diff_genes <- function(x, test_group, is_diff_column, fc_column="log2FoldChange", test_column, test_criterion, log2fc_threshold = 1) {
  require(dplyr)
  possible_test_groups <- c("up", "down", "diff")
  
  if (test_group == "up") {
    up_data <- dplyr::filter_(x, paste(is_diff_column, "==", TRUE)) %>% dplyr::filter_(paste(fc_column,  ">", 0))
    not_up_data <- dplyr::filter_(x, paste(fc_column,  "<", log2fc_threshold))
    
    up_in <- length(which(up_data[,test_column] == test_criterion))
    up_not_in <- length(which(up_data[,test_column] != test_criterion))
    not_up_in <- length(which(not_up_data[,test_column] == test_criterion))
    not_up_not_in <- length(which(not_up_data[,test_column] != test_criterion))
    
    up_table <- data.frame(up = c(up_in, up_not_in), not_up = c(not_up_in, not_up_not_in), row.names = c("in", "not_in"))
    up_ft <- fisher.test(up_table)
    up_pvalue <- up_ft$p.value
    perc_up <- up_in/sum(up_table$up)*100
    results <- list(pvalue = up_pvalue, percent = perc_up, contingency_table = up_table)
    print(paste("fishers exact test for enrichment:"))
    print(paste("gene set tested: upregulated genes"))
    print(paste("tested for enrichment based on the following criterion:", test_column, "==", test_criterion))
    return(results)
    
  } else if (test_group == "down") {
    down_data <- dplyr::filter_(x, paste(is_diff_column, "==", TRUE)) %>% dplyr::filter_(paste(fc_column,  "<", 0))
    not_down_data <- dplyr::filter_(x, paste(fc_column,  ">", -log2fc_threshold))
    
    down_in <- length(which(down_data[,test_column] == test_criterion))
    down_not_in <- length(which(down_data[,test_column] != test_criterion))
    not_down_in <- length(which(not_down_data[,test_column] == test_criterion))
    not_down_not_in <- length(which(not_down_data[,test_column] != test_criterion))
    
    down_table <- data.frame(down = c(down_in, down_not_in), not_down = c(not_down_in, not_down_not_in), row.names = c("in", "not_in"))
    down_ft <- fisher.test(down_table)
    down_pvalue <- down_ft$p.value
    perc_down <- down_in/sum(down_table$down)*100
    results <- list(pvalue = down_pvalue, percent = perc_down, contingency_table = down_table)
    print(paste("fishers exact test for enrichment:"))
    print(paste("gene set tested: downregulated genes"))
    print(paste("tested for enrichment based on the following criterion:", test_column, "==", test_criterion))
    return(results)
    
  } else if (test_group == "diff") {
    diff_data <- dplyr::filter_(x, paste(is_diff_column, "==", TRUE)) 
    not_diff_data <- dplyr::filter_(x, paste(is_diff_column, "==", FALSE))
    
    diff_in <- length(which(diff_data[,test_column] == test_criterion))
    diff_not_in <- length(which(diff_data[,test_column] != test_criterion))
    not_diff_in <- length(which(not_diff_data[,test_column] == test_criterion))
    not_diff_not_in <- length(which(not_diff_data[,test_column] != test_criterion))
    
    diff_table <- data.frame(diff = c(diff_in, diff_not_in), not_diff = c(not_diff_in, not_diff_not_in), row.names = c("in", "not_in"))
    diff_ft <- fisher.test(diff_table)
    diff_pvalue <- diff_ft$p.value
    perc_diff <- diff_in/sum(diff_table$diff)*100
    results <- list(pvalue = diff_pvalue, percent = perc_diff, contingency_table = diff_table)
    
    print(paste("fishers exact test for enrichment:"))
    print(paste("gene set tested: differentially expressed genes"))
    print(paste("tested for enrichment based on the following criterion:", test_column, "==", test_criterion))
    return(results)
    
  } else if (!(test_group %in% possible_test_groups)) {
    stop("invalid test_group argument: should be set to 'up', 'down', or 'diff' ")
  }
}

## convert log2 fc values to linear space
un_log2 <- function(x) {
  if (x > 0) {
    new <- 2^x
    return(new)
  } else if (x < 0) {
    new <- -2^abs(x)
    return(new)
  }
}

## function to filter a results table generated by diff_exp.R to exclude CRY-effect genes and/or include only zld_targets
## nested function to subset a data frame based on whether genes are have CRY-specific expression or if they are zld targets
filter_genes <- function(x, filter_CRY = FALSE, filter_zld = FALSE) {
  if (filter_CRY == TRUE) {
    x <- dplyr::filter(x, is_cry_gene == FALSE)
  }
  
  if (filter_zld == TRUE) {
    x <- dplyr::filter(x, is_zld_target == TRUE)
  }
  return(x)
}

## function to import DEseq results tables into R, subset the data based on provided criteria, and provide a list of gene names that can be used to subset any data to a specific set of genes
filter_genes2 <- function(samples_list, diff_genes_only = "No", change_direction="none") {
  require(dplyr)
  if (!(change_direction %in% c("up", "down", 'none'))) {
    stop("invalid value for 'change direction'. Acceptable values are 'up', 'down', or 'none'")
  }
  
  if (diff_genes_only == "Yes") {
    for (i in samples_list) {
      i_name <- paste(i, "_resOrdered_sig_final.txt", sep = '')
      fn <- paste("./output_data/", i_name, sep = "")
      i_data <- read.table(file = fn, sep = "\t", header = TRUE, as.is = TRUE)
      
      if (change_direction == "up") {
        i_data <- filter(i_data, log2FoldChange > 0)
      }
      
      if (change_direction == "down") {
        i_data <- filter(i_data, log2FoldChange < 0)
      }
      assign(i, i_data)
    }
    
    ## get list of all genes meeting the filtering criteria
    keep_genes <- vector("list", length(samples_list))
    for (i in 1:length(samples_list)) {
      i_data <- eval(as.symbol(samples_list[i]))
      keep_genes[[i]] <- i_data$gene_symbol
    }
    
    keep_genes <- unique(unlist(keep_genes))
    return(keep_genes)
    
  } else if (diff_genes_only == "No") {
    for (i in samples_list) {
      i_name <- paste(i, "_resOrdered.txt", sep = '')
      fn <- paste("./output_data/", i_name, sep = "")
      i_data <- read.table(file = fn, sep = "\t", header = TRUE, as.is = TRUE)
      assign(i, i_data)
    }
    
    ## get list of all genes meeting the filtering criteria
    keep_genes <- vector("list", length(samples_list))
    
    for (i in 1:length(samples_list)) {
      i_data <- eval(as.symbol(samples_list[i]))
      keep_genes[[i]] <- i_data$gene_symbol
    }
    keep_genes <- unique(unlist(keep_genes))
    return(keep_genes)
    
  }
  
  
}

## define function to take a single gene name and return the most up-to-date flybase gene ID (FBgn)
## this function takes as input a synonym table which can be downloaded from flybase. Only the "Primary_FBid", and "symbol_synonym" columns of the synonym table are required
## function arguments: 
## x is a gene name
## gtf is a gtf file used to retrieve IDs for genes where the symbol is not out of date
## synonym table is a table containing synonymous gene names/IDs. this can be retrieved from flybase

get_geneID <- function(x, gtf, synonym_table) {
  
  synonym_table <- as.data.frame(synonym_table)
  
  #Check if the gtf and gtf_genes have been made/modified. If not, then do so
  if(!exists("gtf_genes")){
    gtf$gene_symbol <- GTF.attributes(gtf, "gene_symbol")
    gtf$gene_id <- GTF.attributes(gtf, "gene_id")
    gtf_genes <- dplyr::filter(gtf, feature=="gene")
  }
  if(grepl("FBgn", x)){
    #print("You submitted an FBgn")
    return(x)
  }
  
  #First pass check if the input gene name is in the gtf table
  tmp <- dplyr::select(gtf_genes, gene_symbol, gene_id)
  fbgn_id <- tmp[grep(x, tmp$gene_symbol), ]$gene_id
  
  #Check synonym and fullname tables for the input gene
  synonyms <- vector()
  if (length(fbgn_id) == 0) {
    old_gene_name <- x
    #print(x)
    synonym_index <- grep(old_gene_name, synonym_table$symbol_synonym.s., fixed = TRUE)
    #print(synonym_index)
    if(length(synonym_index) == 0){
      synonym_index <- grep(old_gene_name, synonym_table$fullname_synonym.s., fixed = TRUE)
      #print(synonym_index)
    }
    for (j in synonym_index) {
      synonyms <- strsplit(synonym_table[j,"symbol_synonym.s."], ",")[[1]]
      #print(synonyms)
      if (old_gene_name %in% synonyms) {
        fbgn_id <- synonym_table[j,"primary_FBid"]
        #print(fbgn_id)
      }
    }
    if(length(synonyms) == 0){
      for (j in synonym_index) {
        #print("problem area?")
        synonyms <- strsplit(synonym_table[j,"fullname_synonym.s."], ",")[[1]]
        if (old_gene_name %in% synonyms) {
          fbgn_id <- synonym_table[j,"primary_FBid"]
        }
      }
    }
    if (length(synonym_index >= 0)) {
      sss <- NA
    }
  }
  #print(paste("I found:", fbgn_id, "as the ID for the gene", tmp[grep(fbgn_id, tmp$gene_id), ]$gene_symbol))
  return(fbgn_id)
  
}