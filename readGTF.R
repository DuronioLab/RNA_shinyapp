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


gtf_fn <- "dmel-all-r6.20.ucsc.gtf"
gtf <- readGTF(gtf_fn)
gtf$gene_symbol <- GTF.attributes(gtf, "gene_symbol")
gtf$gene_id <- GTF.attributes(gtf, "gene_id")
