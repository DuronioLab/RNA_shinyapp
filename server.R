

server <- function(input, output, session){
  
  #outputOptions(output, "color_check", suspendWhenHidden = FALSE)
  #observe(input$v_gene_size)
  ## PCA and Scree plot
  PCA_data <- eventReactive(input$PCA_go, {input$PCA_which})
  pca_list_data <- eventReactive(input$PCA_go,
                                 if(length(input$pca_list) == 0){
                                   return(input$pca_list)
                                 }else{
                                   return(unlist(strsplit(input$pca_list, split=", ")))
                                 }
  )
  #saved <<- observe({print(input$pca_list)})
  
  ##Venn diagram reactives
  
  venn_sample_list <- eventReactive(input$venn_go, {unlist(strsplit(input$venn_list, split=", "))})
  
  # venn_sample_list <- eventReactive(input$venn_go, {
  #   validate(
  #     need(if(length(observe({print(input$venn_list)})) < 2){TRUE}, message="You need more than 1 sample!"),
  #     need(if(length(observe({print(input$venn_list)})) > 4){TRUE}, message="You need less then 5 samples!")
  #     )
  #   input$venn_list
  # })
  #venny <<- observe({print(input$venn_list, split=", ")})
  #observe({print(input$venn_list)})
  
  
  ## Gene Ontology reactives
  GO_sample_data <- eventReactive(input$GO_go, {unlist(strsplit(input$GO_sample_list, split=", "))})
  GO_direction <- eventReactive(input$GO_go, {input$direction_go})
  GO_all_bg <- eventReactive(input$GO_go, {input$go_all})
  
  GO_bg_data <- eventReactive(input$GO_go, {unlist(strsplit(input$GO_bg_list, split=", "))})
  GO_bg_direction <- eventReactive(input$GO_go, {input$bg_go})
  
  ## Boxplot cluster reactives
  plot_which_data <- eventReactive(input$cluster_go, {unlist(strsplit(input$plot_which_list, split=", "))})
  box_order_data <- eventReactive(input$cluster_go, {unlist(strsplit(input$box_order_list, split=", "))})
  
  ## Log2 FC heatmap reactives
  fc_heat_comp_data <- eventReactive(input$fc_heat_go, {input$fc_heat_comp})
  fc_heat_sample_sort_data <- eventReactive(input$fc_heat_go, {input$fc_heat_sample_sort})
  fc_cluster_n_data <- eventReactive(input$fc_heat_go, {input$fc_cluster_n})
  fc_heat_method_data <- eventReactive(input$fc_heat_go, {input$fc_heat_method})
  fc_sort_sample_list_data <- eventReactive(input$fc_heat_go, {
    validate(
      need(if(input$fc_heat_comp != input$fc_sort_sample_list){TRUE}, message="Cannot sort by the control!")
    )
    input$fc_sort_sample_list
  })
  fc_heat_diff_only_data <- eventReactive(input$fc_heat_go, {input$fc_heat_diff_only})
  #fc_heat_sample_list_data <- eventReactive(input$fc_heat_go,{unlist(strsplit(input$fc_heat_sample_list, split=", "))})
  fc_heat_sample_list_data <- eventReactive(input$fc_heat_go, {
    validate(
      need(if(length(unlist(strsplit(input$fc_heat_sample_list, split=", "))) < 2 | !any(input$fc_heat_sample_list %in% grep("All", input$fc_heat_sample_list, value = T))){TRUE}, message="Please select two or more comparisons to plot!")
    )
    unlist(strsplit(input$fc_heat_sample_list, split=", "))
  })
  
  #<- eventReactive(input$fc_heat_go, {unlist(strsplit(input$fc_heat_sample_list, split=", "))})
  
  
  ## Clustering reactives
  cluster_method_data <- eventReactive(input$cluster_go, {input$cluster_method})
  # cluster_by_data <- eventReactive(input$cluster_go, {input$cluster_by})
  cluster_n_data <- eventReactive(input$cluster_go, {input$cluster_n})
  direction_of_change_data <- eventReactive(input$cluster_go, {input$direction_of_change})
  direction_bg_change_data <- eventReactive(input$cluster_go, {input$direction_bg_change})
  cluster_sample_data <- eventReactive(input$cluster_go, {unlist(strsplit(input$cluster_sample_list, split=", "))})
  diff_only_data <- eventReactive(input$cluster_go, {input$diff_only})
  
  #Specifically for GO analysis of clusters
  cluster_GO_data <- eventReactive(input$cluster_go, {unlist(strsplit(input$cluster_GO_list, split = ", "))})
  cluster_bg_data <- eventReactive(input$cluster_go, {unlist(strsplit(input$cluster_bg_list, split = ", "))})
  
  ## Count data graph reactives
  count_data <- eventReactive(input$count_go, {unlist(strsplit(input$count_gene_list, split=", "))})
  lc_data <- eventReactive(input$count_go, {input$line_connect})
  sample_data <- eventReactive(input$count_go, {unlist(strsplit(input$count_sample_list, split=", "))})
  
  ## Fold change graph reactives
  fc_data <- eventReactive(input$fc_go, {
    validate(
      need(input$fc_gene_list, "Please input one or more genes to plot!")
    )
    unlist(strsplit(input$fc_gene_list, split=", "))
  })
  fc_sample_data <- eventReactive(input$fc_go, {
    validate(
      need(input$fc_sample_list, "Please select one or more samples to plot!")
    )
    unlist(strsplit(input$fc_sample_list, split=", "))
  })
  
  ## Results table reactives
  first_data <- eventReactive(input$table_go, {input$first_group})
  second_data <- eventReactive(input$table_go, {input$second_group})
  
  ## Volcano plot reactives
  fishers_data <- eventReactive(input$volcano_go, {input$fishers})
  volcano_samples_data <-eventReactive(input$volcano_go, {
    validate(
      need(input$volcano_sample_list, "Please select samples to plot!")
    )
    unlist(strsplit(input$volcano_sample_list, split=", "))
  })
  volcano_padj <- eventReactive(input$volcano_go, {input$padj_n})
  volcano_fc <- eventReactive(input$volcano_go, {input$fc_n})
  volcano_x <- eventReactive(input$volcano_go, {input$top_n})
  volcano_colors <- eventReactive(input$volcano_go, {input$vc_color_list})
  volcano_font_size <- eventReactive(input$volcano_go, {input$v_gene_list})
  volcano_point_size <- eventReactive(input$volcano_go, {input$v_point_size})
  volcano_options <- eventReactive(input$volcano_go, {input$v_checks})
  
  ## MA plot reactives
  
  ma_samples_data <-eventReactive(input$ma_go, {
    validate(
      need(input$ma_sample_list, "Please select samples to plot!")
    )
    unlist(strsplit(input$ma_sample_list, split=", "))
  })
  ma_padj <- eventReactive(input$ma_go, {input$ma_padj_n})
  ma_fc <- eventReactive(input$ma_go, {input$ma_fc_n})
  ma_x <- eventReactive(input$ma_go, {input$ma_top_n})
  ma_colors <- eventReactive(input$ma_go, {input$ma_color_list})
  ma_gene_list <- eventReactive(input$ma_go, {input$ma_gene_list})
  ma_point_size <- eventReactive(input$ma_go, {input$ma_point_size})
  ma_options <- eventReactive(input$ma_go, {input$ma_checks})
  
  
  ## Browser plot reactives
  browser_gene_list_data <- eventReactive(input$browser_go, {
    validate(
      need(input$browser_gene_list, "Please input a single gene to plot!"),
      need(!grepl(',',input$browser_gene_list), "Please input only a single gene!")
    )
    input$browser_gene_list
  })
  track_list_data <-  eventReactive(input$browser_go, {
    validate(
      need(input$track_list, "Please select the tracks you wish to plot.")
    )
    input$track_list
  })
  plusStart_data <-  eventReactive(input$browser_go, {
    validate(
      need(input$plusStart, "Please set 'Extra Upstream' to a number (can be negative)")
    )
    input$plusStart
  })
  plusEnd_data <-  eventReactive(input$browser_go, {
    validate(
      need(input$plusEnd, "Please set 'Extra Downstream' to a number (can be negative)")
    )
    input$plusEnd
  })
  sashimi_data <- eventReactive(input$browser_go, {input$sashimi})
  track_type_data <- eventReactive(input$browser_go, {input$track_type})
  max_height_data <- eventReactive(input$browser_go, {
    validate(
      need(input$max_height, "Please set 'Maximum Height' to a number (including zero)!")
    )
    input$max_height
  })
  
  color_list_data <-eventReactive(input$browser_go, {
    if(nchar(input$color_list) > 0){
      unlist(strsplit(input$color_list, split=", "))
    }else{
      input$color_list
    }
  })
  
  #Ensure that there is a color palette (easily over-ridden within functions)
  paletteLength <- 100
  myColor <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(paletteLength)
  
  ## Test code that outputs the width (dimension[1]) and height (dimension[2]) of the user window
  output$dimension_display <- renderText({
    paste(input$dimension[1], input$dimension[2], input$dimension[2]/input$dimension[1])
  })
  
  get_textSize <- function(plots){
    
    x <- (input$dimension[2]/150)
    x <- x/length(plots)
    return(x)
  }
  
  output$colorText <- renderText({
    x <- c()
    if(!length(color_list_data() == 0)){
      for( i in 1:length(color_list_data())){
        test_color <- tolower(color_list_data()[i])
        if(!test_color %in% colors()){
          x <- append(x, color_list_data()[i])
        }
      }
      if(nchar(x) >= 1){
        paste(c("Could not find the color(s) ", x), collapse = " ")
      }
    }
  })
  
  ## Gene Ontology tables
  output$gene_ont <- DT::renderDataTable({
    withProgress(message = "Generating Gene Ontology Report", value = 0, {
      
      direction <- ""
      if(GO_direction() == "Up"){
        direction <- "up"
      }
      if(GO_direction() == "Down"){
        direction <- "down"
      }
      if(GO_direction() == "Both"){
        direction <- "none"
      }
      
      bg_direction <- ""
      if(GO_bg_direction() == "Up"){
        bg_direction <- "up"
      }
      if(GO_bg_direction() == "Down"){
        bg_direction <- "down"
      }
      if(GO_bg_direction() == "Both"){
        bg_direction <- "none"
      }
      
      incProgress(0.25, detail = "Importing gene sets")
      
      #Grab the filtered gene list
      GO_gene_list <<- filter_genes2(samples_list = GO_sample_data(),
                                     change_direction = direction,
                                     diff_genes_only = "Yes")
      # If the user selects "false" use specific gene set as bg, otherwise use all genes
      if(GO_all_bg() == FALSE){
        bg_gene_list <<- filter_genes2(samples_list = GO_bg_data(),
                                       change_direction = bg_direction,
                                       diff_genes_only = "Yes")
        
        incProgress(0.25, detail = "Performing GO analysis")
        ego <- enrichGO(gene = GO_gene_list,
                        universe = bg_gene_list,
                        keyType = "SYMBOL",
                        OrgDb = org.Dm.eg.db,
                        ont = "CC",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.01,
                        qvalueCutoff = 0.05
        )
        
      }else{
        incProgress(0.25, detail = "Performing GO analysis")
        ego <- enrichGO(gene = GO_gene_list,
                        keyType = "SYMBOL",
                        OrgDb = org.Dm.eg.db,
                        ont = "CC",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.01,
                        qvalueCutoff = 0.05
        )
      }
      ego_df <- as.data.frame(ego)
      if(nrow(ego) == 0){
        ego <- data.frame("V1" = "No enrichment")
      }
      DT::datatable(ego_df)
      
    })
    
  })
  
  ## Log2 Fold change heatmaps
  output$fc_heatmap <- renderPlot({
    withProgress(message = "Making Log2FC heatmap", value = 0, {
      # 
      # fc_heat_comp_data <- eventReactive(input$fc_heat_go, {input$fc_heat_comp})
      # fc_heat_sample_sort_data <- eventReactive(input$fc_heat_go, {input$fc_heat_sample_sort})
      # fc_cluster_n_data <- eventReactive(input$fc_heat_go, {input$fc_cluster_n})
      # fc_heat_method_data <- eventReactive(input$fc_heat_go, {input$fc_heat_method})
      # fc_sort_sample_list_data <- eventReactive(input$fc_heat_go, {input$fc_sort_sample_list})
      # fc_heat_diff_only_data <- eventReactive(input$fc_heat_go, {input$fc_heat_diff_only})
      # fc_heat_sample_list_data <- eventReactive(input$fc_heat_go, {unlist(strsplit(input$fc_heat_sample_list, split=", "))})
      
      controlComp <<- paste("vs_",fc_heat_comp_data(), sep = "")
      
      plotThese <<- fc_heat_sample_list_data()
      #plotThese <- plotThese[! plotThese %in% "All"]
      #plotThese <<- plotThese[! plotThese %in% grep(controlComp, plotThese, value = T)]
      #if(any(plotThese %in% grep(controlComp, plotThese, value = T))){print("found one")}
      
      #Select comparisons for specific control
      selected_L2FC_temp <<- select(tibble::column_to_rownames(all_logfc, var = "gene_symbol"), contains(controlComp), )
      if(!any(plotThese %in% grep("All", plotThese, value = T))){
        selected_L2FC <- select(selected_L2FC_temp, contains(plotThese), )
        if(ncol(selected_L2FC) < 2){
          selected_L2FC <- selected_L2FC_temp
        }
      }
      #global_sL2FC0 <<- selected_L2FC
      
      #Keep only significant genes
      selected_L2FC_samples <- colnames(selected_L2FC)
      selected_L2FC_samples <- sub("_Log2FC*", "", x = selected_L2FC_samples)
      keep_genes <- filter_genes2(samples = selected_L2FC_samples, 
                                  diff_genes_only = fc_heat_diff_only_data(),
                                  change_direction = "none")
      
      selected_L2FC <- selected_L2FC[rownames(selected_L2FC) %in% keep_genes,]
      #global_sL2FC1 <<- selected_L2FC
      
      incProgress(0.2, detail = "Performing clustering")
      #Run clustering
      if(fc_heat_method_data() == "Heirarchical"){
        cluster_data <- perform_hc(input_counts=selected_L2FC, cluster_n=fc_cluster_n_data(), what="anything")
      }else{
        cluster_data <- perform_kc(input_counts=selected_L2FC, cluster_n=fc_cluster_n_data(), what="anything")
      }
      
      #Shorten column names
      controlComp_cols <- colnames(selected_L2FC)
      controlComp_cols_short <- vector()
      for (i in controlComp_cols){
        controlComp_cols_short <- append(controlComp_cols_short,str_split(i,"_vs_")[[1]][1])
      }
      colnames(selected_L2FC) <- controlComp_cols_short
      #global_sL2FC2 <<- selected_L2FC
      
      #Pivot logFC data frames for use in ggplot
      cluster_data <- tibble::rownames_to_column(cluster_data,var = "geneid")
      cluster_data_long <- tidyr::pivot_longer(cluster_data,names_to = "comparison",values_to = "log2FC",
                                               cols = all_of(controlComp_cols))
      selected_L2FC <- tibble::rownames_to_column(selected_L2FC,var = "geneid")
      cluster_data_long <- left_join(cluster_data_long,selected_L2FC,by="geneid")
      
      cluster_data_long$comparison <- sub("*_Log2FC", "", x = cluster_data_long$comparison)
      cluster_data_long$comparison <- gsub("_vs_", "\nvs\n", x = cluster_data_long$comparison)
      cluster_data_long$comparison <- gsub("_", " ", x = cluster_data_long$comparison)
      global_sL2FC3 <<- cluster_data_long 
      
      incProgress(0.2, detail = "Plotting heatmap")
      #test_sort <<- fc_heat_sample_sort_data()
      #test_sortBy <<- fc_sort_sample_list_data()
      
      #Plot heatmaps using ggplot
      #Order genes based on comparison
      if (fc_heat_sample_sort_data() == FALSE){
        fc_heatmap <- ggplot(cluster_data_long,aes(comparison,reorder(geneid,.data[[fc_sort_sample_list_data()]]),fill=log2FC)) + 
          geom_tile() +
          scale_fill_gradient2(low="blue",mid="white",high="red") +
          theme(axis.text.y = element_blank(),axis.ticks.y=element_blank())
      }
      #Order genes based on cluster
      if (fc_heat_sample_sort_data() == TRUE){
        fc_heatmap <- ggplot(cluster_data_long,aes(comparison,reorder(geneid,as.numeric(clusters)))) + 
          geom_tile(aes(fill=log2FC)) +
          scale_fill_gradient2(low="blue",mid="white",high="red") +
          theme(axis.text.y = element_blank(),axis.ticks.y=element_blank()) + 
          geom_ysidebar(aes(yfill=clusters,ycolor=clusters)) +
          ggside(y.pos = "right")
      }
      fc_heatmap <- fc_heatmap +
        ylab(expression(log[2](Fold~Change)))+
        xlab("Genotype comparison")
      
      
      heatmap_download <<- fc_heatmap
      
      print(fc_heatmap)
      
      incProgress(0.2, detail = "Finished!")
    })
  })
  
  ##Sashimi plots
  output$browser <- renderPlot({
    
    withProgress(message = "Making Genome Browser", value = 0, {
      if(sashimi_data()){
        trackType <- c(track_type_data(), "sashimi")
      }else{
        trackType <- c(track_type_data())
      }
      
      fileList = bamFileList
      
      incProgress(0.2, detail = "Indexing BAM")
      
      
      incProgress(0.2, detail = "Finding gene annotation")
      
      #take GTF and subset just the genes (with useful positional info) and subset intron/exon lists for each gene
      info_gtf <- gtf %>% dplyr::select(feature, seqname, start, end, gene_symbol, gene_id)
      gene_info_gtf <- filter(info_gtf, feature == "gene")
      exon_info_gtf <- filter(info_gtf, feature == "exon")
      
      gene_row <- grep(get_geneID(browser_gene_list_data(), gtf, current_fb_ids), gene_info_gtf$gene_id)
      chrom <- gene_info_gtf$seqname[gene_row]
      start <- gene_info_gtf$start[gene_row] - plusStart_data()
      end <- gene_info_gtf$end[gene_row] + plusEnd_data()
      
      gtrack <- GenomeAxisTrack(scale = 0.25,
                                labelPos = "below",
                                cex = 1.5,
                                col = "black")
      grtrack <- GeneRegionTrack(txdb,
                                 genome = "dm6",
                                 chromosome = chrom,
                                 name = browser_gene_list_data(),
                                 showID = TRUE,
                                 size = 5,
                                 stacking = "squish")
      
      displayTracks <- c(gtrack)
      
      #Make the color vector object
      x <- c()
      if(!length(color_list_data()) == 0){
        for( i in 1:length(color_list_data())){
          test_color <- tolower(color_list_data()[i])
          if(!test_color %in% colors()){
            x <- append(x, "grey80")
          }else{
            x <- append(x, test_color)
          }
        }
        if(length(color_list_data()) < length(track_list_data())){
          len <- length(color_list_data()) + 1
          for(i in len:length(track_list_data())){
            x <- append(x, "grey80")
          }
        }
        if(length(color_list_data()) > length(track_list_data())){
          for(i in length(track_list_data()):length(color_list_data())){
            x[-length(x)]
          }
        }
      }else{
        for (i in 1:length(track_list_data())){
          x <- append(x, "grey80")
        }
      }
      assign("color_list", x)
      
      for(i in 1:length(track_list_data())){
        incProgress(0.1, detail = "Making tracks")
        
        indexedBam <- fileList[grep(track_list_data()[i], fileList)]
        
        if(max_height_data() > 0){
          tempTrack <- AlignmentsTrack(range = indexedBam,
                                       genome = "dm6",
                                       type = trackType,
                                       name = track_list_data()[i],
                                       size = 6,
                                       background.title = "grey40",
                                       fill.coverage = color_list[i],
                                       ylim = c(0, max_height_data()))
        }else{
          tempTrack <- AlignmentsTrack(range = indexedBam,
                                       genome = "dm6",
                                       type = trackType,
                                       name = track_list_data()[i],
                                       background.title = "grey40",
                                       fill.coverage = color_list[i],
                                       size = 6)
        }
        
        displayTracks <- append(displayTracks, tempTrack)
      }
      
      displayTracks <- append(displayTracks, grtrack)
      
      incProgress(0.95, detail = "Printing tracks")
      
      # Try-Catch to avoid errors by the sashimiFilter line
      printed <- FALSE
      tryCatch({
        plotTracks(displayTracks,
                   chromosome = chrom,
                   from = start,
                   to = end,
                   sashimiFilter = introns,
                   lwd.sashimiMax = 3,
                   sashimiFilterTolerance = 50L)
        printed <- TRUE
      }, error = function(e) {
        print("caught an error")
      },
      finally = {
        if(printed == FALSE){
          plotTracks(displayTracks,
                     chromosome = chrom,
                     from = start,
                     to = end,
                     lwd.sashimiMax = 3,
                     sashimiFilterTolerance = 50L)
        }
      }
      )
    })
  },
  height = function(){
    session$clientData$output_browser_width * 0.75
    #input$dimension[2] * 0.75
  },
  width = function(){
    session$clientData$output_browser_width
  })
  
  
  
  
  ## Heirarchical cluster plot heatmaps
  output$cluster <- renderPlot({
    
    paletteLength <- 100
    #myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
    myColor <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(paletteLength)
    
    direction <- ""
    if(direction_of_change_data() == "Up"){
      direction <- "up"
    }
    if(direction_of_change_data() == "Down"){
      direction <- "down"
    }
    if(direction_of_change_data() == "Neither"){
      direction <- "none"
    }
    
    
    withProgress(message = "Making Plot", value = 0, {
      
      incProgress(0.10, detail = "Beginning Clustering")
      
      ## filter norm_counts to include only genes meeting the selection criteria
      keep_genes <<- filter_genes2(samples_list = cluster_sample_data(),
                                   change_direction = "none",
                                   diff_genes_only = diff_only_data())
      norm_counts_clust <- norm_counts[rownames(norm_counts) %in% keep_genes,]
      
      ## convert norm_counts to data frame
      norm_counts_clust <- as.data.frame(norm_counts_clust)
      cluster_data <- norm_counts_clust
      
      #change the colnames to the basenames (more readable)
      for (i in sample_table$baseName){
        colnames(cluster_data)[grep(i, colnames(cluster_data))] <- i
      }
      
      cluster_table <<- cluster_data
      
      #### PERFORM CLUSTERING
      
      incProgress(0.20, detail = "Beginning Clustering")
      ## Perform heirarchical clustering
      if (cluster_method_data() == "heirarchical") {
        
        cluster_data <- perform_hc(what = "cd", cluster_n = cluster_n_data(), input_counts = cluster_data)
        hc <- perform_hc(what = "hc", cluster_n = cluster_n_data(), input_counts = cluster_data) 
        
        # ## compute distance matrix
        # sampleDists <- dist(cluster_data)
        # 
        # ## perform hierarchecal clustering
        # hc <- hclust(sampleDists)
        # 
        # 
        # ## extract cluster assignments for each gene
        # hc.cutree <- cutree(hc, cluster_n_data())
        # clusters <- data.frame(clusters = as.factor(hc.cutree))
        # clusters <- tibble::rownames_to_column(clusters, var = "gene_id")
        #cluster_data$clusters <- clusters$clusters
      }
      ## Perform kmeans clustering
      if (cluster_method_data() == "kmeans") {
        
        cluster_data <- perform_kc(what = "cd", cluster_n = cluster_n_data(), input_counts = cluster_data)
        kc <- perform_kc(what = "kc", cluster_n = cluster_n_data(), input_counts = cluster_data) 
        
        # set.seed(20)
        # kc <- kmeans(cluster_data, centers=cluster_n_data(), nstart = 1000, iter.max = 20)
        # kClusters <- as.factor(x = kc$cluster)
        # clusters <- data.frame(clusters = kClusters)
        # clusters <- tibble::rownames_to_column(clusters, var = "gene_id")
        # cluster_data$clusters <- clusters$clusters
      }
      
      #### PERFORM PLOTTING
      incProgress(0.5, detail = "Beginning plotting")
      if (cluster_method_data() == "heirarchical") {
        
        order <- hc$order
        cluster_data <- cluster_data[order,]
        cluster_data_out <<- cluster_data
        hm_data <- as.matrix(dplyr::select(cluster_data, -(clusters)))
        cluster_info <- data.frame(clusters = cluster_data$clusters)
        cluster_info$clusters <- as.character(cluster_info$clusters)
        
        #removed as.character from inline and added extra step
        #hm_info_image_out_pre <<- cluster_info
        cluster_color <- data.frame(clusters = unique(cluster_info$clusters),
                                    color = I(brewer.pal(length(unique(cluster_info$clusters)), name = 'Set3')))
        cluster_info <- left_join(cluster_info, cluster_color, by = "clusters")
        
        #hm_info_image_out <<- cluster_info
        #hm_data_image_out <<- hm_data
        #hm_cluster_image_out <<- cluster_color
        
        
        ### WORKS BEST 7/19/22 MN
        
        ## Simplify cluster info to add the cluster groups
        rownames(cluster_info) <- rownames(hm_data)
        cluster_info <- subset(cluster_info, select = -color)
        
        ph <- pheatmap(hm_data,
                       annotation_names_row = F,
                       color = myColor,
                       show_rownames = F,
                       show_colnames = T,
                       cutree_rows = cluster_n_data(),
                       treeheight_row = 0,
                       treeheight_col = 0,
                       annotation_row = cluster_info)
        
        ph
        
        ## Essential for the ability to download plots
        ph_data <<- hm_data
        ph_color <<- myColor
        ph_rows <<- cluster_n_data()
        ph_anno <<- cluster_info
        
        #rownames(cluster_info) <- make.unique(cluster_info[["clusters"]])
        #annotation <- cluster_info["clusters"]
        #annotation <- sample(cluster_color$clusters, nrow(cluster_color), replace = TRUE)
        
      }
      
      if (cluster_method_data() == "kmeans") {
        cluster_data <- tibble::rownames_to_column(cluster_data, var = "gene_id") %>%
          arrange(clusters) %>%
          tibble::column_to_rownames(var = "gene_id")
        hm_data <- as.matrix(dplyr::select(cluster_data, -(clusters)))
        cluster_info <- data.frame(clusters = as.character(cluster_data$clusters))
        cluster_color <- data.frame(clusters = unique(cluster_info$clusters),
                                    color = I(brewer.pal(length(unique(cluster_info$clusters)), name = 'Set3')))
        cluster_info <- left_join(cluster_info, cluster_color, by = "clusters")
        
        #hm_info_image_out <<- cluster_info
        #hm_data_image_out <<- hm_data
        #hm_cluster_image_out <<- cluster_color
        
        ## Simplify cluster info to add the cluster groups
        rownames(cluster_info) <- rownames(hm_data)
        cluster_info <- subset(cluster_info, select = -color)
        
        ph <- pheatmap(hm_data,
                       annotation_names_row = F,
                       color = myColor,
                       show_rownames = F,
                       show_colnames = T,
                       cutree_rows = cluster_n_data(),
                       treeheight_row = 0,
                       treeheight_col = 0,
                       annotation_row = cluster_info)
        ph
        
        ## Essential for the ability to download plots
        ph_data <<- hm_data
        ph_color <<- myColor
        ph_rows <<- cluster_n_data()
        ph_anno <<- cluster_info
        
      }
    })
    # },
    # height = function(){
    #   session$clientData$output_cluster_width * 0.75
    #   #input$dimension[2] * 0.75
    # })
  })
  
  
  
  
  
  ## Heirarchical cluster plot tables
  output$cluster_table <- DT::renderDataTable({
    
    ## plot heatmap of clustered data
    paletteLength <- 100
    myColor <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(paletteLength)
    
    direction <- ""
    if(direction_of_change_data() == "Up"){
      direction <- "up"
    }
    if(direction_of_change_data() == "Down"){
      direction <- "down"
    }
    if(direction_of_change_data() == "Neither"){
      direction <- "none"
    }
    
    
    withProgress(message = "Making Table", value = 0, {
      
      incProgress(0.33, detail = "Beginning Clustering")
      #if (cluster_by_data() == "counts") {
      
      ## filter norm_counts to include only genes meeting the selection criteria
      
      #keep_genes <- filter_genes2(samples_list =  cluster_sample_data(),
      #                           change_direction = direction)
      #print(head(keep_genes))
      
      # #norm_counts_clust <- norm_counts[rownames(norm_counts) %in% keep_genes,]
      # norm_counts_clust <- norm_counts
      # 
      # ## convert norm_counts to data frame
      # norm_counts_clust <- as.data.frame(norm_counts_clust)
      # cluster_data <- norm_counts_clust
      
      #print(head(cluster_data))
      #}
      # if (cluster_by_data() == "logFC") {
      #   
      #   ## filter all_logfc to contain only genes meeting the selection criteria
      #   keep_genes <- filter_genes2(samples_list =  cluster_sample_data(),
      #                               change_direction = direction)
      #   all_logfc_clust <- all_logfc[rownames(all_logfc) %in% keep_genes,]
      #   cluster_data <- all_logfc_clust
      #   
      #   #print(head(cluster_data))
      # }
      
      #### PERFORM CLUSTERING
      
      incProgress(0.2, detail = "Beginning Clustering")
      
      ## Perform heirarchical clustering
      if (cluster_method_data() == "heirarchical") {
        
        cluster_data <- perform_hc(what = "cd", cluster_n = cluster_n_data())
        hc <- perform_hc(what = "hc", cluster_n = cluster_n_data()) 
        
      }
      
      ## Perform kmeans clustering
      if (cluster_method_data() == "kmeans") {
        
        cluster_data <- perform_kc(what = "cd", cluster_n = cluster_n_data(), input_counts = cluster_data)
        kc <- perform_kc(what = "kc", cluster_n = cluster_n_data(), input_counts = cluster_data) 
        
      }
      
      results_table <- tibble::rownames_to_column(cluster_data, var = "gene_id")
      gene_data_name <- comparisons[1]
      gene_data <- read.table(paste("./output_data/",gene_data_name,"_resOrdered.txt", sep = ""), sep = "\t", header = TRUE, as.is = TRUE)
      colnames(gene_data)[1] <- "gene_id"
      results_table <- left_join(results_table, gene_data, by = "gene_id")
      
      results_table <- results_table %>% dplyr::select(gene_id, clusters)
      
      #write.csv(results_table, file = "raw_cluster_table.csv")
      cluster_table_out <- as.data.frame(results_table)
      DT::datatable(results_table)
      
      
    })
  })
  
  ## Right now only compares to full sets of genes within clusters. No difference in expression
  output$cluster_ont <- DT::renderDataTable({
    
    direction <- ""
    if(direction_of_change_data() == "Up"){
      direction <- "up"
    }
    if(direction_of_change_data() == "Down"){
      direction <- "down"
    }
    if(direction_of_change_data() == "Neither"){
      direction <- "none"
    }
    
    bg_direction <- ""
    if(direction_bg_change_data() == "Up"){
      bg_direction <- "up"
    }
    if(direction_bg_change_data() == "Down"){
      bg_direction <- "down"
    }
    if(direction_bg_change_data() == "Neither"){
      bg_direction <- "none"
    }
    
    #Filter out genes
    keep_genes <<- filter_genes2(samples_list = cluster_sample_data(),
                                 change_direction = "none",
                                 diff_genes_only = diff_only_data())
    norm_counts_clust <- norm_counts[rownames(norm_counts) %in% keep_genes,]
    
    ## convert norm_counts to data frame
    norm_counts_clust <- as.data.frame(norm_counts_clust)
    cluster_data <- norm_counts_clust
    
    cluster_data <- perform_hc(what = "cd", cluster_n = cluster_n_data(), input_counts = cluster_data)
    
    chosen <<- cluster_GO_data()
    background <<- cluster_bg_data()
    GO_raw <<- tibble::rownames_to_column(cluster_data['clusters'], "genes")
    
    GO_select <- GO_raw[GO_raw$clusters %in% chosen, ]
    
    # Now, if the user just wants up or down regulated genes in the cluster.. filter for that.
    if(direction == "up"){
      keep_genes <- filter_genes2(samples_list = cluster_sample_data(),
                                  change_direction = direction,
                                  diff_genes_only = diff_only_data())
      GO_select <- GO_select[GO_select$genes %in% keep_genes,]
    }else if(direction == "down"){
      keep_genes <- filter_genes2(samples_list = cluster_sample_data(),
                                  change_direction = direction,
                                  diff_genes_only = diff_only_data())
      GO_select <- GO_select[GO_select$genes %in% keep_genes,]
    }
    
    go_test <<- GO_select
    
    # If the user is using other clusters as bg, filter up/down regulated genes
    # Might be broken because filter_genes2(samples) may be messing up
    if(background[1] != "all"){
      GO_bg <<- GO_raw[GO_raw$clusters %in% background, ]
      if(direction == "up"){
        keep_genes <- filter_genes2(samples_list = cluster_sample_data(),
                                    change_direction = bg_direction,
                                    diff_genes_only = diff_only_data())
        GO_bg <<- GO_bg[GO_bg$genes %in% keep_genes,]
      }else if(direction == "down"){
        keep_genes <- filter_genes2(samples_list = cluster_sample_data(),
                                    change_direction = bg_direction,
                                    diff_genes_only = diff_only_data())
        GO_bg <<- GO_bg[GO_bg$genes %in% keep_genes,]
      }
    }
    #chosen <- paste(chosen, sep = "_")
    #GOList <- na.omit(c(fbgn_diff_genes[chosen]))
    
    #eg = bitr(GO_select, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Dm.eg.db")
    #ebg = bitr(GO_bg, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Dm.eg.db")
    # ggo <- groupGO(gene = eg$ENTREZID,
    #                OrgDb = org.Dm.eg.db,
    #                ont = "CC",
    #                level = 3,
    #                readable = TRUE)
    
    if(background[1] != "all"){
      ego <- enrichGO(gene = GO_select$genes,
                      universe = GO_bg$genes,
                      keyType = "SYMBOL",
                      OrgDb = org.Dm.eg.db,
                      ont = "CC",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.01,
                      qvalueCutoff = 0.05
      )
      
    }else{
      ego <- enrichGO(gene = GO_select$genes,
                      keyType = "SYMBOL",
                      OrgDb = org.Dm.eg.db,
                      ont = "CC",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.01,
                      qvalueCutoff = 0.05
      )
    }
    # }else{
    #   ego <- data.frame("Please run the heatmap gene ontology first")
    
    ego <- as.data.frame(ego)
    if(nrow(ego) == 0){
      ego <- data.frame("V1" = "No enrichment")
    }
    DT::datatable(ego)
  })
  
  ## Heirarchical cluster box plots 
  output$cluster_box <- renderPlot({
    
    ## plot heatmap of clustered data
    paletteLength <- 100
    myColor <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(paletteLength)
    
    ## Do some formatting of the inputs neccessary for the filter_genes2 function
    direction <- ""
    if(direction_of_change_data() == "Up"){
      direction <- "up"
    }
    if(direction_of_change_data() == "Down"){
      direction <- "down"
    }
    if(direction_of_change_data() == "Neither"){
      direction <- "none"
    }
    
    diff_only <- diff_only_data()
    
    withProgress(message = "Making Box plot(s)", value = 0, {
      
      incProgress(0.25, detail = "Beginning Clustering")
      #if (cluster_by_data() == "counts") {
      ## filter norm_counts to include only genes meeting the selection criteria
      
      keep_genes <- filter_genes2(samples_list =  cluster_sample_data(),
                                  diff_genes_only = diff_only,
                                  change_direction = direction)
      
      norm_counts_clust <- norm_counts[rownames(norm_counts) %in% keep_genes,]
      
      ## convert norm_counts to data frame
      norm_counts_clust <- as.data.frame(norm_counts_clust)
      
      #### PERFORM CLUSTERING
      
      incProgress(0.5, detail = "Beginning Clustering")
      
      ## Perform heirarchical clustering
      if (cluster_method_data() == "heirarchical") {
        
        # ## compute distance matrix
        # sampleDists <- dist(cluster_data)
        # #print("Sample distributions calculated")
        # 
        # ## perform hierarchecal clustering
        # hc <- hclust(sampleDists)
        # #print("Hclust object created")
        # 
        # ## extract cluster assignments for each gene
        # hc.cutree <- cutree(hc, cluster_n_data())
        # clusters <- data.frame(clusters = as.factor(hc.cutree))
        # clusters <- tibble::rownames_to_column(clusters, var = "gene_id")
        # cluster_data$clusters <- clusters$clusters
        
        cluster_data <- perform_hc(what = "cd", cluster_n = cluster_n_data(), input_counts = norm_counts_clust)
        hc <- perform_hc(what = "hc", cluster_n = cluster_n_data(), input_counts = norm_counts_clust)
      }
      
      
      ## Perform kmeans clustering
      if (cluster_method_data() == "kmeans") {
        cluster_data <- perform_kc(what = "cd", cluster_n = cluster_n_data(), input_counts = cluster_data)
        kc <- perform_kc(what = "kc", cluster_n = cluster_n_data(), input_counts = cluster_data) 
        
        # set.seed(20)
        # kc <- kmeans(cluster_data, centers=cluster_n_data(), nstart = 1000, iter.max = 20)
        # kClusters <- as.factor(x = kc$cluster)
        # clusters <- data.frame(clusters = kClusters)
        # clusters <- tibble::rownames_to_column(clusters, var = "gene_id")
        # cluster_data$clusters <- clusters$clusters
      }
      
      ## Get the common names of all the replicates
      incProgress(0.2, detail = "Format Clustering")
      
      ## Using the clustered data, get the mean of each row (gene) in each replicate group
      tibbleList <- list()
      for (i in repList2){
        temp_mean <- as_tibble(cluster_data) %>%dplyr::select(contains(i))
        temp_mean <- as.data.frame(rowMeans(temp_mean))
        colnames(temp_mean) <- i
        assign(i, temp_mean)
        tibbleList <- list(tibbleList, temp_mean)
      }
      
      ## Combine individual tibbles plus the cluster column from the original clustered dataframe
      cluster_tbbl <- bind_cols(tibbleList, clusters[2])
      
      ##Filters by cluster group and puts them in their own tibbles
      for( i in 1:max(as.numeric(cluster_tbbl[,length(cluster_tbbl[1,])]))){
        clust_name <- paste("clust_", i, sep = "")
        temp_clust <- filter(cluster_tbbl, clusters == i) %>%dplyr::select(-one_of("clusters")) %>% tidyr::gather(.)
        assign(clust_name, temp_clust)
      }
      plot_list <- list()
      
      incProgress(1, detail = "Begin Plotting")
      for( i in 1:max(as.numeric(cluster_tbbl[,length(cluster_tbbl[1,])]))){
        plot_name <- paste("p", i, sep = "")
        clust_title <- paste("Cluster ", i, sep = "")
        clust_name <- paste("clust_", i, sep = "")
        df <- eval(as.symbol(clust_name))
        df$key2 <- factor(df$key, levels = unlist(strsplit(box_order_data(), " ")))
        assign(clust_name, df)
        
        #Remove NAs
        df <- df[complete.cases(df), ]
        test <- df
        
        temp_plot <- ggplot(data=df, aes(y=value, x=key2)) +
          geom_boxplot() +
          labs(title=bquote(.(clust_title)~phantom(0)["n= "][.(nrow(filter(cluster_tbbl, clusters == i)))]),x="",y=expression('Log'[2]*' Counts')) +
          theme(plot.title = element_text(hjust = 0.5))
        
        assign(plot_name, temp_plot)
        
        if(i %in% plot_which_data()){
          plot_list[[clust_name]] <- eval(as.symbol(plot_name))
        }
      }
      cluster_box_image_out <<- plot_list
      grid.arrange(grobs = plot_list)
      # }
    })
  },
  height = function(){
    session$clientData$output_cluster_box_width * 0.75
    #input$dimension[2] * 0.75
  })
  
  
  
  output$count <- renderPlot({
    count_list <- list()
    
    #trying to separate genes for boxplot
    #input_gene <- strsplit(input_gene,",")
    
    test_genes <<- count_data()
    for(input_gene in count_data()){
      
      ## Replace the column names of each replicate with the "common names"
      goi <- get_geneID(input_gene, gtf, current_fb_ids)
      
      goi_count <- as.data.frame(fbgn_counts)
      
      ## Find which columns are needed and select only those as well as only the row of the gene in question
      test2 <- vector('numeric')
      
      #for(i in ii){
      for(i in sample_data()){
        test2[length(test2)+1:length(grep(i, colnames(goi_count)))] <- grep(i, colnames(goi_count))
      }
      
      goi_count <- dplyr::filter(goi_count, gene_id == goi)
      goi_count <- dplyr::select(goi_count, test2)
      
      ## Replace the column names of each replicate with the "common names"
      for(i in 1:length(repList2)){
        replace <- grep(repList2[i], colnames(goi_count))
        for(j in replace){
          colnames(goi_count)[j] <- repList2[i]
        }
      }
      
      ## Do some formatting and summarizing to convert the normalized counts into means and SDs
      
      temp_goi <- tibble::rownames_to_column(as.data.frame(t(goi_count)), "key")
      temp_goi$key <- rownames(t(goi_count))
      
      goi_count <- temp_goi
      colnames(goi_count)[2] <- "value"
      grouped_count <- group_by(goi_count, key)
      summ_goi_count <- summarise(grouped_count, mean=mean(value), sd=sd(value), N = sum(!is.na(value)), upper_limit = mean + sd/sqrt(N), lower_limit = mean - sd/sqrt(N))
      summ_goi_count$key <- factor(summ_goi_count$key, levels = sample_data())
      
      ## Make the plot title
      title <- bquote(italic(.(input_gene))~' Log'[2]*' Counts')
      
      # Print out the plot
      c1 <- ggplot(summ_goi_count, aes(x=key, y=mean)) +
        theme_minimal() +
        geom_errorbar(aes(ymax=upper_limit, ymin=lower_limit), color = "grey60", width=0.25) +
        geom_point() +
        labs(title=bquote(.(title)),x="",y=expression('Log'[2]*' Counts')) +
        theme(plot.title = element_text(hjust = 0.5))
      
      if(lc_data() == "Connect"){
        c1 <- c1 +geom_line(group = 1)
      }
      count_list[[match(input_gene, count_data())]] <- c1
      gene_exp_image_out <<- count_list
      # 
    }
    grid.arrange(grobs = count_list)
  },
  height = function(){
    session$clientData$output_count_width * 0.75
    #input$dimension[2] * 0.75
  })
  
  output$fc <- renderPlot({
    fc_list <- list()
    
    test_fc <<- fc_data()
    for(input_gene in fc_data()){
      
      goi <- get_geneID(input_gene, gtf, current_fb_ids)
      
      ## broken because at some point all_logfc had fbgns in a column called "gene_id"
      goi_logfc <- dplyr::filter(all_logfc, gene_id == goi)
      
      temp_fc <- vector('numeric')
      for(i in fc_sample_data()){
        temp_fc[length(temp_fc)+1:length(grep(i, colnames(goi_logfc)))] <- grep(i, colnames(goi_logfc))
      }
      goi_logfc <-dplyr::select(goi_logfc, temp_fc)
      goi_logfc <- tidyr::gather(goi_logfc)
      
      max_fc = ceiling(max(abs(goi_logfc[2])))
      min_fc = -max_fc
      
      ## Make the plot title
      title <- bquote(italic(.(input_gene))~' Log'[2]*' Fold Change')
      
      ## Make plot x labels
      fc_labels <- list()
      for(i in 1:length(fc_sample_data())){
        pre_label <- unlist(strsplit(fc_sample_data()[i], "_vs_"))
        fc_labels[[i]] <- bquote(atop(.(pre_label[1])~' vs',.(pre_label[2])))
      }
      
      ## Plot Log 2 Fold Change
      f1 <- ggplot(goi_logfc, aes(y=value, x=key)) +
        geom_bar(aes(fill = key), stat = "identity")+
        theme_minimal() +
        labs(title=bquote(.(title)),x="",y=expression('Log'[2]*' Fold Change'))+
        ylim(c(min_fc, max_fc)) +
        scale_x_discrete(labels=fc_labels) +
        theme(legend.position="none") +
        theme(plot.title = element_text(hjust = 0.5))
      
      fc_list[[match(input_gene, fc_data())]] <- f1
      fc_image_out <<- fc_list
      grid.arrange(grobs = fc_list)
      
    }
  },
  height = function(){
    session$clientData$output_fc_width * 0.75
    #input$dimension[2] * 0.75
  })
  
  output$res_table <- DT::renderDataTable({
    
    test_res <- paste(first_data(), "_vs_", second_data(), "_resOrdered", sep ="")
    if(!exists(test_res)){
      test_res <- paste(second_data(), "_vs_", first_data(), "_resOrdered", sep = "")
    }
    test_res <- as.data.frame(eval(as.symbol(test_res)))
    raw_res_table_out <- test_res
    DT::datatable(test_res)
    
  })
  
  output$res_table_text <- renderText({
    
    test_res <- paste(first_data(), "_vs_", second_data(), "_resOrdered", sep ="")
    if(!exists(test_res)){
      test_res <- paste(second_data(), "_vs_", first_data(), "_resOrdered", sep = "")
    }
    test_res <- gsub("_resOrdered","", test_res)
    test_res <- gsub("_", " ", test_res)
    print(test_res)
  })
  
  ##PCA plot
  output$PCA <- renderPlot({
    # Created in diff_exp
    test_pca <<- PCA_data()
    
    # Get the indexes of the samples that need to be removed
    exclude <<- pca_list_data()
    
    if(sum(nchar(exclude)) > 0){
      ex <<- lapply(paste(exclude,collapse="|"), grep, rownames(colData(dds2)))
      dds2_ex <- dds2[,-unlist(ex)]
    }else{
      dds2_ex <- dds2
    }
    object <-vst(dds2_ex,blind=FALSE)
    intgroup="sample"
    ntop=500
    
    # calculate the variance for each gene
    rv <- rowVars(assay(object))
    # select the ntop genes by variance
    select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
    # perform a PCA on the data in assay(x) for the selected genes
    pca <- prcomp(t(assay(object)[select,]))
    # the contribution to the total variance for each component
    percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
    if (!all(intgroup %in% names(colData(object)))) {
      stop("the argument 'intgroup' should specify columns of colData(dds)")
    }
    intgroup.df <- as.data.frame(colData(object)[, intgroup, drop=FALSE])
    # add the intgroup factors together to create a new grouping factor
    group <- if (length(intgroup) > 1) {
      factor(apply( intgroup.df, 1, paste, collapse=":"))
    } else {
      colData(object)[[intgroup]]
    }
    # assembly the data for the PCA plot
    d <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], group=group, intgroup.df, name=colnames(object))
    pc1var <- paste0("(", round(percentVar[1]*100,2), "%", ")")
    pc2var <- paste0("(", round(percentVar[2]*100,2), "%", ")")
    pca_xlab <- paste("PC1", pc1var)
    pca_ylab <- paste("PC2", pc2var)
    
    ##only for my test dataset
    #d$name <- gsub("-", "_", d$name)
    
    #This needs to be changed for future
    #d$name <- sub("_dm6_trim_q5_dupsRemoved.bam$", "", sub("Bam*.","",d$name))
    
    d$name <- str_sub(d$name,5,-31)
    
    pca_out <<- ggplot(d,aes(PC1,PC2,color=sample)) +
      geom_point(size=5) +
      theme_minimal(base_size = 14) +
      labs(x = pca_xlab, y = pca_ylab) +
      geom_mark_ellipse(aes(fill = sample,color = sample), expand = unit(0.5,"mm"))
    
    # Scree Plot
    pca_col <- ncol(pca$x)
    var_explained_df <- data.frame(PC= paste0("PC",1:pca_col),
                                   var_explained=(pca$sdev)^2/sum((pca$sdev)^2))
    colors <- rep("#7BAFD4", pca_col)
    var_explained_df$PC <- factor(var_explained_df$PC, levels = var_explained_df$PC)
    scree_out <<- ggplot(var_explained_df, aes(x=PC,y=var_explained, group=1, fill = PC))+
      geom_col()+
      geom_point(size=4)+
      geom_line()+
      theme_minimal(base_size = 14) +
      scale_fill_manual(values = colors)+
      labs(x = pca_xlab, y = pca_ylab) +
      labs(title="Scree plot: PCA on scaled data")+
      theme(legend.position = "none")
    
    
    
    
    if(PCA_data() == "PCA plot"){
      pca_out
    }else{
      scree_out
    }
    
  })
  
  output$volcano <- renderPlot({
    
    textSize <- get_textSize(volcano_samples_data())
    
    withProgress(message = "Making Volcano Plot", value = 0, {    
      
      #set the colors for the plot; more than 3 colors is ignored, fewer sets the rest to black
      volcano_c <- strsplit(volcano_colors(), ", ")
      if(length(volcano_c[[1]]) < 3){
        volcano_c[[1]][(length(volcano_c[[1]])+1):3] <- "black"
      }
      
      
      for (i in volcano_samples_data()) {
        i_name <- paste(i, "_resOrdered.txt", sep = '')
        fn <- paste("./output_data/", i_name, sep = "")
        i_data <- read.table(file = fn, sep = "\t", header = TRUE, as.is = TRUE)
        #i_data <- filter_genes(i_data, filter_CRY = volcano_exclude_CRY_data(), filter_zld = volcano_zld_only_data())
        assign(i_name, i_data)
      }
      
      incProgress(0.1, detail = "Calculating plot scale")
      ## calculate axis range based on max values in comparisons
      ranges <- data.frame(xmin = rep(NA, times = length(volcano_samples_data())),xmax = rep(NA, times = length(volcano_samples_data())), ymin = rep(NA, times = length(volcano_samples_data())), ymax = rep(NA, times = length(volcano_samples_data())) )
      for (i in 1:length(volcano_samples_data())) {
        
        # for each comparison, get the DEseq results table
        i_name <- paste(volcano_samples_data()[i], "_resOrdered", sep = '')
        i_data <- eval(as.symbol(i_name))
        i_data <- as.data.frame(i_data)
        #i_data <- as.data.frame(i_data@listData, row.names = i_data@rownames)
        i_data <- dplyr::filter(i_data, !is.na(padj))
        ranges[i,1:2] <- range(i_data$log2FoldChange)
        ranges[i,3:4] <- range(-log10(i_data$padj))
      }
      
      volcano_list <- list()
      #j <-2
      for (i in volcano_samples_data()) {
        print(length(volcano_samples_data()))
        # for each comparison, get the DEseq results table
        i_name <- paste(i, "_resOrdered", sep = '')
        i_data <- eval(as.symbol(i_name))
        i_data <- as.data.frame(i_data)
        i_data <- dplyr::filter(i_data, !is.na(padj))
        
        ## add color_group column to data frame specifying the different groups of genes to be colored differently on the plot
        i_data$color_group <- "ns"
        
        top_hits <- data.frame(matrix(ncol = length(i_data[1,]), nrow = 0))
        colnames(top_hits) <- colnames(i_data)
        
        for (j in 1:nrow(i_data)) {
          if (i_data$padj[j] <= as.numeric(volcano_padj()) && i_data$log2FoldChange[j] >= volcano_fc()[2]) {
            i_data[j,"color_group"] <- "sig_up"
            top_hits[j,] <- i_data[j,]
          }
          else if (i_data$padj[j] <= as.numeric(volcano_padj()) && i_data$log2FoldChange[j] <= volcano_fc()[1]){
            i_data[j,"color_group"] <- "sig_down"
            top_hits[j,] <- i_data[j,]
          }
          else{
            i_data[j,"color_group"] <- "ns"
            
          }
        }
        
        top_hits <- top_hits[order(top_hits$padj),]
        top_hits <- head(top_hits, volcano_x())
        
        incProgress(j/10, detail = "Individual plots")
        j <- j+1
        
        xmin <- min(ranges$xmin)
        xmax <- max(ranges$xmax)
        ymax <-max(ranges$ymax)
        
        titleIndex <- which(choices[] == i)
        newTitle <- names(titleIndex)
        
        # create volcano plot of differential expression for each comparison
        volcano <- ggplot(i_data, aes(log2FoldChange, -log10(padj), colour=color_group))
        
        # copy this later for MA plots!!
        #volcano <- ggplot(i_data, aes(log2FoldChange, baseMean, colour=color_group))
        
        # Add some formatting to the initial Volcano plot
        volcano <- volcano +
          geom_vline(xintercept=0, color = "gray80",size=0.5) +
          geom_hline(yintercept = 0, color = "gray80", size=0.5) +
          geom_point(size = volcano_point_size()) +
          scale_color_manual(values = c("ns" = volcano_c[[1]][3], "sig_up" = volcano_c[[1]][1], "sig_down" = volcano_c[[1]][2])) +
          labs(title=newTitle) +
          ylab(expression(-log[10](adjusted~p-value)))+
          xlab(expression(log[2](Fold~Change))) +
          ylim(-(0.09 * ymax), ymax + 0.05 * ymax) +
          xlim((xmin + 0.05 * xmin), (xmax + 0.05 * xmax))
        
        # Change the theme elements of the volcano plot
        volcano <- volcano +
          theme_bw() +
          theme(
            plot.background = element_blank(),
            axis.line = element_line(color = "black", size=1, linetype = "solid")
            
          )
        
        if(!"grids" %in% volcano_options()){
          volcano <- volcano +
            theme(
              panel.grid = element_blank(),
              panel.border = element_blank()
            )
        }
        if(!"legend" %in% volcano_options()){
          volcano <- volcano + 
            theme(
              legend.position = "none"
            )
        }
        if("pv" %in% volcano_options()){
          volcano <- volcano + geom_hline(yintercept = -log10(as.numeric(volcano_padj())), linetype="dashed")
        }
        if("fc" %in% volcano_options()){
          volcano <- volcano +
            geom_vline(xintercept = volcano_fc()[1], linetype="dashed") +
            geom_vline(xintercept = volcano_fc()[2], linetype="dashed")
          
        }
        
        ## remove axis labels and numbers
        #axis.title = element_blank(),
        #axis.text = element_blank()
        
        
        # If the user has chosen to label the top X genes, label them.
        if(length(top_hits[,1]) != 0){
          volcano <- volcano +
            geom_text_repel(
              data = top_hits,
              aes(label = rownames(top_hits)),
              color = "grey30",
              size = v_gene_size(),
              max.overlaps = 20,
              box.padding = unit(0.5, "lines"),
              max.iter = 3e3,
              point.padding = unit(0.25, "lines"),
              force = 2,
              segment.color = "grey50",
              show.legend = FALSE
            )
          
        }
        volcano_list[[match(i, volcano_samples_data())]] <- volcano
      }
      incProgress(0.80, detail = "Printing plots...")
      volcano_image_out <<- volcano_list
      grid.arrange(grobs = volcano_list)
    })
  },
  height = function(){
    session$clientData$output_volcano_width * 0.75
    #input$dimension[2] * 0.75
  })
  
  output$venn <- renderPlot({
    
    venn_list <<- list()
    vl <- venn_sample_list()
    
    if(length(vl) == 1){
      venn_list <- list(All_Genes = data.frame(gene_symbol = current_fb_ids$current_symbol)$gene_symbol,
                        sample = rownames_to_column(as.data.frame(get(paste(vl[1], "_resOrdered_sig", sep = ""))), var = "gene_symbol")$gene_symbol
      )
    }else if(length(vl) >= 7){
      #length(vl)
      vl <- vl[-(7:length(vl))]
    }else{
      for(i in 1:length(vl)){
        venn_list[paste(vl[i], sep = "")] <- list(rownames_to_column(as.data.frame(get(paste(vl[i], "_resOrdered_sig", sep = ""))), var = "gene_symbol")$gene_symbol)
      }
    }
    
    saved <<- venn_list
    
    # K20_venn <- list(K20A = rownames_to_column(as.data.frame(germline_l3mbt_mut_vs_soma_l3mbt_mut_resOrdered_sig), var = "gene_symbol")$gene_symbol,
    # K20R = rownames_to_column(as.data.frame(germline_wt_vs_soma_wt_resOrdered_sig), var = "gene_symbol")$gene_symbol)
    venn_list <- venn_list %>% Venn() %>% process_data()
    venn_gg <- ggplot() +
      geom_sf(aes(fill=count), data = venn_region(venn_list)) +
      geom_sf(size = 1, color = "black", data = venn_setedge(venn_list), show.legend = F) +
      geom_sf_text(aes(label = name), data = venn_setlabel(venn_list), size = 4) +
      geom_sf_label(aes(label=count), fontface = "bold", data = venn_region(venn_list),size = 4) +
      theme_void() + scale_fill_gradient(low = "#F4FAFE", high = "#4981BF")
    print(venn_gg)
    
    coln <- lapply(venn_list@region$item, function(col){
      for(i in 1:length(col)){
        return(paste(length(col),"_genes", sep = ""))
      }
    })
    
    max_length <- 0

    for(i in 1:length(venn_list@region$item)){
      if(length(venn_list@region$item[[i]]) > max_length){
        max_length <- length(venn_list@region$item[[i]])
      }
    }
    
    data <- data.frame(matrix(ncol=0,nrow=max_length))
    for(i in 1:length(venn_list@region$item)){
      new <- c(venn_list@region$item[[i]],
               rep(NA, max_length-length(venn_list@region$item[[i]])))
      data[ , ncol(data) +1] <- new
      colnames(data)[ncol(data)] <- paste(length(venn_list@region$item[[i]]), "_genes", sep = "")
    }
    
    #This is the table that is 'caught' by the output$venn_table function
    data_new <<- data[ , colSums(is.na(data)) < nrow(data)]
  })
  
  ## Venn tables
  output$venn_table <- DT::renderDataTable({
    
    #required to get the stupid thing to update :/
    dummy_variable <- venn_sample_list()
    DT::datatable(data_new)
  })
  
  ## MA plot
  output$ma <- renderPlot({
    
    textSize <- get_textSize(ma_samples_data())
    
    withProgress(message = "Making MA Plot", value = 0, {    
      
      #set the colors for the plot; more than 3 colors is ignored, fewer sets the rest to black
      ma_c <- strsplit(ma_colors(), ", ")
      if(length(ma_c[[1]]) < 3){
        ma_c[[1]][(length(ma_c[[1]])+1):3] <- "black"
      }
      
      for (i in ma_samples_data()) {
        i_name <- paste(i, "_resOrdered.txt", sep = '')
        fn <- paste("./output_data/", i_name, sep = "")
        i_data <- read.table(file = fn, sep = "\t", header = TRUE, as.is = TRUE)
        assign(i_name, i_data)
      }
      
      incProgress(0.1, detail = "Calculating plot scale")
      ## calculate axis range based on max values in comparisons
      ranges <- data.frame(xmin = rep(NA, times = length(ma_samples_data())),xmax = rep(NA, times = length(ma_samples_data())), ymin = rep(NA, times = length(ma_samples_data())), ymax = rep(NA, times = length(ma_samples_data())) )
      for (i in 1:length(ma_samples_data())) {
        
        # for each comparison, get the DEseq results table
        i_name <- paste(ma_samples_data()[i], "_resOrdered", sep = '')
        i_data <- eval(as.symbol(i_name))
        i_data <- as.data.frame(i_data)
        #i_data <- as.data.frame(i_data@listData, row.names = i_data@rownames)
        i_data <- dplyr::filter(i_data, !is.na(padj))
        ranges[i,1:2] <- range(i_data$log2FoldChange)
        ranges[i,3:4] <- range(i_data$baseMean)
      }
      
      ma_list <- list()
      #j <-2
      for (i in ma_samples_data()) {
        # for each comparison, get the DEseq results table
        i_name <- paste(i, "_resOrdered", sep = '')
        i_data <- eval(as.symbol(i_name))
        i_data <- as.data.frame(i_data)
        i_data <- dplyr::filter(i_data, !is.na(padj))
        
        ## add color_group column to data frame specifying the different groups of genes to be colored differently on the plot
        i_data$color_group <- "ns"
        
        top_hits <- data.frame(matrix(ncol = length(i_data[1,]), nrow = 0))
        colnames(top_hits) <- colnames(i_data)
        
        for (j in 1:nrow(i_data)) {
          if (i_data$padj[j] <= as.numeric(ma_padj()) && i_data$log2FoldChange[j] >= ma_fc()[2]) {
            i_data[j,"color_group"] <- "sig_up"
            top_hits[j,] <- i_data[j,]
          }
          else if (i_data$padj[j] <= as.numeric(ma_padj()) && i_data$log2FoldChange[j] <= ma_fc()[1]){
            i_data[j,"color_group"] <- "sig_down"
            top_hits[j,] <- i_data[j,]
          }
          else{
            i_data[j,"color_group"] <- "ns"
            
          }
        }
        
        top_hits <- top_hits[order(top_hits$padj),]
        top_hits <- head(top_hits, ma_x())
        
        incProgress(j/10, detail = "Individual plots")
        j <- j+1
        
        xmin <- min(ranges$xmin)
        xmax <- max(ranges$xmax)
        ymax <-max(ranges$ymax)
        
        titleIndex <- which(choices[] == i)
        newTitle <- names(titleIndex)
        
        # create volcano plot of differential expression for each comparison
        #volcano <- ggplot(i_data, aes(log2FoldChange, -log10(padj), colour=color_group))
        
        # copy this later for MA plots!!
        ma <- ggplot(i_data, aes(log2FoldChange, baseMean, colour=color_group))
        
        # Add some formatting to the initial ma plot
        ma <- ma +
          #geom_vline(xintercept=0, color = "gray80",size=0.5) +
          #geom_hline(yintercept = 0, color = "gray80", size=0.5) +
          geom_point(size = ma_point_size()) +
          scale_color_manual(values = c("ns" = ma_c[[1]][3], "sig_up" = ma_c[[1]][1], "sig_down" = ma_c[[1]][2])) +
          labs(title=newTitle) +
          xlab(expression(log[2](Fold~Change)))+
          ylab("mean read depth") +
          ylim(-(0.09 * ymax), ymax + 0.05 * ymax) +
          xlim((xmin + 0.05 * xmin), (xmax + 0.05 * xmax)) +
          coord_flip()
        
        # Change the theme elements of the ma plot
        ma <- ma +
          theme_bw() +
          theme(
            plot.background = element_blank(),
            axis.line = element_line(color = "black", size=1, linetype = "solid")
            
          )
        
        if(!"grids" %in% ma_options()){
          ma <- ma +
            theme(
              panel.grid = element_blank(),
              panel.border = element_blank()
            )
        }
        if(!"legend" %in% ma_options()){
          ma <- ma + 
            theme(
              legend.position = "none"
            )
        }
        if("pv" %in% ma_options()){
          ma <- ma + geom_hline(yintercept = -log10(as.numeric(ma_padj())), linetype="dashed")
        }
        if("fc" %in% ma_options()){
          ma <- ma +
            geom_vline(xintercept = ma_fc()[1], linetype="dashed") +
            geom_vline(xintercept = ma_fc()[2], linetype="dashed")
          
        }
        
        ## remove axis labels and numbers
        #axis.title = element_blank(),
        #axis.text = element_blank()
        
        
        # If the user has chosen to label the top X genes, label them.
        if(length(top_hits[,1]) != 0){
          ma <- ma +
            geom_text_repel(
              data = top_hits,
              aes(label = rownames(top_hits)),
              color = "grey30",
              size = 3,
              max.overlaps = 20,
              box.padding = unit(0.5, "lines"),
              max.iter = 3e3,
              point.padding = unit(0.25, "lines"),
              force = 2,
              segment.color = "grey50",
              show.legend = FALSE
            )
          
        }
        ma_list[[match(i, ma_samples_data())]] <- ma
      }
      incProgress(0.80, detail = "Printing plots...")
      ma_image_out <<- ma_list
      grid.arrange(grobs = ma_list)
    })
  },
  height = function(){
    session$clientData$output_ma_width * 0.75
    #input$dimension[2] * 0.75
  })
  
  
  
  ### DownloadHandler functions that will allow users to download images/tables
  
  output$download_browser_plot <- downloadHandler(
    filename = function(){
      paste("browser.pdf", sep="")
    },
    content = function(file) {
      
      if(sashimi_data()){
        trackType <- c(track_type_data(), "sashimi")
      }else{
        trackType <- c(track_type_data())
      }
      
      fileList = c("./input_data/CRY_10_1_S3_filtered_sorted.bam",
                   "./input_data/WT_10_1_S1_filtered_sorted.bam",
                   "./input_data/CRY_12_3_S7_filtered_sorted.bam",
                   "./input_data/WT_12_3_S3_filtered_sorted.bam",
                   "./input_data/CRY_13_3_S15_filtered_sorted.bam",
                   "./input_data/WT_13_3_S13_filtered_sorted.bam",
                   "./input_data/CRY_14_1_S5_filtered_sorted.bam",
                   "./input_data/WT_14_1_S1_filtered_sorted.bam",
                   "./input_data/CRY_NL_1_S9_filtered_sorted.bam")
      
      
      #take GTF and subset just the genes (with useful positional info) and subset intron/exon lists for each gene
      info_gtf <- gtf %>% dplyr::select(feature, seqname, start, end, gene_symbol, gene_id)
      gene_info_gtf <- filter(info_gtf, feature == "gene")
      exon_info_gtf <- filter(info_gtf, feature == "exon")
      
      print(head(info_gtf))
      print(head(gene_info_gtf))
      print(head(exon_info_gtf))
      
      print(browser_gene_list_data())
      
      gene_row <- grep(get_geneID(browser_gene_list_data(), gtf, current_fb_ids), gene_info_gtf$gene_id)
      chrom <- gene_info_gtf$seqname[gene_row]
      start <- gene_info_gtf$start[gene_row] - plusStart_data()
      end <- gene_info_gtf$end[gene_row] + plusEnd_data()
      
      gtrack <- GenomeAxisTrack(scale = 0.5,
                                labelPos = "below",
                                cex = 1.5,
                                col = "black")
      grtrack <- GeneRegionTrack(txdb,
                                 genome = "dm6",
                                 chromosome = chrom,
                                 name = browser_gene_list_data(),
                                 showID = TRUE,
                                 size = 5,
                                 stacking = "squish")
      
      displayTracks <- c(gtrack)
      
      #Make the color vector object
      x <- c()
      if(!length(color_list_data()) == 0){
        for( i in 1:length(color_list_data())){
          test_color <- tolower(color_list_data()[i])
          if(!test_color %in% colors()){
            x <- append(x, "grey80")
          }else{
            x <- append(x, test_color)
          }
        }
        if(length(color_list_data()) < length(track_list_data())){
          len <- length(color_list_data()) + 1
          for(i in len:length(track_list_data())){
            x <- append(x, "grey80")
          }
        }
        if(length(color_list_data()) > length(track_list_data())){
          for(i in length(track_list_data()):length(color_list_data())){
            x[-length(x)]
          }
        }
      }else{
        for (i in 1:length(track_list_data())){
          x <- append(x, "grey80")
        }
      }
      assign("color_list", x)
      
      for(i in 1:length(track_list_data())){
        indexedBam <- fileList[grep(track_list_data()[i], fileList)]
        
        if(max_height_data() > 0){
          tempTrack <- AlignmentsTrack(range = indexedBam,
                                       genome = "dm6",
                                       type = trackType,
                                       name = track_list_data()[i],
                                       size = 6,
                                       background.title = "grey40",
                                       fill.coverage = color_list[i],
                                       ylim = c(0, max_height_data()))
        }else{
          tempTrack <- AlignmentsTrack(range = indexedBam,
                                       genome = "dm6",
                                       type = trackType,
                                       name = track_list_data()[i],
                                       background.title = "grey40",
                                       fill.coverage = color_list[i],
                                       size = 6)
        }
        
        displayTracks <- append(displayTracks, tempTrack)
        #length(displayTracks)
      }
      displayTracks <- append(displayTracks, grtrack)
      
      # Try-Catch to avoid errors by the sashimiFilter line
      printed <- FALSE
      tryCatch({
        plotTracks(displayTracks,
                   chromosome = chrom,
                   from = start,
                   to = end,
                   sashimiFilter = introns,
                   lwd.sashimiMax = 3)
        printed <- TRUE
      }, error = function(e) {
        print("caught an error")
      },
      finally = {
        if(printed == FALSE){
          plotTracks(displayTracks,
                     chromosome = chrom,
                     from = start,
                     to = end,
                     lwd.sashimiMax = 3)
        }
      }
      )
      dev.off()
    }
  )
  
  output$download_res_table <- downloadHandler(
    filename = function(){
      paste(Sys.Date(),"_results_data.csv", sep = "")
    },
    content = function(file){
      write.csv(raw_res_table_out, file)
    }
  )
  
  output$download_go_table <- downloadHandler(
    filename = function(){
      paste(Sys.Date(),"_GO_results.csv", sep = "")
    },
    content = function(file){
      write.csv(go_table_out, file)
    }
  )
  
  output$download_cluster_table <- downloadHandler(
    
    filename = function(){
      paste(Sys.Date(),"_cluster_data.csv", sep = "")
    },
    content = function(file){
      write.csv(cluster_table_out, file)
    }
  )
  
  output$download_fc_plot <- downloadHandler(
    
    filename = function(){
      paste(fc_data(), "_log2FC_plots.pdf", sep="")
    },
    content = function(file) {
      device <- function(..., width, height) grDevices::png(..., width = width, height = height, res = 300, units = "in")
      ggsave(file = file, plot = grid.arrange(grobs = fc_image_out), device = "pdf", dpi = 300)
    }
  )
  
  output$download_gene_exp_plot <- downloadHandler(
    filename = function(){
      paste(count_data(), "_count_plots.pdf", sep="")
    },
    content = function(file) {
      device <- function(..., width, height) grDevices::png(..., width = width, height = height, res = 300, units = "in")
      ggsave(file = file, plot = grid.arrange(grobs = gene_exp_image_out), device = "pdf", dpi = 300)
    }
  )
  
  output$download_count_box_plot <- downloadHandler(
    filename = function(){
      paste("cluster_box_plots.pdf", sep="")
    },
    content = function(file) {
      device <- function(..., width, height) grDevices::png(..., res = 300)
      ggsave(file = file, plot = grid.arrange(grobs = cluster_box_image_out), device = "pdf", dpi = 300)
    }
  )
  
  output$download_fc_heatmap <- downloadHandler(
    filename = function(){
      paste("fc_heatmap.pdf", sep="")
    },
    content = function(file) {
      device <- function(..., width, height) grDevices::png(..., res = 300)
      ggsave(file = file, plot = heatmap_download, device = "pdf", dpi = 300)
    }
  )
  
  output$download_volcano_plot <- downloadHandler(
    filename = function(){
      paste("volcano_plots.pdf", sep="")
    },
    content = function(file) {
      device <- function(..., width, height) grDevices::png(..., res = 300)
      ggsave(file = file, plot = grid.arrange(grobs = volcano_image_out), device = "pdf", dpi = 300)
    }
  )
  
  output$download_count_hm_plot <- downloadHandler(
    filename = function(){
      paste("cluster_heatmap_plot.pdf", sep="")
    },
    content = function(file) {
      ## Below doesn't work yet
      pheatmap(ph_data,
               filename = file,
               annotation_names_row = F,
               color = ph_color,
               show_rownames = F,
               show_colnames = T,
               cutree_rows = ph_rows,
               treeheight_row = 0,
               treeheight_col = 0,
               annotation_row = ph_anno)
      
    }
    
  )
  
  
}
