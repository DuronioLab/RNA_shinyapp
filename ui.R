ui <- shinyUI(navbarPage("RNA-seq Analysis",
                         tabPanel(title = "Welcome",
                                  shiny::tags$h2("Basic RNA-seq analysis Shiny App"),
                                  shiny::tags$br(),
                                  shiny::tags$b("This app takes a feature counts table, BAM files, and a sample info table and performs a basic DESeq analysis."),
                                  shiny::tags$br(),
                                  shiny::tags$br(),
                                  shiny::tags$b("The app has identified the samples:"),
                                  shiny::tags$br(),
                                  HTML(paste(sample_names, collapse = "<br/>")),
                                  shiny::tags$br(),
                                  shiny::tags$br(),
                                  shiny::tags$b("These data have been analyzed in the following pair-wise comparisons:"),
                                  shiny::tags$br(),
                                  HTML(paste(names(choices), collapse = "<br/>")),
                                  shiny::tags$br(),
                                  shiny::tags$br(),
                                  shiny::tags$b("These data have produced the following output tables:"),
                                  shiny::tags$br(),
                                  HTML(paste("/output_files/",comparisons,"_resOrdered.txt", sep = "", collapse = "<br/>")),
                                  shiny::tags$br(),
                                  shiny::tags$br(),
                                  HTML(paste("E-mail Markus at ", shiny::tags$i("nevil@email.unc.edu"), " if there are problems.", sep = ""))
                                  
                         ),
                         tabPanel(title = "Genome Browser",
                                  sidebarLayout(
                                    sidebarPanel(
                                      width = 3,
                                      fluidRow(
                                        column(12, align = "center",
                                               shiny::tags$b("Genome Browser"))),
                                      br(),
                                      textInput(inputId = "browser_gene_list",
                                                label = "Provide input gene",
                                                value = "Hunchback"),
                                      shiny::tags$p(shiny::tags$h6("If the correct gene is not retrieved, try using the 'FBgn' ID from FlyBase")),
                                      br(),
                                      selectInput(inputId = "track_type",
                                                  label = "Select track type",
                                                  choices = c("Coverage" = "coverage",
                                                              "Pileup" = "pileup")),
                                      checkboxInput(inputId = "sashimi",
                                                    label = "Include sashimi plots?"),
                                      numericInput(inputId ="plusStart",
                                                   label = "Display extra upstream",
                                                   value = "0"),
                                      numericInput(inputId ="plusEnd",
                                                   label = "Display extra downstream",
                                                   value = "0"),
                                      numericInput(inputId ="max_height",
                                                   label = "Maximum height of graphs ('Zero' will auto-select height for each graph)",
                                                   value = "0"),
                                      selectInput(inputId = "track_list",
                                                  label = "Choose samples to plot:",
                                                  choices = sample_table$baseName,
                                                  multiple = TRUE),
                                      hr(),
                                      textInput(inputId = "color_list",
                                                label = "Provide custom colors for each track"),
                                      shiny::tags$h6(shiny::tags$a(href = "http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf", "Click here for available colors")),
                                      actionButton(inputId = "browser_go",
                                                   label = "Update"),
                                      downloadButton("download_browser_plot", label = "Download")
                                    ),
                                    mainPanel(
                                      shiny::tags$head(shiny::tags$script('
                                                            var dimension = [0, 0];
                                                            $(document).on("shiny:connected", function(e) {
                                                            dimension[0] = window.innerWidth;
                                                            dimension[1] = window.innerHeight;
                                                            Shiny.onInputChange("dimension", dimension);
                                                            });
                                                            $(window).resize(function(e) {
                                                            dimension[0] = window.innerWidth;
                                                            dimension[1] = window.innerHeight;
                                                            Shiny.onInputChange("dimension", dimension);
                                                            });
                                                            ')),
                                      span(textOutput("colorText"), style="color:red"),
                                      plotOutput("browser", width = "auto", height = "auto")
                                    )
                                  )),
                         tabPanel(title = "Log2 Fold change heatmaps",
                                  sidebarLayout(
                                    sidebarPanel(
                                      width = 3,
                                      fluidRow(
                                        column(12, align = "center",
                                               shiny::tags$b("Log2 Fold Change heatmaps"))),
                                      br(),
                                      shiny::tags$p("Select a control comparison and plot a heatmap of all Log2 FC values for each comparison."),
                                      hr(),
                                      selectInput(inputId = "fc_heat_comp",
                                                  label = "What is the control",
                                                  choices = unique(sample_names),
                                                  multiple = FALSE),
                                      
                                      br(),
                                      checkboxInput(inputId = "fc_heat_sample_sort",
                                                    label = "Sort by cluster?( False means sort by sample)",
                                                    value = TRUE),
                                      conditionalPanel(
                                        condition = 'input.fc_heat_sample_sort == TRUE',
                                        sliderInput(inputId = "fc_cluster_n",
                                                    label = "Select the number of clusters:",
                                                    value = 5,
                                                    min = 1,
                                                    max = 12),
                                        selectInput(inputId = "fc_heat_method",
                                                    label = "Select clustering method:",
                                                    choices =
                                                      list("heirarchical", "kmeans"))
                                      ),
                                      conditionalPanel(
                                        condition = 'input.fc_heat_sample_sort == FALSE',
                                        selectInput(inputId = "fc_sort_sample_list",
                                                    label = "Sort by sample:",
                                                    choices = unique(sample_names),
                                                    multiple = FALSE)
                                      ),
                                      radioButtons(inputId = "fc_heat_diff_only",
                                                   label = "Filter genes that are differentially expressed?",
                                                   choices = c("Yes", "No"),
                                                   inline = TRUE),
                                      selectInput(inputId = "fc_heat_sample_list",
                                                  label = "Choose samples to plot:",
                                                  choices = choices,
                                                  multiple = TRUE),
                                      br(),
                                      actionButton(inputId = "fc_heat_go",
                                                   label = "Update")
                                      #downloadButton("download_fc_plot", "Download")
                                    ),
                                    mainPanel(
                                      plotOutput("fc_heat")
                                    )
                                  )),
                         tabPanel(title = "Clustering Analysis",
                                  sidebarLayout(
                                    sidebarPanel(
                                      width = 3,
                                      fluidRow(
                                        column(12, align = "center",
                                               shiny::tags$b("Gene Expression Clustering"))),
                                      br(),
                                      shiny::tags$p("Select clustering and filtering parameters to view a clustered heatmap of the data."),
                                      hr(),
                                      selectInput(inputId = "cluster_method",
                                                  label = "Select clustering method:",
                                                  choices =
                                                    list("heirarchical", "kmeans")),
                                      # selectInput(inputId = "cluster_by",
                                      #            label = "Select what values to cluster by:",
                                      #            choices =
                                      #              list("counts", "logFC")),
                                      sliderInput(inputId = "cluster_n",
                                                  label = "Select the number of clusters:",
                                                  value = 5,
                                                  min = 1,
                                                  max = 12),
                                      hr(),
                                      fluidRow(
                                        column(12, align = "center",
                                               shiny::tags$p(shiny::tags$b("Optional Filtering:")),
                                               shiny::tags$p("Only affects results when 'Filter genes that are differentially expressed?' is set to 'yes'.")
                                        )
                                      ),
                                      br(),
                                      radioButtons(inputId = "diff_only",
                                                   label = "Filter genes that are differentially expressed?",
                                                   choices = c("Yes", "No"),
                                                   inline = TRUE),
                                      # selectInput(inputId = "direction_of_change",
                                      #             label = "Expression change:",
                                      #             choices =
                                      #               list("Neither","Up", "Down")),
                                      # 
                                      selectInput(inputId = "cluster_sample_list",
                                                  label = "Differentially expressed in:",
                                                  choices = choices,
                                                  multiple = FALSE),
                                      
                                      actionButton(inputId = "cluster_go",
                                                   label = "Update")
                                    ),
                                    mainPanel(
                                      tabsetPanel(
                                        tabPanel("Heatmap",
                                                 #plot_ly("cluster"),
                                                 #plotly::plotlyOutput("cluster"),
                                                 plotOutput("cluster"),
                                                 #imageOutput("cluster"),
                                                 #renderPlotly("cluster",height ="auto")
                                                 downloadButton("download_count_hm_plot", "Download")),
                                        tabPanel("Gene Ontology",
                                                 textInput(inputId = "cluster_GO_list",
                                                           label = "Which cluster(s) should be used for Gene ontology?",
                                                           value = "1",
                                                           placeholder = "1"),
                                                 selectInput(inputId = "direction_of_change",
                                                             label = "Expression change:",
                                                             choices =
                                                               list("Neither","Up", "Down")),
                                                 textInput(inputId = "cluster_bg_list",
                                                           label = "Which cluster(s) should be used as the background? 'all' indicates all other genes (default).",
                                                           value = "all",
                                                           placeholder = "all"),
                                                 selectInput(inputId = "direction_bg_change",
                                                             label = "Expression change:",
                                                             choices =
                                                               list("Neither","Up", "Down")),
                                                 DT::dataTableOutput("cluster_ont")
                                                 #downloadButton("download_cluster_table", "Download")
                                        ),
                                        tabPanel("Boxplots",
                                                 textInput(inputId = "plot_which_list",
                                                           label = "Which clusters should be plotted?",
                                                           value = "1, 2",
                                                           placeholder = "1, 2"),
                                                 textInput(inputId = "box_order_list",
                                                           label = "Provide order of replicate groups*",
                                                           value = paste(sample_names, collapse = " "),
                                                           placeholder = paste(sample_names, collapse = " ")),
                                                 selectInput(inputId = "direction_of_change",
                                                             label = "Expression change:",
                                                             choices =
                                                               list("Neither","Up", "Down")),
                                                 
                                                 downloadButton("download_count_box_plot", "Download"),
                                                 plotOutput("cluster_box", height = "auto"),
                                                 br(),
                                                 fluidRow(
                                                   column(6,
                                                          shiny::tags$b("For 'counts' data, provide individual replicate groups:"),
                                                          br(),
                                                          HTML(paste(sample_names, collapse = "<br/>")),
                                                          br()
                                                   ),
                                                   column(6,
                                                          shiny::tags$b("For 'logFC' data, provide comparisons:"),
                                                          br(),
                                                          HTML(paste(names(choices), collapse = "<br/>"))
                                                   )
                                                 ),
                                                 
                                        ),
                                        tabPanel("Table",
                                                 DT::dataTableOutput("cluster_table"),
                                                 downloadButton("download_cluster_table", "Download"))
                                      )
                                    )
                                  )
                         ),
                         tabPanel("Gene Count Plots",
                                  sidebarLayout(
                                    sidebarPanel(
                                      width = 3,
                                      fluidRow(
                                        column(12, align = "center",
                                               shiny::tags$b("Log2 Normalized Count Plots"))),
                                      br(),
                                      shiny::tags$p("Select genes and samples to plot the mean Log2 normalized counts for that group of replicates."),
                                      shiny::tags$p("Plots will appear in the order they were listed, and samples will be plotted in the order they were selected."),
                                      hr(),
                                      textInput(inputId = "count_gene_list",
                                                label = "Provide gene(s) separated by commas",
                                                value = "CG11023"),
                                      shiny::tags$p("If the correct gene is not retrieved, try using the 'FBgn' ID from FlyBase"),
                                      br(),
                                      selectInput(inputId = "count_sample_list",
                                                  label = "Choose samples to plot:",
                                                  choices = repList2,
                                                  multiple = TRUE),
                                      radioButtons(inputId = "line_connect",
                                                   label = "Connect points with line?",
                                                   choices = c("Connect", "Do Not Connect"),
                                                   inline = TRUE),
                                      
                                      actionButton(inputId = "count_go",
                                                   label = "Update"),
                                      
                                      downloadButton("download_gene_exp_plot", "Download")
                                    ),
                                    mainPanel(
                                      # verbatimTextOutput("clientdataText"),
                                      plotOutput("count", height = "auto")))),
                         
                         tabPanel("Gene Fold Change Plots",
                                  sidebarLayout(
                                    sidebarPanel(
                                      width = 3,
                                      fluidRow(
                                        column(12, align = "center",
                                               shiny::tags$b("Fold Change plots"))),
                                      br(),
                                      shiny::tags$p("Select genes and 'Mutant vs Wildtype' comparisons to plot the Log2 Fold Change."),
                                      shiny::tags$p("The plots will appear in the order that they are selected."),
                                      hr(),
                                      textInput(inputId = "fc_gene_list",
                                                label = "Provide gene(s) separated by commas",
                                                value = "yellow"),
                                      shiny::tags$p("If the correct gene is not retrieved, try using the 'FBgn' ID from FlyBase"),
                                      br(),
                                      selectInput(inputId = "fc_sample_list",
                                                  label = "Choose samples to plot:",
                                                  choices = choices,
                                                  multiple = TRUE),
                                      actionButton(inputId = "fc_go",
                                                   label = "Update"),
                                      downloadButton("download_fc_plot", "Download")
                                    ),
                                    mainPanel(
                                      plotOutput("fc", height = "auto")
                                    )
                                  )
                         ),
                         tabPanel("Volcano Plots",
                                  sidebarLayout(
                                    sidebarPanel(
                                      width = 3,
                                      fluidRow(
                                        column(12, align = "center",
                                               shiny::tags$b("Volcano plots"))),
                                      br(),
                                      shiny::tags$p("Select sample comparisons to plot the log10 Adjusted p-value vs log2 Fold Change."),
                                      shiny::tags$p("User-defined significant genes will be highlighted in red."),
                                      shiny::tags$p("Insignificant genes are automatically displayed in grey."),
                                      hr(),
                                      shiny::tags$head(shiny::tags$script('
                                                            var dimension = [0, 0];
                                                            $(document).on("shiny:connected", function(e) {
                                                            dimension[0] = window.innerWidth;
                                                            dimension[1] = window.innerHeight;
                                                            Shiny.onInputChange("dimension", dimension);
                                                            });
                                                            $(window).resize(function(e) {
                                                            dimension[0] = window.innerWidth;
                                                            dimension[1] = window.innerHeight;
                                                            Shiny.onInputChange("dimension", dimension);
                                                            });
                                                            ')),
                                      selectInput(inputId = "volcano_sample_list",
                                                  label = "Choose samples to plot:",
                                                  choices = choices,
                                                  multiple = TRUE),
                                      hr(),
                                      shiny::tags$b("Plotting Options"),
                                      br(),
                                      textInput(inputId = "padj_n",
                                                label = "Set the adjusted p-value cut off:",
                                                value = "0.05"),
                                      sliderInput(inputId = "fc_n",
                                                  label = "Set Fold Change cut off:",
                                                  value = c(-1, 1),
                                                  min = -10,
                                                  max = 10,
                                                  step = 0.5),
                                      textInput(inputId = "vc_color_list",
                                                label = "Provide colors for up, down, and insignificant genes (separated by commas)",
                                                value = c("red, blue, grey80")),
                                      shiny::tags$h6(shiny::tags$a(href = "http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf", "Click here for available colors")),
                                      numericInput(inputId = "v_point_size",
                                                   label = "Plot point size:",
                                                   value = 1.5),
                                      checkboxGroupInput(inputId = "v_checks",
                                                         label = "Select options:",
                                                         choices = list("Show Legend" = "legend",
                                                                        "Add grids" = "grids",
                                                                        "Add p-value line" = "pv",
                                                                        "Add Fold Change line" = "fc"
                                                         )
                                      ),
                                      textInput(inputId = 'v_title',
                                                label = "Add custom title(s)",
                                                value = ""),
                                      numericInput(inputId = "top_n",
                                                   label = "Label Top Hits:",
                                                   value = 0),
                                      shiny::tags$p("or"),
                                      textInput(inputId = "v_gene_list",
                                                label = "Highlight specific genes:",
                                                value = ""),
                                      conditionalPanel(
                                        condition = 'input.top_n > 0 || input.v_gene_list != ""',
                                        numericInput(inputId = "v_gene_size",
                                                     label = "Gene name font size:",
                                                     value = 4) 
                                      ),
                                      
                                      
                                      # checkboxInput(inputId = "fishers",
                                      #               label = "Perform Fisher's Exact Test",
                                      #               value = FALSE),
                                      br(),
                                      actionButton(inputId = "volcano_go",
                                                   label = "Update"),
                                      downloadButton("download_volcano_plot", "Download")
                                      
                                    ),
                                    mainPanel(
                                      plotOutput("volcano", height = "auto")
                                    )
                                  )
                         ),
                         
                         tabPanel("MA Plots",
                                  sidebarLayout(
                                    sidebarPanel(
                                      width = 3,
                                      fluidRow(
                                        column(12, align = "center",
                                               shiny::tags$b("MA plots"))),
                                      br(),
                                      shiny::tags$p("Select sample comparisons to plot the mean of normalized counts vs log2 Fold Change."),
                                      shiny::tags$p("User-defined significant genes will be highlighted in red and blue."),
                                      shiny::tags$p("Insignificant genes are automatically displayed in grey."),
                                      hr(),
                                      shiny::tags$head(shiny::tags$script('
                                                            var dimension = [0, 0];
                                                            $(document).on("shiny:connected", function(e) {
                                                            dimension[0] = window.innerWidth;
                                                            dimension[1] = window.innerHeight;
                                                            Shiny.onInputChange("dimension", dimension);
                                                            });
                                                            $(window).resize(function(e) {
                                                            dimension[0] = window.innerWidth;
                                                            dimension[1] = window.innerHeight;
                                                            Shiny.onInputChange("dimension", dimension);
                                                            });
                                                            ')),
                                      selectInput(inputId = "ma_sample_list",
                                                  label = "Choose samples to plot:",
                                                  choices = choices,
                                                  multiple = TRUE),
                                      hr(),
                                      shiny::tags$b("Plotting Options"),
                                      br(),
                                      textInput(inputId = "ma_padj_n",
                                                label = "Set the adjusted p-value cut off:",
                                                value = "0.05"),
                                      sliderInput(inputId = "ma_fc_n",
                                                  label = "Set Fold Change cut off:",
                                                  value = c(-1, 1),
                                                  min = -10,
                                                  max = 10,
                                                  step = 0.5),
                                      textInput(inputId = "ma_color_list",
                                                label = "Provide colors for up, down, and insignificant genes (separated by commas)",
                                                value = c("red, blue, grey80")),
                                      shiny::tags$h6(shiny::tags$a(href = "http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf", "Click here for available colors")),
                                      numericInput(inputId = "ma_point_size",
                                                   label = "Plot point size:",
                                                   value = 1.5),
                                      checkboxGroupInput(inputId = "ma_checks",
                                                         label = "Select options:",
                                                         choices = list("Show Legend" = "legend",
                                                                        "Add grids" = "grids",
                                                                        "Add Fold Change line" = "fc"
                                                         )
                                      ),
                                      textInput(inputId = 'ma_title',
                                                label = "Add custom title(s)",
                                                value = ""),
                                      numericInput(inputId = "ma_top_n",
                                                   label = "Label Top Hits:",
                                                   value = 0),
                                      shiny::tags$p("or"),
                                      textInput(inputId = "ma_gene_list",
                                                label = "Highlight specific genes:",
                                                value = ""),
                                      conditionalPanel(
                                        condition = 'input.ma_top_n > 0 || input.ma_gene_list != ""',
                                        numericInput(inputId = "ma_gene_size",
                                                     label = "Gene name font size:",
                                                     value = 4)
                                      ),
                                      br(),
                                      actionButton(inputId = "ma_go",
                                                   label = "Update"),
                                      downloadButton("download_ma_plot", "Download")
                                      
                                    ),
                                    mainPanel(
                                      plotOutput("ma", height = "auto")
                                    )
                                  )
                         ),
                         tabPanel("Gene Ontology",
                                  shiny::tags$p("Find gene ontology enrichment for a set of differentially expressed genes."),
                                  shiny::tags$p("Select the list of genes and the directionality of change"),
                                  br(),
                                  shiny::tags$p("Note: you can provide a specific set of genes as the background by un-checking the box."),
                                  hr(),
                                  fluidRow(
                                    column(4,
                                           selectInput(inputId = "GO_sample_list",
                                                       label = "Differentially expressed in:",
                                                       choices = choices,
                                                       multiple = FALSE)
                                    ),
                                    column(4, offset = 0.5,
                                           selectInput(inputId = "direction_go",
                                                       label = "Expression change:",
                                                       choices =
                                                         list("Up", "Down", "Both")
                                           )
                                    ),
                                    column(4, offset = 0.5,
                                           checkboxInput(inputId = "go_all",
                                                         label = "Use all other genes as background? (default behavior)",
                                                         value = TRUE)
                                    )),
                                  fluidRow(
                                    column(4,
                                           selectInput(inputId = "GO_bg_list",
                                                       label = "Background gene set:",
                                                       choices = choices,
                                                       multiple = FALSE)
                                    ),
                                    column(4, offset = 0.5,
                                           selectInput(inputId = "bg_go",
                                                       label = "Expression change:",
                                                       choices =
                                                         list("Up", "Down", "Both")
                                           )
                                    )),
                                  fluidRow(
                                    column(5, offset = .75,
                                           actionButton(inputId = "GO_go",
                                                        label = "Update"
                                           )),
                                    column(5, offset = 1,
                                           downloadButton("download_go_table", "Download")
                                    )),
                                  fluidRow(
                                    column(12, align = "center",
                                           shiny::tags$h2(shiny::tags$b(textOutput("gene_ontology_text")))
                                    )),
                                  DT::dataTableOutput("gene_ont")
                         ),
                         
                         tabPanel("Raw Result tables",
                                  fluidRow(
                                    column(2,
                                           selectInput(inputId = "first_group",
                                                       label = "First group:",
                                                       choices = sample_names,
                                                       selected = sample_names[1]
                                           )),
                                    column(2, offset = 0.5,
                                           selectInput(inputId = "second_group",
                                                       label = "Second group:",
                                                       choices = sample_names,
                                                       selected = sample_names[2]
                                           )),
                                    column(2, offset = .75,
                                           actionButton(inputId = "table_go",
                                                        label = "Update"
                                           )),
                                    column(2, offset = 1,
                                           downloadButton("download_res_table", "Download")
                                    )),
                                  fluidRow(
                                    column(12, align = "center",
                                           shiny::tags$h2(shiny::tags$b(textOutput("res_table_text")))
                                    )),
                                  DT::dataTableOutput("res_table")
                         ),
                         
                         
                         ##### Fix so it isn't just using MA plot stuff ####
                         tabPanel("PCA plot",
                                  mainPanel(
                                    plotOutput("PCA")
                                  )
                         )
)

)

