#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("ViSEAGO")


list.of.packages <- c("shiny","readxl", "tibble", "data.table",
                      "DT", "UpSetR")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)){
    install.packages(new.packages)
}
library(shiny)
library(readxl)  # install.packages("readxl") or install.packages("tidyverse")
# library(plyr)
library(tibble)




# load genes background
library(data.table)
library(plotly)
library(stringr)
library(DT)
library(UpSetR)

## connect to Uniprot-GOA ----
Uniprot <- ViSEAGO::Uniprot2GO()
# Display table of available organisms with Uniprot
organisms <- ViSEAGO::available_organisms(Uniprot)
organisms <- organisms$x$data$`!annotation_set`
# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("GO analysis with VISEAGO"),
    tabsetPanel(
        # upload file panel ----
        tabPanel("Upload file",
                 # Sidebar layout with input and output definitions ----
                 sidebarLayout(
                     # Sidebar panel for inputs ----
                     sidebarPanel(
                         # Input: Select a file ----
                         fileInput("file1", "Choose text File",
                                   multiple = FALSE,
                                   accept = c("text/csv",
                                              "text/tsv",
                                              "text/comma-separated-values,text/plain",
                                              "text/tab-separated-values,text/plain",
                                              ".csv",
                                              ".tsv")
                         ),
                         "Input file is a text separated file. Select it here, then, select the appropriate options below to read it while you see the preview on the right.",

                         # Horizontal line ----
                         tags$hr(),

                         # Input: Checkbox if file has header ----
                         checkboxInput("header", "File has header", TRUE),

                         # Input: Select separator ----
                         radioButtons("sep", "Column separator",
                                      choices = c(Comma = ",",
                                                  Semicolon = ";",
                                                  Tab = "\t"),
                                      selected = "\t"),

                         # # Input: Select quotes ----
                         # radioButtons("quote", "Quote",
                         #              choices = c(None = "",
                         #                          "Double Quote" = '"',
                         #                          "Single Quote" = "'"),
                         #              selected = '"'),

                         # # Horizontal line ----
                         # tags$hr(),

                         # Input: Select number of rows to display ----
                         radioButtons("disp", "Display",
                                      choices = c("Only first rows" = "head",
                                                  "All rows" = "all"),
                                      selected = "head")
                     ),
                     # Main panel for displaying outputs ----
                     mainPanel(
                         br(),
                         fluidRow(
                             column(12,textOutput("text_guide")
                             )
                         ),
                         tags$br(),
                         fluidRow(

                             column(12,            # Output: Data file ----
                                    dataTableOutput("contents")
                             )
                         )
                     )
                 )
        ),
        # input data ----
        tabPanel(title = "Select input columns",
                 conditionalPanel(
                     condition = "output.has_file",

                     br(),
                     h4("Columns selection"),
                     wellPanel(
                         style = "background: white",
                         fluidRow(
                             column(width = 12,
                                    "Select the column where the protein identifiers are in the input file:"
                             )
                         ),
                         fluidRow(
                             column(width = 6,
                                    selectizeInput(inputId = "protein_column", label = "Protein columns (one per experiment):", choices = c(""), multiple = FALSE)

                             )
                         ),
                         fluidRow(
                             column(width = 12,
                                    HTML("Select the columns where the associated p-values are for each experiment to compare.<br>
                                     For each column selected, a new experiment will be considered in the comparison:")
                             )
                         ),

                         fluidRow(
                             column(width = 6,
                                    selectizeInput(inputId = "pvalues_columns", label = "P-values (or score) columns (one per experiment):", choices = c(""), multiple = TRUE)

                             )
                         ),
                         fluidRow(
                             column(width = 12,
                                    HTML("This is the p-value which will determine which proteins are considered by the analysis among all the list in the input file<br>
                                         The whole list of proteins (with no threshold) will be considered as the background:")
                             )
                         ),
                         fluidRow(
                             column(width = 6,
                                    numericInput(inputId = "pvalue_threshold", label = "P-value (or score) threshold:", value = 0.05)
                             )
                         )

                     ),

                 )
        ),
        # biological process ----
        tabPanel(title="Analysis",
                 conditionalPanel(
                     condition = "output.has_file",

                     br(),
                     h4("Analysis options"),
                     wellPanel(
                         style = "background: white",
                         fluidRow(
                             column(width = 12,
                                    "Set a name to identify this analysis"
                             )
                         ),

                         fluidRow(
                             column(width = 4,
                                    textInput("experiment_name", "Analysis name")
                             )
                         ),
                         fluidRow(
                             column(width = 12,
                                    "Select the GO annotations domain to use in the analysis"
                             )
                         ),
                         fluidRow(
                             column(width = 6,
                                    selectInput(inputId = "go_type", label = "GO type", selected = "BP", choices = c("", "Biological Process"="BP", "Molecular Function"="MF", "Cellular Component" = "CC" ))
                             )
                         ),
                         fluidRow(
                             column(width = 12,
                                    "Select the species of the proteins"
                             )
                         ),
                         fluidRow(
                             column(width = 6,
                                    selectInput(inputId = "species", label = "Species", selected = "human", choices = c("", organisms ))
                             )
                         ),
                         fluidRow(
                             column(width = 6,
                                    "Use the p-values (or scores) of the differentially expressed proteins in the FGSEA analysis or perform just the enrichment analysis with the DE proteins",
                                    checkboxInput(inputId = "use_ranked_list", label = "Use p-value", value = FALSE)
                             )
                         )
                     ),
                     actionButton("start_analysis_button", label = "Start analysis", icon = icon("play-circle"))
                 )
        ),
        tabPanel(title="Results",

                 # Sidebar with a slider input for number of bins
                 tabsetPanel(

                     # Show a plot of the generated distribution
                     tabPanel(
                         title="Enrichment table",
                         br(),
                         downloadButton(label = "Download enrichment table", outputId = "download_enrichment_table"),

                         br(),
                         fluidRow(
                             column(width = 12,
                                    div(dataTableOutput(outputId = "merged_table"), style = "font-size: 75%; width: 75%")
                             )
                         )

                     ),
                     tabPanel(
                         title="GO count", br(),
                         plotlyOutput(outputId = "go_count_plot", height = "600px")
                     ),
                     tabPanel(
                         title="Upset plot",
                         br(),
                         fluidRow(column(width = 10, offset = 1,
                                         "Visualizations of intersecting GO term sets."
                         )),
                         fluidRow(column(width = 10, offset = 1,
                                         plotOutput(outputId = "upset_plot", height = "800px")
                         ))
                     ),
                     tabPanel(
                         title="GO terms Semantic Similarities",
                         br(),
                         sidebarLayout(
                             sidebarPanel(
                                 selectInput(inputId = "ss_distance", label = "Distance", choices = c("Wang", "Resnik", "Lin", "Rel", "Jiang")),
                                 width = 2
                             ),

                             mainPanel(
                                 br(),
                                 fluidRow(column(width = 9, offset = 1,
                                                 "A Multi Dimensional Scale (MDS) plot provides a representation of distances among a set of enriched GO terms on the two first dimensions. Some patterns could appear at this time and could be investigated in an interactive mode."
                                 )),
                                 fluidRow(column(width = 9, offset = 1,
                                                 plotlyOutput(outputId = "ss_md_plot", height = "800px")
                                 )),
                                 width = 10
                             )
                         )
                     ),
                     tabPanel(
                         title = "Clustering heatmap of GO terms",
                         br(),
                         sidebarLayout(
                             sidebarPanel(
                                 br(),
                                 selectInput(inputId = "go_cluster_heatmap_distance", label = "distance", choices = c("Wang", "Resnik", "Rel", "Lin", "Jiang")),
                                 checkboxInput(inputId = "go_cluster_heatmap_show_ic", label = "Display the Information Content (IC)", value = TRUE),
                                 checkboxInput(inputId = "go_cluster_heatmap_show_labels", label = "Display the GO terms ticks on y axis", value = FALSE),
                                 selectInput(inputId = "go_cluster_heatmap_aggregation_method", label = "Aggregation method", choices = c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid")),
                                 downloadButton(outputId = "go_cluster_table_download", label = "Download table"),
                                 width = 2
                             ),
                             mainPanel(
                                 br(),

                                 fluidRow(column(width = 10, offset = 1,
                                                 "To fully explore the results of this functional analysis, a hierarchical clustering method is performed based on one of SS distances (i.e Wang) between the enriched GO terms and a chosen aggregation criteria (i.e ward.D2)."
                                 )),
                                 fluidRow(column(width = 11, offset = 1,
                                                 "A table of enriched GO terms in at least one of the comparison is displayed."
                                 )),
                                 fluidRow(
                                     column(width = 12,
                                            div(dataTableOutput(outputId = "go_cluster_table"), style = "font-size: 75%; width: 75%")
                                     )
                                 ),
                                 fluidRow(column(width = 10, offset = 1,
                                                 "Enriched GO terms are ranked in a dendrogram as a result of a hierarchical clustering and colored depending on their cluster assignation. Additional illustrations are displayed: the GO description of GO terms (trimmed if more than 50 characters), a heatmap of -log10(pvalue) of enrichment test for each comparison, and the Information Content (IC)."
                                 )),
                                 plotlyOutput(outputId = "go_cluster_heatmap_plot", height = "800px"),

                                 width = 10
                             )
                         )
                     ),
                     tabPanel(
                         title="Multi Dimensional Scaling of GO terms",

                         br(),
                         fluidRow(column(width = 11, offset = 1,
                                         "Multi Dimensional Scale (MDS) plot with the overlay of GO terms clusters. It is a way to check the coherence of GO terms clusters on the MDS plot."
                         )),
                         plotlyOutput(outputId = "multi_dimensional_scaling_plot", height = "800px")
                     ),
                     tabPanel(
                         title = "Similarity between GO clusters",
                         br(),
                         sidebarLayout(
                             sidebarPanel(
                                 br(),
                                 selectInput(inputId = "go_cluster_similarities_distance", label = "Distance for semantic similarities", choices = c("max", "avg","rcmax", "BMA")),
                                 width = 2
                             ),
                             mainPanel(
                                 br(),
                                 fluidRow(column(width = 11, offset = 1,
                                                 "A colored Multi Dimensional Scale (MDS) plot provides a representation of distances between the clusters of GO terms. Each circle represents a cluster of GO terms and its size depends on the number of GO terms that it contains. Clusters of GO terms that are close should share a functional coherence."    )),
                                 plotlyOutput(outputId = "clusters_distances_plot", height = "800px")
                             )
                         )
                     )
                 )
        )
    )
)





# Server logic ----
server <- function(input, output) {

    # has_file reactive variable is true when input$file1 is not null ----
    output$has_file <- reactive({
        !is.null(input$file1)
    })
    # this has to be after previous statement
    outputOptions(output, 'has_file', suspendWhenHidden = FALSE)

    # input data table reactive variable
    input_df <- reactiveVal()

    # receive file and plot table ----
    output$contents <- renderDataTable({
        # input$file1 will be NULL initially. After the user selects
        # and uploads a file, head of that data file by default,
        # or all rows if selected, will be shown.

        req(input$file1)

        # when reading semicolon separated files,
        # having a comma separator causes `read.csv` to error
        tryCatch(
            {
                df <- read.csv(input$file1$datapath,
                               header = input$header,
                               sep = input$sep )
                input_df(df)
            },
            error = function(e) {
                # return a safeError if a parsing error occurs
                stop(safeError(e))
            }
        )
        table <- df
        if(input$disp == "head") {
            table <- head(df)
        }
        datatable(table, rownames = FALSE)

    })




    output$text_guide <- renderText({
        req(input_df())
        paste("Go to next tab ('Select columns to analyze') to choose the columns where the proteins are")
    }    )
    # 1. Genes of interest ----
    df_columns <- reactiveVal()
    ## comparisons: names of experiments as the columns of the pvalues ----
    comparisons <- reactiveVal()

    observeEvent(
        eventExpr = {
            input_df()
        },
        handlerExpr = {
            req(input_df())
            columns <- colnames(input_df())
            df_columns(columns)
        }
    )
    ## background ----
    background <- reactive({
        # get proteins from protein_column\
        table <- input_df()
        proteins <- table %>% dplyr::pull(input$protein_column)
    })

    ## list with all the data ----
    list_all <- reactiveVal()
    ## enrichment results
    enrichmentResults <- reactiveVal()

    can_start <- reactiveVal()

    ## listen to button to check input params to set start_analysis variable to TRUE or FALSE----
    observeEvent(
        eventExpr = {input$start_analysis_button},
        handlerExpr = {
            all_is_ok <- TRUE
            if (experiment_name() == ""){
                all_is_ok <- FALSE
                showNotification("An analysis name is required", type = "error")
            }
            if (is.null(input$protein_column)){
                all_is_ok <- FALSE
                showNotification("You must specify the protein identifier column", type = "error")
            }
            if (is.null(input$pvalues_columns)){
                all_is_ok <- FALSE
                showNotification("You must specify at least one column where the p-values or scores are", type = "error")
            }
            can_start(all_is_ok)
            browser()
        }
    )

    ## add columns for columns of proteins ----
    observe({
        req( df_columns())
        # by default, we select the ones containing 'acc' or 'protein'
        default_selection <- grep("acc", df_columns(), ignore.case=TRUE, value=TRUE)
        default_selection <- c(default_selection, grep("protein", df_columns(), ignore.case=TRUE, value=TRUE))
        updateSelectInput(inputId = "protein_column", choices = df_columns(), selected = unique(default_selection[1]))
    })
    ## add columns for columns of proteins ----
    observe({
        # by default, we select the ones containing 'pvalue' or 'qvalue' or 'prob'
        default_selection <- grep("pvalue", df_columns(), ignore.case=TRUE, value=TRUE)
        default_selection <- c(default_selection, grep("qvalue", df_columns(), ignore.case=TRUE, value=TRUE))
        default_selection <- c(default_selection, grep("prob", df_columns(), ignore.case=TRUE, value=TRUE))
        updateSelectizeInput(inputId = "pvalues_columns", choices = df_columns(), selected = default_selection)
    })
    ## experiments key
    experiments_key <- reactive({
        selected <- input$pvalues_columns
        exp_ids <- which(df_columns()  %in% selected)
        key <- paste(exp_ids, collapse = "_")
        print(paste("experiments selected: ", key))
        key
    })
    ## experiment name
    experiment_name <- reactive({
        name <- input$experiment_name
        name <- str_replace_all(name, " ", "_")
        name <- str_replace_all(name, ",", ".")
        name
    })

    ## ontology ----
    ontology <- reactive({
        ontology = input$go_type
    })

    ## enrichment results ----
    observeEvent(
        eventExpr = {
            input$start_analysis_button
        },
        handlerExpr = {
    browser()
            if(!can_start()){
                return (NULL)
            }

            # set comparisons
            comparisons(input$pvalues_columns)
            # look for the enrichment result if present
            enrichmentResultsfile <- get_rds_path("enrichment_result.rds")
            if (file.exists(enrichmentResultsfile)){
                result <- readRDS(file = enrichmentResultsfile)
            }else{
                myGene2GO <- myGene2GO()
                if (!input$use_ranked_list){
                    background <- background()
                    result <- perform_enrichment(background, myGene2GO, ontology())
                }else{
                    result <- perform_fgsea(myGene2GO, ontology())
                }
                browser()
                saveRDS(result, file = enrichmentResultsfile)
            }
            enrichmentResults(result)
        }
    )

    ## load GO ----
    myGene2GO <- reactive({
        withProgress(message = paste("Loading GO to Uniprot mapping for", input$species, "species"), detail= "Please wait...",
                     expr = {
                         # 2. GO annotation of genes ----

                         ## load GO annotations from Uniprot ----
                         ret <- ViSEAGO::annotate(
                             input$species,
                             Uniprot
                         )
                     }
        )
        ret
    })
    ## function perform_enrichment
    perform_enrichment <- function(background, gene2GO, ontology){
        withProgress(

            min = 1,
            max = length(comparisons())+1,
            message = "Performing enrichment analysis",
            detail = "This may take a minute for the first time",
            {

                i <- 1
                for(comparison in comparisons()){
                    incProgress(amount = 1, message = paste("Enrichment analysis of", comparison))

                    comparison_table <- input_df() %>% dplyr::select(input$protein_column, comparison)
                    comparison_table <- comparison_table %>% dplyr::filter(!is.na(comparison))

                    # comparison_table$NORM_PVALUE_1 <- as.numeric(comparison_table$NORM_PVALUE_1)
                    # table <- comparison[comparison$NORM_PVALUE_1<0.05,c("ACCESSION","NORM_PVALUE_1")]
                    # data.table::setorder(selection, "NORM_PVALUE_1") # order by pvalue

                    # load genes selection
                    selection <- comparison_table %>% dplyr::filter(get(comparison) <= as.numeric(input$pvalue_threshold))
                    selection <- selection %>% dplyr::pull(input$protein_column)

                    topGO_data_file <- get_rds_path(paste0("topGO_", comparison, ".rds" ),FALSE)
                    if (!file.exists(topGO_data_file)){
                        # 3. Functional GO enrichment ----
                        ## 3.1 GO enrichment tests ----
                        ### topGO: create topGOdata for BP ----
                        tryCatch(
                            {
                                BP <- ViSEAGO::create_topGOdata(
                                    geneSel=selection,
                                    allGenes=background,
                                    gene2GO=gene2GO,
                                    ont=ontology,
                                    nodeSize=5
                                )
                            },
                            error = function(e) {
                                showModal(
                                    modalDialog(title = "Some error occurred",
                                                "Are you sure you selected the correct species for your data?",
                                                easyClose = TRUE,
                                                footer = NULL)
                                )
                                return (NULL)
                            }
                        )

                        saveRDS(BP, file = topGO_data_file)
                    }else{
                        BP <- readRDS(file = topGO_data_file)
                    }
                    assign(paste0(ontology,"_",comparison), BP, envir = .GlobalEnv)

                    enrichment_test_result_file <- get_rds_path(paste0("enrichment_result_", comparison, ".rds" ),FALSE)
                    if (!file.exists(enrichment_test_result_file)){
                        ### perform TopGO test using clasic algorithm ----
                        classic<-topGO::runTest(
                            BP,
                            algorithm ="classic",
                            statistic = "fisher",
                            cutoff=0.01
                        )
                        saveRDS(classic, file = enrichment_test_result_file)
                    }else{
                        classic <- readRDS(file = enrichment_test_result_file)
                    }
                    assign(paste0("classic_",comparison), classic, envir = .GlobalEnv)

                    i <- i + 1
                }
                browser()
                input_list <- list()
                for(comparison in comparisons()){
                    input_list[[comparison]] <- c(paste0(ontology,"_",comparison), paste0("classic_",comparison))
                }
                incProgress(amount = 1, message = "Merging all enrichments")
                ViSEAGO::merge_enrich_terms(
                    Input=input_list,
                    cutoff = 0.01,
                    envir = .GlobalEnv
                )
            })
    }
    ## function perform_enrichment
    perform_fgsea <- function(gene2GO, ontology){
        withProgress(
            min = 1,
            max = length(comparisons())+1,
            message = "Performing enrichment analysis",
            detail = "This may take a minute for the first time",
            {
                i <- 1
                for(comparison in comparisons()){
                    comparison_tag <- df_columns()[i]

                    incProgress(amount = 1, message = paste("fgsea analysis of", comparison_tag))
                    comparison_table <- list_all()[[comparison]]
                    comparison_table$NORM_PVALUE_1 <- as.numeric(comparison_table$NORM_PVALUE_1)
                    table <- comparison_table[comparison_table$NORM_PVALUE_1<0.05,c("ACCESSION","NORM_PVALUE_1")]
                    data.table::setorder(table, "NORM_PVALUE_1") # order by pvalue

                    enrichment_test_result_file <- get_rds_path(paste0("fgsea_result_", comparison_tag, ".rds" ), FALSE)
                    if (!file.exists(enrichment_test_result_file)){

                        BP<-ViSEAGO::runfgsea(
                            geneSel=table,
                            ont=ontology,
                            gene2GO=gene2GO,
                            method ="fgseaMultilevel",
                            params = list(
                                scoreType = "pos",
                                minSize=5
                            )
                        )
                        assign(paste0(ontology,"_",comparison_tag), BP)


                        saveRDS(BP, file = enrichment_test_result_file)
                    }else{
                        BP <- readRDS(file = enrichment_test_result_file)
                    }


                    i <- i + 1
                }

                input_list <- list()
                for(comparison_tag in df_columns()){
                    input_list[[comparison_tag]] <- paste0(ontology,"_",comparison_tag)
                }
                incProgress(amount = 1, message = "Merging all enrichments")
                ViSEAGO::merge_enrich_terms(
                    Input=input_list
                )
            })
    }
    ## semantic distances
    ss_all <- reactiveVal()
    ss_all(readRDS(file = get_rds_path("ss_all.rds")))
    ss_wang <- reactiveVal()
    ss_wang(readRDS(file = get_rds_path("ss_wang.rds")))
    ss_resnik <- reactiveVal()
    ss_resnik(readRDS(file = get_rds_path("ss_resnik.rds")))
    ss_rel <- reactiveVal()
    ss_rel(readRDS(file = get_rds_path("ss_rel.rds")))
    ss_lin <- reactiveVal()
    ss_lin(readRDS(file = get_rds_path("ss_lin.rds")))
    ss_jiang <- reactiveVal()
    ss_jiang(readRDS(file = get_rds_path("ss_jiang.rds")))

    ## Show enrichment table ----
    output$merged_table <- renderDT({
        results <- enrichmentResults()
        req(results)
        table <- results@data # table
        # now we filter out some columns
        table %>% dplyr::select(!ends_with("genes") & !ends_with("genes_symbol") & !matches("definition") & !ends_with("log10_pvalue"))
    },
    filter = "top",
    options = list(
        pageLength = 20
    ))

    ## download table ----
    output$download_enrichment_table <- downloadHandler(
        filename = function(){
            "enrichment_table.txt"
        },
        content = function(con){
            results <- enrichmentResults()
            table <- results@data # table
            write.table(table, con, sep = "\t")
        },
        contentType = "text/csv"
    )
    ## show go count ----
    output$go_count_plot <- renderPlotly({
        results <- enrichmentResults()
        ViSEAGO::GOcount(results)
    })

    ## show upset ----
    output$upset_plot <- renderPlot({
        results <- enrichmentResults()
        tmp <- as_tibble(results@data)
        tmp <- tmp %>% dplyr::select(ends_with(".pvalue"))
        names(tmp) <- str_replace(names(tmp), ".pvalue", "")
        tmp2 <- as.data.frame(ifelse(tmp<0.01, 1, 0))
        plot <- UpSetR::upset(tmp2, nsets = ncol(tmp2),
                              mainbar.y.label = "Significant GO term intersections",
                              sets.x.label = "Number of GO termsin each experiment", text.scale = c(2, 1.2, 2, 2, 2, 2)
        )
        plot


    })

    ## semantic similarities MD plot ----
    output$ss_md_plot <- renderPlotly({
        distance <- input$ss_distance
        withProgress({
    browser()
            file <- get_rds_path(paste0('ss_', distance,".rds"))
            if (file.exists(file)){
                myGOs <- readRDS(file)
            }else{
                myGOs <- get_by_distance(distance)
                saveRDS(myGOs, file = file)
            }
            plot_path <- get_rds_path(paste0('mdsplot_', distance,".rds"))
            if (file.exists(plot_path)){
                plot <- readRDS(plot_path)
            }else{
                plot <- ViSEAGO::MDSplot(myGOs)
                saveRDS(plot, file = plot_path)
            }
            plot
            # ViSEAGO::MDSplot(ss_all())
        },message = paste("Loading MD plot using distance",distance), detail = "Please wait a second")
    })



    get_rds_path <- function(file_name, use_experiments_key = TRUE){
        # path using GO ontology
        ont <- ontology()
        folder <- paste0('data/', experiment_name(), '/', ont,'/', experiments_key())
        if(!dir.exists(folder)){
            dir.create(folder, recursive = TRUE, showWarnings = TRUE)
        }
        if (!is.null(use_experiments_key) & use_experiments_key){
            file <- paste0('data/', experiment_name(), '/', ont, '/', experiments_key(),'/',file_name)
        }else{
            file <- paste0('data/', experiment_name(), '/', ont, '/', file_name)
        }
    }

    go_clusters_RV <- reactive({
        withProgress(
            message = "Clustering GO terms",
            detail = "Please wait...",
            {
                show_ic <- input$go_cluster_heatmap_show_ic
                show_labels <- input$go_cluster_heatmap_show_labels
                distance <- input$go_cluster_heatmap_distance
                aggregation_method <- input$go_cluster_heatmap_aggregation_method
                myGOs <- get_by_distance(distance)
                cluster_file <- get_rds_path(paste0('clustering_', show_ic, '_', show_labels, '_',distance, '_',aggregation_method, ".rds"))
                if(file.exists(cluster_file)){
                    clusters <- readRDS(file = cluster_file)
                }else{
                    clusters <-ViSEAGO::GOterms_heatmap(
                        myGOs,
                        showIC=show_ic,
                        showGOlabels=show_labels,
                        GO.tree=list(
                            tree=list(
                                distance=distance,
                                aggreg.method=aggregation_method
                            ),
                            cut=list(
                                dynamic=list(
                                    pamStage=TRUE,
                                    pamRespectsDendro=TRUE,
                                    deepSplit=2,
                                    minClusterSize =2
                                )
                            )
                        ),
                        samples.tree=NULL
                    )
                    saveRDS(clusters, file = cluster_file)
                }
                clusters
            })
    })
    # clustering heatmap of GO terms
    output$go_cluster_heatmap_plot <- renderPlotly({
        clusters <-  go_clusters_RV()
        req(clusters)
        withProgress(message = "Showing heatmap",
                     detail = "Please wait...",
                     {
                         show_ic <- input$go_cluster_heatmap_show_ic
                         show_labels <- input$go_cluster_heatmap_show_labels
                         distance <- input$go_cluster_heatmap_distance
                         aggregation_method <- input$go_cluster_heatmap_aggregation_method
                         heatmap_file <- get_rds_path(paste0('goterms_heatmap_', show_ic, '_', show_labels, '_',distance, '_',aggregation_method, ".rds"))
                         if (file.exists(heatmap_file)){
                             heatmap <- readRDS(file = heatmap_file)
                         }else{
                             heatmap <- ViSEAGO::show_heatmap(
                                 clusters,
                                 "GOterms"
                             )
                             saveRDS(heatmap, file = heatmap_file)
                         }
                         heatmap
                     }
        )
    })

    get_by_distance <- function(distance){
        if (distance == 'Wang'){
            ss_wang()
        }else if (distance == 'Resnik'){
            ss_resnik()
        }else if (distance == 'Rel'){
            ss_rel()
        }else if (distance == 'Lin'){
            ss_lin()
        }else if (distance == 'Jiang'){
            ss_jiang()
        }
    }

    # create the table of the clusters
    clusters_table_obj <- reactive({
        clusters <- go_clusters_RV()
        req(clusters)
        show_ic <- input$go_cluster_heatmap_show_ic
        show_labels <- input$go_cluster_heatmap_show_labels
        distance <- input$go_cluster_heatmap_distance
        aggregation_method <- input$go_cluster_heatmap_aggregation_method
        cluster_table_file <- get_rds_path(paste0('goterms_heatmap_table_', show_ic, '_', show_labels, '_',distance, '_',aggregation_method, ".rds"))
        if (file.exists(cluster_table_file)){
            table_obj <- readRDS(file = cluster_table_file)
        }else{
            table_obj <- ViSEAGO::show_table(clusters, file = NULL)$x$data
            saveRDS(table_obj, file = cluster_table_file)
        }
        table_obj
    })

    # render table of GO clusters----
    output$go_cluster_table <- renderDT({
        table_obj <- clusters_table_obj()

        req(table_obj)
        table <- table_obj
        # now we filter out some columns
        table %>% dplyr::select(!ends_with("genes") & !ends_with("genes_symbol") & !matches("definition") & !ends_with("log10_pvalue"))
        # now we filter out some columns
        # table %>% dplyr::select(!ends_with("genes") & !ends_with("genes_symbol") & !matches("definition") & !ends_with("log10_pvalue"))
    },
    filter = "top",
    options = list(
        pageLength = 5,
        rownames = FALSE
    ))

    # download table of clusters of GO terms
    output$go_cluster_table_download <- downloadHandler(
        filename = function(){
            "GO_clusters.txt"
        },
        content = function(file){
            table_obj <- clusters_table_obj()
            table <- table_obj
            write.table(table, file, sep = "\t")
        },
        contentType = "text/csv"
    )


    # semantic similarities MD plot
    output$multi_dimensional_scaling_plot <- renderPlotly({
        clusters <- go_clusters_RV()
        req(clusters)
        withProgress({
            ViSEAGO::MDSplot(
                clusters,
                "GOterms"
            )
        },message = "Loading plot", detail = "Please wait a second")
    })

    # calculate semantic similarities between GO clusters
    output$clusters_distances_plot <- renderPlotly({
        clusters <- go_clusters_RV()
        req(clusters)
        distance <- input$go_cluster_similarities_distance
        withProgress(message = "Calculating distances between clusters of GO terms", detail = "Please wait a second",
                     {
                         distances<-ViSEAGO::compute_SS_distances(
                             clusters,
                             distance=distance
                         )
                         setProgress(message = "Loading plot")
                         ViSEAGO::MDSplot(
                             distances,
                             "GOclusters"
                         )
                     })


    })
}


# Run the application
shinyApp(ui = ui, server = server)
