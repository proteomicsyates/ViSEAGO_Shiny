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
                             column(width = 6,
                                    selectizeInput(inputId = "protein_column", label = "Protein columns (one per experiment):", choices = c(""), multiple = FALSE)

                             )
                         ),

                         fluidRow(
                             column(width = 12,
                                    "Select the column where the protein identifiers are in the input file:"
                             )
                         ),

                         br(),
                         fluidRow(
                             column(width = 6,
                                    selectizeInput(inputId = "pvalues_columns", label = "P-values columns:", choices = c(""), multiple = TRUE)

                             )
                         ),
                         fluidRow(
                             column(width = 12,
                                    HTML("Select the columns where the associated p-values are for each experiment to compare.<br>
                                     For each column selected, a new experiment will be considered in the comparison:")
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
                             column(width = 4,
                                    textInput("experiment_name", "Analysis name", placeholder = "my_experiment_name_here")
                             )
                         ),
                         fluidRow(
                             column(width = 12,
                                    "Set a name to identify this analysis"
                             )
                         ),
                         br(),
                         fluidRow(
                             column(width = 6,

                                    fluidRow(
                                        column(width = 6,
                                               selectInput(inputId = "go_type", label = "GO type", selected = "BP", choices = c("", "Biological Process"="BP", "Molecular Function"="MF", "Cellular Component" = "CC" ))
                                        )
                                    ),
                                    fluidRow(
                                        column(width = 12,
                                               "Select the GO annotations domain to use in the analysis"
                                        )
                                    )
                             ),
                             column(width = 6,

                                    fluidRow(
                                        column(width = 6,
                                               selectInput(inputId = "species", label = "Species", selected = "human", choices = c("", organisms ))
                                        )
                                    ),
                                    fluidRow(
                                        column(width = 12,
                                               "Select the species of the proteins"
                                        )
                                    )
                             )
                         ),
                         br(),
                         fluidRow(
                             column(width = 6,

                                    fluidRow(
                                        column(width = 6,
                                               numericInput(inputId = "pvalue_threshold", label = "P-value (or score) threshold:", value = 0.05)
                                        )
                                    ),
                                    fluidRow(
                                        column(width = 12,
                                               HTML("This is the p-value which will determine which proteins are considered by the analysis among all the list in the input file<br>
                                         The whole list of proteins (with no threshold) will be considered as the background:")
                                        )
                                    ),
                             ),
                             column(width = 6,

                                    fluidRow(
                                        column(width = 12,
                                               radioButtons(inputId = "use_ranked_list", label = "Type of analysis", choices = c("Regular Enrichment Analysis" = "regular",  "FSGEA" = "FSGEA"))
                                        )
                                    ),
                                    fluidRow(
                                        column(width = 12,
                                               "Regular Enrichment Analysis refers to a GO enrichment analysis using the proteins passing a certain p-value threshold."
                                        )
                                    ),
                                    fluidRow(
                                        column(width = 12,
                                               "With FGSEA (Fast Gene Set Enrichment Analysis) the p-values will be used to rank the proteins and perform the enrichment. In this case, the p-value threshold will be ignore.",
                                        )
                                    )
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

                                 fluidRow(column(width = 12,
                                                 "To fully explore the results of this functional analysis, a hierarchical clustering method is performed based on one of SS distances (i.e Wang) between the enriched GO terms and a chosen aggregation criteria (i.e ward.D2)."
                                 )),
                                 fluidRow(column(width = 12,
                                                 "A table of enriched GO terms in at least one of the comparison is displayed."
                                 )),
                                 br(),
                                 fluidRow(
                                     column(width = 12,
                                            div(dataTableOutput(outputId = "go_cluster_table"), style = "font-size: 75%; width: 75%")
                                     )
                                 ),
                                 fluidRow(column(width = 12,
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
                         fluidRow(column(width = 12,
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
                                 fluidRow(column(width = 12,
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







    ### add columns of the table for selectInput for proteins ----
    observe({
        req( df_columns())
        # by default, we select the ones containing 'acc' or 'protein'
        default_selection <- grep("acc", df_columns(), ignore.case=TRUE, value=TRUE)
        default_selection <- c(default_selection, grep("protein", df_columns(), ignore.case=TRUE, value=TRUE))
        updateSelectInput(inputId = "protein_column", choices = df_columns(), selected = unique(default_selection[1]))
    })
    ### add columns of the table for selectInput for pvalues or scores ----
    observe({
        # by default, we select the ones containing 'pvalue' or 'qvalue' or 'prob'
        default_selection <- grep("pvalue", df_columns(), ignore.case=TRUE, value=TRUE)
        default_selection <- c(default_selection, grep("qvalue", df_columns(), ignore.case=TRUE, value=TRUE))
        default_selection <- c(default_selection, grep("prob", df_columns(), ignore.case=TRUE, value=TRUE))
        updateSelectizeInput(inputId = "pvalues_columns", choices = df_columns(), selected = default_selection)
    })

    ## comparisons: names of experiments as the columns of the pvalues ----
    comparisons <- reactiveVal()
    experiments_key <- reactiveVal()
    experiment_name <- reactiveVal()
    ontology <- reactiveVal()
    species <- reactiveVal()
    perform_FGSEA <- reactiveVal(value = FALSE)

    # read input parameters ----
    read_input_params <- function(){
        print("reading input params")
        ### ontology reactive is set to GO type selected by user ----
        ontology(input$go_type)
        ### species ----
        species(input$species)
        ### experiment name reactive is set to the name set by the user, and replacing spaces and commas ----
        name <- input$experiment_name
        name <- str_replace_all(name, " ", "_")
        name <- str_replace_all(name, ",", ".")
        experiment_name(name)
        ### experiments key reactive is set as the indexes of the columns that are used as the protein and as the scores ----
        selected <- c(input$protein_column, input$pvalues_columns)
        exp_ids <- which(df_columns()  %in% selected)
        key <- paste(exp_ids, collapse = "_")
        print(paste("experiments selected: ", key))
        experiments_key(key)
        # set comparisons, that is the columns of experiments coming from pvalues columns
        comparisons(input$pvalues_columns)
        # perform FGSEA
        perform_FGSEA(input$use_ranked_list == 'FGSEA')
    }

    start_ok <- reactiveVal(value = FALSE)
    are_input_params_ok <- function(){
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
        return(all_is_ok)
    }

    output$all_is_ok <- reactive({
        start_ok()
    })

    # this has to be after previous statement
    # outputOptions(output, 'all_is_ok', suspendWhenHidden = FALSE)

    # START ----
    ## perform enrichment analysys either using cutoff by score or a fgsea, using ranking by score ----
    observeEvent(
        eventExpr = {
            input$start_analysis_button
        },
        handlerExpr = {
            # read input params and set the corresponding reactive variables
            read_input_params()


            print("button start pressed")
            print(paste("ontology:", ontology() ))
            all_is_ok <- are_input_params_ok()
            start_ok(all_is_ok)

            if(!all_is_ok){
                return (NULL)
            }


            # look for the enrichment result if present
            enrichment_file_name <- "enrichment_result.rds"
            if (perform_FGSEA()){
                enrichment_file_name <- "enrichment_result_FGSEA.rds"
            }

            enrichmentResultsfile <- get_rds_path(
                file_name = enrichment_file_name,
                ontology = ontology(),
                species = species(),
                experiment_name = experiment_name(),
                columns_keys = experiments_key()
            )
            if (file.exists(enrichmentResultsfile)){
                print(paste("previous enrichment found at", enrichmentResultsfile))
                result <- readRDS(file = enrichmentResultsfile)
            }else{
                print("previous enrichment not found, doing it now...")
                myGene2GO <- myGene2GO()
                if (!perform_FGSEA()){
                    background <- background()
                    result <- perform_enrichment(background, myGene2GO, ontology())
                }else{
                    result <- perform_fgsea(myGene2GO, ontology())
                }

                if (!is.null(result)){
                    saveRDS(result, file = enrichmentResultsfile)
                }
            }
            enrichmentResults(result)
            if (!is.null(result)) {
                showModal(
                    modalDialog(
                        title = "Enrichment done",
                        div("Now you can go to the results tab and look at the plots"),

                        div("Some plots will require extra computation time...be patient."),
                        easyClose = TRUE,
                        footer = NULL)
                )
            }
        }
    )

    ## load GO ----
    myGene2GO <- reactive({
        withProgress(message = paste("Loading GO to Uniprot mapping for", input$species, "species"), detail= "Please wait...",
                     expr = {
                         # 2. GO annotation of genes ----
                         # using global location since this file is unique by species and shared by different analysis
                         file <- get_rds_path(file_name = paste0("uniprot2GO_", input$species, ".rds"))
                         if(!file.exists(file)){
                             ## load GO annotations from Uniprot ----
                             ret <- ViSEAGO::annotate(
                                 input$species,
                                 Uniprot
                             )
                             saveRDS(ret, file)
                         }else{
                             ret <- readRDS(file)
                         }
                     }
        )
        ret
    })
    ## function perform_enrichment
    perform_enrichment <- function(background, gene2GO, ontology){
        tryCatch(
            {
                withProgress(

                    min = 0,
                    max = length(comparisons()) + 1,
                    message = paste0("Performing GO (", ontology, ") enrichment analysis"),
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

                            topGO_data_file <- get_rds_path(
                                file_name = paste0("enrichment_", comparison, ".rds" ),
                                ontology = ontology(),
                                species = species(),
                                experiment_name =  experiment_name(),
                                columns_keys =  experiments_key()
                            )
                            if (!file.exists(topGO_data_file)){
                                # 3. Functional GO enrichment ----
                                ## 3.1 GO enrichment tests ----
                                ### topGO: create topGOdata for BP ----


                                BP <- ViSEAGO::create_topGOdata(
                                    geneSel=selection,
                                    allGenes=background,
                                    gene2GO=gene2GO,
                                    ont=ontology,
                                    nodeSize=5
                                )


                                saveRDS(BP, file = topGO_data_file)
                            }else{
                                BP <- readRDS(file = topGO_data_file)
                            }
                            assign(paste0(ontology,"_",comparison), BP, envir = .GlobalEnv)
                            enrichment_test_result_file <- get_rds_path(
                                file_name = paste0("enrichment_result_", comparison, ".rds" ),
                                species = species(),
                                ontology = ontology(),
                                experiment_name =  experiment_name(),
                                columns_keys =  experiments_key()
                            )
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
            },
            error = function(e) {
                showModal(
                    modalDialog(title = "Some error occurred",
                                "Are you sure you selected the correct species for your data?",
                                easyClose = TRUE,
                                footer = NULL)
                )
            }
        )
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

                    enrichment_test_result_file <- get_rds_path(
                        file_name = paste0("fgsea_result_", comparison_tag, ".rds" ),
                        species = species(),
                        ontology = ontology(),
                        experiment_name =  experiment_name(),
                        columns_keys =  experiments_key())
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
            if(!is.null(results)) {
                table <- results@data # table
                write.table(table, con, sep = "\t")
            }
        },
        contentType = "text/csv"
    )
    ## show go count ----
    output$go_count_plot <- renderPlotly({
        results <- enrichmentResults()
        if(is.null(results)){
            return (NULL)
        }
        ViSEAGO::GOcount(results)
    })

    ## show upset ----
    output$upset_plot <- renderPlot({
        results <- enrichmentResults()
        if (is.null(results)) {
            return(NULL)
        }
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

    calculate_semantic_similarities <- function(enrichmentResults, gene2GO, distance_type, ontology_type){
        withProgress({

            file <- get_rds_path(
                file_name = paste0('ss_', distance_type, ".rds"),
                ontology =  ontology_type,
                species = species(),
                experiment_name = experiment_name(),
                columns_keys = experiments_key()
            )
            if (file.exists(file)){
                myGOs <- readRDS(file)
            }else{
                myGOs <- compute_semantic_similarity(distance = distance_type,
                                                     myGENE2GO = gene2GO,
                                                     BP_sResults = enrichmentResults)
                saveRDS(myGOs, file = file)
                return(myGOs)
            }
        },message = paste("Calculating enriched GO terms semantic distances using distance:",distance_type), detail = "Please wait some seconds...")

    }
    ## semantic similarities MD plot ----
    output$ss_md_plot <- renderPlotly({
        distance <- input$ss_distance
        if (is.null(enrichmentResults())) {
            return(NULL)
        }
        semantic_similarities <- calculate_semantic_similarities(enrichmentResults(), myGene2GO(), distance, ontology())
        withProgress({
            plot_path <- get_rds_path(
                file_name = paste0('mdsplot_', distance, ".rds"),
                species = species(),
                ontology =  ontology(),
                experiment_name =  experiment_name(),
                columns_keys =  experiments_key())
            if (file.exists(plot_path)){
                plot <- readRDS(plot_path)
            }else{
                plot <- ViSEAGO::MDSplot(semantic_similarities)
                saveRDS(plot, file = plot_path)
            }
            plot
            # ViSEAGO::MDSplot(ss_all())
        },message = "Loading MD plot", detail = "Please wait a second")
    })



    get_rds_path <- function(file_name, ontology = NULL, species = NULL, experiment_name = NULL, columns_keys = NULL){

        folder <- 'data'
        if (!is.null(experiment_name)) {
            folder <- paste0(folder, '/', experiment_name)
        }
        if (!is.null(species)) {
            folder <- paste0(folder, '/', species)
        }
        if (!is.null(ontology)) {
            folder <- paste0(folder, '/', ontology)
        }
        if(!is.null(columns_keys)) {
            folder <- paste0(folder, '/', columns_keys)
        }

        if(!dir.exists(folder)){
            dir.create(folder, recursive = TRUE, showWarnings = TRUE)
        }

        file <- paste0(folder, '/', file_name)

    }

    go_clusters_RV <- reactive({
        withProgress(
            message = "Clustering GO terms",
            detail = "Please wait some seconds...",
            {
                if (is.null(enrichmentResults())) {
                    return(NULL)
                }
                show_ic <- input$go_cluster_heatmap_show_ic
                show_labels <- input$go_cluster_heatmap_show_labels
                distance <- input$go_cluster_heatmap_distance
                aggregation_method <- input$go_cluster_heatmap_aggregation_method
                semantic_similarities <- calculate_semantic_similarities(enrichmentResults(), myGene2GO(), distance, ontology())

                cluster_file <- get_rds_path(
                    file_name = paste0('clustering_', show_ic, '_', show_labels, '_',distance, '_',aggregation_method, ".rds"),
                    species = species(),
                    ontology = ontology(),
                    experiment_name = experiment_name(),
                    columns_keys = experiments_key()
                )
                if(file.exists(cluster_file)){
                    clusters <- readRDS(file = cluster_file)
                }else{
                    clusters <-ViSEAGO::GOterms_heatmap(
                        semantic_similarities,
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
                         heatmap_file <- get_rds_path(
                             file_name = paste0('goterms_heatmap_', show_ic, '_', show_labels, '_',distance, '_',aggregation_method, ".rds"),
                             ontology = ontology(),
                             species = species(),
                             experiment_name = experiment_name(),
                             columns_keys =  experiments_key()
                         )
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

    compute_semantic_similarity <- function(distance, myGENE2GO, BP_sResults){
        # initialyse
        myGOs<-ViSEAGO::build_GO_SS(
            gene2GO=myGENE2GO,
            enrich_GO_terms=BP_sResults
        )

        # compute all available Semantic Similarity (SS) measures
        myGOs<-ViSEAGO::compute_SS_distances(
            myGOs,
            distance = distance
        )

        # if (distance == 'Wang'){
        #     ss_wang()
        # }else if (distance == 'Resnik'){
        #     ss_resnik()
        # }else if (distance == 'Rel'){
        #     ss_rel()
        # }else if (distance == 'Lin'){
        #     ss_lin()
        # }else if (distance == 'Jiang'){
        #     ss_jiang()
        # }
    }

    # create the table of the clusters
    clusters_table_obj <- reactive({
        clusters <- go_clusters_RV()
        req(clusters)
        show_ic <- input$go_cluster_heatmap_show_ic
        show_labels <- input$go_cluster_heatmap_show_labels
        distance <- input$go_cluster_heatmap_distance
        aggregation_method <- input$go_cluster_heatmap_aggregation_method
        cluster_table_file <- get_rds_path(
            file_name = paste0('goterms_heatmap_table_', show_ic, '_', show_labels, '_',distance, '_',aggregation_method, ".rds"),
            experiment_name = experiment_name(),
            columns_keys = experiments_key(),
            ontology = ontology(),
            species = species()
        )
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
        table <- table %>% dplyr::select(!ends_with("genes") & !ends_with("genes_symbol") & !matches("definition") & !ends_with("log10_pvalue"))
        # now we filter out some columns
        # table %>% dplyr::select(!ends_with("genes") & !ends_with("genes_symbol") & !matches("definition") & !ends_with("log10_pvalue"))
        datatable(table, options= list(rownames = FALSE) )
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


    output$clusters_distances_plot <- renderPlotly({
        clusters <- go_clusters_RV()
        req(clusters)
        distance <- input$go_cluster_similarities_distance
        withProgress(
            message = paste("Calculating distances between clusters of GO terms using", distance),
            detail = "Please wait some seconds...",
            {
                cluster_distances_file_name <- paste0("cluster_ss_", distance, ".rds")
                cluster_distances_file <- get_rds_path(
                    file_name = cluster_distances_file_name,
                    ontology = ontology(),
                    experiment_name = experiment_name(),
                    columns_keys = experiments_key()
                )
                if (!file.exists(cluster_distances_file)){
                    # calculate semantic similarities between GO clusters ----
                    distances<-ViSEAGO::compute_SS_distances(
                        clusters,
                        distance=distance
                    )
                    saveRDS(distances, file = cluster_distances_file)
                }else{
                    distances <- readRDS(file = cluster)
                }
                setProgress(message = "Loading plot")

                # MDS plot with semantic similarities between GO clusters ----
                ViSEAGO::MDSplot(
                    distances,
                    "GOclusters"
                )
            })


    })
}


# Run the application
shinyApp(ui = ui, server = server)
