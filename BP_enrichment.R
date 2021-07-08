#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
library(shiny)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
# BiocManager::install("ViSEAGO")

library(readxl)  # install.packages("readxl") or install.packages("tidyverse")
library(plyr)
library(tibble)

# 1. Genes of interest ----
# load excel file ----
setwd("C:\\Users\\salvador\\Desktop\\Jolene\\KLABC-15compare_GOAnalysis")
excel_file <- "KLABC-15compare.xlsx"
comparisons <- excel_sheets(excel_file)
comparisons_tags <- sapply(comparisons, function(x) {
  ret <- gsub("\\s+", "_",x)
  ret <- gsub(",", "_",ret)
})
list_all <- lapply(comparisons, function(x) {
  read_excel(path = excel_file, sheet = x)
}
)
names(list_all) <- comparisons

# get all the proteins as background
accs <- c()
for(i in seq(1,length(list_all))){
  accs <- c(accs,list_all[[i]]$ACCESSION)
}
accs <- unique(accs)
background <- as.character(accs)



# load genes background
library(data.table)
library(plotly)
library(stringr)


# 2. GO annotation of genes ----
## connect to Uniprot-GOA ----
Uniprot<-ViSEAGO::Uniprot2GO()
# Display table of available organisms with Uniprot
ViSEAGO::available_organisms(Uniprot)
## load GO annotations from Uniprot ----
myGENE2GO<-ViSEAGO::annotate(
  "mouse",
  Uniprot
)

comparison <- comparisons[1]

i <- 1
for(comparison in comparisons){
  comparison_tag <- comparisons_tags[i]


  comparison_table <- list_all[[comparison]]
  print(paste(comparison, dim(comparison_table)))

  comparison_table$NORM_PVALUE_1 <- as.numeric(comparison_table$NORM_PVALUE_1)
  # table <- comparison[comparison$NORM_PVALUE_1<0.05,c("ACCESSION","NORM_PVALUE_1")]
  # data.table::setorder(selection, "NORM_PVALUE_1") # order by pvalue

  # load genes selection
  selection <- comparison_table[comparison_table$NORM_PVALUE_1<0.05,]$ACCESSION

  # 3. Functional GO enrichment ----
  ## 3.1 GO enrichment tests ----
  ### topGO: create topGOdata for BP ----
  BP <- ViSEAGO::create_topGOdata(
    geneSel=selection,
    allGenes=background,
    gene2GO=myGENE2GO,
    ont="BP",
    nodeSize=5
  )
  assign(paste0("BP_",comparison_tag), BP)

  ### perform TopGO test using clasic algorithm ----
  classic<-topGO::runTest(
    BP,
    algorithm ="classic",
    statistic = "fisher",
    cutoff=0.01
  )
  assign(paste0("classic_",comparison_tag), classic)

  i <- i + 1
}
## fgsea: perform fgseaMultilevel tests
# BP<-ViSEAGO::runfgsea(
#     geneSel=table,
#     ont="BP",
#     gene2GO=myGENE2GO,
#     method ="fgseaMultilevel",
#     params = list(
#         scoreType = "pos",
#         minSize=5
#     )
# )

## 3.2 Combine enriched GO terms----
### merge results from topGO ----
# build the list of comparisons
input_list <- list()
for(comparison_tag in comparisons_tags){
  input_list[[comparison_tag]] <- c(paste0("BP_",comparison_tag), paste0("classic_",comparison_tag))
}
BP_sResults<-ViSEAGO::merge_enrich_terms(
  Input=input_list
)

### Salva: save data in file
saveRDS(BP_sResults, file = "BP_sResults.rds")

# Start here if enrichment is already done ----
BP_sResults <- readRDS(file = "BP_sResults.rds")



# merge results from fgsea ----
# BP_sResults<-ViSEAGO::merge_enrich_terms(
#     Input=list(
#         condition="BP"
#     )
# )

### print the merged table in a file ----
ViSEAGO::show_table(
  BP_sResults,
  "BP_sResults.xls"
)
r <- BP_sResults@data # table

## 3.3 Graphs of GO enrichment tests ----
### count significant (or not) pvalues by condition ----
p <- ViSEAGO::GOcount(BP_sResults)
p
### display interactions
ViSEAGO::Upset(
  BP_sResults,
  file="OLexport.xls"
)
# make my own upset
library(UpSetR)
movies <- read.csv(system.file("extdata", "movies.csv", package = "UpSetR"),
                   header = T, sep = ";")

tmp <- as.tibble(BP_sResults@data)
tmp <- tmp %>% select(ends_with(".pvalue"))
tmp2 <- as.data.frame(ifelse(tmp<0.01, 1, 0))
plot <- UpSetR::upset(tmp2, nsets = ncol(tmp2))
plot
# 4. GO terms semantic similarity----
## Initialyse----
ss<-ViSEAGO::build_GO_SS(
  gene2GO=myGENE2GO,
  enrich_GO_terms=BP_sResults
)
saveRDS(ss, file = "ss.rds")
## compute all available Semantic Similarity (SS) measures----
ss_all<-ViSEAGO::compute_SS_distances(
  ss,
  distance=c("Resnik","Lin","Rel","Jiang","Wang")
)
saveRDS(ss_all, file = "ss_all.rds")


# 5. Visualization and interpretation of enriched GO terms ----
## 5.1 Multi dimensional scaling of GO terms ----
### display MDSplot ----
ViSEAGO::MDSplot(myGOs)

### print MDSplot ----
ViSEAGO::MDSplot(
  myGOs,
  file="mdsplot1.png"
)
## 5.2 Clustering heatmap of GO terms----
### GOterms heatmap with the default parameters----
Wang_clusters_wardD2<-ViSEAGO::GOterms_heatmap(
  myGOs,
  showIC=TRUE,
  showGOlabels=TRUE,
  GO.tree=list(
    tree=list(
      distance="Wang",
      aggreg.method="ward.D2"
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
### Display the clusters-heatmap ----
ViSEAGO::show_heatmap(
  Wang_clusters_wardD2,
  "GOterms"
)

### print the clusters-heatmap ----
ViSEAGO::show_heatmap(
  Wang_clusters_wardD2,
  "GOterms",
  "cluster_heatmap_Wang_wardD2.png"
)

### Display the clusters-heatmap table ----
t <- ViSEAGO::show_table(Wang_clusters_wardD2)
t <- t$x$data # to get the data table
### Print the clusters-heatmap table ----
ViSEAGO::show_table(
  Wang_clusters_wardD2,
  "cluster_heatmap_Wang_wardD2.xls"
)

## 5.3 Multidimensional scaling of GO terms ----
### display colored MDSplot ----
ViSEAGO::MDSplot(
  Wang_clusters_wardD2,
  "GOterms"
)

### print colored MDSplot ----
ViSEAGO::MDSplot(
  Wang_clusters_wardD2,
  "GOterms",
  file="mdsplot2.png"
)

# 6. Visualization and interpretation of GO clusters ----
## 6.1 Compute semantic similarity between GO clusters ----
### calculate semantic similarites between clusters of GO terms ----
Wang_clusters_wardD2<-ViSEAGO::compute_SS_distances(
  Wang_clusters_wardD2,
  distance=c("max", "avg","rcmax", "BMA")
)
### build and highlight in an interactive MDSplot grouped clusters for one distance object ----
ViSEAGO::MDSplot(
  Wang_clusters_wardD2,
  "GOclusters"
)

### build and highlight in MDSplot grouped clusters for one distance object ----
ViSEAGO::MDSplot(
  Wang_clusters_wardD2,
  "GOclusters",
  file="mdsplot3.png"
)

## 6.2 GO clusters semantic similarities heatmap ----
### GOclusters heatmap ----
Wang_clusters_wardD2<-ViSEAGO::GOclusters_heatmap(
  Wang_clusters_wardD2,
  tree=list(
    distance="BMA",
    aggreg.method="ward.D2"
  )
)
### display the GOClusters heatmap ----
ViSEAGO::show_heatmap(
  Wang_clusters_wardD2,
  "GOclusters"
)

### print the GOClusters heatmap in a file ----
ViSEAGO::show_heatmap(
  Wang_clusters_wardD2,
  "GOclusters",
  "Wang_clusters_wardD2_heatmap_groups.png"
)