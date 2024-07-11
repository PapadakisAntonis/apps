library(shiny)
library(htmltools)
library(topGO)
library(ggplot2)
library('psupertime')
library('Seurat')
library('Signac')

ui <- fluidPage(
  titlePanel("Psupertime analysis"),
  "In this app, users can run of psupertime, ",
  "an R package which uses single cell RNAseq data, where the",
  "cells have labels following a known sequence (e.g. a time series), to identify a small number of genes which place cells in that known order.",
  tags$hr(),
  tags$h3("Choose psupertime parameters"),
  tags$br(),
  "Let's first choose the correct parameters to run psupertime.",
  textInput("y", "Metadata column containing the vector of labels", value = "Sample"),
  textInput("x", "Assay containing gene expression data:", value="RNA"),
  textInput(inputId = "uploadRdsprojectID", label = "Project name : ", value = "psupertime"),
  tags$hr(),
  tags$h3("Seurat object (.rds)"),
  br(),
  "Upload your Seurat object",
  fileInput(inputId = "uploadRdsFile", label = "Choose a Seurat object saved in .rds format", accept = ".rds"),
  actionButton(inputId = "uploadSeuratRdsConfirm", label = "Run psupertime analysis", class="btn-run", icon = icon("check-circle")),
  plotOutput("g1"),
  plotOutput("g5"),
  h4("Gene pseudotime"),
  "Plots profiles of hand-selected genes against psupertime.", br(), br(),
  fluidRow(
    column(
      3, textInput("ps1", "Gene name:", value  = "Nphs1")
    ),
    column
    (6, plotOutput("g2")
    )
  ),
  br(),
  h4("Gene cluster profiles"),
  "Plot the mean profiles of the clustered genes against the pseudotime values",br(),
  fluidRow(
    column(
      3, selectInput("ps2", "Number of clusters:", choices = c(3,4,5), selected=3) 
    ),
    column
    (6, plotOutput("g3"),
    )
  ),
  br(),
  h4("Functional enrichment"), br(),
  fluidRow(
    column(
      3, selectInput("ps3", "Number of clusters:", choices = c(3,4,5), selected=3),
      numericInput("ps4","Minimum number of annotated genes", 5),
      numericInput("ps5", "Maximum p-value", 0.1)
    ),
    column
    (6, plotOutput("g4")
    )
  ),
  tags$h3("Export psupertime object"),
  tags$hr(),
  downloadButton(outputId = "utilitiesConfirmExport", label = "Export app")
)
server <- function(input, output, session) {
  
  options(shiny.maxRequestSize=30*1024^4) 
  #observe_helpers() 
  optCrt="{ option_create: function(data,escape) {return('<div class=\"create\"><strong>' + '</strong></div>');} }" 
  # updateSelectizeInput(session, "ps1", choices = colnames(psuper_obj$x_data), server = TRUE,
  #                      selected = "Nphs1", options = list(
  #                        maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt)))
  
  output$files <- renderTable(input$upload)
  observeEvent(input$uploadSeuratRdsConfirm, {
    tryCatch({
      showModal(modalDialog(div('Psupertime analysis in progress. Please wait...')))
      seurat_object <<- readRDS(input$uploadRdsFile$datapath)
      label<-input$y
      y=factor(seurat_object@meta.data[,input$y])
      assay<-input$x
      x=seurat_object@assays[[input$x]]@data
      var=seurat_object@assays[[input$x]]@var.features
      psuper_obj<-psupertime(x,y, sel_genes='list', gene_list = var)
      go_list     = list()
       for(i in 3:5){go_list[[i-2]]<-psupertime_go_analysis(psuper_obj, org_mapping='org.Mm.eg.db', k=i)}
      # plot gene cluster profiles
      output$g1<- renderPlot({plot_labels_over_psupertime(psuper_obj, label_name='Sample')+ theme(text = element_text(size = 20)) })
      output$g5<- renderPlot({plot_train_results(psuper_obj)+ theme(text = element_text(size = 20)) })
      output$g2<- renderPlot({plot_specified_genes_over_psupertime(psuper_obj,extra_genes = input$ps1, label_name='Sample')+ theme(text = element_text(size = 20)) })
      output$g3<- renderPlot({plot_profiles_of_gene_clusters(go_list[[as.numeric(input$ps2)-2]], label_name='Sample')+ theme(text = element_text(size = 15)) }, width = "auto", height="auto")
      output$g4<- renderPlot({plot_go_results(go_list[[as.numeric(input$ps3)-2]],sig_cutoff = input$ps4, p_cutoff=input$ps5) + theme(text = element_text(size = 12))}, width = "auto", height=1200)

    }, error = function(e) {
      print(paste("Error :  ", e))
      session$sendCustomMessage("handler_alert", "Data input error. Please make sure you are using the correct format.")
    }, finally = { # with or without error
      removeModal()
      Sys.sleep(1) # giving some time for rendering for smoother transition
    })
  })
  
  output$utilitiesConfirmExport <- downloadHandler(
    filename = function() { 
      paste(input$uploadRdsprojectID, Sys.Date(), ".zip", sep="")
    },
    content = function(file) {
      save(psuper_obj,file=file)
    })
}
shinyApp(ui = ui, server = server)

