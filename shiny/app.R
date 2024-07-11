reqPkg = c("data.table", "Matrix", "hdf5r", "reticulate", "ggplot2", 
           "gridExtra", "glue", "readr", "RColorBrewer", "R.utils", "Seurat", "librarian")
newPkg = reqPkg[!(reqPkg %in% installed.packages()[,"Package"])]
if(length(newPkg)){install.packages(newPkg)}
library(librarian)
shelf(shiny, shinyhelper, data.table, Matrix, DT, hdf5r, reticulate, ggplot2, gridExtra, magrittr, ggdendro, htmltools, SGDDNB/ShinyCell, zip)

ui <- fluidPage(
  tags$h3("Shiny app creator"),
  br(),
  "You can create a shiny app to explore any Seurat scRNA-seq object. The input is a Seurat object in rds format and the output",
  "are the files and code required for the data exploration shiny app",
  tags$hr(),
  textInput(inputId = "uploadRdsprojectID", label = "Project name : ", value = "KFO_329"),
  fileInput(inputId = "uploadRdsFile", label = "Choose a Seurat object saved in .rds format", accept = ".rds"),
  tags$hr(),
  tags$br(),
  textInput("assay", "Assay containing gene expression data:", value="RNA"),
  textInput("slot", "Slot in single-cell assay to plot", value="data"),
  actionButton(inputId = "uploadSeuratRdsConfirm", label = "Create shiny app", class="btn-run", icon = icon("check-circle")),
  tags$h3("Export shiny app"),
  tags$hr(),
  downloadButton(outputId = "utilitiesConfirmExport", label = "Export app"), br(),
  "To run the app locally, use RStudio to open either server.R or ui.R in the shiny app folder and click on Run App in the top right corner. "
)
server <- function(input, output, session) {
  
  options(shiny.maxRequestSize=30*1024^4) 
  
  output$files <- renderTable(input$upload)
  observeEvent(input$uploadSeuratRdsConfirm, {
    tryCatch({
      showModal(modalDialog(div('Creation of app in progress. Please wait...')))
      seurat_object <<- readRDS(input$uploadRdsFile$datapath)
      init_seurat_object <<- seurat_object 
      if(is.null(seurat_object@meta.data)) {
        print("Metadata table is missing")
      } else {
        scConf = createConfig(seurat_object)
        makeShinyApp(seurat_object, scConf, gex.assay = input$assay, gex.slot=input$slot,
                     shiny.title = input$uploadRdsprojectID, shiny.footnotes="") 
        tx  <- readLines("shinyApp/ui.R")
        tx2  <- tx[1:(length(tx)-9)]
        tx2[length(tx2)]<-"br()"
        tx2[length(tx2)+1] <- "))) "
        writeLines(tx2, con="shinyApp/ui.R")
      }
      
    }, error = function(e) {
      print(paste("Error :  ", e))
      session$sendCustomMessage("handler_alert", "Data Input error. Please, refer to the help pages for input format.")
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
      zip::zip(file, 
               #files = list.files("shinyApp/"),
               #files=c("shinyApp/sc1def.rds","shinyApp/server.R","shinyApp/ui.R","shinyApp/sc1conf.rds","shinyApp/sc1gene.rds","shinyApp/sc1gexpr.rds","shinyApp/sc1meta.rds"),
               files=dir('shinyApp', full.names = TRUE))
    })
}
shinyApp(ui = ui, server = server)

         