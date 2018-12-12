#' import shiny
#' import shinydashboard
#' import shinyjs
#' import later
#' import igvShiny
#' import DT
#' import VariantAnnotation
#' import cytofSGM
#' importFrom  RColorBrewer brewer.pal
#' import later
#------------------------------------------------------------------------------------------------------------------------
#' @name CytofViz
#' @rdname CytofViz
#' @aliases CytofViz
#------------------------------------------------------------------------------------------------------------------------
state <- new.env(parent=emptyenv())
state$models <- list()
state$tbl.chipSeq <- NULL
model.count <- 0   # for creating default model names
MAX.TF.COUNT <- 50   # we display the top max.tf.cout TFs coming back from cytof

colors <- brewer.pal(8, "Dark2")
totalColorCount <- length(colors)
state$colorNumber <- 0
#------------------------------------------------------------------------------------------------------------------------
.CytofViz <- setClass ("CytofViz",
                    representation = representation(
                       projectName="character",
                       project="CytofProject",
                       quiet="logical")
                    )
#------------------------------------------------------------------------------------------------------------------------
setGeneric('createUI',      signature='obj', function(obj) standardGeneric("createUI"))
setGeneric('createServer',  signature='obj', function(obj, session, input, output) standardGeneric("createServer"))
setGeneric('createApp',     signature='obj', function(obj, port=NA_integer_) standardGeneric("createApp"))
#------------------------------------------------------------------------------------------------------------------------
#' Create an CytofViz object
#'
#' @description
#' a shiny app
#'
#' @rdname CytofViz
#'
#' @param projectName A string indicating a CytofProject subclass
#' @param quiet A logical indicating whether or not the Cytof object should print output
#'
#' @return An object of the CytofViz class
#'
#' @export
#'
CytofViz <- function(projectName, quiet=TRUE)
{
   require(projectName, character.only=TRUE)
   initialization.command <- sprintf("cytofProject <- %s()", projectName)
   eval(parse(text=initialization.command))
   setTargetGene(cytofProject, getSupportedGenes(cytofProject)[1])

   state$tbl.enhancers <- getEnhancers(cytofProject)
   state$tbl.dhs <- getEncodeDHS(cytofProject)
   state$tbl.transcripts <- getTranscriptsTable(cytofProject)

   obj <- .CytofViz(projectName=projectName, project=cytofProject, quiet=quiet)

   obj

} # CytofViz ctor
#------------------------------------------------------------------------------------------------------------------------
#' Create the user interface
#'
#' @rdname createUI
#' @aliases createUI
#'
#' @param obj An object of class CytofViz
#'
#' @export
#'
setMethod('createUI', 'CytofViz',

   function(obj){

      ui <- dashboardPage(
        dashboardHeader(title=sprintf("cytof %s", getTargetGene(obj@project))),
        .createSidebar(obj@project),
        dashboardBody(
           tags$head(tags$style(HTML('
        .main-header .logo {
          font-family: "Georgia", Times, "Times New Roman", serif;
          font-weight: bold;
          font-size: 24px;
          }
        #igvShiny{
          background: white;
          border: 1px solid black;
          border-radius: 5px;
          margin: 5px;
          margin-right: 15px;
          overflow: hidden;
          }
        #table{
          border: 1px solid black;
          border-radius: 5px;
          }
       '))),
       .createBody(obj@project)),
       useShinyjs()
      ) # dashboardPage
     return(ui)
   })
#------------------------------------------------------------------------------------------------------------------------
#' Create the shiny server
#'
#' @rdname createServer
#' @aliases createServer
#'
#' @param obj An object of class CytofViz
#' @param session A shiny session
#' @param input A shiny input object
#' @param output A shiny output object
#'
#' @export
#'
setMethod('createServer', 'CytofViz',

   function(obj, session, input, output){

   observeEvent(input$setOffListGeneButton, ignoreInit=TRUE, {
      newGene <- toupper(isolate(input$offListGene))
      legit <- recognizedGene(obj@project, newGene)
      if(!legit){
        msg <- sprintf("'%s' not found in current datasets.", newGene)
        showModal(modalDialog(title="gene nomination error", msg))
        }
      if(legit){
        shinyjs::html(selector=".logo", html=sprintf("cytof %s", newGene), add=FALSE)
        setTargetGene(obj@project, newGene)
        state$tbl.enhancers <- getEnhancers(obj@project)
        state$tbl.dhs <- getEncodeDHS(obj@project)
        state$tbl.transcripts <- getTranscriptsTable(obj@project)
        showGenomicRegion(session, getGeneRegion(obj@project, flankingPercent=20))
        }
      })

   observeEvent(input$chooseGeneFromList, ignoreInit=TRUE, {
       #  TODO: consolidate this code with that above, setOffListGeneButton code
       newGene <- isolate(input$chooseGeneFromList)
       shinyjs::html(selector=".logo", html=sprintf("cytof %s", newGene), add=FALSE)
       setTargetGene(obj@project, newGene)
       showGenomicRegion(session, getGeneRegion(obj@project, flankingPercent=20))
       printf("new targetGene: %s", getTargetGene(obj@project))
       state$tbl.enhancers <- getEnhancers(obj@project)
       state$tbl.dhs <- getEncodeDHS(obj@project)
       state$tbl.transcripts <- getTranscriptsTable(obj@project)
       loadAndDisplayRelevantVariants(obj@project, session, newGene)
       })


   output$igvShiny <- renderIgvShiny({
      options <- list(genomeName="hg38",
                      initialLocus=getGeneRegion(obj@project, flankingPercent=20),
                      displayMode="EXPANDED",
                      trackHeight=300)
      igvShiny(options) # , height=800)
      })

   output$table = DT::renderDataTable(data.frame(),
                                      width="800px",
                                      class='nowrap display',
                                      selection="single",
                                      extensions="FixedColumns",
                                      options=list(scrollX=TRUE,
                                                   scrollY="500px",
                                                   dom='t',
                                                   paging=FALSE,
                                                   autowWdth=FALSE,
                                                   fixedColumns=list(leftColumns=1)
                                                   ))

   observeEvent(input$table_rows_selected, {
      selectedTableRow <- isolate(input$table_rows_selected)
      dispatch.rowClickInModelTable(obj@project, session, input, output, selectedTableRow)
      }) # observe row selection event

   observeEvent(input$currentGenomicRegion, {
      new.region <- isolate(input$currentGenomicRegion)
      #printf("new region: %s", new.region)
      state[["chromLocRegion"]] <- new.region
      })

   setupIgvAndTableToggling(session, input);
   setupAddTrack(obj@project, session, input, output)
   setupDisplayRegion(obj@project, session, input, output)
   setupBuildModel(obj@project, session, input, output)
   })

#------------------------------------------------------------------------------------------------------------------------
.createBody <- function(project)
{
   body <- #fluidRow(
      tabItems(
         tabItem(tabName="igvAndTable",
            fluidRow(
               column(width=9, id="igvColumn", igvShinyOutput('igvShiny', height="2000px")),
               column(width=3, id="dataTableColumn",
                      selectInput("modelSelector", NULL,  "- no models yet -"),
                      DTOutput("table"),
                      br(),
                      #wellPanel(
                        h4("TF selection will display:"),
                        selectInput("selectRowAction", NULL,  c("(no action)", "Footprints", "ChIP-seq hits",
                                                                sprintf("Plot expression, TF vs target gene")))
                      ) # column
               ) # fluidRow
            ), # tabItem 1
         .createBuildModelTab(project),
         .createVideoTab()
         ) # tabItems

   return(body)

} # createBody
#------------------------------------------------------------------------------------------------------------------------
.createSidebar <- function(cytofProject)
{
  dashboardSidebar(
    sidebarMenu(id="sidebarMenu",
       menuItem("IGV and Current Models", tabName = "igvAndTable"),
       menuItem("Build a new Model",      tabName = "buildModels"),
       menuItem("Introductory video",     tabName = "video")
      ),
    conditionalPanel(id="igvTableWidgets",
        condition = "input.sidebarMenu == 'igvAndTable'",
        actionButton(inputId = "igvHideButton", label = "Toggle IGV"),
        actionButton(inputId = "tableHideButton", label = "Toggle table"),
        selectInput("chooseGeneFromList", "Choose Gene From List:", c("", getSupportedGenes(cytofProject))),
        # div(style="display: inline-block;vertical-align:top; width:100px; margin:0px; !important;",
        textInput("offListGene", label=NULL, placeholder="enter off-list gene here"),
        #div(style="display: inline-block;vertical-align:top; width: 40px; margin-left:0px; margin-top:8px; !important;",
        actionButton(inputId="setOffListGeneButton", label="Set off-list gene"),
        selectInput("addTrack", "Add Track:", .additionalTracksOnOffer(cytofProject)),
        selectInput("displayGenomicRegion", "Display Genomic Region:",
                    c("",
                      "Full Gene" = "fullGeneRegion",
                      "Full Enhancers Region" = "fullEnhancerRegion")),
        actionButton(inputId="removeUserAddedTracks", label="Remove Added Tracks")
        ) # conditionalPanel
     ) # dashboardSidebar

} # .createSidebar
#------------------------------------------------------------------------------------------------------------------------
.additionalTracksOnOffer <- function(cytofProject)
{
   trackOfferings <- list("",
                          "Enhancers"="enhancers",
                          "DNase HS (encode, clustered)"="dhs")

   variantTrackNames <- getVariantDatasetNames(cytofProject)

   for(i in seq_len(length(variantTrackNames))){
      trackName <- variantTrackNames[i]
      if(grepl("gwas", trackName, ignore.case=TRUE)){
         new.list <- c(trackName)
         names(new.list) <- trackName
         trackOfferings <- c(trackOfferings, new.list)
         } # if 'gwas' in the name
      if(grepl(".variants", trackName, ignore.case=TRUE)){
         new.list <- c(trackName)
         names(new.list) <- trackName
         trackOfferings <- c(trackOfferings, new.list)
         } # if '.variants' in the name
      } # for i

   return(trackOfferings)

} # additionalTracksOnOffer
#------------------------------------------------------------------------------------------------------------------------
.createBuildModelTab <- function(project)
{
   footprintDatabases <- getFootprintDatabaseNames(project)
   expressionMatrixNames <- getExpressionMatrixNames(project)

   tab <- tabItem(tabName="buildModels",
             h4("Build cytof regulatory model from DNase footprints in the genomic region currently displayed in IGV, with these constraints:"),
             br(),
             fluidRow(
                column(3, offset=1, checkboxGroupInput("footprintDatabases", "Footprint Databases",
                                                       footprintDatabases,
                                                       selected=footprintDatabases[1])),
                column(4, radioButtons("intersectWithRegions", "Intersect footprints with:",
                                             c("GeneHancer" = "genehancer",
                                               "Encode DHS" = "encodeDHS",
                                               "GeneHancer OR Encode DHS" = "geneHancerOrEncode",
                                               "GeneHancer AND Encode DHS" = "geneHancerAndEncode",
                                               "Use all footprints in region" = "allDNAForFootprints"),
                                       selected="allDNAForFootprints")),
                column(3, radioButtons("motifMapping", "Map motifs to TFs via:",
                                       c("MotifDb", "TFClass", "MotifDb + TFClass"),
                                       selected="MotifDb"))
                ),
             fluidRow(
                column(6, selectInput("expressionSet", "With this gene expression dataset:",
                                      c("", expressionMatrixNames))),
                column(5, textInput("modelNameTextInput", "Model name", width=500))),
             fluidRow(column(width=1, offset=4, actionButton("buildModelButton", "Build model"))),
             br(),br(),
             wellPanel(style="overflow-y:scroll; height:200px", pre(id = "console"))
             ) # tabItem

   return(tab)

} # .createBuildModelTab
#------------------------------------------------------------------------------------------------------------------------
.createVideoTab <- function()
{
   file <- system.file(package="CytofViz", "extdata", "video", "videoTutorial.html")

   tab <- tabItem(tabName="video",
                  h4("Using CytofViz: a video tutorial"),
                  includeHTML(file)
                  )
   return(tab)

} # .createVideoTab
#------------------------------------------------------------------------------------------------------------------------
#' Create a runnable shiny app
#'
#' @rdname createApp
#' @aliases createApp
#'
#' @param obj An object of class CytofViz
#' @param port An integer (e.g., 3838, 60041) NA_integer_ by default
#'
#' @export
#'
setMethod('createApp', 'CytofViz',

  function(obj, port=NA_integer_){

     server <- function(session, input, output){
        x <- createServer(obj, session, input, output)
        }


     if(Sys.info()[["nodename"]] == "riptide.local"){
        shinyOptions <- list(host="0.0.0.0", launch.browser=TRUE)
     }else{
        shinyOptions=list(launch.browser=FALSE, host='0.0.0.0', port=port)
        }

     app <- shinyApp(createUI(obj), server, options=shinyOptions)

     return(app)

  })
#------------------------------------------------------------------------------------------------------------------------
