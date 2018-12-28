# library(RColorBrewer)
library(shiny)
library(shinydashboard)
library(ggplot2)
#------------------------------------------------------------------------------------------------------------------------
if(!exists("tbl"))
   load("tbl.allDays.allProteins.transformed.RData")
   tbl <- tbl.all

conditions <- c("Phenograph", sort(colnames(tbl)[9:35]))

mtx <- as.matrix(tbl.all[, 9:35])
#------------------------------------------------------------------------------------------------------------------------
assignColor <- function(values, min, max, colors)
{
   numberOfBins <- length(colors) - 1
   bin.size <- (max-min)/numberOfBins
   bin.assignments <- 1 + as.numeric(lapply(values, function(value) floor(value/bin.size)))
   colors[bin.assignments]

} # assignColor
#------------------------------------------------------------------------------------------------------------------------
.createSidebar <- function()
{
  dashboardSidebar(
     tags$head(
         tags$style("#plotRenderingPanel{height:82vh !important;}
                     #imageRenderingPanel{height:82vh !important;}")
         ),
     selectInput("chooseProteinFromList", "Choose Protein From List:", conditions,
                 selectize=FALSE, size=length(conditions))
     ) # dashboardSidebar

} # .createSidebar
#------------------------------------------------------------------------------------------------------------------------
ui <- dashboardPage(

   dashboardHeader(title = "cytof erythropoiesis"),
   .createSidebar(),
   dashboardBody(
     fluidRow(
        box(
          conditionalPanel('input.chooseProteinFromList=="Phenograph"', imageOutput("imageRenderingPanel", height=600)),
          conditionalPanel('input.chooseProteinFromList!="Phenograph"', plotOutput("plotRenderingPanel", height=800)),
          width=12)
          #plotOutput("plot1", height = 800), width=12)
      ) # fluidRow
    ) # dashboarBody

) # ui
#------------------------------------------------------------------------------------------------------------------------
server <- function(input, output) {

  output$imageRenderingPanel <- renderImage({
     imagePath <- normalizePath(file.path("./images", "phenograph.png"))
     printf("imagePath: %s (%s)", imagePath, file.exists(imagePath))
     list(src = imagePath, contentType = 'image/png', width = 800, height = 800, alt = "file not found")
     }, deleteFile = FALSE)

  output$plotRenderingPanel <- renderPlot({
     condition <- input$chooseProteinFromList;
     if(condition == "Phenograph") return()
     transform <- input$transformData
     scale <- input$scaleData
     myPalette <- colorRampPalette(c("#5E4FA2", "#3288BD", "#66C2A5", "#ABDDA4",
                                     "#E6F598", "#FFFFBF", "#FEE08B", "#FDAE61",
                                     "#F46D43", "#D53E4F", "#9E0142"))
     condition.wrapped <- sprintf("`%s`", condition)
     ggplot(tbl, aes_string(x="tsne1", y="tsne2", color=condition.wrapped)) +
        geom_point(size=2) + #, shape=23) +
        scale_colour_gradientn(limits = NULL, name = "", colours = myPalette(nrow(tbl) * 2)) +
        ggtitle(condition) +
        theme(plot.title = element_text(size=20, face="bold", hjust = 0.5))
     # } # ! Phenograph
     })

} # server
#------------------------------------------------------------------------------------------------------------------------
app <- shinyApp(ui, server)

