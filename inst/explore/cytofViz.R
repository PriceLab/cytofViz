library(RColorBrewer)
library(shiny)
library(shinydashboard)
#------------------------------------------------------------------------------------------------------------------------
if(!exists("tbl"))
   load("~/github/cytofViz/inst/explore/tbl.allDays.allProteins.Rdata")
   tbl <- tbl.all

conditions <- c("Phenograph", sort(colnames(tbl)[9:35]))
red.colors <- c("#F0F0F0", brewer.pal(9, "Reds")[1:9])
min <- min(asinh(tbl.all[, 9:35]))
max <- max(asinh(tbl.all[, 9:35]))

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
      box(plotOutput("plot1", height = 800), width=12)
      )
    ) # dashboarBody

) # ui
#------------------------------------------------------------------------------------------------------------------------
server <- function(input, output) {

  output$plot1 <- renderPlot({
     condition <- input$chooseProteinFromList;
     if(condition == "Phenograph")
        colors = tbl$color
     else{
        vec <- asinh(tbl[, condition])
        colors <- assignColor(vec, min, max, red.colors)
        }
     plot(tbl$tsne1 , tbl$tsne2, col=colors, main=condition)
     })

} # server
#------------------------------------------------------------------------------------------------------------------------
app <- shinyApp(ui, server)

