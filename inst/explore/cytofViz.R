library(shiny)
library(shinydashboard)
if(!exists("tbl"))
   load("tsnePlot.RData")
#------------------------------------------------------------------------------------------------------------------------
.createSidebar <- function()
{
  dashboardSidebar(
     selectInput("chooseProteinFromList", "Choose Protein From List:", LETTERS)
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
     plot(tbl$tsne1 , tbl$tsne2, col=tbl$color)
    })
}
#------------------------------------------------------------------------------------------------------------------------
app <- shinyApp(ui, server)

