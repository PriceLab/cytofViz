library(RColorBrewer)
library(shiny)
library(shinydashboard)
#------------------------------------------------------------------------------------------------------------------------
if(!exists("tbl"))
   load("tbl.allDays.allProteins.RData")
   tbl <- tbl.all

conditions <- c("Phenograph", sort(colnames(tbl)[9:35]))
tbl.palettes <- brewer.pal.info

red.colors <- c("#F0F0F0", brewer.pal(9, "Reds")[1:9])
spectral.1.20.colors <- c("#5E4FA2", "#466DB0", "#348BBB", "#50A9AF", "#6DC4A4", "#91D3A4", "#B4E0A2", "#D3ED9B", "#EBF7A0",
                          "#F8FCB4", "#FEF6B1", "#FEE695", "#FDD07D", "#FDB567", "#F99655", "#F47346", "#E65948", "#D6404E",
                          "#BA2148", "#9E0142")
spectral.1.10.colors <- c("#5E4FA2", "#378EBA", "#75C8A4", "#BEE4A0", "#F1F9A9", "#FEEDA2", "#FDBE6F", "#F67B49", "#D8434D", "#9E0142")
spectral.1.5.colors  <- c("#5E4FA2", "#88CFA4", "#FFFFBF", "#F88D52", "#9E0142")
spectral.2.10.colors <- c("#00007F", "#0000F0", "#0062FF", "#00D4FF", "#8DFFFF", "#FFFF8D", "#FFD400", "#FF6200", "#F00000", "#7F0000")

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
     radioButtons("transformData", "Transform", choices=c("raw", "asinh"), selected="asinh", inline=TRUE),
     radioButtons("scaleData", "Scale", choices=c("This TF", "All TFs"), selected="This TF", inline=TRUE),
     selectInput("choosePalette", "Select Palette",
                 c("spectral.1-10", "spectral.1-20", "spectral.1-5", "spectral.2-10", rownames(tbl.palettes))),
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
     transform <- input$transformData
     scale <- input$scaleData
       # TODO: refactor this into a tested function
     if(condition == "Phenograph")
        colors = tbl$color
     else{
        palette.name <- input$choosePalette
        if(grepl("spectral.", palette.name, fixed=TRUE)){
          if(palette.name == "spectral.1-20"){
             colors <- spectral.1.20.colors
             color.count <- length(spectral.1.20.colors)
             }
          if(palette.name == "spectral.1-10"){
             colors <- spectral.1.10.colors
             color.count <- length(spectral.1.10.colors)
             }
          if(palette.name == "spectral.1-5"){
             colors <- spectral.1.5.colors
             color.count <- length(spectral.1.5.colors)
             }
          if(palette.name == "spectral.2-10"){
             colors <- spectral.2.10.colors
             color.count <- length(spectral.2.10.colors)
             }

        } else {
           color.count <- tbl.palettes[palette.name, "maxcolors"]
           colors <- brewer.pal(color.count, palette.name)
           }
        working.mtx <- mtx
        if(transform == "asinh")
           working.mtx <- asinh(mtx)
        vec <- working.mtx[, condition]
        if(scale == "This TF"){
           minValue <- min(vec)
           maxValue <- max(vec)
        } else {
           minValue <- min(working.mtx)
           maxValue <- max(working.mtx)
           }
        colors <- assignColor(vec, minValue, maxValue, colors)
        }  # else: not Phenograph
     #browser()
     plot(tbl$tsne1 , tbl$tsne2, col=colors, main=condition, pch=19)
     })

} # server
#------------------------------------------------------------------------------------------------------------------------
app <- shinyApp(ui, server)

