library(RColorBrewer)
library(shiny)
library(shinydashboard)
library(ggplot2)
#------------------------------------------------------------------------------------------------------------------------
if(!exists("tbl"))
   load("tbl.allDays.allProteins.transformed.RData")
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
     #radioButtons("transformData", "Transform", choices=c("raw", "asinh"), selected="asinh", inline=TRUE),
     #radioButtons("scaleData", "Scale", choices=c("This TF", "All TFs"), selected="This TF", inline=TRUE),
     #selectInput("choosePalette", "Select Palette",
     #            c("spectral.1-10", "spectral.1-20", "spectral.1-5", "spectral.2-10", rownames(tbl.palettes))),
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
     if(condition == "Phenograph"){
        colors <- tbl$color
        par (mai=c(1,1,1,1))
        plot(tbl$tsne1 , tbl$tsne2, col=colors, xlab="tsne1", ylab="tsne2", main=condition, pch=19,
             cex.lab=2, cex.axis=1, cex.main=2, cex.sub=2)
     } else {
        myPalette <- colorRampPalette(c("#5E4FA2", "#3288BD", "#66C2A5", "#ABDDA4",
                                        "#E6F598", "#FFFFBF", "#FEE08B", "#FDAE61",
                                        "#F46D43", "#D53E4F", "#9E0142"))
        condition.wrapped <- sprintf("`%s`", condition)
        ggplot(tbl, aes_string(x="tsne1", y="tsne2", color=condition.wrapped)) +
               geom_point(size=2) + #, shape=23) +
           scale_colour_gradientn(limits = NULL, name = "", colours = myPalette(nrow(tbl) * 2)) +
           ggtitle(condition) +
           theme(plot.title = element_text(size=20, face="bold", hjust = 0.5))
        } # ! Phenograph
     })

} # server
#------------------------------------------------------------------------------------------------------------------------
app <- shinyApp(ui, server)

