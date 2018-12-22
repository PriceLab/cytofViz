library(ggplot2)
library(reshape2)
print(load("../extdata/klf1.plotting.data.RData"))
#------------------------------------------------------------------------------------------------------------------------
#melt <- function (data, ..., na.rm = FALSE, value.name = "value")
#{
#    UseMethod("melt", data)
#}
#------------------------------------------------------------------------------------------------------------------------
cytof_wrap_colorPlot <- function(data, xlab, ylab, markers, scaleMarker = FALSE,
                                 colorPalette, # = c("bluered", "spectral1", "spectral2", "heat"),
                                 limits = NA,
                                 pointSize=1,
                                 alpha = 1,
                                 removeOutlier = TRUE){

   browser()
   #removeOutlier <- FALSE
   xyz <- "in cytof_wrap_colorPlot"
    remove_outliers <- function(x, na.rm = TRUE, ...) {
        qnt <- quantile(x, probs=c(.02, .98), na.rm = na.rm, ...)
        x[x <= qnt[1]] <- qnt[1]
        x[x >= qnt[2]] <- qnt[2]
        x
    }

    data <- as.data.frame(data)
    title <- "Marker Expression Level Plot"
    data <- data[,c(xlab, ylab, markers)]

    if(removeOutlier){
        for(m in markers){
            printf("--- removing outliers in column %s", m)
            data[[m]] <- remove_outliers(data[ ,m])
        }
    }

    if(scaleMarker){
        data[ ,markers] <- scale(data[ ,markers], center = TRUE, scale = TRUE)
        ev <- "ScaledExpression"
        data <- melt(data, id.vars = c(xlab, ylab),
                     measure.vars = markers,
                     variable.name = "markers",
                     value.name = ev)
    }else{
        ev <- "Expression"
        data <- melt(data, id.vars = c(xlab, ylab),
                     measure.vars = markers,
                     variable.name = "markers",
                     value.name = ev)
    }


    colorPalette <- "spectral1" # match.arg(colorPalette)
    switch(colorPalette,
           bluered = {
               myPalette <- colorRampPalette(c("blue", "white", "red"))
           },
           spectral1 = {
               myPalette <- colorRampPalette(c("#5E4FA2", "#3288BD", "#66C2A5", "#ABDDA4",
                                               "#E6F598", "#FFFFBF", "#FEE08B", "#FDAE61",
                                               "#F46D43", "#D53E4F", "#9E0142"))
           },
           spectral2 = {
               myPalette <- colorRampPalette(rev(c("#7F0000","red","#FF7F00","yellow","white",
                                                   "cyan", "#007FFF", "blue","#00007F")))
           },
           heat = {
               myPalette <- colorRampPalette(heat.colors(50))
           }
    )
    zlength <- nrow(data)
    grid_row_num <- round(sqrt(length(markers)))
    browser()
    gp <- ggplot(data, aes_string(x = xlab, y = ylab, colour = ev)) +
        facet_wrap(~markers, nrow = grid_row_num, scales = "fixed") +
        scale_colour_gradientn(limits = limits, name = ev, colours = myPalette(zlength * 2)) +
        geom_point(size = pointSize, alpha = alpha) + theme_bw() + coord_fixed() +
        theme(legend.position = "right") + xlab(xlab) + ylab(ylab) + ggtitle(title) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        theme(axis.text=element_text(size=8), axis.title=element_text(size=12,face="bold"))

    return(gp)

} # cytof_wrap_colorPlot
#------------------------------------------------------------------------------------------------------------------------
myPlot <- function(tbl)
{
   myPalette <- colorRampPalette(c("#5E4FA2", "#3288BD", "#66C2A5", "#ABDDA4",
                                   "#E6F598", "#FFFFBF", "#FEE08B", "#FDAE61",
                                   "#F46D43", "#D53E4F", "#9E0142"))
   zlength <- nrow(tbl)
   grid_row_num <- 1
   xlab <- colnames(tbl)[1]
   ylab <- colnames(tbl)[2]
   ev <- "ScaledExpression"
   markers <- colnames(tbl)[3]

   gp <- ggplot(data,
                aes_string(x = xlab, y = ylab, colour = ev)) +
      facet_wrap(~markers, nrow = grid_row_num, scales = "fixed") +
      scale_colour_gradientn(limits = limits, name = ev, colours = myPalette(zlength * 2)) +
      geom_point(size = pointSize, alpha = alpha) + theme_bw() + coord_fixed() +
      theme(legend.position = "right") + xlab(xlab) + ylab(ylab) + ggtitle(title) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      theme(axis.text=element_text(size=8), axis.title=element_text(size=12,face="bold"))



} # myPlot
#------------------------------------------------------------------------------------------------------------------------
scatterPlot <- function(obj,
                        plotMethod,
                        plotFunction,
                        pointSize=1,
                        alpha = 1,
                        addLabel=TRUE,
                        labelSize=1,
                        sampleLabel = TRUE,
                        FlowSOM_k = 40,
                        selectCluster=NULL,
                        selectSamples,
                        facetPlot = FALSE,
                        colorPalette = bluered,
                        labelRepel = FALSE,
                        removeOutlier = TRUE,
                        clusterColor,
                        globalScale = TRUE,
                        centerScale = FALSE)
{

    data <- data.frame(obj$expressionData,
                       obj$dimReducedRes[[plotMethod]],
                       do.call(cbind, obj$clusterRes),
                       check.names = FALSE,
                       stringsAsFactors = FALSE)

    Markers <- obj$allMarkers

    xlab <- colnames(obj$dimReducedRes[[plotMethod]])[1]
    ylab <- colnames(obj$dimReducedRes[[plotMethod]])[2]
    row.names(data) <- row.names(obj$expressionData)

    clusterMethods <- names(obj$clusterRes)
    samples <- sub("_[0-9]*$", "", row.names(obj$expressionData))
    data <- data[samples %in% selectSamples, ,drop=FALSE]
    nsamples <- samples[samples %in% selectSamples]
    data$sample <- nsamples
    sample_num <- length(unique(nsamples))

   browser()
   xyz <- "about to branch on plotFunction"

    if(plotFunction == "Density"){
        colPalette <- colorRampPalette(c("blue", "turquoise", "green",
                                         "yellow", "orange", "red"))
        densCol <- densCols(data[, c(xlab, ylab)], colramp = colPalette)
        data$densCol <- densCol
        gp <- ggplot(data, aes_string(x=xlab, y=ylab)) +
            geom_point(colour=densCol, size = pointSize) + ggtitle("Density Plot") +
            theme(legend.position = "right") + xlab(xlab) + ylab(ylab) + theme_bw() +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
            theme(axis.text=element_text(size=14), axis.title=element_text(size=18,face="bold"))
    }else if(plotFunction == "None"){
        gp <- ggplot(data, aes_string(x=xlab, y=ylab)) +
            geom_point(size = pointSize) + ggtitle("Dot Plot") +
            xlab(xlab) + ylab(ylab) + theme_bw() +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
            theme(axis.text=element_text(size=14), axis.title=element_text(size=18,face="bold"))
    }else if(plotFunction == "Sample"){
        size_legend_row <- ceiling(sample_num/4)
        sample <- "sample"
        gp <- ggplot(data, aes_string(x=xlab, y=ylab, colour = sample)) +
            geom_point(size = pointSize) + ggtitle("Color By Sample") +
            xlab(xlab) + ylab(ylab) + theme_bw() + theme(legend.position = "bottom") +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
            theme(axis.text=element_text(size=14), axis.title=element_text(size=18,face="bold")) +
            guides(colour = guide_legend(nrow = size_legend_row, override.aes = list(size = 4)))
    }else if(plotFunction == "All Markers"){
        gp <- cytof_wrap_colorPlot(data = data,
                              xlab = xlab,
                              ylab = ylab,
                              markers = colnames(obj$expressionData),
                              colorPalette = colorPalette,
                              limits = NULL,
                              pointSize = pointSize,
                              removeOutlier = TRUE)

    }else if(plotFunction == "All Markers(scaled)"){
        gp <- cytof_wrap_colorPlot(data = data,
                                   xlab = xlab,
                                   ylab = ylab,
                                   markers = colnames(obj$expressionData),
                                   scaleMarker = TRUE,
                                   colorPalette = colorPalette,
                                   limits = NULL,
                                   pointSize = pointSize,
                                   removeOutlier = TRUE)

    }else if(plotFunction %in% clusterMethods){

        if(!is.null(selectCluster)){
            clusterIDs <- as.character(data[,plotFunction])
            selectCluster <- as.character(selectCluster)
            data <- data[clusterIDs %in% selectCluster, ,drop=FALSE]
        }
        clusterVec <- obj$clusterRes[[plotFunction]]
        ## make sure they are not factors before transforming to factors
        selectColors <- match(levels(as.factor(data[,plotFunction])), levels(as.factor(clusterVec)))
        clusterColor <- clusterColor[selectColors]

        gp <- cytof_clusterPlot(data = data,
                                xlab = xlab,
                                ylab = ylab,
                                cluster = plotFunction,
                                sample = "sample",
                                title = plotFunction,
                                type = ifelse(facetPlot, 2, 1),
                                point_size = pointSize,
                                addLabel = addLabel,
                                labelSize = labelSize,
                                sampleLabel = sampleLabel,
                                labelRepel = labelRepel,
                                fixCoord = FALSE,
                                clusterColor = clusterColor)
    }else{
        limits <- NULL
        if(globalScale){
          exprData <- obj$expressionData
          markers <- colnames(exprData)
          glimits <- quantile(exprData, probs=c(.02, .98), na.rm = TRUE)
          local.bounds <- as.data.frame(lapply(markers, function(x) quantile(exprData[,x], probs=c(.02, .98), na.rm = TRUE)), col.names = markers)
          gmax <- ifelse(max(local.bounds[2,]) < glimits[2], glimits[2], max(local.bounds[2,]))
          gmin <- ifelse(min(local.bounds[1,]) > glimits[1],min(local.bounds[1,]), glimits[1])
          limits <- c(gmin, gmax)
        }
        if(length(plotFunction > 1)){
          gp <- cytof_wrap_colorPlot(data = data,
                                     xlab = xlab,
                                     ylab = ylab,
                                     markers = plotFunction,
                                     colorPalette = colorPalette,
                                     limits = limits,
                                     scaleMarker = centerScale,
                                     pointSize = pointSize,
                                     alpha = alpha,
                                     removeOutlier = TRUE)
        }else{
          gp <- cytof_colorPlot(data = data,
                                xlab = xlab,
                                ylab = ylab,
                                zlab = plotFunction,
                                colorPalette = colorPalette,
                                limits = limits,
                                pointSize = pointSize,
                                alpha = alpha,
                                removeOutlier = TRUE)
        }
    }

    return(gp)
}
#------------------------------------------------------------------------------------------------------------------------
newf <- function()
{
   #data <- data.frame(obj$expressionData,
   #                    obj$dimReducedRes[[plotMethod]],
   #                    do.call(cbind, obj$clusterRes),
   #                    check.names = FALSE,
   #                    stringsAsFactors = FALSE)

   gp <- cytof_wrap_colorPlot(data = tbl.klf1,
                              xlab = "tsne1",
                              ylab = "tsne2",
                              markers = "KLF1",
                              colorPalette = "spectral1",
                              limits = NULL,
                              pointSize = 1,
                              alpha = 1,
                              removeOutlier = TRUE)

   plot(gp)

} # newf
#------------------------------------------------------------------------------------------------------------------------
if(!exists("tbl.all"))
   load("../docker/tbl.allDays.allProteins.RData")

tbl.klf1 <- tbl.all[, c("tsne1", "tsne2", "KLF1")]
tbl.klf1$KLF1 <- asinh(tbl.klf1$KLF1)

klf1.plotting.data$clusterColor   <- "#000000";
xold <- klf1.plotting.data

#------------------------------------------------------------------------------------------------------------------------
old <- function()
{
   gp <- with(klf1.plotting.data, scatterPlot(obj = xold$obj,
                                              plotMethod = xold$plotMethod,
                                              plotFunction = xold$plotFunction,
                                              pointSize = xold$pointSize,
                                              alpha = xold$alpha,
                                              addLabel = xold$addLabel,
                                              labelSize = xold$labelSize,
                                              sampleLabel = xold$sampleLabel,
                                              FlowSOM_k = xold$FlowSOM_k,
                                              selectSamples = xold$selectSamples,
                                              selectCluster=NULL,
                                              facetPlot = xold$facetPlot,
                                              colorPalette = xold$colorPalette,
                                              labelRepel = xold$labelRepel,
                                              removeOutlier = xold$removeOutlier,
                                              clusterColor = xold$clusterColor,
                                              globalScale = xold$globalScale,
                                              centerScale = xold$centerScale
                                              ))
   plot(gp)

} # old
#------------------------------------------------------------------------------------------------------------------------
tryToGetGoodExpressionData <- function()
{
   dir <- "~/s/data/jeffRanish/fromMarjorieBrand/finalAnalysis/cytofkit-without-Bcl11a-June-2018-USED-FOR-PAPER/Results.orig"
   file <- "cytofkit-without-Bcl11a.RData"
   x <- get(load(file.path(dir, file)))
   names(x) # "expressionData"       "dimReductionMethod"   "visualizationMethods" "dimReducedRes"
            # "clusterRes"           "progressionRes"       "projectName"          "rawFCSdir"            "resultDir"
   class(x$expressionData)
   mtx <- x$expressionData
   dim(mtx)   # 48076 x 27
   klf1.colname <- grep("klf1", colnames(mtx), ignore.case=TRUE, value=TRUE)  # [1] "Ho165Di<165Ho_KLF1_1_>"
   klf1.vec <- mtx[, klf1.colname]
   hist(klf1.vec)
   hist(asinh(klf1.vec))

} # trytoGetGoodExpressionData
#------------------------------------------------------------------------------------------------------------------------
removeOutliers <- function(x, na.rm = TRUE, ...)
{
   qnt <- quantile(x, probs=c(.02, .98), na.rm = na.rm, ...)
   x[x <= qnt[1]] <- qnt[1]
   x[x >= qnt[2]] <- qnt[2]
   x

} # removeOutliers
#------------------------------------------------------------------------------------------------------------------------
stealMatrix <- function()
{
   x.stolen <- get(load("~/github/cytofViz/inst/extdata/klf1.plotting.data.RData"))
   mtx.stolen <- x.stolen$obj$expressionData
   preferred.names <- c("CD90", "TAL1", "GATA2", "CD49F", "KAT3B",
                        "PU.1", "CD123", "CD44", "CD36", "GATA1", "RUNX1", "ATRX",
                        "KLF1", "CD71", "CD45RA", "c-EBPa", "CD235ab", "CD34",
                        "MAFG", "NFE2p45", "IKZF1", "FLI1", "CD41", "BACH1", "CD38",
                        "c-JUN", "c-MYC")

   colnames(mtx.stolen) <- preferred.names
   m2 <- apply(mtx.stolen, 2, removeOutliers)
   rownames(m2) <- NULL
   tbl.all <- cbind(tbl.all[, 1:8], as.data.frame(m2))
   quartz(); hist(tbl.all[, "KLF1"], main="stolen 5:41a")

   save(tbl.all, file="../docker/tbl.allDays.allProteins.transformed.RData")

} # stealMatrix
#------------------------------------------------------------------------------------------------------------------------
experiment.with.ggplot <- function()
{
   if(!exists("tbl")){
      tbl.all <- get(load(file="~/github/cytofViz/inst/docker/tbl.allDays.allProteins.transformed.RData"))
      }

   tbl <- tbl.all[, c("tsne1", "tsne2", "KLF1")]
   #ggplot(mtcars, aes(x=wt, y=mpg)) + geom_point(size=2, shape=23)
   #ggplot(tbl, aes(x=tsne1, y=tsne2)) + geom_point(size=2, shape=23)

   myPalette <- colorRampPalette(c("#5E4FA2", "#3288BD", "#66C2A5", "#ABDDA4",
                                   "#E6F598", "#FFFFBF", "#FEE08B", "#FDAE61",
                                   "#F46D43", "#D53E4F", "#9E0142"))
   #ev <- "ScaledExpression"
   quartz()
   poi <- "KLF1"
   ggplot(tbl, aes_string(x="tsne1", y="tsne2", color=poi)) +
          geom_point(size=1) + #, shape=23) +
          scale_colour_gradientn(limits = NULL, name = "KLF1", colours = myPalette(nrow(tbl) * 2))

} # experiment.with.ggplot
#------------------------------------------------------------------------------------------------------------------------
