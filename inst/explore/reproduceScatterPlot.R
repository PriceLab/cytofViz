library(ggplot2)
library(reshape2)
load("../extdata/klf1.plotting.data.RData")
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
    gp <- ggplot(data, aes_string(x = xlab, y = ylab, colour = ev)) +
        facet_wrap(~markers, nrow = grid_row_num, scales = "fixed") +
        scale_colour_gradientn(limits = limits, name = ev, colours = myPalette(zlength * 2)) +
        geom_point(size = pointSize, alpha = alpha) + theme_bw() + coord_fixed() +
        theme(legend.position = "right") + xlab(xlab) + ylab(ylab) + ggtitle(title) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        theme(axis.text=element_text(size=8), axis.title=element_text(size=12,face="bold"))

    return(gp)
}

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


klf1.plotting.data$clusterColor   <- "#000000";
x <- klf1.plotting.data

gp <- with(klf1.plotting.data, scatterPlot(obj = x$obj,
                                           plotMethod = x$plotMethod,
                                           plotFunction = x$plotFunction,
                                           pointSize = x$pointSize,
                                           alpha = x$alpha,
                                           addLabel = x$addLabel,
                                           labelSize = x$labelSize,
                                           sampleLabel = x$sampleLabel,
                                           FlowSOM_k = x$FlowSOM_k,
                                           selectSamples = x$selectSamples,
                                           selectCluster=NULL,
                                           facetPlot = x$facetPlot,
                                           colorPalette = x$colorPalette,
                                           labelRepel = x$labelRepel,
                                           removeOutlier = x$removeOutlier,
                                           clusterColor = x$clusterColor,
                                           globalScale = x$globalScale,
                                           centerScale = x$centerScale
                                           ))
plot(gp)
