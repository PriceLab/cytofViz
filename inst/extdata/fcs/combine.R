library(flowCore)
print(load("../../explore/tsnePlot.RData"))
tbl.tsne <- tbl
dim(tbl.tsne)  # 48076 8

fcs.files <- grep("fcs$", list.files("."), value=TRUE)
count <- length(fcs.files)

markers <- c("CD34", "CD36", "CD71", "CD38", "CD45RA", "CD123", "CD49f", "CD90", "CD44", "CD41", "CD235ab")
length(markers)
tfs <- c("GATA1", "PU.1", "ATRX", "c-Myc", "KLF1", "FLI1", "TAL1", "GATA2", "RUNX1", "NFE2p45", "BACH1", "IKZF1", "MAFG", "c-JUN", "KAT3B", "C-EBPa")
length(tfs)

#----------------------------------------------------------------------------------------------------
# build a table of identifiers
x.markers <- as.list(sapply(markers, function(marker) as.character(grep(marker, tbl$desc, v=TRUE, ignore.case=TRUE))))
x.tfs <- as.list(sapply(tfs, function(tf) as.character(grep(tf, tbl$desc, v=TRUE, ignore.case=TRUE))))
x.both <- c(x.markers, x.tfs)
tbl.ids <- data.frame(id=as.character(x.both), stringsAsFactors=FALSE)
rownames(tbl.ids) <- names(x.both)

#----------------------------------------------------------------------------------------------------

i <- 1
x <- read.FCS(fcs.files[i], transformation=TRUE, alter.names=TRUE)
tbl.fcs <- pData(parameters(x))
dim(tbl.fcs)
mtx <- exprs(x)  # [1] 3740   57
dim(mtx)

gata1 <- vector("numeric", count)
gata2 <- vector("numeric", count)


for(i in seq_len(count)){
   x <- read.FCS(fcs.files[i], transformation=TRUE, alter.names=TRUE)
   tbl <- pData(parameters(x))
   mtx <- exprs(x)  # [1] 3740   57
   gata1.rowname <-  as.character(tbl$name[grep("GATA1", tbl$desc)]) # [1] "Dy163Di"
   gata2.rowname <-  as.character(tbl$name[grep("GATA2", tbl$desc)]) # [1] "Dy163Di"
   #gata1[i] <- fivenum(asinh(mtx[, gata1.rowname]))[3]
   #gata2[i] <- fivenum(asinh(mtx[, gata2.rowname]))[3]
   #gata1[i] <- mean(asinh(mtx[, gata1.rowname]))
   #gata2[i] <- mean(asinh(mtx[, gata2.rowname]))
   gata1[i] <- sd(asinh(mtx[, gata1.rowname]))
   gata2[i] <- sd(asinh(mtx[, gata2.rowname]))
   #gata1[i] <- sd((mtx[, gata1.rowname]))
   #gata2[i] <- sd((mtx[, gata2.rowname]))
   #printf("%8.2f   %8.2f",, mean(mtx[, gata2.rowname]))
   }

plot(gata1, type="b", ylim=c(0,max(c(gata1, gata2)) * 1.1), col="blue", main="sd asinh count")
lines(gata2, type="b", col="red")
legend(10, 0.99 * max(c(gata1, gata2)), c("GATA1", "GATA2"), c("blue", "red"))
