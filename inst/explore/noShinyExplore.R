# three kinds of files, data structures:
#
#   1) cytofkit-without-Bcl11a_tsne_dimension_reduced_data.csv: has
#      read in as tbl:
#          tsne1     tsne2                         sampleName
#       8.932996 11.170587 export_c01_sample_01_0_Day0_live_1
#       6.301791 12.065454 export_c01_sample_01_0_Day0_live_2
#       5.930183 15.145040 export_c01_sample_01_0_Day0_live_3
#       9.334115  8.528194 export_c01_sample_01_0_Day0_live_4
#       6.425696 -8.455811 export_c01_sample_01_0_Day0_live_5
#       7.110587  8.610645 export_c01_sample_01_0_Day0_live_6
#       ...
#
#  2)  cytofkit-without-Bcl11a_Rphenograph_clusters.csv
#                                       name phenographCluster
#       1 export_c01_sample_01_0_Day0_live_1                 1
#       2 export_c01_sample_01_0_Day0_live_2                 4
#       3 export_c01_sample_01_0_Day0_live_3                 4
#       4 export_c01_sample_01_0_Day0_live_4                 4
#       5 export_c01_sample_01_0_Day0_live_5                11
#       6 export_c01_sample_01_0_Day0_live_6                 4
#       ...
#
# 3) ~/s/data/jeffRanish/fromMarjorieBrand/finalAnalysis/cytofkit-without-Bcl11a-June-2018-USED-FOR-PAPER/*.fcs
#    export_c01_sample_01_0_Day0_live.fcs - Day22
#    library(flowCore)
#    f <- "~/s/data/jeffRanish/fromMarjorieBrand/finalAnalysis/cytofkit-without-Bcl11a-June-2018-USED-FOR-PAPER/export_c01_sample_01_0_Day0_live.fcs"
#    x <- read.FCS(file.name, transformation=FALSE, alter.names=TRUE)
#    tbl.fcs <- pData(parameters(x))
#    mtx <- exprs(x)
#
# 4) a table of identifiers, constructed from the metadata in the fcs file.
#    is this information the same for every fcs file?
#
#
#----------------------------------------------------------------------------------------------------
library(RColorBrewer)

library(flowCore)
library(randomcoloR)
#----------------------------------------------------------------------------------------------------
# build a table of identifiers

f <- "~/s/data/jeffRanish/fromMarjorieBrand/finalAnalysis/cytofkit-without-Bcl11a-June-2018-USED-FOR-PAPER/export_c01_sample_01_0_Day0_live.fcs"
x <- read.FCS(f, transformation=FALSE, alter.names=TRUE)
tbl.fcs <- pData(parameters(x))  # 57 5
mtx <- exprs(x)
colnames(mtx) <- as.character(colnames(mtx))  # they are, weirdly, a list of lists

markers <- c("CD34", "CD36", "CD71", "CD38", "CD45RA", "CD123", "CD49f", "CD90", "CD44", "CD41", "CD235ab")
length(markers)
tfs <- c("GATA1", "PU.1", "ATRX", "c-Myc", "KLF1", "FLI1", "TAL1", "GATA2", "RUNX1", "NFE2p45", "BACH1",
         "IKZF1", "MAFG", "c-JUN", "KAT3B", "C-EBPa")
length(tfs)

x.markers <- as.list(sapply(markers, function(marker) as.character(grep(marker, tbl.fcs$desc, v=TRUE, ignore.case=TRUE))))
x.tfs <- as.list(sapply(tfs, function(tf) as.character(grep(tf, tbl.fcs$desc, v=TRUE, ignore.case=TRUE))))
x.both <- c(x.markers, x.tfs)
desc.ids <- as.character(x.both)
name.ids <- as.character(tbl.fcs$name[match(desc.ids, tbl.fcs$desc)])

tbl.ids <- data.frame(id=as.character(x.both), name=as.character(name.ids), stringsAsFactors=FALSE)   # 27 x 1
rownames(tbl.ids) <- names(x.both)

# head(tbl.ids)
#   tbl.ids
#                     id
#  CD34       149Sm_CD34
#  CD36       155Gd_CD36
#  CD71       175Lu_CD71
#  CD38       172Yb_CD38
#  CD45RA   143Nd_CD45RA
#  CD123     151Eu_CD123
#  CD49f     164Dy_CD49F
#  CD90       161Dy_CD90
#  CD44       153Eu_CD44
#  CD41         89Y_CD41
#  CD235ab 141Pr_CD235ab
#  GATA1     156Gd_GATA1
#  PU.1       167Er_PU.1



directory <- "~/s/data/jeffRanish/fromMarjorieBrand/finalAnalysis/cytofkit-without-Bcl11a-June-2018-USED-FOR-PAPER/Results.orig"
full.path <- file.path(directory, "cytofkit-without-Bcl11a_tsne_dimension_reduced_data.csv")
stopifnot(file.exists(full.path))

tbl <- read.table(full.path, sep=",",  header=TRUE, stringsAsFactors=FALSE)
colnames(tbl) <- c("sampleName", "tsne1", "tsne2")
tbl <- tbl[, c("tsne1", "tsne2", "sampleName")]

sampleInfo <-  strsplit(tbl$sampleName, "_", fixed=TRUE)

tbl.sampleInfo <- data.frame(t(sapply(sampleInfo,c)), stringsAsFactors=FALSE)
colnames(tbl.sampleInfo) <- c("export", "cluster", "ignore", "sample", "day", "dayN", "status", "number")
coi <- c("cluster", "dayN", "number")
tbl.sampleInfo <- tbl.sampleInfo[, coi]
tbl.sampleInfo$cluster <- as.integer(sub("^c", "", tbl.sampleInfo$cluster))
tbl.sampleInfo$dayN <- as.integer(sub("^Day", "", tbl.sampleInfo$dayN))
tbl.sampleInfo$number <- as.integer(tbl.sampleInfo$number)
colnames(tbl.sampleInfo) <- c("cluster", "day", "number")

tbl <- cbind(tbl.sampleInfo, tbl)

f2 <- file.path(directory, "cytofkit-without-Bcl11a_Rphenograph_clusters.csv")
stopifnot(file.exists(f2))

tbl.2 <- read.table(f2, sep=",",  header=TRUE, stringsAsFactors=FALSE)
colnames(tbl.2) <- c("name", "phenographCluster")
tbl$phenographCluster <- tbl.2$phenographCluster
clusterCount <- length(unique(tbl$phenographCluster))
printf("clusterCount: %d", clusterCount)

set.seed(17)
palette <- distinctColorPalette(clusterCount)
tbl$color <- palette[tbl$phenographCluster]
save(tbl, file="tsnePlot-allClusters.RData")
plot(tbl$tsne1 , tbl$tsne2, col=tbl$color)


# identify and color cells living on day 22
# sample name: "export_c14_sample_01_0_Day22_live.fcs"

id <- tbl.ids["CD34", "name"]
cd34.values <- mtx[, id]

tbl.day22 <- tbl[, c("tsne1", "tsne2", "day", "color")]
tbl.day22$color <- "#EEEEEE"   # default value
day22.indices <- which(tbl.day22$day==22)   # 82
length(cd34.values) == length(day22.indices)

#tbl.day22 <- subset(tbl, day==22)[, c("tsne1", "tsne2")]
#length(cd34.values) == nrow(tbl.day22)

red.colors <- brewer.pal(9, "Reds")[4:9]

#----------------------------------------------------------------------------------------------------
assignColor <- function(values, min, max, colors)
{
   numberOfBins <- length(colors) - 1
   bin.size <- (max-min)/numberOfBins
   bin.assignments <- 1 + as.numeric(lapply(values, function(value) floor(value/bin.size)))
   colors[bin.assignments]

} # assignColor
#----------------------------------------------------------------------------------------------------
day22.gene.colors <- assignColor(cd34.values, 0, 6, red.colors)
tbl.day22$color[day22.indices] <- day22.gene.colors

with(tbl.day22, plot(tsne1, tsne2, col=color))


# first entitity in tbl.ids is CD34, which goes by the name 149Sm_CD34
# what are the positions and values of CD34 on day 0?



fcs.directory <- "~/s/data/jeffRanish/fromMarjorieBrand/finalAnalysis/cytofkit-without-Bcl11a-June-2018-USED-FOR-PAPER"
fcs.files <- grep(".fcs$", list.files(fcs.directory), v=TRUE)


fcs.00 <- read.FCS(file.path(fcs.directory, "export_c01_sample_01_0_Day0_live.fcs"), transformation=FALSE, alter.names=TRUE)
fcs.20 <- read.FCS(file.path(fcs.directory, "export_c13_sample_01_0_Day20_live.fcs"), transformation=FALSE, alter.names=TRUE)
fcs.22 <- read.FCS(file.path(fcs.directory, "export_c14_sample_01_0_Day22_live.fcs"), transformation=FALSE, alter.names=TRUE)

mtx.00 <- exprs(fcs.00)
mtx.20 <- exprs(fcs.20)
mtx.22 <- exprs(fcs.22)

tbl.00 <- pData(parameters(fcs.00))
tbl.20 <- pData(parameters(fcs.20))
tbl.22 <- pData(parameters(fcs.22))

tbl.files <- data.frame(day=c(0,2,4,6,8,10,11,12,14,16,18,20,22),
                        file=fcs.files,
                        stringsAsFactors=FALSE)


total.rows <- 0
tbl.big <- data.frame()
tbls.days <- list()
for(r in seq_len(nrow(tbl.files))){
   filename <- tbl.files$file[r]
   dayNumber <- tbl.files$day[r]
   x <- read.FCS(file.path(fcs.directory, filename), transformation=FALSE, alter.names=TRUE)
   tbl.fcs <- pData(parameters(x))  # 57 5
   mtx <- exprs(x)
   colnames(mtx) <- as.character(colnames(mtx))  # they are, weirdly, a list of lists
   if(nrow(mtx) > 5000)
      mtx <- mtx[1:5000,]
   printf("%s matrix entries: %d", filename, nrow(mtx))
   total.rows <- total.rows + nrow(mtx)
   tbl.geom.day <- subset(tbl, day==dayNumber)
   #browser()
   printf("day %d", dayNumber)
   mtx <- mtx[, tbl.ids$name] # just the markers and tfs
   colnames(mtx) <- rownames(tbl.ids)
   tbl.day <- cbind(tbl.geom.day, mtx)
   tbls.days[[r]] <- tbl.day
   } # for filename

printf("grand total: %d", total.rows)
tbl.all <- do.call(rbind, tbls.days)
save(tbl.all, file="tbl.allDays.allProteins.RData")


