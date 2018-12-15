library(RColorBrewer)
print(load("~/github/cytofViz/inst/explore/tbl.allDays.allProteins.Rdata"))
with(tbl.all, plot(tsne1, tsne2, col=color))

#----------------------------------------------------------------------------------------------------
assignColor <- function(values, min, max, colors)
{
   numberOfBins <- length(colors) - 1
   bin.size <- (max-min)/numberOfBins
   bin.assignments <- 1 + as.numeric(lapply(values, function(value) floor(value/bin.size)))
   #browser()
   colors[bin.assignments]

} # assignColor
#----------------------------------------------------------------------------------------------------
red.colors <- c("#F0F0F0", brewer.pal(9, "Reds")[1:9])
all.proteins <- colnames(tbl.all)[9:35]
min <- min(asinh(tbl.all[, 9:35]))
max <- max(asinh(tbl.all[, 9:35]))
for(protein in all.proteins){
   vec <- asinh(tbl.all[, protein])
   colors <- assignColor(vec, min, max, red.colors)
   plot(tbl.all$tsne1, tbl.all$tsne2, col=colors, main=protein)
   Sys.sleep(2)
   }

hist(vec)
with(tbl.all, plot(tsne1, tsne2, col=color))






