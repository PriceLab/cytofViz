library(flowCore)

directory <- "~/s/data/jeffRanish/fromMarjorieBrand/finalAnalysis/cytofkit-without-Bcl11a-June-2018-USED-FOR-PAPER"
fcs.files <- grep("fcs$", list.files(directory), value=TRUE)
stopifnot(length(fcs.files) == 13)
count <- length(fcs.files)
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
