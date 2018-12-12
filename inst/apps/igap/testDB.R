library(RMariaDB)
driver <- RMariaDB::MariaDB()
host <- "genome-mysql.cse.ucsc.edu"
user <- "genome"
dbname <- "hg38"
db <- DBI::dbConnect(driver, user = user, host = host, dbname = dbname)
all.tableNames <- DBI::dbListTables(db);
head(all.tableNames)
d
encode.table.name <- "wgEncodeRegDnaseClustered"

main.clause <- sprintf("select * from %s where", encode.table.name);
chromosome <- "chr19"
start <- 1035007
end <- 1070666

query <- paste(main.clause,
                             sprintf("chrom = '%s'", chromosome),
                             sprintf("and chromStart >= %d", start),
                             sprintf("and chromEnd <= %d", end),
                             collapse = " ")

printf("query: %s", query)

tbl.regions <- dbGetQuery(db, query)
head(tbl.regions)
