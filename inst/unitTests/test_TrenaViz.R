library(RUnit)
library(TrenaViz)
#------------------------------------------------------------------------------------------------------------------------
# create an instance for further tests
if(!exists("tv"))
   tv <- TrenaViz("TrenaProjectIGAP")
#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_constructor()
   test_createUI()


} # runTests
#------------------------------------------------------------------------------------------------------------------------
test_constructor <- function()
{
   printf("--- test_constructor")

   checkTrue("TrenaViz" %in% is(tv))

} # test_constructor
#------------------------------------------------------------------------------------------------------------------------
test_createUI <- function()
{
   printf("--- test_createServer")

   ui <- createUI(tv)
   checkTrue("shiny.tag" %in% is(ui))
   checkTrue(nchar(as.character(ui)) > 10000)

} # test_createUI
#------------------------------------------------------------------------------------------------------------------------
test_createServer <- function()
{
   printf("--- test_createServer")

   ui <- createUI(tv)
   checkTrue("shiny.tag" %in% is(ui))
   checkTrue(nchar(as.character(ui)) > 10000)

} # test_createUI
#------------------------------------------------------------------------------------------------------------------------
run <- function()
{
   createApp(tv)

} # run
#------------------------------------------------------------------------------------------------------------------------
