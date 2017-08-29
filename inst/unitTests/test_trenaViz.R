library(RUnit)
library(trenaViz)
#--------------------------------------------------------------------------------
PORT.RANGE <- 8000:8020
#--------------------------------------------------------------------------------
runTests <- function()
{
  testConstructor();
  testWindowTitle()
  testPing()

} # runTests
#--------------------------------------------------------------------------------
testConstructor <- function()
{
   print("--- testConstructor")
   tv <- trenaViz(PORT.RANGE, quiet=FALSE);
   checkTrue(ready(tv))
   checkTrue(port(tv) %in% PORT.RANGE)
   closeWebSocket(tv)

} # testConstructor
#--------------------------------------------------------------------------------
testWindowTitle <- function()
{
   print("--- testWindowTitle")
   tv <- trenaViz(PORT.RANGE)
   checkTrue(ready(tv))
   checkEquals(getBrowserWindowTitle(tv), "trenaViz")
   setBrowserWindowTitle(tv, "new title");
   checkEquals(getBrowserWindowTitle(tv), "new title")
   closeWebSocket(tv)

} # testWindowTitle
#--------------------------------------------------------------------------------
testPing <- function()
{
   print("--- testPing")
   tv <- trenaViz(PORT.RANGE)
   checkTrue(ready(tv))
   checkEquals(ping(tv), "pong")
   closeWebSocket(tv)

} # testPing
#--------------------------------------------------------------------------------
testIGV <- function()
{
   print("--- testIGV")
   tv <- trenaViz(PORT.RANGE)
   setGenome(tv, "hg38")
   Sys.sleep(5);
   showGenomicRegion(tv, "AQP4")
   Sys.sleep(5);
   chromLocString <- getGenomicRegion(tv)
   checkTrue(grepl("chr18:", chromLocString));

} # testIGV
#--------------------------------------------------------------------------------
if(!interactive())
    runTests()
