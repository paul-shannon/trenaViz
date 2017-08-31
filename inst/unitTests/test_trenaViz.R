library(RUnit)
library(trenaViz)
#--------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#--------------------------------------------------------------------------------
PORT.RANGE <- 8000:8020
tv <- trenaViz(PORT.RANGE, quiet=TRUE);
#--------------------------------------------------------------------------------
runTests <- function()
{
  testConstructor();
  testWindowTitle()
  testPing()
  testIGV()
  testGraph()
  testLoadAndRemoveTracks()

  #closeWebSocket(tv)

} # runTests
#--------------------------------------------------------------------------------
testConstructor <- function()
{
   printf("--- testConstructor")
   checkTrue(ready(tv))
   checkTrue(port(tv) %in% PORT.RANGE)

} # testConstructor
#--------------------------------------------------------------------------------
testWindowTitle <- function()
{
   printf("--- testWindowTitle")
   checkTrue(ready(tv))
   setBrowserWindowTitle(tv, "trenaViz")
   checkEquals(getBrowserWindowTitle(tv), "trenaViz")
   setBrowserWindowTitle(tv, "new title");
   checkEquals(getBrowserWindowTitle(tv), "new title")

} # testWindowTitle
#--------------------------------------------------------------------------------
testPing <- function()
{
   printf("--- testPing")
   checkTrue(ready(tv))
   checkEquals(ping(tv), "pong")

} # testPing
#--------------------------------------------------------------------------------
testIGV <- function()
{
   printf("--- testIGV")
   setGenome(tv, "hg38")
   Sys.sleep(5);
   showGenomicRegion(tv, "AQP4")
   Sys.sleep(5);
   chromLocString <- getGenomicRegion(tv)
   checkTrue(grepl("chr18:", chromLocString));

} # testIGV
#--------------------------------------------------------------------------------
testGraph <- function()
{
   printf("--- testGraph")
   setGraph(tv);
   Sys.sleep(2);
   checkEquals(length(getSelectedNodes(tv)), 0)
   selectNodes(tv, "a")
   Sys.sleep(2)
   checkEquals(getSelectedNodes(tv), "a")

} # testGraph
#--------------------------------------------------------------------------------
testLoadAndRemoveTracks <- function()
{
   printf("--- testLoadAndRemoveTracks")

   raiseTab(tv, "IGV")

   setBrowserWindowTitle(tv, "test load tracks")
   checkEquals(getTrackNames(tv), "Gencode v24")

   segments <- 5
   set.seed(17);
   starts <- 26860646 + round(runif(segments,3,300))
   ends   <- starts + round(runif(segments,3,8))
   tbl.bed <- data.frame(chr=rep("18", segments),
                         start=starts,
                         end=ends,
                         name=LETTERS[1:segments],
                         score=runif(segments, -1, 1),
                         stringsAsFactors=FALSE)

   addBedTrackFromDataFrame(tv, "tbl.bed", tbl.bed, displayMode="EXPANDED", color="darkRed")
   showGenomicRegion(tv, sprintf("chr18:%d-%d", min(tbl.bed$start) - 10, max(tbl.bed$end) + 10))

   tbl.bedGraph <- tbl.bed[, c(1,2,3,5,4)]
   addBedGraphTrackFromDataFrame(tv, "tbl.bedGraph", tbl.bedGraph, color="darkGreen",
                                 minValue=min(tbl.bedGraph$score), maxValue=max(tbl.bedGraph$score),
                                 displayMode="EXPANDED")

   checkEquals(sort(getTrackNames(tv)), c("Gencode v24",  "tbl.bed", "tbl.bedGraph"))
   Sys.sleep(3)
   removeTracksByName(tv, c("tbl.bed", "tbl.bedGraph"))
   checkEquals(getTrackNames(tv), "Gencode v24")

} # testLoadAndRemoveTracks
#--------------------------------------------------------------------------------
if(!interactive())
    runTests()
