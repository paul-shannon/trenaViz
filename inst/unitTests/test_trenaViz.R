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
  testLoadTracks()

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
testLoadTracks <- function()
{
   printf("--- testLoadTracks")

   setBrowserWindowTitle(tv, "test load tracks")
   segments <- 5
   starts <- 26860646 + round(runif(segments,3,30))
   ends   <- starts + round(runif(segments,3,30))
   tbl.bed <- data.frame(chrom=rep("18", segments),
                         start=starts,
                         end=ends,
                         name=LETTERS[1:segments],
                         score=runif(segments, -1, 1),
                         stringsAsFactors=FALSE)

   addBedTrackFromDataFrame(tv, sprintf("data.frame"), tbl.bed, displayMode="EXPANDED", color="darkRed")
   showGenomicRegion(tv, sprintf("chr18:%d-%d", min(tbl.bed$start) - 10, max(tbl.bed$end) + 10))
   #addBedTrackFromHostedFile,
   #addBedGraphTrackFromDataFrame


} # testLoadTracks
#--------------------------------------------------------------------------------
if(!interactive())
    runTests()
