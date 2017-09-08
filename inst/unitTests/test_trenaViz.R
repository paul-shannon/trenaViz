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
test_buildMultiModelGraph_oneModel <- function(display=FALSE)
{
   printf("--- test_buildMultiModelGraph_oneModel")

   load(system.file(package="trenaViz", "extdata", "sampleModelAndRegulatoryRegions.RData"))
   models <- list(tcf7=list(tbl.regulatoryRegions=tbl.regulatoryRegions.strong, tbl.geneModel=tbl.geneModel.strong))
   targetGene <- "TCF7"

   g <- buildMultiModelGraph(tv, targetGene, models)
   nodesInGraph <- nodes(g)
   regionNodes <- unique(models[[1]]$tbl.regulatoryRegions$id)
   tfNodes <- unique(models[[1]]$tbl.regulatoryRegions$geneSymbol)
   checkEquals(length(nodesInGraph), length(regionNodes) + length(tfNodes) + length(targetGene))
   tbl.reg <- models[[1]]$tbl.regulatoryRegions
   checkEquals(length(edgeNames(g)), nrow(tbl.reg) + length(unique(tbl.reg$id)))

   g.lo <- addGeneModelLayout(tv, g, xPos.span=1500)
   min.xPos <- min(as.numeric(nodeData(g.lo, attr="xPos")))
   max.xPos <- max(as.numeric(nodeData(g.lo, attr="xPos")))
   checkEquals(abs(max.xPos - min.xPos), 1500)

   if(display){
     browser()
     setGraph(tv, g.lo) #, names(models))
     setStyle(tv, system.file(package="trenaUtilities", "extdata", "style.js"))
     Sys.sleep(3); fit(tv)
     browser()
     xyz <- 99
     }

} # test_buildMultiModelGraph_oneModel
#------------------------------------------------------------------------------------------------------------------------
test_buildMultiModelGraph_fiveModels <- function(display=FALSE)
{
   printf("--- test_buildMultiModelGraph_fiveModels")
   targetGene <- "AQP4"
   aqp4.tss <- 26865884
   fp.source <- "postgres://whovian/brain_hint_20"
   sources <- list(fp.source)

   prep <- TrenaPrep(targetGene, aqp4.tss, "chr18", aqp4.tss-1000, aqp4.tss+1000, regulatoryRegionSources=sources)
   x <- getRegulatoryRegions(prep)
   closeAllPostgresConnections()
   tbl.regulatoryRegions <- expandRegulatoryRegionsTableByTF(prep, x[[fp.source]])

   tbl.geneModel <- createGeneModel(prep, "randomForest", tbl.regulatoryRegions, mtx)

      # two get multiple models, filter on randomForest score
   tbl.geneModel.rf10 <- subset(tbl.geneModel, randomForest > 10)
   tbl.regulatoryRegions.rf10 <- subset(tbl.regulatoryRegions, tf %in% tbl.geneModel.rf10$tf)

   tbl.geneModel.rf5 <- subset(tbl.geneModel, randomForest > 5)
   tbl.regulatoryRegions.rf5 <- subset(tbl.regulatoryRegions, tf %in% tbl.geneModel.rf5$tf)


   tbl.geneModel.rf3 <- subset(tbl.geneModel, randomForest > 3)
   tbl.regulatoryRegions.rf3 <- subset(tbl.regulatoryRegions, tf %in% tbl.geneModel.rf3$tf)

   tbl.geneModel.rf2 <- subset(tbl.geneModel, randomForest > 2)
   tbl.regulatoryRegions.rf2 <- subset(tbl.regulatoryRegions, tf %in% tbl.geneModel.rf2$tf)

   tbl.geneModel.rf1 <- subset(tbl.geneModel, randomForest > 1)
   tbl.regulatoryRegions.rf1 <- subset(tbl.regulatoryRegions, tf %in% tbl.geneModel.rf1$tf)

   models <- list(rf01=list(tbl.regulatoryRegions=tbl.regulatoryRegions.rf1,  tbl.geneModel=tbl.geneModel.rf1),
                  rf02=list(tbl.regulatoryRegions=tbl.regulatoryRegions.rf2,  tbl.geneModel=tbl.geneModel.rf2),
                  rf03=list(tbl.regulatoryRegions=tbl.regulatoryRegions.rf3,  tbl.geneModel=tbl.geneModel.rf3),
                  rf05=list(tbl.regulatoryRegions=tbl.regulatoryRegions.rf5,  tbl.geneModel=tbl.geneModel.rf5),
                  rf10=list(tbl.regulatoryRegions=tbl.regulatoryRegions.rf10, tbl.geneModel=tbl.geneModel.rf10)
                  )

   #save(models, file="testModel.RData")
   #load("testModel.RData")

   g <- buildMultiModelGraph(prep, models)
   nodesInGraph <- nodes(g)
   regionNodes <- unique(models[[1]]$tbl.regulatoryRegions$id)
   tfNodes <- unique(models[[1]]$tbl.regulatoryRegions$tf)
   checkEquals(length(nodesInGraph), length(regionNodes) + length(tfNodes) + length(targetGene))
   tbl.reg <- models[[1]]$tbl.regulatoryRegions
   #checkEquals(length(edgeNames(g)), nrow(tbl.reg) + length(unique(tbl.reg$id)))

   g.lo <- addGeneModelLayout(prep, g, xPos.span=1500)
   min.xPos <- min(as.numeric(nodeData(g.lo, attr="xPos")))
   max.xPos <- max(as.numeric(nodeData(g.lo, attr="xPos")))
   checkEquals(abs(max.xPos - min.xPos), 1500)

   if(display){
     httpAddGraph(tv, g.lo, names(models))
     loadStyle(tv, system.file(package="trenaUtilities", "extdata", "style.js"))
     Sys.sleep(3); fit(tv)
     browser()
     }

} # test_buildMultiModelGraph_fiveModels
#------------------------------------------------------------------------------------------------------------------------
test_buildMultiModelGraph_twoModels_15k_span <- function(display=FALSE)
{
   printf("--- test_buildMultiModelGraph_twoModels_15k_span")
   targetGene <- "AQP4"
   aqp4.tss <- 26865884
   fp.source <- "postgres://whovian/brain_hint_20"
   sources <- list(fp.source)

   prep <- TrenaPrep(targetGene, aqp4.tss, "chr18", aqp4.tss-5000, aqp4.tss+10000, regulatoryRegionSources=sources)
   x <- getRegulatoryRegions(prep)
   closeAllPostgresConnections()
   tbl.regulatoryRegions <- expandRegulatoryRegionsTableByTF(prep, x[[fp.source]])

   tbl.geneModel <- createGeneModel(prep, "randomForest", tbl.regulatoryRegions, mtx)

      # two get multiple models, filter on randomForest score
   tbl.geneModel.rf10 <- subset(tbl.geneModel, randomForest > 10)
   tbl.regulatoryRegions.rf10 <- subset(tbl.regulatoryRegions, tf %in% tbl.geneModel.rf10$tf)

   tbl.geneModel.rf1 <- subset(tbl.geneModel, randomForest > 1)
   tbl.regulatoryRegions.rf1 <- subset(tbl.regulatoryRegions, tf %in% tbl.geneModel.rf1$tf)


   models <- list(rf01=list(tbl.regulatoryRegions=tbl.regulatoryRegions.rf1,  tbl.geneModel=tbl.geneModel.rf1),
                  rf10=list(tbl.regulatoryRegions=tbl.regulatoryRegions.rf10, tbl.geneModel=tbl.geneModel.rf10)
                  )

   #save(models, file="models.2.big.RData")
   #load("models.2.big.RData")


   g <- buildMultiModelGraph(prep, models)
   nodesInGraph <- nodes(g)
   regionNodes <- unique(models[[1]]$tbl.regulatoryRegions$id)
   tfNodes <- unique(models[[1]]$tbl.regulatoryRegions$tf)
   checkEquals(length(nodesInGraph), length(regionNodes) + length(tfNodes) + length(targetGene))
   tbl.reg <- models[[1]]$tbl.regulatoryRegions
   #checkEquals(length(edgeNames(g)), nrow(tbl.reg) + length(unique(tbl.reg$id)))

   g.lo <- addGeneModelLayout(prep, g, xPos.span=1500)
   min.xPos <- min(as.numeric(nodeData(g.lo, attr="xPos")))
   max.xPos <- max(as.numeric(nodeData(g.lo, attr="xPos")))
   checkEquals(abs(max.xPos - min.xPos), 1500)

   if(display){
     httpAddGraph(tv, g.lo, names(models))
     loadStyle(tv, system.file(package="trenaUtilities", "extdata", "style.js"))
     Sys.sleep(3); fit(tv)
     browser()
     }

} # test_buildMultiModelGraph_fiveModels
#------------------------------------------------------------------------------------------------------------------------
test_geneModelLayout <- function()
{
   printf("--- test_geneModelLayout")
   filename <- "testModel.RData";
   if(!file.exists(filename)){
      printf("   skipping test_geneModelLayout, no file named '%s'", filename)
      return()
      }

   load(filename)
   print(load("ohsu.aqp4.graphLayoutNaNBug.input.RData"))
   # g <- buildMultiModelGraph(prep, models)
   targetGene <- "AQP4"
   aqp4.tss <- 26865884
   fp.source <- "postgres://whovian/brain_hint_20"
   sources <- list(fp.source)

   prep <- TrenaPrep(targetGene, aqp4.tss, "chr18", aqp4.tss-5000, aqp4.tss+10000, regulatoryRegionSources=sources)

   g.lo <- addGeneModelLayout(prep, g, xPos.span=1500)

} # test_geneModelLayout
#------------------------------------------------------------------------------------------------------------------------
if(!interactive())
    runTests()
