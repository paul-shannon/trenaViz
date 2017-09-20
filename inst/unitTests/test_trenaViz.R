library(RUnit)
library(trenaViz)
#--------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#--------------------------------------------------------------------------------
PORT.RANGE <- 8000:8020
if(!exists("tv"))
   tv <- trenaViz(PORT.RANGE, quiet=TRUE);
#--------------------------------------------------------------------------------
runTests <- function()
{
  testConstructor();
  testWindowTitle()
  testPing()
  testIGV()
  test_graphToJSON()
  #testGraph()
  testLoadAndRemoveTracks()
  test_buildMultiModelGraph_oneModel()
  test_buildMultiModelGraph_twoModels()

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
test_graphToJSON <- function()
{
   load(system.file(package="trenaViz", "extdata", "graph.11nodes.14edges.RData"))
   g.json <- trenaViz:::.graphToJSON(g)

} # test_graphToJSON
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

   load(system.file(package="trenaViz", "extdata", "tcf7Model.Rdata"))
   models <- list(tcf7=list(regions=tbl.regRegions, model=tbl.model))
   targetGene <- "TCF7"

   g <- buildMultiModelGraph(tv, targetGene, models)
   nodesInGraph <- nodes(g)
   regionNodes <- unique(models[[1]]$regions$id)
   tfNodes <- unique(models[[1]]$regions$geneSymbol)

   checkEquals(length(nodesInGraph), length(regionNodes) + length(tfNodes) + length(targetGene))
   tbl.reg <- models[[1]]$regions
   checkEquals(length(edgeNames(g)), nrow(tbl.reg) + length(unique(tbl.reg$id)))

   checkTrue(!any(is.null(unlist(nodeData(g), use.names=FALSE))))

   g.lo <- addGeneModelLayout(tv, g, xPos.span=1500)
   min.xPos <- min(as.numeric(nodeData(g.lo, attr="xPos")))
   max.xPos <- max(as.numeric(nodeData(g.lo, attr="xPos")))
   checkEquals(abs(max.xPos - min.xPos), 1500)

   if(display){
     browser()
     setGraph(tv, g.lo, names(models))
     setStyle(tv, system.file(package="trenaUtilities", "extdata", "style.js"))
     Sys.sleep(3); fit(tv)
     browser()
     xyz <- 99
     }

} # test_buildMultiModelGraph_oneModel
#------------------------------------------------------------------------------------------------------------------------
# make sure that the build process does not assume the existence of any one
# node attribute
test_buildMultiModelGraph_oneModel_twoRandomScoresOnly <- function(display=FALSE)
{
   printf("--- test_buildMultiModelGraph_oneModel_twoRandomeScoresOnly")

   load(system.file(package="trenaViz", "extdata", "tcf7Model.Rdata"))

   scores.to.keep <- colnames(tbl.model)[1+sample(1:(ncol(tbl.model)-1), 2)]
   columns.to.keep <- c("gene", scores.to.keep)

   tbl.model.chopped <- tbl.model[, columns.to.keep]
   models <- list(tcf7=list(regions=tbl.regRegions, model=tbl.model.chopped))
   targetGene <- "TCF7"

   g <- buildMultiModelGraph(tv, targetGene, models)
   nodesInGraph <- nodes(g)
   regionNodes <- unique(models[[1]]$regions$id)
   tfNodes <- unique(models[[1]]$regions$geneSymbol)

   checkEquals(length(nodesInGraph), length(regionNodes) + length(tfNodes) + length(targetGene))
   tbl.reg <- models[[1]]$regions
   checkEquals(length(edgeNames(g)), nrow(tbl.reg) + length(unique(tbl.reg$id)))

      # check the attributes
      # keep in mind that their are the real node attributes - whose current values
      # control the cyjs display - and the condition-specific attributes, which
      # are the source for the real ones.  in cyjs, as the user flips from one
      # visual model to the next, the condition-specific attributes are copied
      # into the real ones.

   checkTrue(!any(is.null(unlist(nodeData(g), use.names=FALSE))))

   noa.names <- sort(names(nodeDataDefaults(g)))
   checkTrue(all(scores.to.keep %in% noa.names))
   checkTrue(all(sprintf("tcf7.%s", scores.to.keep) %in% noa.names))
   standard.attributes <- c("distance", "label", "motif", "motifInModel", "type", "xPos", "yPos")
   checkTrue(all(standard.attributes %in% noa.names))
   checkTrue(all(sprintf("tcf7.%s", standard.attributes) %in% noa.names))

   g.lo <- addGeneModelLayout(tv, g, xPos.span=1500)
   min.xPos <- min(as.numeric(nodeData(g.lo, attr="xPos")))
   max.xPos <- max(as.numeric(nodeData(g.lo, attr="xPos")))
   checkEquals(abs(max.xPos - min.xPos), 1500)

   if(display){
     browser()
     setGraph(tv, g.lo, names(models))
     setStyle(tv, system.file(package="trenaUtilities", "extdata", "style.js"))
     Sys.sleep(3); fit(tv)
     browser()
     xyz <- 99
     }

} # test_buildMultiModelGraph_oneModel_twoRandomScoresOnly
#------------------------------------------------------------------------------------------------------------------------
test_buildMultiModelGraph_twoModels <- function(display=FALSE)
{
   printf("--- test_buildMultiModelGraph_twoModels")

   load(system.file(package="trenaViz", "extdata", "sampleModelAndRegulatoryRegions.RData"))

     # copy and modify the gene model slightly but noticeably: remove the large repressor, ARNT2
   tbl.geneModel.2 <- subset(tbl.geneModel.strong, gene != "ARNT2")
   corresponding.motifs <- grep("ARNT2", tbl.regulatoryRegions.strong$geneSymbol)
   stopifnot(length(corresponding.motifs) >= 1)
   tbl.regulatoryRegions.strong.2 <- tbl.regulatoryRegions.strong[-corresponding.motifs,]

   models <- list(tcf7=list(regions=tbl.regulatoryRegions.strong,
                            model=tbl.geneModel.strong),
                  arnt2.deleted=list(regions=tbl.regulatoryRegions.strong.2,
                                     model=tbl.geneModel.2))

   targetGene <- "TCF7"

   g <- buildMultiModelGraph(tv, targetGene, models)
   nodesInGraph <- nodes(g)
   regionNodes <- unique(models[[1]]$regions$id)
   tfNodes <- unique(models[[1]]$regions$geneSymbol)
   checkEquals(length(nodesInGraph), length(regionNodes) + length(tfNodes) + length(targetGene))
   tbl.reg <- models[[1]]$regions
   checkEquals(length(edgeNames(g)), nrow(tbl.reg) + length(unique(tbl.reg$id)))


      #--------------------------------------------------------------------------------
      # has the single ARNT2-related motif, Mmusculus-jaspar2016-Ahr::Arnt-MA0006.1
      # been marked as in model 1, but absent from model 2?
      #--------------------------------------------------------------------------------

   checkEquals(as.list(table(as.logical(nodeData(g, attr="arnt2.deleted.motifInModel"))))$`TRUE`,  51)
   checkEquals(as.list(table(as.logical(nodeData(g, attr="arnt2.deleted.motifInModel"))))$`FALSE`, 1)

   arnt2.motif.node <- "TCF7.fp.upstream.000410.Mmusculus-jaspar2016-Ahr::Arnt-MA0006.1"
   checkEquals(nodeData(g, attr="arnt2.deleted.motifInModel", n=arnt2.motif.node)[[1]], FALSE)

      # the tcf7 model, however, includes all motifs:

   checkTrue(all(as.logical(nodeData(g, attr="tcf7.motifInModel"))))

   g.lo <- addGeneModelLayout(tv, g, xPos.span=1500)
   min.xPos <- min(as.numeric(nodeData(g.lo, attr="xPos")))
   max.xPos <- max(as.numeric(nodeData(g.lo, attr="xPos")))
   checkEquals(abs(max.xPos - min.xPos), 1500)

   if(display){
     setGraph(tv, g.lo, names(models))
     setStyle(tv, system.file(package="trenaUtilities", "extdata", "style.js"))
     Sys.sleep(3); fit(tv)
     browser()
     xyz <- 99
     }

} # test_buildMultiModelGraph_twoModels
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

   models <- list(rf01=list(regions=tbl.regulatoryRegions.rf1,  model=tbl.geneModel.rf1),
                  rf02=list(regions=tbl.regulatoryRegions.rf2,  model=tbl.geneModel.rf2),
                  rf03=list(regions=tbl.regulatoryRegions.rf3,  model=tbl.geneModel.rf3),
                  rf05=list(regions=tbl.regulatoryRegions.rf5,  model=tbl.geneModel.rf5),
                  rf10=list(regions=tbl.regulatoryRegions.rf10, model=tbl.geneModel.rf10)
                  )

   #save(models, file="testModel.RData")
   #load("testModel.RData")

   g <- buildMultiModelGraph(prep, models)
   nodesInGraph <- nodes(g)
   regionNodes <- unique(models[[1]]$regions$id)
   tfNodes <- unique(models[[1]]$regions$tf)
   checkEquals(length(nodesInGraph), length(regionNodes) + length(tfNodes) + length(targetGene))
   tbl.reg <- models[[1]]$regions
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


   models <- list(rf01=list(regions=tbl.regulatoryRegions.rf1,  model=tbl.geneModel.rf1),
                  rf10=list(regions=tbl.regulatoryRegions.rf10, model=tbl.geneModel.rf10)
                  )

   #save(models, file="models.2.big.RData")
   #load("models.2.big.RData")


   g <- buildMultiModelGraph(prep, models)
   nodesInGraph <- nodes(g)
   regionNodes <- unique(models[[1]]$regions$id)
   tfNodes <- unique(models[[1]]$regions$tf)
   checkEquals(length(nodesInGraph), length(regionNodes) + length(tfNodes) + length(targetGene))
   tbl.reg <- models[[1]]$regions
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
