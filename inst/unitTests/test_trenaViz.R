library(RUnit)
library(trenaViz)
#------------------------------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#------------------------------------------------------------------------------------------------------------------------
PORT.RANGE <- 8000:8020
if(!exists("tv"))
   tv <- trenaViz(PORT.RANGE, quiet=TRUE);
#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
  testConstructor();
  testWindowTitle()
  testPing()
  testIGV()
  test_graphToJSON()
  #testGraph()
  testLoadAndRemoveTracks()

  test_buildMultiModelGraph_oneModel_twoRandomScoresOnly()
  test_buildMultiModelGraph_oneModel_allScores()

  test_buildMultiModelGraph_oneModel()
  test_buildMultiModelGraph_twoModels()

  #closeWebSocket(tv)

} # runTests
#------------------------------------------------------------------------------------------------------------------------
testConstructor <- function()
{
   printf("--- testConstructor")
   checkTrue(ready(tv))
   checkTrue(port(tv) %in% PORT.RANGE)

} # testConstructor
#------------------------------------------------------------------------------------------------------------------------
testWindowTitle <- function()
{
   printf("--- testWindowTitle")
   checkTrue(ready(tv))
   setBrowserWindowTitle(tv, "trenaViz")
   checkEquals(getBrowserWindowTitle(tv), "trenaViz")
   setBrowserWindowTitle(tv, "new title");
   checkEquals(getBrowserWindowTitle(tv), "new title")

} # testWindowTitle
#------------------------------------------------------------------------------------------------------------------------
testPing <- function()
{
   printf("--- testPing")
   checkTrue(ready(tv))
   checkEquals(ping(tv), "pong")

} # testPing
#------------------------------------------------------------------------------------------------------------------------
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
#------------------------------------------------------------------------------------------------------------------------
test_graphToJSON <- function()
{
   printf("--- test_graphToJSON")

   load(system.file(package="trenaViz", "extdata", "sampleModelAndRegulatoryRegions.RData"))

     # create a small graph, with just two TFs, 2 regulatory regions
   tbl.model <- head(tbl.geneModel.strong, n=2)[, c("gene", "pearson.coeff", "rf.score")]
   tbl.reg   <- subset(tbl.regulatoryRegions.strong, geneSymbol %in% tbl.model$gene)[1:2,]

     # be sure that there is at least one regulatory region (and motif) per TF
   checkTrue(all(tbl.model$gene %in% tbl.reg$geneSymbol))

   targetGene <- "TCF7"
   model <- list(tcf7=list(model=tbl.model, regions=tbl.reg))
   g <- buildMultiModelGraph(tv, targetGene, model)

   xCoordinate.span <- 1500
   g.lo <- addGeneModelLayout(tv, g, xPos.span=xCoordinate.span)
   g.json <- trenaViz:::.graphToJSON(g.lo)

   checkEquals(class(g.json), "character")
   checkTrue(nchar(g.json) > 3000)

   x <- fromJSON(g.json)
   checkEquals(names(x), "elements")
   x.el <- x$elements
   checkEquals(sort(names(x.el)), c("data", "position"))
   checkTrue(is.data.frame(x.el$data))
   checkTrue(is.data.frame(x.el$position))

     # these back-converted data.frames are a bit odd in that nodes and edges and their attribute
     # are all in the same data.frame.  we can, however, separate them out by ad hoc means

   node.rows <- which(!is.na(x.el$data$type))
   edge.rows <- which(is.na(x.el$data$type))

   tbl.nodes <- x.el$data[node.rows,]
   tbl.edges <- x.el$data[edge.rows,]

      # do some very simple tests on the nodes
   checkEquals(tbl.nodes$label, c("TCF7", "LEF1", "FOXP1", "MA0523.1", "MA0593.1"))
   checkEquals(tbl.nodes$type,  c("targetGene", "TF", "TF", "regulatoryRegion", "regulatoryRegion"))

     # and pm the edges
   checkEquals(tbl.edges$edgeType, c("bindsTo", "bindsTo", "regulatorySiteFor", "regulatorySiteFor"))
   checkEquals(tbl.edges$id,
               c("LEF1->TCF7.fp.upstream.000351.Hsapiens-jaspar2016-TCF7L2-MA0523.1",
                 "FOXP1->TCF7.fp.downstream.000329.Hsapiens-jaspar2016-FOXP2-MA0593.1",
                 "TCF7.fp.upstream.000351.Hsapiens-jaspar2016-TCF7L2-MA0523.1->TCF7",
                 "TCF7.fp.downstream.000329.Hsapiens-jaspar2016-FOXP2-MA0593.1->TCF7"))

     # now check the node positions for plausibility
   tbl.pos <- x.el$position[node.rows,]
   checkEquals(colnames(tbl.pos), c("x", "y"))

     # the layout algorithm ensures that the x coordiantes are centered on zero
     # and that the whole span is not larger than that specified above as "xCoordinate.span"
   checkTrue(all(tbl.pos$x) > -(xCoordinate.span/2))
   checkTrue(all(tbl.pos$x) <  (xCoordinate.span/2))

     # abs(y) are >= 1000
   checkTrue(all(abs(tbl.pos$y) <= 2000))

} # test_graphToJSON
#------------------------------------------------------------------------------------------------------------------------
no_testGraph <- function()
{
   printf("--- testGraph")
   setGraph(tv);
   Sys.sleep(2);
   checkEquals(length(getSelectedNodes(tv)), 0)
   selectNodes(tv, "a")
   Sys.sleep(2)
   checkEquals(getSelectedNodes(tv), "a")

} # no_testGraph
#------------------------------------------------------------------------------------------------------------------------
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

   tbl.bedGraph2 <- tbl.bedGraph
   colnames(tbl.bedGraph2)[1] <- "Chrom"   # make sure this colname tolerates case & extra characters
   addBedGraphTrackFromDataFrame(tv, "tbl.bedGraph2", tbl.bedGraph2, color="purple",
                                 minValue=min(tbl.bedGraph2$score), maxValue=max(tbl.bedGraph2$score),
                                 displayMode="EXPANDED")

   checkEquals(sort(getTrackNames(tv)), c("Gencode v24", "tbl.bedGraph2"))



} # testLoadAndRemoveTracks
#------------------------------------------------------------------------------------------------------------------------
test_buildMultiModelGraph_oneModel <- function(display=FALSE)
{
   printf("--- test_buildMultiModelGraph_oneModel")

   load(system.file(package="trenaViz", "extdata", "tcf7Model.Rdata"))
   colnames(tbl.model) <- gsub(".", "", colnames(tbl.model), fixed=TRUE)
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
# make sure that the build process does not assume the existence of any one node attribute
test_buildMultiModelGraph_oneModel_twoRandomScoresOnly <- function(display=FALSE)
{
   printf("--- test_buildMultiModelGraph_oneModel_twoRandomScoresOnly")

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
# make sure that the build process does not assume the existence of any one node attribute
test_buildMultiModelGraph_oneModel_allScores <- function(display=FALSE)
{
   printf("--- test_buildMultiModelGraph_oneModel_allScores")

   load(system.file(package="trenaViz", "extdata", "tcf7Model.Rdata"))

   scores.to.keep <- colnames(tbl.model)[-1]
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

} # test_buildMultiModelGraph_oneModel_allScores
#------------------------------------------------------------------------------------------------------------------------
test_buildMultiModelGraph_twoModels <- function(display=FALSE)
{
   printf("--- test_buildMultiModelGraph_twoModels")

   load(system.file(package="trenaViz", "extdata", "sampleModelAndRegulatoryRegions.RData"))

   #scores.to.keep <- colnames(tbl.geneModel.strong)[1+sample(1:(ncol(tbl.geneModel.strong)-1), 2)]
   #columns.to.keep <- c("gene", scores.to.keep)
   #tbl.model.chopped <- tbl.geneModel.strong[, columns.to.keep]

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

      # there should be 3 instances of every node attribute.  using "concordance" as an example:
      #    "concordance", "tcf7.concordance", "arnt2.deleted.concordance".
   attributes <- colnames(tbl.geneModel.strong)
   attributes <- attributes[-match("gene", attributes)]

      # check for the three forms, for now ignoring prefixes ("", "tcf7.", "arnt2.deleted.")
   checkTrue(all(unlist(lapply(attributes, function(attribute)
                                  length(grep(attribute, names(nodeDataDefaults(g)))) == 3))))

     # check for "tcf7." attributes
   checkTrue(all(unlist(lapply(attributes, function(attribute)
                                  length(grep(sprintf("tcf7.%s", attribute), names(nodeDataDefaults(g)))) == 1))))

     # check for the "arnt2.deleted." attributes
   checkTrue(all(unlist(lapply(attributes, function(attribute)
                                  length(grep(sprintf("arnt2.deleted.%s", attribute), names(nodeDataDefaults(g)))) == 1))))

     # check for the attributes without prefixes
   checkTrue(all(unlist(lapply(attributes, function(attribute)
                                  length(grep(sprintf("^%s", attribute), names(nodeDataDefaults(g)))) == 1))))

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
