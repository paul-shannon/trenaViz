library(RUnit)
library(trena)
library(trenaViz)
library(MotifDb)
#------------------------------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#------------------------------------------------------------------------------------------------------------------------
PORT.RANGE <- 8000:8020
if(!exists("tv"))
   tv <- trenaViz(PORT.RANGE, quiet=TRUE);
#------------------------------------------------------------------------------------------------------------------------
if(!exists("mtx")){
    load(system.file(package="trena", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
    mtx <- asinh(mtx.sub)
    mtx.var <- apply(mtx, 1, var)
    deleters <- which(mtx.var < 0.01)
    if(length(deleters) > 0)   # 15838 x 638
        mtx <- mtx[-deleters,]
    }
#------------------------------------------------------------------------------------------------------------------------
runTests <- function(display=FALSE)
{
  testConstructor();
  testWindowTitle()
  testPing()

  testIGV()
  testLoadAndRemoveTracks()

  test_graphToJSON()

  test_buildMultiModelGraph_oneModel_twoRandomScoresOnly()
  test_buildMultiModelGraph_oneModel_allScores()
  test_buildMultiModelGraph_oneModel()
  test_buildMultiModelGraph_twoModels()

  test_buildMultiModelGraph_twoModels_two_promoterSpans(display=TRUE)

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

      # back-convert the json for integrity checking
   x <- fromJSON(g.json, flatten=TRUE)
   checkEquals(names(x), "elements")
   x.el <- x$elements
   tbl.nodes <- x.el$nodes
        # check a few of the many columns
   some.expected.columns <- c("data.id", "data.label", "data.distance", "data.xPos", "position.x",
                              "data.type", "data.tcf7.pearson.coeff")
   checkTrue(all(some.expected.columns %in% colnames(tbl.nodes)))
   expected.nodes <- c("TCF7", "LEF1", "FOXP1", "TCF7.fp.upstream.000351.Hsapiens-jaspar2016-TCF7L2-MA0523.1",
                       "TCF7.fp.downstream.000329.Hsapiens-jaspar2016-FOXP2-MA0593.1")
   checkTrue(all(expected.nodes %in% tbl.nodes$data.id))


   tbl.edges <- x.el$edges
   checkEquals(colnames(tbl.edges), c("data.id", "data.source", "data.target", "data.edgeType"))
   checkEquals(sort(unique(tbl.edges$data.edgeType)), c("bindsTo", "regulatorySiteFor"))

     # now check the node positions for plausibility
   checkTrue(all(tbl.nodes$position.x > -(xCoordinate.span/2)))
   checkTrue(all(tbl.nodes$position.x <  (1.2 * (xCoordinate.span/2))))   # whis is this fudge factor needed?

   checkTrue(all(abs(tbl.nodes$position.y) < 2000))

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
   colnames(tbl.model) <- gsub(".", "_", colnames(tbl.model), fixed=TRUE)
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

   colnames(tbl.model) <- gsub(".", "_", colnames(tbl.model), fixed=TRUE)

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

   colnames(tbl.geneModel) <- gsub(".", "_", colnames(tbl.geneModel), fixed=TRUE)
   colnames(tbl.geneModel.strong) <- gsub(".", "_", colnames(tbl.geneModel.strong), fixed=TRUE)

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
     setStyle(tv, "style.js")
     #setStyle(tv, system.file(package="trenaUtilities", "extdata", "style.js"))
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

   trena <- Trena("hg38")
   tbl.geneModel <- createGeneModel(trena, targetGene, solver.names, tbl.motifs.tfs, mtx)

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
test_buildMultiModelGraph_twoModels_two_promoterSpans <- function(display=FALSE)
{
   printf("--- test_buildMultiModelGraph_twoModels_two_promotersSpans")

   #tv <- trenaViz(PORT.RANGE)
   setGenome(tv, "hg38")
   targetGene <- "MEF2C"

   raiseTab(tv, "IGV")
   showGenomicRegion(tv, targetGene)
   targetGene.tss <- 88904257
   chromosome <- "chr5"

   database.filename <- system.file(package="trena", "extdata", "mef2c.neigborhood.hg38.footprints.db")
   database.uri <- sprintf("sqlite://%s", database.filename)
   sources <- list(sqlite=database.uri)

   trena <- Trena("hg38")

       #----------------------------------------------------------------------------------
       # create one gene model with a traditional, conservative promoter
       #----------------------------------------------------------------------------------

   loc.start <- targetGene.tss - 200
   loc.end   <- targetGene.tss + 2000

   showGenomicRegion(tv, sprintf("%s:%d-%d", chromosome, loc.start, loc.end))
   x <- getRegulatoryChromosomalRegions(trena, chromosome, loc.start, loc.end, sources,
                                        targetGene, targetGene.tss)

   names(x) <- names(sources)
   tbl.fp1 <- x[[1]]

   if(display){
     tbl.bed <- unique(tbl.fp1[, 1:3])
     addBedTrackFromDataFrame(tv, "brain fp1", tbl.bed, displayMode="EXPANDED", color="darkRed")
     }

   tbl.fp1 <- associateTranscriptionFactors(MotifDb, tbl.fp1, source="MotifDb", expand.rows=TRUE)

   unmapped <- which(is.na(tbl.fp1$geneSymbol))
   if(length(unmapped) > 0)
      tbl.fp1 <- tbl.fp1[-unmapped,]

   solver.names <- c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman")
   tbl.geneModel.1 <- createGeneModel(trena, targetGene, solver.names, tbl.fp1,  mtx)

       #----------------------------------------------------------------------------------
       # create a second gene regulatory model with a broader putative promoeter
       #----------------------------------------------------------------------------------

   loc.start <- targetGene.tss - 2000
   loc.end   <- targetGene.tss + 5000
      # there will be one data.frame of regulatory regions for every supplied data source
   showGenomicRegion(tv, sprintf("%s:%d-%d", chromosome, loc.start, loc.end))
   x <- getRegulatoryChromosomalRegions(trena, chromosome, loc.start, loc.end, sources,
                                        targetGene, targetGene.tss)
   names(x) <- names(sources)
   tbl.fp2 <- x[[1]]

   if(display){
     tbl.bed <- unique(tbl.fp2[, 1:3])
     addBedTrackFromDataFrame(tv, "brain fp2", tbl.bed, displayMode="EXPANDED", color="darkBlue")
     }


   tbl.fp2 <- associateTranscriptionFactors(MotifDb, tbl.fp2, source="MotifDb", expand.rows=TRUE)

   unmapped <- which(is.na(tbl.fp2$geneSymbol))
   if(length(unmapped) > 0)
      tbl.fp2 <- tbl.fp2[-unmapped,]

   solver.names <- c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman")
   tbl.geneModel.2 <- createGeneModel(trena, targetGene, solver.names, tbl.fp2,  mtx)

       #----------------------------------------------------------------------------------
       # before making a "multimodel", delete all the candidate regulatory regions
       # except those corresponding to TFs in the corresponding gene model
       #----------------------------------------------------------------------------------

   tbl.fp1 <- unique(subset(tbl.fp1, geneSymbol %in% tbl.geneModel.1$gene))    # 4  x 11
   tbl.fp2 <- unique(subset(tbl.fp2, geneSymbol %in% tbl.geneModel.2$gene))    # 21 x 11

       #----------------------------------------------------------------------------------
       # jsavascript (cytoscape.js in the trenaViz webapp) cannot use "." in
       # variable names: they are interpreted as record field delimeters
       # our impefrect solution, for now, is to convert all of these to underscore
       # in the colum titles
       #----------------------------------------------------------------------------------

   colnames(tbl.geneModel.1) <- gsub(".", "_", colnames(tbl.geneModel.1), fixed=TRUE)
   colnames(tbl.geneModel.2) <- gsub(".", "_", colnames(tbl.geneModel.2), fixed=TRUE)

   models <- list(promoter_2000_200 =list(regions=tbl.fp1,   model=tbl.geneModel.1),
                  promoter_5000_2000=list(regions=tbl.fp2,   model=tbl.geneModel.2)
                  )

   g <- buildMultiModelGraph(tv, targetGene, models)
   nodesInGraph <- nodes(g)
   regionNodes <- unique(c(models[[1]]$regions$id, models[[2]]$regions$id))
   tfNodes     <- unique(c(models[[1]]$model$gene, models[[2]]$model$gene))
   checkEquals(length(nodesInGraph), length(regionNodes) + length(tfNodes) + length(targetGene))

   g.lo <- addGeneModelLayout(tv, g, xPos.span=1500)
   min.xPos <- min(as.numeric(nodeData(g.lo, attr="xPos")))
   max.xPos <- max(as.numeric(nodeData(g.lo, attr="xPos")))
   checkEquals(abs(max.xPos - min.xPos), 1500)

   if(display){
     setGraph(tv, g.lo, names(models))
     setStyle(tv, "style.js")
     Sys.sleep(3); fit(tv)
     }

} # test_buildMultiModelGraph_twoModels_two_promoterSpans
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
