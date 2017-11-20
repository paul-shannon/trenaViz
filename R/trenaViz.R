#----------------------------------------------------------------------------------------------------
trenaVizBrowserFile <- system.file(package="trenaViz", "browserCode", "dist", "trenaviz.html")
#----------------------------------------------------------------------------------------------------
.trenaViz <- setClass ("trenaViz",
                            representation = representation (),
                            contains = "BrowserVizClass",
                            prototype = prototype (uri="http://localhost", 9000)
                            )

#----------------------------------------------------------------------------------------------------
setGeneric('buildMultiModelGraph', signature='obj', function(obj, targetGene, models)
                   standardGeneric('buildMultiModelGraph'))
setGeneric('addGeneModelLayout', signature='obj', function(obj, g, xPos.span=1500) standardGeneric('addGeneModelLayout'))
setGeneric('ping',         signature='obj', function (obj) standardGeneric ('ping'))
setGeneric('raiseTab',     signature='obj', function (obj, tabTitle) standardGeneric ('raiseTab'))
setGeneric('getSelection', signature='obj', function (obj) standardGeneric ('getSelection'))
setGeneric('setGenome',    signature='obj', function (obj, genomeName) standardGeneric ('setGenome'))
setGeneric('setGraph',     signature='obj', function (obj, graph=NULL, modelNames=NA) standardGeneric ('setGraph'))
setGeneric('setStyle',     signature='obj', function(obj, filename) standardGeneric ('setStyle'))

setGeneric('showGenomicRegion',   signature='obj', function(obj, regionString) standardGeneric('showGenomicRegion'))
setGeneric('getGenomicRegion',    signature='obj', function(obj) standardGeneric('getGenomicRegion'))

setGeneric('getTrackNames',        signature='obj', function(obj) standardGeneric('getTrackNames'))
setGeneric('removeTracksByName',   signature='obj', function(obj, trackNames) standardGeneric('removeTracksByName'))

setGeneric('addBedTrackFromDataFrame',  signature='obj',
                       function(obj, trackName, tbl.bed, displayMode="COLLAPSED", color="lightgray")
                   standardGeneric('addBedTrackFromDataFrame'))
setGeneric('addBedTrackFromHostedFile',   signature='obj',
                      function(obj, trackName, uri, index.uri=NA, displayMode="COLLAPSED", color="lightgray")
                   standardGeneric('addBedTrackFromHostedFile'))
setGeneric('addBedGraphTrackFromDataFrame', signature='obj',
                      function(obj, trackName, tbl.bed, displayMode="COLLAPSED",minValue=NA, maxValue=NA, color)
                      standardGeneric('addBedGraphTrackFromDataFrame'))

setGeneric('selectNodes',         signature='obj', function(obj, nodeIDs) standardGeneric('selectNodes'))
setGeneric('getSelectedNodes',    signature='obj', function(obj) standardGeneric('getSelectedNodes'))

setGeneric('fit',                 signature='obj', function(obj, margin=50) standardGeneric('fit'))
setGeneric('fitSelected',         signature='obj', function(obj, margin=50) standardGeneric('fitSelected'))

#----------------------------------------------------------------------------------------------------
setMethod('buildMultiModelGraph', 'trenaViz',

  function (obj, targetGene, models){

    stopifnot(is.list(models))
    stopifnot(is.character(names(models)))

      # we require all gene models to have the same scores
      # eg, they all have randomForest and lasso
      # detect disagreement across models, stop if found
    model.colnames <- colnames(models[[1]]$model)

    for(model in models){
       stopifnot(sort(names(model)) == c("model", "regions"))
       stopifnot(is.data.frame(model$model))
       stopifnot(nrow(model$model) >= 2);  # at least two rows
       stopifnot(is.data.frame(model$regions))  # regulatory regions
       stopifnot("gene" %in% colnames(model$model))
       stopifnot(all(model.colnames %in% colnames(model$model)))
       stopifnot(ncol(model$model) >= 2)  # at least "gene" and some score (usually multiple scores)
       missing.essential.colnames <- setdiff(c("motifName", "id", "distance.from.tss", "geneSymbol"),
                                             colnames(model$regions))
       if(length(missing.essential.colnames) > 0)
          stop(sprintf("essential colnames missing in regulatory regions (motif) table: %s",
                       paste(missing.essential.colnames, collapse=",")))
       tfs.in.model <- model$model$gene
       tfs.in.regions <- unique(model$regions$geneSymbol)
       stopifnot(sort(tfs.in.model) == sort(tfs.in.regions))
       } # for model


    g <- graphNEL(edgemode = "directed")
    model.names <- names(models)

    required.node.attribute.specs <- list(type="undefined",
                                          label="default node label",
                                          distance=0,
                                          #pearson=0,
                                          #randomForest=0,
                                          #pcaMax=0,
                                          #concordance=0,
                                          #betaLasso=0,
                                          motif="",
                                          motifInModel=TRUE,
                                          xPos=0,
                                          yPos=0)

       # remove "gene" from the colnames, leaving only the names of the scores we have been given
    score.names <- model.colnames[-match("gene", model.colnames)]
    optional.node.attribute.specs <- lapply(score.names, function(score) return(0))
    names(optional.node.attribute.specs) <- score.names
    node.attribute.specs <- c(required.node.attribute.specs, optional.node.attribute.specs)

    edge.attribute.spec <- list(edgeType="undefined")
    attribute.classes <- c("", model.names)  # "" (no prefix) is the currently displayed set of attibutes

      # create current version of these attributes, and then
      # per-model versions, which get mapped to current
      # in response to user's interactive choice on the cyjs user interface
      # the "current version" is, e.g., "distance".
      # per-model ("wt" and "mut" versions) become "wt.distance" and "mut.distance"
      # and are used by copying e.g. all wt.xxx attributes into the current (non-prefixed)
      # attribute, upon which the cyjs style is defined

    for(class.name in attribute.classes){
       class.name.prefix <- class.name  # with possible "." appended, permits standard and model-specific attributes
       if(nchar(class.name) > 0)
          class.name.prefix <- sprintf("%s.", class.name)
       noa.names.without.prefix <- names(node.attribute.specs)
       noa.names <- sprintf("%s%s", class.name.prefix, noa.names.without.prefix)
       noa.count <- length(node.attribute.specs)
       for(i in 1:noa.count){
          #printf("adding nodeDataDefaults: %s", noa.names[i])
          nodeDataDefaults(g, attr=noa.names[i]) <- node.attribute.specs[[noa.names.without.prefix[i]]]
          }
       } # for class

    #browser()
    edgeDataDefaults(g, attr = "edgeType") <- "undefined"

     #--------------------------------------------------------------------------------
     # 3 kinds of nodes:  1 targetGene, multiple tfs (each a geneSymbol from the
     # model), regulatory regions (binding sites, pfms matched to DNA)
     #--------------------------------------------------------------------------------

    tfs <- c()
    regulatoryRegions <- c()

      # collect all the tf and regulatory region nodes

    for(model in models){
       tbl.model <- model$model
       tfs <- unique(c(tfs, tbl.model$gene))
       tbl.reg <- model$regions
       regulatoryRegions <- unique(c(regulatoryRegions, tbl.reg$id))
       } # for model

    #printf("total tfs: %d   total regulatoryRegions: %d", length(tfs), length(regulatoryRegions))

    all.nodes <- unique(c(targetGene, tfs, regulatoryRegions))
    g <- addNode(all.nodes, g)

    #printf("--- browing after addNode in trenaViz::buildMultiModelGraph")
    #browser()
    nodeData(g, targetGene, "type") <- "targetGene"
    nodeData(g, tfs, "type")         <- "TF"
    nodeData(g, regulatoryRegions, "type")  <- "regulatoryRegion"
    nodeData(g, all.nodes, "label")  <- all.nodes

      # add edges, edge attribute, and the constant attributes for all of the regulatoryRegion nodes

    for(model in models){
       tfs <- model$regions$geneSymbol
       regRegions <- model$regions$id
       suppressWarnings(g <- addEdge(tfs, regRegions, g))
       edgeData(g,  tfs, regRegions, "edgeType") <- "bindsTo"
       suppressWarnings(g <- addEdge(regRegions, targetGene, g))
       edgeData(g, regRegions, targetGene, "edgeType") <- "regulatorySiteFor"
       tokensList <- strsplit(tbl.reg$id, "-")
       motif.labels <- unlist(lapply(tokensList, function(tokens) tokens[length(tokens)]))
       nodeData(g, tbl.reg$id, "label") <- motif.labels
       nodeData(g, tbl.reg$id, "distance") <- tbl.reg$distance.from.tss
       nodeData(g, tbl.reg$id, "motif") <- tbl.reg$motifName
       } # for model

      # now copy in the first model's tf node data

    #tbl.model <- models[[1]]$model

      # now copy in each of the model's tf and regRegion node data in turn
    #browser()
    model.names <- names(models)
    for(model.name in model.names){
       tbl.model <- models[[model.name]]$model
       for(optional.noa.name in names(optional.node.attribute.specs)){
          noa.name <- sprintf("%s.%s", model.name, optional.noa.name)
          nodeData(g, tbl.model$gene, attr=noa.name) <- tbl.model[, optional.noa.name]
          }
       tbl.regRegions <- models[[model.name]]$regions
       regRegionsInThisModel <- unique(tbl.regRegions$id)
       regRegionsNotInThisModel <- setdiff(regulatoryRegions, regRegionsInThisModel)
       attributeName <- sprintf("%s.%s", model.name, "motifInModel")
       nodeData(g, regRegionsInThisModel, attr=attributeName) <- TRUE
       nodeData(g, regRegionsNotInThisModel, attr=attributeName) <- FALSE
       } # for model.name

    g

    }) # buildMultiModelGraph
#----------------------------------------------------------------------------------------------------
setMethod('addGeneModelLayout', 'trenaViz',

  function (obj, g, xPos.span=1500){
    #printf("--- addGeneModelLayout")
    all.distances <- sort(unique(unlist(nodeData(g, attr='distance'), use.names=FALSE)))
    #print(all.distances)

    fp.nodes <- nodes(g)[which(unlist(nodeData(g, attr="type"), use.names=FALSE) == "regulatoryRegion")]
    tf.nodes <- nodes(g)[which(unlist(nodeData(g, attr="type"), use.names=FALSE) == "TF")]
    targetGene.nodes <- nodes(g)[which(unlist(nodeData(g, attr="type"), use.names=FALSE) == "targetGene")]

     # add in a zero in case all of the footprints are up or downstream of the 0 coordinate, the TSS
    span.endpoints <- range(c(0, as.numeric(nodeData(g, fp.nodes, attr="distance"))))
    span <- max(span.endpoints) - min(span.endpoints)
    footprintLayoutFactor <- 1
    #printf("initial:  span: %d  footprintLayoutFactor: %f", span, footprintLayoutFactor)

    footprintLayoutFactor <- xPos.span/span

    #if(span < 600)  #
    #   footprintLayoutFactor <- 600/span
    #if(span > 1000)
    #   footprintLayoutFactor <- span/1000

    #printf("corrected:  span: %d  footprintLayoutFactor: %f", span, footprintLayoutFactor)

    xPos <- as.numeric(nodeData(g, fp.nodes, attr="distance")) * footprintLayoutFactor
    yPos <- 0
    nodeData(g, fp.nodes, "xPos") <- xPos
    nodeData(g, fp.nodes, "yPos") <- yPos

    adjusted.span.endpoints <- range(c(0, as.numeric(nodeData(g, fp.nodes, attr="xPos"))))
    #printf("raw span of footprints: %d   footprintLayoutFactor: %f  new span: %8.0f",
    #       span, footprintLayoutFactor, abs(max(adjusted.span.endpoints) - min(adjusted.span.endpoints)))

    tfs <- names(which(nodeData(g, attr="type") == "TF"))

    for(tf in tfs){
       footprint.neighbors <- edges(g)[[tf]]
       if(length(footprint.neighbors) > 0){
          footprint.positions <- as.integer(nodeData(g, footprint.neighbors, attr="xPos"))
          new.xPos <- mean(footprint.positions)
          #if(is.na(new.xPos)) browser()
          #if(is.nan(new.xPos)) browser()
          #printf("%8s: %5d", tf, new.xPos)
          }
       else{
          new.xPos <- 0
          }
       nodeData(g, tf, "yPos") <- sample(300:1200, 1)
       nodeData(g, tf, "xPos") <- new.xPos
       } # for tf

    nodeData(g, targetGene.nodes, "xPos") <- 0
    nodeData(g, targetGene.nodes, "yPos") <- -200

    g

    }) # addGeneModelLayout

#------------------------------------------------------------------------------------------------------------------------
setupMessageHandlers <- function()
{
   addRMessageHandler("handleResponse", "handleResponse")

} # setupMessageHandlers
#----------------------------------------------------------------------------------------------------
# constructor
trenaViz = function(portRange, host="localhost", title="trenaViz", quiet=TRUE)
{
   if(!quiet){
      printf("want to load %s", trenaVizBrowserFile)
      }

   obj <- .trenaViz(BrowserViz(portRange, host, title, quiet, browserFile=trenaVizBrowserFile,
                               httpQueryProcessingFunction=myQP))
   setBrowserWindowTitle(obj, title)

   obj

} # trenaViz: constructor
#----------------------------------------------------------------------------------------------------
setMethod('ping', 'trenaViz',

  function (obj) {
     send(obj, list(cmd="ping", callback="handleResponse", status="request", payload=""))
     while (!browserResponseReady(obj)){
        if(!obj@quiet) message(sprintf("plot waiting for browser response"));
        Sys.sleep(.1)
        }
     getBrowserResponse(obj)
     }) # ping

#----------------------------------------------------------------------------------------------------
setMethod('raiseTab', 'trenaViz',

  function (obj, tabTitle) {
     send(obj, list(cmd="raiseTab", callback="handleResponse", status="request", payload=tabTitle))
     while (!browserResponseReady(obj)){
        if(!obj@quiet) message(sprintf("plot waiting for browser response"));
        Sys.sleep(.1)
        }
     getBrowserResponse(obj)
     }) # raiseTab

#----------------------------------------------------------------------------------------------------
setMethod('setGenome', 'trenaViz',

  function (obj, genomeName) {
     if(!obj@quiet) printf("trenaViz::addGenome");
     payload <- genomeName
     send(obj, list(cmd="setGenome", callback="handleResponse", status="request", payload=payload))
     while (!browserResponseReady(obj)){
        Sys.sleep(.1)
        }
     getBrowserResponse(obj);
     })

#----------------------------------------------------------------------------------------------------
setMethod('setGraph', 'trenaViz',

  function (obj, graph=NULL, modelNames=NA) {
     if(!obj@quiet){
        printf("trenaViz::setGraph");
        print(graph)
        printf("--- converting graph to JSON");
        }
     g.json <- .graphToJSON(graph)
     if(!obj@quiet)printf("--- g.json conversion complete");
     filename <- "g.json"
     payload <- list(filename=filename, modelNames=modelNames)
     if(!obj@quiet){
        printf("--- about to write file 'g.json' with %d characters", nchar(g.json))
        printf("--- first few characters: %s", substr(g.json, 1, 20))
        }
     write(g.json, file=filename)
     if(!obj@quiet)
         printf("--- file writing complete")
     send(obj, list(cmd="setGraph", callback="handleResponse", status="request", payload=payload))
     while (!browserResponseReady(obj)){
        Sys.sleep(.1)
        }
     getBrowserResponse(obj);
     })

#----------------------------------------------------------------------------------------------------
setMethod('setStyle', 'trenaViz',

  function (obj, filename) {
     send(obj, list(cmd="setStyle", callback="handleResponse", status="request", payload=filename))
     while (!browserResponseReady(obj)){
        Sys.sleep(.1)
        }
     if(!obj@quiet)
        printf("browserResponseReady")
     getBrowserResponse(obj);
     })

#----------------------------------------------------------------------------------------------------
setMethod('showGenomicRegion', 'trenaViz',

   function (obj, regionString) {
     payload <- list(regionString=regionString)
     send(obj, list(cmd="showGenomicRegion", callback="handleResponse", status="request", payload=payload))
     while (!browserResponseReady(obj)){
        Sys.sleep(.1)
        }
     getBrowserResponse(obj);
     })

#----------------------------------------------------------------------------------------------------
setMethod('getGenomicRegion', 'trenaViz',

   function (obj) {
     payload <- ""
     send(obj, list(cmd="getGenomicRegion", callback="handleResponse", status="request", payload=payload))
     while (!browserResponseReady(obj)){
        Sys.sleep(.1)
        }
     getBrowserResponse(obj);
     })

#----------------------------------------------------------------------------------------------------
setMethod('getTrackNames', 'trenaViz',

   function (obj) {
     payload <- ""
     send(obj, list(cmd="getTrackNames", callback="handleResponse", status="request", payload=payload))
     while (!browserResponseReady(obj)){
        Sys.sleep(.1)
        }
     getBrowserResponse(obj);
     })

#----------------------------------------------------------------------------------------------------
setMethod('removeTracksByName', 'trenaViz',

   function (obj, trackNames) {
     payload <- trackNames
     send(obj, list(cmd="removeTracksByName", callback="handleResponse", status="request", payload=payload))
     while (!browserResponseReady(obj)){
        Sys.sleep(.1)
        }
     getBrowserResponse(obj);
     })

#----------------------------------------------------------------------------------------------------
setMethod('addBedTrackFromDataFrame', 'trenaViz',

  function (obj, trackName, tbl.bed, displayMode="COLLAPSED", color) {
     if(!obj@quiet)
        printf("TrenaViz::addBedTrackFromDataFrame");
         # work around the httpuv/igv.js > 1.0.9 bug, in which somehow the last row is dropped.
         # simply duplicate that last row
     tbl.bed <- rbind(tbl.bed, tbl.bed[nrow(tbl.bed),])
     temp.filename <- "tmp.bed"
     write.table(tbl.bed, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE, file=temp.filename)
     payload <- list(name=trackName, bedFileName=temp.filename, displayMode=displayMode, color=color)
     send(obj, list(cmd="addBedTrackFromDataFrame", callback="handleResponse", status="request", payload=payload))
     while (!browserResponseReady(obj)){
        Sys.sleep(.1)
        }
     getBrowserResponse(obj);
     })

#----------------------------------------------------------------------------------------------------
setMethod('addBedGraphTrackFromDataFrame', 'trenaViz',

  function (obj, trackName, tbl.bed, displayMode="COLLAPSED", minValue=NA, maxValue=NA, color="lightgray") {
     if(!obj@quiet)
        printf("TrenaViz::addBedGraphTrackFromDataFrame, color: %s", color);
     found.chromosome.column <- any(grepl("^chr", colnames(tbl.bed), ignore.case=TRUE))
     stopifnot(found.chromosome.column)
     required.colnames <- c("start", "end", "score")
     missing.colnames <- setdiff(required.colnames, colnames(tbl.bed))
     if(length(missing.colnames) > 0){
        printf("addBedGraphTrackFromDataFrame detects missing column name: %s",
               paste(missing.colnames, collapse=", "))
        return()
        }

     if(is.na(minValue))
        minValue <- min(tbl.bed$score)

     if(is.na(maxValue))
        maxValue <- max(tbl.bed$score)

     temp.filename <- "tmp.bedGraph"
     write.table(tbl.bed, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE, file=temp.filename)
     payload <- list(name=trackName,
                     bedFileName=temp.filename,
                     displayMode=displayMode,
                     color=color,
                     min=minValue,
                     max=maxValue)

     send(obj, list(cmd="addBedGraphTrackFromDataFrame", callback="handleResponse", status="request", payload=payload))
     while (!browserResponseReady(obj)){
        Sys.sleep(.1)
        }
     #printf("browserResponseReady")
     getBrowserResponse(obj);
     })

#----------------------------------------------------------------------------------------------------
setMethod('addBedTrackFromHostedFile', 'trenaViz',

  function (obj, trackName, uri, index.uri, displayMode="COLLAPSED", color) {
     if(!obj@quiet)
        printf("TrenaViz::addBedTrackFromHostedFile");
     payload <- list(name=trackName, uri=uri, indexUri=index.uri, displayMode=displayMode, color=color)
     send(obj, list(cmd="addBedTrackFromHostedFile", callback="handleResponse",
                    status="request", payload=payload))
     while (!browserResponseReady(obj)){
        Sys.sleep(.1)
        }
     getBrowserResponse(obj);
     })

#----------------------------------------------------------------------------------------------------
setMethod('selectNodes', 'trenaViz',

  function (obj, nodeIDs) {
     payload <- list(nodeIDs=nodeIDs)
     send(obj, list(cmd="selectNodes", callback="handleResponse", status="request", payload=payload))
     while (!browserResponseReady(obj)){
        Sys.sleep(.1)
        }
     getBrowserResponse(obj);
     })

#----------------------------------------------------------------------------------------------------
setMethod('getSelectedNodes', 'trenaViz',

  function (obj) {
     payload <- ""
     send(obj, list(cmd="getSelectedNodes", callback="handleResponse", status="request", payload=payload))
     while (!browserResponseReady(obj)){
        Sys.sleep(.1)
        }
     result <- fromJSON(getBrowserResponse(obj))$id;
     if(all(is.null(result)))
        return(list())
     else
        return(result)
     })

#----------------------------------------------------------------------------------------------------
setMethod('fit', 'trenaViz',

  function (obj, margin=50) {
     send(obj, list(cmd="fit", callback="handleResponse", status="request", payload=margin))
     while (!browserResponseReady(obj)){
        Sys.sleep(.1)
        }
     getBrowserResponse(obj);
     })

#----------------------------------------------------------------------------------------------------
setMethod('fitSelected', 'trenaViz',

  function (obj, margin=50) {
     send(obj, list(cmd="fitSelected", callback="handleResponse", status="request", payload=margin))
     while (!browserResponseReady(obj)){
        Sys.sleep(.1)
        }
     getBrowserResponse(obj);
     })

#----------------------------------------------------------------------------------------------------
# the basic form
#
#   network =  {"nodes":[{"data":{"id":"o"},"position":{"x":779.5, "y":500 }}],
#               "edges":[{"data":{"source":"o","target":"o","id":"e1"}}]};
#
.graphToJSON <- function(g)
{
   if(length(nodes(g)) == 0)
      return ("{}")

       # allocate more character vectors that we could ever need
    vector.count <- 10 * (length(edgeNames(g)) + length (nodes(g)))
    vec <- vector(mode="character", length=vector.count)
    i <- 1;

    vec[i] <- '{"elements": {"nodes": ['; i <- i + 1;
    nodes <- nodes(g)
    edgeNames <- edgeNames(g)
    edges <- strsplit(edgeNames, "~")  # a list of pairs
    edgeNames <- sub("~", "->", edgeNames)
    names(edges) <- edgeNames

    noa.names <- names(nodeDataDefaults(g))
    eda.names <- names(edgeDataDefaults(g))
    nodeCount <- length(nodes)
    edgeCount <- length(edgeNames)

    for(n in 1:nodeCount){
       node <- nodes[n]
       vec[i] <- '{"data": '; i <- i + 1
       nodeList <- list(id = node)
       this.nodes.data <- nodeData(g, node)[[1]]
       if(length(this.nodes.data) > 0)
          nodeList <- c(nodeList, this.nodes.data)
       nodeList.json <- toJSON(nodeList, auto_unbox=TRUE)
       vec[i] <- nodeList.json; i <- i + 1
       if(all(c("xPos", "yPos") %in% names(nodeDataDefaults(g)))){
          position.markup <- sprintf(', "position": {"x": %f, "y": %f}',
                                     nodeData(g, node, "xPos")[[1]],
                                     nodeData(g, node, "yPos")[[1]])
          vec[i] <- position.markup
          i <- i + 1
          }
        if(n != nodeCount){
           vec [i] <- "},"; i <- i + 1 # sprintf("%s},", x)  # another node coming, add a comma
           }
       } # for n

    vec [i] <- "}]"; i <- i + 1  # close off the last node, the node array ], the nodes element }

    if(edgeCount > 0){
       vec[i] <- ', "edges": [' ; i <- i + 1
       for(e in seq_len(edgeCount)) {
          vec[i] <- '{"data": '; i <- i + 1
          edgeName <- edgeNames[e]
          edge <- edges[[e]]
          sourceNode <- edge[[1]]
          targetNode <- edge[[2]]
          edgeList <- list(id=edgeName, source=sourceNode, target=targetNode)
          this.edges.data <- edgeData(g, sourceNode, targetNode)[[1]]
          if(length(this.edges.data) > 0)
             edgeList <- c(edgeList, this.edges.data)
          edgeList.json <- toJSON(edgeList, auto_unbox=TRUE)
          vec[i] <- edgeList.json; i <- i + 1
          if(e != edgeCount){          # add a comma, ready for the next edge element
             vec [i] <- '},'; i <- i + 1
             }
          } # for e
      vec [i] <- "}]"; i <- i + 1
      } # if edgeCount > 0

   vec [i] <- "}"  # close the edges object
   i <- i + 1;
   vec [i] <- "}"  # close the elements object
   vec.trimmed <- vec [which(vec != "")]
   printf("%d strings used in constructing json", length(vec.trimmed))
   paste0(vec.trimmed, collapse=" ")

} # .graphToJSON
#------------------------------------------------------------------------------------------------------------------------
# {elements: [
#    {data: {id: 'a', score:5}, position: {x: 100, y: 200}},
#    {data: {id: 'b', score:100}, position: {x: 200, y: 200}},
#    {data: {id: 'e1', source: 'a', target: 'b'}}
#    ],  // elements array
# layout: { name: 'preset'},
# style: [{selector: 'node', style: {'content': 'data(id)'}}]
# }
old.graphToJSON <- function(g)
{
    x <- '{"elements": [';
    nodes <- nodes(g)
    edgeNames <- edgeNames(g)
    edges <- strsplit(edgeNames, "~")  # a list of pairs
    edgeNames <- sub("~", "->", edgeNames)
    names(edges) <- edgeNames

    noa.names <- names(nodeDataDefaults(g))
    eda.names <- names(edgeDataDefaults(g))
    nodeCount <- length(nodes)
    edgeCount <- length(edgeNames)

    for(n in 1:nodeCount){
       node <- nodes[n]
       x <- sprintf('%s {"data": {"id": "%s"', x, node);
       nodeAttributeCount <- length(noa.names)
       for(i in seq_len(nodeAttributeCount)){
          noa.name <- noa.names[i];
          value <-  nodeData(g, node, noa.name)[[1]]
          if(is.numeric(value))
             x <- sprintf('%s, "%s": %s', x, noa.name, value)
          else
             x <- sprintf('%s, "%s": "%s"', x, noa.name, value)
          #printf("2: %d", nchar(x))
          } # for i
       #printf("3: %d", nchar(x))
       x <- sprintf('%s}', x)     # close off this node data element
       #printf("4: %d", nchar(x))
       if(all(c("xPos", "yPos") %in% noa.names)){
           xPos <- as.integer(nodeData(g, node, "xPos"))
           yPos <- as.integer(nodeData(g, node, "yPos"))
           x <- sprintf('%s, "position": {"x": %d, "y": %d}', x, xPos, yPos)
           } # add position element
       #printf("5: %d", nchar(x))
       x <- sprintf('%s}', x)     # close off this node data element
       #printf("6: %d", nchar(x))
       if(n != nodeCount)
           x <- sprintf("%s,", x)  # another node coming, add a comma
       #printf("7: %d", nchar(x))
       #browser()
       xyz <- 99
       } # for n

    #printf("--- browser before edge loop")
    #browser()
    for(e in seq_len(edgeCount)) {
       edgeName <- edgeNames[e]
       edge <- edges[[e]]
       sourceNode <- edge[[1]]
       targetNode <- edge[[2]]
       x <- sprintf('%s, {"data": {"id": "%s", "source": "%s", "target": "%s"', x, edgeName, sourceNode, targetNode);
       edgeAttributeCount <- length(eda.names)
       for(i in seq_len(edgeAttributeCount)){
          eda.name <- eda.names[i];
          value <-  edgeData(g, sourceNode, targetNode, eda.name)[[1]]
          if(is.numeric(value))
             x <- sprintf('%s, "%s": %s', x, eda.name, value)
          else
             x <- sprintf('%s, "%s": "%s"', x, eda.name, value)
          } # for each edgeAttribute
       x <- sprintf('%s}}', x)     # close off this edge data element
       } # for e

    #printf("--- browser before closing")
    #browser()
    x <- sprintf("%s]}", x)

    x

} # old.graphToJSON
#------------------------------------------------------------------------------------------------------------------------
myQP <- function(queryString)
{
   #printf("=== TReNA-Viz::myQP");
   #print(queryString)
     # for reasons not quite clear, the query string comes in with extra characters
     # following the expected filename:
     #
     #  "?sampleStyle.js&_=1443650062946"
     #
     # check for that, cleanup the string, then see if the file can be found

   ampersand.loc <- as.integer(regexpr("&", queryString, fixed=TRUE))
   #printf("ampersand.loc: %d", ampersand.loc)

   if(ampersand.loc > 0){
      queryString <- substring(queryString, 1, ampersand.loc - 1);
      }

   questionMark.loc <- as.integer(regexpr("?", queryString, fixed=TRUE));
   #printf("questionMark.loc: %d", questionMark.loc)

   if(questionMark.loc == 1)
      queryString <- substring(queryString, 2, nchar(queryString))

   filename <- queryString;
   #printf("myQP filename: '%s'", filename)
   #printf("       exists?  %s", file.exists(filename));

   stopifnot(file.exists(filename))

   #printf("--- about to scan %s", filename);
      # reconstitute linefeeds though collapsing file into one string, so json
      # structure is intact, and any "//" comment tokens only affect one line
   text <- paste(scan(filename, what=character(0), sep="\n", quiet=TRUE), collapse="\n")
   #printf("%d chars read from %s", nchar(text), filename);

   return(text);

} # myQP
#----------------------------------------------------------------------------------------------------
