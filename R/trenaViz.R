#----------------------------------------------------------------------------------------------------
trenaVizBrowserFile <- system.file(package="trenaViz", "browserCode", "dist", "trenaviz.html")
#----------------------------------------------------------------------------------------------------
.trenaViz <- setClass ("trenaViz",
                            representation = representation (),
                            contains = "BrowserVizClass",
                            prototype = prototype (uri="http://localhost", 9000)
                            )

#----------------------------------------------------------------------------------------------------
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
                       function(obj, trackName, tbl.bed, displayMode="COLLAPSED", color="lightgray", trackHeight=200)
                   standardGeneric('addBedTrackFromDataFrame'))
setGeneric('addBedTrackFromHostedFile',   signature='obj',
                      function(obj, trackName, uri, index.uri=NA, displayMode="COLLAPSED", color="lightgray")
                   standardGeneric('addBedTrackFromHostedFile'))
setGeneric('addBedGraphTrackFromDataFrame', signature='obj',
           function(obj, trackName, tbl.bed, displayMode="COLLAPSED",minValue=NA, maxValue=NA,
                    color="lightgray", trackHeight=200)
                      standardGeneric('addBedGraphTrackFromDataFrame'))

setGeneric('selectNodes',         signature='obj', function(obj, nodeIDs) standardGeneric('selectNodes'))
setGeneric('getSelectedNodes',    signature='obj', function(obj) standardGeneric('getSelectedNodes'))

setGeneric('fit',                 signature='obj', function(obj, margin=50) standardGeneric('fit'))
setGeneric('fitSelected',         signature='obj', function(obj, margin=50) standardGeneric('fitSelected'))

#----------------------------------------------------------------------------------------------------
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
     g.json.as.javascript <- paste("network = ", g.json, sep="")

     if(!obj@quiet)printf("--- g.json conversion complete");
     filename <- file.path(getwd(), "g.json")
     payload <- list(filename=filename, modelNames=modelNames)
     if(!obj@quiet){
        printf("--- about to write file 'g.json' with %d characters", nchar(g.json))
        printf("--- first few characters: %s", substr(g.json, 1, 20))
        }
     write(g.json.as.javascript, file=filename)
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

  function (obj, trackName, tbl.bed, displayMode="COLLAPSED", color, trackHeight=100) {
     if(!obj@quiet)
        printf("TrenaViz::addBedTrackFromDataFrame");
     stopifnot(displayMode %in% c("COLLAPSED", "SQUISHED", "EXPANDED"))
     temp.filename <- "tmp.bed"
     if(!obj@quiet)
        printf("trenaViz.R about to write temporary bed file to %s", temp.filename);
     write.table(tbl.bed, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE, file=temp.filename)
     payload <- list(name=trackName, bedFileName=temp.filename, displayMode=displayMode, color=color,
                     trackHeight=trackHeight)
     send(obj, list(cmd="addBedTrackFromDataFrame", callback="handleResponse", status="request", payload=payload))
     while (!browserResponseReady(obj)){
        Sys.sleep(.1)
        }
     getBrowserResponse(obj);
     })

#----------------------------------------------------------------------------------------------------
setMethod('addBedGraphTrackFromDataFrame', 'trenaViz',

  function (obj, trackName, tbl.bed, displayMode="COLLAPSED", minValue=NA, maxValue=NA, color="lightgray", trackHeight=100) {
     if(!obj@quiet)
        printf("TrenaViz::addBedGraphTrackFromDataFrame, color: %s", color);
     stopifnot(displayMode %in% c("COLLAPSED", "SQUISHED", "EXPANDED"))
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
                     max=maxValue,
                     trackHeight=trackHeight)

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
     stopifnot(displayMode %in% c("COLLAPSED", "SQUISHED", "EXPANDED"))

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

       # allocate more character vectors that we could ever need; unused are deleted at conclusion

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
