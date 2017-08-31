library (httpuv)
library (methods)
#----------------------------------------------------------------------------------------------------
trenaVizBrowserFile <- system.file(package="trenaViz", "browserCode", "dist", "trenaviz.html")
#----------------------------------------------------------------------------------------------------
.trenaViz <- setClass ("trenaViz",
                            representation = representation (),
                            contains = "BrowserVizClass",
                            prototype = prototype (uri="http://localhost", 9000)
                            )

#----------------------------------------------------------------------------------------------------
setGeneric ('ping',         signature='obj', function (obj) standardGeneric ('ping'))
setGeneric ('raiseTab',     signature='obj', function (obj, tabTitle) standardGeneric ('raiseTab'))
setGeneric ('getSelection', signature='obj', function (obj) standardGeneric ('getSelection'))
setGeneric ('setGenome',    signature='obj', function (obj, genomeName) standardGeneric ('setGenome'))
setGeneric ('setGraph',     signature='obj', function (obj, graph=NULL) standardGeneric ('setGraph'))

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
     printf("trenaViz::addGenome");
     payload <- genomeName
     send(obj, list(cmd="setGenome", callback="handleResponse", status="request", payload=payload))
     while (!browserResponseReady(obj)){
        Sys.sleep(.1)
        }
     getBrowserResponse(obj);
     })

#----------------------------------------------------------------------------------------------------
setMethod('setGraph', 'trenaViz',

  function (obj, graph=NULL) {
     printf("trenaViz::setGraph");
     if(is.null(graph))
       graph <- graphNEL()
     payload <- ""
     send(obj, list(cmd="setGraph", callback="handleResponse", status="request", payload=payload))
     while (!browserResponseReady(obj)){
        Sys.sleep(.1)
        }
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
     printf("TrenaViz::addBedTrackFromDataFrame");
     temp.filename <- "tmp.bed"
     write.table(tbl.bed, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE, file=temp.filename)
     payload <- list(name=trackName, bedFileName=temp.filename, displayMode=displayMode, color=color)
     send(obj, list(cmd="addBedTrackFromDataFrame", callback="handleResponse", status="request", payload=payload))
     while (!browserResponseReady(obj)){
        Sys.sleep(.1)
        }
     printf("browserResponseReady")
     getBrowserResponse(obj);
     })

#----------------------------------------------------------------------------------------------------
setMethod('addBedGraphTrackFromDataFrame', 'trenaViz',

  function (obj, trackName, tbl.bed, displayMode="COLLAPSED", minValue=NA, maxValue=NA, color="lightgray") {
     printf("TrenaViz::addBedGraphTrackFromDataFrame, color: %s", color);
     required.colnames <- c("chr", "start", "end", "score")
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
     printf("browserResponseReady")
     getBrowserResponse(obj);
     })

#----------------------------------------------------------------------------------------------------
setMethod('addBedTrackFromHostedFile', 'trenaViz',

  function (obj, trackName, uri, index.uri, displayMode="COLLAPSED", color) {
     printf("TrenaViz::addBedTrackFromHostedFile");
     payload <- list(name=trackName, uri=uri, indexUri=index.uri, displayMode=displayMode, color=color)
     send(obj, list(cmd="addBedTrackFromHostedFile", callback="handleResponse",
                    status="request", payload=payload))
     while (!browserResponseReady(obj)){
        Sys.sleep(.1)
        }
     printf("browserResponseReady")
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
