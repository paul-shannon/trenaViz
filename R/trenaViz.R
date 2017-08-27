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
setGeneric ('ping',  signature='obj', function (obj) standardGeneric ('ping'))
setGeneric ('getSelection',  signature='obj', function (obj) standardGeneric ('getSelection'))
setGeneric ('setGenome', signature='obj', function (obj, genomeName) standardGeneric ('setGenome'))
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

   obj <- .trenaViz(BrowserViz(portRange, host, title, quiet, browserFile=trenaVizBrowserFile))
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
setMethod('setGenome', 'trenaViz',

  function (obj, genomeName) {
     printf("trenaViz::addGenome");
     payload <- genomeName
     send(obj, list(cmd="setGenome", callback="handleResponse", status="request", payload=payload))
     while (!browserResponseReady(obj)){
        Sys.sleep(.1)
        }
     printf("browserResponseReady")
     getBrowserResponse(obj);
     })

#----------------------------------------------------------------------------------------------------
