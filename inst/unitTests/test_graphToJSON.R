library(RUnit)
library(trenaViz)
#------------------------------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#------------------------------------------------------------------------------------------------------------------------
# use ~/github/projects/examples/cyjsMinimal/cyjs.html to test out json strings produced here
#------------------------------------------------------------------------------------------------------------------------
if(!exists("g.big")){
   load(system.file(package="trenaViz", "extdata",  "graph_1669nodes_3260edges_challenge_for_converting_to_json.RData"))
   g.big <- g.lo
   }

if(!exists("g.small")){
   print(load(system.file(package="trenaViz", "extdata",  "graph.11nodes.14edges.RData")))
   g.small <- g
   }


#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
} # runTests
#------------------------------------------------------------------------------------------------------------------------
createTestGraph <- function(nodeCount, edgeCount)
{
   elementCount <- nodeCount^2;
   vec <- rep(0, elementCount)

   set.seed(13);
   vec[sample(1:elementCount, edgeCount)] <- 1
   mtx <- matrix(vec, nrow=nodeCount)

   gam <- graphAM(adjMat=mtx, edgemode="directed")

   as(gam, "graphNEL")

} # createTestGraph
#----------------------------------------------------------------------------------------------------
test_1669_3260 <- function(display=FALSE)
{
   printf("--- test_1669_3260")
   Rprof()
   g.json <- pre.allocate.graphToJSON(g.small)

   if(display){
      writeLines(sprintf("network = %s", g.json), "network.js")
      browseURL("cyjs-readNetworkFromFile.html")
      } # display


   g2 <- fromJSON(g.json, flatten=TRUE)
   checkEquals(lapply(g2, dim), list(nodes=c(11, 25), edges=c(14,4)))

   g.json <- pre.allocate.graphToJSON(g.big)

   if(display){
      writeLines(sprintf("network = %s", g.json), "network.js")
      browseURL("cyjs-readNetworkFromFile.html")
      } # display

   g2 <- fromJSON(g.json, flatten=TRUE)
   checkEquals(lapply(g2, dim), list(nodes=c(11, 25), edges=c(14,4)))

    Rprof(NULL)
   summaryRprof()

 } # test_1669_3260
#------------------------------------------------------------------------------------------------------------------------
test_2_nodes_2_edges_no_attributes <- function(display=FALSE)
{
   printf("--- test_2_nodes_2_edges_no_attributes")

   g <- createTestGraph(2, 2)
   g.json <- graphToJSON(g)

   if(display){
      writeLines(sprintf("network = %s", g.json), "network.js")
      browseURL("cyjs-readNetworkFromFile.html")
      } # display

   g2 <- fromJSON(g.json, flatten=TRUE)
   tbl.nodes <- g2$nodes
   checkEquals(tbl.nodes$data.id, nodes(g))
   tbl.edges <- g2$edges
   checkEquals(dim(tbl.edges), c(2, 3))

 } # test_2_nodes_2_edges_no_attributes
#------------------------------------------------------------------------------------------------------------------------
test_20_nodes_20_edges_no_attributes <- function(display=FALSE)
{
   printf("--- test_20_nodes_20_edges_no_attributes")

   g <- createTestGraph(20, 20)
   g.json <- graphToJSON(g)

   if(display){
      writeLines(sprintf("network = %s", g.json), "network.js")
      browseURL("cyjs-readNetworkFromFile.html")
      } # display

   g2 <- fromJSON(g.json, flatten=TRUE)
   tbl.nodes <- g2$nodes
   checkEquals(tbl.nodes$data.id, nodes(g))
   tbl.edges <- g2$edges
   checkEquals(dim(tbl.edges), c(20, 3))

 } # test_2_nodes_2_edges_no_attributes
#------------------------------------------------------------------------------------------------------------------------
test_200_nodes_200_edges_no_attributes <- function(display=FALSE)
{
   printf("--- test_200_nodes_200_edges_no_attributes")

   g <- createTestGraph(200, 200)
   g.json <- graphToJSON(g)

   if(display){
      writeLines(sprintf("network = %s", g.json), "network.js")
      browseURL("cyjs-readNetworkFromFile.html")
      } # display

   g2 <- fromJSON(g.json, flatten=TRUE)
   tbl.nodes <- g2$nodes
   checkEquals(tbl.nodes$data.id, nodes(g))
   tbl.edges <- g2$edges
   checkEquals(dim(tbl.edges), c(200, 3))

 } # test_200_nodes_200_edges_no_attributes
#------------------------------------------------------------------------------------------------------------------------
test_2000_nodes_2000_edges_no_attributes <- function(display=FALSE)
{
   printf("--- test_2000_nodes_2000_edges_no_attributes")

   print(system.time({   # 4 seconds
      g <- createTestGraph(2000, 2000)
      g.json <- graphToJSON(g)
      }))

   if(display){
      writeLines(sprintf("network = %s", g.json), "network.js")
      browseURL("cyjs-readNetworkFromFile.html")
      } # display

   g2 <- fromJSON(g.json, flatten=TRUE)
   tbl.nodes <- g2$nodes
   checkEquals(tbl.nodes$data.id, nodes(g))
   tbl.edges <- g2$edges
   checkEquals(dim(tbl.edges), c(2000, 3))

 } # test_2000_nodes_2000_edges_no_attributes
#------------------------------------------------------------------------------------------------------------------------
test_1_node <- function(display=FALSE)
{
   printf("--- test_1_node")
   g <- graphNEL(nodes="A", edgemode="directed")
   g.json <- pre.allocate.graphToJSON(g)

   if(display){
      writeLines(sprintf("network = %s", g.json), "network.js")
      browseURL("cyjs-readNetworkFromFile.html")
      } # display

   g2 <- fromJSON(g.json, flatten=TRUE)
   tbl.nodes <- g2$nodes
   checkEquals(tbl.nodes$data.id, nodes(g))

} # test_1_node
#------------------------------------------------------------------------------------------------------------------------
test_2_nodes <- function(display=FALSE)
{
   printf("--- test_2_nodes")

   g <- graphNEL(nodes=c("A", "B"), edgemode="directed")
   g.json <- pre.allocate.graphToJSON(g)

   if(display){
      writeLines(sprintf("network = %s", g.json), "network.js")
      browseURL("cyjs-readNetworkFromFile.html")
      } # display

   g2 <- fromJSON(g.json, flatten=TRUE)
   tbl.nodes <- g2$nodes
   checkEquals(tbl.nodes$data.id, nodes(g))

} # test_2_nodes
#------------------------------------------------------------------------------------------------------------------------
test_2_nodes_1_edge <- function(display=FALSE)
{
   printf("--- test_2_nodes_1_edge")

   g <- graphNEL(nodes=c("X", "Y"), edgemode="directed")
   g <- addEdge("X", "Y", g);
   g.json <- pre.allocate.graphToJSON(g)

   if(display){
      writeLines(sprintf("network = %s", g.json), "network.js")
      browseURL("cyjs-readNetworkFromFile.html")
      } # display

      #  flatten: automatically ‘flatten’ nested data frames into a single non-nested data frame
   g2 <- fromJSON(g.json, flatten=TRUE)
   checkEquals(names(g2), c("nodes", "edges"))
   tbl.nodes <- g2$nodes
   checkEquals(dim(tbl.nodes), c(2,1))
   checkEquals(tbl.nodes$data.id, c("X", "Y"))

   tbl.edges <- g2$edges
   checkEquals(dim(tbl.edges), c(1,3))
   checkEquals(tbl.edges$data.id, "X->Y")

} # test_2_nodes_1_edge
#------------------------------------------------------------------------------------------------------------------------
test_1_node_2_attributes <- function(display=FALSE)
{
   printf("--- test_1_node_2_attributse")

   g <- graphNEL(nodes="A", edgemode="directed")
   nodeDataDefaults(g, "size") <- 0
   nodeData(g, "A", "size") <- 99

   nodeDataDefaults(g, "label") <- ""
   nodeData(g, "A", "label") <- "bigA"

   g.json <- pre.allocate.graphToJSON(g)

   if(display){
      writeLines(sprintf("network = %s", g.json), "network.js")
      browseURL("cyjs-readNetworkFromFile.html")
      } # display

   g2 <- fromJSON(g.json, flatten=TRUE)
   tbl.nodes <- g2$nodes
   checkEquals(tbl.nodes$data.id, nodes(g))
   checkEquals(tbl.nodes$data.size, 99)
   checkEquals(tbl.nodes$data.label, "bigA")

} # test_1_node_2_attributes
#------------------------------------------------------------------------------------------------------------------------
test_2_nodes_1_edge_2_edgeAttribute <- function(display=FALSE)
{
   printf("--- test_2_nodes_2_edgeAttributes")

   g <- graphNEL(nodes=c("X", "Y"), edgemode="directed")
   g <- addEdge("X", "Y", g);
   edgeDataDefaults(g, "weight") <- 0
   edgeDataDefaults(g, "edgeType") <- "generic"
   edgeData(g, "X", "Y", "weight") <- 1.234
   edgeData(g, "X", "Y", "edgeType") <- "regulates"

   g.json <- pre.allocate.graphToJSON(g)

   if(display){
      writeLines(sprintf("network = %s", g.json), "network.js")
      browseURL("cyjs-readNetworkFromFile.html")
      } # display

      #  flatten: automatically ‘flatten’ nested data frames into a single non-nested data frame
   g2 <- fromJSON(g.json, flatten=TRUE)
   checkEquals(names(g2), c("nodes", "edges"))
   tbl.nodes <- g2$nodes
   checkEquals(dim(tbl.nodes), c(2,1))
   checkEquals(tbl.nodes$data.id, c("X", "Y"))

   tbl.edges <- g2$edges
   checkEquals(dim(tbl.edges), c(1,5))
   checkEquals(tbl.edges$data.id, "X->Y")
   checkEquals(tbl.edges$data.source, "X")
   checkEquals(tbl.edges$data.target, "Y")
   checkEquals(tbl.edges$data.weight, 1.234)
   checkEquals(tbl.edges$data.edgeType, "regulates")

} # test_2_nodes_1_edge
#------------------------------------------------------------------------------------------------------------------------
test_smallGraphWithAttributes <- function(display=FALSE)
{
   printf("--- test_smallGraphWithAttributes")
   g <- simpleDemoGraph()
   g.json <- pre.allocate.graphToJSON(g)

   if(display){
      writeLines(sprintf("network = %s", g.json), "network.js")
      browseURL("cyjs-readNetworkFromFile.html")
      } # display

   g2 <- fromJSON(g.json, flatten=TRUE)
   checkEquals(names(g2), c("nodes", "edges"))
   tbl.nodes <- g2$nodes
   tbl.edges <- g2$edges

   checkEquals(dim(tbl.nodes), c(3, 5))
   checkEquals(colnames(tbl.nodes),
               c("data.id", "data.type", "data.lfc", "data.label", "data.count"))
   checkEquals(dim(tbl.edges), c(3, 6))
   checkEquals(colnames(tbl.edges), c("data.id", "data.source", "data.target", "data.edgeType", "data.score", "data.misc"))

} # test_smallGraphWithAttributes
#------------------------------------------------------------------------------------------------------------------------
graphToJSON <- function(g)
{
   if(length(nodes(g)) == 0)
      return ("{}")

    x <- '{"nodes": [';
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
       x <- sprintf('%s{"data": ', x)
       nodeList <- list(id = node)
       this.nodes.data <- nodeData(g, node)[[1]]
       if(length(this.nodes.data) > 0)
          nodeList <- c(nodeList, this.nodes.data)
       nodeList.json <- toJSON(nodeList, auto_unbox=TRUE)
       x <- sprintf("%s %s", x, nodeList.json)
        if(n != nodeCount)
           x <- sprintf("%s},", x)  # another node coming, add a comma
       } # for n

    x <- sprintf("%s}]", x)         # close off the last node, the node array ], the nodes element }

    if(edgeCount > 0){
       x <- sprintf('%s, "edges": [', x)
       for(e in seq_len(edgeCount)) {
          x <- sprintf('%s{"data": ', x)
          edgeName <- edgeNames[e]
          edge <- edges[[e]]
          sourceNode <- edge[[1]]
          targetNode <- edge[[2]]
          edgeList <- list(id=edgeName, source=sourceNode, target=targetNode)
          this.edges.data <- edgeData(g, sourceNode, targetNode)[[1]]
          if(length(this.edges.data) > 0)
             edgeList <- c(edgeList, this.edges.data)
          edgeList.json <- toJSON(edgeList, auto_unbox=TRUE)
          x <- sprintf("%s %s", x, edgeList.json)
          if(e != edgeCount)          # add a comma, ready for the next edge element
             x <- sprintf('%s},', x)
          } # for e
      x <- sprintf("%s}]", x)         # edge elements now complete, close the array
      } # if edgeCount > 0

   x <- sprintf("%s}", x)

   x

} # graphToJSON
#------------------------------------------------------------------------------------------------------------------------
# follow this simplified form, in which neither "node" nor "edge" is neded:
# cy.json({elements:[
#      {data: {id:'a'}},
#      {data: {id: 'e1', source: 'a', target: 'a'}}
#      ]})
#
# or
#
#   network =  {"nodes":[{"data":{"id":"o"},"position":{"x":779.5, "y":500 }}],
#               "edges":[{"data":{"source":"o","target":"o","id":"e1"},"position":{}}]};
#
pre.allocate.graphToJSON <- function(g)
{
   if(length(nodes(g)) == 0)
      return ("{}")

       # allocate more character vectors that we could ever need
    vector.count <- 10 * (length(edgeNames(g)) + length (nodes(g)))
    vec <- vector(mode="character", length=vector.count)
    i <- 1;

    vec[i] <- '{"nodes": ['; i <- i + 1;
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
       #x <- sprintf('%s{"data": ', x)
       vec[i] <- '{"data": '; i <- i + 1
       nodeList <- list(id = node)
       this.nodes.data <- nodeData(g, node)[[1]]
       if(length(this.nodes.data) > 0)
          nodeList <- c(nodeList, this.nodes.data)
       nodeList.json <- toJSON(nodeList, auto_unbox=TRUE)
       #x <- sprintf("%s %s", x, nodeList.json)
       vec[i] <- nodeList.json; i <- i + 1
        if(n != nodeCount){
           #x <- sprintf("%s},", x)  # another node coming, add a comma
           vec [i] <- "},"; i <- i + 1 # sprintf("%s},", x)  # another node coming, add a comma
           }
       } # for n

    #x <- sprintf("%s}]", x)         # close off the last node, the node array ], the nodes element }
    vec [i] <- "}]"; i <- i + 1 # , x)         # close off the last node, the node array ], the nodes element }

    if(edgeCount > 0){
       #x <- sprintf('%s, "edges": [', x)
       vec[i] <- ', "edges": [' ; i <- i + 1
       for(e in seq_len(edgeCount)) {
          #x <- sprintf('%s{"data": ', x)
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
          # x <- sprintf("%s %s", x, edgeList.json)
          vec[i] <- edgeList.json; i <- i + 1
          if(e != edgeCount){          # add a comma, ready for the next edge element
             #x <- sprintf('%s},', x)
             vec [i] <- '},'; i <- i + 1
             }
          } # for e
      #x <- sprintf("%s}]", x)         # edge elements now complete, close the array
      vec [i] <- "}]"; i <- i + 1
      } # if edgeCount > 0

   #x <- sprintf("%s}", x)
   vec [i] <- "}"
   vec.trimmed <- vec [which(vec != "")]
   printf("%d strings used in construction json", length(vec.trimmed))
   paste0(vec.trimmed, collapse=" ")

} # pre.allocate.graphToJSON
#------------------------------------------------------------------------------------------------------------------------
# follow this simplified form, in which neither "node" nor "edge" is neded:
# cy.json({elements:[
#      {data: {id:'a'}},
#      {data: {id: 'e1', source: 'a', target: 'a'}}
#      ]})
#
# or
#
#   network =  {"nodes":[{"data":{"id":"o"},"position":{"x":779.5, "y":500 }}],
#               "edges":[{"data":{"source":"o","target":"o","id":"e1"},"position":{}}]};
#
old.graphToJSON <- function(g)
{
   if(length(nodes(g)) == 0)
      return ("{}")

    x <- '{"nodes": [';
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
       x <- sprintf('%s{"data": {"id": "%s"', x, node);
       this.nodes.data <- nodeData(g, node)[[1]]
       if(length(this.nodes.data) > 0){
          noa.json <- toJSON(nodeData(g, node)[[1]], auto_unbox=TRUE)
          x <- sprintf("%s, %s", x, noa.json)
          }
       x <- sprintf('%s}', x)     # close off this node data element
       if(n != nodeCount)
           x <- sprintf("%s},", x)  # another node coming, add a comma
       } # for n

    x <- sprintf("%s}]", x)         # close off the last node, the node array ], the nodes element }

    if(edgeCount > 0){
       x <- sprintf('%s, "edges": [', x)
       for(e in seq_len(edgeCount)) {
          edgeName <- edgeNames[e]
          edge <- edges[[e]]
          sourceNode <- edge[[1]]
          targetNode <- edge[[2]]
          x <- sprintf('%s {"data": {"id": "%s", "source": "%s", "target": "%s"', x, edgeName, sourceNode, targetNode);
          x <- sprintf('%s}}', x)     # close off this edge data element
          if(e != edgeCount)          # add a comma, ready for the next edge element
             x <- sprintf('%s,', x)
          } # for e
       x <- sprintf("%s]", x)         # edge elements now complete, close the array
        } # if edgeCount > 0

   x <- sprintf("%s}", x)

   x

} # old.graphToJSON
#------------------------------------------------------------------------------------------------------------------------
simpleDemoGraph = function ()
{
  g = new ('graphNEL', edgemode='directed')

  nodeDataDefaults(g, attr='type') <- 'undefined'
  nodeDataDefaults(g, attr='lfc') <-  1.0
  nodeDataDefaults(g, attr='label') <- 'default node label'
  nodeDataDefaults(g, attr='count') <-  0

  edgeDataDefaults(g, attr='edgeType') <- 'undefined'
  edgeDataDefaults(g, attr='score') <-  0.0
  edgeDataDefaults(g, attr= 'misc') <- "default misc"

  g = graph::addNode ('A', g)
  g = graph::addNode ('B', g)
  g = graph::addNode ('C', g)
  nodeData (g, 'A', 'type') = 'kinase'
  nodeData (g, 'B', 'type') = 'transcription factor'
  nodeData (g, 'C', 'type') = 'glycoprotein'

  nodeData (g, 'A', 'lfc') = -3.0
  nodeData (g, 'B', 'lfc') = 0.0
  nodeData (g, 'C', 'lfc') = 3.0

  nodeData (g, 'A', 'count') = 2
  nodeData (g, 'B', 'count') = 30
  nodeData (g, 'C', 'count') = 100

  nodeData (g, 'A', 'label') = 'Gene A'
  nodeData (g, 'B', 'label') = 'Gene B'
  nodeData (g, 'C', 'label') = 'Gene C'

  g = graph::addEdge ('A', 'B', g)
  g = graph::addEdge ('B', 'C', g)
  g = graph::addEdge ('C', 'A', g)

  edgeData (g, 'A', 'B', 'edgeType') = 'phosphorylates'
  edgeData (g, 'B', 'C', 'edgeType') = 'synthetic lethal'

  edgeData (g, 'A', 'B', 'score') =  35.0
  edgeData (g, 'B', 'C', 'score') =  -12

  g

} # simpleDemoGraph
#----------------------------------------------------------------------------------------------------
