buildMultiModelGraph <- function (targetGene, models)
{

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

       #printf("---- tfs.in.model: ")
       #print(tfs.in.model)

       #printf("---- tfs.in.regions: ")
       #print(tfs.in.regions)

       #printf("---- setdiff(tfs.in.model, tfs.in.regions)")
       #print(setdiff(tfs.in.model, tfs.in.regions))

       #printf("---- setdiff(tfs.in.regions, tfs.in.model)")
       #print(setdiff(tfs.in.regions, tfs.in.model))

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

} # buildMultiModelGraph
#----------------------------------------------------------------------------------------------------
addGeneModelLayout <- function (g, xPos.span=1500)
{
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

} # addGeneModelLayout
#------------------------------------------------------------------------------------------------------------------------
