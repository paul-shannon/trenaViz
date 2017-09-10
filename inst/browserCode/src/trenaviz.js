"use strict";
var cytoscape = require('cytoscape');
import css from './css/trenaviz.css';
//----------------------------------------------------------------------------------------------------
var TrenaViz = (function(hub){

  var hub = hub;

//----------------------------------------------------------------------------------------------------
// development aid: used to ensure that the right "this" -- corresponding to the TrenaViz object --
// is avaialable when needed
function checkSignature(obj, callersName)
{
   var success = false;  // be pessimistic
   if(Object.keys(obj).indexOf("signature") >= 0 && obj.signature.indexOf("TrenaViz" == 0)){
      success = true;
      }

   if(!success){
      console.log("--- error: not a TrenaViz object: " + callersName);
      console.log(JSON.stringify(Object.keys(obj)))
      throw new Error("object is not a TrenaViz this!");
      }

} // checkSignature
//----------------------------------------------------------------------------------------------------
function addMessageHandlers()
{
   var self = this;  // the context of the current object, TrenaViz
   checkSignature(self, "addMessageHandlers");

   self.hub.addMessageHandler("ping",               respondToPing.bind(self));
   self.hub.addMessageHandler("raiseTab",           raiseTab.bind(self));
   self.hub.addMessageHandler("setGenome",          setGenome.bind(self));
   self.hub.addMessageHandler("setGraph",           setGraph.bind(self));
   self.hub.addMessageHandler("setStyle",           setStyle.bind(self));

   self.hub.addMessageHandler("getTrackNames",      getTrackNames.bind(self));
   self.hub.addMessageHandler("removeTracksByName", removeTracksByName.bind(self));

   self.hub.addMessageHandler("showGenomicRegion",  showGenomicRegion.bind(self));
   self.hub.addMessageHandler("getGenomicRegion",   getGenomicRegion.bind(self));

   self.hub.addMessageHandler("getSelectedNodes",   getSelectedNodes.bind(self));
   self.hub.addMessageHandler("selectNodes",        selectNodes.bind(self));

   self.hub.addMessageHandler("fit",                fit.bind(self));
   self.hub.addMessageHandler("fitSelected",        fitSelected.bind(self));

   self.hub.addMessageHandler("addBedTrackFromDataFrame",  addBedTrackFromDataFrame.bind(self));
   self.hub.addMessageHandler("addBedTrackFromHostedFile", addBedTrackFromHostedFile.bind(self));

   self.hub.addMessageHandler("addBedGraphTrackFromDataFrame",  addBedGraphTrackFromDataFrame.bind(self));



} // addMessageHandlers
//----------------------------------------------------------------------------------------------------
// called out of the hub once the web page (the DOM) is ready (fully loaded).
// tv(this) is explicitly bound to this function
//   1. create tabs
//   2. window resize handler is bound and assignes
function initializeUI()
{
   var self = this;
   checkSignature(self, "initializeUI");

   var trenaVizDiv = $("#trenaVizDiv");

   var activateFunction = function(event, ui){
      if(ui.newPanel.is("#cyOuterDiv")){
        console.log("cy!");
        self.handleWindowResize();
        if(self.cyjs != null){
           self.cyjs.resize();
	   }
        } // cyOuterDiv
      else if(ui.newPanel.is("#igvOuterDiv")){
         console.log("igv!");
         }
      else{
         console.log("unrecognized tab activated");
	 }
      }; // activateFunction

   var tabOptions = {activate: activateFunction};
   setTimeout(function() {$("#trenaVizDiv").tabs(tabOptions)}, 0);

   var bound_handleWindowResize = this.handleWindowResize.bind(self);
   setTimeout(function(){bound_handleWindowResize();}, 250)
   $(window).resize(bound_handleWindowResize);

}  // initializeUI
//----------------------------------------------------------------------------------------------------
function handleWindowResize ()
{
   var tabsDiv = $("#trenaVizDiv");

   tabsDiv.width(0.98  * $(window).width());
   tabsDiv.height(0.92 * $(window).height());

   $("#cyDiv").width($("#cyMenubarDiv").width()) // Width0.92 * tabsDiv.width());
   $("#cyDiv").height(tabsDiv.height() - 3 * $("#cyMenubarDiv").height()); //tabsDiv.height()-100);

   $("#igvOuterDiv").height($("#trenaVizDiv").height() - (3 * $(".ui-tabs-nav").height()))

} // handleWindowResize
//--------------------------------------------------------------------------------
function raiseTab(msg)
{
  var displayedTabName = msg.payload;

  var status = "success"  // be optimistic
  var returnPayload = "";

  switch(displayedTabName) {
    case "IGV":
       $('a[href="#igvOuterDiv"]').click();
       break;
    case "TRN":
       $('a[href="#cyOuterDiv"]').click();
       break;
    default:
       status = "error";
       returnPayload = "unrecognized tab name: " + displayedTabName;
    } // switch on displayedTabName

  var return_msg = {cmd: msg.callback, status: status, callback: "", payload: returnPayload};

  hub.send(return_msg);

} // raiseTab
//----------------------------------------------------------------------------------------------------
function respondToPing (msg)
{
   var self = this;
   checkSignature(self, "respondToPing")

   var return_msg = {cmd: msg.callback, status: "success", callback: "", payload: "pong"};
   self.hub.send(return_msg);

} // respondToPing
//------------------------------------------------------------------------------------------------------------------------
function setGenome(msg)
{
   var self = this;
   checkSignature(self, "setGenome")

   var supportedGenomes = ["hg38", "mm10"];
   var genomeName = msg.payload;
   var returnPayload = "";

   if(supportedGenomes.indexOf(genomeName) < 0){
      status = "failure"
      returnPayload = "error, unsupported genome: '" + genomeName + "'";
      var return_msg = {cmd: msg.callback, status: status, callback: "", payload: returnPayload};
      hub.send(return_msg);
      } // if unsupported genome

    $('a[href="#igvOuterDiv"]').click();
    setTimeout(function(){self.igvBrowser = initializeIGV(self, genomeName);}, 0);
    self.hub.send({cmd: msg.callback, status: "success", callback: "", payload: ""});

} // setGenome
//----------------------------------------------------------------------------------------------------
function initializeIGV(self, genomeName)
{
   console.log("--- trenaViz, initializeIGV");

   checkSignature(self, "initializeIGV")

    var hg38_options = {
	//locus: "MEF2C",
     flanking: 1000,
     showRuler: true,
     minimumbBases: 10,

     reference: {id: "hg38",
                 //fastaURL: "http://igv.broadinstitute.org/genomes/seq/1kg_v37/human_g1k_v37_decoy.fasta",
                 // cytobandURL: "http://igv.broadinstitute.org/genomes/seq/b37/b37_cytoband.txt"
                 },
     tracks: [
        {name: 'Gencode v24',
         //url: "http://trena.systemsbiology.net/hg38/gencode.v24.annotation.sorted.gtf.gz",
         //indexURL: "http://trena.systemsbiology.net/hg38/gencode.v24.annotation.sorted.gtf.gz.tbi",
         url: "//s3.amazonaws.com/igv.broadinstitute.org/annotations/hg38/genes/gencode.v24.annotation.sorted.gtf.gz",
         indexURL: "//s3.amazonaws.com/igv.broadinstitute.org/annotations/hg38/genes/gencode.v24.annotation.sorted.gtf.gz.tbi",
         format: 'gtf',
         visibilityWindow: 2000000,
         displayMode: 'EXPANDED'
         },
        ]
     }; // hg38_options


   var mm10_options = {//locus: "5:88,621,548-88,999,827", //"22:40,000,000-40,200,000",
         flanking: 2000,
	 showKaryo: false,
         showNavigation: true,
         minimumBases: 5,
         showRuler: true,
         reference: {id: "mm10",
                     fastaURL: "http://trena.systemsbiology.net/mm10/GRCm38.primary_assembly.genome.fa",
                     cytobandURL: "http://trena.systemsbiology.net/mm10/cytoBand.txt"
                     },
         tracks: [
            {name: 'Gencode vM14',
             url: "http://trena.systemsbiology.net/mm10/gencode.vM14.basic.annotation.sorted.gtf.gz",
             indexURL: "http://trena.systemsbiology.net/mm10/gencode.vM14.basic.annotation.sorted.gtf.gz.tbi",
             format: 'gtf',
             visibilityWindow: 2000000,
             displayMode: 'EXPANDED'
             },
            ]
       }; // mm10_options

   var igvOptions = null;

   switch(genomeName) {
      case "hg38":
         igvOptions = hg38_options;
         break;
       case "mm10":
         igvOptions = mm10_options;
         break;
         } // switch on genoneName

    $("#igvDiv").children().remove()

   console.log("--- trenaViz, igv:");
   console.log(igv)
   console.log("about to createBrowser");

   var igvBrowser = igv.createBrowser($("#igvDiv"), igvOptions);

   igvBrowser.on("locuschange",
       function(referenceFrame, chromLocString){
         self.chromLocString = chromLocString;
         });

   return(igvBrowser);

} // initializeIGV
//----------------------------------------------------------------------------------------------------
function showGenomicRegion(msg)
{
   var self = this;
   checkSignature(self, "showGenomicRegion")

   var regionString = msg.payload.regionString;
   self.igvBrowser.search(regionString)

   self.hub.send({cmd: msg.callback, status: "success", callback: "", payload: ""});

} // showGenomicRegion
//----------------------------------------------------------------------------------------------------
function getGenomicRegion(msg)
{
   var self = this;
   checkSignature(self, "getGenomicRegion")

   self.hub.send({cmd: msg.callback, status: "success", callback: "", payload: this.chromLocString});

} // getGenomicRegion
//----------------------------------------------------------------------------------------------------
function getTrackNames(msg)
{
   var self = this;
   checkSignature(self, "getTrackNames");

   var result = [];
   var count = self.igvBrowser.trackViews.length;

   for(var i=0; i < count; i++){
      var trackName = self.igvBrowser.trackViews[i].track.name;
      if(trackName.length > 0){
         result.push(trackName)
	 }
      } // for i

   self.hub.send({cmd: msg.callback, status: "success", callback: "", payload: result});

} // getTrackNames
//----------------------------------------------------------------------------------------------------
function removeTracksByName(msg)
{
   var self = this;
   checkSignature(self, "removeTracksByName")

   var trackNames = msg.payload;
   if(typeof(trackNames) == "string")
      trackNames = [trackNames];

   var count = self.igvBrowser.trackViews.length;

   for(var i=(count-1); i >= 0; i--){
     var trackView = self.igvBrowser.trackViews[i];
     var trackViewName = trackView.track.name;
     var matched = trackNames.indexOf(trackViewName) >= 0;
     console.log(" is " + trackViewName + " in " + JSON.stringify(trackNames) + "? " + matched);
     if (matched){
        self.igvBrowser.removeTrack(trackView.track);
        } // if matched
     } // for i

   self.hub.send({cmd: msg.callback, status: "success", callback: "", payload: ""});


} // removeTracksByName
//----------------------------------------------------------------------------------------------------

function addBedTrackFromDataFrame(msg)
{
   var self = this;
   checkSignature(self, "addBedTrackFromDataFrame")

   console.log("=== addBedTrackFromDataFrame");
   console.log(JSON.stringify(msg));

   var trackName = msg.payload.name;
   var bedFileName = msg.payload.bedFileName;
   var displayMode = msg.payload.displayMode;
   var color = msg.payload.color;
   var url = window.location.href + "?" + bedFileName;

   var config = {format: "bed",
                 name: trackName,
                 url: url,
                 indexed:false,
                 displayMode: displayMode,
                 sourceType: "file",
                 color: color,
                 type: "annotation"};

   console.log(JSON.stringify(config));
   self.igvBrowser.loadTrack(config);

   self.hub.send({cmd: msg.callback, status: "success", callback: "", payload: ""});

} // addBedTrackFromDataFrame
//----------------------------------------------------------------------------------------------------
function addBedGraphTrackFromDataFrame(msg)
{
   var self = this;
   checkSignature(self, "addBedGraphTrackFromDataFrame")

   console.log("--- addBedGraphTrackFromDataFrame");
   console.log(msg.payload)

   var trackName = msg.payload.name;
   var bedFileName = msg.payload.bedFileName;
   var displayMode = msg.payload.displayMode;
   var color = msg.payload.color;
   var minValue = msg.payload.min
   var maxValue = msg.payload.max

   var url = window.location.href + "?" + bedFileName;

   var config = {format: "bedgraph",
                 name: trackName,
                 url: url,
                 min: minValue,
                 max: maxValue,
                 indexed:false,
                 displayMode: displayMode,
                 sourceType: "file",
                 color: color,
                 type: "wig"};

   self.igvBrowser.loadTrack(config);
   self.hub.send({cmd: msg.callback, status: "success", callback: "", payload: ""});

} // addBedGraphTrackFromDataFrame
//----------------------------------------------------------------------------------------------------
function addBedTrackFromHostedFile(msg)
{
   var self = this;
   checkSignature(self, "addBedTrackFromHostedFile")

   console.log("=== addBedTrackFromFile");

   var trackName = msg.payload.name;
   var displayMode = msg.payload.displayMode;
   var color = msg.payload.color;
   var uri       = msg.payload.uri;
   var indexUri  = msg.payload.indexUri;
   var indexed = true;

   if(indexUri==null){
     indexed = false;
     }

   var config = {format: "bed",
                 name: trackName,
                 url: uri,
                 indexed: indexed,
                 displayMode: displayMode,
                 color: color,
                 type: "annotation"};

   if(indexed){
     config.indexURL = indexUri;
     }

   self.igvBrowser.loadTrack(config);

   self.hub.send({cmd: msg.callback, status: "success", callback: "", payload: ""});

} // addBedTrackFromHostedFile
//----------------------------------------------------------------------------------------------------
function selectNodes(msg)
{
   var self = this;
   checkSignature(self, "selectNodes")

   console.log("==== selectNodes");
   console.log(msg.payload);
   var nodeIDs = msg.payload.nodeIDs;

   if(typeof(nodeIDs) == "string")
      nodeIDs = [nodeIDs];

   var filterStrings = [];

   for(var i=0; i < nodeIDs.length; i++){
     var s = '[id="' + nodeIDs[i] + '"]';
     filterStrings.push(s);
     } // for i


   var nodesToSelect = self.cyjs.nodes(filterStrings.join());
   nodesToSelect.select()

   self.hub.send({cmd: msg.callback, status: "success", callback: "", payload: ""});

} // selectNodes
//----------------------------------------------------------------------------------------------------
function getSelectedNodes(msg)
{
   var self = this;
   checkSignature(self, "getSelectedNodes")

   var status = "success";  // be optimistic
   var payload = "";

   if (typeof (this.cyjs) == "undefined"){
      payload = JSON.stringify([]);
      status = "error";
      }
   else if (this.cyjs.nodes().length == 0){
      payload = JSON.stringify([]);
      }
   else {
      payload =  JSON.stringify(self.cyjs.filter("node:selected").map(function(node) {
                                return {id: node.data().id, name: node.data().name}}));
      }

   console.log("getNodes returning payload: " + payload);
   self.hub.send({cmd: msg.callback, status: "success", callback: "", payload: payload});

} // getSelectedNodes
//---------------------------------------------------------------------------------------------------
function fit(msg)
{
   var self = this;
   checkSignature(self, "fit")

   var margin = msg.payload;
   self.cyjs.fit(margin)
   hub.send({cmd: msg.callback, status: "success", callback: "", payload: ""});

} // fit
//----------------------------------------------------------------------------------------------------
function fitSelected(msg)
{
   var self = this;
   checkSignature(self, "fit")

   var selectedNodes = self.cyjs.filter('node:selected');
   var margin = msg.payload;

   if(selectedNodes.length == 0){
     status = "failure";
     payload = "no nodes currently selected"
     }
  else{
    console.log("fitSelected, with margin " + margin);
    self.cyjs.fit(selectedNodes, margin)
    status = "success";
    payload = "";
    }

  hub.send({cmd: msg.callback, status: status, callback: "", payload: payload});

} // fit
//----------------------------------------------------------------------------------------------------
function setGraph(msg)
{
     // soon: graphs = msg.payload, though more complex, since names
     // of graphs will be sent as well

   var self = this;
   checkSignature(self, "setGraph")

   console.log("---> entering setGraph, this: ");
   console.log(this);
   $('a[href="#cyOuterDiv"]').click();
   self.handleWindowResize();
   self.cyjs = initializeTrnCytoscape();
   var temporaryFileName = msg.payload.filename;
   var modelNames = msg.payload.modelNames;
   console.log("modelNames: ");
   console.log(modelNames)
   if(modelNames.length > 1){
     createModelNamesMenu(self, modelNames);
     }


   var status = readNetworkFromFile(temporaryFileName, self.cyjs)
   initializeTrnCytoscapeButtons(self);

   setTimeout(function(){
     console.log("about to call that.fit, self: ");
     console.log(self);
     self.cyjs.fit(100);
     }, 500);

   self.hub.send({cmd: msg.callback, status: "success", callback: "", payload: ""});

} // setGraph
//----------------------------------------------------------------------------------------------------
function createModelNamesMenu(self, modelNames)
{
   if(typeof(modelNames) == "string"){
      modelNames = [modelNames];
      }

   if(modelNames.length < 1){
     return;
     }

  $("#cyModelSelector").remove()

   var html = "<select id='cyModelSelector'>"
   for(var i=0; i < modelNames.length; i++){
      html += "<option value='" + modelNames[i] + "'> " + modelNames[i]  + "</option>";
      } // for i
   html += "</select>"

   $("#cyMenubarDiv").append(html);
   $("#cyModelSelector").change(function(){
     var modelName =$(this).find("option:selected").val();
       self.nextCyModel(self, modelName);
     });

    self.nextCyModel(self, modelNames[0]);

} // createModelNamesMenu
//----------------------------------------------------------------------------------------------------
// not sure why bind does not work on this function, thus necessitating explicit passing of
// this as self.  (pshannon 9 sep 2017)
function nextCyModel(self, modelName)
{
   console.log("--- nextCyModel: " + modelName);
   console.log(self)
   checkSignature(self, "nextCyModel")

   self.cyjs.nodes().show()
   self.cyjs.nodes().filter(function(node){return node.data("type") == "TF"}).map(function(node){node.data({"randomForest": 0})})
   self.cyjs.nodes().filter(function(node){return node.data("type") == "TF"}).map(function(node){node.data({"pearson": 0})})

   var noaName = modelName + "." + "randomForest";
   self.cyjs.nodes("[type='TF']").map(function(node){node.data({"randomForest":  node.data(noaName)})})

   noaName = modelName + "." + "pearson";
   self.cyjs.nodes("[type='TF']").map(function(node){node.data({"pearson":       node.data(noaName)})})

     // now hide all the 0 randomForest TF nodes
   self.cyjs.nodes().filter(function(node){return(node.data("randomForest") == 0 && node.data("type") == "TF")}).hide()

} // nextCyModel
//----------------------------------------------------------------------------------------------------
function setStyle(msg)
{
  console.log("=== entering setStyle");

   var self = this;
   checkSignature(self, "setGraph")

   var filename = msg.payload;
   console.log("setStyle: '" + filename + "'");
   loadStyleFile(filename, self.cyjs);

   var return_msg = {cmd: msg.callback, status: "success", callback: "", payload: ""};
   hub.send(return_msg);

} // setStyle
//----------------------------------------------------------------------------------------------------
function initializeTrnCytoscape()
{
  var options = {container: $("#cyDiv"),
                 elements: {nodes: [{data: {id:'a'}}],
                            edges: [{data:{source:'a', target:'a'}}]},
                 style: cytoscape.stylesheet()
                 .selector('node').style({'background-color': '#ddd',
                                          'label': 'data(id)',
                                          'text-valign': 'center',
                                          'text-halign': 'center',
                                          'border-width': 1})
                 .selector('node:selected').style({'overlay-opacity': 0.2,
                                                   'overlay-color': 'gray'})
                 .selector('edge').style({'line-color': 'black',
                                          'target-arrow-shape': 'triangle',
                                          'target-arrow-color': 'black',
                                          'curve-style': 'bezier'})
                 .selector('edge:selected').style({'overlay-opacity': 0.2,
                                                   'overlay-color': 'gray'})
                };

    console.log("about to call cytoscape with options");
    var cy = cytoscape(options);
    return(cy);

} // initializeTrnCytoscape
//----------------------------------------------------------------------------------------------------
function readNetworkFromFile(filename, targetCy)
{
   var s = window.location.href + "?" + filename;
   fetch(s)
      .then(function(responseObj){
          console.log("fetch in action");
          return responseObj.json();
          })
     .then(function(j){
         targetCy.json(j);
         return "success";
         });

    return "SUCCESS";

} // readNetworkFromFile
//----------------------------------------------------------------------------------------------------
// expected file contents:  vizmap = [{selector:"node",css: {...
function loadStyleFile(filename, cyTarget)
{
   console.log("trenaViz.loadStyleFile, filename: " + filename);

   var s = window.location.href + "?" + filename;
   console.log("=== about to getScript on " + s);

   $.getScript(s)
     .done(function(script, textStatus) {
        console.log(textStatus);
        //console.log("style elements " + layout.length);
        cyTarget.style(vizmap);
       })
    .fail(function( jqxhr, settings, exception ) {
       console.log("getScript error trying to read " + filename);
       console.log("exception: ");
       console.log(exception);
       });

} // loadStyle
//----------------------------------------------------------------------------------------------------
function initializeTrnCytoscapeButtons(self)
{
   checkSignature(self, "intializeTrnCytoscapeButtons")

   $("#cyFitButton").click(function(){self.cyjs.fit(50)});
   $("#cyFitSelectedButton").click(function(){self.cyjs.fit(self.cyjs.nodes(":selected"), 50)});

   $("#cySFNButton").click(function(){self.cyjs.nodes(':selected').neighborhood().nodes().select()});

   $("#cyHideUnselectedButton").click(function(){self.cyjs.nodes(":unselected").hide()});
   $("#cyShowAllButton").click(function(){self.cyjs.nodes().show(); self.cyjs.edges().show()});

   //$("#cyCycleThroughModelsButton").click(function(){nextCyModel("rs3875089")});

} // initializeTrnCytoscapeButtons
//-----------------------------------------------------------------------------------------------------
  return({

    signature: "TrenaViz 0.99.5",

    addMessageHandlers: addMessageHandlers,
    initializeUI: initializeUI,
    handleWindowResize: handleWindowResize.bind(this),
    initializeTrnCytoscape: initializeTrnCytoscape,
    nextCyModel: nextCyModel.bind(this),
    hub: hub,
    cyjs: null,
    igvBrowser: null,
    chromLocString: null
    });

}); // TrenaViz
//----------------------------------------------------------------------------------------------------
var hub = require("browservizjs")
var tv = TrenaViz(hub);
hub.init();
tv.addMessageHandlers()
hub.addOnDocumentReadyFunction(tv.initializeUI.bind(tv));
hub.start();
window.tv = tv;
