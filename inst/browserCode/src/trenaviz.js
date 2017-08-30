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
   self.hub.addMessageHandler("setGenome",          setGenome.bind(self));
   self.hub.addMessageHandler("setGraph",           setGraph.bind(self));

   self.hub.addMessageHandler("showGenomicRegion",  showGenomicRegion.bind(self));
   self.hub.addMessageHandler("getGenomicRegion",   getGenomicRegion.bind(self));

   self.hub.addMessageHandler("getSelectedNodes",   getSelectedNodes.bind(self));
   self.hub.addMessageHandler("selectNodes",        selectNodes.bind(self));

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
         //url: "http://pshannon.systemsbiology.net/hg38/gencode.v24.annotation.sorted.gtf.gz",
         //indexURL: "http://pshannon.systemsbiology.net/hg38/gencode.v24.annotation.sorted.gtf.gz.tbi",
         url: "//s3.amazonaws.com/igv.broadinstitute.org/annotations/hg38/genes/gencode.v24.annotation.sorted.gtf.gz",
         indexURL: "//s3.amazonaws.com/igv.broadinstitute.org/annotations/hg38/genes/gencode.v24.annotation.sorted.gtf.gz.tbi",
         format: 'gtf',
         visibilityWindow: 2000000,
         displayMode: 'EXPANDED'
         },
        ]
     }; // hg38_options


   var mm10_options = {locus: "5:88,621,548-88,999,827", //"22:40,000,000-40,200,000",
         flanking: 2000,
	 showKaryo: true,
         showNavigation: true,
         minimumBases: 5,
         showRuler: true,
         reference: {id: "mm10",
                     fastaURL: "http://pshannon.systemsbiology.net/mm10/GRCm38.primary_assembly.genome.fa",
                     cytobandURL: "http://pshannon.systemsbiology.net/mm10/cytoBand.txt"
                     },
         tracks: [
            {name: 'Gencode vM14',
             url: "http://pshannon.systemsbiology.net/mm10/gencode.vM14.basic.annotation.sorted.gtf.gz",
             indexURL: "http://pshannon.systemsbiology.net/mm10/gencode.vM14.basic.annotation.sorted.gtf.gz.tbi",
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
   initializeTrnCytoscapeButtons(self);

   setTimeout(function(){
     console.log("about to call that.fit, self: ");
     console.log(self);
     self.cyjs.fit(100);
     }, 500);

   self.hub.send({cmd: msg.callback, status: "success", callback: "", payload: ""});

} // setGraph
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
