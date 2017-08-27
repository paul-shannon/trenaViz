"use strict";
import css from './css/trenaviz.css';

//----------------------------------------------------------------------------------------------------
var TrenaViz = (function(){

   var hub;                     // defined in BrowserViz.js, has lots of helpful socket
                                // and message support
   var bvDemo;                  // this simple webapp
   var igvBrowser = null;

//----------------------------------------------------------------------------------------------------
function setHub(newHub)
{
   hub = newHub;

} // setHub
//----------------------------------------------------------------------------------------------------
function addMessageHandlers()
{
   var self = this;  // the context of the current object, Trenaviz

   //var bound_respondToPing = respondToPing.bind(self);
   //hub.addMessageHandler("ping",  bound_respondToPing)
   hub.addMessageHandler("ping",  respondToPing.bind(self));

   hub.addMessageHandler("setGenome",     setGenome.bind(self));


} // addMessageHandlers
//----------------------------------------------------------------------------------------------------
// called out of the hub once the web page (the DOM) is ready (fully loaded)
function initializeUI()
{
  console.log("=== STARTING 527 inst/browserCode/src/trenaviz.js initializeUI");

  var trenaVizDiv = $("#trenaVizDiv");

  console.log("about to call tabs");
  setTimeout(function() {$("#trenaVizDiv").tabs()}, 0);


   var bound_handleWindowResize = this.handleWindowResize.bind(this);
   setTimeout(function(){bound_handleWindowResize();}, 250)
   $(window).resize(bound_handleWindowResize);

   console.log("=== ending inst/browserCode/src/trenaviz.js initializeUI");

}  // initializeUI
//----------------------------------------------------------------------------------------------------
function handleWindowResize ()
{
   console.log("--- starting handleWindowResize");
   console.log(this);
   var tabsDiv = $("#trenaVizDiv");
   //tabsDiv.width(0.95 * $(window).width());
   tabsDiv.height(0.90 * $(window).height());

} // handleWindowResize
//--------------------------------------------------------------------------------
function respondToPing (msg)
{
   console.log("==== entering trenaviz.js::respondToPing, msg:");
   console.log(msg)
   var return_msg = {cmd: msg.callback, status: "success", callback: "", payload: "pong"};
   hub.send(return_msg);

} // respondToPing
//------------------------------------------------------------------------------------------------------------------------
function setGenome(msg)
{
  var supportedGenomes = ["hg38", "mm10"];
  var genomeName = msg.payload;
  var returnPayload = "";

  if(supportedGenomes.indexOf(genomeName) < 0){
     status = "failure"
     returnPayload = "error, unsupported genome: '" + genomeName + "'";
     var return_msg = {cmd: msg.callback, status: status, callback: "", payload: returnPayload};
     hub.send(return_msg);
     } // if unsupported genome

   setTimeout(function(){this.igvBrowser = initializeIGV(genomeName);}, 0);
   hub.send({cmd: msg.callback, status: "success", callback: "", payload: ""});

} // setGenome
//----------------------------------------------------------------------------------------------------
function initializeIGV(genomeName)
{
   console.log("--- trenaViz, initializeIGV");

    var hg38_options = {
	 locus: "MEF2C",
         showRuler: true,

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
         console.log("chromLocString: " + chromLocString);
         //trenaViz.chromLocString = chromLocString;
         });

   console.log("after igv createBrowser");

} // initializeIGV
//-----------------------------------------------------------------------------------------------------
  return({
    setHub: setHub,
    addMessageHandlers: addMessageHandlers,
    initializeUI: initializeUI,
    handleWindowResize: handleWindowResize,
    });

}); // TrenaViz
//--------------------------------------------------------------------------------
var hub = require("browservizjs")
var bvDemo = TrenaViz();
window.bvd = bvDemo;

bvDemo.setHub(hub)
hub.init();

bvDemo.addMessageHandlers()

var bound_initializeUI = bvDemo.initializeUI.bind(bvDemo);
hub.addOnDocumentReadyFunction(bound_initializeUI);

hub.start();

//--------------------------------------------------------------------------------
