"use strict";
//require('jquery')
//require('jquery-ui-dist/jquery-ui.js');
//require("jquery-ui/ui/widgets/tabs");

import css from './css/trenaviz.css';

//----------------------------------------------------------------------------------------------------
var TrenaViz = (function(){

   var hub;                     // defined in BrowserViz.js, has lots of helpful socket
                                // and message support
   var bvDemo;                  // this simple webapp

//----------------------------------------------------------------------------------------------------
function setHub(newHub)
{
   hub = newHub;

} // setHub
//----------------------------------------------------------------------------------------------------
function addMessageHandlers()
{
   var self = this;  // the context of the current object, Trenaviz

   var bound_respondToPing = respondToPing.bind(self);
   hub.addMessageHandler("ping", bound_respondToPing)


} // addMessageHandlers
//----------------------------------------------------------------------------------------------------
// called out of the hub once the web page (the DOM) is ready (fully loaded)
function initializeUI()
{
  console.log("=== STARTING 527 inst/browserCode/src/trenaviz.js initializeUI");

  var trenaVizDiv = $("#trenaVizDiv");

  console.log("about to call tabs");
   setTimeout(function() {$("#trenaVizDiv").tabs()}, 0);

  /************
  trenaVizDiv.tabs({
     activate: function(event, ui){
        if(ui.newPanel.is("#cyOuterDiv")){
          console.log("cy!");
          handleWindowResize();
          cy.resize();
          }
        else if(ui.newPanel.is("#igvDiv")){
          console.log("igv!");
          }
        else{
          console.log("unrecognized tab activated");
          }
        }
     });
    **********/

   var bound_handleWindowResize = this.handleWindowResize.bind(this);
   setTimeout(function(){bound_handleWindowResize();}, 250)
   $(window).resize(bound_handleWindowResize);

   console.log("=== ending inst/browserCode/src/trenaviz.js initializeUI");

}  // initializeUI
//----------------------------------------------------------------------------------------------------
function handleWindowResize ()
{
   var plotDiv = $("#trenavizDiv");
   plotDiv.width(0.95 * $(window).width());
   plotDiv.height(0.90 * $(window).height());

      // an easy way to rescale the canvas when the browser window size changes: just redraw
   if(this.dataReceived)
     this.d3plot(this.dataset, this.xMin, this.xMax, this.yMin, this.yMax);

} // handleWindowResize
//--------------------------------------------------------------------------------
function respondToPing (msg)
{
   console.log("==== entering trenaviz.js::respondToPing, msg:");
   console.log(msg)
   var return_msg = {cmd: msg.callback, status: "success", callback: "", payload: "pong"};
   hub.send(return_msg);

} // respondToPing
//--------------------------------------------------------------------------------

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
