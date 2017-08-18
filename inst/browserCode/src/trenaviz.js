import css from './css/trenaviz.css';
//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------
var Trenaviz = (function(){

   //var plotDiv;
   // var d3plotDiv;
   var selectedRegion;          // assigned out of brushReader function
   var dataReceived = false;    // when true, window resize does replot of data
   var dataset;                 // assigned from payload of incoming plotxy message
   var xMin, xMax, yMin, yMax;  // conveniently calculated in R, part of payload
   var datasetLength;

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
   //hub.addMessageHandler("ping", respondToPing)

   var bound_d3plotPrep = d3plotPrep.bind(self);
   hub.addMessageHandler("plotxy", bound_d3plotPrep)
   //hub.addMessageHandler("plotxy", d3plotPrep)

   //var bound_getSelection = getSelection.bind(self);
   //hub.addMessageHandler("getSelection", bound_getSelection)
   //hub.addMessageHandler("getSelection", getSelection)

} // addMessageHandlers
//----------------------------------------------------------------------------------------------------
// called out of the hub once the web page (the DOM) is ready (fully loaded)
function initializeUI()
{
   console.log("=== STARTING inst/browserCode/src/trenaviz.js initializeUI");
   //var plotDiv = $("#trenavizDiv");
   //var d3plotDiv = d3.select("#trenavizDiv");
   //console.log("div: ");
   //console.log(plotDiv)

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
function brushReader ()
{
  console.log("brushReader");
  selectedRegion = window.d3brush.extent()();
  var x0 = selectedRegion[0][0];
  var x1 = selectedRegion[1][0];
  var width = Math.abs(x0-x1);
  if(width < 0.001)
     selectedRegion = null

}; // d3PlotBrushReader
//--------------------------------------------------------------------------------
function respondToPing (msg)
{
   console.log("==== entering trenaviz.js::respondToPing, msg:");
   console.log(msg)
   var return_msg = {cmd: msg.callback, status: "success", callback: "", payload: "pong"};
   hub.send(return_msg);

} // respondToPing
//--------------------------------------------------------------------------------
function d3plotPrep (msg)
{
   console.log("==== entering trenaviz.js::d3plotPrep, msg:");
   console.log(msg)

   var payload = msg.payload;

      // assign global variables

   this.dataset = [];
   this.datasetLength = payload.x.length

   for(var i=0; i < this.datasetLength; i++){
     this.dataset.push({x: payload.x[i], y: payload.y[i]});
     }

   this.xMin = msg.payload.xMin;
   this.xMax = msg.payload.xMax;
   this.yMin = msg.payload.yMin;
   this.yMax = msg.payload.yMax;
   this.dataReceived = true;

   this.d3plot(this.dataset, this.xMin, this.xMax, this.yMin, this.yMax)

   var return_msg = {cmd: msg.callback, status: "success", callback: "", payload: ""};
   hub.send(return_msg);

} // d3plotPrep
//--------------------------------------------------------------------------------
function d3plot(dataset, xMin, xMax, yMin, yMax)
{
  var plotDiv = $("#trenavizDiv");

  var width = plotDiv.width();
  var height = plotDiv.height();

  var padding = 50;
  var xScale = d3.scaleLinear()
                 .domain([xMin,xMax])
                 .range([padding,width-padding]);

  var yScale = d3.scaleLinear()
                 .domain([yMin, yMax])
                 .range([height-padding, padding]); // note inversion


    // must remove the svg from a d3-selected object, not just a jQuery object
  var d3plotDiv = d3.select("#trenavizDiv");
  d3plotDiv.select("#plotSVG").remove();  // so that append("svg") is not cumulative

    // see https://stackoverflow.com/questions/38237747/how-do-i-apply-a-scale-to-a-d3-v4-0-brush
  /************
  var brushHeight = 1000;
  var d3brush = d3.brush()
                  .extent([[xMin, xMax], [yMin, yMax]])
                  .on("end", brushReader);
                  //.extent([[xScale.range()[0], 0], [xScale.range()[1], brushHeight]])

  window.d3brush = d3brush;
  ***********/
  var svg = d3.select("#trenavizDiv")
      .append("svg")
      .attr("id", "plotSVG")
      .attr("width", width)
      .attr("height", height)
      //.call(d3brush);

  var xAxis = d3.axisBottom()
                .scale(xScale);

  var yAxis = d3.axisLeft()
                .scale(yScale);

  var xTranslationForYAxis = xScale(0);
  var yTranslationForXAxis = yScale(10);

  var tooltip = d3plotDiv.append("div")
                   .attr("class", "tooltip")
                   .style("position", "absolute")
                   .style("z-index", "10")
                   .style("visibility", "hidden")
                   .text("a simple tooltip");

   svg.selectAll("circle")
      .data(dataset)
      .enter()
      .append("circle")
      .attr("cx", function(d){
         return xScale(d.x);
        })
      .attr("cy", function(d){
         return yScale(d.y);
         })
      .attr("r", function(d){
        return 5;
        })
      .style("fill", function(d){
        return "red";
        })
      .on("mouseover", function(d,i){   // no id assigned yet...
         tooltip.text(d.id);
         return tooltip.style("visibility", "visible");
         })
      .on("mousemove", function(){
         return tooltip.style("top",
                (d3.event.pageY-10)+"px").style("left",(d3.event.pageX+10)+"px");})
      .on("mouseout", function(){return tooltip.style("visibility", "hidden");});

   svg.append("g")
      .attr("class", "axis")
      .attr("transform", "translate(0," + (height - padding) + ")")
      .call(xAxis);

   svg.append("g")
      .attr("class", "axis")
      .attr("transform", "translate(" + padding + ",0)")
      .call(yAxis);

} // d3plot
//--------------------------------------------------------------------------------
function getSelection(msg)
{
   var selectedNames = [];
   var selectedRegion = window.d3brush.extent()();
   console.log("---- getSelection, selectedRegion");
   console.log(selectedRegion)

   if(selectedRegion == null){
      var returnMsg = {cmd: msg.callback, callback: "", status: "success",
                        payload: selectedNames};
      console.log("=== returning from getSelection, no selected points");
      console.log(returnMsg);
      hub.send(returnMsg);
      return;
      }

   var x0 = selectedRegion[0][0];
   var x1 = selectedRegion[1][0];
   var y0 = selectedRegion[0][1];
   var y1 = selectedRegion[1][1];

   for(var i=0; i < this.datasetLength; i++){
      var x = this.dataset[i].x;
      var y = this.dataset[i].y;
      if(x >= x0 & x <= x1 & y >= y0 & y <= y1) {
        console.log ("TRUE");
        selectedNames.push("point " + i);
        }
     } // for i

   console.log(" found " + selectedNames.length + " selected points");
   var returnMsg = {cmd: msg.callback, callback: "", status: "success",
                    payload: selectedNames};

   console.log("=== returning from getSelection, 2");
   console.log(returnMsg);
   hub.send(returnMsg);

} // getSelection
//--------------------------------------------------------------------------------

  return({
    setHub: setHub,
    addMessageHandlers: addMessageHandlers,
    initializeUI: initializeUI,
    d3plot: d3plot,
    handleWindowResize: handleWindowResize,
    dataReceived: dataReceived,
    dataset: dataset,
    datasetLength: datasetLength,
    xMin: xMin,
    xMax: xMax,
    yMin: yMin,
    yMax: yMax
    });

}); // Trenaviz
//--------------------------------------------------------------------------------
var hub = require("browservizjs")
var bvDemo = Trenaviz();
window.bvd = bvDemo;

bvDemo.setHub(hub)
hub.init();

bvDemo.addMessageHandlers()

  // mysteriously:  only by assigning the module function to 'f' can it be
  // successfully passed to the hub and called when the document is ready

var bound_initializeUI = bvDemo.initializeUI.bind(bvDemo);
hub.addOnDocumentReadyFunction(bound_initializeUI);

hub.start();

//--------------------------------------------------------------------------------
