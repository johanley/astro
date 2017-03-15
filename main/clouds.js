/*
 Show an image of local cloud cover on a canvas, using a satellite image.
 Always stretch the contrast: otherwise it's too low.
 Current impl uses Web Mapping Service (WMS) servers.
*/
var showClouds = function(input, current_solar_alt_degs, output_image, satellite_name){

  /* these are all required */
  var size_pixels = parseInt(input.pixels_on_a_side);
  var size_degrees = parseFloat(input.degrees_on_a_side);
  var latitude = parseFloat(input.latitude);
  var longitude = parseFloat(input.longitude);
  var locations; //optional
  
  var STRETCH_CONTRAST_ON = true;  
  var STRETCH_CONTRAST_OFF = true;

  var parseInput = function(){
    locations = input.locations ? input.locations : null;
  };
  
  /** The bbox is square, and centered on the input location. */
  var calcBoundingBox = function(){
    var delta = size_degrees/2;
    return {
      sw_long: longitude - delta,
      sw_lat: latitude - delta,
      ne_long: longitude + delta,
      ne_lat: latitude + delta
    };
  };

  /* The parts that vary from one WMS server to the next (so far, at any rate). */  
  var make_wms_server = function(preamble, customizer_fn){
    return {
      preamble: preamble,
      customizer_fn: customizer_fn /* takes a layer string (visible|ir), returns a string (to be added to the URL) */
    };
  };
  var goes_customizer = function(layer){
    var result = '';
    var longitude_of_switchover = -101.38; //CAREFUL! the switchover is really diagonal: http://www.goes.noaa.gov/goes-w.html
    var half = longitude < longitude_of_switchover ? 'west' : 'east';
    if ('visible' === layer){
      result = half +'_vis.map&LAYERS=' + half + '_vis_1km'; 
    }
    else if ('ir' === layer){
      result = half + '_ir.map&LAYERS=' + half + '_ir_4km_gray'; 
    }
    return result;
  };
  var meteosat_customizer = function(layer){
    var result = '&LAYERS=meteosat:msg_';
    return ('visible' === layer) ?  result + 'vis006' : result + 'ir108'; 
  };
  var wms_servers = function(){
    var result = {};
    result.goes = make_wms_server(
      'http://mesonet.agron.iastate.edu/cgi-bin/mapserv/mapserv?map=/mesonet/www/apps/iemwebsite/data/wms/goes/',
       goes_customizer  
    )
    result.meteosat = make_wms_server(
      'http://eumetview.eumetsat.int/geoserv/wms?',
      meteosat_customizer
    );
    return result;
  };
  var wms_server_from = function(satellite_name){
    return wms_servers()[satellite_name];
  }
  var layer = function(){
    var result = input.layer;
    if ('auto_detect' === result){
      //overwrite the value with the 'real' one
      if (current_solar_alt_degs > 0){ //there's ~15m delay in getting images from space
        result = 'visible';
      }
      else {
        result = 'ir';
      }
    }
    return result;
  };
  /*
   GOES example: 
    ir: http://mesonet.agron.iastate.edu/cgi-bin/mapserv/mapserv?map=/mesonet/www/apps/iemwebsite/data/wms/goes/east_ir.map&LAYERS=east_ir_4km_gray&BBOX=-75,40,-73,41&WIDTH=100&HEIGHT=100&REQUEST=GetMap&SERVICE=WMS&VERSION=1.1.1&STYLES=&FORMAT=image/png&BGCOLOR=0xFFFFFF&TRANSPARENT=TRUE&SRS=EPSG:4326
   METEOSAT example: 
    vis: http://eumetview.eumetsat.int/geoserv/wms?SERVICE=WMS&REQUEST=GetMap&SERVICE=WMS&VERSION=1.1.1&LAYERS=meteosat:msg_vis006&STYLES=&FORMAT=image/png&BGCOLOR=0xFFFFFF&TRANSPARENT=TRUE&SRS=EPSG:4326&BBOX=-11.1523918226068,30.6418455813954,0.48155275283938,40.4508576744186&WIDTH=817&HEIGHT=860
    ir: http://eumetview.eumetsat.int/geoserv/wms?SERVICE=WMS&REQUEST=GetMap&SERVICE=WMS&VERSION=1.1.1&LAYERS=meteosat:msg_ir108&STYLES=&FORMAT=image/png&BGCOLOR=0xFFFFFF&TRANSPARENT=TRUE&SRS=EPSG:4326&BBOX=-11.1523918226068,30.6418455813954,0.48155275283938,40.4508576744186&WIDTH=817&HEIGHT=860    
   The only differences between the two are the preamble, and the identifier for the layer.
  */  
  var calcUrlClouds = function(bbox, layer, satellite_name){
    var wms_server = wms_server_from(satellite_name);
    var result = wms_server.preamble;
    result = result + wms_server.customizer_fn(layer());
    result = result + '&BBOX=' + bbox.sw_long+',' + bbox.sw_lat+',' + bbox.ne_long+',' + bbox.ne_lat; 
    result = result + '&WIDTH=' + size_pixels;
    result = result + '&HEIGHT=' + size_pixels;
    result = result + '&REQUEST=GetMap&SERVICE=WMS&VERSION=1.1.1&STYLES=&FORMAT=image/png&BGCOLOR=0xFFFFFF&TRANSPARENT=TRUE&SRS=EPSG:4326';
    console.log("Clouds url: " + result);
    return result;
  };
  
  var setCanvasSizeToMatchImageSize = function(canvas){
    canvas.height = size_pixels;
    canvas.width = size_pixels;
  };
  
  var supports_cors = function(){
    return satellite_name === 'goes';
  };

  /** 
   Showing the image in a canvas instead of an img tag let's you do a contrast stretch on the pixels.
   The server returning the image needs to support CORS, in this impl. 
  */
  var showTheImage = function(url, output, stretch_contrast){
    var canvas = document.getElementById(output);
    setCanvasSizeToMatchImageSize(canvas);
    var ctx = canvas.getContext("2d");
    var img = new Image();
    img.src = '';
    //how to respond later, after the fetch finishes
    img.onload = function() {
      console.log("Success: cloud image.");
      drawImage(this, ctx, canvas, stretch_contrast);
      drawSpecificLocations(ctx);
    };
    img.onerror = function(){
      console.log("Error loading cloud image.");
    };
    //wacky: need to add this function; otherwise the onload doesn't fire!! why??
    img.abort = function(){
      console.log("Aborted loading cloud image.");
    };
    if (supports_cors()) {
      img.crossOrigin = 'anonymous'; //crossOrigin (property), not crossorigin (attr)
      img.style.display = 'none';
      img.src = url; //starts the fetch of the image bytes
    }
    else {
      img.style.display = 'none';
      img.src = UTIL.crossDomainUrl(encodeURIComponent(url) + '&ext=png'); //starts the fetch of the image bytes
    }
  };

  /** Don't do anything until the image comes back from the network. */
  var drawImage = function(img, ctx, canvas, stretch_contrast){
    //console.log('The raw image has been fetched from the network.');
    ctx.drawImage(img, 0, 0);
    if (stretch_contrast){
      //console.log('Applying a linear contrast stretch to the raw image, since it is a bit dark.');
      var imageData = ctx.getImageData(0,0,canvas.width, canvas.height); 
      var data = imageData.data;
      for (var i = 0; i < data.length; i += 4) {
        //this is an 8-bit grey-scale image 0..255; r=g=b
        //var stretched = 2*data[i];
        var stretched = contrastStretch(data[i]);
        data[i] = stretched; // red
        data[i + 1] = stretched; // green
        data[i + 2] = stretched; // blue
      }
      //draws the image a second time, with the altered pixels    
      ctx.putImageData(imageData, 0, 0);
    }
    else {
      //console.log('Not applying a linear contrast stretch.');
    }
  };
  
  /**  Piecewise-linear contrast stretch.  */
  var contrastStretch = function(val){
    var result = val;
    var LOW = 75;
    var HIGH = 150;
    //this seems to show up the light cloud in a useful way!
    if(0 <= val && val <= LOW) {
      //lower values are doubled
      result = val * 2; //this is a fast computation
    }
    else if (LOW < val){
      //higher values
      result = HIGH + (val-LOW)/2; // never saturates; max 240
    }
    return result;
  };

  /** Return an array of objects. Each object has latitude, longitude. */
  var parse_locations = function(){
    var result = []; 
    var all_locations = locations.split(";"); //separate the locations
    var i;
    for (i=0; i<all_locations.length; i++){
      var lat_long = all_locations[i].split(","); //separate lat from long
      var position = {
        latitude: parseFloat(lat_long[0]),
        longitude: parseFloat(lat_long[1])
      };
      result.push(position);
    }
    //console.log('Found this many specific positions: ' + result.length);
    return result;
  };

  /** The NW corner corresponds to the origin of the drawing canvas. */
  var find_north_west_corner = function(){
    var delta = size_degrees/2;
    return {
      longitude: longitude - delta,
      latitude: latitude + delta
    };
  };
  
  var drawSpecificLocations = function(ctx){
    ctx.fillStyle = 'rgb(0,200,0)'; 
    //console.log('Drawing the center of the square.');
    drawLocation(0.5, 0.5, ctx);
    if (locations !== null){
      //console.log('Drawing other locations.');
      var nw = find_north_west_corner();
      var parsed_locations = parse_locations();
      var i;
      for(i=0; i<parsed_locations.length; i++){
        var width_frac = Math.abs(parsed_locations[i].longitude - nw.longitude)/size_degrees;
        var height_frac = Math.abs(nw.latitude - parsed_locations[i].latitude)/size_degrees; 
        drawLocation(width_frac, height_frac, ctx); 
      }
    }
  };
    
  var drawLocation = function(width_frac, height_frac, ctx){
    GRAPH.spot(ctx, width_frac*size_pixels, height_frac*size_pixels, 2);
  };

  var doAll = function(){
    parseInput();
    var bbox = calcBoundingBox();
    
    var cloudsUrl = calcUrlClouds(bbox, layer, satellite_name);
    //console.log('Clouds URL: ' + cloudsUrl);
    //console.log('Showing the clouds image in a canvas.');
    showTheImage(cloudsUrl, output_image, STRETCH_CONTRAST_ON);
  };
  
  return doAll();  
};