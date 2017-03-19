/* 
 Show satellite images. 
 Current impl talks only to WMS servers.  
*/
var SATELLITE = (function(){ 

  var STRETCH_CONTRAST_ON = true;  
  var STRETCH_CONTRAST_OFF = true;

  /* The bounding box is always square, and centered on the input location. */
  var calc_bounding_box = function(size_degrees, latitude, longitude){
    var delta = size_degrees/2;
    return {
      sw_long: longitude - delta,
      sw_lat: latitude - delta,
      ne_long: longitude + delta,
      ne_lat: latitude + delta
    };
  };
  
  /* Return a simple string id for the layer, (visible|ir). */
  var layer_id = function(input_layer, current_solar_alt_degs){
    var result = input_layer;
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

  var is_in_range = function(val, limits){
    return val >= limits.min && val <= limits.max;
  };
  
  /* 
   If returns empty string, then an image cannot be shown, because the position doesn't match the 
   range supported by any server.
   ASSUMES: the satellites have no positional overlap. 
  */
  var satellite_name = function(latitude, longitude){
    var result = '', server, servers, limits; 
    servers = wms_servers(); 
    for (server in servers){
      if (servers.hasOwnProperty(server)){
        limits = servers[server].lat_long_limits;
        if (is_in_range(latitude, limits.φ) && is_in_range(longitude, limits.λ)){
          result = server;
          break;
        }
      }
    }
    return result;
  };
  
  var is_position_supported = function(latitude, longitude){
    return satellite_name(latitude, longitude);
  };
  
  /* Can return an empty array. Each obj in the array has .longitude, .latitude. */  
  var parse_locations = function(input){
    var result = []; 
    var locations = input.locations ? input.locations : null;
    if (locations){
      var all_locations = locations.split(";"); //separate the locations
      var i;
      for (i = 0; i < all_locations.length; i++){
        var lat_long = all_locations[i].split(","); //separate lat from long
        var position = {
          latitude: parseFloat(lat_long[0]),
          longitude: parseFloat(lat_long[1])
        };
        result.push(position);
      }
      //console.log('Found this many specific positions: ' + result.length);
    }
    return result;
  };
  
  /* Collects all the input params, and that data which is simply-derived from the input. */  
  var image_parameters = function(input, current_solar_alt_degs, canvas_id){
    var size_pixels = parseInt(input.pixels_on_a_side);
    var size_degrees = parseFloat(input.degrees_on_a_side);
    var latitude = parseFloat(input.latitude);
    var longitude = parseFloat(input.longitude);
    var delta = size_degrees/2;
    return {
      satellite_name: satellite_name(latitude, longitude),
      canvas_id: canvas_id,
      size_pixels: size_pixels,
      size_degrees: size_degrees,
      latitude: latitude,
      longitude: longitude,
      bbox : calc_bounding_box(size_degrees, latitude, longitude),
      /* The NW corner corresponds to the origin of the drawing canvas. */
      nw_corner : {longitude: longitude-delta, latitude: latitude + delta},
      layer: layer_id(input.layer, current_solar_alt_degs),
      locations: parse_locations(input)
    };
  };

  /* 
   The parts that vary from one WMS server to the next.
     preamble: the start of the server's URL
     layer_fn: input: the result of image_parameters(); output: a string to be added to the server's URL
     lat_long_limits: an object with this data: {λ.min, λ.max, φ.min, φ.max} 
  */  
  var make_wms_server = function(preamble, layer_fn, lat_long_limits){
    return {
      preamble: preamble,
      layer_fn: layer_fn, 
      lat_long_limits: lat_long_limits 
    };
  };
  var goes_layer = function(params){
    var result = '';
    var longitude_of_switchover = -101.38; //CAREFUL! the switchover is really diagonal: http://www.goes.noaa.gov/goes-w.html
    var half = params.longitude < longitude_of_switchover ? 'west' : 'east';
    if ('visible' === params.layer){
      result = half +'_vis.map&LAYERS=' + half + '_vis_1km'; 
    }
    else if ('ir' === params.layer){
      result = half + '_ir.map&LAYERS=' + half + '_ir_4km_gray'; 
    }
    return result;
  };
  var meteosat_layer = function(params){
    var result = '&LAYERS=meteosat:msg_';
    return ('visible' === params.layer) ?  result + 'vis006' : result + 'ir108'; 
  };
  var make_limits = function(λ_min, λ_max, φ_min, φ_max){
    var result = {};
    result.λ = {min:λ_min, max:λ_max};
    result.φ = {min:φ_min, max:φ_max};
    return result;
  };
  var wms_servers = function(){
    var result = {};
    result.goes = make_wms_server(
      'http://mesonet.agron.iastate.edu/cgi-bin/mapserv/mapserv?map=/mesonet/www/apps/iemwebsite/data/wms/goes/',
       goes_layer,
       make_limits(-140,-50,0,70)
    )
    result.meteosat = make_wms_server(
      'http://eumetview.eumetsat.int/geoserv/wms?',
      meteosat_layer,
      make_limits(-56,65,-70,70)
    );
    return result;
  };
  var wms_server_from = function(satellite_name){
    return wms_servers()[satellite_name];
  }
  
  /*
   GOES example: 
    ir: http://mesonet.agron.iastate.edu/cgi-bin/mapserv/mapserv?map=/mesonet/www/apps/iemwebsite/data/wms/goes/east_ir.map&LAYERS=east_ir_4km_gray&BBOX=-75,40,-73,41&WIDTH=100&HEIGHT=100&REQUEST=GetMap&SERVICE=WMS&VERSION=1.1.1&STYLES=&FORMAT=image/png&BGCOLOR=0xFFFFFF&TRANSPARENT=TRUE&SRS=EPSG:4326
   METEOSAT examples: 
    vis: http://eumetview.eumetsat.int/geoserv/wms?SERVICE=WMS&REQUEST=GetMap&SERVICE=WMS&VERSION=1.1.1&LAYERS=meteosat:msg_vis006&STYLES=&FORMAT=image/png&BGCOLOR=0xFFFFFF&TRANSPARENT=TRUE&SRS=EPSG:4326&BBOX=-11.1523918226068,30.6418455813954,0.48155275283938,40.4508576744186&WIDTH=817&HEIGHT=860
    ir: http://eumetview.eumetsat.int/geoserv/wms?SERVICE=WMS&REQUEST=GetMap&SERVICE=WMS&VERSION=1.1.1&LAYERS=meteosat:msg_ir108&STYLES=&FORMAT=image/png&BGCOLOR=0xFFFFFF&TRANSPARENT=TRUE&SRS=EPSG:4326&BBOX=-11.1523918226068,30.6418455813954,0.48155275283938,40.4508576744186&WIDTH=817&HEIGHT=860    
   The only differences between the two are the preamble, and the identifier for the layer.
  */  
  var calc_url_satellite_image = function(params){
    var wms_server = wms_server_from(params.satellite_name);
    var result = wms_server.preamble;
    result = result + wms_server.layer_fn(params);
    result = result + '&BBOX=' + params.bbox.sw_long+',' + params.bbox.sw_lat+',' + params.bbox.ne_long+',' + params.bbox.ne_lat; 
    result = result + '&WIDTH=' + params.size_pixels;
    result = result + '&HEIGHT=' + params.size_pixels;
    result = result + '&REQUEST=GetMap&SERVICE=WMS&VERSION=1.1.1&STYLES=&FORMAT=image/png&BGCOLOR=0xFFFFFF&TRANSPARENT=TRUE&SRS=EPSG:4326';
    console.log("Satellite image url: " + result);
    return result;
  };
  
  var set_canvas_size_to_match_image_size = function(canvas, params){
    canvas.height = params.size_pixels;
    canvas.width = params.size_pixels;
  };
  
  var supports_cors = function(params){
    return params.satellite_name === 'goes';
  };

  /* Showing the image in a canvas instead of an img tag lets you do a contrast stretch on the pixels. */
  var show_the_image = function(url, canvas_id, stretch_contrast, params, more_drawing_fn){
    var canvas = document.getElementById(canvas_id);
    set_canvas_size_to_match_image_size(canvas, params);
    var ctx = canvas.getContext("2d");
    var img = new Image();
    img.src = '';
    //how to respond later, after the fetch finishes
    img.onload = function() {
      console.log("Success: satellite image.");
      draw_image(this, ctx, canvas, stretch_contrast);
      draw_locations_on_satellite_image(params, canvas_id);
      if (more_drawing_fn){
        more_drawing_fn(params, canvas_id);
      }
    };
    img.onerror = function(){
      console.log("Error loading satellite image.");
    };
    //WACKY: need to add this function; otherwise the onload doesn't fire!! why??
    img.abort = function(){
      console.log("Aborted loading satellite image.");
    };
    if (supports_cors(params)) {
      img.crossOrigin = 'anonymous'; //crossOrigin (property), not crossorigin (attr)
      img.style.display = 'none';
      img.src = url; //starts the fetch of the image bytes
    }
    else {
      img.style.display = 'none';
      img.src = UTIL.crossDomainUrl(encodeURIComponent(url) + '&ext=png'); //starts the fetch of the image bytes
    }
  };

  /* Don't do anything until the image comes back from the network. */
  var draw_image = function(img, ctx, canvas, stretch_contrast){
    //console.log('The raw image has been fetched from the network.');
    ctx.drawImage(img, 0, 0);
    if (stretch_contrast){
      //console.log('Applying a linear contrast stretch to the raw image, since it is a bit dark.');
      var image_data = ctx.getImageData(0, 0, canvas.width, canvas.height); 
      var data = image_data.data;
      for (var i = 0; i < data.length; i += 4) {
        //this is an 8-bit grey-scale image 0..255; r=g=b
        //var stretched = 2*data[i];
        var stretched = contrast_stretch(data[i]);
        data[i] = stretched; // red
        data[i + 1] = stretched; // green
        data[i + 2] = stretched; // blue
      }
      //draws the image a second time, with the altered pixels    
      ctx.putImageData(image_data, 0, 0);
    }
    else {
      //console.log('Not applying a linear contrast stretch.');
    }
  };
  
  /*  Piecewise-linear contrast stretch.  */
  var contrast_stretch = function(val){
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
  
  var draw_location = function(width_frac, height_frac, ctx, params){
    GRAPH.spot(ctx, width_frac*params.size_pixels, height_frac*params.size_pixels, 2);
  };
  
  var draw_locations_on_satellite_image = function(params, canvas_id){
    var ctx = document.getElementById(canvas_id).getContext("2d");
    ctx.fillStyle = 'rgb(0,200,0)'; 
    //console.log('Drawing the location at the center of the satellite image.');
    draw_location(0.5, 0.5, ctx, params);
    for(var i = 0; i < params.locations.length; ++i){
      //console.log('Drawing location at longitude ' + params.locations[i].longitude);
      var width_frac = Math.abs(params.locations[i].longitude - params.nw_corner.longitude)/params.size_degrees;
      var height_frac = Math.abs(params.nw_corner.latitude - params.locations[i].latitude)/params.size_degrees; 
      draw_location(width_frac, height_frac, ctx, params); 
    }
  };
  
  /*
   Show a satellite image on a canvas, and allow for its customization.
   Always stretch the contrast: otherwise it's too low.
   Matches the size of the canvas to the size of the image. 
   This impl uses exclusively Web Mapping Service (WMS) servers.
   The input object has:
     .layer (visible|ir)
     .latitude
     .longitude
     .pixels_on_a_side
     .degrees_on_a_side
   The caller can customize the image, if desired, by passing in the last param, which is a function 
   that will draw on top of the returned satellite image. 
     more_drawing_fn(params, canvas_id) 
  */
  var show_image = function(input, current_solar_alt_degs, canvas_id, more_drawing_fn){
    var params = image_parameters(input, current_solar_alt_degs, canvas_id);
    var url_satellite_image = calc_url_satellite_image(params);
    show_the_image(url_satellite_image, canvas_id, STRETCH_CONTRAST_ON, params, more_drawing_fn);
  };
  
  //END OF PRIVATE ITEMS
  
  /* Return the object that contains the items needed by the caller.  */
  
  return {
    is_position_supported: is_position_supported,
    show_image: show_image
  };

}()); // the top level function is invoked here; its return value is stored in SATELLITE, a global variable