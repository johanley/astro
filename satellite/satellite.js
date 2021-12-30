/* 
 Show satellite images. 
 Current impl talks only to WMS servers.
 The impl chooses a server according to the input lat and long.
 The channel (visible or IR) depends on the current altitude of the Sun at the given location (fancy!).
 
 Annoying differences between WMS 1.1 and 1.3:
   - the order of lat-long as it appears in the BBOX parameter changed between WMS 1.1 and WMS 1.3.
   - the name of the projection request param changes from SRS to CRS

 As usual, servers have a habit of changing, which makes support a nuisance. 
 Currently, I'm using a single server, and supporting only North America. 
  - Real Earth SSEC
     2019-12-02: unstable?
     'msLoadMap(): Unable to access file. (/home/wms/data/mapfiles/G16-ABI-FD-BAND13.map)'
     2020-04-16: seeing again, 2020-04-16 Thursday.
       same error as above; also:
       'error on line 1 at column 1: Extra content at the end of the document' 
  - RETIRED: NOAA Nowcoast
      https://nowcoast.noaa.gov/help/#!section=mapservices
      https://nowcoast.noaa.gov/help/#!section=map-service-list
  - RETIRED: Iowa State IEM
     https://mesonet.agron.iastate.edu/
  - RETIRED: UK server, Meteosat
   
*/
var SATELLITE = (function(){ 

  var STRETCH_CONTRAST_ON = true;  

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
  
  /* Return a simple string id for the layer, (visible|ir), according to the current altitude of the Sun. */
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
   Returns the name of the WMS server that fits with the given lat-long.
   If returns an empty string, then an image cannot be shown, because the position doesn't match the 
   range supported by any server.
   ASSUMES: the servers have no positional overlap. 
  */
  var server_name = function(latitude, longitude){
    var result = '', server, servers, limits; 
    servers = wms_servers(); 
    for (server in servers){
      if (servers.hasOwnProperty(server)){
        limits = servers[server].lat_long_limits;
        if (is_in_range(latitude, limits.φ) && is_in_range(longitude, limits.λ)){
          result = server; //the name of the property
          break;
        }
      }
    }
    return result;
  };

  /* Returns an empty string if no server found. Empty strings are falsey. */  
  var is_position_supported = function(latitude, longitude){
    return server_name(latitude, longitude);
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
  
  var pixels_on_a_side_default = function(){
    var width = window.innerWidth;
    var height = window.innerHeight;
    var smallest_dimension = Math.min(width, height);
    return Math.floor(smallest_dimension * 0.90);
  };
  
  var pixels_on_a_side = function(val){
    var result = pixels_on_a_side_default();
    if (val) {
      result = parseInt(val);
    }
    return result;  
  };
  
  /* Collects all the input params, and that data which is simply-derived from the input. */  
  var image_parameters = function(input, current_solar_alt_degs, canvas_id){
    var size_pixels = pixels_on_a_side(input.pixels_on_a_side);
    var size_degrees = parseFloat(input.degrees_on_a_side);
    var latitude = parseFloat(input.latitude);
    var longitude = parseFloat(input.longitude);
    var delta = size_degrees/2;
    return {
      server_name: server_name(latitude, longitude),
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
     wms_version: version of the WMS protocol 
  */  
  var make_wms_server = function(preamble, layer_fn, lat_long_limits, wms_version){
    return {
      preamble: preamble,
      layer_fn: layer_fn, 
      lat_long_limits: lat_long_limits,
      wms_version: wms_version 
    };
  };
  //RETIRED
  var iowa_state_goes_layer = function(params){
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
  //RETIRED
  var nowcoast_layer = function(params){
    var result = '&layers=';
    return ('visible' === params.layer) ?  result + '9' : result + '17'; 
  };
  //RETIRED
  var meteosat_layer = function(params){
    var result = '&LAYERS=meteosat:msg_';
    return ('visible' === params.layer) ?  result + 'vis006' : result + 'ir108'; 
  };
  var real_earth_ssec_layer = function(params){
    //Goes-16, Full Disk, red-visible and IR 
    //http://realearth.ssec.wisc.edu/products/G16-ABI-FD-BAND02  - RED 
    //http://realearth.ssec.wisc.edu/products/G16-ABI-FD-TC  - TRUE COLO
    //http://realearth.ssec.wisc.edu/products/G16-ABI-FD-BAND13  - IR
    return ('visible' === params.layer) ?  '?map=G16-ABI-FD-TC.map&LAYERS=latest' : '?map=G16-ABI-FD-BAND13.map&LAYERS=latest&LAYERS=latest';
  };
  var make_limits = function(λ_min, λ_max, φ_min, φ_max){
    var result = {};
    result.λ = {min:λ_min, max:λ_max};
    result.φ = {min:φ_min, max:φ_max};
    return result;
  };
  var wms_servers = function(){
    var result = {};
    //CURRENTLY USED FOR ALL OF NORTH AMERICA
    result.real_earth_ssec = make_wms_server(
      'http://realearth.ssec.wisc.edu/cgi-bin/mapserv',
      real_earth_ssec_layer,
      make_limits(-140,-50,0,70),
      '1.3.0'
    );
    //NO LONGER USED:
    /*
    result.nowcoast = make_wms_server(
      'https://nowcoast.noaa.gov/arcgis/services/nowcoast/sat_meteo_imagery_time/MapServer/WmsServer?',
      nowcoast_layer,
      make_limits(-140,-50,0,70),
      '1.3.0'
    );
    //currently restricted to the West (higher quality than Nowcoast) 
    result.iowa_state = make_wms_server(
      'http://mesonet.agron.iastate.edu/cgi-bin/mapserv/mapserv?map=/mesonet/www/apps/iemwebsite/data/wms/goes/',
       iowa_state_goes_layer,
       make_limits(-140,-101,0,70),
       '1.1.1'
    );
    result.meteosat = make_wms_server(
      'http://eumetview.eumetsat.int/geoserv/wms?',
      meteosat_layer,
      make_limits(-56,65,-70,70),
      '1.1.1'
    );
    */
    return result;
  };
  var wms_server_from = function(server_name){
    return wms_servers()[server_name];
  }
  
  /*
   Real Earth SSEC examples:
    vis (TRUE COLOR): http://realearth.ssec.wisc.edu/cgi-bin/mapserv?map=G16-ABI-FD-TC.map&REQUEST=GetMap&SERVICE=WMS&VERSION=1.3.0&LAYERS=latest&STYLES=&FORMAT=image/png&BGCOLOR=0xFFFFFF&TRANSPARENT=TRUE&CRS=EPSG:4326&BBOX=35.7656742955095,-117.006018627837,45.8201672081869,-107.454250360794&WIDTH=817&HEIGHT=860
    vis (RED):        http://realearth.ssec.wisc.edu/cgi-bin/mapserv?map=G16-ABI-FD-BAND02.map&REQUEST=GetMap&SERVICE=WMS&VERSION=1.3.0&LAYERS=latest&STYLES=&FORMAT=image/png&BGCOLOR=0xFFFFFF&TRANSPARENT=TRUE&CRS=EPSG:4326&BBOX=35.7656742955095,-117.006018627837,45.8201672081869,-107.454250360794&WIDTH=817&HEIGHT=860
     ir:              http://realearth.ssec.wisc.edu/cgi-bin/mapserv?map=G16-ABI-FD-BAND13.map&LAYERS=latest&LAYERS=latest&VERSION=1.3.0&REQUEST=GetMap&SERVICE=WMS&STYLES=&FORMAT=image/png&WIDTH=480&HEIGHT=480&BBOX=44.58,-66.28,47.58,-63.28&CRS=EPSG
   Iowa State example (GOES): NO LONGER USED 
    ir: http://mesonet.agron.iastate.edu/cgi-bin/mapserv/mapserv?map=/mesonet/www/apps/iemwebsite/data/wms/goes/west_ir.map&LAYERS=west_ir_4km_gray&BBOX=-114.35,48.19,-111.35,51.19&WIDTH=480&HEIGHT=480&REQUEST=GetMap&SERVICE=WMS&VERSION=1.1.1&STYLES=&FORMAT=image/png&BGCOLOR=0xFFFFFF&TRANSPARENT=TRUE&SRS=EPSG
   Nowcoast example (GOES East and West): NO LONGER USED
    ir:  https://nowcoast.noaa.gov/arcgis/services/nowcoast/sat_meteo_imagery_time/MapServer/WmsServer?request=GetMap&format=image/png&version=1.3.0&service=WMS&width=512&height=512&crs=EPSG:4326&bbox=44.76,-64.64,47.76,-61.64&layers=17&styles=
    vis: https://nowcoast.noaa.gov/arcgis/services/nowcoast/sat_meteo_imagery_time/MapServer/WmsServer?request=GetMap&format=image/png&version=1.3.0&service=WMS&width=512&height=512&crs=EPSG:4326&bbox=44.76,-64.64,47.76,-61.64&layers=9&styles=
   METEOSAT examples: NO LONGER USED
    vis: http://eumetview.eumetsat.int/geoserv/wms?SERVICE=WMS&REQUEST=GetMap&SERVICE=WMS&VERSION=1.1.1&LAYERS=meteosat:msg_vis006&STYLES=&FORMAT=image/png&BGCOLOR=0xFFFFFF&TRANSPARENT=TRUE&SRS=EPSG:4326&BBOX=-11.1523918226068,30.6418455813954,0.48155275283938,40.4508576744186&WIDTH=817&HEIGHT=860
    ir: http://eumetview.eumetsat.int/geoserv/wms?SERVICE=WMS&REQUEST=GetMap&SERVICE=WMS&VERSION=1.1.1&LAYERS=meteosat:msg_ir108&STYLES=&FORMAT=image/png&BGCOLOR=0xFFFFFF&TRANSPARENT=TRUE&SRS=EPSG:4326&BBOX=-11.1523918226068,30.6418455813954,0.48155275283938,40.4508576744186&WIDTH=817&HEIGHT=860    
   
   Note two annoying differences between WMS 1.1 and WMS 1.3: 
     - the order of lat-long in the bounding box
     - the change from param name SRS to CRS
  */  
  var calc_url_satellite_image = function(params){
    var wms_server = wms_server_from(params.server_name);
    var result = wms_server.preamble;
    //WARNING: for Iowa State, the layer leaks into the base URL; it needs to come first here: 
    result = result + wms_server.layer_fn(params); 
    result = result + '&VERSION=' + wms_server.wms_version;
    result = result + '&REQUEST=GetMap&SERVICE=WMS&STYLES=&FORMAT=image/png';
    result = result + '&WIDTH=' + params.size_pixels;
    result = result + '&HEIGHT=' + params.size_pixels;
    var PROJECTION = 'EPSG:4326';
    if (wms_server.wms_version === '1.1.1'){
      result = result + '&BGCOLOR=0xFFFFFF&TRANSPARENT=TRUE';
      result = result + '&BBOX=' + params.bbox.sw_long+',' + params.bbox.sw_lat+',' + params.bbox.ne_long+',' + params.bbox.ne_lat; 
      result = result + '&SRS=' + PROJECTION;
    }
    else {
      //result = result + '&BGCOLOR=0xFFFFFF&TRANSPARENT=TRUE';
      result = result + '&BBOX=' + params.bbox.sw_lat+',' + params.bbox.sw_long+',' + params.bbox.ne_lat+',' + params.bbox.ne_long; 
      result = result + '&CRS=' + PROJECTION;
    }
    console.log("Satellite image url: " + result);
    return result;
  };
  
  var set_canvas_size_to_match_image_size = function(canvas, params){
    canvas.height = params.size_pixels;
    canvas.width = params.size_pixels;
  };
  
  var supports_cors = function(params){
    var result = 
      params.server_name === 'real_earth_ssec' ||  
      params.server_name === 'iowa_state' || 
      params.server_name === 'nowcoast'
    ;
    return result;
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
    img.onerror = function(evt){
      console.log("Error loading satellite image." + img.src); //the evt object has no useful error data; sad!
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
      console.log('Applying a linear contrast stretch to the raw image, since it is a bit dark.');
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
      console.log('Not applying a linear contrast stretch.');
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
  
  var draw_center = function(ctx, params){
    var size = 4;
    GRAPH.tickMarkVertical(ctx, 0.5*params.size_pixels, 0.5*params.size_pixels, size);
    GRAPH.tickMarkHorizontal(ctx, 0.5*params.size_pixels, 0.5*params.size_pixels, size);
  };
  
  var draw_location = function(width_frac, height_frac, ctx, params){
    GRAPH.spot(ctx, width_frac*params.size_pixels, height_frac*params.size_pixels, 2);
  };
  
  var draw_with = function(color, ctx){
    ctx.fillStyle = color; 
    ctx.strokeStyle = color; 
  };
  
  var draw_locations_on_satellite_image = function(params, canvas_id){
    var ctx = document.getElementById(canvas_id).getContext("2d");
    var black = 'rgb(0,0,0)';
    var green = 'rgb(0,200,0)';
    draw_with(black, ctx);
    //console.log('Drawing the location at the center of the satellite image.');
    draw_center(ctx, params);
    draw_with(green, ctx);
    for(var i = 0; i < params.locations.length; ++i){
      //console.log('Drawing location at longitude ' + params.locations[i].longitude);
      var width_frac = (params.locations[i].longitude - params.nw_corner.longitude)/params.size_degrees;
      var height_frac = (params.nw_corner.latitude - params.locations[i].latitude)/params.size_degrees;
      if (width_frac >= 0 && height_frac >= 0){
        draw_location(width_frac, height_frac, ctx, params); 
      }
    }
  };
  
  /** Only done for visible, not IR. */
  var stretch_contrast = function(params){
    return params.layer === 'visible' ? true : false;
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
    //show_the_image(url_satellite_image, canvas_id, STRETCH_CONTRAST_ON, params, more_drawing_fn);
    show_the_image(url_satellite_image, canvas_id, stretch_contrast(params), params, more_drawing_fn);
  };
  
  //END OF PRIVATE ITEMS
  
  /* Return the object that contains the items needed by the caller.  */
  
  return {
    is_position_supported: is_position_supported,
    show_image: show_image,
    pixels_on_a_side: pixels_on_a_side 
  };

}()); // the top level function is invoked here; its return value is stored in SATELLITE, a global variable