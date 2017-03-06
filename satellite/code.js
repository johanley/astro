/* 
 Show an image of large-scale cloud cover, over a wide area, using a GOES satellite image.
 Always stretch the contrast: otherwise it's too low.
 No customization showing 'green dots' for specific locations.
 The 'input' object has these props;
   .pixels_on_a_side
   .degrees_on_a_side
   .latitude
   .longitude
   .layer (visible|ir|auto_detect)
*/
var show_large_scale_clouds = function(input, output_image){

  var opts = function(input){
    var result = {};
    result.where = EPH.where(
      parseFloat(input.latitude),
      parseFloat(input.longitude), 
      0,
      true
    );
    result.units = 'degs';
    result.equinox = EPH.when_now(); //equinox of date
    return result;
  };
   
  var sun = EPH.position('sun', EPH.when_now(), opts(input));
  var current_solar_alt_degs = sun.a;

  /* these are all required */
  var size_pixels = parseInt(input.pixels_on_a_side);
  var size_degrees = parseFloat(input.degrees_on_a_side);
  var latitude = parseFloat(input.latitude);
  var longitude = parseFloat(input.longitude);
  
  var STRETCH_CONTRAST_ON = true;  
  var STRETCH_CONTRAST_OFF = true;

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
  
  var calcUrlClouds = function(bbox){
    var url = 'http://mesonet.agron.iastate.edu/cgi-bin/mapserv/mapserv?map=/mesonet/www/apps/iemwebsite/data/wms/goes/';
    var longitude_of_switchover = -101.38; //CAREFUL! the switchover is really diagonal: http://www.goes.noaa.gov/goes-w.html
    var half = longitude < longitude_of_switchover ? 'west' : 'east';
    
    if ('auto_detect' === input.layer){
      //overwrite the value with the 'real' one
      if (current_solar_alt_degs > 0){ //there's ~15m delay in getting images from space
        input.layer = 'visible';
      }
      else {
        input.layer = 'ir';
      }
    }
    
    if ('visible' === input.layer){
      url = url + half +'_vis.map&LAYERS=' + half + '_vis_1km'; 
    }
    else if ('ir' === input.layer){
      url = url + half + '_ir.map&LAYERS=' + half + '_ir_4km_gray'; 
    }
    
    url = url + '&BBOX=' + bbox.sw_long+',' + bbox.sw_lat+',' + bbox.ne_long+',' + bbox.ne_lat; 
    url = url + '&WIDTH=' + size_pixels;
    url = url + '&HEIGHT=' + size_pixels;
    url = url + '&REQUEST=GetMap&SERVICE=WMS&VERSION=1.1.1&STYLES=&FORMAT=image/png&BGCOLOR=0xFFFFFF&TRANSPARENT=TRUE&SRS=EPSG:4326';
    return url;    
  };
  
  var setCanvasSizeToMatchImageSize = function(canvas){
    canvas.height = size_pixels;
    canvas.width = size_pixels;
  };
  
  /** 
   Showing the image in a canvas instead of an img tag let's you do a contrast stretch on the pixels.
   The server returning the image needs to support CORS, in this impl. 
  */
  var showTheImage = function(url, output, stretch_contrast){
    var canvas = document.getElementById(output);
    setCanvasSizeToMatchImageSize(canvas);
    var ctx = canvas.getContext("2d");
    
    //canvas.width = window.innerWidth * 0.50;
    //canvas.width = 480;
    //canvas.height = canvas.width;
    
    var img = new Image();
    img.src = '';
    //how to respond later, after the fetch finishes
    img.onload = function() {
      drawImage(this, ctx, canvas, stretch_contrast);
    };
    img.onerror = function(){
      console.log("Error loading cloud image.");
    };
    //wacky: need to add this function; otherwise the onload doesn't fire!! why??
    img.abort = function(){
      console.log("Aborted loading cloud image.");
    };
    //needed, since the image is on someone else's server; that server supports CORS
    img.crossOrigin = 'anonymous'; //crossOrigin (property), not crossorigin (attr)
    img.style.display = 'none';
    img.src = url; //starts the fetch of the image bytes 
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

  /** The NW corner corresponds to the origin of the drawing canvas. */
  var find_north_west_corner = function(){
    var delta = size_degrees/2;
    return {
      longitude: longitude - delta,
      latitude: latitude + delta
    };
  };
  
  var increment_view_link = function(incr_id, incr_param_name, incr_amount){
    //example, replace '...latitude=10&.. with '...latitude=15&...'
    var anchor = document.getElementById(incr_id);
    var current_url = document.location.href;
    var idx1 = current_url.indexOf(incr_param_name + '=');
    var idx2 = current_url.indexOf('&', idx1);
    var old_val = current_url.substring(idx1 + incr_param_name.length + 1, idx2); // can be negative, '-78'
    var new_val = parseFloat(old_val) + incr_amount; //need to coerce to math, not string concat
    var new_url_start = current_url.substring(0, idx1 + incr_param_name.length);
    var new_url_end = current_url.substring(idx2);
    var new_url = new_url_start + '=' + new_val + new_url_end;
    anchor.href = new_url;
  };

  /* links to change the lat/long by N degs in the cardinal directions. */  
  var increment_view_links = function(){
    increment_view_link('incr_e', 'longitude', 5);
    increment_view_link('incr_w', 'longitude', -5);
    increment_view_link('incr_n', 'latitude', 5);
    increment_view_link('incr_s', 'latitude', -5);
    increment_view_link('incr_in', 'degrees_on_a_side', -2);
    increment_view_link('incr_out', 'degrees_on_a_side', 2);
  };
  
  var click_on = function(id){
    var item = document.getElementById(id);
    item.click();
  };
  var use_scroll_wheel = function(){
    //http://stackoverflow.com/questions/14926366/mousewheel-event-in-modern-browsers
    var change_zoom = function(e){
      // cross-browser wheel delta
      var e = window.event || e; // old IE support
      var delta = Math.max(-1, Math.min(1, (e.wheelDelta || -e.detail)));
      if (delta > 0) {
        click_on('incr_in');
      }
      else {
       click_on('incr_out');
      }
    };
    document.addEventListener("wheel", change_zoom, false);
  };
  
  var is_array = function(value){ //Crockford p61
    return value && typeof value === 'object' && value.constructor === Array;
  };  
  /* Stay in range 0..60. */
  var latitude_constraint = function(val){
	var result = val;
	if (val < 0) {result = 0;}
	if (val > 60) {result = 60;}
    return result;
  };
  /* Stay in range -50...-140 */
  var longitude_constraint = function(val){
	var result = val;
	if (val < -140) {result = -140;}
	if (val > -50) {result = -50;}
    return result;
  };
  var append_item_to_href = function(href, name, value){
    return href + name + '=' + encodeURIComponent(value) + '&';
  };
  var build_href_to_form = function(tweakers){
    var href = '?';
    for (var prop in input){
      if (input.hasOwnProperty(prop)){
	      if (is_array(input[prop])){ //for a multivalued request param
	        for(var i = 0; i < input[prop].length; ++i){
	          href = append_item_to_href(href, prop, input[prop][i]);
	        }
	      }
	      else {
	        href = append_item_to_href(href, prop, input[prop]);
	      }
      }
    }
    return href;
  };
  /* Near the poles this is very sensitive to the motion, and the experience degrades a bit. */
  var use_drag_and_drop_to_recenter = function(canvas_id){
    var canvas = document.getElementById(canvas_id);
    var start = {};
    var is_down = false;
    var current_pos = function(e){
      var pos = {};
      pos.x = parseInt(e.clientX - canvas.offsetLeft);
      pos.y = parseInt(e.clientY - canvas.offsetTop);
      return pos;
    };
    var mouse_down = function(e){
      e.preventDefault();
      if (!is_down) {
        is_down = true;
        start = current_pos(e);
      }
    };
    var mouse_up = function(e){
      if (is_down){
        is_down = false;
        var stop = current_pos(e);
        var displacement = {};
        displacement.x = stop.x - start.x;
        displacement.y = stop.y - start.y;
        console.log('Displacement ' + displacement.x + ' ' + displacement.y);
        var new_center = {}; //calc the new center coords latitude, longitude, and stay in range
        var degs_per_pixel = parseFloat(input.degrees_on_a_side) / parseInt(input.pixels_on_a_side); 
        new_center.latitude = parseFloat(input.latitude) + displacement.y * degs_per_pixel; 
        new_center.latitude = latitude_constraint(new_center.latitude); 
        new_center.longitude = parseFloat(input.longitude) - displacement.x * degs_per_pixel;
        new_center.longitude = longitude_constraint(new_center.longitude);
        //somewhat dangerously, I overwrite the input object; this is reasonably safe since the page is about to be reloaded
        input.latitude = new_center.latitude;
        input.longitude = new_center.longitude;
        var new_url =  'graphic.sky' + build_href_to_form([]);
        document.location = new_url;
      }
    };
    canvas.onmousedown = mouse_down;
    canvas.onmouseup = mouse_up;
    canvas.onmouseout = mouse_up;
  };
  
  
  //New:
  //http://mesonet.agron.iastate.edu/cgi-bin/mapserv/mapserv?map=/mesonet/www/apps/iemwebsite/data/wms/goes/east_ir.map&LAYERS=east_ir_4km_gray&BBOX=-77.16,43.9,-74.16,46.9&WIDTH=480&HEIGHT=480&REQUEST=GetMap&SERVICE=WMS&VERSION=1.1.1&STYLES=&FORMAT=image/png&BGCOLOR=0xFFFFFF&TRANSPARENT=TRUE&SRS=EPSG
  //Orig
  //http://mesonet.agron.iastate.edu/cgi-bin/mapserv/mapserv?map=/mesonet/www/apps/iemwebsite/data/wms/goes/east_ir.map&LAYERS=east_ir_4km_gray&BBOX=-77.16,43.9,-74.16,46.9&WIDTH=480&HEIGHT=480&REQUEST=GetMap&SERVICE=WMS&VERSION=1.1.1&STYLES=&FORMAT=image/png&BGCOLOR=0xFFFFFF&TRANSPARENT=TRUE&SRS=EPSG
  
  /** Return the bbox, since it's needed by the NASA Worldview. */
  var doAll = function(){
    increment_view_links();
    use_scroll_wheel();
    var bbox = calcBoundingBox();
    var cloudsUrl = calcUrlClouds(bbox);
    console.log('Clouds URL: ' + cloudsUrl);
    //console.log('Showing the clouds image in a canvas.');
    showTheImage(cloudsUrl, output_image, STRETCH_CONTRAST_ON);
    use_drag_and_drop_to_recenter('clouds');    
  };
  
  return doAll();  
};