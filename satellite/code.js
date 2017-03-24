/* 
 Show an image of large-scale cloud cover, over a wide area, using a satellite image.
 Always stretch the contrast: otherwise it's too low.
 The 'input' object has these props:
   .pixels_on_a_side
   .degrees_on_a_side
   .latitude
   .longitude
   .layer (visible|ir|auto_detect)
   .locations
*/
var show_large_satellite_image = function(input, output_image){

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
   
  var increment_view_link = function(incr_id, incr_param_name, incr_amount){
    //example, replace '...latitude=10&.. with '...latitude=15&...'
    var anchor = document.getElementById(incr_id);
    var current_url = document.location.href;
    var idx1 = current_url.indexOf(incr_param_name + '=');
    var idx2 = current_url.indexOf('&', idx1);
    var old_val = current_url.substring(idx1 + incr_param_name.length + 1, idx2); // can be negative, '-78'
    var new_val = 0 + parseFloat(old_val) + incr_amount; //need to coerce to math, not string concat
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
    increment_view_link('incr_in', 'degrees_on_a_side', -4);
    increment_view_link('incr_out', 'degrees_on_a_side', 4);
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
        //console.log('Displacement ' + displacement.x + ' ' + displacement.y);
        var new_center = {}; //calc the new center coords latitude, longitude, and stay in range
        var degs_per_pixel = parseFloat(input.degrees_on_a_side) / SATELLITE.pixels_on_a_side(input.pixels_on_a_side); 
        new_center.latitude = parseFloat(input.latitude) + displacement.y * degs_per_pixel; 
        new_center.longitude = parseFloat(input.longitude) - displacement.x * degs_per_pixel;
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
  
  var use_arrow_keys = function(){
      var PAGE_UP = 33, PAGE_DOWN = 34;
      var LEFT_ARROW = 37, UP_ARROW = 38, RIGHT_ARROW = 39, DOWN_ARROW = 40;
      document.onkeydown = function(e) {
        switch (e.keyCode) {
          case PAGE_UP: 
              click_on('incr_in');
              break;
          case PAGE_DOWN: 
              click_on('incr_out');
              break;
          case LEFT_ARROW:
              click_on('incr_w');
              break;
          case UP_ARROW:
              click_on('incr_n');
              break;
          case RIGHT_ARROW:
              click_on('incr_e');
              break;
          case DOWN_ARROW:
              click_on('incr_s');
              break;
        }
    };
  };
  
  var do_all = function(){
    increment_view_links();
    use_scroll_wheel();
    use_drag_and_drop_to_recenter('clouds');
    use_arrow_keys();    
    
    var sun = EPH.position('sun', EPH.when_now(), opts(input));
    SATELLITE.show_image(input, sun.a, 'clouds');
  };
  
  return do_all();  
};