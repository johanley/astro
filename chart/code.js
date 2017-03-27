/* 
 Show a simple star chart.
 If a location is present, then topocentric parallax is calculated for the Moon (and no other object).
*/
var show = function(input){

  //constants
  var NL = '\r\n'; // new line
  var DEG = '&deg;';
  var WEEKDAYS = ['Sun', 'Mon', 'Tue', 'Wed', 'Thu', 'Fri', 'Sat'];
  //indices into the star-messier data structure (array of items)
  var mag = 3;  
  var α = 1  
  var δ = 2;
  var name = 0;
  var gap = 14; //between object and associated text 
  var MARGIN_DEGS = 10; //marginal stars need positions, in order to plot constellation lines that go off-canvas

  var when;
  if (input.time_scale && input.date_time && input.date_time.length > 0){
    when = EPH.when(input.time_scale + ' ' + input.date_time);
  }
  else {
    when = EPH.when_now();
  }

  /* The 'when' as input by the user in the date-time control, and as it appears in the URL. */  
  var when_as_input_str = function(){
    var result = '';
    if (input.time_scale === 'LT'){
      result = when.toStringLT();
    }
    else if (input.time_scale === 'UT'){
      result = when.toStringUT();
    }
    else if (input.time_scale === 'TT'){
      result = when.toStringTT();
    }
    return result.substring(3);
  };
  
  var default_setting = function(name, val){
    if (typeof input[name] === "undefined" || input[name].length === 0){
      input[name] = val;
    }
  };
  var apply_default_settings_if_absent = function(){
    default_setting('place_at_top', 'north');
    default_setting('date_time', when_as_input_str()); 
    default_setting('time_scale', 'LT');
    default_setting('location_name', '');
    default_setting('latitude', '0');
    default_setting('longitude', '0');
    default_setting('limiting_mag', '8.0');
    default_setting('chart_width', '30');
    //default_setting('show_identifiers', ['bayer']);  
    default_setting('style', 'regular');
  };
  apply_default_settings_if_absent();
  
  var star_positions = [];

  var opts_topocentric = {}; //the Moon needs a topocentric position, because of its parallax
  opts_topocentric.where = EPH.where(
    parseFloat(input.latitude),
    parseFloat(input.longitude), 
    parseFloat(input.limiting_mag),
    true
  );
  opts_topocentric.units = 'rads';
  opts_topocentric.equinox = when; //equinox of date

  /* 
   The two center-on styles, ra-decl versus object, are somewhat at odds with each other.
   When centering on an object, we need to calculate the initial center-on ra and decl, in order to play nice with controls
   that tweak the position of the center by a few degrees. Once the tweak-center controls are activated, the user has 
   implicitly (and irreversibly) moved from a center-on-object style to a center-on-ra-decl style.
   SIDE-EFFECT: can override input.center_ra and input.center_dec.  
  */
  var initialize_center = function(){
    var result = {}, object_found = false;
    if (input.center_object){
      ephem = EPH.position(input.center_object, when, opts_topocentric); //fails if the name is unknown
      if (ephem){
        object_found = true;
        if (input.center_object === 'Moon'){
          //needed only because the chart projection relies on ra-dec, not-alt-az
          EPH.apply_parallax_to_αδ(ephem, opts_topocentric.where);          
        }
        result.α = ephem.α;
        result.δ = ephem.δ;
        input.center_ra = EPH.degs(result.α); //side-effect: keep the data 100% consistent; needed for the tweak-center controls
        input.center_dec = EPH.degs(result.δ);
      }
    }
    if (! object_found) {
      result.α = EPH.rads(parseFloat(input.center_ra));
      result.δ = EPH.rads(parseFloat(input.center_dec));
    }
    return result;
  }
  var ctr = initialize_center();
  
  var add_summary_blurb = function(){
    var summary = document.getElementById('summary');
    summary.innerHTML = summary.innerHTML + NL + 
      '   <tr><td>' + UTIL.escapeHtml(input.location_name) +  
      '<td>' + UTIL.escapeHtml(input.time_scale) + ' ' + UTIL.escapeHtml(input.date_time) +   
      '<td>' +  UTIL.escapeHtml(input.latitude) + DEG +  
      '<td>' +  UTIL.escapeHtml(input.longitude) + DEG + 
      '<td>' +  new Date().getTimezoneOffset() + 'm'
    ;  
  };
  add_summary_blurb();

  /* Needed since there's more than 1 form, whose items share names. Return the dom node under the form, whose name matches. */
  var form_child_element_by_name = function(form /*dom node*/, control_name){
    var result, i;
    for(i=0; i < form.children.length; ++i){
      if (form.children[i].name === control_name){
        result = form.children[i];
        break;
      }
    }
    return result;
  };
  
  var array_includes = function(item, items){
    var result = false;
    for(var i = 0; i < items.length; ++i){
      if (items[i] === item){
        result = true;
        break;
      }
    }
    return result;
  };
  
  /* Tweaker-controls/properties are used to tweak the URL. */
  var is_tweaker_param = function(item, tweakers){
    return array_includes(item, tweakers);
  };

  var is_array = function(value){ //Crockford p61
    return value && typeof value === 'object' && value.constructor === Array;
  };  
  var add_hidden_item = function(name, form, value){
    var hidden_input = document.createElement('input');
    hidden_input.type = 'hidden';
    hidden_input.name = name; //use name, not id! there are N such forms on the page, whose items share these names
    hidden_input.value = value; //input[prop];
    form.appendChild(hidden_input);        
  };
  var put_hidden_controls_into = function(form, tweakers){
    for (var prop in input){
      if (input.hasOwnProperty(prop)){
        if (!is_tweaker_param(prop, tweakers)){
          if (is_array(input[prop])){ //for a multivalued request param
            for(var i = 0; i < input[prop].length; ++i){
              add_hidden_item(prop, form, input[prop][i]);
            }
          }
          else {
            add_hidden_item(prop, form, input[prop]);
          }
        }
      }
    }
  };
  
  /*  Tweak the date-time of the URL. */
  var date_time_controls_init = function(){
    var form = document.getElementById('date_time_controls');
    var hidden_input, date_time_unit, num_steps, i, option, when_str, hidden_date_time, sign;
    var tweakers = ['go', 'num_steps', 'date_time_unit'];
    put_hidden_controls_into(form, tweakers);
    //prepopulate date-time controls that modify the URL by recycling req params, if present
    if (input.num_steps){
      document.getElementById('num_steps').value = input.num_steps;
    }
    if (input.date_time_unit){
       date_time_unit = document.getElementById('date_time_unit');
       for(i=0; i < date_time_unit.options.length ; ++i){
         option = date_time_unit.options[i];
         if (option.value === input.date_time_unit){
           option.selected = true;
         }
       }
    }
    form.onsubmit = function(){
      hidden_date_time = form_child_element_by_name(form, 'date_time');
      date_time_unit = document.getElementById('date_time_unit');
      num_steps = document.getElementById('num_steps');
      sign = go_plus ? 1 : -1;
      //alter the hidden date_time control to some new value
      hidden_date_time.value = EPH.date_time_odometer(hidden_date_time.value, date_time_unit.value, sign * num_steps.value);
      return true; //continue with regular form processing
    };
  };
  date_time_controls_init();

  /* Stay in range 0..360. */
  var ra_degs_constraint = function(α){
    var result = α % 360; // range -360..360
    if (result < 0) result = result + 360; // range 0..360
    return result;
  };
  /* Stop at the poles; stay in range -90..90. */
  var dec_degs_constraint = function(δ){
    var result  = Math.min(δ, 90); //stop at the poles; stay in range -90..90
    result = Math.max(result, -90); 
    return result;
  };

  /* 
   Tweak the center of the chart. 
   If the user initially centered on an object, then that is discarded by this action, and it reverts to centering on the object's ra and decl.
   That means that the center-on-object case needs to remember the center as ra, decl, and store those numbers in the forms. 
  */  
  var direction_controls_init = function(){
    var form = document.getElementById('direction_controls');
    var num_degs, hidden_center_ra, hidden_center_dec, sign, op_code, num_degs_input, new_dec, new_ra;
    var tweakers = ['go', 'num_degs'];
    put_hidden_controls_into(form, tweakers);
    //if present, center_object overrides center-on-ra-dec; this form negates center_object, and reverts to center-ra-dec, 
    //so we want to make sure center_object is empty
    form_child_element_by_name(form, 'center_object').value = '';
    //prepopulate controls that modify the URL by recycling req params, if present
    if (input.num_degs){
      document.getElementById('num_degs').value = input.num_degs;
    }
    form.onsubmit = function(){
      hidden_center_ra = form_child_element_by_name(form, 'center_ra');
      hidden_center_dec = form_child_element_by_name(form, 'center_dec');
      num_degs_input = document.getElementById('num_degs').value;
      num_degs = num_degs_input ? parseFloat(num_degs_input) : 0.15 * parseFloat(input.chart_width);
      op_code = go_op;
      sign = op_code % 2 ? 1 : -1;
      if (op_code <= 2){
        new_ra = parseFloat(hidden_center_ra.value)  + sign * num_degs;
        new_ra = ra_degs_constraint(new_ra);
        hidden_center_ra.value = new_ra;
      }
      else {
        new_dec = parseFloat(hidden_center_dec.value)  + sign * num_degs;
        new_dec = dec_degs_constraint(new_dec);
        hidden_center_dec.value = new_dec;
      }
      return true; //continue with regular form processing (a GET to the altered URL)
    };
  };
  direction_controls_init();
  
  var zoom_controls_init = function(){
    var form = document.getElementById('zoom_controls');
    var tweakers = ['go'];
    put_hidden_controls_into(form, tweakers);
    form.onsubmit = function(){
      var val, FACTOR=1.25, MIN=1, MAX=90, ABANDON=false, SUBMIT_FORM=true;
      var hidden_chart_width = form_child_element_by_name(form, 'chart_width');
      if (zoom_go_op === 1){
        val = parseFloat(hidden_chart_width.value) * FACTOR;
        val = Math.min(MAX, val);
      }
      else {
        val = parseFloat(hidden_chart_width.value) / FACTOR;
        val = Math.max(MIN, val);
      }
      if (MIN <= val && val <= MAX){
        hidden_chart_width.value = val;
        return SUBMIT_FORM; //continue with regular form processing (a GET to the altered URL)
      }
      return ABANDON; 
    };
  };
  zoom_controls_init();
  
  var click_on = function(id){
    var item = document.getElementById(id);
    item.click();
  };
  /* Change the center of the chart using the arrow keys, in addition to a click. */
  var use_special_keys = function(){
      document.onkeydown = function(e) {
        switch (e.keyCode) {
          //arrow keys for change of center
          case 37:
              click_on('incr_w');
              break;
          case 38:
              click_on('incr_n');
              break;
          case 39:
              click_on('incr_e');
              break;
          case 40:
              click_on('incr_s');
              break;
          //page up-down for zoom in-out
          case 33:
              click_on('zoom_plus');
              break;
          case 34:
              click_on('zoom_minus');
              break;
        }
    };
  };
  use_special_keys();
  
  var use_scroll_wheel = function(){
    //http://stackoverflow.com/questions/14926366/mousewheel-event-in-modern-browsers
    var change_zoom = function(e){
      // cross-browser wheel delta
      var e = window.event || e; // old IE support
      var delta = Math.max(-1, Math.min(1, (e.wheelDelta || -e.detail)));
      if (delta > 0) {
        click_on('zoom_plus');
      }
      else {
       click_on('zoom_minus');
      }
    };
    document.addEventListener("wheel", change_zoom, false);
  };
  use_scroll_wheel();
  
  /* Near poles this is very sensitive to the motion, and the experience degrades a bit. */
  var use_drag_and_drop_to_recenter = function(chart){
    var canvas = document.getElementById('star_chart');
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
        var new_center = {}; //calc the new center coords ra-dec, and stay in range
        new_center.ra = parseFloat(input.center_ra) + EPH.degs(displacement.x / chart.pixels_per_rad); 
        new_center.ra = ra_degs_constraint(new_center.ra); 
        new_center.dec = parseFloat(input.center_dec) + EPH.degs(displacement.y / chart.pixels_per_rad);
        new_center.dec = dec_degs_constraint(new_center.dec);
        //somewhat dangerously, I overwrite the input object; this is reasonably safe since the page is about to be reloaded
        input.center_ra = new_center.ra;
        input.center_dec = new_center.dec;
        input.center_object = ''; // no longer centered on an object
        var new_url =  'graphic.sky' + build_href_to_form([]);
        document.location = new_url;
      }
    };
    canvas.onmousedown = mouse_down;
    canvas.onmouseup = mouse_up;
    canvas.onmouseout = mouse_up;
  };

  var append_item_to_href = function(href, name, value){
    return href + name + '=' + encodeURIComponent(value) + '&';
  };
  var build_href_to_form = function(tweakers){
    var href = '?';
    for (var prop in input){
      if (input.hasOwnProperty(prop)){
        if (!is_tweaker_param(prop, tweakers)){
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
    }
    return href;
  };
  /* 
   Link back to the main form, such that current settings are preserved.
   Tweaker params are ignored. Multivalued params have special handling. 
  */
  var link_back_to_main_form = function(tweakers){
    var link = document.getElementById('link_back');
    link.href = link.href + build_href_to_form(tweakers);
  };
  link_back_to_main_form(['go', 'num_degs', 'num_steps', 'date_time_unit']);
  
  var standard_canvas_appearance = function(ctx, use_big_font){
    if (use_big_font){
      GRAPH.myBigFont(ctx);
    }
    ctx.textAlign = 'center';
    ctx.textBaseline = 'middle';
  };
  
  var canvas_colors = function(is_regular, ctx, canvas){
    if (is_regular){ //screen-friendly, white on deep blue
      ctx.fillStyle = 'rgb(36,55,114)'; 
      ctx.fillRect(0, 0, canvas.width, canvas.height);    
      ctx.fillStyle = 'white';
      ctx.strokeStyle = 'white';
    }
    else { //print-friendly, black on white
      ctx.fillStyle = 'black';
      ctx.strokeStyle = 'black';
    }
  };

  //SAME    
  var star_size = function(mag){
    var A0 = 3.5; //the max 
    var result;
    var limiting_mag = parseFloat(input.limiting_mag);
    if (mag > limiting_mag){
      result = -1; //won't be drawn; error if even attempted
    }
    else if (mag <= 1.25){
      result = A0; //capped
    }
    else if (mag <= 5.0){
      result = A0*Math.exp(-mag/2.5);
    }
    else {
      result = 0.1;
    }
    return result;
  };
  
  /* Find the pos on Earth where the object is in the zenith (for the given when; geocentric). */
  var find_sub_star_point = function(thing){
    var result = {}; // a 'where' object
    result.φ = thing.δ;
    result.λ = EPH.in2pi(thing.α - when.gmst);
    //corrections for parallax don't apply here! they apply only for the 'real' position.
    //the idea is that for the real position, parallax is already applied before this code is reached, 
    //AND that the parallax is applied to alter ra and dec (not just a, A).
    //the sub-star point then uses that ra and dec to calc it's (non-topocentric) a,A.
    result.is_topocentric = false;
    return result;
  };
  
  /*
   Project from a given center. Use elongation from the center, and the polar angle to the object (wrt the center).
   If the center is near a celestial pole, then use the pole as the center, along with a plain 'rθ projection'.
  */
  var chart_projection = function(ephem, chart /*.ctr.α,ctr.δ in rads*/){
    var result;
      if (Math.abs(ctr.δ) < EPH.rads(88)){
      result = project_from_center(ephem, chart); 
    }
    else {
      result = project_from_pole(ephem, chart);
    }
    return result;
  };
  var project_from_center = function(ephem, chart/*.ctr.α,ctr.δ in rads*/){
    var pos = {};
    var elong = EPH.elongation_between(chart.ctr, ephem);
    var θ = EPH.position_angle_between(chart.ctr, ephem); //NCP-ctr-ephem
    if (input.place_at_top === 'south') {
      θ = θ + Math.PI;
    }
    var r = elong * chart.pixels_per_rad;
    pos.x = chart.ctr.x - r * Math.sin(θ);
    pos.y = chart.ctr.y - r * Math.cos(θ);
    return pos;
  };
  var project_from_pole = function(ephem, chart/*.ctr.α,ctr.δ in rads*/){
    var pos = {}, elong;
    var halfpi = 0.5*Math.PI;
    var θ = ephem.α; 
    if (chart.ctr.δ > 0){ //north celestial pole
      elong = halfpi - ephem.δ;
    }
    else { //south celestial pole
      elong = halfpi + ephem.δ;
    }
    var r = elong * chart.pixels_per_rad;
    pos.x = chart.ctr.x + r * Math.sin(θ);
    pos.y = chart.ctr.y - r * Math.cos(θ);
    return pos;
  };
  
  /* Remember later, so that stars and constellation lines can be plotted in the desired order. */
  var remember_star_position = function(i, x, y, size){
    star_positions[i] = {x:x, y:y, size:size}; //not an array: the index is sparse, non-contiguous
  };
  
  var find_star = function(star, chart, ctx, i){
    var size = star_size(star[mag]);
    var pos = chart_projection({α: star[α], δ:star[δ]}, chart);
    if (size > 0){
      remember_star_position(i, pos.x, pos.y, size);
    }
  };

  var plot_star = function(ctx, star_position){
    GRAPH.spot(ctx, star_position.x, star_position.y, star_position.size);
  };
  
  var plot_star_name = function(ctx, star_position, i){
    var star_name = EPH.stars[i][name];
    if (star_name.length > 0){
      if (star_name.indexOf(' ') !== -1){
        star_name = star_name.split(" ")[0]; //chop off the constellation abbr
      }
      var show_star_name = false;
      if (/^\d+$/.test(star_name)){ //numbers only 
        if (array_includes('flamsteed', input.show_identifiers)){
          show_star_name = true;
        }
      }
      else {
        if (array_includes('bayer', input.show_identifiers)){
          show_star_name = true;
        }
      }
      if (show_star_name){
        GRAPH.text(ctx, star_name, star_position.x, star_position.y + gap); 
      }
    }
  };
  
  var find_stars = function(chart, ctx){
    var i, star;
    var limiting_mag = parseFloat(input.limiting_mag);
    for (i=0; i < EPH.stars.length; ++i){
      star = EPH.stars[i];
      if (star[mag] <= limiting_mag){
        if (EPH.elongation_between(chart.ctr, {α: star[α], δ: star[δ]}) <= chart.width + EPH.rads(MARGIN_DEGS)){
          find_star(star, chart, ctx, i);
        }
      }
    }
  };
  
  var draw_polyline = function(ctx, points){
    var i;
    ctx.save();
    ctx.strokeStyle='rgb(125,125,125)';
    //ctx.strokeStyle='rgb(177,156,217)';
    ctx.beginPath();
    ctx.moveTo(points[0].x,points[0].y);
    for (i=1;i<points.length;i++){
      ctx.lineTo(points[i].x,points[i].y);
    }
    ctx.stroke();
    ctx.closePath();
    ctx.restore();
  };
  
  /* Join up sets of star positions, to show the usual lines denoting a constellation. */
  var connect_the_star_positions = function(ctx){
    var polylines, i, count, polyline, j, points;
    //for all possible polylines
    for (constell in EPH.constellation_lines){
      if (EPH.constellation_lines.hasOwnProperty(constell)){
        polylines = EPH.constellation_lines[constell]; //array of arrays
        for(i=0; i<polylines.length; ++i){
          polyline = polylines[i];
          points = [];
          //find the points that have a corresponding calculated position
          for(j=0; j<polyline.length; ++j){
            if (star_positions[polyline[j]]){
              points.push(star_positions[polyline[j]]);
            }
          }
          if (points.length >= 2){
            //lines(ctx, points);
            draw_polyline(ctx, points);
          }
        }
      }
    }
  };

  /* Do this after the constellation lines, such that the spots overwrite the lines. */  
  var show_stars = function(ctx){
    var i;
    for (i in star_positions){
      if (star_positions.hasOwnProperty(i)){
        plot_star(ctx, star_positions[i]);
        plot_star_name(ctx, star_positions[i], i);
      }
    }
  };
  
  var plot_messier = function(thing, chart, ctx){
    var pos = chart_projection({α: thing[α], δ:thing[δ]}, chart);
    var tiny_gap = 2;
    GRAPH.point(ctx, pos.x+tiny_gap, pos.y+tiny_gap);
    GRAPH.point(ctx, pos.x+tiny_gap, pos.y-tiny_gap);
    GRAPH.point(ctx, pos.x-tiny_gap, pos.y+tiny_gap);
    GRAPH.point(ctx, pos.x-tiny_gap, pos.y-tiny_gap);
    ctx.save();
    ctx.font = '10px Verdana';
    GRAPH.text(ctx, thing[name], pos.x, pos.y+gap);
    ctx.restore();
  };
  
  var show_deep_sky_objects = function(list, chart, ctx){
    var i, dso;
    for (i=0; i < list.length; ++i){
      dso = list[i];
      if (EPH.elongation_between(chart.ctr, {α: dso[α], δ: dso[δ]}) <= chart.width){
          plot_messier(dso, chart, ctx);
      }
    }
  };
  
  var show_comets = function(chart, ctx){
     var comets = EPH.as_array(EPH.comets);
     for(var i=0; i < comets.length; ++i){
       var ephem = EPH.position(comets[i].name, when, opts_topocentric);
        if (EPH.elongation_between(chart.ctr, ephem) <= chart.width){
          var pos = chart_projection(ephem, chart);
          GRAPH.circle(ctx,pos.x, pos.y,3);
          GRAPH.text(ctx, comets[i].alt_name, pos.x, pos.y+gap);
        }
     }
  };
  
  var show_minor_planets = function(chart, ctx){
     var items = EPH.as_array(EPH.minor_planets);
     for(var i=0; i < items.length; ++i){
       var ephem = EPH.position(items[i].name, when, opts_topocentric);
        if (EPH.elongation_between(chart.ctr, ephem) <= chart.width){
          var pos = chart_projection(ephem, chart);
          GRAPH.square(ctx, pos.x, pos.y,6);
          GRAPH.text(ctx, items[i].name, pos.x, pos.y+gap);
        }
     }
  };
  
  /* Use brightness by default, but apparent angular size IF that would result in a bigger dot. */
  var planet_size = function(planet, chart){
    var result = star_size(planet.mag);
    var size_as_fraction_of_chart_width = planet.size / chart.width; // both in rads
    var size_in_pixels = Math.floor(size_as_fraction_of_chart_width * chart.w); // pixels
    if (size_in_pixels > result){
      result = size_in_pixels;
    }
    return result;
  };
  
  var plot_planet = function(planet, chart, ctx){
    var size = planet_size(planet, chart);
    //override the case in which the planet is dimmer than the input limiting stellar magnitude, to always show the planet
    if (size === -1){ 
      size = 0.1;
    }
    var pos = chart_projection(planet, chart);
    GRAPH.spot(ctx, pos.x, pos.y, size);
    GRAPH.text(ctx, planet.symbol, pos.x, pos.y+size+gap);
  };
  
  var show_planets_etc = function(chart, ctx){
      var planets = EPH.as_array(EPH.planets);
      var planet;
      for (var i=0; i < planets.length; ++i){
        planet = EPH.position(planets[i].name, when, opts_topocentric);
        if (planets[i].name === 'Moon'){
          //needed only because the chart projection relies on ra-dec, not-alt-az
          EPH.apply_parallax_to_αδ(planet, opts_topocentric.where);          
        }
        planet.symbol = planets[i].symbol; //add it here for the sake of convenience
        if (EPH.elongation_between(chart.ctr, planet) <= chart.width){
            plot_planet(planet, chart, ctx);
        }
      }
  };
  
  var show_star_chart = function(){
    var chart = {}, ephem;
    chart.width = EPH.rads(parseInt(input.chart_width));
    var canvas = document.getElementById('star_chart');
    canvas.width = window.innerWidth*0.95;
    canvas.height = window.innerHeight*.85;
    var ctx = canvas.getContext('2d');
    //ctx.fillStyle = 'rgb(36,55,114)'; //blue background
    //ctx.fillRect(0, 0, canvas.width, canvas.height);    
    standard_canvas_appearance(ctx, true);
    canvas_colors(input.style === 'regular', ctx, canvas);    
    h = canvas.height;
    w = canvas.width;
    ctr.x = w/2;
    ctr.y = h/2;
    chart.ctr = ctr;
    chart.ctr.α = EPH.rads(input.center_ra);
    chart.ctr.δ = EPH.rads(input.center_dec);
    chart.w = w;
    chart.h = h;
    chart.pixels_per_rad = chart.w / chart.width;
    find_stars(chart, ctx);
    connect_the_star_positions(ctx);
    show_stars(ctx);
    show_deep_sky_objects(EPH.messiers, chart, ctx);
    show_deep_sky_objects(EPH.caldwells, chart, ctx);
    show_planets_etc(chart, ctx);
    show_comets(chart, ctx);
    show_minor_planets(chart, ctx);
    use_drag_and_drop_to_recenter(chart);
  };
  
  /* The top-level call. */
  show_star_chart();

  // animation, over time or space
  //   with optional trailing path, tick marks ever N days, single object
  //   nice for planets, comets, asteroids, moon
  //   printable version for the whole year
  //
  // moon with correct phase, and position angle wrt NCP
  //    this needed for better charting of zoomed-in cases, especially for occultations
  //
  // projection appropriate for horizon views
  //    calc a, A
  //    show the horizon as a straight horizontal line
  //    restrict up-down navigation in some way (max altitude)?  or better: just cut off showing alt < 0 items
  //    probably replace center-on ra, dec with center an alt-az
  //    probably add an option for the projection: regular, near-horizon only
  //    allow zoom
  //
  // center-object: 
  //    constellation names?
  //       the abbr is already in ephem.js 
  //       in that case, I would need to input the ra-decl of each constellation manually somewhere; maybe directly in ephem.js
  //       the IMCEE French data may be useful here
  // 
  // colors control (html color picker not supported universally)
  // colors for printing: need to re-draw the canvas, with different colors!
  // very useful: pop-up window to describe object x; store object coords and details in memory; when within x pixels, show its data
  //   what if 2 things are very close? pick the closest one to the pointer; or maybe show N, not just 1 
};