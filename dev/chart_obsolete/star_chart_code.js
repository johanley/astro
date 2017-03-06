/* 
 Show a simple star chart.
*/
var show = function(input, text_output){

  console.log('Testing...');
  var start = new Date().getTime();
  //EPH.testing();
  console.log('0.');
  EPH.testing();

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
  var star_positions = []; 
    
  /* Completely suppress the display of a dom node.  */
  var dont_show = function(thing){
    if (thing.style) {
      thing.style.display = 'none';
    }
    else {
      console.log('No style: ' + thing);
      thing.style = {};
      thing.style.display = 'none';
    }
  };  

  var when;
  if (input.date_time.length > 0){
    when = EPH.when(input.time_scale + ' ' + input.date_time);
  }
  else {
    when = EPH.when_now();
  }

  var opts_topocentric = {};
  opts_topocentric.where = EPH.where(
    parseFloat(input.latitude),
    parseFloat(input.longitude), 
    parseFloat(input.limiting_mag),
    true
  );
  opts_topocentric.units = 'rads';
  opts_topocentric.equinox = when; //equinox of date
  
  var add_summary_blurb_at_top = function(){
    var summary = document.getElementById('summary');
    var when_text; 
    if (input.date_time.length > 0){
      when_text = input.time_scale + ' ' + input.date_time;
    }
    else {
      when_text = when.toStringLT(WEEKDAYS);
    }
    summary.innerHTML = summary.innerHTML + NL + 
      '   <tr><td>' + input.location_name +  
      '<td>' + when_text +  
      '<td>' +  input.latitude + DEG +  
      '<td>' +  input.longitude + DEG + 
      '<td>' +  new Date().getTimezoneOffset() + 'm'
    ;  
  };
  add_summary_blurb_at_top();

  /* The whole point here is simply to manipulate the date-time of the URL. */
  var date_time_controls_init = function(){
    var form = document.getElementById('date_time_controls');
    var hidden_input, date_time_unit, num_steps, i, option, when_str, hidden_date_time, sign;
    //place original form input params as hidden elements of the form, so that they will be treated the same 
    //as the controls that actually modify the date-time in the URL
    for (var prop in input){
      if (input.hasOwnProperty(prop)){
        if (prop !== 'go' && prop !== 'num_steps' && prop !== 'date_time_unit'){
          hidden_input = document.createElement('input');
          hidden_input.type = 'hidden';
          hidden_input.name = prop;
          hidden_input.id = prop;
          if (prop === 'date_time' && ! input.date_time){
            //special handling for the date_time, since it defaults to 'now'
            //make the default date_time explicit
            if (input.time_scale === 'LT'){
              when_str = when.toStringLT();
            }
            else if (input.time_scale === 'UT'){
              when_str = when.toStringUT();
            }
            else if (input.time_scale === 'TT'){
              when_str = when.toStringTT();
            }
            hidden_input.value = when_str.substring(3);
          }
          else {
            hidden_input.value = input[prop];
          }
          form.appendChild(hidden_input);        
        }
      }
    }
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
      hidden_date_time = document.getElementById('date_time');
      date_time_unit = document.getElementById('date_time_unit');
      num_steps = document.getElementById('num_steps');
      sign = go_plus ? 1 : -1;
      //alter the date_time to some new value
      hidden_date_time.value = EPH.date_time_odometer(hidden_date_time.value, date_time_unit.value, sign * num_steps.value);
      return true; //continue with regular form processing
    };
  };
  date_time_controls_init();

  /* Allow the user to change the center of the chart. */  
  var direction_controls_init = function(){
    var form = document.getElementById('direction_controls');
    var hidden_input, num_degs, option, hidden_center_ra, hidden_center_dec, sign, op_code;
    //place original form input params as hidden elements of the form, so that they will be treated the same 
    //as the controls that actually modify the date-time in the URL
    for (var prop in input){
      if (input.hasOwnProperty(prop)){
        if (prop !== 'go' && prop !== 'num_degs'){
          hidden_input = document.createElement('input');
          hidden_input.type = 'hidden';
          hidden_input.name = prop;
          hidden_input.id = prop;
          hidden_input.value = input[prop];
          form.appendChild(hidden_input);        
        }
      }
    }
    //prepopulate controls that modify the URL by recycling req params, if present
    if (input.num_degs){
      document.getElementById('num_degs').value = input.num_degs;
    }
    form.onsubmit = function(){
      hidden_center_ra = document.getElementById('center_ra');
      hidden_center_dec = document.getElementById('center_dec');
      num_degs = parseFloat(document.getElementById('num_degs').value);
      op_code = go_op;
      sign = op_code % 2 ? 1 : -1;
      if (op_code <= 2){
        hidden_center_ra.value = parseFloat(hidden_center_ra.value)  + sign * num_degs;
      }
      else {
        hidden_center_dec.value = parseFloat(hidden_center_dec.value)  + sign * num_degs;
      }
      return true; //continue with regular form processing (a GET to the altered URL)
    };
  };
  direction_controls_init();

  var standard_canvas_appearance = function(ctx, big_font){
    if (big_font){
      myBigFont(ctx);
    }
    ctx.fillStyle = 'rgb(255,255,255)'; 
    ctx.strokeStyle = 'rgb(255,255,255)';
    ctx.textAlign='center';
    ctx.textBaseline = 'middle';
  };
    
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
  
  /* Find the pos on Earth where the object is in the zenith (for the given when). */
  var find_sub_star_point = function(thing){
    var result = {}; // a 'where' object
    result.φ = thing.δ;
    result.λ = EPH.in2pi(thing.α - when.gmst);
    //corrections for parallax don't apply here! they apply only for the 'real' position
    //the idea is that for the real position, parallax is applied AND it changes ra and dec
    //the sub-star point then uses that ra and dec to calc it's nominal (geocentric) a,A
    result.is_topocentric = false;
    return result;
  };
  
  /* Return the xy position on the canvas, ready for plotting. */
  var chart_projection = function(thing, chart, where_sub_ctr){
    var pos = {};
    //find (a,A) with respect to the where_sub_ctr
    //to avoid having data for two different places in the same object, copy the ephem data into a temp object
    var ephem_nonce = {α: thing.α, δ: thing.δ, equinox: thing.equinox};  
    EPH.convert_αδ_to_aA(ephem_nonce, where_sub_ctr, when); //geocentric, no parallax; adds (a,A,h) and leaves others alone
    //now it's the same as in the regular planisphere
    var halfpi = 0.5*Math.PI;
    var r = (halfpi - ephem_nonce.a) * chart.pixels_per_rad;
    //to rotate the planisphere about the zenith, just add an angle to theta:
    //'polar azimuthal equidistant projection': simple r-theta
    var theta = ephem_nonce.A;  
    pos.x = chart.ctr.x - r * Math.sin(theta);
    pos.y = chart.ctr.y - r * Math.cos(theta);
    return pos;
  };
  
  var find_star = function(star, chart, where_sub_ctr, ctx, i){
    var size = star_size(star[mag]);
    var pos = chart_projection({α: star[α], δ:star[δ]}, chart, where_sub_ctr);
    if (size > 0){
      remember_star_position(i, pos.x, pos.y, size);
    }
  };

  /* Remember later, so that stars and constellation lines can be plotted in the desired order. */
  var remember_star_position = function(i, x, y, size){
    star_positions[i] = {x:x, y:y, size:size}; //not an array: the index is sparse, non-contiguous
  };
  
  var plot_star = function(ctx, star_position){
    spot(ctx, star_position.x, star_position.y, star_position.size);
  };
  
  var plot_star_name = function(ctx, star_position, i){
    //drop the constellation designation
    var star_name = EPH.stars[i][name];
    //text(ctx, i, star_position.x, star_position.y+gap);
    if (star_name.length > 0){
      if (star_name.indexOf(' ') !== -1){
        star_name = star_name.split(" ")[0]; //don't use the constellation abbr
      }
      if (! /^\d+$/.test(star_name)){
        //show bayer letters
        text(ctx, star_name, star_position.x, star_position.y+gap);
      }
      else if (input.include_bayer === '1') {
        //show any other name (Flamsteed numbers don't extend to the southern pole).
        text(ctx, star_name, star_position.x, star_position.y+gap);
      }
    }
  };
  
  var find_stars = function(chart, where_sub_ctr,ctx){
    var i, star;
    var limiting_mag = parseFloat(input.limiting_mag);
    for (i=0; i < EPH.stars.length; ++i){
      star = EPH.stars[i];
      if (star[mag] <= limiting_mag){
        if (EPH.elongation_between(chart.ctr, {α: star[α], δ: star[δ]}) <= chart.width + EPH.rads(MARGIN_DEGS)){
          find_star(star, chart, where_sub_ctr, ctx, i);
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
  
  var plot_messier = function(thing, chart, where_sub_ctr, ctx){
    var pos = chart_projection({α: thing[α], δ:thing[δ]}, chart, where_sub_ctr);
    var tiny_gap = 2;
    point(ctx, pos.x+tiny_gap, pos.y+tiny_gap);
    point(ctx, pos.x+tiny_gap, pos.y-tiny_gap);
    point(ctx, pos.x-tiny_gap, pos.y+tiny_gap);
    point(ctx, pos.x-tiny_gap, pos.y-tiny_gap);
    ctx.save();
    ctx.font = '10px Verdana';
    text(ctx, thing[name], pos.x, pos.y+gap);
    ctx.restore();
  };
  
  var show_messiers = function(chart, where_sub_ctr, ctx){
    var i, messier;
    var limiting_mag = parseFloat(input.limiting_mag_messiers);
    for (i=0; i < EPH.messiers.length; ++i){
      messier = EPH.messiers[i];
      if (messier[mag] <= limiting_mag){
        if (EPH.elongation_between(chart.ctr, {α: messier[α], δ: messier[δ]}) <= chart.width){
            plot_messier(messier, chart, where_sub_ctr, ctx);
        }
      }
    }
  };
  
  var plot_planet = function(planet, chart, where_sub_ctr, ctx){
    var size = star_size(planet.mag);
    var pos = chart_projection(planet, chart, where_sub_ctr);
    if (size > 0){
      spot(ctx, pos.x, pos.y, size);
      text(ctx, planet.symbol, pos.x, pos.y+gap);
    }
  };
  
  var show_planets_etc = function(chart, where_sub_ctr, ctx){
      var planets = EPH.as_array(EPH.planets);
      var planet;
      for (var i=0; i < planets.length; ++i){
        planet = EPH.position(planets[i].name, when, opts_topocentric);
        planet.symbol = planets[i].symbol; //add it here for the sake of convenience
        if (EPH.elongation_between(chart.ctr, planet) <= chart.width){
            plot_planet(planet, chart, where_sub_ctr, ctx);
        }
      }
  };
  
  var show_star_chart = function(){
    var chart = {}, ctr = {}, ephem;
    ctr.α = EPH.rads(parseFloat(input.center_ra));
    ctr.δ = EPH.rads(parseFloat(input.center_dec));
    if (input.center_object){
      //override the ra-dec input, if present
      ephem = EPH.position(input.center_object, when, opts_topocentric); //fails if the name is unknown
      ctr.α = ephem.α;
      ctr.δ = ephem.δ;
    }
    var limiting_mag = parseFloat(input.limiting_mag);
    chart.width = EPH.rads(parseInt(input.chart_width));
    var canvas = document.getElementById('star_chart');
    canvas.width = window.innerWidth*0.94;
    canvas.height = window.innerHeight*.85;
    var ctx = canvas.getContext('2d');
    ctx.fillStyle = 'rgb(36,55,114)';
    ctx.fillRect(0, 0, canvas.width, canvas.height);    
    standard_canvas_appearance(ctx, true);
    h = canvas.height;
    w = canvas.width;
    ctr.x = w/2;
    ctr.y = h/2;
    chart.ctr = ctr;
    chart.w = w;
    chart.h = h;
    chart.pixels_per_rad = chart.w / chart.width;
    var where_sub_ctr = find_sub_star_point(chart.ctr);
    find_stars(chart, where_sub_ctr, ctx);
    connect_the_star_positions(ctx);
    show_stars(ctx);
    show_messiers(chart, where_sub_ctr, ctx);
    show_planets_etc(chart, where_sub_ctr, ctx);
  };
  
  /* The top-level call. */
  show_star_chart();

  var end = new Date().getTime();
  console.log('Done. msecs: ' + (end-start));

  // moon with correct phase, size, and position angle
  // add other, occasional items: minor planets, comets, meteor showers
  // useful: tooltip to describe object x; store object coords and details in memory; when within x pixels, show its data
  //     what if 2 things are very close? pick the closest one to the pointer
  // colors choices
};