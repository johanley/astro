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
  var gap = 12; //between object and associated text
    
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
  
  var opts = {};
  opts.where = EPH.where(
    parseFloat(input.latitude),
    parseFloat(input.longitude), 
    parseFloat(input.limiting_mag)
  );
  opts.units = 'rads';
  opts.equinox = when; //equinox of date
  
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
    //much better to have a continuous size, not discrete chunking
    //var result = -(4/5)*mag+4; //linear function of magnitude; invisible above mag 5;
    var limiting_mag = parseFloat(input.limiting_mag);
    var max = 3.5;
    var result = -(max/limiting_mag)*mag+max; //linear function of magnitude; invisible above mag 5;
    return (result > max ? max : result); //cap the max size
  };
  
  /* Find the pos on Earth where the object is in the zenith (for the given when). */
  var find_sub_star_point = function(thing){
    var result = {};
    result.φ = thing.δ;
    result.λ = EPH.in2pi(thing.α - when.gmst);
    return result;
  };
  /* Return the xy position on the canvas, ready for plotting. */
  var chart_projection = function(α, δ, chart){
    var pos = {};
    pos.x = chart.ctr.x + (chart.ctr.α - α) * chart.dec_factor * chart.pixels_per_rad;
    pos.y = chart.ctr.y + (chart.ctr.δ - δ) * chart.pixels_per_rad;
    return pos;
  };
  var plot_star = function(star, chart, ctx){
    var size = star_size(star[mag]);
    var pos = chart_projection(star[α], star[δ], chart);
    if (size > 0){
      spot(ctx, pos.x, pos.y, size);
    }
    //drop the constellation designation
    var star_name = star[name];
    if (star_name.length > 0){
      if (star_name.indexOf(' ') !== -1){
        star_name = star[name].split(" ")[0]; //don't use the constellation abbr
      }
      text(ctx, star_name, pos.x, pos.y+gap);
    }
  };
  var show_stars = function(chart, box, ctx){
    var i, star;
    var limiting_mag = parseFloat(input.limiting_mag);
    for (i=0; i < EPH.stars.length; ++i){
      star = EPH.stars[i];
      if (star[mag] <= limiting_mag){
        if (star[α] > box.min_α && star[α] < box.max_α &&  
            star[δ] > box.min_δ && star[δ] < box.max_δ){
            plot_star(star, chart, ctx);
        }
      }
    }
  };
  var plot_messier = function(thing, chart, ctx){
    var pos = chart_projection(thing[α], thing[δ], chart);
    var tiny_gap = 2;
    point(ctx, pos.x+tiny_gap, pos.y+tiny_gap);
    point(ctx, pos.x+tiny_gap, pos.y-tiny_gap);
    point(ctx, pos.x-tiny_gap, pos.y+tiny_gap);
    point(ctx, pos.x-tiny_gap, pos.y-tiny_gap);
    text(ctx, thing[name], pos.x, pos.y+gap);
  };
  var show_messiers = function(chart, box, ctx){
    var i, messier;
    var limiting_mag = parseFloat(input.limiting_mag_messiers);
    for (i=0; i < EPH.messiers.length; ++i){
      messier = EPH.messiers[i];
      if (messier[mag] <= limiting_mag){
        if (messier[α] > box.min_α && messier[α] < box.max_α &&  
            messier[δ] > box.min_δ && messier[δ] < box.max_δ){
            plot_messier(messier, chart, ctx);
        }
      }
    }
  };
  var plot_planet = function(planet, chart, ctx){
    var size = star_size(planet.mag);
    var pos = chart_projection(planet.α, planet.δ, chart);
    spot(ctx, pos.x, pos.y, size);
    text(ctx, planet.symbol, pos.x, pos.y+gap);
  };
  var show_planets_etc = function(chart, box, ctx){
      var planets = EPH.as_array(EPH.planets);
      var planet;
      for (var i=0; i < planets.length; ++i){
        planet = EPH.position(planets[i].name, when, opts);
        if (planet.α > box.min_α && planet.α < box.max_α &&  
            planet.δ > box.min_δ && planet.δ < box.max_δ){
            planet.symbol = planets[i].symbol;
            plot_planet(planet, chart, ctx);
        }
      }
  };
  /* 
   Rectangular star chart, with declination not too close to the celestial poles.
  */
  var show_star_chart = function(){
    var i, star;
    var chart = {}, ctr = {}, box = {};
    ctr.α = EPH.rads(parseFloat(input.center_ra));
    ctr.δ = EPH.rads(parseFloat(input.center_dec));
    var limiting_mag = parseFloat(input.limiting_mag);
    var chart_width = EPH.rads(parseInt(input.chart_width));
    var dec_factor = Math.cos(ctr.δ); //plate-carrée projection, center on a circle of declination
    box.min_α = ctr.α - chart_width / dec_factor; //ra-range is wider than dec-range
    box.max_α = ctr.α + chart_width / dec_factor;
    box.min_δ = ctr.δ - chart_width;
    box.max_δ = ctr.δ + chart_width;
    var canvas = document.getElementById('star_chart');
    var ctx = canvas.getContext('2d');
    standard_canvas_appearance(ctx, true);
    h = canvas.height;
    w = canvas.width;
    ctr.x = w/2;
    ctr.y = h/2;
    chart.ctr = ctr;
    chart.w = w;
    chart.h = h;
    chart.dec_factor = dec_factor;
    chart.pixels_per_rad = chart.w / chart_width;
    show_stars(chart, box, ctx);
    show_messiers(chart, box, ctx);
    show_planets_etc(chart, box, ctx);
  };
  show_star_chart();
  
  
  var end = new Date().getTime();
  console.log('Done. msecs: ' + (end-start));


  
  
  
  // PLANISPHERE..........
  
  
  
  
  var planet_rows = function(include_phenomena){
      var planets = EPH.as_array(EPH.planets);
      var rows = [];
      var phen, i, ephem, is_moon, row;
      for (i=0; i < planets.length; ++i){
        ephem = EPH.position(planets[i].name, when, opts);
        row = {thing: planets[i], ephem: ephem};
        if (include_phenomena){
          is_moon = planets[i].name === 'Moon';
          if (planets[i].name === 'Sun'){
            phen = EPH.rise_culmination_set_daily(when, opts.where, planets[i], is_moon);
          }
          else {
            phen = EPH.rise_culmination_set_observation_window(when, opts.where, planets[i], is_moon);
          }
          row.phen = phen;
        }
        rows.push(row);
      }
      return rows;
  };
  /* Simple projection. */
  var plot_on_planisphere = function(ctx, astre, mag, center, R, label, type){
    var r = (90 - astre.ephem.a) * (R/90);
    //to rotate the planisphere about the zenith, just add an angle to theta:
    //'polar azimuthal equidistant projection': simple r-theta
    var theta = EPH.rads(astre.ephem.A);  
    var x = center.x - r * Math.sin(theta);
    var y = center.y - r * Math.cos(theta);
    var size = star_size(mag);
    var gap = 2;
    if (type === 'messier'){
        point(ctx, x+gap, y+gap);
        point(ctx, x+gap, y-gap);
        point(ctx, x-gap, y+gap);
        point(ctx, x-gap, y-gap);
    }
    else if (size > 0){
      spot(ctx, x, y, star_size(mag));
      if (label){
        text(ctx, label, x, y+15);
      }
    }
  };
  var show_planisphere = function(visible_stars, planets_etc){
    var h, w, i, star, R, center, x, y, r, theta, planets, planet, messier;
    var limiting_mag = parseFloat(input.limiting_mag);
    var canvas = document.getElementById('planisphere');
    var ctx = canvas.getContext('2d');
    standard_canvas_appearance(ctx, true);
    h = canvas.height;
    w = canvas.width;
    R = w/2 - 4; // radius of the planisphere
    center = {x: w/2, y:h/2};
    circle(ctx,center.x,center.y,R); //horizon
    for(i=0; i < visible_stars.length; ++i){
      star = visible_stars[i]; //.thing and .ephem
      plot_on_planisphere(ctx, star, star.thing[3], center, R);
    }
    planets = planets_etc !== null ? planets_etc : planet_rows(false);
    for(i=0; i < planets.length; ++i){
      planet = planets[i];
      if (planet.ephem.a > 0 && planet.ephem.mag <= limiting_mag){
        plot_on_planisphere(ctx, planet, planet.ephem.mag, center, R, planet.thing.symbol);
        //console.log('Name:' + planet.thing.name + ' alt:' + planet.ephem.a + ' A:' + planet.ephem.A);
      }
    }
  };
  var visible_stars = EPH.find_visible_stars(when, opts.where, 0, parseFloat(input.limiting_mag), opts);
  //show_planisphere(visible_stars, planet_rows());
  
  
  // tool tip for ra, dec
  // add minor planets
  // add comets
  // add meteor shower radiants
  // center on a specific object, instead of ra, dec; look up its pos, and use its ra, dec as the center
  // fix projection problem near 0h
  // fix ybs data: pi-6 Ori, etc has null 
  // add projection near the celestial poles, for use when the decl in large
  // useful: tooltip to describe object x; store object coords and details in memory; when within x pixels, show its data
  //     what if 2 things are very close? pick the closest one to the pointer
  // drag the canvas updates the display; use to move to an adjacent region
  //     http://rectangleworld.com/blog/archives/129
};