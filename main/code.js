/* 
 Show numerous ephemerides for the given input date.
*/
var show = function(input, MET_OFFICE_API_KEY, is_dev, lang){

  var start = new Date().getTime();

  var NL = '\r\n'; // new line
  var DEG = '&deg;';
  var compass_points = ['N', 'NNE', 'NE', 'ENE', 'E', 'ESE', 'SE', 'SSE', 'S', 'SSW', 'SW', 'WSW', 'W', 'WNW', 'NW', 'NNW', 'N'];
  var WEEKDAYS = ['Su', 'Mo', 'Tu', 'We', 'Th', 'Fr', 'Sa'];
  var WEEKDAYS_LONG = ['Sunday', 'Monday', 'Tuesday', 'Wednesday', 'Thursday', 'Friday', 'Saturday'];
  var RISE = 0;
  var SET = 1;
  
  var default_setting = function(name, val){
    if (typeof input[name] === "undefined"){
      input[name] = val;
    }
  };
  var apply_default_settings_if_absent = function(){
    default_setting('date_time', ''); //now
    default_setting('time_scale', 'LT');
    default_setting('limiting_mag', '5.1');
    default_setting('limiting_mag_messiers', '11.0');
    default_setting('limiting_mag_messiers_planisphere', '8.0');
    default_setting('degrees_on_a_side', '3');
    default_setting('pixels_on_a_side', '480');
    default_setting('layer', 'auto_detect');
    default_setting('twilight', '-12');
    default_setting('planisphere_rotation_angle', '0');
    default_setting('aurora_min_activity_level', '-1');
    default_setting('occultations_num_days_ahead', '10');
    default_setting('occultations_min_mag', '6');
  };
  apply_default_settings_if_absent();
     
  var highest_first = function(row_a, row_b){
    return row_a.ephem.a > row_b.ephem.a ? -1 : 1;
  };
  
  var nearest_minute = function(date_time_str /* 2016-10-31 15:33:35.0 */){
    var result = date_time_str.substring(0,16);
    var seconds = date_time_str.substring(17);
    if (seconds.startsWith('3') || seconds.startsWith('4') || seconds.startsWith('5')){
      result = EPH.date_time_odometer(result, 'min', 1);
    }
    return result.substring(0,16); 
  };

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

  /* First index of desired phenomenon, if present; -1 if not found. */
  var phen_idx = function(phen, type, name){
    var result = -1;
    for (var i = 0; i < phen[type].length; ++i){
      if (phen[type][i].name === name){
        result = i;
        break;
      }
    } 
    return result;
  };
  
  /* NOT BEING USED FOR THE MOMENT - USING ELONGATION INSTEAD.
  var jd_of_earliest = function(target, events){
    var result = Number.MAX_VALUE;
    for (var i = 0; i < events.length; ++i){
      if (events[i].name === target && events[i].when.jd < result){
        result = events[i].when.jd;
      }
    }
    return result < Number.MAX_VALUE ? result : null;
  };  
  var set_transit_or_rise = function(a){
    var result = jd_of_earliest('set', a.phen.horizons);
    if (!result){
      result = jd_of_earliest('upper', a.phen.culminations);
    }
    if (!result){
      result = jd_of_earliest('rise', a.phen.horizons);
    }
    return result;
  };
  var by_time_of_set_transit_rise = function(a, b){
    var a_time = set_transit_or_rise(a);
    var b_time = set_transit_or_rise(b);
    return a_time - b_time;
  };
  */
  var elong_from_sun = function(a){
    var result = a.ephem.elong ? a.ephem.elong : 0; // the Sun itself has no elong property
    if (result < 0){ //elong has a discontinuity at pi, where it jumps to negative numbers
      result = result + 360; //0..360, assumes degs
    }
    return result;
  };
  var by_elongation_from_the_sun = function(a, b){
    return elong_from_sun(a) - elong_from_sun(b);
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
    parseFloat(input.limiting_mag),
    true
  );
  opts.units = 'degs';
  opts.equinox = when; //equinox of date
  
  var opts_rads = {};
  opts_rads.where = EPH.where(
    parseFloat(input.latitude),
    parseFloat(input.longitude), 
    parseFloat(input.limiting_mag),
    true
  );
  opts_rads.equinox = when; //equinox of date

  //this impl suffices for 2 langs (list); it's a simple matter to change it to handle N langs instead (map) 
  var all_text = [];
  var populate_translations = function(){
    all_text['e'] = 'f';
    all_text['Population index; smaller means brighter'] = 'Indice de population; plus petit signifie plus brilliant';
    all_text['Factor by which ZHR is multiplied, for current alt and limiting mag'] = 'Muliplier THZ par ce facteur, pour l&apos;effet de hauteur';
    all_text['Weather'] = 'Météo';

    all_text['Sun'] = 'Soleil';
    all_text['Moon'] = 'Lune';
    all_text['Mercury'] = 'Mercure';
    all_text['Venus'] = 'Vénus';
    all_text['Mars'] = 'Mars';
    all_text['Jupiter'] = 'Jupiter';
    all_text['Saturn'] = 'Saturne';
    all_text['Uranus'] = 'Uranus';
    all_text['Neptune'] = 'Neptune';
        
    all_text['SSW'] = 'SSO';
    all_text['SW'] = 'SO';
    all_text['WSW'] = 'OSO';
    all_text['W'] = 'O';
    all_text['WNW'] = 'ONO';
    all_text['NW'] = 'NO';
    all_text['NNW'] = 'NNO';
    
    all_text['Su'] = 'Dim';
    all_text['Mo'] = 'Lun';
    all_text['Tu'] = 'Mar';
    all_text['We'] = 'Mer';
    all_text['Th'] = 'Jeu';
    all_text['Fr'] = 'Ven';
    all_text['Sa'] = 'Sam';
    
    all_text['bright'] = 'brilliante';
    all_text['fade'] = 'diminuer';
    all_text['steady'] = 'constante';
    all_text['morning'] = 'matin';
    all_text['evening'] = 'soir';
    all_text['early morning'] = 'tôt le matin';
    all_text['early evening'] = 'tôt le soir';
    
    all_text['PA'] = 'AP';
    all_text['DD'] = 'DS'; //sombre, claire
    all_text['DB'] = 'DC';
    all_text['RD'] = 'RS';
    all_text['RB'] = 'RC';
    
    all_text['Europa'] = 'Europe';
    all_text['Ganymede'] = 'Ganymède';
    all_text['Event'] = 'Phénomène';
    all_text['Shadow transit'] = "Passage d'ombre";
    all_text['Satellite transit'] = 'Passage';
    all_text['Egress'] = 'Fin';
    all_text['Ingress'] = 'Début';
    all_text['Disappearance'] = 'Disparition';
    all_text['Reappearance'] = 'Réapparition';
    
    // for the planisphere legend
    all_text['Radiant'] = 'Radiant';
    all_text['Messier'] = 'Messier';
    all_text['Comet'] = 'Comète';
        
    all_text['Chart'] = 'Carte'; //messier table    
  };
  populate_translations();
  var trans = function(english_key){
    var result = english_key;
    if (lang === 'fr'){
      //this data structure handles only one other lang; can be easily changed to handle N langs 
      if (all_text[english_key]){
        result = all_text[english_key];
      }
    }
    return result;
  };
  
  /* +1 for northern hemisphere, -1 for southern.*/
  var sign_latitude = function(){
    return input.latitude < 0 ? -1 : 1; //degs or rads, same condition
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
  var is_tweaker_param = function(item, tweakers){
    return array_includes(item, tweakers);
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
  link_back_to_main_form(['go', 'num_steps', 'date_time_unit']);
  

  var remove_decimal = function(text){
    return text.replace('.', '');
  };
  /* NOT WORKING. Example:   http://www.lightpollutionmap.info/#zoom=8&lat=5649485&lon=-8476439&layers=B0TFFFF */
  var light_pollution_link = function(){
    var lat = remove_decimal(input.latitude);
    var longit = remove_decimal(input.longitude);
    var result = 'http://www.lightpollutionmap.info/#zoom=8&lat=' + lat + '&lon=' + longit + '&layers=B0TFFFF';
    return result; 
  };  
  var google_maps_link = function(){
    return 'https://www.google.ca/maps/@' + input.latitude + ',' + input.longitude + ',11z';
  };
  var add_summary_blurb_at_top = function(){
    var summary = document.getElementById('summary');
    var when_text; 
    if (input.date_time.length > 0){
      when_text = input.time_scale + ' ' + input.date_time;
    }
    else {
      when_text = trans(when.toStringLT(WEEKDAYS));
      when_text = when_text.substring(0, when_text.length-7) + trans(when_text.substring(when_text.length-3));
    }
    var lmst = EPH.lmst(when, opts_rads.where);
    summary.innerHTML = NL + summary.innerHTML + NL + 
      '   <tr><td><a href="' +  google_maps_link() + '">' + UTIL.escapeHtml(input.location_name) + '</a>' +    
      '<td>' + when_text +  
      '<td>' + lmst.hour + ':' + lmst.min + ':' + lmst.sec.substring(0,2) +   
      '<td align="right">' +  UTIL.escapeHtml(input.latitude) + DEG +  
      '<td align="right">' +  UTIL.escapeHtml(input.longitude) + DEG + 
      '<td align="right">' +  new Date().getTimezoneOffset() + 'm'
    ;  
  };
  add_summary_blurb_at_top();

  /* Name of a compass point corresponding to degrees azimuth. */  
  var compass_point = function(degs){
    var width = 360/16;
    var idx = Math.trunc((degs + width*0.5) / width);
    return compass_points[idx];
  };

  /* Example: '23:02:15'. Return empty string if not found. */
  var phen_when = function(phen, type, name){
    var result = '';
    var when, last_idx, when_txt;
    var ph_idx = phen_idx(phen, type, name);
    if (ph_idx > -1){
      when = phen[type][ph_idx].when;
      //result = when.toStringLT().substring(14, 22);
      when_txt = when.toStringLT().substring(3); //chop off the LT
      if (type === 'horizons'){
        result = nearest_minute(when_txt).substring(10);
      }
      else {
        result = when.toStringLT().substring(14, 22);
      }
    }
    return result;
  };
  var chart_link = function(center_on, link_text){
    var result = 
      '<a title=' + trans('Chart') + ' href="../chart/graphic.sky?center_object=' + center_on + 
        '&location_name=' + encodeURIComponent(input.location_name) + 
        '&latitude=' + input.latitude + 
        '&longitude=' + input.longitude + 
        '&date_time=' + encodeURIComponent(input.date_time) + 
        '&time_scale=' + input.time_scale + 
      '">' + 
         link_text + 
       '</a>'  
    ;
    return result; 
  };
  var chart_link_for_date_time = function(center_on, link_text, date_time, time_scale){
    var result = 
      '<a title=' + trans('Chart') + ' href="../chart/graphic.sky?center_object=' + center_on + 
        '&location_name=' + encodeURIComponent(input.location_name) + 
        '&latitude=' + input.latitude + 
        '&longitude=' + input.longitude + 
        '&date_time=' + encodeURIComponent(date_time) + 
        '&time_scale=' + time_scale + 
      '">' + 
         link_text + 
       '</a>'  
    ;
    return result; 
  }
  var populate_sun_moon_planets_table = function(row){
    var thing = row.thing;
    var ephem = row.ephem;
    var phen = row.phen;
    var table = document.getElementById('sun_and_moon_and_planets_table');
    var elong = ephem.elong === undefined ? '' : EPH.round(ephem.elong,0)  + DEG; 
    var mag = ephem.mag === undefined ? ' ' : EPH.round_and_pad(ephem.mag, 1);
    var illum = ephem.illum === undefined ? ' ' : EPH.round(ephem.illum*100,0) + '%';
    var size_seconds = EPH.degs(ephem.size)*60*60;
    var size = size_seconds < 60 ?  EPH.round_and_pad(size_seconds,1) + '"' :  EPH.round_and_pad(size_seconds/60,1) + '\'';
    table.innerHTML = table.innerHTML + NL + 
      '   <tr>' + 
      '<td style="text-align:left; white-space: nowrap;">' +
        chart_link(thing.name, thing.symbol + ' ' + trans(thing.name)) +   
      '<td style="text-align:left;">' + ephem.zodiac.abbr + 
      '<td>' + elong + 
      '<td>' + mag + 
      '<td>' + size +  
      '<td>' + illum +  
      '<td style="text-align:right;">' + EPH.round(ephem.a, 0) + DEG +  
      '<td style="text-align:left;" title="'+ EPH.round(ephem.A,0) + DEG + '">' + trans(compass_point(ephem.A)) +     
      '<td style="text-align:right;">' + EPH.round_and_pad(ephem.α,2) + DEG +  
      '<td style="text-align:right;">' + EPH.round_and_pad(ephem.δ,2) + DEG +  
      '<td style="white-space: nowrap;">' + phen_when(phen, 'horizons', 'rise') +  
      '<td style="white-space: nowrap;">' + phen_when(phen, 'culminations', 'upper') + 
      '<td style="white-space: nowrap;">' + phen_when(phen, 'horizons', 'set')  
    ; 
  };
  var planet_rows = function(exclude_phenomena){
      var planets = EPH.as_array(EPH.planets);
      var rows = [];
      var phen, i, ephem, row;
      for (i=0; i < planets.length; ++i){
        ephem = EPH.position(planets[i].name, when, opts_rads);
        EPH.convert_all_angles_to_degs(ephem); 
        row = {thing: planets[i], ephem: ephem};
        if (!exclude_phenomena){
          if (planets[i].name === 'Sun'){
            phen = EPH.rise_culmination_set_daily(planets[i].name, when, opts_rads);
          }
          else {
            phen = EPH.rise_culmination_set_observation_window(planets[i].name, when, opts_rads);
          }
          row.phen = phen;
        }
        rows.push(row);
      }
      //rows.sort(by_time_of_set_transit_rise);
      rows.sort(by_elongation_from_the_sun);
      return rows;
  };
  var show_sun_moon_and_planets = function(){
    var rows, i;
    var table = document.getElementById('sun_and_moon_and_planets_table');
    if (input.exclude_sun_moon_planets === '1'){
      dont_show(table);
    }
    else {
      rows = planet_rows(false);
      for (i=0; i < rows.length; ++i){
        populate_sun_moon_planets_table(rows[i]);
      }
    }
    return rows;
  };
  var planets_etc = show_sun_moon_and_planets();

  var standard_canvas_appearance = function(ctx, big_font){
    if (big_font){
      GRAPH.myBigFont(ctx);
    }
    ctx.fillStyle = 'rgb(255,255,255)'; 
    ctx.strokeStyle = 'rgb(255,255,255)';
    ctx.textAlign='center';
    ctx.textBaseline = 'middle';
  };


    
  var print_zodiac_stars = function(names, sun, canvas, ctx){
    var i, j, star, ephem, pos, sign;
    EPH.convert_all_angles_to_rads(sun);
    for (i=0; i < names.length; i++){
      for (j=0; j < EPH.stars.length; j++){
        if (EPH.stars[j][0] === names[i]) {
          star = {α:EPH.stars[j][1], δ:EPH.stars[j][2]};
          EPH.convert_αδ_to_λβ(star, when); //adds beta, rads
          star.elong = EPH.elongation_between(sun, star, when); //rads, unsigned
          sign = EPH.delta_longitude_between(star, sun) < 0 ? -1 : 1;
          star.elong = sign * star.elong; // rads, signed
          ephem = {β:EPH.degs(star.β), elong: EPH.degs(star.elong)}; // degs, signed
          pos = ecliptic_rendering_position(ephem, canvas.height, canvas.width);
          GRAPH.tickMarkVertical(ctx,pos.x,pos.y,2);
          GRAPH.tickMarkHorizontal(ctx,pos.x,pos.y,2.5);
        }  
      }
    }
    EPH.convert_all_angles_to_degs(sun); //back to the units used by the caller
  };
  var ecliptic_rendering_position = function(ephem /*.elong .β in degs*/, h, w){
    var x, y, LAT_MAX = 10, sign = sign_latitude(); 
    y = (h/2) - sign * (h/2) * ephem.β/LAT_MAX; //northern hemis: y increases towards negative beta
    var dx = Math.abs((w/2)*(ephem.elong/180));
    if (sign === +1){
      x = ephem.elong < 0 ? dx : w - dx;
    }
    else {
      x = ephem.elong < 0 ? w - dx : dx;
    }
    return {x:x, y:y};
  };
  var print_zodiac = function(idx, x, canvas, ctx){
    var name = EPH.zodiac[idx].abbr;
    var tweak = 10;
    if (canvas.width < 650) {
      name = name.substring(0,1);
      ctx.font = "6px Arial";
      tweak = 5;
    }
    GRAPH.text(ctx, name, x, canvas.height - tweak);
  };
  var show_zodiac_sign_name_edge_case = function(a, b /*the start and end of the sign*/, canvas, ctx, end){
    var midpoint, FRAC = 24;
    //in some cases, the same text can appear at both ends of the chart
    if (a > (canvas.width/FRAC)){
      midpoint = (a+0)/2; //left side
      print_zodiac(end, midpoint, canvas, ctx);
    }
    if ((canvas.width - b) > (canvas.width/FRAC)){
      midpoint = (b+canvas.width)/2; //right side
      print_zodiac(end, midpoint, canvas, ctx);
    }
  };
  var show_zodiac_sign_name_on_ecliptic = function(zodiac_x, canvas){
    var i, a, b, midpoint, end, sign_lat = sign_latitude();
    var ctx = canvas.getContext('2d');
    ctx.textAlign = "center";
    for(i = 0; i < zodiac_x.length; ++i){
      end = (i+1)%12; //goes back to Psc for the last iteration
      a = zodiac_x[i];
      b = zodiac_x[end]; 
      if (Math.abs(a - b) > (canvas.width/2)){
        //widely separated; it's at both ends; edge cases - literally!
        if (sign_lat === +1){
          show_zodiac_sign_name_edge_case(a, b, canvas, ctx, end);
        }
        else {
          show_zodiac_sign_name_edge_case(b, a, canvas, ctx, end); //southern hemis: switch a and b
        }
      }
      else {
        midpoint = (a+b)/2;
        print_zodiac(end, midpoint, canvas, ctx);
      }
    }
  };
  var show_zodiac_hash_marks_on_ecliptic = function(sun_λ, canvas){
    var i, Δλ, x, w, Δx, ctx;
    var result = [];
    ctx = canvas.getContext('2d');
    w = canvas.width;
    var sign = sign_latitude();
    for(i = 0; i < EPH.zodiac.length - 1; ++i){
      Δλ = EPH.zodiac[i].λ_end - sun_λ;
      Δx = Math.abs((Δλ*w)/(2*Math.PI)); //always positive
      if (sign === +1){
        x = (Δλ >= 0) ? w - Δx : Δx; //north 
      }
      else {
        x = (Δλ >= 0) ? Δx : w - Δx; //south 
      }
      result.push(x); //index is parallel to that of zodiac
      GRAPH.tickMarkVertical(ctx,x,canvas.height,8); // ticks at bottom
    }
    return result;
  };
  /* 
   Shows ecliptic latitude and elongation from the Sun.
   Should this show the difference in ecliptic longitude instead?
  */  
  var show_ecliptic = function(planets_etc){
    var planets, i, x, y, point, h, w, tweak, text_tweak, zodiac_x, pos;
    var LAT_MAX = 7; 
    var TICK_SIZE = 10;
    var canvas = document.getElementById('ecliptic');
    var ctx = canvas.getContext('2d');
    canvas.width = window.innerWidth * 0.95;
    canvas.height = canvas.width * 0.12;
    standard_canvas_appearance(ctx, true);
    h = canvas.height;
    w = canvas.width;
    if (input.exclude_ecliptic === '1'){
      dont_show(canvas);
    }
    else {
      GRAPH.line(ctx, 0, h/2, w, h/2); //horizontal ecliptic in the middle
      for(i=1; i <= 7; ++i){ // ticks at top, intervals of 45 degrees
        tweak = 1;
        if (i == 4) tweak = 2;
        if (i == 2 || i == 6) tweak = 1.5;
        GRAPH.tickMarkVertical(ctx,i*(w/8),0,TICK_SIZE*tweak); 
      }
      zodiac_x = show_zodiac_hash_marks_on_ecliptic(EPH.rads(planets_etc[0].ephem.λ), canvas);
      show_zodiac_sign_name_on_ecliptic(zodiac_x, canvas);
      print_zodiac_stars(['α Leo', 'α Vir', 'α Tau', 'α Sco'], planets_etc[0].ephem, canvas, ctx);      
      //title goes here
      planets = planets_etc.length > 0 ? planets_etc : planet_rows(true);
      for(i=0; i < planets.length; ++i){
        pos = ecliptic_rendering_position(planets[i].ephem, h, w);
        GRAPH.square(ctx,pos.x,pos.y,3); //spot for each planet
        text_tweak = planets[i].ephem.β > 0 ? 15 : -15; // text tends towards the center
        GRAPH.text(ctx,planets[i].thing.symbol,pos.x,pos.y+text_tweak);
      }
    }
  };
  show_ecliptic(planets_etc);

  
  var num_hours_after = function(when, start_hour){
    var time = when.toStringLT().substring(14); // LT 2016-05-13 01:01:01.132
    var parts = time.split(":");
    var hour = parseInt(parts[0]);
    if (hour < start_hour) {
      hour = hour + 24; //after midnight, next day
    }
    var result = hour+parseInt(parts[1])/60+parseFloat(parts[2])/3600 - start_hour;
    return result;
  };

  /* 
   Decide which times to use (if any) to indicate when the moon makes the sky brighter. 
   Returns null if no moon_box is to be shown.
   Returns object with two whens: .start, and .end.
  */
  var moon_box = function(obs_window, moon /*rise-set phenomena, not ephem*/) {
    var result, start, end;
    var moon_at_sunset = EPH.position('moon', obs_window.sunset, opts);
    var moon_at_sunrise = EPH.position('moon', obs_window.sunrise, opts);
    var moon_rise_idx = phen_idx(moon, 'horizons', 'rise');
    var moon_set_idx = phen_idx(moon, 'horizons', 'set');
    if (moon_at_sunset.a > 0){
      start = obs_window.sunset;
    }
    else if (moon_rise_idx > -1) { 
      if (obs_window.sunset.jd <= moon.horizons[moon_rise_idx].when.jd &&  moon.horizons[moon_rise_idx].when.jd <= obs_window.sunrise.jd){
        start = moon.horizons[moon_rise_idx].when;
      }
    }
    if (start){
      if (moon_at_sunrise.a > 0 || moon_set_idx === -1){
        end = obs_window.sunrise;
      }
      else {
        end = moon.horizons[moon_set_idx].when;
      }
      result = {
        start: start,
        end: end
      };
    }
    return result;
  };
  var draw_lunar_terminator = function (ctx, x1, y1, x2, y2, start_angle, end_angle) {
      //See: http://stackoverflow.com/questions/21594756/drawing-circle-ellipse-on-html5-canvas-using-mouse-events
      var radiusX = (x2 - x1) * 0.5,   /// radius for x based on input
          radiusY = (y2 - y1) * 0.5,   /// radius for y based on input
          centerX = x1 + radiusX,      /// calc center
          centerY = y1 + radiusY,
          step = 0.10,                 /// resolution of ellipse
          a = start_angle;            /// counter
      /// start a new path
      ctx.beginPath();
      /// begin with start angle
      ctx.moveTo(centerX + radiusX * Math.cos(a),
                 centerY + radiusY * Math.sin(a));
      /// rotate thru to the end angle
      for(; a <= (end_angle+step); a += step) {
          ctx.lineTo(centerX + radiusX * Math.cos(a),
                     centerY + radiusY * Math.sin(a));
      }
      /// close it and stroke it for demo
      //ctx.closePath();
      //ctx.strokeStyle = '#000';
      ctx.stroke();
  };
  /* Draw a simple moon phase chart. */
  var show_moon_phase = function(ctx, ctr, r /* characteristic radius */){
    ctx.save();
    ctx.translate(ctr.x, ctr.y); // draw using the ctr as the natural origin of coord
    if (sign_latitude() === -1){
      ctx.scale(-1, 1); //flip left-right, for southern hemisphere    
    }
    var moon = EPH.position('Moon', when);
    //remember that elongation is a bit strange: the sign indicates evening or morning elongation, so there's 
    //a discontinuity at pi
    var elong = moon.elong > 0 ? moon.elong : moon.elong + 2*Math.PI; //0..2pi
    var north = {x: 0, y: -r};
    var south = {x: 0, y: +r};
    var b = Math.abs(Math.cos(elong) * r); //minor axis of an ellipse
    var waxing = elong < Math.PI; 
    var waning = ! waxing;
    var crescent = elong < 0.5*Math.PI || elong > 1.5*Math.PI;
    var gibbous = ! crescent;
    var x1 = - b;
    var x2 = + b;
    var y1 = - r;
    var y2 = + r;
    ctx.beginPath();
    if ((waxing && crescent) || (waning && gibbous)){
      //ctx.ellipse(0, 0, b, r, 0, -0.5*Math.PI, 0.5*Math.PI); //ellipse on right; not widely supported
      draw_lunar_terminator(ctx, x1, y1, x2, y2, -0.5*Math.PI, 0.5*Math.PI);      
    }
    else if ((waxing && gibbous) || (waning && crescent)){
      //ctx.ellipse(0, 0, b, r, 0, 0.5*Math.PI, 1.5*Math.PI); //ellipse on left; not widely supported
      draw_lunar_terminator(ctx, x1, y1, x2, y2, 0.5*Math.PI, 1.5*Math.PI);      
    }
    if (waxing){
      ctx.moveTo(north.x, north.y);
      ctx.arc(0, 0, r, -0.5*Math.PI, 0.5*Math.PI, false); //half circle on right, starts at north
    }
    else {
      ctx.moveTo(south.x, south.y);
      ctx.arc(0, 0, r, 0.5*Math.PI, -0.5*Math.PI, false); //half circle on left, starts at south
    }
    ctx.stroke();
    ctx.closePath();
    ctx.restore();
  };
  /* Far north has extreme values. In winter the window is often very wide. Return .num_hours */
  var observation_window_config = function(obs_window){
    //defaults for southernish latitudes:
    var start = 15, end = 9;
    //override for the north, increasing the number of charted hours
    if (obs_window.sunset.date.getHours() < 16){
      start = obs_window.sunset.date.getHours() - 1;
    }
    if (obs_window.sunrise.date.getHours() > 7){
      end = obs_window.sunrise.date.getHours() + 2;
    }
    return {
      start_hour: start,
      end_hour: end,
      num_hours: 24+end-start
    };
  };
  var sunset_gradient = function(x, ctx, blue, h, is_set){
    var orange_width = 8;
    var gradient; 
    if (is_set){
      gradient = ctx.createLinearGradient(x-orange_width, 0, x, 0);
    }
    else {
      gradient = ctx.createLinearGradient(x, 0, x+orange_width, 0);
    }
    var sky_sunset = "#FF9900";
    var start_color = is_set ? blue : sky_sunset;
    var end_color = is_set ? sky_sunset : blue;
    gradient.addColorStop(0, start_color);
    gradient.addColorStop(1, end_color);
    if (is_set){
      GRAPH.rect(ctx, x-orange_width, 0, orange_width, h, gradient);
    }
    else {
      GRAPH.rect(ctx, x, 0, orange_width, h, gradient);
    }
  };
  var show_observation_window = function(when, where, planets_etc){
    var obs_window, obs_window_config, planets, i, x, y, point, h, w, tweak, moon_phen, moon_illum, hour, ctr;
    var sunrise_x, sunset_x, twilight_rise_x, twilight_set_x, w_hour;
    var moon_events, moon_x_start, moon_x_end, moon_color;
    //var LAT_MAX = 7; 
    var TICK_SIZE = 10;
    var canvas = document.getElementById('observation_window');
    var ctx = canvas.getContext('2d');
    canvas.width = window.innerWidth * 0.95;
    canvas.height = canvas.width * 0.10;
    standard_canvas_appearance(ctx, true);
    h = canvas.height;
    w = canvas.width;
    if (input.exclude_observation_window === '1'){
      dont_show(canvas);
    }
    else {
      obs_window = EPH.observation_window(when, opts_rads, parseInt(input.twilight)); 
      planets = planets_etc.length > 0 ? planets_etc : planet_rows(false);
      for(i=0; i < planets.length; ++i){
        if (planets[i].thing.name === 'Moon'){ 
          moon_illum = planets[i].ephem.illum;
          moon_phen = planets[i].phen;
          break;
        }
      }
      //indicate daylight hours 
      //calc the x corresponding to each time-diff from a star_hour
      obs_window_config = observation_window_config(obs_window);
      w_hour = w/obs_window_config.num_hours;
      sunset_x = w_hour*num_hours_after(obs_window.sunset, obs_window_config.start_hour);
      sunrise_x = w_hour*num_hours_after(obs_window.sunrise, obs_window_config.start_hour);
      var blue = 'rgb(36,55,114)';
      GRAPH.rect(ctx,0,0,sunset_x,h, blue);
      GRAPH.rect(ctx,sunrise_x,0,w,h,blue);
      //indicate the time near sunset/rise with an orange-y gradient
      sunset_gradient(sunset_x, ctx, blue, h, true);
      sunset_gradient(sunrise_x, ctx, blue, h, false);
      //indicate twilight hours 
      var dusk = 'rgb(125,125,125)';
      if (obs_window.twilight_start && obs_window.twilight_end){
        twilight_set_x = w_hour*num_hours_after(obs_window.twilight_end, obs_window_config.start_hour);
        twilight_rise_x = w_hour*num_hours_after(obs_window.twilight_start, obs_window_config.start_hour);
        GRAPH.rect(ctx, sunset_x, 0, twilight_set_x - sunset_x, h, dusk);
        GRAPH.rect(ctx,twilight_rise_x, 0, sunrise_x - twilight_rise_x, h, dusk);
      }
      else {
        //the entire observation window is in twilight
        GRAPH.rect(ctx, sunset_x, 0, (sunrise_x - sunset_x), h, dusk);
      }
      //box-like thing for moonrise, moonset; height of the box as % of canvas is its percent illumination
      moon_events = moon_box(obs_window, moon_phen);
      if (moon_events !== undefined){
        moon_x_start = w_hour*num_hours_after(moon_events.start, obs_window_config.start_hour);
        moon_x_end = w_hour*num_hours_after(moon_events.end, obs_window_config.start_hour);
        var moon_color = 'rgb(170,170,170)';
        GRAPH.rect(ctx, moon_x_start, 0, moon_x_end - moon_x_start, h*moon_illum, moon_color);
        if (moon_illum > 0.15){
          ctr = {x: moon_x_start + 0.5*(moon_x_end-moon_x_start), y: 0.5*h*moon_illum};
          show_moon_phase(ctx, ctr, 0.20*h*moon_illum);
        }
      }
      //indicate hours as tick marks at the top      
      for(i=1; i <= obs_window_config.num_hours - 1; ++i){
        x = i*w_hour;
        GRAPH.tickMarkVertical(ctx,x,0,TICK_SIZE); // ticks at top, to mark hours (don't do the ends of the interval)
        hour = obs_window_config.start_hour + i;
        if (hour > 23){
          hour = hour - 24;
        }
        if (hour < 10) {
          hour = '0' + hour;
        }
        GRAPH.text(ctx,hour,x,20);
      }
    }
  };
  show_observation_window(when, opts.where, planets_etc);

  var show_galilean_satellites = function(){
    var jupiter, satellites, x, y, h, w, dist, dx, dy, sign;
    var canvas = document.getElementById('galilean_satellites');
    var ctx = canvas.getContext('2d');
    canvas.width = window.innerWidth * 0.45;
    canvas.height = window.innerWidth * 0.114;
    standard_canvas_appearance(ctx, false);
    h = canvas.height;
    w = canvas.width;
    if (input.exclude_galilean_satellites === '1'){
      dont_show(canvas);
    }
    else {
      var JUPITERS_RADIUS = w * 0.015; //the coords of the satellites take this as unit distance 
      jupiter =  {x: w/2, y: h/2, radius: JUPITERS_RADIUS};
      GRAPH.spot(ctx, jupiter.x, jupiter.y, jupiter.radius); //Jupiter's disk in the middle
      satellites = EPH.physical_jupiter(when).satellites;
      sign = sign_latitude();
      for(i=0; i < satellites.length; ++i){
        x = jupiter.x + satellites[i].X*jupiter.radius;
        y = jupiter.y - sign * satellites[i].Y*jupiter.radius; 
        dx = x - jupiter.x;
        dy = y - jupiter.y;
        dist = Math.sqrt(dx*dx + dy*dy);
        if (dist > jupiter.radius){
          GRAPH.square(ctx,x,y,3); 
          GRAPH.text(ctx,satellites[i].symbol,x,y+15);
        }
      }
    }
  };
  show_galilean_satellites();
  
  var add_to_diary_table = function(events){
    var table = document.getElementById('diary_table');
    var when_string;
    if (input.exclude_sky_diary === '1' || events.length === 0){
      dont_show(table);
    }
    else {
      for(var i=0; i < events.length; i++){
        when_string = events[i].when.toStringLT(WEEKDAYS);
        when_string = when_string.substring(3, when_string.length-9) + when_string.substring(when_string.length-4);
        when_string = when_string.substring(0, when_string.length-2) + trans(when_string.substring(when_string.length-2)); 
        table.innerHTML = table.innerHTML + NL + 
          '   <tr><td>' + when_string +    
          '<td style="text-align:left;">' + events[i].text 
        ; 
      }
    }
  };
  add_to_diary_table(EPH.current_events(when, 7));
  
  var add_to_deep_sky_object_table = function(rows, node_name, link_fn /*optional*/){
    var thing, ephem;
    var table = document.getElementById(node_name);
    if (input['exclude_' + node_name] === '1' || rows.length === 0){
      dont_show(table);
    }
    else {
      for(var i = 0; i < rows.length; i++){
        thing = rows[i].thing;
        ephem = rows[i].ephem;
        var link = link_fn ? link_fn(thing) : thing[0];
        var comment = thing[6] ? thing[6] : thing[7];
        var mag = thing[3] ? EPH.round_and_pad(thing[3], 1) : '';
        table.innerHTML = table.innerHTML + NL + 
          '   <tr><td>' + link +
          '<td>' + thing[4] +
          '<td>' + chart_link(thing[0], trans('Chart')) +
          '<td style="text-align:right;">' + mag +  
          '<td style="text-align:right;">' + EPH.round(ephem.a,0) + DEG +
          '<td title="'+ EPH.round(ephem.A,0) + DEG + '">' + trans(compass_point(ephem.A)) +    
          '<td style="text-align:right;">' + EPH.round_and_pad(ephem.α,2) + DEG +  
          '<td style="text-align:right;">' + EPH.round_and_pad(ephem.δ,2) + DEG +   
          '<td>' + thing[5] +
          '<td>' + comment
        ;
      }
    }
  };
  var messier_link = function(thing){
    var messier_num = thing[0].substring(1, thing[0].length); //chop off the 'M'
    return "<a href='https://en.wikipedia.org/wiki/Messier_" + messier_num + "'>" + thing[0] + "</a>";
  };
  var caldwell_link = function(thing){
    return "<a href='https://en.wikipedia.org/wiki/" + thing[8] + "'>" + thing[0] + "</a>";
  };
  var min_alt = 25;
  var visible_messiers = EPH.find_visible_messiers(when, opts.where, min_alt, parseFloat(input.limiting_mag_messiers), opts);
  add_to_deep_sky_object_table(visible_messiers, 'messiers', messier_link);
  var visible_caldwells = EPH.find_visible_caldwells(when, opts.where, min_alt, parseFloat(input.limiting_mag_messiers), opts);
  add_to_deep_sky_object_table(visible_caldwells, 'caldwells', caldwell_link);
  
  var populate_table = function(things, node_name, exclude_flag){
    var table = document.getElementById(node_name);
    if (things.length === 0 || input[exclude_flag] === '1'){
      dont_show(table);
    }
    else {
      var rows = [];
      var i, ephem, thing;
      for (i=0; i < things.length; ++i){
        ephem = EPH.position(things[i].name, when, opts);
        rows.push({thing: things[i], ephem: ephem});
      }
      rows.sort(highest_first);
      for(i=0; i < rows.length; ++i){
        thing = rows[i].thing;
        ephem = rows[i].ephem;
        table.innerHTML = table.innerHTML + NL + 
          '   <tr><td>'  
          + chart_link(thing.name, thing.alt_name) +   
          '<td  style="text-align:right;">' + EPH.round(ephem.elong, 0) + DEG +
          '<td style="text-align:right;">' + EPH.round_and_pad(ephem.mag,1) + 
          '<td style="text-align:right;">' + EPH.round(ephem.a, 0) + DEG +  
          '<td title="'+ EPH.round(ephem.A,0) + DEG + '">' + trans(compass_point(ephem.A)) +     
          '<td style="text-align:right;">' + EPH.round_and_pad(ephem.α,2) + DEG +  
          '<td style="text-align:right;">' + EPH.round_and_pad(ephem.δ,2) + DEG 
        ;
      }
    }
  };
  populate_table(EPH.as_array(EPH.minor_planets), 'minor_planet_table', 'exclude_minor_planets');

  var populate_comets_table = function(things, node_name, exclude_flag){
    var table = document.getElementById(node_name);
    if (things.length === 0 || input[exclude_flag] === '1'){
      dont_show(table);
    }
    else {
      var rows = [];
      var i, ephem, thing;
      for (i=0; i < things.length; ++i){
        ephem = EPH.position(things[i].name, when, opts);
        rows.push({thing: things[i], ephem: ephem});
      }
      rows.sort(highest_first);
      for(i=0; i < rows.length; ++i){
        thing = rows[i].thing;
        ephem = rows[i].ephem;
        table.innerHTML = table.innerHTML + NL + 
          '   <tr><td>' 
          + chart_link(thing.name, thing.alt_name) +   
          '<td  style="text-align:right;">' + EPH.round(ephem.elong, 0) + DEG +
          '<td style="text-align:right;">' + EPH.round_and_pad(thing.mag,1) + 
          '<td style="text-align:right;">' + EPH.round(ephem.a, 0) + DEG +  
          '<td title="'+ EPH.round(ephem.A,0) + DEG + '">' + trans(compass_point(ephem.A)) +     
          '<td style="text-align:right;">' + EPH.round_and_pad(ephem.α,2) + DEG +  
          '<td style="text-align:right;">' + EPH.round_and_pad(ephem.δ,2) + DEG +  
          '<td>' + trans(thing.trend) +  
          '<td>' + trans(thing.when_vis)  
        ;
      }
    }
  };
  populate_comets_table(EPH.as_array(EPH.comets), 'comet_table', 'exclude_comets');
  

  var add_to_meteor_shower_table = function(showers){
    var table = document.getElementById('meteor_shower_table');
    var zhr_factor;
    if (input.exclude_meteor_showers === '1' || showers.length === 0){
      dont_show(table);
    }
    else {
      for(var i=0; i < showers.length; i++){
        zhr_factor = showers[i].zhr_factor ? EPH.round_and_pad(showers[i].zhr_factor, 2) : '';
        var title1 = trans('Population index; smaller means brighter');
        var title2 = trans('Factor by which ZHR is multiplied, for current alt and limiting mag');
        table.innerHTML = table.innerHTML + NL + 
          '   <tr><td>' + showers[i].name +
          '<td>' + showers[i].symbol + 
          '<td>' + showers[i].peak_time.toStringLT().substring(2, 19) + 
          '<td>' + showers[i].v +  
          '<td>' + EPH.round_and_pad(showers[i].radiant.a,1) + DEG +
          '<td>' + EPH.round_and_pad(showers[i].radiant.A,1) + DEG + 
          '<td title="' + title1 + '">' + showers[i].r + 
          '<td>' + showers[i].zhr + 
          '<td title="' + title2 + '">' + zhr_factor
        ;
      }
    }
  };
  add_to_meteor_shower_table(EPH.current_meteor_showers(when, opts));

  /* WMTS tile server; input long-lat. */
  var show_radar_uk_images = function(){
      var parser = new ol.format.WMTSCapabilities();
      //http://datapoint.metoffice.gov.uk/public/data/inspire/view/wmts?REQUEST=GetCapabilities&key=c8258e8b-bf20-45ed-8e95-ca0d76682419
      fetch('http://datapoint.metoffice.gov.uk/public/data/inspire/view/wmts?REQUEST=GetCapabilities&key=' + MET_OFFICE_API_KEY).then(function(response) {
        return response.text();
      }).then(function(text) {
        var result = parser.read(text); //an openlayers object
        
        //Ref: http://stackoverflow.com/questions/41526983/openlayers3-wmts-out-of-memory-error/41564330?noredirect=1#comment70434808_41564330
        
        //ERROR HANDLING - AT LEAST TO THE CONSOLE
        
        //height of the div is too short
        
        //add the eumetsat satellite data
        
        // add the url of the get-capabilities request, to document here what the tiles are etc
        
        //I need to add their 'blurb' on usage somewhere
        
        //the appearance is dark grey - too dark
        
        //1. needed since the GetCapabilities links doesn't include my key:
        result.OperationsMetadata.GetTile.DCP.HTTP.Get[0].href += 'key=' + MET_OFFICE_API_KEY + '&'
        
        //2.fix: correct the bad data coming back from the server:
        result.Contents.TileMatrixSet[0].TileMatrix.forEach(function(m) {
          m.ScaleDenominator *= 111319.49079327358;
          m.TopLeftCorner = m.TopLeftCorner.reverse();
        });
        
        var options = ol.source.WMTS.optionsFromCapabilities(result, {
            layer: 'RADAR_UK_Composite_Highres', 
            matrixSet: 'EPSG:4326'
        });
        var street_map_tile = new ol.layer.Tile({source: new ol.source.OSM(), opacity: 0.7});
        var radar_tile = new ol.layer.Tile({source: new ol.source.WMTS(options), opacity: 0.7});
        var view = new ol.View({
          center: ol.proj.fromLonLat([parseFloat(input.longitude), parseFloat(input.latitude)]),  
          zoom: 7
        }); 
        var map = new ol.Map({
          layers: [street_map_tile, radar_tile],
          target: 'radar_uk',
          view: view,
          controls: []
        });
        var radar_div = document.getElementById('radar_uk');
        var canv = radar_div.children[0].children[0];
        console.log('radar height of canvas: ' + canv.height); //150 in console, 211 in the debugger
        //canv.height = 300; 
      });
  };
  //Example: http://radar.weather.gov/radar.php?rid=lwx&product=N0R&overlay=11101111&loop=no
  var show_radar_us_images = function(radar_station_us){
    var canvas = document.getElementById('radar_us');
    var ctx = canvas.getContext('2d');
    var topo_url, radar_url;
    var topo_img = new Image();
    var radar_img;
    topo_img.onload = function(){
      console.log('Drawing US radar image (topography).');
      canvas.height = topo_img.height;
      canvas.width = topo_img.width;
      ctx.drawImage(topo_img, 0, 0);
      radar_img = new Image();
      radar_img.onload = function(){
        console.log('Drawing US radar image (radar).');
        ctx.drawImage(radar_img, 0, 0);
      };
      radar_img.onabort = function(){
        console.log('Abort load of US radar image (radar).');
      };
      radar_img.onerror = function(){
        console.log('Error load of US radar image (radar).');
      };
      radar_img.style.display = 'none';
      radar_url = 'https://radar.weather.gov/ridge/RadarImg/N0R/' + radar_station_us +'_N0R_0.gif';
      console.log('US radar : ' + radar_url);
      radar_img.src = UTIL.crossDomainUrl(radar_url); 
    };
    topo_img.onabort = function(){
      console.log('Abort load of US radar image (topography).');
    }
    topo_img.onerror = function(){
      console.log('Error load of US radar image (topography).');
    }
    topo_img.style.display = 'none';
    //currently, the NOAA server doesn't support CORS; hence need a workaround
    topo_url = 'https://radar.weather.gov/ridge/Overlays/Topo/Short/' + radar_station_us + '_Topo_Short.jpg';
    console.log('US radar (topography) : ' + topo_url);
    topo_img.src = UTIL.crossDomainUrl(topo_url); 
  };
  var is_covered_by_meteosat = function(){
    return (-70 <= opts.where.φ) && (opts.where.φ <= 70) && (-70 <= opts.where.λ) && (opts.where.λ <= 70);
    //return false; 
  };
  var show_clouds_and_radar = function(){
    //clouds, from a satellite image 
    if (input.exclude_clouds === '1'){
      dont_show(document.getElementById('clouds'));
    }
    else {
      var satellite;
      if (input.country === 'cda' || input.country === 'us'){
         satellite = 'goes';
      }
      else if (is_covered_by_meteosat()){
         satellite = 'meteosat';
      }
      if (satellite){
        var sun = EPH.position('sun', EPH.when_now(), opts);
        showClouds(input, sun.a, 'clouds', satellite);
      }
      else {
        dont_show(document.getElementById('clouds'));
      }
    }
    //radar
    if (input.exclude_radar === '1' || input.country === 'other'){
      //suppress radar
      dont_show(document.getElementById('radar'));
      dont_show(document.getElementById('radar_us'));
    }
    else if ('cda' === input.country){
      //Canadian radar only
      dont_show(document.getElementById('radar_us'));
      dont_show(document.getElementById('radar_uk'));
      showPrecipitation(input.radar_station, is_dev);
    }
    else if ('us' === input.country){
      //American radar only
      dont_show(document.getElementById('radar'));
      dont_show(document.getElementById('radar_uk'));
      show_radar_us_images(input.radar_station);
    }
    else if ('uk' === input.country){
      //UK radar only
      dont_show(document.getElementById('radar'));
      dont_show(document.getElementById('radar_us'));
      //show_radar_uk_images();
    }
  };
  show_clouds_and_radar();

  /* why is a refresh usually needed? what's going on there? */
  var show_clear_sky_clock = function(){
    var csc = document.getElementById('clear_sky_clock');
    var is_north_america = (input.country === 'cda') || (input.country ===  'us');
    if (input.exclude_clear_sky_clock === '1' || ! is_north_america ){
      dont_show(csc);
    }
    else {
      csc.title = trans('Clear Sky Clock, predicted sky conditions');
      var station_id = input.clear_sky_clock_station_id;
      if (station_id){
        csc.src = 'http://www.cleardarksky.com/c/' + station_id + 'csk.gif';
      }
    }
  };
  show_clear_sky_clock();

  var fetch_from_web = function(url, callback, response_type){
    var xhr = new XMLHttpRequest();
    xhr.open('get', url, true);
    if (response_type) {
      xhr.responseType = response_type;
    }
    xhr.onload = function() {
      console.log('Success: fetched from ' + url);
      var status = xhr.status;
      if (status == 200) {
        callback(null, xhr.response);
      } else {
        callback(status);
      }
    };
    xhr.ontimeout = function (e) {
      console.log('Timeout: fetching from ' + url);
    };
    xhr.timeout = 20 * 1000;  
    xhr.send();
  };
  var fetch_json = function(url, callback) {
    return fetch_from_web(url, callback, 'json');
  };
  var fetch_xml = function(url, callback) {
    return fetch_from_web(url, callback);
  };
  var add_to_weather_table = function(name, value, table){
    table.innerHTML = table.innerHTML + NL + '   <tr><td><b><span>' + name + '</span></b><td>' + value; 
  };
  var us_weather_url = function(output_type){
    var result = 'http://forecast.weather.gov/MapClick.php?lat=' + EPH.degs(opts.where.φ) + '&lon=' + EPH.degs(opts.where.λ) + '&FcstType=' + output_type;
    //console.log(result);
    return result;
  };
  // Example URL:  http://forecast.weather.gov/MapClick.php?lat=38.4247341&lon=-86.9624086&FcstType=json
  var show_current_weather_us = function() {
    var table = document.getElementById('weather');
    fetch_json(us_weather_url('json'), function(err, weather){
      if (err != null){
        console.log('Unable to fetch JSON data for US weather: ' + url);
      }
      else {
        table.innerHTML = table.innerHTML + NL + '<caption><a href="' + us_weather_url('html') + '">' + trans('Weather') + '</a>';
        add_to_weather_table('Temperature', weather.currentobservation.Temp + '&deg;F', table);
        add_to_weather_table('Dew Point', weather.currentobservation.Dewp + '&deg;F', table);
        add_to_weather_table('Wind Speed', weather.currentobservation.Winds + 'mph', table);
        add_to_weather_table('Wind Direction', weather.currentobservation.Windd + '&deg;', table);
        add_to_weather_table('Conditions', weather.currentobservation.Weather, table);
        add_to_weather_table('Relative Humidity', weather.currentobservation.Relh + '%', table);
      }
    });
  };
  /* Example: http://datapoint.metoffice.gov.uk/public/data/val/wxobs/all/json/3772?res=hourly&key=c8258e8b-bf20-45ed-8e95-ca0d76682419 */
  var weather_url_uk = function(output_type){
    var result = 'http://datapoint.metoffice.gov.uk/public/data/val/wxobs/all/' + output_type + '/' + input.weather_station + '?res=hourly&key=' + MET_OFFICE_API_KEY;
    return result;
  };
  /* Source: http://www.metoffice.gov.uk/datapoint/support/documentation/code-definitions */
  var uk_weather_codes = function(){
    var result = {};
    result.NA='';
    result[0]='Clear night';
    result[1]='Sunny day';
    result[2]='Partly cloudy';
    result[3]='Partly cloudy';
    result[4]='';
    result[5]='Mist';
    result[6]='Fog';
    result[7]='Cloudy';
    result[8]='Overcast';
    result[9]='Light rain shower';
    result[10]='Light rain shower';
    result[11]='Drizzle';
    result[12]='Light rain';
    result[13]='Heavy rain shower';
    result[14]='Heavy rain shower';
    result[15]='Heavy rain';
    result[16]='Sleet shower';
    result[17]='Sleet shower';
    result[18]='Sleet';
    result[19]='Hail shower';
    result[20]='Hail shower';
    result[21]='Hail';
    result[22]='Light snow shower';
    result[23]='Light snow shower';
    result[24]='Light snow';
    result[25]='Heavy snow shower';
    result[26]='Heavy snow shower';
    result[27]='Heavy snow';
    result[28]='Thunder shower';
    result[29]='Thunder shower';
    result[30]='Thunder';  
    return result;
  };
  var show_current_weather_uk = function(){
    var table = document.getElementById('weather');
    fetch_json(weather_url_uk('json'), function(err, weather){
      if (err != null){
        console.log('Unable to fetch JSON data for UK weather: ' + url);
      }
      else {
        //var details = weather.SiteRep.DV.Location.Period[0].Rep[0];
        var per = weather.SiteRep.DV.Location.Period[0];
        var details = per.Rep[0] ? per.Rep[0] : per.Rep; //inconsistent? changes from array to object, when there's only 1 item
        table.innerHTML = table.innerHTML + NL + '<caption><a href="http://www.metoffice.gov.uk/public/weather/forecast">' + trans('Weather') + '</a>'; //can't map the station id, to link to a corresponding HTML page??
        if (details.T) add_to_weather_table('Temperature', details.T + '&deg;C', table);
        if (details.Dp) add_to_weather_table('Dew Point', details.Dp + '&deg;C', table);
        if (details.S) add_to_weather_table('Wind', details.S + ' m/h ' + details.D, table);
        if (details.W) add_to_weather_table('Conditions', uk_weather_codes()[details.W], table);
        if (details.H) add_to_weather_table('Relative Humidity', details.H + '%', table);
        var tendency = details.Pt === 'F' ? 'falling' : 'rising'; 
        if (details.P) add_to_weather_table('Pressure', details.P + ' hPa ' + tendency, table);
      }
    });
  }; 
  var parse_current_weather_items = function(raw_markup_from_feed){
    var result = [];
    var lines = raw_markup_from_feed.trim().split("<br/>");
    var parts;
    for (var i=0; i < lines.length; ++i) {
      parts = lines[i].split(":</b>");
      if (parts.length === 2){
        result.push({
            name: parts[0].trim().substring(3), //chop off leading '<b>'
            value: parts[1].trim() 
          }
        );
      }
    }
    return result;
  };  
  /* Example: http://weather.gc.ca/rss/city/on-118_e.xml */
  var build_canadian_weather_url_rss = function(){
    var prov = input.prov.toLowerCase();
    var weather_stn_id = input.weather_station;
    return 'http://weather.gc.ca/rss/city/' + prov + '-' + weather_stn_id + '_' + trans('e') + '.xml'; 
  };
  /* Example: 'http://weather.gc.ca/city/pages/on-118_metric_e.html'. */
  var canadian_weather_title_url = function(){
    return 'http://weather.gc.ca/city/pages/' + input.prov.toLowerCase() + '-' + input.weather_station + '_metric_' + trans('e') + '.html';
  };
  /* Weather warnings interfere with the regular layout. Can't go by absolute index. Must find where the current-conditions start, and proceed from there. */
  var current_conditions_idx = function(entries){
    var result = 1; //default; the usual value if no warnings are present
    for (var i = 0; i < entries.length; ++i){
      if (entries[i].title.startsWith("Current Conditions:") || entries[i].title.startsWith("Conditions actuelles")){
        result = i;
        break;
      }
    }
    return result;
  }
  var show_current_weather_canada = function() {
    var table = document.getElementById('weather');
    var url = build_canadian_weather_url_rss();
    fetch_xml(UTIL.crossDomainUrl(url+'&ext=xml'), function(err, weather){
      if (err != null){
        console.log('Unable to fetch XML data for Canadian weather rss feed: ' + url);
      }
      else {
        //extract the first CDATA text: '![CDATA[..this stuff...]]'
        var regex = /\!\[CDATA\[([\s\S]*)\]\]/m;  // advice: http://stackoverflow.com/questions/1979884/how-to-use-javascript-regex-over-multiple-lines
        var lines = parse_current_weather_items(regex.exec(weather)[1]); 
        table.innerHTML = table.innerHTML + NL + '<caption><a href="' + canadian_weather_title_url() + '">' + trans('Weather') + '</a>';
        for(var i = 0; i < lines.length; i++){
          if (i===0 || i === lines.length-1) continue;
          table.innerHTML = table.innerHTML + NL + '   <tr><td><b><span>' + lines[i].name + '</span></b><td>' + lines[i].value; 
        }
      }
    });
  };
  var show_current_weather = function() {
    var table = document.getElementById('weather');
    if (input.exclude_current_weather_conditions === '1' || input.country === 'other'){
      dont_show(table);
    }
    else if (input.country === 'us'){
      show_current_weather_us();
    }
    else if (input.country === 'cda' && input.weather_station.length > 0) {
      show_current_weather_canada();
    }
    else if (input.country === 'uk' && input.weather_station.length > 0) {
      show_current_weather_uk();
    }
  };
  show_current_weather();

  var show_weather_forecast_us = function(table) {
    var i;
    fetch_json(us_weather_url('json'), function(err, weather){
      if (err != null){
        console.log('Unable to fetch JSON data for US weather: ' + url);
      }
      else {
        for (i = 0; i < 4; ++i){
          add_to_weather_table(weather.time.startPeriodName[i], weather.data.weather[i], table);
        }
      }
    });
  };
  var parse_forecast = function(title){
    var name = title.split(":")[0].trim();
    var value = '';
    if (title.indexOf(":") > 0){
      //can be missing in the case of weather alerts
      value = title.split(":")[1].trim();
    }
    return {
      name: name,
      value: value
    }; 
  };
  var show_weather_forecast_canada = function(table) {
    var item, match;
    var table = document.getElementById('weather');
    var url = build_canadian_weather_url_rss();
    fetch_xml(UTIL.crossDomainUrl(url+'&ext=xml'), function(err, weather){
      if (err != null){
        console.log('Unable to fetch XML data for Canadian weather rss feed: ' + url);
      }
      else {
        //extract each title AFTER '![CDATA['
        var start = '![CDATA['; 
        var start_idx = weather.indexOf(start) + start.length; //where to start looking for <title> entries
        var raw_forecast = weather.substring(start_idx);
        var regex = /<title>([\s\S]*?)<\/title>/gi; // ? for lazy matching, to avoid getting it all as one title tag
        //var regex_result = regex.exec(raw_forecast);
        var count = 0;
        while (count < 4 && (match = regex.exec(raw_forecast)) != null){ //find n matches in sequence
          item = parse_forecast(match[1]);
          if (item.value){
            table.innerHTML = table.innerHTML + NL + '   <tr><td><b><span>' + item.name + '</span></b><td>' + item.value; 
          } 
          ++count;
        }
      }
    });
  };
  var weather_url_uk_forecast = function(){
    return 'http://datapoint.metoffice.gov.uk/public/data/val/wxfcs/all/json/' + input.weather_station + '?res=daily&key=' + MET_OFFICE_API_KEY;
  };
  var show_weather_forecast_uk = function(table) {
    fetch_json(weather_url_uk_forecast(), function(err, weather){
      if (err != null){
        console.log('Unable to fetch JSON data for UK weather: ' + url);
      }
      else {
        var periods = weather.SiteRep.DV.Location.Period; //array of 5 items, for 5 days
        for(var i = 0; i < 2; ++i){ //take the first 2 days only
          for (var j = 0; j < periods[i].Rep.length; ++j){ //is Rep always an array?? in current-conditions it isn't
            var forecast_when = EPH.when("UT " + periods[i].value.substring(0,10));
            var forecast_weekday = forecast_when.toStringUT(WEEKDAYS_LONG).substring(25);
            var forecast_day_night = periods[i].Rep[j].$ === 'Night' ?  ' night' : '';
            var day = forecast_weekday + forecast_day_night; 
            var forecast = uk_weather_codes()[periods[i].Rep[j].W] + '.'; //translate the code 'W' into text 
            if (periods[i].Rep[j].$ === 'Day'){
              forecast = forecast + ' High ' + periods[i].Rep[j].Dm + '. PoP ' + periods[i].Rep[j].PPd + '%.';
            }
            else {
              forecast = forecast + ' Low ' + periods[i].Rep[j].Nm + '. PoP ' + periods[i].Rep[j].PPn + '%.';
            }
            forecast = forecast + ' Wind ' + periods[i].Rep[j].S + 'mph ' + periods[i].Rep[j].D + '.';
            add_to_weather_table(day, forecast, table);
          }
        }
      }
    });
  };
  var show_weather_forecast = function() {
    var table = document.getElementById('weather');
    if (input.exclude_weather_forecast === '1' || input.country === 'other'){
      //the same table is shared with current conditions
      //do nothing - the current conditions determines if hidden or not, not the forecast
    }
    else {
      if (input.country === 'us'){
        show_weather_forecast_us(table);
      }
      else if (input.country === 'cda'  && input.weather_station.length > 0){
        show_weather_forecast_canada(table);
      }
      else if (input.country === 'uk' && input.weather_station.length > 0) {
        show_weather_forecast_uk(table);
      }
    }
  };
  show_weather_forecast();

  /* Libration chart for given when, and the next n_days.  */  
  var show_libration_and_phase = function(when_start, n_days){
    var when_libration, i, libration, x, y;
    var points = [];
    var TICK_SIZE = 5;
    var NUM_TICKS = 20;
    var ONE_DAY = 60*60*24;
    var canvas = document.getElementById('libration');
    var ctx = canvas.getContext('2d');
    canvas.width = window.innerWidth * 0.18;
    canvas.height = canvas.width;
    standard_canvas_appearance(ctx, true);
    h = canvas.height;
    w = canvas.width;
    if (input.exclude_libration === '1'){
      dont_show(canvas);
    }
    else {
      var current_dist, percent_diff, mean_dist = 384400; //km
      for(i=0; i < 2*n_days; ++i){
        when_libration = when_start.delta(i * ONE_DAY*0.5);
        libration = EPH.lunar_libration(when_libration);
        if (i === 0) {
          current_dist = (libration.distance * 1.495978707E+8) ; //km
          percent_diff = Math.round(100*(current_dist - mean_dist)/mean_dist);
          var sign = percent_diff >= 0 ? '-' : '+';
          GRAPH.text(ctx, sign + Math.abs(percent_diff) + '%',  w*0.86, h*0.93); 
        }   
        points.push({ //wrt the ctr
          x : 0 + EPH.degs(libration.longitude) * w/NUM_TICKS, 
          y : 0 - EPH.degs(libration.latitude) * h/NUM_TICKS
        });
      }
      //draw the result
      ctx.save();
      ctx.translate(w/2, h/2); //use ctr as natural origin of coords
      if (sign_latitude() === -1){
        ctx.scale(-1, -1); // southern hemisphere: flip up-down and right-left
      }
      GRAPH.line(ctx, -w/2, 0, w/2, 0); //x-axis
      GRAPH.line(ctx, 0, -h/2, 0, h/2); //y-axis
      for(i=-NUM_TICKS/2; i <= NUM_TICKS/2; ++i){
        tweak = 1;
        if (i % 5 === 0) tweak = 2;
        GRAPH.tickMarkVertical(ctx,i*(w/NUM_TICKS),0,TICK_SIZE*tweak); 
        GRAPH.tickMarkHorizontal(ctx,0,i*(h/NUM_TICKS),TICK_SIZE*tweak); 
      }
      GRAPH.lines(ctx, points); 
      GRAPH.spot(ctx, points[0].x, points[0].y, 2);
      ctx.restore();
      show_moon_phase(ctx, {x: w/8, y: h/8}, w/10);
    }
  };
  show_libration_and_phase(when, 27);
  
  
  /* 
   Apply the reverse of the projection.
   Return the altitude and azimth (a,A) of the point under the mouse, in degrees. 
   Returns undefined if r is bigger than R.
  */
  var translate_xy_to_aA = function(pos, R /*radius of the planisphere*/, center /*of the planisphere*/){
    var dx = pos.x - center.x;
    var dy = pos.y - center.y;
    var r = Math.sqrt(dx*dx + dy*dy);
    var sign_lat = sign_latitude();
    var result; 
    if (r <= R){
      var rotation = parseInt(input.planisphere_rotation_angle);
      var A = Math.round(EPH.degs(Math.atan2(-sign_lat*dx, -sign_lat*dy)) - rotation);
      A = A % 360;
      if (A < 0) {
        A = A + 360;
      }
      result = {
        a: Math.round(90 - (r/R)*90),
        A: A 
      };
    }
    return result;
  };
  var find_mouse_pos = function(canvas, evt) {
      var rect = canvas.getBoundingClientRect();
      return {
        x: evt.clientX - rect.left,
        y: evt.clientY - rect.top
      };
  };
  /* Ref: http://stackoverflow.com/questions/6735470/get-pixel-color-from-canvas-on-mouseover */ 
  var show_alt_az_tooltip = function(canvas, R, center){
    canvas.addEventListener('mousemove', function(evt) {
      var mouse_pos = find_mouse_pos(canvas, evt);
      var alt_az = translate_xy_to_aA(mouse_pos, R, center);
      canvas.title = 'alt: ' + alt_az.a + ' az:' + alt_az.A;
    }, false);
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
  var dont_draw_past_the_horizon = function(ctx, ctr, R){
      ctx.beginPath();
      ctx.arc(ctr.x,ctr.y,R,0,Math.PI*2,true);
      ctx.clip();      
  };
  var draw_deep_sky_object = function(pos, ctx){
    var gap = 2;
    GRAPH.point(ctx, pos.x+gap, pos.y+gap);
    GRAPH.point(ctx, pos.x+gap, pos.y-gap);
    GRAPH.point(ctx, pos.x-gap, pos.y+gap);
    GRAPH.point(ctx, pos.x-gap, pos.y-gap);
  };
  var draw_deep_sky_objects = function(visible_objects, center, R, ctx){
    var i, proj, dso;
    var limiting_mag = parseFloat(input.limiting_mag_messiers_planisphere);
    for(i=0; i < visible_objects.length; ++i){
      dso = visible_objects[i];
      if (dso.thing[3] <= limiting_mag){
        proj = project_astre(dso, center, R);
        draw_deep_sky_object(proj, ctx);
      }
    }
  };
  var draw_planets = function(planets_etc, center, R, ctx) {
    var planet, proj, i;
    var planets = planets_etc !== null ? planets_etc : planet_rows(true);
    var limiting_mag = parseFloat(input.limiting_mag);
    for(i=0; i < planets.length; ++i){
      planet = planets[i];
      if (planet.ephem.a > 0 && planet.ephem.mag <= limiting_mag){
        proj = project_astre(planet, center, R);
        GRAPH.spot(ctx, proj.x, proj.y, proj.size);
        GRAPH.text(ctx, planet.thing.symbol, proj.x, proj.y+15);
      }
    }
  };
  var star_color = function(size /*0.1..3.5*/){
    var result = 255;
    var frac = (size - 0.1)/3.5;
    if (frac < 0.20) { 
      result = 175; 
    }
    return 'rgb(' + result +',' + result + ',' + result + ')'; //grey-scale
  };
  var draw_stars = function(star_poss, ctx) {
    var i, star;
    //star_poss is sparse; its indices point to the index in YBS; we need to iterate like this:
    for (i in star_poss){
      if (star_poss.hasOwnProperty(i)){
        star = star_poss[i];
        if (star.size > 0){
          ctx.save();
          var star_col = star_color(star.size)
          ctx.fillStyle = star_col;
          ctx.strokeStyle = star_col;
          GRAPH.spot(ctx, star.x, star.y, star.size);
          ctx.restore();
        }
      }
    }
  };
  var draw_comet = function(pos, ctx){
    GRAPH.circle(ctx,pos.x, pos.y,3);
  };
  var draw_comets = function(center, R, ctx){
     var comets = EPH.as_array(EPH.comets);
     for(var i=0; i < comets.length; ++i){
        var ephem = EPH.position(comets[i].name, when, opts);
        if (ephem.a > 0){
          var proj = project_astre({ephem:ephem}, center, R);
          draw_comet(proj, ctx);
        }
     }
  };
  var draw_meteor_shower_radiant = function(pos, ctx){
    var s = 4;
    GRAPH.line(ctx, pos.x, pos.y-s, pos.x, pos.y+s);
    GRAPH.line(ctx, pos.x-s, pos.y, pos.x+s, pos.y);
  };
  var draw_meteor_shower_radiants = function(center, R, ctx){
     var showers = EPH.current_meteor_showers(when, opts);
     for(var i=0; i < showers.length; ++i){
        if (showers[i].radiant.a > 0){
          var ephem = {a:showers[i].radiant.a , A:showers[i].radiant.A};
          var proj = project_astre({ephem:ephem}, center, R);
          draw_meteor_shower_radiant(proj, ctx);
        }
     }
  };
  var draw_symbol_legend = function(item, x, y, ctx){
    item.drawer({x:x, y:y}, ctx);
    GRAPH.text(ctx, item.title, x+10, y);
  };
  var draw_symbol_legends = function(ctx){
    var x = 20; // top left corner
    var y = 20;
    var items = [
      {title: trans('Messier'), drawer: draw_deep_sky_object},
      {title: trans('Comet'), drawer: draw_comet},
      {title: trans('Radiant'), drawer: draw_meteor_shower_radiant}
    ];
    ctx.save();
    ctx.textAlign="left";
    for (var i = 0; i < items.length; ++i){
      draw_symbol_legend(items[i], x, y, ctx);
      y = y + 15; 
    }
    ctx.restore();
  };
  /* Join up sets of star positions, to show the usual lines denoting a constellation.  */
  var draw_constellation_lines = function(star_poss, ctx){
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
            if (star_poss[polyline[j]]){
              points.push(star_poss[polyline[j]]);
            }
          }
          if (points.length >= 2){
            draw_polyline(ctx, points);
          }
        }
      }
    }
  };
  var project_astre = function(astre /*.ephem (with .mag) */, center, R){
    var r = (90 - astre.ephem.a) * (R/90);
    //to rotate the planisphere about the zenith, just add an angle to theta:
    //'polar azimuthal equidistant projection': simple r-theta
    var rotation = parseInt(input.planisphere_rotation_angle);
    var theta = EPH.rads(astre.ephem.A) + EPH.rads(rotation);
    var sign_lat = sign_latitude();  
    var x = center.x - sign_lat * r * Math.sin(theta);
    var y = center.y - sign_lat * r * Math.cos(theta);
    var size = star_size(astre.ephem.mag);
    return {
      x: x,
      y: y,
      size: size
    };
  };
  var project_star_positions = function(visible_stars, center, R){
    var i, projection;
    var result = [];
    for(i=0; i < visible_stars.length; ++i){
      var projection = project_astre(visible_stars[i], center, R);
      //the array is sparse; the indices are not continuous
      result[visible_stars[i].idx] = projection; //array idx must id the star, in order to draw the constellation lines later  
    }
    return result;
  };
  var draw_the_horizon = function(ctx, center, R){
      GRAPH.circle(ctx,center.x,center.y,R); 
  };
  var blue_background = function(ctx, center, R){
    ctx.save();
    ctx.fillStyle = 'rgb(36,55,114)';
    ctx.beginPath();
    ctx.arc(center.x, center.y, R, 0, Math.PI*2, false);
    ctx.fill();
    ctx.closePath();
    ctx.restore();
  };
  var show_planisphere = function(visible_stars, planets_etc){
    var h, w, i, star, R, center, x, y, r, theta, planets, planet, messier;
    var limiting_mag = parseFloat(input.limiting_mag);
    var canvas = document.getElementById('planisphere');
    var ctx = canvas.getContext('2d');
    canvas.width = window.innerWidth * 0.95;     
    canvas.height = canvas.width;     
    standard_canvas_appearance(ctx, true);
    h = canvas.height;
    w = canvas.width;
    if (input.exclude_planisphere === '1'){
      dont_show(canvas);
    }
    else {
      R = w/2 - 4; // radius of the planisphere
      center = {x: w/2, y:h/2};
      blue_background(ctx, center, R);
      draw_symbol_legends(ctx);
      draw_the_horizon(ctx,center,R); 
      dont_draw_past_the_horizon(ctx, center, R);
      var star_poss = project_star_positions(visible_stars, center, R);
      //draw the constellation lines first, otherwise the lines go inside the star-spots
      draw_constellation_lines(star_poss, ctx);
      draw_stars(star_poss, ctx);
      draw_planets(planets_etc, center, R, ctx);
      //draw_messiers(visible_messiers, center, R, ctx);
      draw_deep_sky_objects(visible_messiers, center, R, ctx);      
      draw_deep_sky_objects(visible_caldwells, center, R, ctx);      
      draw_comets(center, R, ctx);
      draw_meteor_shower_radiants(center, R, ctx);      
      show_alt_az_tooltip(canvas, R, center);
    }
  };
  var visible_stars = EPH.find_visible_stars(when, opts.where, -20, parseFloat(input.limiting_mag), opts);
  show_planisphere(visible_stars, planets_etc);
  
  
  /* The whole point here is simply to maniplate the date-time of the URL. */
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
  
  var min_aurora_kp_level = function (geomagnetic_latitude){
    var result = 0;
    if (input.aurora_min_activity_level === '-1') {
      result = EPH.aurora_min_kp(geomagnetic_latitude);
    }
    else {
      result = parseInt(input.aurora_min_activity_level);
    }
    return result;
  };  
  var max_aurora_level_in_current = function(kp_array){
    var i, result=0, this_level;
    for (i=0; i < kp_array.length; i++){
      this_level = kp_array[i].value;
      if (this_level > result){
        result = this_level;
      }
    }
    console.log('Max auroral activity level: ' + result);
    return result;
  };
  var where_aurora = function(){
    return {
      λ : parseFloat(input.longitude), 
      φ : parseFloat(input.latitude)
    };
  };
  var show_aurora_kp_measurements = function(kp_array, where){
    var when, when_text, value, what_to_show, row;
    kp_array.reverse();
    var table = document.getElementById('aurora_table');
    if (input.exclude_aurora === '1'){
      dont_show(table);
    }
    else {
      var geomagnetic_latitude = EPH.round(EPH.geomagnetic_latitude(where_aurora(), EPH.when_now()),1);
      var min_kp_level = min_aurora_kp_level(geomagnetic_latitude);
      if (max_aurora_level_in_current(kp_array) >= min_kp_level){
        for (var i=0; i < kp_array.length; i++){
          when = EPH.when(kp_array[i].datetime);
          when_text = trans(when.toStringLT(WEEKDAYS));
          row = '   <tr><td>' + when_text.substring(0,19);
          value = kp_array[i].value;
          for (var j = 0; j < 10; ++j){
            what_to_show = (value < j ? '-' : 'x');
            row = row + '<td>' + what_to_show; 
          }
          table.innerHTML = table.innerHTML + NL + row;
        }
      }
      else {
        dont_show(table);
      }
    }
  };
  var show_aurora_data = function(){
    if (input.exclude_aurora !== '1'){
      fetch_aurora_data(show_aurora_kp_measurements);
    }
  };
  show_aurora_data();
  
  var show_occultations = function(){
    var table = document.getElementById('occultation_table');
    var when_occn, occn, rounded_time, look_ahead, min_mag, link, when_text, moon;
    if (input.exclude_occultations !== '1'){
      look_ahead = parseInt(input.occultations_num_days_ahead,10);
      min_mag = parseFloat(input.occultations_min_mag);
      var occns = fetch_occultation_predictions(where_aurora(), when, look_ahead, min_mag);
      if (occns.length > 0){
        for (var i=0; i < occns.length; i++){
          // {UT:'2016-10-20 21:08:18', star:'ZC 123', mag:5.5, ph:'DD', el:92, pa:145}
          occn = occns[i];
          rounded_time = nearest_minute(occn.UT);
          when_occn = EPH.when('UT ' + rounded_time);
          moon = EPH.position('moon', when_occn, opts);
          when_text = when_occn.toStringLT(WEEKDAYS).substring(0,19); 
          link = chart_link_for_date_time('moon', when_text, when_text.substring(3), 'LT');
          row = '   <tr style="text-align:right;"><td>' + link + '<td>' + occn.star + '<td>' + occn.mag + '<td>' + trans(occn.ph) + '<td>' + occn.el + '°<td>' + occn.pa + '°' + '<td>' + Math.round(moon.a) + '°';
          //row = '   <tr style="text-align:right;"><td>' + when_occn.toStringLT(WEEKDAYS).substring(0,19) + '<td>' + occn.star + '<td>' + occn.mag + '<td>' + trans(occn.ph) + '<td>' + occn.el + '°<td>' + occn.pa + '°' + '<td>' + EPH.round(moon.alt, 1) + '°';
          table.innerHTML = table.innerHTML + NL + row;
        }
      }
      else {
        dont_show(table);
      }
    }
    else {
      dont_show(table);
    }
  };
  show_occultations();

  var ph_type = function(ph/*Sh.I*/){
    var parts = ph.split("\.");
    var abbr = parts[0];
    var long_text = {
      'Sh': 'Shadow transit',
      'Tr': 'Satellite transit',
      'Oc': 'Occultation',
      'Ec': 'Eclipse'
    };
    return trans(long_text[abbr]);
  };  
  var ph_action = function(ph/*Sh.I*/){
    var parts = ph.split("\.");
    var abbr = parts[1];
    var long_text = {
      'I': 'Ingress',
      'E': 'Egress',
      'D': 'Disappearance',
      'R': 'Reappearance'
    };
    return trans(long_text[abbr]);
  };  
  var show_jupiter_satellite_phenomena = function(){
    var table = document.getElementById('jupiter_satellite_phenomena');
    var satts = ['Io', 'Europa', 'Ganymede', 'Callisto'];
    var events, event, obs_window, when_text;
    if (input.exclude_jupiter_satellite_phenomena !== '1'){
      obs_window = EPH.observation_window(when, opts_rads, parseInt(input.twilight)); 
      events = fetch_jupiter_satellite_phenomena(obs_window.sunset, obs_window.sunrise, opts.where);
      if (events.length > 0){
        for (var i=0; i < events.length; i++){
          event = events[i];
          when_text = EPH.when(event.when).toStringLT().substring(0,22);
          table.innerHTML  = table.innerHTML + NL + '  <tr><td>' + when_text + '<td>' + trans(satts[event.satt]) + '<td>' + ph_type(event.ph) + ', ' + ph_action(event.ph) + NL;
        }
      }
      else {
        dont_show(table);
      }
    }
    else {
      dont_show(table);
    }
  };
  show_jupiter_satellite_phenomena();
  
  var url_to_large_satellite_image = function(){
    var result = '../satellite/graphic.sky?degrees_on_a_side=10&pixels_on_a_side=800';
    result = result + '&latitude=' + UTIL.escapeHtml(input.latitude);
    result = result + '&longitude=' + UTIL.escapeHtml(input.longitude);
    result = result + '&layer=' + UTIL.escapeHtml(input.layer);
    return result;
  };
  var show_a_larger_satellite_image = function(){
    window.location.href = url_to_large_satellite_image(); 
  };
  var link_to_large_satellite_image = function(){
    if (input.country === 'cda' || input.country === 'us'){
      document.getElementById('clouds').onclick = show_a_larger_satellite_image;
    }
  };
  link_to_large_satellite_image();
  
};