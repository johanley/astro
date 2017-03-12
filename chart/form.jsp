<!doctype html>
<html lang='<s:txt>en</s:txt>'>
<head>
 <meta charset='UTF-8'>
 <meta name="keywords" content="<s:txt>astronomy, custom, star chart</s:txt>">
 <meta name="description" content="<s:txt>Star charts for amateur astronomers.</s:txt>">
 <meta name="viewport" content="width=device-width">
 <link rel="stylesheet" href="../css/styles.css<tags:ver/>" media="all">
 <title><s:txt>Star chart</s:txt></title>
 <script src='../js/ephem.js<tags:ver/>'></script>
 <script src='../js/util.js<tags:ver/>'></script>
 <script>
  var add_option = function(opt_name){
    var select = document.getElementById('center_object');
    var option = document.createElement('option');
    option.text = opt_name; 
    select.appendChild(option);        
  };
  var populate_options_from_array = function(array){
    for(var i = 0; i < array.length; ++i){
      add_option(array[i][0]);
    }
  };
  var populate_options_from_obj = function(obj){
    for (var prop in obj){
      if (obj.hasOwnProperty(prop)){
        add_option(obj[prop].name);
      }
    };
  };
  var populate_names_of_known_objects = function(){
    populate_options_from_obj(EPH.planets);
    populate_options_from_obj(EPH.comets);
    populate_options_from_obj(EPH.minor_planets);
    populate_options_from_array(EPH.messiers);
    populate_options_from_array(EPH.caldwells);
  };
  
  var activate_lat_long_autofill = function(){
     var lat_long_auto = document.getElementById('lat_long_autofill');
     lat_long_auto.onclick = function(){
       if (navigator.geolocation){
         navigator.geolocation.getCurrentPosition(
           function(pos){
             document.getElementById('latitude').value = pos.coords.latitude;
             document.getElementById('longitude').value = pos.coords.longitude;
           }
         );
       }
       else {
         lat_long_auto.innerHTML = 'Disabled by browser';
       }
     };
  };
  
  window.onload = function() {
    populate_names_of_known_objects();
    activate_lat_long_autofill();
    UTIL.prepopulate_form_controls(['show_identifiers']);
  };
 </script>
</head>
<body>
<h2><s:txt>Star chart</s:txt> 
 <span class='home'>
  <form style='display:inline;'>
   <select onChange="if (this.value) window.location.href=this.value" style='vertical-align:text-top;'>
    <option value='form.sky'>Language
     <option value='form.sky?lang=en'>English
     <option value='form.sky?lang=fr'>Français
   </select>
  </form>
 <a href='../main/form.sky'><s:txt>Home</s:txt></a>
 </span>
</h2>

<P>
 <form method='GET' action='graphic.sky' class='user-input-small' id='input_form'>
  <table>
     <tr><td><s:txt>Center on object</s:txt>*:<td>
         <select name='center_object' id='center_object' title='<s:txt>Name of the object. Overrides RA, Dec setting.</s:txt>' autofocus>
           <option>  
         </select>
     <tr><td><s:txt>Center on RA</s:txt>*:<td><input type='text' name='center_ra' id='center_ra' title='<s:txt>Right ascension in degrees</s:txt>' value='0' > 
     <tr><td><s:txt>Center on Dec</s:txt>*:<td><input type='text' name='center_dec' id='center_dec' title='<s:txt>Declination in degrees</s:txt>' value='0' >
     <tr><td><s:txt>Place at top</s:txt>:<td>
         <select name='place_at_top' id='place_at_top'>
           <option value='north' selected><s:txt>North</s:txt>
           <option value='south'><s:txt>South</s:txt>
         </select>
     <tr><td><s:txt>Date-time (defaults to now)</s:txt>:<td><input type='text' name='date_time' id='date_time' title='<s:txt>Example: 2016-03-01 18:00. Leave empty to use the current date and time.</s:txt>' >
     <tr><td><s:txt>Time scale</s:txt>:<td>
       <select name='time_scale' id='time_scale'>
         <option title='<s:txt>Local Time (defined by your browser)</s:txt>' value='LT'><s:txt>Local Time</s:txt>
         <option title='<s:txt>Universal Time (UT1, to be precise)</s:txt>' value='UT'><s:txt>Universal Time</s:txt>
         <option title='<s:txt>Terrestrial Time (fundamental physics time)</s:txt>' value='TT'><s:txt>Terrestrial Time</s:txt>
       </select>
     <tr><td><s:txt>Location</s:txt>:<td><input type='text' name='location_name' id='location_name' value='Ottawa' title='<s:txt>Name of your observing site</s:txt>'>
     <tr><td><s:txt>Latitude</s:txt>:<td><input type='text' id='latitude' name='latitude' value='45.403' title='<s:txt>Latitude in degrees</s:txt>'>
        <td><button type='button' id='lat_long_autofill' title='<s:txt>Let the browser fill in lat/long</s:txt>'><s:txt>Use current location</s:txt></button>
     <tr><td><s:txt>Longitude</s:txt>:<td><input type='text' id='longitude' name='longitude' id='longitude' value='-75.664' title='<s:txt>Longitude in degrees. Negative west of Greenwich</s:txt>'>
     
     <tr><td><s:txt>Limiting stellar magnitude</s:txt>:<td><input type='text' name='limiting_mag' id='limiting_mag' value='8.0' title='<s:txt>Higher number means show dimmer stars</s:txt>'>
     <tr><td><s:txt>Chart width in degrees</s:txt>:<td><input type='text' name='chart_width' id='chart_width'value='30'>
     <tr><td><s:txt>Style</s:txt>:<td>
         <select name='style' id='style' >
           <option value='regular' selected><s:txt>Regular</s:txt>
           <option value='print_friendly'><s:txt>Print-friendly</s:txt> 
         </select>
     <tr><td><s:txt>Show identifiers (multi-select)</s:txt>:<td>
         <select name='show_identifiers' id='show_identifiers' multiple='multiple'>
           <option value='bayer' selected>Bayer
           <option value='flamsteed' selected>Flamsteed 
         </select>
     <tr>
       <td>
       <td><input type='submit' value='<s:txt>Show the star chart</s:txt>'>
       <td><a href='form.sky'><button type='button' id='revert_to_default'><s:txt>Default settings</s:txt></button></a>
  </table>
 </form>

<P><s:txt>Note</s:txt>:
<ul>
 <li><s:txt>all items are optional, except for <em>center on</em>. You can either choose a specific object, or enter a specific right ascension and declination.</s:txt>
 <li><s:txt>the location, latitude, and longitude are used only for lunar parallax.</s:txt> 
 <li><s:txt>the other items show default values. Usually, just accepting the defaults is adequate.</s:txt>
</ul> 
 
<tags:analytics/>

</html>