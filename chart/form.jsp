<!doctype html>
<html lang='en'>
<head>
 <meta charset='UTF-8'>
 <meta name="keywords" content="astronomy">
 <meta name="description" content="Star charts for amateur astronomers.">
 <meta name="viewport" content="width=device-width">
 <link rel="stylesheet" href="../css/styles.css<tags:ver/>" media="all">
 <title>Star chart</title>
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
<h2>Star chart <a href="../main/form.sky" class='home'>Home</a></h2>


<P>
 <form method='GET' action='graphic.sky' class='user-input-small' id='input_form'>
  <table>
     <tr><td>Center on object*:<td>
         <select name='center_object' id='center_object' title='Name of the object. Overrides RA, Dec setting.' autofocus>
           <option>  
         </select>
     <tr><td>Center on RA*:<td><input type='text' name='center_ra' id='center_ra' title='Right ascension in degrees' value='0' > 
     <tr><td>Center on Dec*:<td><input type='text' name='center_dec' id='center_dec' title='Declination in degrees' value='0' >
     <tr><td title='Translation is not complete.'>Preferred language:<td>
         <select name='preferred_lang' id='preferred_lang' title='Output in English or French'>
           <option value='e' selected>English
           <option value='f'>French 
         </select>
     <tr><td>Place at top:<td>
         <select name='place_at_top' id='place_at_top'>
           <option value='north' selected>North
           <option value='south'>South
         </select>
     <tr><td>Date-time (defaults to now):<td><input type='text' name='date_time' id='date_time' title='Example: 2016-03-01 18:00. Leave empty to use the current date and time.' >
     <tr><td>Time scale:<td><select name='time_scale' id='time_scale'><option title='Local Time (defined by your browser)' value='LT'>Local Time<option title='Universal Time (UT1, to be precise)' value='UT'>Universal Time<option title='Terrestrial Time (fundamental physics time)' value='TT'>Terrestrial Time</select>
     <tr><td>Location:<td><input type='text' name='location_name' id='location_name' value='Ottawa' title='Name of your observing site'>
     <tr><td>Latitude:<td><input type='text' id='latitude' name='latitude' value='45.403' title='Latitude in degrees'>
        <td><button type='button' id='lat_long_autofill' title='Let the browser fill in lat/long'>Use current location</button>
     <tr><td>Longitude:<td><input type='text' id='longitude' name='longitude' id='longitude' value='-75.664' title='Longitude in degrees. Negative west of Greenwich'>
     
     <tr><td>Limiting stellar magnitude:<td><input type='text' name='limiting_mag' id='limiting_mag' value='8.0' title='Higher number means show dimmer stars'>
     <tr><td>Chart width in degrees:<td><input type='text' name='chart_width' id='chart_width'value='30'>
     <tr><td>Style:<td>
         <select name='style' id='style' >
           <option value='regular' selected>Regular
           <option value='print_friendly'>Print-friendly 
         </select>
     <tr><td>Show identifiers (multi-select):<td>
         <select name='show_identifiers' id='show_identifiers' multiple='multiple'>
           <option value='bayer' selected>Bayer
           <option value='flamsteed' selected>Flamsteed 
         </select>
     <tr>
       <td>
       <td><input type='submit' value='Show the star chart'>
       <td><a href='form.sky'><button type='button' id='revert_to_default'>Default settings</button></a>
  </table>
 </form>

<P>Note:
<ul>
 <li>all items are optional, except for <em>center on</em>. You can either choose a specific object, or enter a specific 
 right ascension and declination.
 <li>the location, latitude, and longitude are used only for lunar parallax. 
 <li>the other items show default values. 
  Usually, just accepting the defaults is adequate.
</ul> 
 
<tags:analytics/>

</html>