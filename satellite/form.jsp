<!doctype html>
<html lang='<s:txt>en</s:txt>'>
<head>
 <tags:analytics/>
 <meta charset='UTF-8'>
 <meta name="keywords" content="<s:txt>canada, united states, weather, clouds, satellite, goes, visible, infrared, imagery</s:txt>">
 <meta name="description" content="<s:txt>Recent satellite images from the GOES weather satellite (NOAA).</s:txt>">
 <meta name="viewport" content="width=device-width">
 <link rel="stylesheet" href="../css/styles.css<tags:ver/>" media="all">
 <title><s:txt>Satellite images</s:txt></title>
</head>
<script>
  <tags:geolocation/>
  window.onload = function() {
    activate_lat_long_autofill();
  };
  <tags:screen_change/>
</script>
<body>
<h2><s:txt>Satellite images</s:txt> 
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
 <s:txt>The form below lets you browse current imagery from the GOES geosynchronous satellites, for North America.
 These images are especially dramatic when the Sun is low in the sky, and produces longer shadows. 
 The default settings in the form below show the Canadian Maritime provinces.</s:txt> 
 
  <form>
   <select onChange="if (this.value) window.location.href=this.value">
    <tags:locations_sat_image/>
   </select>
  </form>
 
 <P>
 <form method='GET' action='graphic.sky' class='user-input-small' id='input_form_clouds'>
  <table>
     <tr><td><s:txt>Latitude</s:txt>:<td><input type='text' name='latitude' id='latitude' value='47' title='<s:txt>Latitude in degrees</s:txt>'> 
     <tr><td><s:txt>Longitude</s:txt>:<td><input type='text' name='longitude' id='longitude' value='-64' title='<s:txt>Longitude in degrees. Negative west of Greenwich</s:txt>'>
     <tr><td><td><button type='button' id='lat_long_autofill' title='<s:txt>Let the browser fill in lat/long</s:txt>'><s:txt>Current location</s:txt></button> 
     <tr><td><s:txt>Degrees on a side</s:txt><td><input name='degrees_on_a_side' value='10'>
     <tr><td><s:txt>Pixels on a side</s:txt><td><input id='pixels_on_a_side' name='pixels_on_a_side' value=''>
     <tr><td><s:txt>Channel</s:txt>
         <td>
         <select name='layer'>
           <option value='auto_detect' title='<s:txt>Let the system decide</s:txt>' selected><s:txt>Auto-detect</s:txt>
           <option value='visible' title='<s:txt>Day time</s:txt>'><s:txt>Visible (day-time)</s:txt>
           <option value='ir' title='<s:txt>Night time</s:txt>'><s:txt>IR (night-time)</s:txt>
         </select>
     <tr><td><s:txt>Show locations</s:txt><td><input name='locations' id='locations' title='<s:txt>Latitude and longitude; separate with a semi-colon</s:txt>' value='46.25, -63.13' size='20'>
     <tr><td><td><input type='submit' value='<s:txt>Show satellite image</s:txt>'>
  </table>
 </form>

<P><s:txt>Typically, you will use the visible channel during the day, and the infrared channel at night.</s:txt>

<P><s:txt>Credits</s:txt>:
<ul>
 <li><a href='http://radar.weather.gov/'>NOAA/NWS</a>: <s:txt>imagery from the GOES geostationary weather satellite</s:txt>.
 <!--
 <li><a href='https://mesonet.agron.iastate.edu/'>Iowa State University</a>: <s:txt>server for NOAA images captured by the GOES satellite</s:txt>.
 <li><a href='http://www.eumetsat.int/website/home/Satellites/CurrentSatellites/Meteosat/index.html'>EUMETSAT</a>: <s:txt>imagery from the Meteosat geostationary weather satellite</s:txt>.
 -->
</ul>

</body>
</html>