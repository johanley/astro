<!doctype html>
<html lang='en'>
<head>
 <meta charset='UTF-8'>
 <meta name="keywords" content="canada, united states, weather, clouds, satellite, goes, visible, infrared, imagery">
 <meta name="description" content="Recent satellite images from the GOES weather satellite (NOAA).">
 <meta name="viewport" content="width=device-width">
 <link rel="stylesheet" href="../css/styles.css<tags:ver/>" media="all">
 <title>GOES satellite images</title>
</head>
<body>
<h2>GOES satellite images <a href='../main/form.sky' class='home'>Home</a></h2>

 <P>
 The form below lets you browse current imagery from the GOES geosynchronous satellite. 
 These images are especially dramatic when the Sun is low in the sky, and produces longer shadows. 
 The default settings in the form below show the Maritimes. 
 
 <form method='GET' action='graphic.sky' class='user-input-small' id='input_form_clouds'>
  <table>
     <tr><td>Latitude:<td><input type='text' name='latitude' value='47' title='Latitude in degrees'> 
     <tr><td>Longitude:<td><input type='text' name='longitude' value='-64' title='Longitude in degrees. Negative west of Greenwich'>
     <tr><td>Degrees on a side<td><input name='degrees_on_a_side' value='10'>
     <tr><td>Pixels on a side<td><input name='pixels_on_a_side' value='800'>
     <tr><td>Channel: visible or infrared?
         <td>
         <select name='layer'>
           <option value='auto_detect' title='Let the system decide' selected>Auto-detect
           <option value='visible' title='Day time'>Visible (day-time)
           <option value='ir' title='Night time'>IR (night-time)
         </select>
     <tr><td colspan='2' style="text-align:center">
      <input type='submit' value='Show the large-scale satellite image (clouds)'>
  </table>
 </form>

<P>Typically, you will use the visible channel during the day, and the infrared channel at night.

<P>Credits:
<ul>
 <li><a href='http://radar.weather.gov/'>NOAA/NWS</a> : imagery from the GOES geostationary weather satellite.
 <li><a href='https://mesonet.agron.iastate.edu/'>Iowa State University</a>: server for NOAA images captured by the GOES satellite.
</ul>

<tags:analytics/>

</body>
</html>