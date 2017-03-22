<!doctype html>
<html lang='<s:txt>en</s:txt>'>
<head>
 <meta charset='UTF-8'>
 <meta name="keywords" content="<s:txt>canada, united states, weather, clouds, satellite, goes, visible, infrared, imagery</s:txt>">
 <meta name="description" content="<s:txt>Recent satellite images from the GOES weather satellite (NOAA).</s:txt>">
 <meta name="viewport" content="width=device-width">
 <link rel="stylesheet" href="../css/styles.css<tags:ver/>" media="all">
 <title><s:txt>Satellite images</s:txt></title>
</head>
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
 <s:txt>The form below lets you browse current imagery from the GOES and Meteosat geosynchronous satellites.
 GOES has images for North America, and Meteosat has images for Europe and Africa.
 These images are especially dramatic when the Sun is low in the sky, and produces longer shadows. 
 The default settings in the form below show the Canadian Maritime provinces.</s:txt> 
 
  <form>
   <select onChange="if (this.value) window.location.href=this.value">
     <option value=''><s:txt>Choose a location</s:txt>
     <option value=''><s:txt>---Canada---</s:txt>
     <option value='graphic.sky?layer=auto_detect&degrees_on_a_side=15&pixels_on_a_side=800&latitude=52.00&longitude=-121.00'><s:txt>BC, AB</s:txt>
     <option value='graphic.sky?layer=auto_detect&degrees_on_a_side=15&pixels_on_a_side=800&latitude=53.00&longitude=-102.00'><s:txt>SK, MN</s:txt>
     <option value='graphic.sky?layer=auto_detect&degrees_on_a_side=15&pixels_on_a_side=800&latitude=48.00&longitude=-81.80'><s:txt>Ontario</s:txt>
     <option value='graphic.sky?layer=auto_detect&degrees_on_a_side=15&pixels_on_a_side=800&latitude=49.30&longitude=-72.10'><s:txt>Québec</s:txt>
     <option value='graphic.sky?layer=auto_detect&degrees_on_a_side=15&pixels_on_a_side=800&latitude=48.40&longitude=-61.80'><s:txt>East Coast</s:txt>
     <option value=''><s:txt>---U.S.A.---</s:txt>
     <option value='graphic.sky?layer=auto_detect&degrees_on_a_side=15&pixels_on_a_side=800&latitude=42.00&longitude=-75.00'><s:txt>New England</s:txt>
     <option value='graphic.sky?layer=auto_detect&degrees_on_a_side=15&pixels_on_a_side=800&latitude=32.00&longitude=-83.00'><s:txt>Florida</s:txt>
     <option value='graphic.sky?layer=auto_detect&degrees_on_a_side=15&pixels_on_a_side=800&latitude=37.60&longitude=-118.50'><s:txt>California</s:txt>
     <option value='graphic.sky?layer=auto_detect&degrees_on_a_side=15&pixels_on_a_side=800&latitude=44.60&longitude=-118.41'><s:txt>Oregon</s:txt>
     <option value='graphic.sky?layer=auto_detect&degrees_on_a_side=15&pixels_on_a_side=800&latitude=42.92&longitude=-95.03'><s:txt>Iowa</s:txt>
     <option value='graphic.sky?layer=auto_detect&degrees_on_a_side=15&pixels_on_a_side=800&latitude=40.27&longitude=-107.98'><s:txt>Colorado</s:txt>
     <option value=''><s:txt>---Europe---</s:txt>
     <option value='graphic.sky?layer=auto_detect&degrees_on_a_side=15&pixels_on_a_side=800&latitude=54.70&longitude=-3.90'><s:txt>United Kingdom</s:txt>
     <option value='graphic.sky?layer=auto_detect&degrees_on_a_side=15&pixels_on_a_side=800&latitude=42.30&longitude=13.30'><s:txt>Italy</s:txt>
     <option value='graphic.sky?layer=auto_detect&degrees_on_a_side=15&pixels_on_a_side=800&latitude=38.71&longitude=22.26'><s:txt>Greece</s:txt>
     <option value='graphic.sky?layer=auto_detect&degrees_on_a_side=15&pixels_on_a_side=800&latitude=50.93&longitude=10.53'><s:txt>Germany</s:txt>
     <option value='graphic.sky?layer=auto_detect&degrees_on_a_side=15&pixels_on_a_side=800&latitude=47.32&longitude=2.49'><s:txt>France</s:txt>
     <option value='graphic.sky?layer=auto_detect&degrees_on_a_side=15&pixels_on_a_side=800&latitude=40.30&longitude=-3.75'><s:txt>Spain</s:txt>
     <option value=''><s:txt>---Africa---</s:txt>
     <option value='graphic.sky?layer=auto_detect&degrees_on_a_side=15&pixels_on_a_side=800&latitude=26.50&longitude=30.48'><s:txt>Egypt</s:txt>
     <option value='graphic.sky?layer=auto_detect&degrees_on_a_side=15&pixels_on_a_side=800&latitude=-26.79&longitude=24.91'><s:txt>South Africa</s:txt>
     <option value='graphic.sky?layer=auto_detect&degrees_on_a_side=15&pixels_on_a_side=800&latitude=26.38&longitude=18.14'><s:txt>Libya</s:txt>
     <option value='graphic.sky?layer=auto_detect&degrees_on_a_side=15&pixels_on_a_side=800&latitude=-1.43&longitude=23.70'><s:txt>Congo</s:txt>
     <option value='graphic.sky?layer=auto_detect&degrees_on_a_side=15&pixels_on_a_side=800&latitude=12.10&longitude=7.90'><s:txt>Nigeria</s:txt>
     <option value='graphic.sky?layer=auto_detect&degrees_on_a_side=15&pixels_on_a_side=800&latitude=6.70&longitude=40.98'><s:txt>Ethiopia</s:txt>
     <option value='graphic.sky?layer=auto_detect&degrees_on_a_side=15&pixels_on_a_side=800&latitude=30.20&longitude=-5.85'><s:txt>Morocco</s:txt>
     <option value=''><s:txt>---Middle East---</s:txt>
     <option value='graphic.sky?layer=auto_detect&degrees_on_a_side=15&pixels_on_a_side=800&latitude=34.39&longitude=42.11'><s:txt>Iraq</s:txt>
     <option value='graphic.sky?layer=auto_detect&degrees_on_a_side=15&pixels_on_a_side=800&latitude=25.60&longitude=43.40'><s:txt>Saudi Arabia</s:txt>
   </select>
  </form>
 
 <P>
 <form method='GET' action='graphic.sky' class='user-input-small' id='input_form_clouds'>
  <table>
     <tr><td><s:txt>Latitude</s:txt>:<td><input type='text' name='latitude' value='47' title='<s:txt>Latitude in degrees</s:txt>'> 
     <tr><td><s:txt>Longitude</s:txt>:<td><input type='text' name='longitude' value='-64' title='<s:txt>Longitude in degrees. Negative west of Greenwich</s:txt>'>
     <tr><td><s:txt>Degrees on a side</s:txt><td><input name='degrees_on_a_side' value='10'>
     <tr><td><s:txt>Pixels on a side</s:txt><td><input name='pixels_on_a_side' value='800'>
     <tr><td><s:txt>Channel: visible or infrared?</s:txt>
         <td>
         <select name='layer'>
           <option value='auto_detect' title='<s:txt>Let the system decide</s:txt>' selected><s:txt>Auto-detect</s:txt>
           <option value='visible' title='<s:txt>Day time</s:txt>'><s:txt>Visible (day-time)</s:txt>
           <option value='ir' title='<s:txt>Night time</s:txt>'><s:txt>IR (night-time)</s:txt>
         </select>
     <tr><td><s:txt>Show locations</s:txt><td><input name='locations' id='locations' title='<s:txt>Latitude and longitude; separate with a semi-colon</s:txt>' value='46.25, -63.13' size='40'>
     <tr><td colspan='2' style="text-align:center">
      <input type='submit' value='<s:txt>Show the large-scale satellite image (clouds)</s:txt>'>
  </table>
 </form>

<P><s:txt>Typically, you will use the visible channel during the day, and the infrared channel at night.</s:txt>

<P><s:txt>Credits</s:txt>:
<ul>
 <li><a href='http://radar.weather.gov/'>NOAA/NWS</a>: <s:txt>imagery from the GOES geostationary weather satellite</s:txt>.
 <li><a href='https://mesonet.agron.iastate.edu/'>Iowa State University</a>: <s:txt>server for NOAA images captured by the GOES satellite</s:txt>.
 <li><a href='http://www.eumetsat.int/website/home/Satellites/CurrentSatellites/Meteosat/index.html'>EUMETSAT</a>: <s:txt>imagery from the Meteosat geostationary weather satellite</s:txt>.
</ul>

<tags:analytics/>

</body>
</html>