<!doctype html>
<html lang='en'>
<head>
 <meta charset='UTF-8'>
 <meta name="keywords" content="astronomy">
 <meta name="description" content="Astronomical ephemerides for amateur astronomers">
 <meta name="viewport" content="width=device-width">
 <link rel="stylesheet" href="../css/styles.css<tags:ver/>" media="all">
 <style>
  table.report {
    display: inline-block;
    vertical-align: top; 
  }
  body {
    margin:0.5em;
    padding: 0;
    font: 0.75em Verdana, Arial, Helvetica, sans-serif;
  }
  canvas, img {
    border-radius:10px; 
    background-color: rgb(25%, 25%, 25%);
    border-color: rgb(45%,45%,45%); 
    border-style: solid;
    border-width: 2px;
  }
  .cloud {
    width: 30%;
  }
  .radar {
    width: 36%;
    height: 100%; /*no effect*/
  }
  .weather {
    max-width:30%  
  }
  .sky {
    background-color: rgb(36,55,114);
  }
  @media all and (max-width: 750px) {
    .cloud {
      width: 95%;
    }
    .radar {
      width: 95%;
    }
    .weather {
      max-width:95%;  
    }
    #clear_sky_clock {
      width: 95%;
    }
  }
 </style>
 <title>The sky tonight</title>
 <script src='../js/ephem.js<tags:ver/>'></script>
 <script src='../js/util.js<tags:ver/>'></script>
 <script src='../js/general-graphics.js<tags:ver/>'></script>
 <script src='code.js<tags:ver/>'></script>
 <script src='precipitation.js<tags:ver/>'></script>
 <script src='../satellite/satellite.js<tags:ver/>'></script>
 <script src='aurora.js<tags:ver/>'></script>
 <script src='occultations.js<tags:ver/>'></script>
 <script src='jupiter_satellite_phenomena.js<tags:ver/>'></script>
 <c:if test="${param.country eq 'uk'}">
 <!-- 2 scripts for openlayers js mapping tool; tile mapping, used for UK data -->
 <script src="https://cdn.polyfill.io/v2/polyfill.min.js?features=requestAnimationFrame,Element.prototype.classList"></script>
 <script src="https://openlayers.org/en/v3.20.1/build/ol.js" type="text/javascript"></script>
 </c:if> 
 <script> 
  window.onload = function() {
    var formInput = UTIL.requestParams(window);
    show(formInput, "${initParam.metOfficeApiKey}", ${initParam.isDev}, "<s:txt>en</s:txt>");
  };
  <tags:screen_change/>
 </script>
</head>
<body> 
 
 <table id='summary' class='report' style='margin-top:0;'>
    <tr><th><s:txt>Place</s:txt><th><s:txt>Time</s:txt><th title="<s:txt>Local Mean Sidereal Time</s:txt>"><s:txt>LMST</s:txt><th><s:txt>Lat</s:txt><th><s:txt>Long</s:txt><th><s:txt>UT - LT</s:txt>
 </table>
 <span style='white-space:nowrap;'  class='no-print'>
 <form id='date_time_controls' method='GET' action='' class='no-print' style="display:inline; margin:0.25em;">
   <input type='submit' name='go' value='&lt;' onclick="go_plus=0;"  class='no-print' >
   <input type='text' name='num_steps' id='num_steps' size='4'  class='no-print' >
   <select name='date_time_unit' id='date_time_unit'  class='no-print' >
     <option value='year'><s:txt>year</s:txt>
     <option value='month'><s:txt>month</s:txt>
     <option value='day' selected><s:txt>day</s:txt>
     <option value='hour'><s:txt>hour</s:txt>
     <option value='min'><s:txt>min</s:txt>
     <option value='sec'><s:txt>sec</s:txt>
   </select>
   <input type='submit' name='go' value='&gt;' onclick="go_plus=1;"  class='no-print' >
 </form>

  <form style="display:inline; margin:0.25em;">
   <select onChange="if (this.value) window.location.href=this.value">
    <tags:locations/>
   </select>
  </form>
  
  <a href='form.sky' id='link_back' class='no-print'><s:txt>Form</s:txt></a>
 </span>
 

 <P>
 <div>  
  <table id='weather' class='report weather'></table>
  <canvas id='clouds' class='no-print cloud' style='margin:0.25em;' title='<s:txt>Satellite image. Click for a larger image.</s:txt>'></canvas>
  <img id='radar' class='no-print radar'>
  <canvas id='radar_us'  class='no-print radar' style='margin:0.25em;'></canvas> 
  <div id='radar_uk'  class="no-print radar" style='margin:0.25em;float:right;'></div>  
  <!-- <div id='radar_uk'  class="no-print radar" style='float:right;margin:0.25em;'></div> -->  
 </div>
 
<P><canvas id='ecliptic' class='sky no-print' title='<s:txt>Ecliptic, showing the planets and Moon</s:txt>'></canvas>
<P><canvas id='observation_window' class='no-print' title='<s:txt>Observation window: twilight and brightness of the Moon</s:txt>'></canvas>
<P><canvas id='planisphere' class='no-print'></canvas>
 
<P style='text-align:center;'><img id='clear_sky_clock' src='' class='no-print' style='border-radius: 10px;'>
<P>
<canvas id='libration' class='sky no-print' title='<s:txt>Lunar libration, phase, and change in size from mean</s:txt>'></canvas>
<canvas id='galilean_satellites'  class='sky no-print'  style='margin-right:0.25em;' title='<s:txt>Galilean satellites of Jupiter</s:txt>'></canvas>
 
 <table id='sun_and_moon_and_planets_table' class='report' style="text-align:right;">
  <caption><s:txt>Sun, Moon, and Planets</s:txt> 
   <tr style="text-align:center;">
    <th><s:txt>Name</s:txt>
    <th><s:txt>Const</s:txt>
    <th><s:txt>Elong</s:txt>
    <th><s:txt>Mag</s:txt>
    <th><s:txt>Size</s:txt>
    <th><s:txt>Illum</s:txt>
    <th><s:txt>Alt</s:txt>
    <th><s:txt>Az</s:txt>
    <th>α
    <th>δ
    <th><s:txt>Rise</s:txt>
    <th><s:txt>Transit</s:txt>
    <th><s:txt>Set</s:txt>  
 </table>
 
 <table id='diary_table' class='report'>
   <caption><s:txt>Sky Diary</s:txt>
   <tr>
    <th><s:txt>When</s:txt> 
    <th><s:txt>Description</s:txt>
 </table>
 
 <table id='minor_planet_table' class='report'>
   <caption><s:txt>Minor Planets</s:txt> 
   <tr>
    <th><s:txt>Name</s:txt>
    <th><s:txt>Elong</s:txt>
    <th><s:txt>Mag</s:txt>
    <th><s:txt>Alt</s:txt>
    <th><s:txt>Az</s:txt>
    <th>α
    <th>δ
 </table>
 
 <table id='comet_table' class='report'>
  <caption><a href="https://www.ast.cam.ac.uk/~jds/"><s:txt>Bright Comets</s:txt></a> 
  <tr>
    <th><s:txt>Name</s:txt>
    <th><s:txt>Elong</s:txt>
    <th><s:txt>Mag</s:txt>
    <th><s:txt>Alt</s:txt>
    <th><s:txt>Az</s:txt>
    <th>α
    <th>δ
    <th><s:txt>Trend</s:txt>
    <th><s:txt>Visible</s:txt>
 </table>
 
 <table id='jupiter_satellite_phenomena' class='report' title='<s:txt>Times are approximate</s:txt>'>
  <caption><s:txt>Jupiter Satellite Phenomena</s:txt>
   <tr>
    <th><s:txt>When</s:txt>
    <th><s:txt>Satellite</s:txt>
    <th><s:txt>Event</s:txt>
 </table>
 
 <table id='aurora_table' class='report'>
   <caption><a href="http://www.swpc.noaa.gov/products/aurora-30-minute-forecast"><s:txt>Auroral Activity Level</s:txt></a> 
 </table>
 
 <table id='occultation_table' class='report' title='<s:txt>Occultation times are approximate</s:txt>'>
  <caption><a href="http://www.lunar-occultations.com/iota/iotandx.htm"><s:txt>Occultations</s:txt></a>
  <tr>
   <th><s:txt>When</s:txt>
   <th><s:txt>ZC</s:txt>
   <th><s:txt>Mag</s:txt>
   <th><s:txt>Ph</s:txt>
   <th><s:txt>El</s:txt>
   <th><s:txt>PA</s:txt>
   <th><s:txt>Alt</s:txt>
 </table>
 
 <table id='meteor_shower_table' class='report'>
   <caption><s:txt>Meteor Showers</s:txt> 
   <tr>
    <th><s:txt>Name</s:txt> 
    <th><s:txt>Desig</s:txt>
    <th><s:txt>Peak</s:txt>
    <th><s:txt>Km/s</s:txt>
    <th><s:txt>Alt</s:txt>
    <th><s:txt>Az</s:txt>
    <th>r
    <th><s:txt>ZHR</s:txt>
    <th><s:txt>ZHR</s:txt> x
 </table>

 <br><table id='messiers' class='report'>
  <caption><s:txt>Messier Objects</s:txt>
  <tr>
    <th><s:txt>Name</s:txt>
    <th><s:txt>Const</s:txt>
    <th><s:txt>Chart</s:txt>
    <th><s:txt>Mag</s:txt>
    <th><s:txt>Alt</s:txt>
    <th><s:txt>Az</s:txt>
    <th>α
    <th>δ
    <th><s:txt>Type</s:txt>
    <th><s:txt>Comment</s:txt>
 </table>
 
 <br><table id='caldwells' class='report'>
  <caption><s:txt>Caldwell Objects</s:txt>
  <tr>
    <th><s:txt>Name</s:txt>
    <th><s:txt>Const</s:txt>
    <th><s:txt>Chart</s:txt>
    <th><s:txt>Mag</s:txt>
    <th><s:txt>Alt</s:txt>
    <th><s:txt>Az</s:txt>
    <th>α
    <th>δ
    <th><s:txt>Type</s:txt>
    <th><s:txt>Comment</s:txt>
 </table>
 
 <tags:analytics/>
 
</body>
</html>