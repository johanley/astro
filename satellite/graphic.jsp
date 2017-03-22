<!doctype html>
<html lang='<s:txt>en</s:txt>'>
<head>
 <meta charset='UTF-8'>
 <meta name="keywords" content="<s:txt>canada, united states, weather, clouds, satellite, goes, visible, infrared, imagery</s:txt>">
 <meta name="description" content="<s:txt>Recent satellite images from the GOES weather satellite (NOAA).</s:txt>">
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
 </style>
 <title><s:txt>Satellite images</s:txt></title>
 <script src='../js/ephem.js<tags:ver/>'></script>
 <script src='../js/util.js<tags:ver/>'></script>
 <script src='../js/general-graphics.js<tags:ver/>'></script>
 <script src='code.js<tags:ver/>'></script>
 <script src='satellite.js<tags:ver/>'></script>
 <script>
  window.onload = function() {
    var formInput = UTIL.requestParams(window);
    show_large_satellite_image(formInput, 'clouds');
  };
 </script>
</head>
<body> 

 <!-- these links are set dynamically in code --> 
 <table class='report'>
   <td><a id='incr_w' href='' title='<s:txt>West</s:txt>'>◀</a>
   <td><a id='incr_e' href='' title='<s:txt>East</s:txt>'>▶</a>
   <td><a id='incr_n' href='' title='<s:txt>North</s:txt>'>▲</a>
   <td><a id='incr_s' href='' title='<s:txt>South</s:txt>'>▼</a>
   <td><a id='incr_in' href='' title='<s:txt>Zoom in</s:txt>'>+</a>
   <td><a id='incr_out' href='' title='<s:txt>Zoom out</s:txt>'>-</a>
 </table>
 
  <form  style="display:inline; margin:0.25em;">
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
    
 <a href='form.sky'><s:txt>Form</s:txt></a>

 <br><canvas id='clouds' class='no-print' style='border-radius:10px; margin:0.25em;'></canvas>
 
 <P><s:txt>Navigation: arrow keys, page up/down, mouse-drag, mouse-wheel.</s:txt>
 
 <tags:analytics/>
  
</body>
</html>