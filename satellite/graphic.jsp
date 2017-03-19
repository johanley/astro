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
  
 <table class='report' >
   <tr>
     <td><a href='form.sky'><s:txt>Back to Form</s:txt></a>
     <td><a href='graphic.sky?latitude=52&longitude=-121&degrees_on_a_side=15&pixels_on_a_side=800&layer=auto_detect'>BC, AB</a>
     <td><a href='graphic.sky?latitude=53&longitude=-102&degrees_on_a_side=15&pixels_on_a_side=800&layer=auto_detect'>SK, MB</a>
     <td><a href='graphic.sky?latitude=48&longitude=-81.8&degrees_on_a_side=15&pixels_on_a_side=800&layer=auto_detect'>Ontario</a>
     <td><a href='graphic.sky?latitude=49.3&longitude=-72.1&degrees_on_a_side=15&pixels_on_a_side=800&layer=auto_detect'>Québec</a>
     <td><a href='graphic.sky?latitude=48.4&longitude=-61.8&degrees_on_a_side=15&pixels_on_a_side=800&layer=auto_detect'><s:txt>Atlantic Provinces</s:txt></a>
     <td><a href='graphic.sky?latitude=52&longitude=-110&degrees_on_a_side=25&pixels_on_a_side=800&layer=auto_detect'><s:txt>Canada West</s:txt></a>
     <td><a href='graphic.sky?latitude=49.8&longitude=-71.7&degrees_on_a_side=25&pixels_on_a_side=800&layer=auto_detect'><s:txt>Canada East</s:txt></a>
     <td><a href='graphic.sky?latitude=42&longitude=-75&degrees_on_a_side=15&pixels_on_a_side=800&layer=auto_detect'><s:txt>US NE</s:txt></a>
     <td><a href='graphic.sky?latitude=32&longitude=-83&degrees_on_a_side=15&pixels_on_a_side=800&layer=auto_detect'><s:txt>US SE</s:txt></a>
     <td><a href='graphic.sky?latitude=37.6&longitude=-118.5&degrees_on_a_side=15&pixels_on_a_side=800&layer=auto_detect'><s:txt>US California</s:txt></a>
   </tr>
 </table>

 <!-- these links are set dynamically in code --> 
 <table class='report'>
   <td><a id='incr_w' href='' title='<s:txt>West</s:txt>'>◀</a>
   <td><a id='incr_e' href='' title='<s:txt>East</s:txt>'>▶</a>
   <td><a id='incr_n' href='' title='<s:txt>North</s:txt>'>▲</a>
   <td><a id='incr_s' href='' title='<s:txt>South</s:txt>'>▼</a>
   <td><a id='incr_in' href='' title='<s:txt>Zoom in</s:txt>'>+</a>
   <td><a id='incr_out' href='' title='<s:txt>Zoom out</s:txt>'>-</a>
 </table>
 
 <canvas id='clouds' class='no-print' style='border-radius:10px; margin:0.25em;'></canvas>
 
 <P><s:txt>Navigation: arrow keys, page up/down, mouse-drag, mouse-wheel.</s:txt>
 
 <tags:analytics/>
  
</body>
</html>