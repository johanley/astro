<!doctype html>
<html lang='<s:txt>en</s:txt>'>
<head>
 <meta charset='UTF-8'>
 <meta name="keywords" content="<s:txt>canada, united states, weather, clouds, satellite, goes, visible, infrared, imagery</s:txt>">
 <meta name="description" content="<s:txt>Recent satellite images from the GOES weather satellite (NOAA).</s:txt>">
 <meta name="viewport" content="width=device-width">
 <link rel="stylesheet" href="../css/styles.css<tags:ver/>" media="all">
 <style>
  body {
    margin:0.5em;
    padding: 0;
  }
  table.report {
    display: inline-block;
    vertical-align: top; 
  }
  #canvas-container {
     width: 100%;
     text-align:center;
     margin-top: 0.5em;
     margin-bottom: 0.5em;
  }
  canvas {
     display: inline;
  }
  canvas#clouds {
    border-radius:10px; 
    margin:0.25em; 
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
  <tags:screen_change/>
 </script>
</head>
<body> 

<div style='width:100%; text-align:center;'>

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
    <tags:locations_sat_image/>
   </select>
  </form>
  
 <a href='form.sky'><s:txt>Form</s:txt></a>
</div>

<div id="canvas-container">  
 <canvas id='clouds'></canvas>
</div>
 
<div style='width:100%; text-align:center;'>
  <s:txt>Navigation: arrow keys, page up/down, mouse-drag, mouse-wheel.</s:txt>
</div>

<tags:analytics/>
  
</body>
</html>