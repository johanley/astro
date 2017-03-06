<!doctype html>
<html lang='en'>
<head>
 <meta charset='UTF-8'>
 <meta name="keywords" content="astronomy">
 <meta name="description" content="Star chart for amateur astronomers">
 <meta name="viewport" content="width=device-width">
 <link rel="stylesheet" type="text/css" href="../css/styles.css<tags:ver/>" media="all">
 <style>
  body {
    margin:0.5em;
    padding: 0;
    font: 0.75em Verdana, Arial, Helvetica, sans-serif;
  }
  table.report {
    display: inline-block;
    vertical-align: top; 
  }
  #canvas-container {
     width: 100%;
     text-align:center;
     margin-top: 1.0em;
  }
  canvas {
     display: inline;
  }
  canvas#star_chart {
    border-radius:10px; 
    margin:0.25em; 
  }
}    
 </style>
 <title>Star chart</title>
 <script src='../js/ephem.js<tags:ver/>'></script>
 <script src='../js/util.js<tags:ver/>'></script>
 <script src='../js/general-graphics.js<tags:ver/>'></script>
 <script src='code.js<tags:ver/>'></script>
 <script> 
  window.onload = function() {
    var formInput = UTIL.requestParamsMultivalued(window, ['show_identifiers']);
    show(formInput);
  };
 </script>
</head>
<body> 
 
<div style='width:100%; text-align:center;'>

 <form id='zoom_controls' method='GET' action='' style="display:inline; margin:0.25em;" title='Or use the page up/down keys, or mouse scroll'>
   <input type='submit' id='zoom_plus' name='go' value='+' onclick="zoom_go_op=2;">
   <input type='submit' id='zoom_minus' name='go' value='-' onclick="zoom_go_op=1;">
 </form>

 <form id='direction_controls' method='GET' action='' style="display:inline; margin:0.25em;" title='Or use the arrow keys'>
   <input type='submit' id='incr_w' name='go' value='◀' onclick="go_op=1;">
   <input type='submit' id='incr_e' name='go' value='▶' onclick="go_op=2;">
   <input type='submit' id='incr_n' name='go' value='▲' onclick="go_op=3;">
   <input type='submit' id='incr_s' name='go' value='▼' onclick="go_op=4;">
   <input type='text' name='num_degs' id='num_degs' size='2' value=''> degs
 </form>
 
 <table id='summary' class='report' style='margin-top:0;'>
   <tr><th>Place<th>Time<th>Lat<th>Long<th title="Num minutes UT is in advance of browser's local time (LT)">UT - LT
 </table>
 
 <form id='date_time_controls' method='GET' action='' style="display:inline; margin:0.25em;">
   <input type='submit' name='go' value='&lt;' onclick="go_plus=0;">
   <input type='text' name='num_steps' id='num_steps' size='4' value='1'>
   <select name='date_time_unit' id='date_time_unit'>
     <option>year
     <option>month
     <option selected>day
     <option>hour
     <option>min
     <option>sec
   </select>
   <input type='submit' name='go' value='&gt;' onclick="go_plus=1;">
 </form>

 
 <a id='link_back' href='form.sky'>Form</a>
</div>

<div id="canvas-container">  
 <canvas id='star_chart'></canvas>
</div>
 
 <tags:analytics/>
 
</body>
</html>