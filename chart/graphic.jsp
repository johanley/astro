<!doctype html>
<html lang='<s:txt>en</s:txt>'>
<head>
 <meta charset='UTF-8'>
 <meta name="keywords" content="<s:txt>astronomy, custom, star chart</s:txt>">
 <meta name="description" content="<s:txt>Star charts for amateur astronomers.</s:txt>">
 <meta name="viewport" content="width=device-width">
 <link rel="stylesheet" type="text/css" href="../css/styles.css<tags:ver/>" media="all">
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
  canvas#star_chart {
    border-radius:10px; 
    margin:0.25em; 
  }
}    
 </style>
 <title><s:txt>Star chart</s:txt></title>
 <script src='../js/ephem.js<tags:ver/>'></script>
 <script src='../js/util.js<tags:ver/>'></script>
 <script src='../js/general-graphics.js<tags:ver/>'></script>
 <script src='code.js<tags:ver/>'></script>
 <script>
  var show_screen = function(){
    var formInput = UTIL.requestParamsMultivalued(window, ['show_identifiers']);
    show(formInput);
  }; 
  window.onload = show_screen;
  <tags:screen_change/>
 </script>
</head>
<body> 
 
<div style='width:100%; text-align:center;'>

 <form id='zoom_controls' method='GET' action='' style="display:inline; margin:0.25em;">
   <input type='submit' id='zoom_plus' name='go' value='+' onclick="zoom_go_op=2;">
   <input type='submit' id='zoom_minus' name='go' value='-' onclick="zoom_go_op=1;">
 </form>

 <form id='direction_controls' method='GET' action='' style="display:inline; margin:0.25em;">
   <input type='submit' id='incr_w' name='go' value='◀' onclick="go_op=1;">
   <input type='submit' id='incr_e' name='go' value='▶' onclick="go_op=2;">
   <input type='submit' id='incr_n' name='go' value='▲' onclick="go_op=3;">
   <input type='submit' id='incr_s' name='go' value='▼' onclick="go_op=4;">
   <input type='text' name='num_degs' id='num_degs' size='2' value=''> <s:txt>degs</s:txt>
 </form>
 
 <form id='date_time_controls' method='GET' action='' style="display:inline; margin:0.25em;">
   <input type='submit' name='go' value='&lt;' onclick="go_plus=0;">
   <input type='text' name='num_steps' id='num_steps' size='2' value='1'>
   <select name='date_time_unit' id='date_time_unit'>
     <option value='year'><s:txt>year</s:txt>
     <option value='month'><s:txt>month</s:txt>
     <option value='day' selected><s:txt>day</s:txt>
     <option value='hour'><s:txt>hour</s:txt>
     <option value='min'><s:txt>min</s:txt>
     <option value='sec'><s:txt>sec</s:txt>
   </select>
   <input type='submit' name='go' value='&gt;' onclick="go_plus=1;">
 </form>

 
 <a id='link_back' href='form.sky'><s:txt>Form</s:txt></a>
</div>

<div id="canvas-container">  
 <canvas id='star_chart'></canvas>
</div>
 
<div style='width:100%; text-align:center;'>
 <table id='summary' class='report' style='margin-top:0;'>
   <tr>
   <th><s:txt>Place</s:txt>
   <th><s:txt>Time</s:txt>
   <th><s:txt>Lat</s:txt>
   <th><s:txt>Long</s:txt>
   <th title="<s:txt>Num minutes UT is in advance of browser's local time (LT)</s:txt>">UT - LT
 </table>
 <P><s:txt>Navigation: arrow keys, page up/down, mouse-drag, mouse-wheel.</s:txt>
</div>
 
 <tags:analytics/>
 
</body>
</html>