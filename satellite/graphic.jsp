<!doctype html>
<html lang='en'>
<head>
 <meta charset='UTF-8'>
 <meta name="keywords" content="astronomy">
 <meta name="description" content="View NOAA images of cloud cover">
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
 <title>GOES satellite images</title>
 <script src='../js/ephem.js<tags:ver/>'></script>
 <script src='../js/util.js<tags:ver/>'></script>
 <script src='code.js<tags:ver/>'></script>
 <script>
  var click_on = function(id){
    var link = document.getElementById(id);
    link.click();
  };
  var use_arrow_keys = function(){
      var PAGE_UP = 33, PAGE_DOWN = 34;
      var LEFT_ARROW = 37, UP_ARROW = 38, RIGHT_ARROW = 39, DOWN_ARROW = 40;
      document.onkeydown = function(e) {
        switch (e.keyCode) {
          case PAGE_UP: 
              click_on('incr_in');
              break;
          case PAGE_DOWN: 
              click_on('incr_out');
              break;
          case LEFT_ARROW:
              click_on('incr_w');
              break;
          case UP_ARROW:
              click_on('incr_n');
              break;
          case RIGHT_ARROW:
              click_on('incr_e');
              break;
          case DOWN_ARROW:
              click_on('incr_s');
              break;
        }
    };
  };
  window.onload = function() {
    var formInput = UTIL.requestParams(window);
    show_large_scale_clouds(formInput, 'clouds');
    use_arrow_keys();
  };
 </script>
</head>
<body> 
  
 <table class='report' >
   <tr>
     <td><a href='form.sky'>Back to Form</a>
     <td><a href='graphic.sky?latitude=52&longitude=-121&degrees_on_a_side=15&pixels_on_a_side=800&layer=auto_detect'>BC, AB</a>
     <td><a href='graphic.sky?latitude=53&longitude=-102&degrees_on_a_side=15&pixels_on_a_side=800&layer=auto_detect'>SK, MB</a>
     <td><a href='graphic.sky?latitude=48&longitude=-81.8&degrees_on_a_side=15&pixels_on_a_side=800&layer=auto_detect'>Ontario</a>
     <td><a href='graphic.sky?latitude=49.3&longitude=-72.1&degrees_on_a_side=15&pixels_on_a_side=800&layer=auto_detect'>Quebec</a>
     <td><a href='graphic.sky?latitude=48.4&longitude=-61.8&degrees_on_a_side=15&pixels_on_a_side=800&layer=auto_detect'>Atlantic Provinces</a>
     <td><a href='graphic.sky?latitude=52&longitude=-110&degrees_on_a_side=25&pixels_on_a_side=800&layer=auto_detect'>Canada West</a>
     <td><a href='graphic.sky?latitude=49.8&longitude=-71.7&degrees_on_a_side=25&pixels_on_a_side=800&layer=auto_detect'>Canada East</a>
     <td><a href='graphic.sky?latitude=42&longitude=-75&degrees_on_a_side=15&pixels_on_a_side=800&layer=auto_detect'>US NE</a>
     <td><a href='graphic.sky?latitude=32&longitude=-83&degrees_on_a_side=15&pixels_on_a_side=800&layer=auto_detect'>US SE</a>
     <td><a href='graphic.sky?latitude=37.6&longitude=-118.5&degrees_on_a_side=15&pixels_on_a_side=800&layer=auto_detect'>US California</a>
   </tr>
 </table>

 <!-- these links are set dynamically in code --> 
 <table class='report' title='Or use the arrow keys, page up-down keys'>
   <td><a id='incr_w' href='' title='West'>◀</a>
   <td><a id='incr_e' href='' title='East'>▶</a>
   <td><a id='incr_n' href='' title='North'>▲</a>
   <td><a id='incr_s' href='' title='South'>▼</a>
   <td><a id='incr_in' href='' title='Zoom in'>+</a>
   <td><a id='incr_out' href='' title='Zoom out'>-</a>
 </table>
 
 <!-- <canvas id='clouds' class='no-print' style='border-radius:10px; margin:0.25em; width:95%;'></canvas> -->
 <canvas id='clouds' class='no-print' style='border-radius:10px; margin:0.25em;'></canvas>
 
 <P>Navigation: arrow keys, page up/down, mouse-drag, mouse-wheel.
 <tags:analytics/>
  
</body>
</html>