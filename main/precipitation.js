/* 
 Show the most recent radar images from Environment Canada.
 These images show precipitation, both rain and snow.
 To do: add the location indicators to the image.
 Uses server-side code! (CORS problems.) 
*/
var showPrecipitation = function(radar_station, is_dev){

  var image = document.getElementById('radar');
  var imageUrls = [];
  var frameCount = 0;
  var MILLISECONDS_BETWEEN_IMAGES = 1300;
  var NUM_IMAGES = 3;

  /* The 'www' may or may not be present; need to preserve the same style. */
  var url_preamble = function(){
    var result = document.location.href.startsWith('http://www.') ? 'http://www.' : 'http://';
    return result;
  };  
  
  var showMostRecentImages = function(radar_station) {
    //the server returns a json object with url information
    var jsonUrl = '';
    if (is_dev){ 
      jsonUrl = 'http://localhost:8081/astro/fetch/?radarStation=' + radar_station;
    }
    else {
      jsonUrl = UTIL.url_preamble() + 'astronomytonight.net/fetch/?radarStation=' + radar_station;
    }
    
    var xhr = new XMLHttpRequest();
    xhr.onreadystatechange = function() {
      if (xhr.readyState == 4 && xhr.status == 200) {
         console.log('Success: radar station json, ' + jsonUrl);
         var jsonObj = JSON.parse(xhr.responseText);
         imageUrls = jsonObj.radarUrls;
         animateImages();
      }
    };    
    xhr.ontimeout = function (e) {
        console.log('Timeout: radar station json, ' + jsonUrl);
    };
    xhr.open("GET", jsonUrl, true);
    xhr.timeout = 20 * 1000; 
    xhr.send();
  };
  
  var showMostRecentImageOnly = function(){
    var image = document.getElementById('radar');
    image.src = imageUrls[0];
  };
  
  var animateImages = function(){
    setInterval(displayNextImage, MILLISECONDS_BETWEEN_IMAGES);
  };
  
  var displayNextImage = function() {
    image.src = imageUrls[frameCount % NUM_IMAGES];
    ++frameCount;
  };
    
  showMostRecentImages(radar_station);
}; 
