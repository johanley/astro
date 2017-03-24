<%@ include file="/WEB-INF/TagHeader.jspf" %>
<%@ tag pageEncoding="UTF-8" trimDirectiveWhitespaces="true" %>
  var geolocation_success = function(pos){
     document.getElementById('latitude').value = pos.coords.latitude;
     document.getElementById('longitude').value = pos.coords.longitude;
  };
  var geolocation_failure = function(){
    document.getElementById('lat_long_autofill').innerHTML = '<s:txt>Disabled by browser</s:txt>';
  };
  var activate_lat_long_autofill = function(){
     var lat_long_auto = document.getElementById('lat_long_autofill');
     lat_long_auto.onclick = function(){
       if (navigator.geolocation){
         navigator.geolocation.getCurrentPosition(geolocation_success, geolocation_failure);
       }
       else {
         geolocation_failure();
       }
     };
  };
