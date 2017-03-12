<!doctype html>
<html lang='en'>
<head>
 <meta charset='UTF-8'>
 <meta name="keywords" content="astronomy, ephemerides, low precision, amateur, canada, united states, united kingdom, weather, radar, clouds">
 <meta name="description" content="Data of interest to amateur astronomers. Includes weather for North America and the United Kingdom.">
 <meta name="viewport" content="width=device-width">
 <link rel="stylesheet" href="../css/styles.css<tags:ver/>" media="all">
 <title>The sky tonight</title>
</head>
 <script src='find_closest_clear_sky_clock.js<tags:ver/>'></script>
 <script src='../js/util.js<tags:ver/>'></script>
 <script src='../js/ephem.js<tags:ver/>'></script>
 <script> 
  window.onload = function() {
   var is_dev = ${initParam.isDev}; 
   var exclude_none = document.getElementById('exclude_none');
   var exclude_all = document.getElementById('exclude_all');
   var input_form = document.getElementById('input_form');
   var control;
   var toggle = function(on_off){
     for(var i=0; i < input_form.elements.length; ++i){
       control = input_form.elements[i];
       if (control.name.startsWith('exclude_')){
         control.checked = on_off;
       }
     }
   };
   exclude_none.onclick = function(){
     toggle(false);
   };
   exclude_all.onclick = function(){
     toggle(true);
   };
   
   var csc_find_nearest = document.getElementById('csc_find_nearest_station');
   csc_find_nearest.onclick = function(){
     var lat = document.getElementById('latitude').value;
     var long = document.getElementById('longitude').value;
     this.innerHTML = 'Scanning.....';
     if (lat && long){
       var target = {};
       target.φ = lat; 
       target.λ = long;
       var callback = function(id){
         //console.log('The nearest station id: ' + id);
         document.getElementById('clear_sky_clock_station_id').value = id;
         document.getElementById('csc_find_nearest_station').innerHTML = 'Find nearest';
       };
       find_closest_csc_url(target, callback);
     }
   };
   
   var lat_long_auto = document.getElementById('lat_long_autofill');
   lat_long_auto.onclick = function(){
     if (navigator.geolocation){
       navigator.geolocation.getCurrentPosition(
         function(pos){
           document.getElementById('latitude').value = pos.coords.latitude;
           document.getElementById('longitude').value = pos.coords.longitude;
         }
       );
     }
     else {
       lat_long_auto.innerHTML = 'Disabled by browser';
     }
   };

   /*
    IMPORTANT: Naming convention for rows and controls, in the form below. 
    Example
      x=weather_station_id_canada, x'=weather_station_id  
    Allows to hide-disable/show-enable, but share the same 'name', such that the code that receives the form input 
    references the same entity with the same name. 
     row (label+control): id=row_x 
     control: id=x, name=x', where x' denotes the form in which the country-suffix is dropped 
   */ 
   var control_id = function(row_id){
     var result = row_id.substring(4); //chop off 'row_'
     return result;
   };
   var show_hide_rows_per_country = function(){
     //every non-static item has to be handled for every country
     var country = document.getElementById('country').value;
     var show = function(row_id){
       var row = document.getElementById(row_id);
       row.style.display = '';
       var control = document.getElementById(control_id(row_id));
       control.disabled = false; 
     };
     var hide = function(row_id){
       var row = document.getElementById(row_id);
       row.style.display = 'none'; 
       var control = document.getElementById(control_id(row_id));
       control.disabled = true; 
     };
     var show_all = function(ids){
       for (var i = 0; i < ids.length; ++i){
         show(ids[i]);
       }
     };
     var hide_all = function(ids){
       for (var i = 0; i < ids.length; ++i){
         hide(ids[i]);
       }
     };
     var canada_only = ['row_radar_station_cda', 'row_weather_station_cda', 'row_prov_cda'];
     var us_only = ['row_radar_station_us'];
     var uk_only = ['row_weather_station_uk'];
     if ('canada' === country){
       show_all(canada_only);
       hide_all(us_only); 
       hide_all(uk_only); 
       show('row_clear_sky_clock_station_id');
     }
     else if ('us' === country){
       show_all(us_only);
       hide_all(canada_only); 
       hide_all(uk_only); 
       show('row_clear_sky_clock_station_id');
     }
     else if ('uk' === country){
       show_all(uk_only);
       hide_all(canada_only); 
       hide_all(us_only); 
       hide('row_clear_sky_clock_station_id');
     }
     else if ('other' === country){
       hide_all(canada_only); 
       hide_all(us_only); 
       hide('row_clear_sky_clock_station_id');
     } 
   };
   //react to user selections of a new country
   document.getElementById('country').onchange = show_hide_rows_per_country;

   //initial pre-pop, the starting-state of the form
   //inject the URL params into the form, once only; after this, user selections can be out of sync with the URL
   UTIL.prepopulate_form_controls([]);
   show_hide_rows_per_country();
   
  };
 </script>
<body>
<h2>Astronomy Tonight</h2>

This tool collects together, on a single page, basic information of interest to most amateur astronomers. 
It's most useful for residents of Canada, the US, and the UK. Outside of those areas, no weather data is shown.

<P>Below you will find a large form that lets you customize the output data.
Its defaults are for Ottawa, the capital of Canada. 
As a convenience, here are the results for some selected locations (Canada, US, and UK only):

<P>
  <form >
   <select onChange="if (this.value) window.location.href=this.value">
     <option value=''>Choose a location
     <option value=''>---Canada---
     <option value='graphic%2Esky?country=cda&prov=AB&location_name=Banff&latitude=51%2E18&longitude=-115%2E57&locations=51%2E18,-115%2E57&weather_station=49&radar_station=XSM&clear_sky_clock_station_id=BnffAB'>AB Banff
     <option value='graphic%2Esky?country=cda&prov=AB&location_name=Calgary&latitude=51%2E03&longitude=-114%2E09&locations=51%2E03,-114%2E09&weather_station=52&radar_station=XSM&clear_sky_clock_station_id=Calgary'>AB Calgary
     <option value='graphic%2Esky?country=cda&prov=AB&location_name=Edmonton&latitude=53%2E54&longitude=-113%2E50&locations=53%2E54,-113%2E50&weather_station=50&radar_station=WHK&clear_sky_clock_station_id=Edmonton'>AB Edmonton
     <option value='graphic%2Esky?country=cda&prov=AB&location_name=Elk+Island+–+Beaver+Hills&latitude=53%2E58&longitude=-112%2E82&locations=53%2E58,-112%2E82&weather_station=63&radar_station=WHK&clear_sky_clock_station_id=Blackfoot'>AB Elk Island – Beaver Hills
     <option value='graphic%2Esky?country=cda&prov=AB&location_name=Fort+McMurray&latitude=56%2E73&longitude=-111%2E38&locations=56%2E73,-111%2E38&weather_station=20&radar_station=&clear_sky_clock_station_id=FtMcmrryAB'>AB Fort McMurray
     <option value='graphic%2Esky?country=cda&prov=AB&location_name=Grande+Prairie&latitude=55%2E17&longitude=-118%2E80&locations=55%2E17,-118%2E80&weather_station=31&radar_station=WWW&clear_sky_clock_station_id=NDgObAB'>AB Grande Prairie
     <option value='graphic%2Esky?country=cda&prov=AB&location_name=Jasper&latitude=52%2E87&longitude=-118%2E08&locations=52%2E87,-118%2E08&weather_station=70&radar_station=&clear_sky_clock_station_id=MlngLkAB'>AB Jasper
     <option value='graphic%2Esky?country=cda&prov=AB&location_name=Lethbridge&latitude=49%2E69&longitude=-112%2E85&locations=49%2E69,-112%2E85&weather_station=30&radar_station=XSM&clear_sky_clock_station_id=LethAB'>AB Lethbridge
     <option value='graphic%2Esky?country=cda&prov=AB&location_name=Red+Deer&latitude=52%2E26&longitude=-113%2E80&locations=52%2E26,-113%2E80&weather_station=29&radar_station=XSM&clear_sky_clock_station_id=RdDrAB'>AB Red Deer
     <option value='graphic%2Esky?country=cda&prov=BC&location_name=Cranbrook&latitude=49%2E51&longitude=-115%2E76&locations=49%2E51,-115%2E76&weather_station=77&radar_station=&clear_sky_clock_station_id=SprcRdObBC'>BC Cranbrook
     <option value='graphic%2Esky?country=cda&prov=BC&location_name=Kamloops&latitude=50%2E67&longitude=-120%2E29&locations=50%2E67,-120%2E29&weather_station=45&radar_station=XSS&clear_sky_clock_station_id=Kamloops'>BC Kamloops
     <option value='graphic%2Esky?country=cda&prov=BC&location_name=Penticton&latitude=49%2E49&longitude=-119%2E57&locations=49%2E49,-119%2E57&weather_station=84&radar_station=XSS&clear_sky_clock_station_id=PntctonBC'>BC Penticton
     <option value='graphic%2Esky?country=cda&prov=BC&location_name=Prince+George&latitude=53%2E92&longitude=-122%2E75&locations=53%2E92,-122%2E75&weather_station=79&radar_station=XPG&clear_sky_clock_station_id=Prince_George'>BC Prince George
     <option value='graphic%2Esky?country=cda&prov=BC&location_name=Terrace&latitude=54%2E53&longitude=-128%2E61&locations=54%2E53,-128%2E61&weather_station=80&radar_station=&clear_sky_clock_station_id=TrrcBC'>BC Terrace
     <option value='graphic%2Esky?country=cda&prov=BC&location_name=Tofino&latitude=49%2E14&longitude=-125%2E90&locations=49%2E14,-125%2E90&weather_station=17&radar_station=XSI&clear_sky_clock_station_id=TofinoBC'>BC Tofino
     <option value='graphic%2Esky?country=cda&prov=BC&location_name=Vancouver&latitude=49%2E30&longitude=-123%2E14&locations=49%2E30,-123%2E14&weather_station=74&radar_station=WUJ&clear_sky_clock_station_id=Vancouver'>BC Vancouver
     <option value='graphic%2Esky?country=cda&prov=BC&location_name=Victoria&latitude=48%2E41&longitude=-123%2E37&locations=48%2E41,-123%2E37&weather_station=85&radar_station=XSI&clear_sky_clock_station_id=Victoria'>BC Victoria
     <option value='graphic%2Esky?country=cda&prov=BC&location_name=Whistler&latitude=50%2E12&longitude=-122%2E96&locations=50%2E12,-122%2E96&weather_station=86&radar_station=WUJ&clear_sky_clock_station_id=WhstlrBC'>BC Whistler
     <option value='graphic%2Esky?country=cda&prov=MB&location_name=Brandon&latitude=49%2E85&longitude=-99%2E95&locations=49%2E85,-99%2E95&weather_station=52&radar_station=XFW&clear_sky_clock_station_id=BrndnMB'>MB Brandon
     <option value='graphic%2Esky?country=cda&prov=MB&location_name=Churchill&latitude=58%2E77&longitude=-94%2E16&locations=58%2E77,-94%2E16&weather_station=42&radar_station=&clear_sky_clock_station_id=ChrchllMB'>MB Churchill
     <option value='graphic%2Esky?country=cda&prov=MB&location_name=Dauphin&latitude=51%2E15&longitude=-100%2E06&locations=51%2E15,-100%2E06&weather_station=58&radar_station=XFW&clear_sky_clock_station_id=MRossObs'>MB Dauphin
     <option value='graphic%2Esky?country=cda&prov=MB&location_name=The+Pas&latitude=53%2E82&longitude=-101%2E25&locations=53%2E82,-101%2E25&weather_station=30&radar_station=&clear_sky_clock_station_id=TPasMB'>MB The Pas
     <option value='graphic%2Esky?country=cda&prov=MB&location_name=Thompson&latitude=55%2E75&longitude=-97%2E85&locations=55%2E75,-97%2E85&weather_station=34&radar_station=&clear_sky_clock_station_id=ThmpsnMB'>MB Thompson
     <option value='graphic%2Esky?country=cda&prov=MB&location_name=Winnipeg&latitude=49%2E91&longitude=-97%2E14&locations=49%2E91,-97%2E14&weather_station=38&radar_station=XWL&clear_sky_clock_station_id=Winnipeg'>MB Winnipeg
     <option value='graphic%2Esky?country=cda&prov=NB&location_name=Bathurst&latitude=47%2E61&longitude=-65%2E64&locations=47%2E61,-65%2E64&weather_station=28&radar_station=XNC&clear_sky_clock_station_id=BathurstNB'>NB Bathurst
     <option value='graphic%2Esky?country=cda&prov=NB&location_name=Edmunston&latitude=47%2E37&longitude=-68%2E32&locations=47%2E37,-68%2E32&weather_station=32&radar_station=XAM&clear_sky_clock_station_id=RachelObNB'>NB Edmunston
     <option value='graphic%2Esky?country=cda&prov=NB&location_name=Fredericton&latitude=45%2E96&longitude=-66%2E65&locations=45%2E96,-66%2E65&weather_station=29&radar_station=XNC&clear_sky_clock_station_id=FrederictonNB'>NB Fredericton
     <option value='graphic%2Esky?country=cda&prov=NB&location_name=Fundy+National+Park&latitude=45%2E55&longitude=-65%2E02&locations=45%2E55,-65%2E02&weather_station=5&radar_station=XNC&clear_sky_clock_station_id=FndyNPNB'>NB Fundy National Park
     <option value='graphic%2Esky?country=cda&prov=NB&location_name=Kouchibouguac&latitude=46%2E83&longitude=-64%2E93&locations=46%2E83,-64%2E93&weather_station=9&radar_station=XNC&clear_sky_clock_station_id=KchbgNPNB'>NB Kouchibouguac
     <option value='graphic%2Esky?country=cda&prov=NB&location_name=Moncton&latitude=46%2E08&longitude=-64%2E78&locations=46%2E08,-64%2E78&weather_station=36&radar_station=XNC&clear_sky_clock_station_id=Moncton'>NB Moncton
     <option value='graphic%2Esky?country=cda&prov=NB&location_name=Mount+Carleton+Provincial+Park&latitude=47%2E43&longitude=-66%2E91&locations=47%2E43,-66%2E91&weather_station=10&radar_station=XAM&clear_sky_clock_station_id=AASPNB'>NB Mount Carleton Provincial Park
     <option value='graphic%2Esky?country=cda&prov=NB&location_name=St%2E+John&latitude=45%2E27&longitude=-66%2E07&locations=45%2E27,-66%2E07&weather_station=23&radar_station=XNC&clear_sky_clock_station_id=SntJhnNMB'>NB St. John
     <option value='graphic%2Esky?country=cda&prov=NL&location_name=Badger&latitude=49%2E50&longitude=-56%2E07&locations=49%2E50,-56%2E07&weather_station=34&radar_station=XME&clear_sky_clock_station_id=SprngdNL'>NL Badger
     <option value='graphic%2Esky?country=cda&prov=NL&location_name=Cornerbrook&latitude=48%2E95&longitude=-57%2E95&locations=48%2E95,-57%2E95&weather_station=41&radar_station=XME&clear_sky_clock_station_id=CrnrBrkNF'>NL Cornerbrook
     <option value='graphic%2Esky?country=cda&prov=NL&location_name=Happy+Valley-Goose+Bay&latitude=53%2E30&longitude=-60%2E35&locations=53%2E30,-60%2E35&weather_station=23&radar_station=&clear_sky_clock_station_id=GssByNFLD'>NL Happy Valley-Goose Bay
     <option value='graphic%2Esky?country=cda&prov=NL&location_name=St%2E+John%27s&latitude=47%2E56&longitude=-52%2E69&locations=47%2E56,-52%2E69&weather_station=24&radar_station=WTP&clear_sky_clock_station_id=StJohns'>NL St. John's
     <option value='graphic%2Esky?country=cda&prov=NS&location_name=Halifax&latitude=44%2E66&longitude=-63%2E59&locations=44%2E66,-63%2E59&weather_station=19&radar_station=XGO&clear_sky_clock_station_id=Halifax'>NS Halifax
     <option value='graphic%2Esky?country=cda&prov=NS&location_name=Kejimkujik&latitude=44%2E42&longitude=-65%2E26&locations=44%2E42,-65%2E26&weather_station=42&radar_station=XGO&clear_sky_clock_station_id=KjmkjkNS'>NS Kejimkujik
     <option value='graphic%2Esky?country=cda&prov=NS&location_name=Sydney&latitude=46%2E13&longitude=-60%2E19&locations=46%2E13,-60%2E19&weather_station=31&radar_station=XMB&clear_sky_clock_station_id=SydnyNS'>NS Sydney
     <option value='graphic%2Esky?country=cda&prov=NS&location_name=Truro&latitude=45%2E36&longitude=-63%2E29&locations=45%2E36,-63%2E29&weather_station=25&radar_station=XGO&clear_sky_clock_station_id=TruroNS'>NS Truro
     <option value='graphic%2Esky?country=cda&prov=NS&location_name=Yarmouth&latitude=43%2E83&longitude=-66%2E12&locations=43%2E83,-66%2E12&weather_station=29&radar_station=XGO&clear_sky_clock_station_id=ArgyleNS'>NS Yarmouth
     <option value='graphic%2Esky?country=cda&prov=NT&location_name=Yellowknife&latitude=62%2E45&longitude=-114%2E37&locations=62%2E45,-114%2E37&weather_station=24&radar_station=&clear_sky_clock_station_id=YllwknfNWT'>NT Yellowknife
     <option value='graphic%2Esky?country=cda&prov=NU&location_name=Iqaluit&latitude=63%2E75&longitude=-68%2E52&locations=63%2E75,-68%2E52&weather_station=21&radar_station=&clear_sky_clock_station_id=IqltNvT'>NU Iqaluit
     <option value='graphic%2Esky?country=cda&prov=ON&location_name=Algonquin+Park&latitude=45%2E58&longitude=-78%2E41&locations=45%2E58,-78%2E41&weather_station=29&radar_station=WBI&clear_sky_clock_station_id=RckLkAqON'>ON Algonquin Park
     <option value='graphic%2Esky?country=cda&prov=ON&location_name=Barrie&latitude=44%2E39&longitude=-79%2E69&locations=44%2E39,-79%2E69&weather_station=151&radar_station=WKR&clear_sky_clock_station_id=Barrie'>ON Barrie
     <option value='graphic%2Esky?country=cda&prov=ON&location_name=Hamilton&latitude=43%2E25&longitude=-79%2E87&locations=43%2E25,-79%2E87&weather_station=77&radar_station=WKR&clear_sky_clock_station_id=Hamilton'>ON Hamilton
     <option value='graphic%2Esky?country=cda&prov=ON&location_name=Kenora&latitude=49%2E77&longitude=-94%2E50&locations=49%2E77,-94%2E50&weather_station=96&radar_station=XDR&clear_sky_clock_station_id=KenoraON'>ON Kenora
     <option value='graphic%2Esky?country=cda&prov=ON&location_name=Kingston&latitude=44%2E23&longitude=-76%2E54&locations=44%2E23,-76%2E54&weather_station=69&radar_station=XFT&clear_sky_clock_station_id=Kingston'>ON Kingston
     <option value='graphic%2Esky?country=cda&prov=ON&location_name=London&latitude=42%2E99&longitude=-81%2E25&locations=42%2E99,-81%2E25&weather_station=137&radar_station=WSO&clear_sky_clock_station_id=London'>ON London
     <option value='graphic%2Esky?country=cda&prov=ON&location_name=Moosonee&latitude=51%2E28&longitude=-80%2E65&locations=51%2E28,-80%2E65&weather_station=113&radar_station=XTI&clear_sky_clock_station_id=MsneeON'>ON Moosonee
     <option value='graphic%2Esky?country=cda&prov=ON&location_name=North+Frontenac&latitude=44%2E92&longitude=-76%2E94&locations=44%2E92,-76%2E94&weather_station=106&radar_station=XFT&clear_sky_clock_station_id=PlvnAdON'>ON North Frontenac
     <option value='graphic%2Esky?country=cda&prov=ON&location_name=Ottawa&latitude=45%2E40&longitude=-75%2E66&locations=45%2E40,-75%2E66&weather_station=118&radar_station=XFT&clear_sky_clock_station_id=FLO'>ON Ottawa
     <option value='graphic%2Esky?country=cda&prov=ON&location_name=Parry+Sound&latitude=45%2E35&longitude=-80%2E04&locations=45%2E35,-80%2E04&weather_station=103&radar_station=WBI&clear_sky_clock_station_id=PrrySndON'>ON Parry Sound
     <option value='graphic%2Esky?country=cda&prov=ON&location_name=Peterborough&latitude=44%2E31&longitude=-78%2E32&locations=44%2E31,-78%2E32&weather_station=121&radar_station=WKR&clear_sky_clock_station_id=Peterborough'>ON Peterborough
     <option value='graphic%2Esky?country=cda&prov=ON&location_name=Sarnia&latitude=42%2E98&longitude=-82%2E41&locations=42%2E98,-82%2E41&weather_station=147&radar_station=WSO&clear_sky_clock_station_id=Sarnia'>ON Sarnia
     <option value='graphic%2Esky?country=cda&prov=ON&location_name=Sault+Ste%2E+Marie&latitude=46%2E53&longitude=-84%2E36&locations=46%2E53,-84%2E36&weather_station=162&radar_station=WGJ&clear_sky_clock_station_id=SaultON'>ON Sault Ste. Marie
     <option value='graphic%2Esky?country=cda&prov=ON&location_name=Sudbury&latitude=46%2E52&longitude=-80%2E96&locations=46%2E52,-80%2E96&weather_station=40&radar_station=WBI&clear_sky_clock_station_id=Sudbury'>ON Sudbury
     <option value='graphic%2Esky?country=cda&prov=ON&location_name=Thunder+Bay&latitude=48%2E38&longitude=-89%2E24&locations=48%2E38,-89%2E24&weather_station=100&radar_station=XNI&clear_sky_clock_station_id=ThunderBay'>ON Thunder Bay
     <option value='graphic%2Esky?country=cda&prov=ON&location_name=Timmins&latitude=48%2E48&longitude=-81%2E34&locations=48%2E48,-81%2E34&weather_station=127&radar_station=XTI&clear_sky_clock_station_id=Timmins'>ON Timmins
     <option value='graphic%2Esky?country=cda&prov=ON&location_name=Toronto&latitude=43%2E65&longitude=-79%2E38&locations=43%2E65,-79%2E38&weather_station=143&radar_station=WKR&clear_sky_clock_station_id=Toronto'>ON Toronto
     <option value='graphic%2Esky?country=cda&prov=ON&location_name=Windsor&latitude=42%2E31&longitude=-83%2E04&locations=42%2E31,-83%2E04&weather_station=94&radar_station=WSO&clear_sky_clock_station_id=Windsor'>ON Windsor
     <option value='graphic%2Esky?country=cda&prov=PE&location_name=Charlottetown&latitude=46%2E26&longitude=-63%2E14&locations=46%2E26,-63%2E14&weather_station=5&radar_station=XGO&clear_sky_clock_station_id=Charlottetown'>PE Charlottetown
     <option value='graphic%2Esky?country=cda&prov=PE&location_name=Summerside&latitude=46%2E39&longitude=-63%2E79&locations=46%2E39,-63%2E79&weather_station=3&radar_station=XGO&clear_sky_clock_station_id=SlmnPkPEI'>PE Summerside
     <option value='graphic%2Esky?country=cda&prov=QC&location_name=Gaspé&latitude=48%2E83&longitude=-64%2E51&locations=48%2E83,-64%2E51&weather_station=101&radar_station=XAM&clear_sky_clock_station_id=GaspePQ'>QC Gaspé
     <option value='graphic%2Esky?country=cda&prov=QC&location_name=Iles-de-la-Madeleine&latitude=47%2E40&longitude=-61%2E84&locations=47%2E40,-61%2E84&weather_station=103&radar_station=XMB&clear_sky_clock_station_id=BtMntTQC'>QC Iles-de-la-Madeleine
     <option value='graphic%2Esky?country=cda&prov=QC&location_name=La+Vérendrye&latitude=46%2E98&longitude=-76%2E48&locations=46%2E98,-76%2E48&weather_station=30&radar_station=XLA&clear_sky_clock_station_id=RsrvFnqPQ'>QC La Vérendrye
     <option value='graphic%2Esky?country=cda&prov=QC&location_name=Mont+Mégantic&latitude=45%2E46&longitude=-71%2E15&locations=45%2E46,-71%2E15&weather_station=136&radar_station=WVY&clear_sky_clock_station_id=Observatoires_du_Mont_Megantic'>QC Mont Mégantic
     <option value='graphic%2Esky?country=cda&prov=QC&location_name=Montréal&latitude=45%2E50&longitude=-73%2E58&locations=45%2E50,-73%2E58&weather_station=147&radar_station=WMN&clear_sky_clock_station_id=BvObQC'>QC Montréal
     <option value='graphic%2Esky?country=cda&prov=QC&location_name=Québec&latitude=46%2E82&longitude=-71%2E22&locations=46%2E82,-71%2E22&weather_station=133&radar_station=WVY&clear_sky_clock_station_id=Quebec'>QC Québec
     <option value='graphic%2Esky?country=cda&prov=QC&location_name=Rimouski&latitude=48%2E45&longitude=-68%2E53&locations=48%2E45,-68%2E53&weather_station=138&radar_station=XAM&clear_sky_clock_station_id=RmskiQC'>QC Rimouski
     <option value='graphic%2Esky?country=cda&prov=QC&location_name=Rivière-du-Loup&latitude=47%2E83&longitude=-69%2E54&locations=47%2E83,-69%2E54&weather_station=108&radar_station=WMB&clear_sky_clock_station_id=PtQllObPQ'>QC Rivière-du-Loup
     <option value='graphic%2Esky?country=cda&prov=QC&location_name=Saguenay&latitude=48%2E43&longitude=-71%2E07&locations=48%2E43,-71%2E07&weather_station=166&radar_station=WMB&clear_sky_clock_station_id=SrsSgnyPQ'>QC Saguenay
     <option value='graphic%2Esky?country=cda&prov=QC&location_name=Sept-Iles&latitude=50%2E22&longitude=-66%2E37&locations=50%2E22,-66%2E37&weather_station=141&radar_station=XAM&clear_sky_clock_station_id=SeptIlesPQ'>QC Sept-Iles
     <option value='graphic%2Esky?country=cda&prov=QC&location_name=Sherbrooke&latitude=45%2E40&longitude=-71%2E88&locations=45%2E40,-71%2E88&weather_station=136&radar_station=WVY&clear_sky_clock_station_id=BspUObQC'>QC Sherbrooke
     <option value='graphic%2Esky?country=cda&prov=QC&location_name=Trois-Rivières&latitude=46%2E35&longitude=-72%2E59&locations=46%2E35,-72%2E59&weather_station=130&radar_station=WVY&clear_sky_clock_station_id=TroisRivPQ'>QC Trois-Rivières
     <option value='graphic%2Esky?country=cda&prov=QC&location_name=Val-d%27or&latitude=48%2E10&longitude=-77%2E80&locations=48%2E10,-77%2E80&weather_station=149&radar_station=XLA&clear_sky_clock_station_id=ValdOrPQ'>QC Val-d'or
     <option value='graphic%2Esky?country=cda&prov=SK&location_name=Cypress+Hills+–+West+Block&latitude=49%2E60&longitude=-109%2E92&locations=49%2E60,-109%2E92&weather_station=29&radar_station=XBU&clear_sky_clock_station_id=CHDSPWBAB'>SK Cypress Hills – West Block
     <option value='graphic%2Esky?country=cda&prov=SK&location_name=Grasslands+National+Park+–+East&latitude=49%2E07&longitude=-106%2E53&locations=49%2E07,-106%2E53&weather_station=28&radar_station=XBU&clear_sky_clock_station_id=GNPEBSK'>SK Grasslands National Park – East
     <option value='graphic%2Esky?country=cda&prov=SK&location_name=Grasslands+National+Park+–+West&latitude=49%2E18&longitude=-107%2E71&locations=49%2E18,-107%2E71&weather_station=28&radar_station=XBU&clear_sky_clock_station_id=GrssNPESK'>SK Grasslands National Park – West
     <option value='graphic%2Esky?country=cda&prov=SK&location_name=La+Ronge&latitude=55%2E11&longitude=-105%2E28&locations=55%2E11,-105%2E28&weather_station=38&radar_station=&clear_sky_clock_station_id=LRngSK'>SK La Ronge
     <option value='graphic%2Esky?country=cda&prov=SK&location_name=Prince+Albert&latitude=53%2E20&longitude=-105%2E73&locations=53%2E20,-105%2E73&weather_station=27&radar_station=XRA&clear_sky_clock_station_id=PrncAlbrtSK'>SK Prince Albert
     <option value='graphic%2Esky?country=cda&prov=SK&location_name=Regina&latitude=50%2E45&longitude=-104%2E62&locations=50%2E45,-104%2E62&weather_station=32&radar_station=XBE&clear_sky_clock_station_id=Regina'>SK Regina
     <option value='graphic%2Esky?country=cda&prov=SK&location_name=Saskatoon&latitude=52%2E13&longitude=-106%2E69&locations=52%2E13,-106%2E69&weather_station=40&radar_station=XRA&clear_sky_clock_station_id=Saskatoon'>SK Saskatoon
     <option value='graphic%2Esky?country=cda&prov=SK&location_name=Swift+Current&latitude=50%2E29&longitude=-107%2E79&locations=50%2E29,-107%2E79&weather_station=41&radar_station=XBU&clear_sky_clock_station_id=SwtCrnSK'>SK Swift Current
     <option value='graphic%2Esky?country=cda&prov=YT&location_name=Whitehorse&latitude=60%2E72&longitude=-135%2E06&locations=60%2E72,-135%2E06&weather_station=16&radar_station=&clear_sky_clock_station_id=WhthrsYK'>YT Whitehorse
     <option value=''>---U.S.A.---
     <option value='graphic%2Esky?country=us&location_name=Anchorage&latitude=61%2E20&longitude=-149%2E88&locations=61%2E20,-149%2E88&radar_station=AHG&clear_sky_clock_station_id=AnchorAK'>AK Anchorage
     <option value='graphic%2Esky?country=us&location_name=Phoenix&latitude=33%2E43&longitude=-112%2E09&locations=33%2E43,-112%2E09&radar_station=IWA&clear_sky_clock_station_id=Phoenix'>AZ Phoenix
     <option value='graphic%2Esky?country=us&location_name=Los+Angeles&latitude=34%2E03&longitude=-118%2E31&locations=34%2E03,-118%2E31&radar_station=SOX&clear_sky_clock_station_id=LAXCA'>CA Los Angeles
     <option value='graphic%2Esky?country=us&location_name=San+Diego&latitude=32%2E74&longitude=-117%2E21&locations=32%2E74,-117%2E21&radar_station=NKX&clear_sky_clock_station_id=SanDiego'>CA San Diego
     <option value='graphic%2Esky?country=us&location_name=San+Francisco&latitude=37%2E77&longitude=-122%2E43&locations=37%2E77,-122%2E43&radar_station=MUX&clear_sky_clock_station_id=SanFranCA'>CA San Francisco
     <option value='graphic%2Esky?country=us&location_name=Denver&latitude=39%2E73&longitude=-104%2E97&locations=39%2E73,-104%2E97&radar_station=FTG&clear_sky_clock_station_id=DenverCO'>CO Denver
     <option value='graphic%2Esky?country=us&location_name=Washington&latitude=38%2E89&longitude=-77%2E03&locations=38%2E89,-77%2E03&radar_station=LWX&clear_sky_clock_station_id=WashingtonDC'>DC Washington
     <option value='graphic%2Esky?country=us&location_name=Jacksonville&latitude=30%2E32&longitude=-81%2E67&locations=30%2E32,-81%2E67&radar_station=JAX&clear_sky_clock_station_id=JacksonFL'>FL Jacksonville
     <option value='graphic%2Esky?country=us&location_name=Miami&latitude=25%2E76&longitude=-80%2E20&locations=25%2E76,-80%2E20&radar_station=AMX&clear_sky_clock_station_id=MiamiFL'>FL Miami
     <option value='graphic%2Esky?country=us&location_name=Atlanta&latitude=33%2E75&longitude=-84%2E39&locations=33%2E75,-84%2E39&radar_station=FFC&clear_sky_clock_station_id=AtlantaGA'>GA Atlanta
     <option value='graphic%2Esky?country=us&location_name=Des+Moines&latitude=41%2E58&longitude=-93%2E60&locations=41%2E58,-93%2E60&radar_station=DMX&clear_sky_clock_station_id=DeMoinesIA'>IA Des Moines
     <option value='graphic%2Esky?country=us&location_name=Boise&latitude=43%2E60&longitude=-116%2E24&locations=43%2E60,-116%2E24&radar_station=CBX&clear_sky_clock_station_id=BoiseID'>ID Boise
     <option value='graphic%2Esky?country=us&location_name=Chicago&latitude=41%2E88&longitude=-87%2E63&locations=41%2E88,-87%2E63&radar_station=LOT&clear_sky_clock_station_id=Chicago'>IL Chicago
     <option value='graphic%2Esky?country=us&location_name=Indianapolis&latitude=39%2E76&longitude=-86%2E16&locations=39%2E76,-86%2E16&radar_station=IND&clear_sky_clock_station_id=Indianapolis'>IN Indianapolis
     <option value='graphic%2Esky?country=us&location_name=Kansas+City&latitude=39%2E11&longitude=-94%2E63&locations=39%2E11,-94%2E63&radar_station=EAX&clear_sky_clock_station_id=KansasKA'>KS Kansas City
     <option value='graphic%2Esky?country=us&location_name=New+Orleans&latitude=29%2E98&longitude=-90%2E09&locations=29%2E98,-90%2E09&radar_station=LIX&clear_sky_clock_station_id=NewOrlLA'>LA New Orleans
     <option value='graphic%2Esky?country=us&location_name=Boston&latitude=42%2E36&longitude=-71%2E06&locations=42%2E36,-71%2E06&radar_station=BOX&clear_sky_clock_station_id=Boston'>MA Boston
     <option value='graphic%2Esky?country=us&location_name=Bangor&latitude=44%2E80&longitude=-68%2E78&locations=44%2E80,-68%2E78&radar_station=CBW&clear_sky_clock_station_id=BngrME'>ME Bangor
     <option value='graphic%2Esky?country=us&location_name=Detroit&latitude=42%2E35&longitude=-83%2E07&locations=42%2E35,-83%2E07&radar_station=DTX&clear_sky_clock_station_id=Detroit'>MI Detroit
     <option value='graphic%2Esky?country=us&location_name=Minneapolis&latitude=44%2E98&longitude=-93%2E27&locations=44%2E98,-93%2E27&radar_station=MPX&clear_sky_clock_station_id=Minneapolis'>MN Minneapolis
     <option value='graphic%2Esky?country=us&location_name=St%2E+Louis&latitude=38%2E62&longitude=-90%2E23&locations=38%2E62,-90%2E23&radar_station=LSX&clear_sky_clock_station_id=StLouisMO'>MO St. Louis
     <option value='graphic%2Esky?country=us&location_name=Great+Falls&latitude=47%2E51&longitude=-111%2E29&locations=47%2E51,-111%2E29&radar_station=TFX&clear_sky_clock_station_id=GrFallsMT'>MT Great Falls
     <option value='graphic%2Esky?country=us&location_name=Charlotte&latitude=35%2E21&longitude=-80%2E82&locations=35%2E21,-80%2E82&radar_station=GSP&clear_sky_clock_station_id=CharlotteNC'>NC Charlotte
     <option value='graphic%2Esky?country=us&location_name=Bismark&latitude=46%2E81&longitude=-100%2E77&locations=46%2E81,-100%2E77&radar_station=BIS&clear_sky_clock_station_id=BismarkND'>ND Bismark
     <option value='graphic%2Esky?country=us&location_name=Omaha&latitude=41%2E25&longitude=-96%2E01&locations=41%2E25,-96%2E01&radar_station=OAX&clear_sky_clock_station_id=OmahaNE'>NE Omaha
     <option value='graphic%2Esky?country=us&location_name=Albuquerque&latitude=35%2E10&longitude=-106%2E61&locations=35%2E10,-106%2E61&radar_station=ABX&clear_sky_clock_station_id=AlbuqNM'>NM Albuquerque
     <option value='graphic%2Esky?country=us&location_name=Socorro&latitude=34%2E06&longitude=-106%2E91&locations=34%2E06,-106%2E91&radar_station=ABX&clear_sky_clock_station_id=EtscornNM'>NM Socorro
     <option value='graphic%2Esky?country=us&location_name=Buffalo&latitude=42%2E89&longitude=-78%2E87&locations=42%2E89,-78%2E87&radar_station=BUF&clear_sky_clock_station_id=BuffaloNY'>NY Buffalo
     <option value='graphic%2Esky?country=us&location_name=New+York&latitude=40%2E76&longitude=-73%2E98&locations=40%2E76,-73%2E98&radar_station=OKX&clear_sky_clock_station_id=NYCNY'>NY New York
     <option value='graphic%2Esky?country=us&location_name=Cleveland&latitude=41%2E49&longitude=-81%2E70&locations=41%2E49,-81%2E70&radar_station=CLE&clear_sky_clock_station_id=Cleveland'>OH Cleveland
     <option value='graphic%2Esky?country=us&location_name=Columbus&latitude=39%2E95&longitude=-82%2E97&locations=39%2E95,-82%2E97&radar_station=ILN&clear_sky_clock_station_id=Columbus'>OH Columbus
     <option value='graphic%2Esky?country=us&location_name=Oklahoma+City&latitude=35%2E44&longitude=-97%2E53&locations=35%2E44,-97%2E53&radar_station=TLX&clear_sky_clock_station_id=OklahomaOK'>OK Oklahoma City
     <option value='graphic%2Esky?country=us&location_name=Portland&latitude=45%2E52&longitude=-122%2E64&locations=45%2E52,-122%2E64&radar_station=RTX&clear_sky_clock_station_id=PortOR'>OR Portland
     <option value='graphic%2Esky?country=us&location_name=Philadelphia&latitude=39%2E94&longitude=-75%2E16&locations=39%2E94,-75%2E16&radar_station=DIX&clear_sky_clock_station_id=Philadelphia'>PA Philadelphia
     <option value='graphic%2Esky?country=us&location_name=Charleston&latitude=32%2E77&longitude=-79%2E95&locations=32%2E77,-79%2E95&radar_station=CLX&clear_sky_clock_station_id=ChrlstnSC'>SC Charleston
     <option value='graphic%2Esky?country=us&location_name=Sioux+Falls&latitude=43%2E54&longitude=-96%2E73&locations=43%2E54,-96%2E73&radar_station=FSD&clear_sky_clock_station_id=SxFllsSD'>SD Sioux Falls
     <option value='graphic%2Esky?country=us&location_name=Memphis&latitude=35%2E13&longitude=-90%2E03&locations=35%2E13,-90%2E03&radar_station=NQA&clear_sky_clock_station_id=MemphisTN'>TN Memphis
     <option value='graphic%2Esky?country=us&location_name=Nashville&latitude=36%2E17&longitude=-86%2E76&locations=36%2E17,-86%2E76&radar_station=OHX&clear_sky_clock_station_id=NashTE'>TN Nashville
     <option value='graphic%2Esky?country=us&location_name=Austin&latitude=30%2E26&longitude=-97%2E74&locations=30%2E26,-97%2E74&radar_station=GRK&clear_sky_clock_station_id=AustinTX'>TX Austin
     <option value='graphic%2Esky?country=us&location_name=Dallas&latitude=32%2E76&longitude=-96%2E80&locations=32%2E76,-96%2E80&radar_station=FWS&clear_sky_clock_station_id=Dallas'>TX Dallas
     <option value='graphic%2Esky?country=us&location_name=El+Paso&latitude=31%2E77&longitude=-106%2E49&locations=31%2E77,-106%2E49&radar_station=EPZ&clear_sky_clock_station_id=ElPTX'>TX El Paso
     <option value='graphic%2Esky?country=us&location_name=Houston&latitude=29%2E75&longitude=-95%2E38&locations=29%2E75,-95%2E38&radar_station=HGX&clear_sky_clock_station_id=Houston'>TX Houston
     <option value='graphic%2Esky?country=us&location_name=San+Antonio&latitude=29%2E40&longitude=-98%2E50&locations=29%2E40,-98%2E50&radar_station=EWX&clear_sky_clock_station_id=SanAnTX'>TX San Antonio
     <option value='graphic%2Esky?country=us&location_name=Salt+Lake+City&latitude=40%2E76&longitude=-111%2E94&locations=40%2E76,-111%2E94&radar_station=MTX&clear_sky_clock_station_id=SaltLakeUT'>UT Salt Lake City
     <option value='graphic%2Esky?country=us&location_name=Norfolk&latitude=36%2E86&longitude=-76%2E29&locations=36%2E86,-76%2E29&radar_station=AKQ&clear_sky_clock_station_id=TidewaterVA'>VA Norfolk
     <option value='graphic%2Esky?country=us&location_name=Seattle&latitude=47%2E60&longitude=-122%2E32&locations=47%2E60,-122%2E32&radar_station=ATX&clear_sky_clock_station_id=Seattle'>WA Seattle
     <option value='graphic%2Esky?country=us&location_name=Spokane&latitude=47%2E67&longitude=-117%2E42&locations=47%2E67,-117%2E42&radar_station=OTX&clear_sky_clock_station_id=SpknWA'>WA Spokane
     <option value=''>---United Kingdom---
     <option value='graphic%2Esky?country=uk&location_name=Aberdeen&latitude=57%2E14&longitude=-2%2E11&locations=57%2E14,-2%2E11&weather_station=3091&'>Aberdeen
     <option value='graphic%2Esky?country=uk&location_name=Belfast&latitude=54%2E59&longitude=-5%2E93&locations=54%2E59,-5%2E93&weather_station=3917&'>Belfast
     <option value='graphic%2Esky?country=uk&location_name=Cardiff&latitude=51%2E48&longitude=-3%2E18&locations=51%2E48,-3%2E18&weather_station=3716&'>Cardiff
     <option value='graphic%2Esky?country=uk&location_name=Derry&latitude=55%2E03&longitude=-7%2E31&locations=55%2E03,-7%2E31&weather_station=3907&'>Derry
     <option value='graphic%2Esky?country=uk&location_name=Dumfries&latitude=55%2E07&longitude=-3%2E62&locations=55%2E07,-3%2E62&weather_station=3153&'>Dumfries
     <option value='graphic%2Esky?country=uk&location_name=Edinburgh&latitude=55%2E96&longitude=-3%2E19&locations=55%2E96,-3%2E19&weather_station=3166&'>Edinburgh
     <option value='graphic%2Esky?country=uk&location_name=Glasgow&latitude=55%2E86&longitude=-4%2E26&locations=55%2E86,-4%2E26&weather_station=3134&'>Glasgow
     <option value='graphic%2Esky?country=uk&location_name=Inverness&latitude=57%2E47&longitude=-4%2E23&locations=57%2E47,-4%2E23&weather_station=3063&'>Inverness
     <option value='graphic%2Esky?country=uk&location_name=Leeds&latitude=53%2E79&longitude=-1%2E54&locations=53%2E79,-1%2E54&weather_station=3344&'>Leeds
     <option value='graphic%2Esky?country=uk&location_name=Lerwick&latitude=60%2E15&longitude=-1%2E15&locations=60%2E15,-1%2E15&weather_station=3005&'>Lerwick
     <option value='graphic%2Esky?country=uk&location_name=Liverpool&latitude=53%2E40&longitude=-3%2E00&locations=53%2E40,-3%2E00&weather_station=3316&'>Liverpool
     <option value='graphic%2Esky?country=uk&location_name=London&latitude=51%2E50&longitude=-0%2E13&locations=51%2E50,-0%2E13&weather_station=3772&'>London
     <option value='graphic%2Esky?country=uk&location_name=Newcastle&latitude=54%2E98&longitude=-1%2E61&locations=54%2E98,-1%2E61&weather_station=3238&'>Newcastle
     <option value='graphic%2Esky?country=uk&location_name=Oxford&latitude=51%2E75&longitude=-1%2E26&locations=51%2E75,-1%2E26&weather_station=3469&'>Oxford
     <option value='graphic%2Esky?country=uk&location_name=Peterborough&latitude=52%2E57&longitude=-0%2E24&locations=52%2E57,-0%2E24&weather_station=3462&'>Peterborough
     <option value='graphic%2Esky?country=uk&location_name=Plymouth&latitude=50%2E37&longitude=-4%2E14&locations=50%2E37,-4%2E14&weather_station=3827&'>Plymouth
     <option value='graphic%2Esky?country=uk&location_name=Stornoway&latitude=58%2E21&longitude=-6%2E38&locations=58%2E21,-6%2E38&weather_station=3026&'>Stornoway
   </select>
  </form>
  

<P>The information is customized for your location, and includes:
<ul>
 <li>current weather conditions and forecast.
 <li>weather radar images (usually not available in the far north).
 <li>satellite images showing cloud cover (not available in the far north).
 <li>the Clear Sky Clock (Canada and the US only).
 <li>positions of the Sun, Moon, and planets.
 <li>selected asteroids and bright comets. 
 <li>a current diary of sky phenomena.
 <li>Messier/Caldwell objects appearing above the horizon (highest altitude first).
 <li>a planisphere showing the night sky.
 <li>the current observation window, showing the Sun's rise-set times, times of twilight, and the Moon's rise-set times.
 <li>lunar libration, for today and the upcoming month.
 <li>positions of the Galilean satellites of Jupiter.
 <li>upcoming meteor showers (if any).
 <li>auroral activity level (if above your configured threshold).
 <li>upcoming occultations at your location (if any).
 <li>phenomena related to the Galilean satellites of Jupiter (if any).
</ul>

<P>These links will help you look up the data you need when filling out the form:
<ul>
 <li>latitude and longitude: the quickest way is to click the <em>Use current location</em> button, which will use your browser's default location. 
 Otherwise, use <a href='https://www.google.ca/maps/'>Google Maps</a>: right-click on a location, and select "What's here?", to see the corresponding latitude and longitude.
 <li>radar id for US locations: use <a href='http://radar.weather.gov/radar.php?rid=lwx&product=N0R&overlay=11101111&loop=no'>NOAA/NWS</a>.
 <li>clear sky clock: the quickest way is to click the <em>Use lat/long</em> button, 
 and the system will find the station that's nearest to the latitude and longitude location you've 
 already entered in the form. 
 Otherwise, use <a href='http://www.cleardarksky.com/'>cleardarksky.com</a> to infer the station id.
</ul>

<P>If the meaning of an item in the form below is unclear to you, try and hover your mouse over the item; some items have a tooltip.
Required items are marked with *.

<P>
 <form method='GET' action='graphic.sky' class='user-input-small' id='input_form'>
  <table>
    
     <tr><td>Country*:<td>
         <select name='country' id='country' required>
           <option value='cda'>Canada
           <option value='us'>USA 
           <option value='uk'>United Kingdom 
           <option value='other'>Other (no weather data) 
         </select>
     <tr><td>Location name*:<td><input type='text' required id='location_name' name='location_name' value='Ottawa' title='Name of your observing site'>
     <tr><td>Latitude*:<td><input type='text' id='latitude' required name='latitude' value='45.40' title='Latitude in degrees'>
             <button type='button' id='lat_long_autofill' title='Let the browser fill in lat/long'>Use current location</button>
     <tr><td>Longitude*:<td><input type='text' id='longitude' required name='longitude' value='-75.66' title='Longitude in degrees. Negative west of Greenwich'>
     <tr id='row_clear_sky_clock_station_id'>
       <td>Clear Sky Chart station id*
       <td>
         <input name='clear_sky_clock_station_id' id='clear_sky_clock_station_id' value='FLO' size='10' title='Station id, used by the Clear Sky Clock website, for your location'>
         <button type='button' id='csc_find_nearest_station' title='Use the lat/long input above to find the nearest Clear Sky Clock'>Use lat/long</button>
     <tr id='row_radar_station_cda'><td>Radar Station*:<td>
         <select id='radar_station_cda' name='radar_station' title='Nearest radar station. Id used by Environment Canada'>
           <option value='WHK'>AB - Carvel, near Edmonton
           <option value='WHN'>AB - Jimmy Lake, near Cold Lake 
           <option value='XBU'>AB - Schuler, near Medicine Hat 
           <option value='WWW'>AB - Spirit River, near Grande Prairie 
           <option value='XSM'>AB - Strathmore, near Calgary 
           <option value='WUJ'>BC - Aldergrove, near Vancouver 
           <option value='XPG'>BC - Prince George 
           <option value='XSS'>BC - Silver Star Mountain, near Vernon 
           <option value='XSI'>BC - Victoria 
           <option value='XFW'>MB - Foxwarren, near Brandon 
           <option value='XWL'>MB - Woodlands, near Winnipeg 
           <option value='XNC'>NB - Chipman, near Frederiction 
           <option value='WTP'>NL - Holyrood, near St. John's 
           <option value='XME'>NL - Marble Mountain, near Corner Brook 
           <option value='XGO'>NS - Halifax 
           <option value='XMB'>NS - Marion Bridge, near Sydney 
           <option value='WBI'>ON - Britt, near Sudbury 
           <option value='XDR'>ON - Dryden 
           <option value='WSO'>ON - Exeter, near London 
           <option value='XFT' selected>ON - Franktown, near Ottawa 
           <option value='WKR'>ON - King City, near Toronto 
           <option value='WGJ'>ON - Montreal River, near Sault Ste. Marie 
           <option value='XTI'>ON - Northeast Ontario, near Timmins 
           <option value='XNI'>ON - Superior West, near Thunder Bay 
           <option value='WMB'>QC - Lac Castor, near Saguenay 
           <option value='XLA'>QC - Landrienne, near Rouyn-Noranda 
           <option value='WMN'>QC - McGill, near Montréal 
           <option value='XAM'>QC - Val d'Irène, near Mont-Joli 
           <option value='WVY'>QC - Villeroy, Trois-Rivières 
           <option value='XBE'>SK - Bethune, near Regina 
           <option value='XRA'>SK - Radisson, near Saskatoon 
         </select>
     <tr id='row_radar_station_us' style='display:none;'><td>Radar Station*<td><input id='radar_station_us' name='radar_station' title='Nearest US Radar Station Id, eg TYX (NOAA)'  disabled="true">
     <tr id='row_prov_cda'><td>Province*
         <td>
         <select id='prov_cda' name='prov'>
           <option value='NL'>NL
           <option value='NS'>NS
           <option value='NB'>NB
           <option value='PE'>PE
           <option value='QC'>QC
           <option value='ON' selected>ON
           <option value='MB'>MB
           <option value='SK'>SK
           <option value='AB'>AB
           <option value='BC'>BC
           <option value='NU'>NU
           <option value='NT'>NT
           <option value='YT'>YT
         </select>
     <tr id='row_weather_station_cda'><td>Weather Station ID*<td><input id='weather_station_cda' name='weather_station' value='118' title='Environment Canada ID for the nearest weather station'>

     <tr id='row_weather_station_uk' style='display:none;'><td>Weather Station*
       <td>
       <select id='weather_station_uk' name='weather_station' title='Met Office ID for the nearest weather station' disabled="true">
        <option value='3091'>Aberdeen - Aberdeen Airport
        <option value='3080'>Aberdeenshire - Aboyne
        <option value='3088'>Aberdeenshire - Inverbervie
        <option value='3111'>Argyll and Bute - Campbeltown Airport
        <option value='3105'>Argyll and Bute - Islay Airport
        <option value='3100'>Argyll and Bute - Tiree
        <option value='3560'>Bedford - Bedford
        <option value='3660'>Buckinghamshire - High Wycombe
        <option value='3605'>Carmarthenshire - Pembrey Sands Samos
        <option value='99057'>Central Bedfordshire - Woburn
        <option value='3502'>Ceredigion - Aberporth
        <option value='3503'>Ceredigion - Trawsgoed
        <option value='3351'>Cheshire East - Rostherne No 2
        <option value='3305'>Conwy - Capel Curig
        <option value='3808'>Cornwall - Camborne
        <option value='3823'>Cornwall - Cardinham
        <option value='3809'>Cornwall - Culdrose
        <option value='3916'>County Antrim - Ballypatrick Forest
        <option value='3917'>County Antrim - Belfast International Airport
        <option value='3923'>County Armagh - Glenanne
        <option value='3911'>County Londonderry - Lough Fea Samos
        <option value='3907'>County Londonderry - Magilligan No 2
        <option value='3904'>County Tyrone - Castlederg
        <option value='3220'>Cumbria - Carlisle
        <option value='3227'>Cumbria - Great Dun Fell 2
        <option value='3212'>Cumbria - Keswick
        <option value='3225'>Cumbria - Shap
        <option value='3224'>Cumbria - Spadeadam
        <option value='3210'>Cumbria - St. Bees Head
        <option value='3214'>Cumbria - Walney Island
        <option value='3226'>Cumbria - Warcop
        <option value='3313'>Denbighshire - Rhyl
        <option value='3707'>Devon - Chivenor
        <option value='3840'>Devon - Dunkeswell Aerodrome
        <option value='3839'>Devon - Exeter Airport
        <option value='3844'>Devon - Exeter Airport 2
        <option value='99081'>Devon - North Wyke
        <option value='3862'>Dorset - Bournemouth Airport
        <option value='3857'>Dorset - Isle Of Portland
        <option value='3153'>Dumfries and Galloway - Dundrennan
        <option value='3162'>Dumfries and Galloway - Eskdalemuir
        <option value='3132'>Dumfries and Galloway - West Freugh (Esaws)
        <option value='3292'>East Riding of Yorkshire - Bridlington Mrsc
        <option value='3382'>East Riding of Yorkshire - Leconfield Sar
        <option value='3882'>East Sussex - Herstmonceux West End
        <option value='3166'>Edinburgh - Edinburgh/Gogarbank
        <option value='3684'>Essex - Andrewsfield
        <option value='3693'>Essex - Shoeburyness
        <option value='3171'>Fife - Leuchars
        <option value='3321'>Flintshire - Hawarden
        <option value='3647'>Gloucestershire - Little Rissington (Esaws)
        <option value='3772'>Greater London - Heathrow
        <option value='3672'>Greater London - Northolt
        <option value='3894'>Guernsey - Guernsey
        <option value='3405'>Gwynedd - Aberdaron
        <option value='3768'>Hampshire - Farnborough
        <option value='3749'>Hampshire - Middle Wallop
        <option value='3761'>Hampshire - Odiham
        <option value='3522'>Herefordshire - Hereford
        <option value='3520'>Herefordshire - Shobdon Saws
        <option value='3680'>Hertfordshire - Rothamsted
        <option value='3044'>Highland - Altnaharra Saws
        <option value='3041'>Highland - Aonach Mor
        <option value='3034'>Highland - Aultbea
        <option value='3063'>Highland - Aviemore
        <option value='3039'>Highland - Bealach Na Ba
        <option value='3031'>Highland - Loch Glascarnoch Saws
        <option value='3037'>Highland - Skye/Lusa (Samos)
        <option value='3010'>Highland - Sule Skerry (Maws)
        <option value='3047'>Highland - Tulloch Bridge
        <option value='3075'>Highland - Wick John O Groats Airport
        <option value='3302'>Isle of Anglesey - Valley
        <option value='3866'>Isle of Wight - St Catherines Pt.
        <option value='3803'>Isles of Scilly - Scilly St Marys
        <option value='3895'>Jersey - Jersey
        <option value='3784'>Kent - Gravesend-Broadness
        <option value='3796'>Kent - Langdon Bay
        <option value='3797'>Kent - Manston
        <option value='99060'>Lancashire - Stonyhurst
        <option value='3391'>Lincolnshire - Coningsby
        <option value='3379'>Lincolnshire - Cranwell
        <option value='3385'>Lincolnshire - Donna Nook
        <option value='3469'>Lincolnshire - Holbeach
        <option value='3373'>Lincolnshire - Scampton
        <option value='3377'>Lincolnshire - Waddington
        <option value='3392'>Lincolnshire - Wainfleet
        <option value='3316'>Merseyside - Crosby
        <option value='3065'>Moray - Cairn Gorm Summit
        <option value='3066'>Moray - Kinloss
        <option value='3068'>Moray - Lossiemouth
        <option value='3023'>Na h-Eileanan Siar - South Uist Range
        <option value='3026'>Na h-Eileanan Siar - Stornoway
        <option value='3482'>Norfolk - Marham
        <option value='3488'>Norfolk - Weybourne
        <option value='3261'>North Yorkshire - Dishforth Airfield
        <option value='3281'>North Yorkshire - Fylingdales
        <option value='3257'>North Yorkshire - Leeming
        <option value='3266'>North Yorkshire - Linton On Ouse
        <option value='99142'>North Yorkshire - Scarborough
        <option value='3265'>North Yorkshire - Topcliffe
        <option value='3238'>Northumberland - Albemarle
        <option value='3240'>Northumberland - Boulmer
        <option value='3230'>Northumberland - Redesdale Camp (Samos)
        <option value='3354'>Nottinghamshire - Watnall
        <option value='3017'>Orkney Islands - Kirkwall
        <option value='3658'>Oxfordshire - Benson
        <option value='3649'>Oxfordshire - Brize Norton
        <option value='3604'>Pembrokeshire - Milford Haven C.B.
        <option value='3072'>Perth and Kinross - Cairnwell
        <option value='3144'>Perth and Kinross - Strathallan
        <option value='3462'>Peterborough - Wittering
        <option value='3827'>Plymouth - Mount Batten
        <option value='3410'>Powys - Lake Vyrnwy Saws
        <option value='3507'>Powys - Sennybridge
        <option value='3275'>Redcar and Cleveland - Loftus (Samos)
        <option value='3134'>Renfrewshire - Glasgow/Bishopton
        <option value='3158'>Scottish Borders - Charterhall
        <option value='3002'>Shetland Islands - Baltasound
        <option value='3008'>Shetland Islands - Fair Isle
        <option value='3014'>Shetland Islands - Foula
        <option value='3005'>Shetland Islands - Lerwick (S. Screen)
        <option value='3414'>Shropshire - Shawbury
        <option value='3710'>Somerset - Liscombe
        <option value='3853'>Somerset - Yeovilton
        <option value='3136'>South Ayrshire - Prestwick Rnas
        <option value='3628'>South Gloucestershire - Filton
        <option value='3155'>South Lanarkshire - Drumalbin
        <option value='3330'>Staffordshire - Leek
        <option value='3148'>Stirling - Glen Ogle
        <option value='3590'>Suffolk - Wattisham
        <option value='3769'>Surrey - Charlwood
        <option value='3781'>Surrey - Kenley
        <option value='3609'>Swansea - Mumbles Head
        <option value='3716'>Vale of Glamorgan - St-Athan
        <option value='3544'>Warwickshire - Church Lawford
        <option value='3535'>Warwickshire - Coleshill
        <option value='3876'>West Sussex - Shoreham
        <option value='3872'>West Sussex - Thorney Island
        <option value='3344'>West Yorkshire - Bingley Samos
        <option value='3746'>Wiltshire - Boscombe Down
        <option value='3743'>Wiltshire - Larkhill
        <option value='3740'>Wiltshire - Lyneham
        <option value='3529'>Worcestershire - Pershore
        <option value='3976'>Belmullet
        <option value='3980'>Malin Head
        <option value='3952'>Roches Point
        <option value='3204'>Ronaldsway
        <option value='3953'>Valentia Observatory
      </select>     
     
     
     
     <tr><td title='Translation is not complete.'>Preferred language:<td>
         <select name='preferred_lang' id='preferred_lang' title='Output in English or French'>
           <option value='e' selected>English
           <option value='f'>French 
         </select>
     <tr><td>Date-time:<td><input type='text' id='date_time' name='date_time' title='Example: 2016-03-01 18:00. Leave empty to use the current date and time.' >
     <tr><td>Time scale:<td><select name='time_scale' id='time_scale'><option title='Local Time (defined by your browser)' value='LT'>Local Time<option title='Universal Time (UT1, to be precise)' value='UT'>Universal Time<option title='Terrestrial Time (fundamental physics time)' value='TT'>Terrestrial Time</select>
     <tr><td>Show locations<td><input name='locations' id='locations' title='Latitude and longitude; separate with a semi-colon' value='45.510,-73.675;44.228,-76.492;45.255,-76.262' size='40'>
     <tr><td>Limiting Visual Mag:<td><input type='text' name='limiting_mag' id='limiting_mag' value='5.3' title='Higher number means show dimmer stars'>
     <tr><td>Limiting Mag Messier/Caldwell:<td><input type='text' name='limiting_mag_messiers' id='limiting_mag_messiers' value='11.0'>
     <tr><td>Limiting Mag Messier/Caldwell (planisphere):<td><input type='text' name='limiting_mag_messiers_planisphere' id='limiting_mag_messiers_planisphere' value='8.0'>
     <tr><td>Clouds: degrees on a side<td><input name='degrees_on_a_side' id='degrees_on_a_side' value='3'>
     <tr><td>Clouds: pixels on a side<td><input name='pixels_on_a_side' id='pixels_on_a_side' value='480'>
     <tr><td>Clouds: visible or infrared?
         <td>
         <select name='layer' id='layer'>
           <option value='auto_detect' title='Let the system decide' selected>Auto-detect
           <option value='visible' title='Day time'>Visible (day-time)
           <option value='ir' title='Night time'>IR (night-time)
         </select>
     <tr><td>Twilight
         <td>
         <select name='twilight' id='twilight'>
           <option value='-18' title='Sun 18 degrees below the horizon'>Astronomical
           <option value='-12' title='Sun 12 degrees below the horizon' selected>Nautical
           <option value='-6' title='Sun 6 degrees below the horizon'>Civil
         </select>
     <tr><td>Planisphere rotation angle:<td>
         <select name='planisphere_rotation_angle' id='planisphere_rotation_angle' title='Controls the default orientation of the planisphere'>
           <option>0
           <option>45
           <option>90
           <option>135
           <option>180
           <option>225
           <option>270
           <option>315
         </select>
     <tr title='Controls when you see aurora activity data'><td>Minimum auroral activity level:<td>
         <select name='aurora_min_activity_level' id='aurora_min_activity_level' title='Kp, planetary mean'>
           <option value='-1'>Default for your latitude
           <option>0
           <option>1
           <option>2
           <option>3
           <option>4
           <option>5
           <option>6
           <option>7
           <option>8
           <option>9
         </select>
     <tr><td>Occultations: num days to look ahead<td><input name='occultations_num_days_ahead' id='occultations_num_days_ahead' value='10'>
     <tr><td>Occultations: min magnitude<td><input name='occultations_min_mag' id='occultations_min_mag' value='6'>
     
     <tr><td>Exclude clouds?<td><input name='exclude_clouds' id='exclude_clouds'' type='checkbox' value='1'>
     <tr><td>Exclude radar?<td><input name='exclude_radar' id='exclude_radar' type='checkbox' value='1'>
     <tr><td>Exclude ecliptic?<td><input name='exclude_ecliptic' id='exclude_ecliptic' type='checkbox' value='1'>
     <tr><td>Exclude Sun, Moon, and Planets?<td><input name='exclude_sun_moon_planets' id='exclude_sun_moon_planets' type='checkbox' value='1'>
     <tr><td>Exclude sky diary?<td><input name='exclude_sky_diary' id='exclude_sky_diary' type='checkbox' value='1'>
     <tr><td>Exclude minor planets?<td><input name='exclude_minor_planets' id='exclude_minor_planets' type='checkbox' value='1'>
     <tr><td>Exclude comets?<td><input name='exclude_comets' id='exclude_comets' type='checkbox' value='1'>
     <tr><td>Exclude Messier objects?<td><input name='exclude_messiers' id='exclude_messiers' type='checkbox' value='1'>
     <tr><td>Exclude Caldwell objects?<td><input name='exclude_caldwells' id='exclude_caldwells' type='checkbox' value='1'>
     <tr><td>Exclude meteor showers?<td><input name='exclude_meteor_showers' id='exclude_meteor_showers' type='checkbox' value='1'>
     <tr><td>Exclude Galilean satellites?<td><input name='exclude_galilean_satellites' id='exclude_galilean_satellites' type='checkbox' value='1'>
     <tr><td>Exclude observation window?<td><input name='exclude_observation_window' id='exclude_observation_window' type='checkbox' value='1'>
     <tr><td>Exclude current weather?<td><input name='exclude_current_weather_conditions' id='exclude_current_weather_conditions' type='checkbox' value='1'>
     <tr><td>Exclude weather forecast?<td><input name='exclude_weather_forecast' id='exclude_weather_forecast' type='checkbox' value='1'>
     <tr><td>Exclude libration?<td><input name='exclude_libration' id='exclude_libration' type='checkbox' value='1'>
     <tr><td>Exclude planisphere?<td><input name='exclude_planisphere' id='exclude_planisphere' type='checkbox' value='1'>
     <tr><td>Exclude Clear Sky Clock?<td><input name='exclude_clear_sky_clock' id='exclude_clear_sky_clock' type='checkbox' value='1'>
     <tr><td>Exclude aurora?<td><input name='exclude_aurora' id='exclude_aurora' type='checkbox' value='1'>
     <tr><td>Exclude occultations?<td><input name='exclude_occultations' id='exclude_occultations' type='checkbox' value='1'>
     <tr><td>Exclude Jupiter satellite phenomena?<td><input name='exclude_jupiter_satellite_phenomena' id='exclude_jupiter_satellite_phenomena' type='checkbox' value='1'>
     
     <tr><td colspan='2' style="text-align:center">
      <input type='button' value='Exclude none' id='exclude_none'>
      <input type='button' value='Exclude all' id='exclude_all'>
      <input type='submit' value='SHOW the sky tonight'>
      <a href='form.sky'>Reset</a>
  </table>
 </form>

 <P><b>There is a defect related to time zones.</b> 
Example: if you are in Ottawa, then you can view data for Washington DC, since those two cities are in the same time zone.
If you are in Ottawa and you view data for Los Angeles, however, then the times are going to be messed up.
<em>Everything works fine only if you stay in the same time zone used by your browser.</em>
Most of the time, this defect won't bother you, because you're usually concerned with nearby locations.
(This problem is caused by how Javascript treats dates and times.
There are ways around this problem, but I haven't done anything about it yet. 
If you need a work-around, then just temporarily change your computer's time zone setting.) 

 
 <a id='clouds_large_graphic'></a>
 <P> <b>The satellite imagery is particularly useful, because it shows detailed cloud information.</b> 
 The satellite image generated by the above form centers on a given location, and has a small area.
 It's also interesting to see the clouds over a large area of the country. 
 The images are especially dramatic when the Sun is low in the sky, and produces longer shadows. 
 Here's a <a href='../satellite/form.sky'>form for showing larger cloud images</a>.  

<P>Credits:
<ul>
 <li><em><a href='http://www.willbell.com/math/mc1.htm'>Astronomical Algorithms</a></em>, Jean Meeus: many algorithms used here are taken from this book.  
 <li><a href='http://weather.gc.ca/canada_e.html'>Environment Canada</a> : current conditions, forecasts, and radar images.
 <li><a href='http://radar.weather.gov/'>NOAA/NWS</a> : current US weather conditions, forecasts, and radar images.
 <li><a href='https://mesonet.agron.iastate.edu/'>Iowa State University</a>: server for NOAA images captured by the GOES satellite (clouds).
 <li><a href='http://rasc.ca/handbook'>Observer's Handbook</a> of the Royal Astronomical Society of Canada: comments on the Messier objects, and osculating orbital elements for planets.
 <li>Yale Bright Star catalog, <a href='http://cdsarc.u-strasbg.fr/viz-bin/Cat?V/50'>revision 5</a>.
 <li>Allen Rahill (Canadian Meteorological Center) and Attilla Danko: <a href='http://www.cleardarksky.com/csk/'>Clear Sky Clock</a> images.
 <li><a href='http://www.imo.net'>International Meteor Organization</a> : Meteor shower <a href='http://www.imo.net/files/data/vmdb/vmdbrad.txt'>data</a>.
 <li><a href='https://www.ast.cam.ac.uk/~jds/'>British Astronomical Association</a>, Comet section : current comet data.
 <li><a href='http://ssd.jpl.nasa.gov/sbdb.cgi#top'>Jet Propulsion Laboratory</a> : Orbital elements for asteroids and comets.
 <li><a href='http://www.geomag.bgs.ac.uk/education/poles.html'>British Geological Survey</a> : Position of the geomagnetic north pole.
 <li><a href='http://www.swpc.noaa.gov/products/station-k-and-indices'>NOAA Space Weather</a> : Auroral activity (Kp).
 <li><a href='http://www.lunar-occultations.com/iota/iotandx.htm'>International Occultation Timing Association</a> : occultation predictions.
 <li><a href='http://aa.usno.navy.mil/software/mica/micainfo.php'>US Naval Observatory</a>, Multiyear Interactive Computer Almanac (MICA, v2.2.2) : sky diary, phenomena for Galilean satellites.
 <li><a href='http://www.metoffice.gov.uk/datapoint'>Met Office</a> : current UK weather conditions, forecast, and radar images
 <li><a href='http://www.eumetsat.int/website/home/index.html'>EUMETSAT</a> : Meteosat satellite images for Europe (via the Met Office)
</ul>

<P>Code last updated on: ${initParam.lastUpdatedOn}. 
<P>Created by <a href='mailto:webmaster@javapractices.com'>John O'Hanley</a> (Ottawa, Canada).  
<P>Do you want to improve this site? Astronomytonight.net is open source (<a href='https://github.com/johanley/astro'>github</a>). 

<tags:analytics/>

</body>
</html>