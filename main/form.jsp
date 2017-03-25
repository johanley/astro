<!doctype html>
<html lang='<s:txt>en</s:txt>'>
<head>
 <meta charset='UTF-8'>
 <meta name="keywords" content="<s:txt>astronomy, ephemerides, low precision, amateur, canada, united states, united kingdom, weather, radar, clouds</s:txt>">
 <meta name="description" content="<s:txt>Data of interest to amateur astronomers. Includes weather for North America and the United Kingdom.</s:txt>">
 <meta name="viewport" content="width=device-width">
 <link rel="stylesheet" href="../css/styles.css<tags:ver/>" media="all">
 <title><s:txt>The sky tonight</s:txt></title>
</head>
 <script src='find_closest_clear_sky_clock.js<tags:ver/>'></script>
 <script src='../js/util.js<tags:ver/>'></script>
 <script src='../js/ephem.js<tags:ver/>'></script>
 <script>
  <tags:geolocation/> 
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
     this.innerHTML = '<s:txt>Scanning...</s:txt>';
     if (lat && long){
       var target = {};
       target.φ = lat; 
       target.λ = long;
       var callback = function(id){
         //console.log('The nearest station id: ' + id);
         document.getElementById('clear_sky_clock_station_id').value = id;
         document.getElementById('csc_find_nearest_station').innerHTML = '<s:txt>Find nearest</s:txt>';
       };
       find_closest_csc_url(target, callback);
     }
   };
   
   activate_lat_long_autofill(); //in geolocations.tag
   
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
     if ('cda' === country){
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
<h2><s:txt>Astronomy Tonight</s:txt>
 <span class='home'>
  <form style='display:inline;'>
   <select onChange="if (this.value) window.location.href=this.value" style='vertical-align:text-top;'>
    <option value='form.sky'>Language
     <option value='form.sky?lang=en'>English
     <option value='form.sky?lang=fr'>Français
   </select>
  </form>
 </span>
</h2>

<s:txt>This tool collects together, on a single page, basic information of interest to most amateur astronomers. 
It's most useful for residents of Canada, the US, and the UK. Outside of those areas, no weather data is shown.</s:txt>

<P><s:txt>Below you will find a large form that lets you customize the output data.
Its defaults are for Ottawa, the capital of Canada. 
As a convenience, here are the results for some selected locations (Canada, US, and UK only)</s:txt>: 

<P>
  <form >
   <select onChange="if (this.value) window.location.href=this.value">
    <tags:locations/>
   </select>
  </form>

<P><s:txt>The information is customized for your location, and includes</s:txt>:
<ul>
 <li><s:txt>current weather conditions and forecast.</s:txt>
 <li><s:txt>weather radar images (usually not available in the far north).</s:txt>
 <li><s:txt>satellite images showing cloud cover (not available in the far north).</s:txt>
 <li><s:txt>the Clear Sky Clock (Canada and the US only).</s:txt>
 <li><s:txt>positions of the Sun, Moon, and planets.</s:txt>
 <li><s:txt>selected asteroids and bright comets.</s:txt>
 <li><s:txt>a current diary of sky phenomena.</s:txt>
 <li><s:txt>Messier/Caldwell objects appearing above the horizon (highest altitude first).</s:txt>
 <li><s:txt>a planisphere showing the night sky.</s:txt>
 <li><s:txt>the current observation window, showing the Sun's rise-set times, times of twilight, and the Moon's rise-set times.</s:txt>
 <li><s:txt>lunar libration, for today and the upcoming month.</s:txt>
 <li><s:txt>positions of the Galilean satellites of Jupiter.</s:txt>
 <li><s:txt>upcoming meteor showers (if any).</s:txt>
 <li><s:txt>auroral activity level (if above your configured threshold).</s:txt>
 <li><s:txt>upcoming occultations at your location (if any).</s:txt>
 <li><s:txt>phenomena related to the Galilean satellites of Jupiter (if any).</s:txt>
</ul>

<P><s:txt>These links will help you look up the data you need when filling out the form</s:txt>:
<s:txt>
<ul>
 <li>latitude and longitude: the quickest way is to click the <em>Use current location</em> button, which will use your browser's default location. 
 Otherwise, use <a href='https://www.google.ca/maps/'>Google Maps</a>: right-click on a location, and select "What's here?", to see the corresponding latitude and longitude.
 <li>radar id for US locations: use <a href='http://radar.weather.gov/radar.php?rid=lwx&product=N0R&overlay=11101111&loop=no'>NOAA/NWS</a>.
 <li>clear sky clock: the quickest way is to click the <em>Use lat/long</em> button, 
 and the system will find the station that's nearest to the latitude and longitude location you've 
 already entered in the form. 
 Otherwise, use <a href='http://www.cleardarksky.com/'>cleardarksky.com</a> to infer the station id.
</ul>
</s:txt>

<P><s:txt>If the meaning of an item in the form below is unclear to you, try and hover your mouse over the item; some items have a tooltip.
Required items are marked with *.</s:txt>

<P>
 <form method='GET' action='graphic.sky' class='user-input-small' id='input_form'>
  <table>
    
     <tr><td><s:txt>Country</s:txt>*:<td>
         <select name='country' id='country' required>
           <option value='cda'><s:txt>Canada</s:txt>
           <option value='us'><s:txt>United States</s:txt>
           <option value='uk'><s:txt>United Kingdom</s:txt> 
           <option value='other'><s:txt>Other (no weather data)</s:txt>
         </select>
     <tr><td><s:txt>Location name</s:txt>*:<td><input type='text' required id='location_name' name='location_name' value='Ottawa' title='<s:txt>Name of your observing site</s:txt>'>
     <tr><td><s:txt>Latitude</s:txt>*:<td><input type='text' id='latitude' required name='latitude' value='45.40' title='<s:txt>Latitude in degrees</s:txt>'>
             <button type='button' id='lat_long_autofill' title='<s:txt>Let the browser fill in lat/long</s:txt>'><s:txt>Use current location</s:txt></button>
     <tr><td><s:txt>Longitude</s:txt>*:<td><input type='text' id='longitude' required name='longitude' value='-75.66' title='<s:txt>Longitude in degrees. Negative west of Greenwich</s:txt>'>
     <tr id='row_clear_sky_clock_station_id'>
       <td><s:txt>Clear Sky Chart station id</s:txt>*
       <td>
         <input name='clear_sky_clock_station_id' id='clear_sky_clock_station_id' value='FLO' size='10' title='<s:txt>Station id, used by the Clear Sky Clock website, for your location</s:txt>'>
         <button type='button' id='csc_find_nearest_station' title='<s:txt>Use the lat/long input above to find the nearest Clear Sky Clock</s:txt>'><s:txt>Use lat/long</s:txt></button>
     <tr id='row_radar_station_cda'><td><s:txt>Radar Station</s:txt>*:<td>
         <select id='radar_station_cda' name='radar_station' title='<s:txt>Nearest radar station. Id used by Environment Canada</s:txt>'>
           <option value='WHK'>AB - Carvel (Edmonton)
           <option value='WHN'>AB - Jimmy Lake (Cold Lake) 
           <option value='XBU'>AB - Schuler (Medicine Hat) 
           <option value='WWW'>AB - Spirit River (Grande Prairie) 
           <option value='XSM'>AB - Strathmore (Calgary) 
           <option value='WUJ'>BC - Aldergrove (Vancouver) 
           <option value='XPG'>BC - Prince George 
           <option value='XSS'>BC - Silver Star Mountain (Vernon) 
           <option value='XSI'>BC - Victoria 
           <option value='XFW'>MB - Foxwarren (Brandon) 
           <option value='XWL'>MB - Woodlands (Winnipeg) 
           <option value='XNC'>NB - Chipman (Frederiction) 
           <option value='WTP'>NL - Holyrood (St. John's) 
           <option value='XME'>NL - Marble Mountain (Corner Brook) 
           <option value='XGO'>NS - Halifax 
           <option value='XMB'>NS - Marion Bridge (Sydney) 
           <option value='WBI'>ON - Britt (Sudbury) 
           <option value='XDR'>ON - Dryden 
           <option value='WSO'>ON - Exeter (London) 
           <option value='XFT' selected>ON - Franktown (Ottawa) 
           <option value='WKR'>ON - King City (Toronto) 
           <option value='WGJ'>ON - Montreal River (Sault Ste. Marie) 
           <option value='XTI'>ON - Northeast Ontario (Timmins) 
           <option value='XNI'>ON - Superior West (Thunder Bay) 
           <option value='WMB'>QC - Lac Castor (Saguenay) 
           <option value='XLA'>QC - Landrienne (Rouyn-Noranda) 
           <option value='WMN'>QC - McGill (Montréal) 
           <option value='XAM'>QC - Val d'Irène (Mont-Joli) 
           <option value='WVY'>QC - Villeroy (Trois-Rivières) 
           <option value='XBE'>SK - Bethune (Regina) 
           <option value='XRA'>SK - Radisson (Saskatoon) 
         </select>
     <tr id='row_radar_station_us' style='display:none;'><td><s:txt>Radar Station</s:txt>*<td><input id='radar_station_us' name='radar_station' title='<s:txt>Nearest US Radar Station Id, eg TYX (NOAA)</s:txt>' disabled="true">
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
     <tr id='row_weather_station_cda'><td><s:txt>Weather Station ID</s:txt>*<td><input id='weather_station_cda' name='weather_station' value='118' title='<s:txt>Environment Canada ID for the nearest weather station</s:txt>'>

     <tr id='row_weather_station_uk' style='display:none;'><td><s:txt>Weather Station</s:txt>*
       <td>
       <select id='weather_station_uk' name='weather_station' title='<s:txt>Met Office name for the nearest weather station</s:txt>' disabled="true">
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
     
     <tr><td><s:txt>Date-time</s:txt>:<td><input type='text' id='date_time' name='date_time' title='<s:txt>Example: 2016-03-01 18:00. Leave empty to use the current date and time.</s:txt>' >
     <tr><td><s:txt>Time scale</s:txt>:<td>
       <select name='time_scale' id='time_scale'>
        <option title='<s:txt>Local Time (defined by your browser)</s:txt>' value='LT'><s:txt>Local Time</s:txt>
        <option title='<s:txt>Universal Time (UT1, to be precise)</s:txt>' value='UT'><s:txt>Universal Time</s:txt>
        <option title='<s:txt>Terrestrial Time (fundamental physics time)</s:txt>' value='TT'><s:txt>Terrestrial Time</s:txt>
       </select>
     <tr><td><s:txt>Show locations</s:txt><td><input name='locations' id='locations' title='<s:txt>Latitude and longitude; separate with a semi-colon</s:txt>' value='45.510,-73.675;44.228,-76.492;45.255,-76.262' size='40'>
     <tr><td><s:txt>Limiting Visual Mag</s:txt>:<td><input type='text' name='limiting_mag' id='limiting_mag' value='5.3' title='<s:txt>Higher number means show dimmer stars</s:txt>'>
     <tr><td><s:txt>Limiting Mag Messier/Caldwell</s:txt>:<td><input type='text' name='limiting_mag_messiers' id='limiting_mag_messiers' value='11.0'>
     <tr><td><s:txt>Limiting Mag Messier/Caldwell (planisphere)</s:txt>:<td><input type='text' name='limiting_mag_messiers_planisphere' id='limiting_mag_messiers_planisphere' value='8.0'>
     <tr><td><s:txt>Clouds: degrees on a side</s:txt><td><input name='degrees_on_a_side' id='degrees_on_a_side' value='3'>
     <tr><td><s:txt>Clouds: pixels on a side</s:txt><td><input name='pixels_on_a_side' id='pixels_on_a_side' value='480'>
     <tr><td><s:txt>Clouds: visible or infrared?</s:txt>
         <td>
         <select name='layer' id='layer'>
           <option value='auto_detect' title='<s:txt>Let the system decide</s:txt>' selected><s:txt>Auto-detect</s:txt>
           <option value='visible' title='<s:txt>Day time</s:txt>'><s:txt>Visible (day-time)</s:txt>
           <option value='ir' title='<s:txt>Night time</s:txt>'><s:txt>IR (night-time)</s:txt>
         </select>
     <tr><td><s:txt>Twilight</s:txt>
         <td>
         <select name='twilight' id='twilight'>
           <option value='-18' title='<s:txt>Sun 18 degrees below the horizon</s:txt>'><s:txt>Astronomical</s:txt>
           <option value='-12' title='<s:txt>Sun 12 degrees below the horizon</s:txt>' selected><s:txt>Nautical</s:txt>
           <option value='-6' title='<s:txt>Sun 6 degrees below the horizon</s:txt>'><s:txt>Civil</s:txt>
         </select>
     <tr><td><s:txt>Planisphere rotation angle</s:txt>:<td>
         <select name='planisphere_rotation_angle' id='planisphere_rotation_angle' title='<s:txt>Controls the default orientation of the planisphere</s:txt>'>
           <option>0
           <option>45
           <option>90
           <option>135
           <option>180
           <option>225
           <option>270
           <option>315
         </select>
     <tr title='<s:txt>Controls when you see aurora activity data</s:txt>'><td><s:txt>Minimum auroral activity level</s:txt>:<td>
         <select name='aurora_min_activity_level' id='aurora_min_activity_level' title='<s:txt>Kp, planetary mean</s:txt>'>
           <option value='-1'><s:txt>Default for your latitude</s:txt>
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
     <tr><td><s:txt>Occultations: num days to look ahead</s:txt><td><input name='occultations_num_days_ahead' id='occultations_num_days_ahead' value='10'>
     <tr><td><s:txt>Occultations: min magnitude</s:txt><td><input name='occultations_min_mag' id='occultations_min_mag' value='6'>
     
     <tr><td><s:txt>Exclude clouds</s:txt><td><input name='exclude_clouds' id='exclude_clouds'' type='checkbox' value='1'>
     <tr><td><s:txt>Exclude radar</s:txt><td><input name='exclude_radar' id='exclude_radar' type='checkbox' value='1'>
     <tr><td><s:txt>Exclude ecliptic</s:txt><td><input name='exclude_ecliptic' id='exclude_ecliptic' type='checkbox' value='1'>
     <tr><td><s:txt>Exclude Sun, Moon, and Planets</s:txt><td><input name='exclude_sun_moon_planets' id='exclude_sun_moon_planets' type='checkbox' value='1'>
     <tr><td><s:txt>Exclude sky diary</s:txt><td><input name='exclude_sky_diary' id='exclude_sky_diary' type='checkbox' value='1'>
     <tr><td><s:txt>Exclude minor planets</s:txt><td><input name='exclude_minor_planets' id='exclude_minor_planets' type='checkbox' value='1'>
     <tr><td><s:txt>Exclude comets</s:txt><td><input name='exclude_comets' id='exclude_comets' type='checkbox' value='1'>
     <tr><td><s:txt>Exclude Messier objects</s:txt><td><input name='exclude_messiers' id='exclude_messiers' type='checkbox' value='1'>
     <tr><td><s:txt>Exclude Caldwell objects</s:txt><td><input name='exclude_caldwells' id='exclude_caldwells' type='checkbox' value='1'>
     <tr><td><s:txt>Exclude meteor showers</s:txt><td><input name='exclude_meteor_showers' id='exclude_meteor_showers' type='checkbox' value='1'>
     <tr><td><s:txt>Exclude Galilean satellites</s:txt><td><input name='exclude_galilean_satellites' id='exclude_galilean_satellites' type='checkbox' value='1'>
     <tr><td><s:txt>Exclude observation window</s:txt><td><input name='exclude_observation_window' id='exclude_observation_window' type='checkbox' value='1'>
     <tr><td><s:txt>Exclude current weather</s:txt><td><input name='exclude_current_weather_conditions' id='exclude_current_weather_conditions' type='checkbox' value='1'>
     <tr><td><s:txt>Exclude weather forecast</s:txt><td><input name='exclude_weather_forecast' id='exclude_weather_forecast' type='checkbox' value='1'>
     <tr><td><s:txt>Exclude libration</s:txt><td><input name='exclude_libration' id='exclude_libration' type='checkbox' value='1'>
     <tr><td><s:txt>Exclude planisphere</s:txt><td><input name='exclude_planisphere' id='exclude_planisphere' type='checkbox' value='1'>
     <tr><td><s:txt>Exclude Clear Sky Clock</s:txt><td><input name='exclude_clear_sky_clock' id='exclude_clear_sky_clock' type='checkbox' value='1'>
     <tr><td><s:txt>Exclude aurora</s:txt><td><input name='exclude_aurora' id='exclude_aurora' type='checkbox' value='1'>
     <tr><td><s:txt>Exclude occultations</s:txt><td><input name='exclude_occultations' id='exclude_occultations' type='checkbox' value='1'>
     <tr><td><s:txt>Exclude Jupiter satellite phenomena</s:txt><td><input name='exclude_jupiter_satellite_phenomena' id='exclude_jupiter_satellite_phenomena' type='checkbox' value='1'>
     
     <tr><td colspan='2' style="text-align:center">
      <input type='button' value='<s:txt>Exclude none</s:txt>' id='exclude_none'>
      <input type='button' value='<s:txt>Exclude all</s:txt>' id='exclude_all'>
      <input type='submit' value='<s:txt>SHOW the sky tonight</s:txt>'>
      <a href='form.sky'><s:txt>Reset</s:txt></a>
  </table>
 </form>

<s:txt><P><b>There is a defect related to time zones.</b> 
Example: if you are in Ottawa, then you can view data for Washington DC, since those two cities are in the same time zone.
If you are in Ottawa and you view data for Los Angeles, however, then the times are going to be messed up.
<em>Everything works fine only if you stay in the same time zone used by your browser.</em>
Most of the time, this defect won't bother you, because you're usually concerned with nearby locations.
(This problem is caused by how Javascript treats dates and times.
There are ways around this problem, but I haven't done anything about it yet. 
If you need a work-around, then just temporarily change your computer's time zone setting.)</s:txt> 

 
 <a id='clouds_large_graphic'></a>
 <s:txt><P><b>The satellite imagery is particularly useful, because it shows detailed cloud information.</b> 
 The satellite image generated by the above form centers on a given location, and has a small area.
 It's also interesting to see the clouds over a large area of the country. 
 The images are especially dramatic when the Sun is low in the sky, and produces longer shadows. 
 Here's a <a href='../satellite/form.sky'>form for showing larger cloud images</a>.</s:txt>  

<P><s:txt>Credits</s:txt>:
<s:txt><ul>
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
 <li><a href='http://www.geomag.bgs.ac.uk/education/poles.html'>British Geological Survey</a> : position of the geomagnetic north pole.
 <li><a href='http://www.swpc.noaa.gov/products/station-k-and-indices'>NOAA Space Weather</a> : auroral activity (Kp).
 <li><a href='http://www.lunar-occultations.com/iota/iotandx.htm'>International Occultation Timing Association</a> : occultation predictions.
 <li><a href='http://aa.usno.navy.mil/software/mica/micainfo.php'>US Naval Observatory</a>, Multiyear Interactive Computer Almanac (MICA, v2.2.2) : sky diary, phenomena for Galilean satellites.
 <li><a href='http://www.metoffice.gov.uk/datapoint'>Met Office</a> : current UK weather conditions and forecast.
 <li><a href='http://www.eumetsat.int/website/home/index.html'>EUMETSAT</a> : Meteosat satellite images for Europe.
</ul></s:txt>

<P><s:txt>Code last updated on</s:txt>: ${initParam.lastUpdatedOn}.
 
<s:txt><P>Help to improve astronomytonight.net: 
 <ul>  
  <li>send corrections and suggestions to <a href='mailto:webmaster@javapractices.com'>John O'Hanley</a> (Ottawa, Canada).
  <li>contribute via <a href='https://github.com/johanley/astro'>github</a>.
 </ul></s:txt> 

<tags:analytics/>

</body>
</html>