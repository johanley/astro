<!doctype html>
<html lang='en'>
<head>
 <meta charset='UTF-8'>
 <meta name="keywords" content="astronomy">
 <meta name="description" content="Astronomical ephemerides for amateur astronomers">
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
  canvas, img {
    border-radius:10px; 
    background-color: rgb(25%, 25%, 25%);
    border-color: rgb(45%,45%,45%); 
    border-style: solid;
    border-width: 2px;
  }
  .cloud {
    width: 30%;
  }
  .radar {
    width: 36%;
    height: 100%; /*no effect*/
  }
  .weather {
    max-width:30%  
  }
  .sky {
    background-color: rgb(36,55,114);
  }
  @media all and (max-width: 750px) {
    .cloud {
      width: 95%;
    }
    .radar {
      width: 95%;
    }
    .weather {
      max-width:95%;  
    }
    #clear_sky_clock {
      width: 95%;
    }
  }
 </style>
 <title>The sky tonight</title>
 <script src='../js/ephem.js<tags:ver/>'></script>
 <script src='../js/util.js<tags:ver/>'></script>
 <script src='../js/general-graphics.js<tags:ver/>'></script>
 <script src='code.js<tags:ver/>'></script>
 <script src='precipitation.js<tags:ver/>'></script>
 <script src='../satellite/satellite.js<tags:ver/>'></script>
 <script src='aurora.js<tags:ver/>'></script>
 <script src='occultations.js<tags:ver/>'></script>
 <script src='jupiter_satellite_phenomena.js<tags:ver/>'></script>
 <c:if test="${param.country eq 'uk'}">
 <!-- 2 scripts for openlayers js mapping tool; tile mapping, used for UK data -->
 <script src="https://cdn.polyfill.io/v2/polyfill.min.js?features=requestAnimationFrame,Element.prototype.classList"></script>
 <script src="https://openlayers.org/en/v3.20.1/build/ol.js" type="text/javascript"></script>
 </c:if> 
 <script> 
  window.onload = function() {
    var formInput = UTIL.requestParams(window);
    show(formInput, "${initParam.metOfficeApiKey}", ${initParam.isDev}, "<s:txt>en</s:txt>");
  };
  window.onorientationchange = function() { 
    location.reload();
  };
 </script>
</head>
<body> 
 
 <table id='summary' class='report' style='margin-top:0;'>
    <tr><th><s:txt>Place</s:txt><th><s:txt>Time</s:txt><th title="<s:txt>Local Mean Sidereal Time</s:txt>"><s:txt>LMST</s:txt><th><s:txt>Lat</s:txt><th><s:txt>Long</s:txt><th><s:txt>UT - LT</s:txt>
 </table>
 <span style='white-space:nowrap;'  class='no-print'>
 <form id='date_time_controls' method='GET' action='' class='no-print' style="display:inline; margin:0.25em;">
   <input type='submit' name='go' value='&lt;' onclick="go_plus=0;"  class='no-print' >
   <input type='text' name='num_steps' id='num_steps' size='4'  class='no-print' >
   <select name='date_time_unit' id='date_time_unit'  class='no-print' >
     <option value='year'><s:txt>year</s:txt>
     <option value='month'><s:txt>month</s:txt>
     <option value='day' selected><s:txt>day</s:txt>
     <option value='hour'><s:txt>hour</s:txt>
     <option value='min'><s:txt>min</s:txt>
     <option value='sec'><s:txt>sec</s:txt>
   </select>
   <input type='submit' name='go' value='&gt;' onclick="go_plus=1;"  class='no-print' >
 </form>

  <form style="display:inline; margin:0.25em;">
   <select onChange="if (this.value) window.location.href=this.value">
     <option value=''><s:txt>Choose a location</s:txt>
     <option value=''><s:txt>---Canada---</s:txt>
     <option value='graphic%2Esky?country=cda&prov=AB&location_name=Banff&latitude=51%2E18&longitude=-115%2E57&weather_station=49&radar_station=XSM&clear_sky_clock_station_id=BnffAB'>AB Banff
     <option value='graphic%2Esky?country=cda&prov=AB&location_name=Calgary&latitude=51%2E03&longitude=-114%2E09&weather_station=52&radar_station=XSM&clear_sky_clock_station_id=Calgary'>AB Calgary
     <option value='graphic%2Esky?country=cda&prov=AB&location_name=Edmonton&latitude=53%2E54&longitude=-113%2E50&weather_station=50&radar_station=WHK&clear_sky_clock_station_id=Edmonton'>AB Edmonton
     <option value='graphic%2Esky?country=cda&prov=AB&location_name=Elk+Island+–+Beaver+Hills&latitude=53%2E58&longitude=-112%2E82&weather_station=63&radar_station=WHK&clear_sky_clock_station_id=Blackfoot'>AB Elk Island – Beaver Hills
     <option value='graphic%2Esky?country=cda&prov=AB&location_name=Fort+McMurray&latitude=56%2E73&longitude=-111%2E38&weather_station=20&radar_station=&clear_sky_clock_station_id=FtMcmrryAB'>AB Fort McMurray
     <option value='graphic%2Esky?country=cda&prov=AB&location_name=Grande+Prairie&latitude=55%2E17&longitude=-118%2E80&weather_station=31&radar_station=WWW&clear_sky_clock_station_id=NDgObAB'>AB Grande Prairie
     <option value='graphic%2Esky?country=cda&prov=AB&location_name=Jasper&latitude=52%2E87&longitude=-118%2E08&weather_station=70&radar_station=&clear_sky_clock_station_id=MlngLkAB'>AB Jasper
     <option value='graphic%2Esky?country=cda&prov=AB&location_name=Lethbridge&latitude=49%2E69&longitude=-112%2E85&weather_station=30&radar_station=XSM&clear_sky_clock_station_id=LethAB'>AB Lethbridge
     <option value='graphic%2Esky?country=cda&prov=AB&location_name=Red+Deer&latitude=52%2E26&longitude=-113%2E80&weather_station=29&radar_station=XSM&clear_sky_clock_station_id=RdDrAB'>AB Red Deer
     <option value='graphic%2Esky?country=cda&prov=BC&location_name=Cranbrook&latitude=49%2E51&longitude=-115%2E76&weather_station=77&radar_station=&clear_sky_clock_station_id=SprcRdObBC'>BC Cranbrook
     <option value='graphic%2Esky?country=cda&prov=BC&location_name=Kamloops&latitude=50%2E67&longitude=-120%2E29&weather_station=45&radar_station=XSS&clear_sky_clock_station_id=Kamloops'>BC Kamloops
     <option value='graphic%2Esky?country=cda&prov=BC&location_name=Penticton&latitude=49%2E49&longitude=-119%2E57&weather_station=84&radar_station=XSS&clear_sky_clock_station_id=PntctonBC'>BC Penticton
     <option value='graphic%2Esky?country=cda&prov=BC&location_name=Prince+George&latitude=53%2E92&longitude=-122%2E75&weather_station=79&radar_station=XPG&clear_sky_clock_station_id=Prince_George'>BC Prince George
     <option value='graphic%2Esky?country=cda&prov=BC&location_name=Terrace&latitude=54%2E53&longitude=-128%2E61&weather_station=80&radar_station=&clear_sky_clock_station_id=TrrcBC'>BC Terrace
     <option value='graphic%2Esky?country=cda&prov=BC&location_name=Tofino&latitude=49%2E14&longitude=-125%2E90&weather_station=17&radar_station=XSI&clear_sky_clock_station_id=TofinoBC'>BC Tofino
     <option value='graphic%2Esky?country=cda&prov=BC&location_name=Vancouver&latitude=49%2E30&longitude=-123%2E14&weather_station=74&radar_station=WUJ&clear_sky_clock_station_id=Vancouver'>BC Vancouver
     <option value='graphic%2Esky?country=cda&prov=BC&location_name=Victoria&latitude=48%2E41&longitude=-123%2E37&weather_station=85&radar_station=XSI&clear_sky_clock_station_id=Victoria'>BC Victoria
     <option value='graphic%2Esky?country=cda&prov=BC&location_name=Whistler&latitude=50%2E12&longitude=-122%2E96&weather_station=86&radar_station=WUJ&clear_sky_clock_station_id=WhstlrBC'>BC Whistler
     <option value='graphic%2Esky?country=cda&prov=MB&location_name=Brandon&latitude=49%2E85&longitude=-99%2E95&weather_station=52&radar_station=XFW&clear_sky_clock_station_id=BrndnMB'>MB Brandon
     <option value='graphic%2Esky?country=cda&prov=MB&location_name=Churchill&latitude=58%2E77&longitude=-94%2E16&weather_station=42&radar_station=&clear_sky_clock_station_id=ChrchllMB'>MB Churchill
     <option value='graphic%2Esky?country=cda&prov=MB&location_name=Dauphin&latitude=51%2E15&longitude=-100%2E06&weather_station=58&radar_station=XFW&clear_sky_clock_station_id=MRossObs'>MB Dauphin
     <option value='graphic%2Esky?country=cda&prov=MB&location_name=The+Pas&latitude=53%2E82&longitude=-101%2E25&weather_station=30&radar_station=&clear_sky_clock_station_id=TPasMB'>MB The Pas
     <option value='graphic%2Esky?country=cda&prov=MB&location_name=Thompson&latitude=55%2E75&longitude=-97%2E85&weather_station=34&radar_station=&clear_sky_clock_station_id=ThmpsnMB'>MB Thompson
     <option value='graphic%2Esky?country=cda&prov=MB&location_name=Winnipeg&latitude=49%2E91&longitude=-97%2E14&weather_station=38&radar_station=XWL&clear_sky_clock_station_id=Winnipeg'>MB Winnipeg
     <option value='graphic%2Esky?country=cda&prov=NB&location_name=Bathurst&latitude=47%2E61&longitude=-65%2E64&weather_station=28&radar_station=XNC&clear_sky_clock_station_id=BathurstNB'>NB Bathurst
     <option value='graphic%2Esky?country=cda&prov=NB&location_name=Edmunston&latitude=47%2E37&longitude=-68%2E32&weather_station=32&radar_station=XAM&clear_sky_clock_station_id=RachelObNB'>NB Edmunston
     <option value='graphic%2Esky?country=cda&prov=NB&location_name=Fredericton&latitude=45%2E96&longitude=-66%2E65&weather_station=29&radar_station=XNC&clear_sky_clock_station_id=FrederictonNB'>NB Fredericton
     <option value='graphic%2Esky?country=cda&prov=NB&location_name=Fundy+National+Park&latitude=45%2E55&longitude=-65%2E02&weather_station=5&radar_station=XNC&clear_sky_clock_station_id=FndyNPNB'>NB Fundy National Park
     <option value='graphic%2Esky?country=cda&prov=NB&location_name=Kouchibouguac&latitude=46%2E83&longitude=-64%2E93&weather_station=9&radar_station=XNC&clear_sky_clock_station_id=KchbgNPNB'>NB Kouchibouguac
     <option value='graphic%2Esky?country=cda&prov=NB&location_name=Moncton&latitude=46%2E08&longitude=-64%2E78&weather_station=36&radar_station=XNC&clear_sky_clock_station_id=Moncton'>NB Moncton
     <option value='graphic%2Esky?country=cda&prov=NB&location_name=Mount+Carleton+Provincial+Park&latitude=47%2E43&longitude=-66%2E91&weather_station=10&radar_station=XAM&clear_sky_clock_station_id=AASPNB'>NB Mount Carleton Provincial Park
     <option value='graphic%2Esky?country=cda&prov=NB&location_name=St%2E+John&latitude=45%2E27&longitude=-66%2E07&weather_station=23&radar_station=XNC&clear_sky_clock_station_id=SntJhnNMB'>NB St. John
     <option value='graphic%2Esky?country=cda&prov=NL&location_name=Badger&latitude=49%2E50&longitude=-56%2E07&weather_station=34&radar_station=XME&clear_sky_clock_station_id=SprngdNL'>NL Badger
     <option value='graphic%2Esky?country=cda&prov=NL&location_name=Cornerbrook&latitude=48%2E95&longitude=-57%2E95&weather_station=41&radar_station=XME&clear_sky_clock_station_id=CrnrBrkNF'>NL Cornerbrook
     <option value='graphic%2Esky?country=cda&prov=NL&location_name=Happy+Valley-Goose+Bay&latitude=53%2E30&longitude=-60%2E35&weather_station=23&radar_station=&clear_sky_clock_station_id=GssByNFLD'>NL Happy Valley-Goose Bay
     <option value='graphic%2Esky?country=cda&prov=NL&location_name=St%2E+John%27s&latitude=47%2E56&longitude=-52%2E69&weather_station=24&radar_station=WTP&clear_sky_clock_station_id=StJohns'>NL St. John's
     <option value='graphic%2Esky?country=cda&prov=NS&location_name=Halifax&latitude=44%2E66&longitude=-63%2E59&weather_station=19&radar_station=XGO&clear_sky_clock_station_id=Halifax'>NS Halifax
     <option value='graphic%2Esky?country=cda&prov=NS&location_name=Kejimkujik&latitude=44%2E42&longitude=-65%2E26&weather_station=42&radar_station=XGO&clear_sky_clock_station_id=KjmkjkNS'>NS Kejimkujik
     <option value='graphic%2Esky?country=cda&prov=NS&location_name=Sydney&latitude=46%2E13&longitude=-60%2E19&weather_station=31&radar_station=XMB&clear_sky_clock_station_id=SydnyNS'>NS Sydney
     <option value='graphic%2Esky?country=cda&prov=NS&location_name=Truro&latitude=45%2E36&longitude=-63%2E29&weather_station=25&radar_station=XGO&clear_sky_clock_station_id=TruroNS'>NS Truro
     <option value='graphic%2Esky?country=cda&prov=NS&location_name=Yarmouth&latitude=43%2E83&longitude=-66%2E12&weather_station=29&radar_station=XGO&clear_sky_clock_station_id=ArgyleNS'>NS Yarmouth
     <option value='graphic%2Esky?country=cda&prov=NT&location_name=Yellowknife&latitude=62%2E45&longitude=-114%2E37&weather_station=24&radar_station=&clear_sky_clock_station_id=YllwknfNWT'>NT Yellowknife
     <option value='graphic%2Esky?country=cda&prov=NU&location_name=Iqaluit&latitude=63%2E75&longitude=-68%2E52&weather_station=21&radar_station=&clear_sky_clock_station_id=IqltNvT'>NU Iqaluit
     <option value='graphic%2Esky?country=cda&prov=ON&location_name=Algonquin+Park&latitude=45%2E58&longitude=-78%2E41&weather_station=29&radar_station=WBI&clear_sky_clock_station_id=RckLkAqON'>ON Algonquin Park
     <option value='graphic%2Esky?country=cda&prov=ON&location_name=Barrie&latitude=44%2E39&longitude=-79%2E69&weather_station=151&radar_station=WKR&clear_sky_clock_station_id=Barrie'>ON Barrie
     <option value='graphic%2Esky?country=cda&prov=ON&location_name=Hamilton&latitude=43%2E25&longitude=-79%2E87&weather_station=77&radar_station=WKR&clear_sky_clock_station_id=Hamilton'>ON Hamilton
     <option value='graphic%2Esky?country=cda&prov=ON&location_name=Kenora&latitude=49%2E77&longitude=-94%2E50&weather_station=96&radar_station=XDR&clear_sky_clock_station_id=KenoraON'>ON Kenora
     <option value='graphic%2Esky?country=cda&prov=ON&location_name=Kingston&latitude=44%2E23&longitude=-76%2E54&weather_station=69&radar_station=XFT&clear_sky_clock_station_id=Kingston'>ON Kingston
     <option value='graphic%2Esky?country=cda&prov=ON&location_name=London&latitude=42%2E99&longitude=-81%2E25&weather_station=137&radar_station=WSO&clear_sky_clock_station_id=London'>ON London
     <option value='graphic%2Esky?country=cda&prov=ON&location_name=Moosonee&latitude=51%2E28&longitude=-80%2E65&weather_station=113&radar_station=XTI&clear_sky_clock_station_id=MsneeON'>ON Moosonee
     <option value='graphic%2Esky?country=cda&prov=ON&location_name=North+Frontenac&latitude=44%2E92&longitude=-76%2E94&weather_station=106&radar_station=XFT&clear_sky_clock_station_id=PlvnAdON'>ON North Frontenac
     <option value='graphic%2Esky?country=cda&prov=ON&location_name=Ottawa&latitude=45%2E40&longitude=-75%2E66&weather_station=118&radar_station=XFT&clear_sky_clock_station_id=FLO'>ON Ottawa
     <option value='graphic%2Esky?country=cda&prov=ON&location_name=Parry+Sound&latitude=45%2E35&longitude=-80%2E04&weather_station=103&radar_station=WBI&clear_sky_clock_station_id=PrrySndON'>ON Parry Sound
     <option value='graphic%2Esky?country=cda&prov=ON&location_name=Peterborough&latitude=44%2E31&longitude=-78%2E32&weather_station=121&radar_station=WKR&clear_sky_clock_station_id=Peterborough'>ON Peterborough
     <option value='graphic%2Esky?country=cda&prov=ON&location_name=Sarnia&latitude=42%2E98&longitude=-82%2E41&weather_station=147&radar_station=WSO&clear_sky_clock_station_id=Sarnia'>ON Sarnia
     <option value='graphic%2Esky?country=cda&prov=ON&location_name=Sault+Ste%2E+Marie&latitude=46%2E53&longitude=-84%2E36&weather_station=162&radar_station=WGJ&clear_sky_clock_station_id=SaultON'>ON Sault Ste. Marie
     <option value='graphic%2Esky?country=cda&prov=ON&location_name=Sudbury&latitude=46%2E52&longitude=-80%2E96&weather_station=40&radar_station=WBI&clear_sky_clock_station_id=Sudbury'>ON Sudbury
     <option value='graphic%2Esky?country=cda&prov=ON&location_name=Thunder+Bay&latitude=48%2E38&longitude=-89%2E24&weather_station=100&radar_station=XNI&clear_sky_clock_station_id=ThunderBay'>ON Thunder Bay
     <option value='graphic%2Esky?country=cda&prov=ON&location_name=Timmins&latitude=48%2E48&longitude=-81%2E34&weather_station=127&radar_station=XTI&clear_sky_clock_station_id=Timmins'>ON Timmins
     <option value='graphic%2Esky?country=cda&prov=ON&location_name=Toronto&latitude=43%2E65&longitude=-79%2E38&weather_station=143&radar_station=WKR&clear_sky_clock_station_id=Toronto'>ON Toronto
     <option value='graphic%2Esky?country=cda&prov=ON&location_name=Windsor&latitude=42%2E31&longitude=-83%2E04&weather_station=94&radar_station=WSO&clear_sky_clock_station_id=Windsor'>ON Windsor
     <option value='graphic%2Esky?country=cda&prov=PE&location_name=Charlottetown&latitude=46%2E26&longitude=-63%2E14&weather_station=5&radar_station=XGO&clear_sky_clock_station_id=Charlottetown'>PE Charlottetown
     <option value='graphic%2Esky?country=cda&prov=PE&location_name=Summerside&latitude=46%2E39&longitude=-63%2E79&weather_station=3&radar_station=XGO&clear_sky_clock_station_id=SlmnPkPEI'>PE Summerside
     <option value='graphic%2Esky?country=cda&prov=QC&location_name=Gaspé&latitude=48%2E83&longitude=-64%2E51&weather_station=101&radar_station=XAM&clear_sky_clock_station_id=GaspePQ'>QC Gaspé
     <option value='graphic%2Esky?country=cda&prov=QC&location_name=Iles-de-la-Madeleine&latitude=47%2E40&longitude=-61%2E84&weather_station=103&radar_station=XMB&clear_sky_clock_station_id=BtMntTQC'>QC Iles-de-la-Madeleine
     <option value='graphic%2Esky?country=cda&prov=QC&location_name=La+Vérendrye&latitude=46%2E98&longitude=-76%2E48&weather_station=30&radar_station=XLA&clear_sky_clock_station_id=RsrvFnqPQ'>QC La Vérendrye
     <option value='graphic%2Esky?country=cda&prov=QC&location_name=Mont+Mégantic&latitude=45%2E46&longitude=-71%2E15&weather_station=136&radar_station=WVY&clear_sky_clock_station_id=Observatoires_du_Mont_Megantic'>QC Mont Mégantic
     <option value='graphic%2Esky?country=cda&prov=QC&location_name=Montréal&latitude=45%2E50&longitude=-73%2E58&weather_station=147&radar_station=WMN&clear_sky_clock_station_id=BvObQC'>QC Montréal
     <option value='graphic%2Esky?country=cda&prov=QC&location_name=Québec&latitude=46%2E82&longitude=-71%2E22&weather_station=133&radar_station=WVY&clear_sky_clock_station_id=Quebec'>QC Québec
     <option value='graphic%2Esky?country=cda&prov=QC&location_name=Rimouski&latitude=48%2E45&longitude=-68%2E53&weather_station=138&radar_station=XAM&clear_sky_clock_station_id=RmskiQC'>QC Rimouski
     <option value='graphic%2Esky?country=cda&prov=QC&location_name=Rivière-du-Loup&latitude=47%2E83&longitude=-69%2E54&weather_station=108&radar_station=WMB&clear_sky_clock_station_id=PtQllObPQ'>QC Rivière-du-Loup
     <option value='graphic%2Esky?country=cda&prov=QC&location_name=Saguenay&latitude=48%2E43&longitude=-71%2E07&weather_station=166&radar_station=WMB&clear_sky_clock_station_id=SrsSgnyPQ'>QC Saguenay
     <option value='graphic%2Esky?country=cda&prov=QC&location_name=Sept-Iles&latitude=50%2E22&longitude=-66%2E37&weather_station=141&radar_station=XAM&clear_sky_clock_station_id=SeptIlesPQ'>QC Sept-Iles
     <option value='graphic%2Esky?country=cda&prov=QC&location_name=Sherbrooke&latitude=45%2E40&longitude=-71%2E88&weather_station=136&radar_station=WVY&clear_sky_clock_station_id=BspUObQC'>QC Sherbrooke
     <option value='graphic%2Esky?country=cda&prov=QC&location_name=Trois-Rivières&latitude=46%2E35&longitude=-72%2E59&weather_station=130&radar_station=WVY&clear_sky_clock_station_id=TroisRivPQ'>QC Trois-Rivières
     <option value='graphic%2Esky?country=cda&prov=QC&location_name=Val-d%27or&latitude=48%2E10&longitude=-77%2E80&weather_station=149&radar_station=XLA&clear_sky_clock_station_id=ValdOrPQ'>QC Val-d'or
     <option value='graphic%2Esky?country=cda&prov=SK&location_name=Cypress+Hills+–+West+Block&latitude=49%2E60&longitude=-109%2E92&weather_station=29&radar_station=XBU&clear_sky_clock_station_id=CHDSPWBAB'>SK Cypress Hills – West Block
     <option value='graphic%2Esky?country=cda&prov=SK&location_name=Grasslands+National+Park+–+East&latitude=49%2E07&longitude=-106%2E53&weather_station=28&radar_station=XBU&clear_sky_clock_station_id=GNPEBSK'>SK Grasslands National Park – East
     <option value='graphic%2Esky?country=cda&prov=SK&location_name=Grasslands+National+Park+–+West&latitude=49%2E18&longitude=-107%2E71&weather_station=28&radar_station=XBU&clear_sky_clock_station_id=GrssNPESK'>SK Grasslands National Park – West
     <option value='graphic%2Esky?country=cda&prov=SK&location_name=La+Ronge&latitude=55%2E11&longitude=-105%2E28&weather_station=38&radar_station=&clear_sky_clock_station_id=LRngSK'>SK La Ronge
     <option value='graphic%2Esky?country=cda&prov=SK&location_name=Prince+Albert&latitude=53%2E20&longitude=-105%2E73&weather_station=27&radar_station=XRA&clear_sky_clock_station_id=PrncAlbrtSK'>SK Prince Albert
     <option value='graphic%2Esky?country=cda&prov=SK&location_name=Regina&latitude=50%2E45&longitude=-104%2E62&weather_station=32&radar_station=XBE&clear_sky_clock_station_id=Regina'>SK Regina
     <option value='graphic%2Esky?country=cda&prov=SK&location_name=Saskatoon&latitude=52%2E13&longitude=-106%2E69&weather_station=40&radar_station=XRA&clear_sky_clock_station_id=Saskatoon'>SK Saskatoon
     <option value='graphic%2Esky?country=cda&prov=SK&location_name=Swift+Current&latitude=50%2E29&longitude=-107%2E79&weather_station=41&radar_station=XBU&clear_sky_clock_station_id=SwtCrnSK'>SK Swift Current
     <option value='graphic%2Esky?country=cda&prov=YT&location_name=Whitehorse&latitude=60%2E72&longitude=-135%2E06&weather_station=16&radar_station=&clear_sky_clock_station_id=WhthrsYK'>YT Whitehorse
     <option value=''><s:txt>---U.S.A.---</s:txt>
     <option value='graphic%2Esky?country=us&location_name=Anchorage&latitude=61%2E20&longitude=-149%2E88&radar_station=AHG&clear_sky_clock_station_id=AnchorAK'>AK Anchorage
     <option value='graphic%2Esky?country=us&location_name=Phoenix&latitude=33%2E43&longitude=-112%2E09&radar_station=IWA&clear_sky_clock_station_id=Phoenix'>AZ Phoenix
     <option value='graphic%2Esky?country=us&location_name=Los+Angeles&latitude=34%2E03&longitude=-118%2E31&radar_station=SOX&clear_sky_clock_station_id=LAXCA'>CA Los Angeles
     <option value='graphic%2Esky?country=us&location_name=San+Diego&latitude=32%2E74&longitude=-117%2E21&radar_station=NKX&clear_sky_clock_station_id=SanDiego'>CA San Diego
     <option value='graphic%2Esky?country=us&location_name=San+Francisco&latitude=37%2E77&longitude=-122%2E43&radar_station=MUX&clear_sky_clock_station_id=SanFranCA'>CA San Francisco
     <option value='graphic%2Esky?country=us&location_name=Denver&latitude=39%2E73&longitude=-104%2E97&radar_station=FTG&clear_sky_clock_station_id=DenverCO'>CO Denver
     <option value='graphic%2Esky?country=us&location_name=Washington&latitude=38%2E89&longitude=-77%2E03&radar_station=LWX&clear_sky_clock_station_id=WashingtonDC'>DC Washington
     <option value='graphic%2Esky?country=us&location_name=Jacksonville&latitude=30%2E32&longitude=-81%2E67&radar_station=JAX&clear_sky_clock_station_id=JacksonFL'>FL Jacksonville
     <option value='graphic%2Esky?country=us&location_name=Miami&latitude=25%2E76&longitude=-80%2E20&radar_station=AMX&clear_sky_clock_station_id=MiamiFL'>FL Miami
     <option value='graphic%2Esky?country=us&location_name=Atlanta&latitude=33%2E75&longitude=-84%2E39&radar_station=FFC&clear_sky_clock_station_id=AtlantaGA'>GA Atlanta
     <option value='graphic%2Esky?country=us&location_name=Des+Moines&latitude=41%2E58&longitude=-93%2E60&radar_station=DMX&clear_sky_clock_station_id=DeMoinesIA'>IA Des Moines
     <option value='graphic%2Esky?country=us&location_name=Boise&latitude=43%2E60&longitude=-116%2E24&radar_station=CBX&clear_sky_clock_station_id=BoiseID'>ID Boise
     <option value='graphic%2Esky?country=us&location_name=Chicago&latitude=41%2E88&longitude=-87%2E63&radar_station=LOT&clear_sky_clock_station_id=Chicago'>IL Chicago
     <option value='graphic%2Esky?country=us&location_name=Indianapolis&latitude=39%2E76&longitude=-86%2E16&radar_station=IND&clear_sky_clock_station_id=Indianapolis'>IN Indianapolis
     <option value='graphic%2Esky?country=us&location_name=Kansas+City&latitude=39%2E11&longitude=-94%2E63&radar_station=EAX&clear_sky_clock_station_id=KansasKA'>KS Kansas City
     <option value='graphic%2Esky?country=us&location_name=New+Orleans&latitude=29%2E98&longitude=-90%2E09&radar_station=LIX&clear_sky_clock_station_id=NewOrlLA'>LA New Orleans
     <option value='graphic%2Esky?country=us&location_name=Boston&latitude=42%2E36&longitude=-71%2E06&radar_station=BOX&clear_sky_clock_station_id=Boston'>MA Boston
     <option value='graphic%2Esky?country=us&location_name=Bangor&latitude=44%2E80&longitude=-68%2E78&radar_station=CBW&clear_sky_clock_station_id=BngrME'>ME Bangor
     <option value='graphic%2Esky?country=us&location_name=Detroit&latitude=42%2E35&longitude=-83%2E07&radar_station=DTX&clear_sky_clock_station_id=Detroit'>MI Detroit
     <option value='graphic%2Esky?country=us&location_name=Minneapolis&latitude=44%2E98&longitude=-93%2E27&radar_station=MPX&clear_sky_clock_station_id=Minneapolis'>MN Minneapolis
     <option value='graphic%2Esky?country=us&location_name=St%2E+Louis&latitude=38%2E62&longitude=-90%2E23&radar_station=LSX&clear_sky_clock_station_id=StLouisMO'>MO St. Louis
     <option value='graphic%2Esky?country=us&location_name=Great+Falls&latitude=47%2E51&longitude=-111%2E29&radar_station=TFX&clear_sky_clock_station_id=GrFallsMT'>MT Great Falls
     <option value='graphic%2Esky?country=us&location_name=Charlotte&latitude=35%2E21&longitude=-80%2E82&radar_station=GSP&clear_sky_clock_station_id=CharlotteNC'>NC Charlotte
     <option value='graphic%2Esky?country=us&location_name=Bismark&latitude=46%2E81&longitude=-100%2E77&radar_station=BIS&clear_sky_clock_station_id=BismarkND'>ND Bismark
     <option value='graphic%2Esky?country=us&location_name=Omaha&latitude=41%2E25&longitude=-96%2E01&radar_station=OAX&clear_sky_clock_station_id=OmahaNE'>NE Omaha
     <option value='graphic%2Esky?country=us&location_name=Albuquerque&latitude=35%2E10&longitude=-106%2E61&radar_station=ABX&clear_sky_clock_station_id=AlbuqNM'>NM Albuquerque
     <option value='graphic%2Esky?country=us&location_name=Socorro&latitude=34%2E06&longitude=-106%2E91&radar_station=ABX&clear_sky_clock_station_id=EtscornNM'>NM Socorro
     <option value='graphic%2Esky?country=us&location_name=Buffalo&latitude=42%2E89&longitude=-78%2E87&radar_station=BUF&clear_sky_clock_station_id=BuffaloNY'>NY Buffalo
     <option value='graphic%2Esky?country=us&location_name=New+York&latitude=40%2E76&longitude=-73%2E98&radar_station=OKX&clear_sky_clock_station_id=NYCNY'>NY New York
     <option value='graphic%2Esky?country=us&location_name=Cleveland&latitude=41%2E49&longitude=-81%2E70&radar_station=CLE&clear_sky_clock_station_id=Cleveland'>OH Cleveland
     <option value='graphic%2Esky?country=us&location_name=Columbus&latitude=39%2E95&longitude=-82%2E97&radar_station=ILN&clear_sky_clock_station_id=Columbus'>OH Columbus
     <option value='graphic%2Esky?country=us&location_name=Oklahoma+City&latitude=35%2E44&longitude=-97%2E53&radar_station=TLX&clear_sky_clock_station_id=OklahomaOK'>OK Oklahoma City
     <option value='graphic%2Esky?country=us&location_name=Portland&latitude=45%2E52&longitude=-122%2E64&radar_station=RTX&clear_sky_clock_station_id=PortOR'>OR Portland
     <option value='graphic%2Esky?country=us&location_name=Philadelphia&latitude=39%2E94&longitude=-75%2E16&radar_station=DIX&clear_sky_clock_station_id=Philadelphia'>PA Philadelphia
     <option value='graphic%2Esky?country=us&location_name=Charleston&latitude=32%2E77&longitude=-79%2E95&radar_station=CLX&clear_sky_clock_station_id=ChrlstnSC'>SC Charleston
     <option value='graphic%2Esky?country=us&location_name=Sioux+Falls&latitude=43%2E54&longitude=-96%2E73&radar_station=FSD&clear_sky_clock_station_id=SxFllsSD'>SD Sioux Falls
     <option value='graphic%2Esky?country=us&location_name=Memphis&latitude=35%2E13&longitude=-90%2E03&radar_station=NQA&clear_sky_clock_station_id=MemphisTN'>TN Memphis
     <option value='graphic%2Esky?country=us&location_name=Nashville&latitude=36%2E17&longitude=-86%2E76&radar_station=OHX&clear_sky_clock_station_id=NashTE'>TN Nashville
     <option value='graphic%2Esky?country=us&location_name=Austin&latitude=30%2E26&longitude=-97%2E74&radar_station=GRK&clear_sky_clock_station_id=AustinTX'>TX Austin
     <option value='graphic%2Esky?country=us&location_name=Dallas&latitude=32%2E76&longitude=-96%2E80&radar_station=FWS&clear_sky_clock_station_id=Dallas'>TX Dallas
     <option value='graphic%2Esky?country=us&location_name=El+Paso&latitude=31%2E77&longitude=-106%2E49&radar_station=EPZ&clear_sky_clock_station_id=ElPTX'>TX El Paso
     <option value='graphic%2Esky?country=us&location_name=Houston&latitude=29%2E75&longitude=-95%2E38&radar_station=HGX&clear_sky_clock_station_id=Houston'>TX Houston
     <option value='graphic%2Esky?country=us&location_name=San+Antonio&latitude=29%2E40&longitude=-98%2E50&radar_station=EWX&clear_sky_clock_station_id=SanAnTX'>TX San Antonio
     <option value='graphic%2Esky?country=us&location_name=Salt+Lake+City&latitude=40%2E76&longitude=-111%2E94&radar_station=MTX&clear_sky_clock_station_id=SaltLakeUT'>UT Salt Lake City
     <option value='graphic%2Esky?country=us&location_name=Norfolk&latitude=36%2E86&longitude=-76%2E29&radar_station=AKQ&clear_sky_clock_station_id=TidewaterVA'>VA Norfolk
     <option value='graphic%2Esky?country=us&location_name=Seattle&latitude=47%2E60&longitude=-122%2E32&radar_station=ATX&clear_sky_clock_station_id=Seattle'>WA Seattle
     <option value='graphic%2Esky?country=us&location_name=Spokane&latitude=47%2E67&longitude=-117%2E42&radar_station=OTX&clear_sky_clock_station_id=SpknWA'>WA Spokane
     <option value=''><s:txt>---United Kingdom---</s:txt>
     <option value='graphic%2Esky?country=uk&location_name=Aberdeen&latitude=57%2E14&longitude=-2%2E11&weather_station=3091&'>Aberdeen
     <option value='graphic%2Esky?country=uk&location_name=Belfast&latitude=54%2E59&longitude=-5%2E93&weather_station=3917&'>Belfast
     <option value='graphic%2Esky?country=uk&location_name=Cardiff&latitude=51%2E48&longitude=-3%2E18&weather_station=3716&'>Cardiff
     <option value='graphic%2Esky?country=uk&location_name=Derry&latitude=55%2E03&longitude=-7%2E31&weather_station=3907&'>Derry
     <option value='graphic%2Esky?country=uk&location_name=Dumfries&latitude=55%2E07&longitude=-3%2E62&weather_station=3153&'>Dumfries
     <option value='graphic%2Esky?country=uk&location_name=Edinburgh&latitude=55%2E96&longitude=-3%2E19&weather_station=3166&'>Edinburgh
     <option value='graphic%2Esky?country=uk&location_name=Glasgow&latitude=55%2E86&longitude=-4%2E26&weather_station=3134&'>Glasgow
     <option value='graphic%2Esky?country=uk&location_name=Inverness&latitude=57%2E47&longitude=-4%2E23&weather_station=3063&'>Inverness
     <option value='graphic%2Esky?country=uk&location_name=Leeds&latitude=53%2E79&longitude=-1%2E54&weather_station=3344&'>Leeds
     <option value='graphic%2Esky?country=uk&location_name=Lerwick&latitude=60%2E15&longitude=-1%2E15&weather_station=3005&'>Lerwick
     <option value='graphic%2Esky?country=uk&location_name=Liverpool&latitude=53%2E40&longitude=-3%2E00&weather_station=3316&'>Liverpool
     <option value='graphic%2Esky?country=uk&location_name=London&latitude=51%2E50&longitude=-0%2E13&weather_station=3772&'>London
     <option value='graphic%2Esky?country=uk&location_name=Newcastle&latitude=54%2E98&longitude=-1%2E61&weather_station=3238&'>Newcastle
     <option value='graphic%2Esky?country=uk&location_name=Oxford&latitude=51%2E75&longitude=-1%2E26&weather_station=3469&'>Oxford
     <option value='graphic%2Esky?country=uk&location_name=Peterborough&latitude=52%2E57&longitude=-0%2E24&weather_station=3462&'>Peterborough
     <option value='graphic%2Esky?country=uk&location_name=Plymouth&latitude=50%2E37&longitude=-4%2E14&weather_station=3827&'>Plymouth
     <option value='graphic%2Esky?country=uk&location_name=Stornoway&latitude=58%2E21&longitude=-6%2E38&weather_station=3026&'>Stornoway
   </select>
  </form>
  
  <a href='form.sky' id='link_back' class='no-print'><s:txt>Form</s:txt></a>
 </span>
 

 <P>
 <div>  
  <table id='weather' class='report weather'></table>
  <canvas id='clouds' class='no-print cloud' style='margin:0.25em;' title='<s:txt>Satellite image. Click for a larger image.</s:txt>'></canvas>
  <img id='radar' class='no-print radar'>
  <canvas id='radar_us'  class='no-print radar' style='margin:0.25em;'></canvas> 
  <div id='radar_uk'  class="no-print radar" style='margin:0.25em;float:right;'></div>  
  <!-- <div id='radar_uk'  class="no-print radar" style='float:right;margin:0.25em;'></div> -->  
 </div>
 
<P><canvas id='ecliptic' class='sky no-print' title='<s:txt>Ecliptic, showing the planets and Moon</s:txt>'></canvas>
<P><canvas id='observation_window' class='no-print' title='<s:txt>Observation window: twilight and brightness of the Moon</s:txt>'></canvas>
<P><canvas id='planisphere' class='no-print'></canvas>
 
<P style='text-align:center;'><img id='clear_sky_clock' src='' class='no-print' style='border-radius: 10px;'>
<P>
<canvas id='libration' class='sky no-print' title='<s:txt>Lunar libration, phase, and change in size from mean</s:txt>'></canvas>
<canvas id='galilean_satellites'  class='sky no-print'  style='margin-right:0.25em;' title='<s:txt>Galilean satellites of Jupiter</s:txt>'></canvas>
 
 <table id='sun_and_moon_and_planets_table' class='report' style="text-align:right;">
  <caption><s:txt>Sun, Moon, and Planets</s:txt> 
   <tr style="text-align:center;">
    <th><s:txt>Name</s:txt>
    <th><s:txt>Const</s:txt>
    <th><s:txt>Elong</s:txt>
    <th><s:txt>Mag</s:txt>
    <th><s:txt>Size</s:txt>
    <th><s:txt>Illum</s:txt>
    <th><s:txt>Alt</s:txt>
    <th><s:txt>Az</s:txt>
    <th>α
    <th>δ
    <th><s:txt>Rise</s:txt>
    <th><s:txt>Transit</s:txt>
    <th><s:txt>Set</s:txt>  
 </table>
 
 <table id='diary_table' class='report'>
   <caption><s:txt>Sky Diary</s:txt>
   <tr>
    <th><s:txt>When</s:txt> 
    <th><s:txt>Description</s:txt>
 </table>
 
 <table id='minor_planet_table' class='report'>
   <caption><s:txt>Minor Planets</s:txt> 
   <tr>
    <th><s:txt>Name</s:txt>
    <th><s:txt>Elong</s:txt>
    <th><s:txt>Mag</s:txt>
    <th><s:txt>Alt</s:txt>
    <th><s:txt>Az</s:txt>
    <th>α
    <th>δ
 </table>
 
 <table id='comet_table' class='report'>
  <caption><a href="https://www.ast.cam.ac.uk/~jds/"><s:txt>Bright Comets</s:txt></a> 
  <tr>
    <th><s:txt>Name</s:txt>
    <th><s:txt>Elong</s:txt>
    <th><s:txt>Mag</s:txt>
    <th><s:txt>Alt</s:txt>
    <th><s:txt>Az</s:txt>
    <th>α
    <th>δ
    <th><s:txt>Trend</s:txt>
    <th><s:txt>Visible</s:txt>
 </table>
 
 <table id='jupiter_satellite_phenomena' class='report' title='<s:txt>Times are approximate</s:txt>'>
  <caption><s:txt>Jupiter Satellite Phenomena</s:txt>
   <tr>
    <th><s:txt>When</s:txt>
    <th><s:txt>Satellite</s:txt>
    <th><s:txt>Event</s:txt>
 </table>
 
 <table id='aurora_table' class='report'>
   <caption><a href="http://www.swpc.noaa.gov/products/aurora-30-minute-forecast"><s:txt>Auroral Activity Level</s:txt></a> 
 </table>
 
 <table id='occultation_table' class='report' title='<s:txt>Occultation times are approximate</s:txt>'>
  <caption><a href="http://www.lunar-occultations.com/iota/iotandx.htm"><s:txt>Occultations</s:txt></a>
  <tr>
   <th><s:txt>When</s:txt>
   <th><s:txt>ZC</s:txt>
   <th><s:txt>Mag</s:txt>
   <th><s:txt>Ph</s:txt>
   <th><s:txt>El</s:txt>
   <th><s:txt>PA</s:txt>
   <th><s:txt>Alt</s:txt>
 </table>
 
 <table id='meteor_shower_table' class='report'>
   <caption><s:txt>Meteor Showers</s:txt> 
   <tr>
    <th><s:txt>Name</s:txt> 
    <th><s:txt>Desig</s:txt>
    <th><s:txt>Peak</s:txt>
    <th><s:txt>Km/s</s:txt>
    <th><s:txt>Alt</s:txt>
    <th><s:txt>Az</s:txt>
    <th>r
    <th><s:txt>ZHR</s:txt>
    <th><s:txt>ZHR</s:txt> x
 </table>

 <br><table id='messiers' class='report'>
  <caption><s:txt>Messier Objects</s:txt>
  <tr>
    <th><s:txt>Name</s:txt>
    <th><s:txt>Const</s:txt>
    <th><s:txt>Chart</s:txt>
    <th><s:txt>Mag</s:txt>
    <th><s:txt>Alt</s:txt>
    <th><s:txt>Az</s:txt>
    <th>α
    <th>δ
    <th><s:txt>Type</s:txt>
    <th><s:txt>Comment</s:txt>
 </table>
 
 <br><table id='caldwells' class='report'>
  <caption><s:txt>Caldwell Objects</s:txt>
  <tr>
    <th><s:txt>Name</s:txt>
    <th><s:txt>Const</s:txt>
    <th><s:txt>Chart</s:txt>
    <th><s:txt>Mag</s:txt>
    <th><s:txt>Alt</s:txt>
    <th><s:txt>Az</s:txt>
    <th>α
    <th>δ
    <th><s:txt>Type</s:txt>
    <th><s:txt>Comment</s:txt>
 </table>
 
 <tags:analytics/>
 
</body>
</html>