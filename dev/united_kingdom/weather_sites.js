/*
  List of all sites in the UK that have hourly observations. 
  I use this list as well for the forecast. The forecast is actually available for 
  a much longer list of sites (~5000), but for my purposes that seems to be unnecessary.
  Source:
    http://datapoint.metoffice.gov.uk/public/data/val/wxobs/all/json/sitelist?key=c8258e8b-bf20-45ed-8e95-ca0d76682419
  Docs:
    http://www.metoffice.gov.uk/datapoint/support/documentation/uk-observations-site-list-detailed-documentation
  The pertinent data is:
    id, latitude, longitude, name
    object.Locations.Location[i].id etc
 */
var uk_weather_sites = function(){

  var sites =
    	
	{"Locations":{
	 "Location":[
	 {"elevation":"7.0","id":"3066","latitude":"57.6494","longitude":"-3.5606","name":"Kinloss","region":"gr","unitaryAuthArea":"Moray"}
	,{"elevation":"6.0","id":"3068","latitude":"57.712","longitude":"-3.322","obsSource":"LNDSYN","name":"Lossiemouth","region":"gr","unitaryAuthArea":"Moray"}
	,{"elevation":"36.0","id":"3075","latitude":"58.454","longitude":"-3.089","obsSource":"LNDSYN","name":"Wick John O Groats Airport","region":"he","unitaryAuthArea":"Highland"}
	,{"elevation":"15.0","id":"3002","latitude":"60.749","longitude":"-0.854","name":"Baltasound","region":"os","unitaryAuthArea":"Shetland Islands"}
	,{"elevation":"82.0","id":"3005","latitude":"60.139","longitude":"-1.183","obsSource":"LNDSYN","name":"Lerwick (S. Screen)","region":"os","unitaryAuthArea":"Shetland Islands"}
	,{"elevation":"57.0","id":"3008","latitude":"59.527","longitude":"-1.628","name":"Fair Isle","region":"os","unitaryAuthArea":"Shetland Islands"}
	,{"elevation":"22.0","id":"3014","latitude":"60.111","longitude":"-2.063","obsSource":"LNDSYN","name":"Foula","region":"os","unitaryAuthArea":"Shetland Islands"}
	,{"elevation":"11.0","id":"3034","latitude":"57.859","longitude":"-5.636","name":"Aultbea","region":"he","unitaryAuthArea":"Highland"}
	,{"elevation":"18.0","id":"3037","latitude":"57.257","longitude":"-5.809","obsSource":"LNDSYN","name":"Skye\/Lusa (Samos)","region":"he","unitaryAuthArea":"Highland"}
	,{"elevation":"81.0","id":"3044","latitude":"58.288","longitude":"-4.442","name":"Altnaharra Saws","region":"he","unitaryAuthArea":"Highland"}
	,{"elevation":"237.0","id":"3047","latitude":"56.867","longitude":"-4.708","name":"Tulloch Bridge","region":"he","unitaryAuthArea":"Highland"}
	,{"elevation":"117.0","id":"3796","latitude":"51.133","longitude":"1.348","name":"Langdon Bay","region":"se","unitaryAuthArea":"Kent"}
	,{"elevation":"31.0","id":"3803","latitude":"49.913","longitude":"-6.301","name":"Scilly St Marys","region":"sw","unitaryAuthArea":"Isles of Scilly"}
	,{"elevation":"200.0","id":"3823","latitude":"50.502","longitude":"-4.667","name":"Cardinham","region":"sw","unitaryAuthArea":"Cornwall"}
	,{"elevation":"52.0","id":"3857","latitude":"50.522","longitude":"-2.454","name":"Isle Of Portland","region":"sw","unitaryAuthArea":"Dorset"}
	,{"elevation":"13.0","id":"3876","latitude":"50.836","longitude":"-0.292","name":"Shoreham","region":"se","unitaryAuthArea":"West Sussex"}
	,{"elevation":"84.0","id":"3895","latitude":"49.2079","longitude":"-2.1955","obsSource":"LNDSYN","name":"Jersey","region":"sw","unitaryAuthArea":"Jersey"}
	,{"elevation":"225.0","id":"3911","latitude":"54.721","longitude":"-6.814","name":"Lough Fea Samos","region":"ni","unitaryAuthArea":"County Londonderry"}
	,{"elevation":"156.0","id":"3916","latitude":"55.181","longitude":"-6.153","name":"Ballypatrick Forest","region":"ni","unitaryAuthArea":"County Antrim"}
	,{"elevation":"25.0","id":"3953","latitude":"51.933","longitude":"-10.25","obsSource":"LNDSYN","name":"Valentia Observatory"}
	,{"elevation":"177.0","id":"99081","latitude":"50.767","longitude":"-3.9","name":"North Wyke","region":"sw","unitaryAuthArea":"Devon"}
	,{"elevation":"99.0","id":"3520","latitude":"52.242","longitude":"-2.885","name":"Shobdon Saws","region":"wm","unitaryAuthArea":"Herefordshire"}
	,{"elevation":"75.0","id":"3522","latitude":"52.083","longitude":"-2.8","name":"Hereford","region":"wm","unitaryAuthArea":"Herefordshire"}
	,{"elevation":"96.0","id":"3535","latitude":"52.48","longitude":"-1.689","name":"Coleshill","region":"wm","unitaryAuthArea":"Warwickshire"}
	,{"elevation":"85.0","id":"3560","latitude":"52.225","longitude":"-0.464","name":"Bedford","region":"ee","unitaryAuthArea":"Bedford"}
	,{"elevation":"3.0","id":"3605","latitude":"51.713","longitude":"-4.368","name":"Pembrey Sands Samos","region":"wl","unitaryAuthArea":"Carmarthenshire"}
	,{"elevation":"210.0","id":"3647","latitude":"51.86","longitude":"-1.692","name":"Little Rissington (Esaws)","region":"sw","unitaryAuthArea":"Gloucestershire"}
	,{"elevation":"204.0","id":"3660","latitude":"51.68","longitude":"-0.802","name":"High Wycombe","region":"se","unitaryAuthArea":"Buckinghamshire"}
	,{"elevation":"87.0","id":"3684","latitude":"51.896","longitude":"0.453","name":"Andrewsfield","region":"ee","unitaryAuthArea":"Essex"}
	,{"elevation":"49.0","id":"3716","latitude":"51.405","longitude":"-3.44","name":"St-Athan","region":"wl","unitaryAuthArea":"Vale of Glamorgan"}
	,{"elevation":"145.0","id":"3740","latitude":"51.5031","longitude":"-1.9924","name":"Lyneham","region":"sw","unitaryAuthArea":"Wiltshire"}
	,{"elevation":"72.0","id":"3768","latitude":"51.279","longitude":"-0.772","name":"Farnborough","region":"se","unitaryAuthArea":"Hampshire"}
	,{"elevation":"113.0","id":"3153","latitude":"54.803","longitude":"-4.008","name":"Dundrennan","region":"dg","unitaryAuthArea":"Dumfries and Galloway"}
	,{"elevation":"245.0","id":"3155","latitude":"55.627","longitude":"-3.735","name":"Drumalbin","region":"st","unitaryAuthArea":"South Lanarkshire"}
	,{"elevation":"242.0","id":"3162","latitude":"55.311","longitude":"-3.206","name":"Eskdalemuir","region":"dg","unitaryAuthArea":"Dumfries and Galloway"}
	,{"elevation":"57.0","id":"3166","latitude":"55.928","longitude":"-3.343","name":"Edinburgh\/Gogarbank","region":"dg","unitaryAuthArea":"Edinburgh"}
	,{"elevation":"227.0","id":"3226","latitude":"54.572","longitude":"-2.413","name":"Warcop","region":"nw","unitaryAuthArea":"Cumbria"}
	,{"elevation":"211.0","id":"3230","latitude":"55.285","longitude":"-2.279","name":"Redesdale Camp (Samos)","nationalPark":"Northumberland National Park","region":"ne","unitaryAuthArea":"Northumberland"}
	,{"elevation":"146.0","id":"3238","latitude":"55.02","longitude":"-1.88","name":"Albemarle","region":"ne","unitaryAuthArea":"Northumberland"}
	,{"elevation":"32.0","id":"3257","latitude":"54.296","longitude":"-1.53","name":"Leeming","region":"yh","unitaryAuthArea":"North Yorkshire"}
	,{"elevation":"33.0","id":"3261","latitude":"54.134","longitude":"-1.414","name":"Dishforth Airfield","region":"yh","unitaryAuthArea":"North Yorkshire"}
	,{"elevation":"14.0","id":"3266","latitude":"54.045","longitude":"-1.25","name":"Linton On Ouse","region":"yh","unitaryAuthArea":"North Yorkshire"}
	,{"elevation":"73.0","id":"3462","latitude":"52.6111","longitude":"-0.459","name":"Wittering","region":"ee","unitaryAuthArea":"Peterborough"}
	,{"elevation":"21.0","id":"3482","latitude":"52.651","longitude":"0.569","name":"Marham","region":"ee","unitaryAuthArea":"Norfolk"}
	,{"elevation":"6.0","id":"3907","latitude":"55.1604","longitude":"-6.9492","name":"Magilligan No 2","region":"ni","unitaryAuthArea":"County Londonderry"}
	,{"elevation":"228.0","id":"3063","latitude":"57.206","longitude":"-3.827","name":"Aviemore","nationalPark":"Cairngorms National Park","region":"he","unitaryAuthArea":"Highland"}
	,{"elevation":"9.0","id":"3100","latitude":"56.497","longitude":"-6.887","obsSource":"LNDSYN","name":"Tiree","region":"st","unitaryAuthArea":"Argyll and Bute"}
	,{"elevation":"78.0","id":"3809","latitude":"50.085","longitude":"-5.257","obsSource":"LNDSYN","name":"Culdrose","region":"sw","unitaryAuthArea":"Cornwall"}
	,{"elevation":"252.0","id":"3840","latitude":"50.86","longitude":"-3.239","name":"Dunkeswell Aerodrome","region":"sw","unitaryAuthArea":"Devon"}
	,{"elevation":"10.0","id":"3862","latitude":"50.779","longitude":"-1.835","obsSource":"LNDSYN","name":"Bournemouth Airport","region":"sw","unitaryAuthArea":"Dorset"}
	,{"elevation":"9.0","id":"3976","latitude":"54.233","longitude":"-10.001","obsSource":"LNDSYN","name":"Belmullet"}
	,{"elevation":"348.0","id":"3710","latitude":"51.087","longitude":"-3.608","name":"Liscombe","nationalPark":"Exmoor National Park","region":"sw","unitaryAuthArea":"Somerset"}
	,{"elevation":"132.0","id":"3743","latitude":"51.201","longitude":"-1.805","name":"Larkhill","region":"sw","unitaryAuthArea":"Wiltshire"}
	,{"elevation":"126.0","id":"3746","latitude":"51.161","longitude":"-1.754","name":"Boscombe Down","region":"sw","unitaryAuthArea":"Wiltshire"}
	,{"elevation":"67.0","id":"3769","latitude":"51.15","longitude":"-0.2333","name":"Charlwood","region":"se","unitaryAuthArea":"Surrey"}
	,{"elevation":"10.0","id":"3171","latitude":"56.377","longitude":"-2.862","obsSource":"LNDSYN","name":"Leuchars","region":"ta","unitaryAuthArea":"Fife"}
	,{"elevation":"262.0","id":"3281","latitude":"54.362","longitude":"-0.673","name":"Fylingdales","nationalPark":"North York Moors National Park","region":"yh","unitaryAuthArea":"North Yorkshire"}
	,{"elevation":"10.0","id":"3302","latitude":"53.252","longitude":"-4.537","obsSource":"LNDSYN","name":"Valley","region":"wl","unitaryAuthArea":"Isle of Anglesey"}
	,{"elevation":"262.0","id":"3344","latitude":"53.811","longitude":"-1.865","name":"Bingley Samos","region":"yh","unitaryAuthArea":"West Yorkshire"}
	,{"elevation":"62.0","id":"3379","latitude":"53.031","longitude":"-0.502","name":"Cranwell","region":"em","unitaryAuthArea":"Lincolnshire"}
	,{"elevation":"3.0","id":"3392","latitude":"53.088","longitude":"0.274","name":"Wainfleet","region":"em","unitaryAuthArea":"Lincolnshire"}
	,{"elevation":"72.0","id":"3414","latitude":"52.794","longitude":"-2.663","name":"Shawbury","region":"wm","unitaryAuthArea":"Shropshire"}
	,{"elevation":"1245.0","id":"3065","latitude":"57.116","longitude":"-3.642","name":"Cairn Gorm Summit","nationalPark":"Cairngorms National Park","region":"he","unitaryAuthArea":"Moray"}
	,{"elevation":"933.0","id":"3072","latitude":"56.879","longitude":"-3.42","name":"Cairnwell","nationalPark":"Cairngorms National Park","region":"ta","unitaryAuthArea":"Perth and Kinross"}
	,{"elevation":"65.0","id":"3091","latitude":"57.206","longitude":"-2.202","obsSource":"LNDSYN","name":"Aberdeen Airport","region":"gr","unitaryAuthArea":"Aberdeen"}
	,{"elevation":"17.0","id":"3105","latitude":"55.681","longitude":"-6.256","obsSource":"LNDSYN","name":"Islay Airport","region":"st","unitaryAuthArea":"Argyll and Bute"}
	,{"elevation":"12.0","id":"3010","latitude":"59.083","longitude":"-4.404","obsSource":"LNDSYN","name":"Sule Skerry (Maws)","region":"os","unitaryAuthArea":"Highland"}
	,{"elevation":"20.0","id":"3853","latitude":"51.006","longitude":"-2.64","name":"Yeovilton","region":"sw","unitaryAuthArea":"Somerset"}
	,{"elevation":"20.0","id":"3980","latitude":"55.366","longitude":"-7.333","obsSource":"LNDSYN","name":"Malin Head"}
	,{"elevation":"89.0","id":"99057","latitude":"52.017","longitude":"-0.6","name":"Woburn","region":"ee","unitaryAuthArea":"Central Bedfordshire"}
	,{"elevation":"35.0","id":"3529","latitude":"52.148","longitude":"-2.04","name":"Pershore","region":"wm","unitaryAuthArea":"Worcestershire"}
	,{"elevation":"44.0","id":"3604","latitude":"51.708","longitude":"-5.055","name":"Milford Haven C.B.","region":"wl","unitaryAuthArea":"Pembrokeshire"}
	,{"elevation":"59.0","id":"3628","latitude":"51.521","longitude":"-2.576","name":"Filton","region":"sw","unitaryAuthArea":"South Gloucestershire"}
	,{"elevation":"2.0","id":"3693","latitude":"51.554","longitude":"0.83","obsSource":"LNDSYN","name":"Shoeburyness","region":"ee","unitaryAuthArea":"Essex"}
	,{"elevation":"91.0","id":"3749","latitude":"51.149","longitude":"-1.57","name":"Middle Wallop","region":"se","unitaryAuthArea":"Hampshire"}
	,{"elevation":"118.0","id":"3761","latitude":"51.238","longitude":"-0.944","name":"Odiham","region":"se","unitaryAuthArea":"Hampshire"}
	,{"elevation":"11.0","id":"3132","latitude":"54.859","longitude":"-4.936","name":"West Freugh (Esaws)","region":"dg","unitaryAuthArea":"Dumfries and Galloway"}
	,{"elevation":"27.0","id":"3136","latitude":"55.515","longitude":"-4.585","name":"Prestwick Rnas","region":"st","unitaryAuthArea":"South Ayrshire"}
	,{"elevation":"112.0","id":"3158","latitude":"55.709","longitude":"-2.383","name":"Charterhall","region":"dg","unitaryAuthArea":"Scottish Borders"}
	,{"elevation":"124.0","id":"3210","latitude":"54.5181","longitude":"-3.615","obsSource":"LNDSYN","name":"St. Bees Head","region":"nw","unitaryAuthArea":"Cumbria"}
	,{"elevation":"81.0","id":"3212","latitude":"54.614","longitude":"-3.157","name":"Keswick","nationalPark":"Lake District National Park","region":"nw","unitaryAuthArea":"Cumbria"}
	,{"elevation":"28.0","id":"3220","latitude":"54.933","longitude":"-2.963","name":"Carlisle","region":"nw","unitaryAuthArea":"Cumbria"}
	,{"elevation":"285.0","id":"3224","latitude":"55.05","longitude":"-2.553","name":"Spadeadam","region":"nw","unitaryAuthArea":"Cumbria"}
	,{"elevation":"9.0","id":"3316","latitude":"53.497","longitude":"-3.056","obsSource":"LNDSYN","name":"Crosby","region":"nw","unitaryAuthArea":"Merseyside"}
	,{"elevation":"68.0","id":"3377","latitude":"53.175","longitude":"-0.521","name":"Waddington","region":"em","unitaryAuthArea":"Lincolnshire"}
	,{"elevation":"360.0","id":"3410","latitude":"52.757","longitude":"-3.464","name":"Lake Vyrnwy Saws","region":"wl","unitaryAuthArea":"Powys"}
	,{"elevation":"21.0","id":"3488","latitude":"52.949","longitude":"1.127","obsSource":"LNDSYN","name":"Weybourne","region":"ee","unitaryAuthArea":"Norfolk"}
	,{"elevation":"35.0","id":"3351","latitude":"53.3598","longitude":"-2.3805","name":"Rostherne No 2","region":"nw","unitaryAuthArea":"Cheshire East"}
	,{"elevation":"128.0","id":"3680","latitude":"51.8067","longitude":"-0.3602","name":"Rothamsted","region":"ee","unitaryAuthArea":"Hertfordshire"}
	,{"elevation":"140.0","id":"3080","latitude":"57.077","longitude":"-2.836","name":"Aboyne","region":"gr","unitaryAuthArea":"Aberdeenshire"}
	,{"elevation":"134.0","id":"3088","latitude":"56.852","longitude":"-2.264","name":"Inverbervie","region":"gr","unitaryAuthArea":"Aberdeenshire"}
	,{"elevation":"10.0","id":"3111","latitude":"55.441","longitude":"-5.699","obsSource":"LNDSYN","name":"Campbeltown Airport","region":"st","unitaryAuthArea":"Argyll and Bute"}
	,{"elevation":"15.0","id":"3026","latitude":"58.214","longitude":"-6.325","obsSource":"LNDSYN","name":"Stornoway","region":"he","unitaryAuthArea":"Na h-Eileanan Siar"}
	,{"elevation":"265.0","id":"3031","latitude":"57.725","longitude":"-4.896","name":"Loch Glascarnoch Saws","region":"he","unitaryAuthArea":"Highland"}
	,{"elevation":"773.0","id":"3039","latitude":"57.4175","longitude":"-5.689","name":"Bealach Na Ba","region":"he","unitaryAuthArea":"Highland"}
	,{"elevation":"87.0","id":"3808","latitude":"50.218","longitude":"-5.33","name":"Camborne","region":"sw","unitaryAuthArea":"Cornwall"}
	,{"elevation":"16.0","id":"3866","latitude":"50.577","longitude":"-1.297","name":"St Catherines Pt.","region":"se","unitaryAuthArea":"Isle of Wight"}
	,{"elevation":"4.0","id":"3872","latitude":"50.815","longitude":"-0.923","name":"Thorney Island","region":"se","unitaryAuthArea":"West Sussex"}
	,{"elevation":"52.0","id":"3882","latitude":"50.89","longitude":"0.319","name":"Herstmonceux West End","region":"se","unitaryAuthArea":"East Sussex"}
	,{"elevation":"49.0","id":"3904","latitude":"54.707","longitude":"-7.577","name":"Castlederg","region":"ni","unitaryAuthArea":"County Tyrone"}
	,{"elevation":"40.0","id":"3952","latitude":"51.799","longitude":"-8.25","obsSource":"LNDSYN","name":"Roches Point"}
	,{"elevation":"115.0","id":"99060","latitude":"53.85","longitude":"-2.467","name":"Stonyhurst","region":"nw","unitaryAuthArea":"Lancashire"}
	,{"elevation":"63.0","id":"3503","latitude":"52.344","longitude":"-3.947","name":"Trawsgoed","region":"wl","unitaryAuthArea":"Ceredigion"}
	,{"elevation":"307.0","id":"3507","latitude":"52.063","longitude":"-3.614","name":"Sennybridge","region":"wl","unitaryAuthArea":"Powys"}
	,{"elevation":"81.0","id":"3649","latitude":"51.758","longitude":"-1.576","name":"Brize Norton","region":"se","unitaryAuthArea":"Oxfordshire"}
	,{"elevation":"40.0","id":"3672","latitude":"51.548","longitude":"-0.415","name":"Northolt","region":"se","unitaryAuthArea":"Greater London"}
	,{"elevation":"6.0","id":"3707","latitude":"51.089","longitude":"-4.149","name":"Chivenor","region":"sw","unitaryAuthArea":"Devon"}
	,{"elevation":"25.0","id":"3772","latitude":"51.479","longitude":"-0.4491","name":"Heathrow","region":"se","unitaryAuthArea":"Greater London"}
	,{"elevation":"59.0","id":"3134","latitude":"55.907","longitude":"-4.533","name":"Glasgow\/Bishopton","region":"st","unitaryAuthArea":"Renfrewshire"}
	,{"elevation":"15.0","id":"3214","latitude":"54.125","longitude":"-3.257","name":"Walney Island","region":"nw","unitaryAuthArea":"Cumbria"}
	,{"elevation":"255.0","id":"3225","latitude":"54.501","longitude":"-2.684","name":"Shap","region":"nw","unitaryAuthArea":"Cumbria"}
	,{"elevation":"847.0","id":"3227","latitude":"54.684","longitude":"-2.45","name":"Great Dun Fell 2","region":"nw","unitaryAuthArea":"Cumbria"}
	,{"elevation":"298.0","id":"3330","latitude":"53.12755","longitude":"-1.97993","name":"Leek","region":"wm","unitaryAuthArea":"Staffordshire"}
	,{"elevation":"57.0","id":"3373","latitude":"53.307","longitude":"-0.546","name":"Scampton","region":"em","unitaryAuthArea":"Lincolnshire"}
	,{"elevation":"95.0","id":"3405","latitude":"52.789","longitude":"-4.742","name":"Aberdaron","region":"wl","unitaryAuthArea":"Gwynedd"}
	,{"elevation":"27.0","id":"3839","latitude":"50.737","longitude":"-3.405","name":"Exeter Airport","region":"sw","unitaryAuthArea":"Devon"}
	,{"elevation":"4.0","id":"3023","latitude":"57.358","longitude":"-7.397","obsSource":"LNDSYN","name":"South Uist Range","region":"he","unitaryAuthArea":"Na h-Eileanan Siar"}
	,{"elevation":"1130.0","id":"3041","latitude":"56.822","longitude":"-4.97","name":"Aonach Mor","region":"he","unitaryAuthArea":"Highland"}
	,{"elevation":"3.0","id":"3784","latitude":"51.464","longitude":"0.314","name":"Gravesend-Broadness","region":"se","unitaryAuthArea":"Kent"}
	,{"elevation":"54.0","id":"3797","latitude":"51.3422","longitude":"1.3461","obsSource":"LNDSYN","name":"Manston","region":"se","unitaryAuthArea":"Kent"}
	,{"elevation":"50.0","id":"3827","latitude":"50.354","longitude":"-4.121","obsSource":"LNDSYN","name":"Mount Batten","region":"sw","unitaryAuthArea":"Plymouth"}
	,{"elevation":"161.0","id":"3923","latitude":"54.237","longitude":"-6.502","name":"Glenanne","region":"ni","unitaryAuthArea":"County Armagh"}
	,{"elevation":"133.0","id":"3502","latitude":"52.139","longitude":"-4.571","name":"Aberporth","region":"wl","unitaryAuthArea":"Ceredigion"}
	,{"elevation":"107.0","id":"3544","latitude":"52.358","longitude":"-1.33","name":"Church Lawford","region":"wm","unitaryAuthArea":"Warwickshire"}
	,{"elevation":"158.0","id":"3275","latitude":"54.563","longitude":"-0.863","name":"Loftus (Samos)","region":"ne","unitaryAuthArea":"Redcar and Cleveland"}
	,{"elevation":"216.0","id":"3305","latitude":"53.093","longitude":"-3.941","name":"Capel Curig","nationalPark":"Snowdonia National Park","region":"wl","unitaryAuthArea":"Conwy"}
	,{"elevation":"3.0","id":"3469","latitude":"52.873","longitude":"0.141","name":"Holbeach","region":"em","unitaryAuthArea":"Lincolnshire"}
	,{"elevation":"26.0","id":"3017","latitude":"58.954","longitude":"-2.9","obsSource":"LNDSYN","name":"Kirkwall","region":"os","unitaryAuthArea":"Orkney Islands"}
	,{"elevation":"170.0","id":"3781","latitude":"51.303","longitude":"-0.09","name":"Kenley","region":"se","unitaryAuthArea":"Surrey"}
	,{"elevation":"101.0","id":"3894","latitude":"49.44","longitude":"-2.6","obsSource":"LNDSYN","name":"Guernsey","region":"sw","unitaryAuthArea":"Guernsey"}
	,{"elevation":"63.0","id":"3917","latitude":"54.664","longitude":"-6.224","obsSource":"LNDSYN","name":"Belfast International Airport","region":"ni","unitaryAuthArea":"County Antrim"}
	,{"elevation":"89.0","id":"3590","latitude":"52.123","longitude":"0.961","name":"Wattisham","region":"ee","unitaryAuthArea":"Suffolk"}
	,{"elevation":"32.0","id":"3609","latitude":"51.565","longitude":"-3.981","name":"Mumbles Head","region":"wl","unitaryAuthArea":"Swansea"}
	,{"elevation":"57.0","id":"3658","latitude":"51.62","longitude":"-1.097","name":"Benson","region":"se","unitaryAuthArea":"Oxfordshire"}
	,{"elevation":"35.0","id":"3144","latitude":"56.326","longitude":"-3.729","name":"Strathallan","region":"ta","unitaryAuthArea":"Perth and Kinross"}
	,{"elevation":"564.0","id":"3148","latitude":"56.423","longitude":"-4.32","name":"Glen Ogle","nationalPark":"Loch Lomond and the Trossachs National Park","region":"ta","unitaryAuthArea":"Stirling"}
	,{"elevation":"16.0","id":"3204","latitude":"54.0849","longitude":"-4.6321","obsSource":"LNDSYN","name":"Ronaldsway","region":"nw"}
	,{"elevation":"23.0","id":"3240","latitude":"55.421","longitude":"-1.6","obsSource":"LNDSYN","name":"Boulmer","region":"ne","unitaryAuthArea":"Northumberland"}
	,{"elevation":"25.0","id":"3265","latitude":"54.204","longitude":"-1.39","name":"Topcliffe","region":"yh","unitaryAuthArea":"North Yorkshire"}
	,{"elevation":"15.0","id":"3292","latitude":"54.094","longitude":"-0.174","name":"Bridlington Mrsc","region":"yh","unitaryAuthArea":"East Riding of Yorkshire"}
	,{"elevation":"77.0","id":"3313","latitude":"53.259","longitude":"-3.509","name":"Rhyl","region":"wl","unitaryAuthArea":"Denbighshire"}
	,{"elevation":"10.0","id":"3321","latitude":"53.174","longitude":"-2.986","name":"Hawarden","region":"wl","unitaryAuthArea":"Flintshire"}
	,{"elevation":"117.0","id":"3354","latitude":"53.005","longitude":"-1.25","name":"Watnall","region":"em","unitaryAuthArea":"Nottinghamshire"}
	,{"elevation":"7.0","id":"3382","latitude":"53.867","longitude":"-0.433","name":"Leconfield Sar","region":"yh","unitaryAuthArea":"East Riding of Yorkshire"}
	,{"elevation":"8.0","id":"3385","latitude":"53.473","longitude":"0.154","obsSource":"LNDSYN","name":"Donna Nook","region":"em","unitaryAuthArea":"Lincolnshire"}
	,{"elevation":"6.0","id":"3391","latitude":"53.094","longitude":"-0.171","name":"Coningsby","region":"em","unitaryAuthArea":"Lincolnshire"}
	,{"elevation":"110.0","id":"99142","latitude":"54.273","longitude":"-0.421","name":"Scarborough","region":"yh","unitaryAuthArea":"North Yorkshire"}
	,{"elevation":"27.0","id":"3844","latitude":"50.7366","longitude":"-3.40458","name":"Exeter Airport 2","region":"sw","unitaryAuthArea":"Devon"}
	]}
  };
  
  //Using this sort means that both items need to be shown to the user, eg 'Devon, Exeter Airport', going 
  //from the general to the specific. 
  var concatenated_name = function(site){  
	  return site.unitaryAuthArea + site.name;
  }
  
  var by_name = function(site_a, site_b){
		var a = concatenated_name(site_a);
		var b = concatenated_name(site_b);
    if(a < b) return -1;
    if(a > b) return 1;
    return 0;
  }
  
  var result = sites.Locations.Location; // the array of site data
  
  return result.sort(by_name);
};