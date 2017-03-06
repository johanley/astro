var show = function(text_output){
  console.log('Starting');
  
  var sites = uk_weather_sites(); //already has the desired sort-by-names  
  console.log('Num sites: ' + sites.length);
  var output_div = document.getElementById(text_output);
  
  var output_option = function(site){
    var name = site.name; 
    if (site.unitaryAuthArea){
      name = site.unitaryAuthArea + " - " + name;
    }
  	var line = "&lt;option value='" + site.id + "'&gt;" + name;
  	output_div.innerHTML = output_div.innerHTML + line + '\r\n'; 
  };
  
  for (var idx = 0; idx < sites.length; ++idx){
  	output_option(sites[idx]);
  }
	
  console.log('Done');
};