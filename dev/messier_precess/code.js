/*
  Take the contents of the Messier catalog, 
  and apply precession to the star positions, from J2000.0 to the 
  given target epoch.
*/
var show = function(input, text_output){
  console.log('Starting');

  EPH.testing();
  console.log('Size of Messier catalog: ' + messier.length);
  var output_div = document.getElementById(text_output);
  
  var round = function(num){
	  var num_decimals = 10000000;
	  return Math.round(num*num_decimals)/num_decimals;
  };
  
  var output_result = function(nebula, idx, name, mag, constellation, type, comment, commonName){
    //messier[0]=['M1',1.4595316,0.3842633,8.4,"Crab Nebula"];
	//messier[0]=['M1',1.4595316,0.3842633,8.4,"Tau","","","Crab Nebula"];
	var line = 'messier['+idx+']=["'+name+'",' +round(nebula.α)+ ',' +round(nebula.δ)+ ',' +mag+ ',"' +constellation+ '","' +type+ '","' +comment+ '","' +commonName+ '"]';
	//console.log(line);
	output_div.innerHTML = output_div.innerHTML + line + ';' + '\r\n'; 
  };
  
  //var target = parseFloat(input.target_epoch);
  var when_target = EPH.when(input.target_epoch);
  var when_source = EPH.when_j2000;
  var angles = EPH.precession_angles(when_source, when_target);
  
  for (var idx = 0; idx < messier.length; ++idx){
	var nebula = {
	  α: messier[idx][1],
	  δ: messier[idx][2],
	  equinox: when_source
	};
	EPH.apply_precession(nebula, when_target, angles);
	output_result(nebula, idx, messier[idx][0], messier[idx][3], messier[idx][4], messier[idx][5], messier[idx][6],	messier[idx][7]);
  }
  
  console.log('Done');
};