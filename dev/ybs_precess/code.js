/*
  Take the contents of the Yale Bright Star catalog, 
  and apply precession to the star positions, from J2000.0 to the 
  given target epoch.
  This takes a long time; Chrome will barf if it can't get to the end.
  Simple: split the calc into 2 pieces. 
  Better: work with the browser to let it know something is still happening:
    http://stackoverflow.com/questions/6336844/setting-javascript-timeout-limit-in-google-chrome 
*/
var show = function(input, text_output){
  console.log('Starting');

  EPH.testing();
  console.log('1. Size of Yale Bright Star catalog: ' + ybs.length);
  var output_div = document.getElementById(text_output);
  
  var round = function(num){
	  var num_decimals = 10000000;
	  return Math.round(num*num_decimals)/num_decimals;
  };
  
  var output_result = function(star, idx, name, mag){
    //ybs[0]=['',0.0225366,0.7893979,6.70];
	var line = "ybs["+idx+"]=['"+name+"',"+round(star.α)+","+round(star.δ)+","+mag+"]";
	//console.log(line);
	output_div.innerHTML = output_div.innerHTML + line + ';' + '\r\n'; 
  };
  
  var when_target = EPH.when(input.target_epoch); // eg 'J2016.5'
  var when_source = EPH.when_j2000;
  var angles = EPH.precession_angles(when_source, when_target);
  
  for (var idx = 0; idx < ybs.length; ++idx){
	var star = {
	  α: ybs[idx][1],
	  δ: ybs[idx][2],
	  equinox: when_source
	};
	EPH.apply_precession(star, when_target, angles);
	output_result(star, idx, ybs[idx][0], ybs[idx][3]);
  }
  
  console.log('Done');
};