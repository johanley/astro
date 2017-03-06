/*
 *  Clearskyclock.com has index files that store data about each location.
 *  This module downloads that data, and finds the station id that is 
 *  closest to a given target latitude-longitude, and returns it.
 *  The caller uses it to build a url of the form ('FLO' is the station id):
 *     http://www.cleardarksky.com/c/FLOcsk.gif
 *  
 *  Source data url:
 *  http://www.cleardarksky.com/t/chart_keys00.txt
 *  
 *  Example line:
 *  AASPNB|New Brunswick|47.4319444444|-66.9125
 *  
 *   target: 
 *     .φ latitude degs
 *     .λ longitude degs
 *   fx_inform_caller: 
 *     callback function that takes a single string, the id of the CSC station that is closest 
 *     to the given latitude and longitude
 */
var find_closest_csc_url = function(target, fx_inform_caller){

  var to_rads = function(degs){
  	return degs * (Math.PI / 180.0); 
  };
	target.φ = to_rads(target.φ); 
	target.λ = to_rads(target.λ);
	
	var fetch_and_parse_and_tell_the_caller = function (){
      var url = UTIL.crossDomainUrl('http://www.cleardarksky.com/t/chart_keys00.txt&ext=txt');
      //var url = UTIL.crossDomainUrl('https://www.cleardarksky.com/t/chart_keys00.txt');
      var xhr = new XMLHttpRequest();
      xhr.onreadystatechange = function() {
        if (xhr.readyState == 4 && xhr.status == 200) {
           //console.log('Raw response text:' + xhr.responseText);
        	 console.log('Success: fetch station database from cleardarksky.com');
           do_all(xhr.responseText);
        }
      };    
      xhr.ontimeout = function (e) {
      	console.log('Timeout: fetch station database from cleardarksky.com');
          do_all(null);
      };
      xhr.open("GET", url, true);
      xhr.timeout = 20*1000; 
      xhr.send();
    };
    
	/* 
	 * Return an array of objects with .id, .λ, .φ (in rads).
     *  Example line:
     *  AASPNB|New Brunswick|47.4319444444|-66.9125
	 */
	var parsed_objects = function (raw_data_from_csc_server){
		var result = [];
		var lines = raw_data_from_csc_server.split("\n"); // array
		for (var i = 0; i < lines.length; ++i){
			var parts = lines[i].split("|"); 
			var line = {};
			line.id = parts[0];
			line.φ = to_rads(parts[2]);
			line.λ = to_rads(parts[3]);
			result.push(line);
		}
		return result;
    };
    
	/*
	 * Simple, fast max-min boxing, with no exact calculation of real distance. 
	 * Return an array of objects with .id, .λ, .φ (in rads). 
	 */
    var filtered_by_box = function(parsed_objects){
    	var i = 0;
    	var result = [];
    	var place;
    	var max_lat = to_rads(2);
    	var max_long = to_rads(4);
    	for (i = 0; i < parsed_objects.length; ++i){
    		place = parsed_objects[i];
    		if (Math.abs(place.φ - target.φ) < max_lat &&  Math.abs(place.λ - target.λ) < max_long){
    			result.push(place);
    		}
    	}
    	return result;
    };
    
	/*
	 * Calculate the real distance, treating the Earth as a sphere, and select the nearest 
	 * station id.
	 *  
	 * Return the Clear sky clock station id. 
	 */
    var filtered_by_min_distance = function(filtered_by_box){
    	var i = 0;
    	var result = "";
    	var min_dist = 1000; //kms
    	var place;
    	for (i = 0; i < filtered_by_box.length; ++i){
    		place = filtered_by_box[i];
    		var dist = EPH.distance_kms(place, target);
    		if (dist < min_dist){
    			min_dist = dist;
    			result = place.id;
    		}
    	}
    	return result;
    };

    var do_all = function(raw_data){
        var result = "";
        if (raw_data){
        	var base_data = parsed_objects(raw_data); //array of objects
        	var boxed_data = filtered_by_box(base_data); //min-max lat-long
        	if (boxed_data.length > 0){
        		result = filtered_by_min_distance(boxed_data); //station id
        	}
        }
        fx_inform_caller(result);
    };
    
    fetch_and_parse_and_tell_the_caller();
}