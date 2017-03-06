/*
 Fetch aurora data, current Kp, a global measure of activity. 

  callback_Kp
    callback function, to which is passed an array of objects containing:
     .datetime (as a string, format 'UT 2016-11-10 03', for example, denoting 03h00 on Nov 10)
     .values, an array of integers, Kp, in the range 0..9, a global mean measure of magnetic field fluctuations for the past two days
       
  Credit:
    http://www.swpc.noaa.gov/products/station-k-and-indices
    ftp://ftp.swpc.noaa.gov/pub/lists/geomag/AK.txt - Kp estimates, past 2 days
    http://services.swpc.noaa.gov/text/aurora-nowcast-map.txt - nowcast
*/
var fetch_aurora_data = function(callback_Kp){

  /* 
    Example input: '2016 Oct 9'.
    Return: 'UT 2016-10-09'. 
  */
  var canonical_date_Kp = function(raw_date){
    var date = raw_date.trim();
    var year = date.substring(0,4).trim();
    var month_raw = date.substring(4,8).trim();
    var map_months = { 
      Jan: '01',  Feb: '02',  Mar: '03', Apr: '04',  May: '05', Jun: '06',  
      Jul: '07',  Aug: '08',  Sep: '09', Oct: '10', Nov: '11',  Dec: '12'
    };
    var month = map_months[month_raw];
    var day = date.substring(8).trim();
    var result = 'UT ' + year + '-' + month + '-' + pad(day);
    return result;
   };
   
  var pad = function(text){
    return text.length === 1 ? '0' + text : text;  //concatenation, not addition
  };
  
  var extract_Kp_measurements = function(date, line, result /*out-param*/){
    //   'Planetary(estimated Ap)      5     2     1     2     1     2     1     2     1'   -- starts with same text
    var numbers = line.substring('Planetary(estimated Ap)      5'.length); //drop the leading text, and the Ap value
    var values = numbers.trim().split(/\s{5}/); //should be 8 items (ints-as-strings, at this stage)
    for (var i = 0; i < values.length; ++i){
      result.push({
        datetime: date + ' ' + pad(i*3 + ''), //strings, not numbers
        value: parseInt(values[i])
      });
    }
  };
  
  /*
   Return an array of N objects, each having these properties:
     .datetime (as a string, UT, eg '2016-10-11 03' meaning beginning at 03h00)
     .value, 8 integers in the range 0..9
        each value corresponds to a 3 hour period (3*8=24)
        if the value is -1, that means no data is available for that time period, and it should be ignored by the caller
   Data: ftp://ftp.swpc.noaa.gov/pub/lists/geomag/AK.txt
      # Station        Lat Long  Index 00-03 03-06 06-09 09-12 12-15 15-18 18-21 21-24    -- header shows the hours, no need to parse it
     '2016 Oct 9   '    --- starts with 4 numbers
     'Planetary(estimated Ap)      5     2     1     2     1     2     1     2     1'   -- starts with same text
  */
  var parse_Kp_measurements = function(raw_data){
    var result = [];
    var lines = raw_data.split("\n"); // array
    var line, date;
    for (var i = 0; i < lines.length; ++i){
      line = lines[i].trim();
      if (line.match(/^\d{4}\.*/)) {
        date = canonical_date_Kp(line);
      }
      else if (line.startsWith('Planetary(')){
        extract_Kp_measurements(date, line, result);      
      }
    }
    return result;
  };
  
  var do_Kp_measurements = function(){
    var url = UTIL.crossDomainUrl('ftp://ftp.swpc.noaa.gov/pub/lists/geomag/AK.txt&ext=txt');
    var xhr = new XMLHttpRequest();
    xhr.onreadystatechange = function() {
      if (xhr.readyState == 4 && xhr.status == 200) {
        console.log('Success: fetch Kp measurements from NOAA, ' + url);
        var parsed_data = parse_Kp_measurements(xhr.responseText);
        callback_Kp(parsed_data);
      }
    };    
    xhr.ontimeout = function (e) {
      console.log('Timeout: fetch Kp measurements from NOAA, ' + url);
      //don't call anything
    };
    xhr.open("GET", url, true);
    xhr.timeout = 20*1000; 
    xhr.send();
  };

  do_Kp_measurements();
  
  //////////////

  // The 30m prediction (below) is NOT WORKING WELL. SLOW. LITTLE VALUE. BROKEN PIXEL CALC? Leave out for the moment.


  /* 
   Return array of array of ints-as-strings; first level is lat, second level is long. The int value is the auroral prediction data, in the range 0..100.
     array[512][1024] 
     array[512 (-90..+90, latitude)][1024 (0..360, longitude, increasing to the East, I assume)] 
  */
  var parse_global_aurora_predictions_30m = function(raw_data){
    var result = [];
    var lines = raw_data.split("\n"); // array of 512 lines, -90..+90 
    var line = '';
    for (var i = 0; i < lines.length; ++i){
      line = lines[i].trim();
      if (! line.startsWith('#') /*comments*/ && line.length > 0 /*not-empty*/) {
        //the numbers are separatedby 1..3 blank spaces
        var values_for_longitudes = line.split(/\s{1,3}/);  // array of 1024 items, one for each longitude-pixel, 0..1023 eg: '   0   0   19   21  ...'  fixed width format  
        result.push(values_for_longitudes); // each line corresponds to a single latitude level; this is singular at the poles, of course: all then values will necessarily be the same
      }
    }
    return result;  
  };
  
  /* 
   Return pixel.λ, pixel.φ, with each being an INTEGER index into the base data structure parsed above.
   Ambiguous: where the pixel starts on the ground. Does the first pixel-in-longitude start at 0 degrees, or is it centered at 0 degrees?
   In practice, the predictions aren't that precise, so it doesn't matter all that much. But it would be nice to know. 
  */
  var calc_aurora_pixel_for_30m = function(λ, φ /*the input location, in degrees, with longitude positive east of Greenwich. */){
    var result = {};
    var frac_lat = (φ + 90)/180; 
    result.φ = Math.floor(frac_lat * 512); //truncated integer; close enough; which exact pixel probably isn't significant
    if (result.φ === 512){
      result.φ = 511; //avoid overflow
    }
    var longit = λ < 0 ? 360 + λ : λ; 
    var frac_long = longit/360; 
    result.λ = Math.floor(frac_long * 1024); 
    if (result.λ === 1024){
      result.λ = 1023; //avoid overflow
    }
    return result;
  };
  
  /* 
   Return an int in range 0..100, representing the current probability of seeing visible aurora at a 
   given place on the globe. 
  */
  var extract_local_aurora_prediction_30m = function(raw_data){
    if (raw_data){
      var global_aurora_predictions = parse_global_aurora_predictions_30m(raw_data);
      var pixel = calc_aurora_pixel_for_30m(where.λ, where.φ);
      var result = global_aurora_predictions[pixel.φ][pixel.λ];
      console.log('Module, current probability of aurora: ' + result + '%');
      //link to this page for a graphic: http://www.swpc.noaa.gov/products/aurora-30-minute-forecast
      //title: Probability of visible aurora?
      // Kp graphs: http://www.swpc.noaa.gov/products/station-k-and-indices
      // specific Canadian magnetometers: http://www.carisma.ca/about-carisma 
      // table of yearly pole positions: http://wdc.kugi.kyoto-u.ac.jp/poles/polesexp.html  (IGRF-2)
      // table of yearly pole positions (better: one more decimal place) : http://www.geomag.bgs.ac.uk/education/poles.html
      
      var aurora_data = {
        prob_30m : result, 
        geomagnetic_latitude : EPH.geomagnetic_latitude(where, when) 
      };
      callback_30m_prob(aurora_data); 
    }
  };
  
  var fetch_global_30m_predictions = function(callback){
    var url = UTIL.crossDomainUrl('http://services.swpc.noaa.gov/text/aurora-nowcast-map.txt&ext=txt');
    var xhr = new XMLHttpRequest();
    xhr.onreadystatechange = function() {
      if (xhr.readyState == 4 && xhr.status == 200) {
        console.log('Success: fetch aurora predictions from NOAA, ' + url);
        callback(xhr.responseText);
      }
    };    
    xhr.ontimeout = function (e) {
      console.log('Timeout: fetch aurora predictions from NOAA, ' + url);
      callback('');
    };
    xhr.open("GET", url, true);
    xhr.timeout = 20*1000; 
    xhr.send();
  };
  
  // do all the work for the 30m predictions 
  //fetch_global_30m_predictions(extract_local_aurora_prediction_30m);
}