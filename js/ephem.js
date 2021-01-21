/*
 Javascript library of basic, low-precision astronomical calculations.
 
 The intent of this tool is to provide reasonably accurate positions within a year or so of the current date.  
 If you need arc-second precision for the year 2000 B.C., you won't find it here.
 (Higher precision over longer time periods requires significantly more effort.)
 
 The data used by this library requires periodic updating in order to retain its accuracy.
 
 This module places a single javascript object called 'EPH' in global scope.
 (A convention in Javascript is to name global objects using capital letters.)

 Different coordinate systems are used in different cases.
 To make the current sky, you need the mean-equinox-equator-of-date (or even the apparent one).
 When presenting stars in a planisphere, it's likely best to 'pre-precess' catalog positions to some nearby date.
 ----------------------------------------------------------
 mean-equinox     j2000          j2018.5 (or similar)
 ----------------------------------------------------------
 sun              planets        stars
 moon             minor planets  messier
                  meteor showers
 ----------------------------------------------------------

 References
 Astronomical Algorithms, Jean Meeus, 1st edition, 1991
   Unless otherwise noted, algorithms are based on Meeus.
   
 Explanatory Supplement to the Astronomical Almanac and Ephemeris
   an old edition is here:
   https://ia600301.us.archive.org/28/items/astronomicalalmanac1961/131221-explanatory-supplement-1961.pdf
 
 Units
 * When it seems helpful, units can are sometimes appended to variable names, as in blah_km to denote kilometers
 * All angles are in radians, unless otherwise specified.
 * Distances are in different units depending on context (mostly AU, sometimes km).
 * Dates and times are unusual in that different algorithms use different ways of 
   describing a moment in time. Interconversion between these is simple but annoying. 
   In an attempt to make things easier for the caller, this library makes use of a 'when' 
   object, which calculates all of the usual time measures at once, and stores them as properties 
   on a single 'when' object.
   This library uses UT (UT1, to be precise) and TT, but never local time zone data.
   This library does no conversions for local time zone. 
   Any operations regarding local time zone are left to the caller. 
   Note that in Javascript, a Date object already carries both UT and local time zone data.
   Does a JS Date refer to UT1 or UTC? The current date will refer to UTC, civil time.
   But JS has no notion of leap second, so the distinction may be impossible to make. 
 * Longitude is taken negative west of Greenwich. This is contrary to the style adopted by Meeus, but in agreement 
   with many modern tool, such as Google Maps.

---------------------------------------------------------------------------------------------
  Central to this tool are these conventional objects, having important data:
  
  All the aliases for a given moment in time.
  **Items refer to UT**, unless 'tt' or 'lt' appears in the name.
  The given .date property is provided in case you need values in LT.
  There is no function for converting from one offset/timezone to another. 
  when {
    T : number of Julian centuries since J2000
    T_tt : number of Julian centuries since J2000, expressed in Terrestial Time
    d : day of the month
    d_frac : fractional day of the month, with hours-min-sec expressed as a decimal
    date : date object; a back-door, in case LT values are needed
    gmst : Greenwich Mean Sidereal Time, 0..2pi
    hour : hour of the day 0..23
    jd : Julian Date
    jd_tt : Julian Date in Terrestial Time
    m : month of the year 1..12
    min : minute 0..59
    mjd : Modified Julian Date
    msec : milliseconds
    msec_epoch : milliseconds since the epoch used by Javascript, 1970-01-01.0
    sec : seconds 0..59
    weekday : 1..7, in UT, not LT
    y : calendar year, eg 1957
    
    delta(secs) : return a new 'when', that differs from this 'when' by the given number of seconds
    next() : return a new 'when', 24h ahead of this one
    prev() : return a new 'when', 24h behind this one
    startOfDayLT() : return a new 'when', corresponding to the start of the local day
    endOfDayLT() : return a new 'when', corresponding to the end of the local day
    toString(), and variations
  }
  
  The location of an observer on the Earth's surface (height is neglected - low precision)
  where {
    φ - latitude, rads
    λ - longitude, rads
    limiting_mag - limiting magnitude at the observatory site (used for meteor showers, hourly rate)
    is_topocentric - true or false; if false, then taken as geocentric
  }
  // the orbit of an object; different subsets of these items can specify an orbit
  orbit {
    equinox - the coord system used by the coords
    epoch - the moment for which the (osculating) orbit is valid
    a - length of the semi-major axis, in astronomical units (AU); elliptical only
    q - perihelion distance; needed only if parabolic/hyperbolic
    e - eccentricity
    i - inclination of orbit to the plane of the ecliptic
    Ω - longitude of ascending node
    π - longitude of perihelion
    ω - argument of perihelion
    L0 - longitude at epoch
    M0 - mean anomaly at epoch
    n - rads per day, mean motion; elliptical only
    T - time of perihelion passage, expressed as a 'when'; might be the same as the epoch, in some cases
    P - period in fractional days; elliptical only
  }
  //control aspects of the final output in various ways
  options {
    where: where(45,-75), --- no default; if present, (a,A) are added to the position
    equinox: when('...'), -- default is mean equinox of date?? no, default is whatever is calculated; it depends
    units: 'degs', -- default is 'rads'
    time_scale: 'LT' -- (LT|UT|TT), for output of formatted 'when' objects
    --rounding: 2   -- default is no rounding; leave this out??
    precession_angles: blah   -- performance optimization; can this be handled more elegantly?
  }
  //celestial coordinates of an object in the sky
  //the scheme mostly follows that of the Explanatory supplement
  //ξ,ζ,η are not used for planets; just use XYZ, same as for the Sun
  ephem {
   equinox - the equinox to which the coordinates refer
   α,δ,Δ - geocentric equatorial coords (Δ is in AU, including for the Moon) 
   λ,β,Δ - geocentric ecliptic coords
   X,Y,Z - geocentric rectangular coords (both Sun and Planets)
   x,y,z - heliocentric rectangular coords
   l,b,r - heliocentric ecliptic coords ???? still being used? where?
   A,a,h - topocentric: azimuth, altitude, local hour angle; different: these only make sense for the equinox of date
   elong - elongation from the sun, 0..pi 
   phase - phase angle, Sun-thing-Earth, 0..pi
   illum - illuminated fraction of the apparent disk, 0..1
   mag   - apparent magnitude
   size  - apparent angular size of the disk, seconds of arc
  }
  
 Consider: memoization is very useful; should use for precession angles, and likely other things too.
 Consider: if 'ephem' always points to a when, then passing ephem+when is unwanted
 Consider: can N conversion functions be passed as varargs params, to be done in sequence?
 Consider: add heliocentric xyz coord conversion (from lbr?)

 ---------------------------------------------------------------------------------
 Nice example of a well done planisphere (js) : 
   http://www.etwright.org/astro/plani.html
   http://freestarcharts.com/
   
 Another js lib:
  https://github.com/mivion/ephemeris

 The old Explanatory Supplement:
   https://ia600301.us.archive.org/28/items/astronomicalalmanac1961/131221-explanatory-supplement-1961.pdf
   comet mag: page 132
   mag, phase of the planets: pg 311
   moonrise, set: pg 403; can have events 0, 1, or 2 events (high latitudes)! (sun always has 2)

 Ephemerides, for comparison:
 http://aa.usno.navy.mil/data/index.php 
 http://ssd.jpl.nasa.gov/?planet_pos   max err 10' for Saturn, 1800-2050
 http://ssd.jpl.nasa.gov/txt/aprx_pos_planets.pdf
 http://ssd.jpl.nasa.gov/horizons.cgi#top
 http://astropixels.com/ephemeris/planets/mercury2016.html    
 http://theskylive.com/ceres-info - nice!!
 http://aa.usno.navy.mil/data/docs/mrst.php -- good rise and set times for various objects
 
 Yale Bright Star r5 catalog
  http://cdsarc.u-strasbg.fr/viz-bin/Cat?V/50
 Messier catalog
   http://astropixels.com/messier/messiercat.html
 Caldwell catalog
   http://astropixels.com/caldwell/caldwellcat.html 
 Minor planets
   https://www.ast.cam.ac.uk/~jds/
   http://www.minorplanetcenter.net/iau/Ephemerides/EphemOrbEls.html
   http://www.minorplanetcenter.net/iau/info/CometOrbitFormat.html
   For mag estimates, I use M1, K1 from JPL (I don't grok the MPC system; not explicit).
   I could use the JPL for both orbit and mag, if I wanted
     http://ssd.jpl.nasa.gov/sbdb.cgi#top
 Comets
  https://www.ast.cam.ac.uk/~jds/  - what to list, brightest only
  http://ssd.jpl.nasa.gov/sbdb.cgi#top   - orbit data
 Meteor showers
   http://www.imo.net/files/data/vmdb/vmdbrad.txt   --  the radiant catalog, 2012-01-17
   http://www.imo.net/data
   http://imo.net/files/data/calendar/cal2016.pdf
   http://www.imo.net/calendar/2014
 Dominic Ford, BAA
  https://in-the-sky.org/about.php
 Position of the geomagnetic north pole:
    British Geological Survey
    http://www.geomag.bgs.ac.uk/education/poles.html
    
-----------------------------------------------------------------------
  
 This implementation uses an 'immediately-invoked function expression' (IIFE) pattern.
 This is done to put 1 item in global scope, instead of N.
 To see the data visible to a user of the EPH object, go to the bottom of this file.
*/
var EPH = (function(){ 

  //START OF PRIVATE ITEMS
  
  /* 
   Douglas Crockford, Javascript: The Good Parts, page 22.
   Clone an object, and use prototypal inheritance.
   This is used, for example, to share functions among similar objects such as minor planets, with no code repetition.   
  */
  Object.create = function(thing){
    var ConstructorFunc = function(){};
    ConstructorFunc.prototype = thing;
    return new ConstructorFunc();
  };
  
  /* Change an object containing similar things into an array of similar things. */
  var as_array = function(thing){
    var result = [];
    for (prop in thing){
      if (thing.hasOwnProperty(prop)){
        result.push(thing[prop]);
      }
    }
    return result;
  };
  
  //Tested on: Chrome49, FF45, IE11.
  //IE 11 needs these; not tested on any other versions of IE. 
  var add_polyfills = function(){
    if (!String.prototype.startsWith) {
        String.prototype.startsWith = function(searchString, position){
          position = position || 0;
          return this.substr(position, searchString.length) === searchString;
      };
    }
    Math.log10 = Math.log10 || function(x) {
      return Math.log(x) / Math.LN10;
    };
    Math.trunc = Math.trunc || function(x) {
      return x < 0 ? Math.ceil(x) : Math.floor(x);
    }          
    Math.sign = Math.sign || function(x) {
      x = +x; // convert to a number
      if (x === 0 || isNaN(x)) {
        return x;
      }
      return x > 0 ? 1 : -1;
    }    
  };  
  add_polyfills();
  
  var MSEC_PER_DAY = 1000*60*60*24;
  var SEC_PER_DAY = 60*60*24;
  var INTERPOLATION_PROPS = ['α', 'δ', 'Δ', 'size']; //a common choice
  
  /* Julian date of J2000.0. */
  var JD_J2000 = 2451545.0; 
  
  /* Left-pad with a single '0' if 9 or less. This is meant especially for dates. */
  var pad = function(number){
    var padding = number < 10 ? '0' : '';
    return padding + number; 
  };
  
  var rads = function(deg){
     return deg * Math.PI/180;
  };
  
  var degs = function(rad){
     return rad * 180/Math.PI;
  };
  
  /*
   'places' can be negative.  
   No padding is applied if the result of rounding ends in a 0. 
  */
  var round = function(num, places){
    var factor = Math.pow(10, places);
    return Math.round(num*factor)/factor;
  };
  
  /*
   By default, js will round '8.0' to '8'; that is, it will drop any trailing 0's.
   This method will work around that, by returning a string having a fixed number of decimals. 
  */
  var round_and_pad = function(num, places){
    var val = round(num, places) + ''; //coerce to a string
    var decimal_point, num_decimals_present;
    if (places > 0){
      decimal_point = val.indexOf('.');
      if (decimal_point === -1){
        val = val + '.';
        decimal_point = val.indexOf('.');
      }
      num_decimals_present = (val.length - 1) - decimal_point;
      while (num_decimals_present < places) {
        val = val + '0';
        num_decimals_present = (val.length) - 1 - decimal_point;
      }
    }
    return val;
  };
  
  /* To the nearest minute of arc (since this is a low-precision library). */
  var degs_sexagesimal = function(rads) {
    var deg_decimal = degs(rads);
    var sign = deg_decimal < 0 ? -1 : 1;
    var d = Math.abs(deg_decimal);
    var degrees = Math.trunc(d);
    var minutes = Math.trunc((d - degrees)*60);
    var seconds = Math.round((d - (degrees + (minutes/60)))*3600);
    //let's give a nice convenient toString to the result object
    var result = {
      sign: sign,
      deg : degrees,
      min : minutes,
      sec: seconds,
      toString: function(){
        return (sign > 0 ? '+' : '-') + pad(degrees) + "° " + pad(minutes) + "' " + pad(seconds) + "''";
      }
    };
    return result;
  };
  var degs_sexagesimal_OLD = function(rads) {
    var deg_decimal = degs(rads);
    var sign = deg_decimal < 0 ? -1 : 1;
    var d = Math.abs(deg_decimal);
    var degrees = Math.trunc(d);
    var minutes = Math.round((d - degrees)*60);
    //let's give a nice convenient toString to the result object
    var result = {
      sign: sign,
      deg : degrees,
      min : minutes,
      toString: function(){
        return (sign > 0 ? '+' : '-') + pad(degrees) + "° " + pad(minutes) + "'";
      }
    };
    return result;
  };
  
  /* Always positive. */
  var in360 = function(deg){
    var result = deg % 360;
    if (result < 0){
      result = result + 360;
    }
    return result;
  };
  
  /* Always positive. */
  var in2pi = function(rads){
    var twopi = 2*Math.PI;
    var result = rads % twopi;
    if (result < 0){
      result = result + twopi;
    }
    return result;
  };
  
  var zodiac_sign = function (name, abbr, symbol, α_end_hour, α_end_min, λ_end){
    return  {
      name : name,
      abbr: abbr,
      symbol: symbol,
      α_end : rads((α_end_hour + α_end_min/60)*15),
      λ_end : rads(λ_end)
    };
  };
  
  var zodiac = [
    zodiac_sign('Pisces', 'Psc', '♓', 2, 0, 31),
    zodiac_sign('Aries', 'Ari', '♈', 3, 20, 53),
    zodiac_sign('Taurus', 'Tau', '♉', 5, 45, 88),
    zodiac_sign('Gemini', 'Gem', '♊', 8, 0, 118),
    zodiac_sign('Cancer', 'Cnc', '♋', 9, 25, 138),
    zodiac_sign('Leo', 'Leo', '♌', 11, 30, 171),
    zodiac_sign('Virgo', 'Vir', '♍', 14, 20, 218),
    zodiac_sign('Libra', 'Lib', '♎', 15, 40, 239),
    zodiac_sign('Scorpius', 'Sco', '♏', 17, 40, 268),
    zodiac_sign('Sagittarius', 'Sgr', '♐', 20, 0, 295),
    zodiac_sign('Capricorn', 'Cap', '♑', 21, 50, 327),
    zodiac_sign('Aquarius', 'Aqr', '♒', 23, 30, 351),
    zodiac_sign('Pisces', 'Psc', '♓', 23, 59.99999999999, 359.999999999999)
  ];
  
  //planet: name, semidiameter at standard distance, symbol
  
  var where = function(lat_degs, long_degs /* negative west */, limiting_mag, is_topocentric /*boolean*/){
    var result = {
      φ : rads(lat_degs),
      λ : rads(long_degs),
      limiting_mag: limiting_mag,
      is_topocentric: is_topocentric
    };
    return result;
  };
  
  // THE FUNCTIONS BELOW CANNOT TAKE 'when' AS AN ARG, SINCE THE 'when' OBJECT HASN'T BEEN CREATED YET
  
  var fractional_days = function(date_utc){
    return date_utc.getUTCDate() + date_utc.getUTCHours()/24 + date_utc.getUTCMinutes()/(24*60) + date_utc.getUTCSeconds()/(24*60*60) + date_utc.getUTCMilliseconds()/(24*60*60*1000); 
  };

  /*
   Julian date.
   Example: for 1957-10-4.81, JD=2436116.31
   Ref: Astronomical Algorithms, Jean Meeus.
   y, m, d: year, month (1..12), and day (0..31)
   The day is allowed to have a decimal portion.
  */ 
  var find_julian_date = function(y, m, d){
     if (m<3){
       y = y - 1;
       m = m + 12;
     }
     var a = Math.floor(y/100);
     var b = 2 - a + Math.floor(a/4);
     var result = Math.floor(365.25*(y+4716)) + Math.floor(30.6001*(m+1)) + d + b -1524.5;
     return result;
  };
  
  /* Returned units: fractional days.  */
  var convert_jd_utc_to_jd_tt = function(jd_utc, year){
    var delta_t_sec = delta_t(year); 
    return jd_utc + delta_t_sec/(60*60*24);
  };

  /*
   ΔT = TT - UT1
   Returns a result in seconds.
   Uses polynomial fits, for the years 1900..2150.
   Uses a hard-coded value for 2000, to use the measured value at j2000.
   For a low-precision result, the exact value of delta_t is not critical.
   Ref: http://maia.usno.navy.mil/
   Ref: http://eclipse.gsfc.nasa.gov/SEhelp/deltatpoly2004.html
  */  
  var delta_t = function(y){
    var result = 68; //current value at time of writing; a default, for safety

    //first the polynomial approximation    
    if (y >= 2150){
     var  u = (y-1820)/100;    
     result = -20 + 32 * Math.pow(u,2);
    }
    else if (y >= 2050){
      result = -20 + 32 * Math.pow((y-1820)/100, 2) - 0.5628 * (2150 - y);    
    }
    else if (y >= 2005){
     var t = y - 2000;    
     result = 62.92 + 0.32217 * t + 0.005589 * Math.pow(t,2);
    }
    else if (y >= 1986){
      var t = y - 2000;    
      result = 63.86 + 0.3345 * t - 0.060374 * Math.pow(t,2) + 0.0017275 * Math.pow(t,3) + 0.000651814 * Math.pow(t,4) + 0.00002373599 * Math.pow(t,5);
    }
    else if (y >= 1961){
      var t = y - 1975;    
      result = 45.45 + 1.067*t - Math.pow(t,2)/260 - Math.pow(t,3) / 718;
    }
    else if (y >= 1941){
      var t = y - 1950;    
      result = 29.07 + 0.407*t - Math.pow(t,2)/233 + Math.pow(t,3) / 2547;
    }
    else if (y >= 1920){
      var t = y - 1920;    
      result = 21.20 + 0.84493*t - 0.076100 * Math.pow(t,2) + 0.0020936 * Math.pow(t,3);
    }
    else {
      var t = y - 1900;    
      // uses the estimate from 1900..1920 as the value for all other past dates as well, as a simple approximation
      result = -2.79 + 1.494119 * t - 0.0598939 * Math.pow(t,2) + 0.0061966 * Math.pow(t,3) - 0.000197 * Math.pow(t,4);
    }
    
    //for j2000, overwrite with the actual measurement
    if (y === 2000){
      //http://maia.usno.navy.mil/ser7/deltat.data
      result = 63.828; 
    }
    return result;
  };
  
  /* Returns an integer. Sunday=1, Monday=2, etc. */
  var day_of_the_week = function(y, m, d /*no fraction*/){
    var jd = find_julian_date(y, m, d);
    var result = (jd + 1.5) % 7;
    return result + 1;
  };
  
  /* Result is in 0..2pi. jd is in UT. */
  var greenwich_mean_sidereal_time = function(jd, T){
   var result = 280.46061837 + 360.98564736629*(jd - JD_J2000) + 0.000387933*(T*T) - (T*T*T)/38710000;
   result = in2pi(rads(result));
   return result;
  };
  
  // THE FUNCTIONS BELOW CAN TAKE 'when' AS AN ARG

  /* Mean value. No nutation. Rads. */  
  var obliquity_of_ecliptic = function(when){
    var T = when.T;
    var deg = 23 + 26/60 + 21.448/3600 - (46.815/3600)*T - (0.00059/3600)*T*T + (0.001813/3600)*T*T*T;
    return rads(deg);
  };
  
  /* Result in rads. */ 
  var find_mean_sidereal_time_at_longitude = function(when, where){
    var greenwich = greenwich_mean_sidereal_time(when.jd, when.T); //rads
    var result = greenwich + where.λ; //longitude is negative west of Greenwich here; opposite to Meeus
    return in2pi(result);
  };
  
  /* Return .hour .min .sec. Does not return a when object. Numbers below 10 are left-padded with a 0. */
  var rads_to_time = function(angle_rads){
    var degrees = degs(angle_rads);
    var hours = degrees/15;
    var h = Math.trunc(hours);
    var h_frac = hours - h;
    var minutes = h_frac * 60;
    var m = Math.trunc(minutes);
    var s = (minutes - m) * 60;  
    return {
      hour: pad(h),
      min: pad(m), 
      sec: pad(s) //can have a decimal
    };
  };
  
  /* Return local mean sidereal time, as .hour .min .sec. Does not return a when object. Numbers below 10 are left-padded with a 0. */
  var lmst = function(when, where /*rads*/){
    var angle_rads = find_mean_sidereal_time_at_longitude(when, where); //0..2pi
    return rads_to_time(angle_rads);
  };  
    
  /* Result in rads, 0..2pi. Longitude in rads. 'ephem' has an .α in rads. */
  var find_local_hour_angle = function(ephem, when, where){
    var lst = find_mean_sidereal_time_at_longitude(when, where);
    var result = lst - ephem.α;
    return in2pi(result);
  };
  
  /* 
   For a given instant, return an object which stores all of the common ways of describing that instant.
   You can think of this as calculating N aliases for the same moment in time. 
  */
  var when_from_utc = function(date_utc, text /*optional*/){
   //note that the functions called here can't take a 'when' object, since it has not yet been created
   var year = date_utc.getUTCFullYear();
   var month = date_utc.getUTCMonth() + 1;
   var day_frac = fractional_days(date_utc);   
   var julian_date = find_julian_date(year, month, day_frac);
   var julian_date_tt = convert_jd_utc_to_jd_tt(julian_date, year);
   var T_centuries = (julian_date - JD_J2000)/36525;
   var result = {
      y : year,
      m : month,
      d : date_utc.getUTCDate(),
      d_frac : day_frac,
      hour : date_utc.getUTCHours(),
      min : date_utc.getUTCMinutes(),
      sec : date_utc.getUTCSeconds(),
      msec : date_utc.getUTCMilliseconds(),
      msec_epoch: date_utc.getTime(),
      weekday: day_of_the_week(year, month, date_utc.getUTCDate()), //1 is Sunday - note this is UT timescale, not local
      jd : julian_date,
      mjd :  julian_date - 2400000.5,
      jd_tt : julian_date_tt,
      T: T_centuries, // Julian centuries since J2000
      T_tt : (julian_date_tt - JD_J2000)/36525, 
      gmst : greenwich_mean_sidereal_time(julian_date, T_centuries), //0..2pi
      date : date_utc /* back door; eg, in case local values in the current timezone are needed */ 
    };
    result.delta = function(secs){
      var new_msec_epoch = this.msec_epoch + secs*1000;
      var date = new Date(); //any date will do
      date.setTime(new_msec_epoch);
      return when_from_utc(date);
    };
    result.next = function(){
      return this.delta(60*60*24);
    };
    result.prev = function(){
      return this.delta(-1*60*60*24);
    };
    result.startOfDayLT = function(){
      var result = new Date(this.date.getFullYear(), this.date.getMonth(), this.date.getDate());
      return when_from_utc(result);
    };
    result.endOfDayLT = function(){
      var result = new Date(this.date.getFullYear(), this.date.getMonth(), this.date.getDate(), 23, 59, 59, 999);
      return when_from_utc(result);
    };
    //there are various ways to format
    if (text){
      //reuse 'text' if passed to this method
      result.toString = function(){return text;};
    }
    else {
      //otherwise default to UT
      var as_string = when_to_string_ut(result); 
      result.toString = function(){
        return as_string; 
      };
    }
    result.toStringLT = function(weekdays /* array of localized names of days of the week, Sun..Sat */){
      return when_to_string_lt(this, weekdays);
    };
    result.toStringUT = function(weekdays){
      return when_to_string_ut(this, weekdays);
    };
    result.toStringTT = function(weekdays){
      return when_to_string_tt(this, weekdays);
    };
    return result;
  };

  var opt_pad = function(prefix, val, yes_pad){
    var result = '';
    if (val !== undefined){
      result = result + prefix;
      if (yes_pad){
        result = result + pad(val);
      }
      else {
        result = result + val;
      }
    }
    return result;
  };  
  
  /* Format a when in a standard format. */
  var when_to_string_as = function(prefix, when, weekdays){
    var result = prefix + ' ' + opt_pad('', when.y, true) + opt_pad('-', when.m, true) + opt_pad('-', when.d, true);
    result = result + opt_pad(' ', when.hour, true);
    result = result + opt_pad(':', when.min, true);
    result = result + opt_pad(':', when.sec, true);
    result = result + opt_pad('.', when.msec, false);
    if (weekdays !== undefined){ 
      //append the weekday with text from the given array
      result = result + ' ' + weekdays[when.weekday-1]; //1 is Sunday, but arrays are 0-based
    }
    return result;
  };
  
  /* Format a when in the standard way, but using the UT time scale. */
  var when_to_string_ut = function(when, weekdays){
    return when_to_string_as('UT', when, weekdays);
  };
  
  /* Format a when in the standard way, but using the local time zone. */
  var when_to_string_lt = function(when, weekdays){
    return when_to_string_as('LT', pseudo_when_from(when.date), weekdays);
  };
  
  /* This is a 'pseudo-when' because its date is in a local time zone, not UTC. */
  var pseudo_when_from = function(date){
    return {
      y: date.getFullYear(),
      m: date.getMonth() + 1,
      d : date.getDate(),
      hour : date.getHours(),
      min : date.getMinutes(),
      sec : date.getSeconds(),
      msec : date.getMilliseconds(),
      weekday: day_of_the_week(date.getFullYear(), date.getMonth() + 1, date.getDate())
    };
  };
  
  /* Format a when in the standard way, but using the TT time scale. */
  var when_to_string_tt = function(when, weekdays){
    var result = 'TT '; 
    var dt = delta_t(when.y);
    var msecs_tt = when.date.getTime() + dt*1000; // in modern times, TT is slightly ahead of UT
    var pseudo_date = new Date(msecs_tt);
    var pseudo_when = {
      y: pseudo_date.getUTCFullYear(),
      m: pseudo_date.getUTCMonth() + 1,
      d : pseudo_date.getUTCDate(),
      hour : pseudo_date.getUTCHours(),
      min : pseudo_date.getUTCMinutes(),
      sec : pseudo_date.getUTCSeconds(),
      msec : pseudo_date.getUTCMilliseconds(),
      weekday: day_of_the_week(date.getFullYear(), date.getMonth() + 1, date.getDate())
    };
    return when_to_string_as('TT', pseudo_when, weekdays); 
  };

  /* The 'when' derived from the Julian year, eg 2016.5 or similar. */
  var when_from_julian_year = function(julian_year /*eg 2016.5*/, text){
    var msec_j2000 = Date.UTC(2000,0,1,11,58,56,172); //Jan 1, 2000 at 11:58:56.172 UTC; delta-t was 63.8285s at this time
    var msec_per_day = 24*60*60*1000;
    var msec_since_j2000 = (julian_year - 2000) * 365.25 * msec_per_day;
    var date_utc = new Date(msec_j2000 + msec_since_j2000);
    return when_from_utc(date_utc, text);
  };
  
  /* Convenient constant for the most common equinox. */
  var when_j2000 = when_from_julian_year(2000, 'J2000.0');
  
  /* Example: 'J2016.215' */
  var when_parse_julian = function(text){
    var num = parseFloat(text.substring(1)); //chop off the 'J'
    return when_from_julian_year(num, text);
  };
  
  /* 
   Example: 'UT 2016-01-31 02:56:03.123', plus truncations.
   Parse an input string into pieces. Those pieces are suited for building a Javascript Date object. 
   Javascript Date objects know only about UT and LT, not TT. So, when the input has TT, 
   the seconds are tweaked, in order to get the same instant expressed in terms of UT/LT. 
  */
  var when_parse = function(text){
    var original_text = text;
    var style = text.substring(0,2); //first 2 letters
    var text = text.substring(2).trim(); //chop off the 'style'
    var space = text.indexOf(" "); //between date and time
    var dot = text.indexOf("."); //either a fractional day, or a fractional second (but not both)
    var is_fractional_day = (dot === 10); 
    
    //first the date parts only
    var date = (space === -1 ? text : text.substring(0, space));
    var parts = date.split('-');
    var year = parseInt(parts[0], 10);
    var month = parseInt(parts[1], 10);
    var day = is_fractional_day ? Math.floor(parseFloat(parts[2])) : parseInt(parts[2], 10);

    //now for the time parts, which all default to zero, if not present in the input text
    var hour = 0, minute = 0, seconds = 0, msecs = 0; //integers all
    var frac = 0, hour_dec = 0, minute_dec = 0, seconds_dec = 0; // as decimal numbers; used for fractional days
    
    if (! is_fractional_day && text.length > 10){
      //date and time string both present
      var time = text.substring(space+1);
      parts = time.split(':');
      hour = parseInt(parts[0], 10);
      minute = parts.length > 1 ? parseInt(parts[1],10) : 0;
      if (parts.length > 2){
        seconds_dec = parseFloat(parts[2]);
      }
    }
    if (is_fractional_day) {
      frac = parseFloat(parts[2]) - day;
      hour_dec = frac * 24;
      hour = Math.floor(hour_dec);
      minute_dec = (hour_dec - hour) * 60;
      minute = Math.floor(minute_dec);
      seconds_dec = (minute_dec - minute) * 60;
    }

    if ('TT' === style){ //tweak the seconds for ΔT 
      seconds_dec = seconds_dec - delta_t(year); //no longer in 0.0 .. 59.9
      //**according to javascript docs** for the Date object, the runtime will 'odometer' the adjacent values, if 
      //seconds is out of the normal range
    }
    seconds = Math.floor(seconds_dec);
    msecs = Math.round((seconds_dec - seconds)*1000); 
    
    //now we can build a Javascript Date object corresponding to the desired instant
    var result;
    if ('LT' === style){
      result = new Date(year, month-1, day, hour, minute, seconds, msecs); 
    }
    else {
      result = new Date(Date.UTC(year, month-1, day, hour, minute, seconds, msecs));
    }
    return when_from_utc(result, original_text);
  };
  
  /* 
   Example input text: 
   'UT 2016-02-01 13:02:01.123'  (msecs is the finest precision - 3 decimals only).
   'UT 2016-02-01 13:02:01' 
   'LT 2016-02-01 13:02' 
   'UT 2016-02-01 13' 
   'TT 2016-02-01'
   'J2015.0'
   Decimal days are also allowed in the following way:
   'UT 2016-02-01 13.321654...' 
  */
  var when = function(raw_text){
    var text = raw_text.trim().toUpperCase();
    var result;
    if (text.startsWith('J')){
      result = when_parse_julian(text);
    }
    else {
      result = when_parse(text);
    } 
    return result;
  };
  
  var when_now = function(){
    return when_from_utc(new Date());
  };
  
  /*
   Ref: Astronomical Algorithms, Jean Meeus.
   Return .y, .m, .d, with fractional days.
  */
  var find_calendar_date_from_jd = function(jd){
    var temp = jd + 0.5;
    var Z = Math.floor(temp);
    var F = temp - Z;
    var A, alpha, B, C, D, E, y, m, dayFrac;
    if (Z < 2299161){
      A = Z;
    }
    else {
      alpha = Math.floor((Z-1867216.25)/36524.25);
      A = Z + 1 + alpha - Math.floor(alpha/4);
    }
    B = A + 1524;
    C = Math.floor((B-122.1)/365.25);
    D = Math.floor(365.25*C);
    E = Math.floor((B-D)/30.6001);
    dayFrac = B - D - Math.floor(30.6001*E) + F;
    m = E < 14 ? E-1 : E-13;
    y = m > 2 ? C-4716 : C-4715;
    return {
      y:y, 
      m:m,
      d:dayFrac
    };
  };

  var is_leap_year = function(year){
    var result = (year % 4 === 0);
    if (year % 100 === 0){
      //it's a century year; special case
      result = (year % 400 === 0);
    }
    return result;
  };  
  var num_days_in_month = function(year, month){
    var standard = [31,28,31,30,31,30,31,31,30,31,30,31];
    var result = standard[month-1];
    if (month === 2 && is_leap_year(year)){
      result = result + 1;
    }
    return result;
  };
  /*
   Return a modification of the given date-time string, to increment/decrement the given time unit by the 
   given number of steps. When the unit reaches reaches the end of its normal range, it will rollover 
   the next highest unit, like an odometer. 
   date_time_str: '2016-05-22 15:35:20', '2016-05-22 15:20', '2016-05-22 15', '2016-05-22'; ignores msec
   time_unit: 'year, month, day, hour, min, sec'
   num_steps: -9999..+9999
  */
  var date_time_odometer = function(date_time_str, time_unit, num_steps){
    if (Math.abs(num_steps) > 9999){
      var err_message = 'Not allowed to change ' + time_unit + ' by more than 9999 units. Your are trying to change by : ' + num_steps + ' units';
      throw err_message; // early abort
    }
    var when = when_parse('UT ' + date_time_str); //in order to get the parts thereof
    var parts = {
      y: when.y,
      m: when.m,
      d: when.d,
      hour: when.hour,
      min: when.min,
      sec: when.sec
      //ignores msec
    };
    var step = num_steps >=0 ? 1 : -1; 
    var step_year = function(parts){
      parts.y = parts.y + step;
    };
    var step_month = function(parts){
      parts.m = parts.m + step;
      if (parts.m > 12 ){
        parts.m = 1;
        step_year(parts);
      }
      else if (parts.m < 1){
        parts.m = 12;
        step_year(parts);
      }
    };
    var step_day = function(parts){
      parts.d = parts.d + step;
      var days_in_month = num_days_in_month(parts.y, parts.m);
      if (parts.d > days_in_month){
        parts.d = 1;
        step_month(parts);
      }
      else if (parts.d < 1){
        step_month(parts); // first go to the previous month
        parts.d = num_days_in_month(parts.y, parts.m); //the end of the previous month
      }
    };
    var step_hour = function(parts){
      parts.hour = parts.hour + step;
      if (parts.hour > 23 ){
        parts.hour = 0;
        step_day(parts);
      }
      else if (parts.hour < 0){
        parts.hour = 23;
        step_day(parts);
      }
    };
    var step_min = function(parts){
      parts.min = parts.min + step;
      if (parts.min > 59){
        parts.min = 0;
        step_hour(parts);
      }
      else if (parts.min < 0){
        parts.min = 59;
        step_hour(parts);
      }
    };
    var step_sec = function(parts){
      parts.sec = parts.sec + step;
      if (parts.sec > 59){ //ignores/drops fractional seconds 
        parts.sec = 0;
        step_min(parts);
      }
      else if (parts.sec < 0){
        parts.sec = 59;
        step_min(parts);
      }
    };
    var step_func = {
      year: step_year, 
      month: step_month,
      day: step_day,
      hour: step_hour,
      min: step_min,
      sec: step_sec
    };
    if (! step_func[time_unit]) {
       throw new 'Time unit should be one of (year|month|day|hour|min|sec), but yours is: ' + time_unit;
    }
    for (var i = 0; i < Math.abs(num_steps); ++i){
      step_func[time_unit](parts);
    }
    return when_to_string_as('', parts).trim();    
  };
  
  /* 
   Return the elongation in rads between two objects.
   The objects need to have α, δ as properties. Meeus page 105.
  */
  var elongation_between = function(ephem_a, ephem_b){
    var result = Math.sin(ephem_a.δ) * Math.sin(ephem_b.δ) + Math.cos(ephem_a.δ)*Math.cos(ephem_b.δ)*Math.cos(ephem_a.α - ephem_b.α);
    var result = Math.acos(result); //0..pi
    if (Math.abs(result % Math.PI) < rads(0.25)){
      //if near 0 or 180 degs, then use a more accurate formula
      var delta_α = (ephem_a.α - ephem_b.α);
      var delta_δ = (ephem_a.δ - ephem_b.δ);
      result = Math.sqrt(Math.pow(delta_α * Math.cos(ephem_a.δ),2) + Math.pow(delta_δ,2));
    }
    return result; 
  };
  
  /* 
   CAREFUL: this only makes sense when using the EQUINOX OF DATE.
   Adds/resets these ephem props: .a, .A, .h (alt, az, local hour angle).
   If where is topocentric, then a correction for parallax is applied (if the obj is nearby), which decreases the altitude. 
   This is significant for the Moon; for Venus and Mars, the effect is near the limit of this library's resolution.
  */
  var convert_αδ_to_aA = function(ephem, where, when){
    if (ephem.equinox && Math.abs(ephem.equinox.jd - when.jd) > 200){
      console.log("WARNING: when finding a, A (altitude and azimuth) of an object, you should use the mean equinox of date.");
    }
    ephem.h = find_local_hour_angle(ephem, when, where);
    ephem.a = Math.asin(Math.sin(where.φ) * Math.sin(ephem.δ) + Math.cos(where.φ) * Math.cos(ephem.δ) * Math.cos(ephem.h)); //-pi/2 .. + pi/2
    if (where.is_topocentric){
      //apply a correction for parallax (displacement between geocentric and topocentric), if not too far away
      if (ephem.Δ && ephem.Δ < 2.0) { 
        var parallax = Math.asin(4.258750E-5*Math.cos(ephem.a)/ephem.Δ); //-pi/2..pi/2. Meeus p265. Moon in arcmin: 57 avg, range 54..61. 
        ephem.a = ephem.a - parallax;
        /*
        var par = {a: parallax};
        convert_all_angles_to_degs_sexagesimal(par);
        console.log("Parallax: " + par.a.toString());
        */
        //DANGEROUS SIDE EFFECT: this changes the very inputs to this calc; if this func is called again, then result is WRONG
        //An alternative: this func to DO NOTHING if .h, .a, and .A are already set
        //apply_parallax_to_αδ(ephem, where); 
        //apply_parallax_to_λβ(ephem, where, when);
      }
    }
    var numer = Math.sin(ephem.h);
    var denom = Math.cos(ephem.h) * Math.sin(where.φ) - Math.tan(ephem.δ) * Math.cos(where.φ); 
    var az_from_south = Math.atan2(numer, denom); //-pi..+pi
    ephem.A = az_from_south + Math.PI; //0..2pi
  };

  var equatorial_horizontal_parallax = function(ephem){
    return rads(8.794/(3600*ephem.Δ)); 
  };
  
  /* Alters α and δ in place. DON'T call this more than once. */  
  var apply_parallax_to_αδ = function(ephem, where){
    if (ephem.α){
      //Meeus 1991, page 263; spherical Earth
      var pi = equatorial_horizontal_parallax(ephem); 
      var numer = -1*Math.cos(where.φ)*Math.sin(pi)*Math.sin(ephem.h);
      var denom = Math.cos(ephem.δ) - Math.cos(where.φ) * Math.sin(pi) * Math.cos(ephem.h);
      var Δα = Math.atan2(numer, denom); //-pi..+pi
      numer = (Math.sin(ephem.δ) - Math.sin(where.φ) * Math.sin(pi))*Math.cos(Δα);
      denom = Math.cos(ephem.δ) - Math.cos(where.φ) * Math.sin(pi) * Math.cos(ephem.h);
      var δ_prime = Math.atan2(numer, denom); //-pi..+pi
      ephem.δ = δ_prime;
      ephem.α = ephem.α + Δα; 
    }
  };
  
  /* Alters λ and β in place. DON'T call this more than once. */  
  var apply_parallax_to_λβ = function(ephem, where, when){
    if (ephem.λ){
      //Meeus 1991, page 266; spherical Earth
      var S = Math.sin(where.φ);
      var C = Math.cos(where.φ);
      var pi = equatorial_horizontal_parallax(ephem); 
      var λ = ephem.λ;
      var β = ephem.β;
      var θ = find_mean_sidereal_time_at_longitude(when, where); 
      var ε = obliquity_of_ecliptic(when);
      var N = Math.cos(λ)*Math.cos(β) - C * Math.sin(pi)*Math.cos(θ);
      var numer = Math.sin(λ)*Math.cos(β) - Math.sin(pi)*(S*Math.sin(ε) + C*Math.cos(ε)*Math.sin(θ));
      var denom = N;
      var λ_prime = in2pi(Math.atan2(numer, denom));
      numer = Math.cos(λ_prime)*(Math.sin(β) - Math.sin(pi)*(S*Math.cos(ε) - C*Math.sin(ε)*Math.sin(θ)));
      ephem.β = Math.atan(numer/denom); //-pi/2..+pi/2
      ephem.λ = λ_prime; 
    }
  };

  var convert_αδ_to_λβ = function(ephem, when){
    var e = obliquity_of_ecliptic(when);
    var num = Math.sin(ephem.α)*Math.cos(e) + Math.tan(ephem.δ)*Math.sin(e);
    var denom = Math.cos(ephem.α);
    var λ = Math.atan2(num, denom); //-pi..+pi
    var sin_β = Math.sin(ephem.δ) * Math.cos(e) - Math.cos(ephem.δ) * Math.sin(e) * Math.sin(ephem.α);
    ephem.λ = in2pi(λ); //0..2pi
    ephem.β = Math.asin(sin_β); //-pi/2..+pi/2
  };
  
  var convert_λβ_to_αδ = function(ephem, when){
    var e = obliquity_of_ecliptic(when);
    var num = Math.sin(ephem.λ)*Math.cos(e) - Math.tan(ephem.β)*Math.sin(e);
    var denom = Math.cos(ephem.λ);
    var α = Math.atan2(num, denom); //-pi..+pi
    var δ = Math.sin(ephem.β) * Math.cos(e) + Math.cos(ephem.β) * Math.sin(e) * Math.sin(ephem.λ);
    ephem.α = in2pi(α); //0..2pi
    ephem.δ = Math.asin(δ); //-pi/2..+pi/2
  };
  
  var convert_αδ_to_XYZ = function(ephem){
    ephem.X = ephem.Δ*Math.cos(ephem.δ)*Math.cos(ephem.α);
    ephem.Y = ephem.Δ*Math.cos(ephem.δ)*Math.sin(ephem.α);
    ephem.Z = ephem.Δ*Math.sin(ephem.δ);
  };
  
  var convert_XYZ_to_αδ = function(ephem){
    var Δ = Math.sqrt(ephem.X*ephem.X + ephem.Y*ephem.Y + ephem.Z*ephem.Z);
    ephem.α = in2pi(Math.atan2(ephem.Y, ephem.X)); //atan2 is -pi..+pi
    ephem.δ = Math.asin(ephem.Z/Δ); //-pi/2..+pi/2
    ephem.Δ = Δ;
  };
  
  var convert_λβ_to_XYZ = function(ephem, when){
    var e = obliquity_of_ecliptic(when);
    ephem.X = ephem.Δ * Math.cos(ephem.β) * Math.cos(ephem.λ);
    ephem.Y = ephem.Δ * (Math.cos(ephem.β) * Math.sin(ephem.λ) * Math.cos(e) - Math.sin(ephem.β) * Math.sin(e));
    ephem.Z = ephem.Δ * (Math.cos(ephem.β) * Math.sin(ephem.λ) * Math.sin(e) + Math.sin(ephem.β) * Math.cos(e));
  };

  var convert_XYZ_to_λβ = function(ephem){
    convert_XYZ_to_αδ(ephem);
    convert_αδ_to_λβ(ephem);
  };

  /* Convert heliocentric (xyz) to geocentric (XYZ). */  
  var convert_xyz_to_XYZ = function(ephem, sun){
    ephem.X = sun.X + ephem.x;
    ephem.Y = sun.Y + ephem.y;
    ephem.Z = sun.Z + ephem.z;
  };
  
  /* Adds the 'zodiac' property to the ephem. */  
  var convert_ra_to_zodiac_sign = function(ephem){
    var result = zodiac[0]; // default, for safety, to avoid nulls
    for (var idx = 0; idx < zodiac.length; ++idx){
      if (ephem.α < zodiac[idx].α_end){
        result = zodiac[idx];
        break;
      }
    }
    ephem.zodiac = result;
  };
  
  var convert_all_angles_to_degs = function(thing){
    var convert = function(thing, property){
      if (thing[property]){
        thing[property] = degs(thing[property]);
      }
    };
    convert_all_angles(thing, convert);
  };
  
  var convert_all_angles_to_degs_sexagesimal = function(thing){
    var convert = function(thing, property){
      if (thing[property]){
        thing[property] = degs_sexagesimal(thing[property]);
      }
    };
    convert_all_angles(thing, convert);
  };
  
  var convert_all_angles_to_rads = function(thing){
    var convert = function(thing, property){
      if (thing[property]){
        thing[property] = rads(thing[property]);
      }
    };
    convert_all_angles(thing, convert);
  };
  
  /* Apply a conversion function to all angles. */
  var convert_all_angles = function(thing, convert){
    convert(thing, 'α');
    convert(thing, 'δ');
    convert(thing, 'λ');
    convert(thing, 'β');
    convert(thing, 'l');
    convert(thing, 'b');
    convert(thing, 'a');
    convert(thing, 'A');
    convert(thing, 'h');
    
    convert(thing, 'φ');
    //convert(thing, 'λ'); //already done above!
    
    convert(thing, 'i');
    convert(thing, 'Ω');
    convert(thing, 'π');
    convert(thing, 'ω');
    convert(thing, 'L0');
    convert(thing, 'M0');
    convert(thing, 'n');
    convert(thing, 'v');
    
    convert(thing, 'elong');
    //convert(thing, 'size'); // the angular size is left out because it's under an arcmin
    convert(thing, 'phase');
  };
  
  /* 
    Distance in kilometers between two places on the Earth's surface.
    The Earth is modeled here as a simple sphere.
    Params x and y are two objects carrying .φ and .λ , in degrees.
  */
  var distance_kms = function (x, y){
    var a = where(x.φ, x.λ); //converts to rads, for internal calcs
    var b = where(y.φ, y.λ);
    var Δλ = Math.abs(a.λ - b.λ);
    var Δφ = Math.abs(a.φ - b.φ);
    var result = 0.0;
    var value = Math.sin(a.φ) * Math.sin(b.φ) + Math.cos(a.φ) * Math.cos(b.φ) * Math.cos(Δλ);
    if (Math.abs(1 - value) > 0.01) {
      //the long-distance formula can be used; the cosine value is not too close to 1
      result = Math.acos(value); //0..pi rads
    }
    else {
      //the distance is small, and needs a different formula, more accurate at small distances
      value = Math.sqrt(  Math.pow(Math.sin(Δφ/2.0), 2) + Math.cos(a.φ)* Math.cos(b.φ)* Math.pow(Math.sin(Δλ/2.0), 2) );
      result = 2 * Math.abs(Math.asin(value)); // asin is -pi/2..+pi/2 rads, so we need it to be positive here
    }
    result = result * 6371.001; // kilometres
    return result;
  };
  
  /* The angle zenith-object-NCP (north celestial pole), -pi..+pi. Negative before meridian, positive after. Meeus p94. */
  var parallactic_angle = function(ephem, when, where){
    var result = 0;
    if (!ephem.h){
      ephem.h = find_local_hour_angle(ephem, when, where);
    }
    var numer = Math.sin(ephem.h);
    var denom = Math.tan(where.φ) * Math.cos(ephem.δ) - Math.sin(ephem.δ) * Math.cos(ephem.h);
    if (denom !== 0){
      result = Math.atan2(numer, denom); //-pi..pi
    }
    return result;
  };
  
  /* 
   The position angle of the center of the bright limb, 0..2pi. Meeus p316.
   Callers will often want to combine this with a parllactic angle.
  */
  var bright_limb_angle = function(ephem, sun){
    var numer = Math.cos(sun.δ) * Math.sin(sun.α - ephem.α);
    var denom = Math.sin(sun.δ) * Math.cos(ephem.δ) - Math.cos(sun.δ) * Math.sin(ephem.δ) * Math.cos(sun.α - ephem.α);
    var result = Math.atan2(numer, denom); //-pi..pi
    return in2pi(result); //0..2pi
  };
  
  /* Position angle, NCP-from-to. 0..2pi. Measured eastwards. */
  var position_angle_between = function(from, to){
    //this is my own formula; use sine law and cosine law
    var numer = Math.cos(from.δ) * Math.cos(to.δ) * Math.sin(to.α - from.α);
    var denom = Math.sin(to.δ) - Math.cos(elongation_between(from, to)) * Math.sin(from.δ);
    var result = Math.atan2(numer, denom);
    return result; 
  };

  /*
    Returns the geodetic position the geomagnetic North Pole, in an object having .φ, .λ, in degrees. 
    The pole position is important for predicting auroras at a given location. 
    Simple linear extrapolation over a few years from the present.
    British Geological Survey: http://www.geomag.bgs.ac.uk/education/poles.html
  */
  var geomagnetic_north_pole = function(when){
    var num_years = 5;
    var where_pole_2015 = { φ: 80.37, λ: -72.63 };
    var where_pole_2020 = { φ: 80.65, λ: -73.17 };
    var Δφ_annual = (where_pole_2020.φ - where_pole_2015.φ)/num_years;  
    var Δλ_annual = (where_pole_2020.λ - where_pole_2015.λ)/num_years;
    
    var origin = EPH.when("UT 2015-01-01");
    var Δyears = (when.mjd - origin.mjd)/365.25; //close enough
    
    return {
      φ: where_pole_2015.φ + Δφ_annual * Δyears, 
      λ: where_pole_2015.λ + Δλ_annual * Δyears
    }; 
  };
  
  /*
   Where has properties .φ, .λ in degrees, for latitude and longitude. 
   Returns degrees from geomagnetic north pole. Always positive, 0..180. Earth is modeled as a sphere. 
  */
  var geomagnetic_latitude = function(where, when){
    var geomagnetic_pole = geomagnetic_north_pole(when);
    var km = distance_kms(where, geomagnetic_pole);
    var km_per_degree = 111.273; //approximate median value
    var co_latitude = km / km_per_degree;
    return (90 - co_latitude);
  };
  
  /*
   Typical Kp value (auroral activity value) needed before you can see aurora.  
   From: http://www.swpc.noaa.gov/content/tips-viewing-aurora 
  */
  var aurora_min_kp = function(geomagnetic_latitude /*degrees*/){
    var values = [66.5, 64.5, 62.4, 60.4, 58.3, 56.3, 54.2, 52.2, 50.1, 48.1]; //index is Kp 0..9; find the index corresponding to the closest one
    var i, diff = 1000, this_diff, result = 0;
    for (var i = 0; i < values.length; ++i){
      this_diff = Math.abs(geomagnetic_latitude - values[i]);
      if (this_diff < diff){
        diff = this_diff; 
        result = i;
      }
    }
    return result;
  };
  
  /* Overwrites properties of ephem in place: equinox, α, δ. Overwrites (λ,β), (X,Y,Z) as well (if present). */  
  var apply_precession = function(ephem, to_when, precess_angles /*optional*/){
    //do nothing if the difference in equinoxes is less than about half a year. 
    //on that time scale, the amount of precession is on the order of the precision of this library.
    //in addition, it also means that stars and messier objects (in a given calendar year) will not 
    //need to be precessed; this will save significant computation time, I believe, when dealing with 
    //large numbers of objects.
    if (Math.abs(ephem.equinox.jd - to_when.jd) < 200.0) {
      return; 
    }
    
    var from_when = ephem.equinox;    
    var angles = precess_angles;
    if (! precess_angles){
      angles = precession_angles(from_when, to_when);
    }
    var A = Math.cos(ephem.δ) * Math.sin(ephem.α + angles.zeta);
    var B = Math.cos(angles.theta) * Math.cos(ephem.δ) * Math.cos(ephem.α + angles.zeta) - Math.sin(angles.theta) * Math.sin(ephem.δ);
    var C = Math.sin(angles.theta) * Math.cos(ephem.δ) * Math.cos(ephem.α + angles.zeta) + Math.cos(angles.theta) * Math.sin(ephem.δ);
    if (Math.abs(ephem.δ) < rads(85)){
      ephem.δ = Math.asin(C); //-pi/2..+pi/2
    }
    else {
      var temp = Math.acos(Math.sqrt(A*A + B*B)); //0..pi
      ephem.δ = Math.sin(ephem.δ) * temp;
    }
    ephem.α = in2pi(Math.atan2(A, B) + angles.z); // 0..2pi
    
    //keep all the different coord systems in sync, if they are present
    if (ephem.λ){
     convert_αδ_to_λβ(ephem, to_when);
    }
    if (ephem.X){
     convert_λβ_to_XYZ(ephem, to_when);
    }
    ephem.equinox = to_when;
  };
  
  /* Used to transform from one equinox to another. */  
  var precession_angles = function(from_when, to_when){
    var T = from_when.T;
    var t = to_when.T - from_when.T;
    //arcseconds; equatorial coords
    var zeta =  (2306.2181 + 1.39656*T - 0.000139*T*T)*t + (0.30188 - 0.000344*T)*t*t + (0.017998)*t*t*t;
    var z =     (2306.2181 + 1.39656*T - 0.000139*T*T)*t + (1.09468 + 0.000066*T)*t*t + (0.018203)*t*t*t;
    var theta = (2004.3109 - 0.85330*T - 0.000217*T*T)*t - (0.42665 + 0.000217*T)*t*t - (0.041833)*t*t*t;
    var convert_secs = function(arcsec){
      return rads(arcsec/3600);
    };
    return {
     zeta: convert_secs(zeta),
     z: convert_secs(z), 
     theta: convert_secs(theta)
    };
  };
  
  /* Return Δψ, Δε, in rads. The caller decides what to do with it. Low precision, 0.5arcsec. */
  var nutation = function(when){
    var T = when.T;
    var Ω = rads(125.04452 - 1934.136261*T); 
    var L = rads(280.4665 + 36000.7698*T);
    var L_prime = rads(218.3165 + 481267.8813*T);
    var Δψ = -17.20*Math.sin(Ω) - 1.32*Math.sin(2*L) - 0.23*Math.sin(2*L_prime) + 0.21*Math.sin(2*Ω); // arcsec
    var Δε =   9.20*Math.cos(Ω) + 0.57*Math.cos(2*L) + 0.10*Math.cos(2*L_prime) - 0.09*Math.cos(2*Ω); // arcsec
    return {
      Δψ: rads(round(Δψ,1)/3600), 
      Δε: rads(round(Δε,1)/3600)
    };
  };

  var arc_sec_to_rads = function (arc_secs){
    return rads(arc_secs/3600);
  };
  var aberration_memo = []; //N obj's that contain .when, and constants related to the Sun, at a given time
  /* Return Δλ, Δβ, Δα, Δδ, in rads. The caller decides what to do with it.*/
  var annual_aberration = function(ephem, when){
    var i, day_nums, e, pi, θ, ε, T, a, b, λ, β, α, δ;
    var κ = 20.49552; //constant of aberration, arcsecs
    for(i=0; i < aberration_memo.length; i++){
      if (aberration_memo[i].jd === when.jd){
        day_nums = aberration_memo[i]; 
      }
    }
    if (!day_nums){
      T = when.T;
      e = 0.016708617 - 0.000042037*T - 0.0000001236*T*T; //eccentricity
      pi = rads(102.93735 + 1.71953*T + 0.00046*T*T); //longitude of perihelion
      θ = position_sun(when).λ;
      ε = obliquity_of_ecliptic(when);
      day_nums = {jd: when.jd,  e: e, pi: pi, θ: θ, ε: ε};
      aberration_memo.push(day_nums);
    }
    else {
      //keeping these short helps read the formulas below
      e = day_nums.e;  pi = day_nums.pi;  θ = day_nums.θ;  ε = day_nums.ε;
    }
    λ = ephem.λ; β=ephem.β; α=ephem.α; δ=ephem.δ; //just to make the long formulae more legible
    var result = {};
    if (ephem.λ){
      result.Δλ = (-κ*Math.cos(θ - λ) + e*κ*Math.cos(pi - λ))/Math.cos(β),
      result.Δβ = -κ*Math.sin(β) * (Math.sin(θ - λ) - e*Math.sin(pi - λ))
      result.Δλ = arc_sec_to_rads(result.Δλ);  
      result.Δβ = arc_sec_to_rads(result.Δβ);  
    }
    if (ephem.α){
      a = Math.cos(α)*Math.cos(θ)*Math.cos(ε) + Math.sin(α)*Math.sin(θ);
      b = Math.cos(α)*Math.cos(pi)*Math.cos(ε) + Math.sin(α)*Math.sin(pi);
      result.Δα = (-κ*a/Math.cos(δ)) + (e*κ*b/Math.cos(δ));
      a = Math.cos(θ)*Math.cos(ε)*(Math.tan(ε)*Math.cos(δ)  - Math.sin(α)*Math.sin(δ)) + Math.cos(α)*Math.sin(δ)*Math.sin(θ);
      b = Math.cos(pi)*Math.cos(ε)*(Math.tan(ε)*Math.cos(δ) - Math.sin(α)*Math.sin(δ)) + Math.cos(α)*Math.sin(δ)*Math.sin(pi);
      result.Δδ = -κ*a + e*κ*b;
      result.Δα = arc_sec_to_rads(result.Δα);  
      result.Δδ = arc_sec_to_rads(result.Δδ);  
    }
    return result;
  };
  
  /* Overwrite the alt property in place. Assumes rads on input! */
  var add_refraction_to_alt = function(ephem /*rads*/){
    if (ephem.a){
      //alt in degrees is needed by the formula 
      var a_degs = degs(ephem.a);
      var bottom = a_degs + (10.3/(a_degs+5.11)); //degs
      var denom = Math.tan(rads(bottom));
      var refraction_arcmin = 1.02/denom;
      var refraction_rads = rads(refraction_arcmin/60);
      ephem.a = ephem.a + refraction_rads;
    }
  };
  
  /* 
   Geometric position of the Sun for the mean equinox of date.
   Accuracy 0.01 degrees.
   Returns an object having these properties (all angles in rads):
     .λ - ecliptic longitude 
     .β - ecliptic latitude
     .Δ - distance (AU)
     .α - right ascension
     .δ - declination
   Ref: Astronomical Algorithms, Jean Meeus, Chapter 24, low precision 0.01 degrees.
  */
  var position_sun = function(when){
    var T = when.T_tt;
    var L = rads(280.46645 + 36000.76983*T + 0.0003032*T*T); //geometric mean longitude
    var M = rads(357.52910 + 35999.05030*T - 0.0001559*T*T - 0.00000048*T*T*T); //mean anomaly 
    var e = 0.016708617 - 0.000042037*T - 0.0000001236*T*T; //eccentricity of the Earth's orbit; dimensionless
    var C = (1.914600 - 0.004817*T - 0.000014*T*T) * Math.sin(M) + 
            (0.019993 - 0.000101*T) * Math.sin(2*M) + 
             0.000290 * Math.sin(3*M); //equation of center, in deg
    C = rads(C);
    var theta = L + C; //longitude
    var v = M + C; //true anomaly
    var num = 1.000001018 * (1 - e*e);
    var denom = 1 + e * Math.cos(v);
    var R = num/denom; //distance in AU
    var ephem = {
      when: when,
      equinox: when,
      λ: in2pi(theta), 
      β: 0,
      Δ: R
    };
    convert_λβ_to_αδ(ephem, when); 
    convert_λβ_to_XYZ(ephem, when); 
    convert_ra_to_zodiac_sign(ephem); 
    return ephem;
  };
  
  /* Return P, B in rads. A.A., Meeus, chapter 28, page 177, but without L0. */
  var physical_sun = function(when, sun){
    var jd = when.jd_tt;
    var θ = rads((jd - 2398220) * 360/35.38);
    var I = rads(7.25);
    var K = rads(73.6667 + 1.3958333 * (jd - 2396758) / 36525);
    var ε = obliquity_of_ecliptic(when);
    var λ = sun.λ + annual_aberration(sun, when).Δλ;
    var λ_prime = λ + nutation(when).Δψ;
    var x = Math.atan(-1*Math.cos(λ_prime) * Math.tan(ε)) // -pi/2..+pi/2
    var y = Math.atan(-1*Math.cos(λ - K) * Math.tan(I)); // -pi/2..+pi/2
    var P = x + y; //-26..+26, in practice
    var B = Math.asin(Math.sin(λ - K) * Math.sin(I)); // -pi/2..+pi/2
    // η ? how to treat the quadrant??
    return { P: P, B: B };
  };
  
  /* Meeus page 307. Mean equinox of date. */
  var position_moon = function(when){
    var T = when.T;
    //degrees (converted to rads below)
    var L1 = 218.3164591 + T * (481267.88134236 + T * (-0.0013268 + T * ((1 / 538841))));    
    var D = 297.8502042 + T * (445267.1115168 + T * (-0.00163 + T * ((1 / 545868))));        
    var M = 357.5291092 + T * (35999.0502909 + T * (-0.0001536 + T * ((1 / 24490000))));     
    var M1 = 134.9634114 + T * (477198.86763133 + T * (0.008997 + T * ((1 / 69699))));       
    var F = 93.2720993 + T * (483202.0175273 + T * (-0.0034029 - T * ((1 / 3526000))));
    var a1 = 119.75 + 131.849 * T;
    var A2 = 53.09 + 479264.29 * T;
    var A3 = 313.45 + 481266.484 * T;
    //dimensionless   
    var E = 1 + T * (-0.002516 - T * (-0.0000074));
    //convert to rads  
    L1 = rads(in360(L1));
    D = rads(in360(D));
    M = rads(in360(M));
    M1 = rads(in360(M1));
    F = rads(in360(F));
    a1 = rads(in360(a1));
    A2 = rads(in360(A2));
    A3 = rads(in360(A3));
    var calc_longitude = function(){
      var sumL = 0; //unit is 0.000001 degreees, initially
      sumL = sumL + 1 * 1 * 6288774 * Math.sin(0 + 0 * D + 0 * M + 1 * M1 + 0 * F); 
      sumL = sumL + 1 * 1 * 1274027 * Math.sin(0 + 2 * D + 0 * M - 1 * M1 + 0 * F);
      sumL = sumL + 1 * 1 * 658314 * Math.sin(0 + 2 * D + 0 * M + 0 * M1 + 0 * F);
      sumL = sumL + 1 * 1 * 213618 * Math.sin(0 + 0 * D + 0 * M + 2 * M1 + 0 * F);
      sumL = sumL - E * 1 * 185116 * Math.sin(0 + 0 * D + 1 * M + 0 * M1 + 0 * F);
      sumL = sumL - 1 * 1 * 114332 * Math.sin(0 + 0 * D + 0 * M + 0 * M1 + 2 * F);
      sumL = sumL + 1 * 1 * 58793 * Math.sin(0 + 2 * D + 0 * M - 2 * M1 + 0 * F);
      sumL = sumL + E * 1 * 57066 * Math.sin(0 + 2 * D - 1 * M - 1 * M1 + 0 * F);
      sumL = sumL + 1 * 1 * 53322 * Math.sin(0 + 2 * D + 0 * M + 1 * M1 + 0 * F);
      sumL = sumL + E * 1 * 45758 * Math.sin(0 + 2 * D - 1 * M + 0 * M1 + 0 * F);   
      sumL = sumL - E * 1 * 40923 * Math.sin(0 + 0 * D + 1 * M - 1 * M1 + 0 * F);
      sumL = sumL - 1 * 1 * 34720 * Math.sin(0 + 1 * D + 0 * M + 0 * M1 + 0 * F);
      sumL = sumL - E * 1 * 30383 * Math.sin(0 + 0 * D + 1 * M + 1 * M1 + 0 * F);
      sumL = sumL + 1 * 1 * 15327 * Math.sin(0 + 2 * D + 0 * M + 0 * M1 - 2 * F);
      sumL = sumL - 1 * 1 * 12528 * Math.sin(0 + 0 * D + 0 * M + 1 * M1 + 2 * F);
      sumL = sumL + 1 * 1 * 10980 * Math.sin(0 + 0 * D + 0 * M + 1 * M1 - 2 * F);
      sumL = sumL + 1 * 1 * 10675 * Math.sin(0 + 4 * D + 0 * M - 1 * M1 + 0 * F);
      sumL = sumL + 1 * 1 * 10034 * Math.sin(0 + 0 * D + 0 * M + 3 * M1 + 0 * F);
      sumL = sumL + 1 * 1 * 8548 * Math.sin(0 + 4 * D + 0 * M - 2 * M1 + 0 * F);
      sumL = sumL - E * 1 * 7888 * Math.sin(0 + 2 * D + 1 * M - 1 * M1 + 0 * F);    
      sumL = sumL - E * 1 * 6766 * Math.sin(0 + 2 * D + 1 * M + 0 * M1 + 0 * F);
      sumL = sumL - 1 * 1 * 5163 * Math.sin(0 + 1 * D + 0 * M - 1 * M1 + 0 * F);
      sumL = sumL + E * 1 * 4987 * Math.sin(0 + 1 * D + 1 * M + 0 * M1 + 0 * F);
      sumL = sumL + E * 1 * 4036 * Math.sin(0 + 2 * D - 1 * M + 1 * M1 + 0 * F);
      sumL = sumL + 1 * 1 * 3994 * Math.sin(0 + 2 * D + 0 * M + 2 * M1 + 0 * F);
      sumL = sumL + 1 * 1 * 3861 * Math.sin(0 + 4 * D + 0 * M + 0 * M1 + 0 * F);
      sumL = sumL + 1 * 1 * 3665 * Math.sin(0 + 2 * D + 0 * M - 3 * M1 + 0 * F);
      sumL = sumL - E * 1 * 2689 * Math.sin(0 + 0 * D + 1 * M - 2 * M1 + 0 * F);
      sumL = sumL - 1 * 1 * 2602 * Math.sin(0 + 2 * D + 0 * M - 1 * M1 + 2 * F);
      sumL = sumL + E * 1 * 2390 * Math.sin(0 + 2 * D - 1 * M - 2 * M1 + 0 * F);    
      //add more terms here, if needed
      //add corrections due to a1,A2,A3
      sumL = sumL + 3958 * Math.sin(a1) + 1962 * Math.sin(L1 - F) + 318 * Math.sin(A2);
      var result = L1 + rads(sumL/1000000);
      result = in2pi(result);
      return result;        
    };
    var calc_distance = function(){
      //in Meeus the terms are not in order of decreasing amplitude. The largest unretained amplitude
      //is 5.751 km in the following:
      var sumR = 0; // unit is 0.001 km
      sumR = sumR - 1 * 1 * 20905355 * Math.cos(0 + 0 * D + 0 * M + 1 * M1 + 0 * F);
      sumR = sumR - 1 * 1 * 3699111 * Math.cos(0 + 2 * D + 0 * M - 1 * M1 + 0 * F);
      sumR = sumR - 1 * 1 * 2955968 * Math.cos(0 + 2 * D + 0 * M + 0 * M1 + 0 * F);
      sumR = sumR - 1 * 1 * 569925 * Math.cos(0 + 0 * D + 0 * M + 2 * M1 + 0 * F);
      sumR = sumR + E * 1 * 48888 * Math.cos(0 + 0 * D + 1 * M + 0 * M1 + 0 * F);
      sumR = sumR - 1 * 1 * 3149 * Math.cos(0 + 0 * D + 0 * M + 0 * M1 + 2 * F);
      sumR = sumR + 1 * 1 * 246158 * Math.cos(0 + 2 * D + 0 * M - 2 * M1 + 0 * F);
      sumR = sumR - E * 1 * 152138 * Math.cos(0 + 2 * D - 1 * M - 1 * M1 + 0 * F);
      sumR = sumR - 1 * 1 * 170733 * Math.cos(0 + 2 * D + 0 * M + 1 * M1 + 0 * F);
      sumR = sumR - E * 1 * 204586 * Math.cos(0 + 2 * D - 1 * M + 0 * M1 + 0 * F);  
      sumR = sumR - E * 1 * 129620 * Math.cos(0 + 0 * D + 1 * M - 1 * M1 + 0 * F);
      sumR = sumR + 1 * 1 * 108743 * Math.cos(0 + 1 * D + 0 * M + 0 * M1 + 0 * F);
      sumR = sumR + E * 1 * 104755 * Math.cos(0 + 0 * D + 1 * M + 1 * M1 + 0 * F);
      sumR = sumR + 1 * 1 * 10321 * Math.cos(0 + 2 * D + 0 * M + 0 * M1 - 2 * F);
      sumR = sumR + 1 * 1 * 79661 * Math.cos(0 + 0 * D + 0 * M + 1 * M1 - 2 * F);
      sumR = sumR - 1 * 1 * 34782 * Math.cos(0 + 4 * D + 0 * M - 1 * M1 + 0 * F);
      sumR = sumR - 1 * 1 * 23210 * Math.cos(0 + 0 * D + 0 * M + 3 * M1 + 0 * F);
      sumR = sumR - 1 * 1 * 21636 * Math.cos(0 + 4 * D + 0 * M - 2 * M1 + 0 * F);
      sumR = sumR + E * 1 * 24208 * Math.cos(0 + 2 * D + 1 * M - 1 * M1 + 0 * F);
      sumR = sumR + E * 1 * 30284 * Math.cos(0 + 2 * D + 1 * M + 0 * M1 + 0 * F);   
      sumR = sumR - 1 * 1 * 8379 * Math.cos(0 + 1 * D + 0 * M - 1 * M1 + 0 * F);
      sumR = sumR - E * 1 * 16675 * Math.cos(0 + 1 * D + 1 * M + 0 * M1 + 0 * F);
      sumR = sumR - E * 1 * 12831 * Math.cos(0 + 2 * D - 1 * M + 1 * M1 + 0 * F);
      sumR = sumR - 1 * 1 * 10445 * Math.cos(0 + 2 * D + 0 * M + 2 * M1 + 0 * F);
      sumR = sumR - 1 * 1 * 11650 * Math.cos(0 + 4 * D + 0 * M + 0 * M1 + 0 * F);
      sumR = sumR + 1 * 1 * 14403 * Math.cos(0 + 2 * D + 0 * M - 3 * M1 + 0 * F);
      sumR = sumR - E * 1 * 7003 * Math.cos(0 + 0 * D + 1 * M - 2 * M1 + 0 * F);
      sumR = sumR + E * 1 * 10056 * Math.cos(0 + 2 * D - 1 * M - 2 * M1 + 0 * F);
      sumR = sumR + 1 * 1 * 6322 * Math.cos(0 + 1 * D + 0 * M + 1 * M1 + 0 * F);
      sumR = sumR - E * E * 9884 * Math.cos(0 + 2 * D - 2 * M + 0 * M1 + 0 * F);
      sumR = sumR + 1 * 1 * 8752 * Math.cos(0 + 2 * D + 0 * M - 1 * M1 - 2 * F);
      //maintenance HERE if more terms needed
      //no corrections due to A1,A2,A3 are necessary here
      var result = 385000.56 + (sumR/1000);  //km
      result = result / 1.495978707E+8; //AU
      return result;
    };
    var calc_latitude = function(){
      var sumB = 0; // unit is 0.000001 deg
      sumB = sumB + 1 * 1 * 5128122 * Math.sin(0 + 0 * D + 0 * M + 0 * M1 + 1 * F); 
      sumB = sumB + 1 * 1 * 280602 * Math.sin(0 + 0 * D + 0 * M + 1 * M1 + 1 * F);
      sumB = sumB + 1 * 1 * 277693 * Math.sin(0 + 0 * D + 0 * M + 1 * M1 - 1 * F);
      sumB = sumB + 1 * 1 * 173237 * Math.sin(0 + 2 * D + 0 * M + 0 * M1 - 1 * F);
      sumB = sumB + 1 * 1 * 55413 * Math.sin(0 + 2 * D + 0 * M - 1 * M1 + 1 * F);
      sumB = sumB + 1 * 1 * 46271 * Math.sin(0 + 2 * D + 0 * M - 1 * M1 - 1 * F);
      sumB = sumB + 1 * 1 * 32573 * Math.sin(0 + 2 * D + 0 * M + 0 * M1 + 1 * F);
      sumB = sumB + 1 * 1 * 17198 * Math.sin(0 + 0 * D + 0 * M + 2 * M1 + 1 * F);
      sumB = sumB + 1 * 1 * 9266 * Math.sin(0 + 2 * D + 0 * M + 1 * M1 - 1 * F);
      sumB = sumB + 1 * 1 * 8822 * Math.sin(0 + 0 * D + 0 * M + 2 * M1 - 1 * F);   
      sumB = sumB + E * 1 * 8216 * Math.sin(0 + 2 * D - 1 * M + 0 * M1 - 1 * F);
      sumB = sumB + 1 * 1 * 4324 * Math.sin(0 + 2 * D + 0 * M - 2 * M1 - 1 * F);
      sumB = sumB + 1 * 1 * 4200 * Math.sin(0 + 2 * D + 0 * M + 1 * M1 + 1 * F);
      sumB = sumB - E * 1 * 3359 * Math.sin(0 + 2 * D + 1 * M + 0 * M1 - 1 * F);
      sumB = sumB + E * 1 * 2463 * Math.sin(0 + 2 * D - 1 * M - 1 * M1 + 1 * F);
      sumB = sumB + E * 1 * 2211 * Math.sin(0 + 2 * D - 1 * M + 0 * M1 + 1 * F);
      sumB = sumB + E * 1 * 2065 * Math.sin(0 + 2 * D - 1 * M - 1 * M1 - 1 * F);
      sumB = sumB - E * 1 * 1870 * Math.sin(0 + 0 * D + 1 * M - 1 * M1 - 1 * F);
      sumB = sumB + 1 * 1 * 1828 * Math.sin(0 + 4 * D + 0 * M - 1 * M1 - 1 * F);
      sumB = sumB - E * 1 * 1794 * Math.sin(0 + 0 * D + 1 * M + 0 * M1 + 1 * F);   
      sumB = sumB - 1 * 1 * 1749 * Math.sin(0 + 0 * D + 0 * M + 0 * M1 + 3 * F);
      sumB = sumB - E * 1 * 1565 * Math.sin(0 + 0 * D + 1 * M - 1 * M1 + 1 * F);
      sumB = sumB - 1 * 1 * 1491 * Math.sin(0 + 1 * D + 0 * M + 0 * M1 + 1 * F);
      sumB = sumB - E * 1 * 1475 * Math.sin(0 + 0 * D + 1 * M + 1 * M1 + 1 * F);
      sumB = sumB - E * 1 * 1410 * Math.sin(0 + 0 * D + 1 * M + 1 * M1 - 1 * F);
      sumB = sumB - E * 1 * 1344 * Math.sin(0 + 0 * D + 1 * M + 0 * M1 - 1 * F);
      sumB = sumB - 1 * 1 * 1335 * Math.sin(0 + 1 * D + 0 * M + 0 * M1 - 1 * F);
      sumB = sumB + 1 * 1 * 1107 * Math.sin(0 + 0 * D + 0 * M + 3 * M1 + 1 * F);
      sumB = sumB + 1 * 1 * 1021 * Math.sin(0 + 4 * D + 0 * M + 0 * M1 - 1 * F);
      sumB = sumB + 1 * 1 * 833 * Math.sin(0 + 4 * D + 0 * M - 1 * M1 + 1 * F);    
      //add corrections due to A1,A2,A3
      sumB = sumB - 2235 * Math.sin(L1) + 382 * Math.sin(A3) + 175 * Math.sin(a1 - F) + 
                     175 * Math.sin(a1 + F) + 127 * Math.sin(L1 - M1) - 115 * Math.sin(L1 + M1);
      var result = rads(sumB/1000000); //-pi/2..+pi/2
      return result;  
    };
    var ephem = {
      when: when,
      equinox: when, 
      λ: calc_longitude(), 
      Δ: calc_distance(),
      β: calc_latitude()
    };
    convert_λβ_to_αδ(ephem, when); 
    convert_λβ_to_XYZ(ephem, when); 
    convert_ra_to_zodiac_sign(ephem); 
    return ephem;
  }; //end of position_moon
  
  /* 
   Different sources can specify orbits in different ways.
   Core data that's assumed to be always present: 
    equinox,epoch,a,e,i,Ω
   Plus:
    at least 1 of these 2 be present: 
      longitude of perihelion (small-omega-bar, or pi for Meeus)
      argument of perihelion (small-omega)
    at least 1 of these 3 be present; they all serve one purpose: to let you calc the *current* mean anomaly,
    as the first step in finding the object's place in its orbit 
      L0 : mean longitude at epoch
      M0: mean anomaly at epoch
      T : time of perihelion passage
   n: will always be added, if absent
   It's important to keep straight the time to which a mean anomaly M applies: 
     M : current when
     M0: epoch of the orbit
  */
  var add_derived_orbital_items_to = function(orbit, when){
     if (!orbit.n && orbit.a){
       orbit.n = rads(0.9856076686/(orbit.a * Math.sqrt(orbit.a))); //mean motion, rads per day
       orbit.P = 2*Math.PI/orbit.n;
     }
     //1 of these 2 must be present     
     if (!orbit.ω && orbit.π){
       orbit.ω = in2pi(orbit.π - orbit.Ω);
     }
     if (!orbit.π && orbit.ω){
       orbit.pi = in2pi(orbit.ω + orbit.Ω);
     }
     if (!orbit.q && orbit.a){ 
       orbit.q = orbit.a * (1 - orbit.e);
     }
  };

  /*  
   The current mean anomaly is not an orbital param; it's an intermediate param, calculated from the orbit.
   This is the first step in finding the position of the object in its orbit.  
  */  
  var current_mean_anomaly = function(orbit, when){
    var result = 0;
    if (orbit.M0){
      //Minor Planet center - the simplest style
      result = orbit.M0 + orbit.n * (when.jd - orbit.epoch.jd); 
    }
    else if (orbit.T){
      //Meeus example of Enke orbit
      result = 0 + orbit.n * (when.jd - orbit.T.jd);
    }
    else if (orbit.L0){
      //Meeus' Mercury; Observer's Handbook; JPL low-res 
      result = (orbit.L0 - orbit.π) + orbit.n * (when.jd - orbit.epoch.jd);
    }
    return result;
  };
  
  var intermediate_orbit_params = function(orbit, when){
    var eps = obliquity_of_ecliptic(orbit.equinox);
    var sin_eps = Math.sin(eps);
    var cos_eps = Math.cos(eps);
    var F = Math.cos(orbit.Ω);
    var G = Math.sin(orbit.Ω) * cos_eps;
    var H = Math.sin(orbit.Ω) * sin_eps;
    var P = - Math.sin(orbit.Ω) * Math.cos(orbit.i);
    var Q = Math.cos(orbit.Ω) * Math.cos(orbit.i) * cos_eps - Math.sin(orbit.i) * sin_eps;
    var R = Math.cos(orbit.Ω) * Math.cos(orbit.i) * sin_eps + Math.sin(orbit.i) * cos_eps;
    return {
      a: Math.sqrt(F*F + P*P),
      b: Math.sqrt(G*G + Q*Q),
      c: Math.sqrt(H*H + R*R),
      A: Math.atan2(F,P),
      B: Math.atan2(G,Q),
      C: Math.atan2(H,R)
    };
  };

  /* Solve Kepler's equation. */  
  var find_eccentric_anomaly = function(orbit, when){
    var M = current_mean_anomaly(orbit, when);      
    var E = initial_guess_eccentric_anomaly(orbit, M);
    var small_change = rads(0.000001);
    var change = 10; //any big number will do 
    while (Math.abs(change) > small_change){
      change = (M + orbit.e*Math.sin(E) - E)/(1 - orbit.e*Math.cos(E));
      E = E + change;
    }
    return E;
  };
  
  var initial_guess_eccentric_anomaly = function(orbit, M){
    var initial_guess = M;
    var high_eccentricity = (orbit.e > 0.975) && (orbit.e < 1);
    if (high_eccentricity && Math.abs(M) < rads(30)){
      var a = (1-orbit.e)/(4*orbit.e + 0.5);
      var b = M/(8*orbit.e + 1);
      var c = Math.sqrt(b*b + a*a*a);
      var d = Math.sign(b);
      var e = b + d*c;
      if (e >= 0){
        var z = Math.pow(e, 1/3);
        var s_0 = z - a/2;
        var s = s0 - (0.078*Math(s_0,5))/(1+orbit.e);
        initial_guess = M + orbit.e(3*s -4*s*s*s);
      }
    }
    return initial_guess;
  };

  /* Returning the start of a ephem object, to which other coords will be added. */  
  var find_position_in_orbit = function(orbit, when){
    var result;
    if (orbit.e <= 0.99){
      result = find_position_in_orbit_elliptical(orbit, when);
    }
    else {
      //all other cases are approximated as being parabolic
      //the comet with the greatest eccentricity is currently C/1980 E1 (Bowell): 1.057
      result = find_position_in_orbit_parabolic(orbit, when);
    }
    return result;
  };
  
  /* Returns an object with r,v.  */
  var find_position_in_orbit_elliptical = function(orbit, when){
    var E = find_eccentric_anomaly(orbit, when);
    var r = orbit.a*(1 - orbit.e*Math.cos(E));
    var numer = Math.sqrt(1 + orbit.e) * Math.sin(E/2);
    var denom = Math.sqrt(1 - orbit.e) * Math.cos(E/2);
    var v = in2pi(2 * Math.atan2(numer, denom)); //atan2 -pi..+pi
    var starting_ephem = {
      r : r,
      v : v
    };
    return starting_ephem;
  };
  
  /* Returns an object with r,v. Orbit must have q and T (eg, Minor Planet Center, Observer's Handbook. */
  var find_position_in_orbit_parabolic = function(orbit, when){
    //Meeus, page 225
    var t_minus_T = when.jd - orbit.T.jd;
    var W = 0.03649116245 * t_minus_T/(orbit.q * Math.sqrt(orbit.q));
    var G = W/2;
    var Y = Math.pow(G + Math.sqrt(G*G + 1), 1/3);
    var s = Y - 1/Y;
    var v = 2 * Math.atan(s); //-pi..+pi
    var r = orbit.q * (1 + s*s);
    var starting_ephem = {
      r : r,
      v : v
    };
    return starting_ephem;
  };
  
  /* Adds xyz to the ephem. */
  var find_heliocentric_xyz = function(ephem, consts, orbit){
    //these are equatorial rectangular coords, not ecliptical!
    ephem.x = ephem.r * consts.a * Math.sin(consts.A + orbit.ω + ephem.v);
    ephem.y = ephem.r * consts.b * Math.sin(consts.B + orbit.ω + ephem.v);
    ephem.z = ephem.r * consts.c * Math.sin(consts.C + orbit.ω + ephem.v);
  };

  /* Diameter of an object in rads. */  
  var size_rads = function(standard_semi_diam_arcsecs, dist){
    return rads(2*standard_semi_diam_arcsecs/(3600*dist)); 
  };
  
  var sun = {
    name: 'Sun',
    symbol: '☉',
    position: function(when){
      var T = when.T_tt;
      var L = rads(280.46645 + 36000.76983*T + 0.0003032*T*T); //geometric mean longitude
      var M = rads(357.52910 + 35999.05030*T - 0.0001559*T*T - 0.00000048*T*T*T); //mean anomaly 
      var e = 0.016708617 - 0.000042037*T - 0.0000001236*T*T; //eccentricity of the Earth's orbit; dimensionless
      var C = (1.914600 - 0.004817*T - 0.000014*T*T) * Math.sin(M) + 
              (0.019993 - 0.000101*T) * Math.sin(2*M) + 
               0.000290 * Math.sin(3*M); //equation of center, in deg
      C = rads(C);
      var theta = L + C; //longitude
      var v = M + C; //true anomaly
      var num = 1.000001018 * (1 - e*e);
      var denom = 1 + e * Math.cos(v);
      var R = num/denom; //distance in AU
      var ephem = {
        equinox: when,
        λ: in2pi(theta), 
        β: 0,
        Δ: R
      };
      convert_λβ_to_αδ(ephem, when); 
      convert_λβ_to_XYZ(ephem, when); 
      convert_ra_to_zodiac_sign(ephem); 
      return ephem;
    },
    add_physical: function(ephem){
      ephem.size = size_rads(959.63, ephem.Δ);
      ephem.mag = -26.75; //RASC Observer's Handbook 
    },
    ephem: function(when){
      var ephem = this.position(when);
      this.add_physical(ephem);
      return ephem;
    }
  };
  
  var moon = {
    name: 'Moon',
    symbol: '☽',
    position: function(when){
      return position_moon(when);
    },
    add_physical: function(ephem){
      //Meeus page 360
      ephem.size = 2*1.161729E-5/ephem.Δ; // the distance is in AU, not km here
      //ephem.mag = -12.7; //RASC Observer's Handbook
      //http://astronomy.stackexchange.com/questions/10246/is-there-a-simple-analytical-formula-for-the-lunar-phase-brightness-curve
      ephem.mag = -12.73 + 1.49 * Math.abs(ephem.phase) + 0.043 * Math.pow(ephem.phase, 4); 
    },
    ephem: function(when){
      var ephem = this.position(when);
      var sun = position('sun', when);
      add_physical_ephem(ephem, sun);
      this.add_physical(ephem);
      return ephem;
    }
  };
  
  
  /*
   Observer's Handbook, and Astronomical Almanac: to achieve near arc-sec accuracy for planets, for a 
   specific calendar year, interpolate the elements of two osculating orbits given for different times in the year.
   Careful: for fast movers, L0 can go around N times between a and b. That needs to be factored in, since 
   the overall rate of change is calculated.
  */
  var current_osculating_orbit = function(osc_a, osc_b, when){
    if (osc_a.equinox.jd !== osc_b.equinox.jd){
      console.log("Error. Cannon interpolate two orbits having different equinoxes.");
      return null;
    }
    var days_a_to_b = osc_b.epoch.jd - osc_a.epoch.jd;
    var days_a_to_t = when.jd - osc_a.epoch.jd;
    var new_orbital_element = function(name){
       var rate_per_day = (osc_b[name] - osc_a[name])/days_a_to_b; // + -, angle, a, or e
       return osc_a[name] + rate_per_day * days_a_to_t;
    };
    var result = {
      equinox: osc_a.equinox, 
      epoch: when,
      a: new_orbital_element('a'),
      e: new_orbital_element('e'),
      i: new_orbital_element('i'),
      Ω: new_orbital_element('Ω'),
      π: new_orbital_element('π'),
      L0: in2pi(new_orbital_element('L0')) 
    };
    return result;
  };
  var build_osculating_orbit = function(equinox, epoch, a, e, i, Ω, π, L0){
    return {
      equinox: equinox, epoch: epoch,
      a: a, e: e, i: rads(i), Ω: rads(Ω),  π: rads(π), L0: rads(L0)
    };
  };
  
  var planet_orbit_start = when('UT 2019-03-18'); //Observer's Handbook p 23; usually a 240d interval
  var planet_orbit_end = when('UT 2019-11-13');
  
  var mercury = {
    name: 'Mercury',
    symbol: '☿',
    orbit: function(when){
      var osc_a = build_osculating_orbit(when_j2000, planet_orbit_start, 0.387098, 0.205650, 7.0039, 48.3075, 77.4863, 162.0537);
      var osc_b = build_osculating_orbit(when_j2000, planet_orbit_end,   0.387098, 0.205651, 7.0038, 48.3066, 77.4892, (360*3)+64.2144); //4.0914 deg per d; 982 deg per 240d
      var result = current_osculating_orbit(osc_a, osc_b, when);
      return result;
    },
    add_physical: function(ephem){
      var i = degs(ephem.phase);
      ephem.mag = -0.42 + 5 * Math.log10(ephem.r * ephem.Δ) + 0.0380*i - 0.000273*i*i + 0.000002*i*i*i;
      ephem.size = size_rads(3.36, ephem.Δ); 
    },
    ephem: function(when){
      var ephem = position_from_orbit(this.orbit(when), when);
      this.add_physical(ephem);
      return ephem;
    }
  };
  
  var venus = {
    name: 'Venus',
    symbol: '♀',
    orbit: function(when){
      var osc_a = build_osculating_orbit(when_j2000, planet_orbit_start, 0.723329, 0.006734, 3.3946, 76.6249, 131.5550,        261.7249);
      var osc_b = build_osculating_orbit(when_j2000, planet_orbit_end,   0.723331, 0.006733, 3.3946, 76.6247, 131.5080, (360*1)+286.2376); //1.60212 deg per d; 384 deg per 240d
      
      //USING MEEUS' EXPRESSION FOR L AS A BETTER APPROXIMATION
      //THIS IS STILL MISLEADING DATA SINCE ITS THE SAME NUMBERS AS LAST YEAR, OTHER THAN L
      //var osc_a = build_osculating_orbit(when_j2000, planet_orbit_start, 0.723331, 0.006746, 3.3944, 76.6331, 131.6533,        340.87528);
      //var osc_b = build_osculating_orbit(when_j2000, planet_orbit_end,   0.723328, 0.006786, 3.3945, 76.6285, 131.3854, (360*2)+5.386594); //1.60212 deg per d; 384 deg per 240d
      
      var result = current_osculating_orbit(osc_a, osc_b, when);
      return result;
    },
    add_physical: function(ephem){
      var i = degs(ephem.phase);
      ephem.mag = -4.40 + 5 * Math.log10(ephem.r * ephem.Δ) + 0.0009*i + 0.000239*i*i - 0.00000065*i*i*i;
      ephem.size = size_rads(8.34, ephem.Δ); 
    },
    ephem: function(when){
      var ephem = position_from_orbit(this.orbit(when), when);
      this.add_physical(ephem);
      return ephem;
    }
  };
  
  var mars = {
    name: 'Mars',
    symbol: '♂',
    orbit: function(when){
      var osc_a = build_osculating_orbit(when_j2000, planet_orbit_start, 1.523723, 0.093400, 1.8481, 49.5040, 336.1993, 71.7985);
      var osc_b = build_osculating_orbit(when_j2000, planet_orbit_end,   1.523603, 0.093505, 1.8481, 49.5009, 336.1891, 197.5675); //0.52402 deg per d; 126d in 240d 
      var result = current_osculating_orbit(osc_a, osc_b, when);
      return result;
    },
    add_physical: function(ephem){
      var i = degs(ephem.phase);
      ephem.mag = -1.52 + 5 * Math.log10(ephem.r * ephem.Δ) + 0.016*i;
      ephem.size = size_rads(4.68, ephem.Δ); 
    },
    ephem: function(when){
      var ephem = position_from_orbit(this.orbit(when), when);
      this.add_physical(ephem);
      return ephem;
    }
  };
  
  var jupiter = {
    name: 'Jupiter',
    symbol: '♃',
    orbit: function(when){
      var osc_a = build_osculating_orbit(when_j2000, planet_orbit_start, 5.202971, 0.048778, 1.3037, 100.5151, 14.1082, 257.3203);
      var osc_b = build_osculating_orbit(when_j2000, planet_orbit_end,   5.203385, 0.048733, 1.3037, 100.5162, 14.0236, 277.2581); //19.94d per 240d 
      var result = current_osculating_orbit(osc_a, osc_b, when);
      return result;
    },
    add_physical: function(ephem){
      var i = degs(ephem.phase);
      ephem.mag = -9.40 + 5 * Math.log10(ephem.r * ephem.Δ) + 0.005*i;
      ephem.size = size_rads(98.44, ephem.Δ); 
    },
    ephem: function(when){
      var ephem = position_from_orbit(this.orbit(when), when);
      this.add_physical(ephem);
      return ephem;
    }
  };
  
  var saturn = {
    name: 'Saturn',
    symbol: '♄',
    orbit: function(when){
      var osc_a = build_osculating_orbit(when_j2000, planet_orbit_start, 9.571397, 0.051705, 2.4863, 113.5945, 92.7441, 284.6606);
      var osc_b = build_osculating_orbit(when_j2000, planet_orbit_end,   9.572264, 0.051794, 2.4862, 113.5949, 92.0414, 292.7290); //8.0d per 240d 
      var result = current_osculating_orbit(osc_a, osc_b, when);
      return result;
    },
    add_physical: function(ephem){
      var i = degs(ephem.phase);
      var FUDGE_FACTOR_FOR_MISSING_RINGS =  -1.0; //gives the correct answer near opposition, 2019-07-09
      ephem.mag = -8.68 + 5 * Math.log10(ephem.r * ephem.Δ) + FUDGE_FACTOR_FOR_MISSING_RINGS;
      ephem.size = size_rads(82.73, ephem.Δ); 
    },
    ephem: function(when){
      var ephem = position_from_orbit(this.orbit(when), when);
      this.add_physical(ephem);
      return ephem;
    }
  };
  
  var uranus = {
    name: 'Uranus',
    symbol: '⛢',
    orbit: function(when){
      var osc_a = build_osculating_orbit(when_j2000, planet_orbit_start, 19.128240, 0.048830, 0.7708, 74.0658, 174.5193, 35.3703);
      var osc_b = build_osculating_orbit(when_j2000, planet_orbit_end,   19.147300, 0.047739, 0.7706, 74.0832, 174.1226, 38.1445);  
      var result = current_osculating_orbit(osc_a, osc_b, when);
      return result;
    },
    add_physical: function(ephem){
      var i = degs(ephem.phase);
      ephem.mag = -7.19 + 5 * Math.log10(ephem.r * ephem.Δ); 
      ephem.size = size_rads(35.02, ephem.Δ); 
    },
    ephem: function(when){
      var ephem = position_from_orbit(this.orbit(when), when);
      this.add_physical(ephem);
      return ephem;
    }
  };
  
  var neptune = {
    name: 'Neptune',
    symbol: '♆',
    orbit: function(when){
      var osc_a = build_osculating_orbit(when_j2000, planet_orbit_start, 30.090770, 0.007148, 1.7711, 131.7946, 30.7960, 346.5588);
      var osc_b = build_osculating_orbit(when_j2000, planet_orbit_end,   30.142030, 0.008522, 1.7704, 131.7801, 23.5110, 348.0054);  
      var result = current_osculating_orbit(osc_a, osc_b, when);
      return result;
    },
    add_physical: function(ephem){
      var i = degs(ephem.phase);
      ephem.mag = -6.87 + 5 * Math.log10(ephem.r * ephem.Δ); 
      ephem.size = size_rads(33.50, ephem.Δ); 
    },
    ephem: function(when){
      var ephem = position_from_orbit(this.orbit(when), when);
      this.add_physical(ephem);
      return ephem;
    }
  };
  
  var planets = {
    sun: sun,  moon: moon,  mercury: mercury,  venus: venus,  
    mars: mars,  jupiter: jupiter,  saturn: saturn,  uranus: uranus,  
    neptune: neptune
  };

  /* All minor planets share the same ephem function, which calcs the position from a fixed (osculating) orbit. */
  var base_minor_planet = {
    add_physical: function(ephem){
      // Meeus, page 217
      var φ1 = Math.exp(-3.33 * Math.pow(Math.tan(ephem.phase/2), 0.63));
      var φ2 = Math.exp(-1.87 * Math.pow(Math.tan(ephem.phase/2), 1.22));
      ephem.mag = this.H + 5 * Math.log10(ephem.r * ephem.Δ) - 2.5 * Math.log10((1-this.G)*φ1 + this.G*φ2); 
    },
    ephem: function(when){
      var ephem = position_from_orbit(this.orbit, when);
      this.add_physical(ephem);
      return ephem;
    }
  };

  /* Base data for a minor planet*/
  var build_minor_planet = function(name, id, H, G, orbit){
    var result = Object.create(base_minor_planet);
    result.name = name;
    result.alt_name = id + ' ' + name;
    result.orbit = orbit;
    result.H = H; // magnitude at 1 AU
    result.G = G; // 'slope parameter' for magnitude formula
    return result;
  };

  var ceres = build_minor_planet('Ceres', 1, 3.4, 0.12, {
      equinox: when_j2000, 
      epoch: when("UT 2019-04-27"),
      a: 2.769165148633284,
      e: 0.0760090265983052,
      i: rads(10.59406719506626),
      Ω: rads(80.30553090445737),
      ω: rads(73.59769469844186),
      M0: rads(77.37209751948711),
      n: rads(0.2138852265918273) 
    }
  );
  
  var vesta = build_minor_planet('Vesta', 4, 3.0, 0.32, {
      equinox: when_j2000, 
      epoch: when("UT 2019-04-27"),
      a: 2.361417893971132,
      e: 0.08872145978916499,
      i: rads(7.141771087953539),
      Ω: rads(103.8108039679376),
      ω: rads(150.7285410950405),
      M0: rads(95.86193772405923),
      n: rads(0.2716094018601813) 
    }
  );
  
  var eunomia = build_minor_planet('Eunomia', 15, 5.2, 0.23, {
      equinox: when_j2000, 
      epoch: when("UT 2020-05-31"),
      a: 2.643689741651421,
      e: 0.1861776979566804,
      i: rads(11.7536193269134),
      Ω: rads(292.9347962864794),
      ω: rads(98.59525178328343),
      M0: rads(15.02560926150444),
      n: rads(0.2292917057574047) 
    }
  );
  
  var euterpe = build_minor_planet('Euterpe', 27, 7.0, 0.15, {
      equinox: when_j2000, 
      epoch: when("UT 2020-05-31"),
      a: 2.346231773867898,
      e: 0.173130990523773,
      i: rads(1.583696508046424),
      Ω: rads(94.78724571687158),
      ω: rads(356.3779263487263),
      M0: rads(85.06235425066902),
      n: rads(0.274250675533197) 
    }
  );
  
  var melpomene = build_minor_planet('Melpomene', 18, 6.60, 0.25, {
      equinox: when_j2000, 
      epoch: when("UT 2020-05-31"),
      a: 2.295708701265865,
      e: 0.2175053249662481,
      i: rads(10.13159071229021),
      Ω: rads(150.3651871549247),
      ω: rads(228.0776817288102),
      M0: rads(20.46210040363494),
      n: rads(0.2833537100873361) 
    }
  );
  
  var astraea = build_minor_planet('Astraea', 5, 6.9, 0.15, {
      equinox: when_j2000, 
      epoch: when("UT 2019-05-31"),
      a: 2.5740372797851,
      e: 0.1909133446607621,
      i: rads(5.367427194855152),
      Ω: rads(141.5710246346721),
      ω: rads(358.6484193130638),
      M0: rads(17.84634390029505),
      n: rads(0.2386612028566109) 
    }
  );
  
  var minor_planets = {
    ceres: ceres, vesta: vesta, eunomia: eunomia, euterpe: euterpe, melpomene: melpomene, astraea: astraea 
  };
  
  /* 
   Find a thing in an object. Match simply on the name of a property of the given object. 
   If nothing found, inspect each property x of the object, and try to match raw_name to x.name and x.alt_name.
   Case-insensitive. 
  */  
  var find_thing = function(raw_name, object){
    var name = raw_name.toLowerCase();
    var result = object[name];
    if (!result){
      //since the base data structure is an object, and not an array, we need to scan the object's properties
      for (prop in object){
        if (object.hasOwnProperty(prop)){
          if (object[prop].name && object[prop].name.toLowerCase() === name){
            result = object[prop];
            break;
          }
          if (object[prop].alt_name && object[prop].alt_name.toLowerCase() === name){
            result = object[prop];
            break;
          }
        }
      }
    }
    return result;
  };
  
  /* All comets share the same ephem function, which calcs the position from a fixed (osculating) orbit. */
  //BAA 2015: picks out those with peak mag <= 12, and whose elong from the sun at peak > x
  //Having trouble with magnitudes; since they are never exact anyway, I'm going to go with 
  //an approx, hard-coded but fairly recent value from the BAA.
  var base_comet = {
    ephem: function(when){
      return position_from_orbit(this.orbit, when);
    }
  };
  /* 
   Base data for a comet.
   Periodic comets will be elliptical (orbit.a), while non-periodic comets will be quasi-parabolic (orbit.q). 
  */
  var build_comet = function(name, alt_name, mag, trend, when_vis, orbit){
    var result = Object.create(base_comet);
    result.name = name;
    result.alt_name = alt_name;
    result.orbit = orbit;
    result.mag = mag;
    result.trend = trend;
    result.when_vis = when_vis;
    return result;
  };
  /* 
   DON'T CHANGE THIS TO MATCH the Minor Planet Center. 
   For testing only. Meeus pg 217. Excellent agreement: 15" and 17" of arc. 
  */
  var enke_test = build_comet('enke', '2P/Encke', 14.0, 7.0, {
      equinox: when_j2000,
      epoch: when("TT 1990-10-28 13:05"),  
      a: 2.2091404,
      e: 0.8502196,
      i: rads(11.94524),
      Ω: rads(334.75006),
      ω: rads(186.23352),
      T: when("TT 1990-10-28 13:05") 
    }
  );
  
  var comets = {
    /*testing only enke_test: enke_test,*/
  };
  
  /* Match name to an object, compute its ephemeris, then apply the options. */
  var position = function(name_raw, when, options){
    var ephem = null;
    var name = name_raw.toLowerCase();
    var minor_planet, messier, caldwell, comet;
    if (planets[name]){
      ephem = planets[name].ephem(when); 
    }
    else if ((minor_planet = find_thing(name, minor_planets))){
      ephem = minor_planet.ephem(when);
    }
    else if (( messier = find_messier(name))){
      ephem = fixed_ephem(messier);
    }
    else if (( caldwell = find_caldwell(name))){
      ephem = fixed_ephem(caldwell);
    }
    else if ((comet = find_thing(name, comets))){
      ephem = comet.ephem(when); 
    }
    apply_options(ephem, when, options);
    return ephem;
  };

  var apply_options = function(ephem, when, options){
    if(options){
      apply_option_equinox(ephem, options);
      apply_option_where(ephem, when, options); 
      //rounding or sig figs should go here, if present
      apply_option_angular_units(ephem, options);
    }
  };
  
  /* Apply precession.*/  
  var apply_option_equinox = function(ephem, options){
    if (options.equinox){
      apply_precession(ephem, options.equinox, options.precession_angles /*may be absent*/);
    }
  };

  /* Add (a,A,h) to the ephemeris. This only makes sense when the ephemeris is with respect to the mean equinox of date. */
  var apply_option_where = function(ephem, when, options){
    if (options.where){
      convert_αδ_to_aA(ephem, options.where, when);
    }
  };
  
  /* Convert angles to the desired units.  */
  var apply_option_angular_units = function(ephem, options){
    if (options.units){
      if (options.units === 'degs'){
         convert_all_angles_to_degs(ephem);  
      }
      else if (options.units === 'degs_sexagesimal'){
         convert_all_angles_to_degs_sexagesimal(ephem);  
      }
    }
  };
  
  var position_from_orbit = function(orbit, when, options) {
    add_derived_orbital_items_to(orbit, when);
    var consts = intermediate_orbit_params(orbit, when);
    var ephem = find_position_in_orbit(orbit, when); // radius vector, true anomaly only
    ephem.when = when;
    find_heliocentric_xyz(ephem, consts, orbit);
    var sun = position_sun(when); //mean equinox of date
    apply_precession(sun, orbit.equinox);
    convert_xyz_to_XYZ(ephem, sun);
    convert_XYZ_to_αδ(ephem);
    ephem.equinox = orbit.equinox;
    convert_αδ_to_λβ(ephem, when);
    convert_ra_to_zodiac_sign(ephem);
    add_physical_ephem(ephem, sun);
    if (options) {
      apply_options(ephem, when, options);
    }
    return ephem; 
  };

  /* 
   The signed difference in ecliptical longitude (λ) with the sun, in the range -pi..+pi.
   Note that this is not the same as the object's elongation from the sun.
   The sign of the return value states if the object is east/west (+/-) of the Sun.   
   0+ to +pi means east of the Sun, -pi to 0- means west of the Sun. 
   Note the discontinuity between +pi and -pi.
   For objects that stray far from the ecliptic, this number makes less sense.
  */
  var delta_longitude_between = function(ephem, sun){
    var result = ephem.λ - sun.λ; //-2pi..+2pi
    if (result < 0){
      result = result + 2*Math.PI;
    }
    if (result > Math.PI){
      result = result - 2*Math.PI;
    }
    //result is now in -pi..+pi, with a discontinuity at +/-pi
    return result;
  };
  
  /* A sign is given to the elong from the Sun, according to the relative longitude λ of the object and the Sun. */
  var add_physical_ephem = function(ephem, sun){
    var sign = delta_longitude_between(ephem, sun) < 0 ? -1 : 1;
    ephem.elong = sign * elongation_between(ephem, sun); //+0..pi for east, -0..-pi for west 
    
    /* Meeus page 267    
    var numer = Math.pow(ephem.r,2) + Math.pow(ephem.Δ,2) - Math.pow(sun.Δ,2); 
    var denom = 2 * ephem.r * ephem.Δ;
    ephem.phase = Math.acos(numer/denom); // 0..pi
    */
    //Meeus page 316 - this is a better formula for phase because it works for all cases, including the moon, for 
    //which the moon-sun distance is not found.
    var numer = sun.Δ * Math.sin(ephem.elong);
    var denom = ephem.Δ - sun.Δ * Math.cos(ephem.elong);
    ephem.phase = Math.abs(Math.atan2(numer,denom)); // 0..pi
    
    ephem.illum = (1 + Math.cos(ephem.phase))/2; //0..1
    
    //Meeus 316, PA of the bright limb; intended mostly for the Moon, but works for planets too.
    numer = Math.cos(sun.α) * Math.sin(sun.α - ephem.α);
    denom = Math.sin(sun.δ) * Math.cos(ephem.δ) - Math.cos(sun.δ) * Math.sin(ephem.δ) * Math.cos(sun.α - ephem.α);
    ephem.χ = in2pi(Math.atan2(numer, denom));  // 0..2pi position angle of the midpoint of the bright limb; usually near pi/2 or 3pi/2
  }; 

  /* The start-end solar longitudes are approximate. 
    Data: 
     http://www.imo.net/files/data/vmdb/vmdbrad.txt (most data - where is fwhm from?)
     http://www.imo.net/members/imo_showers/working_shower_list  (rates) 
  */  
  var quadrantids = { 
    name: 'Quadrantids',
    symbol: 'QUA',
    λ: {
      start: rads(279.99), //jan 1
      peak: rads(283.16),
      end: rads(284.07) //jan 5
    },
    equinox: when_j2000,
    zhr: 120,  
    fwhm: 0.6, 
    radiant: {
      α:rads(230.1), 
      δ:rads(48.5),
      α_dot:rads(0.40), 
      δ_dot:rads(-0.20)
    },
    v: 42.7, 
    r: 2.1
  };
  var april_lyrids = { 
    name: 'April Lyrids',
    symbol: 'LYR',
    λ: {
      start: rads(29.38), //apr 19
      peak: rads(32.08),
      end: rads(34.26) //apr 24
    },
    equinox: when_j2000,
    zhr: 20,  
    fwhm: 1.3, 
    radiant: {
      α:rads(271.4), 
      δ:rads(33.6),
      α_dot:rads(1.1), 
      δ_dot:rads(0)
    },
    v: 47.6, 
    r: 2.9
  };
  var eta_aquarids = { 
    name: 'η Aquarids',
    symbol: 'ETA',
    λ: {
      start: rads(41.06), //may 1
      peak: rads(45.50),
      end: rads(47.85) //may 8
    },
    equinox: when_j2000,
    zhr: 60,  
    fwhm: 5, 
    radiant: {
      α:rads(338.0), 
      δ:rads(-1.0),
      α_dot:rads(0.9), 
      δ_dot:rads(0.4)
    },
    v: 66.0, 
    r: 2.70
  };
  var s_delta_aquarids = { 
    name: 'S δ Aquarids',
    symbol: 'SDA',
    λ: {
      start: rads(112.96), //jul 15
      peak: rads(125.0),
      end: rads(142.61) //aug 15
    },
    equinox: when_j2000,
    zhr: 20,  
    fwhm: 8, 
    radiant: {
      α:rads(339.0), 
      δ:rads(-16.0),
      α_dot:rads(0.7), 
      δ_dot:rads(0.18)
    },
    v: 41.4,  
    r: 3.2 
  };
  var perseids = { 
    name: 'Perseids',
    symbol: 'PER',
    λ: {
      start: rads(122.50), //jul 25
      peak: rads(140.0),
      end: rads(145.50) //aug 18
    },
    equinox: when_j2000,
    zhr: 90,  
    fwhm: 2.0, 
    radiant: {
      α:rads(46.2), 
      δ:rads(57.4),
      α_dot:rads(1.35), 
      δ_dot:rads(0.12)
    },
    v: 60.0, 
    r: 2.6
  };
  var orionids = { 
    name: 'Orionids',
    symbol: 'ORI',
    λ: {
      start: rads(203.07), //oct 16
      peak: rads(208.40),
      end: rads(214.01) //oct 27
    },
    equinox: when_j2000,
    zhr: 20,  
    fwhm: 2.0, 
    radiant: {
      α:rads(94.5), 
      δ:rads(15.8),
      α_dot:rads(0.7), 
      δ_dot:rads(0.1)
    },
    v: 67.0, 
    r: 2.9
  };
  var leonids = { 
    name: 'Leonids',
    symbol: 'LEO',
    λ: {
      start: rads(233.06), //nov 15
      peak: rads(235.16),
      end: rads(237.09) //nov 19
    },
    equinox: when_j2000,
    zhr: 20,  
    fwhm: 1.0, 
    radiant: {
      α:rads(152.3), 
      δ:rads(22.2),
      α_dot:rads(0.7), 
      δ_dot:rads(-0.42)
    },
    v: 71.1, 
    r: 2.5
  };
  var geminids = { 
    name: 'Geminids',
    symbol: 'GEM',
    λ: {
      start: rads(255.32), //dec 7
      peak: rads(262.0),
      end: rads(264.47) //dec 16
    },
    equinox: when_j2000,
    zhr: 120,  
    fwhm: 1.0, 
    radiant: {
      α:rads(112.3), 
      δ:rads(32.5),
      α_dot:rads(1.02), 
      δ_dot:rads(-0.07)
    },
    v: 35.0, 
    r: 2.6
  };
  var puppids_velids = { 
    name: 'Puppids-Velids',
    symbol: 'PUP',
    λ: {
      start: rads(248.7), //dec 1
      peak: rads(255.0),
      end: rads(262.9) //dec 15
    },
    equinox: when_j2000,
    zhr: 10,  
    fwhm: 1.0, 
    radiant: {
      α:rads(123.0), 
      δ:rads(-45.0),
      α_dot:rads(0.50), 
      δ_dot:rads(0.00)
    },
    v: 40.0, 
    r: 2.9
  };
  
  //ursids: weak
  //taurids: i'm confused: n, s, and not
  
  var meteor_showers = [quadrantids, april_lyrids, eta_aquarids, s_delta_aquarids, perseids, orionids, leonids, geminids, puppids_velids];
  
  /*
   Return a 'when' object for the peak date-time of the shower. 
   The peak time is when the sun's longitude reaches a certain value, characteristic of each shower. 
   'active_showers' should be called first before this function, to make sure the peak time is nearby, 
   and this function converges responsibly to the initial guess. 
  */
  var when_shower_peaks = function(meteor_shower, when_initial_guess){
    var when_guess = when_initial_guess;
    var small_Δλ = rads(0.01);
    var sun = planets.sun.position(when_guess);
    //it would be nice if the sun's position was j2000 to begin with
    apply_precession(sun, meteor_shower.equinox);
    var Δλ = (meteor_shower.λ.peak - sun.λ);
    var Δt_guess = 0; //estimated diff between when_guess and the shower's peak time
    while (Math.abs(Δλ) > small_Δλ){
      //improve the guess
      Δt_guess = Δλ * 365.256363/(2*Math.PI); //guess of fractional days diff; use mean motion of the sun per day  
      Δt_guess = Δt_guess * MSEC_PER_DAY; //msecs
      when_guess = when_from_utc(new Date(when_guess.date.getTime() + Δt_guess));
      sun = planets.sun.position(when_guess);
      apply_precession(sun, meteor_shower.equinox);
      Δλ = (meteor_shower.λ.peak - sun.λ);
    }
    return when_guess;
  };
  
  /* 
   Return the 'current' position of the radiant, as a function of 'when'.
   IMPORTANT: if you next calc (a,A), you must first apply precession, and make sure you are using the 
   current mean equinox of date.
  */
  var current_radiant_position = function(when, meteor_shower, when_peak){
    var Δd = when.jd - when_peak.jd; //fractional days
    var ephem = {
      α: meteor_shower.radiant.α + meteor_shower.radiant.α_dot * Δd,  
      δ: meteor_shower.radiant.δ + meteor_shower.radiant.δ_dot * Δd,
      equinox: meteor_shower.equinox
    };
    return ephem;
  };
  
  //when is it peaking?
  //how many can i expect to see, at what time?
  //   it varies as zhr*sin(a), where a is the radiant's altitude
  //   it also varies with the limiting magnitude
  //   https://en.wikipedia.org/wiki/Zenithal_hourly_rate
  //   in this case, 'where' might have secondary information about sky brightness (limiting mag), and % clear (obstructions) 
  //when is it worth viewing?
  //where is the radiant in my sky? above horizon or not? to calc a, A, you need to use the current equinox
  
  /* Returns an array of meteor showers that are active at the given time. */
  var active_showers = function(when){
    var result = [];
    var sun = planets.sun.position(when);
    for(var i = 0; i < meteor_showers.length; ++i){
      if (meteor_showers[i].λ.start < sun.λ && sun.λ < meteor_showers[i].λ.end){
        result.push(meteor_showers[i]);
      }
    }
    return result;
  };
  
  /* 
   Radiant position in right ascencion and declination.
   The options must always should always include the current equinox, in order to get the proper (a,A,h).
   If the radiant is above the horizon, then sin of the radiant's altitude is added as a property.
  */
  var current_meteor_showers = function(when, options){
    var result = [];
    var current_showers = active_showers(when);
    for (var idx=0; idx < current_showers.length; ++idx){
      //simply copy some base data over into a new object
      //avoid 'v,r' as property names, since they conflict with other ephem items (radius vector and anomaly)
      var shower = {
        name: current_showers[idx].name,
        symbol: current_showers[idx].symbol,
        v: current_showers[idx].v,
        r: current_showers[idx].r,
        zhr: current_showers[idx].zhr
      };
      var when_peak = when_shower_peaks(current_showers[idx], when);
      var ephem = current_radiant_position(when, current_showers[idx], when_peak);
      //we need to apply options piecemeal here, since some intermediate logic needs rads to still be in place
      apply_option_equinox(ephem, options);
      apply_option_where(ephem, when, options); 
      if (ephem.a && ephem.a > 0){
        shower.zhr_factor = Math.sin(ephem.a);
        if (options.where.limiting_mag){
          shower.zhr_factor = shower.zhr_factor / Math.pow(shower.r_idx, 6.5-options.where.limiting_mag);
        }        
      }
      else {
          shower.zhr_factor = 0;
      }
      apply_option_angular_units(ephem, options);
      shower.peak_time = when_peak;
      shower.radiant = ephem;
      result.push(shower);
    }
    return result;
  };
  
  /* Physical ephemeris for Jupiter and its 4 Galilean moons. */
  var physical_jupiter = function(when){
    var d = when.jd_tt - JD_J2000;
    var V = rads(172.74 + 0.00111588*d);
    var M = rads(357.529 + 0.9856003*d);
    var N = rads(20.020 + 0.0830853*d + 0.329*Math.sin(V));
    var J = rads(66.115 + 0.9025179*d - 0.329*Math.sin(V))
    var A = rads(1.915*Math.sin(M) + 0.020*Math.sin(2*M));
    var B = rads(5.555*Math.sin(N) + 0.168*Math.sin(2*N));
    var K = J + A - B;
    var R = 1.00014 - 0.01671*Math.cos(M) - 0.00014*Math.cos(2*M);
    var r = 5.20872 - 0.25208*Math.cos(N) - 0.00611*Math.cos(2*N);
    var Δ = Math.sqrt(r*r + R*R - 2*r*R*Math.cos(K));
    var ψ = Math.asin(Math.sin(K)*R/Δ); // -pi/2...+pi/2
    var d_corr = d - Δ/173;
    var ω1 = in2pi(rads(210.98) + rads(877.8169088)*d_corr + ψ - B);
    var ω2 = in2pi(rads(187.23) + rads(870.1869088)*d_corr + ψ - B);
    var λ = rads(34.35) + rads(0.083091)*d + rads(0.329)*Math.sin(V) + B;
    var D_S = rads(3.12)*Math.sin(λ + rads(42.8));
    var D_E = D_S - rads(2.22)*Math.sin(ψ)*Math.cos(λ + rads(22)) - rads(1.30)*((r - Δ)/Δ)*Math.sin(λ - rads(100.5));
    //satellites
    var u = new Array(5);
    u[1] = rads(163.8067) + rads(203.4058643)*(d_corr) +  ψ - B;
    u[2] = rads(358.4108) + rads(101.2916334)*(d_corr) +  ψ - B;
    u[3] = rads(5.7129) + rads(50.2345179)*(d_corr) +  ψ - B;
    u[4] = rads(224.8151) + rads(21.4879801)*(d_corr) +  ψ - B;
    var G = rads(331.18 + 50.310482*d_corr);
    var H = rads(87.40 + 21.569231*d_corr);
    var dist = new Array(5);
    //use the 'uncorrected' u's
    dist[1] = 5.9073 - 0.0244*Math.cos(2*(u[1] - u[2]));
    dist[2] = 9.3991 - 0.0882*Math.cos(2*(u[2] - u[3]));
    dist[3] = 14.9924 - 0.0216*Math.cos(G);
    dist[4] = 26.3699 - 0.1935*Math.cos(H);
    //the 'corrected' u's
    u[1] = u[1] + rads(0.473)*Math.sin(2*(u[1] - u[2]));
    u[2] = u[2] + rads(1.065)*Math.sin(2*(u[2] - u[3]));
    u[3] = u[3] + rads(0.165)*Math.sin(G);
    u[4] = u[4] + rads(0.841)*Math.sin(H);
    var satellite_names = ['blank', 'Io', 'Europa', 'Ganymede', 'Callisto'];
    var satellite_symbols = ['blank', 'I', 'II', 'III', 'IV'];
    var satellites = [];
    for (var idx = 1; idx < u.length; ++idx){
      satellites.push(
        { 
          name: satellite_names[idx],
          symbol: satellite_symbols[idx],
          X: dist[idx]*Math.sin(u[idx]),
          Y: -1*dist[idx]*Math.cos(u[idx])*Math.sin(D_E)
        }
      );
    }
    return {
      ω1: ω1, 
      ω2: ω2,
      D_S: D_S,
      D_E: D_E,
      satellites: satellites
    };
  };

  /* Optical librations only (~8 deg); no physical (~0.04 deg) or topocentric librations (~1 deg) are included. */  
  var lunar_libration = function(when, moon){
    var I = rads(1.54242);
    if (moon === undefined) {
      moon = position_moon(when); //mean equinox; geometric, not apparent; neglect nutation
    }
    var T = when.T;
    var F = rads(93.2720993 + T * (483202.0175273 + T * (-0.0034029 - T * ((1 / 3526000)))));
    var Ω = rads(125.044555 - 1934.1361849*T + 0.0020762*T*T + T*T*T/467410);
    var W = moon.λ - Ω; //neglect nutation
    var numer = Math.sin(W)*Math.cos(moon.β)*Math.cos(I) - Math.sin(moon.β)*Math.sin(I);
    var denom = Math.cos(W)*Math.cos(moon.β);
    var A = Math.atan2(numer, denom); 
    var l_prime = in2pi(A - F);
    //keep the libration near 0 degrees, + or -
    if (l_prime > rads(90)){
      l_prime = l_prime - 2*Math.PI;
    }
    var b_prime = Math.asin(-1*Math.sin(W)*Math.cos(moon.β)*Math.sin(I) - Math.sin(moon.β)*Math.cos(I)); //-pi/2..+pi/2
    return {
      longitude: l_prime,
      latitude: b_prime,
      distance: moon.Δ // AU
    };
  };

  /*  
    Important: these are ordered by increasing 'when'. All times are in UT.
    The 'when' property here is simple text, not a 'when' object.
    Important: simple text comparison is used to bracket the times.
    Missing from raw output of the MICA tool: eclipses, X at aphelion, X stationary, ad hoc text
  */  
  var all_events = [
     {when:'UT 2021-01-02 14',text:'Earth at perihelion 0.983257060 AU'},
     {when:'UT 2021-01-02 22',text:'Regulus 4.72°S of Moon'},
     {when:'UT 2021-01-06 10',text:'Last Quarter 371387.602 km'},
     {when:'UT 2021-01-06 19',text:'Spica 7.01°S of Moon'},
     {when:'UT 2021-01-09 16',text:'Moon at perigee 367387.455 km'},
     {when:'UT 2021-01-09 21',text:'Saturn 1.67°N of Mercury'},
     {when:'UT 2021-01-10 03',text:'Antares 5.59°S of Moon'},
     {when:'UT 2021-01-11 11',text:'Jupiter 1.48°N of Mercury'},
     {when:'UT 2021-01-11 20',text:'Venus 1.48°N of Moon'},
     {when:'UT 2021-01-13 05',text:'New Moon 374128.716 km'},
     {when:'UT 2021-01-13 21',text:'Saturn 3.22°N of Moon'},
     {when:'UT 2021-01-14 01',text:'Jupiter 3.30°N of Moon'},
     {when:'UT 2021-01-14 08',text:'Mercury 2.32°N of Moon'},
     {when:'UT 2021-01-17 06',text:'Neptune 4.45°N of Moon'},
     {when:'UT 2021-01-20 21',text:'First Quarter 404060.540 km'},
     {when:'UT 2021-01-21 06',text:'Mars 5.05°N of Moon'},
     {when:'UT 2021-01-21 13',text:'Moon at apogee 404359.838 km'},
     {when:'UT 2021-01-22 00',text:'Uranus 1.72°S of Mars'},
     {when:'UT 2021-01-24 02',text:'Mercury at greatest elongation 18.6° East'},
     {when:'UT 2021-01-24 03',text:'Saturn 0.41°S of Sun'},
     {when:'UT 2021-01-24 05',text:'Aldebaran 4.76°S of Moon'},
     {when:'UT 2021-01-27 16',text:'Pollux 3.79°N of Moon'},
     {when:'UT 2021-01-28 19',text:'Full Moon 381520.878 km'},
     {when:'UT 2021-01-29 02',text:'Jupiter 0.52°S of Sun'},
     {when:'UT 2021-01-30 05',text:'Regulus 4.62°S of Moon'},
     {when:'UT 2021-02-03 00',text:'Spica 6.80°S of Moon'},
     {when:'UT 2021-02-03 19',text:'Moon at perigee 370116.150 km'},
     {when:'UT 2021-02-04 18',text:'Last Quarter 370327.722 km'},
     {when:'UT 2021-02-06 05',text:'Saturn 0.39°N of Venus'},
     {when:'UT 2021-02-06 09',text:'Antares 5.43°S of Moon'},
     {when:'UT 2021-02-08 14',text:'Mercury in inferior conjunction 3.63° North'},
     {when:'UT 2021-02-10 11',text:'Saturn 3.41°N of Moon'},
     {when:'UT 2021-02-10 20',text:'Venus 3.20°N of Moon'},
     {when:'UT 2021-02-10 22',text:'Jupiter 3.67°N of Moon'},
     {when:'UT 2021-02-11 03',text:'Mercury 8.28°N of Moon'},
     {when:'UT 2021-02-11 12',text:'Jupiter 0.44°N of Venus'},
     {when:'UT 2021-02-11 19',text:'New Moon 385520.198 km'},
     {when:'UT 2021-02-12 17',text:'Venus 4.81°S of Mercury'},
     {when:'UT 2021-02-13 17',text:'Neptune 4.32°N of Moon'},
     {when:'UT 2021-02-13 19',text:'Jupiter 4.21°S of Mercury'},
     {when:'UT 2021-02-17 16',text:'Uranus 3.03°N of Moon'},
     {when:'UT 2021-02-18 10',text:'Moon at apogee 404467.309 km'},
     {when:'UT 2021-02-18 23',text:'Mars 3.69°N of Moon'},
     {when:'UT 2021-02-19 19',text:'First Quarter 403282.882 km'},
     {when:'UT 2021-02-20 14',text:'Aldebaran 4.99°S of Moon'},
     {when:'UT 2021-02-24 02',text:'Pollux 3.67°N of Moon'},
     {when:'UT 2021-02-26 15',text:'Regulus 4.60°S of Moon'},
     {when:'UT 2021-02-27 08',text:'Full Moon 370594.927 km'},
     {when:'UT 2021-03-02 05',text:'Moon at perigee 365423.177 km'},
     {when:'UT 2021-03-02 07',text:'Spica 6.58°S of Moon'},
     {when:'UT 2021-03-05 07',text:'Jupiter 0.33°S of Mercury'},
     {when:'UT 2021-03-05 14',text:'Antares 5.17°S of Moon'},
     {when:'UT 2021-03-06 02',text:'Last Quarter 372121.853 km'},
     {when:'UT 2021-03-06 11',text:'Mercury at greatest elongation 27.3° West'},
     {when:'UT 2021-03-09 23',text:'Saturn 3.68°N of Moon'},
     {when:'UT 2021-03-10 16',text:'Jupiter 4.05°N of Moon'},
     {when:'UT 2021-03-11 00',text:'Neptune 1.07°S of Sun'},
     {when:'UT 2021-03-11 01',text:'Mercury 3.70°N of Moon'},
     {when:'UT 2021-03-13 00',text:'Venus 3.87°N of Moon'},
     {when:'UT 2021-03-13 03',text:'Neptune 4.26°N of Moon'},
     {when:'UT 2021-03-13 10',text:'New Moon 396124.281 km'},
     {when:'UT 2021-03-14 01',text:'Neptune 0.40°N of Venus'},
     {when:'UT 2021-03-17 02',text:'Uranus 2.72°N of Moon'},
     {when:'UT 2021-03-18 05',text:'Moon at apogee 405252.546 km'},
     {when:'UT 2021-03-19 18',text:'Mars 1.93°N of Moon'},
     {when:'UT 2021-03-19 22',text:'Aldebaran 5.27°S of Moon'},
     {when:'UT 2021-03-20 10',text:'Equinox'},
     {when:'UT 2021-03-21 15',text:'First Quarter 398487.986 km'},
     {when:'UT 2021-03-23 00',text:'Aldebaran 7.02°S of Mars'},
     {when:'UT 2021-03-23 11',text:'Pollux 3.45°N of Moon'},
     {when:'UT 2021-03-26 01',text:'Regulus 4.72°S of Moon'},
     {when:'UT 2021-03-26 07',text:'Venus in superior conjunction 1.35° South'},
     {when:'UT 2021-03-28 19',text:'Full Moon 362173.753 km'},
     {when:'UT 2021-03-29 16',text:'Spica 6.45°S of Moon'},
     {when:'UT 2021-03-29 19',text:'Neptune 1.40°N of Mercury'},
     {when:'UT 2021-03-30 06',text:'Moon at perigee 360309.171 km'},
     {when:'UT 2021-04-01 21',text:'Antares 4.92°S of Moon'},
     {when:'UT 2021-04-04 10',text:'Last Quarter 376583.184 km'},
     {when:'UT 2021-04-06 08',text:'Saturn 3.97°N of Moon'},
     {when:'UT 2021-04-07 07',text:'Jupiter 4.39°N of Moon'},
     {when:'UT 2021-04-09 11',text:'Neptune 4.33°N of Moon'},
     {when:'UT 2021-04-11 06',text:'Mercury 2.99°N of Moon'},
     {when:'UT 2021-04-12 03',text:'New Moon 403642.104 km'},
     {when:'UT 2021-04-12 10',text:'Venus 2.87°N of Moon'},
     {when:'UT 2021-04-13 12',text:'Uranus 2.49°N of Moon'},
     {when:'UT 2021-04-14 18',text:'Moon at apogee 406118.765 km'},
     {when:'UT 2021-04-16 05',text:'Aldebaran 5.47°S of Moon'},
     {when:'UT 2021-04-17 12',text:'Mars 0.13°N of Moon'},
     {when:'UT 2021-04-19 02',text:'Mercury in superior conjunction 0.57° South'},
     {when:'UT 2021-04-19 19',text:'Pollux 3.22°N of Moon'},
     {when:'UT 2021-04-20 07',text:'First Quarter 391306.287 km'},
     {when:'UT 2021-04-22 10',text:'Regulus 4.90°S of Moon'},
     {when:'UT 2021-04-22 23',text:'Uranus 0.25°N of Venus'},
     {when:'UT 2021-04-24 10',text:'Uranus 0.81°S of Mercury'},
     {when:'UT 2021-04-26 03',text:'Spica 6.46°S of Moon'},
     {when:'UT 2021-04-26 09',text:'Venus 1.29°S of Mercury'},
     {when:'UT 2021-04-27 04',text:'Full Moon 357616.113 km'},
     {when:'UT 2021-04-27 15',text:'Moon at perigee 357377.974 km'},
     {when:'UT 2021-04-29 07',text:'Antares 4.78°S of Moon'},
     {when:'UT 2021-04-30 20',text:'Uranus 0.41°S of Sun'},
     {when:'UT 2021-05-03 17',text:'Saturn 4.17°N of Moon'},
     {when:'UT 2021-05-03 20',text:'Last Quarter 383059.100 km'},
     {when:'UT 2021-05-04 21',text:'Jupiter 4.61°N of Moon'},
     {when:'UT 2021-05-06 18',text:'Neptune 4.43°N of Moon'},
     {when:'UT 2021-05-10 21',text:'Uranus 2.36°N of Moon'},
     {when:'UT 2021-05-11 03',text:'Aldebaran 7.99°S of Mercury'},
     {when:'UT 2021-05-11 19',text:'New Moon 406507.120 km'},
     {when:'UT 2021-05-11 22',text:'Moon at apogee 406511.956 km'},
     {when:'UT 2021-05-12 22',text:'Venus 0.71°N of Moon'},
     {when:'UT 2021-05-13 11',text:'Aldebaran 5.54°S of Moon'},
     {when:'UT 2021-05-13 18',text:'Mercury 2.14°N of Moon'},
     {when:'UT 2021-05-16 05',text:'Mars 1.47°S of Moon'},
     {when:'UT 2021-05-17 01',text:'Pollux 3.11°N of Moon'},
     {when:'UT 2021-05-17 06',text:'Mercury at greatest elongation 22.0° East'},
     {when:'UT 2021-05-17 23',text:'Aldebaran 5.87°S of Venus'},
     {when:'UT 2021-05-19 18',text:'Regulus 5.02°S of Moon'},
     {when:'UT 2021-05-19 19',text:'First Quarter 383678.778 km'},
     {when:'UT 2021-05-23 14',text:'Spica 6.52°S of Moon'},
     {when:'UT 2021-05-26 02',text:'Moon at perigee 357310.962 km'},
     {when:'UT 2021-05-26 11',text:'LUNAR ECLIPSE Full Moon 357460.951 km'},
     {when:'UT 2021-05-26 17',text:'Antares 4.76°S of Moon'},
     {when:'UT 2021-05-29 06',text:'Venus 0.42°N of Mercury'},
     {when:'UT 2021-05-31 01',text:'Saturn 4.18°N of Moon'},
     {when:'UT 2021-06-01 09',text:'Jupiter 4.63°N of Moon'},
     {when:'UT 2021-06-02 07',text:'Last Quarter 390480.536 km'},
     {when:'UT 2021-06-02 14',text:'Pollux 5.42°N of Mars'},
     {when:'UT 2021-06-03 01',text:'Neptune 4.46°N of Moon'},
     {when:'UT 2021-06-07 06',text:'Uranus 2.25°N of Moon'},
     {when:'UT 2021-06-08 02',text:'Moon at apogee 406227.881 km'},
     {when:'UT 2021-06-09 17',text:'Aldebaran 5.53°S of Moon'},
     {when:'UT 2021-06-10 11',text:'ANNULAR ECLIPSE New Moon 404245.794 km'},
     {when:'UT 2021-06-10 13',text:'Mercury 3.96°S of Moon'},
     {when:'UT 2021-06-11 01',text:'Mercury in inferior conjunction 3.14° South'},
     {when:'UT 2021-06-12 07',text:'Venus 1.47°S of Moon'},
     {when:'UT 2021-06-13 07',text:'Pollux 3.11°N of Moon'},
     {when:'UT 2021-06-13 20',text:'Mars 2.80°S of Moon'},
     {when:'UT 2021-06-16 00',text:'Regulus 5.01°S of Moon'},
     {when:'UT 2021-06-18 04',text:'First Quarter 377060.898 km'},
     {when:'UT 2021-06-19 22',text:'Spica 6.50°S of Moon'},
     {when:'UT 2021-06-21 04',text:'Solstice'},
     {when:'UT 2021-06-22 15',text:'Pollux 5.26°N of Venus'},
     {when:'UT 2021-06-23 04',text:'Antares 4.75°S of Moon'},
     {when:'UT 2021-06-23 10',text:'Moon at perigee 359956.064 km'},
     {when:'UT 2021-06-24 19',text:'Full Moon 361560.842 km'},
     {when:'UT 2021-06-27 09',text:'Saturn 4.03°N of Moon'},
     {when:'UT 2021-06-28 19',text:'Jupiter 4.45°N of Moon'},
     {when:'UT 2021-06-30 09',text:'Neptune 4.36°N of Moon'},
     {when:'UT 2021-07-01 21',text:'Last Quarter 397469.380 km'},
     {when:'UT 2021-07-04 15',text:'Uranus 2.08°N of Moon'},
     {when:'UT 2021-07-04 20',text:'Mercury at greatest elongation 21.6° West'},
     {when:'UT 2021-07-05 15',text:'Moon at apogee 405341.161 km'},
     {when:'UT 2021-07-05 22',text:'Earth at aphelion 1.016729224 AU'},
     {when:'UT 2021-07-06 23',text:'Aldebaran 5.56°S of Moon'},
     {when:'UT 2021-07-08 05',text:'Mercury 3.75°S of Moon'},
     {when:'UT 2021-07-10 01',text:'New Moon 397520.337 km'},
     {when:'UT 2021-07-10 13',text:'Pollux 3.15°N of Moon'},
     {when:'UT 2021-07-12 09',text:'Venus 3.26°S of Moon'},
     {when:'UT 2021-07-12 10',text:'Mars 3.78°S of Moon'},
     {when:'UT 2021-07-13 06',text:'Regulus 4.90°S of Moon'},
     {when:'UT 2021-07-13 07',text:'Mars 0.49°S of Venus'},
     {when:'UT 2021-07-17 05',text:'Spica 6.34°S of Moon'},
     {when:'UT 2021-07-17 10',text:'First Quarter 372334.538 km'},
     {when:'UT 2021-07-20 13',text:'Antares 4.67°S of Moon'},
     {when:'UT 2021-07-21 10',text:'Moon at perigee 364520.495 km'},
     {when:'UT 2021-07-21 19',text:'Regulus 1.17°S of Venus'},
     {when:'UT 2021-07-24 03',text:'Full Moon 369211.389 km'},
     {when:'UT 2021-07-24 17',text:'Saturn 3.82°N of Moon'},
     {when:'UT 2021-07-25 11',text:'Pollux 5.74°N of Mercury'},
     {when:'UT 2021-07-26 01',text:'Jupiter 4.17°N of Moon'},
     {when:'UT 2021-07-27 18',text:'Neptune 4.18°N of Moon'},
     {when:'UT 2021-07-29 16',text:'Regulus 0.68°S of Mars'},
     {when:'UT 2021-07-31 13',text:'Last Quarter 402464.703 km'},
     {when:'UT 2021-08-01 00',text:'Uranus 1.82°N of Moon'},
     {when:'UT 2021-08-01 14',text:'Mercury in superior conjunction 1.69° North'},
     {when:'UT 2021-08-02 06',text:'Saturn at opposition'},
     {when:'UT 2021-08-02 08',text:'Moon at apogee 404409.687 km'},
     {when:'UT 2021-08-03 07',text:'Aldebaran 5.70°S of Moon'},
     {when:'UT 2021-08-06 20',text:'Pollux 3.13°N of Moon'},
     {when:'UT 2021-08-08 14',text:'New Moon 387818.768 km'},
     {when:'UT 2021-08-09 03',text:'Mercury 3.38°S of Moon'},
     {when:'UT 2021-08-09 12',text:'Regulus 4.81°S of Moon'},
     {when:'UT 2021-08-10 01',text:'Mars 4.30°S of Moon'},
     {when:'UT 2021-08-11 07',text:'Venus 4.29°S of Moon'},
     {when:'UT 2021-08-11 18',text:'Regulus 1.17°S of Mercury'},
     {when:'UT 2021-08-13 10',text:'Spica 6.11°S of Moon'},
     {when:'UT 2021-08-15 15',text:'First Quarter 370026.344 km'},
     {when:'UT 2021-08-16 19',text:'Antares 4.47°S of Moon'},
     {when:'UT 2021-08-17 09',text:'Moon at perigee 369124.445 km'},
     {when:'UT 2021-08-19 04',text:'Mars 0.08°N of Mercury'},
     {when:'UT 2021-08-20 00',text:'Jupiter at opposition'},
     {when:'UT 2021-08-20 22',text:'Saturn 3.70°N of Moon'},
     {when:'UT 2021-08-22 05',text:'Jupiter 3.96°N of Moon'},
     {when:'UT 2021-08-22 12',text:'Full Moon 379232.140 km'},
     {when:'UT 2021-08-24 02',text:'Neptune 4.03°N of Moon'},
     {when:'UT 2021-08-28 09',text:'Uranus 1.53°N of Moon'},
     {when:'UT 2021-08-30 02',text:'Moon at apogee 404099.867 km'},
     {when:'UT 2021-08-30 07',text:'Last Quarter 404073.486 km'},
     {when:'UT 2021-08-30 15',text:'Aldebaran 5.95°S of Moon'},
     {when:'UT 2021-09-03 05',text:'Pollux 2.99°N of Moon'},
     {when:'UT 2021-09-05 06',text:'Spica 1.74°S of Venus'},
     {when:'UT 2021-09-05 20',text:'Regulus 4.82°S of Moon'},
     {when:'UT 2021-09-07 01',text:'New Moon 377018.912 km'},
     {when:'UT 2021-09-07 16',text:'Mars 4.24°S of Moon'},
     {when:'UT 2021-09-08 20',text:'Mercury 6.52°S of Moon'},
     {when:'UT 2021-09-09 17',text:'Spica 5.91°S of Moon'},
     {when:'UT 2021-09-10 02',text:'Venus 4.08°S of Moon'},
     {when:'UT 2021-09-11 10',text:'Moon at perigee 368461.354 km'},
     {when:'UT 2021-09-13 01',text:'Antares 4.21°S of Moon'},
     {when:'UT 2021-09-13 21',text:'First Quarter 370428.263 km'},
     {when:'UT 2021-09-14 04',text:'Mercury at greatest elongation 26.8° East'},
     {when:'UT 2021-09-14 09',text:'Neptune at opposition'},
     {when:'UT 2021-09-17 03',text:'Saturn 3.76°N of Moon'},
     {when:'UT 2021-09-18 07',text:'Jupiter 3.96°N of Moon'},
     {when:'UT 2021-09-20 09',text:'Neptune 4.00°N of Moon'},
     {when:'UT 2021-09-21 00',text:'Full Moon 389988.036 km'},
     {when:'UT 2021-09-22 19',text:'Equinox'},
     {when:'UT 2021-09-23 12',text:'Spica 1.67°N of Mercury'},
     {when:'UT 2021-09-24 16',text:'Uranus 1.35°N of Moon'},
     {when:'UT 2021-09-26 22',text:'Moon at apogee 404640.374 km'},
     {when:'UT 2021-09-26 23',text:'Aldebaran 6.20°S of Moon'},
     {when:'UT 2021-09-29 02',text:'Last Quarter 401694.794 km'},
     {when:'UT 2021-09-30 13',text:'Pollux 2.77°N of Moon'},
     {when:'UT 2021-09-30 15',text:'Spica 1.74°N of Mercury'},
     {when:'UT 2021-10-03 06',text:'Regulus 4.95°S of Moon'},
     {when:'UT 2021-10-06 10',text:'Mars 3.55°S of Moon'},
     {when:'UT 2021-10-06 11',text:'New Moon 367080.665 km'},
     {when:'UT 2021-10-06 18',text:'Mercury 6.91°S of Moon'},
     {when:'UT 2021-10-07 01',text:'Spica 5.82°S of Moon'},
     {when:'UT 2021-10-08 04',text:'Mars 0.65°N of Sun'},
     {when:'UT 2021-10-08 17',text:'Moon at perigee 363385.698 km'},
     {when:'UT 2021-10-09 08',text:'Mars 2.86°N of Mercury'},
     {when:'UT 2021-10-09 16',text:'Mercury in inferior conjunction 1.90° South'},
     {when:'UT 2021-10-09 19',text:'Venus 2.86°S of Moon'},
     {when:'UT 2021-10-10 07',text:'Antares 4.00°S of Moon'},
     {when:'UT 2021-10-13 03',text:'First Quarter 373570.663 km'},
     {when:'UT 2021-10-14 07',text:'Saturn 3.93°N of Moon'},
     {when:'UT 2021-10-15 10',text:'Jupiter 4.14°N of Moon'},
     {when:'UT 2021-10-16 14',text:'Antares 1.47°S of Venus'},
     {when:'UT 2021-10-17 14',text:'Neptune 4.11°N of Moon'},
     {when:'UT 2021-10-20 06',text:'Spica 2.81°S of Mars'},
     {when:'UT 2021-10-20 15',text:'Full Moon 399419.511 km'},
     {when:'UT 2021-10-21 22',text:'Uranus 1.33°N of Moon'},
     {when:'UT 2021-10-24 06',text:'Aldebaran 6.36°S of Moon'},
     {when:'UT 2021-10-24 15',text:'Moon at apogee 405615.069 km'},
     {when:'UT 2021-10-25 05',text:'Mercury at greatest elongation 18.4° West'},
     {when:'UT 2021-10-27 21',text:'Pollux 2.59°N of Moon'},
     {when:'UT 2021-10-28 20',text:'Last Quarter 395946.842 km'},
     {when:'UT 2021-10-29 21',text:'Venus at greatest elongation 47.0° East'},
     {when:'UT 2021-10-30 15',text:'Regulus 5.10°S of Moon'},
     {when:'UT 2021-11-01 02',text:'Spica 4.44°S of Mercury'},
     {when:'UT 2021-11-03 12',text:'Spica 5.84°S of Moon'},
     {when:'UT 2021-11-03 19',text:'Mercury 1.22°S of Moon'},
     {when:'UT 2021-11-04 05',text:'Mars 2.30°S of Moon'},
     {when:'UT 2021-11-04 21',text:'New Moon 359850.596 km'},
     {when:'UT 2021-11-05 00',text:'Uranus at opposition'},
     {when:'UT 2021-11-05 22',text:'Moon at perigee 358843.553 km'},
     {when:'UT 2021-11-06 16',text:'Antares 3.91°S of Moon'},
     {when:'UT 2021-11-08 05',text:'Venus 1.10°S of Moon'},
     {when:'UT 2021-11-10 04',text:'Mars 1.06°S of Mercury'},
     {when:'UT 2021-11-10 14',text:'Saturn 4.11°N of Moon'},
     {when:'UT 2021-11-11 13',text:'First Quarter 379165.191 km'},
     {when:'UT 2021-11-11 17',text:'Jupiter 4.36°N of Moon'},
     {when:'UT 2021-11-13 19',text:'Neptune 4.24°N of Moon'},
     {when:'UT 2021-11-18 02',text:'Uranus 1.45°N of Moon'},
     {when:'UT 2021-11-19 09',text:'PARTIAL ECLIPSE Full Moon 405299.694 km'},
     {when:'UT 2021-11-20 13',text:'Aldebaran 6.39°S of Moon'},
     {when:'UT 2021-11-21 02',text:'Moon at apogee 406279.281 km'},
     {when:'UT 2021-11-24 04',text:'Pollux 2.53°N of Moon'},
     {when:'UT 2021-11-26 23',text:'Regulus 5.16°S of Moon'},
     {when:'UT 2021-11-27 12',text:'Last Quarter 388391.853 km'},
     {when:'UT 2021-11-29 05',text:'Mercury in superior conjunction 0.72° South'},
     {when:'UT 2021-11-30 16',text:'Antares 3.75°S of Mercury'},
     {when:'UT 2021-11-30 23',text:'Spica 5.87°S of Moon'},
     {when:'UT 2021-12-03 00',text:'Mars 0.70°S of Moon'},
     {when:'UT 2021-12-04 03',text:'Antares 3.91°S of Moon'},
     {when:'UT 2021-12-04 08',text:'TOTAL ECLIPSE New Moon 356803.846 km'},
     {when:'UT 2021-12-04 10',text:'Moon at perigee 356794.105 km'},
     {when:'UT 2021-12-04 13',text:'Mercury 0.02°N of Moon'},
     {when:'UT 2021-12-07 01',text:'Venus 1.88°N of Moon'},
     {when:'UT 2021-12-08 02',text:'Saturn 4.19°N of Moon'},
     {when:'UT 2021-12-09 06',text:'Jupiter 4.48°N of Moon'},
     {when:'UT 2021-12-11 01',text:'Neptune 4.24°N of Moon'},
     {when:'UT 2021-12-11 02',text:'First Quarter 386531.530 km'},
     {when:'UT 2021-12-15 06',text:'Uranus 1.54°N of Moon'},
     {when:'UT 2021-12-17 19',text:'Aldebaran 6.38°S of Moon'},
     {when:'UT 2021-12-18 02',text:'Moon at apogee 406319.677 km'},
     {when:'UT 2021-12-19 05',text:'Full Moon 405934.567 km'},
     {when:'UT 2021-12-21 10',text:'Pollux 2.58°N of Moon'},
     {when:'UT 2021-12-21 16',text:'Solstice'},
     {when:'UT 2021-12-24 05',text:'Regulus 5.06°S of Moon'},
     {when:'UT 2021-12-26 18',text:'Antares 4.55°S of Mars'},
     {when:'UT 2021-12-27 02',text:'Last Quarter 380818.523 km'},
     {when:'UT 2021-12-28 08',text:'Spica 5.76°S of Moon'},
     {when:'UT 2021-12-29 01',text:'Venus 4.23°N of Mercury'},
     {when:'UT 2021-12-31 14',text:'Antares 3.87°S of Moon'},
     {when:'UT 2021-12-31 20',text:'Mars 0.95°N of Moon'},
     {when:'UT 2022-01-01 23',text:'Moon at perigee 358032.559 km'},
     {when:'UT 2022-01-02 19',text:'New Moon 358676.638 km'},
     {when:'UT 2022-01-03 08',text:'Venus 7.53°N of Moon'},
     {when:'UT 2022-01-04 01',text:'Mercury 3.12°N of Moon'},
     {when:'UT 2022-01-04 07',text:'Earth at perihelion 0.983336540 AU'},
     {when:'UT 2022-01-04 17',text:'Saturn 4.19°N of Moon'},
     {when:'UT 2022-01-06 00',text:'Jupiter 4.45°N of Moon'},
     {when:'UT 2022-01-07 10',text:'Neptune 4.08°N of Moon'},
     {when:'UT 2022-01-07 11',text:'Mercury at greatest elongation 19.2° East'},
     {when:'UT 2022-01-09 01',text:'Venus in inferior conjunction 4.85° North'},
     {when:'UT 2022-01-09 18',text:'First Quarter 394416.641 km'},
     {when:'UT 2022-01-11 11',text:'Uranus 1.45°N of Moon'},
     {when:'UT 2022-01-14 02',text:'Aldebaran 6.46°S of Moon'},
     {when:'UT 2022-01-14 09',text:'Moon at apogee 405804.852 km'},
     {when:'UT 2022-01-17 16',text:'Pollux 2.63°N of Moon'},
     {when:'UT 2022-01-18 00',text:'Full Moon 401023.649 km'},
     {when:'UT 2022-01-20 11',text:'Regulus 4.91°S of Moon'},
     {when:'UT 2022-01-23 10',text:'Mercury in inferior conjunction 3.30° North'},
     {when:'UT 2022-01-24 14',text:'Spica 5.52°S of Moon'},
     {when:'UT 2022-01-25 14',text:'Last Quarter 374709.289 km'},
     {when:'UT 2022-01-27 23',text:'Antares 3.71°S of Moon'},
     {when:'UT 2022-01-29 15',text:'Mars 2.41°N of Moon'},
     {when:'UT 2022-01-30 02',text:'Venus 10.15°N of Moon'},
     {when:'UT 2022-01-30 07',text:'Moon at perigee 362251.771 km'},
     {when:'UT 2022-01-31 00',text:'Mercury 7.57°N of Moon'}
  ];
  
  /* Events within the next n days, in ascending order. The date-time returned uses LT.*/
  var current_events = function(when_start, n_days /*where?, maybe in the future*/){
    var result = [];
    //to avoid all problems with time zones, and related errors, start with the previous day or two
    //this way the UT/LT distinction doesn't matter
    var start = find_calendar_date_from_jd(when_start.jd - 2); 
    var end = find_calendar_date_from_jd(when_start.jd + n_days); 
    //this impl uses simple text compare 
    var start_text = 'UT ' + start.y + '-' + pad(start.m) + '-' + pad(Math.floor(start.d)) + " 00:00"; 
    var end_text = 'UT ' + end.y + '-' + pad(end.m) + '-' + pad(Math.floor(end.d)) + " 2359";
    var fixed_locale = 'en'; 
    var event;
    for(var idx = 0; idx < all_events.length; ++idx){
      var event_time = all_events[idx].when; 
      if (event_time.localeCompare(start_text, fixed_locale) >= 0 && event_time.localeCompare(end_text, fixed_locale) <= 0){
        event = {
          when: when(all_events[idx].when),
          when_text: when(all_events[idx].when).toStringLT().substring(0, 19), /* SHOULD I REALLY HAVE THIS POLICY? LIKELY NOT. */
          text: all_events[idx].text 
        };
        result.push(event);
      }
    }
    return result;
  };

  /* Coords have equinox 2021.5. */
  var build_messiers = function(){
    var messier = [110];
    messier[0]=["M1",1.4651798,0.3844894,8.4,"Tau","NB","!! famous Crab Neb. supernova remnant","Crab Nebula"];
    messier[1]=["M2",5.6487478,0.0159331,6.5,"Aqr","GC","200-mm telescope needed to resolve",""];
    messier[2]=["M3",3.5918444,0.4934997,6.2,"CVn","GC","!! contains many variable stars",""];
    messier[3]=["M4",4.2975271,-0.4639415,5.6,"Sco","GC","bright globular near Antares",""];
    messier[4]=["M5",4.0128997,0.0350122,5.6,"Ser","GC","!! one of the sky's finest globulars",""];
    messier[5]=["M6",4.631679,-0.5624617,4.2,"Sco","OC","!! Butterfly Cluster; best at low power","Butterfly Cluster"];
    messier[6]=["M7",4.6920333,-0.6077146,3.3,"Sco","OC","!! excellent in binocs or rich-field scope","Ptolemy's Cluster"];
    messier[7]=["M8",4.7347244,-0.4255288,6,"Sgr","NB","!! Lagoon Nebula w/open cl. NGC 6530","Lagoon Nebula"];
    messier[8]=["M9",4.5398629,-0.3235411,7.7,"Oph","GC","smallest of Ophiuchus globulars",""];
    messier[9]=["M10",4.4428889,-0.0721197,6.6,"Oph","GC","rich globular cluster; M12 is 3°NW",""];
    messier[10]=["M11",4.940386,-0.1089069,6.3,"Sct","OC","!! Wild Duck Cl.; the best open cluster?","Wild Duck Cluster"];
    messier[11]=["M12",4.3996155,-0.0346816,6.7,"Oph","GC","loose globular cluster near M10",""];
    messier[12]=["M13",4.3740943,0.6357668,5.8,"Her","GC","!! Hercules Cluster; NGC 6207 0.5°NE","Great Hercules Globular"];
    messier[13]=["M14",4.6195769,-0.0569219,7.6,"Oph","GC","200-mm telescope needed to resolve",""];
    messier[14]=["M15",5.6332204,0.2140087,6.2,"Peg","GC","rich, compact globular","Great Pegasus Globular"];
    messier[15]=["M16",4.799738,-0.2403878,6.4,"Ser","NB","Eagle Neb. w/open cl.; use neb. filter","Eagle Nebula"];
    messier[16]=["M17",4.8085576,-0.2822575,7,"Sgr","NB","!! Swan or Omega Nebula; use neb. filter","Omega Nebula"];
    messier[17]=["M18",4.8046684,-0.2988463,7.5,"Sgr","OC","sparse cluster; 1°S of M17",""];
    messier[18]=["M19",4.4677425,-0.4589517,6.8,"Oph","GC","oblate globular; M62 4°S",""];
    messier[19]=["M20",4.7294298,-0.4019779,9,"Sgr","NB","!! Trifid Nebula; look for dark lanes","Trifid Nebula"];
    messier[20]=["M21",4.7381334,-0.3926512,6.5,"Sgr","OC","0.7°NE of M20; sparse cluster",""];
    messier[21]=["M22",4.8769353,-0.4167974,5.1,"Sgr","GC","spectacular from southern latitude","Sagittarius Cluster"];
    messier[22]=["M23",4.7039544,-0.3319268,6.9,"Sgr","OC","bright, loose open cluster",""];
    messier[23]=["M24",4.791634,-0.3227263,4.6,"Sgr","OC","rich star cloud; best in big binoculars","Sagittarius Star Cloud"];
    messier[24]=["M25",4.8558001,-0.335683,6.5,"Sgr","OC","bright but sparse open cluster",""];
    messier[25]=["M26",4.9147578,-0.1636463,8,"Sct","OC","bright, coarse cluster",""];
    messier[26]=["M27",5.2382921,0.3975257,7.4,"Vul","NB","!! Dumbbell Nebula; a superb object","Dumbbell Nebula"];
    messier[27]=["M28",4.8250606,-0.4337763,6.8,"Sgr","GC","compact globular near M22",""];
    messier[28]=["M29",5.3437329,0.6737636,7.1,"Cyg","OC","small, poor open cluster 2°S of γ Cygni",""];
    messier[29]=["M30",5.6793822,-0.402909,7.2,"Cap","GC","toughest in one-night Messier marathon",""];
    messier[30]=["M31",0.1875328,0.7222927,3.4,"And","GY","!! Andromeda Gal.; look for dust lanes","Andromeda Galaxy"];
    messier[31]=["M32",0.1918992,0.7153097,8.1,"And","GY","closest companion to M31",""];
    messier[32]=["M33",0.4150213,0.5368574,5.7,"Tri","GY","large diffuse spiral; requires dark sky","Triangulum Galaxy"];
    messier[33]=["M34",0.7129285,0.7482944,5.5,"Per","OC","best at low power",""];
    messier[34]=["M35",1.6153818,0.4246097,5.3,"Gem","OC","!! look for sm. cluster NGC 2158 0.25°S",""];
    messier[35]=["M36",1.4727303,0.5959501,6.3,"Aur","OC","bright but scattered group; use low pow.",""];
    messier[36]=["M37",1.5437761,0.5681676,6.2,"Aur","OC","!! finest of three Auriga clusters; very rich",""];
    messier[37]=["M38",1.440528,0.6256875,7.4,"Aur","OC","look for small cluster NGC 1907 0.5°S",""];
    messier[38]=["M39",5.6416789,0.8469928,4.6,"Cyg","OC","very sparse cluster; use low power",""];
    messier[39]=["M40",3.2438052,1.0116667,8.4,"UMa","OC","double star Winneke 4; separation 50arcsec","Winnecke 4"];
    messier[40]=["M41",1.7799062,-0.3622945,4.6,"CMa","OC","4°S of Sirius; bright but coarse",""];
    messier[41]=["M42",1.4680687,-0.0949014,4,"Ori","NB","!! Orion Nebula; finest in northern sky","Great Nebula in Orion"];
    messier[42]=["M43",1.468948,-0.0917035,9,"Ori","NB","detached part of Orion Nebula","De Mairan's Nebula"];
    messier[43]=["M44",2.2747517,0.3474272,3.7,"Cnc","OC","!! Beehive or Praesepe; use low power","Beehive Cluster"];
    messier[44]=["M45",0.9960672,0.4220557,1.6,"Tau","OC","!! Pleiades; look for subtle nebulosity","Pleiades"];
    messier[45]=["M46",2.0192913,-0.2595014,6,"Pup","OC","!! contains planetary nebula NGC 2438",""];
    messier[46]=["M47",1.996608,-0.2539315,5.2,"Pup","OC","coarse cluster 1.5°W of M46",""];
    messier[47]=["M48",2.1592393,-0.1023847,5.5,"Hya","OC","former lost Messier; large sparse cl.",""];
    messier[48]=["M49",3.2763893,0.1375555,8.4,"Vir","GY","very bright elliptical",""];
    messier[49]=["M50",1.8510715,-0.1460175,6.3,"Mon","OC","between Sirius & Procyon; use low mag",""];
    messier[50]=["M51",3.5382345,0.821576,8.4,"CVn","GY","!! Whirlpool Galaxy; superb in big scope","Whirlpool Galaxy"];
    messier[51]=["M52",6.1311923,1.0768962,7.3,"Cas","OC","young, rich cl.; faint Bubble Neb. nearby",""];
    messier[52]=["M53",3.4642718,0.3150853,7.6,"Com","GC","150-mm telescope needed to resolve",""];
    messier[53]=["M54",4.958809,-0.531531,7.6,"Sgr","GC","not easily resolved",""];
    messier[54]=["M55",5.1546628,-0.5395818,6.3,"Sgr","GC","bright, loose globular cluster",""];
    messier[55]=["M56",5.0502796,0.5274874,8.3,"Lyr","GC","within a rich starfield",""];
    messier[56]=["M57",4.9497496,0.5770281,8.8,"Lyr","NB","!! Ring Nebula; an amazing smoke ring","Ring Nebula"];
    messier[57]=["M58",3.3108258,0.2041796,9.7,"Vir","GY","bright barred spiral; M59 and M60 1°E",""];
    messier[58]=["M59",3.3295812,0.2012777,9.6,"Vir","GY","bright elliptical paired with M60",""];
    messier[59]=["M60",3.3369965,0.1995352,8.8,"Vir","GY","bright elliptical with M59 and NGC 4647",""];
    messier[60]=["M61",3.2419417,0.0758789,9.7,"Vir","GY","face-on two-armed spiral",""];
    messier[61]=["M62",4.4618075,-0.5261591,6.5,"Oph","GC","asymmetrical; in rich field",""];
    messier[62]=["M63",3.4765265,0.7316457,8.6,"CVn","GY","!! Sunflower Galaxy; bright, elongated","Sunflower Galaxy"];
    messier[63]=["M64",3.3935965,0.3764213,8.5,"Com","GY","!! Black Eye Gal; eye needs big scope","Black Eye Galaxy"];
    messier[64]=["M65",2.9671532,0.2262907,9.3,"Leo","GY","!! bright elongated spiral",""];
    messier[65]=["M66",2.9728222,0.2245434,8.9,"Leo","GY","!! M65 and NGC 3628 in same field",""];
    messier[66]=["M67",2.3194345,0.2048217,6.1,"Cnc","OC","one of the oldest star clusters known",""];
    messier[67]=["M68",3.3189356,-0.4689328,7.8,"Hya","GC","150-mm telescope needed to resolve",""];
    messier[68]=["M69",4.8555152,-0.5643223,7.6,"Sgr","GC","small, poor globular cluster",""];
    messier[69]=["M70",4.9069885,-0.5633436,7.9,"Sgr","GC","small globular 2°E of M69",""];
    messier[70]=["M71",5.2131183,0.3288301,8.2,"Sge","GC","loose globular; looks like an open cluster",""];
    messier[71]=["M72",5.4745692,-0.2173093,9.3,"Aqr","GC","near the Saturn Nebula, NGC 7009",""];
    messier[72]=["M73",5.4985624,-0.2190187,9,"Aqr","OC","group of four stars only; an asterism",""];
    messier[73]=["M74",0.4269855,0.2773748,9.4,"Psc","GY","faint, elusive spiral; tough in small scope",""];
    messier[74]=["M75",5.2681262,-0.3814208,8.5,"Sgr","GC","small and distant; 59,000 ly away",""];
    messier[75]=["M76",0.4527592,0.9018894,10.1,"Per","NB","Little Dumbell; faint but distinct","Little Dumbbell Nebula"];
    messier[76]=["M77",0.7147227,0.0021629,8.9,"Cet","GY","a Seyfert galaxy; with starlike nucleus",""];
    messier[77]=["M78",1.5175741,0.0009889,8.3,"Ori","NB","bright featureless reflection nebula",""];
    messier[78]=["M79",1.4197638,-0.42816,7.7,"Lep","GC","200-mm telescope needed to resolve",""];
    messier[79]=["M80",4.2685749,-0.4020371,7.3,"Sco","GC","very compressed globular",""];
    messier[80]=["M81",2.606399,1.2036478,6.9,"UMa","GY","!! bright spiral visible in binoculars","Bode's Galaxy"];
    messier[81]=["M82",2.6073598,1.2144097,8.4,"UMa","GY","!! the exploding galaxy; M81 0.5°S","Cigar Galaxy"];
    messier[82]=["M83",3.5701398,-0.5231742,7.6,"Hya","GY","large and diffuse; superb from far south","Southern Pinwheel"];
    messier[83]=["M84",3.2558671,0.2227806,9.1,"Vir","GY","!! w/M86 in Markarian's Chain",""];
    messier[84]=["M85",3.2575879,0.3155743,9.1,"Com","GY","bright elliptical shape",""];
    messier[85]=["M86",3.2606641,0.2239452,8.9,"Vir","GY","!! w/many NGC galaxies in Chain",""];
    messier[86]=["M87",3.2807288,0.2143512,8.6,"Vir","GY","famous jet and black hole",""];
    messier[87]=["M88",3.2863874,0.2498413,9.6,"Com","GY","bright multiple-arm spiral",""];
    messier[88]=["M89",3.3020985,0.2169758,9.8,"Vir","GY","elliptical; resembles M87 but smaller",""];
    messier[89]=["M90",3.3068921,0.2277403,9.5,"Vir","GY","bright barred spiral near M89",""];
    messier[90]=["M91",3.3012144,0.2510094,10.2,"Com","GY","some lists say M91=M58, not NGC 4548",""];
    messier[91]=["M92",4.5280877,0.7524329,6.4,"Her","GC","9°NE of M13; fine but often overlooked",""];
    messier[92]=["M93",2.0311781,-0.4174763,6,"Pup","OC","compact, bright cluster; fairly rich",""];
    messier[93]=["M94",3.3680889,0.7158753,8.2,"CVn","GY","very bright and comet-like",""];
    messier[94]=["M95",2.8149273,0.2022265,9.7,"Leo","GY","bright barred spiral",""];
    messier[95]=["M96",2.827141,0.2042546,9.2,"Leo","GY","M95 in same field",""];
    messier[96]=["M97",2.9497543,0.9581723,9.9,"UMa","NB","Owl Nebula; distinct grey oval","Owl Nebula"];
    messier[97]=["M98",3.207016,0.25826,10.1,"Com","GY","nearly edge-on spiral near star 6 Com. B.",""];
    messier[98]=["M99",3.2288222,0.2498276,9.9,"Com","GY","nearly face-on spiral near M98",""];
    messier[99]=["M100",3.2466967,0.2742657,9.3,"Com","GY","face-on spiral with starlike nucleus",""];
    messier[100]=["M101",3.6824695,0.9467937,7.9,"UMa","GY","!! Pinwheel Gal; diffuse face-on spiral","Pinwheel Galaxy"];
    messier[101]=["M102",3.9579296,0.9718791,9.9,"Dra","GY","or is M102=M101? (look for NGC 5907)",""];
    messier[102]=["M103",0.4129564,1.061331,7.4,"Cas","OC","three NGC open clusters nearby",""];
    messier[103]=["M104",3.3210097,-0.2048056,8,"Vir","GY","!! Sombrero Galaxy; look for dust lane","Sombrero Galaxy"];
    messier[104]=["M105",2.8315115,0.2176327,9.3,"Leo","GY","bright elliptical near M95 and M96",""];
    messier[105]=["M106",3.2286762,0.82375,8.4,"CVn","GY","!! superb large, bright spiral",""];
    messier[106]=["M107",4.3358569,-0.2285388,7.9,"Oph","GC","small, faint globular",""];
    messier[107]=["M108",2.9354125,0.9695229,10,"UMa","GY","nearly edge-on; paired with M97 0.75°SE",""];
    messier[108]=["M109",3.1359514,0.9296259,9.8,"UMa","GY","barred spiral near γ UMA",""];
    messier[109]=["M110",0.181418,0.7295672,8.5,"And","GY","more distant companion to M31",""];
    return messier;
  };  
  var messiers = build_messiers();

  var find_deep_sky_object = function(name_or_alias /*case-insensitive*/, target){
    var result;
    var name = name_or_alias.toLowerCase(); 
    for(var idx=0; idx < target.length; ++idx){ 
      if (target[idx][FIXED.name].toLowerCase() === name || target[idx][FIXED.alt_name].toLowerCase() === name){
        result = target[idx];
        break;
      }
    }
    return result;
  };
  
  /* Find a messier object either by 'M1' or 'Crab Nebula' (case-insensitive). */  
  var find_messier = function(name_raw){
    return find_deep_sky_object(name_raw, messiers);
  };
  
  /* Coords have equinox 2021.5. */
  var build_caldwells = function(){
    var caldwell = [42];
    caldwell[0]=["C68",4.988799,-0.6443354,9.7,"CrA","Bn","","","NGC_6729"];
    caldwell[1]=["C69",4.5167247,-0.6479299,12.8,"Sco","Pl","","Bug Nebula","NGC_6302"];
    caldwell[2]=["C70",0.243969,-0.6556699,8.1,"Scl","Sp","","","NGC_300"];
    caldwell[3]=["C71",2.0641366,-0.6738107,5.8,"Pup","Oc","","","NGC_2477"];
    caldwell[4]=["C72",0.0697072,-0.6817938,8.2,"Scl","Sb","","","NGC_55"];
    caldwell[5]=["C73",1.3736071,-0.698592,7.3,"Col","Gc","","","NGC_1851"];
    caldwell[6]=["C74",2.6555634,-0.70754,8.2,"Vel","Pl","","Eight Burst Nebula","NGC_3132"];
    caldwell[7]=["C75",4.3069478,-0.7105974,5.8,"Sco","Oc","","","NGC_6124"];
    caldwell[8]=["C76",4.4310115,-0.7301343,2.6,"Sco","Oc","","","NGC_6231"];
    caldwell[9]=["C77",3.5201817,-0.7527258,7,"Cen","Px","","Centaurus A","Centaurus_A"];
    caldwell[10]=["C78",4.7540984,-0.7626289,6.6,"CrA","Gc","","","NGC_6541"];
    caldwell[11]=["C79",2.6986501,-0.8120095,6.7,"Vel","Gc","","","NGC_3201"];
    caldwell[12]=["C80",3.5259873,-0.8306794,3.6,"Cen","Gc","","Omega Centauri","Omega_Centauri"];
    caldwell[13]=["C81",4.5689921,-0.8453361,8.1,"Ara","Gc","","","NGC_6352"];
    caldwell[14]=["C82",4.3760524,-0.8518353,5.2,"Ara","Oc","","","NGC_6193"];
    caldwell[15]=["C83",3.4324577,-0.8653592,9.5,"Cen","Sp","","","NGC_4945"];
    caldwell[16]=["C84",3.6118379,-0.8983827,7.6,"Cen","Gc","","","NGC_5286"];
    caldwell[17]=["C85",2.2724809,-0.9275345,2.5,"Vel","Oc","","Omicron Vel Cluster","IC_2391"];
    caldwell[18]=["C86",4.6358167,-0.9368278,5.6,"Ara","Gc","","","NGC_6397"];
    caldwell[19]=["C87",0.8416381,-0.9623188,8.4,"Hor","Gc","","","NGC_1261"];
    caldwell[20]=["C88",3.9588911,-0.9718378,7.9,"Cir","Oc","","","NGC_5823"];
    caldwell[21]=["C89",4.2790851,-1.0114302,5.4,"Nor","Oc","","S Norma Cluster","NGC_6087"];
    caldwell[22]=["C90",2.4522175,-1.0194281,9.7,"Car","Pl","","","NGC_2867"];
    caldwell[23]=["C91",2.9117363,-1.0259597,3,"Car","Oc","","","NGC_3532"];
    caldwell[24]=["C92",2.8127443,-1.0468463,6.2,"Car","Bn","","Eta Carinae Nebula","Carina_Nebula"];
    caldwell[25]=["C93",5.0299939,-1.0462625,5.4,"Pav","Gc","","","NGC_6752"];
    caldwell[26]=["C94",3.381137,-1.0550461,4.2,"Cru","Oc","","","Jewel Box"];
    caldwell[27]=["C95",4.2129807,-1.0569321,5.1,"TrA","Oc","","","NGC_6025"];
    caldwell[28]=["C96",2.0885232,-1.0633562,3.8,"Car","Oc","","","NGC_2516"];
    caldwell[29]=["C97",3.0417223,-1.0774919,5.3,"Cen","Oc","","","NGC_3766"];
    caldwell[30]=["C98",3.3317338,-1.1010282,6.9,"Cru","Oc","","","NGC_4609"];
    caldwell[31]=["C99",3.3786105,-1.1015895,undefined,"Cru","Dn","","Coalsack Nebula","Coalsack_Nebula"];
    caldwell[32]=["C100",3.0438884,-1.1022179,4.5,"Cen","Oc","","Lambda Centauri Nebula","IC_2944"];
    caldwell[33]=["C101",5.0258072,-1.1137574,9,"Pav","Sb","","","NGC_6744"];
    caldwell[34]=["C102",2.8098669,-1.125966,1.9,"Car","Oc","","Theta Car Cluster","IC_2602"];
    caldwell[35]=["C103",1.4772201,-1.205828,1,"Dor","Bn","","Tarantula Nebula","Tarantula_Nebula"];
    caldwell[36]=["C104",0.2789282,-1.2345565,6.6,"Tuc","Gc","","","NGC_362"];
    caldwell[37]=["C105",3.4080283,-1.2391647,7.3,"Mus","Gc","","","NGC_4833"];
    caldwell[38]=["C106",0.1092752,-1.2560144,4,"Tuc","Gc","","47 Tucanae","47_Tucanae"];
    caldwell[39]=["C107",4.3121592,-1.2609521,9.3,"Aps","Gc","","","NGC_6101"];
    caldwell[40]=["C108",3.2597478,-1.2703478,7.8,"Mus","Gc","","","NGC_4372"];
    caldwell[41]=["C109",2.6581852,-1.4132399,11.6,"Cha","Pl","","","NGC_3195"];    
    return caldwell;
  };
  var caldwells = build_caldwells();

  /* Find a Caldwell object either by 'C103' or 'Tarantula Nebula' (case-insensitive). */  
  var find_caldwell = function(name_raw){
    return find_deep_sky_object(name_raw, caldwells);
  };
  
  /* This exist in order to avoid creating a large number of identical objects. */
  var when_fixed_equinox = when("J2021.5");
  
  /* Array indices, for arrays of stars and messier-caldwell objects. */
  var FIXED = {
    name: 0,  α: 1, δ: 2, mag: 3, constellation: 4, type: 5, description: 6, alt_name: 7
  };
  
  /* Ephem object for an object with a fixed position (messier object, star). */
  var fixed_ephem = function(thing) {
    //this impl depends on a conventional data structure (array), and indices, to represent all 'fixed' things
    var result = {
      equinox: when_fixed_equinox, 
      name: thing[FIXED.name],      
      α: thing[FIXED.α],
      δ: thing[FIXED.δ],
      mag: thing[FIXED.mag]
    };
    if (thing[FIXED.alt_name]){
      result.alt_name = thing[FIXED.alt_name];
    }
    return result;
  };
  
  /* 
   Return data for all items above a given altitude (for a given location and time).
   The items must share the FIXED structure defined above.
   Return an object: .thing, .ephem, and .idx (the index into the given items array).
   For convenience, the alt is in degrees here.  
   Limited to those objects brighter than the given limiting magnitude. 
  */
  var fixed_ephem_find_visible = function(items, when, where, minimum_alt /*degrees*/, limiting_mag, opts, sort_by_fn /*optional*/){
    var rows = [];
    var options = {where: where, equinox: when};
    for (var i = 0; i < items.length; ++i){
      if (items[i][FIXED.mag] <= limiting_mag){
        var ephem = fixed_ephem(items[i]);
        apply_options(ephem, when, options); //note that these are not the options passed by the caller!
        if (ephem.a > rads(minimum_alt)) {
          apply_options(ephem, when, opts); //the opts passed by the caller
          rows.push({thing: items[i], ephem: ephem, idx: i});
        }
      }
    }
    if (sort_by_fn){
      rows.sort(sort_by_fn);
    }
    return rows;
  };
  
  var highest_first = function(row_a, row_b){
      return row_a.ephem.a > row_b.ephem.a ? -1 : 1;
  };
  
  /* 
   Return ephem's for all Messier objects above a given altitude (for a given location and time).
   For convenience, the alt is in degrees here. 
   Limited to those objects brighter than the given limiting magnitude. 
   Sorts the retured array, with the highest in altitude coming first.
  */
  var find_visible_messiers = function(when, where, minimum_alt /*degrees*/, limiting_mag, opts){
    return fixed_ephem_find_visible(messiers, when, where, minimum_alt, limiting_mag, opts, highest_first);
  };

  /* 
   Return ephem's for all Caldwell objects (declination < -35) above a given altitude (for a given location and time).
   For convenience, the alt is in degrees here. 
   Limited to those objects brighter than the given limiting magnitude. 
   Sorts the retured array, with the highest in altitude coming first.
  */
  var find_visible_caldwells = function(when, where, minimum_alt /*degrees*/, limiting_mag, opts){
    return fixed_ephem_find_visible(caldwells, when, where, minimum_alt, limiting_mag, opts, highest_first);
  };
  
  /* 
   Return ephem's for all stars (in the Yale Bright Star catalog) above a given altitude (for a given location and time).
   Return a star object, its ephem, and the index into my version of the star catalog.
   For convenience, the alt is in degrees here.  
   Limited to those objects brighter than the given limiting magnitude. Not sorted in any specific way. 
  */
  var find_visible_stars = function(when, where, minimum_alt /*degrees*/, limiting_mag, opts){
    return fixed_ephem_find_visible(stars, when, where, minimum_alt, limiting_mag, opts);
  };
    
  /*
   We most often need to interpolate a position, for example right asenscion and declination; that means interpolating on  
   two numeric values in parallel. (Meeus, 1st ed, equation 3.3).
   table: an array of n row-objects, having .when and .thing properties. The whens are evenly spaced, and go forward in time.
   The table must have enough values to ensure that 3-point quadratic interpolation can take place; this means that, 
   at the endpoints, the caller needs to supply an extra data point outside the range of nominal interest.  
   props: the names if the simple numeric properties of thing that are to be interpolated.
   when: the time for which the interpolation is done.
   Returns a new row in the table, with the .thing having the interpolated values of each specified prop.
  */
  var interpolate = function(table, props, when){
    var result = {}, prop;
    if (table.length < 3) {
      console.log('Cannot interpolate. Table has less than 3 rows.');
      return result;
    }
    result.when = when;
    result.thing = {};
    result.thing.equinox = table[0].thing.equinox; //always copy over the same equinox
    var i, j, delta, n, x2 = {delta:1000000 /*a too-large number*/};
    var a, b, c, tabular_interval;
    //find the nearest when
    for (i = 0; i < table.length; ++i){
      delta = table[i].when.jd - when.jd;
      if (Math.abs(delta) <= Math.abs(x2.delta)){
        x2.i = i;
        x2.delta = delta;
      }
    }
    tabular_interval = table[1].when.jd - table[0].when.jd; //assumes even spacing 
    n = (when.jd - table[x2.i].when.jd)/tabular_interval; //in units of the tabular interval; can be either sign; best if abs(n)<= 0.5
    //interpolate each numeric property separately
    for (j=0; j < props.length; ++j){
      prop = props[j];
      a = table[x2.i].thing[prop] - table[x2.i - 1].thing[prop];
      b = table[x2.i + 1].thing[prop] - table[x2.i].thing[prop];
      c = b - a;
      if (a === undefined || b === undefined || c === undefined){
        console.log('Cannot interpolate. Table does not have enough rows. Needs 3 data-points in region of interest. You likely need to add 2 more rows to the table, 1 on each end.');
      }
      else {
        result.thing[prop] = table[x2.i].thing[prop] + (n/2) * (a + b + n*c); 
      }
    }
    return result;
  };
  
  /* Create a table of positions from which values can be interpolated. */
  var interpolation_table = function(when_start, where, thing, bin_gross, num_bin_gross, bin_fine, num_bin_fine, fine_offset, props){
    var i, gross_table = [], fine_table = [], row, when_x, when_fine_start;
    //fill up the gross table with accurate positions, for the mean equinox of date
    for(i = 0; i <= num_bin_gross; i++){
      when_x = when_start.delta(i*bin_gross*SEC_PER_DAY);
      gross_table.push({
        when: when_x,
        thing: position(thing.name, when_x, {equinox: when_x}) //always mean equinox of date; ra, dec; a,A not desired yet
      });
    }
    //proceed to find the fine position, by interpolating from the gross; the interval in the fine table is chosen 
    //such that changes between rows are always sufficient to detect the presence or absence of any phenomenon;
    //the fine table's start-end don't match the gross-table: at the ends of the gross table there aren't enough 
    //points to do the 3-point interpolation; in addition, such points would never be used later in the calc.
    when_fine_start = when_start.delta((bin_gross - fine_offset*bin_fine)*SEC_PER_DAY);
    for(i = 0; i <= num_bin_fine; i++){
      when_x = when_fine_start.delta(i*bin_fine*SEC_PER_DAY);
      //simply re-use the gross table entry data in the fine table, when they match up; otherwise, interpolate
      row = matching_thing_in_gross_table(when_x, gross_table);
      if (!row){
        row = interpolate(gross_table, props, when_x);
      } 
      fine_table.push(row);
    }
    return fine_table;
  };

  /* Return a match to an existing row in the table. */
  var matching_thing_in_gross_table = function(when_x, gross_table){
    var i, result;
    var ONE_MINUTE = 1/(60*24);
    for (i = 0; i < gross_table.length; ++i){
      if (Math.abs(gross_table[i].when.jd - when_x.jd) <= ONE_MINUTE){
        result = gross_table[i];
        break;
      }
    }
    return result;
  };
  
  /* 
   Return an array of .when, .name. 
   The second is controlled by 'names', and provides a tag to distinguish events where the func is increasing/decreasing.
  */
  var find_zeros = function(fine_table, where, zero_func, epsilon, props, event_names){
    var result = []; 
    //find the zero_func for all rows in fine_table, and add it to the row temporarily
    var i, start, end, when, name;
    for (i = 0; i < fine_table.length; ++i){
      fine_table[i].thing.temp = zero_func(fine_table[i], where); 
    }
    //scan for a bin which holds the target; ignore the start-end bins that are adjacent to the region of real interest,
    //since they exist only for interpolation; there's no need or desire to find roots outside the region of real interest.
    for (i = 1; i < fine_table.length-2; ++i){
      //if the zero_func changes sign, then it contains a zero
      start = zero_func(fine_table[i], where);
      end = zero_func(fine_table[i+1], where);
      if (changes_sign_between(start, end)){
        //start a binary search for the time when LHA takes that value
        when = binary_search_with_table(zero_func, fine_table, i, epsilon, where, props);
        name = start < end ? event_names.increasing : event_names.decreasing;
        result.push({
          when: when, 
          name: name
        });              
      }
    }
    for (i = 0; i < fine_table.length; ++i){
      delete fine_table[i].thing.temp;
    }
    return result;
  };
  
  /* 
   Search between i and i+1 for the 'when' corresponding to a 0 value for the function.
   Return the when corresponding to a zero value of the func, within the given epsilon.
   Function values come from interpolating the given table.
   zero_func: the function whose zeroes are to be find. Functin of (row, where).
   fine_table: the table to be used to interpolate (row has .when, .thing)
   start_i: index into the table; the zero needs to come between this index and the next; that is, this index 
   identifies a bin in the table where the zero must be
   epsilon_days: the small interval (in days) below which the binary search will cease
   where: location of observation
   props: the names of .thing properties that are interpolated
  */
  var binary_search_with_table = function(zero_func, fine_table, start_i, epsilon_days, where, props){
    var beginning, start, end, mid; //these are all rows in the table (.when, .thing)
    var width, mid_when, delta, iter_count, result, s_val, m_val;
    iter_count = 0;
    beginning = fine_table[start_i]; 
    start = fine_table[start_i]; //initial values for dynamic pointers, that are reset below
    end = fine_table[start_i+1];
    width = end.when.jd - start.when.jd;
    while (width > epsilon_days && iter_count < 1000 /*avoid infinite loop, just in case*/){
      mid_when = midpoint_when(start, width, beginning);
      mid = interpolate(fine_table, props, mid_when);
      s_val = zero_func(start, where);
      m_val = zero_func(mid, where);
      if (changes_sign_between(s_val, m_val)){
        end = mid;
      }
      else {
        start = mid;
      }
      width = end.when.jd - start.when.jd;
      ++iter_count;
    }
    delta = (mid.when.jd - beginning.when.jd) * SEC_PER_DAY; 
    result = beginning.when.delta(delta);
    if (Math.abs(result.jd -2457703.3604166666) < 0.001){
      console.log("Trapped LT 2016-11-10 15:39");
    }
    return result; 
  };   
  var changes_sign_between = function(a, b){
    var product = Math.sign(a) * Math.sign(b);
    return product <= 0;
  }
  var midpoint_when = function(start, width, beginning){
    var mid_jd = start.when.jd + (width)/2;
    var delta_sec = (mid_jd - beginning.when.jd) * SEC_PER_DAY; 
    result = beginning.when.delta(delta_sec);
    return result;
  };

  // events, rise, set, etc:
    
  /* The start_jd..end_jd interval may not be evenly split by the bin_width. If not, the end has to be a bit beyond end_jd, to make sure we find events near the end. */
  var num_bins = function(start_jd, end_jd, bin_width /*jd*/){
    var is_evenly_divisible = ((end_jd - start_jd) % bin_width === 0);
    var result = Math.trunc((end_jd - start_jd)/bin_width);
    if (!is_evenly_divisible){
      result = result + 1;
    }
    return result;
  };
  /* Can be reused to find different kinds of events. Each bin can have up to 1 target event, not more. */
  var event_detection_table = function(name, start_wh, end_wh, bin_width /*jd*/, opts){
    var i, wh, result = []; //gross bins; change of sign of the zero_func between the start and end of a bin signals a zero exists in the bin
    var n_bins = num_bins(start_wh.jd, end_wh.jd, bin_width);
    for (i = 0; i <= n_bins; ++i){ //note the equal sign here is needed!
      wh = start_wh.delta(i * bin_width * 3600*24); 
      result.push({
        when : wh,
        ephem: position(name, wh, opts)  /* opts.where will add a, A to the ephem; that includes parallax, but no refraction or semi-diameter. */
      });
    }
    return result;
  };
  /* Binary search between start and end. We know there's a zero somewhere in there. Return .when .name .ephem .val */
  var binary_search_for_zeros_core = function(st, en, name, zero_func /*(ephem, when, where)*/, eps /*jd*/, opts, event_names){
    var start = st; //starting values, changed below 
    var end = en;
    var mid, SECS_PER_DAY = 60*60*24;
    var width = 1; //day
    var num_iters = 0;
    var middle_of = function(start, end){
      var result = {};
      result.when = start.when.delta(SECS_PER_DAY*(end.when.jd - start.when.jd)/2);
      result.ephem = position(name, result.when, opts);
      result.val = zero_func(result.ephem, result.when, opts.where);
      return result;
    };
    while (width > eps && num_iters < 500 /*safety first*/){
      mid = middle_of(start, end);
      if (changes_sign_between(st.val, mid.val)){
        end = mid;
      }
      else {
        start = mid;
      }
      width = end.when.jd - start.when.jd;
      ++num_iters;
    }
    //best guess at the zero is the midpoint
    var result = middle_of(start, end);
    //apply the callers names for the events
    var is_increasing = (st.val < end.val);
    result.name = is_increasing ? event_names.increasing : event_names.decreasing; 
    return result;
  };
  /* Scan the table looking for bins with a zero. At most 1 zero per bin. Return [] of .when .name .ephem .val */
  var binary_search_for_zeros = function(name, detection_table, zero_func /*(ephem, when, where)*/, eps /*jd*/, opts, event_names, end_wh){
    var start, end, result = [], possible;
    var point = function(i){
      return {
        when : detection_table[i].when,
        ephem : detection_table[i].ephem, 
        val : zero_func(detection_table[i].ephem, detection_table[i].when, opts.where) //add this, but don't intefere with the table data 
      };
    }
    for(var i = 0; i < detection_table.length - 1; ++i){
      start = point(i); //since these are taken pairwise, the index doesn't go all the way to the end
      end = point(i + 1);
      if (changes_sign_between(start.val, end.val)) {
        possible = binary_search_for_zeros_core(start, end, name, zero_func, eps, opts, event_names);
        if (possible.when.jd <= end_wh.jd){
          result.push(possible);
        }
      }
    }
    return result;
  };
  /* Return [] of .when .name .ephem .val */
  var find_events = function(name, start_wh, end_wh, opts, bin_width /*jd*/, zero_func/*(ephem, when, where)*/, eps /*jd*/, event_names){
    var detection_table = event_detection_table(name, start_wh, end_wh, bin_width, opts); // Nx .ephem .when
    var events = binary_search_for_zeros(name, detection_table, zero_func, eps, opts, event_names, end_wh); // Nx .when .name .ephem .val
    return events;
  }; 
  /* Return .culminations, .horizons, two arrays, possibly empty. Each is an [] of .when .name .ephem .val */
  var rise_culmination_set = function(name, when_start, when_end, opts /*.where, rads*/){
    //I decided not to care about grazes, and detecting if the object ever gets above the horizon: simplifies the algo.
    var result = {culminations: [], horizons: []}; //default empty data
    var EPSILON = 5/(60*60*24); //5 seconds
    var HOUR = 1/24;
    var event_names = {increasing: 'upper', decreasing: 'lower'};
    result.culminations = find_events(name, when_start, when_end, opts, HOUR, zero_culmination, EPSILON, event_names);
    event_names = {increasing: 'rise', decreasing: 'set'};
    result.horizons = find_events(name, when_start, when_end, opts, HOUR, zero_altitude, EPSILON, event_names);
    return result;
  };
  var zero_culmination = function(ephem, when, where){
    return lha_for_culmination(ephem, when, where); 
  };
  var zero_altitude = function(ephem, when, where){
    var result;
    convert_αδ_to_aA(ephem, where, when); //adds a, A, lha; takes parallax into account, but not semi-diam or refraction
    var semi_diam = ephem.size ? ephem.size/2 : 0; 
    var h0 = alt_at_rise_set(semi_diam);    
    result = ephem.a - h0;
    return result;
  };
  /* Return .culminations, .horizons, two arrays, possibly empty. Each array is [] of .when .name .ephem .val */
  var rise_culmination_set_daily = function(name, when, opts) {
    var wh_start = when.startOfDayLT(); //the when param acts as an indicator for a specific day
    var wh_end = when.endOfDayLT();
    return rise_culmination_set(name, wh_start, wh_end, opts); 
  };
  /* 
   Return all twilight phenomena for the given local day (for the Sun). 
   Return [] of .when .name ('rise'|'set') .ephem .val. 
   This is not for the observation window, which always spans two days.
   twilight_degs: -6,-12,-18 for civil, nautical, and astronomical twilight, respectively. 
  */
  var twilight = function(when, opts /*.where, rads*/, twilight_degs) {
    var when_start = when.startOfDayLT(); //the when param acts as an indicator for a specific day
    var when_end = when.endOfDayLT();
    event_names = {increasing: 'rise', decreasing: 'set'};
    var EPSILON = 5/(60*60*24); //5 seconds
    var HOUR = 1/24;
    var result = find_events('sun', when_start, when_end, opts, HOUR, twilight_altitude(twilight_degs), EPSILON, event_names);
    return result;
  };
  /* Returns a function, not a number. Uses closure around the arg. No correction for refraction, or apparent size. */
  var twilight_altitude = function(twilight_degs /* eg -6 */){
    var result = function(ephem, when, where){
      convert_αδ_to_aA(ephem, where, when); //adds a, A, lha; takes parallax into account
      var h0 = rads(twilight_degs);    
      return ephem.a - h0; 
    };
    return result;
  }; 
  /*
   Default: today-sunset..tomorrow-sunrise.
   If the Sun is down, AND it's past midnight, then use yesterday-sunset..today-sunrise.
   This is the usual time of interest to an amateur astronomer.  
   Return an object with .sunset, .sunrise; optional: .twilight_end, .twilight_start.
  */
  var observation_window_memo = []; 
  /* From today-sunset..tomorrow-sunrise. Return whens for: .sunrise .sunset .twilight_end .twilight_start. */
  var observation_window_for_nominal_day = function(when, opts, twilight_degs /*optional*/){
    var RISE = 0, SET = 1;
    var tomorrow = when.next();
    var todays = rise_culmination_set_daily('sun', when, opts);  // return .culminations, .horizons, two [] of .when .name .ephem .val */
    var tomorrows = rise_culmination_set_daily('sun', tomorrow, opts);
    var result = {
      sunset: todays.horizons[SET].when,
      sunrise: tomorrows.horizons[RISE].when
    };
    if (twilight_degs){
      var todays_twil = twilight(when, opts, twilight_degs); //[] of .when .name ('rise'|'set') .ephem .val
      var tomorrows_twil = twilight(tomorrow, opts, twilight_degs);
      if (todays_twil.length > 0 && tomorrows_twil.length > 0){ //no twilight sometimes, if you're far enough north
        result.twilight_end = todays_twil[SET].when; 
        result.twilight_start = tomorrows_twil[RISE].when; 
      }
    }
    return result;
  };
  /* Return an object with .sunset, .sunrise; optional: .twilight_end, .twilight_start. */
  var observation_window = function(when_nominal, opts, twilight_degs /*optional*/) {
    var when, i, memo, result;
    //if it's past local midnight AND the Sun is still below the horizon, then take 'when' as yesterday instead of today
    if (when_nominal.date.getHours() < 12 && position('sun', when_nominal, opts).a < 0){
      when = when_nominal.prev();
    }
    else {
      when = when_nominal;
    }
    //early return iff the result has been computed before
    for (i=0; i < observation_window_memo.length; ++i){
      memo = observation_window_memo[i];
      if (memo.when.jd === when.jd && memo.where.λ === where.λ && memo.where.φ === where.φ){
        if (twilight_degs){ 
          if (twilight_degs === memo.twilight_degs){
            return memo.result;
          }
        }
        else {
          return memo.result;
        }
      } 
    }
    return observation_window_for_nominal_day(when, opts, twilight_degs /*optional*/);
    observation_window_memo.push({when: when, where: where, twilight_degs: twilight_degs, result: result}); //remember the result in case needed again
    return result;
  };
  /* 
   Return only those phenomena in the current observation window.   
   Return .culminations, .horizons, two arrays, possibly empty. Each is an [] of .when .name .ephem .val 
  */
  var rise_culmination_set_observation_window = function(name, when, opts) {
    var obs_window = observation_window(when, opts);
    var result = rise_culmination_set(name, obs_window.sunset, obs_window.sunrise, opts);
    return result;
  }
  /* The 0s of this function correspond to both upper and lower culmination. */
  var lha_for_culmination = function(thing, when, where){
     var result = find_local_hour_angle(thing, when, where); //0..2pi
     return decircularize(result);
  };
  /* 
   Time and right ascenscion have discontinuities. This removes them, but at the same time degrades the data by 
   mapping 2 values of x to the same output value. This is done in order to find zeroes, instead of multiples of pi.
   The 0s correspond to 0/2pi and pi, which is what's needed for culmination.
   This is appropriate only for finding culminations. 
  */
  var decircularize = function(rads /*0..2pi*/){
    //this is a kind of see-saw function: 0(up)..pi/2(down)..0..-pi/2(up)..0
    var result = rads;
    if (0.5*Math.PI < rads && rads <= 1.5*Math.PI){ 
      result = Math.PI - rads;
    }
    else if (1.5*Math.PI < rads && rads <= 2*Math.PI){
      result = rads - 2*Math.PI;
    }
    return result;
  };
    var alt_at_rise_set = function(semi_diam /*rads*/){
    var result;
    var REFRACTION_AT_HORIZON = rads(34/60);
    var result = -1*(semi_diam + REFRACTION_AT_HORIZON);
    //note: parallax in altitude is always calculated anyway (for aA), so it's not added in here
    //parallax has the opposite sign to the other 2 effects, [refraction + semi-diameter]
    return result;
  };

  /* Coords have equinox 2021.5. */
  var build_stars = function(){
    /* Source: Yale Bright Star Catalog r5.  Name, Right Ascension, Declination (J2019.5), and Magnitude.*/
    var ybs = [9096];
    //removal of leading whitespace cuts down on the overall size of this js file
ybs[0]=['',0.0273974,0.7914863,6.7];
ybs[1]=['',0.0269007,-0.0066915,6.29];
ybs[2]=['33 Psc',0.0280811,-0.0975263,4.61];
ybs[3]=['86 Peg',0.0296926,0.2358945,5.51];
ybs[4]=['',0.0322531,1.0220004,5.96];
ybs[5]=['',0.0322978,-0.8544321,5.7];
ybs[6]=['10 Cas',0.0330476,1.1225216,5.59];
ybs[7]=['',0.0337005,0.5086069,6.13];
ybs[8]=['',0.0346028,-0.401214,6.18];
ybs[9]=['',0.0366526,-0.3013618,6.19];
ybs[10]=['',0.0385551,-0.0423988,6.43];
ybs[11]=['',0.0387234,-0.3907665,5.94];
ybs[12]=['',0.0399173,-0.5831116,5.68];
ybs[13]=['',0.0405913,-0.0406342,6.07];
ybs[14]=['α And',0.0414545,0.5098135,2.06];
ybs[15]=['',0.0409676,-0.1519184,5.99];
ybs[16]=['',0.042759,0.6413433,6.19];
ybs[17]=['',0.0421174,-0.3046978,6.06];
ybs[18]=['',0.0435518,0.4464966,6.23];
ybs[19]=['',0.0460494,1.3933715,6.01];
ybs[20]=['β Cas',0.0450051,1.0344446,2.27];
ybs[21]=['87 Peg',0.0442814,0.3199457,5.53];
ybs[22]=['',0.0441326,-0.9404244,6.33];
ybs[23]=['κ1 Scl',0.0455574,-0.4863918,5.42];
ybs[24]=['ε Phe',0.0457813,-0.7963574,3.88];
ybs[25]=['34 Psc',0.0486277,0.1966135,5.51];
ybs[26]=['22 And',0.0499481,0.8061987,5.03];
ybs[27]=['',0.0507578,0.9998139,6.74];
ybs[28]=['',0.0497995,-0.0895188,5.84];
ybs[29]=['γ3 Oct',0.0478983,-1.4329907,5.28];
ybs[30]=['',0.051531,-0.2174758,5.85];
ybs[31]=['',0.0509125,-1.275921,6.64];
ybs[32]=['6 Cet',0.0539312,-0.2678822,4.89];
ybs[33]=['κ2 Scl',0.0552482,-0.4831105,5.41];
ybs[34]=['θ Scl',0.0559259,-0.6111014,5.25];
ybs[35]=['',0.0572304,0.8425057,6.16];
ybs[36]=['',0.0578578,-0.3109971,5.25];
ybs[37]=['',0.0609277,0.6599583,6.73];
ybs[38]=['γ Peg',0.0625982,0.2670893,2.83];
ybs[39]=['',0.0633415,0.4731011,6.3];
ybs[40]=['23 And',0.0638829,0.7182859,5.72];
ybs[41]=['',0.0645297,-0.4520835,5.94];
ybs[42]=['',0.0646815,-0.4566699,6.31];
ybs[43]=['',0.0661491,0.5816409,6.25];
ybs[44]=['χ Peg',0.0685782,0.3547574,4.8];
ybs[45]=['',0.0678832,-0.1337117,5.12];
ybs[46]=['',0.0615328,-1.4813427,5.77];
ybs[47]=['7 Cet',0.0686399,-0.3283548,4.44];
ybs[48]=['',0.0700323,0.3910164,6.24];
ybs[49]=['35 Psc',0.0701928,0.1560369,5.79];
ybs[50]=['',0.0698344,-0.1649389,5.75];
ybs[51]=['',0.0708548,0.5524883,6.45];
ybs[52]=['',0.0711029,0.4782633,6.35];
ybs[53]=['',0.0700286,-0.6071132,6.17];
ybs[54]=['',0.076305,1.3451288,6.35];
ybs[55]=['',0.0763392,0.7629548,6.15];
ybs[56]=['',0.0751758,-0.5467594,5.67];
ybs[57]=['',0.0736795,-1.32282,6.49];
ybs[58]=['36 Psc',0.0771238,0.1458984,6.11];
ybs[59]=['',0.0790689,1.0760423,5.74];
ybs[60]=['',0.0776543,-0.3506575,6.47];
ybs[61]=['',0.0798257,0.8389246,5.89];
ybs[62]=['θ And',0.0795136,0.6772053,4.61];
ybs[63]=['',0.0773958,-1.3728969,6.77];
ybs[64]=['',0.0823209,0.8997586,6.14];
ybs[65]=['',0.0812986,-0.330422,6.45];
ybs[66]=['',0.0824585,0.0315591,6.17];
ybs[67]=['σ And',0.0849094,0.6441062,4.52];
ybs[68]=['',0.0846329,0.1976607,6.05];
ybs[69]=['26 And',0.086578,0.7663808,6.11];
ybs[70]=['',0.0862406,0.5521611,5.87];
ybs[71]=['',0.0863631,-0.1384658,6.46];
ybs[72]=['',0.0862815,-0.7525163,6.33];
ybs[73]=['ι Cet',0.0895521,-0.1519247,3.56];
ybs[74]=['',0.0908956,0.7129487,6.33];
ybs[75]=['',0.0926685,0.8549406,6.52];
ybs[76]=['ζ Tuc',0.0919889,-1.1301968,4.23];
ybs[77]=['',0.0939637,0.5420125,5.9];
ybs[78]=['',0.0955093,0.5764922,5.79];
ybs[79]=['41 Psc',0.0947132,0.1450275,5.37];
ybs[80]=['',0.096076,0.1936637,6.56];
ybs[81]=['ρ And',0.0971232,0.664757,5.18];
ybs[82]=['π Tuc',0.0943942,-1.2131053,5.51];
ybs[83]=['ι Scl',0.0985959,-0.5037461,5.18];
ybs[84]=['',0.099731,-0.3479951,5.12];
ybs[85]=['42 Psc',0.102706,0.2373926,6.23];
ybs[86]=['',0.0976273,-1.3492756,5.97];
ybs[87]=['9 Cet',0.1045223,-0.2110168,6.39];
ybs[88]=['',0.1059511,-0.5396044,6.55];
ybs[89]=['',0.1098522,0.6753766,7.39];
ybs[90]=['',0.1109524,0.9099971,5.57];
ybs[91]=['12 Cas',0.1134146,1.0812328,5.4];
ybs[92]=['',0.1116791,-0.0366552,6.07];
ybs[93]=['',0.1146679,0.9279198,5.74];
ybs[94]=['44 Psc',0.1156593,0.0359302,5.77];
ybs[95]=['β Hyi',0.1161231,-1.3462641,2.8];
ybs[96]=['α Phe',0.1192693,-0.7363061,2.39];
ybs[97]=['κ Phe',0.11891,-0.7602849,3.94];
ybs[98]=['10 Cet',0.1209743,0.0012066,6.19];
ybs[99]=['',0.1235664,-0.4438093,5.98];
ybs[100]=['47 Psc',0.1272763,0.3143655,5.06];
ybs[101]=['',0.1282344,0.7769018,5.17];
ybs[102]=['η Scl',0.1265007,-0.5740117,4.81];
ybs[103]=['48 Psc',0.1279822,0.2890921,6.06];
ybs[104]=['',0.1284904,0.1799167,6.04];
ybs[105]=['',0.1284186,-0.3528402,6.43];
ybs[106]=['',0.1286812,-0.6945758,5.43];
ybs[107]=['',0.1312991,0.6460983,6.26];
ybs[108]=['',0.1297934,-0.8798913,6.26];
ybs[109]=['',0.1409595,1.3463122,6.21];
ybs[110]=['',0.1376427,1.0488701,5.94];
ybs[111]=['28 And',0.1363985,0.5213349,5.23];
ybs[112]=['',0.1350456,-0.2573579,6.14];
ybs[113]=['',0.1347279,-0.5584708,6.57];
ybs[114]=['12 Cet',0.1358633,-0.0669961,5.72];
ybs[115]=['',0.1372351,-0.4131049,5.19];
ybs[116]=['',0.1374857,-0.7124581,6.19];
ybs[117]=['',0.1372917,-0.8394404,5.69];
ybs[118]=['13 Cas',0.1425832,1.163052,6.18];
ybs[119]=['',0.1421263,0.5881795,5.87];
ybs[120]=['λ Cas',0.1438589,0.9536606,4.73];
ybs[121]=['',0.1434551,0.9242907,5.6];
ybs[122]=['λ1 Phe',0.1415584,-0.8497148,4.77];
ybs[123]=['β1 Tuc',0.1418804,-1.0967567,4.37];
ybs[124]=['β2 Tuc',0.1419455,-1.0968924,4.54];
ybs[125]=['',0.1466687,0.7611935,6.7];
ybs[126]=['',0.1510995,1.24093,6.42];
ybs[127]=['κ Cas',0.1493968,1.1004315,4.16];
ybs[128]=['52 Psc',0.1471278,0.3562722,5.38];
ybs[129]=['51 Psc',0.146202,0.1234648,5.67];
ybs[130]=['',0.147101,0.4834388,6.67];
ybs[131]=['',0.1481686,0.495651,6.3];
ybs[132]=['',0.1499914,0.960165,5.93];
ybs[133]=['β3 Tuc',0.1470283,-1.0980332,5.09];
ybs[134]=['16 Cas',0.1557148,1.1670768,6.48];
ybs[135]=['',0.1516108,-0.5138243,5.55];
ybs[136]=['θ Tuc',0.1495891,-1.2417619,6.13];
ybs[137]=['',0.1547718,-0.9120175,5.57];
ybs[138]=['',0.1572595,0.235434,6.4];
ybs[139]=['13 Cet',0.1585879,-0.0606421,5.2];
ybs[140]=['14 Cet',0.1599071,-0.0067604,5.93];
ybs[141]=['',0.1629541,0.9474829,5.08];
ybs[142]=['',0.1615803,0.2325625,6.41];
ybs[143]=['',0.1644649,1.0549511,5.79];
ybs[144]=['λ2 Phe',0.1601491,-0.8357095,5.51];
ybs[145]=['',0.1594899,-0.9547601,6.06];
ybs[146]=['',0.1635152,0.4777467,6.5];
ybs[147]=['',0.1620174,-0.2592763,6.45];
ybs[148]=['',0.1622498,-0.3966144,6.06];
ybs[149]=['',0.1656104,0.7785341,5.13];
ybs[150]=['ζ Cas',0.166596,0.9427402,3.66];
ybs[151]=['π And',0.1659621,0.5905766,4.36];
ybs[152]=['53 Psc',0.1654202,0.2679041,5.89];
ybs[153]=['',0.1669275,0.4211872,6.47];
ybs[154]=['',0.1680306,0.6198974,5.48];
ybs[155]=['',0.1812382,1.4418463,6.4];
ybs[156]=['',0.1675982,-0.430209,5.57];
ybs[157]=['',0.1638834,-1.134579,6.42];
ybs[158]=['',0.1684884,0.0567813,6.39];
ybs[159]=['',0.1670894,-0.9472966,6.41];
ybs[160]=['ε And',0.1732357,0.5136438,4.37];
ybs[161]=['',0.1761187,0.8634553,5.43];
ybs[162]=['δ And',0.1766271,0.5406807,3.27];
ybs[163]=['54 Psc',0.1767042,0.3729497,5.87];
ybs[164]=['55 Psc',0.1791656,0.3762261,5.36];
ybs[165]=['α Cas',0.1821244,0.9888163,2.23];
ybs[166]=['',0.1724156,-1.2744265,6.85];
ybs[167]=['',0.1789415,-0.5906863,6.69];
ybs[168]=['',0.1783964,-0.7797925,6.01];
ybs[169]=['',0.181311,-0.2862194,6.49];
ybs[170]=['',0.1815697,-0.4134103,6.14];
ybs[171]=['',0.1823963,-0.0739005,5.91];
ybs[172]=['32 And',0.1845397,0.6907373,5.33];
ybs[173]=['',0.1806336,-1.0356199,5.89];
ybs[174]=['',0.1891927,1.1565446,5.83];
ybs[175]=['',0.186498,0.4319139,6.04];
ybs[176]=['ξ Cas',0.1888216,0.8836624,4.8];
ybs[177]=['μ Phe',0.1847366,-0.8022806,4.59];
ybs[178]=['',0.1909751,1.0274913,6.17];
ybs[179]=['ξ Phe',0.1865013,-0.9840864,5.7];
ybs[180]=['π Cas',0.1949039,0.8227868,4.94];
ybs[181]=['λ1 Scl',0.1908769,-0.6692598,6.06];
ybs[182]=['',0.1904243,-1.0497269,5.98];
ybs[183]=['ρ Tuc',0.1892829,-1.1405805,5.39];
ybs[184]=['β Cet',0.1948761,-0.3118761,2.04];
ybs[185]=['',0.199166,0.8374362,5.67];
ybs[186]=['',0.1959965,-0.207593,6.02];
ybs[187]=['η Phe',0.1933516,-1.0008685,4.36];
ybs[188]=['21 Cas',0.205566,1.310835,5.66];
ybs[189]=['ο Cas',0.2004191,0.8447709,4.54];
ybs[190]=['φ1 Cet',0.1975478,-0.1831203,4.76];
ybs[191]=['λ2 Scl',0.1973537,-0.6685351,5.9];
ybs[192]=['',0.2030083,0.9658472,5.42];
ybs[193]=['',0.199858,-0.3820306,5.24];
ybs[194]=['',0.2005679,-0.7428002,5.94];
ybs[195]=['',0.1983623,-1.088743,6.07];
ybs[196]=['',0.2094958,1.2119942,6.33];
ybs[197]=['',0.2028769,-0.078747,6.15];
ybs[198]=['',0.2005977,-0.9354555,6.15];
ybs[199]=['18 Cet',0.20315,-0.2227659,6.15];
ybs[200]=['',0.2072333,0.9673048,6.52];
ybs[201]=['',0.2067359,0.7850246,6.05];
ybs[202]=['',0.2040673,-0.2846091,6.47];
ybs[203]=['',0.2093357,1.0418149,6.39];
ybs[204]=['23 Cas',0.2148639,1.3083778,5.41];
ybs[205]=['',0.2040168,-0.8278864,5.8];
ybs[206]=['',0.2062042,-0.3910363,5.5];
ybs[207]=['57 Psc',0.2080393,0.2721445,5.38];
ybs[208]=['',0.2163316,1.2704598,5.87];
ybs[209]=['58 Psc',0.2100852,0.211028,5.5];
ybs[210]=['59 Psc',0.2110276,0.3437599,6.13];
ybs[211]=['ζ And',0.2115565,0.4255865,4.06];
ybs[212]=['60 Psc',0.2116522,0.1196932,5.99];
ybs[213]=['61 Psc',0.2140377,0.3672575,6.54];
ybs[214]=['',0.2128913,-0.3131877,5.7];
ybs[215]=['η Cas',0.2197636,1.0111168,3.44];
ybs[216]=['',0.2141531,-0.3770867,5.57];
ybs[217]=['62 Psc',0.2155701,0.1294508,5.93];
ybs[218]=['',0.2159616,0.0942047,5.75];
ybs[219]=['ν Cas',0.2184438,0.8916059,4.89];
ybs[220]=['δ Psc',0.217289,0.1344243,4.43];
ybs[221]=['64 Psc',0.2186534,0.297709,5.07];
ybs[222]=['ν And',0.2225579,0.7190007,4.53];
ybs[223]=['',0.2203641,-0.2346513,5.59];
ybs[224]=['',0.2194207,-0.4192195,5.9];
ybs[225]=['',0.2178974,-0.8129893,6.27];
ybs[226]=['65 Psc',0.2226909,0.485684,7];
ybs[227]=['65 Psc',0.2227199,0.4856743,7.1];
ybs[228]=['',0.2208447,-0.4056986,6.28];
ybs[229]=['',0.2271103,1.1233672,5.39];
ybs[230]=['',0.2247663,0.7874747,6.15];
ybs[231]=['φ2 Cet',0.2234413,-0.1837424,5.19];
ybs[232]=['λ Hyi',0.2151847,-1.3056173,5.07];
ybs[233]=['',0.2294077,1.080751,6.07];
ybs[234]=['',0.2277359,0.9010216,6.39];
ybs[235]=['',0.2228119,-0.7553423,6.48];
ybs[236]=['',0.2489125,1.4629937,5.62];
ybs[237]=['',0.2303839,0.9021209,6.21];
ybs[238]=['ρ Phe',0.2254006,-0.8878529,5.22];
ybs[239]=['',0.2286966,0.0611152,6.37];
ybs[240]=['',0.2372443,1.06885,4.82];
ybs[241]=['',0.230675,-0.7608341,6.9];
ybs[242]=['',0.2359701,0.6748327,6.69];
ybs[243]=['',0.2344399,-0.4169477,5.46];
ybs[244]=['20 Cet',0.2360914,-0.0179372,4.77];
ybs[245]=['',0.2384893,0.6550995,6.06];
ybs[246]=['',0.2401719,0.9216299,6.27];
ybs[247]=['',0.2367426,-0.4304074,6.46];
ybs[248]=['λ1 Tuc',0.2321955,-1.2110475,6.22];
ybs[249]=['υ1 Cas',0.2456351,1.0312969,4.83];
ybs[250]=['66 Psc',0.2431612,0.3369285,5.74];
ybs[251]=['21 Cet',0.2416319,-0.1505267,6.16];
ybs[252]=['',0.245742,0.8516297,6.27];
ybs[253]=['',0.2378696,-1.0952815,5.7];
ybs[254]=['36 And',0.2448653,0.4144202,5.47];
ybs[255]=['',0.2460907,0.4306269,6.2];
ybs[256]=['',0.250915,1.0142579,6.21];
ybs[257]=['',0.25454,1.202393,6.37];
ybs[258]=['67 Psc',0.2493081,0.4769202,6.09];
ybs[259]=['',0.2478094,-0.1262067,5.85];
ybs[260]=['γ Cas',0.2531698,1.0617297,2.47];
ybs[261]=['υ2 Cas',0.2529171,1.0349293,4.63];
ybs[262]=['',0.2534922,1.055553,5.55];
ybs[263]=['φ3 Cet',0.2491621,-0.1946146,5.31];
ybs[264]=['',0.2485593,-0.4827488,6.1];
ybs[265]=['μ And',0.2528535,0.6739661,3.87];
ybs[266]=['λ2 Tuc',0.2434773,-1.2114458,5.45];
ybs[267]=['η And',0.2546456,0.4107355,4.42];
ybs[268]=['',0.2569468,0.802076,6.12];
ybs[269]=['',0.2613647,1.1600845,5.97];
ybs[270]=['68 Psc',0.2574608,0.5080313,5.42];
ybs[271]=['',0.2592713,0.5945744,5.98];
ybs[272]=['',0.2576099,0.2410589,6.32];
ybs[273]=['',0.259464,0.3755985,6.37];
ybs[274]=['',0.2704679,1.2409029,6.39];
ybs[275]=['φ4 Cet',0.2609663,-0.1965989,5.61];
ybs[276]=['α Scl',0.2602288,-0.510365,4.31];
ybs[277]=['',0.2585672,-1.0573311,6.23];
ybs[278]=['',0.2673969,0.7823724,6.84];
ybs[279]=['',0.2674116,0.7824112,6.04];
ybs[280]=['',0.2659209,0.1151677,6.11];
ybs[281]=['',0.3145728,1.5042433,4.25];
ybs[282]=['',0.4689787,1.5552764,6.46];
ybs[283]=['',0.2773741,0.8927396,6.47];
ybs[284]=['ξ Scl',0.2718531,-0.6772104,5.59];
ybs[285]=['',0.2804301,0.8288782,6.45];
ybs[286]=['39 And',0.2797868,0.7236158,5.98];
ybs[287]=['σ Psc',0.2792595,0.5571019,5.5];
ybs[288]=['',0.2834373,1.0679672,5.92];
ybs[289]=['σ Scl',0.2769068,-0.5486745,5.5];
ybs[290]=['ε Psc',0.27953,0.1397159,4.28];
ybs[291]=['ω Phe',0.2746006,-0.9928694,6.11];
ybs[292]=['25 Cet',0.2798386,-0.0824066,5.43];
ybs[293]=['',0.286569,1.0767842,5.84];
ybs[294]=['',0.2849951,0.9183431,5.99];
ybs[295]=['',0.2783151,-0.8077792,5.36];
ybs[296]=['',0.280653,-0.5133143,6.29];
ybs[297]=['26 Cet',0.2832752,0.02586,6.04];
ybs[298]=['',0.2881923,0.8922971,6.54];
ybs[299]=['',0.2864018,0.5196459,6.19];
ybs[300]=['',0.277213,-1.1404143,6.21];
ybs[301]=['',0.2872008,0.6999817,6.72];
ybs[302]=['',0.3517431,1.5210501,6.25];
ybs[303]=['73 Psc',0.2879449,0.1007271,6];
ybs[304]=['72 Psc',0.288975,0.2628627,5.68];
ybs[305]=['',0.2956191,1.097398,6.54];
ybs[306]=['ψ1 Psc',0.2916334,0.3767828,5.34];
ybs[307]=['ψ1 Psc',0.2916915,0.3766422,5.56];
ybs[308]=['',0.3103681,1.3984589,6.29];
ybs[309]=['77 Psc',0.2920538,0.0876687,6.35];
ybs[310]=['77 Psc',0.2922139,0.087688,7.25];
ybs[311]=['27 Cet',0.2909963,-0.1721666,6.12];
ybs[312]=['',0.298101,0.9957019,6.43];
ybs[313]=['28 Cet',0.293055,-0.1697293,5.58];
ybs[314]=['',0.2986664,0.9357204,6.38];
ybs[315]=['75 Psc',0.2953703,0.2281269,6.12];
ybs[316]=['',0.2930815,-0.4167467,6.14];
ybs[317]=['μ Cas',0.3035898,0.960535,5.17];
ybs[318]=['β Phe',0.292517,-0.813392,3.31];
ybs[319]=['',0.2942838,-0.6203984,6.61];
ybs[320]=['41 And',0.3021746,0.7689277,5.03];
ybs[321]=['',0.2978332,-0.4168175,6.37];
ybs[322]=['',0.3049427,1.0188863,5.79];
ybs[323]=['78 Psc',0.3019946,0.5607149,6.25];
ybs[324]=['ψ2 Psc',0.3015442,0.3639631,5.55];
ybs[325]=['30 Cet',0.3004053,-0.1687932,5.82];
ybs[326]=['80 Psc',0.3031903,0.1006016,5.52];
ybs[327]=['υ Phe',0.3000928,-0.7220868,5.21];
ybs[328]=['ι Tuc',0.2973798,-1.0761835,5.37];
ybs[329]=['',0.3238045,1.392555,5.64];
ybs[330]=['η Cet',0.3039776,-0.1757185,3.45];
ybs[331]=['φ And',0.3087556,0.8265195,4.25];
ybs[332]=['31 Cas',0.3147526,1.2024017,5.29];
ybs[333]=['β And',0.309523,0.6236875,2.06];
ybs[334]=['ζ Phe',0.3023051,-0.9622261,3.92];
ybs[335]=['ψ3 Psc',0.3096816,0.3450988,5.55];
ybs[336]=['44 And',0.3121836,0.7364486,5.65];
ybs[337]=['',0.3119551,0.4463119,5.8];
ybs[338]=['',0.317807,1.1225364,5.55];
ybs[339]=['θ Cas',0.315981,0.9645317,4.33];
ybs[340]=['',0.3112558,0.2755561,6.06];
ybs[341]=['32 Cas',0.3190114,1.1367794,5.57];
ybs[342]=['32 Cet',0.3110146,-0.1534506,6.4];
ybs[343]=['33 Cet',0.3127121,0.0446723,5.95];
ybs[344]=['45 And',0.315852,0.6603984,5.81];
ybs[345]=['82 Psc',0.3154924,0.5504526,5.16];
ybs[346]=['',0.3097799,-1.0049671,6.41];
ybs[347]=['χ Psc',0.3168312,0.369112,4.66];
ybs[348]=['τ Psc',0.31786,0.5271508,4.51];
ybs[349]=['34 Cet',0.3177427,-0.0373032,5.94];
ybs[350]=['',0.3252826,1.0789516,6.41];
ybs[351]=['',0.3221049,0.793277,6.11];
ybs[352]=['',0.3236768,0.526701,6.19];
ybs[353]=['',0.3425709,1.3966634,6.26];
ybs[354]=['',0.3202805,-0.5356159,6.52];
ybs[355]=['',0.3217597,-0.6587353,5.92];
ybs[356]=['φ Psc',0.3269013,0.4310452,4.65];
ybs[357]=['ζ Psc',0.3266123,0.1341938,5.24];
ybs[358]=['ζ Psc',0.3267142,0.1342471,6.3];
ybs[359]=['',0.328422,0.4999168,6.43];
ybs[360]=['87 Psc',0.3284409,0.2835637,5.98];
ybs[361]=['',0.3393953,1.2541394,7.83];
ybs[362]=['37 Cet',0.3293466,-0.136305,5.13];
ybs[363]=['88 Psc',0.3308607,0.124068,6.03];
ybs[364]=['38 Cet',0.3312611,-0.0150205,5.7];
ybs[365]=['',0.3389714,0.8411652,6.61];
ybs[366]=['ν Phe',0.332176,-0.7926963,4.96];
ybs[367]=['',0.3382375,0.5799334,6.02];
ybs[368]=['',0.3418488,0.7856569,6.34];
ybs[369]=['39 Cet',0.339031,-0.0416663,5.41];
ybs[370]=['',0.342969,0.5560192,6.73];
ybs[371]=['',0.358598,1.3558208,6.31];
ybs[372]=['',0.3466571,0.8295971,6.25];
ybs[373]=['κ Tuc',0.3336519,-1.20014,4.86];
ybs[374]=['89 Psc',0.3443195,0.0650522,5.16];
ybs[375]=['',0.3491083,0.6544757,6.46];
ybs[376]=['',0.3394652,-1.1568937,6.24];
ybs[377]=['',0.365566,1.3325736,6.38];
ybs[378]=['φ Cas',0.3553972,1.018295,4.98];
ybs[379]=['υ Psc',0.3519162,0.4778124,4.76];
ybs[380]=['35 Cas',0.3601608,1.1304582,6.34];
ybs[381]=['42 Cet',0.3530174,-0.0069198,5.87];
ybs[382]=['',0.3741749,1.3759728,6.07];
ybs[383]=['',0.3563425,-0.0547103,6.23];
ybs[384]=['',0.3557526,-0.1941956,6.15];
ybs[385]=['91 Psc',0.3591762,0.5035314,5.23];
ybs[386]=['ξ And',0.3648393,0.7965827,4.88];
ybs[387]=['',0.3697428,1.016738,6.45];
ybs[388]=['',0.3653144,0.0320842,6.2];
ybs[389]=['43 Cet',0.3651262,-0.0058959,6.49];
ybs[390]=['',0.3645633,-0.3310796,6.35];
ybs[391]=['47 And',0.3704984,0.6602003,5.58];
ybs[392]=['',0.3702031,0.5996521,6.29];
ybs[393]=['',0.369055,0.3591998,5.97];
ybs[394]=['',0.3813146,1.2407765,6.49];
ybs[395]=['ψ Cas',0.3816908,1.1910342,4.74];
ybs[396]=['',0.3687703,-0.5381515,5.84];
ybs[397]=['44 Cet',0.3714037,-0.1378088,6.21];
ybs[398]=['θ Cet',0.3713215,-0.1408776,3.6];
ybs[399]=['δ Cas',0.3806035,1.0532459,2.68];
ybs[400]=['',0.3727343,-0.1187373,5.91];
ybs[401]=['',0.3740097,-0.271377,6.14];
ybs[402]=['',0.3748317,-0.0477717,6.15];
ybs[403]=['',0.3786214,0.4122991,6.18];
ybs[404]=['',0.3736247,-0.7222342,5.42];
ybs[405]=['',0.382141,0.7604219,5.96];
ybs[406]=['',0.3812199,0.6054712,6.31];
ybs[407]=['',0.3736296,-0.7752195,6.26];
ybs[408]=['46 Cet',0.3781968,-0.2528555,4.9];
ybs[409]=['ρ Psc',0.381436,0.3365593,5.38];
ybs[410]=['94 Psc',0.3833582,0.3377456,5.5];
ybs[411]=['',0.3854029,0.6019385,6.27];
ybs[412]=['',0.3820345,-0.0050167,6.41];
ybs[413]=['ω And',0.3880804,0.7944318,4.83];
ybs[414]=['',0.3870362,0.7192767,6.46];
ybs[415]=['',0.3839933,0.0636411,6.58];
ybs[416]=['',0.3744937,-1.1215132,5.93];
ybs[417]=['47 Cet',0.3836268,-0.2259428,5.66];
ybs[418]=['',0.3885039,0.7059239,6.6];
ybs[419]=['',0.3837836,-0.5660447,5.79];
ybs[420]=['α UMi',0.7817861,1.5593953,2.02];
ybs[421]=['',0.3876556,-0.1883341,6.13];
ybs[422]=['',0.3905572,0.1408861,6.2];
ybs[423]=['38 Cas',0.4051573,1.2282737,5.81];
ybs[424]=['',0.4031502,1.155553,6.14];
ybs[425]=['γ Phe',0.38963,-0.7541134,3.41];
ybs[426]=['49 And',0.3988241,0.8223582,5.27];
ybs[427]=['',0.3914005,-0.5873534,6.58];
ybs[428]=['97 Psc',0.3972583,0.3222934,6.02];
ybs[429]=['48 Cet',0.3954526,-0.3755753,5.12];
ybs[430]=['μ Psc',0.3984019,0.1091586,4.84];
ybs[431]=['',0.3945065,-0.8141227,6.31];
ybs[432]=['',0.3987763,-0.455485,5.93];
ybs[433]=['η Psc',0.4042044,0.2697581,3.62];
ybs[434]=['',0.4073593,0.609295,6.39];
ybs[435]=['',0.4138284,1.0199222,5.7];
ybs[436]=['δ Phe',0.4020319,-0.8545574,3.95];
ybs[437]=['',0.4045427,-0.5266217,5.82];
ybs[438]=['χ Cas',0.4160736,1.035706,4.71];
ybs[439]=['',0.4038817,-0.7935209,6.17];
ybs[440]=['',0.4107209,-0.1554193,6.59];
ybs[441]=['',0.4096865,-0.6415025,5.51];
ybs[442]=['',0.4168078,0.6518246,5.88];
ybs[443]=['',0.4079017,-0.8659942,6.28];
ybs[444]=['',0.4136148,-0.1206993,5.76];
ybs[445]=['',0.4327885,1.2986941,6.58];
ybs[446]=['',0.418814,0.3241082,5.89];
ybs[447]=['49 Cet',0.4174738,-0.2716882,5.63];
ybs[448]=['',0.4238875,0.7188248,6.38];
ybs[449]=['',0.4181233,-0.5547134,6.12];
ybs[450]=['',0.4266456,0.8522772,5.92];
ybs[451]=['101 Psc',0.422922,0.2577967,6.22];
ybs[452]=['40 Cas',0.4375494,1.2766842,5.28];
ybs[453]=['',0.423578,0.3061805,5.8];
ybs[454]=['υ And',0.4279244,0.7245664,4.09];
ybs[455]=['50 Cet',0.42338,-0.266879,5.42];
ybs[456]=['',0.4190696,-1.0128149,6.01];
ybs[457]=['',0.4343658,1.0137961,5.56];
ybs[458]=['τ Scl',0.4238073,-0.5200782,5.69];
ybs[459]=['π Psc',0.4286661,0.2138143,5.57];
ybs[460]=['51 And',0.4333764,0.850623,3.57];
ybs[461]=['',0.4356092,0.794276,6.36];
ybs[462]=['',0.4306492,-0.1622284,6.24];
ybs[463]=['',0.4093824,-1.3682491,6.11];
ybs[464]=['',0.4255585,-1.0151137,6.18];
ybs[465]=['χ And',0.4391701,0.7765772,4.98];
ybs[466]=['',0.4433151,0.9420696,6.39];
ybs[467]=['',0.4337602,-0.6356422,5.94];
ybs[468]=['α Eri',0.4298254,-0.9970677,0.46];
ybs[469]=['',0.43584,-0.3694278,5.58];
ybs[470]=['',0.4356413,-0.4348194,6.7];
ybs[471]=['105 Psc',0.4400056,0.2882282,5.97];
ybs[472]=['',0.4448785,0.7575771,5.61];
ybs[473]=['τ And',0.4444373,0.7100899,4.94];
ybs[474]=['43 Cas',0.4536184,1.1894564,5.59];
ybs[475]=['',0.4347395,-0.9307882,6.84];
ybs[476]=['42 Cas',0.4565387,1.2344737,5.18];
ybs[477]=['',0.4517256,1.0672024,6.71];
ybs[478]=['',0.4526534,1.0251294,6.37];
ybs[479]=['',0.4497684,0.7456317,4.95];
ybs[480]=['',0.4472762,0.4512355,6.17];
ybs[481]=['',0.448877,0.5263075,5.99];
ybs[482]=['',0.4389057,-0.9789485,5.87];
ybs[483]=['',0.4389349,-0.9788903,5.76];
ybs[484]=['',0.4557947,1.0738891,6.34];
ybs[485]=['ν Psc',0.4474744,0.0976605,4.44];
ybs[486]=['',0.4507623,0.6170339,5.64];
ybs[487]=['44 Cas',0.4572895,1.0586937,5.78];
ybs[488]=['',0.4485747,-0.1957692,5.75];
ybs[489]=['107 Psc',0.4523706,0.3556353,5.24];
ybs[490]=['',0.4467848,-0.6636615,6.17];
ybs[491]=['',0.4563576,0.7928999,6.34];
ybs[492]=['φ Per',0.4582342,0.8865595,4.07];
ybs[493]=['π Scl',0.4499212,-0.5623285,5.25];
ybs[494]=['',0.4494097,-0.6409649,5.72];
ybs[495]=['',0.461351,1.006073,6.21];
ybs[496]=['',0.4529727,-0.0625269,4.99];
ybs[497]=['',0.4474202,-0.8714583,6.64];
ybs[498]=['',0.4633829,0.9982656,6.25];
ybs[499]=['',0.4584465,0.5637265,6.34];
ybs[500]=['',0.461485,0.8071633,6.35];
ybs[501]=['',0.4473869,-1.0590911,5.71];
ybs[502]=['',0.450764,-0.9360676,5.52];
ybs[503]=['',0.4581398,-0.0812988,6.19];
ybs[504]=['109 Psc',0.4629919,0.352387,6.27];
ybs[505]=['τ Cet',0.4586299,-0.2762865,3.5];
ybs[506]=['ο Psc',0.4648244,0.1617032,4.26];
ybs[507]=['',0.4768867,1.1162907,5.63];
ybs[508]=['',0.4251394,-1.4462847,5.87];
ybs[509]=['',0.4671762,-0.0981981,5.34];
ybs[510]=['ε Scl',0.4653434,-0.4353796,5.31];
ybs[511]=['',0.4701837,0.30578,6.55];
ybs[512]=['τ1 Hyi',0.4424115,-1.379511,6.33];
ybs[513]=['',0.4669099,-0.4754655,6.39];
ybs[514]=['',0.4761701,0.8087204,6.32];
ybs[515]=['',0.4665955,-0.8850458,5.49];
ybs[516]=['',0.4665146,-0.9322667,5.04];
ybs[517]=['',0.4796252,0.6642569,5.94];
ybs[518]=['4 Ari',0.4771317,0.2977885,5.84];
ybs[519]=['',0.4796889,0.5724089,5.79];
ybs[520]=['',0.4720621,-0.726987,6.18];
ybs[521]=['',0.4208147,-1.4776059,5.69];
ybs[522]=['',0.482621,0.8378127,5.82];
ybs[523]=['',0.4780003,0.0661824,5.91];
ybs[524]=['',0.4744407,-0.6466992,6.32];
ybs[525]=['',0.4901771,0.9082538,5.9];
ybs[526]=['1 Ari',0.4857904,0.3906268,5.86];
ybs[527]=['χ Cet',0.4827816,-0.1846601,4.67];
ybs[528]=['',0.4812492,-0.5404684,6.34];
ybs[529]=['1 Per',0.4948687,0.9643469,5.52];
ybs[530]=['',0.4887467,0.1945895,5.94];
ybs[531]=['',0.4831942,-0.6684224,6.37];
ybs[532]=['2 Per',0.4953967,0.8883421,5.79];
ybs[533]=['',0.4851733,-0.8327036,6.14];
ybs[534]=['',0.4984444,0.9002413,6.26];
ybs[535]=['ζ Cet',0.4909662,-0.1785352,3.73];
ybs[536]=['',0.5028554,0.9722027,6.45];
ybs[537]=['',0.4875683,-0.8744145,5.94];
ybs[538]=['ε Cas',0.5059884,1.1130819,3.38];
ybs[539]=['55 And',0.4999812,0.712704,5.4];
ybs[540]=['α Tri',0.4987856,0.5180862,3.41];
ybs[541]=['γ1 Ari',0.5005266,0.3386112,4.83];
ybs[542]=['γ2 Ari',0.5005265,0.3385724,4.75];
ybs[543]=['',0.4969882,-0.2936311,5.8];
ybs[544]=['ω Cas',0.5135738,1.2006076,4.99];
ybs[545]=['ξ Psc',0.5003418,0.0574678,4.62];
ybs[546]=['τ2 Hyi',0.4696318,-1.3974842,6.06];
ybs[547]=['',0.5070091,0.7122121,6.24];
ybs[548]=['',0.5071743,0.6498406,6.26];
ybs[549]=['β Ari',0.5054031,0.3649996,2.64];
ybs[550]=['',0.4987577,-0.6717684,6.1];
ybs[551]=['ψ Phe',0.4996432,-0.8062955,4.41];
ybs[552]=['',0.511322,0.6524447,5.89];
ybs[553]=['56 And',0.5124064,0.6519878,5.67];
ybs[554]=['φ Phe',0.5029092,-0.7398792,5.11];
ybs[555]=['7 Ari',0.5107438,0.4133254,5.74];
ybs[556]=['',0.5105364,0.0341089,6.01];
ybs[557]=['',0.524036,1.0786465,6.02];
ybs[558]=['',0.5203407,0.7295209,6.78];
ybs[559]=['ι Ari',0.5171828,0.3127925,5.1];
ybs[560]=['',0.5190472,0.4870958,5.82];
ybs[561]=['56 Cet',0.5134541,-0.3913474,4.85];
ybs[562]=['χ Eri',0.5094923,-0.8989194,3.7];
ybs[563]=['',0.5290214,1.1296632,5.26];
ybs[564]=['3 Per',0.5233293,0.8605873,5.69];
ybs[565]=['λ Ari',0.5198198,0.4136457,4.79];
ybs[566]=['η2 Hyi',0.5038653,-1.1788361,4.69];
ybs[567]=['',0.50815,-1.060405,6.06];
ybs[568]=['',0.5463442,1.3616878,6.04];
ybs[569]=['',0.5140269,-0.9016682,6.1];
ybs[570]=['',0.5149315,-0.8252042,4.83];
ybs[571]=['48 Cas',0.5400361,1.2393556,4.54];
ybs[572]=['',0.5209487,-0.5753081,6.35];
ybs[573]=['',0.5270435,0.3693456,5.87];
ybs[574]=['',0.5261557,0.2163926,6.09];
ybs[575]=['',0.5459928,1.2907253,6.23];
ybs[576]=['50 Cas',0.5468111,1.2657806,3.98];
ybs[577]=['47 Cas',0.555626,1.3505948,5.38];
ybs[578]=['112 Psc',0.5291337,0.0558627,5.88];
ybs[579]=['57 Cet',0.5269989,-0.3616471,5.41];
ybs[580]=['',0.5169675,-1.1400594,6.37];
ybs[581]=['υ Cet',0.5280256,-0.3660697,4];
ybs[582]=['52 Cas',0.5432633,1.134535,6];
ybs[583]=['',0.5302063,-0.1469604,5.51];
ybs[584]=['',0.5259247,-0.7317628,5.57];
ybs[585]=['53 Cas',0.5437655,1.1256091,5.58];
ybs[586]=['4 Per',0.5399497,0.9527816,5.04];
ybs[587]=['α Hyi',0.5211278,-1.0727811,2.86];
ybs[588]=['49 Cas',0.5569433,1.3302359,5.22];
ybs[589]=['σ Hyi',0.5053801,-1.3656085,6.16];
ybs[590]=['π For',0.5332295,-0.5218266,5.35];
ybs[591]=['α Psc',0.5373895,0.0500313,3.82];
ybs[592]=['α Psc',0.5373895,0.0500313,4.33];
ybs[593]=['',0.5769333,1.4206379,6.05];
ybs[594]=['',0.5511251,1.1380512,6.52];
ybs[595]=['ε Tri',0.542057,0.582706,5.5];
ybs[596]=['',0.5246853,-1.1512666,6.1];
ybs[597]=['',0.5399506,0.2370068,5.94];
ybs[598]=['χ Phe',0.5348047,-0.7786003,5.14];
ybs[599]=['γ1 And',0.5464096,0.7405811,2.26];
ybs[600]=['γ2 And',0.5464607,0.7406004,4.84];
ybs[601]=['10 Ari',0.5448803,0.4544503,5.63];
ybs[602]=['',0.5385701,-0.5159563,6.42];
ybs[603]=['60 Cet',0.5423511,0.0040317,5.43];
ybs[604]=['',0.5411026,-0.2653441,5.86];
ybs[605]=['',0.5449515,0.3203701,6.21];
ybs[606]=['61 Cet',0.5449965,-0.00415,5.93];
ybs[607]=['',0.5443657,-0.0698318,5.62];
ybs[608]=['ν For',0.5473915,-0.509542,4.69];
ybs[609]=['κ Ari',0.5575125,0.3970637,5.03];
ybs[610]=['',0.5556414,0.1457235,6.31];
ybs[611]=['11 Ari',0.5586967,0.4504064,6.15];
ybs[612]=['',0.5567192,0.0023872,6.28];
ybs[613]=['α Ari',0.5601873,0.4112706,2];
ybs[614]=['',0.5680814,1.0214491,5.67];
ybs[615]=['',0.5668545,0.7777294,6.42];
ybs[616]=['58 And',0.5663131,0.6625332,4.82];
ybs[617]=['',0.5741552,0.9414963,6.31];
ybs[618]=['β Tri',0.5708355,0.6124032,3];
ybs[619]=['14 Ari',0.570064,0.4544954,4.98];
ybs[620]=['',0.5697057,0.3023853,6.43];
ybs[621]=['',0.5662798,-0.3085443,6.1];
ybs[622]=['',0.5907234,1.2937687,6.29];
ybs[623]=['5 Per',0.5803148,1.0078605,6.36];
ybs[624]=['59 And',0.5768012,0.6831212,5.63];
ybs[625]=['59 And',0.5768669,0.6831842,6.1];
ybs[626]=['',0.5697005,-0.4231534,6.48];
ybs[627]=['15 Ari',0.5751764,0.3421,5.7];
ybs[628]=['',0.5672922,-0.7577451,5.85];
ybs[629]=['16 Ari',0.57783,0.4544381,6.02];
ybs[630]=['5 Tri',0.578921,0.5519911,6.23];
ybs[631]=['64 Cet',0.5781103,0.1513223,5.63];
ybs[632]=['',0.5713175,-0.7629663,6.32];
ybs[633]=['',0.5725191,-0.885296,6.12];
ybs[634]=['',0.5778174,-0.1736918,6.01];
ybs[635]=['63 Cet',0.5789716,-0.0301058,5.93];
ybs[636]=['55 Cas',0.5942944,1.1628059,6.07];
ybs[637]=['',0.5900435,1.023824,6.44];
ybs[638]=['6 Tri',0.5830592,0.5306352,4.94];
ybs[639]=['60 And',0.587221,0.7737308,4.83];
ybs[640]=['',0.5840098,0.4235532,5.96];
ybs[641]=['',0.5892027,0.8930074,5.31];
ybs[642]=['η Ari',0.5847116,0.371944,5.27];
ybs[643]=['',0.5909484,0.8304934,6.06];
ybs[644]=['19 Ari',0.5856854,0.2684254,5.71];
ybs[645]=['ξ1 Cet',0.5853101,0.1561477,4.37];
ybs[646]=['66 Cet',0.5841741,-0.040031,5.54];
ybs[647]=['',0.5847554,-0.3647795,5.86];
ybs[648]=['μ For',0.584049,-0.5344878,5.28];
ybs[649]=['',0.5993639,0.8361947,6.33];
ybs[650]=['',0.6038048,0.9975261,6.48];
ybs[651]=['7 Tri',0.5987171,0.5839514,5.28];
ybs[652]=['20 Ari',0.5977688,0.4517292,5.79];
ybs[653]=['21 Ari',0.5975174,0.4388141,5.58];
ybs[654]=['',0.5957211,-0.1634732,6.55];
ybs[655]=['',0.5908076,-0.7167567,5.91];
ybs[656]=['δ Tri',0.6036209,0.5990476,4.87];
ybs[657]=['8 Per',0.6088385,1.0122585,5.75];
ybs[658]=['7 Per',0.6091459,1.0055725,5.98];
ybs[659]=['',0.606171,0.7750225,6.7];
ybs[660]=['γ Tri',0.6047526,0.5924674,4.01];
ybs[661]=['',0.6038613,0.4165488,6.55];
ybs[662]=['67 Cet',0.6023788,-0.1103647,5.51];
ybs[663]=['π1 Hyi',0.5877293,-1.1823208,5.55];
ybs[664]=['',0.6191273,1.1246022,6.6];
ybs[665]=['θ Ari',0.6079236,0.3490578,5.62];
ybs[666]=['62 And',0.6138359,0.8286484,5.3];
ybs[667]=['',0.6133654,0.81281,6.21];
ybs[668]=['',0.607086,0.0323977,5.58];
ybs[669]=['',0.6143402,0.8561416,6.37];
ybs[670]=['φ Eri',0.5989702,-0.8973305,3.56];
ybs[671]=['10 Tri',0.6117458,0.5016194,5.03];
ybs[672]=['',0.6116772,0.4060675,6.46];
ybs[673]=['',0.6150256,0.6969617,6.63];
ybs[674]=['π2 Hyi',0.5930869,-1.180664,5.69];
ybs[675]=['',0.6199965,0.8274338,6.11];
ybs[676]=['',0.6166948,0.5285934,6.47];
ybs[677]=['ο Cet',0.6127541,-0.0502553,3.04];
ybs[678]=['63 And',0.6213591,0.8770093,5.59];
ybs[679]=['',0.6106211,-0.4511213,6.34];
ybs[680]=['',0.6141869,-0.0741342,6.5];
ybs[681]=['9 Per',0.6277598,0.9763836,5.17];
ybs[682]=['',0.6120359,-0.7286791,6.37];
ybs[683]=['',0.629139,0.724196,5.82];
ybs[684]=['',0.6134865,-0.9747096,5.81];
ybs[685]=['69 Cet',0.6241621,0.0086067,5.28];
ybs[686]=['',0.6343177,0.9679786,6.28];
ybs[687]=['70 Cet',0.6252839,-0.0137494,5.42];
ybs[688]=['',0.6242712,-0.1864098,5.46];
ybs[689]=['',0.624377,-0.3065663,5.87];
ybs[690]=['64 And',0.6364136,0.874465,5.19];
ybs[691]=['κ For',0.6262331,-0.4139791,5.2];
ybs[692]=['10 Per',0.6405426,0.98971,6.25];
ybs[693]=['',0.6281983,-0.3186525,6.22];
ybs[694]=['',0.6241165,-0.7522847,6.31];
ybs[695]=['65 And',0.641711,0.8792048,4.71];
ybs[696]=['',0.6282953,-0.6541392,6.53];
ybs[697]=['',0.6268583,-0.8900336,5.92];
ybs[698]=['ξ Ari',0.636923,0.1868717,5.47];
ybs[699]=['',0.6339918,-0.4494383,6.44];
ybs[700]=['71 Cet',0.6373145,-0.0468383,6.33];
ybs[701]=['δ Hyi',0.6202057,-1.1966323,4.09];
ybs[702]=['',0.6345186,-0.7111174,6.18];
ybs[703]=['ι Cas',0.6582963,1.1780531,4.52];
ybs[704]=['ρ Cet',0.6413647,-0.2128339,4.89];
ybs[705]=['66 And',0.6515204,0.8842734,6.12];
ybs[706]=['',0.6415343,-0.2660764,5.83];
ybs[707]=['',0.647367,0.4731414,6.18];
ybs[708]=['11 Tri',0.6490125,0.5567068,5.54];
ybs[709]=['',0.643958,-0.3481391,5.88];
ybs[710]=['λ Hor',0.6348812,-1.0509583,5.35];
ybs[711]=['κ Hyi',0.6240582,-1.2836666,5.01];
ybs[712]=['',0.6586194,0.9709492,6.51];
ybs[713]=['12 Tri',0.652028,0.5194935,5.29];
ybs[714]=['ξ2 Cet',0.6514596,0.1493193,4.28];
ybs[715]=['',0.6506241,0.0358884,6.45];
ybs[716]=['13 Tri',0.6548382,0.5240715,5.89];
ybs[717]=['κ Eri',0.6447766,-0.8309181,4.25];
ybs[718]=['',0.6365575,-1.1588706,6.41];
ybs[719]=['',0.6564852,0.4112677,6.19];
ybs[720]=['φ For',0.6498607,-0.5884494,5.14];
ybs[721]=['',0.6577259,0.1686118,6.07];
ybs[722]=['',0.6613723,0.592165,6.25];
ybs[723]=['',0.6523926,-0.5411784,6.11];
ybs[724]=['',0.6622674,0.4420847,5.92];
ybs[725]=['26 Ari',0.6625628,0.3481905,6.15];
ybs[726]=['',0.6584404,-0.3942341,6.77];
ybs[727]=['27 Ari',0.6636732,0.3106402,6.23];
ybs[728]=['',0.6626004,0.0061055,6];
ybs[729]=['',0.6597109,-0.437932,6.51];
ybs[730]=['',0.6482882,-1.1205753,6.37];
ybs[731]=['',0.661162,-0.3918426,6.1];
ybs[732]=['14 Tri',0.6694292,0.6325299,5.15];
ybs[733]=['',0.6659108,0.0412164,5.25];
ybs[734]=['',0.6727456,0.604518,5.83];
ybs[735]=['75 Cet',0.6686945,-0.0164219,5.35];
ybs[736]=['σ Cet',0.6680608,-0.2644277,4.75];
ybs[737]=['29 Ari',0.672317,0.2640432,6.04];
ybs[738]=['',0.6681583,-0.6341374,6.3];
ybs[739]=['',0.6986931,1.2725254,5.16];
ybs[740]=['λ1 For',0.6720105,-0.6031192,5.9];
ybs[741]=['',0.6748474,-0.3474658,6.21];
ybs[742]=['',0.6842448,0.6938979,6.36];
ybs[743]=['',0.6954384,1.1490855,5.78];
ybs[744]=['',0.6849501,0.6528428,5.71];
ybs[745]=['ω For',0.6753854,-0.491117,4.9];
ybs[746]=['15 Tri',0.6854402,0.6070321,5.35];
ybs[747]=['',0.6815945,0.1320259,6.18];
ybs[748]=['77 Cet',0.6796854,-0.1355453,5.75];
ybs[749]=['',0.6860032,0.1218197,5.82];
ybs[750]=['ν Cet',0.6850712,0.0992431,4.86];
ybs[751]=['',0.6747561,-0.8901183,6.24];
ybs[752]=['',0.6907117,0.6776287,5.9];
ybs[753]=['',0.6894227,0.5532706,6.1];
ybs[754]=['',0.6909369,0.5996315,5.3];
ybs[755]=['80 Cet',0.6853062,-0.1350679,5.53];
ybs[756]=['',0.6924754,0.6979255,6.54];
ybs[757]=['',0.6911751,0.5756861,6.25];
ybs[758]=['',0.6723399,-1.0907123,6.77];
ybs[759]=['31 Ari',0.6885359,0.2188664,5.68];
ybs[760]=['30 Ari',0.6902923,0.431809,7.09];
ybs[761]=['30 Ari',0.690496,0.4317941,6.5];
ybs[762]=['',0.6882197,0.136526,5.81];
ybs[763]=['ι1 For',0.6854013,-0.5227593,5.75];
ybs[764]=['',0.6965426,0.6600609,6.18];
ybs[765]=['',0.6972843,0.6663916,6.3];
ybs[766]=['',0.6944525,0.1359165,6.39];
ybs[767]=['81 Cet',0.6928117,-0.0576627,5.65];
ybs[768]=['λ2 For',0.6888364,-0.6018905,5.79];
ybs[769]=['ν Ari',0.6983176,0.3849022,5.43];
ybs[770]=['',0.746354,1.4230854,5.78];
ybs[771]=['',0.6969781,0.0616979,6.21];
ybs[772]=['μ Hyi',0.6599796,-1.3790711,5.28];
ybs[773]=['ι2 For',0.6947981,-0.5253801,5.83];
ybs[774]=['η Hor',0.6898948,-0.9154359,5.31];
ybs[775]=['δ Cet',0.7006942,0.0073354,4.07];
ybs[776]=['',0.6949759,-0.6614532,6.49];
ybs[777]=['ε Cet',0.7007532,-0.2056095,4.84];
ybs[778]=['33 Ari',0.7066207,0.4738932,5.3];
ybs[779]=['',0.704227,0.1082691,6.25];
ybs[780]=['',0.703611,-0.1633908,5.78];
ybs[781]=['11 Per',0.7182,0.9633559,5.77];
ybs[782]=['',0.7023257,-0.5330648,6.52];
ybs[783]=['',0.7178628,0.935785,5.84];
ybs[784]=['12 Per',0.7139036,0.7030987,4.91];
ybs[785]=['',0.7008203,-0.7470017,4.75];
ybs[786]=['84 Cet',0.7082986,-0.0105547,5.71];
ybs[787]=['',0.7274074,1.1853307,5.95];
ybs[788]=['',0.7177049,0.843971,6.48];
ybs[789]=['μ Ari',0.7137641,0.3508523,5.69];
ybs[790]=['ι Eri',0.7047232,-0.6940168,4.11];
ybs[791]=['',0.7107404,-0.0544967,6.05];
ybs[792]=['',0.7094246,-0.2523476,5.98];
ybs[793]=['',0.714028,0.1890598,6.3];
ybs[794]=['',0.698103,-1.1203299,6.55];
ybs[795]=['θ Per',0.7228652,0.8607676,4.12];
ybs[796]=['14 Per',0.7221159,0.7746994,5.43];
ybs[797]=['35 Ari',0.718722,0.4851585,4.66];
ybs[798]=['ζ Hor',0.703927,-0.9504826,5.21];
ybs[799]=['',0.7204138,0.4490423,6.35];
ybs[800]=['γ Cet',0.7174173,0.0580534,3.47];
ybs[801]=['',0.7110705,-0.6683398,6.01];
ybs[802]=['ε Hyi',0.6977941,-1.1898812,4.11];
ybs[803]=['',0.7108521,-0.8104193,6.1];
ybs[804]=['36 Ari',0.7222244,0.3116094,6.46];
ybs[805]=['ο Ari',0.7231641,0.2688088,5.77];
ybs[806]=['ι Hor',0.7124339,-0.885049,5.41];
ybs[807]=['π Cet',0.7205931,-0.2403055,4.25];
ybs[808]=['38 Ari',0.7248878,0.2187882,5.18];
ybs[809]=['μ Cet',0.7247493,0.178093,4.27];
ybs[810]=['',0.7163376,-0.7057602,6.36];
ybs[811]=['',0.7456919,1.2168863,6.18];
ybs[812]=['',0.7263914,0.0837993,6.03];
ybs[813]=['',0.7210079,-0.5660964,6.22];
ybs[814]=['τ1 Eri',0.7247443,-0.3225843,4.47];
ybs[815]=['',0.7343743,0.6295872,6.25];
ybs[816]=['',0.734737,0.622106,6.3];
ybs[817]=['',0.7193787,-0.9159558,6.15];
ybs[818]=['',0.7245131,-0.8062977,6.85];
ybs[819]=['',0.7147923,-1.1628078,6.26];
ybs[820]=['39 Ari',0.7382327,0.5120094,4.51];
ybs[821]=['',0.7466368,0.997845,6.25];
ybs[822]=['',0.7318555,-0.3761273,6.49];
ybs[823]=['',0.7337211,-0.3908925,6.47];
ybs[824]=['40 Ari',0.7406464,0.3206547,5.82];
ybs[825]=['',0.7589201,1.203855,5.8];
ybs[826]=['',0.7418475,0.4411585,5.86];
ybs[827]=['',0.745259,0.653003,6.45];
ybs[828]=['',0.7372471,-0.2159279,6.9];
ybs[829]=['γ Hor',0.7239647,-1.1102858,5.74];
ybs[830]=['η Per',0.7517152,0.9770925,3.76];
ybs[831]=['η1 For',0.7349365,-0.6189266,6.51];
ybs[832]=['π Ari',0.7439273,0.3063481,5.22];
ybs[833]=['ζ Hyi',0.7237782,-1.1785671,4.84];
ybs[834]=['41 Ari',0.747231,0.4773229,3.63];
ybs[835]=['',0.7565635,1.0193082,6.45];
ybs[836]=['16 Per',0.7502494,0.6703183,4.23];
ybs[837]=['β For',0.74171,-0.5640454,4.46];
ybs[838]=['',0.7554936,0.8190715,5.88];
ybs[839]=['17 Per',0.7541785,0.6134344,4.53];
ybs[840]=['γ1 For',0.7452738,-0.4271195,6.14];
ybs[841]=['γ2 For',0.7454027,-0.4861361,5.39];
ybs[842]=['',0.7609892,0.9265034,6.36];
ybs[843]=['σ Ari',0.7534742,0.2647569,5.49];
ybs[844]=['η2 For',0.7466281,-0.624053,5.92];
ybs[845]=['',0.7628372,0.8492115,6.26];
ybs[846]=['τ2 Eri',0.750559,-0.3650611,4.75];
ybs[847]=['η3 For',0.748494,-0.6211322,5.47];
ybs[848]=['ν Hor',0.7395866,-1.0946383,5.26];
ybs[849]=['',0.74887,-0.6954064,6.36];
ybs[850]=['τ Per',0.7670574,0.9223882,3.95];
ybs[851]=['20 Per',0.763902,0.6706285,5.33];
ybs[852]=['',0.7609401,0.2892051,6.31];
ybs[853]=['',0.7573107,-0.2213475,6.04];
ybs[854]=['',0.7541081,-0.536288,6.4];
ybs[855]=['',0.7587348,-0.1632591,6.32];
ybs[856]=['',0.7752172,1.0752435,5.59];
ybs[857]=['',0.7775999,1.1243083,6.24];
ybs[858]=['',0.7616391,-0.3890268,5.95];
ybs[859]=['ψ For',0.7610251,-0.6693358,5.92];
ybs[860]=['',0.7782574,0.8961629,6.22];
ybs[861]=['',0.7767585,0.8246597,6.02];
ybs[862]=['',0.7539086,-1.0964574,6.03];
ybs[863]=['ρ2 Ari',0.7723996,0.321448,5.91];
ybs[864]=['',0.7618036,-0.8692357,4];
ybs[865]=['ρ3 Ari',0.7751266,0.3160579,5.63];
ybs[866]=['',0.7739722,0.1477854,5.97];
ybs[867]=['',0.7627303,-0.8863607,6.21];
ybs[868]=['ν Hyi',0.7433593,-1.3086276,4.75];
ybs[869]=['21 Per',0.7792883,0.5588468,5.11];
ybs[870]=['η Eri',0.7743948,-0.1538037,3.89];
ybs[871]=['',0.7753792,-0.0632951,5.17];
ybs[872]=['',0.782822,0.6754443,6.04];
ybs[873]=['',0.777567,0.0800515,6.11];
ybs[874]=['47 Ari',0.7824131,0.3622209,5.8];
ybs[875]=['π Per',0.7860267,0.6937268,4.7];
ybs[876]=['',0.7625382,-1.1231007,6.56];
ybs[877]=['',0.825152,1.3875436,5.49];
ybs[878]=['24 Per',0.7871538,0.6155391,4.93];
ybs[879]=['4 Eri',0.7781937,-0.4149786,5.45];
ybs[880]=['',0.7772312,-0.5195808,6.29];
ybs[881]=['',0.7910745,0.8256326,5.47];
ybs[882]=['',0.7900342,0.7176368,5.89];
ybs[883]=['ε Ari',0.7873444,0.3739364,4.63];
ybs[884]=['ε Ari',0.7873444,0.3739364,4.63];
ybs[885]=['6 Eri',0.7812539,-0.410518,5.84];
ybs[886]=['',0.7959338,0.9151756,5.28];
ybs[887]=['',0.7960212,0.9151851,6.74];
ybs[888]=['',0.7844634,-0.0470817,5.23];
ybs[889]=['',0.7783443,-0.6650703,6.41];
ybs[890]=['',0.7922302,0.6669946,6.11];
ybs[891]=['',0.7846739,-0.1691486,6.14];
ybs[892]=['λ Cet',0.7891956,0.1569405,4.7];
ybs[893]=['θ1 Eri',0.7813776,-0.7019643,3.24];
ybs[894]=['θ2 Eri',0.7814211,-0.7019596,4.35];
ybs[895]=['5 Eri',0.7887765,-0.0415467,5.56];
ybs[896]=['',0.7855101,-0.5030414,6.14];
ybs[897]=['ζ For',0.7877718,-0.4396407,5.71];
ybs[898]=['',0.7936997,0.1911908,5.95];
ybs[899]=['',0.787688,-0.5658814,6.31];
ybs[900]=['7 Eri',0.7938414,-0.048773,6.11];
ybs[901]=['49 Ari',0.7992485,0.4633137,5.9];
ybs[902]=['',0.8517757,1.4233174,5.95];
ybs[903]=['ρ1 Eri',0.7950977,-0.1322745,5.75];
ybs[904]=['',0.7985134,0.0945941,6.25];
ybs[905]=['β Hor',0.7819363,-1.1167731,4.99];
ybs[906]=['93 Cet',0.8006842,0.0774284,5.61];
ybs[907]=['α Cet',0.8002626,0.0728379,2.53];
ybs[908]=['',0.7983882,-0.1723977,5.83];
ybs[909]=['',0.7994402,-0.1118945,6.19];
ybs[910]=['ε For',0.796517,-0.4888284,5.89];
ybs[911]=['γ Per',0.8131836,0.9353035,2.93];
ybs[912]=['',0.8063097,0.4948499,6.36];
ybs[913]=['ρ2 Eri',0.8018081,-0.1326772,5.32];
ybs[914]=['',0.8166936,0.9911392,4.76];
ybs[915]=['τ3 Eri',0.799989,-0.4108657,4.09];
ybs[916]=['',0.8171849,0.9800168,6.11];
ybs[917]=['ρ Per',0.814015,0.6793296,3.39];
ybs[918]=['',0.8252804,1.1194423,5.89];
ybs[919]=['',0.8148357,0.7097359,6.05];
ybs[920]=['',0.8110496,0.2781841,6.49];
ybs[921]=['ρ3 Eri',0.8086521,-0.1312137,5.26];
ybs[922]=['',0.8104805,0.0339695,6.05];
ybs[923]=['52 Ari',0.814681,0.4422253,6.8];
ybs[924]=['52 Ari',0.814681,0.4422253,7];
ybs[925]=['',0.8013951,-0.8184126,5.82];
ybs[926]=['',0.8273768,0.9127137,6.31];
ybs[927]=['',0.8184674,0.2315919,5.62];
ybs[928]=['',0.8478948,1.2998037,4.87];
ybs[929]=['',0.8258587,0.8271074,6.41];
ybs[930]=['μ Hor',0.8034025,-1.0411688,5.11];
ybs[931]=['',0.8186611,-0.1048355,5.27];
ybs[932]=['β Per',0.82718,0.7162282,2.12];
ybs[933]=['ι Per',0.831579,0.8673285,4.05];
ybs[934]=['53 Ari',0.8231127,0.3134894,6.11];
ybs[935]=['θ Hyi',0.795499,-1.253473,5.53];
ybs[936]=['54 Ari',0.8271778,0.3294529,6.27];
ybs[937]=['κ Per',0.8331794,0.7843161,3.8];
ybs[938]=['',0.828157,0.1492605,6.28];
ybs[939]=['',0.8236449,-0.4843218,6.19];
ybs[940]=['55 Ari',0.8330039,0.508898,5.72];
ybs[941]=['',0.8352961,0.4869566,6.42];
ybs[942]=['',0.836589,0.4708345,6.02];
ybs[943]=['ω Per',0.8407548,0.692752,4.63];
ybs[944]=['',0.8369949,0.2086172,5.98];
ybs[945]=['',0.8462041,0.8343627,6.33];
ybs[946]=['',0.8446876,0.7409946,6.15];
ybs[947]=['δ Ari',0.8415101,0.3456915,4.35];
ybs[948]=['',0.8401563,0.2291249,6.12];
ybs[949]=['',0.835735,-0.4129079,6.38];
ybs[950]=['56 Ari',0.8444027,0.4771153,5.79];
ybs[951]=['',0.8394677,-0.0651272,6.05];
ybs[952]=['',0.8504268,0.8422295,5.9];
ybs[953]=['',0.8389801,-0.2782945,6.26];
ybs[954]=['',0.8446691,0.1176445,5.56];
ybs[955]=['',0.8203016,-1.2074866,6.15];
ybs[956]=['',0.8340736,-0.8491658,6.12];
ybs[957]=['',0.8863369,1.3580575,5.45];
ybs[958]=['94 Cet',0.8459089,-0.0194872,5.06];
ybs[959]=['α For',0.8420183,-0.5045233,3.87];
ybs[960]=['',0.8615973,0.9986619,5.79];
ybs[961]=['',0.9503505,1.4832133,5.61];
ybs[962]=['',0.8568612,0.7432057,6.07];
ybs[963]=['',0.8702049,1.1473123,6.36];
ybs[964]=['',0.8429183,-0.773878,5.93];
ybs[965]=['',0.8628592,0.8903957,5.03];
ybs[966]=['',0.8459081,-0.6259511,6.27];
ybs[967]=['',0.8580792,0.534685,5.52];
ybs[968]=['ζ Ari',0.8558339,0.3686687,4.89];
ybs[969]=['',0.8619641,0.7927991,6.16];
ybs[970]=['',0.8487978,-0.5187971,6.16];
ybs[971]=['',0.8600955,0.5748197,6.31];
ybs[972]=['',0.8612534,0.6067962,6.25];
ybs[973]=['',0.8425561,-0.9990595,5.74];
ybs[974]=['',0.8635704,0.5630719,6.06];
ybs[975]=['',0.8665691,0.707925,6.45];
ybs[976]=['',0.8549012,-0.4541616,6.25];
ybs[977]=['',0.8152797,-1.3771957,5.57];
ybs[978]=['30 Per',0.8693702,0.7697344,5.47];
ybs[979]=['',0.8599222,-0.1019325,6.17];
ybs[980]=['ζ Eri',0.8590483,-0.1525651,4.8];
ybs[981]=['',0.8809842,1.1471836,4.84];
ybs[982]=['',0.8690169,0.686977,5.96];
ybs[983]=['29 Per',0.8734157,0.8778901,5.15];
ybs[984]=['14 Eri',0.8623615,-0.1584123,6.14];
ybs[985]=['31 Per',0.8755776,0.8756663,5.03];
ybs[986]=['',0.8598995,-0.5366752,6.65];
ybs[987]=['',0.873019,0.5986472,4.82];
ybs[988]=['95 Cet',0.8703503,-0.014886,5.38];
ybs[989]=['',0.8680763,-0.5012482,5.91];
ybs[990]=['15 Eri',0.869694,-0.391547,4.88];
ybs[991]=['59 Ari',0.8779883,0.4738188,5.9];
ybs[992]=['κ1 Cet',0.8747828,0.0601658,4.83];
ybs[993]=['',0.8712077,-0.3225797,5.71];
ybs[994]=['',0.8645719,-0.8320656,5.85];
ybs[995]=['',0.8798498,0.508325,4.47];
ybs[996]=['60 Ari',0.8801076,0.4492354,6.12];
ybs[997]=['',0.887527,0.8577719,5.93];
ybs[998]=['32 Per',0.8852868,0.7575692,4.95];
ybs[999]=['τ4 Eri',0.8747263,-0.378402,3.69];
ybs[1000]=['',0.8749246,-0.4196842,5.61];
ybs[1001]=['τ1 Ari',0.8834495,0.3704137,5.28];
ybs[1002]=['ζ1 Ret',0.8646865,-1.0907876,5.54];
ybs[1003]=['κ2 Cet',0.8824351,0.0654816,5.69];
ybs[1004]=['',0.8756643,-0.7503678,4.27];
ybs[1005]=['',0.9013136,1.1285435,5.23];
ybs[1006]=['ζ2 Ret',0.8666321,-1.0895884,5.24];
ybs[1007]=['',0.8934068,0.8602495,5.29];
ybs[1008]=['62 Ari',0.8879117,0.4831647,5.52];
ybs[1009]=['',0.879956,-0.4630347,6.39];
ybs[1010]=['',0.8649497,-1.1667395,6.05];
ybs[1011]=['τ2 Ari',0.8901011,0.3633343,5.09];
ybs[1012]=['',0.8828775,-0.4111837,5.52];
ybs[1013]=['α Per',0.8982721,0.8715475,1.79];
ybs[1014]=['',0.886612,-0.4452673,6.35];
ybs[1015]=['',0.8981679,0.5866171,5.61];
ybs[1016]=['',0.9050608,0.9424067,6.51];
ybs[1017]=['',0.8824843,-0.8325355,6.39];
ybs[1018]=['64 Ari',0.8970227,0.4328261,5.5];
ybs[1019]=['',0.8935459,0.0865192,6.38];
ybs[1020]=['',0.8916361,-0.1347179,6.2];
ybs[1021]=['ι Hyi',0.8528285,-1.3493086,5.52];
ybs[1022]=['',0.9014083,0.7213758,6.51];
ybs[1023]=['65 Ari',0.8974444,0.3643987,6.08];
ybs[1024]=['',0.8960264,0.2217347,6.04];
ybs[1025]=['',0.9053596,0.8586156,6.09];
ybs[1026]=['ο Tau',0.8987352,0.1588886,3.6];
ybs[1027]=['',0.8927644,-0.5695352,6.5];
ybs[1028]=['',0.9276195,1.2555227,6.32];
ybs[1029]=['',0.9169933,1.0529348,6.49];
ybs[1030]=['',0.9145104,0.8575873,4.98];
ybs[1031]=['',0.9199083,1.0474273,4.21];
ybs[1032]=['',0.9088083,0.3286493,6.57];
ybs[1033]=['',0.9181547,0.8712918,5.58];
ybs[1034]=['ξ Tau',0.909041,0.1711569,3.74];
ybs[1035]=['',0.9097489,0.2235546,6.28];
ybs[1036]=['',0.9234914,1.0288917,4.54];
ybs[1037]=['',0.9149996,0.591331,5.61];
ybs[1038]=['χ1 For',0.902163,-0.6256388,6.39];
ybs[1039]=['',0.9247404,1.0373982,6.13];
ybs[1040]=['34 Per',0.9202944,0.8653638,4.67];
ybs[1041]=['',0.9044433,-0.4754857,5.93];
ybs[1042]=['',0.9235376,0.9690847,5.09];
ybs[1043]=['',0.9204249,0.8204892,6.24];
ybs[1044]=['66 Ari',0.9150096,0.3992862,6.03];
ybs[1045]=['',0.903049,-0.7254054,6.32];
ybs[1046]=['',0.9121239,-0.1957071,5.73];
ybs[1047]=['',0.925655,0.8408282,5.82];
ybs[1048]=['σ Per',0.9254659,0.8389377,4.36];
ybs[1049]=['',0.8907438,-1.2138665,6.15];
ybs[1050]=['χ2 For',0.9092645,-0.6214713,5.71];
ybs[1051]=['',0.9494945,1.2813706,6.57];
ybs[1052]=['',0.929545,0.860127,6.29];
ybs[1053]=['χ3 For',0.9120251,-0.6244769,6.5];
ybs[1054]=['',0.9309928,0.8634601,6.39];
ybs[1055]=['',0.9193893,-0.1174992,5.99];
ybs[1056]=['4 Tau',0.9232221,0.1991218,5.14];
ybs[1057]=['',0.9189888,-0.2199446,5.59];
ybs[1058]=['',0.9323204,0.8394207,5.47];
ybs[1059]=['',0.8976005,-1.2088454,5.96];
ybs[1060]=['',0.9278544,0.4824785,5.96];
ybs[1061]=['5 Tau',0.9253003,0.2270485,4.11];
ybs[1062]=['',0.9245891,0.1092737,5.94];
ybs[1063]=['',0.9393069,1.0268824,6.4];
ybs[1064]=['36 Per',0.9334852,0.8050939,5.31];
ybs[1065]=['17 Eri',0.9236572,-0.0873169,4.73];
ybs[1066]=['',0.9398743,1.0112412,6.37];
ybs[1067]=['',0.9343457,0.7841242,6.41];
ybs[1068]=['',0.939435,0.960729,5.98];
ybs[1069]=['',0.9339367,0.6201702,5.9];
ybs[1070]=['',0.9192169,-0.742837,5.78];
ybs[1071]=['',0.9206342,-0.7207754,6.12];
ybs[1072]=['',0.9459162,1.0491438,6.46];
ybs[1073]=['',0.9381557,0.6976172,5.81];
ybs[1074]=['6 Tau',0.9327281,0.164849,5.77];
ybs[1075]=['',0.9689446,1.3231002,6.27];
ybs[1076]=['',0.9219924,-0.8255898,5.99];
ybs[1077]=['',0.9285895,-0.4457969,6.38];
ybs[1078]=['κ Ret',0.9151592,-1.0971916,4.72];
ybs[1079]=['ε Eri',0.9336123,-0.1638324,3.73];
ybs[1080]=['',0.9397055,0.3124778,6.17];
ybs[1081]=['7 Tau',0.9412617,0.4282198,5.92];
ybs[1082]=['ψ Per',0.9513235,0.8423412,4.23];
ybs[1083]=['τ5 Eri',0.9369703,-0.3763225,4.27];
ybs[1084]=['',0.942328,0.1132438,6.49];
ybs[1085]=['',0.9303447,-0.8780219,5.68];
ybs[1086]=['',0.9409869,-0.1710055,6.25];
ybs[1087]=['',0.9210427,-1.1591999,5.83];
ybs[1088]=['',0.9373144,-0.5412134,6.2];
ybs[1089]=['',0.9600704,0.9948688,6.3];
ybs[1090]=['',0.9399492,-0.5550835,6.4];
ybs[1091]=['',0.930567,-1.0636971,6.41];
ybs[1092]=['',0.9575956,0.7444222,6.42];
ybs[1093]=['',0.9467847,-0.194141,5.57];
ybs[1094]=['',0.950744,0.0114767,5.71];
ybs[1095]=['20 Eri',0.9480193,-0.3036335,5.23];
ybs[1096]=['10 Tau',0.9511094,0.0082278,4.28];
ybs[1097]=['',0.9555964,0.270529,6.39];
ybs[1098]=['',0.9610375,0.3662512,6.5];
ybs[1099]=['',0.9366247,-1.1465674,6.75];
ybs[1100]=['',0.9775691,1.104514,5.1];
ybs[1101]=['',0.9506267,-0.7017095,4.58];
ybs[1102]=['',1.1277323,1.5102191,5.86];
ybs[1103]=['',0.9579157,-0.1278033,5.85];
ybs[1104]=['',0.9130149,-1.3662248,5.7];
ybs[1105]=['',0.9627541,0.2898173,6.16];
ybs[1106]=['21 Eri',0.9602886,-0.0969926,5.96];
ybs[1107]=['',0.9795678,1.0478355,5.76];
ybs[1108]=['',0.9710029,0.6570792,5.57];
ybs[1109]=['τ For',0.9585771,-0.4864944,6.01];
ybs[1110]=['12 Tau',0.9641845,0.0545489,5.57];
ybs[1111]=['',0.9630606,-0.0580231,6.23];
ybs[1112]=['',0.9619089,-0.1809652,6.19];
ybs[1113]=['11 Tau',0.9689209,0.4432699,6.11];
ybs[1114]=['',0.9646703,-0.0183632,6.12];
ybs[1115]=['',0.9651025,-0.2645624,6.33];
ybs[1116]=['22 Eri',0.9673688,-0.0897519,5.53];
ybs[1117]=['δ Per',0.9794116,0.8352197,3.01];
ybs[1118]=['40 Per',0.976267,0.5939764,4.97];
ybs[1119]=['',0.9951483,1.1740354,5.8];
ybs[1120]=['',0.9697477,-0.204817,6.49];
ybs[1121]=['13 Tau',0.9754596,0.3450109,5.69];
ybs[1122]=['',0.9846245,0.8480583,6.06];
ybs[1123]=['',0.9701204,-0.3406335,6.59];
ybs[1124]=['',0.9945717,1.1067243,4.8];
ybs[1125]=['',0.9869762,0.8057492,6.11];
ybs[1126]=['ο Per',0.9846805,0.5646984,3.83];
ybs[1127]=['14 Tau',0.9818824,0.3443841,6.14];
ybs[1128]=['',0.9857618,0.637506,5.59];
ybs[1129]=['δ For',0.9734754,-0.5562509,5];
ybs[1130]=['ν Per',0.9890016,0.7442906,3.77];
ybs[1131]=['δ Eri',0.9786159,-0.1692323,3.54];
ybs[1132]=['',0.9849016,0.3664331,6.1];
ybs[1133]=['',1.0101069,1.2380539,5.44];
ybs[1134]=['',0.9799676,-0.1818399,5.6];
ybs[1135]=['16 Tau',0.9864838,0.4250881,5.46];
ybs[1136]=['',0.9926581,0.7984477,5.66];
ybs[1137]=['17 Tau',0.9867903,0.4220138,3.7];
ybs[1138]=['',0.9757937,-0.6500714,4.59];
ybs[1139]=['18 Tau',0.9880683,0.4346797,5.64];
ybs[1140]=['19 Tau',0.9882584,0.4281878,4.3];
ybs[1141]=['24 Eri',0.9843762,-0.0191389,5.25];
ybs[1142]=['',1.0002201,0.9771667,6.1];
ybs[1143]=['γ Cam',1.0151823,1.2460929,4.63];
ybs[1144]=['20 Tau',0.9909541,0.4264475,3.87];
ybs[1145]=['25 Eri',0.9862933,-0.0040209,5.55];
ybs[1146]=['21 Tau',0.9913103,0.4297096,5.76];
ybs[1147]=['22 Tau',0.9919277,0.4292431,6.43];
ybs[1148]=['29 Tau',0.9896796,0.1067436,5.35];
ybs[1149]=['',0.941312,-1.3657681,6.29];
ybs[1150]=['',1.0101602,1.1447649,4.47];
ybs[1151]=['23 Tau',0.9931216,0.419123,4.18];
ybs[1152]=['',0.9811554,-0.7084912,6.45];
ybs[1153]=['',1.0101911,1.105863,5.85];
ybs[1154]=['',0.9918193,0.1198832,5.91];
ybs[1155]=['',1.0031291,0.8866513,6.14];
ybs[1156]=['',1.0081646,0.9980191,6.46];
ybs[1157]=['π Eri',0.9911641,-0.2100658,4.42];
ybs[1158]=['',1.0002749,0.5875641,6.57];
ybs[1159]=['',0.9999446,0.5630428,6.25];
ybs[1160]=['η Tau',0.9981841,0.4218484,2.87];
ybs[1161]=['',1.0202915,1.1967825,6.32];
ybs[1162]=['',0.9839401,-0.8376702,6.49];
ybs[1163]=['',0.9822396,-0.9460962,6.3];
ybs[1164]=['',0.9858174,-0.825427,5.73];
ybs[1165]=['',1.006306,0.7684235,6.02];
ybs[1166]=['σ For',0.9919314,-0.5108995,5.9];
ybs[1167]=['',1.0019279,0.4099057,5.45];
ybs[1168]=['τ6 Eri',0.9938706,-0.4046412,4.23];
ybs[1169]=['30 Tau',1.0011781,0.1956191,5.07];
ybs[1170]=['β Ret',0.9793854,-1.1299289,3.85];
ybs[1171]=['',1.0104656,0.7859519,5.66];
ybs[1172]=['42 Per',1.0075301,0.5786744,5.11];
ybs[1173]=['27 Tau',1.0055017,0.4209339,3.63];
ybs[1174]=['',0.9957386,-0.5207478,6.55];
ybs[1175]=['28 Tau',1.0056139,0.4223881,5.09];
ybs[1176]=['τ7 Eri',0.9973876,-0.4155556,5.24];
ybs[1177]=['',1.0024827,0.0051041,5.91];
ybs[1178]=['',1.0079483,0.4149664,6.17];
ybs[1179]=['ρ For',0.9983368,-0.5253921,5.54];
ybs[1180]=['',1.0087317,0.389357,6.07];
ybs[1181]=['',0.9976123,-0.6290297,6.21];
ybs[1182]=['',1.0015791,-0.3636975,5.81];
ybs[1183]=['',1.0105937,0.4475606,5.26];
ybs[1184]=['',1.0008676,-0.6555013,5.4];
ybs[1185]=['',1.000904,-0.6554724,4.73];
ybs[1186]=['',1.0178555,0.6007831,5.77];
ybs[1187]=['',1.0274671,1.0129415,5.8];
ybs[1188]=['',1.0161159,0.3856303,6.83];
ybs[1189]=['',1.0142962,0.2288008,6.3];
ybs[1190]=['',1.0047069,-0.6306904,4.17];
ybs[1191]=['',1.0422358,1.2545874,6.34];
ybs[1192]=['',1.0184869,0.5450913,6.25];
ybs[1193]=['',1.0263052,0.8502006,5.76];
ybs[1194]=['31 Tau',1.0173178,0.115155,5.67];
ybs[1195]=['',1.0098044,-0.6346264,6.86];
ybs[1196]=['',1.022746,0.3035054,5.97];
ybs[1197]=['30 Eri',1.0199584,-0.0924763,5.48];
ybs[1198]=['ζ Per',1.0275119,0.5575592,2.85];
ybs[1199]=['',1.044322,1.1018752,5.03];
ybs[1200]=['',1.0427858,1.0676111,5];
ybs[1201]=['',1.0218156,-0.3206479,6.22];
ybs[1202]=['',1.0364073,0.8365835,5.37];
ybs[1203]=['γ Hyi',0.9901444,-1.2945682,3.24];
ybs[1204]=['',1.032941,0.5429275,6.1];
ybs[1205]=['43 Per',1.0394046,0.8858644,5.28];
ybs[1206]=['32 Eri',1.027,-0.0504506,6.14];
ybs[1207]=['32 Eri',1.0270071,-0.0504845,4.79];
ybs[1208]=['τ8 Eri',1.0237451,-0.428479,4.65];
ybs[1209]=['',1.0230577,-0.6051005,5.11];
ybs[1210]=['',1.0379024,0.6133477,5.49];
ybs[1211]=['',1.0219836,-0.8173555,5.93];
ybs[1212]=['',1.0309917,-0.2100926,6];
ybs[1213]=['32 Tau',1.0390861,0.3933803,5.63];
ybs[1214]=['',1.0259978,-0.7032808,5.71];
ybs[1215]=['ε Per',1.0441539,0.6993669,2.89];
ybs[1216]=['33 Tau',1.0399631,0.4055524,6.06];
ybs[1217]=['',1.0416557,0.4280011,6.16];
ybs[1218]=['',1.0447599,0.6086811,6.53];
ybs[1219]=['',1.0392305,0.1064813,6.09];
ybs[1220]=['',1.0370008,-0.1691171,6.19];
ybs[1221]=['',1.0468506,0.6789415,6.3];
ybs[1222]=['',1.0259582,-0.9185388,6.46];
ybs[1223]=['ξ Per',1.0487943,0.6257199,4.04];
ybs[1224]=['',1.0520099,0.678588,6.38];
ybs[1225]=['',1.1072392,1.4094057,5.1];
ybs[1226]=['γ Eri',1.0429779,-0.2347136,2.95];
ybs[1227]=['',1.0469104,-0.0944202,5.83];
ybs[1228]=['',1.0509331,0.1813495,6.37];
ybs[1229]=['',1.0588089,0.6466264,6.41];
ybs[1230]=['',1.0494287,-0.2184161,5.6];
ybs[1231]=['',1.0312602,-1.1065747,6.14];
ybs[1232]=['',1.0552552,0.3029186,6.32];
ybs[1233]=['',1.0561523,0.3185765,5.89];
ybs[1234]=['λ Tau',1.0553757,0.2190309,3.47];
ybs[1235]=['τ9 Eri',1.0508726,-0.4181236,4.66];
ybs[1236]=['',1.0831346,1.1996798,5.87];
ybs[1237]=['',1.0745102,1.0334612,5.06];
ybs[1238]=['',1.0600433,0.17552,5.67];
ybs[1239]=['35 Eri',1.0586478,-0.0260197,5.28];
ybs[1240]=['',1.0436126,-0.9955737,6.05];
ybs[1241]=['',1.0538984,-0.5311296,5.93];
ybs[1242]=['δ Ret',1.0432246,-1.0705839,4.56];
ybs[1243]=['',1.0850753,1.1445378,6.17];
ybs[1244]=['',1.0633936,-0.0036735,5.38];
ybs[1245]=['',1.0508669,-0.898934,6.51];
ybs[1246]=['ν Tau',1.0659717,0.1055457,3.91];
ybs[1247]=['36 Tau',1.0718576,0.421731,5.47];
ybs[1248]=['40 Tau',1.0685139,0.0958786,5.33];
ybs[1249]=['',1.0694762,0.1440771,5.46];
ybs[1250]=['',1.0833862,0.9436133,6.31];
ybs[1251]=['37 Tau',1.0732362,0.3864048,4.36];
ybs[1252]=['',1.0702702,0.0503465,5.36];
ybs[1253]=['',1.0662229,-0.3505683,6.46];
ybs[1254]=['',1.0671019,-0.3508172,7.01];
ybs[1255]=['',1.0897936,1.088838,6.99];
ybs[1256]=['λ Per',1.0829558,0.8797832,4.29];
ybs[1257]=['39 Tau',1.0760344,0.3851246,5.9];
ybs[1258]=['',1.0695474,-0.2885229,6.39];
ybs[1259]=['γ Ret',1.0524877,-1.0838508,4.51];
ybs[1260]=['',1.0706949,-0.2222654,5.61];
ybs[1261]=['ι Ret',1.0544098,-1.0649948,4.97];
ybs[1262]=['',1.0717608,-0.3547236,6.13];
ybs[1263]=['41 Tau',1.0817965,0.4826976,5.2];
ybs[1264]=['ψ Tau',1.0836082,0.5071531,5.23];
ybs[1265]=['',1.0964867,1.0465544,6.28];
ybs[1266]=['',0.9552546,-1.4818331,6.41];
ybs[1267]=['',1.0776456,-0.1535752,6.26];
ybs[1268]=['48 Per',1.0918361,0.8337093,4.04];
ybs[1269]=['',1.0765321,-0.357011,6.34];
ybs[1270]=['',1.0755812,-0.4816163,5.59];
ybs[1271]=['',1.0955382,0.9579074,6.18];
ybs[1272]=['49 Per',1.0894566,0.6594469,6.09];
ybs[1273]=['50 Per',1.0910229,0.6648885,5.51];
ybs[1274]=['',1.0861044,0.2656187,6.01];
ybs[1275]=['',1.0874463,0.3036111,5.89];
ybs[1276]=['',1.117813,1.2597673,6.03];
ybs[1277]=['',1.1128707,1.1965121,6.32];
ybs[1278]=['ω1 Tau',1.0926637,0.3432108,5.5];
ybs[1279]=['',1.0918339,0.2348126,5.95];
ybs[1280]=['',1.0826617,-0.7480594,6.59];
ybs[1281]=['',1.101166,0.587149,5.72];
ybs[1282]=['44 Tau',1.100195,0.4631303,5.41];
ybs[1283]=['',1.0920263,-0.2850203,5.37];
ybs[1284]=['',1.1930057,1.4635143,5.57];
ybs[1285]=['37 Eri',1.0970505,-0.1198873,5.44];
ybs[1286]=['',1.0874743,-0.7995169,6.59];
ybs[1287]=['45 Tau',1.1016593,0.0973447,5.72];
ybs[1288]=['',1.0988202,-0.1529792,5.7];
ybs[1289]=['',1.0803097,-1.1199091,6.38];
ybs[1290]=['',1.1072304,0.3024885,6.09];
ybs[1291]=['',1.1205317,1.0037875,6.08];
ybs[1292]=['',1.1088602,0.3921276,6.12];
ybs[1293]=['ο1 Eri',1.103554,-0.1183916,4.04];
ybs[1294]=['',1.0976567,-0.6146903,6.44];
ybs[1295]=['',1.1019449,-0.3543334,5.79];
ybs[1296]=['',1.1145257,0.6635745,6.45];
ybs[1297]=['δ Hor',1.0976475,-0.7319719,4.93];
ybs[1298]=['μ Per',1.119128,0.8458225,4.14];
ybs[1299]=['',1.1995726,1.4553458,5.46];
ybs[1300]=['',1.1292434,1.0803916,5.7];
ybs[1301]=['52 Per',1.1185712,0.7074912,4.71];
ybs[1302]=['',1.1115974,0.1791677,6.23];
ybs[1303]=['',1.1112985,0.1561008,6.51];
ybs[1304]=['46 Tau',1.1113903,0.1356026,5.29];
ybs[1305]=['',1.1127822,0.2235163,6.25];
ybs[1306]=['47 Tau',1.1131366,0.1626084,4.84];
ybs[1307]=['',1.1114701,-0.0191357,6.44];
ybs[1308]=['',1.129783,1.0107563,5.71];
ybs[1309]=['',1.1275118,0.9366079,5.19];
ybs[1310]=['',1.1160638,0.1756493,5.22];
ybs[1311]=['',1.1048409,-0.7734321,6.71];
ybs[1312]=['',1.1819759,1.411456,5.43];
ybs[1313]=['39 Eri',1.1144776,-0.1780831,4.87];
ybs[1314]=['48 Tau',1.121341,0.2697031,6.32];
ybs[1315]=['μ Tau',1.1200846,0.1561134,4.29];
ybs[1316]=['',1.1195365,0.1091214,6.93];
ybs[1317]=['',1.1197834,0.108893,6.31];
ybs[1318]=['',1.1097378,-0.7034437,6.37];
ybs[1319]=['',1.1338851,0.8787137,4.61];
ybs[1320]=['ο2 Eri',1.1183896,-0.1326487,4.43];
ybs[1321]=['α Hor',1.1113988,-0.7372481,3.86];
ybs[1322]=['',1.1463079,1.1377861,5.27];
ybs[1323]=['',1.1328527,0.7363932,6.22];
ybs[1324]=['ω2 Tau',1.1280253,0.3600648,4.94];
ybs[1325]=['',1.1381385,0.8744005,5.45];
ybs[1326]=['51 Tau',1.1329809,0.3775185,5.65];
ybs[1327]=['',1.1273654,-0.1120561,5.94];
ybs[1328]=['',1.1424468,0.8896107,5.55];
ybs[1329]=['',1.1326438,0.1664648,6.54];
ybs[1330]=['',1.1504978,1.0608956,5.39];
ybs[1331]=['α Ret',1.1113588,-1.0894476,3.35];
ybs[1332]=['',1.1420175,0.7305629,5.92];
ybs[1333]=['γ Dor',1.1195764,-0.8976986,4.25];
ybs[1334]=['53 Tau',1.13754,0.3698837,5.35];
ybs[1335]=['',1.1130855,-1.0845299,5.45];
ybs[1336]=['56 Tau',1.1383353,0.380902,5.38];
ybs[1337]=['',1.1502821,0.9870876,5.88];
ybs[1338]=['54 Per',1.1423698,0.6041758,4.93];
ybs[1339]=['',1.1411834,0.5585667,6.16];
ybs[1340]=['',1.1309885,-0.3606614,6];
ybs[1341]=['γ Tau',1.1389013,0.2736308,3.65];
ybs[1342]=['υ4 Eri',1.1288255,-0.5889953,3.56];
ybs[1343]=['φ Tau',1.1418035,0.4782366,4.95];
ybs[1344]=['',1.1379749,0.1775327,6.31];
ybs[1345]=['53 Per',1.148056,0.8124223,4.85];
ybs[1346]=['57 Tau',1.1395796,0.24584,5.59];
ybs[1347]=['',1.1554669,1.0413529,6.19];
ybs[1348]=['',1.1324653,-0.4000165,6.07];
ybs[1349]=['',1.1417421,0.3279926,6.12];
ybs[1350]=['ε Ret',1.1207567,-1.0341039,4.44];
ybs[1351]=['58 Tau',1.1424247,0.2643351,5.26];
ybs[1352]=['',1.119973,-1.0628423,6.37];
ybs[1353]=['',1.1435814,0.242846,6.17];
ybs[1354]=['',1.1338571,-0.5908665,6.37];
ybs[1355]=['',1.142473,0.1078757,5.77];
ybs[1356]=['',1.1431447,0.1618876,6.53];
ybs[1357]=['',1.1418796,-0.1081323,6.27];
ybs[1358]=['',1.1421324,-0.1316414,5.85];
ybs[1359]=['',1.1342735,-0.7717372,5.34];
ybs[1360]=['',1.1309622,-0.9216893,6.09];
ybs[1361]=['',1.1456038,-0.0008451,5.86];
ybs[1362]=['',1.1413944,-0.3593575,5.38];
ybs[1363]=['60 Tau',1.1487315,0.2465547,5.72];
ybs[1364]=['χ Tau',1.151452,0.4481694,5.37];
ybs[1365]=['',1.1503821,0.3642597,5.91];
ybs[1366]=['',1.1567199,0.7413561,6.23];
ybs[1367]=['θ Ret',1.1253695,-1.1031166,5.87];
ybs[1368]=['δ1 Tau',1.1526822,0.3070278,3.76];
ybs[1369]=['',1.1449966,-0.4481775,6.01];
ybs[1370]=['',1.1554514,0.3670572,5.99];
ybs[1371]=['63 Tau',1.1547564,0.2936672,5.64];
ybs[1372]=['55 Per',1.1601466,0.5965304,5.73];
ybs[1373]=['62 Tau',1.1575747,0.4249787,6.36];
ybs[1374]=['56 Per',1.160735,0.5935477,5.76];
ybs[1375]=['δ2 Tau',1.1577486,0.305297,4.8];
ybs[1376]=['66 Tau',1.1564555,0.1659686,5.12];
ybs[1377]=['',1.1729115,1.0058697,6.32];
ybs[1378]=['ξ Eri',1.1552049,-0.0645244,5.17];
ybs[1379]=['',1.1518928,-0.4335977,5.83];
ybs[1380]=['',1.1615394,0.3331763,5.98];
ybs[1381]=['',1.1515609,-0.6195236,6.39];
ybs[1382]=['κ1 Tau',1.1634823,0.3899347,4.22];
ybs[1383]=['κ2 Tau',1.1636896,0.3882907,5.28];
ybs[1384]=['δ3 Tau',1.1638476,0.3137357,4.29];
ybs[1385]=['',1.1670842,0.5495385,5.28];
ybs[1386]=['70 Tau',1.1643496,0.2790511,6.46];
ybs[1387]=['υ Tau',1.1676056,0.3989977,4.28];
ybs[1388]=['43 Eri',1.1555974,-0.5928617,3.96];
ybs[1389]=['71 Tau',1.1675021,0.2734163,4.49];
ybs[1390]=['η Ret',1.1437216,-1.1054349,5.24];
ybs[1391]=['π Tau',1.1686042,0.2576239,4.69];
ybs[1392]=['',1.1672772,0.1507538,6.06];
ybs[1393]=['',1.1595003,-0.6057992,6.55];
ybs[1394]=['72 Tau',1.1719049,0.4021795,5.53];
ybs[1395]=['',1.1699368,0.0371131,6.23];
ybs[1396]=['',1.2044197,1.2666222,5.94];
ybs[1397]=['',1.1722919,0.1965106,5.88];
ybs[1398]=['',1.1750009,0.378151,5.72];
ybs[1399]=['',1.1606199,-0.769916,6.39];
ybs[1400]=['',1.1546661,-0.9952374,6.29];
ybs[1401]=['',1.1790837,0.5307095,6.4];
ybs[1402]=['75 Tau',1.176665,0.2863384,4.97];
ybs[1403]=['76 Tau',1.1763881,0.258084,5.9];
ybs[1404]=['ε Tau',1.1775407,0.3355648,3.53];
ybs[1405]=['',1.1687396,-0.4194782,6.11];
ybs[1406]=['θ1 Tau',1.1772397,0.2793996,3.84];
ybs[1407]=['θ2 Tau',1.1776146,0.2778039,3.4];
ybs[1408]=['',1.1745037,0.03325,6.15];
ybs[1409]=['79 Tau',1.1782769,0.228526,5.03];
ybs[1410]=['',1.1765603,0.0249072,5.55];
ybs[1411]=['',1.1579923,-1.0679711,5.94];
ybs[1412]=['1 Cam',1.1944265,0.9416966,5.77];
ybs[1413]=['',1.1681972,-0.8185673,6.1];
ybs[1414]=['',1.1869222,0.5672881,6.21];
ybs[1415]=['',1.1820213,0.1844346,6.79];
ybs[1416]=['',1.1763349,-0.3388052,5.96];
ybs[1417]=['80 Tau',1.1840726,0.2737287,5.58];
ybs[1418]=['',1.1785976,-0.2269337,5.6];
ybs[1419]=['',1.1906416,0.6990925,6.26];
ybs[1420]=['',1.1834305,0.1799085,6.48];
ybs[1421]=['δ Men',1.1195934,-1.3990913,5.69];
ybs[1422]=['',1.1859186,0.2834262,4.78];
ybs[1423]=['81 Tau',1.1862785,0.2746649,5.48];
ybs[1424]=['',1.173133,-0.7315282,6.44];
ybs[1425]=['83 Tau',1.1860911,0.2403258,5.4];
ybs[1426]=['',1.1831443,-0.2364351,6.24];
ybs[1427]=['85 Tau',1.191587,0.2774424,6.02];
ybs[1428]=['',1.1779644,-0.8110423,6.16];
ybs[1429]=['57 Per',1.1996246,0.7523708,6.09];
ybs[1430]=['',1.1694603,-1.0903772,5.75];
ybs[1431]=['',1.1921655,0.0951994,6.39];
ybs[1432]=['45 Eri',1.1911004,0.000013,4.91];
ybs[1433]=['',1.1886828,-0.237367,6.21];
ybs[1434]=['',1.18445,-0.6214824,5.96];
ybs[1435]=['',1.2149041,1.1223143,5.94];
ybs[1436]=['',1.1942437,-0.0552377,5.81];
ybs[1437]=['',1.1990211,0.3152144,6.25];
ybs[1438]=['δ Cae',1.1846192,-0.7838038,5.07];
ybs[1439]=['ρ Tau',1.2002133,0.2598462,4.65];
ybs[1440]=['',1.2042011,0.5062213,5.88];
ybs[1441]=['',1.1998164,0.1650512,6.01];
ybs[1442]=['',1.1972257,-0.1874816,6.06];
ybs[1443]=['',1.2011527,0.0979502,5.68];
ybs[1444]=['46 Eri',1.1997514,-0.1168539,5.72];
ybs[1445]=['',1.2011587,-0.1185825,6.09];
ybs[1446]=['47 Eri',1.2009214,-0.1429052,5.11];
ybs[1447]=['',1.200903,-0.1558013,5.26];
ybs[1448]=['υ1 Eri',1.19711,-0.5187601,4.51];
ybs[1449]=['58 Per',1.2138126,0.7209417,4.25];
ybs[1450]=['',1.2085337,0.3477462,6.36];
ybs[1451]=['ν Men',1.1307291,-1.4229571,5.79];
ybs[1452]=['α Tau',1.2093155,0.2888834,0.85];
ybs[1453]=['88 Tau',1.2079303,0.1780866,4.25];
ybs[1454]=['',1.2120528,0.4081134,6.02];
ybs[1455]=['',1.2054063,-0.1691861,6.37];
ybs[1456]=['',1.2040529,-0.3469262,6.13];
ybs[1457]=['',1.2090791,-0.0622964,6.33];
ybs[1458]=['ν Eri',1.2103606,-0.0577708,3.93];
ybs[1459]=['υ2 Eri',1.2059711,-0.5326625,3.82];
ybs[1460]=['α Dor',1.1975647,-0.9599528,3.27];
ybs[1461]=['2 Cam',1.2290562,0.9339884,5.35];
ybs[1462]=['3 Cam',1.2287709,0.9271239,5.05];
ybs[1463]=['',1.2610937,1.3377659,6.49];
ybs[1464]=['',1.2144797,0.0181577,5.31];
ybs[1465]=['',1.2209617,0.4709134,6.47];
ybs[1466]=['',1.2197016,0.3617404,5.92];
ybs[1467]=['89 Tau',1.2190593,0.2805595,5.79];
ybs[1468]=['90 Tau',1.2189383,0.2190803,4.27];
ybs[1469]=['51 Eri',1.2159901,-0.0424374,5.23];
ybs[1470]=['',1.1946876,-1.0957106,5.79];
ybs[1471]=['',1.2116265,-0.5353691,6.3];
ybs[1472]=['',1.2247805,0.4408571,6.22];
ybs[1473]=['σ1 Tau',1.2234003,0.2764737,5.07];
ybs[1474]=['σ2 Tau',1.2239357,0.278538,4.69];
ybs[1475]=['',1.2228981,0.1380892,5.39];
ybs[1476]=['53 Eri',1.2180986,-0.2489241,3.87];
ybs[1477]=['',1.2348671,0.8437042,5.67];
ybs[1478]=['',1.2212892,-0.2108675,5.01];
ybs[1479]=['93 Tau',1.2272113,0.2136003,5.46];
ybs[1480]=['',1.1365041,-1.4459942,6.76];
ybs[1481]=['',1.2443047,1.0395126,6.5];
ybs[1482]=['',1.2231062,-0.2498987,5.45];
ybs[1483]=['',1.2255723,-0.0176627,6.1];
ybs[1484]=['',1.2361164,0.6688094,5.99];
ybs[1485]=['',1.233417,0.5001233,5.78];
ybs[1486]=['',1.2730732,1.326048,6.06];
ybs[1487]=['',1.2087192,-1.0827157,5.4];
ybs[1488]=['',1.2435534,0.8728875,5.87];
ybs[1489]=['59 Per',1.24108,0.757545,5.29];
ybs[1490]=['',1.2261392,-0.4265905,5.58];
ybs[1491]=['54 Eri',1.2277637,-0.3426287,4.32];
ybs[1492]=['τ Tau',1.2371703,0.401364,4.28];
ybs[1493]=['',1.2200091,-0.9011399,6.44];
ybs[1494]=['95 Tau',1.2415159,0.4211116,6.13];
ybs[1495]=['',1.2466378,0.7125383,6.08];
ybs[1496]=['',1.2444191,0.5742831,6.45];
ybs[1497]=['α Cae',1.2272288,-0.729956,4.45];
ybs[1498]=['β Cae',1.234028,-0.6475993,5.05];
ybs[1499]=['',1.2246096,-1.0280497,6.53];
ybs[1500]=['55 Eri',1.2418398,-0.1527982,6.82];
ybs[1501]=['55 Eri',1.2418761,-0.1528418,6.7];
ybs[1502]=['',1.2462586,0.1952076,5.4];
ybs[1503]=['56 Eri',1.2440824,-0.1477411,5.9];
ybs[1504]=['',1.2391307,-0.5362764,5.68];
ybs[1505]=['',1.2786062,1.2387781,6.37];
ybs[1506]=['4 Cam',1.2645042,0.9912381,5.3];
ybs[1507]=['',1.2523138,0.4130472,6.35];
ybs[1508]=['',1.2439064,-0.3251199,5.53];
ybs[1509]=['',1.257633,0.7042407,5.97];
ybs[1510]=['',1.2648623,0.9710836,6.26];
ybs[1511]=['λ Pic',1.23625,-0.8803781,5.31];
ybs[1512]=['',1.2546137,0.3276424,6.01];
ybs[1513]=['',1.2411228,-0.7160352,6.25];
ybs[1514]=['',1.2532536,0.204958,5.37];
ybs[1515]=['μ Eri',1.2504321,-0.0561431,4.02];
ybs[1516]=['',1.2478901,-0.3708022,5.72];
ybs[1517]=['',1.2543693,-0.0509101,6.33];
ybs[1518]=['',1.3283944,1.4176203,5.07];
ybs[1519]=['',1.2506262,-0.5928383,6.86];
ybs[1520]=['',1.2535369,-0.489564,6.19];
ybs[1521]=['',1.2507596,-0.686243,6.05];
ybs[1522]=['',1.2832934,1.1089774,5.44];
ybs[1523]=['',1.2684647,0.5694018,5.86];
ybs[1524]=['',1.2679579,0.549312,5.58];
ybs[1525]=['κ Dor',1.2421416,-1.0418579,5.27];
ybs[1526]=['',1.2104654,-1.3546222,6.05];
ybs[1527]=['58 Eri',1.2591173,-0.294917,5.51];
ybs[1528]=['',1.2713096,0.6549175,4.88];
ybs[1529]=['',1.264814,0.0632624,6.03];
ybs[1530]=['',1.27749,0.8512991,5.66];
ybs[1531]=['',1.2638879,-0.0983923,5.78];
ybs[1532]=['96 Tau',1.2695841,0.2782053,6.08];
ybs[1533]=['59 Eri',1.2632263,-0.2843659,5.77];
ybs[1534]=['ζ Cae',1.2595413,-0.5233093,6.37];
ybs[1535]=['',1.2442841,-1.1028959,6.46];
ybs[1536]=['μ Men',1.2342132,-1.2372923,5.54];
ybs[1537]=['α Cam',1.2924259,1.1584834,4.29];
ybs[1538]=['π3 Ori',1.2697177,0.1221238,3.19];
ybs[1539]=['π2 Ori',1.2731538,0.155957,4.36];
ybs[1540]=['',1.2683896,-0.2397005,6.26];
ybs[1541]=['',1.286619,0.9228397,6.41];
ybs[1542]=['97 Tau',1.2768537,0.329426,5.1];
ybs[1543]=['',1.2619847,-0.766958,6.72];
ybs[1544]=['60 Eri',1.270436,-0.2824216,5.03];
ybs[1545]=['',1.2842157,0.7438748,5.71];
ybs[1546]=['2 Aur',1.2831548,0.6411882,4.78];
ybs[1547]=['π4 Ori',1.2756335,0.0984384,3.69];
ybs[1548]=['',1.2780434,0.1747046,6.11];
ybs[1549]=['',1.2833842,0.4875013,5.97];
ybs[1550]=['5 Cam',1.2951112,0.9650308,5.52];
ybs[1551]=['ο1 Ori',1.2817344,0.2493199,4.74];
ybs[1552]=['',1.2695978,-0.7205618,6.07];
ybs[1553]=['',1.2932993,0.7695857,6.08];
ybs[1554]=['',1.275193,-0.6086194,5.86];
ybs[1555]=['ω Eri',1.2826131,-0.0945705,4.39];
ybs[1556]=['',1.2995252,0.9233132,5.75];
ybs[1557]=['5 Ori',1.2850081,0.0443677,5.33];
ybs[1558]=['ι Pic',1.2714948,-0.9324592,5.61];
ybs[1559]=['ι Pic',1.2715749,-0.9324303,6.42];
ybs[1560]=['',1.2873752,0.0279811,6.61];
ybs[1561]=['',1.2925756,0.3406616,6.37];
ybs[1562]=['π5 Ori',1.2888092,0.0431819,3.72];
ybs[1563]=['7 Cam',1.3047138,0.9387102,4.47];
ybs[1564]=['6 Ori',1.2914419,0.2000046,5.19];
ybs[1565]=['π1 Ori',1.2918974,0.1777458,4.65];
ybs[1566]=['',1.2913759,0.1363534,5.33];
ybs[1567]=['',1.3311947,1.2967495,6.06];
ybs[1568]=['',1.2992698,0.6318328,6.07];
ybs[1569]=['',1.291329,0.0087405,5.99];
ybs[1570]=['',1.298407,0.429783,6.37];
ybs[1571]=['',1.2961868,0.2630694,5.81];
ybs[1572]=['ι Aur',1.302002,0.5794188,2.69];
ybs[1573]=['',1.2964259,0.0948042,6.5];
ybs[1574]=['',1.2918798,-0.2915984,5.7];
ybs[1575]=['ο2 Ori',1.2984572,0.2364389,4.07];
ybs[1576]=['',1.2927502,-0.2859667,5.72];
ybs[1577]=['62 Eri',1.2979302,-0.0896902,5.51];
ybs[1578]=['',1.2932177,-0.4484581,6.72];
ybs[1579]=['',1.2899498,-0.6910675,6.1];
ybs[1580]=['',1.3029588,0.2999453,5.48];
ybs[1581]=['99 Tau',1.3051519,0.4185363,5.79];
ybs[1582]=['',1.3391825,1.2879142,6.66];
ybs[1583]=['8 Cam',1.3155052,0.9282746,6.08];
ybs[1584]=['',1.3412751,1.2931994,5.96];
ybs[1585]=['98 Tau',1.3067045,0.4377609,5.81];
ybs[1586]=['',1.301929,-0.0180668,6.23];
ybs[1587]=['ω Aur',1.3121328,0.6618509,4.94];
ybs[1588]=['',1.3244442,1.0665312,6.03];
ybs[1589]=['',1.3309345,1.1667835,6.19];
ybs[1590]=['',1.3034629,-0.2478284,6.15];
ybs[1591]=['',1.3057865,-0.0380635,6.35];
ybs[1592]=['',1.2882023,-1.0212577,6.12];
ybs[1593]=['',1.2808408,-1.1631105,6.41];
ybs[1594]=['5 Aur',1.316796,0.6880991,5.95];
ybs[1595]=['',1.3099219,0.2543636,6.09];
ybs[1596]=['π6 Ori',1.3075318,0.0304665,4.47];
ybs[1597]=['6 Aur',1.3171679,0.692641,6.58];
ybs[1598]=['β Cam',1.3322982,1.0554178,4.03];
ybs[1599]=['',1.3089442,-0.2852671,5.66];
ybs[1600]=['ε Aur',1.3243372,0.765378,2.99];
ybs[1601]=['',1.2773598,-1.2631465,6.28];
ybs[1602]=['',1.3115632,-0.2578707,7.71];
ybs[1603]=['63 Eri',1.3127419,-0.1785914,5.38];
ybs[1604]=['',1.316304,0.0636294,7.03];
ybs[1605]=['',1.3163985,0.0636437,6.66];
ybs[1606]=['64 Eri',1.3130512,-0.2182838,4.79];
ybs[1607]=['ζ Aur',1.326385,0.7174207,3.75];
ybs[1608]=['',1.316627,-0.0355206,6.32];
ybs[1609]=['',1.3171653,-0.0998857,6.22];
ybs[1610]=['',1.3300383,0.7237983,6.14];
ybs[1611]=['',1.4818629,1.4963655,6.51];
ybs[1612]=['ψ Eri',1.319827,-0.1246846,4.81];
ybs[1613]=['',1.3218523,0.0131248,5.92];
ybs[1614]=['',1.3225891,0.0285985,6.24];
ybs[1615]=['ι Tau',1.328112,0.3773243,4.64];
ybs[1616]=['',1.3192924,-0.3494485,4.91];
ybs[1617]=['11 Cam',1.3439858,1.0297424,5.08];
ybs[1618]=['12 Cam',1.3442616,1.0305903,6.08];
ybs[1619]=['',1.3458429,1.0680926,6.04];
ybs[1620]=['',1.3256919,-0.0729618,5.85];
ybs[1621]=['',1.3335084,0.5327305,6.14];
ybs[1622]=['',1.3352267,0.5645891,6.62];
ybs[1623]=['',1.3222459,-0.4580675,5.02];
ybs[1624]=['η Men',1.2853603,-1.3073108,5.47];
ybs[1625]=['',1.344428,0.9500376,7.24];
ybs[1626]=['',1.3189981,-0.6926871,6.03];
ybs[1627]=['',1.3350807,0.4838821,6.6];
ybs[1628]=['',1.3336199,0.3718686,6.19];
ybs[1629]=['1 Lep',1.3249464,-0.3973353,5.75];
ybs[1630]=['',1.3229367,-0.5539992,5.94];
ybs[1631]=['',1.3612476,1.2158827,6.41];
ybs[1632]=['9 Aur',1.3455055,0.9010253,5];
ybs[1633]=['11 Ori',1.3342979,0.2693483,4.68];
ybs[1634]=['',1.3415247,0.6276895,6.52];
ybs[1635]=['',1.3301576,-0.2502919,6.41];
ybs[1636]=['η Aur',1.3440158,0.7201532,3.17];
ybs[1637]=['',1.3386883,0.3461729,6.44];
ybs[1638]=['',1.3749121,1.2910316,5.43];
ybs[1639]=['',1.3455019,0.7540146,6.2];
ybs[1640]=['',1.3298524,-0.4251446,5.61];
ybs[1641]=['',1.3351144,-0.0525606,6.05];
ybs[1642]=['',1.3606748,1.1335032,6.41];
ybs[1643]=['',1.3373874,0.0210394,6.17];
ybs[1644]=['η1 Pic',1.3237252,-0.8573402,5.38];
ybs[1645]=['',1.3860128,1.3350993,6.37];
ybs[1646]=['',1.3290135,-0.7280845,6.31];
ybs[1647]=['γ1 Cae',1.3315872,-0.6188026,4.55];
ybs[1648]=['γ2 Cae',1.331699,-0.6226765,6.34];
ybs[1649]=['ε Lep',1.3368007,-0.3899611,3.19];
ybs[1650]=['',1.3358032,-0.4559569,5.73];
ybs[1651]=['104 Tau',1.3469992,0.3258858,5];
ybs[1652]=['66 Eri',1.3431432,-0.0807689,5.12];
ybs[1653]=['106 Tau',1.3486266,0.3568331,5.3];
ybs[1654]=['103 Tau',1.3501105,0.4239721,5.5];
ybs[1655]=['105 Tau',1.3491956,0.3792838,5.89];
ybs[1656]=['',1.3421803,-0.2285432,6.05];
ybs[1657]=['13 Ori',1.3474738,0.1657845,6.17];
ybs[1658]=['η2 Pic',1.3330943,-0.8648011,5.03];
ybs[1659]=['14 Ori',1.3485001,0.1487897,5.34];
ybs[1660]=['',1.345716,-0.2175306,5.97];
ybs[1661]=['β Eri',1.3478763,-0.0883076,2.79];
ybs[1662]=['',1.3328309,-0.9490955,6.27];
ybs[1663]=['',1.3627464,0.8200841,5.68];
ybs[1664]=['',1.3603688,0.6514846,6.02];
ybs[1665]=['',1.3574416,0.4896739,6.01];
ybs[1666]=['',1.3498706,-0.1507704,5.78];
ybs[1667]=['16 Ori',1.3548541,0.1720091,5.43];
ybs[1668]=['68 Eri',1.3517239,-0.077315,5.12];
ybs[1669]=['ζ Dor',1.3346643,-1.0025988,4.72];
ybs[1670]=['',1.374589,1.0799022,6.17];
ybs[1671]=['15 Ori',1.3566993,0.2726722,4.82];
ybs[1672]=['β Men',1.3196795,-1.244154,5.31];
ybs[1673]=['14 Cam',1.3767731,1.0945829,6.5];
ybs[1674]=['λ Eri',1.3534014,-0.1523338,4.27];
ybs[1675]=['',1.348324,-0.6229382,6.52];
ybs[1676]=['',1.3576513,-0.0094192,6.1];
ybs[1677]=['',1.3050378,-1.366054,6.29];
ybs[1678]=['',1.400158,1.2791355,5.74];
ybs[1679]=['',1.3654151,0.2804794,5.18];
ybs[1680]=['',1.3615761,-0.0388991,6.25];
ybs[1681]=['',1.4230951,1.3831674,5.05];
ybs[1682]=['',1.3631095,-0.0430376,5.9];
ybs[1683]=['',1.3835469,1.0372199,6.15];
ybs[1684]=['μ Aur',1.3740256,0.6720953,4.86];
ybs[1685]=['',1.3648311,0.0094158,6.67];
ybs[1686]=['',1.3651333,0.0185297,5.89];
ybs[1687]=['',1.3808559,0.9291598,6.2];
ybs[1688]=['',1.3630315,-0.2063716,5.68];
ybs[1689]=['',1.359683,-0.4517634,6.41];
ybs[1690]=['',1.3427599,-1.1060609,5.2];
ybs[1691]=['ι Lep',1.3670373,-0.2067288,4.45];
ybs[1692]=['',1.3694534,-0.105296,5.91];
ybs[1693]=['ρ Ori',1.3719036,0.0503536,4.46];
ybs[1694]=['',1.362851,-0.6522361,6.57];
ybs[1695]=['',1.3340075,-1.2742615,6.27];
ybs[1696]=['',1.3728971,0.0347647,6.09];
ybs[1697]=['μ Lep',1.3696359,-0.2824186,3.31];
ybs[1698]=['',1.3739811,0.0101922,6.32];
ybs[1699]=['',1.3726572,-0.1417896,6.37];
ybs[1700]=['κ Lep',1.3710692,-0.2254509,4.36];
ybs[1701]=['14 Aur',1.3823459,0.5709069,5.02];
ybs[1702]=['',1.3920622,0.9356333,6.5];
ybs[1703]=['α Aur',1.3887558,0.8032028,0.08];
ybs[1704]=['',1.3782839,0.0903959,5.5];
ybs[1705]=['',1.374351,-0.2545223,6.21];
ybs[1706]=['108 Tau',1.3821103,0.3893394,6.27];
ybs[1707]=['',1.3863431,0.5992459,5.96];
ybs[1708]=['β Ori',1.3769456,-0.142739,0.12];
ybs[1709]=['',1.5325287,1.4910315,6.6];
ybs[1710]=['',1.3724324,-0.6248588,6.98];
ybs[1711]=['ξ Men',1.2935882,-1.4388214,5.85];
ybs[1712]=['',1.380543,-0.0241947,6.15];
ybs[1713]=['18 Ori',1.3843283,0.1983373,5.56];
ybs[1714]=['15 Cam',1.4020397,1.0146961,6.13];
ybs[1715]=['',1.406699,1.0938621,5.61];
ybs[1716]=['',1.3754994,-0.6275122,5.76];
ybs[1717]=['',1.395405,0.7472366,5.48];
ybs[1718]=['',1.3799803,-0.4698498,5.07];
ybs[1719]=['',1.3866773,0.0343729,6.42];
ybs[1720]=['',1.3970376,0.7066154,6.18];
ybs[1721]=['16 Aur',1.3944775,0.5828183,4.54];
ybs[1722]=['',1.3717819,-0.9076989,6.05];
ybs[1723]=['17 Aur',1.3950944,0.5897207,6.14];
ybs[1724]=['λ Aur',1.3990592,0.7002262,4.71];
ybs[1725]=['',1.3812412,-0.6091882,6.66];
ybs[1726]=['',1.3864833,-0.2987913,6.56];
ybs[1727]=['',1.398083,0.5893849,5.41];
ybs[1728]=['',1.4032573,0.7757276,6.62];
ybs[1729]=['18 Aur',1.399812,0.5935216,6.49];
ybs[1730]=['τ Ori',1.3903824,-0.1190786,3.6];
ybs[1731]=['',1.4061361,0.8200242,6.54];
ybs[1732]=['',1.3904114,-0.2355845,5.5];
ybs[1733]=['',1.4039437,0.7174417,5.52];
ybs[1734]=['109 Tau',1.398751,0.3860182,4.94];
ybs[1735]=['19 Aur',1.4025238,0.5930362,5.03];
ybs[1736]=['',1.3985318,0.351781,6.08];
ybs[1737]=['',1.3794431,-0.9103472,6.49];
ybs[1738]=['ο Col',1.3886657,-0.6086556,4.83];
ybs[1739]=['θ Dor',1.3689672,-1.1721856,4.83];
ybs[1740]=['',1.4519483,1.3612269,6.56];
ybs[1741]=['21 Ori',1.3976163,0.0456708,5.34];
ybs[1742]=['',1.3953371,-0.3160593,5.96];
ybs[1743]=['',1.3992176,-0.0242863,6.34];
ybs[1744]=['ρ Aur',1.4107985,0.7299649,5.23];
ybs[1745]=['',1.4064777,0.4882934,6.33];
ybs[1746]=['16 Cam',1.4194295,1.0046633,5.28];
ybs[1747]=['',1.4075278,0.5164396,5.76];
ybs[1748]=['',1.3972921,-0.3228701,6.36];
ybs[1749]=['',1.3973508,-0.3226908,6.54];
ybs[1750]=['',1.4059301,0.346171,6.18];
ybs[1751]=['λ Lep',1.3987358,-0.2296141,4.29];
ybs[1752]=['ν Lep',1.4005571,-0.2145886,5.3];
ybs[1753]=['',1.3973673,-0.4773129,5.99];
ybs[1754]=['',1.4027983,-0.0933119,6.39];
ybs[1755]=['',1.4152517,0.7164294,5.54];
ybs[1756]=['',1.4069832,0.0703674,6.57];
ybs[1757]=['',1.4022279,-0.3703436,4.71];
ybs[1758]=['',1.4089112,0.147449,5.8];
ybs[1759]=['',1.4077328,-0.0069233,5.68];
ybs[1760]=['22 Ori',1.4087448,-0.0063339,4.73];
ybs[1761]=['',1.4011448,-0.6052537,6.34];
ybs[1762]=['ζ Pic',1.3958117,-0.8828772,5.45];
ybs[1763]=['22 Aur',1.416968,0.5053663,6.46];
ybs[1764]=['',1.4086394,-0.2397477,6.56];
ybs[1765]=['23 Ori',1.4135624,0.0621944,5];
ybs[1766]=['',1.407844,-0.4320285,5.06];
ybs[1767]=['',1.4052567,-0.5990904,6.09];
ybs[1768]=['σ Aur',1.4229542,0.6528153,4.99];
ybs[1769]=['110 Tau',1.4175226,0.2917849,6.08];
ybs[1770]=['',1.4225624,0.5452745,6.28];
ybs[1771]=['',1.422573,0.5438637,5.94];
ybs[1772]=['',1.4166157,0.0932211,6.35];
ybs[1773]=['',1.4152019,-0.1465556,5.9];
ybs[1774]=['',1.4252729,0.6086488,6.55];
ybs[1775]=['111 Tau',1.4210187,0.3037138,4.99];
ybs[1776]=['',1.4172322,-0.0024632,5.7];
ybs[1777]=['',1.4178685,-0.0148127,6.11];
ybs[1778]=['8 Lep',1.4158464,-0.2427491,5.25];
ybs[1779]=['29 Ori',1.4180094,-0.1359537,4.14];
ybs[1780]=['',1.4139967,-0.4657746,6.49];
ybs[1781]=['',1.4212425,0.04138,6.32];
ybs[1782]=['27 Ori',1.4205947,-0.0152402,5.08];
ybs[1783]=['η Ori',1.4205185,-0.041517,3.36];
ybs[1784]=['ψ1 Ori',1.4218496,0.0325406,4.95];
ybs[1785]=['γ Ori',1.4236928,0.111135,1.64];
ybs[1786]=['β Tau',1.429652,0.4995951,1.65];
ybs[1787]=['',1.4199601,-0.2959708,5.65];
ybs[1788]=['',1.4141959,-0.6921934,5.71];
ybs[1789]=['',1.4326736,0.6191394,6.15];
ybs[1790]=['',1.4322236,0.6005428,5.94];
ybs[1791]=['',1.4323392,0.5808397,6.15];
ybs[1792]=['',1.4154241,-0.6513212,6.82];
ybs[1793]=['113 Tau',1.4282865,0.2917771,6.25];
ybs[1794]=['',1.4226284,-0.179965,5.61];
ybs[1795]=['',1.4251377,-0.0091893,6.57];
ybs[1796]=['κ Pic',1.4083349,-0.9793912,6.11];
ybs[1797]=['17 Cam',1.4495274,1.1009926,5.42];
ybs[1798]=['',1.4263252,0.0093961,6.16];
ybs[1799]=['',1.4334231,0.527532,5.74];
ybs[1800]=['φ Aur',1.435865,0.6020043,5.07];
ybs[1801]=['',1.4272266,-0.0960094,6.23];
ybs[1802]=['',1.4303225,0.1201873,6.42];
ybs[1803]=['115 Tau',1.4330202,0.3137925,5.42];
ybs[1804]=['',1.4331829,0.2665906,6.16];
ybs[1805]=['114 Tau',1.4352192,0.3831601,4.88];
ybs[1806]=['ψ2 Ori',1.4310142,0.0543238,4.59];
ybs[1807]=['',1.4264974,-0.3434477,5.89];
ybs[1808]=['',1.4205587,-0.7715708,6.08];
ybs[1809]=['116 Tau',1.4355196,0.2773438,5.5];
ybs[1810]=['',1.3543207,-1.422731,6.51];
ybs[1811]=['117 Tau',1.4367369,0.3011603,5.77];
ybs[1812]=['',1.4315282,-0.2074142,6.35];
ybs[1813]=['θ Pic',1.4192123,-0.9127755,6.27];
ybs[1814]=['',1.439013,0.2390216,6.35];
ybs[1815]=['',1.4361416,0.0229457,6.41];
ybs[1816]=['118 Tau',1.4425141,0.4392332,5.47];
ybs[1817]=['',1.4444514,0.509668,6.24];
ybs[1818]=['',1.4334597,-0.3727837,6.07];
ybs[1819]=['',1.4500711,0.7239059,6];
ybs[1820]=['',1.4497137,0.695351,6.37];
ybs[1821]=['',1.439982,-0.0574494,6.39];
ybs[1822]=['',1.4302066,-0.7143049,5.87];
ybs[1823]=['18 Cam',1.4591128,0.998938,6.48];
ybs[1824]=['β Lep',1.4362629,-0.3620363,2.84];
ybs[1825]=['',1.4419332,-0.0598775,5.79];
ybs[1826]=['',1.4487182,0.3923049,6.29];
ybs[1827]=['',1.4471724,0.2683506,5.94];
ybs[1828]=['',1.4443916,0.0314952,5.78];
ybs[1829]=['31 Ori',1.443502,-0.0187927,4.71];
ybs[1830]=['',1.4355185,-0.6495107,5.57];
ybs[1831]=['λ Dor',1.4252283,-1.0279126,5.14];
ybs[1832]=['',1.4462972,0.0736412,6.21];
ybs[1833]=['',1.4396278,-0.525358,6.75];
ybs[1834]=['32 Ori',1.4483461,0.1040736,4.2];
ybs[1835]=['',1.44594,-0.1294955,6.33];
ybs[1836]=['33 Ori',1.450242,0.0577164,5.46];
ybs[1837]=['χ Aur',1.457916,0.5620971,4.76];
ybs[1838]=['',1.4949495,1.3099343,6.17];
ybs[1839]=['119 Tau',1.4550604,0.3247813,4.38];
ybs[1840]=['',1.4617585,0.7351731,6.55];
ybs[1841]=['',1.4551008,0.297966,5.46];
ybs[1842]=['',1.4503443,-0.1168267,6.22];
ybs[1843]=['10 Lep',1.4488309,-0.3638804,5.55];
ybs[1844]=['',1.4611328,0.5727224,6.48];
ybs[1845]=['δ Ori',1.4534576,-0.004715,6.85];
ybs[1846]=['δ Ori',1.4534499,-0.0049718,2.23];
ybs[1847]=['',1.4812557,1.1642635,6.26];
ybs[1848]=['',1.4619973,0.6063134,6.27];
ybs[1849]=['υ Ori',1.4528604,-0.1271828,4.62];
ybs[1850]=['',1.4431684,-0.8213936,5.46];
ybs[1851]=['19 Cam',1.4806406,1.1199087,6.15];
ybs[1852]=['120 Tau',1.4607891,0.323824,5.69];
ybs[1853]=['',1.4263355,-1.1973931,6.03];
ybs[1854]=['',1.4613917,0.3575755,6.18];
ybs[1855]=['',1.4563774,-0.0275413,5.35];
ybs[1856]=['ε Col',1.4485149,-0.6188197,3.87];
ybs[1857]=['',1.4582563,-0.0297511,6.46];
ybs[1858]=['35 Ori',1.4622732,0.2499108,5.64];
ybs[1859]=['α Lep',1.4559498,-0.3108128,2.58];
ybs[1860]=['',1.4763521,0.9501683,5.73];
ybs[1861]=['',1.437662,-1.087314,6.59];
ybs[1862]=['',1.4600365,-0.019942,5.34];
ybs[1863]=['',1.4743271,0.8329973,6.11];
ybs[1864]=['',1.4495523,-0.8012919,5.86];
ybs[1865]=['',1.4620347,0.0248022,6.59];
ybs[1866]=['38 Ori',1.4635094,0.0659744,5.36];
ybs[1867]=['',1.4624043,-0.0178429,6.22];
ybs[1868]=['',1.4623958,-0.025435,5.93];
ybs[1869]=['121 Tau',1.4694191,0.4197849,5.38];
ybs[1870]=['φ1 Ori',1.4660834,0.1658457,4.41];
ybs[1871]=['',1.4555112,-0.671936,5.48];
ybs[1872]=['',1.471647,0.4830098,6.27];
ybs[1873]=['λ Ori',1.467489,0.1736047,3.54];
ybs[1874]=['λ Ori',1.4675036,0.1736193,5.61];
ybs[1875]=['',1.4568726,-0.6130629,5.78];
ybs[1876]=['',1.4416285,-1.1154805,6.19];
ybs[1877]=['',1.4678713,0.1789417,5.6];
ybs[1878]=['',1.4764516,0.7015157,6.09];
ybs[1879]=['',1.6062118,1.4814145,6.11];
ybs[1880]=['',1.4663684,-0.1046572,5.67];
ybs[1881]=['',1.4664995,-0.1045314,4.78];
ybs[1882]=['',1.4603925,-0.5207275,6.53];
ybs[1883]=['',1.4740923,0.4529365,6.49];
ybs[1884]=['',1.4679437,-0.0782042,6.56];
ybs[1885]=['',1.467997,-0.0770165,6.24];
ybs[1886]=['42 Ori',1.4680329,-0.0842257,4.59];
ybs[1887]=['θ1 Ori',1.4674819,-0.0938046,6.73];
ybs[1888]=['θ1 Ori',1.4674966,-0.0937706,7.96];
ybs[1889]=['θ1 Ori',1.4675255,-0.0938482,5.13];
ybs[1890]=['θ1 Ori',1.4675837,-0.0938145,6.7];
ybs[1891]=['θ2 Ori',1.4679899,-0.0943098,5.08];
ybs[1892]=['',1.4686246,-0.0759609,6.38];
ybs[1893]=['ι Ori',1.4681972,-0.1029303,2.77];
ybs[1894]=['',1.4690069,-0.0565545,6.4];
ybs[1895]=['45 Ori',1.4692176,-0.0845336,5.26];
ybs[1896]=['',1.4769442,0.4701173,5.83];
ybs[1897]=['ε Ori',1.4717721,-0.0207664,1.7];
ybs[1898]=['',1.4799514,0.585914,6.33];
ybs[1899]=['122 Tau',1.4761618,0.2976121,5.54];
ybs[1900]=['',1.4717702,-0.0983659,6.54];
ybs[1901]=['φ2 Ori',1.4751738,0.1623556,4.09];
ybs[1902]=['',1.4759737,0.1928003,5.94];
ybs[1903]=['',1.4662885,-0.5771334,5.78];
ybs[1904]=['ζ Tau',1.4788668,0.3692038,3];
ybs[1905]=['',1.4732601,-0.105646,5.72];
ybs[1906]=['',1.4580708,-0.9579876,6.43];
ybs[1907]=['',1.4769793,0.1564419,6.12];
ybs[1908]=['26 Aur',1.4836074,0.5323827,5.4];
ybs[1909]=['',1.4704958,-0.5008322,6.26];
ybs[1910]=['',1.5035998,1.1467926,5.6];
ybs[1911]=['',1.45347,-1.1207362,5.34];
ybs[1912]=['',1.4770245,-0.1034431,6.05];
ybs[1913]=['',1.4754547,-0.2053188,6.11];
ybs[1914]=['',1.4799669,0.1318169,5.88];
ybs[1915]=['',1.4848286,0.4647534,6.37];
ybs[1916]=['β Dor',1.4565379,-1.0904124,3.76];
ybs[1917]=['',1.4789564,-0.083817,6.19];
ybs[1918]=['',1.4864768,0.510085,5.96];
ybs[1919]=['',1.4969908,0.9335835,6.23];
ybs[1920]=['ν1 Col',1.4753486,-0.4862445,6.16];
ybs[1921]=['',1.4688285,-0.8255678,6.11];
ybs[1922]=['125 Tau',1.4881998,0.4521653,5.18];
ybs[1923]=['',1.4867692,0.3800133,6.34];
ybs[1924]=['',1.4632578,-1.0272691,6.75];
ybs[1925]=['σ Ori',1.482775,-0.04519,3.81];
ybs[1926]=['',1.4829424,-0.0450885,6.65];
ybs[1927]=['',1.4821206,-0.1145462,5.96];
ybs[1928]=['ω Ori',1.484932,0.0721161,4.57];
ybs[1929]=['ν2 Col',1.4773536,-0.5005266,5.31];
ybs[1930]=['',1.462574,-1.067493,6.32];
ybs[1931]=['49 Ori',1.4832098,-0.1257041,4.8];
ybs[1932]=['',1.4922178,0.5474716,6.04];
ybs[1933]=['',1.4926969,0.5572929,6.11];
ybs[1934]=['',1.4861141,-0.0620345,6];
ybs[1935]=['24 Cam',1.5047031,0.9876827,6.05];
ybs[1936]=['',1.4858587,-0.1692314,6.5];
ybs[1937]=['23 Cam',1.5102531,1.0731056,6.15];
ybs[1938]=['',1.4844901,-0.3113472,6.38];
ybs[1939]=['',1.4954064,0.5148176,6.43];
ybs[1940]=['126 Tau',1.4946069,0.2887355,4.86];
ybs[1941]=['',1.4809924,-0.7102894,5.82];
ybs[1942]=['ζ Ori',1.4915763,-0.0337376,2.05];
ybs[1943]=['ζ Ori',1.4915836,-0.0337377,4.21];
ybs[1944]=['',1.4909479,-0.0491341,6.22];
ybs[1945]=['',1.4975703,0.4072811,6.59];
ybs[1946]=['',1.4919768,-0.0195333,4.95];
ybs[1947]=['γ Men',1.4443899,-1.3321442,5.19];
ybs[1948]=['',1.4982181,0.3956539,6.36];
ybs[1949]=['',1.4931209,0.0060624,5.93];
ybs[1950]=['α Col',1.4853961,-0.5945247,2.64];
ybs[1951]=['',1.4913009,-0.1815086,6.52];
ybs[1952]=['',1.4862639,-0.5693064,5.45];
ybs[1953]=['',1.4955268,-0.0503847,6.42];
ybs[1954]=['',1.4700699,-1.1614859,6.31];
ybs[1955]=['',1.5037395,0.405135,6.21];
ybs[1956]=['',1.4950936,-0.2917536,6.21];
ybs[1957]=['51 Ori',1.4991983,0.0258933,4.91];
ybs[1958]=['',1.4582928,-1.2867979,5.78];
ybs[1959]=['',1.4974466,-0.3058037,6.15];
ybs[1960]=['',1.4932918,-0.5827844,6.34];
ybs[1961]=['',1.5007289,-0.1184635,6.02];
ybs[1962]=['12 Lep',1.4972181,-0.3903355,5.87];
ybs[1963]=['26 Cam',1.5198351,0.9795159,5.94];
ybs[1964]=['',1.5020457,-0.0280046,6.31];
ybs[1965]=['ο Aur',1.5165522,0.8697554,5.47];
ybs[1966]=['',1.4966796,-0.5327876,6.19];
ybs[1967]=['',1.4967223,-0.6049088,5.29];
ybs[1968]=['',1.5155362,0.7071067,6.58];
ybs[1969]=['',1.5022993,-0.3237422,5.73];
ybs[1970]=['',1.5320272,1.0962976,6.13];
ybs[1971]=['',1.5138077,0.3613207,6.95];
ybs[1972]=['',1.510439,0.070085,6.09];
ybs[1973]=['',1.5218639,0.7423396,6.29];
ybs[1974]=['',1.5070928,-0.3511346,6.34];
ybs[1975]=['',1.5019122,-0.6876339,6.25];
ybs[1976]=['',1.5068584,-0.3911943,6.15];
ybs[1977]=['γ Lep',1.5069518,-0.3916599,3.6];
ybs[1978]=['',1.5022716,-0.7997919,6.39];
ybs[1979]=['129 Tau',1.5184181,0.2762697,6];
ybs[1980]=['',1.514566,-0.0743742,6.34];
ybs[1981]=['',1.5186567,0.1663084,5.79];
ybs[1982]=['',1.5171058,0.0205036,5.95];
ybs[1983]=['131 Tau',1.5203804,0.25298,5.72];
ybs[1984]=['130 Tau',1.5214536,0.3095411,5.49];
ybs[1985]=['ι Men',1.4586596,-1.3754552,6.05];
ybs[1986]=['29 Cam',1.5376487,0.9934996,6.54];
ybs[1987]=['133 Tau',1.5225174,0.2427023,5.29];
ybs[1988]=['',1.5500283,1.1951051,6.2];
ybs[1989]=['τ Aur',1.5300727,0.6839312,4.52];
ybs[1990]=['μ Col',1.5131921,-0.563729,5.17];
ybs[1991]=['',1.5256692,0.3643406,6.07];
ybs[1992]=['ζ Lep',1.518133,-0.2585773,3.55];
ybs[1993]=['52 Ori',1.5234952,0.1127505,5.27];
ybs[1994]=['',1.5188265,-0.2832898,6.17];
ybs[1995]=['',1.5204351,-0.1837267,6.03];
ybs[1996]=['132 Tau',1.5286348,0.4288779,4.86];
ybs[1997]=['',1.5386989,0.8991763,6.29];
ybs[1998]=['κ Ori',1.5218274,-0.1686616,2.06];
ybs[1999]=['',1.518084,-0.4997339,6.22];
ybs[2000]=['30 Cam',1.5454338,1.0291804,6.14];
ybs[2001]=['',1.5256332,-0.0713672,5.97];
ybs[2002]=['',1.5143055,-0.8131543,5.31];
ybs[2003]=['',1.5187365,-0.6225292,6.32];
ybs[2004]=['134 Tau',1.530469,0.2208932,4.91];
ybs[2005]=['υ Aur',1.5381,0.6511798,4.74];
ybs[2006]=['ν Aur',1.5401723,0.683343,3.97];
ybs[2007]=['',1.5373049,0.4882059,5.56];
ybs[2008]=['',1.5325236,0.1723687,5.8];
ybs[2009]=['δ Dor',1.5045421,-1.1471634,4.35];
ybs[2010]=['135 Tau',1.5346053,0.2497602,5.52];
ybs[2011]=['',1.521315,-0.7094135,6.61];
ybs[2012]=['',1.5395146,0.5607587,6.25];
ybs[2013]=['',1.5330854,0.0772857,5.97];
ybs[2014]=['β Pic',1.517543,-0.8911631,3.85];
ybs[2015]=['',1.5297231,-0.2526965,5.49];
ybs[2016]=['π Men',1.4635925,-1.4042362,5.65];
ybs[2017]=['',1.5169175,-0.948661,6.18];
ybs[2018]=['',1.5342268,0.0354147,5.98];
ybs[2019]=['',1.5452965,0.6907645,6.45];
ybs[2020]=['',1.5306138,-0.4008432,5.87];
ybs[2021]=['31 Cam',1.5572296,1.0452857,5.2];
ybs[2022]=['',1.5450182,0.5920323,5.98];
ybs[2023]=['ξ Aur',1.5561817,0.9723084,4.99];
ybs[2024]=['',1.5431541,0.3468265,6.06];
ybs[2025]=['55 Ori',1.5376589,-0.1311409,5.35];
ybs[2026]=['',1.5280099,-0.7831291,6.38];
ybs[2027]=['137 Tau',1.5428468,0.2474062,5.59];
ybs[2028]=['136 Tau',1.5475789,0.4819789,4.58];
ybs[2029]=['δ Lep',1.5369418,-0.3643353,3.81];
ybs[2030]=['',1.5375314,-0.4000625,6.17];
ybs[2031]=['56 Ori',1.5426854,0.0324397,4.78];
ybs[2032]=['',1.5472227,0.3543424,6.71];
ybs[2033]=['',1.5409183,-0.1577349,5.97];
ybs[2034]=['β Col',1.5346562,-0.6241963,3.12];
ybs[2035]=['',1.5697731,1.1536069,6.25];
ybs[2036]=['γ Pic',1.5281087,-0.9802024,4.51];
ybs[2037]=['',1.5394827,-0.513906,6.45];
ybs[2038]=['',1.5313044,-0.9208869,6.35];
ybs[2039]=['',1.5618533,0.9041749,6.49];
ybs[2040]=['',1.5550052,0.5533379,5.9];
ybs[2041]=['χ1 Ori',1.5518614,0.3539303,4.41];
ybs[2042]=['',1.5507894,0.1848242,6.12];
ybs[2043]=['',1.5331578,-0.9093909,5.17];
ybs[2044]=['',1.552201,0.2053387,6.59];
ybs[2045]=['',1.5506839,0.0563389,6.31];
ybs[2046]=['57 Ori',1.5542977,0.344738,5.92];
ybs[2047]=['',1.5415024,-0.6567223,5.63];
ybs[2048]=['',1.5652765,0.8557443,6.47];
ybs[2049]=['',1.5425043,-0.6723402,6.7];
ybs[2050]=['λ Col',1.5441649,-0.5898863,4.87];
ybs[2051]=['',1.5526668,0.0169435,6];
ybs[2052]=['',1.5517995,-0.0708836,6.57];
ybs[2053]=['',1.4230825,-1.4794883,6.2];
ybs[2054]=['',1.5484902,-0.3427028,6.69];
ybs[2055]=['α Ori',1.5548084,0.1293143,0.5];
ybs[2056]=['λ Men',1.5156821,-1.2687801,6.53];
ybs[2057]=['',1.5581405,0.3521525,5.4];
ybs[2058]=['ε Dor',1.5266131,-1.1675525,5.11];
ybs[2059]=['',1.5521599,-0.2054496,5.66];
ybs[2060]=['',1.5617643,0.5051622,6.32];
ybs[2061]=['',1.5589236,0.2430723,6.6];
ybs[2062]=['',1.5492852,-0.508676,6.36];
ybs[2063]=['',1.5447846,-0.7490623,6.55];
ybs[2064]=['',1.5558227,-0.0805351,5.87];
ybs[2065]=['',1.5561872,-0.0835368,6.28];
ybs[2066]=['',1.5389399,-0.9974941,5.94];
ybs[2067]=['',1.5337208,-1.1175193,6.36];
ybs[2068]=['',1.5631718,0.4232594,6.02];
ybs[2069]=['',1.5605373,0.1660028,5.99];
ybs[2070]=['',1.5621767,0.2011048,5.87];
ybs[2071]=['δ Aur',1.5764447,0.9474434,3.72];
ybs[2072]=['',1.6062211,1.3191612,6.4];
ybs[2073]=['',1.5775838,0.9655247,6.44];
ybs[2074]=['',1.5776729,0.9520223,6.14];
ybs[2075]=['',1.5753094,0.8713441,5.89];
ybs[2076]=['',1.5514925,-0.6973561,5.57];
ybs[2077]=['',1.5476816,-0.878931,6.52];
ybs[2078]=['139 Tau',1.5678728,0.452993,4.82];
ybs[2079]=['η Lep',1.559391,-0.2472461,3.71];
ybs[2080]=['',1.5583113,-0.3986078,5.96];
ybs[2081]=['ξ Col',1.5543813,-0.6478431,4.97];
ybs[2082]=['β Aur',1.5756316,0.784479,1.9];
ybs[2083]=['',1.5499571,-0.8661076,6.1];
ybs[2084]=['',1.559764,-0.4051608,6.36];
ybs[2085]=['π Aur',1.5774793,0.8017442,4.26];
ybs[2086]=['σ Col',1.5583968,-0.5476984,5.5];
ybs[2087]=['',1.5665007,0.0213846,6.22];
ybs[2088]=['',1.5503399,-0.918614,5.29];
ybs[2089]=['θ Aur',1.5759764,0.6494765,2.62];
ybs[2090]=['',1.5790456,0.7782662,6.22];
ybs[2091]=['',1.5676923,-0.01734,6.22];
ybs[2092]=['',1.5604104,-0.558063,6.44];
ybs[2093]=['',1.5712215,0.223557,5.7];
ybs[2094]=['59 Ori',1.5687192,0.0320701,5.9];
ybs[2095]=['36 Aur',1.582178,0.8360303,5.73];
ybs[2096]=['',1.5457536,-1.1010752,4.65];
ybs[2097]=['60 Ori',1.570505,0.0096583,5.22];
ybs[2098]=['',1.5459152,-1.1253747,6.63];
ybs[2099]=['',1.5855017,0.8544803,5.96];
ybs[2100]=['γ Col',1.5633779,-0.6157913,4.36];
ybs[2101]=['1 Mon',1.5709686,-0.1637464,6.12];
ybs[2102]=['2 Mon',1.571202,-0.1668206,5.03];
ybs[2103]=['',1.57393,-0.0252119,6.63];
ybs[2104]=['',1.5819664,0.5416411,5.98];
ybs[2105]=['',1.5810876,0.4812108,6.05];
ybs[2106]=['',1.5903536,0.8709831,6.05];
ybs[2107]=['',1.5757395,-0.0536597,4.53];
ybs[2108]=['',1.5607463,-0.9324384,6.45];
ybs[2109]=['',1.5902027,0.7570661,6.42];
ybs[2110]=['',1.583854,0.3909469,6.37];
ybs[2111]=['',1.5675923,-0.7685364,5.81];
ybs[2112]=['',1.5764131,-0.2251498,6.22];
ybs[2113]=['38 Aur',1.5919451,0.7489128,6.1];
ybs[2114]=['η Col',1.5699455,-0.7472628,3.96];
ybs[2115]=['',1.6015546,1.0365489,6.34];
ybs[2116]=['',1.5896757,0.5695697,6.24];
ybs[2117]=['',1.5978063,0.9000759,6.45];
ybs[2118]=['μ Ori',1.5863588,0.1683535,4.12];
ybs[2119]=['κ Men',1.5220786,-1.3850223,5.47];
ybs[2120]=['',1.608814,1.1074044,6.39];
ybs[2121]=['',1.5856364,0.0295477,6.59];
ybs[2122]=['3 Mon',1.583242,-0.1849924,4.95];
ybs[2123]=['',1.5799276,-0.443639,6.05];
ybs[2124]=['64 Ori',1.5914271,0.3436277,5.14];
ybs[2125]=['',1.579756,-0.5918854,5.55];
ybs[2126]=['39 Aur',1.5996144,0.7501185,5.87];
ybs[2127]=['',1.5909224,0.2038324,6.08];
ybs[2128]=['1 Gem',1.5944792,0.4059783,4.16];
ybs[2129]=['χ2 Ori',1.5934746,0.3514386,4.63];
ybs[2130]=['',1.5862489,-0.2530521,6.2];
ybs[2131]=['',1.5992396,0.662547,6.34];
ybs[2132]=['',1.5765826,-0.8939044,5.67];
ybs[2133]=['',1.6012809,0.5863589,6.23];
ybs[2134]=['',1.5887899,-0.4587837,5.04];
ybs[2135]=['',1.6038857,0.6175659,6.12];
ybs[2136]=['',1.5937937,-0.1171403,5.21];
ybs[2137]=['40 Aur',1.6059966,0.6715844,5.36];
ybs[2138]=['63 Ori',1.5974882,0.0945463,5.67];
ybs[2139]=['66 Ori',1.5974565,0.072531,5.63];
ybs[2140]=['',1.6046025,0.5150259,6.08];
ybs[2141]=['',1.609974,0.7304181,6.12];
ybs[2142]=['17 Lep',1.5967375,-0.2877576,4.93];
ybs[2143]=['',1.59322,-0.5615593,5.65];
ybs[2144]=['',1.5990071,-0.1788245,5.87];
ybs[2145]=['',1.5813676,-1.0489105,6.45];
ybs[2146]=['37 Cam',1.6226359,1.0285247,5.36];
ybs[2147]=['',1.6140016,0.7164713,6.36];
ybs[2148]=['',1.6044455,-0.0732626,5.38];
ybs[2149]=['θ Lep',1.6019036,-0.2607303,4.67];
ybs[2150]=['',1.599806,-0.4223487,6.95];
ybs[2151]=['',1.5930316,-0.7860817,6.35];
ybs[2152]=['',1.5938866,-0.7868204,5.93];
ybs[2153]=['ν Ori',1.6091924,0.2576814,4.42];
ybs[2154]=['',1.5979013,-0.6198826,5.8];
ybs[2155]=['',1.6051461,-0.1950835,6.66];
ybs[2156]=['',1.594111,-0.8458084,6.58];
ybs[2157]=['',1.6032203,-0.4034189,5.47];
ybs[2158]=['',1.6009903,-0.519445,5.81];
ybs[2159]=['36 Cam',1.6363018,1.1468744,5.32];
ybs[2160]=['',1.6051374,-0.3807724,5.78];
ybs[2161]=['',1.6142617,0.1512346,6.55];
ybs[2162]=['19 Lep',1.6084472,-0.3345813,5.31];
ybs[2163]=['',1.6180816,0.3871958,5.93];
ybs[2164]=['',1.6049917,-0.5989243,5.83];
ybs[2165]=['π1 Col',1.6028659,-0.7383139,6.12];
ybs[2166]=['',1.6296788,0.9187523,6.3];
ybs[2167]=['3 Gem',1.6189647,0.4033092,5.75];
ybs[2168]=['',1.6148126,0.0435367,5.73];
ybs[2169]=['41 Aur',1.6286306,0.8500899,6.82];
ybs[2170]=['41 Aur',1.6286377,0.850056,6.09];
ybs[2171]=['θ Col',1.6068577,-0.6502605,5.02];
ybs[2172]=['',1.604191,-0.7870601,6.51];
ybs[2173]=['',1.6172981,-0.09977,6.17];
ybs[2174]=['',1.61386,-0.3915195,5.5];
ybs[2175]=['π2 Col',1.6081045,-0.7357991,5.5];
ybs[2176]=['',1.6156672,-0.3164497,6.35];
ybs[2177]=['',1.6168256,-0.2546334,5.56];
ybs[2178]=['',1.6244149,0.3163123,6.33];
ybs[2179]=['5 Gem',1.6268969,0.4261032,5.8];
ybs[2180]=['',1.6174814,-0.3975776,5.71];
ybs[2181]=['',1.6109859,-0.7742412,6.27];
ybs[2182]=['',1.638212,0.8929956,6.04];
ybs[2183]=['',1.6307645,0.5704875,5.78];
ybs[2184]=['',1.628176,0.3815654,6.56];
ybs[2185]=['',1.626136,0.2379287,6.04];
ybs[2186]=['',1.617257,-0.4661106,6.27];
ybs[2187]=['68 Ori',1.6288092,0.3452951,5.75];
ybs[2188]=['η1 Dor',1.5977684,-1.1526668,5.71];
ybs[2189]=['',1.6234483,-0.1179876,6.15];
ybs[2190]=['',1.6024418,-1.0848697,5.05];
ybs[2191]=['6 Gem',1.6302346,0.3997076,6.39];
ybs[2192]=['69 Ori',1.6288074,0.2814158,4.95];
ybs[2193]=['ξ Ori',1.6282305,0.2478775,4.48];
ybs[2194]=['',1.6206834,-0.4740299,5.72];
ybs[2195]=['40 Cam',1.6476149,1.0470314,5.35];
ybs[2196]=['',1.6266085,-0.081541,6.18];
ybs[2197]=['',1.6182133,-0.7043993,5.58];
ybs[2198]=['',1.6141193,-0.8651217,6.49];
ybs[2199]=['',1.6271283,-0.1144368,5.05];
ybs[2200]=['',1.6235431,-0.4623082,6.09];
ybs[2201]=['',1.6354613,0.3259031,6.58];
ybs[2202]=['',1.633635,0.1853591,6.45];
ybs[2203]=['',1.6633502,1.2096751,4.8];
ybs[2204]=['',1.6311018,-0.0438318,6.62];
ybs[2205]=['',1.6200317,-0.790419,6.31];
ybs[2206]=['δ Pic',1.6175622,-0.959479,4.81];
ybs[2207]=['',1.630663,-0.3101445,6.52];
ybs[2208]=['',1.6394442,0.3123878,5.88];
ybs[2209]=['1 Lyn',1.6576023,1.073472,4.98];
ybs[2210]=['η Gem',1.6413795,0.392674,3.28];
ybs[2211]=['',1.6454124,0.6307632,6.92];
ybs[2212]=['',1.6361399,-0.0654311,5.83];
ybs[2213]=['κ Aur',1.6438838,0.5146919,4.35];
ybs[2214]=['71 Ori',1.6411164,0.3342011,5.2];
ybs[2215]=['ν Dor',1.6083304,-1.2016218,5.06];
ybs[2216]=['',1.6421861,0.241604,5.91];
ybs[2217]=['72 Ori',1.6434826,0.2816034,5.3];
ybs[2218]=['',1.6391931,-0.0798704,5.83];
ybs[2219]=['',1.6346849,-0.4165988,6.39];
ybs[2220]=['',1.6335743,-0.5131862,6.54];
ybs[2221]=['γ Mon',1.6401922,-0.1096547,3.98];
ybs[2222]=['42 Aur',1.6544926,0.8100872,6.52];
ybs[2223]=['73 Ori',1.6447905,0.2189093,5.33];
ybs[2224]=['8 Gem',1.6477254,0.4182008,6.08];
ybs[2225]=['',1.6441917,0.1057256,6.07];
ybs[2226]=['',1.6446281,0.0746142,6.64];
ybs[2227]=['',1.6435298,-0.0090868,5.65];
ybs[2228]=['',1.6430347,-0.0859338,5.99];
ybs[2229]=['',1.6477924,0.2997168,6.39];
ybs[2230]=['',1.6450237,0.0202559,6.37];
ybs[2231]=['',1.6426134,-0.1578405,6.1];
ybs[2232]=['2 Lyn',1.664691,1.029746,4.48];
ybs[2233]=['43 Aur',1.657556,0.8089706,6.38];
ybs[2234]=['9 Gem',1.6506023,0.4141951,6.25];
ybs[2235]=['74 Ori',1.647805,0.2140355,5.04];
ybs[2236]=['',1.6408951,-0.3539591,5.91];
ybs[2237]=['',1.6416303,-0.3226223,5.99];
ybs[2238]=['',1.6438106,-0.239578,5.01];
ybs[2239]=['η2 Dor',1.6200936,-1.1448545,5.01];
ybs[2240]=['',1.6469984,0.0187004,6.63];
ybs[2241]=['75 Ori',1.6506258,0.1733682,5.39];
ybs[2242]=['',1.6499221,0.1229391,6.57];
ybs[2243]=['',1.6453556,-0.2901862,5.92];
ybs[2244]=['',1.6527239,0.2451988,6.59];
ybs[2245]=['',1.6511374,0.0888541,5.71];
ybs[2246]=['',1.644021,-0.5200536,6.67];
ybs[2247]=['',1.6550852,0.2508565,6.16];
ybs[2248]=['',1.6491632,-0.3966109,6.07];
ybs[2249]=['6 Mon',1.651946,-0.1873561,6.75];
ybs[2250]=['κ Col',1.646358,-0.6134726,4.37];
ybs[2251]=['4 Lyn',1.6753692,1.0360314,5.94];
ybs[2252]=['',1.6592876,0.3021994,6.32];
ybs[2253]=['',1.6574208,0.1577284,6.24];
ybs[2254]=['',1.6521838,-0.2936572,5.14];
ybs[2255]=['α Men',1.6126263,-1.3047772,5.09];
ybs[2256]=['',1.6463031,-0.6854482,6];
ybs[2257]=['',1.6482557,-0.6588019,5.53];
ybs[2258]=['45 Aur',1.6733907,0.9327113,5.36];
ybs[2259]=['',1.6488874,-0.6503481,5.87];
ybs[2260]=['',1.6543838,-0.3486591,5.52];
ybs[2261]=['',1.6574795,-0.1640675,5.36];
ybs[2262]=['',1.6571341,-0.2624066,6.06];
ybs[2263]=['',1.66372,0.2555218,5.69];
ybs[2264]=['',1.6587676,-0.1500397,6.22];
ybs[2265]=['',1.6576386,-0.3654065,5.81];
ybs[2266]=['',1.6692923,0.5153905,6.43];
ybs[2267]=['7 Mon',1.661334,-0.1367223,5.27];
ybs[2268]=['',1.643265,-1.0336224,6.43];
ybs[2269]=['',1.6627346,-0.0515772,4.9];
ybs[2270]=['',1.6671069,0.2049923,6.54];
ybs[2271]=['',1.6697836,0.3098327,6.35];
ybs[2272]=['',1.6508032,-0.9205302,6.41];
ybs[2273]=['',1.6600641,-0.6005178,5.78];
ybs[2274]=['',1.6691926,0.0393946,6.31];
ybs[2275]=['',1.6550342,-0.8791067,7.04];
ybs[2276]=['ζ CMa',1.6630335,-0.5248929,3.02];
ybs[2277]=['',1.6351774,-1.2515856,6.64];
ybs[2278]=['',1.6685967,-0.2056828,5.64];
ybs[2279]=['',1.704699,1.2308098,5.97];
ybs[2280]=['μ Gem',1.6766473,0.3927218,2.88];
ybs[2281]=['',1.6747151,0.2191767,6];
ybs[2282]=['',1.6640996,-0.5961144,5.53];
ybs[2283]=['ψ1 Aur',1.6866558,0.8600049,4.91];
ybs[2284]=['',1.660942,-0.8508784,6.6];
ybs[2285]=['',1.6940347,0.98211,5.64];
ybs[2286]=['',1.6774427,0.0654848,6.4];
ybs[2287]=['5 Lyn',1.6959791,1.0193205,5.21];
ybs[2288]=['β CMa',1.6739781,-0.3135993,1.98];
ybs[2289]=['',1.677441,-0.0820251,6.67];
ybs[2290]=['δ Col',1.6707189,-0.5837799,3.85];
ybs[2291]=['',1.6853336,0.5182563,6.71];
ybs[2292]=['ε Mon',1.6794799,0.0799377,4.44];
ybs[2293]=['',1.6795091,0.0799861,6.72];
ybs[2294]=['',1.6808158,0.1548437,6.26];
ybs[2295]=['',1.678217,-0.172561,6.19];
ybs[2296]=['',1.6847685,0.2800145,6.33];
ybs[2297]=['',1.6787464,-0.2632709,6.24];
ybs[2298]=['',1.6879748,0.4068937,6.06];
ybs[2299]=['',1.6806492,-0.2014658,5.22];
ybs[2300]=['',1.6786748,-0.3455389,6.6];
ybs[2301]=['',1.6757199,-0.5550553,6.34];
ybs[2302]=['',1.6872761,0.2567092,6.24];
ybs[2303]=['',1.6813367,-0.2264643,6.12];
ybs[2304]=['',1.685898,0.1234364,5.98];
ybs[2305]=['',1.6790314,-0.4466381,5.63];
ybs[2306]=['',1.6860797,0.025964,6.66];
ybs[2307]=['',1.6858456,-0.0167477,5.87];
ybs[2308]=['',1.6993841,0.8271176,6.56];
ybs[2309]=['',1.6881585,0.0394183,6.51];
ybs[2310]=['',1.6788477,-0.6408935,5.62];
ybs[2311]=['',1.6879713,-0.0681182,6.35];
ybs[2312]=['',1.6823757,-0.5025346,6.39];
ybs[2313]=['',1.6973194,0.5680752,6.43];
ybs[2314]=['ν Pic',1.6725299,-0.9840525,5.61];
ybs[2315]=['',1.6886756,-0.1380202,6.4];
ybs[2316]=['',1.6760251,-0.9109494,5.98];
ybs[2317]=['',1.6817998,-0.7033196,6.31];
ybs[2318]=['',1.6918757,-0.0265534,5.87];
ybs[2319]=['',1.6913926,-0.0804832,6.15];
ybs[2320]=['α Car',1.677386,-0.9199359,-0.72];
ybs[2321]=['',1.6933498,0.0144249,6.71];
ybs[2322]=['',1.6920422,-0.1313465,6.27];
ybs[2323]=['',1.6854046,-0.6122157,6.25];
ybs[2324]=['16 Gem',1.6983047,0.3574647,6.22];
ybs[2325]=['6 Lyn',1.7132678,1.0148438,5.88];
ybs[2326]=['48 Aur',1.7014771,0.5319382,5.55];
ybs[2327]=['',1.6950028,0.0505014,5.55];
ybs[2328]=['',1.6944283,0.0049688,5.2];
ybs[2329]=['',1.6945383,-0.0050719,5.55];
ybs[2330]=['',1.6759743,-1.0220015,6.48];
ybs[2331]=['',1.6718519,-1.1116891,6.27];
ybs[2332]=['47 Aur',1.7089174,0.8145363,5.9];
ybs[2333]=['',1.7029472,0.4704073,6.47];
ybs[2334]=['',1.7004169,0.283148,6.23];
ybs[2335]=['',1.6810922,-0.9218732,6.51];
ybs[2336]=['',1.6995211,0.179574,6.15];
ybs[2337]=['ν Gem',1.7027432,0.3525007,4.15];
ybs[2338]=['10 Mon',1.6974299,-0.0833755,5.06];
ybs[2339]=['',1.6776804,-1.0523256,5.8];
ybs[2340]=['',1.7625375,1.3888907,6.54];
ybs[2341]=['',1.6990681,0.0331124,6.48];
ybs[2342]=['',1.685538,-0.8410878,5.76];
ybs[2343]=['',1.6932231,-0.4515352,6.07];
ybs[2344]=['',1.7845362,1.4327588,6.65];
ybs[2345]=['',1.7025439,0.1920566,6.59];
ybs[2346]=['π1 Dor',1.6686665,-1.2216592,5.56];
ybs[2347]=['',1.6923448,-0.6616522,6.48];
ybs[2348]=['',1.6780824,-1.107266,6.46];
ybs[2349]=['',1.7033198,0.0459124,6.16];
ybs[2350]=['β Mon',1.701085,-0.1230118,4.6];
ybs[2351]=['β Mon',1.7011212,-0.123041,5.4];
ybs[2352]=['β Mon',1.7011212,-0.123041,5.6];
ybs[2353]=['',1.6998376,-0.3051056,5.77];
ybs[2354]=['',1.6801407,-1.1142323,6.27];
ybs[2355]=['λ CMa',1.697187,-0.568888,4.48];
ybs[2356]=['',1.7072328,0.1573099,6.57];
ybs[2357]=['',1.761882,1.3609022,5.73];
ybs[2358]=['',1.6993141,-0.5652466,5.74];
ybs[2359]=['',1.7480726,1.2858739,6.24];
ybs[2360]=['',1.7122317,0.2953456,6.2];
ybs[2361]=['',1.7069575,-0.1762324,5.93];
ybs[2362]=['',1.6990469,-0.7171484,6.32];
ybs[2363]=['',1.6903841,-1.0125775,5.82];
ybs[2364]=['',1.7119696,0.1960756,6.14];
ybs[2365]=['19 Gem',1.7141762,0.2772726,6.4];
ybs[2366]=['',1.7185231,0.5661406,5.87];
ybs[2367]=['',1.7085438,-0.2297593,6.16];
ybs[2368]=['',1.7141496,0.2055201,6.65];
ybs[2369]=['',1.714802,0.2011942,5.23];
ybs[2370]=['7 Lyn',1.7293286,0.9657713,6.45];
ybs[2371]=['π2 Dor',1.6811554,-1.2165557,5.38];
ybs[2372]=['',1.717352,0.2034432,6.03];
ybs[2373]=['',1.7120852,-0.2165651,5.15];
ybs[2374]=['',1.708781,-0.4849517,5.93];
ybs[2375]=['',1.714214,-0.1426788,5.43];
ybs[2376]=['12 Mon',1.7168023,0.0844515,5.84];
ybs[2377]=['',1.7240415,0.5760679,6.42];
ybs[2378]=['',1.7032215,-0.8771123,5.27];
ybs[2379]=['13 Mon',1.7194379,0.1276819,4.5];
ybs[2380]=['',1.7166977,-0.1027304,5.6];
ybs[2381]=['ξ1 CMa',1.7137019,-0.4090205,4.33];
ybs[2382]=['',1.7103494,-0.6156756,5.84];
ybs[2383]=['',1.7010403,-0.9925378,5.22];
ybs[2384]=['',1.7090586,-0.7144105,6.2];
ybs[2385]=['',1.7227401,0.2467455,5.53];
ybs[2386]=['',1.718233,-0.1951926,6.24];
ybs[2387]=['',1.7118566,-0.6450149,6.34];
ybs[2388]=['8 Lyn',1.7438465,1.0726969,5.94];
ybs[2389]=['',1.722306,-0.0216083,5.1];
ybs[2390]=['',1.7587049,1.2518754,5.92];
ybs[2391]=['',1.7167731,-0.5593389,5.69];
ybs[2392]=['49 Aur',1.7302985,0.4887543,5.27];
ybs[2393]=['',1.7151812,-0.6582283,5.24];
ybs[2394]=['',1.7095654,-0.904823,5.6];
ybs[2395]=['',1.788424,1.3882315,5.45];
ybs[2396]=['11 Lyn',1.7430003,0.9920008,5.85];
ybs[2397]=['',1.7207381,-0.3654987,6.4];
ybs[2398]=['14 Mon',1.7276058,0.1318441,6.45];
ybs[2399]=['',1.7367126,0.6706583,5.29];
ybs[2400]=['',1.7299637,0.1740035,5.88];
ybs[2401]=['',1.7186912,-0.6744429,6.44];
ybs[2402]=['',1.7021615,-1.1446567,6.29];
ybs[2403]=['',1.7295014,0.0152082,5.8];
ybs[2404]=['',1.7077641,-1.0802892,6.15];
ybs[2405]=['',1.7216795,-0.6326822,5.42];
ybs[2406]=['μ Pic',1.7117038,-1.0257457,5.7];
ybs[2407]=['',1.732846,0.0781643,6.55];
ybs[2408]=['ξ2 CMa',1.7276929,-0.4011324,4.54];
ybs[2409]=['',1.7251976,-0.5713264,5.62];
ybs[2410]=['',1.7188146,-0.9136173,6.19];
ybs[2411]=['',1.7399678,0.4288452,6.44];
ybs[2412]=['',1.7350627,-0.0912879,5.52];
ybs[2413]=['51 Aur',1.7459727,0.6871423,5.69];
ybs[2414]=['ψ3 Aur',1.7467089,0.6960711,5.2];
ybs[2415]=['γ Gem',1.7407584,0.2858717,1.93];
ybs[2416]=['',1.7390211,0.1067362,6.06];
ybs[2417]=['ν1 CMa',1.7336453,-0.3260129,5.7];
ybs[2418]=['',1.72853,-0.6422569,5.59];
ybs[2419]=['53 Aur',1.7442233,0.5055147,5.79];
ybs[2420]=['',1.7401254,0.1890749,6.38];
ybs[2421]=['ψ2 Aur',1.7491048,0.7412073,4.79];
ybs[2422]=['',1.7355917,-0.2328307,5.97];
ybs[2423]=['ν2 CMa',1.7349451,-0.3364149,3.95];
ybs[2424]=['',1.740075,0.0468497,6.17];
ybs[2425]=['',1.7307434,-0.6301992,6.35];
ybs[2426]=['',1.7410581,0.086171,6.15];
ybs[2427]=['',1.7347991,-0.3950432,6.35];
ybs[2428]=['',1.7519777,0.7678178,6.41];
ybs[2429]=['',1.7254807,-0.9249176,4.39];
ybs[2430]=['',1.7469917,0.3841502,6.04];
ybs[2431]=['',1.7395462,-0.2269774,6.12];
ybs[2432]=['54 Aur',1.7492865,0.4929186,6.03];
ybs[2433]=['',1.7489988,0.4289866,6.38];
ybs[2434]=['',1.7428096,-0.0447471,6.14];
ybs[2435]=['',1.7451793,0.0816829,6.57];
ybs[2436]=['',1.744239,0.0278073,6.21];
ybs[2437]=['ν3 CMa',1.7402515,-0.3186525,4.43];
ybs[2438]=['',1.735566,-0.6661243,6.04];
ybs[2439]=['',1.7345864,-0.7256431,6.34];
ybs[2440]=['',1.7364978,-0.6459482,5.71];
ybs[2441]=['',1.7392041,-0.5647812,5.27];
ybs[2442]=['',1.7433602,-0.2948545,6.03];
ybs[2443]=['',1.7497094,0.2262308,5.97];
ybs[2444]=['',1.7464695,-0.2472521,4.82];
ybs[2445]=['ν Pup',1.7384356,-0.75426,3.17];
ybs[2446]=['',1.7587299,0.6267469,6.46];
ybs[2447]=['25 Gem',1.7571285,0.4917388,6.42];
ybs[2448]=['',1.752679,0.1108339,6.51];
ybs[2449]=['',1.7475104,-0.4139287,6.05];
ybs[2450]=['15 Mon',1.7547641,0.1723331,4.66];
ybs[2451]=['',1.75669,0.2858098,6.28];
ybs[2452]=['',1.7561505,0.1916647,6.11];
ybs[2453]=['ψ4 Aur',1.7656073,0.7767008,5.02];
ybs[2454]=['',1.7476685,-0.5321706,5.71];
ybs[2455]=['',1.7549108,0.0082667,5.79];
ybs[2456]=['',1.7418397,-0.8419556,4.93];
ybs[2457]=['',1.7711814,0.9297894,6.27];
ybs[2458]=['',1.7657839,0.6479382,6.19];
ybs[2459]=['',1.7482951,-0.6663639,6.58];
ybs[2460]=['26 Gem',1.7612835,0.3075783,5.21];
ybs[2461]=['',1.759041,0.1103553,6.37];
ybs[2462]=['',1.737653,-1.0743003,6.18];
ybs[2463]=['',1.7582704,-0.160383,5.19];
ybs[2464]=['12 Lyn',1.7808046,1.0370257,4.87];
ybs[2465]=['',1.770001,0.6298265,6.31];
ybs[2466]=['ε Gem',1.7682535,0.4382167,2.98];
ybs[2467]=['',1.7638083,0.052546,6.19];
ybs[2468]=['',1.7537798,-0.7046125,6.12];
ybs[2469]=['',1.7514646,-0.8324537,6.65];
ybs[2470]=['13 Lyn',1.7830833,0.9973582,5.35];
ybs[2471]=['30 Gem',1.768021,0.2304644,4.49];
ybs[2472]=['',1.7661822,0.0682297,5.9];
ybs[2473]=['28 Gem',1.7720329,0.5052249,5.44];
ybs[2474]=['',1.7613464,-0.3922035,6.13];
ybs[2475]=['',1.7584288,-0.6705686,6.29];
ybs[2476]=['ψ5 Aur',1.7814833,0.7601409,5.25];
ybs[2477]=['ξ Gem',1.7736876,0.2246543,3.36];
ybs[2478]=['',1.7889736,0.9717818,6.33];
ybs[2479]=['',1.7889299,0.9717819,6.28];
ybs[2480]=['ψ6 Aur',1.7858917,0.8510978,5.22];
ybs[2481]=['',1.7632485,-0.6844491,6.3];
ybs[2482]=['32 Gem',1.7763559,0.2211243,6.46];
ybs[2483]=['42 Cam',1.8028511,1.1788824,5.14];
ybs[2484]=['α CMa',1.7719857,-0.2921644,-1.46];
ybs[2485]=['10 CMa',1.7684207,-0.54269,5.2];
ybs[2486]=['',1.7703052,-0.4776075,6.45];
ybs[2487]=['16 Mon',1.7789819,0.1494487,5.93];
ybs[2488]=['',1.7727592,-0.4099032,6.05];
ybs[2489]=['',1.7709173,-0.5342399,6.54];
ybs[2490]=['',1.7757254,-0.2586617,5.32];
ybs[2491]=['',1.7830607,0.3170991,6.2];
ybs[2492]=['',1.7723422,-0.5553178,5.92];
ybs[2493]=['',1.773002,-0.5405759,5.8];
ybs[2494]=['',1.7787885,-0.1768311,5.66];
ybs[2495]=['17 Mon',1.7824086,0.1398424,4.77];
ybs[2496]=['11 CMa',1.7795069,-0.2522068,5.29];
ybs[2497]=['',1.7480922,-1.2530897,6.51];
ybs[2498]=['18 Mon',1.7845187,0.0416631,4.47];
ybs[2499]=['',1.7748672,-0.6905234,6.62];
ybs[2500]=['',1.7830546,-0.157486,5.07];
ybs[2501]=['12 CMa',1.7800038,-0.3672204,6.08];
ybs[2502]=['',1.7756108,-0.6597294,6.21];
ybs[2503]=['43 Cam',1.815182,1.2018329,5.12];
ybs[2504]=['',1.7937142,0.5686381,5.71];
ybs[2505]=['',1.7712108,-0.911495,6.57];
ybs[2506]=['',1.7863782,-0.0234658,5.75];
ybs[2507]=['',1.773199,-0.9151449,5.8];
ybs[2508]=['ψ7 Aur',1.7989273,0.7287571,5.02];
ybs[2509]=['',1.7897115,0.0170384,6.15];
ybs[2510]=['',1.7806433,-0.6624304,5.26];
ybs[2511]=['33 Gem',1.7936207,0.2823356,5.85];
ybs[2512]=['14 Lyn',1.8106709,1.037086,5.33];
ybs[2513]=['',1.7905185,-0.0401034,5.74];
ybs[2514]=['',1.7886945,-0.2647726,5.39];
ybs[2515]=['',1.7776063,-0.8951843,5.4];
ybs[2516]=['',1.7764456,-0.9550276,6.46];
ybs[2517]=['35 Gem',1.7961104,0.2336455,5.65];
ybs[2518]=['',1.7790601,-0.969786,5.61];
ybs[2519]=['',1.8463919,1.3429559,4.55];
ybs[2520]=['',1.7916953,-0.4206564,6.33];
ybs[2521]=['36 Gem',1.8013455,0.3793313,5.27];
ybs[2522]=['',1.7973797,-0.0099037,5.77];
ybs[2523]=['',1.7591379,-1.276544,6.37];
ybs[2524]=['',1.8094324,0.782109,6.26];
ybs[2525]=['',1.8033854,0.4114511,5.65];
ybs[2526]=['',1.7965584,-0.1408069,6.29];
ybs[2527]=['',1.7947364,-0.2986298,5.79];
ybs[2528]=['',1.7659009,-1.2297093,6.11];
ybs[2529]=['',1.7931526,-0.4775233,7.04];
ybs[2530]=['κ CMa',1.7917796,-0.5678366,3.96];
ybs[2531]=['59 Aur',1.8086065,0.6779093,6.12];
ybs[2532]=['θ Gem',1.8073052,0.59225,3.6];
ybs[2533]=['60 Aur',1.8094465,0.6703833,6.3];
ybs[2534]=['',1.8085718,0.6241434,6.01];
ybs[2535]=['',1.8010995,0.0526152,6.38];
ybs[2536]=['',1.7954694,-0.4503734,6.33];
ybs[2537]=['',1.7942053,-0.5538352,5.7];
ybs[2538]=['',1.7915383,-0.7937067,6.55];
ybs[2539]=['ψ8 Aur',1.8126198,0.6715452,6.48];
ybs[2540]=['',1.7912194,-0.8140345,5.14];
ybs[2541]=['',1.7961874,-0.6002846,4.99];
ybs[2542]=['α Pic',1.7820383,-1.0815182,3.27];
ybs[2543]=['',1.8063886,0.145781,5.77];
ybs[2544]=['',1.8039729,-0.0932616,6.3];
ybs[2545]=['τ Pup',1.7910092,-0.8838476,2.93];
ybs[2546]=['',1.7903681,-0.9363373,4.4];
ybs[2547]=['',1.808891,0.1914358,6.24];
ybs[2548]=['',1.8187778,0.7993156,6.34];
ybs[2549]=['',1.8186065,0.7658685,6.13];
ybs[2550]=['',1.7997179,-0.6328083,5.96];
ybs[2551]=['ζ Men',1.7375892,-1.4108186,5.64];
ybs[2552]=['15 Lyn',1.8288016,1.0191401,4.35];
ybs[2553]=['',1.828454,1.0041454,6.05];
ybs[2554]=['',1.7902809,-1.0519999,6.11];
ybs[2555]=['',1.7982337,-0.8433316,6.42];
ybs[2556]=['38 Gem',1.8145126,0.2294968,4.65];
ybs[2557]=['',1.8075262,-0.3326704,5.64];
ybs[2558]=['',1.8077411,-0.3309352,6.14];
ybs[2559]=['',1.8058334,-0.4709798,6.4];
ybs[2560]=['ψ9 Aur',1.8244144,0.8071194,5.87];
ybs[2561]=['37 Gem',1.8179006,0.4423819,5.73];
ybs[2562]=['',1.8116334,-0.102639,6.41];
ybs[2563]=['15 CMa',1.8085042,-0.3534661,4.83];
ybs[2564]=['',1.8129724,-0.0201651,5.45];
ybs[2565]=['',1.8261692,0.8146353,5.86];
ybs[2566]=['θ CMa',1.8116195,-0.2106072,4.07];
ybs[2567]=['',1.8035198,-0.7423213,6.52];
ybs[2568]=['',1.8082204,-0.4985997,6.04];
ybs[2569]=['',1.8142228,-0.0311535,6.21];
ybs[2570]=['',1.8099544,-0.4287802,6.21];
ybs[2571]=['',1.8039596,-0.7680029,6.46];
ybs[2572]=['ο1 CMa',1.8108861,-0.4225813,3.87];
ybs[2573]=['',1.8491049,1.2352704,5.68];
ybs[2574]=['',1.8154001,-0.0494333,6.04];
ybs[2575]=['',1.8112679,-0.4181218,6.91];
ybs[2576]=['',1.8183994,0.144787,6.29];
ybs[2577]=['16 Lyn',1.829039,0.7865151,4.9];
ybs[2578]=['',1.8256978,0.5873258,5.89];
ybs[2579]=['',1.8030984,-0.9445275,6.57];
ybs[2580]=['17 CMa',1.8150285,-0.3566307,5.74];
ybs[2581]=['',1.8221817,0.1732574,5.92];
ybs[2582]=['π CMa',1.8175632,-0.3519525,4.68];
ybs[2583]=['',1.8113135,-0.739913,6.32];
ybs[2584]=['',1.8023662,-1.0361759,6.41];
ybs[2585]=['μ CMa',1.8199232,-0.245618,5];
ybs[2586]=['',1.8089166,-0.8838306,6.26];
ybs[2587]=['',1.8181322,-0.4009103,5.3];
ybs[2588]=['ι CMa',1.8199248,-0.2981622,4.37];
ybs[2589]=['',1.8266087,0.2073018,6.27];
ybs[2590]=['',1.8183158,-0.5553484,6.36];
ybs[2591]=['',1.8240297,-0.1432673,6.34];
ybs[2592]=['62 Aur',1.8348231,0.6635979,6];
ybs[2593]=['39 Gem',1.8331119,0.4546654,6.1];
ybs[2594]=['ι Vol',1.7941895,-1.2390077,5.4];
ybs[2595]=['',1.8245637,-0.3880417,6.61];
ybs[2596]=['',1.8218375,-0.6173439,6.29];
ybs[2597]=['40 Gem',1.836049,0.4517457,6.4];
ybs[2598]=['',1.831783,0.1324941,6.27];
ybs[2599]=['',1.8258512,-0.4304075,5.46];
ybs[2600]=['',1.8188052,-0.8508541,4.95];
ybs[2601]=['',2.0498193,1.5158076,5.07];
ybs[2602]=['',1.8329488,0.0623342,5.97];
ybs[2603]=['',1.8263422,-0.4811443,6.23];
ybs[2604]=['',1.8241494,-0.620243,6.23];
ybs[2605]=['',1.8347608,0.1271648,6.35];
ybs[2606]=['',1.828192,-0.4746418,6.37];
ybs[2607]=['41 Gem',1.8391332,0.2800812,5.68];
ybs[2608]=['',1.8303274,-0.4440882,5.59];
ybs[2609]=['',1.868704,1.2339026,6.5];
ybs[2610]=['ε CMa',1.830285,-0.506193,1.5];
ybs[2611]=['',1.829133,-0.5958911,5.06];
ybs[2612]=['',1.8442955,0.5651807,6.59];
ybs[2613]=['',1.8306481,-0.5415464,6.42];
ybs[2614]=['',1.8385228,-0.0942189,6.3];
ybs[2615]=['',1.8350989,-0.377591,6.26];
ybs[2616]=['',1.8388288,-0.1472776,5.96];
ybs[2617]=['',1.837252,-0.352385,6.31];
ybs[2618]=['',1.8296407,-0.7993353,6.22];
ybs[2619]=['',1.8399347,-0.1611746,6.49];
ybs[2620]=['',1.8379873,-0.3866047,6.53];
ybs[2621]=['',1.8449468,0.0835303,6.63];
ybs[2622]=['ω Gem',1.8488375,0.4220686,5.18];
ybs[2623]=['',1.848628,0.3093254,5.94];
ybs[2624]=['',1.8479467,0.2670994,5.74];
ybs[2625]=['',1.8459621,0.0964341,6.59];
ybs[2626]=['',1.82859,-0.9731931,6.27];
ybs[2627]=['',1.8491681,0.2904505,5.82];
ybs[2628]=['',1.8455663,-0.0240464,6.17];
ybs[2629]=['',1.8394074,-0.4977853,6.27];
ybs[2630]=['',1.8282671,-0.9848038,6.45];
ybs[2631]=['',1.8456665,-0.100434,5.2];
ybs[2632]=['',1.841247,-0.4406438,5.63];
ybs[2633]=['',1.8396927,-0.5846308,6.4];
ybs[2634]=['',1.8671135,1.043139,6.44];
ybs[2635]=['',1.8538321,0.5114536,5.93];
ybs[2636]=['',1.8647463,0.920204,6.12];
ybs[2637]=['',1.8620808,0.8332382,6.38];
ybs[2638]=['σ CMa',1.8438335,-0.4881125,3.47];
ybs[2639]=['',1.852118,0.1589192,5.97];
ybs[2640]=['19 Mon',1.8499661,-0.0745584,4.99];
ybs[2641]=['',1.8536451,0.1905648,5.13];
ybs[2642]=['ζ Gem',1.856082,0.3584367,3.79];
ybs[2643]=['',1.854694,0.2192346,5.98];
ybs[2644]=['',1.8386229,-0.8976935,5.14];
ybs[2645]=['ο2 CMa',1.8497145,-0.4165413,3.02];
ybs[2646]=['',1.8563776,0.0253926,6.57];
ybs[2647]=['',1.8550477,-0.0934957,5.62];
ybs[2648]=['',1.8543018,-0.1772799,6.45];
ybs[2649]=['γ CMa',1.8532402,-0.2734311,4.12];
ybs[2650]=['',1.8453612,-0.7581091,6.43];
ybs[2651]=['44 Gem',1.8613932,0.3945012,6.02];
ybs[2652]=['',1.8658001,0.6010817,5.55];
ybs[2653]=['',1.8387904,-1.0292488,6.02];
ybs[2654]=['',1.8317346,-1.1858989,5.17];
ybs[2655]=['',1.8623871,0.1597276,5.78];
ybs[2656]=['',1.8574702,-0.3851215,6.09];
ybs[2657]=['',1.8709158,0.5929653,5.91];
ybs[2658]=['',1.8532308,-0.7395031,5.2];
ybs[2659]=['',1.8527475,-0.7616824,5.54];
ybs[2660]=['',1.8528563,-0.7617458,6.79];
ybs[2661]=['',1.8708344,0.4911737,6.48];
ybs[2662]=['',1.8624572,-0.1866678,6.49];
ybs[2663]=['',1.8703384,0.395642,7.68];
ybs[2664]=['',1.8520399,-0.8659795,4.93];
ybs[2665]=['',1.8746141,0.5898648,6.28];
ybs[2666]=['',1.8482578,-1.0334227,5.5];
ybs[2667]=['',1.8764764,0.6529161,6.16];
ybs[2668]=['',1.8685836,0.0850925,6.11];
ybs[2669]=['',1.8601555,-0.6075794,6.14];
ybs[2670]=['',1.8661434,-0.1977241,5.39];
ybs[2671]=['',1.8657543,-0.2169171,6.48];
ybs[2672]=['',1.86244,-0.5356374,6.34];
ybs[2673]=['',1.9043662,1.252764,6.35];
ybs[2674]=['',1.871808,0.1297812,5.75];
ybs[2675]=['',1.853125,-0.9910498,5.17];
ybs[2676]=['45 Gem',1.8744791,0.2774261,5.44];
ybs[2677]=['',1.8621642,-0.6705027,6.11];
ybs[2678]=['',1.8664556,-0.4362487,6.08];
ybs[2679]=['',1.8579996,-0.8795421,6.46];
ybs[2680]=['',1.8669498,-0.4658719,6.62];
ybs[2681]=['θ Men',1.8115829,-1.3866496,5.45];
ybs[2682]=['',1.8687068,-0.4167006,5.71];
ybs[2683]=['',1.866731,-0.7143295,5.79];
ybs[2684]=['',1.8822974,0.3701944,6.43];
ybs[2685]=['δ CMa',1.8730283,-0.4612686,1.84];
ybs[2686]=['',1.8777855,-0.18122,6.21];
ybs[2687]=['',1.875005,-0.4202718,6.65];
ybs[2688]=['63 Aur',1.8898837,0.6856242,4.9];
ybs[2689]=['τ Gem',1.8871691,0.5272357,4.41];
ybs[2690]=['',1.8663568,-0.9076151,5.96];
ybs[2691]=['',1.8785073,-0.2839779,6.03];
ybs[2692]=['47 Gem',1.8880854,0.4680913,5.78];
ybs[2693]=['20 Mon',1.8818852,-0.0745883,4.92];
ybs[2694]=['',1.8743721,-0.6927463,4.83];
ybs[2695]=['',1.8983105,0.8969386,5.47];
ybs[2696]=['',1.878854,-0.4409956,5.69];
ybs[2697]=['',1.8810394,-0.3267533,6.23];
ybs[2698]=['48 Gem',1.8925714,0.4204639,5.85];
ybs[2699]=['21 Mon',1.8871055,-0.005915,5.45];
ybs[2700]=['',1.8814114,-0.4804501,5.46];
ybs[2701]=['',1.960657,1.4174337,6.31];
ybs[2702]=['',1.8893269,0.0980442,6.09];
ybs[2703]=['',1.8943463,0.4745026,6.43];
ybs[2704]=['',1.8594316,-1.2020312,6.47];
ybs[2705]=['',1.8904914,0.0949003,6.16];
ybs[2706]=['δ Mon',1.8891569,-0.0092498,4.15];
ybs[2707]=['18 Lyn',1.9102109,1.0401832,5.2];
ybs[2708]=['',1.887667,-0.365125,5.84];
ybs[2709]=['51 Gem',1.896322,0.281363,5];
ybs[2710]=['26 CMa',1.8896844,-0.4534331,5.92];
ybs[2711]=['',1.8822067,-0.854666,5.14];
ybs[2712]=['',1.8888764,-0.5385894,6.1];
ybs[2713]=['',1.9086296,0.8238079,5.58];
ybs[2714]=['',1.9013341,0.430613,6.89];
ybs[2715]=['',1.8942556,-0.1970334,5.78];
ybs[2716]=['',1.8904835,-0.4801675,6.59];
ybs[2717]=['52 Gem',1.9024537,0.4336506,5.82];
ybs[2718]=['',1.8901607,-0.6384735,5.96];
ybs[2719]=['',1.8892162,-0.7074899,5.31];
ybs[2720]=['',1.901285,0.2107884,5.62];
ybs[2721]=['',1.9000517,0.0536334,5.35];
ybs[2722]=['',1.8950365,-0.3963859,6.01];
ybs[2723]=['',1.8991547,-0.0687612,5.75];
ybs[2724]=['',1.8992697,-0.1742861,5.9];
ybs[2725]=['',1.8968019,-0.400457,6.36];
ybs[2726]=['',1.895741,-0.4781223,6.12];
ybs[2727]=['γ1 Vol',1.8697376,-1.2310247,5.69];
ybs[2728]=['γ2 Vol',1.8699337,-1.2310542,3.78];
ybs[2729]=['',1.91656,0.9091539,5.92];
ybs[2730]=['53 Gem',1.90805,0.4862177,5.71];
ybs[2731]=['',1.9001803,-0.1807311,6.03];
ybs[2732]=['',1.8900956,-0.8167594,4.49];
ybs[2733]=['',1.8963573,-0.5431808,6.6];
ybs[2734]=['',1.9873417,1.4375232,4.96];
ybs[2735]=['',1.8971269,-0.530199,6.33];
ybs[2736]=['24 Mon',1.9042509,-0.0034958,6.41];
ybs[2737]=['27 CMa',1.8986145,-0.4606068,4.66];
ybs[2738]=['',1.893105,-0.789252,4.89];
ybs[2739]=['',1.9059949,0.1385563,5.82];
ybs[2740]=['',1.8945251,-0.7797718,5.1];
ybs[2741]=['ω CMa',1.9010333,-0.4679468,3.85];
ybs[2742]=['',1.9011964,-0.4725771,5.58];
ybs[2743]=['',1.9205599,0.8626182,5.05];
ybs[2744]=['',1.9056183,-0.1854058,5.95];
ybs[2745]=['64 Aur',1.9178051,0.7128447,5.78];
ybs[2746]=['',1.8859751,-1.1035203,6.02];
ybs[2747]=['',1.9054387,-0.4150331,6.32];
ybs[2748]=['',1.9032073,-0.5362566,5.36];
ybs[2749]=['',1.9174207,0.5395774,6.24];
ybs[2750]=['',1.9077204,-0.2727105,5.46];
ybs[2751]=['',1.9008968,-0.7236913,5.94];
ybs[2752]=['',1.9131048,0.1159014,6.65];
ybs[2753]=['',1.8997168,-0.8183541,5.72];
ybs[2754]=['',1.8990485,-0.8431704,4.76];
ybs[2755]=['λ Gem',1.9169343,0.2879788,3.58];
ybs[2756]=['',1.9090424,-0.4076225,4.79];
ybs[2757]=['',1.9136543,-0.1172858,6.29];
ybs[2758]=['',1.908718,-0.4873061,4.64];
ybs[2759]=['',1.9018032,-0.9169697,5.97];
ybs[2760]=['',1.9101962,-0.5399405,6.32];
ybs[2761]=['',1.9079695,-0.6694787,5.8];
ybs[2762]=['',1.9093517,-0.639355,5.03];
ybs[2763]=['',1.906243,-0.8170531,5.66];
ybs[2764]=['47 Cam',1.9380139,1.0447441,6.35];
ybs[2765]=['π Pup',1.9107143,-0.6481668,2.7];
ybs[2766]=['',1.9140685,-0.468404,6.46];
ybs[2767]=['',1.9310611,0.7437499,6.35];
ybs[2768]=['',1.932277,0.7886463,5.77];
ybs[2769]=['δ Gem',1.925999,0.3829411,3.53];
ybs[2770]=['',1.9220292,0.0471177,5.89];
ybs[2771]=['',1.9240136,0.1239473,5.91];
ybs[2772]=['',1.925702,0.2635706,6.45];
ybs[2773]=['29 CMa',1.9179751,-0.4293405,4.98];
ybs[2774]=['τ CMa',1.9181114,-0.4362397,4.4];
ybs[2775]=['19 Lyn',1.9399126,0.9641493,6.53];
ybs[2776]=['19 Lyn',1.9399994,0.9640958,5.45];
ybs[2777]=['',1.9197644,-0.3372146,6.09];
ybs[2778]=['',1.9186826,-0.4647187,5.28];
ybs[2779]=['',1.9158125,-0.6418355,4.66];
ybs[2780]=['',1.921772,-0.2868608,5.7];
ybs[2781]=['',1.9143481,-0.768413,5.85];
ybs[2782]=['',1.9172527,-0.6419885,5.11];
ybs[2783]=['',1.9167816,-0.6850537,5.25];
ybs[2784]=['',1.9359521,0.6798707,6.4];
ybs[2785]=['65 Aur',1.9350448,0.6408546,5.13];
ybs[2786]=['',1.9199908,-0.5893624,6.3];
ybs[2787]=['56 Gem',1.9338921,0.3560717,5.1];
ybs[2788]=['',1.9284005,-0.2513563,5.45];
ybs[2789]=['',2.0005089,1.4110587,6.41];
ybs[2790]=['',1.9299561,-0.155686,6.55];
ybs[2791]=['',1.9277127,-0.3995628,6.61];
ybs[2792]=['',1.9276726,-0.4713299,6.01];
ybs[2793]=['',1.9336568,0.0023563,5.99];
ybs[2794]=['',1.9283948,-0.4526175,5.87];
ybs[2795]=['δ Vol',1.9059612,-1.1867645,3.98];
ybs[2796]=['',1.94876,0.904839,5.8];
ybs[2797]=['66 Aur',1.9444138,0.709108,5.19];
ybs[2798]=['',1.9332327,-0.1574523,6.43];
ybs[2799]=['',1.9346399,-0.0527302,6.23];
ybs[2800]=['57 Gem',1.9407432,0.4364649,5.03];
ybs[2801]=['',1.9615005,1.1569193,6.47];
ybs[2802]=['58 Gem',1.9406356,0.399721,6.02];
ybs[2803]=['',1.9350386,-0.1051589,5.82];
ybs[2804]=['',1.9337042,-0.332641,4.96];
ybs[2805]=['',1.9236861,-0.9137306,6.05];
ybs[2806]=['',1.9237081,-0.9136967,6.6];
ybs[2807]=['',1.9249659,-0.9097964,5.39];
ybs[2808]=['59 Gem',1.9455707,0.481616,5.76];
ybs[2809]=['',1.9446777,0.2700688,6.41];
ybs[2810]=['21 Lyn',1.9562077,0.8581223,4.64];
ybs[2811]=['',1.9365866,-0.5579207,5.43];
ybs[2812]=['1 CMi',1.946757,0.2029131,5.3];
ybs[2813]=['ι Gem',1.9506807,0.4843986,3.79];
ybs[2814]=['',1.9388456,-0.4865459,5.38];
ybs[2815]=['',1.9388506,-0.562783,5.39];
ybs[2816]=['',1.9405728,-0.5281366,6.6];
ybs[2817]=['',1.9444733,-0.2835212,5.33];
ybs[2818]=['',1.9425502,-0.4006584,6.19];
ybs[2819]=['η CMa',1.9414435,-0.5121879,2.45];
ybs[2820]=['ε CMi',1.9496324,0.161131,4.99];
ybs[2821]=['',1.940591,-0.6262389,6.31];
ybs[2822]=['',1.9770036,1.1941332,5.64];
ybs[2823]=['',1.9451392,-0.3325858,6.24];
ybs[2824]=['',1.9466137,-0.2407793,5.78];
ybs[2825]=['',1.9499985,-0.1015617,5.97];
ybs[2826]=['',1.9441004,-0.5559282,5.35];
ybs[2827]=['',1.9552651,0.375093,6.54];
ybs[2828]=['',1.9532418,0.1843757,6.37];
ybs[2829]=['61 Gem',1.955659,0.3527762,5.93];
ybs[2830]=['',1.9509498,-0.079965,6.76];
ybs[2831]=['',1.9471483,-0.3844358,6.05];
ybs[2832]=['',1.9542301,0.1913697,6.41];
ybs[2833]=['',1.9474103,-0.4408978,5.78];
ybs[2834]=['',1.9440726,-0.6515919,6.97];
ybs[2835]=['',1.9440798,-0.6516064,6.84];
ybs[2836]=['',1.9654863,0.8401709,5.72];
ybs[2837]=['β CMi',1.9561496,0.1438978,2.9];
ybs[2838]=['63 Gem',1.9592018,0.3735001,5.22];
ybs[2839]=['',1.9484108,-0.5547101,6.31];
ybs[2840]=['',1.7409416,-1.5172149,6.47];
ybs[2841]=['22 Lyn',1.9702813,0.866143,5.36];
ybs[2842]=['',1.952957,-0.4146316,6.56];
ybs[2843]=['η CMi',1.9599644,0.120372,5.25];
ybs[2844]=['ρ Gem',1.9656227,0.5539454,4.18];
ybs[2845]=['',1.9551697,-0.3125727,5.63];
ybs[2846]=['γ CMi',1.9605925,0.1549914,4.32];
ybs[2847]=['',1.9543426,-0.4037066,5.61];
ybs[2848]=['',1.9526098,-0.596645,5.9];
ybs[2849]=['64 Gem',1.966454,0.4899531,5.05];
ybs[2850]=['',1.9635374,0.2629152,6.22];
ybs[2851]=['',1.9585767,-0.2024924,5.79];
ybs[2852]=['',1.957516,-0.3997614,5.95];
ybs[2853]=['65 Gem',1.9685024,0.4864246,5.01];
ybs[2854]=['',1.9500417,-0.891209,5.1];
ybs[2855]=['',1.9584379,-0.5096514,5.54];
ybs[2856]=['6 CMi',1.9678257,0.2087531,4.54];
ybs[2857]=['',1.9652351,-0.0340516,5.59];
ybs[2858]=['',1.9655524,-0.1325908,5.86];
ybs[2859]=['',1.9651947,-0.1810327,5.75];
ybs[2860]=['',1.9650147,-0.2625831,6.05];
ybs[2861]=['',1.9597327,-0.6607028,6.58];
ybs[2862]=['',1.9621059,-0.5566516,6.38];
ybs[2863]=['',1.9621205,-0.5566273,7.13];
ybs[2864]=['',1.9782657,0.6780483,6.54];
ybs[2865]=['',1.9631136,-0.5498128,5.77];
ybs[2866]=['',1.9668574,-0.4026545,4.85];
ybs[2867]=['',1.9627994,-0.6781961,5.43];
ybs[2868]=['',1.9718427,-0.0920288,6.24];
ybs[2869]=['',1.9767763,0.2973891,5.42];
ybs[2870]=['σ Pup',1.9631219,-0.7565477,3.25];
ybs[2871]=['',1.9815185,0.3986384,6.54];
ybs[2872]=['δ1 CMi',1.9775226,0.0325916,5.25];
ybs[2873]=['',1.9702375,-0.5412016,4.65];
ybs[2874]=['',1.9699161,-0.6525097,6.65];
ybs[2875]=['',1.9771509,-0.1558212,5.9];
ybs[2876]=['',1.9657443,-0.9197368,5.87];
ybs[2877]=['',1.9731364,-0.6317998,6.68];
ybs[2878]=['68 Gem',1.9845898,0.2753925,5.25];
ybs[2879]=['δ2 CMi',1.9823465,0.0565952,5.59];
ybs[2880]=['',1.9592473,-1.1267024,6.39];
ybs[2881]=['',1.9743869,-0.6271771,6.61];
ybs[2882]=['α Gem',1.9895627,0.5557175,2.88];
ybs[2883]=['α Gem',1.9895627,0.5557126,1.98];
ybs[2884]=['',1.9678562,-0.9502552,5.96];
ybs[2885]=['',1.9864838,0.1836135,6.28];
ybs[2886]=['',2.0006933,0.9722497,5.92];
ybs[2887]=['',1.9772575,-0.6284673,6.3];
ybs[2888]=['',1.991897,0.539525,5.33];
ybs[2889]=['',1.9825099,-0.251083,6.21];
ybs[2890]=['',1.9959703,0.7501792,6.3];
ybs[2891]=['',1.9821353,-0.3396433,5.66];
ybs[2892]=['',1.9812157,-0.4321152,5.85];
ybs[2893]=['δ3 CMi',1.9870177,0.0580019,5.81];
ybs[2894]=['',1.98438,-0.2543252,4.97];
ybs[2895]=['',1.9987635,0.8051373,5.65];
ybs[2896]=['',1.9891921,0.0467161,6.55];
ybs[2897]=['υ Gem',1.9951076,0.4685662,4.06];
ybs[2898]=['',1.9852037,-0.3899778,4.45];
ybs[2899]=['',1.9807617,-0.6999891,6.26];
ybs[2900]=['',1.9805658,-0.7528287,6.52];
ybs[2901]=['',1.9862777,-0.4105312,5.83];
ybs[2902]=['',1.986314,-0.4105507,5.87];
ybs[2903]=['',1.9836924,-0.6350586,5.54];
ybs[2904]=['',1.9869126,-0.4566625,6.65];
ybs[2905]=['',1.9854216,-0.5848835,6.11];
ybs[2906]=['',2.004931,0.850388,5.92];
ybs[2907]=['',2.0017362,0.6977063,6.38];
ybs[2908]=['',1.987319,-0.4722889,5.77];
ybs[2909]=['',1.9840452,-0.6973241,6.76];
ybs[2910]=['',1.9972015,0.1014461,5.91];
ybs[2911]=['ε Men',1.9391247,-1.381211,5.53];
ybs[2912]=['',1.9954114,-0.1459174,6.27];
ybs[2913]=['',1.9942734,-0.2538011,5.7];
ybs[2914]=['',1.9907462,-0.4959884,4.64];
ybs[2915]=['',1.9942733,-0.3876293,6.34];
ybs[2916]=['70 Gem',2.0069243,0.610837,5.56];
ybs[2917]=['',1.9862252,-0.8992442,6.28];
ybs[2918]=['',2.0051228,0.4242934,6.27];
ybs[2919]=['25 Mon',1.9999237,-0.0726172,5.13];
ybs[2920]=['',1.9967823,-0.3447279,5.74];
ybs[2921]=['23 Lyn',2.0184481,0.9953854,6.06];
ybs[2922]=['ο Gem',2.0095979,0.6027258,4.9];
ybs[2923]=['',2.0092899,0.4218808,6.17];
ybs[2924]=['',2.0011854,-0.2529124,6.53];
ybs[2925]=['',1.9992296,-0.4158161,6.37];
ybs[2926]=['',1.9905144,-0.9177384,4.94];
ybs[2927]=['',2.0144918,0.6683459,5.73];
ybs[2928]=['',2.0126886,0.5577873,6.17];
ybs[2929]=['',1.9991225,-0.6111818,4.53];
ybs[2930]=['74 Gem',2.0102556,0.3075984,5.05];
ybs[2931]=['',2.0193037,0.8391567,5.56];
ybs[2932]=['',1.995496,-0.8531074,5.72];
ybs[2933]=['',1.9917884,-0.9762727,6.39];
ybs[2934]=['',2.000759,-0.6165712,6.6];
ybs[2935]=['α CMi',2.009062,0.0903117,0.38];
ybs[2936]=['',2.0036182,-0.4435704,4.7];
ybs[2937]=['',2.0006471,-0.6642767,6.38];
ybs[2938]=['24 Lyn',2.0281427,1.0237726,4.99];
ybs[2939]=['',2.0074481,-0.3268925,6.72];
ybs[2940]=['',2.0058363,-0.4686541,4.5];
ybs[2941]=['',2.0058727,-0.4686882,4.62];
ybs[2942]=['',2.012618,0.0904067,6.02];
ybs[2943]=['',2.0169917,0.4008544,5.89];
ybs[2944]=['',2.0033764,-0.6988541,6.59];
ybs[2945]=['',2.0158383,0.239452,6.24];
ybs[2946]=['',2.0049979,-0.6378674,5.8];
ybs[2947]=['',2.0040548,-0.6777321,6.19];
ybs[2948]=['',2.0085693,-0.4697308,6.5];
ybs[2949]=['',2.002378,-0.8491208,5.68];
ybs[2950]=['',2.0142458,-0.1437618,6.01];
ybs[2951]=['',2.0131072,-0.2672904,4.94];
ybs[2952]=['',2.0122493,-0.344035,5.93];
ybs[2953]=['',2.0080684,-0.669488,4.84];
ybs[2954]=['',2.0250936,0.5925057,6.02];
ybs[2955]=['',2.0092709,-0.6665427,5.73];
ybs[2956]=['',2.0095555,-0.6686618,5.76];
ybs[2957]=['',2.0205178,0.2343769,5.77];
ybs[2958]=['',2.0189784,0.0623627,5.94];
ybs[2959]=['',2.0213719,0.2470774,5.56];
ybs[2960]=['',2.010334,-0.6567708,6];
ybs[2961]=['',2.0319615,0.8793144,5.27];
ybs[2962]=['α Mon',2.017058,-0.1675957,3.93];
ybs[2963]=['',2.0050587,-0.9306719,6.06];
ybs[2964]=['',2.0140898,-0.4886393,6.76];
ybs[2965]=['σ Gem',2.0274213,0.5031985,4.28];
ybs[2966]=['',2.0146035,-0.5534794,6.56];
ybs[2967]=['51 Cam',2.0451067,1.1414739,5.92];
ybs[2968]=['',2.0172403,-0.3907563,6.18];
ybs[2969]=['49 Cam',2.0437364,1.0956564,6.49];
ybs[2970]=['',2.027414,0.3900282,6.21];
ybs[2971]=['',1.9848973,-1.2971955,7.16];
ybs[2972]=['',1.9849046,-1.2971955,7.26];
ybs[2973]=['',2.0159449,-0.6734348,5.42];
ybs[2974]=['',2.0254247,0.0023935,6.19];
ybs[2975]=['76 Gem',2.0307963,0.4490965,5.31];
ybs[2976]=['',2.016021,-0.7798761,6.41];
ybs[2977]=['κ Gem',2.0321876,0.4249016,3.57];
ybs[2978]=['',2.0190162,-0.6733582,6.54];
ybs[2979]=['',2.0308351,0.223517,6.43];
ybs[2980]=['',2.0232286,-0.4608233,5.64];
ybs[2981]=['',2.0300061,0.0410538,6.47];
ybs[2982]=['β Gem',2.0361219,0.4882159,1.14];
ybs[2983]=['79 Gem',2.0351212,0.3536575,6.33];
ybs[2984]=['',2.0269735,-0.4460434,6.55];
ybs[2985]=['1 Pup',2.026366,-0.496783,4.59];
ybs[2986]=['',2.0245302,-0.6301085,5.6];
ybs[2987]=['',2.0240196,-0.6792194,6.89];
ybs[2988]=['3 Pup',2.0275142,-0.506273,3.96];
ybs[2989]=['',2.093573,1.3998691,6.56];
ybs[2990]=['',2.0229294,-0.7893334,5.06];
ybs[2991]=['',2.0424045,0.6538606,5.18];
ybs[2992]=['',1.9860598,-1.3558184,6.18];
ybs[2993]=['',2.0266677,-0.6676663,6.4];
ybs[2994]=['',2.0264441,-0.7153474,5.17];
ybs[2995]=['81 Gem',2.0392787,0.3221222,4.88];
ybs[2996]=['',2.0310008,-0.4315647,5.62];
ybs[2997]=['',2.0232865,-0.8734496,6.57];
ybs[2998]=['',2.0183116,-1.0242035,6.43];
ybs[2999]=['',2.0287211,-0.6303345,5.8];
ybs[3000]=['11 CMi',2.0396488,0.1870037,5.3];
ybs[3001]=['2 Pup',2.0353495,-0.2572529,6.89];
ybs[3002]=['2 Pup',2.0353785,-0.2573354,6.07];
ybs[3003]=['',2.0304147,-0.6631548,5.88];
ybs[3004]=['',2.0214986,-1.0172135,6.21];
ybs[3005]=['π Gem',2.0459098,0.5822615,5.14];
ybs[3006]=['',2.038054,-0.1191391,5.49];
ybs[3007]=['4 Pup',2.037405,-0.2551236,5.04];
ybs[3008]=['',2.03263,-0.6621942,6.54];
ybs[3009]=['',2.0334044,-0.6636065,3.61];
ybs[3010]=['',2.035028,-0.5973644,5.37];
ybs[3011]=['',2.0409616,-0.2221676,6.39];
ybs[3012]=['',2.0332767,-0.7645496,6.03];
ybs[3013]=['82 Gem',2.0500791,0.4029304,6.18];
ybs[3014]=['',2.0374161,-0.663008,5.88];
ybs[3015]=['',2.0426156,-0.3939891,5.9];
ybs[3016]=['ζ Vol',2.0138453,-1.2681124,3.95];
ybs[3017]=['',2.0389752,-0.7001138,6.57];
ybs[3018]=['',2.0447456,-0.2800422,6.34];
ybs[3019]=['',2.0452321,-0.2804552,6.43];
ybs[3020]=['',2.0628964,0.9437519,6.02];
ybs[3021]=['5 Pup',2.0462005,-0.2137611,5.48];
ybs[3022]=['',2.051792,0.2324034,6.04];
ybs[3023]=['',2.0334896,-0.9909251,6.12];
ybs[3024]=['',2.0413733,-0.6874063,6.31];
ybs[3025]=['',2.0512663,0.0746602,6.53];
ybs[3026]=['ο Pup',2.0463161,-0.4536427,4.5];
ybs[3027]=['',2.0428148,-0.6730925,5.08];
ybs[3028]=['',2.0283528,-1.1540953,6.38];
ybs[3029]=['',2.0427992,-0.8144209,5.23];
ybs[3030]=['',2.0252209,-1.2195304,6.18];
ybs[3031]=['',2.0696012,0.9625941,6.38];
ybs[3032]=['',2.0613087,0.5790573,6.03];
ybs[3033]=['',2.0458634,-0.7104677,6.14];
ybs[3034]=['',2.0528457,-0.2340241,6.23];
ybs[3035]=['',2.0504676,-0.4357607,5.33];
ybs[3036]=['6 Pup',2.0536274,-0.3016571,5.18];
ybs[3037]=['ξ Pup',2.0516336,-0.4348466,3.34];
ybs[3038]=['',2.0463191,-0.822616,4.71];
ybs[3039]=['',2.0560408,-0.1612496,5.61];
ybs[3040]=['',2.053811,-0.3536442,6.56];
ybs[3041]=['',2.050964,-0.616074,5.93];
ybs[3042]=['',2.0591219,0.0562227,6.18];
ybs[3043]=['',2.0553273,-0.3417205,6.12];
ybs[3044]=['',2.052563,-0.5819653,5.6];
ybs[3045]=['',2.0647016,0.3363043,5.99];
ybs[3046]=['',2.0592199,-0.1952071,6.16];
ybs[3047]=['',2.0502987,-0.8103285,4.11];
ybs[3048]=['',2.0454309,-0.9865597,6.33];
ybs[3049]=['',2.0514177,-0.7820319,6.32];
ybs[3050]=['',2.0501635,-0.8187834,5.84];
ybs[3051]=['ζ CMi',2.0630435,0.0298562,5.14];
ybs[3052]=['',2.0590894,-0.4290766,6.45];
ybs[3053]=['',2.0649247,0.056212,6.31];
ybs[3054]=['',2.0488966,-0.9855092,5.59];
ybs[3055]=['8 Pup',2.0624871,-0.2247238,6.36];
ybs[3056]=['9 Pup',2.0628432,-0.2435497,5.17];
ybs[3057]=['25 Lyn',2.0771415,0.8260367,6.25];
ybs[3058]=['26 Lyn',2.0781274,0.8291522,5.45];
ybs[3059]=['φ Gem',2.0717516,0.4661538,4.97];
ybs[3060]=['',2.0623379,-0.3705364,5.63];
ybs[3061]=['',2.0568299,-0.779036,6.45];
ybs[3062]=['',2.0488901,-1.0531072,5.78];
ybs[3063]=['',2.0550606,-0.8825312,5.91];
ybs[3064]=['',2.0676043,-0.0957289,5.76];
ybs[3065]=['10 Pup',2.0651821,-0.2601057,5.69];
ybs[3066]=['',2.0596953,-0.7531376,6.32];
ybs[3067]=['',2.1062941,1.2890573,5.41];
ybs[3068]=['',2.0519727,-1.049055,6.72];
ybs[3069]=['',2.086448,0.9851653,6.72];
ybs[3070]=['',2.061149,-0.7495236,6.04];
ybs[3071]=['',2.0641615,-0.6067075,5.01];
ybs[3072]=['',2.0636629,-0.7091674,3.73];
ybs[3073]=['',2.0499798,-1.1562978,5.79];
ybs[3074]=['',2.1296333,1.3860879,5.42];
ybs[3075]=['',2.0816501,0.6170538,6.23];
ybs[3076]=['',2.0656266,-0.6792773,4.49];
ybs[3077]=['',2.0675596,-0.6356621,5.43];
ybs[3078]=['85 Gem',2.0809475,0.3460242,5.35];
ybs[3079]=['',2.079954,0.153671,5.86];
ybs[3080]=['',2.063885,-0.9498738,5.7];
ybs[3081]=['',2.0667654,-0.8669029,4.63];
ybs[3082]=['',2.0679332,-0.8405505,4.24];
ybs[3083]=['',2.0724957,-0.627182,5.49];
ybs[3084]=['',2.0746485,-0.6091993,6.15];
ybs[3085]=['',2.0836301,0.0772721,6.17];
ybs[3086]=['',2.0934319,0.7665153,6.34];
ybs[3087]=['1 Cnc',2.0865833,0.2745668,5.78];
ybs[3088]=['',2.0772864,-0.5406223,6.44];
ybs[3089]=['',2.087545,0.1497932,6.05];
ybs[3090]=['',2.0873261,0.0186416,6.35];
ybs[3091]=['',2.0823418,-0.5295971,6.33];
ybs[3092]=['',2.0750279,-0.9187546,6.38];
ybs[3093]=['',2.079012,-0.7662534,6.02];
ybs[3094]=['11 Pup',2.0847254,-0.4003546,4.2];
ybs[3095]=['',2.0911271,0.1248672,6.41];
ybs[3096]=['',2.093304,0.2872665,5.99];
ybs[3097]=['',2.0740364,-1.0011327,5.63];
ybs[3098]=['',2.1080678,1.0295111,5.77];
ybs[3099]=['',2.0819448,-0.7120031,6.78];
ybs[3100]=['',2.1893092,1.4658918,6.49];
ybs[3101]=['53 Cam',2.1097999,1.051795,6.01];
ybs[3102]=['14 CMi',2.0920441,0.0377928,5.29];
ybs[3103]=['',2.0842954,-0.7411497,6.09];
ybs[3104]=['',2.1136981,1.1000613,6.4];
ybs[3105]=['',2.0879644,-0.5304703,4.79];
ybs[3106]=['',2.0842303,-0.7602464,5.35];
ybs[3107]=['',2.0978169,0.2300744,6.02];
ybs[3108]=['',2.0856854,-0.7708859,5.09];
ybs[3109]=['χ Car',2.0827264,-0.9257354,3.47];
ybs[3110]=['',2.0855523,-0.836869,6.22];
ybs[3111]=['',2.1133204,0.9985413,6.49];
ybs[3112]=['',2.0798695,-1.0574015,5.74];
ybs[3113]=['',2.0880239,-0.7965126,5.17];
ybs[3114]=['27 Mon',2.0979297,-0.06527,4.93];
ybs[3115]=['12 Pup',2.0944728,-0.4078869,5.11];
ybs[3116]=['ω1 Cnc',2.1041226,0.4421306,5.83];
ybs[3117]=['',2.1033424,0.3448006,6.25];
ybs[3118]=['',2.0823849,-1.0329713,6.25];
ybs[3119]=['',2.1044027,0.4105444,6.34];
ybs[3120]=['3 Cnc',2.1032037,0.3010367,5.55];
ybs[3121]=['',2.0894146,-0.8605205,4.41];
ybs[3122]=['',2.1088498,0.6170093,6.34];
ybs[3123]=['',2.098019,-0.3221733,4.61];
ybs[3124]=['ω2 Cnc',2.1075926,0.436835,6.31];
ybs[3125]=['',2.0897524,-0.8989815,6.44];
ybs[3126]=['5 Cnc',2.1063005,0.2861376,5.99];
ybs[3127]=['',2.1023185,-0.0513492,6.51];
ybs[3128]=['',2.1047229,0.0841085,5.65];
ybs[3129]=['',2.0931431,-0.7902096,5.99];
ybs[3130]=['',2.0863316,-1.0535202,5.6];
ybs[3131]=['',2.0833971,-1.1057636,6.14];
ybs[3132]=['',2.0954224,-0.6869093,5.24];
ybs[3133]=['28 Mon',2.1044889,-0.0253622,4.68];
ybs[3134]=['',2.0935751,-0.8732981,6.32];
ybs[3135]=['',2.0936554,-0.8732449,6.34];
ybs[3136]=['',2.1075345,0.154513,6.22];
ybs[3137]=['',2.1091583,0.0396769,4.39];
ybs[3138]=['',2.0987915,-0.7944231,6.61];
ybs[3139]=['',2.0908966,-1.0626236,5.81];
ybs[3140]=['',2.0982076,-0.8559355,6.02];
ybs[3141]=['χ Gem',2.115497,0.4840224,4.94];
ybs[3142]=['',2.1096202,-0.1116732,6.33];
ybs[3143]=['',2.099235,-0.8540174,6.12];
ybs[3144]=['',2.0945983,-1.0518625,6.33];
ybs[3145]=['',2.0943596,-1.0584847,5.17];
ybs[3146]=['',2.1049144,-0.6517822,5.95];
ybs[3147]=['',2.107022,-0.6477184,6.34];
ybs[3148]=['',2.1003321,-0.9461732,5.87];
ybs[3149]=['',2.1027018,-0.9525237,6.1];
ybs[3150]=['',2.1205513,0.3277721,6.15];
ybs[3151]=['',2.0970176,-1.1105104,4.82];
ybs[3152]=['',2.1114492,-0.5676737,5.82];
ybs[3153]=['',2.1032382,-0.9689308,6.28];
ybs[3154]=['',2.1096087,-0.7220646,5.52];
ybs[3155]=['8 Cnc',2.1217615,0.2278643,5.12];
ybs[3156]=['',2.124637,0.4793906,6.21];
ybs[3157]=['ζ Pup',2.1133411,-0.6992656,2.25];
ybs[3158]=['',2.1127692,-0.7506695,6.29];
ybs[3159]=['28 Lyn',2.1321322,0.753928,6.26];
ybs[3160]=['14 Pup',2.1190332,-0.3454046,6.13];
ybs[3161]=['μ1 Cnc',2.1274609,0.3939661,5.99];
ybs[3162]=['',2.1166868,-0.5713676,5.31];
ybs[3163]=['',2.0899768,-1.2793993,6.34];
ybs[3164]=['',2.124608,-0.0111059,6.41];
ybs[3165]=['27 Lyn',2.13832,0.8978441,4.84];
ybs[3166]=['',2.1270928,-0.1624548,6.23];
ybs[3167]=['',2.1459488,1.0154907,5.93];
ybs[3168]=['μ2 Cnc',2.133776,0.375561,5.3];
ybs[3169]=['',2.123102,-0.5869853,6.14];
ybs[3170]=['',2.117564,-0.8840556,5.95];
ybs[3171]=['',2.1205908,-0.8210253,6.19];
ybs[3172]=['',2.1189308,-0.9279969,5.53];
ybs[3173]=['',2.1417612,0.7394295,6.27];
ybs[3174]=['',2.1595228,1.1939476,5.32];
ybs[3175]=['',2.1303896,-0.3598529,5.38];
ybs[3176]=['12 Cnc',2.1376205,0.2369604,6.27];
ybs[3177]=['ρ Pup',2.1313159,-0.4252948,2.81];
ybs[3178]=['',2.1162898,-1.0977797,6.3];
ybs[3179]=['',2.1265261,-0.791147,5.05];
ybs[3180]=['ζ Mon',2.1366056,-0.0531945,4.34];
ybs[3181]=['',2.1378928,-0.1990338,6.32];
ybs[3182]=['',2.1366167,-0.3565187,6.36];
ybs[3183]=['ψ Cnc',2.1456514,0.4440541,5.73];
ybs[3184]=['16 Pup',2.1379724,-0.3370073,4.4];
ybs[3185]=['',2.1501739,0.6748519,6.58];
ybs[3186]=['',2.140031,-0.2847188,5.68];
ybs[3187]=['',2.135478,-0.6587792,6.37];
ybs[3188]=['',2.1379271,-0.5303465,6.65];
ybs[3189]=['',2.2188715,1.4374428,6.32];
ybs[3190]=['',2.1475693,0.2541974,6.23];
ybs[3191]=['',2.1379572,-0.6199257,6.2];
ybs[3192]=['',2.1621958,0.9841188,5.85];
ybs[3193]=['',2.1487093,0.170274,6.07];
ybs[3194]=['18 Pup',2.1452984,-0.2419723,5.54];
ybs[3195]=['',2.1371616,-0.8508224,5.7];
ybs[3196]=['',2.1393729,-0.7712098,5.21];
ybs[3197]=['',2.1403184,-0.7453418,6.26];
ybs[3198]=['γ1 Vel',2.1386883,-0.8274618,4.27];
ybs[3199]=['γ2 Vel',2.1388855,-0.8273021,1.78];
ybs[3200]=['ζ1 Cnc',2.1530417,0.3068677,5.63];
ybs[3201]=['ζ1 Cnc',2.1530417,0.3068677,6.02];
ybs[3202]=['ζ2 Cnc',2.1530853,0.3068676,6.2];
ybs[3203]=['19 Pup',2.1479816,-0.2267538,4.72];
ybs[3204]=['',2.1493595,-0.1367941,5.36];
ybs[3205]=['',2.1396602,-0.83779,5.23];
ybs[3206]=['',2.1536049,0.2432688,6.54];
ybs[3207]=['15 Cnc',2.1575647,0.5164549,5.64];
ybs[3208]=['',2.191155,1.3210035,5.54];
ybs[3209]=['',2.1322912,-1.1146508,6.28];
ybs[3210]=['',2.1382914,-0.9799986,5.66];
ybs[3211]=['',2.145978,-0.6520055,6.44];
ybs[3212]=['',2.1352939,-1.0710468,4.76];
ybs[3213]=['',2.171361,1.0526657,6.45];
ybs[3214]=['',2.1564279,0.2870765,6.01];
ybs[3215]=['ε Vol',2.1292718,-1.1987033,4.35];
ybs[3216]=['',2.1597014,0.4026748,6.56];
ybs[3217]=['',2.1473102,-0.692611,4.45];
ybs[3218]=['',2.1474474,-0.7514048,4.75];
ybs[3219]=['',2.1460243,-0.8469545,5.82];
ybs[3220]=['',2.1616496,0.3073424,6.47];
ybs[3221]=['20 Pup',2.1568869,-0.2767101,4.99];
ybs[3222]=['',2.1539025,-0.5231895,6.52];
ybs[3223]=['',2.1622185,0.2265762,6.38];
ybs[3224]=['',2.1497071,-0.8152347,5.76];
ybs[3225]=['',2.1539464,-0.6630538,6.43];
ybs[3226]=['',2.1519814,-0.8086063,6.03];
ybs[3227]=['29 Lyn',2.1799657,1.0385233,5.64];
ybs[3228]=['',2.1947557,1.2625325,5.98];
ybs[3229]=['',2.1568161,-0.6277207,4.78];
ybs[3230]=['',2.1577581,-0.5870463,6.37];
ybs[3231]=['',2.1599957,-0.5621211,6.06];
ybs[3232]=['',2.1588925,-0.6351031,5.08];
ybs[3233]=['',2.1589208,-0.635428,6.11];
ybs[3234]=['',2.1600083,-0.6205848,5.78];
ybs[3235]=['',2.1590208,-0.7053627,4.44];
ybs[3236]=['',2.1566888,-0.8213168,5.13];
ybs[3237]=['',2.1866458,1.0897569,5.71];
ybs[3238]=['',2.1812634,0.9437929,6.27];
ybs[3239]=['',2.1563124,-0.8772395,5.51];
ybs[3240]=['',2.1718389,0.2034872,7.13];
ybs[3241]=['β Cnc',2.1715416,0.1591417,3.52];
ybs[3242]=['',2.160236,-0.8011207,5.83];
ybs[3243]=['',2.1674329,-0.5409281,6.21];
ybs[3244]=['',2.1759524,0.1535587,6.29];
ybs[3245]=['',2.1676824,-0.6277928,6.16];
ybs[3246]=['30 Lyn',2.1910648,1.0066034,5.89];
ybs[3247]=['',2.1722835,-0.3732876,6.6];
ybs[3248]=['',2.1642459,-0.8816748,6.44];
ybs[3249]=['21 Pup',2.1745553,-0.2854093,6.16];
ybs[3250]=['',2.1908927,0.9338426,6.49];
ybs[3251]=['',2.1790952,-0.2216591,5.98];
ybs[3252]=['',2.1624119,-1.0992523,5.16];
ybs[3253]=['',2.1766251,-0.5248433,6.45];
ybs[3254]=['χ Cnc',2.1876302,0.4738362,5.14];
ybs[3255]=['',2.2014059,1.0569873,6.41];
ybs[3256]=['',2.1886419,0.3609115,5.83];
ybs[3257]=['',2.1828963,-0.1786238,6.32];
ybs[3258]=['',2.1777836,-0.6199368,5.58];
ybs[3259]=['',2.177345,-0.6534902,6.7];
ybs[3260]=['λ Cnc',2.1895627,0.4180599,5.98];
ybs[3261]=['',2.1858521,0.0677005,6.05];
ybs[3262]=['',2.1788854,-0.6410185,4.45];
ybs[3263]=['',2.1873942,-0.0170768,6.18];
ybs[3264]=['',2.1875529,-0.0942159,6.13];
ybs[3265]=['',2.1830624,-0.6049118,6.43];
ybs[3266]=['',2.1745506,-1.0338424,6.42];
ybs[3267]=['31 Lyn',2.2004264,0.752549,4.25];
ybs[3268]=['',2.1877393,-0.4013171,6.13];
ybs[3269]=['',2.2053402,0.9276269,5.51];
ybs[3270]=['',2.1922534,-0.0291764,6.5];
ybs[3271]=['',2.1917515,-0.3516595,5.58];
ybs[3272]=['',2.1753186,-1.1463552,5.07];
ybs[3273]=['',2.194264,-0.3081565,5.75];
ybs[3274]=['',2.1913967,-0.5781206,4.83];
ybs[3275]=['',2.1910999,-0.637985,5.2];
ybs[3276]=['20 Cnc',2.2017039,0.3187298,5.95];
ybs[3277]=['',2.1972079,-0.1090676,6.15];
ybs[3278]=['',2.1911822,-0.6927256,6.16];
ybs[3279]=['',2.2085487,0.7318871,6.02];
ybs[3280]=['',2.198898,-0.1328797,5.96];
ybs[3281]=['22 Pup',2.1982058,-0.2290706,6.11];
ybs[3282]=['21 Cnc',2.2038896,0.1843307,6.08];
ybs[3283]=['',2.1979848,-0.461083,5.9];
ybs[3284]=['',2.2098183,0.6098231,6.06];
ybs[3285]=['',2.1889783,-1.0130297,5.97];
ybs[3286]=['',2.1955754,-0.8475345,4.82];
ybs[3287]=['',2.20643,-0.0835625,6.01];
ybs[3288]=['',2.1994734,-0.6694394,6.32];
ybs[3289]=['1 Hya',2.2063567,-0.0667054,5.61];
ybs[3290]=['',2.1878761,-1.1200705,6.12];
ybs[3291]=['25 Cnc',2.2124281,0.2962649,6.14];
ybs[3292]=['',2.1970296,-0.9109557,5.85];
ybs[3293]=['κ1 Vol',2.1805371,-1.2493688,5.37];
ybs[3294]=['κ2 Vol',2.1813939,-1.2492006,5.65];
ybs[3295]=['',2.2330436,1.1732857,5.88];
ybs[3296]=['φ1 Cnc',2.2155484,0.4855846,5.57];
ybs[3297]=['',2.2109285,0.035447,5.73];
ybs[3298]=['',2.212493,0.1307783,5.13];
ybs[3299]=['ε Car',2.1945503,-1.0398592,1.86];
ybs[3300]=['',2.2072091,-0.4053448,5.68];
ybs[3301]=['',2.2213918,0.7955362,6.32];
ybs[3302]=['φ2 Cnc',2.216901,0.4688416,6.32];
ybs[3303]=['φ2 Cnc',2.2169155,0.4688562,6.3];
ybs[3304]=['24 Cnc',2.2163059,0.4269498,7.02];
ybs[3305]=['24 Cnc',2.2163278,0.4269692,7.81];
ybs[3306]=['',2.2110504,-0.0694234,3.9];
ybs[3307]=['',2.2078037,-0.420923,5.28];
ybs[3308]=['',2.2090275,-0.3685602,6.01];
ybs[3309]=['',2.210623,-0.3056194,6.44];
ybs[3310]=['α Cha',2.172616,-1.3436874,4.07];
ybs[3311]=['27 Cnc',2.2162164,0.2196094,5.5];
ybs[3312]=['',2.2118813,-0.2618186,5.98];
ybs[3313]=['2 Hya',2.2145098,-0.0708449,5.59];
ybs[3314]=['',2.206524,-0.7477053,5.98];
ybs[3315]=['ο UMa',2.2341989,1.0584499,3.36];
ybs[3316]=['',2.2153237,-0.2200188,5.54];
ybs[3317]=['',2.2180751,-0.1131267,6.59];
ybs[3318]=['',2.2105396,-0.7369589,5.47];
ybs[3319]=['',2.212568,-0.6829636,6.53];
ybs[3320]=['',2.2126116,-0.6829783,7.25];
ybs[3321]=['28 Cnc',2.2247959,0.4201386,6.1];
ybs[3322]=['',2.2084082,-0.9040662,5.17];
ybs[3323]=['',2.2154073,-0.5111549,6.73];
ybs[3324]=['',2.2470367,1.2085622,6.31];
ybs[3325]=['29 Cnc',2.2245085,0.2467597,5.95];
ybs[3326]=['η Vol',2.189798,-1.2822846,5.29];
ybs[3327]=['',2.2187982,-0.365052,6.56];
ybs[3328]=['',2.2171738,-0.5540542,6.33];
ybs[3329]=['',2.2234255,-0.0451985,6.39];
ybs[3330]=['',2.2225498,-0.1551336,6.43];
ybs[3331]=['',2.2200847,-0.457358,6.62];
ybs[3332]=['θ Cha',2.1815456,-1.3535593,4.35];
ybs[3333]=['',2.2122714,-0.9229126,6.05];
ybs[3334]=['',2.2247914,-0.1714077,6];
ybs[3335]=['',2.2201521,-0.6141132,5.75];
ybs[3336]=['',2.2232776,-0.4039415,6.51];
ybs[3337]=['',2.2246152,-0.3669185,6.67];
ybs[3338]=['',2.2084963,-1.1287399,5.97];
ybs[3339]=['β Vol',2.2076961,-1.155549,3.77];
ybs[3340]=['',2.2371641,0.6491251,6.18];
ybs[3341]=['',2.2166164,-0.96139,6.53];
ybs[3342]=['',2.2174425,-0.9278277,5.09];
ybs[3343]=['',2.2434463,0.9257309,6.24];
ybs[3344]=['',2.2657102,1.3028439,6.31];
ybs[3345]=['',2.2268885,-0.4783132,6.7];
ybs[3346]=['2 UMa',2.253693,1.1356831,5.47];
ybs[3347]=['υ1 Cnc',2.2374183,0.4190075,5.75];
ybs[3348]=['',2.224678,-0.7720152,5.79];
ybs[3349]=['θ Cnc',2.2375987,0.3145199,5.35];
ybs[3350]=['',2.2242393,-0.8377893,5.33];
ybs[3351]=['',2.2261021,-0.781869,4.99];
ybs[3352]=['',2.2441074,0.6622135,5.9];
ybs[3353]=['',2.23872,0.1700047,6.83];
ybs[3354]=['',2.2311423,-0.5625665,5.65];
ybs[3355]=['',2.2273178,-0.8099176,5.99];
ybs[3356]=['',2.2310213,-0.6421826,6.69];
ybs[3357]=['32 Lyn',2.2459833,0.6346341,6.24];
ybs[3358]=['η Cnc',2.2425298,0.355469,5.33];
ybs[3359]=['',2.2361266,-0.342978,5.42];
ybs[3360]=['',2.2259913,-0.9645375,6.36];
ybs[3361]=['υ2 Cnc',2.2439306,0.4190598,6.36];
ybs[3362]=['',2.213615,-1.2246117,5.53];
ybs[3363]=['',2.2313148,-0.7820909,6.3];
ybs[3364]=['34 Cnc',2.2420211,0.1743917,6.46];
ybs[3365]=['',2.2349164,-0.6830831,6.31];
ybs[3366]=['',2.2408096,-0.263607,6.38];
ybs[3367]=['',2.2334265,-0.8367136,6.39];
ybs[3368]=['',2.2468562,0.2300792,6.28];
ybs[3369]=['33 Lyn',2.25194,0.6343286,5.78];
ybs[3370]=['',2.2464911,0.0817169,5.87];
ybs[3371]=['',2.2778876,1.2837321,6.15];
ybs[3372]=['',2.2487646,0.1462081,6.03];
ybs[3373]=['',2.2427907,-0.4307597,6.19];
ybs[3374]=['',2.234315,-0.9506419,6.34];
ybs[3375]=['',2.2476101,-0.0388583,5.81];
ybs[3376]=['',2.2415688,-0.5510887,6.38];
ybs[3377]=['',2.2419579,-0.6057716,6.36];
ybs[3378]=['',2.2369769,-0.9300174,5.69];
ybs[3379]=['35 Cnc',2.253907,0.3405957,6.58];
ybs[3380]=['',2.2433426,-0.6710008,6.49];
ybs[3381]=['',2.2446597,-0.6793418,5.96];
ybs[3382]=['',2.2436489,-0.8211,6.24];
ybs[3383]=['π1 UMa',2.2736511,1.1334838,5.64];
ybs[3384]=['',2.253807,0.0465705,6.33];
ybs[3385]=['',2.1947414,-1.4134441,5.69];
ybs[3386]=['',2.2572851,0.265953,6.32];
ybs[3387]=['',2.2558153,0.1142231,5.99];
ybs[3388]=['',2.2558373,0.1142667,7.25];
ybs[3389]=['',2.2488179,-0.5702604,6.43];
ybs[3390]=['3 Hya',2.2537415,-0.1406308,5.72];
ybs[3391]=['',2.2484295,-0.6577495,6.3];
ybs[3392]=['',2.2687824,0.9306931,5.66];
ybs[3393]=['',2.272875,1.0447976,6.48];
ybs[3394]=['',2.2531907,-0.4698237,5.96];
ybs[3395]=['π2 UMa',2.2779767,1.1213807,4.6];
ybs[3396]=['',2.2514712,-0.69892,6.47];
ybs[3397]=['',2.2714894,0.9223741,6.42];
ybs[3398]=['36 Cnc',2.2613414,0.1671946,5.88];
ybs[3399]=['',2.2487869,-0.8729982,5.01];
ybs[3400]=['',2.2727507,0.9186487,5.91];
ybs[3401]=['',2.2674249,0.5711661,5.94];
ybs[3402]=['δ Hya',2.2636712,0.0982164,4.16];
ybs[3403]=['',2.2624766,-0.0874365,6.19];
ybs[3404]=['37 Cnc',2.2656577,0.1657769,6.53];
ybs[3405]=['',2.2537042,-0.8909104,5.8];
ybs[3406]=['',2.2507381,-1.0137627,4.86];
ybs[3407]=['',2.2504099,-1.0175292,5.26];
ybs[3408]=['',2.2662969,-0.1176175,6.51];
ybs[3409]=['',2.2363902,-1.2816061,6.12];
ybs[3410]=['σ Hya',2.2684037,0.0569801,4.44];
ybs[3411]=['',2.2617272,-0.5903042,6.48];
ybs[3412]=['η Pyx',2.2636467,-0.4595673,5.27];
ybs[3413]=['',2.2607304,-0.702033,6.55];
ybs[3414]=['34 Lyn',2.2798141,0.7985972,5.37];
ybs[3415]=['',2.2760596,0.5561425,6.1];
ybs[3416]=['',2.2713861,0.1385841,6.45];
ybs[3417]=['',2.267363,-0.3458116,6.33];
ybs[3418]=['',2.261956,-0.7516315,4.14];
ybs[3419]=['39 Cnc',2.2747814,0.3478537,6.39];
ybs[3420]=['',2.2759121,0.3419566,6.44];
ybs[3421]=['ε Cnc',2.2762643,0.3397743,6.3];
ybs[3422]=['',2.2692842,-0.3968656,5.05];
ybs[3423]=['6 Hya',2.2734907,-0.2190813,4.98];
ybs[3424]=['',2.2588587,-1.0983278,5.47];
ybs[3425]=['ζ Pyx',2.2715549,-0.5172826,4.89];
ybs[3426]=['',2.2697886,-0.640253,6.13];
ybs[3427]=['',2.2661336,-0.9279412,6.47];
ybs[3428]=['',2.2885248,0.81721,6.22];
ybs[3429]=['',2.2779612,-0.1593398,6.63];
ybs[3430]=['β Pyx',2.2730541,-0.6175932,3.97];
ybs[3431]=['',2.2737859,-0.7040901,5.2];
ybs[3432]=['',2.2689402,-0.9340398,5.48];
ybs[3433]=['9 Hya',2.2807937,-0.279622,4.88];
ybs[3434]=['',2.2714327,-0.9273291,5.19];
ybs[3435]=['',2.267515,-1.0540731,6.36];
ybs[3436]=['',2.2746956,-0.790088,5.71];
ybs[3437]=['',2.2747791,-0.8155264,3.84];
ybs[3438]=['',2.2828384,-0.2102094,6.45];
ybs[3439]=['',2.2729,-0.9250092,3.62];
ybs[3440]=['',2.2728783,-0.9266382,5.61];
ybs[3441]=['γ Cnc',2.2886886,0.373328,4.66];
ybs[3442]=['45 Cnc',2.2880745,0.2199531,5.64];
ybs[3443]=['',2.2931029,0.6429652,6.33];
ybs[3444]=['',2.2773245,-0.8271903,4.77];
ybs[3445]=['',2.2766561,-0.8552115,5.9];
ybs[3446]=['η Hya',2.2879004,0.0579479,4.3];
ybs[3447]=['',2.2743958,-1.0057043,6.34];
ybs[3448]=['',2.2806255,-0.7939275,5.23];
ybs[3449]=['',2.2736858,-1.044377,4.33];
ybs[3450]=['',2.2912858,0.0742808,6.37];
ybs[3451]=['',2.289563,-0.1276223,4.62];
ybs[3452]=['θ Vol',2.265236,-1.2298206,5.2];
ybs[3453]=['δ Cnc',2.2946915,0.3154705,3.94];
ybs[3454]=['',2.2818575,-0.8408498,5.51];
ybs[3455]=['',2.2854604,-0.6286957,6.42];
ybs[3456]=['46 Cnc',2.2980364,0.5343929,6.13];
ybs[3457]=['49 Cnc',2.2947397,0.1745786,5.66];
ybs[3458]=['',2.2817234,-0.9281309,5.52];
ybs[3459]=['',2.2822032,-0.928374,4.86];
ybs[3460]=['α Pyx',2.288375,-0.5805825,3.68];
ybs[3461]=['10 Hya',2.2958021,0.0977629,6.13];
ybs[3462]=['',2.3158091,1.1628653,6.2];
ybs[3463]=['',2.2816495,-0.9748039,6.29];
ybs[3464]=['',2.2969932,-0.0467766,6.41];
ybs[3465]=['',2.2945952,-0.3708276,6.11];
ybs[3466]=['ι Cnc',2.303678,0.5006556,6.57];
ybs[3467]=['ι Cnc',2.3038086,0.5005633,4.02];
ybs[3468]=['',2.287887,-0.8709421,5.16];
ybs[3469]=['',2.2914847,-0.7457444,4.07];
ybs[3470]=['',2.300041,-0.0371482,5.7];
ybs[3471]=['',2.2937733,-0.6497207,5.76];
ybs[3472]=['',2.3001131,-0.1934863,6.25];
ybs[3473]=['50 Cnc',2.3043211,0.2099648,5.87];
ybs[3474]=['ε Hya',2.3034789,0.1106372,3.38];
ybs[3475]=['',2.2983949,-0.4444817,6.1];
ybs[3476]=['12 Hya',2.3011754,-0.2378436,4.32];
ybs[3477]=['δ Vel',2.292036,-0.9562179,1.96];
ybs[3478]=['',2.3053176,-0.0345093,5.29];
ybs[3479]=['',2.2984173,-0.8049657,3.91];
ybs[3480]=['',2.300283,-0.7191661,6.21];
ybs[3481]=['',2.2933517,-1.0263244,6.21];
ybs[3482]=['',2.3024173,-0.6056743,6.37];
ybs[3483]=['',2.2868228,-1.1918888,6.32];
ybs[3484]=['ρ Hya',2.3106905,0.1004837,4.36];
ybs[3485]=['',2.3088199,-0.1158714,6.09];
ybs[3486]=['',2.3005293,-0.8027193,5.46];
ybs[3487]=['',2.2898609,-1.1502477,6.05];
ybs[3488]=['',2.3040331,-0.8069621,5.75];
ybs[3489]=['',2.305838,-0.7298502,6.36];
ybs[3490]=['',2.300624,-0.9922096,4.49];
ybs[3491]=['',2.3207157,0.5795183,6.25];
ybs[3492]=['14 Hya',2.3144903,-0.0615034,5.31];
ybs[3493]=['',2.3078629,-0.7425362,6.43];
ybs[3494]=['η Cha',2.2713112,-1.3795196,5.47];
ybs[3495]=['',2.3066078,-0.9238162,6.3];
ybs[3496]=['',2.3211703,0.3272638,6.16];
ybs[3497]=['5 UMa',2.3349404,1.0800051,5.73];
ybs[3498]=['',2.3334125,1.029286,6.25];
ybs[3499]=['',2.3156778,-0.3687805,6.47];
ybs[3500]=['35 Lyn',2.3273187,0.7617451,5.15];
ybs[3501]=['',2.3284762,0.7894214,5.99];
ybs[3502]=['54 Cnc',2.3222607,0.2664954,6.38];
ybs[3503]=['',2.3281921,0.7316513,5.99];
ybs[3504]=['',2.3157589,-0.573542,5.21];
ybs[3505]=['',2.3166594,-0.5156419,5.87];
ybs[3506]=['',2.3145479,-0.7051382,5.48];
ybs[3507]=['',2.3181011,-0.5008961,6.17];
ybs[3508]=['',2.3155627,-0.6845643,6.39];
ybs[3509]=['γ Pyx',2.318882,-0.4850487,4.01];
ybs[3510]=['σ1 Cnc',2.3295781,0.5653481,5.66];
ybs[3511]=['',2.3149168,-0.7921872,4.93];
ybs[3512]=['53 Cnc',2.328992,0.4917833,6.23];
ybs[3513]=['ρ1 Cnc',2.3295176,0.4930333,5.95];
ybs[3514]=['15 Hya',2.3240482,-0.1266917,5.54];
ybs[3515]=['',2.2794969,-1.3813894,6.05];
ybs[3516]=['',2.3175074,-0.7360255,6];
ybs[3517]=['',2.327997,0.0917694,6.33];
ybs[3518]=['',2.3181845,-0.8135048,5.1];
ybs[3519]=['',2.335588,0.6188189,6.14];
ybs[3520]=['',2.3279696,-0.2324017,6.13];
ybs[3521]=['',2.3223585,-0.7432663,6.55];
ybs[3522]=['6 UMa',2.3494153,1.1260894,5.58];
ybs[3523]=['57 Cnc',2.3367823,0.532268,5.39];
ybs[3524]=['',2.3270223,-0.5688226,6.5];
ybs[3525]=['',2.3277739,-0.6392721,6.42];
ybs[3526]=['',2.3283658,-0.677297,5.82];
ybs[3527]=['',2.3219787,-1.0073203,5.59];
ybs[3528]=['',2.316311,-1.1671749,5.35];
ybs[3529]=['',2.3359796,-0.0962924,6];
ybs[3530]=['',2.3271892,-0.8454581,5.91];
ybs[3531]=['ρ2 Cnc',2.3428672,0.4859736,5.22];
ybs[3532]=['',2.3413162,0.2992934,6.64];
ybs[3533]=['',2.3271036,-0.9112571,6.39];
ybs[3534]=['',2.2910722,-1.3889952,5.79];
ybs[3535]=['',2.3117483,-1.2676611,6.11];
ybs[3536]=['',2.3487093,0.7949663,5.74];
ybs[3537]=['',2.347029,0.7001923,5.89];
ybs[3538]=['ζ Hya',2.3410576,0.1023185,3.11];
ybs[3539]=['',2.3328539,-0.7073817,6.47];
ybs[3540]=['',2.3284234,-0.990153,6.03];
ybs[3541]=['60 Cnc',2.3435367,0.2014593,5.41];
ybs[3542]=['',2.3324904,-0.8308344,5.33];
ybs[3543]=['17 Hya',2.3411198,-0.1405589,6.91];
ybs[3544]=['17 Hya',2.341127,-0.1405734,6.67];
ybs[3545]=['',2.3395906,-0.3198215,5.75];
ybs[3546]=['σ2 Cnc',2.3486298,0.572931,5.45];
ybs[3547]=['δ Pyx',2.3406869,-0.4845922,4.89];
ybs[3548]=['',2.3463424,0.0724849,6.14];
ybs[3549]=['',2.3489679,0.2977547,6.17];
ybs[3550]=['',2.3425875,-0.4171622,6.39];
ybs[3551]=['',2.3313345,-1.0548173,5.78];
ybs[3552]=['ο1 Cnc',2.3494035,0.2659696,5.2];
ybs[3553]=['',2.3390804,-0.7875746,6.26];
ybs[3554]=['61 Cnc',2.353046,0.5262078,6.29];
ybs[3555]=['',2.3455788,-0.293093,5.96];
ybs[3556]=['ο2 Cnc',2.3508864,0.270481,5.67];
ybs[3557]=['',2.3553437,0.6233999,6.51];
ybs[3558]=['',2.3512118,0.1623815,6.19];
ybs[3559]=['',2.3363097,-1.0179256,6.38];
ybs[3560]=['ι UMa',2.3591819,0.8370084,3.14];
ybs[3561]=['',2.3379002,-0.9607778,5.71];
ybs[3562]=['',2.3367052,-1.0598966,3.84];
ybs[3563]=['α Cnc',2.3547097,0.2054861,4.25];
ybs[3564]=['',2.3529112,0.0254384,6.59];
ybs[3565]=['',2.3429882,-0.9216561,4.69];
ybs[3566]=['σ3 Cnc',2.3599458,0.564333,5.2];
ybs[3567]=['ρ UMa',2.3756337,1.1788617,4.76];
ybs[3568]=['',2.3579077,0.3150348,6.38];
ybs[3569]=['',2.3550385,-0.2830424,5.86];
ybs[3570]=['',2.3651054,0.7277612,3.97];
ybs[3571]=['',2.3643725,0.6548365,6.44];
ybs[3572]=['',2.4412427,1.4676524,6.33];
ybs[3573]=['',2.3452814,-1.0352081,4.92];
ybs[3574]=['',2.3502815,-0.8492307,5.87];
ybs[3575]=['',2.3590254,-0.336722,6.18];
ybs[3576]=['',2.3569663,-0.5042369,6.25];
ybs[3577]=['',2.3695308,0.691636,6.36];
ybs[3578]=['66 Cnc',2.3680407,0.5614219,5.82];
ybs[3579]=['',2.3544744,-0.8258737,5.18];
ybs[3580]=['67 Cnc',2.3696944,0.4855024,6.07];
ybs[3581]=['',2.3677929,0.0969605,6.07];
ybs[3582]=['',2.3601005,-0.7214966,4.45];
ybs[3583]=['',2.380493,0.9459249,5.75];
ybs[3584]=['',2.3612336,-0.7549989,6.07];
ybs[3585]=['κ UMa',2.3783795,0.8215341,3.6];
ybs[3586]=['ν Cnc',2.3736113,0.4252828,5.45];
ybs[3587]=['',2.3695704,-0.0099194,5.67];
ybs[3588]=['',2.3654571,-0.4668606,6.2];
ybs[3589]=['',2.3559175,-1.0326786,5.16];
ybs[3590]=['',2.3731725,0.1258767,5.85];
ybs[3591]=['',2.3655629,-0.7321608,5.55];
ybs[3592]=['70 Cnc',2.3799399,0.48541,6.38];
ybs[3593]=['',2.3689928,-0.6891968,6.27];
ybs[3594]=['',2.3861954,0.8454969,5.95];
ybs[3595]=['',2.3616721,-1.0655043,5.79];
ybs[3596]=['',2.3667201,-0.9123487,5.23];
ybs[3597]=['',2.3833821,0.5635716,6.46];
ybs[3598]=['',2.3740376,-0.4466368,6.74];
ybs[3599]=['',2.3927215,1.0342309,6.45];
ybs[3600]=['σ1 UMa',2.4009336,1.1656237,5.14];
ybs[3601]=['',2.362211,-1.2002453,5.88];
ybs[3602]=['',2.3724976,-0.936118,6.4];
ybs[3603]=['',2.3906288,0.6695949,4.56];
ybs[3604]=['ω Hya',2.3871929,0.0873572,4.97];
ybs[3605]=['',2.377572,-0.8235174,3.75];
ybs[3606]=['α Vol',2.3683296,-1.1603247,4];
ybs[3607]=['σ2 UMa',2.4096416,1.1701736,4.8];
ybs[3608]=['',2.3941048,0.3995678,6.4];
ybs[3609]=['',2.391574,0.0240052,6.17];
ybs[3610]=['15 UMa',2.4014938,0.8991344,4.48];
ybs[3611]=['',2.3971118,0.5664074,6.5];
ybs[3612]=['τ Cnc',2.3967237,0.5160309,5.43];
ybs[3613]=['',2.3796448,-1.0112262,6.44];
ybs[3614]=['κ Cnc',2.3950706,0.1846628,5.24];
ybs[3615]=['τ UMa',2.4114456,1.1069704,4.67];
ybs[3616]=['',2.4005737,0.5898191,5.93];
ybs[3617]=['75 Cnc',2.4000571,0.4632299,5.98];
ybs[3618]=['ξ Cnc',2.4024063,0.3832274,5.14];
ybs[3619]=['κ Pyx',2.3954299,-0.4528442,4.58];
ybs[3620]=['',2.3875439,-0.9754728,6.11];
ybs[3621]=['19 Hya',2.3987629,-0.1514496,5.6];
ybs[3622]=['',2.3908359,-0.8953423,6.73];
ybs[3623]=['',2.384719,-1.12725,6.37];
ybs[3624]=['',2.4264635,1.2490591,6.55];
ybs[3625]=['λ Vel',2.3945461,-0.7595704,2.21];
ybs[3626]=['',2.403934,0.2002952,6.48];
ybs[3627]=['',2.400798,-0.2172223,5.77];
ybs[3628]=['',2.3983566,-0.4687212,6.15];
ybs[3629]=['',2.4001155,-0.3214322,5.73];
ybs[3630]=['',2.4082903,0.5388591,5.95];
ybs[3631]=['79 Cnc',2.4067219,0.3823632,6.01];
ybs[3632]=['20 Hya',2.4026349,-0.1549167,5.46];
ybs[3633]=['',2.3815198,-1.2326495,4.71];
ybs[3634]=['',2.3788399,-1.2686677,4.48];
ybs[3635]=['ε Pyx',2.4035472,-0.5315168,5.59];
ybs[3636]=['',2.4347119,1.2715676,5.96];
ybs[3637]=['',2.4057054,-0.4060547,6.53];
ybs[3638]=['',2.4019025,-0.864165,6.48];
ybs[3639]=['16 UMa',2.426109,1.0704676,5.13];
ybs[3640]=['',2.4131749,0.0938849,6.35];
ybs[3641]=['π1 Cnc',2.4150141,0.2601736,6.51];
ybs[3642]=['',2.4143936,0.0659386,6.14];
ybs[3643]=['36 Lyn',2.4225265,0.7527248,5.32];
ybs[3644]=['',2.412766,-0.3462191,5.73];
ybs[3645]=['',2.4079219,-0.7846445,5];
ybs[3646]=['21 Hya',2.4150785,-0.1256465,6.11];
ybs[3647]=['',2.410837,-0.68675,6];
ybs[3648]=['',2.4209757,0.3698981,6.48];
ybs[3649]=['',2.4099333,-0.8145943,5.79];
ybs[3650]=['',2.4065177,-1.0307152,3.44];
ybs[3651]=['17 UMa',2.432161,0.9887437,5.27];
ybs[3652]=['',2.4142598,-0.7627591,5.57];
ybs[3653]=['18 UMa',2.4335133,0.9412784,4.83];
ybs[3654]=['',2.4075306,-1.0891903,3.97];
ybs[3655]=['',2.428438,0.6028945,5.97];
ybs[3656]=['θ Hya',2.4237363,0.0388196,3.88];
ybs[3657]=['',2.4526128,1.2902234,6.5];
ybs[3658]=['',2.418501,-0.6755469,6.31];
ybs[3659]=['',2.4178162,-0.7393766,6.29];
ybs[3660]=['π2 Cnc',2.4278209,0.2592008,5.34];
ybs[3661]=['',2.4187298,-0.827779,5.92];
ybs[3662]=['',2.4213432,-0.772058,5.85];
ybs[3663]=['',2.4150504,-1.0385375,5.54];
ybs[3664]=['',2.4225728,-0.7560317,5.25];
ybs[3665]=['',2.427894,-0.2638071,6.35];
ybs[3666]=['',2.4388896,0.8155247,5.97];
ybs[3667]=['',2.4251882,-0.6578604,5.86];
ybs[3668]=['ζ Oct',2.3264412,-1.4922614,5.42];
ybs[3669]=['',2.4213815,-0.971443,5.27];
ybs[3670]=['',2.426112,-0.7966689,6.25];
ybs[3671]=['23 Hya',2.4336955,-0.1124657,5.24];
ybs[3672]=['',2.4280268,-0.6747504,4.94];
ybs[3673]=['24 Hya',2.4336085,-0.1542081,5.47];
ybs[3674]=['',2.428683,-0.6545637,4.62];
ybs[3675]=['β Car',2.4148312,-1.2183556,1.68];
ybs[3676]=['',2.4423891,0.6156261,5.75];
ybs[3677]=['',2.4353699,-0.2559439,5.84];
ybs[3678]=['',2.4297565,-0.785208,6.04];
ybs[3679]=['',2.4391913,0.1991411,6.41];
ybs[3680]=['38 Lyn',2.4442349,0.6407273,3.82];
ybs[3681]=['',2.4254897,-1.0206477,6.02];
ybs[3682]=['',2.4311628,-0.7741658,5.12];
ybs[3683]=['',2.4268348,-1.0065026,6.32];
ybs[3684]=['',2.4338556,-0.6892689,5.33];
ybs[3685]=['',2.4083422,-1.3395757,6.14];
ybs[3686]=['',2.4295499,-1.0058664,4.34];
ybs[3687]=['',2.4531128,0.8931535,6.13];
ybs[3688]=['',2.4577951,0.9879722,5.47];
ybs[3689]=['ι Car',2.4332776,-1.0361337,2.25];
ybs[3690]=['',2.4363506,-0.9527108,6.33];
ybs[3691]=['',2.4536269,0.6649022,6.12];
ybs[3692]=['',2.4460868,-0.1990703,6.62];
ybs[3693]=['',2.4382868,-0.8926013,5.26];
ybs[3694]=['',2.4459382,-0.2779639,5.78];
ybs[3695]=['α Lyn',2.4537818,0.5986521,3.13];
ybs[3696]=['26 Hya',2.4469945,-0.2106052,4.79];
ybs[3697]=['',2.4554681,0.5726347,6.16];
ybs[3698]=['',2.4409084,-0.9014963,5.87];
ybs[3699]=['27 Hya',2.450152,-0.168387,4.8];
ybs[3700]=['',2.4464794,-0.5968173,6.39];
ybs[3701]=['',2.4541167,0.2666654,6.53];
ybs[3702]=['',2.4329337,-1.2004422,5.39];
ybs[3703]=['',2.4357282,-1.1718466,6.11];
ybs[3704]=['',2.4519314,-0.2741903,6.33];
ybs[3705]=['',2.4506557,-0.5559336,6.82];
ybs[3706]=['',2.4493914,-0.6575248,6.05];
ybs[3707]=['',2.4443408,-0.9647886,6.28];
ybs[3708]=['θ Pyx',2.4541363,-0.4547962,4.72];
ybs[3709]=['',2.4873561,1.3090616,6.29];
ybs[3710]=['',2.4319683,-1.3087396,5.29];
ybs[3711]=['',2.4321704,-1.3059522,5.86];
ybs[3712]=['',2.4759453,1.1143398,6.28];
ybs[3713]=['',2.464288,0.4379028,6.41];
ybs[3714]=['',2.4604674,-0.1733409,6.53];
ybs[3715]=['',2.4714125,0.8985011,6.31];
ybs[3716]=['',2.4551289,-0.7380552,5.58];
ybs[3717]=['',2.4683274,0.6369331,6.67];
ybs[3718]=['',2.4498417,-1.0907753,4.81];
ybs[3719]=['',2.4585627,-0.6958179,6.54];
ybs[3720]=['',2.4573662,-0.8052971,5.75];
ybs[3721]=['κ Leo',2.4692209,0.4553352,4.46];
ybs[3722]=['',2.4543323,-0.9705324,5.63];
ybs[3723]=['λ Pyx',2.461527,-0.504868,4.69];
ybs[3724]=['κ Vel',2.4555929,-0.9617348,2.5];
ybs[3725]=['',2.463598,-0.6606124,6.48];
ybs[3726]=['',2.472835,0.2878369,6.29];
ybs[3727]=['',2.4658291,-0.6897381,6.06];
ybs[3728]=['28 Hya',2.4717134,-0.0909518,5.59];
ybs[3729]=['',2.4640194,-0.9046104,6.08];
ybs[3730]=['',2.4610389,-1.0540993,6.3];
ybs[3731]=['',2.4760366,-0.0271898,6.01];
ybs[3732]=['',2.4637005,-1.0776017,5.99];
ybs[3733]=['',2.4873857,0.7942405,5.41];
ybs[3734]=['29 Hya',2.4796639,-0.1626273,6.54];
ybs[3735]=['',2.4769957,-0.5040785,6.1];
ybs[3736]=['',2.4754144,-0.7085323,6.2];
ybs[3737]=['',2.4928588,0.971283,6.45];
ybs[3738]=['α Hya',2.4811753,-0.1527682,1.98];
ybs[3739]=['',2.4796199,-0.3916195,4.69];
ybs[3740]=['',2.4820781,-0.1076089,5.38];
ybs[3741]=['',2.5307048,1.4177097,4.29];
ybs[3742]=['',2.4696071,-1.0828697,5.77];
ybs[3743]=['',2.4740188,-0.9332808,5.11];
ybs[3744]=['ω Leo',2.4853782,0.1564167,5.41];
ybs[3745]=['3 Leo',2.485482,0.1412612,5.71];
ybs[3746]=['',2.4807028,-0.6126477,6.65];
ybs[3747]=['23 UMa',2.5010347,1.0989682,3.67];
ybs[3748]=['',2.4876699,-0.0235929,6.27];
ybs[3749]=['τ1 Hya',2.4881235,-0.0499819,4.6];
ybs[3750]=['',2.4892706,-0.0401465,6.14];
ybs[3751]=['',2.4748922,-1.1348778,6.05];
ybs[3752]=['',2.4897997,-0.0757858,6.26];
ybs[3753]=['',2.4879715,-0.3637874,5.66];
ybs[3754]=['7 LMi',2.4958825,0.5857354,5.85];
ybs[3755]=['ε Ant',2.4876824,-0.6291258,4.51];
ybs[3756]=['',2.4877128,-0.6719301,6.19];
ybs[3757]=['',2.4906186,-0.4091112,6.24];
ybs[3758]=['22 UMa',2.517071,1.2585352,5.72];
ybs[3759]=['8 LMi',2.4995029,0.6109944,5.37];
ybs[3760]=['',2.4908658,-0.4657377,5.48];
ybs[3761]=['24 UMa',2.5148031,1.2170813,4.56];
ybs[3762]=['',2.4931919,-0.2735361,5.85];
ybs[3763]=['λ Leo',2.4999378,0.399198,4.31];
ybs[3764]=['',2.5229058,1.2953936,6.46];
ybs[3765]=['θ UMa',2.5059401,0.9002606,3.17];
ybs[3766]=['',2.4841561,-1.088522,5.92];
ybs[3767]=['',2.4753831,-1.2513363,5.47];
ybs[3768]=['',2.5069666,0.861188,6.76];
ybs[3769]=['6 Leo',2.5006679,0.1679019,5.07];
ybs[3770]=['ζ1 Ant',2.4944166,-0.5582738,7];
ybs[3771]=['ζ1 Ant',2.4944676,-0.55824,6.18];
ybs[3772]=['ξ Leo',2.5006382,0.1955461,4.97];
ybs[3773]=['',2.4824194,-1.1658188,5.91];
ybs[3774]=['',2.490682,-0.900805,5.45];
ybs[3775]=['',2.4988594,-0.1858404,6.14];
ybs[3776]=['ψ Vel',2.4938768,-0.7079402,3.6];
ybs[3777]=['τ2 Hya',2.5005219,-0.0223535,4.57];
ybs[3778]=['',2.5000929,-0.1826712,6.13];
ybs[3779]=['ζ2 Ant',2.4978241,-0.5579388,5.93];
ybs[3780]=['',2.4977526,-0.6250127,5.87];
ybs[3781]=['9 LMi',2.5081123,0.6351371,6.18];
ybs[3782]=['',2.506993,0.493437,6.53];
ybs[3783]=['',2.4915443,-1.0202645,5.88];
ybs[3784]=['',2.5036792,0.0308606,6.11];
ybs[3785]=['ι Cha',2.4582265,-1.4116204,5.36];
ybs[3786]=['',2.5016707,-0.3402717,5.74];
ybs[3787]=['',2.5120986,0.8169133,6.52];
ybs[3788]=['',2.5012898,-0.5013265,6.46];
ybs[3789]=['26 UMa',2.5145313,0.9067804,4.5];
ybs[3790]=['10 LMi',2.5112398,0.6335721,4.55];
ybs[3791]=['',2.5049576,-0.1501221,6.12];
ybs[3792]=['',2.504386,-0.2375915,5.94];
ybs[3793]=['',2.4952856,-0.9971048,3.13];
ybs[3794]=['',2.5098261,0.407665,6.25];
ybs[3795]=['',2.506296,-0.1271678,6.24];
ybs[3796]=['',2.5304855,1.2737906,6.42];
ybs[3797]=['',2.5009559,-0.7111391,5.35];
ybs[3798]=['',2.5054194,-0.3702185,5.01];
ybs[3799]=['',2.5150117,0.689835,4.81];
ybs[3800]=['',2.5063721,-0.400729,5.91];
ybs[3801]=['',2.5163749,0.6958014,6.76];
ybs[3802]=['',2.5045442,-0.6846049,6.43];
ybs[3803]=['',2.4957254,-1.1661411,6.27];
ybs[3804]=['33 Hya',2.5116046,-0.1049214,5.56];
ybs[3805]=['11 LMi',2.5174814,0.6233155,5.41];
ybs[3806]=['',2.4992478,-1.097544,6.1];
ybs[3807]=['',2.5067973,-0.8569786,5.12];
ybs[3808]=['7 Leo',2.5178787,0.2492809,6.36];
ybs[3809]=['',2.508449,-0.8962556,5.01];
ybs[3810]=['',2.5219371,0.5421764,5.56];
ybs[3811]=['',2.4947522,-1.2771679,5.47];
ybs[3812]=['',2.5157317,-0.3434889,6.31];
ybs[3813]=['',2.5136871,-0.626933,6.49];
ybs[3814]=['',2.5360488,1.1724087,5.94];
ybs[3815]=['',2.5092156,-1.0354322,4.08];
ybs[3816]=['8 Leo',2.5229922,0.2851945,5.69];
ybs[3817]=['10 Leo',2.5235137,0.1176082,5];
ybs[3818]=['',2.5199704,-0.4328405,6.53];
ybs[3819]=['42 Lyn',2.5294053,0.7006094,5.25];
ybs[3820]=['',2.5218837,-0.4432082,5.7];
ybs[3821]=['',2.5185263,-0.8525668,6.17];
ybs[3822]=['34 Hya',2.5259887,-0.1661904,6.4];
ybs[3823]=['',2.522399,-0.5633215,5.63];
ybs[3824]=['',2.5288912,0.0794372,4.68];
ybs[3825]=['',2.5236171,-0.6316914,5.98];
ybs[3826]=['',2.5202668,-0.8631087,4.35];
ybs[3827]=['',2.5198313,-0.9257462,6.19];
ybs[3828]=['',2.5484286,1.2066947,5.69];
ybs[3829]=['27 UMa',2.5520629,1.2593125,5.17];
ybs[3830]=['',2.5216831,-0.9383925,5.45];
ybs[3831]=['',2.5158252,-1.1352929,6.56];
ybs[3832]=['',2.5257803,-0.7555302,5.5];
ybs[3833]=['',2.5650316,1.3619627,6.23];
ybs[3834]=['',2.5287697,-0.6931042,6.7];
ybs[3835]=['ι Hya',2.5348852,-0.0216586,3.91];
ybs[3836]=['37 Hya',2.534395,-0.1861991,6.31];
ybs[3837]=['',2.573268,1.3794407,6.17];
ybs[3838]=['',2.5367694,-0.1896732,6.37];
ybs[3839]=['κ Hya',2.5365679,-0.2518601,5.06];
ybs[3840]=['',2.5431751,0.5441823,5.89];
ybs[3841]=['43 Lyn',2.5452625,0.692179,5.62];
ybs[3842]=['ο Leo',2.5407589,0.1709315,3.52];
ybs[3843]=['13 Leo',2.5432713,0.4505403,6.24];
ybs[3844]=['',2.5487163,0.8435533,6.39];
ybs[3845]=['',2.5507626,0.9470928,6.47];
ybs[3846]=['',2.5304989,-1.0720859,4.52];
ybs[3847]=['13 LMi',2.5481901,0.6107655,6.14];
ybs[3848]=['',2.5406164,-0.4134728,4.77];
ybs[3849]=['',2.558132,1.1324435,6.17];
ybs[3850]=['ζ Cha',2.5010033,-1.4143705,5.11];
ybs[3851]=['15 Leo',2.5517187,0.5214198,5.64];
ybs[3852]=['',2.5447859,-0.4191306,4.94];
ybs[3853]=['',2.5367186,-1.0137216,5.32];
ybs[3854]=['',2.5382077,-1.0010891,5.8];
ybs[3855]=['28 UMa',2.5636927,1.1092144,6.34];
ybs[3856]=['ψ Leo',2.5521081,0.2429907,5.35];
ybs[3857]=['',2.5464334,-0.6213486,6.41];
ybs[3858]=['',2.5416701,-0.9653915,6];
ybs[3859]=['',2.5555657,0.3274946,6.5];
ybs[3860]=['',2.5657865,0.9953242,5.2];
ybs[3861]=['θ Ant',2.553255,-0.4864036,4.79];
ybs[3862]=['',2.5491668,-0.8958342,6.15];
ybs[3863]=['ε Leo',2.5615734,0.4131931,2.98];
ybs[3864]=['',2.5531751,-0.6923817,6.82];
ybs[3865]=['',2.5500907,-0.9423193,5.56];
ybs[3866]=['',2.5625766,0.1153417,5.79];
ybs[3867]=['18 Leo',2.5636484,0.2043765,5.63];
ybs[3868]=['',2.5582643,-0.5288791,6.45];
ybs[3869]=['',2.5634664,0.029417,5.65];
ybs[3870]=['19 Leo',2.5681941,0.2001534,6.45];
ybs[3871]=['',2.5741834,0.8014615,5.09];
ybs[3872]=['',2.5687438,0.197719,6.02];
ybs[3873]=['',2.558506,-0.9998184,6.46];
ybs[3874]=['',2.5562013,-1.0927063,3.69];
ybs[3875]=['',2.5833394,1.1430517,6.31];
ybs[3876]=['',2.562788,-0.7828688,5.55];
ybs[3877]=['',2.5594769,-1.0278952,6.22];
ybs[3878]=['υ UMa',2.5853331,1.0286476,3.8];
ybs[3879]=['20 Leo',2.5788813,0.3678869,6.09];
ybs[3880]=['υ Car',2.5640584,-1.1374686,3.01];
ybs[3881]=['υ Car',2.5641021,-1.1374784,6.26];
ybs[3882]=['',2.5759953,-0.6507864,5.97];
ybs[3883]=['4 Sex',2.5814413,0.0740432,6.24];
ybs[3884]=['φ UMa',2.5898764,0.9418269,4.59];
ybs[3885]=['',2.5716462,-0.9863313,6.06];
ybs[3886]=['23 Leo',2.5839341,0.2262769,6.46];
ybs[3887]=['',2.5777127,-0.6347701,6.37];
ybs[3888]=['',2.5778072,-0.7999512,5.08];
ybs[3889]=['6 Sex',2.5844669,-0.0758307,6.01];
ybs[3890]=['22 Leo',2.5878855,0.4240039,5.32];
ybs[3891]=['',2.5849819,-0.1096616,6.42];
ybs[3892]=['ν Cha',2.5582981,-1.3417396,5.45];
ybs[3893]=['υ1 Hya',2.5853247,-0.2608948,4.12];
ybs[3894]=['',2.5810292,-0.820928,5.73];
ybs[3895]=['μ Leo',2.5917598,0.4521285,3.88];
ybs[3896]=['7 Sex',2.5888288,0.041058,6.02];
ybs[3897]=['',2.5887686,-0.0004566,6.35];
ybs[3898]=['',2.5875377,-0.2903594,6.08];
ybs[3899]=['γ Sex',2.5899486,-0.1432356,5.05];
ybs[3900]=['',2.5838137,-0.8080059,5.62];
ybs[3901]=['',2.6031813,1.0648875,6.27];
ybs[3902]=['',2.5853293,-0.8141841,4.58];
ybs[3903]=['',2.5825283,-1.0389459,5.79];
ybs[3904]=['',2.5810541,-1.0968796,5.57];
ybs[3905]=['',2.5954912,0.1022099,5.95];
ybs[3906]=['',2.5915457,-0.4788159,6.3];
ybs[3907]=['31 UMa',2.6053804,0.8677304,5.27];
ybs[3908]=['',2.6191309,1.2701802,5.83];
ybs[3909]=['',2.5969871,-0.454392,4.88];
ybs[3910]=['',2.5906678,-0.9682252,6.48];
ybs[3911]=['',2.5984776,-0.3942816,6.24];
ybs[3912]=['',2.6123572,1.0003391,5.93];
ybs[3913]=['',2.6000446,-0.3335652,4.94];
ybs[3914]=['',2.5945541,-0.894465,5.93];
ybs[3915]=['',2.5968031,-0.7921377,5.71];
ybs[3916]=['',2.6074064,0.1541158,5.85];
ybs[3917]=['',2.5990499,-0.8787085,5.72];
ybs[3918]=['19 LMi',2.6136204,0.714753,5.14];
ybs[3919]=['',2.6149187,0.7908287,6.3];
ybs[3920]=['',2.6048138,-0.714319,6.41];
ybs[3921]=['',2.6081964,-0.4651864,6.28];
ybs[3922]=['',2.6072244,-0.5850605,5.84];
ybs[3923]=['',2.6087206,-0.4813264,6.32];
ybs[3924]=['',2.6691664,1.4627975,6.37];
ybs[3925]=['',2.6056015,-0.8977785,6.37];
ybs[3926]=['',2.6165259,0.4826791,6.3];
ybs[3927]=['ν Leo',2.6152814,0.2153977,5.26];
ybs[3928]=['',2.6147812,0.1433064,6.04];
ybs[3929]=['',2.623785,0.9897436,5.48];
ybs[3930]=['φ Vel',2.6076143,-0.9541839,3.54];
ybs[3931]=['',2.6091089,-0.9205199,6.12];
ybs[3932]=['',2.6216634,0.5155975,5.73];
ybs[3933]=['',2.6116068,-0.8467921,6.05];
ybs[3934]=['',2.6028574,-1.2477732,6.35];
ybs[3935]=['12 Sex',2.6216337,0.057264,6.7];
ybs[3936]=['',2.6184158,-0.4198186,6.21];
ybs[3937]=['η Ant',2.6171176,-0.6282243,5.23];
ybs[3938]=['',2.6085829,-1.127351,6.58];
ybs[3939]=['',2.6068898,-1.2078529,6.2];
ybs[3940]=['π Leo',2.6238786,0.1385844,4.7];
ybs[3941]=['20 LMi',2.6278576,0.5553555,5.36];
ybs[3942]=['',2.6354936,0.3812606,5.66];
ybs[3943]=['',2.6237067,-0.9957203,6.52];
ybs[3944]=['',2.6442666,0.9387541,5.74];
ybs[3945]=['',2.6287252,-0.9332039,6.2];
ybs[3946]=['',2.6345108,-0.5355022,6.54];
ybs[3947]=['',2.6297786,-1.0027613,6.2];
ybs[3948]=['',2.6466741,0.9122081,6.14];
ybs[3949]=['',2.6387021,-0.168924,6.12];
ybs[3950]=['',2.6297171,-1.0563623,5.94];
ybs[3951]=['13 Sex',2.6409215,0.0540397,6.45];
ybs[3952]=['',2.638423,-0.4436872,6.7];
ybs[3953]=['',2.6401356,-0.3177633,5.86];
ybs[3954]=['',2.6363064,-0.8157798,6.12];
ybs[3955]=['',2.6413273,-0.4256938,5.7];
ybs[3956]=['',2.6341049,-1.0521392,6.19];
ybs[3957]=['',2.6331866,-1.086657,6.42];
ybs[3958]=['',2.6411119,-0.6995408,6.43];
ybs[3959]=['',2.6478726,0.2731832,6.37];
ybs[3960]=['υ2 Hya',2.644931,-0.2298568,4.6];
ybs[3961]=['',2.6364851,-1.0819043,6.14];
ybs[3962]=['',2.6449864,-0.6368534,6.27];
ybs[3963]=['14 Sex',2.6525248,0.0960954,6.21];
ybs[3964]=['21 LMi',2.6559128,0.6132916,4.48];
ybs[3965]=['η Leo',2.6550943,0.2907215,3.52];
ybs[3966]=['',2.6487241,-0.8286011,5.08];
ybs[3967]=['',2.6537318,-0.3010217,5.6];
ybs[3968]=['',2.6482166,-0.9126916,6.52];
ybs[3969]=['',2.6594626,0.5497485,6.24];
ybs[3970]=['31 Leo',2.6574653,0.1726427,4.37];
ybs[3971]=['α Sex',2.6574322,-0.0083334,4.49];
ybs[3972]=['α Leo',2.6595356,0.2070188,1.35];
ybs[3973]=['μ1 Cha',2.6183113,-1.4367286,5.52];
ybs[3974]=['',2.6570857,-0.6534411,6.36];
ybs[3975]=['',2.6608445,-0.1918243,6.53];
ybs[3976]=['',2.6600177,-0.2743243,6.27];
ybs[3977]=['',2.6715342,0.7078153,6.32];
ybs[3978]=['',2.6659738,-0.2129671,6.24];
ybs[3979]=['17 Sex',2.6668379,-0.1486089,5.91];
ybs[3980]=['',2.6605674,-0.9061248,4.86];
ybs[3981]=['',2.6666451,-0.225539,5.31];
ybs[3982]=['',2.6636859,-0.62767,6.13];
ybs[3983]=['',2.6724541,0.6509263,5.85];
ybs[3984]=['λ Hya',2.6687921,-0.2174786,3.61];
ybs[3985]=['',2.6586354,-1.1505422,5.28];
ybs[3986]=['18 Sex',2.6703512,-0.1487868,5.65];
ybs[3987]=['μ2 Cha',2.6339281,-1.425419,6.6];
ybs[3988]=['34 Leo',2.6738003,0.2312264,6.44];
ybs[3989]=['',2.6618759,-1.0760874,5.6];
ybs[3990]=['',2.6719701,-0.1295606,6.25];
ybs[3991]=['',2.6683377,-0.7299217,5.98];
ybs[3992]=['',2.6618064,-1.2006024,5.81];
ybs[3993]=['',2.6748567,-0.5011393,6.28];
ybs[3994]=['19 Sex',2.6787491,0.0786751,5.77];
ybs[3995]=['',2.6775835,-0.3361597,6.44];
ybs[3996]=['',2.6836215,0.4717382,6.04];
ybs[3997]=['',2.6717866,-1.0286045,6.4];
ybs[3998]=['',2.6903956,1.0450683,6.25];
ybs[3999]=['',2.6726575,-1.0152099,5.72];
ybs[4000]=['',2.6756194,-0.9122866,6.16];
ybs[4001]=['',2.6804586,-0.473612,6.25];
ybs[4002]=['',2.6864048,0.3675734,6.02];
ybs[4003]=['',2.6807208,-0.5783854,6.38];
ybs[4004]=['22 LMi',2.6892707,0.5473447,6.46];
ybs[4005]=['',2.6820713,-0.7060381,5.9];
ybs[4006]=['',2.7043427,1.2734812,6.4];
ybs[4007]=['',2.6800182,-0.8960592,5.28];
ybs[4008]=['',2.6779984,-1.0476345,6.1];
ybs[4009]=['',2.6828517,-0.7054231,6.35];
ybs[4010]=['',2.6803752,-0.9051836,5.78];
ybs[4011]=['',2.7032546,1.2383524,6.66];
ybs[4012]=['',2.6793303,-1.078019,6.41];
ybs[4013]=['',2.6862671,-0.737041,3.85];
ybs[4014]=['23 LMi',2.6941714,0.5096847,5.35];
ybs[4015]=['',2.6796107,-1.1602973,5.16];
ybs[4016]=['32 UMa',2.7034064,1.1344661,5.82];
ybs[4017]=['24 LMi',2.6951537,0.4987222,6.49];
ybs[4018]=['',2.6940758,0.3077452,6.55];
ybs[4019]=['',2.6890912,-0.6392373,6.19];
ybs[4020]=['35 Leo',2.6953573,0.4083235,5.97];
ybs[4021]=['ζ Leo',2.6960168,0.4068248,3.44];
ybs[4022]=['',2.6960904,0.4409315,5.84];
ybs[4023]=['λ UMa',2.6982368,0.7471138,3.45];
ybs[4024]=['',2.6930957,-0.1974155,6.08];
ybs[4025]=['37 Leo',2.6957952,0.2377219,5.41];
ybs[4026]=['',2.6896834,-0.7543327,5.6];
ybs[4027]=['ω Car',2.6801655,-1.2242643,3.32];
ybs[4028]=['',2.688146,-0.9613566,6.16];
ybs[4029]=['39 Leo',2.6984235,0.4013927,5.82];
ybs[4030]=['',2.6955887,-0.362652,6.57];
ybs[4031]=['',2.7025529,0.4765983,6.52];
ybs[4032]=['ε Sex',2.6995992,-0.1427149,5.24];
ybs[4033]=['',2.6912631,-1.0473897,6.22];
ybs[4034]=['',2.7065726,0.8142385,6.43];
ybs[4035]=['',2.6944083,-0.8955778,6.3];
ybs[4036]=['',2.7086455,0.8427922,6];
ybs[4037]=['',2.7169031,1.1979697,5.96];
ybs[4038]=['',2.7061628,0.429408,6.4];
ybs[4039]=['',2.7013978,-0.5078929,5.34];
ybs[4040]=['',2.6956851,-1.0723326,3.4];
ybs[4041]=['',2.7123372,0.9367266,6.45];
ybs[4042]=['',2.713546,0.9443662,6];
ybs[4043]=['',2.7034232,-0.6442535,6.3];
ybs[4044]=['40 Leo',2.7092224,0.3379355,4.79];
ybs[4045]=['',2.7067294,-0.2205484,6];
ybs[4046]=['',2.7025975,-0.7291388,5.96];
ybs[4047]=['γ1 Leo',2.7102608,0.3444069,2.61];
ybs[4048]=['γ2 Leo',2.7102825,0.3443875,3.8];
ybs[4049]=['',2.7079671,-0.0910073,6.37];
ybs[4050]=['',2.7098839,-0.1600028,6.32];
ybs[4051]=['',2.702747,-0.9811938,5.81];
ybs[4052]=['',2.7601285,1.4685447,5.5];
ybs[4053]=['',2.7071167,-0.9623384,4.57];
ybs[4054]=['23 Sex',2.7146111,0.0380638,6.66];
ybs[4055]=['',2.7041736,-1.130707,5.67];
ybs[4056]=['',2.7103173,-0.8344036,5.65];
ybs[4057]=['',2.7203105,0.7176855,5.76];
ybs[4058]=['',2.7147159,-0.3157972,6.51];
ybs[4059]=['μ UMa',2.7209856,0.7223973,3.05];
ybs[4060]=['42 Leo',2.7183196,0.2594702,6.12];
ybs[4061]=['',2.7161374,-0.415733,6.5];
ybs[4062]=['',2.7299453,1.1424376,4.97];
ybs[4063]=['',2.7166898,-0.395095,6.51];
ybs[4064]=['',2.712756,-0.9800342,4.5];
ybs[4065]=['27 LMi',2.7241881,0.5898998,5.9];
ybs[4066]=['',2.7194214,-0.3486473,6.13];
ybs[4067]=['43 Leo',2.7232845,0.1122813,6.07];
ybs[4068]=['',2.7266787,0.5149842,6.39];
ybs[4069]=['',2.7243043,0.0974742,6.54];
ybs[4070]=['',2.7194542,-0.7288336,4.83];
ybs[4071]=['28 LMi',2.7287088,0.5865895,5.5];
ybs[4072]=['25 Sex',2.7250241,-0.073016,5.97];
ybs[4073]=['',2.7236148,-0.5283375,6.27];
ybs[4074]=['',2.7643266,1.4389816,5.26];
ybs[4075]=['',2.7285085,0.0394191,6.32];
ybs[4076]=['',2.7246235,-0.6653079,5.33];
ybs[4077]=['',2.7253321,-0.7341328,6.27];
ybs[4078]=['44 Leo',2.7331185,0.1514071,5.61];
ybs[4079]=['',2.7210087,-1.1695602,4.99];
ybs[4080]=['30 LMi',2.736431,0.5879357,4.74];
ybs[4081]=['',2.725503,-1.0133956,6.35];
ybs[4082]=['',2.7350029,-0.1251322,5.57];
ybs[4083]=['',2.73237,-0.7431223,6.18];
ybs[4084]=['μ Hya',2.7363897,-0.2957684,3.81];
ybs[4085]=['',2.7304665,-1.0242645,5.95];
ybs[4086]=['',2.7433705,0.7241481,6.02];
ybs[4087]=['',2.7409353,0.3360517,6.15];
ybs[4088]=['',2.7461746,0.8495284,6.44];
ybs[4089]=['',2.7361728,-0.7478524,6.13];
ybs[4090]=['β LMi',2.7450702,0.6387371,4.21];
ybs[4091]=['45 Leo',2.7435872,0.168464,6.04];
ybs[4092]=['',2.7262827,-1.294007,4];
ybs[4093]=['',2.7484396,0.7871746,6.35];
ybs[4094]=['α Ant',2.7407788,-0.5441568,4.25];
ybs[4095]=['',2.7278025,-1.2929611,6.19];
ybs[4096]=['35 UMa',2.755035,1.1434593,6.32];
ybs[4097]=['',2.7386398,-0.9597134,5.58];
ybs[4098]=['',2.7572694,1.1195707,6.12];
ybs[4099]=['',2.7481213,-0.0672464,6.05];
ybs[4100]=['',2.7410926,-1.0079107,4.66];
ybs[4101]=['',2.7441634,-0.8642143,6.1];
ybs[4102]=['36 UMa',2.7576016,0.9751104,4.84];
ybs[4103]=['32 LMi',2.7548065,0.6774417,5.77];
ybs[4104]=['',2.7431005,-1.0271207,3.82];
ybs[4105]=['',2.7406343,-1.1486859,6.01];
ybs[4106]=['δ Sex',2.7513865,-0.0497377,5.21];
ybs[4107]=['',2.7509908,-0.5196577,5.58];
ybs[4108]=['δ Ant',2.7514391,-0.5361272,5.56];
ybs[4109]=['β Sex',2.7549647,-0.0130498,5.09];
ybs[4110]=['',2.7471311,-1.121944,5.29];
ybs[4111]=['',2.784395,1.4029392,6.52];
ybs[4112]=['',2.7578642,-0.1352348,6.2];
ybs[4113]=['',2.7578596,-0.2390966,5.58];
ybs[4114]=['33 LMi',2.7622956,0.5631893,5.9];
ybs[4115]=['',2.757046,-0.4641659,6.51];
ybs[4116]=['',2.778847,1.3194918,4.84];
ybs[4117]=['46 Leo',2.7634811,0.2448015,5.46];
ybs[4118]=['',2.7551014,-1.0727999,6.43];
ybs[4119]=['',2.7524568,-1.1710405,6.19];
ybs[4120]=['',2.7611795,-0.4947753,6.05];
ybs[4121]=['',2.7709903,0.9317625,6.45];
ybs[4122]=['',2.7684544,0.7036157,4.75];
ybs[4123]=['ρ Leo',2.7660955,0.1604904,3.85];
ybs[4124]=['',2.7585797,-0.9394496,4.89];
ybs[4125]=['',2.7614641,-0.7885001,5.74];
ybs[4126]=['',2.7613985,-0.7885486,6.09];
ybs[4127]=['34 LMi',2.7695728,0.6087223,5.58];
ybs[4128]=['',2.7527115,-1.2584481,4.74];
ybs[4129]=['',2.7641066,-0.7806869,5.91];
ybs[4130]=['',2.7610806,-1.0785496,3.32];
ybs[4131]=['37 UMa',2.7773797,0.9943326,5.16];
ybs[4132]=['',2.7555855,-1.2798937,4.93];
ybs[4133]=['',2.765738,-0.8223046,5.02];
ybs[4134]=['',2.7646183,-1.0258723,6];
ybs[4135]=['44 Hya',2.7708846,-0.4163788,5.08];
ybs[4136]=['48 Leo',2.7747368,0.1194152,5.08];
ybs[4137]=['',2.7673921,-1.0175551,6.14];
ybs[4138]=['49 Leo',2.7757919,0.1490268,5.67];
ybs[4139]=['',2.7750277,-0.406448,6.1];
ybs[4140]=['35 LMi',2.7819812,0.6320713,6.28];
ybs[4141]=['',2.7707187,-1.0663833,6.23];
ybs[4142]=['',2.7780954,-0.326044,6.49];
ybs[4143]=['',2.7758345,-0.69245,5.38];
ybs[4144]=['',2.775563,-0.7640424,6.08];
ybs[4145]=['',2.7810077,-0.186667,6.57];
ybs[4146]=['φ2 Hya',2.7808778,-0.2872174,6.03];
ybs[4147]=['',2.7798415,-0.4675188,6.29];
ybs[4148]=['',2.7820767,-0.2154124,5.7];
ybs[4149]=['',2.7769051,-1.0065231,4.45];
ybs[4150]=['',2.7849351,-0.2070078,6.52];
ybs[4151]=['',2.7562076,-1.4317293,7.07];
ybs[4152]=['',2.7848597,-0.4803944,4.89];
ybs[4153]=['',2.7864846,-0.2355597,4.82];
ybs[4154]=['',2.7800916,-1.0415533,5.08];
ybs[4155]=['',2.7943585,0.9347268,5.52];
ybs[4156]=['37 LMi',2.7921985,0.5561274,4.71];
ybs[4157]=['',2.7847382,-0.8436556,3.84];
ybs[4158]=['38 LMi',2.7940813,0.659692,5.85];
ybs[4159]=['',2.7849844,-1.0270463,5.45];
ybs[4160]=['',2.7742217,-1.3337953,6.3];
ybs[4161]=['φ3 Hya',2.7909336,-0.2965137,4.91];
ybs[4162]=['',2.7921138,-0.2191431,6.04];
ybs[4163]=['',2.7876561,-1.0012708,5.91];
ybs[4164]=['γ Cha',2.77381,-1.3739135,4.11];
ybs[4165]=['',2.7915989,-0.7481523,6.11];
ybs[4166]=['',2.8069588,1.1925905,5.75];
ybs[4167]=['',2.7906755,-1.0348997,4.66];
ybs[4168]=['38 UMa',2.8073432,1.1449961,5.12];
ybs[4169]=['',2.7917362,-1.0285105,5.92];
ybs[4170]=['',2.7932605,-0.9724235,4.28];
ybs[4171]=['',2.8125006,1.2036308,5];
ybs[4172]=['33 Sex',2.8034362,-0.0323669,6.26];
ybs[4173]=['',2.8005805,-0.6257771,6.37];
ybs[4174]=['',2.8073088,0.5512443,6.02];
ybs[4175]=['',2.7966238,-1.1381839,5.52];
ybs[4176]=['',2.791579,-1.3021204,6.07];
ybs[4177]=['39 UMa',2.8146215,0.9963374,5.8];
ybs[4178]=['',2.8017788,-1.0435276,6.42];
ybs[4179]=['40 LMi',2.8108931,0.4574935,5.51];
ybs[4180]=['',2.8081657,-0.2458823,6.24];
ybs[4181]=['',2.8135237,0.8044342,5.18];
ybs[4182]=['41 LMi',2.8125334,0.4027375,5.08];
ybs[4183]=['35 Sex',2.8120007,0.0808894,5.79];
ybs[4184]=['',2.8087609,-0.572972,5.64];
ybs[4185]=['',2.8210354,1.1745701,6];
ybs[4186]=['',2.8056273,-1.1271218,4.82];
ybs[4187]=['',2.8160837,0.3428751,6.27];
ybs[4188]=['',2.8078567,-1.0354839,5.38];
ybs[4189]=['θ Car',2.8088011,-1.1258684,2.76];
ybs[4190]=['',2.8115503,-1.0590629,4.57];
ybs[4191]=['36 Sex',2.8198633,0.0414445,6.28];
ybs[4192]=['41 UMa',2.8261699,0.9992386,6.34];
ybs[4193]=['42 LMi',2.8233158,0.5335233,5.24];
ybs[4194]=['',2.812739,-1.1233307,5.77];
ybs[4195]=['',2.8139034,-1.1183088,4.82];
ybs[4196]=['',2.8015357,-1.394451,5.97];
ybs[4197]=['',2.8240021,0.1092478,6.37];
ybs[4198]=['51 Leo',2.8255197,0.3277329,5.49];
ybs[4199]=['52 Leo',2.8255195,0.2457606,5.48];
ybs[4200]=['η Car',2.8182693,-1.0436648,6.21];
ybs[4201]=['',2.8142492,-1.2387175,6.26];
ybs[4202]=['',2.8151787,-1.2386308,6.46];
ybs[4203]=['',2.8145785,-1.2663619,6.27];
ybs[4204]=['',2.827093,-0.3038689,5.42];
ybs[4205]=['',2.837238,1.1347806,6.39];
ybs[4206]=['μ Vel',2.8261116,-0.8645264,2.69];
ybs[4207]=['',2.823568,-1.0597108,6.25];
ybs[4208]=['',2.8304651,-0.2683585,6.67];
ybs[4209]=['',2.8233327,-1.1279823,5.34];
ybs[4210]=['',2.8243121,-1.1235905,5.23];
ybs[4211]=['',2.8267013,-0.9925855,5.23];
ybs[4212]=['',2.8258747,-1.1256859,4.85];
ybs[4213]=['43 LMi',2.8367568,0.5114119,6.15];
ybs[4214]=['',2.8351716,-0.0361794,5.93];
ybs[4215]=['',2.8328787,-0.5550547,5.88];
ybs[4216]=['',2.8296793,-1.004989,6.36];
ybs[4217]=['53 Leo',2.8378413,0.1820579,5.34];
ybs[4218]=['',2.8315248,-1.047775,6];
ybs[4219]=['40 Sex',2.837818,-0.0722269,6.61];
ybs[4220]=['44 LMi',2.8408406,0.4862428,6.04];
ybs[4221]=['δ1 Cha',2.8162873,-1.4064409,5.47];
ybs[4222]=['ν Hya',2.8391492,-0.2846246,3.11];
ybs[4223]=['',2.8396588,-0.1739565,5.86];
ybs[4224]=['δ2 Cha',2.8185453,-1.4076738,4.45];
ybs[4225]=['43 UMa',2.8470644,0.9855487,5.67];
ybs[4226]=['42 UMa',2.8480812,1.0333313,5.58];
ybs[4227]=['41 Sex',2.8421865,-0.1572902,5.79];
ybs[4228]=['',2.8403271,-0.5964188,5.61];
ybs[4229]=['',2.8372935,-1.0373842,5.91];
ybs[4230]=['',2.8456904,-0.0559711,5.95];
ybs[4231]=['',2.8527177,0.9154363,6.65];
ybs[4232]=['',2.8527957,0.91436,6.44];
ybs[4233]=['',2.8578972,1.2171767,5.93];
ybs[4234]=['',2.8507013,0.0158946,6.38];
ybs[4235]=['',2.8523173,-0.0055158,6.31];
ybs[4236]=['44 UMa',2.8573969,0.9506844,5.1];
ybs[4237]=['46 LMi',2.8558211,0.5951616,3.83];
ybs[4238]=['ω UMa',2.8588779,0.7518031,4.71];
ybs[4239]=['',2.8558449,-0.041365,6.12];
ybs[4240]=['',2.8510036,-1.0010365,5.25];
ybs[4241]=['',2.8559852,-0.353493,5.24];
ybs[4242]=['',2.8562863,-0.2715791,6.38];
ybs[4243]=['',2.8572135,-0.0391648,5.45];
ybs[4244]=['48 LMi',2.8617671,0.4428926,6.2];
ybs[4245]=['',2.8595718,-0.2421286,5.66];
ybs[4246]=['',2.8630445,0.5920109,5.72];
ybs[4247]=['',2.8552257,-1.0291873,3.78];
ybs[4248]=['46 UMa',2.8663921,0.5827975,5.03];
ybs[4249]=['54 Leo',2.8657249,0.4299554,4.5];
ybs[4250]=['54 Leo',2.8657612,0.4299409,6.3];
ybs[4251]=['',2.8634095,-0.3626798,6.44];
ybs[4252]=['',2.8554121,-1.2363049,5.99];
ybs[4253]=['',2.8623294,-0.7394279,6.11];
ybs[4254]=['',2.8687134,0.7311734,6.03];
ybs[4255]=['55 Leo',2.865875,0.0108533,5.91];
ybs[4256]=['',2.8594785,-1.0810843,5.93];
ybs[4257]=['56 Leo',2.8673185,0.1059439,5.81];
ybs[4258]=['',2.8483912,-1.3905738,6.33];
ybs[4259]=['',2.8686101,0.3880999,6.14];
ybs[4260]=['50 LMi',2.8699188,0.443048,6.35];
ybs[4261]=['',2.8630104,-1.0582274,5.92];
ybs[4262]=['',2.8868196,1.3553229,6.2];
ybs[4263]=['ι Ant',2.8698535,-0.6501877,4.6];
ybs[4264]=['',2.8713965,-0.8880285,5.91];
ybs[4265]=['',2.8822283,0.9034979,6.17];
ybs[4266]=['',2.874076,-1.0445328,6.11];
ybs[4267]=['47 UMa',2.8827335,0.7036234,5.05];
ybs[4268]=['',2.8830157,0.6279245,6];
ybs[4269]=['',2.8705321,-1.3127495,6.13];
ybs[4270]=['',2.8862123,0.7925606,5.47];
ybs[4271]=['',2.8833379,0.2022868,6.55];
ybs[4272]=['',2.8808855,-0.5908429,5.71];
ybs[4273]=['',2.8871298,0.8968581,6.43];
ybs[4274]=['',2.8823261,-0.2874472,5.89];
ybs[4275]=['',2.8865932,0.7469249,6.02];
ybs[4276]=['',2.8904334,1.1048852,6.39];
ybs[4277]=['α Crt',2.8834411,-0.3213945,4.08];
ybs[4278]=['49 UMa',2.8886965,0.6823611,5.08];
ybs[4279]=['',2.8853103,-0.2478201,5.88];
ybs[4280]=['',2.8802601,-1.0722578,6.16];
ybs[4281]=['58 Leo',2.8870778,0.0611168,4.84];
ybs[4282]=['',2.8840421,-0.7665993,5.81];
ybs[4283]=['',2.8847912,-0.7389992,4.39];
ybs[4284]=['59 Leo',2.8879153,0.1044683,4.99];
ybs[4285]=['β UMa',2.8934158,0.9820367,2.37];
ybs[4286]=['',2.8845447,-0.9064102,6.15];
ybs[4287]=['',2.8886113,-0.2776574,6.34];
ybs[4288]=['',2.8872384,-0.5577238,6.07];
ybs[4289]=['61 Leo',2.8925559,-0.04539,4.74];
ybs[4290]=['60 Leo',2.8949566,0.350178,4.42];
ybs[4291]=['α UMa',2.9018013,1.0757275,1.79];
ybs[4292]=['',2.8948413,-0.4703208,6.23];
ybs[4293]=['',2.8987459,-0.0151602,6.14];
ybs[4294]=['',2.8774691,-1.4254391,6.71];
ybs[4295]=['',2.898673,-0.1993118,5.5];
ybs[4296]=['62 Leo',2.9003525,-0.0020419,5.95];
ybs[4297]=['',2.8985449,-0.5598484,6.46];
ybs[4298]=['',2.9002247,-0.2365027,6.34];
ybs[4299]=['51 UMa',2.9047133,0.6654086,6];
ybs[4300]=['χ Leo',2.9065534,0.1260088,4.63];
ybs[4301]=['',2.9037773,-0.8341876,5.67];
ybs[4302]=['η Oct',2.8753823,-1.4784577,6.19];
ybs[4303]=['',2.9056399,-0.6269404,5.43];
ybs[4304]=['χ1 Hya',2.9076122,-0.4783945,4.94];
ybs[4305]=['',2.9087946,-0.1955692,6.09];
ybs[4306]=['',2.9061488,-0.8640922,6.13];
ybs[4307]=['χ2 Hya',2.9103568,-0.478294,5.71];
ybs[4308]=['',2.9106014,-0.8958593,6.3];
ybs[4309]=['65 Leo',2.9147385,0.0320965,5.52];
ybs[4310]=['',2.9133333,-0.5034281,6.77];
ybs[4311]=['',2.9121816,-0.891395,6.32];
ybs[4312]=['64 Leo',2.9182328,0.4050378,6.46];
ybs[4313]=['',2.9121228,-1.0261102,6.02];
ybs[4314]=['',2.9154525,-0.5707891,6.59];
ybs[4315]=['',2.9122181,-1.0915407,4.61];
ybs[4316]=['',2.911516,-1.1336997,6.41];
ybs[4317]=['',2.9159164,-0.7462193,5.15];
ybs[4318]=['',2.918829,-0.5286847,6.54];
ybs[4319]=['',2.9130303,-1.2390893,5.57];
ybs[4320]=['',2.9278389,1.1710004,6.06];
ybs[4321]=['',2.9203893,-0.5251608,6.49];
ybs[4322]=['67 Leo',2.9232881,0.4283307,5.68];
ybs[4323]=['',2.9255924,0.63168,5.74];
ybs[4324]=['',2.922455,-0.4921363,5.44];
ybs[4325]=['ψ UMa',2.9272067,0.7746072,3.01];
ybs[4326]=['',2.9270929,0.7520731,5.89];
ybs[4327]=['',2.9213147,-1.0313456,3.91];
ybs[4328]=['',2.9211135,-1.0832206,5.13];
ybs[4329]=['',2.9274697,-0.5669598,5.81];
ybs[4330]=['',2.9388287,1.1895251,6.4];
ybs[4331]=['',2.9358863,0.2492883,6.3];
ybs[4332]=['',2.9314943,-1.0222794,6.88];
ybs[4333]=['β Crt',2.9352878,-0.4004297,4.48];
ybs[4334]=['',2.9407967,0.9560379,6.63];
ybs[4335]=['',2.9396084,0.6230199,6.41];
ybs[4336]=['',2.9377652,-0.5681231,6.38];
ybs[4337]=['ψ Crt',2.9390293,-0.3249313,6.13];
ybs[4338]=['',2.9393075,-0.3816401,6.4];
ybs[4339]=['',2.9334426,-1.2488434,6.35];
ybs[4340]=['',2.9388766,-0.8590215,5.36];
ybs[4341]=['',2.9446084,0.7150839,6.33];
ybs[4342]=['',2.938832,-1.0547844,4.6];
ybs[4343]=['',2.9406004,-0.8701098,6.11];
ybs[4344]=['',2.9419831,-0.7764881,5.8];
ybs[4345]=['',2.9393709,-1.1220186,5.23];
ybs[4346]=['69 Leo',2.9446398,-0.0032646,5.42];
ybs[4347]=['δ Leo',2.9463138,0.3561562,2.56];
ybs[4348]=['',2.9458766,0.138635,5.79];
ybs[4349]=['θ Leo',2.9468476,0.267246,3.34];
ybs[4350]=['',2.9436256,-0.9311152,5.76];
ybs[4351]=['',2.942852,-1.0426027,5.74];
ybs[4352]=['72 Leo',2.9511088,0.4010432,4.63];
ybs[4353]=['',2.955221,0.9190117,6.5];
ybs[4354]=['',2.9492276,-0.7653549,6.21];
ybs[4355]=['73 Leo',2.9539184,0.2302082,5.32];
ybs[4356]=['',2.9543439,0.2221311,6.67];
ybs[4357]=['',2.9579134,0.861473,5.88];
ybs[4358]=['φ Leo',2.9572764,-0.0657864,4.47];
ybs[4359]=['',2.9585984,-0.1265777,6.14];
ybs[4360]=['',2.9560259,-0.8028095,6.31];
ybs[4361]=['75 Leo',2.9600562,0.0330369,5.18];
ybs[4362]=['',2.9593357,-0.6655308,6.27];
ybs[4363]=['',2.9613651,-0.6083333,6.45];
ybs[4364]=['ξ UMa',2.9641627,0.5482324,4.87];
ybs[4365]=['ξ UMa',2.9641699,0.5482324,4.41];
ybs[4366]=['',2.9616238,-0.6397008,6.68];
ybs[4367]=['ν UMa',2.9654697,0.5755464,3.48];
ybs[4368]=['',2.964747,0.2071173,6.66];
ybs[4369]=['',2.9592175,-1.185799,6.06];
ybs[4370]=['55 UMa',2.9683657,0.6644068,4.78];
ybs[4371]=['76 Leo',2.9671511,0.0267511,5.91];
ybs[4372]=['δ Crt',2.9688989,-0.2599926,3.56];
ybs[4373]=['',2.9766039,1.1690658,6.21];
ybs[4374]=['',2.967934,-1.1292342,5.99];
ybs[4375]=['',2.9635787,-1.3925352,6.35];
ybs[4376]=['σ Leo',2.9768638,0.1031737,4.05];
ybs[4377]=['',2.9687711,-1.3135414,6.27];
ybs[4378]=['',2.9803416,0.9940856,6.43];
ybs[4379]=['',2.9710449,-1.2585982,6.41];
ybs[4380]=['π Cen',2.9757696,-0.953109,3.89];
ybs[4381]=['',2.9850129,1.1207174,6.02];
ybs[4382]=['56 UMa',2.9845158,0.7568551,4.99];
ybs[4383]=['',2.9819411,-0.7812786,6.12];
ybs[4384]=['',2.9862672,0.0002348,6.05];
ybs[4385]=['λ Crt',2.9864386,-0.329836,5.09];
ybs[4386]=['',2.9856402,-0.6332565,5];
ybs[4387]=['',2.9788198,-1.3565818,6.43];
ybs[4388]=['',2.9850415,-0.9930511,5.79];
ybs[4389]=['ι Leo',2.9890535,0.1817045,3.94];
ybs[4390]=['79 Leo',2.989496,0.0225062,5.39];
ybs[4391]=['',2.9858386,-1.1357417,5.11];
ybs[4392]=['ε Crt',2.9919216,-0.1915981,4.83];
ybs[4393]=['',2.9906338,-0.746782,6.12];
ybs[4394]=['',2.9936673,0.1974305,5.8];
ybs[4395]=['γ Crt',2.9930676,-0.3107075,4.08];
ybs[4396]=['',2.9891208,-1.263181,5.59];
ybs[4397]=['',2.9982838,0.9727092,5.75];
ybs[4398]=['81 Leo',2.9964219,0.2851518,5.57];
ybs[4399]=['',2.9955972,-0.6314853,5.22];
ybs[4400]=['80 Leo',2.9973408,0.065303,6.37];
ybs[4401]=['',2.9958524,-0.6608892,5.89];
ybs[4402]=['',3.0000998,0.5817548,6.32];
ybs[4403]=['',2.9961957,-1.118602,5.17];
ybs[4404]=['83 Leo',3.0013574,0.0505198,6.5];
ybs[4405]=['',3.0000714,-1.0687304,5.3];
ybs[4406]=['κ Crt',3.0030373,-0.2177329,5.94];
ybs[4407]=['',3.00109,-0.9298849,5.81];
ybs[4408]=['τ Leo',3.0065121,0.0477791,4.95];
ybs[4409]=['',3.0063075,-0.0317399,6.25];
ybs[4410]=['',3.006462,-0.61867,6.45];
ybs[4411]=['',3.0119844,1.0761645,5.83];
ybs[4412]=['57 UMa',3.0116688,0.6844884,5.31];
ybs[4413]=['',3.0090676,-0.7468748,5.08];
ybs[4414]=['',3.0147114,0.9881846,6.28];
ybs[4415]=['',3.0071955,-1.2669874,6.09];
ybs[4416]=['85 Leo',3.014259,0.2669419,5.74];
ybs[4417]=['',3.0168085,0.9467179,6.41];
ybs[4418]=['',3.0138211,-0.4290371,5.76];
ybs[4419]=['',3.0251165,1.4138629,6.15];
ybs[4420]=['',3.0175979,0.8122546,6.35];
ybs[4421]=['58 UMa',3.0180087,0.7514443,5.94];
ybs[4422]=['87 Leo',3.0168613,-0.0544951,4.77];
ybs[4423]=['86 Leo',3.0176973,0.3192379,5.52];
ybs[4424]=['λ Dra',3.0222963,1.2079826,3.84];
ybs[4425]=['',3.0196349,0.8344489,6.42];
ybs[4426]=['',3.0208987,0.8494584,6.56];
ybs[4427]=['88 Leo',3.0231938,0.248633,6.2];
ybs[4428]=['',3.0204852,-1.0715819,6.38];
ybs[4429]=['',3.0261845,1.0640161,5.48];
ybs[4430]=['',3.0232302,-0.3646951,6.24];
ybs[4431]=['ο1 Cen',3.0227888,-1.0395363,5.13];
ybs[4432]=['ο2 Cen',3.0229845,-1.0408211,5.15];
ybs[4433]=['',3.0252596,-0.5128159,5.81];
ybs[4434]=['',3.0252742,-0.5127771,5.64];
ybs[4435]=['',3.0257977,-0.4688919,6.16];
ybs[4436]=['',3.0276484,-0.1386906,5.95];
ybs[4437]=['',3.0275184,-0.7078231,5.64];
ybs[4438]=['',3.0250919,-1.1707857,5.9];
ybs[4439]=['',3.028016,-0.5446495,5.04];
ybs[4440]=['ξ Hya',3.0284484,-0.5580983,3.54];
ybs[4441]=['',3.0295829,-0.2862247,6.05];
ybs[4442]=['',3.0328585,0.6404765,6.4];
ybs[4443]=['',3.0311017,-0.7104516,5.39];
ybs[4444]=['',3.0337261,0.1903219,6.55];
ybs[4445]=['89 Leo',3.0345661,0.0513305,5.77];
ybs[4446]=['90 Leo',3.0361127,0.2910851,5.95];
ybs[4447]=['',3.0379869,0.9541062,5.63];
ybs[4448]=['',3.0349528,-0.5750925,5.98];
ybs[4449]=['',3.0376764,0.3546922,6.45];
ybs[4450]=['',3.0359644,-0.9491654,4.62];
ybs[4451]=['2 Dra',3.0424455,1.2078324,5.2];
ybs[4452]=['',3.0368264,-0.8596737,5.5];
ybs[4453]=['',3.0380373,-0.8288835,5.71];
ybs[4454]=['',3.0405152,0.1883569,6.56];
ybs[4455]=['',3.0430933,0.4827934,5.8];
ybs[4456]=['',3.0411255,-0.833582,5.25];
ybs[4457]=['λ Cen',3.0402961,-1.1019795,3.13];
ybs[4458]=['θ Crt',3.0446193,-0.1731599,4.7];
ybs[4459]=['',3.0440878,-0.5879857,5.74];
ybs[4460]=['',3.0444901,-0.6520006,6.31];
ybs[4461]=['υ Leo',3.0458157,-0.0164586,4.3];
ybs[4462]=['',3.0429214,-1.0676408,5.83];
ybs[4463]=['',3.0459987,-0.5778293,6.29];
ybs[4464]=['',3.0501369,0.8813767,6.14];
ybs[4465]=['',3.0457136,-1.071675,5.15];
ybs[4466]=['',3.0482896,-0.8354258,5.44];
ybs[4467]=['59 UMa',3.0520882,0.7593293,5.59];
ybs[4468]=['',3.0511504,0.1529779,6.17];
ybs[4469]=['π Cha',3.0463681,-1.3267259,5.65];
ybs[4470]=['60 UMa',3.0530459,0.81533,6.1];
ybs[4471]=['',3.0543747,1.1209853,6.46];
ybs[4472]=['',3.0528755,0.5847963,6.27];
ybs[4473]=['ω Vir',3.0524419,0.1398877,5.36];
ybs[4474]=['',3.0521518,-0.0445985,6.22];
ybs[4475]=['',3.049083,-1.1822762,5.96];
ybs[4476]=['',3.0538614,0.7852132,6.44];
ybs[4477]=['',3.0505744,-1.0811541,5.15];
ybs[4478]=['ι Crt',3.0532792,-0.2324979,5.48];
ybs[4479]=['',3.0547141,-0.4335456,6.42];
ybs[4480]=['',3.0584408,-0.2546064,6.21];
ybs[4481]=['',3.0583827,-0.29216,6.19];
ybs[4482]=['',3.0565101,-1.1434876,5.17];
ybs[4483]=['',3.0614022,1.0096952,6.37];
ybs[4484]=['ο Hya',3.0599433,-0.6084915,4.7];
ybs[4485]=['92 Leo',3.0626257,0.3705941,5.26];
ybs[4486]=['61 UMa',3.0638293,0.5948493,5.33];
ybs[4487]=['',3.0619969,-0.9440121,5.96];
ybs[4488]=['',3.064015,-0.5116555,6.44];
ybs[4489]=['',3.0627118,-1.0857571,4.94];
ybs[4490]=['',3.0668994,0.9608589,6.27];
ybs[4491]=['62 UMa',3.0660924,0.5519914,5.73];
ybs[4492]=['',3.0647826,-0.7542468,5.55];
ybs[4493]=['',3.0665943,-0.56931,5.22];
ybs[4494]=['3 Dra',3.070284,1.1628366,5.3];
ybs[4495]=['',3.0683034,0.3855691,6.59];
ybs[4496]=['',3.0680567,-0.3562783,6.22];
ybs[4497]=['',3.0621737,-1.4524508,6.33];
ybs[4498]=['',3.0740914,-0.6511768,5.98];
ybs[4499]=['',3.0710836,-1.3862412,6.39];
ybs[4500]=['',3.0762145,-0.1186238,6.07];
ybs[4501]=['',3.0742133,-1.0927306,5.03];
ybs[4502]=['',3.0776126,0.4380584,6.02];
ybs[4503]=['',3.0757932,-1.0995182,6.1];
ybs[4504]=['ζ Crt',3.0798737,-0.3223673,4.73];
ybs[4505]=['ξ Vir',3.0822131,0.14205,4.85];
ybs[4506]=['',3.0817171,-0.8585133,6.26];
ybs[4507]=['ν Vir',3.0847174,0.1118749,4.03];
ybs[4508]=['χ UMa',3.0856665,0.8318231,3.71];
ybs[4509]=['',3.084015,-0.7995262,5.29];
ybs[4510]=['λ Mus',3.0833032,-1.1667193,3.64];
ybs[4511]=['',3.0895238,0.9688116,5.27];
ybs[4512]=['',3.087339,-1.0698491,4.11];
ybs[4513]=['',3.087475,-0.7089538,4.91];
ybs[4514]=['',3.090105,-0.6287805,6.17];
ybs[4515]=['',3.0907542,-0.530693,6.48];
ybs[4516]=['',3.090891,-1.0090782,5.41];
ybs[4517]=['93 Leo',3.0940137,0.3507997,4.53];
ybs[4518]=['4 Vir',3.0936851,0.1418304,5.32];
ybs[4519]=['',3.0957313,-0.1820883,6.26];
ybs[4520]=['μ Mus',3.094841,-1.1682235,4.72];
ybs[4521]=['',3.0968801,0.2472189,5.88];
ybs[4522]=['',3.0972713,-0.4689575,5.11];
ybs[4523]=['',3.0984909,-0.0076477,6.15];
ybs[4524]=['β Leo',3.0986905,0.2522415,2.14];
ybs[4525]=['',3.0995148,0.281403,6.04];
ybs[4526]=['',3.1014986,0.6075854,5.7];
ybs[4527]=['',3.1012115,-1.1154036,4.32];
ybs[4528]=['',3.1022778,-1.2277593,4.97];
ybs[4529]=['',3.1041619,-0.2789645,6.13];
ybs[4530]=['β Vir',3.1058024,0.0287126,3.61];
ybs[4531]=['',3.1045864,-1.0955266,5.7];
ybs[4532]=['',3.1054239,-0.4781745,6.48];
ybs[4533]=['',3.1068059,0.2122194,6.35];
ybs[4534]=['',3.1072835,-0.0951719,5.64];
ybs[4535]=['',3.107864,0.5804159,6.27];
ybs[4536]=['',3.1076871,-0.7905159,4.46];
ybs[4537]=['',3.1087075,-0.2148047,6.35];
ybs[4538]=['',3.1101138,-0.5402602,5.85];
ybs[4539]=['',3.1107036,-1.1401492,4.9];
ybs[4540]=['',3.1158157,0.6562257,6.45];
ybs[4541]=['',3.1121262,-0.9967124,5.57];
ybs[4542]=['β Hya',3.1154246,-0.5938954,4.28];
ybs[4543]=['',3.1177679,-0.6141172,6.17];
ybs[4544]=['γ UMa',3.1195484,0.9350612,2.44];
ybs[4545]=['',3.1195158,0.0075447,6.3];
ybs[4546]=['',3.1209818,-1.004082,6.06];
ybs[4547]=['',3.1220615,-0.660931,6.46];
ybs[4548]=['',3.1232906,-0.4508806,5.3];
ybs[4549]=['6 Vir',3.1248155,0.145285,5.58];
ybs[4550]=['65 UMa',3.1250403,0.809087,6.54];
ybs[4551]=['65 UMa',3.1254394,0.808961,7.03];
ybs[4552]=['',3.1256381,0.6394313,6.49];
ybs[4553]=['',3.1244886,-1.1065137,5.91];
ybs[4554]=['95 Leo',3.127539,0.2709971,5.53];
ybs[4555]=['',3.1274815,-0.4991052,5.93];
ybs[4556]=['66 UMa',3.1288789,0.9857432,5.84];
ybs[4557]=['η Crt',3.1290032,-0.3014274,5.18];
ybs[4558]=['',3.1285351,-0.6947954,6.13];
ybs[4559]=['',3.1328592,1.0721466,6.22];
ybs[4560]=['',3.1321133,-0.823659,6.26];
ybs[4561]=['',3.1335655,-0.5835503,6.21];
ybs[4562]=['',3.1343895,0.7020398,6.62];
ybs[4563]=['',3.1361956,-1.0920277,5.57];
ybs[4564]=['',3.1382053,0.5611965,6.42];
ybs[4565]=['',3.1391906,1.0706727,6.76];
ybs[4566]=['',3.1387631,-0.9850101,5.44];
ybs[4567]=['',3.1391416,-0.7167529,6.79];
ybs[4568]=['',3.1411304,-1.1250194,5.61];
ybs[4569]=['',3.1416278,-0.4542845,6.43];
ybs[4570]=['',3.1422847,0.0071708,6.17];
ybs[4571]=['',3.1433182,0.576793,5.96];
ybs[4572]=['',3.142827,-0.9043661,6.05];
ybs[4573]=['ε Cha',3.1447576,-1.3673196,4.91];
ybs[4574]=['',3.146194,0.5919337,6.5];
ybs[4575]=['7 Vir',3.146175,0.0617075,5.37];
ybs[4576]=['',3.1477053,1.4090629,6.17];
ybs[4577]=['',3.1496391,-0.1844081,5.55];
ybs[4578]=['',3.149496,-0.3832205,6.28];
ybs[4579]=['π Vir',3.1502099,0.1133499,4.66];
ybs[4580]=['',3.1501287,-0.3452014,5.26];
ybs[4581]=['',3.1508954,-0.0329474,6.31];
ybs[4582]=['',3.1529022,-1.0057163,6.16];
ybs[4583]=['',3.153622,0.6269616,5.59];
ybs[4584]=['67 UMa',3.1555993,0.7491977,5.21];
ybs[4585]=['',3.1569432,-1.4922957,6.05];
ybs[4586]=['',3.1572903,-1.2498054,6.42];
ybs[4587]=['',3.1579457,-1.209721,5.89];
ybs[4588]=['',3.1588841,-0.1361932,6.22];
ybs[4589]=['θ1 Cru',3.1596649,-1.1071052,4.33];
ybs[4590]=['',3.1624057,-0.7427046,5.15];
ybs[4591]=['',3.1628516,-1.2973654,6.44];
ybs[4592]=['2 Com',3.1650439,0.3724445,5.87];
ybs[4593]=['θ2 Cru',3.1653385,-1.1045355,4.72];
ybs[4594]=['',3.1667882,-1.1946575,5.35];
ybs[4595]=['κ Cha',3.1674399,-1.3375999,5.04];
ybs[4596]=['',3.1653203,1.4872683,6.27];
ybs[4597]=['',3.1680972,-1.0662012,5.96];
ybs[4598]=['ο Vir',3.1691184,0.1503322,4.12];
ybs[4599]=['',3.1690902,1.3401716,5.8];
ybs[4600]=['',3.1709941,1.0963007,6.13];
ybs[4601]=['',3.1722159,-1.1461031,6.33];
ybs[4602]=['',3.1723833,-0.6250642,6.23];
ybs[4603]=['',3.1725694,-0.0567461,6.37];
ybs[4604]=['',3.174182,-1.2002712,6.23];
ybs[4605]=['',3.1744014,-1.1489294,6.06];
ybs[4606]=['η Cru',3.1765711,-1.1298083,4.15];
ybs[4607]=['',3.1808606,-1.317489,5.18];
ybs[4608]=['',3.1817819,-0.8862956,4.47];
ybs[4609]=['',3.181753,-0.8880749,6.37];
ybs[4610]=['',3.1824679,-0.8519369,5.34];
ybs[4611]=['δ Cen',3.1829705,-0.8873621,2.6];
ybs[4612]=['',3.1832433,-1.0640719,6.22];
ybs[4613]=['α Crv',3.1831485,-0.433688,4.02];
ybs[4614]=['',3.1853041,-0.7757239,5.75];
ybs[4615]=['',3.1853465,-0.7217108,5.48];
ybs[4616]=['10 Vir',3.1886711,0.0310355,5.95];
ybs[4617]=['',3.1887696,1.3010001,6.35];
ybs[4618]=['',3.1902828,-0.6078034,6.17];
ybs[4619]=['11 Vir',3.1902715,0.0992634,5.72];
ybs[4620]=['ε Crv',3.1906201,-0.3968754,3];
ybs[4621]=['',3.192571,-0.6630476,6.06];
ybs[4622]=['3 Com',3.1923016,0.2912886,6.39];
ybs[4623]=['',3.1933335,0.4740636,6.01];
ybs[4624]=['',3.1949479,-1.0715805,6.08];
ybs[4625]=['3 Crv',3.1947273,-0.4140277,5.46];
ybs[4626]=['',3.1947158,-0.7948634,6.61];
ybs[4627]=['',3.1968189,-0.8984776,6.23];
ybs[4628]=['ρ Cen',3.1973854,-0.9160908,3.96];
ybs[4629]=['',3.1936853,1.424022,6];
ybs[4630]=['4 Com',3.1980661,0.4494355,5.66];
ybs[4631]=['68 UMa',3.1974902,0.9937018,6.43];
ybs[4632]=['',3.1987858,0.4959631,6.49];
ybs[4633]=['5 Com',3.1993937,0.3564387,5.57];
ybs[4634]=['',3.200592,-1.1007851,5.92];
ybs[4635]=['',3.2025025,-1.2264679,6.17];
ybs[4636]=['',3.1991021,1.3525756,5.14];
ybs[4637]=['',3.2041545,-0.5976886,6.5];
ybs[4638]=['',3.2050669,-0.6815273,5.76];
ybs[4639]=['',3.2078341,-1.3734531,6.35];
ybs[4640]=['12 Vir',3.2049846,0.1770244,5.85];
ybs[4641]=['',3.2058795,-0.5918803,6.33];
ybs[4642]=['',3.2078132,-0.8001173,5.31];
ybs[4643]=['',3.2089924,-1.126227,6.22];
ybs[4644]=['',3.2104565,0.9305273,6.16];
ybs[4645]=['',3.2118756,-0.3658836,5.83];
ybs[4646]=['δ Cru',3.2127198,-1.0274457,2.8];
ybs[4647]=['',3.2126477,-0.1820713,6.11];
ybs[4648]=['',3.2142003,-0.7336048,6.26];
ybs[4649]=['',3.2120741,1.2231369,5.71];
ybs[4650]=['δ UMa',3.2134888,0.9933209,3.31];
ybs[4651]=['',3.2153329,-0.4096812,6.54];
ybs[4652]=['γ Crv',3.2154176,-0.3082485,2.59];
ybs[4653]=['6 Com',3.2161886,0.257951,5.1];
ybs[4654]=['',3.2184192,-1.2694494,6.22];
ybs[4655]=['',3.2143939,1.264167,6.29];
ybs[4656]=['2 CVn',3.2166365,0.7075721,5.66];
ybs[4657]=['7 Com',3.2176366,0.4158405,4.95];
ybs[4658]=['',3.2183021,0.5749468,5];
ybs[4659]=['',3.2213713,-1.1486382,6.06];
ybs[4660]=['',3.2208655,-0.2934414,6.05];
ybs[4661]=['ε Mus',3.2234669,-1.1882228,4.11];
ybs[4662]=['',3.2225037,0.9262774,5.81];
ybs[4663]=['',3.2227046,0.5029672,5.7];
ybs[4664]=['β Cha',3.2273678,-1.3863413,4.26];
ybs[4665]=['',3.2241391,-0.6320396,6.15];
ybs[4666]=['',3.2237539,0.2622332,6.34];
ybs[4667]=['',3.2256142,-0.0711002,6.99];
ybs[4668]=['',3.2256506,-0.0709984,6.54];
ybs[4669]=['ζ Cru',3.2271948,-1.119146,4.04];
ybs[4670]=['',3.2271376,0.5258657,6.23];
ybs[4671]=['13 Vir',3.2278739,-0.0158214,5.9];
ybs[4672]=['',3.2295385,-0.9645094,5];
ybs[4673]=['',3.217489,1.503598,6.33];
ybs[4674]=['',3.2293699,0.4518399,6.48];
ybs[4675]=['8 Com',3.2306162,0.3999504,6.27];
ybs[4676]=['',3.2099015,1.5273381,6.28];
ybs[4677]=['',3.2279087,1.3097175,5.38];
ybs[4678]=['9 Com',3.2313591,0.4893503,6.33];
ybs[4679]=['η Vir',3.2332623,-0.0137212,3.89];
ybs[4680]=['3 CVn',3.2326331,0.8528541,5.29];
ybs[4681]=['',3.2345225,-0.389117,5.97];
ybs[4682]=['',3.2361321,-1.1512536,6.21];
ybs[4683]=['',3.2350048,0.4625166,5.54];
ybs[4684]=['',3.2348621,0.451739,6.15];
ybs[4685]=['16 Vir',3.2351836,0.0557336,4.96];
ybs[4686]=['ζ Crv',3.2361966,-0.3898196,5.21];
ybs[4687]=['11 Com',3.2367324,0.3084625,4.74];
ybs[4688]=['',3.2365722,0.4701139,7.13];
ybs[4689]=['',3.2377652,-0.2388435,5.14];
ybs[4690]=['ε Cru',3.2399538,-1.0562778,3.59];
ybs[4691]=['70 UMa',3.2370527,1.0078354,5.55];
ybs[4692]=['',3.2425213,-0.9860035,5.92];
ybs[4693]=['ζ2 Mus',3.2434261,-1.180559,5.15];
ybs[4694]=['ζ1 Mus',3.2437827,-1.1942696,5.74];
ybs[4695]=['',3.2430841,0.4303071,6.19];
ybs[4696]=['',3.2463165,-1.0087162,5.39];
ybs[4697]=['12 Com',3.244496,0.4490212,4.81];
ybs[4698]=['17 Vir',3.2447017,0.0905209,6.4];
ybs[4699]=['',3.2618621,-1.5022889,6.33];
ybs[4700]=['',3.2482894,-1.1824731,6.36];
ybs[4701]=['6 Crv',3.248429,-0.4356272,5.68];
ybs[4702]=['',3.2494881,-0.6201471,5.32];
ybs[4703]=['',3.2496211,-0.6880452,6.4];
ybs[4704]=['',3.2502014,-0.6812043,5.79];
ybs[4705]=['4 CVn',3.2499725,0.7404341,6.06];
ybs[4706]=['5 CVn',3.2509488,0.8978534,4.8];
ybs[4707]=['13 Com',3.2523554,0.4534299,5.18];
ybs[4708]=['',3.2545746,-0.7243664,6.25];
ybs[4709]=['',3.2529537,0.444427,6.42];
ybs[4710]=['',3.2572668,-1.1499885,6.3];
ybs[4711]=['',3.2563169,-0.744093,6.11];
ybs[4712]=['',3.2563831,-0.2047135,5.95];
ybs[4713]=['',3.256946,-0.4863901,6.09];
ybs[4714]=['',3.2572359,-0.616194,5.73];
ybs[4715]=['',3.2564781,0.4155135,6.03];
ybs[4716]=['71 UMa',3.2553622,0.9888782,5.81];
ybs[4717]=['',3.2554771,1.1114924,6.32];
ybs[4718]=['6 CVn',3.2589914,0.6789279,5.02];
ybs[4719]=['',3.2625822,-1.10377,4.86];
ybs[4720]=['α1 Cru',3.2629468,-1.1033626,1.33];
ybs[4721]=['α2 Cru',3.2629906,-1.1033675,1.73];
ybs[4722]=['',3.2624553,-0.9000609,4.82];
ybs[4723]=['14 Com',3.2614738,0.4738475,4.95];
ybs[4724]=['',3.2636386,-0.8557729,6.26];
ybs[4725]=['',3.2637682,-0.5750657,5.55];
ybs[4726]=['',3.2665183,-1.1154045,6];
ybs[4727]=['γ Com',3.2638074,0.4913014,4.36];
ybs[4728]=['16 Com',3.2640333,0.4661203,5];
ybs[4729]=['',3.2667233,-1.0316771,5.5];
ybs[4730]=['',3.2608607,1.2533357,6.24];
ybs[4731]=['',3.2672335,0.1482045,6.37];
ybs[4732]=['',3.2678803,-0.2923553,6.35];
ybs[4733]=['σ Cen',3.2690617,-0.8787614,3.91];
ybs[4734]=['',3.2704971,-1.1250416,6.04];
ybs[4735]=['73 UMa',3.2663892,0.970298,5.7];
ybs[4736]=['',3.267984,-0.0826249,6.22];
ybs[4737]=['',3.2709216,-1.0806035,6.22];
ybs[4738]=['',3.2704239,-0.6834732,5.44];
ybs[4739]=['',3.2714026,-0.9865736,6.15];
ybs[4740]=['',3.2712142,0.4556695,6.54];
ybs[4741]=['',3.2716883,0.4499536,6.65];
ybs[4742]=['17 Com',3.272422,0.4501914,5.29];
ybs[4743]=['18 Com',3.2747717,0.4187083,5.48];
ybs[4744]=['',3.2772769,-0.9886131,5.8];
ybs[4745]=['',3.2773958,-0.7305032,6.02];
ybs[4746]=['20 Com',3.275974,0.362635,5.69];
ybs[4747]=['δ Crv',3.2767937,-0.2903215,2.95];
ybs[4748]=['',3.2777161,-0.2358234,6.35];
ybs[4749]=['',3.2786965,-0.4156551,5.63];
ybs[4750]=['74 UMa',3.2766551,1.0173034,5.35];
ybs[4751]=['7 CVn',3.2771624,0.8973945,6.21];
ybs[4752]=['75 UMa',3.2771561,1.0236159,6.08];
ybs[4753]=['γ Cru',3.2828311,-0.9988848,1.63];
ybs[4754]=['γ Cru',3.2833266,-0.9983222,6.42];
ybs[4755]=['4 Dra',3.2770583,1.2057167,4.95];
ybs[4756]=['21 Com',3.2815767,0.4267096,5.46];
ybs[4757]=['',3.2805656,0.924293,6.21];
ybs[4758]=['',3.285092,-1.0392109,5.48];
ybs[4759]=['',3.2877335,-1.2761871,5.88];
ybs[4760]=['',3.2831816,0.1306488,6.05];
ybs[4761]=['',3.2863165,-1.1104588,5.95];
ybs[4762]=['',3.2845043,-0.0902513,6.19];
ybs[4763]=['γ Mus',3.289,-1.2610266,3.87];
ybs[4764]=['',3.2865439,-0.5698865,6.46];
ybs[4765]=['η Crv',3.2864189,-0.2847435,4.31];
ybs[4766]=['',3.2887197,-0.2439554,5.74];
ybs[4767]=['20 Vir',3.2905462,0.1776246,6.26];
ybs[4768]=['',3.2921307,-0.3475008,6.26];
ybs[4769]=['',3.2929557,-0.2259966,5.58];
ybs[4770]=['22 Com',3.2927381,0.4217533,6.29];
ybs[4771]=['21 Vir',3.293846,-0.1670333,5.48];
ybs[4772]=['',3.2950694,-0.8731495,6.38];
ybs[4773]=['',3.2930163,0.5782124,5.42];
ybs[4774]=['',3.2936325,0.5806076,6.24];
ybs[4775]=['β CVn',3.2933532,0.7197587,4.26];
ybs[4776]=['β Crv',3.2965786,-0.4104137,2.65];
ybs[4777]=['κ Dra',3.2916648,1.21597,3.87];
ybs[4778]=['',3.2981473,-0.7817612,5.77];
ybs[4779]=['23 Com',3.2983362,0.3928893,4.81];
ybs[4780]=['',3.3018398,-1.0814087,6.22];
ybs[4781]=['24 Com',3.2994681,0.3186791,6.56];
ybs[4782]=['24 Com',3.2995771,0.3186744,5.02];
ybs[4783]=['',3.2995764,0.3798385,5.85];
ybs[4784]=['',3.3027135,-0.7180309,5.13];
ybs[4785]=['6 Dra',3.2970776,1.220049,4.94];
ybs[4786]=['',3.3038458,-0.6979253,5.8];
ybs[4787]=['',3.3035032,-0.3603302,6.2];
ybs[4788]=['α Mus',3.3095473,-1.2087038,2.69];
ybs[4789]=['25 Vir',3.3069624,-0.103848,5.87];
ybs[4790]=['',3.3046082,1.0361809,5.5];
ybs[4791]=['25 Com',3.3076158,0.2962059,5.68];
ybs[4792]=['τ Cen',3.3113067,-0.8492622,3.86];
ybs[4793]=['',3.311091,-0.475723,5.45];
ybs[4794]=['',3.3190187,-1.3175024,6.49];
ybs[4795]=['',3.3125073,0.0552309,6.33];
ybs[4796]=['',3.3168799,-1.1747982,6.25];
ybs[4797]=['',3.3138249,0.030312,5.71];
ybs[4798]=['',3.3143456,0.1199106,7.08];
ybs[4799]=['',3.3155684,-0.3205858,6];
ybs[4800]=['',3.3170365,-0.5330258,5.89];
ybs[4801]=['9 CVn',3.3152663,0.7113351,6.37];
ybs[4802]=['',3.3165737,0.3934239,6.38];
ybs[4803]=['χ Vir',3.3176977,-0.1416065,4.66];
ybs[4804]=['',3.3214627,-1.162904,6.26];
ybs[4805]=['26 Com',3.3169631,0.3655521,5.46];
ybs[4806]=['',3.3175381,0.6254221,6.45];
ybs[4807]=['',3.3206967,-0.6999701,4.64];
ybs[4808]=['',3.3273665,-0.807446,5.84];
ybs[4809]=['γ Cen',3.3279908,-0.8565623,2.17];
ybs[4810]=['',3.3310657,-1.2134422,6.33];
ybs[4811]=['',3.3265486,-0.2291848,6.08];
ybs[4812]=['',3.3265631,-0.229209,5.98];
ybs[4813]=['',3.3300745,-1.0437675,4.93];
ybs[4814]=['27 Vir',3.3277287,0.1799209,6.19];
ybs[4815]=['γ Vir',3.3281869,-0.0273514,3.65];
ybs[4816]=['γ Vir',3.3281869,-0.0273514,3.68];
ybs[4817]=['',3.3290133,-0.3469063,6.03];
ybs[4818]=['ρ Vir',3.3290895,0.1765907,4.88];
ybs[4819]=['31 Vir',3.3294038,0.1167454,5.59];
ybs[4820]=['',3.3340953,-1.102632,5.31];
ybs[4821]=['',3.3326829,-0.8540006,4.66];
ybs[4822]=['',3.333851,-0.9785149,6.08];
ybs[4823]=['76 UMa',3.3270256,1.0924952,6.07];
ybs[4824]=['',3.3352785,-0.9825092,6];
ybs[4825]=['',3.3367382,-1.0301028,6.4];
ybs[4826]=['',3.3362737,-0.7032852,6.44];
ybs[4827]=['',3.3368058,-0.0295733,5.93];
ybs[4828]=['',3.3385901,-0.6364623,6.39];
ybs[4829]=['',3.3386415,-0.4963948,5.48];
ybs[4830]=['',3.3336146,1.0653143,6.38];
ybs[4831]=['',3.3439458,-1.2033772,6.16];
ybs[4832]=['ι Cru',3.3462634,-1.0663679,4.69];
ybs[4833]=['',3.3399639,0.7676945,6.33];
ybs[4834]=['β Mus',3.3494034,-1.1907552,3.05];
ybs[4835]=['10 CVn',3.3423776,0.6834978,5.95];
ybs[4836]=['',3.3428987,0.7910347,4.99];
ybs[4837]=['32 Vir',3.3453924,0.131878,5.22];
ybs[4838]=['',3.3494083,-0.9879624,4.65];
ybs[4839]=['33 Vir',3.348679,0.1644589,5.67];
ybs[4840]=['',3.3507411,-0.5835108,5.86];
ybs[4841]=['27 Com',3.3498013,0.287287,5.12];
ybs[4842]=['',3.3378417,1.4050543,6.4];
ybs[4843]=['β Cru',3.3553686,-1.0438056,1.25];
ybs[4844]=['',3.3515997,0.1018174,6.34];
ybs[4845]=['34 Vir',3.3523749,0.2066635,6.07];
ybs[4846]=['',3.3539544,-0.1120329,6.26];
ybs[4847]=['',3.3555857,-0.4357859,6.44];
ybs[4848]=['35 Vir',3.355188,0.060314,6.41];
ybs[4849]=['',3.352013,1.0936883,5.89];
ybs[4850]=['',3.3579924,-0.4837087,5.66];
ybs[4851]=['28 Com',3.3567743,0.2345034,6.56];
ybs[4852]=['',3.3648768,-1.2584382,5.55];
ybs[4853]=['7 Dra',3.3529683,1.1636669,5.43];
ybs[4854]=['',3.3590523,0.4315037,6.31];
ybs[4855]=['29 Com',3.3596698,0.2444434,5.7];
ybs[4856]=['11 CVn',3.3583787,0.8438666,6.27];
ybs[4857]=['',3.3579277,1.0507413,5.85];
ybs[4858]=['',3.3662507,-1.0562313,6.75];
ybs[4859]=['30 Com',3.3612345,0.4788371,5.78];
ybs[4860]=['ι Oct',3.3922997,-1.4823235,5.46];
ybs[4861]=['',3.3665131,-0.8478194,6.24];
ybs[4862]=['',3.3693928,-0.9233521,5.73];
ybs[4863]=['',3.3656391,0.3970025,6.43];
ybs[4864]=['',3.3678764,-0.5954392,4.91];
ybs[4865]=['',3.3649939,0.652756,5.89];
ybs[4866]=['',3.3710506,-1.054988,5.72];
ybs[4867]=['',3.3706817,-0.1824737,6.41];
ybs[4868]=['37 Vir',3.371589,0.0513137,6.02];
ybs[4869]=['',3.3734555,-0.6945956,5.98];
ybs[4870]=['',3.3742107,-0.8414357,6.33];
ybs[4871]=['',3.3733806,-0.4687015,6.15];
ybs[4872]=['',3.3757393,-0.9415344,6.24];
ybs[4873]=['31 Com',3.3717318,0.4786383,4.94];
ybs[4874]=['32 Com',3.3740422,0.2959616,6.32];
ybs[4875]=['',3.3786334,-0.9611341,5.93];
ybs[4876]=['',3.3751628,0.2793572,6.3];
ybs[4877]=['',3.3801004,-1.0549642,5.76];
ybs[4878]=['',3.378717,-0.8562542,4.33];
ybs[4879]=['',3.3799751,-0.7032852,4.27];
ybs[4880]=['κ Cru',3.3820945,-1.0558069,5.9];
ybs[4881]=['38 Vir',3.378502,-0.0640444,6.11];
ybs[4882]=['',3.3568591,1.4538786,5.85];
ybs[4883]=['',3.3573623,1.4537867,5.28];
ybs[4884]=['35 Com',3.3787628,0.3687635,4.9];
ybs[4885]=['',3.3844302,-1.0218348,6.58];
ybs[4886]=['',3.3804642,-0.0757518,6.44];
ybs[4887]=['λ Cru',3.3857085,-1.0343327,4.62];
ybs[4888]=['μ1 Cru',3.3853838,-0.9999692,4.03];
ybs[4889]=['μ2 Cru',3.3854711,-0.9998043,5.17];
ybs[4890]=['41 Vir',3.3811636,0.2147151,6.25];
ybs[4891]=['',3.3834827,-0.2053361,6];
ybs[4892]=['ψ Vir',3.3836457,-0.1685144,4.79];
ybs[4893]=['',3.3867624,-0.7726248,5.89];
ybs[4894]=['',3.3826461,0.5832567,6.26];
ybs[4895]=['ε UMa',3.3814169,0.9746511,1.77];
ybs[4896]=['',3.3882646,-0.7510498,5.47];
ybs[4897]=['',3.3946506,-1.261895,5.93];
ybs[4898]=['',3.3913118,-0.994003,5.32];
ybs[4899]=['',3.3855898,0.8217089,5.84];
ybs[4900]=['δ Vir',3.3889865,0.0572709,3.38];
ybs[4901]=['',3.3903998,-0.2695316,6.17];
ybs[4902]=['',3.3931927,-0.4638436,6.62];
ybs[4903]=['',3.3960787,-0.8956075,5.16];
ybs[4904]=['α1 CVn',3.3903738,0.6666922,5.6];
ybs[4905]=['α2 CVn',3.3904682,0.6667552,2.9];
ybs[4906]=['8 Dra',3.3873551,1.1400919,5.24];
ybs[4907]=['',3.3913208,0.9421881,5.82];
ybs[4908]=['',3.397745,-0.3991526,6.31];
ybs[4909]=['',3.3951297,0.8039163,6.12];
ybs[4910]=['36 Com',3.4033348,0.3018329,4.78];
ybs[4911]=['44 Vir',3.406746,-0.0685484,5.79];
ybs[4912]=['',3.4109363,-0.5867926,6.02];
ybs[4913]=['δ Mus',3.4198106,-1.2507743,3.62];
ybs[4914]=['37 Com',3.4090746,0.5352836,4.9];
ybs[4915]=['46 Vir',3.4108439,-0.0608086,5.99];
ybs[4916]=['',3.4108394,0.3186552,6.2];
ybs[4917]=['',3.4008656,1.3152236,6.01];
ybs[4918]=['9 Dra',3.4065913,1.1603237,5.32];
ybs[4919]=['38 Com',3.413091,0.2968398,5.96];
ybs[4920]=['',3.4233883,-1.2495021,6.03];
ybs[4921]=['78 UMa',3.4105587,0.981764,4.93];
ybs[4922]=['ε Vir',3.4175889,0.1892621,2.83];
ybs[4923]=['ξ1 Cen',3.4243826,-0.8664208,4.85];
ybs[4924]=['',3.414841,1.1081962,6];
ybs[4925]=['',3.4248606,-0.3612494,5.58];
ybs[4926]=['',3.4188943,1.0402324,6.53];
ybs[4927]=['48 Vir',3.4252839,-0.0659442,6.59];
ybs[4928]=['',3.4296673,-0.7210221,6.26];
ybs[4929]=['',3.4330225,-0.911581,6.43];
ybs[4930]=['',3.4362746,-0.8478503,4.71];
ybs[4931]=['',3.4374695,-0.7278581,5.59];
ybs[4932]=['ξ2 Cen',3.4390729,-0.873025,4.27];
ybs[4933]=['14 CVn',3.4328243,0.6228061,5.25];
ybs[4934]=['',3.4415588,-1.0467614,5.99];
ybs[4935]=['',3.4332192,0.7880842,5.63];
ybs[4936]=['39 Com',3.4356899,0.3671945,5.99];
ybs[4937]=['',3.4387682,-0.6279081,6.54];
ybs[4938]=['',3.43479,0.5046581,6.54];
ybs[4939]=['40 Com',3.4357738,0.3927249,5.6];
ybs[4940]=['',3.4273548,1.2725263,6.31];
ybs[4941]=['',3.442358,-0.9350453,5.71];
ybs[4942]=['θ Mus',3.4449632,-1.1418071,5.51];
ybs[4943]=['',3.4349026,1.0808352,6.14];
ybs[4944]=['41 Com',3.439204,0.4801438,4.8];
ybs[4945]=['49 Vir',3.442773,-0.1894498,5.19];
ybs[4946]=['',3.4423215,0.4789433,6.19];
ybs[4947]=['',3.4455689,-0.1588029,5.55];
ybs[4948]=['ψ Hya',3.4479776,-0.4054796,4.95];
ybs[4949]=['',3.4486231,-0.1684682,6.32];
ybs[4950]=['',3.4482618,0.1729277,5.78];
ybs[4951]=['50 Vir',3.45088,-0.1822743,5.94];
ybs[4952]=['',3.4507558,0.2920722,5.91];
ybs[4953]=['θ Vir',3.4516772,-0.0986628,4.38];
ybs[4954]=['',3.4498037,0.6511635,6.02];
ybs[4955]=['',3.4569211,-0.9194542,6.06];
ybs[4956]=['',3.4617374,-1.2227023,5.91];
ybs[4957]=['15 CVn',3.4500238,0.6705513,6.28];
ybs[4958]=['α Com',3.4515833,0.3039555,5.22];
ybs[4959]=['α Com',3.4515833,0.3039555,5.22];
ybs[4960]=['',3.4574219,-0.7390934,5.79];
ybs[4961]=['17 CVn',3.4515637,0.6699415,5.91];
ybs[4962]=['',3.4613552,-1.1068271,6.33];
ybs[4963]=['',3.4584949,-0.7589167,5.25];
ybs[4964]=['',3.4499265,1.0841121,6.54];
ybs[4965]=['',3.4629571,-1.0477999,4.6];
ybs[4966]=['',3.4738485,-1.3711399,5.85];
ybs[4967]=['',3.4656159,-1.1578607,5.9];
ybs[4968]=['',3.4593718,-0.4654002,6.5];
ybs[4969]=['',3.461299,-0.6617728,4.85];
ybs[4970]=['',3.465769,-1.0459801,6.16];
ybs[4971]=['53 Vir',3.4610043,-0.2847042,5.04];
ybs[4972]=['',3.4648702,-0.7472334,6.22];
ybs[4973]=['β Com',3.4596649,0.4845781,4.26];
ybs[4974]=['',3.4608782,0.4213979,6.33];
ybs[4975]=['',3.4674442,-0.8868629,5.89];
ybs[4976]=['',3.4628193,0.1997084,5.77];
ybs[4977]=['',3.4629489,0.3252946,6.53];
ybs[4978]=['',3.4712711,-1.0262057,5.89];
ybs[4979]=['',3.4714866,-1.0335262,4.92];
ybs[4980]=['54 Vir',3.4670991,-0.3305684,6.28];
ybs[4981]=['',3.469732,-0.7548951,6.16];
ybs[4982]=['',3.4656019,0.3248649,6.11];
ybs[4983]=['η Mus',3.4764126,-1.1869569,4.8];
ybs[4984]=['',3.4773656,-1.2181153,6.37];
ybs[4985]=['55 Vir',3.4703228,-0.3498376,5.33];
ybs[4986]=['',3.4731974,-0.8564322,5.89];
ybs[4987]=['',3.4674908,0.6988177,4.92];
ybs[4988]=['',3.4714289,0.1957968,5.67];
ybs[4989]=['',3.4748559,-0.6367716,6.19];
ybs[4990]=['',3.4827703,-1.1388493,6.07];
ybs[4991]=['57 Vir',3.4781753,-0.3500456,5.22];
ybs[4992]=['',3.4849513,-1.1675633,4.87];
ybs[4993]=['',3.465121,1.2685987,6.59];
ybs[4994]=['19 CVn',3.4753899,0.7110839,5.79];
ybs[4995]=['',3.4798851,-0.0262422,6.68];
ybs[4996]=['',3.4823006,-0.5518562,5.1];
ybs[4997]=['',3.4788168,0.3305413,6.45];
ybs[4998]=['',3.4840606,-0.7695558,5.84];
ybs[4999]=['',3.4585571,1.4025054,6.25];
ybs[5000]=['',3.480115,0.3433461,6.45];
ybs[5001]=['59 Vir',3.4812808,0.1625114,5.22];
ybs[5002]=['',3.4946919,-1.2592203,6.04];
ybs[5003]=['',3.4833426,0.2367136,5.33];
ybs[5004]=['',3.4845599,-0.0137792,6.37];
ybs[5005]=['σ Vir',3.4849502,0.0934959,4.8];
ybs[5006]=['',3.490141,-0.897077,6.19];
ybs[5007]=['20 CVn',3.484144,0.7061545,4.73];
ybs[5008]=['',3.4783616,1.191973,6.2];
ybs[5009]=['61 Vir',3.4887414,-0.3215603,4.74];
ybs[5010]=['γ Hya',3.4910663,-0.4063865,3];
ybs[5011]=['',3.4904111,0.0623989,6.62];
ybs[5012]=['',3.4882849,0.593157,5.82];
ybs[5013]=['21 CVn',3.4869672,0.8651464,5.15];
ybs[5014]=['',3.4992525,-1.0452006,6.18];
ybs[5015]=['',3.4909097,0.6111357,6.02];
ybs[5016]=['',3.4991715,-0.9225863,5.48];
ybs[5017]=['',3.5000514,-0.9758619,6.02];
ybs[5018]=['ι Cen',3.4986117,-0.6427085,2.75];
ybs[5019]=['',3.5004419,-0.8201781,5.77];
ybs[5020]=['',3.5103424,-1.2611482,6.05];
ybs[5021]=['',3.4984552,0.0493826,6.26];
ybs[5022]=['23 CVn',3.496241,0.6987988,5.6];
ybs[5023]=['',3.5022642,-0.3421018,6.21];
ybs[5024]=['',3.5081308,-1.0661186,6.18];
ybs[5025]=['',3.5082923,-1.0663997,4.53];
ybs[5026]=['',3.5063261,-0.9127199,5.83];
ybs[5027]=['',3.5028292,0.0344729,5.69];
ybs[5028]=['',3.5088526,-0.8387161,6.16];
ybs[5029]=['',3.5095925,-0.8495317,6.38];
ybs[5030]=['64 Vir',3.5048329,0.0880124,5.87];
ybs[5031]=['',3.5145457,-1.1283107,4.53];
ybs[5032]=['ι1 Mus',3.5206456,-1.308982,5.05];
ybs[5033]=['',3.5096787,-0.5812259,6.22];
ybs[5034]=['63 Vir',3.5088757,-0.3114907,5.37];
ybs[5035]=['',3.5037625,0.7642978,6.35];
ybs[5036]=['',3.5132536,-0.871525,6.48];
ybs[5037]=['65 Vir',3.5099961,-0.0878986,5.89];
ybs[5038]=['',3.5199017,-1.1274242,5.31];
ybs[5039]=['',3.5231263,-1.234624,5.67];
ybs[5040]=['66 Vir',3.5154035,-0.0920736,5.75];
ybs[5041]=['ι2 Mus',3.5302164,-1.3055567,6.63];
ybs[5042]=['',3.5119174,0.6444142,6.07];
ybs[5043]=['',3.5149717,0.2150315,6.44];
ybs[5044]=['ζ UMa',3.5115243,0.9566777,2.27];
ybs[5045]=['ζ UMa',3.5115897,0.9566148,3.95];
ybs[5046]=['α Vir',3.5182782,-0.1967475,0.98];
ybs[5047]=['',3.5174342,0.4143936,5.78];
ybs[5048]=['',3.5228573,-0.6958017,5.09];
ybs[5049]=['',3.5224924,-0.0227543,5.97];
ybs[5050]=['',3.5264166,-0.7262162,5.69];
ybs[5051]=['',3.5273722,-0.8596605,6.31];
ybs[5052]=['80 UMa',3.5171785,0.9577777,4.01];
ybs[5053]=['',3.5284367,-0.8637951,6.28];
ybs[5054]=['68 Vir',3.524964,-0.223732,5.25];
ybs[5055]=['',3.5277397,-0.702915,6.4];
ybs[5056]=['',3.535897,-1.2171704,6.2];
ybs[5057]=['',3.5220551,0.8013998,5.88];
ybs[5058]=['69 Vir',3.5282128,-0.280729,4.76];
ybs[5059]=['',3.5369805,-1.1307367,6.11];
ybs[5060]=['',3.5201169,1.1021722,6.5];
ybs[5061]=['',3.537565,-0.8949324,5.06];
ybs[5062]=['70 Vir',3.5320563,0.2385533,4.98];
ybs[5063]=['',3.5198204,1.2615257,5.79];
ybs[5064]=['',3.5247003,1.1279097,6.66];
ybs[5065]=['',3.5251433,1.1276289,7.04];
ybs[5066]=['',3.5293051,0.918653,6.34];
ybs[5067]=['',3.5315994,0.7089339,6.47];
ybs[5068]=['',3.5358399,-0.0257448,6.43];
ybs[5069]=['',3.530267,0.8809788,6.8];
ybs[5070]=['',3.538195,-0.4082659,4.97];
ybs[5071]=['71 Vir',3.5355303,0.1868846,5.65];
ybs[5072]=['',3.5570711,-1.3557378,6.48];
ybs[5073]=['',3.5327293,0.8832691,6.43];
ybs[5074]=['κ Oct',3.5992369,-1.4950795,5.58];
ybs[5075]=['',3.5309821,1.0443181,5.4];
ybs[5076]=['',3.5390066,0.123367,6.17];
ybs[5077]=['',3.5388414,0.1030242,6.51];
ybs[5078]=['72 Vir',3.5410612,-0.1148543,6.09];
ybs[5079]=['',3.5443292,-0.6897148,3.88];
ybs[5080]=['',3.5463226,-0.492583,6.47];
ybs[5081]=['',3.5219323,1.3706546,5.77];
ybs[5082]=['',3.5488649,-0.6721124,6.16];
ybs[5083]=['',3.5566365,-1.1474177,6.37];
ybs[5084]=['73 Vir',3.5483095,-0.3288016,6.01];
ybs[5085]=['74 Vir',3.5477645,-0.111106,4.69];
ybs[5086]=['',3.5438779,0.7329663,6.08];
ybs[5087]=['',3.5508905,-0.5027021,5.69];
ybs[5088]=['',3.5508051,-0.5179302,6.45];
ybs[5089]=['75 Vir',3.5518146,-0.2700538,5.55];
ybs[5090]=['76 Vir',3.5522011,-0.1793302,5.21];
ybs[5091]=['',3.5523462,-0.1274938,6.68];
ybs[5092]=['',3.550951,0.4230063,6.11];
ybs[5093]=['',3.5596013,-0.8444209,6.33];
ybs[5094]=['',3.5602797,-0.5832946,6.44];
ybs[5095]=['78 Vir',3.5570748,0.0619463,4.94];
ybs[5096]=['',3.5596973,-0.2325469,5.91];
ybs[5097]=['ζ Vir',3.559588,-0.0123106,3.37];
ybs[5098]=['',3.5574647,0.6750859,6.37];
ybs[5099]=['81 UMa',3.5558748,0.9641016,5.6];
ybs[5100]=['',3.5593905,0.6470458,4.98];
ybs[5101]=['80 Vir',3.5632738,-0.0960881,5.73];
ybs[5102]=['24 CVn',3.5575728,0.8535799,4.7];
ybs[5103]=['',3.5721414,-1.0786288,5.63];
ybs[5104]=['',3.5631858,0.1761978,6.49];
ybs[5105]=['',3.5827133,-1.3228259,6.34];
ybs[5106]=['',3.5611196,0.7694726,6.84];
ybs[5107]=['',3.569545,-0.6034792,6.5];
ybs[5108]=['',3.5709198,-0.7723484,5.98];
ybs[5109]=['',3.5798135,-1.2313921,6.1];
ybs[5110]=['',3.5692299,-0.4643282,5.78];
ybs[5111]=['',3.5722794,-0.8122281,5.9];
ybs[5112]=['',3.5759825,-1.0214319,6.42];
ybs[5113]=['',3.5691836,0.4276808,5.74];
ybs[5114]=['',3.5789663,-1.0076072,6.01];
ybs[5115]=['',3.5853384,-1.2373796,6.59];
ybs[5116]=['',3.5671709,0.8618009,6.49];
ybs[5117]=['25 CVn',3.5710156,0.6315661,4.82];
ybs[5118]=['',3.5775593,-0.5178299,5.83];
ybs[5119]=['',3.5743595,0.2477127,6.52];
ybs[5120]=['',3.5854042,-1.12897,5.79];
ybs[5121]=['',3.5561679,1.3340787,6.57];
ybs[5122]=['ε Cen',3.5834391,-0.9350556,2.3];
ybs[5123]=['',3.5717116,0.8832385,6.48];
ybs[5124]=['',3.5837701,-0.8736876,6];
ybs[5125]=['',3.5820746,-0.6956265,6.27];
ybs[5126]=['',3.5826507,-0.7009299,5.6];
ybs[5127]=['',3.5782477,0.3168941,6.48];
ybs[5128]=['',3.580719,0.1856621,5.57];
ybs[5129]=['',3.5679135,1.2415082,5.5];
ybs[5130]=['',3.5930372,-1.0279133,5.38];
ybs[5131]=['',3.5916184,-0.9541355,5.01];
ybs[5132]=['82 UMa',3.5794153,0.9217588,5.46];
ybs[5133]=['',3.5833342,0.53937,6.21];
ybs[5134]=['1 Boo',3.5853551,0.3464013,5.75];
ybs[5135]=['',3.5851018,0.4879426,6.23];
ybs[5136]=['',3.589735,-0.4111601,6.59];
ybs[5137]=['',3.5910206,-0.5882614,6.05];
ybs[5138]=['',3.5833423,0.8798404,6.32];
ybs[5139]=['2 Boo',3.5868938,0.390739,5.62];
ybs[5140]=['82 Vir',3.5899114,-0.1537819,5.01];
ybs[5141]=['',3.5969356,-0.9926686,6];
ybs[5142]=['',3.5965508,-0.8883417,6.41];
ybs[5143]=['',3.5829049,0.9965687,6.29];
ybs[5144]=['83 UMa',3.5846978,0.9524862,4.66];
ybs[5145]=['',3.5962651,-0.7244651,5.98];
ybs[5146]=['',3.5922506,0.1445214,6.16];
ybs[5147]=['',3.5995623,-0.7360928,5.98];
ybs[5148]=['',3.6024858,-0.8922197,6.47];
ybs[5149]=['84 Vir',3.5960365,0.0598715,5.36];
ybs[5150]=['',3.5927482,0.7254695,6.3];
ybs[5151]=['',3.5939807,0.6087905,5.98];
ybs[5152]=['',3.5873652,1.1294798,5.85];
ybs[5153]=['',3.5998612,-0.0978495,6.51];
ybs[5154]=['',3.5987271,0.394318,6.13];
ybs[5155]=['83 Vir',3.6026227,-0.284253,5.6];
ybs[5156]=['',3.603953,-0.4469457,6.21];
ybs[5157]=['',3.607692,-0.4576808,5.81];
ybs[5158]=['1 Cen',3.608155,-0.5785931,4.23];
ybs[5159]=['',3.5986248,0.9068196,6.02];
ybs[5160]=['85 Vir',3.6073667,-0.2770637,6.19];
ybs[5161]=['',3.6158866,-1.0942632,6.51];
ybs[5162]=['',3.612956,-0.8995354,4.65];
ybs[5163]=['86 Vir',3.6088511,-0.2187537,5.51];
ybs[5164]=['',3.6137087,-0.634579,5.15];
ybs[5165]=['',3.6164274,-0.8788791,5.91];
ybs[5166]=['',3.6172175,-0.8801292,5.45];
ybs[5167]=['',3.6041413,0.9734091,6.5];
ybs[5168]=['',3.614421,-0.1713191,6.05];
ybs[5169]=['',3.6090786,0.7152647,5.87];
ybs[5170]=['',3.6095495,0.6701532,5.94];
ybs[5171]=['87 Vir',3.6154301,-0.3135772,5.43];
ybs[5172]=['3 Boo',3.6116095,0.4467238,5.95];
ybs[5173]=['',3.6129612,0.1089746,6.33];
ybs[5174]=['',3.5900548,1.3605988,5.91];
ybs[5175]=['τ Boo',3.6141214,0.302814,4.5];
ybs[5176]=['',3.6125111,0.6708348,5.5];
ybs[5177]=['84 UMa',3.6101968,0.9481657,5.7];
ybs[5178]=['',3.658993,-1.444618,5.95];
ybs[5179]=['',3.6223386,-0.6250053,6.53];
ybs[5180]=['ν Cen',3.6250683,-0.7294414,3.41];
ybs[5181]=['η UMa',3.6145317,0.8588185,1.86];
ybs[5182]=['2 Cen',3.6246082,-0.6031332,4.19];
ybs[5183]=['μ Cen',3.6255804,-0.7431611,3.04];
ybs[5184]=['',3.6368072,-1.2131245,5.75];
ybs[5185]=['',3.6198755,0.5425163,5.62];
ybs[5186]=['89 Vir',3.6261244,-0.318352,4.97];
ybs[5187]=['',3.6273799,-0.509416,6.18];
ybs[5188]=['',3.6285917,-0.6982548,6.44];
ybs[5189]=['',3.6210084,0.6882962,7.4];
ybs[5190]=['υ Boo',3.6238108,0.2738701,4.07];
ybs[5191]=['6 Boo',3.6247404,0.3692776,4.91];
ybs[5192]=['',3.6292288,-0.3491201,6.53];
ybs[5193]=['',3.5861119,1.4424185,5.98];
ybs[5194]=['',3.6245583,0.6375105,6.38];
ybs[5195]=['',3.6280698,0.0940956,6.01];
ybs[5196]=['',3.635216,-0.8203875,5.77];
ybs[5197]=['',3.6367454,-0.9235787,5.25];
ybs[5198]=['',3.634118,-0.6377251,6.35];
ybs[5199]=['',3.6326588,-0.4275451,6.45];
ybs[5200]=['3 Cen',3.634977,-0.5777043,4.56];
ybs[5201]=['3 Cen',3.6350135,-0.5777091,6.06];
ybs[5202]=['',3.6357666,-0.5537052,6.12];
ybs[5203]=['',3.6235383,1.0713358,5.96];
ybs[5204]=['',3.6303798,0.6050481,6.65];
ybs[5205]=['',3.6307239,0.6031625,5.87];
ybs[5206]=['',3.6267964,1.0198565,6.46];
ybs[5207]=['',3.6439387,-0.9333793,5.89];
ybs[5208]=['',3.6498389,-1.1825876,5.71];
ybs[5209]=['',3.6335111,0.5993207,4.74];
ybs[5210]=['',3.6362191,0.2104832,6.04];
ybs[5211]=['4 Cen',3.6409869,-0.5590815,4.73];
ybs[5212]=['',3.6425594,-0.6242922,5.54];
ybs[5213]=['',3.6446753,-0.8243729,6.1];
ybs[5214]=['',3.6439628,-0.618187,6.19];
ybs[5215]=['7 Boo',3.6400735,0.3111489,5.7];
ybs[5216]=['10 Dra',3.6305453,1.1277897,4.65];
ybs[5217]=['',3.6282265,1.1904788,6.4];
ybs[5218]=['',3.6455747,-0.5004677,6.04];
ybs[5219]=['',3.6396634,0.4981655,5.9];
ybs[5220]=['',3.6503754,-0.9122107,5.71];
ybs[5221]=['ζ Cen',3.6516401,-0.8271633,2.55];
ybs[5222]=['90 Vir',3.6469085,-0.0280637,5.15];
ybs[5223]=['',3.6481958,-0.1424833,6.19];
ybs[5224]=['',3.6554011,-0.9466032,6.14];
ybs[5225]=['η Boo',3.6464748,0.3192711,2.68];
ybs[5226]=['',3.6563938,-0.9565942,6];
ybs[5227]=['',3.6520375,-0.5478519,6.51];
ybs[5228]=['86 UMa',3.6418066,0.9359063,5.7];
ybs[5229]=['',3.6550454,-0.8150152,5.83];
ybs[5230]=['',3.6776514,-1.3734557,6.09];
ybs[5231]=['',3.6618304,-1.1133583,4.71];
ybs[5232]=['',3.66586,-1.1502485,6.2];
ybs[5233]=['',3.6515657,0.2435047,6.16];
ybs[5234]=['92 Vir',3.6545574,0.016513,5.91];
ybs[5235]=['',3.6526741,0.5572483,6.32];
ybs[5236]=['',3.6593625,-0.4036413,6.14];
ybs[5237]=['9 Boo',3.654503,0.4780024,5.01];
ybs[5238]=['φ Cen',3.6633962,-0.7366123,3.83];
ybs[5239]=['υ1 Cen',3.6652741,-0.7837827,3.87];
ybs[5240]=['47 Hya',3.6640189,-0.4376607,5.15];
ybs[5241]=['',3.6681578,-0.8809316,5.91];
ybs[5242]=['',3.6731995,-1.0748571,6.49];
ybs[5243]=['',3.6762079,-1.1584069,5.97];
ybs[5244]=['',3.6639041,0.2538681,6];
ybs[5245]=['10 Boo',3.6636899,0.3768555,5.76];
ybs[5246]=['',3.6573558,1.0714326,6.37];
ybs[5247]=['48 Hya',3.6704975,-0.4383181,5.77];
ybs[5248]=['',3.6692873,-0.0637618,6.4];
ybs[5249]=['',3.6766428,-0.7038104,6.13];
ybs[5250]=['υ2 Cen',3.678607,-0.7977314,4.34];
ybs[5251]=['θ Aps',3.6979412,-1.342134,5.5];
ybs[5252]=['',3.6756816,0.1534416,5.99];
ybs[5253]=['11 Boo',3.6745797,0.4761859,6.23];
ybs[5254]=['τ Vir',3.6771568,0.0251565,4.26];
ybs[5255]=['',3.6809398,-0.4805392,5.48];
ybs[5256]=['',3.686608,-0.9829025,5.92];
ybs[5257]=['β Cen',3.688586,-1.0554965,0.61];
ybs[5258]=['',3.6838769,-0.5547805,6.18];
ybs[5259]=['',3.6860413,-0.7247637,6.11];
ybs[5260]=['',3.6808573,0.1672643,6.2];
ybs[5261]=['',3.6785217,0.7967541,6.27];
ybs[5262]=['',3.6873976,-0.3931203,6.3];
ybs[5263]=['',3.6852343,0.1864724,6.3];
ybs[5264]=['',3.685624,0.1299193,6.26];
ybs[5265]=['',3.6870558,0.0837472,6.24];
ybs[5266]=['',3.6886175,-0.0957099,6.39];
ybs[5267]=['',3.6897072,-0.2630907,6.28];
ybs[5268]=['',3.696746,-0.9559406,6.17];
ybs[5269]=['',3.7110223,-1.3081483,6.02];
ybs[5270]=['',3.6817486,0.8878348,6.15];
ybs[5271]=['',3.6998935,-1.0440087,6.42];
ybs[5272]=['',3.6753316,1.1968683,6.34];
ybs[5273]=['',3.6901373,0.0383137,6.28];
ybs[5274]=['',3.6931549,-0.2868962,6.56];
ybs[5275]=['χ Cen',3.6973445,-0.7204996,4.36];
ybs[5276]=['',3.6980013,-0.7538735,6.2];
ybs[5277]=['π Hya',3.6983548,-0.467474,3.27];
ybs[5278]=['θ Cen',3.6999741,-0.6365511,2.06];
ybs[5279]=['',3.7081591,-1.1049553,6.4];
ybs[5280]=['95 Vir',3.6994736,-0.1643234,5.46];
ybs[5281]=['α Dra',3.686898,1.1217828,3.65];
ybs[5282]=['',3.7108891,-1.0363364,6.34];
ybs[5283]=['',3.7190552,-1.2288183,6.05];
ybs[5284]=['',3.709743,-0.7604782,6.17];
ybs[5285]=['',3.7212204,-1.2185911,6.06];
ybs[5286]=['',3.7132315,-0.9006875,6];
ybs[5287]=['',3.7147776,-0.9344482,4.75];
ybs[5288]=['96 Vir',3.7095182,-0.1821341,6.47];
ybs[5289]=['',3.7035363,0.7636345,5.27];
ybs[5290]=['13 Boo',3.704866,0.8614376,5.25];
ybs[5291]=['',3.7176377,-0.2862775,4.91];
ybs[5292]=['',3.7063734,1.0338733,6.46];
ybs[5293]=['η Aps',3.7571407,-1.4155655,4.91];
ybs[5294]=['12 Boo',3.7148432,0.4361746,4.83];
ybs[5295]=['3 UMi',3.6963045,1.3001278,6.45];
ybs[5296]=['',3.7492412,-1.3572118,6.47];
ybs[5297]=['',3.7202389,0.0220236,6.43];
ybs[5298]=['',3.729486,-0.9383875,5.56];
ybs[5299]=['',3.7246611,-0.4269819,6.34];
ybs[5300]=['',3.7183778,0.5619104,6.11];
ybs[5301]=['',3.7312556,-0.9551406,6.11];
ybs[5302]=['50 Hya',3.7262976,-0.4775413,5.08];
ybs[5303]=['',3.7234613,0.0403047,5.01];
ybs[5304]=['',3.728261,-0.4662138,6.24];
ybs[5305]=['κ Vir',3.7264812,-0.181053,4.19];
ybs[5306]=['',3.737036,-0.9980691,5.07];
ybs[5307]=['',3.7297078,-0.0164986,5.91];
ybs[5308]=['',3.7352275,-0.7319372,5.61];
ybs[5309]=['',3.7385741,-0.9356523,6.39];
ybs[5310]=['',3.7453505,-1.1639002,5.75];
ybs[5311]=['4 UMi',3.7035867,1.3516915,4.82];
ybs[5312]=['',3.7327569,-0.1055409,6.36];
ybs[5313]=['14 Boo',3.7311988,0.224446,5.54];
ybs[5314]=['',3.7361988,-0.5128,6.08];
ybs[5315]=['',3.7394438,-0.7871429,6.31];
ybs[5316]=['',3.7443379,-1.0474196,6.39];
ybs[5317]=['',3.7863904,-1.4476599,6.42];
ybs[5318]=['κ1 Boo',3.7272766,0.9021244,6.69];
ybs[5319]=['κ2 Boo',3.7273708,0.9021682,4.54];
ybs[5320]=['15 Boo',3.7345741,0.1745529,5.29];
ybs[5321]=['',3.7348735,0.0564912,6.45];
ybs[5322]=['',3.7375872,-0.3193965,5.43];
ybs[5323]=['',3.7336024,0.3800256,6.39];
ybs[5324]=['',3.7196205,1.210075,5.24];
ybs[5325]=['',3.7317704,0.7229034,6.24];
ybs[5326]=['ε Aps',3.7747099,-1.3998554,5.06];
ybs[5327]=['',3.7419153,-0.5818988,6.55];
ybs[5328]=['ι Vir',3.7400024,-0.1064585,4.08];
ybs[5329]=['δ Oct',3.7989062,-1.4619424,4.32];
ybs[5330]=['α Boo',3.7379307,0.3330667,-0.04];
ybs[5331]=['',3.7415136,-0.1173019,6.44];
ybs[5332]=['',3.7420683,-0.057514,6.15];
ybs[5333]=['',3.7397175,0.3283467,5.98];
ybs[5334]=['',3.7448505,-0.3260977,6.22];
ybs[5335]=['',3.7351601,0.9151896,6.58];
ybs[5336]=['',3.7417682,0.3494579,6.25];
ybs[5337]=['',3.7405976,0.6919487,6.38];
ybs[5338]=['',3.7510516,-0.5815244,6.54];
ybs[5339]=['',3.7588476,-1.0711244,5.23];
ybs[5340]=['ι Boo',3.7390683,0.8947981,4.75];
ybs[5341]=['λ Boo',3.7402677,0.8026653,4.18];
ybs[5342]=['',3.7459195,0.2646737,5.8];
ybs[5343]=['',3.7487413,-0.1333602,6.47];
ybs[5344]=['ι Lup',3.7559087,-0.8055707,3.55];
ybs[5345]=['',3.7517298,-0.3283731,5.9];
ybs[5346]=['',3.7535404,-0.4522796,5.87];
ybs[5347]=['',3.7555457,-0.6475458,5.94];
ybs[5348]=['',3.7605022,-0.9858386,4.33];
ybs[5349]=['λ Vir',3.7536679,-0.2350828,4.52];
ybs[5350]=['',3.7442382,0.8937569,6.2];
ybs[5351]=['',3.7476796,0.6180374,4.81];
ybs[5352]=['',3.7590989,-0.7532263,5.56];
ybs[5353]=['',3.7464401,0.8360666,6.32];
ybs[5354]=['',3.7615754,-0.7903697,4.77];
ybs[5355]=['18 Boo',3.7538127,0.2252531,5.41];
ybs[5356]=['υ Vir',3.7553145,-0.0412522,5.14];
ybs[5357]=['ψ Cen',3.7606359,-0.6629277,4.05];
ybs[5358]=['',3.75587,0.0049949,6.19];
ybs[5359]=['',3.7516331,0.6749059,6.86];
ybs[5360]=['20 Boo',3.7558397,0.2829,4.86];
ybs[5361]=['',3.7706806,-1.0220031,4.92];
ybs[5362]=['',3.7509052,0.9558453,6.53];
ybs[5363]=['',3.7554085,0.675371,6.33];
ybs[5364]=['',3.7571935,0.529381,6.44];
ybs[5365]=['',3.7701832,-0.8450415,6.09];
ybs[5366]=['',3.7682746,-0.6088423,5.56];
ybs[5367]=['',3.7733395,-0.8878322,6.02];
ybs[5368]=['',3.7715285,-0.6913101,4.42];
ybs[5369]=['',3.7826547,-1.1919114,5.61];
ybs[5370]=['',3.775524,-0.9297951,6];
ybs[5371]=['51 Hya',3.7714241,-0.4860883,4.77];
ybs[5372]=['',3.7847842,-1.1566189,6.36];
ybs[5373]=['2 Lib',3.7724733,-0.2061409,6.21];
ybs[5374]=['',3.771443,0.01998,6.27];
ybs[5375]=['',3.7718183,0.1457025,6.86];
ybs[5376]=['',3.7718256,0.1457315,5.12];
ybs[5377]=['',3.7702717,0.4405402,6.22];
ybs[5378]=['',3.7746073,0.142191,5.95];
ybs[5379]=['',3.8046485,-1.3408296,6.07];
ybs[5380]=['',3.7788273,-0.4346356,5.32];
ybs[5381]=['',3.791132,-1.1504733,5.85];
ybs[5382]=['',3.7754163,0.099892,5.1];
ybs[5383]=['',3.7779502,-0.2053584,6.49];
ybs[5384]=['',3.7758758,0.1394243,6.19];
ybs[5385]=['τ1 Lup',3.7853028,-0.7909369,4.56];
ybs[5386]=['τ2 Lup',3.7854991,-0.7936953,4.35];
ybs[5387]=['',3.7817027,-0.3502162,6.61];
ybs[5388]=['',3.7855598,-0.7402832,6.32];
ybs[5389]=['',3.7831834,-0.4703416,6.48];
ybs[5390]=['',3.7881238,-0.6976017,6.35];
ybs[5391]=['',3.790006,-0.806867,5.83];
ybs[5392]=['',3.7802236,0.6684055,6.27];
ybs[5393]=['',3.7974674,-1.0348562,6.45];
ybs[5394]=['θ Boo',3.7783655,0.9032861,4.05];
ybs[5395]=['22 Boo',3.7850038,0.3338994,5.39];
ybs[5396]=['104 Vir',3.7897195,-0.1084876,6.17];
ybs[5397]=['52 Hya',3.7936447,-0.5163907,4.97];
ybs[5398]=['',3.8096118,-1.1835337,5.83];
ybs[5399]=['φ Vir',3.7931027,-0.0405512,4.81];
ybs[5400]=['106 Vir',3.7953594,-0.1220989,5.42];
ybs[5401]=['',3.7887033,0.7143522,6.63];
ybs[5402]=['',3.802819,-0.792665,5.5];
ybs[5403]=['',3.8039207,-0.8659239,5.37];
ybs[5404]=['',3.7937841,0.4920762,7.62];
ybs[5405]=['',3.7939148,0.4921054,7.12];
ybs[5406]=['',3.7924429,0.6300914,6.1];
ybs[5407]=['',3.8061194,-0.7145281,6.39];
ybs[5408]=['',3.800191,0.0128116,5.94];
ybs[5409]=['',3.8070855,-0.6800517,5.97];
ybs[5410]=['24 Boo',3.7934252,0.8682916,5.59];
ybs[5411]=['',3.8140086,-0.9945178,6.93];
ybs[5412]=['',3.7993621,0.5532038,6.06];
ybs[5413]=['',3.7980768,0.7278177,6.35];
ybs[5414]=['',3.8040949,0.0816408,6.02];
ybs[5415]=['σ Lup',3.8138815,-0.8822784,4.42];
ybs[5416]=['',3.8182089,-0.9615352,5.87];
ybs[5417]=['',3.8178747,-0.921073,5.87];
ybs[5418]=['',3.8154605,-0.5376994,6.09];
ybs[5419]=['ρ Boo',3.8081305,0.5284361,3.58];
ybs[5420]=['5 UMi',3.7852,1.3194753,4.25];
ybs[5421]=['',3.8201166,-0.7364091,6.6];
ybs[5422]=['',3.8262508,-1.0490969,6.4];
ybs[5423]=['',3.810447,0.4639637,6.01];
ybs[5424]=['26 Boo',3.811461,0.3868699,5.92];
ybs[5425]=['γ Boo',3.8089504,0.6669632,3.03];
ybs[5426]=['',3.8017261,1.1011443,6.09];
ybs[5427]=['',3.8061312,1.049488,6.27];
ybs[5428]=['',3.8225292,-0.3583575,6.5];
ybs[5429]=['',3.8261737,-0.7262345,5.87];
ybs[5430]=['η Cen',3.8261201,-0.7374144,2.31];
ybs[5431]=['',3.8144898,0.6434278,6.43];
ybs[5432]=['',3.8100032,0.9652322,5.76];
ybs[5433]=['',3.838107,-1.1872488,6.04];
ybs[5434]=['',3.8298438,-0.80875,5.55];
ybs[5435]=['',3.8183723,0.566202,6.33];
ybs[5436]=['',3.8299346,-0.6927192,6.13];
ybs[5437]=['σ Boo',3.8205732,0.5175198,4.46];
ybs[5438]=['',3.8201815,0.6376126,6.03];
ybs[5439]=['',3.8314153,-0.7034414,5.74];
ybs[5440]=['',3.8342892,-0.8068,5.41];
ybs[5441]=['',3.8174984,0.9943454,6.48];
ybs[5442]=['',3.8197106,0.860011,5.74];
ybs[5443]=['ρ Lup',3.8368707,-0.8642519,4.05];
ybs[5444]=['',3.8270155,0.4041738,6.38];
ybs[5445]=['',3.8317182,-0.216382,6.2];
ybs[5446]=['',3.8383068,-0.6786921,6.02];
ybs[5447]=['',3.8423898,-0.814648,6.07];
ybs[5448]=['',3.8435168,-0.8577805,6.39];
ybs[5449]=['α1 Cen',3.845196,-1.0633737,-0.01];
ybs[5450]=['α2 Cen',3.8452105,-1.0633785,1.33];
ybs[5451]=['',3.8489595,-0.9866709,6.3];
ybs[5452]=['',3.8363839,0.3177585,5.91];
ybs[5453]=['α Cir',3.8584029,-1.1356128,3.19];
ybs[5454]=['',3.8354454,0.7600872,5.7];
ybs[5455]=['',3.8551828,-1.0246283,6.22];
ybs[5456]=['',3.8499883,-0.6322651,5.67];
ybs[5457]=['',3.835078,0.9412765,5.85];
ybs[5458]=['33 Boo',3.8381486,0.773399,5.39];
ybs[5459]=['α Lup',3.8544358,-0.8286671,2.3];
ybs[5460]=['α Aps',3.8861235,-1.3811356,3.83];
ybs[5461]=['',3.8541419,-0.6612078,4];
ybs[5462]=['',3.8455821,0.3819505,6.1];
ybs[5463]=['',3.8472917,0.2346226,5.91];
ybs[5464]=['',3.8534285,-0.5414741,6.37];
ybs[5465]=['π1 Boo',3.847306,0.2849609,4.94];
ybs[5466]=['π2 Boo',3.8473279,0.2849512,5.88];
ybs[5467]=['ζ Boo',3.8492132,0.238014,4.83];
ybs[5468]=['ζ Boo',3.8492132,0.238014,4.43];
ybs[5469]=['',3.8097034,1.3886956,6.26];
ybs[5470]=['31 Boo',3.8515237,0.1408605,4.86];
ybs[5471]=['32 Boo',3.85178,0.201928,5.56];
ybs[5472]=['',3.8703167,-1.098954,5.36];
ybs[5473]=['',3.8523133,0.3670904,6.38];
ybs[5474]=['4 Lib',3.8592513,-0.4378663,5.73];
ybs[5475]=['',3.8614567,-0.6154701,4.05];
ybs[5476]=['',3.8682798,-1.022196,6.11];
ybs[5477]=['μ Vir',3.8580208,-0.1003354,3.88];
ybs[5478]=['',3.8691666,-0.9720019,6.1];
ybs[5479]=['',3.867267,-0.6157821,4.92];
ybs[5480]=['34 Boo',3.8587871,0.4614198,4.81];
ybs[5481]=['',4.1080371,-1.5386313,6.48];
ybs[5482]=['',3.8510181,1.0676359,6.25];
ybs[5483]=['',3.8596833,0.70457,5.73];
ybs[5484]=['',3.8743475,-0.8295609,5.74];
ybs[5485]=['',3.8769854,-0.9158203,5.21];
ybs[5486]=['',3.8672352,-0.026311,6.07];
ybs[5487]=['54 Hya',3.8713816,-0.4456259,4.94];
ybs[5488]=['',3.8777752,-0.9127115,6.07];
ybs[5489]=['',3.8718023,-0.4056573,5.81];
ybs[5490]=['',3.8859301,-1.1638249,5.91];
ybs[5491]=['108 Vir',3.868529,0.0109536,5.69];
ybs[5492]=['ο Boo',3.8669829,0.2945192,4.6];
ybs[5493]=['5 Lib',3.8709375,-0.2713844,6.33];
ybs[5494]=['',3.8720447,-0.3711527,6.4];
ybs[5495]=['ε Boo',3.8655865,0.47098,5.12];
ybs[5496]=['ε Boo',3.8655866,0.4709655,2.7];
ybs[5497]=['',3.8673824,0.328035,6.13];
ybs[5498]=['',3.8765507,-0.669855,5.94];
ybs[5499]=['',3.8787429,-0.7617728,6.3];
ybs[5500]=['',3.8664566,0.5706978,6.28];
ybs[5501]=['109 Vir',3.8717511,0.0314754,3.72];
ybs[5502]=['',3.8707744,0.2625413,5.63];
ybs[5503]=['',3.8766185,-0.3737401,6.06];
ybs[5504]=['55 Hya',3.8773838,-0.4487834,5.63];
ybs[5505]=['',3.8864528,-0.9905801,6.23];
ybs[5506]=['56 Hya',3.8790207,-0.456863,5.24];
ybs[5507]=['57 Hya',3.879962,-0.4666161,5.77];
ybs[5508]=['',3.8793887,-0.2256448,6.35];
ybs[5509]=['',3.8832574,-0.640941,6.04];
ybs[5510]=['',3.9070321,-1.2789198,5.6];
ybs[5511]=['',3.8857992,-0.424812,5.68];
ybs[5512]=['',3.8833952,-0.0163401,6.14];
ybs[5513]=['μ Lib',3.8855482,-0.2484855,5.31];
ybs[5514]=['',3.880499,0.4237314,6.14];
ybs[5515]=['π1 Oct',3.9525826,-1.454051,5.65];
ybs[5516]=['58 Hya',3.8901777,-0.4895334,4.41];
ybs[5517]=['',3.9023598,-1.1152133,5.87];
ybs[5518]=['ο Lup',3.8966814,-0.7620626,4.32];
ybs[5519]=['',3.8831985,0.6583853,6.16];
ybs[5520]=['α1 Lib',3.8915701,-0.2807365,5.15];
ybs[5521]=['α2 Lib',3.8924079,-0.281511,2.75];
ybs[5522]=['',3.8872791,0.4979029,5.8];
ybs[5523]=['38 Boo',3.8836992,0.8033359,5.74];
ybs[5524]=['',3.888688,0.4158065,5.85];
ybs[5525]=['11 Lib',3.8926599,-0.0416585,4.94];
ybs[5526]=['',3.8925438,-0.00602,6.18];
ybs[5527]=['',3.8843886,0.8951172,6.51];
ybs[5528]=['39 Boo',3.8852002,0.8487942,5.69];
ybs[5529]=['ζ Cir',3.9119617,-1.1532721,6.09];
ybs[5530]=['',3.9287842,-1.3395006,5.34];
ybs[5531]=['',3.8892412,0.6489835,5.48];
ybs[5532]=['',3.9001474,-0.5351888,6.29];
ybs[5533]=['',3.9017233,-0.6613109,5.03];
ybs[5534]=['ξ Boo',3.8937315,0.3318487,4.55];
ybs[5535]=['π2 Oct',3.9651725,-1.4507253,5.65];
ybs[5536]=['',3.9150328,-1.0506902,5.2];
ybs[5537]=['',3.9392026,-1.3481684,5.93];
ybs[5538]=['12 Lib',3.9077446,-0.4315972,5.3];
ybs[5539]=['',3.9093342,-0.5827115,5.82];
ybs[5540]=['',3.9025471,0.2725783,6.4];
ybs[5541]=['θ Cir',3.9203925,-1.0972247,5.11];
ybs[5542]=['',3.8920572,1.033344,5.46];
ybs[5543]=['',3.9024812,0.3327629,6.01];
ybs[5544]=['ξ1 Lib',3.9075898,-0.2091744,5.8];
ybs[5545]=['',3.937077,-1.3110342,6.2];
ybs[5546]=['',3.9175284,-0.9231995,5.38];
ybs[5547]=['ω Oct',3.9976606,-1.4812076,5.91];
ybs[5548]=['',3.9142099,-0.592396,5.32];
ybs[5549]=['',3.9182892,-0.8371438,5.64];
ybs[5550]=['',3.9206492,-0.89941,6.64];
ybs[5551]=['',3.9181505,-0.6894355,6.36];
ybs[5552]=['',3.9175285,-0.5711126,6.06];
ybs[5553]=['β UMi',3.8862672,1.2927226,2.08];
ybs[5554]=['ξ2 Lib',3.9179939,-0.2006314,5.46];
ybs[5555]=['',3.9205159,-0.5103901,6.29];
ybs[5556]=['',3.9253977,-0.8543055,6.35];
ybs[5557]=['',3.9149317,0.2506389,5.77];
ybs[5558]=['',3.921321,-0.3752615,5.74];
ybs[5559]=['',3.9133334,0.5622461,6.12];
ybs[5560]=['16 Lib',3.9196211,-0.0773505,4.49];
ybs[5561]=['β Lup',3.9267735,-0.7543105,2.68];
ybs[5562]=['',3.9269819,-0.6979844,6.15];
ybs[5563]=['',3.921136,-0.0044128,5.53];
ybs[5564]=['',3.9183947,0.3747176,6.49];
ybs[5565]=['',3.9191309,0.2845384,5.71];
ybs[5566]=['κ Cen',3.9294771,-0.7363344,3.13];
ybs[5567]=['59 Hya',3.9266961,-0.4841914,5.65];
ybs[5568]=['17 Lib',3.9243376,-0.1961764,6.6];
ybs[5569]=['',3.9295969,-0.6626328,6.47];
ybs[5570]=['',3.9307947,-0.7547603,6.1];
ybs[5571]=['',3.9143082,0.8646846,5.63];
ybs[5572]=['18 Lib',3.9272617,-0.195983,5.87];
ybs[5573]=['',3.9270417,-0.0885582,6.09];
ybs[5574]=['',3.9289986,0.0782451,5.93];
ybs[5575]=['',3.9382755,-0.6657081,5.89];
ybs[5576]=['δ Lib',3.9362698,-0.1501499,4.92];
ybs[5577]=['',3.9414105,-0.6011358,6.22];
ybs[5578]=['40 Boo',3.9289131,0.6838314,5.64];
ybs[5579]=['',3.9179858,1.1492476,4.6];
ybs[5580]=['',3.9376748,-0.0495488,5.52];
ybs[5581]=['60 Hya',3.9417892,-0.4912085,5.85];
ybs[5582]=['',3.9350091,0.383299,6.38];
ybs[5583]=['η Cir',3.9558707,-1.1190033,5.17];
ybs[5584]=['',3.939723,-0.0039151,5.71];
ybs[5585]=['',3.9458008,-0.5711871,5.44];
ybs[5586]=['',3.8789919,1.4385628,5.64];
ybs[5587]=['',3.933009,0.8236822,6.37];
ybs[5588]=['',3.9676666,-1.256407,6.52];
ybs[5589]=['',3.9438712,-0.0543634,6.61];
ybs[5590]=['ω Boo',3.9403028,0.4350124,4.81];
ybs[5591]=['110 Vir',3.9443989,0.0350468,4.4];
ybs[5592]=['β Boo',3.9390251,0.7034862,3.5];
ybs[5593]=['σ Lib',3.9502708,-0.4426998,3.29];
ybs[5594]=['',3.9536816,-0.7146077,6.41];
ybs[5595]=['π Lup',3.9557684,-0.8226358,4.72];
ybs[5596]=['π Lup',3.9557684,-0.8226358,4.82];
ybs[5597]=['',3.9563269,-0.7181961,5.15];
ybs[5598]=['',3.9355384,1.0492996,5.93];
ybs[5599]=['',3.9442759,0.6130034,5.51];
ybs[5600]=['',3.9495738,0.0944152,6.5];
ybs[5601]=['',3.969809,-1.1406923,6.17];
ybs[5602]=['',3.94389,0.777738,6.65];
ybs[5603]=['',3.9465092,0.6018415,6.59];
ybs[5604]=['',3.9578193,-0.4515508,6.67];
ybs[5605]=['',3.9601072,-0.6343611,6.27];
ybs[5606]=['ψ Boo',3.9504284,0.4688774,4.54];
ybs[5607]=['',3.9659926,-0.8581765,5.77];
ybs[5608]=['44 Boo',3.9466883,0.8302768,4.76];
ybs[5609]=['',3.9613141,-0.5410664,5.96];
ybs[5610]=['',3.9605668,-0.3859609,6.17];
ybs[5611]=['',3.9769006,-1.1722478,5.76];
ybs[5612]=['ν Lib',3.9611585,-0.2851671,5.2];
ybs[5613]=['',3.9760415,-1.1121851,6.28];
ybs[5614]=['',3.9689009,-0.7097412,5.79];
ybs[5615]=['',3.970977,-0.7495995,5.85];
ybs[5616]=['λ Lup',3.9719401,-0.7916944,4.05];
ybs[5617]=['47 Boo',3.9538038,0.838956,5.57];
ybs[5618]=['',3.9916685,-1.2714673,6.01];
ybs[5619]=['',3.9457343,1.1490658,6.13];
ybs[5620]=['',3.9594089,0.6348381,6.35];
ybs[5621]=['',3.9651268,0.0945358,6.16];
ybs[5622]=['',3.9815262,-1.0734253,6.3];
ybs[5623]=['',3.9633178,0.3204419,6.02];
ybs[5624]=['45 Boo',3.9629522,0.4326225,4.93];
ybs[5625]=['',3.9570625,0.9507544,5.25];
ybs[5626]=['',3.9772215,-0.6784628,5.98];
ybs[5627]=['',3.983209,-0.9673693,5.54];
ybs[5628]=['46 Boo',3.9676798,0.457622,5.67];
ybs[5629]=['',3.9702286,0.2295789,6.1];
ybs[5630]=['',3.9685693,0.4368102,5.81];
ybs[5631]=['',3.9775514,-0.4609987,5.76];
ybs[5632]=['',3.9839052,-0.7916373,6.44];
ybs[5633]=['',3.9836868,-0.7916715,7.39];
ybs[5634]=['',3.9986239,-1.2244919,5.81];
ybs[5635]=['',3.9915097,-1.0790191,6.32];
ybs[5636]=['κ1 Lup',3.985652,-0.8520279,3.87];
ybs[5637]=['κ2 Lup',3.9857616,-0.8521295,5.69];
ybs[5638]=['',3.9662969,0.8722043,6.39];
ybs[5639]=['ζ Lup',3.9874077,-0.9106926,3.41];
ybs[5640]=['',3.9881836,-0.8429627,6.33];
ybs[5641]=['',3.9892956,-0.7780685,4.82];
ybs[5642]=['ι1 Lib',3.9856877,-0.3468219,4.54];
ybs[5643]=['',3.9902017,-0.6312992,6.1];
ybs[5644]=['',3.9839377,0.3297967,5.89];
ybs[5645]=['',3.9904928,-0.4204093,6.47];
ybs[5646]=['ι2 Lib',3.9904781,-0.3442983,6.08];
ybs[5647]=['23 Lib',3.9913507,-0.4431118,6.45];
ybs[5648]=['',3.9931708,-0.4585455,5.84];
ybs[5649]=['',3.9867776,0.3352117,6.68];
ybs[5650]=['1 Lup',3.9965637,-0.5514887,4.91];
ybs[5651]=['',4.0071128,-1.0643337,5.73];
ybs[5652]=['26 Lib',3.995841,-0.3114971,6.17];
ybs[5653]=['',4.0029235,-0.8404186,5.95];
ybs[5654]=['δ Cir',4.0085989,-1.0652671,5.09];
ybs[5655]=['',3.9901796,0.3997506,6.3];
ybs[5656]=['ε Cir',4.0120181,-1.1115664,4.86];
ybs[5657]=['',4.0033035,-0.7255217,5.16];
ybs[5658]=['',4.0038782,-0.760316,6.04];
ybs[5659]=['',3.996725,-0.0974162,6.28];
ybs[5660]=['β Cir',4.0108459,-1.0276272,4.07];
ybs[5661]=['γ TrA',4.0184168,-1.2000257,2.89];
ybs[5662]=['',3.9748632,1.1816028,6.17];
ybs[5663]=['',3.9898786,0.6664611,6.2];
ybs[5664]=['',3.9923588,0.5534256,5.99];
ybs[5665]=['3 Ser',3.9979428,0.0848371,5.33];
ybs[5666]=['χ Boo',3.994134,0.5076328,5.26];
ybs[5667]=['',3.9922161,0.7346492,6.13];
ybs[5668]=['',4.0039383,-0.3923077,5.5];
ybs[5669]=['4 Ser',4.0008102,0.0051285,5.63];
ybs[5670]=['',4.0167242,-1.0572062,5.46];
ybs[5671]=['δ Boo',3.9984105,0.5800806,3.47];
ybs[5672]=['',4.0124125,-0.7179974,6.28];
ybs[5673]=['μ Lup',4.0144331,-0.8369242,4.27];
ybs[5674]=['',4.0259161,-1.1791035,6.28];
ybs[5675]=['β Lib',4.0062677,-0.1651248,2.61];
ybs[5676]=['2 Lup',4.01053,-0.5275508,4.34];
ybs[5677]=['',4.0158212,-0.713236,5.59];
ybs[5678]=['',4.0143109,-0.546055,6.18];
ybs[5679]=['',4.0182129,-0.6488003,6.2];
ybs[5680]=['',4.012258,-0.0094027,5.89];
ybs[5681]=['',3.991918,1.174042,5.13];
ybs[5682]=['',4.0115234,0.3577121,5.7];
ybs[5683]=['',3.9925888,1.2019443,6.51];
ybs[5684]=['5 Ser',4.016021,0.0294661,5.06];
ybs[5685]=['δ Lup',4.0264366,-0.710761,3.22];
ybs[5686]=['',4.0273881,-0.7125435,6.2];
ybs[5687]=['',4.0268896,-0.6683777,6.48];
ybs[5688]=['ν1 Lup',4.0301895,-0.83782,5];
ybs[5689]=['ν2 Lup',4.0287504,-0.8446292,5.65];
ybs[5690]=['',4.0358047,-1.0599777,5.67];
ybs[5691]=['28 Lib',4.0234995,-0.3182599,6.17];
ybs[5692]=['',4.0158799,0.5661554,6.32];
ybs[5693]=['ο Lib',4.0239795,-0.2727011,6.3];
ybs[5694]=['γ Cir',4.0365491,-1.0366569,4.51];
ybs[5695]=['φ1 Lup',4.0281347,-0.634206,3.56];
ybs[5696]=['',4.0225594,-0.043454,6.35];
ybs[5697]=['',4.0241467,-0.1029916,5.54];
ybs[5698]=['ε Lup',4.0323712,-0.7812967,3.37];
ybs[5699]=['ο CrB',4.0187732,0.51556,5.51];
ybs[5700]=['6 Ser',4.0235551,0.0111522,5.35];
ybs[5701]=['',4.0231827,0.4342636,6.39];
ybs[5702]=['φ2 Lup',4.034057,-0.6446199,4.54];
ybs[5703]=['',4.0504505,-1.1935113,5.89];
ybs[5704]=['11 UMi',4.0015934,1.2522003,5.02];
ybs[5705]=['',4.017398,0.9055088,5.66];
ybs[5706]=['',4.0205176,0.7741871,6.19];
ybs[5707]=['7 Ser',4.0291198,0.2180218,6.28];
ybs[5708]=['50 Boo',4.0259194,0.5734777,5.37];
ybs[5709]=['υ Lup',4.041141,-0.6943795,5.37];
ybs[5710]=['',4.0363091,-0.2171989,5.72];
ybs[5711]=['8 Ser',4.0353634,-0.0191587,6.12];
ybs[5712]=['',4.042648,-0.6674844,7.03];
ybs[5713]=['ε Lib',4.0376817,-0.181466,4.94];
ybs[5714]=['',4.0436643,-0.6773292,4.6];
ybs[5715]=['',4.0554836,-1.1275678,5.71];
ybs[5716]=['',4.029176,0.689504,5.5];
ybs[5717]=['η CrB',4.0321041,0.5273042,5.58];
ybs[5718]=['η CrB',4.0321041,0.5273042,6.08];
ybs[5719]=['ρ Oct',4.1387252,-1.475351,5.57];
ybs[5720]=['κ1 Aps',4.0749205,-1.2821395,5.49];
ybs[5721]=['',4.0274627,1.0816053,5.98];
ybs[5722]=['',4.0352502,0.7888184,6.01];
ybs[5723]=['μ1 Boo',4.0374142,0.6510473,4.31];
ybs[5724]=['μ2 Boo',4.0375245,0.6505335,6.5];
ybs[5725]=['γ UMi',4.0173643,1.2524001,3.05];
ybs[5726]=['',4.0521639,-0.643005,5.45];
ybs[5727]=['',4.0273695,1.1041928,5.79];
ybs[5728]=['',4.0580383,-0.9018235,6.1];
ybs[5729]=['τ1 Ser',4.0438785,0.267972,5.17];
ybs[5730]=['',4.0441802,0.3387069,6.27];
ybs[5731]=['',4.0453927,0.5979779,5.46];
ybs[5732]=['',4.0618664,-0.8169114,5.24];
ybs[5733]=['ζ1 Lib',4.0555889,-0.293036,5.64];
ybs[5734]=['ι Dra',4.0378711,1.0278464,3.29];
ybs[5735]=['',4.0516673,0.4368213,6.02];
ybs[5736]=['10 Ser',4.0566978,0.0308751,5.17];
ybs[5737]=['β CrB',4.0523069,0.5067084,3.68];
ybs[5738]=['',4.0453534,0.9415324,6.45];
ybs[5739]=['',4.0659694,-0.3630406,6.22];
ybs[5740]=['ζ3 Lib',4.0661342,-0.291152,5.82];
ybs[5741]=['',4.0730833,-0.6753511,6.25];
ybs[5742]=['',4.0554259,0.8225413,6.15];
ybs[5743]=['',4.0718018,-0.5751372,6.46];
ybs[5744]=['',4.0494496,1.0856264,6.5];
ybs[5745]=['',4.0504203,1.0576104,5.9];
ybs[5746]=['',4.0708327,-0.3531955,6.22];
ybs[5747]=['',4.1113364,-1.3611191,6.18];
ybs[5748]=['',4.0664772,0.1484734,6.57];
ybs[5749]=['',4.0557278,0.9620571,6.43];
ybs[5750]=['',4.0633409,0.5447797,6.46];
ybs[5751]=['',4.0634855,0.6410885,6.37];
ybs[5752]=['',4.0746942,-0.3445643,5.52];
ybs[5753]=['ν1 Boo',4.0653197,0.7114091,5.02];
ybs[5754]=['ζ4 Lib',4.0759486,-0.2953826,5.5];
ybs[5755]=['',4.0772161,-0.4286656,7];
ybs[5756]=['',4.0938931,-1.1463873,6.51];
ybs[5757]=['',4.0816859,-0.7005228,5.82];
ybs[5758]=['',4.0567563,1.0825696,6.38];
ybs[5759]=['',4.0674536,0.6378177,6.38];
ybs[5760]=['τ2 Ser',4.0716506,0.2789795,6.22];
ybs[5761]=['ε TrA',4.0959006,-1.1586643,4.11];
ybs[5762]=['11 Ser',4.0756717,-0.0219526,5.51];
ybs[5763]=['',4.0830552,-0.6880123,6.36];
ybs[5764]=['ν2 Boo',4.0690357,0.712574,5.02];
ybs[5765]=['36 Lib',4.0837641,-0.4907448,5.15];
ybs[5766]=['γ Lup',4.0866124,-0.7197277,2.78];
ybs[5767]=['37 Lib',4.0812299,-0.1768947,4.62];
ybs[5768]=['θ CrB',4.0744632,0.5460734,4.14];
ybs[5769]=['',4.0818335,-0.1006324,6.51];
ybs[5770]=['',4.0823597,-0.1615145,5.17];
ybs[5771]=['',4.0900771,-0.7858991,4.54];
ybs[5772]=['κ2 Aps',4.1136711,-1.2830727,5.65];
ybs[5773]=['',4.079111,0.2978707,6.45];
ybs[5774]=['',4.0914204,-0.7760939,5.43];
ybs[5775]=['',4.0633125,1.1193878,5.79];
ybs[5776]=['',4.1216571,-1.3290542,5.95];
ybs[5777]=['γ Lib',4.0872607,-0.2593515,3.91];
ybs[5778]=['δ Ser',4.0833372,0.1826811,3.8];
ybs[5779]=['δ Ser',4.0833371,0.1827102,3.8];
ybs[5780]=['',4.0908135,-0.5787994,6.24];
ybs[5781]=['',4.0848015,0.0278969,6.56];
ybs[5782]=['',4.1119563,-1.2268945,6.44];
ybs[5783]=['α CrB',4.0823083,0.4650256,2.23];
ybs[5784]=['υ Lib',4.0942609,-0.4922639,3.58];
ybs[5785]=['τ3 Ser',4.0863927,0.30692,6.12];
ybs[5786]=['',4.0880626,0.195396,6.07];
ybs[5787]=['ω Lup',4.0994057,-0.7441503,4.33];
ybs[5788]=['',4.1034251,-0.9152784,5.44];
ybs[5789]=['14 Ser',4.0913469,-0.0110226,6.51];
ybs[5790]=['μ CrB',4.0842402,0.6796223,5.11];
ybs[5791]=['',4.0961559,-0.4598847,6.19];
ybs[5792]=['16 Ser',4.0907329,0.173487,5.26];
ybs[5793]=['',4.109042,-1.0467897,5.95];
ybs[5794]=['τ5 Ser',4.0905205,0.280107,5.93];
ybs[5795]=['',4.1013755,-0.6846893,6.57];
ybs[5796]=['',4.0974689,-0.4051082,5.78];
ybs[5797]=['',4.1020726,-0.684116,6.04];
ybs[5798]=['',4.0867628,0.6685244,6.42];
ybs[5799]=['',4.099663,-0.4935055,6.32];
ybs[5800]=['',4.099447,-0.3680068,5.84];
ybs[5801]=['',4.0833872,0.9398843,5.97];
ybs[5802]=['τ Lib',4.1014489,-0.5209236,3.66];
ybs[5803]=['',4.0917837,0.5222255,6.52];
ybs[5804]=['41 Lib',4.1021681,-0.3380843,5.38];
ybs[5805]=['',4.1007793,-0.1546959,6.5];
ybs[5806]=['',4.1007865,-0.1546377,6.48];
ybs[5807]=['',4.0870082,0.9075629,6.74];
ybs[5808]=['',4.0862881,0.9522569,5.74];
ybs[5809]=['',4.1042504,-0.4052469,6.34];
ybs[5810]=['ψ1 Lup',4.1064891,-0.6017966,4.67];
ybs[5811]=['',4.1124667,-0.8343278,6.23];
ybs[5812]=['',4.1084931,-0.5459716,6.34];
ybs[5813]=['φ Boo',4.0954046,0.7030869,5.24];
ybs[5814]=['42 Lib',4.1083195,-0.4168949,4.96];
ybs[5815]=['',4.11322,-0.7806672,4.64];
ybs[5816]=['θ UMi',4.0614828,1.3487388,4.96];
ybs[5817]=['',4.0999817,0.6039889,6.11];
ybs[5818]=['',4.0931883,0.9501451,5.87];
ybs[5819]=['',4.0493694,1.4028113,6.58];
ybs[5820]=['',4.096971,0.8155667,5.75];
ybs[5821]=['',4.1067231,0.2091723,6.25];
ybs[5822]=['',4.1197896,-0.8649265,6.04];
ybs[5823]=['ζ1 CrB',4.1023127,0.6382306,6];
ybs[5824]=['ζ2 CrB',4.1023491,0.638216,5.07];
ybs[5825]=['',4.0980389,0.8788466,5.84];
ybs[5826]=['',4.126475,-1.0533725,6.48];
ybs[5827]=['',4.1191681,-0.6543627,5.24];
ybs[5828]=['κ Lib',4.1154437,-0.3446404,4.74];
ybs[5829]=['ψ2 Lup',4.1192391,-0.6069865,4.75];
ybs[5830]=['τ6 Ser',4.1101444,0.278497,6.01];
ybs[5831]=['',4.0999374,1.0097695,6.45];
ybs[5832]=['ι Ser',4.1124879,0.3421281,4.52];
ybs[5833]=['χ Ser',4.1137502,0.2230501,5.33];
ybs[5834]=['',4.0916009,1.2080068,5.62];
ybs[5835]=['τ7 Ser',4.1140987,0.3210754,5.81];
ybs[5836]=['',4.1269851,-0.7310468,5.94];
ybs[5837]=['',4.1216992,-0.2637238,6.31];
ybs[5838]=['η Lib',4.1245931,-0.2747047,5.41];
ybs[5839]=['γ CrB',4.1174494,0.4577697,3.84];
ybs[5840]=['',4.119773,0.2373771,6.48];
ybs[5841]=['',4.1446227,-1.1433181,6.18];
ybs[5842]=['',4.14463,-1.1433181,6.39];
ybs[5843]=['ψ Ser',4.1238408,0.0427311,5.88];
ybs[5844]=['α Ser',4.1247615,0.1109849,2.65];
ybs[5845]=['π CrB',4.1226302,0.5663433,5.56];
ybs[5846]=['',4.1343749,-0.4909149,6.51];
ybs[5847]=['',4.1165081,0.9126942,5.51];
ybs[5848]=['τ8 Ser',4.1263077,0.3001574,6.14];
ybs[5849]=['',4.1296923,0.0939135,5.58];
ybs[5850]=['',4.136938,-0.606466,5.61];
ybs[5851]=['',4.1310094,0.0144062,6.33];
ybs[5852]=['',4.1402048,-0.7026572,6.42];
ybs[5853]=['25 Ser',4.132975,-0.0326416,5.4];
ybs[5854]=['',4.1403511,-0.6629021,6.01];
ybs[5855]=['',4.1471878,-0.9163418,6.07];
ybs[5856]=['',4.1360019,-0.1079619,6.24];
ybs[5857]=['β Ser',4.1328531,0.2680159,3.67];
ybs[5858]=['λ Ser',4.1342218,0.1271893,4.43];
ybs[5859]=['',4.1528318,-0.9297953,5.77];
ybs[5860]=['υ Ser',4.1376939,0.2452187,5.71];
ybs[5861]=['',4.1518103,-0.8547913,5.84];
ybs[5862]=['',4.1529443,-0.7935232,6.12];
ybs[5863]=['',4.157361,-0.962013,5.73];
ybs[5864]=['',4.1417758,0.2395244,6];
ybs[5865]=['',4.1454874,-0.0677735,5.53];
ybs[5866]=['',4.1666665,-1.1382176,6.54];
ybs[5867]=['',4.1403002,0.5527605,6.44];
ybs[5868]=['',4.1325071,0.9670699,5.92];
ybs[5869]=['κ Ser',4.1438923,0.3155034,4.09];
ybs[5870]=['',4.1428009,0.4902965,5.85];
ybs[5871]=['μ Ser',4.1484136,-0.0609907,3.53];
ybs[5872]=['',4.1585266,-0.8224664,6.01];
ybs[5873]=['χ Lup',4.1553252,-0.5880155,3.95];
ybs[5874]=['',4.1681581,-1.0937863,6.19];
ybs[5875]=['1 Sco',4.1550888,-0.4505564,4.64];
ybs[5876]=['',4.1320566,1.0914197,5.19];
ybs[5877]=['',4.1371048,0.9653666,5.86];
ybs[5878]=['ω Ser',4.1511709,0.0372182,5.23];
ybs[5879]=['δ CrB',4.147337,0.4538561,4.63];
ybs[5880]=['',4.164607,-0.8844974,6.6];
ybs[5881]=['κ TrA',4.1785184,-1.1984207,5.09];
ybs[5882]=['ε Ser',4.153391,0.0770399,3.71];
ybs[5883]=['',4.1606445,-0.5227209,6.4];
ybs[5884]=['',4.1525206,0.263018,5.2];
ybs[5885]=['36 Ser',4.1555595,-0.0550488,5.11];
ybs[5886]=['',4.1575683,-0.2477832,6.19];
ybs[5887]=['β TrA',4.1759863,-1.1081472,2.85];
ybs[5888]=['',4.1744425,-1.0612488,6.15];
ybs[5889]=['ρ Ser',4.154807,0.3650222,4.76];
ybs[5890]=['',4.1772621,-1.0513729,5.77];
ybs[5891]=['κ CrB',4.1540709,0.6212309,4.82];
ybs[5892]=['λ Lib',4.165172,-0.3530763,5.03];
ybs[5893]=['ζ UMi',4.1160627,1.3565987,4.32];
ybs[5894]=['2 Sco',4.1665687,-0.4431329,4.59];
ybs[5895]=['',4.1797565,-1.0566916,5.76];
ybs[5896]=['',4.1677828,-0.4292699,5.39];
ybs[5897]=['',4.1679075,-0.4195831,5.42];
ybs[5898]=['θ Lib',4.1671919,-0.293072,4.15];
ybs[5899]=['',4.1622223,0.3026543,6.36];
ybs[5900]=['',4.1705249,-0.4782313,6.14];
ybs[5901]=['39 Ser',4.1635187,0.2292315,6.1];
ybs[5902]=['3 Sco',4.1711348,-0.4416655,5.87];
ybs[5903]=['',4.1650806,0.2794707,6.09];
ybs[5904]=['χ Her',4.1600135,0.7398224,4.62];
ybs[5905]=['47 Lib',4.172434,-0.339377,5.94];
ybs[5906]=['',4.1750737,-0.5435859,6.21];
ybs[5907]=['4 Sco',4.1748558,-0.4595001,5.62];
ybs[5908]=['',4.1781179,-0.6968305,6.03];
ybs[5909]=['40 Ser',4.1700805,0.1486718,6.29];
ybs[5910]=['',4.1929883,-1.1361682,5.75];
ybs[5911]=['',4.1827966,-0.8416509,6.31];
ybs[5912]=['',4.1572887,0.9732561,5.81];
ybs[5913]=['',4.1782694,-0.5558362,6.29];
ybs[5914]=['',4.1692745,0.3534075,5.44];
ybs[5915]=['ξ1 Lup',4.1812483,-0.5938889,5.12];
ybs[5916]=['ξ2 Lup',4.1812991,-0.59385,5.62];
ybs[5917]=['',4.177654,-0.2523871,6.37];
ybs[5918]=['ρ Sco',4.1810133,-0.5109472,3.88];
ybs[5919]=['',4.1833756,-0.6326072,5.8];
ybs[5920]=['',4.1790504,-0.2598847,6.13];
ybs[5921]=['',4.1740735,0.3239152,6.26];
ybs[5922]=['2 Her',4.1685028,0.751827,5.37];
ybs[5923]=['γ Ser',4.1776208,0.272279,3.85];
ybs[5924]=['',4.1841383,-0.3672813,5.85];
ybs[5925]=['',4.1884914,-0.6556025,6.31];
ybs[5926]=['λ CrB',4.1738478,0.6612246,5.45];
ybs[5927]=['',4.1956614,-0.943885,6.1];
ybs[5928]=['4 Her',4.1723646,0.7418418,5.75];
ybs[5929]=['',4.2024374,-1.1141403,6.41];
ybs[5930]=['φ Ser',4.1811088,0.2505172,5.54];
ybs[5931]=['48 Lib',4.1861608,-0.2502774,4.88];
ybs[5932]=['',4.1882391,-0.4344401,5.43];
ybs[5933]=['',4.1930556,-0.7296206,4.99];
ybs[5934]=['π Sco',4.1894747,-0.4568267,2.89];
ybs[5935]=['',4.1950108,-0.7105638,6.49];
ybs[5936]=['',4.2009859,-0.953591,6.13];
ybs[5937]=['ε CrB',4.1821642,0.4680457,4.15];
ybs[5938]=['η Lup',4.1955675,-0.671191,3.41];
ybs[5939]=['',4.1724294,1.0271269,6.31];
ybs[5940]=['',4.1811906,0.691752,6.31];
ybs[5941]=['',4.2096547,-1.0925722,6.25];
ybs[5942]=['',4.199053,-0.7067604,6.21];
ybs[5943]=['δ Sco',4.1958095,-0.3958595,2.32];
ybs[5944]=['49 Lib',4.195563,-0.2895982,5.47];
ybs[5945]=['',4.2252863,-1.2646205,5.7];
ybs[5946]=['',4.2005123,-0.5576045,6.33];
ybs[5947]=['',4.1877255,0.6385069,5.62];
ybs[5948]=['',4.1905532,0.4513493,6.33];
ybs[5949]=['50 Lib',4.1973293,-0.1478401,5.55];
ybs[5950]=['',4.181411,0.954503,4.95];
ybs[5951]=['ι1 Nor',4.211913,-1.0093783,4.63];
ybs[5952]=['η Nor',4.2097476,-0.8602334,4.65];
ybs[5953]=['',4.1971748,0.0762409,5.83];
ybs[5954]=['',4.1874136,0.8695402,6.05];
ybs[5955]=['',4.2062082,-0.5095344,6.03];
ybs[5956]=['5 Her',4.1984186,0.3099577,5.12];
ybs[5957]=['',4.2099088,-0.6747525,4.89];
ybs[5958]=['ρ CrB',4.1969664,0.5802247,5.41];
ybs[5959]=['',4.2090733,-0.4524471,5];
ybs[5960]=['',4.2103263,-0.5595259,6.01];
ybs[5961]=['ι CrB',4.1988541,0.5199705,4.99];
ybs[5962]=['π Ser',4.2028473,0.3969899,4.83];
ybs[5963]=['',4.2115099,-0.4325653,6.21];
ybs[5964]=['',4.2135465,-0.5807064,6.1];
ybs[5965]=['',4.215147,-0.6618372,5.9];
ybs[5966]=['43 Ser',4.2098531,0.0860233,6.08];
ybs[5967]=['ξ Sco',4.2130281,-0.1995024,5.07];
ybs[5968]=['ξ Sco',4.2130281,-0.1995024,4.77];
ybs[5969]=['',4.228646,-0.9817033,6.16];
ybs[5970]=['δ Nor',4.2237706,-0.7894102,4.72];
ybs[5971]=['',4.2003197,0.9225297,5.93];
ybs[5972]=['υ Her',4.2039197,0.8024717,4.76];
ybs[5973]=['',4.2067429,0.6383281,5.83];
ybs[5974]=['β1 Sco',4.2179825,-0.3466685,2.62];
ybs[5975]=['β2 Sco',4.2180043,-0.3466054,4.92];
ybs[5976]=['θ Dra',4.1988657,1.021129,4.01];
ybs[5977]=['θ Lup',4.2237391,-0.6433063,4.23];
ybs[5978]=['',4.2210412,-0.413,5.92];
ybs[5979]=['',4.2188546,-0.1108046,6.53];
ybs[5980]=['',4.2199624,-0.1081458,6.41];
ybs[5981]=['',4.2267132,-0.6424864,5.73];
ybs[5982]=['',4.2179035,0.1403082,6.29];
ybs[5983]=['ω1 Sco',4.223994,-0.3617303,3.96];
ybs[5984]=['ι2 Nor',4.2371833,-1.0121098,5.57];
ybs[5985]=['',4.2042821,1.0358967,6.19];
ybs[5986]=['',4.2248517,-0.2465659,6.32];
ybs[5987]=['ω2 Sco',4.226613,-0.3652065,4.32];
ybs[5988]=['',4.228757,-0.4279181,6.33];
ybs[5989]=['',4.2324973,-0.6834862,7.05];
ybs[5990]=['',4.2325111,-0.6832728,6.65];
ybs[5991]=['',4.2299731,-0.4604615,5.38];
ybs[5992]=['11 Sco',4.227207,-0.223431,5.78];
ybs[5993]=['',4.2324956,-0.4143606,5.88];
ybs[5994]=['45 Ser',4.2265478,0.1716625,5.63];
ybs[5995]=['',4.2250187,0.3798924,6.14];
ybs[5996]=['',4.2363627,-0.5708081,6.19];
ybs[5997]=['',4.2379249,-0.5864452,5.54];
ybs[5998]=['κ Her',4.2282667,0.296549,5];
ybs[5999]=['κ Her',4.2282956,0.29668,6.25];
ybs[6000]=['47 Ser',4.2302646,0.1479765,5.73];
ybs[6001]=['',4.2326773,0.0593229,5.91];
ybs[6002]=['',4.2374985,-0.3210682,6.47];
ybs[6003]=['8 Her',4.2313218,0.2993278,6.14];
ybs[6004]=['',4.2334766,0.1103657,5.97];
ybs[6005]=['',4.2445086,-0.7186227,5.86];
ybs[6006]=['',4.2366539,-0.061471,5.37];
ybs[6007]=['',4.2427903,-0.5143636,5.13];
ybs[6008]=['τ CrB',4.2313772,0.6359154,4.76];
ybs[6009]=['ζ Nor',4.2546936,-0.9703006,5.81];
ybs[6010]=['δ1 Aps',4.2919087,-1.3743678,4.68];
ybs[6011]=['δ2 Aps',4.292322,-1.3738676,5.27];
ybs[6012]=['',4.2540889,-0.9376784,5.83];
ybs[6013]=['φ Her',4.2300213,0.7832919,4.26];
ybs[6014]=['κ Nor',4.2550468,-0.9544125,4.94];
ybs[6015]=['',4.2167123,1.1825187,5.44];
ybs[6016]=['ν Sco',4.2465156,-0.3404052,6.3];
ybs[6017]=['ν Sco',4.246596,-0.3405942,4.01];
ybs[6018]=['13 Sco',4.2482716,-0.488348,4.59];
ybs[6019]=['12 Sco',4.2481249,-0.4969198,5.67];
ybs[6020]=['δ TrA',4.2647652,-1.1124349,3.85];
ybs[6021]=['ψ Sco',4.2462899,-0.1765964,4.94];
ybs[6022]=['',4.243437,0.1685667,6.53];
ybs[6023]=['16 Sco',4.2467702,-0.1501247,5.43];
ybs[6024]=['',4.2012114,1.3392821,5.56];
ybs[6025]=['',4.2431258,0.2899201,6.08];
ybs[6026]=['',4.2301289,1.0102344,6.33];
ybs[6027]=['',4.2728358,-1.1866988,5.75];
ybs[6028]=['',4.2320342,0.9734306,6.49];
ybs[6029]=['10 Her',4.24355,0.4091126,5.7];
ybs[6030]=['',4.2656647,-1.0116689,5.63];
ybs[6031]=['',4.2502057,-0.0746036,6.25];
ybs[6032]=['',4.2544956,-0.4271721,6.41];
ybs[6033]=['',4.2432583,0.5809886,6.29];
ybs[6034]=['',4.2575252,-0.5770809,5.92];
ybs[6035]=['θ Nor',4.2622028,-0.8277168,5.14];
ybs[6036]=['',4.2437127,0.6347893,5.63];
ybs[6037]=['9 Her',4.2512781,0.086701,5.48];
ybs[6038]=['χ Sco',4.2544162,-0.2075319,5.22];
ybs[6039]=['',4.2625405,-0.749656,6.14];
ybs[6040]=['',4.2433593,0.7386264,5.87];
ybs[6041]=['',4.2575036,-0.3693185,6.41];
ybs[6042]=['',4.2483161,0.4645552,6.5];
ybs[6043]=['',4.2581578,-0.3244233,6.32];
ybs[6044]=['',4.2594776,-0.4455761,6.05];
ybs[6045]=['',4.2691233,-0.9400841,5.44];
ybs[6046]=['δ Oph',4.256312,-0.0654049,2.74];
ybs[6047]=['',4.2554739,0.1020824,6.31];
ybs[6048]=['γ1 Nor',4.270094,-0.8747581,4.99];
ybs[6049]=['',4.2718082,-0.9274349,6.33];
ybs[6050]=['18 Sco',4.2620379,-0.1469885,5.5];
ybs[6051]=['',4.2632925,-0.2600788,6.09];
ybs[6052]=['',4.2789321,-1.0114257,6.49];
ybs[6053]=['σ CrB',4.2563958,0.5900209,5.64];
ybs[6054]=['σ CrB',4.2563958,0.5900161,6.66];
ybs[6055]=['16 Her',4.2604968,0.3273511,5.69];
ybs[6056]=['',4.2684232,-0.3727255,6.61];
ybs[6057]=['',4.267564,-0.0699023,6.18];
ybs[6058]=['',4.2615216,0.4776942,6.14];
ybs[6059]=['',4.2433687,1.1709421,6.21];
ybs[6060]=['',4.2744716,-0.500298,4.78];
ybs[6061]=['λ Nor',4.2795358,-0.7456824,5.45];
ybs[6062]=['γ2 Nor',4.2824405,-0.8762571,4.02];
ybs[6063]=['',4.2854251,-0.9632468,5.77];
ybs[6064]=['υ CrB',4.2656206,0.5078622,5.78];
ybs[6065]=['ε Oph',4.2736976,-0.0827916,3.24];
ybs[6066]=['',4.2777593,-0.3537517,6.29];
ybs[6067]=['',4.2800146,-0.5403041,5.49];
ybs[6068]=['',4.2770338,-0.2604602,5.94];
ybs[6069]=['19 UMi',4.2334723,1.3233519,5.48];
ybs[6070]=['',4.2847978,-0.6890703,6.12];
ybs[6071]=['ο Sco',4.2844956,-0.4227087,4.55];
ybs[6072]=['20 UMi',4.2412585,1.3117256,6.39];
ybs[6073]=['',4.2938666,-0.8660543,5.33];
ybs[6074]=['σ Sco',4.286961,-0.4475458,2.89];
ybs[6075]=['',4.2935293,-0.7672689,5.88];
ybs[6076]=['',4.2656606,1.0420175,5.4];
ybs[6077]=['',4.2804458,0.3679533,6.05];
ybs[6078]=['',4.2508568,1.2800554,5.98];
ybs[6079]=['',4.308068,-1.1025643,6.15];
ybs[6080]=['',4.2751391,0.8549884,5.91];
ybs[6081]=['',4.2789365,0.6921655,5.46];
ybs[6082]=['τ Her',4.2777496,0.8074378,3.89];
ybs[6083]=['σ Ser',4.289878,0.0171011,4.82];
ybs[6084]=['',4.2999723,-0.6848914,5.4];
ybs[6085]=['γ Her',4.2885824,0.3334208,3.75];
ybs[6086]=['',4.2924903,-0.0371543,6.23];
ybs[6087]=['',4.2993296,-0.580284,6.47];
ybs[6088]=['ζ TrA',4.3231451,-1.2240069,4.91];
ybs[6089]=['',4.3041993,-0.792328,6.33];
ybs[6090]=['',4.302096,-0.6564868,5.42];
ybs[6091]=['',4.2680262,1.1956027,6.41];
ybs[6092]=['γ Aps',4.3494954,-1.3777723,3.89];
ybs[6093]=['ξ CrB',4.2888758,0.5383041,4.85];
ybs[6094]=['ψ Oph',4.2994668,-0.3505639,4.5];
ybs[6095]=['',4.3022905,-0.5192596,6.63];
ybs[6096]=['',4.3022977,-0.5192838,5.84];
ybs[6097]=['ν1 CrB',4.2898744,0.5890467,5.2];
ybs[6098]=['ν2 CrB',4.2904458,0.5873801,5.39];
ybs[6099]=['ι TrA',4.319539,-1.1188322,5.27];
ybs[6100]=['',4.2924948,0.5634632,6.4];
ybs[6101]=['21 Her',4.2988714,0.1204226,5.85];
ybs[6102]=['ρ Oph',4.3060734,-0.4100623,5.02];
ybs[6103]=['ρ Oph',4.3060661,-0.4100429,5.92];
ybs[6104]=['',4.3200367,-1.0235646,5.69];
ybs[6105]=['ε Nor',4.3143195,-0.8308078,4.47];
ybs[6106]=['η UMi',4.2625772,1.3212731,4.95];
ybs[6107]=['ω Her',4.3040215,0.2440941,4.57];
ybs[6108]=['χ Oph',4.312152,-0.322944,4.42];
ybs[6109]=['',4.3054951,0.3289055,6.7];
ybs[6110]=['',4.3264731,-1.0088235,6.06];
ybs[6111]=['',4.3074947,0.1982712,6.11];
ybs[6112]=['',4.3183028,-0.6497119,5.79];
ybs[6113]=['25 Her',4.3029787,0.6518117,5.54];
ybs[6114]=['',4.3106101,0.0401501,6.07];
ybs[6115]=['',4.3316799,-1.0764938,5.2];
ybs[6116]=['',4.2837975,1.2053192,5.25];
ybs[6117]=['',4.2974104,0.9626648,5.74];
ybs[6118]=['',4.3148288,-0.1334248,5.23];
ybs[6119]=['υ Oph',4.31519,-0.1469261,4.63];
ybs[6120]=['',4.2938328,1.0759597,5.67];
ybs[6121]=['',4.3252299,-0.8078938,5.35];
ybs[6122]=['η Dra',4.2947672,1.0727762,2.74];
ybs[6123]=['',4.5738472,-1.5272869,6.57];
ybs[6124]=['α Sco',4.3228696,-0.4621234,0.96];
ybs[6125]=['',4.3490175,-1.239728,5.5];
ybs[6126]=['',4.3182223,0.0107995,5.39];
ybs[6127]=['',4.3196037,-0.1426804,6.48];
ybs[6128]=['',4.4106598,-1.4534351,6.57];
ybs[6129]=['',4.4917228,-1.5048043,6.04];
ybs[6130]=['',4.3240474,-0.2547562,5.68];
ybs[6131]=['22 Sco',4.3263068,-0.4391316,4.79];
ybs[6132]=['',4.3336294,-0.7306222,5.33];
ybs[6133]=['',4.3318695,-0.6064887,4.23];
ybs[6134]=['',4.3269283,-0.1319519,6.5];
ybs[6135]=['',4.3314877,-0.4639538,6.1];
ybs[6136]=['30 Her',4.3168451,0.7301651,5.04];
ybs[6137]=['φ Oph',4.3300434,-0.2907324,4.28];
ybs[6138]=['β Her',4.324698,0.3742727,2.77];
ybs[6139]=['λ Oph',4.3284169,0.0338382,3.82];
ybs[6140]=['',4.3165155,0.8964271,6.29];
ybs[6141]=['θ TrA',4.3538611,-1.1438501,5.52];
ybs[6142]=['',4.3262202,0.3566381,5.25];
ybs[6143]=['ω Oph',4.3345846,-0.3754352,4.45];
ybs[6144]=['',4.3290468,0.3865954,5.76];
ybs[6145]=['μ Nor',4.3441991,-0.7694935,4.94];
ybs[6146]=['34 Her',4.3227172,0.8537316,6.45];
ybs[6147]=['',4.3276997,0.6140049,6.25];
ybs[6148]=['28 Her',4.3356339,0.0955885,5.63];
ybs[6149]=['29 Her',4.3354708,0.1997312,4.84];
ybs[6150]=['',4.3488573,-0.7904188,6.46];
ybs[6151]=['15 Dra',4.3107552,1.1994125,5];
ybs[6152]=['',4.3303244,0.7950595,5.65];
ybs[6153]=['β Aps',4.3905015,-1.3536101,4.24];
ybs[6154]=['',4.3541287,-0.7487677,5.47];
ybs[6155]=['τ Sco',4.3512173,-0.4932079,2.82];
ybs[6156]=['',4.3536964,-0.6160649,4.16];
ybs[6157]=['',4.3667792,-1.0651971,6.18];
ybs[6158]=['σ Her',4.340625,0.7399026,4.2];
ybs[6159]=['',4.3476297,0.2969554,6.41];
ybs[6160]=['',4.3316241,1.0607897,5.94];
ybs[6161]=['12 Oph',4.3523215,-0.0413149,5.75];
ybs[6162]=['η1 TrA',4.3791454,-1.192685,5.91];
ybs[6163]=['',4.2960305,1.3773397,5.56];
ybs[6164]=['',4.3631652,-0.7581701,5.83];
ybs[6165]=['ζ Oph',4.3560977,-0.1851665,2.56];
ybs[6166]=['',4.3532639,0.2697537,6.3];
ybs[6167]=['',4.3752784,-1.0556829,6.18];
ybs[6168]=['',4.3656378,-0.6502841,5.91];
ybs[6169]=['',4.3597461,-0.1148371,6.09];
ybs[6170]=['',4.3247354,1.2665291,6.3];
ybs[6171]=['',4.3580559,0.2381531,6.31];
ybs[6172]=['',4.3875891,-1.1775952,6.03];
ybs[6173]=['',4.3494288,0.8128117,5.79];
ybs[6174]=['16 Dra',4.3489354,0.9225391,5.53];
ybs[6175]=['17 Dra',4.349093,0.9229613,5.08];
ybs[6176]=['17 Dra',4.3491221,0.9229565,6.53];
ybs[6177]=['',4.3762194,-0.851772,5.65];
ybs[6178]=['',4.3777395,-0.8672783,5.65];
ybs[6179]=['',4.3669426,-0.167469,6.35];
ybs[6180]=['',4.3713727,-0.3569016,6.26];
ybs[6181]=['',4.3186648,1.3509016,6.34];
ybs[6182]=['',4.3770913,-0.579207,5.87];
ybs[6183]=['',4.3760244,-0.4277432,6.09];
ybs[6184]=['36 Her',4.3705398,0.072725,6.93];
ybs[6185]=['37 Her',4.3708012,0.0729436,5.77];
ybs[6186]=['',4.3756272,-0.3103559,4.96];
ybs[6187]=['',4.3835199,-0.8047644,6.23];
ybs[6188]=['',4.3508348,1.1000877,6.16];
ybs[6189]=['',4.3565321,0.9769262,5.29];
ybs[6190]=['42 Her',4.3604223,0.8532378,4.9];
ybs[6191]=['',4.373366,-0.0181625,6.24];
ybs[6192]=['',4.3771149,-0.34844,5.57];
ybs[6193]=['',4.3714374,0.2156307,6.08];
ybs[6194]=['',4.4019285,-1.1719333,5.13];
ybs[6195]=['14 Oph',4.3755453,0.0199191,5.74];
ybs[6196]=['',4.3862569,-0.7183358,6.2];
ybs[6197]=['',4.3911225,-0.9283532,5.96];
ybs[6198]=['',4.3716274,0.4331626,6.06];
ybs[6199]=['',4.3868677,-0.7182328,6.12];
ybs[6200]=['',4.3862371,-0.6666351,6.05];
ybs[6201]=['',4.3852764,-0.5610345,6.46];
ybs[6202]=['ζ Her',4.3725353,0.5508775,2.81];
ybs[6203]=['39 Her',4.3741645,0.4690924,5.92];
ybs[6204]=['',4.3903947,-0.7134552,5.71];
ybs[6205]=['',4.3990946,-1.0217324,5.74];
ybs[6206]=['',4.3878712,-0.4798714,6.58];
ybs[6207]=['α TrA',4.4111409,-1.2053918,1.92];
ybs[6208]=['',4.3910388,-0.4982542,6.02];
ybs[6209]=['',4.4033182,-1.0188927,5.58];
ybs[6210]=['η Her',4.3791781,0.6786344,3.53];
ybs[6211]=['',4.399418,-0.6879118,5.48];
ybs[6212]=['',4.383647,0.5934128,5.99];
ybs[6213]=['18 Dra',4.3680032,1.1265876,4.83];
ybs[6214]=['16 Oph',4.3920731,0.0171447,6.03];
ybs[6215]=['25 Sco',4.3989908,-0.4462081,6.71];
ybs[6216]=['',4.3782158,0.9712916,6.16];
ybs[6217]=['',4.3910263,0.2741428,5.56];
ybs[6218]=['43 Her',4.3932782,0.149133,5.15];
ybs[6219]=['η Ara',4.4141532,-1.0310886,3.76];
ybs[6220]=['',4.3889397,0.7536159,6.05];
ybs[6221]=['',4.4266364,-1.1818714,6.32];
ybs[6222]=['19 Oph',4.3993162,0.0353833,6.1];
ybs[6223]=['',4.4244151,-1.1416212,6.13];
ybs[6224]=['45 Her',4.4018669,0.0909287,5.24];
ybs[6225]=['',4.4055318,-0.2608553,6.03];
ybs[6226]=['',4.4167596,-0.8740755,6.47];
ybs[6227]=['',4.3882213,0.9903647,4.85];
ybs[6228]=['',4.3489285,1.3766471,6.32];
ybs[6229]=['',4.4032049,0.2365599,6.35];
ybs[6230]=['',4.4099823,-0.274077,6.1];
ybs[6231]=['ε Sco',4.4138384,-0.5991522,2.29];
ybs[6232]=['',4.3983124,0.7365593,5.87];
ybs[6233]=['20 Oph',4.4114175,-0.1888243,4.65];
ybs[6234]=['',4.4176474,-0.6553638,6.11];
ybs[6235]=['',4.4203311,-0.720217,5.22];
ybs[6236]=['',4.4094483,0.2308273,5.91];
ybs[6237]=['μ1 Sco',4.4214896,-0.6646597,3.08];
ybs[6238]=['',4.4134793,-0.0469392,6.32];
ybs[6239]=['',4.4236734,-0.7310993,6.49];
ybs[6240]=['47 Her',4.4129222,0.1258767,5.49];
ybs[6241]=['',4.4324543,-1.0112957,5.94];
ybs[6242]=['μ2 Sco',4.4235178,-0.664132,3.57];
ybs[6243]=['',4.4393698,-1.1048372,6.02];
ybs[6244]=['52 Her',4.4063744,0.8019285,4.82];
ybs[6245]=['21 Oph',4.4178968,0.020614,5.51];
ybs[6246]=['',4.4084614,0.7573781,6.13];
ybs[6247]=['',4.4298185,-0.751968,5.96];
ybs[6248]=['50 Her',4.4134584,0.5196056,5.72];
ybs[6249]=['',4.4136258,0.5675493,6.13];
ybs[6250]=['',4.4311429,-0.7302456,5.45];
ybs[6251]=['',4.4309365,-0.7335331,6.32];
ybs[6252]=['ζ1 Sco',4.4310254,-0.739947,4.73];
ybs[6253]=['',4.4318733,-0.7310101,6.45];
ybs[6254]=['',4.412601,0.7306148,6.29];
ybs[6255]=['',4.4324388,-0.7304806,6.59];
ybs[6256]=['',4.433012,-0.7419793,5.88];
ybs[6257]=['',4.3727702,1.3521855,5.98];
ybs[6258]=['49 Her',4.4203133,0.2607427,6.52];
ybs[6259]=['',4.4274333,-0.3569115,5.88];
ybs[6260]=['51 Her',4.4185055,0.4297262,5.04];
ybs[6261]=['ζ2 Sco',4.4335938,-0.7399273,3.62];
ybs[6262]=['',4.4352201,-0.7188006,5.77];
ybs[6263]=['',4.4330224,-0.5344298,6.35];
ybs[6264]=['',4.4410473,-0.8850128,6.33];
ybs[6265]=['',4.442634,-0.9130902,5.94];
ybs[6266]=['',4.4588325,-1.2094948,5.79];
ybs[6267]=['',4.4300456,-0.0287255,6.25];
ybs[6268]=['',4.432568,-0.2064002,6.57];
ybs[6269]=['53 Her',4.4234809,0.5526997,5.32];
ybs[6270]=['23 Oph',4.4320308,-0.1079887,5.25];
ybs[6271]=['ι Oph',4.4288953,0.1768288,4.38];
ybs[6272]=['',4.4391159,-0.5853814,6.37];
ybs[6273]=['',4.4423011,-0.7130704,6.15];
ybs[6274]=['',4.4386829,-0.2938921,6.37];
ybs[6275]=['ζ Ara',4.4523669,-0.9777597,3.13];
ybs[6276]=['',4.4239571,0.8269849,6];
ybs[6277]=['',4.4324641,0.3652156,5.41];
ybs[6278]=['27 Sco',4.4444358,-0.5810461,5.48];
ybs[6279]=['',4.4504391,-0.8844024,5.55];
ybs[6280]=['',4.4342588,0.2371311,6.34];
ybs[6281]=['24 Oph',4.4422959,-0.4046068,5.58];
ybs[6282]=['56 Her',4.4327674,0.4485025,6.08];
ybs[6283]=['54 Her',4.4345269,0.3211452,5.35];
ybs[6284]=['',4.4433134,-0.3415982,6.27];
ybs[6285]=['ε1 Ara',4.4562833,-0.9283635,4.06];
ybs[6286]=['',4.4445964,-0.1919043,6.19];
ybs[6287]=['',4.4586926,-0.9534285,5.65];
ybs[6288]=['',4.452037,-0.6571565,6.09];
ybs[6289]=['κ Oph',4.444892,0.1630679,3.2];
ybs[6290]=['',4.4596658,-0.8495934,6];
ybs[6291]=['',4.4441378,0.2417664,6.37];
ybs[6292]=['',4.4502322,-0.2600724,6.59];
ybs[6293]=['',4.4598008,-0.7938102,6.65];
ybs[6294]=['',4.4665732,-1.0295337,6.11];
ybs[6295]=['57 Her',4.443609,0.4419307,6.28];
ybs[6296]=['',4.4360135,0.8727709,6.56];
ybs[6297]=['',4.444472,0.4249785,6.32];
ybs[6298]=['',4.4561694,-0.4384722,5.86];
ybs[6299]=['26 Oph',4.4570306,-0.4366767,5.75];
ybs[6300]=['',4.4595467,-0.6276984,5.97];
ybs[6301]=['',4.4656259,-0.8929191,6.45];
ybs[6302]=['',4.4441135,0.7414264,6.34];
ybs[6303]=['ε2 Ara',4.4718423,-0.9296653,5.29];
ybs[6304]=['19 Dra',4.4337377,1.1362403,4.89];
ybs[6305]=['',4.4648659,-0.5615298,5.03];
ybs[6306]=['',4.4573031,0.1143739,6.59];
ybs[6307]=['30 Oph',4.4601726,-0.0742228,4.82];
ybs[6308]=['20 Dra',4.4354634,1.1345759,6.41];
ybs[6309]=['',4.4778615,-1.0077619,5.73];
ybs[6310]=['29 Oph',4.4641773,-0.3301339,6.26];
ybs[6311]=['ε UMi',4.380067,1.4311472,4.23];
ybs[6312]=['',4.4737084,-0.8235983,6.06];
ybs[6313]=['ε Her',4.4554539,0.5392328,3.92];
ybs[6314]=['',4.4587807,0.3944786,5.65];
ybs[6315]=['',4.4616213,0.2603943,6.31];
ybs[6316]=['',4.4737761,-0.6663822,5.91];
ybs[6317]=['',4.4594209,0.4741398,6.55];
ybs[6318]=['',4.4637588,0.1469714,6.33];
ybs[6319]=['',4.4495307,0.9888583,6.03];
ybs[6320]=['',4.4796677,-0.7946427,6.28];
ybs[6321]=['59 Her',4.4610667,0.5853549,5.25];
ybs[6322]=['',4.4645191,0.4446395,5.75];
ybs[6323]=['',4.4778201,-0.5960466,4.87];
ybs[6324]=['',4.4325511,1.2757453,6.3];
ybs[6325]=['',4.4641168,0.5559765,6.36];
ybs[6326]=['',4.4685468,0.2454421,4.98];
ybs[6327]=['',4.4827196,-0.77026,6.19];
ybs[6328]=['',4.4687201,0.2527583,6.52];
ybs[6329]=['',4.4769047,-0.3581935,6.3];
ybs[6330]=['',4.4708555,0.2369528,5.93];
ybs[6331]=['',4.4722166,0.2362962,6.08];
ybs[6332]=['',4.471595,0.3431627,6.35];
ybs[6333]=['',4.4845998,-0.6502206,5.98];
ybs[6334]=['',4.4459257,1.2069807,6.4];
ybs[6335]=['61 Her',4.4692435,0.6175874,6.69];
ybs[6336]=['',4.4850918,-0.6192157,6.13];
ybs[6337]=['',4.4573834,1.0579994,6.13];
ybs[6338]=['',4.4784189,0.0117717,6.01];
ybs[6339]=['',4.4832399,-0.3768556,6.3];
ybs[6340]=['',4.4709703,0.606702,6.04];
ybs[6341]=['',4.4751323,0.3415748,6.17];
ybs[6342]=['',4.4795953,-0.0160543,5.64];
ybs[6343]=['',4.486462,-0.463214,6.29];
ybs[6344]=['60 Her',4.4784064,0.2218807,4.91];
ybs[6345]=['',4.5032791,-1.076884,6.39];
ybs[6346]=['',4.5150518,-1.2347367,6.22];
ybs[6347]=['',4.4819348,0.169397,6.37];
ybs[6348]=['',4.4821557,0.1819784,6.37];
ybs[6349]=['',4.4609787,1.1269721,6.1];
ybs[6350]=['',4.4854839,-0.0293843,6.38];
ybs[6351]=['',4.4756312,0.7641746,6.43];
ybs[6352]=['',4.4741551,0.8512979,6.09];
ybs[6353]=['',4.4820696,0.3849604,5.56];
ybs[6354]=['',4.4920349,-0.3078001,5.99];
ybs[6355]=['',4.4949554,-0.5310999,5.97];
ybs[6356]=['',4.4913321,-0.0193029,6.06];
ybs[6357]=['',4.5182672,-1.173216,5.89];
ybs[6358]=['μ Dra',4.4758122,0.9501941,5.83];
ybs[6359]=['μ Dra',4.4758049,0.9501941,5.8];
ybs[6360]=['',4.5041191,-0.778114,5.08];
ybs[6361]=['',4.4944064,-0.0682241,6.36];
ybs[6362]=['',4.5352288,-1.3012279,6.25];
ybs[6363]=['',4.5085654,-0.8534453,5.84];
ybs[6364]=['',4.4985378,-0.1841155,5.56];
ybs[6365]=['',4.4875982,0.7066707,6.34];
ybs[6366]=['',4.4889744,0.6267226,5.39];
ybs[6367]=['η Oph',4.501257,-0.2748915,2.43];
ybs[6368]=['',4.4549904,1.3136555,6.21];
ybs[6369]=['η Sco',4.5103506,-0.7550919,3.33];
ybs[6370]=['',4.5106229,-0.6899516,5.67];
ybs[6371]=['',4.510604,-0.6780009,6.3];
ybs[6372]=['',4.489048,0.8868991,6.46];
ybs[6373]=['',4.5205875,-0.9932951,6.09];
ybs[6374]=['',4.5019106,0.2171532,6.57];
ybs[6375]=['',4.5097119,-0.4412045,6.54];
ybs[6376]=['',4.5106529,-0.4849619,6.14];
ybs[6377]=['',4.4953311,0.7112438,5.08];
ybs[6378]=['',4.5133277,-0.56658,6.01];
ybs[6379]=['',4.5063983,0.137357,6.33];
ybs[6380]=['63 Her',4.502708,0.4225902,6.19];
ybs[6381]=['',4.5202047,-0.6944698,6.6];
ybs[6382]=['37 Oph',4.5093974,0.1843222,5.33];
ybs[6383]=['',4.5117014,0.0057212,6.65];
ybs[6384]=['',4.4986062,0.9142623,6.29];
ybs[6385]=['ζ Dra',4.4892245,1.1464756,3.17];
ybs[6386]=['',4.5236115,-0.5859273,5.53];
ybs[6387]=['',4.5250961,-0.673986,5.96];
ybs[6388]=['',4.5039053,0.8678083,6.04];
ybs[6389]=['',4.5492707,-1.2228757,6.53];
ybs[6390]=['36 Oph',4.5234023,-0.4647046,5.11];
ybs[6391]=['36 Oph',4.5233877,-0.4646803,5.07];
ybs[6392]=['',4.5257878,-0.5276674,6.21];
ybs[6393]=['',4.5228578,-0.254936,5.99];
ybs[6394]=['',4.5282526,-0.6243345,6.12];
ybs[6395]=['α1 Her',4.5187878,0.2507514,3.48];
ybs[6396]=['α2 Her',4.5188096,0.2507466,5.39];
ybs[6397]=['',4.5427244,-1.0422259,5.91];
ybs[6398]=['',4.5311604,-0.5704558,5.55];
ybs[6399]=['δ Her',4.5200373,0.4331219,3.14];
ybs[6400]=['ι Aps',4.5575274,-1.2242161,5.41];
ybs[6401]=['',4.5261656,0.0377632,6.17];
ybs[6402]=['',4.5285484,-0.1093829,6.09];
ybs[6403]=['',4.5274732,0.0207392,5.88];
ybs[6404]=['41 Oph',4.5278961,-0.0081598,4.73];
ybs[6405]=['',4.5407122,-0.814279,5.48];
ybs[6406]=['ζ Aps',4.5563952,-1.183154,4.78];
ybs[6407]=['π Her',4.5195179,0.6420374,3.16];
ybs[6408]=['',4.5229711,0.4139922,5.96];
ybs[6409]=['',4.5394241,-0.7705755,5.76];
ybs[6410]=['',4.5061392,1.0969374,5.56];
ybs[6411]=['',4.5367339,-0.5685342,6.36];
ybs[6412]=['',4.5429632,-0.8741297,6.27];
ybs[6413]=['ο Oph',4.5349163,-0.4242618,5.2];
ybs[6414]=['ο Oph',4.5349016,-0.4242134,6.8];
ybs[6415]=['',4.5395374,-0.6110516,5.91];
ybs[6416]=['',4.5420858,-0.772199,6.65];
ybs[6417]=['',4.5359358,-0.2850694,6.43];
ybs[6418]=['',4.6055526,-1.4114999,5.88];
ybs[6419]=['',4.5313095,0.4026308,6.45];
ybs[6420]=['68 Her',4.5296541,0.5773208,4.82];
ybs[6421]=['',4.5336607,0.3018814,6];
ybs[6422]=['',4.536234,0.1892497,5.03];
ybs[6423]=['',4.5375583,0.10584,6.51];
ybs[6424]=['',4.5428364,-0.3102655,6.02];
ybs[6425]=['69 Her',4.5309411,0.6504821,4.65];
ybs[6426]=['',4.5263266,0.8668846,7.48];
ybs[6427]=['',4.5586947,-1.0127985,5.88];
ybs[6428]=['',4.5428415,-0.1036327,6.32];
ybs[6429]=['',4.5642278,-1.0975042,5.7];
ybs[6430]=['',4.545874,-0.3377726,6.52];
ybs[6431]=['',4.5593775,-0.9868788,5.8];
ybs[6432]=['',4.5363335,0.5026875,5.65];
ybs[6433]=['',4.5339793,0.6770125,5.94];
ybs[6434]=['ξ Oph',4.5478378,-0.3688355,4.39];
ybs[6435]=['ν Ser',4.5467485,-0.2245714,4.33];
ybs[6436]=['',4.5651549,-1.0592695,5.77];
ybs[6437]=['',4.5236987,1.058508,6.32];
ybs[6438]=['',4.546886,-0.1870319,6.46];
ybs[6439]=['',4.5558568,-0.660154,6.41];
ybs[6440]=['ι Ara',4.559175,-0.8288048,5.25];
ybs[6441]=['',4.5433678,0.3148023,5];
ybs[6442]=['θ Oph',4.5523961,-0.4366614,3.27];
ybs[6443]=['',4.5556482,-0.6270803,6.47];
ybs[6444]=['',4.5423938,0.4453561,5.38];
ybs[6445]=['',4.5569507,-0.6499512,5.93];
ybs[6446]=['70 Her',4.5456676,0.4272453,5.12];
ybs[6447]=['72 Her',4.5442343,0.5663164,5.39];
ybs[6448]=['43 Oph',4.5584291,-0.4915155,5.35];
ybs[6449]=['',4.5630687,-0.7710988,5.12];
ybs[6450]=['β Ara',4.5688008,-0.9694883,2.85];
ybs[6451]=['γ Ara',4.5693058,-0.9842791,3.34];
ybs[6452]=['',4.5488374,0.2916636,6.35];
ybs[6453]=['74 Her',4.5420498,0.8066979,5.59];
ybs[6454]=['',4.5552078,-0.0420164,6.29];
ybs[6455]=['',4.5481661,0.5015775,6.35];
ybs[6456]=['',4.5428068,0.8406899,6.43];
ybs[6457]=['κ Ara',4.5713646,-0.8840244,5.23];
ybs[6458]=['',4.5484713,0.6973416,5.51];
ybs[6459]=['',4.5661078,-0.6058772,6.16];
ybs[6460]=['',4.5822015,-1.1004729,6.24];
ybs[6461]=['',4.5639833,-0.3745424,5.85];
ybs[6462]=['',4.5635046,-0.3222561,6.21];
ybs[6463]=['',4.5658623,-0.4234417,6.19];
ybs[6464]=['',4.5755818,-0.9069766,6.19];
ybs[6465]=['',4.5596217,0.1541876,5.77];
ybs[6466]=['',4.5747201,-0.8004159,5.29];
ybs[6467]=['',4.5766313,-0.8839553,5.92];
ybs[6468]=['',4.5475536,0.9320197,5.67];
ybs[6469]=['73 Her',4.5597233,0.4004106,5.74];
ybs[6470]=['',4.5618047,0.2841904,5.71];
ybs[6471]=['',4.5619991,0.2720606,6.35];
ybs[6472]=['',4.5800747,-0.9130421,5.75];
ybs[6473]=['ρ Her',4.5571512,0.6480052,5.47];
ybs[6474]=['ρ Her',4.557173,0.6479908,4.52];
ybs[6475]=['44 Oph',4.571387,-0.4222377,4.17];
ybs[6476]=['',4.5833619,-0.9631702,5.94];
ybs[6477]=['',4.5586313,0.6730733,6.49];
ybs[6478]=['',4.5687513,-0.0291311,6.44];
ybs[6479]=['',4.573865,-0.4530911,6.44];
ybs[6480]=['',4.5605357,0.6446137,6.28];
ybs[6481]=['45 Oph',4.5759448,-0.5215669,4.29];
ybs[6482]=['',4.5717847,-0.089077,4.54];
ybs[6483]=['',4.5771163,-0.5190773,6];
ybs[6484]=['',4.5678086,0.2949608,5.98];
ybs[6485]=['',4.5738194,-0.2186783,6.21];
ybs[6486]=['',4.5699504,0.1322663,6.06];
ybs[6487]=['σ Oph',4.5709418,0.0719622,4.34];
ybs[6488]=['',4.5678623,0.4688203,6.41];
ybs[6489]=['δ Ara',4.5947836,-1.0593876,3.62];
ybs[6490]=['',4.5831988,-0.6421787,6.02];
ybs[6491]=['',4.5716591,0.3501795,5.54];
ybs[6492]=['',4.5854521,-0.6725139,6.39];
ybs[6493]=['',4.5780368,-0.1435475,6.37];
ybs[6494]=['',4.5955179,-0.993703,5.95];
ybs[6495]=['',4.5707728,0.6052581,5.94];
ybs[6496]=['',4.5811733,0.005491,5.44];
ybs[6497]=['υ Sco',4.5912074,-0.6511943,2.69];
ybs[6498]=['77 Her',4.5697418,0.8419963,5.85];
ybs[6499]=['α Ara',4.596795,-0.8707509,2.95];
ybs[6500]=['',4.5639002,1.0477308,5.65];
ybs[6501]=['',4.5855962,-0.103588,6.37];
ybs[6502]=['',4.5963826,-0.8037356,6.03];
ybs[6503]=['',4.5658092,1.0233629,6.51];
ybs[6504]=['',4.5880669,-0.0188082,5.31];
ybs[6505]=['',4.5954914,-0.5884745,6.44];
ybs[6506]=['',4.55956,1.1744002,6.43];
ybs[6507]=['51 Oph',4.5934007,-0.4184833,4.81];
ybs[6508]=['',4.5949136,-0.458744,6.05];
ybs[6509]=['',4.587489,0.2078657,6.39];
ybs[6510]=['',4.5982203,-0.5985384,6.17];
ybs[6511]=['',4.6017412,-0.7188527,5.84];
ybs[6512]=['',4.5921113,0.047295,5.59];
ybs[6513]=['',4.6142301,-1.0447252,6.28];
ybs[6514]=['λ Her',4.5885035,0.4554531,4.41];
ybs[6515]=['λ Sco',4.6036124,-0.6478184,1.63];
ybs[6516]=['',4.5890733,0.5435548,5.61];
ybs[6517]=['',4.5291829,1.3982704,5.72];
ybs[6518]=['',4.612362,-0.9313981,6.1];
ybs[6519]=['',4.5875574,0.6783594,6.43];
ybs[6520]=['',4.5956769,0.2079748,6.42];
ybs[6521]=['78 Her',4.5931464,0.4955521,5.62];
ybs[6522]=['',4.6017708,-0.1005001,5.62];
ybs[6523]=['',4.6081691,-0.5688811,5.7];
ybs[6524]=['β Dra',4.5855079,0.9125648,2.79];
ybs[6525]=['σ Ara',4.6131845,-0.8118893,4.59];
ybs[6526]=['',4.5936911,0.597888,6.56];
ybs[6527]=['',4.6128325,-0.6536656,6.48];
ybs[6528]=['',4.5861831,1.009869,6.4];
ybs[6529]=['',4.6003206,0.3358543,5.64];
ybs[6530]=['',4.6016452,0.2845589,5.69];
ybs[6531]=['',4.6019501,0.2588012,6.48];
ybs[6532]=['',4.6075381,-0.196433,5.55];
ybs[6533]=['52 Oph',4.6103009,-0.3849572,6.57];
ybs[6534]=['',4.6165246,-0.6745195,4.29];
ybs[6535]=['',4.621303,-0.8739094,5.93];
ybs[6536]=['53 Oph',4.6060691,0.1670926,5.81];
ybs[6537]=['π Ara',4.6245283,-0.9514007,5.25];
ybs[6538]=['',4.5980992,0.7195955,5.74];
ybs[6539]=['',4.6051137,0.2878192,6.4];
ybs[6540]=['',4.7489834,-1.4820481,6.45];
ybs[6541]=['θ Sco',4.6201763,-0.7506522,1.87];
ybs[6542]=['ν1 Dra',4.5928134,0.9628943,4.88];
ybs[6543]=['ν2 Dra',4.5932072,0.9627012,4.87];
ybs[6544]=['α Oph',4.6073679,0.2189899,2.08];
ybs[6545]=['',4.6204185,-0.6645679,6.26];
ybs[6546]=['',4.6237542,-0.7485988,6.1];
ybs[6547]=['',4.6116509,0.366237,6.1];
ybs[6548]=['',4.5984284,1.0043529,6.17];
ybs[6549]=['ξ Ser',4.6199735,-0.2689549,3.54];
ybs[6550]=['',4.6200529,-0.2719654,5.94];
ybs[6551]=['',4.6096148,0.6508192,6.1];
ybs[6552]=['',4.6119385,0.4917029,6.38];
ybs[6553]=['',4.6553201,-1.2606223,6.49];
ybs[6554]=['27 Dra',4.5897062,1.1889248,5.05];
ybs[6555]=['μ Oph',4.6208246,-0.1418976,4.62];
ybs[6556]=['',4.6222966,-0.1908949,5.75];
ybs[6557]=['λ Ara',4.6340779,-0.8626351,4.77];
ybs[6558]=['',4.6139081,0.5370954,6.02];
ybs[6559]=['79 Her',4.6181632,0.4240889,5.77];
ybs[6560]=['',4.637707,-0.8191056,5.79];
ybs[6561]=['26 Dra',4.6041933,1.0796959,5.23];
ybs[6562]=['82 Her',4.6128558,0.8477678,5.37];
ybs[6563]=['',4.6261119,0.0352113,6.26];
ybs[6564]=['',4.6414695,-0.8817312,6.24];
ybs[6565]=['',4.6249145,0.2324508,6.12];
ybs[6566]=['',4.6308671,-0.0377434,6.19];
ybs[6567]=['',4.6234809,0.571222,6.37];
ybs[6568]=['κ Sco',4.6424772,-0.6813547,2.41];
ybs[6569]=['ο Ser',4.6365809,-0.2248797,4.26];
ybs[6570]=['η Pav',4.6593646,-1.1297653,3.62];
ybs[6571]=['',4.6439335,-0.6449759,5.54];
ybs[6572]=['',4.6284884,0.5444076,6.03];
ybs[6573]=['μ Ara',4.6506693,-0.9048135,5.15];
ybs[6574]=['',4.6547209,-1.0044834,6.01];
ybs[6575]=['',4.6448785,-0.5769981,6.4];
ybs[6576]=['ι Her',4.6254412,0.8027788,3.8];
ybs[6577]=['',4.6345291,0.2647401,6.34];
ybs[6578]=['',4.6364127,0.1100155,5.95];
ybs[6579]=['',4.6316618,0.5458977,6.28];
ybs[6580]=['',4.6337445,0.4276703,6.36];
ybs[6581]=['',4.6454179,-0.4868165,6.36];
ybs[6582]=['',4.6379675,0.2782543,5.52];
ybs[6583]=['58 Oph',4.6457256,-0.3785906,4.87];
ybs[6584]=['ω Dra',4.6112844,1.1998442,4.8];
ybs[6585]=['',4.6523643,-0.7458921,5.87];
ybs[6586]=['',4.6097868,1.2140269,6.42];
ybs[6587]=['',4.6306917,0.7585358,6.59];
ybs[6588]=['',4.6470558,-0.2359116,6.39];
ybs[6589]=['',4.6466975,-0.123702,6.3];
ybs[6590]=['83 Her',4.6397706,0.42857,5.52];
ybs[6591]=['β Oph',4.6449196,0.0795674,2.77];
ybs[6592]=['',4.6440897,0.2493478,6.24];
ybs[6593]=['',4.6293199,1.000078,6.77];
ybs[6594]=['',4.6109084,1.2643831,5.86];
ybs[6595]=['',4.6332318,0.9042283,5.99];
ybs[6596]=['84 Her',4.6436494,0.4244523,5.71];
ybs[6597]=['61 Oph',4.6497627,0.0448841,6.17];
ybs[6598]=['',4.6498645,0.0448747,6.56];
ybs[6599]=['',4.6481065,0.2513681,6.19];
ybs[6600]=['',4.64141,0.7692676,6.34];
ybs[6601]=['',4.6626414,-0.6652847,6.43];
ybs[6602]=['',4.6706335,-0.9670369,6.11];
ybs[6603]=['ι1 Sco',4.6647852,-0.7004536,3.03];
ybs[6604]=['3 Sgr',4.6640189,-0.4858469,4.54];
ybs[6605]=['',4.6646534,-0.3924217,6.18];
ybs[6606]=['',4.6444856,0.9388724,5.75];
ybs[6607]=['',4.6534,0.5497343,6.23];
ybs[6608]=['',4.6636986,-0.2571216,5.94];
ybs[6609]=['',4.6679211,-0.4709016,6.35];
ybs[6610]=['',4.6784476,-0.9357887,5.92];
ybs[6611]=['μ Her',4.6570146,0.4836955,3.42];
ybs[6612]=['',4.6841497,-1.0501354,5.78];
ybs[6613]=['',4.6539711,0.6784831,6.52];
ybs[6614]=['',4.6542935,0.6861826,6.68];
ybs[6615]=['',4.6603902,0.3087619,5.72];
ybs[6616]=['',4.6712534,-0.5534199,4.83];
ybs[6617]=['γ Oph',4.6642732,0.0471445,3.75];
ybs[6618]=['',4.6745212,-0.6466139,3.21];
ybs[6619]=['ι2 Sco',4.6761357,-0.6997948,4.81];
ybs[6620]=['',4.6815111,-0.9273755,6.09];
ybs[6621]=['',4.6661675,0.0662938,6.22];
ybs[6622]=['',4.6925726,-1.1430528,6.49];
ybs[6623]=['',4.7156379,-1.3295553,6.07];
ybs[6624]=['ψ1 Dra',4.6319237,1.2590695,4.58];
ybs[6625]=['ψ1 Dra',4.6320441,1.2592104,5.79];
ybs[6626]=['',4.6658579,0.3588353,5.69];
ybs[6627]=['',4.6705036,0.0341354,6.47];
ybs[6628]=['',4.6832957,-0.7959478,6.11];
ybs[6629]=['',4.6587778,0.8308754,6.43];
ybs[6630]=['',4.6675917,0.3359702,6.12];
ybs[6631]=['',4.6821131,-0.7116845,5.96];
ybs[6632]=['87 Her',4.6674144,0.4471039,5.12];
ybs[6633]=['',4.6800689,-0.5333979,6.66];
ybs[6634]=['',4.7549003,-1.4221366,6.35];
ybs[6635]=['',4.6847308,-0.6074243,5.9];
ybs[6636]=['',4.6851538,-0.6007476,5.84];
ybs[6637]=['',4.6880031,-0.733038,6.2];
ybs[6638]=['',4.6762857,0.2084287,6.17];
ybs[6639]=['',4.6872831,-0.5954634,6.06];
ybs[6640]=['',4.687819,-0.611248,6.45];
ybs[6641]=['',4.6879976,-0.6218165,6.03];
ybs[6642]=['',4.6740568,0.5116854,5.5];
ybs[6643]=['',4.6762268,0.3894148,5.98];
ybs[6644]=['30 Dra',4.6669565,0.8862004,5.02];
ybs[6645]=['',4.689527,-0.6062168,6.17];
ybs[6646]=['',4.6898049,-0.6090912,5.6];
ybs[6647]=['',4.6822992,-0.0216518,6.35];
ybs[6648]=['',4.6914135,-0.6071777,6.38];
ybs[6649]=['',4.6853369,-0.107288,6.21];
ybs[6650]=['',4.6920954,-0.6065944,5.96];
ybs[6651]=['',4.6923323,-0.6079757,6.42];
ybs[6652]=['88 Her',4.6714536,0.8445494,6.68];
ybs[6653]=['',4.6816092,0.2674176,6.46];
ybs[6654]=['',4.6873178,-0.1902938,6.18];
ybs[6655]=['',4.6848174,0.0227139,5.95];
ybs[6656]=['',4.6944289,-0.601596,5.96];
ybs[6657]=['',4.6771891,0.6993204,6.46];
ybs[6658]=['',4.6874633,0.1064324,5.77];
ybs[6659]=['',4.6974993,-0.6366611,6.06];
ybs[6660]=['',4.6959131,-0.4344045,6.2];
ybs[6661]=['',4.6808815,0.6977524,6.04];
ybs[6662]=['',4.6801521,0.8140097,6.38];
ybs[6663]=['',4.7052324,-0.7739399,4.86];
ybs[6664]=['',4.6916389,0.1942169,6.38];
ybs[6665]=['90 Her',4.6862174,0.6982144,5.16];
ybs[6666]=['',4.7055737,-0.7034858,6.43];
ybs[6667]=['',4.7000914,-0.3281922,6.52];
ybs[6668]=['',4.7039046,-0.4898554,5.8];
ybs[6669]=['',4.7017171,-0.2760081,5.89];
ybs[6670]=['',4.7094383,-0.7280966,4.88];
ybs[6671]=['',4.7100104,-0.6830803,6.29];
ybs[6672]=['',4.7010575,0.0116699,5.82];
ybs[6673]=['89 Her',4.6961921,0.4546205,5.46];
ybs[6674]=['',4.7033617,-0.0712674,5.47];
ybs[6675]=['',4.698211,0.39204,5.58];
ybs[6676]=['ξ Dra',4.6857591,0.9925599,3.75];
ybs[6677]=['',4.7044174,0.001137,5.97];
ybs[6678]=['',4.7035787,0.1132099,6.29];
ybs[6679]=['',4.71408,-0.6433024,5.74];
ybs[6680]=['',4.7124604,-0.5019482,6.01];
ybs[6681]=['',4.7144304,-0.5280174,5.16];
ybs[6682]=['',4.7144522,-0.5280223,7.04];
ybs[6683]=['θ Her',4.6992608,0.650114,3.86];
ybs[6684]=['',4.7056556,0.1927384,6.36];
ybs[6685]=['',4.7042171,0.4187852,6.3];
ybs[6686]=['ν Oph',4.7133099,-0.1705852,3.34];
ybs[6687]=['',4.6940108,0.9768448,6.1];
ybs[6688]=['4 Sgr',4.7172174,-0.4156655,4.76];
ybs[6689]=['35 Dra',4.6621543,1.3431534,5.04];
ybs[6690]=['',4.7011489,0.7914951,6.02];
ybs[6691]=['ξ Her',4.7062753,0.5104534,3.7];
ybs[6692]=['',4.7179788,-0.3549796,6.21];
ybs[6693]=['γ Dra',4.6997659,0.8986219,2.23];
ybs[6694]=['',4.7156789,-0.0841474,5.87];
ybs[6695]=['ν Her',4.7094512,0.5268953,4.41];
ybs[6696]=['',4.7266117,-0.6348889,6.3];
ybs[6697]=['',4.7183013,0.0109933,6.37];
ybs[6698]=['ζ Ser',4.7194407,-0.0643979,4.62];
ybs[6699]=['',4.7100127,0.6333328,6];
ybs[6700]=['66 Oph',4.7181865,0.0762539,4.64];
ybs[6701]=['93 Her',4.7168155,0.2923621,4.67];
ybs[6702]=['67 Oph',4.7199044,0.051178,3.97];
ybs[6703]=['6 Sgr',4.7238852,-0.2994269,6.28];
ybs[6704]=['',4.7263937,-0.3975724,5.77];
ybs[6705]=['',4.6642464,1.3666142,6.24];
ybs[6706]=['',4.7101497,0.7937004,6.48];
ybs[6707]=['',4.7208146,0.1094159,6.34];
ybs[6708]=['',4.7184715,0.3404495,6.5];
ybs[6709]=['χ Oct',5.0041338,-1.5271279,5.28];
ybs[6710]=['',4.7207934,0.2634415,6.26];
ybs[6711]=['68 Oph',4.7247998,0.0228024,4.45];
ybs[6712]=['7 Sgr',4.7305821,-0.4237727,5.34];
ybs[6713]=['ψ2 Dra',4.6897655,1.2566787,5.45];
ybs[6714]=['',4.7184764,0.5797008,5.99];
ybs[6715]=['',4.7312852,-0.3964762,6.74];
ybs[6716]=['',4.7147947,0.7941513,5.67];
ybs[6717]=['95 Her',4.7229079,0.3769265,5.18];
ybs[6718]=['95 Her',4.7229442,0.3769314,4.96];
ybs[6719]=['',4.7746322,-1.3244383,5.86];
ybs[6720]=['',4.7294867,-0.0934949,6.76];
ybs[6721]=['τ Oph',4.7309509,-0.1427394,5.94];
ybs[6722]=['τ Oph',4.7309436,-0.1427442,5.24];
ybs[6723]=['',4.6850643,1.3119246,6.36];
ybs[6724]=['9 Sgr',4.7350434,-0.4251306,5.97];
ybs[6725]=['',4.7227982,0.5814116,6.15];
ybs[6726]=['96 Her',4.7268086,0.363641,5.28];
ybs[6727]=['',4.7398274,-0.6265467,6];
ybs[6728]=['',4.7556388,-1.1265293,6.41];
ybs[6729]=['97 Her',4.7272293,0.4001097,6.21];
ybs[6730]=['γ1 Sgr',4.7402936,-0.5162164,4.69];
ybs[6731]=['θ Ara',4.7486295,-0.8741964,3.66];
ybs[6732]=['',4.7306117,0.3423462,6.5];
ybs[6733]=['π Pav',4.7588513,-1.1111344,4.35];
ybs[6734]=['γ2 Sgr',4.743767,-0.5309427,2.99];
ybs[6735]=['',4.7372929,0.0335429,6.14];
ybs[6736]=['',4.7466186,-0.6285978,5.95];
ybs[6737]=['',4.7489746,-0.7578351,5.77];
ybs[6738]=['',4.7489746,-0.7578351,5.77];
ybs[6739]=['',4.7791821,-1.2856861,5.85];
ybs[6740]=['70 Oph',4.7409078,0.0436781,4.03];
ybs[6741]=['',4.7285833,0.8458954,6.21];
ybs[6742]=['',4.7366463,0.4179221,6.34];
ybs[6743]=['',4.7442206,-0.1452182,5.85];
ybs[6744]=['',4.7446559,-0.0828652,5.77];
ybs[6745]=['',4.7439314,-0.007735,6.34];
ybs[6746]=['',4.7417185,0.2095641,7.04];
ybs[6747]=['',4.7564364,-0.798704,6.15];
ybs[6748]=['',4.7641338,-1.030343,6.38];
ybs[6749]=['ι Pav',4.7666524,-1.0820387,5.49];
ybs[6750]=['',4.7493894,-0.3741951,6.28];
ybs[6751]=['',4.740381,0.3778648,6.15];
ybs[6752]=['',4.7360341,0.6996469,6.52];
ybs[6753]=['98 Her',4.7426621,0.3878519,5.06];
ybs[6754]=['',4.7535986,-0.4965923,4.57];
ybs[6755]=['',4.7371948,0.7321562,6.34];
ybs[6756]=['',4.741304,0.5625861,5.71];
ybs[6757]=['',4.7519044,-0.2993199,5.52];
ybs[6758]=['71 Oph',4.7487576,0.1525064,4.64];
ybs[6759]=['72 Oph',4.7489156,0.166993,3.73];
ybs[6760]=['',4.7596499,-0.6399638,6.58];
ybs[6761]=['',4.7570321,-0.4444918,6.61];
ybs[6762]=['',4.7860036,-1.2347023,6.73];
ybs[6763]=['99 Her',4.7466163,0.5334744,5.04];
ybs[6764]=['',4.7507753,0.2282096,6.63];
ybs[6765]=['',4.7621636,-0.5709694,6.43];
ybs[6766]=['',4.7677985,-0.829151,6.07];
ybs[6767]=['ο Her',4.748965,0.5020729,3.83];
ybs[6768]=['',4.7624851,-0.5362171,5.53];
ybs[6769]=['100 Her',4.7503172,0.4556305,5.86];
ybs[6770]=['100 Her',4.7503173,0.4555626,5.9];
ybs[6771]=['ε Tel',4.7683539,-0.8019468,4.53];
ybs[6772]=['',4.7540228,0.2493979,6.37];
ybs[6773]=['',4.7601408,-0.2431077,6.39];
ybs[6774]=['',4.7674374,-0.7217456,5.86];
ybs[6775]=['102 Her',4.7546189,0.3633646,4.36];
ybs[6776]=['',4.7662411,-0.5898104,6.16];
ybs[6777]=['δ UMi',4.561344,1.5082532,4.36];
ybs[6778]=['',4.7447051,0.88709,6.29];
ybs[6779]=['',4.7478559,0.7586203,5];
ybs[6780]=['',4.7457353,0.8676801,6.32];
ybs[6781]=['',4.7507243,0.6354008,5.48];
ybs[6782]=['101 Her',4.755189,0.3499413,5.1];
ybs[6783]=['73 Oph',4.7587792,0.0697888,5.73];
ybs[6784]=['',4.7836765,-1.1114511,6.47];
ybs[6785]=['',4.7602801,0.0545445,5.69];
ybs[6786]=['',4.7670226,-0.3462038,6.36];
ybs[6787]=['',4.7559807,0.5318794,6.38];
ybs[6788]=['',4.7636397,0.0581198,5.51];
ybs[6789]=['11 Sgr',4.7692654,-0.4135496,4.98];
ybs[6790]=['',4.7705775,-0.5043091,6.51];
ybs[6791]=['',4.7608456,0.2876689,6.09];
ybs[6792]=['',4.7766706,-0.7213239,5.47];
ybs[6793]=['',4.7897195,-1.100375,5.6];
ybs[6794]=['',4.7575363,0.671301,6.4];
ybs[6795]=['',4.7592151,0.636553,5.58];
ybs[6796]=['',4.8010135,-1.1906491,6.33];
ybs[6797]=['40 Dra',4.7055957,1.3962711,6.04];
ybs[6798]=['41 Dra',4.7060134,1.3963301,5.68];
ybs[6799]=['24 UMi',4.5508657,1.5154538,5.79];
ybs[6800]=['μ Sgr',4.7780536,-0.3674157,3.86];
ybs[6801]=['',4.7747936,-0.0698917,6.59];
ybs[6802]=['',4.7670952,0.58387,5.88];
ybs[6803]=['104 Her',4.7678614,0.5482376,4.97];
ybs[6804]=['14 Sgr',4.7802699,-0.3788285,5.44];
ybs[6805]=['',4.7602253,0.947579,5.95];
ybs[6806]=['',4.7885057,-0.7714002,5.46];
ybs[6807]=['',4.7950088,-0.9776274,5.33];
ybs[6808]=['',4.7742826,0.382008,6.12];
ybs[6809]=['',4.7940168,-0.891148,6.06];
ybs[6810]=['15 Sgr',4.7843735,-0.3616333,5.38];
ybs[6811]=['16 Sgr',4.7843594,-0.3556943,5.95];
ybs[6812]=['',4.7708322,0.7182686,6.36];
ybs[6813]=['',4.7855904,-0.3255557,6.07];
ybs[6814]=['',4.7725939,0.6768496,6.04];
ybs[6815]=['',4.7620351,1.0544462,6.49];
ybs[6816]=['',4.8072671,-1.114849,6.18];
ybs[6817]=['φ Oct',4.8279493,-1.30954,5.47];
ybs[6818]=['',4.7869964,-0.0629867,6.36];
ybs[6819]=['',4.7803178,0.5099002,6.56];
ybs[6820]=['η Sgr',4.7956633,-0.641445,3.11];
ybs[6821]=['',4.7954088,-0.5951165,6.16];
ybs[6822]=['',4.7873308,0.0416516,6.01];
ybs[6823]=['',4.7942566,-0.4999159,6.19];
ybs[6824]=['',4.7942177,-0.4935745,6.4];
ybs[6825]=['',4.8571647,-1.4000421,5.95];
ybs[6826]=['',4.7928686,-0.3030693,5.75];
ybs[6827]=['',4.8005388,-0.7378938,6.3];
ybs[6828]=['',4.790981,-0.052327,6];
ybs[6829]=['',4.7941413,-0.3220811,6.54];
ybs[6830]=['',4.7970322,-0.4718102,4.65];
ybs[6831]=['',4.7934914,-0.170156,6.31];
ybs[6832]=['',4.7916858,0.0177156,6.63];
ybs[6833]=['',4.7835812,0.7359667,5.59];
ybs[6834]=['',4.7997667,-0.4467104,6.51];
ybs[6835]=['',4.782918,0.7891981,6.29];
ybs[6836]=['',4.7995869,-0.3247944,6.84];
ybs[6837]=['',4.7781118,0.9877882,6.37];
ybs[6838]=['36 Dra',4.7734803,1.1240707,5.03];
ybs[6839]=['',4.795437,0.2406218,6.3];
ybs[6840]=['',4.7956149,0.3166218,5.99];
ybs[6841]=['',4.7900611,0.7146386,6.11];
ybs[6842]=['',4.7954001,0.4067726,6.63];
ybs[6843]=['ξ Pav',4.8223667,-1.0730505,4.36];
ybs[6844]=['',4.8100799,-0.6540832,6.45];
ybs[6845]=['',4.8005257,0.1268853,5.39];
ybs[6846]=['',4.8056932,-0.2761257,5.39];
ybs[6847]=['δ Sgr',4.8099972,-0.5204004,2.7];
ybs[6848]=['105 Her',4.7999319,0.4268437,5.27];
ybs[6849]=['',4.8120759,-0.4346469,6.25];
ybs[6850]=['',4.8162051,-0.6744812,5.1];
ybs[6851]=['',4.8112168,-0.3289687,5.75];
ybs[6852]=['',4.8143299,-0.4959907,6.16];
ybs[6853]=['37 Dra',4.7785202,1.2001543,5.95];
ybs[6854]=['74 Oph',4.8081294,0.0591384,4.86];
ybs[6855]=['',4.8027034,0.517956,5.99];
ybs[6856]=['106 Her',4.8049261,0.3834874,4.95];
ybs[6857]=['η Ser',4.8102846,-0.0503961,3.26];
ybs[6858]=['',4.8185985,-0.6397876,5.34];
ybs[6859]=['',4.832646,-1.0996893,6.14];
ybs[6860]=['κ Lyr',4.8023441,0.6296275,4.33];
ybs[6861]=['',4.8107013,0.0950735,6.13];
ybs[6862]=['',4.8211774,-0.632258,5.55];
ybs[6863]=['',4.825253,-0.7696414,5.25];
ybs[6864]=['108 Her',4.8074216,0.5213304,5.63];
ybs[6865]=['107 Her',4.8077525,0.5040717,5.12];
ybs[6866]=['',4.8180878,-0.1781334,6.33];
ybs[6867]=['ε Sgr',4.8240867,-0.5999002,1.85];
ybs[6868]=['',4.8015774,0.8963716,6.3];
ybs[6869]=['',4.8188827,-0.2094799,5.73];
ybs[6870]=['',4.8129279,0.4066104,5.41];
ybs[6871]=['',4.8153142,0.2101538,5.89];
ybs[6872]=['ζ Sct',4.8207523,-0.15571,4.68];
ybs[6873]=['',4.816085,0.3113459,5.25];
ybs[6874]=['',4.8068878,0.8680694,6.4];
ybs[6875]=['',4.8171412,0.2914756,6.22];
ybs[6876]=['18 Sgr',4.8276242,-0.5365712,5.6];
ybs[6877]=['',4.8293653,-0.6279358,6.15];
ybs[6878]=['',4.8223012,-0.062317,6.38];
ybs[6879]=['',4.8088018,0.8575334,5.05];
ybs[6880]=['',4.8252359,-0.1232569,6.31];
ybs[6881]=['',4.8316472,-0.5922148,6.3];
ybs[6882]=['',4.8368753,-0.8395471,5.46];
ybs[6883]=['109 Her',4.8197707,0.3801731,3.84];
ybs[6884]=['21 Sgr',4.8285848,-0.3582833,4.81];
ybs[6885]=['α Tel',4.8370349,-0.8020463,3.51];
ybs[6886]=['',4.826119,-0.0273345,6.15];
ybs[6887]=['',4.8680286,-1.290631,5.89];
ybs[6888]=['',4.8267351,0.0889787,6.74];
ybs[6889]=['',4.8200608,0.6763472,6.36];
ybs[6890]=['',4.8288087,0.1404219,5.65];
ybs[6891]=['μ Lyr',4.8212073,0.6897548,5.12];
ybs[6892]=['',4.8250948,0.4783689,6.27];
ybs[6893]=['ζ Tel',4.8453874,-0.856178,4.13];
ybs[6894]=['',4.829754,0.2614578,6.37];
ybs[6895]=['',4.8397941,-0.520135,5.92];
ybs[6896]=['',4.8511082,-1.0036862,5.76];
ybs[6897]=['',4.8392313,-0.4646054,6.31];
ybs[6898]=['',4.8430182,-0.6803354,5.64];
ybs[6899]=['',4.8182415,0.9304936,6.32];
ybs[6900]=['',4.9157014,-1.4274128,6.27];
ybs[6901]=['λ Sgr',4.8402241,-0.4434315,2.81];
ybs[6902]=['',4.8408657,-0.46674,6.27];
ybs[6903]=['',4.8466605,-0.7649816,6.36];
ybs[6904]=['ν Pav',4.8580236,-1.0866679,4.64];
ybs[6905]=['',4.829366,0.5208524,5.83];
ybs[6906]=['59 Ser',4.8359084,0.0036752,5.21];
ybs[6907]=['',4.8397808,-0.3104089,6.2];
ybs[6908]=['φ Dra',4.8016101,1.2452667,4.22];
ybs[6909]=['',4.846623,-0.6778069,6.63];
ybs[6910]=['',4.850028,-0.82387,5.7];
ybs[6911]=['39 Dra',4.8180855,1.0264823,4.98];
ybs[6912]=['',4.8325854,0.4618765,6.53];
ybs[6913]=['',4.8385286,0.0656836,6.07];
ybs[6914]=['',4.8445806,-0.4636683,6.5];
ybs[6915]=['χ Dra',4.8023787,1.2696161,3.57];
ybs[6916]=['',4.8390576,0.1083678,5.73];
ybs[6917]=['',4.8463101,-0.4405342,6.59];
ybs[6918]=['γ Sct',4.8451367,-0.2539507,4.7];
ybs[6919]=['',4.854535,-0.7312414,6.04];
ybs[6920]=['',4.8476752,-0.2542219,5.96];
ybs[6921]=['',4.8496638,-0.3266006,5.66];
ybs[6922]=['δ1 Tel',4.8578962,-0.8010722,4.96];
ybs[6923]=['60 Ser',4.8467866,-0.0343747,5.39];
ybs[6924]=['',4.8541519,-0.5754808,5.34];
ybs[6925]=['',4.8585094,-0.7590568,5.72];
ybs[6926]=['δ2 Tel',4.8590913,-0.798316,5.07];
ybs[6927]=['',4.8667293,-1.0243556,6.44];
ybs[6928]=['',4.8493439,-0.0996256,6.28];
ybs[6929]=['',4.8483202,0.0712308,6.69];
ybs[6930]=['',4.8600957,-0.6926629,5.16];
ybs[6931]=['',4.8454129,0.4168153,5.9];
ybs[6932]=['',4.8550606,-0.3208977,5.14];
ybs[6933]=['42 Dra',4.8260085,1.1445375,4.82];
ybs[6934]=['',4.8547232,-0.1881318,5.72];
ybs[6935]=['',4.857053,-0.3334987,6.68];
ybs[6936]=['',4.862979,-0.6959439,6.22];
ybs[6937]=['',4.8345545,1.0395823,6.43];
ybs[6938]=['',4.8503348,0.3635782,6.5];
ybs[6939]=['θ CrA',4.8652628,-0.7381812,4.64];
ybs[6940]=['κ1 CrA',4.864537,-0.6754863,6.32];
ybs[6941]=['κ2 CrA',4.8645228,-0.6755882,5.65];
ybs[6942]=['',4.8705467,-0.9228173,6.22];
ybs[6943]=['',4.8521499,0.2957467,5.77];
ybs[6944]=['',4.8588759,-0.2552895,6.37];
ybs[6945]=['61 Ser',4.8566412,-0.0172113,5.94];
ybs[6946]=['',4.8571999,0.0641708,6.43];
ybs[6947]=['',4.8605206,-0.2591501,5.5];
ybs[6948]=['',4.8667317,-0.5759347,5.28];
ybs[6949]=['24 Sgr',4.865998,-0.4191326,5.49];
ybs[6950]=['',4.8645704,-0.2589332,5.76];
ybs[6951]=['',4.8630657,-0.1028648,6.36];
ybs[6952]=['',4.9612266,-1.4536532,7.16];
ybs[6953]=['25 Sgr',4.8688638,-0.4224428,6.51];
ybs[6954]=['',4.859279,0.4124951,5.84];
ybs[6955]=['',4.8625806,0.1446175,6.42];
ybs[6956]=['',4.8592316,0.5335728,5.48];
ybs[6957]=['',4.8722475,-0.3634144,6.48];
ybs[6958]=['',4.8704879,-0.1912652,5.14];
ybs[6959]=['',4.8616298,0.5394827,6.59];
ybs[6960]=['',4.8754308,-0.5180153,6.37];
ybs[6961]=['α Sct',4.8711137,-0.143563,3.85];
ybs[6962]=['',4.8549935,0.9098827,6.56];
ybs[6963]=['',4.8662043,0.3575218,6.57];
ybs[6964]=['',4.8686068,0.1904159,6.4];
ybs[6965]=['',4.8701506,0.3180321,5.78];
ybs[6966]=['45 Dra',4.8561416,0.9959304,4.77];
ybs[6967]=['',4.8490066,1.1423598,6.59];
ybs[6968]=['',4.871222,0.4123211,5.61];
ybs[6969]=['',4.8731521,0.2966094,6.21];
ybs[6970]=['ζ Pav',4.9110706,-1.2462535,4.01];
ybs[6971]=['',4.8626311,0.9140534,5.36];
ybs[6972]=['',4.8694783,0.601725,6.1];
ybs[6973]=['',4.8759669,0.1595533,5.39];
ybs[6974]=['',4.8906834,-0.8358192,5.86];
ybs[6975]=['',4.8768786,0.1167848,5.45];
ybs[6976]=['',4.8834027,-0.373112,5.94];
ybs[6977]=['',4.8838511,-0.2440775,6.47];
ybs[6978]=['',4.8861308,-0.4098843,5.81];
ybs[6979]=['',4.8918418,-0.7533738,5.37];
ybs[6980]=['',4.8791396,0.1996879,6.42];
ybs[6981]=['',4.8812688,-0.0050546,5.75];
ybs[6982]=['',4.9353476,-1.358592,6.39];
ybs[6983]=['',4.8786954,0.2830558,6.29];
ybs[6984]=['',4.90642,-1.1278408,6.37];
ybs[6985]=['',4.8756183,0.5844782,5.42];
ybs[6986]=['',4.8876784,-0.3670671,5.86];
ybs[6987]=['',4.8848414,-0.0553856,6.49];
ybs[6988]=['',4.884432,-0.0190786,6.66];
ybs[6989]=['α Lyr',4.8767139,0.6772402,0.03];
ybs[6990]=['',4.8842101,0.154533,6.4];
ybs[6991]=['',4.8756545,0.7547018,6.2];
ybs[6992]=['',4.9118298,-1.1262248,5.78];
ybs[6993]=['',4.900606,-0.8390325,6.49];
ybs[6994]=['',4.8376096,1.3537152,5.64];
ybs[6995]=['',4.8920473,-0.1356027,5.84];
ybs[6996]=['',4.8898602,0.0922411,6.38];
ybs[6997]=['',4.8817674,0.6926872,6.04];
ybs[6998]=['',4.8908531,0.1287933,6.28];
ybs[6999]=['26 Sgr',4.9007523,-0.4155848,6.23];
ybs[7000]=['',4.919863,-1.1317983,4.79];
ybs[7001]=['',4.8707174,1.143321,6.06];
ybs[7002]=['',4.8997171,-0.2538091,6.42];
ybs[7003]=['',4.9180882,-1.0658909,6.04];
ybs[7004]=['',4.8906465,0.5387911,6.36];
ybs[7005]=['',4.8879887,0.7148124,6.25];
ybs[7006]=['',4.877111,1.0916379,5.74];
ybs[7007]=['',4.8909886,0.6700023,6.45];
ybs[7008]=['δ Sct',4.9019759,-0.1576075,4.72];
ybs[7009]=['λ CrA',4.9098501,-0.6684699,5.13];
ybs[7010]=['',4.9184059,-0.9923532,6.22];
ybs[7011]=['',4.9051882,-0.3361728,6.35];
ybs[7012]=['',4.9033361,-0.1230617,6.15];
ybs[7013]=['',4.8052005,1.4518892,6.17];
ybs[7014]=['',4.9112711,-0.6404544,6.32];
ybs[7015]=['',4.94082,-1.2735368,6.06];
ybs[7016]=['',4.8885537,0.9113579,6];
ybs[7017]=['',4.9120626,-0.6216614,4.87];
ybs[7018]=['',4.8978317,0.5522159,6.41];
ybs[7019]=['',4.9150415,-0.6922444,5.43];
ybs[7020]=['ε Sct',4.907394,-0.1440312,4.9];
ybs[7021]=['',4.8996209,0.6068291,6.47];
ybs[7022]=['',4.9088026,-0.1186047,6.31];
ybs[7023]=['',4.9137453,-0.4361143,5.83];
ybs[7024]=['θ Pav',4.9337737,-1.1353722,5.73];
ybs[7025]=['',4.9246506,-0.8738803,6.54];
ybs[7026]=['',4.9156919,-0.3661273,6.36];
ybs[7027]=['φ Sgr',4.9174537,-0.4706595,3.17];
ybs[7028]=['4 Aql',4.912738,0.0363647,5.02];
ybs[7029]=['',4.9043547,0.6863146,6.45];
ybs[7030]=['',4.8918311,1.0955613,6.09];
ybs[7031]=['',4.9059239,0.6384326,6.01];
ybs[7032]=['',4.9072943,0.5576265,5.7];
ybs[7033]=['',4.9187255,-0.3417737,6.42];
ybs[7034]=['28 Sgr',4.9202502,-0.3903927,5.37];
ybs[7035]=['',4.9112116,0.4121269,6.31];
ybs[7036]=['',4.9154143,0.0964048,5.83];
ybs[7037]=['46 Dra',4.9002204,0.9697344,5.04];
ybs[7038]=['μ CrA',4.9272546,-0.704781,5.24];
ybs[7039]=['ε1 Lyr',4.9089664,0.6927769,5.06];
ybs[7040]=['ε1 Lyr',4.9089591,0.6927963,6.02];
ybs[7041]=['ε2 Lyr',4.9091517,0.6917834,5.14];
ybs[7042]=['ε2 Lyr',4.9091517,0.6917786,5.37];
ybs[7043]=['',4.9214235,-0.1762864,5.71];
ybs[7044]=['ζ1 Lyr',4.9109789,0.6567399,4.36];
ybs[7045]=['ζ2 Lyr',4.9111105,0.6565559,5.73];
ybs[7046]=['',4.9153158,0.3841276,6.51];
ybs[7047]=['5 Aql',4.9200164,-0.0163585,5.9];
ybs[7048]=['',4.9041177,0.9406388,6.11];
ybs[7049]=['110 Her',4.915666,0.3590197,4.19];
ybs[7050]=['η1 CrA',4.9322566,-0.7619111,5.49];
ybs[7051]=['β Sct',4.9232061,-0.0824323,4.22];
ybs[7052]=['',4.9172092,0.4657646,4.83];
ybs[7053]=['',4.9350891,-0.7990858,5.81];
ybs[7054]=['',4.9245859,-0.0991361,5.2];
ybs[7055]=['',4.9202277,0.3269053,6.17];
ybs[7056]=['η2 CrA',4.9354747,-0.7576091,5.61];
ybs[7057]=['111 Her',4.9216959,0.3177549,4.36];
ybs[7058]=['',4.9336648,-0.6060255,6.62];
ybs[7059]=['',4.9102962,0.9585364,6.23];
ybs[7060]=['',4.9306172,-0.3242088,6.47];
ybs[7061]=['',4.917048,0.723715,6.07];
ybs[7062]=['λ Pav',4.948887,-1.0848959,4.22];
ybs[7063]=['',4.9067995,1.0658921,5.99];
ybs[7064]=['',4.926681,0.0744657,6.21];
ybs[7065]=['',4.9342888,-0.3336406,6.75];
ybs[7066]=['29 Sgr',4.9346708,-0.3542784,5.24];
ybs[7067]=['',4.9269401,0.4108404,6.15];
ybs[7068]=['',4.9200587,0.8087772,6.52];
ybs[7069]=['',4.9251897,0.5547008,6.06];
ybs[7070]=['',4.8996668,1.2359572,6.44];
ybs[7071]=['',4.9341922,-0.1027429,5.99];
ybs[7072]=['',4.9183299,0.9252411,5.88];
ybs[7073]=['',4.9336677,0.0150416,6.25];
ybs[7074]=['',4.9298035,0.3377944,5.88];
ybs[7075]=['κ Tel',4.9495798,-0.9089641,5.17];
ybs[7076]=['30 Sgr',4.9398647,-0.3863383,6.61];
ybs[7077]=['',4.9371002,-0.1375516,6.8];
ybs[7078]=['',4.922824,0.8569542,6.4];
ybs[7079]=['',4.9310932,0.4375913,6.59];
ybs[7080]=['',4.9482091,-0.81276,5.54];
ybs[7081]=['',4.9519253,-0.9058808,6.31];
ybs[7082]=['',4.9399678,-0.1701253,5.83];
ybs[7083]=['',4.950919,-0.8435548,6.19];
ybs[7084]=['',4.925476,0.8515927,6.12];
ybs[7085]=['',4.9505706,-0.8125904,6.19];
ybs[7086]=['',4.932943,0.5524865,6.64];
ybs[7087]=['',4.9382838,0.1920375,6.55];
ybs[7088]=['ν1 Lyr',4.9330224,0.5731446,5.91];
ybs[7089]=['8 Aql',4.9414512,-0.0574368,6.1];
ybs[7090]=['ν2 Lyr',4.9335447,0.5685739,5.25];
ybs[7091]=['',4.9471822,-0.464655,6.29];
ybs[7092]=['',4.9479244,-0.5122865,6.13];
ybs[7093]=['',4.9483303,-0.5359302,6.63];
ybs[7094]=['β Lyr',4.9343693,0.5827467,3.45];
ybs[7095]=['κ Pav',4.9705038,-1.1729244,4.44];
ybs[7096]=['',4.9575713,-0.8700462,6.6];
ybs[7097]=['',4.9437216,0.2442195,6.14];
ybs[7098]=['',4.9489337,-0.1666504,6.34];
ybs[7099]=['',4.9694546,-1.0955639,6.48];
ybs[7100]=['',4.9412184,0.5028389,6.18];
ybs[7101]=['112 Her',4.9444838,0.374418,5.48];
ybs[7102]=['33 Sgr',4.9536171,-0.3723041,5.69];
ybs[7103]=['',4.9408718,0.6381889,6.09];
ybs[7104]=['ν1 Sgr',4.9544086,-0.3964802,4.83];
ybs[7105]=['',4.9097542,1.293449,5.27];
ybs[7106]=['',4.9428126,0.7227495,6.28];
ybs[7107]=['',4.9565167,-0.2718252,5.1];
ybs[7108]=['ν2 Sgr',4.9585425,-0.3951871,4.99];
ybs[7109]=['σ Sgr',4.9593372,-0.4584586,2.02];
ybs[7110]=['',4.9646403,-0.7449251,5.36];
ybs[7111]=['',4.9395712,0.9250566,5.51];
ybs[7112]=['50 Dra',4.9116395,1.3169864,5.35];
ybs[7113]=['ο Dra',4.9371621,1.0369862,4.66];
ybs[7114]=['',4.9600292,-0.2853202,5.79];
ybs[7115]=['ω Pav',4.97644,-1.0501611,5.14];
ybs[7116]=['',4.9624532,-0.4039446,5.93];
ybs[7117]=['',4.9660318,-0.6512463,5.38];
ybs[7118]=['',4.9839162,-1.1627693,6.01];
ybs[7119]=['δ1 Lyr',4.9500941,0.6457659,5.58];
ybs[7120]=['',4.9527006,0.4876051,5.62];
ybs[7121]=['113 Her',4.9552347,0.3957281,4.59];
ybs[7122]=['λ Tel',4.9749567,-0.9234184,4.87];
ybs[7123]=['',4.9589441,0.1159636,5.57];
ybs[7124]=['',4.970115,-0.6945224,6.31];
ybs[7125]=['',4.9469571,0.8855106,4.92];
ybs[7126]=['',4.9520768,0.7200146,7.3];
ybs[7127]=['δ2 Lyr',4.9534882,0.6445026,4.3];
ybs[7128]=['',4.9552669,0.593363,6.02];
ybs[7129]=['θ1 Ser',4.9623537,0.073879,4.62];
ybs[7130]=['θ2 Ser',4.9624556,0.0738501,4.98];
ybs[7131]=['',4.9632571,-0.0309022,6.22];
ybs[7132]=['',4.963317,0.043643,6.15];
ybs[7133]=['ξ1 Sgr',4.9681592,-0.3599991,5.08];
ybs[7134]=['',4.9548103,0.7266039,5.44];
ybs[7135]=['',4.9611682,0.3145822,6.63];
ybs[7136]=['',4.9613312,0.3165072,5.69];
ybs[7137]=['η Sct',4.9663822,-0.101514,4.83];
ybs[7138]=['ξ2 Sgr',4.9698715,-0.3678545,3.51];
ybs[7139]=['',4.9730336,-0.54115,6.12];
ybs[7140]=['ε CrA',4.9749527,-0.6471123,4.87];
ybs[7141]=['',4.9486313,1.0038238,6.22];
ybs[7142]=['',4.9539178,0.8532602,5.77];
ybs[7143]=['',4.9726967,-0.4336478,6.62];
ybs[7144]=['',4.9771051,-0.689471,6.49];
ybs[7145]=['13 Lyr',4.9566861,0.7675067,4.04];
ybs[7146]=['64 Ser',4.9670236,0.0447704,5.57];
ybs[7147]=['',4.9728961,-0.3926806,6.14];
ybs[7148]=['',4.9047498,1.3956661,6.39];
ybs[7149]=['',4.9994031,-1.1994246,5.88];
ybs[7150]=['',4.9647128,0.5747556,5.22];
ybs[7151]=['',4.9717791,0.1094446,6.21];
ybs[7152]=['',4.9772592,-0.3235129,6.37];
ybs[7153]=['',4.9707064,0.3035331,5.38];
ybs[7154]=['',4.9768233,-0.2235693,5.53];
ybs[7155]=['10 Aql',4.9731799,0.2432514,5.89];
ybs[7156]=['',4.9817366,-0.4347738,6.36];
ybs[7157]=['',4.9851124,-0.6462773,6.69];
ybs[7158]=['',4.9851924,-0.6462965,6.4];
ybs[7159]=['',4.9728221,0.3460072,6.5];
ybs[7160]=['11 Aql',4.9745579,0.2382946,5.23];
ybs[7161]=['',4.9755442,0.1775298,6.75];
ybs[7162]=['',4.9688125,0.6683962,5.89];
ybs[7163]=['48 Dra',4.9615964,1.0095757,5.66];
ybs[7164]=['ε Aql',4.9768091,0.2635337,4.02];
ybs[7165]=['',4.9901597,-0.7309013,6.23];
ybs[7166]=['γ Lyr',4.9730892,0.5710733,3.24];
ybs[7167]=['',4.9719216,0.7105184,6.22];
ybs[7168]=['υ Dra',4.9485509,1.2448613,4.82];
ybs[7169]=['',4.9769473,0.458352,5.27];
ybs[7170]=['',4.9869497,-0.3955515,6.24];
ybs[7171]=['',4.9780092,0.3987364,6.29];
ybs[7172]=['',4.9647033,1.0167428,6.46];
ybs[7173]=['',4.973879,0.6850162,6.41];
ybs[7174]=['',4.9863453,-0.2661701,6.32];
ybs[7175]=['',4.9590056,1.1394775,5.63];
ybs[7176]=['ζ CrA',4.9944016,-0.7341265,4.75];
ybs[7177]=['',4.9987383,-0.8898599,5.93];
ybs[7178]=['',4.9632986,1.089545,6.45];
ybs[7179]=['λ Lyr',4.9777937,0.5615902,4.93];
ybs[7180]=['12 Aql',4.9865286,-0.099602,4.02];
ybs[7181]=['ζ Sgr',4.9915461,-0.5209395,2.6];
ybs[7182]=['',4.9906677,-0.4330929,5.65];
ybs[7183]=['',4.9721256,0.8873263,6.3];
ybs[7184]=['',4.9949559,-0.6670705,5.74];
ybs[7185]=['',4.9830534,0.3375726,6.39];
ybs[7186]=['',4.9428321,1.3232219,6.22];
ybs[7187]=['',4.9842364,0.364172,6.69];
ybs[7188]=['',4.9786439,0.7106193,6.65];
ybs[7189]=['',4.9836296,0.4594272,5.69];
ybs[7190]=['',4.9930635,-0.3353203,6.05];
ybs[7191]=['',4.981661,0.5905124,6.01];
ybs[7192]=['',4.9932906,-0.3328425,6.37];
ybs[7193]=['',4.9849573,0.4373417,6.72];
ybs[7194]=['',4.9861349,0.389139,6.4];
ybs[7195]=['',4.9889973,0.1467228,6.3];
ybs[7196]=['14 Aql',4.9918158,-0.0639866,5.42];
ybs[7197]=['',4.977542,0.882523,5.38];
ybs[7198]=['',4.9994815,-0.5412859,5.5];
ybs[7199]=['',4.9855332,0.587364,6.39];
ybs[7200]=['ρ Tel',5.0092135,-0.9129163,5.16];
ybs[7201]=['',4.9943637,0.0323221,5.83];
ybs[7202]=['16 Lyr',4.9831239,0.8197214,5.01];
ybs[7203]=['',4.9908296,0.3437212,6.09];
ybs[7204]=['ο Sgr',5.0002302,-0.3788763,3.77];
ybs[7205]=['49 Dra',4.9792006,0.9719701,5.48];
ybs[7206]=['',4.9971105,0.0587113,6.73];
ybs[7207]=['',4.9984087,-0.0986376,6.9];
ybs[7208]=['',5.0271257,-1.1935998,5.33];
ybs[7209]=['',4.9943949,0.37177,6.52];
ybs[7210]=['',5.0114614,-0.842371,5.97];
ybs[7211]=['',4.968677,1.214077,6.52];
ybs[7212]=['15 Aql',5.0007792,-0.0697718,5.42];
ybs[7213]=['γ CrA',5.0085114,-0.6462739,4.93];
ybs[7214]=['γ CrA',5.0085114,-0.6462739,4.99];
ybs[7215]=['σ Oct',5.6124218,-1.5507379,5.47];
ybs[7216]=['',4.9856308,0.91269,6.31];
ybs[7217]=['',5.0043699,-0.2727274,5.97];
ybs[7218]=['',5.002218,-0.0258108,6.53];
ybs[7219]=['',5.0105444,-0.6593065,6.16];
ybs[7220]=['',5.0206139,-0.9718762,6.49];
ybs[7221]=['τ Sgr',5.0103252,-0.4823349,3.32];
ybs[7222]=['ζ Aql',5.0021069,0.2425533,2.99];
ybs[7223]=['λ Aql',5.0064303,-0.0846152,3.44];
ybs[7224]=['',4.9994193,0.554633,5.56];
ybs[7225]=['',4.9994968,0.5369859,6.06];
ybs[7226]=['',5.009561,-0.2826412,6.03];
ybs[7227]=['',5.0128766,-0.4991965,6.04];
ybs[7228]=['',5.0108351,-0.3264271,6.29];
ybs[7229]=['δ CrA',5.0171254,-0.7061799,4.59];
ybs[7230]=['',5.0065008,0.1442417,6.09];
ybs[7231]=['',5.003085,0.5228267,6.31];
ybs[7232]=['',5.0101785,0.0118074,6.56];
ybs[7233]=['',5.01588,-0.4297263,6.3];
ybs[7234]=['',4.9612545,1.3453093,6.54];
ybs[7235]=['18 Aql',5.0090461,0.1938385,5.09];
ybs[7236]=['',5.0158152,-0.3360553,5.54];
ybs[7237]=['',5.0070671,0.4238597,5.77];
ybs[7238]=['51 Dra',4.9977626,0.9325336,5.38];
ybs[7239]=['',4.9991492,0.871915,6.43];
ybs[7240]=['',5.0068253,0.500266,5.55];
ybs[7241]=['α CrA',5.0218742,-0.6609275,4.11];
ybs[7242]=['',5.0228131,-0.694499,6.46];
ybs[7243]=['',5.0223682,-0.6305625,6.56];
ybs[7244]=['',5.024254,-0.7305277,5.88];
ybs[7245]=['',5.0046463,0.7234076,6.49];
ybs[7246]=['β CrA',5.0243836,-0.6859922,4.11];
ybs[7247]=['',5.0131012,0.2947608,6.07];
ybs[7248]=['17 Lyr',5.0101274,0.5678704,5.23];
ybs[7249]=['ι Lyr',5.0093976,0.6306768,5.28];
ybs[7250]=['',5.0133697,0.3793324,6.23];
ybs[7251]=['π Sgr',5.0223616,-0.3662995,2.89];
ybs[7252]=['',5.0224806,-0.3450062,6.13];
ybs[7253]=['19 Aql',5.0180463,0.1066238,5.22];
ybs[7254]=['',5.0162216,0.294733,6.48];
ybs[7255]=['',5.0287326,-0.6801222,6.36];
ybs[7256]=['',5.022033,-0.0068392,6.34];
ybs[7257]=['',5.0294909,-0.5142654,6.3];
ybs[7258]=['',5.0371107,-0.8804943,6.13];
ybs[7259]=['',5.01721,0.6045172,6.74];
ybs[7260]=['',5.0335951,-0.6552899,6.57];
ybs[7261]=['τ Pav',5.0560614,-1.2069089,6.27];
ybs[7262]=['',5.0131834,0.9156153,5.81];
ybs[7263]=['',5.0341796,-0.3773492,6.41];
ybs[7264]=['',5.0376776,-0.4514947,5.8];
ybs[7265]=['',5.0586034,-1.162761,5.53];
ybs[7266]=['20 Aql',5.034592,-0.1379129,5.34];
ybs[7267]=['',5.0282384,0.4672734,6.36];
ybs[7268]=['',5.0449591,-0.7880972,5.92];
ybs[7269]=['',5.0372776,-0.2137084,5.51];
ybs[7270]=['19 Lyr',5.0291304,0.5466443,5.98];
ybs[7271]=['',5.0269785,0.7062654,6.18];
ybs[7272]=['',5.0332566,0.2946797,6.73];
ybs[7273]=['',5.0332408,0.3768509,5.93];
ybs[7274]=['21 Aql',5.0387449,0.0406962,5.15];
ybs[7275]=['',5.0387348,0.0969347,6.49];
ybs[7276]=['',5.0523905,-0.7928482,5.4];
ybs[7277]=['55 Dra',5.0171215,1.1521705,6.25];
ybs[7278]=['',5.0477465,-0.4213242,6.25];
ybs[7279]=['ψ Sgr',5.0477333,-0.4401302,4.85];
ybs[7280]=['',5.0293349,0.8707681,6.75];
ybs[7281]=['',5.0293639,0.870802,6.57];
ybs[7282]=['53 Dra',5.0268925,0.9930242,5.12];
ybs[7283]=['',5.0524999,-0.5843726,6.25];
ybs[7284]=['',5.0608641,-0.9310671,6.38];
ybs[7285]=['η Lyr',5.037415,0.6838925,4.39];
ybs[7286]=['',5.0439072,0.3532906,6];
ybs[7287]=['',5.0453748,0.2639373,5.57];
ybs[7288]=['1 Sge',5.0449432,0.3712502,5.64];
ybs[7289]=['',5.0450838,0.5334647,5.85];
ybs[7290]=['22 Aql',5.0508965,0.085071,5.59];
ybs[7291]=['43 Sgr',5.0566189,-0.3300936,4.96];
ybs[7292]=['',5.0475649,0.4798732,6.54];
ybs[7293]=['1 Vul',5.0489813,0.3740168,4.77];
ybs[7294]=['',5.0502461,0.2545415,5.63];
ybs[7295]=['',5.0477187,0.4881009,6.16];
ybs[7296]=['54 Dra',5.0365976,1.007806,4.99];
ybs[7297]=['δ Dra',5.0289421,1.1815692,3.07];
ybs[7298]=['',5.0434869,0.8745776,6.27];
ybs[7299]=['59 Dra',5.0106157,1.3368512,5.13];
ybs[7300]=['',5.0566077,0.0361596,6.19];
ybs[7301]=['θ Lyr',5.0488667,0.6662436,4.36];
ybs[7302]=['ω1 Aql',5.0563317,0.2030759,5.28];
ybs[7303]=['',5.0662028,-0.6175022,5.59];
ybs[7304]=['',5.0624522,-0.2704499,6.06];
ybs[7305]=['2 Vul',5.0555061,0.4025707,5.43];
ybs[7306]=['23 Aql',5.0598618,0.0196484,5.1];
ybs[7307]=['',5.0890122,-1.1925421,6.34];
ybs[7308]=['24 Aql',5.0612255,0.006629,6.41];
ybs[7309]=['',5.0504317,0.8209804,6];
ybs[7310]=['',5.069383,-0.5546008,6.58];
ybs[7311]=['',5.05641,0.5421409,6.68];
ybs[7312]=['',5.061035,0.168576,6.32];
ybs[7313]=['',5.0603623,0.3429771,6.58];
ybs[7314]=['',5.0698403,-0.3902719,5.58];
ybs[7315]=['κ Cyg',5.0509706,0.9321498,3.77];
ybs[7316]=['η Tel',5.0814314,-0.949125,5.05];
ybs[7317]=['',5.0741684,-0.6098506,6.48];
ybs[7318]=['28 Aql',5.0643262,0.2166955,5.53];
ybs[7319]=['ω2 Aql',5.0653528,0.2020415,6.02];
ybs[7320]=['26 Aql',5.0688406,-0.0938,5.01];
ybs[7321]=['',5.0774466,-0.7325801,6.34];
ybs[7322]=['',5.06088,0.5834559,6.6];
ybs[7323]=['27 Aql',5.0688889,-0.0148479,5.49];
ybs[7324]=['β1 Sgr',5.0796891,-0.7752104,4.01];
ybs[7325]=['',5.060481,0.6542528,6.22];
ybs[7326]=['',5.0740058,-0.3349658,6.26];
ybs[7327]=['ρ1 Sgr',5.0741931,-0.3107587,3.93];
ybs[7328]=['',5.0579681,0.8658601,6.31];
ybs[7329]=['υ Sgr',5.0743552,-0.2777328,4.61];
ybs[7330]=['β2 Sgr',5.0822481,-0.7811541,4.29];
ybs[7331]=['ρ2 Sgr',5.074974,-0.318805,5.87];
ybs[7332]=['',5.0632303,0.6522558,6.31];
ybs[7333]=['',5.0672869,0.6148361,6.31];
ybs[7334]=['',5.0767757,-0.1423969,6.31];
ybs[7335]=['α Sgr',5.0848912,-0.7081309,3.97];
ybs[7336]=['',5.076561,-0.0036676,5.83];
ybs[7337]=['',5.0871328,-0.7623531,6.17];
ybs[7338]=['',5.0618116,0.9497554,6.26];
ybs[7339]=['τ Dra',5.0402242,1.2809705,4.45];
ybs[7340]=['',5.0799401,-0.1284135,6.32];
ybs[7341]=['',5.0781673,0.1737584,6.35];
ybs[7342]=['',5.0869334,-0.4855871,6.04];
ybs[7343]=['',5.0643348,1.0068184,5.91];
ybs[7344]=['',5.0794284,0.261168,6.64];
ybs[7345]=['3 Vul',5.0777267,0.4591097,5.18];
ybs[7346]=['',5.0761234,0.5857402,6.06];
ybs[7347]=['',5.0894684,-0.510778,5.93];
ybs[7348]=['',5.0611517,1.1245452,6.52];
ybs[7349]=['χ1 Sgr',5.0901649,-0.426991,5.03];
ybs[7350]=['χ3 Sgr',5.0911021,-0.4174528,5.43];
ybs[7351]=['',5.082043,0.3544321,6.4];
ybs[7352]=['',5.0693656,1.0089518,6.43];
ybs[7353]=['',5.0883613,-0.0844826,6.52];
ybs[7354]=['',5.0901304,-0.241782,5.69];
ybs[7355]=['',5.0805469,0.5805856,6.37];
ybs[7356]=['2 Sge',5.0847295,0.2963759,6.25];
ybs[7357]=['',5.1029968,-0.9473668,5.69];
ybs[7358]=['π Dra',5.0648291,1.147659,4.59];
ybs[7359]=['2 Cyg',5.08316,0.5177441,4.97];
ybs[7360]=['31 Aql',5.0875363,0.2092311,5.16];
ybs[7361]=['',5.0843048,0.4909797,6.53];
ybs[7362]=['50 Sgr',5.0946134,-0.3793007,5.59];
ybs[7363]=['',5.0827192,0.6369592,6.36];
ybs[7364]=['δ Aql',5.0901474,0.0551282,3.36];
ybs[7365]=['',5.0937713,-0.261953,5.72];
ybs[7366]=['',5.0947346,-0.2531905,6.7];
ybs[7367]=['',5.0976582,-0.5183397,5.67];
ybs[7368]=['',5.0787338,0.8781472,6.51];
ybs[7369]=['',5.0816176,0.7580155,5.84];
ybs[7370]=['',5.1199092,-1.193578,5.96];
ybs[7371]=['',5.0889895,0.3545668,6.31];
ybs[7372]=['4 Vul',5.0894587,0.3463162,5.16];
ybs[7373]=['',5.0890519,0.4355746,6.19];
ybs[7374]=['ν Aql',5.0946923,0.0066846,4.66];
ybs[7375]=['',5.1121539,-0.966829,6.13];
ybs[7376]=['',5.0937458,0.228083,5.74];
ybs[7377]=['5 Vul',5.0926916,0.3515439,5.63];
ybs[7378]=['',5.093827,0.3479439,5.81];
ybs[7379]=['',5.109097,-0.757472,5.71];
ybs[7380]=['μ Tel',5.1151632,-0.9610394,6.3];
ybs[7381]=['λ UMi',4.4103958,1.5532759,6.38];
ybs[7382]=['4 Cyg',5.091676,0.634635,5.15];
ybs[7383]=['',5.0987776,0.2500598,6.32];
ybs[7384]=['',5.1025823,0.0519331,5.85];
ybs[7385]=['',5.1103098,-0.4701828,5.52];
ybs[7386]=['',5.1121724,-0.5593076,6.6];
ybs[7387]=['35 Aql',5.105539,0.0348346,5.8];
ybs[7388]=['',5.088359,1.0135315,6.6];
ybs[7389]=['',5.1073341,-0.1221402,6.61];
ybs[7390]=['',5.0979491,0.6629798,6.34];
ybs[7391]=['',5.1068327,0.0050936,6.25];
ybs[7392]=['α Vul',5.1033571,0.4312778,4.44];
ybs[7393]=['8 Vul',5.1044223,0.4330882,5.81];
ybs[7394]=['',5.1066364,0.2555436,5.56];
ybs[7395]=['ι1 Cyg',5.0961783,0.913946,5.75];
ybs[7396]=['7 Vul',5.1063384,0.3547458,6.33];
ybs[7397]=['',5.1145732,-0.371156,6.13];
ybs[7398]=['',5.1250939,-0.927437,5.75];
ybs[7399]=['',5.1105611,0.0514928,6.09];
ybs[7400]=['',5.0906271,1.0925999,6.38];
ybs[7401]=['36 Aql',5.1128835,-0.0478655,5.03];
ybs[7402]=['',5.1121866,0.0609255,6.05];
ybs[7403]=['',5.1264889,-0.7893104,5.61];
ybs[7404]=['β1 Cyg',5.1120216,0.4887984,3.08];
ybs[7405]=['β2 Cyg',5.112167,0.4888957,5.11];
ybs[7406]=['',5.1118945,0.6331178,6.25];
ybs[7407]=['ι2 Cyg',5.1061604,0.9036533,3.79];
ybs[7408]=['',5.1148798,0.4653728,5.87];
ybs[7409]=['',5.1295731,-0.6978974,5.89];
ybs[7410]=['',5.0628458,1.3900536,6.05];
ybs[7411]=['ι Tel',5.1347833,-0.838639,4.9];
ybs[7412]=['',5.0276665,1.4573606,6.53];
ybs[7413]=['8 Cyg',5.1163064,0.602137,4.74];
ybs[7414]=['',5.1133425,0.87883,5.53];
ybs[7415]=['',5.1124207,0.9735176,6.37];
ybs[7416]=['μ Aql',5.1274936,0.129624,4.45];
ybs[7417]=['37 Aql',5.1325987,-0.1834643,5.12];
ybs[7418]=['51 Sgr',5.1370761,-0.4305754,5.65];
ybs[7419]=['',5.1341289,-0.1293561,6.34];
ybs[7420]=['',5.1345697,-0.2130003,6.27];
ybs[7421]=['',5.1497154,-1.0111228,6.18];
ybs[7422]=['',5.1573345,-1.162992,6.39];
ybs[7423]=['',5.1240935,0.6773613,6.61];
ybs[7424]=['9 Vul',5.1291992,0.3459516,5];
ybs[7425]=['',5.1334474,0.0516967,6.38];
ybs[7426]=['',5.1386235,-0.3281845,6.11];
ybs[7427]=['52 Sgr',5.1400414,-0.4334399,4.6];
ybs[7428]=['9 Cyg',5.1299703,0.515071,5.38];
ybs[7429]=['',5.1237851,0.8606257,5.96];
ybs[7430]=['',5.141305,-0.3173292,5.64];
ybs[7431]=['',5.1285974,0.7410843,5.35];
ybs[7432]=['',5.1362813,0.1954592,6.68];
ybs[7433]=['κ Aql',5.1402004,-0.1217911,4.95];
ybs[7434]=['ι Aql',5.1392669,-0.0215914,4.36];
ybs[7435]=['',5.1203736,1.0507933,6.29];
ybs[7436]=['',5.1367352,0.252038,6.38];
ybs[7437]=['',5.1086739,1.2398066,6.07];
ybs[7438]=['',5.1264042,0.8950866,5.73];
ybs[7439]=['',5.1358858,0.3950519,6.32];
ybs[7440]=['',5.1281132,0.8414741,6.67];
ybs[7441]=['',5.1434247,-0.2487434,5.47];
ybs[7442]=['',5.1648152,-1.1484672,6.09];
ybs[7443]=['',5.1395137,0.197618,5.98];
ybs[7444]=['11 Cyg',5.1337893,0.6456535,6.05];
ybs[7445]=['',5.1381111,0.3557327,7.14];
ybs[7446]=['',5.1574978,-0.9488721,6.26];
ybs[7447]=['42 Aql',5.1440328,-0.0802449,5.46];
ybs[7448]=['',5.1541287,-0.7893643,6.25];
ybs[7449]=['σ Dra',5.1150023,1.2166346,4.68];
ybs[7450]=['ε Sge',5.1411419,0.2881942,5.66];
ybs[7451]=['',5.1547703,-0.6873532,6.61];
ybs[7452]=['',5.1334831,0.8776807,6.52];
ybs[7453]=['',5.1400528,0.512831,6.43];
ybs[7454]=['',5.138683,0.6707859,6.5];
ybs[7455]=['',5.136946,0.7809327,5.17];
ybs[7456]=['θ Cyg',5.1357121,0.8773795,4.48];
ybs[7457]=['53 Sgr',5.1535771,-0.4080051,6.34];
ybs[7458]=['',5.1482532,0.0598988,6.35];
ybs[7459]=['',5.1453541,0.3636005,6.48];
ybs[7460]=['',5.1548639,-0.4080172,5.97];
ybs[7461]=['σ Aql',5.1498297,0.0950896,5.17];
ybs[7462]=['',5.1504488,0.2901074,6.38];
ybs[7463]=['54 Sgr',5.1572365,-0.2834784,6.2];
ybs[7464]=['',5.142353,0.8610441,6.47];
ybs[7465]=['φ Cyg',5.1497081,0.5271561,4.69];
ybs[7466]=['α Sge',5.1533353,0.3152894,4.37];
ybs[7467]=['45 Aql',5.1566979,-0.009947,5.67];
ybs[7468]=['',5.1511617,0.5939325,6.1];
ybs[7469]=['',5.1548808,0.358276,6.5];
ybs[7470]=['14 Cyg',5.1493365,0.7482021,5.4];
ybs[7471]=['',5.1450891,0.9603494,5.82];
ybs[7472]=['',5.1555856,0.4148406,6.64];
ybs[7473]=['',5.1578277,0.2420229,6.01];
ybs[7474]=['',5.1497078,0.8030015,6.2];
ybs[7475]=['β Sge',5.1575084,0.3059112,4.37];
ybs[7476]=['55 Sgr',5.1650606,-0.2805062,5.06];
ybs[7477]=['',5.1581956,0.3927719,6.36];
ybs[7478]=['',5.1707941,-0.6542586,6.16];
ybs[7479]=['',5.1547573,0.7527405,6.16];
ybs[7480]=['46 Aql',5.1627792,0.2137191,6.34];
ybs[7481]=['',5.2353639,-1.4187922,6.39];
ybs[7482]=['',5.1552545,0.7954536,5.06];
ybs[7483]=['',5.1695741,-0.2690852,5.49];
ybs[7484]=['χ Aql',5.1643337,0.2073225,5.27];
ybs[7485]=['',5.200499,-1.2644521,5.41];
ybs[7486]=['',5.1604773,0.703465,6.23];
ybs[7487]=['',5.1511264,1.0569363,6.51];
ybs[7488]=['',5.1647688,0.5128438,6.49];
ybs[7489]=['',5.1643032,0.566861,5.94];
ybs[7490]=['16 Cyg',5.1591577,0.8827327,5.96];
ybs[7491]=['',5.1593841,0.8825974,6.2];
ybs[7492]=['',5.1661941,0.5363552,6.05];
ybs[7493]=['10 Vul',5.1688318,0.4507224,5.49];
ybs[7494]=['',5.1809577,-0.5559724,5.52];
ybs[7495]=['',5.1697218,0.4745236,6.28];
ybs[7496]=['',5.1598135,0.9689197,6.48];
ybs[7497]=['ν Tel',5.1913124,-0.9827555,5.35];
ybs[7498]=['ψ Aql',5.1730186,0.2331018,6.26];
ybs[7499]=['',5.1690825,0.597166,6.05];
ybs[7500]=['',5.2009865,-1.1651308,6.45];
ybs[7501]=['',5.1682135,0.7299941,5.84];
ybs[7502]=['56 Sgr',5.1819566,-0.3439563,4.86];
ybs[7503]=['',5.1792353,-0.049388,6.48];
ybs[7504]=['15 Cyg',5.1707573,0.6528792,4.89];
ybs[7505]=['',5.1734737,0.5116918,6.82];
ybs[7506]=['υ Aql',5.1779969,0.1338114,5.91];
ybs[7507]=['',5.1724759,0.60156,6.57];
ybs[7508]=['',5.194891,-0.9221081,6.25];
ybs[7509]=['',5.1646677,1.0134883,6.22];
ybs[7510]=['',5.1729327,0.7115654,6.34];
ybs[7511]=['',5.2056821,-1.1440422,6.05];
ybs[7512]=['γ Aql',5.1804931,0.186176,2.72];
ybs[7513]=['',5.1665995,0.9964942,6.27];
ybs[7514]=['',5.2020784,-1.0647472,6.21];
ybs[7515]=['δ Cyg',5.1733542,0.7886081,2.87];
ybs[7516]=['',5.1768617,0.6308413,6.43];
ybs[7517]=['',5.1777736,0.6120225,6.09];
ybs[7518]=['',5.2035131,-1.0321357,5.42];
ybs[7519]=['',5.1891064,-0.2382146,6.11];
ybs[7520]=['',5.1817058,0.4396103,6.62];
ybs[7521]=['17 Cyg',5.1803238,0.5895998,4.99];
ybs[7522]=['',5.181047,0.5749549,6.18];
ybs[7523]=['δ Sge',5.1851419,0.3244296,3.82];
ybs[7524]=['',5.2001858,-0.8290573,5.94];
ybs[7525]=['',5.1946603,-0.5014974,6.05];
ybs[7526]=['',5.1867231,0.4439829,5.95];
ybs[7527]=['',5.1933145,-0.1887702,6.04];
ybs[7528]=['',5.190295,0.1876051,6.44];
ybs[7529]=['',5.1846166,0.6712845,5.77];
ybs[7530]=['π Aql',5.1911083,0.2071834,5.72];
ybs[7531]=['',5.1673453,1.2110762,5.92];
ybs[7532]=['ζ Sge',5.1920598,0.335055,5];
ybs[7533]=['',5.1839579,0.837095,6.12];
ybs[7534]=['',5.2112432,-0.9584343,5.74];
ybs[7535]=['',5.2113454,-0.958531,6.5];
ybs[7536]=['',5.190312,0.6172576,6.53];
ybs[7537]=['',5.1908887,0.5845482,6.44];
ybs[7538]=['',5.206743,-0.6949549,5.33];
ybs[7539]=['51 Aql',5.2009162,-0.1868848,5.39];
ybs[7540]=['',5.1981776,0.1388958,6.51];
ybs[7541]=['',5.1933121,0.6765803,6.11];
ybs[7542]=['',5.1957659,0.4973487,6.38];
ybs[7543]=['α Aql',5.2002914,0.1557567,0.77];
ybs[7544]=['',5.2208727,-1.0666227,6.24];
ybs[7545]=['',5.2024116,-0.0419709,6.13];
ybs[7546]=['ο Aql',5.2013019,0.1827628,5.11];
ybs[7547]=['57 Sgr',5.2073964,-0.3314106,5.92];
ybs[7548]=['',5.202499,0.1690593,6.25];
ybs[7549]=['',5.1782299,1.1954126,6.34];
ybs[7550]=['χ Cyg',5.1984298,0.5754331,4.23];
ybs[7551]=['12 Vul',5.201054,0.3955959,4.95];
ybs[7552]=['19 Cyg',5.1981511,0.6768074,5.12];
ybs[7553]=['',5.198289,0.7095714,5.69];
ybs[7554]=['',5.1991442,0.6611691,6.06];
ybs[7555]=['',5.2057637,0.2039477,6.13];
ybs[7556]=['η Aql',5.2079214,0.0185393,3.9];
ybs[7557]=['',5.2111957,-0.2538769,6.48];
ybs[7558]=['',5.2066864,0.1816528,6.54];
ybs[7559]=['',5.2051456,0.4371812,5.57];
ybs[7560]=['9 Sge',5.2068515,0.3268744,6.23];
ybs[7561]=['',5.2117115,-0.0533615,5.65];
ybs[7562]=['20 Cyg',5.1974487,0.9257879,5.03];
ybs[7563]=['',5.2009288,0.8278664,6.2];
ybs[7564]=['',5.2167153,-0.4168469,6.18];
ybs[7565]=['',5.2398215,-1.2060948,5.75];
ybs[7566]=['',5.2117538,0.0778003,6.53];
ybs[7567]=['ι Sgr',5.2217569,-0.7297275,4.13];
ybs[7568]=['ε Dra',5.1839965,1.2273536,3.83];
ybs[7569]=['',5.2057138,0.6368483,6.1];
ybs[7570]=['56 Aql',5.2154876,-0.1486448,5.79];
ybs[7571]=['',5.2205381,-0.5757573,6.46];
ybs[7572]=['',5.2304378,-1.0099642,6.53];
ybs[7573]=['',5.2311676,-1.0269945,5.26];
ybs[7574]=['',5.240572,-1.1990829,6.39];
ybs[7575]=['',5.2038433,0.821768,5.62];
ybs[7576]=['ε Pav',5.2492242,-1.2714705,3.96];
ybs[7577]=['',5.2043692,0.8375545,5.91];
ybs[7578]=['13 Vul',5.2114453,0.4212666,4.58];
ybs[7579]=['57 Aql',5.217614,-0.1425856,5.71];
ybs[7580]=['57 Aql',5.2176507,-0.14276,6.49];
ybs[7581]=['ξ Aql',5.2154259,0.1486821,4.71];
ybs[7582]=['58 Aql',5.2178642,0.0057826,5.61];
ybs[7583]=['ω Sgr',5.2235452,-0.4579951,4.7];
ybs[7584]=['',5.2173163,0.1256278,6.15];
ybs[7585]=['',5.2206194,-0.1165214,6.51];
ybs[7586]=['',5.2083094,0.8353898,6.29];
ybs[7587]=['',5.2160487,0.425459,5.52];
ybs[7588]=['β Aql',5.2201402,0.1128289,3.71];
ybs[7589]=['μ1 Pav',5.2466988,-1.1674325,5.76];
ybs[7590]=['59 Sgr',5.2284062,-0.4731804,4.52];
ybs[7591]=['',5.2321299,-0.6632113,6.55];
ybs[7592]=['',5.2167476,0.6467103,5.76];
ybs[7593]=['',5.2183862,0.5280161,6.57];
ybs[7594]=['23 Cyg',5.208627,1.0049693,5.14];
ybs[7595]=['10 Sge',5.2228908,0.2913476,5.36];
ybs[7596]=['φ Aql',5.224013,0.2004033,5.28];
ybs[7597]=['',5.2096788,1.0431118,6.06];
ybs[7598]=['μ2 Pav',5.2531828,-1.1673287,5.31];
ybs[7599]=['22 Cyg',5.2212856,0.6727339,4.94];
ybs[7600]=['61 Sgr',5.2323532,-0.2693426,5.02];
ybs[7601]=['η Cyg',5.2233981,0.6133382,3.89];
ybs[7602]=['',5.2270205,0.3675097,6.48];
ybs[7603]=['',5.2372378,-0.5319481,6.28];
ybs[7604]=['60 Sgr',5.2371187,-0.4561573,4.83];
ybs[7605]=['ψ Cyg',5.2193491,0.9162434,4.92];
ybs[7606]=['',5.2252194,0.6337182,6.02];
ybs[7607]=['',5.2447342,-0.8602853,6.17];
ybs[7608]=['11 Sge',5.2304582,0.2940569,5.53];
ybs[7609]=['θ1 Sgr',5.2409238,-0.6146411,4.37];
ybs[7610]=['θ2 Sgr',5.2414127,-0.6045416,5.3];
ybs[7611]=['',5.2514419,-1.0352433,5.13];
ybs[7612]=['',5.2176274,1.0176686,6.09];
ybs[7613]=['',5.2444058,-0.750194,6.14];
ybs[7614]=['',5.2271679,0.7055811,5.45];
ybs[7615]=['',5.2433455,-0.6569757,5.95];
ybs[7616]=['',5.2461155,-0.7863145,5.81];
ybs[7617]=['',5.2434743,-0.5871912,5.66];
ybs[7618]=['',5.2243771,0.8894373,6.43];
ybs[7619]=['',5.219983,1.0280722,4.96];
ybs[7620]=['',5.2219447,0.990391,6.12];
ybs[7621]=['γ Sge',5.2347281,0.341242,3.47];
ybs[7622]=['',5.2380388,0.0250907,6.17];
ybs[7623]=['',5.240196,-0.1727582,5.88];
ybs[7624]=['',5.2301402,0.7386218,6.43];
ybs[7625]=['',5.2486312,-0.71128,6.29];
ybs[7626]=['',5.2337424,0.5418031,5.49];
ybs[7627]=['14 Vul',5.2364225,0.404237,5.67];
ybs[7628]=['',5.2331467,0.6661037,6.32];
ybs[7629]=['',5.2476503,-0.3957788,6.01];
ybs[7630]=['',5.2692471,-1.1738741,6.07];
ybs[7631]=['13 Sge',5.2404645,0.3067723,5.37];
ybs[7632]=['',5.2360612,0.799918,5.92];
ybs[7633]=['25 Cyg',5.2390809,0.6475654,5.19];
ybs[7634]=['',5.244807,0.1504226,5.91];
ybs[7635]=['63 Sgr',5.2498557,-0.2369448,5.71];
ybs[7636]=['62 Sgr',5.253336,-0.4825553,4.58];
ybs[7637]=['',5.235227,0.9095866,6.15];
ybs[7638]=['',5.2577155,-0.6611144,4.77];
ybs[7639]=['15 Vul',5.2446533,0.4854486,4.64];
ybs[7640]=['',5.2305048,1.109914,5.96];
ybs[7641]=['',5.2449078,0.6485553,6.2];
ybs[7642]=['',5.2475774,0.4339084,5.88];
ybs[7643]=['16 Vul',5.2487872,0.4363153,5.22];
ybs[7644]=['',5.2578508,-0.3932879,6.45];
ybs[7645]=['',5.2607916,-0.5584057,4.99];
ybs[7646]=['26 Cyg',5.2445726,0.8755501,5.05];
ybs[7647]=['',5.2585698,-0.1292906,6.72];
ybs[7648]=['',5.2544779,0.3239697,5.96];
ybs[7649]=['',5.2812185,-1.1569909,6.45];
ybs[7650]=['',5.2555519,0.2808764,5.67];
ybs[7651]=['δ Pav',5.2828614,-1.1539724,3.56];
ybs[7652]=['',5.2300194,1.2291689,6.33];
ybs[7653]=['62 Aql',5.2599574,-0.0112988,5.68];
ybs[7654]=['',5.2661015,-0.5748655,6.53];
ybs[7655]=['τ Aql',5.2586234,0.1281071,5.52];
ybs[7656]=['',5.2555749,0.5228717,5.71];
ybs[7657]=['',5.2633701,-0.2013594,6.34];
ybs[7658]=['15 Sge',5.2581499,0.2990083,5.8];
ybs[7659]=['ξ Tel',5.2753632,-0.9218361,4.94];
ybs[7660]=['',5.2764175,-0.9591068,6.26];
ybs[7661]=['65 Sgr',5.2649315,-0.219959,6.55];
ybs[7662]=['64 Dra',5.2433935,1.1323988,5.27];
ybs[7663]=['',5.2617446,0.406183,6.45];
ybs[7664]=['',5.2597476,0.5634047,5.64];
ybs[7665]=['η Sge',5.2626535,0.3499994,5.1];
ybs[7666]=['',5.2640445,0.271622,6.34];
ybs[7667]=['',5.2679887,-0.0700828,6.47];
ybs[7668]=['65 Dra',5.2472013,1.1291478,6.57];
ybs[7669]=['',5.261906,0.6726617,6.19];
ybs[7670]=['',5.2583406,0.8428497,6.16];
ybs[7671]=['ρ Dra',5.2486692,1.1856851,4.51];
ybs[7672]=['69 Dra',5.2315384,1.3358911,6.2];
ybs[7673]=['',5.2608268,0.9058559,6.14];
ybs[7674]=['17 Vul',5.2700821,0.4132518,5.07];
ybs[7675]=['27 Cyg',5.2672694,0.6289361,5.36];
ybs[7676]=['64 Aql',5.2758532,-0.0107276,5.99];
ybs[7677]=['',5.2920654,-1.0028436,6.37];
ybs[7678]=['',5.2614949,0.9844312,6.21];
ybs[7679]=['',5.274703,0.1651659,6.43];
ybs[7680]=['',5.278291,-0.174513,6.18];
ybs[7681]=['',5.2578422,1.1161834,6.26];
ybs[7682]=['',5.2756446,0.2919563,6.42];
ybs[7683]=['',5.2656018,0.9290143,5.85];
ybs[7684]=['',5.3636631,-1.4527929,6.17];
ybs[7685]=['',5.2731346,0.6019035,6.11];
ybs[7686]=['',5.2781521,0.188317,6.31];
ybs[7687]=['66 Dra',5.2616414,1.0831159,5.39];
ybs[7688]=['',5.2700415,0.8777675,6.54];
ybs[7689]=['',5.2909343,-0.6289463,5.32];
ybs[7690]=['',5.2576909,1.1883822,6.28];
ybs[7691]=['θ Sge',5.2835081,0.3661662,6.48];
ybs[7692]=['',5.296508,-0.7455151,6.22];
ybs[7693]=['',5.3073007,-1.1056515,6.09];
ybs[7694]=['28 Cyg',5.2806052,0.6440955,4.93];
ybs[7695]=['',5.2897986,-0.15319,6.49];
ybs[7696]=['θ Aql',5.2901473,-0.0131992,3.23];
ybs[7697]=['18 Vul',5.2859725,0.4706965,5.52];
ybs[7698]=['ξ1 Cap',5.2933863,-0.2151478,6.34];
ybs[7699]=['',5.2883676,0.3700047,6.22];
ybs[7700]=['',5.3055198,-0.914186,5.65];
ybs[7701]=['ξ2 Cap',5.2954292,-0.2190713,5.85];
ybs[7702]=['',5.2896218,0.3829369,6.26];
ybs[7703]=['',5.295688,0.0162872,6.27];
ybs[7704]=['19 Vul',5.2913971,0.4690431,5.49];
ybs[7705]=['20 Vul',5.2923337,0.4632851,5.92];
ybs[7706]=['66 Aql',5.2985596,-0.0164667,5.47];
ybs[7707]=['',5.291504,0.8343077,6.92];
ybs[7708]=['',5.3083937,-0.4706432,5.73];
ybs[7709]=['',5.2996864,0.4242025,6.56];
ybs[7710]=['ρ Aql',5.3026163,0.2664053,4.95];
ybs[7711]=['',5.3109227,-0.5225189,6.3];
ybs[7712]=['',5.2932711,0.8993536,6.01];
ybs[7713]=['68 Dra',5.2880174,1.0846121,5.75];
ybs[7714]=['',5.3136062,-0.6350737,6.39];
ybs[7715]=['',5.3137522,-0.6131404,6.53];
ybs[7716]=['30 Cyg',5.2969679,0.8182407,4.83];
ybs[7717]=['21 Vul',5.3019838,0.5019756,5.18];
ybs[7718]=['',5.3273046,-1.1023929,6.27];
ybs[7719]=['',5.2989838,0.7582629,6.14];
ybs[7720]=['',5.3009522,0.6400393,6.45];
ybs[7721]=['31 Cyg',5.2984216,0.8169439,3.79];
ybs[7722]=['29 Cyg',5.3029079,0.6435528,4.97];
ybs[7723]=['',5.3018804,0.7360052,6.71];
ybs[7724]=['3 Cap',5.3126437,-0.2141447,6.32];
ybs[7725]=['',5.3065706,0.4478298,4.78];
ybs[7726]=['33 Cyg',5.2966123,0.9884444,4.3];
ybs[7727]=['22 Vul',5.3076872,0.4114707,5.15];
ybs[7728]=['',5.2964294,1.0595278,5.79];
ybs[7729]=['',5.3068107,0.5898567,5.66];
ybs[7730]=['23 Vul',5.3086838,0.4866186,4.52];
ybs[7731]=['',5.3252889,-0.8315151,6.31];
ybs[7732]=['18 Sge',5.3113565,0.3781412,6.13];
ybs[7733]=['α1 Cap',5.3181818,-0.2171265,4.24];
ybs[7734]=['4 Cap',5.3201239,-0.3794681,5.87];
ybs[7735]=['',5.3268709,-0.8292338,6.13];
ybs[7736]=['κ Cep',5.2714347,1.3574306,4.39];
ybs[7737]=['32 Cyg',5.306396,0.8339409,3.98];
ybs[7738]=['',5.3094591,0.6800659,6.27];
ybs[7739]=['24 Vul',5.3132403,0.4317697,5.32];
ybs[7740]=['α2 Cap',5.319957,-0.2177585,3.57];
ybs[7741]=['',5.3073184,0.8778959,6.31];
ybs[7742]=['',5.3088852,0.7966823,5.91];
ybs[7743]=['',5.3113538,0.6479308,6.48];
ybs[7744]=['',5.3328464,-0.9596098,6.27];
ybs[7745]=['',5.3131623,0.7056802,5.24];
ybs[7746]=['',5.3163015,0.5099125,6.22];
ybs[7747]=['σ Cap',5.3260072,-0.3324844,5.28];
ybs[7748]=['',5.3154961,0.7468208,6.29];
ybs[7749]=['34 Cyg',5.3170572,0.6649867,4.81];
ybs[7750]=['',5.3310574,-0.508381,6.3];
ybs[7751]=['',5.3330576,-0.6214122,6.46];
ybs[7752]=['',5.3374286,-0.8714383,6.27];
ybs[7753]=['',5.3183626,0.7120984,5.84];
ybs[7754]=['',5.3268791,-0.017625,6.06];
ybs[7755]=['36 Cyg',5.3201199,0.6469617,5.58];
ybs[7756]=['35 Cyg',5.3209767,0.611756,5.17];
ybs[7757]=['',5.3254195,0.2318774,6.21];
ybs[7758]=['',5.3301495,-0.1098262,6.63];
ybs[7759]=['ν Cap',5.3313416,-0.2214819,4.76];
ybs[7760]=['',5.3276641,0.2376602,5.95];
ybs[7761]=['',5.3318998,-0.2568384,6.1];
ybs[7762]=['β Cap',5.3329247,-0.2567736,3.08];
ybs[7763]=['',5.3211438,0.8096771,6.45];
ybs[7764]=['',5.3291084,0.2554844,6.13];
ybs[7765]=['κ1 Sgr',5.3403153,-0.7326842,5.59];
ybs[7766]=['',5.3290697,0.311752,5.8];
ybs[7767]=['',5.3186456,0.9680523,5.76];
ybs[7768]=['',5.3259149,0.6492842,6.57];
ybs[7769]=['',5.3132097,1.1680008,5.93];
ybs[7770]=['',5.3277629,0.688921,6.23];
ybs[7771]=['',5.396261,-1.4117981,5.77];
ybs[7772]=['',5.3259533,0.8186689,6.5];
ybs[7773]=['κ2 Sgr',5.3465607,-0.7391847,5.64];
ybs[7774]=['',5.3414973,-0.1672817,6.3];
ybs[7775]=['25 Vul',5.336263,0.4278821,5.54];
ybs[7776]=['α Pav',5.3552602,-0.9889664,1.94];
ybs[7777]=['',5.327955,0.9366328,6.18];
ybs[7778]=['71 Dra',5.3231074,1.0877951,5.72];
ybs[7779]=['',5.3401516,0.2551929,6.17];
ybs[7780]=['',5.3417711,0.0944797,5.31];
ybs[7781]=['',5.3355219,0.7190945,6.39];
ybs[7782]=['γ Cyg',5.3363461,0.7038291,2.2];
ybs[7783]=['',5.3384781,0.5468982,6.09];
ybs[7784]=['',5.3354333,0.8004899,5.58];
ybs[7785]=['',5.3548074,-0.7107849,6.09];
ybs[7786]=['',5.3386069,0.7172623,5.93];
ybs[7787]=['',5.3527442,-0.4990261,5.85];
ybs[7788]=['',5.3392444,0.7514234,6.2];
ybs[7789]=['',5.3482102,0.0198874,6.15];
ybs[7790]=['',5.3240623,1.2033869,5.55];
ybs[7791]=['',5.3297581,1.1178747,5.69];
ybs[7792]=['39 Cyg',5.3438476,0.5630516,4.43];
ybs[7793]=['',5.3430843,0.6553154,5.9];
ybs[7794]=['',5.3593723,-0.6515523,6.25];
ybs[7795]=['',5.3530506,-0.0476295,6.11];
ybs[7796]=['',5.35278,0.1767616,6.33];
ybs[7797]=['',5.3521648,0.3749138,5.66];
ybs[7798]=['',5.4183509,-1.4174156,5.91];
ybs[7799]=['',5.3537225,0.347956,6.41];
ybs[7800]=['π Cap',5.3605497,-0.3165967,5.25];
ybs[7801]=['',5.3455865,0.9358918,6.51];
ybs[7802]=['',5.3554062,0.3034626,6.22];
ybs[7803]=['',5.3675519,-0.6199966,6.1];
ybs[7804]=['',5.347377,1.0414538,6.44];
ybs[7805]=['',5.3666067,-0.273477,6.41];
ybs[7806]=['',5.3632663,0.1485242,6.25];
ybs[7807]=['68 Aql',5.3648762,-0.05734,6.13];
ybs[7808]=['ρ Cap',5.3672537,-0.3096383,4.78];
ybs[7809]=['',5.3580229,0.6004061,6.39];
ybs[7810]=['',5.3641042,0.0525226,6.21];
ybs[7811]=['',5.3702897,-0.3895355,6.16];
ybs[7812]=['40 Cyg',5.3597738,0.6721664,5.62];
ybs[7813]=['',5.3534025,0.9897825,6.36];
ybs[7814]=['43 Cyg',5.3568208,0.8631544,5.69];
ybs[7815]=['ο Cap',5.3717053,-0.3231232,6.74];
ybs[7816]=['ο Cap',5.371807,-0.3230649,5.94];
ybs[7817]=['69 Aql',5.3702508,-0.0490891,4.91];
ybs[7818]=['',5.3767428,-0.5068257,6.39];
ybs[7819]=['',5.368259,0.3518685,6.55];
ybs[7820]=['41 Cyg',5.3680826,0.5313029,4.01];
ybs[7821]=['42 Cyg',5.3675888,0.6375249,5.88];
ybs[7822]=['1 Del',5.3726852,0.1914458,6.08];
ybs[7823]=['',5.3768138,-0.2614998,6.12];
ybs[7824]=['',5.4016171,-1.2136219,6.11];
ybs[7825]=['',5.3752993,0.3609218,6.18];
ybs[7826]=['',5.3766751,0.197818,7.11];
ybs[7827]=['',5.3699761,0.8028798,6.41];
ybs[7828]=['',5.3849932,-0.434056,6.36];
ybs[7829]=['',5.3668338,0.9798419,5.91];
ybs[7830]=['ω1 Cyg',5.3700468,0.8556424,4.95];
ybs[7831]=['',5.3824296,-0.1706799,5.65];
ybs[7832]=['ν Mic',5.3903943,-0.7756476,5.11];
ybs[7833]=['44 Cyg',5.3747588,0.6459337,6.19];
ybs[7834]=['φ1 Pav',5.3989129,-1.0560315,4.76];
ybs[7835]=['',5.3794952,0.4516618,6.34];
ybs[7836]=['θ Cep',5.3666105,1.1007256,4.22];
ybs[7837]=['ω2 Cyg',5.3755127,0.8603395,5.44];
ybs[7838]=['ε Del',5.3853882,0.1985789,4.03];
ybs[7839]=['',5.3944574,-0.6634791,6.44];
ybs[7840]=['',5.3754572,0.9142606,6.18];
ybs[7841]=['',5.3903967,-0.2381725,6.13];
ybs[7842]=['',5.3935516,-0.530554,6.4];
ybs[7843]=['',5.3883927,0.1768787,6.56];
ybs[7844]=['η Del',5.3885519,0.2286717,5.38];
ybs[7845]=['ρ Pav',5.4077687,-1.0725689,4.88];
ybs[7846]=['',5.3769234,0.9922843,6.14];
ybs[7847]=['',5.3826847,0.7551319,6.6];
ybs[7848]=['',5.3892492,0.3675673,6.48];
ybs[7849]=['μ1 Oct',5.4306732,-1.3282355,6];
ybs[7850]=['μ2 Oct',5.4289083,-1.3137518,6.55];
ybs[7851]=['',5.3963341,-0.2871145,6.19];
ybs[7852]=['47 Cyg',5.3875713,0.6165459,4.61];
ybs[7853]=['',5.386844,0.7303646,6.49];
ybs[7854]=['',5.3664673,1.2671879,6.27];
ybs[7855]=['α Ind',5.4064534,-0.8240593,3.11];
ybs[7856]=['',5.3870427,0.8162645,5.78];
ybs[7857]=['ζ Del',5.3944309,0.2574259,4.68];
ybs[7858]=['',5.4178178,-1.0965996,6.22];
ybs[7859]=['70 Aql',5.4011171,-0.0431821,4.89];
ybs[7860]=['26 Vul',5.3976917,0.4530537,6.41];
ybs[7861]=['φ2 Pav',5.4183381,-1.0554283,5.12];
ybs[7862]=['',5.3907379,0.9063346,6.11];
ybs[7863]=['',5.4067803,-0.4369004,6.36];
ybs[7864]=['',5.4035658,0.0030198,6.22];
ybs[7865]=['73 Dra',5.3721159,1.3094882,5.2];
ybs[7866]=['27 Vul',5.4017755,0.4631737,5.59];
ybs[7867]=['υ Pav',5.4275242,-1.1638332,5.15];
ybs[7868]=['β Del',5.4042178,0.2560649,3.63];
ybs[7869]=['ι Del',5.4054844,0.1999109,5.43];
ybs[7870]=['71 Aql',5.4081083,-0.0179556,4.32];
ybs[7871]=['48 Cyg',5.4035592,0.5523727,6.32];
ybs[7872]=['',5.4056768,0.3201887,6.25];
ybs[7873]=['',5.4036194,0.5514904,6.49];
ybs[7874]=['',5.402678,0.6702879,6.2];
ybs[7875]=['τ Cap',5.4125846,-0.2596671,5.22];
ybs[7876]=['',5.4119918,-0.0407696,6.22];
ybs[7877]=['29 Vul',5.408254,0.3713651,4.82];
ybs[7878]=['θ Del',5.4094137,0.2337281,5.72];
ybs[7879]=['',5.4178198,-0.5821477,5.47];
ybs[7880]=['28 Vul',5.4082017,0.4222413,5.04];
ybs[7881]=['',5.4084491,0.4146399,5.91];
ybs[7882]=['κ Del',5.4112465,0.1773762,5.05];
ybs[7883]=['1 Aqr',5.4127615,0.0098317,5.16];
ybs[7884]=['',5.4168811,-0.413584,6.37];
ybs[7885]=['',5.4108732,0.2777661,5.97];
ybs[7886]=['υ Cap',5.4160687,-0.315231,5.1];
ybs[7887]=['75 Dra',5.3529041,1.4223494,5.46];
ybs[7888]=['',5.4187447,-0.4636915,6.51];
ybs[7889]=['',5.4110936,0.3821228,6.08];
ybs[7890]=['',5.4099893,0.5307748,5.68];
ybs[7891]=['',5.4181517,-0.280069,5.8];
ybs[7892]=['α Del',5.4132935,0.2790596,3.77];
ybs[7893]=['',5.4144135,0.1976902,6.42];
ybs[7894]=['74 Dra',5.3586718,1.4165747,5.96];
ybs[7895]=['',5.4223911,-0.5501378,5.76];
ybs[7896]=['',5.422218,-0.4524286,6.28];
ybs[7897]=['',5.4120143,0.7095874,6.06];
ybs[7898]=['',5.4110035,0.7983796,6.58];
ybs[7899]=['β Pav',5.4405061,-1.1540776,3.42];
ybs[7900]=['',5.4180368,0.3492877,6.45];
ybs[7901]=['',5.4292123,-0.6890603,6.29];
ybs[7902]=['',5.4085973,0.9788094,6.48];
ybs[7903]=['',5.4170467,0.5215504,6.08];
ybs[7904]=['10 Del',5.4204554,0.2558775,5.99];
ybs[7905]=['',5.4140381,0.7598417,5.95];
ybs[7906]=['η Ind',5.4349525,-0.9048181,4.51];
ybs[7907]=['49 Cyg',5.4188732,0.5652205,5.51];
ybs[7908]=['',5.4184264,0.6834662,6.51];
ybs[7909]=['',5.4234226,0.3071659,6.22];
ybs[7910]=['α Cyg',5.419968,0.7916453,1.25];
ybs[7911]=['',5.4138023,1.0573628,6.01];
ybs[7912]=['',5.4223845,0.729457,5.67];
ybs[7913]=['',5.4245403,0.6201881,6.66];
ybs[7914]=['δ Del',5.4299927,0.264469,4.43];
ybs[7915]=['51 Cyg',5.4230568,0.8799592,5.39];
ybs[7916]=['',5.3525317,1.4607977,6.19];
ybs[7917]=['',5.4389101,-0.4741704,6.5];
ybs[7918]=['',5.4290479,0.6224933,6.47];
ybs[7919]=['',5.4442395,-0.6827632,5.5];
ybs[7920]=['σ Pav',5.459868,-1.198961,5.41];
ybs[7921]=['',5.4439971,-0.6290268,6.49];
ybs[7922]=['ψ Cap',5.4426584,-0.4396699,4.14];
ybs[7923]=['17 Cap',5.4428495,-0.3741034,5.93];
ybs[7924]=['',5.424128,1.0590567,6.15];
ybs[7925]=['30 Vul',5.4358582,0.4424342,4.91];
ybs[7926]=['',5.4269537,0.9981974,6.32];
ybs[7927]=['',5.4386838,0.3171189,6.38];
ybs[7928]=['52 Cyg',5.4391017,0.5375453,4.22];
ybs[7929]=['ι Mic',5.4538483,-0.7663401,5.11];
ybs[7930]=['',5.4320033,0.9872776,5.78];
ybs[7931]=['4 Cep',5.4255488,1.1647586,5.58];
ybs[7932]=['',5.4462002,-0.0420052,6.27];
ybs[7933]=['γ1 Del',5.4438716,0.2828167,5.14];
ybs[7934]=['γ2 Del',5.4439297,0.282812,4.27];
ybs[7935]=['ε Cyg',5.4413791,0.5942819,2.46];
ybs[7936]=['ε Aqr',5.4490829,-0.164334,3.77];
ybs[7937]=['3 Aqr',5.4492217,-0.0863514,4.42];
ybs[7938]=['ζ Ind',5.4583116,-0.8053996,4.89];
ybs[7939]=['13 Del',5.4492202,0.1062654,5.58];
ybs[7940]=['',5.4492575,0.0591124,6.4];
ybs[7941]=['',5.4362088,1.0063374,4.51];
ybs[7942]=['',5.4455922,0.6013376,4.92];
ybs[7943]=['η Cep',5.4354775,1.0806731,3.43];
ybs[7944]=['',5.442662,0.8135219,6.3];
ybs[7945]=['',5.4690254,-1.0881663,6.28];
ybs[7946]=['',5.4690617,-1.0881662,6.59];
ybs[7947]=['',5.4566214,-0.4485592,5.86];
ybs[7948]=['',5.4409768,0.9263309,6.33];
ybs[7949]=['λ Cyg',5.4465028,0.6382819,4.53];
ybs[7950]=['',5.4565899,-0.3133737,6.21];
ybs[7951]=['α Mic',5.4598499,-0.5881518,4.9];
ybs[7952]=['',5.4457965,0.7969122,6.4];
ybs[7953]=['',5.4309159,1.2187759,6.41];
ybs[7954]=['ι Ind',5.4674376,-0.8993087,5.05];
ybs[7955]=['',5.4477419,0.836224,5.57];
ybs[7956]=['',5.4633446,-0.5580345,6.36];
ybs[7957]=['',5.4645659,-0.6602897,5.52];
ybs[7958]=['',5.447694,0.9160778,6.27];
ybs[7959]=['15 Del',5.4570039,0.2203687,5.98];
ybs[7960]=['14 Del',5.4578976,0.1386692,6.33];
ybs[7961]=['',5.4587463,0.0981886,6.21];
ybs[7962]=['',5.4623416,-0.2175316,5.88];
ybs[7963]=['55 Cyg',5.4527215,0.8062507,4.84];
ybs[7964]=['',5.4513651,0.907415,6.29];
ybs[7965]=['β Mic',5.4685961,-0.5776227,6.04];
ybs[7966]=['ω Cap',5.4676834,-0.4684003,4.11];
ybs[7967]=['',5.4611517,0.316475,6.52];
ybs[7968]=['4 Aqr',5.4653437,-0.0967743,5.99];
ybs[7969]=['',5.4569438,0.8158032,6.33];
ybs[7970]=['56 Cyg',5.4578285,0.7703969,5.04];
ybs[7971]=['5 Aqr',5.4684671,-0.0946848,5.55];
ybs[7972]=['β Ind',5.4823972,-1.0187688,3.65];
ybs[7973]=['',5.4762346,-0.6933751,5.35];
ybs[7974]=['',5.4645535,0.4944895,5.77];
ybs[7975]=['',5.4728067,-0.4136569,6.33];
ybs[7976]=['μ Aqr',5.4707786,-0.1553557,4.73];
ybs[7977]=['',5.4747687,-0.5347024,6.35];
ybs[7978]=['',5.4808008,-0.8839198,6.24];
ybs[7979]=['',5.4526879,1.1191555,6.45];
ybs[7980]=['',5.4727688,-0.2005617,6.38];
ybs[7981]=['31 Vul',5.4674673,0.4743596,4.59];
ybs[7982]=['',5.4667242,0.5747538,6.44];
ybs[7983]=['',5.4777091,-0.4859499,6.41];
ybs[7984]=['',5.4764805,-0.1188067,6.44];
ybs[7985]=['',5.4717245,0.5189156,6.34];
ybs[7986]=['19 Cap',5.4803847,-0.311369,5.78];
ybs[7987]=['57 Cyg',5.4716406,0.7761388,4.78];
ybs[7988]=['76 Dra',5.4144222,1.4417946,5.75];
ybs[7989]=['',5.4718754,0.7900097,5.45];
ybs[7990]=['',5.4725846,0.7416359,6.66];
ybs[7991]=['',5.4749709,0.5850395,5.47];
ybs[7992]=['',5.4814027,-0.02252,6.55];
ybs[7993]=['',5.4772222,0.4992453,6.56];
ybs[7994]=['32 Vul',5.4780531,0.4911406,5.01];
ybs[7995]=['',5.4767357,0.7118454,6.7];
ybs[7996]=['',5.4836183,0.0805645,6.05];
ybs[7997]=['17 Del',5.4830796,0.2409354,5.17];
ybs[7998]=['16 Del',5.4832497,0.2208158,5.58];
ybs[7999]=['',5.4893185,-0.4574981,5.7];
ybs[8000]=['',5.4865649,-0.0607011,6.57];
ybs[8001]=['7 Aqr',5.489323,-0.1677924,5.51];
ybs[8002]=['',5.4388955,1.4072928,5.39];
ybs[8003]=['',5.4902629,0.009554,6.05];
ybs[8004]=['',5.4928843,-0.2783393,5.87];
ybs[8005]=['',5.5126498,-1.1889914,6.37];
ybs[8006]=['',5.4827639,0.8290488,5.67];
ybs[8007]=['α Oct',5.5293953,-1.3428052,5.15];
ybs[8008]=['',5.4842209,0.8928817,6.63];
ybs[8009]=['',5.4861704,0.7855467,5.96];
ybs[8010]=['',5.4972973,-0.2512996,6.01];
ybs[8011]=['',5.4851617,0.8868375,5.81];
ybs[8012]=['',5.485288,0.8600856,5.9];
ybs[8013]=['',5.5059886,-0.8932635,5.76];
ybs[8014]=['ν Cyg',5.4889543,0.7199651,3.94];
ybs[8015]=['',5.4840759,0.9943291,6.23];
ybs[8016]=['18 Del',5.4954668,0.1906495,5.48];
ybs[8017]=['',5.5036469,-0.6291011,6.11];
ybs[8018]=['33 Vul',5.49445,0.3911285,5.31];
ybs[8019]=['20 Cap',5.5013652,-0.3307497,6.25];
ybs[8020]=['ε Equ',5.4984397,0.0764124,5.23];
ybs[8021]=['',5.4938262,0.7776459,5.55];
ybs[8022]=['',5.4947728,0.7334661,6.16];
ybs[8023]=['',5.501479,0.2951165,6.66];
ybs[8024]=['',5.5026834,0.1326667,5.99];
ybs[8025]=['γ Mic',5.509154,-0.5615147,4.67];
ybs[8026]=['',5.4942562,0.8822017,5.61];
ybs[8027]=['11 Aqr',5.5051737,-0.0810744,6.21];
ybs[8028]=['',5.5135961,-0.7490296,6.64];
ybs[8029]=['',5.473589,1.3265929,6.05];
ybs[8030]=['',5.5040924,0.3388457,5.65];
ybs[8031]=['',5.5109925,-0.4676714,6.05];
ybs[8032]=['',5.5144582,-0.6709878,5.94];
ybs[8033]=['59 Cyg',5.5002246,0.8308784,4.74];
ybs[8034]=['ζ Mic',5.5166996,-0.6727492,5.3];
ybs[8035]=['',5.4975698,1.0388747,5.51];
ybs[8036]=['',5.5171892,-0.4825121,6.25];
ybs[8037]=['',5.5068285,0.6302621,5.97];
ybs[8038]=['',5.5467174,-1.3286188,6.58];
ybs[8039]=['60 Cyg',5.5062213,0.8070585,5.37];
ybs[8040]=['',5.5156785,-0.0146395,6.5];
ybs[8041]=['μ Ind',5.5274741,-0.9536547,5.16];
ybs[8042]=['',5.5158634,0.0282377,6.25];
ybs[8043]=['',5.5154326,0.2585868,6.31];
ybs[8044]=['12 Aqr',5.5205296,-0.1001295,7.31];
ybs[8045]=['12 Aqr',5.5205367,-0.1001246,5.89];
ybs[8046]=['η Cap',5.5223361,-0.3450259,4.84];
ybs[8047]=['',5.5481306,-1.2755677,5.68];
ybs[8048]=['',5.5116238,0.7832475,6.19];
ybs[8049]=['',5.514869,0.6762003,6.07];
ybs[8050]=['',5.5133463,0.8017117,6.48];
ybs[8051]=['',5.5097579,0.9905661,5.83];
ybs[8052]=['3 Equ',5.5224318,0.0975514,5.61];
ybs[8053]=['',5.5230057,0.0528572,6.42];
ybs[8054]=['',5.5232918,0.0411251,6.33];
ybs[8055]=['η Mic',5.5318883,-0.7208015,5.53];
ybs[8056]=['δ Mic',5.5296928,-0.5242609,5.68];
ybs[8057]=['',5.5181839,0.7280511,6.33];
ybs[8058]=['',5.5158204,0.8803086,6.37];
ybs[8059]=['',5.5427812,-1.1142329,5.76];
ybs[8060]=['',5.5172837,0.8193986,6.32];
ybs[8061]=['θ Cap',5.5289862,-0.2992499,4.07];
ybs[8062]=['',5.5314759,-0.5629465,5.18];
ybs[8063]=['4 Equ',5.5262023,0.1055077,5.94];
ybs[8064]=['',5.5171988,0.9315216,5.9];
ybs[8065]=['ξ Cyg',5.5227221,0.7681955,3.72];
ybs[8066]=['24 Cap',5.5343625,-0.4349078,4.5];
ybs[8067]=['',5.5565218,-1.2645796,6.2];
ybs[8068]=['',5.5297554,0.471441,6.12];
ybs[8069]=['',5.536828,-0.3031221,6.17];
ybs[8070]=['',5.5301106,0.5457975,5.82];
ybs[8071]=['61 Cyg',5.5315971,0.6777661,5.21];
ybs[8072]=['61 Cyg',5.5316482,0.6777225,6.03];
ybs[8073]=['χ Cap',5.5404934,-0.3683632,5.3];
ybs[8074]=['',5.5351816,0.2748226,6.34];
ybs[8075]=['63 Cyg',5.5298275,0.8331418,4.55];
ybs[8076]=['',5.5393776,0.1235229,6.15];
ybs[8077]=['27 Cap',5.5447926,-0.3572355,6.25];
ybs[8078]=['ο Pav',5.5646231,-1.2223698,5.02];
ybs[8079]=['ν Aqr',5.5447438,-0.1969317,4.51];
ybs[8080]=['',5.5395041,0.5287259,5.59];
ybs[8081]=['',5.5460313,0.0529143,6.45];
ybs[8082]=['',5.5498689,-0.1617076,6.27];
ybs[8083]=['γ Equ',5.547467,0.1783765,4.69];
ybs[8084]=['6 Equ',5.5482474,0.1769328,6.07];
ybs[8085]=['',5.5261867,1.2482408,5.87];
ybs[8086]=['',5.5571212,-0.7012763,5.83];
ybs[8087]=['',5.5479738,0.3934553,6.68];
ybs[8088]=['',5.5539543,-0.2510337,6.48];
ybs[8089]=['',5.5446906,0.7957109,6.63];
ybs[8090]=['',5.5606816,-0.6865379,5.26];
ybs[8091]=['',5.5498444,0.6350894,6.54];
ybs[8092]=['',5.5454534,0.9364004,5.73];
ybs[8093]=['',5.5469364,0.8339272,6.46];
ybs[8094]=['',5.5617128,-0.6341523,5.96];
ybs[8095]=['',5.5411519,1.1062543,6.54];
ybs[8096]=['',5.5612997,-0.480486,5.42];
ybs[8097]=['',5.5874406,-1.3134505,6.63];
ybs[8098]=['',5.5196536,1.3650734,5.91];
ybs[8099]=['',5.5405974,1.1969191,7.33];
ybs[8100]=['',5.5732114,-0.9280359,5.75];
ybs[8101]=['ζ Cyg',5.5582317,0.5291207,3.2];
ybs[8102]=['',5.5610151,0.2805118,6.27];
ybs[8103]=['',5.5702835,-0.7053937,6.21];
ybs[8104]=['',5.5651532,-0.1835272,6.77];
ybs[8105]=['',5.5516737,1.048513,5.64];
ybs[8106]=['',5.5601993,0.640941,6.05];
ybs[8107]=['',5.5663692,0.0031813,6.38];
ybs[8108]=['',5.5689526,-0.3011525,6.04];
ybs[8109]=['δ Equ',5.5655392,0.1762248,4.49];
ybs[8110]=['',5.5724481,-0.630419,6.12];
ybs[8111]=['',5.5840059,-1.1273142,6.31];
ybs[8112]=['',5.5636358,0.5234412,6.17];
ybs[8113]=['φ Cap',5.571316,-0.3588615,5.24];
ybs[8114]=['29 Cap',5.5716801,-0.263212,5.28];
ybs[8115]=['',5.6560677,-1.4785334,6.45];
ybs[8116]=['τ Cyg',5.5660576,0.6655922,3.72];
ybs[8117]=['α Equ',5.5715108,0.0931698,3.92];
ybs[8118]=['',5.5753253,-0.0264771,6.48];
ybs[8119]=['',5.5595198,1.1256241,6.39];
ybs[8120]=['',5.5780729,-0.2301728,6.4];
ybs[8121]=['ε Mic',5.5817154,-0.5599241,4.71];
ybs[8122]=['',5.569198,0.8388741,6.46];
ybs[8123]=['30 Cap',5.5813767,-0.3123105,5.43];
ybs[8124]=['',5.5733197,0.7390078,6.43];
ybs[8125]=['31 Cap',5.5827005,-0.3031797,7.05];
ybs[8126]=['θ Ind',5.5910824,-0.9312696,4.39];
ybs[8127]=['15 Aqr',5.582048,-0.0772862,5.82];
ybs[8128]=['',5.5858286,-0.500456,6.4];
ybs[8129]=['σ Cyg',5.5774728,0.6891552,4.23];
ybs[8130]=['',5.5772043,0.746552,6.19];
ybs[8131]=['',5.5918883,-0.7841806,6];
ybs[8132]=['υ Cyg',5.5798324,0.6106573,4.43];
ybs[8133]=['',5.5750408,0.944019,6.13];
ybs[8134]=['',5.5894916,-0.458345,6.56];
ybs[8135]=['',5.5846486,0.1971318,5.96];
ybs[8136]=['',5.5758075,0.9754458,5.98];
ybs[8137]=['θ1 Mic',5.5943267,-0.7106552,4.82];
ybs[8138]=['',5.5969952,-0.8699666,6.38];
ybs[8139]=['',5.5759067,1.0245528,6.42];
ybs[8140]=['68 Cyg',5.581809,0.7685931,5];
ybs[8141]=['',5.5839817,0.7178942,6.15];
ybs[8142]=['',5.6120576,-1.2154602,6.41];
ybs[8143]=['',5.5860518,0.6689695,5.83];
ybs[8144]=['',5.5903373,0.3860376,6.29];
ybs[8145]=['',5.6169325,-1.2515,6.09];
ybs[8146]=['16 Aqr',5.5946428,-0.0779771,5.87];
ybs[8147]=['',5.5860162,0.8657169,5.76];
ybs[8148]=['α Cep',5.581056,1.0939175,2.44];
ybs[8149]=['9 Equ',5.5944011,0.129969,5.82];
ybs[8150]=['',5.5844364,1.0247728,5.66];
ybs[8151]=['',5.5939547,0.4179722,5.57];
ybs[8152]=['',5.5926554,0.5680157,5.68];
ybs[8153]=['ι Cap',5.6000637,-0.2921996,4.28];
ybs[8154]=['',5.5651373,1.3456909,5.95];
ybs[8155]=['',5.5949723,0.5708114,6.04];
ybs[8156]=['',5.5931929,0.7057716,6.4];
ybs[8157]=['6 Cep',5.5842377,1.1338269,5.18];
ybs[8158]=['',5.603537,-0.3940253,5.6];
ybs[8159]=['1 Peg',5.5984882,0.3472682,4.08];
ybs[8160]=['',5.5518113,1.4193033,6.15];
ybs[8161]=['17 Aqr',5.6028978,-0.1610343,5.99];
ybs[8162]=['',5.6600664,-1.4414039,6.38];
ybs[8163]=['',5.6102084,-0.8119556,6.31];
ybs[8164]=['β Equ',5.6023277,0.1204965,5.16];
ybs[8165]=['',5.5899053,1.0620139,6.11];
ybs[8166]=['θ2 Mic',5.610253,-0.7140715,5.77];
ybs[8167]=['γ Pav',5.6207856,-1.1392115,4.22];
ybs[8168]=['',5.6008619,0.5306231,6.05];
ybs[8169]=['33 Cap',5.6085106,-0.3623072,5.41];
ybs[8170]=['',5.6084372,-0.3953813,6.38];
ybs[8171]=['',5.5970682,0.8636129,5.69];
ybs[8172]=['',5.600942,0.6759124,6.63];
ybs[8173]=['18 Aqr',5.6084495,-0.2231365,5.49];
ybs[8174]=['γ Ind',5.6189822,-0.952366,6.12];
ybs[8175]=['',5.6036162,0.654492,6.58];
ybs[8176]=['',5.606635,0.4252903,5.71];
ybs[8177]=['',5.6088681,0.1792016,6.35];
ybs[8178]=['20 Aqr',5.6111511,-0.0576804,6.36];
ybs[8179]=['',5.6054604,0.6535296,6.47];
ybs[8180]=['',5.607232,0.4434086,6.15];
ybs[8181]=['19 Aqr',5.6128536,-0.1685116,5.7];
ybs[8182]=['',5.6314334,-1.2114402,5.34];
ybs[8183]=['',5.6083965,0.4297384,6.32];
ybs[8184]=['',5.6091436,0.4584597,5.68];
ybs[8185]=['21 Aqr',5.6129944,-0.0604414,5.49];
ybs[8186]=['',5.6187093,-0.6586075,5.63];
ybs[8187]=['',5.6550949,-1.3952639,6.47];
ybs[8188]=['',5.6217018,-0.7409542,5.51];
ybs[8189]=['',5.6154103,0.0109651,6.46];
ybs[8190]=['ζ Cap',5.6194818,-0.3895103,3.74];
ybs[8191]=['',5.6180591,0.0208975,6.13];
ybs[8192]=['',5.6098296,0.8624856,6.58];
ybs[8193]=['35 Cap',5.6219792,-0.3682965,5.78];
ybs[8194]=['',5.6117185,0.8169542,5.6];
ybs[8195]=['69 Cyg',5.6141276,0.6416048,5.94];
ybs[8196]=['',5.6175271,0.3398075,6.07];
ybs[8197]=['',5.6308845,-0.9356876,6.39];
ybs[8198]=['',5.6260401,-0.2002547,6.61];
ybs[8199]=['36 Cap',5.6284339,-0.3789541,4.51];
ybs[8200]=['5 PsA',5.6301792,-0.5435609,6.5];
ybs[8201]=['70 Cyg',5.6209851,0.6494531,5.31];
ybs[8202]=['',5.6183144,0.8539735,5.31];
ybs[8203]=['35 Vul',5.6226476,0.4835082,5.41];
ybs[8204]=['',5.6176021,0.9248961,6.03];
ybs[8205]=['',5.6263864,0.1446909,6.4];
ybs[8206]=['',5.6245624,0.5640868,5.8];
ybs[8207]=['',5.628711,0.3141703,6.44];
ybs[8208]=['',5.6339047,-0.3325311,6.57];
ybs[8209]=['',5.6286031,0.3887589,5.93];
ybs[8210]=['',5.6200317,1.0444789,6.1];
ybs[8211]=['2 Peg',5.632713,0.4142362,4.57];
ybs[8212]=['',5.6267444,0.9688902,6.12];
ybs[8213]=['7 Cep',5.620748,1.167686,5.44];
ybs[8214]=['71 Cyg',5.6297496,0.8139424,5.24];
ybs[8215]=['ξ Gru',5.6437436,-0.7170393,5.29];
ybs[8216]=['6 PsA',5.6441237,-0.5907739,5.97];
ybs[8217]=['',5.6382844,0.213506,6.08];
ybs[8218]=['β Aqr',5.6404159,-0.0955651,2.91];
ybs[8219]=['',5.6495023,-0.9187635,6.41];
ybs[8220]=['',5.6788785,-1.3848208,6.18];
ybs[8221]=['',5.6452134,-0.4275114,6.43];
ybs[8222]=['',5.6495273,-0.7810761,5.57];
ybs[8223]=['',5.6332972,0.9259537,6.02];
ybs[8224]=['β Cep',5.6240104,1.2331696,3.23];
ybs[8225]=['',5.6030326,1.4070478,5.97];
ybs[8226]=['',5.6436502,0.4099837,6.7];
ybs[8227]=['',5.6533317,-0.747498,6.32];
ybs[8228]=['',5.6382079,0.9200597,6.16];
ybs[8229]=['',5.6355753,1.0568809,5.53];
ybs[8230]=['',5.6555034,-0.5166074,6.41];
ybs[8231]=['37 Cap',5.6551064,-0.3488525,5.69];
ybs[8232]=['',5.6448495,0.8739525,5.75];
ybs[8233]=['',5.6569997,-0.4076629,6.4];
ybs[8234]=['',5.6466019,0.8019839,6.25];
ybs[8235]=['',5.6711735,-1.1296896,6.2];
ybs[8236]=['',5.6529016,0.3988299,6.47];
ybs[8237]=['',5.656676,-0.0678281,5.77];
ybs[8238]=['ρ Cyg',5.6496004,0.7974109,4.02];
ybs[8239]=['8 PsA',5.661074,-0.4550825,5.73];
ybs[8240]=['ν Oct',5.6888111,-1.3489854,3.76];
ybs[8241]=['72 Cyg',5.6533528,0.674234,4.9];
ybs[8242]=['7 PsA',5.6640214,-0.5750995,6.11];
ybs[8243]=['',5.6560328,0.4938234,6.31];
ybs[8244]=['',5.6567152,0.4284616,6.11];
ybs[8245]=['',5.6513797,0.9039899,6.15];
ybs[8246]=['ε Cap',5.6648158,-0.3380486,4.68];
ybs[8247]=['',5.6599768,0.5262623,6.36];
ybs[8248]=['',5.6586065,0.7936308,5.53];
ybs[8249]=['',5.6665034,-0.0051102,6.25];
ybs[8250]=['ξ Aqr',5.667484,-0.1353786,4.69];
ybs[8251]=['3 Peg',5.6670746,0.1172139,6.18];
ybs[8252]=['74 Cyg',5.6627813,0.7070481,5.01];
ybs[8253]=['5 Peg',5.666914,0.3388755,5.45];
ybs[8254]=['',5.6740068,-0.5861024,6.28];
ybs[8255]=['',5.6786503,-0.9121248,6.21];
ybs[8256]=['4 Peg',5.6705985,0.102441,5.67];
ybs[8257]=['',5.681313,-0.9710848,6.33];
ybs[8258]=['',5.6648638,0.7818041,6.2];
ybs[8259]=['',5.6750309,-0.182891,6.08];
ybs[8260]=['',5.6710985,0.4467469,6.16];
ybs[8261]=['',5.6651772,0.9449155,6.15];
ybs[8262]=['',5.6723999,0.3554046,5.85];
ybs[8263]=['25 Aqr',5.6751385,0.0408702,5.1];
ybs[8264]=['γ Cap',5.6778841,-0.2890957,3.68];
ybs[8265]=['9 Cep',5.6657612,1.0852362,4.73];
ybs[8266]=['λ Oct',5.733288,-1.441948,5.29];
ybs[8267]=['',5.6706984,1.0050828,5.62];
ybs[8268]=['',5.6853957,-0.4363879,6.49];
ybs[8269]=['42 Cap',5.6841787,-0.2434527,5.18];
ybs[8270]=['75 Cyg',5.6768096,0.7569863,5.11];
ybs[8271]=['41 Cap',5.6864185,-0.4042872,5.24];
ybs[8272]=['',5.704371,-1.2375948,6.01];
ybs[8273]=['26 Aqr',5.6865618,0.0241577,5.67];
ybs[8274]=['κ Cap',5.6891273,-0.3275525,4.73];
ybs[8275]=['7 Peg',5.6868639,0.1008605,5.3];
ybs[8276]=['',5.6785822,0.9594178,6.2];
ybs[8277]=['76 Cyg',5.6829626,0.7139081,6.11];
ybs[8278]=['',5.6880288,0.1906542,6.09];
ybs[8279]=['',5.6916157,-0.3407171,6.22];
ybs[8280]=['',5.992125,-1.5478522,6.57];
ybs[8281]=['44 Cap',5.6908376,-0.2495924,5.88];
ybs[8282]=['',5.6850914,0.6214954,6.07];
ybs[8283]=['',5.6852525,0.800489,6.17];
ybs[8284]=['',5.6976492,-0.6711303,6.3];
ybs[8285]=['77 Cyg',5.6864923,0.7186587,5.69];
ybs[8286]=['π1 Cyg',5.6848,0.8951533,4.67];
ybs[8287]=['45 Cap',5.6949595,-0.2556914,5.99];
ybs[8288]=['',5.7016739,-0.8621715,6.45];
ybs[8289]=['',5.6873001,0.8674152,6.09];
ybs[8290]=['ι PsA',5.6994622,-0.5746695,4.34];
ybs[8291]=['',5.6896577,0.7200199,5.49];
ybs[8292]=['79 Cyg',5.6911626,0.6699112,5.65];
ybs[8293]=['ε Peg',5.6951924,0.1740869,2.39];
ybs[8294]=['μ1 Cyg',5.6945674,0.5033912,4.73];
ybs[8295]=['μ2 Cyg',5.6945456,0.5033961,6.08];
ybs[8296]=['46 Cap',5.6991504,-0.1567796,5.09];
ybs[8297]=['',5.6871756,1.0362032,6.08];
ybs[8298]=['9 Peg',5.6964503,0.3045517,4.34];
ybs[8299]=['',5.6965506,0.2595563,5.94];
ybs[8300]=['κ Peg',5.6968378,0.4493274,4.13];
ybs[8301]=['μ Cep',5.6904956,1.0276356,4.08];
ybs[8302]=['11 Cep',5.6820065,1.2463405,4.56];
ybs[8303]=['47 Cap',5.7046797,-0.1601475,6];
ybs[8304]=['λ Cap',5.7058715,-0.1966236,5.58];
ybs[8305]=['',5.7013471,0.6275697,6.4];
ybs[8306]=['12 Peg',5.7031405,0.4022786,5.29];
ybs[8307]=['δ Cap',5.7081747,-0.2797229,2.87];
ybs[8308]=['',5.7144069,-0.823847,5.58];
ybs[8309]=['',5.6868139,1.2639548,5.17];
ybs[8310]=['',5.7045038,0.4479108,6.28];
ybs[8311]=['θ PsA',5.7115636,-0.5375238,5.01];
ybs[8312]=['',5.6962296,1.0918803,5.95];
ybs[8313]=['11 Peg',5.708635,0.0486325,5.64];
ybs[8314]=['',5.7034409,0.7532989,6.54];
ybs[8315]=['',5.707659,0.3018448,6.21];
ybs[8316]=['',5.7231258,-1.1276803,5.62];
ybs[8317]=['',5.7105588,-0.1015219,6.17];
ybs[8318]=['ο Ind',5.7271667,-1.213493,5.53];
ybs[8319]=['ν Cep',5.6988028,1.0685007,4.29];
ybs[8320]=['π2 Cyg',5.7054369,0.8623602,4.23];
ybs[8321]=['',5.7117957,0.6402061,6.47];
ybs[8322]=['',5.7196389,-0.220296,6.31];
ybs[8323]=['',5.713267,0.6763022,6.12];
ybs[8324]=['12 Cep',5.7074766,1.0610396,5.52];
ybs[8325]=['',5.722051,-0.29223,6.38];
ybs[8326]=['',5.7179342,0.3588997,6.29];
ybs[8327]=['',5.7045668,1.2261111,6.29];
ybs[8328]=['14 Peg',5.7194326,0.5284021,5.04];
ybs[8329]=['13 Peg',5.7210453,0.303455,5.29];
ybs[8330]=['',5.7183323,0.719946,6.48];
ybs[8331]=['',5.7285357,-0.3232605,6.16];
ybs[8332]=['',5.7157224,1.0711718,6.17];
ybs[8333]=['',5.7272117,0.3478126,5.77];
ybs[8334]=['',5.7245626,0.6918144,6.17];
ybs[8335]=['',5.7303819,0.3730604,6.89];
ybs[8336]=['μ Cap',5.7354081,-0.2347405,5.08];
ybs[8337]=['',5.7454243,-1.0783306,5.9];
ybs[8338]=['γ Gru',5.7387305,-0.6503583,3.01];
ybs[8339]=['15 Peg',5.7310568,0.5043149,5.53];
ybs[8340]=['',5.7366667,-0.1781904,6.59];
ybs[8341]=['16 Peg',5.7335942,0.4542557,5.08];
ybs[8342]=['',5.7279332,0.9756139,5.71];
ybs[8343]=['',5.7361803,0.3450591,5.68];
ybs[8344]=['',5.7379228,0.1215908,6.15];
ybs[8345]=['',5.7390512,-0.0728475,5.71];
ybs[8346]=['',5.7253717,1.1493736,6.37];
ybs[8347]=['π Ind',5.749659,-1.0087406,6.19];
ybs[8348]=['',5.7408869,-0.0558286,6.2];
ybs[8349]=['',5.73909,0.3459348,6.39];
ybs[8350]=['',5.747254,-0.532389,6.41];
ybs[8351]=['',5.7494101,-0.6484027,5.46];
ybs[8352]=['',5.752286,-0.6570099,6.18];
ybs[8353]=['δ Ind',5.7568153,-0.9579972,4.4];
ybs[8354]=['κ1 Ind',5.7596024,-1.0281518,6.12];
ybs[8355]=['',5.7771953,-1.3536438,6.41];
ybs[8356]=['13 Cep',5.7404389,0.9898423,5.8];
ybs[8357]=['',5.7482701,0.372498,6.4];
ybs[8358]=['17 Peg',5.7508141,0.2125702,5.54];
ybs[8359]=['',5.7420892,1.0758987,6.13];
ybs[8360]=['',5.7424892,1.1418534,5.86];
ybs[8361]=['',5.7567343,-0.0928757,6.33];
ybs[8362]=['',5.750252,0.851225,6.42];
ybs[8363]=['',5.7592589,-0.3679031,6.12];
ybs[8364]=['',5.7621596,-0.668315,5.5];
ybs[8365]=['',5.7818405,-1.3266929,5.95];
ybs[8366]=['',5.7676728,-0.9735242,6.01];
ybs[8367]=['',5.7597469,-0.0745174,6.22];
ybs[8368]=['',5.747625,1.1122706,4.91];
ybs[8369]=['',5.7496974,1.1564395,6.43];
ybs[8370]=['18 Peg',5.7648457,0.1190547,6];
ybs[8371]=['η PsA',5.7686039,-0.4947934,5.42];
ybs[8372]=['ε Ind',5.7805982,-0.989277,4.69];
ybs[8373]=['',5.757511,1.096098,5.93];
ybs[8374]=['',5.7600477,1.0081358,6.59];
ybs[8375]=['28 Aqr',5.7691098,0.0123758,5.58];
ybs[8376]=['',5.7656685,0.5778787,6.46];
ybs[8377]=['20 Peg',5.7689092,0.2307988,5.6];
ybs[8378]=['19 Peg',5.7692764,0.1459325,5.65];
ybs[8379]=['',5.774316,-0.3106552,6.28];
ybs[8380]=['',5.7510414,1.3107386,6.35];
ybs[8381]=['29 Aqr',5.775366,-0.2942577,6.37];
ybs[8382]=['',5.7730235,0.1933512,6.37];
ybs[8383]=['',5.7793171,-0.5200995,7.1];
ybs[8384]=['',5.7652575,1.0924359,6.66];
ybs[8385]=['16 Cep',5.7576398,1.2790384,5.03];
ybs[8386]=['30 Aqr',5.7787922,-0.1120127,5.54];
ybs[8387]=['ο Aqr',5.7788962,-0.0357901,4.69];
ybs[8388]=['',5.7710783,0.9247882,5.78];
ybs[8389]=['21 Peg',5.7786611,0.2005565,5.8];
ybs[8390]=['13 PsA',5.7841629,-0.5203127,6.47];
ybs[8391]=['14 Cep',5.7718117,1.0141209,5.56];
ybs[8392]=['',5.7762374,0.7811139,5.6];
ybs[8393]=['',5.7850301,-0.4663084,5.96];
ybs[8394]=['κ2 Ind',5.7916026,-1.039008,5.62];
ybs[8395]=['32 Aqr',5.7853098,-0.0139912,5.3];
ybs[8396]=['λ Gru',5.7918928,-0.6883222,4.46];
ybs[8397]=['',5.7836974,0.5767772,6.38];
ybs[8398]=['ν Peg',5.7890893,0.0901263,4.84];
ybs[8399]=['α Aqr',5.7896338,-0.0037428,2.96];
ybs[8400]=['',5.7865371,0.4673818,5.78];
ybs[8401]=['18 Cep',5.7793322,1.1034751,5.29];
ybs[8402]=['ξ Cep',5.7787925,1.1297951,4.29];
ybs[8403]=['ι Aqr',5.792723,-0.240232,4.27];
ybs[8404]=['23 Peg',5.7881821,0.5073514,5.7];
ybs[8405]=['',5.8148842,-1.3225054,6.55];
ybs[8406]=['',5.7863399,0.8176839,6.13];
ybs[8407]=['',5.7888907,0.7891939,6.44];
ybs[8408]=['',5.7480394,1.4481485,6.98];
ybs[8409]=['',5.7897246,0.7874883,5.14];
ybs[8410]=['α Gru',5.8013599,-0.8177776,1.74];
ybs[8411]=['20 Cep',5.7842944,1.0976477,5.27];
ybs[8412]=['',5.7888189,0.8436385,6.27];
ybs[8413]=['19 Cep',5.7849492,1.0888248,5.11];
ybs[8414]=['',5.7904741,0.7915759,6.19];
ybs[8415]=['ι Peg',5.7945213,0.4441962,3.76];
ybs[8416]=['μ PsA',5.8016039,-0.573911,4.5];
ybs[8417]=['',5.8202047,-1.3266117,6.15];
ybs[8418]=['υ PsA',5.8018475,-0.5923289,4.99];
ybs[8419]=['',5.7900635,0.9852104,6.39];
ybs[8420]=['',5.7966698,0.3417571,5.75];
ybs[8421]=['',5.7967999,0.3160136,6.35];
ybs[8422]=['',5.8030381,-0.5762998,6.37];
ybs[8423]=['25 Peg',5.7982055,0.3806309,5.78];
ybs[8424]=['35 Aqr',5.8039151,-0.3213789,5.81];
ybs[8425]=['',5.80895,-0.8377737,6.43];
ybs[8426]=['',5.8000852,0.447668,6.11];
ybs[8427]=['',5.7940042,1.0288088,6.32];
ybs[8428]=['',5.7954717,0.9322303,6.14];
ybs[8429]=['',5.8083619,-0.5918135,5.37];
ybs[8430]=['',5.7993451,0.8709584,6.42];
ybs[8431]=['',5.8085505,-0.4919417,6.44];
ybs[8432]=['τ PsA',5.8092783,-0.5662193,4.92];
ybs[8433]=['',5.8012837,0.8001968,6.11];
ybs[8434]=['π1 Peg',5.8040202,0.5808163,5.58];
ybs[8435]=['θ Peg',5.8087955,0.1100278,3.53];
ybs[8436]=['',5.8096266,-0.0661091,6.27];
ybs[8437]=['38 Aqr',5.8109494,-0.1999892,5.46];
ybs[8438]=['',5.8105563,-0.0726193,6.01];
ybs[8439]=['π2 Peg',5.8073402,0.5809262,4.29];
ybs[8440]=['',5.8090672,0.3442368,6.18];
ybs[8441]=['',5.8093914,0.2571985,6.33];
ybs[8442]=['',5.8129343,-0.3687171,6.09];
ybs[8443]=['',5.8105502,0.2047427,5.78];
ybs[8444]=['28 Peg',5.8098563,0.3679934,6.46];
ybs[8445]=['',5.8112161,0.5351101,6.32];
ybs[8446]=['',5.815849,0.2818234,5.95];
ybs[8447]=['39 Aqr',5.8188676,-0.2458645,6.03];
ybs[8448]=['',5.8119408,0.8888942,5.4];
ybs[8449]=['',5.8213979,-0.4576386,6.17];
ybs[8450]=['ζ Cep',5.8102155,1.0176592,3.35];
ybs[8451]=['',5.8168964,0.4373236,5.92];
ybs[8452]=['',5.8200163,-0.0805274,6.39];
ybs[8453]=['24 Cep',5.8041447,1.2644437,4.79];
ybs[8454]=['λ Cep',5.8130141,1.0388386,5.04];
ybs[8455]=['',5.8247826,-0.4376175,5.58];
ybs[8456]=['ψ Oct',5.8462547,-1.3509449,5.51];
ybs[8457]=['',5.814488,0.9938977,5.24];
ybs[8458]=['',5.8061541,1.2604313,6.37];
ybs[8459]=['',5.8081996,1.2259047,5.5];
ybs[8460]=['',5.819582,0.6058329,5.33];
ybs[8461]=['',5.8149463,1.0330856,6.3];
ybs[8462]=['',5.8290987,-0.7203715,6.23];
ybs[8463]=['λ PsA',5.827334,-0.4827513,5.43];
ybs[8464]=['',5.8152004,1.0623154,5.35];
ybs[8465]=['41 Aqr',5.8271456,-0.3659404,5.32];
ybs[8466]=['ε Oct',5.8569251,-1.4020402,5.1];
ybs[8467]=['',5.823416,0.5011798,5.89];
ybs[8468]=['',5.8164982,1.1065074,5.79];
ybs[8469]=['',5.8332918,-0.7739542,6.1];
ybs[8470]=['',5.8241768,0.6950285,4.49];
ybs[8471]=['μ1 Gru',5.8333304,-0.7197568,4.79];
ybs[8472]=['',5.8237568,0.7949629,5.53];
ybs[8473]=['μ2 Gru',5.8369466,-0.7246549,5.1];
ybs[8474]=['',5.8278496,0.7515612,5.71];
ybs[8475]=['',5.8228736,1.1042639,6.11];
ybs[8476]=['',5.8340554,0.1510958,6.21];
ybs[8477]=['',5.8373663,-0.4501287,6.15];
ybs[8478]=['',5.8174572,1.2813182,6.08];
ybs[8479]=['ε Cep',5.8285671,0.9974741,4.19];
ybs[8480]=['',5.8366757,-0.0259801,6.15];
ybs[8481]=['42 Aqr',5.8379108,-0.2220669,5.34];
ybs[8482]=['',5.8389414,-0.4019852,6.17];
ybs[8483]=['1 Lac',5.8333692,0.6607218,4.13];
ybs[8484]=['θ Aqr',5.8379666,-0.1339616,4.16];
ybs[8485]=['',5.8381759,-0.1558944,5.79];
ybs[8486]=['',5.845277,-0.934092,5.37];
ybs[8487]=['α Tuc',5.8466764,-1.0498402,2.86];
ybs[8488]=['',5.8358878,0.4871559,6.37];
ybs[8489]=['44 Aqr',5.839128,-0.0921405,5.75];
ybs[8490]=['υ Oct',5.9131594,-1.4947591,5.77];
ybs[8491]=['',5.8347204,1.0005632,5.88];
ybs[8492]=['',5.8432501,-0.002262,6.39];
ybs[8493]=['45 Aqr',5.8475646,-0.2303243,5.95];
ybs[8494]=['',5.8556601,-1.0018405,6.34];
ybs[8495]=['',5.8463322,0.6610923,6.17];
ybs[8496]=['25 Cep',5.8421126,1.0980321,5.75];
ybs[8497]=['ρ Aqr',5.852646,-0.1346078,5.37];
ybs[8498]=['30 Peg',5.8535787,0.1029421,5.37];
ybs[8499]=['',5.8555929,0.1447833,6.17];
ybs[8500]=['ν Ind',5.8743992,-1.2591834,5.29];
ybs[8501]=['47 Aqr',5.8589547,-0.3750603,5.13];
ybs[8502]=['',5.8555884,0.4720084,6.47];
ybs[8503]=['γ Aqr',5.8589098,-0.0223098,3.84];
ybs[8504]=['',5.8534604,0.8916809,6.42];
ybs[8505]=['31 Peg',5.8580979,0.2149235,5.01];
ybs[8506]=['π1 Gru',5.8644621,-0.8000338,6.62];
ybs[8507]=['32 Peg',5.856959,0.4963618,4.81];
ybs[8508]=['2 Lac',5.8552203,0.8141171,4.57];
ybs[8509]=['π2 Gru',5.8662107,-0.7996978,5.62];
ybs[8510]=['',5.8406755,1.3368558,6.66];
ybs[8511]=['',5.8802634,-1.3073498,6.04];
ybs[8512]=['',5.8765778,-1.2273487,5.78];
ybs[8513]=['',5.8589446,0.7363077,6.41];
ybs[8514]=['49 Aqr',5.8673887,-0.4302783,5.53];
ybs[8515]=['',5.8671917,-0.1236579,5.93];
ybs[8516]=['',5.8745405,-1.0068373,5.32];
ybs[8517]=['33 Peg',5.8673069,0.3657812,6.04];
ybs[8518]=['51 Aqr',5.8696868,-0.0825096,5.78];
ybs[8519]=['50 Aqr',5.8712867,-0.2342211,5.76];
ybs[8520]=['',5.8634323,1.0017085,6.16];
ybs[8521]=['',5.8680161,0.6751464,6.22];
ybs[8522]=['',5.8631249,1.0913406,6.04];
ybs[8523]=['β Lac',5.8660956,0.9134793,4.43];
ybs[8524]=['π Aqr',5.874664,0.0259571,4.66];
ybs[8525]=['δ Tuc',5.8854007,-1.1319541,4.48];
ybs[8526]=['4 Lac',5.8703821,0.8654379,4.57];
ybs[8527]=['',5.8789808,-0.411419,6.29];
ybs[8528]=['',5.8761593,0.3238328,6.26];
ybs[8529]=['53 Aqr',5.8805749,-0.2902724,6.57];
ybs[8530]=['53 Aqr',5.8805895,-0.2902918,6.35];
ybs[8531]=['',5.8077665,1.5012608,5.27];
ybs[8532]=['',5.8912452,-1.1759802,5.55];
ybs[8533]=['34 Peg',5.8804967,0.078603,5.75];
ybs[8534]=['',5.8805332,0.6554394,6.46];
ybs[8535]=['',5.8636428,1.3675114,6.76];
ybs[8536]=['35 Peg',5.8858819,0.0838773,4.79];
ybs[8537]=['ν Gru',5.8900725,-0.6810538,5.47];
ybs[8538]=['',5.8834492,0.6967334,6.14];
ybs[8539]=['',5.8809058,0.9868684,6.57];
ybs[8540]=['',5.8850575,0.5576417,5.98];
ybs[8541]=['δ1 Gru',5.8928669,-0.7572109,3.97];
ybs[8542]=['',5.8755106,1.237101,5.47];
ybs[8543]=['ζ1 Aqr',5.8901815,0.001574,4.59];
ybs[8544]=['ζ2 Aqr',5.8902106,0.0015788,4.42];
ybs[8545]=['δ2 Gru',5.8950005,-0.7616404,4.11];
ybs[8546]=['26 Cep',5.880814,1.1386927,5.46];
ybs[8547]=['36 Peg',5.8913831,0.1612581,5.58];
ybs[8548]=['',5.8946829,-0.471179,5.95];
ybs[8549]=['',5.8912672,0.4690323,5.79];
ybs[8550]=['',5.8955847,-0.2234772,6.4];
ybs[8551]=['37 Peg',5.8950865,0.079279,5.48];
ybs[8552]=['56 Aqr',5.896765,-0.2526377,6.37];
ybs[8553]=['',5.8863281,1.1204294,6.29];
ybs[8554]=['',5.8935486,0.6254594,6.56];
ybs[8555]=['ζ PsA',5.8995836,-0.4531352,6.43];
ybs[8556]=['δ Cep',5.890371,1.0214675,3.75];
ybs[8557]=['5 Lac',5.8923629,0.8345733,4.36];
ybs[8558]=['σ Aqr',5.8982639,-0.1844329,4.82];
ybs[8559]=['38 Peg',5.8949162,0.5704293,5.65];
ybs[8560]=['',5.8948393,0.8633586,6.4];
ybs[8561]=['β PsA',5.902355,-0.5626087,4.29];
ybs[8562]=['',5.922765,-1.3728732,6.15];
ybs[8563]=['ρ1 Cep',5.8767419,1.3769909,5.83];
ybs[8564]=['6 Lac',5.8966824,0.7545776,4.51];
ybs[8565]=['',5.9010348,-0.048872,6.16];
ybs[8566]=['',5.901085,-0.1124698,6.14];
ybs[8567]=['ν Tuc',5.9098316,-1.0798511,4.81];
ybs[8568]=['58 Aqr',5.9028107,-0.1884001,6.38];
ybs[8569]=['',5.9017005,0.5175559,6.35];
ybs[8570]=['α Lac',5.8999837,0.8795313,3.77];
ybs[8571]=['39 Peg',5.9063168,0.3550208,6.42];
ybs[8572]=['',5.9072115,0.2788088,6.32];
ybs[8573]=['',5.9052944,0.6962272,5.88];
ybs[8574]=['',5.9043159,0.9450718,6.35];
ybs[8575]=['60 Aqr',5.9129786,-0.0255287,5.89];
ybs[8576]=['ρ2 Cep',5.8907094,1.3776712,5.5];
ybs[8577]=['υ Aqr',5.9160574,-0.3594806,5.2];
ybs[8578]=['',5.9221403,-1.0083074,6.23];
ybs[8579]=['',5.9101733,0.9902367,5.71];
ybs[8580]=['',5.9064801,1.2221643,6.6];
ybs[8581]=['',5.9200969,-0.4167729,5.97];
ybs[8582]=['η Aqr',5.918668,-0.0001008,4.02];
ybs[8583]=['',5.907457,1.2301984,6.34];
ybs[8584]=['',5.9020089,1.33234,5.68];
ybs[8585]=['σ1 Gru',5.9242369,-0.7063492,6.28];
ybs[8586]=['',5.9245031,-0.5506849,5.82];
ybs[8587]=['σ2 Gru',5.9263788,-0.706493,5.86];
ybs[8588]=['8 Lac',5.9202954,0.6936982,5.73];
ybs[8589]=['',5.9215159,0.6228921,6.1];
ybs[8590]=['',5.9239675,0.2061042,6.4];
ybs[8591]=['',5.9201015,0.8758571,6.29];
ybs[8592]=['',5.9197663,0.9805574,6.38];
ybs[8593]=['',5.9260147,0.2214695,6.3];
ybs[8594]=['',5.9244856,0.6242081,6.3];
ybs[8595]=['κ Aqr',5.9291924,-0.0718357,5.03];
ybs[8596]=['',5.9361131,-0.9176903,6.65];
ybs[8597]=['',5.9319149,-0.1358824,6.23];
ybs[8598]=['9 Lac',5.9265415,0.901591,4.63];
ybs[8599]=['',5.9338457,-0.4997824,6.47];
ybs[8600]=['31 Cep',5.9179041,1.2872642,5.08];
ybs[8601]=['',5.9344135,-0.5754178,5.66];
ybs[8602]=['',5.9307423,0.7905523,6.4];
ybs[8603]=['40 Peg',5.9337701,0.3426883,5.82];
ybs[8604]=['',5.9381467,-0.4924053,6.31];
ybs[8605]=['',5.9435795,-1.0002391,5.97];
ybs[8606]=['',5.9318507,0.9932344,5.21];
ybs[8607]=['10 Lac',5.9351234,0.6835182,4.88];
ybs[8608]=['',5.940967,-0.5331325,5.87];
ybs[8609]=['41 Peg',5.9377268,0.3454643,6.21];
ybs[8610]=['',5.9239504,1.3174387,5.79];
ybs[8611]=['',5.9365078,0.658081,6.03];
ybs[8612]=['30 Cep',5.9315785,1.111718,5.19];
ybs[8613]=['ε PsA',5.9421512,-0.4700331,4.17];
ybs[8614]=['',5.9424611,-0.0600646,6.31];
ybs[8615]=['β Oct',5.9696298,-1.4183937,4.15];
ybs[8616]=['',5.9425769,0.2559032,5.71];
ybs[8617]=['11 Lac',5.9404849,0.774735,4.46];
ybs[8618]=['',5.9392944,0.9417573,5.93];
ybs[8619]=['ζ Peg',5.9451711,0.1910127,3.4];
ybs[8620]=['',5.9510778,-0.8220066,5.98];
ybs[8621]=['β Gru',5.9513018,-0.8163196,2.1];
ybs[8622]=['19 PsA',5.9496481,-0.510471,6.17];
ybs[8623]=['',5.9451551,0.542425,6.34];
ybs[8624]=['',5.9514552,-0.7702961,6.07];
ybs[8625]=['12 Lac',5.9447794,0.7040375,5.25];
ybs[8626]=['ο Peg',5.9462014,0.5134825,4.79];
ybs[8627]=['',5.9472849,0.2553295,5.9];
ybs[8628]=['',5.9452974,0.7271441,5.94];
ybs[8629]=['ρ Gru',5.954797,-0.7208428,4.85];
ybs[8630]=['',5.9523712,-0.1430918,6.45];
ybs[8631]=['',5.9587655,-1.0539364,6.3];
ybs[8632]=['67 Aqr',5.9531401,-0.1195487,6.41];
ybs[8633]=['',5.9482175,0.9428593,6.12];
ybs[8634]=['66 Aqr',5.9548153,-0.3266745,4.69];
ybs[8635]=['η Peg',5.9516252,0.5294365,2.94];
ybs[8636]=['',5.9511577,0.6617564,6.43];
ybs[8637]=['',5.9516054,0.8252215,6.39];
ybs[8638]=['',5.954991,0.1929005,6.51];
ybs[8639]=['',5.956202,0.6907809,5.95];
ybs[8640]=['η Gru',5.964391,-0.931774,4.85];
ybs[8641]=['13 Lac',5.9561755,0.7318591,5.08];
ybs[8642]=['',5.9643995,-0.8104251,5.51];
ybs[8643]=['',5.9664421,-0.8528595,6.62];
ybs[8644]=['',5.9679339,-0.8651971,6.48];
ybs[8645]=['45 Peg',5.9625611,0.3399932,6.25];
ybs[8646]=['',5.9590844,0.9185774,6.55];
ybs[8647]=['',5.9689808,-0.8172629,6.56];
ybs[8648]=['ξ Oct',5.9877469,-1.3964345,5.35];
ybs[8649]=['',5.983871,-1.342792,6.73];
ybs[8650]=['ξ Peg',5.9679913,0.2144398,4.19];
ybs[8651]=['',5.965201,0.7794593,5.76];
ybs[8652]=['λ Peg',5.9671413,0.4132807,3.95];
ybs[8653]=['',5.9712974,-0.5942373,6.28];
ybs[8654]=['',5.9765689,-1.074602,6.37];
ybs[8655]=['68 Aqr',5.9721068,-0.34033,5.26];
ybs[8656]=['',5.9734007,-0.6651108,6.71];
ybs[8657]=['',5.9812154,-1.2258078,6.34];
ybs[8658]=['τ1 Aqr',5.9727436,-0.2433427,5.66];
ybs[8659]=['',5.9738684,-0.4502604,6.3];
ybs[8660]=['ε Gru',5.9770475,-0.8936594,3.49];
ybs[8661]=['70 Aqr',5.9761482,-0.1822394,6.19];
ybs[8662]=['',5.9701003,1.0227033,6.36];
ybs[8663]=['',5.9741364,0.6550327,5.9];
ybs[8664]=['τ2 Aqr',5.98093,-0.2352411,4.01];
ybs[8665]=['',5.9828961,-0.5705661,6.33];
ybs[8666]=['',5.9804297,0.1848837,6.54];
ybs[8667]=['',5.9764213,0.9517112,6.12];
ybs[8668]=['',5.9758012,1.1004711,6.06];
ybs[8669]=['μ Peg',5.982289,0.4313739,3.48];
ybs[8670]=['',5.987583,-0.6814207,5.42];
ybs[8671]=['',5.9912335,-1.0431285,6.46];
ybs[8672]=['',5.9766293,1.1987678,6.19];
ybs[8673]=['',5.98063,0.9776805,5.43];
ybs[8674]=['',5.9931983,-1.1008492,6.12];
ybs[8675]=['14 Lac',5.9835869,0.7342233,5.92];
ybs[8676]=['',5.985189,0.3360662,6.4];
ybs[8677]=['',5.9825409,0.8864736,6.21];
ybs[8678]=['21 PsA',5.9887896,-0.5135048,5.97];
ybs[8679]=['ι Cep',5.9797376,1.1574103,3.52];
ybs[8680]=['γ PsA',5.9939731,-0.5717859,4.46];
ybs[8681]=['',5.987442,1.0788121,5.6];
ybs[8682]=['σ Peg',5.9929349,0.1736631,5.16];
ybs[8683]=['λ Aqr',5.9940502,-0.1302902,3.74];
ybs[8684]=['15 Lac',5.9908602,0.7579449,4.94];
ybs[8685]=['τ1 Gru',5.9990759,-0.8461924,6.04];
ybs[8686]=['ρ Ind',6.0044774,-1.2210086,6.05];
ybs[8687]=['',5.9661141,1.4532943,4.74];
ybs[8688]=['',5.9956366,0.2959347,5.64];
ybs[8689]=['74 Aqr',5.9978592,-0.2007459,5.8];
ybs[8690]=['',5.9943648,0.8818557,6.46];
ybs[8691]=['',5.9959683,0.7030525,6.34];
ybs[8692]=['',5.9948851,1.0509641,6.01];
ybs[8693]=['',5.9979789,0.7830238,5.81];
ybs[8694]=['δ Aqr',6.0030147,-0.2741194,3.27];
ybs[8695]=['78 Aqr',6.002567,-0.1237401,6.19];
ybs[8696]=['77 Aqr',6.003492,-0.2819925,5.56];
ybs[8697]=['',6.000022,0.7067152,5.81];
ybs[8698]=['',6.0058854,-0.6330933,6.4];
ybs[8699]=['',6.0024319,0.2976939,6.12];
ybs[8700]=['1 Psc',6.0043296,0.0205899,6.11];
ybs[8701]=['',6.0052272,-0.0850456,5.72];
ybs[8702]=['ρ Peg',6.0052835,0.1558729,4.9];
ybs[8703]=['',6.0041274,0.6491218,5.91];
ybs[8704]=['',6.0084668,-0.5500917,6.1];
ybs[8705]=['δ PsA',6.0088789,-0.5659158,4.21];
ybs[8706]=['',6.0108337,-0.5489123,6.48];
ybs[8707]=['τ3 Gru',6.0128406,-0.8352082,5.7];
ybs[8708]=['',6.0071904,0.6364651,5.74];
ybs[8709]=['',6.0123671,0.208804,6.51];
ybs[8710]=['16 Lac',6.0099515,0.7281352,5.59];
ybs[8711]=['',6.0099597,0.8700257,4.95];
ybs[8712]=['',6.0144014,-0.0819376,6.31];
ybs[8713]=['α PsA',6.016263,-0.5149916,1.16];
ybs[8714]=['51 Peg',6.0149206,0.3644986,5.49];
ybs[8715]=['',6.0154516,0.0685152,6.28];
ybs[8716]=['',6.0127904,0.851711,5.43];
ybs[8717]=['',6.0204029,-0.6179783,6.13];
ybs[8718]=['',6.015607,0.688083,6.18];
ybs[8719]=['',6.0186171,-0.0397904,6.16];
ybs[8720]=['',6.0192039,-0.0225986,6.37];
ybs[8721]=['',5.9792992,1.4871836,5.9];
ybs[8722]=['',6.019929,0.1653253,6.43];
ybs[8723]=['',6.0204942,0.1301184,6.33];
ybs[8724]=['52 Peg',6.0225755,0.206725,5.75];
ybs[8725]=['',6.0247379,-0.5121945,5.51];
ybs[8726]=['',6.0245515,-0.2261108,6.07];
ybs[8727]=['2 Psc',6.0238137,0.0188215,5.43];
ybs[8728]=['',6.0268663,-0.4371781,5.65];
ybs[8729]=['',6.0218655,0.9210105,6.29];
ybs[8730]=['',6.0215473,1.0459807,6.43];
ybs[8731]=['',6.0282375,-0.4452495,6.29];
ybs[8732]=['ζ Gru',6.0307259,-0.9187125,4.12];
ybs[8733]=['',5.9957868,1.4741212,4.71];
ybs[8734]=['',6.031749,-0.8872233,5.68];
ybs[8735]=['3 Psc',6.0289479,0.0052641,6.21];
ybs[8736]=['',6.0292854,0.0545844,5.83];
ybs[8737]=['',6.0257395,0.9959018,5];
ybs[8738]=['',6.0289645,0.5445225,6.6];
ybs[8739]=['',6.0322561,-0.5015682,5.55];
ybs[8740]=['',6.0281563,0.7939635,6.5];
ybs[8741]=['',6.0324492,-0.3957526,6.28];
ybs[8742]=['81 Aqr',6.0323379,-0.1212172,6.21];
ybs[8743]=['',6.0297477,0.6776042,6.54];
ybs[8744]=['',6.0329051,-0.0802065,5.94];
ybs[8745]=['',6.0377708,-0.6336383,6.47];
ybs[8746]=['',6.031979,0.9987025,6.2];
ybs[8747]=['ο And',6.0341049,0.7407535,3.62];
ybs[8748]=['82 Aqr',6.03735,-0.1127159,6.15];
ybs[8749]=['',6.0383367,-0.3622344,5.97];
ybs[8750]=['',6.0370088,0.5567003,6.57];
ybs[8751]=['2 And',6.0370849,0.748289,5.1];
ybs[8752]=['π PsA',6.0418004,-0.604465,5.11];
ybs[8753]=['',6.0377112,0.7709981,6.39];
ybs[8754]=['',6.0487082,-1.19911,5.52];
ybs[8755]=['',6.0373659,0.9660822,6.5];
ybs[8756]=['',6.0440522,-0.7219052,5.79];
ybs[8757]=['',6.0434925,-0.0816652,6.68];
ybs[8758]=['β Psc',6.0430752,0.0686996,4.53];
ybs[8759]=['κ Gru',6.0472047,-0.9398371,5.37];
ybs[8760]=['β Peg',6.0423967,0.4921646,2.42];
ybs[8761]=['',6.0436614,0.1175109,6.41];
ybs[8762]=['',6.0401116,1.0569958,6.74];
ybs[8763]=['',6.0400234,1.0241738,6.43];
ybs[8764]=['',6.0404694,1.1750482,5.24];
ybs[8765]=['3 And',6.0438497,0.8756046,4.65];
ybs[8766]=['α Peg',6.0468359,0.2674121,2.49];
ybs[8767]=['83 Aqr',6.048789,-0.132248,5.43];
ybs[8768]=['',6.0490914,-0.2960568,6.14];
ybs[8769]=['',6.0483222,0.2911105,6.44];
ybs[8770]=['',6.0492788,0.0248415,6.39];
ybs[8771]=['',6.0652822,-1.3851642,6.12];
ybs[8772]=['θ Gru',6.0566634,-0.7575425,4.28];
ybs[8773]=['',6.0535361,0.3252244,6.13];
ybs[8774]=['86 Aqr',6.0555573,-0.4123605,4.47];
ybs[8775]=['υ Gru',6.0566538,-0.6767628,5.61];
ybs[8776]=['',6.0579816,-0.8637646,6.33];
ybs[8777]=['',6.05452,0.3495432,6.3];
ybs[8778]=['',6.0583805,-0.8826091,5.83];
ybs[8779]=['',6.0652464,-1.2822866,6.15];
ybs[8780]=['55 Peg',6.05668,0.1662604,4.52];
ybs[8781]=['56 Peg',6.0570058,0.4465411,4.76];
ybs[8782]=['1 Cas',6.0542463,1.0391034,4.85];
ybs[8783]=['',6.0584477,0.5749544,6.02];
ybs[8784]=['',6.0586422,0.3708964,5.99];
ybs[8785]=['',6.0575629,0.8060743,6.66];
ybs[8786]=['',6.0568468,0.9238547,6.11];
ybs[8787]=['',6.0628812,-0.5010247,5.6];
ybs[8788]=['',6.056681,1.0444763,6.4];
ybs[8789]=['4 And',6.0591024,0.8116457,5.33];
ybs[8790]=['5 And',6.059494,0.8624107,5.7];
ybs[8791]=['',6.0615375,0.7797848,6.56];
ybs[8792]=['5 Psc',6.0640574,0.0391747,5.4];
ybs[8793]=['',6.0592598,1.1126472,6.26];
ybs[8794]=['',6.071715,-1.1648422,6.47];
ybs[8795]=['',6.0820661,-1.410149,6.41];
ybs[8796]=['',6.0599286,1.1229305,6.21];
ybs[8797]=['88 Aqr',6.0675874,-0.3674902,3.66];
ybs[8798]=['',6.0689465,-0.4881985,5.87];
ybs[8799]=['',6.0700451,-0.7460171,5.81];
ybs[8800]=['57 Peg',6.0676853,0.1534858,5.12];
ybs[8801]=['',6.0691865,-0.2512167,6.42];
ybs[8802]=['89 Aqr',6.0696337,-0.3899168,4.69];
ybs[8803]=['',6.0709208,-0.7064171,5.83];
ybs[8804]=['π Cep',6.0588552,1.3177961,4.41];
ybs[8805]=['ι Gru',6.0718445,-0.7876618,3.9];
ybs[8806]=['58 Peg',6.0698584,0.173466,5.39];
ybs[8807]=['2 Cas',6.0679097,1.0375972,5.7];
ybs[8808]=['',6.0734647,-0.5132663,6.51];
ybs[8809]=['',6.0727848,0.309123,5.71];
ybs[8810]=['6 And',6.0713828,0.7620305,5.94];
ybs[8811]=['59 Peg',6.0773382,0.1542367,5.16];
ybs[8812]=['60 Peg',6.0775496,0.4706165,6.17];
ybs[8813]=['',6.0844983,-0.8639661,6.8];
ybs[8814]=['',6.0885524,-1.0922729,6.12];
ybs[8815]=['7 And',6.0804566,0.8643496,4.52];
ybs[8816]=['',6.0829704,0.5159004,6.35];
ybs[8817]=['',6.0835032,0.9998223,5.56];
ybs[8818]=['',6.0847625,0.1951678,5.82];
ybs[8819]=['φ Aqr',6.0887343,-0.1035243,4.22];
ybs[8820]=['',6.0918925,-0.7153773,5.77];
ybs[8821]=['',6.0902801,-0.1845021,6.12];
ybs[8822]=['',6.08782,0.8854954,6.31];
ybs[8823]=['',6.0886243,0.5216623,6.41];
ybs[8824]=['',6.0897526,0.4227269,6.36];
ybs[8825]=['',6.0941621,-0.0589726,5.55];
ybs[8826]=['ψ1 Aqr',6.0955972,-0.1565601,4.21];
ybs[8827]=['61 Peg',6.0947976,0.4950679,6.49];
ybs[8828]=['',6.1009244,-1.0800701,5.66];
ybs[8829]=['',6.0885443,1.2976263,5.84];
ybs[8830]=['',6.0956721,0.434389,6.6];
ybs[8831]=['',6.099274,-0.7744295,5.92];
ybs[8832]=['',6.0999658,-0.7169256,6.47];
ybs[8833]=['γ Tuc',6.1028594,-1.0143528,3.99];
ybs[8834]=['',6.1116316,-1.3850044,6.33];
ybs[8835]=['χ Aqr',6.0997607,-0.1328027,5.06];
ybs[8836]=['',6.0932399,1.2392808,5.56];
ybs[8837]=['γ Psc',6.1010681,0.0593392,3.69];
ybs[8838]=['',6.0985634,0.9308055,5.54];
ybs[8839]=['',6.0972211,1.0835116,6.53];
ybs[8840]=['',6.1070822,-1.1755373,6.13];
ybs[8841]=['',6.1033575,-0.2023769,6.34];
ybs[8842]=['',6.101192,0.7903172,6.43];
ybs[8843]=['ψ2 Aqr',6.1043725,-0.15821,4.39];
ybs[8844]=['φ Gru',6.1057762,-0.7104656,5.53];
ybs[8845]=['8 And',6.1031846,0.8575325,4.85];
ybs[8846]=['',6.1040662,0.7959857,6.48];
ybs[8847]=['τ Oct',6.1548305,-1.5233104,5.49];
ybs[8848]=['γ Scl',6.1085606,-0.5657331,4.41];
ybs[8849]=['9 And',6.1061019,0.7311427,6.02];
ybs[8850]=['ψ3 Aqr',6.1089917,-0.1656841,4.98];
ybs[8851]=['94 Aqr',6.1096715,-0.2328451,5.08];
ybs[8852]=['',6.1002774,1.316272,6.38];
ybs[8853]=['96 Aqr',6.1108747,-0.0873811,5.55];
ybs[8854]=['',6.1109679,-0.3134158,5.93];
ybs[8855]=['',6.1089023,0.7898497,6.5];
ybs[8856]=['',6.1124785,-0.5862588,6.37];
ybs[8857]=['ο Cep',6.1065344,1.1908287,4.75];
ybs[8858]=['',6.1108369,0.6093155,6.32];
ybs[8859]=['11 And',6.1108516,0.8507286,5.44];
ybs[8860]=['',6.1117152,0.8464625,6.32];
ybs[8861]=['10 And',6.1125827,0.7364586,5.79];
ybs[8862]=['',6.1175173,-0.8759574,6.05];
ybs[8863]=['7 Psc',6.1149248,0.0959818,5.05];
ybs[8864]=['',6.1164708,-0.1010557,6.17];
ybs[8865]=['τ Peg',6.1160831,0.4164052,4.6];
ybs[8866]=['',6.1138293,1.083639,6.45];
ybs[8867]=['63 Peg',6.1168612,0.5329014,5.59];
ybs[8868]=['',6.1191269,-0.468946,5.64];
ybs[8869]=['',6.1163187,0.7720355,6.13];
ybs[8870]=['12 And',6.1170607,0.668465,5.77];
ybs[8871]=['',6.1153045,1.0878817,6.39];
ybs[8872]=['64 Peg',6.1216043,0.557294,5.32];
ybs[8873]=['',6.1218853,0.4664739,6.62];
ybs[8874]=['',6.1268977,-1.0461092,6.09];
ybs[8875]=['97 Aqr',6.1251271,-0.2604208,5.2];
ybs[8876]=['65 Peg',6.125005,0.36559,6.29];
ybs[8877]=['98 Aqr',6.1265403,-0.3487582,3.97];
ybs[8878]=['66 Peg',6.1268125,0.2169807,5.08];
ybs[8879]=['',6.1239648,1.0515915,5.56];
ybs[8880]=['',6.1309388,-0.9370685,6.15];
ybs[8881]=['',6.1301546,-0.7505996,6.1];
ybs[8882]=['',6.1288683,0.0071492,6.31];
ybs[8883]=['',6.1322728,-0.9036111,5.75];
ybs[8884]=['',6.1297976,0.5698436,6.69];
ybs[8885]=['',6.131589,-0.324094,6.19];
ybs[8886]=['',6.1371606,-0.990139,5.59];
ybs[8887]=['',6.1331819,0.7196182,6.72];
ybs[8888]=['67 Peg',6.1344082,0.5672902,5.57];
ybs[8889]=['4 Cas',6.1339708,1.0891047,4.98];
ybs[8890]=['υ Peg',6.1368008,0.4105458,4.4];
ybs[8891]=['99 Aqr',6.1399579,-0.3582029,4.39];
ybs[8892]=['ο Gru',6.1426925,-0.918099,5.52];
ybs[8893]=['',6.1452019,-1.1599912,6.45];
ybs[8894]=['',6.1455696,-1.0185321,5.63];
ybs[8895]=['',6.145016,-0.8733403,6.2];
ybs[8896]=['κ Psc',6.1437061,0.0239817,4.94];
ybs[8897]=['9 Psc',6.1450739,0.0216598,6.25];
ybs[8898]=['13 And',6.144268,0.751023,5.75];
ybs[8899]=['',6.1486215,-0.6182981,6.32];
ybs[8900]=['69 Peg',6.1468056,0.4413199,5.98];
ybs[8901]=['θ Psc',6.1481966,0.113402,4.28];
ybs[8902]=['',6.1488026,-0.1977658,6.37];
ybs[8903]=['',6.1443867,1.2300772,5.6];
ybs[8904]=['',6.153346,-1.099421,5.68];
ybs[8905]=['',6.1530596,-0.774562,6.43];
ybs[8906]=['',6.1528188,-0.1596534,6.18];
ybs[8907]=['',6.1530198,0.4043304,6.35];
ybs[8908]=['70 Peg',6.153344,0.2247846,4.55];
ybs[8909]=['',6.1550861,-0.0770406,6.25];
ybs[8910]=['',6.1573228,0.8596056,6.17];
ybs[8911]=['',6.1568002,1.0239428,4.91];
ybs[8912]=['',6.1597708,0.6768508,6.05];
ybs[8913]=['',6.1615652,-0.1076791,6.39];
ybs[8914]=['',6.1636729,-0.7805951,6.02];
ybs[8915]=['14 And',6.1625127,0.6868775,5.22];
ybs[8916]=['',6.1637657,-0.0692619,6.49];
ybs[8917]=['100 Aqr',6.164617,-0.3708934,6.29];
ybs[8918]=['',6.1644546,0.4978103,6.41];
ybs[8919]=['13 Psc',6.1656503,-0.0168773,6.38];
ybs[8920]=['',6.1726535,-1.3485522,5.81];
ybs[8921]=['',6.1674369,0.6121108,6.65];
ybs[8922]=['β Scl',6.1702465,-0.6579792,4.37];
ybs[8923]=['',6.1375073,1.524186,5.58];
ybs[8924]=['101 Aqr',6.1714815,-0.3629504,4.71];
ybs[8925]=['71 Peg',6.172129,0.3947554,5.32];
ybs[8926]=['',6.1730536,0.7884873,6.24];
ybs[8927]=['',6.1741311,0.3658173,6.06];
ybs[8928]=['72 Peg',6.1742018,0.5488054,4.98];
ybs[8929]=['14 Psc',6.1752063,-0.0196966,5.87];
ybs[8930]=['',6.1802979,-1.1269663,7.4];
ybs[8931]=['',6.1782003,-0.2640129,5.96];
ybs[8932]=['15 And',6.1770824,0.7043343,5.59];
ybs[8933]=['73 Peg',6.1771752,0.5867136,5.63];
ybs[8934]=['ι Phe',6.1794483,-0.7416948,4.71];
ybs[8935]=['',6.1777671,0.6657191,6.18];
ybs[8936]=['',6.181273,-0.1282014,6.39];
ybs[8937]=['',6.1781623,1.2524698,5.84];
ybs[8938]=['',6.1828774,0.4307504,6.45];
ybs[8939]=['16 Psc',6.1849602,0.0387692,5.68];
ybs[8940]=['',6.1853562,0.5763646,6.35];
ybs[8941]=['',6.1881558,-0.5541718,6.52];
ybs[8942]=['',6.1945668,-1.3395542,6];
ybs[8943]=['',6.1905626,-0.2258651,5.65];
ybs[8944]=['',6.1915452,-0.7919141,4.74];
ybs[8945]=['74 Peg',6.1904711,0.295741,6.26];
ybs[8946]=['λ And',6.189892,0.8129256,3.82];
ybs[8947]=['',6.1897679,0.7775147,5.8];
ybs[8948]=['75 Peg',6.1917022,0.3232302,5.53];
ybs[8949]=['',6.191694,0.8084172,6.58];
ybs[8950]=['ι And',6.1924134,0.7572501,4.29];
ybs[8951]=['θ Phe',6.1985918,-0.8119016,6.09];
ybs[8952]=['18 And',6.1967424,0.8829777,5.3];
ybs[8953]=['ω1 Aqr',6.1998339,-0.2461335,5];
ybs[8954]=['ι Psc',6.2004911,0.1002805,4.13];
ybs[8955]=['',6.2003401,0.1709809,5.97];
ybs[8956]=['',6.1964109,1.3161878,5.95];
ybs[8957]=['',6.197255,1.2936731,5.98];
ybs[8958]=['',6.2007939,0.6592417,6.53];
ybs[8959]=['γ Cep',6.1970334,1.3570237,3.21];
ybs[8960]=['μ Scl',6.2036119,-0.5576983,5.31];
ybs[8961]=['κ And',6.2023387,0.7758542,4.14];
ybs[8962]=['',6.2035516,0.6429815,6.23];
ybs[8963]=['',6.2056739,-0.419594,6.6];
ybs[8964]=['',6.2057717,-0.2017818,5.89];
ybs[8965]=['103 Aqr',6.207652,-0.3125517,5.34];
ybs[8966]=['',6.2068538,0.8662339,6.26];
ybs[8967]=['104 Aqr',6.2084724,-0.3088717,4.82];
ybs[8968]=['',6.2091931,0.1286291,5.89];
ybs[8969]=['λ Psc',6.2096521,0.03315,4.5];
ybs[8970]=['',6.2088042,1.0014585,6.24];
ybs[8971]=['',6.2103728,0.7873408,6.57];
ybs[8972]=['',6.211525,-0.2675312,5.28];
ybs[8973]=['ω2 Aqr',6.2126418,-0.2517746,4.49];
ybs[8974]=['',6.2106378,1.1280922,6.56];
ybs[8975]=['',6.211457,1.0785928,6.4];
ybs[8976]=['77 Peg',6.215419,0.1824008,5.06];
ybs[8977]=['',6.2174554,-0.2646797,6.36];
ybs[8978]=['',6.2184225,-0.7847682,6.09];
ybs[8979]=['',6.2204113,-1.2282028,6.07];
ybs[8980]=['',6.2218118,-1.3730844,5.75];
ybs[8981]=['',6.2193414,-1.1219851,5.72];
ybs[8982]=['78 Peg',6.2180644,0.5145421,4.93];
ybs[8983]=['106 Aqr',6.2191059,-0.3169084,5.24];
ybs[8984]=['',6.2203489,-0.4560013,6.17];
ybs[8985]=['',6.2214957,0.9759737,6.51];
ybs[8986]=['',6.2270971,-0.6992314,6.31];
ybs[8987]=['107 Aqr',6.2270135,-0.3239081,5.29];
ybs[8988]=['ψ And',6.2269308,0.8122722,4.95];
ybs[8989]=['19 Psc',6.2286085,0.0629395,5.04];
ybs[8990]=['',6.2293017,1.1676555,5.95];
ybs[8991]=['σ Phe',6.2325667,-0.8745345,5.18];
ybs[8992]=['',6.2332405,-1.1916172,6.89];
ybs[8993]=['τ Cas',6.2313386,1.0257555,4.87];
ybs[8994]=['',6.2324498,-0.2057972,5.73];
ybs[8995]=['',6.2312304,1.0048019,5.51];
ybs[8996]=['',6.2335687,0.8194676,6.07];
ybs[8997]=['20 Psc',6.2353839,-0.0461137,5.49];
ybs[8998]=['',6.2349963,1.1855408,5.04];
ybs[8999]=['',6.238008,-0.109275,6.07];
ybs[9000]=['',6.2392148,0.0407313,6.46];
ybs[9001]=['δ Scl',6.239728,-0.4888791,4.57];
ybs[9002]=['',6.2382588,1.1343934,6.41];
ybs[9003]=['6 Cas',6.2390998,1.0879337,5.43];
ybs[9004]=['',6.2393861,1.048916,6.34];
ybs[9005]=['',6.2407136,1.0311865,6.33];
ybs[9006]=['',6.2423205,-0.2747415,6.24];
ybs[9007]=['21 Psc',6.2419949,0.0208688,5.77];
ybs[9008]=['',6.2434191,-1.094668,6.59];
ybs[9009]=['',6.2429126,0.6378282,5.9];
ybs[9010]=['79 Peg',6.2428128,0.5054838,5.97];
ybs[9011]=['',6.2436453,-0.4400288,6.42];
ybs[9012]=['',6.2454438,-0.1719946,5.94];
ybs[9013]=['',6.2458771,0.9030555,6.44];
ybs[9014]=['',6.2468025,-0.2492738,5.72];
ybs[9015]=['80 Peg',6.250253,0.1646361,5.79];
ybs[9016]=['108 Aqr',6.2502977,-0.3279346,5.18];
ybs[9017]=['γ1 Oct',6.2540224,-1.4294117,5.11];
ybs[9018]=['22 Psc',6.2529305,0.053231,5.55];
ybs[9019]=['',6.2525988,1.3564538,6.55];
ybs[9020]=['',6.2547628,0.3803155,6.11];
ybs[9021]=['φ Peg',6.2551955,0.3357999,5.08];
ybs[9022]=['',6.2552845,-0.2466407,5.87];
ybs[9023]=['',6.2546614,1.3205922,6.39];
ybs[9024]=['82 Peg',6.2557727,0.1931581,5.3];
ybs[9025]=['',6.2567688,-0.1549333,5.75];
ybs[9026]=['24 Psc',6.2571334,-0.0529866,5.93];
ybs[9027]=['25 Psc',6.257797,0.0385754,6.28];
ybs[9028]=['',6.2589877,-0.4207903,6.24];
ybs[9029]=['',6.2633933,-0.4698872,6.35];
ybs[9030]=['ρ Cas',6.2634133,1.0056432,4.54];
ybs[9031]=['',6.2646575,-0.7012791,6.03];
ybs[9032]=['',6.2652022,0.003994,5.61];
ybs[9033]=['26 Psc',6.2667391,0.1255029,6.21];
ybs[9034]=['',6.2674076,-0.5550495,6.1];
ybs[9035]=['',6.2674076,-0.554395,6.83];
ybs[9036]=['',6.2678313,0.455089,6.54];
ybs[9037]=['',6.2685718,1.0041211,6];
ybs[9038]=['',6.2685817,0.828604,6];
ybs[9039]=['',6.2727268,-0.4296571,6.31];
ybs[9040]=['',6.2735476,0.397372,6.15];
ybs[9041]=['',6.2723223,1.4540477,6.59];
ybs[9042]=['',6.2751451,0.7466174,5.97];
ybs[9043]=['',6.2755175,-0.4625807,6.26];
ybs[9044]=['',6.2754906,0.9743392,5.55];
ybs[9045]=['',6.2763882,-1.0967073,5.97];
ybs[9046]=['γ2 Oct',6.2774057,-1.432048,5.73];
ybs[9047]=['η Tuc',6.2774984,-1.1201286,5];
ybs[9048]=['',6.2773096,1.0496986,6.47];
ybs[9049]=['ψ Peg',6.278205,0.440889,4.66];
ybs[9050]=['1 Cet',6.2808113,-0.274502,6.26];
ybs[9051]=['',6.2810584,0.8989896,4.8];
ybs[9052]=['27 Psc',6.2822051,-0.0599767,4.86];
ybs[9053]=['',6.282841,0.5672558,6.52];
ybs[9054]=['π Phe',0.0001456,-0.9184994,5.13];
ybs[9055]=['',6.2826418,0.8121497,6.54];
ybs[9056]=['σ Cas',0.0004752,0.9751974,4.88];
ybs[9057]=['ω Psc',0.0018045,0.1218769,4.01];
ybs[9058]=['',0.0024736,-0.5125212,5.62];
ybs[9059]=['',0.0025685,0.5906917,6.58];
ybs[9060]=['',0.0025685,0.5906917,6.58];
ybs[9061]=['ε Tuc',0.0044351,-1.1424493,4.5];
ybs[9062]=['',0.0061967,-0.7709269,6.29];
ybs[9063]=['',0.0065506,0.4719026,6.46];
ybs[9064]=['',0.0070718,1.0416024,6.19];
ybs[9065]=['',0.0079979,0.7919088,6.38];
ybs[9066]=['τ Phe',0.0094818,-0.8498062,5.71];
ybs[9067]=['',0.0106052,-0.8764613,5.53];
ybs[9068]=['',0.0105954,0.8744336,6.22];
ybs[9069]=['θ Oct',0.0116832,-1.3429635,4.78];
ybs[9070]=['',0.011891,1.0706329,5.55];
ybs[9071]=['',0.0123756,0.7415365,6.25];
ybs[9072]=['29 Psc',0.0127628,-0.0507508,5.1];
ybs[9073]=['85 Peg',0.0142892,0.474758,5.75];
ybs[9074]=['30 Psc',0.0133578,-0.102878,4.41];
ybs[9075]=['',0.0140592,-0.2540575,7.1];
ybs[9076]=['ζ Scl',0.0149669,-0.5166278,5.01];
ybs[9077]=['31 Psc',0.0152989,0.1584171,6.32];
ybs[9078]=['32 Psc',0.0156988,0.1501898,5.63];
ybs[9079]=['',0.0162252,1.1557321,5.86];
ybs[9080]=['',0.0177119,-0.3477818,6.25];
ybs[9081]=['',0.018443,-0.4193258,6.44];
ybs[9082]=['',0.0198405,1.1128503,6.24];
ybs[9083]=['2 Cet',0.0211147,-0.3004835,4.55];
ybs[9084]=['',0.0217663,1.1664366,6.29];
ybs[9085]=['9 Cas',0.0233337,1.0892154,5.88];
ybs[9086]=['',0.0236808,-0.2863949,5.78];
ybs[9087]=['',0.0237128,-0.508745,6.4];
ybs[9088]=['3 Cet',0.0244418,-0.1813358,4.94];
ybs[9089]=['',0.02543,1.174368,5.67];
ybs[9090]=['',0.0249729,0.7367365,6.01];
ybs[9091]=['',0.0243451,-1.2702176,7.31];
ybs[9092]=['',0.0262082,0.6070148,6.12];
ybs[9093]=['',0.0251235,-1.2447213,5.59];
ybs[9094]=['',0.0263589,0.4671994,6.25];
ybs[9095]=['',0.0271702,1.0722226,5.8];
    return ybs;
  };  
  var stars = build_stars();
  
  var build_constellation_lines = function(){
    var constell = {};
    constell.Psc = [[348,379,356,347],[356,433,459,506,591,545,485,430,357,290,220,9057,8954,8901,8837,8896,8969,8954],[8837,8758]];    
    constell.Ari = [[820,797,613,549,541],[968,883,613],[947,883],[834,797]];
    constell.Tau = [[1786,1492,1404,1368,1341,1406,1452,1904],[1492,1387,1382,1251,1160],[1341,1234,1169,1034,1026],[1452,1468,1453,1315,1246,1096],[2028,2078,1904,1839,1940,2004]];
    constell.Gem = [[2882,2689,2466,2337],[2466,2280,2210,2128],[2982,2769,2642,2415],[2769,2755,2477],[2977,2897,2813,2689,2532]];
    constell.Cnc = [[3467,3441,3453,3563],[3241,3453,3349,3358,3441],[3349,3200]];
    constell.Leo = [[3721,3895,4021,4047,3965,3972,4349,4524,4347,4290,4047],[4517,4347,4349,4389,4376,4408,4461,4422,4358,4289],[3895,3863,3763],[3842,3972,4123],[4249,4347,4517]];
    constell.Vir = [[5046,5328,5477,5501,5254,5095,4900,4598,4507,4530,4679,4815,4953,5046],[4922,4900],[4892,5046,5058,5009]];
    constell.Lib = [[5521,5675,5777,5593,5521],[5560,5675,5767],[5802,5784,5593]];
    constell.Sco = [[6017,5974,5943,5934,5918],[5943,6074,6124,6155,6231,6237,6252,6369,6541,6603,6515]];
    constell.Sgr = [[6734,6847,6901,7027,7109,7221,7181,6867,6734],[7181,7027,6847,6867],[6800,6901,7138,7204,7251,7327,7329],[6867,6820]];
    constell.Cap = [[7733,7762,7808,7922,7966,8066,8190,8246,8274,8307,8264,8153,8061,7875,7733]];
    constell.Aqr = [[8524,8543,8503,8399,8387,8250,8079,7936,7937],[8218,8250],[8399,8484,8558,8683,8826,8694,8664,8683],[8582,8543],[8826,8953],[8826,8967],[8826,8877],[8826,8797]];
    constell.UMi = [[420,6777,6311,5893,6106,5725,5553,5893]];
    constell.Cas = [[538,399,260,165,20]];
    constell.Cep = [[8959,8679,8450,8148,8224,8679],[8224,8959]];
    constell.Dra = [[4424,4777,5281,5734,5976,6122,6385,6915,7339,7568,7671,7113,6911,6676,6542,6524,6693,6676],[6676,7358,7297,7168,6908,6385]];
    constell.Cam = [[1563,1598,1537,1143,1150,1031],[1537,2519,2203],[1598,1150]];
    constell.UMa = [[5181,5044,4895,4650,4544,4285,4291],[4544,4508,4325,4238,4059,4023],[4325,4367,4364],[3315,3615,3747,3878,3884,3765,3560],[3765,3585]];
    constell.Peg = [[8760,14,38,8766,8780,8702,8536,8435,8293,8298,8159,8211,8300,8439,8635,8760],[8760,8669,8652,8415],[8908,8766,8650,8619,8505],[8760,8766]];
    constell.Equ = [[8117,8164,8109,8083]];                    
    constell.Del = [[7934,7914,7868,7892,7934],[7868,7838]];
    constell.Sge = [[7621,7523,7475],[7523,7466]];
    constell.Vul = [[7730,7639,7578,7392,7293]];
    constell.Aql = [[7164,7222,7512,7543,7588,7696,7870],[7696,7364,7222],[7512,7364,7223]];
    constell.Sct = [[7137,7051,6961,7008,7137],[7008,6918]];
    constell.Lyr = [[6989,7041,7044,7127,7166,7094,7044]];
    constell.Lyr = [[6989,7039,7044,7127,7166,7094,7044,6989]];
    constell.Cyg = [[7910,7782,7601,7404],[7515,7782,7935]];
    constell.Her = [[5904,6013,6082,6158,6210,6407,6313,6202,6210],[6576,6407,6683],[6313,6399,6514,6611,6691,6767],[6202,6138,6395,6775,6883,7049,7121],[6085,6138,6399]];
    constell.Oph = [[6065,6046,6139,6289,6544,6591,6617,6555,6367,6434,6442,6475],[6367,6165,6065,6046,6139]];
    constell.Ser = [[5962,5832,5869,5923,5857,5832],[5857,5778,5844,5882,5871,6046],[7129,6857,6686,6569,6549,6435]];
    constell.CrB = [[5768,5737,5783,5839,5879,5937,5961,6093,6008,5824,5768]];
    constell.Boo = [[5330,5419,5425,5592,5723,5671,5496,5330],[5764,5592,5341,5340,5394],[5467,5465,5492,5534,5330,5225,5175,5190],[5496,5606,5590],[5671,5425]];
    constell.Com = [[4973,4727,4697,4687,4910,4958,4973]];
    constell.CVn = [[5117,5007,4905,4775,4836]];
    constell.Crv = [[4613,4620,4652,4747,4776,4620]];
    constell.Crt = [[4395,4333,4277,4372],[4458,4392,4372,4395,4504]];
    constell.And = [[460,386,265,151,14,162,333,473,599],[162,160,211,267],[14,67,62,8950],[8747,8950,8961,8946,8845,8815,8765]];
    constell.Lac = [[8523,8598,8570,8526,8557,8508,8617,8564,8607,8470,8483]];
    constell.Cet = [[73,184,330,73],[330,398,535,505,330],[535,677,775,800,750,714,809,892,907,800]];
    constell.Eri = [[1674,1661,1555,1515,1458,1293,1226,1157,1131,1079,980,870,807,814,846,915,999,1083,1168,1208,1235,1459,1388,1342,1190,1004,893,790,717,670,562,468]];
    constell.Ori = [[2055,1785,1846,1708,1998,1942,2055],[1942,1897,1846],[2041,2153,2118,2055],[2118,2193,2129],[1785,1596,1562,1547,1538,1539,1565,1575,1551],[1575,1633,1671],[2055,1873,1785]];
    constell.CMi = [[3137,2935,2837,2856]];
    constell.Mon = [[2450,2379,2292],[2379,2498,2706,2962,3180],[2706,2350,2221]];
    constell.CMa = [[2566,2649,2588,2585,2566],[2588,2484,2288],[2484,2582,2645,2685,2819],[2685,2638,2610,2530,2355],[2276,2610],[2645,2572,2381]];
    constell.Lep = [[2149,2079,1992,1859,2029,1977,1824,1649,1697],[1859,1697,1751,1700,1691]];
    constell.Hya = [[3402,3474,3538,3446,3410,3402],[3538,3656,3835,3738,3893,3960,3984,4084,4222,4440,4542,5010,5277,5371,5516],[3738,3696,3476,3433]];
    constell.Sex = [[3971,4109]];
    constell.Aur = [[2082,2071,1703],[1703,2082,2089,1786,1572,1636,1703],[2213,2089]];
    constell.Per = [[1126,1198,1223,1215,1130,1117,1013,911,830,850,933,937,932,917],[1117,1268,1298,1256],[1449,1301,1215],[933,795]];
    constell.Tri = [[540,618,660,540]];
    constell.LMi = [[4237,4090,3964,3790],[4237,4156,4080,3964]];
    constell.Lyn = [[3695,3680,3603,3570,3267,2810,2938,2552,2464,2232]];    
    constell.PsA = [[8613,8713,8705,8680,8561,8416,8290]];
    constell.Mic = [[7951,8025,8121,8137,7951]];
    constell.Scl = [[276,102,9001,8848,8922]];
    constell.Gru = [[8338,8396,8410,8621,8805,8772,8629,8541,8410],[8732,8660,8621],[8640,8660]];
    constell.Ind = [[7855,7906,7972,8372,8126,7906]];
    constell.Tuc = [[8487,8833,123,76,9061,8525,8487]];
    constell.Phe = [[96,24,318,425,436,334,318,96]];
    constell.For = [[1129,959,837,745,608]];
    constell.Cae = [[1647,1498,1497,1438]];
    constell.Col = [[1856,1950,2034,1856],[2034,2114,2100,2250,2290,2100,2034]];
    constell.Pic = [[2542,2036,2014,2206,2542]];
    constell.Dor = [[1333,1460,1669,1916,1460],[1916,2009]];
    constell.Ret = [[1350,1331,1170,1242,1350]];
    constell.Hor = [[1321,798,930,905]];
    constell.Pyx = [[3708,3619,3547,3509,3460,3430],[3723,3547]];
    constell.Ant = [[4094,3937,3755,3861,4094,4263]];
    constell.Pup = [[2545,2445,2765,3037,3094,3177,3157,2870,2545]];
    constell.Vel = [[3199,3477,3724,3930,4206,3776,3625,3199]];
    constell.Cru = [[4843,4646],[4720,4753]];
    constell.Cen = [[5122,4809,5221],[5449,5122,5257],[4809,4733,4611,4380],[5182,5158,5211,5200,5182,5018],[4733,4628,4457],[5122,5221,5250,5238,5275,5357,5278,5180,5183,5221],[5018,5180,5238,5430,5566]];
    constell.Lup = [[5385,5459,5344],[5415,5443,5459,5518,5561,5695,5702,5685,5616],[5459,5595,5636,5639],[5595,5616,5673,5688,5698],[5695,5810,5873,5977],[5702,5766,5787],[5702,5938]];
    constell.Cir = [[5694,5453,5660]];
    constell.TrA = [[6207,5887,5661,6207]];
    constell.Nor = [[5970,6105,6062,5952,5970]];
    constell.Ara = [[6219,6275,6285,6499,6450,6451,6489,6219]];
    constell.Aps = [[5460,6010,6153,6092,6010]];
    constell.Mus = [[4510,4661,4788,4834,4913,4763,4788]];
    constell.Cha = [[3332,3310],[3310,4164,4573,4664,4224,4164]];
    constell.Oct = [[8615,8240,5329,8615]];
    constell.Pav = [[7776,8167,7899,7651,7776],[7651,7576,6970,7095,7651],[7095,7062,6843,6733,7062],[6733,6570]];
    constell.Tuc = [[8487,8525,9061,76,123,8833,8487]];
    constell.Hyi = [[587,95,1203,701,587],[566,701,802,833]];
    constell.CrA = [[7140,7213,7241,7246,7229,6939]];
    constell.Tel = [[6771,6885,6893,7659]];
    constell.Car = [[2320,3109,3299,3689,4189,4027,3675,2320]];
    constell.Vol = [[3606,3339,3215,2795,2728,3215,3606]];
    return constell;
  };
  var constellation_lines = build_constellation_lines();
  
  //END OF PRIVATE ITEMS
  
  /* Return the object that contains the items needed by the caller  */
  
  return {
  
    testing: function(){
      console.log('ephemjs is indeed visible.');
    },
    
    /* main items: */
    when: when,
    when_now: when_now,
    when_j2000: when_j2000,
    delta_t: delta_t,
    date_time_odometer: date_time_odometer,
    lmst : lmst,
    
    where : where,
    
    position : position,
    position_from_orbit : position_from_orbit,
    find_visible_messiers: find_visible_messiers,
    find_visible_caldwells: find_visible_caldwells,
    find_visible_stars: find_visible_stars,
    current_meteor_showers : current_meteor_showers,
    
    current_events : current_events,

    rise_culmination_set : rise_culmination_set,
    rise_culmination_set_daily : rise_culmination_set_daily,
    rise_culmination_set_observation_window : rise_culmination_set_observation_window,
    observation_window : observation_window,
    twilight : twilight,
    
    physical_jupiter: physical_jupiter, 
    lunar_libration: lunar_libration,
    physical_sun: physical_sun,
    
    /* in-memory databases */
    planets : planets, /*object*/ 
    minor_planets : minor_planets, /*object*/ 
    messiers: messiers, /*array*/
    caldwells: caldwells, /*array*/
    comets: comets, /*object*/
    stars: stars, /*array*/
    constellation_lines: constellation_lines, /* array of arrays */
    zodiac: zodiac,

    /* utility functions */
    find_calendar_date_from_jd : find_calendar_date_from_jd,
    elongation_between : elongation_between, 
    delta_longitude_between : delta_longitude_between, 
    convert_αδ_to_λβ : convert_αδ_to_λβ,
    convert_λβ_to_αδ : convert_λβ_to_αδ,
    convert_αδ_to_XYZ : convert_αδ_to_XYZ,
    convert_XYZ_to_αδ : convert_XYZ_to_αδ,
    convert_λβ_to_XYZ : convert_λβ_to_XYZ,
    convert_αδ_to_aA : convert_αδ_to_aA,
    convert_ra_to_zodiac_sign : convert_ra_to_zodiac_sign,
    apply_parallax_to_αδ : apply_parallax_to_αδ, 
    apply_parallax_to_λβ : apply_parallax_to_λβ,
    distance_kms : distance_kms,
    parallactic_angle : parallactic_angle,
    position_angle_between : position_angle_between,
    bright_limb_angle : bright_limb_angle,
    geomagnetic_north_pole : geomagnetic_north_pole,
    geomagnetic_latitude: geomagnetic_latitude,
    aurora_min_kp : aurora_min_kp,
    apply_precession : apply_precession,
    precession_angles : precession_angles,
    nutation : nutation,
    annual_aberration: annual_aberration,
    add_refraction_to_alt : add_refraction_to_alt,
    convert_all_angles_to_degs : convert_all_angles_to_degs,  
    convert_all_angles_to_degs_sexagesimal : convert_all_angles_to_degs_sexagesimal,
    convert_all_angles_to_rads : convert_all_angles_to_rads,  
    rads : rads,
    degs : degs,
    round: round,
    round_and_pad: round_and_pad,
    as_array: as_array,
    in360: in360,
    in2pi: in2pi
  }; 
}()); // the top level function is invoked here; its return value is stored in EPH, a global variable