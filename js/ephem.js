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
   
   PROBLEM: sadly, the Observer's Handbook no longer (as of 2023, maybe earlier) has the data in the same form! 
   It has elements only for one day of the year (Feb 25, not mid-year!), not two.
   I will retain the old algorithm, but this is not really appropriate for the current data.
   PROBLEM: the Handbook data for 2023 has bad values of L for Mercury and Venus, and likely the other planets as well.
   Maybe I should abandon the Handbook and just use the Astronomical Almanac directly.
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
  
  var planet_orbit_start = when('UT 2024-03-31'); //Observer's Handbook p 23; usually a 240d interval
  var planet_orbit_end = when('UT 2024-10-17'); //In 2024, Uranus and Neptune only use the spring date
  
  var planet_orbit_fixed_date = when('UT 2023-02-25'); //AS OF 2023, USES ONLY ONE DATE!
  
  var mercury = {
    name: 'Mercury',
    symbol: '☿',
    orbit: function(when){
      var osc_a = build_osculating_orbit(when_j2000, planet_orbit_start, 0.387099, 0.205630, 7.0036, 48.3005, 77.4951, 131.9541);
      var osc_b = build_osculating_orbit(when_j2000, planet_orbit_end,   0.387098, 0.205645, 7.0035, 48.3001, 77.4954, (360*2)+230.4244); //4.0914 deg per d; 818 deg per 200d
      var result = current_osculating_orbit(osc_a, osc_b, when);
      return result;
      //return build_osculating_orbit(when_j2000, planet_orbit_fixed_date, 0.386734, 0.205623, 7.0036, 48.3027, 77.4937, 98.0867 /*293.7 or so?*/);
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
      var osc_a = build_osculating_orbit(when_j2000, planet_orbit_start, 0.723325, 0.006753, 3.3944, 76.6119, 131.8760,        329.6470);
      var osc_b = build_osculating_orbit(when_j2000, planet_orbit_end,   0.723327, 0.006740, 3.3944, 76.6117, 131.8490, (360*1)+290.0732); //1.60212 deg per d; 320 deg per 200d
      
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
      var osc_a = build_osculating_orbit(when_j2000, planet_orbit_start, 1.523618, 0.093272, 1.8479, 49.4892, 336.1640, 316.0197);
      var osc_b = build_osculating_orbit(when_j2000, planet_orbit_end,   1.523780, 0.093374, 1.8477, 49.4894, 336.2096, (360*1)+60.8245); //0.52402 deg per d; 105deg in 200d 
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
      var osc_a = build_osculating_orbit(when_j2000, planet_orbit_start, 5.202152, 0.048274, 1.3036, 100.5179, 13.9749, 50.1707);
      var osc_b = build_osculating_orbit(when_j2000, planet_orbit_end,   5.202223, 0.048281, 1.3036, 100.5180, 13.9899, 66.7923); //17deg per 200d 
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
      var osc_a = build_osculating_orbit(when_j2000, planet_orbit_start, 9.557734, 0.055079, 2.4874, 113.6131, 90.5280, 346.4908);
      var osc_b = build_osculating_orbit(when_j2000, planet_orbit_end,   9.551585, 0.055316, 2.4876, 113.6203, 91.1664, 353.1777); //7.0d per 200d 
      var result = current_osculating_orbit(osc_a, osc_b, when);
      return result;
    },
    add_physical: function(ephem){
      var i = degs(ephem.phase);
      var FUDGE_FACTOR_FOR_MISSING_RINGS =  -0.6+0.24; //gives the correct answer 2024-01-10
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
      return build_osculating_orbit(when_j2000, planet_orbit_start, 19.27939, 0.045165, 0.7718, 74.0400, 166.1995, 57.0765);
      
      //var osc_a = build_osculating_orbit(when_j2000, planet_orbit_start, 19.128240, 0.048830, 0.7708, 74.0658, 174.5193, 35.3703);
      //var osc_b = build_osculating_orbit(when_j2000, planet_orbit_end,   19.147300, 0.047739, 0.7706, 74.0832, 174.1226, 38.1445);  
      //var result = current_osculating_orbit(osc_a, osc_b, when);
      //return result;
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
      return build_osculating_orbit(when_j2000, planet_orbit_start, 30.28281, 0.012845, 1.7693, 131.7551, 36.624, 358.0542);
      
      //var osc_a = build_osculating_orbit(when_j2000, planet_orbit_start, 30.090770, 0.007148, 1.7711, 131.7946, 30.7960, 346.5588);
      //var osc_b = build_osculating_orbit(when_j2000, planet_orbit_end,   30.142030, 0.008522, 1.7704, 131.7801, 23.5110, 348.0054);  
      //var result = current_osculating_orbit(osc_a, osc_b, when);
      //return result;
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

  var ceres = build_minor_planet('Ceres', 1, 3.33, 0.12, {
      equinox: when_j2000, 
      epoch: when("UT 2023-09-13"),
      a: 2.767254360873952,
      e: 0.0789125317658808,
      i: rads(10.5868796009696),
      Ω: rads(80.25497772273573),
      ω: rads(73.42179714001003),
      M0: rads(60.0787728227207),
      n: rads(0.2141067962925504) 
    }
  );
  
  var vesta = build_minor_planet('Vesta', 4, 3.21, 0.32, {
      equinox: when_j2000, 
      epoch: when("UT 2023-09-13"),
      a: 2.361922083328795,
      e: 0.08944909117827099,
      i: rads(7.142176968213055),
      Ω: rads(103.7099923672372),
      ω: rads(151.6622800356213),
      M0: rads(169.3518013373338),
      n: rads(0.2715224375630777) 
    }
  );
  
  var eunomia = build_minor_planet('Eunomia', 15, 5.38, 0.23, {
      equinox: when_j2000, 
      epoch: when("UT 2023-09-13"),
      a: 2.642676431373159,
      e: 0.1871006981311412,
      i: rads(11.7548511018996),
      Ω: rads(292.9065238235527),
      ω: rads(98.75997502602053),
      M0: rads(289.9857010739341),
      n: rads(0.2294235981314534) 
    }
  );
  
  var euterpe = build_minor_planet('Euterpe', 27, 7.08, 0.15, {
      equinox: when_j2000, 
      epoch: when("UT 2023-09-13"),
      a: 2.347209309913783,
      e: 0.1713306154270552,
      i: rads(1.583331369630248),
      Ω: rads(94.76828627329486),
      ω: rads(356.3962983014516),
      M0: rads(53.9326458221896),
      n: rads(0.2740793686926032) 
    }
  );
  
  var melpomene = build_minor_planet('Melpomene', 18, 6.57, 0.25, {
      equinox: when_j2000, 
      epoch: when("UT 2023-09-13"),
      a: 2.295862691743098,
      e: 0.2183933697004635,
      i: rads(10.13154294601918),
      Ω: rads(150.3556147719062),
      ω: rads(228.1499914936565),
      M0: rads(0.3995156373718076),
      n: rads(0.2833252024749643) 
    }
  );
  
  var astraea = build_minor_planet('Astraea', 5, 7.01, 0.15, {
      equinox: when_j2000, 
      epoch: when("UT 2023-09-13"),
      a: 2.576885950717378,
      e: 0.1874390271740186,
      i: rads(5.358657374808717),
      Ω: rads(141.4689968514151),
      ω: rads(359.1367438916351),
      M0: rads(303.3980248437335),
      n: rads(0.2382655629346103) 
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
  
     {when:'UT 2024-01-01 15',text:'Moon at apogee 404909.384 km'},
     {when:'UT 2024-01-03 01',text:'Earth at perihelion 0.983306994 AU'},
     {when:'UT 2024-01-04 04',text:'Last Quarter 401013.909 km'},
     {when:'UT 2024-01-05 00',text:'Spica 1.98°S of Moon'},
     {when:'UT 2024-01-06 08',text:'Antares 6.42°S of Venus'},
     {when:'UT 2024-01-08 15',text:'Antares 0.78°S of Moon'},
     {when:'UT 2024-01-08 20',text:'Venus 5.70°N of Moon'},
     {when:'UT 2024-01-09 19',text:'Mercury 6.59°N of Moon'},
     {when:'UT 2024-01-10 09',text:'Mars 4.16°N of Moon'},
     {when:'UT 2024-01-11 12',text:'New Moon 365203.399 km'},
     {when:'UT 2024-01-12 15',text:'Mercury at greatest elongation 23.5° West'},
     {when:'UT 2024-01-13 11',text:'Moon at perigee 362266.647 km'},
     {when:'UT 2024-01-14 10',text:'Saturn 2.14°N of Moon'},
     {when:'UT 2024-01-15 20',text:'Neptune 0.95°N of Moon'},
     {when:'UT 2024-01-18 04',text:'First Quarter 374667.969 km'},
     {when:'UT 2024-01-18 21',text:'Jupiter 2.77°S of Moon'},
     {when:'UT 2024-01-19 20',text:'Uranus 2.96°S of Moon'},
     {when:'UT 2024-01-21 11',text:'Aldebaran 9.53°S of Moon'},
     {when:'UT 2024-01-24 20',text:'Pollux 1.73°N of Moon'},
     {when:'UT 2024-01-25 18',text:'Full Moon 400994.594 km'},
     {when:'UT 2024-01-27 16',text:'Mars 0.24°S of Mercury'},
     {when:'UT 2024-01-27 17',text:'Regulus 3.58°S of Moon'},
     {when:'UT 2024-01-29 08',text:'Moon at apogee 405777.270 km'},
     {when:'UT 2024-02-01 08',text:'Spica 1.67°S of Moon'},
     {when:'UT 2024-02-02 23',text:'Last Quarter 394437.592 km'},
     {when:'UT 2024-02-05 01',text:'Antares 0.56°S of Moon'},
     {when:'UT 2024-02-07 19',text:'Venus 5.43°N of Moon'},
     {when:'UT 2024-02-08 06',text:'Mars 4.21°N of Moon'},
     {when:'UT 2024-02-08 22',text:'Mercury 3.21°N of Moon'},
     {when:'UT 2024-02-09 23',text:'New Moon 358746.317 km'},
     {when:'UT 2024-02-10 19',text:'Moon at perigee 358088.072 km'},
     {when:'UT 2024-02-11 01',text:'Saturn 1.81°N of Moon'},
     {when:'UT 2024-02-12 07',text:'Neptune 0.68°N of Moon'},
     {when:'UT 2024-02-15 08',text:'Jupiter 3.16°S of Moon'},
     {when:'UT 2024-02-16 02',text:'Uranus 3.23°S of Moon'},
     {when:'UT 2024-02-16 15',text:'First Quarter 380696.948 km'},
     {when:'UT 2024-02-17 17',text:'Aldebaran 9.75°S of Moon'},
     {when:'UT 2024-02-21 02',text:'Pollux 1.62°N of Moon'},
     {when:'UT 2024-02-22 16',text:'Mars 0.63°S of Venus'},
     {when:'UT 2024-02-23 23',text:'Regulus 3.55°S of Moon'},
     {when:'UT 2024-02-24 13',text:'Full Moon 405917.286 km'},
     {when:'UT 2024-02-25 15',text:'Moon at apogee 406311.572 km'},
     {when:'UT 2024-02-28 09',text:'Mercury in superior conjunction 1.83° South'},
     {when:'UT 2024-02-28 14',text:'Saturn 0.21°N of Mercury'},
     {when:'UT 2024-02-28 21',text:'Saturn 1.62°S of Sun'},
     {when:'UT 2024-03-03 09',text:'Antares 0.35°S of Moon'},
     {when:'UT 2024-03-03 15',text:'Last Quarter 386520.139 km'},
     {when:'UT 2024-03-08 05',text:'Mars 3.52°N of Moon'},
     {when:'UT 2024-03-08 17',text:'Venus 3.28°N of Moon'},
     {when:'UT 2024-03-08 18',text:'Neptune 0.49°S of Mercury'},
     {when:'UT 2024-03-09 17',text:'Saturn 1.52°N of Moon'},
     {when:'UT 2024-03-10 07',text:'Moon at perigee 356894.957 km'},
     {when:'UT 2024-03-10 09',text:'New Moon 356901.504 km'},
     {when:'UT 2024-03-10 19',text:'Neptune 0.52°N of Moon'},
     {when:'UT 2024-03-11 03',text:'Mercury 1.02°N of Moon'},
     {when:'UT 2024-03-14 01',text:'Jupiter 3.59°S of Moon'},
     {when:'UT 2024-03-14 12',text:'Uranus 3.44°S of Moon'},
     {when:'UT 2024-03-16 00',text:'Aldebaran 9.92°S of Moon'},
     {when:'UT 2024-03-17 04',text:'First Quarter 388177.033 km'},
     {when:'UT 2024-03-17 11',text:'Neptune 1.22°S of Sun'},
     {when:'UT 2024-03-19 07',text:'Pollux 1.50°N of Moon'},
     {when:'UT 2024-03-20 03',text:'Equinox'},
     {when:'UT 2024-03-22 02',text:'Saturn 0.35°S of Venus'},
     {when:'UT 2024-03-22 05',text:'Regulus 3.61°S of Moon'},
     {when:'UT 2024-03-23 16',text:'Moon at apogee 406294.218 km'},
     {when:'UT 2024-03-24 23',text:'Mercury at greatest elongation 18.7° East'},
     {when:'UT 2024-03-25 07',text:'ECLIPSE Full Moon 405393.769 km'},
     {when:'UT 2024-03-26 20',text:'Spica 1.41°S of Moon'},
     {when:'UT 2024-03-30 15',text:'Antares 0.26°S of Moon'},
     {when:'UT 2024-04-02 03',text:'Last Quarter 379185.697 km'},
     {when:'UT 2024-04-03 11',text:'Neptune 0.29°N of Venus'},
     {when:'UT 2024-04-06 04',text:'Mars 1.98°N of Moon'},
     {when:'UT 2024-04-06 09',text:'Saturn 1.22°N of Moon'},
     {when:'UT 2024-04-07 08',text:'Neptune 0.42°N of Moon'},
     {when:'UT 2024-04-07 17',text:'Venus 0.39°S of Moon'},
     {when:'UT 2024-04-07 18',text:'Moon at perigee 358849.735 km'},
     {when:'UT 2024-04-08 18',text:'ECLIPSE New Moon 359807.205 km'},
     {when:'UT 2024-04-09 01',text:'Mercury 2.19°N of Moon'},
     {when:'UT 2024-04-10 21',text:'Jupiter 3.98°S of Moon'},
     {when:'UT 2024-04-11 00',text:'Uranus 3.55°S of Moon'},
     {when:'UT 2024-04-11 03',text:'Saturn 0.48°S of Mars'},
     {when:'UT 2024-04-11 23',text:'Mercury in inferior conjunction 2.22° North'},
     {when:'UT 2024-04-12 09',text:'Aldebaran 9.95°S of Moon'},
     {when:'UT 2024-04-15 14',text:'Pollux 1.48°N of Moon'},
     {when:'UT 2024-04-15 19',text:'First Quarter 395631.087 km'},
     {when:'UT 2024-04-18 12',text:'Regulus 3.62°S of Moon'},
     {when:'UT 2024-04-18 23',text:'Venus 1.96°S of Mercury'},
     {when:'UT 2024-04-20 02',text:'Moon at apogee 405623.444 km'},
     {when:'UT 2024-04-20 08',text:'Uranus 0.53°N of Jupiter'},
     {when:'UT 2024-04-23 03',text:'Spica 1.42°S of Moon'},
     {when:'UT 2024-04-24 00',text:'Full Moon 399781.076 km'},
     {when:'UT 2024-04-26 21',text:'Antares 0.30°S of Moon'},
     {when:'UT 2024-04-29 04',text:'Neptune 0.04°N of Mars'},
     {when:'UT 2024-05-01 11',text:'Last Quarter 373619.969 km'},
     {when:'UT 2024-05-03 23',text:'Saturn 0.84°N of Moon'},
     {when:'UT 2024-05-04 19',text:'Neptune 0.28°N of Moon'},
     {when:'UT 2024-05-05 02',text:'Mars 0.20°S of Moon'},
     {when:'UT 2024-05-05 22',text:'Moon at perigee 363163.207 km'},
     {when:'UT 2024-05-06 08',text:'Mercury 3.82°S of Moon'},
     {when:'UT 2024-05-07 16',text:'Venus 3.50°S of Moon'},
     {when:'UT 2024-05-08 03',text:'New Moon 366738.583 km'},
     {when:'UT 2024-05-08 13',text:'Uranus 3.62°S of Moon'},
     {when:'UT 2024-05-08 18',text:'Jupiter 4.33°S of Moon'},
     {when:'UT 2024-05-09 19',text:'Aldebaran 9.89°S of Moon'},
     {when:'UT 2024-05-09 21',text:'Mercury at greatest elongation 26.4° West'},
     {when:'UT 2024-05-12 23',text:'Pollux 1.59°N of Moon'},
     {when:'UT 2024-05-13 09',text:'Uranus 0.26°S of Sun'},
     {when:'UT 2024-05-15 12',text:'First Quarter 401392.001 km'},
     {when:'UT 2024-05-15 19',text:'Regulus 3.49°S of Moon'},
     {when:'UT 2024-05-17 19',text:'Moon at apogee 404639.793 km'},
     {when:'UT 2024-05-18 09',text:'Uranus 0.47°N of Venus'},
     {when:'UT 2024-05-18 19',text:'Jupiter 0.73°S of Sun'},
     {when:'UT 2024-05-20 10',text:'Spica 1.37°S of Moon'},
     {when:'UT 2024-05-23 10',text:'Jupiter 0.20°S of Venus'},
     {when:'UT 2024-05-23 14',text:'Full Moon 390646.380 km'},
     {when:'UT 2024-05-24 03',text:'Antares 0.36°S of Moon'},
     {when:'UT 2024-05-30 17',text:'Last Quarter 370419.839 km'},
     {when:'UT 2024-05-31 01',text:'Uranus 1.35°N of Mercury'},
     {when:'UT 2024-05-31 08',text:'Saturn 0.38°N of Moon'},
     {when:'UT 2024-06-01 03',text:'Neptune 0.02°N of Moon'},
     {when:'UT 2024-06-01 17',text:'Aldebaran 5.34°S of Venus'},
     {when:'UT 2024-06-02 07',text:'Moon at perigee 368102.060 km'},
     {when:'UT 2024-06-03 00',text:'Mars 2.41°S of Moon'},
     {when:'UT 2024-06-04 10',text:'Jupiter 0.12°N of Mercury'},
     {when:'UT 2024-06-04 16',text:'Venus in superior conjunction 0.06° South'},
     {when:'UT 2024-06-05 01',text:'Uranus 3.74°S of Moon'},
     {when:'UT 2024-06-05 14',text:'Jupiter 4.67°S of Moon'},
     {when:'UT 2024-06-05 18',text:'Mercury 4.66°S of Moon'},
     {when:'UT 2024-06-06 04',text:'Aldebaran 9.85°S of Moon'},
     {when:'UT 2024-06-06 13',text:'New Moon 376368.756 km'},
     {when:'UT 2024-06-06 14',text:'Venus 4.54°S of Moon'},
     {when:'UT 2024-06-08 15',text:'Aldebaran 5.45°S of Mercury'},
     {when:'UT 2024-06-09 08',text:'Pollux 1.74°N of Moon'},
     {when:'UT 2024-06-12 04',text:'Regulus 3.26°S of Moon'},
     {when:'UT 2024-06-14 05',text:'First Quarter 404000.346 km'},
     {when:'UT 2024-06-14 14',text:'Moon at apogee 404076.734 km'},
     {when:'UT 2024-06-14 17',text:'Mercury in superior conjunction 0.95° North'},
     {when:'UT 2024-06-16 18',text:'Spica 1.18°S of Moon'},
     {when:'UT 2024-06-17 13',text:'Venus 0.88°S of Mercury'},
     {when:'UT 2024-06-20 11',text:'Antares 0.33°S of Moon'},
     {when:'UT 2024-06-20 21',text:'Solstice'},
     {when:'UT 2024-06-22 01',text:'Full Moon 380041.766 km'},
     {when:'UT 2024-06-27 12',text:'Moon at perigee 369286.249 km'},
     {when:'UT 2024-06-27 15',text:'Saturn 0.08°S of Moon'},
     {when:'UT 2024-06-28 09',text:'Neptune 0.30°S of Moon'},
     {when:'UT 2024-06-28 22',text:'Last Quarter 369872.432 km'},
     {when:'UT 2024-06-29 10',text:'Pollux 4.86°N of Mercury'},
     {when:'UT 2024-07-01 18',text:'Mars 4.09°S of Moon'},
     {when:'UT 2024-07-02 10',text:'Uranus 3.95°S of Moon'},
     {when:'UT 2024-07-03 08',text:'Jupiter 5.02°S of Moon'},
     {when:'UT 2024-07-03 12',text:'Aldebaran 9.92°S of Moon'},
     {when:'UT 2024-07-05 05',text:'Earth at aphelion 1.016725489 AU'},
     {when:'UT 2024-07-05 23',text:'New Moon 387023.263 km'},
     {when:'UT 2024-07-06 15',text:'Venus 3.86°S of Moon'},
     {when:'UT 2024-07-06 16',text:'Pollux 1.84°N of Moon'},
     {when:'UT 2024-07-07 06',text:'Pollux 5.68°N of Venus'},
     {when:'UT 2024-07-07 19',text:'Mercury 3.22°S of Moon'},
     {when:'UT 2024-07-09 12',text:'Regulus 3.04°S of Moon'},
     {when:'UT 2024-07-12 08',text:'Moon at apogee 404362.069 km'},
     {when:'UT 2024-07-13 07',text:'Aldebaran 4.82°S of Jupiter'},
     {when:'UT 2024-07-13 23',text:'First Quarter 402735.027 km'},
     {when:'UT 2024-07-14 03',text:'Spica 0.91°S of Moon'},
     {when:'UT 2024-07-15 09',text:'Uranus 0.55°N of Mars'},
     {when:'UT 2024-07-17 20',text:'Antares 0.19°S of Moon'},
     {when:'UT 2024-07-21 10',text:'Full Moon 369928.915 km'},
     {when:'UT 2024-07-22 07',text:'Mercury at greatest elongation 26.9° East'},
     {when:'UT 2024-07-24 06',text:'Moon at perigee 364917.305 km'},
     {when:'UT 2024-07-24 21',text:'Saturn 0.39°S of Moon'},
     {when:'UT 2024-07-25 15',text:'Neptune 0.57°S of Moon'},
     {when:'UT 2024-07-27 12',text:'Regulus 2.65°N of Mercury'},
     {when:'UT 2024-07-28 03',text:'Last Quarter 372035.666 km'},
     {when:'UT 2024-07-29 18',text:'Uranus 4.22°S of Moon'},
     {when:'UT 2024-07-30 11',text:'Mars 5.03°S of Moon'},
     {when:'UT 2024-07-30 18',text:'Aldebaran 10.09°S of Moon'},
     {when:'UT 2024-07-31 00',text:'Jupiter 5.38°S of Moon'},
     {when:'UT 2024-08-03 00',text:'Pollux 1.82°N of Moon'},
     {when:'UT 2024-08-04 11',text:'New Moon 396838.914 km'},
     {when:'UT 2024-08-04 22',text:'Regulus 1.09°S of Venus'},
     {when:'UT 2024-08-05 19',text:'Aldebaran 4.99°S of Mars'},
     {when:'UT 2024-08-05 20',text:'Regulus 2.92°S of Moon'},
     {when:'UT 2024-08-05 22',text:'Venus 1.74°S of Moon'},
     {when:'UT 2024-08-06 00',text:'Mercury 7.47°S of Moon'},
     {when:'UT 2024-08-06 15',text:'Venus 5.93°N of Mercury'},
     {when:'UT 2024-08-09 02',text:'Moon at apogee 405297.132 km'},
     {when:'UT 2024-08-10 10',text:'Spica 0.65°S of Moon'},
     {when:'UT 2024-08-11 22',text:'Regulus 5.53°N of Mercury'},
     {when:'UT 2024-08-12 15',text:'First Quarter 397979.646 km'},
     {when:'UT 2024-08-14 05',text:'Antares 0.00°N of Moon'},
     {when:'UT 2024-08-14 17',text:'Jupiter 0.31°S of Mars'},
     {when:'UT 2024-08-19 02',text:'Mercury in inferior conjunction 4.51° South'},
     {when:'UT 2024-08-19 18',text:'Full Moon 361969.146 km'},
     {when:'UT 2024-08-21 03',text:'Saturn 0.46°S of Moon'},
     {when:'UT 2024-08-21 05',text:'Moon at perigee 360195.683 km'},
     {when:'UT 2024-08-21 22',text:'Neptune 0.69°S of Moon'},
     {when:'UT 2024-08-26 00',text:'Uranus 4.44°S of Moon'},
     {when:'UT 2024-08-26 09',text:'Last Quarter 376703.625 km'},
     {when:'UT 2024-08-26 23',text:'Aldebaran 10.27°S of Moon'},
     {when:'UT 2024-08-27 13',text:'Jupiter 5.67°S of Moon'},
     {when:'UT 2024-08-28 00',text:'Mars 5.27°S of Moon'},
     {when:'UT 2024-08-30 05',text:'Pollux 1.72°N of Moon'},
     {when:'UT 2024-09-01 09',text:'Mercury 5.03°S of Moon'},
     {when:'UT 2024-09-02 02',text:'Regulus 2.92°S of Moon'},
     {when:'UT 2024-09-03 02',text:'New Moon 403894.934 km'},
     {when:'UT 2024-09-05 03',text:'Mercury at greatest elongation 18.1° West'},
     {when:'UT 2024-09-05 10',text:'Venus 1.18°N of Moon'},
     {when:'UT 2024-09-05 15',text:'Moon at apogee 406211.127 km'},
     {when:'UT 2024-09-06 17',text:'Spica 0.52°S of Moon'},
     {when:'UT 2024-09-08 05',text:'Saturn at opposition'},
     {when:'UT 2024-09-09 07',text:'Regulus 0.50°S of Mercury'},
     {when:'UT 2024-09-10 13',text:'Antares 0.15°N of Moon'},
     {when:'UT 2024-09-11 06',text:'First Quarter 391038.639 km'},
     {when:'UT 2024-09-17 10',text:'Saturn 0.31°S of Moon'},
     {when:'UT 2024-09-17 13',text:'Spica 2.63°S of Venus'},
     {when:'UT 2024-09-18 03',text:'ECLIPSE Full Moon 357484.738 km'},
     {when:'UT 2024-09-18 08',text:'Neptune 0.67°S of Moon'},
     {when:'UT 2024-09-18 13',text:'Moon at perigee 357285.860 km'},
     {when:'UT 2024-09-21 00',text:'Neptune at opposition'},
     {when:'UT 2024-09-22 07',text:'Uranus 4.53°S of Moon'},
     {when:'UT 2024-09-22 13',text:'Equinox'},
     {when:'UT 2024-09-23 06',text:'Aldebaran 10.37°S of Moon'},
     {when:'UT 2024-09-23 23',text:'Jupiter 5.83°S of Moon'},
     {when:'UT 2024-09-24 19',text:'Last Quarter 383362.457 km'},
     {when:'UT 2024-09-25 12',text:'Mars 4.90°S of Moon'},
     {when:'UT 2024-09-26 11',text:'Pollux 1.65°N of Moon'},
     {when:'UT 2024-09-29 08',text:'Regulus 2.96°S of Moon'},
     {when:'UT 2024-09-30 21',text:'Mercury in superior conjunction 1.30° North'},
     {when:'UT 2024-10-02 19',text:'ECLIPSE New Moon 406515.237 km'},
     {when:'UT 2024-10-02 20',text:'Moon at apogee 406515.628 km'},
     {when:'UT 2024-10-03 00',text:'Mercury 1.80°N of Moon'},
     {when:'UT 2024-10-03 23',text:'Spica 0.51°S of Moon'},
     {when:'UT 2024-10-05 20',text:'Venus 3.00°N of Moon'},
     {when:'UT 2024-10-07 19',text:'Antares 0.16°N of Moon'},
     {when:'UT 2024-10-09 15',text:'Spica 2.65°S of Mercury'},
     {when:'UT 2024-10-10 19',text:'First Quarter 383541.171 km'},
     {when:'UT 2024-10-14 18',text:'Saturn 0.11°S of Moon'},
     {when:'UT 2024-10-15 18',text:'Neptune 0.59°S of Moon'},
     {when:'UT 2024-10-17 01',text:'Moon at perigee 357174.538 km'},
     {when:'UT 2024-10-17 11',text:'Full Moon 357367.473 km'},
     {when:'UT 2024-10-19 16',text:'Uranus 4.47°S of Moon'},
     {when:'UT 2024-10-20 15',text:'Aldebaran 10.34°S of Moon'},
     {when:'UT 2024-10-21 06',text:'Pollux 5.72°N of Mars'},
     {when:'UT 2024-10-21 08',text:'Jupiter 5.81°S of Moon'},
     {when:'UT 2024-10-23 18',text:'Pollux 1.70°N of Moon'},
     {when:'UT 2024-10-23 20',text:'Mars 3.91°S of Moon'},
     {when:'UT 2024-10-24 08',text:'Last Quarter 391090.884 km'},
     {when:'UT 2024-10-25 19',text:'Antares 3.11°S of Venus'},
     {when:'UT 2024-10-26 14',text:'Regulus 2.91°S of Moon'},
     {when:'UT 2024-10-29 23',text:'Moon at apogee 406161.454 km'},
     {when:'UT 2024-10-31 05',text:'Spica 0.52°S of Moon'},
     {when:'UT 2024-11-01 13',text:'New Moon 403829.297 km'},
     {when:'UT 2024-11-03 08',text:'Mercury 2.11°N of Moon'},
     {when:'UT 2024-11-04 01',text:'Antares 0.08°N of Moon'},
     {when:'UT 2024-11-05 00',text:'Venus 3.10°N of Moon'},
     {when:'UT 2024-11-09 06',text:'First Quarter 376943.730 km'},
     {when:'UT 2024-11-10 04',text:'Antares 2.02°S of Mercury'},
     {when:'UT 2024-11-11 02',text:'Saturn 0.09°S of Moon'},
     {when:'UT 2024-11-12 02',text:'Neptune 0.62°S of Moon'},
     {when:'UT 2024-11-14 11',text:'Moon at perigee 360109.291 km'},
     {when:'UT 2024-11-15 21',text:'Full Moon 361873.389 km'},
     {when:'UT 2024-11-16 01',text:'Uranus 4.36°S of Moon'},
     {when:'UT 2024-11-16 08',text:'Mercury at greatest elongation 22.6° East'},
     {when:'UT 2024-11-17 02',text:'Aldebaran 10.25°S of Moon'},
     {when:'UT 2024-11-17 03',text:'Uranus at opposition'},
     {when:'UT 2024-11-17 15',text:'Jupiter 5.64°S of Moon'},
     {when:'UT 2024-11-20 03',text:'Pollux 1.87°N of Moon'},
     {when:'UT 2024-11-20 21',text:'Mars 2.44°S of Moon'},
     {when:'UT 2024-11-22 21',text:'Regulus 2.71°S of Moon'},
     {when:'UT 2024-11-23 01',text:'Last Quarter 398384.758 km'},
     {when:'UT 2024-11-26 12',text:'Moon at apogee 405314.013 km'},
     {when:'UT 2024-11-27 12',text:'Spica 0.42°S of Moon'},
     {when:'UT 2024-12-01 06',text:'New Moon 396279.998 km'},
     {when:'UT 2024-12-01 07',text:'Antares 0.02°N of Moon'},
     {when:'UT 2024-12-02 02',text:'Mercury 4.95°N of Moon'},
     {when:'UT 2024-12-04 23',text:'Venus 2.26°N of Moon'},
     {when:'UT 2024-12-06 02',text:'Mercury in inferior conjunction 1.39° North'},
     {when:'UT 2024-12-07 21',text:'Jupiter at opposition'},
     {when:'UT 2024-12-08 09',text:'Saturn 0.31°S of Moon'},
     {when:'UT 2024-12-08 15',text:'First Quarter 372341.900 km'},
     {when:'UT 2024-12-09 09',text:'Neptune 0.83°S of Moon'},
     {when:'UT 2024-12-10 11',text:'Antares 7.13°S of Mercury'},
     {when:'UT 2024-12-12 13',text:'Moon at perigee 365360.719 km'},
     {when:'UT 2024-12-13 10',text:'Uranus 4.35°S of Moon'},
     {when:'UT 2024-12-14 13',text:'Aldebaran 10.23°S of Moon'},
     {when:'UT 2024-12-14 20',text:'Jupiter 5.47°S of Moon'},
     {when:'UT 2024-12-15 09',text:'Full Moon 370401.804 km'},
     {when:'UT 2024-12-17 13',text:'Pollux 2.05°N of Moon'},
     {when:'UT 2024-12-18 09',text:'Mars 0.91°S of Moon'},
     {when:'UT 2024-12-20 06',text:'Regulus 2.42°S of Moon'},
     {when:'UT 2024-12-21 09',text:'Solstice'},
     {when:'UT 2024-12-22 00',text:'Antares 7.09°S of Mercury'},
     {when:'UT 2024-12-22 22',text:'Last Quarter 403246.198 km'},
     {when:'UT 2024-12-24 07',text:'Moon at apogee 404484.761 km'},
     {when:'UT 2024-12-24 20',text:'Spica 0.17°S of Moon'},
     {when:'UT 2024-12-25 03',text:'Mercury at greatest elongation 22.0° West'},
     {when:'UT 2024-12-28 15',text:'Antares 0.09°N of Moon'},
     {when:'UT 2024-12-29 04',text:'Mercury 6.39°N of Moon'},
     {when:'UT 2024-12-30 22',text:'New Moon 385604.317 km'},
     {when:'UT 2025-01-03 15',text:'Venus 1.44°N of Moon'},
     {when:'UT 2025-01-04 13',text:'Earth at perihelion 0.983327407 AU'},
     {when:'UT 2025-01-04 17',text:'Saturn 0.68°S of Moon'},
     {when:'UT 2025-01-05 15',text:'Neptune 1.15°S of Moon'},
     {when:'UT 2025-01-07 00',text:'First Quarter 370409.743 km'},
     {when:'UT 2025-01-08 00',text:'Moon at perigee 370170.692 km'},
     {when:'UT 2025-01-09 16',text:'Uranus 4.49°S of Moon'},
     {when:'UT 2025-01-10 05',text:'Venus at greatest elongation 47.2° East'},
     {when:'UT 2025-01-10 21',text:'Aldebaran 10.35°S of Moon'},
     {when:'UT 2025-01-10 23',text:'Jupiter 5.43°S of Moon'},
     {when:'UT 2025-01-13 22',text:'Pollux 2.12°N of Moon'},
     {when:'UT 2025-01-14 04',text:'Mars 0.23°S of Moon'},
     {when:'UT 2025-01-16 03',text:'Mars at opposition'},
     {when:'UT 2025-01-16 16',text:'Regulus 2.21°S of Moon'},
     {when:'UT 2025-01-20 05',text:'Saturn 2.52°S of Venus'},
     {when:'UT 2025-01-21 05',text:'Spica 0.13°N of Moon'},
     {when:'UT 2025-01-21 17',text:'Pollux 2.40°N of Mars'},
     {when:'UT 2025-01-21 21',text:'Last Quarter 404019.114 km'},
     {when:'UT 2025-01-25 00',text:'Antares 0.27°N of Moon'},
     {when:'UT 2025-01-28 21',text:'Mercury 2.54°N of Moon'},
     {when:'UT 2025-01-29 13',text:'New Moon 374263.240 km'},
     {when:'UT 2025-02-01 05',text:'Saturn 1.09°S of Moon'},
     {when:'UT 2025-02-01 20',text:'Venus 2.35°N of Moon'},
     {when:'UT 2025-02-01 23',text:'Neptune 1.42°S of Moon'}
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

  /* Coords have equinox 2024.5. */
  var build_messiers = function(){
    var messier = [110];
    messier[0]=["M1",1.465968,0.38452,8.4,"Tau","NB","!! famous Crab Neb. supernova remnant","Crab Nebula"];
    messier[1]=["M2",5.649416,0.0161679,6.5,"Aqr","GC","200-mm telescope needed to resolve",""];
    messier[2]=["M3",3.5924471,0.4932374,6.2,"CVn","GC","!! contains many variable stars",""];
    messier[3]=["M4",4.2983316,-0.4640589,5.6,"Sco","GC","bright globular near Antares",""];
    messier[4]=["M5",4.0135629,0.0348246,5.6,"Ser","GC","!! one of the sky's finest globulars",""];
    messier[5]=["M6",4.6325332,-0.5624851,4.2,"Sco","OC","!! Butterfly Cluster; best at low power","Butterfly Cluster"];
    messier[6]=["M7",4.692907,-0.6077204,3.3,"Sco","OC","!! excellent in binocs or rich-field scope","Ptolemy's Cluster"];
    messier[7]=["M8",4.7355274,-0.4255221,6,"Sgr","NB","!! Lagoon Nebula w/open cl. NGC 6530","Lagoon Nebula"];
    messier[8]=["M9",4.5406301,-0.323591,7.7,"Oph","GC","smallest of Ophiuchus globulars",""];
    messier[9]=["M10",4.4435802,-0.0721972,6.6,"Oph","GC","rich globular cluster; M12 is 3°NW",""];
    messier[10]=["M11",4.941088,-0.1088409,6.3,"Sct","OC","!! Wild Duck Cl.; the best open cluster?","Wild Duck Cluster"];
    messier[11]=["M12",4.400296,-0.0347712,6.7,"Oph","GC","loose globular cluster near M10",""];
    messier[12]=["M13",4.3745624,0.6356701,5.8,"Her","GC","!! Hercules Cluster; NGC 6207 0.5°NE","Great Hercules Globular"];
    messier[13]=["M14",4.6202644,-0.0569489,7.6,"Oph","GC","200-mm telescope needed to resolve",""];
    messier[14]=["M15",5.633853,0.2142408,6.2,"Peg","GC","rich, compact globular","Great Pegasus Globular"];
    messier[15]=["M16",4.8004802,-0.2403623,6.4,"Ser","NB","Eagle Neb. w/open cl.; use neb. filter","Eagle Nebula"];
    messier[16]=["M17",4.8093127,-0.2822294,7,"Sgr","NB","!! Swan or Omega Nebula; use neb. filter","Omega Nebula"];
    messier[17]=["M18",4.8054288,-0.2988193,7.5,"Sgr","OC","sparse cluster; 1°S of M17",""];
    messier[18]=["M19",4.4685532,-0.4590222,6.8,"Oph","GC","oblate globular; M62 4°S",""];
    messier[19]=["M20",4.7302246,-0.4019728,9,"Sgr","NB","!! Trifid Nebula; look for dark lanes","Trifid Nebula"];
    messier[20]=["M21",4.7389251,-0.3926436,6.5,"Sgr","OC","0.7°NE of M20; sparse cluster",""];
    messier[21]=["M22",4.8777335,-0.4167495,5.1,"Sgr","GC","spectacular from southern latitude","Sagittarius Cluster"];
    messier[22]=["M23",4.7047258,-0.3319291,6.9,"Sgr","OC","bright, loose open cluster",""];
    messier[23]=["M24",4.7924021,-0.3227031,4.6,"Sgr","OC","rich star cloud; best in big binoculars","Sagittarius Star Cloud"];
    messier[24]=["M25",4.8565717,-0.3356413,6.5,"Sgr","OC","bright but sparse open cluster",""];
    messier[25]=["M26",4.9154759,-0.1635876,8,"Sct","OC","bright, coarse cluster",""];
    messier[26]=["M27",5.2388572,0.397672,7.4,"Vul","NB","!! Dumbbell Nebula; a superb object","Dumbbell Nebula"];
    messier[27]=["M28",4.8258657,-0.4337434,6.8,"Sgr","GC","compact globular near M22",""];
    messier[28]=["M29",5.344216,0.6739357,7.1,"Cyg","OC","small, poor open cluster 2°S of γ Cygni",""];
    messier[29]=["M30",5.6801236,-0.4026689,7.2,"Cap","GC","toughest in one-night Messier marathon",""];
    messier[30]=["M31",0.1882517,0.722579,3.4,"And","GY","!! Andromeda Gal.; look for dust lanes","Andromeda Galaxy"];
    messier[31]=["M32",0.1926185,0.7155958,8.1,"And","GY","closest companion to M31",""];
    messier[32]=["M33",0.4157623,0.5371241,5.7,"Tri","GY","large diffuse spiral; requires dark sky","Triangulum Galaxy"];
    messier[33]=["M34",0.7137766,0.7485148,5.5,"Per","OC","best at low power",""];
    messier[34]=["M35",1.6161844,0.4245966,5.3,"Gem","OC","!! look for sm. cluster NGC 2158 0.25°S",""];
    messier[35]=["M36",1.473598,0.5959785,6.3,"Aur","OC","bright but scattered group; use low pow.",""];
    messier[36]=["M37",1.5446331,0.5681753,6.2,"Aur","OC","!! finest of three Auriga clusters; very rich",""];
    messier[37]=["M38",1.4414078,0.6257252,7.4,"Aur","OC","look for small cluster NGC 1907 0.5°S",""];
    messier[38]=["M39",5.6421525,0.8472263,4.6,"Cyg","OC","very sparse cluster; use low power",""];
    messier[39]=["M40",3.2444285,1.0113768,8.4,"UMa","OC","double star Winneke 4; separation 50arcsec","Winnecke 4"];
    messier[40]=["M41",1.7804691,-0.362355,4.6,"CMa","OC","4°S of Sirius; bright but coarse",""];
    messier[41]=["M42",1.468712,-0.0948716,4,"Ori","NB","!! Orion Nebula; finest in northern sky","Great Nebula in Orion"];
    messier[42]=["M43",1.4695923,-0.091674,9,"Ori","NB","detached part of Orion Nebula","De Mairan's Nebula"];
    messier[43]=["M44",2.275503,0.3472384,3.7,"Cnc","OC","!! Beehive or Praesepe; use low power","Beehive Cluster"];
    messier[44]=["M45",0.9968481,0.4222141,1.6,"Tau","OC","!! Pleiades; look for subtle nebulosity","Pleiades"];
    messier[45]=["M46",2.0198925,-0.2596279,6,"Pup","OC","!! contains planetary nebula NGC 2438",""];
    messier[46]=["M47",1.99721,-0.254052,5.2,"Pup","OC","coarse cluster 1.5°W of M46",""];
    messier[47]=["M48",2.1598853,-0.1025465,5.5,"Hya","OC","former lost Messier; large sparse cl.",""];
    messier[48]=["M49",3.2770548,0.1372667,8.4,"Vir","GY","very bright elliptical",""];
    messier[49]=["M50",1.8517013,-0.1460982,6.3,"Mon","OC","between Sirius & Procyon; use low mag",""];
    messier[50]=["M51",3.5387844,0.8213072,8.4,"CVn","GY","!! Whirlpool Galaxy; superb in big scope","Whirlpool Galaxy"];
    messier[51]=["M52",6.1317814,1.0771844,7.3,"Cas","OC","young, rich cl.; faint Bubble Neb. nearby",""];
    messier[52]=["M53",3.4649126,0.3148089,7.6,"Com","GC","150-mm telescope needed to resolve",""];
    messier[53]=["M54",4.9596461,-0.5314598,7.6,"Sgr","GC","not easily resolved",""];
    messier[54]=["M55",5.1554915,-0.5394569,6.3,"Sgr","GC","bright, loose globular cluster",""];
    messier[55]=["M56",5.0507903,0.5275841,8.3,"Lyr","GC","within a rich starfield",""];
    messier[56]=["M57",4.9502361,0.5770967,8.8,"Lyr","NB","!! Ring Nebula; an amazing smoke ring","Ring Nebula"];
    messier[57]=["M58",3.3114865,0.2038923,9.7,"Vir","GY","bright barred spiral; M59 and M60 1°E",""];
    messier[58]=["M59",3.330241,0.2009914,9.6,"Vir","GY","bright elliptical paired with M60",""];
    messier[59]=["M60",3.337656,0.1992493,8.8,"Vir","GY","bright elliptical with M59 and NGC 4647",""];
    messier[60]=["M61",3.2426105,0.0755889,9.7,"Vir","GY","face-on two-armed spiral",""];
    messier[61]=["M62",4.4626424,-0.5262312,6.5,"Oph","GC","asymmetrical; in rich field",""];
    messier[62]=["M63",3.4771114,0.7313704,8.6,"CVn","GY","!! Sunflower Galaxy; bright, elongated","Sunflower Galaxy"];
    messier[63]=["M64",3.3942387,0.376139,8.5,"Com","GY","!! Black Eye Gal; eye needs big scope","Black Eye Galaxy"];
    messier[64]=["M65",2.9678358,0.2260036,9.3,"Leo","GY","!! bright elongated spiral",""];
    messier[65]=["M66",2.9735043,0.224256,8.9,"Leo","GY","!! M65 and NGC 3628 in same field",""];
    messier[66]=["M67",2.3201498,0.2046232,6.1,"Cnc","OC","one of the oldest star clusters known",""];
    messier[67]=["M68",3.3196327,-0.4692197,7.8,"Hya","GC","150-mm telescope needed to resolve",""];
    messier[68]=["M69",4.8563688,-0.5642806,7.6,"Sgr","GC","small, poor globular cluster",""];
    messier[69]=["M70",4.9078401,-0.5632871,7.9,"Sgr","GC","small globular 2°E of M69",""];
    messier[70]=["M71",5.213702,0.3289701,8.2,"Sge","GC","loose globular; looks like an open cluster",""];
    messier[71]=["M72",5.4752867,-0.2171079,9.3,"Aqr","GC","near the Saturn Nebula, NGC 7009",""];
    messier[72]=["M73",5.4992792,-0.2188124,9,"Aqr","OC","group of four stars only; an asterism",""];
    messier[73]=["M74",0.4276908,0.2776401,9.4,"Psc","GY","faint, elusive spiral; tough in small scope",""];
    messier[74]=["M75",5.2688964,-0.3812669,8.5,"Sgr","GC","small and distant; 59,000 ly away",""];
    messier[75]=["M76",0.4535916,0.9021515,10.1,"Per","NB","Little Dumbell; faint but distinct","Little Dumbbell Nebula"];
    messier[76]=["M77",0.715394,0.002383,8.9,"Cet","GY","a Seyfert galaxy; with starlike nucleus",""];
    messier[77]=["M78",1.5182454,0.0010043,8.3,"Ori","NB","bright featureless reflection nebula",""];
    messier[78]=["M79",1.4203032,-0.4281162,7.7,"Lep","GC","200-mm telescope needed to resolve",""];
    messier[79]=["M80",4.2693578,-0.4021622,7.3,"Sco","GC","very compressed globular",""];
    messier[80]=["M81",2.607456,1.203397,6.9,"UMa","GY","!! bright spiral visible in binoculars","Bode's Galaxy"];
    messier[81]=["M82",2.6084289,1.2141588,8.4,"UMa","GY","!! the exploding galaxy; M81 0.5°S","Cigar Galaxy"];
    messier[82]=["M83",3.5708807,-0.5234393,7.6,"Hya","GY","large and diffuse; superb from far south","Southern Pinwheel"];
    messier[83]=["M84",3.2565305,0.222491,9.1,"Vir","GY","!! w/M86 in Markarian's Chain",""];
    messier[84]=["M85",3.2582478,0.3152848,9.1,"Com","GY","bright elliptical shape",""];
    messier[85]=["M86",3.2613272,0.2236558,8.9,"Vir","GY","!! w/many NGC galaxies in Chain",""];
    messier[86]=["M87",3.2813909,0.2140626,8.6,"Vir","GY","famous jet and black hole",""];
    messier[87]=["M88",3.2870476,0.2495528,9.6,"Com","GY","bright multiple-arm spiral",""];
    messier[88]=["M89",3.3027592,0.216688,9.8,"Vir","GY","elliptical; resembles M87 but smaller",""];
    messier[89]=["M90",3.3075519,0.2274528,9.5,"Vir","GY","bright barred spiral near M89",""];
    messier[90]=["M91",3.3018735,0.2507216,10.2,"Com","GY","some lists say M91=M58, not NGC 4548",""];
    messier[91]=["M92",4.5284904,0.7523795,6.4,"Her","GC","9°NE of M13; fine but often overlooked",""];
    messier[92]=["M93",2.0317332,-0.4176059,6,"Pup","OC","compact, bright cluster; fairly rich",""];
    messier[93]=["M94",3.3687028,0.7155913,8.2,"CVn","GY","very bright and comet-like",""];
    messier[94]=["M95",2.8156174,0.2019504,9.7,"Leo","GY","bright barred spiral",""];
    messier[95]=["M96",2.8278306,0.2039774,9.2,"Leo","GY","M95 in same field",""];
    messier[96]=["M97",2.9505042,0.9578861,9.9,"UMa","NB","Owl Nebula; distinct grey oval","Owl Nebula"];
    messier[97]=["M98",3.2076819,0.2579691,10.1,"Com","GY","nearly edge-on spiral near star 6 Com. B.",""];
    messier[98]=["M99",3.2294867,0.2495373,9.9,"Com","GY","nearly face-on spiral near M98",""];
    messier[99]=["M100",3.247359,0.2739759,9.3,"Com","GY","face-on spiral with starlike nucleus",""];
    messier[100]=["M101",3.6829319,0.9465439,7.9,"UMa","GY","!! Pinwheel Gal; diffuse face-on spiral","Pinwheel Galaxy"];
    messier[101]=["M102",3.9582894,0.9716795,9.9,"Dra","GY","or is M102=M101? (look for NGC 5907)",""];
    messier[102]=["M103",0.413837,1.0615979,7.4,"Cas","OC","three NGC open clusters nearby",""];
    messier[103]=["M104",3.3216915,-0.2050924,8,"Vir","GY","!! Sombrero Galaxy; look for dust lane","Sombrero Galaxy"];
    messier[104]=["M105",2.8322021,0.2173551,9.3,"Leo","GY","bright elliptical near M95 and M96",""];
    messier[105]=["M106",3.2293197,0.8234596,8.4,"CVn","GY","!! superb large, bright spiral",""];
    messier[106]=["M107",4.3365909,-0.2286458,7.9,"Oph","GC","small, faint globular",""];
    messier[107]=["M108",2.9361703,0.9692376,10,"UMa","GY","nearly edge-on; paired with M97 0.75°SE",""];
    messier[108]=["M109",3.1366244,0.9293344,9.8,"UMa","GY","barred spiral near γ UMA",""];
    messier[109]=["M110",0.1821361,0.7298539,8.5,"And","GY","more distant companion to M31",""];
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
  
  /* Coords have equinox 2024.5. */
  var build_caldwells = function(){
    var caldwell = [42];
    caldwell[0]=["C68",4.9896806,-0.6442558,9.7,"CrA","Bn","","","NGC_6729"];
    caldwell[1]=["C69",4.5176121,-0.6479864,12.8,"Sco","Pl","","Bug Nebula","NGC_6302"];
    caldwell[2]=["C70",0.2445857,-0.655387,8.1,"Scl","Sp","","","NGC_300"];
    caldwell[3]=["C71",2.0646025,-0.6739488,5.8,"Pup","Oc","","","NGC_2477"];
    caldwell[4]=["C72",0.0703616,-0.6815031,8.2,"Scl","Sb","","","NGC_55"];
    caldwell[5]=["C73",1.374038,-0.6985349,7.3,"Col","Gc","","","NGC_1851"];
    caldwell[6]=["C74",2.6561179,-0.7077978,8.2,"Vel","Pl","","Eight Burst Nebula","NGC_3132"];
    caldwell[7]=["C75",4.3078493,-0.7107122,5.8,"Sco","Oc","","","NGC_6124"];
    caldwell[8]=["C76",4.4319332,-0.7302151,2.6,"Sco","Oc","","","NGC_6231"];
    caldwell[9]=["C77",3.5209537,-0.7529966,7,"Cen","Px","","Centaurus A","Centaurus_A"];
    caldwell[10]=["C78",4.7550476,-0.7626166,6.6,"CrA","Gc","","","NGC_6541"];
    caldwell[11]=["C79",2.6991893,-0.8122728,6.7,"Vel","Gc","","","NGC_3201"];
    caldwell[12]=["C80",3.5267781,-0.8309495,3.6,"Cen","Gc","","Omega Centauri","Omega_Centauri"];
    caldwell[13]=["C81",4.5699884,-0.8453776,8.1,"Ara","Gc","","","NGC_6352"];
    caldwell[14]=["C82",4.3770378,-0.8519314,5.2,"Ara","Oc","","","NGC_6193"];
    caldwell[15]=["C83",3.433227,-0.8656384,9.5,"Cen","Sp","","","NGC_4945"];
    caldwell[16]=["C84",3.6126749,-0.8986425,7.6,"Cen","Gc","","","NGC_5286"];
    caldwell[17]=["C85",2.2728548,-0.9277227,2.5,"Vel","Oc","","Omicron Vel Cluster","IC_2391"];
    caldwell[18]=["C86",4.636883,-0.9368499,5.6,"Ara","Gc","","","NGC_6397"];
    caldwell[19]=["C87",0.841997,-0.9621246,8.4,"Hor","Gc","","","NGC_1261"];
    caldwell[20]=["C88",3.9598737,-0.9720371,7.9,"Cir","Oc","","","NGC_5823"];
    caldwell[21]=["C89",4.2801787,-1.0115524,5.4,"Nor","Oc","","S Norma Cluster","NGC_6087"];
    caldwell[22]=["C90",2.4525869,-1.019653,9.7,"Car","Pl","","","NGC_2867"];
    caldwell[23]=["C91",2.9122977,-1.0262435,3,"Car","Oc","","","NGC_3532"];
    caldwell[24]=["C92",2.8132524,-1.0471222,6.2,"Car","Bn","","Eta Carinae Nebula","Carina_Nebula"];
    caldwell[25]=["C93",5.0311433,-1.0461713,5.4,"Pav","Gc","","","NGC_6752"];
    caldwell[26]=["C94",3.3819301,-1.0553293,4.2,"Cru","Oc","","","Jewel Box"];
    caldwell[27]=["C95",4.2141052,-1.0570715,5.1,"TrA","Oc","","","NGC_6025"];
    caldwell[28]=["C96",2.0887385,-1.0635005,3.8,"Car","Oc","","","NGC_2516"];
    caldwell[29]=["C97",3.0423394,-1.077782,5.3,"Cen","Oc","","","NGC_3766"];
    caldwell[30]=["C98",3.3325135,-1.1013144,6.9,"Cru","Oc","","","NGC_4609"];
    caldwell[31]=["C99",3.3794167,-1.1018728,undefined,"Cru","Dn","","Coalsack Nebula","Coalsack_Nebula"];
    caldwell[32]=["C100",3.0445033,-1.102508,4.5,"Cen","Oc","","Lambda Centauri Nebula","IC_2944"];
    caldwell[33]=["C101",5.0270418,-1.1136674,9,"Pav","Sb","","","NGC_6744"];
    caldwell[34]=["C102",2.8103387,-1.1262417,1.9,"Car","Oc","","Theta Car Cluster","IC_2602"];
    caldwell[35]=["C103",1.4771316,-1.2058007,1,"Dor","Bn","","Tarantula Nebula","Tarantula_Nebula"];
    caldwell[36]=["C104",0.2793695,-1.2342763,6.6,"Tuc","Gc","","","NGC_362"];
    caldwell[37]=["C105",3.4089226,-1.2394459,7.3,"Mus","Gc","","","NGC_4833"];
    caldwell[38]=["C106",0.1098483,-1.2557247,4,"Tuc","Gc","","47 Tucanae","47_Tucanae"];
    caldwell[39]=["C107",4.3136691,-1.2610654,9.3,"Aps","Gc","","","NGC_6101"];
    caldwell[40]=["C108",3.2605301,-1.2706372,7.8,"Mus","Gc","","","NGC_4372"];
    caldwell[41]=["C109",2.6580025,-1.413498,11.6,"Cha","Pl","","","NGC_3195"];    
   return caldwell;
  };
  var caldwells = build_caldwells();

  /* Find a Caldwell object either by 'C103' or 'Tarantula Nebula' (case-insensitive). */  
  var find_caldwell = function(name_raw){
    return find_deep_sky_object(name_raw, caldwells);
  };
  
  /* This exist in order to avoid creating a large number of identical objects. */
  var when_fixed_equinox = when("J2024.5");
  
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

  /* Coords have equinox 2024.5. */
  var build_stars = function(){
    /* Source: Yale Bright Star Catalog r5.  Name, Right Ascension, Declination (J2024.5), and Magnitude.*/
    var ybs = [9096];
    //removal of leading whitespace cuts down on the overall size of this js file
ybs[0]=['',0.0280765,0.7917777,6.7];
ybs[1]=['',0.0275716,-0.0064002,6.29];
ybs[2]=['33 Psc',0.0287512,-0.0972349,4.61];
ybs[3]=['86 Peg',0.0303657,0.2361859,5.51];
ybs[4]=['',0.0329396,1.0222917,5.96];
ybs[5]=['',0.0329579,-0.8541408,5.7];
ybs[6]=['10 Cas',0.0337388,1.1228129,5.59];
ybs[7]=['',0.034377,0.5088982,6.13];
ybs[8]=['',0.0352694,-0.4009227,6.18];
ybs[9]=['',0.0373202,-0.3010705,6.19];
ybs[10]=['',0.0392256,-0.0421075,6.43];
ybs[11]=['',0.0393897,-0.3904752,5.94];
ybs[12]=['',0.0405806,-0.5828203,5.68];
ybs[13]=['',0.0412617,-0.040343,6.07];
ybs[14]=['α And',0.0421323,0.5101047,2.06];
ybs[15]=['',0.0416367,-0.1516271,5.99];
ybs[16]=['',0.0434394,0.6416345,6.19];
ybs[17]=['',0.0427845,-0.3044065,6.06];
ybs[18]=['',0.0442289,0.4467878,6.23];
ybs[19]=['',0.0467958,1.3936627,6.01];
ybs[20]=['β Cas',0.0456983,1.0347358,2.27];
ybs[21]=['87 Peg',0.0449566,0.3202369,5.53];
ybs[22]=['',0.0447858,-0.9401332,6.33];
ybs[23]=['κ1 Scl',0.0462213,-0.4861006,5.42];
ybs[24]=['ε Phe',0.0464385,-0.7960662,3.88];
ybs[25]=['34 Psc',0.0493015,0.1969046,5.51];
ybs[26]=['22 And',0.0506343,0.8064899,5.03];
ybs[27]=['',0.0514519,1.000105,6.74];
ybs[28]=['',0.0504691,-0.0892276,5.84];
ybs[29]=['γ3 Oct',0.0484682,-1.4326996,5.28];
ybs[30]=['',0.0521986,-0.2171847,5.85];
ybs[31]=['',0.0515343,-1.2756299,6.64];
ybs[32]=['6 Cet',0.0545978,-0.2675911,4.89];
ybs[33]=['κ2 Scl',0.0559107,-0.4828195,5.41];
ybs[34]=['θ Scl',0.0565853,-0.6108104,5.25];
ybs[35]=['',0.0579202,0.8427967,6.16];
ybs[36]=['',0.0585233,-0.3107061,5.25];
ybs[37]=['',0.0616125,0.6602493,6.73];
ybs[38]=['γ Peg',0.0632742,0.2673802,2.83];
ybs[39]=['',0.0640219,0.473392,6.3];
ybs[40]=['23 And',0.0645702,0.7185767,5.72];
ybs[41]=['',0.0651915,-0.4517927,5.94];
ybs[42]=['',0.0653431,-0.456379,6.31];
ybs[43]=['',0.0668328,0.5819317,6.25];
ybs[44]=['χ Peg',0.0692566,0.3550482,4.8];
ybs[45]=['',0.0685515,-0.1334209,5.12];
ybs[46]=['',0.0620035,-1.4810518,5.77];
ybs[47]=['7 Cet',0.069304,-0.328064,4.44];
ybs[48]=['',0.0707117,0.3913072,6.24];
ybs[49]=['35 Psc',0.070867,0.1563277,5.79];
ybs[50]=['',0.0705019,-0.1646481,5.75];
ybs[51]=['',0.0715386,0.5527791,6.45];
ybs[52]=['',0.0717846,0.478554,6.35];
ybs[53]=['',0.0706854,-0.6068224,6.17];
ybs[54]=['',0.0770733,1.3454195,6.35];
ybs[55]=['',0.0770315,0.7632454,6.15];
ybs[56]=['',0.0758333,-0.5464688,5.67];
ybs[57]=['',0.0742655,-1.3225293,6.49];
ybs[58]=['36 Psc',0.0777981,0.146189,6.11];
ybs[59]=['',0.0797828,1.0763328,5.74];
ybs[60]=['',0.078317,-0.3503669,6.47];
ybs[61]=['',0.0805226,0.8392151,5.89];
ybs[62]=['θ And',0.0802032,0.6774959,4.61];
ybs[63]=['',0.077954,-1.3726063,6.77];
ybs[64]=['',0.0830222,0.9000491,6.14];
ybs[65]=['',0.0819614,-0.3301315,6.45];
ybs[66]=['',0.0831302,0.0318496,6.17];
ybs[67]=['σ And',0.085599,0.6443966,4.52];
ybs[68]=['',0.0853088,0.1979512,6.05];
ybs[69]=['26 And',0.0872734,0.7666712,6.11];
ybs[70]=['',0.0869271,0.5524515,5.87];
ybs[71]=['',0.0870305,-0.1381754,6.46];
ybs[72]=['',0.0869289,-0.7522259,6.33];
ybs[73]=['ι Cet',0.0902191,-0.1516344,3.56];
ybs[74]=['',0.0915896,0.713239,6.33];
ybs[75]=['',0.0933706,0.8552308,6.52];
ybs[76]=['ζ Tuc',0.0926029,-1.1299066,4.23];
ybs[77]=['',0.0946511,0.5423026,5.9];
ybs[78]=['',0.0961984,0.5767823,5.79];
ybs[79]=['41 Psc',0.0953882,0.1453177,5.37];
ybs[80]=['',0.0967524,0.1939539,6.56];
ybs[81]=['ρ And',0.0978164,0.6650471,5.18];
ybs[82]=['π Tuc',0.0949914,-1.2128151,5.51];
ybs[83]=['ι Scl',0.0992509,-0.503456,5.18];
ybs[84]=['',0.1003914,-0.3477051,5.12];
ybs[85]=['42 Psc',0.1033842,0.2376825,6.23];
ybs[86]=['',0.0981718,-1.3489855,5.97];
ybs[87]=['9 Cet',0.1051867,-0.2107269,6.39];
ybs[88]=['',0.1066035,-0.5393146,6.55];
ybs[89]=['',0.1105489,0.6756663,7.39];
ybs[90]=['',0.111665,0.9102868,5.57];
ybs[91]=['12 Cas',0.1141477,1.0815224,5.4];
ybs[92]=['',0.1123488,-0.0363655,6.07];
ybs[93]=['',0.1153836,0.9282093,5.74];
ybs[94]=['44 Psc',0.1163315,0.0362197,5.77];
ybs[95]=['β Hyi',0.116646,-1.3459746,2.8];
ybs[96]=['α Phe',0.1199088,-0.7360167,2.39];
ybs[97]=['κ Phe',0.119548,-0.7599955,3.94];
ybs[98]=['10 Cet',0.1216453,0.001496,6.19];
ybs[99]=['',0.1242202,-0.4435201,5.98];
ybs[100]=['47 Psc',0.1279593,0.3146547,5.06];
ybs[101]=['',0.1289421,0.7771909,5.17];
ybs[102]=['η Scl',0.1271478,-0.5737225,4.81];
ybs[103]=['48 Psc',0.1286642,0.2893811,6.06];
ybs[104]=['',0.1291681,0.1802058,6.04];
ybs[105]=['',0.1290757,-0.3525511,6.43];
ybs[106]=['',0.1293209,-0.6942867,5.43];
ybs[107]=['',0.1319989,0.6463873,6.26];
ybs[108]=['',0.1304186,-0.8796023,6.26];
ybs[109]=['',0.1418104,1.3466007,6.21];
ybs[110]=['',0.1383834,1.0491588,5.94];
ybs[111]=['28 And',0.1370923,0.5216237,5.23];
ybs[112]=['',0.1357062,-0.257069,6.14];
ybs[113]=['',0.1353743,-0.558182,6.57];
ybs[114]=['12 Cet',0.1365316,-0.0667073,5.72];
ybs[115]=['',0.1378885,-0.4128162,5.19];
ybs[116]=['',0.1381221,-0.7121693,6.19];
ybs[117]=['',0.1379181,-0.8391517,5.69];
ybs[118]=['13 Cas',0.1433503,1.1633405,6.18];
ybs[119]=['',0.1428248,0.588468,5.87];
ybs[120]=['λ Cas',0.1445889,0.953949,4.73];
ybs[121]=['',0.1441814,0.9245792,5.6];
ybs[122]=['λ1 Phe',0.1421825,-0.8494263,4.77];
ybs[123]=['β1 Tuc',0.1424709,-1.0964681,4.37];
ybs[124]=['β2 Tuc',0.1425359,-1.0966039,4.54];
ybs[125]=['',0.1473803,0.7614819,6.7];
ybs[126]=['',0.151899,1.2412181,6.42];
ybs[127]=['κ Cas',0.1501534,1.1007197,4.16];
ybs[128]=['52 Psc',0.1478147,0.3565605,5.38];
ybs[129]=['51 Psc',0.1468783,0.1237532,5.67];
ybs[130]=['',0.1477945,0.4837271,6.67];
ybs[131]=['',0.1488629,0.4959392,6.3];
ybs[132]=['',0.1507248,0.9604531,5.93];
ybs[133]=['β3 Tuc',0.1476156,-1.0977449,5.09];
ybs[134]=['16 Cas',0.1564919,1.1673648,6.48];
ybs[135]=['',0.1522568,-0.5135362,5.55];
ybs[136]=['θ Tuc',0.1501327,-1.2414737,6.13];
ybs[137]=['',0.1553846,-0.9117295,5.57];
ybs[138]=['',0.1579414,0.2357219,6.4];
ybs[139]=['13 Cet',0.1592561,-0.0603543,5.2];
ybs[140]=['14 Cet',0.1605778,-0.0064726,5.93];
ybs[141]=['',0.163691,0.9477705,5.08];
ybs[142]=['',0.1622624,0.2328502,6.41];
ybs[143]=['',0.1652202,1.0552387,5.79];
ybs[144]=['λ2 Phe',0.1607686,-0.8354218,5.51];
ybs[145]=['',0.1600954,-0.9544723,6.06];
ybs[146]=['',0.1642107,0.4780343,6.5];
ybs[147]=['',0.1626759,-0.2589886,6.45];
ybs[148]=['',0.162901,-0.3963267,6.06];
ybs[149]=['',0.1663288,0.7788216,5.13];
ybs[150]=['ζ Cas',0.1673336,0.9430276,3.66];
ybs[151]=['π And',0.1666654,0.590864,4.36];
ybs[152]=['53 Psc',0.1661043,0.2681916,5.89];
ybs[153]=['',0.1676202,0.4214746,6.47];
ybs[154]=['',0.1687364,0.6201848,5.48];
ybs[155]=['',0.182316,1.442133,6.4];
ybs[156]=['',0.1682468,-0.4299216,5.57];
ybs[157]=['',0.1644522,-1.1342914,6.42];
ybs[158]=['',0.1691621,0.0570686,6.39];
ybs[159]=['',0.1676928,-0.9470091,6.41];
ybs[160]=['ε And',0.1739351,0.5139309,4.37];
ybs[161]=['',0.1768495,0.8637423,5.43];
ybs[162]=['δ And',0.1773288,0.5409676,3.27];
ybs[163]=['54 Psc',0.1773953,0.3732366,5.87];
ybs[164]=['55 Psc',0.1798571,0.3765129,5.36];
ybs[165]=['α Cas',0.1828757,0.9891029,2.23];
ybs[166]=['',0.1729227,-1.2741394,6.85];
ybs[167]=['',0.1795777,-0.5903995,6.69];
ybs[168]=['',0.1790161,-0.7795057,6.01];
ybs[169]=['',0.1819665,-0.2859327,6.49];
ybs[170]=['',0.1822176,-0.4131236,6.14];
ybs[171]=['',0.1830633,-0.0736138,5.91];
ybs[172]=['32 And',0.185255,0.6910238,5.33];
ybs[173]=['',0.1812162,-1.0353332,5.89];
ybs[174]=['',0.1899886,1.1568309,5.83];
ybs[175]=['',0.1871939,0.4322003,6.04];
ybs[176]=['ξ Cas',0.1895594,0.8839487,4.8];
ybs[177]=['μ Phe',0.1853521,-0.8019941,4.59];
ybs[178]=['',0.1917379,1.0277775,6.17];
ybs[179]=['ξ Phe',0.1870908,-0.9838,5.7];
ybs[180]=['π Cas',0.1956358,0.8230727,4.94];
ybs[181]=['λ1 Scl',0.191504,-0.6689736,6.06];
ybs[182]=['',0.1909991,-1.0494407,5.98];
ybs[183]=['ρ Tuc',0.1898342,-1.1402943,5.39];
ybs[184]=['β Cet',0.1955288,-0.3115902,2.04];
ybs[185]=['',0.1999011,0.8377219,5.67];
ybs[186]=['',0.1966555,-0.2073072,6.02];
ybs[187]=['η Phe',0.1939351,-1.0005825,4.36];
ybs[188]=['21 Cas',0.2064613,1.3111203,5.66];
ybs[189]=['ο Cas',0.2011556,0.8450565,4.54];
ybs[190]=['φ1 Cet',0.1982082,-0.1828345,4.76];
ybs[191]=['λ2 Scl',0.1979794,-0.6682493,5.9];
ybs[192]=['',0.2037645,0.9661327,5.42];
ybs[193]=['',0.2005056,-0.381745,5.24];
ybs[194]=['',0.2011855,-0.7425145,5.94];
ybs[195]=['',0.1989233,-1.0884573,6.07];
ybs[196]=['',0.2103288,1.2122793,6.33];
ybs[197]=['',0.2035432,-0.0784615,6.15];
ybs[198]=['',0.2011898,-0.9351699,6.15];
ybs[199]=['18 Cet',0.2038076,-0.2224804,6.15];
ybs[200]=['',0.2079914,0.96759,6.52];
ybs[201]=['',0.2074667,0.7853099,6.05];
ybs[202]=['',0.2047209,-0.2843236,6.47];
ybs[203]=['',0.2101105,1.0421,6.39];
ybs[204]=['23 Cas',0.2157669,1.3086625,5.41];
ybs[205]=['',0.2046234,-0.827601,5.8];
ybs[206]=['',0.2068505,-0.390751,5.5];
ybs[207]=['57 Psc',0.2087271,0.2724297,5.38];
ybs[208]=['',0.2172051,1.2707445,5.87];
ybs[209]=['58 Psc',0.2107692,0.2113131,5.5];
ybs[210]=['59 Psc',0.2117205,0.3440449,6.13];
ybs[211]=['ζ And',0.2122552,0.4258715,4.06];
ybs[212]=['60 Psc',0.2123306,0.1199782,5.99];
ybs[213]=['61 Psc',0.2147325,0.3675423,6.54];
ybs[214]=['',0.2135423,-0.3129029,5.7];
ybs[215]=['η Cas',0.2205362,1.0114013,3.44];
ybs[216]=['',0.2147995,-0.3768019,5.57];
ybs[217]=['62 Psc',0.2162491,0.1297355,5.93];
ybs[218]=['',0.2166384,0.0944894,5.75];
ybs[219]=['ν Cas',0.2191932,0.8918904,4.89];
ybs[220]=['δ Psc',0.2179684,0.1347089,4.43];
ybs[221]=['64 Psc',0.2193438,0.2979935,5.07];
ybs[222]=['ν And',0.2232853,0.7192849,4.53];
ybs[223]=['',0.2210198,-0.2343668,5.59];
ybs[224]=['',0.2200633,-0.418935,5.9];
ybs[225]=['',0.2185017,-0.8127048,6.27];
ybs[226]=['65 Psc',0.2233959,0.4859682,7];
ybs[227]=['65 Psc',0.2234249,0.4859585,7.1];
ybs[228]=['',0.2214882,-0.4054142,6.28];
ybs[229]=['',0.2279183,1.1236512,5.39];
ybs[230]=['',0.2255026,0.7877588,6.15];
ybs[231]=['φ2 Cet',0.2241002,-0.1834581,5.19];
ybs[232]=['λ Hyi',0.2156264,-1.3053325,5.07];
ybs[233]=['',0.2302032,1.0810348,6.07];
ybs[234]=['',0.2284901,0.9013056,6.39];
ybs[235]=['',0.2234221,-0.7550581,6.48];
ybs[236]=['',0.2502496,1.4632762,5.62];
ybs[237]=['',0.2311392,0.9024046,6.21];
ybs[238]=['ρ Phe',0.2259914,-0.8875688,5.22];
ybs[239]=['',0.2293716,0.0613991,6.37];
ybs[240]=['',0.2380404,1.0691333,4.82];
ybs[241]=['',0.2312824,-0.7605504,6.9];
ybs[242]=['',0.2366956,0.675116,6.69];
ybs[243]=['',0.2350808,-0.4166643,5.46];
ybs[244]=['20 Cet',0.2367611,-0.0176538,4.77];
ybs[245]=['',0.2392133,0.6553827,6.06];
ybs[246]=['',0.2409344,0.921913,6.27];
ybs[247]=['',0.2373821,-0.430124,6.46];
ybs[248]=['λ1 Tuc',0.232688,-1.2107639,6.22];
ybs[249]=['υ1 Cas',0.2464246,1.0315796,4.83];
ybs[250]=['66 Psc',0.2438568,0.3372114,5.74];
ybs[251]=['21 Cet',0.2422922,-0.1502437,6.16];
ybs[252]=['',0.2464941,0.8519124,6.27];
ybs[253]=['',0.2384071,-1.0949982,5.7];
ybs[254]=['36 And',0.2455674,0.414703,5.47];
ybs[255]=['',0.2467943,0.4309096,6.2];
ybs[256]=['',0.2517025,1.0145403,6.21];
ybs[257]=['',0.2554015,1.2026751,6.37];
ybs[258]=['67 Psc',0.2500163,0.4772027,6.09];
ybs[259]=['',0.2484712,-0.1259242,5.85];
ybs[260]=['γ Cas',0.2539718,1.0620119,2.47];
ybs[261]=['υ2 Cas',0.2537111,1.0352115,4.63];
ybs[262]=['',0.2542925,1.0558352,5.55];
ybs[263]=['φ3 Cet',0.2498189,-0.1943321,5.31];
ybs[264]=['',0.2491926,-0.4824663,6.1];
ybs[265]=['μ And',0.2535828,0.6742483,3.87];
ybs[266]=['λ2 Tuc',0.2439611,-1.2111629,5.45];
ybs[267]=['η And',0.2553486,0.4110176,4.42];
ybs[268]=['',0.2576945,0.8023579,6.12];
ybs[269]=['',0.262209,1.160366,5.97];
ybs[270]=['68 Psc',0.2581731,0.5083131,5.42];
ybs[271]=['',0.2599929,0.5948562,5.98];
ybs[272]=['',0.2582992,0.2413407,6.32];
ybs[273]=['',0.2601645,0.3758802,6.37];
ybs[274]=['',0.2713667,1.2411838,6.39];
ybs[275]=['φ4 Cet',0.2616223,-0.1963173,5.61];
ybs[276]=['α Scl',0.2608577,-0.5100834,4.31];
ybs[277]=['',0.2591059,-1.0570493,6.23];
ybs[278]=['',0.2681445,0.7826535,6.84];
ybs[279]=['',0.2681592,0.7826923,6.04];
ybs[280]=['',0.2666007,0.1154489,6.11];
ybs[281]=['',0.3166735,1.5045198,4.25];
ybs[282]=['',0.4784198,1.5555358,6.46];
ybs[283]=['',0.2781443,0.8930199,6.47];
ybs[284]=['ξ Scl',0.2724611,-0.6769297,5.59];
ybs[285]=['',0.2811892,0.8291582,6.45];
ybs[286]=['39 And',0.280529,0.7238959,5.98];
ybs[287]=['σ Psc',0.2799805,0.5573821,5.5];
ybs[288]=['',0.2842567,1.0682471,5.92];
ybs[289]=['σ Scl',0.277529,-0.5483942,5.5];
ybs[290]=['ε Psc',0.2802122,0.1399961,4.28];
ybs[291]=['ω Phe',0.2751503,-0.9925888,6.11];
ybs[292]=['25 Cet',0.2805029,-0.0821265,5.43];
ybs[293]=['',0.2873932,1.0770638,5.84];
ybs[294]=['',0.2857735,0.9186228,5.99];
ybs[295]=['',0.2789022,-0.8074989,5.36];
ybs[296]=['',0.2812784,-0.5130343,6.29];
ybs[297]=['26 Cet',0.2839482,0.0261399,6.04];
ybs[298]=['',0.2889662,0.8925765,6.54];
ybs[299]=['',0.2871199,0.5199255,6.19];
ybs[300]=['',0.2777101,-1.140134,6.21];
ybs[301]=['',0.2879414,0.7002612,6.72];
ybs[302]=['',0.3545249,1.5213233,6.25];
ybs[303]=['73 Psc',0.2886242,0.1010066,6];
ybs[304]=['72 Psc',0.2896683,0.2631421,5.68];
ybs[305]=['',0.2964561,1.0976768,6.54];
ybs[306]=['ψ1 Psc',0.2923375,0.3770619,5.34];
ybs[307]=['ψ1 Psc',0.2923957,0.3769213,5.56];
ybs[308]=['',0.3115519,1.3987365,6.29];
ybs[309]=['77 Psc',0.2927322,0.0879478,6.35];
ybs[310]=['77 Psc',0.2928922,0.0879671,7.25];
ybs[311]=['27 Cet',0.2916527,-0.1718874,6.12];
ybs[312]=['',0.2989042,0.9959805,6.43];
ybs[313]=['28 Cet',0.2937116,-0.1694502,5.58];
ybs[314]=['',0.2994539,0.935999,6.38];
ybs[315]=['75 Psc',0.296061,0.2284058,6.12];
ybs[316]=['',0.2937151,-0.4164677,6.14];
ybs[317]=['μ Cas',0.3043855,0.9608132,5.17];
ybs[318]=['β Phe',0.293099,-0.8131129,3.31];
ybs[319]=['',0.2948943,-0.6201195,6.61];
ybs[320]=['41 And',0.3029296,0.769206,5.03];
ybs[321]=['',0.2984663,-0.4165389,6.37];
ybs[322]=['',0.305756,1.0191643,5.79];
ybs[323]=['78 Psc',0.3027201,0.5609931,6.25];
ybs[324]=['ψ2 Psc',0.3022482,0.3642414,5.55];
ybs[325]=['30 Cet',0.3010616,-0.1685148,5.82];
ybs[326]=['80 Psc',0.3038701,0.1008798,5.52];
ybs[327]=['υ Phe',0.3006878,-0.7218083,5.21];
ybs[328]=['ι Tuc',0.2978923,-1.0759048,5.37];
ybs[329]=['',0.3249915,1.3928313,5.64];
ybs[330]=['η Cet',0.304633,-0.1754404,3.45];
ybs[331]=['φ And',0.3095229,0.8267972,4.25];
ybs[332]=['31 Cas',0.3156577,1.2026789,5.29];
ybs[333]=['β And',0.3102579,0.6239651,2.06];
ybs[334]=['ζ Phe',0.3028514,-0.9619479,3.92];
ybs[335]=['ψ3 Psc',0.3103846,0.3453764,5.55];
ybs[336]=['44 And',0.3129358,0.736726,5.65];
ybs[337]=['',0.3126689,0.4465893,5.8];
ybs[338]=['',0.3186677,1.1228132,5.55];
ybs[339]=['θ Cas',0.3167828,0.9648087,4.33];
ybs[340]=['',0.311952,0.2758336,6.06];
ybs[341]=['32 Cas',0.3198799,1.1370562,5.57];
ybs[342]=['32 Cet',0.3116717,-0.1531732,6.4];
ybs[343]=['33 Cet',0.313387,0.0449496,5.95];
ybs[344]=['45 And',0.3165934,0.6606754,5.81];
ybs[345]=['82 Psc',0.3162189,0.5507296,5.16];
ybs[346]=['',0.3103108,-1.0046895,6.41];
ybs[347]=['χ Psc',0.3175373,0.3693889,4.66];
ybs[348]=['τ Psc',0.318584,0.5274277,4.51];
ybs[349]=['34 Cet',0.3184102,-0.0370264,5.94];
ybs[350]=['',0.3261277,1.0792277,6.41];
ybs[351]=['',0.3228697,0.7935534,6.11];
ybs[352]=['',0.3244017,0.5269773,6.19];
ybs[353]=['',0.3437998,1.3969378,6.26];
ybs[354]=['',0.3208969,-0.5353393,6.52];
ybs[355]=['',0.3223593,-0.6584588,5.92];
ybs[356]=['φ Psc',0.3276153,0.4313212,4.65];
ybs[357]=['ζ Psc',0.3272959,0.1344698,5.24];
ybs[358]=['ζ Psc',0.3273978,0.1345232,6.3];
ybs[359]=['',0.3291443,0.5001926,6.43];
ybs[360]=['87 Psc',0.3291393,0.2838396,5.98];
ybs[361]=['',0.3403629,1.2544142,7.83];
ybs[362]=['37 Cet',0.3300046,-0.1360292,5.13];
ybs[363]=['88 Psc',0.3315435,0.1243437,6.03];
ybs[364]=['38 Cet',0.3319306,-0.0147449,5.7];
ybs[365]=['',0.3397509,0.8414401,6.61];
ybs[366]=['ν Phe',0.3327504,-0.7924208,4.96];
ybs[367]=['',0.3389719,0.5802084,6.02];
ybs[368]=['',0.3426176,0.7859315,6.34];
ybs[369]=['39 Cet',0.3396979,-0.0413914,5.41];
ybs[370]=['',0.3437009,0.5562936,6.73];
ybs[371]=['',0.3597384,1.3560937,6.31];
ybs[372]=['',0.3474363,0.8298712,6.25];
ybs[373]=['κ Tuc',0.3340771,-1.1998646,4.86];
ybs[374]=['89 Psc',0.3449969,0.0653265,5.16];
ybs[375]=['',0.3498559,0.6547495,6.46];
ybs[376]=['',0.3399152,-1.1566189,6.24];
ybs[377]=['',0.3666669,1.3328457,6.38];
ybs[378]=['φ Cas',0.3562329,1.0185682,4.98];
ybs[379]=['υ Psc',0.3526392,0.478086,4.76];
ybs[380]=['35 Cas',0.3610501,1.1307309,6.34];
ybs[381]=['42 Cet',0.3536877,-0.0066463,5.87];
ybs[382]=['',0.375387,1.3762441,6.07];
ybs[383]=['',0.3570079,-0.0544372,6.23];
ybs[384]=['',0.3564036,-0.1939224,6.15];
ybs[385]=['91 Psc',0.3599036,0.5038042,5.23];
ybs[386]=['ξ And',0.3656167,0.7968549,4.88];
ybs[387]=['',0.3705842,1.0170097,6.45];
ybs[388]=['',0.3659887,0.0323565,6.2];
ybs[389]=['43 Cet',0.3657965,-0.0056237,6.49];
ybs[390]=['',0.3651985,-0.3308073,6.35];
ybs[391]=['47 And',0.3712514,0.6604719,5.58];
ybs[392]=['',0.3709462,0.5999238,6.29];
ybs[393]=['',0.3697654,0.3594717,5.97];
ybs[394]=['',0.3823027,1.241047,6.49];
ybs[395]=['ψ Cas',0.3826342,1.1913046,4.74];
ybs[396]=['',0.3693785,-0.5378797,5.84];
ybs[397]=['44 Cet',0.37206,-0.1375372,6.21];
ybs[398]=['θ Cet',0.3719774,-0.140606,3.6];
ybs[399]=['δ Cas',0.3814649,1.0535165,2.68];
ybs[400]=['',0.3733926,-0.1184658,5.91];
ybs[401]=['',0.374651,-0.2711057,6.14];
ybs[402]=['',0.3754976,-0.0475005,6.15];
ybs[403]=['',0.3793395,0.4125699,6.18];
ybs[404]=['',0.3742019,-0.7219628,5.42];
ybs[405]=['',0.3829155,0.7606923,5.96];
ybs[406]=['',0.381966,0.6057417,6.31];
ybs[407]=['',0.3741962,-0.7749481,6.26];
ybs[408]=['46 Cet',0.37884,-0.2525846,4.9];
ybs[409]=['ρ Psc',0.382145,0.3368298,5.38];
ybs[410]=['94 Psc',0.3840675,0.3380159,5.5];
ybs[411]=['',0.3861493,0.6022086,6.27];
ybs[412]=['',0.382705,-0.0047463,6.41];
ybs[413]=['ω And',0.3888638,0.7947015,4.83];
ybs[414]=['',0.3878036,0.7195466,6.46];
ybs[415]=['',0.3846712,0.0639113,6.58];
ybs[416]=['',0.3749434,-1.1212419,5.93];
ybs[417]=['47 Cet',0.3842726,-0.2256725,5.66];
ybs[418]=['',0.389269,0.7061936,6.6];
ybs[419]=['',0.3843851,-0.5657744,5.79];
ybs[420]=['α UMi',0.8010193,1.5596001,2.02];
ybs[421]=['',0.3883056,-0.1880643,6.13];
ybs[422]=['',0.3912439,0.1411556,6.2];
ybs[423]=['38 Cas',0.406151,1.2285415,5.81];
ybs[424]=['',0.4040809,1.1558211,6.14];
ybs[425]=['γ Phe',0.3901969,-0.7538438,3.41];
ybs[426]=['49 And',0.3996171,0.8226268,5.27];
ybs[427]=['',0.3919974,-0.587084,6.58];
ybs[428]=['97 Psc',0.3979669,0.3225621,6.02];
ybs[429]=['48 Cet',0.3960792,-0.3753064,5.12];
ybs[430]=['μ Psc',0.3990853,0.1094272,4.84];
ybs[431]=['',0.3950588,-0.8138536,6.31];
ybs[432]=['',0.3993917,-0.4552165,5.93];
ybs[433]=['η Psc',0.4049071,0.2700261,3.62];
ybs[434]=['',0.4081109,0.6095626,6.39];
ybs[435]=['',0.4146904,1.0201891,5.7];
ybs[436]=['δ Phe',0.4025718,-0.8542892,3.95];
ybs[437]=['',0.4051469,-0.5263537,5.82];
ybs[438]=['χ Cas',0.4169435,1.0359725,4.71];
ybs[439]=['',0.4044362,-0.7932529,6.17];
ybs[440]=['',0.4113736,-0.1551521,6.59];
ybs[441]=['',0.4102707,-0.6412352,5.51];
ybs[442]=['',0.4175689,0.6520911,5.88];
ybs[443]=['',0.4084366,-0.8657267,6.28];
ybs[444]=['',0.4142715,-0.1204324,5.76];
ybs[445]=['',0.4338984,1.2989587,6.58];
ybs[446]=['',0.4195248,0.3243744,5.89];
ybs[447]=['49 Cet',0.4181118,-0.2714218,5.63];
ybs[448]=['',0.4246634,0.7190904,6.38];
ybs[449]=['',0.4187209,-0.554447,6.12];
ybs[450]=['',0.4274547,0.8525425,5.92];
ybs[451]=['101 Psc',0.4236246,0.2580624,6.22];
ybs[452]=['40 Cas',0.4386287,1.2769481,5.28];
ybs[453]=['',0.4242869,0.3064462,5.8];
ybs[454]=['υ And',0.4287026,0.7248316,4.09];
ybs[455]=['50 Cet',0.4240182,-0.2666133,5.42];
ybs[456]=['',0.4195504,-1.0125487,6.01];
ybs[457]=['',0.4352339,1.0140605,5.56];
ybs[458]=['τ Scl',0.4244096,-0.5198125,5.69];
ybs[459]=['π Psc',0.4293634,0.2140794,5.57];
ybs[460]=['51 And',0.434187,0.8508875,3.57];
ybs[461]=['',0.4364055,0.7945402,6.36];
ybs[462]=['',0.4313003,-0.1619636,6.24];
ybs[463]=['',0.4094887,-1.3679817,6.11];
ybs[464]=['',0.4260356,-1.0148482,6.18];
ybs[465]=['χ And',0.4399629,0.7768409,4.98];
ybs[466]=['',0.4441582,0.9423329,6.39];
ybs[467]=['',0.4343407,-0.6353778,5.94];
ybs[468]=['α Eri',0.4303083,-0.9968028,0.46];
ybs[469]=['',0.4364633,-0.3691636,5.58];
ybs[470]=['',0.4362551,-0.4345552,6.7];
ybs[471]=['105 Psc',0.4407134,0.2884918,5.97];
ybs[472]=['',0.4456682,0.7578401,5.61];
ybs[473]=['τ And',0.4452161,0.7103531,4.94];
ybs[474]=['43 Cas',0.4546083,1.1897184,5.59];
ybs[475]=['',0.4352455,-0.9305239,6.84];
ybs[476]=['42 Cas',0.4575778,1.2347353,5.18];
ybs[477]=['',0.4526278,1.0674646,6.71];
ybs[478]=['',0.4535346,1.0253914,6.37];
ybs[479]=['',0.4505565,0.7458941,4.95];
ybs[480]=['',0.4480083,0.4514983,6.17];
ybs[481]=['',0.4496215,0.5265701,5.99];
ybs[482]=['',0.4393923,-0.9786846,5.87];
ybs[483]=['',0.4394215,-0.9786265,5.76];
ybs[484]=['',0.4567025,1.0741508,6.34];
ybs[485]=['ν Psc',0.4481578,0.0979233,4.44];
ybs[486]=['',0.4515235,0.6172962,5.64];
ybs[487]=['44 Cas',0.4581897,1.0589552,5.78];
ybs[488]=['',0.4492206,-0.1955066,5.75];
ybs[489]=['107 Psc',0.4530889,0.3558974,5.24];
ybs[490]=['',0.4473572,-0.6633987,6.17];
ybs[491]=['',0.4571591,0.7931615,6.34];
ybs[492]=['φ Per',0.4590634,0.8868209,4.07];
ybs[493]=['π Scl',0.4505123,-0.562066,5.25];
ybs[494]=['',0.4499861,-0.6407024,5.72];
ybs[495]=['',0.462227,1.0063339,6.21];
ybs[496]=['',0.4536356,-0.0622649,4.99];
ybs[497]=['',0.4479412,-0.8711955,6.64];
ybs[498]=['',0.4642562,0.9985263,6.25];
ybs[499]=['',0.4591991,0.5639878,6.34];
ybs[500]=['',0.4622916,0.8074243,6.35];
ybs[501]=['',0.4478333,-1.0588283,5.71];
ybs[502]=['',0.4512624,-0.9358052,5.52];
ybs[503]=['',0.4588002,-0.0810374,6.19];
ybs[504]=['109 Psc',0.4637108,0.3526477,6.27];
ybs[505]=['τ Cet',0.4592643,-0.2760252,3.5];
ybs[506]=['ο Psc',0.4655167,0.1619637,4.26];
ybs[507]=['',0.4778318,1.1165496,5.63];
ybs[508]=['',0.4248511,-1.4460191,5.87];
ybs[509]=['',0.4678343,-0.0979378,5.34];
ybs[510]=['ε Scl',0.4659535,-0.4351191,5.31];
ybs[511]=['',0.4708964,0.3060398,6.55];
ybs[512]=['τ1 Hyi',0.4424385,-1.3792475,6.33];
ybs[513]=['',0.4675133,-0.4752052,6.39];
ybs[514]=['',0.4769812,0.8089794,6.32];
ybs[515]=['',0.4671062,-0.8847855,5.49];
ybs[516]=['',0.4670089,-0.9320064,5.04];
ybs[517]=['',0.4804015,0.6645155,5.94];
ybs[518]=['4 Ari',0.4778438,0.2980474,5.84];
ybs[519]=['',0.4804467,0.5726675,5.79];
ybs[520]=['',0.4726151,-0.7267274,6.18];
ybs[521]=['',0.4202143,-1.4773398,5.69];
ybs[522]=['',0.4834423,0.8380709,5.82];
ybs[523]=['',0.4786801,0.0664411,5.91];
ybs[524]=['',0.4750111,-0.6464399,6.32];
ybs[525]=['',0.4910242,0.9085109,5.9];
ybs[526]=['1 Ari',0.4865175,0.3908845,5.86];
ybs[527]=['χ Cet',0.4834273,-0.184402,4.67];
ybs[528]=['',0.4818392,-0.5402101,6.34];
ybs[529]=['1 Per',0.4957394,0.9646034,5.52];
ybs[530]=['',0.4894446,0.1948468,5.94];
ybs[531]=['',0.4837582,-0.6681643,6.37];
ybs[532]=['2 Per',0.4962384,0.8885985,5.79];
ybs[533]=['',0.4856948,-0.8324458,6.14];
ybs[534]=['',0.4992912,0.9004973,6.26];
ybs[535]=['ζ Cet',0.4916123,-0.1782782,3.73];
ybs[536]=['',0.5037325,0.9724581,6.45];
ybs[537]=['',0.4880759,-0.874157,5.94];
ybs[538]=['ε Cas',0.5069465,1.1133368,3.38];
ybs[539]=['55 And',0.500773,0.7129598,5.4];
ybs[540]=['α Tri',0.4995361,0.5183422,3.41];
ybs[541]=['γ1 Ari',0.5012468,0.3388669,4.83];
ybs[542]=['γ2 Ari',0.5012468,0.3388281,4.75];
ybs[543]=['',0.4976171,-0.2933749,5.8];
ybs[544]=['ω Cas',0.5146142,1.2008614,4.99];
ybs[545]=['ξ Psc',0.5010208,0.0577236,4.62];
ybs[546]=['τ2 Hyi',0.4695499,-1.3972242,6.06];
ybs[547]=['',0.5078024,0.7124668,6.24];
ybs[548]=['',0.5079529,0.6500954,6.26];
ybs[549]=['β Ari',0.506128,0.3652546,2.64];
ybs[550]=['',0.4993177,-0.6715125,6.1];
ybs[551]=['ψ Phe',0.5001685,-0.8060397,4.41];
ybs[552]=['',0.512102,0.6526989,5.89];
ybs[553]=['56 And',0.5131865,0.6522418,5.67];
ybs[554]=['φ Phe',0.5034519,-0.7396239,5.11];
ybs[555]=['7 Ari',0.5114773,0.4135796,5.74];
ybs[556]=['',0.5112123,0.0343631,6.01];
ybs[557]=['',0.5249793,1.0788988,6.02];
ybs[558]=['',0.5211413,0.7297737,6.78];
ybs[559]=['ι Ari',0.5179004,0.3130458,5.1];
ybs[560]=['',0.5197948,0.4873489,5.82];
ybs[561]=['56 Cet',0.514066,-0.3910936,4.85];
ybs[562]=['χ Eri',0.5099844,-0.898665,3.7];
ybs[563]=['',0.5300042,1.1299148,5.26];
ybs[564]=['3 Per',0.5241698,0.8608397,5.69];
ybs[565]=['λ Ari',0.5205544,0.4138986,4.79];
ybs[566]=['η2 Hyi',0.5041958,-1.1785809,4.69];
ybs[567]=['',0.5085676,-1.0601503,6.06];
ybs[568]=['',0.5477301,1.3619368,6.04];
ybs[569]=['',0.5145166,-0.9014144,6.1];
ybs[570]=['',0.515447,-0.8249506,4.83];
ybs[571]=['48 Cas',0.5411431,1.2396055,4.54];
ybs[572]=['',0.5215256,-0.5750553,6.35];
ybs[573]=['',0.5277713,0.3695975,5.87];
ybs[574]=['',0.5268589,0.2166446,6.09];
ybs[575]=['',0.5471907,1.2909743,6.23];
ybs[576]=['50 Cas',0.5479641,1.2660295,3.98];
ybs[577]=['47 Cas',0.556985,1.3508424,5.38];
ybs[578]=['112 Psc',0.5298129,0.0561142,5.88];
ybs[579]=['57 Cet',0.5276144,-0.3613952,5.41];
ybs[580]=['',0.5173249,-1.139806,6.37];
ybs[581]=['υ Cet',0.5286403,-0.365818,4];
ybs[582]=['52 Cas',0.5442578,1.1347845,6];
ybs[583]=['',0.5308555,-0.146709,5.51];
ybs[584]=['',0.5264642,-0.7315107,5.57];
ybs[585]=['53 Cas',0.5447529,1.1258585,5.58];
ybs[586]=['4 Per',0.5408316,0.9530315,5.04];
ybs[587]=['α Hyi',0.5215319,-1.0725283,2.86];
ybs[588]=['49 Cas',0.5582434,1.3304832,5.22];
ybs[589]=['σ Hyi',0.5053734,-1.3653534,6.16];
ybs[590]=['π For',0.5338152,-0.5215756,5.35];
ybs[591]=['α Psc',0.5380679,0.0502817,3.82];
ybs[592]=['α Psc',0.5380679,0.0502817,4.33];
ybs[593]=['',0.5786574,1.4208821,6.05];
ybs[594]=['',0.5521268,1.1382994,6.52];
ybs[595]=['ε Tri',0.5428272,0.5829556,5.5];
ybs[596]=['',0.5250288,-1.1510144,6.1];
ybs[597]=['',0.5406578,0.2372568,5.94];
ybs[598]=['χ Phe',0.535329,-0.7783496,5.14];
ybs[599]=['γ1 And',0.5472192,0.74083,2.26];
ybs[600]=['γ2 And',0.5472702,0.7408494,4.84];
ybs[601]=['10 Ari',0.5456251,0.4546995,5.63];
ybs[602]=['',0.5391562,-0.5157061,6.42];
ybs[603]=['60 Cet',0.5430226,0.0042813,5.43];
ybs[604]=['',0.5417328,-0.2650943,5.86];
ybs[605]=['',0.5456726,0.3206193,6.21];
ybs[606]=['61 Cet',0.5456669,-0.0039008,5.93];
ybs[607]=['',0.5450261,-0.0695825,5.62];
ybs[608]=['ν For',0.5479776,-0.5092931,4.69];
ybs[609]=['κ Ari',0.5582481,0.3973109,5.03];
ybs[610]=['',0.5563349,0.1459711,6.31];
ybs[611]=['11 Ari',0.5594424,0.4506535,6.15];
ybs[612]=['',0.5573906,0.0026346,6.28];
ybs[613]=['α Ari',0.5609259,0.4115175,2];
ybs[614]=['',0.5690088,1.0216947,5.67];
ybs[615]=['',0.5676798,0.7779752,6.42];
ybs[616]=['58 And',0.5671062,0.6627791,4.82];
ybs[617]=['',0.5750438,0.941741,6.31];
ybs[618]=['β Tri',0.5716172,0.6126484,3];
ybs[619]=['14 Ari',0.5708119,0.4547407,4.98];
ybs[620]=['',0.5704257,0.3026307,6.43];
ybs[621]=['',0.5669009,-0.3082984,6.1];
ybs[622]=['',0.5919661,1.2940107,6.29];
ybs[623]=['5 Per',0.5812393,1.0081042,6.36];
ybs[624]=['59 And',0.5776017,0.6833654,5.63];
ybs[625]=['59 And',0.5776674,0.6834284,6.1];
ybs[626]=['',0.5703007,-0.422908,6.48];
ybs[627]=['15 Ari',0.5759039,0.3423445,5.7];
ybs[628]=['',0.5678149,-0.7574994,5.85];
ybs[629]=['16 Ari',0.5785788,0.4546822,6.02];
ybs[630]=['5 Tri',0.5796903,0.5522351,6.23];
ybs[631]=['64 Cet',0.5788055,0.1515664,5.63];
ybs[632]=['',0.5718378,-0.7627211,6.32];
ybs[633]=['',0.5729969,-0.885051,6.12];
ybs[634]=['',0.5784604,-0.1734477,6.01];
ybs[635]=['63 Cet',0.5796377,-0.0298618,5.93];
ybs[636]=['55 Cas',0.5953433,1.1630473,6.07];
ybs[637]=['',0.5909811,1.0240661,6.44];
ybs[638]=['6 Tri',0.5838244,0.5308785,4.94];
ybs[639]=['60 And',0.5880499,0.7739734,4.83];
ybs[640]=['',0.5847533,0.4237963,5.96];
ybs[641]=['',0.590075,0.8932496,5.31];
ybs[642]=['η Ari',0.5854453,0.372187,5.27];
ybs[643]=['',0.5917973,0.8307354,6.06];
ybs[644]=['19 Ari',0.5864007,0.2686682,5.71];
ybs[645]=['ξ1 Cet',0.5860064,0.1563906,4.37];
ybs[646]=['66 Cet',0.5848386,-0.0397879,5.54];
ybs[647]=['',0.5853649,-0.3645365,5.86];
ybs[648]=['μ For',0.5846247,-0.5342447,5.28];
ybs[649]=['',0.6002171,0.8364354,6.33];
ybs[650]=['',0.6047323,0.997766,6.48];
ybs[651]=['7 Tri',0.5994967,0.5841921,5.28];
ybs[652]=['20 Ari',0.5985194,0.4519701,5.79];
ybs[653]=['21 Ari',0.5982654,0.439055,5.58];
ybs[654]=['',0.596365,-0.1632319,6.55];
ybs[655]=['',0.591337,-0.7165147,5.91];
ybs[656]=['δ Tri',0.6044049,0.5992876,4.87];
ybs[657]=['8 Per',0.6097764,1.0124975,5.75];
ybs[658]=['7 Per',0.61008,1.0058115,5.98];
ybs[659]=['',0.6070048,0.775262,6.7];
ybs[660]=['γ Tri',0.6055352,0.5927071,4.01];
ybs[661]=['',0.6046055,0.4167887,6.55];
ybs[662]=['67 Cet',0.6030315,-0.1101246,5.51];
ybs[663]=['π1 Hyi',0.5880054,-1.1820783,5.55];
ybs[664]=['',0.6201522,1.1248395,6.6];
ybs[665]=['θ Ari',0.6086552,0.349297,5.62];
ybs[666]=['62 And',0.6146901,0.8288866,5.3];
ybs[667]=['',0.6142137,0.8130483,6.21];
ybs[668]=['',0.6077623,0.032637,5.58];
ybs[669]=['',0.6152049,0.8563797,6.37];
ybs[670]=['φ Eri',0.5994351,-0.8970898,3.56];
ybs[671]=['10 Tri',0.6125087,0.501858,5.03];
ybs[672]=['',0.6124202,0.4063061,6.46];
ybs[673]=['',0.6158375,0.6971997,6.63];
ybs[674]=['π2 Hyi',0.5933617,-1.1804223,5.69];
ybs[675]=['',0.6208518,0.8276709,6.11];
ybs[676]=['',0.6174643,0.5288311,6.47];
ybs[677]=['ο Cet',0.6134167,-0.0500169,3.04];
ybs[678]=['63 And',0.6222342,0.8772462,5.59];
ybs[679]=['',0.6112111,-0.4508825,6.34];
ybs[680]=['',0.6148454,-0.073896,6.5];
ybs[681]=['9 Per',0.6286842,0.9766195,5.17];
ybs[682]=['',0.6125574,-0.7284405,6.37];
ybs[683]=['',0.6299618,0.7244316,5.82];
ybs[684]=['',0.61391,-0.9744713,5.81];
ybs[685]=['69 Cet',0.6248346,0.0088432,5.28];
ybs[686]=['',0.6352399,0.9682133,6.28];
ybs[687]=['70 Cet',0.6259525,-0.0135132,5.42];
ybs[688]=['',0.62491,-0.1861734,5.46];
ybs[689]=['',0.6249941,-0.3063298,5.87];
ybs[690]=['64 And',0.6372919,0.8746994,5.19];
ybs[691]=['κ For',0.6268289,-0.413743,5.2];
ybs[692]=['10 Per',0.641479,0.9899437,6.25];
ybs[693]=['',0.6288128,-0.3184168,6.22];
ybs[694]=['',0.624628,-0.7520482,6.31];
ybs[695]=['65 And',0.6425929,0.8794382,4.71];
ybs[696]=['',0.6288349,-0.6539034,6.53];
ybs[697]=['',0.6273182,-0.8897975,5.92];
ybs[698]=['ξ Ari',0.6376268,0.187106,5.47];
ybs[699]=['',0.6345794,-0.4492035,6.44];
ybs[700]=['71 Cet',0.6379773,-0.046604,6.33];
ybs[701]=['δ Hyi',0.6204452,-1.1963952,4.09];
ybs[702]=['',0.6350407,-0.7108827,6.18];
ybs[703]=['ι Cas',0.6593981,1.1782835,4.52];
ybs[704]=['ρ Cet',0.641998,-0.2126004,4.89];
ybs[705]=['66 And',0.6524072,0.8845051,6.12];
ybs[706]=['',0.6421577,-0.2658429,5.83];
ybs[707]=['',0.648128,0.4733739,6.18];
ybs[708]=['11 Tri',0.6497931,0.5569389,5.54];
ybs[709]=['',0.6445655,-0.347906,5.88];
ybs[710]=['λ Hor',0.6352501,-1.0507237,5.35];
ybs[711]=['κ Hyi',0.6241526,-1.28343,5.01];
ybs[712]=['',0.6595514,0.9711796,6.51];
ybs[713]=['12 Tri',0.6528002,0.5197251,5.29];
ybs[714]=['ξ2 Cet',0.6521572,0.149551,4.28];
ybs[715]=['',0.6513014,0.0361202,6.45];
ybs[716]=['13 Tri',0.6556118,0.5243026,5.89];
ybs[717]=['κ Eri',0.6452556,-0.8306851,4.25];
ybs[718]=['',0.636832,-1.1586363,6.41];
ybs[719]=['',0.6572338,0.4114985,6.19];
ybs[720]=['φ For',0.650414,-0.5882174,5.14];
ybs[721]=['',0.6584272,0.1688424,6.07];
ybs[722]=['',0.6621637,0.592395,6.25];
ybs[723]=['',0.6529572,-0.5409468,6.11];
ybs[724]=['',0.6630233,0.4423145,5.92];
ybs[725]=['26 Ari',0.6632989,0.3484202,6.15];
ybs[726]=['',0.6590372,-0.3940036,6.77];
ybs[727]=['27 Ari',0.6644018,0.3108698,6.23];
ybs[728]=['',0.6632724,0.0063353,6];
ybs[729]=['',0.6602982,-0.4377017,6.51];
ybs[730]=['',0.648595,-1.1203429,6.37];
ybs[731]=['',0.661759,-0.3916126,6.1];
ybs[732]=['14 Tri',0.6702329,0.6327584,5.15];
ybs[733]=['',0.6665892,0.0414456,5.25];
ybs[734]=['',0.6735421,0.6047459,5.83];
ybs[735]=['75 Cet',0.6693625,-0.0161933,5.35];
ybs[736]=['σ Cet',0.6686829,-0.264199,4.75];
ybs[737]=['29 Ari',0.673037,0.2642712,6.04];
ybs[738]=['',0.6686964,-0.6339087,6.3];
ybs[739]=['',0.6999746,1.2727485,5.16];
ybs[740]=['λ1 For',0.6725564,-0.6028911,5.9];
ybs[741]=['',0.6754524,-0.3472383,6.21];
ybs[742]=['',0.6850692,0.6941237,6.36];
ybs[743]=['',0.696526,1.1493092,5.78];
ybs[744]=['',0.6857622,0.6530685,5.71];
ybs[745]=['ω For',0.6759588,-0.4908895,4.9];
ybs[746]=['15 Tri',0.6862394,0.6072577,5.35];
ybs[747]=['',0.6822898,0.1322522,6.18];
ybs[748]=['77 Cet',0.6803314,-0.1353187,5.75];
ybs[749]=['',0.6866968,0.1220452,5.82];
ybs[750]=['ν Cet',0.6857606,0.0994688,4.86];
ybs[751]=['',0.6752022,-0.8898908,6.24];
ybs[752]=['',0.6915322,0.6778533,5.9];
ybs[753]=['',0.6902083,0.5534954,6.1];
ybs[754]=['',0.691735,0.599856,5.3];
ybs[755]=['80 Cet',0.6859521,-0.1348423,5.53];
ybs[756]=['',0.6933025,0.6981497,6.54];
ybs[757]=['',0.6919667,0.5759106,6.25];
ybs[758]=['',0.6726622,-1.0904842,6.77];
ybs[759]=['31 Ari',0.6892481,0.2190914,5.68];
ybs[760]=['30 Ari',0.6910488,0.4320337,7.09];
ybs[761]=['30 Ari',0.6912526,0.4320188,6.5];
ybs[762]=['',0.6889161,0.136751,5.81];
ybs[763]=['ι1 For',0.6859659,-0.5225337,5.75];
ybs[764]=['',0.6973588,0.6602844,6.18];
ybs[765]=['',0.6981025,0.666615,6.3];
ybs[766]=['',0.695149,0.1361404,6.39];
ybs[767]=['81 Cet',0.6934719,-0.0574385,5.65];
ybs[768]=['λ2 For',0.6893801,-0.6016655,5.79];
ybs[769]=['ν Ari',0.6990645,0.3851254,5.43];
ybs[770]=['',0.7483574,1.4232992,5.78];
ybs[771]=['',0.6976606,0.0619213,6.21];
ybs[772]=['μ Hyi',0.6597306,-1.3788408,5.28];
ybs[773]=['ι2 For',0.6953609,-0.5251562,5.83];
ybs[774]=['η Hor',0.6903244,-0.9152111,5.31];
ybs[775]=['δ Cet',0.7013665,0.0075582,4.07];
ybs[776]=['',0.6955016,-0.6612294,6.49];
ybs[777]=['ε Cet',0.7013849,-0.2053868,4.84];
ybs[778]=['33 Ari',0.7073888,0.4741148,5.3];
ybs[779]=['',0.7049185,0.1084912,6.25];
ybs[780]=['',0.7042508,-0.1631686,5.78];
ybs[781]=['11 Per',0.7191471,0.9635753,5.77];
ybs[782]=['',0.7028856,-0.5328424,6.52];
ybs[783]=['',0.7187941,0.9360045,5.84];
ybs[784]=['12 Per',0.7147364,0.703319,4.91];
ybs[785]=['',0.7013171,-0.7467789,4.75];
ybs[786]=['84 Cet',0.7089676,-0.0103334,5.71];
ybs[787]=['',0.7285565,1.1855483,5.95];
ybs[788]=['',0.7185916,0.8441905,6.48];
ybs[789]=['μ Ari',0.7145049,0.3510725,5.69];
ybs[790]=['ι Eri',0.705237,-0.6937948,4.11];
ybs[791]=['',0.7114009,-0.0542758,6.05];
ybs[792]=['',0.7100466,-0.2521265,5.98];
ybs[793]=['',0.7147355,0.18928,6.3];
ybs[794]=['',0.6983866,-1.1201066,6.55];
ybs[795]=['θ Per',0.7237606,0.8609861,4.12];
ybs[796]=['14 Per',0.7229756,0.7749181,5.43];
ybs[797]=['35 Ari',0.7194942,0.4853778,4.66];
ybs[798]=['ζ Hor',0.7043338,-0.9502604,5.21];
ybs[799]=['',0.7211775,0.4492613,6.35];
ybs[800]=['γ Cet',0.7180994,0.058273,3.47];
ybs[801]=['',0.7115912,-0.668119,6.01];
ybs[802]=['ε Hyi',0.6979975,-1.1896579,4.11];
ybs[803]=['',0.711323,-0.8101984,6.1];
ybs[804]=['36 Ari',0.7229575,0.3118281,6.46];
ybs[805]=['ο Ari',0.7238883,0.2690273,5.77];
ybs[806]=['ι Hor',0.712872,-0.8848284,5.41];
ybs[807]=['π Cet',0.7212169,-0.2400866,4.25];
ybs[808]=['38 Ari',0.7256018,0.2190063,5.18];
ybs[809]=['μ Cet',0.7254551,0.1783112,4.27];
ybs[810]=['',0.7168454,-0.7055404,6.36];
ybs[811]=['',0.7468987,1.2171003,6.18];
ybs[812]=['',0.7270787,0.0840171,6.03];
ybs[813]=['',0.7215566,-0.5658775,6.22];
ybs[814]=['τ1 Eri',0.7253506,-0.3223661,4.47];
ybs[815]=['',0.7351876,0.6298035,6.25];
ybs[816]=['',0.7355482,0.6223222,6.3];
ybs[817]=['',0.7197995,-0.9157366,6.15];
ybs[818]=['',0.7249826,-0.8060795,6.85];
ybs[819]=['',0.7150214,-1.1625877,6.26];
ybs[820]=['39 Ari',0.7390139,0.512225,4.51];
ybs[821]=['',0.7476148,0.9980589,6.25];
ybs[822]=['',0.7324495,-0.3759105,6.49];
ybs[823]=['',0.7343116,-0.3906761,6.47];
ybs[824]=['40 Ari',0.7413827,0.3208698,5.82];
ybs[825]=['',0.7601134,1.2040664,5.8];
ybs[826]=['',0.7426115,0.4413733,5.86];
ybs[827]=['',0.7460813,0.6532171,6.45];
ybs[828]=['',0.737875,-0.2157122,6.9];
ybs[829]=['γ Hor',0.7242465,-1.1100674,5.74];
ybs[830]=['η Per',0.7526813,0.9773054,3.76];
ybs[831]=['η1 For',0.7354682,-0.6187104,6.51];
ybs[832]=['π Ari',0.7446607,0.3065625,5.22];
ybs[833]=['ζ Hyi',0.7239826,-1.1783487,4.84];
ybs[834]=['41 Ari',0.7480045,0.4775367,3.63];
ybs[835]=['',0.75756,1.0195201,6.45];
ybs[836]=['16 Per',0.751078,0.6705315,4.23];
ybs[837]=['β For',0.7422564,-0.5638306,4.46];
ybs[838]=['',0.7563785,0.8192835,5.88];
ybs[839]=['17 Per',0.75499,0.6136467,4.53];
ybs[840]=['γ1 For',0.7458548,-0.4269053,6.14];
ybs[841]=['γ2 For',0.7459691,-0.485922,5.39];
ybs[842]=['',0.7619279,0.9267143,6.36];
ybs[843]=['σ Ari',0.7541992,0.2649694,5.49];
ybs[844]=['η2 For',0.7471565,-0.6238391,5.92];
ybs[845]=['',0.7637373,0.8494221,6.26];
ybs[846]=['τ2 Eri',0.751154,-0.364848,4.75];
ybs[847]=['η3 For',0.749023,-0.6209187,5.47];
ybs[848]=['ν Hor',0.7398767,-1.094423,5.26];
ybs[849]=['',0.7493753,-0.6951929,6.36];
ybs[850]=['τ Per',0.7679955,0.922598,3.95];
ybs[851]=['20 Per',0.764733,0.6708389,5.33];
ybs[852]=['',0.7616709,0.2894161,6.31];
ybs[853]=['',0.7579366,-0.2211357,6.04];
ybs[854]=['',0.7546604,-0.5360756,6.4];
ybs[855]=['',0.7593727,-0.1630477,6.32];
ybs[856]=['',0.7762658,1.0754516,5.59];
ybs[857]=['',0.7786984,1.1245159,6.24];
ybs[858]=['',0.7622276,-0.3888159,5.95];
ybs[859]=['ψ For',0.761537,-0.6691248,5.92];
ybs[860]=['',0.7791844,0.8963704,6.22];
ybs[861]=['',0.7776506,0.8248675,6.02];
ybs[862]=['',0.7541909,-1.0962449,6.03];
ybs[863]=['ρ2 Ari',0.7731384,0.3216567,5.91];
ybs[864]=['',0.7622364,-0.8690248,4];
ybs[865]=['ρ3 Ari',0.7758643,0.316266,5.63];
ybs[866]=['',0.7746735,0.1479937,5.97];
ybs[867]=['',0.7631545,-0.88615,6.21];
ybs[868]=['ν Hyi',0.7432954,-1.308413,4.75];
ybs[869]=['21 Per',0.7800874,0.5590541,5.11];
ybs[870]=['η Eri',0.7750342,-0.1535954,3.89];
ybs[871]=['',0.7760373,-0.063087,5.17];
ybs[872]=['',0.7836577,0.6756509,6.04];
ybs[873]=['',0.7782543,0.0802592,6.11];
ybs[874]=['47 Ari',0.783162,0.3624275,5.8];
ybs[875]=['π Per',0.7868693,0.6939327,4.7];
ybs[876]=['',0.7627899,-1.12289,6.56];
ybs[877]=['',0.82698,1.3877412,5.49];
ybs[878]=['24 Per',0.7879709,0.6157448,4.93];
ybs[879]=['4 Eri',0.7787745,-0.4147711,5.45];
ybs[880]=['',0.7777852,-0.519373,6.29];
ybs[881]=['',0.7919703,0.8258374,5.47];
ybs[882]=['',0.790886,0.7178418,5.89];
ybs[883]=['ε Ari',0.7880965,0.374142,4.63];
ybs[884]=['ε Ari',0.7880965,0.374142,4.63];
ybs[885]=['6 Eri',0.7818355,-0.4103111,5.84];
ybs[886]=['',0.7968758,0.9153794,5.28];
ybs[887]=['',0.7969632,0.9153889,6.74];
ybs[888]=['',0.7851247,-0.0468755,5.23];
ybs[889]=['',0.7788548,-0.6648628,6.41];
ybs[890]=['',0.7930646,0.6671992,6.11];
ybs[891]=['',0.7853097,-0.1689424,6.14];
ybs[892]=['λ Cet',0.7898993,0.1571458,4.7];
ybs[893]=['θ1 Eri',0.7818749,-0.7017575,3.24];
ybs[894]=['θ2 Eri',0.7819185,-0.7017527,4.35];
ybs[895]=['5 Eri',0.7894389,-0.0413413,5.56];
ybs[896]=['',0.7860677,-0.5028354,6.14];
ybs[897]=['ζ For',0.7883456,-0.4394351,5.71];
ybs[898]=['',0.7944109,0.1913951,5.95];
ybs[899]=['',0.7882278,-0.5656758,6.31];
ybs[900]=['7 Eri',0.7945023,-0.0485687,6.11];
ybs[901]=['49 Ari',0.8000239,0.4635168,5.9];
ybs[902]=['',0.8539254,1.4235091,5.95];
ybs[903]=['ρ1 Eri',0.795741,-0.1320704,5.75];
ybs[904]=['',0.7992042,0.0947974,6.25];
ybs[905]=['β Hor',0.7821864,-1.1165663,4.99];
ybs[906]=['93 Cet',0.8013714,0.0776313,5.61];
ybs[907]=['α Cet',0.8009488,0.0730408,2.53];
ybs[908]=['',0.7990228,-0.1721944,5.83];
ybs[909]=['',0.8000876,-0.1116914,6.19];
ybs[910]=['ε For',0.7970771,-0.4886247,5.89];
ybs[911]=['γ Per',0.8141419,0.9355037,2.93];
ybs[912]=['',0.8070943,0.4950516,6.36];
ybs[913]=['ρ2 Eri',0.8024511,-0.1324745,5.32];
ybs[914]=['',0.8176893,0.9913386,4.76];
ybs[915]=['τ3 Eri',0.8005689,-0.4106627,4.09];
ybs[916]=['',0.8181729,0.9802161,6.11];
ybs[917]=['ρ Per',0.8148572,0.6795297,3.39];
ybs[918]=['',0.8263935,1.1196399,5.89];
ybs[919]=['',0.8156889,0.7099358,6.05];
ybs[920]=['',0.8117809,0.2783848,6.49];
ybs[921]=['ρ3 Eri',0.8092952,-0.1310125,5.26];
ybs[922]=['',0.8111587,0.0341703,6.05];
ybs[923]=['52 Ari',0.8154524,0.4424252,6.8];
ybs[924]=['52 Ari',0.8154524,0.4424252,7];
ybs[925]=['',0.8018423,-0.8182099,5.82];
ybs[926]=['',0.8283255,0.9129109,6.31];
ybs[927]=['',0.8191885,0.231791,5.62];
ybs[928]=['',0.8493533,1.2999964,4.87];
ybs[929]=['',0.8267628,0.8273049,6.41];
ybs[930]=['μ Hor',0.8037151,-1.0409665,5.11];
ybs[931]=['',0.8193097,-0.1046364,5.27];
ybs[932]=['β Per',0.8280378,0.7164255,2.12];
ybs[933]=['ι Per',0.8325041,0.8675248,4.05];
ybs[934]=['53 Ari',0.8238529,0.3136875,6.11];
ybs[935]=['θ Hyi',0.7955363,-1.253269,5.53];
ybs[936]=['54 Ari',0.8279222,0.3296502,6.27];
ybs[937]=['κ Per',0.8340657,0.784512,3.8];
ybs[938]=['',0.8288602,0.1494575,6.28];
ybs[939]=['',0.8242033,-0.4841238,6.19];
ybs[940]=['55 Ari',0.8337952,0.509094,5.72];
ybs[941]=['',0.8360816,0.4871521,6.42];
ybs[942]=['',0.8373702,0.4710297,6.02];
ybs[943]=['ω Per',0.8416061,0.6929463,4.63];
ybs[944]=['',0.8377117,0.2088124,5.98];
ybs[945]=['',0.8471159,0.8345558,6.33];
ybs[946]=['',0.845558,0.7411881,6.15];
ybs[947]=['δ Ari',0.8422594,0.3458857,4.35];
ybs[948]=['',0.8408779,0.2293193,6.12];
ybs[949]=['',0.8363112,-0.4127125,6.38];
ybs[950]=['56 Ari',0.8451864,0.4773088,5.79];
ybs[951]=['',0.8401245,-0.0649326,6.05];
ybs[952]=['',0.8513434,0.8424217,5.9];
ybs[953]=['',0.8395891,-0.2780998,6.26];
ybs[954]=['',0.8453658,0.117838,5.56];
ybs[955]=['',0.820412,-1.2072879,6.15];
ybs[956]=['',0.8344992,-0.84897,6.12];
ybs[957]=['',0.8880545,1.3582416,5.45];
ybs[958]=['94 Cet',0.8465756,-0.019294,5.06];
ybs[959]=['α For',0.8425692,-0.5043292,3.87];
ybs[960]=['',0.8626119,0.9988516,5.79];
ybs[961]=['',0.9537283,1.4833824,5.61];
ybs[962]=['',0.8577347,0.7433965,6.07];
ybs[963]=['',0.8713706,1.1475001,6.36];
ybs[964]=['',0.8433766,-0.7736841,5.93];
ybs[965]=['',0.8638039,0.8905852,5.03];
ybs[966]=['',0.8464213,-0.6257579,6.27];
ybs[967]=['',0.8588808,0.5348755,5.52];
ybs[968]=['ζ Ari',0.8565899,0.3688597,4.89];
ybs[969]=['',0.8628597,0.7929887,6.16];
ybs[970]=['',0.8493439,-0.5186045,6.16];
ybs[971]=['',0.8609097,0.5750097,6.31];
ybs[972]=['',0.8620779,0.606986,6.25];
ybs[973]=['',0.8428889,-0.9988656,5.74];
ybs[974]=['',0.8643813,0.5632612,6.06];
ybs[975]=['',0.8674303,0.7081136,6.45];
ybs[976]=['',0.8554648,-0.4539703,6.25];
ybs[977]=['',0.8148692,-1.3769958,5.57];
ybs[978]=['30 Per',0.8702571,0.7699224,5.47];
ybs[979]=['',0.8605706,-0.1017424,6.17];
ybs[980]=['ζ Eri',0.8596853,-0.1523748,4.8];
ybs[981]=['',0.8821541,1.147369,4.84];
ybs[982]=['',0.8698705,0.6871651,5.96];
ybs[983]=['29 Per',0.8743559,0.8780772,5.15];
ybs[984]=['14 Eri',0.8629971,-0.1582227,6.14];
ybs[985]=['31 Per',0.8765171,0.8758529,5.03];
ybs[986]=['',0.860439,-0.536485,6.65];
ybs[987]=['',0.8738424,0.5988344,4.82];
ybs[988]=['95 Cet',0.8710179,-0.0146982,5.38];
ybs[989]=['',0.8686253,-0.5010599,5.91];
ybs[990]=['15 Eri',0.870273,-0.391359,4.88];
ybs[991]=['59 Ari',0.8787743,0.4740049,5.9];
ybs[992]=['κ1 Cet',0.8754672,0.0603526,4.83];
ybs[993]=['',0.8718041,-0.3223921,5.71];
ybs[994]=['',0.8649994,-0.8318765,5.85];
ybs[995]=['',0.8806459,0.5085107,4.47];
ybs[996]=['60 Ari',0.8808869,0.449421,6.12];
ybs[997]=['',0.8884595,0.8579558,5.93];
ybs[998]=['32 Per',0.8861712,0.7577536,4.95];
ybs[999]=['τ4 Eri',0.8753083,-0.3782152,3.69];
ybs[1000]=['',0.8754957,-0.4194974,5.61];
ybs[1001]=['τ1 Ari',0.884208,0.3705985,5.28];
ybs[1002]=['ζ1 Ret',0.8649315,-1.0905985,5.54];
ybs[1003]=['κ2 Cet',0.8831208,0.0656667,5.69];
ybs[1004]=['',0.8761266,-0.7501812,4.27];
ybs[1005]=['',0.9024675,1.1287243,5.23];
ybs[1006]=['ζ2 Ret',0.8668776,-1.0893997,5.24];
ybs[1007]=['',0.8943418,0.8604321,5.29];
ybs[1008]=['62 Ari',0.8887014,0.4833485,5.52];
ybs[1009]=['',0.8805148,-0.4628491,6.39];
ybs[1010]=['',0.8651019,-1.1665504,6.05];
ybs[1011]=['τ2 Ari',0.8908582,0.3635177,5.09];
ybs[1012]=['',0.8834503,-0.4109987,5.52];
ybs[1013]=['α Per',0.8992144,0.871729,1.79];
ybs[1014]=['',0.8871752,-0.4450832,6.35];
ybs[1015]=['',0.8989905,0.5867986,5.61];
ybs[1016]=['',0.9060474,0.9425866,6.51];
ybs[1017]=['',0.8829079,-0.8323504,6.39];
ybs[1018]=['64 Ari',0.897799,0.4330079,5.5];
ybs[1019]=['',0.8942366,0.0867017,6.38];
ybs[1020]=['',0.8922763,-0.1345349,6.2];
ybs[1021]=['ι Hyi',0.8525251,-1.3491168,5.52];
ybs[1022]=['',0.9022804,0.7215566,6.51];
ybs[1023]=['65 Ari',0.8982024,0.3645804,6.08];
ybs[1024]=['',0.8967487,0.2219167,6.04];
ybs[1025]=['',0.9062963,0.8587955,6.09];
ybs[1026]=['ο Tau',0.8994428,0.15907,3.6];
ybs[1027]=['',0.89329,-0.5693525,6.5];
ybs[1028]=['',0.9290061,1.2556973,6.32];
ybs[1029]=['',0.9180706,1.0531119,6.49];
ybs[1030]=['',0.9154484,0.857765,4.98];
ybs[1031]=['',0.9209813,1.0476038,4.21];
ybs[1032]=['',0.9095577,0.3288283,6.57];
ybs[1033]=['',0.919101,0.8714687,5.58];
ybs[1034]=['ξ Tau',0.9097517,0.171336,3.74];
ybs[1035]=['',0.9104722,0.2237334,6.28];
ybs[1036]=['',0.9245488,1.0290674,4.54];
ybs[1037]=['',0.9158258,0.5915087,5.61];
ybs[1038]=['χ1 For',0.9026687,-0.6254581,6.39];
ybs[1039]=['',0.9258057,1.0375735,6.13];
ybs[1040]=['34 Per',0.9212378,0.8655402,4.67];
ybs[1041]=['',0.9049963,-0.4753056,5.93];
ybs[1042]=['',0.9245474,0.9692604,5.09];
ybs[1043]=['',0.9213448,0.8206656,6.24];
ybs[1044]=['66 Ari',0.9157781,0.3994639,6.03];
ybs[1045]=['',0.903517,-0.725225,6.32];
ybs[1046]=['',0.9127492,-0.1955288,5.73];
ybs[1047]=['',0.9265863,0.8410033,5.82];
ybs[1048]=['σ Per',0.9263962,0.8391129,4.36];
ybs[1049]=['',0.8908071,-1.2136833,6.15];
ybs[1050]=['χ2 For',0.9097707,-0.6212923,5.71];
ybs[1051]=['',0.950962,1.2815401,6.57];
ybs[1052]=['',0.9304874,0.8603013,6.29];
ybs[1053]=['χ3 For',0.9125299,-0.6242985,6.5];
ybs[1054]=['',0.9319374,0.863634,6.39];
ybs[1055]=['',0.9200329,-0.1173226,5.99];
ybs[1056]=['4 Tau',0.92394,0.1992976,5.14];
ybs[1057]=['',0.919608,-0.2197679,5.59];
ybs[1058]=['',0.9332523,0.8395943,5.47];
ybs[1059]=['',0.8976697,-1.2086636,5.96];
ybs[1060]=['',0.9286476,0.4826531,5.96];
ybs[1061]=['5 Tau',0.926025,0.2272238,4.11];
ybs[1062]=['',0.9252856,0.1094492,5.94];
ybs[1063]=['',0.9403671,1.0270543,6.4];
ybs[1064]=['36 Per',0.9344,0.8052672,5.31];
ybs[1065]=['17 Eri',0.9243078,-0.0871412,4.73];
ybs[1066]=['',0.9409213,1.0114131,6.37];
ybs[1067]=['',0.9352506,0.7842973,6.41];
ybs[1068]=['',0.9404428,0.960901,5.98];
ybs[1069]=['',0.9347751,0.6203434,5.9];
ybs[1070]=['',0.919675,-0.7426602,5.78];
ybs[1071]=['',0.9211014,-0.720599,6.12];
ybs[1072]=['',0.9469987,1.0493142,6.46];
ybs[1073]=['',0.9390238,0.6977894,5.81];
ybs[1074]=['6 Tau',0.9334381,0.1650226,5.77];
ybs[1075]=['',0.9705665,1.3232651,6.27];
ybs[1076]=['',0.9224116,-0.8254137,5.99];
ybs[1077]=['',0.9291489,-0.4456223,6.38];
ybs[1078]=['κ Ret',0.9153794,-1.097014,4.72];
ybs[1079]=['ε Eri',0.9342445,-0.1636591,3.73];
ybs[1080]=['',0.9404525,0.3126497,6.17];
ybs[1081]=['7 Tau',0.9420402,0.4283913,5.92];
ybs[1082]=['ψ Per',0.9522606,0.8425104,4.23];
ybs[1083]=['τ5 Eri',0.9375485,-0.3761499,4.27];
ybs[1084]=['',0.9430258,0.1134151,6.49];
ybs[1085]=['',0.930734,-0.8778478,5.68];
ybs[1086]=['',0.9416172,-0.1708339,6.25];
ybs[1087]=['',0.9211821,-1.1590236,5.83];
ybs[1088]=['',0.9378441,-0.5410409,6.2];
ybs[1089]=['',0.9611093,0.9950358,6.3];
ybs[1090]=['',0.9404742,-0.5549117,6.4];
ybs[1091]=['',0.9308172,-1.063523,6.41];
ybs[1092]=['',0.9584862,0.7445898,6.42];
ybs[1093]=['',0.9474091,-0.1939707,5.57];
ybs[1094]=['',0.9514177,0.011646,5.71];
ybs[1095]=['20 Eri',0.9486161,-0.3034636,5.23];
ybs[1096]=['10 Tau',0.9517823,0.008397,4.28];
ybs[1097]=['',0.9563334,0.2706971,6.39];
ybs[1098]=['',0.9618002,0.3664181,6.5];
ybs[1099]=['',0.9367759,-1.1463947,6.75];
ybs[1100]=['',0.9787205,1.1046768,5.1];
ybs[1101]=['',0.9510972,-0.7015402,4.58];
ybs[1102]=['',1.1329526,1.5103431,5.86];
ybs[1103]=['',0.958556,-0.1276357,5.85];
ybs[1104]=['',0.9125748,-1.3660466,5.7];
ybs[1105]=['',0.9634964,0.2899837,6.16];
ybs[1106]=['21 Eri',0.9609363,-0.0968256,5.96];
ybs[1107]=['',0.9806587,1.0479978,5.76];
ybs[1108]=['',0.9718595,0.6572437,5.57];
ybs[1109]=['τ For',0.9591219,-0.4863269,6.01];
ybs[1110]=['12 Tau',0.9648686,0.054715,5.57];
ybs[1111]=['',0.9637177,-0.0578567,6.23];
ybs[1112]=['',0.9625361,-0.1807985,6.19];
ybs[1113]=['11 Tau',0.969706,0.4434348,6.11];
ybs[1114]=['',0.9653368,-0.0181972,6.12];
ybs[1115]=['',0.9657085,-0.2643965,6.33];
ybs[1116]=['22 Eri',0.9680182,-0.0895866,5.53];
ybs[1117]=['δ Per',0.9803501,0.8353821,3.01];
ybs[1118]=['40 Per',0.9771011,0.5941396,4.97];
ybs[1119]=['',0.9964032,1.1741939,5.8];
ybs[1120]=['',0.9703687,-0.2046522,6.49];
ybs[1121]=['13 Tau',0.9762173,0.3451743,5.69];
ybs[1122]=['',0.9855709,0.8482195,6.06];
ybs[1123]=['',0.9707061,-0.3404689,6.59];
ybs[1124]=['',0.9957313,1.106883,4.8];
ybs[1125]=['',0.9879006,0.8059098,6.11];
ybs[1126]=['ο Per',0.9855053,0.5648595,3.83];
ybs[1127]=['14 Tau',0.9826403,0.3445459,6.14];
ybs[1128]=['',0.9866128,0.6376668,5.59];
ybs[1129]=['δ For',0.9739965,-0.556087,5];
ybs[1130]=['ν Per',0.989897,0.7444507,3.77];
ybs[1131]=['δ Eri',0.9792455,-0.1690697,3.54];
ybs[1132]=['',0.9856658,0.3665942,6.1];
ybs[1133]=['',1.0114926,1.2382087,5.44];
ybs[1134]=['',0.9805941,-0.1816776,5.6];
ybs[1135]=['16 Tau',0.9872649,0.4252488,5.46];
ybs[1136]=['',0.9935797,0.7986069,5.66];
ybs[1137]=['17 Tau',0.9875704,0.4221745,3.7];
ybs[1138]=['',0.9762811,-0.6499081,4.59];
ybs[1139]=['18 Tau',0.9888523,0.43484,5.64];
ybs[1140]=['19 Tau',0.9890405,0.4283481,4.3];
ybs[1141]=['24 Eri',0.9850425,-0.0189777,5.25];
ybs[1142]=['',1.0012547,0.977324,6.1];
ybs[1143]=['γ Cam',1.0165894,1.2462465,4.63];
ybs[1144]=['20 Tau',0.9917359,0.4266071,3.87];
ybs[1145]=['25 Eri',0.9869633,-0.0038601,5.55];
ybs[1146]=['21 Tau',0.992093,0.4298691,5.76];
ybs[1147]=['22 Tau',0.9927104,0.4294025,6.43];
ybs[1148]=['29 Tau',0.9903766,0.1069035,5.35];
ybs[1149]=['',0.9408505,-1.3655964,6.29];
ybs[1150]=['',1.0113754,1.1449197,4.47];
ybs[1151]=['23 Tau',0.9939014,0.4192821,4.18];
ybs[1152]=['',0.9816188,-0.7083291,6.45];
ybs[1153]=['',1.0113545,1.1060179,5.85];
ybs[1154]=['',0.9925196,0.1200426,5.91];
ybs[1155]=['',1.0041016,0.8868079,6.14];
ybs[1156]=['',1.009218,0.9981745,6.46];
ybs[1157]=['π Eri',0.9917831,-0.2099062,4.42];
ybs[1158]=['',1.0011093,0.5877214,6.57];
ybs[1159]=['',1.0007705,0.5632002,6.25];
ybs[1160]=['η Tau',0.998965,0.4220063,2.87];
ybs[1161]=['',1.0215958,1.1969349,6.32];
ybs[1162]=['',0.9843416,-0.8375088,6.49];
ybs[1163]=['',0.9825743,-0.9459345,6.3];
ybs[1164]=['',0.986225,-0.8252661,5.73];
ybs[1165]=['',1.0072151,0.7685793,6.02];
ybs[1166]=['σ For',0.9924655,-0.5107401,5.9];
ybs[1167]=['',1.0027056,0.4100627,5.45];
ybs[1168]=['τ6 Eri',0.9944369,-0.4044822,4.23];
ybs[1169]=['30 Tau',1.0018977,0.1957762,5.07];
ybs[1170]=['β Ret',0.9795436,-1.1297664,3.85];
ybs[1171]=['',1.0113839,0.7861067,5.66];
ybs[1172]=['42 Per',1.0083621,0.57883,5.11];
ybs[1173]=['27 Tau',1.0062829,0.4210899,3.63];
ybs[1174]=['',0.9962693,-0.5205894,6.55];
ybs[1175]=['28 Tau',1.0063955,0.4225441,5.09];
ybs[1176]=['τ7 Eri',0.9979505,-0.4153976,5.24];
ybs[1177]=['',1.0031549,0.0052609,5.91];
ybs[1178]=['',1.0087279,0.4151218,6.17];
ybs[1179]=['ρ For',0.9988657,-0.5252343,5.54];
ybs[1180]=['',1.0095039,0.3895123,6.07];
ybs[1181]=['',0.998105,-0.6288716,6.21];
ybs[1182]=['',1.0021566,-0.3635405,5.81];
ybs[1183]=['',1.0113833,0.4477153,5.26];
ybs[1184]=['',1.0013498,-0.6553441,5.4];
ybs[1185]=['',1.0013863,-0.6553152,4.73];
ybs[1186]=['',1.0186965,0.6009361,5.77];
ybs[1187]=['',1.0285381,1.013092,5.8];
ybs[1188]=['',1.0168875,0.3857837,6.83];
ybs[1189]=['',1.0150248,0.2289547,6.3];
ybs[1190]=['',1.0051982,-0.6305342,4.17];
ybs[1191]=['',1.0436765,1.2547342,6.34];
ybs[1192]=['',1.0193084,0.5452441,6.25];
ybs[1193]=['',1.0272603,0.8503515,5.76];
ybs[1194]=['31 Tau',1.0180175,0.1153081,5.67];
ybs[1195]=['',1.0102936,-0.6344714,6.86];
ybs[1196]=['',1.0234949,0.3036572,5.97];
ybs[1197]=['30 Eri',1.0206064,-0.0923239,5.48];
ybs[1198]=['ζ Per',1.0283385,0.5577098,2.85];
ybs[1199]=['',1.0454907,1.1020215,5.03];
ybs[1200]=['',1.0439145,1.0677578,5];
ybs[1201]=['',1.022404,-0.3204959,6.22];
ybs[1202]=['',1.0373563,0.8367318,5.37];
ybs[1203]=['γ Hyi',0.9899559,-1.2944083,3.24];
ybs[1204]=['',1.033763,0.5430767,6.1];
ybs[1205]=['43 Per',1.0403833,0.886012,5.28];
ybs[1206]=['32 Eri',1.0276583,-0.0502999,6.14];
ybs[1207]=['32 Eri',1.0276655,-0.0503338,4.79];
ybs[1208]=['τ8 Eri',1.0243023,-0.4283274,4.65];
ybs[1209]=['',1.0235565,-0.6049488,5.11];
ybs[1210]=['',1.0387502,0.6134957,5.49];
ybs[1211]=['',1.0223894,-0.8172035,5.93];
ybs[1212]=['',1.0316094,-0.2099429,6];
ybs[1213]=['32 Tau',1.0398614,0.393528,5.63];
ybs[1214]=['',1.0264574,-0.7031298,5.71];
ybs[1215]=['ε Per',1.0450369,0.6995133,2.89];
ybs[1216]=['33 Tau',1.040742,0.4056999,6.06];
ybs[1217]=['',1.0424415,0.4281482,6.16];
ybs[1218]=['',1.0456066,0.6088274,6.53];
ybs[1219]=['',1.0399283,0.106629,6.09];
ybs[1220]=['',1.0376289,-0.1689689,6.19];
ybs[1221]=['',1.0477252,0.6790873,6.3];
ybs[1222]=['',1.0263028,-0.9183877,6.46];
ybs[1223]=['ξ Per',1.0496479,0.6258651,4.04];
ybs[1224]=['',1.052885,0.6787324,6.38];
ybs[1225]=['',1.1095131,1.4095358,5.1];
ybs[1226]=['γ Eri',1.0435886,-0.2345669,2.95];
ybs[1227]=['',1.0475574,-0.0942745,5.83];
ybs[1228]=['',1.0516504,0.1814942,6.37];
ybs[1229]=['',1.0596718,0.6467691,6.41];
ybs[1230]=['',1.0500436,-0.218271,5.6];
ybs[1231]=['',1.0314317,-1.1064249,6.14];
ybs[1232]=['',1.0560055,0.3030622,6.32];
ybs[1233]=['',1.0569069,0.3187199,5.89];
ybs[1234]=['λ Tau',1.0561032,0.2191745,3.47];
ybs[1235]=['τ9 Eri',1.0514312,-0.4179788,4.66];
ybs[1236]=['',1.0844677,1.1998162,5.87];
ybs[1237]=['',1.0756115,1.0335998,5.06];
ybs[1238]=['',1.0607593,0.1756624,5.67];
ybs[1239]=['35 Eri',1.0593122,-0.0258769,5.28];
ybs[1240]=['',1.043895,-0.9954271,6.05];
ybs[1241]=['',1.0544205,-0.5309856,5.93];
ybs[1242]=['δ Ret',1.0434348,-1.0704372,4.56];
ybs[1243]=['',1.0863142,1.1446737,6.17];
ybs[1244]=['',1.0640636,-0.003532,5.38];
ybs[1245]=['',1.0512198,-0.8987892,6.51];
ybs[1246]=['ν Tau',1.0666697,0.1056866,3.91];
ybs[1247]=['36 Tau',1.0726434,0.4218704,5.47];
ybs[1248]=['40 Tau',1.0692094,0.0960188,5.33];
ybs[1249]=['',1.0701842,0.1442171,5.46];
ybs[1250]=['',1.0844126,0.9437497,6.31];
ybs[1251]=['37 Tau',1.0740114,0.3865438,4.36];
ybs[1252]=['',1.070954,0.0504863,5.36];
ybs[1253]=['',1.0668005,-0.3504274,6.46];
ybs[1254]=['',1.0676795,-0.3506766,7.01];
ybs[1255]=['',1.0909588,1.0889727,6.99];
ybs[1256]=['λ Per',1.0839382,0.8799197,4.29];
ybs[1257]=['39 Tau',1.0768094,0.385263,5.9];
ybs[1258]=['',1.0701425,-0.288383,6.39];
ybs[1259]=['γ Ret',1.0526805,-1.0837064,4.51];
ybs[1260]=['',1.071308,-0.2221257,5.61];
ybs[1261]=['ι Ret',1.0546232,-1.0648509,4.97];
ybs[1262]=['',1.0723369,-0.3545842,6.13];
ybs[1263]=['41 Tau',1.0826024,0.4828344,5.2];
ybs[1264]=['ψ Tau',1.0844223,0.5072894,5.23];
ybs[1265]=['',1.0976063,1.0466874,6.28];
ybs[1266]=['',0.9531001,-1.4816651,6.41];
ybs[1267]=['',1.0782769,-0.1534372,6.26];
ybs[1268]=['48 Per',1.0927921,0.8338435,4.04];
ybs[1269]=['',1.0771074,-0.3568728,6.34];
ybs[1270]=['',1.0761181,-0.4814778,5.59];
ybs[1271]=['',1.0965779,0.9580406,6.18];
ybs[1272]=['49 Per',1.0903279,0.6595817,6.09];
ybs[1273]=['50 Per',1.0918967,0.6650229,5.51];
ybs[1274]=['',1.0868455,0.2657544,6.01];
ybs[1275]=['',1.0881981,0.3037465,5.89];
ybs[1276]=['',1.1192997,1.2598946,6.03];
ybs[1277]=['',1.1142076,1.1966408,6.32];
ybs[1278]=['ω1 Tau',1.0934271,0.3433448,5.5];
ybs[1279]=['',1.0925667,0.2349469,5.95];
ybs[1280]=['',1.0830937,-0.7479228,6.59];
ybs[1281]=['',1.1020099,0.5872808,5.72];
ybs[1282]=['44 Tau',1.1009957,0.4632623,5.41];
ybs[1283]=['',1.0926215,-0.2848861,5.37];
ybs[1284]=['',1.1961952,1.4636214,5.57];
ybs[1285]=['37 Eri',1.0976902,-0.1197544,5.44];
ybs[1286]=['',1.0878798,-0.7993815,6.59];
ybs[1287]=['45 Tau',1.1023557,0.0974764,5.72];
ybs[1288]=['',1.0994511,-0.1528468,5.7];
ybs[1289]=['',1.0804496,-1.1197718,6.38];
ybs[1290]=['',1.1079827,0.3026187,6.09];
ybs[1291]=['',1.1216149,1.0039142,6.08];
ybs[1292]=['',1.1096391,0.3922574,6.12];
ybs[1293]=['ο1 Eri',1.104194,-0.1182604,4.04];
ybs[1294]=['',1.0981445,-0.6145575,6.44];
ybs[1295]=['',1.1025197,-0.3542017,5.79];
ybs[1296]=['',1.1154013,0.6637029,6.45];
ybs[1297]=['δ Hor',1.0980854,-0.7318391,4.93];
ybs[1298]=['μ Per',1.1200951,0.8459496,4.14];
ybs[1299]=['',1.2025883,1.4554511,5.46];
ybs[1300]=['',1.1304082,1.080516,5.7];
ybs[1301]=['52 Per',1.1194665,0.7076184,4.71];
ybs[1302]=['',1.1123157,0.1792968,6.23];
ybs[1303]=['',1.1120106,0.15623,6.51];
ybs[1304]=['46 Tau',1.1120969,0.1357317,5.29];
ybs[1305]=['',1.1135126,0.2236451,6.25];
ybs[1306]=['47 Tau',1.1138504,0.1627371,4.84];
ybs[1307]=['',1.1121361,-0.0190066,6.44];
ybs[1308]=['',1.1308745,1.0108806,5.71];
ybs[1309]=['',1.1285409,0.9367328,5.19];
ybs[1310]=['',1.1167812,0.1757772,5.22];
ybs[1311]=['',1.1052576,-0.7733012,6.71];
ybs[1312]=['',1.1843267,1.4115661,5.43];
ybs[1313]=['39 Eri',1.1151015,-0.1779547,4.87];
ybs[1314]=['48 Tau',1.1220846,0.2698296,6.32];
ybs[1315]=['μ Tau',1.1207969,0.1562403,4.29];
ybs[1316]=['',1.1202362,0.1092484,6.93];
ybs[1317]=['',1.1204831,0.10902,6.31];
ybs[1318]=['',1.1101874,-0.703314,6.37];
ybs[1319]=['',1.1348748,0.8788369,4.61];
ybs[1320]=['ο2 Eri',1.1190256,-0.1325214,4.43];
ybs[1321]=['α Hor',1.1118325,-0.7371189,3.86];
ybs[1322]=['',1.1475537,1.137906,5.27];
ybs[1323]=['',1.133763,0.7365167,6.22];
ybs[1324]=['ω2 Tau',1.1287955,0.3601896,4.94];
ybs[1325]=['',1.139126,0.8745226,5.45];
ybs[1326]=['51 Tau',1.1337566,0.3776419,5.65];
ybs[1327]=['',1.1280067,-0.1119311,5.94];
ybs[1328]=['',1.1434449,0.8897317,5.55];
ybs[1329]=['',1.1333591,0.1665884,6.54];
ybs[1330]=['',1.1516448,1.0610144,5.39];
ybs[1331]=['α Ret',1.1115296,-1.0893183,3.35];
ybs[1332]=['',1.142926,0.730684,5.92];
ybs[1333]=['γ Dor',1.1199184,-0.8975716,4.25];
ybs[1334]=['53 Tau',1.1383136,0.370006,5.35];
ybs[1335]=['',1.1132618,-1.0844011,5.45];
ybs[1336]=['56 Tau',1.1391123,0.3810241,5.38];
ybs[1337]=['',1.1513561,0.9872064,5.88];
ybs[1338]=['54 Per',1.1432238,0.6042968,4.93];
ybs[1339]=['',1.1420201,0.558688,6.16];
ybs[1340]=['',1.1315599,-0.3605373,6];
ybs[1341]=['γ Tau',1.1396466,0.2737528,3.65];
ybs[1342]=['υ4 Eri',1.1293204,-0.5888707,3.56];
ybs[1343]=['φ Tau',1.1426119,0.4783577,4.95];
ybs[1344]=['',1.1386933,0.1776549,6.31];
ybs[1345]=['53 Per',1.1490076,0.8125417,4.85];
ybs[1346]=['57 Tau',1.140317,0.2459617,5.59];
ybs[1347]=['',1.1565938,1.0414704,6.19];
ybs[1348]=['',1.1330246,-0.3998929,6.07];
ybs[1349]=['',1.1425033,0.3281137,6.12];
ybs[1350]=['ε Ret',1.1209865,-1.0339771,4.44];
ybs[1351]=['58 Tau',1.1431675,0.2644561,5.26];
ybs[1352]=['',1.1201727,-1.0627153,6.37];
ybs[1353]=['',1.1443181,0.2429667,6.17];
ybs[1354]=['',1.1343509,-0.5907432,6.37];
ybs[1355]=['',1.1431727,0.1079966,5.77];
ybs[1356]=['',1.143859,0.1620084,6.53];
ybs[1357]=['',1.1425218,-0.1080112,6.27];
ybs[1358]=['',1.1427683,-0.1315204,5.85];
ybs[1359]=['',1.1346874,-0.7716141,5.34];
ybs[1360]=['',1.1312856,-0.9215652,6.09];
ybs[1361]=['',1.1462745,-0.0007249,5.86];
ybs[1362]=['',1.1419658,-0.3592362,5.38];
ybs[1363]=['60 Tau',1.1494694,0.246674,5.72];
ybs[1364]=['χ Tau',1.152251,0.4482879,5.37];
ybs[1365]=['',1.1511545,0.3643786,5.91];
ybs[1366]=['',1.1576352,0.7414733,6.23];
ybs[1367]=['θ Ret',1.1255197,-1.1029911,5.87];
ybs[1368]=['δ1 Tau',1.1534376,0.3071461,3.76];
ybs[1369]=['',1.1455399,-0.4480572,6.01];
ybs[1370]=['',1.1562249,0.3671747,5.99];
ybs[1371]=['63 Tau',1.155508,0.2937849,5.64];
ybs[1372]=['55 Per',1.1609991,0.5966467,5.73];
ybs[1373]=['62 Tau',1.1583665,0.4250957,6.36];
ybs[1374]=['56 Per',1.1615863,0.5936638,5.76];
ybs[1375]=['δ2 Tau',1.1585037,0.3054139,4.8];
ybs[1376]=['66 Tau',1.1571712,0.1660859,5.12];
ybs[1377]=['',1.1740066,1.0059825,6.32];
ybs[1378]=['ξ Eri',1.1558586,-0.0644068,5.17];
ybs[1379]=['',1.1524404,-0.4334792,5.83];
ybs[1380]=['',1.1623029,0.3332922,5.98];
ybs[1381]=['',1.152042,-0.619405,6.39];
ybs[1382]=['κ1 Tau',1.1642633,0.3900501,4.22];
ybs[1383]=['κ2 Tau',1.1644701,0.3884061,5.28];
ybs[1384]=['δ3 Tau',1.1646055,0.313851,4.29];
ybs[1385]=['',1.1679194,0.5496529,5.28];
ybs[1386]=['70 Tau',1.1650973,0.2791662,6.46];
ybs[1387]=['υ Tau',1.1683897,0.3991119,4.28];
ybs[1388]=['43 Eri',1.1560886,-0.5927441,3.96];
ybs[1389]=['71 Tau',1.1682482,0.2735306,4.49];
ybs[1390]=['η Ret',1.1438643,-1.1053142,5.24];
ybs[1391]=['π Tau',1.1693458,0.2577379,4.69];
ybs[1392]=['',1.1679889,0.1508682,6.06];
ybs[1393]=['',1.1599862,-0.6056827,6.55];
ybs[1394]=['72 Tau',1.1726902,0.4022926,5.53];
ybs[1395]=['',1.1706178,0.0372267,6.23];
ybs[1396]=['',1.205958,1.2667264,5.94];
ybs[1397]=['',1.1730163,0.1966237,5.88];
ybs[1398]=['',1.1757787,0.3782633,5.72];
ybs[1399]=['',1.1610317,-0.7697998,6.39];
ybs[1400]=['',1.1549262,-0.9951196,6.29];
ybs[1401]=['',1.1799128,0.5308207,6.4];
ybs[1402]=['75 Tau',1.1774153,0.2864502,4.97];
ybs[1403]=['76 Tau',1.1771301,0.2581959,5.9];
ybs[1404]=['ε Tau',1.1783056,0.3356764,3.53];
ybs[1405]=['',1.169291,-0.4193642,6.11];
ybs[1406]=['θ1 Tau',1.1779879,0.2795113,3.84];
ybs[1407]=['θ2 Tau',1.1783624,0.2779154,3.4];
ybs[1408]=['',1.1751836,0.0333624,6.15];
ybs[1409]=['79 Tau',1.1790105,0.2286374,5.03];
ybs[1410]=['',1.177238,0.0250191,5.55];
ybs[1411]=['',1.1581778,-1.0678542,5.94];
ybs[1412]=['1 Cam',1.1954701,0.9418036,5.77];
ybs[1413]=['',1.1685816,-0.8184532,6.1];
ybs[1414]=['',1.1877654,0.5673972,6.21];
ybs[1415]=['',1.1827425,0.184545,6.79];
ybs[1416]=['',1.1769111,-0.3386932,5.96];
ybs[1417]=['80 Tau',1.1848194,0.2738385,5.58];
ybs[1418]=['',1.1792064,-0.2268224,5.6];
ybs[1419]=['',1.1915402,0.6992006,6.26];
ybs[1420]=['',1.1841506,0.1800186,6.48];
ybs[1421]=['δ Men',1.1187526,-1.398964,5.69];
ybs[1422]=['',1.1866683,0.2835356,4.78];
ybs[1423]=['81 Tau',1.1870256,0.2747742,5.48];
ybs[1424]=['',1.1735627,-0.7314154,6.44];
ybs[1425]=['83 Tau',1.1868283,0.2404351,5.4];
ybs[1426]=['',1.1837502,-0.236325,6.24];
ybs[1427]=['85 Tau',1.192335,0.2775502,6.02];
ybs[1428]=['',1.1783519,-0.8109308,6.16];
ybs[1429]=['57 Per',1.2005499,0.7524764,6.09];
ybs[1430]=['',1.1696165,-1.0902634,5.75];
ybs[1431]=['',1.1928623,0.095307,6.39];
ybs[1432]=['45 Eri',1.1917713,0.0001209,4.91];
ybs[1433]=['',1.1892884,-0.2372584,6.21];
ybs[1434]=['',1.1849276,-0.6213726,5.96];
ybs[1435]=['',1.2161431,1.1224157,5.94];
ybs[1436]=['',1.1948997,-0.0551306,5.81];
ybs[1437]=['',1.1997806,0.3153202,6.25];
ybs[1438]=['δ Cae',1.185021,-0.783694,5.07];
ybs[1439]=['ρ Tau',1.2009565,0.2599516,4.65];
ybs[1440]=['',1.205023,0.5063256,5.88];
ybs[1441]=['',1.2005326,0.1651567,6.01];
ybs[1442]=['',1.1978452,-0.1873753,6.06];
ybs[1443]=['',1.2018503,0.0980554,5.68];
ybs[1444]=['46 Eri',1.2003905,-0.1167483,5.72];
ybs[1445]=['',1.2017972,-0.1184773,6.09];
ybs[1446]=['47 Eri',1.2015533,-0.1427999,5.11];
ybs[1447]=['',1.2015313,-0.155696,5.26];
ybs[1448]=['υ1 Eri',1.197626,-0.5186537,4.51];
ybs[1449]=['58 Per',1.2147236,0.7210435,4.25];
ybs[1450]=['',1.2093035,0.3478494,6.36];
ybs[1451]=['ν Men',1.1296305,-1.4228328,5.79];
ybs[1452]=['α Tau',1.2100676,0.2889864,0.85];
ybs[1453]=['88 Tau',1.2086503,0.1781899,4.25];
ybs[1454]=['',1.2128418,0.4082156,6.02];
ybs[1455]=['',1.2060308,-0.169082,6.37];
ybs[1456]=['',1.2046255,-0.3468218,6.13];
ybs[1457]=['',1.209733,-0.0621934,6.33];
ybs[1458]=['ν Eri',1.2110158,-0.0576681,3.93];
ybs[1459]=['υ2 Eri',1.2064816,-0.5325585,3.82];
ybs[1460]=['α Dor',1.197848,-0.9598465,3.27];
ybs[1461]=['2 Cam',1.2300986,0.9340859,5.35];
ybs[1462]=['3 Cam',1.2298079,0.9272215,5.05];
ybs[1463]=['',1.2629349,1.3378545,6.49];
ybs[1464]=['',1.2151556,0.0182593,5.31];
ybs[1465]=['',1.2217721,0.4710132,6.47];
ybs[1466]=['',1.2204761,0.3618406,5.92];
ybs[1467]=['89 Tau',1.2198091,0.2806598,5.79];
ybs[1468]=['90 Tau',1.2196702,0.2191807,4.27];
ybs[1469]=['51 Eri',1.2166494,-0.0423362,5.23];
ybs[1470]=['',1.1948315,-1.0956035,5.79];
ybs[1471]=['',1.2121356,-0.5352667,6.3];
ybs[1472]=['',1.2255809,0.4409558,6.22];
ybs[1473]=['σ1 Tau',1.2241491,0.2765728,5.07];
ybs[1474]=['σ2 Tau',1.2246851,0.2786369,4.69];
ybs[1475]=['',1.2236071,0.1381885,5.39];
ybs[1476]=['53 Eri',1.2187,-0.2488235,3.87];
ybs[1477]=['',1.2358474,0.8438001,5.67];
ybs[1478]=['',1.2219016,-0.2107678,5.01];
ybs[1479]=['93 Tau',1.2279418,0.2136984,5.46];
ybs[1480]=['',1.135069,-1.4458714,6.76];
ybs[1481]=['',1.2454456,1.0396059,6.5];
ybs[1482]=['',1.2237072,-0.2497994,5.45];
ybs[1483]=['',1.2262384,-0.0175642,6.1];
ybs[1484]=['',1.237005,0.668905,5.99];
ybs[1485]=['',1.2342383,0.5002197,5.78];
ybs[1486]=['',1.2748604,1.3261332,6.06];
ybs[1487]=['',1.2088767,-1.0826125,5.4];
ybs[1488]=['',1.2445535,0.872981,5.87];
ybs[1489]=['59 Per',1.2420118,0.7576392,5.29];
ybs[1490]=['',1.2266855,-0.4264921,5.58];
ybs[1491]=['54 Eri',1.2283367,-0.3425307,4.32];
ybs[1492]=['τ Tau',1.2379582,0.4014593,4.28];
ybs[1493]=['',1.2203343,-0.9010398,6.44];
ybs[1494]=['95 Tau',1.2423104,0.4212057,6.13];
ybs[1495]=['',1.2475475,0.7126311,6.08];
ybs[1496]=['',1.2452687,0.5743764,6.45];
ybs[1497]=['α Cae',1.2276542,-0.7298579,4.45];
ybs[1498]=['β Cae',1.2344909,-0.6475031,5.05];
ybs[1499]=['',1.224826,-1.0279508,6.53];
ybs[1500]=['55 Eri',1.2424682,-0.1527041,6.82];
ybs[1501]=['55 Eri',1.2425045,-0.1527478,6.7];
ybs[1502]=['',1.2469841,0.1953004,5.4];
ybs[1503]=['56 Eri',1.2447122,-0.1476476,5.9];
ybs[1504]=['',1.2396378,-0.5361815,5.68];
ybs[1505]=['',1.280087,1.2388618,6.37];
ybs[1506]=['4 Cam',1.2655999,0.9913259,5.3];
ybs[1507]=['',1.2531061,0.4131383,6.35];
ybs[1508]=['',1.2444843,-0.3250264,5.53];
ybs[1509]=['',1.2585396,0.7043304,5.97];
ybs[1510]=['',1.2659398,0.9711712,6.26];
ybs[1511]=['λ Pic',1.2365876,-0.8802824,5.31];
ybs[1512]=['',1.2553789,0.327733,6.01];
ybs[1513]=['',1.2415538,-0.7159409,6.25];
ybs[1514]=['',1.2539821,0.2050489,5.37];
ybs[1515]=['μ Eri',1.2510875,-0.0560514,4.02];
ybs[1516]=['',1.2484536,-0.3707098,5.72];
ybs[1517]=['',1.2550262,-0.0508194,6.33];
ybs[1518]=['',1.3308992,1.4176899,5.07];
ybs[1519]=['',1.2511107,-0.5927466,6.86];
ybs[1520]=['',1.2540603,-0.4894731,6.19];
ybs[1521]=['',1.2512039,-0.6861514,6.05];
ybs[1522]=['',1.2845261,1.1090599,5.44];
ybs[1523]=['',1.2693138,0.5694885,5.86];
ybs[1524]=['',1.2687992,0.5493989,5.58];
ybs[1525]=['κ Dor',1.2423405,-1.0417638,5.27];
ybs[1526]=['',1.2098948,-1.3545194,6.05];
ybs[1527]=['58 Eri',1.2597039,-0.2948277,5.51];
ybs[1528]=['',1.2721945,0.6550033,4.88];
ybs[1529]=['',1.2655026,0.0633501,6.03];
ybs[1530]=['',1.2784795,0.8513832,5.66];
ybs[1531]=['',1.2645314,-0.0983043,5.78];
ybs[1532]=['96 Tau',1.2703346,0.2782916,6.08];
ybs[1533]=['59 Eri',1.2638161,-0.2842778,5.77];
ybs[1534]=['ζ Cae',1.2600521,-0.5232201,6.37];
ybs[1535]=['',1.2444087,-1.1028024,6.46];
ybs[1536]=['μ Men',1.2340901,-1.237196,5.54];
ybs[1537]=['α Cam',1.2937378,1.1585633,4.29];
ybs[1538]=['π3 Ori',1.2704229,0.1222101,3.19];
ybs[1539]=['π2 Ori',1.2738686,0.1560424,4.36];
ybs[1540]=['',1.2689926,-0.2396138,6.26];
ybs[1541]=['',1.2876597,0.9229213,6.41];
ybs[1542]=['97 Tau',1.2776201,0.3295103,5.1];
ybs[1543]=['',1.262388,-0.7668695,6.72];
ybs[1544]=['60 Eri',1.2710262,-0.2823354,5.03];
ybs[1545]=['',1.285144,0.7439571,5.71];
ybs[1546]=['2 Aur',1.2840344,0.6412707,4.78];
ybs[1547]=['π4 Ori',1.276332,0.0985231,3.69];
ybs[1548]=['',1.2787637,0.1747887,6.11];
ybs[1549]=['',1.2842033,0.4875838,5.97];
ybs[1550]=['5 Cam',1.2961872,0.96511,5.52];
ybs[1551]=['ο1 Ori',1.2824765,0.2494029,4.74];
ybs[1552]=['',1.2700244,-0.7204754,6.07];
ybs[1553]=['',1.2942419,0.7696654,6.08];
ybs[1554]=['',1.2756697,-0.6085346,5.86];
ybs[1555]=['ω Eri',1.2832576,-0.0944878,4.39];
ybs[1556]=['',1.3005676,0.9233912,5.75];
ybs[1557]=['5 Ori',1.2856915,0.0444497,5.33];
ybs[1558]=['ι Pic',1.2717903,-0.9323733,5.61];
ybs[1559]=['ι Pic',1.2718704,-0.9323444,6.42];
ybs[1560]=['',1.288054,0.0280625,6.61];
ybs[1561]=['',1.293346,0.3407415,6.37];
ybs[1562]=['π5 Ori',1.2894923,0.0432629,3.72];
ybs[1563]=['7 Cam',1.3057688,0.9387867,4.47];
ybs[1564]=['6 Ori',1.2921697,0.2000848,5.19];
ybs[1565]=['π1 Ori',1.2926187,0.177826,4.65];
ybs[1566]=['',1.2920853,0.1364337,5.33];
ybs[1567]=['',1.3328732,1.2968184,6.06];
ybs[1568]=['',1.3001463,0.6319109,6.07];
ybs[1569]=['',1.2920024,0.0088208,5.99];
ybs[1570]=['',1.2992067,0.4298613,6.37];
ybs[1571]=['',1.2969333,0.2631483,5.81];
ybs[1572]=['ι Aur',1.3028568,0.5794961,2.69];
ybs[1573]=['',1.2971235,0.0948831,6.5];
ybs[1574]=['',1.2924666,-0.2915183,5.7];
ybs[1575]=['ο2 Ori',1.2991958,0.2365172,4.07];
ybs[1576]=['',1.2933387,-0.2858867,5.72];
ybs[1577]=['62 Eri',1.2985759,-0.0896117,5.51];
ybs[1578]=['',1.2937538,-0.4483783,6.72];
ybs[1579]=['',1.2903891,-0.6909868,6.1];
ybs[1580]=['',1.3037167,0.3000224,5.48];
ybs[1581]=['99 Tau',1.305948,0.4186127,5.79];
ybs[1582]=['',1.3408298,1.2879809,6.66];
ybs[1583]=['8 Cam',1.3165531,0.9283481,6.08];
ybs[1584]=['',1.3429425,1.2932655,5.96];
ybs[1585]=['98 Tau',1.3075071,0.4378369,5.81];
ybs[1586]=['',1.3025948,-0.0179894,6.23];
ybs[1587]=['ω Aur',1.3130233,0.6619253,4.94];
ybs[1588]=['',1.3256275,1.0666022,6.03];
ybs[1589]=['',1.3322679,1.1668525,6.19];
ybs[1590]=['',1.3040627,-0.2477515,6.15];
ybs[1591]=['',1.3064467,-0.0379872,6.35];
ybs[1592]=['',1.2884163,-1.0211765,6.12];
ybs[1593]=['',1.2808651,-1.1630272,6.41];
ybs[1594]=['5 Aur',1.317699,0.6881722,5.95];
ybs[1595]=['',1.3106661,0.2544387,6.09];
ybs[1596]=['π6 Ori',1.3082113,0.0305422,4.47];
ybs[1597]=['6 Aur',1.318073,0.6927141,6.58];
ybs[1598]=['β Cam',1.3334693,1.0554865,4.03];
ybs[1599]=['',1.3095326,-0.2851918,5.66];
ybs[1600]=['ε Aur',1.3252798,0.765449,2.99];
ybs[1601]=['',1.2771528,-1.2630622,6.28];
ybs[1602]=['',1.3121599,-0.257796,7.71];
ybs[1603]=['63 Eri',1.313362,-0.1785171,5.38];
ybs[1604]=['',1.3169929,0.0637027,7.03];
ybs[1605]=['',1.3170875,0.063717,6.66];
ybs[1606]=['64 Eri',1.3136596,-0.2182096,4.79];
ybs[1607]=['ζ Aur',1.3273028,0.7174911,3.75];
ybs[1608]=['',1.317288,-0.0354474,6.32];
ybs[1609]=['',1.317808,-0.0998127,6.22];
ybs[1610]=['',1.3309595,0.7238677,6.14];
ybs[1611]=['',1.4866372,1.4963906,6.51];
ybs[1612]=['ψ Eri',1.3204625,-0.1246123,4.81];
ybs[1613]=['',1.322527,0.0131965,5.92];
ybs[1614]=['',1.3232681,0.0286701,6.24];
ybs[1615]=['ι Tau',1.3288951,0.3773943,4.64];
ybs[1616]=['',1.3198605,-0.349376,4.91];
ybs[1617]=['11 Cam',1.3451295,1.0298078,5.08];
ybs[1618]=['12 Cam',1.3454063,1.0306556,6.08];
ybs[1619]=['',1.3470308,1.0681574,6.04];
ybs[1620]=['',1.3263422,-0.0728912,5.85];
ybs[1621]=['',1.3343464,0.5327989,6.14];
ybs[1622]=['',1.3360772,0.564657,6.62];
ybs[1623]=['',1.3227775,-0.4579959,5.02];
ybs[1624]=['η Men',1.2849947,-1.3072287,5.47];
ybs[1625]=['',1.3454962,0.9501028,7.24];
ybs[1626]=['',1.3194348,-0.6926146,6.03];
ybs[1627]=['',1.3359006,0.4839501,6.6];
ybs[1628]=['',1.3344014,0.371937,6.19];
ybs[1629]=['1 Lep',1.3254987,-0.3972645,5.75];
ybs[1630]=['',1.3234328,-0.5539277,5.94];
ybs[1631]=['',1.362688,1.2159431,6.41];
ybs[1632]=['9 Aur',1.3465353,0.9010903,5];
ybs[1633]=['11 Ori',1.3350471,0.2694165,4.68];
ybs[1634]=['',1.3424016,0.6277556,6.52];
ybs[1635]=['',1.3307562,-0.2502225,6.41];
ybs[1636]=['η Aur',1.3449359,0.7202186,3.17];
ybs[1637]=['',1.3394616,0.3462398,6.44];
ybs[1638]=['',1.3765785,1.2910881,5.43];
ybs[1639]=['',1.3464397,0.7540796,6.2];
ybs[1640]=['',1.3303952,-0.4250751,5.61];
ybs[1641]=['',1.3357704,-0.0524926,6.05];
ybs[1642]=['',1.3619557,1.1335638,6.41];
ybs[1643]=['',1.3380644,0.0211067,6.17];
ybs[1644]=['η1 Pic',1.3240696,-0.8572689,5.38];
ybs[1645]=['',1.3878772,1.3351526,6.37];
ybs[1646]=['',1.3294322,-0.7280147,6.31];
ybs[1647]=['γ1 Cae',1.3320565,-0.6187336,4.55];
ybs[1648]=['γ2 Cae',1.3321666,-0.6226076,6.34];
ybs[1649]=['ε Lep',1.3373551,-0.3898936,3.19];
ybs[1650]=['',1.3363352,-0.4558891,5.73];
ybs[1651]=['104 Tau',1.3477662,0.3259504,5];
ybs[1652]=['66 Eri',1.3437912,-0.0807032,5.12];
ybs[1653]=['106 Tau',1.3494036,0.3568972,5.3];
ybs[1654]=['103 Tau',1.3509098,0.4240358,5.5];
ybs[1655]=['105 Tau',1.3499799,0.3793477,5.89];
ybs[1656]=['',1.3427852,-0.2284773,6.05];
ybs[1657]=['13 Ori',1.3481923,0.165849,6.17];
ybs[1658]=['η2 Pic',1.333433,-0.8647325,5.03];
ybs[1659]=['14 Ori',1.3492136,0.1488539,5.34];
ybs[1660]=['',1.3463242,-0.2174656,5.97];
ybs[1661]=['β Eri',1.3485221,-0.0882433,2.79];
ybs[1662]=['',1.3331065,-0.9490269,6.27];
ybs[1663]=['',1.3637231,0.8201441,5.68];
ybs[1664]=['',1.3612572,0.6515453,6.02];
ybs[1665]=['',1.3582644,0.4897355,6.01];
ybs[1666]=['',1.3504983,-0.1507066,5.78];
ybs[1667]=['16 Ori',1.3555745,0.1720714,5.43];
ybs[1668]=['68 Eri',1.3523728,-0.0772518,5.12];
ybs[1669]=['ζ Dor',1.3348914,-1.0025307,4.72];
ybs[1670]=['',1.3757949,1.0799589,6.17];
ybs[1671]=['15 Ori',1.3574499,0.2727341,4.82];
ybs[1672]=['β Men',1.3195172,-1.2440815,5.31];
ybs[1673]=['14 Cam',1.3779986,1.0946389,6.5];
ybs[1674]=['λ Eri',1.3540286,-0.152271,4.27];
ybs[1675]=['',1.3487907,-0.6228739,6.52];
ybs[1676]=['',1.3583196,-0.0093576,6.1];
ybs[1677]=['',1.3043547,-1.3659773,6.29];
ybs[1678]=['',1.4017859,1.2791848,5.74];
ybs[1679]=['',1.3661683,0.2805387,5.18];
ybs[1680]=['',1.362236,-0.0388386,6.25];
ybs[1681]=['',1.4252851,1.38321,5.05];
ybs[1682]=['',1.3637681,-0.0429776,5.9];
ybs[1683]=['',1.3847027,1.037274,6.15];
ybs[1684]=['μ Aur',1.374924,0.6721522,4.86];
ybs[1685]=['',1.3655047,0.0094753,6.67];
ybs[1686]=['',1.3658096,0.0185891,5.89];
ybs[1687]=['',1.3819101,0.9292146,6.2];
ybs[1688]=['',1.3636427,-0.2063116,5.68];
ybs[1689]=['',1.3602157,-0.4517024,6.41];
ybs[1690]=['',1.3428645,-1.105995,5.2];
ybs[1691]=['ι Lep',1.3676483,-0.2066699,4.45];
ybs[1692]=['',1.3700942,-0.1052378,5.91];
ybs[1693]=['ρ Ori',1.372589,0.0504111,4.46];
ybs[1694]=['',1.3633041,-0.652176,6.57];
ybs[1695]=['',1.3337512,-1.2741931,6.27];
ybs[1696]=['',1.373578,0.0348219,6.09];
ybs[1697]=['μ Lep',1.370224,-0.2823604,3.31];
ybs[1698]=['',1.374655,0.0102491,6.32];
ybs[1699]=['',1.3732874,-0.1417324,6.37];
ybs[1700]=['κ Lep',1.3716747,-0.2253931,4.36];
ybs[1701]=['14 Aur',1.3832008,0.5709614,5.02];
ybs[1702]=['',1.3931224,0.935685,6.5];
ybs[1703]=['α Aur',1.3897238,0.8032555,0.08];
ybs[1704]=['',1.3789808,0.0904516,5.5];
ybs[1705]=['',1.3749476,-0.2544655,6.21];
ybs[1706]=['108 Tau',1.3828987,0.389394,6.27];
ybs[1707]=['',1.3872098,0.5992993,5.96];
ybs[1708]=['β Ori',1.3775754,-0.142683,0.12];
ybs[1709]=['',1.5370512,1.4910419,6.6];
ybs[1710]=['',1.3728972,-0.6248014,6.98];
ybs[1711]=['ξ Men',1.2921483,-1.4387414,5.85];
ybs[1712]=['',1.381207,-0.0241396,6.15];
ybs[1713]=['18 Ori',1.3850569,0.1983912,5.56];
ybs[1714]=['15 Cam',1.403173,1.0147449,6.13];
ybs[1715]=['',1.4079266,1.0939096,5.61];
ybs[1716]=['',1.3759629,-0.6274557,5.76];
ybs[1717]=['',1.3963419,0.7472874,5.48];
ybs[1718]=['',1.3805059,-0.4697946,5.07];
ybs[1719]=['',1.3873581,0.0344261,6.42];
ybs[1720]=['',1.3979537,0.7066657,6.18];
ybs[1721]=['16 Aur',1.3953376,0.5828693,4.54];
ybs[1722]=['',1.372087,-0.9076414,6.05];
ybs[1723]=['17 Aur',1.3959574,0.5897716,6.14];
ybs[1724]=['λ Aur',1.3999722,0.7002759,4.71];
ybs[1725]=['',1.3817124,-0.6091333,6.66];
ybs[1726]=['',1.387066,-0.298738,6.56];
ybs[1727]=['',1.398946,0.5894349,5.41];
ybs[1728]=['',1.4042102,0.7757761,6.62];
ybs[1729]=['18 Aur',1.4006767,0.5935711,6.49];
ybs[1730]=['τ Ori',1.391019,-0.1190263,3.6];
ybs[1731]=['',1.4071153,0.8200718,6.54];
ybs[1732]=['',1.3910136,-0.2355323,5.5];
ybs[1733]=['',1.4048654,0.7174899,5.52];
ybs[1734]=['109 Tau',1.3995387,0.3860679,4.94];
ybs[1735]=['19 Aur',1.4033884,0.5930849,5.03];
ybs[1736]=['',1.3993082,0.3518309,6.08];
ybs[1737]=['',1.3797456,-0.9102918,6.49];
ybs[1738]=['ο Col',1.3891369,-0.6086029,4.83];
ybs[1739]=['θ Dor',1.3689602,-1.1721271,4.83];
ybs[1740]=['',1.4539803,1.3612612,6.56];
ybs[1741]=['21 Ori',1.3983004,0.045721,5.34];
ybs[1742]=['',1.3959142,-0.3160085,5.96];
ybs[1743]=['',1.3998816,-0.0242366,6.34];
ybs[1744]=['ρ Aur',1.411727,0.7300112,5.23];
ybs[1745]=['',1.4073014,0.488341,6.33];
ybs[1746]=['16 Cam',1.4205539,1.0047071,5.28];
ybs[1747]=['',1.4083621,0.5164869,5.76];
ybs[1748]=['',1.397867,-0.3228199,6.36];
ybs[1749]=['',1.3979257,-0.3226406,6.54];
ybs[1750]=['',1.4067048,0.3462187,6.18];
ybs[1751]=['λ Lep',1.3993396,-0.2295642,4.29];
ybs[1752]=['ν Lep',1.4011655,-0.2145393,5.3];
ybs[1753]=['',1.3978897,-0.4772626,5.99];
ybs[1754]=['',1.4034424,-0.0932633,6.39];
ybs[1755]=['',1.4161734,0.7164744,5.54];
ybs[1756]=['',1.4076745,0.0704148,6.57];
ybs[1757]=['',1.4027873,-0.3702948,4.71];
ybs[1758]=['',1.4096249,0.1474959,5.8];
ybs[1759]=['',1.4084017,-0.0068761,5.68];
ybs[1760]=['22 Ori',1.4094139,-0.006287,4.73];
ybs[1761]=['',1.401617,-0.6052045,6.34];
ybs[1762]=['ζ Pic',1.3961334,-0.8828265,5.45];
ybs[1763]=['22 Aur',1.4177984,0.5054109,6.46];
ybs[1764]=['',1.40924,-0.2397007,6.56];
ybs[1765]=['23 Ori',1.4142513,0.06224,5];
ybs[1766]=['',1.4083823,-0.4319813,5.06];
ybs[1767]=['',1.4057314,-0.5990425,6.09];
ybs[1768]=['σ Aur',1.4238457,0.6528581,4.99];
ybs[1769]=['110 Tau',1.4182801,0.2918292,6.08];
ybs[1770]=['',1.4234082,0.5453174,6.28];
ybs[1771]=['',1.4234183,0.5439066,5.94];
ybs[1772]=['',1.4173135,0.0932657,6.35];
ybs[1773]=['',1.4158303,-0.1465105,5.9];
ybs[1774]=['',1.4261448,0.6086909,6.55];
ybs[1775]=['111 Tau',1.42178,0.3037571,4.99];
ybs[1776]=['',1.4179024,-0.0024187,5.7];
ybs[1777]=['',1.4185352,-0.0147684,6.11];
ybs[1778]=['8 Lep',1.416446,-0.2427042,5.25];
ybs[1779]=['29 Ori',1.4186409,-0.1359094,4.14];
ybs[1780]=['',1.414523,-0.4657291,6.49];
ybs[1781]=['',1.4219253,0.0414234,6.32];
ybs[1782]=['27 Ori',1.4212612,-0.0151966,5.08];
ybs[1783]=['η Ori',1.4211775,-0.0414734,3.36];
ybs[1784]=['ψ1 Ori',1.4225299,0.0325838,4.95];
ybs[1785]=['γ Ori',1.4243959,0.1111776,1.64];
ybs[1786]=['β Tau',1.4304805,0.499636,1.65];
ybs[1787]=['',1.4205432,-0.2959271,5.65];
ybs[1788]=['',1.4146282,-0.692148,5.71];
ybs[1789]=['',1.4335503,0.6191794,6.15];
ybs[1790]=['',1.4330923,0.6005829,5.94];
ybs[1791]=['',1.4331996,0.5808798,6.15];
ybs[1792]=['',1.4158756,-0.6512761,6.82];
ybs[1793]=['113 Tau',1.4290441,0.2918184,6.25];
ybs[1794]=['',1.4232469,-0.1799221,5.61];
ybs[1795]=['',1.425806,-0.0091471,6.57];
ybs[1796]=['κ Pic',1.4085775,-0.9793441,6.11];
ybs[1797]=['17 Cam',1.4507683,1.1010277,5.42];
ybs[1798]=['',1.4269988,0.0094379,6.16];
ybs[1799]=['',1.4342623,0.5275718,5.74];
ybs[1800]=['φ Aur',1.4367344,0.6020434,5.07];
ybs[1801]=['',1.4278698,-0.0959678,6.23];
ybs[1802]=['',1.4310283,0.120228,6.42];
ybs[1803]=['115 Tau',1.4337849,0.3138324,5.42];
ybs[1804]=['',1.4339327,0.2666305,6.16];
ybs[1805]=['114 Tau',1.4360065,0.3831994,4.88];
ybs[1806]=['ψ2 Ori',1.4317009,0.0543643,4.59];
ybs[1807]=['',1.4270652,-0.3434059,5.89];
ybs[1808]=['',1.4209493,-0.7715273,6.08];
ybs[1809]=['116 Tau',1.4362728,0.277383,5.5];
ybs[1810]=['',1.3530837,-1.4226682,6.51];
ybs[1811]=['117 Tau',1.4374976,0.3011991,5.77];
ybs[1812]=['',1.4321384,-0.2073738,6.35];
ybs[1813]=['θ Pic',1.4195105,-0.9127316,6.27];
ybs[1814]=['',1.4397543,0.2390598,6.35];
ybs[1815]=['',1.4368192,0.0229847,6.41];
ybs[1816]=['118 Tau',1.4433208,0.4392704,5.47];
ybs[1817]=['',1.445284,0.5097046,6.24];
ybs[1818]=['',1.4340177,-0.3727438,6.07];
ybs[1819]=['',1.4509979,0.7239408,6];
ybs[1820]=['',1.450626,0.6953861,6.37];
ybs[1821]=['',1.4406363,-0.0574115,6.39];
ybs[1822]=['',1.4306273,-0.7142642,5.87];
ybs[1823]=['18 Cam',1.4602339,0.9989703,6.48];
ybs[1824]=['β Lep',1.4368245,-0.3619973,2.84];
ybs[1825]=['',1.4425869,-0.0598401,5.79];
ybs[1826]=['',1.4495089,0.3923403,6.29];
ybs[1827]=['',1.4479229,0.2683864,5.94];
ybs[1828]=['',1.4450716,0.0315319,5.78];
ybs[1829]=['31 Ori',1.4441675,-0.0187558,4.71];
ybs[1830]=['',1.4359701,-0.6494715,5.57];
ybs[1831]=['λ Dor',1.4254213,-1.0278704,5.14];
ybs[1832]=['',1.4469895,0.0736773,6.21];
ybs[1833]=['',1.4401313,-0.52532,6.75];
ybs[1834]=['32 Ori',1.4490473,0.1041091,4.2];
ybs[1835]=['',1.4465732,-0.1294593,6.33];
ybs[1836]=['33 Ori',1.4509297,0.0577514,5.46];
ybs[1837]=['χ Aur',1.4587694,0.5621298,4.76];
ybs[1838]=['',1.4967093,1.3099562,6.17];
ybs[1839]=['119 Tau',1.4558288,0.3248148,4.38];
ybs[1840]=['',1.4626915,0.7352047,6.55];
ybs[1841]=['',1.4558607,0.2979995,5.46];
ybs[1842]=['',1.4509813,-0.1167918,6.22];
ybs[1843]=['10 Lep',1.4493916,-0.363845,5.55];
ybs[1844]=['',1.4619906,0.5727542,6.48];
ybs[1845]=['δ Ori',1.4541272,-0.0046809,6.85];
ybs[1846]=['δ Ori',1.4541194,-0.0049378,2.23];
ybs[1847]=['',1.4826011,1.1642894,6.26];
ybs[1848]=['',1.4628692,0.606345,6.27];
ybs[1849]=['υ Ori',1.4534943,-0.1271486,4.62];
ybs[1850]=['',1.4435287,-0.8213566,5.46];
ybs[1851]=['19 Cam',1.4819113,1.1199347,6.15];
ybs[1852]=['120 Tau',1.4615573,0.3238559,5.69];
ybs[1853]=['',1.4262703,-1.1973512,6.03];
ybs[1854]=['',1.4621709,0.3576072,6.18];
ybs[1855]=['',1.4570403,-0.0275081,5.35];
ybs[1856]=['ε Col',1.4489798,-0.6187842,3.87];
ybs[1857]=['',1.4589186,-0.0297184,6.46];
ybs[1858]=['35 Ori',1.4630181,0.2499423,5.64];
ybs[1859]=['α Lep',1.4565278,-0.3107795,2.58];
ybs[1860]=['',1.4774291,0.9501957,5.73];
ybs[1861]=['',1.4377827,-1.0872753,6.59];
ybs[1862]=['',1.4607017,-0.0199099,5.34];
ybs[1863]=['',1.4753172,0.8330252,6.11];
ybs[1864]=['',1.4499246,-0.8012567,5.86];
ybs[1865]=['',1.4627129,0.0248338,6.59];
ybs[1866]=['38 Ori',1.4641995,0.0660056,5.36];
ybs[1867]=['',1.4630701,-0.0178115,6.22];
ybs[1868]=['',1.4630594,-0.0254036,5.93];
ybs[1869]=['121 Tau',1.4702195,0.4198143,5.38];
ybs[1870]=['φ1 Ori',1.4668029,0.1658761,4.41];
ybs[1871]=['',1.4559519,-0.6719025,5.48];
ybs[1872]=['',1.4724701,0.4830386,6.27];
ybs[1873]=['λ Ori',1.4682108,0.1736347,3.54];
ybs[1874]=['λ Ori',1.4682254,0.1736492,5.61];
ybs[1875]=['',1.4573398,-0.6130299,5.78];
ybs[1876]=['',1.4417091,-1.115443,6.19];
ybs[1877]=['',1.4685947,0.1789716,5.6];
ybs[1878]=['',1.4773677,0.7015431,6.09];
ybs[1879]=['',1.6103367,1.4814036,6.11];
ybs[1880]=['',1.4670089,-0.1046269,5.67];
ybs[1881]=['',1.4671401,-0.1045012,4.78];
ybs[1882]=['',1.4608973,-0.5206954,6.53];
ybs[1883]=['',1.4749045,0.4529645,6.49];
ybs[1884]=['',1.4685919,-0.0781744,6.56];
ybs[1885]=['',1.4686456,-0.0769867,6.24];
ybs[1886]=['42 Ori',1.4686793,-0.0841959,4.59];
ybs[1887]=['θ1 Ori',1.4681256,-0.0937746,6.73];
ybs[1888]=['θ1 Ori',1.4681403,-0.0937406,7.96];
ybs[1889]=['θ1 Ori',1.4681692,-0.0938183,5.13];
ybs[1890]=['θ1 Ori',1.4682274,-0.0937845,6.7];
ybs[1891]=['θ2 Ori',1.4686335,-0.09428,5.08];
ybs[1892]=['',1.4692735,-0.0759313,6.38];
ybs[1893]=['ι Ori',1.4688382,-0.1029005,2.77];
ybs[1894]=['',1.4696615,-0.056525,6.4];
ybs[1895]=['45 Ori',1.469864,-0.0845042,5.26];
ybs[1896]=['',1.4777626,0.4701445,5.83];
ybs[1897]=['ε Ori',1.472437,-0.0207377,1.7];
ybs[1898]=['',1.480815,0.5859403,6.33];
ybs[1899]=['122 Tau',1.4769217,0.2976395,5.54];
ybs[1900]=['',1.4724125,-0.0983372,6.54];
ybs[1901]=['φ2 Ori',1.4758923,0.1623833,4.09];
ybs[1902]=['',1.4767013,0.1928278,5.94];
ybs[1903]=['',1.4667707,-0.577103,5.78];
ybs[1904]=['ζ Tau',1.4796501,0.3692305,3];
ybs[1905]=['',1.4739003,-0.1056177,5.72];
ybs[1906]=['',1.4583298,-0.9579548,6.43];
ybs[1907]=['',1.477696,0.1564691,6.12];
ybs[1908]=['26 Aur',1.4844495,0.532408,5.4];
ybs[1909]=['',1.471008,-0.5008031,6.26];
ybs[1910]=['',1.5049151,1.146812,5.6];
ybs[1911]=['',1.4535418,-1.1207021,5.34];
ybs[1912]=['',1.4776653,-0.1034159,6.05];
ybs[1913]=['',1.4760652,-0.2052911,6.11];
ybs[1914]=['',1.4806763,0.1318432,5.88];
ybs[1915]=['',1.4856452,0.4647783,6.37];
ybs[1916]=['β Dor',1.4566532,-1.0903792,3.76];
ybs[1917]=['',1.479603,-0.0837903,6.19];
ybs[1918]=['',1.4873102,0.5101094,5.96];
ybs[1919]=['',1.4980545,0.9336048,6.23];
ybs[1920]=['ν1 Col',1.4758662,-0.4862168,6.16];
ybs[1921]=['',1.4691852,-0.8255382,6.11];
ybs[1922]=['125 Tau',1.4890118,0.4521892,5.18];
ybs[1923]=['',1.4875561,0.3800377,6.34];
ybs[1924]=['',1.4634491,-1.0272378,6.75];
ybs[1925]=['σ Ori',1.4834328,-0.0451645,3.81];
ybs[1926]=['',1.4836003,-0.045063,6.65];
ybs[1927]=['',1.4827582,-0.1145205,5.96];
ybs[1928]=['ω Ori',1.485624,0.072141,4.57];
ybs[1929]=['ν2 Col',1.4778658,-0.5004994,5.31];
ybs[1930]=['',1.4627186,-1.0674615,6.32];
ybs[1931]=['49 Ori',1.483844,-0.1256787,4.8];
ybs[1932]=['',1.4930659,0.5474944,6.04];
ybs[1933]=['',1.493549,0.5573156,6.11];
ybs[1934]=['',1.486767,-0.06201,6];
ybs[1935]=['24 Cam',1.505815,0.9877018,6.05];
ybs[1936]=['',1.48648,-0.1692068,6.5];
ybs[1937]=['23 Cam',1.5114596,1.0731231,6.15];
ybs[1938]=['',1.4850676,-0.3113222,6.38];
ybs[1939]=['',1.4962418,0.5148394,6.43];
ybs[1940]=['126 Tau',1.4953642,0.2887575,4.86];
ybs[1941]=['',1.4814137,-0.7102633,5.82];
ybs[1942]=['ζ Ori',1.4922375,-0.0337147,2.05];
ybs[1943]=['ζ Ori',1.4922448,-0.0337147,4.21];
ybs[1944]=['',1.4916046,-0.0491109,6.22];
ybs[1945]=['',1.4983667,0.4073023,6.59];
ybs[1946]=['',1.4926421,-0.0195105,4.95];
ybs[1947]=['γ Men',1.4438724,-1.3321074,5.19];
ybs[1948]=['',1.4990105,0.3956749,6.36];
ybs[1949]=['',1.4937936,0.0060849,5.93];
ybs[1950]=['α Col',1.4858706,-0.5944999,2.64];
ybs[1951]=['',1.4919185,-0.1814855,6.52];
ybs[1952]=['',1.486749,-0.5692819,5.45];
ybs[1953]=['',1.4961831,-0.0503629,6.42];
ybs[1954]=['',1.4700724,-1.1614566,6.31];
ybs[1955]=['',1.5045352,0.4051544,6.21];
ybs[1956]=['',1.4956773,-0.2917317,6.21];
ybs[1957]=['51 Ori',1.4998768,0.0259141,4.91];
ybs[1958]=['',1.4579715,-1.2867651,5.78];
ybs[1959]=['',1.4980258,-0.3057824,6.15];
ybs[1960]=['',1.4937712,-0.5827619,6.34];
ybs[1961]=['',1.5013653,-0.1184432,6.02];
ybs[1962]=['12 Lep',1.4977694,-0.3903142,5.87];
ybs[1963]=['26 Cam',1.5209396,0.9795306,5.94];
ybs[1964]=['',1.5027085,-0.0279847,6.31];
ybs[1965]=['ο Aur',1.517568,0.869771,5.47];
ybs[1966]=['',1.4971792,-0.5327661,6.19];
ybs[1967]=['',1.4971922,-0.6048873,5.29];
ybs[1968]=['',1.5164558,0.7071226,6.58];
ybs[1969]=['',1.5028726,-0.3237224,5.73];
ybs[1970]=['',1.5332652,1.0963087,6.13];
ybs[1971]=['',1.5145886,0.3613372,6.95];
ybs[1972]=['',1.5111304,0.0701025,6.09];
ybs[1973]=['',1.5228019,0.7423537,6.29];
ybs[1974]=['',1.5076572,-0.3511161,6.34];
ybs[1975]=['',1.5023443,-0.6876139,6.25];
ybs[1976]=['',1.5074094,-0.3911758,6.15];
ybs[1977]=['γ Lep',1.5075026,-0.3916414,3.6];
ybs[1978]=['',1.5026432,-0.799772,6.39];
ybs[1979]=['129 Tau',1.5191715,0.2762849,6];
ybs[1980]=['',1.5152152,-0.074358,6.34];
ybs[1981]=['',1.5193765,0.1663234,5.79];
ybs[1982]=['',1.5177828,0.0205191,5.95];
ybs[1983]=['131 Tau',1.5211266,0.2529945,5.72];
ybs[1984]=['130 Tau',1.5222177,0.3095553,5.49];
ybs[1985]=['ι Men',1.4578668,-1.3754225,6.05];
ybs[1986]=['29 Cam',1.538767,0.9935091,6.54];
ybs[1987]=['133 Tau',1.5232605,0.2427162,5.29];
ybs[1988]=['',1.5514381,1.195111,6.2];
ybs[1989]=['τ Aur',1.530981,0.683943,4.52];
ybs[1990]=['μ Col',1.5136791,-0.5637123,5.17];
ybs[1991]=['',1.5264512,0.3643536,6.07];
ybs[1992]=['ζ Lep',1.518727,-0.258562,3.55];
ybs[1993]=['52 Ori',1.5241991,0.1127642,5.27];
ybs[1994]=['',1.5194127,-0.2832747,6.17];
ybs[1995]=['',1.521052,-0.1837121,6.03];
ybs[1996]=['132 Tau',1.5294389,0.42889,4.86];
ybs[1997]=['',1.5397364,0.8991855,6.29];
ybs[1998]=['κ Ori',1.5224488,-0.1686474,2.06];
ybs[1999]=['',1.518596,-0.4997186,6.22];
ybs[2000]=['30 Cam',1.5465891,1.0291876,6.14];
ybs[2001]=['',1.5262833,-0.0713542,5.97];
ybs[2002]=['',1.5146688,-0.8131379,5.31];
ybs[2003]=['',1.5191985,-0.6225141,6.32];
ybs[2004]=['134 Tau',1.5312053,0.2209049,4.91];
ybs[2005]=['υ Aur',1.5389929,0.6511892,4.74];
ybs[2006]=['ν Aur',1.5410805,0.6833518,3.97];
ybs[2007]=['',1.5381306,0.4882156,5.56];
ybs[2008]=['',1.5332453,0.1723798,5.8];
ybs[2009]=['δ Dor',1.5045681,-1.1471441,4.35];
ybs[2010]=['135 Tau',1.5353505,0.2497706,5.52];
ybs[2011]=['',1.521736,-0.7093992,6.61];
ybs[2012]=['',1.5403685,0.5607677,6.25];
ybs[2013]=['',1.5337789,0.0772965,5.97];
ybs[2014]=['β Pic',1.5178537,-0.8911476,3.85];
ybs[2015]=['',1.5303189,-0.2526846,5.49];
ybs[2016]=['π Men',1.4625398,-1.4042049,5.65];
ybs[2017]=['',1.5171826,-0.9486454,6.18];
ybs[2018]=['',1.534908,0.0354252,5.98];
ybs[2019]=['',1.5462083,0.6907718,6.45];
ybs[2020]=['',1.5311613,-0.4008316,5.87];
ybs[2021]=['31 Cam',1.5584031,1.0452895,5.2];
ybs[2022]=['',1.5458851,0.5920397,5.98];
ybs[2023]=['ξ Aur',1.55728,0.9723125,4.99];
ybs[2024]=['',1.5439304,0.3468345,6.06];
ybs[2025]=['55 Ori',1.5382915,-0.1311313,5.35];
ybs[2026]=['',1.5283909,-0.7831167,6.38];
ybs[2027]=['137 Tau',1.5435914,0.2474142,5.59];
ybs[2028]=['136 Tau',1.5484023,0.4819855,4.58];
ybs[2029]=['δ Lep',1.5375017,-0.3643255,3.81];
ybs[2030]=['',1.5380792,-0.4000529,6.17];
ybs[2031]=['56 Ori',1.5433658,0.0324478,4.78];
ybs[2032]=['',1.5480015,0.3543491,6.71];
ybs[2033]=['',1.5415429,-0.1577263,5.97];
ybs[2034]=['β Col',1.5351174,-0.6241858,3.12];
ybs[2035]=['',1.5711018,1.153607,6.25];
ybs[2036]=['γ Pic',1.5283453,-0.98019,4.51];
ybs[2037]=['',1.5399892,-0.513897,6.45];
ybs[2038]=['',1.5315922,-0.9208754,6.35];
ybs[2039]=['',1.5628947,0.9041773,6.49];
ybs[2040]=['',1.5558562,0.5533423,5.9];
ybs[2041]=['χ1 Ori',1.5526401,0.3539357,4.41];
ybs[2042]=['',1.5515148,0.1848299,6.12];
ybs[2043]=['',1.5334545,-0.9093799,5.17];
ybs[2044]=['',1.5529327,0.205344,6.59];
ybs[2045]=['',1.5513713,0.0563446,6.31];
ybs[2046]=['57 Ori',1.5550733,0.3447427,5.92];
ybs[2047]=['',1.5419488,-0.6567138,5.63];
ybs[2048]=['',1.5662831,0.8557457,6.47];
ybs[2049]=['',1.5429432,-0.672332,6.7];
ybs[2050]=['λ Col',1.5446408,-0.5898786,4.87];
ybs[2051]=['',1.5533427,0.0169487,6];
ybs[2052]=['',1.5524498,-0.0708782,6.57];
ybs[2053]=['',1.4206059,-1.4794451,6.2];
ybs[2054]=['',1.5490572,-0.3426963,6.69];
ybs[2055]=['α Ori',1.5555172,0.1293189,0.5];
ybs[2056]=['λ Men',1.5154189,-1.268764,6.53];
ybs[2057]=['',1.5589186,0.352156,5.4];
ybs[2058]=['ε Dor',1.5266015,-1.1675396,5.11];
ybs[2059]=['',1.5527701,-0.2054443,5.66];
ybs[2060]=['',1.5625964,0.5051647,6.32];
ybs[2061]=['',1.5596669,0.2430756,6.6];
ybs[2062]=['',1.5497937,-0.5086698,6.36];
ybs[2063]=['',1.5451846,-0.7490547,6.55];
ybs[2064]=['',1.5564701,-0.0805308,5.87];
ybs[2065]=['',1.5568338,-0.0835327,6.28];
ybs[2066]=['',1.5391596,-0.9974849,5.94];
ybs[2067]=['',1.5337938,-1.1175085,6.36];
ybs[2068]=['',1.5639741,0.4232615,6.02];
ybs[2069]=['',1.5612571,0.1660057,5.99];
ybs[2070]=['',1.5629071,0.2011072,5.87];
ybs[2071]=['δ Aur',1.5775211,0.9474416,3.72];
ybs[2072]=['',1.6080251,1.3191506,6.4];
ybs[2073]=['',1.5786761,0.9655226,6.44];
ybs[2074]=['',1.5787532,0.9520202,6.14];
ybs[2075]=['',1.5763268,0.8713426,5.89];
ybs[2076]=['',1.5519193,-0.6973505,5.57];
ybs[2077]=['',1.5480008,-0.8789243,6.52];
ybs[2078]=['139 Tau',1.5686856,0.4529937,4.82];
ybs[2079]=['η Lep',1.5599884,-0.2472429,3.71];
ybs[2080]=['',1.5588595,-0.3986043,5.96];
ybs[2081]=['ξ Col',1.5548317,-0.6478384,4.97];
ybs[2082]=['β Aur',1.5765935,0.7844774,1.9];
ybs[2083]=['',1.5502853,-0.8661016,6.1];
ybs[2084]=['',1.56031,-0.4051576,6.36];
ybs[2085]=['π Aur',1.5784514,0.8017421,4.26];
ybs[2086]=['σ Col',1.55889,-0.5476949,5.5];
ybs[2087]=['',1.5671778,0.0213858,6.22];
ybs[2088]=['',1.5506293,-0.9186081,5.29];
ybs[2089]=['θ Aur',1.5768687,0.6494748,2.62];
ybs[2090]=['',1.5800039,0.7782637,6.22];
ybs[2091]=['',1.5683582,-0.0173392,6.22];
ybs[2092]=['',1.5608994,-0.5580601,6.44];
ybs[2093]=['',1.5719588,0.2235568,5.7];
ybs[2094]=['59 Ori',1.5693995,0.0320706,5.9];
ybs[2095]=['36 Aur',1.5831715,0.8360268,5.73];
ybs[2096]=['',1.5458505,-1.1010679,4.65];
ybs[2097]=['60 Ori',1.5711787,0.0096582,5.22];
ybs[2098]=['',1.5459758,-1.1253674,6.63];
ybs[2099]=['',1.5865074,0.8544759,5.96];
ybs[2100]=['γ Col',1.5638426,-0.6157892,4.36];
ybs[2101]=['1 Mon',1.5715914,-0.1637465,6.12];
ybs[2102]=['2 Mon',1.5718239,-0.1668208,5.03];
ybs[2103]=['',1.5745936,-0.0252129,6.63];
ybs[2104]=['',1.5828127,0.5416377,5.98];
ybs[2105]=['',1.5819107,0.4812076,6.05];
ybs[2106]=['',1.5913706,0.8709772,6.05];
ybs[2107]=['',1.5763948,-0.0536613,4.53];
ybs[2108]=['',1.5610244,-0.9324355,6.45];
ybs[2109]=['',1.5911491,0.7570603,6.42];
ybs[2110]=['',1.5846451,0.390943,6.37];
ybs[2111]=['',1.5679814,-0.7685355,5.81];
ybs[2112]=['',1.5770173,-0.2251515,6.22];
ybs[2113]=['38 Aur',1.5928869,0.7489065,6.1];
ybs[2114]=['η Col',1.5703464,-0.7472626,3.96];
ybs[2115]=['',1.602718,1.0365397,6.34];
ybs[2116]=['',1.5905333,0.5695641,6.24];
ybs[2117]=['',1.5988445,0.9000678,6.45];
ybs[2118]=['μ Ori',1.5870792,0.1683488,4.12];
ybs[2119]=['κ Men',1.5212005,-1.385008,5.47];
ybs[2120]=['',1.6100678,1.1073931,6.39];
ybs[2121]=['',1.5863159,0.0295433,6.59];
ybs[2122]=['3 Mon',1.5838584,-0.1849961,4.95];
ybs[2123]=['',1.58046,-0.4436417,6.05];
ybs[2124]=['64 Ori',1.5922024,0.3436216,5.14];
ybs[2125]=['',1.580231,-0.591888,5.55];
ybs[2126]=['39 Aur',1.6005569,0.7501099,5.87];
ybs[2127]=['',1.5916536,0.2038265,6.08];
ybs[2128]=['1 Gem',1.5952754,0.4059713,4.16];
ybs[2129]=['χ2 Ori',1.5942524,0.3514319,4.63];
ybs[2130]=['',1.5868445,-0.2530567,6.2];
ybs[2131]=['',1.6001378,0.6625386,6.34];
ybs[2132]=['',1.5768907,-0.8939061,5.67];
ybs[2133]=['',1.6021454,0.5863499,6.23];
ybs[2134]=['',1.5893169,-0.4587891,5.04];
ybs[2135]=['',1.6047635,0.6175561,6.12];
ybs[2136]=['',1.5944304,-0.1171471,5.21];
ybs[2137]=['40 Aur',1.6068991,0.671574,5.36];
ybs[2138]=['63 Ori',1.5981868,0.0945384,5.67];
ybs[2139]=['66 Ori',1.5981486,0.0725231,5.63];
ybs[2140]=['',1.6054383,0.515016,6.08];
ybs[2141]=['',1.6109058,0.7304066,6.12];
ybs[2142]=['17 Lep',1.5973222,-0.2877653,4.93];
ybs[2143]=['',1.5937076,-0.5615659,5.65];
ybs[2144]=['',1.5996253,-0.1788328,5.87];
ybs[2145]=['',1.5815317,-1.0489136,6.45];
ybs[2146]=['37 Cam',1.6237899,1.0285094,5.36];
ybs[2147]=['',1.6149261,0.7164586,6.36];
ybs[2148]=['',1.6050951,-0.0732725,5.38];
ybs[2149]=['θ Lep',1.6024968,-0.2607395,4.67];
ybs[2150]=['',1.600346,-0.4223572,6.95];
ybs[2151]=['',1.5934107,-0.7860882,6.35];
ybs[2152]=['',1.5942653,-0.7868272,5.93];
ybs[2153]=['ν Ori',1.6099401,0.2576701,4.42];
ybs[2154]=['',1.5983643,-0.6198905,5.8];
ybs[2155]=['',1.6057595,-0.1950936,6.66];
ybs[2156]=['',1.594453,-0.8458153,6.58];
ybs[2157]=['',1.6037669,-0.4034285,5.47];
ybs[2158]=['',1.6014946,-0.5194539,5.81];
ybs[2159]=['36 Cam',1.6376173,1.1468551,5.32];
ybs[2160]=['',1.6056918,-0.3807825,5.78];
ybs[2161]=['',1.614977,0.1512218,6.55];
ybs[2162]=['19 Lep',1.6090168,-0.3345923,5.31];
ybs[2163]=['',1.6188713,0.3871819,5.93];
ybs[2164]=['',1.6054638,-0.5989343,5.83];
ybs[2165]=['π1 Col',1.6032717,-0.7383234,6.12];
ybs[2166]=['',1.6307309,0.918735,6.3];
ybs[2167]=['3 Gem',1.6197599,0.403295,5.75];
ybs[2168]=['',1.6154962,0.0435238,5.73];
ybs[2169]=['41 Aur',1.6296328,0.8500729,6.82];
ybs[2170]=['41 Aur',1.6296399,0.850039,6.09];
ybs[2171]=['θ Col',1.6073071,-0.650271,5.02];
ybs[2172]=['',1.6045697,-0.7870699,6.51];
ybs[2173]=['',1.6179399,-0.0997837,6.17];
ybs[2174]=['',1.6144107,-0.3915321,5.5];
ybs[2175]=['π2 Col',1.6085118,-0.73581,5.5];
ybs[2176]=['',1.6162428,-0.3164629,6.35];
ybs[2177]=['',1.6174207,-0.2546469,5.56];
ybs[2178]=['',1.6251811,0.3162965,6.33];
ybs[2179]=['5 Gem',1.6277,0.4260867,5.8];
ybs[2180]=['',1.6180301,-0.3975913,5.71];
ybs[2181]=['',1.611372,-0.774253,6.27];
ybs[2182]=['',1.6392442,0.8929758,6.04];
ybs[2183]=['',1.6316221,0.5704699,5.78];
ybs[2184]=['',1.6289637,0.3815486,6.56];
ybs[2185]=['',1.6268775,0.2379125,6.04];
ybs[2186]=['',1.6177814,-0.4661242,6.27];
ybs[2187]=['68 Ori',1.6295848,0.3452781,5.75];
ybs[2188]=['η1 Dor',1.5977835,-1.1526747,5.71];
ybs[2189]=['',1.6240847,-0.118003,6.15];
ybs[2190]=['',1.6025611,-1.0848789,5.05];
ybs[2191]=['6 Gem',1.6310284,0.3996902,6.39];
ybs[2192]=['69 Ori',1.6295625,0.2813988,4.95];
ybs[2193]=['ξ Ori',1.6289751,0.2478607,4.48];
ybs[2194]=['',1.6212049,-0.4740445,5.72];
ybs[2195]=['40 Cam',1.648789,1.0470089,5.35];
ybs[2196]=['',1.6272557,-0.0815574,6.18];
ybs[2197]=['',1.6186368,-0.7044131,5.58];
ybs[2198]=['',1.6144484,-0.8651344,6.49];
ybs[2199]=['',1.6277658,-0.1144533,5.05];
ybs[2200]=['',1.624069,-0.4623237,6.09];
ybs[2201]=['',1.6362305,0.3258842,6.58];
ybs[2202]=['',1.6343605,0.1853407,6.45];
ybs[2203]=['',1.6647896,1.2096479,4.8];
ybs[2204]=['',1.63176,-0.0438495,6.62];
ybs[2205]=['',1.6204086,-0.7904334,6.31];
ybs[2206]=['δ Pic',1.6178177,-0.9594926,4.81];
ybs[2207]=['',1.6312407,-0.310162,6.52];
ybs[2208]=['',1.640209,0.3123677,5.88];
ybs[2209]=['1 Lyn',1.6588081,1.0734465,4.98];
ybs[2210]=['η Gem',1.6421709,0.3926533,3.28];
ybs[2211]=['',1.6462956,0.6307413,6.92];
ybs[2212]=['',1.6367918,-0.0654503,5.83];
ybs[2213]=['κ Aur',1.6447192,0.5146705,4.35];
ybs[2214]=['71 Ori',1.6418883,0.3341805,5.2];
ybs[2215]=['ν Dor',1.6082485,-1.2016327,5.06];
ybs[2216]=['',1.6429287,0.2415831,5.91];
ybs[2217]=['72 Ori',1.6442377,0.2815821,5.3];
ybs[2218]=['',1.6398407,-0.0798905,5.83];
ybs[2219]=['',1.6352271,-0.4166175,6.39];
ybs[2220]=['',1.6340813,-0.5132045,6.54];
ybs[2221]=['γ Mon',1.6408312,-0.109675,3.98];
ybs[2222]=['42 Aur',1.6554687,0.8100627,6.52];
ybs[2223]=['73 Ori',1.6455261,0.2188876,5.33];
ybs[2224]=['8 Gem',1.6485255,0.4181783,6.08];
ybs[2225]=['',1.6448935,0.1057042,6.07];
ybs[2226]=['',1.6453208,0.0745926,6.64];
ybs[2227]=['',1.6441981,-0.0091081,5.65];
ybs[2228]=['',1.6436806,-0.0859549,5.99];
ybs[2229]=['',1.6485531,0.2996943,6.39];
ybs[2230]=['',1.6457005,0.0202342,6.37];
ybs[2231]=['',1.6432381,-0.1578615,6.1];
ybs[2232]=['2 Lyn',1.6658449,1.0297185,4.48];
ybs[2233]=['43 Aur',1.6585313,0.8089452,6.38];
ybs[2234]=['9 Gem',1.651401,0.4141718,6.25];
ybs[2235]=['74 Ori',1.6485391,0.2140129,5.04];
ybs[2236]=['',1.6414586,-0.3539796,5.91];
ybs[2237]=['',1.6422041,-0.322643,5.99];
ybs[2238]=['',1.6444106,-0.2395994,5.01];
ybs[2239]=['η2 Dor',1.6201229,-1.1448689,5.01];
ybs[2240]=['',1.6476748,0.0186781,6.63];
ybs[2241]=['75 Ori',1.6513476,0.1733448,5.39];
ybs[2242]=['',1.650629,0.122916,6.57];
ybs[2243]=['',1.6459397,-0.290208,5.92];
ybs[2244]=['',1.6534675,0.2451748,6.59];
ybs[2245]=['',1.6518342,0.0888306,5.71];
ybs[2246]=['',1.6445254,-0.520075,6.67];
ybs[2247]=['',1.6558306,0.2508318,6.16];
ybs[2248]=['',1.6497125,-0.3966338,6.07];
ybs[2249]=['6 Mon',1.6525619,-0.1873799,6.75];
ybs[2250]=['κ Col',1.6468243,-0.6134947,4.37];
ybs[2251]=['4 Lyn',1.6765295,1.0360008,5.94];
ybs[2252]=['',1.660049,0.3021735,6.32];
ybs[2253]=['',1.6581379,0.1577031,6.24];
ybs[2254]=['',1.6527669,-0.293681,5.14];
ybs[2255]=['α Men',1.6122284,-1.3047894,5.09];
ybs[2256]=['',1.6467364,-0.6854702,6];
ybs[2257]=['',1.6487017,-0.6588246,5.53];
ybs[2258]=['45 Aur',1.6744526,0.9326813,5.36];
ybs[2259]=['',1.6493373,-0.6503709,5.87];
ybs[2260]=['',1.6549492,-0.3486835,5.52];
ybs[2261]=['',1.6581024,-0.1640928,5.36];
ybs[2262]=['',1.6577271,-0.2624318,6.06];
ybs[2263]=['',1.6644668,0.2554947,5.69];
ybs[2264]=['',1.6593947,-0.1500654,6.22];
ybs[2265]=['',1.6581985,-0.3654319,5.81];
ybs[2266]=['',1.6701276,0.5153617,6.43];
ybs[2267]=['7 Mon',1.661965,-0.1367487,5.27];
ybs[2268]=['',1.6434479,-1.0336435,6.43];
ybs[2269]=['',1.6633905,-0.0516041,4.9];
ybs[2270]=['',1.6678382,0.2049641,6.54];
ybs[2271]=['',1.6705474,0.3098038,6.35];
ybs[2272]=['',1.6510921,-0.9205535,6.41];
ybs[2273]=['',1.6605362,-0.6005439,5.78];
ybs[2274]=['',1.6698749,0.0393658,6.31];
ybs[2275]=['',1.6553544,-0.8791313,7.04];
ybs[2276]=['ζ CMa',1.6635364,-0.5249198,3.02];
ybs[2277]=['',1.6349682,-1.2516043,6.64];
ybs[2278]=['',1.6692072,-0.2057114,5.64];
ybs[2279]=['',1.7061865,1.2307707,5.97];
ybs[2280]=['μ Gem',1.6774384,0.3926909,2.88];
ybs[2281]=['',1.6754506,0.2191463,6];
ybs[2282]=['',1.6645736,-0.5961416,5.53];
ybs[2283]=['ψ1 Aur',1.687663,0.8599711,4.91];
ybs[2284]=['',1.6612819,-0.8509047,6.6];
ybs[2285]=['',1.6951389,0.982074,5.64];
ybs[2286]=['',1.6781326,0.0654536,6.4];
ybs[2287]=['5 Lyn',1.6971202,1.019284,5.21];
ybs[2288]=['β CMa',1.674555,-0.3136294,1.98];
ybs[2289]=['',1.6780881,-0.0820562,6.67];
ybs[2290]=['δ Col',1.6711983,-0.5838091,3.85];
ybs[2291]=['',1.6861697,0.5182228,6.71];
ybs[2292]=['ε Mon',1.6801741,0.0799059,4.44];
ybs[2293]=['',1.6802033,0.0799544,6.72];
ybs[2294]=['',1.681532,0.1548115,6.26];
ybs[2295]=['',1.6788374,-0.1725923,6.19];
ybs[2296]=['',1.6855227,0.2799813,6.33];
ybs[2297]=['',1.6793392,-0.2633024,6.24];
ybs[2298]=['',1.6887705,0.4068595,6.06];
ybs[2299]=['',1.681261,-0.2014978,5.22];
ybs[2300]=['',1.6792414,-0.3455704,6.6];
ybs[2301]=['',1.6762111,-0.5550859,6.34];
ybs[2302]=['',1.688023,0.2566752,6.24];
ybs[2303]=['',1.6819409,-0.2264965,6.12];
ybs[2304]=['',1.6866048,0.1234029,5.98];
ybs[2305]=['',1.6795636,-0.4466697,5.63];
ybs[2306]=['',1.6867582,0.0259304,6.66];
ybs[2307]=['',1.6865117,-0.0167812,5.87];
ybs[2308]=['',1.7003693,0.8270801,6.56];
ybs[2309]=['',1.6888408,0.039384,6.51];
ybs[2310]=['',1.6793025,-0.640925,5.62];
ybs[2311]=['',1.6886225,-0.0681524,6.35];
ybs[2312]=['',1.6828875,-0.5025671,6.39];
ybs[2313]=['',1.6981749,0.5680383,6.43];
ybs[2314]=['ν Pic',1.6727646,-0.9840822,5.61];
ybs[2315]=['',1.6893063,-0.1380546,6.4];
ybs[2316]=['',1.6763224,-0.9109801,5.98];
ybs[2317]=['',1.6822251,-0.7033519,6.31];
ybs[2318]=['',1.692539,-0.0265887,5.87];
ybs[2319]=['',1.6920402,-0.0805184,6.15];
ybs[2320]=['α Car',1.6776764,-0.9199669,-0.72];
ybs[2321]=['',1.6940249,0.0143892,6.71];
ybs[2322]=['',1.6926749,-0.1313818,6.27];
ybs[2323]=['',1.6858722,-0.6122491,6.25];
ybs[2324]=['16 Gem',1.6990837,0.3574276,6.22];
ybs[2325]=['6 Lyn',1.7144031,1.0148022,5.88];
ybs[2326]=['48 Aur',1.7023181,0.5319001,5.55];
ybs[2327]=['',1.6956883,0.0504652,5.55];
ybs[2328]=['',1.6951006,0.0049327,5.2];
ybs[2329]=['',1.6952078,-0.0051079,5.55];
ybs[2330]=['',1.6761712,-1.0220321,6.48];
ybs[2331]=['',1.6719362,-1.1117185,6.27];
ybs[2332]=['47 Aur',1.7098944,0.814496,5.9];
ybs[2333]=['',1.7037651,0.4703688,6.47];
ybs[2334]=['',1.701172,0.2831102,6.23];
ybs[2335]=['',1.6813812,-0.9219053,6.51];
ybs[2336]=['',1.7002445,0.1795365,6.15];
ybs[2337]=['ν Gem',1.7035204,0.3524623,4.15];
ybs[2338]=['10 Mon',1.6980767,-0.0834124,5.06];
ybs[2339]=['',1.6778434,-1.0523567,5.8];
ybs[2340]=['',1.7647635,1.3888348,6.54];
ybs[2341]=['',1.6997486,0.033075,6.48];
ybs[2342]=['',1.6858852,-0.8411212,5.76];
ybs[2343]=['',1.6937538,-0.4515709,6.07];
ybs[2344]=['',1.7872566,1.4326966,6.65];
ybs[2345]=['',1.703271,0.1920182,6.59];
ybs[2346]=['π1 Dor',1.6685406,-1.2216876,5.56];
ybs[2347]=['',1.6927904,-0.6616876,6.48];
ybs[2348]=['',1.6781735,-1.1072972,6.46];
ybs[2349]=['',1.7040041,0.0458738,6.16];
ybs[2350]=['β Mon',1.7017202,-0.1230498,4.6];
ybs[2351]=['β Mon',1.7017564,-0.1230789,5.4];
ybs[2352]=['β Mon',1.7017564,-0.1230789,5.6];
ybs[2353]=['',1.7004175,-0.3051432,5.77];
ybs[2354]=['',1.6802217,-1.1142641,6.27];
ybs[2355]=['λ CMa',1.697673,-0.5689248,4.48];
ybs[2356]=['',1.7079495,0.1572701,6.57];
ybs[2357]=['',1.7638959,1.3608465,5.73];
ybs[2358]=['',1.6998017,-0.5652841,5.74];
ybs[2359]=['',1.749723,1.2858222,6.24];
ybs[2360]=['',1.7129904,0.2953044,6.2];
ybs[2361]=['',1.707577,-0.176272,5.93];
ybs[2362]=['',1.6994658,-0.7171857,6.32];
ybs[2363]=['',1.6905916,-1.0126123,5.82];
ybs[2364]=['',1.7126978,0.1960344,6.14];
ybs[2365]=['19 Gem',1.7149292,0.2772308,6.4];
ybs[2366]=['',1.7193772,0.5660976,5.87];
ybs[2367]=['',1.7091472,-0.2297994,6.16];
ybs[2368]=['',1.7148807,0.2054783,6.65];
ybs[2369]=['',1.7155318,0.2011523,5.23];
ybs[2370]=['7 Lyn',1.7304157,0.9657252,6.45];
ybs[2371]=['π2 Dor',1.6810429,-1.2165878,5.38];
ybs[2372]=['',1.7180824,0.2034006,6.03];
ybs[2373]=['',1.7126927,-0.2166062,5.15];
ybs[2374]=['',1.7092998,-0.4849919,5.93];
ybs[2375]=['',1.7148435,-0.1427206,5.43];
ybs[2376]=['12 Mon',1.7174977,0.084409,5.84];
ybs[2377]=['',1.7248995,0.5760232,6.42];
ybs[2378]=['',1.703545,-0.8771509,5.27];
ybs[2379]=['13 Mon',1.7201459,0.1276386,4.5];
ybs[2380]=['',1.717339,-0.1027728,5.6];
ybs[2381]=['ξ1 CMa',1.7142478,-0.4090621,4.33];
ybs[2382]=['',1.7108162,-0.6157162,5.84];
ybs[2383]=['',1.7012684,-0.9925757,5.22];
ybs[2384]=['',1.7094792,-0.7144507,6.2];
ybs[2385]=['',1.7234836,0.2467013,5.53];
ybs[2386]=['',1.718847,-0.1952356,6.24];
ybs[2387]=['',1.7123105,-0.645056,6.34];
ybs[2388]=['8 Lyn',1.7450453,1.0726466,5.94];
ybs[2389]=['',1.7229707,-0.0216524,5.1];
ybs[2390]=['',1.7602429,1.2518207,5.92];
ybs[2391]=['',1.7172635,-0.5593814,5.69];
ybs[2392]=['49 Aur',1.7311225,0.4887079,5.27];
ybs[2393]=['',1.7156291,-0.6582703,5.24];
ybs[2394]=['',1.7098689,-0.9048634,5.6];
ybs[2395]=['',1.7906359,1.3881683,5.45];
ybs[2396]=['11 Lyn',1.7441106,0.9919507,5.85];
ybs[2397]=['',1.7212987,-0.3655423,6.4];
ybs[2398]=['14 Mon',1.7283149,0.1317985,6.45];
ybs[2399]=['',1.7376116,0.67061,5.29];
ybs[2400]=['',1.7306852,0.1739572,5.88];
ybs[2401]=['',1.7191316,-0.6744859,6.44];
ybs[2402]=['',1.7021959,-1.1446949,6.29];
ybs[2403]=['',1.7301767,0.0151621,5.8];
ybs[2404]=['',1.7078943,-1.080329,6.15];
ybs[2405]=['',1.7221392,-0.6327261,5.42];
ybs[2406]=['μ Pic',1.7118987,-1.0257866,5.7];
ybs[2407]=['',1.7335394,0.0781171,6.55];
ybs[2408]=['ξ2 CMa',1.7282418,-0.401178,4.54];
ybs[2409]=['',1.7256834,-0.5713713,5.62];
ybs[2410]=['',1.7191119,-0.9136603,6.19];
ybs[2411]=['',1.7407701,0.428796,6.44];
ybs[2412]=['',1.7357073,-0.0913356,5.52];
ybs[2413]=['51 Aur',1.7468792,0.6870914,5.69];
ybs[2414]=['ψ3 Aur',1.7476196,0.69602,5.2];
ybs[2415]=['γ Gem',1.7415138,0.2858223,1.93];
ybs[2416]=['',1.7397228,0.1066873,6.06];
ybs[2417]=['ν1 CMa',1.734219,-0.3260602,5.7];
ybs[2418]=['',1.7289856,-0.6423027,5.59];
ybs[2419]=['53 Aur',1.7450531,0.5054643,5.79];
ybs[2420]=['',1.7408513,0.1890257,6.38];
ybs[2421]=['ψ2 Aur',1.7500383,0.7411555,4.79];
ybs[2422]=['',1.7361945,-0.2328786,5.97];
ybs[2423]=['ν2 CMa',1.7355155,-0.3364626,3.95];
ybs[2424]=['',1.7407594,0.0468005,6.17];
ybs[2425]=['',1.7312045,-0.6302457,6.35];
ybs[2426]=['',1.7417538,0.0861215,6.15];
ybs[2427]=['',1.7353502,-0.3950909,6.35];
ybs[2428]=['',1.7529254,0.7677652,6.41];
ybs[2429]=['',1.7257696,-0.9249626,4.39];
ybs[2430]=['',1.7477787,0.384099,6.04];
ybs[2431]=['',1.7401508,-0.2270264,6.12];
ybs[2432]=['54 Aur',1.7501115,0.4928667,6.03];
ybs[2433]=['',1.749801,0.4289348,6.38];
ybs[2434]=['',1.7434677,-0.0447971,6.14];
ybs[2435]=['',1.7458737,0.0816322,6.57];
ybs[2436]=['',1.7449179,0.0277569,6.21];
ybs[2437]=['ν3 CMa',1.7408277,-0.3187017,4.43];
ybs[2438]=['',1.736011,-0.6661722,6.04];
ybs[2439]=['',1.7350022,-0.7256907,6.34];
ybs[2440]=['',1.736952,-0.6459964,5.71];
ybs[2441]=['',1.739693,-0.5648301,5.27];
ybs[2442]=['',1.7439439,-0.2949046,6.03];
ybs[2443]=['',1.7504463,0.2261788,5.97];
ybs[2444]=['',1.747068,-0.2473031,4.82];
ybs[2445]=['ν Pup',1.7388365,-0.7543087,3.17];
ybs[2446]=['',1.7596082,0.6266923,6.46];
ybs[2447]=['25 Gem',1.7579528,0.4916847,6.42];
ybs[2448]=['',1.7533819,0.1107811,6.51];
ybs[2449]=['',1.7480553,-0.4139801,6.05];
ybs[2450]=['15 Mon',1.755485,0.1722797,4.66];
ybs[2451]=['',1.7574451,0.2857559,6.28];
ybs[2452]=['',1.756877,0.1916109,6.11];
ybs[2453]=['ψ4 Aur',1.7665592,0.7766442,5.02];
ybs[2454]=['',1.7481704,-0.5322219,5.71];
ybs[2455]=['',1.7555841,0.0082132,5.79];
ybs[2456]=['',1.7421889,-0.8420053,4.93];
ybs[2457]=['',1.7722351,0.9297312,6.27];
ybs[2458]=['',1.7666713,0.6478816,6.19];
ybs[2459]=['',1.7487404,-0.6664154,6.58];
ybs[2460]=['26 Gem',1.7620453,0.307523,5.21];
ybs[2461]=['',1.7597437,0.1103007,6.37];
ybs[2462]=['',1.7377933,-1.0743488,6.18];
ybs[2463]=['',1.758895,-0.1604374,5.19];
ybs[2464]=['12 Lyn',1.7819578,1.0369648,4.87];
ybs[2465]=['',1.7708802,0.6297687,6.31];
ybs[2466]=['ε Gem',1.7690583,0.4381594,2.98];
ybs[2467]=['',1.7644943,0.05249,6.19];
ybs[2468]=['',1.754207,-0.7046656,6.12];
ybs[2469]=['',1.7518205,-0.8325061,6.65];
ybs[2470]=['13 Lyn',1.7841954,0.9972966,5.35];
ybs[2471]=['30 Gem',1.768759,0.2304071,4.49];
ybs[2472]=['',1.7668727,0.068173,5.9];
ybs[2473]=['28 Gem',1.7728618,0.5051665,5.44];
ybs[2474]=['',1.7618989,-0.3922588,6.13];
ybs[2475]=['',1.7588726,-0.6706231,6.29];
ybs[2476]=['ψ5 Aur',1.7824252,0.7600798,5.25];
ybs[2477]=['ξ Gem',1.7744237,0.2245955,3.36];
ybs[2478]=['',1.7900613,0.9717185,6.33];
ybs[2479]=['',1.7900177,0.9717186,6.28];
ybs[2480]=['ψ6 Aur',1.7868874,0.8510354,5.22];
ybs[2481]=['',1.763686,-0.6845049,6.3];
ybs[2482]=['32 Gem',1.777091,0.2210647,6.46];
ybs[2483]=['42 Cam',1.8042083,1.1788151,5.14];
ybs[2484]=['α CMa',1.7725708,-0.2922227,-1.46];
ybs[2485]=['10 CMa',1.7689193,-0.5427473,5.2];
ybs[2486]=['',1.7708283,-0.4776653,6.45];
ybs[2487]=['16 Mon',1.7796957,0.1493884,5.93];
ybs[2488]=['',1.7733061,-0.4099618,6.05];
ybs[2489]=['',1.7714192,-0.5342979,6.54];
ybs[2490]=['',1.7763209,-0.2587211,5.32];
ybs[2491]=['',1.7838251,0.3170376,6.2];
ybs[2492]=['',1.7728359,-0.5553762,5.92];
ybs[2493]=['',1.7735015,-0.5406345,5.8];
ybs[2494]=['',1.7794085,-0.1768914,5.66];
ybs[2495]=['17 Mon',1.7831196,0.1397811,4.77];
ybs[2496]=['11 CMa',1.7801044,-0.2522673,5.29];
ybs[2497]=['',1.7478905,-1.2531411,6.51];
ybs[2498]=['18 Mon',1.7852015,0.0416012,4.47];
ybs[2499]=['',1.7753023,-0.6905825,6.62];
ybs[2500]=['',1.7836803,-0.1575475,5.07];
ybs[2501]=['12 CMa',1.780565,-0.367281,6.08];
ybs[2502]=['',1.7760603,-0.6597887,6.21];
ybs[2503]=['43 Cam',1.8165842,1.2017622,5.12];
ybs[2504]=['',1.7945668,0.5685735,5.71];
ybs[2505]=['',1.7715132,-0.911553,6.57];
ybs[2506]=['',1.7870425,-0.0235283,5.75];
ybs[2507]=['',1.7734987,-0.9152035,5.8];
ybs[2508]=['ψ7 Aur',1.7998517,0.728691,5.02];
ybs[2509]=['',1.7903873,0.0169751,6.15];
ybs[2510]=['',1.7810919,-0.6624912,5.26];
ybs[2511]=['33 Gem',1.7943741,0.2822711,5.85];
ybs[2512]=['14 Lyn',1.8118209,1.0370166,5.33];
ybs[2513]=['',1.791178,-0.040167,5.74];
ybs[2514]=['',1.7892883,-0.2648357,5.39];
ybs[2515]=['',1.7779213,-0.8952442,5.4];
ybs[2516]=['',1.7767132,-0.9550872,6.46];
ybs[2517]=['35 Gem',1.7968489,0.2335802,5.65];
ybs[2518]=['',1.7793151,-0.9698463,5.61];
ybs[2519]=['',1.848272,1.3428763,4.55];
ybs[2520]=['',1.792239,-0.4207203,6.33];
ybs[2521]=['36 Gem',1.8021296,0.3792646,5.27];
ybs[2522]=['',1.7980478,-0.0099693,5.77];
ybs[2523]=['',1.7588639,-1.2765985,6.37];
ybs[2524]=['',1.8103847,0.78204,6.26];
ybs[2525]=['',1.8041801,0.4113838,5.65];
ybs[2526]=['',1.7971891,-0.1408723,6.29];
ybs[2527]=['',1.7953198,-0.2986946,5.79];
ybs[2528]=['',1.7657661,-1.2297658,6.11];
ybs[2529]=['',1.7936764,-0.4775876,7.04];
ybs[2530]=['κ CMa',1.7922691,-0.5679006,3.96];
ybs[2531]=['59 Aur',1.8095055,0.6778405,6.12];
ybs[2532]=['θ Gem',1.8081667,0.5921816,3.6];
ybs[2533]=['60 Aur',1.810342,0.6703142,6.3];
ybs[2534]=['',1.8094468,0.6240746,6.01];
ybs[2535]=['',1.8017854,0.0525486,6.38];
ybs[2536]=['',1.796003,-0.4504384,6.33];
ybs[2537]=['',1.7947005,-0.5538999,5.7];
ybs[2538]=['',1.7919201,-0.7937706,6.55];
ybs[2539]=['ψ8 Aur',1.8135156,0.6714753,6.48];
ybs[2540]=['',1.7915892,-0.8140982,5.14];
ybs[2541]=['',1.7966639,-0.6003498,4.99];
ybs[2542]=['α Pic',1.782174,-1.0815794,3.27];
ybs[2543]=['',1.8071012,0.1457128,5.77];
ybs[2544]=['',1.8046173,-0.0933291,6.3];
ybs[2545]=['τ Pup',1.7913333,-0.8839113,2.93];
ybs[2546]=['',1.7906525,-0.9364008,4.4];
ybs[2547]=['',1.8096168,0.1913669,6.24];
ybs[2548]=['',1.8197392,0.7992439,6.34];
ybs[2549]=['',1.8195492,0.7657969,6.13];
ybs[2550]=['',1.8001807,-0.6328745,5.96];
ybs[2551]=['ζ Men',1.7364783,-1.4108668,5.64];
ybs[2552]=['15 Lyn',1.8299304,1.0190656,4.35];
ybs[2553]=['',1.8295678,1.004071,6.05];
ybs[2554]=['',1.7904535,-1.0520634,6.11];
ybs[2555]=['',1.7985857,-0.8433974,6.42];
ybs[2556]=['38 Gem',1.8152497,0.2294264,4.65];
ybs[2557]=['',1.8080993,-0.3327388,5.64];
ybs[2558]=['',1.8083147,-0.3310037,6.14];
ybs[2559]=['',1.80636,-0.4710477,6.4];
ybs[2560]=['ψ9 Aur',1.82538,0.8070462,5.87];
ybs[2561]=['37 Gem',1.8187054,0.4423104,5.73];
ybs[2562]=['',1.8122752,-0.1027086,6.41];
ybs[2563]=['15 CMa',1.8090706,-0.3535348,4.83];
ybs[2564]=['',1.8136376,-0.0202351,5.45];
ybs[2565]=['',1.8271391,0.8145615,5.86];
ybs[2566]=['θ CMa',1.8122299,-0.2106768,4.07];
ybs[2567]=['',1.8039305,-0.7423886,6.52];
ybs[2568]=['',1.808737,-0.4986683,6.04];
ybs[2569]=['',1.8148849,-0.0312239,6.21];
ybs[2570]=['',1.8104959,-0.4288494,6.21];
ybs[2571]=['',1.8043567,-0.7680703,6.46];
ybs[2572]=['ο1 CMa',1.8114297,-0.4226507,3.87];
ybs[2573]=['',1.8505793,1.2351901,5.68];
ybs[2574]=['',1.816057,-0.049504,6.04];
ybs[2575]=['',1.8118131,-0.4181913,6.91];
ybs[2576]=['',1.8191116,0.1447154,6.29];
ybs[2577]=['16 Lyn',1.8299923,0.7864405,4.9];
ybs[2578]=['',1.8265564,0.5872521,5.89];
ybs[2579]=['',1.8033772,-0.9445947,6.57];
ybs[2580]=['17 CMa',1.815594,-0.3567012,5.74];
ybs[2581]=['',1.822902,0.1731848,5.92];
ybs[2582]=['π CMa',1.8181304,-0.3520238,4.68];
ybs[2583]=['',1.811726,-0.7399825,6.32];
ybs[2584]=['',1.802558,-1.0362428,6.41];
ybs[2585]=['μ CMa',1.8205234,-0.24569,5];
ybs[2586]=['',1.8092422,-0.8838994,6.26];
ybs[2587]=['',1.8186834,-0.4009817,5.3];
ybs[2588]=['ι CMa',1.8205089,-0.2982342,4.37];
ybs[2589]=['',1.827339,0.2072279,6.27];
ybs[2590]=['',1.8188114,-0.5554199,6.36];
ybs[2591]=['',1.8246599,-0.1433404,6.34];
ybs[2592]=['62 Aur',1.8357141,0.6635217,6];
ybs[2593]=['39 Gem',1.8339204,0.4545897,6.1];
ybs[2594]=['ι Vol',1.7940353,-1.2390723,5.4];
ybs[2595]=['',1.8251193,-0.388115,6.61];
ybs[2596]=['',1.822308,-0.6174164,6.29];
ybs[2597]=['40 Gem',1.8368564,0.4516691,6.4];
ybs[2598]=['',1.8324915,0.1324188,6.27];
ybs[2599]=['',1.8263926,-0.4304811,5.46];
ybs[2600]=['',1.819154,-0.8509257,4.95];
ybs[2601]=['',2.0553591,1.5156728,5.07];
ybs[2602]=['',1.8336374,0.0622585,5.97];
ybs[2603]=['',1.8268659,-0.481218,6.23];
ybs[2604]=['',1.8246188,-0.6203162,6.23];
ybs[2605]=['',1.8354677,0.1270887,6.35];
ybs[2606]=['',1.8287181,-0.4747161,6.37];
ybs[2607]=['41 Gem',1.839885,0.2800038,5.68];
ybs[2608]=['',1.8308643,-0.4441631,5.59];
ybs[2609]=['',1.8701702,1.2338169,6.5];
ybs[2610]=['ε CMa',1.8307997,-0.5062679,1.5];
ybs[2611]=['',1.8296128,-0.5959657,5.06];
ybs[2612]=['',1.8451444,0.5651018,6.59];
ybs[2613]=['',1.8311496,-0.5416214,6.42];
ybs[2614]=['',1.8391672,-0.0942961,6.3];
ybs[2615]=['',1.8356582,-0.3776672,6.26];
ybs[2616]=['',1.839458,-0.1473549,5.96];
ybs[2617]=['',1.8378195,-0.3524618,6.31];
ybs[2618]=['',1.8300219,-0.7994099,6.22];
ybs[2619]=['',1.84056,-0.1612522,6.49];
ybs[2620]=['',1.8385438,-0.3866817,6.53];
ybs[2621]=['',1.8456413,0.0834513,6.63];
ybs[2622]=['ω Gem',1.8496343,0.4219885,5.18];
ybs[2623]=['',1.8493886,0.3092453,5.94];
ybs[2624]=['',1.8486943,0.2670195,5.74];
ybs[2625]=['',1.8466601,0.0963548,6.59];
ybs[2626]=['',1.8288468,-0.9732674,6.27];
ybs[2627]=['',1.8499228,0.2903703,5.82];
ybs[2628]=['',1.8462304,-0.0241256,6.17];
ybs[2629]=['',1.8399256,-0.4978627,6.27];
ybs[2630]=['',1.8285134,-0.984878,6.45];
ybs[2631]=['',1.8463092,-0.1005132,5.2];
ybs[2632]=['',1.8417855,-0.4407218,5.63];
ybs[2633]=['',1.8401776,-0.5847083,6.4];
ybs[2634]=['',1.8682626,1.0430537,6.44];
ybs[2635]=['',1.8546601,0.5113721,5.93];
ybs[2636]=['',1.8657837,0.9201194,6.12];
ybs[2637]=['',1.863059,0.8331543,6.38];
ybs[2638]=['σ CMa',1.8443554,-0.4881912,3.47];
ybs[2639]=['',1.8528338,0.1588381,5.97];
ybs[2640]=['19 Mon',1.8506161,-0.0746388,4.99];
ybs[2641]=['',1.85437,0.1904833,5.13];
ybs[2642]=['ζ Gem',1.8568577,0.3583546,3.79];
ybs[2643]=['',1.8554273,0.2191529,5.98];
ybs[2644]=['',1.8389413,-0.8977706,5.14];
ybs[2645]=['ο2 CMa',1.8502615,-0.4166217,3.02];
ybs[2646]=['',1.8570557,0.0253104,6.57];
ybs[2647]=['',1.8556924,-0.0935776,5.62];
ybs[2648]=['',1.8549226,-0.1773615,6.45];
ybs[2649]=['γ CMa',1.8538326,-0.2735124,4.12];
ybs[2650]=['',1.8457665,-0.7581882,6.43];
ybs[2651]=['44 Gem',1.8621803,0.3944175,6.02];
ybs[2652]=['',1.8666623,0.6009968,5.55];
ybs[2653]=['',1.8389941,-1.029326,6.02];
ybs[2654]=['',1.8317103,-1.1859741,5.17];
ybs[2655]=['',1.863103,0.1596437,5.78];
ybs[2656]=['',1.8580278,-0.385204,6.09];
ybs[2657]=['',1.8717743,0.592879,5.91];
ybs[2658]=['',1.8536464,-0.7395844,5.2];
ybs[2659]=['',1.8531514,-0.7617636,5.54];
ybs[2660]=['',1.8532602,-0.761827,6.79];
ybs[2661]=['',1.8716542,0.4910875,6.48];
ybs[2662]=['',1.8630754,-0.1867517,6.49];
ybs[2663]=['',1.8711257,0.3955559,7.68];
ybs[2664]=['',1.8523816,-0.8660604,4.93];
ybs[2665]=['',1.8754711,0.5897775,6.28];
ybs[2666]=['',1.8484582,-1.0335026,5.5];
ybs[2667]=['',1.8773599,0.6528283,6.16];
ybs[2668]=['',1.8692783,0.0850069,6.11];
ybs[2669]=['',1.8606321,-0.6076626,6.14];
ybs[2670]=['',1.8667585,-0.197809,5.39];
ybs[2671]=['',1.8663638,-0.217002,6.48];
ybs[2672]=['',1.8629452,-0.5357213,6.34];
ybs[2673]=['',1.9058734,1.2526683,6.35];
ybs[2674]=['',1.8725153,0.1296947,5.75];
ybs[2675]=['',1.8533684,-0.991131,5.17];
ybs[2676]=['45 Gem',1.8752292,0.2773389,5.44];
ybs[2677]=['',1.8626137,-0.6705865,6.11];
ybs[2678]=['',1.8669965,-0.4363337,6.08];
ybs[2679]=['',1.8583326,-0.8796247,6.46];
ybs[2680]=['',1.8674806,-0.465957,6.62];
ybs[2681]=['θ Men',1.8107336,-1.386719,5.45];
ybs[2682]=['',1.8692544,-0.4167862,5.71];
ybs[2683]=['',1.8671601,-0.7144145,5.79];
ybs[2684]=['',1.883076,0.3701049,6.43];
ybs[2685]=['δ CMa',1.8735609,-0.4613555,1.84];
ybs[2686]=['',1.8784055,-0.1813081,6.21];
ybs[2687]=['',1.8755517,-0.4203591,6.65];
ybs[2688]=['63 Aur',1.890781,0.6855327,4.9];
ybs[2689]=['τ Gem',1.8880013,0.5271449,4.41];
ybs[2690]=['',1.8666708,-0.9077,5.96];
ybs[2691]=['',1.8790971,-0.2840663,6.03];
ybs[2692]=['47 Gem',1.8888963,0.4680002,5.78];
ybs[2693]=['20 Mon',1.8825354,-0.0746776,4.92];
ybs[2694]=['',1.8748122,-0.6928335,4.83];
ybs[2695]=['',1.8993269,0.8968447,5.47];
ybs[2696]=['',1.8793939,-0.441084,5.69];
ybs[2697]=['',1.8816163,-0.3268423,6.23];
ybs[2698]=['48 Gem',1.8933659,0.4203716,5.85];
ybs[2699]=['21 Mon',1.8877748,-0.0060057,5.45];
ybs[2700]=['',1.8819377,-0.4805392,5.46];
ybs[2701]=['',1.9630706,1.4173226,6.31];
ybs[2702]=['',1.8900251,0.0979528,6.09];
ybs[2703]=['',1.8951591,0.4744098,6.43];
ybs[2704]=['',1.8593794,-1.2021142,6.47];
ybs[2705]=['',1.8911886,0.0948086,6.16];
ybs[2706]=['δ Mon',1.8898252,-0.0093411,4.15];
ybs[2707]=['18 Lyn',1.9113502,1.040086,5.2];
ybs[2708]=['',1.8882321,-0.3652159,5.84];
ybs[2709]=['51 Gem',1.8970728,0.2812697,5];
ybs[2710]=['26 CMa',1.8902204,-0.4535246,5.92];
ybs[2711]=['',1.8825588,-0.8547554,5.14];
ybs[2712]=['',1.8893819,-0.5386806,6.1];
ybs[2713]=['',1.9095974,0.8237112,5.58];
ybs[2714]=['',1.9021316,0.4305183,6.89];
ybs[2715]=['',1.8948713,-0.1971262,5.78];
ybs[2716]=['',1.8910104,-0.4802592,6.59];
ybs[2717]=['52 Gem',1.9032523,0.4335556,5.82];
ybs[2718]=['',1.8906263,-0.6385651,5.96];
ybs[2719]=['',1.8896504,-0.7075813,5.31];
ybs[2720]=['',1.9020149,0.2106938,5.62];
ybs[2721]=['',1.9007374,0.053539,5.35];
ybs[2722]=['',1.8955918,-0.3964789,6.01];
ybs[2723]=['',1.8998066,-0.0688553,5.75];
ybs[2724]=['',1.8998921,-0.1743802,5.9];
ybs[2725]=['',1.8973559,-0.4005505,6.36];
ybs[2726]=['',1.8962688,-0.4782154,6.12];
ybs[2727]=['γ1 Vol',1.8696203,-1.2311105,5.69];
ybs[2728]=['γ2 Vol',1.8698164,-1.2311401,3.78];
ybs[2729]=['',1.917583,0.909055,5.92];
ybs[2730]=['53 Gem',1.9088663,0.4861211,5.71];
ybs[2731]=['',1.9008008,-0.1808255,6.03];
ybs[2732]=['',1.8904719,-0.8168509,4.49];
ybs[2733]=['',1.8968615,-0.5432741,6.6];
ybs[2734]=['',1.9899988,1.4374049,4.96];
ybs[2735]=['',1.897636,-0.5302925,6.33];
ybs[2736]=['24 Mon',1.9049209,-0.0035913,6.41];
ybs[2737]=['27 CMa',1.8991485,-0.4607007,4.66];
ybs[2738]=['',1.8934973,-0.7893444,4.89];
ybs[2739]=['',1.9067042,0.1384603,5.82];
ybs[2740]=['',1.8949228,-0.7798646,5.1];
ybs[2741]=['ω CMa',1.9015649,-0.4680414,3.85];
ybs[2742]=['',1.9017264,-0.4726717,5.58];
ybs[2743]=['',1.9215505,0.8625182,5.05];
ybs[2744]=['',1.9062376,-0.1855017,5.95];
ybs[2745]=['64 Aur',1.918713,0.7127455,5.78];
ybs[2746]=['',1.8860967,-1.1036107,6.02];
ybs[2747]=['',1.9059884,-0.4151289,6.32];
ybs[2748]=['',1.9037144,-0.5363518,5.36];
ybs[2749]=['',1.9182557,0.5394782,6.24];
ybs[2750]=['',1.9083144,-0.272807,5.46];
ybs[2751]=['',1.901324,-0.7237859,5.94];
ybs[2752]=['',1.9138077,0.1158035,6.65];
ybs[2753]=['',1.9000931,-0.8184483,5.72];
ybs[2754]=['',1.8994096,-0.8432645,4.76];
ybs[2755]=['λ Gem',1.9176865,0.2878798,3.58];
ybs[2756]=['',1.9095946,-0.4077193,4.79];
ybs[2757]=['',1.9142929,-0.1173839,6.29];
ybs[2758]=['',1.9092432,-0.4874028,4.64];
ybs[2759]=['',1.9021144,-0.9170645,5.97];
ybs[2760]=['',1.9107024,-0.5400376,6.32];
ybs[2761]=['',1.9084227,-0.6695752,5.8];
ybs[2762]=['',1.9098182,-0.6394519,5.03];
ybs[2763]=['',1.9066207,-0.8171491,5.66];
ybs[2764]=['47 Cam',1.9391533,1.0446393,6.35];
ybs[2765]=['π Pup',1.9111771,-0.648264,2.7];
ybs[2766]=['',1.9146006,-0.4685022,6.46];
ybs[2767]=['',1.9319829,0.743647,6.35];
ybs[2768]=['',1.9332223,0.7885431,5.77];
ybs[2769]=['δ Gem',1.92678,0.3828397,3.53];
ybs[2770]=['',1.922713,0.0470174,5.89];
ybs[2771]=['',1.9247186,0.1238464,5.91];
ybs[2772]=['',1.9264466,0.2634692,6.45];
ybs[2773]=['29 CMa',1.9185206,-0.4294397,4.98];
ybs[2774]=['τ CMa',1.9186546,-0.4363389,4.4];
ybs[2775]=['19 Lyn',1.9409751,0.964044,6.53];
ybs[2776]=['19 Lyn',1.9410619,0.9639905,5.45];
ybs[2777]=['',1.9203393,-0.3373143,6.09];
ybs[2778]=['',1.9192161,-0.4648182,5.28];
ybs[2779]=['',1.9162784,-0.6419341,4.66];
ybs[2780]=['',1.9223622,-0.2869611,5.7];
ybs[2781]=['',1.9147538,-0.7685113,5.85];
ybs[2782]=['',1.9177187,-0.6420875,5.11];
ybs[2783]=['',1.9172285,-0.6851526,5.25];
ybs[2784]=['',1.9368431,0.6797665,6.4];
ybs[2785]=['65 Aur',1.9359188,0.6407506,5.13];
ybs[2786]=['',1.9204786,-0.5894622,6.3];
ybs[2787]=['56 Gem',1.9346644,0.3559681,5.1];
ybs[2788]=['',1.9290013,-0.2514584,5.45];
ybs[2789]=['',2.002823,1.410937,6.41];
ybs[2790]=['',1.9305842,-0.1557885,6.55];
ybs[2791]=['',1.9282683,-0.3996647,6.61];
ybs[2792]=['',1.9282043,-0.4714318,6.01];
ybs[2793]=['',1.9343284,0.0022527,5.99];
ybs[2794]=['',1.928933,-0.4527196,5.87];
ybs[2795]=['δ Vol',1.9059509,-1.1868604,3.98];
ybs[2796]=['',1.9497757,0.9047313,5.8];
ybs[2797]=['66 Aur',1.9453175,0.7090015,5.19];
ybs[2798]=['',1.9338604,-0.1575557,6.43];
ybs[2799]=['',1.9352964,-0.0528341,6.23];
ybs[2800]=['57 Gem',1.9415409,0.4363594,5.03];
ybs[2801]=['',1.9627847,1.1568081,6.47];
ybs[2802]=['58 Gem',1.9414214,0.3996155,6.02];
ybs[2803]=['',1.9356808,-0.1052629,5.82];
ybs[2804]=['',1.934281,-0.3327445,4.96];
ybs[2805]=['',1.9240024,-0.9138314,6.05];
ybs[2806]=['',1.9240244,-0.9137975,6.6];
ybs[2807]=['',1.9252853,-0.9098975,5.39];
ybs[2808]=['59 Gem',1.9463834,0.4815092,5.76];
ybs[2809]=['',1.9454238,0.2699622,6.41];
ybs[2810]=['21 Lyn',1.9571911,0.8580125,4.64];
ybs[2811]=['',1.9370877,-0.5580251,5.43];
ybs[2812]=['1 CMi',1.9474837,0.202806,5.3];
ybs[2813]=['ι Gem',1.951494,0.4842904,3.79];
ybs[2814]=['',1.9393727,-0.4866508,5.38];
ybs[2815]=['',1.93935,-0.5628879,5.39];
ybs[2816]=['',1.9410852,-0.5282421,6.6];
ybs[2817]=['',1.9450652,-0.2836277,5.33];
ybs[2818]=['',1.9431061,-0.4007643,6.19];
ybs[2819]=['η CMa',1.9419616,-0.5122935,2.45];
ybs[2820]=['ε CMi',1.9503473,0.1610231,4.99];
ybs[2821]=['',1.9410654,-0.6263444,6.31];
ybs[2822]=['',1.9783512,1.1940179,5.64];
ybs[2823]=['',1.9457165,-0.3326924,6.24];
ybs[2824]=['',1.9472181,-0.2408863,5.78];
ybs[2825]=['',1.9506418,-0.1016697,5.97];
ybs[2826]=['',1.9446027,-0.5560345,5.35];
ybs[2827]=['',1.9560424,0.3749836,6.54];
ybs[2828]=['',1.9539631,0.1842668,6.37];
ybs[2829]=['61 Gem',1.9564294,0.3526667,5.93];
ybs[2830]=['',1.9515991,-0.0800732,6.76];
ybs[2831]=['',1.9477096,-0.384543,6.05];
ybs[2832]=['',1.9549534,0.1912605,6.41];
ybs[2833]=['',1.9479533,-0.4410051,5.78];
ybs[2834]=['',1.9445365,-0.6516983,6.97];
ybs[2835]=['',1.9445437,-0.6517128,6.84];
ybs[2836]=['',1.9664575,0.8400587,5.72];
ybs[2837]=['β CMi',1.9568596,0.1437881,2.9];
ybs[2838]=['63 Gem',1.9599785,0.3733896,5.22];
ybs[2839]=['',1.9489139,-0.5548176,6.31];
ybs[2840]=['',1.736038,-1.5172635,6.47];
ybs[2841]=['22 Lyn',1.971268,0.8660295,5.36];
ybs[2842]=['',1.9535089,-0.4147404,6.56];
ybs[2843]=['η CMi',1.9606679,0.1202613,5.25];
ybs[2844]=['ρ Gem',1.96646,0.5538332,4.18];
ybs[2845]=['',1.9557533,-0.3126821,5.63];
ybs[2846]=['γ CMi',1.9613056,0.1548805,4.32];
ybs[2847]=['',1.954898,-0.4038158,5.61];
ybs[2848]=['',1.953097,-0.5967537,5.9];
ybs[2849]=['64 Gem',1.9672684,0.4898406,5.05];
ybs[2850]=['',1.9642808,0.2628035,6.22];
ybs[2851]=['',1.9591922,-0.2026027,5.79];
ybs[2852]=['',1.9580729,-0.3998715,5.95];
ybs[2853]=['65 Gem',1.9693154,0.4863116,5.01];
ybs[2854]=['',1.9503775,-0.8913169,5.1];
ybs[2855]=['',1.958958,-0.5097616,5.54];
ybs[2856]=['6 CMi',1.9685535,0.2086403,4.54];
ybs[2857]=['',1.9658969,-0.0341637,5.59];
ybs[2858]=['',1.9661875,-0.132703,5.86];
ybs[2859]=['',1.9658164,-0.1811448,5.75];
ybs[2860]=['',1.9656133,-0.2626951,6.05];
ybs[2861]=['',1.960194,-0.6608134,6.58];
ybs[2862]=['',1.9626092,-0.5567628,6.38];
ybs[2863]=['',1.9626237,-0.5567386,7.13];
ybs[2864]=['',1.9791521,0.6779327,6.54];
ybs[2865]=['',1.9636195,-0.5499243,5.77];
ybs[2866]=['',1.9674138,-0.402767,4.85];
ybs[2867]=['',1.9632533,-0.6783075,5.43];
ybs[2868]=['',1.9724889,-0.0921427,6.24];
ybs[2869]=['',1.9775293,0.2972738,5.42];
ybs[2870]=['σ Pup',1.9635386,-0.7566592,3.25];
ybs[2871]=['',1.982302,0.3985219,6.54];
ybs[2872]=['δ1 CMi',1.9782023,0.0324762,5.25];
ybs[2873]=['',1.9707471,-0.5413151,4.65];
ybs[2874]=['',1.9703818,-0.6526231,6.65];
ybs[2875]=['',1.9777798,-0.1559365,5.9];
ybs[2876]=['',1.9660621,-0.919849,5.87];
ybs[2877]=['',1.9736111,-0.631914,6.68];
ybs[2878]=['68 Gem',1.9853361,0.2752752,5.25];
ybs[2879]=['δ2 CMi',1.9830325,0.0564785,5.59];
ybs[2880]=['',1.9593512,-1.1268129,6.39];
ybs[2881]=['',1.9748636,-0.6272916,6.61];
ybs[2882]=['α Gem',1.990399,0.5555989,2.88];
ybs[2883]=['α Gem',1.990399,0.555594,1.98];
ybs[2884]=['',1.968151,-0.950368,5.96];
ybs[2885]=['',1.9872043,0.1834957,6.28];
ybs[2886]=['',2.0017526,0.9721281,5.92];
ybs[2887]=['',1.9777338,-0.6285826,6.3];
ybs[2888]=['',1.9927272,0.5394058,5.33];
ybs[2889]=['',1.9831123,-0.2511997,6.21];
ybs[2890]=['',1.9968886,0.7500588,6.3];
ybs[2891]=['',1.9827118,-0.33976,5.66];
ybs[2892]=['',1.9817633,-0.4322315,5.85];
ybs[2893]=['δ3 CMi',1.9877041,0.0578839,5.81];
ybs[2894]=['',1.9849815,-0.2544424,4.97];
ybs[2895]=['',1.9997102,0.8050162,5.65];
ybs[2896]=['',1.9898755,0.0465976,6.55];
ybs[2897]=['υ Gem',1.995913,0.4684461,4.06];
ybs[2898]=['',1.985765,-0.3900953,4.45];
ybs[2899]=['',1.9812074,-0.7001053,6.26];
ybs[2900]=['',1.9809862,-0.7529449,6.52];
ybs[2901]=['',1.9868325,-0.4106489,5.83];
ybs[2902]=['',1.9868689,-0.4106684,5.87];
ybs[2903]=['',1.9841666,-0.6351757,5.54];
ybs[2904]=['',1.9874526,-0.4567804,6.65];
ybs[2905]=['',1.9859159,-0.585001,6.11];
ybs[2906]=['',2.0059031,0.8502653,5.92];
ybs[2907]=['',2.0026291,0.6975844,6.38];
ybs[2908]=['',1.9878538,-0.472407,5.77];
ybs[2909]=['',1.9844925,-0.6974412,6.76];
ybs[2910]=['',1.9978995,0.1013255,5.91];
ybs[2911]=['ε Men',1.9383778,-1.3813158,5.53];
ybs[2912]=['',1.9960433,-0.1460376,6.27];
ybs[2913]=['',1.9948754,-0.2539209,5.7];
ybs[2914]=['',1.9912731,-0.4961073,4.64];
ybs[2915]=['',1.9948358,-0.3877492,6.34];
ybs[2916]=['70 Gem',2.0077802,0.6107138,5.56];
ybs[2917]=['',1.9865605,-0.8993619,6.28];
ybs[2918]=['',2.0059131,0.4241706,6.27];
ybs[2919]=['25 Mon',2.0005753,-0.0727386,5.13];
ybs[2920]=['',1.9973579,-0.3448485,5.74];
ybs[2921]=['23 Lyn',2.019524,0.9952591,6.06];
ybs[2922]=['ο Gem',2.0104503,0.6026019,4.9];
ybs[2923]=['',2.0100792,0.4217569,6.17];
ybs[2924]=['',2.0017879,-0.2530341,6.53];
ybs[2925]=['',1.9997835,-0.4159373,6.37];
ybs[2926]=['',1.9908374,-0.9178572,4.94];
ybs[2927]=['',2.0153706,0.6682207,5.73];
ybs[2928]=['',2.0135239,0.5576625,6.17];
ybs[2929]=['',1.9996077,-0.6113029,4.53];
ybs[2930]=['74 Gem',2.0110103,0.3074743,5.05];
ybs[2931]=['',2.0202671,0.8390302,5.56];
ybs[2932]=['',1.9958627,-0.8532276,5.72];
ybs[2933]=['',1.9920659,-0.9763919,6.39];
ybs[2934]=['',2.0012421,-0.6166928,6.6];
ybs[2935]=['α CMi',2.0097569,0.0901879,0.38];
ybs[2936]=['',2.0041634,-0.4436928,4.7];
ybs[2937]=['',2.0011105,-0.6643983,6.38];
ybs[2938]=['24 Lyn',2.0292428,1.0236438,4.99];
ybs[2939]=['',2.0080295,-0.3270159,6.72];
ybs[2940]=['',2.0063735,-0.468777,4.5];
ybs[2941]=['',2.0064098,-0.4688111,4.62];
ybs[2942]=['',2.0133128,0.090282,6.02];
ybs[2943]=['',2.0177741,0.4007285,5.89];
ybs[2944]=['',2.0038249,-0.6989764,6.59];
ybs[2945]=['',2.0165734,0.2393264,6.24];
ybs[2946]=['',2.0054729,-0.6379901,5.8];
ybs[2947]=['',2.0045128,-0.6778545,6.19];
ybs[2948]=['',2.0091062,-0.4698544,6.5];
ybs[2949]=['',2.002748,-0.8492427,5.68];
ybs[2950]=['',2.0148786,-0.1438869,6.01];
ybs[2951]=['',2.013706,-0.2674153,4.94];
ybs[2952]=['',2.0128259,-0.3441596,5.93];
ybs[2953]=['',2.0085304,-0.6696115,4.84];
ybs[2954]=['',2.0259408,0.5923777,6.02];
ybs[2955]=['',2.0097342,-0.6666665,5.73];
ybs[2956]=['',2.0100179,-0.6687857,5.76];
ybs[2957]=['',2.0212514,0.2342501,5.77];
ybs[2958]=['',2.0196658,0.0622363,5.94];
ybs[2959]=['',2.022109,0.2469504,5.56];
ybs[2960]=['',2.0108016,-0.6568949,6];
ybs[2961]=['',2.0329477,0.8791846,5.27];
ybs[2962]=['α Mon',2.0176844,-0.1677216,3.93];
ybs[2963]=['',2.0053746,-0.9307946,6.06];
ybs[2964]=['',2.0146207,-0.4887643,6.76];
ybs[2965]=['σ Gem',2.0282362,0.5030698,4.28];
ybs[2966]=['',2.0151118,-0.5536046,6.56];
ybs[2967]=['51 Cam',2.0463438,1.1413406,5.92];
ybs[2968]=['',2.0178029,-0.3908822,6.18];
ybs[2969]=['49 Cam',2.0449115,1.0955234,6.49];
ybs[2970]=['',2.0281924,0.3898996,6.21];
ybs[2971]=['',1.9846171,-1.2973127,7.16];
ybs[2972]=['',1.9846244,-1.2973128,7.26];
ybs[2973]=['',2.016406,-0.6735604,5.42];
ybs[2974]=['',2.0260963,0.0022655,6.19];
ybs[2975]=['76 Gem',2.0315931,0.448967,5.31];
ybs[2976]=['',2.0164318,-0.7800017,6.41];
ybs[2977]=['κ Gem',2.0329766,0.4247718,3.57];
ybs[2978]=['',2.0194775,-0.6734846,6.54];
ybs[2979]=['',2.0315653,0.2233875,6.43];
ybs[2980]=['',2.0237694,-0.4609508,5.64];
ybs[2981]=['',2.0306877,0.0409245,6.47];
ybs[2982]=['β Gem',2.0369311,0.488085,1.14];
ybs[2983]=['79 Gem',2.0358883,0.3535269,6.33];
ybs[2984]=['',2.0275193,-0.4461719,6.55];
ybs[2985]=['1 Pup',2.026895,-0.4969113,4.59];
ybs[2986]=['',2.02501,-0.6302363,5.6];
ybs[2987]=['',2.024479,-0.679347,6.89];
ybs[2988]=['3 Pup',2.0280401,-0.5064016,3.96];
ybs[2989]=['',2.0957055,1.3997233,6.56];
ybs[2990]=['',2.0233361,-0.7894608,5.06];
ybs[2991]=['',2.0432743,0.6537281,5.18];
ybs[2992]=['',1.9855087,-1.3559359,6.18];
ybs[2993]=['',2.0271323,-0.6677947,6.4];
ybs[2994]=['',2.0268876,-0.7154757,5.17];
ybs[2995]=['81 Gem',2.0400364,0.3219905,4.88];
ybs[2996]=['',2.0315515,-0.4316943,5.62];
ybs[2997]=['',2.0236445,-0.873577,6.57];
ybs[2998]=['',2.0185506,-1.0243297,6.43];
ybs[2999]=['',2.0292013,-0.6304635,5.8];
ybs[3000]=['11 CMi',2.0403689,0.1868719,5.3];
ybs[3001]=['2 Pup',2.0359519,-0.2573836,6.89];
ybs[3002]=['2 Pup',2.0359808,-0.2574661,6.07];
ybs[3003]=['',2.0308816,-0.6632842,5.88];
ybs[3004]=['',2.021745,-1.0173405,6.21];
ybs[3005]=['π Gem',2.0467513,0.582128,5.14];
ybs[3006]=['',2.0386938,-0.1192705,5.49];
ybs[3007]=['4 Pup',2.038008,-0.2552548,5.04];
ybs[3008]=['',2.0330975,-0.6623242,6.54];
ybs[3009]=['',2.0338714,-0.6637366,3.61];
ybs[3010]=['',2.0355216,-0.5974949,5.37];
ybs[3011]=['',2.0415739,-0.2222998,6.39];
ybs[3012]=['',2.0336974,-0.7646797,6.03];
ybs[3013]=['82 Gem',2.0508603,0.4027959,6.18];
ybs[3014]=['',2.0378838,-0.6631392,5.88];
ybs[3015]=['',2.0431786,-0.3941216,5.9];
ybs[3016]=['ζ Vol',2.0136728,-1.2682374,3.95];
ybs[3017]=['',2.039427,-0.7002454,6.57];
ybs[3018]=['',2.0453419,-0.2801753,6.34];
ybs[3019]=['',2.0458283,-0.2805884,6.43];
ybs[3020]=['',2.0639218,0.9436141,6.02];
ybs[3021]=['5 Pup',2.0468152,-0.2138946,5.48];
ybs[3022]=['',2.052524,0.2322685,6.04];
ybs[3023]=['',2.0337623,-0.9910552,6.12];
ybs[3024]=['',2.0418309,-0.6875386,6.31];
ybs[3025]=['',2.0519566,0.0745254,6.53];
ybs[3026]=['ο Pup',2.0468607,-0.4537762,4.5];
ybs[3027]=['',2.0432787,-0.6732251,5.08];
ybs[3028]=['',2.0284329,-1.1542241,6.38];
ybs[3029]=['',2.043195,-0.8145535,5.23];
ybs[3030]=['',2.025177,-1.2196583,6.18];
ybs[3031]=['',2.0706397,0.9624545,6.38];
ybs[3032]=['',2.0621476,0.5789199,6.03];
ybs[3033]=['',2.0463113,-0.710601,6.14];
ybs[3034]=['',2.053455,-0.2341593,6.23];
ybs[3035]=['',2.0510181,-0.4358953,5.33];
ybs[3036]=['6 Pup',2.054218,-0.3017925,5.18];
ybs[3037]=['ξ Pup',2.0521845,-0.4349815,3.34];
ybs[3038]=['',2.0467108,-0.8227495,4.71];
ybs[3039]=['',2.0566698,-0.1613857,5.61];
ybs[3040]=['',2.0543866,-0.3537796,6.56];
ybs[3041]=['',2.0514519,-0.6162087,5.93];
ybs[3042]=['',2.0598073,0.0560859,6.18];
ybs[3043]=['',2.0559065,-0.3418564,6.12];
ybs[3044]=['',2.053064,-0.5821005,5.6];
ybs[3045]=['',2.0654622,0.3361661,5.99];
ybs[3046]=['',2.05984,-0.1953439,6.16];
ybs[3047]=['',2.0506978,-0.8104631,4.11];
ybs[3048]=['',2.0457098,-0.9866929,6.33];
ybs[3049]=['',2.051832,-0.7821667,6.32];
ybs[3050]=['',2.050558,-0.8189179,5.84];
ybs[3051]=['ζ CMi',2.0637221,0.0297183,5.14];
ybs[3052]=['',2.0596426,-0.4292134,6.45];
ybs[3053]=['',2.0656101,0.0560737,6.31];
ybs[3054]=['',2.0491771,-0.9856434,5.59];
ybs[3055]=['8 Pup',2.0630993,-0.2248614,6.36];
ybs[3056]=['9 Pup',2.0634503,-0.2436875,5.17];
ybs[3057]=['25 Lyn',2.0780889,0.8258952,6.25];
ybs[3058]=['26 Lyn',2.0790764,0.8290105,5.45];
ybs[3059]=['φ Gem',2.0725512,0.4660137,4.97];
ybs[3060]=['',2.062909,-0.370674,5.63];
ybs[3061]=['',2.0572464,-0.7791722,6.45];
ybs[3062]=['',2.0491066,-1.0532414,5.78];
ybs[3063]=['',2.0554179,-0.882667,5.91];
ybs[3064]=['',2.0682507,-0.0958679,5.76];
ybs[3065]=['10 Pup',2.0657847,-0.2602441,5.69];
ybs[3066]=['',2.060125,-0.7532746,6.32];
ybs[3067]=['',2.1078305,1.2889084,5.41];
ybs[3068]=['',2.0521941,-1.0491899,6.72];
ybs[3069]=['',2.0875011,0.9850215,6.72];
ybs[3070]=['',2.0615806,-0.7496609,6.04];
ybs[3071]=['',2.0646543,-0.6068456,5.01];
ybs[3072]=['',2.0641135,-0.7093054,3.73];
ybs[3073]=['',2.0500628,-1.1564322,5.79];
ybs[3074]=['',2.1316257,1.3859331,5.42];
ybs[3075]=['',2.0825014,0.6169111,6.23];
ybs[3076]=['',2.0660904,-0.6794158,4.49];
ybs[3077]=['',2.0680414,-0.635801,5.43];
ybs[3078]=['85 Gem',2.0817101,0.3458817,5.35];
ybs[3079]=['',2.0806643,0.1535289,5.86];
ybs[3080]=['',2.0641969,-0.9500118,5.7];
ybs[3081]=['',2.0671343,-0.8670417,4.63];
ybs[3082]=['',2.068318,-0.8406895,4.24];
ybs[3083]=['',2.0729814,-0.6273222,5.49];
ybs[3084]=['',2.0751414,-0.6093401,6.15];
ybs[3085]=['',2.0843206,0.077129,6.17];
ybs[3086]=['',2.0943459,0.7663697,6.34];
ybs[3087]=['1 Cnc',2.0873256,0.2744229,5.78];
ybs[3088]=['',2.0778044,-0.5407637,6.44];
ybs[3089]=['',2.0882541,0.1496491,6.05];
ybs[3090]=['',2.0880018,0.0184975,6.35];
ybs[3091]=['',2.0828639,-0.5297398,6.33];
ybs[3092]=['',2.0753645,-0.9188955,6.38];
ybs[3093]=['',2.0794378,-0.7663953,6.02];
ybs[3094]=['11 Pup',2.0852889,-0.400498,4.2];
ybs[3095]=['',2.0918298,0.1247222,6.41];
ybs[3096]=['',2.0940495,0.2871209,5.99];
ybs[3097]=['',2.0743086,-1.0012733,5.63];
ybs[3098]=['',2.1091551,1.0293618,5.77];
ybs[3099]=['',2.0823963,-0.7121457,6.78];
ybs[3100]=['',2.1922315,1.4657224,6.49];
ybs[3101]=['53 Cam',2.1109085,1.0516453,6.01];
ybs[3102]=['14 CMi',2.0927246,0.0376476,5.29];
ybs[3103]=['',2.0847339,-0.7412929,6.09];
ybs[3104]=['',2.1148592,1.0999106,6.4];
ybs[3105]=['',2.0884867,-0.5306145,4.79];
ybs[3106]=['',2.0846598,-0.7603896,5.35];
ybs[3107]=['',2.0985469,0.2299277,6.02];
ybs[3108]=['',2.0861099,-0.7710295,5.09];
ybs[3109]=['χ Car',2.0830596,-0.9258782,3.47];
ybs[3110]=['',2.0859419,-0.8370125,6.22];
ybs[3111]=['',2.1143787,0.9983907,6.49];
ybs[3112]=['',2.080089,-1.0575436,5.74];
ybs[3113]=['',2.0884358,-0.7966568,5.17];
ybs[3114]=['27 Mon',2.0985841,-0.0654167,4.93];
ybs[3115]=['12 Pup',2.0950347,-0.4080328,5.11];
ybs[3116]=['ω1 Cnc',2.1049123,0.4419823,5.83];
ybs[3117]=['',2.1041035,0.3446525,6.25];
ybs[3118]=['',2.0826297,-1.0331141,6.25];
ybs[3119]=['',2.1051829,0.4103961,6.34];
ybs[3120]=['3 Cnc',2.1039526,0.3008886,5.55];
ybs[3121]=['',2.0897912,-0.860665,4.41];
ybs[3122]=['',2.1096983,0.6168598,6.34];
ybs[3123]=['',2.0986059,-0.32232,4.61];
ybs[3124]=['ω2 Cnc',2.1083805,0.4366859,6.31];
ybs[3125]=['',2.090105,-0.8991262,6.44];
ybs[3126]=['5 Cnc',2.1070451,0.2859888,5.99];
ybs[3127]=['',2.1029765,-0.051497,6.51];
ybs[3128]=['',2.105415,0.0839601,5.65];
ybs[3129]=['',2.093559,-0.7903551,5.99];
ybs[3130]=['',2.0865568,-1.0536639,5.6];
ybs[3131]=['',2.0835617,-1.1059066,6.14];
ybs[3132]=['',2.0958864,-0.6870554,5.24];
ybs[3133]=['28 Mon',2.1051535,-0.0255105,4.68];
ybs[3134]=['',2.0939447,-0.8734436,6.32];
ybs[3135]=['',2.094025,-0.8733905,6.34];
ybs[3136]=['',2.1082445,0.1543639,6.22];
ybs[3137]=['',2.1098392,0.0395274,4.39];
ybs[3138]=['',2.099206,-0.79457,6.61];
ybs[3139]=['',2.0911134,-1.0627685,5.81];
ybs[3140]=['',2.0985884,-0.8560822,6.02];
ybs[3141]=['χ Gem',2.116299,0.4838712,4.94];
ybs[3142]=['',2.1102631,-0.1118229,6.33];
ybs[3143]=['',2.0996171,-0.8541644,6.12];
ybs[3144]=['',2.0948273,-1.0520084,6.33];
ybs[3145]=['',2.0945816,-1.0586305,5.17];
ybs[3146]=['',2.1053939,-0.6519307,5.95];
ybs[3147]=['',2.1075034,-0.6478674,6.34];
ybs[3148]=['',2.1006541,-0.9463205,5.87];
ybs[3149]=['',2.1030195,-0.9526716,6.1];
ybs[3150]=['',2.1213067,0.3276197,6.15];
ybs[3151]=['',2.0971802,-1.1106569,4.82];
ybs[3152]=['',2.1119608,-0.5678238,5.82];
ybs[3153]=['',2.1035435,-0.9690788,6.28];
ybs[3154]=['',2.1100593,-0.7222142,5.52];
ybs[3155]=['8 Cnc',2.12249,0.2277116,5.12];
ybs[3156]=['',2.1254368,0.4792372,6.21];
ybs[3157]=['ζ Pup',2.1138021,-0.6994161,2.25];
ybs[3158]=['',2.1132072,-0.7508199,6.29];
ybs[3159]=['28 Lyn',2.1330347,0.7537727,6.26];
ybs[3160]=['14 Pup',2.1196146,-0.3455566,6.13];
ybs[3161]=['μ1 Cnc',2.1282347,0.393812,5.99];
ybs[3162]=['',2.1171976,-0.571519,5.31];
ybs[3163]=['',2.0898037,-1.2795439,6.34];
ybs[3164]=['',2.1252762,-0.0112593,6.41];
ybs[3165]=['27 Lyn',2.1392991,0.8976873,4.84];
ybs[3166]=['',2.1277232,-0.1626088,6.23];
ybs[3167]=['',2.1470137,1.015332,5.93];
ybs[3168]=['μ2 Cnc',2.1345441,0.3754053,5.3];
ybs[3169]=['',2.1236079,-0.5871383,6.14];
ybs[3170]=['',2.1179312,-0.8842072,5.95];
ybs[3171]=['',2.1209948,-0.8211777,6.19];
ybs[3172]=['',2.1192695,-0.9281488,5.53];
ybs[3173]=['',2.1426557,0.7392719,6.27];
ybs[3174]=['',2.1608058,1.1937856,5.32];
ybs[3175]=['',2.1309676,-0.3600077,5.38];
ybs[3176]=['12 Cnc',2.1383508,0.2368038,6.27];
ybs[3177]=['ρ Pup',2.131875,-0.4254498,2.81];
ybs[3178]=['',2.1164738,-1.097931,6.3];
ybs[3179]=['',2.1269466,-0.7913008,5.05];
ybs[3180]=['ζ Mon',2.1372635,-0.0533508,4.34];
ybs[3181]=['',2.1385141,-0.1991905,6.32];
ybs[3182]=['',2.137196,-0.3566751,6.36];
ybs[3183]=['ψ Cnc',2.1464386,0.4438955,5.73];
ybs[3184]=['16 Pup',2.1385572,-0.337164,4.4];
ybs[3185]=['',2.1510399,0.6746923,6.58];
ybs[3186]=['',2.1406301,-0.284876,5.68];
ybs[3187]=['',2.1359583,-0.6589353,6.37];
ybs[3188]=['',2.1384539,-0.5305031,6.65];
ybs[3189]=['',2.221272,1.4372666,6.32];
ybs[3190]=['',2.1483037,0.2540383,6.23];
ybs[3191]=['',2.1384527,-0.6200824,6.2];
ybs[3192]=['',2.1632305,0.9839562,5.85];
ybs[3193]=['',2.1494222,0.1701146,6.07];
ybs[3194]=['18 Pup',2.1459089,-0.2421308,5.54];
ybs[3195]=['',2.137552,-0.8509788,5.7];
ybs[3196]=['',2.139805,-0.7713668,5.21];
ybs[3197]=['',2.1407628,-0.745499,6.26];
ybs[3198]=['γ1 Vel',2.1390919,-0.8276186,4.27];
ybs[3199]=['γ2 Vel',2.1392892,-0.827459,1.78];
ybs[3200]=['ζ1 Cnc',2.1537897,0.3067073,5.63];
ybs[3201]=['ζ1 Cnc',2.1537897,0.3067073,6.02];
ybs[3202]=['ζ2 Cnc',2.1538333,0.3067072,6.2];
ybs[3203]=['19 Pup',2.1485962,-0.2269129,4.72];
ybs[3204]=['',2.1499968,-0.1369535,5.36];
ybs[3205]=['',2.1400584,-0.8379471,5.23];
ybs[3206]=['',2.1543362,0.2431083,6.54];
ybs[3207]=['15 Cnc',2.1583734,0.5162934,5.64];
ybs[3208]=['',2.1927547,1.3208338,5.54];
ybs[3209]=['',2.1324592,-1.114806,6.28];
ybs[3210]=['',2.1385959,-0.9801553,5.66];
ybs[3211]=['',2.1464622,-0.6521641,6.44];
ybs[3212]=['',2.1355138,-1.0712027,4.76];
ybs[3213]=['',2.1724536,1.0525008,6.45];
ybs[3214]=['',2.1571705,0.2869153,6.01];
ybs[3215]=['ε Vol',2.1293092,-1.1988577,4.35];
ybs[3216]=['',2.1604755,0.4025128,6.56];
ybs[3217]=['',2.1477783,-0.69277,4.45];
ybs[3218]=['',2.14789,-0.7515638,4.75];
ybs[3219]=['',2.1464185,-0.8471132,5.82];
ybs[3220]=['',2.1623974,0.3071799,6.47];
ybs[3221]=['20 Pup',2.1574889,-0.2768714,4.99];
ybs[3222]=['',2.1544331,-0.52335,6.52];
ybs[3223]=['',2.1629451,0.2264136,6.38];
ybs[3224]=['',2.1501191,-0.8153942,5.76];
ybs[3225]=['',2.1544273,-0.6632143,6.43];
ybs[3226]=['',2.1523971,-0.8087664,6.03];
ybs[3227]=['29 Lyn',2.1810423,1.0383564,5.64];
ybs[3228]=['',2.196169,1.262362,5.98];
ybs[3229]=['',2.1573108,-0.627882,4.78];
ybs[3230]=['',2.1582676,-0.5872078,6.37];
ybs[3231]=['',2.160514,-0.5622831,6.06];
ybs[3232]=['',2.1593847,-0.6352649,5.08];
ybs[3233]=['',2.1594129,-0.6355897,6.11];
ybs[3234]=['',2.160506,-0.6207469,5.78];
ybs[3235]=['',2.1594853,-0.7055245,4.44];
ybs[3236]=['',2.1570988,-0.821478,5.13];
ybs[3237]=['',2.1877723,1.0895884,5.71];
ybs[3238]=['',2.1822639,0.9436257,6.27];
ybs[3239]=['',2.1566911,-0.8774007,5.51];
ybs[3240]=['',2.1725594,0.2033222,7.13];
ybs[3241]=['β Cnc',2.1722511,0.1589769,3.52];
ybs[3242]=['',2.1606569,-0.8012828,5.83];
ybs[3243]=['',2.167959,-0.541092,6.21];
ybs[3244]=['',2.1766604,0.1533928,6.29];
ybs[3245]=['',2.1681784,-0.6279567,6.16];
ybs[3246]=['30 Lyn',2.1921104,1.0064338,5.89];
ybs[3247]=['',2.1728603,-0.3734526,6.6];
ybs[3248]=['',2.1646235,-0.8818379,6.44];
ybs[3249]=['21 Pup',2.1751559,-0.2855748,6.16];
ybs[3250]=['',2.1918841,0.933673,6.49];
ybs[3251]=['',2.1797122,-0.2218257,5.98];
ybs[3252]=['',2.1626083,-1.0994149,5.16];
ybs[3253]=['',2.1771573,-0.5250094,6.45];
ybs[3254]=['χ Cnc',2.188423,0.4736675,5.14];
ybs[3255]=['',2.2024937,1.0568153,6.41];
ybs[3256]=['',2.1894025,0.3607425,5.83];
ybs[3257]=['',2.1835241,-0.1787914,6.32];
ybs[3258]=['',2.1782836,-0.6201032,5.58];
ybs[3259]=['',2.1778325,-0.6536564,6.7];
ybs[3260]=['λ Cnc',2.1903391,0.4178907,5.98];
ybs[3261]=['',2.1865392,0.0675322,6.05];
ybs[3262]=['',2.1793778,-0.6411851,4.45];
ybs[3263]=['',2.188061,-0.0172454,6.18];
ybs[3264]=['',2.1882014,-0.0943845,6.13];
ybs[3265]=['',2.1835684,-0.6050793,6.43];
ybs[3266]=['',2.1748185,-1.0340079,6.42];
ybs[3267]=['31 Lyn',2.2013178,0.7523773,4.25];
ybs[3268]=['',2.1883093,-0.4014858,6.13];
ybs[3269]=['',2.2063242,0.927454,5.51];
ybs[3270]=['',2.1929174,-0.0293462,6.5];
ybs[3271]=['',2.1923355,-0.3518292,5.58];
ybs[3272]=['',2.1754588,-1.1465209,5.07];
ybs[3273]=['',2.1948597,-0.3083268,5.75];
ybs[3274]=['',2.1919129,-0.5782902,4.83];
ybs[3275]=['',2.1915949,-0.6381545,5.2];
ybs[3276]=['20 Cnc',2.2024524,0.3185578,5.95];
ybs[3277]=['',2.197853,-0.1092386,6.15];
ybs[3278]=['',2.1916563,-0.6928951,6.16];
ybs[3279]=['',2.2094299,0.7317134,6.02];
ybs[3280]=['',2.1995374,-0.1330511,5.96];
ybs[3281]=['22 Pup',2.1988217,-0.2292418,6.11];
ybs[3282]=['21 Cnc',2.2046043,0.1841582,6.08];
ybs[3283]=['',2.1985385,-0.4612541,5.9];
ybs[3284]=['',2.2106526,0.6096491,6.06];
ybs[3285]=['',2.1892685,-1.0131986,5.97];
ybs[3286]=['',2.1959785,-0.847705,4.82];
ybs[3287]=['',2.2070813,-0.0837357,6.01];
ybs[3288]=['',2.1999578,-0.6696109,6.32];
ybs[3289]=['1 Hya',2.207012,-0.0668785,5.61];
ybs[3290]=['',2.1880557,-1.1202392,6.12];
ybs[3291]=['25 Cnc',2.2131703,0.2960904,6.14];
ybs[3292]=['',2.1973961,-0.9111266,5.85];
ybs[3293]=['κ1 Vol',2.1804902,-1.2495357,5.37];
ybs[3294]=['κ2 Vol',2.1813478,-1.2493677,5.65];
ybs[3295]=['',2.2342616,1.1731063,5.88];
ybs[3296]=['φ1 Cnc',2.2163422,0.4854093,5.57];
ybs[3297]=['',2.2116077,0.0352728,5.73];
ybs[3298]=['',2.2131946,0.1306037,5.13];
ybs[3299]=['ε Car',2.1948183,-1.0400295,1.86];
ybs[3300]=['',2.2077794,-0.4055181,5.68];
ybs[3301]=['',2.2222993,0.7953596,6.32];
ybs[3302]=['φ2 Cnc',2.2176897,0.468666,6.32];
ybs[3303]=['φ2 Cnc',2.2177043,0.4686806,6.3];
ybs[3304]=['24 Cnc',2.2170827,0.4267744,7.02];
ybs[3305]=['24 Cnc',2.2171047,0.4267937,7.81];
ybs[3306]=['',2.2117051,-0.0695976,3.9];
ybs[3307]=['',2.2083697,-0.4210964,5.28];
ybs[3308]=['',2.209608,-0.368734,6.01];
ybs[3309]=['',2.2112202,-0.3057935,6.44];
ybs[3310]=['α Cha',2.1722467,-1.3438524,4.07];
ybs[3311]=['27 Cnc',2.2169393,0.2194339,5.5];
ybs[3312]=['',2.2124897,-0.261993,5.98];
ybs[3313]=['2 Hya',2.2151642,-0.0710199,5.59];
ybs[3314]=['',2.2069774,-0.7478784,5.98];
ybs[3315]=['ο UMa',2.235278,1.0582703,3.36];
ybs[3316]=['',2.2159425,-0.220194,5.54];
ybs[3317]=['',2.2187196,-0.1133026,6.59];
ybs[3318]=['',2.2109984,-0.737133,5.47];
ybs[3319]=['',2.213049,-0.6831381,6.53];
ybs[3320]=['',2.2130926,-0.6831528,7.25];
ybs[3321]=['28 Cnc',2.2255702,0.4199612,6.1];
ybs[3322]=['',2.2087815,-0.9042398,5.17];
ybs[3323]=['',2.2159476,-0.5113301,6.73];
ybs[3324]=['',2.248307,1.2083796,6.31];
ybs[3325]=['29 Cnc',2.2252377,0.2465823,5.95];
ybs[3326]=['η Vol',2.1896687,-1.2824537,5.29];
ybs[3327]=['',2.2193803,-0.365228,6.56];
ybs[3328]=['',2.2177008,-0.5542298,6.33];
ybs[3329]=['',2.224086,-0.0453756,6.39];
ybs[3330]=['',2.2231845,-0.1553105,6.43];
ybs[3331]=['',2.2206414,-0.4575343,6.62];
ybs[3332]=['θ Cha',2.1811341,-1.3537264,4.35];
ybs[3333]=['',2.2126338,-0.9230871,6.05];
ybs[3334]=['',2.2254222,-0.1715851,6];
ybs[3335]=['',2.2206593,-0.6142895,5.75];
ybs[3336]=['',2.2238495,-0.4041186,6.51];
ybs[3337]=['',2.2251972,-0.3670958,6.67];
ybs[3338]=['',2.2086723,-1.1289134,5.97];
ybs[3339]=['β Vol',2.2078355,-1.1557223,3.77];
ybs[3340]=['',2.2380089,0.6489449,6.18];
ybs[3341]=['',2.2169539,-0.9615655,6.53];
ybs[3342]=['',2.2178029,-0.9280033,5.09];
ybs[3343]=['',2.2444201,0.9255492,6.24];
ybs[3344]=['',2.2671959,1.3026571,6.31];
ybs[3345]=['',2.2274397,-0.4784911,6.7];
ybs[3346]=['2 UMa',2.25485,1.1354991,5.47];
ybs[3347]=['υ1 Cnc',2.2381912,0.4188272,5.75];
ybs[3348]=['',2.2251237,-0.7721925,5.79];
ybs[3349]=['θ Cnc',2.2383441,0.3143396,5.35];
ybs[3350]=['',2.2246532,-0.8379666,5.33];
ybs[3351]=['',2.2265435,-0.7820467,4.99];
ybs[3352]=['',2.2449559,0.6620316,5.9];
ybs[3353]=['',2.2394302,0.1698241,6.83];
ybs[3354]=['',2.2316681,-0.5627454,5.65];
ybs[3355]=['',2.2277463,-0.8100955,5.99];
ybs[3356]=['',2.23152,-0.6423614,6.69];
ybs[3357]=['32 Lyn',2.2468217,0.6344518,6.24];
ybs[3358]=['η Cnc',2.2432854,0.3552875,5.33];
ybs[3359]=['',2.2367157,-0.343158,5.42];
ybs[3360]=['',2.2263289,-0.9647152,6.36];
ybs[3361]=['υ2 Cnc',2.244703,0.418878,6.36];
ybs[3362]=['',2.213639,-1.2247864,5.53];
ybs[3363]=['',2.2317571,-0.7822698,6.3];
ybs[3364]=['34 Cnc',2.2427322,0.1742103,6.46];
ybs[3365]=['',2.2354006,-0.6832628,6.31];
ybs[3366]=['',2.2414189,-0.2637881,6.38];
ybs[3367]=['',2.2338427,-0.836893,6.39];
ybs[3368]=['',2.2475804,0.2298968,6.28];
ybs[3369]=['33 Lyn',2.2527775,0.634145,5.78];
ybs[3370]=['',2.2471807,0.0815345,5.87];
ybs[3371]=['',2.2793084,1.2835426,6.15];
ybs[3372]=['',2.2494689,0.1460252,6.03];
ybs[3373]=['',2.2433568,-0.4309413,6.19];
ybs[3374]=['',2.2346644,-0.9508215,6.34];
ybs[3375]=['',2.2482722,-0.0390409,5.81];
ybs[3376]=['',2.2420995,-0.55127,6.38];
ybs[3377]=['',2.2424707,-0.6059529,6.36];
ybs[3378]=['',2.2373406,-0.9301976,5.69];
ybs[3379]=['35 Cnc',2.254658,0.3404116,6.58];
ybs[3380]=['',2.2438325,-0.6711824,6.49];
ybs[3381]=['',2.2451467,-0.6795237,5.96];
ybs[3382]=['',2.244075,-0.8212817,6.24];
ybs[3383]=['π1 UMa',2.2747975,1.1332952,5.64];
ybs[3384]=['',2.2544885,0.0463864,6.33];
ybs[3385]=['',2.1939201,-1.4136143,5.69];
ybs[3386]=['',2.2580174,0.2657682,6.32];
ybs[3387]=['',2.2565122,0.1140386,5.99];
ybs[3388]=['',2.2565341,0.1140822,7.25];
ybs[3389]=['',2.2493433,-0.5704433,6.43];
ybs[3390]=['3 Hya',2.2543805,-0.1408148,5.72];
ybs[3391]=['',2.248925,-0.6579323,6.3];
ybs[3392]=['',2.269753,0.9305056,5.66];
ybs[3393]=['',2.2739291,1.0446092,6.48];
ybs[3394]=['',2.2537468,-0.4700076,5.96];
ybs[3395]=['π2 UMa',2.2791067,1.1211912,4.6];
ybs[3396]=['',2.2519518,-0.6991035,6.47];
ybs[3397]=['',2.2724542,0.9221861,6.42];
ybs[3398]=['36 Cnc',2.2620502,0.1670089,5.88];
ybs[3399]=['',2.2491871,-0.8731811,5.01];
ybs[3400]=['',2.273713,0.9184603,5.91];
ybs[3401]=['',2.2682394,0.570979,5.94];
ybs[3402]=['δ Hya',2.2643642,0.0980301,4.16];
ybs[3403]=['',2.2631279,-0.0876225,6.19];
ybs[3404]=['37 Cnc',2.2663661,0.1655902,6.53];
ybs[3405]=['',2.2540954,-0.8910943,5.8];
ybs[3406]=['',2.2510451,-1.013946,4.86];
ybs[3407]=['',2.2507137,-1.0177124,5.26];
ybs[3408]=['',2.2669414,-0.1178044,6.51];
ybs[3409]=['',2.2362903,-1.281786,6.12];
ybs[3410]=['σ Hya',2.2690873,0.0567928,4.44];
ybs[3411]=['',2.2622476,-0.59049,6.48];
ybs[3412]=['η Pyx',2.2642066,-0.4597536,5.27];
ybs[3413]=['',2.2612112,-0.7022186,6.55];
ybs[3414]=['34 Lyn',2.2807121,0.7984073,5.37];
ybs[3415]=['',2.2768684,0.5559535,6.1];
ybs[3416]=['',2.2720881,0.1383962,6.45];
ybs[3417]=['',2.2679534,-0.3459987,6.33];
ybs[3418]=['',2.262417,-0.7518173,4.14];
ybs[3419]=['39 Cnc',2.2755329,0.347665,6.39];
ybs[3420]=['',2.276662,0.3417676,6.44];
ybs[3421]=['ε Cnc',2.2770136,0.3395852,6.3];
ybs[3422]=['',2.2698616,-0.3970531,5.05];
ybs[3423]=['6 Hya',2.2741121,-0.2192697,4.98];
ybs[3424]=['',2.259089,-1.098513,5.47];
ybs[3425]=['ζ Pyx',2.2720991,-0.5174706,4.89];
ybs[3426]=['',2.2702933,-0.6404406,6.13];
ybs[3427]=['',2.2665057,-0.928128,6.47];
ybs[3428]=['',2.2894296,0.8170182,6.22];
ybs[3429]=['',2.2785965,-0.1595293,6.63];
ybs[3430]=['β Pyx',2.273567,-0.6177815,3.97];
ybs[3431]=['',2.274268,-0.7042786,5.2];
ybs[3432]=['',2.2693092,-0.9342273,5.48];
ybs[3433]=['9 Hya',2.2814011,-0.2798121,4.88];
ybs[3434]=['',2.2718065,-0.9275171,5.19];
ybs[3435]=['',2.2677925,-1.0542602,6.36];
ybs[3436]=['',2.2751422,-0.7902767,5.71];
ybs[3437]=['',2.2752141,-0.8157151,3.84];
ybs[3438]=['',2.2834623,-0.2104,6.45];
ybs[3439]=['',2.2732756,-0.9251975,3.62];
ybs[3440]=['',2.2732529,-0.9268265,5.61];
ybs[3441]=['γ Cnc',2.2894455,0.3731361,4.66];
ybs[3442]=['45 Cnc',2.2887945,0.2197614,5.64];
ybs[3443]=['',2.2939376,0.6427724,6.33];
ybs[3444]=['',2.2777543,-0.8273796,4.77];
ybs[3445]=['',2.2770718,-0.8554006,5.9];
ybs[3446]=['η Hya',2.2885841,0.0577563,4.3];
ybs[3447]=['',2.2747162,-1.005893,6.34];
ybs[3448]=['',2.2810715,-0.7941176,5.23];
ybs[3449]=['',2.273974,-1.0445655,4.33];
ybs[3450]=['',2.291973,0.0740884,6.37];
ybs[3451]=['',2.2902058,-0.1278143,4.62];
ybs[3452]=['θ Vol',2.2652755,-1.2300071,5.2];
ybs[3453]=['δ Cnc',2.2954337,0.3152774,3.94];
ybs[3454]=['',2.2822816,-0.8410401,5.51];
ybs[3455]=['',2.2859712,-0.6288868,6.42];
ybs[3456]=['46 Cnc',2.2988361,0.5341991,6.13];
ybs[3457]=['49 Cnc',2.2954491,0.1743854,5.66];
ybs[3458]=['',2.2820993,-0.9283211,5.52];
ybs[3459]=['',2.282579,-0.9285644,4.86];
ybs[3460]=['α Pyx',2.2889019,-0.5807743,3.68];
ybs[3461]=['10 Hya',2.2964944,0.0975696,6.13];
ybs[3462]=['',2.3169755,1.1626675,6.2];
ybs[3463]=['',2.2819948,-0.9749942,6.29];
ybs[3464]=['',2.2976539,-0.0469702,6.41];
ybs[3465]=['',2.2951812,-0.3710207,6.11];
ybs[3466]=['ι Cnc',2.3044674,0.5004605,6.57];
ybs[3467]=['ι Cnc',2.304598,0.5003681,4.02];
ybs[3468]=['',2.288297,-0.8711337,5.16];
ybs[3469]=['',2.2919534,-0.7459368,4.07];
ybs[3470]=['',2.3007038,-0.0373425,5.7];
ybs[3471]=['',2.2942782,-0.6499136,5.76];
ybs[3472]=['',2.3007414,-0.1936806,6.25];
ybs[3473]=['50 Cnc',2.3050381,0.2097695,5.87];
ybs[3474]=['ε Hya',2.3041739,0.1104422,3.38];
ybs[3475]=['',2.2989621,-0.4446757,6.1];
ybs[3476]=['12 Hya',2.3017937,-0.2380381,4.32];
ybs[3477]=['δ Vel',2.2923968,-0.9564104,1.96];
ybs[3478]=['',2.3059811,-0.0347047,5.29];
ybs[3479]=['',2.2988619,-0.8051596,3.91];
ybs[3480]=['',2.3007637,-0.7193604,6.21];
ybs[3481]=['',2.2936616,-1.0265172,6.21];
ybs[3482]=['',2.3029381,-0.6058691,6.37];
ybs[3483]=['',2.2869413,-1.1920801,6.32];
ybs[3484]=['ρ Hya',2.3113831,0.1002871,4.36];
ybs[3485]=['',2.3094658,-0.1160676,6.09];
ybs[3486]=['',2.3009753,-0.8029137,5.46];
ybs[3487]=['',2.2900413,-1.1504397,6.05];
ybs[3488]=['',2.3044779,-0.8071572,5.75];
ybs[3489]=['',2.3063155,-0.7300458,6.36];
ybs[3490]=['',2.3009623,-0.992404,4.49];
ybs[3491]=['',2.3215262,0.5793195,6.25];
ybs[3492]=['14 Hya',2.315148,-0.0617008,5.31];
ybs[3493]=['',2.3083357,-0.7427322,6.43];
ybs[3494]=['η Cha',2.2708305,-1.3797074,5.47];
ybs[3495]=['',2.3069927,-0.9240119,6.3];
ybs[3496]=['',2.3219136,0.3270649,6.16];
ybs[3497]=['5 UMa',2.3360048,1.0798033,5.73];
ybs[3498]=['',2.3344336,1.0290845,6.25];
ybs[3499]=['',2.316266,-0.3689782,6.47];
ybs[3500]=['35 Lyn',2.3281917,0.7615449,5.15];
ybs[3501]=['',2.3293605,0.789221,5.99];
ybs[3502]=['54 Cnc',2.3229897,0.2662964,6.38];
ybs[3503]=['',2.3290531,0.731451,5.99];
ybs[3504]=['',2.3162915,-0.5737396,5.21];
ybs[3505]=['',2.317209,-0.5158398,5.87];
ybs[3506]=['',2.3150363,-0.7053356,5.48];
ybs[3507]=['',2.318655,-0.5010943,6.17];
ybs[3508]=['',2.3160587,-0.6847619,6.39];
ybs[3509]=['γ Pyx',2.3194404,-0.4852471,4.01];
ybs[3510]=['σ1 Cnc',2.3303832,0.5651475,5.66];
ybs[3511]=['',2.3153704,-0.7923847,4.93];
ybs[3512]=['53 Cnc',2.3297762,0.4915828,6.23];
ybs[3513]=['ρ1 Cnc',2.3303022,0.4928327,5.95];
ybs[3514]=['15 Hya',2.324692,-0.1268911,5.54];
ybs[3515]=['',2.2790126,-1.3815791,6.05];
ybs[3516]=['',2.3179846,-0.7362235,6];
ybs[3517]=['',2.3286874,0.0915691,6.33];
ybs[3518]=['',2.3186292,-0.8137029,5.1];
ybs[3519]=['',2.3364086,0.618617,6.14];
ybs[3520]=['',2.3285904,-0.2326019,6.13];
ybs[3521]=['',2.3228337,-0.7434654,6.55];
ybs[3522]=['6 UMa',2.3505213,1.1258845,5.58];
ybs[3523]=['57 Cnc',2.3375769,0.5320659,5.39];
ybs[3524]=['',2.3275577,-0.5690226,6.5];
ybs[3525]=['',2.3282873,-0.6394724,6.42];
ybs[3526]=['',2.3288665,-0.6774974,5.82];
ybs[3527]=['',2.3223124,-1.0075193,5.59];
ybs[3528]=['',2.3164803,-1.1673727,5.35];
ybs[3529]=['',2.3366302,-0.0964944,6];
ybs[3530]=['',2.327621,-0.8456582,5.91];
ybs[3531]=['ρ2 Cnc',2.3436484,0.4857702,5.22];
ybs[3532]=['',2.3420516,0.2990903,6.64];
ybs[3533]=['',2.3275011,-0.9114571,6.39];
ybs[3534]=['',2.2905504,-1.3891874,5.79];
ybs[3535]=['',2.3117314,-1.2678579,6.11];
ybs[3536]=['',2.3495918,0.7947616,5.74];
ybs[3537]=['',2.3478751,0.699988,5.89];
ybs[3538]=['ζ Hya',2.34175,0.1021155,3.11];
ybs[3539]=['',2.3333445,-0.707583,6.47];
ybs[3540]=['',2.3287715,-0.9903534,6.03];
ybs[3541]=['60 Cnc',2.3442502,0.2012558,5.41];
ybs[3542]=['',2.3329303,-0.8310356,5.33];
ybs[3543]=['17 Hya',2.3417611,-0.140762,6.91];
ybs[3544]=['17 Hya',2.3417683,-0.1407765,6.67];
ybs[3545]=['',2.3401921,-0.3200243,5.75];
ybs[3546]=['σ2 Cnc',2.3494346,0.5727264,5.45];
ybs[3547]=['δ Pyx',2.3412476,-0.4847952,4.89];
ybs[3548]=['',2.3470284,0.0722807,6.14];
ybs[3549]=['',2.3497025,0.29755,6.17];
ybs[3550]=['',2.3431659,-0.4173655,6.39];
ybs[3551]=['',2.3316331,-1.0550183,5.78];
ybs[3552]=['ο1 Cnc',2.350131,0.2657649,5.2];
ybs[3553]=['',2.3395409,-0.7877772,6.26];
ybs[3554]=['61 Cnc',2.3538369,0.5260022,6.29];
ybs[3555]=['',2.3461869,-0.293297,5.96];
ybs[3556]=['ο2 Cnc',2.3516147,0.2702759,5.67];
ybs[3557]=['',2.3561628,0.6231939,6.51];
ybs[3558]=['',2.3519167,0.1621764,6.19];
ybs[3559]=['',2.33664,-1.0181276,6.38];
ybs[3560]=['ι UMa',2.3600806,0.8368016,3.14];
ybs[3561]=['',2.3382709,-0.9609802,5.71];
ybs[3562]=['',2.3370013,-1.0600987,3.84];
ybs[3563]=['α Cnc',2.3554236,0.2052802,4.25];
ybs[3564]=['',2.3535874,0.0252329,6.59];
ybs[3565]=['',2.3433839,-0.9218595,4.69];
ybs[3566]=['σ3 Cnc',2.3607466,0.564126,5.2];
ybs[3567]=['ρ UMa',2.3767931,1.1786515,4.76];
ybs[3568]=['',2.3586457,0.3148282,6.38];
ybs[3569]=['',2.3556494,-0.2832483,5.86];
ybs[3570]=['',2.3659582,0.7275531,3.97];
ybs[3571]=['',2.3652003,0.6546286,6.44];
ybs[3572]=['',2.4437239,1.4674293,6.33];
ybs[3573]=['',2.3456013,-1.035412,4.92];
ybs[3574]=['',2.3507168,-0.8494357,5.87];
ybs[3575]=['',2.3596244,-0.3369288,6.18];
ybs[3576]=['',2.3575236,-0.5044432,6.25];
ybs[3577]=['',2.37037,0.6914271,6.36];
ybs[3578]=['66 Cnc',2.3688396,0.5612133,5.82];
ybs[3579]=['',2.3549215,-0.8260795,5.18];
ybs[3580]=['67 Cnc',2.3704725,0.4852934,6.07];
ybs[3581]=['',2.3684836,0.0967519,6.07];
ybs[3582]=['',2.3605908,-0.7217036,4.45];
ybs[3583]=['',2.3814424,0.9457137,5.75];
ybs[3584]=['',2.3617116,-0.7552061,6.07];
ybs[3585]=['κ UMa',2.3792669,0.8213234,3.6];
ybs[3586]=['ν Cnc',2.3743739,0.425073,5.45];
ybs[3587]=['',2.3702393,-0.0101284,5.67];
ybs[3588]=['',2.3660251,-0.4670686,6.2];
ybs[3589]=['',2.356243,-1.0328847,5.16];
ybs[3590]=['',2.3738691,0.125667,5.85];
ybs[3591]=['',2.3660503,-0.7323689,5.55];
ybs[3592]=['70 Cnc',2.3807169,0.4851989,6.38];
ybs[3593]=['',2.3694961,-0.6894056,6.27];
ybs[3594]=['',2.3870916,0.8452846,5.95];
ybs[3595]=['',2.3619724,-1.0657116,5.79];
ybs[3596]=['',2.3671274,-0.9125571,5.23];
ybs[3597]=['',2.3841797,0.5633599,6.46];
ybs[3598]=['',2.3746116,-0.4468467,6.74];
ybs[3599]=['',2.3937259,1.0340173,6.45];
ybs[3600]=['σ1 UMa',2.4020627,1.1654085,5.14];
ybs[3601]=['',2.3623545,-1.2004527,5.88];
ybs[3602]=['',2.3728932,-0.9363275,6.4];
ybs[3603]=['',2.3914571,0.6693818,4.56];
ybs[3604]=['ω Hya',2.3878813,0.0871448,4.97];
ybs[3605]=['',2.3780253,-0.8237279,3.75];
ybs[3606]=['α Vol',2.3685326,-1.1605333,4];
ybs[3607]=['σ2 UMa',2.4107721,1.1699566,4.8];
ybs[3608]=['',2.3948594,0.399354,6.4];
ybs[3609]=['',2.3922496,0.0237919,6.17];
ybs[3610]=['15 UMa',2.4024119,0.898919,4.48];
ybs[3611]=['',2.3979082,0.566193,6.5];
ybs[3612]=['τ Cnc',2.3975066,0.5158166,5.43];
ybs[3613]=['',2.3799945,-1.0114372,6.44];
ybs[3614]=['κ Cnc',2.3957785,0.1844488,5.24];
ybs[3615]=['τ UMa',2.4125048,1.1067531,4.67];
ybs[3616]=['',2.4013762,0.5896039,5.93];
ybs[3617]=['75 Cnc',2.4008263,0.4630149,5.98];
ybs[3618]=['ξ Cnc',2.4031564,0.3830119,5.14];
ybs[3619]=['κ Pyx',2.3960045,-0.4530583,4.58];
ybs[3620]=['',2.3879202,-0.9756853,6.11];
ybs[3621]=['19 Hya',2.3994038,-0.1516643,5.6];
ybs[3622]=['',2.3912587,-0.8955555,6.73];
ybs[3623]=['',2.3849687,-1.1274619,6.37];
ybs[3624]=['',2.4277072,1.2488389,6.55];
ybs[3625]=['λ Vel',2.395029,-0.7597843,2.21];
ybs[3626]=['',2.4046447,0.2000794,6.48];
ybs[3627]=['',2.4014255,-0.2174375,5.77];
ybs[3628]=['',2.3989276,-0.4689358,6.15];
ybs[3629]=['',2.4007209,-0.3216472,5.73];
ybs[3630]=['',2.4090778,0.5386425,5.95];
ybs[3631]=['79 Cnc',2.4074714,0.3821469,6.01];
ybs[3632]=['20 Hya',2.4032752,-0.1551322,5.46];
ybs[3633]=['',2.3816194,-1.2328608,4.71];
ybs[3634]=['',2.3788644,-1.2688784,4.48];
ybs[3635]=['ε Pyx',2.4041028,-0.5317325,5.59];
ybs[3636]=['',2.4359959,1.2713458,5.96];
ybs[3637]=['',2.4062922,-0.4062708,6.53];
ybs[3638]=['',2.4023433,-0.8643804,6.48];
ybs[3639]=['16 UMa',2.4271294,1.0702475,5.13];
ybs[3640]=['',2.4138641,0.0936673,6.35];
ybs[3641]=['π1 Cnc',2.4157365,0.2599557,6.51];
ybs[3642]=['',2.4150773,0.0657208,6.14];
ybs[3643]=['36 Lyn',2.4233771,0.7525054,5.32];
ybs[3644]=['',2.4133669,-0.3464366,5.73];
ybs[3645]=['',2.4083979,-0.7848611,5];
ybs[3646]=['21 Hya',2.415725,-0.1258644,6.11];
ybs[3647]=['',2.4113484,-0.6869671,6];
ybs[3648]=['',2.4217211,0.3696791,6.48];
ybs[3649]=['',2.4103978,-0.8148112,5.79];
ybs[3650]=['',2.4068626,-1.0309315,3.44];
ybs[3651]=['17 UMa',2.4331202,0.9885225,5.27];
ybs[3652]=['',2.4147455,-0.7629769,5.57];
ybs[3653]=['18 UMa',2.4344443,0.9410569,4.83];
ybs[3654]=['',2.4078279,-1.0894067,3.97];
ybs[3655]=['',2.4292401,0.6026739,5.97];
ybs[3656]=['θ Hya',2.4244147,0.0386,3.88];
ybs[3657]=['',2.4539261,1.2899983,6.5];
ybs[3658]=['',2.4190174,-0.6757655,6.31];
ybs[3659]=['',2.4183111,-0.7395951,6.29];
ybs[3660]=['π2 Cnc',2.4285424,0.2589804,5.34];
ybs[3661]=['',2.4191909,-0.8279976,5.92];
ybs[3662]=['',2.421827,-0.7722771,5.85];
ybs[3663]=['',2.4153926,-1.0387555,5.54];
ybs[3664]=['',2.4230627,-0.7562511,5.25];
ybs[3665]=['',2.4285134,-0.2640275,6.35];
ybs[3666]=['',2.4397605,0.8153022,5.97];
ybs[3667]=['',2.4257113,-0.6580802,5.86];
ybs[3668]=['ζ Oct',2.3242533,-1.4924605,5.42];
ybs[3669]=['',2.421771,-0.9716621,5.27];
ybs[3670]=['',2.4265874,-0.7968889,6.25];
ybs[3671]=['23 Hya',2.4343451,-0.1126872,5.24];
ybs[3672]=['',2.4285451,-0.6749708,4.94];
ybs[3673]=['24 Hya',2.43425,-0.1544296,5.47];
ybs[3674]=['',2.4292077,-0.6547842,4.62];
ybs[3675]=['β Car',2.4149755,-1.2185735,1.68];
ybs[3676]=['',2.4431926,0.6154029,5.75];
ybs[3677]=['',2.4359914,-0.2561657,5.84];
ybs[3678]=['',2.4302371,-0.7854288,6.04];
ybs[3679]=['',2.4399002,0.1989186,6.41];
ybs[3680]=['38 Lyn',2.4450454,0.6405038,3.82];
ybs[3681]=['',2.4258486,-1.0208676,6.02];
ybs[3682]=['',2.4316479,-0.7743869,5.12];
ybs[3683]=['',2.4272039,-1.0067228,6.32];
ybs[3684]=['',2.4343704,-0.6894904,5.33];
ybs[3685]=['',2.408184,-1.3397922,6.14];
ybs[3686]=['',2.4299204,-1.0060871,4.34];
ybs[3687]=['',2.4540137,0.8929283,6.13];
ybs[3688]=['',2.4587452,0.9877461,5.47];
ybs[3689]=['ι Car',2.4336283,-1.0363551,2.25];
ybs[3690]=['',2.4367558,-0.9529328,6.33];
ybs[3691]=['',2.4544428,0.6646769,6.12];
ybs[3692]=['',2.4467201,-0.1992942,6.62];
ybs[3693]=['',2.4387238,-0.8928236,5.26];
ybs[3694]=['',2.4465558,-0.2781877,5.78];
ybs[3695]=['α Lyn',2.4545789,0.5984268,3.13];
ybs[3696]=['26 Hya',2.4476255,-0.2108293,4.79];
ybs[3697]=['',2.456258,0.5724091,6.16];
ybs[3698]=['',2.4413418,-0.9017192,5.87];
ybs[3699]=['27 Hya',2.4507913,-0.1686116,4.8];
ybs[3700]=['',2.4470235,-0.5970412,6.39];
ybs[3701]=['',2.4548382,0.26644,6.53];
ybs[3702]=['',2.4331159,-1.2006635,5.39];
ybs[3703]=['',2.4359505,-1.1720685,6.11];
ybs[3704]=['',2.4525501,-0.2744152,6.33];
ybs[3705]=['',2.4512113,-0.5561583,6.82];
ybs[3706]=['',2.4499187,-0.6577492,6.05];
ybs[3707]=['',2.4447416,-0.9650121,6.28];
ybs[3708]=['θ Pyx',2.4547168,-0.4550216,4.72];
ybs[3709]=['',2.4886884,1.3088302,6.29];
ybs[3710]=['',2.4319309,-1.3089607,5.29];
ybs[3711]=['',2.432141,-1.3061734,5.86];
ybs[3712]=['',2.4769825,1.1141105,6.28];
ybs[3713]=['',2.4650444,0.4376756,6.41];
ybs[3714]=['',2.4611062,-0.1735674,6.53];
ybs[3715]=['',2.4723107,0.8982726,6.31];
ybs[3716]=['',2.4556318,-0.7382807,5.58];
ybs[3717]=['',2.4691327,0.6367052,6.67];
ybs[3718]=['',2.4501554,-1.0909998,4.81];
ybs[3719]=['',2.45908,-0.696044,6.54];
ybs[3720]=['',2.4578455,-0.805523,5.75];
ybs[3721]=['κ Leo',2.4699807,0.455107,4.46];
ybs[3722]=['',2.4547331,-0.9707578,5.63];
ybs[3723]=['λ Pyx',2.4620967,-0.5050946,4.69];
ybs[3724]=['κ Vel',2.4559991,-0.9619604,2.5];
ybs[3725]=['',2.4641269,-0.6608394,6.48];
ybs[3726]=['',2.4735594,0.2876081,6.29];
ybs[3727]=['',2.4663496,-0.6899656,6.06];
ybs[3728]=['28 Hya',2.4723678,-0.0911803,5.59];
ybs[3729]=['',2.4644579,-0.9048376,6.08];
ybs[3730]=['',2.461387,-1.0543258,6.3];
ybs[3731]=['',2.4767026,-0.0274192,6.01];
ybs[3732]=['',2.4640313,-1.0778288,5.99];
ybs[3733]=['',2.4882371,0.7940092,5.41];
ybs[3734]=['29 Hya',2.4803054,-0.1628573,6.54];
ybs[3735]=['',2.4775675,-0.504308,6.1];
ybs[3736]=['',2.475931,-0.7087615,6.2];
ybs[3737]=['',2.4937872,0.9710507,6.45];
ybs[3738]=['α Hya',2.4818187,-0.1529984,1.98];
ybs[3739]=['',2.4802169,-0.3918495,4.69];
ybs[3740]=['',2.4827297,-0.1078393,5.38];
ybs[3741]=['',2.5324571,1.4174708,4.29];
ybs[3742]=['',2.4699361,-1.0830978,5.77];
ybs[3743]=['',2.4744461,-0.9335098,5.11];
ybs[3744]=['ω Leo',2.4860772,0.1561857,5.41];
ybs[3745]=['3 Leo',2.4861782,0.1410302,5.71];
ybs[3746]=['',2.481248,-0.6128779,6.65];
ybs[3747]=['23 UMa',2.5020467,1.0987344,3.67];
ybs[3748]=['',2.4883366,-0.0238243,6.27];
ybs[3749]=['τ1 Hya',2.4887856,-0.0502134,4.6];
ybs[3750]=['',2.4899345,-0.0403782,6.14];
ybs[3751]=['',2.4751761,-1.1351069,6.05];
ybs[3752]=['',2.4904572,-0.0760176,6.26];
ybs[3753]=['',2.4885749,-0.3640188,5.66];
ybs[3754]=['7 LMi',2.4966697,0.5855025,5.85];
ybs[3755]=['ε Ant',2.4882243,-0.6293572,4.51];
ybs[3756]=['',2.4882427,-0.6721615,6.19];
ybs[3757]=['',2.491213,-0.4093431,6.24];
ybs[3758]=['22 UMa',2.5182692,1.2582986,5.72];
ybs[3759]=['8 LMi',2.5002961,0.6107609,5.37];
ybs[3760]=['',2.491448,-0.4659697,5.48];
ybs[3761]=['24 UMa',2.5159365,1.2168451,4.56];
ybs[3762]=['',2.4938135,-0.2737685,5.85];
ybs[3763]=['λ Leo',2.5006822,0.3989644,4.31];
ybs[3764]=['',2.5241742,1.295156,6.46];
ybs[3765]=['θ UMa',2.506829,0.900026,3.17];
ybs[3766]=['',2.4844868,-1.0887528,5.92];
ybs[3767]=['',2.4755093,-1.2515655,5.47];
ybs[3768]=['',2.5078386,0.8609532,6.76];
ybs[3769]=['6 Leo',2.5013684,0.1676682,5.07];
ybs[3770]=['ζ1 Ant',2.4949778,-0.5585064,7];
ybs[3771]=['ζ1 Ant',2.4950288,-0.5584726,6.18];
ybs[3772]=['ξ Leo',2.5013437,0.1953124,4.97];
ybs[3773]=['',2.4826738,-1.1660492,5.91];
ybs[3774]=['',2.49113,-0.9010369,5.45];
ybs[3775]=['',2.4994975,-0.1860737,6.14];
ybs[3776]=['ψ Vel',2.4943972,-0.7081727,3.6];
ybs[3777]=['τ2 Hya',2.5011889,-0.0225872,4.57];
ybs[3778]=['',2.5007316,-0.1829048,6.13];
ybs[3779]=['ζ2 Ant',2.4983859,-0.558172,5.93];
ybs[3780]=['',2.4982974,-0.6252459,5.87];
ybs[3781]=['9 LMi',2.5089103,0.6349021,6.18];
ybs[3782]=['',2.5077568,0.4932022,6.53];
ybs[3783]=['',2.4919279,-1.0204966,5.88];
ybs[3784]=['',2.5043554,0.0306263,6.11];
ybs[3785]=['ι Cha',2.4577498,-1.4118464,5.36];
ybs[3786]=['',2.50228,-0.3405056,5.74];
ybs[3787]=['',2.5129522,0.8166776,6.52];
ybs[3788]=['',2.5018653,-0.5015603,6.46];
ybs[3789]=['26 UMa',2.5154206,0.9065443,4.5];
ybs[3790]=['10 LMi',2.5120369,0.6333366,4.55];
ybs[3791]=['',2.5056023,-0.1503565,6.12];
ybs[3792]=['',2.5050149,-0.2378259,5.94];
ybs[3793]=['',2.4956849,-0.9973375,3.13];
ybs[3794]=['',2.5105713,0.4074297,6.25];
ybs[3795]=['',2.5069449,-0.1274025,6.24];
ybs[3796]=['',2.5317022,1.2735517,6.42];
ybs[3797]=['',2.5014767,-0.7113729,5.35];
ybs[3798]=['',2.5060231,-0.370453,5.01];
ybs[3799]=['',2.5158235,0.6895988,4.81];
ybs[3800]=['',2.5069698,-0.4009637,5.91];
ybs[3801]=['',2.5171882,0.6955649,6.76];
ybs[3802]=['',2.5050737,-0.6848392,6.43];
ybs[3803]=['',2.4959867,-1.1663739,6.27];
ybs[3804]=['33 Hya',2.5122575,-0.105157,5.56];
ybs[3805]=['11 LMi',2.5182747,0.6230789,5.41];
ybs[3806]=['',2.4995777,-1.0977775,6.1];
ybs[3807]=['',2.5072687,-0.8572134,5.12];
ybs[3808]=['7 Leo',2.518593,0.2490442,6.36];
ybs[3809]=['',2.5089043,-0.8964906,5.01];
ybs[3810]=['',2.5227099,0.541939,5.56];
ybs[3811]=['',2.4948419,-1.2774005,5.47];
ybs[3812]=['',2.5163416,-0.3437252,6.31];
ybs[3813]=['',2.514234,-0.6271689,6.49];
ybs[3814]=['',2.5371135,1.172169,5.94];
ybs[3815]=['',2.5095961,-1.0356674,4.08];
ybs[3816]=['8 Leo',2.5237126,0.2849569,5.69];
ybs[3817]=['10 Leo',2.5242045,0.1173706,5];
ybs[3818]=['',2.5205629,-0.4330775,6.53];
ybs[3819]=['42 Lyn',2.5302174,0.7003708,5.25];
ybs[3820]=['',2.5224743,-0.4434455,5.7];
ybs[3821]=['',2.5190026,-0.8528035,6.17];
ybs[3822]=['34 Hya',2.5266314,-0.1664284,6.4];
ybs[3823]=['',2.5229631,-0.5635589,5.63];
ybs[3824]=['',2.5295755,0.0791986,4.68];
ybs[3825]=['',2.5241645,-0.6319291,5.98];
ybs[3826]=['',2.5207394,-0.8633457,4.35];
ybs[3827]=['',2.5202766,-0.9259832,6.19];
ybs[3828]=['',2.5495266,1.2064529,5.69];
ybs[3829]=['27 UMa',2.5532365,1.2590701,5.17];
ybs[3830]=['',2.5221229,-0.9386298,5.45];
ybs[3831]=['',2.5161292,-1.1355292,6.56];
ybs[3832]=['',2.5262927,-0.7557682,5.5];
ybs[3833]=['',2.5664511,1.3617182,6.23];
ybs[3834]=['',2.5293014,-0.6933427,6.7];
ybs[3835]=['ι Hya',2.5355525,-0.0218981,3.91];
ybs[3836]=['37 Hya',2.5350346,-0.1864386,6.31];
ybs[3837]=['',2.5747473,1.3791949,6.17];
ybs[3838]=['',2.5374086,-0.189913,6.37];
ybs[3839]=['κ Hya',2.5371962,-0.2520999,5.06];
ybs[3840]=['',2.5439453,0.5439414,5.89];
ybs[3841]=['43 Lyn',2.546069,0.6919377,5.62];
ybs[3842]=['ο Leo',2.5414583,0.1706911,3.52];
ybs[3843]=['13 Leo',2.5440216,0.4502994,6.24];
ybs[3844]=['',2.5495701,0.8433115,6.39];
ybs[3845]=['',2.551659,0.9468507,6.47];
ybs[3846]=['',2.5308628,-1.0723247,4.52];
ybs[3847]=['13 LMi',2.5489751,0.6105238,6.14];
ybs[3848]=['',2.541215,-0.4137133,4.77];
ybs[3849]=['',2.5591451,1.1322002,6.17];
ybs[3850]=['ζ Cha',2.5005685,-1.4146042,5.11];
ybs[3851]=['15 Leo',2.5524827,0.5211775,5.64];
ybs[3852]=['',2.5453838,-0.4193718,4.94];
ybs[3853]=['',2.5371234,-1.0139614,5.32];
ybs[3854]=['',2.5386204,-1.0013292,5.8];
ybs[3855]=['28 UMa',2.5646834,1.1089702,6.34];
ybs[3856]=['ψ Leo',2.5528191,0.2427483,5.35];
ybs[3857]=['',2.5469874,-0.62159,6.41];
ybs[3858]=['',2.5421033,-0.9656321,6];
ybs[3859]=['',2.5562913,0.3272517,6.5];
ybs[3860]=['',2.5667018,0.9950797,5.2];
ybs[3861]=['θ Ant',2.5538404,-0.4866461,4.79];
ybs[3862]=['',2.5496344,-0.8960761,6.15];
ybs[3863]=['ε Leo',2.5623143,0.4129492,2.98];
ybs[3864]=['',2.5537118,-0.6926242,6.82];
ybs[3865]=['',2.550538,-0.9425613,5.56];
ybs[3866]=['',2.563266,0.1150976,5.79];
ybs[3867]=['18 Leo',2.5643523,0.2041323,5.63];
ybs[3868]=['',2.5588414,-0.5291225,6.45];
ybs[3869]=['',2.564142,0.0291728,5.65];
ybs[3870]=['19 Leo',2.5688971,0.1999085,6.45];
ybs[3871]=['',2.575016,0.8012157,5.09];
ybs[3872]=['',2.5694463,0.197474,6.02];
ybs[3873]=['',2.5589271,-1.0000617,6.46];
ybs[3874]=['',2.5565614,-1.0929492,3.69];
ybs[3875]=['',2.5843487,1.1428044,6.31];
ybs[3876]=['',2.5633003,-0.7831129,5.55];
ybs[3877]=['',2.5598822,-1.0281387,6.22];
ybs[3878]=['υ UMa',2.5862593,1.0284,3.8];
ybs[3879]=['20 Leo',2.5796121,0.3676403,6.09];
ybs[3880]=['υ Car',2.5643853,-1.1377129,3.01];
ybs[3881]=['υ Car',2.564429,-1.1377226,6.26];
ybs[3882]=['',2.5765473,-0.6510326,5.97];
ybs[3883]=['4 Sex',2.5821237,0.0737962,6.24];
ybs[3884]=['φ UMa',2.5907572,0.9415786,4.59];
ybs[3885]=['',2.5720794,-0.9865768,6.06];
ybs[3886]=['23 Leo',2.5846405,0.2260295,6.46];
ybs[3887]=['',2.578269,-0.6350166,6.37];
ybs[3888]=['',2.5783178,-0.8001977,5.08];
ybs[3889]=['6 Sex',2.5851262,-0.0760781,6.01];
ybs[3890]=['22 Leo',2.5886255,0.4237559,5.32];
ybs[3891]=['',2.5856359,-0.1099091,6.42];
ybs[3892]=['ν Cha',2.5582801,-1.3419829,5.45];
ybs[3893]=['υ1 Hya',2.5859545,-0.2611424,4.12];
ybs[3894]=['',2.5815338,-0.8211749,5.73];
ybs[3895]=['μ Leo',2.5925046,0.45188,3.88];
ybs[3896]=['7 Sex',2.589506,0.0408099,6.02];
ybs[3897]=['',2.5894394,-0.0007047,6.35];
ybs[3898]=['',2.5881629,-0.2906073,6.08];
ybs[3899]=['γ Sex',2.5905975,-0.1434839,5.05];
ybs[3900]=['',2.5843233,-0.8082533,5.62];
ybs[3901]=['',2.6041217,1.0646372,6.27];
ybs[3902]=['',2.5858372,-0.8144317,4.58];
ybs[3903]=['',2.5829365,-1.0391931,5.79];
ybs[3904]=['',2.5814229,-1.0971265,5.57];
ybs[3905]=['',2.5961776,0.1019608,5.95];
ybs[3906]=['',2.5921375,-0.4790645,6.3];
ybs[3907]=['31 UMa',2.6062269,0.8674797,5.27];
ybs[3908]=['',2.6202704,1.2699275,5.83];
ybs[3909]=['',2.5975843,-0.4546414,4.88];
ybs[3910]=['',2.591117,-0.9684736,6.48];
ybs[3911]=['',2.5990859,-0.3945312,6.24];
ybs[3912]=['',2.6132573,1.0000874,5.93];
ybs[3913]=['',2.6006635,-0.333815,4.94];
ybs[3914]=['',2.5950362,-0.8947139,5.93];
ybs[3915]=['',2.597321,-0.792387,5.71];
ybs[3916]=['',2.6081004,0.1538649,5.85];
ybs[3917]=['',2.5995393,-0.8789582,5.72];
ybs[3918]=['19 LMi',2.6144186,0.7145011,5.14];
ybs[3919]=['',2.6157377,0.7905767,6.3];
ybs[3920]=['',2.6053555,-0.7145695,6.41];
ybs[3921]=['',2.608793,-0.4654375,6.28];
ybs[3922]=['',2.607797,-0.5853114,5.84];
ybs[3923]=['',2.6093142,-0.4815775,6.32];
ybs[3924]=['',2.671057,1.4625378,6.37];
ybs[3925]=['',2.6060858,-0.8980291,6.37];
ybs[3926]=['',2.6172733,0.4824269,6.3];
ybs[3927]=['ν Leo',2.6159843,0.2151456,5.26];
ybs[3928]=['',2.6154732,0.1430544,6.04];
ybs[3929]=['',2.6246755,0.9894902,5.48];
ybs[3930]=['φ Vel',2.6080759,-0.9544349,3.54];
ybs[3931]=['',2.6095854,-0.9207711,6.12];
ybs[3932]=['',2.6224164,0.5153445,5.73];
ybs[3933]=['',2.6121111,-0.8470436,6.05];
ybs[3934]=['',2.6030815,-1.2480234,6.35];
ybs[3935]=['12 Sex',2.6223129,0.057011,6.7];
ybs[3936]=['',2.6190217,-0.4200711,6.21];
ybs[3937]=['η Ant',2.6176825,-0.6284766,5.23];
ybs[3938]=['',2.608942,-1.1276021,6.58];
ybs[3939]=['',2.6071696,-1.2081038,6.2];
ybs[3940]=['π Leo',2.6245696,0.1383311,4.7];
ybs[3941]=['20 LMi',2.6286173,0.5551016,5.36];
ybs[3942]=['',2.6362212,0.3810056,5.66];
ybs[3943]=['',2.624155,-0.9959736,6.52];
ybs[3944]=['',2.6451272,0.9384978,5.74];
ybs[3945]=['',2.6292031,-0.9334579,6.2];
ybs[3946]=['',2.6350978,-0.535757,6.54];
ybs[3947]=['',2.6302259,-1.0030154,6.2];
ybs[3948]=['',2.6475238,0.9119515,6.14];
ybs[3949]=['',2.6393491,-0.1691795,6.12];
ybs[3950]=['',2.6301355,-1.0566165,5.94];
ybs[3951]=['13 Sex',2.6416,0.0537839,6.45];
ybs[3952]=['',2.6390272,-0.4439426,6.7];
ybs[3953]=['',2.6407604,-0.318019,5.86];
ybs[3954]=['',2.6368274,-0.8160349,6.12];
ybs[3955]=['',2.6419349,-0.4259496,5.7];
ybs[3956]=['',2.6345277,-1.0523939,6.19];
ybs[3957]=['',2.6335877,-1.0869117,6.42];
ybs[3958]=['',2.6416652,-0.6997966,6.43];
ybs[3959]=['',2.6485822,0.2729265,6.37];
ybs[3960]=['υ2 Hya',2.6455695,-0.2301131,4.6];
ybs[3961]=['',2.636891,-1.0821595,6.14];
ybs[3962]=['',2.6455547,-0.6371097,6.27];
ybs[3963]=['14 Sex',2.653209,0.0958381,6.21];
ybs[3964]=['21 LMi',2.6566794,0.6130338,4.48];
ybs[3965]=['η Leo',2.655806,0.2904638,3.52];
ybs[3966]=['',2.6492447,-0.8288579,5.08];
ybs[3967]=['',2.6543603,-0.3012792,5.6];
ybs[3968]=['',2.648709,-0.9129484,6.52];
ybs[3969]=['',2.6602163,0.5494902,6.24];
ybs[3970]=['31 Leo',2.6581598,0.1723847,4.37];
ybs[3971]=['α Sex',2.658102,-0.0085915,4.49];
ybs[3972]=['α Leo',2.6602349,0.2067605,1.35];
ybs[3973]=['μ1 Cha',2.6179008,-1.4369811,5.52];
ybs[3974]=['',2.6576528,-0.653699,6.36];
ybs[3975]=['',2.6614892,-0.1920827,6.53];
ybs[3976]=['',2.6606506,-0.2745827,6.27];
ybs[3977]=['',2.672318,0.7075553,6.32];
ybs[3978]=['',2.6666159,-0.2132263,6.24];
ybs[3979]=['17 Sex',2.6674889,-0.1488682,5.91];
ybs[3980]=['',2.6610663,-0.9063832,4.86];
ybs[3981]=['',2.6672854,-0.2257982,5.31];
ybs[3982]=['',2.6642596,-0.6279288,6.13];
ybs[3983]=['',2.6732253,0.6506663,5.85];
ybs[3984]=['λ Hya',2.6694337,-0.2177382,3.61];
ybs[3985]=['',2.6590034,-1.1508003,5.28];
ybs[3986]=['18 Sex',2.6710023,-0.1490465,5.65];
ybs[3987]=['μ2 Cha',2.6336301,-1.4256737,6.6];
ybs[3988]=['34 Leo',2.6745022,0.2309662,6.44];
ybs[3989]=['',2.6622975,-1.076346,5.6];
ybs[3990]=['',2.6726239,-0.1298206,6.25];
ybs[3991]=['',2.6688898,-0.7301812,5.98];
ybs[3992]=['',2.6621306,-1.200861,5.81];
ybs[3993]=['',2.6754558,-0.5013997,6.28];
ybs[3994]=['19 Sex',2.6794302,0.0784142,5.77];
ybs[3995]=['',2.6782089,-0.3364204,6.44];
ybs[3996]=['',2.6843581,0.4714767,6.04];
ybs[3997]=['',2.6722386,-1.0288645,6.4];
ybs[3998]=['',2.6912854,1.0448059,6.25];
ybs[3999]=['',2.6731163,-1.01547,5.72];
ybs[4000]=['',2.6761211,-0.912547,6.16];
ybs[4001]=['',2.6810631,-0.473873,6.25];
ybs[4002]=['',2.687125,0.3673115,6.02];
ybs[4003]=['',2.6813071,-0.5786465,6.38];
ybs[4004]=['22 LMi',2.6900192,0.5470824,6.46];
ybs[4005]=['',2.6826321,-0.7062994,5.9];
ybs[4006]=['',2.7054158,1.273217,6.4];
ybs[4007]=['',2.6805269,-0.8963202,5.28];
ybs[4008]=['',2.6784433,-1.0478953,6.1];
ybs[4009]=['',2.6834127,-0.7056845,6.35];
ybs[4010]=['',2.680881,-0.9054447,5.78];
ybs[4011]=['',2.7042833,1.2380884,6.66];
ybs[4012]=['',2.6797592,-1.0782799,6.41];
ybs[4013]=['',2.6868217,-0.7373029,3.85];
ybs[4014]=['23 LMi',2.6949128,0.5094218,5.35];
ybs[4015]=['',2.6799832,-1.1605582,5.16];
ybs[4016]=['32 UMa',2.7043423,1.1342021,5.82];
ybs[4017]=['24 LMi',2.6958931,0.4984592,6.49];
ybs[4018]=['',2.6947868,0.3074824,6.55];
ybs[4019]=['',2.6896674,-0.6394994,6.19];
ybs[4020]=['35 Leo',2.6960826,0.4080605,5.97];
ybs[4021]=['ζ Leo',2.6967418,0.4065618,3.44];
ybs[4022]=['',2.6968205,0.4406684,5.84];
ybs[4023]=['λ UMa',2.6990235,0.7468505,3.45];
ybs[4024]=['',2.6937414,-0.1976782,6.08];
ybs[4025]=['37 Leo',2.6964966,0.2374589,5.41];
ybs[4026]=['',2.6902348,-0.7545949,5.6];
ybs[4027]=['ω Car',2.680477,-1.2245253,3.32];
ybs[4028]=['',2.688634,-0.9616186,6.16];
ybs[4029]=['39 Leo',2.6991474,0.4011294,5.82];
ybs[4030]=['',2.696212,-0.362915,6.57];
ybs[4031]=['',2.7032877,0.4763344,6.52];
ybs[4032]=['ε Sex',2.7002522,-0.1429785,5.24];
ybs[4033]=['',2.6917142,-1.0476522,6.22];
ybs[4034]=['',2.7073736,0.8139741,6.43];
ybs[4035]=['',2.6949219,-0.8958407,6.3];
ybs[4036]=['',2.7094535,0.8425276,6];
ybs[4037]=['',2.7178807,1.1977041,5.96];
ybs[4038]=['',2.7068899,0.4291437,6.4];
ybs[4039]=['',2.7019997,-0.5081567,5.34];
ybs[4040]=['',2.6961251,-1.0725956,3.4];
ybs[4041]=['',2.7131729,0.9364615,6.45];
ybs[4042]=['',2.7143839,0.944101,6];
ybs[4043]=['',2.7040013,-0.6445175,6.3];
ybs[4044]=['40 Leo',2.7099363,0.3376708,4.79];
ybs[4045]=['',2.7073728,-0.2208128,6];
ybs[4046]=['',2.7031578,-0.7294027,5.96];
ybs[4047]=['γ1 Leo',2.7109754,0.3441421,2.61];
ybs[4048]=['γ2 Leo',2.7109971,0.3441226,3.8];
ybs[4049]=['',2.7086269,-0.0912718,6.37];
ybs[4050]=['',2.7105351,-0.1602676,6.32];
ybs[4051]=['',2.7032329,-0.9814577,5.81];
ybs[4052]=['',2.7618533,1.4682741,5.5];
ybs[4053]=['',2.7076115,-0.9626029,4.57];
ybs[4054]=['23 Sex',2.7152867,0.0377984,6.66];
ybs[4055]=['',2.7045824,-1.1309711,5.67];
ybs[4056]=['',2.7108539,-0.8346685,5.65];
ybs[4057]=['',2.7210854,0.7174195,5.76];
ybs[4058]=['',2.7153474,-0.3160625,6.51];
ybs[4059]=['μ UMa',2.7217613,0.7221311,3.05];
ybs[4060]=['42 Leo',2.7190223,0.2592044,6.12];
ybs[4061]=['',2.7167552,-0.4159986,6.5];
ybs[4062]=['',2.7308713,1.1421704,4.97];
ybs[4063]=['',2.7173107,-0.3953606,6.51];
ybs[4064]=['',2.7132463,-0.9802993,4.5];
ybs[4065]=['27 LMi',2.7249381,0.5896332,5.9];
ybs[4066]=['',2.7200489,-0.3489133,6.13];
ybs[4067]=['43 Leo',2.7239688,0.1120149,6.07];
ybs[4068]=['',2.7274161,0.5147174,6.39];
ybs[4069]=['',2.7249867,0.0972077,6.54];
ybs[4070]=['',2.7200186,-0.7290995,4.83];
ybs[4071]=['28 LMi',2.7294574,0.5863225,5.5];
ybs[4072]=['25 Sex',2.7256864,-0.0732826,5.97];
ybs[4073]=['',2.7242167,-0.5286039,6.27];
ybs[4074]=['',2.7658051,1.4387105,5.26];
ybs[4075]=['',2.729184,0.0391521,6.32];
ybs[4076]=['',2.7252019,-0.6655745,5.33];
ybs[4077]=['',2.7258967,-0.7343994,6.27];
ybs[4078]=['44 Leo',2.7338071,0.1511396,5.61];
ybs[4079]=['',2.7213992,-1.1698263,4.99];
ybs[4080]=['30 LMi',2.7371784,0.5876677,4.74];
ybs[4081]=['',2.7259849,-1.0136622,6.35];
ybs[4082]=['',2.7356593,-0.1254,5.57];
ybs[4083]=['',2.7329344,-0.7433898,6.18];
ybs[4084]=['μ Hya',2.7370256,-0.2960364,3.81];
ybs[4085]=['',2.7309461,-1.0245317,5.95];
ybs[4086]=['',2.7441413,0.7238794,6.02];
ybs[4087]=['',2.7416459,0.3357833,6.15];
ybs[4088]=['',2.7469731,0.8492594,6.44];
ybs[4089]=['',2.7367371,-0.7481203,6.13];
ybs[4090]=['β LMi',2.7458247,0.6384682,4.21];
ybs[4091]=['45 Leo',2.7442773,0.1681953,6.04];
ybs[4092]=['',2.7265395,-1.2942737,4];
ybs[4093]=['',2.7492224,0.7869054,6.35];
ybs[4094]=['α Ant',2.741381,-0.5444252,4.25];
ybs[4095]=['',2.7280625,-1.293228,6.19];
ybs[4096]=['35 UMa',2.7559469,1.1431893,6.32];
ybs[4097]=['',2.7391476,-0.9599816,5.58];
ybs[4098]=['',2.7581656,1.1193005,6.12];
ybs[4099]=['',2.7487847,-0.0675157,6.05];
ybs[4100]=['',2.7415835,-1.0081792,4.66];
ybs[4101]=['',2.7447023,-0.8644831,6.1];
ybs[4102]=['36 UMa',2.7584334,0.9748401,4.84];
ybs[4103]=['32 LMi',2.7555658,0.6771717,5.77];
ybs[4104]=['',2.7435843,-1.0273894,3.82];
ybs[4105]=['',2.741052,-1.1489543,6.01];
ybs[4106]=['δ Sex',2.7520519,-0.0500073,5.21];
ybs[4107]=['',2.7515983,-0.5199273,5.58];
ybs[4108]=['δ Ant',2.7520442,-0.5363968,5.56];
ybs[4109]=['β Sex',2.7556342,-0.0133198,5.09];
ybs[4110]=['',2.7475695,-1.1222131,5.29];
ybs[4111]=['',2.7856659,1.4026661,6.52];
ybs[4112]=['',2.7585203,-0.1355051,6.2];
ybs[4113]=['',2.7585039,-0.2393669,5.58];
ybs[4114]=['33 LMi',2.7630346,0.5629185,5.9];
ybs[4115]=['',2.7576623,-0.4644361,6.51];
ybs[4116]=['',2.77992,1.3192193,4.84];
ybs[4117]=['46 Leo',2.7641788,0.2445306,5.46];
ybs[4118]=['',2.7555703,-1.0730699,6.43];
ybs[4119]=['',2.7528661,-1.1713103,6.19];
ybs[4120]=['',2.761792,-0.495046,6.05];
ybs[4121]=['',2.7718031,0.9314907,6.45];
ybs[4122]=['',2.7692154,0.7033443,4.75];
ybs[4123]=['ρ Leo',2.7667837,0.1602191,3.85];
ybs[4124]=['',2.7591018,-0.93972,4.89];
ybs[4125]=['',2.7620263,-0.7887708,5.74];
ybs[4126]=['',2.7619606,-0.7888193,6.09];
ybs[4127]=['34 LMi',2.7703175,0.6084508,5.58];
ybs[4128]=['',2.7530402,-1.2587178,4.74];
ybs[4129]=['',2.7646712,-0.7809579,5.91];
ybs[4130]=['',2.7615497,-1.0788202,3.32];
ybs[4131]=['37 UMa',2.7782102,0.9940602,5.16];
ybs[4132]=['',2.7558899,-1.2801638,4.93];
ybs[4133]=['',2.7662938,-0.8225758,5.02];
ybs[4134]=['',2.7651123,-1.0261434,6];
ybs[4135]=['44 Hya',2.7715089,-0.4166505,5.08];
ybs[4136]=['48 Leo',2.7754203,0.119143,5.08];
ybs[4137]=['',2.7678906,-1.0178264,6.14];
ybs[4138]=['49 Leo',2.7764785,0.1487545,5.67];
ybs[4139]=['',2.7756537,-0.4067202,6.1];
ybs[4140]=['35 LMi',2.7827272,0.6317984,6.28];
ybs[4141]=['',2.7711983,-1.066655,6.23];
ybs[4142]=['',2.7787313,-0.3263165,6.49];
ybs[4143]=['',2.776419,-0.6927223,5.38];
ybs[4144]=['',2.7761341,-0.7643146,6.08];
ybs[4145]=['',2.7816592,-0.1869398,6.57];
ybs[4146]=['φ2 Hya',2.7815183,-0.2874901,6.03];
ybs[4147]=['',2.7804604,-0.4677915,6.29];
ybs[4148]=['',2.7827252,-0.2156853,5.7];
ybs[4149]=['',2.7774119,-1.0067954,4.45];
ybs[4150]=['',2.7855847,-0.207281,6.52];
ybs[4151]=['',2.7560949,-1.4319994,7.07];
ybs[4152]=['',2.7854776,-0.4806675,4.89];
ybs[4153]=['',2.7871312,-0.235833,4.82];
ybs[4154]=['',2.7805863,-1.0418259,5.08];
ybs[4155]=['',2.7951636,0.9344526,5.52];
ybs[4156]=['37 LMi',2.7929314,0.5558535,4.71];
ybs[4157]=['',2.7852948,-0.8439288,3.84];
ybs[4158]=['38 LMi',2.7948292,0.6594179,5.85];
ybs[4159]=['',2.785487,-1.0273195,5.45];
ybs[4160]=['',2.7744591,-1.3340673,6.3];
ybs[4161]=['φ3 Hya',2.791574,-0.2967875,4.91];
ybs[4162]=['',2.7927625,-0.219417,6.04];
ybs[4163]=['',2.7881693,-1.0015442,5.91];
ybs[4164]=['γ Cha',2.7739553,-1.3741855,4.11];
ybs[4165]=['',2.7921771,-0.7484261,6.11];
ybs[4166]=['',2.8078702,1.1923152,5.75];
ybs[4167]=['',2.7911778,-1.0351734,4.66];
ybs[4168]=['38 UMa',2.8082246,1.1447207,5.12];
ybs[4169]=['',2.7922414,-1.0287844,5.92];
ybs[4170]=['',2.7937856,-0.9726975,4.28];
ybs[4171]=['',2.8134161,1.203355,5];
ybs[4172]=['33 Sex',2.804104,-0.0326419,6.26];
ybs[4173]=['',2.801181,-0.6260518,6.37];
ybs[4174]=['',2.8080385,0.550969,6.02];
ybs[4175]=['',2.7970814,-1.1384582,5.52];
ybs[4176]=['',2.7918869,-1.3023943,6.07];
ybs[4177]=['39 UMa',2.8154369,0.9960613,5.8];
ybs[4178]=['',2.8022829,-1.0438024,6.42];
ybs[4179]=['40 LMi',2.8116106,0.4572177,5.51];
ybs[4180]=['',2.8088127,-0.2461577,6.24];
ybs[4181]=['',2.8142921,0.8041582,5.18];
ybs[4182]=['41 LMi',2.8132444,0.4024616,5.08];
ybs[4183]=['35 Sex',2.8126792,0.0806136,5.79];
ybs[4184]=['',2.8093704,-0.5732475,5.64];
ybs[4185]=['',2.8219256,1.1742935,6];
ybs[4186]=['',2.8060961,-1.1273971,4.82];
ybs[4187]=['',2.8167879,0.3425989,6.27];
ybs[4188]=['',2.8083668,-1.0357593,5.38];
ybs[4189]=['θ Car',2.8092724,-1.1261439,2.76];
ybs[4190]=['',2.8120531,-1.0593387,4.57];
ybs[4191]=['36 Sex',2.8205381,0.0411679,6.28];
ybs[4192]=['41 UMa',2.8269812,0.9989614,6.34];
ybs[4193]=['42 LMi',2.8240406,0.5332465,5.24];
ybs[4194]=['',2.8132138,-1.1236066,5.77];
ybs[4195]=['',2.8143814,-1.1185848,4.82];
ybs[4196]=['',2.8016608,-1.3947258,5.97];
ybs[4197]=['',2.824683,0.1089708,6.37];
ybs[4198]=['51 Leo',2.8262214,0.3274559,5.49];
ybs[4199]=['52 Leo',2.8262131,0.2454835,5.48];
ybs[4200]=['η Car',2.8187812,-1.0439412,6.21];
ybs[4201]=['',2.8146484,-1.2389935,6.26];
ybs[4202]=['',2.8155787,-1.2389069,6.46];
ybs[4203]=['',2.8149515,-1.2666379,6.27];
ybs[4204]=['',2.8277357,-0.3041462,5.42];
ybs[4205]=['',2.8380961,1.1345024,6.39];
ybs[4206]=['μ Vel',2.8266766,-0.8648035,2.69];
ybs[4207]=['',2.8240765,-1.0599877,6.25];
ybs[4208]=['',2.8311116,-0.268636,6.67];
ybs[4209]=['',2.8238114,-1.1282591,5.34];
ybs[4210]=['',2.8247935,-1.1238675,5.23];
ybs[4211]=['',2.827234,-0.9928627,5.23];
ybs[4212]=['',2.826356,-1.125963,4.85];
ybs[4213]=['43 LMi',2.8374768,0.5111339,6.15];
ybs[4214]=['',2.8358394,-0.0364573,5.93];
ybs[4215]=['',2.8334947,-0.5553325,5.88];
ybs[4216]=['',2.8302094,-1.0052664,6.36];
ybs[4217]=['53 Leo',2.8385282,0.1817797,5.34];
ybs[4218]=['',2.8320416,-1.0480526,6];
ybs[4219]=['40 Sex',2.8384826,-0.0725051,6.61];
ybs[4220]=['44 LMi',2.8415571,0.4859643,6.04];
ybs[4221]=['δ1 Cha',2.8163962,-1.4067171,5.47];
ybs[4222]=['ν Hya',2.8397947,-0.2849029,3.11];
ybs[4223]=['',2.8403146,-0.1742349,5.86];
ybs[4224]=['δ2 Cha',2.8186537,-1.4079502,4.45];
ybs[4225]=['43 UMa',2.8478629,0.9852697,5.67];
ybs[4226]=['42 UMa',2.8488934,1.0330523,5.58];
ybs[4227]=['41 Sex',2.8428439,-0.1575687,5.79];
ybs[4228]=['',2.8409393,-0.5966972,5.61];
ybs[4229]=['',2.8378166,-1.0376624,5.91];
ybs[4230]=['',2.8463566,-0.05625,5.95];
ybs[4231]=['',2.8534965,0.9151569,6.65];
ybs[4232]=['',2.8535743,0.9140805,6.44];
ybs[4233]=['',2.8587887,1.2168968,5.93];
ybs[4234]=['',2.8513736,0.0156154,6.38];
ybs[4235]=['',2.8529877,-0.0057952,6.31];
ybs[4236]=['44 UMa',2.8581822,0.9504045,5.1];
ybs[4237]=['46 LMi',2.8565476,0.5948819,3.83];
ybs[4238]=['ω UMa',2.8596248,0.7515232,4.71];
ybs[4239]=['',2.8565124,-0.0416447,6.12];
ybs[4240]=['',2.8515442,-1.0013158,5.25];
ybs[4241]=['',2.8566258,-0.3537727,5.24];
ybs[4242]=['',2.8569344,-0.2718588,6.38];
ybs[4243]=['',2.8578813,-0.0394446,5.45];
ybs[4244]=['48 LMi',2.8624762,0.4426125,6.2];
ybs[4245]=['',2.8602227,-0.2424086,5.66];
ybs[4246]=['',2.8637693,0.5917306,5.72];
ybs[4247]=['',2.8557599,-1.0294669,3.78];
ybs[4248]=['46 UMa',2.8671151,0.5825169,5.03];
ybs[4249]=['54 Leo',2.8664322,0.429675,4.5];
ybs[4250]=['54 Leo',2.8664685,0.4296604,6.3];
ybs[4251]=['',2.8640501,-0.3629601,6.44];
ybs[4252]=['',2.8558464,-1.2365846,5.99];
ybs[4253]=['',2.8629272,-0.7397081,6.11];
ybs[4254]=['',2.8694547,0.7308927,6.03];
ybs[4255]=['55 Leo',2.8665468,0.0105728,5.91];
ybs[4256]=['',2.8599973,-1.0813643,5.93];
ybs[4257]=['56 Leo',2.8679978,0.1056633,5.81];
ybs[4258]=['',2.8485995,-1.3908528,6.33];
ybs[4259]=['',2.8693131,0.3878192,6.14];
ybs[4260]=['50 LMi',2.8706268,0.4427671,6.35];
ybs[4261]=['',2.863539,-1.0585076,5.92];
ybs[4262]=['',2.8878253,1.3550407,6.2];
ybs[4263]=['ι Ant',2.870465,-0.6504685,4.6];
ybs[4264]=['',2.8719719,-0.8883095,5.91];
ybs[4265]=['',2.882994,0.9032161,6.17];
ybs[4266]=['',2.8746144,-1.044814,6.11];
ybs[4267]=['47 UMa',2.8834676,0.7033416,5.05];
ybs[4268]=['',2.8837407,0.6276426,6];
ybs[4269]=['',2.8709074,-1.3130303,6.13];
ybs[4270]=['',2.8869578,0.7922785,5.47];
ybs[4271]=['',2.8840241,0.202005,6.55];
ybs[4272]=['',2.8815061,-0.5911246,5.71];
ybs[4273]=['',2.8878925,0.896576,6.43];
ybs[4274]=['',2.8829749,-0.287729,5.89];
ybs[4275]=['',2.8873321,0.7466428,6.02];
ybs[4276]=['',2.8912482,1.1046028,6.39];
ybs[4277]=['α Crt',2.8840873,-0.3216764,4.08];
ybs[4278]=['49 UMa',2.8894266,0.6820789,5.08];
ybs[4279]=['',2.8859625,-0.2481021,5.88];
ybs[4280]=['',2.8807928,-1.0725394,6.16];
ybs[4281]=['58 Leo',2.8877532,0.0608347,4.84];
ybs[4282]=['',2.8846416,-0.7668812,5.81];
ybs[4283]=['',2.8853947,-0.7392812,4.39];
ybs[4284]=['59 Leo',2.8885939,0.1041862,4.99];
ybs[4285]=['β UMa',2.8941938,0.9817541,2.37];
ybs[4286]=['',2.8851211,-0.9066921,6.15];
ybs[4287]=['',2.8892615,-0.2779396,6.34];
ybs[4288]=['',2.8878636,-0.5580059,6.07];
ybs[4289]=['61 Leo',2.8932236,-0.0456725,4.74];
ybs[4290]=['60 Leo',2.8956535,0.3498953,4.42];
ybs[4291]=['α UMa',2.9026002,1.0754443,1.79];
ybs[4292]=['',2.8954761,-0.4706034,6.23];
ybs[4293]=['',2.8994158,-0.0154431,6.14];
ybs[4294]=['',2.8776199,-1.4257204,6.71];
ybs[4295]=['',2.8993297,-0.1995947,5.5];
ybs[4296]=['62 Leo',2.9010233,-0.002325,5.95];
ybs[4297]=['',2.899172,-0.5601313,6.46];
ybs[4298]=['',2.9008789,-0.2367857,6.34];
ybs[4299]=['51 UMa',2.9054378,0.6651253,6];
ybs[4300]=['χ Leo',2.907233,0.1257253,4.63];
ybs[4301]=['',2.9043726,-0.8344709,5.67];
ybs[4302]=['η Oct',2.8752237,-1.4787389,6.19];
ybs[4303]=['',2.9062615,-0.6272238,5.43];
ybs[4304]=['χ1 Hya',2.9082482,-0.478678,4.94];
ybs[4305]=['',2.9094522,-0.1958528,6.09];
ybs[4306]=['',2.9067401,-0.8643757,6.13];
ybs[4307]=['χ2 Hya',2.9109931,-0.4785777,5.71];
ybs[4308]=['',2.911189,-0.8961431,6.3];
ybs[4309]=['65 Leo',2.9154115,0.0318124,5.52];
ybs[4310]=['',2.913968,-0.5037121,6.77];
ybs[4311]=['',2.9127706,-0.8916788,6.32];
ybs[4312]=['64 Leo',2.9189314,0.4047535,6.46];
ybs[4313]=['',2.9126844,-1.026394,6.02];
ybs[4314]=['',2.9160815,-0.5710732,6.59];
ybs[4315]=['',2.9127616,-1.0918245,4.61];
ybs[4316]=['',2.9120448,-1.1339836,6.41];
ybs[4317]=['',2.9165271,-0.7465034,5.15];
ybs[4318]=['',2.9194624,-0.528969,6.54];
ybs[4319]=['',2.9135096,-1.2393733,5.57];
ybs[4320]=['',2.9286558,1.1707155,6.06];
ybs[4321]=['',2.9210233,-0.5254453,6.49];
ybs[4322]=['67 Leo',2.9239878,0.4280461,5.68];
ybs[4323]=['',2.926309,0.6313953,5.74];
ybs[4324]=['',2.923092,-0.4924208,5.44];
ybs[4325]=['ψ UMa',2.9279382,0.7743224,3.01];
ybs[4326]=['',2.9278218,0.7517882,5.89];
ybs[4327]=['',2.9218794,-1.0316301,3.91];
ybs[4328]=['',2.9216643,-1.0835051,5.13];
ybs[4329]=['',2.9281012,-0.5672446,5.81];
ybs[4330]=['',2.9396458,1.1892396,6.4];
ybs[4331]=['',2.9365723,0.2490029,6.3];
ybs[4332]=['',2.9320659,-1.0225645,6.88];
ybs[4333]=['β Crt',2.9359335,-0.400715,4.48];
ybs[4334]=['',2.9415498,0.9557523,6.63];
ybs[4335]=['',2.9403213,0.6227343,6.41];
ybs[4336]=['',2.9383985,-0.5684086,6.38];
ybs[4337]=['ψ Crt',2.9396805,-0.3252168,6.13];
ybs[4338]=['',2.939955,-0.3819257,6.4];
ybs[4339]=['',2.9339331,-1.2491287,6.35];
ybs[4340]=['',2.9394796,-0.859307,5.36];
ybs[4341]=['',2.9453288,0.714798,6.33];
ybs[4342]=['',2.9393996,-1.05507,4.6];
ybs[4343]=['',2.9412025,-0.8703955,6.11];
ybs[4344]=['',2.9425973,-0.7767738,5.8];
ybs[4345]=['',2.9399204,-1.1223042,5.23];
ybs[4346]=['69 Leo',2.9453105,-0.0035505,5.42];
ybs[4347]=['δ Leo',2.9470057,0.3558703,2.56];
ybs[4348]=['',2.9465554,0.138349,5.79];
ybs[4349]=['θ Leo',2.947534,0.26696,3.34];
ybs[4350]=['',2.9442196,-0.931401,5.76];
ybs[4351]=['',2.9434244,-1.0428885,5.74];
ybs[4352]=['72 Leo',2.951803,0.4007569,4.63];
ybs[4353]=['',2.9559625,0.9187253,6.5];
ybs[4354]=['',2.9498451,-0.765641,6.21];
ybs[4355]=['73 Leo',2.9546021,0.2299219,5.32];
ybs[4356]=['',2.9550271,0.2218447,6.67];
ybs[4357]=['',2.9586462,0.8611864,5.88];
ybs[4358]=['φ Leo',2.9579439,-0.066073,4.47];
ybs[4359]=['',2.9592626,-0.1268643,6.14];
ybs[4360]=['',2.9566413,-0.803096,6.31];
ybs[4361]=['75 Leo',2.9607289,0.0327502,5.18];
ybs[4362]=['',2.9599653,-0.6658175,6.27];
ybs[4363]=['',2.9619997,-0.6086201,6.45];
ybs[4364]=['ξ UMa',2.964865,0.5479455,4.87];
ybs[4365]=['ξ UMa',2.9648722,0.5479455,4.41];
ybs[4366]=['',2.962256,-0.6399876,6.68];
ybs[4367]=['ν UMa',2.9661737,0.5752594,3.48];
ybs[4368]=['',2.9654287,0.2068303,6.66];
ybs[4369]=['',2.9597582,-1.1860857,6.06];
ybs[4370]=['55 UMa',2.9690759,0.6641197,4.78];
ybs[4371]=['76 Leo',2.9678234,0.026464,5.91];
ybs[4372]=['δ Crt',2.9695565,-0.2602797,3.56];
ybs[4373]=['',2.9773873,1.1687783,6.21];
ybs[4374]=['',2.9684985,-1.1295214,5.99];
ybs[4375]=['',2.9639633,-1.3928221,6.35];
ybs[4376]=['σ Leo',2.9775397,0.1028862,4.05];
ybs[4377]=['',2.9692516,-1.3138285,6.27];
ybs[4378]=['',2.9810843,0.9937979,6.43];
ybs[4379]=['',2.9715627,-1.2588855,6.41];
ybs[4380]=['π Cen',2.9763729,-0.9533965,3.89];
ybs[4381]=['',2.9857777,1.1204294,6.02];
ybs[4382]=['56 UMa',2.9852297,0.7565672,4.99];
ybs[4383]=['',2.9825662,-0.7815664,6.12];
ybs[4384]=['',2.9869381,-0.0000532,6.05];
ybs[4385]=['λ Crt',2.9870941,-0.330124,5.09];
ybs[4386]=['',2.986278,-0.6335445,5];
ybs[4387]=['',2.9792738,-1.3568695,6.43];
ybs[4388]=['',2.9856429,-0.9933391,5.79];
ybs[4389]=['ι Leo',2.9897325,0.1814164,3.94];
ybs[4390]=['79 Leo',2.9901679,0.0222181,5.39];
ybs[4391]=['',2.9864124,-1.1360297,5.11];
ybs[4392]=['ε Crt',2.9925842,-0.1918863,4.83];
ybs[4393]=['',2.9912642,-0.7470702,6.12];
ybs[4394]=['',2.9943468,0.1971422,5.8];
ybs[4395]=['γ Crt',2.9937247,-0.3109958,4.08];
ybs[4396]=['',2.9896526,-1.2634691,5.59];
ybs[4397]=['',2.9990156,0.9724207,5.75];
ybs[4398]=['81 Leo',2.9971052,0.2848634,5.57];
ybs[4399]=['',2.9962372,-0.6317737,5.22];
ybs[4400]=['80 Leo',2.9980144,0.0650146,6.37];
ybs[4401]=['',2.9964905,-0.6611777,5.89];
ybs[4402]=['',3.0007977,0.5814662,6.32];
ybs[4403]=['',2.9967798,-1.1188904,5.17];
ybs[4404]=['83 Leo',3.0020304,0.0502312,6.5];
ybs[4405]=['',3.0006676,-1.069019,5.3];
ybs[4406]=['κ Crt',3.0036993,-0.2180216,5.94];
ybs[4407]=['',3.0017064,-0.9301735,5.81];
ybs[4408]=['τ Leo',3.0071849,0.0474903,4.95];
ybs[4409]=['',3.0069772,-0.0320288,6.25];
ybs[4410]=['',3.007105,-0.6189588,6.45];
ybs[4411]=['',3.012725,1.0758754,5.83];
ybs[4412]=['57 UMa',3.0123705,0.6841994,5.31];
ybs[4413]=['',3.009703,-0.7471638,5.08];
ybs[4414]=['',3.0154382,0.9878954,6.28];
ybs[4415]=['',3.0077421,-1.2672763,6.09];
ybs[4416]=['85 Leo',3.0149401,0.2666527,5.74];
ybs[4417]=['',3.0175297,0.9464287,6.41];
ybs[4418]=['',3.0144751,-0.4293263,5.76];
ybs[4419]=['',3.0260005,1.4135734,6.15];
ybs[4420]=['',3.0183067,0.8119653,6.35];
ybs[4421]=['58 UMa',3.0187131,0.751155,5.94];
ybs[4422]=['87 Leo',3.0175303,-0.0547844,4.77];
ybs[4423]=['86 Leo',3.0183801,0.3189486,5.52];
ybs[4424]=['λ Dra',3.0230583,1.2076932,3.84];
ybs[4425]=['',3.0203448,0.8341595,6.42];
ybs[4426]=['',3.0216094,0.849169,6.56];
ybs[4427]=['88 Leo',3.0238735,0.2483436,6.2];
ybs[4428]=['',3.0210917,-1.0718713,6.38];
ybs[4429]=['',3.0269157,1.0637266,5.48];
ybs[4430]=['',3.0238881,-0.3649845,6.24];
ybs[4431]=['ο1 Cen',3.0234011,-1.0398258,5.13];
ybs[4432]=['ο2 Cen',3.0235967,-1.0411105,5.15];
ybs[4433]=['',3.0259115,-0.5131054,5.81];
ybs[4434]=['',3.0259261,-0.5130666,5.64];
ybs[4435]=['',3.0264516,-0.4691815,6.16];
ybs[4436]=['',3.0283147,-0.1389802,5.95];
ybs[4437]=['',3.0281611,-0.7081127,5.64];
ybs[4438]=['',3.0256829,-1.1710752,5.9];
ybs[4439]=['',3.028667,-0.5449391,5.04];
ybs[4440]=['ξ Hya',3.0290988,-0.5583879,3.54];
ybs[4441]=['',3.0302443,-0.2865144,6.05];
ybs[4442]=['',3.0335529,0.6401867,6.4];
ybs[4443]=['',3.0317451,-0.7107413,5.39];
ybs[4444]=['',3.0344031,0.1900321,6.55];
ybs[4445]=['89 Leo',3.0352387,0.0510407,5.77];
ybs[4446]=['90 Leo',3.0367928,0.2907952,5.95];
ybs[4447]=['',3.0387003,0.9538162,5.63];
ybs[4448]=['',3.0356036,-0.5753823,5.98];
ybs[4449]=['',3.0383585,0.3544023,6.45];
ybs[4450]=['',3.0365926,-0.9494552,4.62];
ybs[4451]=['2 Dra',3.0431921,1.2075424,5.2];
ybs[4452]=['',3.0374621,-0.8599636,5.5];
ybs[4453]=['',3.0386755,-0.8291734,5.71];
ybs[4454]=['',3.0411917,0.1880669,6.56];
ybs[4455]=['',3.0437792,0.4825034,5.8];
ybs[4456]=['',3.0417643,-0.833872,5.25];
ybs[4457]=['λ Cen',3.040909,-1.1022695,3.13];
ybs[4458]=['θ Crt',3.0452853,-0.17345,4.7];
ybs[4459]=['',3.0447399,-0.5882758,5.74];
ybs[4460]=['',3.0451396,-0.6522907,6.31];
ybs[4461]=['υ Leo',3.0464862,-0.0167488,4.3];
ybs[4462]=['',3.0435403,-1.0679309,5.83];
ybs[4463]=['',3.0466516,-0.5781194,6.29];
ybs[4464]=['',3.05084,0.8810864,6.14];
ybs[4465]=['',3.0463335,-1.0719651,5.15];
ybs[4466]=['',3.0489306,-0.835716,5.44];
ybs[4467]=['59 UMa',3.0527838,0.759039,5.59];
ybs[4468]=['',3.0518253,0.1526876,6.17];
ybs[4469]=['π Cha',3.046928,-1.3270161,5.65];
ybs[4470]=['60 UMa',3.0537441,0.8150396,6.1];
ybs[4471]=['',3.055098,1.1206949,6.46];
ybs[4472]=['',3.0535635,0.5845059,6.27];
ybs[4473]=['ω Vir',3.0531165,0.1395973,5.36];
ybs[4474]=['',3.0528216,-0.0448888,6.22];
ybs[4475]=['',3.0496884,-1.1825665,5.96];
ybs[4476]=['',3.0545578,0.7849228,6.44];
ybs[4477]=['',3.0511958,-1.0814444,5.15];
ybs[4478]=['ι Crt',3.0539441,-0.2327883,5.48];
ybs[4479]=['',3.0553734,-0.433836,6.42];
ybs[4480]=['',3.0591055,-0.2548968,6.21];
ybs[4481]=['',3.0590464,-0.2924505,6.19];
ybs[4482]=['',3.0571269,-1.1437781,5.17];
ybs[4483]=['',3.0621101,1.0094046,6.37];
ybs[4484]=['ο Hya',3.0605977,-0.6087821,4.7];
ybs[4485]=['92 Leo',3.0633055,0.3703035,5.26];
ybs[4486]=['61 UMa',3.0645155,0.5945587,5.33];
ybs[4487]=['',3.0626359,-0.9443027,5.96];
ybs[4488]=['',3.0646734,-0.5119461,6.44];
ybs[4489]=['',3.0633393,-1.0860477,4.94];
ybs[4490]=['',3.0676013,0.9605682,6.27];
ybs[4491]=['62 UMa',3.0667768,0.5517008,5.73];
ybs[4492]=['',3.0654327,-0.7545374,5.55];
ybs[4493]=['',3.0672513,-0.5696007,5.22];
ybs[4494]=['3 Dra',3.0710028,1.1625459,5.3];
ybs[4495]=['',3.068983,0.3852784,6.59];
ybs[4496]=['',3.0687197,-0.356569,6.22];
ybs[4497]=['',3.0626505,-1.4527414,6.33];
ybs[4498]=['',3.0747474,-0.6514676,5.98];
ybs[4499]=['',3.0716449,-1.386532,6.39];
ybs[4500]=['',3.0768832,-0.1189147,6.07];
ybs[4501]=['',3.0748466,-1.0930215,5.03];
ybs[4502]=['',3.0782922,0.4377675,6.02];
ybs[4503]=['',3.0764267,-1.099809,6.1];
ybs[4504]=['ζ Crt',3.0805387,-0.3226582,4.73];
ybs[4505]=['ξ Vir',3.0828865,0.141759,4.85];
ybs[4506]=['',3.082368,-0.8588042,6.26];
ybs[4507]=['ν Vir',3.0853902,0.1115839,4.03];
ybs[4508]=['χ UMa',3.0863552,0.831532,3.71];
ybs[4509]=['',3.0846688,-0.7998173,5.29];
ybs[4510]=['λ Mus',3.0839346,-1.1670103,3.64];
ybs[4511]=['',3.0902166,0.9685205,5.27];
ybs[4512]=['',3.0879812,-1.0701401,4.11];
ybs[4513]=['',3.0881325,-0.7092448,4.91];
ybs[4514]=['',3.0907651,-0.6290716,6.17];
ybs[4515]=['',3.0914165,-0.5309842,6.48];
ybs[4516]=['',3.0915386,-1.0093693,5.41];
ybs[4517]=['93 Leo',3.0946897,0.3505086,4.53];
ybs[4518]=['4 Vir',3.094358,0.1415393,5.32];
ybs[4519]=['',3.0963998,-0.1823794,6.26];
ybs[4520]=['μ Mus',3.0954802,-1.1685147,4.72];
ybs[4521]=['',3.0975543,0.2469277,5.88];
ybs[4522]=['',3.0979357,-0.4692487,5.11];
ybs[4523]=['',3.0991617,-0.007939,6.15];
ybs[4524]=['β Leo',3.0993646,0.2519502,2.14];
ybs[4525]=['',3.1001892,0.2811117,6.04];
ybs[4526]=['',3.1021776,0.6072941,5.7];
ybs[4527]=['',3.1018586,-1.1156948,4.32];
ybs[4528]=['',3.102917,-1.2280506,4.97];
ybs[4529]=['',3.1048297,-0.2792558,6.13];
ybs[4530]=['β Vir',3.1064736,0.0284213,3.61];
ybs[4531]=['',3.1052365,-1.0958179,5.7];
ybs[4532]=['',3.1060894,-0.4784658,6.48];
ybs[4533]=['',3.107479,0.211928,6.35];
ybs[4534]=['',3.1079535,-0.0954632,5.64];
ybs[4535]=['',3.1085413,0.5801245,6.27];
ybs[4536]=['',3.1083482,-0.7908072,4.46];
ybs[4537]=['',3.1093764,-0.215096,6.35];
ybs[4538]=['',3.1107793,-0.5405515,5.85];
ybs[4539]=['',3.1113551,-1.1404406,4.9];
ybs[4540]=['',3.1164923,0.6559343,6.45];
ybs[4541]=['',3.112784,-0.9970038,5.57];
ybs[4542]=['β Hya',3.1160905,-0.5941868,4.28];
ybs[4543]=['',3.118434,-0.6144086,6.17];
ybs[4544]=['γ UMa',3.120228,0.9347698,2.44];
ybs[4545]=['',3.1201868,0.0072533,6.3];
ybs[4546]=['',3.1216435,-1.0043735,6.06];
ybs[4547]=['',3.1227281,-0.6612224,6.46];
ybs[4548]=['',3.123959,-0.4511721,5.3];
ybs[4549]=['6 Vir',3.1254871,0.1449935,5.58];
ybs[4550]=['65 UMa',3.1257162,0.8087955,6.54];
ybs[4551]=['65 UMa',3.1261151,0.8086695,7.03];
ybs[4552]=['',3.1263124,0.6391398,6.49];
ybs[4553]=['',3.1251498,-1.1068051,5.91];
ybs[4554]=['95 Leo',3.128211,0.2707056,5.53];
ybs[4555]=['',3.1281503,-0.4993967,5.93];
ybs[4556]=['66 UMa',3.1295553,0.9854518,5.84];
ybs[4557]=['η Crt',3.129673,-0.3017188,5.18];
ybs[4558]=['',3.129203,-0.6950869,6.13];
ybs[4559]=['',3.1335346,1.0718552,6.22];
ybs[4560]=['',3.1327813,-0.8239505,6.26];
ybs[4561]=['',3.134235,-0.5838418,6.21];
ybs[4562]=['',3.1350621,0.7017483,6.62];
ybs[4563]=['',3.1368637,-1.0923192,5.57];
ybs[4564]=['',3.1388768,0.5609051,6.42];
ybs[4565]=['',3.1398627,1.0703812,6.76];
ybs[4566]=['',3.1394329,-0.9853016,5.44];
ybs[4567]=['',3.139812,-0.7170444,6.79];
ybs[4568]=['',3.1418013,-1.1253109,5.61];
ybs[4569]=['',3.1422988,-0.454576,6.43];
ybs[4570]=['',3.1429557,0.0068793,6.17];
ybs[4571]=['',3.1439888,0.5765015,5.96];
ybs[4572]=['',3.1434985,-0.9046576,6.05];
ybs[4573]=['ε Cha',3.1454335,-1.3676111,4.91];
ybs[4574]=['',3.146864,0.5916422,6.5];
ybs[4575]=['7 Vir',3.1468459,0.061416,5.37];
ybs[4576]=['',3.1483647,1.4087714,6.17];
ybs[4577]=['',3.1503105,-0.1846995,5.55];
ybs[4578]=['',3.1501679,-0.3835119,6.28];
ybs[4579]=['π Vir',3.1508805,0.1130585,4.66];
ybs[4580]=['',3.1508006,-0.3454928,5.26];
ybs[4581]=['',3.1515665,-0.0332389,6.31];
ybs[4582]=['',3.1535785,-1.0060078,6.16];
ybs[4583]=['',3.1542903,0.6266701,5.59];
ybs[4584]=['67 UMa',3.1562663,0.7489063,5.21];
ybs[4585]=['',3.1576759,-1.4925863,6.05];
ybs[4586]=['',3.1579753,-1.2500969,6.42];
ybs[4587]=['',3.1586295,-1.2100124,5.89];
ybs[4588]=['',3.1595558,-0.1364846,6.22];
ybs[4589]=['θ1 Cru',3.1603466,-1.1073967,4.33];
ybs[4590]=['',3.1630823,-0.7429961,5.15];
ybs[4591]=['',3.1635451,-1.2976568,6.44];
ybs[4592]=['2 Com',3.1657122,0.3721531,5.87];
ybs[4593]=['θ2 Cru',3.1660234,-1.1048269,4.72];
ybs[4594]=['',3.1674781,-1.1949488,5.35];
ybs[4595]=['κ Cha',3.168143,-1.3378913,5.04];
ybs[4596]=['',3.1659031,1.4869777,6.27];
ybs[4597]=['',3.1687824,-1.0664926,5.96];
ybs[4598]=['ο Vir',3.1697881,0.1500408,4.12];
ybs[4599]=['',3.1697266,1.3398802,5.8];
ybs[4600]=['',3.1716482,1.0960093,6.13];
ybs[4601]=['',3.1729068,-1.1463944,6.33];
ybs[4602]=['',3.1730608,-0.6253555,6.23];
ybs[4603]=['',3.1732409,-0.0570375,6.37];
ybs[4604]=['',3.1748776,-1.2005626,6.23];
ybs[4605]=['',3.1750938,-1.1492208,6.06];
ybs[4606]=['η Cru',3.1772639,-1.1300996,4.15];
ybs[4607]=['',3.1815761,-1.3177802,5.18];
ybs[4608]=['',3.1824673,-0.8865869,4.47];
ybs[4609]=['',3.1824385,-0.8883662,6.37];
ybs[4610]=['',3.1831526,-0.8522281,5.34];
ybs[4611]=['δ Cen',3.1836564,-0.8876533,2.6];
ybs[4612]=['',3.1839363,-1.0643631,6.22];
ybs[4613]=['α Crv',3.1838251,-0.4339792,4.02];
ybs[4614]=['',3.1859877,-0.7760151,5.75];
ybs[4615]=['',3.1860287,-0.722002,5.48];
ybs[4616]=['10 Vir',3.1893416,0.0307443,5.95];
ybs[4617]=['',3.1893905,1.300709,6.35];
ybs[4618]=['',3.1909637,-0.6080945,6.17];
ybs[4619]=['11 Vir',3.190941,0.0989723,5.72];
ybs[4620]=['ε Crv',3.191297,-0.3971666,3];
ybs[4621]=['',3.1932536,-0.6633387,6.06];
ybs[4622]=['3 Com',3.1929681,0.2909975,6.39];
ybs[4623]=['',3.1939966,0.4737725,6.01];
ybs[4624]=['',3.1956475,-1.0718716,6.08];
ybs[4625]=['3 Crv',3.1954051,-0.4143188,5.46];
ybs[4626]=['',3.1954026,-0.7951545,6.61];
ybs[4627]=['',3.1975102,-0.8987686,6.23];
ybs[4628]=['ρ Cen',3.1980777,-0.9163819,3.96];
ybs[4629]=['',3.1942531,1.4237309,6];
ybs[4630]=['4 Com',3.198729,0.4491444,5.66];
ybs[4631]=['68 UMa',3.198136,0.9934108,6.43];
ybs[4632]=['',3.1994477,0.4956721,6.49];
ybs[4633]=['5 Com',3.2000583,0.3561477,5.57];
ybs[4634]=['',3.201297,-1.101076,5.92];
ybs[4635]=['',3.2032233,-1.2267589,6.17];
ybs[4636]=['',3.1996971,1.3522846,5.14];
ybs[4637]=['',3.204838,-0.5979795,6.5];
ybs[4638]=['',3.2057529,-0.6818182,5.76];
ybs[4639]=['',3.2086022,-1.3737439,6.35];
ybs[4640]=['12 Vir',3.2056522,0.1767335,5.85];
ybs[4641]=['',3.2065631,-0.5921712,6.33];
ybs[4642]=['',3.2085041,-0.8004081,5.31];
ybs[4643]=['',3.2097048,-1.1265178,6.22];
ybs[4644]=['',3.2111004,0.9302365,6.16];
ybs[4645]=['',3.2125544,-0.3661743,5.83];
ybs[4646]=['δ Cru',3.2134253,-1.0277364,2.8];
ybs[4647]=['',3.2133225,-0.182362,6.11];
ybs[4648]=['',3.2148904,-0.7338955,6.26];
ybs[4649]=['',3.2126882,1.2228461,5.71];
ybs[4650]=['δ UMa',3.2141275,0.9930301,3.31];
ybs[4651]=['',3.2160132,-0.4099719,6.54];
ybs[4652]=['γ Crv',3.2160954,-0.3085392,2.59];
ybs[4653]=['6 Com',3.2168537,0.2576604,5.1];
ybs[4654]=['',3.2191625,-1.26974,6.22];
ybs[4655]=['',3.2149977,1.2638763,6.29];
ybs[4656]=['2 CVn',3.2172887,0.7072815,5.66];
ybs[4657]=['7 Com',3.2182977,0.4155498,4.95];
ybs[4658]=['',3.2189585,0.5746562,5];
ybs[4659]=['',3.2220942,-1.1489287,6.06];
ybs[4660]=['',3.2215434,-0.293732,6.05];
ybs[4661]=['ε Mus',3.2241974,-1.1885133,4.11];
ybs[4662]=['',3.2231432,0.9259868,5.81];
ybs[4663]=['',3.2233625,0.5026767,5.7];
ybs[4664]=['β Cha',3.2281734,-1.3866317,4.26];
ybs[4665]=['',3.2248277,-0.6323301,6.15];
ybs[4666]=['',3.2244184,0.2619427,6.34];
ybs[4667]=['',3.2262869,-0.0713907,6.99];
ybs[4668]=['',3.2263233,-0.0712889,6.54];
ybs[4669]=['ζ Cru',3.2279173,-1.1194364,4.04];
ybs[4670]=['',3.227794,0.5255753,6.23];
ybs[4671]=['13 Vir',3.2285452,-0.0161117,5.9];
ybs[4672]=['',3.2302466,-0.9647997,5];
ybs[4673]=['',3.2178167,1.5033079,6.33];
ybs[4674]=['',3.2300284,0.4515496,6.48];
ybs[4675]=['8 Com',3.2312762,0.3996601,6.27];
ybs[4676]=['',3.2101025,1.5270475,6.28];
ybs[4677]=['',3.2284853,1.3094271,5.38];
ybs[4678]=['9 Com',3.2320161,0.48906,6.33];
ybs[4679]=['η Vir',3.2339337,-0.0140114,3.89];
ybs[4680]=['3 CVn',3.2332736,0.8525639,5.29];
ybs[4681]=['',3.2352045,-0.3894073,5.97];
ybs[4682]=['',3.236865,-1.1515437,6.21];
ybs[4683]=['',3.2356622,0.4622264,5.54];
ybs[4684]=['',3.2355198,0.4514488,6.15];
ybs[4685]=['16 Vir',3.235853,0.0554434,4.96];
ybs[4686]=['ζ Crv',3.2368789,-0.3901098,5.21];
ybs[4687]=['11 Com',3.2373945,0.3081723,4.74];
ybs[4688]=['',3.237229,0.4698237,7.13];
ybs[4689]=['',3.238443,-0.2391337,5.14];
ybs[4690]=['ε Cru',3.2406756,-1.0565679,3.59];
ybs[4691]=['70 UMa',3.2376795,1.0075452,5.55];
ybs[4692]=['',3.2432368,-0.9862935,5.92];
ybs[4693]=['ζ2 Mus',3.2441694,-1.180849,5.15];
ybs[4694]=['ζ1 Mus',3.2445291,-1.1945595,5.74];
ybs[4695]=['',3.2437415,0.4300171,6.19];
ybs[4696]=['',3.247036,-1.0090061,5.39];
ybs[4697]=['12 Com',3.2451524,0.4487312,4.81];
ybs[4698]=['17 Vir',3.2453699,0.090231,6.4];
ybs[4699]=['',3.2630733,-1.5025776,6.33];
ybs[4700]=['',3.2490365,-1.1827629,6.36];
ybs[4701]=['6 Crv',3.2491145,-0.435917,5.68];
ybs[4702]=['',3.2501816,-0.6204369,5.32];
ybs[4703]=['',3.250318,-0.6883349,6.4];
ybs[4704]=['',3.250898,-0.6814941,5.79];
ybs[4705]=['4 CVn',3.2506146,0.7401444,6.06];
ybs[4706]=['5 CVn',3.2515798,0.8975636,4.8];
ybs[4707]=['13 Com',3.2530106,0.4531402,5.18];
ybs[4708]=['',3.2552747,-0.724656,6.25];
ybs[4709]=['',3.2536091,0.4441373,6.42];
ybs[4710]=['',3.2580132,-1.150278,6.3];
ybs[4711]=['',3.2570187,-0.7443825,6.11];
ybs[4712]=['',3.257061,-0.2050031,5.95];
ybs[4713]=['',3.2576348,-0.4866796,6.09];
ybs[4714]=['',3.2579307,-0.6164836,5.73];
ybs[4715]=['',3.2571343,0.415224,6.03];
ybs[4716]=['71 UMa',3.2559828,0.9885886,5.81];
ybs[4717]=['',3.256081,1.1112028,6.32];
ybs[4718]=['6 CVn',3.2596347,0.6786385,5.02];
ybs[4719]=['',3.2633231,-1.1040593,4.86];
ybs[4720]=['α1 Cru',3.2636879,-1.1036519,1.33];
ybs[4721]=['α2 Cru',3.2637317,-1.1036568,1.73];
ybs[4722]=['',3.2631707,-0.9003502,4.82];
ybs[4723]=['14 Com',3.2621269,0.4735581,4.95];
ybs[4724]=['',3.2643506,-0.8560622,6.26];
ybs[4725]=['',3.2644622,-0.5753551,5.55];
ybs[4726]=['',3.2672636,-1.1156937,6];
ybs[4727]=['γ Com',3.2644593,0.4910121,4.36];
ybs[4728]=['16 Com',3.2646863,0.465831,5];
ybs[4729]=['',3.2674553,-1.0319663,5.5];
ybs[4730]=['',3.2614259,1.2530463,6.24];
ybs[4731]=['',3.2678989,0.1479153,6.37];
ybs[4732]=['',3.2685623,-0.2926444,6.35];
ybs[4733]=['σ Cen',3.2697775,-0.8790505,3.91];
ybs[4734]=['',3.2712467,-1.1253306,6.04];
ybs[4735]=['73 UMa',3.267007,0.9700088,5.7];
ybs[4736]=['',3.268658,-0.082914,6.22];
ybs[4737]=['',3.2716632,-1.0808925,6.22];
ybs[4738]=['',3.2711254,-0.6837623,5.44];
ybs[4739]=['',3.2721308,-0.9868626,6.15];
ybs[4740]=['',3.2718666,0.4553805,6.54];
ybs[4741]=['',3.2723409,0.4496646,6.65];
ybs[4742]=['17 Com',3.2730745,0.4499024,5.29];
ybs[4743]=['18 Com',3.2754253,0.4184194,5.48];
ybs[4744]=['',3.2780079,-0.9889019,5.8];
ybs[4745]=['',3.2781022,-0.730792,6.02];
ybs[4746]=['20 Com',3.2766301,0.3623461,5.69];
ybs[4747]=['δ Crv',3.2774764,-0.2906103,2.95];
ybs[4748]=['',3.2783966,-0.2361121,6.35];
ybs[4749]=['',3.2793851,-0.4159438,5.63];
ybs[4750]=['74 UMa',3.2772624,1.0170146,5.35];
ybs[4751]=['7 CVn',3.2777839,0.8971057,6.21];
ybs[4752]=['75 UMa',3.2777622,1.0233271,6.08];
ybs[4753]=['γ Cru',3.283566,-0.9991733,1.63];
ybs[4754]=['γ Cru',3.2840617,-0.9986108,6.42];
ybs[4755]=['4 Dra',3.2776261,1.2054279,4.95];
ybs[4756]=['21 Com',3.2822291,0.4264209,5.46];
ybs[4757]=['',3.2811829,0.9240043,6.21];
ybs[4758]=['',3.2858341,-1.0394994,5.48];
ybs[4759]=['',3.2885448,-1.2764755,5.88];
ybs[4760]=['',3.2838471,0.1303603,6.05];
ybs[4761]=['',3.2870724,-1.1107472,5.95];
ybs[4762]=['',3.285179,-0.0905398,6.19];
ybs[4763]=['γ Mus',3.2898052,-1.2613149,3.87];
ybs[4764]=['',3.2872419,-0.5701749,6.46];
ybs[4765]=['η Crv',3.2871022,-0.2850319,4.31];
ybs[4766]=['',3.2894013,-0.2442437,5.74];
ybs[4767]=['20 Vir',3.2912094,0.1773364,6.26];
ybs[4768]=['',3.2928175,-0.347789,6.26];
ybs[4769]=['',3.2936368,-0.2262847,5.58];
ybs[4770]=['22 Com',3.2933893,0.4214651,6.29];
ybs[4771]=['21 Vir',3.2945245,-0.1673214,5.48];
ybs[4772]=['',3.2957936,-0.8734375,6.38];
ybs[4773]=['',3.2936585,0.5779243,5.42];
ybs[4774]=['',3.2942744,0.5803195,6.24];
ybs[4775]=['β CVn',3.2939854,0.7194706,4.26];
ybs[4776]=['β Crv',3.2972692,-0.4107017,2.65];
ybs[4777]=['κ Dra',3.292218,1.2156818,3.87];
ybs[4778]=['',3.2988635,-0.7820491,5.77];
ybs[4779]=['23 Com',3.2989883,0.3926014,4.81];
ybs[4780]=['',3.3025983,-1.0816965,6.22];
ybs[4781]=['24 Com',3.3001239,0.3183913,6.56];
ybs[4782]=['24 Com',3.3002329,0.3183865,5.02];
ybs[4783]=['',3.300229,0.3795506,5.85];
ybs[4784]=['',3.3034254,-0.7183186,5.13];
ybs[4785]=['6 Dra',3.297625,1.219761,4.94];
ybs[4786]=['',3.3045563,-0.6982129,5.8];
ybs[4787]=['',3.3041919,-0.3606179,6.2];
ybs[4788]=['α Mus',3.3103472,-1.2089912,2.69];
ybs[4789]=['25 Vir',3.3076383,-0.1041355,5.87];
ybs[4790]=['',3.3051992,1.0358933,5.5];
ybs[4791]=['25 Com',3.308272,0.2959185,5.68];
ybs[4792]=['τ Cen',3.3120338,-0.8495495,3.86];
ybs[4793]=['',3.3117874,-0.4760103,5.45];
ybs[4794]=['',3.319889,-1.3177893,6.49];
ybs[4795]=['',3.3131755,0.0549437,6.33];
ybs[4796]=['',3.3176727,-1.1750852,6.25];
ybs[4797]=['',3.3144943,0.0300248,5.71];
ybs[4798]=['',3.3150105,0.1196234,7.08];
ybs[4799]=['',3.3162561,-0.3208729,6];
ybs[4800]=['',3.3177375,-0.5333128,5.89];
ybs[4801]=['9 CVn',3.3158938,0.711048,6.37];
ybs[4802]=['',3.3172235,0.3931368,6.38];
ybs[4803]=['χ Vir',3.3183759,-0.1418934,4.66];
ybs[4804]=['',3.3222547,-1.1631908,6.26];
ybs[4805]=['26 Com',3.3176145,0.3652651,5.46];
ybs[4806]=['',3.3181722,0.6251352,6.45];
ybs[4807]=['',3.3214114,-0.7002569,4.64];
ybs[4808]=['',3.3280938,-0.8077325,5.84];
ybs[4809]=['γ Cen',3.3287242,-0.8568487,2.17];
ybs[4810]=['',3.331884,-1.2137285,6.33];
ybs[4811]=['',3.3272321,-0.2294713,6.08];
ybs[4812]=['',3.3272466,-0.2294955,5.98];
ybs[4813]=['',3.3308395,-1.0440538,4.93];
ybs[4814]=['27 Vir',3.3283899,0.1796344,6.19];
ybs[4815]=['γ Vir',3.3288593,-0.0276378,3.65];
ybs[4816]=['γ Vir',3.3288593,-0.0276378,3.68];
ybs[4817]=['',3.329704,-0.3471926,6.03];
ybs[4818]=['ρ Vir',3.3297507,0.1763043,4.88];
ybs[4819]=['31 Vir',3.3300684,0.1164591,5.59];
ybs[4820]=['',3.3348768,-1.1029181,5.31];
ybs[4821]=['',3.3334175,-0.8542868,4.66];
ybs[4822]=['',3.3346049,-0.978801,6.08];
ybs[4823]=['76 UMa',3.3275927,1.0922088,6.07];
ybs[4824]=['',3.3360338,-0.9827952,6];
ybs[4825]=['',3.3375035,-1.0303888,6.4];
ybs[4826]=['',3.3369926,-0.7035711,6.44];
ybs[4827]=['',3.3374784,-0.0298592,5.93];
ybs[4828]=['',3.3393033,-0.6367482,6.39];
ybs[4829]=['',3.3393434,-0.4966806,5.48];
ybs[4830]=['',3.3341849,1.0650282,6.38];
ybs[4831]=['',3.3447693,-1.2036627,6.16];
ybs[4832]=['ι Cru',3.3470419,-1.0666532,4.69];
ybs[4833]=['',3.3405793,0.7674087,6.33];
ybs[4834]=['β Mus',3.3502253,-1.1910404,3.05];
ybs[4835]=['10 CVn',3.3430011,0.6832122,5.95];
ybs[4836]=['',3.3435107,0.7907491,4.99];
ybs[4837]=['32 Vir',3.3460555,0.1315926,5.22];
ybs[4838]=['',3.3501707,-0.9882476,4.65];
ybs[4839]=['33 Vir',3.3493399,0.1641737,5.67];
ybs[4840]=['',3.3514521,-0.5837959,5.86];
ybs[4841]=['27 Com',3.3504545,0.2870018,5.12];
ybs[4842]=['',3.3381728,1.4047684,6.4];
ybs[4843]=['β Cru',3.3561461,-1.0440904,1.25];
ybs[4844]=['',3.3522645,0.1015323,6.34];
ybs[4845]=['34 Vir',3.353033,0.2063785,6.07];
ybs[4846]=['',3.3546323,-0.1123178,6.26];
ybs[4847]=['',3.3562856,-0.4360708,6.44];
ybs[4848]=['35 Vir',3.3558552,0.0600292,6.41];
ybs[4849]=['',3.3525661,1.0934033,5.89];
ybs[4850]=['',3.3586963,-0.4839934,5.66];
ybs[4851]=['28 Com',3.3574304,0.2342187,6.56];
ybs[4852]=['',3.3657481,-1.2587224,5.55];
ybs[4853]=['7 Dra',3.3534973,1.1633819,5.43];
ybs[4854]=['',3.3596942,0.4312191,6.31];
ybs[4855]=['29 Com',3.360325,0.2441589,5.7];
ybs[4856]=['11 CVn',3.358979,0.843582,6.27];
ybs[4857]=['',3.3584893,1.0504566,5.85];
ybs[4858]=['',3.3670367,-1.0565154,6.75];
ybs[4859]=['30 Com',3.3618725,0.4785526,5.78];
ybs[4860]=['ι Oct',3.3938431,-1.4826048,5.46];
ybs[4861]=['',3.3672579,-0.8481035,6.24];
ybs[4862]=['',3.370151,-0.923636,5.73];
ybs[4863]=['',3.3662829,0.3967183,6.43];
ybs[4864]=['',3.3685917,-0.5957233,4.91];
ybs[4865]=['',3.3656154,0.6524718,5.89];
ybs[4866]=['',3.3718387,-1.0552718,5.72];
ybs[4867]=['',3.3713649,-0.1827575,6.41];
ybs[4868]=['37 Vir',3.3722565,0.0510299,6.02];
ybs[4869]=['',3.3741823,-0.6948793,5.98];
ybs[4870]=['',3.3749569,-0.8417193,6.33];
ybs[4871]=['',3.3740855,-0.4689852,6.15];
ybs[4872]=['',3.3765033,-0.9418179,6.24];
ybs[4873]=['31 Com',3.3723682,0.4783545,4.94];
ybs[4874]=['32 Com',3.3746926,0.2956779,6.32];
ybs[4875]=['',3.3794025,-0.9614174,5.93];
ybs[4876]=['',3.3758144,0.2790737,6.3];
ybs[4877]=['',3.380893,-1.0552475,5.76];
ybs[4878]=['',3.379467,-0.8565376,4.33];
ybs[4879]=['',3.3807045,-0.7035684,4.27];
ybs[4880]=['κ Cru',3.3828884,-1.0560899,5.9];
ybs[4881]=['38 Vir',3.3791774,-0.0643277,6.11];
ybs[4882]=['',3.3570004,1.4535938,5.85];
ybs[4883]=['',3.3575028,1.453502,5.28];
ybs[4884]=['35 Com',3.3794073,0.3684802,4.9];
ybs[4885]=['',3.3852159,-1.0221177,6.58];
ybs[4886]=['',3.3811404,-0.076035,6.44];
ybs[4887]=['λ Cru',3.3864981,-1.0346155,4.62];
ybs[4888]=['μ1 Cru',3.3861646,-1.0002521,4.03];
ybs[4889]=['μ2 Cru',3.3862519,-1.0000871,5.17];
ybs[4890]=['41 Vir',3.3818195,0.214432,6.25];
ybs[4891]=['',3.3841683,-0.2056191,6];
ybs[4892]=['ψ Vir',3.3843286,-0.1687974,4.79];
ybs[4893]=['',3.3875025,-0.7729075,5.89];
ybs[4894]=['',3.3832711,0.5829737,6.26];
ybs[4895]=['ε UMa',3.3819857,0.974368,1.77];
ybs[4896]=['',3.3890021,-0.7513324,5.47];
ybs[4897]=['',3.3955508,-1.2621772,5.93];
ybs[4898]=['',3.3920937,-0.9942854,5.32];
ybs[4899]=['',3.386185,0.8214261,5.84];
ybs[4900]=['δ Vir',3.3896534,0.0569883,3.38];
ybs[4901]=['',3.3910906,-0.2698141,6.17];
ybs[4902]=['',3.3939001,-0.4641259,6.62];
ybs[4903]=['',3.3968414,-0.8958896,5.16];
ybs[4904]=['α1 CVn',3.3909882,0.6664097,5.6];
ybs[4905]=['α2 CVn',3.3910826,0.6664728,2.9];
ybs[4906]=['8 Dra',3.3878716,1.1398092,5.24];
ybs[4907]=['',3.3918926,0.9419056,5.82];
ybs[4908]=['',3.3984471,-0.3994345,6.31];
ybs[4909]=['',3.3957247,0.8036341,6.12];
ybs[4910]=['36 Com',3.4039823,0.3015514,4.78];
ybs[4911]=['44 Vir',3.4074222,-0.0688297,5.79];
ybs[4912]=['',3.4116589,-0.5870736,6.02];
ybs[4913]=['δ Mus',3.4207236,-1.2510546,3.62];
ybs[4914]=['37 Com',3.4096998,0.5350025,4.9];
ybs[4915]=['46 Vir',3.4115195,-0.0610895,5.99];
ybs[4916]=['',3.4114847,0.3183742,6.2];
ybs[4917]=['',3.4012504,1.3149419,6.01];
ybs[4918]=['9 Dra',3.4070867,1.1600424,5.32];
ybs[4919]=['38 Com',3.4137381,0.2965591,5.96];
ybs[4920]=['',3.4243033,-1.249782,6.03];
ybs[4921]=['78 UMa',3.4111137,0.981483,4.93];
ybs[4922]=['ε Vir',3.4182446,0.1889817,2.83];
ybs[4923]=['ξ1 Cen',3.4251494,-0.8667007,4.85];
ybs[4924]=['',3.4153541,1.1079155,6];
ybs[4925]=['',3.4255624,-0.3615292,5.58];
ybs[4926]=['',3.4194291,1.0399521,6.53];
ybs[4927]=['48 Vir',3.4259602,-0.066224,6.59];
ybs[4928]=['',3.4304111,-0.7213015,6.26];
ybs[4929]=['',3.4338017,-0.9118601,6.43];
ybs[4930]=['',3.4370416,-0.8481292,4.71];
ybs[4931]=['',3.4382163,-0.7281369,5.59];
ybs[4932]=['ξ2 Cen',3.4398459,-0.8733037,4.27];
ybs[4933]=['14 CVn',3.4334351,0.6225269,5.25];
ybs[4934]=['',3.4423791,-1.0470398,5.99];
ybs[4935]=['',3.4338059,0.787805,5.63];
ybs[4936]=['39 Com',3.4363283,0.3669156,5.99];
ybs[4937]=['',3.4395012,-0.6281868,6.54];
ybs[4938]=['',3.4354143,0.5043791,6.54];
ybs[4939]=['40 Com',3.4364098,0.3924459,5.6];
ybs[4940]=['',3.4277585,1.2722466,6.31];
ybs[4941]=['',3.4431461,-0.9353236,5.71];
ybs[4942]=['θ Mus',3.4458249,-1.1420852,5.51];
ybs[4943]=['',3.4354154,1.0805562,6.14];
ybs[4944]=['41 Com',3.4398304,0.4798651,4.8];
ybs[4945]=['49 Vir',3.4434606,-0.1897281,5.19];
ybs[4946]=['',3.4429476,0.4786649,6.19];
ybs[4947]=['',3.4462539,-0.159081,5.55];
ybs[4948]=['ψ Hya',3.4486863,-0.4057575,4.95];
ybs[4949]=['',3.4493091,-0.1687461,6.32];
ybs[4950]=['',3.4489174,0.1726498,5.78];
ybs[4951]=['50 Vir',3.4515674,-0.182552,5.94];
ybs[4952]=['',3.4514,0.2917946,5.91];
ybs[4953]=['θ Vir',3.4523569,-0.0989404,4.38];
ybs[4954]=['',3.4504073,0.6508858,6.02];
ybs[4955]=['',3.4577108,-0.9197313,6.06];
ybs[4956]=['',3.4626616,-1.2229789,5.91];
ybs[4957]=['15 CVn',3.4506245,0.6702736,6.28];
ybs[4958]=['α Com',3.4522264,0.303678,5.22];
ybs[4959]=['α Com',3.4522264,0.303678,5.22];
ybs[4960]=['',3.4581755,-0.7393704,5.79];
ybs[4961]=['17 CVn',3.4521641,0.6696639,5.91];
ybs[4962]=['',3.4622095,-1.1071038,6.33];
ybs[4963]=['',3.4592521,-0.7591937,5.25];
ybs[4964]=['',3.4504302,1.0838344,6.54];
ybs[4965]=['',3.463788,-1.0480764,4.6];
ybs[4966]=['',3.4749904,-1.3714154,5.85];
ybs[4967]=['',3.466499,-1.1581369,5.9];
ybs[4968]=['',3.4600886,-0.465677,6.5];
ybs[4969]=['',3.4620414,-0.6620495,4.85];
ybs[4970]=['',3.4666006,-1.0462563,6.16];
ybs[4971]=['53 Vir',3.4617021,-0.2849809,5.04];
ybs[4972]=['',3.4656271,-0.7475097,6.22];
ybs[4973]=['β Com',3.4602878,0.4843013,4.26];
ybs[4974]=['',3.4615081,0.4211212,6.33];
ybs[4975]=['',3.4682298,-0.887139,5.89];
ybs[4976]=['',3.4634717,0.1994319,5.77];
ybs[4977]=['',3.4635888,0.3250181,6.53];
ybs[4978]=['',3.4720981,-1.0264814,5.89];
ybs[4979]=['',3.4723163,-1.033802,4.92];
ybs[4980]=['54 Vir',3.4678021,-0.3308445,6.28];
ybs[4981]=['',3.4704914,-0.7551709,6.16];
ybs[4982]=['',3.4662416,0.3245886,6.11];
ybs[4983]=['η Mus',3.4773211,-1.1872322,4.8];
ybs[4984]=['',3.4782979,-1.2183905,6.37];
ybs[4985]=['55 Vir',3.4710281,-0.3501135,5.33];
ybs[4986]=['',3.4739779,-0.8567077,5.89];
ybs[4987]=['',3.4680833,0.6985415,4.92];
ybs[4988]=['',3.4720811,0.1955211,5.67];
ybs[4989]=['',3.4755974,-0.637047,6.19];
ybs[4990]=['',3.4836531,-1.1391239,6.07];
ybs[4991]=['57 Vir',3.4788814,-0.3503207,5.22];
ybs[4992]=['',3.4858526,-1.1678377,4.87];
ybs[4993]=['',3.4654947,1.2683224,6.59];
ybs[4994]=['19 CVn',3.4759785,0.7108085,5.79];
ybs[4995]=['',3.4805586,-0.0265171,6.68];
ybs[4996]=['',3.4830315,-0.5521309,5.1];
ybs[4997]=['',3.4794546,0.3302663,6.45];
ybs[4998]=['',3.4848265,-0.7698303,5.84];
ybs[4999]=['',3.4586937,1.4022284,6.25];
ybs[5000]=['',3.4807513,0.3430712,6.45];
ybs[5001]=['59 Vir',3.4819358,0.1622366,5.22];
ybs[5002]=['',3.4956764,-1.2594938,6.04];
ybs[5003]=['',3.4839899,0.236439,5.33];
ybs[5004]=['',3.4852322,-0.0140537,6.37];
ybs[5005]=['σ Vir',3.4856119,0.0932214,4.8];
ybs[5006]=['',3.4909369,-0.8973509,6.19];
ybs[5007]=['20 CVn',3.4847314,0.70588,4.73];
ybs[5008]=['',3.4787905,1.1916979,6.2];
ybs[5009]=['61 Vir',3.4894455,-0.3218343,4.74];
ybs[5010]=['γ Hya',3.4917803,-0.4066604,3];
ybs[5011]=['',3.4910758,0.062125,6.62];
ybs[5012]=['',3.488889,0.5928829,5.82];
ybs[5013]=['21 CVn',3.4875223,0.8648721,5.15];
ybs[5014]=['',3.5000996,-1.0454736,6.18];
ybs[5015]=['',3.4915107,0.6108619,6.02];
ybs[5016]=['',3.4999773,-0.9228593,5.48];
ybs[5017]=['',3.5008736,-0.9761348,6.02];
ybs[5018]=['ι Cen',3.499359,-0.6429816,2.75];
ybs[5019]=['',3.5012227,-0.8204509,5.77];
ybs[5020]=['',3.5113424,-1.2614201,6.05];
ybs[5021]=['',3.4991212,0.0491095,6.26];
ybs[5022]=['23 CVn',3.4968268,0.6985255,5.6];
ybs[5023]=['',3.5029719,-0.3423745,6.21];
ybs[5024]=['',3.5089911,-1.0663907,6.18];
ybs[5025]=['',3.5091528,-1.0666718,4.53];
ybs[5026]=['',3.5071317,-0.9129922,5.83];
ybs[5027]=['',3.5034966,0.0342002,5.69];
ybs[5028]=['',3.5096402,-0.8389881,6.16];
ybs[5029]=['',3.5103828,-0.8498037,6.38];
ybs[5030]=['64 Vir',3.5054947,0.0877399,5.87];
ybs[5031]=['',3.5154411,-1.1285821,4.53];
ybs[5032]=['ι1 Mus',3.5217198,-1.3092528,5.05];
ybs[5033]=['',3.5104186,-0.5814979,6.22];
ybs[5034]=['63 Vir',3.5095804,-0.3117627,5.37];
ybs[5035]=['',3.5043344,0.7640253,6.35];
ybs[5036]=['',3.5140506,-0.8717966,6.48];
ybs[5037]=['65 Vir',3.5106764,-0.0881704,5.89];
ybs[5038]=['',3.5207997,-1.1276951,5.31];
ybs[5039]=['',3.5241084,-1.2348945,5.67];
ybs[5040]=['66 Vir',3.5160843,-0.0923449,5.75];
ybs[5041]=['ι2 Mus',3.5312947,-1.3058264,6.63];
ybs[5042]=['',3.512509,0.6441425,6.07];
ybs[5043]=['',3.5156194,0.2147601,6.44];
ybs[5044]=['ζ UMa',3.5120457,0.956406,2.27];
ybs[5045]=['ζ UMa',3.5121111,0.956343,3.95];
ybs[5046]=['α Vir',3.5189706,-0.1970185,0.98];
ybs[5047]=['',3.5180581,0.4141225,5.78];
ybs[5048]=['',3.5236189,-0.6960722,5.09];
ybs[5049]=['',3.5231658,-0.0230249,5.97];
ybs[5050]=['',3.5271849,-0.7264863,5.69];
ybs[5051]=['',3.5281707,-0.8599305,6.31];
ybs[5052]=['80 UMa',3.5176974,0.9575066,4.01];
ybs[5053]=['',3.5292366,-0.864065,6.28];
ybs[5054]=['68 Vir',3.5256598,-0.2240023,5.25];
ybs[5055]=['',3.5285038,-0.703185,6.4];
ybs[5056]=['',3.5368718,-1.2174394,6.2];
ybs[5057]=['',3.5226143,0.8011292,5.88];
ybs[5058]=['69 Vir',3.5289155,-0.2809989,4.76];
ybs[5059]=['',3.5378903,-1.1310057,6.11];
ybs[5060]=['',3.520575,1.1019014,6.5];
ybs[5061]=['',3.5383764,-0.8952013,5.06];
ybs[5062]=['70 Vir',3.5327003,0.2382838,4.98];
ybs[5063]=['',3.5201545,1.2612549,5.79];
ybs[5064]=['',3.5251415,1.1276394,6.66];
ybs[5065]=['',3.5255844,1.1273586,7.04];
ybs[5066]=['',3.5298317,0.9183832,6.34];
ybs[5067]=['',3.5321752,0.7086644,6.47];
ybs[5068]=['',3.5365137,-0.0260139,6.43];
ybs[5069]=['',3.530804,0.8807091,6.8];
ybs[5070]=['',3.5389147,-0.4085347,4.97];
ybs[5071]=['71 Vir',3.5361801,0.1866154,5.65];
ybs[5072]=['',3.5582818,-1.3560044,6.48];
ybs[5073]=['',3.5332649,0.8829996,6.43];
ybs[5074]=['κ Oct',3.6017097,-1.4953402,5.58];
ybs[5075]=['',3.5314626,1.0440485,5.4];
ybs[5076]=['',3.5396636,0.1230983,6.17];
ybs[5077]=['',3.5395007,0.1027554,6.51];
ybs[5078]=['72 Vir',3.5417452,-0.1151228,6.09];
ybs[5079]=['',3.5450945,-0.6899829,3.88];
ybs[5080]=['',3.5470552,-0.4928509,6.47];
ybs[5081]=['',3.5220701,1.370384,5.77];
ybs[5082]=['',3.5496278,-0.67238,6.16];
ybs[5083]=['',3.5575686,-1.1476844,6.37];
ybs[5084]=['73 Vir',3.5490198,-0.3290692,6.01];
ybs[5085]=['74 Vir',3.5484484,-0.1113738,4.69];
ybs[5086]=['',3.5444461,0.7326981,6.08];
ybs[5087]=['',3.5516253,-0.5029695,5.69];
ybs[5088]=['',3.5515422,-0.5181976,6.45];
ybs[5089]=['75 Vir',3.5525178,-0.270321,5.55];
ybs[5090]=['76 Vir',3.5528931,-0.1795974,5.21];
ybs[5091]=['',3.5530321,-0.127761,6.68];
ybs[5092]=['',3.5515697,0.4227389,6.11];
ybs[5093]=['',3.5604056,-0.8446872,6.33];
ybs[5094]=['',3.561029,-0.5835609,6.44];
ybs[5095]=['78 Vir',3.5577384,0.0616797,4.94];
ybs[5096]=['',3.5603963,-0.2328132,5.91];
ybs[5097]=['ζ Vir',3.5602604,-0.0125769,3.37];
ybs[5098]=['',3.5580413,0.6748193,6.37];
ybs[5099]=['81 UMa',3.5563766,0.9638348,5.6];
ybs[5100]=['',3.559972,0.6467795,4.98];
ybs[5101]=['80 Vir',3.5639563,-0.096354,5.73];
ybs[5102]=['24 CVn',3.5581087,0.8533133,4.7];
ybs[5103]=['',3.5730396,-1.0788937,5.63];
ybs[5104]=['',3.5638355,0.1759319,6.49];
ybs[5105]=['',3.5838766,-1.3230894,6.34];
ybs[5106]=['',3.5616755,0.7692064,6.84];
ybs[5107]=['',3.5702994,-0.6037444,6.5];
ybs[5108]=['',3.571709,-0.7726134,5.98];
ybs[5109]=['',3.5808353,-1.231656,6.1];
ybs[5110]=['',3.5699615,-0.4645934,5.78];
ybs[5111]=['',3.5730789,-0.8124929,5.9];
ybs[5112]=['',3.576854,-1.0216963,6.42];
ybs[5113]=['',3.5697995,0.4274156,5.74];
ybs[5114]=['',3.579833,-1.0078712,6.01];
ybs[5115]=['',3.5863712,-1.2376428,6.59];
ybs[5116]=['',3.5677015,0.8615355,6.49];
ybs[5117]=['25 CVn',3.5715977,0.6313011,4.82];
ybs[5118]=['',3.5783004,-0.5180941,5.83];
ybs[5119]=['',3.5749995,0.2474481,6.52];
ybs[5120]=['',3.58634,-1.1292332,5.79];
ybs[5121]=['',3.5563524,1.3338119,6.57];
ybs[5122]=['ε Cen',3.5842791,-0.9353191,2.3];
ybs[5123]=['',3.5722345,0.8829736,6.48];
ybs[5124]=['',3.5845902,-0.873951,6];
ybs[5125]=['',3.5828494,-0.6958902,6.27];
ybs[5126]=['',3.5834268,-0.7011935,5.6];
ybs[5127]=['',3.5788782,0.31663,6.48];
ybs[5128]=['',3.5813667,0.1853983,5.57];
ybs[5129]=['',3.5682317,1.2412428,5.5];
ybs[5130]=['',3.5939192,-1.0281755,5.38];
ybs[5131]=['',3.5924684,-0.9543979,5.01];
ybs[5132]=['82 UMa',3.5799233,0.9214949,5.46];
ybs[5133]=['',3.5839305,0.5391066,6.21];
ybs[5134]=['1 Boo',3.5859809,0.3461381,5.75];
ybs[5135]=['',3.5857063,0.4876793,6.23];
ybs[5136]=['',3.5904611,-0.4114227,6.59];
ybs[5137]=['',3.5917761,-0.5885239,6.05];
ybs[5138]=['',3.5838625,0.8795769,6.32];
ybs[5139]=['2 Boo',3.587513,0.390476,5.62];
ybs[5140]=['82 Vir',3.5906019,-0.1540445,5.01];
ybs[5141]=['',3.5978032,-0.9929303,6];
ybs[5142]=['',3.5973795,-0.8886035,6.41];
ybs[5143]=['',3.5833833,0.9963052,6.29];
ybs[5144]=['83 UMa',3.585193,0.9522229,4.66];
ybs[5145]=['',3.5970495,-0.7247269,5.98];
ybs[5146]=['',3.5929031,0.144259,6.16];
ybs[5147]=['',3.6003501,-0.7363542,5.98];
ybs[5148]=['',3.6033177,-0.8924807,6.47];
ybs[5149]=['84 Vir',3.5966998,0.0596096,5.36];
ybs[5150]=['',3.5933064,0.7252072,6.3];
ybs[5151]=['',3.5945628,0.6085283,5.98];
ybs[5152]=['',3.5877701,1.1292169,5.85];
ybs[5153]=['',3.6005448,-0.0981109,6.51];
ybs[5154]=['',3.5993445,0.3940565,6.13];
ybs[5155]=['83 Vir',3.6033316,-0.284514,5.6];
ybs[5156]=['',3.6046863,-0.4472065,6.21];
ybs[5157]=['',3.6084276,-0.4579412,5.81];
ybs[5158]=['1 Cen',3.6089117,-0.5788533,4.23];
ybs[5159]=['',3.5991313,0.9065581,6.02];
ybs[5160]=['85 Vir',3.6080749,-0.2773241,6.19];
ybs[5161]=['',3.6168158,-1.0945225,6.51];
ybs[5162]=['',3.6137938,-0.8997951,4.65];
ybs[5163]=['86 Vir',3.6095513,-0.2190139,5.51];
ybs[5164]=['',3.6144774,-0.6348385,5.15];
ybs[5165]=['',3.6172593,-0.8791383,5.91];
ybs[5166]=['',3.6180501,-0.8803882,5.45];
ybs[5167]=['',3.604621,0.9731483,6.5];
ybs[5168]=['',3.6151149,-0.1715786,6.05];
ybs[5169]=['',3.6096354,0.7150045,5.87];
ybs[5170]=['',3.6101162,0.6698931,5.94];
ybs[5171]=['87 Vir',3.6161442,-0.3138365,5.43];
ybs[5172]=['3 Boo',3.6122172,0.446464,5.95];
ybs[5173]=['',3.6136177,0.1087149,6.33];
ybs[5174]=['',3.5901337,1.3603362,5.91];
ybs[5175]=['τ Boo',3.6147509,0.3025545,4.5];
ybs[5176]=['',3.6130771,0.6705751,5.5];
ybs[5177]=['84 UMa',3.6106844,0.9479056,5.7];
ybs[5178]=['',3.6608034,-1.4448712,5.95];
ybs[5179]=['',3.6231069,-0.6252637,6.53];
ybs[5180]=['ν Cen',3.6258605,-0.7296994,3.41];
ybs[5181]=['η UMa',3.6150488,0.858559,1.86];
ybs[5182]=['2 Cen',3.6253724,-0.6033913,4.19];
ybs[5183]=['μ Cen',3.6263761,-0.7434191,3.04];
ybs[5184]=['',3.6378492,-1.2133809,5.75];
ybs[5185]=['',3.6204655,0.5422576,5.62];
ybs[5186]=['89 Vir',3.6268402,-0.3186099,4.97];
ybs[5187]=['',3.6281269,-0.5096737,6.18];
ybs[5188]=['',3.6293773,-0.6985123,6.44];
ybs[5189]=['',3.6215687,0.6880376,7.4];
ybs[5190]=['υ Boo',3.6244438,0.2736119,4.07];
ybs[5191]=['6 Boo',3.6253589,0.3690195,4.91];
ybs[5192]=['',3.6299495,-0.3493776,6.53];
ybs[5193]=['',3.5858131,1.4421553,5.98];
ybs[5194]=['',3.625129,0.6372524,6.38];
ybs[5195]=['',3.6287279,0.0938379,6.01];
ybs[5196]=['',3.6360352,-0.8206441,5.77];
ybs[5197]=['',3.6375998,-0.9238351,5.25];
ybs[5198]=['',3.6348912,-0.6379819,6.35];
ybs[5199]=['',3.6333924,-0.4278021,6.45];
ybs[5200]=['3 Cen',3.635738,-0.577961,4.56];
ybs[5201]=['3 Cen',3.6357745,-0.5779658,6.06];
ybs[5202]=['',3.6365231,-0.5539618,6.12];
ybs[5203]=['',3.6239616,1.0710775,5.96];
ybs[5204]=['',3.6309561,0.6047908,6.65];
ybs[5205]=['',3.6313005,0.6029053,5.87];
ybs[5206]=['',3.6272461,1.0195986,6.46];
ybs[5207]=['',3.6447993,-0.9336347,5.89];
ybs[5208]=['',3.6508572,-1.1828422,5.71];
ybs[5209]=['',3.634088,0.5990638,4.74];
ybs[5210]=['',3.6368604,0.2102267,6.04];
ybs[5211]=['4 Cen',3.6417452,-0.5593373,4.73];
ybs[5212]=['',3.6433313,-0.6245479,5.54];
ybs[5213]=['',3.6454984,-0.8246282,6.1];
ybs[5214]=['',3.6447336,-0.6184425,6.19];
ybs[5215]=['7 Boo',3.6406996,0.3108929,5.7];
ybs[5216]=['10 Dra',3.6309277,1.1275323,4.65];
ybs[5217]=['',3.6285565,1.1902211,6.4];
ybs[5218]=['',3.6463227,-0.5007229,6.04];
ybs[5219]=['',3.6402586,0.4979095,5.9];
ybs[5220]=['',3.6512301,-0.9124652,5.71];
ybs[5221]=['ζ Cen',3.6524659,-0.8274176,2.55];
ybs[5222]=['90 Vir',3.6475835,-0.0283188,5.15];
ybs[5223]=['',3.6488871,-0.1427381,6.19];
ybs[5224]=['',3.6562711,-0.946857,6.14];
ybs[5225]=['η Boo',3.6470991,0.319016,2.68];
ybs[5226]=['',3.6572685,-0.9568479,6];
ybs[5227]=['',3.6527954,-0.5481061,6.51];
ybs[5228]=['86 UMa',3.6422877,0.9356506,5.7];
ybs[5229]=['',3.6558684,-0.8152691,5.83];
ybs[5230]=['',3.6790683,-1.3737062,6.09];
ybs[5231]=['',3.6627961,-1.1136111,4.71];
ybs[5232]=['',3.6668576,-1.1505008,6.2];
ybs[5233]=['',3.6522013,0.2432504,6.16];
ybs[5234]=['92 Vir',3.655226,0.0162591,5.91];
ybs[5235]=['',3.6532562,0.5569941,6.32];
ybs[5236]=['',3.6600951,-0.4038945,6.14];
ybs[5237]=['9 Boo',3.6550998,0.4777485,5.01];
ybs[5238]=['φ Cen',3.664199,-0.7368649,3.83];
ybs[5239]=['υ1 Cen',3.6660904,-0.784035,3.87];
ybs[5240]=['47 Hya',3.6647579,-0.4379132,5.15];
ybs[5241]=['',3.6690064,-0.8811835,5.91];
ybs[5242]=['',3.6741438,-1.0751082,6.49];
ybs[5243]=['',3.6772187,-1.1586577,5.97];
ybs[5244]=['',3.6645373,0.2536155,6];
ybs[5245]=['10 Boo',3.6643034,0.3766029,5.76];
ybs[5246]=['',3.6577632,1.071179,6.37];
ybs[5247]=['48 Hya',3.6712374,-0.4385697,5.77];
ybs[5248]=['',3.6699677,-0.0640136,6.4];
ybs[5249]=['',3.67744,-0.7040611,6.13];
ybs[5250]=['υ2 Cen',3.6794309,-0.7979818,4.34];
ybs[5251]=['θ Aps',3.6992746,-1.3423814,5.5];
ybs[5252]=['',3.6763296,0.1531907,5.99];
ybs[5253]=['11 Boo',3.6751742,0.4759349,6.23];
ybs[5254]=['τ Vir',3.677824,0.0249059,4.26];
ybs[5255]=['',3.6816888,-0.4807893,5.48];
ybs[5256]=['',3.6875059,-0.9831517,5.92];
ybs[5257]=['β Cen',3.689525,-1.0557454,0.61];
ybs[5258]=['',3.6846411,-0.5550301,6.18];
ybs[5259]=['',3.686846,-0.725013,6.11];
ybs[5260]=['',3.681503,0.1670142,6.2];
ybs[5261]=['',3.6790401,0.7965037,6.27];
ybs[5262]=['',3.6881314,-0.3933694,6.3];
ybs[5263]=['',3.6858768,0.186223,6.3];
ybs[5264]=['',3.6862752,0.1296699,6.26];
ybs[5265]=['',3.687714,0.0834981,6.24];
ybs[5266]=['',3.689303,-0.0959588,6.39];
ybs[5267]=['',3.6904191,-0.2633395,6.28];
ybs[5268]=['',3.6976347,-0.9561883,6.17];
ybs[5269]=['',3.7122786,-1.3083937,6.02];
ybs[5270]=['',3.6822352,0.8875848,6.15];
ybs[5271]=['',3.7008302,-1.0442558,6.42];
ybs[5272]=['',3.6756246,1.1966173,6.34];
ybs[5273]=['',3.6908025,0.038065,6.28];
ybs[5274]=['',3.6938709,-0.2871444,6.56];
ybs[5275]=['χ Cen',3.6981506,-0.7207471,4.36];
ybs[5276]=['',3.6988169,-0.7541209,6.2];
ybs[5277]=['π Hya',3.6991036,-0.4677214,3.27];
ybs[5278]=['θ Cen',3.7007593,-0.6367982,2.06];
ybs[5279]=['',3.7091416,-1.1052012,6.4];
ybs[5280]=['95 Vir',3.7001702,-0.1645706,5.46];
ybs[5281]=['α Dra',3.6872551,1.1215336,3.65];
ybs[5282]=['',3.7118257,-1.0365818,6.34];
ybs[5283]=['',3.7201737,-1.2290624,6.05];
ybs[5284]=['',3.7105633,-0.7607238,6.17];
ybs[5285]=['',3.7223263,-1.2188349,6.06];
ybs[5286]=['',3.7141016,-0.9009325,6];
ybs[5287]=['',3.7156627,-0.934693,4.75];
ybs[5288]=['96 Vir',3.7102181,-0.1823797,6.47];
ybs[5289]=['',3.7040585,0.7633879,5.27];
ybs[5290]=['13 Boo',3.7053556,0.8611912,5.25];
ybs[5291]=['',3.7183555,-0.2865219,4.91];
ybs[5292]=['',3.7067823,1.0336271,6.46];
ybs[5293]=['η Aps',3.7588894,-1.4158034,4.91];
ybs[5294]=['12 Boo',3.7154404,0.4359298,4.83];
ybs[5295]=['3 UMi',3.6964224,1.29988,6.45];
ybs[5296]=['',3.7506806,-1.357451,6.47];
ybs[5297]=['',3.7209064,0.0217796,6.43];
ybs[5298]=['',3.7303778,-0.93863,5.56];
ybs[5299]=['',3.7254051,-0.4272252,6.34];
ybs[5300]=['',3.7189486,0.5616661,6.11];
ybs[5301]=['',3.732156,-0.9553828,6.11];
ybs[5302]=['50 Hya',3.7270519,-0.4777843,5.08];
ybs[5303]=['',3.7241258,0.0400612,5.01];
ybs[5304]=['',3.7290132,-0.4664565,6.24];
ybs[5305]=['κ Vir',3.7271816,-0.181296,4.19];
ybs[5306]=['',3.7379607,-0.9983104,5.07];
ybs[5307]=['',3.7303815,-0.0167411,5.91];
ybs[5308]=['',3.736045,-0.7321787,5.61];
ybs[5309]=['',3.7394676,-0.9358933,6.39];
ybs[5310]=['',3.7464059,-1.1641401,5.75];
ybs[5311]=['4 UMi',3.7035605,1.3514449,4.82];
ybs[5312]=['',3.733445,-0.1057829,6.36];
ybs[5313]=['14 Boo',3.7318328,0.2242038,5.54];
ybs[5314]=['',3.7369618,-0.5130414,6.08];
ybs[5315]=['',3.7402795,-0.7873838,6.31];
ybs[5316]=['',3.7452954,-1.0476596,6.39];
ybs[5317]=['',3.7884803,-1.4478927,6.42];
ybs[5318]=['κ1 Boo',3.7277436,0.9018815,6.69];
ybs[5319]=['κ2 Boo',3.7278378,0.9019253,4.54];
ybs[5320]=['15 Boo',3.7352163,0.1743112,5.29];
ybs[5321]=['',3.7355352,0.0562496,6.45];
ybs[5322]=['',3.7383123,-0.3196376,5.43];
ybs[5323]=['',3.7342083,0.3797838,6.39];
ybs[5324]=['',3.7198694,1.2098309,5.24];
ybs[5325]=['',3.7322982,0.7226613,6.24];
ybs[5326]=['ε Aps',3.7763817,-1.4000903,5.06];
ybs[5327]=['',3.7426946,-0.5821393,6.55];
ybs[5328]=['ι Vir',3.7406909,-0.1066993,4.08];
ybs[5329]=['δ Oct',3.8012109,-1.4621729,4.32];
ybs[5330]=['α Boo',3.738545,0.3328256,-0.04];
ybs[5331]=['',3.742204,-0.1175424,6.44];
ybs[5332]=['',3.7427487,-0.0577545,6.15];
ybs[5333]=['',3.7403325,0.3281059,5.98];
ybs[5334]=['',3.7455775,-0.3263377,6.22];
ybs[5335]=['',3.735619,0.914948,6.58];
ybs[5336]=['',3.7423792,0.3492174,6.25];
ybs[5337]=['',3.7411324,0.691708,6.38];
ybs[5338]=['',3.7518323,-0.5817633,6.54];
ybs[5339]=['',3.7598279,-1.071362,5.23];
ybs[5340]=['ι Boo',3.7395348,0.8945571,4.75];
ybs[5341]=['λ Boo',3.7407686,0.8024246,4.18];
ybs[5342]=['',3.7465456,0.2644338,5.8];
ybs[5343]=['',3.7494346,-0.1335995,6.47];
ybs[5344]=['ι Lup',3.7567547,-0.8058088,3.55];
ybs[5345]=['',3.7524577,-0.3286119,5.9];
ybs[5346]=['',3.7542927,-0.4525181,5.87];
ybs[5347]=['',3.7563437,-0.647784,5.94];
ybs[5348]=['',3.7614287,-0.9860759,4.33];
ybs[5349]=['λ Vir',3.754379,-0.2353213,4.52];
ybs[5350]=['',3.7447036,0.8935168,6.2];
ybs[5351]=['',3.7482325,0.6177978,4.81];
ybs[5352]=['',3.7599282,-0.7534639,5.56];
ybs[5353]=['',3.7469275,0.8358268,6.32];
ybs[5354]=['',3.7624175,-0.7906069,4.77];
ybs[5355]=['18 Boo',3.7544453,0.2250146,5.41];
ybs[5356]=['υ Vir',3.7559924,-0.0414904,5.14];
ybs[5357]=['ψ Cen',3.761439,-0.663165,4.05];
ybs[5358]=['',3.7565401,0.0047568,6.19];
ybs[5359]=['',3.7521704,0.6746671,6.86];
ybs[5360]=['20 Boo',3.7564618,0.2826618,4.86];
ybs[5361]=['',3.7716323,-1.0222387,4.92];
ybs[5362]=['',3.7513399,0.9556063,6.53];
ybs[5363]=['',3.755945,0.6751328,6.33];
ybs[5364]=['',3.757766,0.5291431,6.44];
ybs[5365]=['',3.7710474,-0.8452772,6.09];
ybs[5366]=['',3.7690648,-0.6090783,5.56];
ybs[5367]=['',3.7742222,-0.8880673,6.02];
ybs[5368]=['',3.7723417,-0.6915456,4.42];
ybs[5369]=['',3.783764,-1.1921449,5.61];
ybs[5370]=['',3.7764265,-0.9300299,6];
ybs[5371]=['51 Hya',3.7721858,-0.4863238,4.77];
ybs[5372]=['',3.7858533,-1.1568521,6.36];
ybs[5373]=['2 Lib',3.7731802,-0.2063762,6.21];
ybs[5374]=['',3.7721106,0.0197445,6.27];
ybs[5375]=['',3.7724641,0.1454671,6.86];
ybs[5376]=['',3.7724713,0.1454961,5.12];
ybs[5377]=['',3.7708618,0.4403045,6.22];
ybs[5378]=['',3.7752535,0.1419561,5.95];
ybs[5379]=['',3.8060869,-1.3410592,6.07];
ybs[5380]=['',3.7795788,-0.4348699,5.32];
ybs[5381]=['',3.7921978,-1.1507053,5.85];
ybs[5382]=['',3.77607,0.0996572,5.1];
ybs[5383]=['',3.7786573,-0.2055928,6.49];
ybs[5384]=['',3.7765226,0.1391895,6.19];
ybs[5385]=['τ1 Lup',3.7861508,-0.79117,4.56];
ybs[5386]=['τ2 Lup',3.7863481,-0.7939283,4.35];
ybs[5387]=['',3.7824373,-0.3504499,6.61];
ybs[5388]=['',3.7863908,-0.7405162,6.32];
ybs[5389]=['',3.7839431,-0.470575,6.48];
ybs[5390]=['',3.788942,-0.6978343,6.35];
ybs[5391]=['',3.7908608,-0.8070992,5.83];
ybs[5392]=['',3.7807573,0.6681715,6.27];
ybs[5393]=['',3.7984379,-1.0350871,6.45];
ybs[5394]=['θ Boo',3.7788165,0.9030518,4.05];
ybs[5395]=['22 Boo',3.7856141,0.3336662,5.39];
ybs[5396]=['104 Vir',3.7904097,-0.1087199,6.17];
ybs[5397]=['52 Hya',3.7944161,-0.5166223,4.97];
ybs[5398]=['',3.8107259,-1.1837624,5.83];
ybs[5399]=['φ Vir',3.7937808,-0.0407829,4.81];
ybs[5400]=['106 Vir',3.7960521,-0.1223302,5.42];
ybs[5401]=['',3.7892219,0.7141197,6.63];
ybs[5402]=['',3.8036717,-0.792895,5.5];
ybs[5403]=['',3.8048025,-0.8661537,5.37];
ybs[5404]=['',3.7943602,0.4918446,7.62];
ybs[5405]=['',3.7944909,0.4918738,7.12];
ybs[5406]=['',3.7929851,0.6298596,6.1];
ybs[5407]=['',3.8069464,-0.7147574,6.39];
ybs[5408]=['',3.8008597,0.0125811,5.94];
ybs[5409]=['',3.8079021,-0.6802809,5.97];
ybs[5410]=['24 Boo',3.7938873,0.8680599,5.59];
ybs[5411]=['',3.8149592,-0.9947457,6.93];
ybs[5412]=['',3.799923,0.5529732,6.06];
ybs[5413]=['',3.7985893,0.7275868,6.35];
ybs[5414]=['',3.8047512,0.0814111,6.02];
ybs[5415]=['σ Lup',3.8147732,-0.8825064,4.42];
ybs[5416]=['',3.8191416,-0.9617624,5.87];
ybs[5417]=['',3.8187859,-0.9213002,5.87];
ybs[5418]=['',3.81624,-0.5379271,6.09];
ybs[5419]=['ρ Boo',3.8086963,0.5282071,3.58];
ybs[5420]=['5 UMi',3.78519,1.3192421,4.25];
ybs[5421]=['',3.8209535,-0.7366359,6.6];
ybs[5422]=['',3.8272427,-1.0493226,6.4];
ybs[5423]=['',3.8110275,0.4637351,6.01];
ybs[5424]=['26 Boo',3.8120583,0.3866414,5.92];
ybs[5425]=['γ Boo',3.8094793,0.6667342,3.03];
ybs[5426]=['',3.8020449,1.1009141,6.09];
ybs[5427]=['',3.8064891,1.0492586,6.27];
ybs[5428]=['',3.8232689,-0.3585839,6.5];
ybs[5429]=['',3.8270084,-0.7264603,5.87];
ybs[5430]=['η Cen',3.8269586,-0.7376402,2.31];
ybs[5431]=['',3.8150245,0.6431999,6.43];
ybs[5432]=['',3.8104132,0.9650034,5.76];
ybs[5433]=['',3.8392419,-1.1874723,6.04];
ybs[5434]=['',3.8307089,-0.8089751,5.55];
ybs[5435]=['',3.8189272,0.5659748,6.33];
ybs[5436]=['',3.8307594,-0.6929443,6.13];
ybs[5437]=['σ Boo',3.8211399,0.517293,4.46];
ybs[5438]=['',3.8207169,0.6373857,6.03];
ybs[5439]=['',3.8322437,-0.7036661,5.74];
ybs[5440]=['',3.8351546,-0.8070243,5.41];
ybs[5441]=['',3.8178888,0.994118,6.48];
ybs[5442]=['',3.8201691,0.859784,5.74];
ybs[5443]=['ρ Lup',3.8377605,-0.8644757,4.05];
ybs[5444]=['',3.8276075,0.4039482,6.38];
ybs[5445]=['',3.83243,-0.2166067,6.2];
ybs[5446]=['',3.8391287,-0.6789156,6.02];
ybs[5447]=['',3.8432602,-0.8148708,6.07];
ybs[5448]=['',3.8444056,-0.858003,6.39];
ybs[5449]=['α1 Cen',3.8462064,-1.0635959,-0.01];
ybs[5450]=['α2 Cen',3.846221,-1.0636007,1.33];
ybs[5451]=['',3.8499172,-0.9868924,6.3];
ybs[5452]=['',3.8369935,0.3175346,5.91];
ybs[5453]=['α Cir',3.8594861,-1.1358324,3.19];
ybs[5454]=['',3.8359392,0.7598631,5.7];
ybs[5455]=['',3.8561679,-1.0248486,6.22];
ybs[5456]=['',3.8507983,-0.6324864,5.67];
ybs[5457]=['',3.8354931,0.9410523,5.85];
ybs[5458]=['33 Boo',3.8386369,0.7731755,5.39];
ybs[5459]=['α Lup',3.8553148,-0.8288875,2.3];
ybs[5460]=['α Aps',3.8878249,-1.3813498,3.83];
ybs[5461]=['',3.8549612,-0.6614283,4];
ybs[5462]=['',3.8461772,0.3817283,6.1];
ybs[5463]=['',3.8479175,0.2344008,5.91];
ybs[5464]=['',3.854214,-0.5416947,6.37];
ybs[5465]=['π1 Boo',3.8479216,0.2847391,4.94];
ybs[5466]=['π2 Boo',3.8479435,0.2847294,5.88];
ybs[5467]=['ζ Boo',3.8498382,0.2377926,4.83];
ybs[5468]=['ζ Boo',3.8498382,0.2377926,4.43];
ybs[5469]=['',3.8093946,1.3884668,6.26];
ybs[5470]=['31 Boo',3.8521677,0.1406395,4.86];
ybs[5471]=['32 Boo',3.852412,0.201707,5.56];
ybs[5472]=['',3.8713684,-1.0991714,5.36];
ybs[5473]=['',3.8529111,0.3668695,6.38];
ybs[5474]=['4 Lib',3.860012,-0.4380858,5.73];
ybs[5475]=['',3.8622637,-0.6156891,4.05];
ybs[5476]=['',3.8692678,-1.0224138,6.11];
ybs[5477]=['μ Vir',3.8587111,-0.1005552,3.88];
ybs[5478]=['',3.8701218,-0.9722194,6.1];
ybs[5479]=['',3.8680749,-0.616,4.92];
ybs[5480]=['34 Boo',3.8593628,0.4612002,4.81];
ybs[5481]=['',4.1163992,-1.5387958,6.48];
ybs[5482]=['',3.8513441,1.0674148,6.25];
ybs[5483]=['',3.8601912,0.7043505,5.73];
ybs[5484]=['',3.8752316,-0.8297775,5.74];
ybs[5485]=['',3.8779111,-0.9160364,5.21];
ybs[5486]=['',3.8679112,-0.026529,6.07];
ybs[5487]=['54 Hya',3.8721455,-0.4458431,4.94];
ybs[5488]=['',3.8786996,-0.9129274,6.07];
ybs[5489]=['',3.8725568,-0.4058744,5.81];
ybs[5490]=['',3.8870596,-1.1640392,5.91];
ybs[5491]=['108 Vir',3.8691979,0.0107358,5.69];
ybs[5492]=['ο Boo',3.8675952,0.2943012,4.6];
ybs[5493]=['5 Lib',3.8716625,-0.2716016,6.33];
ybs[5494]=['',3.8727914,-0.3713698,6.4];
ybs[5495]=['ε Boo',3.8661592,0.4707617,5.12];
ybs[5496]=['ε Boo',3.8661592,0.4707472,2.7];
ybs[5497]=['',3.8679875,0.327817,6.13];
ybs[5498]=['',3.8773765,-0.6700712,5.94];
ybs[5499]=['',3.8796009,-0.7619886,6.3];
ybs[5500]=['',3.8670035,0.5704797,6.28];
ybs[5501]=['109 Vir',3.872416,0.0312583,3.72];
ybs[5502]=['',3.8713931,0.262324,5.63];
ybs[5503]=['',3.8773662,-0.3739562,6.06];
ybs[5504]=['55 Hya',3.878149,-0.4489994,5.63];
ybs[5505]=['',3.8874254,-0.9907943,6.23];
ybs[5506]=['56 Hya',3.8797881,-0.4570787,5.24];
ybs[5507]=['57 Hya',3.8807318,-0.4668316,5.77];
ybs[5508]=['',3.8801047,-0.2258604,6.35];
ybs[5509]=['',3.8840753,-0.6411559,6.04];
ybs[5510]=['',3.908376,-1.2791299,5.6];
ybs[5511]=['',3.8865595,-0.4250264,5.68];
ybs[5512]=['',3.8840694,-0.0165549,6.14];
ybs[5513]=['μ Lib',3.8862693,-0.2486999,5.31];
ybs[5514]=['',3.8810814,0.423516,6.14];
ybs[5515]=['π1 Oct',3.9550591,-1.4542515,5.65];
ybs[5516]=['58 Hya',3.8909544,-0.4897469,4.41];
ybs[5517]=['',3.9034413,-1.1154243,5.87];
ybs[5518]=['ο Lup',3.8975431,-0.7622748,4.32];
ybs[5519]=['',3.8837172,0.6581704,6.16];
ybs[5520]=['α1 Lib',3.8922983,-0.2809497,5.15];
ybs[5521]=['α2 Lib',3.8931364,-0.2817241,2.75];
ybs[5522]=['',3.8878425,0.4976889,5.8];
ybs[5523]=['38 Boo',3.884166,0.8031211,5.74];
ybs[5524]=['',3.8892715,0.4155927,5.85];
ybs[5525]=['11 Lib',3.8933391,-0.0418715,4.94];
ybs[5526]=['',3.8932159,-0.0062331,6.18];
ybs[5527]=['',3.8848136,0.8949025,6.51];
ybs[5528]=['39 Boo',3.885647,0.8485797,5.69];
ybs[5529]=['ζ Cir',3.9130906,-1.1534812,6.09];
ybs[5530]=['',3.9303329,-1.3397062,5.34];
ybs[5531]=['',3.8897618,0.6487698,5.48];
ybs[5532]=['',3.9009373,-0.5354003,6.29];
ybs[5533]=['',3.9025507,-0.6615221,5.03];
ybs[5534]=['ξ Boo',3.8943338,0.3316359,4.55];
ybs[5535]=['π2 Oct',3.9676193,-1.4509231,5.65];
ybs[5536]=['',3.9160596,-1.0508987,5.2];
ybs[5537]=['',3.9407962,-1.3483718,5.93];
ybs[5538]=['12 Lib',3.9085087,-0.4318072,5.3];
ybs[5539]=['',3.9101386,-0.5829212,5.82];
ybs[5540]=['',3.9031619,0.2723672,6.4];
ybs[5541]=['θ Cir',3.9214633,-1.0974321,5.11];
ybs[5542]=['',3.8923947,1.0331308,5.46];
ybs[5543]=['',3.9030827,0.3325519,6.01];
ybs[5544]=['ξ1 Lib',3.9083036,-0.2093844,5.8];
ybs[5545]=['',3.9385321,-1.311238,6.2];
ybs[5546]=['',3.9184694,-0.9234075,5.38];
ybs[5547]=['ω Oct',4.0007882,-1.4813983,5.91];
ybs[5548]=['',3.9150179,-0.5926046,5.32];
ybs[5549]=['',3.9191869,-0.8373516,5.64];
ybs[5550]=['',3.9215781,-0.8996173,6.64];
ybs[5551]=['',3.91899,-0.6896433,6.36];
ybs[5552]=['',3.9183307,-0.5713206,6.06];
ybs[5553]=['β UMi',3.8862464,1.2925082,2.08];
ybs[5554]=['ξ2 Lib',3.9187064,-0.2008393,5.46];
ybs[5555]=['',3.9213015,-0.5105974,6.29];
ybs[5556]=['',3.9263051,-0.8545119,6.35];
ybs[5557]=['',3.9155505,0.2504304,5.77];
ybs[5558]=['',3.9220727,-0.3754687,5.74];
ybs[5559]=['',3.9138763,0.5620373,6.12];
ybs[5560]=['16 Lib',3.920308,-0.077558,4.49];
ybs[5561]=['β Lup',3.9276382,-0.7545166,2.68];
ybs[5562]=['',3.9278259,-0.6981904,6.15];
ybs[5563]=['',3.9218078,-0.00462,5.53];
ybs[5564]=['',3.9189853,0.3745098,6.49];
ybs[5565]=['',3.9197421,0.2843308,5.71];
ybs[5566]=['κ Cen',3.9303355,-0.7365399,3.13];
ybs[5567]=['59 Hya',3.9274755,-0.4843975,5.65];
ybs[5568]=['17 Lib',3.9250494,-0.196383,6.6];
ybs[5569]=['',3.9304292,-0.6628383,6.47];
ybs[5570]=['',3.9316603,-0.7549655,6.1];
ybs[5571]=['',3.9147405,0.8644759,5.63];
ybs[5572]=['18 Lib',3.9279736,-0.196189,5.87];
ybs[5573]=['',3.927731,-0.0887642,6.09];
ybs[5574]=['',3.9296534,0.0780395,5.93];
ybs[5575]=['',3.9391102,-0.6659118,5.89];
ybs[5576]=['δ Lib',3.9369722,-0.150354,4.92];
ybs[5577]=['',3.9422249,-0.6013389,6.22];
ybs[5578]=['40 Boo',3.9294157,0.6836257,5.64];
ybs[5579]=['',3.9182014,1.1490397,4.6];
ybs[5580]=['',3.9383561,-0.0497526,5.52];
ybs[5581]=['60 Hya',3.9425721,-0.4914114,5.85];
ybs[5582]=['',3.9355963,0.3830947,6.38];
ybs[5583]=['η Cir',3.9569788,-1.1192033,5.17];
ybs[5584]=['',3.9403948,-0.0041185,5.71];
ybs[5585]=['',3.9466067,-0.5713892,5.44];
ybs[5586]=['',3.8781913,1.438347,5.64];
ybs[5587]=['',3.9334561,0.8234774,6.37];
ybs[5588]=['',3.9689973,-1.2566044,6.52];
ybs[5589]=['',3.9445536,-0.054566,6.61];
ybs[5590]=['ω Boo',3.9408767,0.4348091,4.81];
ybs[5591]=['110 Vir',3.9450625,0.0348443,4.4];
ybs[5592]=['β Boo',3.9395192,0.7032827,3.5];
ybs[5593]=['σ Lib',3.9510418,-0.442901,3.29];
ybs[5594]=['',3.9545362,-0.7148081,6.41];
ybs[5595]=['π Lup',3.9566679,-0.8228358,4.72];
ybs[5596]=['π Lup',3.9566679,-0.8228358,4.82];
ybs[5597]=['',3.9571833,-0.718396,5.15];
ybs[5598]=['',3.9358476,1.0490953,5.93];
ybs[5599]=['',3.9447994,0.612801,5.51];
ybs[5600]=['',3.9502248,0.0942139,6.5];
ybs[5601]=['',3.9709484,-1.1408893,6.17];
ybs[5602]=['',3.9443545,0.7775354,6.65];
ybs[5603]=['',3.9470359,0.6016395,6.59];
ybs[5604]=['',3.9585933,-0.4517504,6.67];
ybs[5605]=['',3.9609349,-0.6345602,6.27];
ybs[5606]=['ψ Boo',3.9509925,0.4686762,4.54];
ybs[5607]=['',3.9669114,-0.8583743,5.77];
ybs[5608]=['44 Boo',3.9471294,0.8300748,4.76];
ybs[5609]=['',3.9621132,-0.5412652,5.96];
ybs[5610]=['',3.9613243,-0.3861599,6.17];
ybs[5611]=['',3.9780853,-1.1724433,5.76];
ybs[5612]=['ν Lib',3.9618919,-0.285366,5.2];
ybs[5613]=['',3.9771502,-1.1123808,6.28];
ybs[5614]=['',3.9697563,-0.7099384,5.79];
ybs[5615]=['',3.9718481,-0.7497963,5.85];
ybs[5616]=['λ Lup',3.972829,-0.791891,4.05];
ybs[5617]=['47 Boo',3.9542392,0.8387555,5.57];
ybs[5618]=['',3.9930497,-1.2716595,6.01];
ybs[5619]=['',3.9459374,1.1488636,6.13];
ybs[5620]=['',3.9599232,0.6346389,6.35];
ybs[5621]=['',3.9657775,0.0943378,6.16];
ybs[5622]=['',3.9825973,-1.0736198,6.3];
ybs[5623]=['',3.9639179,0.3202435,6.02];
ybs[5624]=['45 Boo',3.9635246,0.432424,4.93];
ybs[5625]=['',3.9574362,0.9505547,5.25];
ybs[5626]=['',3.9780668,-0.6786582,5.98];
ybs[5627]=['',3.9841955,-0.9675634,5.54];
ybs[5628]=['46 Boo',3.9682452,0.4574245,5.67];
ybs[5629]=['',3.9708494,0.229382,6.1];
ybs[5630]=['',3.9691401,0.4366129,5.81];
ybs[5631]=['',3.9783298,-0.461194,5.76];
ybs[5632]=['',3.9847965,-0.7918312,6.44];
ybs[5633]=['',3.9845781,-0.7918655,7.39];
ybs[5634]=['',3.9999059,-1.2246826,5.81];
ybs[5635]=['',3.9925897,-1.0792114,6.32];
ybs[5636]=['κ1 Lup',3.9865721,-0.8522215,3.87];
ybs[5637]=['κ2 Lup',3.9866818,-0.8523231,5.69];
ybs[5638]=['',3.966713,0.8720065,6.39];
ybs[5639]=['ζ Lup',3.9883599,-0.9108858,3.41];
ybs[5640]=['',3.9890997,-0.8431558,6.33];
ybs[5641]=['',3.990182,-0.7782613,4.82];
ybs[5642]=['ι1 Lib',3.9864375,-0.3470155,4.54];
ybs[5643]=['',3.9910326,-0.6314918,6.1];
ybs[5644]=['',3.9845342,0.3296027,5.89];
ybs[5645]=['',3.9912616,-0.4206018,6.47];
ybs[5646]=['ι2 Lib',3.9912275,-0.3444909,6.08];
ybs[5647]=['23 Lib',3.9921257,-0.4433041,6.45];
ybs[5648]=['',3.99395,-0.4587374,5.84];
ybs[5649]=['',3.9873726,0.3350183,6.68];
ybs[5650]=['1 Lup',3.9973701,-0.5516799,4.91];
ybs[5651]=['',4.0081841,-1.0645226,5.73];
ybs[5652]=['26 Lib',3.9965828,-0.3116884,6.17];
ybs[5653]=['',4.0038416,-0.8406084,5.95];
ybs[5654]=['δ Cir',4.0096716,-1.0654556,5.09];
ybs[5655]=['',3.9907582,0.3995579,6.3];
ybs[5656]=['ε Cir',4.0131401,-1.1117542,4.86];
ybs[5657]=['',4.0041708,-0.7257114,5.16];
ybs[5658]=['',4.0047598,-0.7605056,6.04];
ybs[5659]=['',3.9974175,-0.0976074,6.28];
ybs[5660]=['β Cir',4.0118859,-1.0278152,4.07];
ybs[5661]=['γ TrA',4.0196646,-1.200212,2.89];
ybs[5662]=['',3.9750082,1.1814068,6.17];
ybs[5663]=['',3.9903776,0.6662684,6.2];
ybs[5664]=['',3.9928943,0.5532335,5.99];
ybs[5665]=['3 Ser',3.998595,0.0846462,5.33];
ybs[5666]=['χ Boo',3.9946828,0.5074411,5.26];
ybs[5667]=['',3.9926891,0.734457,6.13];
ybs[5668]=['',4.0047008,-0.3924973,5.5];
ybs[5669]=['4 Ser',4.00148,0.0049382,5.63];
ybs[5670]=['',4.0177921,-1.0573929,5.46];
ybs[5671]=['δ Boo',3.9989371,0.5798898,3.47];
ybs[5672]=['',4.0132783,-0.7181851,6.28];
ybs[5673]=['μ Lup',4.0153518,-0.8371114,4.27];
ybs[5674]=['',4.0271333,-1.1792881,6.28];
ybs[5675]=['β Lib',4.0069756,-0.1653139,2.61];
ybs[5676]=['2 Lup',4.0113307,-0.5277389,4.34];
ybs[5677]=['',4.0166857,-0.7134229,5.59];
ybs[5678]=['',4.0151176,-0.5462422,6.18];
ybs[5679]=['',4.0190538,-0.6489867,6.2];
ybs[5680]=['',4.0129311,-0.0095904,5.89];
ybs[5681]=['',3.9920663,1.1738497,5.13];
ybs[5682]=['',4.0121111,0.3575242,5.7];
ybs[5683]=['',3.9926928,1.2017522,6.51];
ybs[5684]=['5 Ser',4.0166853,0.0292792,5.06];
ybs[5685]=['δ Lup',4.0273018,-0.7109455,3.22];
ybs[5686]=['',4.0282542,-0.7127278,6.2];
ybs[5687]=['',4.0277388,-0.6685622,6.48];
ybs[5688]=['ν1 Lup',4.0311119,-0.8380037,5];
ybs[5689]=['ν2 Lup',4.0296759,-0.8448132,5.65];
ybs[5690]=['',4.0368815,-1.0601601,5.67];
ybs[5691]=['28 Lib',4.0242447,-0.3184451,6.17];
ybs[5692]=['',4.0164087,0.5659685,6.32];
ybs[5693]=['ο Lib',4.0247135,-0.2728862,6.3];
ybs[5694]=['γ Cir',4.0376048,-1.0368391,4.51];
ybs[5695]=['φ1 Lup',4.0289719,-0.6343901,3.56];
ybs[5696]=['',4.0232402,-0.0436395,6.35];
ybs[5697]=['',4.024841,-0.1031766,5.54];
ybs[5698]=['ε Lup',4.0332671,-0.7814799,3.37];
ybs[5699]=['ο CrB',4.0193172,0.5153737,5.51];
ybs[5700]=['6 Ser',4.0242235,0.010967,5.35];
ybs[5701]=['',4.0237493,0.4340783,6.39];
ybs[5702]=['φ2 Lup',4.0348987,-0.6448027,4.54];
ybs[5703]=['',4.0517021,-1.1936903,5.89];
ybs[5704]=['11 UMi',4.0015948,1.2520101,5.02];
ybs[5705]=['',4.0177836,0.9053222,5.66];
ybs[5706]=['',4.020969,0.7740012,6.19];
ybs[5707]=['7 Ser',4.0297407,0.2178378,6.28];
ybs[5708]=['50 Boo',4.0264447,0.573293,5.37];
ybs[5709]=['υ Lup',4.0420021,-0.6945607,5.37];
ybs[5710]=['',4.0370303,-0.2173812,5.72];
ybs[5711]=['8 Ser',4.0360387,-0.0193412,6.12];
ybs[5712]=['',4.0434992,-0.6676652,7.03];
ybs[5713]=['ε Lib',4.0383944,-0.181648,4.94];
ybs[5714]=['',4.0445193,-0.6775098,4.6];
ybs[5715]=['',4.0566411,-1.1277457,5.71];
ybs[5716]=['',4.0296606,0.68932,5.5];
ybs[5717]=['η CrB',4.0326431,0.5271209,5.58];
ybs[5718]=['η CrB',4.0326431,0.5271209,6.08];
ybs[5719]=['ρ Oct',4.1419582,-1.4755088,5.57];
ybs[5720]=['κ1 Aps',4.0763809,-1.2823128,5.49];
ybs[5721]=['',4.0277096,1.0814209,5.98];
ybs[5722]=['',4.0356924,0.7886358,6.01];
ybs[5723]=['μ1 Boo',4.0379117,0.6508652,4.31];
ybs[5724]=['μ2 Boo',4.0380222,0.6503515,6.5];
ybs[5725]=['γ UMi',4.0173562,1.2522134,3.05];
ybs[5726]=['',4.0530074,-0.6431837,5.45];
ybs[5727]=['',4.0275924,1.1040084,5.79];
ybs[5728]=['',4.059002,-0.9020008,6.1];
ybs[5729]=['τ1 Ser',4.0444866,0.2677914,5.17];
ybs[5730]=['',4.0447705,0.3385263,6.27];
ybs[5731]=['',4.0459076,0.5977976,5.46];
ybs[5732]=['',4.0627845,-0.8170878,5.24];
ybs[5733]=['ζ1 Lib',4.0563295,-0.2932138,5.64];
ybs[5734]=['ι Dra',4.0381648,1.0276644,3.29];
ybs[5735]=['',4.0522308,0.4366425,6.02];
ybs[5736]=['10 Ser',4.0573616,0.0306975,5.17];
ybs[5737]=['β CrB',4.0528501,0.5065297,3.68];
ybs[5738]=['',4.0457098,0.9413521,6.45];
ybs[5739]=['',4.0667288,-0.363216,6.22];
ybs[5740]=['ζ3 Lib',4.066875,-0.2913274,5.82];
ybs[5741]=['',4.0739417,-0.6755249,6.25];
ybs[5742]=['',4.0558483,0.8223633,6.15];
ybs[5743]=['',4.0726243,-0.5753113,6.46];
ybs[5744]=['',4.0496848,1.085447,6.5];
ybs[5745]=['',4.0506833,1.0574313,5.9];
ybs[5746]=['',4.0715897,-0.3533699,6.22];
ybs[5747]=['',4.1131382,-1.3612837,6.18];
ybs[5748]=['',4.0671134,0.148298,6.57];
ybs[5749]=['',4.0560675,0.9618792,6.43];
ybs[5750]=['',4.0638712,0.5446036,6.46];
ybs[5751]=['',4.0639832,0.6409124,6.37];
ybs[5752]=['',4.0754492,-0.3447377,5.52];
ybs[5753]=['ν1 Boo',4.0657902,0.7112334,5.02];
ybs[5754]=['ζ4 Lib',4.0766909,-0.2955557,5.5];
ybs[5755]=['',4.0779943,-0.4288384,7];
ybs[5756]=['',4.0950899,-1.1465561,6.51];
ybs[5757]=['',4.0825554,-0.7006946,5.82];
ybs[5758]=['',4.0569923,1.0823919,6.38];
ybs[5759]=['',4.0679519,0.6376425,6.38];
ybs[5760]=['τ2 Ser',4.0722546,0.2788054,6.22];
ybs[5761]=['ε TrA',4.0971159,-1.1588326,4.11];
ybs[5762]=['11 Ser',4.0763478,-0.0221258,5.51];
ybs[5763]=['',4.0839199,-0.6881837,6.36];
ybs[5764]=['ν2 Boo',4.0695052,0.7123992,5.02];
ybs[5765]=['36 Lib',4.0845611,-0.4909161,5.15];
ybs[5766]=['γ Lup',4.0874906,-0.7198983,2.78];
ybs[5767]=['37 Lib',4.0819429,-0.1770667,4.62];
ybs[5768]=['θ CrB',4.0749918,0.5458999,4.14];
ybs[5769]=['',4.0825282,-0.1008042,6.51];
ybs[5770]=['',4.0830691,-0.1616862,5.17];
ybs[5771]=['',4.0909852,-0.7860689,4.54];
ybs[5772]=['κ2 Aps',4.1151563,-1.2832368,5.65];
ybs[5773]=['',4.0797098,0.2976983,6.45];
ybs[5774]=['',4.0923241,-0.7762634,5.43];
ybs[5775]=['',4.0635045,1.1192116,5.79];
ybs[5776]=['',4.1233108,-1.3292163,5.95];
ybs[5777]=['γ Lib',4.0879944,-0.259522,3.91];
ybs[5778]=['δ Ser',4.0839646,0.1825097,3.8];
ybs[5779]=['δ Ser',4.0839645,0.1825388,3.8];
ybs[5780]=['',4.0916394,-0.5789691,6.24];
ybs[5781]=['',4.0854659,0.0277258,6.56];
ybs[5782]=['',4.1132992,-1.227059,6.44];
ybs[5783]=['α CrB',4.0828611,0.464854,2.23];
ybs[5784]=['υ Lib',4.0950593,-0.4924327,3.58];
ybs[5785]=['τ3 Ser',4.0869888,0.3067493,6.12];
ybs[5786]=['',4.0886867,0.1952257,6.07];
ybs[5787]=['ω Lup',4.1002963,-0.7443179,4.33];
ybs[5788]=['',4.1044071,-0.915445,5.44];
ybs[5789]=['14 Ser',4.0920205,-0.0111921,6.51];
ybs[5790]=['μ CrB',4.0847206,0.6794511,5.11];
ybs[5791]=['',4.0969447,-0.4600531,6.19];
ybs[5792]=['16 Ser',4.0913623,0.1733173,5.26];
ybs[5793]=['',4.1101285,-1.046955,5.95];
ybs[5794]=['τ5 Ser',4.0911233,0.2799373,5.93];
ybs[5795]=['',4.1022414,-0.6848564,6.57];
ybs[5796]=['',4.098242,-0.4052763,5.78];
ybs[5797]=['',4.1029384,-0.684283,6.04];
ybs[5798]=['',4.0872471,0.6683537,6.42];
ybs[5799]=['',4.1004623,-0.4936731,6.32];
ybs[5800]=['',4.1002099,-0.3681744,5.84];
ybs[5801]=['',4.0837356,0.9397128,5.97];
ybs[5802]=['τ Lib',4.1022569,-0.5210907,3.66];
ybs[5803]=['',4.0923182,0.522056,6.52];
ybs[5804]=['41 Lib',4.1029231,-0.3382512,5.38];
ybs[5805]=['',4.1014875,-0.1548631,6.5];
ybs[5806]=['',4.1014947,-0.1548049,6.48];
ybs[5807]=['',4.0873767,0.9073923,6.74];
ybs[5808]=['',4.0866272,0.9520861,5.74];
ybs[5809]=['',4.105024,-0.4054134,6.34];
ybs[5810]=['ψ1 Lup',4.1073246,-0.6019625,4.67];
ybs[5811]=['',4.1134032,-0.8344923,6.23];
ybs[5812]=['',4.1093098,-0.546137,6.34];
ybs[5813]=['φ Boo',4.0958741,0.7029183,5.24];
ybs[5814]=['42 Lib',4.1090968,-0.4170604,4.96];
ybs[5815]=['',4.1141295,-0.7808314,4.64];
ybs[5816]=['θ UMi',4.0611273,1.3485622,4.96];
ybs[5817]=['',4.1004881,0.6038214,6.11];
ybs[5818]=['',4.0935272,0.9499759,5.87];
ybs[5819]=['',4.0486867,1.4026318,6.58];
ybs[5820]=['',4.0973892,0.8153985,5.75];
ybs[5821]=['',4.1073432,0.2090065,6.25];
ybs[5822]=['',4.1207443,-0.8650892,6.04];
ybs[5823]=['ζ1 CrB',4.1028065,0.6380636,6];
ybs[5824]=['ζ2 CrB',4.1028429,0.6380491,5.07];
ybs[5825]=['',4.0984224,0.8786787,5.84];
ybs[5826]=['',4.1275729,-1.0535336,6.48];
ybs[5827]=['',4.1200245,-0.6545255,5.24];
ybs[5828]=['κ Lib',4.1162012,-0.3448042,4.74];
ybs[5829]=['ψ2 Lup',4.120078,-0.6071493,4.75];
ybs[5830]=['τ6 Ser',4.1107466,0.278332,6.01];
ybs[5831]=['',4.1002289,1.009602,6.45];
ybs[5832]=['ι Ser',4.1130732,0.3419636,4.52];
ybs[5833]=['χ Ser',4.1143666,0.222886,5.33];
ybs[5834]=['',4.0916474,1.2078373,5.62];
ybs[5835]=['τ7 Ser',4.1146895,0.3209113,5.81];
ybs[5836]=['',4.127874,-0.7312078,5.94];
ybs[5837]=['',4.1224355,-0.263886,6.31];
ybs[5838]=['η Lib',4.1253324,-0.2748662,5.41];
ybs[5839]=['γ CrB',4.1180014,0.4576064,3.84];
ybs[5840]=['',4.1203854,0.2372143,6.48];
ybs[5841]=['',4.1458334,-1.1434747,6.18];
ybs[5842]=['',4.1458407,-1.1434747,6.39];
ybs[5843]=['ψ Ser',4.1245014,0.0425694,5.88];
ybs[5844]=['α Ser',4.1254054,0.1108234,2.65];
ybs[5845]=['π CrB',4.1231471,0.5661813,5.56];
ybs[5846]=['',4.1351764,-0.4910741,6.51];
ybs[5847]=['',4.116867,0.9125306,5.51];
ybs[5848]=['τ8 Ser',4.1269034,0.2999963,6.14];
ybs[5849]=['',4.1303404,0.0937532,5.58];
ybs[5850]=['',4.1377786,-0.6066245,5.61];
ybs[5851]=['',4.1316768,0.0142462,6.33];
ybs[5852]=['',4.1410834,-0.7028149,6.42];
ybs[5853]=['25 Ser',4.1336539,-0.0328011,5.4];
ybs[5854]=['',4.1412135,-0.6630598,6.01];
ybs[5855]=['',4.1481797,-0.9164978,6.07];
ybs[5856]=['',4.1366994,-0.1081207,6.24];
ybs[5857]=['β Ser',4.133457,0.2678563,3.67];
ybs[5858]=['λ Ser',4.1348615,0.1270301,4.43];
ybs[5859]=['',4.153834,-0.9299499,5.77];
ybs[5860]=['υ Ser',4.1383037,0.2450603,5.71];
ybs[5861]=['',4.1527651,-0.8549461,5.84];
ybs[5862]=['',4.1538664,-0.7936777,6.12];
ybs[5863]=['',4.1583875,-0.9621665,5.73];
ybs[5864]=['',4.1423868,0.2393671,6];
ybs[5865]=['',4.1461751,-0.06793,5.53];
ybs[5866]=['',4.1678774,-1.1383688,6.54];
ybs[5867]=['',4.14082,0.5526028,6.44];
ybs[5868]=['',4.1328245,0.9669103,5.92];
ybs[5869]=['κ Ser',4.1444831,0.3153465,4.09];
ybs[5870]=['',4.1433409,0.4901394,5.85];
ybs[5871]=['μ Ser',4.1490997,-0.0611464,3.53];
ybs[5872]=['',4.1594647,-0.8226196,6.01];
ybs[5873]=['χ Lup',4.1561611,-0.5881695,3.95];
ybs[5874]=['',4.1693118,-1.0939371,6.19];
ybs[5875]=['1 Sco',4.1558795,-0.4507105,4.64];
ybs[5876]=['',4.1322587,1.0912599,5.19];
ybs[5877]=['',4.1374224,0.965208,5.86];
ybs[5878]=['ω Ser',4.1518327,0.0370631,5.23];
ybs[5879]=['δ CrB',4.1478879,0.4537001,4.63];
ybs[5880]=['',4.1655819,-0.8846491,6.6];
ybs[5881]=['κ TrA',4.1798321,-1.1985689,5.09];
ybs[5882]=['ε Ser',4.1540429,0.0768854,3.71];
ybs[5883]=['',4.1614586,-0.5228736,6.4];
ybs[5884]=['',4.153125,0.2628633,5.2];
ybs[5885]=['36 Ser',4.1562441,-0.0552028,5.11];
ybs[5886]=['',4.158302,-0.2479366,6.19];
ybs[5887]=['β TrA',4.1771599,-1.108296,2.85];
ybs[5888]=['',4.1755616,-1.061398,6.15];
ybs[5889]=['ρ Ser',4.1553834,0.3648681,4.76];
ybs[5890]=['',4.1783718,-1.0515214,5.77];
ybs[5891]=['κ CrB',4.1545649,0.6210766,4.82];
ybs[5892]=['λ Lib',4.1659348,-0.3532278,5.03];
ybs[5893]=['ζ UMi',4.1156255,1.3564349,4.32];
ybs[5894]=['2 Sco',4.1673579,-0.4432842,4.59];
ybs[5895]=['',4.1808723,-1.0568395,5.76];
ybs[5896]=['',4.168568,-0.4294209,5.39];
ybs[5897]=['',4.1686897,-0.419734,5.42];
ybs[5898]=['θ Lib',4.1679381,-0.2932231,4.15];
ybs[5899]=['',4.1628157,0.302502,6.36];
ybs[5900]=['',4.1713254,-0.4783816,6.14];
ybs[5901]=['39 Ser',4.1641316,0.2290795,6.1];
ybs[5902]=['3 Sco',4.1719239,-0.4418156,5.87];
ybs[5903]=['',4.1656802,0.2793191,6.09];
ybs[5904]=['χ Her',4.160458,0.7396695,4.62];
ybs[5905]=['47 Lib',4.1731932,-0.3395267,5.94];
ybs[5906]=['',4.1758961,-0.543735,6.21];
ybs[5907]=['4 Sco',4.1756507,-0.4596492,5.62];
ybs[5908]=['',4.1789989,-0.6969788,6.03];
ybs[5909]=['40 Ser',4.170714,0.1485215,6.29];
ybs[5910]=['',4.1942046,-1.1363128,5.75];
ybs[5911]=['',4.1837492,-0.841798,6.31];
ybs[5912]=['',4.1575957,0.9731026,5.81];
ybs[5913]=['',4.1790963,-0.5559845,6.29];
ybs[5914]=['',4.1698533,0.353257,5.44];
ybs[5915]=['ξ1 Lup',4.182089,-0.5940364,5.12];
ybs[5916]=['ξ2 Lup',4.1821398,-0.5939975,5.62];
ybs[5917]=['',4.1783896,-0.2525355,6.37];
ybs[5918]=['ρ Sco',4.1818251,-0.5110948,3.88];
ybs[5919]=['',4.1842311,-0.6327542,5.8];
ybs[5920]=['',4.1797882,-0.2600328,6.13];
ybs[5921]=['',4.1746604,0.3237658,6.26];
ybs[5922]=['2 Her',4.1689405,0.7516762,5.37];
ybs[5923]=['γ Ser',4.1782218,0.2721306,3.85];
ybs[5924]=['',4.1849061,-0.3674281,5.85];
ybs[5925]=['',4.1893565,-0.6557482,6.31];
ybs[5926]=['λ CrB',4.174324,0.6610751,5.45];
ybs[5927]=['',4.1966824,-0.9440289,6.1];
ybs[5928]=['4 Her',4.1728065,0.741692,5.75];
ybs[5929]=['',4.2036264,-1.1142824,6.41];
ybs[5930]=['φ Ser',4.1817154,0.2503696,5.54];
ybs[5931]=['48 Lib',4.1868962,-0.2504237,4.88];
ybs[5932]=['',4.1890272,-0.4345859,5.43];
ybs[5933]=['',4.193953,-0.7297652,4.99];
ybs[5934]=['π Sco',4.1902698,-0.4569721,2.89];
ybs[5935]=['',4.1958998,-0.7107079,6.49];
ybs[5936]=['',4.2020152,-0.9537335,6.13];
ybs[5937]=['ε CrB',4.182708,0.4678983,4.15];
ybs[5938]=['η Lup',4.1964398,-0.6713349,3.41];
ybs[5939]=['',4.1726868,1.0269771,6.31];
ybs[5940]=['',4.1816534,0.6916044,6.31];
ybs[5941]=['',4.2108186,-1.0927125,6.25];
ybs[5942]=['',4.1999409,-0.7069034,6.21];
ybs[5943]=['δ Sco',4.1965864,-0.3960033,2.32];
ybs[5944]=['49 Lib',4.1963095,-0.2897422,5.47];
ybs[5945]=['',4.2267726,-1.2647568,5.7];
ybs[5946]=['',4.2013418,-0.5577471,6.33];
ybs[5947]=['',4.1882092,0.6383609,5.62];
ybs[5948]=['',4.1911016,0.4512041,6.33];
ybs[5949]=['50 Lib',4.198038,-0.1479836,5.55];
ybs[5950]=['',4.1817271,0.9543554,4.95];
ybs[5951]=['ι1 Nor',4.2129907,-1.009518,4.63];
ybs[5952]=['η Nor',4.2107155,-0.8603737,4.65];
ybs[5953]=['',4.1978263,0.0760974,5.83];
ybs[5954]=['',4.1877858,0.8693941,6.05];
ybs[5955]=['',4.2070217,-0.5096756,6.03];
ybs[5956]=['5 Her',4.1990083,0.3098144,5.12];
ybs[5957]=['',4.2107842,-0.6748928,4.89];
ybs[5958]=['ρ CrB',4.1974711,0.5800811,5.41];
ybs[5959]=['',4.2098684,-0.4525876,5];
ybs[5960]=['',4.2111574,-0.5596661,6.01];
ybs[5961]=['ι CrB',4.1993797,0.5198274,4.99];
ybs[5962]=['π Ser',4.2034115,0.3968478,4.83];
ybs[5963]=['',4.2122989,-0.4327051,6.21];
ybs[5964]=['',4.2143854,-0.5808457,6.1];
ybs[5965]=['',4.2160176,-0.6619762,5.9];
ybs[5966]=['43 Ser',4.210502,0.085883,6.08];
ybs[5967]=['ξ Sco',4.2137508,-0.1996419,5.07];
ybs[5968]=['ξ Sco',4.2137508,-0.1996419,4.77];
ybs[5969]=['',4.2297032,-0.9818388,6.16];
ybs[5970]=['δ Nor',4.2247011,-0.7895469,4.72];
ybs[5971]=['',4.2006552,0.9223869,5.93];
ybs[5972]=['υ Her',4.2043272,0.8023298,4.76];
ybs[5973]=['',4.2072247,0.638187,5.83];
ybs[5974]=['β1 Sco',4.2187462,-0.3468067,2.62];
ybs[5975]=['β2 Sco',4.2187679,-0.3467436,4.92];
ybs[5976]=['θ Dra',4.1991222,1.0209858,4.01];
ybs[5977]=['θ Lup',4.224603,-0.643443,4.23];
ybs[5978]=['',4.2218248,-0.4131375,5.92];
ybs[5979]=['',4.2195541,-0.1109426,6.53];
ybs[5980]=['',4.2206612,-0.1082835,6.41];
ybs[5981]=['',4.2275772,-0.6426224,5.73];
ybs[5982]=['',4.2185382,0.1401699,6.29];
ybs[5983]=['ω1 Sco',4.2247624,-0.361867,3.96];
ybs[5984]=['ι2 Nor',4.238269,-1.012243,5.57];
ybs[5985]=['',4.2045233,1.0357549,6.19];
ybs[5986]=['',4.2255875,-0.2467024,6.32];
ybs[5987]=['ω2 Sco',4.2273825,-0.3653425,4.32];
ybs[5988]=['',4.2295457,-0.4280536,6.33];
ybs[5989]=['',4.2333789,-0.6836207,7.05];
ybs[5990]=['',4.2333926,-0.6834073,6.65];
ybs[5991]=['',4.2307721,-0.4605966,5.38];
ybs[5992]=['11 Sco',4.2279365,-0.2235669,5.78];
ybs[5993]=['',4.2332803,-0.4144951,5.88];
ybs[5994]=['45 Ser',4.227174,0.1715265,5.63];
ybs[5995]=['',4.2255868,0.379756,6.14];
ybs[5996]=['',4.2372,-0.5709416,6.19];
ybs[5997]=['',4.2387682,-0.5865783,5.54];
ybs[5998]=['κ Her',4.2288589,0.2964135,5];
ybs[5999]=['κ Her',4.2288877,0.2965444,6.25];
ybs[6000]=['47 Ser',4.230897,0.1478414,5.73];
ybs[6001]=['',4.2333329,0.0591884,5.91];
ybs[6002]=['',4.2382557,-0.3212013,6.47];
ybs[6003]=['8 Her',4.231913,0.299193,6.14];
ybs[6004]=['',4.2341189,0.1102314,5.97];
ybs[6005]=['',4.2454072,-0.718754,5.86];
ybs[6006]=['',4.2373409,-0.0616044,5.37];
ybs[6007]=['',4.2436082,-0.5144954,5.13];
ybs[6008]=['τ CrB',4.2318574,0.6357806,4.76];
ybs[6009]=['ζ Nor',4.2557465,-0.9704293,5.81];
ybs[6010]=['δ1 Aps',4.2939178,-1.3744865,4.68];
ybs[6011]=['δ2 Aps',4.2943279,-1.3739862,5.27];
ybs[6012]=['',4.2551161,-0.9378073,5.83];
ybs[6013]=['φ Her',4.2304351,0.7831567,4.26];
ybs[6014]=['κ Nor',4.2560871,-0.9545411,4.94];
ybs[6015]=['',4.2167565,1.18238,5.44];
ybs[6016]=['ν Sco',4.2472788,-0.3405361,6.3];
ybs[6017]=['ν Sco',4.2473593,-0.340725,4.01];
ybs[6018]=['13 Sco',4.2490811,-0.4884783,4.59];
ybs[6019]=['12 Sco',4.2489372,-0.4970502,5.67];
ybs[6020]=['δ TrA',4.265969,-1.1125609,3.85];
ybs[6021]=['ψ Sco',4.2470073,-0.1767273,4.94];
ybs[6022]=['',4.2440637,0.1684351,6.53];
ybs[6023]=['16 Sco',4.2474806,-0.1502555,5.43];
ybs[6024]=['',4.2008044,1.3391394,5.56];
ybs[6025]=['',4.2437192,0.2897883,6.08];
ybs[6026]=['',4.2303885,1.0100992,6.33];
ybs[6027]=['',4.2741597,-1.1868227,5.75];
ybs[6028]=['',4.2323252,0.973296,6.49];
ybs[6029]=['10 Her',4.2441082,0.408981,5.7];
ybs[6030]=['',4.2667559,-1.0117947,5.63];
ybs[6031]=['',4.2508962,-0.0747334,6.25];
ybs[6032]=['',4.2552856,-0.4273009,6.41];
ybs[6033]=['',4.2437585,0.5808568,6.29];
ybs[6034]=['',4.2583666,-0.5772088,5.92];
ybs[6035]=['θ Nor',4.2631595,-0.8278435,5.14];
ybs[6036]=['',4.2441921,0.6346577,5.63];
ybs[6037]=['9 Her',4.2519263,0.0865714,5.48];
ybs[6038]=['χ Sco',4.2551423,-0.2076607,5.22];
ybs[6039]=['',4.2634559,-0.7497826,6.14];
ybs[6040]=['',4.2437935,0.7384947,5.87];
ybs[6041]=['',4.258276,-0.3694464,6.41];
ybs[6042]=['',4.2488564,0.4644248,6.5];
ybs[6043]=['',4.2589169,-0.3245511,6.32];
ybs[6044]=['',4.2602737,-0.4457036,6.05];
ybs[6045]=['',4.2701549,-0.9402089,5.44];
ybs[6046]=['δ Oph',4.2570001,-0.0655332,2.74];
ybs[6047]=['',4.256118,0.1019539,6.31];
ybs[6048]=['γ1 Nor',4.2710804,-0.8748828,4.99];
ybs[6049]=['',4.2728309,-0.9275591,6.33];
ybs[6050]=['18 Sco',4.2627478,-0.1471153,5.5];
ybs[6051]=['',4.2640334,-0.2602052,6.09];
ybs[6052]=['',4.2800258,-1.0115479,6.49];
ybs[6053]=['σ CrB',4.2568915,0.5898926,5.64];
ybs[6054]=['σ CrB',4.2568916,0.5898878,6.66];
ybs[6055]=['16 Her',4.2610787,0.3272239,5.69];
ybs[6056]=['',4.2691971,-0.3728506,6.61];
ybs[6057]=['',4.2682534,-0.0700276,6.18];
ybs[6058]=['',4.2620567,0.4775672,6.14];
ybs[6059]=['',4.2434245,1.1708104,6.21];
ybs[6060]=['',4.2752869,-0.5004215,4.78];
ybs[6061]=['λ Nor',4.2804512,-0.7458046,5.45];
ybs[6062]=['γ2 Nor',4.2834296,-0.8763785,4.02];
ybs[6063]=['',4.2864778,-0.9633673,5.77];
ybs[6064]=['υ CrB',4.2661453,0.5077363,5.78];
ybs[6065]=['ε Oph',4.2743905,-0.0829154,3.24];
ybs[6066]=['',4.2785279,-0.3538743,6.29];
ybs[6067]=['',4.2808443,-0.5404262,5.49];
ybs[6068]=['',4.2777752,-0.260583,5.94];
ybs[6069]=['19 UMi',4.2331195,1.3232175,5.48];
ybs[6070]=['',4.2856873,-0.689191,6.12];
ybs[6071]=['ο Sco',4.2852859,-0.4228295,4.55];
ybs[6072]=['20 UMi',4.2409498,1.3115933,6.39];
ybs[6073]=['',4.2948508,-0.8661726,5.33];
ybs[6074]=['σ Sco',4.2877595,-0.447666,2.89];
ybs[6075]=['',4.2944571,-0.7673873,5.88];
ybs[6076]=['',4.2658817,1.0418916,5.4];
ybs[6077]=['',4.2810147,0.3678313,6.05];
ybs[6078]=['',4.2506558,1.2799256,5.98];
ybs[6079]=['',4.3092691,-1.1026788,6.15];
ybs[6080]=['',4.2755064,0.854865,5.91];
ybs[6081]=['',4.2793881,0.6920431,5.46];
ybs[6082]=['τ Her',4.2781442,0.8073152,3.89];
ybs[6083]=['σ Ser',4.2905444,0.0169817,4.82];
ybs[6084]=['',4.3008615,-0.6850081,5.4];
ybs[6085]=['γ Her',4.2891614,0.3333011,3.75];
ybs[6086]=['',4.2931712,-0.037273,6.23];
ybs[6087]=['',4.3001757,-0.5804009,6.47];
ybs[6088]=['ζ TrA',4.3245627,-1.2241173,4.91];
ybs[6089]=['',4.3051416,-0.7924436,6.33];
ybs[6090]=['',4.302973,-0.656603,5.42];
ybs[6091]=['',4.2680291,1.1954774,6.41];
ybs[6092]=['γ Aps',4.3515615,-1.3778755,3.89];
ybs[6093]=['ξ CrB',4.289388,0.5381843,4.85];
ybs[6094]=['ψ Oph',4.3002354,-0.3506808,4.5];
ybs[6095]=['',4.3031143,-0.5193757,6.63];
ybs[6096]=['',4.3031215,-0.5193999,5.84];
ybs[6097]=['ν1 CrB',4.2903677,0.5889272,5.2];
ybs[6098]=['ν2 CrB',4.2909398,0.5872608,5.39];
ybs[6099]=['ι TrA',4.3207648,-1.1189436,5.27];
ybs[6100]=['',4.2929976,0.5633445,6.4];
ybs[6101]=['21 Her',4.2995101,0.1203055,5.85];
ybs[6102]=['ρ Oph',4.3068608,-0.4101774,5.02];
ybs[6103]=['ρ Oph',4.3068534,-0.410158,5.92];
ybs[6104]=['',4.3211498,-1.0236759,5.69];
ybs[6105]=['ε Nor',4.3152848,-0.8309206,4.47];
ybs[6106]=['η UMi',4.2622185,1.3211463,4.95];
ybs[6107]=['ω Her',4.3046258,0.2439785,4.57];
ybs[6108]=['χ Oph',4.3129129,-0.3230575,4.42];
ybs[6109]=['',4.3060747,0.3287903,6.7];
ybs[6110]=['',4.3275731,-1.008933,6.06];
ybs[6111]=['',4.3081119,0.1981565,6.11];
ybs[6112]=['',4.3191783,-0.6498237,5.79];
ybs[6113]=['25 Her',4.3034456,0.6516957,5.54];
ybs[6114]=['',4.3112703,0.0400362,6.07];
ybs[6115]=['',4.3328532,-1.076602,5.2];
ybs[6116]=['',4.2837758,1.2051981,5.25];
ybs[6117]=['',4.2976982,0.9625473,5.74];
ybs[6118]=['',4.3155359,-0.1335376,5.23];
ybs[6119]=['υ Oph',4.3159007,-0.1470388,4.63];
ybs[6120]=['',4.2940103,1.0758412,5.67];
ybs[6121]=['',4.3261832,-0.8080037,5.35];
ybs[6122]=['η Dra',4.2949482,1.072658,2.74];
ybs[6123]=['',4.5813732,-1.5273261,6.57];
ybs[6124]=['α Sco',4.3236749,-0.462234,0.96];
ybs[6125]=['',4.3504815,-1.2398314,5.5];
ybs[6126]=['',4.3188903,0.0106876,5.39];
ybs[6127]=['',4.3203133,-0.1427919,6.48];
ybs[6128]=['',4.4136933,-1.4535213,6.57];
ybs[6129]=['',4.4969103,-1.5048672,6.04];
ybs[6130]=['',4.3247887,-0.2548664,5.68];
ybs[6131]=['22 Sco',4.3271047,-0.4392413,4.79];
ybs[6132]=['',4.3345431,-0.7307299,5.33];
ybs[6133]=['',4.3327282,-0.6065968,4.23];
ybs[6134]=['',4.3276351,-0.1320614,6.5];
ybs[6135]=['',4.3322941,-0.4640621,6.1];
ybs[6136]=['30 Her',4.3172753,0.7300528,5.04];
ybs[6137]=['φ Oph',4.3307953,-0.290841,4.28];
ybs[6138]=['β Her',4.3252629,0.3741626,2.77];
ybs[6139]=['λ Oph',4.3290787,0.0337291,3.82];
ybs[6140]=['',4.31685,0.8963147,6.29];
ybs[6141]=['θ TrA',4.3551323,-1.1439523,5.52];
ybs[6142]=['',4.3267905,0.3565284,5.25];
ybs[6143]=['ω Oph',4.3353623,-0.3755427,4.45];
ybs[6144]=['',4.3296077,0.3864865,5.76];
ybs[6145]=['μ Nor',4.3451335,-0.7695983,4.94];
ybs[6146]=['34 Her',4.323079,0.8536209,6.45];
ybs[6147]=['',4.3281802,0.6138956,6.25];
ybs[6148]=['28 Her',4.3362789,0.0954814,5.63];
ybs[6149]=['29 Her',4.3360869,0.199624,4.84];
ybs[6150]=['',4.3498035,-0.7905224,6.46];
ybs[6151]=['15 Dra',4.3107374,1.1992985,5];
ybs[6152]=['',4.3307196,0.7949508,5.65];
ybs[6153]=['β Aps',4.3924262,-1.3537021,4.24];
ybs[6154]=['',4.3550534,-0.7488698,5.47];
ybs[6155]=['τ Sco',4.3520348,-0.4933108,2.82];
ybs[6156]=['',4.3545606,-0.6161671,4.16];
ybs[6157]=['',4.3679457,-1.0652957,6.18];
ybs[6158]=['σ Her',4.341048,0.7397967,4.2];
ybs[6159]=['',4.3482173,0.2968515,6.41];
ybs[6160]=['',4.3318113,1.0606814,5.94];
ybs[6161]=['12 Oph',4.3530037,-0.0414176,5.75];
ybs[6162]=['η1 TrA',4.3805101,-1.1927802,5.91];
ybs[6163]=['',4.2953413,1.3772218,5.56];
ybs[6164]=['',4.3640956,-0.7582697,5.83];
ybs[6165]=['ζ Oph',4.3568199,-0.1852681,2.56];
ybs[6166]=['',4.3538594,0.2696513,6.3];
ybs[6167]=['',4.3764355,-1.0557791,6.18];
ybs[6168]=['',4.3665173,-0.650383,5.91];
ybs[6169]=['',4.3604486,-0.1149376,6.09];
ybs[6170]=['',4.3245472,1.2664189,6.3];
ybs[6171]=['',4.3586605,0.2380521,6.31];
ybs[6172]=['',4.3889263,-1.1776881,6.03];
ybs[6173]=['',4.3498119,0.8127083,5.79];
ybs[6174]=['16 Dra',4.3492467,0.9224356,5.53];
ybs[6175]=['17 Dra',4.3494039,0.9228577,5.08];
ybs[6176]=['17 Dra',4.349433,0.922853,6.53];
ybs[6177]=['',4.3772048,-0.851868,5.65];
ybs[6178]=['',4.378735,-0.8673739,5.65];
ybs[6179]=['',4.3676599,-0.1675676,6.35];
ybs[6180]=['',4.3721461,-0.356999,6.26];
ybs[6181]=['',4.3181318,1.3507897,6.34];
ybs[6182]=['',4.3779423,-0.5793028,5.87];
ybs[6183]=['',4.3768208,-0.4278393,6.09];
ybs[6184]=['36 Her',4.3711908,0.0726274,6.93];
ybs[6185]=['37 Her',4.3714521,0.072846,5.77];
ybs[6186]=['',4.3763864,-0.3104521,4.96];
ybs[6187]=['',4.3844777,-0.8048584,6.23];
ybs[6188]=['',4.35097,1.0999846,6.16];
ybs[6189]=['',4.3567984,0.9768247,5.29];
ybs[6190]=['42 Her',4.3607797,0.8531373,4.9];
ybs[6191]=['',4.3740419,-0.0182594,6.24];
ybs[6192]=['',4.3778858,-0.3485358,5.57];
ybs[6193]=['',4.3720482,0.2155333,6.08];
ybs[6194]=['',4.4032582,-1.1720221,5.13];
ybs[6195]=['14 Oph',4.3762107,0.0198229,5.74];
ybs[6196]=['',4.3871692,-0.718429,6.2];
ybs[6197]=['',4.3921631,-0.9284451,5.96];
ybs[6198]=['',4.3721713,0.4330652,6.06];
ybs[6199]=['',4.38778,-0.7183259,6.12];
ybs[6200]=['',4.3871253,-0.6667284,6.05];
ybs[6201]=['',4.3861209,-0.5611281,6.46];
ybs[6202]=['ζ Her',4.3730374,0.5507804,2.81];
ybs[6203]=['39 Her',4.3746961,0.4689957,5.92];
ybs[6204]=['',4.3913051,-0.7135473,5.71];
ybs[6205]=['',4.4002189,-1.0218221,5.74];
ybs[6206]=['',4.388686,-0.4799642,6.58];
ybs[6207]=['α TrA',4.4125396,-1.2054781,1.92];
ybs[6208]=['',4.3918603,-0.4983461,6.02];
ybs[6209]=['',4.4044403,-1.0189812,5.58];
ybs[6210]=['η Her',4.3796269,0.6785391,3.53];
ybs[6211]=['',4.4003169,-0.6880014,5.48];
ybs[6212]=['',4.3841319,0.5933188,5.99];
ybs[6213]=['18 Dra',4.3680977,1.1264892,4.83];
ybs[6214]=['16 Oph',4.3927393,0.017053,6.03];
ybs[6215]=['25 Sco',4.3997944,-0.4462978,6.71];
ybs[6216]=['',4.3784839,0.971196,6.16];
ybs[6217]=['',4.3916195,0.2740508,5.56];
ybs[6218]=['43 Her',4.3939076,0.1490417,5.15];
ybs[6219]=['η Ara',4.4152894,-1.0311741,3.76];
ybs[6220]=['',4.3893513,0.7535233,6.05];
ybs[6221]=['',4.4279901,-1.1819534,6.32];
ybs[6222]=['19 Oph',4.3999774,0.0352936,6.1];
ybs[6223]=['',4.425697,-1.1417038,6.13];
ybs[6224]=['45 Her',4.4025126,0.0908397,5.24];
ybs[6225]=['',4.4062769,-0.2609432,6.03];
ybs[6226]=['',4.4177639,-0.8741603,6.47];
ybs[6227]=['',4.3884709,0.9902719,4.85];
ybs[6228]=['',4.3482144,1.3765434,6.32];
ybs[6229]=['',4.4038089,0.2364713,6.35];
ybs[6230]=['',4.4107315,-0.2741637,6.1];
ybs[6231]=['ε Sco',4.4146996,-0.5992378,2.29];
ybs[6232]=['',4.3987319,0.7364693,5.87];
ybs[6233]=['20 Oph',4.4121417,-0.1889106,4.65];
ybs[6234]=['',4.4185328,-0.6554484,6.11];
ybs[6235]=['',4.4212471,-0.7203008,5.22];
ybs[6236]=['',4.4100539,0.2307405,5.91];
ybs[6237]=['μ1 Sco',4.4223794,-0.6647431,3.08];
ybs[6238]=['',4.4141633,-0.0470249,6.32];
ybs[6239]=['',4.424595,-0.7311822,6.49];
ybs[6240]=['47 Her',4.4135579,0.1257908,5.49];
ybs[6241]=['',4.4335727,-1.0113761,5.94];
ybs[6242]=['μ2 Sco',4.4244075,-0.6642149,3.57];
ybs[6243]=['',4.440599,-1.1049156,6.02];
ybs[6244]=['52 Her',4.406758,0.8018407,4.82];
ybs[6245]=['21 Oph',4.418562,0.0205295,5.51];
ybs[6246]=['',4.4088694,0.7572909,6.13];
ybs[6247]=['',4.4307514,-0.7520491,5.96];
ybs[6248]=['50 Her',4.41397,0.5195198,5.72];
ybs[6249]=['',4.4141191,0.5674635,6.13];
ybs[6250]=['',4.4320646,-0.7303263,5.45];
ybs[6251]=['',4.4318599,-0.7336139,6.32];
ybs[6252]=['ζ1 Sco',4.431952,-0.7400278,4.73];
ybs[6253]=['',4.4327954,-0.7310907,6.45];
ybs[6254]=['',4.4130224,0.7305288,6.29];
ybs[6255]=['',4.4333608,-0.730561,6.59];
ybs[6256]=['',4.4339399,-0.7420595,5.88];
ybs[6257]=['',4.3722044,1.3520883,5.98];
ybs[6258]=['49 Her',4.4209098,0.2606588,6.52];
ybs[6259]=['',4.4282085,-0.3569934,5.88];
ybs[6260]=['51 Her',4.4190486,0.4296418,5.04];
ybs[6261]=['ζ2 Sco',4.4345206,-0.7400074,3.62];
ybs[6262]=['',4.4361364,-0.7188803,5.77];
ybs[6263]=['',4.4338592,-0.53451,6.35];
ybs[6264]=['',4.4420615,-0.8850908,6.33];
ybs[6265]=['',4.4436688,-0.9131677,5.94];
ybs[6266]=['',4.4602504,-1.2095677,5.79];
ybs[6267]=['',4.4307246,-0.0288066,6.25];
ybs[6268]=['',4.4332976,-0.2064806,6.57];
ybs[6269]=['53 Her',4.4239795,0.5526168,5.32];
ybs[6270]=['23 Oph',4.4327321,-0.1080693,5.25];
ybs[6271]=['ι Oph',4.4295162,0.1767473,4.38];
ybs[6272]=['',4.4399729,-0.5854599,6.37];
ybs[6273]=['',4.4432151,-0.713148,6.15];
ybs[6274]=['',4.4394388,-0.2939708,6.37];
ybs[6275]=['ζ Ara',4.4534559,-0.9778345,3.13];
ybs[6276]=['',4.4243244,0.8269021,6];
ybs[6277]=['',4.4330279,0.3651351,5.41];
ybs[6278]=['27 Sco',4.4452914,-0.5811231,5.48];
ybs[6279]=['',4.4514538,-0.8844778,5.55];
ybs[6280]=['',4.434862,0.2370511,6.34];
ybs[6281]=['24 Oph',4.4430872,-0.4046845,5.58];
ybs[6282]=['56 Her',4.4333035,0.4484221,6.08];
ybs[6283]=['54 Her',4.4351046,0.3210653,5.35];
ybs[6284]=['',4.4440843,-0.3416756,6.27];
ybs[6285]=['ε1 Ara',4.4573312,-0.9284372,4.06];
ybs[6286]=['',4.4453219,-0.1919813,6.19];
ybs[6287]=['',4.459761,-0.9535015,5.65];
ybs[6288]=['',4.4529253,-0.6572314,6.09];
ybs[6289]=['κ Oph',4.4455167,0.1629909,3.2];
ybs[6290]=['',4.4606578,-0.8496662,6];
ybs[6291]=['',4.4447394,0.2416892,6.37];
ybs[6292]=['',4.4509781,-0.2601478,6.59];
ybs[6293]=['',4.4607588,-0.7938829,6.65];
ybs[6294]=['',4.4677146,-1.0296044,6.11];
ybs[6295]=['57 Her',4.4441469,0.4418534,6.28];
ybs[6296]=['',4.4363502,0.8726914,6.56];
ybs[6297]=['',4.4450157,0.4249015,6.32];
ybs[6298]=['',4.4569726,-0.438546,5.86];
ybs[6299]=['26 Oph',4.4578332,-0.4367502,5.75];
ybs[6300]=['',4.4604225,-0.6277712,5.97];
ybs[6301]=['',4.4666479,-0.8929902,6.45];
ybs[6302]=['',4.4445271,0.7413492,6.34];
ybs[6303]=['ε2 Ara',4.4728926,-0.9297346,5.29];
ybs[6304]=['19 Dra',4.4338049,1.1361601,4.89];
ybs[6305]=['',4.4657147,-0.5616011,5.03];
ybs[6306]=['',4.4579416,0.1143005,6.59];
ybs[6307]=['30 Oph',4.4608645,-0.0742955,4.82];
ybs[6308]=['20 Dra',4.435533,1.1344963,6.41];
ybs[6309]=['',4.4789817,-1.0078295,5.73];
ybs[6310]=['29 Oph',4.464945,-0.3302054,6.26];
ybs[6311]=['ε UMi',4.3787788,1.431052,4.23];
ybs[6312]=['',4.4746851,-0.8236671,6.06];
ybs[6313]=['ε Her',4.4559561,0.5391588,3.92];
ybs[6314]=['',4.4593342,0.3944056,5.65];
ybs[6315]=['',4.462217,0.260322,6.31];
ybs[6316]=['',4.4746698,-0.6664509,5.91];
ybs[6317]=['',4.4599471,0.474067,6.55];
ybs[6318]=['',4.4643879,0.1468997,6.33];
ybs[6319]=['',4.4497738,0.9887826,6.03];
ybs[6320]=['',4.4806276,-0.7947098,6.28];
ybs[6321]=['59 Her',4.4615505,0.5852824,5.25];
ybs[6322]=['',4.4650554,0.444568,5.75];
ybs[6323]=['',4.4786834,-0.5961143,4.87];
ybs[6324]=['',4.4323005,1.2756648,6.3];
ybs[6325]=['',4.4646121,0.5559049,6.36];
ybs[6326]=['',4.4691469,0.2453718,4.98];
ybs[6327]=['',4.483666,-0.7703262,6.19];
ybs[6328]=['',4.469318,0.2526881,6.52];
ybs[6329]=['',4.4776817,-0.3582614,6.3];
ybs[6330]=['',4.4714581,0.2368832,5.93];
ybs[6331]=['',4.4728194,0.2362269,6.08];
ybs[6332]=['',4.4721648,0.3430932,6.35];
ybs[6333]=['',4.4854867,-0.6502863,5.98];
ybs[6334]=['',4.4458582,1.2069039,6.4];
ybs[6335]=['61 Her',4.4697135,0.6175173,6.69];
ybs[6336]=['',4.4859651,-0.6192813,6.13];
ybs[6337]=['',4.4575534,1.0579259,6.13];
ybs[6338]=['',4.4790865,0.0117042,6.01];
ybs[6339]=['',4.4840232,-0.3769217,6.3];
ybs[6340]=['',4.4714448,0.6066324,6.04];
ybs[6341]=['',4.4757025,0.3415064,6.17];
ybs[6342]=['',4.4802708,-0.0161214,5.64];
ybs[6343]=['',4.4872749,-0.4632792,6.29];
ybs[6344]=['60 Her',4.4790133,0.2218132,4.91];
ybs[6345]=['',4.5044797,-1.0769443,6.39];
ybs[6346]=['',4.5165412,-1.2347936,6.22];
ybs[6347]=['',4.4825572,0.1693305,6.37];
ybs[6348]=['',4.4827744,0.1819119,6.37];
ybs[6349]=['',4.4610559,1.1268996,6.1];
ybs[6350]=['',4.4861632,-0.0294498,6.38];
ybs[6351]=['',4.4760306,0.7641063,6.43];
ybs[6352]=['',4.4745027,0.8512291,6.09];
ybs[6353]=['',4.4826256,0.384894,5.56];
ybs[6354]=['',4.4927963,-0.3078637,5.99];
ybs[6355]=['',4.4957935,-0.5311627,5.97];
ybs[6356]=['',4.4920086,-0.0193667,6.06];
ybs[6357]=['',4.5196194,-1.1732721,5.89];
ybs[6358]=['μ Dra',4.4760868,0.9501258,5.83];
ybs[6359]=['μ Dra',4.4760795,0.9501258,5.8];
ybs[6360]=['',4.5050712,-0.7781742,5.08];
ybs[6361]=['',4.4950968,-0.0682871,6.36];
ybs[6362]=['',4.5369385,-1.3012791,6.25];
ybs[6363]=['',4.5095636,-0.8535041,5.84];
ybs[6364]=['',4.4992618,-0.1841772,5.56];
ybs[6365]=['',4.4880265,0.7066058,6.34];
ybs[6366]=['',4.4894395,0.6266581,5.39];
ybs[6367]=['η Oph',4.5020083,-0.2749525,2.43];
ybs[6368]=['',4.4545897,1.3135812,6.21];
ybs[6369]=['η Sco',4.5112903,-0.7551503,3.33];
ybs[6370]=['',4.5115295,-0.6900099,5.67];
ybs[6371]=['',4.511505,-0.6780592,6.3];
ybs[6372]=['',4.4893702,0.8868346,6.46];
ybs[6373]=['',4.5216977,-0.9933505,6.09];
ybs[6374]=['',4.5025186,0.2170924,6.57];
ybs[6375]=['',4.5105177,-0.4412631,6.54];
ybs[6376]=['',4.5114743,-0.4850202,6.14];
ybs[6377]=['',4.4957567,0.7111811,5.08];
ybs[6378]=['',4.5141805,-0.5666375,6.01];
ybs[6379]=['',4.5070298,0.1372974,6.33];
ybs[6380]=['63 Her',4.5032507,0.4225296,6.19];
ybs[6381]=['',4.521114,-0.6945253,6.6];
ybs[6382]=['37 Oph',4.5100151,0.1842635,5.33];
ybs[6383]=['',4.5123708,0.0056632,6.65];
ybs[6384]=['',4.4989075,0.9142005,6.29];
ybs[6385]=['ζ Dra',4.4892663,1.1464111,3.17];
ybs[6386]=['',4.5244725,-0.5859818,5.53];
ybs[6387]=['',4.5259959,-0.6740401,5.96];
ybs[6388]=['',4.5042398,0.867748,6.04];
ybs[6389]=['',4.5507349,-1.2229228,6.53];
ybs[6390]=['36 Oph',4.5242168,-0.4647592,5.11];
ybs[6391]=['36 Oph',4.5242022,-0.4647349,5.07];
ybs[6392]=['',4.5266257,-0.5277214,6.21];
ybs[6393]=['',4.5236033,-0.2549908,5.99];
ybs[6394]=['',4.52913,-0.6243877,6.12];
ybs[6395]=['α1 Her',4.5193854,0.2506954,3.48];
ybs[6396]=['α2 Her',4.5194072,0.2506906,5.39];
ybs[6397]=['',4.5438874,-1.042275,5.91];
ybs[6398]=['',4.5320153,-0.5705082,5.55];
ybs[6399]=['δ Her',4.520576,0.4330662,3.14];
ybs[6400]=['ι Aps',4.5589959,-1.2242608,5.41];
ybs[6401]=['',4.5268257,0.0377093,6.17];
ybs[6402]=['',4.5292509,-0.1094361,6.09];
ybs[6403]=['',4.5281382,0.0206857,5.88];
ybs[6404]=['41 Oph',4.5285694,-0.0082132,4.73];
ybs[6405]=['',4.5416875,-0.8143286,5.48];
ybs[6406]=['ζ Aps',4.5577715,-1.1831991,4.78];
ybs[6407]=['π Her',4.519975,0.6419816,3.16];
ybs[6408]=['',4.5235163,0.4139374,5.96];
ybs[6409]=['',4.5403738,-0.7706256,5.76];
ybs[6410]=['',4.5062538,1.0968778,5.56];
ybs[6411]=['',4.5375883,-0.568585,6.36];
ybs[6412]=['',4.5439776,-0.8741787,6.27];
ybs[6413]=['ο Oph',4.5357168,-0.4243132,5.2];
ybs[6414]=['ο Oph',4.5357022,-0.4242648,6.8];
ybs[6415]=['',4.5404095,-0.6111016,5.91];
ybs[6416]=['',4.5430366,-0.7722482,6.65];
ybs[6417]=['',4.5366908,-0.2851204,6.43];
ybs[6418]=['',4.608028,-1.4115306,5.88];
ybs[6419]=['',4.5318584,0.4025784,6.45];
ybs[6420]=['68 Her',4.5301383,0.5772679,4.82];
ybs[6421]=['',4.5342424,0.3018297,6];
ybs[6422]=['',4.53685,0.1891987,5.03];
ybs[6423]=['',4.5381988,0.1057894,6.51];
ybs[6424]=['',4.5435994,-0.3103146,6.02];
ybs[6425]=['69 Her',4.5313939,0.6504295,4.65];
ybs[6426]=['',4.5266601,0.8668307,7.48];
ybs[6427]=['',4.5598273,-1.012843,5.88];
ybs[6428]=['',4.5435423,-0.1036818,6.32];
ybs[6429]=['',4.5654617,-1.0975471,5.7];
ybs[6430]=['',4.546646,-0.3378208,6.52];
ybs[6431]=['',4.5604844,-0.9869231,5.8];
ybs[6432]=['',4.5368467,0.5026365,5.65];
ybs[6433]=['',4.5344197,0.6769609,5.94];
ybs[6434]=['ξ Oph',4.5486199,-0.3688831,4.39];
ybs[6435]=['ν Ser',4.5474851,-0.2246194,4.33];
ybs[6436]=['',4.5663395,-1.0593121,5.77];
ybs[6437]=['',4.5238606,1.0584533,6.32];
ybs[6438]=['',4.5476114,-0.1870798,6.46];
ybs[6439]=['',4.5567513,-0.6601993,6.41];
ybs[6440]=['ι Ara',4.5601602,-0.8288492,5.25];
ybs[6441]=['',4.5439452,0.3147534,5];
ybs[6442]=['θ Oph',4.5532013,-0.4367077,3.27];
ybs[6443]=['',4.5565278,-0.6271257,6.47];
ybs[6444]=['',4.5429276,0.4453069,5.38];
ybs[6445]=['',4.5578406,-0.6499962,5.93];
ybs[6446]=['70 Her',4.5462077,0.427197,5.12];
ybs[6447]=['72 Her',4.5447225,0.5662676,5.39];
ybs[6448]=['43 Oph',4.5592542,-0.4915601,5.35];
ybs[6449]=['',4.5640198,-0.771142,5.12];
ybs[6450]=['β Ara',4.5698923,-0.9695299,2.85];
ybs[6451]=['γ Ara',4.5704109,-0.9843205,3.34];
ybs[6452]=['',4.549422,0.2916163,6.35];
ybs[6453]=['74 Her',4.542421,0.8066486,5.59];
ybs[6454]=['',4.5558908,-0.0420619,6.29];
ybs[6455]=['',4.5486794,0.5015299,6.35];
ybs[6456]=['',4.5431568,0.8406408,6.43];
ybs[6457]=['κ Ara',4.5723876,-0.8840652,5.23];
ybs[6458]=['',4.5489014,0.6972941,5.51];
ybs[6459]=['',4.5669786,-0.6059195,6.16];
ybs[6460]=['',4.583441,-1.1005106,6.24];
ybs[6461]=['',4.5647675,-0.3745854,5.85];
ybs[6462]=['',4.5642718,-0.3222992,6.21];
ybs[6463]=['',4.5666632,-0.4234842,6.19];
ybs[6464]=['',4.5766219,-0.9070162,6.19];
ybs[6465]=['',4.5602479,0.1541433,5.77];
ybs[6466]=['',4.5756886,-0.8004557,5.29];
ybs[6467]=['',4.5776544,-0.8839946,5.92];
ybs[6468]=['',4.5478374,0.9319719,5.67];
ybs[6469]=['73 Her',4.5602723,0.4003664,5.74];
ybs[6470]=['',4.5623915,0.2841467,5.71];
ybs[6471]=['',4.5625896,0.272017,6.35];
ybs[6472]=['',4.5811197,-0.9130804,5.75];
ybs[6473]=['ρ Her',4.5576041,0.6479602,5.47];
ybs[6474]=['ρ Her',4.557626,0.6479458,4.52];
ybs[6475]=['44 Oph',4.5721876,-0.4222786,4.17];
ybs[6476]=['',4.5844486,-0.9632075,5.94];
ybs[6477]=['',4.5590726,0.6730287,6.49];
ybs[6478]=['',4.5694307,-0.0291727,6.44];
ybs[6479]=['',4.5746765,-0.4531312,6.44];
ybs[6480]=['',4.56099,0.6445697,6.28];
ybs[6481]=['45 Oph',4.5767817,-0.5216064,4.29];
ybs[6482]=['',4.5724814,-0.0891178,4.54];
ybs[6483]=['',4.5779523,-0.5191165,6];
ybs[6484]=['',4.568392,0.2949189,5.98];
ybs[6485]=['',4.5745545,-0.2187185,6.21];
ybs[6486]=['',4.570583,0.132225,6.06];
ybs[6487]=['σ Oph',4.571592,0.0719212,4.34];
ybs[6488]=['',4.5683872,0.4687784,6.41];
ybs[6489]=['δ Ara',4.5959704,-1.0594217,3.62];
ybs[6490]=['',4.5840859,-0.6422161,6.02];
ybs[6491]=['',4.5722246,0.3501387,5.54];
ybs[6492]=['',4.5863533,-0.6725506,6.39];
ybs[6493]=['',4.5787495,-0.1435865,6.37];
ybs[6494]=['',4.5966336,-0.9937368,5.95];
ybs[6495]=['',4.5712441,0.6052171,5.94];
ybs[6496]=['',4.5818427,0.005453,5.44];
ybs[6497]=['υ Sco',4.5920989,-0.6512294,2.69];
ybs[6498]=['77 Her',4.5700896,0.8419549,5.85];
ybs[6499]=['α Ara',4.5978097,-0.8707844,2.95];
ybs[6500]=['',4.5640712,1.0476877,5.65];
ybs[6501]=['',4.5862972,-0.1036247,6.37];
ybs[6502]=['',4.5973539,-0.8037692,6.03];
ybs[6503]=['',4.5660071,1.0233204,6.51];
ybs[6504]=['',4.5887433,-0.0188442,5.31];
ybs[6505]=['',4.5963556,-0.5885083,6.44];
ybs[6506]=['',4.5595426,1.1743558,6.43];
ybs[6507]=['51 Oph',4.5942004,-0.4185178,4.81];
ybs[6508]=['',4.5957275,-0.458778,6.05];
ybs[6509]=['',4.588099,0.2078295,6.39];
ybs[6510]=['',4.5990887,-0.5985715,6.17];
ybs[6511]=['',4.6026657,-0.7188847,5.84];
ybs[6512]=['',4.5927685,0.0472601,5.59];
ybs[6513]=['',4.6154007,-1.0447535,6.28];
ybs[6514]=['λ Her',4.5890328,0.4554172,4.41];
ybs[6515]=['λ Sco',4.6045026,-0.6478499,1.63];
ybs[6516]=['',4.5895694,0.5435191,5.61];
ybs[6517]=['',4.5282096,1.3982171,5.72];
ybs[6518]=['',4.613423,-0.9314271,6.1];
ybs[6519]=['',4.5879952,0.6783232,6.43];
ybs[6520]=['',4.5962867,0.2079409,6.42];
ybs[6521]=['78 Her',4.5936609,0.4955175,5.62];
ybs[6522]=['',4.6024709,-0.1005322,5.62];
ybs[6523]=['',4.6090254,-0.5689113,5.7];
ybs[6524]=['β Dra',4.5858049,0.912528,2.79];
ybs[6525]=['σ Ara',4.6141613,-0.811918,4.59];
ybs[6526]=['',4.594165,0.5978535,6.56];
ybs[6527]=['',4.6137257,-0.6536944,6.48];
ybs[6528]=['',4.5863938,1.0098323,6.4];
ybs[6529]=['',4.6008904,0.3358218,5.64];
ybs[6530]=['',4.6022314,0.2845267,5.69];
ybs[6531]=['',4.6025443,0.2587692,6.48];
ybs[6532]=['',4.6082667,-0.1964634,5.55];
ybs[6533]=['52 Oph',4.6110894,-0.3849868,6.57];
ybs[6534]=['',4.6174275,-0.6745473,4.29];
ybs[6535]=['',4.6223208,-0.8739358,5.93];
ybs[6536]=['53 Oph',4.6066911,0.1670617,5.81];
ybs[6537]=['π Ara',4.6256065,-0.9514261,5.25];
ybs[6538]=['',4.5985164,0.7195623,5.74];
ybs[6539]=['',4.6056988,0.287788,6.4];
ybs[6540]=['',4.7531316,-1.4820369,6.45];
ybs[6541]=['θ Sco',4.6211181,-0.7506789,1.87];
ybs[6542]=['ν1 Dra',4.5930684,0.9628596,4.88];
ybs[6543]=['ν2 Dra',4.5934624,0.9626666,4.87];
ybs[6544]=['α Oph',4.6079743,0.2189594,2.08];
ybs[6545]=['',4.6213169,-0.6645945,6.26];
ybs[6546]=['',4.6246949,-0.7486245,6.1];
ybs[6547]=['',4.6122106,0.3662078,6.1];
ybs[6548]=['',4.598644,1.0043198,6.17];
ybs[6549]=['ξ Ser',4.6207244,-0.2689817,3.54];
ybs[6550]=['',4.6208048,-0.2719922,5.94];
ybs[6551]=['',4.610065,0.6507894,6.1];
ybs[6552]=['',4.6124541,0.4916737,6.38];
ybs[6553]=['',4.656899,-1.2606386,6.49];
ybs[6554]=['27 Dra',4.5896568,1.1888891,5.05];
ybs[6555]=['μ Oph',4.621537,-0.1419242,4.62];
ybs[6556]=['',4.6230236,-0.190921,5.75];
ybs[6557]=['λ Ara',4.6350882,-0.8626577,4.77];
ybs[6558]=['',4.6144063,0.5370668,6.02];
ybs[6559]=['79 Her',4.6187031,0.4240616,5.77];
ybs[6560]=['',4.6386889,-0.8191272,5.79];
ybs[6561]=['26 Dra',4.6043224,1.0796645,5.23];
ybs[6562]=['82 Her',4.6131981,0.8477388,5.37];
ybs[6563]=['',4.6267726,0.0351862,6.26];
ybs[6564]=['',4.6424934,-0.8817517,6.24];
ybs[6565]=['',4.6255167,0.2324254,6.12];
ybs[6566]=['',4.631549,-0.037767,6.19];
ybs[6567]=['',4.6239652,0.5711962,6.37];
ybs[6568]=['κ Sco',4.643384,-0.6813749,2.41];
ybs[6569]=['ο Ser',4.6373183,-0.2249017,4.26];
ybs[6570]=['η Pav',4.6606522,-1.1297806,3.62];
ybs[6571]=['',4.6448232,-0.6449958,5.54];
ybs[6572]=['',4.6289834,0.5443833,6.03];
ybs[6573]=['μ Ara',4.6517105,-0.9048314,5.15];
ybs[6574]=['',4.6558496,-1.0045001,6.01];
ybs[6575]=['',4.6457388,-0.5770176,6.4];
ybs[6576]=['ι Her',4.6258115,0.8027536,3.8];
ybs[6577]=['',4.6351213,0.2647175,6.34];
ybs[6578]=['',4.6370515,0.1099934,5.95];
ybs[6579]=['',4.6321562,0.5458743,6.28];
ybs[6580]=['',4.634283,0.4276474,6.36];
ybs[6581]=['',4.6462428,-0.4868359,6.36];
ybs[6582]=['',4.6385554,0.2782327,5.52];
ybs[6583]=['58 Oph',4.6465122,-0.3786099,4.87];
ybs[6584]=['ω Dra',4.6112098,1.1998148,4.8];
ybs[6585]=['',4.6533041,-0.7459095,5.87];
ybs[6586]=['',4.6096798,1.2139971,6.42];
ybs[6587]=['',4.6310873,0.758512,6.59];
ybs[6588]=['',4.6477966,-0.2359305,6.39];
ybs[6589]=['',4.6474046,-0.123721,6.3];
ybs[6590]=['83 Her',4.6403087,0.4285489,5.52];
ybs[6591]=['β Oph',4.6455673,0.0795479,2.77];
ybs[6592]=['',4.6446866,0.249328,6.24];
ybs[6593]=['',4.6295383,1.0000539,6.77];
ybs[6594]=['',4.6106628,1.2643535,5.86];
ybs[6595]=['',4.6335334,0.9042053,5.99];
ybs[6596]=['84 Her',4.644189,0.4244323,5.71];
ybs[6597]=['61 Oph',4.6504206,0.044866,6.17];
ybs[6598]=['',4.6505224,0.0448565,6.56];
ybs[6599]=['',4.6487028,0.2513495,6.19];
ybs[6600]=['',4.6417994,0.769247,6.34];
ybs[6601]=['',4.6635408,-0.6652991,6.43];
ybs[6602]=['',4.6717268,-0.9670489,6.11];
ybs[6603]=['ι1 Sco',4.6657016,-0.7004673,3.03];
ybs[6604]=['3 Sgr',4.6648436,-0.4858608,4.54];
ybs[6605]=['',4.6654449,-0.3924355,6.18];
ybs[6606]=['',4.6447593,0.9388527,5.75];
ybs[6607]=['',4.6538926,0.5497171,6.23];
ybs[6608]=['',4.6644461,-0.2571356,5.94];
ybs[6609]=['',4.6687403,-0.4709144,6.35];
ybs[6610]=['',4.679514,-0.9357984,5.92];
ybs[6611]=['μ Her',4.6575326,0.4836794,3.42];
ybs[6612]=['',4.6853287,-1.0501435,5.78];
ybs[6613]=['',4.6544074,0.6784661,6.52];
ybs[6614]=['',4.6547261,0.6861657,6.68];
ybs[6615]=['',4.6609683,0.3087468,5.72];
ybs[6616]=['',4.6721043,-0.5534317,4.83];
ybs[6617]=['γ Oph',4.6649304,0.0471306,3.75];
ybs[6618]=['',4.6754121,-0.6466248,3.21];
ybs[6619]=['ι2 Sco',4.6770519,-0.6998052,4.81];
ybs[6620]=['',4.6825706,-0.9273844,6.09];
ybs[6621]=['',4.6668192,0.0662804,6.22];
ybs[6622]=['',4.6938828,-1.1430584,6.49];
ybs[6623]=['',4.7174935,-1.3295541,6.07];
ybs[6624]=['ψ1 Dra',4.631693,1.259046,4.58];
ybs[6625]=['ψ1 Dra',4.631813,1.2591869,5.79];
ybs[6626]=['',4.6664197,0.3588219,5.69];
ybs[6627]=['',4.6711646,0.0341233,6.47];
ybs[6628]=['',4.6842642,-0.7959562,6.11];
ybs[6629]=['',4.6591299,0.8308599,6.43];
ybs[6630]=['',4.668161,0.3359572,6.12];
ybs[6631]=['',4.6830353,-0.7116932,5.96];
ybs[6632]=['87 Her',4.6679457,0.4470909,5.12];
ybs[6633]=['',4.6809119,-0.5334072,6.66];
ybs[6634]=['',4.7575156,-1.4221238,6.35];
ybs[6635]=['',4.6856043,-0.6074322,5.9];
ybs[6636]=['',4.6860245,-0.6007554,5.84];
ybs[6637]=['',4.6889364,-0.733045,6.2];
ybs[6638]=['',4.676895,0.2084183,6.17];
ybs[6639]=['',4.6881515,-0.5954706,6.06];
ybs[6640]=['',4.6886942,-0.611255,6.45];
ybs[6641]=['',4.6888774,-0.6218235,6.03];
ybs[6642]=['',4.6745642,0.5116743,5.5];
ybs[6643]=['',4.6767782,0.3894044,5.98];
ybs[6644]=['30 Dra',4.6672707,0.8861872,5.02];
ybs[6645]=['',4.6904,-0.6062233,6.17];
ybs[6646]=['',4.6906791,-0.6090977,5.6];
ybs[6647]=['',4.6829765,-0.0216605,6.35];
ybs[6648]=['',4.6922869,-0.6071836,6.38];
ybs[6649]=['',4.6860392,-0.1072957,6.21];
ybs[6650]=['',4.6929685,-0.6066002,5.96];
ybs[6651]=['',4.6932061,-0.6079815,6.42];
ybs[6652]=['88 Her',4.6717966,0.8445376,6.68];
ybs[6653]=['',4.6822004,0.2674087,6.46];
ybs[6654]=['',4.6880449,-0.190301,6.18];
ybs[6655]=['',4.6854817,0.022706,5.95];
ybs[6656]=['',4.6952999,-0.6016011,5.96];
ybs[6657]=['',4.6776151,0.6993102,6.46];
ybs[6658]=['',4.6881031,0.1064253,5.77];
ybs[6659]=['',4.6983857,-0.6366654,6.06];
ybs[6660]=['',4.6967193,-0.4344091,6.2];
ybs[6661]=['',4.6813082,0.6977433,6.04];
ybs[6662]=['',4.6805145,0.8140003,6.38];
ybs[6663]=['',4.7061882,-0.7739419,4.86];
ybs[6664]=['',4.6922525,0.1942109,6.38];
ybs[6665]=['90 Her',4.6866438,0.6982069,5.16];
ybs[6666]=['',4.7064919,-0.7034877,6.43];
ybs[6667]=['',4.7008616,-0.3281956,6.52];
ybs[6668]=['',4.704731,-0.4898578,5.8];
ybs[6669]=['',4.7024706,-0.2760111,5.89];
ybs[6670]=['',4.7103691,-0.7280974,4.88];
ybs[6671]=['',4.7109185,-0.6830808,6.29];
ybs[6672]=['',4.701725,0.0116667,5.82];
ybs[6673]=['89 Her',4.6967206,0.4546159,5.46];
ybs[6674]=['',4.7040534,-0.07127,5.47];
ybs[6675]=['',4.6987615,0.3920359,5.58];
ybs[6676]=['ξ Dra',4.6859836,0.9925522,3.75];
ybs[6677]=['',4.705088,0.0011348,5.97];
ybs[6678]=['',4.7042165,0.1132074,6.29];
ybs[6679]=['',4.7149694,-0.6433018,5.74];
ybs[6680]=['',4.7132914,-0.501948,6.01];
ybs[6681]=['',4.7152713,-0.5280167,5.16];
ybs[6682]=['',4.7152931,-0.5280216,7.04];
ybs[6683]=['θ Her',4.6997101,0.6501103,3.86];
ybs[6684]=['',4.7062696,0.1927366,6.36];
ybs[6685]=['',4.7047583,0.4187829,6.3];
ybs[6686]=['ν Oph',4.7140311,-0.1705848,3.34];
ybs[6687]=['',4.6942502,0.9768395,6.1];
ybs[6688]=['4 Sgr',4.718017,-0.415664,4.76];
ybs[6689]=['35 Dra',4.6615687,1.3431387,5.04];
ybs[6690]=['',4.7015248,0.7914919,6.02];
ybs[6691]=['ξ Her',4.706783,0.5104517,3.7];
ybs[6692]=['',4.7187578,-0.3549778,6.21];
ybs[6693]=['γ Dra',4.7000706,0.8986183,2.23];
ybs[6694]=['',4.7163744,-0.0841464,5.87];
ybs[6695]=['ν Her',4.7099526,0.5268945,4.41];
ybs[6696]=['',4.7274974,-0.6348846,6.3];
ybs[6697]=['',4.7189691,0.0109951,6.37];
ybs[6698]=['ζ Ser',4.7201305,-0.0643958,4.62];
ybs[6699]=['',4.7104696,0.6333322,6];
ybs[6700]=['66 Oph',4.7188352,0.0762556,4.64];
ybs[6701]=['93 Her',4.7173987,0.2923635,4.67];
ybs[6702]=['67 Oph',4.7205604,0.0511803,3.97];
ybs[6703]=['6 Sgr',4.7246461,-0.2994234,6.28];
ybs[6704]=['',4.7271871,-0.3975682,5.77];
ybs[6705]=['',4.6635114,1.3666,6.24];
ybs[6706]=['',4.7105243,0.7936998,6.48];
ybs[6707]=['',4.7214535,0.1094185,6.34];
ybs[6708]=['',4.7190391,0.3404513,6.5];
ybs[6709]=['χ Oct',5.0113815,-1.5270431,5.28];
ybs[6710]=['',4.7213857,0.2634441,6.26];
ybs[6711]=['68 Oph',4.7254641,0.0228061,4.45];
ybs[6712]=['7 Sgr',4.7313845,-0.4237673,5.34];
ybs[6713]=['ψ2 Dra',4.6895395,1.2566721,5.45];
ybs[6714]=['',4.7189565,0.5797027,5.99];
ybs[6715]=['',4.7320782,-0.3964705,6.74];
ybs[6716]=['',4.715169,0.7941521,5.67];
ybs[6717]=['95 Her',4.7234634,0.3769297,5.18];
ybs[6718]=['95 Her',4.7234998,0.3769345,4.96];
ybs[6719]=['',4.77646,-1.3244199,5.86];
ybs[6720]=['',4.730185,-0.0934898,6.76];
ybs[6721]=['τ Oph',4.7316637,-0.1427338,5.94];
ybs[6722]=['τ Oph',4.7316565,-0.1427387,5.24];
ybs[6723]=['',4.684635,1.3119166,6.36];
ybs[6724]=['9 Sgr',4.7358462,-0.4251239,5.97];
ybs[6725]=['',4.7232776,0.5814147,6.15];
ybs[6726]=['96 Her',4.7273687,0.3636453,5.28];
ybs[6727]=['',4.7407092,-0.6265386,6];
ybs[6728]=['',4.7569215,-1.1265165,6.41];
ybs[6729]=['97 Her',4.7277769,0.4001141,6.21];
ybs[6730]=['γ1 Sgr',4.7411299,-0.5162081,4.69];
ybs[6731]=['θ Ara',4.7496486,-0.8741857,3.66];
ybs[6732]=['',4.7311788,0.3423516,6.5];
ybs[6733]=['π Pav',4.7601104,-1.1111207,4.35];
ybs[6734]=['γ2 Sgr',4.744609,-0.5309334,2.99];
ybs[6735]=['',4.737954,0.0335502,6.14];
ybs[6736]=['',4.7475013,-0.6285877,5.95];
ybs[6737]=['',4.7499212,-0.7578243,5.77];
ybs[6738]=['',4.7499212,-0.7578243,5.77];
ybs[6739]=['',4.7808452,-1.2856664,5.85];
ybs[6740]=['70 Oph',4.741566,0.0436866,4.03];
ybs[6741]=['',4.7289252,0.8459001,6.21];
ybs[6742]=['',4.7371879,0.4179293,6.34];
ybs[6743]=['',4.7449342,-0.1452088,5.85];
ybs[6744]=['',4.745351,-0.0828557,5.77];
ybs[6745]=['',4.7446046,-0.0077257,6.34];
ybs[6746]=['',4.7423275,0.2095727,7.04];
ybs[6747]=['',4.7574064,-0.798691,6.15];
ybs[6748]=['',4.7652898,-1.0303277,6.38];
ybs[6749]=['ι Pav',4.7678706,-1.0820227,5.49];
ybs[6750]=['',4.7501747,-0.3741842,6.28];
ybs[6751]=['',4.7409363,0.3778731,6.15];
ybs[6752]=['',4.7364598,0.6996539,6.52];
ybs[6753]=['98 Her',4.743214,0.3878608,5.06];
ybs[6754]=['',4.7544274,-0.4965802,4.57];
ybs[6755]=['',4.7376038,0.7321634,6.34];
ybs[6756]=['',4.7417912,0.5625946,5.71];
ybs[6757]=['',4.7526652,-0.2993082,5.52];
ybs[6758]=['71 Oph',4.7493837,0.1525171,4.64];
ybs[6759]=['72 Oph',4.7495375,0.1670038,3.73];
ybs[6760]=['',4.7605376,-0.6399499,6.58];
ybs[6761]=['',4.7578418,-0.4444787,6.61];
ybs[6762]=['',4.7875066,-1.2346806,6.73];
ybs[6763]=['99 Her',4.7471152,0.5334844,5.04];
ybs[6764]=['',4.7513786,0.2282208,6.63];
ybs[6765]=['',4.7630215,-0.5709547,6.43];
ybs[6766]=['',4.7687871,-0.8291347,6.07];
ybs[6767]=['ο Her',4.749476,0.5020836,3.83];
ybs[6768]=['',4.763329,-0.5362024,5.53];
ybs[6769]=['100 Her',4.7508454,0.4556416,5.86];
ybs[6770]=['100 Her',4.7508456,0.4555737,5.9];
ybs[6771]=['ε Tel',4.7693256,-0.8019304,4.53];
ybs[6772]=['',4.7546196,0.2494101,6.37];
ybs[6773]=['',4.7608839,-0.2430937,6.39];
ybs[6774]=['',4.7683645,-0.7217294,5.86];
ybs[6775]=['102 Her',4.7551791,0.363377,4.36];
ybs[6776]=['',4.7671069,-0.5897946,6.16];
ybs[6777]=['δ UMi',4.5572102,1.5082088,4.36];
ybs[6778]=['',4.7450185,0.8870994,6.29];
ybs[6779]=['',4.7482507,0.7586307,5];
ybs[6780]=['',4.7460625,0.8676899,6.32];
ybs[6781]=['',4.7511805,0.635412,5.48];
ybs[6782]=['101 Her',4.7557536,0.3499538,5.1];
ybs[6783]=['73 Oph',4.7594298,0.0698024,5.73];
ybs[6784]=['',4.7849352,-1.1114302,6.47];
ybs[6785]=['',4.7609352,0.0545586,5.69];
ybs[6786]=['',4.7677985,-0.3461878,6.36];
ybs[6787]=['',4.7564802,0.5318922,6.38];
ybs[6788]=['',4.7642937,0.0581349,5.51];
ybs[6789]=['11 Sgr',4.7700641,-0.4135329,4.98];
ybs[6790]=['',4.771409,-0.504292,6.51];
ybs[6791]=['',4.7614304,0.2876831,6.09];
ybs[6792]=['',4.7775973,-0.7213051,5.47];
ybs[6793]=['',4.790962,-1.1003523,5.6];
ybs[6794]=['',4.7579759,0.6713142,6.4];
ybs[6795]=['',4.7596708,0.6365667,5.58];
ybs[6796]=['',4.802411,-1.1906231,6.33];
ybs[6797]=['40 Dra',4.7046135,1.3962689,6.04];
ybs[6798]=['41 Dra',4.7050306,1.3963281,5.68];
ybs[6799]=['24 UMi',4.5461403,1.5154063,5.79];
ybs[6800]=['μ Sgr',4.7788365,-0.3673964,3.86];
ybs[6801]=['',4.7754849,-0.0698734,6.59];
ybs[6802]=['',4.7675738,0.583886,5.88];
ybs[6803]=['104 Her',4.7683546,0.5482539,4.97];
ybs[6804]=['14 Sgr',4.7810566,-0.3788086,5.44];
ybs[6805]=['',4.7604912,0.947593,5.95];
ybs[6806]=['',4.7894593,-0.7713779,5.46];
ybs[6807]=['',4.7961106,-0.9776032,5.33];
ybs[6808]=['',4.7748366,0.3820261,6.12];
ybs[6809]=['',4.7950473,-0.891124,6.06];
ybs[6810]=['15 Sgr',4.7851544,-0.3616122,5.38];
ybs[6811]=['16 Sgr',4.7851383,-0.3556732,5.95];
ybs[6812]=['',4.7712488,0.7182857,6.36];
ybs[6813]=['',4.7863595,-0.3255342,6.07];
ybs[6814]=['',4.773031,0.6768672,6.04];
ybs[6815]=['',4.7621932,1.0544607,6.49];
ybs[6816]=['',4.8085297,-1.1148212,6.18];
ybs[6817]=['φ Oct',4.829703,-1.3095061,5.47];
ybs[6818]=['',4.7876856,-0.0629649,6.36];
ybs[6819]=['',4.7808262,0.5099201,6.56];
ybs[6820]=['η Sgr',4.7965511,-0.6414206,3.11];
ybs[6821]=['',4.7962764,-0.5950922,6.16];
ybs[6822]=['',4.7879896,0.0416735,6.01];
ybs[6823]=['',4.7950862,-0.4998919,6.19];
ybs[6824]=['',4.7950449,-0.4935506,6.4];
ybs[6825]=['',4.8595078,-1.3999997,5.95];
ybs[6826]=['',4.7936304,-0.3030458,5.75];
ybs[6827]=['',4.8014737,-0.737868,6.3];
ybs[6828]=['',4.7916672,-0.052304,6];
ybs[6829]=['',4.7949092,-0.3220572,6.54];
ybs[6830]=['',4.7978513,-0.4717854,4.65];
ybs[6831]=['',4.7942123,-0.1701323,6.31];
ybs[6832]=['',4.7923516,0.0177388,6.63];
ybs[6833]=['',4.7839888,0.7359875,5.59];
ybs[6834]=['',4.8005768,-0.4466849,6.51];
ybs[6835]=['',4.783296,0.7892187,6.29];
ybs[6836]=['',4.8003556,-0.3247689,6.84];
ybs[6837]=['',4.7783417,0.9878074,6.37];
ybs[6838]=['36 Dra',4.7735439,1.1240885,5.03];
ybs[6839]=['',4.7960367,0.2406461,6.3];
ybs[6840]=['',4.7961906,0.3166461,5.99];
ybs[6841]=['',4.7904799,0.7146613,6.11];
ybs[6842]=['',4.7959459,0.4067969,6.63];
ybs[6843]=['ξ Pav',4.8235708,-1.0730183,4.36];
ybs[6844]=['',4.8109732,-0.6540546,6.45];
ybs[6845]=['',4.8011597,0.126911,5.39];
ybs[6846]=['',4.8064464,-0.2760984,5.39];
ybs[6847]=['δ Sgr',4.8108343,-0.5203719,2.7];
ybs[6848]=['105 Her',4.8004707,0.4268693,5.27];
ybs[6849]=['',4.8128815,-0.4346178,6.25];
ybs[6850]=['',4.8171078,-0.6744509,5.1];
ybs[6851]=['',4.8119867,-0.3289398,5.75];
ybs[6852]=['',4.8151577,-0.4959609,6.16];
ybs[6853]=['37 Dra',4.7784427,1.2001736,5.95];
ybs[6854]=['74 Oph',4.8087832,0.0591664,4.86];
ybs[6855]=['',4.8032089,0.5179823,5.99];
ybs[6856]=['106 Her',4.8054799,0.3835144,4.95];
ybs[6857]=['η Ser',4.8109702,-0.0503675,3.26];
ybs[6858]=['',4.8194852,-0.6397566,5.34];
ybs[6859]=['',4.8338851,-1.0996541,6.14];
ybs[6860]=['κ Lyr',4.8028036,0.6296538,4.33];
ybs[6861]=['',4.8113446,0.0951022,6.13];
ybs[6862]=['',4.8220606,-0.6322262,5.55];
ybs[6863]=['',4.8262046,-0.7696084,5.25];
ybs[6864]=['108 Her',4.8079259,0.5213581,5.63];
ybs[6865]=['107 Her',4.8082634,0.5040995,5.12];
ybs[6866]=['',4.818811,-0.1781025,6.33];
ybs[6867]=['ε Sgr',4.8249558,-0.5998676,1.85];
ybs[6868]=['',4.8018852,0.8963976,6.3];
ybs[6869]=['',4.8196152,-0.2094488,5.73];
ybs[6870]=['',4.813474,0.4066398,5.41];
ybs[6871]=['',4.8159233,0.2101838,5.89];
ybs[6872]=['ζ Sct',4.8214688,-0.1556784,4.68];
ybs[6873]=['',4.8166626,0.3113762,5.25];
ybs[6874]=['',4.8072161,0.8680969,6.4];
ybs[6875]=['',4.8177251,0.2915062,6.22];
ybs[6876]=['18 Sgr',4.8284674,-0.5365375,5.6];
ybs[6877]=['',4.8302464,-0.6279017,6.15];
ybs[6878]=['',4.8229902,-0.0622849,6.38];
ybs[6879]=['',4.8091374,0.8575615,5.05];
ybs[6880]=['',4.8259427,-0.123224,6.31];
ybs[6881]=['',4.8325128,-0.59218,6.3];
ybs[6882]=['',4.8378686,-0.8395108,5.46];
ybs[6883]=['109 Her',4.8203258,0.3802044,3.84];
ybs[6884]=['21 Sgr',4.8293641,-0.3582494,4.81];
ybs[6885]=['α Tel',4.8380049,-0.8020099,3.51];
ybs[6886]=['',4.8267978,-0.0273013,6.15];
ybs[6887]=['',4.8697001,-1.2905856,5.89];
ybs[6888]=['',4.8273802,0.089012,6.74];
ybs[6889]=['',4.8204991,0.6763786,6.36];
ybs[6890]=['',4.8294387,0.1404558,5.65];
ybs[6891]=['μ Lyr',4.8216392,0.6897865,5.12];
ybs[6892]=['',4.8256155,0.4784017,6.27];
ybs[6893]=['ζ Tel',4.8463913,-0.8561392,4.13];
ybs[6894]=['',4.8303474,0.261492,6.37];
ybs[6895]=['',4.8406306,-0.5200978,5.92];
ybs[6896]=['',4.8522323,-1.0036457,5.76];
ybs[6897]=['',4.8400472,-0.4645684,6.31];
ybs[6898]=['',4.843923,-0.6802973,5.64];
ybs[6899]=['',4.8185234,0.9305245,6.32];
ybs[6900]=['',4.9183488,-1.4273535,6.27];
ybs[6901]=['λ Sgr',4.8410323,-0.4433942,2.81];
ybs[6902]=['',4.8416823,-0.4667025,6.27];
ybs[6903]=['',4.8476087,-0.7649424,6.36];
ybs[6904]=['ν Pav',4.8592429,-1.0866254,4.64];
ybs[6905]=['',4.8298708,0.5208865,5.83];
ybs[6906]=['59 Ser',4.8365783,0.0037112,5.21];
ybs[6907]=['',4.8405445,-0.3103717,6.2];
ybs[6908]=['φ Dra',4.8014208,1.2452926,4.22];
ybs[6909]=['',4.8475265,-0.6777678,6.63];
ybs[6910]=['',4.8510108,-0.8238298,5.7];
ybs[6911]=['39 Dra',4.8182775,1.026513,4.98];
ybs[6912]=['',4.8331123,0.4619115,6.53];
ybs[6913]=['',4.8391805,0.0657203,6.07];
ybs[6914]=['',4.845396,-0.4636298,6.5];
ybs[6915]=['χ Dra',4.8021151,1.2696423,3.57];
ybs[6916]=['',4.8396971,0.1084047,5.73];
ybs[6917]=['',4.8471172,-0.4404951,6.59];
ybs[6918]=['γ Sct',4.8458827,-0.253912,4.7];
ybs[6919]=['',4.8554648,-0.7312,6.04];
ybs[6920]=['',4.8484211,-0.2541824,5.96];
ybs[6921]=['',4.8504325,-0.3265606,5.66];
ybs[6922]=['δ1 Tel',4.8588647,-0.8010298,4.96];
ybs[6923]=['60 Ser',4.8474675,-0.0343356,5.39];
ybs[6924]=['',4.85501,-0.5754395,5.34];
ybs[6925]=['',4.8594538,-0.7590142,5.72];
ybs[6926]=['δ2 Tel',4.8600582,-0.7982732,5.07];
ybs[6927]=['',4.8678737,-1.0243106,6.44];
ybs[6928]=['',4.8500437,-0.0995857,6.28];
ybs[6929]=['',4.8489706,0.0712704,6.69];
ybs[6930]=['',4.8610059,-0.6926199,5.16];
ybs[6931]=['',4.8459559,0.416854,5.9];
ybs[6932]=['',4.8558275,-0.3208562,5.14];
ybs[6933]=['42 Dra',4.8260417,1.1445706,4.82];
ybs[6934]=['',4.8554491,-0.1880904,5.72];
ybs[6935]=['',4.8578239,-0.3334566,6.68];
ybs[6936]=['',4.8638907,-0.6959001,6.22];
ybs[6937]=['',4.834733,1.0396178,6.43];
ybs[6938]=['',4.8508959,0.3636184,6.5];
ybs[6939]=['θ CrA',4.8661958,-0.7381367,4.64];
ybs[6940]=['κ1 CrA',4.8654388,-0.675442,6.32];
ybs[6941]=['κ2 CrA',4.8654246,-0.6755439,5.65];
ybs[6942]=['',4.8715978,-0.9227713,6.22];
ybs[6943]=['',4.8527329,0.2957874,5.77];
ybs[6944]=['',4.8596221,-0.2552468,6.37];
ybs[6945]=['61 Ser',4.8573171,-0.0171693,5.94];
ybs[6946]=['',4.8578523,0.064213,6.43];
ybs[6947]=['',4.861268,-0.259107,5.5];
ybs[6948]=['',4.8675896,-0.5758898,5.28];
ybs[6949]=['24 Sgr',4.8667972,-0.4190879,5.49];
ybs[6950]=['',4.8653177,-0.2588889,5.76];
ybs[6951]=['',4.8637664,-0.1028209,6.36];
ybs[6952]=['',4.9642964,-1.453581,7.16];
ybs[6953]=['25 Sgr',4.8696641,-0.4223972,6.51];
ybs[6954]=['',4.8598237,0.4125379,5.84];
ybs[6955]=['',4.8632096,0.1446612,6.42];
ybs[6956]=['',4.8597322,0.5336155,5.48];
ybs[6957]=['',4.8730279,-0.3633679,6.48];
ybs[6958]=['',4.8712145,-0.1912192,5.14];
ybs[6959]=['',4.8621282,0.5395261,6.59];
ybs[6960]=['',4.8762656,-0.5179678,6.37];
ybs[6961]=['α Sct',4.8718263,-0.1435168,3.85];
ybs[6962]=['',4.8552933,0.9099242,6.56];
ybs[6963]=['',4.8667676,0.3575665,6.57];
ybs[6964]=['',4.8692223,0.1904613,6.4];
ybs[6965]=['',4.8707268,0.3180779,5.78];
ybs[6966]=['45 Dra',4.8563672,0.9959722,4.77];
ybs[6967]=['',4.8490452,1.1423995,6.59];
ybs[6968]=['',4.871767,0.4123673,5.61];
ybs[6969]=['',4.8737351,0.2966561,6.21];
ybs[6970]=['ζ Pav',4.9125907,-1.2461958,4.01];
ybs[6971]=['',4.8629282,0.9140971,5.36];
ybs[6972]=['',4.8699516,0.6017707,6.1];
ybs[6973]=['',4.8765916,0.1596009,5.39];
ybs[6974]=['',4.8916716,-0.8357673,5.86];
ybs[6975]=['',4.8775158,0.1168326,5.45];
ybs[6976]=['',4.8841861,-0.3730622,5.94];
ybs[6977]=['',4.8845936,-0.2440277,6.47];
ybs[6978]=['',4.8869265,-0.4098338,5.81];
ybs[6979]=['',4.8927817,-0.7533217,5.37];
ybs[6980]=['',4.8797523,0.1997364,6.42];
ybs[6981]=['',4.8819411,-0.0050055,5.75];
ybs[6982]=['',4.9373375,-1.3585273,6.39];
ybs[6983]=['',4.8792827,0.2831041,6.29];
ybs[6984]=['',4.9076937,-1.1277844,6.37];
ybs[6985]=['',4.876099,0.5845257,5.42];
ybs[6986]=['',4.8884597,-0.3670162,5.86];
ybs[6987]=['',4.8855283,-0.0553355,6.49];
ybs[6988]=['',4.8851085,-0.0190286,6.66];
ybs[6989]=['α Lyr',4.8771536,0.6772879,0.03];
ybs[6990]=['',4.8848363,0.1545829,6.4];
ybs[6991]=['',4.876055,0.7547492,6.2];
ybs[6992]=['',4.9131004,-1.1261669,5.78];
ybs[6993]=['',4.9015958,-0.8389778,6.49];
ybs[6994]=['',4.8369691,1.3537515,5.64];
ybs[6995]=['',4.8927574,-0.1355505,5.84];
ybs[6996]=['',4.8905046,0.0922926,6.38];
ybs[6997]=['',4.8821999,0.6927364,6.04];
ybs[6998]=['',4.8914869,0.1288451,6.28];
ybs[6999]=['26 Sgr',4.9015496,-0.4155301,6.23];
ybs[7000]=['',4.9211413,-1.1317381,4.79];
ybs[7001]=['',4.8707565,1.143367,6.06];
ybs[7002]=['',4.9004623,-0.2537547,6.42];
ybs[7003]=['',4.9192753,-1.0658312,6.04];
ybs[7004]=['',4.8911459,0.5388429,6.36];
ybs[7005]=['',4.8884105,0.7148633,6.25];
ybs[7006]=['',4.8772284,1.0916857,5.74];
ybs[7007]=['',4.8914322,0.6700541,6.45];
ybs[7008]=['δ Sct',4.9026923,-0.1575524,4.72];
ybs[7009]=['λ CrA',4.9107468,-0.6684126,5.13];
ybs[7010]=['',4.9195138,-0.9922934,6.22];
ybs[7011]=['',4.9059591,-0.3361168,6.35];
ybs[7012]=['',4.9040424,-0.1230063,6.15];
ybs[7013]=['',4.8034416,1.451916,6.17];
ybs[7014]=['',4.912155,-0.6403967,6.32];
ybs[7015]=['',4.9424175,-1.2734706,6.06];
ybs[7016]=['',4.8888544,0.911409,6];
ybs[7017]=['',4.9129382,-0.6216034,4.87];
ybs[7018]=['',4.8983261,0.5522697,6.41];
ybs[7019]=['',4.9159491,-0.6921856,5.43];
ybs[7020]=['ε Sct',4.9081064,-0.1439746,4.9];
ybs[7021]=['',4.900093,0.6068834,6.47];
ybs[7022]=['',4.9095076,-0.1185477,6.31];
ybs[7023]=['',4.9145493,-0.4360559,5.83];
ybs[7024]=['θ Pav',4.9350558,-1.135308,5.73];
ybs[7025]=['',4.9256619,-0.8738187,6.54];
ybs[7026]=['',4.9164723,-0.3660684,6.36];
ybs[7027]=['φ Sgr',4.9182698,-0.4706,3.17];
ybs[7028]=['4 Aql',4.9133985,0.0364228,5.02];
ybs[7029]=['',4.9047912,0.6863702,6.45];
ybs[7030]=['',4.8919447,1.0956134,6.09];
ybs[7031]=['',4.9063826,0.6384888,6.01];
ybs[7032]=['',4.9077869,0.557683,5.7];
ybs[7033]=['',4.9194979,-0.3417138,6.42];
ybs[7034]=['28 Sgr',4.9210385,-0.3903324,5.37];
ybs[7035]=['',4.9117576,0.4121846,6.31];
ybs[7036]=['',4.9160577,0.0964637,5.83];
ybs[7037]=['46 Dra',4.9004737,0.9697889,5.04];
ybs[7038]=['μ CrA',4.9281677,-0.7047187,5.24];
ybs[7039]=['ε1 Lyr',4.9094,0.6928339,5.06];
ybs[7040]=['ε1 Lyr',4.9093928,0.6928533,6.02];
ybs[7041]=['ε2 Lyr',4.9095858,0.6918405,5.14];
ybs[7042]=['ε2 Lyr',4.9095858,0.6918357,5.37];
ybs[7043]=['',4.9221452,-0.1762258,5.71];
ybs[7044]=['ζ1 Lyr',4.9114296,0.6567975,4.36];
ybs[7045]=['ζ2 Lyr',4.9115612,0.6566135,5.73];
ybs[7046]=['',4.9158714,0.3841864,6.51];
ybs[7047]=['5 Aql',4.920692,-0.0162983,5.9];
ybs[7048]=['',4.9043963,0.9406944,6.11];
ybs[7049]=['110 Her',4.9162298,0.3590786,4.19];
ybs[7050]=['η1 CrA',4.9331989,-0.7618473,5.49];
ybs[7051]=['β Sct',4.9239006,-0.0823712,4.22];
ybs[7052]=['',4.9177367,0.465824,4.83];
ybs[7053]=['',4.9360522,-0.7990213,5.81];
ybs[7054]=['',4.9252851,-0.0990746,5.2];
ybs[7055]=['',4.9208019,0.3269655,6.17];
ybs[7056]=['η2 CrA',4.9364144,-0.7575445,5.61];
ybs[7057]=['111 Her',4.922273,0.3178156,4.36];
ybs[7058]=['',4.9345328,-0.6059614,6.62];
ybs[7059]=['',4.9105601,0.9585938,6.23];
ybs[7060]=['',4.9313838,-0.3241456,6.47];
ybs[7061]=['',4.9174667,0.7237743,6.07];
ybs[7062]=['λ Pav',4.9500944,-1.0848275,4.22];
ybs[7063]=['',4.9069529,1.0659484,5.99];
ybs[7064]=['',4.9273307,0.0745278,6.21];
ybs[7065]=['',4.9350582,-0.3335764,6.75];
ybs[7066]=['29 Sgr',4.9354469,-0.354214,5.24];
ybs[7067]=['',4.927487,0.4109025,6.15];
ybs[7068]=['',4.9204308,0.8088374,6.52];
ybs[7069]=['',4.9256841,0.5547624,6.06];
ybs[7070]=['',4.8995146,1.2360114,6.44];
ybs[7071]=['',4.9348924,-0.1026787,5.99];
ybs[7072]=['',4.918622,0.9253007,5.88];
ybs[7073]=['',4.9343343,0.0151057,6.25];
ybs[7074]=['',4.9303745,0.3378573,5.88];
ybs[7075]=['κ Tel',4.9506144,-0.9088954,5.17];
ybs[7076]=['30 Sgr',4.9406511,-0.3862725,6.61];
ybs[7077]=['',4.9378105,-0.1374865,6.8];
ybs[7078]=['',4.9231659,0.8570151,6.4];
ybs[7079]=['',4.931631,0.4376546,6.59];
ybs[7080]=['',4.9491794,-0.8126917,5.54];
ybs[7081]=['',4.9529573,-0.9058115,6.31];
ybs[7082]=['',4.9406875,-0.1700594,5.83];
ybs[7083]=['',4.9519081,-0.8434858,6.19];
ybs[7084]=['',4.9258216,0.8516544,6.12];
ybs[7085]=['',4.9515406,-0.8125215,6.19];
ybs[7086]=['',4.9334386,0.5525504,6.64];
ybs[7087]=['',4.9388995,0.1921029,6.55];
ybs[7088]=['ν1 Lyr',4.9335098,0.5732085,5.91];
ybs[7089]=['8 Aql',4.9421385,-0.0573705,6.1];
ybs[7090]=['ν2 Lyr',4.9340339,0.5686379,5.25];
ybs[7091]=['',4.9479953,-0.4645871,6.29];
ybs[7092]=['',4.9487547,-0.5122184,6.13];
ybs[7093]=['',4.9491695,-0.5358619,6.63];
ybs[7094]=['β Lyr',4.9348528,0.582811,3.45];
ybs[7095]=['κ Pav',4.9718452,-1.1728498,4.44];
ybs[7096]=['',4.9585774,-0.8699753,6.6];
ybs[7097]=['',4.9443219,0.2442864,6.14];
ybs[7098]=['',4.9496523,-0.166582,6.34];
ybs[7099]=['',4.9706733,-1.0954896,6.48];
ybs[7100]=['',4.9417332,0.5029051,6.18];
ybs[7101]=['112 Her',4.9450432,0.3744852,5.48];
ybs[7102]=['33 Sgr',4.9543986,-0.3722344,5.69];
ybs[7103]=['',4.9413321,0.638255,6.09];
ybs[7104]=['ν1 Sgr',4.955198,-0.3964102,4.83];
ybs[7105]=['',4.909421,1.2935061,5.27];
ybs[7106]=['',4.9432333,0.7228161,6.28];
ybs[7107]=['',4.9572665,-0.2717546,5.1];
ybs[7108]=['ν2 Sgr',4.9593314,-0.3951159,4.99];
ybs[7109]=['σ Sgr',4.9601476,-0.4583872,2.02];
ybs[7110]=['',4.9655714,-0.7448523,5.36];
ybs[7111]=['',4.9398652,0.9251223,5.51];
ybs[7112]=['50 Dra',4.9112088,1.317044,5.35];
ybs[7113]=['ο Dra',4.9373523,1.0370512,4.66];
ybs[7114]=['',4.9607831,-0.2852486,5.79];
ybs[7115]=['ω Pav',4.9776016,-1.0500848,5.14];
ybs[7116]=['',4.9632449,-0.4038723,5.93];
ybs[7117]=['',4.9669177,-0.651173,5.38];
ybs[7118]=['',4.9852365,-1.162691,6.01];
ybs[7119]=['δ1 Lyr',4.9505516,0.6458346,5.58];
ybs[7120]=['',4.9532214,0.4876746,5.62];
ybs[7121]=['113 Her',4.9557875,0.3957983,4.59];
ybs[7122]=['λ Tel',4.9759999,-0.9233426,4.87];
ybs[7123]=['',4.9595821,0.1160349,5.57];
ybs[7124]=['',4.9710207,-0.694448,6.31];
ybs[7125]=['',4.9472812,0.8855784,4.92];
ybs[7126]=['',4.9524994,0.7200838,7.3];
ybs[7127]=['δ2 Lyr',4.9539464,0.6445722,4.3];
ybs[7128]=['',4.9557471,0.5934332,6.02];
ybs[7129]=['θ1 Ser',4.9630037,0.0739512,4.62];
ybs[7130]=['θ2 Ser',4.9631056,0.0739223,4.98];
ybs[7131]=['',4.9639368,-0.0308298,6.22];
ybs[7132]=['',4.9639756,0.0437154,6.15];
ybs[7133]=['ξ1 Sgr',4.9689363,-0.3599253,5.08];
ybs[7134]=['',4.9552297,0.726674,5.44];
ybs[7135]=['',4.9617472,0.314654,6.63];
ybs[7136]=['',4.9619096,0.3165791,5.69];
ybs[7137]=['η Sct',4.9670819,-0.1014407,4.83];
ybs[7138]=['ξ2 Sgr',4.970651,-0.3677801,3.51];
ybs[7139]=['',4.9738738,-0.5410748,6.12];
ybs[7140]=['ε CrA',4.9758363,-0.6470365,4.87];
ybs[7141]=['',4.9488571,1.0038921,6.22];
ybs[7142]=['',4.9542645,0.8533299,5.77];
ybs[7143]=['',4.973498,-0.4335727,6.62];
ybs[7144]=['',4.9780079,-0.6893947,6.49];
ybs[7145]=['13 Lyr',4.9570841,0.7675773,4.04];
ybs[7146]=['64 Ser',4.9676819,0.0448439,5.57];
ybs[7147]=['',4.9736837,-0.3926054,6.14];
ybs[7148]=['',4.9038033,1.3957217,6.39];
ybs[7149]=['',5.0007917,-1.1993419,5.88];
ybs[7150]=['',4.9652009,0.5748285,5.22];
ybs[7151]=['',4.9724191,0.1095195,6.21];
ybs[7152]=['',4.9780245,-0.3234365,6.37];
ybs[7153]=['',4.9712891,0.3036077,5.38];
ybs[7154]=['',4.9775582,-0.223493,5.53];
ybs[7155]=['10 Aql',4.9737809,0.2433267,5.89];
ybs[7156]=['',4.982538,-0.4346961,6.36];
ybs[7157]=['',4.985995,-0.6461986,6.69];
ybs[7158]=['',4.986075,-0.6462178,6.4];
ybs[7159]=['',4.9733915,0.3460824,6.5];
ybs[7160]=['11 Aql',4.9751604,0.2383702,5.23];
ybs[7161]=['',4.9761647,0.1776058,6.75];
ybs[7162]=['',4.9692608,0.6684702,5.89];
ybs[7163]=['48 Dra',4.961818,1.0096477,5.66];
ybs[7164]=['ε Aql',4.9774041,0.2636099,4.02];
ybs[7165]=['',4.9910819,-0.7308212,6.23];
ybs[7166]=['γ Lyr',4.9735792,0.5711485,3.24];
ybs[7167]=['',4.9723501,0.7105933,6.22];
ybs[7168]=['υ Dra',4.9483833,1.2449295,4.82];
ybs[7169]=['',4.9774794,0.4584283,5.27];
ybs[7170]=['',4.9877377,-0.3954723,6.24];
ybs[7171]=['',4.9785616,0.398813,6.29];
ybs[7172]=['',4.964918,1.0168156,6.46];
ybs[7173]=['',4.9743199,0.6850916,6.41];
ybs[7174]=['',4.9870928,-0.2660911,6.32];
ybs[7175]=['',4.9590623,1.1395487,5.63];
ybs[7176]=['ζ CrA',4.9953252,-0.7340452,4.75];
ybs[7177]=['',4.9997542,-0.8897775,5.93];
ybs[7178]=['',4.9634288,1.0896174,6.45];
ybs[7179]=['λ Lyr',4.9782876,0.5616667,4.93];
ybs[7180]=['12 Aql',4.9872275,-0.099523,4.02];
ybs[7181]=['ζ Sgr',4.9923778,-0.5208591,2.6];
ybs[7182]=['',4.9914682,-0.4330127,5.65];
ybs[7183]=['',4.9724506,0.8874012,6.3];
ybs[7184]=['',4.9958472,-0.6669891,5.74];
ybs[7185]=['',4.9836258,0.3376507,6.39];
ybs[7186]=['',4.9423801,1.3232884,6.22];
ybs[7187]=['',4.9848003,0.3642503,6.69];
ybs[7188]=['',4.9790728,0.710696,6.65];
ybs[7189]=['',4.9841616,0.4595054,5.69];
ybs[7190]=['',4.993832,-0.3352395,6.05];
ybs[7191]=['',4.9821436,0.59059,6.01];
ybs[7192]=['',4.9940583,-0.3327616,6.37];
ybs[7193]=['',4.985497,0.4374203,6.72];
ybs[7194]=['',4.9866908,0.3892178,6.4];
ybs[7195]=['',4.9896268,0.1468025,6.3];
ybs[7196]=['14 Aql',4.9925047,-0.0639061,5.42];
ybs[7197]=['',4.9778709,0.8825994,5.38];
ybs[7198]=['',5.0003204,-0.5412032,5.5];
ybs[7199]=['',4.9860173,0.5874427,6.39];
ybs[7200]=['ρ Tel',5.010245,-0.9128309,5.16];
ybs[7201]=['',4.9950256,0.0324033,5.83];
ybs[7202]=['16 Lyr',4.983494,0.8197994,5.01];
ybs[7203]=['',4.9914002,0.3438014,6.09];
ybs[7204]=['ο Sgr',5.0010124,-0.3787934,3.77];
ybs[7205]=['49 Dra',4.9794595,0.972047,5.48];
ybs[7206]=['',4.997765,0.0587933,6.73];
ybs[7207]=['',4.9991073,-0.0985553,6.9];
ybs[7208]=['',5.028496,-1.1935093,5.33];
ybs[7209]=['',4.9949567,0.3718512,6.52];
ybs[7210]=['',5.0124445,-0.842285,5.97];
ybs[7211]=['',4.9685913,1.2141509,6.52];
ybs[7212]=['15 Aql',5.0014697,-0.0696888,5.42];
ybs[7213]=['γ CrA',5.0093926,-0.6461887,4.93];
ybs[7214]=['γ CrA',5.0093926,-0.6461887,4.99];
ybs[7215]=['σ Oct',5.6221333,-1.5505087,5.47];
ybs[7216]=['',4.9859387,0.9127687,6.31];
ybs[7217]=['',5.0051189,-0.2726434,5.97];
ybs[7218]=['',5.0028961,-0.0257274,6.53];
ybs[7219]=['',5.0114312,-0.6592207,6.16];
ybs[7220]=['',5.0216917,-0.9717877,6.49];
ybs[7221]=['τ Sgr',5.011142,-0.4822492,3.32];
ybs[7222]=['ζ Aql',5.0027088,0.2426367,2.99];
ybs[7223]=['λ Aql',5.0071249,-0.0845306,3.44];
ybs[7224]=['',4.9999171,0.5547156,5.56];
ybs[7225]=['',5.0000013,0.5370685,6.06];
ybs[7226]=['',5.0103129,-0.2825557,6.03];
ybs[7227]=['',5.0136994,-0.4991102,6.04];
ybs[7228]=['',5.0116003,-0.3263413,6.29];
ybs[7229]=['δ CrA',5.0180334,-0.7060923,4.59];
ybs[7230]=['',5.0071312,0.1443263,6.09];
ybs[7231]=['',5.003595,0.5229103,6.31];
ybs[7232]=['',5.0108461,0.011893,6.56];
ybs[7233]=['',5.0166784,-0.4296391,6.3];
ybs[7234]=['',4.9606936,1.345381,6.54];
ybs[7235]=['18 Aql',5.0096623,0.1939238,5.09];
ybs[7236]=['',5.0165833,-0.3359681,5.54];
ybs[7237]=['',5.0076122,0.4239445,5.77];
ybs[7238]=['51 Dra',4.9980565,0.9326157,5.38];
ybs[7239]=['',4.9994874,0.8719975,6.43];
ybs[7240]=['',5.0073438,0.5003507,5.55];
ybs[7241]=['α CrA',5.022761,-0.6608386,4.11];
ybs[7242]=['',5.0237151,-0.6944098,6.46];
ybs[7243]=['',5.0232417,-0.6304734,6.56];
ybs[7244]=['',5.0251734,-0.7304381,5.88];
ybs[7245]=['',5.0050707,0.7234916,6.49];
ybs[7246]=['β CrA',5.0252816,-0.6859026,4.11];
ybs[7247]=['',5.0136876,0.2948472,6.07];
ybs[7248]=['17 Lyr',5.0106206,0.567956,5.23];
ybs[7249]=['ι Lyr',5.009865,0.6307621,5.28];
ybs[7250]=['',5.0139297,0.3794189,6.23];
ybs[7251]=['π Sgr',5.023139,-0.3662105,2.89];
ybs[7252]=['',5.0232513,-0.3449172,6.13];
ybs[7253]=['19 Aql',5.0186875,0.1067116,5.22];
ybs[7254]=['',5.0168081,0.2948203,6.48];
ybs[7255]=['',5.0296275,-0.6800314,6.36];
ybs[7256]=['',5.0227059,-0.0067503,6.34];
ybs[7257]=['',5.0303182,-0.5141744,6.3];
ybs[7258]=['',5.0381161,-0.8804012,6.13];
ybs[7259]=['',5.0176888,0.6046047,6.74];
ybs[7260]=['',5.0344786,-0.6551978,6.57];
ybs[7261]=['τ Pav',5.0574526,-1.2068105,6.27];
ybs[7262]=['',5.013492,0.9157018,5.81];
ybs[7263]=['',5.0349601,-0.3772569,6.41];
ybs[7264]=['',5.0384824,-0.4514014,5.8];
ybs[7265]=['',5.0599084,-1.1626619,5.53];
ybs[7266]=['20 Aql',5.0353013,-0.1378205,5.34];
ybs[7267]=['',5.0287695,0.4673641,6.36];
ybs[7268]=['',5.0459069,-0.7880019,5.92];
ybs[7269]=['',5.0380085,-0.2136153,5.51];
ybs[7270]=['19 Lyr',5.0296328,0.5467352,5.98];
ybs[7271]=['',5.027413,0.7063557,6.18];
ybs[7272]=['',5.0338435,0.2947717,6.73];
ybs[7273]=['',5.0338023,0.3769429,5.93];
ybs[7274]=['21 Aql',5.0394046,0.0407897,5.15];
ybs[7275]=['',5.0393789,0.0970282,6.49];
ybs[7276]=['',5.0533403,-0.7927509,5.4];
ybs[7277]=['55 Dra',5.0171674,1.1522579,6.25];
ybs[7278]=['',5.0485407,-0.4212281,6.25];
ybs[7279]=['ψ Sgr',5.0485338,-0.4400341,4.85];
ybs[7280]=['',5.029677,0.870859,6.75];
ybs[7281]=['',5.029706,0.8708929,6.57];
ybs[7282]=['53 Dra',5.0271383,0.9931144,5.12];
ybs[7283]=['',5.0533526,-0.5842753,6.25];
ybs[7284]=['',5.0619031,-0.9309674,6.38];
ybs[7285]=['η Lyr',5.0378608,0.6839856,4.39];
ybs[7286]=['',5.0444765,0.3533856,6];
ybs[7287]=['',5.0459713,0.2640327,5.57];
ybs[7288]=['1 Sge',5.0455069,0.3713454,5.64];
ybs[7289]=['',5.0455921,0.53356,5.85];
ybs[7290]=['22 Aql',5.051544,0.0851679,5.59];
ybs[7291]=['43 Sgr',5.0573838,-0.3299951,4.96];
ybs[7292]=['',5.0480926,0.4799692,6.54];
ybs[7293]=['1 Vul',5.0495442,0.3741131,4.77];
ybs[7294]=['',5.0508455,0.2546382,5.63];
ybs[7295]=['',5.0482435,0.4881969,6.16];
ybs[7296]=['54 Dra',5.0368308,1.0078989,4.99];
ybs[7297]=['δ Dra',5.0289376,1.18166,3.07];
ybs[7298]=['',5.043828,0.8746724,6.27];
ybs[7299]=['59 Dra',5.0101172,1.3369368,5.13];
ybs[7300]=['',5.0572687,0.0362581,6.19];
ybs[7301]=['θ Lyr',5.0493214,0.6663399,4.36];
ybs[7302]=['ω1 Aql',5.0569461,0.2031743,5.28];
ybs[7303]=['',5.0670679,-0.617401,5.59];
ybs[7304]=['',5.0631991,-0.2703498,6.06];
ybs[7305]=['2 Vul',5.0560601,0.4026689,5.43];
ybs[7306]=['23 Aql',5.0605274,0.0197477,5.1];
ybs[7307]=['',5.090365,-1.1924347,6.34];
ybs[7308]=['24 Aql',5.0618946,0.0067287,6.41];
ybs[7309]=['',5.0508073,0.8210771,6];
ybs[7310]=['',5.0702231,-0.5544989,6.58];
ybs[7311]=['',5.0569157,0.5422393,6.68];
ybs[7312]=['',5.0616593,0.1686757,6.32];
ybs[7313]=['',5.0609354,0.3430766,6.58];
ybs[7314]=['',5.0706236,-0.3901699,5.58];
ybs[7315]=['κ Cyg',5.0512712,0.9322466,3.77];
ybs[7316]=['η Tel',5.0824817,-0.9490197,5.05];
ybs[7317]=['',5.0750298,-0.6097473,6.48];
ybs[7318]=['28 Aql',5.0649369,0.216796,5.53];
ybs[7319]=['ω2 Aql',5.0659677,0.2021424,6.02];
ybs[7320]=['26 Aql',5.0695372,-0.0936982,5.01];
ybs[7321]=['',5.0783625,-0.7324759,6.34];
ybs[7322]=['',5.0613701,0.5835555,6.6];
ybs[7323]=['27 Aql',5.0695639,-0.014746,5.49];
ybs[7324]=['β1 Sgr',5.0806265,-0.7751056,4.01];
ybs[7325]=['',5.0609418,0.6543522,6.22];
ybs[7326]=['',5.0747717,-0.3348625,6.26];
ybs[7327]=['ρ1 Sgr',5.0749516,-0.3106554,3.93];
ybs[7328]=['',5.0583167,0.8659589,6.31];
ybs[7329]=['υ Sgr',5.0751039,-0.2776294,4.61];
ybs[7330]=['β2 Sgr',5.0831885,-0.7810486,4.29];
ybs[7331]=['ρ2 Sgr',5.0757348,-0.3187016,5.87];
ybs[7332]=['',5.0636922,0.6523561,6.31];
ybs[7333]=['',5.0677648,0.6149374,6.31];
ybs[7334]=['',5.0774857,-0.1422929,6.31];
ybs[7335]=['α Sgr',5.0857946,-0.7080247,3.97];
ybs[7336]=['',5.0772329,-0.0035637,5.83];
ybs[7337]=['',5.0880627,-0.7622463,6.17];
ybs[7338]=['',5.0620998,0.9498552,6.26];
ybs[7339]=['τ Dra',5.0399696,1.2810643,4.45];
ybs[7340]=['',5.0806462,-0.1283086,6.32];
ybs[7341]=['',5.0787905,0.1738627,6.35];
ybs[7342]=['',5.0877474,-0.4854803,6.04];
ybs[7343]=['',5.0645732,1.0069189,5.91];
ybs[7344]=['',5.0800266,0.2612727,6.64];
ybs[7345]=['3 Vul',5.0782631,0.4592139,5.18];
ybs[7346]=['',5.0766136,0.585844,6.06];
ybs[7347]=['',5.0902911,-0.5106706,5.93];
ybs[7348]=['',5.06125,1.1246448,6.52];
ybs[7349]=['χ1 Sgr',5.090959,-0.4268833,5.03];
ybs[7350]=['χ3 Sgr',5.0918932,-0.4173449,5.43];
ybs[7351]=['',5.0826133,0.3545375,6.4];
ybs[7352]=['',5.0696027,1.0090536,6.43];
ybs[7353]=['',5.0890552,-0.0843754,6.52];
ybs[7354]=['',5.0908682,-0.2416744,5.69];
ybs[7355]=['',5.0810395,0.5806906,6.37];
ybs[7356]=['2 Sge',5.0853175,0.296482,6.25];
ybs[7357]=['',5.1040424,-0.9472557,5.69];
ybs[7358]=['π Dra',5.0648925,1.1477596,4.59];
ybs[7359]=['2 Cyg',5.0836762,0.5178498,4.97];
ybs[7360]=['31 Aql',5.0881496,0.209338,5.16];
ybs[7361]=['',5.0848306,0.4910857,6.53];
ybs[7362]=['50 Sgr',5.0953922,-0.3791919,5.59];
ybs[7363]=['',5.0831891,0.6370648,6.36];
ybs[7364]=['δ Aql',5.0908034,0.0552358,3.36];
ybs[7365]=['',5.0945148,-0.2618444,5.72];
ybs[7366]=['',5.0954755,-0.2530817,6.7];
ybs[7367]=['',5.0984832,-0.5182301,5.67];
ybs[7368]=['',5.0790768,0.8782517,6.51];
ybs[7369]=['',5.0820312,0.7581207,5.84];
ybs[7370]=['',5.1212553,-1.1934623,5.96];
ybs[7371]=['',5.0895601,0.3546741,6.31];
ybs[7372]=['4 Vul',5.0900319,0.3464236,5.16];
ybs[7373]=['',5.0895967,0.4356819,6.19];
ybs[7374]=['ν Aql',5.0953614,0.0067934,4.66];
ybs[7375]=['',5.1132139,-0.9667155,6.13];
ybs[7376]=['',5.0943539,0.2281916,5.74];
ybs[7377]=['5 Vul',5.0932632,0.3516522,5.63];
ybs[7378]=['',5.0943998,0.3480525,5.81];
ybs[7379]=['',5.1100221,-0.7573593,5.71];
ybs[7380]=['μ Tel',5.1162178,-0.960925,6.3];
ybs[7381]=['λ UMi',4.395063,1.5531871,6.38];
ybs[7382]=['4 Cyg',5.0921476,0.634743,5.15];
ybs[7383]=['',5.0993796,0.2501697,6.32];
ybs[7384]=['',5.1032392,0.052044,5.85];
ybs[7385]=['',5.1111172,-0.4700698,5.52];
ybs[7386]=['',5.1130114,-0.559194,6.6];
ybs[7387]=['35 Aql',5.1062006,0.0349463,5.8];
ybs[7388]=['',5.0885948,1.0136386,6.6];
ybs[7389]=['',5.1080381,-0.1220279,6.61];
ybs[7390]=['',5.0984091,0.6630895,6.34];
ybs[7391]=['',5.1075023,0.0052057,6.25];
ybs[7392]=['α Vul',5.103904,0.431389,4.44];
ybs[7393]=['8 Vul',5.1049686,0.4331996,5.81];
ybs[7394]=['',5.1072371,0.2556556,5.56];
ybs[7395]=['ι1 Cyg',5.0964987,0.9140552,5.75];
ybs[7396]=['7 Vul',5.1069096,0.3548578,6.33];
ybs[7397]=['',5.1153485,-0.3710418,6.13];
ybs[7398]=['',5.1261209,-0.92732,5.75];
ybs[7399]=['',5.1112182,0.0516059,6.09];
ybs[7400]=['',5.0907754,1.0927076,6.38];
ybs[7401]=['36 Aql',5.1135673,-0.0477518,5.03];
ybs[7402]=['',5.1128411,0.0610391,6.05];
ybs[7403]=['',5.1274287,-0.789193,5.61];
ybs[7404]=['β1 Cyg',5.1125498,0.4889119,3.08];
ybs[7405]=['β2 Cyg',5.1126951,0.4890092,5.11];
ybs[7406]=['',5.1123684,0.6332313,6.25];
ybs[7407]=['ι2 Cyg',5.1064895,0.9037651,3.79];
ybs[7408]=['',5.115416,0.4654871,5.87];
ybs[7409]=['',5.1304675,-0.6977791,5.89];
ybs[7410]=['',5.0620179,1.3901536,6.05];
ybs[7411]=['ι Tel',5.1357499,-0.8385194,4.9];
ybs[7412]=['',5.0259033,1.4574507,6.53];
ybs[7413]=['8 Cyg',5.1167931,0.6022516,4.74];
ybs[7414]=['',5.1136896,0.8789438,5.53];
ybs[7415]=['',5.1126969,0.9736311,6.37];
ybs[7416]=['μ Aql',5.1281297,0.1297417,4.45];
ybs[7417]=['37 Aql',5.133319,-0.1833453,5.12];
ybs[7418]=['51 Sgr',5.137869,-0.4304552,5.65];
ybs[7419]=['',5.1348345,-0.1292366,6.34];
ybs[7420]=['',5.1352981,-0.2128807,6.27];
ybs[7421]=['',5.1508077,-1.0109992,6.18];
ybs[7422]=['',5.1586142,-1.1628664,6.39];
ybs[7423]=['',5.1245496,0.677478,6.61];
ybs[7424]=['9 Vul',5.1297741,0.3460697,5];
ybs[7425]=['',5.1341046,0.0518159,6.38];
ybs[7426]=['',5.1393848,-0.3280639,6.11];
ybs[7427]=['52 Sgr',5.140835,-0.4333189,4.6];
ybs[7428]=['9 Cyg',5.1304904,0.5151893,5.38];
ybs[7429]=['',5.1241453,0.8607423,5.96];
ybs[7430]=['',5.142063,-0.3172079,5.64];
ybs[7431]=['',5.1290244,0.7412022,5.35];
ybs[7432]=['',5.1368997,0.1955792,6.68];
ybs[7433]=['κ Aql',5.1409038,-0.12167,4.95];
ybs[7434]=['ι Aql',5.1399436,-0.0214706,4.36];
ybs[7435]=['',5.1205772,1.0509089,6.29];
ybs[7436]=['',5.1373378,0.2521581,6.38];
ybs[7437]=['',5.1085622,1.2399191,6.07];
ybs[7438]=['',5.1267422,0.8952039,5.73];
ybs[7439]=['',5.1364459,0.3951718,6.32];
ybs[7440]=['',5.1284857,0.8415919,6.67];
ybs[7441]=['',5.1441628,-0.2486215,5.47];
ybs[7442]=['',5.1660693,-1.1483396,6.09];
ybs[7443]=['',5.1401315,0.1977388,5.98];
ybs[7444]=['11 Cyg',5.1342599,0.6457728,6.05];
ybs[7445]=['',5.1386834,0.3558532,7.14];
ybs[7446]=['',5.1585356,-0.9487464,6.26];
ybs[7447]=['42 Aql',5.144725,-0.0801229,5.46];
ybs[7448]=['',5.1550651,-0.7892396,6.25];
ybs[7449]=['σ Dra',5.1149478,1.2167488,4.68];
ybs[7450]=['ε Sge',5.1417342,0.2883154,5.66];
ybs[7451]=['',5.1556574,-0.6872283,6.61];
ybs[7452]=['',5.1338338,0.8777999,6.52];
ybs[7453]=['',5.1405744,0.512952,6.43];
ybs[7454]=['',5.1391433,0.6709065,6.5];
ybs[7455]=['',5.1373536,0.7810528,5.17];
ybs[7456]=['θ Cyg',5.1360633,0.8774993,4.48];
ybs[7457]=['53 Sgr',5.154362,-0.4078806,6.34];
ybs[7458]=['',5.1489083,0.0600219,6.35];
ybs[7459]=['',5.1459243,0.3637229,6.48];
ybs[7460]=['',5.1556486,-0.4078923,5.97];
ybs[7461]=['σ Aql',5.1504755,0.0952132,5.17];
ybs[7462]=['',5.1510409,0.2902312,6.38];
ybs[7463]=['54 Sgr',5.157984,-0.2833528,6.2];
ybs[7464]=['',5.1427155,0.8611657,6.47];
ybs[7465]=['φ Cyg',5.1502254,0.5272797,4.69];
ybs[7466]=['α Sge',5.1539202,0.3154138,4.37];
ybs[7467]=['45 Aql',5.1573715,-0.0098216,5.67];
ybs[7468]=['',5.1516544,0.5940564,6.1];
ybs[7469]=['',5.1554531,0.3584008,6.5];
ybs[7470]=['14 Cyg',5.1497623,0.7483255,5.4];
ybs[7471]=['',5.1453818,0.9604717,5.82];
ybs[7472]=['',5.1561405,0.4149657,6.64];
ybs[7473]=['',5.1584337,0.2421486,6.01];
ybs[7474]=['',5.1501052,0.803125,6.2];
ybs[7475]=['β Sge',5.1580963,0.3060368,4.37];
ybs[7476]=['55 Sgr',5.165807,-0.2803786,5.06];
ybs[7477]=['',5.1587576,0.3928976,6.36];
ybs[7478]=['',5.1716654,-0.6541295,6.16];
ybs[7479]=['',5.1551814,0.7528654,6.16];
ybs[7480]=['46 Aql',5.1633931,0.213846,6.34];
ybs[7481]=['',5.2376815,-1.4186464,6.39];
ybs[7482]=['',5.1556567,0.7955786,5.06];
ybs[7483]=['',5.1703172,-0.2689565,5.49];
ybs[7484]=['χ Aql',5.1649495,0.2074499,5.27];
ybs[7485]=['',5.2019833,-1.2643152,5.41];
ybs[7486]=['',5.1609254,0.7035913,6.23];
ybs[7487]=['',5.1513298,1.0570601,6.51];
ybs[7488]=['',5.1652921,0.5129713,6.49];
ybs[7489]=['',5.1648072,0.5669883,5.94];
ybs[7490]=['16 Cyg',5.1595089,0.8828587,5.96];
ybs[7491]=['',5.1597354,0.8827234,6.2];
ybs[7492]=['',5.1667093,0.536483,6.05];
ybs[7493]=['10 Vul',5.1693761,0.450851,5.49];
ybs[7494]=['',5.1817902,-0.5558407,5.52];
ybs[7495]=['',5.1702584,0.4746524,6.28];
ybs[7496]=['',5.1601018,0.9690458,6.48];
ybs[7497]=['ν Tel',5.1923712,-0.9826211,5.35];
ybs[7498]=['ψ Aql',5.1736276,0.2332315,6.26];
ybs[7499]=['',5.1695755,0.5972946,6.05];
ybs[7500]=['',5.2022564,-1.1649938,6.45];
ybs[7501]=['',5.1686502,0.7301225,5.84];
ybs[7502]=['56 Sgr',5.1827206,-0.3438243,4.86];
ybs[7503]=['',5.1799191,-0.0492568,6.48];
ybs[7504]=['15 Cyg',5.1712284,0.6530083,4.89];
ybs[7505]=['',5.173998,0.5118215,6.82];
ybs[7506]=['υ Aql',5.1786328,0.1339423,5.91];
ybs[7507]=['',5.1729676,0.6016895,6.57];
ybs[7508]=['',5.1959024,-0.9219727,6.25];
ybs[7509]=['',5.1649179,1.0136157,6.22];
ybs[7510]=['',5.1733785,0.711695,6.34];
ybs[7511]=['',5.2069174,-1.143904,6.05];
ybs[7512]=['γ Aql',5.1811151,0.1863076,2.72];
ybs[7513]=['',5.1668656,0.9966221,6.27];
ybs[7514]=['',5.2032133,-1.06461,6.21];
ybs[7515]=['δ Cyg',5.1737624,0.7887378,2.87];
ybs[7516]=['',5.1773423,0.6309719,6.43];
ybs[7517]=['',5.1782617,0.6121534,6.09];
ybs[7518]=['',5.2046139,-1.0319981,5.42];
ybs[7519]=['',5.1898402,-0.2380807,6.11];
ybs[7520]=['',5.1822545,0.4397422,6.62];
ybs[7521]=['17 Cyg',5.1808207,0.5897313,4.99];
ybs[7522]=['',5.1815494,0.5750866,6.18];
ybs[7523]=['δ Sge',5.1857255,0.3245624,3.82];
ybs[7524]=['',5.2011377,-0.8289205,5.94];
ybs[7525]=['',5.1954728,-0.5013621,6.05];
ybs[7526]=['',5.1872707,0.4441161,5.95];
ybs[7527]=['',5.1940348,-0.1886352,6.04];
ybs[7528]=['',5.1909168,0.1877392,6.44];
ybs[7529]=['',5.1850814,0.6714171,5.77];
ybs[7530]=['π Aql',5.1917249,0.2073178,5.72];
ybs[7531]=['',5.1673199,1.2112042,5.92];
ybs[7532]=['ζ Sge',5.1926407,0.3351896,5];
ybs[7533]=['',5.1843408,0.8372275,6.12];
ybs[7534]=['',5.2122784,-0.9582947,5.74];
ybs[7535]=['',5.2123807,-0.9583914,6.5];
ybs[7536]=['',5.1907992,0.6173917,6.53];
ybs[7537]=['',5.1913885,0.5846825,6.44];
ybs[7538]=['',5.2076278,-0.6948165,5.33];
ybs[7539]=['51 Aql',5.2016358,-0.1867479,5.39];
ybs[7540]=['',5.1988125,0.139032,6.51];
ybs[7541]=['',5.1937755,0.6767152,6.11];
ybs[7542]=['',5.1962967,0.4974843,6.38];
ybs[7543]=['α Aql',5.2009219,0.1558935,0.77];
ybs[7544]=['',5.2220049,-1.0664807,6.24];
ybs[7545]=['',5.2030933,-0.0418336,6.13];
ybs[7546]=['ο Aql',5.2019253,0.1828997,5.11];
ybs[7547]=['57 Sgr',5.2081556,-0.331272,5.92];
ybs[7548]=['',5.2031261,0.1691966,6.25];
ybs[7549]=['',5.1782399,1.1955435,6.34];
ybs[7550]=['χ Cyg',5.1989336,0.5755694,4.23];
ybs[7551]=['12 Vul',5.2016175,0.3957328,4.95];
ybs[7552]=['19 Cyg',5.198615,0.6769435,5.12];
ybs[7553]=['',5.1987386,0.7097076,5.69];
ybs[7554]=['',5.1996148,0.6613055,6.06];
ybs[7555]=['',5.2063816,0.2040858,6.13];
ybs[7556]=['η Aql',5.2085875,0.018678,3.9];
ybs[7557]=['',5.211933,-0.2537373,6.48];
ybs[7558]=['',5.2073102,0.1817911,6.54];
ybs[7559]=['',5.2056965,0.4373192,5.57];
ybs[7560]=['9 Sge',5.2074354,0.3270129,6.23];
ybs[7561]=['',5.2123961,-0.0532218,5.65];
ybs[7562]=['20 Cyg',5.1977769,0.9259238,5.03];
ybs[7563]=['',5.2013195,0.8280033,6.2];
ybs[7564]=['',5.2174992,-0.4167059,6.18];
ybs[7565]=['',5.2411518,-1.2059479,5.75];
ybs[7566]=['',5.2124048,0.07794,6.53];
ybs[7567]=['ι Sgr',5.2226554,-0.7295852,4.13];
ybs[7568]=['ε Dra',5.1839412,1.227486,3.83];
ybs[7569]=['',5.2061948,0.6369864,6.1];
ybs[7570]=['56 Aql',5.2161967,-0.1485042,5.79];
ybs[7571]=['',5.2213743,-0.5756154,6.46];
ybs[7572]=['',5.2315117,-1.0098197,6.53];
ybs[7573]=['',5.232257,-1.0268498,5.26];
ybs[7574]=['',5.2418884,-1.1989358,6.39];
ybs[7575]=['',5.2042378,0.8219056,5.62];
ybs[7576]=['ε Pav',5.2507063,-1.2713212,3.96];
ybs[7577]=['',5.2047549,0.8376922,5.91];
ybs[7578]=['13 Vul',5.2120015,0.4214062,4.58];
ybs[7579]=['57 Aql',5.2183216,-0.1424444,5.71];
ybs[7580]=['57 Aql',5.2183582,-0.1426189,6.49];
ybs[7581]=['ξ Aql',5.2160586,0.1488227,4.71];
ybs[7582]=['58 Aql',5.2185337,0.0059238,5.61];
ybs[7583]=['ω Sgr',5.2243414,-0.4578525,4.7];
ybs[7584]=['',5.217955,0.1257689,6.15];
ybs[7585]=['',5.2213202,-0.1163794,6.51];
ybs[7586]=['',5.208697,0.8355286,6.29];
ybs[7587]=['',5.216604,0.4255998,5.52];
ybs[7588]=['β Aql',5.2207823,0.1129707,3.71];
ybs[7589]=['μ1 Pav',5.2479572,-1.1672839,5.76];
ybs[7590]=['59 Sgr',5.2292069,-0.4730365,4.52];
ybs[7591]=['',5.2329984,-0.6630664,6.55];
ybs[7592]=['',5.2172259,0.6468512,5.76];
ybs[7593]=['',5.2189085,0.5281575,6.57];
ybs[7594]=['23 Cyg',5.2088943,1.0051081,5.14];
ybs[7595]=['10 Sge',5.2234855,0.2914901,5.36];
ybs[7596]=['φ Aql',5.2246323,0.2005461,5.28];
ybs[7597]=['',5.2099102,1.0432509,6.06];
ybs[7598]=['μ2 Pav',5.2544388,-1.1671785,5.31];
ybs[7599]=['22 Cyg',5.2217538,0.672876,4.94];
ybs[7600]=['61 Sgr',5.233094,-0.2691976,5.02];
ybs[7601]=['η Cyg',5.22389,0.6134808,3.89];
ybs[7602]=['',5.2275938,0.3676533,6.48];
ybs[7603]=['',5.2380571,-0.5318019,6.28];
ybs[7604]=['60 Sgr',5.2379133,-0.4560111,4.83];
ybs[7605]=['ψ Cyg',5.219688,0.9163849,4.92];
ybs[7606]=['',5.2257037,0.6338613,6.02];
ybs[7607]=['',5.2456969,-0.8601372,6.17];
ybs[7608]=['11 Sge',5.2310525,0.2942013,5.53];
ybs[7609]=['θ1 Sgr',5.2417724,-0.614494,4.37];
ybs[7610]=['θ2 Sgr',5.2422574,-0.6043944,5.3];
ybs[7611]=['',5.2525342,-1.0350935,5.13];
ybs[7612]=['',5.2178851,1.0178097,6.09];
ybs[7613]=['',5.2453108,-0.750046,6.14];
ybs[7614]=['',5.2276227,0.7057247,5.45];
ybs[7615]=['',5.2442103,-0.656828,5.95];
ybs[7616]=['',5.2470378,-0.7861661,5.81];
ybs[7617]=['',5.2443125,-0.5870434,5.66];
ybs[7618]=['',5.2247347,0.8895802,6.43];
ybs[7619]=['',5.2202315,1.0282139,4.96];
ybs[7620]=['',5.2222276,0.9905332,6.12];
ybs[7621]=['γ Sge',5.2353093,0.3413875,3.47];
ybs[7622]=['',5.2387034,0.0252371,6.17];
ybs[7623]=['',5.2409108,-0.1726113,5.88];
ybs[7624]=['',5.2305805,0.7387661,6.43];
ybs[7625]=['',5.2495179,-0.7111309,6.29];
ybs[7626]=['',5.2342612,0.5419484,5.49];
ybs[7627]=['14 Vul',5.2369855,0.4043829,5.67];
ybs[7628]=['',5.2336189,0.6662488,6.32];
ybs[7629]=['',5.248426,-0.39563,6.01];
ybs[7630]=['',5.270508,-1.1737199,6.07];
ybs[7631]=['13 Sge',5.2410557,0.3069193,5.37];
ybs[7632]=['',5.2364723,0.8000638,5.92];
ybs[7633]=['25 Cyg',5.2395613,0.647712,5.19];
ybs[7634]=['',5.2454399,0.1505706,5.91];
ybs[7635]=['63 Sgr',5.2505871,-0.2367955,5.71];
ybs[7636]=['62 Sgr',5.2541378,-0.4824051,4.58];
ybs[7637]=['',5.2355733,0.9097322,6.15];
ybs[7638]=['',5.2585802,-0.6609631,4.77];
ybs[7639]=['15 Vul',5.2451917,0.4855966,4.64];
ybs[7640]=['',5.2306657,1.1100583,5.96];
ybs[7641]=['',5.2453884,0.6487034,6.2];
ybs[7642]=['',5.2481322,0.4340571,5.88];
ybs[7643]=['16 Vul',5.2493413,0.4364643,5.22];
ybs[7644]=['',5.258625,-0.3931366,6.45];
ybs[7645]=['',5.2616179,-0.5582536,4.99];
ybs[7646]=['26 Cyg',5.2449425,0.875698,5.05];
ybs[7647]=['',5.2592731,-0.1291391,6.72];
ybs[7648]=['',5.2550649,0.3241201,5.96];
ybs[7649]=['',5.2824483,-1.1568337,6.45];
ybs[7650]=['',5.2561508,0.2810271,5.67];
ybs[7651]=['δ Pav',5.2840861,-1.1538148,3.56];
ybs[7652]=['',5.2299778,1.2293131,6.33];
ybs[7653]=['62 Aql',5.2606312,-0.011147,5.68];
ybs[7654]=['',5.266933,-0.5747121,6.53];
ybs[7655]=['τ Aql',5.2592623,0.1282586,5.52];
ybs[7656]=['',5.256102,0.5230224,5.71];
ybs[7657]=['',5.2640917,-0.2012067,6.34];
ybs[7658]=['15 Sge',5.2587441,0.2991597,5.8];
ybs[7659]=['ξ Tel',5.2763589,-0.9216804,4.94];
ybs[7660]=['',5.2774394,-0.9589508,6.26];
ybs[7661]=['65 Sgr',5.2656579,-0.2198059,6.55];
ybs[7662]=['64 Dra',5.2435283,1.1325465,5.27];
ybs[7663]=['',5.2623086,0.4063353,6.45];
ybs[7664]=['',5.2602613,0.5635565,5.64];
ybs[7665]=['η Sge',5.2632338,0.3501519,5.1];
ybs[7666]=['',5.2646463,0.2717748,6.34];
ybs[7667]=['',5.2686771,-0.069929,6.47];
ybs[7668]=['65 Dra',5.2473417,1.1292963,6.57];
ybs[7669]=['',5.262379,0.672814,6.19];
ybs[7670]=['',5.2587321,0.8430011,6.16];
ybs[7671]=['ρ Dra',5.2487219,1.185834,4.51];
ybs[7672]=['69 Dra',5.2311514,1.3360357,6.2];
ybs[7673]=['',5.2611804,0.9060079,6.14];
ybs[7674]=['17 Vul',5.2706446,0.4134062,5.07];
ybs[7675]=['27 Cyg',5.2677601,0.6290897,5.36];
ybs[7676]=['64 Aql',5.2765268,-0.0105719,5.99];
ybs[7677]=['',5.2931183,-1.0026838,6.37];
ybs[7678]=['',5.2617916,0.9845833,6.21];
ybs[7679]=['',5.2753328,0.1653214,6.43];
ybs[7680]=['',5.2790053,-0.1743566,6.18];
ybs[7681]=['',5.2580033,1.1163346,6.26];
ybs[7682]=['',5.2762414,0.292112,6.42];
ybs[7683]=['',5.2659409,0.9291675,5.85];
ybs[7684]=['',5.366286,-1.4526159,6.17];
ybs[7685]=['',5.2736359,0.6020586,6.11];
ybs[7686]=['',5.2787762,0.1884734,6.31];
ybs[7687]=['66 Dra',5.2618436,1.0832681,5.39];
ybs[7688]=['',5.2704146,0.8779218,6.54];
ybs[7689]=['',5.2917827,-0.6287868,5.32];
ybs[7690]=['',5.2577421,1.1885334,6.28];
ybs[7691]=['θ Sge',5.284085,0.3663238,6.48];
ybs[7692]=['',5.2974034,-0.7453543,6.22];
ybs[7693]=['',5.3084524,-1.105488,6.09];
ybs[7694]=['28 Cyg',5.2810917,0.6442525,4.93];
ybs[7695]=['',5.2905072,-0.1530308,6.49];
ybs[7696]=['θ Aql',5.2908215,-0.0130399,3.23];
ybs[7697]=['18 Vul',5.2865189,0.4708547,5.52];
ybs[7698]=['ξ1 Cap',5.2941105,-0.2149877,6.34];
ybs[7699]=['',5.2889437,0.3701635,6.22];
ybs[7700]=['',5.3065042,-0.914023,5.65];
ybs[7701]=['ξ2 Cap',5.2961543,-0.2189107,5.85];
ybs[7702]=['',5.2901944,0.3830961,6.26];
ybs[7703]=['',5.296355,0.0164478,6.27];
ybs[7704]=['19 Vul',5.2919444,0.4692027,5.49];
ybs[7705]=['20 Vul',5.2928829,0.4634449,5.92];
ybs[7706]=['66 Aql',5.2992345,-0.0163054,5.47];
ybs[7707]=['',5.2919059,0.8344673,6.92];
ybs[7708]=['',5.3091873,-0.4704795,5.73];
ybs[7709]=['',5.3002477,0.4243641,6.56];
ybs[7710]=['ρ Aql',5.3032212,0.2665676,4.95];
ybs[7711]=['',5.3117323,-0.5223545,6.3];
ybs[7712]=['',5.2936354,0.8995136,6.01];
ybs[7713]=['68 Dra',5.2882257,1.0847708,5.75];
ybs[7714]=['',5.3144542,-0.6349087,6.39];
ybs[7715]=['',5.3145922,-0.6129753,6.53];
ybs[7716]=['30 Cyg',5.2973792,0.8184016,4.83];
ybs[7717]=['21 Vul',5.3025218,0.5021378,5.18];
ybs[7718]=['',5.3284458,-1.1022246,6.27];
ybs[7719]=['',5.2994249,0.7584243,6.14];
ybs[7720]=['',5.3014426,0.6402012,6.45];
ybs[7721]=['31 Cyg',5.2988338,0.8171052,3.79];
ybs[7722]=['29 Cyg',5.3033972,0.6437152,4.97];
ybs[7723]=['',5.3023319,0.7361673,6.71];
ybs[7724]=['3 Cap',5.3133669,-0.21398,6.32];
ybs[7725]=['',5.3071255,0.447993,4.78];
ybs[7726]=['33 Cyg',5.296914,0.9886052,4.3];
ybs[7727]=['22 Vul',5.3082529,0.4116342,5.15];
ybs[7728]=['',5.2966669,1.0596886,5.79];
ybs[7729]=['',5.30732,0.59002,5.66];
ybs[7730]=['23 Vul',5.3092271,0.4867824,4.52];
ybs[7731]=['',5.3262213,-0.8313473,6.31];
ybs[7732]=['18 Sge',5.3119318,0.3783056,6.13];
ybs[7733]=['α1 Cap',5.3189056,-0.2169604,4.24];
ybs[7734]=['4 Cap',5.3208902,-0.3793016,5.87];
ybs[7735]=['',5.3278017,-0.8290656,6.13];
ybs[7736]=['κ Cep',5.2709645,1.3575851,4.39];
ybs[7737]=['32 Cyg',5.3068007,0.8341041,3.98];
ybs[7738]=['',5.3099351,0.6802298,6.27];
ybs[7739]=['24 Vul',5.3138005,0.4319346,5.32];
ybs[7740]=['α2 Cap',5.3206808,-0.2175921,3.57];
ybs[7741]=['',5.3076985,0.8780593,6.31];
ybs[7742]=['',5.3093094,0.7968461,5.91];
ybs[7743]=['',5.3118425,0.6480952,6.48];
ybs[7744]=['',5.3338556,-0.9594403,6.27];
ybs[7745]=['',5.3136284,0.7058451,5.24];
ybs[7746]=['',5.3168382,0.5100781,6.22];
ybs[7747]=['σ Cap',5.3267604,-0.3323164,5.28];
ybs[7748]=['',5.3159448,0.7469862,6.29];
ybs[7749]=['34 Cyg',5.3175401,0.6651524,4.81];
ybs[7750]=['',5.3318606,-0.5082119,6.3];
ybs[7751]=['',5.3338983,-0.6212425,6.46];
ybs[7752]=['',5.3383804,-0.8712676,6.27];
ybs[7753]=['',5.3188267,0.7122645,5.84];
ybs[7754]=['',5.3275542,-0.0174568,6.06];
ybs[7755]=['36 Cyg',5.3206101,0.6471282,5.58];
ybs[7756]=['35 Cyg',5.3214798,0.6119227,5.17];
ybs[7757]=['',5.3260341,0.2320452,6.21];
ybs[7758]=['',5.3308467,-0.1096573,6.63];
ybs[7759]=['ν Cap',5.332066,-0.2213127,4.76];
ybs[7760]=['',5.3282774,0.2378286,5.95];
ybs[7761]=['',5.3326331,-0.256669,6.1];
ybs[7762]=['β Cap',5.3336579,-0.256604,3.08];
ybs[7763]=['',5.3215637,0.8098438,6.45];
ybs[7764]=['',5.3297172,0.2556531,6.13];
ybs[7765]=['κ1 Sgr',5.3411984,-0.7325128,5.59];
ybs[7766]=['',5.329664,0.3119206,5.8];
ybs[7767]=['',5.3189685,0.9682184,5.76];
ybs[7768]=['',5.326405,0.6494521,6.57];
ybs[7769]=['',5.3133163,1.1681656,5.93];
ybs[7770]=['',5.3282378,0.6890893,6.23];
ybs[7771]=['',5.3983389,-1.4116137,5.77];
ybs[7772]=['',5.3263695,0.8188368,6.5];
ybs[7773]=['κ2 Sgr',5.3474456,-0.7390119,5.64];
ybs[7774]=['',5.342208,-0.1671101,6.3];
ybs[7775]=['25 Vul',5.3368261,0.4280524,5.54];
ybs[7776]=['α Pav',5.3562856,-0.9887915,1.94];
ybs[7777]=['',5.3283024,0.9368011,6.18];
ybs[7778]=['71 Dra',5.3233229,1.0879623,5.72];
ybs[7779]=['',5.340761,0.2553642,6.17];
ybs[7780]=['',5.3424197,0.0946514,5.31];
ybs[7781]=['',5.3359857,0.7192647,6.39];
ybs[7782]=['γ Cyg',5.3368163,0.7039995,2.2];
ybs[7783]=['',5.3390053,0.547069,6.09];
ybs[7784]=['',5.3358603,0.80066,5.58];
ybs[7785]=['',5.3556792,-0.7106101,6.09];
ybs[7786]=['',5.3390718,0.7174332,5.93];
ybs[7787]=['',5.3535425,-0.4988519,5.85];
ybs[7788]=['',5.3396948,0.7515945,6.2];
ybs[7789]=['',5.3488765,0.0200605,6.15];
ybs[7790]=['',5.3241131,1.2035543,5.55];
ybs[7791]=['',5.3299406,1.1180435,5.69];
ybs[7792]=['39 Cyg',5.34437,0.5632237,4.43];
ybs[7793]=['',5.3435743,0.6554873,5.9];
ybs[7794]=['',5.3602206,-0.6513765,6.25];
ybs[7795]=['',5.3537327,-0.0474552,6.11];
ybs[7796]=['',5.3534092,0.1769358,6.33];
ybs[7797]=['',5.3527438,0.3750879,5.66];
ybs[7798]=['',5.4204545,-1.4172262,5.91];
ybs[7799]=['',5.3543088,0.3481305,6.41];
ybs[7800]=['π Cap',5.3612967,-0.3164207,5.25];
ybs[7801]=['',5.3459384,0.9360644,6.51];
ybs[7802]=['',5.3560041,0.3036374,6.22];
ybs[7803]=['',5.3683878,-0.6198189,6.1];
ybs[7804]=['',5.3476467,1.0416267,6.44];
ybs[7805]=['',5.3673425,-0.2732995,6.41];
ybs[7806]=['',5.3639025,0.1487008,6.25];
ybs[7807]=['68 Aql',5.3655604,-0.0571629,6.13];
ybs[7808]=['ρ Cap',5.3679986,-0.3094606,4.78];
ybs[7809]=['',5.3585344,0.6005815,6.39];
ybs[7810]=['',5.3647629,0.0526995,6.21];
ybs[7811]=['',5.3710552,-0.3893571,6.16];
ybs[7812]=['40 Cyg',5.3602597,0.6723423,5.62];
ybs[7813]=['',5.3537176,0.9899569,6.36];
ybs[7814]=['43 Cyg',5.3572194,0.8633296,5.69];
ybs[7815]=['ο Cap',5.3724534,-0.3229446,6.74];
ybs[7816]=['ο Cap',5.372555,-0.3228863,5.94];
ybs[7817]=['69 Aql',5.3709331,-0.0489108,4.91];
ybs[7818]=['',5.3775411,-0.5066459,6.39];
ybs[7819]=['',5.3688452,0.3520463,6.55];
ybs[7820]=['41 Cyg',5.3686177,0.5314806,4.01];
ybs[7821]=['42 Cyg',5.3680885,0.6377025,5.88];
ybs[7822]=['1 Del',5.3733115,0.1916247,6.08];
ybs[7823]=['',5.3775461,-0.26132,6.12];
ybs[7824]=['',5.4028904,-1.2134364,6.11];
ybs[7825]=['',5.3758836,0.3611013,6.18];
ybs[7826]=['',5.3773,0.1979978,7.11];
ybs[7827]=['',5.3704082,0.803058,6.41];
ybs[7828]=['',5.3857698,-0.4338743,6.36];
ybs[7829]=['',5.36716,0.9800194,5.91];
ybs[7830]=['ω1 Cyg',5.3704521,0.8558206,4.95];
ybs[7831]=['',5.3831399,-0.1704988,5.65];
ybs[7832]=['ν Mic',5.3912878,-0.7754646,5.11];
ybs[7833]=['44 Cyg',5.3752565,0.646113,6.19];
ybs[7834]=['φ1 Pav',5.3999822,-1.0558467,4.76];
ybs[7835]=['',5.380055,0.4518422,6.34];
ybs[7836]=['θ Cep',5.3668261,1.100903,4.22];
ybs[7837]=['ω2 Cyg',5.3759166,0.860519,5.44];
ybs[7838]=['ε Del',5.3860132,0.1987607,4.03];
ybs[7839]=['',5.3953051,-0.6632952,6.44];
ybs[7840]=['',5.3758301,0.9144401,6.18];
ybs[7841]=['',5.3911227,-0.2379896,6.13];
ybs[7842]=['',5.3943553,-0.5303703,6.4];
ybs[7843]=['',5.3890229,0.1770612,6.56];
ybs[7844]=['η Del',5.3891699,0.2288541,5.38];
ybs[7845]=['ρ Pav',5.4088508,-1.072382,4.88];
ybs[7846]=['',5.3772429,0.9924641,6.14];
ybs[7847]=['',5.3831407,0.755313,6.6];
ybs[7848]=['',5.3898326,0.3677499,6.48];
ybs[7849]=['μ1 Oct',5.4322302,-1.3280435,6];
ybs[7850]=['μ2 Oct',5.4304146,-1.3135602,6.55];
ybs[7851]=['',5.3970717,-0.2869302,6.19];
ybs[7852]=['47 Cyg',5.388081,0.6167281,4.61];
ybs[7853]=['',5.387311,0.7305467,6.49];
ybs[7854]=['',5.3663996,1.2673652,6.27];
ybs[7855]=['α Ind',5.4073663,-0.8238728,3.11];
ybs[7856]=['',5.3874716,0.8164466,5.78];
ybs[7857]=['ζ Del',5.3950423,0.2576097,4.68];
ybs[7858]=['',5.4189208,-1.0964105,6.22];
ybs[7859]=['70 Aql',5.4017977,-0.0429968,4.89];
ybs[7860]=['26 Vul',5.3982528,0.4532383,6.41];
ybs[7861]=['φ2 Pav',5.4194004,-1.055239,5.12];
ybs[7862]=['',5.391119,0.9065176,6.11];
ybs[7863]=['',5.4075558,-0.4367138,6.36];
ybs[7864]=['',5.404236,0.0032057,6.22];
ybs[7865]=['73 Dra',5.3719252,1.3096669,5.2];
ybs[7866]=['27 Vul',5.4023342,0.4633591,5.59];
ybs[7867]=['υ Pav',5.4287053,-1.163642,5.15];
ybs[7868]=['β Del',5.4048299,0.2562509,3.63];
ybs[7869]=['ι Del',5.4061099,0.2000972,5.43];
ybs[7870]=['71 Aql',5.4087833,-0.0177687,4.32];
ybs[7871]=['48 Cyg',5.4040917,0.5525586,6.32];
ybs[7872]=['',5.4062734,0.3203751,6.25];
ybs[7873]=['',5.4041522,0.5516763,6.49];
ybs[7874]=['',5.4031708,0.6704736,6.2];
ybs[7875]=['τ Cap',5.4133147,-0.2594792,5.22];
ybs[7876]=['',5.4126718,-0.0405818,6.22];
ybs[7877]=['29 Vul',5.4088378,0.371552,4.82];
ybs[7878]=['θ Del',5.4100314,0.2339153,5.72];
ybs[7879]=['',5.4186367,-0.5819586,5.47];
ybs[7880]=['28 Vul',5.4087721,0.4224283,5.04];
ybs[7881]=['',5.4090216,0.4148269,5.91];
ybs[7882]=['κ Del',5.4118775,0.1775638,5.05];
ybs[7883]=['1 Aqr',5.4134303,0.0100196,5.16];
ybs[7884]=['',5.4176495,-0.4133951,6.37];
ybs[7885]=['',5.4114805,0.2779536,5.97];
ybs[7886]=['υ Cap',5.4168121,-0.3150424,5.1];
ybs[7887]=['75 Dra',5.3520108,1.4225235,5.46];
ybs[7888]=['',5.4195265,-0.4635022,6.51];
ybs[7889]=['',5.4116748,0.3823104,6.08];
ybs[7890]=['',5.4105291,0.5309621,5.68];
ybs[7891]=['',5.4188864,-0.2798798,5.8];
ybs[7892]=['α Del',5.4139006,0.2792477,3.77];
ybs[7893]=['',5.4150398,0.1978785,6.42];
ybs[7894]=['74 Dra',5.3578445,1.4167501,5.96];
ybs[7895]=['',5.4231976,-0.5499478,5.76];
ybs[7896]=['',5.4229964,-0.4522385,6.28];
ybs[7897]=['',5.4124937,0.7097751,6.06];
ybs[7898]=['',5.4114454,0.7985671,6.58];
ybs[7899]=['β Pav',5.4416682,-1.1538835,3.42];
ybs[7900]=['',5.4186269,0.3494768,6.45];
ybs[7901]=['',5.4300642,-0.6888687,6.29];
ybs[7902]=['',5.4089356,0.9789964,6.48];
ybs[7903]=['',5.41759,0.5217392,6.08];
ybs[7904]=['10 Del',5.4210684,0.2560671,5.99];
ybs[7905]=['',5.4144975,0.7600299,5.95];
ybs[7906]=['η Ind',5.4359015,-0.9046252,4.51];
ybs[7907]=['49 Cyg',5.4194035,0.5654098,5.51];
ybs[7908]=['',5.4189167,0.6836553,6.51];
ybs[7909]=['',5.4240235,0.3073562,6.22];
ybs[7910]=['α Cyg',5.4204147,0.7918348,1.25];
ybs[7911]=['',5.4140783,1.0575509,6.01];
ybs[7912]=['',5.4228578,0.7296471,5.67];
ybs[7913]=['',5.4250537,0.6203787,6.66];
ybs[7914]=['δ Del',5.4306042,0.2646607,4.43];
ybs[7915]=['51 Cyg',5.4234605,0.8801494,5.39];
ybs[7916]=['',5.3510831,1.4609717,6.19];
ybs[7917]=['',5.4396928,-0.4739766,6.5];
ybs[7918]=['',5.4295611,0.6226848,6.47];
ybs[7919]=['',5.4450867,-0.6825683,5.5];
ybs[7920]=['σ Pav',5.4610866,-1.1987627,5.41];
ybs[7921]=['',5.4448258,-0.6288319,6.49];
ybs[7922]=['ψ Cap',5.4434314,-0.4394754,4.14];
ybs[7923]=['17 Cap',5.4436057,-0.3739088,5.93];
ybs[7924]=['',5.4244059,1.0592471,6.15];
ybs[7925]=['30 Vul',5.4364256,0.4426272,4.91];
ybs[7926]=['',5.4272831,0.9983884,6.32];
ybs[7927]=['',5.4392832,0.3173126,6.38];
ybs[7928]=['52 Cyg',5.4396427,0.537739,4.22];
ybs[7929]=['ι Mic',5.454726,-0.7661431,5.11];
ybs[7930]=['',5.4323422,0.9874698,5.78];
ybs[7931]=['4 Cep',5.4257069,1.1649494,5.58];
ybs[7932]=['',5.4468802,-0.0418099,6.27];
ybs[7933]=['γ1 Del',5.4444795,0.2830114,5.14];
ybs[7934]=['γ2 Del',5.4445376,0.2830067,4.27];
ybs[7935]=['ε Cyg',5.4419032,0.5944761,2.46];
ybs[7936]=['ε Aqr',5.4497896,-0.1641381,3.77];
ybs[7937]=['3 Aqr',5.4499113,-0.0861554,4.42];
ybs[7938]=['ζ Ind',5.4592052,-0.8052016,4.89];
ybs[7939]=['13 Del',5.4498681,0.1064613,5.58];
ybs[7940]=['',5.4499157,0.0593083,6.4];
ybs[7941]=['',5.4365348,1.0065305,4.51];
ybs[7942]=['',5.4461145,0.6015328,4.92];
ybs[7943]=['η Cep',5.4357388,1.080866,3.43];
ybs[7944]=['',5.4431032,0.8137164,6.3];
ybs[7945]=['',5.4701006,-1.087966,6.28];
ybs[7946]=['',5.4701369,-1.087966,6.59];
ybs[7947]=['',5.4573955,-0.4483617,5.86];
ybs[7948]=['',5.4413583,0.926525,6.33];
ybs[7949]=['λ Cyg',5.4470132,0.6384772,4.53];
ybs[7950]=['',5.4573303,-0.3131761,6.21];
ybs[7951]=['α Mic',5.4606633,-0.5879535,4.9];
ybs[7952]=['',5.4462458,0.7971073,6.4];
ybs[7953]=['',5.4309893,1.2189678,6.41];
ybs[7954]=['ι Ind',5.4683755,-0.8991089,5.05];
ybs[7955]=['',5.4481735,0.8364196,5.57];
ybs[7956]=['',5.4641484,-0.5578356,6.36];
ybs[7957]=['',5.465402,-0.6600905,5.52];
ybs[7958]=['',5.4480833,0.9162734,6.27];
ybs[7959]=['15 Del',5.4576268,0.2205663,5.98];
ybs[7960]=['14 Del',5.4585387,0.138867,6.33];
ybs[7961]=['',5.4593961,0.0983865,6.21];
ybs[7962]=['',5.4630597,-0.2173328,5.88];
ybs[7963]=['55 Cyg',5.4531681,0.8064474,4.84];
ybs[7964]=['',5.4517603,0.9076114,6.29];
ybs[7965]=['β Mic',5.4694052,-0.5774227,6.04];
ybs[7966]=['ω Cap',5.4684617,-0.4682004,4.11];
ybs[7967]=['',5.4617527,0.3166735,6.52];
ybs[7968]=['4 Aqr',5.4660353,-0.0965749,5.99];
ybs[7969]=['',5.457387,0.8160008,6.33];
ybs[7970]=['56 Cyg',5.4582916,0.7705946,5.04];
ybs[7971]=['5 Aqr',5.4691581,-0.0944847,5.55];
ybs[7972]=['β Ind',5.4834077,-1.0185658,3.65];
ybs[7973]=['',5.4770804,-0.6931734,5.35];
ybs[7974]=['',5.4651097,0.4946887,5.77];
ybs[7975]=['',5.4735703,-0.4134559,6.33];
ybs[7976]=['μ Aqr',5.4714826,-0.1551552,4.73];
ybs[7977]=['',5.4755644,-0.534501,6.35];
ybs[7978]=['',5.4817271,-0.8837171,6.24];
ybs[7979]=['',5.4529152,1.1193521,6.45];
ybs[7980]=['',5.4734826,-0.2003607,6.38];
ybs[7981]=['31 Vul',5.4680293,0.4745594,4.59];
ybs[7982]=['',5.4672576,0.5749535,6.44];
ybs[7983]=['',5.478491,-0.4857479,6.41];
ybs[7984]=['',5.4771765,-0.1186049,6.44];
ybs[7985]=['',5.4722747,0.5191164,6.34];
ybs[7986]=['19 Cap',5.4811231,-0.3111664,5.78];
ybs[7987]=['57 Cyg',5.472104,0.7763395,4.78];
ybs[7988]=['76 Dra',5.4133754,1.4419827,5.75];
ybs[7989]=['',5.472333,0.7902104,5.45];
ybs[7990]=['',5.473062,0.7418368,6.66];
ybs[7991]=['',5.4755022,0.585241,5.47];
ybs[7992]=['',5.4820783,-0.0223172,6.55];
ybs[7993]=['',5.4777784,0.4994472,6.56];
ybs[7994]=['32 Vul',5.4786117,0.4913426,5.01];
ybs[7995]=['',5.4772251,0.7120472,6.7];
ybs[7996]=['',5.4842723,0.0807678,6.05];
ybs[7997]=['17 Del',5.4836991,0.2411385,5.17];
ybs[7998]=['16 Del',5.4838737,0.221019,5.58];
ybs[7999]=['',5.4900917,-0.4572936,5.7];
ybs[8000]=['',5.4872485,-0.0604973,6.57];
ybs[8001]=['7 Aqr',5.4900291,-0.167588,5.51];
ybs[8002]=['',5.4382445,1.4074864,5.39];
ybs[8003]=['',5.4909319,0.0097586,6.05];
ybs[8004]=['',5.4936143,-0.2781342,5.87];
ybs[8005]=['',5.5138259,-1.1887821,6.37];
ybs[8006]=['',5.4832065,0.8292518,5.67];
ybs[8007]=['α Oct',5.5309249,-1.3425926,5.15];
ybs[8008]=['',5.4846324,0.8930851,6.63];
ybs[8009]=['',5.4866328,0.7857504,5.96];
ybs[8010]=['',5.4980211,-0.2510935,6.01];
ybs[8011]=['',5.4855766,0.8870411,5.81];
ybs[8012]=['',5.4857165,0.8602892,5.9];
ybs[8013]=['',5.5069134,-0.8930556,5.76];
ybs[8014]=['ν Cyg',5.4894429,0.7201694,3.94];
ybs[8015]=['',5.4844255,0.9945324,6.23];
ybs[8016]=['18 Del',5.4960979,0.1908552,5.48];
ybs[8017]=['',5.5044669,-0.6288937,6.11];
ybs[8018]=['33 Vul',5.4950357,0.3913339,5.31];
ybs[8019]=['20 Cap',5.5021066,-0.3305428,6.25];
ybs[8020]=['ε Equ',5.4990949,0.0766187,5.23];
ybs[8021]=['',5.4942934,0.7778512,5.55];
ybs[8022]=['',5.4952575,0.7336716,6.16];
ybs[8023]=['',5.5020875,0.2953235,6.66];
ybs[8024]=['',5.5033269,0.1328739,5.99];
ybs[8025]=['γ Mic',5.5099531,-0.5613062,4.67];
ybs[8026]=['',5.4946758,0.8824071,5.61];
ybs[8027]=['11 Aqr',5.5058613,-0.0808667,6.21];
ybs[8028]=['',5.5144555,-0.7488202,6.64];
ybs[8029]=['',5.4734125,1.3267939,6.05];
ybs[8030]=['',5.5046911,0.3390532,5.65];
ybs[8031]=['',5.5117661,-0.4674625,6.05];
ybs[8032]=['',5.5152899,-0.6707782,5.94];
ybs[8033]=['59 Cyg',5.5006703,0.831085,4.74];
ybs[8034]=['ζ Mic',5.5175316,-0.6725392,5.3];
ybs[8035]=['',5.4978904,1.0390808,5.51];
ybs[8036]=['',5.5179659,-0.482302,6.25];
ybs[8037]=['',5.5073504,0.6304701,5.97];
ybs[8038]=['',5.5481799,-1.3284027,6.58];
ybs[8039]=['60 Cyg',5.5066789,0.8072664,5.37];
ybs[8040]=['',5.5163524,-0.0144297,6.5];
ybs[8041]=['μ Ind',5.5284265,-0.9534425,5.16];
ybs[8042]=['',5.5165286,0.0284476,6.25];
ybs[8043]=['',5.51605,0.2587966,6.31];
ybs[8044]=['12 Aqr',5.5212207,-0.0999186,7.31];
ybs[8045]=['12 Aqr',5.5212279,-0.0999138,5.89];
ybs[8046]=['η Cap',5.5230792,-0.3448147,4.84];
ybs[8047]=['',5.5494436,-1.2753513,5.68];
ybs[8048]=['',5.5120924,0.7834565,6.19];
ybs[8049]=['',5.5153774,0.6764099,6.07];
ybs[8050]=['',5.5138076,0.8019211,6.48];
ybs[8051]=['',5.5101182,0.9907747,5.83];
ybs[8052]=['3 Equ',5.523083,0.0977626,5.61];
ybs[8053]=['',5.523666,0.0530685,6.42];
ybs[8054]=['',5.5239545,0.0413365,6.33];
ybs[8055]=['η Mic',5.5327339,-0.7205884,5.53];
ybs[8056]=['δ Mic',5.530479,-0.5240482,5.68];
ybs[8057]=['',5.5186749,0.7282614,6.33];
ybs[8058]=['',5.5162464,0.8805184,6.37];
ybs[8059]=['',5.5438522,-1.1140176,5.76];
ybs[8060]=['',5.5177384,0.8196087,6.32];
ybs[8061]=['θ Cap',5.5297187,-0.2990374,4.07];
ybs[8062]=['',5.5322724,-0.5627335,5.18];
ybs[8063]=['4 Equ',5.526852,0.1057197,5.94];
ybs[8064]=['',5.517598,0.9317318,5.9];
ybs[8065]=['ξ Cyg',5.5231989,0.7684068,3.72];
ybs[8066]=['24 Cap',5.5351255,-0.4346942,4.5];
ybs[8067]=['',5.5578046,-1.2643617,6.2];
ybs[8068]=['',5.5303247,0.4716536,6.12];
ybs[8069]=['',5.5375608,-0.302908,6.17];
ybs[8070]=['',5.5306605,0.5460102,5.82];
ybs[8071]=['61 Cyg',5.5321078,0.6779791,5.21];
ybs[8072]=['61 Cyg',5.5321589,0.6779355,6.03];
ybs[8073]=['χ Cap',5.5412404,-0.3681484,5.3];
ybs[8074]=['',5.5357966,0.2750363,6.34];
ybs[8075]=['63 Cyg',5.5302791,0.8333544,4.55];
ybs[8076]=['',5.5400241,0.1237374,6.15];
ybs[8077]=['27 Cap',5.5455367,-0.3570198,6.25];
ybs[8078]=['ο Pav',5.5658218,-1.2221502,5.02];
ybs[8079]=['ν Aqr',5.5454538,-0.196716,4.51];
ybs[8080]=['',5.5400598,0.5289405,5.59];
ybs[8081]=['',5.5466919,0.0531302,6.45];
ybs[8082]=['',5.5505717,-0.161491,6.27];
ybs[8083]=['γ Equ',5.5481027,0.1785926,4.69];
ybs[8084]=['6 Equ',5.5488834,0.1771492,6.07];
ybs[8085]=['',5.5262585,1.2484527,5.87];
ybs[8086]=['',5.5579554,-0.7010583,5.83];
ybs[8087]=['',5.5485636,0.3936715,6.68];
ybs[8088]=['',5.554675,-0.2508163,6.48];
ybs[8089]=['',5.5451613,0.7959265,6.63];
ybs[8090]=['',5.5615104,-0.6863191,5.26];
ybs[8091]=['',5.5503716,0.635306,6.54];
ybs[8092]=['',5.5458579,0.9366161,5.73];
ybs[8093]=['',5.5473917,0.8341432,6.46];
ybs[8094]=['',5.5625252,-0.6339333,5.96];
ybs[8095]=['',5.5414297,1.1064691,6.54];
ybs[8096]=['',5.5620709,-0.4802672,5.42];
ybs[8097]=['',5.5888205,-1.3132267,6.63];
ybs[8098]=['',5.519358,1.3652839,5.91];
ybs[8099]=['',5.5407658,1.1971339,7.33];
ybs[8100]=['',5.5741359,-0.9278148,5.75];
ybs[8101]=['ζ Cyg',5.5587896,0.5293389,3.2];
ybs[8102]=['',5.5616305,0.2807305,6.27];
ybs[8103]=['',5.5711166,-0.7051731,6.21];
ybs[8104]=['',5.5658597,-0.1833076,6.77];
ybs[8105]=['',5.5520063,1.0487299,5.64];
ybs[8106]=['',5.5607264,0.6411596,6.05];
ybs[8107]=['',5.5670395,0.0034011,6.38];
ybs[8108]=['',5.5696828,-0.3009322,6.04];
ybs[8109]=['δ Equ',5.566176,0.1764444,4.49];
ybs[8110]=['',5.5732577,-0.630198,6.12];
ybs[8111]=['',5.5850714,-1.127091,6.31];
ybs[8112]=['',5.5641959,0.5236605,6.17];
ybs[8113]=['φ Cap',5.5720583,-0.3586407,5.24];
ybs[8114]=['29 Cap',5.5724023,-0.2629912,5.28];
ybs[8115]=['',5.6585817,-1.4782972,6.45];
ybs[8116]=['τ Cyg',5.5665781,0.6658119,3.72];
ybs[8117]=['α Equ',5.572164,0.0933906,3.92];
ybs[8118]=['',5.5760013,-0.0262556,6.48];
ybs[8119]=['',5.5597862,1.1258426,6.39];
ybs[8120]=['',5.578788,-0.2299508,6.4];
ybs[8121]=['ε Mic',5.5825042,-0.5597013,4.71];
ybs[8122]=['',5.5696565,0.8390944,6.46];
ybs[8123]=['30 Cap',5.5821084,-0.3120878,5.43];
ybs[8124]=['',5.5738175,0.739229,6.43];
ybs[8125]=['31 Cap',5.5834302,-0.3029567,7.05];
ybs[8126]=['θ Ind',5.5920032,-0.9310451,4.39];
ybs[8127]=['15 Aqr',5.5827335,-0.0770634,5.82];
ybs[8128]=['',5.5866019,-0.5002325,6.4];
ybs[8129]=['σ Cyg',5.577988,0.6893771,4.23];
ybs[8130]=['',5.5777003,0.7467739,6.19];
ybs[8131]=['',5.5927445,-0.783956,6];
ybs[8132]=['υ Cyg',5.5803715,0.6108797,4.43];
ybs[8133]=['',5.57545,0.9442405,6.13];
ybs[8134]=['',5.5902545,-0.4581208,6.56];
ybs[8135]=['',5.5852821,0.1973551,5.96];
ybs[8136]=['',5.5761988,0.9756673,5.98];
ybs[8137]=['θ1 Mic',5.595157,-0.7104301,4.82];
ybs[8138]=['',5.5978849,-0.869741,6.38];
ybs[8139]=['',5.5762661,1.0247744,6.42];
ybs[8140]=['68 Cyg',5.5822981,0.7688158,5];
ybs[8141]=['',5.5844888,0.7181173,6.15];
ybs[8142]=['',5.6132165,-1.2152318,6.41];
ybs[8143]=['',5.5865748,0.669193,5.83];
ybs[8144]=['',5.5909326,0.3862619,6.29];
ybs[8145]=['',5.6181477,-1.2512708,6.09];
ybs[8146]=['16 Aqr',5.5953282,-0.077752,5.87];
ybs[8147]=['',5.5864672,0.8659404,5.76];
ybs[8148]=['α Cep',5.5813625,1.09414,2.44];
ybs[8149]=['9 Equ',5.5950478,0.1301941,5.82];
ybs[8150]=['',5.5847988,1.024996,5.66];
ybs[8151]=['',5.5945433,0.4181972,5.57];
ybs[8152]=['',5.5932079,0.5682404,5.68];
ybs[8153]=['ι Cap',5.6007899,-0.2919735,4.28];
ybs[8154]=['',5.5649702,1.3459104,5.95];
ybs[8155]=['',5.5955244,0.5710365,6.04];
ybs[8156]=['',5.5937057,0.7059964,6.4];
ybs[8157]=['6 Cep',5.5845071,1.1340501,5.18];
ybs[8158]=['',5.604284,-0.3937986,5.6];
ybs[8159]=['1 Peg',5.5990925,0.347494,4.08];
ybs[8160]=['',5.5512057,1.4195202,6.15];
ybs[8161]=['17 Aqr',5.6035985,-0.1608076,5.99];
ybs[8162]=['',5.6620416,-1.441167,6.38];
ybs[8163]=['',5.6110708,-0.8117276,6.31];
ybs[8164]=['β Equ',5.6029764,0.1207231,5.16];
ybs[8165]=['',5.5902423,1.0622382,6.11];
ybs[8166]=['θ2 Mic',5.6110813,-0.7138435,5.77];
ybs[8167]=['γ Pav',5.6218454,-1.1389816,4.22];
ybs[8168]=['',5.601425,0.5308494,6.05];
ybs[8169]=['33 Cap',5.6092505,-0.3620795,5.41];
ybs[8170]=['',5.6091841,-0.3951536,6.38];
ybs[8171]=['',5.5975231,0.8638385,5.69];
ybs[8172]=['',5.6014655,0.6761387,6.63];
ybs[8173]=['18 Aqr',5.6091617,-0.2229088,5.49];
ybs[8174]=['γ Ind',5.6199055,-0.9521364,6.12];
ybs[8175]=['',5.6041466,0.6547188,6.58];
ybs[8176]=['',5.6072233,0.4255176,5.71];
ybs[8177]=['',5.6095061,0.1794293,6.35];
ybs[8178]=['20 Aqr',5.6118325,-0.0574523,6.36];
ybs[8179]=['',5.6059914,0.6537567,6.47];
ybs[8180]=['',5.6078163,0.443636,6.15];
ybs[8181]=['19 Aqr',5.6135553,-0.1682831,5.7];
ybs[8182]=['',5.6325744,-1.2112084,5.34];
ybs[8183]=['',5.6089839,0.4299661,6.32];
ybs[8184]=['',5.6097248,0.4586875,5.68];
ybs[8185]=['21 Aqr',5.6136763,-0.0602129,5.49];
ybs[8186]=['',5.6195192,-0.658378,5.63];
ybs[8187]=['',5.6567298,-1.3950279,6.47];
ybs[8188]=['',5.6225364,-0.7407241,5.51];
ybs[8189]=['',5.6160793,0.011194,6.46];
ybs[8190]=['ζ Cap',5.6202264,-0.3892807,3.74];
ybs[8191]=['',5.6187263,0.0211269,6.13];
ybs[8192]=['',5.6102883,0.8627135,6.58];
ybs[8193]=['35 Cap',5.6227192,-0.3680664,5.78];
ybs[8194]=['',5.6121963,0.8171825,5.6];
ybs[8195]=['69 Cyg',5.6146635,0.6418335,5.94];
ybs[8196]=['',5.6181344,0.3400368,6.07];
ybs[8197]=['',5.6317953,-0.9354559,6.39];
ybs[8198]=['',5.6267471,-0.2000239,6.61];
ybs[8199]=['36 Cap',5.6291755,-0.3787228,4.51];
ybs[8200]=['5 PsA',5.6309571,-0.5433293,6.5];
ybs[8201]=['70 Cyg',5.62152,0.649683,5.31];
ybs[8202]=['',5.618779,0.8542029,5.31];
ybs[8203]=['35 Vul',5.6232247,0.4837385,5.41];
ybs[8204]=['',5.6180343,0.9251254,6.03];
ybs[8205]=['',5.6270314,0.1449218,6.4];
ybs[8206]=['',5.6251205,0.5643173,5.8];
ybs[8207]=['',5.6293243,0.3144016,6.44];
ybs[8208]=['',5.6346364,-0.3322989,6.57];
ybs[8209]=['',5.6292013,0.3889901,5.93];
ybs[8210]=['',5.6203938,1.0447086,6.1];
ybs[8211]=['2 Peg',5.6333063,0.4144682,4.57];
ybs[8212]=['',5.6271564,0.9691212,6.12];
ybs[8213]=['7 Cep',5.6209985,1.1679158,5.44];
ybs[8214]=['71 Cyg',5.630233,0.8141739,5.24];
ybs[8215]=['ξ Gru',5.6445661,-0.7168054,5.29];
ybs[8216]=['6 PsA',5.6449111,-0.5905399,5.97];
ybs[8217]=['',5.6389173,0.213739,6.08];
ybs[8218]=['β Aqr',5.6411035,-0.0953317,2.91];
ybs[8219]=['',5.6503992,-0.9185285,6.41];
ybs[8220]=['',5.6804281,-1.3845808,6.18];
ybs[8221]=['',5.6459633,-0.4272772,6.43];
ybs[8222]=['',5.6503692,-0.7808411,5.57];
ybs[8223]=['',5.6337336,0.9261858,6.02];
ybs[8224]=['β Cep',5.6241727,1.2334,3.23];
ybs[8225]=['',5.602593,1.4072743,5.97];
ybs[8226]=['',5.6442455,0.4102176,6.7];
ybs[8227]=['',5.6541617,-0.7472624,6.32];
ybs[8228]=['',5.6386487,0.9202927,6.16];
ybs[8229]=['',5.6359348,1.0571134,5.53];
ybs[8230]=['',5.6562716,-0.5163714,6.41];
ybs[8231]=['37 Cap',5.6558396,-0.3486166,5.69];
ybs[8232]=['',5.6453129,0.8741866,5.75];
ybs[8233]=['',5.6577444,-0.4074267,6.4];
ybs[8234]=['',5.6470938,0.8022184,6.25];
ybs[8235]=['',5.6721988,-1.129451,6.2];
ybs[8236]=['',5.6535002,0.3990654,6.47];
ybs[8237]=['',5.6573585,-0.0675919,5.77];
ybs[8238]=['ρ Cyg',5.6500946,0.7976459,4.02];
ybs[8239]=['8 PsA',5.661828,-0.4548456,5.73];
ybs[8240]=['ν Oct',5.6902047,-1.3487438,3.76];
ybs[8241]=['72 Cyg',5.6538866,0.6744696,4.9];
ybs[8242]=['7 PsA',5.6648019,-0.5748621,6.11];
ybs[8243]=['',5.6566117,0.4940595,6.31];
ybs[8244]=['',5.657308,0.4286978,6.11];
ybs[8245]=['',5.6518319,0.9042252,6.15];
ybs[8246]=['ε Cap',5.665546,-0.3378111,4.68];
ybs[8247]=['',5.6605489,0.526499,6.36];
ybs[8248]=['',5.6591042,0.7938673,5.53];
ybs[8249]=['',5.6671752,-0.0048724,6.25];
ybs[8250]=['ξ Aqr',5.6681778,-0.1351406,4.69];
ybs[8251]=['3 Peg',5.6677257,0.1174518,6.18];
ybs[8252]=['74 Cyg',5.6633075,0.7072853,5.01];
ybs[8253]=['5 Peg',5.6675256,0.3391134,5.45];
ybs[8254]=['',5.6747884,-0.5858632,6.28];
ybs[8255]=['',5.6795351,-0.9118849,6.21];
ybs[8256]=['4 Peg',5.6712522,0.1026795,5.67];
ybs[8257]=['',5.6822251,-0.9708444,6.33];
ybs[8258]=['',5.665367,0.7820417,6.2];
ybs[8259]=['',5.6757326,-0.1826517,6.08];
ybs[8260]=['',5.6716892,0.4469855,6.16];
ybs[8261]=['',5.6656145,0.9451531,6.15];
ybs[8262]=['',5.6730088,0.3556434,5.85];
ybs[8263]=['25 Aqr',5.6758027,0.0411095,5.1];
ybs[8264]=['γ Cap',5.6786043,-0.288856,3.68];
ybs[8265]=['9 Cep',5.6661124,1.0854739,4.73];
ybs[8266]=['λ Oct',5.7351317,-1.4416994,5.29];
ybs[8267]=['',5.6711054,1.0053213,5.62];
ybs[8268]=['',5.6861431,-0.436147,6.49];
ybs[8269]=['42 Cap',5.6848904,-0.2432119,5.18];
ybs[8270]=['75 Cyg',5.6773236,0.7572258,5.11];
ybs[8271]=['41 Cap',5.6871594,-0.4040461,5.24];
ybs[8272]=['',5.7055021,-1.2373507,6.01];
ybs[8273]=['26 Aqr',5.6872288,0.0243989,5.67];
ybs[8274]=['κ Cap',5.6898536,-0.3273109,4.73];
ybs[8275]=['7 Peg',5.6875182,0.1011017,5.3];
ybs[8276]=['',5.6790168,0.9596577,6.2];
ybs[8277]=['76 Cyg',5.6834909,0.7141487,6.11];
ybs[8278]=['',5.6886682,0.1908956,6.09];
ybs[8279]=['',5.6923442,-0.340475,6.22];
ybs[8280]=['',5.9964448,-1.5475729,6.57];
ybs[8281]=['44 Cap',5.69155,-0.2493505,5.88];
ybs[8282]=['',5.6856449,0.6217363,6.07];
ybs[8283]=['',5.6857544,0.80073,6.17];
ybs[8284]=['',5.698448,-0.6708873,6.3];
ybs[8285]=['77 Cyg',5.68702,0.7188998,5.69];
ybs[8286]=['π1 Cyg',5.6852661,0.8953942,4.67];
ybs[8287]=['45 Cap',5.6956726,-0.2554488,5.99];
ybs[8288]=['',5.7025315,-0.8619279,6.45];
ybs[8289]=['',5.6877782,0.8676565,6.09];
ybs[8290]=['ι PsA',5.7002371,-0.5744262,4.34];
ybs[8291]=['',5.6901857,0.7202616,5.49];
ybs[8292]=['79 Cyg',5.6917047,0.6701531,5.65];
ybs[8293]=['ε Peg',5.6958349,0.1743295,2.39];
ybs[8294]=['μ1 Cyg',5.6951493,0.5036337,4.73];
ybs[8295]=['μ2 Cyg',5.6951274,0.5036385,6.08];
ybs[8296]=['46 Cap',5.6998467,-0.1565363,5.09];
ybs[8297]=['',5.6875702,1.0364445,6.08];
ybs[8298]=['9 Peg',5.6970705,0.3047945,4.34];
ybs[8299]=['',5.6971787,0.2597991,5.94];
ybs[8300]=['κ Peg',5.697431,0.4495702,4.13];
ybs[8301]=['μ Cep',5.6908969,1.0278774,4.08];
ybs[8302]=['11 Cep',5.6821871,1.2465809,4.56];
ybs[8303]=['47 Cap',5.7053764,-0.1599034,6];
ybs[8304]=['λ Cap',5.7065741,-0.1963793,5.58];
ybs[8305]=['',5.7019018,0.6278132,6.4];
ybs[8306]=['12 Peg',5.7037435,0.4025225,5.29];
ybs[8307]=['δ Cap',5.7088912,-0.2794782,2.87];
ybs[8308]=['',5.7152472,-0.8236014,5.58];
ybs[8309]=['',5.686968,1.264196,5.17];
ybs[8310]=['',5.7050982,0.4481549,6.28];
ybs[8311]=['θ PsA',5.7123284,-0.5372786,5.01];
ybs[8312]=['',5.6965896,1.092123,5.95];
ybs[8313]=['11 Peg',5.7092982,0.0488772,5.64];
ybs[8314]=['',5.7039621,0.7535428,6.54];
ybs[8315]=['',5.7082806,0.3020894,6.21];
ybs[8316]=['',5.7241227,-1.1274332,5.62];
ybs[8317]=['',5.7112459,-0.1012769,6.17];
ybs[8318]=['ο Ind',5.7282492,-1.2132453,5.53];
ybs[8319]=['ν Cep',5.699181,1.0687438,4.29];
ybs[8320]=['π2 Cyg',5.7059221,0.8626045,4.23];
ybs[8321]=['',5.7123492,0.6404513,6.47];
ybs[8322]=['',5.7203447,-0.2200495,6.31];
ybs[8323]=['',5.7138118,0.6765477,6.12];
ybs[8324]=['12 Cep',5.7078637,1.0612841,5.52];
ybs[8325]=['',5.7227685,-0.2919832,6.38];
ybs[8326]=['',5.7185466,0.3591459,6.29];
ybs[8327]=['',5.7047936,1.2263552,6.29];
ybs[8328]=['14 Peg',5.7200126,0.5286485,5.04];
ybs[8329]=['13 Peg',5.7216676,0.3037017,5.29];
ybs[8330]=['',5.7188664,0.7201923,6.48];
ybs[8331]=['',5.729258,-0.3230127,6.16];
ybs[8332]=['',5.7161063,1.0714177,6.17];
ybs[8333]=['',5.7278269,0.3480602,5.77];
ybs[8334]=['',5.7251056,0.6920616,6.17];
ybs[8335]=['',5.7309929,0.3733085,6.89];
ybs[8336]=['μ Cap',5.7361153,-0.2344916,5.08];
ybs[8337]=['',5.7463732,-1.0780802,5.9];
ybs[8338]=['γ Gru',5.7395162,-0.6501089,3.01];
ybs[8339]=['15 Peg',5.7316434,0.5045631,5.53];
ybs[8340]=['',5.7373649,-0.1779413,6.59];
ybs[8341]=['16 Peg',5.7341908,0.4545043,5.08];
ybs[8342]=['',5.7283773,0.9758616,5.71];
ybs[8343]=['',5.7367968,0.3453081,5.68];
ybs[8344]=['',5.7385752,0.1218401,6.15];
ybs[8345]=['',5.7397331,-0.072598,5.71];
ybs[8346]=['',5.7256984,1.1496209,6.37];
ybs[8347]=['π Ind',5.7505651,-1.0084896,6.19];
ybs[8348]=['',5.7415662,-0.0555789,6.2];
ybs[8349]=['',5.7397066,0.3461843,6.39];
ybs[8350]=['',5.7480126,-0.5321383,6.41];
ybs[8351]=['',5.7501933,-0.6481517,5.46];
ybs[8352]=['',5.7530706,-0.6567585,6.18];
ybs[8353]=['δ Ind',5.7576944,-0.9577451,4.4];
ybs[8354]=['κ1 Ind',5.7605148,-1.0278993,6.12];
ybs[8355]=['',5.7785054,-1.3533887,6.41];
ybs[8356]=['13 Cep',5.7408806,0.9900919,5.8];
ybs[8357]=['',5.748883,0.3727488,6.4];
ybs[8358]=['17 Peg',5.7514531,0.2128214,5.54];
ybs[8359]=['',5.742482,1.0761486,6.13];
ybs[8360]=['',5.7428321,1.1421033,5.86];
ybs[8361]=['',5.7574189,-0.0926236,6.33];
ybs[8362]=['',5.750754,0.8514761,6.42];
ybs[8363]=['',5.759986,-0.3676506,6.12];
ybs[8364]=['',5.7629449,-0.6680622,5.5];
ybs[8365]=['',5.783073,-1.3264372,5.95];
ybs[8366]=['',5.7685548,-0.9732705,6.01];
ybs[8367]=['',5.7604287,-0.0742649,6.22];
ybs[8368]=['',5.7479946,1.1125213,4.91];
ybs[8369]=['',5.7500313,1.1566905,6.43];
ybs[8370]=['18 Peg',5.7654994,0.119308,6];
ybs[8371]=['η PsA',5.7693522,-0.4945396,5.42];
ybs[8372]=['ε Ind',5.7814825,-0.9890215,4.69];
ybs[8373]=['',5.7578973,1.0963501,5.93];
ybs[8374]=['',5.7604878,1.0083884,6.59];
ybs[8375]=['28 Aqr',5.7697789,0.0126297,5.58];
ybs[8376]=['',5.7662454,0.578132,6.46];
ybs[8377]=['20 Peg',5.7695465,0.2310526,5.6];
ybs[8378]=['19 Peg',5.7699263,0.1461864,5.65];
ybs[8379]=['',5.7750325,-0.3104006,6.28];
ybs[8380]=['',5.7511563,1.3109898,6.35];
ybs[8381]=['29 Aqr',5.7760798,-0.294003,6.37];
ybs[8382]=['',5.7736666,0.1936056,6.37];
ybs[8383]=['',5.7800685,-0.5198441,7.1];
ybs[8384]=['',5.7656501,1.0926892,6.66];
ybs[8385]=['16 Cep',5.7578237,1.2792905,5.03];
ybs[8386]=['30 Aqr',5.779479,-0.1117574,5.54];
ybs[8387]=['ο Aqr',5.7795722,-0.0355349,4.69];
ybs[8388]=['',5.7715598,0.9250424,5.78];
ybs[8389]=['21 Peg',5.7793034,0.2008117,5.8];
ybs[8390]=['13 PsA',5.7849137,-0.5200567,6.47];
ybs[8391]=['14 Cep',5.7722535,1.0143752,5.56];
ybs[8392]=['',5.776768,0.7813688,5.6];
ybs[8393]=['',5.7857711,-0.4660522,5.96];
ybs[8394]=['κ2 Ind',5.7925072,-1.038751,5.62];
ybs[8395]=['32 Aqr',5.7859827,-0.0137351,5.3];
ybs[8396]=['λ Gru',5.7926768,-0.6880651,4.46];
ybs[8397]=['',5.7842776,0.5770331,6.38];
ybs[8398]=['ν Peg',5.7897477,0.090383,4.84];
ybs[8399]=['α Aqr',5.7903053,-0.0034861,2.96];
ybs[8400]=['',5.7871379,0.4676381,5.78];
ybs[8401]=['18 Cep',5.7797242,1.1037304,5.29];
ybs[8402]=['ξ Cep',5.779165,1.1300503,4.29];
ybs[8403]=['ι Aqr',5.7934276,-0.2399748,4.27];
ybs[8404]=['23 Peg',5.788776,0.507608,5.7];
ybs[8405]=['',5.8160732,-1.3222452,6.55];
ybs[8406]=['',5.7868627,0.8179402,6.13];
ybs[8407]=['',5.7894223,0.7894506,6.44];
ybs[8408]=['',5.7475027,1.4483992,6.98];
ybs[8409]=['',5.7902569,0.787745,5.14];
ybs[8410]=['α Gru',5.8021748,-0.8175192,1.74];
ybs[8411]=['20 Cep',5.784693,1.0979037,5.27];
ybs[8412]=['',5.7893344,0.8438952,6.27];
ybs[8413]=['19 Cep',5.7853539,1.0890808,5.11];
ybs[8414]=['',5.7910055,0.7918328,6.19];
ybs[8415]=['ι Peg',5.7951271,0.4444536,3.76];
ybs[8416]=['μ PsA',5.802362,-0.5736527,4.5];
ybs[8417]=['',5.8213972,-1.3263508,6.15];
ybs[8418]=['υ PsA',5.8026092,-0.5920705,4.99];
ybs[8419]=['',5.7905264,0.9854671,6.39];
ybs[8420]=['',5.7972923,0.3420148,5.75];
ybs[8421]=['',5.7974263,0.3162713,6.35];
ybs[8422]=['',5.8037964,-0.5760412,6.37];
ybs[8423]=['25 Peg',5.7988221,0.3808888,5.78];
ybs[8424]=['35 Aqr',5.8046307,-0.3211202,5.81];
ybs[8425]=['',5.8097686,-0.8375143,6.43];
ybs[8426]=['',5.8006911,0.4479262,6.11];
ybs[8427]=['',5.7944477,1.0290662,6.32];
ybs[8428]=['',5.7959587,0.9324878,6.14];
ybs[8429]=['',5.8091223,-0.5915542,5.37];
ybs[8430]=['',5.799855,0.8712165,6.42];
ybs[8431]=['',5.8092928,-0.4916824,6.44];
ybs[8432]=['τ PsA',5.8100337,-0.5659599,4.92];
ybs[8433]=['',5.8018155,0.8004552,6.11];
ybs[8434]=['π1 Peg',5.8046029,0.581075,5.58];
ybs[8435]=['θ Peg',5.8094517,0.1102871,3.53];
ybs[8436]=['',5.8103063,-0.0658497,6.27];
ybs[8437]=['38 Aqr',5.8116472,-0.1997296,5.46];
ybs[8438]=['',5.8112369,-0.0723598,6.01];
ybs[8439]=['π2 Peg',5.8079235,0.5811853,4.29];
ybs[8440]=['',5.8096904,0.3444962,6.18];
ybs[8441]=['',5.8100274,0.2574579,6.33];
ybs[8442]=['',5.8136562,-0.3684572,6.09];
ybs[8443]=['',5.8111936,0.2050022,5.78];
ybs[8444]=['28 Peg',5.8104761,0.3682529,6.46];
ybs[8445]=['',5.8118085,0.5353697,6.32];
ybs[8446]=['',5.816482,0.2820837,5.95];
ybs[8447]=['39 Aqr',5.8195712,-0.2456039,6.03];
ybs[8448]=['',5.8124488,0.8891539,5.4];
ybs[8449]=['',5.8221328,-0.4573776,6.17];
ybs[8450]=['ζ Cep',5.8106714,1.0179187,3.35];
ybs[8451]=['',5.8175061,0.437584,5.92];
ybs[8452]=['',5.8206977,-0.0802666,6.39];
ybs[8453]=['24 Cep',5.8043908,1.2647024,4.79];
ybs[8454]=['λ Cep',5.8134607,1.0390985,5.04];
ybs[8455]=['',5.8255138,-0.4373561,5.58];
ybs[8456]=['ψ Oct',5.8474766,-1.3506807,5.51];
ybs[8457]=['',5.8149566,0.9941578,5.24];
ybs[8458]=['',5.8064077,1.2606903,6.37];
ybs[8459]=['',5.8084994,1.2261639,5.5];
ybs[8460]=['',5.8201627,0.6060937,5.33];
ybs[8461]=['',5.8153967,1.0333458,6.3];
ybs[8462]=['',5.8298817,-0.7201095,6.23];
ybs[8463]=['λ PsA',5.8280721,-0.4824895,5.43];
ybs[8464]=['',5.8156355,1.0625756,5.35];
ybs[8465]=['41 Aqr',5.8278657,-0.3656786,5.32];
ybs[8466]=['ε Oct',5.8583018,-1.4017747,5.1];
ybs[8467]=['',5.8240161,0.501441,5.89];
ybs[8468]=['',5.8169073,1.1067678,5.79];
ybs[8469]=['',5.8340865,-0.7736916,6.1];
ybs[8470]=['',5.8247401,0.6952898,4.49];
ybs[8471]=['μ1 Gru',5.8341123,-0.7194943,4.79];
ybs[8472]=['',5.824296,0.7952242,5.53];
ybs[8473]=['μ2 Gru',5.8377288,-0.7243919,5.1];
ybs[8474]=['',5.8284008,0.7518231,5.71];
ybs[8475]=['',5.8232875,1.104525,6.11];
ybs[8476]=['',5.8347071,0.1513584,6.21];
ybs[8477]=['',5.8380979,-0.4498656,6.15];
ybs[8478]=['',5.8176885,1.2815787,6.08];
ybs[8479]=['ε Cep',5.8290398,0.9977361,4.19];
ybs[8480]=['',5.8373499,-0.0257172,6.15];
ybs[8481]=['42 Aqr',5.8386101,-0.2218038,5.34];
ybs[8482]=['',5.8396655,-0.401722,6.17];
ybs[8483]=['1 Lac',5.8339417,0.6609843,4.13];
ybs[8484]=['θ Aqr',5.8386545,-0.1336985,4.16];
ybs[8485]=['',5.8388666,-0.1556313,5.79];
ybs[8486]=['',5.8461149,-0.933828,5.37];
ybs[8487]=['α Tuc',5.8475619,-1.049576,2.86];
ybs[8488]=['',5.836492,0.4874188,6.37];
ybs[8489]=['44 Aqr',5.8398105,-0.0918773,5.75];
ybs[8490]=['υ Oct',5.9152785,-1.4944879,5.77];
ybs[8491]=['',5.8351943,1.0008259,5.88];
ybs[8492]=['',5.8439213,-0.0019982,6.39];
ybs[8493]=['45 Aqr',5.8482643,-0.23006,5.95];
ybs[8494]=['',5.8565198,-1.0015752,6.34];
ybs[8495]=['',5.8469072,0.6613564,6.17];
ybs[8496]=['25 Cep',5.8425403,1.0982957,5.75];
ybs[8497]=['ρ Aqr',5.8533334,-0.1343429,5.37];
ybs[8498]=['30 Peg',5.8542371,0.1032071,5.37];
ybs[8499]=['',5.8562462,0.1450486,6.17];
ybs[8500]=['ν Ind',5.8754293,-1.2589159,5.29];
ybs[8501]=['47 Aqr',5.8596728,-0.3747946,5.13];
ybs[8502]=['',5.8561976,0.4722737,6.47];
ybs[8503]=['γ Aqr',5.8595834,-0.0220441,3.84];
ybs[8504]=['',5.8539809,0.8919459,6.42];
ybs[8505]=['31 Peg',5.8587426,0.2151891,5.01];
ybs[8506]=['π1 Gru',5.8652549,-0.7997675,6.62];
ybs[8507]=['32 Peg',5.8575647,0.4966273,4.81];
ybs[8508]=['2 Lac',5.8557632,0.8143823,4.57];
ybs[8509]=['π2 Gru',5.8670029,-0.7994313,5.62];
ybs[8510]=['',5.8408225,1.3371192,6.66];
ybs[8511]=['',5.8813573,-1.3070816,6.04];
ybs[8512]=['',5.8775706,-1.2270809,5.78];
ybs[8513]=['',5.8595069,0.7365733,6.41];
ybs[8514]=['49 Aqr',5.8681136,-0.4300116,5.53];
ybs[8515]=['',5.8678773,-0.1233912,5.93];
ybs[8516]=['',5.8753943,-1.0065697,5.32];
ybs[8517]=['33 Peg',5.8679328,0.3660479,6.04];
ybs[8518]=['51 Aqr',5.8703674,-0.0822427,5.78];
ybs[8519]=['50 Aqr',5.8719855,-0.233954,5.76];
ybs[8520]=['',5.8639176,1.0019747,6.16];
ybs[8521]=['',5.8685929,0.6754131,6.22];
ybs[8522]=['',5.8635673,1.0916068,6.04];
ybs[8523]=['β Lac',5.8666136,0.9137459,4.43];
ybs[8524]=['π Aqr',5.8753319,0.0262246,4.66];
ybs[8525]=['δ Tuc',5.8863119,-1.1316853,4.48];
ybs[8526]=['4 Lac',5.8709158,0.8657049,4.57];
ybs[8527]=['',5.8797017,-0.4111509,6.29];
ybs[8528]=['',5.8767915,0.3241006,6.26];
ybs[8529]=['53 Aqr',5.8812799,-0.2900042,6.57];
ybs[8530]=['53 Aqr',5.8812945,-0.2900236,6.35];
ybs[8531]=['',5.8064143,1.5015193,5.27];
ybs[8532]=['',5.892183,-1.1757108,5.55];
ybs[8533]=['34 Peg',5.8811587,0.0788713,5.75];
ybs[8534]=['',5.8811163,0.6557076,6.46];
ybs[8535]=['',5.8637374,1.3677776,6.76];
ybs[8536]=['35 Peg',5.8865434,0.0841462,4.79];
ybs[8537]=['ν Gru',5.8908338,-0.6807845,5.47];
ybs[8538]=['',5.8840253,0.697002,6.14];
ybs[8539]=['',5.8814041,0.9871366,6.57];
ybs[8540]=['',5.885658,0.5579104,5.98];
ybs[8541]=['δ1 Gru',5.8936425,-0.7569413,3.97];
ybs[8542]=['',5.8758481,1.2373686,5.47];
ybs[8543]=['ζ1 Aqr',5.8908522,0.0018433,4.59];
ybs[8544]=['ζ2 Aqr',5.8908813,0.0018481,4.42];
ybs[8545]=['δ2 Gru',5.8957765,-0.7613706,4.11];
ybs[8546]=['26 Cep',5.8812375,1.1389609,5.46];
ybs[8547]=['36 Peg',5.8920359,0.1615276,5.58];
ybs[8548]=['',5.89541,-0.4709092,5.95];
ybs[8549]=['',5.8918818,0.4693018,5.79];
ybs[8550]=['',5.8962807,-0.2232073,6.4];
ybs[8551]=['37 Peg',5.8957487,0.0795489,5.48];
ybs[8552]=['56 Aqr',5.8974643,-0.2523677,6.37];
ybs[8553]=['',5.8867661,1.1206982,6.29];
ybs[8554]=['',5.8941396,0.6257291,6.56];
ybs[8555]=['ζ PsA',5.9003076,-0.4528648,6.43];
ybs[8556]=['δ Cep',5.8908598,1.0217368,3.75];
ybs[8557]=['5 Lac',5.8929113,0.8348428,4.36];
ybs[8558]=['σ Aqr',5.8989552,-0.1841627,4.82];
ybs[8559]=['38 Peg',5.8955164,0.5706991,5.65];
ybs[8560]=['',5.8953812,0.8636284,6.4];
ybs[8561]=['β PsA',5.9030942,-0.5623381,4.29];
ybs[8562]=['',5.9239473,-1.3726004,6.15];
ybs[8563]=['ρ1 Cep',5.8768253,1.3772586,5.83];
ybs[8564]=['6 Lac',5.8972501,0.7548476,4.51];
ybs[8565]=['',5.901711,-0.0486015,6.16];
ybs[8566]=['',5.9017682,-0.1121993,6.14];
ybs[8567]=['ν Tuc',5.9107012,-1.0795797,4.81];
ybs[8568]=['58 Aqr',5.9035022,-0.1881294,6.38];
ybs[8569]=['',5.9023097,0.5178265,6.35];
ybs[8570]=['α Lac',5.900523,0.8798016,3.77];
ybs[8571]=['39 Peg',5.906948,0.3552919,6.42];
ybs[8572]=['',5.9078518,0.27908,6.32];
ybs[8573]=['',5.9058755,0.6964981,5.88];
ybs[8574]=['',5.9048377,0.9453426,6.35];
ybs[8575]=['60 Aqr',5.9136523,-0.0252569,5.89];
ybs[8576]=['ρ2 Cep',5.8908099,1.3779405,5.5];
ybs[8577]=['υ Aqr',5.9167676,-0.3592085,5.2];
ybs[8578]=['',5.9229744,-1.0080347,6.23];
ybs[8579]=['',5.9106823,0.9905082,5.71];
ybs[8580]=['',5.9068561,1.2224353,6.6];
ybs[8581]=['',5.9208136,-0.4165004,5.97];
ybs[8582]=['η Aqr',5.9193389,0.0001716,4.02];
ybs[8583]=['',5.9078262,1.2304696,6.34];
ybs[8584]=['',5.9022337,1.3326106,5.68];
ybs[8585]=['σ1 Gru',5.9249951,-0.7060762,6.28];
ybs[8586]=['',5.9252368,-0.5504119,5.82];
ybs[8587]=['σ2 Gru',5.9271365,-0.7062198,5.86];
ybs[8588]=['8 Lac',5.9208804,0.6939707,5.73];
ybs[8589]=['',5.9221128,0.6231647,6.1];
ybs[8590]=['',5.9246171,0.2063772,6.4];
ybs[8591]=['',5.9206483,0.8761297,6.29];
ybs[8592]=['',5.9202826,0.9808298,6.38];
ybs[8593]=['',5.9266627,0.2217426,6.3];
ybs[8594]=['',5.9250829,0.624481,6.3];
ybs[8595]=['κ Aqr',5.9298706,-0.0715623,5.03];
ybs[8596]=['',5.9369135,-0.9174162,6.65];
ybs[8597]=['',5.9325996,-0.1356087,6.23];
ybs[8598]=['9 Lac',5.9270839,0.9018642,4.63];
ybs[8599]=['',5.9345711,-0.4995084,6.47];
ybs[8600]=['31 Cep',5.9182177,1.2875365,5.08];
ybs[8601]=['',5.935149,-0.5751438,5.66];
ybs[8602]=['',5.9313116,0.7908259,6.4];
ybs[8603]=['40 Peg',5.9344055,0.3429622,5.82];
ybs[8604]=['',5.9388705,-0.4921309,6.31];
ybs[8605]=['',5.9444015,-0.9999642,5.97];
ybs[8606]=['',5.9323677,0.9935082,5.21];
ybs[8607]=['10 Lac',5.9357134,0.6837923,4.88];
ybs[8608]=['',5.9416956,-0.5328579,5.87];
ybs[8609]=['41 Peg',5.9383622,0.3457386,6.21];
ybs[8610]=['',5.9242255,1.3177116,5.79];
ybs[8611]=['',5.9371022,0.6583552,6.03];
ybs[8612]=['30 Cep',5.9320464,1.1119917,5.19];
ybs[8613]=['ε PsA',5.9428716,-0.4697584,4.17];
ybs[8614]=['',5.9431379,-0.0597898,6.31];
ybs[8615]=['β Oct',5.9708845,-1.4181164,4.15];
ybs[8616]=['',5.9432224,0.2561779,5.71];
ybs[8617]=['11 Lac',5.94106,0.7750095,4.46];
ybs[8618]=['',5.9398304,0.9420318,5.93];
ybs[8619]=['ζ Peg',5.9458234,0.1912877,3.4];
ybs[8620]=['',5.9518509,-0.821731,5.98];
ybs[8621]=['β Gru',5.9520736,-0.816044,2.1];
ybs[8622]=['19 PsA',5.9503725,-0.5101956,6.17];
ybs[8623]=['',5.9457679,0.5427001,6.34];
ybs[8624]=['',5.9522181,-0.7700204,6.07];
ybs[8625]=['12 Lac',5.9453682,0.7043125,5.25];
ybs[8626]=['ο Peg',5.946818,0.5137576,4.79];
ybs[8627]=['',5.9479308,0.2556048,5.9];
ybs[8628]=['',5.9458825,0.7274191,5.94];
ybs[8629]=['ρ Gru',5.9555504,-0.7205668,4.85];
ybs[8630]=['',5.9530558,-0.142816,6.45];
ybs[8631]=['',5.9595997,-1.0536601,6.3];
ybs[8632]=['67 Aqr',5.9538223,-0.1192729,6.41];
ybs[8633]=['',5.9487565,0.9431346,6.12];
ybs[8634]=['66 Aqr',5.955518,-0.3263986,4.69];
ybs[8635]=['η Peg',5.9522407,0.5297121,2.94];
ybs[8636]=['',5.9517546,0.662032,6.43];
ybs[8637]=['',5.9521736,0.8254971,6.39];
ybs[8638]=['',5.9556436,0.1931765,6.51];
ybs[8639]=['',5.9567956,0.6910569,5.95];
ybs[8640]=['η Gru',5.9651847,-0.9314972,4.85];
ybs[8641]=['13 Lac',5.9567624,0.7321351,5.08];
ybs[8642]=['',5.9651663,-0.8101482,5.51];
ybs[8643]=['',5.9672168,-0.8525825,6.62];
ybs[8644]=['',5.9687107,-0.8649199,6.48];
ybs[8645]=['45 Peg',5.9631996,0.3402699,6.25];
ybs[8646]=['',5.9596338,0.9188538,6.55];
ybs[8647]=['',5.9697476,-0.8169856,6.56];
ybs[8648]=['ξ Oct',5.9888984,-1.3961556,5.35];
ybs[8649]=['',5.9849115,-1.3425134,6.73];
ybs[8650]=['ξ Peg',5.9686426,0.2147169,4.19];
ybs[8651]=['',5.965782,0.7797362,5.76];
ybs[8652]=['λ Peg',5.9677726,0.4135577,3.95];
ybs[8653]=['',5.9720287,-0.5939598,6.28];
ybs[8654]=['',5.9774021,-1.0743241,6.37];
ybs[8655]=['68 Aqr',5.9728093,-0.3400525,5.26];
ybs[8656]=['',5.9741412,-0.6648332,6.71];
ybs[8657]=['',5.9821271,-1.2255294,6.34];
ybs[8658]=['τ1 Aqr',5.9734366,-0.2430651,5.66];
ybs[8659]=['',5.9745822,-0.4499827,6.3];
ybs[8660]=['ε Gru',5.9778275,-0.8933814,3.49];
ybs[8661]=['70 Aqr',5.9768353,-0.1819615,6.19];
ybs[8662]=['',5.9706242,1.0229807,6.36];
ybs[8663]=['',5.9747393,0.6553104,5.9];
ybs[8664]=['τ2 Aqr',5.9816217,-0.2349628,4.01];
ybs[8665]=['',5.9836223,-0.5702876,6.33];
ybs[8666]=['',5.9810844,0.185162,6.54];
ybs[8667]=['',5.9769688,0.9519891,6.12];
ybs[8668]=['',5.9762988,1.100749,6.06];
ybs[8669]=['μ Peg',5.9829202,0.4316523,3.48];
ybs[8670]=['',5.9883227,-0.6811419,5.42];
ybs[8671]=['',5.9920481,-1.0428493,6.46];
ybs[8672]=['',5.9770749,1.1990457,6.19];
ybs[8673]=['',5.9811722,0.9779588,5.43];
ybs[8674]=['',5.994033,-1.1005699,6.12];
ybs[8675]=['14 Lac',5.9841802,0.7345018,5.92];
ybs[8676]=['',5.9858301,0.3363449,6.4];
ybs[8677]=['',5.9831061,0.8867521,6.21];
ybs[8678]=['21 PsA',5.9895081,-0.5132258,5.97];
ybs[8679]=['ι Cep',5.9802101,1.1576885,3.52];
ybs[8680]=['γ PsA',5.9946974,-0.5715065,4.46];
ybs[8681]=['',5.9879545,1.0790909,5.6];
ybs[8682]=['σ Peg',5.9935912,0.1739424,5.16];
ybs[8683]=['λ Aqr',5.994732,-0.1300108,3.74];
ybs[8684]=['15 Lac',5.9914517,0.758224,4.94];
ybs[8685]=['τ1 Gru',5.999839,-0.8459125,6.04];
ybs[8686]=['ρ Ind',6.0053678,-1.2207283,6.05];
ybs[8687]=['',5.9660142,1.4535712,4.74];
ybs[8688]=['',5.9962823,0.2962143,5.64];
ybs[8689]=['74 Aqr',5.9985468,-0.2004662,5.8];
ybs[8690]=['',5.994935,0.8821352,6.46];
ybs[8691]=['',5.9965693,0.7033321,6.34];
ybs[8692]=['',5.9954113,1.0512435,6.01];
ybs[8693]=['',5.9985682,0.7833035,5.81];
ybs[8694]=['δ Aqr',6.0037082,-0.2738393,3.27];
ybs[8695]=['78 Aqr',6.003248,-0.12346,6.19];
ybs[8696]=['77 Aqr',6.0041863,-0.2817123,5.56];
ybs[8697]=['',6.0006235,0.7069951,5.81];
ybs[8698]=['',6.0066149,-0.6328129,6.4];
ybs[8699]=['',6.0030781,0.297974,6.12];
ybs[8700]=['1 Psc',6.0049989,0.0208702,6.11];
ybs[8701]=['',6.0059049,-0.0847653,5.72];
ybs[8702]=['ρ Peg',6.0059419,0.1561532,4.9];
ybs[8703]=['',6.0047375,0.6494021,5.91];
ybs[8704]=['',6.0091861,-0.5498111,6.1];
ybs[8705]=['δ PsA',6.0095999,-0.5656352,4.21];
ybs[8706]=['',6.0115526,-0.5486315,6.48];
ybs[8707]=['τ3 Gru',6.0135974,-0.8349273,5.7];
ybs[8708]=['',6.0078027,0.6367456,5.74];
ybs[8709]=['',6.0130215,0.2090849,6.51];
ybs[8710]=['16 Lac',6.0105524,0.7284159,5.59];
ybs[8711]=['',6.0105375,0.8703064,4.95];
ybs[8712]=['',6.0150787,-0.0816565,6.31];
ybs[8713]=['α PsA',6.0169774,-0.5147104,1.16];
ybs[8714]=['51 Peg',6.0155621,0.3647797,5.49];
ybs[8715]=['',6.0161172,0.0687964,6.28];
ybs[8716]=['',6.0133725,0.8519919,5.43];
ybs[8717]=['',6.0211276,-0.6176968,6.13];
ybs[8718]=['',6.0162147,0.6883642,6.18];
ybs[8719]=['',6.019291,-0.039509,6.16];
ybs[8720]=['',6.0198765,-0.0223172,6.37];
ybs[8721]=['',5.9788621,1.4874608,5.9];
ybs[8722]=['',6.0205872,0.1656068,6.43];
ybs[8723]=['',6.0211553,0.1303999,6.33];
ybs[8724]=['52 Peg',6.0232307,0.2070067,5.75];
ybs[8725]=['',6.0254507,-0.5119126,5.51];
ybs[8726]=['',6.0252396,-0.225829,6.07];
ybs[8727]=['2 Psc',6.0244832,0.0191033,5.43];
ybs[8728]=['',6.0275718,-0.4368961,5.65];
ybs[8729]=['',6.0224374,0.9212921,6.29];
ybs[8730]=['',6.0220881,1.0462623,6.43];
ybs[8731]=['',6.0289435,-0.4449674,6.29];
ybs[8732]=['ζ Gru',6.031492,-0.9184302,4.12];
ybs[8733]=['',5.9956043,1.4744007,4.71];
ybs[8734]=['',6.0325088,-0.886941,5.68];
ybs[8735]=['3 Psc',6.0296185,0.0055462,6.21];
ybs[8736]=['',6.0299523,0.0548665,5.83];
ybs[8737]=['',6.026296,0.9961837,5];
ybs[8738]=['',6.0295911,0.5448046,6.6];
ybs[8739]=['',6.0329667,-0.5012858,5.55];
ybs[8740]=['',6.0287525,0.7942456,6.5];
ybs[8741]=['',6.0331503,-0.3954702,6.28];
ybs[8742]=['81 Aqr',6.0330176,-0.1209348,6.21];
ybs[8743]=['',6.0303599,0.6778864,6.54];
ybs[8744]=['',6.0335818,-0.079924,5.94];
ybs[8745]=['',6.0384936,-0.6333555,6.47];
ybs[8746]=['',6.0325375,0.9989848,6.2];
ybs[8747]=['ο And',6.0347102,0.741036,3.62];
ybs[8748]=['82 Aqr',6.038029,-0.1124331,6.15];
ybs[8749]=['',6.0390344,-0.3619516,5.97];
ybs[8750]=['',6.0376356,0.556983,6.57];
ybs[8751]=['2 And',6.0376899,0.7485718,5.1];
ybs[8752]=['π PsA',6.0425194,-0.6041819,5.11];
ybs[8753]=['',6.0383134,0.7712808,6.39];
ybs[8754]=['',6.0495525,-1.1988264,5.52];
ybs[8755]=['',6.0379343,0.9663649,6.5];
ybs[8756]=['',6.0447839,-0.7216219,5.79];
ybs[8757]=['',6.0441691,-0.0813821,6.68];
ybs[8758]=['β Psc',6.0437414,0.0689827,4.53];
ybs[8759]=['κ Gru',6.0479687,-0.9395536,5.37];
ybs[8760]=['β Peg',6.0430304,0.4924477,2.42];
ybs[8761]=['',6.0443242,0.1177941,6.41];
ybs[8762]=['',6.0406584,1.0572787,6.74];
ybs[8763]=['',6.0405791,1.0244568,6.43];
ybs[8764]=['',6.0409728,1.1753311,5.24];
ybs[8765]=['3 And',6.0444378,0.8758878,4.65];
ybs[8766]=['α Peg',6.0474882,0.2676955,2.49];
ybs[8767]=['83 Aqr',6.0494689,-0.1319645,5.43];
ybs[8768]=['',6.0497829,-0.2957732,6.14];
ybs[8769]=['',6.0489729,0.291394,6.44];
ybs[8770]=['',6.0499481,0.0251251,6.39];
ybs[8771]=['',6.0662877,-1.3848795,6.12];
ybs[8772]=['θ Gru',6.0573962,-0.7572585,4.28];
ybs[8773]=['',6.0541847,0.3255083,6.13];
ybs[8774]=['86 Aqr',6.056257,-0.4120765,4.47];
ybs[8775]=['υ Gru',6.0573773,-0.6764787,5.61];
ybs[8776]=['',6.0587286,-0.8634805,6.33];
ybs[8777]=['',6.0551669,0.3498271,6.3];
ybs[8778]=['',6.0591303,-0.882325,5.83];
ybs[8779]=['',6.0661291,-1.2820019,6.15];
ybs[8780]=['55 Peg',6.05734,0.1665445,4.52];
ybs[8781]=['56 Peg',6.0576455,0.4468252,4.76];
ybs[8782]=['1 Cas',6.0548049,1.0393873,4.85];
ybs[8783]=['',6.0590766,0.5752386,6.02];
ybs[8784]=['',6.0592879,0.3711806,5.99];
ybs[8785]=['',6.058166,0.8063585,6.66];
ybs[8786]=['',6.0574313,0.9241388,6.11];
ybs[8787]=['',6.063587,-0.5007402,5.6];
ybs[8788]=['',6.0572394,1.0447604,6.4];
ybs[8789]=['4 And',6.0597052,0.8119299,5.33];
ybs[8790]=['5 And',6.0600896,0.8626949,5.7];
ybs[8791]=['',6.0621451,0.7800692,6.56];
ybs[8792]=['5 Psc',6.0647258,0.0394593,5.4];
ybs[8793]=['',6.0597996,1.1129315,6.26];
ybs[8794]=['',6.072528,-1.1645572,6.47];
ybs[8795]=['',6.0830952,-1.4098633,6.41];
ybs[8796]=['',6.0604653,1.1232147,6.21];
ybs[8797]=['88 Aqr',6.0682823,-0.3672054,3.66];
ybs[8798]=['',6.0696503,-0.4879136,5.87];
ybs[8799]=['',6.0707729,-0.7457322,5.81];
ybs[8800]=['57 Peg',6.0683466,0.1537706,5.12];
ybs[8801]=['',6.0698733,-0.2509318,6.42];
ybs[8802]=['89 Aqr',6.0703299,-0.3896319,4.69];
ybs[8803]=['',6.071644,-0.7061322,5.83];
ybs[8804]=['π Cep',6.0592754,1.3180803,4.41];
ybs[8805]=['ι Gru',6.0725768,-0.7873768,3.9];
ybs[8806]=['58 Peg',6.0705185,0.1737509,5.39];
ybs[8807]=['2 Cas',6.0684752,1.0378819,5.7];
ybs[8808]=['',6.0741698,-0.5129812,6.51];
ybs[8809]=['',6.0734363,0.3094081,5.71];
ybs[8810]=['6 And',6.0719953,0.7623155,5.94];
ybs[8811]=['59 Peg',6.0779999,0.154522,5.16];
ybs[8812]=['60 Peg',6.0781903,0.4709019,6.17];
ybs[8813]=['',6.0852364,-0.8636803,6.8];
ybs[8814]=['',6.0893318,-1.0919869,6.12];
ybs[8815]=['7 And',6.0810589,0.8646351,4.52];
ybs[8816]=['',6.0836085,0.5161861,6.35];
ybs[8817]=['',6.0840842,1.000108,5.56];
ybs[8818]=['',6.0854221,0.1954536,5.82];
ybs[8819]=['φ Aqr',6.089411,-0.1032382,4.22];
ybs[8820]=['',6.0926115,-0.7150912,5.77];
ybs[8821]=['',6.0909614,-0.184216,6.12];
ybs[8822]=['',6.0884218,0.8857813,6.31];
ybs[8823]=['',6.0892629,0.5219484,6.41];
ybs[8824]=['',6.0903984,0.423013,6.36];
ybs[8825]=['',6.0948362,-0.0586863,5.55];
ybs[8826]=['ψ1 Aqr',6.0962767,-0.1562737,4.21];
ybs[8827]=['61 Peg',6.0954391,0.4953543,6.49];
ybs[8828]=['',6.1016939,-1.0797834,5.66];
ybs[8829]=['',6.0890142,1.2979123,5.84];
ybs[8830]=['',6.0963179,0.4346754,6.6];
ybs[8831]=['',6.099997,-0.7741429,5.92];
ybs[8832]=['',6.1006829,-0.7166389,6.47];
ybs[8833]=['γ Tuc',6.1036142,-1.014066,3.99];
ybs[8834]=['',6.1125664,-1.3847172,6.33];
ybs[8835]=['χ Aqr',6.1004387,-0.1325161,5.06];
ybs[8836]=['',6.0937511,1.2395671,5.56];
ybs[8837]=['γ Psc',6.1017359,0.0596259,3.69];
ybs[8838]=['',6.0991626,0.931092,5.54];
ybs[8839]=['',6.0977905,1.0837981,6.53];
ybs[8840]=['',6.1078752,-1.1752503,6.13];
ybs[8841]=['',6.1040391,-0.2020901,6.34];
ybs[8842]=['',6.1018098,0.7906038,6.43];
ybs[8843]=['ψ2 Aqr',6.1050517,-0.1579232,4.39];
ybs[8844]=['φ Gru',6.1064913,-0.7101787,5.53];
ybs[8845]=['8 And',6.1037953,0.8578193,4.85];
ybs[8846]=['',6.1046841,0.7962726,6.48];
ybs[8847]=['τ Oct',6.1563046,-1.5230216,5.49];
ybs[8848]=['γ Scl',6.1092636,-0.565446,4.41];
ybs[8849]=['9 And',6.1067269,0.7314296,6.02];
ybs[8850]=['ψ3 Aqr',6.109671,-0.165397,4.98];
ybs[8851]=['94 Aqr',6.1103544,-0.2325579,5.08];
ybs[8852]=['',6.1007447,1.3165586,6.38];
ybs[8853]=['96 Aqr',6.11155,-0.0870939,5.55];
ybs[8854]=['',6.111655,-0.3131286,5.93];
ybs[8855]=['',6.1095224,0.7901368,6.5];
ybs[8856]=['',6.1131822,-0.5859715,6.37];
ybs[8857]=['ο Cep',6.1070772,1.1911156,4.75];
ybs[8858]=['',6.111473,0.6096027,6.32];
ybs[8859]=['11 And',6.1114656,0.8510157,5.44];
ybs[8860]=['',6.11233,0.8467497,6.32];
ybs[8861]=['10 And',6.1132089,0.7367458,5.79];
ybs[8862]=['',6.1182457,-0.8756699,6.05];
ybs[8863]=['7 Psc',6.115591,0.0962691,5.05];
ybs[8864]=['',6.1171466,-0.1007683,6.17];
ybs[8865]=['τ Peg',6.1167327,0.4166926,4.6];
ybs[8866]=['',6.1144077,1.0839263,6.45];
ybs[8867]=['63 Peg',6.1175037,0.5331888,5.59];
ybs[8868]=['',6.1198219,-0.4686584,5.64];
ybs[8869]=['',6.1169426,0.772323,6.13];
ybs[8870]=['12 And',6.1176936,0.6687525,5.77];
ybs[8871]=['',6.1158827,1.0881691,6.39];
ybs[8872]=['64 Peg',6.122246,0.5575817,5.32];
ybs[8873]=['',6.1225327,0.4667616,6.62];
ybs[8874]=['',6.1276469,-1.0458213,6.09];
ybs[8875]=['97 Aqr',6.1258103,-0.2601329,5.2];
ybs[8876]=['65 Peg',6.1256584,0.3658779,6.29];
ybs[8877]=['98 Aqr',6.1272278,-0.3484703,3.97];
ybs[8878]=['66 Peg',6.1274735,0.2172687,5.08];
ybs[8879]=['',6.124555,1.0518793,5.56];
ybs[8880]=['',6.1316697,-0.9367804,6.15];
ybs[8881]=['',6.1308668,-0.7503115,6.1];
ybs[8882]=['',6.1295389,0.0074372,6.31];
ybs[8883]=['',6.1329992,-0.9033229,5.75];
ybs[8884]=['',6.1304401,0.5701317,6.69];
ybs[8885]=['',6.1322747,-0.3238059,6.19];
ybs[8886]=['',6.137896,-0.9898506,5.59];
ybs[8887]=['',6.1338147,0.7199064,6.72];
ybs[8888]=['67 Peg',6.1350517,0.5675785,5.57];
ybs[8889]=['4 Cas',6.134559,1.089393,4.98];
ybs[8890]=['υ Peg',6.1374533,0.4108342,4.4];
ybs[8891]=['99 Aqr',6.1406444,-0.3579144,4.39];
ybs[8892]=['ο Gru',6.1434167,-0.9178103,5.52];
ybs[8893]=['',6.1459646,-1.1597024,6.45];
ybs[8894]=['',6.1463053,-1.0182434,5.63];
ybs[8895]=['',6.1457347,-0.8730515,6.2];
ybs[8896]=['κ Psc',6.1443761,0.0242704,4.94];
ybs[8897]=['9 Psc',6.145744,0.0219485,6.25];
ybs[8898]=['13 And',6.1449014,0.7513117,5.75];
ybs[8899]=['',6.1493202,-0.6180092,6.32];
ybs[8900]=['69 Peg',6.1474579,0.4416087,5.98];
ybs[8901]=['θ Psc',6.1488631,0.1136909,4.28];
ybs[8902]=['',6.1494813,-0.197477,6.37];
ybs[8903]=['',6.1449441,1.2303659,5.6];
ybs[8904]=['',6.1540908,-1.099132,5.68];
ybs[8905]=['',6.1537674,-0.7742729,6.43];
ybs[8906]=['',6.1534958,-0.1593644,6.18];
ybs[8907]=['',6.1536746,0.4046194,6.35];
ybs[8908]=['70 Peg',6.1540063,0.2250736,4.55];
ybs[8909]=['',6.1557599,-0.0767515,6.25];
ybs[8910]=['',6.1579513,0.8598948,6.17];
ybs[8911]=['',6.1574109,1.0242319,4.91];
ybs[8912]=['',6.160413,0.6771401,6.05];
ybs[8913]=['',6.1622399,-0.1073897,6.39];
ybs[8914]=['',6.1643782,-0.7803057,6.02];
ybs[8915]=['14 And',6.1631549,0.6871669,5.22];
ybs[8916]=['',6.164439,-0.0689725,6.49];
ybs[8917]=['100 Aqr',6.1653013,-0.3706039,6.29];
ybs[8918]=['',6.1651068,0.4980997,6.41];
ybs[8919]=['13 Psc',6.1663218,-0.0165878,6.38];
ybs[8920]=['',6.1734661,-1.3482625,5.81];
ybs[8921]=['',6.1680843,0.6124003,6.65];
ybs[8922]=['β Scl',6.1709427,-0.6576896,4.37];
ybs[8923]=['',6.1372331,1.5244741,5.58];
ybs[8924]=['101 Aqr',6.1721648,-0.3626607,4.71];
ybs[8925]=['71 Peg',6.1727865,0.3950451,5.32];
ybs[8926]=['',6.1736924,0.788777,6.24];
ybs[8927]=['',6.1747899,0.3661071,6.06];
ybs[8928]=['72 Peg',6.1748534,0.5490951,4.98];
ybs[8929]=['14 Psc',6.1758778,-0.0194068,5.87];
ybs[8930]=['',6.1810316,-1.1266764,7.4];
ybs[8931]=['',6.1788795,-0.2637231,5.96];
ybs[8932]=['15 And',6.1777272,0.7046242,5.59];
ybs[8933]=['73 Peg',6.1778256,0.5870035,5.63];
ybs[8934]=['ι Phe',6.1801468,-0.7414049,4.71];
ybs[8935]=['',6.178414,0.666009,6.18];
ybs[8936]=['',6.1819477,-0.1279114,6.39];
ybs[8937]=['',6.1787408,1.2527597,5.84];
ybs[8938]=['',6.1835349,0.4310404,6.45];
ybs[8939]=['16 Psc',6.18563,0.0390593,5.68];
ybs[8940]=['',6.1860087,0.5766547,6.35];
ybs[8941]=['',6.1888438,-0.5538816,6.52];
ybs[8942]=['',6.1953467,-1.3392639,6];
ybs[8943]=['',6.1912397,-0.2255749,5.65];
ybs[8944]=['',6.192243,-0.7916238,4.74];
ybs[8945]=['74 Peg',6.1911338,0.2960312,6.26];
ybs[8946]=['λ And',6.1905344,0.8132158,3.82];
ybs[8947]=['',6.1904122,0.777805,5.8];
ybs[8948]=['75 Peg',6.1923643,0.3235205,5.53];
ybs[8949]=['',6.1923372,0.8087075,6.58];
ybs[8950]=['ι And',6.1930595,0.7575403,4.29];
ybs[8951]=['θ Phe',6.1992886,-0.8116112,6.09];
ybs[8952]=['18 And',6.1973828,0.8832681,5.3];
ybs[8953]=['ω1 Aqr',6.200511,-0.245843,5];
ybs[8954]=['ι Psc',6.2011596,0.100571,4.13];
ybs[8955]=['',6.2010069,0.1712714,5.97];
ybs[8956]=['',6.196985,1.3164782,5.95];
ybs[8957]=['',6.1978383,1.2939635,5.98];
ybs[8958]=['',6.2014463,0.6595322,6.53];
ybs[8959]=['γ Cep',6.1975891,1.3573141,3.21];
ybs[8960]=['μ Scl',6.2042972,-0.5574078,5.31];
ybs[8961]=['κ And',6.2029866,0.7761447,4.14];
ybs[8962]=['',6.2042052,0.6432721,6.23];
ybs[8963]=['',6.2063549,-0.4193034,6.6];
ybs[8964]=['',6.2064473,-0.2014911,5.89];
ybs[8965]=['103 Aqr',6.20833,-0.312261,5.34];
ybs[8966]=['',6.2074987,0.8665246,6.26];
ybs[8967]=['104 Aqr',6.2091503,-0.308581,4.82];
ybs[8968]=['',6.2098613,0.1289198,5.89];
ybs[8969]=['λ Psc',6.2103223,0.0334407,4.5];
ybs[8970]=['',6.2094414,1.0017491,6.24];
ybs[8971]=['',6.2110225,0.7876315,6.57];
ybs[8972]=['',6.2122017,-0.2672405,5.28];
ybs[8973]=['ω2 Aqr',6.213318,-0.2514838,4.49];
ybs[8974]=['',6.2112644,1.1283829,6.56];
ybs[8975]=['',6.2120892,1.0788835,6.4];
ybs[8976]=['77 Peg',6.2160864,0.1826916,5.06];
ybs[8977]=['',6.2181315,-0.2643888,6.36];
ybs[8978]=['',6.2191121,-0.7844773,6.09];
ybs[8979]=['',6.2211332,-1.2279118,6.07];
ybs[8980]=['',6.2225713,-1.3727935,5.75];
ybs[8981]=['',6.2200508,-1.1216942,5.72];
ybs[8982]=['78 Peg',6.2187247,0.514833,4.93];
ybs[8983]=['106 Aqr',6.2197829,-0.3166175,5.24];
ybs[8984]=['',6.2210288,-0.4557103,6.17];
ybs[8985]=['',6.2221402,0.9762647,6.51];
ybs[8986]=['',6.2277817,-0.6989404,6.31];
ybs[8987]=['107 Aqr',6.2276899,-0.3236171,5.29];
ybs[8988]=['ψ And',6.2275845,0.8125632,4.95];
ybs[8989]=['19 Psc',6.2292785,0.0632306,5.04];
ybs[8990]=['',6.229936,1.1679466,5.95];
ybs[8991]=['σ Phe',6.2332551,-0.8742434,5.18];
ybs[8992]=['',6.2339477,-1.1913261,6.89];
ybs[8993]=['τ Cas',6.2319848,1.0260466,4.87];
ybs[8994]=['',6.2331238,-0.205506,5.73];
ybs[8995]=['',6.2318777,1.005093,5.51];
ybs[8996]=['',6.2342242,0.8197587,6.07];
ybs[8997]=['20 Psc',6.2360555,-0.0458226,5.49];
ybs[8998]=['',6.2356328,1.185832,5.04];
ybs[8999]=['',6.2386804,-0.1089838,6.07];
ybs[9000]=['',6.2398853,0.0410226,6.46];
ybs[9001]=['δ Scl',6.2404056,-0.4885879,4.57];
ybs[9002]=['',6.2389019,1.1346846,6.41];
ybs[9003]=['6 Cas',6.2397464,1.088225,5.43];
ybs[9004]=['',6.240035,1.0492072,6.34];
ybs[9005]=['',6.2413641,1.0314777,6.33];
ybs[9006]=['',6.2429947,-0.2744502,6.24];
ybs[9007]=['21 Psc',6.2426656,0.0211601,5.77];
ybs[9008]=['',6.2441123,-1.0943767,6.59];
ybs[9009]=['',6.2435749,0.6381194,5.9];
ybs[9010]=['79 Peg',6.2434773,0.505775,5.97];
ybs[9011]=['',6.2443216,-0.4397376,6.42];
ybs[9012]=['',6.2461167,-0.1717033,5.94];
ybs[9013]=['',6.2465343,0.9033467,6.44];
ybs[9014]=['',6.2474761,-0.2489825,5.72];
ybs[9015]=['80 Peg',6.2509224,0.1649274,5.79];
ybs[9016]=['108 Aqr',6.2509718,-0.3276433,5.18];
ybs[9017]=['γ1 Oct',6.2547522,-1.4291203,5.11];
ybs[9018]=['22 Psc',6.253601,0.0535223,5.55];
ybs[9019]=['',6.2532292,1.3567451,6.55];
ybs[9020]=['',6.2554305,0.3806069,6.11];
ybs[9021]=['φ Peg',6.2558636,0.3360913,5.08];
ybs[9022]=['',6.2559575,-0.2463493,5.87];
ybs[9023]=['',6.2553002,1.3208836,6.39];
ybs[9024]=['82 Peg',6.2564421,0.1934494,5.3];
ybs[9025]=['',6.2574409,-0.1546419,5.75];
ybs[9026]=['24 Psc',6.2578047,-0.0526952,5.93];
ybs[9027]=['25 Psc',6.2584677,0.0388668,6.28];
ybs[9028]=['',6.2596618,-0.4204989,6.24];
ybs[9029]=['',6.2640671,-0.4695958,6.35];
ybs[9030]=['ρ Cas',6.2640753,1.0059346,4.54];
ybs[9031]=['',6.2653329,-0.7009876,6.03];
ybs[9032]=['',6.2658731,0.0042854,5.61];
ybs[9033]=['26 Psc',6.2674094,0.1257944,6.21];
ybs[9034]=['',6.2680813,-0.554758,6.1];
ybs[9035]=['',6.2680813,-0.5541035,6.83];
ybs[9036]=['',6.2685001,0.4553804,6.54];
ybs[9037]=['',6.2692362,1.0044126,6];
ybs[9038]=['',6.2692481,0.8288955,6];
ybs[9039]=['',6.2733991,-0.4293656,6.31];
ybs[9040]=['',6.2742174,0.3976635,6.15];
ybs[9041]=['',6.272967,1.4543392,6.59];
ybs[9042]=['',6.2758139,0.7469089,5.97];
ybs[9043]=['',6.2761895,-0.4622892,6.26];
ybs[9044]=['',6.2761584,0.9746307,5.55];
ybs[9045]=['',6.2770628,-1.0964158,5.97];
ybs[9046]=['γ2 Oct',6.278088,-1.4317565,5.73];
ybs[9047]=['η Tuc',6.2781726,-1.1198371,5];
ybs[9048]=['',6.2779777,1.0499901,6.47];
ybs[9049]=['ψ Peg',6.2788753,0.4411805,4.66];
ybs[9050]=['1 Cet',6.2814824,-0.2742105,6.26];
ybs[9051]=['',6.2817287,0.8992811,4.8];
ybs[9052]=['27 Psc',6.2828761,-0.0596852,4.86];
ybs[9053]=['',0.0003267,0.5675473,6.52];
ybs[9054]=['π Phe',0.0008163,-0.9182079,5.13];
ybs[9055]=['',0.0001274,0.8124412,6.54];
ybs[9056]=['σ Cas',0.0011464,0.9754889,4.88];
ybs[9057]=['ω Psc',0.0024755,0.1221684,4.01];
ybs[9058]=['',0.0031441,-0.5122297,5.62];
ybs[9059]=['',0.00324,0.5909832,6.58];
ybs[9060]=['',0.00324,0.5909832,6.58];
ybs[9061]=['ε Tuc',0.005103,-1.1421578,4.5];
ybs[9062]=['',0.0068658,-0.7706354,6.29];
ybs[9063]=['',0.0072226,0.4721941,6.46];
ybs[9064]=['',0.0077464,1.0418939,6.19];
ybs[9065]=['',0.0086713,0.7922002,6.38];
ybs[9066]=['τ Phe',0.0101495,-0.8495147,5.71];
ybs[9067]=['',0.0112724,-0.8761698,5.53];
ybs[9068]=['',0.0112701,0.8747251,6.22];
ybs[9069]=['θ Oct',0.0123391,-1.342672,4.78];
ybs[9070]=['',0.0125685,1.0709244,5.55];
ybs[9071]=['',0.0130499,0.741828,6.25];
ybs[9072]=['29 Psc',0.0134335,-0.0504594,5.1];
ybs[9073]=['85 Peg',0.0149623,0.4750495,5.75];
ybs[9074]=['30 Psc',0.0140283,-0.1025866,4.41];
ybs[9075]=['',0.0147291,-0.2537661,7.1];
ybs[9076]=['ζ Scl',0.0156353,-0.5163363,5.01];
ybs[9077]=['31 Psc',0.0159705,0.1587086,6.32];
ybs[9078]=['32 Psc',0.0163704,0.1504813,5.63];
ybs[9079]=['',0.0169071,1.1560235,5.86];
ybs[9080]=['',0.0183809,-0.3474903,6.25];
ybs[9081]=['',0.0191115,-0.4190343,6.44];
ybs[9082]=['',0.0205234,1.1131417,6.24];
ybs[9083]=['2 Cet',0.0217837,-0.3001921,4.55];
ybs[9084]=['',0.0224523,1.166728,6.29];
ybs[9085]=['9 Cas',0.0240178,1.0895068,5.88];
ybs[9086]=['',0.0243497,-0.2861035,5.78];
ybs[9087]=['',0.0243798,-0.5084536,6.4];
ybs[9088]=['3 Cet',0.0251114,-0.1810444,4.94];
ybs[9089]=['',0.0261189,1.1746594,5.67];
ybs[9090]=['',0.0256506,0.7370279,6.01];
ybs[9091]=['',0.0249928,-1.2699262,7.31];
ybs[9092]=['',0.0268845,0.6073062,6.12];
ybs[9093]=['',0.0257725,-1.24443,5.59];
ybs[9094]=['',0.0270338,0.4674908,6.25];
ybs[9095]=['',0.0278559,1.0725139,5.8];    
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