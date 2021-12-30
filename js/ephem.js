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
      var FUDGE_FACTOR_FOR_MISSING_RINGS =  -0.8; //gives the correct answer near opposition, 2022-07
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
      epoch: when("UT 2022-01-21"),
      a: 2.766043062222408,
      e: 0.07850100198908602,
      i: rads(10.58769305845201),
      Ω: rads(80.26859547732911),
      ω: rads(73.63703979153577),
      M0: rads(291.3755993017663),
      n: rads(0.2142474533574742) 
    }
  );
  
  var vesta = build_minor_planet('Vesta', 4, 3.4, 0.32, {
      equinox: when_j2000, 
      epoch: when("UT 2022-01-21"),
      a: 2.361266458114362,
      e: 0.08823417531213737,
      i: rads(7.141717168552266),
      Ω: rads(103.8039247181175),
      ω: rads(151.0909385501822),
      M0: rads(7.03152246557125),
      n: rads(0.2716355310980331) 
    }
  );
  
  var eunomia = build_minor_planet('Eunomia', 15, 5.41, 0.23, {
      equinox: when_j2000, 
      epoch: when("UT 2022-01-21"),
      a: 2.644075438339401,
      e: 0.1866253417318819,
      i: rads(11.75199956073668),
      Ω: rads(292.9260977033457),
      ω: rads(98.63168561763565),
      M0: rads(152.5002318505929),
      n: rads(0.2292415367124858) 
    }
  );
  
  var euterpe = build_minor_planet('Euterpe', 27, 7.08, 0.15, {
      equinox: when_j2000, 
      epoch: when("UT 2022-01-21"),
      a: 2.34807608298257,
      e: 0.172072813394719,
      i: rads(1.58339470560508),
      Ω: rads(94.78287259702536),
      ω: rads(356.2569067954329),
      M0: rads(249.6942624999856),
      n: rads(0.273927621465773) 
    }
  );
  
  var melpomene = build_minor_planet('Melpomene', 18, 6.53, 0.25, {
      equinox: when_j2000, 
      epoch: when("UT 2022-01-21"),
      a: 2.295788850979223,
      e: 0.2179091998946986,
      i: rads(10.13249376222656),
      Ω: rads(150.3617349643433),
      ω: rads(228.119234422882),
      M0: rads(190.3739342479254),
      n: rads(0.2833388717102693) 
    }
  );
  
  var astraea = build_minor_planet('Astraea', 5, 6.99, 0.15, {
      equinox: when_j2000, 
      epoch: when("UT 2022-01-21"),
      a: 2.575176611830832,
      e: 0.1900993645991224,
      i: rads(5.367623287373819),
      Ω: rads(141.5703596793466),
      ω: rads(358.7403904206299),
      M0: rads(160.9820880322897),
      n: rads(0.2385028345272263) 
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
     {when:'UT 2022-01-31 00',text:'Mercury 7.57°N of Moon'},
     {when:'UT 2022-02-01 06',text:'New Moon 365182.588 km'},
     {when:'UT 2022-02-01 09',text:'Saturn 4.21°N of Moon'},
     {when:'UT 2022-02-02 21',text:'Jupiter 4.32°N of Moon'},
     {when:'UT 2022-02-03 21',text:'Neptune 3.86°N of Moon'},
     {when:'UT 2022-02-04 19',text:'Saturn 0.85°S of Sun'},
     {when:'UT 2022-02-07 20',text:'Uranus 1.18°N of Moon'},
     {when:'UT 2022-02-08 14',text:'First Quarter 400920.807 km'},
     {when:'UT 2022-02-10 09',text:'Aldebaran 6.68°S of Moon'},
     {when:'UT 2022-02-11 03',text:'Moon at apogee 404896.549 km'},
     {when:'UT 2022-02-13 01',text:'Mars 6.58°S of Venus'},
     {when:'UT 2022-02-13 23',text:'Pollux 2.55°N of Moon'},
     {when:'UT 2022-02-16 17',text:'Full Moon 391887.790 km'},
     {when:'UT 2022-02-16 18',text:'Regulus 4.84°S of Moon'},
     {when:'UT 2022-02-16 21',text:'Mercury at greatest elongation 26.3° West'},
     {when:'UT 2022-02-20 20',text:'Spica 5.25°S of Moon'},
     {when:'UT 2022-02-23 23',text:'Last Quarter 371026.934 km'},
     {when:'UT 2022-02-24 06',text:'Antares 3.44°S of Moon'},
     {when:'UT 2022-02-26 22',text:'Moon at perigee 367788.952 km'},
     {when:'UT 2022-02-27 06',text:'Venus 8.74°N of Moon'},
     {when:'UT 2022-02-27 09',text:'Mars 3.52°N of Moon'},
     {when:'UT 2022-02-28 20',text:'Mercury 3.73°N of Moon'},
     {when:'UT 2022-03-01 00',text:'Saturn 4.30°N of Moon'},
     {when:'UT 2022-03-02 13',text:'Saturn 0.69°N of Mercury'},
     {when:'UT 2022-03-02 18',text:'New Moon 375031.477 km'},
     {when:'UT 2022-03-02 19',text:'Jupiter 4.14°N of Moon'},
     {when:'UT 2022-03-03 09',text:'Neptune 3.71°N of Moon'},
     {when:'UT 2022-03-05 14',text:'Jupiter 0.98°S of Sun'},
     {when:'UT 2022-03-07 06',text:'Uranus 0.84°N of Moon'},
     {when:'UT 2022-03-09 17',text:'Aldebaran 6.96°S of Moon'},
     {when:'UT 2022-03-10 11',text:'First Quarter 404095.540 km'},
     {when:'UT 2022-03-10 23',text:'Moon at apogee 404267.872 km'},
     {when:'UT 2022-03-12 14',text:'Mars 3.99°S of Venus'},
     {when:'UT 2022-03-13 08',text:'Pollux 2.36°N of Moon'},
     {when:'UT 2022-03-13 12',text:'Neptune 1.12°S of Sun'},
     {when:'UT 2022-03-16 02',text:'Regulus 4.91°S of Moon'},
     {when:'UT 2022-03-18 07',text:'Full Moon 380830.618 km'},
     {when:'UT 2022-03-20 02',text:'Spica 5.08°S of Moon'},
     {when:'UT 2022-03-20 09',text:'Venus at greatest elongation 46.6° West'},
     {when:'UT 2022-03-20 16',text:'Equinox'},
     {when:'UT 2022-03-20 22',text:'Jupiter 1.29°N of Mercury'},
     {when:'UT 2022-03-23 11',text:'Antares 3.20°S of Moon'},
     {when:'UT 2022-03-23 12',text:'Neptune 1.03°N of Mercury'},
     {when:'UT 2022-03-24 00',text:'Moon at perigee 369760.001 km'},
     {when:'UT 2022-03-25 06',text:'Last Quarter 370166.069 km'},
     {when:'UT 2022-03-28 03',text:'Mars 4.10°N of Moon'},
     {when:'UT 2022-03-28 10',text:'Venus 6.69°N of Moon'},
     {when:'UT 2022-03-28 12',text:'Saturn 4.43°N of Moon'},
     {when:'UT 2022-03-29 13',text:'Saturn 2.16°S of Venus'},
     {when:'UT 2022-03-30 15',text:'Jupiter 3.93°N of Moon'},
     {when:'UT 2022-03-30 19',text:'Neptune 3.69°N of Moon'},
     {when:'UT 2022-04-01 00',text:'Mercury 2.55°N of Moon'},
     {when:'UT 2022-04-01 06',text:'New Moon 386233.744 km'},
     {when:'UT 2022-04-02 23',text:'Mercury in superior conjunction 1.04° South'},
     {when:'UT 2022-04-03 17',text:'Uranus 0.57°N of Moon'},
     {when:'UT 2022-04-04 22',text:'Saturn 0.32°N of Mars'},
     {when:'UT 2022-04-06 01',text:'Aldebaran 7.16°S of Moon'},
     {when:'UT 2022-04-07 19',text:'Moon at apogee 404437.938 km'},
     {when:'UT 2022-04-09 07',text:'First Quarter 403029.927 km'},
     {when:'UT 2022-04-09 16',text:'Pollux 2.16°N of Moon'},
     {when:'UT 2022-04-12 11',text:'Regulus 5.05°S of Moon'},
     {when:'UT 2022-04-12 20',text:'Neptune 0.11°S of Jupiter'},
     {when:'UT 2022-04-16 11',text:'Spica 5.06°S of Moon'},
     {when:'UT 2022-04-16 19',text:'Full Moon 370264.512 km'},
     {when:'UT 2022-04-18 14',text:'Uranus 2.14°S of Mercury'},
     {when:'UT 2022-04-19 15',text:'Moon at perigee 365143.266 km'},
     {when:'UT 2022-04-19 18',text:'Antares 3.07°S of Moon'},
     {when:'UT 2022-04-23 12',text:'Last Quarter 372072.679 km'},
     {when:'UT 2022-04-24 21',text:'Saturn 4.51°N of Moon'},
     {when:'UT 2022-04-25 22',text:'Mars 3.90°N of Moon'},
     {when:'UT 2022-04-27 02',text:'Venus 3.78°N of Moon'},
     {when:'UT 2022-04-27 03',text:'Neptune 3.71°N of Moon'},
     {when:'UT 2022-04-27 08',text:'Jupiter 3.64°N of Moon'},
     {when:'UT 2022-04-27 19',text:'Neptune 0.01°N of Venus'},
     {when:'UT 2022-04-29 08',text:'Mercury at greatest elongation 20.6° East'},
     {when:'UT 2022-04-30 19',text:'Jupiter 0.25°N of Venus'},
     {when:'UT 2022-04-30 20',text:'ECLIPSE New Moon 396520.139 km'},
     {when:'UT 2022-05-01 04',text:'Uranus 0.40°N of Moon'},
     {when:'UT 2022-05-02 14',text:'Mercury 1.85°N of Moon'},
     {when:'UT 2022-05-03 09',text:'Aldebaran 7.23°S of Moon'},
     {when:'UT 2022-05-05 07',text:'Uranus 0.36°S of Sun'},
     {when:'UT 2022-05-05 13',text:'Moon at apogee 405285.349 km'},
     {when:'UT 2022-05-07 00',text:'Pollux 2.08°N of Moon'},
     {when:'UT 2022-05-09 00',text:'First Quarter 398275.651 km'},
     {when:'UT 2022-05-09 20',text:'Regulus 5.13°S of Moon'},
     {when:'UT 2022-05-13 21',text:'Spica 5.09°S of Moon'},
     {when:'UT 2022-05-16 04',text:'ECLIPSE Full Moon 362126.539 km'},
     {when:'UT 2022-05-17 03',text:'Antares 3.05°S of Moon'},
     {when:'UT 2022-05-17 15',text:'Moon at perigee 360298.232 km'},
     {when:'UT 2022-05-17 23',text:'Neptune 0.57°N of Mars'},
     {when:'UT 2022-05-21 19',text:'Mercury in inferior conjunction 1.24° South'},
     {when:'UT 2022-05-22 05',text:'Saturn 4.46°N of Moon'},
     {when:'UT 2022-05-22 19',text:'Last Quarter 376453.794 km'},
     {when:'UT 2022-05-24 10',text:'Neptune 3.67°N of Moon'},
     {when:'UT 2022-05-24 19',text:'Mars 2.78°N of Moon'},
     {when:'UT 2022-05-25 00',text:'Jupiter 3.25°N of Moon'},
     {when:'UT 2022-05-27 03',text:'Venus 0.20°N of Moon'},
     {when:'UT 2022-05-28 14',text:'Uranus 0.25°N of Moon'},
     {when:'UT 2022-05-29 00',text:'Jupiter 0.63°N of Mars'},
     {when:'UT 2022-05-29 13',text:'Mercury 3.72°S of Moon'},
     {when:'UT 2022-05-30 12',text:'New Moon 403791.204 km'},
     {when:'UT 2022-05-30 16',text:'Aldebaran 7.22°S of Moon'},
     {when:'UT 2022-06-02 01',text:'Moon at apogee 406191.987 km'},
     {when:'UT 2022-06-03 06',text:'Pollux 2.11°N of Moon'},
     {when:'UT 2022-06-06 03',text:'Regulus 5.07°S of Moon'},
     {when:'UT 2022-06-07 15',text:'First Quarter 391288.432 km'},
     {when:'UT 2022-06-10 07',text:'Spica 5.04°S of Moon'},
     {when:'UT 2022-06-11 13',text:'Uranus 1.61°N of Venus'},
     {when:'UT 2022-06-13 14',text:'Antares 3.05°S of Moon'},
     {when:'UT 2022-06-14 12',text:'Full Moon 357656.377 km'},
     {when:'UT 2022-06-14 23',text:'Moon at perigee 357432.471 km'},
     {when:'UT 2022-06-16 15',text:'Mercury at greatest elongation 23.2° West'},
     {when:'UT 2022-06-18 12',text:'Saturn 4.27°N of Moon'},
     {when:'UT 2022-06-20 17',text:'Neptune 3.52°N of Moon'},
     {when:'UT 2022-06-21 03',text:'Last Quarter 382826.040 km'},
     {when:'UT 2022-06-21 09',text:'Solstice'},
     {when:'UT 2022-06-21 14',text:'Jupiter 2.74°N of Moon'},
     {when:'UT 2022-06-22 18',text:'Mars 0.95°N of Moon'},
     {when:'UT 2022-06-23 14',text:'Aldebaran 3.02°S of Mercury'},
     {when:'UT 2022-06-24 22',text:'Uranus 0.05°N of Moon'},
     {when:'UT 2022-06-26 08',text:'Venus 2.69°S of Moon'},
     {when:'UT 2022-06-26 22',text:'Aldebaran 7.24°S of Moon'},
     {when:'UT 2022-06-27 08',text:'Mercury 3.94°S of Moon'},
     {when:'UT 2022-06-29 03',text:'New Moon 406573.819 km'},
     {when:'UT 2022-06-29 06',text:'Moon at apogee 406580.103 km'},
     {when:'UT 2022-06-30 12',text:'Pollux 2.19°N of Moon'},
     {when:'UT 2022-07-02 00',text:'Aldebaran 4.17°S of Venus'},
     {when:'UT 2022-07-03 09',text:'Regulus 4.92°S of Moon'},
     {when:'UT 2022-07-04 07',text:'Earth at aphelion 1.016715374 AU'},
     {when:'UT 2022-07-07 02',text:'First Quarter 383719.795 km'},
     {when:'UT 2022-07-07 16',text:'Spica 4.86°S of Moon'},
     {when:'UT 2022-07-11 00',text:'Antares 2.97°S of Moon'},
     {when:'UT 2022-07-13 09',text:'Moon at perigee 357263.701 km'},
     {when:'UT 2022-07-13 19',text:'Full Moon 357417.822 km'},
     {when:'UT 2022-07-15 20',text:'Saturn 4.05°N of Moon'},
     {when:'UT 2022-07-16 20',text:'Mercury in superior conjunction 1.52° North'},
     {when:'UT 2022-07-16 22',text:'Pollux 5.24°N of Mercury'},
     {when:'UT 2022-07-18 01',text:'Neptune 3.29°N of Moon'},
     {when:'UT 2022-07-19 01',text:'Jupiter 2.22°N of Moon'},
     {when:'UT 2022-07-20 14',text:'Last Quarter 390331.133 km'},
     {when:'UT 2022-07-21 17',text:'Mars 1.06°S of Moon'},
     {when:'UT 2022-07-22 06',text:'Uranus 0.24°S of Moon'},
     {when:'UT 2022-07-24 04',text:'Aldebaran 7.37°S of Moon'},
     {when:'UT 2022-07-26 10',text:'Moon at apogee 406274.499 km'},
     {when:'UT 2022-07-26 14',text:'Venus 4.17°S of Moon'},
     {when:'UT 2022-07-27 18',text:'Pollux 2.20°N of Moon'},
     {when:'UT 2022-07-28 18',text:'New Moon 404343.797 km'},
     {when:'UT 2022-07-29 21',text:'Mercury 3.60°S of Moon'},
     {when:'UT 2022-07-30 15',text:'Regulus 4.78°S of Moon'},
     {when:'UT 2022-08-01 09',text:'Uranus 1.37°N of Mars'},
     {when:'UT 2022-08-03 22',text:'Spica 4.59°S of Moon'},
     {when:'UT 2022-08-04 05',text:'Regulus 0.74°S of Mercury'},
     {when:'UT 2022-08-05 11',text:'First Quarter 377031.067 km'},
     {when:'UT 2022-08-07 09',text:'Antares 2.77°S of Moon'},
     {when:'UT 2022-08-07 10',text:'Pollux 6.56°N of Venus'},
     {when:'UT 2022-08-10 17',text:'Moon at perigee 359827.980 km'},
     {when:'UT 2022-08-12 02',text:'Full Moon 361411.858 km'},
     {when:'UT 2022-08-12 04',text:'Saturn 3.91°N of Moon'},
     {when:'UT 2022-08-14 10',text:'Neptune 3.09°N of Moon'},
     {when:'UT 2022-08-14 17',text:'Saturn at opposition'},
     {when:'UT 2022-08-15 10',text:'Jupiter 1.86°N of Moon'},
     {when:'UT 2022-08-18 15',text:'Uranus 0.55°S of Moon'},
     {when:'UT 2022-08-19 05',text:'Last Quarter 397557.125 km'},
     {when:'UT 2022-08-19 12',text:'Mars 2.68°S of Moon'},
     {when:'UT 2022-08-20 10',text:'Aldebaran 7.60°S of Moon'},
     {when:'UT 2022-08-22 22',text:'Moon at apogee 405417.800 km'},
     {when:'UT 2022-08-24 01',text:'Pollux 2.09°N of Moon'},
     {when:'UT 2022-08-25 21',text:'Venus 4.29°S of Moon'},
     {when:'UT 2022-08-26 21',text:'Regulus 4.74°S of Moon'},
     {when:'UT 2022-08-27 08',text:'New Moon 397575.973 km'},
     {when:'UT 2022-08-27 16',text:'Mercury at greatest elongation 27.3° East'},
     {when:'UT 2022-08-29 11',text:'Mercury 6.64°S of Moon'},
     {when:'UT 2022-08-31 04',text:'Spica 4.35°S of Moon'},
     {when:'UT 2022-09-03 15',text:'Antares 2.52°S of Moon'},
     {when:'UT 2022-09-03 18',text:'First Quarter 372300.486 km'},
     {when:'UT 2022-09-05 01',text:'Regulus 0.78°S of Venus'},
     {when:'UT 2022-09-07 18',text:'Moon at perigee 364492.051 km'},
     {when:'UT 2022-09-08 11',text:'Saturn 3.94°N of Moon'},
     {when:'UT 2022-09-09 01',text:'Aldebaran 4.34°S of Mars'},
     {when:'UT 2022-09-10 10',text:'Full Moon 369131.268 km'},
     {when:'UT 2022-09-10 19',text:'Neptune 3.03°N of Moon'},
     {when:'UT 2022-09-11 15',text:'Jupiter 1.81°N of Moon'},
     {when:'UT 2022-09-14 23',text:'Uranus 0.79°S of Moon'},
     {when:'UT 2022-09-16 18',text:'Aldebaran 7.84°S of Moon'},
     {when:'UT 2022-09-16 22',text:'Neptune at opposition'},
     {when:'UT 2022-09-17 02',text:'Mars 3.61°S of Moon'},
     {when:'UT 2022-09-17 22',text:'Last Quarter 402716.607 km'},
     {when:'UT 2022-09-19 15',text:'Moon at apogee 404555.725 km'},
     {when:'UT 2022-09-20 08',text:'Pollux 1.91°N of Moon'},
     {when:'UT 2022-09-23 01',text:'Equinox'},
     {when:'UT 2022-09-23 05',text:'Regulus 4.82°S of Moon'},
     {when:'UT 2022-09-23 07',text:'Mercury in inferior conjunction 2.86° South'},
     {when:'UT 2022-08-27 08',text:'New Moon 397575.971 km'},
     {when:'UT 2022-09-25 05',text:'Venus 2.75°S of Moon'},
     {when:'UT 2022-09-25 08',text:'Mercury 6.66°S of Moon'},
     {when:'UT 2022-09-26 01',text:'Venus 3.76°N of Mercury'},
     {when:'UT 2022-09-26 20',text:'Jupiter at opposition'},
     {when:'UT 2022-09-27 10',text:'Spica 4.23°S of Moon'},
     {when:'UT 2022-09-30 21',text:'Antares 2.33°S of Moon'},
     {when:'UT 2022-10-03 00',text:'First Quarter 370124.856 km'},
     {when:'UT 2022-10-04 17',text:'Moon at perigee 369325.025 km'},
     {when:'UT 2022-10-05 16',text:'Saturn 4.08°N of Moon'},
     {when:'UT 2022-10-08 03',text:'Neptune 3.10°N of Moon'},
     {when:'UT 2022-10-08 18',text:'Jupiter 2.07°N of Moon'},
     {when:'UT 2022-10-08 21',text:'Mercury at greatest elongation 18.0° West'},
     {when:'UT 2022-10-09 21',text:'Full Moon 379482.926 km'},
     {when:'UT 2022-10-12 07',text:'Uranus 0.85°S of Moon'},
     {when:'UT 2022-10-14 03',text:'Aldebaran 7.99°S of Moon'},
     {when:'UT 2022-10-15 05',text:'Mars 3.62°S of Moon'},
     {when:'UT 2022-10-17 10',text:'Moon at apogee 404327.832 km'},
     {when:'UT 2022-10-17 15',text:'Spica 3.48°S of Venus'},
     {when:'UT 2022-10-17 16',text:'Pollux 1.76°N of Moon'},
     {when:'UT 2022-10-17 17',text:'Last Quarter 404273.295 km'},
     {when:'UT 2022-10-20 13',text:'Regulus 4.93°S of Moon'},
     {when:'UT 2022-10-22 21',text:'Venus in superior conjunction 1.05° North'},
     {when:'UT 2022-10-24 16',text:'Mercury 0.39°S of Moon'},
     {when:'UT 2022-10-24 18',text:'Spica 4.23°S of Moon'},
     {when:'UT 2022-10-25 11',text:'ECLIPSE New Moon 376366.905 km'},
     {when:'UT 2022-10-25 12',text:'Venus 0.00°N of Moon'},
     {when:'UT 2022-10-28 03',text:'Antares 2.25°S of Moon'},
     {when:'UT 2022-10-29 15',text:'Moon at perigee 368290.502 km'},
     {when:'UT 2022-11-01 07',text:'First Quarter 370690.904 km'},
     {when:'UT 2022-11-01 21',text:'Saturn 4.19°N of Moon'},
     {when:'UT 2022-11-04 08',text:'Neptune 3.20°N of Moon'},
     {when:'UT 2022-11-04 20',text:'Jupiter 2.39°N of Moon'},
     {when:'UT 2022-11-08 11',text:'ECLIPSE Full Moon 390659.790 km'},
     {when:'UT 2022-11-08 13',text:'Uranus 0.75°S of Moon'},
     {when:'UT 2022-11-08 17',text:'Mercury in superior conjunction 0.09° North'},
     {when:'UT 2022-11-09 08',text:'Uranus at opposition'},
     {when:'UT 2022-11-10 11',text:'Aldebaran 8.02°S of Moon'},
     {when:'UT 2022-11-11 14',text:'Mars 2.47°S of Moon'},
     {when:'UT 2022-11-14 00',text:'Pollux 1.74°N of Moon'},
     {when:'UT 2022-11-14 07',text:'Moon at apogee 404921.387 km'},
     {when:'UT 2022-11-16 13',text:'Last Quarter 401640.891 km'},
     {when:'UT 2022-11-16 22',text:'Regulus 4.95°S of Moon'},
     {when:'UT 2022-11-21 04',text:'Spica 4.23°S of Moon'},
     {when:'UT 2022-11-22 17',text:'Venus 1.35°N of Mercury'},
     {when:'UT 2022-11-23 13',text:'Antares 3.16°S of Mercury'},
     {when:'UT 2022-11-23 18',text:'Antares 4.54°S of Venus'},
     {when:'UT 2022-11-23 23',text:'New Moon 366158.439 km'},
     {when:'UT 2022-11-24 12',text:'Antares 2.26°S of Moon'},
     {when:'UT 2022-11-24 14',text:'Venus 2.33°N of Moon'},
     {when:'UT 2022-11-24 15',text:'Mercury 0.94°N of Moon'},
     {when:'UT 2022-11-26 02',text:'Moon at perigee 362826.006 km'},
     {when:'UT 2022-11-29 05',text:'Saturn 4.17°N of Moon'},
     {when:'UT 2022-11-30 15',text:'First Quarter 373976.351 km'},
     {when:'UT 2022-12-01 13',text:'Neptune 3.18°N of Moon'},
     {when:'UT 2022-12-02 01',text:'Jupiter 2.51°N of Moon'},
     {when:'UT 2022-12-05 18',text:'Uranus 0.65°S of Moon'},
     {when:'UT 2022-12-07 19',text:'Aldebaran 7.99°S of Moon'},
     {when:'UT 2022-12-08 04',text:'Full Moon 400238.742 km'},
     {when:'UT 2022-12-08 06',text:'Mars at opposition'},
     {when:'UT 2022-12-11 08',text:'Pollux 1.83°N of Moon'},
     {when:'UT 2022-12-12 00',text:'Moon at apogee 405868.535 km'},
     {when:'UT 2022-12-14 06',text:'Regulus 4.81°S of Moon'},
     {when:'UT 2022-12-16 09',text:'Last Quarter 395489.753 km'},
     {when:'UT 2022-12-18 14',text:'Spica 4.11°S of Moon'},
     {when:'UT 2022-12-21 16',text:'Mercury at greatest elongation 20.1° East'},
     {when:'UT 2022-12-21 22',text:'Solstice'},
     {when:'UT 2022-12-21 23',text:'Antares 2.25°S of Moon'},
     {when:'UT 2022-12-22 04',text:'Aldebaran 8.24°S of Mars'},
     {when:'UT 2022-12-23 10',text:'New Moon 359080.704 km'},
     {when:'UT 2022-12-24 08',text:'Moon at perigee 358270.156 km'},
     {when:'UT 2022-12-24 11',text:'Venus 3.47°N of Moon'},
     {when:'UT 2022-12-24 19',text:'Mercury 3.76°N of Moon'},
     {when:'UT 2022-12-26 16',text:'Saturn 4.02°N of Moon'},
     {when:'UT 2022-12-28 20',text:'Neptune 2.98°N of Moon'},
     {when:'UT 2022-12-29 09',text:'Venus 1.41°S of Mercury'},
     {when:'UT 2022-12-29 11',text:'Jupiter 2.30°N of Moon'},
     {when:'UT 2022-12-30 01',text:'First Quarter 379784.405 km'},
     {when:'UT 2023-01-01 22',text:'Uranus 0.71°S of Moon'},
     {when:'UT 2023-01-03 20',text:'Mars 0.54°N of Moon'},
     {when:'UT 2023-01-04 01',text:'Aldebaran 8.05°S of Moon'},
     {when:'UT 2023-01-04 16',text:'Earth at perihelion 0.983295578 AU'},
     {when:'UT 2023-01-06 23',text:'Full Moon 405790.853 km'},
     {when:'UT 2023-01-07 13',text:'Mercury in inferior conjunction 2.79° North'},
     {when:'UT 2023-01-07 14',text:'Pollux 1.91°N of Moon'},
     {when:'UT 2023-01-08 09',text:'Moon at apogee 406458.329 km'},
     {when:'UT 2023-01-10 12',text:'Regulus 4.61°S of Moon'},
     {when:'UT 2023-01-14 23',text:'Spica 3.84°S of Moon'},
     {when:'UT 2023-01-15 02',text:'Last Quarter 387528.963 km'},
     {when:'UT 2023-01-18 10',text:'Antares 2.10°S of Moon'},
     {when:'UT 2023-01-20 08',text:'Mercury 6.94°N of Moon'},
     {when:'UT 2023-01-21 21',text:'New Moon 356568.957 km'},
     {when:'UT 2023-01-22 20',text:'Saturn 0.37°N of Venus'},
     {when:'UT 2023-01-23 07',text:'Saturn 3.83°N of Moon'},
     {when:'UT 2023-01-23 08',text:'Venus 3.45°N of Moon'},
     {when:'UT 2023-01-25 06',text:'Neptune 2.69°N of Moon'},
     {when:'UT 2023-01-26 02',text:'Jupiter 1.81°N of Moon'},
     {when:'UT 2023-01-28 15',text:'First Quarter 387404.920 km'},
     {when:'UT 2023-01-29 04',text:'Uranus 0.95°S of Moon'},
     {when:'UT 2023-01-30 06',text:'Mercury at greatest elongation 25.0° West'},
     {when:'UT 2023-01-31 04',text:'Mars 0.11°N of Moon'},
     {when:'UT 2023-01-31 07',text:'Aldebaran 8.25°S of Moon'}
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

  /* Coords have equinox 2022.5. */
  var build_messiers = function(){
    var messier = [110];
    messier[0]=["M1",1.4654425,0.3844996,8.4,"Tau","NB","!! famous Crab Neb. supernova remnant","Crab Nebula"];
    messier[1]=["M2",5.6489706,0.0160114,6.5,"Aqr","GC","200-mm telescope needed to resolve",""];
    messier[2]=["M3",3.5920453,0.4934123,6.2,"CVn","GC","!! contains many variable stars",""];
    messier[3]=["M4",4.2977952,-0.4639807,5.6,"Sco","GC","bright globular near Antares",""];
    messier[4]=["M5",4.0131208,0.0349496,5.6,"Ser","GC","!! one of the sky's finest globulars",""];
    messier[5]=["M6",4.6319637,-0.5624695,4.2,"Sco","OC","!! Butterfly Cluster; best at low power","Butterfly Cluster"];
    messier[6]=["M7",4.6923245,-0.6077165,3.3,"Sco","OC","!! excellent in binocs or rich-field scope","Ptolemy's Cluster"];
    messier[7]=["M8",4.734992,-0.4255266,6,"Sgr","NB","!! Lagoon Nebula w/open cl. NGC 6530","Lagoon Nebula"];
    messier[8]=["M9",4.5401186,-0.3235578,7.7,"Oph","GC","smallest of Ophiuchus globulars",""];
    messier[9]=["M10",4.4431193,-0.0721456,6.6,"Oph","GC","rich globular cluster; M12 is 3°NW",""];
    messier[10]=["M11",4.94062,-0.108885,6.3,"Sct","OC","!! Wild Duck Cl.; the best open cluster?","Wild Duck Cluster"];
    messier[11]=["M12",4.3998423,-0.0347114,6.7,"Oph","GC","loose globular cluster near M10",""];
    messier[12]=["M13",4.3742504,0.6357345,5.8,"Her","GC","!! Hercules Cluster; NGC 6207 0.5°NE","Great Hercules Globular"];
    messier[13]=["M14",4.619806,-0.0569309,7.6,"Oph","GC","200-mm telescope needed to resolve",""];
    messier[14]=["M15",5.6334313,0.214086,6.2,"Peg","GC","rich, compact globular","Great Pegasus Globular"];
    messier[15]=["M16",4.7999854,-0.2403793,6.4,"Ser","NB","Eagle Neb. w/open cl.; use neb. filter","Eagle Nebula"];
    messier[16]=["M17",4.8088093,-0.2822482,7,"Sgr","NB","!! Swan or Omega Nebula; use neb. filter","Omega Nebula"];
    messier[17]=["M18",4.8049219,-0.2988373,7.5,"Sgr","OC","sparse cluster; 1°S of M17",""];
    messier[18]=["M19",4.4680127,-0.4589752,6.8,"Oph","GC","oblate globular; M62 4°S",""];
    messier[19]=["M20",4.7296947,-0.4019762,9,"Sgr","NB","!! Trifid Nebula; look for dark lanes","Trifid Nebula"];
    messier[20]=["M21",4.7383973,-0.3926487,6.5,"Sgr","OC","0.7°NE of M20; sparse cluster",""];
    messier[21]=["M22",4.8772013,-0.4167815,5.1,"Sgr","GC","spectacular from southern latitude","Sagittarius Cluster"];
    messier[22]=["M23",4.7042115,-0.3319276,6.9,"Sgr","OC","bright, loose open cluster",""];
    messier[23]=["M24",4.79189,-0.3227186,4.6,"Sgr","OC","rich star cloud; best in big binoculars","Sagittarius Star Cloud"];
    messier[24]=["M25",4.8560573,-0.3356691,6.5,"Sgr","OC","bright but sparse open cluster",""];
    messier[25]=["M26",4.9149972,-0.1636267,8,"Sct","OC","bright, coarse cluster",""];
    messier[26]=["M27",5.2384804,0.3975744,7.4,"Vul","NB","!! Dumbbell Nebula; a superb object","Dumbbell Nebula"];
    messier[27]=["M28",4.825329,-0.4337654,6.8,"Sgr","GC","compact globular near M22",""];
    messier[28]=["M29",5.3438939,0.673821,7.1,"Cyg","OC","small, poor open cluster 2°S of γ Cygni",""];
    messier[29]=["M30",5.6796294,-0.402829,7.2,"Cap","GC","toughest in one-night Messier marathon",""];
    messier[30]=["M31",0.1877724,0.7223881,3.4,"And","GY","!! Andromeda Gal.; look for dust lanes","Andromeda Galaxy"];
    messier[31]=["M32",0.1921389,0.715405,8.1,"And","GY","closest companion to M31",""];
    messier[32]=["M33",0.4152683,0.5369463,5.7,"Tri","GY","large diffuse spiral; requires dark sky","Triangulum Galaxy"];
    messier[33]=["M34",0.7132111,0.7483679,5.5,"Per","OC","best at low power",""];
    messier[34]=["M35",1.6156493,0.4246053,5.3,"Gem","OC","!! look for sm. cluster NGC 2158 0.25°S",""];
    messier[35]=["M36",1.4730195,0.5959596,6.3,"Aur","OC","bright but scattered group; use low pow.",""];
    messier[36]=["M37",1.5440618,0.5681702,6.2,"Aur","OC","!! finest of three Auriga clusters; very rich",""];
    messier[37]=["M38",1.4408212,0.6257001,7.4,"Aur","OC","look for small cluster NGC 1907 0.5°S",""];
    messier[38]=["M39",5.6418368,0.8470706,4.6,"Cyg","OC","very sparse cluster; use low power",""];
    messier[39]=["M40",3.244013,1.0115701,8.4,"UMa","OC","double star Winneke 4; separation 50arcsec","Winnecke 4"];
    messier[40]=["M41",1.7800938,-0.3623146,4.6,"CMa","OC","4°S of Sirius; bright but coarse",""];
    messier[41]=["M42",1.4682831,-0.0948914,4,"Ori","NB","!! Orion Nebula; finest in northern sky","Great Nebula in Orion"];
    messier[42]=["M43",1.4691627,-0.0916936,9,"Ori","NB","detached part of Orion Nebula","De Mairan's Nebula"];
    messier[43]=["M44",2.2750021,0.3473643,3.7,"Cnc","OC","!! Beehive or Praesepe; use low power","Beehive Cluster"];
    messier[44]=["M45",0.9963275,0.4221086,1.6,"Tau","OC","!! Pleiades; look for subtle nebulosity","Pleiades"];
    messier[45]=["M46",2.0194917,-0.2595435,6,"Pup","OC","!! contains planetary nebula NGC 2438",""];
    messier[46]=["M47",1.9968087,-0.2539717,5.2,"Pup","OC","coarse cluster 1.5°W of M46",""];
    messier[47]=["M48",2.1594546,-0.1024386,5.5,"Hya","OC","former lost Messier; large sparse cl.",""];
    messier[48]=["M49",3.2766112,0.1374592,8.4,"Vir","GY","very bright elliptical",""];
    messier[49]=["M50",1.8512815,-0.1460443,6.3,"Mon","OC","between Sirius & Procyon; use low mag",""];
    messier[50]=["M51",3.5384178,0.8214864,8.4,"CVn","GY","!! Whirlpool Galaxy; superb in big scope","Whirlpool Galaxy"];
    messier[51]=["M52",6.1313886,1.0769923,7.3,"Cas","OC","young, rich cl.; faint Bubble Neb. nearby",""];
    messier[52]=["M53",3.4644854,0.3149932,7.6,"Com","GC","150-mm telescope needed to resolve",""];
    messier[53]=["M54",4.959088,-0.5315073,7.6,"Sgr","GC","not easily resolved",""];
    messier[54]=["M55",5.154939,-0.5395402,6.3,"Sgr","GC","bright, loose globular cluster",""];
    messier[55]=["M56",5.0504498,0.5275196,8.3,"Lyr","GC","within a rich starfield",""];
    messier[56]=["M57",4.9499118,0.5770509,8.8,"Lyr","NB","!! Ring Nebula; an amazing smoke ring","Ring Nebula"];
    messier[57]=["M58",3.311046,0.2040838,9.7,"Vir","GY","bright barred spiral; M59 and M60 1°E",""];
    messier[58]=["M59",3.3298011,0.2011822,9.6,"Vir","GY","bright elliptical paired with M60",""];
    messier[59]=["M60",3.3372163,0.1994399,8.8,"Vir","GY","bright elliptical with M59 and NGC 4647",""];
    messier[60]=["M61",3.2421646,0.0757822,9.7,"Vir","GY","face-on two-armed spiral",""];
    messier[61]=["M62",4.4620858,-0.5261832,6.5,"Oph","GC","asymmetrical; in rich field",""];
    messier[62]=["M63",3.4767215,0.7315539,8.6,"CVn","GY","!! Sunflower Galaxy; bright, elongated","Sunflower Galaxy"];
    messier[63]=["M64",3.3938106,0.3763272,8.5,"Com","GY","!! Black Eye Gal; eye needs big scope","Black Eye Galaxy"];
    messier[64]=["M65",2.9673807,0.226195,9.3,"Leo","GY","!! bright elongated spiral",""];
    messier[65]=["M66",2.9730495,0.2244476,8.9,"Leo","GY","!! M65 and NGC 3628 in same field",""];
    messier[66]=["M67",2.3196729,0.2047556,6.1,"Cnc","OC","one of the oldest star clusters known",""];
    messier[67]=["M68",3.319168,-0.4690285,7.8,"Hya","GC","150-mm telescope needed to resolve",""];
    messier[68]=["M69",4.8557997,-0.5643085,7.6,"Sgr","GC","small, poor globular cluster",""];
    messier[69]=["M70",4.9072724,-0.5633248,7.9,"Sgr","GC","small globular 2°E of M69",""];
    messier[70]=["M71",5.2133128,0.3288767,8.2,"Sge","GC","loose globular; looks like an open cluster",""];
    messier[71]=["M72",5.4748084,-0.2172422,9.3,"Aqr","GC","near the Saturn Nebula, NGC 7009",""];
    messier[72]=["M73",5.4988014,-0.21895,9,"Aqr","OC","group of four stars only; an asterism",""];
    messier[73]=["M74",0.4272206,0.2774632,9.4,"Psc","GY","faint, elusive spiral; tough in small scope",""];
    messier[74]=["M75",5.2683829,-0.3813695,8.5,"Sgr","GC","small and distant; 59,000 ly away",""];
    messier[75]=["M76",0.4530366,0.9019768,10.1,"Per","NB","Little Dumbell; faint but distinct","Little Dumbbell Nebula"];
    messier[76]=["M77",0.7149464,0.0022363,8.9,"Cet","GY","a Seyfert galaxy; with starlike nucleus",""];
    messier[77]=["M78",1.5177979,0.000994,8.3,"Ori","NB","bright featureless reflection nebula",""];
    messier[78]=["M79",1.4199436,-0.4281454,7.7,"Lep","GC","200-mm telescope needed to resolve",""];
    messier[79]=["M80",4.2688359,-0.4020788,7.3,"Sco","GC","very compressed globular",""];
    messier[80]=["M81",2.6067514,1.2035642,6.9,"UMa","GY","!! bright spiral visible in binoculars","Bode's Galaxy"];
    messier[81]=["M82",2.6077162,1.2143261,8.4,"UMa","GY","!! the exploding galaxy; M81 0.5°S","Cigar Galaxy"];
    messier[82]=["M83",3.5703867,-0.5232625,7.6,"Hya","GY","large and diffuse; superb from far south","Southern Pinwheel"];
    messier[83]=["M84",3.2560882,0.222684,9.1,"Vir","GY","!! w/M86 in Markarian's Chain",""];
    messier[84]=["M85",3.2578079,0.3154778,9.1,"Com","GY","bright elliptical shape",""];
    messier[85]=["M86",3.2608851,0.2238487,8.9,"Vir","GY","!! w/many NGC galaxies in Chain",""];
    messier[86]=["M87",3.2809495,0.214255,8.6,"Vir","GY","famous jet and black hole",""];
    messier[87]=["M88",3.2866075,0.2497451,9.6,"Com","GY","bright multiple-arm spiral",""];
    messier[88]=["M89",3.3023187,0.2168799,9.8,"Vir","GY","elliptical; resembles M87 but smaller",""];
    messier[89]=["M90",3.307112,0.2276444,9.5,"Vir","GY","bright barred spiral near M89",""];
    messier[90]=["M91",3.3014341,0.2509135,10.2,"Com","GY","some lists say M91=M58, not NGC 4548",""];
    messier[91]=["M92",4.5282219,0.7524151,6.4,"Her","GC","9°NE of M13; fine but often overlooked",""];
    messier[92]=["M93",2.0313631,-0.4175195,6,"Pup","OC","compact, bright cluster; fairly rich",""];
    messier[93]=["M94",3.3682936,0.7157807,8.2,"CVn","GY","very bright and comet-like",""];
    messier[94]=["M95",2.8151574,0.2021345,9.7,"Leo","GY","bright barred spiral",""];
    messier[95]=["M96",2.8273709,0.2041622,9.2,"Leo","GY","M95 in same field",""];
    messier[96]=["M97",2.9500043,0.9580769,9.9,"UMa","NB","Owl Nebula; distinct grey oval","Owl Nebula"];
    messier[97]=["M98",3.2072379,0.258163,10.1,"Com","GY","nearly edge-on spiral near star 6 Com. B.",""];
    messier[98]=["M99",3.2290437,0.2497308,9.9,"Com","GY","nearly face-on spiral near M98",""];
    messier[99]=["M100",3.2469174,0.2741691,9.3,"Com","GY","face-on spiral with starlike nucleus",""];
    messier[100]=["M101",3.6826236,0.9467104,7.9,"UMa","GY","!! Pinwheel Gal; diffuse face-on spiral","Pinwheel Galaxy"];
    messier[101]=["M102",3.9580495,0.9718126,9.9,"Dra","GY","or is M102=M101? (look for NGC 5907)",""];
    messier[102]=["M103",0.4132499,1.06142,7.4,"Cas","OC","three NGC open clusters nearby",""];
    messier[103]=["M104",3.3212369,-0.2049012,8,"Vir","GY","!! Sombrero Galaxy; look for dust lane","Sombrero Galaxy"];
    messier[104]=["M105",2.8317417,0.2175402,9.3,"Leo","GY","bright elliptical near M95 and M96",""];
    messier[105]=["M106",3.2288907,0.8236532,8.4,"CVn","GY","!! superb large, bright spiral",""];
    messier[106]=["M107",4.3361016,-0.2285745,7.9,"Oph","GC","small, faint globular",""];
    messier[107]=["M108",2.9356651,0.9694278,10,"UMa","GY","nearly edge-on; paired with M97 0.75°SE",""];
    messier[108]=["M109",3.1361758,0.9295287,9.8,"UMa","GY","barred spiral near γ UMA",""];
    messier[109]=["M110",0.1816573,0.7296627,8.5,"And","GY","more distant companion to M31",""];    
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
  
  /* Coords have equinox 2022.5. */
  var build_caldwells = function(){
    var caldwell = [42];
    caldwell[0]=["C68",4.9890929,-0.6443089,9.7,"CrA","Bn","","","NGC_6729"];
    caldwell[1]=["C69",4.5170205,-0.6479487,12.8,"Sco","Pl","","Bug Nebula","NGC_6302"];
    caldwell[2]=["C70",0.2441746,-0.6555756,8.1,"Scl","Sp","","","NGC_300"];
    caldwell[3]=["C71",2.0642919,-0.6738567,5.8,"Pup","Oc","","","NGC_2477"];
    caldwell[4]=["C72",0.0699254,-0.6816969,8.2,"Scl","Sb","","","NGC_55"];
    caldwell[5]=["C73",1.3737507,-0.6985729,7.3,"Col","Gc","","","NGC_1851"];
    caldwell[6]=["C74",2.6557482,-0.7076259,8.2,"Vel","Pl","","Eight Burst Nebula","NGC_3132"];
    caldwell[7]=["C75",4.3072483,-0.7106357,5.8,"Sco","Oc","","","NGC_6124"];
    caldwell[8]=["C76",4.4313187,-0.7301613,2.6,"Sco","Oc","","","NGC_6231"];
    caldwell[9]=["C77",3.520439,-0.7528161,7,"Cen","Px","","Centaurus A","Centaurus_A"];
    caldwell[10]=["C78",4.7544148,-0.7626248,6.6,"CrA","Gc","","","NGC_6541"];
    caldwell[11]=["C79",2.6988298,-0.8120972,6.7,"Vel","Gc","","","NGC_3201"];
    caldwell[12]=["C80",3.5262509,-0.8307694,3.6,"Cen","Gc","","Omega Centauri","Omega_Centauri"];
    caldwell[13]=["C81",4.5693242,-0.84535,8.1,"Ara","Gc","","","NGC_6352"];
    caldwell[14]=["C82",4.3763808,-0.8518674,5.2,"Ara","Oc","","","NGC_6193"];
    caldwell[15]=["C83",3.4327141,-0.8654523,9.5,"Cen","Sp","","","NGC_4945"];
    caldwell[16]=["C84",3.6121169,-0.8984693,7.6,"Cen","Gc","","","NGC_5286"];
    caldwell[17]=["C85",2.2726055,-0.9275972,2.5,"Vel","Oc","","Omicron Vel Cluster","IC_2391"];
    caldwell[18]=["C86",4.6361721,-0.9368352,5.6,"Ara","Gc","","","NGC_6397"];
    caldwell[19]=["C87",0.8417577,-0.962254,8.4,"Hor","Gc","","","NGC_1261"];
    caldwell[20]=["C88",3.9592186,-0.9719043,7.9,"Cir","Oc","","","NGC_5823"];
    caldwell[21]=["C89",4.2794496,-1.0114709,5.4,"Nor","Oc","","S Norma Cluster","NGC_6087"];
    caldwell[22]=["C90",2.4523406,-1.0195031,9.7,"Car","Pl","","","NGC_2867"];
    caldwell[23]=["C91",2.9119234,-1.0260543,3,"Car","Oc","","","NGC_3532"];
    caldwell[24]=["C92",2.8129137,-1.0469383,6.2,"Car","Bn","","Eta Carinae Nebula","Carina_Nebula"];
    caldwell[25]=["C93",5.0303771,-1.0462321,5.4,"Pav","Gc","","","NGC_6752"];
    caldwell[26]=["C94",3.3814013,-1.0551405,4.2,"Cru","Oc","","","Jewel Box"];
    caldwell[27]=["C95",4.2133555,-1.0569786,5.1,"TrA","Oc","","","NGC_6025"];
    caldwell[28]=["C96",2.088595,-1.0634043,3.8,"Car","Oc","","","NGC_2516"];
    caldwell[29]=["C97",3.041928,-1.0775886,5.3,"Cen","Oc","","","NGC_3766"];
    caldwell[30]=["C98",3.3319937,-1.1011236,6.9,"Cru","Oc","","","NGC_4609"];
    caldwell[31]=["C99",3.3788792,-1.1016839,undefined,"Cru","Dn","","Coalsack Nebula","Coalsack_Nebula"];
    caldwell[32]=["C100",3.0440934,-1.1023146,4.5,"Cen","Oc","","Lambda Centauri Nebula","IC_2944"];
    caldwell[33]=["C101",5.0262188,-1.1137274,9,"Pav","Sb","","","NGC_6744"];
    caldwell[34]=["C102",2.8100241,-1.1260579,1.9,"Car","Oc","","Theta Car Cluster","IC_2602"];
    caldwell[35]=["C103",1.4771906,-1.2058189,1,"Dor","Bn","","Tarantula Nebula","Tarantula_Nebula"];
    caldwell[36]=["C104",0.2790753,-1.2344631,6.6,"Tuc","Gc","","","NGC_362"];
    caldwell[37]=["C105",3.4083263,-1.2392585,7.3,"Mus","Gc","","","NGC_4833"];
    caldwell[38]=["C106",0.1094663,-1.2559178,4,"Tuc","Gc","","47 Tucanae","47_Tucanae"];
    caldwell[39]=["C107",4.3126624,-1.2609899,9.3,"Aps","Gc","","","NGC_6101"];
    caldwell[40]=["C108",3.2600085,-1.2704443,7.8,"Mus","Gc","","","NGC_4372"];
    caldwell[41]=["C109",2.6581245,-1.413326,11.6,"Cha","Pl","","","NGC_3195"];    
   return caldwell;
  };
  var caldwells = build_caldwells();

  /* Find a Caldwell object either by 'C103' or 'Tarantula Nebula' (case-insensitive). */  
  var find_caldwell = function(name_raw){
    return find_deep_sky_object(name_raw, caldwells);
  };
  
  /* This exist in order to avoid creating a large number of identical objects. */
  var when_fixed_equinox = when("J2022.5");
  
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

  /* Coords have equinox 2022.5. */
  var build_stars = function(){
    /* Source: Yale Bright Star Catalog r5.  Name, Right Ascension, Declination (J2022.5), and Magnitude.*/
    var ybs = [9096];
    //removal of leading whitespace cuts down on the overall size of this js file
ybs[0]=['',0.0276237,0.7915835,6.7];
ybs[1]=['',0.0271243,-0.0065944,6.29];
ybs[2]=['33 Psc',0.0283045,-0.0974292,4.61];
ybs[3]=['86 Peg',0.029917,0.2359916,5.51];
ybs[4]=['',0.0324819,1.0220975,5.96];
ybs[5]=['',0.0325179,-0.854335,5.7];
ybs[6]=['10 Cas',0.033278,1.1226187,5.59];
ybs[7]=['',0.033926,0.508704,6.13];
ybs[8]=['',0.034825,-0.4011169,6.18];
ybs[9]=['',0.0368751,-0.3012647,6.19];
ybs[10]=['',0.0387786,-0.0423017,6.43];
ybs[11]=['',0.0389455,-0.3906694,5.94];
ybs[12]=['',0.0401384,-0.5830145,5.68];
ybs[13]=['',0.0408147,-0.0405371,6.07];
ybs[14]=['α And',0.0416805,0.5099106,2.06];
ybs[15]=['',0.0411906,-0.1518213,5.99];
ybs[16]=['',0.0429858,0.6414404,6.19];
ybs[17]=['',0.0423398,-0.3046007,6.06];
ybs[18]=['',0.0437775,0.4465937,6.23];
ybs[19]=['',0.046298,1.3934686,6.01];
ybs[20]=['β Cas',0.0452361,1.0345417,2.27];
ybs[21]=['87 Peg',0.0445064,0.3200427,5.53];
ybs[22]=['',0.0443504,-0.9403274,6.33];
ybs[23]=['κ1 Scl',0.0457787,-0.4862947,5.42];
ybs[24]=['ε Phe',0.0460004,-0.7962603,3.88];
ybs[25]=['34 Psc',0.0488523,0.1967105,5.51];
ybs[26]=['22 And',0.0501768,0.8062958,5.03];
ybs[27]=['',0.0509891,0.9999109,6.74];
ybs[28]=['',0.0500227,-0.0894217,5.84];
ybs[29]=['γ3 Oct',0.0480884,-1.4328936,5.28];
ybs[30]=['',0.0517535,-0.2173788,5.85];
ybs[31]=['',0.0511198,-1.2758239,6.64];
ybs[32]=['6 Cet',0.0541534,-0.2677852,4.89];
ybs[33]=['κ2 Scl',0.055469,-0.4830135,5.41];
ybs[34]=['θ Scl',0.0561457,-0.6110044,5.25];
ybs[35]=['',0.0574603,0.8426027,6.16];
ybs[36]=['',0.0580797,-0.3109001,5.25];
ybs[37]=['',0.061156,0.6600553,6.73];
ybs[38]=['γ Peg',0.0628235,0.2671863,2.83];
ybs[39]=['',0.0635683,0.4731981,6.3];
ybs[40]=['23 And',0.064112,0.7183828,5.72];
ybs[41]=['',0.0647503,-0.4519866,5.94];
ybs[42]=['',0.0649021,-0.4565729,6.31];
ybs[43]=['',0.066377,0.5817378,6.25];
ybs[44]=['χ Peg',0.0688043,0.3548544,4.8];
ybs[45]=['',0.068106,-0.1336147,5.12];
ybs[46]=['',0.0616898,-1.4812458,5.77];
ybs[47]=['7 Cet',0.0688613,-0.3282579,4.44];
ybs[48]=['',0.0702587,0.3911133,6.24];
ybs[49]=['35 Psc',0.0704175,0.1561338,5.79];
ybs[50]=['',0.0700569,-0.1648419,5.75];
ybs[51]=['',0.0710827,0.5525852,6.45];
ybs[52]=['',0.0713301,0.4783602,6.35];
ybs[53]=['',0.0702476,-0.6070163,6.17];
ybs[54]=['',0.076561,1.3452257,6.35];
ybs[55]=['',0.0765699,0.7630517,6.15];
ybs[56]=['',0.075395,-0.5466626,5.67];
ybs[57]=['',0.0738749,-1.3227231,6.49];
ybs[58]=['36 Psc',0.0773486,0.1459952,6.11];
ybs[59]=['',0.0793068,1.0761391,5.74];
ybs[60]=['',0.0778752,-0.3505607,6.47];
ybs[61]=['',0.080058,0.8390214,5.89];
ybs[62]=['θ And',0.0797435,0.6773022,4.61];
ybs[63]=['',0.0775819,-1.3728001,6.77];
ybs[64]=['',0.0825546,0.8998555,6.14];
ybs[65]=['',0.0815195,-0.3303252,6.45];
ybs[66]=['',0.0826824,0.0316559,6.17];
ybs[67]=['σ And',0.0851392,0.644203,4.52];
ybs[68]=['',0.0848582,0.1977576,6.05];
ybs[69]=['26 And',0.0868098,0.7664776,6.11];
ybs[70]=['',0.0864695,0.5522579,5.87];
ybs[71]=['',0.0865856,-0.138369,6.46];
ybs[72]=['',0.0864973,-0.7524195,6.33];
ybs[73]=['ι Cet',0.0897745,-0.151828,3.56];
ybs[74]=['',0.0911269,0.7130455,6.33];
ybs[75]=['',0.0929025,0.8550374,6.52];
ybs[76]=['ζ Tuc',0.0921936,-1.1301001,4.23];
ybs[77]=['',0.0941928,0.5421092,5.9];
ybs[78]=['',0.095739,0.5765889,5.79];
ybs[79]=['41 Psc',0.0949382,0.1451242,5.37];
ybs[80]=['',0.0963015,0.1937605,6.56];
ybs[81]=['ρ And',0.0973542,0.6648537,5.18];
ybs[82]=['π Tuc',0.0945933,-1.2130086,5.51];
ybs[83]=['ι Scl',0.0988142,-0.5036494,5.18];
ybs[84]=['',0.0999511,-0.3478984,5.12];
ybs[85]=['42 Psc',0.102932,0.2374892,6.23];
ybs[86]=['',0.0978089,-1.3491789,5.97];
ybs[87]=['9 Cet',0.1047438,-0.2109202,6.39];
ybs[88]=['',0.1061686,-0.5395078,6.55];
ybs[89]=['',0.1100844,0.6754731,7.39];
ybs[90]=['',0.1111899,0.9100937,5.57];
ybs[91]=['12 Cas',0.1136589,1.0813293,5.4];
ybs[92]=['',0.1119023,-0.0365586,6.07];
ybs[93]=['',0.1149064,0.9280163,5.74];
ybs[94]=['44 Psc',0.1158834,0.0360267,5.77];
ybs[95]=['β Hyi',0.1162975,-1.3461676,2.8];
ybs[96]=['α Phe',0.1194825,-0.7362096,2.39];
ybs[97]=['κ Phe',0.1191227,-0.7601884,3.94];
ybs[98]=['10 Cet',0.121198,0.0013031,6.19];
ybs[99]=['',0.1237843,-0.4437129,5.98];
ybs[100]=['47 Psc',0.1275039,0.3144619,5.06];
ybs[101]=['',0.1284702,0.7769982,5.17];
ybs[102]=['η Scl',0.1267164,-0.5739153,4.81];
ybs[103]=['48 Psc',0.1282095,0.2891884,6.06];
ybs[104]=['',0.1287163,0.1800131,6.04];
ybs[105]=['',0.1286376,-0.3527438,6.43];
ybs[106]=['',0.1288944,-0.6944794,5.43];
ybs[107]=['',0.1315323,0.6461946,6.26];
ybs[108]=['',0.1300018,-0.8797949,6.26];
ybs[109]=['',0.141243,1.3464084,6.21];
ybs[110]=['',0.1378895,1.0489663,5.94];
ybs[111]=['28 And',0.1366297,0.5214311,5.23];
ybs[112]=['',0.1352658,-0.2572616,6.14];
ybs[113]=['',0.1349433,-0.5583745,6.57];
ybs[114]=['12 Cet',0.1360861,-0.0668998,5.72];
ybs[115]=['',0.1374529,-0.4130086,5.19];
ybs[116]=['',0.1376979,-0.7123618,6.19];
ybs[117]=['',0.1375005,-0.8393442,5.69];
ybs[118]=['13 Cas',0.1428389,1.1631482,6.18];
ybs[119]=['',0.1423591,0.5882756,5.87];
ybs[120]=['λ Cas',0.1441022,0.9537567,4.73];
ybs[121]=['',0.1436971,0.9243869,5.6];
ybs[122]=['λ1 Phe',0.1417665,-0.8496187,4.77];
ybs[123]=['β1 Tuc',0.1420773,-1.0966605,4.37];
ybs[124]=['β2 Tuc',0.1421423,-1.0967962,4.54];
ybs[125]=['',0.1469058,0.7612897,6.7];
ybs[126]=['',0.1513659,1.241026,6.42];
ybs[127]=['κ Cas',0.149649,1.1005275,4.16];
ybs[128]=['52 Psc',0.1473568,0.3563683,5.38];
ybs[129]=['51 Psc',0.1464274,0.123561,5.67];
ybs[130]=['',0.1473322,0.4835349,6.67];
ybs[131]=['',0.1484001,0.4957471,6.3];
ybs[132]=['',0.1502358,0.960261,5.93];
ybs[133]=['β3 Tuc',0.1472241,-1.0979371,5.09];
ybs[134]=['16 Cas',0.1559737,1.1671728,6.48];
ybs[135]=['',0.1518261,-0.5137283,5.55];
ybs[136]=['θ Tuc',0.1497703,-1.2416659,6.13];
ybs[137]=['',0.1549761,-0.9119215,5.57];
ybs[138]=['',0.1574868,0.23553,6.4];
ybs[139]=['13 Cet',0.1588106,-0.0605462,5.2];
ybs[140]=['14 Cet',0.1601307,-0.0066644,5.93];
ybs[141]=['',0.1631997,0.9475788,5.08];
ybs[142]=['',0.1618077,0.2326584,6.41];
ybs[143]=['',0.1647166,1.055047,5.79];
ybs[144]=['λ2 Phe',0.1603556,-0.8356136,5.51];
ybs[145]=['',0.1596918,-0.9546642,6.06];
ybs[146]=['',0.163747,0.4778425,6.5];
ybs[147]=['',0.1622369,-0.2591804,6.45];
ybs[148]=['',0.1624668,-0.3965185,6.06];
ybs[149]=['',0.1658498,0.7786299,5.13];
ybs[150]=['ζ Cas',0.1668418,0.942836,3.66];
ybs[151]=['π And',0.1661965,0.5906724,4.36];
ybs[152]=['53 Psc',0.1656482,0.2679999,5.89];
ybs[153]=['',0.1671584,0.421283,6.47];
ybs[154]=['',0.1682659,0.6199932,5.48];
ybs[155]=['',0.1815971,1.4419419,6.4];
ybs[156]=['',0.1678144,-0.4301132,5.57];
ybs[157]=['',0.164073,-1.1344831,6.42];
ybs[158]=['',0.1687129,0.056877,6.39];
ybs[159]=['',0.1672905,-0.9472007,6.41];
ybs[160]=['ε And',0.1734688,0.5137395,4.37];
ybs[161]=['',0.1763622,0.863551,5.43];
ybs[162]=['δ And',0.176861,0.5407764,3.27];
ybs[163]=['54 Psc',0.1769346,0.3730453,5.87];
ybs[164]=['55 Psc',0.1793961,0.3763217,5.36];
ybs[165]=['α Cas',0.1823748,0.9889118,2.23];
ybs[166]=['',0.1725847,-1.2743308,6.85];
ybs[167]=['',0.1791536,-0.5905907,6.69];
ybs[168]=['',0.178603,-0.7796969,6.01];
ybs[169]=['',0.1815295,-0.2861238,6.49];
ybs[170]=['',0.1817857,-0.4133147,6.14];
ybs[171]=['',0.1826186,-0.0738049,5.91];
ybs[172]=['32 And',0.1847781,0.6908328,5.33];
ybs[173]=['',0.1808278,-1.0355243,5.89];
ybs[174]=['',0.189458,1.1566401,5.83];
ybs[175]=['',0.18673,0.4320094,6.04];
ybs[176]=['ξ Cas',0.1890675,0.8837578,4.8];
ybs[177]=['μ Phe',0.1849418,-0.8021851,4.59];
ybs[178]=['',0.1912293,1.0275867,6.17];
ybs[179]=['ξ Phe',0.1866978,-0.9839909,5.7];
ybs[180]=['π Cas',0.1951478,0.8228821,4.94];
ybs[181]=['λ1 Scl',0.1910859,-0.6691644,6.06];
ybs[182]=['',0.1906159,-1.0496315,5.98];
ybs[183]=['ρ Tuc',0.1894667,-1.1404851,5.39];
ybs[184]=['β Cet',0.1950936,-0.3117808,2.04];
ybs[185]=['',0.199411,0.8375314,5.67];
ybs[186]=['',0.1962161,-0.2074977,6.02];
ybs[187]=['η Phe',0.1935461,-1.0007731,4.36];
ybs[188]=['21 Cas',0.2058643,1.3109301,5.66];
ybs[189]=['ο Cas',0.2006646,0.8448661,4.54];
ybs[190]=['φ1 Cet',0.1977679,-0.183025,4.76];
ybs[191]=['λ2 Scl',0.1975623,-0.6684399,5.9];
ybs[192]=['',0.2032603,0.9659424,5.42];
ybs[193]=['',0.2000739,-0.3819354,5.24];
ybs[194]=['',0.2007738,-0.742705,5.94];
ybs[195]=['',0.1985493,-1.0886478,6.07];
ybs[196]=['',0.2097734,1.2120892,6.33];
ybs[197]=['',0.203099,-0.0786518,6.15];
ybs[198]=['',0.2007951,-0.9353603,6.15];
ybs[199]=['18 Cet',0.2033692,-0.2226707,6.15];
ybs[200]=['',0.207486,0.9673998,6.52];
ybs[201]=['',0.2069795,0.7851197,6.05];
ybs[202]=['',0.2042852,-0.2845139,6.47];
ybs[203]=['',0.2095939,1.0419099,6.39];
ybs[204]=['23 Cas',0.2151648,1.3084727,5.41];
ybs[205]=['',0.2042191,-0.8277913,5.8];
ybs[206]=['',0.2064196,-0.3909412,5.5];
ybs[207]=['57 Psc',0.2082686,0.2722396,5.38];
ybs[208]=['',0.2166227,1.2705547,5.87];
ybs[209]=['58 Psc',0.2103132,0.211123,5.5];
ybs[210]=['59 Psc',0.2112586,0.3438549,6.13];
ybs[211]=['ζ And',0.2117894,0.4256815,4.06];
ybs[212]=['60 Psc',0.2118783,0.1197882,5.99];
ybs[213]=['61 Psc',0.2142693,0.3673525,6.54];
ybs[214]=['',0.2131083,-0.3130928,5.7];
ybs[215]=['η Cas',0.2200211,1.0112116,3.44];
ybs[216]=['',0.2143686,-0.3769918,5.57];
ybs[217]=['62 Psc',0.2157964,0.1295457,5.93];
ybs[218]=['',0.2161872,0.0942996,5.75];
ybs[219]=['ν Cas',0.2186936,0.8917007,4.89];
ybs[220]=['δ Psc',0.2175154,0.1345191,4.43];
ybs[221]=['64 Psc',0.2188835,0.2978038,5.07];
ybs[222]=['ν And',0.2228004,0.7190954,4.53];
ybs[223]=['',0.2205826,-0.2345564,5.59];
ybs[224]=['',0.2196349,-0.4191246,5.9];
ybs[225]=['',0.2180989,-0.8128945,6.27];
ybs[226]=['65 Psc',0.2229259,0.4857787,7];
ybs[227]=['65 Psc',0.2229549,0.485769,7.1];
ybs[228]=['',0.2210592,-0.4056038,6.28];
ybs[229]=['',0.2273796,1.1234619,5.39];
ybs[230]=['',0.2250117,0.7875694,6.15];
ybs[231]=['φ2 Cet',0.2236609,-0.1836476,5.19];
ybs[232]=['λ Hyi',0.215332,-1.3055224,5.07];
ybs[233]=['',0.2296728,1.0808456,6.07];
ybs[234]=['',0.2279873,0.9011163,6.39];
ybs[235]=['',0.2230153,-0.7552476,6.48];
ybs[236]=['',0.2493576,1.4630879,5.62];
ybs[237]=['',0.2306356,0.9022155,6.21];
ybs[238]=['ρ Phe',0.2255975,-0.8877582,5.22];
ybs[239]=['',0.2289216,0.0612099,6.37];
ybs[240]=['',0.2375096,1.0689444,4.82];
ybs[241]=['',0.2308774,-0.7607395,6.9];
ybs[242]=['',0.2362119,0.6749271,6.69];
ybs[243]=['',0.2346536,-0.4168532,5.46];
ybs[244]=['20 Cet',0.2363147,-0.0178427,4.77];
ybs[245]=['',0.2387306,0.6551939,6.06];
ybs[246]=['',0.2404261,0.9217243,6.27];
ybs[247]=['',0.2369558,-0.4303129,6.46];
ybs[248]=['λ1 Tuc',0.2323597,-1.210953,6.22];
ybs[249]=['υ1 Cas',0.2458982,1.0313912,4.83];
ybs[250]=['66 Psc',0.2433931,0.3370228,5.74];
ybs[251]=['21 Cet',0.241852,-0.1504324,6.16];
ybs[252]=['',0.2459927,0.8517239,6.27];
ybs[253]=['',0.2380488,-1.0951871,5.7];
ybs[254]=['36 And',0.2450993,0.4145145,5.47];
ybs[255]=['',0.2463252,0.4307212,6.2];
ybs[256]=['',0.2511775,1.014352,6.21];
ybs[257]=['',0.2548271,1.202487,6.37];
ybs[258]=['67 Psc',0.2495442,0.4770144,6.09];
ybs[259]=['',0.24803,-0.1261125,5.85];
ybs[260]=['γ Cas',0.2534371,1.0618238,2.47];
ybs[261]=['υ2 Cas',0.2531817,1.0350234,4.63];
ybs[262]=['',0.2537589,1.0556471,5.55];
ybs[263]=['φ3 Cet',0.2493811,-0.1945205,5.31];
ybs[264]=['',0.2487704,-0.4826547,6.1];
ybs[265]=['μ And',0.2530966,0.6740602,3.87];
ybs[266]=['λ2 Tuc',0.2436386,-1.2113515,5.45];
ybs[267]=['η And',0.2548799,0.4108296,4.42];
ybs[268]=['',0.257196,0.80217,6.12];
ybs[269]=['',0.2616461,1.1601783,5.97];
ybs[270]=['68 Psc',0.2576982,0.5081252,5.42];
ybs[271]=['',0.2595118,0.5946684,5.98];
ybs[272]=['',0.2578397,0.2411528,6.32];
ybs[273]=['',0.2596975,0.3756924,6.37];
ybs[274]=['',0.2707674,1.2409966,6.39];
ybs[275]=['φ4 Cet',0.261185,-0.196505,5.61];
ybs[276]=['α Scl',0.2604384,-0.5102711,4.31];
ybs[277]=['',0.2587468,-1.0572372,6.23];
ybs[278]=['',0.2676461,0.7824661,6.84];
ybs[279]=['',0.2676607,0.7825049,6.04];
ybs[280]=['',0.2661475,0.1152614,6.11];
ybs[281]=['',0.3152713,1.5043355,4.25];
ybs[282]=['',0.4720912,1.555363,6.46];
ybs[283]=['',0.2776308,0.8928331,6.47];
ybs[284]=['ξ Scl',0.2720558,-0.6771168,5.59];
ybs[285]=['',0.2806831,0.8289715,6.45];
ybs[286]=['39 And',0.2800342,0.7237092,5.98];
ybs[287]=['σ Psc',0.2794998,0.5571953,5.5];
ybs[288]=['',0.2837104,1.0680605,5.92];
ybs[289]=['σ Scl',0.2771142,-0.5485811,5.5];
ybs[290]=['ε Psc',0.2797574,0.1398093,4.28];
ybs[291]=['ω Phe',0.2747838,-0.9927758,6.11];
ybs[292]=['25 Cet',0.2800601,-0.0823132,5.43];
ybs[293]=['',0.2868437,1.0768774,5.84];
ybs[294]=['',0.2852545,0.9184363,5.99];
ybs[295]=['',0.2785108,-0.8076858,5.36];
ybs[296]=['',0.2808614,-0.513221,6.29];
ybs[297]=['26 Cet',0.2834995,0.0259533,6.04];
ybs[298]=['',0.2884502,0.8923902,6.54];
ybs[299]=['',0.2866412,0.5197391,6.19];
ybs[300]=['',0.2773787,-1.1403209,6.21];
ybs[301]=['',0.2874477,0.7000749,6.72];
ybs[302]=['',0.3526673,1.5211412,6.25];
ybs[303]=['73 Psc',0.2881713,0.1008203,6];
ybs[304]=['72 Psc',0.2892061,0.2629559,5.68];
ybs[305]=['',0.295898,1.0974909,6.54];
ybs[306]=['ψ1 Psc',0.2918681,0.3768758,5.34];
ybs[307]=['ψ1 Psc',0.2919262,0.3767352,5.56];
ybs[308]=['',0.3107624,1.3985515,6.29];
ybs[309]=['77 Psc',0.2922799,0.0877617,6.35];
ybs[310]=['77 Psc',0.29244,0.087781,7.25];
ybs[311]=['27 Cet',0.2912151,-0.1720736,6.12];
ybs[312]=['',0.2983687,0.9957948,6.43];
ybs[313]=['28 Cet',0.2932739,-0.1696362,5.58];
ybs[314]=['',0.2989288,0.9358133,6.38];
ybs[315]=['75 Psc',0.2956005,0.2282199,6.12];
ybs[316]=['',0.2932927,-0.4166537,6.14];
ybs[317]=['μ Cas',0.303855,0.9606277,5.17];
ybs[318]=['β Phe',0.292711,-0.813299,3.31];
ybs[319]=['',0.2944873,-0.6203054,6.61];
ybs[320]=['41 And',0.3024263,0.7690205,5.03];
ybs[321]=['',0.2980443,-0.4167246,6.37];
ybs[322]=['',0.3052137,1.018979,5.79];
ybs[323]=['78 Psc',0.3022365,0.5608076,6.25];
ybs[324]=['ψ2 Psc',0.3017789,0.3640559,5.55];
ybs[325]=['30 Cet',0.3006241,-0.1687004,5.82];
ybs[326]=['80 Psc',0.3034169,0.1006943,5.52];
ybs[327]=['υ Phe',0.3002912,-0.7219939,5.21];
ybs[328]=['ι Tuc',0.2975506,-1.0760906,5.37];
ybs[329]=['',0.3241998,1.3926471,5.64];
ybs[330]=['η Cet',0.304196,-0.1756258,3.45];
ybs[331]=['φ And',0.3090113,0.8266121,4.25];
ybs[332]=['31 Cas',0.3150542,1.2024941,5.29];
ybs[333]=['β And',0.309768,0.62378,2.06];
ybs[334]=['ζ Phe',0.3024872,-0.9621334,3.92];
ybs[335]=['ψ3 Psc',0.3099159,0.3451914,5.55];
ybs[336]=['44 And',0.3124343,0.7365411,5.65];
ybs[337]=['',0.312193,0.4464044,5.8];
ybs[338]=['',0.3180938,1.1226287,5.55];
ybs[339]=['θ Cas',0.3162482,0.9646241,4.33];
ybs[340]=['',0.3114879,0.2756486,6.06];
ybs[341]=['32 Cas',0.3193008,1.1368717,5.57];
ybs[342]=['32 Cet',0.3112336,-0.1533581,6.4];
ybs[343]=['33 Cet',0.312937,0.0447648,5.95];
ybs[344]=['45 And',0.3160991,0.6604907,5.81];
ybs[345]=['82 Psc',0.3157345,0.5505449,5.16];
ybs[346]=['',0.3099569,-1.0048746,6.41];
ybs[347]=['χ Psc',0.3170666,0.3692043,4.66];
ybs[348]=['τ Psc',0.3181013,0.5272431,4.51];
ybs[349]=['34 Cet',0.3179652,-0.0372109,5.94];
ybs[350]=['',0.3255642,1.0790436,6.41];
ybs[351]=['',0.3223598,0.7933691,6.11];
ybs[352]=['',0.3239184,0.5267931,6.19];
ybs[353]=['',0.3429802,1.3967549,6.26];
ybs[354]=['',0.320486,-0.5355237,6.52];
ybs[355]=['',0.3219596,-0.6586431,5.92];
ybs[356]=['φ Psc',0.3271393,0.4311372,4.65];
ybs[357]=['ζ Psc',0.3268402,0.1342858,5.24];
ybs[358]=['ζ Psc',0.3269421,0.1343391,6.3];
ybs[359]=['',0.3286628,0.5000087,6.43];
ybs[360]=['87 Psc',0.3286737,0.2836557,5.98];
ybs[361]=['',0.3397177,1.254231,7.83];
ybs[362]=['37 Cet',0.3295659,-0.1362131,5.13];
ybs[363]=['88 Psc',0.3310883,0.1241599,6.03];
ybs[364]=['38 Cet',0.3314843,-0.0149286,5.7];
ybs[365]=['',0.3392312,0.8412569,6.61];
ybs[366]=['ν Phe',0.3323675,-0.7926045,4.96];
ybs[367]=['',0.3384823,0.5800251,6.02];
ybs[368]=['',0.342105,0.7857485,6.34];
ybs[369]=['39 Cet',0.3392533,-0.0415746,5.41];
ybs[370]=['',0.343213,0.5561107,6.73];
ybs[371]=['',0.3589779,1.3559117,6.31];
ybs[372]=['',0.3469168,0.8296885,6.25];
ybs[373]=['κ Tuc',0.3337936,-1.2000482,4.86];
ybs[374]=['89 Psc',0.3445453,0.0651436,5.16];
ybs[375]=['',0.3493575,0.654567,6.46];
ybs[376]=['',0.3396152,-1.1568021,6.24];
ybs[377]=['',0.3659328,1.3326643,6.38];
ybs[378]=['φ Cas',0.3556758,1.0183861,4.98];
ybs[379]=['υ Psc',0.3521572,0.4779036,4.76];
ybs[380]=['35 Cas',0.3604571,1.1305491,6.34];
ybs[381]=['42 Cet',0.3532408,-0.0068286,5.87];
ybs[382]=['',0.3745786,1.3760633,6.07];
ybs[383]=['',0.3565643,-0.0546193,6.23];
ybs[384]=['',0.3559696,-0.1941045,6.15];
ybs[385]=['91 Psc',0.3594187,0.5036223,5.23];
ybs[386]=['ξ And',0.3650984,0.7966734,4.88];
ybs[387]=['',0.3700232,1.0168286,6.45];
ybs[388]=['',0.3655391,0.032175,6.2];
ybs[389]=['43 Cet',0.3653496,-0.0058052,6.49];
ybs[390]=['',0.364775,-0.3309888,6.35];
ybs[391]=['47 And',0.3707494,0.6602908,5.58];
ybs[392]=['',0.3704508,0.5997426,6.29];
ybs[393]=['',0.3692918,0.3592905,5.97];
ybs[394]=['',0.3816439,1.2408667,6.49];
ybs[395]=['ψ Cas',0.3820052,1.1911243,4.74];
ybs[396]=['',0.3689731,-0.5380609,5.84];
ybs[397]=['44 Cet',0.3716225,-0.1377182,6.21];
ybs[398]=['θ Cet',0.3715401,-0.1407871,3.6];
ybs[399]=['δ Cas',0.3808906,1.0533361,2.68];
ybs[400]=['',0.3729538,-0.1186468,5.91];
ybs[401]=['',0.3742234,-0.2712865,6.14];
ybs[402]=['',0.3750536,-0.0476813,6.15];
ybs[403]=['',0.3788607,0.4123894,6.18];
ybs[404]=['',0.3738171,-0.7221437,5.42];
ybs[405]=['',0.3823991,0.760512,5.96];
ybs[406]=['',0.3814686,0.6055614,6.31];
ybs[407]=['',0.3738185,-0.775129,6.26];
ybs[408]=['46 Cet',0.3784112,-0.2527652,4.9];
ybs[409]=['ρ Psc',0.3816723,0.3366495,5.38];
ybs[410]=['94 Psc',0.3835947,0.3378357,5.5];
ybs[411]=['',0.3856517,0.6020286,6.27];
ybs[412]=['',0.382258,-0.0049266,6.41];
ybs[413]=['ω And',0.3883415,0.7945217,4.83];
ybs[414]=['',0.3872919,0.7193667,6.46];
ybs[415]=['',0.3842193,0.0637311,6.58];
ybs[416]=['',0.3746436,-1.1214227,5.93];
ybs[417]=['47 Cet',0.3838421,-0.2258527,5.66];
ybs[418]=['',0.3887589,0.7060138,6.6];
ybs[419]=['',0.3839841,-0.5659546,5.79];
ybs[420]=['α UMi',0.7881203,1.559464,2.02];
ybs[421]=['',0.3878723,-0.1882441,6.13];
ybs[422]=['',0.3907861,0.140976,6.2];
ybs[423]=['38 Cas',0.4054884,1.228363,5.81];
ybs[424]=['',0.4034603,1.1556424,6.14];
ybs[425]=['γ Phe',0.3898189,-0.7540235,3.41];
ybs[426]=['49 And',0.3990884,0.8224478,5.27];
ybs[427]=['',0.3915994,-0.5872636,6.58];
ybs[428]=['97 Psc',0.3974945,0.3223829,6.02];
ybs[429]=['48 Cet',0.3956615,-0.3754857,5.12];
ybs[430]=['μ Psc',0.3986297,0.1092481,4.84];
ybs[431]=['',0.3946906,-0.814033,6.31];
ybs[432]=['',0.3989814,-0.4553955,5.93];
ybs[433]=['η Psc',0.4044386,0.2698474,3.62];
ybs[434]=['',0.4076098,0.6093842,6.39];
ybs[435]=['',0.4141156,1.0200112,5.7];
ybs[436]=['δ Phe',0.4022119,-0.854468,3.95];
ybs[437]=['',0.4047441,-0.5265323,5.82];
ybs[438]=['χ Cas',0.4163635,1.0357948,4.71];
ybs[439]=['',0.4040666,-0.7934315,6.17];
ybs[440]=['',0.4109384,-0.1553302,6.59];
ybs[441]=['',0.4098813,-0.6414134,5.51];
ybs[442]=['',0.4170615,0.6519135,5.88];
ybs[443]=['',0.40808,-0.865905,6.28];
ybs[444]=['',0.4138337,-0.1206103,5.76];
ybs[445]=['',0.4331583,1.2987823,6.58];
ybs[446]=['',0.4190509,0.324197,5.89];
ybs[447]=['49 Cet',0.4176864,-0.2715994,5.63];
ybs[448]=['',0.4241461,0.7189133,6.38];
ybs[449]=['',0.4183225,-0.5546246,6.12];
ybs[450]=['',0.4269153,0.8523657,5.92];
ybs[451]=['101 Psc',0.4231562,0.2578853,6.22];
ybs[452]=['40 Cas',0.437909,1.2767722,5.28];
ybs[453]=['',0.4238143,0.3062691,5.8];
ybs[454]=['υ And',0.4281838,0.7246548,4.09];
ybs[455]=['50 Cet',0.4235928,-0.2667904,5.42];
ybs[456]=['',0.4192299,-1.0127262,6.01];
ybs[457]=['',0.4346551,1.0138842,5.56];
ybs[458]=['τ Scl',0.424008,-0.5199896,5.69];
ybs[459]=['π Psc',0.4288986,0.2139027,5.57];
ybs[460]=['51 And',0.4336466,0.8507112,3.57];
ybs[461]=['',0.4358746,0.7943641,6.36];
ybs[462]=['',0.4308662,-0.1621401,6.24];
ybs[463]=['',0.4094178,-1.36816,6.11];
ybs[464]=['',0.4257175,-1.0150252,6.18];
ybs[465]=['χ And',0.4394343,0.7766651,4.98];
ybs[466]=['',0.4435961,0.9421574,6.39];
ybs[467]=['',0.4339537,-0.6355541,5.94];
ybs[468]=['α Eri',0.4299864,-0.9969794,0.46];
ybs[469]=['',0.4360477,-0.3693398,5.58];
ybs[470]=['',0.4358459,-0.4347313,6.7];
ybs[471]=['105 Psc',0.4402415,0.2883161,5.97];
ybs[472]=['',0.4451417,0.7576648,5.61];
ybs[473]=['τ And',0.4446968,0.7101777,4.94];
ybs[474]=['43 Cas',0.4539482,1.1895438,5.59];
ybs[475]=['',0.4349082,-0.9307001,6.84];
ybs[476]=['42 Cas',0.456885,1.2345609,5.18];
ybs[477]=['',0.4520263,1.0672898,6.71];
ybs[478]=['',0.4529471,1.0252167,6.37];
ybs[479]=['',0.4500311,0.7457192,4.95];
ybs[480]=['',0.4475202,0.4513231,6.17];
ybs[481]=['',0.4491252,0.5263951,5.99];
ybs[482]=['',0.4390679,-0.9788605,5.87];
ybs[483]=['',0.4390971,-0.9788023,5.76];
ybs[484]=['',0.4560972,1.0739763,6.34];
ybs[485]=['ν Psc',0.4477022,0.0977481,4.44];
ybs[486]=['',0.451016,0.6171214,5.64];
ybs[487]=['44 Cas',0.4575895,1.0587809,5.78];
ybs[488]=['',0.44879,-0.1956817,5.75];
ybs[489]=['107 Psc',0.45261,0.3557227,5.24];
ybs[490]=['',0.4469756,-0.6635739,6.17];
ybs[491]=['',0.4566247,0.7929872,6.34];
ybs[492]=['φ Per',0.4585105,0.8866466,4.07];
ybs[493]=['π Scl',0.4501183,-0.562241,5.25];
ybs[494]=['',0.4496018,-0.6408774,5.72];
ybs[495]=['',0.4616429,1.00616,6.21];
ybs[496]=['',0.4531937,-0.0624396,4.99];
ybs[497]=['',0.4475939,-0.8713707,6.64];
ybs[498]=['',0.463674,0.9983525,6.25];
ybs[499]=['',0.4586974,0.5638136,6.34];
ybs[500]=['',0.4617538,0.8072503,6.35];
ybs[501]=['',0.4475357,-1.0590035,5.71];
ybs[502]=['',0.4509301,-0.9359801,5.52];
ybs[503]=['',0.4583599,-0.0812116,6.19];
ybs[504]=['109 Psc',0.4632316,0.3524739,6.27];
ybs[505]=['τ Cet',0.4588414,-0.2761994,3.5];
ybs[506]=['ο Psc',0.4650552,0.16179,4.26];
ybs[507]=['',0.4772016,1.1163771,5.63];
ybs[508]=['',0.425043,-1.4461962,5.87];
ybs[509]=['',0.4673956,-0.0981113,5.34];
ybs[510]=['ε Scl',0.4655467,-0.4352928,5.31];
ybs[511]=['',0.4704212,0.3058666,6.55];
ybs[512]=['τ1 Hyi',0.4424204,-1.3794232,6.33];
ybs[513]=['',0.467111,-0.4753787,6.39];
ybs[514]=['',0.4764404,0.8088067,6.32];
ybs[515]=['',0.4667657,-0.884959,5.49];
ybs[516]=['',0.4666793,-0.9321799,5.04];
ybs[517]=['',0.4798839,0.6643431,5.94];
ybs[518]=['4 Ari',0.4773691,0.2978748,5.84];
ybs[519]=['',0.4799415,0.5724951,5.79];
ybs[520]=['',0.4722464,-0.7269005,6.18];
ybs[521]=['',0.420614,-1.4775172,5.69];
ybs[522]=['',0.4828947,0.8378988,5.82];
ybs[523]=['',0.4782269,0.0662686,5.91];
ybs[524]=['',0.4746308,-0.6466128,6.32];
ybs[525]=['',0.4904594,0.9083395,5.9];
ybs[526]=['1 Ari',0.4860328,0.3907127,5.86];
ybs[527]=['χ Cet',0.4829968,-0.1845741,4.67];
ybs[528]=['',0.4814459,-0.5403823,6.34];
ybs[529]=['1 Per',0.4951589,0.9644324,5.52];
ybs[530]=['',0.4889793,0.1946752,5.94];
ybs[531]=['',0.4833822,-0.6683364,6.37];
ybs[532]=['2 Per',0.4956772,0.8884276,5.79];
ybs[533]=['',0.4853471,-0.8326176,6.14];
ybs[534]=['',0.4987266,0.9003267,6.26];
ybs[535]=['ζ Cet',0.4911816,-0.1784495,3.73];
ybs[536]=['',0.5031477,0.9722879,6.45];
ybs[537]=['',0.4877375,-0.8743286,5.94];
ybs[538]=['ε Cas',0.5063077,1.1131669,3.38];
ybs[539]=['55 And',0.5002451,0.7127893,5.4];
ybs[540]=['α Tri',0.4990358,0.5181716,3.41];
ybs[541]=['γ1 Ari',0.5007666,0.3386964,4.83];
ybs[542]=['γ2 Ari',0.5007666,0.3386576,4.75];
ybs[543]=['',0.4971978,-0.2935457,5.8];
ybs[544]=['ω Cas',0.5139205,1.2006922,4.99];
ybs[545]=['ξ Psc',0.5005681,0.0575531,4.62];
ybs[546]=['τ2 Hyi',0.4696043,-1.3973975,6.06];
ybs[547]=['',0.5072735,0.712297,6.24];
ybs[548]=['',0.5074338,0.6499256,6.26];
ybs[549]=['β Ari',0.5056447,0.3650846,2.64];
ybs[550]=['',0.4989443,-0.6716831,6.1];
ybs[551]=['ψ Phe',0.4998183,-0.8062102,4.41];
ybs[552]=['',0.511582,0.6525295,5.89];
ybs[553]=['56 And',0.5126664,0.6520725,5.67];
ybs[554]=['φ Phe',0.5030901,-0.7397941,5.11];
ybs[555]=['7 Ari',0.5109883,0.4134101,5.74];
ybs[556]=['',0.5107617,0.0341936,6.01];
ybs[557]=['',0.5243504,1.0787306,6.02];
ybs[558]=['',0.5206075,0.7296052,6.78];
ybs[559]=['ι Ari',0.517422,0.312877,5.1];
ybs[560]=['',0.5192964,0.4871802,5.82];
ybs[561]=['56 Cet',0.5136581,-0.3912628,4.85];
ybs[562]=['χ Eri',0.5096563,-0.8988346,3.7];
ybs[563]=['',0.5293489,1.1297471,5.26];
ybs[564]=['3 Per',0.5236095,0.8606715,5.69];
ybs[565]=['λ Ari',0.5200646,0.41373,4.79];
ybs[566]=['η2 Hyi',0.5039755,-1.1787511,4.69];
ybs[567]=['',0.5082892,-1.0603201,6.06];
ybs[568]=['',0.5468059,1.3617708,6.04];
ybs[569]=['',0.5141902,-0.9015836,6.1];
ybs[570]=['',0.5151033,-0.8251196,4.83];
ybs[571]=['48 Cas',0.5404049,1.239439,4.54];
ybs[572]=['',0.521141,-0.5752238,6.35];
ybs[573]=['',0.5272861,0.3694296,5.87];
ybs[574]=['',0.5263901,0.2164766,6.09];
ybs[575]=['',0.5463919,1.2908083,6.23];
ybs[576]=['50 Cas',0.5471953,1.2658636,3.98];
ybs[577]=['47 Cas',0.5560787,1.3506774,5.38];
ybs[578]=['112 Psc',0.5293601,0.0559465,5.88];
ybs[579]=['57 Cet',0.5272041,-0.3615632,5.41];
ybs[580]=['',0.5170866,-1.1399749,6.37];
ybs[581]=['υ Cet',0.5282305,-0.3659858,4];
ybs[582]=['52 Cas',0.5435947,1.1346182,6];
ybs[583]=['',0.5304227,-0.1468766,5.51];
ybs[584]=['',0.5261046,-0.7316788,5.57];
ybs[585]=['53 Cas',0.5440946,1.1256922,5.58];
ybs[586]=['4 Per',0.5402436,0.9528649,5.04];
ybs[587]=['α Hyi',0.5212625,-1.0726968,2.86];
ybs[588]=['49 Cas',0.5573765,1.3303183,5.22];
ybs[589]=['σ Hyi',0.5053778,-1.3655235,6.16];
ybs[590]=['π For',0.5334247,-0.5217429,5.35];
ybs[591]=['α Psc',0.5376156,0.0501148,3.82];
ybs[592]=['α Psc',0.5376156,0.0501148,4.33];
ybs[593]=['',0.5775075,1.4207193,6.05];
ybs[594]=['',0.5514589,1.1381339,6.52];
ybs[595]=['ε Tri',0.5423137,0.5827892,5.5];
ybs[596]=['',0.5247998,-1.1511825,6.1];
ybs[597]=['',0.5401863,0.2370901,5.94];
ybs[598]=['χ Phe',0.5349794,-0.7785167,5.14];
ybs[599]=['γ1 And',0.5466795,0.7406641,2.26];
ybs[600]=['γ2 And',0.5467305,0.7406834,4.84];
ybs[601]=['10 Ari',0.5451285,0.4545333,5.63];
ybs[602]=['',0.5387655,-0.5158729,6.42];
ybs[603]=['60 Cet',0.5425749,0.0041149,5.43];
ybs[604]=['',0.5413127,-0.2652608,5.86];
ybs[605]=['',0.5451919,0.3204532,6.21];
ybs[606]=['61 Cet',0.54522,-0.0040669,5.93];
ybs[607]=['',0.5445859,-0.0697487,5.62];
ybs[608]=['ν For',0.5475869,-0.509459,4.69];
ybs[609]=['κ Ari',0.5577577,0.3971461,5.03];
ybs[610]=['',0.5558725,0.1458061,6.31];
ybs[611]=['11 Ari',0.5589452,0.4504888,6.15];
ybs[612]=['',0.556943,0.0024697,6.28];
ybs[613]=['α Ari',0.5604335,0.4113529,2];
ybs[614]=['',0.5683905,1.021531,5.67];
ybs[615]=['',0.5671296,0.7778113,6.42];
ybs[616]=['58 And',0.5665775,0.6626152,4.82];
ybs[617]=['',0.5744514,0.9415779,6.31];
ybs[618]=['β Tri',0.5710961,0.612485,3];
ybs[619]=['14 Ari',0.5703133,0.4545771,4.98];
ybs[620]=['',0.5699457,0.3024672,6.43];
ybs[621]=['',0.5664869,-0.3084623,6.1];
ybs[622]=['',0.5911375,1.2938494,6.29];
ybs[623]=['5 Per',0.5806229,1.0079417,6.36];
ybs[624]=['59 And',0.577068,0.6832026,5.63];
ybs[625]=['59 And',0.5771337,0.6832656,6.1];
ybs[626]=['',0.5699006,-0.4230716,6.48];
ybs[627]=['15 Ari',0.5754189,0.3421815,5.7];
ybs[628]=['',0.5674664,-0.7576632,5.85];
ybs[629]=['16 Ari',0.5780796,0.4545195,6.02];
ybs[630]=['5 Tri',0.5791774,0.5520724,6.23];
ybs[631]=['64 Cet',0.578342,0.1514037,5.63];
ybs[632]=['',0.5714909,-0.7628845,6.32];
ybs[633]=['',0.5726783,-0.8852143,6.12];
ybs[634]=['',0.5780317,-0.1736104,6.01];
ybs[635]=['63 Cet',0.5791936,-0.0300244,5.93];
ybs[636]=['55 Cas',0.5946439,1.1628864,6.07];
ybs[637]=['',0.590356,1.0239047,6.44];
ybs[638]=['6 Tri',0.5833143,0.5307163,4.94];
ybs[639]=['60 And',0.5874973,0.7738116,4.83];
ybs[640]=['',0.5842576,0.4236343,5.96];
ybs[641]=['',0.5894934,0.8930881,5.31];
ybs[642]=['η Ari',0.5849561,0.372025,5.27];
ybs[643]=['',0.5912314,0.8305741,6.06];
ybs[644]=['19 Ari',0.5859238,0.2685063,5.71];
ybs[645]=['ξ1 Cet',0.5855422,0.1562287,4.37];
ybs[646]=['66 Cet',0.5843956,-0.03995,5.54];
ybs[647]=['',0.5849586,-0.3646985,5.86];
ybs[648]=['μ For',0.5842409,-0.5344068,5.28];
ybs[649]=['',0.5996483,0.836275,6.33];
ybs[650]=['',0.6041139,0.9976061,6.48];
ybs[651]=['7 Tri',0.5989769,0.5840317,5.28];
ybs[652]=['20 Ari',0.598019,0.4518095,5.79];
ybs[653]=['21 Ari',0.5977667,0.4388944,5.58];
ybs[654]=['',0.5959357,-0.1633927,6.55];
ybs[655]=['',0.5909841,-0.716676,5.91];
ybs[656]=['δ Tri',0.6038822,0.5991276,4.87];
ybs[657]=['8 Per',0.6091511,1.0123382,5.75];
ybs[658]=['7 Per',0.6094572,1.0056522,5.98];
ybs[659]=['',0.6064489,0.7751024,6.7];
ybs[660]=['γ Tri',0.6050134,0.5925473,4.01];
ybs[661]=['',0.6041093,0.4166288,6.55];
ybs[662]=['67 Cet',0.6025964,-0.1102847,5.51];
ybs[663]=['π1 Hyi',0.5878213,-1.18224,5.55];
ybs[664]=['',0.6194688,1.1246813,6.6];
ybs[665]=['θ Ari',0.6081674,0.3491376,5.62];
ybs[666]=['62 And',0.6141206,0.8287278,5.3];
ybs[667]=['',0.6136481,0.8128895,6.21];
ybs[668]=['',0.6073114,0.0324775,5.58];
ybs[669]=['',0.6146284,0.856221,6.37];
ybs[670]=['φ Eri',0.5991251,-0.8972503,3.56];
ybs[671]=['10 Tri',0.6120001,0.5016989,5.03];
ybs[672]=['',0.6119249,0.406147,6.46];
ybs[673]=['',0.6152962,0.697041,6.63];
ybs[674]=['π2 Hyi',0.5931785,-1.1805835,5.69];
ybs[675]=['',0.6202816,0.8275129,6.11];
ybs[676]=['',0.6169513,0.5286727,6.47];
ybs[677]=['ο Cet',0.612975,-0.0501759,3.04];
ybs[678]=['63 And',0.6216507,0.8770883,5.59];
ybs[679]=['',0.6108178,-0.4510417,6.34];
ybs[680]=['',0.6144064,-0.0740548,6.5];
ybs[681]=['9 Per',0.6280679,0.9764623,5.17];
ybs[682]=['',0.6122097,-0.7285996,6.37];
ybs[683]=['',0.6294132,0.7242746,5.82];
ybs[684]=['',0.6136277,-0.9746302,5.81];
ybs[685]=['69 Cet',0.6243863,0.0086856,5.28];
ybs[686]=['',0.6346251,0.9680568,6.28];
ybs[687]=['70 Cet',0.6255068,-0.0136707,5.42];
ybs[688]=['',0.6244841,-0.186331,5.46];
ybs[689]=['',0.6245827,-0.3064874,5.87];
ybs[690]=['64 And',0.6367064,0.8745431,5.19];
ybs[691]=['κ For',0.6264317,-0.4139004,5.2];
ybs[692]=['10 Per',0.6408547,0.9897879,6.25];
ybs[693]=['',0.6284032,-0.3185739,6.22];
ybs[694]=['',0.624287,-0.7522059,6.31];
ybs[695]=['65 And',0.6420049,0.8792826,4.71];
ybs[696]=['',0.6284752,-0.6540606,6.53];
ybs[697]=['',0.6270116,-0.8899549,5.92];
ybs[698]=['ξ Ari',0.6371576,0.1869498,5.47];
ybs[699]=['',0.6341877,-0.44936,6.44];
ybs[700]=['71 Cet',0.6375354,-0.0467602,6.33];
ybs[701]=['δ Hyi',0.6202855,-1.1965533,4.09];
ybs[702]=['',0.6346927,-0.7110392,6.18];
ybs[703]=['ι Cas',0.6586635,1.1781299,4.52];
ybs[704]=['ρ Cet',0.6415758,-0.2127561,4.89];
ybs[705]=['66 And',0.6518159,0.8843506,6.12];
ybs[706]=['',0.6417421,-0.2659985,5.83];
ybs[707]=['',0.6476207,0.4732189,6.18];
ybs[708]=['11 Tri',0.6492727,0.5567842,5.54];
ybs[709]=['',0.6441605,-0.3480614,5.88];
ybs[710]=['λ Hor',0.6350041,-1.0508801,5.35];
ybs[711]=['κ Hyi',0.6240896,-1.2835877,5.01];
ybs[712]=['',0.65893,0.971026,6.51];
ybs[713]=['12 Tri',0.6522854,0.5195707,5.29];
ybs[714]=['ξ2 Cet',0.6516922,0.1493966,4.28];
ybs[715]=['',0.6508499,0.0359657,6.45];
ybs[716]=['13 Tri',0.6550961,0.5241485,5.89];
ybs[717]=['κ Eri',0.6449363,-0.8308404,4.25];
ybs[718]=['',0.636649,-1.1587925,6.41];
ybs[719]=['',0.6567347,0.4113446,6.19];
ybs[720]=['φ For',0.6500452,-0.5883721,5.14];
ybs[721]=['',0.6579597,0.1686887,6.07];
ybs[722]=['',0.6616361,0.5922417,6.25];
ybs[723]=['',0.6525808,-0.5411012,6.11];
ybs[724]=['',0.6625194,0.4421613,5.92];
ybs[725]=['26 Ari',0.6628081,0.3482671,6.15];
ybs[726]=['',0.6586393,-0.3941573,6.77];
ybs[727]=['27 Ari',0.663916,0.3107167,6.23];
ybs[728]=['',0.6628244,0.0061821,6];
ybs[729]=['',0.6599066,-0.4378552,6.51];
ybs[730]=['',0.6483905,-1.1204978,6.37];
ybs[731]=['',0.661361,-0.3917659,6.1];
ybs[732]=['14 Tri',0.6696971,0.6326061,5.15];
ybs[733]=['',0.6661369,0.0412928,5.25];
ybs[734]=['',0.673011,0.604594,5.83];
ybs[735]=['75 Cet',0.6689172,-0.0163457,5.35];
ybs[736]=['σ Cet',0.6682682,-0.2643515,4.75];
ybs[737]=['29 Ari',0.672557,0.2641192,6.04];
ybs[738]=['',0.6683376,-0.6340612,6.3];
ybs[739]=['',0.6991201,1.2725998,5.16];
ybs[740]=['λ1 For',0.6721925,-0.6030432,5.9];
ybs[741]=['',0.675049,-0.3473899,6.21];
ybs[742]=['',0.6845196,0.6939732,6.36];
ybs[743]=['',0.6958009,1.1491601,5.78];
ybs[744]=['',0.6852208,0.6529181,5.71];
ybs[745]=['ω For',0.6755765,-0.4910411,4.9];
ybs[746]=['15 Tri',0.6857066,0.6071073,5.35];
ybs[747]=['',0.6818262,0.1321013,6.18];
ybs[748]=['77 Cet',0.6799007,-0.1354698,5.75];
ybs[749]=['',0.6862344,0.1218949,5.82];
ybs[750]=['ν Cet',0.685301,0.0993183,4.86];
ybs[751]=['',0.6749048,-0.8900425,6.24];
ybs[752]=['',0.6909852,0.6777035,5.9];
ybs[753]=['',0.6896846,0.5533456,6.1];
ybs[754]=['',0.6912029,0.5997063,5.3];
ybs[755]=['80 Cet',0.6855215,-0.1349927,5.53];
ybs[756]=['',0.6927511,0.6980002,6.54];
ybs[757]=['',0.6914389,0.575761,6.25];
ybs[758]=['',0.6724473,-1.0906362,6.77];
ybs[759]=['31 Ari',0.6887733,0.2189415,5.68];
ybs[760]=['30 Ari',0.6905444,0.4318839,7.09];
ybs[761]=['30 Ari',0.6907482,0.431869,6.5];
ybs[762]=['',0.6884519,0.136601,5.81];
ybs[763]=['ι1 For',0.6855895,-0.5226841,5.75];
ybs[764]=['',0.6968146,0.6601354,6.18];
ybs[765]=['',0.697557,0.6664661,6.3];
ybs[766]=['',0.6946846,0.1359912,6.39];
ybs[767]=['81 Cet',0.6930318,-0.0575879,5.65];
ybs[768]=['λ2 For',0.6890176,-0.6018155,5.79];
ybs[769]=['ν Ari',0.6985666,0.3849766,5.43];
ybs[770]=['',0.7470213,1.4231567,5.78];
ybs[771]=['',0.6972056,0.0617724,6.21];
ybs[772]=['μ Hyi',0.6598964,-1.3789943,5.28];
ybs[773]=['ι2 For',0.6949857,-0.5253055,5.83];
ybs[774]=['η Hor',0.690038,-0.915361,5.31];
ybs[775]=['δ Cet',0.7009183,0.0074097,4.07];
ybs[776]=['',0.6951511,-0.6613786,6.49];
ybs[777]=['ε Cet',0.7009638,-0.2055353,4.84];
ybs[778]=['33 Ari',0.7068767,0.4739671,5.3];
ybs[779]=['',0.7044575,0.1083431,6.25];
ybs[780]=['',0.7038242,-0.1633168,5.78];
ybs[781]=['11 Per',0.7185157,0.9634291,5.77];
ybs[782]=['',0.7025123,-0.5329906,6.52];
ybs[783]=['',0.7181732,0.9358582,5.84];
ybs[784]=['12 Per',0.7141812,0.7031722,4.91];
ybs[785]=['',0.7009859,-0.7469274,4.75];
ybs[786]=['84 Cet',0.7085216,-0.010481,5.71];
ybs[787]=['',0.7277903,1.1854033,5.95];
ybs[788]=['',0.7180004,0.8440442,6.48];
ybs[789]=['μ Ari',0.714011,0.3509257,5.69];
ybs[790]=['ι Eri',0.7048944,-0.6939427,4.11];
ybs[791]=['',0.7109606,-0.0544231,6.05];
ybs[792]=['',0.7096319,-0.2522739,5.98];
ybs[793]=['',0.7142638,0.1891332,6.3];
ybs[794]=['',0.6981975,-1.1202554,6.55];
ybs[795]=['θ Per',0.7231636,0.8608404,4.12];
ybs[796]=['14 Per',0.7224025,0.7747723,5.43];
ybs[797]=['35 Ari',0.7189793,0.4852316,4.66];
ybs[798]=['ζ Hor',0.7040626,-0.9504085,5.21];
ybs[799]=['',0.7206683,0.4491153,6.35];
ybs[800]=['γ Cet',0.7176447,0.0581266,3.47];
ybs[801]=['',0.711244,-0.6682662,6.01];
ybs[802]=['ε Hyi',0.6978619,-1.1898068,4.11];
ybs[803]=['',0.711009,-0.8103456,6.1];
ybs[804]=['36 Ari',0.7224687,0.3116823,6.46];
ybs[805]=['ο Ari',0.7234055,0.2688816,5.77];
ybs[806]=['ι Hor',0.7125799,-0.8849755,5.41];
ybs[807]=['π Cet',0.7208011,-0.2402325,4.25];
ybs[808]=['38 Ari',0.7251258,0.2188609,5.18];
ybs[809]=['μ Cet',0.7249846,0.1781658,4.27];
ybs[810]=['',0.7165069,-0.7056869,6.36];
ybs[811]=['',0.746094,1.2169576,6.18];
ybs[812]=['',0.7266205,0.0838719,6.03];
ybs[813]=['',0.7211908,-0.5660234,6.22];
ybs[814]=['τ1 Eri',0.7249464,-0.3225115,4.47];
ybs[815]=['',0.7346454,0.6296593,6.25];
ybs[816]=['',0.7350074,0.6221781,6.3];
ybs[817]=['',0.719519,-0.9158827,6.15];
ybs[818]=['',0.7246696,-0.806225,6.85];
ybs[819]=['',0.7148687,-1.1627345,6.26];
ybs[820]=['39 Ari',0.7384931,0.5120813,4.51];
ybs[821]=['',0.7469627,0.9979163,6.25];
ybs[822]=['',0.7320535,-0.376055,6.49];
ybs[823]=['',0.7339179,-0.3908203,6.47];
ybs[824]=['40 Ari',0.7408918,0.3207264,5.82];
ybs[825]=['',0.7593178,1.2039255,5.8];
ybs[826]=['',0.7421021,0.4412301,5.86];
ybs[827]=['',0.7455331,0.6530744,6.45];
ybs[828]=['',0.7374564,-0.215856,6.9];
ybs[829]=['γ Hor',0.7240586,-1.110213,5.74];
ybs[830]=['η Per',0.7520372,0.9771635,3.76];
ybs[831]=['η1 For',0.7351137,-0.6188545,6.51];
ybs[832]=['π Ari',0.7441718,0.3064196,5.22];
ybs[833]=['ζ Hyi',0.7238463,-1.1784943,4.84];
ybs[834]=['41 Ari',0.7474888,0.4773942,3.63];
ybs[835]=['',0.7568956,1.0193788,6.45];
ybs[836]=['16 Per',0.7505256,0.6703894,4.23];
ybs[837]=['β For',0.7418921,-0.5639738,4.46];
ybs[838]=['',0.7557886,0.8191422,5.88];
ybs[839]=['17 Per',0.754449,0.6135052,4.53];
ybs[840]=['γ1 For',0.7454675,-0.4270481,6.14];
ybs[841]=['γ2 For',0.7455915,-0.4860647,5.39];
ybs[842]=['',0.7613021,0.9265737,6.36];
ybs[843]=['σ Ari',0.7537158,0.2648278,5.49];
ybs[844]=['η2 For',0.7468042,-0.6239816,5.92];
ybs[845]=['',0.7631372,0.8492817,6.26];
ybs[846]=['τ2 Eri',0.7507573,-0.36499,4.75];
ybs[847]=['η3 For',0.7486703,-0.621061,5.47];
ybs[848]=['ν Hor',0.7396833,-1.0945665,5.26];
ybs[849]=['',0.7490385,-0.6953352,6.36];
ybs[850]=['τ Per',0.76737,0.9224582,3.95];
ybs[851]=['20 Per',0.764179,0.6706987,5.33];
ybs[852]=['',0.7611837,0.2892754,6.31];
ybs[853]=['',0.7575194,-0.2212769,6.04];
ybs[854]=['',0.7542922,-0.5362172,6.4];
ybs[855]=['',0.7589474,-0.1631886,6.32];
ybs[856]=['',0.7755667,1.0753129,5.59];
ybs[857]=['',0.777966,1.1243776,6.24];
ybs[858]=['',0.7618353,-0.3889565,5.95];
ybs[859]=['ψ For',0.7611957,-0.6692654,5.92];
ybs[860]=['',0.7785664,0.8962321,6.22];
ybs[861]=['',0.7770559,0.824729,6.02];
ybs[862]=['',0.7540027,-1.0963865,6.03];
ybs[863]=['ρ2 Ari',0.7726459,0.3215176,5.91];
ybs[864]=['',0.7619479,-0.8691654,4];
ybs[865]=['ρ3 Ari',0.7753725,0.3161273,5.63];
ybs[866]=['',0.774206,0.1478548,5.97];
ybs[867]=['',0.7628717,-0.8862905,6.21];
ybs[868]=['ν Hyi',0.7433379,-1.3085561,4.75];
ybs[869]=['21 Per',0.7795547,0.5589159,5.11];
ybs[870]=['η Eri',0.7746079,-0.1537342,3.89];
ybs[871]=['',0.7755986,-0.0632257,5.17];
ybs[872]=['',0.7831005,0.6755132,6.04];
ybs[873]=['',0.7777961,0.0801208,6.11];
ybs[874]=['47 Ari',0.7826627,0.3622898,5.8];
ybs[875]=['π Per',0.7863075,0.6937955,4.7];
ybs[876]=['',0.7626221,-1.1230305,6.56];
ybs[877]=['',0.8257609,1.3876095,5.49];
ybs[878]=['24 Per',0.7874261,0.6156077,4.93];
ybs[879]=['4 Eri',0.7783873,-0.4149094,5.45];
ybs[880]=['',0.7774159,-0.5195115,6.29];
ybs[881]=['',0.7913731,0.8257009,5.47];
ybs[882]=['',0.7903181,0.7177051,5.89];
ybs[883]=['ε Ari',0.7875951,0.374005,4.63];
ybs[884]=['ε Ari',0.7875951,0.374005,4.63];
ybs[885]=['6 Eri',0.7814477,-0.410449,5.84];
ybs[886]=['',0.7962478,0.9152436,5.28];
ybs[887]=['',0.7963352,0.915253,6.74];
ybs[888]=['',0.7846838,-0.047013,5.23];
ybs[889]=['',0.7785145,-0.6650011,6.41];
ybs[890]=['',0.7925083,0.6670628,6.11];
ybs[891]=['',0.7848859,-0.1690798,6.14];
ybs[892]=['λ Cet',0.7894302,0.157009,4.7];
ybs[893]=['θ1 Eri',0.7815434,-0.7018954,3.24];
ybs[894]=['θ2 Eri',0.7815869,-0.7018906,4.35];
ybs[895]=['5 Eri',0.7889973,-0.0414782,5.56];
ybs[896]=['',0.785696,-0.5029727,6.14];
ybs[897]=['ζ For',0.7879631,-0.4395721,5.71];
ybs[898]=['',0.7939368,0.1912589,5.95];
ybs[899]=['',0.7878679,-0.5658129,6.31];
ybs[900]=['7 Eri',0.7940617,-0.0487049,6.11];
ybs[901]=['49 Ari',0.799507,0.4633814,5.9];
ybs[902]=['',0.8524918,1.4233813,5.95];
ybs[903]=['ρ1 Eri',0.7953121,-0.1322064,5.75];
ybs[904]=['',0.7987437,0.0946619,6.25];
ybs[905]=['β Hor',0.7820197,-1.1167041,4.99];
ybs[906]=['93 Cet',0.8009132,0.0774961,5.61];
ybs[907]=['α Cet',0.8004913,0.0729055,2.53];
ybs[908]=['',0.7985997,-0.1723299,5.83];
ybs[909]=['',0.799656,-0.1118268,6.19];
ybs[910]=['ε For',0.7967037,-0.4887605,5.89];
ybs[911]=['γ Per',0.813503,0.9353703,2.93];
ybs[912]=['',0.8065712,0.4949172,6.36];
ybs[913]=['ρ2 Eri',0.8020224,-0.1326096,5.32];
ybs[914]=['',0.8170254,0.9912057,4.76];
ybs[915]=['τ3 Eri',0.8001823,-0.410798,4.09];
ybs[916]=['',0.8175142,0.9800833,6.11];
ybs[917]=['ρ Per',0.8142957,0.6793963,3.39];
ybs[918]=['',0.8256514,1.1195082,5.89];
ybs[919]=['',0.8151201,0.7098026,6.05];
ybs[920]=['',0.8112934,0.278251,6.49];
ybs[921]=['ρ3 Eri',0.8088665,-0.1311466,5.26];
ybs[922]=['',0.8107066,0.0340364,6.05];
ybs[923]=['52 Ari',0.8149381,0.442292,6.8];
ybs[924]=['52 Ari',0.8149381,0.442292,7];
ybs[925]=['',0.8015442,-0.818345,5.82];
ybs[926]=['',0.827693,0.9127795,6.31];
ybs[927]=['',0.8187077,0.2316583,5.62];
ybs[928]=['',0.8483808,1.299868,4.87];
ybs[929]=['',0.8261601,0.8271733,6.41];
ybs[930]=['μ Hor',0.8035067,-1.0411014,5.11];
ybs[931]=['',0.8188773,-0.1047691,5.27];
ybs[932]=['β Per',0.8274659,0.716294,2.12];
ybs[933]=['ι Per',0.8318874,0.867394,4.05];
ybs[934]=['53 Ari',0.8233594,0.3135555,6.11];
ybs[935]=['θ Hyi',0.7955114,-1.253405,5.53];
ybs[936]=['54 Ari',0.8274259,0.3295187,6.27];
ybs[937]=['κ Per',0.8334748,0.7843814,3.8];
ybs[938]=['',0.8283914,0.1493262,6.28];
ybs[939]=['',0.823831,-0.4842558,6.19];
ybs[940]=['55 Ari',0.8332676,0.5089634,5.72];
ybs[941]=['',0.8355579,0.4870218,6.42];
ybs[942]=['',0.8368494,0.4708996,6.02];
ybs[943]=['ω Per',0.8410386,0.6928167,4.63];
ybs[944]=['',0.8372338,0.2086823,5.98];
ybs[945]=['',0.846508,0.8344271,6.33];
ybs[946]=['',0.8449777,0.7410591,6.15];
ybs[947]=['δ Ari',0.8417599,0.3457563,4.35];
ybs[948]=['',0.8403968,0.2291897,6.12];
ybs[949]=['',0.835927,-0.4128427,6.38];
ybs[950]=['56 Ari',0.8446639,0.4771798,5.79];
ybs[951]=['',0.8396866,-0.0650623,6.05];
ybs[952]=['',0.8507323,0.8422936,5.9];
ybs[953]=['',0.8391831,-0.2782296,6.26];
ybs[954]=['',0.8449013,0.117709,5.56];
ybs[955]=['',0.8203383,-1.2074204,6.15];
ybs[956]=['',0.8342155,-0.8491005,6.12];
ybs[957]=['',0.8869092,1.3581189,5.45];
ybs[958]=['94 Cet',0.8461311,-0.0194228,5.06];
ybs[959]=['α For',0.842202,-0.5044586,3.87];
ybs[960]=['',0.8619354,0.9987252,5.79];
ybs[961]=['',0.9514751,1.4832698,5.61];
ybs[962]=['',0.8571523,0.7432693,6.07];
ybs[963]=['',0.8705934,1.1473749,6.36];
ybs[964]=['',0.8430711,-0.7738134,5.93];
ybs[965]=['',0.8631741,0.8904589,5.03];
ybs[966]=['',0.8460792,-0.6258867,6.27];
ybs[967]=['',0.8583464,0.5347485,5.52];
ybs[968]=['ζ Ari',0.8560859,0.3687323,4.89];
ybs[969]=['',0.8622626,0.7928623,6.16];
ybs[970]=['',0.8489798,-0.5187329,6.16];
ybs[971]=['',0.8603669,0.574883,6.31];
ybs[972]=['',0.8615282,0.6068595,6.25];
ybs[973]=['',0.842667,-0.9989949,5.74];
ybs[974]=['',0.8638407,0.5631351,6.06];
ybs[975]=['',0.8668561,0.7079879,6.45];
ybs[976]=['',0.855089,-0.4540978,6.25];
ybs[977]=['',0.8151427,-1.3771291,5.57];
ybs[978]=['30 Per',0.8696658,0.7697971,5.47];
ybs[979]=['',0.8601383,-0.1018691,6.17];
ybs[980]=['ζ Eri',0.8592607,-0.1525016,4.8];
ybs[981]=['',0.8813741,1.1472454,4.84];
ybs[982]=['',0.8693014,0.6870398,5.96];
ybs[983]=['29 Per',0.8737291,0.8779525,5.15];
ybs[984]=['14 Eri',0.8625734,-0.1583491,6.14];
ybs[985]=['31 Per',0.8758908,0.8757285,5.03];
ybs[986]=['',0.8600793,-0.5366118,6.65];
ybs[987]=['',0.8732935,0.5987096,4.82];
ybs[988]=['95 Cet',0.8705728,-0.0148234,5.38];
ybs[989]=['',0.8682593,-0.5011854,5.91];
ybs[990]=['15 Eri',0.869887,-0.3914843,4.88];
ybs[991]=['59 Ari',0.8782503,0.4738809,5.9];
ybs[992]=['κ1 Cet',0.8750109,0.0602281,4.83];
ybs[993]=['',0.8714065,-0.3225172,5.71];
ybs[994]=['',0.8647144,-0.8320026,5.85];
ybs[995]=['',0.8801151,0.5083869,4.47];
ybs[996]=['60 Ari',0.8803674,0.4492973,6.12];
ybs[997]=['',0.8878378,0.8578332,5.93];
ybs[998]=['32 Per',0.8855816,0.7576307,4.95];
ybs[999]=['τ4 Eri',0.8749203,-0.3783397,3.69];
ybs[1000]=['',0.875115,-0.4196219,5.61];
ybs[1001]=['τ1 Ari',0.8837023,0.3704753,5.28];
ybs[1002]=['ζ1 Ret',0.8647681,-1.0907246,5.54];
ybs[1003]=['κ2 Cet',0.8826637,0.0655433,5.69];
ybs[1004]=['',0.8758184,-0.7503056,4.27];
ybs[1005]=['',0.9016981,1.1286038,5.23];
ybs[1006]=['ζ2 Ret',0.8667139,-1.0895255,5.24];
ybs[1007]=['',0.8937184,0.8603104,5.29];
ybs[1008]=['62 Ari',0.8881749,0.483226,5.52];
ybs[1009]=['',0.8801422,-0.4629728,6.39];
ybs[1010]=['',0.8650004,-1.1666765,6.05];
ybs[1011]=['τ2 Ari',0.8903534,0.3633955,5.09];
ybs[1012]=['',0.8830685,-0.411122,5.52];
ybs[1013]=['α Per',0.8985862,0.871608,1.79];
ybs[1014]=['',0.8867998,-0.4452059,6.35];
ybs[1015]=['',0.8984421,0.5866776,5.61];
ybs[1016]=['',0.9053896,0.9424667,6.51];
ybs[1017]=['',0.8826255,-0.8324738,6.39];
ybs[1018]=['64 Ari',0.8972815,0.4328867,5.5];
ybs[1019]=['',0.8937762,0.08658,6.38];
ybs[1020]=['',0.8918495,-0.1346569,6.2];
ybs[1021]=['ι Hyi',0.8527272,-1.3492447,5.52];
ybs[1022]=['',0.901699,0.7214361,6.51];
ybs[1023]=['65 Ari',0.8976971,0.3644593,6.08];
ybs[1024]=['',0.8962671,0.2217954,6.04];
ybs[1025]=['',0.9056718,0.8586756,6.09];
ybs[1026]=['ο Tau',0.8989711,0.1589491,3.6];
ybs[1027]=['',0.8929396,-0.5694743,6.5];
ybs[1028]=['',0.9280816,1.2555809,6.32];
ybs[1029]=['',0.9173524,1.0529939,6.49];
ybs[1030]=['',0.9148231,0.8576465,4.98];
ybs[1031]=['',0.9202659,1.0474862,4.21];
ybs[1032]=['',0.909058,0.328709,6.57];
ybs[1033]=['',0.9184701,0.8713508,5.58];
ybs[1034]=['ξ Tau',0.9092779,0.1712166,3.74];
ybs[1035]=['',0.90999,0.2236142,6.28];
ybs[1036]=['',0.9238438,1.0289503,4.54];
ybs[1037]=['',0.915275,0.5913902,5.61];
ybs[1038]=['χ1 For',0.9023315,-0.6255785,6.39];
ybs[1039]=['',0.9250954,1.0374567,6.13];
ybs[1040]=['34 Per',0.9206088,0.8654227,4.67];
ybs[1041]=['',0.9046276,-0.4754257,5.93];
ybs[1042]=['',0.9238742,0.9691433,5.09];
ybs[1043]=['',0.9207315,0.820548,6.24];
ybs[1044]=['66 Ari',0.9152658,0.3993455,6.03];
ybs[1045]=['',0.903205,-0.7253453,6.32];
ybs[1046]=['',0.9123324,-0.1956477,5.73];
ybs[1047]=['',0.9259654,0.8408866,5.82];
ybs[1048]=['σ Per',0.9257759,0.8389961,4.36];
ybs[1049]=['',0.8907649,-1.2138054,6.15];
ybs[1050]=['χ2 For',0.9094332,-0.6214116,5.71];
ybs[1051]=['',0.9499835,1.2814272,6.57];
ybs[1052]=['',0.9298591,0.8601851,6.29];
ybs[1053]=['χ3 For',0.9121934,-0.6244174,6.5];
ybs[1054]=['',0.9313076,0.8635181,6.39];
ybs[1055]=['',0.9196038,-0.1174403,5.99];
ybs[1056]=['4 Tau',0.9234613,0.1991804,5.14];
ybs[1057]=['',0.9191952,-0.2198857,5.59];
ybs[1058]=['',0.932631,0.8394786,5.47];
ybs[1059]=['',0.8976235,-1.2087848,5.96];
ybs[1060]=['',0.9281188,0.4825367,5.96];
ybs[1061]=['5 Tau',0.9255418,0.227107,4.11];
ybs[1062]=['',0.9248212,0.1093322,5.94];
ybs[1063]=['',0.9396602,1.0269397,6.4];
ybs[1064]=['36 Per',0.9337901,0.8051517,5.31];
ybs[1065]=['17 Eri',0.923874,-0.0872583,4.73];
ybs[1066]=['',0.9402233,1.0112985,6.37];
ybs[1067]=['',0.9346473,0.7841819,6.41];
ybs[1068]=['',0.9397709,0.9607864,5.98];
ybs[1069]=['',0.9342162,0.6202279,5.9];
ybs[1070]=['',0.9193696,-0.742778,5.78];
ybs[1071]=['',0.9207899,-0.7207165,6.12];
ybs[1072]=['',0.946277,1.0492006,6.46];
ybs[1073]=['',0.9384451,0.6976746,5.81];
ybs[1074]=['6 Tau',0.9329648,0.1649069,5.77];
ybs[1075]=['',0.9694851,1.3231552,6.27];
ybs[1076]=['',0.9221321,-0.8255311,5.99];
ybs[1077]=['',0.9287759,-0.4457387,6.38];
ybs[1078]=['κ Ret',0.9152326,-1.0971324,4.72];
ybs[1079]=['ε Eri',0.933823,-0.1637746,3.73];
ybs[1080]=['',0.9399545,0.3125351,6.17];
ybs[1081]=['7 Tau',0.9415212,0.428277,5.92];
ybs[1082]=['ψ Per',0.9516359,0.8423976,4.23];
ybs[1083]=['τ5 Eri',0.937163,-0.376265,4.27];
ybs[1084]=['',0.9425606,0.1133009,6.49];
ybs[1085]=['',0.9304744,-0.8779639,5.68];
ybs[1086]=['',0.941197,-0.1709483,6.25];
ybs[1087]=['',0.9210892,-1.1591411,5.83];
ybs[1088]=['',0.937491,-0.5411558,6.2];
ybs[1089]=['',0.9604167,0.9949245,6.3];
ybs[1090]=['',0.9401242,-0.5550262,6.4];
ybs[1091]=['',0.9306504,-1.0636391,6.41];
ybs[1092]=['',0.9578924,0.7444781,6.42];
ybs[1093]=['',0.9469928,-0.1940842,5.57];
ybs[1094]=['',0.9509686,0.0115332,5.71];
ybs[1095]=['20 Eri',0.9482182,-0.3035769,5.23];
ybs[1096]=['10 Tau',0.9513337,0.0082842,4.28];
ybs[1097]=['',0.9558421,0.270585,6.39];
ybs[1098]=['',0.9612917,0.3663069,6.5];
ybs[1099]=['',0.9366751,-1.1465098,6.75];
ybs[1100]=['',0.9779528,1.1045683,5.1];
ybs[1101]=['',0.9507835,-0.701653,4.58];
ybs[1102]=['',1.1294701,1.5102606,5.86];
ybs[1103]=['',0.9581291,-0.1277474,5.85];
ybs[1104]=['',0.9128681,-1.3661654,5.7];
ybs[1105]=['',0.9630015,0.2898728,6.16];
ybs[1106]=['21 Eri',0.9605045,-0.0969369,5.96];
ybs[1107]=['',0.9799314,1.0478896,5.76];
ybs[1108]=['',0.9712884,0.6571341,5.57];
ybs[1109]=['τ For',0.9587587,-0.4864386,6.01];
ybs[1110]=['12 Tau',0.9644125,0.0546043,5.57];
ybs[1111]=['',0.9632796,-0.0579676,6.23];
ybs[1112]=['',0.9621179,-0.1809096,6.19];
ybs[1113]=['11 Tau',0.9691826,0.4433249,6.11];
ybs[1114]=['',0.9648925,-0.0183078,6.12];
ybs[1115]=['',0.9653045,-0.2645071,6.33];
ybs[1116]=['22 Eri',0.9675852,-0.0896968,5.53];
ybs[1117]=['δ Per',0.9797244,0.8352739,3.01];
ybs[1118]=['40 Per',0.976545,0.5940308,4.97];
ybs[1119]=['',0.9955665,1.1740883,5.8];
ybs[1120]=['',0.9699547,-0.204762,6.49];
ybs[1121]=['13 Tau',0.9757122,0.3450654,5.69];
ybs[1122]=['',0.9849399,0.8481121,6.06];
ybs[1123]=['',0.9703156,-0.3405786,6.59];
ybs[1124]=['',0.9949582,1.1067772,4.8];
ybs[1125]=['',0.9872843,0.8058027,6.11];
ybs[1126]=['ο Per',0.9849554,0.5647521,3.83];
ybs[1127]=['14 Tau',0.982135,0.3444381,6.14];
ybs[1128]=['',0.9860455,0.6375596,5.59];
ybs[1129]=['δ For',0.9736491,-0.5561963,5];
ybs[1130]=['ν Per',0.9893,0.744344,3.77];
ybs[1131]=['δ Eri',0.9788258,-0.1691781,3.54];
ybs[1132]=['',0.9851563,0.3664868,6.1];
ybs[1133]=['',1.0105687,1.2381055,5.44];
ybs[1134]=['',0.9801764,-0.1817858,5.6];
ybs[1135]=['16 Tau',0.9867442,0.4251417,5.46];
ybs[1136]=['',0.9929653,0.7985008,5.66];
ybs[1137]=['17 Tau',0.9870503,0.4220674,3.7];
ybs[1138]=['',0.9759562,-0.650017,4.59];
ybs[1139]=['18 Tau',0.9883297,0.4347332,5.64];
ybs[1140]=['19 Tau',0.988519,0.4282412,4.3];
ybs[1141]=['24 Eri',0.9845983,-0.0190851,5.25];
ybs[1142]=['',1.0005649,0.9772191,6.1];
ybs[1143]=['γ Cam',1.0156512,1.2461442,4.63];
ybs[1144]=['20 Tau',0.9912147,0.4265007,3.87];
ybs[1145]=['25 Eri',0.9865166,-0.0039673,5.55];
ybs[1146]=['21 Tau',0.9915712,0.4297628,5.76];
ybs[1147]=['22 Tau',0.9921886,0.4292962,6.43];
ybs[1148]=['29 Tau',0.9899119,0.1067969,5.35];
ybs[1149]=['',0.941158,-1.3657109,6.29];
ybs[1150]=['',1.0105652,1.1448166,4.47];
ybs[1151]=['23 Tau',0.9933815,0.419176,4.18];
ybs[1152]=['',0.9813099,-0.7084371,6.45];
ybs[1153]=['',1.0105789,1.1059147,5.85];
ybs[1154]=['',0.9920527,0.1199363,5.91];
ybs[1155]=['',1.0034533,0.8867035,6.14];
ybs[1156]=['',1.0085157,0.9980709,6.46];
ybs[1157]=['π Eri',0.9913704,-0.2100126,4.42];
ybs[1158]=['',1.0005531,0.5876166,6.57];
ybs[1159]=['',1.0002199,0.5630953,6.25];
ybs[1160]=['η Tau',0.9984444,0.4219011,2.87];
ybs[1161]=['',1.0207261,1.1968333,6.32];
ybs[1162]=['',0.984074,-0.8376164,6.49];
ybs[1163]=['',0.9823512,-0.9460423,6.3];
ybs[1164]=['',0.9859532,-0.8253734,5.73];
ybs[1165]=['',1.006609,0.7684754,6.02];
ybs[1166]=['σ For',0.9921094,-0.5108463,5.9];
ybs[1167]=['',1.0021872,0.4099581,5.45];
ybs[1168]=['τ6 Eri',0.9940594,-0.4045882,4.23];
ybs[1169]=['30 Tau',1.0014179,0.1956715,5.07];
ybs[1170]=['β Ret',0.9794381,-1.1298747,3.85];
ybs[1171]=['',1.0107717,0.7860035,5.66];
ybs[1172]=['42 Per',1.0078074,0.5787263,5.11];
ybs[1173]=['27 Tau',1.0057621,0.4209859,3.63];
ybs[1174]=['',0.9959155,-0.520695,6.55];
ybs[1175]=['28 Tau',1.0058744,0.4224401,5.09];
ybs[1176]=['τ7 Eri',0.9975753,-0.4155029,5.24];
ybs[1177]=['',1.0027067,0.0051564,5.91];
ybs[1178]=['',1.0082082,0.4150182,6.17];
ybs[1179]=['ρ For',0.998513,-0.5253395,5.54];
ybs[1180]=['',1.008989,0.3894088,6.07];
ybs[1181]=['',0.9977765,-0.628977,6.21];
ybs[1182]=['',1.0017716,-0.3636451,5.81];
ybs[1183]=['',1.0108569,0.4476122,5.26];
ybs[1184]=['',1.0010283,-0.6554489,5.4];
ybs[1185]=['',1.0010648,-0.65542,4.73];
ybs[1186]=['',1.0181358,0.6008341,5.77];
ybs[1187]=['',1.0278241,1.0129917,5.8];
ybs[1188]=['',1.0163731,0.3856814,6.83];
ybs[1189]=['',1.014539,0.2288521,6.3];
ybs[1190]=['',1.0048706,-0.6306383,4.17];
ybs[1191]=['',1.0427159,1.2546364,6.34];
ybs[1192]=['',1.0187607,0.5451422,6.25];
ybs[1193]=['',1.0266236,0.8502509,5.76];
ybs[1194]=['31 Tau',1.017551,0.115206,5.67];
ybs[1195]=['',1.0099675,-0.6345748,6.86];
ybs[1196]=['',1.0229957,0.3035561,5.97];
ybs[1197]=['30 Eri',1.0201744,-0.0924255,5.48];
ybs[1198]=['ζ Per',1.0277875,0.5576095,2.85];
ybs[1199]=['',1.0447115,1.101924,5.03];
ybs[1200]=['',1.043162,1.06766,5];
ybs[1201]=['',1.0220117,-0.3205972,6.22];
ybs[1202]=['',1.0367236,0.836633,5.37];
ybs[1203]=['γ Hyi',0.9900815,-1.2945149,3.24];
ybs[1204]=['',1.033215,0.5429773,6.1];
ybs[1205]=['43 Per',1.0397308,0.8859136,5.28];
ybs[1206]=['32 Eri',1.0272194,-0.0504004,6.14];
ybs[1207]=['32 Eri',1.0272266,-0.0504343,4.79];
ybs[1208]=['τ8 Eri',1.0239308,-0.4284284,4.65];
ybs[1209]=['',1.0232239,-0.6050499,5.11];
ybs[1210]=['',1.038185,0.613397,5.49];
ybs[1211]=['',1.0221189,-0.8173048,5.93];
ybs[1212]=['',1.0311976,-0.2100427,6];
ybs[1213]=['32 Tau',1.0393446,0.3934295,5.63];
ybs[1214]=['',1.026151,-0.7032304,5.71];
ybs[1215]=['ε Per',1.0444482,0.6994157,2.89];
ybs[1216]=['33 Tau',1.0402227,0.4056016,6.06];
ybs[1217]=['',1.0419176,0.4280502,6.16];
ybs[1218]=['',1.0450421,0.6087299,6.53];
ybs[1219]=['',1.0394631,0.1065306,6.09];
ybs[1220]=['',1.0372102,-0.1690677,6.19];
ybs[1221]=['',1.0471421,0.6789901,6.3];
ybs[1222]=['',1.0260731,-0.9184884,6.46];
ybs[1223]=['ξ Per',1.0490789,0.6257683,4.04];
ybs[1224]=['',1.0523016,0.6786362,6.38];
ybs[1225]=['',1.1079968,1.4094491,5.1];
ybs[1226]=['γ Eri',1.0431814,-0.2346647,2.95];
ybs[1227]=['',1.047126,-0.0943716,5.83];
ybs[1228]=['',1.0511722,0.1813978,6.37];
ybs[1229]=['',1.0590965,0.646674,6.41];
ybs[1230]=['',1.0496337,-0.2183677,5.6];
ybs[1231]=['',1.0313173,-1.1065247,6.14];
ybs[1232]=['',1.0555053,0.3029665,6.32];
ybs[1233]=['',1.0564038,0.3186243,5.89];
ybs[1234]=['λ Tau',1.0556182,0.2190788,3.47];
ybs[1235]=['τ9 Eri',1.0510588,-0.4180753,4.66];
ybs[1236]=['',1.0835789,1.1997253,5.87];
ybs[1237]=['',1.0748772,1.0335074,5.06];
ybs[1238]=['',1.0602819,0.1755675,5.67];
ybs[1239]=['35 Eri',1.0588693,-0.0259721,5.28];
ybs[1240]=['',1.0437067,-0.9955249,6.05];
ybs[1241]=['',1.0540724,-0.5310816,5.93];
ybs[1242]=['δ Ret',1.0432946,-1.070535,4.56];
ybs[1243]=['',1.0854882,1.1445831,6.17];
ybs[1244]=['',1.0636169,-0.0036263,5.38];
ybs[1245]=['',1.0509845,-0.8988857,6.51];
ybs[1246]=['ν Tau',1.0662043,0.1055927,3.91];
ybs[1247]=['36 Tau',1.0721196,0.4217775,5.47];
ybs[1248]=['40 Tau',1.0687458,0.0959254,5.33];
ybs[1249]=['',1.0697122,0.1441238,5.46];
ybs[1250]=['',1.0837283,0.9436588,6.31];
ybs[1251]=['37 Tau',1.0734946,0.3864511,4.36];
ybs[1252]=['',1.0704981,0.0503931,5.36];
ybs[1253]=['',1.0664154,-0.3505213,6.46];
ybs[1254]=['',1.0672944,-0.3507703,7.01];
ybs[1255]=['',1.0901819,1.0888829,6.99];
ybs[1256]=['λ Per',1.0832832,0.8798287,4.29];
ybs[1257]=['39 Tau',1.0762927,0.3851708,5.9];
ybs[1258]=['',1.0697458,-0.2884763,6.39];
ybs[1259]=['γ Ret',1.052552,-1.0838027,4.51];
ybs[1260]=['',1.0708993,-0.2222188,5.61];
ybs[1261]=['ι Ret',1.0544809,-1.0649469,4.97];
ybs[1262]=['',1.0719528,-0.3546771,6.13];
ybs[1263]=['41 Tau',1.0820651,0.4827432,5.2];
ybs[1264]=['ψ Tau',1.0838796,0.5071985,5.23];
ybs[1265]=['',1.0968599,1.0465987,6.28];
ybs[1266]=['',0.9545353,-1.4817771,6.41];
ybs[1267]=['',1.077856,-0.1535292,6.26];
ybs[1268]=['48 Per',1.0921547,0.8337541,4.04];
ybs[1269]=['',1.0767239,-0.3569649,6.34];
ybs[1270]=['',1.0757602,-0.4815701,5.59];
ybs[1271]=['',1.0958847,0.9579519,6.18];
ybs[1272]=['49 Per',1.089747,0.6594918,6.09];
ybs[1273]=['50 Per',1.0913142,0.6649333,5.51];
ybs[1274]=['',1.0863514,0.2656639,6.01];
ybs[1275]=['',1.0876969,0.3036562,5.89];
ybs[1276]=['',1.1183084,1.2598098,6.03];
ybs[1277]=['',1.1133163,1.196555,6.32];
ybs[1278]=['ω1 Tau',1.0929181,0.3432555,5.5];
ybs[1279]=['',1.0920781,0.2348574,5.95];
ybs[1280]=['',1.0828057,-0.7480138,6.59];
ybs[1281]=['',1.1014473,0.5871929,5.72];
ybs[1282]=['44 Tau',1.1004619,0.4631743,5.41];
ybs[1283]=['',1.0922247,-0.2849756,5.37];
ybs[1284]=['',1.1940682,1.4635501,5.57];
ybs[1285]=['37 Eri',1.0972637,-0.119843,5.44];
ybs[1286]=['',1.0876094,-0.7994717,6.59];
ybs[1287]=['45 Tau',1.1018914,0.0973886,5.72];
ybs[1288]=['',1.0990305,-0.152935,5.7];
ybs[1289]=['',1.0803563,-1.1198633,6.38];
ybs[1290]=['',1.1074811,0.3025319,6.09];
ybs[1291]=['',1.1208927,1.0038297,6.08];
ybs[1292]=['',1.1091198,0.3921709,6.12];
ybs[1293]=['ο1 Eri',1.1037673,-0.1183479,4.04];
ybs[1294]=['',1.0978193,-0.614646,6.44];
ybs[1295]=['',1.1021365,-0.3542895,5.79];
ybs[1296]=['',1.1148175,0.6636173,6.45];
ybs[1297]=['δ Hor',1.0977935,-0.7319276,4.93];
ybs[1298]=['μ Per',1.1194504,0.8458649,4.14];
ybs[1299]=['',1.2005773,1.455381,5.46];
ybs[1300]=['',1.1296316,1.0804331,5.7];
ybs[1301]=['52 Per',1.1188696,0.7075336,4.71];
ybs[1302]=['',1.1118368,0.1792107,6.23];
ybs[1303]=['',1.1115359,0.1561439,6.51];
ybs[1304]=['46 Tau',1.1116258,0.1356456,5.29];
ybs[1305]=['',1.1130256,0.2235593,6.25];
ybs[1306]=['47 Tau',1.1133745,0.1626513,4.84];
ybs[1307]=['',1.1116921,-0.0190927,6.44];
ybs[1308]=['',1.1301468,1.0107978,5.71];
ybs[1309]=['',1.1278548,0.9366495,5.19];
ybs[1310]=['',1.1163029,0.175692,5.22];
ybs[1311]=['',1.1049798,-0.7733885,6.71];
ybs[1312]=['',1.1827592,1.4114928,5.43];
ybs[1313]=['39 Eri',1.1146856,-0.1780403,4.87];
ybs[1314]=['48 Tau',1.1215889,0.2697453,6.32];
ybs[1315]=['μ Tau',1.120322,0.1561557,4.29];
ybs[1316]=['',1.1197698,0.1091637,6.93];
ybs[1317]=['',1.1200166,0.1089353,6.31];
ybs[1318]=['',1.1098877,-0.7034004,6.37];
ybs[1319]=['',1.1342149,0.8787548,4.61];
ybs[1320]=['ο2 Eri',1.1186016,-0.1326062,4.43];
ybs[1321]=['α Hor',1.1115433,-0.737205,3.86];
ybs[1322]=['',1.1467231,1.1378261,5.27];
ybs[1323]=['',1.1331561,0.7364344,6.22];
ybs[1324]=['ω2 Tau',1.128282,0.3601064,4.94];
ybs[1325]=['',1.1384677,0.8744413,5.45];
ybs[1326]=['51 Tau',1.1332394,0.3775596,5.65];
ybs[1327]=['',1.1275792,-0.1120144,5.94];
ybs[1328]=['',1.1427795,0.8896511,5.55];
ybs[1329]=['',1.1328822,0.166506,6.54];
ybs[1330]=['',1.1508801,1.0609353,5.39];
ybs[1331]=['α Ret',1.1114157,-1.0894045,3.35];
ybs[1332]=['',1.1423203,0.7306033,5.92];
ybs[1333]=['γ Dor',1.1196904,-0.8976563,4.25];
ybs[1334]=['53 Tau',1.1377979,0.3699245,5.35];
ybs[1335]=['',1.1131443,-1.0844869,5.45];
ybs[1336]=['56 Tau',1.1385943,0.3809427,5.38];
ybs[1337]=['',1.1506401,0.9871272,5.88];
ybs[1338]=['54 Per',1.1426544,0.6042162,4.93];
ybs[1339]=['',1.1414623,0.5586072,6.16];
ybs[1340]=['',1.131179,-0.36062,6];
ybs[1341]=['γ Tau',1.1391498,0.2736715,3.65];
ybs[1342]=['υ4 Eri',1.1289905,-0.5889537,3.56];
ybs[1343]=['φ Tau',1.142073,0.478277,4.95];
ybs[1344]=['',1.1382144,0.1775735,6.31];
ybs[1345]=['53 Per',1.1483732,0.8124621,4.85];
ybs[1346]=['57 Tau',1.1398254,0.2458806,5.59];
ybs[1347]=['',1.1558425,1.0413921,6.19];
ybs[1348]=['',1.1326517,-0.3999753,6.07];
ybs[1349]=['',1.1419958,0.328033,6.12];
ybs[1350]=['ε Ret',1.1208333,-1.0340616,4.44];
ybs[1351]=['58 Tau',1.1426723,0.2643755,5.26];
ybs[1352]=['',1.1200395,-1.0628,6.37];
ybs[1353]=['',1.1438269,0.2428863,6.17];
ybs[1354]=['',1.1340217,-0.5908254,6.37];
ybs[1355]=['',1.1427062,0.107916,5.77];
ybs[1356]=['',1.1433828,0.1619279,6.53];
ybs[1357]=['',1.1420937,-0.1080919,6.27];
ybs[1358]=['',1.1423444,-0.1316011,5.85];
ybs[1359]=['',1.1344114,-0.7716962,5.34];
ybs[1360]=['',1.13107,-0.9216479,6.09];
ybs[1361]=['',1.1458273,-0.000805,5.86];
ybs[1362]=['',1.1415849,-0.359317,5.38];
ybs[1363]=['60 Tau',1.1489774,0.2465945,5.72];
ybs[1364]=['χ Tau',1.1517183,0.4482089,5.37];
ybs[1365]=['',1.1506396,0.3642994,5.91];
ybs[1366]=['',1.1570249,0.7413952,6.23];
ybs[1367]=['θ Ret',1.1254196,-1.1030748,5.87];
ybs[1368]=['δ1 Tau',1.152934,0.3070673,3.76];
ybs[1369]=['',1.1451777,-0.4481374,6.01];
ybs[1370]=['',1.1557092,0.3670963,5.99];
ybs[1371]=['63 Tau',1.1550069,0.2937065,5.64];
ybs[1372]=['55 Per',1.1604308,0.5965692,5.73];
ybs[1373]=['62 Tau',1.1578386,0.4250177,6.36];
ybs[1374]=['56 Per',1.1610187,0.5935864,5.76];
ybs[1375]=['δ2 Tau',1.1580003,0.305336,4.8];
ybs[1376]=['66 Tau',1.1566941,0.1660077,5.12];
ybs[1377]=['',1.1732765,1.0059073,6.32];
ybs[1378]=['ξ Eri',1.1554228,-0.0644852,5.17];
ybs[1379]=['',1.1520753,-0.4335582,5.83];
ybs[1380]=['',1.1617939,0.333215,5.98];
ybs[1381]=['',1.1517213,-0.6194841,6.39];
ybs[1382]=['κ1 Tau',1.1637426,0.3899732,4.22];
ybs[1383]=['κ2 Tau',1.1639498,0.3883292,5.28];
ybs[1384]=['δ3 Tau',1.1641002,0.3137741,4.29];
ybs[1385]=['',1.1673626,0.5495767,5.28];
ybs[1386]=['70 Tau',1.1645988,0.2790895,6.46];
ybs[1387]=['υ Tau',1.167867,0.3990358,4.28];
ybs[1388]=['43 Eri',1.1557611,-0.5928225,3.96];
ybs[1389]=['71 Tau',1.1677508,0.2734544,4.49];
ybs[1390]=['η Ret',1.1437691,-1.1053947,5.24];
ybs[1391]=['π Tau',1.1688514,0.2576619,4.69];
ybs[1392]=['',1.1675144,0.1507919,6.06];
ybs[1393]=['',1.1596623,-0.6057603,6.55];
ybs[1394]=['72 Tau',1.1721667,0.4022172,5.53];
ybs[1395]=['',1.1701638,0.037151,6.23];
ybs[1396]=['',1.2049324,1.266657,5.94];
ybs[1397]=['',1.1725333,0.1965483,5.88];
ybs[1398]=['',1.1752601,0.3781885,5.72];
ybs[1399]=['',1.1607572,-0.7698773,6.39];
ybs[1400]=['',1.1547528,-0.9951981,6.29];
ybs[1401]=['',1.1793601,0.5307466,6.4];
ybs[1402]=['75 Tau',1.1769151,0.2863757,4.97];
ybs[1403]=['76 Tau',1.1766354,0.2581213,5.9];
ybs[1404]=['ε Tau',1.1777957,0.335602,3.53];
ybs[1405]=['',1.1689234,-0.4194402,6.11];
ybs[1406]=['θ1 Tau',1.1774891,0.2794368,3.84];
ybs[1407]=['θ2 Tau',1.1778639,0.2778411,3.4];
ybs[1408]=['',1.1747304,0.0332875,6.15];
ybs[1409]=['79 Tau',1.1785215,0.2285631,5.03];
ybs[1410]=['',1.1767862,0.0249445,5.55];
ybs[1411]=['',1.1580541,-1.0679322,5.94];
ybs[1412]=['1 Cam',1.1947743,0.9417323,5.77];
ybs[1413]=['',1.1683253,-0.8185293,6.1];
ybs[1414]=['',1.1872032,0.5673245,6.21];
ybs[1415]=['',1.1822617,0.1844714,6.79];
ybs[1416]=['',1.176527,-0.3387678,5.96];
ybs[1417]=['80 Tau',1.1843215,0.2737653,5.58];
ybs[1418]=['',1.1788005,-0.2268966,5.6];
ybs[1419]=['',1.1909411,0.6991286,6.26];
ybs[1420]=['',1.1836705,0.1799452,6.48];
ybs[1421]=['δ Men',1.1193129,-1.3990489,5.69];
ybs[1422]=['',1.1861685,0.2834627,4.78];
ybs[1423]=['81 Tau',1.1865275,0.2747014,5.48];
ybs[1424]=['',1.1732762,-0.7314906,6.44];
ybs[1425]=['83 Tau',1.1863368,0.2403623,5.4];
ybs[1426]=['',1.1833463,-0.2363984,6.24];
ybs[1427]=['85 Tau',1.1918363,0.2774783,6.02];
ybs[1428]=['',1.1780936,-0.8110051,6.16];
ybs[1429]=['57 Per',1.199933,0.752406,6.09];
ybs[1430]=['',1.1695124,-1.0903393,5.75];
ybs[1431]=['',1.1923978,0.0952353,6.39];
ybs[1432]=['45 Eri',1.191324,0.000049,4.91];
ybs[1433]=['',1.1888847,-0.2373308,6.21];
ybs[1434]=['',1.1846092,-0.6214458,5.96];
ybs[1435]=['',1.2153171,1.1223481,5.94];
ybs[1436]=['',1.1944624,-0.055202,5.81];
ybs[1437]=['',1.1992742,0.3152497,6.25];
ybs[1438]=['δ Cae',1.1847531,-0.7837672,5.07];
ybs[1439]=['ρ Tau',1.200461,0.2598813,4.65];
ybs[1440]=['',1.2044751,0.5062561,5.88];
ybs[1441]=['',1.2000551,0.1650864,6.01];
ybs[1442]=['',1.1974322,-0.1874461,6.06];
ybs[1443]=['',1.2013852,0.0979853,5.68];
ybs[1444]=['46 Eri',1.1999644,-0.1168186,5.72];
ybs[1445]=['',1.2013715,-0.1185474,6.09];
ybs[1446]=['47 Eri',1.201132,-0.1428701,5.11];
ybs[1447]=['',1.2011125,-0.1557662,5.26];
ybs[1448]=['υ1 Eri',1.197282,-0.5187246,4.51];
ybs[1449]=['58 Per',1.2141163,0.7209757,4.25];
ybs[1450]=['',1.2087903,0.3477807,6.36];
ybs[1451]=['ν Men',1.1303626,-1.4229157,5.79];
ybs[1452]=['α Tau',1.2095662,0.2889177,0.85];
ybs[1453]=['88 Tau',1.2081703,0.178121,4.25];
ybs[1454]=['',1.2123158,0.4081475,6.02];
ybs[1455]=['',1.2056145,-0.1691514,6.37];
ybs[1456]=['',1.2042438,-0.3468914,6.13];
ybs[1457]=['',1.2092971,-0.0622621,6.33];
ybs[1458]=['ν Eri',1.210579,-0.0577366,3.93];
ybs[1459]=['υ2 Eri',1.2061413,-0.5326278,3.82];
ybs[1460]=['α Dor',1.1976591,-0.9599174,3.27];
ybs[1461]=['2 Cam',1.2294036,0.9340209,5.35];
ybs[1462]=['3 Cam',1.2291165,0.9271564,5.05];
ybs[1463]=['',1.2617073,1.3377955,6.49];
ybs[1464]=['',1.214705,0.0181915,5.31];
ybs[1465]=['',1.2212318,0.4709467,6.47];
ybs[1466]=['',1.2199598,0.3617738,5.92];
ybs[1467]=['89 Tau',1.2193092,0.280593,5.79];
ybs[1468]=['90 Tau',1.2191822,0.2191138,4.27];
ybs[1469]=['51 Eri',1.2162098,-0.0424037,5.23];
ybs[1470]=['',1.1947355,-1.0956749,5.79];
ybs[1471]=['',1.2117962,-0.535335,6.3];
ybs[1472]=['',1.2250473,0.44089,6.22];
ybs[1473]=['σ1 Tau',1.2236499,0.2765068,5.07];
ybs[1474]=['σ2 Tau',1.2241855,0.278571,4.69];
ybs[1475]=['',1.2231344,0.1381223,5.39];
ybs[1476]=['53 Eri',1.218299,-0.2488906,3.87];
ybs[1477]=['',1.2351938,0.8437362,5.67];
ybs[1478]=['',1.2214934,-0.2108343,5.01];
ybs[1479]=['93 Tau',1.2274548,0.213633,5.46];
ybs[1480]=['',1.1360253,-1.4459533,6.76];
ybs[1481]=['',1.244685,1.0395437,6.5];
ybs[1482]=['',1.2233066,-0.2498656,5.45];
ybs[1483]=['',1.2257943,-0.0176299,6.1];
ybs[1484]=['',1.2364126,0.6688413,5.99];
ybs[1485]=['',1.2336908,0.5001555,5.78];
ybs[1486]=['',1.2736689,1.3260765,6.06];
ybs[1487]=['',1.2087717,-1.0826813,5.4];
ybs[1488]=['',1.2438868,0.8729187,5.87];
ybs[1489]=['59 Per',1.2413906,0.7575764,5.29];
ybs[1490]=['',1.2263213,-0.4265577,5.58];
ybs[1491]=['54 Eri',1.2279547,-0.342596,4.32];
ybs[1492]=['τ Tau',1.2374329,0.4013958,4.28];
ybs[1493]=['',1.2201175,-0.9011065,6.44];
ybs[1494]=['95 Tau',1.2417807,0.421143,6.13];
ybs[1495]=['',1.246941,0.7125693,6.08];
ybs[1496]=['',1.2447023,0.5743143,6.45];
ybs[1497]=['α Cae',1.2273706,-0.7299233,4.45];
ybs[1498]=['β Cae',1.2341823,-0.6475672,5.05];
ybs[1499]=['',1.2246817,-1.0280168,6.53];
ybs[1500]=['55 Eri',1.2420493,-0.1527668,6.82];
ybs[1501]=['55 Eri',1.2420856,-0.1528105,6.7];
ybs[1502]=['',1.2465004,0.1952386,5.4];
ybs[1503]=['56 Eri',1.2442923,-0.1477099,5.9];
ybs[1504]=['',1.2392997,-0.5362447,5.68];
ybs[1505]=['',1.2790997,1.238806,6.37];
ybs[1506]=['4 Cam',1.2648694,0.9912674,5.3];
ybs[1507]=['',1.2525779,0.4130776,6.35];
ybs[1508]=['',1.244099,-0.3250887,5.53];
ybs[1509]=['',1.2579352,0.7042707,5.97];
ybs[1510]=['',1.2652214,0.9711128,6.26];
ybs[1511]=['λ Pic',1.2363625,-0.8803462,5.31];
ybs[1512]=['',1.2548688,0.3276726,6.01];
ybs[1513]=['',1.2412665,-0.7160038,6.25];
ybs[1514]=['',1.2534965,0.2049883,5.37];
ybs[1515]=['μ Eri',1.2506505,-0.0561125,4.02];
ybs[1516]=['',1.2480779,-0.3707714,5.72];
ybs[1517]=['',1.2545883,-0.0508798,6.33];
ybs[1518]=['',1.3292291,1.4176436,5.07];
ybs[1519]=['',1.2507877,-0.5928077,6.86];
ybs[1520]=['',1.2537114,-0.4895337,6.19];
ybs[1521]=['',1.2509077,-0.6862125,6.05];
ybs[1522]=['',1.2837042,1.1090049,5.44];
ybs[1523]=['',1.2687477,0.5694307,5.86];
ybs[1524]=['',1.2682384,0.549341,5.58];
ybs[1525]=['κ Dor',1.2422079,-1.0418265,5.27];
ybs[1526]=['',1.2102751,-1.354588,6.05];
ybs[1527]=['58 Eri',1.2593128,-0.2948872,5.51];
ybs[1528]=['',1.2716045,0.6549461,4.88];
ybs[1529]=['',1.2650436,0.0632917,6.03];
ybs[1530]=['',1.2778198,0.8513272,5.66];
ybs[1531]=['',1.2641024,-0.0983629,5.78];
ybs[1532]=['96 Tau',1.2698343,0.2782341,6.08];
ybs[1533]=['59 Eri',1.2634229,-0.2843365,5.77];
ybs[1534]=['ζ Cae',1.2597116,-0.5232795,6.37];
ybs[1535]=['',1.2443256,-1.1028647,6.46];
ybs[1536]=['μ Men',1.2341721,-1.2372602,5.54];
ybs[1537]=['α Cam',1.2928632,1.1585101,4.29];
ybs[1538]=['π3 Ori',1.2699528,0.1221526,3.19];
ybs[1539]=['π2 Ori',1.2733921,0.1559855,4.36];
ybs[1540]=['',1.2685906,-0.2396716,6.26];
ybs[1541]=['',1.2869659,0.9228669,6.41];
ybs[1542]=['97 Tau',1.2771092,0.3294541,5.1];
ybs[1543]=['',1.2621191,-0.7669285,6.72];
ybs[1544]=['60 Eri',1.2706327,-0.2823929,5.03];
ybs[1545]=['',1.2845251,0.7439022,5.71];
ybs[1546]=['2 Aur',1.283448,0.6412157,4.78];
ybs[1547]=['π4 Ori',1.2758663,0.0984667,3.69];
ybs[1548]=['',1.2782835,0.1747327,6.11];
ybs[1549]=['',1.2836572,0.4875288,5.97];
ybs[1550]=['5 Cam',1.2954698,0.9650572,5.52];
ybs[1551]=['ο1 Ori',1.2819818,0.2493476,4.74];
ybs[1552]=['',1.26974,-0.720533,6.07];
ybs[1553]=['',1.2936135,0.7696123,6.08];
ybs[1554]=['',1.2753519,-0.6085911,5.86];
ybs[1555]=['ω Eri',1.2828279,-0.0945429,4.39];
ybs[1556]=['',1.2998726,0.9233392,5.75];
ybs[1557]=['5 Ori',1.2852359,0.044395,5.33];
ybs[1558]=['ι Pic',1.2715933,-0.9324306,5.61];
ybs[1559]=['ι Pic',1.2716734,-0.9324016,6.42];
ybs[1560]=['',1.2876015,0.0280082,6.61];
ybs[1561]=['',1.2928324,0.3406883,6.37];
ybs[1562]=['π5 Ori',1.2890369,0.043209,3.72];
ybs[1563]=['7 Cam',1.3050654,0.9387358,4.47];
ybs[1564]=['6 Ori',1.2916845,0.2000314,5.19];
ybs[1565]=['π1 Ori',1.2921378,0.1777725,4.65];
ybs[1566]=['',1.2916124,0.1363802,5.33];
ybs[1567]=['',1.3317541,1.2967725,6.06];
ybs[1568]=['',1.2995619,0.6318589,6.07];
ybs[1569]=['',1.2915535,0.0087673,5.99];
ybs[1570]=['',1.2986735,0.4298091,6.37];
ybs[1571]=['',1.2964356,0.2630957,5.81];
ybs[1572]=['ι Aur',1.3022869,0.5794446,2.69];
ybs[1573]=['',1.2966584,0.0948306,6.5];
ybs[1574]=['',1.2920754,-0.2915717,5.7];
ybs[1575]=['ο2 Ori',1.2987034,0.236465,4.07];
ybs[1576]=['',1.2929463,-0.28594,5.72];
ybs[1577]=['62 Eri',1.2981454,-0.089664,5.51];
ybs[1578]=['',1.2933964,-0.4484315,6.72];
ybs[1579]=['',1.2900962,-0.6910406,6.1];
ybs[1580]=['',1.3032114,0.299971,5.48];
ybs[1581]=['99 Tau',1.3054172,0.4185618,5.79];
ybs[1582]=['',1.3397316,1.2879365,6.66];
ybs[1583]=['8 Cam',1.3158545,0.9282992,6.08];
ybs[1584]=['',1.3418308,1.2932215,5.96];
ybs[1585]=['98 Tau',1.306972,0.4377863,5.81];
ybs[1586]=['',1.3021509,-0.018041,6.23];
ybs[1587]=['ω Aur',1.3124296,0.6618757,4.94];
ybs[1588]=['',1.3248386,1.0665549,6.03];
ybs[1589]=['',1.3313789,1.1668065,6.19];
ybs[1590]=['',1.3036628,-0.2478028,6.15];
ybs[1591]=['',1.3060066,-0.038038,6.35];
ybs[1592]=['',1.2882736,-1.0212306,6.12];
ybs[1593]=['',1.2808489,-1.1630828,6.41];
ybs[1594]=['5 Aur',1.317097,0.6881235,5.95];
ybs[1595]=['',1.31017,0.2543887,6.09];
ybs[1596]=['π6 Ori',1.3077583,0.0304917,4.47];
ybs[1597]=['6 Aur',1.3174696,0.6926654,6.58];
ybs[1598]=['β Cam',1.3326886,1.0554408,4.03];
ybs[1599]=['',1.3091404,-0.285242,5.66];
ybs[1600]=['ε Aur',1.3246514,0.7654017,2.99];
ybs[1601]=['',1.2772908,-1.2631184,6.28];
ybs[1602]=['',1.3117621,-0.2578458,7.71];
ybs[1603]=['63 Eri',1.3129486,-0.1785666,5.38];
ybs[1604]=['',1.3165336,0.0636539,7.03];
ybs[1605]=['',1.3166282,0.0636682,6.66];
ybs[1606]=['64 Eri',1.313254,-0.218259,4.79];
ybs[1607]=['ζ Aur',1.3266909,0.7174442,3.75];
ybs[1608]=['',1.3168473,-0.0354962,6.32];
ybs[1609]=['',1.3173795,-0.0998614,6.22];
ybs[1610]=['',1.3303454,0.7238215,6.14];
ybs[1611]=['',1.483454,1.496374,6.51];
ybs[1612]=['ψ Eri',1.3200388,-0.1246605,4.81];
ybs[1613]=['',1.3220772,0.0131487,5.92];
ybs[1614]=['',1.3228154,0.0286224,6.24];
ybs[1615]=['ι Tau',1.328373,0.3773477,4.64];
ybs[1616]=['',1.3194818,-0.3494243,4.91];
ybs[1617]=['11 Cam',1.344367,1.0297642,5.08];
ybs[1618]=['12 Cam',1.3446432,1.0306121,6.08];
ybs[1619]=['',1.3462389,1.0681142,6.04];
ybs[1620]=['',1.3259087,-0.0729383,5.85];
ybs[1621]=['',1.3337878,0.5327533,6.14];
ybs[1622]=['',1.3355102,0.5646117,6.62];
ybs[1623]=['',1.3224231,-0.4580436,5.02];
ybs[1624]=['η Men',1.2852384,-1.3072834,5.47];
ybs[1625]=['',1.344784,0.9500594,7.24];
ybs[1626]=['',1.3191436,-0.6926629,6.03];
ybs[1627]=['',1.335354,0.4839048,6.6];
ybs[1628]=['',1.3338804,0.3718915,6.19];
ybs[1629]=['1 Lep',1.3251305,-0.3973117,5.75];
ybs[1630]=['',1.323102,-0.5539754,5.94];
ybs[1631]=['',1.3617277,1.2159029,6.41];
ybs[1632]=['9 Aur',1.3458488,0.901047,5];
ybs[1633]=['11 Ori',1.3345477,0.2693711,4.68];
ybs[1634]=['',1.341817,0.6277115,6.52];
ybs[1635]=['',1.3303571,-0.2502687,6.41];
ybs[1636]=['η Aur',1.3443225,0.7201751,3.17];
ybs[1637]=['',1.338946,0.3461952,6.44];
ybs[1638]=['',1.3754675,1.2910505,5.43];
ybs[1639]=['',1.3458145,0.7540363,6.2];
ybs[1640]=['',1.3300333,-0.4251214,5.61];
ybs[1641]=['',1.3353331,-0.0525379,6.05];
ybs[1642]=['',1.3611017,1.1335235,6.41];
ybs[1643]=['',1.3376131,0.0210619,6.17];
ybs[1644]=['η1 Pic',1.32384,-0.8573164,5.38];
ybs[1645]=['',1.3866342,1.3351172,6.37];
ybs[1646]=['',1.3291531,-0.7280612,6.31];
ybs[1647]=['γ1 Cae',1.3317437,-0.6187796,4.55];
ybs[1648]=['γ2 Cae',1.3318549,-0.6226535,6.34];
ybs[1649]=['ε Lep',1.3369855,-0.3899386,3.19];
ybs[1650]=['',1.3359805,-0.4559343,5.73];
ybs[1651]=['104 Tau',1.3472549,0.3259074,5];
ybs[1652]=['66 Eri',1.3433592,-0.080747,5.12];
ybs[1653]=['106 Tau',1.3488856,0.3568545,5.3];
ybs[1654]=['103 Tau',1.3503769,0.4239934,5.5];
ybs[1655]=['105 Tau',1.3494571,0.3793051,5.89];
ybs[1656]=['',1.342382,-0.2285212,6.05];
ybs[1657]=['13 Ori',1.3477133,0.165806,6.17];
ybs[1658]=['η2 Pic',1.3332072,-0.8647782,5.03];
ybs[1659]=['14 Ori',1.3487379,0.1488111,5.34];
ybs[1660]=['',1.3459187,-0.2175089,5.97];
ybs[1661]=['β Eri',1.3480916,-0.0882862,2.79];
ybs[1662]=['',1.3329227,-0.9490726,6.27];
ybs[1663]=['',1.3630719,0.8201041,5.68];
ybs[1664]=['',1.3606649,0.6515048,6.02];
ybs[1665]=['',1.3577158,0.4896944,6.01];
ybs[1666]=['',1.3500798,-0.1507491,5.78];
ybs[1667]=['16 Ori',1.3550942,0.1720299,5.43];
ybs[1668]=['68 Eri',1.3519402,-0.0772939,5.12];
ybs[1669]=['ζ Dor',1.33474,-1.0025761,4.72];
ybs[1670]=['',1.3749909,1.0799211,6.17];
ybs[1671]=['15 Ori',1.3569495,0.2726929,4.82];
ybs[1672]=['β Men',1.3196254,-1.2441298,5.31];
ybs[1673]=['14 Cam',1.3771816,1.0946016,6.5];
ybs[1674]=['λ Eri',1.3536105,-0.1523129,4.27];
ybs[1675]=['',1.3484796,-0.6229167,6.52];
ybs[1676]=['',1.3578741,-0.0093987,6.1];
ybs[1677]=['',1.30481,-1.3660284,6.29];
ybs[1678]=['',1.4007006,1.279152,5.74];
ybs[1679]=['',1.3656662,0.2804992,5.18];
ybs[1680]=['',1.3617961,-0.0388789,6.25];
ybs[1681]=['',1.423825,1.3831817,5.05];
ybs[1682]=['',1.363329,-0.0430176,5.9];
ybs[1683]=['',1.3839322,1.037238,6.15];
ybs[1684]=['μ Aur',1.3743251,0.6721143,4.86];
ybs[1685]=['',1.3650556,0.0094356,6.67];
ybs[1686]=['',1.3653587,0.0185495,5.89];
ybs[1687]=['',1.3812073,0.9291781,6.2];
ybs[1688]=['',1.3632352,-0.2063516,5.68];
ybs[1689]=['',1.3598606,-0.4517431,6.41];
ybs[1690]=['',1.3427947,-1.1060389,5.2];
ybs[1691]=['ι Lep',1.3672409,-0.2067091,4.45];
ybs[1692]=['',1.369667,-0.1052766,5.91];
ybs[1693]=['ρ Ori',1.3721321,0.0503728,4.46];
ybs[1694]=['',1.363002,-0.652216,6.57];
ybs[1695]=['',1.333922,-1.2742387,6.27];
ybs[1696]=['',1.373124,0.0347838,6.09];
ybs[1697]=['μ Lep',1.3698319,-0.2823992,3.31];
ybs[1698]=['',1.3742058,0.0102112,6.32];
ybs[1699]=['',1.3728673,-0.1417705,6.37];
ybs[1700]=['κ Lep',1.371271,-0.2254316,4.36];
ybs[1701]=['14 Aur',1.3826308,0.5709251,5.02];
ybs[1702]=['',1.3924156,0.9356506,6.5];
ybs[1703]=['α Aur',1.3890784,0.8032204,0.08];
ybs[1704]=['',1.3785162,0.0904145,5.5];
ybs[1705]=['',1.3745499,-0.2545033,6.21];
ybs[1706]=['108 Tau',1.3823731,0.3893577,6.27];
ybs[1707]=['',1.386632,0.5992637,5.96];
ybs[1708]=['β Ori',1.3771555,-0.1427203,0.12];
ybs[1709]=['',1.5340361,1.4910351,6.6];
ybs[1710]=['',1.3725873,-0.6248397,6.98];
ybs[1711]=['ξ Men',1.293108,-1.4387948,5.85];
ybs[1712]=['',1.3807643,-0.0241763,6.15];
ybs[1713]=['18 Ori',1.3845712,0.1983553,5.56];
ybs[1714]=['15 Cam',1.4024175,1.0147124,6.13];
ybs[1715]=['',1.4071082,1.093878,5.61];
ybs[1716]=['',1.3756539,-0.6274933,5.76];
ybs[1717]=['',1.3957173,0.7472536,5.48];
ybs[1718]=['',1.3801555,-0.4698314,5.07];
ybs[1719]=['',1.3869042,0.0343907,6.42];
ybs[1720]=['',1.397343,0.7066322,6.18];
ybs[1721]=['16 Aur',1.3947642,0.5828353,4.54];
ybs[1722]=['',1.3718836,-0.9076797,6.05];
ybs[1723]=['17 Aur',1.395382,0.5897377,6.14];
ybs[1724]=['λ Aur',1.3993635,0.7002428,4.71];
ybs[1725]=['',1.3813982,-0.6091699,6.66];
ybs[1726]=['',1.3866776,-0.2987735,6.56];
ybs[1727]=['',1.3983707,0.5894016,5.41];
ybs[1728]=['',1.4035749,0.7757438,6.62];
ybs[1729]=['18 Aur',1.4001002,0.5935382,6.49];
ybs[1730]=['τ Ori',1.3905946,-0.1190611,3.6];
ybs[1731]=['',1.4064625,0.8200401,6.54];
ybs[1732]=['',1.3906121,-0.235567,5.5];
ybs[1733]=['',1.4042509,0.7174578,5.52];
ybs[1734]=['109 Tau',1.3990136,0.3860348,4.94];
ybs[1735]=['19 Aur',1.402812,0.5930524,5.03];
ybs[1736]=['',1.3987906,0.3517977,6.08];
ybs[1737]=['',1.3795439,-0.9103287,6.49];
ybs[1738]=['ο Col',1.3888228,-0.6086381,4.83];
ybs[1739]=['θ Dor',1.3689649,-1.1721661,4.83];
ybs[1740]=['',1.4526255,1.3612384,6.56];
ybs[1741]=['21 Ori',1.3978443,0.0456876,5.34];
ybs[1742]=['',1.3955295,-0.3160423,5.96];
ybs[1743]=['',1.399439,-0.0242697,6.34];
ybs[1744]=['ρ Aur',1.411108,0.7299803,5.23];
ybs[1745]=['',1.4067522,0.4883093,6.33];
ybs[1746]=['16 Cam',1.4198043,1.004678,5.28];
ybs[1747]=['',1.4078059,0.5164554,5.76];
ybs[1748]=['',1.3974838,-0.3228533,6.36];
ybs[1749]=['',1.3975424,-0.3226741,6.54];
ybs[1750]=['',1.4061884,0.3461869,6.18];
ybs[1751]=['λ Lep',1.3989371,-0.2295974,4.29];
ybs[1752]=['ν Lep',1.4007599,-0.2145721,5.3];
ybs[1753]=['',1.3975414,-0.4772961,5.99];
ybs[1754]=['',1.403013,-0.0932957,6.39];
ybs[1755]=['',1.4155589,0.7164444,5.54];
ybs[1756]=['',1.4072136,0.0703832,6.57];
ybs[1757]=['',1.4024144,-0.3703273,4.71];
ybs[1758]=['',1.4091491,0.1474646,5.8];
ybs[1759]=['',1.4079557,-0.0069075,5.68];
ybs[1760]=['22 Ori',1.4089678,-0.0063182,4.73];
ybs[1761]=['',1.4013022,-0.6052373,6.34];
ybs[1762]=['ζ Pic',1.3959189,-0.8828603,5.45];
ybs[1763]=['22 Aur',1.4172448,0.5053812,6.46];
ybs[1764]=['',1.4088396,-0.239732,6.56];
ybs[1765]=['23 Ori',1.413792,0.0622096,5];
ybs[1766]=['',1.4080234,-0.4320128,5.06];
ybs[1767]=['',1.4054149,-0.5990744,6.09];
ybs[1768]=['σ Aur',1.4232514,0.6528296,4.99];
ybs[1769]=['110 Tau',1.4177751,0.2917997,6.08];
ybs[1770]=['',1.4228443,0.5452888,6.28];
ybs[1771]=['',1.4228547,0.543878,5.94];
ybs[1772]=['',1.4168483,0.093236,6.35];
ybs[1773]=['',1.4154114,-0.1465406,5.9];
ybs[1774]=['',1.4255635,0.6086629,6.55];
ybs[1775]=['111 Tau',1.4212724,0.3037283,4.99];
ybs[1776]=['',1.4174556,-0.0024483,5.7];
ybs[1777]=['',1.4180907,-0.0147979,6.11];
ybs[1778]=['8 Lep',1.4160462,-0.2427341,5.25];
ybs[1779]=['29 Ori',1.4182199,-0.1359389,4.14];
ybs[1780]=['',1.4141722,-0.4657594,6.49];
ybs[1781]=['',1.4214701,0.0413945,6.32];
ybs[1782]=['27 Ori',1.4208169,-0.0152256,5.08];
ybs[1783]=['η Ori',1.4207381,-0.0415024,3.36];
ybs[1784]=['ψ1 Ori',1.4220764,0.0325551,4.95];
ybs[1785]=['γ Ori',1.4239272,0.1111492,1.64];
ybs[1786]=['β Tau',1.4299282,0.4996088,1.65];
ybs[1787]=['',1.4201545,-0.2959562,5.65];
ybs[1788]=['',1.41434,-0.6921782,5.71];
ybs[1789]=['',1.4329658,0.6191528,6.15];
ybs[1790]=['',1.4325132,0.6005562,5.94];
ybs[1791]=['',1.432626,0.5808531,6.15];
ybs[1792]=['',1.4155746,-0.6513062,6.82];
ybs[1793]=['113 Tau',1.428539,0.2917909,6.25];
ybs[1794]=['',1.4228345,-0.1799507,5.61];
ybs[1795]=['',1.4253604,-0.0091752,6.57];
ybs[1796]=['κ Pic',1.4084158,-0.9793755,6.11];
ybs[1797]=['17 Cam',1.449941,1.1010044,5.42];
ybs[1798]=['',1.4265497,0.00941,6.16];
ybs[1799]=['',1.4337029,0.5275453,5.74];
ybs[1800]=['φ Aur',1.4361548,0.6020173,5.07];
ybs[1801]=['',1.427441,-0.0959955,6.23];
ybs[1802]=['',1.4305578,0.1202009,6.42];
ybs[1803]=['115 Tau',1.4332751,0.3138058,5.42];
ybs[1804]=['',1.4334328,0.266604,6.16];
ybs[1805]=['114 Tau',1.4354816,0.3831732,4.88];
ybs[1806]=['ψ2 Ori',1.4312431,0.0543373,4.59];
ybs[1807]=['',1.4266867,-0.3434337,5.89];
ybs[1808]=['',1.4206889,-0.7715563,6.08];
ybs[1809]=['116 Tau',1.4357707,0.2773569,5.5];
ybs[1810]=['',1.3539082,-1.4227101,6.51];
ybs[1811]=['117 Tau',1.4369904,0.3011732,5.77];
ybs[1812]=['',1.4317316,-0.2074007,6.35];
ybs[1813]=['θ Pic',1.4193117,-0.9127609,6.27];
ybs[1814]=['',1.4392601,0.2390344,6.35];
ybs[1815]=['',1.4363674,0.0229587,6.41];
ybs[1816]=['118 Tau',1.442783,0.4392457,5.47];
ybs[1817]=['',1.4447289,0.5096803,6.24];
ybs[1818]=['',1.4336457,-0.3727704,6.07];
ybs[1819]=['',1.45038,0.7239176,6];
ybs[1820]=['',1.4500178,0.6953627,6.37];
ybs[1821]=['',1.4402001,-0.0574368,6.39];
ybs[1822]=['',1.4303468,-0.7142913,5.87];
ybs[1823]=['18 Cam',1.4594865,0.9989488,6.48];
ybs[1824]=['β Lep',1.4364501,-0.3620233,2.84];
ybs[1825]=['',1.4421511,-0.059865,5.79];
ybs[1826]=['',1.4489818,0.3923167,6.29];
ybs[1827]=['',1.4474225,0.2683625,5.94];
ybs[1828]=['',1.4446183,0.0315075,5.78];
ybs[1829]=['31 Ori',1.4437238,-0.0187804,4.71];
ybs[1830]=['',1.4356691,-0.6494976,5.57];
ybs[1831]=['λ Dor',1.4252926,-1.0278985,5.14];
ybs[1832]=['',1.446528,0.0736532,6.21];
ybs[1833]=['',1.4397957,-0.5253453,6.75];
ybs[1834]=['32 Ori',1.4485798,0.1040854,4.2];
ybs[1835]=['',1.446151,-0.1294834,6.33];
ybs[1836]=['33 Ori',1.4504712,0.0577281,5.46];
ybs[1837]=['χ Aur',1.4582004,0.562108,4.76];
ybs[1838]=['',1.4955361,1.3099417,6.17];
ybs[1839]=['119 Tau',1.4553165,0.3247925,4.38];
ybs[1840]=['',1.4620695,0.7351836,6.55];
ybs[1841]=['',1.4553541,0.2979772,5.46];
ybs[1842]=['',1.4505566,-0.1168151,6.22];
ybs[1843]=['10 Lep',1.4490178,-0.3638686,5.55];
ybs[1844]=['',1.4614187,0.572733,6.48];
ybs[1845]=['δ Ori',1.4536808,-0.0047036,6.85];
ybs[1846]=['δ Ori',1.4536731,-0.0049605,2.23];
ybs[1847]=['',1.4817042,1.1642722,6.26];
ybs[1848]=['',1.462288,0.606324,6.27];
ybs[1849]=['υ Ori',1.4530717,-0.1271714,4.62];
ybs[1850]=['',1.4432885,-0.8213813,5.46];
ybs[1851]=['19 Cam',1.4810642,1.1199174,6.15];
ybs[1852]=['120 Tau',1.4610452,0.3238346,5.69];
ybs[1853]=['',1.4263138,-1.1973792,6.03];
ybs[1854]=['',1.4616514,0.3575861,6.18];
ybs[1855]=['',1.4565984,-0.0275302,5.35];
ybs[1856]=['ε Col',1.4486699,-0.6188079,3.87];
ybs[1857]=['',1.458477,-0.0297402,6.46];
ybs[1858]=['35 Ori',1.4625215,0.2499213,5.64];
ybs[1859]=['α Lep',1.4561425,-0.3108017,2.58];
ybs[1860]=['',1.4767111,0.9501775,5.73];
ybs[1861]=['',1.4377022,-1.0873011,6.59];
ybs[1862]=['',1.4602583,-0.0199313,5.34];
ybs[1863]=['',1.4746571,0.8330066,6.11];
ybs[1864]=['',1.4496764,-0.8012801,5.86];
ybs[1865]=['',1.4622608,0.0248127,6.59];
ybs[1866]=['38 Ori',1.4637395,0.0659848,5.36];
ybs[1867]=['',1.4626262,-0.0178324,6.22];
ybs[1868]=['',1.462617,-0.0254245,5.93];
ybs[1869]=['121 Tau',1.4696859,0.4197947,5.38];
ybs[1870]=['φ1 Ori',1.4663232,0.1658559,4.41];
ybs[1871]=['',1.4556581,-0.6719248,5.48];
ybs[1872]=['',1.4719214,0.4830194,6.27];
ybs[1873]=['λ Ori',1.4677296,0.1736147,3.54];
ybs[1874]=['λ Ori',1.4677442,0.1736293,5.61];
ybs[1875]=['',1.4570283,-0.6130519,5.78];
ybs[1876]=['',1.4416554,-1.115468,6.19];
ybs[1877]=['',1.4681124,0.1789517,5.6];
ybs[1878]=['',1.4767569,0.7015249,6.09];
ybs[1879]=['',1.6075869,1.481411,6.11];
ybs[1880]=['',1.4665819,-0.1046471,5.67];
ybs[1881]=['',1.466713,-0.1045213,4.78];
ybs[1882]=['',1.4605608,-0.5207168,6.53];
ybs[1883]=['',1.474363,0.4529458,6.49];
ybs[1884]=['',1.4681597,-0.0781942,6.56];
ybs[1885]=['',1.4682132,-0.0770065,6.24];
ybs[1886]=['42 Ori',1.4682483,-0.0842157,4.59];
ybs[1887]=['θ1 Ori',1.4676965,-0.0937945,6.73];
ybs[1888]=['θ1 Ori',1.4677111,-0.0937606,7.96];
ybs[1889]=['θ1 Ori',1.4677401,-0.0938382,5.13];
ybs[1890]=['θ1 Ori',1.4677982,-0.0938045,6.7];
ybs[1891]=['θ2 Ori',1.4682044,-0.0942999,5.08];
ybs[1892]=['',1.4688409,-0.075951,6.38];
ybs[1893]=['ι Ori',1.4684109,-0.1029203,2.77];
ybs[1894]=['',1.4692251,-0.0565447,6.4];
ybs[1895]=['45 Ori',1.469433,-0.0845238,5.26];
ybs[1896]=['',1.477217,0.4701264,5.83];
ybs[1897]=['ε Ori',1.4719937,-0.0207568,1.7];
ybs[1898]=['',1.4802393,0.5859228,6.33];
ybs[1899]=['122 Tau',1.4764151,0.2976212,5.54];
ybs[1900]=['',1.4719843,-0.0983563,6.54];
ybs[1901]=['φ2 Ori',1.4754133,0.1623649,4.09];
ybs[1902]=['',1.4762162,0.1928095,5.94];
ybs[1903]=['',1.4664492,-0.5771233,5.78];
ybs[1904]=['ζ Tau',1.4791279,0.3692127,3];
ybs[1905]=['',1.4734735,-0.1056365,5.72];
ybs[1906]=['',1.4581571,-0.9579767,6.43];
ybs[1907]=['',1.4772182,0.156451,6.12];
ybs[1908]=['26 Aur',1.4838881,0.5323911,5.4];
ybs[1909]=['',1.4706665,-0.5008225,6.26];
ybs[1910]=['',1.5040382,1.1467991,5.6];
ybs[1911]=['',1.4534939,-1.1207249,5.34];
ybs[1912]=['',1.4772381,-0.103434,6.05];
ybs[1913]=['',1.4756582,-0.2053095,6.11];
ybs[1914]=['',1.4802033,0.1318257,5.88];
ybs[1915]=['',1.4851008,0.4647617,6.37];
ybs[1916]=['β Dor',1.4565763,-1.0904013,3.76];
ybs[1917]=['',1.479172,-0.0838081,6.19];
ybs[1918]=['',1.4867546,0.5100931,5.96];
ybs[1919]=['',1.4973453,0.9335906,6.23];
ybs[1920]=['ν1 Col',1.4755211,-0.4862353,6.16];
ybs[1921]=['',1.4689474,-0.8255579,6.11];
ybs[1922]=['125 Tau',1.4884705,0.4521733,5.18];
ybs[1923]=['',1.4870315,0.3800214,6.34];
ybs[1924]=['',1.4633216,-1.0272586,6.75];
ybs[1925]=['σ Ori',1.4829943,-0.0451815,3.81];
ybs[1926]=['',1.4831617,-0.04508,6.65];
ybs[1927]=['',1.4823331,-0.1145376,5.96];
ybs[1928]=['ω Ori',1.4851627,0.0721244,4.57];
ybs[1929]=['ν2 Col',1.4775243,-0.5005175,5.31];
ybs[1930]=['',1.4626222,-1.0674825,6.32];
ybs[1931]=['49 Ori',1.4834212,-0.1256956,4.8];
ybs[1932]=['',1.4925005,0.5474792,6.04];
ybs[1933]=['',1.492981,0.5573005,6.11];
ybs[1934]=['',1.4863318,-0.0620263,6];
ybs[1935]=['24 Cam',1.5050737,0.9876891,6.05];
ybs[1936]=['',1.4860658,-0.1692232,6.5];
ybs[1937]=['23 Cam',1.5106553,1.0731115,6.15];
ybs[1938]=['',1.4846826,-0.3113388,6.38];
ybs[1939]=['',1.4956849,0.5148249,6.43];
ybs[1940]=['126 Tau',1.4948593,0.2887428,4.86];
ybs[1941]=['',1.4811329,-0.7102807,5.82];
ybs[1942]=['ζ Ori',1.4917967,-0.03373,2.05];
ybs[1943]=['ζ Ori',1.491804,-0.03373,4.21];
ybs[1944]=['',1.4911668,-0.0491263,6.22];
ybs[1945]=['',1.4978358,0.4072882,6.59];
ybs[1946]=['',1.4921986,-0.0195257,4.95];
ybs[1947]=['γ Men',1.4442173,-1.3321319,5.19];
ybs[1948]=['',1.4984822,0.3956609,6.36];
ybs[1949]=['',1.4933452,0.00607,5.93];
ybs[1950]=['α Col',1.4855543,-0.5945164,2.64];
ybs[1951]=['',1.4915067,-0.1815009,6.52];
ybs[1952]=['',1.4864256,-0.5692982,5.45];
ybs[1953]=['',1.4957455,-0.0503774,6.42];
ybs[1954]=['',1.4700708,-1.1614761,6.31];
ybs[1955]=['',1.5040048,0.4051415,6.21];
ybs[1956]=['',1.4952881,-0.2917463,6.21];
ybs[1957]=['51 Ori',1.4994244,0.0259003,4.91];
ybs[1958]=['',1.4581856,-1.2867869,5.78];
ybs[1959]=['',1.4976397,-0.3057966,6.15];
ybs[1960]=['',1.4934516,-0.5827769,6.34];
ybs[1961]=['',1.500941,-0.1184567,6.02];
ybs[1962]=['12 Lep',1.4974019,-0.3903284,5.87];
ybs[1963]=['26 Cam',1.5202033,0.9795208,5.94];
ybs[1964]=['',1.5022666,-0.027998,6.31];
ybs[1965]=['ο Aur',1.5168908,0.8697606,5.47];
ybs[1966]=['',1.4968461,-0.5327804,6.19];
ybs[1967]=['',1.4968789,-0.6049016,5.29];
ybs[1968]=['',1.5158427,0.707112,6.58];
ybs[1969]=['',1.5024904,-0.3237356,5.73];
ybs[1970]=['',1.5324399,1.0963014,6.13];
ybs[1971]=['',1.514068,0.3613262,6.95];
ybs[1972]=['',1.5106695,0.0700908,6.09];
ybs[1973]=['',1.5221765,0.7423443,6.29];
ybs[1974]=['',1.5072809,-0.3511284,6.34];
ybs[1975]=['',1.5020562,-0.6876272,6.25];
ybs[1976]=['',1.5070421,-0.3911881,6.15];
ybs[1977]=['γ Lep',1.5071354,-0.3916537,3.6];
ybs[1978]=['',1.5023955,-0.7997852,6.39];
ybs[1979]=['129 Tau',1.5186692,0.2762748,6];
ybs[1980]=['',1.5147824,-0.0743688,6.34];
ybs[1981]=['',1.5188966,0.1663134,5.79];
ybs[1982]=['',1.5173315,0.0205088,5.95];
ybs[1983]=['131 Tau',1.5206291,0.2529848,5.72];
ybs[1984]=['130 Tau',1.5217083,0.3095458,5.49];
ybs[1985]=['ι Men',1.4583953,-1.3754443,6.05];
ybs[1986]=['29 Cam',1.5380215,0.9935028,6.54];
ybs[1987]=['133 Tau',1.5227651,0.242707,5.29];
ybs[1988]=['',1.5504982,1.1951071,6.2];
ybs[1989]=['τ Aur',1.5303755,0.6839352,4.52];
ybs[1990]=['μ Col',1.5133544,-0.5637234,5.17];
ybs[1991]=['',1.5259298,0.364345,6.07];
ybs[1992]=['ζ Lep',1.518331,-0.2585722,3.55];
ybs[1993]=['52 Ori',1.5237299,0.1127551,5.27];
ybs[1994]=['',1.5190219,-0.2832848,6.17];
ybs[1995]=['',1.5206407,-0.1837218,6.03];
ybs[1996]=['132 Tau',1.5289028,0.4288819,4.86];
ybs[1997]=['',1.5390448,0.8991794,6.29];
ybs[1998]=['κ Ori',1.5220345,-0.1686568,2.06];
ybs[1999]=['',1.5182546,-0.4997288,6.22];
ybs[2000]=['30 Cam',1.5458189,1.0291829,6.14];
ybs[2001]=['',1.5258499,-0.0713628,5.97];
ybs[2002]=['',1.5144266,-0.8131489,5.31];
ybs[2003]=['',1.5188905,-0.6225242,6.32];
ybs[2004]=['134 Tau',1.5307144,0.2208971,4.91];
ybs[2005]=['υ Aur',1.5383976,0.6511829,4.74];
ybs[2006]=['ν Aur',1.5404751,0.6833459,3.97];
ybs[2007]=['',1.5375802,0.4882092,5.56];
ybs[2008]=['',1.5327641,0.1723724,5.8];
ybs[2009]=['δ Dor',1.5045508,-1.147157,4.35];
ybs[2010]=['135 Tau',1.5348537,0.2497637,5.52];
ybs[2011]=['',1.5214553,-0.7094087,6.61];
ybs[2012]=['',1.5397992,0.5607618,6.25];
ybs[2013]=['',1.5333166,0.0772893,5.97];
ybs[2014]=['β Pic',1.5176465,-0.8911579,3.85];
ybs[2015]=['',1.5299217,-0.2526925,5.49];
ybs[2016]=['π Men',1.4632415,-1.4042258,5.65];
ybs[2017]=['',1.5170059,-0.9486558,6.18];
ybs[2018]=['',1.5344538,0.0354182,5.98];
ybs[2019]=['',1.5456004,0.690767,6.45];
ybs[2020]=['',1.5307963,-0.4008393,5.87];
ybs[2021]=['31 Cam',1.5576207,1.045287,5.2];
ybs[2022]=['',1.5453072,0.5920348,5.98];
ybs[2023]=['ξ Aur',1.5565478,0.9723098,4.99];
ybs[2024]=['',1.5434129,0.3468292,6.06];
ybs[2025]=['55 Ori',1.5378698,-0.1311376,5.35];
ybs[2026]=['',1.5281369,-0.7831249,6.38];
ybs[2027]=['137 Tau',1.543095,0.2474089,5.59];
ybs[2028]=['136 Tau',1.5478534,0.4819811,4.58];
ybs[2029]=['δ Lep',1.5371284,-0.364332,3.81];
ybs[2030]=['',1.537714,-0.4000593,6.17];
ybs[2031]=['56 Ori',1.5429122,0.0324424,4.78];
ybs[2032]=['',1.5474823,0.3543446,6.71];
ybs[2033]=['',1.5411265,-0.157732,5.97];
ybs[2034]=['β Col',1.5348099,-0.6241928,3.12];
ybs[2035]=['',1.570216,1.153607,6.25];
ybs[2036]=['γ Pic',1.5281876,-0.9801982,4.51];
ybs[2037]=['',1.5396516,-0.513903,6.45];
ybs[2038]=['',1.5314003,-0.920883,6.35];
ybs[2039]=['',1.5622004,0.9041757,6.49];
ybs[2040]=['',1.5552889,0.5533394,5.9];
ybs[2041]=['χ1 Ori',1.552121,0.3539321,4.41];
ybs[2042]=['',1.5510312,0.1848262,6.12];
ybs[2043]=['',1.5332567,-0.9093872,5.17];
ybs[2044]=['',1.5524449,0.2053405,6.59];
ybs[2045]=['',1.5509131,0.0563408,6.31];
ybs[2046]=['57 Ori',1.5545563,0.3447396,5.92];
ybs[2047]=['',1.5416512,-0.6567194,5.63];
ybs[2048]=['',1.565612,0.8557448,6.47];
ybs[2049]=['',1.5426506,-0.6723375,6.7];
ybs[2050]=['λ Col',1.5443236,-0.5898837,4.87];
ybs[2051]=['',1.5528921,0.0169453,6];
ybs[2052]=['',1.5520163,-0.0708818,6.57];
ybs[2053]=['',1.4222567,-1.479474,6.2];
ybs[2054]=['',1.5486792,-0.3427006,6.69];
ybs[2055]=['α Ori',1.5550446,0.1293158,0.5];
ybs[2056]=['λ Men',1.5155944,-1.2687747,6.53];
ybs[2057]=['',1.5583999,0.3521537,5.4];
ybs[2058]=['ε Dor',1.5266092,-1.1675482,5.11];
ybs[2059]=['',1.5523633,-0.2054478,5.66];
ybs[2060]=['',1.5620417,0.5051631,6.32];
ybs[2061]=['',1.5591714,0.2430734,6.6];
ybs[2062]=['',1.5494547,-0.5086739,6.36];
ybs[2063]=['',1.5449179,-0.7490598,6.55];
ybs[2064]=['',1.5560385,-0.0805336,5.87];
ybs[2065]=['',1.5564028,-0.0835354,6.28];
ybs[2066]=['',1.5390131,-0.997491,5.94];
ybs[2067]=['',1.5337451,-1.1175157,6.36];
ybs[2068]=['',1.5634392,0.4232602,6.02];
ybs[2069]=['',1.5607773,0.1660038,5.99];
ybs[2070]=['',1.5624202,0.2011056,5.87];
ybs[2071]=['δ Aur',1.5768035,0.9474428,3.72];
ybs[2072]=['',1.6068225,1.3191577,6.4];
ybs[2073]=['',1.5779479,0.965524,6.44];
ybs[2074]=['',1.578033,0.9520216,6.14];
ybs[2075]=['',1.5756485,0.8713436,5.89];
ybs[2076]=['',1.5516347,-0.6973542,5.57];
ybs[2077]=['',1.547788,-0.8789288,6.52];
ybs[2078]=['139 Tau',1.5681438,0.4529933,4.82];
ybs[2079]=['η Lep',1.5595901,-0.247245,3.71];
ybs[2080]=['',1.558494,-0.3986066,5.96];
ybs[2081]=['ξ Col',1.5545315,-0.6478416,4.97];
ybs[2082]=['β Aur',1.5759523,0.7844785,1.9];
ybs[2083]=['',1.5500665,-0.8661056,6.1];
ybs[2084]=['',1.559946,-0.4051597,6.36];
ybs[2085]=['π Aur',1.5778034,0.8017436,4.26];
ybs[2086]=['σ Col',1.5585612,-0.5476972,5.5];
ybs[2087]=['',1.5667264,0.021385,6.22];
ybs[2088]=['',1.5504364,-0.918612,5.29];
ybs[2089]=['θ Aur',1.5762738,0.6494759,2.62];
ybs[2090]=['',1.579365,0.7782654,6.22];
ybs[2091]=['',1.5679143,-0.0173397,6.22];
ybs[2092]=['',1.5605734,-0.558062,6.44];
ybs[2093]=['',1.5714673,0.223557,5.7];
ybs[2094]=['59 Ori',1.568946,0.0320703,5.9];
ybs[2095]=['36 Aur',1.5825092,0.8360291,5.73];
ybs[2096]=['',1.5457859,-1.1010727,4.65];
ybs[2097]=['60 Ori',1.5707295,0.0096583,5.22];
ybs[2098]=['',1.5459354,-1.1253723,6.63];
ybs[2099]=['',1.5858369,0.8544789,5.96];
ybs[2100]=['γ Col',1.5635328,-0.6157906,4.36];
ybs[2101]=['1 Mon',1.5711762,-0.1637464,6.12];
ybs[2102]=['2 Mon',1.5714093,-0.1668206,5.03];
ybs[2103]=['',1.5741512,-0.0252122,6.63];
ybs[2104]=['',1.5822485,0.54164,5.98];
ybs[2105]=['',1.581362,0.4812097,6.05];
ybs[2106]=['',1.5906926,0.8709811,6.05];
ybs[2107]=['',1.5759579,-0.0536602,4.53];
ybs[2108]=['',1.560839,-0.9324375,6.45];
ybs[2109]=['',1.5905182,0.7570642,6.42];
ybs[2110]=['',1.5841177,0.3909457,6.37];
ybs[2111]=['',1.567722,-0.7685361,5.81];
ybs[2112]=['',1.5766145,-0.2251504,6.22];
ybs[2113]=['38 Aur',1.5922591,0.7489107,6.1];
ybs[2114]=['η Col',1.5700792,-0.7472627,3.96];
ybs[2115]=['',1.6019424,1.0365459,6.34];
ybs[2116]=['',1.5899616,0.5695678,6.24];
ybs[2117]=['',1.5981524,0.9000732,6.45];
ybs[2118]=['μ Ori',1.5865989,0.168352,4.12];
ybs[2119]=['κ Men',1.5217858,-1.3850176,5.47];
ybs[2120]=['',1.6092319,1.1074007,6.39];
ybs[2121]=['',1.5858629,0.0295462,6.59];
ybs[2122]=['3 Mon',1.5834475,-0.1849936,4.95];
ybs[2123]=['',1.5801051,-0.4436399,6.05];
ybs[2124]=['64 Ori',1.5916855,0.3436257,5.14];
ybs[2125]=['',1.5799143,-0.5918862,5.55];
ybs[2126]=['39 Aur',1.5999286,0.7501156,5.87];
ybs[2127]=['',1.5911661,0.2038305,6.08];
ybs[2128]=['1 Gem',1.5947446,0.405976,4.16];
ybs[2129]=['χ2 Ori',1.5937338,0.3514364,4.63];
ybs[2130]=['',1.5864475,-0.2530536,6.2];
ybs[2131]=['',1.599539,0.6625442,6.34];
ybs[2132]=['',1.5766853,-0.8939049,5.67];
ybs[2133]=['',1.6015691,0.5863559,6.23];
ybs[2134]=['',1.5889655,-0.4587855,5.04];
ybs[2135]=['',1.6041783,0.6175626,6.12];
ybs[2136]=['',1.594006,-0.1171425,5.21];
ybs[2137]=['40 Aur',1.6062975,0.671581,5.36];
ybs[2138]=['63 Ori',1.597721,0.0945437,5.67];
ybs[2139]=['66 Ori',1.5976872,0.0725284,5.63];
ybs[2140]=['',1.6048811,0.5150226,6.08];
ybs[2141]=['',1.6102846,0.7304143,6.12];
ybs[2142]=['17 Lep',1.5969324,-0.2877601,4.93];
ybs[2143]=['',1.5933826,-0.5615615,5.65];
ybs[2144]=['',1.5992131,-0.1788273,5.87];
ybs[2145]=['',1.5814223,-1.0489115,6.45];
ybs[2146]=['37 Cam',1.6230206,1.0285196,5.36];
ybs[2147]=['',1.6143097,0.7164671,6.36];
ybs[2148]=['',1.604662,-0.0732659,5.38];
ybs[2149]=['θ Lep',1.6021013,-0.2607334,4.67];
ybs[2150]=['',1.599986,-0.4223515,6.95];
ybs[2151]=['',1.593158,-0.7860839,6.35];
ybs[2152]=['',1.5940129,-0.7868226,5.93];
ybs[2153]=['ν Ori',1.6094416,0.2576777,4.42];
ybs[2154]=['',1.5980556,-0.6198852,5.8];
ybs[2155]=['',1.6053505,-0.1950868,6.66];
ybs[2156]=['',1.594225,-0.8458107,6.58];
ybs[2157]=['',1.6034025,-0.4034221,5.47];
ybs[2158]=['',1.6011584,-0.5194479,5.81];
ybs[2159]=['36 Cam',1.6367403,1.146868,5.32];
ybs[2160]=['',1.6053222,-0.3807757,5.78];
ybs[2161]=['',1.6145001,0.1512303,6.55];
ybs[2162]=['19 Lep',1.608637,-0.3345849,5.31];
ybs[2163]=['',1.6183449,0.3871912,5.93];
ybs[2164]=['',1.6051491,-0.5989276,5.83];
ybs[2165]=['π1 Col',1.6030012,-0.7383171,6.12];
ybs[2166]=['',1.6300295,0.9187466,6.3];
ybs[2167]=['3 Gem',1.6192298,0.4033045,5.75];
ybs[2168]=['',1.6150405,0.0435324,5.73];
ybs[2169]=['41 Aur',1.6289647,0.8500843,6.82];
ybs[2170]=['41 Aur',1.6289718,0.8500504,6.09];
ybs[2171]=['θ Col',1.6070075,-0.650264,5.02];
ybs[2172]=['',1.6043173,-0.7870634,6.51];
ybs[2173]=['',1.617512,-0.0997745,6.17];
ybs[2174]=['',1.6140435,-0.3915237,5.5];
ybs[2175]=['π2 Col',1.6082403,-0.7358027,5.5];
ybs[2176]=['',1.6158591,-0.3164541,6.35];
ybs[2177]=['',1.6170239,-0.2546379,5.56];
ybs[2178]=['',1.6246703,0.316307,6.33];
ybs[2179]=['5 Gem',1.6271646,0.4260977,5.8];
ybs[2180]=['',1.6176643,-0.3975821,5.71];
ybs[2181]=['',1.6111146,-0.7742452,6.27];
ybs[2182]=['',1.6385561,0.892989,6.04];
ybs[2183]=['',1.6310504,0.5704817,5.78];
ybs[2184]=['',1.6284385,0.3815598,6.56];
ybs[2185]=['',1.6263831,0.2379233,6.04];
ybs[2186]=['',1.6174318,-0.4661151,6.27];
ybs[2187]=['68 Ori',1.6290678,0.3452894,5.75];
ybs[2188]=['η1 Dor',1.5977734,-1.1526694,5.71];
ybs[2189]=['',1.6236604,-0.1179927,6.15];
ybs[2190]=['',1.6024815,-1.0848728,5.05];
ybs[2191]=['6 Gem',1.6304992,0.3997018,6.39];
ybs[2192]=['69 Ori',1.6290591,0.2814102,4.95];
ybs[2193]=['ξ Ori',1.6284787,0.247872,4.48];
ybs[2194]=['',1.6208572,-0.4740347,5.72];
ybs[2195]=['40 Cam',1.6480063,1.047024,5.35];
ybs[2196]=['',1.6268242,-0.0815464,6.18];
ybs[2197]=['',1.6183545,-0.7044039,5.58];
ybs[2198]=['',1.614229,-0.8651259,6.49];
ybs[2199]=['',1.6273408,-0.1144422,5.05];
ybs[2200]=['',1.6237184,-0.4623133,6.09];
ybs[2201]=['',1.6357177,0.3258969,6.58];
ybs[2202]=['',1.6338769,0.185353,6.45];
ybs[2203]=['',1.66383,1.2096661,4.8];
ybs[2204]=['',1.6313212,-0.0438377,6.62];
ybs[2205]=['',1.6201574,-0.7904238,6.31];
ybs[2206]=['δ Pic',1.6176473,-0.9594835,4.81];
ybs[2207]=['',1.6308556,-0.3101503,6.52];
ybs[2208]=['',1.6396991,0.3123811,5.88];
ybs[2209]=['1 Lyn',1.6580043,1.0734635,4.98];
ybs[2210]=['η Gem',1.6416433,0.3926671,3.28];
ybs[2211]=['',1.6457068,0.6307559,6.92];
ybs[2212]=['',1.6363572,-0.0654375,5.83];
ybs[2213]=['κ Aur',1.6441623,0.5146848,4.35];
ybs[2214]=['71 Ori',1.6413737,0.3341942,5.2];
ybs[2215]=['ν Dor',1.6083031,-1.2016255,5.06];
ybs[2216]=['',1.6424336,0.2415971,5.91];
ybs[2217]=['72 Ori',1.6437343,0.2815964,5.3];
ybs[2218]=['',1.639409,-0.0798771,5.83];
ybs[2219]=['',1.6348656,-0.416605,6.39];
ybs[2220]=['',1.6337433,-0.5131923,6.54];
ybs[2221]=['γ Mon',1.6404052,-0.1096614,3.98];
ybs[2222]=['42 Aur',1.654818,0.8100791,6.52];
ybs[2223]=['73 Ori',1.6450357,0.2189021,5.33];
ybs[2224]=['8 Gem',1.6479921,0.4181933,6.08];
ybs[2225]=['',1.6444256,0.1057185,6.07];
ybs[2226]=['',1.644859,0.074607,6.64];
ybs[2227]=['',1.6437526,-0.0090939,5.65];
ybs[2228]=['',1.64325,-0.0859408,5.99];
ybs[2229]=['',1.6480459,0.2997093,6.39];
ybs[2230]=['',1.6452493,0.0202487,6.37];
ybs[2231]=['',1.6428217,-0.1578475,6.1];
ybs[2232]=['2 Lyn',1.6650756,1.0297369,4.48];
ybs[2233]=['43 Aur',1.6578811,0.8089621,6.38];
ybs[2234]=['9 Gem',1.6508685,0.4141874,6.25];
ybs[2235]=['74 Ori',1.6480497,0.214028,5.04];
ybs[2236]=['',1.6410829,-0.3539659,5.91];
ybs[2237]=['',1.6418216,-0.3226292,5.99];
ybs[2238]=['',1.6440106,-0.2395851,5.01];
ybs[2239]=['η2 Dor',1.6201033,-1.1448593,5.01];
ybs[2240]=['',1.6472239,0.018693,6.63];
ybs[2241]=['75 Ori',1.6508664,0.1733604,5.39];
ybs[2242]=['',1.6501577,0.1229315,6.57];
ybs[2243]=['',1.6455503,-0.2901934,5.92];
ybs[2244]=['',1.6529718,0.2451908,6.59];
ybs[2245]=['',1.6513697,0.0888463,5.71];
ybs[2246]=['',1.6441891,-0.5200607,6.67];
ybs[2247]=['',1.6553337,0.2508483,6.16];
ybs[2248]=['',1.6493463,-0.3966186,6.07];
ybs[2249]=['6 Mon',1.6521513,-0.187364,6.75];
ybs[2250]=['κ Col',1.6465134,-0.61348,4.37];
ybs[2251]=['4 Lyn',1.675756,1.0360212,5.94];
ybs[2252]=['',1.6595414,0.3021908,6.32];
ybs[2253]=['',1.6576598,0.15772,6.24];
ybs[2254]=['',1.6523782,-0.2936651,5.14];
ybs[2255]=['α Men',1.6124937,-1.3047813,5.09];
ybs[2256]=['',1.6464476,-0.6854555,6];
ybs[2257]=['',1.6484044,-0.6588095,5.53];
ybs[2258]=['45 Aur',1.6737447,0.9327013,5.36];
ybs[2259]=['',1.6490374,-0.6503557,5.87];
ybs[2260]=['',1.6545723,-0.3486672,5.52];
ybs[2261]=['',1.6576872,-0.1640759,5.36];
ybs[2262]=['',1.6573318,-0.262415,6.06];
ybs[2263]=['',1.6639689,0.2555128,5.69];
ybs[2264]=['',1.6589766,-0.1500482,6.22];
ybs[2265]=['',1.6578252,-0.365415,5.81];
ybs[2266]=['',1.6695708,0.5153809,6.43];
ybs[2267]=['7 Mon',1.6615443,-0.1367311,5.27];
ybs[2268]=['',1.643326,-1.0336294,6.43];
ybs[2269]=['',1.6629532,-0.0515861,4.9];
ybs[2270]=['',1.6673507,0.2049829,6.54];
ybs[2271]=['',1.6700382,0.3098231,6.35];
ybs[2272]=['',1.6508995,-0.920538,6.41];
ybs[2273]=['',1.6602215,-0.6005265,5.78];
ybs[2274]=['',1.66942,0.039385,6.31];
ybs[2275]=['',1.655141,-0.8791149,7.04];
ybs[2276]=['ζ CMa',1.6632012,-0.5249018,3.02];
ybs[2277]=['',1.6351076,-1.2515918,6.64];
ybs[2278]=['',1.6688002,-0.2056923,5.64];
ybs[2279]=['',1.7051948,1.2307969,5.97];
ybs[2280]=['μ Gem',1.676911,0.3927115,2.88];
ybs[2281]=['',1.6749603,0.2191666,6];
ybs[2282]=['',1.6642576,-0.5961235,5.53];
ybs[2283]=['ψ1 Aur',1.6869916,0.8599936,4.91];
ybs[2284]=['',1.6610553,-0.8508872,6.6];
ybs[2285]=['',1.6944028,0.982098,5.64];
ybs[2286]=['',1.6776727,0.0654744,6.4];
ybs[2287]=['5 Lyn',1.6963595,1.0193084,5.21];
ybs[2288]=['β CMa',1.6741704,-0.3136093,1.98];
ybs[2289]=['',1.6776567,-0.0820354,6.67];
ybs[2290]=['δ Col',1.6708787,-0.5837896,3.85];
ybs[2291]=['',1.6856123,0.5182451,6.71];
ybs[2292]=['ε Mon',1.6797113,0.0799271,4.44];
ybs[2293]=['',1.6797405,0.0799755,6.72];
ybs[2294]=['',1.6810546,0.154833,6.26];
ybs[2295]=['',1.6784238,-0.1725714,6.19];
ybs[2296]=['',1.6850199,0.2800035,6.33];
ybs[2297]=['',1.678944,-0.2632813,6.24];
ybs[2298]=['',1.68824,0.4068823,6.06];
ybs[2299]=['',1.6808532,-0.2014764,5.22];
ybs[2300]=['',1.6788637,-0.3455494,6.6];
ybs[2301]=['',1.6758836,-0.5550655,6.34];
ybs[2302]=['',1.6875251,0.2566979,6.24];
ybs[2303]=['',1.6815381,-0.226475,6.12];
ybs[2304]=['',1.6861336,0.1234253,5.98];
ybs[2305]=['',1.6792088,-0.4466486,5.63];
ybs[2306]=['',1.6863059,0.0259529,6.66];
ybs[2307]=['',1.6860676,-0.0167588,5.87];
ybs[2308]=['',1.6997125,0.8271051,6.56];
ybs[2309]=['',1.6883859,0.0394069,6.51];
ybs[2310]=['',1.6789993,-0.640904,5.62];
ybs[2311]=['',1.6881884,-0.0681296,6.35];
ybs[2312]=['',1.6825463,-0.5025454,6.39];
ybs[2313]=['',1.6976046,0.568063,6.43];
ybs[2314]=['ν Pic',1.6726081,-0.9840624,5.61];
ybs[2315]=['',1.6888858,-0.1380316,6.4];
ybs[2316]=['',1.6761242,-0.9109596,5.98];
ybs[2317]=['',1.6819416,-0.7033303,6.31];
ybs[2318]=['',1.6920968,-0.0265651,5.87];
ybs[2319]=['',1.6916085,-0.0804949,6.15];
ybs[2320]=['α Car',1.6774828,-0.9199462,-0.72];
ybs[2321]=['',1.6935748,0.014413,6.71];
ybs[2322]=['',1.6922531,-0.1313582,6.27];
ybs[2323]=['',1.6855605,-0.6122268,6.25];
ybs[2324]=['16 Gem',1.6985644,0.3574524,6.22];
ybs[2325]=['6 Lyn',1.7136463,1.01483,5.88];
ybs[2326]=['48 Aur',1.7017575,0.5319255,5.55];
ybs[2327]=['',1.6952313,0.0504893,5.55];
ybs[2328]=['',1.6946524,0.0049568,5.2];
ybs[2329]=['',1.6947615,-0.0050839,5.55];
ybs[2330]=['',1.6760399,-1.0220117,6.48];
ybs[2331]=['',1.67188,-1.1116989,6.27];
ybs[2332]=['47 Aur',1.7092431,0.8145229,5.9];
ybs[2333]=['',1.7032199,0.4703945,6.47];
ybs[2334]=['',1.7006686,0.2831354,6.23];
ybs[2335]=['',1.6811886,-0.9218839,6.51];
ybs[2336]=['',1.6997622,0.1795615,6.15];
ybs[2337]=['ν Gem',1.7030023,0.3524879,4.15];
ybs[2338]=['10 Mon',1.6976455,-0.0833878,5.06];
ybs[2339]=['',1.6777347,-1.0523359,5.8];
ybs[2340]=['',1.7632796,1.3888721,6.54];
ybs[2341]=['',1.6992949,0.0331,6.48];
ybs[2342]=['',1.6856537,-0.8410989,5.76];
ybs[2343]=['',1.6934,-0.4515471,6.07];
ybs[2344]=['',1.7854432,1.4327382,6.65];
ybs[2345]=['',1.7027863,0.1920438,6.59];
ybs[2346]=['π1 Dor',1.6686246,-1.2216687,5.56];
ybs[2347]=['',1.6924933,-0.661664,6.48];
ybs[2348]=['',1.6781128,-1.1072764,6.46];
ybs[2349]=['',1.7035479,0.0458996,6.16];
ybs[2350]=['β Mon',1.7012967,-0.1230244,4.6];
ybs[2351]=['β Mon',1.701333,-0.1230536,5.4];
ybs[2352]=['β Mon',1.701333,-0.1230536,5.6];
ybs[2353]=['',1.7000309,-0.3051181,5.77];
ybs[2354]=['',1.6801677,-1.1142429,6.27];
ybs[2355]=['λ CMa',1.697349,-0.5689003,4.48];
ybs[2356]=['',1.7074717,0.1572966,6.57];
ybs[2357]=['',1.7625534,1.3608837,5.73];
ybs[2358]=['',1.6994767,-0.5652591,5.74];
ybs[2359]=['',1.7486228,1.2858567,6.24];
ybs[2360]=['',1.7124846,0.2953319,6.2];
ybs[2361]=['',1.707164,-0.1762456,5.93];
ybs[2362]=['',1.6991865,-0.7171608,6.32];
ybs[2363]=['',1.6904533,-1.0125891,5.82];
ybs[2364]=['',1.7122123,0.1960619,6.14];
ybs[2365]=['19 Gem',1.7144272,0.2772587,6.4];
ybs[2366]=['',1.7188078,0.5661263,5.87];
ybs[2367]=['',1.7087449,-0.2297726,6.16];
ybs[2368]=['',1.7143933,0.2055062,6.65];
ybs[2369]=['',1.7150453,0.2011802,5.23];
ybs[2370]=['7 Lyn',1.729691,0.965756,6.45];
ybs[2371]=['π2 Dor',1.6811179,-1.2165664,5.38];
ybs[2372]=['',1.7175955,0.203429,6.03];
ybs[2373]=['',1.7122877,-0.2165788,5.15];
ybs[2374]=['',1.7089539,-0.4849651,5.93];
ybs[2375]=['',1.7144238,-0.1426927,5.43];
ybs[2376]=['12 Mon',1.7170341,0.0844374,5.84];
ybs[2377]=['',1.7243275,0.576053,6.42];
ybs[2378]=['',1.7033293,-0.8771252,5.27];
ybs[2379]=['13 Mon',1.7196739,0.1276675,4.5];
ybs[2380]=['',1.7169115,-0.1027445,5.6];
ybs[2381]=['ξ1 CMa',1.7138838,-0.4090343,4.33];
ybs[2382]=['',1.710505,-0.6156892,5.84];
ybs[2383]=['',1.7011163,-0.9925505,5.22];
ybs[2384]=['',1.7091988,-0.7144239,6.2];
ybs[2385]=['',1.7229879,0.2467308,5.53];
ybs[2386]=['',1.7184377,-0.1952069,6.24];
ybs[2387]=['',1.7120079,-0.6450286,6.34];
ybs[2388]=['8 Lyn',1.7442461,1.0726802,5.94];
ybs[2389]=['',1.7225275,-0.021623,5.1];
ybs[2390]=['',1.7592176,1.2518572,5.92];
ybs[2391]=['',1.7169365,-0.5593531,5.69];
ybs[2392]=['49 Aur',1.7305731,0.4887388,5.27];
ybs[2393]=['',1.7153305,-0.6582423,5.24];
ybs[2394]=['',1.7096666,-0.9048365,5.6];
ybs[2395]=['',1.7891615,1.3882105,5.45];
ybs[2396]=['11 Lyn',1.7433704,0.9919842,5.85];
ybs[2397]=['',1.720925,-0.3655132,6.4];
ybs[2398]=['14 Mon',1.7278421,0.1318289,6.45];
ybs[2399]=['',1.7370123,0.6706423,5.29];
ybs[2400]=['',1.7302042,0.1739881,5.88];
ybs[2401]=['',1.718838,-0.6744572,6.44];
ybs[2402]=['',1.702173,-1.1446694,6.29];
ybs[2403]=['',1.7297265,0.0151929,5.8];
ybs[2404]=['',1.7078075,-1.0803024,6.15];
ybs[2405]=['',1.7218327,-0.6326968,5.42];
ybs[2406]=['μ Pic',1.7117687,-1.0257593,5.7];
ybs[2407]=['',1.7330771,0.0781486,6.55];
ybs[2408]=['ξ2 CMa',1.7278759,-0.4011476,4.54];
ybs[2409]=['',1.7253595,-0.5713413,5.62];
ybs[2410]=['',1.7189137,-0.9136316,6.19];
ybs[2411]=['',1.7402353,0.4288288,6.44];
ybs[2412]=['',1.7352776,-0.0913038,5.52];
ybs[2413]=['51 Aur',1.7462749,0.6871253,5.69];
ybs[2414]=['ψ3 Aur',1.7470125,0.6960541,5.2];
ybs[2415]=['γ Gem',1.7410102,0.2858553,1.93];
ybs[2416]=['',1.739255,0.1067199,6.06];
ybs[2417]=['ν1 CMa',1.7338365,-0.3260286,5.7];
ybs[2418]=['',1.7286819,-0.6422722,5.59];
ybs[2419]=['53 Aur',1.7444999,0.505498,5.79];
ybs[2420]=['',1.7403674,0.1890585,6.38];
ybs[2421]=['ψ2 Aur',1.749416,0.7411901,4.79];
ybs[2422]=['',1.7357926,-0.2328466,5.97];
ybs[2423]=['ν2 CMa',1.7351352,-0.3364308,3.95];
ybs[2424]=['',1.7403032,0.0468333,6.17];
ybs[2425]=['',1.7308971,-0.6302147,6.35];
ybs[2426]=['',1.74129,0.0861545,6.15];
ybs[2427]=['',1.7349828,-0.3950591,6.35];
ybs[2428]=['',1.7522936,0.7678003,6.41];
ybs[2429]=['',1.725577,-0.9249326,4.39];
ybs[2430]=['',1.7472541,0.3841332,6.04];
ybs[2431]=['',1.7397478,-0.2269937,6.12];
ybs[2432]=['54 Aur',1.7495615,0.4929013,6.03];
ybs[2433]=['',1.7492662,0.4289694,6.38];
ybs[2434]=['',1.743029,-0.0447638,6.14];
ybs[2435]=['',1.7454107,0.081666,6.57];
ybs[2436]=['',1.7444653,0.0277905,6.21];
ybs[2437]=['ν3 CMa',1.7404436,-0.3186689,4.43];
ybs[2438]=['',1.7357143,-0.6661402,6.04];
ybs[2439]=['',1.734725,-0.7256589,6.34];
ybs[2440]=['',1.7366492,-0.6459643,5.71];
ybs[2441]=['',1.7393671,-0.5647974,5.27];
ybs[2442]=['',1.7435548,-0.2948712,6.03];
ybs[2443]=['',1.749955,0.2262135,5.97];
ybs[2444]=['',1.746669,-0.2472691,4.82];
ybs[2445]=['ν Pup',1.7385693,-0.7542762,3.17];
ybs[2446]=['',1.7590227,0.6267287,6.46];
ybs[2447]=['25 Gem',1.7574033,0.4917208,6.42];
ybs[2448]=['',1.7529133,0.1108163,6.51];
ybs[2449]=['',1.7476921,-0.4139458,6.05];
ybs[2450]=['15 Mon',1.7550044,0.1723154,4.66];
ybs[2451]=['',1.7569417,0.2857919,6.28];
ybs[2452]=['',1.7563927,0.1916468,6.11];
ybs[2453]=['ψ4 Aur',1.7659246,0.776682,5.02];
ybs[2454]=['',1.7478358,-0.5321877,5.71];
ybs[2455]=['',1.7551352,0.0082489,5.79];
ybs[2456]=['',1.7419561,-0.8419722,4.93];
ybs[2457]=['',1.7715327,0.9297701,6.27];
ybs[2458]=['',1.7660797,0.6479194,6.19];
ybs[2459]=['',1.7484435,-0.666381,6.58];
ybs[2460]=['26 Gem',1.7615374,0.3075598,5.21];
ybs[2461]=['',1.7592752,0.1103371,6.37];
ybs[2462]=['',1.7376998,-1.0743165,6.18];
ybs[2463]=['',1.7584786,-0.1604011,5.19];
ybs[2464]=['12 Lyn',1.781189,1.0370055,4.87];
ybs[2465]=['',1.770294,0.6298073,6.31];
ybs[2466]=['ε Gem',1.7685218,0.4381976,2.98];
ybs[2467]=['',1.7640369,0.0525274,6.19];
ybs[2468]=['',1.7539222,-0.7046302,6.12];
ybs[2469]=['',1.7515832,-0.8324711,6.65];
ybs[2470]=['13 Lyn',1.783454,0.9973377,5.35];
ybs[2471]=['30 Gem',1.768267,0.2304453,4.49];
ybs[2472]=['',1.7664124,0.0682108,5.9];
ybs[2473]=['28 Gem',1.7723092,0.5052055,5.44];
ybs[2474]=['',1.7615305,-0.3922219,6.13];
ybs[2475]=['',1.7585767,-0.6705868,6.29];
ybs[2476]=['ψ5 Aur',1.7817973,0.7601206,5.25];
ybs[2477]=['ξ Gem',1.773933,0.2246348,3.36];
ybs[2478]=['',1.7893362,0.9717607,6.33];
ybs[2479]=['',1.7892925,0.9717608,6.28];
ybs[2480]=['ψ6 Aur',1.7862236,0.851077,5.22];
ybs[2481]=['',1.7633943,-0.6844676,6.3];
ybs[2482]=['32 Gem',1.776601,0.2211044,6.46];
ybs[2483]=['42 Cam',1.8033035,1.17886,5.14];
ybs[2484]=['α CMa',1.7721807,-0.2921838,-1.46];
ybs[2485]=['10 CMa',1.7685869,-0.5427091,5.2];
ybs[2486]=['',1.7704796,-0.4776268,6.45];
ybs[2487]=['16 Mon',1.7792198,0.1494287,5.93];
ybs[2488]=['',1.7729415,-0.4099227,6.05];
ybs[2489]=['',1.7710846,-0.5342592,6.54];
ybs[2490]=['',1.7759239,-0.2586814,5.32];
ybs[2491]=['',1.7833155,0.3170786,6.2];
ybs[2492]=['',1.7725068,-0.5553372,5.92];
ybs[2493]=['',1.7731685,-0.5405954,5.8];
ybs[2494]=['',1.7789952,-0.1768512,5.66];
ybs[2495]=['17 Mon',1.7826456,0.139822,4.77];
ybs[2496]=['11 CMa',1.7797061,-0.2522269,5.29];
ybs[2497]=['',1.748025,-1.2531069,6.51];
ybs[2498]=['18 Mon',1.7847463,0.0416425,4.47];
ybs[2499]=['',1.7750123,-0.6905431,6.62];
ybs[2500]=['',1.7832631,-0.1575065,5.07];
ybs[2501]=['12 CMa',1.7801909,-0.3672406,6.08];
ybs[2502]=['',1.7757606,-0.6597492,6.21];
ybs[2503]=['43 Cam',1.8156494,1.2018094,5.12];
ybs[2504]=['',1.7939984,0.5686166,5.71];
ybs[2505]=['',1.7713116,-0.9115143,6.57];
ybs[2506]=['',1.7865996,-0.0234866,5.75];
ybs[2507]=['',1.7732989,-0.9151644,5.8];
ybs[2508]=['ψ7 Aur',1.7992354,0.7287351,5.02];
ybs[2509]=['',1.7899368,0.0170173,6.15];
ybs[2510]=['',1.7807928,-0.6624506,5.26];
ybs[2511]=['33 Gem',1.7938718,0.2823142,5.85];
ybs[2512]=['14 Lyn',1.8110542,1.0370629,5.33];
ybs[2513]=['',1.7907383,-0.0401246,5.74];
ybs[2514]=['',1.7888924,-0.2647936,5.39];
ybs[2515]=['',1.7777113,-0.8952042,5.4];
ybs[2516]=['',1.7765348,-0.9550475,6.46];
ybs[2517]=['35 Gem',1.7963566,0.2336237,5.65];
ybs[2518]=['',1.7791451,-0.9698061,5.61];
ybs[2519]=['',1.8470188,1.3429294,4.55];
ybs[2520]=['',1.7918765,-0.4206777,6.33];
ybs[2521]=['36 Gem',1.8016069,0.3793091,5.27];
ybs[2522]=['',1.7976024,-0.0099256,5.77];
ybs[2523]=['',1.7590466,-1.2765621,6.37];
ybs[2524]=['',1.8097498,0.782086,6.26];
ybs[2525]=['',1.8036503,0.4114287,5.65];
ybs[2526]=['',1.7967686,-0.1408287,6.29];
ybs[2527]=['',1.7949308,-0.2986514,5.79];
ybs[2528]=['',1.765856,-1.2297281,6.11];
ybs[2529]=['',1.7933272,-0.4775447,7.04];
ybs[2530]=['κ CMa',1.7919427,-0.5678579,3.96];
ybs[2531]=['59 Aur',1.8089062,0.6778864,6.12];
ybs[2532]=['θ Gem',1.8075924,0.5922272,3.6];
ybs[2533]=['60 Aur',1.809745,0.6703603,6.3];
ybs[2534]=['',1.8088635,0.6241205,6.01];
ybs[2535]=['',1.8013282,0.052593,6.38];
ybs[2536]=['',1.7956473,-0.4503951,6.33];
ybs[2537]=['',1.7943704,-0.5538568,5.7];
ybs[2538]=['',1.7916656,-0.793728,6.55];
ybs[2539]=['ψ8 Aur',1.8129184,0.671522,6.48];
ybs[2540]=['',1.7913427,-0.8140557,5.14];
ybs[2541]=['',1.7963462,-0.6003063,4.99];
ybs[2542]=['α Pic',1.7820836,-1.0815386,3.27];
ybs[2543]=['',1.8066261,0.1457583,5.77];
ybs[2544]=['',1.8041877,-0.0932841,6.3];
ybs[2545]=['τ Pup',1.7911172,-0.8838688,2.93];
ybs[2546]=['',1.7904629,-0.9363584,4.4];
ybs[2547]=['',1.809133,0.1914128,6.24];
ybs[2548]=['',1.8190983,0.7992918,6.34];
ybs[2549]=['',1.8189207,0.7658447,6.13];
ybs[2550]=['',1.7998722,-0.6328304,5.96];
ybs[2551]=['ζ Men',1.737219,-1.4108347,5.64];
ybs[2552]=['15 Lyn',1.8291779,1.0191153,4.35];
ybs[2553]=['',1.8288253,1.0041207,6.05];
ybs[2554]=['',1.7903384,-1.052021,6.11];
ybs[2555]=['',1.798351,-0.8433535,6.42];
ybs[2556]=['38 Gem',1.8147583,0.2294734,4.65];
ybs[2557]=['',1.8077172,-0.3326932,5.64];
ybs[2558]=['',1.8079323,-0.330958,6.14];
ybs[2559]=['',1.806009,-0.4710024,6.4];
ybs[2560]=['ψ9 Aur',1.8247363,0.807095,5.87];
ybs[2561]=['37 Gem',1.8181689,0.4423581,5.73];
ybs[2562]=['',1.8118473,-0.1026622,6.41];
ybs[2563]=['15 CMa',1.808693,-0.353489,4.83];
ybs[2564]=['',1.8131941,-0.0201884,5.45];
ybs[2565]=['',1.8264925,0.8146108,5.86];
ybs[2566]=['θ CMa',1.8118229,-0.2106304,4.07];
ybs[2567]=['',1.8036567,-0.7423437,6.52];
ybs[2568]=['',1.8083926,-0.4986226,6.04];
ybs[2569]=['',1.8144435,-0.031177,6.21];
ybs[2570]=['',1.8101349,-0.4288033,6.21];
ybs[2571]=['',1.804092,-0.7680254,6.46];
ybs[2572]=['ο1 CMa',1.8110673,-0.4226044,3.87];
ybs[2573]=['',1.8495965,1.2352437,5.68];
ybs[2574]=['',1.815619,-0.0494569,6.04];
ybs[2575]=['',1.8114497,-0.4181449,6.91];
ybs[2576]=['',1.8186368,0.1447631,6.29];
ybs[2577]=['16 Lyn',1.8293568,0.7864903,4.9];
ybs[2578]=['',1.825984,0.5873012,5.89];
ybs[2579]=['',1.8031913,-0.9445499,6.57];
ybs[2580]=['17 CMa',1.815217,-0.3566542,5.74];
ybs[2581]=['',1.8224218,0.1732332,5.92];
ybs[2582]=['π CMa',1.8177523,-0.3519762,4.68];
ybs[2583]=['',1.811451,-0.7399362,6.32];
ybs[2584]=['',1.8024302,-1.0361982,6.41];
ybs[2585]=['μ CMa',1.8201233,-0.245642,5];
ybs[2586]=['',1.8090251,-0.8838535,6.26];
ybs[2587]=['',1.8183159,-0.4009341,5.3];
ybs[2588]=['ι CMa',1.8201195,-0.2981862,4.37];
ybs[2589]=['',1.8268521,0.2072772,6.27];
ybs[2590]=['',1.818481,-0.5553722,6.36];
ybs[2591]=['',1.8242398,-0.1432917,6.34];
ybs[2592]=['62 Aur',1.8351201,0.6635725,6];
ybs[2593]=['39 Gem',1.8333814,0.4546402,6.1];
ybs[2594]=['ι Vol',1.7941382,-1.2390293,5.4];
ybs[2595]=['',1.8247489,-0.3880661,6.61];
ybs[2596]=['',1.8219943,-0.6173681,6.29];
ybs[2597]=['40 Gem',1.8363182,0.4517202,6.4];
ybs[2598]=['',1.8320192,0.132469,6.27];
ybs[2599]=['',1.8260316,-0.430432,5.46];
ybs[2600]=['',1.8189215,-0.8508779,4.95];
ybs[2601]=['',2.0516689,1.5157628,5.07];
ybs[2602]=['',1.8331784,0.062309,5.97];
ybs[2603]=['',1.8265168,-0.4811688,6.23];
ybs[2604]=['',1.8243059,-0.6202674,6.23];
ybs[2605]=['',1.8349964,0.1271395,6.35];
ybs[2606]=['',1.8283674,-0.4746665,6.37];
ybs[2607]=['41 Gem',1.8393838,0.2800554,5.68];
ybs[2608]=['',1.8305063,-0.4441132,5.59];
ybs[2609]=['',1.8691928,1.2338741,6.5];
ybs[2610]=['ε CMa',1.8304566,-0.5062179,1.5];
ybs[2611]=['',1.8292929,-0.595916,5.06];
ybs[2612]=['',1.8445785,0.5651544,6.59];
ybs[2613]=['',1.8308152,-0.5415714,6.42];
ybs[2614]=['',1.8387376,-0.0942446,6.3];
ybs[2615]=['',1.8352853,-0.3776164,6.26];
ybs[2616]=['',1.8390386,-0.1473034,5.96];
ybs[2617]=['',1.8374412,-0.3524106,6.31];
ybs[2618]=['',1.8297678,-0.7993601,6.22];
ybs[2619]=['',1.8401431,-0.1612004,6.49];
ybs[2620]=['',1.8381728,-0.3866303,6.53];
ybs[2621]=['',1.8451783,0.083504,6.63];
ybs[2622]=['ω Gem',1.8491031,0.422042,5.18];
ybs[2623]=['',1.8488816,0.3092987,5.94];
ybs[2624]=['',1.8481959,0.2670728,5.74];
ybs[2625]=['',1.8461948,0.0964077,6.59];
ybs[2626]=['',1.8286756,-0.9732178,6.27];
ybs[2627]=['',1.8494196,0.2904238,5.82];
ybs[2628]=['',1.8457877,-0.0240728,6.17];
ybs[2629]=['',1.8395801,-0.4978111,6.27];
ybs[2630]=['',1.8283492,-0.9848285,6.45];
ybs[2631]=['',1.8458807,-0.1004604,5.2];
ybs[2632]=['',1.8414265,-0.4406698,5.63];
ybs[2633]=['',1.8398543,-0.5846566,6.4];
ybs[2634]=['',1.8674966,1.0431106,6.44];
ybs[2635]=['',1.8541081,0.5114265,5.93];
ybs[2636]=['',1.8650921,0.9201758,6.12];
ybs[2637]=['',1.8624069,0.8332103,6.38];
ybs[2638]=['σ CMa',1.8440075,-0.4881387,3.47];
ybs[2639]=['',1.8523566,0.1588922,5.97];
ybs[2640]=['19 Mon',1.8501827,-0.0745852,4.99];
ybs[2641]=['',1.8538867,0.1905376,5.13];
ybs[2642]=['ζ Gem',1.8563406,0.3584094,3.79];
ybs[2643]=['',1.8549384,0.2192074,5.98];
ybs[2644]=['',1.838729,-0.8977192,5.14];
ybs[2645]=['ο2 CMa',1.8498968,-0.4165681,3.02];
ybs[2646]=['',1.8566036,0.0253653,6.57];
ybs[2647]=['',1.8552626,-0.093523,5.62];
ybs[2648]=['',1.8545087,-0.1773071,6.45];
ybs[2649]=['γ CMa',1.8534377,-0.2734582,4.12];
ybs[2650]=['',1.8454963,-0.7581354,6.43];
ybs[2651]=['44 Gem',1.8616556,0.3944733,6.02];
ybs[2652]=['',1.8660875,0.6010534,5.55];
ybs[2653]=['',1.8388583,-1.0292745,6.02];
ybs[2654]=['',1.8317265,-1.185924,5.17];
ybs[2655]=['',1.8626257,0.1596996,5.78];
ybs[2656]=['',1.8576561,-0.385149,6.09];
ybs[2657]=['',1.871202,0.5929366,5.91];
ybs[2658]=['',1.8533694,-0.7395302,5.2];
ybs[2659]=['',1.8528821,-0.7617095,5.54];
ybs[2660]=['',1.8529909,-0.7617728,6.79];
ybs[2661]=['',1.8711077,0.491145,6.48];
ybs[2662]=['',1.8626632,-0.1866957,6.49];
ybs[2663]=['',1.8706008,0.3956133,7.68];
ybs[2664]=['',1.8521538,-0.8660064,4.93];
ybs[2665]=['',1.8748998,0.5898358,6.28];
ybs[2666]=['',1.8483246,-1.0334493,5.5];
ybs[2667]=['',1.8767709,0.6528869,6.16];
ybs[2668]=['',1.8688151,0.085064,6.11];
ybs[2669]=['',1.8603144,-0.6076071,6.14];
ybs[2670]=['',1.8663485,-0.1977524,5.39];
ybs[2671]=['',1.8659575,-0.2169454,6.48];
ybs[2672]=['',1.8626084,-0.5356654,6.34];
ybs[2673]=['',1.9048687,1.2527321,6.35];
ybs[2674]=['',1.8720438,0.1297524,5.75];
ybs[2675]=['',1.8532061,-0.9910768,5.17];
ybs[2676]=['45 Gem',1.8747292,0.2773971,5.44];
ybs[2677]=['',1.8623141,-0.6705306,6.11];
ybs[2678]=['',1.8666359,-0.436277,6.08];
ybs[2679]=['',1.8581106,-0.8795696,6.46];
ybs[2680]=['',1.8671267,-0.4659003,6.62];
ybs[2681]=['θ Men',1.8112999,-1.3866728,5.45];
ybs[2682]=['',1.8688894,-0.4167291,5.71];
ybs[2683]=['',1.866874,-0.7143578,5.79];
ybs[2684]=['',1.882557,0.3701646,6.43];
ybs[2685]=['δ CMa',1.8732058,-0.4612975,1.84];
ybs[2686]=['',1.8779922,-0.1812493,6.21];
ybs[2687]=['',1.8751873,-0.4203009,6.65];
ybs[2688]=['63 Aur',1.8901828,0.6855937,4.9];
ybs[2689]=['τ Gem',1.8874465,0.5272054,4.41];
ybs[2690]=['',1.8664615,-0.9076434,5.96];
ybs[2691]=['',1.8787039,-0.2840074,6.03];
ybs[2692]=['47 Gem',1.8883557,0.468061,5.78];
ybs[2693]=['20 Mon',1.8821019,-0.0746181,4.92];
ybs[2694]=['',1.8745188,-0.6927753,4.83];
ybs[2695]=['',1.8986493,0.8969073,5.47];
ybs[2696]=['',1.879034,-0.441025,5.69];
ybs[2697]=['',1.8812317,-0.3267829,6.23];
ybs[2698]=['48 Gem',1.8928362,0.4204331,5.85];
ybs[2699]=['21 Mon',1.8873286,-0.0059452,5.45];
ybs[2700]=['',1.8815868,-0.4804798,5.46];
ybs[2701]=['',1.9614618,1.4173968,6.31];
ybs[2702]=['',1.8895597,0.0980138,6.09];
ybs[2703]=['',1.8946172,0.4744717,6.43];
ybs[2704]=['',1.8594142,-1.2020589,6.47];
ybs[2705]=['',1.8907238,0.0948698,6.16];
ybs[2706]=['δ Mon',1.8893797,-0.0092802,4.15];
ybs[2707]=['18 Lyn',1.9105907,1.0401509,5.2];
ybs[2708]=['',1.8878553,-0.3651553,5.84];
ybs[2709]=['51 Gem',1.8965723,0.2813319,5];
ybs[2710]=['26 CMa',1.8898631,-0.4534636,5.92];
ybs[2711]=['',1.8823241,-0.8546958,5.14];
ybs[2712]=['',1.8890449,-0.5386198,6.1];
ybs[2713]=['',1.9089522,0.8237757,5.58];
ybs[2714]=['',1.9015999,0.4305815,6.89];
ybs[2715]=['',1.8944608,-0.1970643,5.78];
ybs[2716]=['',1.8906591,-0.480198,6.59];
ybs[2717]=['52 Gem',1.9027199,0.433619,5.82];
ybs[2718]=['',1.8903159,-0.638504,5.96];
ybs[2719]=['',1.8893609,-0.7075204,5.31];
ybs[2720]=['',1.9015283,0.2107569,5.62];
ybs[2721]=['',1.9002802,0.0536019,5.35];
ybs[2722]=['',1.8952216,-0.3964169,6.01];
ybs[2723]=['',1.899372,-0.0687925,5.75];
ybs[2724]=['',1.8994772,-0.1743175,5.9];
ybs[2725]=['',1.8969866,-0.4004882,6.36];
ybs[2726]=['',1.8959169,-0.4781533,6.12];
ybs[2727]=['γ1 Vol',1.8696985,-1.2310533,5.69];
ybs[2728]=['γ2 Vol',1.8698946,-1.2310828,3.78];
ybs[2729]=['',1.916901,0.9091209,5.92];
ybs[2730]=['53 Gem',1.9083221,0.4861855,5.71];
ybs[2731]=['',1.9003871,-0.1807626,6.03];
ybs[2732]=['',1.8902211,-0.8167899,4.49];
ybs[2733]=['',1.8965254,-0.5432119,6.6];
ybs[2734]=['',1.9882279,1.4374838,4.96];
ybs[2735]=['',1.8972966,-0.5302301,6.33];
ybs[2736]=['24 Mon',1.9044743,-0.0035277,6.41];
ybs[2737]=['27 CMa',1.8987925,-0.460638,4.66];
ybs[2738]=['',1.8932358,-0.7892828,4.89];
ybs[2739]=['',1.9062314,0.1385243,5.82];
ybs[2740]=['',1.8946577,-0.7798027,5.1];
ybs[2741]=['ω CMa',1.9012105,-0.4679783,3.85];
ybs[2742]=['',1.901373,-0.4726086,5.58];
ybs[2743]=['',1.9208901,0.8625849,5.05];
ybs[2744]=['',1.9058248,-0.1854377,5.95];
ybs[2745]=['64 Aur',1.9181078,0.7128117,5.78];
ybs[2746]=['',1.8860156,-1.1035504,6.02];
ybs[2747]=['',1.905622,-0.4150651,6.32];
ybs[2748]=['',1.9033763,-0.5362883,5.36];
ybs[2749]=['',1.917699,0.5395443,6.24];
ybs[2750]=['',1.9079184,-0.2727427,5.46];
ybs[2751]=['',1.9010392,-0.7237228,5.94];
ybs[2752]=['',1.9133391,0.1158688,6.65];
ybs[2753]=['',1.8998423,-0.8183854,5.72];
ybs[2754]=['',1.8991688,-0.8432018,4.76];
ybs[2755]=['λ Gem',1.917185,0.2879459,3.58];
ybs[2756]=['',1.9092265,-0.4076548,4.79];
ybs[2757]=['',1.9138671,-0.1173185,6.29];
ybs[2758]=['',1.9088931,-0.4873383,4.64];
ybs[2759]=['',1.901907,-0.9170013,5.97];
ybs[2760]=['',1.9103649,-0.5399729,6.32];
ybs[2761]=['',1.9081206,-0.6695109,5.8];
ybs[2762]=['',1.9095072,-0.6393873,5.03];
ybs[2763]=['',1.9063689,-0.8170851,5.66];
ybs[2764]=['47 Cam',1.9383938,1.0447092,6.35];
ybs[2765]=['π Pup',1.9108686,-0.6481992,2.7];
ybs[2766]=['',1.9142459,-0.4684367,6.46];
ybs[2767]=['',1.9313684,0.7437156,6.35];
ybs[2768]=['',1.9325921,0.7886119,5.77];
ybs[2769]=['δ Gem',1.9262594,0.3829073,3.53];
ybs[2770]=['',1.9222571,0.0470843,5.89];
ybs[2771]=['',1.9242486,0.1239137,5.91];
ybs[2772]=['',1.9259502,0.2635368,6.45];
ybs[2773]=['29 CMa',1.9181569,-0.4293735,4.98];
ybs[2774]=['τ CMa',1.9182925,-0.4362727,4.4];
ybs[2775]=['19 Lyn',1.9402668,0.9641142,6.53];
ybs[2776]=['19 Lyn',1.9403536,0.9640608,5.45];
ybs[2777]=['',1.919956,-0.3372478,6.09];
ybs[2778]=['',1.9188604,-0.4647519,5.28];
ybs[2779]=['',1.9159678,-0.6418683,4.66];
ybs[2780]=['',1.9219687,-0.2868942,5.7];
ybs[2781]=['',1.9144833,-0.7684458,5.85];
ybs[2782]=['',1.917408,-0.6420215,5.11];
ybs[2783]=['',1.9169305,-0.6850867,5.25];
ybs[2784]=['',1.9362491,0.679836,6.4];
ybs[2785]=['65 Aur',1.9353362,0.64082,5.13];
ybs[2786]=['',1.9201534,-0.5893957,6.3];
ybs[2787]=['56 Gem',1.9341495,0.3560372,5.1];
ybs[2788]=['',1.9286008,-0.2513904,5.45];
ybs[2789]=['',2.0012806,1.4110182,6.41];
ybs[2790]=['',1.9301655,-0.1557201,6.55];
ybs[2791]=['',1.9278979,-0.3995968,6.61];
ybs[2792]=['',1.9278498,-0.4713638,6.01];
ybs[2793]=['',1.9338807,0.0023218,5.99];
ybs[2794]=['',1.9285742,-0.4526515,5.87];
ybs[2795]=['δ Vol',1.9059578,-1.1867965,3.98];
ybs[2796]=['',1.9490986,0.9048032,5.8];
ybs[2797]=['66 Aur',1.9447151,0.7090725,5.19];
ybs[2798]=['',1.9334419,-0.1574868,6.43];
ybs[2799]=['',1.9348587,-0.0527648,6.23];
ybs[2800]=['57 Gem',1.9410091,0.4364298,5.03];
ybs[2801]=['',1.9619286,1.1568823,6.47];
ybs[2802]=['58 Gem',1.9408976,0.3996858,6.02];
ybs[2803]=['',1.9352527,-0.1051936,5.82];
ybs[2804]=['',1.9338965,-0.3326755,4.96];
ybs[2805]=['',1.9237915,-0.9137642,6.05];
ybs[2806]=['',1.9238135,-0.9137303,6.6];
ybs[2807]=['',1.9250724,-0.9098301,5.39];
ybs[2808]=['59 Gem',1.9458416,0.4815804,5.76];
ybs[2809]=['',1.9449264,0.2700333,6.41];
ybs[2810]=['21 Lyn',1.9565355,0.8580857,4.64];
ybs[2811]=['',1.9367536,-0.5579555,5.43];
ybs[2812]=['1 CMi',1.9469992,0.2028774,5.3];
ybs[2813]=['ι Gem',1.9509518,0.4843626,3.79];
ybs[2814]=['',1.9390213,-0.4865808,5.38];
ybs[2815]=['',1.9390171,-0.5628179,5.39];
ybs[2816]=['',1.9407436,-0.5281718,6.6];
ybs[2817]=['',1.9446706,-0.2835567,5.33];
ybs[2818]=['',1.9427355,-0.4006937,6.19];
ybs[2819]=['η CMa',1.9416162,-0.5122231,2.45];
ybs[2820]=['ε CMi',1.9498707,0.1610951,4.99];
ybs[2821]=['',1.9407492,-0.6262741,6.31];
ybs[2822]=['',1.9774528,1.1940948,5.64];
ybs[2823]=['',1.9453316,-0.3326213,6.24];
ybs[2824]=['',1.9468152,-0.2408149,5.78];
ybs[2825]=['',1.9502129,-0.1015976,5.97];
ybs[2826]=['',1.9442679,-0.5559636,5.35];
ybs[2827]=['',1.9555242,0.3750566,6.54];
ybs[2828]=['',1.9534822,0.1843394,6.37];
ybs[2829]=['61 Gem',1.9559158,0.3527397,5.93];
ybs[2830]=['',1.9511662,-0.080001,6.76];
ybs[2831]=['',1.9473354,-0.3844715,6.05];
ybs[2832]=['',1.9544712,0.1913333,6.41];
ybs[2833]=['',1.9475913,-0.4409336,5.78];
ybs[2834]=['',1.9442273,-0.6516273,6.97];
ybs[2835]=['',1.9442344,-0.6516418,6.84];
ybs[2836]=['',1.9658101,0.8401335,5.72];
ybs[2837]=['β CMi',1.9563862,0.1438612,2.9];
ybs[2838]=['63 Gem',1.9594607,0.3734633,5.22];
ybs[2839]=['',1.9485785,-0.5547459,6.31];
ybs[2840]=['',1.7393081,-1.5172312,6.47];
ybs[2841]=['22 Lyn',1.9706102,0.8661052,5.36];
ybs[2842]=['',1.953141,-0.4146679,6.56];
ybs[2843]=['η CMi',1.9601989,0.1203352,5.25];
ybs[2844]=['ρ Gem',1.9659018,0.553908,4.18];
ybs[2845]=['',1.9553643,-0.3126092,5.63];
ybs[2846]=['γ CMi',1.9608302,0.1549545,4.32];
ybs[2847]=['',1.9545277,-0.403743,5.61];
ybs[2848]=['',1.9527722,-0.5966812,5.9];
ybs[2849]=['64 Gem',1.9667255,0.4899156,5.05];
ybs[2850]=['',1.9637852,0.262878,6.22];
ybs[2851]=['',1.9587818,-0.2025291,5.79];
ybs[2852]=['',1.9577017,-0.3997981,5.95];
ybs[2853]=['65 Gem',1.9687734,0.486387,5.01];
ybs[2854]=['',1.9501536,-0.891245,5.1];
ybs[2855]=['',1.9586113,-0.5096881,5.54];
ybs[2856]=['6 CMi',1.9680683,0.2087155,4.54];
ybs[2857]=['',1.9654557,-0.034089,5.59];
ybs[2858]=['',1.9657641,-0.1326282,5.86];
ybs[2859]=['',1.965402,-0.1810701,5.75];
ybs[2860]=['',1.9652142,-0.2626204,6.05];
ybs[2861]=['',1.9598864,-0.6607396,6.58];
ybs[2862]=['',1.9622736,-0.5566887,6.38];
ybs[2863]=['',1.9622882,-0.5566644,7.13];
ybs[2864]=['',1.9785612,0.6780098,6.54];
ybs[2865]=['',1.9632823,-0.54985,5.77];
ybs[2866]=['',1.9670429,-0.402692,4.85];
ybs[2867]=['',1.9629507,-0.6782332,5.43];
ybs[2868]=['',1.9720581,-0.0920667,6.24];
ybs[2869]=['',1.9770273,0.2973507,5.42];
ybs[2870]=['σ Pup',1.9632608,-0.7565848,3.25];
ybs[2871]=['',1.9817797,0.3985996,6.54];
ybs[2872]=['δ1 CMi',1.9777491,0.0325532,5.25];
ybs[2873]=['',1.9704074,-0.5412394,4.65];
ybs[2874]=['',1.9700714,-0.6525475,6.65];
ybs[2875]=['',1.9773606,-0.1558596,5.9];
ybs[2876]=['',1.9658503,-0.9197742,5.87];
ybs[2877]=['',1.9732946,-0.6318378,6.68];
ybs[2878]=['68 Gem',1.9848386,0.2753534,5.25];
ybs[2879]=['δ2 CMi',1.9825752,0.0565563,5.59];
ybs[2880]=['',1.959282,-1.1267392,6.39];
ybs[2881]=['',1.9745458,-0.6272152,6.61];
ybs[2882]=['α Gem',1.9898415,0.555678,2.88];
ybs[2883]=['α Gem',1.9898415,0.5556731,1.98];
ybs[2884]=['',1.9679544,-0.9502928,5.96];
ybs[2885]=['',1.986724,0.1835743,6.28];
ybs[2886]=['',2.0010464,0.9722092,5.92];
ybs[2887]=['',1.9774163,-0.6285057,6.3];
ybs[2888]=['',1.9921738,0.5394853,5.33];
ybs[2889]=['',1.9827107,-0.2511219,6.21];
ybs[2890]=['',1.9962764,0.7501391,6.3];
ybs[2891]=['',1.9823275,-0.3396822,5.66];
ybs[2892]=['',1.9813982,-0.432154,5.85];
ybs[2893]=['δ3 CMi',1.9872465,0.0579626,5.81];
ybs[2894]=['',1.9845805,-0.2543642,4.97];
ybs[2895]=['',1.9990791,0.805097,5.65];
ybs[2896]=['',1.9894199,0.0466766,6.55];
ybs[2897]=['υ Gem',1.9953761,0.4685262,4.06];
ybs[2898]=['',1.9853908,-0.390017,4.45];
ybs[2899]=['',1.9809102,-0.7000278,6.26];
ybs[2900]=['',1.9807059,-0.7528674,6.52];
ybs[2901]=['',1.9864626,-0.4105704,5.83];
ybs[2902]=['',1.986499,-0.4105899,5.87];
ybs[2903]=['',1.9838505,-0.6350976,5.54];
ybs[2904]=['',1.9870926,-0.4567018,6.65];
ybs[2905]=['',1.9855864,-0.5849226,6.11];
ybs[2906]=['',2.0052551,0.8503471,5.92];
ybs[2907]=['',2.0020339,0.6976657,6.38];
ybs[2908]=['',1.9874973,-0.4723283,5.77];
ybs[2909]=['',1.9841943,-0.6973631,6.76];
ybs[2910]=['',1.9974342,0.1014059,5.91];
ybs[2911]=['ε Men',1.9388758,-1.381246,5.53];
ybs[2912]=['',1.995622,-0.1459575,6.27];
ybs[2913]=['',1.9944741,-0.253841,5.7];
ybs[2914]=['',1.9909218,-0.496028,4.64];
ybs[2915]=['',1.9944608,-0.3876693,6.34];
ybs[2916]=['70 Gem',2.0072096,0.610796,5.56];
ybs[2917]=['',1.9863369,-0.8992834,6.28];
ybs[2918]=['',2.0053863,0.4242525,6.27];
ybs[2919]=['25 Mon',2.0001409,-0.0726576,5.13];
ybs[2920]=['',1.9969742,-0.3447681,5.74];
ybs[2921]=['23 Lyn',2.0188068,0.9953433,6.06];
ybs[2922]=['ο Gem',2.0098821,0.6026845,4.9];
ybs[2923]=['',2.009553,0.4218395,6.17];
ybs[2924]=['',2.0013862,-0.252953,6.53];
ybs[2925]=['',1.9994142,-0.4158565,6.37];
ybs[2926]=['',1.9906221,-0.917778,4.94];
ybs[2927]=['',2.0147848,0.6683042,5.73];
ybs[2928]=['',2.0129671,0.5577457,6.17];
ybs[2929]=['',1.9992842,-0.6112222,4.53];
ybs[2930]=['74 Gem',2.0105072,0.3075571,5.05];
ybs[2931]=['',2.0196249,0.8391146,5.56];
ybs[2932]=['',1.9956183,-0.8531474,5.72];
ybs[2933]=['',1.9918809,-0.9763125,6.39];
ybs[2934]=['',2.00092,-0.6166117,6.6];
ybs[2935]=['α CMi',2.0092936,0.0902704,0.38];
ybs[2936]=['',2.0037999,-0.4436112,4.7];
ybs[2937]=['',2.0008016,-0.6643172,6.38];
ybs[2938]=['24 Lyn',2.0285094,1.0237297,4.99];
ybs[2939]=['',2.0076419,-0.3269336,6.72];
ybs[2940]=['',2.0060154,-0.4686951,4.5];
ybs[2941]=['',2.0060517,-0.4687291,4.62];
ybs[2942]=['',2.0128496,0.0903652,6.02];
ybs[2943]=['',2.0172525,0.4008124,5.89];
ybs[2944]=['',2.0035259,-0.6988949,6.59];
ybs[2945]=['',2.0160834,0.2394102,6.24];
ybs[2946]=['',2.0051563,-0.6379083,5.8];
ybs[2947]=['',2.0042075,-0.6777729,6.19];
ybs[2948]=['',2.0087483,-0.469772,6.5];
ybs[2949]=['',2.0025013,-0.8491614,5.68];
ybs[2950]=['',2.0144567,-0.1438035,6.01];
ybs[2951]=['',2.0133068,-0.267332,4.94];
ybs[2952]=['',2.0124415,-0.3440765,5.93];
ybs[2953]=['',2.0082224,-0.6695292,4.84];
ybs[2954]=['',2.025376,0.5924631,6.02];
ybs[2955]=['',2.0094253,-0.666584,5.73];
ybs[2956]=['',2.0097096,-0.6687031,5.76];
ybs[2957]=['',2.0207623,0.2343346,5.77];
ybs[2958]=['',2.0192076,0.0623206,5.94];
ybs[2959]=['',2.0216176,0.2470351,5.56];
ybs[2960]=['',2.0104898,-0.6568121,6];
ybs[2961]=['',2.0322903,0.8792711,5.27];
ybs[2962]=['α Mon',2.0172668,-0.1676377,3.93];
ybs[2963]=['',2.005164,-0.9307128,6.06];
ybs[2964]=['',2.0142667,-0.4886809,6.76];
ybs[2965]=['σ Gem',2.027693,0.5031556,4.28];
ybs[2966]=['',2.014773,-0.5535211,6.56];
ybs[2967]=['51 Cam',2.0455191,1.1414295,5.92];
ybs[2968]=['',2.0174278,-0.3907982,6.18];
ybs[2969]=['49 Cam',2.0441281,1.0956121,6.49];
ybs[2970]=['',2.0276735,0.3899854,6.21];
ybs[2971]=['',1.9848039,-1.2972346,7.16];
ybs[2972]=['',1.9848113,-1.2972346,7.26];
ybs[2973]=['',2.0160986,-0.6734767,5.42];
ybs[2974]=['',2.0256486,0.0023509,6.19];
ybs[2975]=['76 Gem',2.0310619,0.4490534,5.31];
ybs[2976]=['',2.0161579,-0.779918,6.41];
ybs[2977]=['κ Gem',2.0324506,0.4248584,3.57];
ybs[2978]=['',2.01917,-0.6734003,6.54];
ybs[2979]=['',2.0310785,0.2234738,6.43];
ybs[2980]=['',2.0234089,-0.4608658,5.64];
ybs[2981]=['',2.0302333,0.0410107,6.47];
ybs[2982]=['β Gem',2.0363916,0.4881723,1.14];
ybs[2983]=['79 Gem',2.0353769,0.353614,6.33];
ybs[2984]=['',2.0271554,-0.4460862,6.55];
ybs[2985]=['1 Pup',2.0265423,-0.4968257,4.59];
ybs[2986]=['',2.0246901,-0.6301511,5.6];
ybs[2987]=['',2.0241727,-0.6792619,6.89];
ybs[2988]=['3 Pup',2.0276895,-0.5063158,3.96];
ybs[2989]=['',2.0942842,1.3998205,6.56];
ybs[2990]=['',2.023065,-0.7893758,5.06];
ybs[2991]=['',2.0426945,0.6538164,5.18];
ybs[2992]=['',1.9858762,-1.3558575,6.18];
ybs[2993]=['',2.0268226,-0.6677091,6.4];
ybs[2994]=['',2.0265919,-0.7153902,5.17];
ybs[2995]=['81 Gem',2.0395313,0.3220783,4.88];
ybs[2996]=['',2.0311844,-0.4316079,5.62];
ybs[2997]=['',2.0234058,-0.873492,6.57];
ybs[2998]=['',2.0183913,-1.0242455,6.43];
ybs[2999]=['',2.0288812,-0.6303775,5.8];
ybs[3000]=['11 CMi',2.0398888,0.1869598,5.3];
ybs[3001]=['2 Pup',2.0355503,-0.2572965,6.89];
ybs[3002]=['2 Pup',2.0355793,-0.2573789,6.07];
ybs[3003]=['',2.0305703,-0.6631979,5.88];
ybs[3004]=['',2.0215808,-1.0172559,6.21];
ybs[3005]=['π Gem',2.0461903,0.582217,5.14];
ybs[3006]=['',2.0382673,-0.1191829,5.49];
ybs[3007]=['4 Pup',2.037606,-0.2551673,5.04];
ybs[3008]=['',2.0327858,-0.6622375,6.54];
ybs[3009]=['',2.0335601,-0.6636499,3.61];
ybs[3010]=['',2.0351925,-0.5974079,5.37];
ybs[3011]=['',2.0411657,-0.2222117,6.39];
ybs[3012]=['',2.0334169,-0.7645929,6.03];
ybs[3013]=['82 Gem',2.0503395,0.4028856,6.18];
ybs[3014]=['',2.037572,-0.6630517,5.88];
ybs[3015]=['',2.0428033,-0.3940332,5.9];
ybs[3016]=['ζ Vol',2.0137878,-1.2681541,3.95];
ybs[3017]=['',2.0391258,-0.7001577,6.57];
ybs[3018]=['',2.0449443,-0.2800866,6.34];
ybs[3019]=['',2.0454308,-0.2804996,6.43];
ybs[3020]=['',2.0632382,0.943706,6.02];
ybs[3021]=['5 Pup',2.0464054,-0.2138056,5.48];
ybs[3022]=['',2.052036,0.2323585,6.04];
ybs[3023]=['',2.0335805,-0.9909684,6.12];
ybs[3024]=['',2.0415258,-0.6874504,6.31];
ybs[3025]=['',2.0514964,0.0746153,6.53];
ybs[3026]=['ο Pup',2.0464977,-0.4536872,4.5];
ybs[3027]=['',2.0429694,-0.6731367,5.08];
ybs[3028]=['',2.0283795,-1.1541383,6.38];
ybs[3029]=['',2.0429311,-0.8144651,5.23];
ybs[3030]=['',2.0252063,-1.219573,6.18];
ybs[3031]=['',2.0699474,0.9625476,6.38];
ybs[3032]=['',2.0615883,0.5790115,6.03];
ybs[3033]=['',2.0460127,-0.7105121,6.14];
ybs[3034]=['',2.0530488,-0.2340691,6.23];
ybs[3035]=['',2.0506511,-0.4358056,5.33];
ybs[3036]=['6 Pup',2.0538243,-0.3017022,5.18];
ybs[3037]=['ξ Pup',2.0518173,-0.4348915,3.34];
ybs[3038]=['',2.0464497,-0.8226605,4.71];
ybs[3039]=['',2.0562504,-0.161295,5.61];
ybs[3040]=['',2.0540029,-0.3536893,6.56];
ybs[3041]=['',2.0511266,-0.6161189,5.93];
ybs[3042]=['',2.0593503,0.0561771,6.18];
ybs[3043]=['',2.0555203,-0.3417658,6.12];
ybs[3044]=['',2.05273,-0.5820104,5.6];
ybs[3045]=['',2.0649551,0.3362583,5.99];
ybs[3046]=['',2.0594266,-0.1952527,6.16];
ybs[3047]=['',2.0504317,-0.8103734,4.11];
ybs[3048]=['',2.0455239,-0.9866041,6.33];
ybs[3049]=['',2.0515558,-0.7820768,6.32];
ybs[3050]=['',2.050295,-0.8188282,5.84];
ybs[3051]=['ζ CMi',2.0632697,0.0298102,5.14];
ybs[3052]=['',2.0592738,-0.4291222,6.45];
ybs[3053]=['',2.0651532,0.0561659,6.31];
ybs[3054]=['',2.0489901,-0.9855539,5.59];
ybs[3055]=['8 Pup',2.0626911,-0.2247696,6.36];
ybs[3056]=['9 Pup',2.0630456,-0.2435957,5.17];
ybs[3057]=['25 Lyn',2.0774574,0.8259896,6.25];
ybs[3058]=['26 Lyn',2.0784438,0.829105,5.45];
ybs[3059]=['φ Gem',2.0720181,0.4661072,4.97];
ybs[3060]=['',2.0625282,-0.3705822,5.63];
ybs[3061]=['',2.0569687,-0.7790814,6.45];
ybs[3062]=['',2.0489623,-1.0531519,5.78];
ybs[3063]=['',2.0551797,-0.8825765,5.91];
ybs[3064]=['',2.0678198,-0.0957752,5.76];
ybs[3065]=['10 Pup',2.065383,-0.2601518,5.69];
ybs[3066]=['',2.0598386,-0.7531833,6.32];
ybs[3067]=['',2.1068063,1.2890077,5.41];
ybs[3068]=['',2.0520465,-1.0491,6.72];
ybs[3069]=['',2.0867991,0.9851174,6.72];
ybs[3070]=['',2.0612929,-0.7495693,6.04];
ybs[3071]=['',2.0643258,-0.6067536,5.01];
ybs[3072]=['',2.0638131,-0.7092134,3.73];
ybs[3073]=['',2.0500075,-1.1563426,5.79];
ybs[3074]=['',2.1302977,1.3860364,5.42];
ybs[3075]=['',2.0819339,0.6170062,6.23];
ybs[3076]=['',2.0657812,-0.6793235,4.49];
ybs[3077]=['',2.0677202,-0.6357084,5.43];
ybs[3078]=['85 Gem',2.0812017,0.3459767,5.35];
ybs[3079]=['',2.0801908,0.1536237,5.86];
ybs[3080]=['',2.063989,-0.9499198,5.7];
ybs[3081]=['',2.0668883,-0.8669492,4.63];
ybs[3082]=['',2.0680615,-0.8405968,4.24];
ybs[3083]=['',2.0726576,-0.6272287,5.49];
ybs[3084]=['',2.0748128,-0.6092462,6.15];
ybs[3085]=['',2.0838603,0.0772244,6.17];
ybs[3086]=['',2.0937366,0.7664668,6.34];
ybs[3087]=['1 Cnc',2.0868307,0.2745188,5.78];
ybs[3088]=['',2.0774591,-0.5406694,6.44];
ybs[3089]=['',2.0877814,0.1497452,6.05];
ybs[3090]=['',2.0875514,0.0185936,6.35];
ybs[3091]=['',2.0825158,-0.5296446,6.33];
ybs[3092]=['',2.0751401,-0.9188016,6.38];
ybs[3093]=['',2.0791539,-0.7663007,6.02];
ybs[3094]=['11 Pup',2.0849132,-0.4004024,4.2];
ybs[3095]=['',2.0913613,0.1248189,6.41];
ybs[3096]=['',2.0935525,0.287218,5.99];
ybs[3097]=['',2.0741271,-1.0011796,5.63];
ybs[3098]=['',2.1084303,1.0294614,5.77];
ybs[3099]=['',2.0820953,-0.7120506,6.78];
ybs[3100]=['',2.1902842,1.4658354,6.49];
ybs[3101]=['53 Cam',2.1101695,1.0517452,6.01];
ybs[3102]=['14 CMi',2.092271,0.0377444,5.29];
ybs[3103]=['',2.0844416,-0.7411974,6.09];
ybs[3104]=['',2.1140852,1.1000111,6.4];
ybs[3105]=['',2.0881385,-0.5305183,4.79];
ybs[3106]=['',2.0843735,-0.7602941,5.35];
ybs[3107]=['',2.0980603,0.2300255,6.02];
ybs[3108]=['',2.0858269,-0.7709337,5.09];
ybs[3109]=['χ Car',2.0828375,-0.925783,3.47];
ybs[3110]=['',2.0856822,-0.8369168,6.22];
ybs[3111]=['',2.1136732,0.9984911,6.49];
ybs[3112]=['',2.0799427,-1.0574489,5.74];
ybs[3113]=['',2.0881612,-0.7965607,5.17];
ybs[3114]=['27 Mon',2.0981478,-0.0653189,4.93];
ybs[3115]=['12 Pup',2.0946601,-0.4079355,5.11];
ybs[3116]=['ω1 Cnc',2.1043858,0.4420812,5.83];
ybs[3117]=['',2.1035961,0.3447513,6.25];
ybs[3118]=['',2.0824666,-1.0330189,6.25];
ybs[3119]=['',2.1046628,0.410495,6.34];
ybs[3120]=['3 Cnc',2.1034534,0.3009874,5.55];
ybs[3121]=['',2.0895402,-0.8605686,4.41];
ybs[3122]=['',2.1091327,0.6169595,6.34];
ybs[3123]=['',2.0982147,-0.3222222,4.61];
ybs[3124]=['ω2 Cnc',2.1078552,0.4367853,6.31];
ybs[3125]=['',2.0898699,-0.8990297,6.44];
ybs[3126]=['5 Cnc',2.1065487,0.286088,5.99];
ybs[3127]=['',2.1025378,-0.0513985,6.51];
ybs[3128]=['',2.1049536,0.084059,5.65];
ybs[3129]=['',2.0932817,-0.7902581,5.99];
ybs[3130]=['',2.0864067,-1.0535681,5.6];
ybs[3131]=['',2.083452,-1.1058113,6.14];
ybs[3132]=['',2.0955771,-0.686958,5.24];
ybs[3133]=['28 Mon',2.1047105,-0.0254116,4.68];
ybs[3134]=['',2.0936983,-0.8733466,6.32];
ybs[3135]=['',2.0937786,-0.8732934,6.34];
ybs[3136]=['',2.1077712,0.1544634,6.22];
ybs[3137]=['',2.1093853,0.0396271,4.39];
ybs[3138]=['',2.0989296,-0.7944721,6.61];
ybs[3139]=['',2.0909689,-1.0626719,5.81];
ybs[3140]=['',2.0983346,-0.8559844,6.02];
ybs[3141]=['χ Gem',2.1157643,0.483972,4.94];
ybs[3142]=['',2.1098345,-0.1117231,6.33];
ybs[3143]=['',2.0993624,-0.8540664,6.12];
ybs[3144]=['',2.0946746,-1.0519111,6.33];
ybs[3145]=['',2.0944336,-1.0585333,5.17];
ybs[3146]=['',2.1050742,-0.6518317,5.95];
ybs[3147]=['',2.1071825,-0.6477681,6.34];
ybs[3148]=['',2.1004394,-0.9462223,5.87];
ybs[3149]=['',2.1028077,-0.952573,6.1];
ybs[3150]=['',2.1208031,0.3277213,6.15];
ybs[3151]=['',2.0970718,-1.1105592,4.82];
ybs[3152]=['',2.1116197,-0.5677237,5.82];
ybs[3153]=['',2.1033399,-0.9689802,6.28];
ybs[3154]=['',2.1097589,-0.7221144,5.52];
ybs[3155]=['8 Cnc',2.1220044,0.2278134,5.12];
ybs[3156]=['',2.1249036,0.4793395,6.21];
ybs[3157]=['ζ Pup',2.1134948,-0.6993158,2.25];
ybs[3158]=['',2.1129152,-0.7507196,6.29];
ybs[3159]=['28 Lyn',2.132433,0.7538762,6.26];
ybs[3160]=['14 Pup',2.119227,-0.3454552,6.13];
ybs[3161]=['μ1 Cnc',2.1277189,0.3939148,5.99];
ybs[3162]=['',2.1168571,-0.5714181,5.31];
ybs[3163]=['',2.0899192,-1.2794475,6.34];
ybs[3164]=['',2.1248307,-0.011157,6.41];
ybs[3165]=['27 Lyn',2.1386464,0.8977919,4.84];
ybs[3166]=['',2.1273029,-0.1625061,6.23];
ybs[3167]=['',2.1463038,1.0154379,5.93];
ybs[3168]=['μ2 Cnc',2.134032,0.3755091,5.3];
ybs[3169]=['',2.1232706,-0.5870363,6.14];
ybs[3170]=['',2.1176864,-0.8841062,5.95];
ybs[3171]=['',2.1207255,-0.8210761,6.19];
ybs[3172]=['',2.1190437,-0.9280475,5.53];
ybs[3173]=['',2.1420594,0.739377,6.27];
ybs[3174]=['',2.1599505,1.1938937,5.32];
ybs[3175]=['',2.1305823,-0.3599045,5.38];
ybs[3176]=['12 Cnc',2.137864,0.2369082,6.27];
ybs[3177]=['ρ Pup',2.1315023,-0.4253464,2.81];
ybs[3178]=['',2.1163512,-1.0978301,6.3];
ybs[3179]=['',2.1266663,-0.7911982,5.05];
ybs[3180]=['ζ Mon',2.1368249,-0.0532466,4.34];
ybs[3181]=['',2.1380999,-0.199086,6.32];
ybs[3182]=['',2.1368098,-0.3565708,6.36];
ybs[3183]=['ψ Cnc',2.1459138,0.4440012,5.73];
ybs[3184]=['16 Pup',2.1381674,-0.3370595,4.4];
ybs[3185]=['',2.1504626,0.6747987,6.58];
ybs[3186]=['',2.1402307,-0.2847712,5.68];
ybs[3187]=['',2.1356381,-0.6588312,6.37];
ybs[3188]=['',2.1381027,-0.5303987,6.65];
ybs[3189]=['',2.2196723,1.4373841,6.32];
ybs[3190]=['',2.1478141,0.2541444,6.23];
ybs[3191]=['',2.1381224,-0.6199779,6.2];
ybs[3192]=['',2.1625407,0.9840646,5.85];
ybs[3193]=['',2.1489469,0.1702209,6.07];
ybs[3194]=['18 Pup',2.1455019,-0.2420251,5.54];
ybs[3195]=['',2.1372917,-0.8508745,5.7];
ybs[3196]=['',2.1395169,-0.7712621,5.21];
ybs[3197]=['',2.1404665,-0.7453942,6.26];
ybs[3198]=['γ1 Vel',2.1388229,-0.827514,4.27];
ybs[3199]=['γ2 Vel',2.1390201,-0.8273544,1.78];
ybs[3200]=['ζ1 Cnc',2.153291,0.3068142,5.63];
ybs[3201]=['ζ1 Cnc',2.153291,0.3068142,6.02];
ybs[3202]=['ζ2 Cnc',2.1533346,0.3068142,6.2];
ybs[3203]=['19 Pup',2.1481865,-0.2268068,4.72];
ybs[3204]=['',2.1495719,-0.1368472,5.36];
ybs[3205]=['',2.1397929,-0.8378424,5.23];
ybs[3206]=['',2.1538487,0.2432153,6.54];
ybs[3207]=['15 Cnc',2.1578343,0.5164011,5.64];
ybs[3208]=['',2.1916884,1.3209469,5.54];
ybs[3209]=['',2.1323472,-1.1147025,6.28];
ybs[3210]=['',2.1383929,-0.9800508,5.66];
ybs[3211]=['',2.1461394,-0.6520584,6.44];
ybs[3212]=['',2.1353672,-1.0710987,4.76];
ybs[3213]=['',2.1717253,1.0526108,6.45];
ybs[3214]=['',2.1566754,0.2870228,6.01];
ybs[3215]=['ε Vol',2.1292843,-1.1987547,4.35];
ybs[3216]=['',2.1599595,0.4026209,6.56];
ybs[3217]=['',2.1474662,-0.692664,4.45];
ybs[3218]=['',2.1475949,-0.7514578,4.75];
ybs[3219]=['',2.1461557,-0.8470074,5.82];
ybs[3220]=['',2.1618989,0.3072882,6.47];
ybs[3221]=['20 Pup',2.1570875,-0.2767639,4.99];
ybs[3222]=['',2.1540793,-0.523243,6.52];
ybs[3223]=['',2.1624607,0.226522,6.38];
ybs[3224]=['',2.1498444,-0.8152879,5.76];
ybs[3225]=['',2.1541067,-0.6631073,6.43];
ybs[3226]=['',2.15212,-0.8086597,6.03];
ybs[3227]=['29 Lyn',2.1803246,1.0384677,5.64];
ybs[3228]=['',2.1952269,1.2624757,5.98];
ybs[3229]=['',2.156981,-0.6277744,4.78];
ybs[3230]=['',2.1579279,-0.5871002,6.37];
ybs[3231]=['',2.1601684,-0.5621751,6.06];
ybs[3232]=['',2.1590566,-0.635157,5.08];
ybs[3233]=['',2.1590848,-0.6354819,6.11];
ybs[3234]=['',2.1601742,-0.6206388,5.78];
ybs[3235]=['',2.1591756,-0.7054166,4.44];
ybs[3236]=['',2.1568254,-0.8213705,5.13];
ybs[3237]=['',2.1870214,1.0897008,5.71];
ybs[3238]=['',2.181597,0.9437372,6.27];
ybs[3239]=['',2.1564386,-0.8772932,5.51];
ybs[3240]=['',2.1720791,0.2034322,7.13];
ybs[3241]=['β Cnc',2.1717781,0.1590868,3.52];
ybs[3242]=['',2.1603763,-0.8011747,5.83];
ybs[3243]=['',2.1676083,-0.5409827,6.21];
ybs[3244]=['',2.1761884,0.1535034,6.29];
ybs[3245]=['',2.1678478,-0.6278475,6.16];
ybs[3246]=['30 Lyn',2.1914134,1.0065469,5.89];
ybs[3247]=['',2.1724758,-0.3733426,6.6];
ybs[3248]=['',2.1643718,-0.8817292,6.44];
ybs[3249]=['21 Pup',2.1747555,-0.2854644,6.16];
ybs[3250]=['',2.1912232,0.9337861,6.49];
ybs[3251]=['',2.1793009,-0.2217146,5.98];
ybs[3252]=['',2.1624774,-1.0993065,5.16];
ybs[3253]=['',2.1768025,-0.5248987,6.45];
ybs[3254]=['χ Cnc',2.1878945,0.47378,5.14];
ybs[3255]=['',2.2017685,1.05693,6.41];
ybs[3256]=['',2.1888955,0.3608552,5.83];
ybs[3257]=['',2.1831055,-0.1786797,6.32];
ybs[3258]=['',2.1779502,-0.6199923,5.58];
ybs[3259]=['',2.1775075,-0.6535456,6.7];
ybs[3260]=['λ Cnc',2.1898215,0.4180035,5.98];
ybs[3261]=['',2.1860811,0.0676444,6.05];
ybs[3262]=['',2.1790495,-0.641074,4.45];
ybs[3263]=['',2.1876164,-0.017133,6.18];
ybs[3264]=['',2.1877691,-0.0942721,6.13];
ybs[3265]=['',2.1832311,-0.6049676,6.43];
ybs[3266]=['',2.1746399,-1.0338976,6.42];
ybs[3267]=['31 Lyn',2.2007235,0.7524918,4.25];
ybs[3268]=['',2.1879293,-0.4013733,6.13];
ybs[3269]=['',2.2056682,0.9275693,5.51];
ybs[3270]=['',2.1924747,-0.0292329,6.5];
ybs[3271]=['',2.1919462,-0.351716,5.58];
ybs[3272]=['',2.1753654,-1.1464104,5.07];
ybs[3273]=['',2.1944626,-0.3082132,5.75];
ybs[3274]=['',2.1915687,-0.5781771,4.83];
ybs[3275]=['',2.1912649,-0.6380415,5.2];
ybs[3276]=['20 Cnc',2.2019534,0.3186725,5.95];
ybs[3277]=['',2.1974229,-0.1091246,6.15];
ybs[3278]=['',2.1913402,-0.6927821,6.16];
ybs[3279]=['',2.2088424,0.7318292,6.02];
ybs[3280]=['',2.1991112,-0.1329368,5.96];
ybs[3281]=['22 Pup',2.1984111,-0.2291276,6.11];
ybs[3282]=['21 Cnc',2.2041278,0.1842732,6.08];
ybs[3283]=['',2.1981694,-0.46114,5.9];
ybs[3284]=['',2.2100964,0.6097651,6.06];
ybs[3285]=['',2.189075,-1.013086,5.97];
ybs[3286]=['',2.1957097,-0.8475913,4.82];
ybs[3287]=['',2.2066471,-0.0836202,6.01];
ybs[3288]=['',2.1996348,-0.6694966,6.32];
ybs[3289]=['1 Hya',2.2065751,-0.0667631,5.61];
ybs[3290]=['',2.187936,-1.1201267,6.12];
ybs[3291]=['25 Cnc',2.2126755,0.2962068,6.14];
ybs[3292]=['',2.1971517,-0.9110126,5.85];
ybs[3293]=['κ1 Vol',2.1805215,-1.2494245,5.37];
ybs[3294]=['κ2 Vol',2.1813786,-1.2492563,5.65];
ybs[3295]=['',2.2334497,1.1732259,5.88];
ybs[3296]=['φ1 Cnc',2.215813,0.4855262,5.57];
ybs[3297]=['',2.2111549,0.0353889,5.73];
ybs[3298]=['',2.2127269,0.1307201,5.13];
ybs[3299]=['ε Car',2.1946397,-1.0399159,1.86];
ybs[3300]=['',2.2073992,-0.4054026,5.68];
ybs[3301]=['',2.2216943,0.7954774,6.32];
ybs[3302]=['φ2 Cnc',2.2171639,0.4687831,6.32];
ybs[3303]=['φ2 Cnc',2.2171784,0.4687977,6.3];
ybs[3304]=['24 Cnc',2.2165649,0.4268914,7.02];
ybs[3305]=['24 Cnc',2.2165868,0.4269107,7.81];
ybs[3306]=['',2.2112687,-0.0694814,3.9];
ybs[3307]=['',2.2079924,-0.4209808,5.28];
ybs[3308]=['',2.209221,-0.3686181,6.01];
ybs[3309]=['',2.2108221,-0.3056774,6.44];
ybs[3310]=['α Cha',2.172493,-1.3437424,4.07];
ybs[3311]=['27 Cnc',2.2164574,0.2195509,5.5];
ybs[3312]=['',2.2120841,-0.2618767,5.98];
ybs[3313]=['2 Hya',2.214728,-0.0709032,5.59];
ybs[3314]=['',2.2066751,-0.747763,5.98];
ybs[3315]=['ο UMa',2.2345587,1.0583901,3.36];
ybs[3316]=['',2.21553,-0.2200772,5.54];
ybs[3317]=['',2.2182899,-0.1131853,6.59];
ybs[3318]=['',2.2106925,-0.7370169,5.47];
ybs[3319]=['',2.2127283,-0.6830217,6.53];
ybs[3320]=['',2.2127719,-0.6830364,7.25];
ybs[3321]=['28 Cnc',2.225054,0.4200795,6.1];
ybs[3322]=['',2.2085326,-0.9041241,5.17];
ybs[3323]=['',2.2155874,-0.5112133,6.73];
ybs[3324]=['',2.2474602,1.2085013,6.31];
ybs[3325]=['29 Cnc',2.2247515,0.2467006,5.95];
ybs[3326]=['η Vol',2.1897549,-1.282341,5.29];
ybs[3327]=['',2.2189923,-0.3651106,6.56];
ybs[3328]=['',2.2173494,-0.5541127,6.33];
ybs[3329]=['',2.2236457,-0.0452575,6.39];
ybs[3330]=['',2.2227613,-0.1551926,6.43];
ybs[3331]=['',2.2202703,-0.4574167,6.62];
ybs[3332]=['θ Cha',2.1814085,-1.353615,4.35];
ybs[3333]=['',2.2123922,-0.9229707,6.05];
ybs[3334]=['',2.2250017,-0.1714668,6];
ybs[3335]=['',2.2203212,-0.614172,5.75];
ybs[3336]=['',2.2234682,-0.4040005,6.51];
ybs[3337]=['',2.2248092,-0.3669776,6.67];
ybs[3338]=['',2.208555,-1.1287977,5.97];
ybs[3339]=['β Vol',2.2077426,-1.1556067,3.77];
ybs[3340]=['',2.2374457,0.6490651,6.18];
ybs[3341]=['',2.2167289,-0.9614485,6.53];
ybs[3342]=['',2.2175626,-0.9278862,5.09];
ybs[3343]=['',2.2437709,0.9256703,6.24];
ybs[3344]=['',2.2662056,1.3027817,6.31];
ybs[3345]=['',2.2270722,-0.4783725,6.7];
ybs[3346]=['2 UMa',2.2540787,1.1356218,5.47];
ybs[3347]=['υ1 Cnc',2.2376759,0.4189474,5.75];
ybs[3348]=['',2.2248265,-0.7720743,5.79];
ybs[3349]=['θ Cnc',2.2378472,0.3144598,5.35];
ybs[3350]=['',2.2243773,-0.8378484,5.33];
ybs[3351]=['',2.2262492,-0.7819282,4.99];
ybs[3352]=['',2.2443903,0.6621529,5.9];
ybs[3353]=['',2.2389567,0.1699445,6.83];
ybs[3354]=['',2.2313176,-0.5626261,5.65];
ybs[3355]=['',2.2274606,-0.8099769,5.99];
ybs[3356]=['',2.2311875,-0.6422422,6.69];
ybs[3357]=['32 Lyn',2.2462628,0.6345734,6.24];
ybs[3358]=['η Cnc',2.2427816,0.3554085,5.33];
ybs[3359]=['',2.236323,-0.343038,5.42];
ybs[3360]=['',2.2261039,-0.9645967,6.36];
ybs[3361]=['υ2 Cnc',2.244188,0.4189992,6.36];
ybs[3362]=['',2.213623,-1.2246699,5.53];
ybs[3363]=['',2.2314623,-0.7821505,6.3];
ybs[3364]=['34 Cnc',2.2422581,0.1743312,6.46];
ybs[3365]=['',2.2350778,-0.683143,6.31];
ybs[3366]=['',2.2410127,-0.2636674,6.38];
ybs[3367]=['',2.2335652,-0.8367734,6.39];
ybs[3368]=['',2.2470976,0.2300184,6.28];
ybs[3369]=['33 Lyn',2.2522192,0.6342674,5.78];
ybs[3370]=['',2.246721,0.0816562,5.87];
ybs[3371]=['',2.2783614,1.283669,6.15];
ybs[3372]=['',2.2489994,0.1461472,6.03];
ybs[3373]=['',2.2429794,-0.4308202,6.19];
ybs[3374]=['',2.2344315,-0.9507017,6.34];
ybs[3375]=['',2.2478308,-0.0389191,5.81];
ybs[3376]=['',2.2417457,-0.5511492,6.38];
ybs[3377]=['',2.2421288,-0.605832,6.36];
ybs[3378]=['',2.2370981,-0.9300774,5.69];
ybs[3379]=['35 Cnc',2.2541573,0.3405344,6.58];
ybs[3380]=['',2.2435059,-0.6710613,6.49];
ybs[3381]=['',2.244822,-0.6794024,5.96];
ybs[3382]=['',2.243791,-0.8211605,6.24];
ybs[3383]=['π1 UMa',2.2740333,1.133421,5.64];
ybs[3384]=['',2.2540342,0.0465091,6.33];
ybs[3385]=['',2.1944679,-1.4135008,5.69];
ybs[3386]=['',2.2575292,0.2658914,6.32];
ybs[3387]=['',2.2560476,0.1141616,5.99];
ybs[3388]=['',2.2560696,0.1142052,7.25];
ybs[3389]=['',2.2489931,-0.5703214,6.43];
ybs[3390]=['3 Hya',2.2539545,-0.1406921,5.72];
ybs[3391]=['',2.2485947,-0.6578104,6.3];
ybs[3392]=['',2.269106,0.9306306,5.66];
ybs[3393]=['',2.2732264,1.0447348,6.48];
ybs[3394]=['',2.2533761,-0.4698849,5.96];
ybs[3395]=['π2 UMa',2.2783534,1.1213175,4.6];
ybs[3396]=['',2.2516314,-0.6989812,6.47];
ybs[3397]=['',2.271811,0.9223114,6.42];
ybs[3398]=['36 Cnc',2.2615776,0.1671327,5.88];
ybs[3399]=['',2.2489203,-0.8730592,5.01];
ybs[3400]=['',2.2730715,0.9185859,5.91];
ybs[3401]=['',2.2676964,0.5711038,5.94];
ybs[3402]=['δ Hya',2.2639022,0.0981543,4.16];
ybs[3403]=['',2.2626937,-0.0874985,6.19];
ybs[3404]=['37 Cnc',2.2658939,0.1657147,6.53];
ybs[3405]=['',2.2538346,-0.8909717,5.8];
ybs[3406]=['',2.2508404,-1.0138238,4.86];
ybs[3407]=['',2.2505112,-1.0175903,5.26];
ybs[3408]=['',2.2665117,-0.1176798,6.51];
ybs[3409]=['',2.2363569,-1.2816661,6.12];
ybs[3410]=['σ Hya',2.2686316,0.0569177,4.44];
ybs[3411]=['',2.2619007,-0.5903661,6.48];
ybs[3412]=['η Pyx',2.2638333,-0.4596294,5.27];
ybs[3413]=['',2.2608907,-0.7020949,6.55];
ybs[3414]=['34 Lyn',2.2801135,0.798534,5.37];
ybs[3415]=['',2.2763292,0.5560795,6.1];
ybs[3416]=['',2.2716201,0.1385215,6.45];
ybs[3417]=['',2.2675598,-0.345874,6.33];
ybs[3418]=['',2.2621096,-0.7516934,4.14];
ybs[3419]=['39 Cnc',2.2750319,0.3477908,6.39];
ybs[3420]=['',2.276162,0.3418936,6.44];
ybs[3421]=['ε Cnc',2.2765141,0.3397113,6.3];
ybs[3422]=['',2.2694766,-0.3969281,5.05];
ybs[3423]=['6 Hya',2.2736978,-0.2191441,4.98];
ybs[3424]=['',2.2589355,-1.0983895,5.47];
ybs[3425]=['ζ Pyx',2.2717363,-0.5173453,4.89];
ybs[3426]=['',2.2699568,-0.6403156,6.13];
ybs[3427]=['',2.2662577,-0.9280035,6.47];
ybs[3428]=['',2.2888264,0.8171461,6.22];
ybs[3429]=['',2.278173,-0.159403,6.63];
ybs[3430]=['β Pyx',2.2732251,-0.617656,3.97];
ybs[3431]=['',2.2739466,-0.704153,5.2];
ybs[3432]=['',2.2690632,-0.9341023,5.48];
ybs[3433]=['9 Hya',2.2809962,-0.2796853,4.88];
ybs[3434]=['',2.2715573,-0.9273918,5.19];
ybs[3435]=['',2.2676075,-1.0541354,6.36];
ybs[3436]=['',2.2748445,-0.7901509,5.71];
ybs[3437]=['',2.2749241,-0.8155893,3.84];
ybs[3438]=['',2.2830464,-0.2102729,6.45];
ybs[3439]=['',2.2730252,-0.925072,3.62];
ybs[3440]=['',2.2730031,-0.926701,5.61];
ybs[3441]=['γ Cnc',2.2889409,0.373264,4.66];
ybs[3442]=['45 Cnc',2.2883145,0.2198893,5.64];
ybs[3443]=['',2.2933812,0.6429009,6.33];
ybs[3444]=['',2.2774678,-0.8272534,4.77];
ybs[3445]=['',2.2767947,-0.8552745,5.9];
ybs[3446]=['η Hya',2.2881283,0.057884,4.3];
ybs[3447]=['',2.2745026,-1.0057672,6.34];
ybs[3448]=['',2.2807742,-0.7939909,5.23];
ybs[3449]=['',2.2737818,-1.0444398,4.33];
ybs[3450]=['',2.2915148,0.0742167,6.37];
ybs[3451]=['',2.2897773,-0.1276862,4.62];
ybs[3452]=['θ Vol',2.2652492,-1.2298828,5.2];
ybs[3453]=['δ Cnc',2.2949389,0.3154062,3.94];
ybs[3454]=['',2.2819989,-0.8409132,5.51];
ybs[3455]=['',2.2856307,-0.6287594,6.42];
ybs[3456]=['46 Cnc',2.298303,0.5343283,6.13];
ybs[3457]=['49 Cnc',2.2949762,0.1745142,5.66];
ybs[3458]=['',2.2818487,-0.9281943,5.52];
ybs[3459]=['',2.2823285,-0.9284375,4.86];
ybs[3460]=['α Pyx',2.2885507,-0.5806464,3.68];
ybs[3461]=['10 Hya',2.2960328,0.0976985,6.13];
ybs[3462]=['',2.316198,1.1627994,6.2];
ybs[3463]=['',2.2817646,-0.9748673,6.29];
ybs[3464]=['',2.2972134,-0.0468411,6.41];
ybs[3465]=['',2.2947905,-0.3708919,6.11];
ybs[3466]=['ι Cnc',2.3039412,0.5005906,6.57];
ybs[3467]=['ι Cnc',2.3040718,0.5004982,4.02];
ybs[3468]=['',2.2880237,-0.871006,5.16];
ybs[3469]=['',2.2916409,-0.7458085,4.07];
ybs[3470]=['',2.3002619,-0.0372129,5.7];
ybs[3471]=['',2.2939416,-0.6497849,5.76];
ybs[3472]=['',2.3003225,-0.193551,6.25];
ybs[3473]=['50 Cnc',2.3045601,0.2098997,5.87];
ybs[3474]=['ε Hya',2.3037106,0.1105723,3.38];
ybs[3475]=['',2.2985839,-0.4445464,6.1];
ybs[3476]=['12 Hya',2.3013815,-0.2379084,4.32];
ybs[3477]=['δ Vel',2.2921563,-0.9562821,1.96];
ybs[3478]=['',2.3055388,-0.0345744,5.29];
ybs[3479]=['',2.2985655,-0.8050303,3.91];
ybs[3480]=['',2.3004432,-0.7192308,6.21];
ybs[3481]=['',2.293455,-1.0263886,6.21];
ybs[3482]=['',2.3025909,-0.6057392,6.37];
ybs[3483]=['',2.2868623,-1.1919526,6.32];
ybs[3484]=['ρ Hya',2.3109214,0.1004182,4.36];
ybs[3485]=['',2.3090352,-0.1159368,6.09];
ybs[3486]=['',2.300678,-0.8027841,5.46];
ybs[3487]=['',2.289921,-1.1503117,6.05];
ybs[3488]=['',2.3041813,-0.8070271,5.75];
ybs[3489]=['',2.3059971,-0.7299154,6.36];
ybs[3490]=['',2.3007368,-0.9922744,4.49];
ybs[3491]=['',2.3209859,0.579452,6.25];
ybs[3492]=['14 Hya',2.3147095,-0.0615692,5.31];
ybs[3493]=['',2.3080205,-0.7426015,6.43];
ybs[3494]=['η Cha',2.2711511,-1.3795822,5.47];
ybs[3495]=['',2.3067361,-0.9238815,6.3];
ybs[3496]=['',2.3214181,0.3271975,6.16];
ybs[3497]=['5 UMa',2.3352952,1.0799378,5.73];
ybs[3498]=['',2.3337529,1.0292188,6.25];
ybs[3499]=['',2.3158739,-0.3688464,6.47];
ybs[3500]=['35 Lyn',2.3276097,0.7616784,5.15];
ybs[3501]=['',2.328771,0.7893546,5.99];
ybs[3502]=['54 Cnc',2.3225037,0.2664291,6.38];
ybs[3503]=['',2.3284792,0.7315846,5.99];
ybs[3504]=['',2.3159364,-0.5736078,5.21];
ybs[3505]=['',2.3168426,-0.5157079,5.87];
ybs[3506]=['',2.3147107,-0.705204,5.48];
ybs[3507]=['',2.3182857,-0.5009622,6.17];
ybs[3508]=['',2.315728,-0.6846301,6.39];
ybs[3509]=['γ Pyx',2.3190681,-0.4851148,4.01];
ybs[3510]=['σ1 Cnc',2.3298465,0.5652812,5.66];
ybs[3511]=['',2.315068,-0.792253,4.93];
ybs[3512]=['53 Cnc',2.3292534,0.4917165,6.23];
ybs[3513]=['ρ1 Cnc',2.3297792,0.4929664,5.95];
ybs[3514]=['15 Hya',2.3242628,-0.1267581,5.54];
ybs[3515]=['',2.2793357,-1.3814526,6.05];
ybs[3516]=['',2.3176665,-0.7360915,6];
ybs[3517]=['',2.3282272,0.0917026,6.33];
ybs[3518]=['',2.3183327,-0.8135708,5.1];
ybs[3519]=['',2.3358615,0.6187516,6.14];
ybs[3520]=['',2.3281766,-0.2324684,6.13];
ybs[3521]=['',2.3225169,-0.7433326,6.55];
ybs[3522]=['6 UMa',2.349784,1.1260211,5.58];
ybs[3523]=['57 Cnc',2.3370472,0.5322007,5.39];
ybs[3524]=['',2.3272008,-0.5688893,6.5];
ybs[3525]=['',2.327945,-0.6393389,6.42];
ybs[3526]=['',2.3285327,-0.6773638,5.82];
ybs[3527]=['',2.3220899,-1.0073866,5.59];
ybs[3528]=['',2.3163674,-1.1672408,5.35];
ybs[3529]=['',2.3361965,-0.0963597,6];
ybs[3530]=['',2.3273331,-0.8455248,5.91];
ybs[3531]=['ρ2 Cnc',2.3431276,0.4859058,5.22];
ybs[3532]=['',2.3415614,0.2992257,6.64];
ybs[3533]=['',2.3272361,-0.9113237,6.39];
ybs[3534]=['',2.2908984,-1.3890592,5.79];
ybs[3535]=['',2.3117427,-1.2677267,6.11];
ybs[3536]=['',2.3490035,0.7948981,5.74];
ybs[3537]=['',2.347311,0.7001242,5.89];
ybs[3538]=['ζ Hya',2.3412884,0.1022509,3.11];
ybs[3539]=['',2.3330174,-0.7074488,6.47];
ybs[3540]=['',2.3285394,-0.9902198,6.03];
ybs[3541]=['60 Cnc',2.3437745,0.2013915,5.41];
ybs[3542]=['',2.332637,-0.8309015,5.33];
ybs[3543]=['17 Hya',2.3413336,-0.1406266,6.91];
ybs[3544]=['17 Hya',2.3413407,-0.1406411,6.67];
ybs[3545]=['',2.3397911,-0.3198891,5.75];
ybs[3546]=['σ2 Cnc',2.3488981,0.5728628,5.45];
ybs[3547]=['δ Pyx',2.3408738,-0.4846598,4.89];
ybs[3548]=['',2.346571,0.0724169,6.14];
ybs[3549]=['',2.3492127,0.2976865,6.17];
ybs[3550]=['',2.3427803,-0.4172299,6.39];
ybs[3551]=['',2.331434,-1.0548843,5.78];
ybs[3552]=['ο1 Cnc',2.349646,0.2659014,5.2];
ybs[3553]=['',2.3392339,-0.7876422,6.26];
ybs[3554]=['61 Cnc',2.3533097,0.5261393,6.29];
ybs[3555]=['',2.3457815,-0.293161,5.96];
ybs[3556]=['ο2 Cnc',2.3511292,0.2704127,5.67];
ybs[3557]=['',2.3556167,0.6233312,6.51];
ybs[3558]=['',2.3514468,0.1623132,6.19];
ybs[3559]=['',2.3364198,-1.0179929,6.38];
ybs[3560]=['ι UMa',2.3594815,0.8369395,3.14];
ybs[3561]=['',2.3380237,-0.9608453,5.71];
ybs[3562]=['',2.3368039,-1.059964,3.84];
ybs[3563]=['α Cnc',2.3549476,0.2054175,4.25];
ybs[3564]=['',2.3531366,0.02537,6.59];
ybs[3565]=['',2.3431201,-0.9217238,4.69];
ybs[3566]=['σ3 Cnc',2.3602127,0.564264,5.2];
ybs[3567]=['ρ UMa',2.3760203,1.1787917,4.76];
ybs[3568]=['',2.3581537,0.3149659,6.38];
ybs[3569]=['',2.3552421,-0.283111,5.86];
ybs[3570]=['',2.3653897,0.7276918,3.97];
ybs[3571]=['',2.3646485,0.6547672,6.44];
ybs[3572]=['',2.4420708,1.4675781,6.33];
ybs[3573]=['',2.345388,-1.035276,4.92];
ybs[3574]=['',2.3504266,-0.849299,5.87];
ybs[3575]=['',2.359225,-0.3367909,6.18];
ybs[3576]=['',2.357152,-0.5043057,6.25];
ybs[3577]=['',2.3698106,0.6915664,6.36];
ybs[3578]=['66 Cnc',2.368307,0.5613524,5.82];
ybs[3579]=['',2.3546234,-0.8259423,5.18];
ybs[3580]=['67 Cnc',2.3699538,0.4854327,6.07];
ybs[3581]=['',2.3680231,0.096891,6.07];
ybs[3582]=['',2.3602639,-0.7215656,4.45];
ybs[3583]=['',2.3808095,0.9458545,5.75];
ybs[3584]=['',2.361393,-0.755068,6.07];
ybs[3585]=['κ UMa',2.3786753,0.8214639,3.6];
ybs[3586]=['ν Cnc',2.3738655,0.4252129,5.45];
ybs[3587]=['',2.3697934,-0.0099891,5.67];
ybs[3588]=['',2.3656464,-0.4669299,6.2];
ybs[3589]=['',2.3560261,-1.0327473,5.16];
ybs[3590]=['',2.3734047,0.1258068,5.85];
ybs[3591]=['',2.3657253,-0.7322302,5.55];
ybs[3592]=['70 Cnc',2.3801989,0.4853396,6.38];
ybs[3593]=['',2.3691605,-0.6892663,6.27];
ybs[3594]=['',2.3864942,0.8454261,5.95];
ybs[3595]=['',2.3617722,-1.0655734,5.79];
ybs[3596]=['',2.3668559,-0.9124182,5.23];
ybs[3597]=['',2.383648,0.563501,6.46];
ybs[3598]=['',2.374229,-0.4467068,6.74];
ybs[3599]=['',2.3930564,1.0341597,6.45];
ybs[3600]=['σ1 UMa',2.4013101,1.165552,5.14];
ybs[3601]=['',2.3622588,-1.2003144,5.88];
ybs[3602]=['',2.3726295,-0.9361879,6.4];
ybs[3603]=['',2.390905,0.6695239,4.56];
ybs[3604]=['ω Hya',2.3874224,0.0872864,4.97];
ybs[3605]=['',2.3777231,-0.8235876,3.75];
ybs[3606]=['α Vol',2.3683973,-1.1603942,4];
ybs[3607]=['σ2 UMa',2.4100185,1.1701013,4.8];
ybs[3608]=['',2.3943564,0.3994965,6.4];
ybs[3609]=['',2.3917992,0.0239341,6.17];
ybs[3610]=['15 UMa',2.4017999,0.8990626,4.48];
ybs[3611]=['',2.3973773,0.566336,6.5];
ybs[3612]=['τ Cnc',2.3969847,0.5159595,5.43];
ybs[3613]=['',2.3797614,-1.0112965,6.44];
ybs[3614]=['κ Cnc',2.3953066,0.1845915,5.24];
ybs[3615]=['τ UMa',2.4117987,1.106898,4.67];
ybs[3616]=['',2.4008412,0.5897474,5.93];
ybs[3617]=['75 Cnc',2.4003136,0.4631583,5.98];
ybs[3618]=['ξ Cnc',2.4026564,0.3831556,5.14];
ybs[3619]=['κ Pyx',2.3956214,-0.4529156,4.58];
ybs[3620]=['',2.3876693,-0.9755436,6.11];
ybs[3621]=['19 Hya',2.3989765,-0.1515212,5.6];
ybs[3622]=['',2.3909768,-0.8954134,6.73];
ybs[3623]=['',2.3848022,-1.1273206,6.37];
ybs[3624]=['',2.4268782,1.2489857,6.55];
ybs[3625]=['λ Vel',2.3947071,-0.7596417,2.21];
ybs[3626]=['',2.4041709,0.2002232,6.48];
ybs[3627]=['',2.4010072,-0.217294,5.77];
ybs[3628]=['',2.3985469,-0.4687927,6.15];
ybs[3629]=['',2.4003173,-0.3215038,5.73];
ybs[3630]=['',2.4085528,0.5387869,5.95];
ybs[3631]=['79 Cnc',2.4069718,0.3822911,6.01];
ybs[3632]=['20 Hya',2.4028483,-0.1549885,5.46];
ybs[3633]=['',2.381553,-1.23272,4.71];
ybs[3634]=['',2.3788481,-1.2687379,4.48];
ybs[3635]=['ε Pyx',2.4037324,-0.5315887,5.59];
ybs[3636]=['',2.4351401,1.2714937,5.96];
ybs[3637]=['',2.405901,-0.4061267,6.53];
ybs[3638]=['',2.4020494,-0.8642368,6.48];
ybs[3639]=['16 UMa',2.4264492,1.0703943,5.13];
ybs[3640]=['',2.4134046,0.0938124,6.35];
ybs[3641]=['π1 Cnc',2.4152549,0.260101,6.51];
ybs[3642]=['',2.4146215,0.065866,6.14];
ybs[3643]=['36 Lyn',2.42281,0.7526517,5.32];
ybs[3644]=['',2.4129663,-0.3462916,5.73];
ybs[3645]=['',2.4080805,-0.7847167,5];
ybs[3646]=['21 Hya',2.415294,-0.1257191,6.11];
ybs[3647]=['',2.4110074,-0.6868224,6];
ybs[3648]=['',2.4212242,0.3698251,6.48];
ybs[3649]=['',2.4100882,-0.8146666,5.79];
ybs[3650]=['',2.4066327,-1.0307873,3.44];
ybs[3651]=['17 UMa',2.4324808,0.98867,5.27];
ybs[3652]=['',2.4144217,-0.7628317,5.57];
ybs[3653]=['18 UMa',2.4338236,0.9412046,4.83];
ybs[3654]=['',2.4076297,-1.0892624,3.97];
ybs[3655]=['',2.4287054,0.602821,5.97];
ybs[3656]=['θ Hya',2.4239624,0.0387464,3.88];
ybs[3657]=['',2.4530508,1.2901484,6.5];
ybs[3658]=['',2.4186731,-0.6756197,6.31];
ybs[3659]=['',2.4179812,-0.7394494,6.29];
ybs[3660]=['π2 Cnc',2.4280614,0.2591273,5.34];
ybs[3661]=['',2.4188835,-0.8278518,5.92];
ybs[3662]=['',2.4215045,-0.772131,5.85];
ybs[3663]=['',2.4151645,-1.0386102,5.54];
ybs[3664]=['',2.4227361,-0.7561048,5.25];
ybs[3665]=['',2.4281005,-0.2638806,6.35];
ybs[3666]=['',2.4391799,0.8154506,5.97];
ybs[3667]=['',2.4253626,-0.6579336,5.86];
ybs[3668]=['ζ Oct',2.3257134,-1.4923278,5.42];
ybs[3669]=['',2.4215113,-0.971516,5.27];
ybs[3670]=['',2.4262705,-0.7967422,6.25];
ybs[3671]=['23 Hya',2.433912,-0.1125395,5.24];
ybs[3672]=['',2.4281996,-0.6748238,4.94];
ybs[3673]=['24 Hya',2.4338224,-0.1542819,5.47];
ybs[3674]=['',2.4288579,-0.6546372,4.62];
ybs[3675]=['β Car',2.4148793,-1.2184282,1.68];
ybs[3676]=['',2.442657,0.6155517,5.75];
ybs[3677]=['',2.4355771,-0.2560178,5.84];
ybs[3678]=['',2.4299167,-0.7852816,6.04];
ybs[3679]=['',2.4394276,0.199067,6.41];
ybs[3680]=['38 Lyn',2.4445051,0.6406528,3.82];
ybs[3681]=['',2.4256093,-1.020721,6.02];
ybs[3682]=['',2.4313245,-0.7742395,5.12];
ybs[3683]=['',2.4269579,-1.006576,6.32];
ybs[3684]=['',2.4340272,-0.6893427,5.33];
ybs[3685]=['',2.4082896,-1.3396479,6.14];
ybs[3686]=['',2.4296734,-1.00594,4.34];
ybs[3687]=['',2.4534132,0.8930784,6.13];
ybs[3688]=['',2.4581119,0.9878969,5.47];
ybs[3689]=['ι Car',2.4333945,-1.0362075,2.25];
ybs[3690]=['',2.4364856,-0.9527848,6.33];
ybs[3691]=['',2.4538989,0.6648271,6.12];
ybs[3692]=['',2.4462979,-0.1991449,6.62];
ybs[3693]=['',2.4384325,-0.8926754,5.26];
ybs[3694]=['',2.4461441,-0.2780385,5.78];
ybs[3695]=['α Lyn',2.4540475,0.598577,3.13];
ybs[3696]=['26 Hya',2.4472048,-0.2106799,4.79];
ybs[3697]=['',2.4557314,0.5725595,6.16];
ybs[3698]=['',2.4410528,-0.9015706,5.87];
ybs[3699]=['27 Hya',2.4503651,-0.1684619,4.8];
ybs[3700]=['',2.4466608,-0.5968919,6.39];
ybs[3701]=['',2.4543572,0.2665903,6.53];
ybs[3702]=['',2.4329945,-1.200516,5.39];
ybs[3703]=['',2.4358023,-1.1719206,6.11];
ybs[3704]=['',2.4521376,-0.2742653,6.33];
ybs[3705]=['',2.4508409,-0.5560085,6.82];
ybs[3706]=['',2.4495672,-0.6575996,6.05];
ybs[3707]=['',2.4444744,-0.9648631,6.28];
ybs[3708]=['θ Pyx',2.4543298,-0.4548713,4.72];
ybs[3709]=['',2.4878004,1.3089845,6.29];
ybs[3710]=['',2.4319559,-1.3088133,5.29];
ybs[3711]=['',2.4321607,-1.3060259,5.86];
ybs[3712]=['',2.4762911,1.1142634,6.28];
ybs[3713]=['',2.4645401,0.4378271,6.41];
ybs[3714]=['',2.4606803,-0.1734164,6.53];
ybs[3715]=['',2.4717119,0.898425,6.31];
ybs[3716]=['',2.4552965,-0.7381304,5.58];
ybs[3717]=['',2.4685959,0.6368572,6.67];
ybs[3718]=['',2.4499463,-1.0908501,4.81];
ybs[3719]=['',2.4587351,-0.6958933,6.54];
ybs[3720]=['',2.457526,-0.8053724,5.75];
ybs[3721]=['κ Leo',2.4694742,0.4552591,4.46];
ybs[3722]=['',2.4544659,-0.9706075,5.63];
ybs[3723]=['λ Pyx',2.4617169,-0.5049435,4.69];
ybs[3724]=['κ Vel',2.4557283,-0.96181,2.5];
ybs[3725]=['',2.4637743,-0.660688,6.48];
ybs[3726]=['',2.4730765,0.2877606,6.29];
ybs[3727]=['',2.4660026,-0.6898139,6.06];
ybs[3728]=['28 Hya',2.4719316,-0.091028,5.59];
ybs[3729]=['',2.4641656,-0.9046862,6.08];
ybs[3730]=['',2.4611549,-1.0541748,6.3];
ybs[3731]=['',2.4762586,-0.0272663,6.01];
ybs[3732]=['',2.4638108,-1.0776774,5.99];
ybs[3733]=['',2.4876696,0.7941634,5.41];
ybs[3734]=['29 Hya',2.4798777,-0.162704,6.54];
ybs[3735]=['',2.4771863,-0.504155,6.1];
ybs[3736]=['',2.4755866,-0.7086087,6.2];
ybs[3737]=['',2.4931683,0.9712056,6.45];
ybs[3738]=['α Hya',2.4813897,-0.1528449,1.98];
ybs[3739]=['',2.4798189,-0.3916961,4.69];
ybs[3740]=['',2.4822953,-0.1076857,5.38];
ybs[3741]=['',2.5312894,1.4176301,4.29];
ybs[3742]=['',2.4697168,-1.0829457,5.77];
ybs[3743]=['',2.4741612,-0.9333571,5.11];
ybs[3744]=['ω Leo',2.4856112,0.1563397,5.41];
ybs[3745]=['3 Leo',2.4857141,0.1411843,5.71];
ybs[3746]=['',2.4808845,-0.6127244,6.65];
ybs[3747]=['23 UMa',2.5013721,1.0988903,3.67];
ybs[3748]=['',2.4878921,-0.02367,6.27];
ybs[3749]=['τ1 Hya',2.4883442,-0.050059,4.6];
ybs[3750]=['',2.4894919,-0.0402238,6.14];
ybs[3751]=['',2.4749868,-1.1349542,6.05];
ybs[3752]=['',2.4900189,-0.0758631,6.26];
ybs[3753]=['',2.4881726,-0.3638645,5.66];
ybs[3754]=['7 LMi',2.4961449,0.5856577,5.85];
ybs[3755]=['ε Ant',2.487863,-0.6292029,4.51];
ybs[3756]=['',2.4878894,-0.6720072,6.19];
ybs[3757]=['',2.4908168,-0.4091885,6.24];
ybs[3758]=['22 UMa',2.5174705,1.2584564,5.72];
ybs[3759]=['8 LMi',2.4997673,0.6109166,5.37];
ybs[3760]=['',2.4910599,-0.465815,5.48];
ybs[3761]=['24 UMa',2.515181,1.2170026,4.56];
ybs[3762]=['',2.4933991,-0.2736136,5.85];
ybs[3763]=['λ Leo',2.5001859,0.3991201,4.31];
ybs[3764]=['',2.5233288,1.2953144,6.46];
ybs[3765]=['θ UMa',2.5062364,0.9001824,3.17];
ybs[3766]=['',2.4842664,-1.0885989,5.92];
ybs[3767]=['',2.4754252,-1.2514127,5.47];
ybs[3768]=['',2.5072573,0.8611098,6.76];
ybs[3769]=['6 Leo',2.5009014,0.1678241,5.07];
ybs[3770]=['ζ1 Ant',2.4946036,-0.5583513,7];
ybs[3771]=['ζ1 Ant',2.4946547,-0.5583175,6.18];
ybs[3772]=['ξ Leo',2.5008734,0.1954682,4.97];
ybs[3773]=['',2.4825042,-1.1658956,5.91];
ybs[3774]=['',2.4908313,-0.9008823,5.45];
ybs[3775]=['',2.4990721,-0.1859182,6.14];
ybs[3776]=['ψ Vel',2.4940503,-0.7080177,3.6];
ybs[3777]=['τ2 Hya',2.5007443,-0.0224314,4.57];
ybs[3778]=['',2.5003058,-0.1827491,6.13];
ybs[3779]=['ζ2 Ant',2.4980114,-0.5580165,5.93];
ybs[3780]=['',2.4979342,-0.6250904,5.87];
ybs[3781]=['9 LMi',2.5083783,0.6350588,6.18];
ybs[3782]=['',2.5072476,0.4933588,6.53];
ybs[3783]=['',2.4916722,-1.0203418,5.88];
ybs[3784]=['',2.5039046,0.0307825,6.11];
ybs[3785]=['ι Cha',2.4580679,-1.4116958,5.36];
ybs[3786]=['',2.5018738,-0.3403497,5.74];
ybs[3787]=['',2.5123832,0.8168348,6.52];
ybs[3788]=['',2.5014816,-0.5014044,6.46];
ybs[3789]=['26 UMa',2.5148278,0.9067017,4.5];
ybs[3790]=['10 LMi',2.5115055,0.6334936,4.55];
ybs[3791]=['',2.5051725,-0.1502002,6.12];
ybs[3792]=['',2.5045956,-0.2376696,5.94];
ybs[3793]=['',2.4954187,-0.9971823,3.13];
ybs[3794]=['',2.5100745,0.4075866,6.25];
ybs[3795]=['',2.5065123,-0.127246,6.24];
ybs[3796]=['',2.5308912,1.273711,6.42];
ybs[3797]=['',2.5011295,-0.711217,5.35];
ybs[3798]=['',2.5056206,-0.3702967,5.01];
ybs[3799]=['',2.5152823,0.6897563,4.81];
ybs[3800]=['',2.5065714,-0.4008073,5.91];
ybs[3801]=['',2.516646,0.6957226,6.76];
ybs[3802]=['',2.5047207,-0.684683,6.43];
ybs[3803]=['',2.4958125,-1.1662187,6.27];
ybs[3804]=['33 Hya',2.5118222,-0.1049999,5.56];
ybs[3805]=['11 LMi',2.5177459,0.6232367,5.41];
ybs[3806]=['',2.4993578,-1.0976218,6.1];
ybs[3807]=['',2.5069544,-0.8570569,5.12];
ybs[3808]=['7 Leo',2.5181168,0.249202,6.36];
ybs[3809]=['',2.5086007,-0.8963339,5.01];
ybs[3810]=['',2.5221947,0.5420973,5.56];
ybs[3811]=['',2.4947821,-1.2772454,5.47];
ybs[3812]=['',2.515935,-0.3435677,6.31];
ybs[3813]=['',2.5138694,-0.6270116,6.49];
ybs[3814]=['',2.5364038,1.1723288,5.94];
ybs[3815]=['',2.5093425,-1.0355106,4.08];
ybs[3816]=['8 Leo',2.5232323,0.2851153,5.69];
ybs[3817]=['10 Leo',2.523744,0.117529,5];
ybs[3818]=['',2.5201679,-0.4329195,6.53];
ybs[3819]=['42 Lyn',2.529676,0.7005299,5.25];
ybs[3820]=['',2.5220805,-0.4432873,5.7];
ybs[3821]=['',2.5186851,-0.8526457,6.17];
ybs[3822]=['34 Hya',2.526203,-0.1662697,6.4];
ybs[3823]=['',2.522587,-0.5634006,5.63];
ybs[3824]=['',2.5291193,0.0793577,4.68];
ybs[3825]=['',2.5237995,-0.6317706,5.98];
ybs[3826]=['',2.5204243,-0.8631877,4.35];
ybs[3827]=['',2.5199798,-0.9258252,6.19];
ybs[3828]=['',2.5487947,1.2066141,5.69];
ybs[3829]=['27 UMa',2.5524543,1.2592317,5.17];
ybs[3830]=['',2.5218297,-0.9384716,5.45];
ybs[3831]=['',2.5159265,-1.1353717,6.56];
ybs[3832]=['',2.5259511,-0.7556095,5.5];
ybs[3833]=['',2.565505,1.3618812,6.23];
ybs[3834]=['',2.5289469,-0.6931837,6.7];
ybs[3835]=['ι Hya',2.5351077,-0.0217384,3.91];
ybs[3836]=['37 Hya',2.5346082,-0.1862789,6.31];
ybs[3837]=['',2.5737614,1.3793588,6.17];
ybs[3838]=['',2.5369825,-0.1897531,6.37];
ybs[3839]=['κ Hya',2.5367773,-0.25194,5.06];
ybs[3840]=['',2.5434319,0.544102,5.89];
ybs[3841]=['43 Lyn',2.5455314,0.6920986,5.62];
ybs[3842]=['ο Leo',2.5409921,0.1708514,3.52];
ybs[3843]=['13 Leo',2.5435214,0.45046,6.24];
ybs[3844]=['',2.549001,0.8434727,6.39];
ybs[3845]=['',2.5510614,0.9470121,6.47];
ybs[3846]=['',2.5306202,-1.0721655,4.52];
ybs[3847]=['13 LMi',2.5484518,0.610685,6.14];
ybs[3848]=['',2.5408159,-0.413553,4.77];
ybs[3849]=['',2.5584698,1.1323624,6.17];
ybs[3850]=['ζ Cha',2.5008586,-1.4144484,5.11];
ybs[3851]=['15 Leo',2.5519734,0.5213391,5.64];
ybs[3852]=['',2.5449852,-0.419211,4.94];
ybs[3853]=['',2.5368535,-1.0138015,5.32];
ybs[3854]=['',2.5383453,-1.0011692,5.8];
ybs[3855]=['28 UMa',2.564023,1.109133,6.34];
ybs[3856]=['ψ Leo',2.5523451,0.2429099,5.35];
ybs[3857]=['',2.5466181,-0.621429,6.41];
ybs[3858]=['',2.5418145,-0.9654717,6];
ybs[3859]=['',2.5558076,0.3274136,6.5];
ybs[3860]=['',2.5660916,0.9952427,5.2];
ybs[3861]=['θ Ant',2.5534501,-0.4864844,4.79];
ybs[3862]=['',2.5493227,-0.8959148,6.15];
ybs[3863]=['ε Leo',2.5618204,0.4131118,2.98];
ybs[3864]=['',2.553354,-0.6924625,6.82];
ybs[3865]=['',2.5502398,-0.9424,5.56];
ybs[3866]=['',2.5628064,0.1152603,5.79];
ybs[3867]=['18 Leo',2.563883,0.2042951,5.63];
ybs[3868]=['',2.5584566,-0.5289602,6.45];
ybs[3869]=['',2.5636916,0.0293356,5.65];
ybs[3870]=['19 Leo',2.5684285,0.2000718,6.45];
ybs[3871]=['',2.574461,0.8013796,5.09];
ybs[3872]=['',2.568978,0.1976373,6.02];
ybs[3873]=['',2.5586464,-0.9998995,6.46];
ybs[3874]=['',2.5563213,-1.0927872,3.69];
ybs[3875]=['',2.5836759,1.1429693,6.31];
ybs[3876]=['',2.5629587,-0.7829502,5.55];
ybs[3877]=['',2.559612,-1.0279763,6.22];
ybs[3878]=['υ UMa',2.5856419,1.0285651,3.8];
ybs[3879]=['20 Leo',2.5791249,0.3678048,6.09];
ybs[3880]=['υ Car',2.5641673,-1.13755,3.01];
ybs[3881]=['υ Car',2.564211,-1.1375598,6.26];
ybs[3882]=['',2.5761793,-0.6508685,5.97];
ybs[3883]=['4 Sex',2.5816688,0.0739609,6.24];
ybs[3884]=['φ UMa',2.59017,0.9417442,4.59];
ybs[3885]=['',2.5717906,-0.9864131,6.06];
ybs[3886]=['23 Leo',2.5841696,0.2261945,6.46];
ybs[3887]=['',2.5778981,-0.6348523,6.37];
ybs[3888]=['',2.5779774,-0.8000334,5.08];
ybs[3889]=['6 Sex',2.5846867,-0.0759131,6.01];
ybs[3890]=['22 Leo',2.5881322,0.4239212,5.32];
ybs[3891]=['',2.5851999,-0.1097441,6.42];
ybs[3892]=['ν Cha',2.5582922,-1.3418207,5.45];
ybs[3893]=['υ1 Hya',2.5855346,-0.2609774,4.12];
ybs[3894]=['',2.5811974,-0.8210103,5.73];
ybs[3895]=['μ Leo',2.5920081,0.4520457,3.88];
ybs[3896]=['7 Sex',2.5890545,0.0409753,6.02];
ybs[3897]=['',2.5889922,-0.0005393,6.35];
ybs[3898]=['',2.5877461,-0.290442,6.08];
ybs[3899]=['γ Sex',2.5901649,-0.1433183,5.05];
ybs[3900]=['',2.5839836,-0.8080884,5.62];
ybs[3901]=['',2.6034948,1.0648041,6.27];
ybs[3902]=['',2.5854986,-0.8142666,4.58];
ybs[3903]=['',2.5826643,-1.0390283,5.79];
ybs[3904]=['',2.5811771,-1.0969619,5.57];
ybs[3905]=['',2.59572,0.1021269,5.95];
ybs[3906]=['',2.591743,-0.4788988,6.3];
ybs[3907]=['31 UMa',2.6056626,0.8676468,5.27];
ybs[3908]=['',2.6195109,1.270096,5.83];
ybs[3909]=['',2.5971862,-0.4544752,4.88];
ybs[3910]=['',2.5908176,-0.968308,6.48];
ybs[3911]=['',2.5986803,-0.3943648,6.24];
ybs[3912]=['',2.6126573,1.0002552,5.93];
ybs[3913]=['',2.6002509,-0.3336484,4.94];
ybs[3914]=['',2.5947148,-0.8945479,5.93];
ybs[3915]=['',2.5969758,-0.7922208,5.71];
ybs[3916]=['',2.6076377,0.1540322,5.85];
ybs[3917]=['',2.599213,-0.8787917,5.72];
ybs[3918]=['19 LMi',2.6138865,0.7146691,5.14];
ybs[3919]=['',2.6151917,0.7907447,6.3];
ybs[3920]=['',2.6049943,-0.7144025,6.41];
ybs[3921]=['',2.6083953,-0.4652701,6.28];
ybs[3922]=['',2.6074153,-0.5851441,5.84];
ybs[3923]=['',2.6089185,-0.4814101,6.32];
ybs[3924]=['',2.6697975,1.462711,6.37];
ybs[3925]=['',2.6057629,-0.897862,6.37];
ybs[3926]=['',2.616775,0.4825951,6.3];
ybs[3927]=['ν Leo',2.6155157,0.2153137,5.26];
ybs[3928]=['',2.6150119,0.1432224,6.04];
ybs[3929]=['',2.6240819,0.9896592,5.48];
ybs[3930]=['φ Vel',2.6077682,-0.9542676,3.54];
ybs[3931]=['',2.6092677,-0.9206037,6.12];
ybs[3932]=['',2.6219144,0.5155132,5.73];
ybs[3933]=['',2.6117749,-0.8468759,6.05];
ybs[3934]=['',2.6029322,-1.2478566,6.35];
ybs[3935]=['12 Sex',2.6218601,0.0571797,6.7];
ybs[3936]=['',2.6186177,-0.4199028,6.21];
ybs[3937]=['η Ant',2.6173059,-0.6283084,5.23];
ybs[3938]=['',2.6087026,-1.1274347,6.58];
ybs[3939]=['',2.6069831,-1.2079365,6.2];
ybs[3940]=['π Leo',2.6241089,0.1385,4.7];
ybs[3941]=['20 LMi',2.6281109,0.5552709,5.36];
ybs[3942]=['',2.6357362,0.3811756,5.66];
ybs[3943]=['',2.6238561,-0.9958047,6.52];
ybs[3944]=['',2.6445535,0.9386687,5.74];
ybs[3945]=['',2.6288845,-0.9332885,6.2];
ybs[3946]=['',2.6347065,-0.5355871,6.54];
ybs[3947]=['',2.6299277,-1.002846,6.2];
ybs[3948]=['',2.6469574,0.9121226,6.14];
ybs[3949]=['',2.6389178,-0.1690091,6.12];
ybs[3950]=['',2.6298566,-1.056447,5.94];
ybs[3951]=['13 Sex',2.6411477,0.0539544,6.45];
ybs[3952]=['',2.6386244,-0.4437723,6.7];
ybs[3953]=['',2.6403438,-0.3178485,5.86];
ybs[3954]=['',2.6364801,-0.8158649,6.12];
ybs[3955]=['',2.6415299,-0.4257791,5.7];
ybs[3956]=['',2.6342458,-1.0522241,6.19];
ybs[3957]=['',2.6333203,-1.0867419,6.42];
ybs[3958]=['',2.6412964,-0.699626,6.43];
ybs[3959]=['',2.6481091,0.2730977,6.37];
ybs[3960]=['υ2 Hya',2.6451439,-0.2299422,4.6];
ybs[3961]=['',2.6366204,-1.0819894,6.14];
ybs[3962]=['',2.6451758,-0.6369388,6.27];
ybs[3963]=['14 Sex',2.6527529,0.0960096,6.21];
ybs[3964]=['21 LMi',2.6561684,0.6132056,4.48];
ybs[3965]=['η Leo',2.6553315,0.2906356,3.52];
ybs[3966]=['',2.6488977,-0.8286867,5.08];
ybs[3967]=['',2.6539413,-0.3011075,5.6];
ybs[3968]=['',2.6483807,-0.9127772,6.52];
ybs[3969]=['',2.6597139,0.5496624,6.24];
ybs[3970]=['31 Leo',2.6576968,0.1725567,4.37];
ybs[3971]=['α Sex',2.6576555,-0.0084194,4.49];
ybs[3972]=['α Leo',2.6597687,0.2069327,1.35];
ybs[3973]=['μ1 Cha',2.6181748,-1.4368128,5.52];
ybs[3974]=['',2.6572747,-0.6535271,6.36];
ybs[3975]=['',2.6610594,-0.1919104,6.53];
ybs[3976]=['',2.6602286,-0.2744104,6.27];
ybs[3977]=['',2.6717955,0.7077286,6.32];
ybs[3978]=['',2.6661879,-0.2130535,6.24];
ybs[3979]=['17 Sex',2.6670549,-0.1486953,5.91];
ybs[3980]=['',2.6607337,-0.9062109,4.86];
ybs[3981]=['',2.6668585,-0.2256254,5.31];
ybs[3982]=['',2.6638771,-0.6277562,6.13];
ybs[3983]=['',2.6727112,0.6508397,5.85];
ybs[3984]=['λ Hya',2.6690059,-0.2175651,3.61];
ybs[3985]=['',2.6587581,-1.1506282,5.28];
ybs[3986]=['18 Sex',2.6705682,-0.1488734,5.65];
ybs[3987]=['μ2 Cha',2.633829,-1.4255039,6.6];
ybs[3988]=['34 Leo',2.6740343,0.2311397,6.44];
ybs[3989]=['',2.6620165,-1.0761736,5.6];
ybs[3990]=['',2.6721881,-0.1296473,6.25];
ybs[3991]=['',2.6685217,-0.7300082,5.98];
ybs[3992]=['',2.6619145,-1.2006886,5.81];
ybs[3993]=['',2.6750564,-0.5012261,6.28];
ybs[3994]=['19 Sex',2.6789761,0.0785881,5.77];
ybs[3995]=['',2.677792,-0.3362466,6.44];
ybs[3996]=['',2.6838671,0.4716511,6.04];
ybs[3997]=['',2.6719373,-1.0286912,6.4];
ybs[3998]=['',2.6906923,1.0449808,6.25];
ybs[3999]=['',2.6728104,-1.0152966,5.72];
ybs[4000]=['',2.6757866,-0.9123734,6.16];
ybs[4001]=['',2.6806601,-0.473699,6.25];
ybs[4002]=['',2.6866449,0.3674861,6.02];
ybs[4003]=['',2.6809162,-0.5784724,6.38];
ybs[4004]=['22 LMi',2.6895202,0.5472573,6.46];
ybs[4005]=['',2.6822582,-0.7061252,5.9];
ybs[4006]=['',2.7047005,1.2733932,6.4];
ybs[4007]=['',2.6801877,-0.8961462,5.28];
ybs[4008]=['',2.6781467,-1.0477214,6.1];
ybs[4009]=['',2.6830387,-0.7055103,6.35];
ybs[4010]=['',2.6805438,-0.9052706,5.78];
ybs[4011]=['',2.7035976,1.2382644,6.66];
ybs[4012]=['',2.6794733,-1.078106,6.41];
ybs[4013]=['',2.686452,-0.7371283,3.85];
ybs[4014]=['23 LMi',2.6944185,0.509597,5.35];
ybs[4015]=['',2.6797349,-1.1603842,5.16];
ybs[4016]=['32 UMa',2.7037185,1.1343781,5.82];
ybs[4017]=['24 LMi',2.6954002,0.4986345,6.49];
ybs[4018]=['',2.6943128,0.3076576,6.55];
ybs[4019]=['',2.6892832,-0.6393246,6.19];
ybs[4020]=['35 Leo',2.6955991,0.4082359,5.97];
ybs[4021]=['ζ Leo',2.6962585,0.4067371,3.44];
ybs[4022]=['',2.6963338,0.4408438,5.84];
ybs[4023]=['λ UMa',2.6984991,0.747026,3.45];
ybs[4024]=['',2.6933109,-0.197503,6.08];
ybs[4025]=['37 Leo',2.696029,0.2376343,5.41];
ybs[4026]=['',2.6898672,-0.7544201,5.6];
ybs[4027]=['ω Car',2.6802693,-1.2243513,3.32];
ybs[4028]=['',2.6883086,-0.9614439,6.16];
ybs[4029]=['39 Leo',2.6986648,0.401305,5.82];
ybs[4030]=['',2.6957965,-0.3627397,6.57];
ybs[4031]=['',2.7027979,0.4765103,6.52];
ybs[4032]=['ε Sex',2.6998169,-0.1428028,5.24];
ybs[4033]=['',2.6914135,-1.0474772,6.22];
ybs[4034]=['',2.7068397,0.8141504,6.43];
ybs[4035]=['',2.6945795,-0.8956654,6.3];
ybs[4036]=['',2.7089148,0.842704,6];
ybs[4037]=['',2.7172291,1.1978812,5.96];
ybs[4038]=['',2.7064052,0.4293199,6.4];
ybs[4039]=['',2.7015984,-0.5079808,5.34];
ybs[4040]=['',2.6958317,-1.0724203,3.4];
ybs[4041]=['',2.7126158,0.9366382,6.45];
ybs[4042]=['',2.7138254,0.9442778,6];
ybs[4043]=['',2.7036159,-0.6443415,6.3];
ybs[4044]=['40 Leo',2.7094604,0.3378473,4.79];
ybs[4045]=['',2.7069438,-0.2206365,6];
ybs[4046]=['',2.7027842,-0.7292268,5.96];
ybs[4047]=['γ1 Leo',2.710499,0.3443186,2.61];
ybs[4048]=['γ2 Leo',2.7105207,0.3442992,3.8];
ybs[4049]=['',2.708187,-0.0910954,6.37];
ybs[4050]=['',2.7101009,-0.160091,6.32];
ybs[4051]=['',2.702909,-0.9812818,5.81];
ybs[4052]=['',2.7607043,1.4684545,5.5];
ybs[4053]=['',2.7072816,-0.9624266,4.57];
ybs[4054]=['23 Sex',2.7148363,0.0379753,6.66];
ybs[4055]=['',2.7043099,-1.130795,5.67];
ybs[4056]=['',2.7104962,-0.8344919,5.65];
ybs[4057]=['',2.7205689,0.7175969,5.76];
ybs[4058]=['',2.7149264,-0.3158856,6.51];
ybs[4059]=['μ UMa',2.7212442,0.7223086,3.05];
ybs[4060]=['42 Leo',2.7185539,0.2593816,6.12];
ybs[4061]=['',2.7163433,-0.4158215,6.5];
ybs[4062]=['',2.7302541,1.1423486,4.97];
ybs[4063]=['',2.7168968,-0.3951835,6.51];
ybs[4064]=['',2.7129194,-0.9801226,4.5];
ybs[4065]=['27 LMi',2.7244381,0.5898109,5.9];
ybs[4066]=['',2.7196306,-0.348736,6.13];
ybs[4067]=['43 Leo',2.7235126,0.1121925,6.07];
ybs[4068]=['',2.7269245,0.5148953,6.39];
ybs[4069]=['',2.7245317,0.0973854,6.54];
ybs[4070]=['',2.7196423,-0.7289222,4.83];
ybs[4071]=['28 LMi',2.7289583,0.5865005,5.5];
ybs[4072]=['25 Sex',2.7252449,-0.0731049,5.97];
ybs[4073]=['',2.7238154,-0.5284263,6.27];
ybs[4074]=['',2.76482,1.4388912,5.26];
ybs[4075]=['',2.7287336,0.0393301,6.32];
ybs[4076]=['',2.7248163,-0.6653968,5.33];
ybs[4077]=['',2.7255203,-0.7342217,6.27];
ybs[4078]=['44 Leo',2.7333481,0.151318,5.61];
ybs[4079]=['',2.7211389,-1.1696489,4.99];
ybs[4080]=['30 LMi',2.7366801,0.5878464,4.74];
ybs[4081]=['',2.7256636,-1.0134845,6.35];
ybs[4082]=['',2.7352217,-0.1252215,5.57];
ybs[4083]=['',2.7325582,-0.7432115,6.18];
ybs[4084]=['μ Hya',2.7366017,-0.2958577,3.81];
ybs[4085]=['',2.7306264,-1.0243535,5.95];
ybs[4086]=['',2.7436274,0.7240586,6.02];
ybs[4087]=['',2.7411722,0.3359623,6.15];
ybs[4088]=['',2.7464408,0.8494388,6.44];
ybs[4089]=['',2.7363609,-0.7479417,6.13];
ybs[4090]=['β LMi',2.7453217,0.6386475,4.21];
ybs[4091]=['45 Leo',2.7438173,0.1683744,6.04];
ybs[4092]=['',2.7263683,-1.2940959,4];
ybs[4093]=['',2.7487006,0.7870849,6.35];
ybs[4094]=['α Ant',2.7409796,-0.5442462,4.25];
ybs[4095]=['',2.7278892,-1.2930501,6.19];
ybs[4096]=['35 UMa',2.7553391,1.1433693,6.32];
ybs[4097]=['',2.738809,-0.9598028,5.58];
ybs[4098]=['',2.7575682,1.1194807,6.12];
ybs[4099]=['',2.7483425,-0.0673362,6.05];
ybs[4100]=['',2.7412563,-1.0080002,4.66];
ybs[4101]=['',2.744343,-0.8643039,6.1];
ybs[4102]=['36 UMa',2.7578789,0.9750203,4.84];
ybs[4103]=['32 LMi',2.7550596,0.6773517,5.77];
ybs[4104]=['',2.7432617,-1.0272103,3.82];
ybs[4105]=['',2.7407735,-1.1487753,6.01];
ybs[4106]=['δ Sex',2.7516083,-0.0498275,5.21];
ybs[4107]=['',2.7511933,-0.5197476,5.58];
ybs[4108]=['δ Ant',2.7516408,-0.536217,5.56];
ybs[4109]=['β Sex',2.7551879,-0.0131398,5.09];
ybs[4110]=['',2.7472772,-1.1220337,5.29];
ybs[4111]=['',2.7848189,1.4028482,6.52];
ybs[4112]=['',2.7580829,-0.1353249,6.2];
ybs[4113]=['',2.7580743,-0.2391867,5.58];
ybs[4114]=['33 LMi',2.7625419,0.5630991,5.9];
ybs[4115]=['',2.7572514,-0.464256,6.51];
ybs[4116]=['',2.7792048,1.319401,4.84];
ybs[4117]=['46 Leo',2.7637137,0.2447112,5.46];
ybs[4118]=['',2.7552577,-1.0728899,6.43];
ybs[4119]=['',2.7525932,-1.1711304,6.19];
ybs[4120]=['',2.7613836,-0.4948656,6.05];
ybs[4121]=['',2.7712613,0.9316719,6.45];
ybs[4122]=['',2.7687081,0.7035253,4.75];
ybs[4123]=['ρ Leo',2.7663249,0.1604,3.85];
ybs[4124]=['',2.7587537,-0.9395397,4.89];
ybs[4125]=['',2.7616515,-0.7885903,5.74];
ybs[4126]=['',2.7615859,-0.7886388,6.09];
ybs[4127]=['34 LMi',2.7698211,0.6086318,5.58];
ybs[4128]=['',2.7528211,-1.258538,4.74];
ybs[4129]=['',2.7642948,-0.7807773,5.91];
ybs[4130]=['',2.761237,-1.0786398,3.32];
ybs[4131]=['37 UMa',2.7776566,0.9942418,5.16];
ybs[4132]=['',2.755687,-1.2799837,4.93];
ybs[4133]=['',2.7659233,-0.822395,5.02];
ybs[4134]=['',2.764783,-1.0259627,6];
ybs[4135]=['44 Hya',2.7710927,-0.4164694,5.08];
ybs[4136]=['48 Leo',2.7749646,0.1193245,5.08];
ybs[4137]=['',2.7675582,-1.0176455,6.14];
ybs[4138]=['49 Leo',2.7760208,0.148936,5.67];
ybs[4139]=['',2.7752364,-0.4065387,6.1];
ybs[4140]=['35 LMi',2.7822299,0.6319804,6.28];
ybs[4141]=['',2.7708785,-1.0664738,6.23];
ybs[4142]=['',2.7783073,-0.3261348,6.49];
ybs[4143]=['',2.7760293,-0.6925408,5.38];
ybs[4144]=['',2.7757534,-0.7641331,6.08];
ybs[4145]=['',2.7812248,-0.1867579,6.57];
ybs[4146]=['φ2 Hya',2.7810913,-0.2873083,6.03];
ybs[4147]=['',2.7800478,-0.4676097,6.29];
ybs[4148]=['',2.7822929,-0.2155034,5.7];
ybs[4149]=['',2.777074,-1.0066138,4.45];
ybs[4150]=['',2.7851516,-0.2070989,6.52];
ybs[4151]=['',2.7561702,-1.4318193,7.07];
ybs[4152]=['',2.7850656,-0.4804854,4.89];
ybs[4153]=['',2.7867002,-0.2356508,4.82];
ybs[4154]=['',2.7802565,-1.0416441,5.08];
ybs[4155]=['',2.7946269,0.9346354,5.52];
ybs[4156]=['37 LMi',2.7924428,0.5560361,4.71];
ybs[4157]=['',2.7849238,-0.8437467,3.84];
ybs[4158]=['38 LMi',2.7943306,0.6596006,5.85];
ybs[4159]=['',2.7851519,-1.0271374,5.45];
ybs[4160]=['',2.7743009,-1.3338859,6.3];
ybs[4161]=['φ3 Hya',2.7911471,-0.296605,4.91];
ybs[4162]=['',2.79233,-0.2192344,6.04];
ybs[4163]=['',2.7878271,-1.0013619,5.91];
ybs[4164]=['γ Cha',2.7738585,-1.3740042,4.11];
ybs[4165]=['',2.7917916,-0.7482435,6.11];
ybs[4166]=['',2.8072627,1.1924988,5.75];
ybs[4167]=['',2.7908429,-1.0349909,4.66];
ybs[4168]=['38 UMa',2.807637,1.1449043,5.12];
ybs[4169]=['',2.7919046,-1.0286018,5.92];
ybs[4170]=['',2.7934356,-0.9725148,4.28];
ybs[4171]=['',2.8128058,1.2035389,5];
ybs[4172]=['33 Sex',2.8036588,-0.0324586,6.26];
ybs[4173]=['',2.8007806,-0.6258687,6.37];
ybs[4174]=['',2.807552,0.5511526,6.02];
ybs[4175]=['',2.7967763,-1.1382753,5.52];
ybs[4176]=['',2.7916817,-1.3022117,6.07];
ybs[4177]=['39 UMa',2.8148933,0.9962453,5.8];
ybs[4178]=['',2.8019468,-1.0436192,6.42];
ybs[4179]=['40 LMi',2.8111323,0.4574016,5.51];
ybs[4180]=['',2.8083813,-0.2459741,6.24];
ybs[4181]=['',2.8137799,0.8043422,5.18];
ybs[4182]=['41 LMi',2.8127704,0.4026456,5.08];
ybs[4183]=['35 Sex',2.8122269,0.0807975,5.79];
ybs[4184]=['',2.8089641,-0.5730638,5.64];
ybs[4185]=['',2.8213322,1.1744779,6];
ybs[4186]=['',2.8057836,-1.1272136,4.82];
ybs[4187]=['',2.8163184,0.342783,6.27];
ybs[4188]=['',2.8080267,-1.0355757,5.38];
ybs[4189]=['θ Car',2.8089582,-1.1259603,2.76];
ybs[4190]=['',2.8117179,-1.0591548,4.57];
ybs[4191]=['36 Sex',2.8200883,0.0413523,6.28];
ybs[4192]=['41 UMa',2.8264404,0.9991462,6.34];
ybs[4193]=['42 LMi',2.8235574,0.533431,5.24];
ybs[4194]=['',2.8128972,-1.1234227,5.77];
ybs[4195]=['',2.8140627,-1.1184008,4.82];
ybs[4196]=['',2.8015775,-1.3945426,5.97];
ybs[4197]=['',2.8242291,0.1091555,6.37];
ybs[4198]=['51 Leo',2.8257536,0.3276406,5.49];
ybs[4199]=['52 Leo',2.8257507,0.2456683,5.48];
ybs[4200]=['η Car',2.81844,-1.043757,6.21];
ybs[4201]=['',2.8143823,-1.2388095,6.26];
ybs[4202]=['',2.815312,-1.2387229,6.46];
ybs[4203]=['',2.8147028,-1.2664539,6.27];
ybs[4204]=['',2.8273072,-0.3039613,5.42];
ybs[4205]=['',2.8375241,1.1346879,6.39];
ybs[4206]=['μ Vel',2.8262999,-0.8646187,2.69];
ybs[4207]=['',2.8237375,-1.0598031,6.25];
ybs[4208]=['',2.8306806,-0.268451,6.67];
ybs[4209]=['',2.8234922,-1.1280745,5.34];
ybs[4210]=['',2.8244725,-1.1236828,5.23];
ybs[4211]=['',2.8268789,-0.9926779,5.23];
ybs[4212]=['',2.8260351,-1.1257783,4.85];
ybs[4213]=['43 LMi',2.8369968,0.5113193,6.15];
ybs[4214]=['',2.8353942,-0.036272,5.93];
ybs[4215]=['',2.833084,-0.5551473,5.88];
ybs[4216]=['',2.829856,-1.0050814,6.36];
ybs[4217]=['53 Leo',2.8380703,0.1819652,5.34];
ybs[4218]=['',2.8316971,-1.0478675,6];
ybs[4219]=['40 Sex',2.8380395,-0.0723197,6.61];
ybs[4220]=['44 LMi',2.8410795,0.48615,6.04];
ybs[4221]=['δ1 Cha',2.8163237,-1.4065329,5.47];
ybs[4222]=['ν Hya',2.8393643,-0.2847174,3.11];
ybs[4223]=['',2.8398774,-0.1740493,5.86];
ybs[4224]=['δ2 Cha',2.8185815,-1.4077659,4.45];
ybs[4225]=['43 UMa',2.8473306,0.9854557,5.67];
ybs[4226]=['42 UMa',2.848352,1.0332383,5.58];
ybs[4227]=['41 Sex',2.8424056,-0.157383,5.79];
ybs[4228]=['',2.8405312,-0.5965116,5.61];
ybs[4229]=['',2.8374678,-1.0374769,5.91];
ybs[4230]=['',2.8459125,-0.0560641,5.95];
ybs[4231]=['',2.8529773,0.9153432,6.65];
ybs[4232]=['',2.8530553,0.9142668,6.44];
ybs[4233]=['',2.8581944,1.2170834,5.93];
ybs[4234]=['',2.8509254,0.0158015,6.38];
ybs[4235]=['',2.8525408,-0.0056089,6.31];
ybs[4236]=['44 UMa',2.8576587,0.9505911,5.1];
ybs[4237]=['46 LMi',2.8560633,0.5950683,3.83];
ybs[4238]=['ω UMa',2.8591269,0.7517098,4.71];
ybs[4239]=['',2.8560674,-0.0414582,6.12];
ybs[4240]=['',2.8511838,-1.0011296,5.25];
ybs[4241]=['',2.8561987,-0.3535862,5.24];
ybs[4242]=['',2.8565023,-0.2716723,6.38];
ybs[4243]=['',2.8574361,-0.039258,5.45];
ybs[4244]=['48 LMi',2.8620035,0.4427993,6.2];
ybs[4245]=['',2.8597887,-0.2422219,5.66];
ybs[4246]=['',2.8632861,0.5919175,5.72];
ybs[4247]=['',2.8554037,-1.0292805,3.78];
ybs[4248]=['46 UMa',2.8666331,0.582704,5.03];
ybs[4249]=['54 Leo',2.8659606,0.429862,4.5];
ybs[4250]=['54 Leo',2.865997,0.4298474,6.3];
ybs[4251]=['',2.863623,-0.3627732,6.44];
ybs[4252]=['',2.8555568,-1.2363981,5.99];
ybs[4253]=['',2.8625287,-0.7395213,6.11];
ybs[4254]=['',2.8689605,0.7310798,6.03];
ybs[4255]=['55 Leo',2.8660989,0.0107598,5.91];
ybs[4256]=['',2.8596515,-1.0811777,5.93];
ybs[4257]=['56 Leo',2.8675449,0.1058503,5.81];
ybs[4258]=['',2.8484607,-1.3906668,6.33];
ybs[4259]=['',2.8688444,0.3880063,6.14];
ybs[4260]=['50 LMi',2.8701548,0.4429544,6.35];
ybs[4261]=['',2.8631865,-1.0583208,5.92];
ybs[4262]=['',2.887155,1.3552288,6.2];
ybs[4263]=['ι Ant',2.8700573,-0.6502813,4.6];
ybs[4264]=['',2.8715883,-0.8881222,5.91];
ybs[4265]=['',2.8824836,0.9034039,6.17];
ybs[4266]=['',2.8742554,-1.0446265,6.11];
ybs[4267]=['47 UMa',2.8829782,0.7035295,5.05];
ybs[4268]=['',2.8832574,0.6278305,6];
ybs[4269]=['',2.8706572,-1.3128431,6.13];
ybs[4270]=['',2.8864608,0.7924666,5.47];
ybs[4271]=['',2.8835666,0.2021929,6.55];
ybs[4272]=['',2.8810923,-0.5909368,5.71];
ybs[4273]=['',2.8873841,0.8967641,6.43];
ybs[4274]=['',2.8825423,-0.2875411,5.89];
ybs[4275]=['',2.8868396,0.7468308,6.02];
ybs[4276]=['',2.8907051,1.1047911,6.39];
ybs[4277]=['α Crt',2.8836565,-0.3214885,4.08];
ybs[4278]=['49 UMa',2.8889399,0.6822671,5.08];
ybs[4279]=['',2.8855277,-0.2479141,5.88];
ybs[4280]=['',2.8804376,-1.0723517,6.16];
ybs[4281]=['58 Leo',2.887303,0.0610227,4.84];
ybs[4282]=['',2.8842419,-0.7666933,5.81];
ybs[4283]=['',2.8849923,-0.7390932,4.39];
ybs[4284]=['59 Leo',2.8881415,0.1043743,4.99];
ybs[4285]=['β UMa',2.8936752,0.9819425,2.37];
ybs[4286]=['',2.8847368,-0.9065041,6.15];
ybs[4287]=['',2.888828,-0.2777514,6.34];
ybs[4288]=['',2.8874468,-0.5578178,6.07];
ybs[4289]=['61 Leo',2.8927785,-0.0454842,4.74];
ybs[4290]=['60 Leo',2.8951889,0.3500838,4.42];
ybs[4291]=['α UMa',2.9020676,1.0756331,1.79];
ybs[4292]=['',2.8950529,-0.470415,6.23];
ybs[4293]=['',2.8989692,-0.0152545,6.14];
ybs[4294]=['',2.8775194,-1.4255328,6.71];
ybs[4295]=['',2.8988919,-0.1994061,5.5];
ybs[4296]=['62 Leo',2.9005761,-0.0021362,5.95];
ybs[4297]=['',2.8987539,-0.5599427,6.46];
ybs[4298]=['',2.9004428,-0.236597,6.34];
ybs[4299]=['51 UMa',2.9049548,0.6653142,6];
ybs[4300]=['χ Leo',2.90678,0.1259143,4.63];
ybs[4301]=['',2.9039757,-0.8342821,5.67];
ybs[4302]=['η Oct',2.8753298,-1.4785514,6.19];
ybs[4303]=['',2.9058471,-0.6270349,5.43];
ybs[4304]=['χ1 Hya',2.9078242,-0.478489,4.94];
ybs[4305]=['',2.9090138,-0.1956637,6.09];
ybs[4306]=['',2.9063459,-0.8641867,6.13];
ybs[4307]=['χ2 Hya',2.9105689,-0.4783886,5.71];
ybs[4308]=['',2.9107972,-0.8959539,6.3];
ybs[4309]=['65 Leo',2.9149629,0.0320018,5.52];
ybs[4310]=['',2.9135449,-0.5035228,6.77];
ybs[4311]=['',2.9123779,-0.8914896,6.32];
ybs[4312]=['64 Leo',2.9184657,0.404943,6.46];
ybs[4313]=['',2.91231,-1.0262048,6.02];
ybs[4314]=['',2.9156621,-0.5708838,6.59];
ybs[4315]=['',2.9123993,-1.0916353,4.61];
ybs[4316]=['',2.9116922,-1.1337943,6.41];
ybs[4317]=['',2.9161199,-0.746314,5.15];
ybs[4318]=['',2.9190401,-0.5287795,6.54];
ybs[4319]=['',2.91319,-1.239184,5.57];
ybs[4320]=['',2.9281113,1.1709054,6.06];
ybs[4321]=['',2.9206006,-0.5252556,6.49];
ybs[4322]=['67 Leo',2.9235213,0.4282359,5.68];
ybs[4323]=['',2.9258313,0.6315851,5.74];
ybs[4324]=['',2.9226673,-0.4922311,5.44];
ybs[4325]=['ψ UMa',2.9274505,0.7745123,3.01];
ybs[4326]=['',2.9273359,0.7519781,5.89];
ybs[4327]=['',2.9215029,-1.0314404,3.91];
ybs[4328]=['',2.921297,-1.0833154,5.13];
ybs[4329]=['',2.9276802,-0.5670547,5.81];
ybs[4330]=['',2.9391012,1.18943,6.4];
ybs[4331]=['',2.936115,0.2491932,6.3];
ybs[4332]=['',2.9316848,-1.0223744,6.88];
ybs[4333]=['β Crt',2.935503,-0.4005248,4.48];
ybs[4334]=['',2.9410478,0.9559427,6.63];
ybs[4335]=['',2.9398461,0.6229247,6.41];
ybs[4336]=['',2.9379763,-0.5682182,6.38];
ybs[4337]=['ψ Crt',2.9392463,-0.3250265,6.13];
ybs[4338]=['',2.9395234,-0.3817353,6.4];
ybs[4339]=['',2.9336061,-1.2489385,6.35];
ybs[4340]=['',2.9390776,-0.8591167,5.36];
ybs[4341]=['',2.9448485,0.7149886,6.33];
ybs[4342]=['',2.9390212,-1.0548796,4.6];
ybs[4343]=['',2.9408011,-0.870205,6.11];
ybs[4344]=['',2.9421878,-0.7765833,5.8];
ybs[4345]=['',2.939554,-1.1221138,5.23];
ybs[4346]=['69 Leo',2.9448633,-0.0033599,5.42];
ybs[4347]=['δ Leo',2.9465444,0.3560609,2.56];
ybs[4348]=['',2.9461029,0.1385397,5.79];
ybs[4349]=['θ Leo',2.9470764,0.2671507,3.34];
ybs[4350]=['',2.9438236,-0.9312105,5.76];
ybs[4351]=['',2.9430428,-1.042698,5.74];
ybs[4352]=['72 Leo',2.9513402,0.4009478,4.63];
ybs[4353]=['',2.9554682,0.9189162,6.5];
ybs[4354]=['',2.9494334,-0.7654502,6.21];
ybs[4355]=['73 Leo',2.9541463,0.2301128,5.32];
ybs[4356]=['',2.9545716,0.2220356,6.67];
ybs[4357]=['',2.9581577,0.8613774,5.88];
ybs[4358]=['φ Leo',2.9574989,-0.0658819,4.47];
ybs[4359]=['',2.9588198,-0.1266732,6.14];
ybs[4360]=['',2.956231,-0.802905,6.31];
ybs[4361]=['75 Leo',2.9602804,0.0329414,5.18];
ybs[4362]=['',2.9595456,-0.6656264,6.27];
ybs[4363]=['',2.9615766,-0.6084289,6.45];
ybs[4364]=['ξ UMa',2.9643969,0.5481368,4.87];
ybs[4365]=['ξ UMa',2.964404,0.5481368,4.41];
ybs[4366]=['',2.9618345,-0.6397964,6.68];
ybs[4367]=['ν UMa',2.9657044,0.5754507,3.48];
ybs[4368]=['',2.9649742,0.2070216,6.66];
ybs[4369]=['',2.9593977,-1.1858946,6.06];
ybs[4370]=['55 UMa',2.9686025,0.6643111,4.78];
ybs[4371]=['76 Leo',2.9673752,0.0266554,5.91];
ybs[4372]=['δ Crt',2.9691181,-0.2600883,3.56];
ybs[4373]=['',2.9768651,1.16897,6.21];
ybs[4374]=['',2.9681221,-1.1293299,5.99];
ybs[4375]=['',2.9637069,-1.3926308,6.35];
ybs[4376]=['σ Leo',2.9770891,0.1030779,4.05];
ybs[4377]=['',2.9689312,-1.3136371,6.27];
ybs[4378]=['',2.9805892,0.9939897,6.43];
ybs[4379]=['',2.9712175,-1.258694,6.41];
ybs[4380]=['π Cen',2.9759707,-0.9532048,3.89];
ybs[4381]=['',2.9852679,1.1206214,6.02];
ybs[4382]=['56 UMa',2.9847538,0.7567591,4.99];
ybs[4383]=['',2.9821494,-0.7813745,6.12];
ybs[4384]=['',2.9864909,0.0001388,6.05];
ybs[4385]=['λ Crt',2.9866571,-0.329932,5.09];
ybs[4386]=['',2.9858528,-0.6333525,5];
ybs[4387]=['',2.9789711,-1.3566777,6.43];
ybs[4388]=['',2.9852419,-0.9931471,5.79];
ybs[4389]=['ι Leo',2.9892798,0.1816085,3.94];
ybs[4390]=['79 Leo',2.9897199,0.0224101,5.39];
ybs[4391]=['',2.9860298,-1.1358377,5.11];
ybs[4392]=['ε Crt',2.9921425,-0.1916942,4.83];
ybs[4393]=['',2.9908439,-0.7468781,6.12];
ybs[4394]=['',2.9938938,0.1973344,5.8];
ybs[4395]=['γ Crt',2.9932866,-0.3108036,4.08];
ybs[4396]=['',2.989298,-1.263277,5.59];
ybs[4397]=['',2.9985278,0.972613,5.75];
ybs[4398]=['81 Leo',2.9966497,0.2850557,5.57];
ybs[4399]=['',2.9958105,-0.6315814,5.22];
ybs[4400]=['80 Leo',2.9975653,0.0652069,6.37];
ybs[4401]=['',2.9960651,-0.6609854,5.89];
ybs[4402]=['',3.0003324,0.5816586,6.32];
ybs[4403]=['',2.9963903,-1.1186982,5.17];
ybs[4404]=['83 Leo',3.0015817,0.0504236,6.5];
ybs[4405]=['',3.0002701,-1.0688266,5.3];
ybs[4406]=['κ Crt',3.003258,-0.2178291,5.94];
ybs[4407]=['',3.0012955,-0.9299811,5.81];
ybs[4408]=['τ Leo',3.0067363,0.0476828,4.95];
ybs[4409]=['',3.0065307,-0.0318362,6.25];
ybs[4410]=['',3.0066763,-0.6187663,6.45];
ybs[4411]=['',3.0122313,1.0760681,5.83];
ybs[4412]=['57 UMa',3.0119027,0.6843921,5.31];
ybs[4413]=['',3.0092794,-0.7469712,5.08];
ybs[4414]=['',3.0149537,0.9880882,6.28];
ybs[4415]=['',3.0073776,-1.2670837,6.09];
ybs[4416]=['85 Leo',3.014486,0.2668455,5.74];
ybs[4417]=['',3.0170489,0.9466215,6.41];
ybs[4418]=['',3.0140391,-0.4291335,5.76];
ybs[4419]=['',3.0254114,1.4137664,6.15];
ybs[4420]=['',3.0178342,0.8121582,6.35];
ybs[4421]=['58 UMa',3.0182435,0.7513479,5.94];
ybs[4422]=['87 Leo',3.0170843,-0.0545915,4.77];
ybs[4423]=['86 Leo',3.0179249,0.3191414,5.52];
ybs[4424]=['λ Dra',3.0225504,1.2078862,3.84];
ybs[4425]=['',3.0198716,0.8343524,6.42];
ybs[4426]=['',3.0211357,0.8493619,6.56];
ybs[4427]=['88 Leo',3.0234204,0.2485366,6.2];
ybs[4428]=['',3.0206874,-1.0716784,6.38];
ybs[4429]=['',3.0264283,1.0639196,5.48];
ybs[4430]=['',3.0234495,-0.3647916,6.24];
ybs[4431]=['ο1 Cen',3.0229929,-1.0396328,5.13];
ybs[4432]=['ο2 Cen',3.0231886,-1.0409176,5.15];
ybs[4433]=['',3.0254769,-0.5129124,5.81];
ybs[4434]=['',3.0254915,-0.5128736,5.64];
ybs[4435]=['',3.0260156,-0.4689884,6.16];
ybs[4436]=['',3.0278705,-0.1387871,5.95];
ybs[4437]=['',3.0277326,-0.7079196,5.64];
ybs[4438]=['',3.0252888,-1.1708822,5.9];
ybs[4439]=['',3.028233,-0.544746,5.04];
ybs[4440]=['ξ Hya',3.0286652,-0.5581948,3.54];
ybs[4441]=['',3.0298033,-0.2863213,6.05];
ybs[4442]=['',3.03309,0.6403799,6.4];
ybs[4443]=['',3.0313162,-0.7105482,5.39];
ybs[4444]=['',3.0339518,0.1902253,6.55];
ybs[4445]=['89 Leo',3.0347903,0.0512339,5.77];
ybs[4446]=['90 Leo',3.0363394,0.2909884,5.95];
ybs[4447]=['',3.0382247,0.9540095,5.63];
ybs[4448]=['',3.0351697,-0.5751891,5.98];
ybs[4449]=['',3.0379038,0.3545956,6.45];
ybs[4450]=['',3.0361738,-0.949262,4.62];
ybs[4451]=['2 Dra',3.0426944,1.2077357,5.2];
ybs[4452]=['',3.0370383,-0.8597704,5.5];
ybs[4453]=['',3.03825,-0.8289801,5.71];
ybs[4454]=['',3.0407407,0.1882602,6.56];
ybs[4455]=['',3.043322,0.4826968,5.8];
ybs[4456]=['',3.0413384,-0.8336786,5.25];
ybs[4457]=['λ Cen',3.0405003,-1.1020762,3.13];
ybs[4458]=['θ Crt',3.0448413,-0.1732566,4.7];
ybs[4459]=['',3.0443051,-0.5880824,5.74];
ybs[4460]=['',3.0447066,-0.6520973,6.31];
ybs[4461]=['υ Leo',3.0460392,-0.0165554,4.3];
ybs[4462]=['',3.0431276,-1.0677375,5.83];
ybs[4463]=['',3.0462163,-0.577926,6.29];
ybs[4464]=['',3.0503713,0.88128,6.14];
ybs[4465]=['',3.0459202,-1.0717717,5.15];
ybs[4466]=['',3.0485033,-0.8355225,5.44];
ybs[4467]=['59 UMa',3.0523201,0.7592326,5.59];
ybs[4468]=['',3.0513754,0.1528811,6.17];
ybs[4469]=['π Cha',3.0465547,-1.3268227,5.65];
ybs[4470]=['60 UMa',3.0532786,0.8152332,6.1];
ybs[4471]=['',3.0546159,1.1208885,6.46];
ybs[4472]=['',3.0531048,0.5846995,6.27];
ybs[4473]=['ω Vir',3.0526667,0.1397909,5.36];
ybs[4474]=['',3.0523751,-0.0446953,6.22];
ybs[4475]=['',3.0492848,-1.182373,5.96];
ybs[4476]=['',3.0540936,0.7851164,6.44];
ybs[4477]=['',3.0507815,-1.0812509,5.15];
ybs[4478]=['ι Crt',3.0535008,-0.2325947,5.48];
ybs[4479]=['',3.0549339,-0.4336424,6.42];
ybs[4480]=['',3.0586624,-0.2547032,6.21];
ybs[4481]=['',3.0586039,-0.2922569,6.19];
ybs[4482]=['',3.0567157,-1.1435845,5.17];
ybs[4483]=['',3.0616382,1.0095983,6.37];
ybs[4484]=['ο Hya',3.0601614,-0.6085884,4.7];
ybs[4485]=['92 Leo',3.0628523,0.3704973,5.26];
ybs[4486]=['61 UMa',3.0640581,0.5947525,5.33];
ybs[4487]=['',3.0622099,-0.944109,5.96];
ybs[4488]=['',3.0642345,-0.5117524,6.44];
ybs[4489]=['',3.0629209,-1.085854,4.94];
ybs[4490]=['',3.0671334,0.960762,6.27];
ybs[4491]=['62 UMa',3.0663206,0.5518946,5.73];
ybs[4492]=['',3.0649993,-0.7543436,5.55];
ybs[4493]=['',3.0668133,-0.5694069,5.22];
ybs[4494]=['3 Dra',3.0705237,1.1627397,5.3];
ybs[4495]=['',3.06853,0.3854722,6.59];
ybs[4496]=['',3.0682777,-0.3563752,6.22];
ybs[4497]=['',3.0623326,-1.4525477,6.33];
ybs[4498]=['',3.0743101,-0.6512737,5.98];
ybs[4499]=['',3.0712706,-1.3863381,6.39];
ybs[4500]=['',3.0764374,-0.1187208,6.07];
ybs[4501]=['',3.0744244,-1.0928276,5.03];
ybs[4502]=['',3.0778392,0.4379614,6.02];
ybs[4503]=['',3.0760043,-1.0996151,6.1];
ybs[4504]=['ζ Crt',3.0800954,-0.3224643,4.73];
ybs[4505]=['ξ Vir',3.0824375,0.141953,4.85];
ybs[4506]=['',3.081934,-0.8586102,6.26];
ybs[4507]=['ν Vir',3.0849417,0.1117779,4.03];
ybs[4508]=['χ UMa',3.0858961,0.8317261,3.71];
ybs[4509]=['',3.0842329,-0.7996232,5.29];
ybs[4510]=['λ Mus',3.0835136,-1.1668163,3.64];
ybs[4511]=['',3.0897548,0.9687146,5.27];
ybs[4512]=['',3.087553,-1.0699461,4.11];
ybs[4513]=['',3.0876941,-0.7090508,4.91];
ybs[4514]=['',3.090325,-0.6288775,6.17];
ybs[4515]=['',3.090975,-0.5307901,6.48];
ybs[4516]=['',3.0911068,-1.0091752,5.41];
ybs[4517]=['93 Leo',3.0942391,0.3507027,4.53];
ybs[4518]=['4 Vir',3.0939094,0.1417334,5.32];
ybs[4519]=['',3.0959542,-0.1821853,6.26];
ybs[4520]=['μ Mus',3.095054,-1.1683205,4.72];
ybs[4521]=['',3.0971049,0.2471219,5.88];
ybs[4522]=['',3.0974928,-0.4690546,5.11];
ybs[4523]=['',3.0987145,-0.0077448,6.15];
ybs[4524]=['β Leo',3.0989152,0.2521444,2.14];
ybs[4525]=['',3.0997396,0.2813059,6.04];
ybs[4526]=['',3.1017249,0.6074883,5.7];
ybs[4527]=['',3.1014271,-1.1155007,4.32];
ybs[4528]=['',3.1024908,-1.2278564,4.97];
ybs[4529]=['',3.1043845,-0.2790616,6.13];
ybs[4530]=['β Vir',3.1060261,0.0286155,3.61];
ybs[4531]=['',3.1048031,-1.0956237,5.7];
ybs[4532]=['',3.1056457,-0.4782716,6.48];
ybs[4533]=['',3.1070303,0.2121222,6.35];
ybs[4534]=['',3.1075069,-0.095269,5.64];
ybs[4535]=['',3.1080898,0.5803188,6.27];
ybs[4536]=['',3.1079074,-0.790613,4.46];
ybs[4537]=['',3.1089304,-0.2149018,6.35];
ybs[4538]=['',3.1103356,-0.5403573,5.85];
ybs[4539]=['',3.1109207,-1.1402464,4.9];
ybs[4540]=['',3.1160412,0.6561286,6.45];
ybs[4541]=['',3.1123454,-0.9968096,5.57];
ybs[4542]=['β Hya',3.1156465,-0.5939926,4.28];
ybs[4543]=['',3.1179899,-0.6142143,6.17];
ybs[4544]=['γ UMa',3.119775,0.9349641,2.44];
ybs[4545]=['',3.1197395,0.0074476,6.3];
ybs[4546]=['',3.1212023,-1.0041792,6.06];
ybs[4547]=['',3.1222836,-0.6610281,6.46];
ybs[4548]=['',3.1235134,-0.4509778,5.3];
ybs[4549]=['6 Vir',3.1250394,0.1451878,5.58];
ybs[4550]=['65 UMa',3.1252657,0.8089898,6.54];
ybs[4551]=['65 UMa',3.1256646,0.8088638,7.03];
ybs[4552]=['',3.1258629,0.6393341,6.49];
ybs[4553]=['',3.1247089,-1.1066108,5.91];
ybs[4554]=['95 Leo',3.127763,0.2708999,5.53];
ybs[4555]=['',3.1277044,-0.4992024,5.93];
ybs[4556]=['66 UMa',3.1291044,0.9856461,5.84];
ybs[4557]=['η Crt',3.1292265,-0.3015245,5.18];
ybs[4558]=['',3.1287577,-0.6948926,6.13];
ybs[4559]=['',3.1330844,1.0720495,6.22];
ybs[4560]=['',3.1323359,-0.8237562,6.26];
ybs[4561]=['',3.1337887,-0.5836474,6.21];
ybs[4562]=['',3.1346137,0.7019426,6.62];
ybs[4563]=['',3.1364183,-1.0921249,5.57];
ybs[4564]=['',3.1384292,0.5610994,6.42];
ybs[4565]=['',3.1394147,1.0705756,6.76];
ybs[4566]=['',3.1389863,-0.9851072,5.44];
ybs[4567]=['',3.139365,-0.71685,6.79];
ybs[4568]=['',3.141354,-1.1251165,5.61];
ybs[4569]=['',3.1418514,-0.4543817,6.43];
ybs[4570]=['',3.1425084,0.0070736,6.17];
ybs[4571]=['',3.1435417,0.5766958,5.96];
ybs[4572]=['',3.1430508,-0.9044633,6.05];
ybs[4573]=['ε Cha',3.1449828,-1.3674168,4.91];
ybs[4574]=['',3.1464173,0.5918365,6.5];
ybs[4575]=['7 Vir',3.1463986,0.0616103,5.37];
ybs[4576]=['',3.1479252,1.4089658,6.17];
ybs[4577]=['',3.1498629,-0.1845052,5.55];
ybs[4578]=['',3.14972,-0.3833176,6.28];
ybs[4579]=['π Vir',3.1504334,0.1132528,4.66];
ybs[4580]=['',3.1503526,-0.3452985,5.26];
ybs[4581]=['',3.1511191,-0.0330446,6.31];
ybs[4582]=['',3.1531276,-1.0058135,6.16];
ybs[4583]=['',3.1538448,0.6268644,5.59];
ybs[4584]=['67 UMa',3.1558216,0.7491006,5.21];
ybs[4585]=['',3.1571871,-1.4923925,6.05];
ybs[4586]=['',3.1575185,-1.2499026,6.42];
ybs[4587]=['',3.1581736,-1.2098181,5.89];
ybs[4588]=['',3.159108,-0.1362903,6.22];
ybs[4589]=['θ1 Cru',3.1598921,-1.1072024,4.33];
ybs[4590]=['',3.1626312,-0.7428018,5.15];
ybs[4591]=['',3.1630827,-1.2974626,6.44];
ybs[4592]=['2 Com',3.1652667,0.3723473,5.87];
ybs[4593]=['θ2 Cru',3.1655667,-1.1046327,4.72];
ybs[4594]=['',3.1670181,-1.1947546,5.35];
ybs[4595]=['κ Cha',3.1676742,-1.3376971,5.04];
ybs[4596]=['',3.1655148,1.4871714,6.27];
ybs[4597]=['',3.1683256,-1.0662983,5.96];
ybs[4598]=['ο Vir',3.1693416,0.150235,4.12];
ybs[4599]=['',3.1693024,1.3400744,5.8];
ybs[4600]=['',3.1712122,1.0962035,6.13];
ybs[4601]=['',3.1724461,-1.1462002,6.33];
ybs[4602]=['',3.1726092,-0.6251613,6.23];
ybs[4603]=['',3.1727932,-0.0568433,6.37];
ybs[4604]=['',3.1744138,-1.2003684,6.23];
ybs[4605]=['',3.1746321,-1.1490265,6.06];
ybs[4606]=['η Cru',3.176802,-1.1299054,4.15];
ybs[4607]=['',3.181099,-1.3175861,5.18];
ybs[4608]=['',3.1820103,-0.8863927,4.47];
ybs[4609]=['',3.1819815,-0.888172,6.37];
ybs[4610]=['',3.1826961,-0.8520339,5.34];
ybs[4611]=['δ Cen',3.1831991,-0.8874592,2.6];
ybs[4612]=['',3.1834743,-1.064169,6.22];
ybs[4613]=['α Crv',3.1833741,-0.4337851,4.02];
ybs[4614]=['',3.185532,-0.775821,5.75];
ybs[4615]=['',3.1855739,-0.7218079,5.48];
ybs[4616]=['10 Vir',3.1888946,0.0309384,5.95];
ybs[4617]=['',3.1889766,1.3009031,6.35];
ybs[4618]=['',3.1905098,-0.6079004,6.17];
ybs[4619]=['11 Vir',3.1904946,0.0991664,5.72];
ybs[4620]=['ε Crv',3.1908457,-0.3969725,3];
ybs[4621]=['',3.1927985,-0.6631447,6.06];
ybs[4622]=['3 Com',3.1925238,0.2911916,6.39];
ybs[4623]=['',3.1935545,0.4739665,6.01];
ybs[4624]=['',3.195181,-1.0716775,6.08];
ybs[4625]=['3 Crv',3.1949532,-0.4141247,5.46];
ybs[4626]=['',3.1949447,-0.7949604,6.61];
ybs[4627]=['',3.1970493,-0.8985746,6.23];
ybs[4628]=['ρ Cen',3.1976161,-0.9161878,3.96];
ybs[4629]=['',3.1938746,1.423925,6];
ybs[4630]=['4 Com',3.1982871,0.4493385,5.66];
ybs[4631]=['68 UMa',3.1977055,0.9936048,6.43];
ybs[4632]=['',3.1990064,0.4958661,6.49];
ybs[4633]=['5 Com',3.1996152,0.3563417,5.57];
ybs[4634]=['',3.200827,-1.1008821,5.92];
ybs[4635]=['',3.2027427,-1.2265649,6.17];
ybs[4636]=['',3.1993005,1.3524786,5.14];
ybs[4637]=['',3.2043823,-0.5977856,6.5];
ybs[4638]=['',3.2052956,-0.6816243,5.76];
ybs[4639]=['',3.20809,-1.37355,6.35];
ybs[4640]=['12 Vir',3.2052071,0.1769274,5.85];
ybs[4641]=['',3.2061073,-0.5919773,6.33];
ybs[4642]=['',3.2080434,-0.8002142,5.31];
ybs[4643]=['',3.2092298,-1.1263239,6.22];
ybs[4644]=['',3.2106711,0.9304304,6.16];
ybs[4645]=['',3.2121018,-0.3659805,5.83];
ybs[4646]=['δ Cru',3.2129549,-1.0275426,2.8];
ybs[4647]=['',3.2128726,-0.1821682,6.11];
ybs[4648]=['',3.2144303,-0.7337017,6.26];
ybs[4649]=['',3.2122788,1.22304,5.71];
ybs[4650]=['δ UMa',3.2137017,0.993224,3.31];
ybs[4651]=['',3.2155597,-0.4097781,6.54];
ybs[4652]=['γ Crv',3.2156435,-0.3083454,2.59];
ybs[4653]=['6 Com',3.2164103,0.2578542,5.1];
ybs[4654]=['',3.2186668,-1.2695462,6.22];
ybs[4655]=['',3.2145952,1.2640701,6.29];
ybs[4656]=['2 CVn',3.2168539,0.7074752,5.66];
ybs[4657]=['7 Com',3.217857,0.4157436,4.95];
ybs[4658]=['',3.2185209,0.5748499,5];
ybs[4659]=['',3.2216122,-1.148735,6.06];
ybs[4660]=['',3.2210914,-0.2935383,6.05];
ybs[4661]=['ε Mus',3.2237104,-1.1883197,4.11];
ybs[4662]=['',3.2227169,0.9261805,5.81];
ybs[4663]=['',3.222924,0.5028703,5.7];
ybs[4664]=['β Cha',3.2276362,-1.3864381,4.26];
ybs[4665]=['',3.2243686,-0.6321364,6.15];
ybs[4666]=['',3.2239754,0.2621363,6.34];
ybs[4667]=['',3.2258385,-0.0711971,6.99];
ybs[4668]=['',3.2258749,-0.0710953,6.54];
ybs[4669]=['ζ Cru',3.2274356,-1.1192428,4.04];
ybs[4670]=['',3.2273564,0.5257689,6.23];
ybs[4671]=['13 Vir',3.2280977,-0.0159182,5.9];
ybs[4672]=['',3.2297745,-0.9646062,5];
ybs[4673]=['',3.2175983,1.5035013,6.33];
ybs[4674]=['',3.2295894,0.4517431,6.48];
ybs[4675]=['8 Com',3.2308362,0.3998537,6.27];
ybs[4676]=['',3.2099683,1.5272412,6.28];
ybs[4677]=['',3.228101,1.3096207,5.38];
ybs[4678]=['9 Com',3.2315781,0.4892535,6.33];
ybs[4679]=['η Vir',3.2334861,-0.0138179,3.89];
ybs[4680]=['3 CVn',3.2328466,0.8527574,5.29];
ybs[4681]=['',3.2347498,-0.3892138,5.97];
ybs[4682]=['',3.2363764,-1.1513503,6.21];
ybs[4683]=['',3.2352239,0.4624198,5.54];
ybs[4684]=['',3.2350814,0.4516423,6.15];
ybs[4685]=['16 Vir',3.2354067,0.0556368,4.96];
ybs[4686]=['ζ Crv',3.2364241,-0.3899164,5.21];
ybs[4687]=['11 Com',3.2369531,0.3083658,4.74];
ybs[4688]=['',3.2367911,0.4700171,7.13];
ybs[4689]=['',3.2379911,-0.2389402,5.14];
ybs[4690]=['ε Cru',3.2401944,-1.0563745,3.59];
ybs[4691]=['70 UMa',3.2372616,1.0077387,5.55];
ybs[4692]=['',3.2427597,-0.9861002,5.92];
ybs[4693]=['ζ2 Mus',3.2436738,-1.1806557,5.15];
ybs[4694]=['ζ1 Mus',3.2440314,-1.1943662,5.74];
ybs[4695]=['',3.2433032,0.4302104,6.19];
ybs[4696]=['',3.2465563,-1.0088128,5.39];
ybs[4697]=['12 Com',3.2447148,0.4489245,4.81];
ybs[4698]=['17 Vir',3.2449244,0.0904242,6.4];
ybs[4699]=['',3.262265,-1.5023851,6.33];
ybs[4700]=['',3.2485384,-1.1825697,6.36];
ybs[4701]=['6 Crv',3.2486575,-0.4357238,5.68];
ybs[4702]=['',3.2497192,-0.6202437,5.32];
ybs[4703]=['',3.2498534,-0.6881418,6.4];
ybs[4704]=['',3.2504336,-0.6813009,5.79];
ybs[4705]=['4 CVn',3.2501866,0.7403376,6.06];
ybs[4706]=['5 CVn',3.2511592,0.8977568,4.8];
ybs[4707]=['13 Com',3.2525739,0.4533333,5.18];
ybs[4708]=['',3.2548079,-0.7244629,6.25];
ybs[4709]=['',3.2531722,0.4443304,6.42];
ybs[4710]=['',3.2575155,-1.150085,6.3];
ybs[4711]=['',3.2565508,-0.7441895,6.11];
ybs[4712]=['',3.256609,-0.20481,5.95];
ybs[4713]=['',3.2571756,-0.4864866,6.09];
ybs[4714]=['',3.2574675,-0.6162905,5.73];
ybs[4715]=['',3.2566968,0.415417,6.03];
ybs[4716]=['71 UMa',3.2555691,0.9887817,5.81];
ybs[4717]=['',3.2556785,1.1113959,6.32];
ybs[4718]=['6 CVn',3.2592059,0.6788314,5.02];
ybs[4719]=['',3.2628291,-1.1038664,4.86];
ybs[4720]=['α1 Cru',3.2631938,-1.103459,1.33];
ybs[4721]=['α2 Cru',3.2632376,-1.1034639,1.73];
ybs[4722]=['',3.2626937,-0.9001573,4.82];
ybs[4723]=['14 Com',3.2616915,0.4737511,4.95];
ybs[4724]=['',3.2638759,-0.8558693,6.26];
ybs[4725]=['',3.2639995,-0.5751622,5.55];
ybs[4726]=['',3.2667667,-1.1155009,6];
ybs[4727]=['γ Com',3.2640247,0.491205,4.36];
ybs[4728]=['16 Com',3.264251,0.4660238,5];
ybs[4729]=['',3.2669673,-1.0317735,5.5];
ybs[4730]=['',3.2610491,1.2532393,6.24];
ybs[4731]=['',3.2674553,0.1481081,6.37];
ybs[4732]=['',3.2681076,-0.2924517,6.35];
ybs[4733]=['σ Cen',3.2693002,-0.8788578,3.91];
ybs[4734]=['',3.2707469,-1.1251379,6.04];
ybs[4735]=['73 UMa',3.2665952,0.9702016,5.7];
ybs[4736]=['',3.2682087,-0.0827212,6.22];
ybs[4737]=['',3.2711687,-1.0806998,6.22];
ybs[4738]=['',3.2706577,-0.6835696,5.44];
ybs[4739]=['',3.2716453,-0.98667,6.15];
ybs[4740]=['',3.2714317,0.4555732,6.54];
ybs[4741]=['',3.2719058,0.4498573,6.65];
ybs[4742]=['17 Com',3.2726395,0.4500951,5.29];
ybs[4743]=['18 Com',3.2749896,0.418612,5.48];
ybs[4744]=['',3.2775205,-0.9887094,5.8];
ybs[4745]=['',3.2776313,-0.7305994,6.02];
ybs[4746]=['20 Com',3.2761927,0.3625387,5.69];
ybs[4747]=['δ Crv',3.2770213,-0.2904178,2.95];
ybs[4748]=['',3.2779429,-0.2359196,6.35];
ybs[4749]=['',3.278926,-0.4157513,5.63];
ybs[4750]=['74 UMa',3.2768576,1.0172071,5.35];
ybs[4751]=['7 CVn',3.2773696,0.8972983,6.21];
ybs[4752]=['75 UMa',3.2773582,1.0235196,6.08];
ybs[4753]=['γ Cru',3.283076,-0.998981,1.63];
ybs[4754]=['γ Cru',3.2835716,-0.9984184,6.42];
ybs[4755]=['4 Dra',3.2772476,1.2056204,4.95];
ybs[4756]=['21 Com',3.2817942,0.4266134,5.46];
ybs[4757]=['',3.2807714,0.9241968,6.21];
ybs[4758]=['',3.2853393,-1.0393071,5.48];
ybs[4759]=['',3.2880038,-1.2762832,5.88];
ybs[4760]=['',3.2834034,0.1305526,6.05];
ybs[4761]=['',3.2865684,-1.1105549,5.95];
ybs[4762]=['',3.2847292,-0.0903475,6.19];
ybs[4763]=['γ Mus',3.2892683,-1.2611227,3.87];
ybs[4764]=['',3.2867766,-0.5699827,6.46];
ybs[4765]=['η Crv',3.2866467,-0.2848396,4.31];
ybs[4766]=['',3.2889469,-0.2440515,5.74];
ybs[4767]=['20 Vir',3.2907672,0.1775285,6.26];
ybs[4768]=['',3.2923596,-0.3475969,6.26];
ybs[4769]=['',3.2931827,-0.2260926,5.58];
ybs[4770]=['22 Com',3.2929552,0.4216572,6.29];
ybs[4771]=['21 Vir',3.2940722,-0.1671293,5.48];
ybs[4772]=['',3.2953108,-0.8732455,6.38];
ybs[4773]=['',3.2932303,0.5781163,5.42];
ybs[4774]=['',3.2938465,0.5805115,6.24];
ybs[4775]=['β CVn',3.2935639,0.7196627,4.26];
ybs[4776]=['β Crv',3.2968088,-0.4105097,2.65];
ybs[4777]=['κ Dra',3.2918492,1.2158739,3.87];
ybs[4778]=['',3.298386,-0.7818571,5.77];
ybs[4779]=['23 Com',3.2985536,0.3927933,4.81];
ybs[4780]=['',3.3020926,-1.0815046,6.22];
ybs[4781]=['24 Com',3.2996867,0.3185832,6.56];
ybs[4782]=['24 Com',3.2997957,0.3185784,5.02];
ybs[4783]=['',3.299794,0.3797425,5.85];
ybs[4784]=['',3.3029507,-0.7181268,5.13];
ybs[4785]=['6 Dra',3.2972601,1.219953,4.94];
ybs[4786]=['',3.3040826,-0.6980212,5.8];
ybs[4787]=['',3.3037327,-0.3604261,6.2];
ybs[4788]=['α Mus',3.3098138,-1.2087996,2.69];
ybs[4789]=['25 Vir',3.3071877,-0.1039439,5.87];
ybs[4790]=['',3.3048053,1.036085,5.5];
ybs[4791]=['25 Com',3.3078346,0.2961101,5.68];
ybs[4792]=['τ Cen',3.3115491,-0.849358,3.86];
ybs[4793]=['',3.3113231,-0.4758188,5.45];
ybs[4794]=['',3.3193086,-1.3175981,6.49];
ybs[4795]=['',3.31273,0.0551352,6.33];
ybs[4796]=['',3.3171441,-1.1748939,6.25];
ybs[4797]=['',3.3140481,0.0302162,5.71];
ybs[4798]=['',3.3145673,0.1198148,7.08];
ybs[4799]=['',3.3157976,-0.3206815,6];
ybs[4800]=['',3.3172701,-0.5331215,5.89];
ybs[4801]=['9 CVn',3.3154755,0.7112394,6.37];
ybs[4802]=['',3.3167903,0.3933282,6.38];
ybs[4803]=['χ Vir',3.3179237,-0.1417021,4.66];
ybs[4804]=['',3.3217266,-1.1629996,6.26];
ybs[4805]=['26 Com',3.3171802,0.3654564,5.46];
ybs[4806]=['',3.3177495,0.6253265,6.45];
ybs[4807]=['',3.3209349,-0.7000657,4.64];
ybs[4808]=['',3.3276089,-0.8075415,5.84];
ybs[4809]=['γ Cen',3.3282353,-0.8566578,2.17];
ybs[4810]=['',3.3313384,-1.2135377,6.33];
ybs[4811]=['',3.3267764,-0.2292803,6.08];
ybs[4812]=['',3.326791,-0.2293045,5.98];
ybs[4813]=['',3.3303295,-1.0438629,4.93];
ybs[4814]=['27 Vir',3.3279491,0.1798254,6.19];
ybs[4815]=['γ Vir',3.328411,-0.0274468,3.65];
ybs[4816]=['γ Vir',3.328411,-0.0274468,3.68];
ybs[4817]=['',3.3292435,-0.3470017,6.03];
ybs[4818]=['ρ Vir',3.3293099,0.1764952,4.88];
ybs[4819]=['31 Vir',3.3296253,0.11665,5.59];
ybs[4820]=['',3.3343557,-1.1027274,5.31];
ybs[4821]=['',3.3329277,-0.854096,4.66];
ybs[4822]=['',3.3341023,-0.9786102,6.08];
ybs[4823]=['76 UMa',3.3272147,1.0923997,6.07];
ybs[4824]=['',3.3355303,-0.9826046,6];
ybs[4825]=['',3.3369933,-1.0301982,6.4];
ybs[4826]=['',3.3365133,-0.7033805,6.44];
ybs[4827]=['',3.33703,-0.0296686,5.93];
ybs[4828]=['',3.3388278,-0.6365576,6.39];
ybs[4829]=['',3.3388754,-0.4964901,5.48];
ybs[4830]=['',3.3338047,1.0652189,6.38];
ybs[4831]=['',3.3442202,-1.2034724,6.16];
ybs[4832]=['ι Cru',3.3465229,-1.066463,4.69];
ybs[4833]=['',3.340169,0.7675992,6.33];
ybs[4834]=['β Mus',3.3496773,-1.1908503,3.05];
ybs[4835]=['10 CVn',3.3425854,0.6834026,5.95];
ybs[4836]=['',3.3431027,0.7909395,4.99];
ybs[4837]=['32 Vir',3.3456134,0.1317829,5.22];
ybs[4838]=['',3.3496624,-0.9880575,4.65];
ybs[4839]=['33 Vir',3.3488993,0.1643638,5.67];
ybs[4840]=['',3.3509781,-0.5836058,5.86];
ybs[4841]=['27 Com',3.350019,0.287192,5.12];
ybs[4842]=['',3.337952,1.404959,6.4];
ybs[4843]=['β Cru',3.3556277,-1.0439005,1.25];
ybs[4844]=['',3.3518213,0.1017223,6.34];
ybs[4845]=['34 Vir',3.3525943,0.2065685,6.07];
ybs[4846]=['',3.3541804,-0.1121279,6.26];
ybs[4847]=['',3.355819,-0.4358809,6.44];
ybs[4848]=['35 Vir',3.3554104,0.0602191,6.41];
ybs[4849]=['',3.3521974,1.0935933,5.89];
ybs[4850]=['',3.3582271,-0.4838036,5.66];
ybs[4851]=['28 Com',3.356993,0.2344085,6.56];
ybs[4852]=['',3.3651671,-1.2585329,5.55];
ybs[4853]=['7 Dra',3.3531447,1.1635719,5.43];
ybs[4854]=['',3.3592663,0.4314088,6.31];
ybs[4855]=['29 Com',3.3598882,0.2443486,5.7];
ybs[4856]=['11 CVn',3.3585788,0.8437718,6.27];
ybs[4857]=['',3.3581149,1.0506464,5.85];
ybs[4858]=['',3.3665127,-1.056326,6.75];
ybs[4859]=['30 Com',3.3614472,0.4787422,5.78];
ybs[4860]=['ι Oct',3.3928133,-1.4824172,5.46];
ybs[4861]=['',3.3667613,-0.8479141,6.24];
ybs[4862]=['',3.3696455,-0.9234467,5.73];
ybs[4863]=['',3.3658537,0.3969077,6.43];
ybs[4864]=['',3.3681148,-0.5955339,4.91];
ybs[4865]=['',3.365201,0.6526612,5.89];
ybs[4866]=['',3.3713133,-1.0550826,5.72];
ybs[4867]=['',3.3709094,-0.1825683,6.41];
ybs[4868]=['37 Vir',3.3718115,0.0512191,6.02];
ybs[4869]=['',3.3736977,-0.6946902,5.98];
ybs[4870]=['',3.3744594,-0.8415302,6.33];
ybs[4871]=['',3.3736155,-0.4687961,6.15];
ybs[4872]=['',3.3759939,-0.9416289,6.24];
ybs[4873]=['31 Com',3.371944,0.4785437,4.94];
ybs[4874]=['32 Com',3.374259,0.295867,6.32];
ybs[4875]=['',3.3788897,-0.9612286,5.93];
ybs[4876]=['',3.37538,0.2792627,6.3];
ybs[4877]=['',3.3803646,-1.0550587,5.76];
ybs[4878]=['',3.378967,-0.8563487,4.33];
ybs[4879]=['',3.3802182,-0.7033796,4.27];
ybs[4880]=['κ Cru',3.3823591,-1.0559012,5.9];
ybs[4881]=['38 Vir',3.3787271,-0.0641389,6.11];
ybs[4882]=['',3.3569061,1.4537837,5.85];
ybs[4883]=['',3.357409,1.4536918,5.28];
ybs[4884]=['35 Com',3.3789777,0.368669,4.9];
ybs[4885]=['',3.384692,-1.0219291,6.58];
ybs[4886]=['',3.3806896,-0.0758462,6.44];
ybs[4887]=['λ Cru',3.3859716,-1.034427,4.62];
ybs[4888]=['μ1 Cru',3.385644,-1.0000635,4.03];
ybs[4889]=['μ2 Cru',3.3857313,-0.9998986,5.17];
ybs[4890]=['41 Vir',3.3813822,0.2146207,6.25];
ybs[4891]=['',3.3837112,-0.2054304,6];
ybs[4892]=['ψ Vir',3.3838733,-0.1686087,4.79];
ybs[4893]=['',3.3870091,-0.772719,5.89];
ybs[4894]=['',3.3828544,0.5831623,6.26];
ybs[4895]=['ε UMa',3.3816066,0.9745567,1.77];
ybs[4896]=['',3.3885104,-0.751144,5.47];
ybs[4897]=['',3.3949505,-1.2619891,5.93];
ybs[4898]=['',3.3915724,-0.9940972,5.32];
ybs[4899]=['',3.3857882,0.8216146,5.84];
ybs[4900]=['δ Vir',3.3892088,0.0571767,3.38];
ybs[4901]=['',3.39063,-0.2696258,6.17];
ybs[4902]=['',3.3934285,-0.4639377,6.62];
ybs[4903]=['',3.3963329,-0.8957016,5.16];
ybs[4904]=['α1 CVn',3.3905786,0.666598,5.6];
ybs[4905]=['α2 CVn',3.390673,0.6666611,2.9];
ybs[4906]=['8 Dra',3.3875273,1.1399977,5.24];
ybs[4907]=['',3.3915114,0.9420939,5.82];
ybs[4908]=['',3.397979,-0.3992465,6.31];
ybs[4909]=['',3.395328,0.8038222,6.12];
ybs[4910]=['36 Com',3.4035507,0.3017391,4.78];
ybs[4911]=['44 Vir',3.4069714,-0.0686422,5.79];
ybs[4912]=['',3.4111771,-0.5868863,6.02];
ybs[4913]=['δ Mus',3.4201148,-1.2508677,3.62];
ybs[4914]=['37 Com',3.409283,0.5351899,4.9];
ybs[4915]=['46 Vir',3.4110691,-0.0609022,5.99];
ybs[4916]=['',3.4110545,0.3185615,6.2];
ybs[4917]=['',3.4009939,1.3151297,6.01];
ybs[4918]=['9 Dra',3.4067565,1.1602299,5.32];
ybs[4919]=['38 Com',3.4133067,0.2967462,5.96];
ybs[4920]=['',3.4236932,-1.2495954,6.03];
ybs[4921]=['78 UMa',3.4107437,0.9816703,4.93];
ybs[4922]=['ε Vir',3.4178074,0.1891687,2.83];
ybs[4923]=['ξ1 Cen',3.4246382,-0.8665141,4.85];
ybs[4924]=['',3.415012,1.1081026,6];
ybs[4925]=['',3.4250945,-0.3613427,5.58];
ybs[4926]=['',3.4190726,1.040139,6.53];
ybs[4927]=['48 Vir',3.4255093,-0.0660375,6.59];
ybs[4928]=['',3.4299152,-0.7211152,6.26];
ybs[4929]=['',3.4332822,-0.911674,6.43];
ybs[4930]=['',3.4365302,-0.8479433,4.71];
ybs[4931]=['',3.4377184,-0.7279511,5.59];
ybs[4932]=['ξ2 Cen',3.4393305,-0.8731179,4.27];
ybs[4933]=['14 CVn',3.433028,0.622713,5.25];
ybs[4934]=['',3.4418322,-1.0468542,5.99];
ybs[4935]=['',3.4334148,0.7879911,5.63];
ybs[4936]=['39 Com',3.4359027,0.3671015,5.99];
ybs[4937]=['',3.4390125,-0.628001,6.54];
ybs[4938]=['',3.4349981,0.5045651,6.54];
ybs[4939]=['40 Com',3.4359858,0.3926319,5.6];
ybs[4940]=['',3.4274894,1.2724331,6.31];
ybs[4941]=['',3.4426207,-0.935138,5.71];
ybs[4942]=['θ Mus',3.4452504,-1.1418998,5.51];
ybs[4943]=['',3.4350735,1.0807422,6.14];
ybs[4944]=['41 Com',3.4394128,0.4800509,4.8];
ybs[4945]=['49 Vir',3.4430022,-0.1895426,5.19];
ybs[4946]=['',3.4425302,0.4788505,6.19];
ybs[4947]=['',3.4457972,-0.1588957,5.55];
ybs[4948]=['ψ Hya',3.4482138,-0.4055722,4.95];
ybs[4949]=['',3.4488517,-0.1685609,6.32];
ybs[4950]=['',3.4484803,0.1728351,5.78];
ybs[4951]=['50 Vir',3.4511091,-0.1823669,5.94];
ybs[4952]=['',3.4509705,0.2919796,5.91];
ybs[4953]=['θ Vir',3.4519038,-0.0987553,4.38];
ybs[4954]=['',3.4500049,0.6510709,6.02];
ybs[4955]=['',3.4571843,-0.9195466,6.06];
ybs[4956]=['',3.4620454,-1.2227945,5.91];
ybs[4957]=['15 CVn',3.450224,0.6704588,6.28];
ybs[4958]=['α Com',3.4517977,0.303863,5.22];
ybs[4959]=['α Com',3.4517977,0.303863,5.22];
ybs[4960]=['',3.4576731,-0.7391858,5.79];
ybs[4961]=['17 CVn',3.4517638,0.669849,5.91];
ybs[4962]=['',3.4616399,-1.1069193,6.33];
ybs[4963]=['',3.4587473,-0.759009,5.25];
ybs[4964]=['',3.4500944,1.0840195,6.54];
ybs[4965]=['',3.463234,-1.0478921,4.6];
ybs[4966]=['',3.4742288,-1.3712317,5.85];
ybs[4967]=['',3.4659102,-1.1579528,5.9];
ybs[4968]=['',3.4596107,-0.4654925,6.5];
ybs[4969]=['',3.4615464,-0.661865,4.85];
ybs[4970]=['',3.4660462,-1.0460721,6.16];
ybs[4971]=['53 Vir',3.4612369,-0.2847964,5.04];
ybs[4972]=['',3.4651225,-0.7473255,6.22];
ybs[4973]=['β Com',3.4598725,0.4844859,4.26];
ybs[4974]=['',3.4610882,0.4213057,6.33];
ybs[4975]=['',3.467706,-0.886955,5.89];
ybs[4976]=['',3.4630368,0.1996162,5.77];
ybs[4977]=['',3.4631622,0.3252025,6.53];
ybs[4978]=['',3.4715467,-1.0262976,5.89];
ybs[4979]=['',3.4717631,-1.0336182,4.92];
ybs[4980]=['54 Vir',3.4673334,-0.3306604,6.28];
ybs[4981]=['',3.4699851,-0.754987,6.16];
ybs[4982]=['',3.4658151,0.3247728,6.11];
ybs[4983]=['η Mus',3.4767153,-1.1870487,4.8];
ybs[4984]=['',3.4776763,-1.218207,6.37];
ybs[4985]=['55 Vir',3.4705579,-0.3499296,5.33];
ybs[4986]=['',3.4734575,-0.856524,5.89];
ybs[4987]=['',3.4676883,0.6987256,4.92];
ybs[4988]=['',3.4716463,0.1957049,5.67];
ybs[4989]=['',3.475103,-0.6368634,6.19];
ybs[4990]=['',3.4830645,-1.1389408,6.07];
ybs[4991]=['57 Vir',3.4784107,-0.3501373,5.22];
ybs[4992]=['',3.4852516,-1.1676548,4.87];
ybs[4993]=['',3.4652456,1.2685066,6.59];
ybs[4994]=['19 CVn',3.4755861,0.7109921,5.79];
ybs[4995]=['',3.4801096,-0.0263338,6.68];
ybs[4996]=['',3.4825442,-0.5519478,5.1];
ybs[4997]=['',3.4790294,0.3304496,6.45];
ybs[4998]=['',3.4843159,-0.7696473,5.84];
ybs[4999]=['',3.4586026,1.402413,6.25];
ybs[5000]=['',3.4803271,0.3432544,6.45];
ybs[5001]=['59 Vir',3.4814992,0.1624197,5.22];
ybs[5002]=['',3.49502,-1.2593115,6.04];
ybs[5003]=['',3.4835584,0.2366221,5.33];
ybs[5004]=['',3.484784,-0.0138707,6.37];
ybs[5005]=['σ Vir',3.4851708,0.0934044,4.8];
ybs[5006]=['',3.4904063,-0.8971683,6.19];
ybs[5007]=['20 CVn',3.4843398,0.706063,4.73];
ybs[5008]=['',3.4785046,1.1918813,6.2];
ybs[5009]=['61 Vir',3.4889761,-0.3216516,4.74];
ybs[5010]=['γ Hya',3.4913043,-0.4064778,3];
ybs[5011]=['',3.4906326,0.0623076,6.62];
ybs[5012]=['',3.4884863,0.5930656,5.82];
ybs[5013]=['21 CVn',3.4871522,0.8650549,5.15];
ybs[5014]=['',3.4995348,-1.0452916,6.18];
ybs[5015]=['',3.49111,0.6110444,6.02];
ybs[5016]=['',3.4994401,-0.9226773,5.48];
ybs[5017]=['',3.5003254,-0.9759528,6.02];
ybs[5018]=['ι Cen',3.4988608,-0.6427995,2.75];
ybs[5019]=['',3.5007021,-0.820269,5.77];
ybs[5020]=['',3.5106756,-1.2612389,6.05];
ybs[5021]=['',3.4986772,0.0492916,6.26];
ybs[5022]=['23 CVn',3.4964363,0.6987077,5.6];
ybs[5023]=['',3.5025001,-0.3421928,6.21];
ybs[5024]=['',3.5084175,-1.0662093,6.18];
ybs[5025]=['',3.508579,-1.0664904,4.53];
ybs[5026]=['',3.5065946,-0.9128107,5.83];
ybs[5027]=['',3.5030516,0.034382,5.69];
ybs[5028]=['',3.5091151,-0.8388068,6.16];
ybs[5029]=['',3.5098559,-0.8496224,6.38];
ybs[5030]=['64 Vir',3.5050535,0.0879215,5.87];
ybs[5031]=['',3.5148441,-1.1284012,4.53];
ybs[5032]=['ι1 Mus',3.5210035,-1.3090723,5.05];
ybs[5033]=['',3.5099253,-0.5813166,6.22];
ybs[5034]=['63 Vir',3.5091106,-0.3115813,5.37];
ybs[5035]=['',3.5039531,0.764207,6.35];
ybs[5036]=['',3.5135192,-0.8716155,6.48];
ybs[5037]=['65 Vir',3.5102229,-0.0879892,5.89];
ybs[5038]=['',3.5202009,-1.1275145,5.31];
ybs[5039]=['',3.5234536,-1.2347142,5.67];
ybs[5040]=['66 Vir',3.5156305,-0.092164,5.75];
ybs[5041]=['ι2 Mus',3.5305757,-1.3056467,6.63];
ybs[5042]=['',3.5121146,0.6443237,6.07];
ybs[5043]=['',3.5151876,0.214941,6.44];
ybs[5044]=['ζ UMa',3.5116981,0.9565871,2.27];
ybs[5045]=['ζ UMa',3.5117636,0.9565242,3.95];
ybs[5046]=['α Vir',3.518509,-0.1968379,0.98];
ybs[5047]=['',3.5176422,0.4143032,5.78];
ybs[5048]=['',3.5231111,-0.6958919,5.09];
ybs[5049]=['',3.5227169,-0.0228445,5.97];
ybs[5050]=['',3.5266727,-0.7263062,5.69];
ybs[5051]=['',3.5276384,-0.8597505,6.31];
ybs[5052]=['80 UMa',3.5173514,0.9576873,4.01];
ybs[5053]=['',3.5287033,-0.863885,6.28];
ybs[5054]=['68 Vir',3.5251959,-0.2238221,5.25];
ybs[5055]=['',3.5279944,-0.703005,6.4];
ybs[5056]=['',3.5362218,-1.2172601,6.2];
ybs[5057]=['',3.5222415,0.8013096,5.88];
ybs[5058]=['69 Vir',3.528447,-0.280819,4.76];
ybs[5059]=['',3.5372837,-1.1308264,6.11];
ybs[5060]=['',3.5202696,1.1020819,6.5];
ybs[5061]=['',3.5378354,-0.895022,5.06];
ybs[5062]=['70 Vir',3.5322709,0.2384635,4.98];
ybs[5063]=['',3.5199318,1.2614354,5.79];
ybs[5064]=['',3.5248474,1.1278196,6.66];
ybs[5065]=['',3.5252903,1.1275388,7.04];
ybs[5066]=['',3.5294807,0.918563,6.34];
ybs[5067]=['',3.5317913,0.7088441,6.47];
ybs[5068]=['',3.5360645,-0.0258345,6.43];
ybs[5069]=['',3.530446,0.8808889,6.8];
ybs[5070]=['',3.5384349,-0.4083555,4.97];
ybs[5071]=['71 Vir',3.5357469,0.1867948,5.65];
ybs[5072]=['',3.5574744,-1.3558267,6.48];
ybs[5073]=['',3.5329079,0.8831792,6.43];
ybs[5074]=['κ Oct',3.6000594,-1.4951665,5.58];
ybs[5075]=['',3.5311423,1.0442282,5.4];
ybs[5076]=['',3.5392256,0.1232774,6.17];
ybs[5077]=['',3.5390612,0.1029346,6.51];
ybs[5078]=['72 Vir',3.5412892,-0.1149438,6.09];
ybs[5079]=['',3.5445843,-0.6898042,3.88];
ybs[5080]=['',3.5465668,-0.4926723,6.47];
ybs[5081]=['',3.5219782,1.3705644,5.77];
ybs[5082]=['',3.5491192,-0.6722016,6.16];
ybs[5083]=['',3.5569471,-1.1475066,6.37];
ybs[5084]=['73 Vir',3.5485463,-0.3288908,6.01];
ybs[5085]=['74 Vir',3.5479925,-0.1111953,4.69];
ybs[5086]=['',3.5440673,0.7328769,6.08];
ybs[5087]=['',3.5511354,-0.5027912,5.69];
ybs[5088]=['',3.5510508,-0.5180193,6.45];
ybs[5089]=['75 Vir',3.552049,-0.2701429,5.55];
ybs[5090]=['76 Vir',3.5524318,-0.1794193,5.21];
ybs[5091]=['',3.5525748,-0.1275828,6.68];
ybs[5092]=['',3.5511573,0.4229171,6.11];
ybs[5093]=['',3.5598694,-0.8445096,6.33];
ybs[5094]=['',3.5605294,-0.5833834,6.44];
ybs[5095]=['78 Vir',3.557296,0.0618575,4.94];
ybs[5096]=['',3.5599303,-0.2326356,5.91];
ybs[5097]=['ζ Vir',3.5598121,-0.0123994,3.37];
ybs[5098]=['',3.5576569,0.6749971,6.37];
ybs[5099]=['81 UMa',3.5560421,0.9640127,5.6];
ybs[5100]=['',3.5595843,0.646957,4.98];
ybs[5101]=['80 Vir',3.5635013,-0.0961767,5.73];
ybs[5102]=['24 CVn',3.5577514,0.8534911,4.7];
ybs[5103]=['',3.5724407,-1.0787171,5.63];
ybs[5104]=['',3.5634024,0.1761092,6.49];
ybs[5105]=['',3.5831009,-1.3229137,6.34];
ybs[5106]=['',3.5613049,0.7693838,6.84];
ybs[5107]=['',3.5697965,-0.6035676,6.5];
ybs[5108]=['',3.5711828,-0.7724368,5.98];
ybs[5109]=['',3.580154,-1.2314801,6.1];
ybs[5110]=['',3.5694737,-0.4644166,5.78];
ybs[5111]=['',3.5725459,-0.8123164,5.9];
ybs[5112]=['',3.5762729,-1.02152,6.42];
ybs[5113]=['',3.5693889,0.4275924,5.74];
ybs[5114]=['',3.5792552,-1.0076952,6.01];
ybs[5115]=['',3.5856825,-1.2374673,6.59];
ybs[5116]=['',3.5673478,0.8617124,6.49];
ybs[5117]=['25 CVn',3.5712096,0.6314777,4.82];
ybs[5118]=['',3.5778063,-0.517918,5.83];
ybs[5119]=['',3.5745729,0.2476245,6.52];
ybs[5120]=['',3.585716,-1.1290577,5.79];
ybs[5121]=['',3.5562294,1.3339898,6.57];
ybs[5122]=['ε Cen',3.583719,-0.9351435,2.3];
ybs[5123]=['',3.571886,0.8831501,6.48];
ybs[5124]=['',3.5840434,-0.8737754,6];
ybs[5125]=['',3.5823329,-0.6957144,6.27];
ybs[5126]=['',3.5829094,-0.7010178,5.6];
ybs[5127]=['',3.5784579,0.3168061,6.48];
ybs[5128]=['',3.5809349,0.1855741,5.57];
ybs[5129]=['',3.5680195,1.2414197,5.5];
ybs[5130]=['',3.5933311,-1.0280007,5.38];
ybs[5131]=['',3.5919017,-0.9542229,5.01];
ybs[5132]=['82 UMa',3.5795847,0.9216709,5.46];
ybs[5133]=['',3.5835329,0.5392822,6.21];
ybs[5134]=['1 Boo',3.5855637,0.3463136,5.75];
ybs[5135]=['',3.5853033,0.4878548,6.23];
ybs[5136]=['',3.589977,-0.4112476,6.59];
ybs[5137]=['',3.5912724,-0.5883489,6.05];
ybs[5138]=['',3.5835157,0.8797526,6.32];
ybs[5139]=['2 Boo',3.5871002,0.3906514,5.62];
ybs[5140]=['82 Vir',3.5901415,-0.1538694,5.01];
ybs[5141]=['',3.5972248,-0.9927558,6];
ybs[5142]=['',3.596827,-0.888429,6.41];
ybs[5143]=['',3.5830644,0.9964809,6.29];
ybs[5144]=['83 UMa',3.5848628,0.9523984,4.66];
ybs[5145]=['',3.5965266,-0.7245524,5.98];
ybs[5146]=['',3.5924681,0.1444339,6.16];
ybs[5147]=['',3.5998249,-0.73618,5.98];
ybs[5148]=['',3.6027631,-0.8923067,6.47];
ybs[5149]=['84 Vir',3.5962576,0.0597842,5.36];
ybs[5150]=['',3.5929342,0.7253821,6.3];
ybs[5151]=['',3.5941747,0.6087031,5.98];
ybs[5152]=['',3.5875002,1.1293922,5.85];
ybs[5153]=['',3.6000891,-0.0979366,6.51];
ybs[5154]=['',3.5989329,0.3942308,6.13];
ybs[5155]=['83 Vir',3.602859,-0.28434,5.6];
ybs[5156]=['',3.6041974,-0.4470326,6.21];
ybs[5157]=['',3.6079372,-0.4577676,5.81];
ybs[5158]=['1 Cen',3.6084072,-0.5786798,4.23];
ybs[5159]=['',3.5987937,0.9067325,6.02];
ybs[5160]=['85 Vir',3.6076028,-0.2771505,6.19];
ybs[5161]=['',3.6161963,-1.0943497,6.51];
ybs[5162]=['',3.6132352,-0.899622,4.65];
ybs[5163]=['86 Vir',3.6090845,-0.2188405,5.51];
ybs[5164]=['',3.6139649,-0.6346655,5.15];
ybs[5165]=['',3.6167047,-0.8789655,5.91];
ybs[5166]=['',3.617495,-0.8802155,5.45];
ybs[5167]=['',3.6043012,0.9733222,6.5];
ybs[5168]=['',3.6146523,-0.1714056,6.05];
ybs[5169]=['',3.6092642,0.7151779,5.87];
ybs[5170]=['',3.6097384,0.6700665,5.94];
ybs[5171]=['87 Vir',3.6156681,-0.3136636,5.43];
ybs[5172]=['3 Boo',3.6118121,0.4466372,5.95];
ybs[5173]=['',3.61318,0.108888,6.33];
ybs[5174]=['',3.5900811,1.3605113,5.91];
ybs[5175]=['τ Boo',3.6143312,0.3027275,4.5];
ybs[5176]=['',3.6126998,0.6707482,5.5];
ybs[5177]=['84 UMa',3.6103594,0.948079,5.7];
ybs[5178]=['',3.6595958,-1.4447025,5.95];
ybs[5179]=['',3.6225947,-0.6250915,6.53];
ybs[5180]=['ν Cen',3.6253323,-0.7295274,3.41];
ybs[5181]=['η UMa',3.6147041,0.858732,1.86];
ybs[5182]=['2 Cen',3.6248629,-0.6032193,4.19];
ybs[5183]=['μ Cen',3.6258456,-0.7432471,3.04];
ybs[5184]=['',3.6371544,-1.21321,5.75];
ybs[5185]=['',3.6200722,0.54243,5.62];
ybs[5186]=['89 Vir',3.626363,-0.318438,4.97];
ybs[5187]=['',3.6276289,-0.5095019,6.18];
ybs[5188]=['',3.6288536,-0.6983406,6.44];
ybs[5189]=['',3.6211952,0.68821,7.4];
ybs[5190]=['υ Boo',3.6240218,0.2737841,4.07];
ybs[5191]=['6 Boo',3.6249466,0.3691915,4.91];
ybs[5192]=['',3.629469,-0.3492059,6.53];
ybs[5193]=['',3.586012,1.4423308,5.98];
ybs[5194]=['',3.6247485,0.6374244,6.38];
ybs[5195]=['',3.6282892,0.0940097,6.01];
ybs[5196]=['',3.6354891,-0.8204731,5.77];
ybs[5197]=['',3.6370302,-0.9236642,5.25];
ybs[5198]=['',3.6343757,-0.6378107,6.35];
ybs[5199]=['',3.6329033,-0.4276308,6.45];
ybs[5200]=['3 Cen',3.6352307,-0.5777899,4.56];
ybs[5201]=['3 Cen',3.6352671,-0.5777947,6.06];
ybs[5202]=['',3.6360188,-0.5537908,6.12];
ybs[5203]=['',3.6236794,1.0712497,5.96];
ybs[5204]=['',3.6305719,0.6049623,6.65];
ybs[5205]=['',3.6309161,0.6030768,5.87];
ybs[5206]=['',3.6269463,1.0197705,6.46];
ybs[5207]=['',3.6442255,-0.9334644,5.89];
ybs[5208]=['',3.6501783,-1.1826725,5.71];
ybs[5209]=['',3.6337034,0.599235,4.74];
ybs[5210]=['',3.6364329,0.2103977,6.04];
ybs[5211]=['4 Cen',3.6412396,-0.5591668,4.73];
ybs[5212]=['',3.6428166,-0.6243775,5.54];
ybs[5213]=['',3.6449496,-0.824458,6.1];
ybs[5214]=['',3.6442197,-0.6182722,6.19];
ybs[5215]=['7 Boo',3.6402822,0.3110636,5.7];
ybs[5216]=['10 Dra',3.6306728,1.1277039,4.65];
ybs[5217]=['',3.6283365,1.1903929,6.4];
ybs[5218]=['',3.645824,-0.5005528,6.04];
ybs[5219]=['',3.6398618,0.4980801,5.9];
ybs[5220]=['',3.6506603,-0.9122955,5.71];
ybs[5221]=['ζ Cen',3.6519153,-0.8272481,2.55];
ybs[5222]=['90 Vir',3.6471335,-0.0281488,5.15];
ybs[5223]=['',3.6484262,-0.1425683,6.19];
ybs[5224]=['',3.655691,-0.9466878,6.14];
ybs[5225]=['η Boo',3.6466829,0.3191861,2.68];
ybs[5226]=['',3.6566853,-0.9566788,6];
ybs[5227]=['',3.6522901,-0.5479366,6.51];
ybs[5228]=['86 UMa',3.641967,0.9358211,5.7];
ybs[5229]=['',3.6553197,-0.8150999,5.83];
ybs[5230]=['',3.6781234,-1.3735392,6.09];
ybs[5231]=['',3.6621522,-1.1134426,4.71];
ybs[5232]=['',3.6661925,-1.1503326,6.2];
ybs[5233]=['',3.6517775,0.2434199,6.16];
ybs[5234]=['92 Vir',3.6547803,0.0164284,5.91];
ybs[5235]=['',3.6528682,0.5571636,6.32];
ybs[5236]=['',3.6596066,-0.4037257,6.14];
ybs[5237]=['9 Boo',3.654702,0.4779178,5.01];
ybs[5238]=['φ Cen',3.6636638,-0.7366965,3.83];
ybs[5239]=['υ1 Cen',3.6655462,-0.7838668,3.87];
ybs[5240]=['47 Hya',3.6642652,-0.4377449,5.15];
ybs[5241]=['',3.6684406,-0.8810156,5.91];
ybs[5242]=['',3.6735142,-1.0749408,6.49];
ybs[5243]=['',3.6765448,-1.1584905,5.97];
ybs[5244]=['',3.6641152,0.2537839,6];
ybs[5245]=['10 Boo',3.6638944,0.3767713,5.76];
ybs[5246]=['',3.6574916,1.071348,6.37];
ybs[5247]=['48 Hya',3.6707441,-0.438402,5.77];
ybs[5248]=['',3.6695141,-0.0638457,6.4];
ybs[5249]=['',3.6769085,-0.7038939,6.13];
ybs[5250]=['υ2 Cen',3.6788816,-0.7978149,4.34];
ybs[5251]=['θ Aps',3.6983854,-1.3422165,5.5];
ybs[5252]=['',3.6758976,0.153358,5.99];
ybs[5253]=['11 Boo',3.6747779,0.4761022,6.23];
ybs[5254]=['τ Vir',3.6773792,0.0250729,4.26];
ybs[5255]=['',3.6811894,-0.4806226,5.48];
ybs[5256]=['',3.6869073,-0.9829856,5.92];
ybs[5257]=['β Cen',3.688899,-1.0555795,0.61];
ybs[5258]=['',3.6841316,-0.5548637,6.18];
ybs[5259]=['',3.6863095,-0.7248469,6.11];
ybs[5260]=['',3.6810725,0.1671809,6.2];
ybs[5261]=['',3.6786945,0.7966706,6.27];
ybs[5262]=['',3.6876422,-0.3932034,6.3];
ybs[5263]=['',3.6854485,0.1863893,6.3];
ybs[5264]=['',3.6858411,0.1298362,6.26];
ybs[5265]=['',3.6872752,0.0836642,6.24];
ybs[5266]=['',3.688846,-0.0957929,6.39];
ybs[5267]=['',3.6899445,-0.2631737,6.28];
ybs[5268]=['',3.6970422,-0.9560232,6.17];
ybs[5269]=['',3.7114409,-1.3082301,6.02];
ybs[5270]=['',3.6819108,0.8877514,6.15];
ybs[5271]=['',3.7002057,-1.0440911,6.42];
ybs[5272]=['',3.6754293,1.1967846,6.34];
ybs[5273]=['',3.690359,0.0382308,6.28];
ybs[5274]=['',3.6933935,-0.2869789,6.56];
ybs[5275]=['χ Cen',3.6976132,-0.7205821,4.36];
ybs[5276]=['',3.6982731,-0.753956,6.2];
ybs[5277]=['π Hya',3.6986044,-0.4675565,3.27];
ybs[5278]=['θ Cen',3.7002358,-0.6366335,2.06];
ybs[5279]=['',3.7084866,-1.1050373,6.4];
ybs[5280]=['95 Vir',3.6997058,-0.1644058,5.46];
ybs[5281]=['α Dra',3.687017,1.1216997,3.65];
ybs[5282]=['',3.7112012,-1.0364182,6.34];
ybs[5283]=['',3.7194279,-1.2288997,6.05];
ybs[5284]=['',3.7100164,-0.7605601,6.17];
ybs[5285]=['',3.7215889,-1.2186724,6.06];
ybs[5286]=['',3.7135215,-0.9007692,6];
ybs[5287]=['',3.7150725,-0.9345298,4.75];
ybs[5288]=['96 Vir',3.7097515,-0.182216,6.47];
ybs[5289]=['',3.7037103,0.7635523,5.27];
ybs[5290]=['13 Boo',3.7050292,0.8613554,5.25];
ybs[5291]=['',3.717877,-0.286359,4.91];
ybs[5292]=['',3.7065097,1.0337912,6.46];
ybs[5293]=['η Aps',3.7577231,-1.4156448,4.91];
ybs[5294]=['12 Boo',3.7150423,0.436093,4.83];
ybs[5295]=['3 UMi',3.6963437,1.3000452,6.45];
ybs[5296]=['',3.7497207,-1.3572916,6.47];
ybs[5297]=['',3.7204614,0.0219422,6.43];
ybs[5298]=['',3.7297832,-0.9384684,5.56];
ybs[5299]=['',3.7249091,-0.427063,6.34];
ybs[5300]=['',3.718568,0.5618289,6.11];
ybs[5301]=['',3.7315557,-0.9552214,6.11];
ybs[5302]=['50 Hya',3.726549,-0.4776223,5.08];
ybs[5303]=['',3.7236827,0.0402235,5.01];
ybs[5304]=['',3.7285118,-0.4662947,6.24];
ybs[5305]=['κ Vir',3.7267147,-0.181134,4.19];
ybs[5306]=['',3.7373442,-0.9981496,5.07];
ybs[5307]=['',3.7299324,-0.0165794,5.91];
ybs[5308]=['',3.7355,-0.7320177,5.61];
ybs[5309]=['',3.7388719,-0.9357327,6.39];
ybs[5310]=['',3.7457022,-1.1639802,5.75];
ybs[5311]=['4 UMi',3.7035779,1.3516093,4.82];
ybs[5312]=['',3.7329862,-0.1056216,6.36];
ybs[5313]=['14 Boo',3.7314102,0.2243652,5.54];
ybs[5314]=['',3.7364531,-0.5128805,6.08];
ybs[5315]=['',3.7397223,-0.7872232,6.31];
ybs[5316]=['',3.744657,-1.0474996,6.39];
ybs[5317]=['',3.7870863,-1.4477375,6.42];
ybs[5318]=['κ1 Boo',3.7274322,0.9020434,6.69];
ybs[5319]=['κ2 Boo',3.7275265,0.9020872,4.54];
ybs[5320]=['15 Boo',3.7347881,0.1744723,5.29];
ybs[5321]=['',3.735094,0.0564107,6.45];
ybs[5322]=['',3.7378289,-0.3194769,5.43];
ybs[5323]=['',3.7338043,0.379945,6.39];
ybs[5324]=['',3.7197035,1.2099936,5.24];
ybs[5325]=['',3.7319464,0.7228227,6.24];
ybs[5326]=['ε Aps',3.7752668,-1.3999337,5.06];
ybs[5327]=['',3.742175,-0.581979,6.55];
ybs[5328]=['ι Vir',3.7402319,-0.1065387,4.08];
ybs[5329]=['δ Oct',3.7996735,-1.4620193,4.32];
ybs[5330]=['α Boo',3.7381354,0.3329863,-0.04];
ybs[5331]=['',3.7417438,-0.1173821,6.44];
ybs[5332]=['',3.7422951,-0.0575942,6.15];
ybs[5333]=['',3.7399225,0.3282664,5.98];
ybs[5334]=['',3.7450928,-0.3261777,6.22];
ybs[5335]=['',3.7353131,0.9151091,6.58];
ybs[5336]=['',3.7419719,0.3493777,6.25];
ybs[5337]=['',3.7407759,0.6918684,6.38];
ybs[5338]=['',3.7513118,-0.5816041,6.54];
ybs[5339]=['',3.7591743,-1.0712036,5.23];
ybs[5340]=['ι Boo',3.7392238,0.8947177,4.75];
ybs[5341]=['λ Boo',3.7404347,0.8025851,4.18];
ybs[5342]=['',3.7461282,0.2645937,5.8];
ybs[5343]=['',3.7489724,-0.13344,6.47];
ybs[5344]=['ι Lup',3.7561907,-0.8056501,3.55];
ybs[5345]=['',3.7519724,-0.3284527,5.9];
ybs[5346]=['',3.7537911,-0.4523591,5.87];
ybs[5347]=['',3.7558117,-0.6476252,5.94];
ybs[5348]=['',3.760811,-0.9859177,4.33];
ybs[5349]=['λ Vir',3.753905,-0.2351623,4.52];
ybs[5350]=['',3.7443934,0.8936769,6.2];
ybs[5351]=['',3.7478639,0.6179575,4.81];
ybs[5352]=['',3.7593753,-0.7533055,5.56];
ybs[5353]=['',3.7466026,0.8359866,6.32];
ybs[5354]=['',3.761856,-0.7904488,4.77];
ybs[5355]=['18 Boo',3.7540235,0.2251736,5.41];
ybs[5356]=['υ Vir',3.7555405,-0.0413316,5.14];
ybs[5357]=['ψ Cen',3.7609036,-0.6630068,4.05];
ybs[5358]=['',3.7560933,0.0049155,6.19];
ybs[5359]=['',3.7518122,0.6748263,6.86];
ybs[5360]=['20 Boo',3.756047,0.2828206,4.86];
ybs[5361]=['',3.7709977,-1.0220817,4.92];
ybs[5362]=['',3.7510501,0.9557656,6.53];
ybs[5363]=['',3.7555874,0.6752916,6.33];
ybs[5364]=['',3.7573844,0.5293017,6.44];
ybs[5365]=['',3.7704712,-0.8451201,6.09];
ybs[5366]=['',3.768538,-0.608921,5.56];
ybs[5367]=['',3.7736337,-0.8879106,6.02];
ybs[5368]=['',3.7717995,-0.6913886,4.42];
ybs[5369]=['',3.7830243,-1.1919892,5.61];
ybs[5370]=['',3.7758248,-0.9298734,6];
ybs[5371]=['51 Hya',3.771678,-0.4861668,4.77];
ybs[5372]=['',3.7851405,-1.1566967,6.36];
ybs[5373]=['2 Lib',3.7727089,-0.2062193,6.21];
ybs[5374]=['',3.7716655,0.0199015,6.27];
ybs[5375]=['',3.7720336,0.145624,6.86];
ybs[5376]=['',3.7720408,0.145653,5.12];
ybs[5377]=['',3.7704684,0.4404617,6.22];
ybs[5378]=['',3.7748227,0.1421127,5.95];
ybs[5379]=['',3.8051277,-1.3409062,6.07];
ybs[5380]=['',3.7790777,-0.4347137,5.32];
ybs[5381]=['',3.7914872,-1.1505506,5.85];
ybs[5382]=['',3.7756342,0.0998137,5.1];
ybs[5383]=['',3.7781859,-0.2054366,6.49];
ybs[5384]=['',3.7760914,0.139346,6.19];
ybs[5385]=['τ1 Lup',3.7855854,-0.7910146,4.56];
ybs[5386]=['τ2 Lup',3.7857821,-0.793773,4.35];
ybs[5387]=['',3.7819475,-0.3502941,6.61];
ybs[5388]=['',3.7858368,-0.7403609,6.32];
ybs[5389]=['',3.7834366,-0.4704194,6.48];
ybs[5390]=['',3.7883965,-0.6976793,6.35];
ybs[5391]=['',3.7902909,-0.8069444,5.83];
ybs[5392]=['',3.7804015,0.6683275,6.27];
ybs[5393]=['',3.7977908,-1.0349331,6.45];
ybs[5394]=['θ Boo',3.7785158,0.903208,4.05];
ybs[5395]=['22 Boo',3.7852073,0.3338216,5.39];
ybs[5396]=['104 Vir',3.7899496,-0.108565,6.17];
ybs[5397]=['52 Hya',3.7939018,-0.5164679,4.97];
ybs[5398]=['',3.8099831,-1.1836099,5.83];
ybs[5399]=['φ Vir',3.7933287,-0.0406284,4.81];
ybs[5400]=['106 Vir',3.7955903,-0.122176,5.42];
ybs[5401]=['',3.7888762,0.7142747,6.63];
ybs[5402]=['',3.8031032,-0.7927417,5.5];
ybs[5403]=['',3.8042146,-0.8660005,5.37];
ybs[5404]=['',3.7939761,0.491999,7.62];
ybs[5405]=['',3.7941068,0.4920282,7.12];
ybs[5406]=['',3.7926236,0.6300141,6.1];
ybs[5407]=['',3.8063951,-0.7146045,6.39];
ybs[5408]=['',3.8004139,0.0127347,5.94];
ybs[5409]=['',3.8073576,-0.6801281,5.97];
ybs[5410]=['24 Boo',3.7935792,0.8682144,5.59];
ybs[5411]=['',3.8143254,-0.9945938,6.93];
ybs[5412]=['',3.7995491,0.5531269,6.06];
ybs[5413]=['',3.7982476,0.7277407,6.35];
ybs[5414]=['',3.8043137,0.0815642,6.02];
ybs[5415]=['σ Lup',3.8141787,-0.8823545,4.42];
ybs[5416]=['',3.8185198,-0.961611,5.87];
ybs[5417]=['',3.8181784,-0.9211487,5.87];
ybs[5418]=['',3.8157203,-0.5377753,6.09];
ybs[5419]=['ρ Boo',3.8083191,0.5283598,3.58];
ybs[5420]=['5 UMi',3.7851966,1.3193976,4.25];
ybs[5421]=['',3.8203955,-0.7364847,6.6];
ybs[5422]=['',3.8265814,-1.0491721,6.4];
ybs[5423]=['',3.8106405,0.4638875,6.01];
ybs[5424]=['26 Boo',3.8116601,0.3867937,5.92];
ybs[5425]=['γ Boo',3.8091267,0.6668868,3.03];
ybs[5426]=['',3.8018323,1.1010676,6.09];
ybs[5427]=['',3.8062505,1.0494115,6.27];
ybs[5428]=['',3.8227758,-0.358433,6.5];
ybs[5429]=['',3.8264519,-0.7263098,5.87];
ybs[5430]=['η Cen',3.8263996,-0.7374897,2.31];
ybs[5431]=['',3.814668,0.6433518,6.43];
ybs[5432]=['',3.8101398,0.9651559,5.76];
ybs[5433]=['',3.8384852,-1.1873234,6.04];
ybs[5434]=['',3.8301322,-0.8088251,5.55];
ybs[5435]=['',3.8185573,0.5661262,6.33];
ybs[5436]=['',3.8302095,-0.6927942,6.13];
ybs[5437]=['σ Boo',3.8207621,0.5174442,4.46];
ybs[5438]=['',3.8203599,0.637537,6.03];
ybs[5439]=['',3.8316914,-0.7035163,5.74];
ybs[5440]=['',3.8345777,-0.8068748,5.41];
ybs[5441]=['',3.8176285,0.9942696,6.48];
ybs[5442]=['',3.8198634,0.8599353,5.74];
ybs[5443]=['ρ Lup',3.8371673,-0.8643265,4.05];
ybs[5444]=['',3.8272128,0.4040986,6.38];
ybs[5445]=['',3.8319555,-0.2164569,6.2];
ybs[5446]=['',3.8385807,-0.6787666,6.02];
ybs[5447]=['',3.8426799,-0.8147223,6.07];
ybs[5448]=['',3.843813,-0.8578547,6.39];
ybs[5449]=['α1 Cen',3.8455327,-1.0634478,-0.01];
ybs[5450]=['α2 Cen',3.8455473,-1.0634526,1.33];
ybs[5451]=['',3.8492786,-0.9867448,6.3];
ybs[5452]=['',3.8365871,0.3176838,5.91];
ybs[5453]=['α Cir',3.8587639,-1.135686,3.19];
ybs[5454]=['',3.83561,0.7600125,5.7];
ybs[5455]=['',3.8555112,-1.0247017,6.22];
ybs[5456]=['',3.8502583,-0.6323389,5.67];
ybs[5457]=['',3.8352164,0.9412017,5.85];
ybs[5458]=['33 Boo',3.8383114,0.7733245,5.39];
ybs[5459]=['α Lup',3.8547288,-0.8287406,2.3];
ybs[5460]=['α Aps',3.8866903,-1.381207,3.83];
ybs[5461]=['',3.854415,-0.6612813,4];
ybs[5462]=['',3.8457804,0.3818764,6.1];
ybs[5463]=['',3.8475003,0.2345487,5.91];
ybs[5464]=['',3.8536903,-0.5415477,6.37];
ybs[5465]=['π1 Boo',3.8475112,0.284887,4.94];
ybs[5466]=['π2 Boo',3.8475331,0.2848773,5.88];
ybs[5467]=['ζ Boo',3.8494216,0.2379402,4.83];
ybs[5468]=['ζ Boo',3.8494216,0.2379402,4.43];
ybs[5469]=['',3.8096003,1.3886194,6.26];
ybs[5470]=['31 Boo',3.8517384,0.1407868,4.86];
ybs[5471]=['32 Boo',3.8519906,0.2018543,5.56];
ybs[5472]=['',3.8706672,-1.0990265,5.36];
ybs[5473]=['',3.8525126,0.3670167,6.38];
ybs[5474]=['4 Lib',3.8595048,-0.4379394,5.73];
ybs[5475]=['',3.8617257,-0.6155431,4.05];
ybs[5476]=['',3.8686091,-1.0222686,6.11];
ybs[5477]=['μ Vir',3.8582509,-0.1004087,3.88];
ybs[5478]=['',3.8694849,-0.9720744,6.1];
ybs[5479]=['',3.8675363,-0.6158547,4.92];
ybs[5480]=['34 Boo',3.858979,0.4613466,4.81];
ybs[5481]=['',4.1108151,-1.5386863,6.48];
ybs[5482]=['',3.8511268,1.0675622,6.25];
ybs[5483]=['',3.8598526,0.7044968,5.73];
ybs[5484]=['',3.8746422,-0.8296331,5.74];
ybs[5485]=['',3.8772939,-0.9158924,5.21];
ybs[5486]=['',3.8674605,-0.0263837,6.07];
ybs[5487]=['54 Hya',3.8716362,-0.4456983,4.94];
ybs[5488]=['',3.8780833,-0.9127835,6.07];
ybs[5489]=['',3.8720538,-0.4057297,5.81];
ybs[5490]=['',3.8863065,-1.1638964,5.91];
ybs[5491]=['108 Vir',3.868752,0.010881,5.69];
ybs[5492]=['ο Boo',3.867187,0.2944465,4.6];
ybs[5493]=['5 Lib',3.8711791,-0.2714568,6.33];
ybs[5494]=['',3.8722936,-0.3712251,6.4];
ybs[5495]=['ε Boo',3.8657774,0.4709072,5.12];
ybs[5496]=['ε Boo',3.8657774,0.4708927,2.7];
ybs[5497]=['',3.8675841,0.3279623,6.13];
ybs[5498]=['',3.876826,-0.6699271,5.94];
ybs[5499]=['',3.8790289,-0.7618448,6.3];
ybs[5500]=['',3.8666389,0.5706251,6.28];
ybs[5501]=['109 Vir',3.8719728,0.031403,3.72];
ybs[5502]=['',3.8709806,0.2624689,5.63];
ybs[5503]=['',3.8768677,-0.3738121,6.06];
ybs[5504]=['55 Hya',3.8776389,-0.4488554,5.63];
ybs[5505]=['',3.8867769,-0.9906515,6.23];
ybs[5506]=['56 Hya',3.8792765,-0.4569349,5.24];
ybs[5507]=['57 Hya',3.8802186,-0.4666879,5.77];
ybs[5508]=['',3.8796274,-0.2257167,6.35];
ybs[5509]=['',3.88353,-0.6410126,6.04];
ybs[5510]=['',3.9074799,-1.2789899,5.6];
ybs[5511]=['',3.8860527,-0.4248835,5.68];
ybs[5512]=['',3.8836199,-0.0164117,6.14];
ybs[5513]=['μ Lib',3.8857885,-0.248557,5.31];
ybs[5514]=['',3.8806931,0.4236596,6.14];
ybs[5515]=['π1 Oct',3.9534073,-1.4541179,5.65];
ybs[5516]=['58 Hya',3.8904365,-0.4896046,4.41];
ybs[5517]=['',3.9027202,-1.1152837,5.87];
ybs[5518]=['ο Lup',3.8969686,-0.7621334,4.32];
ybs[5519]=['',3.8833714,0.6583137,6.16];
ybs[5520]=['α1 Lib',3.8918128,-0.2808076,5.15];
ybs[5521]=['α2 Lib',3.8926507,-0.2815821,2.75];
ybs[5522]=['',3.8874669,0.4978316,5.8];
ybs[5523]=['38 Boo',3.8838548,0.8032643,5.74];
ybs[5524]=['',3.8888825,0.4157353,5.85];
ybs[5525]=['11 Lib',3.8928863,-0.0417295,4.94];
ybs[5526]=['',3.8927678,-0.0060911,6.18];
ybs[5527]=['',3.8845303,0.8950456,6.51];
ybs[5528]=['39 Boo',3.8853491,0.8487227,5.69];
ybs[5529]=['ζ Cir',3.9123379,-1.1533418,6.09];
ybs[5530]=['',3.9293002,-1.3395692,5.34];
ybs[5531]=['',3.8894147,0.6489123,5.48];
ybs[5532]=['',3.9004107,-0.5352593,6.29];
ybs[5533]=['',3.9019991,-0.6613813,5.03];
ybs[5534]=['ξ Boo',3.8939322,0.3317778,4.55];
ybs[5535]=['π2 Oct',3.9659873,-1.4507913,5.65];
ybs[5536]=['',3.915375,-1.0507597,5.2];
ybs[5537]=['',3.9397335,-1.3482363,5.93];
ybs[5538]=['12 Lib',3.9079993,-0.4316672,5.3];
ybs[5539]=['',3.9096023,-0.5827814,5.82];
ybs[5540]=['',3.902752,0.2725079,6.4];
ybs[5541]=['θ Cir',3.9207494,-1.0972939,5.11];
ybs[5542]=['',3.8921697,1.0332729,5.46];
ybs[5543]=['',3.9026817,0.3326926,6.01];
ybs[5544]=['ξ1 Lib',3.9078277,-0.2092444,5.8];
ybs[5545]=['',3.9375618,-1.3111022,6.2];
ybs[5546]=['',3.917842,-0.9232689,5.38];
ybs[5547]=['ω Oct',3.9987018,-1.4812712,5.91];
ybs[5548]=['',3.9144792,-0.5924656,5.32];
ybs[5549]=['',3.9185884,-0.8372131,5.64];
ybs[5550]=['',3.9209588,-0.8994791,6.64];
ybs[5551]=['',3.9184303,-0.6895048,6.36];
ybs[5552]=['',3.9177959,-0.571182,6.06];
ybs[5553]=['β UMi',3.8862602,1.2926511,2.08];
ybs[5554]=['ξ2 Lib',3.9182314,-0.2007007,5.46];
ybs[5555]=['',3.9207778,-0.5104592,6.29];
ybs[5556]=['',3.9257001,-0.8543743,6.35];
ybs[5557]=['',3.9151379,0.2505694,5.77];
ybs[5558]=['',3.9215715,-0.3753306,5.74];
ybs[5559]=['',3.9135144,0.5621765,6.12];
ybs[5560]=['16 Lib',3.9198501,-0.0774197,4.49];
ybs[5561]=['β Lup',3.9270617,-0.7543792,2.68];
ybs[5562]=['',3.9272632,-0.6980531,6.15];
ybs[5563]=['',3.9213599,-0.0044819,5.53];
ybs[5564]=['',3.9185916,0.3746484,6.49];
ybs[5565]=['',3.9193346,0.2844692,5.71];
ybs[5566]=['κ Cen',3.9297632,-0.7364029,3.13];
ybs[5567]=['59 Hya',3.9269559,-0.4842601,5.65];
ybs[5568]=['17 Lib',3.9245749,-0.1962453,6.6];
ybs[5569]=['',3.9298743,-0.6627013,6.47];
ybs[5570]=['',3.9310832,-0.7548287,6.1];
ybs[5571]=['',3.9144523,0.864615,5.63];
ybs[5572]=['18 Lib',3.9274989,-0.1960517,5.87];
ybs[5573]=['',3.9272715,-0.0886269,6.09];
ybs[5574]=['',3.9292169,0.0781765,5.93];
ybs[5575]=['',3.9385537,-0.665776,5.89];
ybs[5576]=['δ Lib',3.9365039,-0.1502179,4.92];
ybs[5577]=['',3.9416819,-0.6012035,6.22];
ybs[5578]=['40 Boo',3.9290806,0.6837628,5.64];
ybs[5579]=['',3.9180577,1.1491783,4.6];
ybs[5580]=['',3.9379019,-0.0496167,5.52];
ybs[5581]=['60 Hya',3.9420502,-0.4912762,5.85];
ybs[5582]=['',3.9352048,0.3832309,6.38];
ybs[5583]=['η Cir',3.95624,-1.11907,5.17];
ybs[5584]=['',3.939947,-0.0039829,5.71];
ybs[5585]=['',3.9460694,-0.5712545,5.44];
ybs[5586]=['',3.8787246,1.4384909,5.64];
ybs[5587]=['',3.933158,0.8236139,6.37];
ybs[5588]=['',3.96811,-1.2564728,6.52];
ybs[5589]=['',3.9440987,-0.054431,6.61];
ybs[5590]=['ω Boo',3.9404941,0.4349446,4.81];
ybs[5591]=['110 Vir',3.9446201,0.0349793,4.4];
ybs[5592]=['β Boo',3.9391898,0.7034184,3.5];
ybs[5593]=['σ Lib',3.9505278,-0.4427669,3.29];
ybs[5594]=['',3.9539664,-0.7146745,6.41];
ybs[5595]=['π Lup',3.9560682,-0.8227025,4.72];
ybs[5596]=['π Lup',3.9560682,-0.8227025,4.82];
ybs[5597]=['',3.9566124,-0.7182628,5.15];
ybs[5598]=['',3.9356415,1.0492315,5.93];
ybs[5599]=['',3.9444504,0.6129359,5.51];
ybs[5600]=['',3.9497908,0.0943481,6.5];
ybs[5601]=['',3.9701887,-1.140758,6.17];
ybs[5602]=['',3.9440448,0.7776704,6.65];
ybs[5603]=['',3.9466848,0.6017741,6.59];
ybs[5604]=['',3.9580772,-0.4516173,6.67];
ybs[5605]=['',3.9603831,-0.6344275,6.27];
ybs[5606]=['ψ Boo',3.9506164,0.4688103,4.54];
ybs[5607]=['',3.9662989,-0.8582425,5.77];
ybs[5608]=['44 Boo',3.9468354,0.8302095,4.76];
ybs[5609]=['',3.9615804,-0.5411327,5.96];
ybs[5610]=['',3.9608193,-0.3860273,6.17];
ybs[5611]=['',3.9772954,-1.172313,5.76];
ybs[5612]=['ν Lib',3.961403,-0.2852334,5.2];
ybs[5613]=['',3.976411,-1.1122504,6.28];
ybs[5614]=['',3.969186,-0.709807,5.79];
ybs[5615]=['',3.9712673,-0.7496651,5.85];
ybs[5616]=['λ Lup',3.9722364,-0.79176,4.05];
ybs[5617]=['47 Boo',3.9539489,0.8388891,5.57];
ybs[5618]=['',3.9921287,-1.2715314,6.01];
ybs[5619]=['',3.945802,1.1489984,6.13];
ybs[5620]=['',3.9595803,0.6347717,6.35];
ybs[5621]=['',3.9653437,0.0944698,6.16];
ybs[5622]=['',3.9818832,-1.0734902,6.3];
ybs[5623]=['',3.9635178,0.3203758,6.02];
ybs[5624]=['45 Boo',3.963143,0.4325563,4.93];
ybs[5625]=['',3.957187,0.9506878,5.25];
ybs[5626]=['',3.9775032,-0.6785279,5.98];
ybs[5627]=['',3.9835378,-0.967434,5.54];
ybs[5628]=['46 Boo',3.9678683,0.4575561,5.67];
ybs[5629]=['',3.9704356,0.2295132,6.1];
ybs[5630]=['',3.9687596,0.4367444,5.81];
ybs[5631]=['',3.9778108,-0.4610638,5.76];
ybs[5632]=['',3.9842023,-0.7917019,6.44];
ybs[5633]=['',3.9839839,-0.7917362,7.39];
ybs[5634]=['',3.9990511,-1.2245555,5.81];
ybs[5635]=['',3.9918697,-1.0790833,6.32];
ybs[5636]=['κ1 Lup',3.9859587,-0.8520924,3.87];
ybs[5637]=['κ2 Lup',3.9860683,-0.8521941,5.69];
ybs[5638]=['',3.9664356,0.8721383,6.39];
ybs[5639]=['ζ Lup',3.9877251,-0.9107571,3.41];
ybs[5640]=['',3.9884889,-0.8430271,6.33];
ybs[5641]=['',3.989591,-0.7781328,4.82];
ybs[5642]=['ι1 Lib',3.9859376,-0.3468864,4.54];
ybs[5643]=['',3.9904786,-0.6313634,6.1];
ybs[5644]=['',3.9841365,0.3297321,5.89];
ybs[5645]=['',3.990749,-0.4204735,6.47];
ybs[5646]=['ι2 Lib',3.9907279,-0.3443625,6.08];
ybs[5647]=['23 Lib',3.991609,-0.4431759,6.45];
ybs[5648]=['',3.9934305,-0.4586095,5.84];
ybs[5649]=['',3.9869759,0.3351472,6.68];
ybs[5650]=['1 Lup',3.9968325,-0.5515524,4.91];
ybs[5651]=['',4.0074698,-1.0643967,5.73];
ybs[5652]=['26 Lib',3.9960883,-0.3115609,6.17];
ybs[5653]=['',4.0032295,-0.8404819,5.95];
ybs[5654]=['δ Cir',4.0089564,-1.0653299,5.09];
ybs[5655]=['',3.9903725,0.3996863,6.3];
ybs[5656]=['ε Cir',4.012392,-1.111629,4.86];
ybs[5657]=['',4.0035926,-0.725585,5.16];
ybs[5658]=['',4.004172,-0.7603792,6.04];
ybs[5659]=['',3.9969558,-0.0974799,6.28];
ybs[5660]=['β Cir',4.0111925,-1.0276899,4.07];
ybs[5661]=['γ TrA',4.0188326,-1.2000879,2.89];
ybs[5662]=['',3.9749115,1.1815374,6.17];
ybs[5663]=['',3.9900449,0.6663969,6.2];
ybs[5664]=['',3.9925373,0.5533616,5.99];
ybs[5665]=['3 Ser',3.9981602,0.0847735,5.33];
ybs[5666]=['χ Boo',3.9943169,0.5075689,5.26];
ybs[5667]=['',3.9923738,0.7345851,6.13];
ybs[5668]=['',4.0041924,-0.3923709,5.5];
ybs[5669]=['4 Ser',4.0010335,0.005065,5.63];
ybs[5670]=['',4.0170801,-1.0572685,5.46];
ybs[5671]=['δ Boo',3.9985861,0.580017,3.47];
ybs[5672]=['',4.0127011,-0.71806,6.28];
ybs[5673]=['μ Lup',4.0147393,-0.8369866,4.27];
ybs[5674]=['',4.0263217,-1.1791651,6.28];
ybs[5675]=['β Lib',4.0065037,-0.1651879,2.61];
ybs[5676]=['2 Lup',4.0107969,-0.5276135,4.34];
ybs[5677]=['',4.0161093,-0.7132983,5.59];
ybs[5678]=['',4.0145798,-0.5461174,6.18];
ybs[5679]=['',4.0184932,-0.6488625,6.2];
ybs[5680]=['',4.0124823,-0.0094653,5.89];
ybs[5681]=['',3.9919674,1.1739779,5.13];
ybs[5682]=['',4.0117193,0.3576495,5.7];
ybs[5683]=['',3.9926234,1.2018803,6.51];
ybs[5684]=['5 Ser',4.0162424,0.0294038,5.06];
ybs[5685]=['δ Lup',4.026725,-0.7108225,3.22];
ybs[5686]=['',4.0276768,-0.712605,6.2];
ybs[5687]=['',4.0271727,-0.6684392,6.48];
ybs[5688]=['ν1 Lup',4.0304969,-0.8378813,5];
ybs[5689]=['ν2 Lup',4.0290588,-0.8446905,5.65];
ybs[5690]=['',4.0361636,-1.0600385,5.67];
ybs[5691]=['28 Lib',4.0237479,-0.3183217,6.17];
ybs[5692]=['',4.0160562,0.5660931,6.32];
ybs[5693]=['ο Lib',4.0242242,-0.2727628,6.3];
ybs[5694]=['γ Cir',4.036901,-1.0367177,4.51];
ybs[5695]=['φ1 Lup',4.0284137,-0.6342674,3.56];
ybs[5696]=['',4.0227863,-0.0435159,6.35];
ybs[5697]=['',4.0243781,-0.1030533,5.54];
ybs[5698]=['ε Lup',4.0326698,-0.7813578,3.37];
ybs[5699]=['ο CrB',4.0189545,0.5154979,5.51];
ybs[5700]=['6 Ser',4.0237779,0.0110905,5.35];
ybs[5701]=['',4.0233716,0.4342018,6.39];
ybs[5702]=['φ2 Lup',4.0343376,-0.6446809,4.54];
ybs[5703]=['',4.0508676,-1.193571,5.89];
ybs[5704]=['11 UMi',4.0015938,1.2521369,5.02];
ybs[5705]=['',4.0175265,0.9054466,5.66];
ybs[5706]=['',4.020668,0.7741251,6.19];
ybs[5707]=['7 Ser',4.0293267,0.2179605,6.28];
ybs[5708]=['50 Boo',4.0260945,0.5734161,5.37];
ybs[5709]=['υ Lup',4.041428,-0.6944399,5.37];
ybs[5710]=['',4.0365495,-0.2172597,5.72];
ybs[5711]=['8 Ser',4.0355885,-0.0192196,6.12];
ybs[5712]=['',4.0429317,-0.6675447,7.03];
ybs[5713]=['ε Lib',4.0379193,-0.1815266,4.94];
ybs[5714]=['',4.0439493,-0.6773894,4.6];
ybs[5715]=['',4.0558693,-1.1276271,5.71];
ybs[5716]=['',4.0293375,0.6894426,5.5];
ybs[5717]=['η CrB',4.0322838,0.5272431,5.58];
ybs[5718]=['η CrB',4.0322838,0.5272431,6.08];
ybs[5719]=['ρ Oct',4.1398018,-1.4754037,5.57];
ybs[5720]=['κ1 Aps',4.0754071,-1.2821973,5.49];
ybs[5721]=['',4.027545,1.0815438,5.98];
ybs[5722]=['',4.0353976,0.7887575,6.01];
ybs[5723]=['μ1 Boo',4.03758,0.6509866,4.31];
ybs[5724]=['μ2 Boo',4.0376904,0.6504728,6.5];
ybs[5725]=['γ UMi',4.0173615,1.2523379,3.05];
ybs[5726]=['',4.0524451,-0.6430646,5.45];
ybs[5727]=['',4.0274438,1.1041313,5.79];
ybs[5728]=['',4.0583595,-0.9018826,6.1];
ybs[5729]=['τ1 Ser',4.0440812,0.2679118,5.17];
ybs[5730]=['',4.044377,0.3386467,6.27];
ybs[5731]=['',4.0455643,0.5979178,5.46];
ybs[5732]=['',4.0621724,-0.8169702,5.24];
ybs[5733]=['ζ1 Lib',4.0558358,-0.2930953,5.64];
ybs[5734]=['ι Dra',4.037969,1.0277857,3.29];
ybs[5735]=['',4.0518551,0.4367617,6.02];
ybs[5736]=['10 Ser',4.0569191,0.0308159,5.17];
ybs[5737]=['β CrB',4.052488,0.5066488,3.68];
ybs[5738]=['',4.0454722,0.9414723,6.45];
ybs[5739]=['',4.0662225,-0.3630991,6.22];
ybs[5740]=['ζ3 Lib',4.0663811,-0.2912105,5.82];
ybs[5741]=['',4.0733694,-0.675409,6.25];
ybs[5742]=['',4.0555667,0.8224819,6.15];
ybs[5743]=['',4.0720759,-0.5751953,6.46];
ybs[5744]=['',4.049528,1.0855666,6.5];
ybs[5745]=['',4.0505079,1.0575507,5.9];
ybs[5746]=['',4.071085,-0.3532537,6.22];
ybs[5747]=['',4.1119367,-1.361174,6.18];
ybs[5748]=['',4.0666893,0.1484149,6.57];
ybs[5749]=['',4.055841,0.9619978,6.43];
ybs[5750]=['',4.0635177,0.544721,6.46];
ybs[5751]=['',4.0636514,0.6410298,6.37];
ybs[5752]=['',4.0749459,-0.3446221,5.52];
ybs[5753]=['ν1 Boo',4.0654766,0.7113505,5.02];
ybs[5754]=['ζ4 Lib',4.076196,-0.2954403,5.5];
ybs[5755]=['',4.0774755,-0.4287232,7];
ybs[5756]=['',4.0942919,-1.1464436,6.51];
ybs[5757]=['',4.0819757,-0.7005801,5.82];
ybs[5758]=['',4.056835,1.0825103,6.38];
ybs[5759]=['',4.0676197,0.6377593,6.38];
ybs[5760]=['τ2 Ser',4.0718519,0.2789215,6.22];
ybs[5761]=['ε TrA',4.0963057,-1.1587204,4.11];
ybs[5762]=['11 Ser',4.075897,-0.0220103,5.51];
ybs[5763]=['',4.0833434,-0.6880695,6.36];
ybs[5764]=['ν2 Boo',4.0691922,0.7125157,5.02];
ybs[5765]=['36 Lib',4.0840297,-0.4908019,5.15];
ybs[5766]=['γ Lup',4.0869051,-0.7197846,2.78];
ybs[5767]=['37 Lib',4.0814675,-0.1769521,4.62];
ybs[5768]=['θ CrB',4.0746394,0.5460155,4.14];
ybs[5769]=['',4.0820651,-0.1006897,6.51];
ybs[5770]=['',4.0825962,-0.1615717,5.17];
ybs[5771]=['',4.0903798,-0.7859557,4.54];
ybs[5772]=['κ2 Aps',4.114166,-1.2831275,5.65];
ybs[5773]=['',4.0793106,0.2978132,6.45];
ybs[5774]=['',4.0917216,-0.7761504,5.43];
ybs[5775]=['',4.0633765,1.1193291,5.79];
ybs[5776]=['',4.1222082,-1.3291083,5.95];
ybs[5777]=['γ Lib',4.0875052,-0.2594083,3.91];
ybs[5778]=['δ Ser',4.0835463,0.182624,3.8];
ybs[5779]=['δ Ser',4.0835462,0.1826531,3.8];
ybs[5780]=['',4.0910888,-0.578856,6.24];
ybs[5781]=['',4.085023,0.0278398,6.56];
ybs[5782]=['',4.1124038,-1.2269494,6.44];
ybs[5783]=['α CrB',4.0824926,0.4649684,2.23];
ybs[5784]=['υ Lib',4.094527,-0.4923202,3.58];
ybs[5785]=['τ3 Ser',4.0865914,0.3068631,6.12];
ybs[5786]=['',4.0882706,0.1953392,6.07];
ybs[5787]=['ω Lup',4.0997026,-0.7442062,4.33];
ybs[5788]=['',4.1037524,-0.915334,5.44];
ybs[5789]=['14 Ser',4.0915715,-0.0110791,6.51];
ybs[5790]=['μ CrB',4.0844004,0.6795652,5.11];
ybs[5791]=['',4.0964188,-0.4599409,6.19];
ybs[5792]=['16 Ser',4.0909427,0.1734304,5.26];
ybs[5793]=['',4.1094041,-1.0468448,5.95];
ybs[5794]=['τ5 Ser',4.0907214,0.2800504,5.93];
ybs[5795]=['',4.1016641,-0.684745,6.57];
ybs[5796]=['',4.0977266,-0.4051642,5.78];
ybs[5797]=['',4.1023612,-0.6841717,6.04];
ybs[5798]=['',4.0869243,0.6684675,6.42];
ybs[5799]=['',4.0999294,-0.4935614,6.32];
ybs[5800]=['',4.0997013,-0.3680627,5.84];
ybs[5801]=['',4.0835033,0.9398271,5.97];
ybs[5802]=['τ Lib',4.1017182,-0.5209794,3.66];
ybs[5803]=['',4.0919619,0.522169,6.52];
ybs[5804]=['41 Lib',4.1024198,-0.3381399,5.38];
ybs[5805]=['',4.1010153,-0.1547516,6.5];
ybs[5806]=['',4.1010225,-0.1546934,6.48];
ybs[5807]=['',4.087131,0.907506,6.74];
ybs[5808]=['',4.0864012,0.9521999,5.74];
ybs[5809]=['',4.1045082,-0.4053024,6.34];
ybs[5810]=['ψ1 Lup',4.1067676,-0.6018519,4.67];
ybs[5811]=['',4.1127788,-0.8343827,6.23];
ybs[5812]=['',4.1087653,-0.5460268,6.34];
ybs[5813]=['φ Boo',4.0955611,0.7030307,5.24];
ybs[5814]=['42 Lib',4.1085786,-0.4169501,4.96];
ybs[5815]=['',4.1135231,-0.780722,4.64];
ybs[5816]=['θ UMi',4.0613642,1.34868,4.96];
ybs[5817]=['',4.1001505,0.6039331,6.11];
ybs[5818]=['',4.0933012,0.9500887,5.87];
ybs[5819]=['',4.0491416,1.4027515,6.58];
ybs[5820]=['',4.0971104,0.8155106,5.75];
ybs[5821]=['',4.1069298,0.209117,6.25];
ybs[5822]=['',4.1201078,-0.8649808,6.04];
ybs[5823]=['ζ1 CrB',4.1024773,0.6381749,6];
ybs[5824]=['ζ2 CrB',4.1025137,0.6381604,5.07];
ybs[5825]=['',4.0981667,0.8787906,5.84];
ybs[5826]=['',4.1268409,-1.0534262,6.48];
ybs[5827]=['',4.1194535,-0.654417,5.24];
ybs[5828]=['κ Lib',4.1156962,-0.344695,4.74];
ybs[5829]=['ψ2 Lup',4.1195187,-0.6070408,4.75];
ybs[5830]=['τ6 Ser',4.1103451,0.278442,6.01];
ybs[5831]=['',4.1000346,1.0097137,6.45];
ybs[5832]=['ι Ser',4.112683,0.3420732,4.52];
ybs[5833]=['χ Ser',4.1139557,0.2229954,5.33];
ybs[5834]=['',4.0916163,1.2079503,5.62];
ybs[5835]=['τ7 Ser',4.1142956,0.3210207,5.81];
ybs[5836]=['',4.1272814,-0.7311005,5.94];
ybs[5837]=['',4.1219446,-0.2637779,6.31];
ybs[5838]=['η Lib',4.1248395,-0.2747585,5.41];
ybs[5839]=['γ CrB',4.1176334,0.4577153,3.84];
ybs[5840]=['',4.1199771,0.2373228,6.48];
ybs[5841]=['',4.1450262,-1.1433704,6.18];
ybs[5842]=['',4.1450335,-1.1433703,6.39];
ybs[5843]=['ψ Ser',4.124061,0.0426772,5.88];
ybs[5844]=['α Ser',4.1249761,0.1109311,2.65];
ybs[5845]=['π CrB',4.1228025,0.5662893,5.56];
ybs[5846]=['',4.134642,-0.490968,6.51];
ybs[5847]=['',4.1166277,0.9126397,5.51];
ybs[5848]=['τ8 Ser',4.1265063,0.3001037,6.14];
ybs[5849]=['',4.1299083,0.09386,5.58];
ybs[5850]=['',4.1372182,-0.6065188,5.61];
ybs[5851]=['',4.1312318,0.0143529,6.33];
ybs[5852]=['',4.1404976,-0.7027098,6.42];
ybs[5853]=['25 Ser',4.1332013,-0.0326948,5.4];
ybs[5854]=['',4.1406385,-0.6629547,6.01];
ybs[5855]=['',4.1475184,-0.9163939,6.07];
ybs[5856]=['',4.1362344,-0.1080149,6.24];
ybs[5857]=['β Ser',4.1330544,0.2679627,3.67];
ybs[5858]=['λ Ser',4.134435,0.1271362,4.43];
ybs[5859]=['',4.1531658,-0.9298468,5.77];
ybs[5860]=['υ Ser',4.1378972,0.2451659,5.71];
ybs[5861]=['',4.1521286,-0.8548429,5.84];
ybs[5862]=['',4.1532516,-0.7935747,6.12];
ybs[5863]=['',4.1577031,-0.9620642,5.73];
ybs[5864]=['',4.1419795,0.239472,6];
ybs[5865]=['',4.1457166,-0.0678257,5.53];
ybs[5866]=['',4.1670701,-1.1382681,6.54];
ybs[5867]=['',4.1404735,0.5527079,6.44];
ybs[5868]=['',4.1326129,0.9670167,5.92];
ybs[5869]=['κ Ser',4.1440892,0.3154511,4.09];
ybs[5870]=['',4.1429809,0.4902441,5.85];
ybs[5871]=['μ Ser',4.1486423,-0.0610426,3.53];
ybs[5872]=['',4.1588393,-0.8225175,6.01];
ybs[5873]=['χ Lup',4.1556038,-0.5880669,3.95];
ybs[5874]=['',4.1685426,-1.0938366,6.19];
ybs[5875]=['1 Sco',4.1553524,-0.4506078,4.64];
ybs[5876]=['',4.132124,1.0913664,5.19];
ybs[5877]=['',4.1372107,0.9653137,5.86];
ybs[5878]=['ω Ser',4.1513915,0.0371665,5.23];
ybs[5879]=['δ CrB',4.1475207,0.4538041,4.63];
ybs[5880]=['',4.1649319,-0.884548,6.6];
ybs[5881]=['κ TrA',4.1789562,-1.1984701,5.09];
ybs[5882]=['ε Ser',4.1536083,0.0769884,3.71];
ybs[5883]=['',4.1609159,-0.5227718,6.4];
ybs[5884]=['',4.152722,0.2629664,5.2];
ybs[5885]=['36 Ser',4.1557877,-0.0551002,5.11];
ybs[5886]=['',4.1578129,-0.2478344,6.19];
ybs[5887]=['β TrA',4.1763774,-1.1081968,2.85];
ybs[5888]=['',4.1748155,-1.0612985,6.15];
ybs[5889]=['ρ Ser',4.1549992,0.3649708,4.76];
ybs[5890]=['',4.1776319,-1.0514224,5.77];
ybs[5891]=['κ CrB',4.1542355,0.6211795,4.82];
ybs[5892]=['λ Lib',4.1654263,-0.3531268,5.03];
ybs[5893]=['ζ UMi',4.1159168,1.3565441,4.32];
ybs[5894]=['2 Sco',4.1668318,-0.4431834,4.59];
ybs[5895]=['',4.1801283,-1.0567409,5.76];
ybs[5896]=['',4.1680445,-0.4293203,5.39];
ybs[5897]=['',4.1681682,-0.4196334,5.42];
ybs[5898]=['θ Lib',4.1674406,-0.2931224,4.15];
ybs[5899]=['',4.1624201,0.3026035,6.36];
ybs[5900]=['',4.1707917,-0.4782814,6.14];
ybs[5901]=['39 Ser',4.163723,0.2291808,6.1];
ybs[5902]=['3 Sco',4.1713979,-0.4417155,5.87];
ybs[5903]=['',4.1652805,0.2794202,6.09];
ybs[5904]=['χ Her',4.1601617,0.7397714,4.62];
ybs[5905]=['47 Lib',4.172687,-0.3394269,5.94];
ybs[5906]=['',4.1753478,-0.5436357,6.21];
ybs[5907]=['4 Sco',4.1751208,-0.4595498,5.62];
ybs[5908]=['',4.1784116,-0.6968799,6.03];
ybs[5909]=['40 Ser',4.1702916,0.1486217,6.29];
ybs[5910]=['',4.1933937,-1.1362164,5.75];
ybs[5911]=['',4.1831141,-0.8417,6.31];
ybs[5912]=['',4.157391,0.9732049,5.81];
ybs[5913]=['',4.178545,-0.5558857,6.29];
ybs[5914]=['',4.1694674,0.3533573,5.44];
ybs[5915]=['ξ1 Lup',4.1815285,-0.5939381,5.12];
ybs[5916]=['ξ2 Lup',4.1815793,-0.5938992,5.62];
ybs[5917]=['',4.1778992,-0.2524366,6.37];
ybs[5918]=['ρ Sco',4.1812839,-0.5109964,3.88];
ybs[5919]=['',4.1836608,-0.6326563,5.8];
ybs[5920]=['',4.1792963,-0.2599341,6.13];
ybs[5921]=['',4.1742691,0.3238654,6.26];
ybs[5922]=['2 Her',4.1686487,0.7517767,5.37];
ybs[5923]=['γ Ser',4.1778211,0.2722295,3.85];
ybs[5924]=['',4.1843942,-0.3673303,5.85];
ybs[5925]=['',4.1887797,-0.6556511,6.31];
ybs[5926]=['λ CrB',4.1740065,0.6611748,5.45];
ybs[5927]=['',4.1960017,-0.943933,6.1];
ybs[5928]=['4 Her',4.1725119,0.7417919,5.75];
ybs[5929]=['',4.2028336,-1.1141877,6.41];
ybs[5930]=['φ Ser',4.181311,0.250468,5.54];
ybs[5931]=['48 Lib',4.1864059,-0.2503261,4.88];
ybs[5932]=['',4.1885018,-0.4344888,5.43];
ybs[5933]=['',4.1933547,-0.7296689,4.99];
ybs[5934]=['π Sco',4.1897398,-0.4568752,2.89];
ybs[5935]=['',4.1953071,-0.7106119,6.49];
ybs[5936]=['',4.201329,-0.9536385,6.13];
ybs[5937]=['ε CrB',4.1823454,0.4679965,4.15];
ybs[5938]=['η Lup',4.1958583,-0.671239,3.41];
ybs[5939]=['',4.1725152,1.027077,6.31];
ybs[5940]=['',4.1813449,0.6917028,6.31];
ybs[5941]=['',4.2100426,-1.092619,6.25];
ybs[5942]=['',4.199349,-0.7068081,6.21];
ybs[5943]=['δ Sco',4.1960685,-0.3959075,2.32];
ybs[5944]=['49 Lib',4.1958119,-0.2896462,5.47];
ybs[5945]=['',4.2257816,-1.264666,5.7];
ybs[5946]=['',4.2007888,-0.5576521,6.33];
ybs[5947]=['',4.1878868,0.6384582,5.62];
ybs[5948]=['',4.190736,0.4513009,6.33];
ybs[5949]=['50 Lib',4.1975655,-0.1478879,5.55];
ybs[5950]=['',4.1815163,0.9544538,4.95];
ybs[5951]=['ι1 Nor',4.2122722,-1.0094249,4.63];
ybs[5952]=['η Nor',4.2100702,-0.8602802,4.65];
ybs[5953]=['',4.1973919,0.0761931,5.83];
ybs[5954]=['',4.1875376,0.8694915,6.05];
ybs[5955]=['',4.2064794,-0.5095815,6.03];
ybs[5956]=['5 Her',4.1986151,0.3099099,5.12];
ybs[5957]=['',4.2102005,-0.6747993,4.89];
ybs[5958]=['ρ CrB',4.1971346,0.5801768,5.41];
ybs[5959]=['',4.2093383,-0.452494,5];
ybs[5960]=['',4.2106033,-0.5595726,6.01];
ybs[5961]=['ι CrB',4.1990293,0.5199228,4.99];
ybs[5962]=['π Ser',4.2030354,0.3969425,4.83];
ybs[5963]=['',4.2117729,-0.4326119,6.21];
ybs[5964]=['',4.2138261,-0.5807528,6.1];
ybs[5965]=['',4.2154372,-0.6618836,5.9];
ybs[5966]=['43 Ser',4.2100694,0.0859765,6.08];
ybs[5967]=['ξ Sco',4.213269,-0.199549,5.07];
ybs[5968]=['ξ Sco',4.213269,-0.199549,4.77];
ybs[5969]=['',4.2289984,-0.9817485,6.16];
ybs[5970]=['δ Nor',4.2240808,-0.7894558,4.72];
ybs[5971]=['',4.2004316,0.9224821,5.93];
ybs[5972]=['υ Her',4.2040555,0.8024244,4.76];
ybs[5973]=['',4.2069035,0.6382811,5.83];
ybs[5974]=['β1 Sco',4.218237,-0.3467146,2.62];
ybs[5975]=['β2 Sco',4.2182588,-0.3466515,4.92];
ybs[5976]=['θ Dra',4.1989512,1.0210812,4.01];
ybs[5977]=['θ Lup',4.224027,-0.6433519,4.23];
ybs[5978]=['',4.2213024,-0.4130459,5.92];
ybs[5979]=['',4.2190877,-0.1108506,6.53];
ybs[5980]=['',4.2201953,-0.1081917,6.41];
ybs[5981]=['',4.2270012,-0.6425318,5.73];
ybs[5982]=['',4.218115,0.1402621,6.29];
ybs[5983]=['ω1 Sco',4.2242501,-0.3617759,3.96];
ybs[5984]=['ι2 Nor',4.2375451,-1.0121543,5.57];
ybs[5985]=['',4.2043625,1.0358495,6.19];
ybs[5986]=['',4.225097,-0.2466114,6.32];
ybs[5987]=['ω2 Sco',4.2268695,-0.3652518,4.32];
ybs[5988]=['',4.2290199,-0.4279633,6.33];
ybs[5989]=['',4.2327911,-0.6835311,7.05];
ybs[5990]=['',4.2328049,-0.6833177,6.65];
ybs[5991]=['',4.2302394,-0.4605065,5.38];
ybs[5992]=['11 Sco',4.2274501,-0.2234763,5.78];
ybs[5993]=['',4.2327571,-0.4144055,5.88];
ybs[5994]=['45 Ser',4.2267565,0.1716172,5.63];
ybs[5995]=['',4.2252081,0.3798469,6.14];
ybs[5996]=['',4.2366418,-0.5708526,6.19];
ybs[5997]=['',4.238206,-0.5864896,5.54];
ybs[5998]=['κ Her',4.2284641,0.2965038,5];
ybs[5999]=['κ Her',4.2284929,0.2966348,6.25];
ybs[6000]=['47 Ser',4.2304754,0.1479314,5.73];
ybs[6001]=['',4.2328959,0.0592781,5.91];
ybs[6002]=['',4.2377509,-0.3211126,6.47];
ybs[6003]=['8 Her',4.2315189,0.2992829,6.14];
ybs[6004]=['',4.2336907,0.1103209,5.97];
ybs[6005]=['',4.2448081,-0.7186665,5.86];
ybs[6006]=['',4.2368829,-0.0615154,5.37];
ybs[6007]=['',4.2430629,-0.5144076,5.13];
ybs[6008]=['τ CrB',4.2315373,0.6358705,4.76];
ybs[6009]=['ζ Nor',4.2550445,-0.9703435,5.81];
ybs[6010]=['δ1 Aps',4.2925782,-1.3744075,4.68];
ybs[6011]=['δ2 Aps',4.2929904,-1.3739072,5.27];
ybs[6012]=['',4.2544313,-0.9377214,5.83];
ybs[6013]=['φ Her',4.2301592,0.7832468,4.26];
ybs[6014]=['κ Nor',4.2553936,-0.9544554,4.94];
ybs[6015]=['',4.216727,1.1824725,5.44];
ybs[6016]=['ν Sco',4.24677,-0.3404489,6.3];
ybs[6017]=['ν Sco',4.2468504,-0.3406378,4.01];
ybs[6018]=['13 Sco',4.2485414,-0.4883914,4.59];
ybs[6019]=['12 Sco',4.2483957,-0.4969632,5.67];
ybs[6020]=['δ TrA',4.2651664,-1.112477,3.85];
ybs[6021]=['ψ Sco',4.246529,-0.1766401,4.94];
ybs[6022]=['',4.2436459,0.1685228,6.53];
ybs[6023]=['16 Sco',4.247007,-0.1501683,5.43];
ybs[6024]=['',4.2010756,1.3392345,5.56];
ybs[6025]=['',4.2433236,0.2898762,6.08];
ybs[6026]=['',4.2302154,1.0101893,6.33];
ybs[6027]=['',4.273277,-1.1867402,5.75];
ybs[6028]=['',4.2321312,0.9733857,6.49];
ybs[6029]=['10 Her',4.2437361,0.4090687,5.7];
ybs[6030]=['',4.2660284,-1.0117108,5.63];
ybs[6031]=['',4.2504358,-0.0746469,6.25];
ybs[6032]=['',4.254759,-0.427215,6.41];
ybs[6033]=['',4.243425,0.5809446,6.29];
ybs[6034]=['',4.2578056,-0.5771235,5.92];
ybs[6035]=['θ Nor',4.2625217,-0.827759,5.14];
ybs[6036]=['',4.2438725,0.6347455,5.63];
ybs[6037]=['9 Her',4.2514942,0.0866578,5.48];
ybs[6038]=['χ Sco',4.2546582,-0.2075749,5.22];
ybs[6039]=['',4.2628456,-0.7496982,6.14];
ybs[6040]=['',4.243504,0.7385825,5.87];
ybs[6041]=['',4.2577611,-0.3693611,6.41];
ybs[6042]=['',4.2484962,0.4645117,6.5];
ybs[6043]=['',4.2584108,-0.324466,6.32];
ybs[6044]=['',4.2597429,-0.4456186,6.05];
ybs[6045]=['',4.2694671,-0.9401257,5.44];
ybs[6046]=['δ Oph',4.2565414,-0.0654477,2.74];
ybs[6047]=['',4.2556886,0.1020395,6.31];
ybs[6048]=['γ1 Nor',4.2704228,-0.8747997,4.99];
ybs[6049]=['',4.272149,-0.9274764,6.33];
ybs[6050]=['18 Sco',4.2622745,-0.1470308,5.5];
ybs[6051]=['',4.2635395,-0.2601209,6.09];
ybs[6052]=['',4.2792966,-1.0114664,6.49];
ybs[6053]=['σ CrB',4.2565611,0.5899781,5.64];
ybs[6054]=['σ CrB',4.2565611,0.5899733,6.66];
ybs[6055]=['16 Her',4.2606908,0.3273087,5.69];
ybs[6056]=['',4.2686811,-0.3727673,6.61];
ybs[6057]=['',4.2677938,-0.0699441,6.18];
ybs[6058]=['',4.2617,0.4776519,6.14];
ybs[6059]=['',4.2433873,1.1708982,6.21];
ybs[6060]=['',4.2747433,-0.5003392,4.78];
ybs[6061]=['λ Nor',4.279841,-0.7457232,5.45];
ybs[6062]=['γ2 Nor',4.2827702,-0.8762976,4.02];
ybs[6063]=['',4.285776,-0.963287,5.77];
ybs[6064]=['υ CrB',4.2657955,0.5078202,5.78];
ybs[6065]=['ε Oph',4.2739285,-0.0828329,3.24];
ybs[6066]=['',4.2780155,-0.3537926,6.29];
ybs[6067]=['',4.2802912,-0.5403448,5.49];
ybs[6068]=['',4.2772809,-0.2605011,5.94];
ybs[6069]=['19 UMi',4.2333546,1.3233071,5.48];
ybs[6070]=['',4.2850943,-0.6891105,6.12];
ybs[6071]=['ο Sco',4.284759,-0.422749,4.55];
ybs[6072]=['20 UMi',4.2411555,1.3116815,6.39];
ybs[6073]=['',4.2941947,-0.8660937,5.33];
ybs[6074]=['σ Sco',4.2872272,-0.4475859,2.89];
ybs[6075]=['',4.2938386,-0.7673084,5.88];
ybs[6076]=['',4.2657343,1.0419755,5.4];
ybs[6077]=['',4.2806354,0.3679126,6.05];
ybs[6078]=['',4.2507897,1.2800122,5.98];
ybs[6079]=['',4.3084683,-1.1026025,6.15];
ybs[6080]=['',4.2752616,0.8549472,5.91];
ybs[6081]=['',4.279087,0.6921247,5.46];
ybs[6082]=['τ Her',4.2778811,0.8073969,3.89];
ybs[6083]=['σ Ser',4.2901001,0.0170613,4.82];
ybs[6084]=['',4.3002687,-0.6849303,5.4];
ybs[6085]=['γ Her',4.2887754,0.3333809,3.75];
ybs[6086]=['',4.2927172,-0.0371939,6.23];
ybs[6087]=['',4.2996116,-0.580323,6.47];
ybs[6088]=['ζ TrA',4.3236176,-1.2240437,4.91];
ybs[6089]=['',4.3045134,-0.7923666,6.33];
ybs[6090]=['',4.3023883,-0.6565255,5.42];
ybs[6091]=['',4.2680271,1.195561,6.41];
ybs[6092]=['γ Aps',4.3501839,-1.3778067,3.89];
ybs[6093]=['ξ CrB',4.2890465,0.5382641,4.85];
ybs[6094]=['ψ Oph',4.299723,-0.3506029,4.5];
ybs[6095]=['',4.3025651,-0.5192983,6.63];
ybs[6096]=['',4.3025723,-0.5193225,5.84];
ybs[6097]=['ν1 CrB',4.2900388,0.5890069,5.2];
ybs[6098]=['ν2 CrB',4.2906105,0.5873403,5.39];
ybs[6099]=['ι TrA',4.3199476,-1.1188694,5.27];
ybs[6100]=['',4.2926624,0.5634236,6.4];
ybs[6101]=['21 Her',4.2990843,0.1203835,5.85];
ybs[6102]=['ρ Oph',4.3063359,-0.4101007,5.02];
ybs[6103]=['ρ Oph',4.3063285,-0.4100813,5.92];
ybs[6104]=['',4.3204077,-1.0236017,5.69];
ybs[6105]=['ε Nor',4.3146412,-0.8308454,4.47];
ybs[6106]=['η UMi',4.2624576,1.3212309,4.95];
ybs[6107]=['ω Her',4.3042229,0.2440556,4.57];
ybs[6108]=['χ Oph',4.3124056,-0.3229819,4.42];
ybs[6109]=['',4.3056883,0.3288671,6.7];
ybs[6110]=['',4.3268398,-1.00886,6.06];
ybs[6111]=['',4.3077004,0.198233,6.11];
ybs[6112]=['',4.3185946,-0.6497492,5.79];
ybs[6113]=['25 Her',4.3031343,0.651773,5.54];
ybs[6114]=['',4.3108302,0.0401121,6.07];
ybs[6115]=['',4.3320709,-1.0765299,5.2];
ybs[6116]=['',4.2837903,1.2052789,5.25];
ybs[6117]=['',4.2975063,0.9626256,5.74];
ybs[6118]=['',4.3150645,-0.1334624,5.23];
ybs[6119]=['υ Oph',4.3154269,-0.1469637,4.63];
ybs[6120]=['',4.293892,1.0759202,5.67];
ybs[6121]=['',4.3255476,-0.8079304,5.35];
ybs[6122]=['η Dra',4.2948275,1.0727367,2.74];
ybs[6123]=['',4.5763544,-1.5273002,6.57];
ybs[6124]=['α Sco',4.323138,-0.4621603,0.96];
ybs[6125]=['',4.3495054,-1.2397626,5.5];
ybs[6126]=['',4.318445,0.0107622,5.39];
ybs[6127]=['',4.3198402,-0.1427176,6.48];
ybs[6128]=['',4.4116705,-1.4534639,6.57];
ybs[6129]=['',4.4934509,-1.5048255,6.04];
ybs[6130]=['',4.3242945,-0.254793,5.68];
ybs[6131]=['22 Sco',4.3265728,-0.4391682,4.79];
ybs[6132]=['',4.3339339,-0.7306581,5.33];
ybs[6133]=['',4.3321557,-0.6065247,4.23];
ybs[6134]=['',4.3271639,-0.1319884,6.5];
ybs[6135]=['',4.3317565,-0.4639899,6.1];
ybs[6136]=['30 Her',4.3169885,0.7301276,5.04];
ybs[6137]=['φ Oph',4.330294,-0.2907686,4.28];
ybs[6138]=['β Her',4.3248863,0.374236,2.77];
ybs[6139]=['λ Oph',4.3286375,0.0338018,3.82];
ybs[6140]=['',4.316627,0.8963896,6.29];
ybs[6141]=['θ TrA',4.3542848,-1.1438842,5.52];
ybs[6142]=['',4.3264103,0.3566015,5.25];
ybs[6143]=['ω Oph',4.3348438,-0.3754711,4.45];
ybs[6144]=['',4.3292338,0.3865591,5.76];
ybs[6145]=['μ Nor',4.3445105,-0.7695285,4.94];
ybs[6146]=['34 Her',4.3228378,0.8536947,6.45];
ybs[6147]=['',4.3278599,0.6139684,6.25];
ybs[6148]=['28 Her',4.3358489,0.0955528,5.63];
ybs[6149]=['29 Her',4.3356762,0.1996955,4.84];
ybs[6150]=['',4.3491727,-0.7904534,6.46];
ybs[6151]=['15 Dra',4.3107492,1.1993745,5];
ybs[6152]=['',4.3304561,0.7950232,5.65];
ybs[6153]=['β Aps',4.3911429,-1.3536408,4.24];
ybs[6154]=['',4.3544369,-0.7488018,5.47];
ybs[6155]=['τ Sco',4.3514898,-0.4932423,2.82];
ybs[6156]=['',4.3539845,-0.616099,4.16];
ybs[6157]=['',4.367168,-1.06523,6.18];
ybs[6158]=['σ Her',4.340766,0.7398673,4.2];
ybs[6159]=['',4.3478255,0.2969207,6.41];
ybs[6160]=['',4.3316865,1.0607536,5.94];
ybs[6161]=['12 Oph',4.3525489,-0.0413492,5.75];
ybs[6162]=['η1 TrA',4.3796003,-1.1927168,5.91];
ybs[6163]=['',4.2958006,1.3773004,5.56];
ybs[6164]=['',4.3634753,-0.7582033,5.83];
ybs[6165]=['ζ Oph',4.3563384,-0.1852004,2.56];
ybs[6166]=['',4.3534624,0.2697196,6.3];
ybs[6167]=['',4.3756641,-1.055715,6.18];
ybs[6168]=['',4.3659309,-0.6503171,5.91];
ybs[6169]=['',4.3599803,-0.1148706,6.09];
ybs[6170]=['',4.3246726,1.2664924,6.3];
ybs[6171]=['',4.3582574,0.2381194,6.31];
ybs[6172]=['',4.3880348,-1.1776262,6.03];
ybs[6173]=['',4.3495565,0.8127772,5.79];
ybs[6174]=['16 Dra',4.3490392,0.9225046,5.53];
ybs[6175]=['17 Dra',4.3491966,0.9229267,5.08];
ybs[6176]=['17 Dra',4.3492257,0.922922,6.53];
ybs[6177]=['',4.3765479,-0.851804,5.65];
ybs[6178]=['',4.3780713,-0.8673102,5.65];
ybs[6179]=['',4.3671817,-0.1675018,6.35];
ybs[6180]=['',4.3716305,-0.3569341,6.26];
ybs[6181]=['',4.3184871,1.3508643,6.34];
ybs[6182]=['',4.3773749,-0.579239,5.87];
ybs[6183]=['',4.3762899,-0.4277753,6.09];
ybs[6184]=['36 Her',4.3707568,0.0726924,6.93];
ybs[6185]=['37 Her',4.3710181,0.0729111,5.77];
ybs[6186]=['',4.3758803,-0.310388,4.96];
ybs[6187]=['',4.3838391,-0.8047958,6.23];
ybs[6188]=['',4.3508799,1.1000533,6.16];
ybs[6189]=['',4.3566209,0.9768924,5.29];
ybs[6190]=['42 Her',4.3605414,0.8532043,4.9];
ybs[6191]=['',4.3735913,-0.0181948,6.24];
ybs[6192]=['',4.3773718,-0.348472,5.57];
ybs[6193]=['',4.371641,0.2155983,6.08];
ybs[6194]=['',4.4023717,-1.1719629,5.13];
ybs[6195]=['14 Oph',4.3757671,0.019887,5.74];
ybs[6196]=['',4.386561,-0.7183669,6.2];
ybs[6197]=['',4.3914694,-0.9283838,5.96];
ybs[6198]=['',4.3718087,0.4331301,6.06];
ybs[6199]=['',4.3871718,-0.7182638,6.12];
ybs[6200]=['',4.3865331,-0.6666663,6.05];
ybs[6201]=['',4.3855579,-0.5610657,6.46];
ybs[6202]=['ζ Her',4.3727026,0.5508451,2.81];
ybs[6203]=['39 Her',4.3743417,0.4690601,5.92];
ybs[6204]=['',4.3906982,-0.7134859,5.71];
ybs[6205]=['',4.3994694,-1.0217624,5.74];
ybs[6206]=['',4.3881428,-0.4799024,6.58];
ybs[6207]=['α TrA',4.4116071,-1.2054206,1.92];
ybs[6208]=['',4.3913126,-0.4982848,6.02];
ybs[6209]=['',4.4036922,-1.0189223,5.58];
ybs[6210]=['η Her',4.3793277,0.6786027,3.53];
ybs[6211]=['',4.3997176,-0.6879417,5.48];
ybs[6212]=['',4.3838087,0.5933814,5.99];
ybs[6213]=['18 Dra',4.3680347,1.1265548,4.83];
ybs[6214]=['16 Oph',4.3922951,0.0171141,6.03];
ybs[6215]=['25 Sco',4.3992587,-0.446238,6.71];
ybs[6216]=['',4.3783052,0.9712597,6.16];
ybs[6217]=['',4.3912241,0.2741121,5.56];
ybs[6218]=['43 Her',4.393488,0.1491026,5.15];
ybs[6219]=['η Ara',4.4145319,-1.0311171,3.76];
ybs[6220]=['',4.3890769,0.753585,6.05];
ybs[6221]=['',4.4270876,-1.1818987,6.32];
ybs[6222]=['19 Oph',4.3995366,0.0353534,6.1];
ybs[6223]=['',4.4248424,-1.1416488,6.13];
ybs[6224]=['45 Her',4.4020821,0.090899,5.24];
ybs[6225]=['',4.4057801,-0.2608846,6.03];
ybs[6226]=['',4.4170943,-0.8741038,6.47];
ybs[6227]=['',4.3883045,0.9903338,4.85];
ybs[6228]=['',4.3486904,1.3766126,6.32];
ybs[6229]=['',4.4034062,0.2365304,6.35];
ybs[6230]=['',4.4102321,-0.2741059,6.1];
ybs[6231]=['ε Sco',4.4141255,-0.5991808,2.29];
ybs[6232]=['',4.3984522,0.7365293,5.87];
ybs[6233]=['20 Oph',4.4116589,-0.1888531,4.65];
ybs[6234]=['',4.4179425,-0.655392,6.11];
ybs[6235]=['',4.4206364,-0.7202449,5.22];
ybs[6236]=['',4.4096502,0.2307984,5.91];
ybs[6237]=['μ1 Sco',4.4217862,-0.6646875,3.08];
ybs[6238]=['',4.4137073,-0.0469678,6.32];
ybs[6239]=['',4.4239806,-0.731127,6.49];
ybs[6240]=['47 Her',4.4131341,0.1258481,5.49];
ybs[6241]=['',4.4328271,-1.0113226,5.94];
ybs[6242]=['μ2 Sco',4.4238144,-0.6641597,3.57];
ybs[6243]=['',4.4397795,-1.1048633,6.02];
ybs[6244]=['52 Her',4.4065023,0.8018992,4.82];
ybs[6245]=['21 Oph',4.4181185,0.0205858,5.51];
ybs[6246]=['',4.4085974,0.757349,6.13];
ybs[6247]=['',4.4301295,-0.7519951,5.96];
ybs[6248]=['50 Her',4.4136289,0.519577,5.72];
ybs[6249]=['',4.4137902,0.5675207,6.13];
ybs[6250]=['',4.4314501,-0.7302725,5.45];
ybs[6251]=['',4.4312443,-0.7335601,6.32];
ybs[6252]=['ζ1 Sco',4.4313342,-0.739974,4.73];
ybs[6253]=['',4.4321806,-0.731037,6.45];
ybs[6254]=['',4.4127415,0.7305862,6.29];
ybs[6255]=['',4.4327461,-0.7305074,6.59];
ybs[6256]=['',4.4333213,-0.742006,5.88];
ybs[6257]=['',4.3725815,1.3521531,5.98];
ybs[6258]=['49 Her',4.4205121,0.2607147,6.52];
ybs[6259]=['',4.4276917,-0.3569388,5.88];
ybs[6260]=['51 Her',4.4186866,0.429698,5.04];
ybs[6261]=['ζ2 Sco',4.4339027,-0.739954,3.62];
ybs[6262]=['',4.4355255,-0.7188272,5.77];
ybs[6263]=['',4.4333014,-0.5344566,6.35];
ybs[6264]=['',4.4413854,-0.8850389,6.33];
ybs[6265]=['',4.4429789,-0.9131161,5.94];
ybs[6266]=['',4.4593051,-1.2095192,5.79];
ybs[6267]=['',4.4302719,-0.0287526,6.25];
ybs[6268]=['',4.4328112,-0.206427,6.57];
ybs[6269]=['53 Her',4.4236471,0.5526721,5.32];
ybs[6270]=['23 Oph',4.4322645,-0.1080156,5.25];
ybs[6271]=['ι Oph',4.4291023,0.1768016,4.38];
ybs[6272]=['',4.4394015,-0.5854076,6.37];
ybs[6273]=['',4.4426058,-0.7130963,6.15];
ybs[6274]=['',4.4389349,-0.2939184,6.37];
ybs[6275]=['ζ Ara',4.4527299,-0.9777846,3.13];
ybs[6276]=['',4.4240795,0.8269573,6];
ybs[6277]=['',4.432652,0.3651887,5.41];
ybs[6278]=['27 Sco',4.444721,-0.5810718,5.48];
ybs[6279]=['',4.4507773,-0.8844276,5.55];
ybs[6280]=['',4.4344599,0.2371044,6.34];
ybs[6281]=['24 Oph',4.4425597,-0.4046327,5.58];
ybs[6282]=['56 Her',4.4329461,0.4484757,6.08];
ybs[6283]=['54 Her',4.4347195,0.3211186,5.35];
ybs[6284]=['',4.4435704,-0.341624,6.27];
ybs[6285]=['ε1 Ara',4.4566326,-0.9283881,4.06];
ybs[6286]=['',4.4448382,-0.19193,6.19];
ybs[6287]=['',4.4590487,-0.9534528,5.65];
ybs[6288]=['',4.4523331,-0.6571815,6.09];
ybs[6289]=['κ Oph',4.4451003,0.1630422,3.2];
ybs[6290]=['',4.4599964,-0.8496177,6];
ybs[6291]=['',4.4443383,0.2417406,6.37];
ybs[6292]=['',4.4504808,-0.2600976,6.59];
ybs[6293]=['',4.4601201,-0.7938345,6.65];
ybs[6294]=['',4.4669537,-1.0295573,6.11];
ybs[6295]=['57 Her',4.4437883,0.441905,6.28];
ybs[6296]=['',4.4361257,0.8727444,6.56];
ybs[6297]=['',4.4446532,0.4249528,6.32];
ybs[6298]=['',4.4564372,-0.4384969,5.86];
ybs[6299]=['26 Oph',4.4572982,-0.4367012,5.75];
ybs[6300]=['',4.4598386,-0.6277227,5.97];
ybs[6301]=['',4.4659665,-0.8929428,6.45];
ybs[6302]=['',4.4442514,0.7414006,6.34];
ybs[6303]=['ε2 Ara',4.4721924,-0.9296885,5.29];
ybs[6304]=['19 Dra',4.4337601,1.1362136,4.89];
ybs[6305]=['',4.4651488,-0.5615536,5.03];
ybs[6306]=['',4.4575159,0.1143494,6.59];
ybs[6307]=['30 Oph',4.4604032,-0.0742471,4.82];
ybs[6308]=['20 Dra',4.4354866,1.1345494,6.41];
ybs[6309]=['',4.4782349,-1.0077845,5.73];
ybs[6310]=['29 Oph',4.4644332,-0.3301578,6.26];
ybs[6311]=['ε UMi',4.3796374,1.4311155,4.23];
ybs[6312]=['',4.474034,-0.8236213,6.06];
ybs[6313]=['ε Her',4.4556213,0.5392081,3.92];
ybs[6314]=['',4.4589652,0.3944543,5.65];
ybs[6315]=['',4.4618199,0.2603702,6.31];
ybs[6316]=['',4.474074,-0.6664051,5.91];
ybs[6317]=['',4.4595963,0.4741155,6.55];
ybs[6318]=['',4.4639685,0.1469475,6.33];
ybs[6319]=['',4.4496117,0.9888331,6.03];
ybs[6320]=['',4.4799876,-0.7946651,6.28];
ybs[6321]=['59 Her',4.4612279,0.5853307,5.25];
ybs[6322]=['',4.4646978,0.4446156,5.75];
ybs[6323]=['',4.4781078,-0.5960692,4.87];
ybs[6324]=['',4.4324675,1.2757185,6.3];
ybs[6325]=['',4.4642819,0.5559526,6.36];
ybs[6326]=['',4.4687468,0.2454186,4.98];
ybs[6327]=['',4.4830351,-0.7702821,6.19];
ybs[6328]=['',4.4689194,0.2527349,6.52];
ybs[6329]=['',4.4771637,-0.3582161,6.3];
ybs[6330]=['',4.4710563,0.2369296,5.93];
ybs[6331]=['',4.4724175,0.2362731,6.08];
ybs[6332]=['',4.4717849,0.3431395,6.35];
ybs[6333]=['',4.4848954,-0.6502426,5.98];
ybs[6334]=['',4.4459032,1.2069551,6.4];
ybs[6335]=['61 Her',4.4694002,0.617564,6.69];
ybs[6336]=['',4.4853829,-0.6192376,6.13];
ybs[6337]=['',4.4574401,1.0579749,6.13];
ybs[6338]=['',4.4786414,0.0117492,6.01];
ybs[6339]=['',4.483501,-0.3768777,6.3];
ybs[6340]=['',4.4711285,0.6066788,6.04];
ybs[6341]=['',4.4753224,0.341552,6.17];
ybs[6342]=['',4.4798205,-0.0160767,5.64];
ybs[6343]=['',4.4867329,-0.4632358,6.29];
ybs[6344]=['60 Her',4.4786087,0.2218582,4.91];
ybs[6345]=['',4.5036793,-1.0769041,6.39];
ybs[6346]=['',4.5155482,-1.2347557,6.22];
ybs[6347]=['',4.4821423,0.1693748,6.37];
ybs[6348]=['',4.4823619,0.1819562,6.37];
ybs[6349]=['',4.4610044,1.126948,6.1];
ybs[6350]=['',4.4857103,-0.0294062,6.38];
ybs[6351]=['',4.4757643,0.7641518,6.43];
ybs[6352]=['',4.4742709,0.8512749,6.09];
ybs[6353]=['',4.4822549,0.3849383,5.56];
ybs[6354]=['',4.4922887,-0.3078213,5.99];
ybs[6355]=['',4.4952348,-0.5311208,5.97];
ybs[6356]=['',4.4915576,-0.0193242,6.06];
ybs[6357]=['',4.5187179,-1.1732347,5.89];
ybs[6358]=['μ Dra',4.4759037,0.9501713,5.83];
ybs[6359]=['μ Dra',4.4758964,0.9501713,5.8];
ybs[6360]=['',4.5044365,-0.7781341,5.08];
ybs[6361]=['',4.4946365,-0.0682451,6.36];
ybs[6362]=['',4.5357986,-1.301245,6.25];
ybs[6363]=['',4.5088981,-0.8534649,5.84];
ybs[6364]=['',4.4987791,-0.1841361,5.56];
ybs[6365]=['',4.4877409,0.7066491,6.34];
ybs[6366]=['',4.4891294,0.6267011,5.39];
ybs[6367]=['η Oph',4.5015074,-0.2749119,2.43];
ybs[6368]=['',4.4548568,1.3136308,6.21];
ybs[6369]=['η Sco',4.5106638,-0.7551114,3.33];
ybs[6370]=['',4.5109251,-0.6899711,5.67];
ybs[6371]=['',4.5109043,-0.6780204,6.3];
ybs[6372]=['',4.4891554,0.8868776,6.46];
ybs[6373]=['',4.5209576,-0.9933136,6.09];
ybs[6374]=['',4.5021133,0.2171329,6.57];
ybs[6375]=['',4.5099805,-0.4412241,6.54];
ybs[6376]=['',4.5109267,-0.4849814,6.14];
ybs[6377]=['',4.4954729,0.7112229,5.08];
ybs[6378]=['',4.513612,-0.5665992,6.01];
ybs[6379]=['',4.5066088,0.1373371,6.33];
ybs[6380]=['63 Her',4.5028889,0.42257,6.19];
ybs[6381]=['',4.5205078,-0.6944883,6.6];
ybs[6382]=['37 Oph',4.5096033,0.1843026,5.33];
ybs[6383]=['',4.5119246,0.0057019,6.65];
ybs[6384]=['',4.4987066,0.9142417,6.29];
ybs[6385]=['ζ Dra',4.4892384,1.1464541,3.17];
ybs[6386]=['',4.5238985,-0.5859455,5.53];
ybs[6387]=['',4.525396,-0.6740041,5.96];
ybs[6388]=['',4.5040168,0.8677882,6.04];
ybs[6389]=['',4.5497587,-1.2228914,6.53];
ybs[6390]=['36 Oph',4.5236737,-0.4647228,5.11];
ybs[6391]=['36 Oph',4.5236592,-0.4646985,5.07];
ybs[6392]=['',4.5260671,-0.5276854,6.21];
ybs[6393]=['',4.5231063,-0.2549543,5.99];
ybs[6394]=['',4.5285451,-0.6243522,6.12];
ybs[6395]=['α1 Her',4.518987,0.2507327,3.48];
ybs[6396]=['α2 Her',4.5190088,0.2507279,5.39];
ybs[6397]=['',4.5431121,-1.0422423,5.91];
ybs[6398]=['',4.5314454,-0.5704733,5.55];
ybs[6399]=['δ Her',4.5202169,0.4331033,3.14];
ybs[6400]=['ι Aps',4.5580168,-1.224231,5.41];
ybs[6401]=['',4.5263856,0.0377452,6.17];
ybs[6402]=['',4.5287826,-0.1094006,6.09];
ybs[6403]=['',4.5276949,0.0207213,5.88];
ybs[6404]=['41 Oph',4.5281206,-0.0081776,4.73];
ybs[6405]=['',4.5410373,-0.8142956,5.48];
ybs[6406]=['ζ Aps',4.556854,-1.1831691,4.78];
ybs[6407]=['π Her',4.5196703,0.6420188,3.16];
ybs[6408]=['',4.5231529,0.413974,5.96];
ybs[6409]=['',4.5397406,-0.7705922,5.76];
ybs[6410]=['',4.5061774,1.0969175,5.56];
ybs[6411]=['',4.5370187,-0.5685511,6.36];
ybs[6412]=['',4.5433013,-0.8741461,6.27];
ybs[6413]=['ο Oph',4.5351831,-0.424279,5.2];
ybs[6414]=['ο Oph',4.5351685,-0.4242306,6.8];
ybs[6415]=['',4.5398281,-0.6110683,5.91];
ybs[6416]=['',4.5424027,-0.7722154,6.65];
ybs[6417]=['',4.5361875,-0.2850864,6.43];
ybs[6418]=['',4.6063776,-1.4115102,5.88];
ybs[6419]=['',4.5314925,0.4026133,6.45];
ybs[6420]=['68 Her',4.5298155,0.5773032,4.82];
ybs[6421]=['',4.5338546,0.3018642,6];
ybs[6422]=['',4.5364393,0.1892326,5.03];
ybs[6423]=['',4.5377718,0.1058231,6.51];
ybs[6424]=['',4.5430907,-0.3102819,6.02];
ybs[6425]=['69 Her',4.531092,0.6504646,4.65];
ybs[6426]=['',4.5264378,0.8668666,7.48];
ybs[6427]=['',4.5590722,-1.0128133,5.88];
ybs[6428]=['',4.5430751,-0.1036491,6.32];
ybs[6429]=['',4.5646391,-1.0975185,5.7];
ybs[6430]=['',4.5461314,-0.3377886,6.52];
ybs[6431]=['',4.5597464,-0.9868936,5.8];
ybs[6432]=['',4.5365046,0.5026705,5.65];
ybs[6433]=['',4.5341261,0.6769953,5.94];
ybs[6434]=['ξ Oph',4.5480985,-0.3688514,4.39];
ybs[6435]=['ν Ser',4.546994,-0.2245874,4.33];
ybs[6436]=['',4.5655497,-1.0592838,5.77];
ybs[6437]=['',4.5237527,1.0584897,6.32];
ybs[6438]=['',4.5471278,-0.1870479,6.46];
ybs[6439]=['',4.5561549,-0.6601691,6.41];
ybs[6440]=['ι Ara',4.5595034,-0.8288196,5.25];
ybs[6441]=['',4.5435603,0.314786,5];
ybs[6442]=['θ Oph',4.5526645,-0.4366768,3.27];
ybs[6443]=['',4.5559414,-0.6270955,6.47];
ybs[6444]=['',4.5425718,0.4453397,5.38];
ybs[6445]=['',4.5572473,-0.6499662,5.93];
ybs[6446]=['70 Her',4.5458476,0.4272292,5.12];
ybs[6447]=['72 Her',4.544397,0.5663001,5.39];
ybs[6448]=['43 Oph',4.5587041,-0.4915304,5.35];
ybs[6449]=['',4.5633857,-0.7711133,5.12];
ybs[6450]=['β Ara',4.5691646,-0.9695022,2.85];
ybs[6451]=['γ Ara',4.5696742,-0.9842929,3.34];
ybs[6452]=['',4.5490322,0.2916478,6.35];
ybs[6453]=['74 Her',4.5421735,0.8066815,5.59];
ybs[6454]=['',4.5554354,-0.0420316,6.29];
ybs[6455]=['',4.5483372,0.5015616,6.35];
ybs[6456]=['',4.5429235,0.8406735,6.43];
ybs[6457]=['κ Ara',4.5717056,-0.8840381,5.23];
ybs[6458]=['',4.5486147,0.6973258,5.51];
ybs[6459]=['',4.5663981,-0.6058913,6.16];
ybs[6460]=['',4.5826146,-1.1004855,6.24];
ybs[6461]=['',4.5642447,-0.3745568,5.85];
ybs[6462]=['',4.5637603,-0.3222705,6.21];
ybs[6463]=['',4.5661293,-0.4234559,6.19];
ybs[6464]=['',4.5759285,-0.9069899,6.19];
ybs[6465]=['',4.5598305,0.1541728,5.77];
ybs[6466]=['',4.5750429,-0.8004292,5.29];
ybs[6467]=['',4.5769723,-0.8839685,5.92];
ybs[6468]=['',4.5476482,0.9320038,5.67];
ybs[6469]=['73 Her',4.5599063,0.4003959,5.74];
ybs[6470]=['',4.5620003,0.2841758,5.71];
ybs[6471]=['',4.5621959,0.2720461,6.35];
ybs[6472]=['',4.580423,-0.9130549,5.75];
ybs[6473]=['ρ Her',4.5573022,0.6479902,5.47];
ybs[6474]=['ρ Her',4.557324,0.6479758,4.52];
ybs[6475]=['44 Oph',4.5716538,-0.4222514,4.17];
ybs[6476]=['',4.5837241,-0.9631826,5.94];
ybs[6477]=['',4.5587784,0.6730584,6.49];
ybs[6478]=['',4.5689777,-0.029145,6.44];
ybs[6479]=['',4.5741355,-0.4531045,6.44];
ybs[6480]=['',4.5606871,0.644599,6.28];
ybs[6481]=['45 Oph',4.5762238,-0.5215801,4.29];
ybs[6482]=['',4.5720169,-0.0890906,4.54];
ybs[6483]=['',4.5773949,-0.5190904,6];
ybs[6484]=['',4.5680031,0.2949468,5.98];
ybs[6485]=['',4.5740644,-0.2186917,6.21];
ybs[6486]=['',4.5701613,0.1322525,6.06];
ybs[6487]=['σ Oph',4.5711585,0.0719485,4.34];
ybs[6488]=['',4.5680373,0.4688063,6.41];
ybs[6489]=['δ Ara',4.5951792,-1.059399,3.62];
ybs[6490]=['',4.5834945,-0.6421912,6.02];
ybs[6491]=['',4.5718476,0.3501659,5.54];
ybs[6492]=['',4.5857525,-0.6725261,6.39];
ybs[6493]=['',4.5782744,-0.1435605,6.37];
ybs[6494]=['',4.5958898,-0.9937143,5.95];
ybs[6495]=['',4.5709299,0.6052444,5.94];
ybs[6496]=['',4.5813964,0.0054783,5.44];
ybs[6497]=['υ Sco',4.5915046,-0.651206,2.69];
ybs[6498]=['77 Her',4.5698577,0.8419825,5.85];
ybs[6499]=['α Ara',4.5971332,-0.8707621,2.95];
ybs[6500]=['',4.5639572,1.0477164,5.65];
ybs[6501]=['',4.5858298,-0.1036003,6.37];
ybs[6502]=['',4.5967063,-0.8037469,6.03];
ybs[6503]=['',4.5658752,1.0233487,6.51];
ybs[6504]=['',4.5882924,-0.0188202,5.31];
ybs[6505]=['',4.5957795,-0.5884858,6.44];
ybs[6506]=['',4.5595542,1.1743854,6.43];
ybs[6507]=['51 Oph',4.5936673,-0.4184949,4.81];
ybs[6508]=['',4.5951849,-0.4587554,6.05];
ybs[6509]=['',4.5876924,0.2078536,6.39];
ybs[6510]=['',4.5985098,-0.5985495,6.17];
ybs[6511]=['',4.6020494,-0.7188634,5.84];
ybs[6512]=['',4.5923304,0.0472833,5.59];
ybs[6513]=['',4.6146203,-1.0447347,6.28];
ybs[6514]=['λ Her',4.5886799,0.4554411,4.41];
ybs[6515]=['λ Sco',4.6039091,-0.647829,1.63];
ybs[6516]=['',4.5892387,0.5435429,5.61];
ybs[6517]=['',4.5288584,1.3982527,5.72];
ybs[6518]=['',4.6127157,-0.9314078,6.1];
ybs[6519]=['',4.5877033,0.6783473,6.43];
ybs[6520]=['',4.5958802,0.2079635,6.42];
ybs[6521]=['78 Her',4.5933179,0.4955405,5.62];
ybs[6522]=['',4.6020041,-0.1005108,5.62];
ybs[6523]=['',4.6084546,-0.5688912,5.7];
ybs[6524]=['β Dra',4.5856069,0.9125525,2.79];
ybs[6525]=['σ Ara',4.6135101,-0.8118989,4.59];
ybs[6526]=['',4.5938491,0.5978765,6.56];
ybs[6527]=['',4.6131303,-0.6536752,6.48];
ybs[6528]=['',4.5862533,1.0098568,6.4];
ybs[6529]=['',4.6005105,0.3358435,5.64];
ybs[6530]=['',4.6018406,0.2845481,5.69];
ybs[6531]=['',4.6021481,0.2587905,6.48];
ybs[6532]=['',4.6077809,-0.1964431,5.55];
ybs[6533]=['52 Oph',4.6105637,-0.3849671,6.57];
ybs[6534]=['',4.6168255,-0.6745288,4.29];
ybs[6535]=['',4.6216422,-0.8739182,5.93];
ybs[6536]=['53 Oph',4.6062764,0.1670823,5.81];
ybs[6537]=['π Ara',4.6248877,-0.9514092,5.25];
ybs[6538]=['',4.5982383,0.7195844,5.74];
ybs[6539]=['',4.6053087,0.2878088,6.4];
ybs[6540]=['',4.7503663,-1.4820445,6.45];
ybs[6541]=['θ Sco',4.6204902,-0.7506611,1.87];
ybs[6542]=['ν1 Dra',4.5928984,0.9628827,4.88];
ybs[6543]=['ν2 Dra',4.5932923,0.9626897,4.87];
ybs[6544]=['α Oph',4.60757,0.2189797,2.08];
ybs[6545]=['',4.620718,-0.6645768,6.26];
ybs[6546]=['',4.6240678,-0.7486074,6.1];
ybs[6547]=['',4.6118375,0.3662273,6.1];
ybs[6548]=['',4.5985002,1.0043419,6.17];
ybs[6549]=['ξ Ser',4.6202238,-0.2689638,3.54];
ybs[6550]=['',4.6203036,-0.2719744,5.94];
ybs[6551]=['',4.6097649,0.6508093,6.1];
ybs[6552]=['',4.6121103,0.4916931,6.38];
ybs[6553]=['',4.6558464,-1.2606278,6.49];
ybs[6554]=['27 Dra',4.5896897,1.1889129,5.05];
ybs[6555]=['μ Oph',4.6210621,-0.1419065,4.62];
ybs[6556]=['',4.6225389,-0.1909036,5.75];
ybs[6557]=['λ Ara',4.6344146,-0.8626427,4.77];
ybs[6558]=['',4.6140742,0.5370858,6.02];
ybs[6559]=['79 Her',4.6183432,0.4240798,5.77];
ybs[6560]=['',4.6380343,-0.8191128,5.79];
ybs[6561]=['26 Dra',4.6042363,1.0796855,5.23];
ybs[6562]=['82 Her',4.6129699,0.8477581,5.37];
ybs[6563]=['',4.6263321,0.0352029,6.26];
ybs[6564]=['',4.6418108,-0.881738,6.24];
ybs[6565]=['',4.6251152,0.2324423,6.12];
ybs[6566]=['',4.6310944,-0.0377513,6.19];
ybs[6567]=['',4.6236423,0.5712134,6.37];
ybs[6568]=['κ Sco',4.6427795,-0.6813615,2.41];
ybs[6569]=['ο Ser',4.6368267,-0.2248871,4.26];
ybs[6570]=['η Pav',4.6597938,-1.1297705,3.62];
ybs[6571]=['',4.64423,-0.6449826,5.54];
ybs[6572]=['',4.6286534,0.5443995,6.03];
ybs[6573]=['μ Ara',4.6510164,-0.9048195,5.15];
ybs[6574]=['',4.6550972,-1.004489,6.01];
ybs[6575]=['',4.6451653,-0.5770046,6.4];
ybs[6576]=['ι Her',4.6255646,0.8027704,3.8];
ybs[6577]=['',4.6347265,0.2647325,6.34];
ybs[6578]=['',4.6366256,0.1100081,5.95];
ybs[6579]=['',4.6318266,0.5458899,6.28];
ybs[6580]=['',4.633924,0.4276626,6.36];
ybs[6581]=['',4.6456928,-0.486823,6.36];
ybs[6582]=['',4.6381634,0.2782471,5.52];
ybs[6583]=['58 Oph',4.6459878,-0.3785971,4.87];
ybs[6584]=['ω Dra',4.6112595,1.1998344,4.8];
ybs[6585]=['',4.6526775,-0.745898,5.87];
ybs[6586]=['',4.6097511,1.214017,6.42];
ybs[6587]=['',4.6308236,0.7585278,6.59];
ybs[6588]=['',4.6473027,-0.235918,6.39];
ybs[6589]=['',4.6469332,-0.1237084,6.3];
ybs[6590]=['83 Her',4.6399499,0.428563,5.52];
ybs[6591]=['β Oph',4.6451355,0.0795609,2.77];
ybs[6592]=['',4.6442886,0.2493411,6.24];
ybs[6593]=['',4.6293927,1.00007,6.77];
ybs[6594]=['',4.6108265,1.2643732,5.86];
ybs[6595]=['',4.6333323,0.9042206,5.99];
ybs[6596]=['84 Her',4.6438293,0.4244456,5.71];
ybs[6597]=['61 Oph',4.649982,0.0448781,6.17];
ybs[6598]=['',4.6500838,0.0448686,6.56];
ybs[6599]=['',4.6483053,0.2513619,6.19];
ybs[6600]=['',4.6415398,0.7692608,6.34];
ybs[6601]=['',4.6629412,-0.6652895,6.43];
ybs[6602]=['',4.6709979,-0.9670409,6.11];
ybs[6603]=['ι1 Sco',4.6650906,-0.7004582,3.03];
ybs[6604]=['3 Sgr',4.6642938,-0.4858516,4.54];
ybs[6605]=['',4.6649172,-0.3924263,6.18];
ybs[6606]=['',4.6445769,0.9388658,5.75];
ybs[6607]=['',4.6535642,0.5497285,6.23];
ybs[6608]=['',4.6639477,-0.2571263,5.94];
ybs[6609]=['',4.6681942,-0.4709059,6.35];
ybs[6610]=['',4.6788031,-0.9357919,5.92];
ybs[6611]=['μ Her',4.6571873,0.4836901,3.42];
ybs[6612]=['',4.6845427,-1.0501381,5.78];
ybs[6613]=['',4.6541165,0.6784774,6.52];
ybs[6614]=['',4.6544377,0.6861769,6.68];
ybs[6615]=['',4.6605829,0.3087569,5.72];
ybs[6616]=['',4.671537,-0.5534239,4.83];
ybs[6617]=['γ Oph',4.6644923,0.0471398,3.75];
ybs[6618]=['',4.6748182,-0.6466175,3.21];
ybs[6619]=['ι2 Sco',4.6764411,-0.6997983,4.81];
ybs[6620]=['',4.6818643,-0.9273785,6.09];
ybs[6621]=['',4.6663847,0.0662893,6.22];
ybs[6622]=['',4.6930094,-1.1430547,6.49];
ybs[6623]=['',4.7162564,-1.329555,6.07];
ybs[6624]=['ψ1 Dra',4.6318468,1.2590617,4.58];
ybs[6625]=['ψ1 Dra',4.6319671,1.2592026,5.79];
ybs[6626]=['',4.6660452,0.3588308,5.69];
ybs[6627]=['',4.670724,0.0341313,6.47];
ybs[6628]=['',4.6836185,-0.7959506,6.11];
ybs[6629]=['',4.6588951,0.8308702,6.43];
ybs[6630]=['',4.6677815,0.3359658,6.12];
ybs[6631]=['',4.6824205,-0.7116875,5.96];
ybs[6632]=['87 Her',4.6675915,0.4470995,5.12];
ybs[6633]=['',4.6803499,-0.533401,6.66];
ybs[6634]=['',4.7557721,-1.4221324,6.35];
ybs[6635]=['',4.685022,-0.607427,5.9];
ybs[6636]=['',4.6854441,-0.6007502,5.84];
ybs[6637]=['',4.6883142,-0.7330404,6.2];
ybs[6638]=['',4.6764888,0.2084252,6.17];
ybs[6639]=['',4.6875725,-0.5954659,6.06];
ybs[6640]=['',4.6881107,-0.6112503,6.45];
ybs[6641]=['',4.6882908,-0.6218189,6.03];
ybs[6642]=['',4.674226,0.5116817,5.5];
ybs[6643]=['',4.6764106,0.3894113,5.98];
ybs[6644]=['30 Dra',4.6670612,0.886196,5.02];
ybs[6645]=['',4.689818,-0.606219,6.17];
ybs[6646]=['',4.6900963,-0.6090934,5.6];
ybs[6647]=['',4.682525,-0.0216547,6.35];
ybs[6648]=['',4.6917047,-0.6071797,6.38];
ybs[6649]=['',4.685571,-0.1072906,6.21];
ybs[6650]=['',4.6923864,-0.6065964,5.96];
ybs[6651]=['',4.6926236,-0.6079777,6.42];
ybs[6652]=['88 Her',4.6715679,0.8445455,6.68];
ybs[6653]=['',4.6818063,0.2674146,6.46];
ybs[6654]=['',4.6875601,-0.1902962,6.18];
ybs[6655]=['',4.6850388,0.0227113,5.95];
ybs[6656]=['',4.6947193,-0.6015978,5.96];
ybs[6657]=['',4.6773311,0.699317,6.46];
ybs[6658]=['',4.6876766,0.10643,5.77];
ybs[6659]=['',4.6977948,-0.6366626,6.06];
ybs[6660]=['',4.6961818,-0.434406,6.2];
ybs[6661]=['',4.6810237,0.6977493,6.04];
ybs[6662]=['',4.6802729,0.8140065,6.38];
ybs[6663]=['',4.705551,-0.7739406,4.86];
ybs[6664]=['',4.6918435,0.1942149,6.38];
ybs[6665]=['90 Her',4.6863595,0.6982119,5.16];
ybs[6666]=['',4.7058797,-0.7034865,6.43];
ybs[6667]=['',4.7003481,-0.3281933,6.52];
ybs[6668]=['',4.70418,-0.4898562,5.8];
ybs[6669]=['',4.7019683,-0.2760092,5.89];
ybs[6670]=['',4.7097486,-0.7280969,4.88];
ybs[6671]=['',4.7103131,-0.6830805,6.29];
ybs[6672]=['',4.70128,0.0116688,5.82];
ybs[6673]=['89 Her',4.6963683,0.4546189,5.46];
ybs[6674]=['',4.7035923,-0.0712683,5.47];
ybs[6675]=['',4.6983945,0.3920386,5.58];
ybs[6676]=['ξ Dra',4.685834,0.9925573,3.75];
ybs[6677]=['',4.7046409,0.0011363,5.97];
ybs[6678]=['',4.7037913,0.1132091,6.29];
ybs[6679]=['',4.7143764,-0.6433022,5.74];
ybs[6680]=['',4.7127374,-0.5019482,6.01];
ybs[6681]=['',4.7147107,-0.5280172,5.16];
ybs[6682]=['',4.7147325,-0.5280221,7.04];
ybs[6683]=['θ Her',4.6994106,0.6501127,3.86];
ybs[6684]=['',4.7058603,0.1927378,6.36];
ybs[6685]=['',4.7043975,0.4187844,6.3];
ybs[6686]=['ν Oph',4.7135503,-0.1705851,3.34];
ybs[6687]=['',4.6940906,0.976843,6.1];
ybs[6688]=['4 Sgr',4.7174839,-0.415665,4.76];
ybs[6689]=['35 Dra',4.6619591,1.3431485,5.04];
ybs[6690]=['',4.7012742,0.791494,6.02];
ybs[6691]=['ξ Her',4.7064445,0.5104528,3.7];
ybs[6692]=['',4.7182384,-0.354979,6.21];
ybs[6693]=['γ Dra',4.6998675,0.8986207,2.23];
ybs[6694]=['',4.7159107,-0.0841471,5.87];
ybs[6695]=['ν Her',4.7096184,0.526895,4.41];
ybs[6696]=['',4.7269069,-0.6348875,6.3];
ybs[6697]=['',4.7185239,0.0109938,6.37];
ybs[6698]=['ζ Ser',4.7196706,-0.0643972,4.62];
ybs[6699]=['',4.710165,0.6333326,6];
ybs[6700]=['66 Oph',4.7184027,0.0762544,4.64];
ybs[6701]=['93 Her',4.7170099,0.2923625,4.67];
ybs[6702]=['67 Oph',4.7201231,0.0511787,3.97];
ybs[6703]=['6 Sgr',4.7241388,-0.2994258,6.28];
ybs[6704]=['',4.7266582,-0.397571,5.77];
ybs[6705]=['',4.6640014,1.3666095,6.24];
ybs[6706]=['',4.7102746,0.7937002,6.48];
ybs[6707]=['',4.7210276,0.1094167,6.34];
ybs[6708]=['',4.7186607,0.3404501,6.5];
ybs[6709]=['χ Oct',5.0065528,-1.5270999,5.28];
ybs[6710]=['',4.7209908,0.2634423,6.26];
ybs[6711]=['68 Oph',4.7250212,0.0228036,4.45];
ybs[6712]=['7 Sgr',4.7308496,-0.4237709,5.34];
ybs[6713]=['ψ2 Dra',4.6896902,1.2566765,5.45];
ybs[6714]=['',4.7186364,0.5797014,5.99];
ybs[6715]=['',4.7315495,-0.3964743,6.74];
ybs[6716]=['',4.7149195,0.7941516,5.67];
ybs[6717]=['95 Her',4.723093,0.3769275,5.18];
ybs[6718]=['95 Her',4.7231294,0.3769324,4.96];
ybs[6719]=['',4.7752415,-1.3244322,5.86];
ybs[6720]=['',4.7297194,-0.0934932,6.76];
ybs[6721]=['τ Oph',4.7311885,-0.1427375,5.94];
ybs[6722]=['τ Oph',4.7311812,-0.1427424,5.24];
ybs[6723]=['',4.6849212,1.311922,6.36];
ybs[6724]=['9 Sgr',4.735311,-0.4251284,5.97];
ybs[6725]=['',4.722958,0.5814126,6.15];
ybs[6726]=['96 Her',4.7269953,0.3636424,5.28];
ybs[6727]=['',4.7401213,-0.626544,6];
ybs[6728]=['',4.7560664,-1.126525,6.41];
ybs[6729]=['97 Her',4.7274118,0.4001112,6.21];
ybs[6730]=['γ1 Sgr',4.7405723,-0.5162136,4.69];
ybs[6731]=['θ Ara',4.7489692,-0.8741929,3.66];
ybs[6732]=['',4.7308007,0.342348,6.5];
ybs[6733]=['π Pav',4.759271,-1.1111299,4.35];
ybs[6734]=['γ2 Sgr',4.7440476,-0.5309396,2.99];
ybs[6735]=['',4.7375132,0.0335453,6.14];
ybs[6736]=['',4.7469128,-0.6285945,5.95];
ybs[6737]=['',4.7492902,-0.7578315,5.77];
ybs[6738]=['',4.7492902,-0.7578315,5.77];
ybs[6739]=['',4.7797365,-1.2856796,5.85];
ybs[6740]=['70 Oph',4.7411272,0.0436809,4.03];
ybs[6741]=['',4.7286973,0.845897,6.21];
ybs[6742]=['',4.7368269,0.4179245,6.34];
ybs[6743]=['',4.7444585,-0.145215,5.85];
ybs[6744]=['',4.7448876,-0.0828621,5.77];
ybs[6745]=['',4.7441558,-0.0077319,6.34];
ybs[6746]=['',4.7419215,0.209567,7.04];
ybs[6747]=['',4.7567598,-0.7986997,6.15];
ybs[6748]=['',4.7645191,-1.0303379,6.38];
ybs[6749]=['ι Pav',4.7670585,-1.0820334,5.49];
ybs[6750]=['',4.7496511,-0.3741915,6.28];
ybs[6751]=['',4.7405661,0.3778675,6.15];
ybs[6752]=['',4.736176,0.6996492,6.52];
ybs[6753]=['98 Her',4.7428461,0.3878549,5.06];
ybs[6754]=['',4.7538749,-0.4965883,4.57];
ybs[6755]=['',4.7373311,0.7321586,6.34];
ybs[6756]=['',4.7414664,0.5625889,5.71];
ybs[6757]=['',4.752158,-0.299316,5.52];
ybs[6758]=['71 Oph',4.7489663,0.1525099,4.64];
ybs[6759]=['72 Oph',4.7491229,0.1669966,3.73];
ybs[6760]=['',4.7599458,-0.6399592,6.58];
ybs[6761]=['',4.757302,-0.4444875,6.61];
ybs[6762]=['',4.7865046,-1.2346951,6.73];
ybs[6763]=['99 Her',4.7467826,0.5334777,5.04];
ybs[6764]=['',4.7509764,0.2282133,6.63];
ybs[6765]=['',4.7624495,-0.5709645,6.43];
ybs[6766]=['',4.768128,-0.8291456,6.07];
ybs[6767]=['ο Her',4.7491353,0.5020764,3.83];
ybs[6768]=['',4.7627664,-0.5362122,5.53];
ybs[6769]=['100 Her',4.7504932,0.4556342,5.86];
ybs[6770]=['100 Her',4.7504934,0.4555663,5.9];
ybs[6771]=['ε Tel',4.7686778,-0.8019414,4.53];
ybs[6772]=['',4.7542217,0.2494019,6.37];
ybs[6773]=['',4.7603885,-0.2431031,6.39];
ybs[6774]=['',4.7677464,-0.7217402,5.86];
ybs[6775]=['102 Her',4.7548056,0.3633687,4.36];
ybs[6776]=['',4.7665297,-0.5898052,6.16];
ybs[6777]=['δ UMi',4.5599653,1.5082385,4.36];
ybs[6778]=['',4.7448096,0.8870931,6.29];
ybs[6779]=['',4.7479875,0.7586238,5];
ybs[6780]=['',4.7458444,0.8676833,6.32];
ybs[6781]=['',4.7508764,0.6354045,5.48];
ybs[6782]=['101 Her',4.7553772,0.3499454,5.1];
ybs[6783]=['73 Oph',4.7589961,0.0697933,5.73];
ybs[6784]=['',4.7840961,-1.1114442,6.47];
ybs[6785]=['',4.7604985,0.0545492,5.69];
ybs[6786]=['',4.7672812,-0.3461985,6.36];
ybs[6787]=['',4.7561472,0.5318836,6.38];
ybs[6788]=['',4.7638577,0.0581248,5.51];
ybs[6789]=['11 Sgr',4.7695316,-0.4135441,4.98];
ybs[6790]=['',4.7708547,-0.5043035,6.51];
ybs[6791]=['',4.7610406,0.2876736,6.09];
ybs[6792]=['',4.7769795,-0.7213177,5.47];
ybs[6793]=['',4.7901337,-1.1003675,5.6];
ybs[6794]=['',4.7576828,0.6713054,6.4];
ybs[6795]=['',4.759367,0.6365575,5.58];
ybs[6796]=['',4.8014794,-1.1906405,6.33];
ybs[6797]=['40 Dra',4.7052683,1.3962704,6.04];
ybs[6798]=['41 Dra',4.7056858,1.3963295,5.68];
ybs[6799]=['24 UMi',4.5492896,1.5154381,5.79];
ybs[6800]=['μ Sgr',4.7783146,-0.3674093,3.86];
ybs[6801]=['',4.775024,-0.0698856,6.59];
ybs[6802]=['',4.7672547,0.5838753,5.88];
ybs[6803]=['104 Her',4.7680258,0.548243,4.97];
ybs[6804]=['14 Sgr',4.7805321,-0.3788219,5.44];
ybs[6805]=['',4.7603139,0.9475837,5.95];
ybs[6806]=['',4.7888236,-0.7713928,5.46];
ybs[6807]=['',4.7953761,-0.9776194,5.33];
ybs[6808]=['',4.7744672,0.382014,6.12];
ybs[6809]=['',4.7943603,-0.89114,6.06];
ybs[6810]=['15 Sgr',4.7846338,-0.3616263,5.38];
ybs[6811]=['16 Sgr',4.784619,-0.3556873,5.95];
ybs[6812]=['',4.770971,0.7182743,6.36];
ybs[6813]=['',4.7858468,-0.3255485,6.07];
ybs[6814]=['',4.7727396,0.6768555,6.04];
ybs[6815]=['',4.7620878,1.054451,6.49];
ybs[6816]=['',4.807688,-1.1148398,6.18];
ybs[6817]=['φ Oct',4.8285339,-1.3095287,5.47];
ybs[6818]=['',4.7872261,-0.0629795,6.36];
ybs[6819]=['',4.7804873,0.5099068,6.56];
ybs[6820]=['η Sgr',4.7959592,-0.6414369,3.11];
ybs[6821]=['',4.795698,-0.5951085,6.16];
ybs[6822]=['',4.7875504,0.0416589,6.01];
ybs[6823]=['',4.7945331,-0.4999079,6.19];
ybs[6824]=['',4.7944935,-0.4935666,6.4];
ybs[6825]=['',4.8579458,-1.400028,5.95];
ybs[6826]=['',4.7931225,-0.3030615,5.75];
ybs[6827]=['',4.8008504,-0.7378852,6.3];
ybs[6828]=['',4.7912098,-0.0523194,6];
ybs[6829]=['',4.7943973,-0.3220732,6.54];
ybs[6830]=['',4.7973052,-0.471802,4.65];
ybs[6831]=['',4.7937317,-0.1701482,6.31];
ybs[6832]=['',4.7919078,0.0177233,6.63];
ybs[6833]=['',4.7837171,0.7359736,5.59];
ybs[6834]=['',4.8000367,-0.4467019,6.51];
ybs[6835]=['',4.783044,0.789205,6.29];
ybs[6836]=['',4.7998432,-0.3247859,6.84];
ybs[6837]=['',4.7781884,0.9877946,6.37];
ybs[6838]=['36 Dra',4.7735015,1.1240766,5.03];
ybs[6839]=['',4.7956369,0.2406299,6.3];
ybs[6840]=['',4.7958068,0.3166299,5.99];
ybs[6841]=['',4.7902007,0.7146461,6.11];
ybs[6842]=['',4.795582,0.4067807,6.63];
ybs[6843]=['ξ Pav',4.8227681,-1.0730398,4.36];
ybs[6844]=['',4.8103777,-0.6540737,6.45];
ybs[6845]=['',4.800737,0.1268938,5.39];
ybs[6846]=['',4.8059443,-0.2761166,5.39];
ybs[6847]=['δ Sgr',4.8102762,-0.520391,2.7];
ybs[6848]=['105 Her',4.8001115,0.4268522,5.27];
ybs[6849]=['',4.8123444,-0.4346372,6.25];
ybs[6850]=['',4.816506,-0.6744711,5.1];
ybs[6851]=['',4.8114734,-0.3289591,5.75];
ybs[6852]=['',4.8146058,-0.4959808,6.16];
ybs[6853]=['37 Dra',4.7784944,1.2001608,5.95];
ybs[6854]=['74 Oph',4.8083474,0.0591477,4.86];
ybs[6855]=['',4.8028719,0.5179647,5.99];
ybs[6856]=['106 Her',4.8051107,0.3834964,4.95];
ybs[6857]=['η Ser',4.8105131,-0.0503866,3.26];
ybs[6858]=['',4.8188941,-0.6397773,5.34];
ybs[6859]=['',4.8330591,-1.0996776,6.14];
ybs[6860]=['κ Lyr',4.8024973,0.6296363,4.33];
ybs[6861]=['',4.8109157,0.095083,6.13];
ybs[6862]=['',4.8214718,-0.6322474,5.55];
ybs[6863]=['',4.8255702,-0.7696304,5.25];
ybs[6864]=['108 Her',4.8075897,0.5213396,5.63];
ybs[6865]=['107 Her',4.8079228,0.504081,5.12];
ybs[6866]=['',4.8183289,-0.1781231,6.33];
ybs[6867]=['ε Sgr',4.8243764,-0.5998894,1.85];
ybs[6868]=['',4.80168,0.8963802,6.3];
ybs[6869]=['',4.8191269,-0.2094696,5.73];
ybs[6870]=['',4.8131099,0.4066202,5.41];
ybs[6871]=['',4.8155172,0.2101638,5.89];
ybs[6872]=['ζ Sct',4.8209911,-0.1556995,4.68];
ybs[6873]=['',4.8162775,0.311356,5.25];
ybs[6874]=['',4.8069973,0.8680786,6.4];
ybs[6875]=['',4.8173358,0.2914858,6.22];
ybs[6876]=['18 Sgr',4.8279053,-0.53656,5.6];
ybs[6877]=['',4.829659,-0.6279245,6.15];
ybs[6878]=['',4.8225309,-0.0623063,6.38];
ybs[6879]=['',4.8089136,0.8575428,5.05];
ybs[6880]=['',4.8254715,-0.123246,6.31];
ybs[6881]=['',4.8319357,-0.5922032,6.3];
ybs[6882]=['',4.8372064,-0.839535,5.46];
ybs[6883]=['109 Her',4.8199557,0.3801835,3.84];
ybs[6884]=['21 Sgr',4.8288446,-0.358272,4.81];
ybs[6885]=['α Tel',4.8373582,-0.8020342,3.51];
ybs[6886]=['',4.8263453,-0.0273234,6.15];
ybs[6887]=['',4.8685858,-1.290616,5.89];
ybs[6888]=['',4.8269502,0.0889898,6.74];
ybs[6889]=['',4.8202069,0.6763577,6.36];
ybs[6890]=['',4.8290187,0.1404332,5.65];
ybs[6891]=['μ Lyr',4.8213513,0.6897653,5.12];
ybs[6892]=['',4.8252683,0.4783798,6.27];
ybs[6893]=['ζ Tel',4.8457221,-0.8561651,4.13];
ybs[6894]=['',4.8299518,0.2614692,6.37];
ybs[6895]=['',4.8400729,-0.5201226,5.92];
ybs[6896]=['',4.8514829,-1.0036727,5.76];
ybs[6897]=['',4.8395033,-0.4645931,6.31];
ybs[6898]=['',4.8433198,-0.6803227,5.64];
ybs[6899]=['',4.8183355,0.9305039,6.32];
ybs[6900]=['',4.9165841,-1.4273931,6.27];
ybs[6901]=['λ Sgr',4.8404935,-0.4434191,2.81];
ybs[6902]=['',4.8411379,-0.4667275,6.27];
ybs[6903]=['',4.8469765,-0.7649686,6.36];
ybs[6904]=['ν Pav',4.85843,-1.0866537,4.64];
ybs[6905]=['',4.8295342,0.5208637,5.83];
ybs[6906]=['59 Ser',4.8361317,0.0036872,5.21];
ybs[6907]=['',4.8400354,-0.3103965,6.2];
ybs[6908]=['φ Dra',4.801547,1.2452753,4.22];
ybs[6909]=['',4.8469242,-0.6777939,6.63];
ybs[6910]=['',4.8503556,-0.8238566,5.7];
ybs[6911]=['39 Dra',4.8181495,1.0264925,4.98];
ybs[6912]=['',4.8327611,0.4618881,6.53];
ybs[6913]=['',4.8387459,0.0656958,6.07];
ybs[6914]=['',4.8448524,-0.4636555,6.5];
ybs[6915]=['χ Dra',4.8022909,1.2696248,3.57];
ybs[6916]=['',4.8392708,0.1083801,5.73];
ybs[6917]=['',4.8465791,-0.4405212,6.59];
ybs[6918]=['γ Sct',4.8453854,-0.2539379,4.7];
ybs[6919]=['',4.854845,-0.7312277,6.04];
ybs[6920]=['',4.8479238,-0.2542088,5.96];
ybs[6921]=['',4.84992,-0.3265873,5.66];
ybs[6922]=['δ1 Tel',4.8582191,-0.8010581,4.96];
ybs[6923]=['60 Ser',4.8470135,-0.0343617,5.39];
ybs[6924]=['',4.8544379,-0.5754671,5.34];
ybs[6925]=['',4.8588242,-0.7590427,5.72];
ybs[6926]=['δ2 Tel',4.8594136,-0.7983018,5.07];
ybs[6927]=['',4.8671108,-1.0243407,6.44];
ybs[6928]=['',4.8495772,-0.0996123,6.28];
ybs[6929]=['',4.848537,0.071244,6.69];
ybs[6930]=['',4.8603991,-0.6926486,5.16];
ybs[6931]=['',4.8455939,0.4168281,5.9];
ybs[6932]=['',4.8553162,-0.3208839,5.14];
ybs[6933]=['42 Dra',4.8260196,1.1445485,4.82];
ybs[6934]=['',4.8549652,-0.188118,5.72];
ybs[6935]=['',4.85731,-0.3334847,6.68];
ybs[6936]=['',4.8632829,-0.6959294,6.22];
ybs[6937]=['',4.834614,1.0395941,6.43];
ybs[6938]=['',4.8505218,0.3635916,6.5];
ybs[6939]=['θ CrA',4.8655738,-0.7381664,4.64];
ybs[6940]=['κ1 CrA',4.8648376,-0.6754716,6.32];
ybs[6941]=['κ2 CrA',4.8648234,-0.6755734,5.65];
ybs[6942]=['',4.8708971,-0.922802,6.22];
ybs[6943]=['',4.8523443,0.2957602,5.77];
ybs[6944]=['',4.8591246,-0.2552753,6.37];
ybs[6945]=['61 Ser',4.8568665,-0.0171973,5.94];
ybs[6946]=['',4.8574174,0.0641849,6.43];
ybs[6947]=['',4.8607697,-0.2591358,5.5];
ybs[6948]=['',4.8670176,-0.5759197,5.28];
ybs[6949]=['24 Sgr',4.8662644,-0.4191177,5.49];
ybs[6950]=['',4.8648195,-0.2589185,5.76];
ybs[6951]=['',4.8632993,-0.1028502,6.36];
ybs[6952]=['',4.9622502,-1.4536292,7.16];
ybs[6953]=['25 Sgr',4.8691306,-0.4224276,6.51];
ybs[6954]=['',4.8594605,0.4125094,5.84];
ybs[6955]=['',4.8627902,0.1446321,6.42];
ybs[6956]=['',4.8593984,0.533587,5.48];
ybs[6957]=['',4.8725076,-0.363399,6.48];
ybs[6958]=['',4.8707301,-0.1912499,5.14];
ybs[6959]=['',4.8617959,0.5394972,6.59];
ybs[6960]=['',4.875709,-0.5179995,6.37];
ybs[6961]=['α Sct',4.8713513,-0.1435476,3.85];
ybs[6962]=['',4.8550934,0.9098965,6.56];
ybs[6963]=['',4.8663921,0.3575367,6.57];
ybs[6964]=['',4.868812,0.190431,6.4];
ybs[6965]=['',4.8703427,0.3180473,5.78];
ybs[6966]=['45 Dra',4.8562168,0.9959443,4.77];
ybs[6967]=['',4.8490195,1.1423731,6.59];
ybs[6968]=['',4.8714036,0.4123365,5.61];
ybs[6969]=['',4.8733464,0.296625,6.21];
ybs[6970]=['ζ Pav',4.9115774,-1.2462343,4.01];
ybs[6971]=['',4.8627301,0.9140679,5.36];
ybs[6972]=['',4.8696361,0.6017402,6.1];
ybs[6973]=['',4.8761751,0.1595691,5.39];
ybs[6974]=['',4.8910128,-0.8358019,5.86];
ybs[6975]=['',4.877091,0.1168007,5.45];
ybs[6976]=['',4.8836639,-0.3730954,5.94];
ybs[6977]=['',4.8840986,-0.2440609,6.47];
ybs[6978]=['',4.886396,-0.4098675,5.81];
ybs[6979]=['',4.8921551,-0.7533565,5.37];
ybs[6980]=['',4.8793438,0.1997041,6.42];
ybs[6981]=['',4.8814929,-0.0050383,5.75];
ybs[6982]=['',4.936011,-1.3585705,6.39];
ybs[6983]=['',4.8788912,0.2830719,6.29];
ybs[6984]=['',4.9068446,-1.127822,6.37];
ybs[6985]=['',4.8757786,0.584494,5.42];
ybs[6986]=['',4.8879389,-0.3670502,5.86];
ybs[6987]=['',4.8850704,-0.0553689,6.49];
ybs[6988]=['',4.8846575,-0.019062,6.66];
ybs[6989]=['α Lyr',4.8768605,0.6772561,0.03];
ybs[6990]=['',4.8844188,0.1545496,6.4];
ybs[6991]=['',4.875788,0.7547176,6.2];
ybs[6992]=['',4.9122533,-1.1262056,5.78];
ybs[6993]=['',4.900936,-0.8390143,6.49];
ybs[6994]=['',4.8373961,1.3537273,5.64];
ybs[6995]=['',4.892284,-0.1355854,5.84];
ybs[6996]=['',4.890075,0.0922582,6.38];
ybs[6997]=['',4.8819116,0.6927036,6.04];
ybs[6998]=['',4.8910644,0.1288105,6.28];
ybs[6999]=['26 Sgr',4.9010181,-0.4155666,6.23];
ybs[7000]=['',4.9202891,-1.1317783,4.79];
ybs[7001]=['',4.8707305,1.1433363,6.06];
ybs[7002]=['',4.8999655,-0.253791,6.42];
ybs[7003]=['',4.9184839,-1.065871,6.04];
ybs[7004]=['',4.8908129,0.5388084,6.36];
ybs[7005]=['',4.8881293,0.7148293,6.25];
ybs[7006]=['',4.8771501,1.0916538,5.74];
ybs[7007]=['',4.8911365,0.6700195,6.45];
ybs[7008]=['δ Sct',4.9022147,-0.1575892,4.72];
ybs[7009]=['λ CrA',4.910149,-0.6684509,5.13];
ybs[7010]=['',4.9187752,-0.9923333,6.22];
ybs[7011]=['',4.9054452,-0.3361542,6.35];
ybs[7012]=['',4.9035716,-0.1230433,6.15];
ybs[7013]=['',4.8046143,1.4518982,6.17];
ybs[7014]=['',4.9115657,-0.6404352,6.32];
ybs[7015]=['',4.9413525,-1.2735148,6.06];
ybs[7016]=['',4.8886539,0.9113749,6];
ybs[7017]=['',4.9123544,-0.6216421,4.87];
ybs[7018]=['',4.8979965,0.5522338,6.41];
ybs[7019]=['',4.915344,-0.6922248,5.43];
ybs[7020]=['ε Sct',4.9076315,-0.1440124,4.9];
ybs[7021]=['',4.8997783,0.6068472,6.47];
ybs[7022]=['',4.9090376,-0.1185857,6.31];
ybs[7023]=['',4.9140133,-0.4360948,5.83];
ybs[7024]=['θ Pav',4.9342011,-1.1353508,5.73];
ybs[7025]=['',4.9249877,-0.8738598,6.54];
ybs[7026]=['',4.915952,-0.3661077,6.36];
ybs[7027]=['φ Sgr',4.9177258,-0.4706397,3.17];
ybs[7028]=['4 Aql',4.9129582,0.0363841,5.02];
ybs[7029]=['',4.9045002,0.6863331,6.45];
ybs[7030]=['',4.891869,1.0955787,6.09];
ybs[7031]=['',4.9060768,0.6384513,6.01];
ybs[7032]=['',4.9074585,0.5576453,5.7];
ybs[7033]=['',4.918983,-0.3417537,6.42];
ybs[7034]=['28 Sgr',4.9205129,-0.3903726,5.37];
ybs[7035]=['',4.9113936,0.4121461,6.31];
ybs[7036]=['',4.9156288,0.0964244,5.83];
ybs[7037]=['46 Dra',4.9003049,0.9697526,5.04];
ybs[7038]=['μ CrA',4.927559,-0.7047602,5.24];
ybs[7039]=['ε1 Lyr',4.9091109,0.6927959,5.06];
ybs[7040]=['ε1 Lyr',4.9091037,0.6928153,6.02];
ybs[7041]=['ε2 Lyr',4.9092964,0.6918024,5.14];
ybs[7042]=['ε2 Lyr',4.9092964,0.6917976,5.37];
ybs[7043]=['',4.9216641,-0.1762662,5.71];
ybs[7044]=['ζ1 Lyr',4.9111291,0.6567591,4.36];
ybs[7045]=['ζ2 Lyr',4.9112607,0.6565751,5.73];
ybs[7046]=['',4.915501,0.3841471,6.51];
ybs[7047]=['5 Aql',4.9202416,-0.0163385,5.9];
ybs[7048]=['',4.9042106,0.9406574,6.11];
ybs[7049]=['110 Her',4.9158539,0.3590393,4.19];
ybs[7050]=['η1 CrA',4.9325707,-0.7618899,5.49];
ybs[7051]=['β Sct',4.9234376,-0.082412,4.22];
ybs[7052]=['',4.917385,0.4657844,4.83];
ybs[7053]=['',4.9354102,-0.7990644,5.81];
ybs[7054]=['',4.924819,-0.0991157,5.2];
ybs[7055]=['',4.9204191,0.3269253,6.17];
ybs[7056]=['η2 CrA',4.9357879,-0.7575876,5.61];
ybs[7057]=['111 Her',4.9218883,0.3177751,4.36];
ybs[7058]=['',4.9339542,-0.6060042,6.62];
ybs[7059]=['',4.9103842,0.9585555,6.23];
ybs[7060]=['',4.9308728,-0.3241878,6.47];
ybs[7061]=['',4.9171875,0.7237348,6.07];
ybs[7062]=['λ Pav',4.9492895,-1.0848731,4.22];
ybs[7063]=['',4.9068506,1.0659108,5.99];
ybs[7064]=['',4.9268975,0.0744864,6.21];
ybs[7065]=['',4.9345453,-0.3336193,6.75];
ybs[7066]=['29 Sgr',4.9349295,-0.354257,5.24];
ybs[7067]=['',4.9271224,0.4108611,6.15];
ybs[7068]=['',4.9201827,0.8087972,6.52];
ybs[7069]=['',4.9253545,0.5547213,6.06];
ybs[7070]=['',4.8996161,1.2359753,6.44];
ybs[7071]=['',4.9344256,-0.1027216,5.99];
ybs[7072]=['',4.9184273,0.9252609,5.88];
ybs[7073]=['',4.9338899,0.015063,6.25];
ybs[7074]=['',4.9299939,0.3378153,5.88];
ybs[7075]=['κ Tel',4.9499247,-0.9089412,5.17];
ybs[7076]=['30 Sgr',4.9401268,-0.3863164,6.61];
ybs[7077]=['',4.937337,-0.1375299,6.8];
ybs[7078]=['',4.922938,0.8569745,6.4];
ybs[7079]=['',4.9312724,0.4376124,6.59];
ybs[7080]=['',4.9485325,-0.8127373,5.54];
ybs[7081]=['',4.9522693,-0.9058577,6.31];
ybs[7082]=['',4.9402077,-0.1701034,5.83];
ybs[7083]=['',4.9512487,-0.8435318,6.19];
ybs[7084]=['',4.9255912,0.8516132,6.12];
ybs[7085]=['',4.9508939,-0.8125674,6.19];
ybs[7086]=['',4.9331082,0.5525078,6.64];
ybs[7087]=['',4.938489,0.1920593,6.55];
ybs[7088]=['ν1 Lyr',4.9331848,0.5731659,5.91];
ybs[7089]=['8 Aql',4.9416803,-0.0574147,6.1];
ybs[7090]=['ν2 Lyr',4.9337078,0.5685952,5.25];
ybs[7091]=['',4.9474533,-0.4646324,6.29];
ybs[7092]=['',4.9482012,-0.5122638,6.13];
ybs[7093]=['',4.94861,-0.5359074,6.63];
ybs[7094]=['β Lyr',4.9345305,0.5827681,3.45];
ybs[7095]=['κ Pav',4.970951,-1.1728995,4.44];
ybs[7096]=['',4.9579067,-0.8700226,6.6];
ybs[7097]=['',4.9439217,0.2442418,6.14];
ybs[7098]=['',4.9491733,-0.1666276,6.34];
ybs[7099]=['',4.9698608,-1.0955392,6.48];
ybs[7100]=['',4.94139,0.502861,6.18];
ybs[7101]=['112 Her',4.9446702,0.3744404,5.48];
ybs[7102]=['33 Sgr',4.9538776,-0.3722809,5.69];
ybs[7103]=['',4.9410252,0.6382109,6.09];
ybs[7104]=['ν1 Sgr',4.9546717,-0.3964569,4.83];
ybs[7105]=['',4.9096431,1.293468,5.27];
ybs[7106]=['',4.9429529,0.7227717,6.28];
ybs[7107]=['',4.9567666,-0.2718017,5.1];
ybs[7108]=['ν2 Sgr',4.9588055,-0.3951634,4.99];
ybs[7109]=['σ Sgr',4.9596073,-0.4584349,2.02];
ybs[7110]=['',4.9649507,-0.7449009,5.36];
ybs[7111]=['',4.9396692,0.9250785,5.51];
ybs[7112]=['50 Dra',4.911496,1.3170056,5.35];
ybs[7113]=['ο Dra',4.9372255,1.0370079,4.66];
ybs[7114]=['',4.9602805,-0.2852964,5.79];
ybs[7115]=['ω Pav',4.9768273,-1.0501357,5.14];
ybs[7116]=['',4.9627171,-0.4039205,5.93];
ybs[7117]=['',4.9663271,-0.6512219,5.38];
ybs[7118]=['',4.9843564,-1.1627433,6.01];
ybs[7119]=['δ1 Lyr',4.9502466,0.6457888,5.58];
ybs[7120]=['',4.9528742,0.4876283,5.62];
ybs[7121]=['113 Her',4.955419,0.3957515,4.59];
ybs[7122]=['λ Tel',4.9753045,-0.9233932,4.87];
ybs[7123]=['',4.9591567,0.1159874,5.57];
ybs[7124]=['',4.9704169,-0.6944976,6.31];
ybs[7125]=['',4.9470651,0.8855332,4.92];
ybs[7126]=['',4.9522176,0.7200376,7.3];
ybs[7127]=['δ2 Lyr',4.9536409,0.6445258,4.3];
ybs[7128]=['',4.955427,0.5933864,6.02];
ybs[7129]=['θ1 Ser',4.9625704,0.073903,4.62];
ybs[7130]=['θ2 Ser',4.9626722,0.0738741,4.98];
ybs[7131]=['',4.9634837,-0.0308781,6.22];
ybs[7132]=['',4.9635365,0.0436671,6.15];
ybs[7133]=['ξ1 Sgr',4.9684182,-0.3599745,5.08];
ybs[7134]=['',4.9549501,0.7266273,5.44];
ybs[7135]=['',4.9613612,0.3146061,6.63];
ybs[7136]=['',4.961524,0.3165312,5.69];
ybs[7137]=['η Sct',4.9666155,-0.1014896,4.83];
ybs[7138]=['ξ2 Sgr',4.9701313,-0.3678297,3.51];
ybs[7139]=['',4.9733137,-0.541125,6.12];
ybs[7140]=['ε CrA',4.9752472,-0.647087,4.87];
ybs[7141]=['',4.9487066,1.0038466,6.22];
ybs[7142]=['',4.9540334,0.8532834,5.77];
ybs[7143]=['',4.9729638,-0.4336228,6.62];
ybs[7144]=['',4.9774061,-0.6894456,6.49];
ybs[7145]=['13 Lyr',4.9568188,0.7675302,4.04];
ybs[7146]=['64 Ser',4.967243,0.0447948,5.57];
ybs[7147]=['',4.9731587,-0.3926556,6.14];
ybs[7148]=['',4.9044344,1.3956847,6.39];
ybs[7149]=['',4.999866,-1.1993971,5.88];
ybs[7150]=['',4.9648755,0.5747799,5.22];
ybs[7151]=['',4.9719924,0.1094695,6.21];
ybs[7152]=['',4.9775143,-0.3234875,6.37];
ybs[7153]=['',4.9709006,0.303558,5.38];
ybs[7154]=['',4.9770682,-0.2235439,5.53];
ybs[7155]=['10 Aql',4.9733802,0.2432765,5.89];
ybs[7156]=['',4.9820038,-0.4347479,6.36];
ybs[7157]=['',4.9854066,-0.6462511,6.69];
ybs[7158]=['',4.9854866,-0.6462703,6.4];
ybs[7159]=['',4.9730119,0.3460322,6.5];
ybs[7160]=['11 Aql',4.9747587,0.2383198,5.23];
ybs[7161]=['',4.975751,0.1775551,6.75];
ybs[7162]=['',4.9689619,0.6684208,5.89];
ybs[7163]=['48 Dra',4.9616702,1.0095997,5.66];
ybs[7164]=['ε Aql',4.9770074,0.2635591,4.02];
ybs[7165]=['',4.9904671,-0.7308746,6.23];
ybs[7166]=['γ Lyr',4.9732525,0.5710984,3.24];
ybs[7167]=['',4.9720644,0.7105434,6.22];
ybs[7168]=['υ Dra',4.948495,1.244884,4.82];
ybs[7169]=['',4.9771246,0.4583774,5.27];
ybs[7170]=['',4.9872124,-0.3955251,6.24];
ybs[7171]=['',4.9781933,0.3987619,6.29];
ybs[7172]=['',4.9647749,1.016767,6.46];
ybs[7173]=['',4.974026,0.6850413,6.41];
ybs[7174]=['',4.9865945,-0.2661438,6.32];
ybs[7175]=['',4.9590245,1.1395012,5.63];
ybs[7176]=['ζ CrA',4.9947095,-0.7340994,4.75];
ybs[7177]=['',4.9990769,-0.8898325,5.93];
ybs[7178]=['',4.963342,1.0895692,6.45];
ybs[7179]=['λ Lyr',4.9779583,0.5616157,4.93];
ybs[7180]=['12 Aql',4.9867616,-0.0995757,4.02];
ybs[7181]=['ζ Sgr',4.9918233,-0.5209128,2.6];
ybs[7182]=['',4.9909346,-0.4330662,5.65];
ybs[7183]=['',4.9722339,0.8873512,6.3];
ybs[7184]=['',4.995253,-0.6670434,5.74];
ybs[7185]=['',4.9832442,0.3375986,6.39];
ybs[7186]=['',4.9426815,1.323244,6.22];
ybs[7187]=['',4.9844244,0.3641981,6.69];
ybs[7188]=['',4.9787869,0.7106449,6.65];
ybs[7189]=['',4.9838069,0.4594532,5.69];
ybs[7190]=['',4.9933197,-0.3352934,6.05];
ybs[7191]=['',4.9818219,0.5905382,6.01];
ybs[7192]=['',4.9935465,-0.3328155,6.37];
ybs[7193]=['',4.9851372,0.4373679,6.72];
ybs[7194]=['',4.9863202,0.3891652,6.4];
ybs[7195]=['',4.9892071,0.1467493,6.3];
ybs[7196]=['14 Aql',4.9920454,-0.0639598,5.42];
ybs[7197]=['',4.9776516,0.8825485,5.38];
ybs[7198]=['',4.9997611,-0.5412583,5.5];
ybs[7199]=['',4.9856946,0.5873902,6.39];
ybs[7200]=['ρ Tel',5.0095573,-0.9128878,5.16];
ybs[7201]=['',4.9945844,0.0323492,5.83];
ybs[7202]=['16 Lyr',4.9832473,0.8197474,5.01];
ybs[7203]=['',4.9910198,0.3437479,6.09];
ybs[7204]=['ο Sgr',5.0004909,-0.3788487,3.77];
ybs[7205]=['49 Dra',4.9792869,0.9719958,5.48];
ybs[7206]=['',4.9973287,0.0587386,6.73];
ybs[7207]=['',4.9986415,-0.0986102,6.9];
ybs[7208]=['',5.0275825,-1.1935697,5.33];
ybs[7209]=['',4.9945822,0.3717971,6.52];
ybs[7210]=['',5.0117892,-0.8423424,5.97];
ybs[7211]=['',4.9686485,1.2141016,6.52];
ybs[7212]=['15 Aql',5.0010094,-0.0697442,5.42];
ybs[7213]=['γ CrA',5.0088051,-0.6462455,4.93];
ybs[7214]=['γ CrA',5.0088051,-0.6462455,4.99];
ybs[7215]=['σ Oct',5.615683,-1.5506617,5.47];
ybs[7216]=['',4.9857334,0.9127163,6.31];
ybs[7217]=['',5.0046195,-0.2726994,5.97];
ybs[7218]=['',5.002444,-0.025783,6.53];
ybs[7219]=['',5.01084,-0.6592779,6.16];
ybs[7220]=['',5.0209732,-0.9718467,6.49];
ybs[7221]=['τ Sgr',5.0105975,-0.4823064,3.32];
ybs[7222]=['ζ Aql',5.0023076,0.2425811,2.99];
ybs[7223]=['λ Aql',5.0066618,-0.084587,3.44];
ybs[7224]=['',4.9995852,0.5546605,5.56];
ybs[7225]=['',4.999665,0.5370134,6.06];
ybs[7226]=['',5.0098117,-0.2826127,6.03];
ybs[7227]=['',5.0131509,-0.4991678,6.04];
ybs[7228]=['',5.0110902,-0.3263985,6.29];
ybs[7229]=['δ CrA',5.0174281,-0.7061507,4.59];
ybs[7230]=['',5.0067109,0.1442699,6.09];
ybs[7231]=['',5.003255,0.5228546,6.31];
ybs[7232]=['',5.010401,0.0118359,6.56];
ybs[7233]=['',5.0161461,-0.4296972,6.3];
ybs[7234]=['',4.9610676,1.3453332,6.54];
ybs[7235]=['18 Aql',5.0092515,0.1938669,5.09];
ybs[7236]=['',5.0160712,-0.3360262,5.54];
ybs[7237]=['',5.0072488,0.423888,5.77];
ybs[7238]=['51 Dra',4.9978606,0.932561,5.38];
ybs[7239]=['',4.9992619,0.8719425,6.43];
ybs[7240]=['',5.0069981,0.5002942,5.55];
ybs[7241]=['α CrA',5.0221698,-0.6608979,4.11];
ybs[7242]=['',5.0231138,-0.6944693,6.46];
ybs[7243]=['',5.0226594,-0.6305328,6.56];
ybs[7244]=['',5.0245605,-0.7304978,5.88];
ybs[7245]=['',5.0047877,0.7234356,6.49];
ybs[7246]=['β CrA',5.024683,-0.6859624,4.11];
ybs[7247]=['',5.0132966,0.2947896,6.07];
ybs[7248]=['17 Lyr',5.0102918,0.5678989,5.23];
ybs[7249]=['ι Lyr',5.0095534,0.6307052,5.28];
ybs[7250]=['',5.0135564,0.3793612,6.23];
ybs[7251]=['π Sgr',5.0226207,-0.3662699,2.89];
ybs[7252]=['',5.0227375,-0.3449766,6.13];
ybs[7253]=['19 Aql',5.0182601,0.106653,5.22];
ybs[7254]=['',5.0164171,0.2947621,6.48];
ybs[7255]=['',5.0290309,-0.6800919,6.36];
ybs[7256]=['',5.0222573,-0.0068096,6.34];
ybs[7257]=['',5.0297667,-0.5142351,6.3];
ybs[7258]=['',5.0374459,-0.8804633,6.13];
ybs[7259]=['',5.0173696,0.6045463,6.74];
ybs[7260]=['',5.0338896,-0.6552592,6.57];
ybs[7261]=['τ Pav',5.0565252,-1.2068762,6.27];
ybs[7262]=['',5.0132863,0.9156441,5.81];
ybs[7263]=['',5.0344397,-0.3773185,6.41];
ybs[7264]=['',5.0379459,-0.4514636,5.8];
ybs[7265]=['',5.0590384,-1.162728,5.53];
ybs[7266]=['20 Aql',5.0348285,-0.1378821,5.34];
ybs[7267]=['',5.0284154,0.4673036,6.36];
ybs[7268]=['',5.045275,-0.7880654,5.92];
ybs[7269]=['',5.0375212,-0.2136774,5.51];
ybs[7270]=['19 Lyr',5.0292978,0.5466746,5.98];
ybs[7271]=['',5.0271233,0.7062955,6.18];
ybs[7272]=['',5.0334522,0.2947104,6.73];
ybs[7273]=['',5.033428,0.3768815,5.93];
ybs[7274]=['21 Aql',5.0389648,0.0407273,5.15];
ybs[7275]=['',5.0389495,0.0969658,6.49];
ybs[7276]=['',5.0527071,-0.7928158,5.4];
ybs[7277]=['55 Dra',5.0171368,1.1521996,6.25];
ybs[7278]=['',5.0480112,-0.4212922,6.25];
ybs[7279]=['ψ Sgr',5.0480001,-0.4400982,4.85];
ybs[7280]=['',5.0294489,0.8707984,6.75];
ybs[7281]=['',5.0294779,0.8708323,6.57];
ybs[7282]=['53 Dra',5.0269745,0.9930543,5.12];
ybs[7283]=['',5.0527841,-0.5843402,6.25];
ybs[7284]=['',5.0612105,-0.9310339,6.38];
ybs[7285]=['η Lyr',5.0375636,0.6839235,4.39];
ybs[7286]=['',5.044097,0.3533223,6];
ybs[7287]=['',5.0455736,0.2639691,5.57];
ybs[7288]=['1 Sge',5.0451311,0.3712819,5.64];
ybs[7289]=['',5.0452532,0.5334964,5.85];
ybs[7290]=['22 Aql',5.0511123,0.0851032,5.59];
ybs[7291]=['43 Sgr',5.0568739,-0.3300608,4.96];
ybs[7292]=['',5.0477408,0.4799052,6.54];
ybs[7293]=['1 Vul',5.0491689,0.3740489,4.77];
ybs[7294]=['',5.0504459,0.2545737,5.63];
ybs[7295]=['',5.0478937,0.4881329,6.16];
ybs[7296]=['54 Dra',5.0366753,1.007837,4.99];
ybs[7297]=['δ Dra',5.0289406,1.1815995,3.07];
ybs[7298]=['',5.0436006,0.8746092,6.27];
ybs[7299]=['59 Dra',5.0104496,1.3368797,5.13];
ybs[7300]=['',5.056828,0.0361924,6.19];
ybs[7301]=['θ Lyr',5.0490183,0.6662757,4.36];
ybs[7302]=['ω1 Aql',5.0565365,0.2031087,5.28];
ybs[7303]=['',5.0664912,-0.6174685,5.59];
ybs[7304]=['',5.0627012,-0.2704166,6.06];
ybs[7305]=['2 Vul',5.0556907,0.4026034,5.43];
ybs[7306]=['23 Aql',5.0600837,0.0196815,5.1];
ybs[7307]=['',5.0894632,-1.1925064,6.34];
ybs[7308]=['24 Aql',5.0614485,0.0066622,6.41];
ybs[7309]=['',5.0505569,0.8210126,6];
ybs[7310]=['',5.0696631,-0.5545669,6.58];
ybs[7311]=['',5.0565786,0.5421737,6.68];
ybs[7312]=['',5.0612431,0.1686092,6.32];
ybs[7313]=['',5.0605533,0.3430103,6.58];
ybs[7314]=['',5.0701014,-0.3902379,5.58];
ybs[7315]=['κ Cyg',5.0510708,0.932182,3.77];
ybs[7316]=['η Tel',5.0817815,-0.9490899,5.05];
ybs[7317]=['',5.0744555,-0.6098162,6.48];
ybs[7318]=['28 Aql',5.0645298,0.216729,5.53];
ybs[7319]=['ω2 Aql',5.0655577,0.2020751,6.02];
ybs[7320]=['26 Aql',5.0690728,-0.0937661,5.01];
ybs[7321]=['',5.0777519,-0.7325454,6.34];
ybs[7322]=['',5.0610434,0.5834891,6.6];
ybs[7323]=['27 Aql',5.0691139,-0.0148139,5.49];
ybs[7324]=['β1 Sgr',5.0800016,-0.7751755,4.01];
ybs[7325]=['',5.0606346,0.6542859,6.22];
ybs[7326]=['',5.0742611,-0.3349314,6.26];
ybs[7327]=['ρ1 Sgr',5.074446,-0.3107243,3.93];
ybs[7328]=['',5.0580843,0.865893,6.31];
ybs[7329]=['υ Sgr',5.0746048,-0.2776983,4.61];
ybs[7330]=['β2 Sgr',5.0825616,-0.781119,4.29];
ybs[7331]=['ρ2 Sgr',5.0752276,-0.3187706,5.87];
ybs[7332]=['',5.0633843,0.6522892,6.31];
ybs[7333]=['',5.0674462,0.6148699,6.31];
ybs[7334]=['',5.0770124,-0.1423622,6.31];
ybs[7335]=['α Sgr',5.0851924,-0.7080955,3.97];
ybs[7336]=['',5.076785,-0.003633,5.83];
ybs[7337]=['',5.0874428,-0.7623175,6.17];
ybs[7338]=['',5.0619077,0.9497887,6.26];
ybs[7339]=['τ Dra',5.0401394,1.2810018,4.45];
ybs[7340]=['',5.0801755,-0.1283785,6.32];
ybs[7341]=['',5.078375,0.1737931,6.35];
ybs[7342]=['',5.0872047,-0.4855515,6.04];
ybs[7343]=['',5.0644143,1.0068519,5.91];
ybs[7344]=['',5.0796278,0.2612029,6.64];
ybs[7345]=['3 Vul',5.0779055,0.4591444,5.18];
ybs[7346]=['',5.0762868,0.5857748,6.06];
ybs[7347]=['',5.0897426,-0.5107422,5.93];
ybs[7348]=['',5.0611845,1.1245784,6.52];
ybs[7349]=['χ1 Sgr',5.0904296,-0.4269551,5.03];
ybs[7350]=['χ3 Sgr',5.0913658,-0.4174169,5.43];
ybs[7351]=['',5.0822331,0.3544672,6.4];
ybs[7352]=['',5.0694446,1.0089857,6.43];
ybs[7353]=['',5.0885926,-0.0844469,6.52];
ybs[7354]=['',5.0903764,-0.2417462,5.69];
ybs[7355]=['',5.0807111,0.5806206,6.37];
ybs[7356]=['2 Sge',5.0849255,0.2964113,6.25];
ybs[7357]=['',5.1033454,-0.9473298,5.69];
ybs[7358]=['π Dra',5.0648503,1.1476925,4.59];
ybs[7359]=['2 Cyg',5.0833321,0.5177794,4.97];
ybs[7360]=['31 Aql',5.0877407,0.2092667,5.16];
ybs[7361]=['',5.0844801,0.491015,6.53];
ybs[7362]=['50 Sgr',5.094873,-0.3792645,5.59];
ybs[7363]=['',5.0828758,0.6369944,6.36];
ybs[7364]=['δ Aql',5.090366,0.055164,3.36];
ybs[7365]=['',5.0940191,-0.2619168,5.72];
ybs[7366]=['',5.0949816,-0.2531543,6.7];
ybs[7367]=['',5.0979332,-0.5183032,5.67];
ybs[7368]=['',5.0788482,0.878182,6.51];
ybs[7369]=['',5.0817555,0.7580505,5.84];
ybs[7370]=['',5.1203579,-1.1935395,5.96];
ybs[7371]=['',5.0891797,0.3546026,6.31];
ybs[7372]=['4 Vul',5.0896498,0.346352,5.16];
ybs[7373]=['',5.0892335,0.4356104,6.19];
ybs[7374]=['ν Aql',5.0949153,0.0067209,4.66];
ybs[7375]=['',5.1125072,-0.9667912,6.13];
ybs[7376]=['',5.0939485,0.2281192,5.74];
ybs[7377]=['5 Vul',5.0928821,0.35158,5.63];
ybs[7378]=['',5.0940179,0.3479801,5.81];
ybs[7379]=['',5.1094054,-0.7574344,5.71];
ybs[7380]=['μ Tel',5.1155147,-0.9610013,6.3];
ybs[7381]=['λ UMi',4.405267,1.5532467,6.38];
ybs[7382]=['4 Cyg',5.0918332,0.634671,5.15];
ybs[7383]=['',5.0989783,0.2500964,6.32];
ybs[7384]=['',5.1028013,0.05197,5.85];
ybs[7385]=['',5.1105789,-0.4701452,5.52];
ybs[7386]=['',5.1124521,-0.5592697,6.6];
ybs[7387]=['35 Aql',5.1057596,0.0348718,5.8];
ybs[7388]=['',5.0884376,1.0135672,6.6];
ybs[7389]=['',5.1075688,-0.1221028,6.61];
ybs[7390]=['',5.0981024,0.6630163,6.34];
ybs[7391]=['',5.1070559,0.0051309,6.25];
ybs[7392]=['α Vul',5.1035394,0.4313149,4.44];
ybs[7393]=['8 Vul',5.1046044,0.4331253,5.81];
ybs[7394]=['',5.1068366,0.2555809,5.56];
ybs[7395]=['ι1 Cyg',5.0962851,0.9139824,5.75];
ybs[7396]=['7 Vul',5.1065288,0.3547831,6.33];
ybs[7397]=['',5.1148316,-0.371118,6.13];
ybs[7398]=['',5.1254363,-0.9273981,5.75];
ybs[7399]=['',5.1107801,0.0515305,6.09];
ybs[7400]=['',5.0906765,1.0926358,6.38];
ybs[7401]=['36 Aql',5.1131114,-0.0478276,5.03];
ybs[7402]=['',5.1124047,0.0609634,6.05];
ybs[7403]=['',5.1268022,-0.7892713,5.61];
ybs[7404]=['β1 Cyg',5.1121977,0.4888362,3.08];
ybs[7405]=['β2 Cyg',5.112343,0.4889335,5.11];
ybs[7406]=['',5.1120525,0.6331556,6.25];
ybs[7407]=['ι2 Cyg',5.1062701,0.9036905,3.79];
ybs[7408]=['',5.1150585,0.4654109,5.87];
ybs[7409]=['',5.1298712,-0.697858,5.89];
ybs[7410]=['',5.06257,1.390087,6.05];
ybs[7411]=['ι Tel',5.1351055,-0.8385991,4.9];
ybs[7412]=['',5.0270792,1.4573907,6.53];
ybs[7413]=['8 Cyg',5.1164687,0.6021752,4.74];
ybs[7414]=['',5.1134582,0.8788679,5.53];
ybs[7415]=['',5.1125128,0.9735554,6.37];
ybs[7416]=['μ Aql',5.1277056,0.1296632,4.45];
ybs[7417]=['37 Aql',5.1328388,-0.1834246,5.12];
ybs[7418]=['51 Sgr',5.1373404,-0.4305354,5.65];
ybs[7419]=['',5.1343641,-0.1293163,6.34];
ybs[7420]=['',5.1348125,-0.2129604,6.27];
ybs[7421]=['',5.1500795,-1.0110816,6.18];
ybs[7422]=['',5.1577612,-1.1629502,6.39];
ybs[7423]=['',5.1242456,0.6774002,6.61];
ybs[7424]=['9 Vul',5.1293909,0.345991,5];
ybs[7425]=['',5.1336665,0.0517364,6.38];
ybs[7426]=['',5.1388773,-0.3281443,6.11];
ybs[7427]=['52 Sgr',5.1403059,-0.4333996,4.6];
ybs[7428]=['9 Cyg',5.1301436,0.5151104,5.38];
ybs[7429]=['',5.1239051,0.8606646,5.96];
ybs[7430]=['',5.1415577,-0.3172888,5.64];
ybs[7431]=['',5.1287397,0.7411236,5.35];
ybs[7432]=['',5.1364875,0.1954992,6.68];
ybs[7433]=['κ Aql',5.1404349,-0.1217507,4.95];
ybs[7434]=['ι Aql',5.1394925,-0.0215511,4.36];
ybs[7435]=['',5.1204414,1.0508318,6.29];
ybs[7436]=['',5.1369361,0.252078,6.38];
ybs[7437]=['',5.1086367,1.2398441,6.07];
ybs[7438]=['',5.1265168,0.8951257,5.73];
ybs[7439]=['',5.1360725,0.3950918,6.32];
ybs[7440]=['',5.1282373,0.8415134,6.67];
ybs[7441]=['',5.1436707,-0.2487028,5.47];
ybs[7442]=['',5.1652333,-1.1484247,6.09];
ybs[7443]=['',5.1397197,0.1976583,5.98];
ybs[7444]=['11 Cyg',5.1339462,0.6456933,6.05];
ybs[7445]=['',5.1383018,0.3557729,7.14];
ybs[7446]=['',5.1578437,-0.9488302,6.26];
ybs[7447]=['42 Aql',5.1442635,-0.0802042,5.46];
ybs[7448]=['',5.1544409,-0.7893228,6.25];
ybs[7449]=['σ Dra',5.1149842,1.2166727,4.68];
ybs[7450]=['ε Sge',5.1413393,0.2882346,5.66];
ybs[7451]=['',5.155066,-0.6873116,6.61];
ybs[7452]=['',5.1336,0.8777205,6.52];
ybs[7453]=['',5.1402266,0.5128713,6.43];
ybs[7454]=['',5.1388364,0.6708261,6.5];
ybs[7455]=['',5.1370819,0.7809727,5.17];
ybs[7456]=['θ Cyg',5.1358292,0.8774194,4.48];
ybs[7457]=['53 Sgr',5.1538387,-0.4079636,6.34];
ybs[7458]=['',5.1484716,0.0599398,6.35];
ybs[7459]=['',5.1455442,0.3636413,6.48];
ybs[7460]=['',5.1551255,-0.4079756,5.97];
ybs[7461]=['σ Aql',5.150045,0.0951308,5.17];
ybs[7462]=['',5.1506461,0.2901487,6.38];
ybs[7463]=['54 Sgr',5.1574857,-0.2834366,6.2];
ybs[7464]=['',5.1424739,0.8610847,6.47];
ybs[7465]=['φ Cyg',5.1498805,0.5271973,4.69];
ybs[7466]=['α Sge',5.1535303,0.3153308,4.37];
ybs[7467]=['45 Aql',5.1569224,-0.0099052,5.67];
ybs[7468]=['',5.1513259,0.5939738,6.1];
ybs[7469]=['',5.1550715,0.3583176,6.5];
ybs[7470]=['14 Cyg',5.1494784,0.7482432,5.4];
ybs[7471]=['',5.1451867,0.9603902,5.82];
ybs[7472]=['',5.1557706,0.4148823,6.64];
ybs[7473]=['',5.1580297,0.2420648,6.01];
ybs[7474]=['',5.1498403,0.8030426,6.2];
ybs[7475]=['β Sge',5.1577044,0.3059531,4.37];
ybs[7476]=['55 Sgr',5.1653094,-0.2804637,5.06];
ybs[7477]=['',5.1583829,0.3928138,6.36];
ybs[7478]=['',5.1710845,-0.6542156,6.16];
ybs[7479]=['',5.1548986,0.7527821,6.16];
ybs[7480]=['46 Aql',5.1629838,0.2137614,6.34];
ybs[7481]=['',5.2361369,-1.4187437,6.39];
ybs[7482]=['',5.1553886,0.7954953,5.06];
ybs[7483]=['',5.1698218,-0.2690423,5.49];
ybs[7484]=['χ Aql',5.164539,0.2073649,5.27];
ybs[7485]=['',5.2009939,-1.2644065,5.41];
ybs[7486]=['',5.1606267,0.7035071,6.23];
ybs[7487]=['',5.1511943,1.0569775,6.51];
ybs[7488]=['',5.1649433,0.5128863,6.49];
ybs[7489]=['',5.1644712,0.5669034,5.94];
ybs[7490]=['16 Cyg',5.1592748,0.8827747,5.96];
ybs[7491]=['',5.1595012,0.8826394,6.2];
ybs[7492]=['',5.1663658,0.5363978,6.05];
ybs[7493]=['10 Vul',5.1690132,0.4507653,5.49];
ybs[7494]=['',5.1812352,-0.5559285,5.52];
ybs[7495]=['',5.1699007,0.4745666,6.28];
ybs[7496]=['',5.1599096,0.9689617,6.48];
ybs[7497]=['ν Tel',5.1916654,-0.9827108,5.35];
ybs[7498]=['ψ Aql',5.1732216,0.233145,6.26];
ybs[7499]=['',5.1692468,0.5972089,6.05];
ybs[7500]=['',5.2014099,-1.1650852,6.45];
ybs[7501]=['',5.1683591,0.7300369,5.84];
ybs[7502]=['56 Sgr',5.1822113,-0.3439123,4.86];
ybs[7503]=['',5.1794632,-0.0493443,6.48];
ybs[7504]=['15 Cyg',5.1709143,0.6529222,4.89];
ybs[7505]=['',5.1736484,0.511735,6.82];
ybs[7506]=['υ Aql',5.1782089,0.133855,5.91];
ybs[7507]=['',5.1726398,0.6016032,6.57];
ybs[7508]=['',5.1952282,-0.922063,6.25];
ybs[7509]=['',5.1647511,1.0135307,6.22];
ybs[7510]=['',5.1730813,0.7116086,6.34];
ybs[7511]=['',5.2060939,-1.1439962,6.05];
ybs[7512]=['γ Aql',5.1807004,0.1862199,2.72];
ybs[7513]=['',5.1666882,0.9965368,6.27];
ybs[7514]=['',5.2024567,-1.0647015,6.21];
ybs[7515]=['δ Cyg',5.1734902,0.7886513,2.87];
ybs[7516]=['',5.1770219,0.6308848,6.43];
ybs[7517]=['',5.1779363,0.6120662,6.09];
ybs[7518]=['',5.2038801,-1.0320899,5.42];
ybs[7519]=['',5.189351,-0.23817,6.11];
ybs[7520]=['',5.1818887,0.4396542,6.62];
ybs[7521]=['17 Cyg',5.1804894,0.5896436,4.99];
ybs[7522]=['',5.1812145,0.5749988,6.18];
ybs[7523]=['δ Sge',5.1853364,0.3244738,3.82];
ybs[7524]=['',5.2005031,-0.8290117,5.94];
ybs[7525]=['',5.1949311,-0.5014523,6.05];
ybs[7526]=['',5.1869056,0.4440273,5.95];
ybs[7527]=['',5.1935546,-0.1887252,6.04];
ybs[7528]=['',5.1905023,0.1876498,6.44];
ybs[7529]=['',5.1847715,0.6713287,5.77];
ybs[7530]=['π Aql',5.1913138,0.2072282,5.72];
ybs[7531]=['',5.1673369,1.2111189,5.92];
ybs[7532]=['ζ Sge',5.1922534,0.3350999,5];
ybs[7533]=['',5.1840856,0.8371391,6.12];
ybs[7534]=['',5.2115883,-0.9583878,5.74];
ybs[7535]=['',5.2116905,-0.9584845,6.5];
ybs[7536]=['',5.1904744,0.6173023,6.53];
ybs[7537]=['',5.1910553,0.5845929,6.44];
ybs[7538]=['',5.2070379,-0.6949088,5.33];
ybs[7539]=['51 Aql',5.2011561,-0.1868392,5.39];
ybs[7540]=['',5.1983893,0.1389412,6.51];
ybs[7541]=['',5.1934665,0.6766253,6.11];
ybs[7542]=['',5.1959428,0.4973939,6.38];
ybs[7543]=['α Aql',5.2005016,0.1558023,0.77];
ybs[7544]=['',5.2212501,-1.0665754,6.24];
ybs[7545]=['',5.2026388,-0.0419251,6.13];
ybs[7546]=['ο Aql',5.2015097,0.1828084,5.11];
ybs[7547]=['57 Sgr',5.2076495,-0.3313644,5.92];
ybs[7548]=['',5.202708,0.1691051,6.25];
ybs[7549]=['',5.1782333,1.1954562,6.34];
ybs[7550]=['χ Cyg',5.1985978,0.5754785,4.23];
ybs[7551]=['12 Vul',5.2012418,0.3956415,4.95];
ybs[7552]=['19 Cyg',5.1983058,0.6768528,5.12];
ybs[7553]=['',5.1984388,0.7096168,5.69];
ybs[7554]=['',5.1993011,0.6612145,6.06];
ybs[7555]=['',5.2059697,0.2039937,6.13];
ybs[7556]=['η Aql',5.2081434,0.0185855,3.9];
ybs[7557]=['',5.2114414,-0.2538304,6.48];
ybs[7558]=['',5.2068944,0.1816989,6.54];
ybs[7559]=['',5.2053292,0.4372272,5.57];
ybs[7560]=['9 Sge',5.2070461,0.3269206,6.23];
ybs[7561]=['',5.2119397,-0.0533149,5.65];
ybs[7562]=['20 Cyg',5.1975581,0.9258332,5.03];
ybs[7563]=['',5.2010591,0.827912,6.2];
ybs[7564]=['',5.2169766,-0.4167999,6.18];
ybs[7565]=['',5.240265,-1.2060458,5.75];
ybs[7566]=['',5.2119708,0.0778469,6.53];
ybs[7567]=['ι Sgr',5.2220564,-0.7296801,4.13];
ybs[7568]=['ε Dra',5.1839781,1.2273977,3.83];
ybs[7569]=['',5.2058741,0.6368944,6.1];
ybs[7570]=['56 Aql',5.215724,-0.148598,5.79];
ybs[7571]=['',5.2208168,-0.57571,6.46];
ybs[7572]=['',5.2307958,-1.0099161,6.53];
ybs[7573]=['',5.2315308,-1.0269463,5.26];
ybs[7574]=['',5.2410109,-1.1990339,6.39];
ybs[7575]=['',5.2039748,0.8218138,5.62];
ybs[7576]=['ε Pav',5.2497183,-1.2714208,3.96];
ybs[7577]=['',5.2044978,0.8376004,5.91];
ybs[7578]=['13 Vul',5.2116307,0.4213131,4.58];
ybs[7579]=['57 Aql',5.2178499,-0.1425386,5.71];
ybs[7580]=['57 Aql',5.2178865,-0.142713,6.49];
ybs[7581]=['ξ Aql',5.2156368,0.1487289,4.71];
ybs[7582]=['58 Aql',5.2180874,0.0058297,5.61];
ybs[7583]=['ω Sgr',5.2238106,-0.4579476,4.7];
ybs[7584]=['',5.2175292,0.1256748,6.15];
ybs[7585]=['',5.220853,-0.1164741,6.51];
ybs[7586]=['',5.2084386,0.8354361,6.29];
ybs[7587]=['',5.2162338,0.4255059,5.52];
ybs[7588]=['β Aql',5.2203542,0.1128762,3.71];
ybs[7589]=['μ1 Pav',5.2471184,-1.167383,5.76];
ybs[7590]=['59 Sgr',5.2286731,-0.4731325,4.52];
ybs[7591]=['',5.2324194,-0.663163,6.55];
ybs[7592]=['',5.216907,0.6467573,5.76];
ybs[7593]=['',5.2185603,0.5280632,6.57];
ybs[7594]=['23 Cyg',5.2087161,1.0050156,5.14];
ybs[7595]=['10 Sge',5.2230891,0.2913951,5.36];
ybs[7596]=['φ Aql',5.2242194,0.2004508,5.28];
ybs[7597]=['',5.2097559,1.0431581,6.06];
ybs[7598]=['μ2 Pav',5.2536015,-1.1672786,5.31];
ybs[7599]=['22 Cyg',5.2214417,0.6727812,4.94];
ybs[7600]=['61 Sgr',5.2326002,-0.2692943,5.02];
ybs[7601]=['η Cyg',5.223562,0.6133857,3.89];
ybs[7602]=['',5.2272116,0.3675576,6.48];
ybs[7603]=['',5.2375109,-0.5318994,6.28];
ybs[7604]=['60 Sgr',5.2373836,-0.4561086,4.83];
ybs[7605]=['ψ Cyg',5.2194621,0.9162906,4.92];
ybs[7606]=['',5.2253808,0.6337659,6.02];
ybs[7607]=['',5.2450551,-0.860236,6.17];
ybs[7608]=['11 Sge',5.2306563,0.294105,5.53];
ybs[7609]=['θ1 Sgr',5.2412067,-0.6145921,4.37];
ybs[7610]=['θ2 Sgr',5.2416943,-0.6044926,5.3];
ybs[7611]=['',5.2518061,-1.0351934,5.13];
ybs[7612]=['',5.2177133,1.0177156,6.09];
ybs[7613]=['',5.2447075,-0.7501447,6.14];
ybs[7614]=['',5.2273195,0.7056289,5.45];
ybs[7615]=['',5.2436338,-0.6569264,5.95];
ybs[7616]=['',5.246423,-0.7862651,5.81];
ybs[7617]=['',5.2437537,-0.5871419,5.66];
ybs[7618]=['',5.2244963,0.8894849,6.43];
ybs[7619]=['',5.2200658,1.0281195,4.96];
ybs[7620]=['',5.222039,0.9904384,6.12];
ybs[7621]=['γ Sge',5.2349218,0.3412905,3.47];
ybs[7622]=['',5.2382604,0.0251395,6.17];
ybs[7623]=['',5.2404343,-0.1727092,5.88];
ybs[7624]=['',5.2302869,0.7386699,6.43];
ybs[7625]=['',5.2489268,-0.7112303,6.29];
ybs[7626]=['',5.2339153,0.5418515,5.49];
ybs[7627]=['14 Vul',5.2366102,0.4042856,5.67];
ybs[7628]=['',5.2333041,0.6661521,6.32];
ybs[7629]=['',5.2479089,-0.3957293,6.01];
ybs[7630]=['',5.2696675,-1.1738227,6.07];
ybs[7631]=['13 Sge',5.2406616,0.3068213,5.37];
ybs[7632]=['',5.2361982,0.7999666,5.92];
ybs[7633]=['25 Cyg',5.239241,0.6476143,5.19];
ybs[7634]=['',5.2450179,0.1504719,5.91];
ybs[7635]=['63 Sgr',5.2500995,-0.236895,5.71];
ybs[7636]=['62 Sgr',5.2536032,-0.4825052,4.58];
ybs[7637]=['',5.2353424,0.9096351,6.15];
ybs[7638]=['',5.2580038,-0.661064,4.77];
ybs[7639]=['15 Vul',5.2448328,0.4854979,4.64];
ybs[7640]=['',5.2305584,1.1099621,5.96];
ybs[7641]=['',5.245068,0.6486047,6.2];
ybs[7642]=['',5.2477623,0.4339579,5.88];
ybs[7643]=['16 Vul',5.2489719,0.4363649,5.22];
ybs[7644]=['',5.2581089,-0.3932375,6.45];
ybs[7645]=['',5.261067,-0.558355,4.99];
ybs[7646]=['26 Cyg',5.2446959,0.8755994,5.05];
ybs[7647]=['',5.2588043,-0.1292401,6.72];
ybs[7648]=['',5.2546735,0.3240198,5.96];
ybs[7649]=['',5.2816285,-1.1569385,6.45];
ybs[7650]=['',5.2557515,0.2809266,5.67];
ybs[7651]=['δ Pav',5.2832697,-1.1539199,3.56];
ybs[7652]=['',5.2300056,1.229217,6.33];
ybs[7653]=['62 Aql',5.260182,-0.0112482,5.68];
ybs[7654]=['',5.2663787,-0.5748144,6.53];
ybs[7655]=['τ Aql',5.2588364,0.1281576,5.52];
ybs[7656]=['',5.2557506,0.5229219,5.71];
ybs[7657]=['',5.2636106,-0.2013085,6.34];
ybs[7658]=['15 Sge',5.258348,0.2990587,5.8];
ybs[7659]=['ξ Tel',5.2756951,-0.9217842,4.94];
ybs[7660]=['',5.2767582,-0.9590548,6.26];
ybs[7661]=['65 Sgr',5.2651736,-0.219908,6.55];
ybs[7662]=['64 Dra',5.2434385,1.1324481,5.27];
ybs[7663]=['',5.2619326,0.4062337,6.45];
ybs[7664]=['',5.2599188,0.5634553,5.64];
ybs[7665]=['η Sge',5.262847,0.3500502,5.1];
ybs[7666]=['',5.2642451,0.2716729,6.34];
ybs[7667]=['',5.2682182,-0.0700315,6.47];
ybs[7668]=['65 Dra',5.2472481,1.1291973,6.57];
ybs[7669]=['',5.2620637,0.6727124,6.19];
ybs[7670]=['',5.2584711,0.8429002,6.16];
ybs[7671]=['ρ Dra',5.2486868,1.1857347,4.51];
ybs[7672]=['69 Dra',5.2314095,1.3359393,6.2];
ybs[7673]=['',5.2609447,0.9059066,6.14];
ybs[7674]=['17 Vul',5.2702696,0.4133033,5.07];
ybs[7675]=['27 Cyg',5.267433,0.6289873,5.36];
ybs[7676]=['64 Aql',5.2760778,-0.0106757,5.99];
ybs[7677]=['',5.2924164,-1.0027904,6.37];
ybs[7678]=['',5.2615938,0.9844819,6.21];
ybs[7679]=['',5.2749129,0.1652177,6.43];
ybs[7680]=['',5.2785291,-0.1744609,6.18];
ybs[7681]=['',5.2578959,1.1162338,6.26];
ybs[7682]=['',5.2758435,0.2920082,6.42];
ybs[7683]=['',5.2657148,0.9290654,5.85];
ybs[7684]=['',5.3645382,-1.452734,6.17];
ybs[7685]=['',5.2733017,0.6019552,6.11];
ybs[7686]=['',5.2783602,0.1883691,6.31];
ybs[7687]=['66 Dra',5.2617088,1.0831666,5.39];
ybs[7688]=['',5.2701658,0.8778189,6.54];
ybs[7689]=['',5.2912171,-0.6288932,5.32];
ybs[7690]=['',5.257708,1.1884326,6.28];
ybs[7691]=['θ Sge',5.2837004,0.3662187,6.48];
ybs[7692]=['',5.2968065,-0.7454615,6.22];
ybs[7693]=['',5.3076846,-1.105597,6.09];
ybs[7694]=['28 Cyg',5.2807674,0.6441478,4.93];
ybs[7695]=['',5.2900348,-0.153137,6.49];
ybs[7696]=['θ Aql',5.290372,-0.0131461,3.23];
ybs[7697]=['18 Vul',5.2861546,0.4707492,5.52];
ybs[7698]=['ξ1 Cap',5.2936277,-0.2150945,6.34];
ybs[7699]=['',5.2885596,0.3700576,6.22];
ybs[7700]=['',5.305848,-0.9141317,5.65];
ybs[7701]=['ξ2 Cap',5.2956709,-0.2190178,5.85];
ybs[7702]=['',5.2898127,0.38299,6.26];
ybs[7703]=['',5.2959104,0.0163407,6.27];
ybs[7704]=['19 Vul',5.2915795,0.4690963,5.49];
ybs[7705]=['20 Vul',5.2925168,0.4633384,5.92];
ybs[7706]=['66 Aql',5.2987846,-0.0164129,5.47];
ybs[7707]=['',5.291638,0.8343609,6.92];
ybs[7708]=['',5.3086582,-0.4705887,5.73];
ybs[7709]=['',5.2998735,0.4242563,6.56];
ybs[7710]=['ρ Aql',5.302818,0.2664594,4.95];
ybs[7711]=['',5.3111926,-0.5224641,6.3];
ybs[7712]=['',5.2933925,0.899407,6.01];
ybs[7713]=['68 Dra',5.2880869,1.084665,5.75];
ybs[7714]=['',5.3138889,-0.6350187,6.39];
ybs[7715]=['',5.3140322,-0.6130854,6.53];
ybs[7716]=['30 Cyg',5.297105,0.8182943,4.83];
ybs[7717]=['21 Vul',5.3021631,0.5020297,5.18];
ybs[7718]=['',5.3276851,-1.1023369,6.27];
ybs[7719]=['',5.2991309,0.7583167,6.14];
ybs[7720]=['',5.3011157,0.6400933,6.45];
ybs[7721]=['31 Cyg',5.298559,0.8169977,3.79];
ybs[7722]=['29 Cyg',5.303071,0.643607,4.97];
ybs[7723]=['',5.3020309,0.7360592,6.71];
ybs[7724]=['3 Cap',5.3128848,-0.2140898,6.32];
ybs[7725]=['',5.3067555,0.4478842,4.78];
ybs[7726]=['33 Cyg',5.2967129,0.988498,4.3];
ybs[7727]=['22 Vul',5.3078758,0.4115252,5.15];
ybs[7728]=['',5.2965086,1.0595814,5.79];
ybs[7729]=['',5.3069805,0.5899111,5.66];
ybs[7730]=['23 Vul',5.3088649,0.4866732,4.52];
ybs[7731]=['',5.3255997,-0.8314592,6.31];
ybs[7732]=['18 Sge',5.3115482,0.3781959,6.13];
ybs[7733]=['α1 Cap',5.3184231,-0.2170712,4.24];
ybs[7734]=['4 Cap',5.3203794,-0.3794126,5.87];
ybs[7735]=['',5.3271812,-0.8291777,6.13];
ybs[7736]=['κ Cep',5.2712781,1.3574821,4.39];
ybs[7737]=['32 Cyg',5.3065309,0.8339953,3.98];
ybs[7738]=['',5.3096177,0.6801205,6.27];
ybs[7739]=['24 Vul',5.3134271,0.4318246,5.32];
ybs[7740]=['α2 Cap',5.3201983,-0.2177031,3.57];
ybs[7741]=['',5.3074451,0.8779503,6.31];
ybs[7742]=['',5.3090266,0.7967369,5.91];
ybs[7743]=['',5.3115167,0.6479856,6.48];
ybs[7744]=['',5.3331828,-0.9595533,6.27];
ybs[7745]=['',5.3133177,0.7057352,5.24];
ybs[7746]=['',5.3164804,0.5099677,6.22];
ybs[7747]=['σ Cap',5.3262583,-0.3324284,5.28];
ybs[7748]=['',5.3156456,0.7468759,6.29];
ybs[7749]=['34 Cyg',5.3172181,0.6650419,4.81];
ybs[7750]=['',5.3313252,-0.5083247,6.3];
ybs[7751]=['',5.3333379,-0.6213556,6.46];
ybs[7752]=['',5.3377459,-0.8713814,6.27];
ybs[7753]=['',5.3185173,0.7121538,5.84];
ybs[7754]=['',5.3271041,-0.017569,6.06];
ybs[7755]=['36 Cyg',5.3202833,0.6470172,5.58];
ybs[7756]=['35 Cyg',5.3211444,0.6118115,5.17];
ybs[7757]=['',5.3256244,0.2319333,6.21];
ybs[7758]=['',5.3303819,-0.1097699,6.63];
ybs[7759]=['ν Cap',5.3315831,-0.2214255,4.76];
ybs[7760]=['',5.3278685,0.2377163,5.95];
ybs[7761]=['',5.3321443,-0.2567819,6.1];
ybs[7762]=['β Cap',5.3331691,-0.2567171,3.08];
ybs[7763]=['',5.3212838,0.8097326,6.45];
ybs[7764]=['',5.3293113,0.2555406,6.13];
ybs[7765]=['κ1 Sgr',5.3406097,-0.7326271,5.59];
ybs[7766]=['',5.3292678,0.3118082,5.8];
ybs[7767]=['',5.3187533,0.9681076,5.76];
ybs[7768]=['',5.3260783,0.6493402,6.57];
ybs[7769]=['',5.3132453,1.1680558,5.93];
ybs[7770]=['',5.3279212,0.688977,6.23];
ybs[7771]=['',5.3969541,-1.4117367,5.77];
ybs[7772]=['',5.326092,0.8187249,6.5];
ybs[7773]=['κ2 Sgr',5.3468557,-0.7391272,5.64];
ybs[7774]=['',5.3417342,-0.1672246,6.3];
ybs[7775]=['25 Vul',5.3364507,0.4279389,5.54];
ybs[7776]=['α Pav',5.355602,-0.9889081,1.94];
ybs[7777]=['',5.3280708,0.9366889,6.18];
ybs[7778]=['71 Dra',5.3231793,1.0878508,5.72];
ybs[7779]=['',5.3403547,0.25525,6.17];
ybs[7780]=['',5.3419873,0.0945369,5.31];
ybs[7781]=['',5.3356765,0.7191512,6.39];
ybs[7782]=['γ Cyg',5.3365028,0.7038859,2.2];
ybs[7783]=['',5.3386538,0.5469551,6.09];
ybs[7784]=['',5.3355757,0.8005466,5.58];
ybs[7785]=['',5.355098,-0.7107267,6.09];
ybs[7786]=['',5.3387619,0.7173192,5.93];
ybs[7787]=['',5.3530103,-0.498968,5.85];
ybs[7788]=['',5.3393945,0.7514804,6.2];
ybs[7789]=['',5.3484323,0.0199451,6.15];
ybs[7790]=['',5.3240793,1.2034427,5.55];
ybs[7791]=['',5.329819,1.117931,5.69];
ybs[7792]=['39 Cyg',5.3440217,0.5631089,4.43];
ybs[7793]=['',5.3432476,0.6553727,5.9];
ybs[7794]=['',5.3596551,-0.6514937,6.25];
ybs[7795]=['',5.353278,-0.0475714,6.11];
ybs[7796]=['',5.3529898,0.1768196,6.33];
ybs[7797]=['',5.3523578,0.3749719,5.66];
ybs[7798]=['',5.4190526,-1.4173525,5.91];
ybs[7799]=['',5.3539179,0.3480141,6.41];
ybs[7800]=['π Cap',5.3607987,-0.3165381,5.25];
ybs[7801]=['',5.3457038,0.9359493,6.51];
ybs[7802]=['',5.3556055,0.3035208,6.22];
ybs[7803]=['',5.3678305,-0.6199374,6.1];
ybs[7804]=['',5.3474669,1.0415114,6.44];
ybs[7805]=['',5.3668519,-0.2734178,6.41];
ybs[7806]=['',5.3634783,0.148583,6.25];
ybs[7807]=['68 Aql',5.3651043,-0.0572809,6.13];
ybs[7808]=['ρ Cap',5.367502,-0.3095791,4.78];
ybs[7809]=['',5.3581934,0.6004645,6.39];
ybs[7810]=['',5.3643237,0.0525816,6.21];
ybs[7811]=['',5.3705449,-0.389476,6.16];
ybs[7812]=['40 Cyg',5.3599357,0.672225,5.62];
ybs[7813]=['',5.3535075,0.9898406,6.36];
ybs[7814]=['43 Cyg',5.3569537,0.8632128,5.69];
ybs[7815]=['ο Cap',5.3719547,-0.3230637,6.74];
ybs[7816]=['ο Cap',5.3720563,-0.3230054,5.94];
ybs[7817]=['69 Aql',5.3704782,-0.0490297,4.91];
ybs[7818]=['',5.3770089,-0.5067658,6.39];
ybs[7819]=['',5.3684544,0.3519277,6.55];
ybs[7820]=['41 Cyg',5.368261,0.5313621,4.01];
ybs[7821]=['42 Cyg',5.3677553,0.6375841,5.88];
ybs[7822]=['1 Del',5.372894,0.1915054,6.08];
ybs[7823]=['',5.3770579,-0.2614399,6.12];
ybs[7824]=['',5.4020417,-1.2135601,6.11];
ybs[7825]=['',5.3754941,0.3609816,6.18];
ybs[7826]=['',5.3768834,0.1978779,7.11];
ybs[7827]=['',5.3701201,0.8029392,6.41];
ybs[7828]=['',5.3852521,-0.4339955,6.36];
ybs[7829]=['',5.3669426,0.9799011,5.91];
ybs[7830]=['ω1 Cyg',5.3701819,0.8557018,4.95];
ybs[7831]=['',5.3826663,-0.1706196,5.65];
ybs[7832]=['ν Mic',5.3906921,-0.7755866,5.11];
ybs[7833]=['44 Cyg',5.3749247,0.6459935,6.19];
ybs[7834]=['φ1 Pav',5.3992694,-1.0559699,4.76];
ybs[7835]=['',5.3796818,0.4517219,6.34];
ybs[7836]=['θ Cep',5.3666823,1.1007847,4.22];
ybs[7837]=['ω2 Cyg',5.3756473,0.8603993,5.44];
ybs[7838]=['ε Del',5.3855965,0.1986395,4.03];
ybs[7839]=['',5.39474,-0.6634178,6.44];
ybs[7840]=['',5.3755815,0.9143204,6.18];
ybs[7841]=['',5.3906387,-0.2381115,6.13];
ybs[7842]=['',5.3938195,-0.5304928,6.4];
ybs[7843]=['',5.3886028,0.1769395,6.56];
ybs[7844]=['η Del',5.3887579,0.2287325,5.38];
ybs[7845]=['ρ Pav',5.4081294,-1.0725066,4.88];
ybs[7846]=['',5.3770299,0.9923443,6.14];
ybs[7847]=['',5.3828367,0.7551923,6.6];
ybs[7848]=['',5.3894436,0.3676282,6.48];
ybs[7849]=['μ1 Oct',5.4311924,-1.3281716,6];
ybs[7850]=['μ2 Oct',5.4294106,-1.313688,6.55];
ybs[7851]=['',5.39658,-0.2870531,6.19];
ybs[7852]=['47 Cyg',5.3877412,0.6166066,4.61];
ybs[7853]=['',5.3869996,0.7304253,6.49];
ybs[7854]=['',5.3664448,1.267247,6.27];
ybs[7855]=['α Ind',5.4067578,-0.8239972,3.11];
ybs[7856]=['',5.3871857,0.8163252,5.78];
ybs[7857]=['ζ Del',5.3946347,0.2574871,4.68];
ybs[7858]=['',5.4181856,-1.0965366,6.22];
ybs[7859]=['70 Aql',5.401344,-0.0431203,4.89];
ybs[7860]=['26 Vul',5.3978787,0.4531153,6.41];
ybs[7861]=['φ2 Pav',5.4186923,-1.0553652,5.12];
ybs[7862]=['',5.3908649,0.9063956,6.11];
ybs[7863]=['',5.4070388,-0.4368382,6.36];
ybs[7864]=['',5.4037892,0.0030818,6.22];
ybs[7865]=['73 Dra',5.3720524,1.3095478,5.2];
ybs[7866]=['27 Vul',5.4019617,0.4632355,5.59];
ybs[7867]=['υ Pav',5.427918,-1.1637695,5.15];
ybs[7868]=['β Del',5.4044218,0.2561269,3.63];
ybs[7869]=['ι Del',5.4056929,0.1999729,5.43];
ybs[7870]=['71 Aql',5.4083333,-0.0178933,4.32];
ybs[7871]=['48 Cyg',5.4037367,0.5524346,6.32];
ybs[7872]=['',5.4058757,0.3202508,6.25];
ybs[7873]=['',5.403797,0.5515523,6.49];
ybs[7874]=['',5.4028423,0.6703498,6.2];
ybs[7875]=['τ Cap',5.412828,-0.2596045,5.22];
ybs[7876]=['',5.4122185,-0.040707,6.22];
ybs[7877]=['29 Vul',5.4084486,0.3714273,4.82];
ybs[7878]=['θ Del',5.4096196,0.2337905,5.72];
ybs[7879]=['',5.4180921,-0.5820847,5.47];
ybs[7880]=['28 Vul',5.4083918,0.4223036,5.04];
ybs[7881]=['',5.40864,0.4147022,5.91];
ybs[7882]=['κ Del',5.4114569,0.1774387,5.05];
ybs[7883]=['1 Aqr',5.4129844,0.0098943,5.16];
ybs[7884]=['',5.4171373,-0.4135211,6.37];
ybs[7885]=['',5.4110757,0.2778286,5.97];
ybs[7886]=['υ Cap',5.4163165,-0.3151682,5.1];
ybs[7887]=['75 Dra',5.3526066,1.4224074,5.46];
ybs[7888]=['',5.4190053,-0.4636284,6.51];
ybs[7889]=['',5.4112873,0.3821853,6.08];
ybs[7890]=['',5.4101692,0.5308372,5.68];
ybs[7891]=['',5.4183966,-0.2800059,5.8];
ybs[7892]=['α Del',5.4134959,0.2791223,3.77];
ybs[7893]=['',5.4146223,0.1977529,6.42];
ybs[7894]=['74 Dra',5.3583963,1.4166332,5.96];
ybs[7895]=['',5.4226599,-0.5500745,5.76];
ybs[7896]=['',5.4224775,-0.4523652,6.28];
ybs[7897]=['',5.4121741,0.70965,6.06];
ybs[7898]=['',5.4111508,0.7984421,6.58];
ybs[7899]=['β Pav',5.4408936,-1.1540129,3.42];
ybs[7900]=['',5.4182335,0.3493507,6.45];
ybs[7901]=['',5.4294963,-0.6889965,6.29];
ybs[7902]=['',5.40871,0.9788717,6.48];
ybs[7903]=['',5.4172278,0.5216133,6.08];
ybs[7904]=['10 Del',5.4206597,0.2559407,5.99];
ybs[7905]=['',5.4141912,0.7599045,5.95];
ybs[7906]=['η Ind',5.4352689,-0.9047538,4.51];
ybs[7907]=['49 Cyg',5.41905,0.5652836,5.51];
ybs[7908]=['',5.4185898,0.6835292,6.51];
ybs[7909]=['',5.4236229,0.3072293,6.22];
ybs[7910]=['α Cyg',5.4201169,0.7917084,1.25];
ybs[7911]=['',5.4138943,1.0574255,6.01];
ybs[7912]=['',5.4225423,0.7295204,5.67];
ybs[7913]=['',5.4247114,0.6202516,6.66];
ybs[7914]=['δ Del',5.4301965,0.2645329,4.43];
ybs[7915]=['51 Cyg',5.4231913,0.8800226,5.39];
ybs[7916]=['',5.3520494,1.4608557,6.19];
ybs[7917]=['',5.439171,-0.4741058,6.5];
ybs[7918]=['',5.429219,0.6225571,6.47];
ybs[7919]=['',5.444522,-0.6826983,5.5];
ybs[7920]=['σ Pav',5.4602743,-1.1988949,5.41];
ybs[7921]=['',5.4442734,-0.6289618,6.49];
ybs[7922]=['ψ Cap',5.4429161,-0.4396051,4.14];
ybs[7923]=['17 Cap',5.4431016,-0.3740385,5.93];
ybs[7924]=['',5.4242206,1.0591202,6.15];
ybs[7925]=['30 Vul',5.4360473,0.4424985,4.91];
ybs[7926]=['',5.4270635,0.9982611,6.32];
ybs[7927]=['',5.4388836,0.3171835,6.38];
ybs[7928]=['52 Cyg',5.439282,0.5376099,4.22];
ybs[7929]=['ι Mic',5.4541409,-0.7662744,5.11];
ybs[7930]=['',5.4321162,0.9873416,5.78];
ybs[7931]=['4 Cep',5.4256015,1.1648222,5.58];
ybs[7932]=['',5.4464269,-0.0419401,6.27];
ybs[7933]=['γ1 Del',5.4440742,0.2828816,5.14];
ybs[7934]=['γ2 Del',5.4441323,0.2828769,4.27];
ybs[7935]=['ε Cyg',5.4415538,0.5943466,2.46];
ybs[7936]=['ε Aqr',5.4493185,-0.1642687,3.77];
ybs[7937]=['3 Aqr',5.4494516,-0.0862861,4.42];
ybs[7938]=['ζ Ind',5.4586095,-0.8053336,4.89];
ybs[7939]=['13 Del',5.4494362,0.1063307,5.58];
ybs[7940]=['',5.4494769,0.0591777,6.4];
ybs[7941]=['',5.4363175,1.0064017,4.51];
ybs[7942]=['',5.4457663,0.6014027,4.92];
ybs[7943]=['η Cep',5.4355646,1.0807374,3.43];
ybs[7944]=['',5.4428091,0.8135867,6.3];
ybs[7945]=['',5.4693839,-1.0880995,6.28];
ybs[7946]=['',5.4694202,-1.0880995,6.59];
ybs[7947]=['',5.4568795,-0.4484934,5.86];
ybs[7948]=['',5.4411039,0.9263956,6.33];
ybs[7949]=['λ Cyg',5.446673,0.638347,4.53];
ybs[7950]=['',5.4568367,-0.3133078,6.21];
ybs[7951]=['α Mic',5.4601211,-0.5880857,4.9];
ybs[7952]=['',5.4459462,0.7969772,6.4];
ybs[7953]=['',5.4309404,1.2188398,6.41];
ybs[7954]=['ι Ind',5.4677503,-0.8992421,5.05];
ybs[7955]=['',5.4478857,0.8362892,5.57];
ybs[7956]=['',5.4636125,-0.5579682,6.36];
ybs[7957]=['',5.4648446,-0.6602234,5.52];
ybs[7958]=['',5.4478237,0.916143,6.27];
ybs[7959]=['15 Del',5.4572115,0.2204346,5.98];
ybs[7960]=['14 Del',5.4581113,0.1387351,6.33];
ybs[7961]=['',5.4589629,0.0982545,6.21];
ybs[7962]=['',5.462581,-0.2174654,5.88];
ybs[7963]=['55 Cyg',5.4528704,0.8063162,4.84];
ybs[7964]=['',5.4514968,0.9074804,6.29];
ybs[7965]=['β Mic',5.4688659,-0.5775561,6.04];
ybs[7966]=['ω Cap',5.4679429,-0.4683337,4.11];
ybs[7967]=['',5.4613521,0.3165412,6.52];
ybs[7968]=['4 Aqr',5.4655742,-0.0967078,5.99];
ybs[7969]=['',5.4570915,0.8158691,6.33];
ybs[7970]=['56 Cyg',5.4579829,0.7704628,5.04];
ybs[7971]=['5 Aqr',5.4686974,-0.0946181,5.55];
ybs[7972]=['β Ind',5.4827341,-1.0187012,3.65];
ybs[7973]=['',5.4765166,-0.6933079,5.35];
ybs[7974]=['',5.4647389,0.4945559,5.77];
ybs[7975]=['',5.4730612,-0.4135899,6.33];
ybs[7976]=['μ Aqr',5.4710133,-0.1552889,4.73];
ybs[7977]=['',5.4750339,-0.5346353,6.35];
ybs[7978]=['',5.4811096,-0.8838523,6.24];
ybs[7979]=['',5.4527637,1.119221,6.45];
ybs[7980]=['',5.4730068,-0.2004947,6.38];
ybs[7981]=['31 Vul',5.4676547,0.4744262,4.59];
ybs[7982]=['',5.466902,0.5748204,6.44];
ybs[7983]=['',5.4779697,-0.4858826,6.41];
ybs[7984]=['',5.4767125,-0.1187394,6.44];
ybs[7985]=['',5.4719079,0.5189825,6.34];
ybs[7986]=['19 Cap',5.4806309,-0.3113015,5.78];
ybs[7987]=['57 Cyg',5.471795,0.7762057,4.78];
ybs[7988]=['76 Dra',5.4140737,1.4418573,5.75];
ybs[7989]=['',5.4720279,0.7900766,5.45];
ybs[7990]=['',5.4727437,0.7417029,6.66];
ybs[7991]=['',5.475148,0.5851067,5.47];
ybs[7992]=['',5.4816279,-0.0224524,6.55];
ybs[7993]=['',5.4774076,0.4993126,6.56];
ybs[7994]=['32 Vul',5.4782393,0.4912079,5.01];
ybs[7995]=['',5.4768988,0.7119127,6.7];
ybs[7996]=['',5.4838363,0.0806323,6.05];
ybs[7997]=['17 Del',5.4832861,0.2410031,5.17];
ybs[7998]=['16 Del',5.4834577,0.2208835,5.58];
ybs[7999]=['',5.4895762,-0.45743,5.7];
ybs[8000]=['',5.4867927,-0.0606332,6.57];
ybs[8001]=['7 Aqr',5.4895583,-0.1677243,5.51];
ybs[8002]=['',5.4386787,1.4073573,5.39];
ybs[8003]=['',5.4904859,0.0096222,6.05];
ybs[8004]=['',5.4931276,-0.278271,5.87];
ybs[8005]=['',5.5130419,-1.1889216,6.37];
ybs[8006]=['',5.4829114,0.8291164,5.67];
ybs[8007]=['α Oct',5.5299054,-1.3427344,5.15];
ybs[8008]=['',5.4843581,0.8929495,6.63];
ybs[8009]=['',5.4863246,0.7856146,5.96];
ybs[8010]=['',5.4975385,-0.2512309,6.01];
ybs[8011]=['',5.4853,0.8869054,5.81];
ybs[8012]=['',5.4854308,0.8601535,5.9];
ybs[8013]=['',5.5062969,-0.8931942,5.76];
ybs[8014]=['ν Cyg',5.4891172,0.7200332,3.94];
ybs[8015]=['',5.4841925,0.9943969,6.23];
ybs[8016]=['18 Del',5.4956772,0.1907181,5.48];
ybs[8017]=['',5.5039203,-0.629032,6.11];
ybs[8018]=['33 Vul',5.4946452,0.3911969,5.31];
ybs[8019]=['20 Cap',5.5016123,-0.3306808,6.25];
ybs[8020]=['ε Equ',5.4986581,0.0764812,5.23];
ybs[8021]=['',5.4939819,0.7777143,5.55];
ybs[8022]=['',5.4949344,0.7335346,6.16];
ybs[8023]=['',5.5016818,0.2951855,6.66];
ybs[8024]=['',5.5028979,0.1327358,5.99];
ybs[8025]=['γ Mic',5.5094204,-0.5614452,4.67];
ybs[8026]=['',5.4943961,0.8822701,5.61];
ybs[8027]=['11 Aqr',5.5054029,-0.0810052,6.21];
ybs[8028]=['',5.5138826,-0.7489598,6.64];
ybs[8029]=['',5.4735302,1.3266599,6.05];
ybs[8030]=['',5.5042919,0.3389149,5.65];
ybs[8031]=['',5.5112504,-0.4676017,6.05];
ybs[8032]=['',5.5147355,-0.670918,5.94];
ybs[8033]=['59 Cyg',5.5003731,0.8309472,4.74];
ybs[8034]=['ζ Mic',5.516977,-0.6726792,5.3];
ybs[8035]=['',5.4976766,1.0389434,5.51];
ybs[8036]=['',5.5174481,-0.4824421,6.25];
ybs[8037]=['',5.5070024,0.6303314,5.97];
ybs[8038]=['',5.5472052,-1.3285468,6.58];
ybs[8039]=['60 Cyg',5.5063739,0.8071278,5.37];
ybs[8040]=['',5.5159031,-0.0145696,6.5];
ybs[8041]=['μ Ind',5.5277916,-0.953584,5.16];
ybs[8042]=['',5.5160852,0.0283077,6.25];
ybs[8043]=['',5.5156384,0.2586567,6.31];
ybs[8044]=['12 Aqr',5.5207599,-0.1000592,7.31];
ybs[8045]=['12 Aqr',5.5207671,-0.1000544,5.89];
ybs[8046]=['η Cap',5.5225838,-0.3449555,4.84];
ybs[8047]=['',5.5485684,-1.2754956,5.68];
ybs[8048]=['',5.51178,0.7833171,6.19];
ybs[8049]=['',5.5150385,0.6762702,6.07];
ybs[8050]=['',5.5135001,0.8017815,6.48];
ybs[8051]=['',5.509878,0.9906356,5.83];
ybs[8052]=['3 Equ',5.5226489,0.0976218,5.61];
ybs[8053]=['',5.5232258,0.0529276,6.42];
ybs[8054]=['',5.5235127,0.0411955,6.33];
ybs[8055]=['η Mic',5.5321702,-0.7207305,5.53];
ybs[8056]=['δ Mic',5.5299549,-0.52419,5.68];
ybs[8057]=['',5.5183476,0.7281212,6.33];
ybs[8058]=['',5.5159624,0.8803785,6.37];
ybs[8059]=['',5.5431383,-1.1141611,5.76];
ybs[8060]=['',5.5174352,0.8194686,6.32];
ybs[8061]=['θ Cap',5.5292304,-0.2991791,4.07];
ybs[8062]=['',5.5317415,-0.5628756,5.18];
ybs[8063]=['4 Equ',5.5264189,0.1055784,5.94];
ybs[8064]=['',5.5173319,0.9315917,5.9];
ybs[8065]=['ξ Cyg',5.522881,0.7682659,3.72];
ybs[8066]=['24 Cap',5.5346168,-0.4348366,4.5];
ybs[8067]=['',5.5569496,-1.264507,6.2];
ybs[8068]=['',5.5299451,0.4715118,6.12];
ybs[8069]=['',5.5370722,-0.3030507,6.17];
ybs[8070]=['',5.5302939,0.5458684,5.82];
ybs[8071]=['61 Cyg',5.5317673,0.6778371,5.21];
ybs[8072]=['61 Cyg',5.5318184,0.6777935,6.03];
ybs[8073]=['χ Cap',5.5407424,-0.3682916,5.3];
ybs[8074]=['',5.5353866,0.2748938,6.34];
ybs[8075]=['63 Cyg',5.529978,0.8332126,4.55];
ybs[8076]=['',5.5395931,0.1235944,6.15];
ybs[8077]=['27 Cap',5.5450407,-0.3571636,6.25];
ybs[8078]=['ο Pav',5.5650228,-1.2222966,5.02];
ybs[8079]=['ν Aqr',5.5449805,-0.1968598,4.51];
ybs[8080]=['',5.5396893,0.5287974,5.59];
ybs[8081]=['',5.5462515,0.0529862,6.45];
ybs[8082]=['',5.5501032,-0.1616354,6.27];
ybs[8083]=['γ Equ',5.5476789,0.1784485,4.69];
ybs[8084]=['6 Equ',5.5484594,0.1770049,6.07];
ybs[8085]=['',5.5262107,1.2483114,5.87];
ybs[8086]=['',5.5573993,-0.7012037,5.83];
ybs[8087]=['',5.5481704,0.3935274,6.68];
ybs[8088]=['',5.5541946,-0.2509612,6.48];
ybs[8089]=['',5.5448475,0.7957827,6.63];
ybs[8090]=['',5.5609579,-0.686465,5.26];
ybs[8091]=['',5.5500201,0.6351616,6.54];
ybs[8092]=['',5.5455883,0.9364723,5.73];
ybs[8093]=['',5.5470882,0.8339992,6.46];
ybs[8094]=['',5.5619836,-0.6340793,5.96];
ybs[8095]=['',5.5412445,1.1063259,6.54];
ybs[8096]=['',5.5615568,-0.4804131,5.42];
ybs[8097]=['',5.5879008,-1.3133759,6.63];
ybs[8098]=['',5.5195552,1.3651436,5.91];
ybs[8099]=['',5.5406536,1.1969907,7.33];
ybs[8100]=['',5.5735196,-0.9279622,5.75];
ybs[8101]=['ζ Cyg',5.5584177,0.5291934,3.2];
ybs[8102]=['',5.5612202,0.2805847,6.27];
ybs[8103]=['',5.5705612,-0.7053201,6.21];
ybs[8104]=['',5.5653887,-0.183454,6.77];
ybs[8105]=['',5.5517846,1.0485853,5.64];
ybs[8106]=['',5.560375,0.6410138,6.05];
ybs[8107]=['',5.5665926,0.0032545,6.38];
ybs[8108]=['',5.569196,-0.3010791,6.04];
ybs[8109]=['δ Equ',5.5657515,0.176298,4.49];
ybs[8110]=['',5.572718,-0.6303453,6.12];
ybs[8111]=['',5.5843612,-1.1272398,6.31];
ybs[8112]=['',5.5638225,0.5235143,6.17];
ybs[8113]=['φ Cap',5.5715634,-0.3587879,5.24];
ybs[8114]=['29 Cap',5.5719209,-0.2631384,5.28];
ybs[8115]=['',5.6569069,-1.4784547,6.45];
ybs[8116]=['τ Cyg',5.5662311,0.6656654,3.72];
ybs[8117]=['α Equ',5.5717286,0.0932434,3.92];
ybs[8118]=['',5.5755507,-0.0264033,6.48];
ybs[8119]=['',5.5596086,1.1256969,6.39];
ybs[8120]=['',5.5783112,-0.2300988,6.4];
ybs[8121]=['ε Mic',5.5819784,-0.5598498,4.71];
ybs[8122]=['',5.5693509,0.8389475,6.46];
ybs[8123]=['30 Cap',5.5816206,-0.3122362,5.43];
ybs[8124]=['',5.5734856,0.7390815,6.43];
ybs[8125]=['31 Cap',5.5829438,-0.3031054,7.05];
ybs[8126]=['θ Ind',5.5913894,-0.9311948,4.39];
ybs[8127]=['15 Aqr',5.5822765,-0.077212,5.82];
ybs[8128]=['',5.5860864,-0.5003816,6.4];
ybs[8129]=['σ Cyg',5.5776445,0.6892292,4.23];
ybs[8130]=['',5.5773696,0.746626,6.19];
ybs[8131]=['',5.5921737,-0.7841058,6];
ybs[8132]=['υ Cyg',5.5800121,0.6107314,4.43];
ybs[8133]=['',5.5751772,0.9440928,6.13];
ybs[8134]=['',5.5897459,-0.4582703,6.56];
ybs[8135]=['',5.5848598,0.1972063,5.96];
ybs[8136]=['',5.5759379,0.9755196,5.98];
ybs[8137]=['θ1 Mic',5.5946035,-0.7105802,4.82];
ybs[8138]=['',5.5972918,-0.8698914,6.38];
ybs[8139]=['',5.5760265,1.0246267,6.42];
ybs[8140]=['68 Cyg',5.581972,0.7686673,5];
ybs[8141]=['',5.5841507,0.7179685,6.15];
ybs[8142]=['',5.612444,-1.2153841,6.41];
ybs[8143]=['',5.5862261,0.669044,5.83];
ybs[8144]=['',5.5905357,0.3861123,6.29];
ybs[8145]=['',5.6173377,-1.2514236,6.09];
ybs[8146]=['16 Aqr',5.5948713,-0.0779021,5.87];
ybs[8147]=['',5.5861665,0.8657914,5.76];
ybs[8148]=['α Cep',5.5811582,1.0939917,2.44];
ybs[8149]=['9 Equ',5.5946167,0.1300441,5.82];
ybs[8150]=['',5.5845572,1.0248472,5.66];
ybs[8151]=['',5.5941509,0.4180472,5.57];
ybs[8152]=['',5.5928395,0.5680906,5.68];
ybs[8153]=['ι Cap',5.6003058,-0.2921243,4.28];
ybs[8154]=['',5.5650817,1.3457641,5.95];
ybs[8155]=['',5.5951563,0.5708864,6.04];
ybs[8156]=['',5.5933638,0.7058465,6.4];
ybs[8157]=['6 Cep',5.5843275,1.1339013,5.18];
ybs[8158]=['',5.603786,-0.3939498,5.6];
ybs[8159]=['1 Peg',5.5986896,0.3473434,4.08];
ybs[8160]=['',5.5516098,1.4193756,6.15];
ybs[8161]=['17 Aqr',5.6031314,-0.1609587,5.99];
ybs[8162]=['',5.6607255,-1.441325,6.38];
ybs[8163]=['',5.6104959,-0.8118796,6.31];
ybs[8164]=['β Equ',5.6025439,0.120572,5.16];
ybs[8165]=['',5.5900176,1.0620887,6.11];
ybs[8166]=['θ2 Mic',5.6105291,-0.7139955,5.77];
ybs[8167]=['γ Pav',5.6211389,-1.1391349,4.22];
ybs[8168]=['',5.6010496,0.5306985,6.05];
ybs[8169]=['33 Cap',5.6087572,-0.3622313,5.41];
ybs[8170]=['',5.6086862,-0.3953054,6.38];
ybs[8171]=['',5.5972198,0.8636881,5.69];
ybs[8172]=['',5.6011165,0.6759878,6.63];
ybs[8173]=['18 Aqr',5.6086869,-0.2230606,5.49];
ybs[8174]=['γ Ind',5.61929,-0.9522895,6.12];
ybs[8175]=['',5.603793,0.6545676,6.58];
ybs[8176]=['',5.6068311,0.4253661,5.71];
ybs[8177]=['',5.6090808,0.1792775,6.35];
ybs[8178]=['20 Aqr',5.6113782,-0.0576044,6.36];
ybs[8179]=['',5.6056374,0.6536053,6.47];
ybs[8180]=['',5.6074268,0.4434844,6.15];
ybs[8181]=['19 Aqr',5.6130875,-0.1684355,5.7];
ybs[8182]=['',5.6318138,-1.211363,5.34];
ybs[8183]=['',5.6085923,0.4298143,6.32];
ybs[8184]=['',5.6093373,0.4585356,5.68];
ybs[8185]=['21 Aqr',5.6132217,-0.0603653,5.49];
ybs[8186]=['',5.6189793,-0.658531,5.63];
ybs[8187]=['',5.6556402,-1.3951852,6.47];
ybs[8188]=['',5.62198,-0.7408775,5.51];
ybs[8189]=['',5.6156333,0.0110414,6.46];
ybs[8190]=['ζ Cap',5.61973,-0.3894338,3.74];
ybs[8191]=['',5.6182815,0.020974,6.13];
ybs[8192]=['',5.6099825,0.8625615,6.58];
ybs[8193]=['35 Cap',5.6222259,-0.3682198,5.78];
ybs[8194]=['',5.6118778,0.8170303,5.6];
ybs[8195]=['69 Cyg',5.6143062,0.641681,5.94];
ybs[8196]=['',5.6177295,0.3398839,6.07];
ybs[8197]=['',5.6311881,-0.9356104,6.39];
ybs[8198]=['',5.6262758,-0.2001778,6.61];
ybs[8199]=['36 Cap',5.6286811,-0.378877,4.51];
ybs[8200]=['5 PsA',5.6304385,-0.5434837,6.5];
ybs[8201]=['70 Cyg',5.6211634,0.6495297,5.31];
ybs[8202]=['',5.6184693,0.8540499,5.31];
ybs[8203]=['35 Vul',5.62284,0.483585,5.41];
ybs[8204]=['',5.6177462,0.9249725,6.03];
ybs[8205]=['',5.6266014,0.1447679,6.4];
ybs[8206]=['',5.6247484,0.5641636,5.8];
ybs[8207]=['',5.6289155,0.3142474,6.44];
ybs[8208]=['',5.6341486,-0.3324537,6.57];
ybs[8209]=['',5.6288025,0.3888359,5.93];
ybs[8210]=['',5.6201524,1.0445554,6.1];
ybs[8211]=['2 Peg',5.6329107,0.4143135,4.57];
ybs[8212]=['',5.6268817,0.9689672,6.12];
ybs[8213]=['7 Cep',5.6208315,1.1677626,5.44];
ybs[8214]=['71 Cyg',5.6299107,0.8140196,5.24];
ybs[8215]=['ξ Gru',5.6440178,-0.7169614,5.29];
ybs[8216]=['6 PsA',5.6443862,-0.5906959,5.97];
ybs[8217]=['',5.6384954,0.2135836,6.08];
ybs[8218]=['β Aqr',5.6406451,-0.0954873,2.91];
ybs[8219]=['',5.6498013,-0.9186852,6.41];
ybs[8220]=['',5.6793954,-1.3847408,6.18];
ybs[8221]=['',5.6454634,-0.4274333,6.43];
ybs[8222]=['',5.649808,-0.7809978,5.57];
ybs[8223]=['',5.6334426,0.926031,6.02];
ybs[8224]=['β Cep',5.6240645,1.2332464,3.23];
ybs[8225]=['',5.6028863,1.4071233,5.97];
ybs[8226]=['',5.6438486,0.4100616,6.7];
ybs[8227]=['',5.6536084,-0.7474195,6.32];
ybs[8228]=['',5.6383548,0.9201374,6.16];
ybs[8229]=['',5.6356952,1.0569584,5.53];
ybs[8230]=['',5.6557595,-0.5165287,6.41];
ybs[8231]=['37 Cap',5.6553508,-0.3487739,5.69];
ybs[8232]=['',5.645004,0.8740305,5.75];
ybs[8233]=['',5.6572479,-0.4075842,6.4];
ybs[8234]=['',5.6467659,0.8020621,6.25];
ybs[8235]=['',5.6715153,-1.1296101,6.2];
ybs[8236]=['',5.6531011,0.3989084,6.47];
ybs[8237]=['',5.6569035,-0.0677494,5.77];
ybs[8238]=['ρ Cyg',5.6497651,0.7974892,4.02];
ybs[8239]=['8 PsA',5.6613254,-0.4550035,5.73];
ybs[8240]=['ν Oct',5.6892759,-1.3489049,3.76];
ybs[8241]=['72 Cyg',5.6535307,0.6743125,4.9];
ybs[8242]=['7 PsA',5.6642816,-0.5750204,6.11];
ybs[8243]=['',5.6562257,0.4939021,6.31];
ybs[8244]=['',5.6569128,0.4285403,6.11];
ybs[8245]=['',5.6515305,0.9040683,6.15];
ybs[8246]=['ε Cap',5.6650592,-0.3379695,4.68];
ybs[8247]=['',5.6601675,0.5263412,6.36];
ybs[8248]=['',5.6587724,0.7937096,5.53];
ybs[8249]=['',5.6667273,-0.005031,6.25];
ybs[8250]=['ξ Aqr',5.6677152,-0.1352993,4.69];
ybs[8251]=['3 Peg',5.6672916,0.1172932,6.18];
ybs[8252]=['74 Cyg',5.6629567,0.7071271,5.01];
ybs[8253]=['5 Peg',5.6671179,0.3389548,5.45];
ybs[8254]=['',5.6742673,-0.5860227,6.28];
ybs[8255]=['',5.6789453,-0.9120449,6.21];
ybs[8256]=['4 Peg',5.6708164,0.1025205,5.67];
ybs[8257]=['',5.6816171,-0.9710047,6.33];
ybs[8258]=['',5.6650315,0.7818833,6.2];
ybs[8259]=['',5.6752648,-0.1828112,6.08];
ybs[8260]=['',5.6712954,0.4468264,6.16];
ybs[8261]=['',5.665323,0.9449947,6.15];
ybs[8262]=['',5.6726028,0.3554842,5.85];
ybs[8263]=['25 Aqr',5.6753599,0.04095,5.1];
ybs[8264]=['γ Cap',5.6781242,-0.2890158,3.68];
ybs[8265]=['9 Cep',5.6658783,1.0853154,4.73];
ybs[8266]=['λ Oct',5.7339032,-1.4418652,5.29];
ybs[8267]=['',5.670834,1.0051623,5.62];
ybs[8268]=['',5.6856449,-0.4363076,6.49];
ybs[8269]=['42 Cap',5.684416,-0.2433725,5.18];
ybs[8270]=['75 Cyg',5.6769809,0.7570661,5.11];
ybs[8271]=['41 Cap',5.6866655,-0.4042068,5.24];
ybs[8272]=['',5.7047482,-1.2375134,6.01];
ybs[8273]=['26 Aqr',5.6867842,0.024238,5.67];
ybs[8274]=['κ Cap',5.6893694,-0.3274719,4.73];
ybs[8275]=['7 Peg',5.687082,0.1009409,5.3];
ybs[8276]=['',5.678727,0.9594977,6.2];
ybs[8277]=['76 Cyg',5.6831387,0.7139883,6.11];
ybs[8278]=['',5.688242,0.1907347,6.09];
ybs[8279]=['',5.6918585,-0.3406364,6.22];
ybs[8280]=['',5.9935758,-1.5477592,6.57];
ybs[8281]=['44 Cap',5.6910751,-0.2495118,5.88];
ybs[8282]=['',5.6852759,0.6215757,6.07];
ybs[8283]=['',5.6854198,0.8005693,6.17];
ybs[8284]=['',5.6979155,-0.6710493,6.3];
ybs[8285]=['77 Cyg',5.6866682,0.7187391,5.69];
ybs[8286]=['π1 Cyg',5.6849554,0.8952336,4.67];
ybs[8287]=['45 Cap',5.6951972,-0.2556106,5.99];
ybs[8288]=['',5.7019598,-0.8620903,6.45];
ybs[8289]=['',5.6874595,0.8674957,6.09];
ybs[8290]=['ι PsA',5.6997206,-0.5745885,4.34];
ybs[8291]=['',5.6898337,0.7201004,5.49];
ybs[8292]=['79 Cyg',5.6913433,0.6699918,5.65];
ybs[8293]=['ε Peg',5.6954066,0.1741677,2.39];
ybs[8294]=['μ1 Cyg',5.6947614,0.503472,4.73];
ybs[8295]=['μ2 Cyg',5.6947395,0.5034769,6.08];
ybs[8296]=['46 Cap',5.6993825,-0.1566985,5.09];
ybs[8297]=['',5.6873071,1.0362837,6.08];
ybs[8298]=['9 Peg',5.6966571,0.3046326,4.34];
ybs[8299]=['',5.69676,0.2596372,5.94];
ybs[8300]=['κ Peg',5.6970355,0.4494083,4.13];
ybs[8301]=['μ Cep',5.6906294,1.0277162,4.08];
ybs[8302]=['11 Cep',5.6820667,1.2464207,4.56];
ybs[8303]=['47 Cap',5.704912,-0.1600661,6];
ybs[8304]=['λ Cap',5.7061057,-0.1965421,5.58];
ybs[8305]=['',5.701532,0.6276508,6.4];
ybs[8306]=['12 Peg',5.7033415,0.4023599,5.29];
ybs[8307]=['δ Cap',5.7084135,-0.2796413,2.87];
ybs[8308]=['',5.714687,-0.8237651,5.58];
ybs[8309]=['',5.6868653,1.2640352,5.17];
ybs[8310]=['',5.704702,0.4479922,6.28];
ybs[8311]=['θ PsA',5.7118185,-0.5374421,5.01];
ybs[8312]=['',5.6963496,1.0919612,5.95];
ybs[8313]=['11 Peg',5.7088561,0.048714,5.64];
ybs[8314]=['',5.7036146,0.7533802,6.54];
ybs[8315]=['',5.7078662,0.3019263,6.21];
ybs[8316]=['',5.7234582,-1.1275979,5.62];
ybs[8317]=['',5.7107879,-0.1014403,6.17];
ybs[8318]=['ο Ind',5.7275276,-1.2134105,5.53];
ybs[8319]=['ν Cep',5.6989289,1.0685817,4.29];
ybs[8320]=['π2 Cyg',5.7055986,0.8624416,4.23];
ybs[8321]=['',5.7119802,0.6402878,6.47];
ybs[8322]=['',5.7198742,-0.2202138,6.31];
ybs[8323]=['',5.7134486,0.676384,6.12];
ybs[8324]=['12 Cep',5.7076056,1.0611211,5.52];
ybs[8325]=['',5.7222902,-0.2921478,6.38];
ybs[8326]=['',5.7181383,0.3589817,6.29];
ybs[8327]=['',5.7046424,1.2261924,6.29];
ybs[8328]=['14 Peg',5.7196259,0.5284842,5.04];
ybs[8329]=['13 Peg',5.7212527,0.3035372,5.29];
ybs[8330]=['',5.7185103,0.7200281,6.48];
ybs[8331]=['',5.7287765,-0.3231779,6.16];
ybs[8332]=['',5.7158504,1.0712538,6.17];
ybs[8333]=['',5.7274168,0.3478951,5.77];
ybs[8334]=['',5.7247436,0.6918968,6.17];
ybs[8335]=['',5.7305855,0.3731431,6.89];
ybs[8336]=['μ Cap',5.7356438,-0.2346576,5.08];
ybs[8337]=['',5.7457407,-1.0782472,5.9];
ybs[8338]=['γ Gru',5.7389925,-0.6502752,3.01];
ybs[8339]=['15 Peg',5.7312524,0.5043976,5.53];
ybs[8340]=['',5.7368995,-0.1781074,6.59];
ybs[8341]=['16 Peg',5.733793,0.4543386,5.08];
ybs[8342]=['',5.7280812,0.9756965,5.71];
ybs[8343]=['',5.7363858,0.3451421,5.68];
ybs[8344]=['',5.7381403,0.1216739,6.15];
ybs[8345]=['',5.7392785,-0.0727643,5.71];
ybs[8346]=['',5.7254806,1.1494561,6.37];
ybs[8347]=['π Ind',5.7499611,-1.0086569,6.19];
ybs[8348]=['',5.7411133,-0.0557453,6.2];
ybs[8349]=['',5.7392955,0.346018,6.39];
ybs[8350]=['',5.7475069,-0.5323055,6.41];
ybs[8351]=['',5.7496712,-0.6483191,5.46];
ybs[8352]=['',5.7525476,-0.6569262,6.18];
ybs[8353]=['δ Ind',5.7571084,-0.9579132,4.4];
ybs[8354]=['κ1 Ind',5.7599066,-1.0280677,6.12];
ybs[8355]=['',5.7776322,-1.3535588,6.41];
ybs[8356]=['13 Cep',5.7405861,0.9899255,5.8];
ybs[8357]=['',5.7484744,0.3725816,6.4];
ybs[8358]=['17 Peg',5.7510271,0.212654,5.54];
ybs[8359]=['',5.7422201,1.075982,6.13];
ybs[8360]=['',5.7426035,1.1419367,5.86];
ybs[8361]=['',5.7569625,-0.0927916,6.33];
ybs[8362]=['',5.7504193,0.8513087,6.42];
ybs[8363]=['',5.7595013,-0.3678189,6.12];
ybs[8364]=['',5.7624214,-0.6682308,5.5];
ybs[8365]=['',5.7822516,-1.3266077,5.95];
ybs[8366]=['',5.7679669,-0.9734396,6.01];
ybs[8367]=['',5.7599742,-0.0744332,6.22];
ybs[8368]=['',5.7477482,1.1123541,4.91];
ybs[8369]=['',5.7498087,1.1565232,6.43];
ybs[8370]=['18 Peg',5.7650636,0.1191392,6];
ybs[8371]=['η PsA',5.7688533,-0.4947088,5.42];
ybs[8372]=['ε Ind',5.780893,-0.9891919,4.69];
ybs[8373]=['',5.7576398,1.096182,5.93];
ybs[8374]=['',5.7601944,1.00822,6.59];
ybs[8375]=['28 Aqr',5.7693328,0.0124604,5.58];
ybs[8376]=['',5.7658608,0.5779631,6.46];
ybs[8377]=['20 Peg',5.7691216,0.2308834,5.6];
ybs[8378]=['19 Peg',5.769493,0.1460172,5.65];
ybs[8379]=['',5.7745548,-0.3105704,6.28];
ybs[8380]=['',5.7510798,1.3108223,6.35];
ybs[8381]=['29 Aqr',5.7756039,-0.2941728,6.37];
ybs[8382]=['',5.7732378,0.193436,6.37];
ybs[8383]=['',5.7795676,-0.5200144,7.1];
ybs[8384]=['',5.7653883,1.0925204,6.66];
ybs[8385]=['16 Cep',5.7577011,1.2791224,5.03];
ybs[8386]=['30 Aqr',5.7790211,-0.1119276,5.54];
ybs[8387]=['ο Aqr',5.7791215,-0.035705,4.69];
ybs[8388]=['',5.7712388,0.9248729,5.78];
ybs[8389]=['21 Peg',5.7788752,0.2006415,5.8];
ybs[8390]=['13 PsA',5.7844132,-0.5202273,6.47];
ybs[8391]=['14 Cep',5.7719589,1.0142057,5.56];
ybs[8392]=['',5.7764142,0.7811989,5.6];
ybs[8393]=['',5.7852771,-0.466223,5.96];
ybs[8394]=['κ2 Ind',5.7919042,-1.0389223,5.62];
ybs[8395]=['32 Aqr',5.7855341,-0.0139059,5.3];
ybs[8396]=['λ Gru',5.7921542,-0.6882365,4.46];
ybs[8397]=['',5.7838908,0.5768625,6.38];
ybs[8398]=['ν Peg',5.7893088,0.0902118,4.84];
ybs[8399]=['α Aqr',5.7898576,-0.0036572,2.96];
ybs[8400]=['',5.7867374,0.4674672,5.78];
ybs[8401]=['18 Cep',5.7794629,1.1035602,5.29];
ybs[8402]=['ξ Cep',5.7789167,1.1298802,4.29];
ybs[8403]=['ι Aqr',5.7929579,-0.2401463,4.27];
ybs[8404]=['23 Peg',5.78838,0.5074369,5.7];
ybs[8405]=['',5.8152807,-1.3224187,6.55];
ybs[8406]=['',5.7865142,0.8177694,6.13];
ybs[8407]=['',5.7890679,0.7892795,6.44];
ybs[8408]=['',5.7478609,1.4482321,6.98];
ybs[8409]=['',5.789902,0.7875738,5.14];
ybs[8410]=['α Gru',5.8016315,-0.8176915,1.74];
ybs[8411]=['20 Cep',5.7844273,1.0977331,5.27];
ybs[8412]=['',5.7889907,0.8437241,6.27];
ybs[8413]=['19 Cep',5.7850841,1.0889101,5.11];
ybs[8414]=['',5.7906512,0.7916615,6.19];
ybs[8415]=['ι Peg',5.7947232,0.444282,3.76];
ybs[8416]=['μ PsA',5.8018566,-0.5738249,4.5];
ybs[8417]=['',5.8206024,-1.3265247,6.15];
ybs[8418]=['υ PsA',5.8021014,-0.5922428,4.99];
ybs[8419]=['',5.7902178,0.9852959,6.39];
ybs[8420]=['',5.7968773,0.341843,5.75];
ybs[8421]=['',5.7970087,0.3160995,6.35];
ybs[8422]=['',5.8032909,-0.5762136,6.37];
ybs[8423]=['25 Peg',5.798411,0.3807169,5.78];
ybs[8424]=['35 Aqr',5.8041536,-0.3212926,5.81];
ybs[8425]=['',5.8092229,-0.8376872,6.43];
ybs[8426]=['',5.8002872,0.4477541,6.11];
ybs[8427]=['',5.794152,1.0288946,6.32];
ybs[8428]=['',5.795634,0.9323161,6.14];
ybs[8429]=['',5.8086154,-0.5917271,5.37];
ybs[8430]=['',5.7995151,0.8710444,6.42];
ybs[8431]=['',5.8087979,-0.4918553,6.44];
ybs[8432]=['τ PsA',5.8095301,-0.5661329,4.92];
ybs[8433]=['',5.801461,0.8002829,6.11];
ybs[8434]=['π1 Peg',5.8042144,0.5809025,5.58];
ybs[8435]=['θ Peg',5.8090142,0.1101142,3.53];
ybs[8436]=['',5.8098531,-0.0660227,6.27];
ybs[8437]=['38 Aqr',5.811182,-0.1999027,5.46];
ybs[8438]=['',5.8107832,-0.0725328,6.01];
ybs[8439]=['π2 Peg',5.8075346,0.5810126,4.29];
ybs[8440]=['',5.8092749,0.3443233,6.18];
ybs[8441]=['',5.8096034,0.257285,6.33];
ybs[8442]=['',5.813175,-0.3686305,6.09];
ybs[8443]=['',5.8107647,0.2048292,5.78];
ybs[8444]=['28 Peg',5.8100629,0.3680799,6.46];
ybs[8445]=['',5.8114136,0.5351966,6.32];
ybs[8446]=['',5.81606,0.2819102,5.95];
ybs[8447]=['39 Aqr',5.8191021,-0.2457777,6.03];
ybs[8448]=['',5.8121102,0.8889807,5.4];
ybs[8449]=['',5.8216429,-0.4575516,6.17];
ybs[8450]=['ζ Cep',5.8103674,1.0177457,3.35];
ybs[8451]=['',5.8170996,0.4374104,5.92];
ybs[8452]=['',5.8202434,-0.0804405,6.39];
ybs[8453]=['24 Cep',5.8042267,1.2645299,4.79];
ybs[8454]=['λ Cep',5.8131629,1.0389252,5.04];
ybs[8455]=['',5.8250263,-0.4375304,5.58];
ybs[8456]=['ψ Oct',5.8466622,-1.3508568,5.51];
ybs[8457]=['',5.8146442,0.9939844,5.24];
ybs[8458]=['',5.8062386,1.2605176,6.37];
ybs[8459]=['',5.8082995,1.2259911,5.5];
ybs[8460]=['',5.8197756,0.6059198,5.33];
ybs[8461]=['',5.8150964,1.0331723,6.3];
ybs[8462]=['',5.8293597,-0.7202842,6.23];
ybs[8463]=['λ PsA',5.8275801,-0.482664,5.43];
ybs[8464]=['',5.8153454,1.0624021,5.35];
ybs[8465]=['41 Aqr',5.8273856,-0.3658531,5.32];
ybs[8466]=['ε Oct',5.8573843,-1.4019517,5.1];
ybs[8467]=['',5.823616,0.5012668,5.89];
ybs[8468]=['',5.8166346,1.1065942,5.79];
ybs[8469]=['',5.8335568,-0.7738667,6.1];
ybs[8470]=['',5.8243646,0.6951156,4.49];
ybs[8471]=['μ1 Gru',5.833591,-0.7196693,4.79];
ybs[8472]=['',5.8239365,0.79505,5.53];
ybs[8473]=['μ2 Gru',5.8372074,-0.7245673,5.1];
ybs[8474]=['',5.8280333,0.7516485,5.71];
ybs[8475]=['',5.8230116,1.1043509,6.11];
ybs[8476]=['',5.8342726,0.1511833,6.21];
ybs[8477]=['',5.8376102,-0.450041,6.15];
ybs[8478]=['',5.8175343,1.281405,6.08];
ybs[8479]=['ε Cep',5.8287247,0.9975614,4.19];
ybs[8480]=['',5.8369004,-0.0258925,6.15];
ybs[8481]=['42 Aqr',5.8381439,-0.2219792,5.34];
ybs[8482]=['',5.8391828,-0.4018975,6.17];
ybs[8483]=['1 Lac',5.8335601,0.6608093,4.13];
ybs[8484]=['θ Aqr',5.8381959,-0.1338739,4.16];
ybs[8485]=['',5.8384061,-0.1558067,5.79];
ybs[8486]=['',5.8455563,-0.934004,5.37];
ybs[8487]=['α Tuc',5.8469716,-1.0497522,2.86];
ybs[8488]=['',5.8360892,0.4872435,6.37];
ybs[8489]=['44 Aqr',5.8393555,-0.0920528,5.75];
ybs[8490]=['υ Oct',5.9138673,-1.4946687,5.77];
ybs[8491]=['',5.8348784,1.0006507,5.88];
ybs[8492]=['',5.8434738,-0.0021741,6.39];
ybs[8493]=['45 Aqr',5.8477978,-0.2302362,5.95];
ybs[8494]=['',5.8559467,-1.0017521,6.34];
ybs[8495]=['',5.8465238,0.6611803,6.17];
ybs[8496]=['25 Cep',5.8422552,1.09812,5.75];
ybs[8497]=['ρ Aqr',5.8528752,-0.1345195,5.37];
ybs[8498]=['30 Peg',5.8537981,0.1030304,5.37];
ybs[8499]=['',5.8558106,0.1448717,6.17];
ybs[8500]=['ν Ind',5.8747427,-1.2590943,5.29];
ybs[8501]=['47 Aqr',5.8591941,-0.3749718,5.13];
ybs[8502]=['',5.8557914,0.4720968,6.47];
ybs[8503]=['γ Aqr',5.8591343,-0.0222213,3.84];
ybs[8504]=['',5.8536339,0.8917693,6.42];
ybs[8505]=['31 Peg',5.8583128,0.215012,5.01];
ybs[8506]=['π1 Gru',5.8647264,-0.799945,6.62];
ybs[8507]=['32 Peg',5.8571609,0.4964503,4.81];
ybs[8508]=['2 Lac',5.8554013,0.8142055,4.57];
ybs[8509]=['π2 Gru',5.8664748,-0.799609,5.62];
ybs[8510]=['',5.8407246,1.3369436,6.66];
ybs[8511]=['',5.8806282,-1.3072605,6.04];
ybs[8512]=['',5.8769089,-1.2272594,5.78];
ybs[8513]=['',5.859132,0.7363962,6.41];
ybs[8514]=['49 Aqr',5.8676304,-0.4301894,5.53];
ybs[8515]=['',5.8674202,-0.123569,5.93];
ybs[8516]=['',5.8748252,-1.0067481,5.32];
ybs[8517]=['33 Peg',5.8675155,0.3658701,6.04];
ybs[8518]=['51 Aqr',5.8699137,-0.0824206,5.78];
ybs[8519]=['50 Aqr',5.8715197,-0.2341321,5.76];
ybs[8520]=['',5.863594,1.0017972,6.16];
ybs[8521]=['',5.8682084,0.6752353,6.22];
ybs[8522]=['',5.8632724,1.0914293,6.04];
ybs[8523]=['β Lac',5.8662682,0.9135682,4.43];
ybs[8524]=['π Aqr',5.8748867,0.0260463,4.66];
ybs[8525]=['δ Tuc',5.8857045,-1.1318645,4.48];
ybs[8526]=['4 Lac',5.87056,0.8655269,4.57];
ybs[8527]=['',5.8792211,-0.4113296,6.29];
ybs[8528]=['',5.8763701,0.3239221,6.26];
ybs[8529]=['53 Aqr',5.8808099,-0.290183,6.57];
ybs[8530]=['53 Aqr',5.8808245,-0.2902024,6.35];
ybs[8531]=['',5.8073173,1.501347,5.27];
ybs[8532]=['',5.8915579,-1.1758904,5.55];
ybs[8533]=['34 Peg',5.8807174,0.0786924,5.75];
ybs[8534]=['',5.8807276,0.6555288,6.46];
ybs[8535]=['',5.8636744,1.3676002,6.76];
ybs[8536]=['35 Peg',5.8861024,0.0839669,4.79];
ybs[8537]=['ν Gru',5.8903263,-0.680964,5.47];
ybs[8538]=['',5.8836412,0.6968229,6.14];
ybs[8539]=['',5.8810719,0.9869578,6.57];
ybs[8540]=['',5.8852576,0.5577312,5.98];
ybs[8541]=['δ1 Gru',5.8931255,-0.7571211,3.97];
ybs[8542]=['',5.8756231,1.2371902,5.47];
ybs[8543]=['ζ1 Aqr',5.8904051,0.0016638,4.59];
ybs[8544]=['ζ2 Aqr',5.8904342,0.0016686,4.42];
ybs[8545]=['δ2 Gru',5.8952592,-0.7615505,4.11];
ybs[8546]=['26 Cep',5.8809551,1.1387821,5.46];
ybs[8547]=['36 Peg',5.8916007,0.1613479,5.58];
ybs[8548]=['',5.8949253,-0.471089,5.95];
ybs[8549]=['',5.8914721,0.4691221,5.79];
ybs[8550]=['',5.8958167,-0.2233872,6.4];
ybs[8551]=['37 Peg',5.8953072,0.079369,5.48];
ybs[8552]=['56 Aqr',5.8969981,-0.2525477,6.37];
ybs[8553]=['',5.8864741,1.120519,6.29];
ybs[8554]=['',5.8937456,0.6255493,6.56];
ybs[8555]=['ζ PsA',5.8998249,-0.4530451,6.43];
ybs[8556]=['δ Cep',5.8905339,1.0215572,3.75];
ybs[8557]=['5 Lac',5.8925457,0.8346631,4.36];
ybs[8558]=['σ Aqr',5.8984944,-0.1843429,4.82];
ybs[8559]=['38 Peg',5.8951162,0.5705192,5.65];
ybs[8560]=['',5.8950199,0.8634485,6.4];
ybs[8561]=['β PsA',5.9026014,-0.5625185,4.29];
ybs[8562]=['',5.9231594,-1.3727823,6.15];
ybs[8563]=['ρ1 Cep',5.8767698,1.3770801,5.83];
ybs[8564]=['6 Lac',5.8968717,0.7546676,4.51];
ybs[8565]=['',5.9012602,-0.0487818,6.16];
ybs[8566]=['',5.9013127,-0.1123796,6.14];
ybs[8567]=['ν Tuc',5.9101215,-1.0797606,4.81];
ybs[8568]=['58 Aqr',5.9030412,-0.1883098,6.38];
ybs[8569]=['',5.9019035,0.5176461,6.35];
ybs[8570]=['α Lac',5.9001634,0.8796214,3.77];
ybs[8571]=['39 Peg',5.9065272,0.3551112,6.42];
ybs[8572]=['',5.907425,0.2788992,6.32];
ybs[8573]=['',5.9054881,0.6963175,5.88];
ybs[8574]=['',5.9044898,0.9451621,6.35];
ybs[8575]=['60 Aqr',5.9132032,-0.0254381,5.89];
ybs[8576]=['ρ2 Cep',5.890743,1.3777609,5.5];
ybs[8577]=['υ Aqr',5.9162941,-0.3593899,5.2];
ybs[8578]=['',5.9224184,-1.0082165,6.23];
ybs[8579]=['',5.9103429,0.9903272,5.71];
ybs[8580]=['',5.9066054,1.2222546,6.6];
ybs[8581]=['',5.9203358,-0.4166821,5.97];
ybs[8582]=['η Aqr',5.9188916,-0.00001,4.02];
ybs[8583]=['',5.90758,1.2302888,6.34];
ybs[8584]=['',5.9020839,1.3324302,5.68];
ybs[8585]=['σ1 Gru',5.9244897,-0.7062582,6.28];
ybs[8586]=['',5.9247476,-0.5505939,5.82];
ybs[8587]=['σ2 Gru',5.9266314,-0.7064019,5.86];
ybs[8588]=['8 Lac',5.9204904,0.693789,5.73];
ybs[8589]=['',5.9217149,0.6229829,6.1];
ybs[8590]=['',5.924184,0.2061952,6.4];
ybs[8591]=['',5.9202837,0.875948,6.29];
ybs[8592]=['',5.9199384,0.9806482,6.38];
ybs[8593]=['',5.9262307,0.2215605,6.3];
ybs[8594]=['',5.9246847,0.624299,6.3];
ybs[8595]=['κ Aqr',5.9294184,-0.0717446,5.03];
ybs[8596]=['',5.9363799,-0.9175989,6.65];
ybs[8597]=['',5.9321432,-0.1357912,6.23];
ybs[8598]=['9 Lac',5.9267223,0.9016821,4.63];
ybs[8599]=['',5.9340875,-0.4996911,6.47];
ybs[8600]=['31 Cep',5.9180086,1.287355,5.08];
ybs[8601]=['',5.9346587,-0.5753265,5.66];
ybs[8602]=['',5.930932,0.7906435,6.4];
ybs[8603]=['40 Peg',5.9339819,0.3427796,5.82];
ybs[8604]=['',5.938388,-0.4923138,6.31];
ybs[8605]=['',5.9438536,-1.0001475,5.97];
ybs[8606]=['',5.932023,0.9933257,5.21];
ybs[8607]=['10 Lac',5.9353201,0.6836096,4.88];
ybs[8608]=['',5.9412098,-0.533041,5.87];
ybs[8609]=['41 Peg',5.9379386,0.3455557,6.21];
ybs[8610]=['',5.9240421,1.3175297,5.79];
ybs[8611]=['',5.9367059,0.6581724,6.03];
ybs[8612]=['30 Cep',5.9317344,1.1118092,5.19];
ybs[8613]=['ε PsA',5.9423913,-0.4699416,4.17];
ybs[8614]=['',5.9426867,-0.059973,6.31];
ybs[8615]=['β Oct',5.9700484,-1.4183013,4.15];
ybs[8616]=['',5.9427921,0.2559947,5.71];
ybs[8617]=['11 Lac',5.9406766,0.7748265,4.46];
ybs[8618]=['',5.9394731,0.9418488,5.93];
ybs[8619]=['ζ Peg',5.9453886,0.1911043,3.4];
ybs[8620]=['',5.9513356,-0.8219147,5.98];
ybs[8621]=['β Gru',5.9515591,-0.8162277,2.1];
ybs[8622]=['19 PsA',5.9498896,-0.5103792,6.17];
ybs[8623]=['',5.9453594,0.5425167,6.34];
ybs[8624]=['',5.9517095,-0.7702042,6.07];
ybs[8625]=['12 Lac',5.9449756,0.7041292,5.25];
ybs[8626]=['ο Peg',5.9464069,0.5135742,4.79];
ybs[8627]=['',5.9475002,0.2554213,5.9];
ybs[8628]=['',5.9454924,0.7272358,5.94];
ybs[8629]=['ρ Gru',5.9550482,-0.7207508,4.85];
ybs[8630]=['',5.9525994,-0.1429999,6.45];
ybs[8631]=['',5.9590436,-1.0538443,6.3];
ybs[8632]=['67 Aqr',5.9533675,-0.1194568,6.41];
ybs[8633]=['',5.9483972,0.942951,6.12];
ybs[8634]=['66 Aqr',5.9550495,-0.3265826,4.69];
ybs[8635]=['η Peg',5.9518304,0.5295284,2.94];
ybs[8636]=['',5.9513566,0.6618483,6.43];
ybs[8637]=['',5.9517948,0.8253134,6.39];
ybs[8638]=['',5.9552086,0.1929925,6.51];
ybs[8639]=['',5.9563998,0.6908729,5.95];
ybs[8640]=['η Gru',5.9646556,-0.9316818,4.85];
ybs[8641]=['13 Lac',5.9563711,0.7319511,5.08];
ybs[8642]=['',5.9646551,-0.8103328,5.51];
ybs[8643]=['',5.9667004,-0.8527672,6.62];
ybs[8644]=['',5.9681929,-0.8651047,6.48];
ybs[8645]=['45 Peg',5.9627739,0.3400854,6.25];
ybs[8646]=['',5.9592675,0.9186695,6.55];
ybs[8647]=['',5.9692364,-0.8171705,6.56];
ybs[8648]=['ξ Oct',5.988131,-1.3963415,5.35];
ybs[8649]=['',5.984218,-1.3426992,6.73];
ybs[8650]=['ξ Peg',5.9682084,0.2145321,4.19];
ybs[8651]=['',5.9653947,0.7795516,5.76];
ybs[8652]=['λ Peg',5.9673517,0.413373,3.95];
ybs[8653]=['',5.9715412,-0.5941448,6.28];
ybs[8654]=['',5.9768467,-1.0745094,6.37];
ybs[8655]=['68 Aqr',5.972341,-0.3402375,5.26];
ybs[8656]=['',5.9736475,-0.6650183,6.71];
ybs[8657]=['',5.9815194,-1.225715,6.34];
ybs[8658]=['τ1 Aqr',5.9729746,-0.2432501,5.66];
ybs[8659]=['',5.9741064,-0.4501678,6.3];
ybs[8660]=['ε Gru',5.9773075,-0.8935668,3.49];
ybs[8661]=['70 Aqr',5.9763773,-0.1821467,6.19];
ybs[8662]=['',5.9702749,1.0227958,6.36];
ybs[8663]=['',5.9743373,0.6551253,5.9];
ybs[8664]=['τ2 Aqr',5.9811606,-0.2351483,4.01];
ybs[8665]=['',5.9831382,-0.5704733,6.33];
ybs[8666]=['',5.9806479,0.1849765,6.54];
ybs[8667]=['',5.9766038,0.9518039,6.12];
ybs[8668]=['',5.9759671,1.1005638,6.06];
ybs[8669]=['μ Peg',5.9824994,0.4314667,3.48];
ybs[8670]=['',5.9878296,-0.6813278,5.42];
ybs[8671]=['',5.9915051,-1.0430354,6.46];
ybs[8672]=['',5.9767778,1.1988604,6.19];
ybs[8673]=['',5.9808107,0.9777732,5.43];
ybs[8674]=['',5.9934766,-1.1007561,6.12];
ybs[8675]=['14 Lac',5.9837847,0.7343162,5.92];
ybs[8676]=['',5.9854027,0.3361591,6.4];
ybs[8677]=['',5.9827293,0.8865664,6.21];
ybs[8678]=['21 PsA',5.9890291,-0.5134118,5.97];
ybs[8679]=['ι Cep',5.9798951,1.157503,3.52];
ybs[8680]=['γ PsA',5.9942145,-0.5716928,4.46];
ybs[8681]=['',5.9876128,1.078905,5.6];
ybs[8682]=['σ Peg',5.9931537,0.1737562,5.16];
ybs[8683]=['λ Aqr',5.9942775,-0.1301971,3.74];
ybs[8684]=['15 Lac',5.9910573,0.7580379,4.94];
ybs[8685]=['τ1 Gru',5.9993303,-0.8460991,6.04];
ybs[8686]=['ρ Ind',6.0047743,-1.2209152,6.05];
ybs[8687]=['',5.9660811,1.4533866,4.74];
ybs[8688]=['',5.9958518,0.2960279,5.64];
ybs[8689]=['74 Aqr',5.9980884,-0.2006527,5.8];
ybs[8690]=['',5.9945549,0.8819489,6.46];
ybs[8691]=['',5.9961686,0.7031457,6.34];
ybs[8692]=['',5.9950605,1.0510572,6.01];
ybs[8693]=['',5.9981753,0.783117,5.81];
ybs[8694]=['δ Aqr',6.0032459,-0.274026,3.27];
ybs[8695]=['78 Aqr',6.002794,-0.1236468,6.19];
ybs[8696]=['77 Aqr',6.0037235,-0.2818991,5.56];
ybs[8697]=['',6.0002225,0.7068085,5.81];
ybs[8698]=['',6.0061286,-0.6329999,6.4];
ybs[8699]=['',6.0026473,0.2977872,6.12];
ybs[8700]=['1 Psc',6.0045527,0.0206833,6.11];
ybs[8701]=['',6.0054531,-0.0849522,5.72];
ybs[8702]=['ρ Peg',6.0055029,0.1559663,4.9];
ybs[8703]=['',6.0043308,0.6492152,5.91];
ybs[8704]=['',6.0087066,-0.5499982,6.1];
ybs[8705]=['δ PsA',6.0091192,-0.5658222,4.21];
ybs[8706]=['',6.0110734,-0.5488187,6.48];
ybs[8707]=['τ3 Gru',6.0130929,-0.8351146,5.7];
ybs[8708]=['',6.0073945,0.6365586,5.74];
ybs[8709]=['',6.0125852,0.2088977,6.51];
ybs[8710]=['16 Lac',6.0101518,0.7282287,5.59];
ybs[8711]=['',6.0101523,0.8701193,4.95];
ybs[8712]=['',6.0146272,-0.0818439,6.31];
ybs[8713]=['α PsA',6.0165011,-0.5148979,1.16];
ybs[8714]=['51 Peg',6.0151344,0.3645923,5.49];
ybs[8715]=['',6.0156735,0.0686089,6.28];
ybs[8716]=['',6.0129844,0.8518047,5.43];
ybs[8717]=['',6.0206445,-0.6178845,6.13];
ybs[8718]=['',6.0158096,0.6881768,6.18];
ybs[8719]=['',6.0188417,-0.0396966,6.16];
ybs[8720]=['',6.0194281,-0.0225048,6.37];
ybs[8721]=['',5.9791541,1.487276,5.9];
ybs[8722]=['',6.0201484,0.1654192,6.43];
ybs[8723]=['',6.0207146,0.1302123,6.33];
ybs[8724]=['52 Peg',6.0227939,0.2068189,5.75];
ybs[8725]=['',6.0249755,-0.5121005,5.51];
ybs[8726]=['',6.0247809,-0.2260169,6.07];
ybs[8727]=['2 Psc',6.0240368,0.0189154,5.43];
ybs[8728]=['',6.0271015,-0.4370841,5.65];
ybs[8729]=['',6.0220561,0.9211043,6.29];
ybs[8730]=['',6.0217275,1.0460746,6.43];
ybs[8731]=['',6.0284728,-0.4451554,6.29];
ybs[8732]=['ζ Gru',6.0309813,-0.9186184,4.12];
ybs[8733]=['',5.9957263,1.4742144,4.71];
ybs[8734]=['',6.0320023,-0.8871292,5.68];
ybs[8735]=['3 Psc',6.0291714,0.0053581,6.21];
ybs[8736]=['',6.0295077,0.0546784,5.83];
ybs[8737]=['',6.025925,0.9959957,5];
ybs[8738]=['',6.0291734,0.5446165,6.6];
ybs[8739]=['',6.032493,-0.501474,5.55];
ybs[8740]=['',6.028355,0.7940575,6.5];
ybs[8741]=['',6.0326829,-0.3956585,6.28];
ybs[8742]=['81 Aqr',6.0325645,-0.121123,6.21];
ybs[8743]=['',6.0299517,0.6776983,6.54];
ybs[8744]=['',6.0331307,-0.0801123,5.94];
ybs[8745]=['',6.0380117,-0.6335441,6.47];
ybs[8746]=['',6.0321651,0.9987966,6.2];
ybs[8747]=['ο And',6.0343066,0.7408477,3.62];
ybs[8748]=['82 Aqr',6.0375763,-0.1126216,6.15];
ybs[8749]=['',6.0385693,-0.3621401,5.97];
ybs[8750]=['',6.0372177,0.5567945,6.57];
ybs[8751]=['2 And',6.0372865,0.7483833,5.1];
ybs[8752]=['π PsA',6.0420401,-0.6043706,5.11];
ybs[8753]=['',6.0379119,0.7710923,6.39];
ybs[8754]=['',6.0489897,-1.1990154,5.52];
ybs[8755]=['',6.0375553,0.9661764,6.5];
ybs[8756]=['',6.0442961,-0.7218108,5.79];
ybs[8757]=['',6.043718,-0.0815709,6.68];
ybs[8758]=['β Psc',6.0432972,0.0687939,4.53];
ybs[8759]=['κ Gru',6.0474594,-0.9397426,5.37];
ybs[8760]=['β Peg',6.0426079,0.492259,2.42];
ybs[8761]=['',6.0438823,0.1176053,6.41];
ybs[8762]=['',6.0402938,1.0570901,6.74];
ybs[8763]=['',6.0402086,1.0242681,6.43];
ybs[8764]=['',6.0406372,1.1751425,5.24];
ybs[8765]=['3 And',6.0440457,0.875699,4.65];
ybs[8766]=['α Peg',6.0470533,0.2675065,2.49];
ybs[8767]=['83 Aqr',6.0490156,-0.1321535,5.43];
ybs[8768]=['',6.0493219,-0.2959623,6.14];
ybs[8769]=['',6.0485391,0.291205,6.44];
ybs[8770]=['',6.0495019,0.0249361,6.39];
ybs[8771]=['',6.0656176,-1.3850693,6.12];
ybs[8772]=['θ Gru',6.0569077,-0.7574478,4.28];
ybs[8773]=['',6.0537523,0.3253191,6.13];
ybs[8774]=['86 Aqr',6.0557906,-0.4122658,4.47];
ybs[8775]=['υ Gru',6.056895,-0.6766681,5.61];
ybs[8776]=['',6.0582307,-0.8636699,6.33];
ybs[8777]=['',6.0547356,0.3496378,6.3];
ybs[8778]=['',6.0586304,-0.8825144,5.83];
ybs[8779]=['',6.0655407,-1.2821917,6.15];
ybs[8780]=['55 Peg',6.0569,0.1663551,4.52];
ybs[8781]=['56 Peg',6.057219,0.4466358,4.76];
ybs[8782]=['1 Cas',6.0544325,1.0391981,4.85];
ybs[8783]=['',6.0586573,0.5750491,6.02];
ybs[8784]=['',6.0588574,0.3709911,5.99];
ybs[8785]=['',6.0577639,0.806169,6.66];
ybs[8786]=['',6.0570416,0.9239494,6.11];
ybs[8787]=['',6.0631165,-0.5009298,5.6];
ybs[8788]=['',6.0568671,1.044571,6.4];
ybs[8789]=['4 And',6.0593033,0.8117404,5.33];
ybs[8790]=['5 And',6.0596925,0.8625054,5.7];
ybs[8791]=['',6.06174,0.7798796,6.56];
ybs[8792]=['5 Psc',6.0642802,0.0392696,5.4];
ybs[8793]=['',6.0594397,1.112742,6.26];
ybs[8794]=['',6.0719861,-1.1647472,6.47];
ybs[8795]=['',6.0824094,-1.4100538,6.41];
ybs[8796]=['',6.0601075,1.1230252,6.21];
ybs[8797]=['88 Aqr',6.0678191,-0.3673953,3.66];
ybs[8798]=['',6.0691811,-0.4881035,5.87];
ybs[8799]=['',6.0702877,-0.7459222,5.81];
ybs[8800]=['57 Peg',6.0679057,0.1535807,5.12];
ybs[8801]=['',6.0694155,-0.2511217,6.42];
ybs[8802]=['89 Aqr',6.0698658,-0.3898218,4.69];
ybs[8803]=['',6.0711619,-0.7063222,5.83];
ybs[8804]=['π Cep',6.0589953,1.3178909,4.41];
ybs[8805]=['ι Gru',6.0720886,-0.7875668,3.9];
ybs[8806]=['58 Peg',6.0700785,0.1735609,5.39];
ybs[8807]=['2 Cas',6.0680982,1.0376921,5.7];
ybs[8808]=['',6.0736998,-0.5131712,6.51];
ybs[8809]=['',6.073002,0.309218,5.71];
ybs[8810]=['6 And',6.0715869,0.7621255,5.94];
ybs[8811]=['59 Peg',6.0775588,0.1543318,5.16];
ybs[8812]=['60 Peg',6.0777632,0.4707116,6.17];
ybs[8813]=['',6.0847444,-0.8638709,6.8];
ybs[8814]=['',6.0888123,-1.0921776,6.12];
ybs[8815]=['7 And',6.0806574,0.8644448,4.52];
ybs[8816]=['',6.0831831,0.5159956,6.35];
ybs[8817]=['',6.0836969,0.9999176,5.56];
ybs[8818]=['',6.0849823,0.1952631,5.82];
ybs[8819]=['φ Aqr',6.0889599,-0.1034289,4.22];
ybs[8820]=['',6.0921322,-0.715282,5.77];
ybs[8821]=['',6.0905072,-0.1844068,6.12];
ybs[8822]=['',6.0880206,0.8855907,6.31];
ybs[8823]=['',6.0888372,0.5217577,6.41];
ybs[8824]=['',6.0899678,0.4228223,6.36];
ybs[8825]=['',6.0943868,-0.0588771,5.55];
ybs[8826]=['ψ1 Aqr',6.0958237,-0.1564646,4.21];
ybs[8827]=['61 Peg',6.0950115,0.4951634,6.49];
ybs[8828]=['',6.1011809,-1.0799745,5.66];
ybs[8829]=['',6.0887009,1.2977216,5.84];
ybs[8830]=['',6.0958873,0.4344844,6.6];
ybs[8831]=['',6.099515,-0.774334,5.92];
ybs[8832]=['',6.1002048,-0.71683,6.47];
ybs[8833]=['γ Tuc',6.1031111,-1.0142572,3.99];
ybs[8834]=['',6.1119434,-1.3849086,6.33];
ybs[8835]=['χ Aqr',6.0999867,-0.1327071,5.06];
ybs[8836]=['',6.0934102,1.2393762,5.56];
ybs[8837]=['γ Psc',6.1012907,0.0594348,3.69];
ybs[8838]=['',6.0987631,0.930901,5.54];
ybs[8839]=['',6.0974109,1.0836071,6.53];
ybs[8840]=['',6.1073466,-1.1754416,6.13];
ybs[8841]=['',6.1035847,-0.2022813,6.34];
ybs[8842]=['',6.1013979,0.7904127,6.43];
ybs[8843]=['ψ2 Aqr',6.1045989,-0.1581144,4.39];
ybs[8844]=['φ Gru',6.1060146,-0.71037,5.53];
ybs[8845]=['8 And',6.1033882,0.8576281,4.85];
ybs[8846]=['',6.1042721,0.7960814,6.48];
ybs[8847]=['τ Oct',6.1553234,-1.5232141,5.49];
ybs[8848]=['γ Scl',6.1087949,-0.5656374,4.41];
ybs[8849]=['9 And',6.1063103,0.7312383,6.02];
ybs[8850]=['ψ3 Aqr',6.1092181,-0.1655884,4.98];
ybs[8851]=['94 Aqr',6.1098991,-0.2327494,5.08];
ybs[8852]=['',6.1004332,1.3163676,6.38];
ybs[8853]=['96 Aqr',6.1110998,-0.0872854,5.55];
ybs[8854]=['',6.1111969,-0.3133201,5.93];
ybs[8855]=['',6.109109,0.7899454,6.5];
ybs[8856]=['',6.1127131,-0.586163,6.37];
ybs[8857]=['ο Cep',6.1067153,1.1909243,4.75];
ybs[8858]=['',6.111049,0.6094112,6.32];
ybs[8859]=['11 And',6.1110562,0.8508243,5.44];
ybs[8860]=['',6.1119201,0.8465582,6.32];
ybs[8861]=['10 And',6.1127914,0.7365543,5.79];
ybs[8862]=['',6.1177601,-0.8758616,6.05];
ybs[8863]=['7 Psc',6.1151469,0.0960775,5.05];
ybs[8864]=['',6.1166961,-0.1009599,6.17];
ybs[8865]=['τ Peg',6.1162996,0.416501,4.6];
ybs[8866]=['',6.1140221,1.0837347,6.45];
ybs[8867]=['63 Peg',6.1170753,0.5329972,5.59];
ybs[8868]=['',6.1193586,-0.4688502,5.64];
ybs[8869]=['',6.1165266,0.7721313,6.13];
ybs[8870]=['12 And',6.1172716,0.6685609,5.77];
ybs[8871]=['',6.1154972,1.0879775,6.39];
ybs[8872]=['64 Peg',6.1218182,0.5573899,5.32];
ybs[8873]=['',6.1221011,0.4665698,6.62];
ybs[8874]=['',6.1271475,-1.0460133,6.09];
ybs[8875]=['97 Aqr',6.1253548,-0.2603248,5.2];
ybs[8876]=['65 Peg',6.1252228,0.365686,6.29];
ybs[8877]=['98 Aqr',6.1267695,-0.3486622,3.97];
ybs[8878]=['66 Peg',6.1270329,0.2170767,5.08];
ybs[8879]=['',6.1241615,1.0516874,5.56];
ybs[8880]=['',6.1311824,-0.9369725,6.15];
ybs[8881]=['',6.130392,-0.7505036,6.1];
ybs[8882]=['',6.1290918,0.0072452,6.31];
ybs[8883]=['',6.132515,-0.903515,5.75];
ybs[8884]=['',6.1300117,0.5699397,6.69];
ybs[8885]=['',6.1318176,-0.323998,6.19];
ybs[8886]=['',6.1374057,-0.9900429,5.59];
ybs[8887]=['',6.1333928,0.7197143,6.72];
ybs[8888]=['67 Peg',6.1346227,0.5673863,5.57];
ybs[8889]=['4 Cas',6.1341668,1.0892008,4.98];
ybs[8890]=['υ Peg',6.1370183,0.410642,4.4];
ybs[8891]=['99 Aqr',6.1401867,-0.3581068,4.39];
ybs[8892]=['ο Gru',6.1429339,-0.9180027,5.52];
ybs[8893]=['',6.1454562,-1.1598949,6.45];
ybs[8894]=['',6.1458149,-1.0184359,5.63];
ybs[8895]=['',6.1452556,-0.873244,6.2];
ybs[8896]=['κ Psc',6.1439294,0.0240779,4.94];
ybs[8897]=['9 Psc',6.1452973,0.021756,6.25];
ybs[8898]=['13 And',6.1444791,0.7511192,5.75];
ybs[8899]=['',6.1488544,-0.6182018,6.32];
ybs[8900]=['69 Peg',6.147023,0.4414162,5.98];
ybs[8901]=['θ Psc',6.1484188,0.1134983,4.28];
ybs[8902]=['',6.1490288,-0.1976695,6.37];
ybs[8903]=['',6.1445725,1.2301735,5.6];
ybs[8904]=['',6.1535943,-1.0993247,5.68];
ybs[8905]=['',6.1532955,-0.7744656,6.43];
ybs[8906]=['',6.1530445,-0.1595571,6.18];
ybs[8907]=['',6.1532381,0.4044267,6.35];
ybs[8908]=['70 Peg',6.1535648,0.2248809,4.55];
ybs[8909]=['',6.1553107,-0.0769442,6.25];
ybs[8910]=['',6.1575323,0.859702,6.17];
ybs[8911]=['',6.1570038,1.0240391,4.91];
ybs[8912]=['',6.1599848,0.6769472,6.05];
ybs[8913]=['',6.1617901,-0.1075826,6.39];
ybs[8914]=['',6.163908,-0.7804987,6.02];
ybs[8915]=['14 And',6.1627267,0.686974,5.22];
ybs[8916]=['',6.1639901,-0.0691654,6.49];
ybs[8917]=['100 Aqr',6.1648451,-0.3707969,6.29];
ybs[8918]=['',6.164672,0.4979068,6.41];
ybs[8919]=['13 Psc',6.1658741,-0.0167808,6.38];
ybs[8920]=['',6.1729245,-1.3484557,5.81];
ybs[8921]=['',6.1676527,0.6122073,6.65];
ybs[8922]=['β Scl',6.1704786,-0.6578827,4.37];
ybs[8923]=['',6.1374168,1.5242821,5.58];
ybs[8924]=['101 Aqr',6.1717093,-0.3628538,4.71];
ybs[8925]=['71 Peg',6.1723481,0.3948519,5.32];
ybs[8926]=['',6.1732665,0.7885839,6.24];
ybs[8927]=['',6.1743507,0.3659139,6.06];
ybs[8928]=['72 Peg',6.174419,0.548902,4.98];
ybs[8929]=['14 Psc',6.1754301,-0.0196,5.87];
ybs[8930]=['',6.1805425,-1.1268697,7.4];
ybs[8931]=['',6.1784267,-0.2639163,5.96];
ybs[8932]=['15 And',6.1772973,0.7044309,5.59];
ybs[8933]=['73 Peg',6.177392,0.5868103,5.63];
ybs[8934]=['ι Phe',6.1796811,-0.7415981,4.71];
ybs[8935]=['',6.1779827,0.6658157,6.18];
ybs[8936]=['',6.1814979,-0.1281047,6.39];
ybs[8937]=['',6.1783551,1.2525665,5.84];
ybs[8938]=['',6.1830966,0.4308471,6.45];
ybs[8939]=['16 Psc',6.1851834,0.0388659,5.68];
ybs[8940]=['',6.1855737,0.5764613,6.35];
ybs[8941]=['',6.1883852,-0.5540751,6.52];
ybs[8942]=['',6.1948269,-1.3394574,6];
ybs[8943]=['',6.1907883,-0.2257684,5.65];
ybs[8944]=['',6.1917778,-0.7918173,4.74];
ybs[8945]=['74 Peg',6.190692,0.2958377,6.26];
ybs[8946]=['λ And',6.1901061,0.8130223,3.82];
ybs[8947]=['',6.1899827,0.7776115,5.8];
ybs[8948]=['75 Peg',6.1919229,0.323327,5.53];
ybs[8949]=['',6.1919084,0.808514,6.58];
ybs[8950]=['ι And',6.1926288,0.7573468,4.29];
ybs[8951]=['θ Phe',6.1988241,-0.8118048,6.09];
ybs[8952]=['18 And',6.1969558,0.8830745,5.3];
ybs[8953]=['ω1 Aqr',6.2000596,-0.2460367,5];
ybs[8954]=['ι Psc',6.2007139,0.1003774,4.13];
ybs[8955]=['',6.2005623,0.1710777,5.97];
ybs[8956]=['',6.1966022,1.3162846,5.95];
ybs[8957]=['',6.1974494,1.2937699,5.98];
ybs[8958]=['',6.2010113,0.6593385,6.53];
ybs[8959]=['γ Cep',6.1972185,1.3571205,3.21];
ybs[8960]=['μ Scl',6.2038404,-0.5576015,5.31];
ybs[8961]=['κ And',6.2025546,0.775951,4.14];
ybs[8962]=['',6.2037695,0.6430783,6.23];
ybs[8963]=['',6.2059009,-0.4194971,6.6];
ybs[8964]=['',6.2059969,-0.2016849,5.89];
ybs[8965]=['103 Aqr',6.207878,-0.3124548,5.34];
ybs[8966]=['',6.2070688,0.8663308,6.26];
ybs[8967]=['104 Aqr',6.2086984,-0.3087748,4.82];
ybs[8968]=['',6.2094159,0.128726,5.89];
ybs[8969]=['λ Psc',6.2098755,0.0332469,4.5];
ybs[8970]=['',6.2090166,1.0015554,6.24];
ybs[8971]=['',6.2105893,0.7874377,6.57];
ybs[8972]=['',6.2117506,-0.2674343,5.28];
ybs[8973]=['ω2 Aqr',6.2128672,-0.2516776,4.49];
ybs[8974]=['',6.2108466,1.1281891,6.56];
ybs[8975]=['',6.2116677,1.0786897,6.4];
ybs[8976]=['77 Peg',6.2156415,0.1824977,5.06];
ybs[8977]=['',6.2176808,-0.2645827,6.36];
ybs[8978]=['',6.2186524,-0.7846713,6.09];
ybs[8979]=['',6.220652,-1.2281058,6.07];
ybs[8980]=['',6.2220651,-1.3729874,5.75];
ybs[8981]=['',6.2195779,-1.1218881,5.72];
ybs[8982]=['78 Peg',6.2182845,0.5146391,4.93];
ybs[8983]=['106 Aqr',6.2193316,-0.3168115,5.24];
ybs[8984]=['',6.2205755,-0.4559043,6.17];
ybs[8985]=['',6.2217105,0.9760707,6.51];
ybs[8986]=['',6.2273253,-0.6991344,6.31];
ybs[8987]=['107 Aqr',6.227239,-0.3238111,5.29];
ybs[8988]=['ψ And',6.2271486,0.8123692,4.95];
ybs[8989]=['19 Psc',6.2288318,0.0630365,5.04];
ybs[8990]=['',6.2295131,1.1677525,5.95];
ybs[8991]=['σ Phe',6.2327962,-0.8744375,5.18];
ybs[8992]=['',6.2334763,-1.1915201,6.89];
ybs[8993]=['τ Cas',6.231554,1.0258526,4.87];
ybs[8994]=['',6.2326744,-0.2057001,5.73];
ybs[8995]=['',6.2314462,1.004899,5.51];
ybs[8996]=['',6.2337872,0.8195646,6.07];
ybs[8997]=['20 Psc',6.2356078,-0.0460167,5.49];
ybs[8998]=['',6.2352084,1.1856379,5.04];
ybs[8999]=['',6.2382321,-0.1091779,6.07];
ybs[9000]=['',6.2394383,0.0408284,6.46];
ybs[9001]=['δ Scl',6.2399539,-0.488782,4.57];
ybs[9002]=['',6.2384731,1.1344904,6.41];
ybs[9003]=['6 Cas',6.2393153,1.0880308,5.43];
ybs[9004]=['',6.2396024,1.049013,6.34];
ybs[9005]=['',6.2409304,1.0312836,6.33];
ybs[9006]=['',6.2425452,-0.2746444,6.24];
ybs[9007]=['21 Psc',6.2422184,0.0209659,5.77];
ybs[9008]=['',6.2436502,-1.0945709,6.59];
ybs[9009]=['',6.2431334,0.6379253,5.9];
ybs[9010]=['79 Peg',6.2430343,0.5055809,5.97];
ybs[9011]=['',6.2438707,-0.4399318,6.42];
ybs[9012]=['',6.2456681,-0.1718975,5.94];
ybs[9013]=['',6.2460961,0.9031525,6.44];
ybs[9014]=['',6.247027,-0.2491767,5.72];
ybs[9015]=['80 Peg',6.2504761,0.1647332,5.79];
ybs[9016]=['108 Aqr',6.2505224,-0.3278375,5.18];
ybs[9017]=['γ1 Oct',6.2542659,-1.4293145,5.11];
ybs[9018]=['22 Psc',6.253154,0.0533281,5.55];
ybs[9019]=['',6.2528089,1.3565509,6.55];
ybs[9020]=['',6.2549854,0.3804126,6.11];
ybs[9021]=['φ Peg',6.2554182,0.3358971,5.08];
ybs[9022]=['',6.2555088,-0.2465435,5.87];
ybs[9023]=['',6.2548743,1.3206893,6.39];
ybs[9024]=['82 Peg',6.2559958,0.1932552,5.3];
ybs[9025]=['',6.2569928,-0.1548361,5.75];
ybs[9026]=['24 Psc',6.2573572,-0.0528894,5.93];
ybs[9027]=['25 Psc',6.2580206,0.0386725,6.28];
ybs[9028]=['',6.2592124,-0.4206932,6.24];
ybs[9029]=['',6.2636179,-0.4697901,6.35];
ybs[9030]=['ρ Cas',6.263634,1.0057403,4.54];
ybs[9031]=['',6.2648827,-0.7011819,6.03];
ybs[9032]=['',6.2654258,0.0040911,5.61];
ybs[9033]=['26 Psc',6.2669625,0.1256001,6.21];
ybs[9034]=['',6.2676322,-0.5549523,6.1];
ybs[9035]=['',6.2676322,-0.5542978,6.83];
ybs[9036]=['',6.2680542,0.4551861,6.54];
ybs[9037]=['',6.2687933,1.0042182,6];
ybs[9038]=['',6.2688038,0.8287011,6];
ybs[9039]=['',6.2729509,-0.4295599,6.31];
ybs[9040]=['',6.2737708,0.3974692,6.15];
ybs[9041]=['',6.2725371,1.4541449,6.59];
ybs[9042]=['',6.275368,0.7467145,5.97];
ybs[9043]=['',6.2757415,-0.4624836,6.26];
ybs[9044]=['',6.2757132,0.9744363,5.55];
ybs[9045]=['',6.2766131,-1.0966101,5.97];
ybs[9046]=['γ2 Oct',6.2776333,-1.4319508,5.73];
ybs[9047]=['η Tuc',6.2777232,-1.1200314,5];
ybs[9048]=['',6.2775322,1.0497958,6.47];
ybs[9049]=['ψ Peg',6.2784285,0.4409862,4.66];
ybs[9050]=['1 Cet',6.281035,-0.2744049,6.26];
ybs[9051]=['',6.2812818,0.8990867,4.8];
ybs[9052]=['27 Psc',6.2824288,-0.0598796,4.86];
ybs[9053]=['',6.2830647,0.5673529,6.52];
ybs[9054]=['π Phe',0.0003692,-0.9184022,5.13];
ybs[9055]=['',6.2828654,0.8122468,6.54];
ybs[9056]=['σ Cas',0.0006989,0.9752946,4.88];
ybs[9057]=['ω Psc',0.0020282,0.1219741,4.01];
ybs[9058]=['',0.0026971,-0.512424,5.62];
ybs[9059]=['',0.0027923,0.5907889,6.58];
ybs[9060]=['',0.0027923,0.5907889,6.58];
ybs[9061]=['ε Tuc',0.0046578,-1.1423521,4.5];
ybs[9062]=['',0.0064197,-0.7708298,6.29];
ybs[9063]=['',0.0067746,0.4719997,6.46];
ybs[9064]=['',0.0072966,1.0416995,6.19];
ybs[9065]=['',0.0082223,0.7920059,6.38];
ybs[9066]=['τ Phe',0.0097044,-0.849709,5.71];
ybs[9067]=['',0.0108276,-0.8763641,5.53];
ybs[9068]=['',0.0108203,0.8745308,6.22];
ybs[9069]=['θ Oct',0.0119019,-1.3428663,4.78];
ybs[9070]=['',0.0121168,1.0707301,5.55];
ybs[9071]=['',0.0126004,0.7416336,6.25];
ybs[9072]=['29 Psc',0.0129863,-0.0506537,5.1];
ybs[9073]=['85 Peg',0.0145136,0.4748552,5.75];
ybs[9074]=['30 Psc',0.0135813,-0.1027809,4.41];
ybs[9075]=['',0.0142825,-0.2539604,7.1];
ybs[9076]=['ζ Scl',0.0151897,-0.5165306,5.01];
ybs[9077]=['31 Psc',0.0155227,0.1585143,6.32];
ybs[9078]=['32 Psc',0.0159226,0.150287,5.63];
ybs[9079]=['',0.0164524,1.1558292,5.86];
ybs[9080]=['',0.0179349,-0.3476846,6.25];
ybs[9081]=['',0.0186658,-0.4192286,6.44];
ybs[9082]=['',0.0200681,1.1129474,6.24];
ybs[9083]=['2 Cet',0.0213377,-0.3003863,4.55];
ybs[9084]=['',0.0219949,1.1665337,6.29];
ybs[9085]=['9 Cas',0.0235617,1.0893126,5.88];
ybs[9086]=['',0.0239038,-0.2862977,5.78];
ybs[9087]=['',0.0239351,-0.5086478,6.4];
ybs[9088]=['3 Cet',0.024665,-0.1812387,4.94];
ybs[9089]=['',0.0256595,1.1744652,5.67];
ybs[9090]=['',0.0251988,0.7368336,6.01];
ybs[9091]=['',0.0245611,-1.2701205,7.31];
ybs[9092]=['',0.0264336,0.6071119,6.12];
ybs[9093]=['',0.0253399,-1.2446242,5.59];
ybs[9094]=['',0.0265839,0.4672965,6.25];
ybs[9095]=['',0.0273987,1.0723197,5.8];    
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