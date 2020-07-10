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

  var ceres = build_minor_planet('Ceres', 1, 3.34, 0.12, {
      equinox: when_j2000, 
      epoch: when("UT 2019-04-27"),
      a: 2.76916515450648,
      e: 0.07600902910070946,
      i: rads(10.59406704424526),
      Ω: rads(80.30553156826473),
      ω: rads(73.597694115971),
      M0: rads(77.37209588584763),
      n: rads(.213885225911375) 
    }
  );
  
  var vesta = build_minor_planet('Vesta', 4, 3.20, 0.32, {
      equinox: when_j2000, 
      epoch: when("UT 2019-04-27"),
      a: 2.361417895877277,
      e: 0.08872145950956178,
      i: rads(7.141770811873426),
      Ω: rads(103.810804427096),
      ω: rads(150.7285412870121),
      M0: rads(95.86193620017683),
      n: rads(0.2716094015313151) 
    }
  );
  
  var eunomia = build_minor_planet('Eunomia', 15, 5.28, 0.23, {
      equinox: when_j2000, 
      epoch: when("UT 2019-04-27"),
      a: 2.644100301686379,
      e: 0.1860843467561324,
      i: rads(11.75243014616266),
      Ω: rads(292.9343400475747),
      ω: rads(98.49867350140117),
      M0: rads(283.3877048843139),
      n: rads(0.2292383032694505) 
    }
  );
  
  var euterpe = build_minor_planet('Euterpe', 27, 7.0, 0.15, {
      equinox: when_j2000, 
      epoch: when("UT 2019-04-27"),
      a: 2.346664413776671,
      e: 0.173225849497636,
      i: rads(1.583713843590394),
      Ω: rads(94.78794181847543),
      ω: rads(356.4498803824397),
      M0: rads(335.316057551255),
      n: rads(0.2741748361944277) 
    }
  );
  
  var melpomene = build_minor_planet('Melpomene', 18, 6.51, 0.25, {
      equinox: when_j2000, 
      epoch: when("UT 2019-04-27"),
      a: 2.29665350984091,
      e: 0.2176743622521839,
      i: rads(10.12873129564968),
      Ω: rads(150.3838618286583),
      ω: rads(227.9508469405316),
      M0: rads(267.2543806758699),
      n: rads(0.2831788769147128) 
    }
  );
  
  var astraea = build_minor_planet('Astraea', 5, 6.85, 0.15, {
      equinox: when_j2000, 
      epoch: when("UT 2019-04-27"),
      a: 2.574248918456554,
      e: 0.191094508984634,
      i: rads(5.366988248754509),
      Ω: rads(141.5766048356492),
      ω: rads(358.6876073473166),
      M0: rads(282.3662885068957),
      n: rads(0.238631771610838) 
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
  var panstarrs_2017_T2 =  build_comet('PanSTARRS (2017 T2)', 'PanSTARRS (2017 T2)', 9.0, 'fade', 'all night', {
      equinox: when_j2000,
      epoch: when("TT 2019-03-03"),
      q: 1.614716047120907, 
      e: 0.999817074351556,
      i: rads(57.23616498855034),
      Ω: rads(64.37852784007528),
      ω: rads(93.00023485020597),
      T: when("TT 2020-05-04.91659834") 
    }
  );
  var neowise_2020_f3 =  build_comet('C/2020 F3 (NEOWISE)', 'C/2020 F3 (NEOWISE)', 1.0, 'fading', 'early morning', {
      equinox: when_j2000,
      epoch: when("TT 2020-04-19"),
      q: 0.2946766901808353, 
      e: 0.9991762343449366,
      i: rads(128.9375082086831),
      Ω: rads(61.00968684540688),
      ω: rads(37.27670558969492),
      T: when("TT 2020-07-03.68062954") 
    }
  );
  var lemmon_2019_u6 =  build_comet('C/2019 U6 (Lemmon)', 'C/2019 U6 (Lemmon)', 7.5, 'fade', 'evening', {
      equinox: when_j2000,
      epoch: when("TT 2019-12-23"),
      q: 0.9142791243582419, 
      e: 0.9979256826345774,
      i: rads(61.00264167203654),
      Ω: rads(235.7075243379916),
      ω: rads(329.6563185311068),
      T: when("TT 2020-06-18.81966081") 
    }
  );
  
  var comets = {
    /*testing only enke_test: enke_test,*/
    lemmon_2019_u6: lemmon_2019_u6,
    neowise_2020_f3: neowise_2020_f3, 
    panstarrs_2017_T2: panstarrs_2017_T2
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
     {when:'UT 2018-12-03 04',text:'Spica 7.70°S of Moon'},
     {when:'UT 2018-12-03 19',text:'Venus 3.63°S of Moon'},
     {when:'UT 2018-12-05 21',text:'Mercury 1.88°S of Moon'},
     {when:'UT 2018-12-06 13',text:'Jupiter 3.45°S of Moon'},
     {when:'UT 2018-12-06 19',text:'Antares 8.55°S of Moon'},
     {when:'UT 2018-12-07 07',text:'New Moon 394111.993 km'},
     {when:'UT 2018-12-07 15',text:'Neptune 0.04°S of Mars'},
     {when:'UT 2018-12-09 05',text:'Saturn 1.13°S of Moon'},
     {when:'UT 2018-12-12 12',text:'Moon at apogee 405176.775 km'},
     {when:'UT 2018-12-14 14',text:'Neptune 2.98°N of Moon'},
     {when:'UT 2018-12-14 23',text:'Mars 3.56°N of Moon'},
     {when:'UT 2018-12-15 12',text:'Mercury at greatest elongation 21.3° West'},
     {when:'UT 2018-12-18 04',text:'Uranus 5.01°N of Moon'},
     {when:'UT 2018-12-20 02',text:'Antares 5.27°S of Jupiter'},
     {when:'UT 2018-12-21 08',text:'Aldebaran 1.67°S of Moon'},
     {when:'UT 2018-12-21 15',text:'Jupiter 0.87°S of Mercury'},
     {when:'UT 2018-12-21 22',text:'Solstice'},
     {when:'UT 2018-12-22 18',text:'Full Moon 363368.422 km'},
     {when:'UT 2018-12-24 08',text:'Pollux 6.99°N of Moon'},
     {when:'UT 2018-12-24 10',text:'Moon at perigee 361061.498 km'},
     {when:'UT 2018-12-26 17',text:'Regulus 2.52°S of Moon'},
     {when:'UT 2018-12-29 10',text:'Last Quarter 375797.473 km'},
     {when:'UT 2018-12-30 10',text:'Spica 7.86°S of Moon'},
     {when:'UT 2019-01-01 22',text:'Venus 1.28°S of Moon'},
     {when:'UT 2019-01-02 06',text:'Saturn 0.48°N of Sun'},
     {when:'UT 2019-01-03 02',text:'Antares 8.60°S of Moon'},
     {when:'UT 2019-01-03 05',text:'Earth at perihelion 0.983301165 AU'},
     {when:'UT 2019-01-03 08',text:'Jupiter 3.13°S of Moon'},
     {when:'UT 2019-01-04 18',text:'Mercury 2.77°S of Moon'},
     {when:'UT 2019-01-05 19',text:'Saturn 0.87°S of Moon'},
     {when:'UT 2019-01-06 01',text:'PARTIAL ECLIPSE New Moon 402609.472 km'},
     {when:'UT 2019-01-06 05',text:'Venus at greatest elongation 47.0° West'},
     {when:'UT 2019-01-09 04',text:'Moon at apogee 406117.417 km'},
     {when:'UT 2019-01-10 22',text:'Neptune 3.12°N of Moon'},
     {when:'UT 2019-01-12 20',text:'Mars 5.34°N of Moon'},
     {when:'UT 2019-01-13 11',text:'Saturn 1.72°N of Mercury'},
     {when:'UT 2019-01-14 07',text:'First Quarter 392908.626 km'},
     {when:'UT 2019-01-14 12',text:'Uranus 5.11°N of Moon'},
     {when:'UT 2019-01-15 21',text:'Antares 7.93°S of Venus'},
     {when:'UT 2019-01-17 19',text:'Aldebaran 1.63°S of Moon'},
     {when:'UT 2019-01-20 20',text:'Pollux 6.99°N of Moon'},
     {when:'UT 2019-01-21 05',text:'TOTAL ECLIPSE Full Moon 357714.618 km'},
     {when:'UT 2019-01-21 20',text:'Moon at perigee 357342.248 km'},
     {when:'UT 2019-01-22 06',text:'Jupiter 2.43°S of Venus'},
     {when:'UT 2019-01-23 02',text:'Regulus 2.54°S of Moon'},
     {when:'UT 2019-01-26 16',text:'Spica 7.89°S of Moon'},
     {when:'UT 2019-01-27 21',text:'Last Quarter 382160.403 km'},
     {when:'UT 2019-01-30 03',text:'Mercury in superior conjunction 2.08° South'},
     {when:'UT 2019-01-30 07',text:'Antares 8.61°S of Moon'},
     {when:'UT 2019-01-31 00',text:'Jupiter 2.76°S of Moon'},
     {when:'UT 2019-01-31 18',text:'Venus 0.09°S of Moon'},
     {when:'UT 2019-02-02 07',text:'Saturn 0.62°S of Moon'},
     {when:'UT 2019-02-04 21',text:'New Moon 406469.645 km'},
     {when:'UT 2019-02-05 07',text:'Mercury 0.19°N of Moon'},
     {when:'UT 2019-02-05 09',text:'Moon at apogee 406555.016 km'},
     {when:'UT 2019-02-07 06',text:'Neptune 3.14°N of Moon'},
     {when:'UT 2019-02-10 16',text:'Mars 6.08°N of Moon'},
     {when:'UT 2019-02-10 20',text:'Uranus 5.06°N of Moon'},
     {when:'UT 2019-02-12 22',text:'First Quarter 384994.797 km'},
     {when:'UT 2019-02-13 20',text:'Uranus 1.05°S of Mars'},
     {when:'UT 2019-02-14 04',text:'Aldebaran 1.70°S of Moon'},
     {when:'UT 2019-02-17 07',text:'Pollux 6.96°N of Moon'},
     {when:'UT 2019-02-18 14',text:'Saturn 1.09°S of Venus'},
     {when:'UT 2019-02-19 09',text:'Moon at perigee 356760.688 km'},
     {when:'UT 2019-02-19 11',text:'Neptune 0.77°S of Mercury'},
     {when:'UT 2019-02-19 14',text:'Regulus 2.52°S of Moon'},
     {when:'UT 2019-02-19 16',text:'Full Moon 356842.872 km'},
     {when:'UT 2019-02-23 01',text:'Spica 7.78°S of Moon'},
     {when:'UT 2019-02-26 11',text:'Last Quarter 389808.992 km'},
     {when:'UT 2019-02-26 14',text:'Antares 8.48°S of Moon'},
     {when:'UT 2019-02-27 01',text:'Mercury at greatest elongation 18.1° East'},
     {when:'UT 2019-02-27 14',text:'Jupiter 2.32°S of Moon'},
     {when:'UT 2019-03-01 18',text:'Saturn 0.31°S of Moon'},
     {when:'UT 2019-03-02 21',text:'Venus 1.20°N of Moon'},
     {when:'UT 2019-03-04 11',text:'Moon at apogee 406390.596 km'},
     {when:'UT 2019-03-06 14',text:'Neptune 3.16°N of Moon'},
     {when:'UT 2019-03-06 16',text:'New Moon 404730.352 km'},
     {when:'UT 2019-03-07 01',text:'Neptune 0.96°S of Sun'},
     {when:'UT 2019-03-07 13',text:'Mercury 8.43°N of Moon'},
     {when:'UT 2019-03-10 04',text:'Uranus 4.90°N of Moon'},
     {when:'UT 2019-03-11 12',text:'Mars 5.78°N of Moon'},
     {when:'UT 2019-03-13 11',text:'Aldebaran 1.90°S of Moon'},
     {when:'UT 2019-03-14 10',text:'First Quarter 377824.517 km'},
     {when:'UT 2019-03-15 02',text:'Mercury in inferior conjunction 3.50° North'},
     {when:'UT 2019-03-16 16',text:'Pollux 6.81°N of Moon'},
     {when:'UT 2019-03-19 00',text:'Regulus 2.57°S of Moon'},
     {when:'UT 2019-03-19 20',text:'Moon at perigee 359377.000 km'},
     {when:'UT 2019-03-20 22',text:'Equinox'},
     {when:'UT 2019-03-21 02',text:'Full Moon 360767.844 km'},
     {when:'UT 2019-03-22 07',text:'Neptune 3.40°S of Mercury'},
     {when:'UT 2019-03-22 12',text:'Spica 7.64°S of Moon'},
     {when:'UT 2019-03-25 22',text:'Antares 8.25°S of Moon'},
     {when:'UT 2019-03-27 02',text:'Jupiter 1.90°S of Moon'},
     {when:'UT 2019-03-28 04',text:'Last Quarter 397206.134 km'},
     {when:'UT 2019-03-29 05',text:'Saturn 0.05°N of Moon'},
     {when:'UT 2019-04-01 00',text:'Moon at apogee 405577.361 km'},
     {when:'UT 2019-04-02 04',text:'Venus 2.68°N of Moon'},
     {when:'UT 2019-04-02 19',text:'Neptune 0.38°S of Mercury'},
     {when:'UT 2019-04-02 23',text:'Neptune 3.27°N of Moon'},
     {when:'UT 2019-04-05 09',text:'New Moon 398063.954 km'},
     {when:'UT 2019-04-06 13',text:'Uranus 4.77°N of Moon'},
     {when:'UT 2019-04-09 07',text:'Mars 4.72°N of Moon'},
     {when:'UT 2019-04-09 16',text:'Aldebaran 2.14°S of Moon'},
     {when:'UT 2019-04-10 04',text:'Neptune 0.31°N of Venus'},
     {when:'UT 2019-04-11 20',text:'Mercury at greatest elongation 27.7° West'},
     {when:'UT 2019-04-12 19',text:'First Quarter 372634.104 km'},
     {when:'UT 2019-04-12 23',text:'Pollux 6.56°N of Moon'},
     {when:'UT 2019-04-15 09',text:'Regulus 2.75°S of Moon'},
     {when:'UT 2019-04-16 22',text:'Aldebaran 6.54°S of Mars'},
     {when:'UT 2019-04-18 22',text:'Spica 7.61°S of Moon'},
     {when:'UT 2019-04-19 11',text:'Full Moon 368584.124 km'},
     {when:'UT 2019-04-22 08',text:'Antares 8.04°S of Moon'},
     {when:'UT 2019-04-22 23',text:'Uranus 0.49°S of Sun'},
     {when:'UT 2019-04-23 12',text:'Jupiter 1.64°S of Moon'},
     {when:'UT 2019-04-25 14',text:'Saturn 0.37°N of Moon'},
     {when:'UT 2019-04-26 22',text:'Last Quarter 402459.720 km'},
     {when:'UT 2019-04-28 18',text:'Moon at apogee 404581.976 km'},
     {when:'UT 2019-04-30 08',text:'Neptune 3.49°N of Moon'},
     {when:'UT 2019-05-02 12',text:'Venus 3.62°N of Moon'},
     {when:'UT 2019-05-03 06',text:'Mercury 2.91°N of Moon'},
     {when:'UT 2019-05-03 23',text:'Uranus 4.73°N of Moon'},
     {when:'UT 2019-05-04 23',text:'New Moon 388270.744 km'},
     {when:'UT 2019-05-06 22',text:'Aldebaran 2.29°S of Moon'},
     {when:'UT 2019-05-08 00',text:'Mars 3.23°N of Moon'},
     {when:'UT 2019-05-08 08',text:'Uranus 1.38°N of Mercury'},
     {when:'UT 2019-05-10 04',text:'Pollux 6.32°N of Moon'},
     {when:'UT 2019-05-12 01',text:'First Quarter 370050.384 km'},
     {when:'UT 2019-05-12 15',text:'Regulus 2.99°S of Moon'},
     {when:'UT 2019-05-13 22',text:'Moon at perigee 369008.823 km'},
     {when:'UT 2019-05-16 07',text:'Spica 7.69°S of Moon'},
     {when:'UT 2019-05-18 08',text:'Uranus 1.15°N of Venus'},
     {when:'UT 2019-05-18 21',text:'Full Moon 378792.915 km'},
     {when:'UT 2019-05-19 17',text:'Antares 7.94°S of Moon'},
     {when:'UT 2019-05-20 17',text:'Jupiter 1.69°S of Moon'},
     {when:'UT 2019-05-21 13',text:'Mercury in superior conjunction 0.33° North'},
     {when:'UT 2019-05-22 22',text:'Saturn 0.52°N of Moon'},
     {when:'UT 2019-05-26 12',text:'Aldebaran 6.67°S of Mercury'},
     {when:'UT 2019-05-26 13',text:'Moon at apogee 404137.606 km'},
     {when:'UT 2019-05-26 17',text:'Last Quarter 404126.779 km'},
     {when:'UT 2019-05-27 17',text:'Neptune 3.71°N of Moon'},
     {when:'UT 2019-05-31 10',text:'Uranus 4.79°N of Moon'},
     {when:'UT 2019-06-01 18',text:'Venus 3.24°N of Moon'},
     {when:'UT 2019-06-03 06',text:'Aldebaran 2.32°S of Moon'},
     {when:'UT 2019-06-03 10',text:'New Moon 377501.210 km'},
     {when:'UT 2019-06-04 16',text:'Mercury 3.66°N of Moon'},
     {when:'UT 2019-06-05 15',text:'Mars 1.58°N of Moon'},
     {when:'UT 2019-06-06 10',text:'Pollux 6.18°N of Moon'},
     {when:'UT 2019-06-07 23',text:'Moon at perigee 368503.889 km'},
     {when:'UT 2019-06-08 20',text:'Regulus 3.17°S of Moon'},
     {when:'UT 2019-06-10 06',text:'First Quarter 370204.768 km'},
     {when:'UT 2019-06-10 15',text:'Jupiter at opposition'},
     {when:'UT 2019-06-12 13',text:'Spica 7.83°S of Moon'},
     {when:'UT 2019-06-16 01',text:'Antares 7.95°S of Moon'},
     {when:'UT 2019-06-16 19',text:'Jupiter 1.99°S of Moon'},
     {when:'UT 2019-06-17 09',text:'Full Moon 389563.000 km'},
     {when:'UT 2019-06-17 21',text:'Aldebaran 4.80°S of Venus'},
     {when:'UT 2019-06-18 15',text:'Mars 0.24°S of Mercury'},
     {when:'UT 2019-06-19 04',text:'Saturn 0.44°N of Moon'},
     {when:'UT 2019-06-21 05',text:'Pollux 5.74°N of Mercury'},
     {when:'UT 2019-06-21 16',text:'Solstice'},
     {when:'UT 2019-06-23 07',text:'Pollux 5.60°N of Mars'},
     {when:'UT 2019-06-23 08',text:'Moon at apogee 404548.046 km'},
     {when:'UT 2019-06-23 23',text:'Mercury at greatest elongation 25.2° East'},
     {when:'UT 2019-06-24 01',text:'Neptune 3.84°N of Moon'},
     {when:'UT 2019-06-25 10',text:'Last Quarter 401871.683 km'},
     {when:'UT 2019-06-27 22',text:'Uranus 4.84°N of Moon'},
     {when:'UT 2019-06-30 16',text:'Aldebaran 2.29°S of Moon'},
     {when:'UT 2019-07-01 22',text:'Venus 1.64°N of Moon'},
     {when:'UT 2019-07-02 19',text:'TOTAL ECLIPSE New Moon 367744.913 km'},
     {when:'UT 2019-07-03 18',text:'Pollux 6.14°N of Moon'},
     {when:'UT 2019-07-04 06',text:'Mars 0.09°S of Moon'},
     {when:'UT 2019-07-04 09',text:'Mercury 3.25°S of Moon'},
     {when:'UT 2019-07-04 22',text:'Earth at aphelion 1.016754345 AU'},
     {when:'UT 2019-07-05 05',text:'Moon at perigee 363725.908 km'},
     {when:'UT 2019-07-06 03',text:'Regulus 3.24°S of Moon'},
     {when:'UT 2019-07-07 14',text:'Mars 3.84°N of Mercury'},
     {when:'UT 2019-07-09 11',text:'First Quarter 372965.473 km'},
     {when:'UT 2019-07-09 17',text:'Saturn at opposition'},
     {when:'UT 2019-07-09 18',text:'Spica 7.90°S of Moon'},
     {when:'UT 2019-07-13 07',text:'Antares 7.98°S of Moon'},
     {when:'UT 2019-07-13 20',text:'Jupiter 2.34°S of Moon'},
     {when:'UT 2019-07-16 07',text:'Saturn 0.22°N of Moon'},
     {when:'UT 2019-07-16 22',text:'PARTIAL ECLIPSE Full Moon 398907.789 km'},
     {when:'UT 2019-07-21 00',text:'Moon at apogee 405480.551 km'},
     {when:'UT 2019-07-21 08',text:'Neptune 3.81°N of Moon'},
     {when:'UT 2019-07-21 13',text:'Mercury in inferior conjunction 4.96° South'},
     {when:'UT 2019-07-23 16',text:'Pollux 6.09°N of Venus'},
     {when:'UT 2019-07-24 11',text:'Venus 5.72°N of Mercury'},
     {when:'UT 2019-07-25 01',text:'Last Quarter 396428.522 km'},
     {when:'UT 2019-07-25 07',text:'Uranus 4.80°N of Moon'},
     {when:'UT 2019-07-26 06',text:'Pollux 11.61°N of Mercury'},
     {when:'UT 2019-07-28 01',text:'Aldebaran 2.30°S of Moon'},
     {when:'UT 2019-07-31 02',text:'Mercury 4.52°S of Moon'},
     {when:'UT 2019-07-31 04',text:'Pollux 6.14°N of Moon'},
     {when:'UT 2019-07-31 21',text:'Venus 0.59°S of Moon'},
     {when:'UT 2019-08-01 03',text:'New Moon 360606.299 km'},
     {when:'UT 2019-08-01 20',text:'Mars 1.68°S of Moon'},
     {when:'UT 2019-08-02 07',text:'Moon at perigee 359397.901 km'},
     {when:'UT 2019-08-02 12',text:'Regulus 3.21°S of Moon'},
     {when:'UT 2019-08-05 22',text:'Pollux 9.33°N of Mercury'},
     {when:'UT 2019-08-06 01',text:'Spica 7.84°S of Moon'},
     {when:'UT 2019-08-07 18',text:'First Quarter 378084.205 km'},
     {when:'UT 2019-08-09 13',text:'Antares 7.93°S of Moon'},
     {when:'UT 2019-08-09 23',text:'Jupiter 2.47°S of Moon'},
     {when:'UT 2019-08-12 10',text:'Saturn 0.04°N of Moon'},
     {when:'UT 2019-08-14 06',text:'Venus in superior conjunction 1.27° North'},
     {when:'UT 2019-08-15 12',text:'Full Moon 404934.615 km'},
     {when:'UT 2019-08-17 11',text:'Moon at apogee 406244.477 km'},
     {when:'UT 2019-08-17 13',text:'Neptune 3.69°N of Moon'},
     {when:'UT 2019-08-17 23',text:'Regulus 0.70°S of Mars'},
     {when:'UT 2019-08-21 04',text:'Regulus 0.96°S of Venus'},
     {when:'UT 2019-08-21 15',text:'Uranus 4.65°N of Moon'},
     {when:'UT 2019-08-23 15',text:'Last Quarter 389163.356 km'},
     {when:'UT 2019-08-24 10',text:'Aldebaran 2.43°S of Moon'},
     {when:'UT 2019-08-24 13',text:'Mars 0.31°S of Venus'},
     {when:'UT 2019-08-27 15',text:'Pollux 6.08°N of Moon'},
     {when:'UT 2019-08-29 03',text:'Regulus 1.36°S of Mercury'},
     {when:'UT 2019-08-29 22',text:'Regulus 3.20°S of Moon'},
     {when:'UT 2019-08-30 01',text:'Mercury 1.94°S of Moon'},
     {when:'UT 2019-08-30 10',text:'Mars 3.06°S of Moon'},
     {when:'UT 2019-08-30 11',text:'New Moon 357223.789 km'},
     {when:'UT 2019-08-30 16',text:'Moon at perigee 357176.280 km'},
     {when:'UT 2019-09-02 09',text:'Spica 7.69°S of Moon'},
     {when:'UT 2019-09-02 11',text:'Mars 1.08°N of Sun'},
     {when:'UT 2019-09-03 11',text:'Mars 0.70°S of Mercury'},
     {when:'UT 2019-09-04 02',text:'Mercury in superior conjunction 1.71° North'},
     {when:'UT 2019-09-05 19',text:'Antares 7.75°S of Moon'},
     {when:'UT 2019-09-06 03',text:'First Quarter 385082.374 km'},
     {when:'UT 2019-09-06 07',text:'Jupiter 2.30°S of Moon'},
     {when:'UT 2019-09-08 14',text:'Saturn 0.04°N of Moon'},
     {when:'UT 2019-09-10 07',text:'Neptune at opposition'},
     {when:'UT 2019-09-13 14',text:'Moon at apogee 406377.344 km'},
     {when:'UT 2019-09-13 18',text:'Neptune 3.61°N of Moon'},
     {when:'UT 2019-09-13 22',text:'Venus 0.34°N of Mercury'},
     {when:'UT 2019-09-14 05',text:'Full Moon 406247.221 km'},
     {when:'UT 2019-09-17 20',text:'Uranus 4.47°N of Moon'},
     {when:'UT 2019-09-20 17',text:'Aldebaran 2.66°S of Moon'},
     {when:'UT 2019-09-22 03',text:'Last Quarter 381672.275 km'},
     {when:'UT 2019-09-23 08',text:'Equinox'},
     {when:'UT 2019-09-24 00',text:'Pollux 5.89°N of Moon'},
     {when:'UT 2019-09-26 09',text:'Regulus 3.28°S of Moon'},
     {when:'UT 2019-08-30 11',text:'New Moon 357223.789 km'},
     {when:'UT 2019-09-28 01',text:'Mars 4.06°S of Moon'},
     {when:'UT 2019-09-28 02',text:'Moon at perigee 357802.250 km'},
     {when:'UT 2019-09-28 23',text:'Spica 1.42°S of Mercury'},
     {when:'UT 2019-09-29 13',text:'Venus 4.37°S of Moon'},
     {when:'UT 2019-09-29 20',text:'Spica 7.56°S of Moon'},
     {when:'UT 2019-09-29 22',text:'Mercury 6.24°S of Moon'},
     {when:'UT 2019-10-03 01',text:'Spica 3.13°S of Venus'},
     {when:'UT 2019-10-03 03',text:'Antares 7.50°S of Moon'},
     {when:'UT 2019-10-03 20',text:'Jupiter 1.87°S of Moon'},
     {when:'UT 2019-10-05 17',text:'First Quarter 392930.717 km'},
     {when:'UT 2019-10-05 21',text:'Saturn 0.26°N of Moon'},
     {when:'UT 2019-10-10 18',text:'Moon at apogee 405898.952 km'},
     {when:'UT 2019-10-10 23',text:'Neptune 3.66°N of Moon'},
     {when:'UT 2019-10-13 21',text:'Full Moon 402366.196 km'},
     {when:'UT 2019-10-15 00',text:'Uranus 4.38°N of Moon'},
     {when:'UT 2019-10-17 22',text:'Aldebaran 2.90°S of Moon'},
     {when:'UT 2019-10-20 04',text:'Mercury at greatest elongation 24.6° East'},
     {when:'UT 2019-10-21 07',text:'Pollux 5.63°N of Moon'},
     {when:'UT 2019-10-21 13',text:'Last Quarter 375440.610 km'},
     {when:'UT 2019-10-23 18',text:'Regulus 3.49°S of Moon'},
     {when:'UT 2019-10-26 11',text:'Moon at perigee 361310.865 km'},
     {when:'UT 2019-10-26 17',text:'Mars 4.53°S of Moon'},
     {when:'UT 2019-10-27 06',text:'Spica 7.55°S of Moon'},
     {when:'UT 2019-10-28 04',text:'New Moon 363684.922 km'},
     {when:'UT 2019-10-28 08',text:'Uranus at opposition'},
     {when:'UT 2019-10-29 14',text:'Venus 3.91°S of Moon'},
     {when:'UT 2019-10-29 15',text:'Mercury 6.67°S of Moon'},
     {when:'UT 2019-10-30 08',text:'Venus 2.72°N of Mercury'},
     {when:'UT 2019-10-30 13',text:'Antares 7.31°S of Moon'},
     {when:'UT 2019-10-31 14',text:'Jupiter 1.31°S of Moon'},
     {when:'UT 2019-11-02 07',text:'Saturn 0.59°N of Moon'},
     {when:'UT 2019-11-04 10',text:'First Quarter 399886.039 km'},
     {when:'UT 2019-11-07 05',text:'Neptune 3.85°N of Moon'},
     {when:'UT 2019-11-07 09',text:'Moon at apogee 405057.885 km'},
     {when:'UT 2019-11-08 15',text:'Spica 3.06°S of Mars'},
     {when:'UT 2019-11-09 11',text:'Antares 3.95°S of Venus'},
     {when:'UT 2019-11-11 04',text:'Uranus 4.44°N of Moon'},
     {when:'UT 2019-11-11 15',text:'Mercury in inferior conjunction 0.02° North'},
     {when:'UT 2019-11-12 14',text:'Full Moon 393972.407 km'},
     {when:'UT 2019-11-14 04',text:'Aldebaran 3.03°S of Moon'},
     {when:'UT 2019-11-17 12',text:'Pollux 5.41°N of Moon'},
     {when:'UT 2019-11-19 21',text:'Last Quarter 371501.386 km'},
     {when:'UT 2019-11-20 00',text:'Regulus 3.71°S of Moon'},
     {when:'UT 2019-11-23 08',text:'Moon at perigee 366716.233 km'},
     {when:'UT 2019-11-23 16',text:'Spica 7.65°S of Moon'},
     {when:'UT 2019-11-24 09',text:'Mars 4.34°S of Moon'},
     {when:'UT 2019-11-24 14',text:'Jupiter 1.41°N of Venus'},
     {when:'UT 2019-11-25 03',text:'Mercury 1.91°S of Moon'},
     {when:'UT 2019-11-26 15',text:'New Moon 372887.395 km'},
     {when:'UT 2019-11-26 23',text:'Antares 7.24°S of Moon'},
     {when:'UT 2019-11-28 10',text:'Mercury at greatest elongation 20.1° West'},
     {when:'UT 2019-11-28 11',text:'Jupiter 0.73°S of Moon'},
     {when:'UT 2019-11-28 19',text:'Venus 1.87°S of Moon'},
     {when:'UT 2019-11-29 21',text:'Saturn 0.93°N of Moon'},
     {when:'UT 2019-12-04 07',text:'First Quarter 403933.279 km'},
     {when:'UT 2019-12-04 12',text:'Neptune 4.05°N of Moon'},
     {when:'UT 2019-12-05 04',text:'Moon at apogee 404445.826 km'},
     {when:'UT 2019-12-08 11',text:'Uranus 4.59°N of Moon'},
     {when:'UT 2019-12-11 05',text:'Saturn 1.81°N of Venus'},
     {when:'UT 2019-12-11 12',text:'Aldebaran 3.03°S of Moon'},
     {when:'UT 2019-12-12 05',text:'Full Moon 382863.920 km'},
     {when:'UT 2019-12-14 18',text:'Pollux 5.31°N of Moon'},
     {when:'UT 2019-12-15 16',text:'Antares 5.14°S of Mercury'},
     {when:'UT 2019-12-17 05',text:'Regulus 3.84°S of Moon'},
     {when:'UT 2019-12-18 20',text:'Moon at perigee 370264.727 km'},
     {when:'UT 2019-12-19 05',text:'Last Quarter 370293.435 km'},
     {when:'UT 2019-12-20 22',text:'Spica 7.75°S of Moon'},
     {when:'UT 2019-12-22 04',text:'Solstice'},
     {when:'UT 2019-12-23 02',text:'Mars 3.53°S of Moon'},
     {when:'UT 2019-12-24 08',text:'Antares 7.26°S of Moon'},
     {when:'UT 2019-12-25 11',text:'Mercury 1.94°S of Moon'},
     {when:'UT 2019-12-26 05',text:'ANNULAR ECLIPSE New Moon 384225.938 km'},
     {when:'UT 2019-12-26 07',text:'Jupiter 0.18°S of Moon'},
     {when:'UT 2019-12-27 12',text:'Saturn 1.20°N of Moon'},
     {when:'UT 2019-12-27 18',text:'Jupiter 0.10°N of Sun'},
     {when:'UT 2019-12-29 02',text:'Venus 1.00°N of Moon'},
     {when:'UT 2019-12-31 21',text:'Neptune 4.14°N of Moon'},
     {when:'UT 2020-01-02 02',text:'Moon at apogee 404580.190 km'},
     {when:'UT 2020-01-02 15',text:'Jupiter 1.50°N of Mercury'},
     {when:'UT 2020-01-03 05',text:'First Quarter 403730.503 km'},
     {when:'UT 2020-01-04 18',text:'Uranus 4.67°N of Moon'},
     {when:'UT 2020-01-05 08',text:'Earth at perihelion 0.983243564 AU'},
     {when:'UT 2020-01-07 22',text:'Aldebaran 3.02°S of Moon'},
     {when:'UT 2020-01-10 15',text:'Mercury in superior conjunction 1.93° South'},
     {when:'UT 2020-01-10 19',text:'Full Moon 371542.737 km'},
     {when:'UT 2020-01-11 03',text:'Pollux 5.31°N of Moon'},
     {when:'UT 2020-01-12 05',text:'Saturn 2.06°N of Mercury'},
     {when:'UT 2020-01-13 12',text:'Regulus 3.84°S of Moon'},
     {when:'UT 2020-01-13 15',text:'Saturn 0.04°N of Sun'},
     {when:'UT 2020-01-13 20',text:'Moon at perigee 365958.483 km'},
     {when:'UT 2020-01-17 03',text:'Spica 7.73°S of Moon'},
     {when:'UT 2020-01-17 04',text:'Antares 4.82°S of Mars'},
     {when:'UT 2020-01-17 13',text:'Last Quarter 371872.008 km'},
     {when:'UT 2020-01-20 15',text:'Antares 7.24°S of Moon'},
     {when:'UT 2020-01-20 19',text:'Mars 2.25°S of Moon'},
     {when:'UT 2020-01-23 03',text:'Jupiter 0.36°N of Moon'},
     {when:'UT 2020-01-24 02',text:'Saturn 1.45°N of Moon'},
     {when:'UT 2020-01-24 22',text:'New Moon 395268.915 km'},
     {when:'UT 2020-01-25 18',text:'Mercury 1.33°N of Moon'},
     {when:'UT 2020-01-27 19',text:'Neptune 0.08°N of Venus'},
     {when:'UT 2020-01-28 06',text:'Neptune 4.11°N of Moon'},
     {when:'UT 2020-01-28 07',text:'Venus 4.08°N of Moon'},
     {when:'UT 2020-01-29 21',text:'Moon at apogee 405392.731 km'}
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

  /* Coords have equinox 2020.5. */
  var build_messiers = function(){
    var messier = [110];
    messier[0]=["M1",1.464917,0.3844791,8.4,"Tau","NB","!! famous Crab Neb. supernova remnant","Crab Nebula"];
    messier[1]=["M2",5.6485251,0.0158548,6.5,"Aqr","GC","200-mm telescope needed to resolve",""];
    messier[2]=["M3",3.5916435,0.4935872,6.2,"CVn","GC","!! contains many variable stars",""];
    messier[3]=["M4",4.297259,-0.4639024,5.6,"Sco","GC","bright globular near Antares",""];
    messier[4]=["M5",4.0126787,0.0350748,5.6,"Ser","GC","!! one of the sky's finest globulars",""];
    messier[5]=["M6",4.6313944,-0.5624539,4.2,"Sco","OC","!! Butterfly Cluster; best at low power","Butterfly Cluster"];
    messier[6]=["M7",4.6917421,-0.6077126,3.3,"Sco","OC","!! excellent in binocs or rich-field scope","Ptolemy's Cluster"];
    messier[7]=["M8",4.7344567,-0.4255309,6,"Sgr","NB","!! Lagoon Nebula w/open cl. NGC 6530","Lagoon Nebula"];
    messier[8]=["M9",4.5396071,-0.3235244,7.7,"Oph","GC","smallest of Ophiuchus globulars",""];
    messier[9]=["M10",4.4426585,-0.0720938,6.6,"Oph","GC","rich globular cluster; M12 is 3°NW",""];
    messier[10]=["M11",4.940152,-0.1089289,6.3,"Sct","OC","!! Wild Duck Cl.; the best open cluster?","Wild Duck Cluster"];
    messier[11]=["M12",4.3993886,-0.0346516,6.7,"Oph","GC","loose globular cluster near M10",""];
    messier[12]=["M13",4.3739383,0.635799,5.8,"Her","GC","!! Hercules Cluster; NGC 6207 0.5°NE","Great Hercules Globular"];
    messier[13]=["M14",4.6193477,-0.0569129,7.6,"Oph","GC","200-mm telescope needed to resolve",""];
    messier[14]=["M15",5.6330095,0.2139313,6.2,"Peg","GC","rich, compact globular","Great Pegasus Globular"];
    messier[15]=["M16",4.7994907,-0.2403963,6.4,"Ser","NB","Eagle Neb. w/open cl.; use neb. filter","Eagle Nebula"];
    messier[16]=["M17",4.8083059,-0.2822668,7,"Sgr","NB","!! Swan or Omega Nebula; use neb. filter","Omega Nebula"];
    messier[17]=["M18",4.804415,-0.2988552,7.5,"Sgr","OC","sparse cluster; 1°S of M17",""];
    messier[18]=["M19",4.4674723,-0.4589281,6.8,"Oph","GC","oblate globular; M62 4°S",""];
    messier[19]=["M20",4.7291648,-0.4019795,9,"Sgr","NB","!! Trifid Nebula; look for dark lanes","Trifid Nebula"];
    messier[20]=["M21",4.7378696,-0.3926537,6.5,"Sgr","OC","0.7°NE of M20; sparse cluster",""];
    messier[21]=["M22",4.8766692,-0.4168133,5.1,"Sgr","GC","spectacular from southern latitude","Sagittarius Cluster"];
    messier[22]=["M23",4.7036972,-0.331926,6.9,"Sgr","OC","bright, loose open cluster",""];
    messier[23]=["M24",4.7913779,-0.3227339,4.6,"Sgr","OC","rich star cloud; best in big binoculars","Sagittarius Star Cloud"];
    messier[24]=["M25",4.8555429,-0.3356969,6.5,"Sgr","OC","bright but sparse open cluster",""];
    messier[25]=["M26",4.9145185,-0.1636658,8,"Sct","OC","bright, coarse cluster",""];
    messier[26]=["M27",5.2381037,0.3974769,7.4,"Vul","NB","!! Dumbbell Nebula; a superb object","Dumbbell Nebula"];
    messier[27]=["M28",4.8247922,-0.4337872,6.8,"Sgr","GC","compact globular near M22",""];
    messier[28]=["M29",5.3435718,0.6737063,7.1,"Cyg","OC","small, poor open cluster 2°S of γ Cygni",""];
    messier[29]=["M30",5.679135,-0.4029889,7.2,"Cap","GC","toughest in one-night Messier marathon",""];
    messier[30]=["M31",0.1872932,0.7221972,3.4,"And","GY","!! Andromeda Gal.; look for dust lanes","Andromeda Galaxy"];
    messier[31]=["M32",0.1916594,0.7152143,8.1,"And","GY","closest companion to M31",""];
    messier[32]=["M33",0.4147744,0.5367685,5.7,"Tri","GY","large diffuse spiral; requires dark sky","Triangulum Galaxy"];
    messier[33]=["M34",0.7126459,0.7482209,5.5,"Per","OC","best at low power",""];
    messier[34]=["M35",1.6151142,0.424614,5.3,"Gem","OC","!! look for sm. cluster NGC 2158 0.25°S",""];
    messier[35]=["M36",1.472441,0.5959406,6.3,"Aur","OC","bright but scattered group; use low pow.",""];
    messier[36]=["M37",1.5434905,0.5681649,6.2,"Aur","OC","!! finest of three Auriga clusters; very rich",""];
    messier[37]=["M38",1.4402347,0.6256749,7.4,"Aur","OC","look for small cluster NGC 1907 0.5°S",""];
    messier[38]=["M39",5.641521,0.8469149,4.6,"Cyg","OC","very sparse cluster; use low power",""];
    messier[39]=["M40",3.2435974,1.0117634,8.4,"UMa","OC","double star Winneke 4; separation 50arcsec","Winnecke 4"];
    messier[40]=["M41",1.7797186,-0.3622743,4.6,"CMa","OC","4°S of Sirius; bright but coarse",""];
    messier[41]=["M42",1.4678543,-0.0949114,4,"Ori","NB","!! Orion Nebula; finest in northern sky","Great Nebula in Orion"];
    messier[42]=["M43",1.4687332,-0.0917134,9,"Ori","NB","detached part of Orion Nebula","De Mairan's Nebula"];
    messier[43]=["M44",2.2745012,0.34749,3.7,"Cnc","OC","!! Beehive or Praesepe; use low power","Beehive Cluster"];
    messier[44]=["M45",0.995807,0.4220029,1.6,"Tau","OC","!! Pleiades; look for subtle nebulosity","Pleiades"];
    messier[45]=["M46",2.0190909,-0.2594593,6,"Pup","OC","!! contains planetary nebula NGC 2438",""];
    messier[46]=["M47",1.9964073,-0.2538914,5.2,"Pup","OC","coarse cluster 1.5°W of M46",""];
    messier[47]=["M48",2.1590239,-0.1023307,5.5,"Hya","OC","former lost Messier; large sparse cl.",""];
    messier[48]=["M49",3.2761675,0.1376518,8.4,"Vir","GY","very bright elliptical",""];
    messier[49]=["M50",1.8508616,-0.1459906,6.3,"Mon","OC","between Sirius & Procyon; use low mag",""];
    messier[50]=["M51",3.5380512,0.8216656,8.4,"CVn","GY","!! Whirlpool Galaxy; superb in big scope","Whirlpool Galaxy"];
    messier[51]=["M52",6.130996,1.0768002,7.3,"Cas","OC","young, rich cl.; faint Bubble Neb. nearby",""];
    messier[52]=["M53",3.4640582,0.3151775,7.6,"Com","GC","150-mm telescope needed to resolve",""];
    messier[53]=["M54",4.9585299,-0.5315547,7.6,"Sgr","GC","not easily resolved",""];
    messier[54]=["M55",5.1543866,-0.5396234,6.3,"Sgr","GC","bright, loose globular cluster",""];
    messier[55]=["M56",5.0501093,0.5274552,8.3,"Lyr","GC","within a rich starfield",""];
    messier[56]=["M57",4.9495874,0.5770052,8.8,"Lyr","NB","!! Ring Nebula; an amazing smoke ring","Ring Nebula"];
    messier[57]=["M58",3.3106055,0.2042754,9.7,"Vir","GY","bright barred spiral; M59 and M60 1°E",""];
    messier[58]=["M59",3.3293613,0.2013731,9.6,"Vir","GY","bright elliptical paired with M60",""];
    messier[59]=["M60",3.3367767,0.1996305,8.8,"Vir","GY","bright elliptical with M59 and NGC 4647",""];
    messier[60]=["M61",3.2417188,0.0759756,9.7,"Vir","GY","face-on two-armed spiral",""];
    messier[61]=["M62",4.4615292,-0.526135,6.5,"Oph","GC","asymmetrical; in rich field",""];
    messier[62]=["M63",3.4763315,0.7317374,8.6,"CVn","GY","!! Sunflower Galaxy; bright, elongated","Sunflower Galaxy"];
    messier[63]=["M64",3.3933825,0.3765154,8.5,"Com","GY","!! Black Eye Gal; eye needs big scope","Black Eye Galaxy"];
    messier[64]=["M65",2.9669257,0.2263864,9.3,"Leo","GY","!! bright elongated spiral",""];
    messier[65]=["M66",2.9725948,0.2246391,8.9,"Leo","GY","!! M65 and NGC 3628 in same field",""];
    messier[66]=["M67",2.3191961,0.2048878,6.1,"Cnc","OC","one of the oldest star clusters known",""];
    messier[67]=["M68",3.3187033,-0.4688372,7.8,"Hya","GC","150-mm telescope needed to resolve",""];
    messier[68]=["M69",4.8552307,-0.5643362,7.6,"Sgr","GC","small, poor globular cluster",""];
    messier[69]=["M70",4.9067047,-0.5633624,7.9,"Sgr","GC","small globular 2°E of M69",""];
    messier[70]=["M71",5.2129237,0.3287834,8.2,"Sge","GC","loose globular; looks like an open cluster",""];
    messier[71]=["M72",5.47433,-0.2173764,9.3,"Aqr","GC","near the Saturn Nebula, NGC 7009",""];
    messier[72]=["M73",5.4983235,-0.2190875,9,"Aqr","OC","group of four stars only; an asterism",""];
    messier[73]=["M74",0.4267504,0.2772864,9.4,"Psc","GY","faint, elusive spiral; tough in small scope",""];
    messier[74]=["M75",5.2678694,-0.381472,8.5,"Sgr","GC","small and distant; 59,000 ly away",""];
    messier[75]=["M76",0.4524818,0.901802,10.1,"Per","NB","Little Dumbell; faint but distinct","Little Dumbbell Nebula"];
    messier[76]=["M77",0.7144989,0.0020895,8.9,"Cet","GY","a Seyfert galaxy; with starlike nucleus",""];
    messier[77]=["M78",1.5173504,0.0009837,8.3,"Ori","NB","bright featureless reflection nebula",""];
    messier[78]=["M79",1.419584,-0.4281746,7.7,"Lep","GC","200-mm telescope needed to resolve",""];
    messier[79]=["M80",4.268314,-0.4019954,7.3,"Sco","GC","very compressed globular",""];
    messier[80]=["M81",2.6060464,1.2037314,6.9,"UMa","GY","!! bright spiral visible in binoculars","Bode's Galaxy"];
    messier[81]=["M82",2.6070032,1.2144933,8.4,"UMa","GY","!! the exploding galaxy; M81 0.5°S","Cigar Galaxy"];
    messier[82]=["M83",3.5698929,-0.5230858,7.6,"Hya","GY","large and diffuse; superb from far south","Southern Pinwheel"];
    messier[83]=["M84",3.2556459,0.2228771,9.1,"Vir","GY","!! w/M86 in Markarian's Chain",""];
    messier[84]=["M85",3.2573679,0.3156708,9.1,"Com","GY","bright elliptical shape",""];
    messier[85]=["M86",3.2604431,0.2240417,8.9,"Vir","GY","!! w/many NGC galaxies in Chain",""];
    messier[86]=["M87",3.2805081,0.2144474,8.6,"Vir","GY","famous jet and black hole",""];
    messier[87]=["M88",3.2861673,0.2499374,9.6,"Com","GY","bright multiple-arm spiral",""];
    messier[88]=["M89",3.3018783,0.2170717,9.8,"Vir","GY","elliptical; resembles M87 but smaller",""];
    messier[89]=["M90",3.3066721,0.2278361,9.5,"Vir","GY","bright barred spiral near M89",""];
    messier[90]=["M91",3.3009947,0.2511053,10.2,"Com","GY","some lists say M91=M58, not NGC 4548",""];
    messier[91]=["M92",4.5279534,0.7524507,6.4,"Her","GC","9°NE of M13; fine but often overlooked",""];
    messier[92]=["M93",2.030993,-0.4174332,6,"Pup","OC","compact, bright cluster; fairly rich",""];
    messier[93]=["M94",3.3678842,0.71597,8.2,"CVn","GY","very bright and comet-like",""];
    messier[94]=["M95",2.8146973,0.2023186,9.7,"Leo","GY","bright barred spiral",""];
    messier[95]=["M96",2.8269111,0.204347,9.2,"Leo","GY","M95 in same field",""];
    messier[96]=["M97",2.9495043,0.9582677,9.9,"UMa","NB","Owl Nebula; distinct grey oval","Owl Nebula"];
    messier[97]=["M98",3.206794,0.2583569,10.1,"Com","GY","nearly edge-on spiral near star 6 Com. B.",""];
    messier[98]=["M99",3.2286007,0.2499244,9.9,"Com","GY","nearly face-on spiral near M98",""];
    messier[99]=["M100",3.2464759,0.2743623,9.3,"Com","GY","face-on spiral with starlike nucleus",""];
    messier[100]=["M101",3.6823153,0.946877,7.9,"UMa","GY","!! Pinwheel Gal; diffuse face-on spiral","Pinwheel Galaxy"];
    messier[101]=["M102",3.9578097,0.9719457,9.9,"Dra","GY","or is M102=M101? (look for NGC 5907)",""];
    messier[102]=["M103",0.412663,1.061242,7.4,"Cas","OC","three NGC open clusters nearby",""];
    messier[103]=["M104",3.3207824,-0.20471,8,"Vir","GY","!! Sombrero Galaxy; look for dust lane","Sombrero Galaxy"];
    messier[104]=["M105",2.8312813,0.2177252,9.3,"Leo","GY","bright elliptical near M95 and M96",""];
    messier[105]=["M106",3.2284616,0.8238468,8.4,"CVn","GY","!! superb large, bright spiral",""];
    messier[106]=["M107",4.3356122,-0.228503,7.9,"Oph","GC","small, faint globular",""];
    messier[107]=["M108",2.9351599,0.969618,10,"UMa","GY","nearly edge-on; paired with M97 0.75°SE",""];
    messier[108]=["M109",3.135727,0.929723,9.8,"UMa","GY","barred spiral near γ UMA",""];
    messier[109]=["M110",0.1811787,0.7294716,8.5,"And","GY","more distant companion to M31",""];
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
  
  /* Coords have equinox 2020.5. */
  var build_caldwells = function(){
    var caldwell = [42];
    caldwell[0]=["C68",4.9885051,-0.6443619,9.7,"CrA","Bn","","","NGC_6729"];
    caldwell[1]=["C69",4.5164289,-0.647911,12.8,"Sco","Pl","","Bug Nebula","NGC_6302"];
    caldwell[2]=["C70",0.2437634,-0.6557641,8.1,"Scl","Sp","","","NGC_300"];
    caldwell[3]=["C71",2.0639813,-0.6737647,5.8,"Pup","Oc","","","NGC_2477"];
    caldwell[4]=["C72",0.0694891,-0.6818908,8.2,"Scl","Sb","","","NGC_55"];
    caldwell[5]=["C73",1.3734635,-0.698611,7.3,"Col","Gc","","","NGC_1851"];
    caldwell[6]=["C74",2.6553786,-0.7074541,8.2,"Vel","Pl","","Eight Burst Nebula","NGC_3132"];
    caldwell[7]=["C75",4.3066473,-0.710559,5.8,"Sco","Oc","","","NGC_6124"];
    caldwell[8]=["C76",4.4307043,-0.7301073,2.6,"Sco","Oc","","","NGC_6231"];
    caldwell[9]=["C77",3.5199245,-0.7526355,7,"Cen","Px","","Centaurus A","Centaurus_A"];
    caldwell[10]=["C78",4.753782,-0.7626329,6.6,"CrA","Gc","","","NGC_6541"];
    caldwell[11]=["C79",2.6984704,-0.8119217,6.7,"Vel","Gc","","","NGC_3201"];
    caldwell[12]=["C80",3.5257238,-0.8305893,3.6,"Cen","Gc","","Omega Centauri","Omega_Centauri"];
    caldwell[13]=["C81",4.56866,-0.8453222,8.1,"Ara","Gc","","","NGC_6352"];
    caldwell[14]=["C82",4.3757239,-0.8518032,5.2,"Ara","Oc","","","NGC_6193"];
    caldwell[15]=["C83",3.4322013,-0.8652661,9.5,"Cen","Sp","","","NGC_4945"];
    caldwell[16]=["C84",3.611559,-0.898296,7.6,"Cen","Gc","","","NGC_5286"];
    caldwell[17]=["C85",2.2723562,-0.9274718,2.5,"Vel","Oc","","Omicron Vel Cluster","IC_2391"];
    caldwell[18]=["C86",4.6354613,-0.9368203,5.6,"Ara","Gc","","","NGC_6397"];
    caldwell[19]=["C87",0.8415184,-0.9623835,8.4,"Hor","Gc","","","NGC_1261"];
    caldwell[20]=["C88",3.9585636,-0.9717713,7.9,"Cir","Oc","","","NGC_5823"];
    caldwell[21]=["C89",4.2787206,-1.0113894,5.4,"Nor","Oc","","S Norma Cluster","NGC_6087"];
    caldwell[22]=["C90",2.4520943,-1.0193531,9.7,"Car","Pl","","","NGC_2867"];
    caldwell[23]=["C91",2.9115491,-1.0258651,3,"Car","Oc","","","NGC_3532"];
    caldwell[24]=["C92",2.812575,-1.0467544,6.2,"Car","Bn","","Eta Carinae Nebula","Carina_Nebula"];
    caldwell[25]=["C93",5.0296107,-1.0462928,5.4,"Pav","Gc","","","NGC_6752"];
    caldwell[26]=["C94",3.3808727,-1.0549517,4.2,"Cru","Oc","","","Jewel Box"];
    caldwell[27]=["C95",4.212606,-1.0568855,5.1,"TrA","Oc","","","NGC_6025"];
    caldwell[28]=["C96",2.0884514,-1.0633081,3.8,"Car","Oc","","","NGC_2516"];
    caldwell[29]=["C97",3.0415167,-1.0773952,5.3,"Cen","Oc","","","NGC_3766"];
    caldwell[30]=["C98",3.331474,-1.1009328,6.9,"Cru","Oc","","","NGC_4609"];
    caldwell[31]=["C99",3.3783419,-1.101495,undefined,"Cru","Dn","","Coalsack Nebula","Coalsack_Nebula"];
    caldwell[32]=["C100",3.0436835,-1.1021212,4.5,"Cen","Oc","","Lambda Centauri Nebula","IC_2944"];
    caldwell[33]=["C101",5.0253956,-1.1137874,9,"Pav","Sb","","","NGC_6744"];
    caldwell[34]=["C102",2.8097096,-1.1258742,1.9,"Car","Oc","","Theta Car Cluster","IC_2602"];
    caldwell[35]=["C103",1.4772497,-1.205837,1,"Dor","Bn","","Tarantula Nebula","Tarantula_Nebula"];
    caldwell[36]=["C104",0.2787811,-1.2346499,6.6,"Tuc","Gc","","","NGC_362"];
    caldwell[37]=["C105",3.4077304,-1.239071,7.3,"Mus","Gc","","","NGC_4833"];
    caldwell[38]=["C106",0.1090841,-1.256111,4,"Tuc","Gc","","47 Tucanae","47_Tucanae"];
    caldwell[39]=["C107",4.3116561,-1.2609142,9.3,"Aps","Gc","","","NGC_6101"];
    caldwell[40]=["C108",3.2594872,-1.2702513,7.8,"Mus","Gc","","","NGC_4372"];
    caldwell[41]=["C109",2.6582457,-1.4131539,11.6,"Cha","Pl","","","NGC_3195"];
    return caldwell;
  };
  var caldwells = build_caldwells();

  /* Find a Caldwell object either by 'C103' or 'Tarantula Nebula' (case-insensitive). */  
  var find_caldwell = function(name_raw){
    return find_deep_sky_object(name_raw, caldwells);
  };
  
  /* This exist in order to avoid creating a large number of identical objects. */
  var when_fixed_equinox = when("J2020.5");
  
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

  /* Coords have equinox 2020.5. */
  var build_stars = function(){
    /* Source: Yale Bright Star Catalog r5.  Name, Right Ascension, Declination (J2019.5), and Magnitude.*/
    var ybs = [9096];
    //removal of leading whitespace cuts down on the overall size of this js file
ybs[0]=['',0.027171,0.7913892,6.7];
ybs[1]=['',0.0266771,-0.0067887,6.29];
ybs[2]=['33 Psc',0.0278577,-0.0976234,4.61];
ybs[3]=['86 Peg',0.0294683,0.2357974,5.51];
ybs[4]=['',0.0320244,1.0219033,5.96];
ybs[5]=['',0.0320778,-0.8545293,5.7];
ybs[6]=['10 Cas',0.0328173,1.1224245,5.59];
ybs[7]=['',0.0334751,0.5085098,6.13];
ybs[8]=['',0.0343805,-0.4013111,6.18];
ybs[9]=['',0.03643,-0.3014589,6.19];
ybs[10]=['',0.0383316,-0.0424959,6.43];
ybs[11]=['',0.0385013,-0.3908636,5.94];
ybs[12]=['',0.0396963,-0.5832087,5.68];
ybs[13]=['',0.0403678,-0.0407313,6.07];
ybs[14]=['α And',0.0412287,0.5097164,2.06];
ybs[15]=['',0.0407445,-0.1520154,5.99];
ybs[16]=['',0.0425323,0.6412462,6.19];
ybs[17]=['',0.041895,-0.3047948,6.06];
ybs[18]=['',0.0433261,0.4463996,6.23];
ybs[19]=['',0.0458009,1.3932745,6.01];
ybs[20]=['β Cas',0.0447742,1.0343475,2.27];
ybs[21]=['87 Peg',0.0440563,0.3198486,5.53];
ybs[22]=['',0.0439148,-0.9405215,6.33];
ybs[23]=['κ1 Scl',0.0453361,-0.4864888,5.42];
ybs[24]=['ε Phe',0.0455622,-0.7964544,3.88];
ybs[25]=['34 Psc',0.0484031,0.1965164,5.51];
ybs[26]=['22 And',0.0497194,0.8061017,5.03];
ybs[27]=['',0.0505265,0.9997168,6.74];
ybs[28]=['',0.0495763,-0.0896158,5.84];
ybs[29]=['γ3 Oct',0.0477082,-1.4330878,5.28];
ybs[30]=['',0.0513085,-0.2175729,5.85];
ybs[31]=['',0.0507051,-1.276018,6.64];
ybs[32]=['6 Cet',0.053709,-0.2679792,4.89];
ybs[33]=['κ2 Scl',0.0550274,-0.4832075,5.41];
ybs[34]=['θ Scl',0.055706,-0.6111984,5.25];
ybs[35]=['',0.0570006,0.8424087,6.16];
ybs[36]=['',0.057636,-0.3110941,5.25];
ybs[37]=['',0.0606995,0.6598613,6.73];
ybs[38]=['γ Peg',0.0623729,0.2669923,2.83];
ybs[39]=['',0.0631147,0.4730042,6.3];
ybs[40]=['23 And',0.0636539,0.7181889,5.72];
ybs[41]=['',0.0643091,-0.4521805,5.94];
ybs[42]=['',0.0644609,-0.4567668,6.31];
ybs[43]=['',0.0659213,0.5815439,6.25];
ybs[44]=['χ Peg',0.0683521,0.3546605,4.8];
ybs[45]=['',0.0676605,-0.1338086,5.12];
ybs[46]=['',0.0613757,-1.4814397,5.77];
ybs[47]=['7 Cet',0.0684185,-0.3284517,4.44];
ybs[48]=['',0.0698058,0.3909195,6.24];
ybs[49]=['35 Psc',0.0699681,0.15594,5.79];
ybs[50]=['',0.0696119,-0.1650358,5.75];
ybs[51]=['',0.0706269,0.5523914,6.45];
ybs[52]=['',0.0708757,0.4781664,6.35];
ybs[53]=['',0.0698097,-0.6072101,6.17];
ybs[54]=['',0.0760492,1.345032,6.35];
ybs[55]=['',0.0761085,0.7628579,6.15];
ybs[56]=['',0.0749565,-0.5468563,5.67];
ybs[57]=['',0.0734841,-1.3229169,6.49];
ybs[58]=['36 Psc',0.0768991,0.1458015,6.11];
ybs[59]=['',0.0788311,1.0759454,5.74];
ybs[60]=['',0.0774334,-0.3507544,6.47];
ybs[61]=['',0.0795934,0.8388277,5.89];
ybs[62]=['θ And',0.0792838,0.6771084,4.61];
ybs[63]=['',0.0772096,-1.3729938,6.77];
ybs[64]=['',0.0820872,0.8996618,6.14];
ybs[65]=['',0.0810777,-0.3305189,6.45];
ybs[66]=['',0.0822346,0.0314623,6.17];
ybs[67]=['σ And',0.0846796,0.6440094,4.52];
ybs[68]=['',0.0844076,0.1975639,6.05];
ybs[69]=['26 And',0.0863463,0.766284,6.11];
ybs[70]=['',0.0860118,0.5520643,5.87];
ybs[71]=['',0.0861406,-0.1385626,6.46];
ybs[72]=['',0.0860657,-0.7526131,6.33];
ybs[73]=['ι Cet',0.0893298,-0.1520215,3.56];
ybs[74]=['',0.0906644,0.712852,6.33];
ybs[75]=['',0.0924345,0.8548439,6.52];
ybs[76]=['ζ Tuc',0.0917842,-1.1302936,4.23];
ybs[77]=['',0.0937345,0.5419157,5.9];
ybs[78]=['',0.0952796,0.5763954,5.79];
ybs[79]=['41 Psc',0.0944882,0.1449308,5.37];
ybs[80]=['',0.0958505,0.193567,6.56];
ybs[81]=['ρ And',0.0968922,0.6646603,5.18];
ybs[82]=['π Tuc',0.094195,-1.213202,5.51];
ybs[83]=['ι Scl',0.0983775,-0.5038428,5.18];
ybs[84]=['',0.0995108,-0.3480918,5.12];
ybs[85]=['42 Psc',0.1024799,0.2372959,6.23];
ybs[86]=['',0.0974457,-1.3493723,5.97];
ybs[87]=['9 Cet',0.1043008,-0.2111134,6.39];
ybs[88]=['',0.1057336,-0.539701,6.55];
ybs[89]=['',0.10962,0.67528,7.39];
ybs[90]=['',0.1107149,0.9099006,5.57];
ybs[91]=['12 Cas',0.1131704,1.0811362,5.4];
ybs[92]=['',0.1114558,-0.0367517,6.07];
ybs[93]=['',0.1144295,0.9278233,5.74];
ybs[94]=['44 Psc',0.1154353,0.0358337,5.77];
ybs[95]=['β Hyi',0.1159488,-1.3463606,2.8];
ybs[96]=['α Phe',0.1190562,-0.7364026,2.39];
ybs[97]=['κ Phe',0.1186973,-0.7603814,3.94];
ybs[98]=['10 Cet',0.1207507,0.0011102,6.19];
ybs[99]=['',0.1233484,-0.4439058,5.98];
ybs[100]=['47 Psc',0.1270486,0.3142692,5.06];
ybs[101]=['',0.1279985,0.7768054,5.17];
ybs[102]=['η Scl',0.126285,-0.5741081,4.81];
ybs[103]=['48 Psc',0.1277549,0.2889957,6.06];
ybs[104]=['',0.1282645,0.1798204,6.04];
ybs[105]=['',0.1281995,-0.3529365,6.43];
ybs[106]=['',0.1284679,-0.6946721,5.43];
ybs[107]=['',0.1310659,0.646002,6.26];
ybs[108]=['',0.129585,-0.8799876,6.26];
ybs[109]=['',0.1406761,1.346216,6.21];
ybs[110]=['',0.1373959,1.0487738,5.94];
ybs[111]=['28 And',0.1361673,0.5212386,5.23];
ybs[112]=['',0.1348254,-0.2574541,6.14];
ybs[113]=['',0.1345124,-0.5585671,6.57];
ybs[114]=['12 Cet',0.1356406,-0.0670924,5.72];
ybs[115]=['',0.1370173,-0.4132011,5.19];
ybs[116]=['',0.1372736,-0.7125543,6.19];
ybs[117]=['',0.1370829,-0.8395367,5.69];
ybs[118]=['13 Cas',0.1423276,1.1629558,6.18];
ybs[119]=['',0.1418934,0.5880833,5.87];
ybs[120]=['λ Cas',0.1436157,0.9535644,4.73];
ybs[121]=['',0.143213,0.9241945,5.6];
ybs[122]=['λ1 Phe',0.1413504,-0.849811,4.77];
ybs[123]=['β1 Tuc',0.1416835,-1.0968529,4.37];
ybs[124]=['β2 Tuc',0.1417486,-1.0969886,4.54];
ybs[125]=['',0.1464315,0.7610974,6.7];
ybs[126]=['',0.1508332,1.2408339,6.42];
ybs[127]=['κ Cas',0.1491448,1.1003354,4.16];
ybs[128]=['52 Psc',0.1468989,0.3561761,5.38];
ybs[129]=['51 Psc',0.1459766,0.1233687,5.67];
ybs[130]=['',0.1468699,0.4833427,6.67];
ybs[131]=['',0.1479373,0.4955549,6.3];
ybs[132]=['',0.1497471,0.9600689,5.93];
ybs[133]=['β3 Tuc',0.1468325,-1.0981293,5.09];
ybs[134]=['16 Cas',0.1554559,1.1669808,6.48];
ybs[135]=['',0.1513954,-0.5139204,5.55];
ybs[136]=['θ Tuc',0.1494079,-1.241858,6.13];
ybs[137]=['',0.1545675,-0.9121135,5.57];
ybs[138]=['',0.1570322,0.235338,6.4];
ybs[139]=['13 Cet',0.1583652,-0.0607381,5.2];
ybs[140]=['14 Cet',0.1596836,-0.0068563,5.93];
ybs[141]=['',0.1627086,0.947387,5.08];
ybs[142]=['',0.161353,0.2324666,6.41];
ybs[143]=['',0.1642132,1.0548553,5.79];
ybs[144]=['λ2 Phe',0.1599426,-0.8358054,5.51];
ybs[145]=['',0.1592881,-0.954856,6.06];
ybs[146]=['',0.1632833,0.4776508,6.5];
ybs[147]=['',0.1617979,-0.2593722,6.45];
ybs[148]=['',0.1620327,-0.3967102,6.06];
ybs[149]=['',0.1653709,0.7784383,5.13];
ybs[150]=['ζ Cas',0.1663501,0.9426444,3.66];
ybs[151]=['π And',0.1657277,0.5904807,4.36];
ybs[152]=['53 Psc',0.1651921,0.2678083,5.89];
ybs[153]=['',0.1666966,0.4210914,6.47];
ybs[154]=['',0.1677954,0.6198017,5.48];
ybs[155]=['',0.1808797,1.4417507,6.4];
ybs[156]=['',0.167382,-0.4303048,5.57];
ybs[157]=['',0.1636938,-1.1346749,6.42];
ybs[158]=['',0.1682638,0.0566855,6.39];
ybs[159]=['',0.1668882,-0.9473924,6.41];
ybs[160]=['ε And',0.1730027,0.5135481,4.37];
ybs[161]=['',0.1758751,0.8633597,5.43];
ybs[162]=['δ And',0.1763932,0.5405851,3.27];
ybs[163]=['54 Psc',0.1764739,0.372854,5.87];
ybs[164]=['55 Psc',0.1789351,0.3761305,5.36];
ybs[165]=['α Cas',0.181874,0.9887207,2.23];
ybs[166]=['',0.1722466,-1.2745223,6.85];
ybs[167]=['',0.1787295,-0.5907819,6.69];
ybs[168]=['',0.1781898,-0.7798882,6.01];
ybs[169]=['',0.1810925,-0.2863149,6.49];
ybs[170]=['',0.1813538,-0.4135058,6.14];
ybs[171]=['',0.182174,-0.073996,5.91];
ybs[172]=['32 And',0.1843014,0.6906418,5.33];
ybs[173]=['',0.1804394,-1.0357155,5.89];
ybs[174]=['',0.1889275,1.1564492,5.83];
ybs[175]=['',0.1862661,0.4318184,6.04];
ybs[176]=['ξ Cas',0.1885757,0.883567,4.8];
ybs[177]=['μ Phe',0.1845315,-0.8023761,4.59];
ybs[178]=['',0.190721,1.0273959,6.17];
ybs[179]=['ξ Phe',0.1863047,-0.9841819,5.7];
ybs[180]=['π Cas',0.1946599,0.8226915,4.94];
ybs[181]=['λ1 Scl',0.1906678,-0.6693552,6.06];
ybs[182]=['',0.1902327,-1.0498223,5.98];
ybs[183]=['ρ Tuc',0.1890991,-1.140676,5.39];
ybs[184]=['β Cet',0.1946585,-0.3119714,2.04];
ybs[185]=['',0.1989211,0.8373409,5.67];
ybs[186]=['',0.1957768,-0.2076884,6.02];
ybs[187]=['η Phe',0.1931571,-1.0009639,4.36];
ybs[188]=['21 Cas',0.2052679,1.3107398,5.66];
ybs[189]=['ο Cas',0.2001737,0.8446756,4.54];
ybs[190]=['φ1 Cet',0.1973277,-0.1832155,4.76];
ybs[191]=['λ2 Scl',0.1971451,-0.6686304,5.9];
ybs[192]=['',0.2027564,0.9657521,5.42];
ybs[193]=['',0.1996421,-0.3821259,5.24];
ybs[194]=['',0.200362,-0.7428954,5.94];
ybs[195]=['',0.1981752,-1.0888383,6.07];
ybs[196]=['',0.2092183,1.2118992,6.33];
ybs[197]=['',0.2026548,-0.0788421,6.15];
ybs[198]=['',0.2004003,-0.9355507,6.15];
ybs[199]=['18 Cet',0.2029307,-0.222861,6.15];
ybs[200]=['',0.2069807,0.9672097,6.52];
ybs[201]=['',0.2064923,0.7849295,6.05];
ybs[202]=['',0.2038494,-0.2847042,6.47];
ybs[203]=['',0.2090775,1.0417198,6.39];
ybs[204]=['23 Cas',0.2145632,1.3082828,5.41];
ybs[205]=['',0.2038146,-0.8279816,5.8];
ybs[206]=['',0.2059887,-0.3911314,5.5];
ybs[207]=['57 Psc',0.2078101,0.2720494,5.38];
ybs[208]=['',0.2160407,1.2703649,5.87];
ybs[209]=['58 Psc',0.2098572,0.210933,5.5];
ybs[210]=['59 Psc',0.2107967,0.3436649,6.13];
ybs[211]=['ζ And',0.2113236,0.4254915,4.06];
ybs[212]=['60 Psc',0.2114261,0.1195982,5.99];
ybs[213]=['61 Psc',0.2138061,0.3671626,6.54];
ybs[214]=['',0.2126743,-0.3132827,5.7];
ybs[215]=['η Cas',0.2195062,1.011022,3.44];
ybs[216]=['',0.2139376,-0.3771816,5.57];
ybs[217]=['62 Psc',0.2153437,0.1293559,5.93];
ybs[218]=['',0.2157359,0.0941098,5.75];
ybs[219]=['ν Cas',0.2181941,0.891511,4.89];
ybs[220]=['δ Psc',0.2170625,0.1343294,4.43];
ybs[221]=['64 Psc',0.2184233,0.2976141,5.07];
ybs[222]=['ν And',0.2223155,0.7189059,4.53];
ybs[223]=['',0.2201455,-0.2347461,5.59];
ybs[224]=['',0.2192065,-0.4193143,5.9];
ybs[225]=['',0.217696,-0.8130842,6.27];
ybs[226]=['65 Psc',0.222456,0.4855892,7];
ybs[227]=['65 Psc',0.222485,0.4855795,7.1];
ybs[228]=['',0.2206302,-0.4057934,6.28];
ybs[229]=['',0.2268411,1.1232725,5.39];
ybs[230]=['',0.2245209,0.7873799,6.15];
ybs[231]=['φ2 Cet',0.2232216,-0.1838371,5.19];
ybs[232]=['λ Hyi',0.2150375,-1.3057122,5.07];
ybs[233]=['',0.2291427,1.0806564,6.07];
ybs[234]=['',0.2274846,0.900927,6.39];
ybs[235]=['',0.2226085,-0.7554371,6.48];
ybs[236]=['',0.248468,1.4628995,5.62];
ybs[237]=['',0.2301322,0.9020263,6.21];
ybs[238]=['ρ Phe',0.2252036,-0.8879476,5.22];
ybs[239]=['',0.2284716,0.0610206,6.37];
ybs[240]=['',0.2369791,1.0687555,4.82];
ybs[241]=['',0.2304724,-0.7609287,6.9];
ybs[242]=['',0.2357282,0.6747382,6.69];
ybs[243]=['',0.2342263,-0.4170423,5.46];
ybs[244]=['20 Cet',0.2358682,-0.0180317,4.77];
ybs[245]=['',0.238248,0.6550051,6.06];
ybs[246]=['',0.2399179,0.9215355,6.27];
ybs[247]=['',0.2365294,-0.4305018,6.46];
ybs[248]=['λ1 Tuc',0.2320313,-1.2111421,6.22];
ybs[249]=['υ1 Cas',0.245372,1.0312027,4.83];
ybs[250]=['66 Psc',0.2429294,0.3368342,5.74];
ybs[251]=['21 Cet',0.2414118,-0.1506211,6.16];
ybs[252]=['',0.2454914,0.8515354,6.27];
ybs[253]=['',0.2376904,-1.0953759,5.7];
ybs[254]=['36 And',0.2446313,0.414326,5.47];
ybs[255]=['',0.2458562,0.4305327,6.2];
ybs[256]=['',0.2506526,1.0141638,6.21];
ybs[257]=['',0.254253,1.202299,6.37];
ybs[258]=['67 Psc',0.2490721,0.476826,6.09];
ybs[259]=['',0.2475887,-0.1263009,5.85];
ybs[260]=['γ Cas',0.2529026,1.0616356,2.47];
ybs[261]=['υ2 Cas',0.2526526,1.0348352,4.63];
ybs[262]=['',0.2532256,1.055459,5.55];
ybs[263]=['φ3 Cet',0.2489432,-0.1947088,5.31];
ybs[264]=['',0.2483482,-0.482843,6.1];
ybs[265]=['μ And',0.2526105,0.673872,3.87];
ybs[266]=['λ2 Tuc',0.243316,-1.21154,5.45];
ybs[267]=['η And',0.2544113,0.4106415,4.42];
ybs[268]=['',0.2566976,0.801982,6.12];
ybs[269]=['',0.2610835,1.1599906,5.97];
ybs[270]=['68 Psc',0.2572234,0.5079373,5.42];
ybs[271]=['',0.2590309,0.5944805,5.98];
ybs[272]=['',0.2573802,0.2409649,6.32];
ybs[273]=['',0.2592305,0.3755046,6.37];
ybs[274]=['',0.2701685,1.2408093,6.39];
ybs[275]=['φ4 Cet',0.2607477,-0.1966928,5.61];
ybs[276]=['α Scl',0.2600191,-0.5104589,4.31];
ybs[277]=['',0.2583876,-1.0574251,6.23];
ybs[278]=['',0.2671477,0.7822787,6.84];
ybs[279]=['',0.2671624,0.7823175,6.04];
ybs[280]=['',0.2656943,0.1150739,6.11];
ybs[281]=['',0.3138759,1.5041511,4.25];
ybs[282]=['',0.4659,1.5551897,6.46];
ybs[283]=['',0.2771174,0.8926462,6.47];
ybs[284]=['ξ Scl',0.2716505,-0.677304,5.59];
ybs[285]=['',0.2801771,0.8287848,6.45];
ybs[286]=['39 And',0.2795395,0.7235224,5.98];
ybs[287]=['σ Psc',0.2790191,0.5570085,5.5];
ybs[288]=['',0.2831643,1.067874,5.92];
ybs[289]=['σ Scl',0.2766994,-0.548768,5.5];
ybs[290]=['ε Psc',0.2793026,0.1396225,4.28];
ybs[291]=['ω Phe',0.2744173,-0.9929629,6.11];
ybs[292]=['25 Cet',0.2796172,-0.0825,5.43];
ybs[293]=['',0.2862944,1.076691,5.84];
ybs[294]=['',0.2847357,0.9182498,5.99];
ybs[295]=['',0.2781193,-0.8078726,5.36];
ybs[296]=['',0.2804445,-0.5134077,6.29];
ybs[297]=['26 Cet',0.2830508,0.0257667,6.04];
ybs[298]=['',0.2879344,0.8922039,6.54];
ybs[299]=['',0.2861625,0.5195527,6.19];
ybs[300]=['',0.2770473,-1.1405078,6.21];
ybs[301]=['',0.286954,0.6998885,6.72];
ybs[302]=['',0.350822,1.520959,6.25];
ybs[303]=['73 Psc',0.2877184,0.100634,6];
ybs[304]=['72 Psc',0.2887439,0.2627696,5.68];
ybs[305]=['',0.2953402,1.097305,6.54];
ybs[306]=['ψ1 Psc',0.2913987,0.3766897,5.34];
ybs[307]=['ψ1 Psc',0.2914568,0.3765491,5.56];
ybs[308]=['',0.3099741,1.3983664,6.29];
ybs[309]=['77 Psc',0.2918277,0.0875756,6.35];
ybs[310]=['77 Psc',0.2919878,0.0875949,7.25];
ybs[311]=['27 Cet',0.2907775,-0.1722597,6.12];
ybs[312]=['',0.2978334,0.995609,6.43];
ybs[313]=['28 Cet',0.2928362,-0.1698223,5.58];
ybs[314]=['',0.298404,0.9356276,6.38];
ybs[315]=['75 Psc',0.2951401,0.228034,6.12];
ybs[316]=['',0.2928703,-0.4168397,6.14];
ybs[317]=['μ Cas',0.3033246,0.9604423,5.17];
ybs[318]=['β Phe',0.292323,-0.813485,3.31];
ybs[319]=['',0.2940803,-0.6204914,6.61];
ybs[320]=['41 And',0.301923,0.7688349,5.03];
ybs[321]=['',0.2976222,-0.4169104,6.37];
ybs[322]=['',0.3046717,1.0187936,5.79];
ybs[323]=['78 Psc',0.3017529,0.5606221,6.25];
ybs[324]=['ψ2 Psc',0.3013096,0.3638703,5.55];
ybs[325]=['30 Cet',0.3001866,-0.168886,5.82];
ybs[326]=['80 Psc',0.3029638,0.1005089,5.52];
ybs[327]=['υ Phe',0.2998945,-0.7221796,5.21];
ybs[328]=['ι Tuc',0.2972089,-1.0762764,5.37];
ybs[329]=['',0.3234093,1.3924629,5.64];
ybs[330]=['η Cet',0.3037591,-0.1758112,3.45];
ybs[331]=['φ And',0.3084999,0.826427,4.25];
ybs[332]=['31 Cas',0.3144511,1.2023093,5.29];
ybs[333]=['β And',0.3092781,0.6235949,2.06];
ybs[334]=['ζ Phe',0.302123,-0.9623189,3.92];
ybs[335]=['ψ3 Psc',0.3094474,0.3450063,5.55];
ybs[336]=['44 And',0.3119329,0.7363562,5.65];
ybs[337]=['',0.3117172,0.4462194,5.8];
ybs[338]=['',0.3175203,1.1224441,5.55];
ybs[339]=['θ Cas',0.3157138,0.9644394,4.33];
ybs[340]=['',0.3110237,0.2754636,6.06];
ybs[341]=['32 Cas',0.318722,1.1366871,5.57];
ybs[342]=['32 Cet',0.3107955,-0.1535431,6.4];
ybs[343]=['33 Cet',0.3124871,0.0445799,5.95];
ybs[344]=['45 And',0.3156049,0.660306,5.81];
ybs[345]=['82 Psc',0.3152502,0.5503602,5.16];
ybs[346]=['',0.3096029,-1.0050597,6.41];
ybs[347]=['χ Psc',0.3165958,0.3690196,4.66];
ybs[348]=['τ Psc',0.3176187,0.5270585,4.51];
ybs[349]=['34 Cet',0.3175201,-0.0373955,5.94];
ybs[350]=['',0.325001,1.0788595,6.41];
ybs[351]=['',0.32185,0.7931848,6.11];
ybs[352]=['',0.3234352,0.5266089,6.19];
ybs[353]=['',0.3421618,1.3965718,6.26];
ybs[354]=['',0.320075,-0.5357081,6.52];
ybs[355]=['',0.3215599,-0.6588275,5.92];
ybs[356]=['φ Psc',0.3266633,0.4309531,4.65];
ybs[357]=['ζ Psc',0.3263845,0.1341018,5.24];
ybs[358]=['ζ Psc',0.3264864,0.1341551,6.3];
ybs[359]=['',0.3281812,0.4998248,6.43];
ybs[360]=['87 Psc',0.3282081,0.2834717,5.98];
ybs[361]=['',0.339073,1.2540478,7.83];
ybs[362]=['37 Cet',0.3291272,-0.1363969,5.13];
ybs[363]=['88 Psc',0.3306332,0.1239761,6.03];
ybs[364]=['38 Cet',0.3310379,-0.0151123,5.7];
ybs[365]=['',0.3387117,0.8410736,6.61];
ybs[366]=['ν Phe',0.3319845,-0.7927882,4.96];
ybs[367]=['',0.3379928,0.5798418,6.02];
ybs[368]=['',0.3415926,0.7855654,6.34];
ybs[369]=['39 Cet',0.3388087,-0.0417579,5.41];
ybs[370]=['',0.3427251,0.5559277,6.73];
ybs[371]=['',0.3582182,1.3557298,6.31];
ybs[372]=['',0.3463974,0.8295057,6.25];
ybs[373]=['κ Tuc',0.3335101,-1.2002318,4.86];
ybs[374]=['89 Psc',0.3440937,0.0649607,5.16];
ybs[375]=['',0.3488592,0.6543844,6.46];
ybs[376]=['',0.3393152,-1.1569853,6.24];
ybs[377]=['',0.3651994,1.3324828,6.38];
ybs[378]=['φ Cas',0.3551188,1.0182039,4.98];
ybs[379]=['υ Psc',0.3516752,0.4777212,4.76];
ybs[380]=['35 Cas',0.3598645,1.1303673,6.34];
ybs[381]=['42 Cet',0.352794,-0.007011,5.87];
ybs[382]=['',0.3737714,1.3758824,6.07];
ybs[383]=['',0.3561207,-0.0548014,6.23];
ybs[384]=['',0.3555357,-0.1942867,6.15];
ybs[385]=['91 Psc',0.3589337,0.5034404,5.23];
ybs[386]=['ξ And',0.3645802,0.7964919,4.88];
ybs[387]=['',0.3694624,1.0166474,6.45];
ybs[388]=['',0.3650896,0.0319935,6.2];
ybs[389]=['43 Cet',0.3649028,-0.0059867,6.49];
ybs[390]=['',0.3643515,-0.3311704,6.35];
ybs[391]=['47 And',0.3702475,0.6601097,5.58];
ybs[392]=['',0.3699554,0.5995615,6.29];
ybs[393]=['',0.3688182,0.3591092,5.97];
ybs[394]=['',0.3809854,1.2406863,6.49];
ybs[395]=['ψ Cas',0.3813765,1.190944,4.74];
ybs[396]=['',0.3685676,-0.5382421,5.84];
ybs[397]=['44 Cet',0.371185,-0.1378993,6.21];
ybs[398]=['θ Cet',0.3711028,-0.1409681,3.6];
ybs[399]=['δ Cas',0.3803165,1.0531557,2.68];
ybs[400]=['',0.3725149,-0.1188278,5.91];
ybs[401]=['',0.3737959,-0.2714674,6.14];
ybs[402]=['',0.3746098,-0.0478622,6.15];
ybs[403]=['',0.378382,0.4122088,6.18];
ybs[404]=['',0.3734323,-0.7223246,5.42];
ybs[405]=['',0.3818829,0.7603317,5.96];
ybs[406]=['',0.3809712,0.605381,6.31];
ybs[407]=['',0.3734407,-0.77531,6.26];
ybs[408]=['46 Cet',0.3779825,-0.2529458,4.9];
ybs[409]=['ρ Psc',0.3811997,0.3364691,5.38];
ybs[410]=['94 Psc',0.3831218,0.3376555,5.5];
ybs[411]=['',0.3851542,0.6018485,6.27];
ybs[412]=['',0.3818111,-0.0051069,6.41];
ybs[413]=['ω And',0.3878193,0.7943418,4.83];
ybs[414]=['',0.3867804,0.7191867,6.46];
ybs[415]=['',0.3837673,0.063551,6.58];
ybs[416]=['',0.3743437,-1.1216036,5.93];
ybs[417]=['47 Cet',0.3834115,-0.2260329,5.66];
ybs[418]=['',0.3882489,0.7058339,6.6];
ybs[419]=['',0.383583,-0.5661348,5.79];
ybs[420]=['α UMi',0.7755275,1.5593261,2.02];
ybs[421]=['',0.387439,-0.188424,6.13];
ybs[422]=['',0.3903283,0.1407963,6.2];
ybs[423]=['38 Cas',0.4048263,1.2281844,5.81];
ybs[424]=['',0.4028401,1.1554636,6.14];
ybs[425]=['γ Phe',0.389441,-0.7542033,3.41];
ybs[426]=['49 And',0.3985598,0.8222687,5.27];
ybs[427]=['',0.3912015,-0.5874432,6.58];
ybs[428]=['97 Psc',0.3970221,0.3222038,6.02];
ybs[429]=['48 Cet',0.3952437,-0.375665,5.12];
ybs[430]=['μ Psc',0.3981741,0.109069,4.84];
ybs[431]=['',0.3943224,-0.8142124,6.31];
ybs[432]=['',0.3985711,-0.4555746,5.93];
ybs[433]=['η Psc',0.4039702,0.2696688,3.62];
ybs[434]=['',0.4071088,0.6092057,6.39];
ybs[435]=['',0.4135411,1.0198333,5.7];
ybs[436]=['δ Phe',0.401852,-0.8546468,3.95];
ybs[437]=['',0.4043413,-0.526711,5.82];
ybs[438]=['χ Cas',0.4157837,1.0356171,4.71];
ybs[439]=['',0.4036969,-0.7936102,6.17];
ybs[440]=['',0.4105033,-0.1555084,6.59];
ybs[441]=['',0.4094918,-0.6415917,5.51];
ybs[442]=['',0.4165542,0.6517358,5.88];
ybs[443]=['',0.4077233,-0.8660834,6.28];
ybs[444]=['',0.4133959,-0.1207883,5.76];
ybs[445]=['',0.4324189,1.2986059,6.58];
ybs[446]=['',0.4185771,0.3240194,5.89];
ybs[447]=['49 Cet',0.4172611,-0.271777,5.63];
ybs[448]=['',0.4236289,0.7187362,6.38];
ybs[449]=['',0.4179241,-0.5548022,6.12];
ybs[450]=['',0.426376,0.8521888,5.92];
ybs[451]=['101 Psc',0.4226879,0.2577081,6.22];
ybs[452]=['40 Cas',0.4371899,1.2765962,5.28];
ybs[453]=['',0.4233417,0.306092,5.8];
ybs[454]=['υ And',0.4276651,0.724478,4.09];
ybs[455]=['50 Cet',0.4231673,-0.2669676,5.42];
ybs[456]=['',0.4189093,-1.0129037,6.01];
ybs[457]=['',0.4340765,1.013708,5.56];
ybs[458]=['τ Scl',0.4236065,-0.5201667,5.69];
ybs[459]=['π Psc',0.4284337,0.213726,5.57];
ybs[460]=['51 And',0.4331063,0.8505348,3.57];
ybs[461]=['',0.4353438,0.7941879,6.36];
ybs[462]=['',0.4304322,-0.1623167,6.24];
ybs[463]=['',0.4093471,-1.3683383,6.11];
ybs[464]=['',0.4253994,-1.0152022,6.18];
ybs[465]=['χ And',0.4389059,0.7764892,4.98];
ybs[466]=['',0.4430342,0.9419819,6.39];
ybs[467]=['',0.4335667,-0.6357304,5.94];
ybs[468]=['α Eri',0.4296644,-0.9971561,0.46];
ybs[469]=['',0.4356322,-0.3695159,5.58];
ybs[470]=['',0.4354367,-0.4349075,6.7];
ybs[471]=['105 Psc',0.4397697,0.2881402,5.97];
ybs[472]=['',0.4446153,0.7574894,5.61];
ybs[473]=['τ And',0.4441777,0.7100022,4.94];
ybs[474]=['43 Cas',0.4532886,1.1893691,5.59];
ybs[475]=['',0.4345708,-0.9308764,6.84];
ybs[476]=['42 Cas',0.4561926,1.2343865,5.18];
ybs[477]=['',0.451425,1.067115,6.71];
ybs[478]=['',0.4523598,1.025042,6.37];
ybs[479]=['',0.4495058,0.7455442,4.95];
ybs[480]=['',0.4470322,0.4511479,6.17];
ybs[481]=['',0.4486289,0.52622,5.99];
ybs[482]=['',0.4387435,-0.9790364,5.87];
ybs[483]=['',0.4387726,-0.9789783,5.76];
ybs[484]=['',0.4554922,1.0738018,6.34];
ybs[485]=['ν Psc',0.4472467,0.0975729,4.44];
ybs[486]=['',0.4505087,0.6169464,5.64];
ybs[487]=['44 Cas',0.4569896,1.0586065,5.78];
ybs[488]=['',0.4483594,-0.1958568,5.75];
ybs[489]=['107 Psc',0.4521312,0.3555479,5.24];
ybs[490]=['',0.446594,-0.6637491,6.17];
ybs[491]=['',0.4560905,0.7928127,6.34];
ybs[492]=['φ Per',0.4579578,0.8864724,4.07];
ybs[493]=['π Scl',0.4497242,-0.562416,5.25];
ybs[494]=['',0.4492175,-0.6410524,5.72];
ybs[495]=['',0.4610591,1.005986,6.21];
ybs[496]=['',0.4527517,-0.0626143,4.99];
ybs[497]=['',0.4472465,-0.8715459,6.64];
ybs[498]=['',0.4630919,0.9981787,6.25];
ybs[499]=['',0.4581957,0.5636393,6.34];
ybs[500]=['',0.4612161,0.8070763,6.35];
ybs[501]=['',0.4472381,-1.0591787,5.71];
ybs[502]=['',0.4505978,-0.936155,5.52];
ybs[503]=['',0.4579197,-0.0813859,6.19];
ybs[504]=['109 Psc',0.4627523,0.3523,6.27];
ybs[505]=['τ Cet',0.4584185,-0.2763736,3.5];
ybs[506]=['ο Psc',0.4645937,0.1616163,4.26];
ybs[507]=['',0.4765718,1.1162044,5.63];
ybs[508]=['',0.4252361,-1.4463732,5.87];
ybs[509]=['',0.4669569,-0.0982848,5.34];
ybs[510]=['ε Scl',0.46514,-0.4354664,5.31];
ybs[511]=['',0.4699462,0.3056934,6.55];
ybs[512]=['τ1 Hyi',0.4424027,-1.3795988,6.33];
ybs[513]=['',0.4667088,-0.4755523,6.39];
ybs[514]=['',0.4758998,0.808634,6.32];
ybs[515]=['',0.4664253,-0.8851326,5.49];
ybs[516]=['',0.4663498,-0.9323535,5.04];
ybs[517]=['',0.4793664,0.6641707,5.94];
ybs[518]=['4 Ari',0.4768944,0.2977022,5.84];
ybs[519]=['',0.4794364,0.5723227,5.79];
ybs[520]=['',0.4718777,-0.7270735,6.18];
ybs[521]=['',0.4210161,-1.4776945,5.69];
ybs[522]=['',0.4823472,0.8377267,5.82];
ybs[523]=['',0.4777737,0.0660961,5.91];
ybs[524]=['',0.4742506,-0.6467856,6.32];
ybs[525]=['',0.4898949,0.908168,5.9];
ybs[526]=['1 Ari',0.4855481,0.3905409,5.86];
ybs[527]=['χ Cet',0.4825664,-0.1847462,4.67];
ybs[528]=['',0.4810526,-0.5405546,6.34];
ybs[529]=['1 Per',0.4945785,0.9642614,5.52];
ybs[530]=['',0.488514,0.1945037,5.94];
ybs[531]=['',0.4830062,-0.6685085,6.37];
ybs[532]=['2 Per',0.4951163,0.8882566,5.79];
ybs[533]=['',0.4849995,-0.8327895,6.14];
ybs[534]=['',0.4981622,0.900156,6.26];
ybs[535]=['ζ Cet',0.4907508,-0.1786209,3.73];
ybs[536]=['',0.5025631,0.9721176,6.45];
ybs[537]=['',0.4873991,-0.8745003,5.94];
ybs[538]=['ε Cas',0.5056692,1.1129969,3.38];
ybs[539]=['55 And',0.4997173,0.7126187,5.4];
ybs[540]=['α Tri',0.4985355,0.5180009,3.41];
ybs[541]=['γ1 Ari',0.5002865,0.338526,4.83];
ybs[542]=['γ2 Ari',0.5002865,0.3384872,4.75];
ybs[543]=['',0.4967785,-0.2937165,5.8];
ybs[544]=['ω Cas',0.5132272,1.200523,4.99];
ybs[545]=['ξ Psc',0.5001154,0.0573826,4.62];
ybs[546]=['τ2 Hyi',0.4696594,-1.3975708,6.06];
ybs[547]=['',0.5067447,0.7121271,6.24];
ybs[548]=['',0.5069148,0.6497557,6.26];
ybs[549]=['β Ari',0.5051615,0.3649146,2.64];
ybs[550]=['',0.498571,-0.6718537,6.1];
ybs[551]=['ψ Phe',0.4994681,-0.8063808,4.41];
ybs[552]=['',0.511062,0.65236,5.89];
ybs[553]=['56 And',0.5121464,0.6519031,5.67];
ybs[554]=['φ Phe',0.5027283,-0.7399643,5.11];
ybs[555]=['7 Ari',0.5104993,0.4132406,5.74];
ybs[556]=['',0.5103112,0.0340241,6.01];
ybs[557]=['',0.5237217,1.0785623,6.02];
ybs[558]=['',0.5200738,0.7294366,6.78];
ybs[559]=['ι Ari',0.5169436,0.3127081,5.1];
ybs[560]=['',0.518798,0.4870114,5.82];
ybs[561]=['56 Cet',0.5132502,-0.3914321,4.85];
ybs[562]=['χ Eri',0.5093282,-0.8990043,3.7];
ybs[563]=['',0.5286939,1.1295793,5.26];
ybs[564]=['3 Per',0.5230492,0.8605032,5.69];
ybs[565]=['λ Ari',0.519575,0.4135613,4.79];
ybs[566]=['η2 Hyi',0.5037552,-1.1789212,4.69];
ybs[567]=['',0.5080108,-1.0604898,6.06];
ybs[568]=['',0.5458828,1.3616048,6.04];
ybs[569]=['',0.5138637,-0.9017528,6.1];
ybs[570]=['',0.5147597,-0.8252888,4.83];
ybs[571]=['48 Cas',0.5396673,1.2392723,4.54];
ybs[572]=['',0.5207564,-0.5753924,6.35];
ybs[573]=['',0.526801,0.3692616,5.87];
ybs[574]=['',0.5259214,0.2163085,6.09];
ybs[575]=['',0.5455938,1.2906422,6.23];
ybs[576]=['50 Cas',0.546427,1.2656976,3.98];
ybs[577]=['47 Cas',0.5551735,1.3505123,5.38];
ybs[578]=['112 Psc',0.5289073,0.0557788,5.88];
ybs[579]=['57 Cet',0.5267938,-0.3617311,5.41];
ybs[580]=['',0.5168483,-1.1401439,6.37];
ybs[581]=['υ Cet',0.5278208,-0.3661537,4];
ybs[582]=['52 Cas',0.542932,1.1344519,6];
ybs[583]=['',0.52999,-0.1470442,5.51];
ybs[584]=['',0.5257449,-0.7318468,5.57];
ybs[585]=['53 Cas',0.5434365,1.1255259,5.58];
ybs[586]=['4 Per',0.5396558,0.9526982,5.04];
ybs[587]=['α Hyi',0.5209931,-1.0728654,2.86];
ybs[588]=['49 Cas',0.5565104,1.3301534,5.22];
ybs[589]=['σ Hyi',0.5053826,-1.3656935,6.16];
ybs[590]=['π For',0.5330342,-0.5219103,5.35];
ybs[591]=['α Psc',0.5371633,0.0499479,3.82];
ybs[592]=['α Psc',0.5371633,0.0499479,4.33];
ybs[593]=['',0.5763597,1.4205564,6.05];
ybs[594]=['',0.5507914,1.1379684,6.52];
ybs[595]=['ε Tri',0.5418004,0.5826227,5.5];
ybs[596]=['',0.5245708,-1.1513507,6.1];
ybs[597]=['',0.5397149,0.2369235,5.94];
ybs[598]=['χ Phe',0.5346299,-0.7786839,5.14];
ybs[599]=['γ1 And',0.5461399,0.740498,2.26];
ybs[600]=['γ2 And',0.5461909,0.7405174,4.84];
ybs[601]=['10 Ari',0.544632,0.4543672,5.63];
ybs[602]=['',0.5383747,-0.5160397,6.42];
ybs[603]=['60 Cet',0.5421272,0.0039485,5.43];
ybs[604]=['',0.5408926,-0.2654274,5.86];
ybs[605]=['',0.5447111,0.320287,6.21];
ybs[606]=['61 Cet',0.5447731,-0.0042331,5.93];
ybs[607]=['',0.5441456,-0.0699149,5.62];
ybs[608]=['ν For',0.5471961,-0.5096249,4.69];
ybs[609]=['κ Ari',0.5572673,0.3969812,5.03];
ybs[610]=['',0.5554102,0.145641,6.31];
ybs[611]=['11 Ari',0.5584481,0.450324,6.15];
ybs[612]=['',0.5564955,0.0023047,6.28];
ybs[613]=['α Ari',0.5599412,0.4111883,2];
ybs[614]=['',0.5677724,1.0213672,5.67];
ybs[615]=['',0.5665795,0.7776474,6.42];
ybs[616]=['58 And',0.5660488,0.6624512,4.82];
ybs[617]=['',0.5738591,0.9414147,6.31];
ybs[618]=['β Tri',0.570575,0.6123215,3];
ybs[619]=['14 Ari',0.5698148,0.4544135,4.98];
ybs[620]=['',0.5694657,0.3023035,6.43];
ybs[621]=['',0.5660728,-0.3086263,6.1];
ybs[622]=['',0.5903095,1.293688,6.29];
ybs[623]=['5 Per',0.5800068,1.0077792,6.36];
ybs[624]=['59 And',0.5765345,0.6830397,5.63];
ybs[625]=['59 And',0.5766002,0.6831027,6.1];
ybs[626]=['',0.5695005,-0.4232353,6.48];
ybs[627]=['15 Ari',0.574934,0.3420185,5.7];
ybs[628]=['',0.5671179,-0.7578271,5.85];
ybs[629]=['16 Ari',0.5775805,0.4543567,6.02];
ybs[630]=['5 Tri',0.5786647,0.5519098,6.23];
ybs[631]=['64 Cet',0.5778785,0.151241,5.63];
ybs[632]=['',0.5711441,-0.763048,6.32];
ybs[633]=['',0.5723598,-0.8853777,6.12];
ybs[634]=['',0.577603,-0.1737732,6.01];
ybs[635]=['63 Cet',0.5787495,-0.0301871,5.93];
ybs[636]=['55 Cas',0.5939449,1.1627254,6.07];
ybs[637]=['',0.5897311,1.0237432,6.44];
ybs[638]=['6 Tri',0.5828042,0.5305541,4.94];
ybs[639]=['60 And',0.5869448,0.7736499,4.83];
ybs[640]=['',0.583762,0.4234721,5.96];
ybs[641]=['',0.588912,0.8929266,5.31];
ybs[642]=['η Ari',0.584467,0.3718629,5.27];
ybs[643]=['',0.5906656,0.8304127,6.06];
ybs[644]=['19 Ari',0.585447,0.2683444,5.71];
ybs[645]=['ξ1 Cet',0.585078,0.1560667,4.37];
ybs[646]=['66 Cet',0.5839526,-0.0401121,5.54];
ybs[647]=['',0.5845523,-0.3648605,5.86];
ybs[648]=['μ For',0.583857,-0.5345689,5.28];
ybs[649]=['',0.5990796,0.8361145,6.33];
ybs[650]=['',0.6034957,0.9974461,6.48];
ybs[651]=['7 Tri',0.5984573,0.5838711,5.28];
ybs[652]=['20 Ari',0.5975187,0.4516489,5.79];
ybs[653]=['21 Ari',0.5972681,0.4387337,5.58];
ybs[654]=['',0.5955064,-0.1635536,6.55];
ybs[655]=['',0.5906311,-0.7168374,5.91];
ybs[656]=['δ Tri',0.6033596,0.5989676,4.87];
ybs[657]=['8 Per',0.6085259,1.0121788,5.75];
ybs[658]=['7 Per',0.6088346,1.0054928,5.98];
ybs[659]=['',0.6058932,0.7749427,6.7];
ybs[660]=['γ Tri',0.6044918,0.5923875,4.01];
ybs[661]=['',0.6036132,0.4164688,6.55];
ybs[662]=['67 Cet',0.6021613,-0.1104448,5.51];
ybs[663]=['π1 Hyi',0.5876373,-1.1824017,5.55];
ybs[664]=['',0.6187858,1.124523,6.6];
ybs[665]=['θ Ari',0.6076797,0.3489781,5.62];
ybs[666]=['62 And',0.6135513,0.828569,5.3];
ybs[667]=['',0.6130827,0.8127306,6.21];
ybs[668]=['',0.6068605,0.0323179,5.58];
ybs[669]=['',0.614052,0.8560622,6.37];
ybs[670]=['φ Eri',0.5988152,-0.8974108,3.56];
ybs[671]=['10 Tri',0.6114916,0.5015399,5.03];
ybs[672]=['',0.6114296,0.4059879,6.46];
ybs[673]=['',0.614755,0.6968823,6.63];
ybs[674]=['π2 Hyi',0.5929953,-1.1807446,5.69];
ybs[675]=['',0.6197115,0.8273547,6.11];
ybs[676]=['',0.6164383,0.5285141,6.47];
ybs[677]=['ο Cet',0.6125333,-0.0503348,3.04];
ybs[678]=['63 And',0.6210674,0.8769303,5.59];
ybs[679]=['',0.6104244,-0.4512009,6.34];
ybs[680]=['',0.6139675,-0.0742136,6.5];
ybs[681]=['9 Per',0.6274517,0.976305,5.17];
ybs[682]=['',0.6118621,-0.7287586,6.37];
ybs[683]=['',0.6288648,0.7241175,5.82];
ybs[684]=['',0.6133453,-0.9747891,5.81];
ybs[685]=['69 Cet',0.623938,0.0085279,5.28];
ybs[686]=['',0.6340105,0.9679003,6.28];
ybs[687]=['70 Cet',0.625061,-0.0138282,5.42];
ybs[688]=['',0.6240582,-0.1864887,5.46];
ybs[689]=['',0.6241714,-0.3066451,5.87];
ybs[690]=['64 And',0.6361209,0.8743869,5.19];
ybs[691]=['κ For',0.6260344,-0.4140578,5.2];
ybs[692]=['10 Per',0.6402306,0.9896321,6.25];
ybs[693]=['',0.6279935,-0.3187312,6.22];
ybs[694]=['',0.623946,-0.7523636,6.31];
ybs[695]=['65 And',0.6414172,0.8791269,4.71];
ybs[696]=['',0.6281155,-0.6542178,6.53];
ybs[697]=['',0.6267051,-0.8901123,5.92];
ybs[698]=['ξ Ari',0.6366885,0.1867936,5.47];
ybs[699]=['',0.6337959,-0.4495166,6.44];
ybs[700]=['71 Cet',0.6370935,-0.0469164,6.33];
ybs[701]=['δ Hyi',0.6201258,-1.1967114,4.09];
ybs[702]=['',0.6343446,-0.7111957,6.18];
ybs[703]=['ι Cas',0.6579292,1.1779762,4.52];
ybs[704]=['ρ Cet',0.6411536,-0.2129118,4.89];
ybs[705]=['66 And',0.6512248,0.8841961,6.12];
ybs[706]=['',0.6413265,-0.2661542,5.83];
ybs[707]=['',0.6471134,0.4730639,6.18];
ybs[708]=['11 Tri',0.6487523,0.5566293,5.54];
ybs[709]=['',0.6437556,-0.3482168,5.88];
ybs[710]=['λ Hor',0.6347582,-1.0510366,5.35];
ybs[711]=['κ Hyi',0.6240269,-1.2837454,5.01];
ybs[712]=['',0.6583089,0.9708723,6.51];
ybs[713]=['12 Tri',0.6517707,0.5194163,5.29];
ybs[714]=['ξ2 Cet',0.6512271,0.149242,4.28];
ybs[715]=['',0.6503984,0.0358111,6.45];
ybs[716]=['13 Tri',0.6545804,0.5239944,5.89];
ybs[717]=['κ Eri',0.6446169,-0.8309957,4.25];
ybs[718]=['',0.6364661,-1.1589488,6.41];
ybs[719]=['',0.6562357,0.4111907,6.19];
ybs[720]=['φ For',0.6496763,-0.5885268,5.14];
ybs[721]=['',0.6574921,0.1685349,6.07];
ybs[722]=['',0.6611085,0.5920884,6.25];
ybs[723]=['',0.6522044,-0.5412556,6.11];
ybs[724]=['',0.6620155,0.4420081,5.92];
ybs[725]=['26 Ari',0.6623174,0.3481139,6.15];
ybs[726]=['',0.6582415,-0.394311,6.77];
ybs[727]=['27 Ari',0.6634303,0.3105637,6.23];
ybs[728]=['',0.6623763,0.0060289,6];
ybs[729]=['',0.6595151,-0.4380088,6.51];
ybs[730]=['',0.6481859,-1.1206527,6.37];
ybs[731]=['',0.6609631,-0.3919193,6.1];
ybs[732]=['14 Tri',0.6691614,0.6324537,5.15];
ybs[733]=['',0.6656847,0.04114,5.25];
ybs[734]=['',0.6724801,0.604442,5.83];
ybs[735]=['75 Cet',0.6684719,-0.0164982,5.35];
ybs[736]=['σ Cet',0.6678535,-0.264504,4.75];
ybs[737]=['29 Ari',0.672077,0.2639672,6.04];
ybs[738]=['',0.6679789,-0.6342137,6.3];
ybs[739]=['',0.6982662,1.272451,5.16];
ybs[740]=['λ1 For',0.6718285,-0.6031952,5.9];
ybs[741]=['',0.6746457,-0.3475417,6.21];
ybs[742]=['',0.6839701,0.6938226,6.36];
ybs[743]=['',0.695076,1.1490109,5.78];
ybs[744]=['',0.6846795,0.6527676,5.71];
ybs[745]=['ω For',0.6751942,-0.4911928,4.9];
ybs[746]=['15 Tri',0.6851739,0.6069569,5.35];
ybs[747]=['',0.6813627,0.1319504,6.18];
ybs[748]=['77 Cet',0.6794701,-0.1356209,5.75];
ybs[749]=['',0.685772,0.1217445,5.82];
ybs[750]=['ν Cet',0.6848415,0.0991679,4.86];
ybs[751]=['',0.6746074,-0.8901942,6.24];
ybs[752]=['',0.6904382,0.6775538,5.9];
ybs[753]=['',0.6891609,0.5531956,6.1];
ybs[754]=['',0.690671,0.5995566,5.3];
ybs[755]=['80 Cet',0.6850909,-0.1351431,5.53];
ybs[756]=['',0.6921997,0.6978507,6.54];
ybs[757]=['',0.6909113,0.5756112,6.25];
ybs[758]=['',0.6722325,-1.0907883,6.77];
ybs[759]=['31 Ari',0.6882985,0.2187914,5.68];
ybs[760]=['30 Ari',0.6900401,0.4317341,7.09];
ybs[761]=['30 Ari',0.6902438,0.4317192,6.5];
ybs[762]=['',0.6879876,0.1364509,5.81];
ybs[763]=['ι1 For',0.6852131,-0.5228345,5.75];
ybs[764]=['',0.6962706,0.6599863,6.18];
ybs[765]=['',0.6970116,0.6663171,6.3];
ybs[766]=['',0.6942203,0.1358418,6.39];
ybs[767]=['81 Cet',0.6925916,-0.0577375,5.65];
ybs[768]=['λ2 For',0.6886552,-0.6019655,5.79];
ybs[769]=['ν Ari',0.6980687,0.3848277,5.43];
ybs[770]=['',0.7456872,1.4230141,5.78];
ybs[771]=['',0.6967506,0.0616234,6.21];
ybs[772]=['μ Hyi',0.6600629,-1.3791478,5.28];
ybs[773]=['ι2 For',0.6946106,-0.5254547,5.83];
ybs[774]=['η Hor',0.6897516,-0.9155109,5.31];
ybs[775]=['δ Cet',0.7004701,0.0072612,4.07];
ybs[776]=['',0.6948007,-0.6615278,6.49];
ybs[777]=['ε Cet',0.7005426,-0.2056838,4.84];
ybs[778]=['33 Ari',0.7063647,0.4738193,5.3];
ybs[779]=['',0.7039965,0.108195,6.25];
ybs[780]=['',0.7033977,-0.1634649,5.78];
ybs[781]=['11 Per',0.7178844,0.9632827,5.77];
ybs[782]=['',0.7021391,-0.533139,6.52];
ybs[783]=['',0.7175524,0.9357118,5.84];
ybs[784]=['12 Per',0.713626,0.7030253,4.91];
ybs[785]=['',0.7006546,-0.7470759,4.75];
ybs[786]=['84 Cet',0.7080756,-0.0106285,5.71];
ybs[787]=['',0.7270246,1.1852582,5.95];
ybs[788]=['',0.7174094,0.8438978,6.48];
ybs[789]=['μ Ari',0.7135172,0.3507788,5.69];
ybs[790]=['ι Eri',0.7045519,-0.6940908,4.11];
ybs[791]=['',0.7105202,-0.0545703,6.05];
ybs[792]=['',0.7092172,-0.2524213,5.98];
ybs[793]=['',0.7137922,0.1889863,6.3];
ybs[794]=['',0.6980085,-1.1204043,6.55];
ybs[795]=['θ Per',0.7225668,0.8606947,4.12];
ybs[796]=['14 Per',0.7218295,0.7746265,5.43];
ybs[797]=['35 Ari',0.7184646,0.4850853,4.66];
ybs[798]=['ζ Hor',0.7037914,-0.9505567,5.21];
ybs[799]=['',0.7201593,0.4489693,6.35];
ybs[800]=['γ Cet',0.7171899,0.0579802,3.47];
ybs[801]=['',0.7108969,-0.6684134,6.01];
ybs[802]=['ε Hyi',0.6977264,-1.1899557,4.11];
ybs[803]=['',0.7106951,-0.8104929,6.1];
ybs[804]=['36 Ari',0.7219801,0.3115365,6.46];
ybs[805]=['ο Ari',0.7229228,0.268736,5.77];
ybs[806]=['ι Hor',0.7122879,-0.8851225,5.41];
ybs[807]=['π Cet',0.7203852,-0.2403786,4.25];
ybs[808]=['38 Ari',0.7246499,0.2187154,5.18];
ybs[809]=['μ Cet',0.7245141,0.1780203,4.27];
ybs[810]=['',0.7161683,-0.7058335,6.36];
ybs[811]=['',0.7452899,1.2168149,6.18];
ybs[812]=['',0.7261624,0.0837266,6.03];
ybs[813]=['',0.7208251,-0.5661694,6.22];
ybs[814]=['τ1 Eri',0.7245422,-0.322657,4.47];
ybs[815]=['',0.7341032,0.6295151,6.25];
ybs[816]=['',0.7344667,0.6220339,6.3];
ybs[817]=['',0.7192384,-0.9160289,6.15];
ybs[818]=['',0.7243566,-0.8063705,6.85];
ybs[819]=['',0.7147161,-1.1628812,6.26];
ybs[820]=['39 Ari',0.7379723,0.5119376,4.51];
ybs[821]=['',0.7463109,0.9977737,6.25];
ybs[822]=['',0.7316575,-0.3761996,6.49];
ybs[823]=['',0.7335242,-0.3909647,6.47];
ybs[824]=['40 Ari',0.740401,0.320583,5.82];
ybs[825]=['',0.7585226,1.2037845,5.8];
ybs[826]=['',0.7415928,0.4410868,5.86];
ybs[827]=['',0.744985,0.6529315,6.45];
ybs[828]=['',0.7370377,-0.2159999,6.9];
ybs[829]=['γ Hor',0.7238708,-1.1103586,5.74];
ybs[830]=['η Per',0.7513933,0.9770216,3.76];
ybs[831]=['η1 For',0.7347592,-0.6189987,6.51];
ybs[832]=['π Ari',0.7436828,0.3062766,5.22];
ybs[833]=['ζ Hyi',0.7237101,-1.1786399,4.84];
ybs[834]=['41 Ari',0.7469732,0.4772516,3.63];
ybs[835]=['',0.7562315,1.0192375,6.45];
ybs[836]=['16 Per',0.7499733,0.6702472,4.23];
ybs[837]=['β For',0.7415279,-0.5641171,4.46];
ybs[838]=['',0.7551988,0.8190007,5.88];
ybs[839]=['17 Per',0.753908,0.6133635,4.53];
ybs[840]=['γ1 For',0.7450801,-0.4271909,6.14];
ybs[841]=['γ2 For',0.7452138,-0.4862075,5.39];
ybs[842]=['',0.7606764,0.926433,6.36];
ybs[843]=['σ Ari',0.7532325,0.264686,5.49];
ybs[844]=['η2 For',0.746452,-0.6241243,5.92];
ybs[845]=['',0.7625373,0.8491413,6.26];
ybs[846]=['τ2 Eri',0.7503607,-0.3651322,4.75];
ybs[847]=['η3 For',0.7483177,-0.6212034,5.47];
ybs[848]=['ν Hor',0.7394899,-1.0947101,5.26];
ybs[849]=['',0.7487016,-0.6954775,6.36];
ybs[850]=['τ Per',0.7667447,0.9223183,3.95];
ybs[851]=['20 Per',0.7636251,0.6705583,5.33];
ybs[852]=['',0.7606965,0.2891347,6.31];
ybs[853]=['',0.7571021,-0.2214181,6.04];
ybs[854]=['',0.7539239,-0.5363588,6.4];
ybs[855]=['',0.7585221,-0.1633297,6.32];
ybs[856]=['',0.7748678,1.0751741,5.59];
ybs[857]=['',0.7772339,1.1242391,6.24];
ybs[858]=['',0.761443,-0.3890971,5.95];
ybs[859]=['ψ For',0.7608545,-0.6694061,5.92];
ybs[860]=['',0.7779485,0.8960937,6.22];
ybs[861]=['',0.7764612,0.8245904,6.02];
ybs[862]=['',0.7538144,-1.0965282,6.03];
ybs[863]=['ρ2 Ari',0.7721534,0.3213784,5.91];
ybs[864]=['',0.7616593,-0.869306,4];
ybs[865]=['ρ3 Ari',0.7748807,0.3159885,5.63];
ybs[866]=['',0.7737384,0.1477159,5.97];
ybs[867]=['',0.7625889,-0.8864309,6.21];
ybs[868]=['ν Hyi',0.7433807,-1.3086992,4.75];
ybs[869]=['21 Per',0.779022,0.5587777,5.11];
ybs[870]=['η Eri',0.7741817,-0.1538731,3.89];
ybs[871]=['',0.7751599,-0.0633645,5.17];
ybs[872]=['',0.7825434,0.6753754,6.04];
ybs[873]=['',0.7773379,0.0799823,6.11];
ybs[874]=['47 Ari',0.7821635,0.3621519,5.8];
ybs[875]=['π Per',0.7857459,0.6936582,4.7];
ybs[876]=['',0.7624543,-1.123171,6.56];
ybs[877]=['',0.8245433,1.3874777,5.49];
ybs[878]=['24 Per',0.7868815,0.6154705,4.93];
ybs[879]=['4 Eri',0.7780002,-0.4150478,5.45];
ybs[880]=['',0.7770465,-0.51965,6.29];
ybs[881]=['',0.790776,0.8255642,5.47];
ybs[882]=['',0.7897503,0.7175684,5.89];
ybs[883]=['ε Ari',0.7870938,0.3738678,4.63];
ybs[884]=['ε Ari',0.7870938,0.3738678,4.63];
ybs[885]=['6 Eri',0.78106,-0.410587,5.84];
ybs[886]=['',0.79562,0.9151076,5.28];
ybs[887]=['',0.7957073,0.9151171,6.74];
ybs[888]=['',0.784243,-0.0471505,5.23];
ybs[889]=['',0.7781742,-0.6651395,6.41];
ybs[890]=['',0.7919521,0.6669264,6.11];
ybs[891]=['',0.784462,-0.1692173,6.14];
ybs[892]=['λ Cet',0.7889611,0.1568721,4.7];
ybs[893]=['θ1 Eri',0.7812118,-0.7020333,3.24];
ybs[894]=['θ2 Eri',0.7812554,-0.7020286,4.35];
ybs[895]=['5 Eri',0.7885557,-0.0416152,5.56];
ybs[896]=['',0.7853243,-0.5031101,6.14];
ybs[897]=['ζ For',0.7875806,-0.4397092,5.71];
ybs[898]=['',0.7934627,0.1911226,5.95];
ybs[899]=['',0.7875081,-0.56595,6.31];
ybs[900]=['7 Eri',0.7936212,-0.0488411,6.11];
ybs[901]=['49 Ari',0.7989901,0.4632459,5.9];
ybs[902]=['',0.8510602,1.4232533,5.95];
ybs[903]=['ρ1 Eri',0.7948833,-0.1323425,5.75];
ybs[904]=['',0.7982832,0.0945263,6.25];
ybs[905]=['β Hor',0.781853,-1.116842,4.99];
ybs[906]=['93 Cet',0.8004551,0.0773608,5.61];
ybs[907]=['α Cet',0.8000339,0.0727702,2.53];
ybs[908]=['',0.7981767,-0.1724655,5.83];
ybs[909]=['',0.7992243,-0.1119622,6.19];
ybs[910]=['ε For',0.7963303,-0.4888964,5.89];
ybs[911]=['γ Per',0.8128643,0.9352367,2.93];
ybs[912]=['',0.8060483,0.4947827,6.36];
ybs[913]=['ρ2 Eri',0.8015937,-0.1327447,5.32];
ybs[914]=['',0.8163618,0.9910726,4.76];
ybs[915]=['τ3 Eri',0.7997957,-0.4109334,4.09];
ybs[916]=['',0.8168556,0.9799503,6.11];
ybs[917]=['ρ Per',0.8137343,0.6792629,3.39];
ybs[918]=['',0.8249095,1.1193764,5.89];
ybs[919]=['',0.8145513,0.7096693,6.05];
ybs[920]=['',0.8108058,0.2781172,6.49];
ybs[921]=['ρ3 Eri',0.8084377,-0.1312808,5.26];
ybs[922]=['',0.8102545,0.0339025,6.05];
ybs[923]=['52 Ari',0.8144239,0.4421586,6.8];
ybs[924]=['52 Ari',0.8144239,0.4421586,7];
ybs[925]=['',0.801246,-0.8184802,5.82];
ybs[926]=['',0.8270606,0.912648,6.31];
ybs[927]=['',0.818227,0.2315255,5.62];
ybs[928]=['',0.847409,1.2997394,4.87];
ybs[929]=['',0.8255575,0.8270415,6.41];
ybs[930]=['μ Hor',0.8032983,-1.0412363,5.11];
ybs[931]=['',0.8184449,-0.1049019,5.27];
ybs[932]=['β Per',0.8268941,0.7161625,2.12];
ybs[933]=['ι Per',0.8312708,0.8672631,4.05];
ybs[934]=['53 Ari',0.8228659,0.3134233,6.11];
ybs[935]=['θ Hyi',0.7954867,-1.253541,5.53];
ybs[936]=['54 Ari',0.8269297,0.3293871,6.27];
ybs[937]=['κ Per',0.832884,0.7842507,3.8];
ybs[938]=['',0.8279226,0.1491948,6.28];
ybs[939]=['',0.8234587,-0.4843878,6.19];
ybs[940]=['55 Ari',0.8327401,0.5088327,5.72];
ybs[941]=['',0.8350343,0.4868914,6.42];
ybs[942]=['',0.8363287,0.4707694,6.02];
ybs[943]=['ω Per',0.8404711,0.6926871,4.63];
ybs[944]=['',0.836756,0.2085522,5.98];
ybs[945]=['',0.8459002,0.8342983,6.33];
ybs[946]=['',0.8443975,0.7409301,6.15];
ybs[947]=['δ Ari',0.8412604,0.3456268,4.35];
ybs[948]=['',0.8399157,0.22906,6.12];
ybs[949]=['',0.8355429,-0.4129731,6.38];
ybs[950]=['56 Ari',0.8441416,0.4770508,5.79];
ybs[951]=['',0.8392488,-0.0651921,6.05];
ybs[952]=['',0.8501213,0.8421654,5.9];
ybs[953]=['',0.8387771,-0.2783595,6.26];
ybs[954]=['',0.8444369,0.11758,5.56];
ybs[955]=['',0.8202649,-1.2075529,6.15];
ybs[956]=['',0.8339317,-0.849231,6.12];
ybs[957]=['',0.8857649,1.3579961,5.45];
ybs[958]=['94 Cet',0.8456867,-0.0195516,5.06];
ybs[959]=['α For',0.8418347,-0.504588,3.87];
ybs[960]=['',0.8612591,0.9985986,5.79];
ybs[961]=['',0.9492272,1.4831568,5.61];
ybs[962]=['',0.85657,0.7431421,6.07];
ybs[963]=['',0.8698165,1.1472497,6.36];
ybs[964]=['',0.8427656,-0.7739426,5.93];
ybs[965]=['',0.8625444,0.8903325,5.03];
ybs[966]=['',0.8457371,-0.6260156,6.27];
ybs[967]=['',0.857812,0.5346214,5.52];
ybs[968]=['ζ Ari',0.8555819,0.3686049,4.89];
ybs[969]=['',0.8616656,0.7927358,6.16];
ybs[970]=['',0.8486158,-0.5188613,6.16];
ybs[971]=['',0.8598242,0.5747563,6.31];
ybs[972]=['',0.8609786,0.6067329,6.25];
ybs[973]=['',0.8424451,-0.9991242,5.74];
ybs[974]=['',0.8633002,0.5630088,6.06];
ybs[975]=['',0.8662821,0.7078621,6.45];
ybs[976]=['',0.8547133,-0.4542253,6.25];
ybs[977]=['',0.8154169,-1.3772624,5.57];
ybs[978]=['30 Per',0.8690747,0.7696717,5.47];
ybs[979]=['',0.8597061,-0.1019959,6.17];
ybs[980]=['ζ Eri',0.858836,-0.1526286,4.8];
ybs[981]=['',0.8805943,1.1471218,4.84];
ybs[982]=['',0.8687324,0.6869143,5.96];
ybs[983]=['29 Per',0.8731024,0.8778277,5.15];
ybs[984]=['14 Eri',0.8621497,-0.1584756,6.14];
ybs[985]=['31 Per',0.8752645,0.875604,5.03];
ybs[986]=['',0.8597197,-0.5367386,6.65];
ybs[987]=['',0.8727446,0.5985848,4.82];
ybs[988]=['95 Cet',0.8701277,-0.0149487,5.38];
ybs[989]=['',0.8678932,-0.501311,5.91];
ybs[990]=['15 Eri',0.869501,-0.3916097,4.88];
ybs[991]=['59 Ari',0.8777263,0.4737568,5.9];
ybs[992]=['κ1 Cet',0.8745546,0.0601035,4.83];
ybs[993]=['',0.8710089,-0.3226423,5.71];
ybs[994]=['',0.8644295,-0.8321287,5.85];
ybs[995]=['',0.8795844,0.5082631,4.47];
ybs[996]=['60 Ari',0.8798479,0.4491735,6.12];
ybs[997]=['',0.8872162,0.8577106,5.93];
ybs[998]=['32 Per',0.884992,0.7575076,4.95];
ybs[999]=['τ4 Eri',0.8745323,-0.3784643,3.69];
ybs[1000]=['',0.8747342,-0.4197465,5.61];
ybs[1001]=['τ1 Ari',0.8831967,0.370352,5.28];
ybs[1002]=['ζ1 Ret',0.8646048,-1.0908507,5.54];
ybs[1003]=['κ2 Cet',0.8822065,0.0654199,5.69];
ybs[1004]=['',0.8755102,-0.75043,4.27];
ybs[1005]=['',0.9009291,1.1284832,5.23];
ybs[1006]=['ζ2 Ret',0.8665503,-1.0896513,5.24];
ybs[1007]=['',0.8930952,0.8601886,5.29];
ybs[1008]=['62 Ari',0.8876485,0.4831033,5.52];
ybs[1009]=['',0.8797697,-0.4630967,6.39];
ybs[1010]=['',0.864899,-1.1668025,6.05];
ybs[1011]=['τ2 Ari',0.8898487,0.3632732,5.09];
ybs[1012]=['',0.8826866,-0.4112454,5.52];
ybs[1013]=['α Per',0.8979581,0.871487,1.79];
ybs[1014]=['',0.8864243,-0.4453288,6.35];
ybs[1015]=['',0.8978938,0.5865566,5.61];
ybs[1016]=['',0.904732,0.9423467,6.51];
ybs[1017]=['',0.8823432,-0.8325972,6.39];
ybs[1018]=['64 Ari',0.896764,0.4327655,5.5];
ybs[1019]=['',0.8933157,0.0864583,6.38];
ybs[1020]=['',0.8914227,-0.134779,6.2];
ybs[1021]=['ι Hyi',0.8529299,-1.3493725,5.52];
ybs[1022]=['',0.9011177,0.7213155,6.51];
ybs[1023]=['65 Ari',0.8971918,0.3643381,6.08];
ybs[1024]=['',0.8957856,0.221674,6.04];
ybs[1025]=['',0.9050475,0.8585556,6.09];
ybs[1026]=['ο Tau',0.8984994,0.1588281,3.6];
ybs[1027]=['',0.8925892,-0.5695962,6.5];
ybs[1028]=['',0.9271575,1.2554644,6.32];
ybs[1029]=['',0.9166343,1.0528757,6.49];
ybs[1030]=['',0.9141978,0.857528,4.98];
ybs[1031]=['',0.9195507,1.0473685,4.21];
ybs[1032]=['',0.9085585,0.3285895,6.57];
ybs[1033]=['',0.9178394,0.8712328,5.58];
ybs[1034]=['ξ Tau',0.9088041,0.1710972,3.74];
ybs[1035]=['',0.9095078,0.2234949,6.28];
ybs[1036]=['',0.923139,1.0288331,4.54];
ybs[1037]=['',0.9147243,0.5912717,5.61];
ybs[1038]=['χ1 For',0.9019944,-0.625699,6.39];
ybs[1039]=['',0.9243854,1.0373397,6.13];
ybs[1040]=['34 Per',0.91998,0.865305,4.67];
ybs[1041]=['',0.904259,-0.4755458,5.93];
ybs[1042]=['',0.9232011,0.9690261,5.09];
ybs[1043]=['',0.9201183,0.8204303,6.24];
ybs[1044]=['66 Ari',0.9147535,0.399227,6.03];
ybs[1045]=['',0.902893,-0.7254656,6.32];
ybs[1046]=['',0.9119155,-0.1957666,5.73];
ybs[1047]=['',0.9253446,0.8407697,5.82];
ybs[1048]=['σ Per',0.9251558,0.8388792,4.36];
ybs[1049]=['',0.8907228,-1.2139276,6.15];
ybs[1050]=['χ2 For',0.9090957,-0.621531,5.71];
ybs[1051]=['',0.9490056,1.2813141,6.57];
ybs[1052]=['',0.9292309,0.8600689,6.29];
ybs[1053]=['χ3 For',0.9118568,-0.6245364,6.5];
ybs[1054]=['',0.930678,0.8634021,6.39];
ybs[1055]=['',0.9191747,-0.1175582,5.99];
ybs[1056]=['4 Tau',0.9229828,0.1990632,5.14];
ybs[1057]=['',0.9187824,-0.2200036,5.59];
ybs[1058]=['',0.9320098,0.8393628,5.47];
ybs[1059]=['',0.8975775,-1.208906,5.96];
ybs[1060]=['',0.92759,0.4824202,5.96];
ybs[1061]=['5 Tau',0.9250587,0.2269901,4.11];
ybs[1062]=['',0.9243569,0.1092152,5.94];
ybs[1063]=['',0.9389536,1.026825,6.4];
ybs[1064]=['36 Per',0.9331804,0.8050361,5.31];
ybs[1065]=['17 Eri',0.9234403,-0.0873755,4.73];
ybs[1066]=['',0.9395254,1.0111839,6.37];
ybs[1067]=['',0.9340441,0.7840664,6.41];
ybs[1068]=['',0.9390992,0.9606717,5.98];
ybs[1069]=['',0.9336573,0.6201124,5.9];
ybs[1070]=['',0.9190642,-0.7428959,5.78];
ybs[1071]=['',0.9204785,-0.7208342,6.12];
ybs[1072]=['',0.9455555,1.0490869,6.46];
ybs[1073]=['',0.9378664,0.6975597,5.81];
ybs[1074]=['6 Tau',0.9324915,0.1647911,5.77];
ybs[1075]=['',0.9684044,1.3230452,6.27];
ybs[1076]=['',0.9218526,-0.8256485,5.99];
ybs[1077]=['',0.928403,-0.4458551,6.38];
ybs[1078]=['κ Ret',0.9150858,-1.0972509,4.72];
ybs[1079]=['ε Eri',0.9334016,-0.1638902,3.73];
ybs[1080]=['',0.9394565,0.3124205,6.17];
ybs[1081]=['7 Tau',0.9410022,0.4281626,5.92];
ybs[1082]=['ψ Per',0.9510112,0.8422848,4.23];
ybs[1083]=['τ5 Eri',0.9367776,-0.37638,4.27];
ybs[1084]=['',0.9420954,0.1131867,6.49];
ybs[1085]=['',0.9302149,-0.87808,5.68];
ybs[1086]=['',0.9407769,-0.1710628,6.25];
ybs[1087]=['',0.9209963,-1.1592587,5.83];
ybs[1088]=['',0.9371378,-0.5412709,6.2];
ybs[1089]=['',0.9597242,0.994813,6.3];
ybs[1090]=['',0.9397742,-0.5551408,6.4];
ybs[1091]=['',0.9304836,-1.0637552,6.41];
ybs[1092]=['',0.9572987,0.7443663,6.42];
ybs[1093]=['',0.9465766,-0.1941977,5.57];
ybs[1094]=['',0.9505195,0.0114203,5.71];
ybs[1095]=['20 Eri',0.9478204,-0.3036902,5.23];
ybs[1096]=['10 Tau',0.9508851,0.0081714,4.28];
ybs[1097]=['',0.9553508,0.2704729,6.39];
ybs[1098]=['',0.9607833,0.3661956,6.5];
ybs[1099]=['',0.9365744,-1.146625,6.75];
ybs[1100]=['',0.9771854,1.1044597,5.1];
ybs[1101]=['',0.9504699,-0.701766,4.58];
ybs[1102]=['',1.1259969,1.5101774,5.86];
ybs[1103]=['',0.9577023,-0.1278592,5.85];
ybs[1104]=['',0.913162,-1.3662842,5.7];
ybs[1105]=['',0.9625066,0.2897618,6.16];
ybs[1106]=['21 Eri',0.9600727,-0.0970483,5.96];
ybs[1107]=['',0.9792042,1.0477813,5.76];
ybs[1108]=['',0.9707174,0.6570244,5.57];
ybs[1109]=['τ For',0.9583955,-0.4865502,6.01];
ybs[1110]=['12 Tau',0.9639565,0.0544935,5.57];
ybs[1111]=['',0.9628416,-0.0580786,6.23];
ybs[1112]=['',0.9616998,-0.1810208,6.19];
ybs[1113]=['11 Tau',0.9686593,0.4432148,6.11];
ybs[1114]=['',0.9644481,-0.0184185,6.12];
ybs[1115]=['',0.9649005,-0.2646177,6.33];
ybs[1116]=['22 Eri',0.9671524,-0.0898071,5.53];
ybs[1117]=['δ Per',0.9790989,0.8351655,3.01];
ybs[1118]=['40 Per',0.975989,0.593922,4.97];
ybs[1119]=['',0.9947301,1.1739825,5.8];
ybs[1120]=['',0.9695407,-0.2048719,6.49];
ybs[1121]=['13 Tau',0.9752071,0.3449564,5.69];
ybs[1122]=['',0.9843091,0.8480046,6.06];
ybs[1123]=['',0.9699251,-0.3406885,6.59];
ybs[1124]=['',0.9941854,1.1066714,4.8];
ybs[1125]=['',0.9866681,0.8056956,6.11];
ybs[1126]=['ο Per',0.9844056,0.5646446,3.83];
ybs[1127]=['14 Tau',0.9816298,0.3443301,6.14];
ybs[1128]=['',0.9854782,0.6374523,5.59];
ybs[1129]=['δ For',0.9733017,-0.5563056,5];
ybs[1130]=['ν Per',0.9887032,0.7442372,3.77];
ybs[1131]=['δ Eri',0.978406,-0.1692866,3.54];
ybs[1132]=['',0.9846469,0.3663794,6.1];
ybs[1133]=['',1.0096452,1.2380022,5.44];
ybs[1134]=['',0.9797588,-0.181894,5.6];
ybs[1135]=['16 Tau',0.9862235,0.4250345,5.46];
ybs[1136]=['',0.9923509,0.7983946,5.66];
ybs[1137]=['17 Tau',0.9865302,0.4219602,3.7];
ybs[1138]=['',0.9756313,-0.6501259,4.59];
ybs[1139]=['18 Tau',0.987807,0.4346263,5.64];
ybs[1140]=['19 Tau',0.9879977,0.4281343,4.3];
ybs[1141]=['24 Eri',0.9841541,-0.0191927,5.25];
ybs[1142]=['',0.9998753,0.9771142,6.1];
ybs[1143]=['γ Cam',1.0147135,1.2460417,4.63];
ybs[1144]=['20 Tau',0.9906935,0.4263942,3.87];
ybs[1145]=['25 Eri',0.98607,-0.0040745,5.55];
ybs[1146]=['21 Tau',0.9910494,0.4296564,5.76];
ybs[1147]=['22 Tau',0.9916668,0.4291899,6.43];
ybs[1148]=['29 Tau',0.9894472,0.1066902,5.35];
ybs[1149]=['',0.9414661,-1.3658253,6.29];
ybs[1150]=['',1.0097552,1.1447132,4.47];
ybs[1151]=['23 Tau',0.9928617,0.4190699,4.18];
ybs[1152]=['',0.981001,-0.7085452,6.45];
ybs[1153]=['',1.0098035,1.1058113,5.85];
ybs[1154]=['',0.9915858,0.11983,5.91];
ybs[1155]=['',1.0028051,0.886599,6.14];
ybs[1156]=['',1.0078135,0.9979673,6.46];
ybs[1157]=['π Eri',0.9909578,-0.210119,4.42];
ybs[1158]=['',0.9999968,0.5875116,6.57];
ybs[1159]=['',0.9996694,0.5629903,6.25];
ybs[1160]=['η Tau',0.9979238,0.4217958,2.87];
ybs[1161]=['',1.0198569,1.1967317,6.32];
ybs[1162]=['',0.9838063,-0.837724,6.49];
ybs[1163]=['',0.982128,-0.9461502,6.3];
ybs[1164]=['',0.9856815,-0.8254807,5.73];
ybs[1165]=['',1.006003,0.7683715,6.02];
ybs[1166]=['σ For',0.9917533,-0.5109526,5.9];
ybs[1167]=['',1.0016687,0.4098534,5.45];
ybs[1168]=['τ6 Eri',0.9936818,-0.4046942,4.23];
ybs[1169]=['30 Tau',1.0009382,0.1955667,5.07];
ybs[1170]=['β Ret',0.9793327,-1.1299831,3.85];
ybs[1171]=['',1.0101596,0.7859002,5.66];
ybs[1172]=['42 Per',1.0072528,0.5786225,5.11];
ybs[1173]=['27 Tau',1.0052413,0.4208818,3.63];
ybs[1174]=['',0.9955617,-0.5208007,6.55];
ybs[1175]=['28 Tau',1.0053534,0.422336,5.09];
ybs[1176]=['τ7 Eri',0.9972,-0.4156083,5.24];
ybs[1177]=['',1.0022586,0.0050518,5.91];
ybs[1178]=['',1.0076885,0.4149145,6.17];
ybs[1179]=['ρ For',0.9981605,-0.5254447,5.54];
ybs[1180]=['',1.0084743,0.3893052,6.07];
ybs[1181]=['',0.997448,-0.6290824,6.21];
ybs[1182]=['',1.0013867,-0.3637499,5.81];
ybs[1183]=['',1.0103306,0.4475089,5.26];
ybs[1184]=['',1.0007068,-0.6555538,5.4];
ybs[1185]=['',1.0007433,-0.6555248,4.73];
ybs[1186]=['',1.0175752,0.6007321,5.77];
ybs[1187]=['',1.0271102,1.0128912,5.8];
ybs[1188]=['',1.0158588,0.3855791,6.83];
ybs[1189]=['',1.0140533,0.2287495,6.3];
ybs[1190]=['',1.0045431,-0.6307426,4.17];
ybs[1191]=['',1.0417557,1.2545384,6.34];
ybs[1192]=['',1.0182131,0.5450403,6.25];
ybs[1193]=['',1.025987,0.8501503,5.76];
ybs[1194]=['31 Tau',1.0170846,0.1151039,5.67];
ybs[1195]=['',1.0096413,-0.6346781,6.86];
ybs[1196]=['',1.0224964,0.3034548,5.97];
ybs[1197]=['30 Eri',1.0197425,-0.0925272,5.48];
ybs[1198]=['ζ Per',1.0272364,0.557509,2.85];
ybs[1199]=['',1.0439326,1.1018263,5.03];
ybs[1200]=['',1.0424097,1.0675621,5];
ybs[1201]=['',1.0216195,-0.3206986,6.22];
ybs[1202]=['',1.036091,0.836534,5.37];
ybs[1203]=['γ Hyi',0.9902073,-1.2946215,3.24];
ybs[1204]=['',1.032667,0.5428777,6.1];
ybs[1205]=['43 Per',1.0390784,0.8858152,5.28];
ybs[1206]=['32 Eri',1.0267805,-0.0505009,6.14];
ybs[1207]=['32 Eri',1.0267877,-0.0505348,4.79];
ybs[1208]=['τ8 Eri',1.0235594,-0.4285295,4.65];
ybs[1209]=['',1.0228914,-0.6051511,5.11];
ybs[1210]=['',1.0376199,0.6132983,5.49];
ybs[1211]=['',1.0218483,-0.8174062,5.93];
ybs[1212]=['',1.0307859,-0.2101426,6];
ybs[1213]=['32 Tau',1.0388277,0.393331,5.63];
ybs[1214]=['',1.0258446,-0.7033312,5.71];
ybs[1215]=['ε Per',1.0438596,0.699318,2.89];
ybs[1216]=['33 Tau',1.0397035,0.4055032,6.06];
ybs[1217]=['',1.0413938,0.4279521,6.16];
ybs[1218]=['',1.0444777,0.6086323,6.53];
ybs[1219]=['',1.0389979,0.1064321,6.09];
ybs[1220]=['',1.0367915,-0.1691666,6.19];
ybs[1221]=['',1.046559,0.6788929,6.3];
ybs[1222]=['',1.0258434,-0.9185891,6.46];
ybs[1223]=['ξ Per',1.0485098,0.6256714,4.04];
ybs[1224]=['',1.0517182,0.6785398,6.38];
ybs[1225]=['',1.1064819,1.4093623,5.1];
ybs[1226]=['γ Eri',1.0427743,-0.2347625,2.95];
ybs[1227]=['',1.0466947,-0.0944689,5.83];
ybs[1228]=['',1.050694,0.1813013,6.37];
ybs[1229]=['',1.0585214,0.6465788,6.41];
ybs[1230]=['',1.0492238,-0.2184645,5.6];
ybs[1231]=['',1.031203,-1.1066246,6.14];
ybs[1232]=['',1.0550052,0.3028707,6.32];
ybs[1233]=['',1.0559007,0.3185287,5.89];
ybs[1234]=['λ Tau',1.0551333,0.218983,3.47];
ybs[1235]=['τ9 Eri',1.0506864,-0.4181719,4.66];
ybs[1236]=['',1.0826905,1.1996342,5.87];
ybs[1237]=['',1.0741431,1.0334149,5.06];
ybs[1238]=['',1.0598046,0.1754725,5.67];
ybs[1239]=['35 Eri',1.0584264,-0.0260673,5.28];
ybs[1240]=['',1.0435184,-0.9956226,6.05];
ybs[1241]=['',1.0537244,-0.5311776,5.93];
ybs[1242]=['δ Ret',1.0431545,-1.0706328,4.56];
ybs[1243]=['',1.0846625,1.1444924,6.17];
ybs[1244]=['',1.0631702,-0.0037207,5.38];
ybs[1245]=['',1.0507493,-0.8989823,6.51];
ybs[1246]=['ν Tau',1.065739,0.1054987,3.91];
ybs[1247]=['36 Tau',1.0715957,0.4216845,5.47];
ybs[1248]=['40 Tau',1.0682821,0.0958318,5.33];
ybs[1249]=['',1.0692402,0.1440304,5.46];
ybs[1250]=['',1.0830441,0.9435678,6.31];
ybs[1251]=['37 Tau',1.0729778,0.3863584,4.36];
ybs[1252]=['',1.0700423,0.0502999,5.36];
ybs[1253]=['',1.0660303,-0.3506153,6.46];
ybs[1254]=['',1.0669094,-0.3508641,7.01];
ybs[1255]=['',1.0894053,1.088793,6.99];
ybs[1256]=['λ Per',1.0826284,0.8797376,4.29];
ybs[1257]=['39 Tau',1.0757761,0.3850785,5.9];
ybs[1258]=['',1.0693491,-0.2885696,6.39];
ybs[1259]=['γ Ret',1.0524235,-1.0838989,4.51];
ybs[1260]=['',1.0704905,-0.222312,5.61];
ybs[1261]=['ι Ret',1.0543388,-1.0650428,4.97];
ybs[1262]=['',1.0715687,-0.3547701,6.13];
ybs[1263]=['41 Tau',1.0815279,0.4826519,5.2];
ybs[1264]=['ψ Tau',1.0833369,0.5071076,5.23];
ybs[1265]=['',1.0961136,1.04651,6.28];
ybs[1266]=['',0.955975,-1.481889,6.41];
ybs[1267]=['',1.0774353,-0.1536212,6.26];
ybs[1268]=['48 Per',1.0915174,0.8336645,4.04];
ybs[1269]=['',1.0763404,-0.3570571,6.34];
ybs[1270]=['',1.0754022,-0.4816625,5.59];
ybs[1271]=['',1.0951917,0.9578629,6.18];
ybs[1272]=['49 Per',1.0891662,0.6594019,6.09];
ybs[1273]=['50 Per',1.0907317,0.6648436,5.51];
ybs[1274]=['',1.0858574,0.2655734,6.01];
ybs[1275]=['',1.0871957,0.3035659,5.89];
ybs[1276]=['',1.1173176,1.2597247,6.03];
ybs[1277]=['',1.1124252,1.1964691,6.32];
ybs[1278]=['ω1 Tau',1.0924092,0.3431661,5.5];
ybs[1279]=['',1.0915896,0.2347679,5.95];
ybs[1280]=['',1.0825177,-0.748105,6.59];
ybs[1281]=['',1.1008847,0.587105,5.72];
ybs[1282]=['44 Tau',1.0999281,0.4630862,5.41];
ybs[1283]=['',1.091828,-0.2850651,5.37];
ybs[1284]=['',1.1919438,1.4634784,5.57];
ybs[1285]=['37 Eri',1.0968373,-0.1199317,5.44];
ybs[1286]=['',1.0873392,-0.799562,6.59];
ybs[1287]=['45 Tau',1.1014272,0.0973007,5.72];
ybs[1288]=['',1.0986099,-0.1530234,5.7];
ybs[1289]=['',1.080263,-1.1199549,6.38];
ybs[1290]=['',1.1069796,0.302445,6.09];
ybs[1291]=['',1.1201707,1.0037452,6.08];
ybs[1292]=['',1.1086006,0.3920843,6.12];
ybs[1293]=['ο1 Eri',1.1033407,-0.1184354,4.04];
ybs[1294]=['',1.0974941,-0.6147346,6.44];
ybs[1295]=['',1.1017534,-0.3543773,5.79];
ybs[1296]=['',1.1142339,0.6635317,6.45];
ybs[1297]=['δ Hor',1.0975016,-0.7320162,4.93];
ybs[1298]=['μ Per',1.1188057,0.8457801,4.14];
ybs[1299]=['',1.1985684,1.4553105,5.46];
ybs[1300]=['',1.1288553,1.0803501,5.7];
ybs[1301]=['52 Per',1.1182729,0.7074487,4.71];
ybs[1302]=['',1.111358,0.1791246,6.23];
ybs[1303]=['',1.1110612,0.1560577,6.51];
ybs[1304]=['46 Tau',1.1111547,0.1355595,5.29];
ybs[1305]=['',1.1125387,0.2234734,6.25];
ybs[1306]=['47 Tau',1.1128986,0.1625654,4.84];
ybs[1307]=['',1.1112481,-0.0191788,6.44];
ybs[1308]=['',1.1294192,1.0107148,5.71];
ybs[1309]=['',1.1271689,0.9365662,5.19];
ybs[1310]=['',1.1158246,0.1756066,5.22];
ybs[1311]=['',1.104702,-0.7734758,6.71];
ybs[1312]=['',1.1811929,1.4114191,5.43];
ybs[1313]=['39 Eri',1.1142697,-0.1781259,4.87];
ybs[1314]=['48 Tau',1.1210932,0.2696608,6.32];
ybs[1315]=['μ Tau',1.1198472,0.1560711,4.29];
ybs[1316]=['',1.1193033,0.109079,6.93];
ybs[1317]=['',1.1195502,0.1088506,6.31];
ybs[1318]=['',1.109588,-0.7034869,6.37];
ybs[1319]=['',1.1335552,0.8786726,4.61];
ybs[1320]=['ο2 Eri',1.1181776,-0.1326912,4.43];
ybs[1321]=['α Hor',1.1112542,-0.7372912,3.86];
ybs[1322]=['',1.1458928,1.137746,5.27];
ybs[1323]=['',1.1325493,0.736352,6.22];
ybs[1324]=['ω2 Tau',1.1277686,0.3600231,4.94];
ybs[1325]=['',1.1378094,0.8743598,5.45];
ybs[1326]=['51 Tau',1.1327224,0.3774772,5.65];
ybs[1327]=['',1.1271516,-0.1120978,5.94];
ybs[1328]=['',1.1421141,0.8895704,5.55];
ybs[1329]=['',1.1324054,0.1664235,6.54];
ybs[1330]=['',1.1501156,1.060856,5.39];
ybs[1331]=['α Ret',1.1113019,-1.0894906,3.35];
ybs[1332]=['',1.1417147,0.7305225,5.92];
ybs[1333]=['γ Dor',1.1194625,-0.897741,4.25];
ybs[1334]=['53 Tau',1.1372822,0.3698429,5.35];
ybs[1335]=['',1.1130268,-1.0845728,5.45];
ybs[1336]=['56 Tau',1.1380763,0.3808613,5.38];
ybs[1337]=['',1.1499242,0.9870479,5.88];
ybs[1338]=['54 Per',1.1420851,0.6041354,4.93];
ybs[1339]=['',1.1409046,0.5585263,6.16];
ybs[1340]=['',1.130798,-0.3607027,6];
ybs[1341]=['γ Tau',1.1386529,0.2735902,3.65];
ybs[1342]=['υ4 Eri',1.1286605,-0.5890369,3.56];
ybs[1343]=['φ Tau',1.1415341,0.4781961,4.95];
ybs[1344]=['',1.1377354,0.177492,6.31];
ybs[1345]=['53 Per',1.1477388,0.8123824,4.85];
ybs[1346]=['57 Tau',1.1393338,0.2457994,5.59];
ybs[1347]=['',1.1550914,1.0413137,6.19];
ybs[1348]=['',1.1322788,-0.4000578,6.07];
ybs[1349]=['',1.1414884,0.3279521,6.12];
ybs[1350]=['ε Ret',1.1206801,-1.0341462,4.44];
ybs[1351]=['58 Tau',1.1421772,0.2642948,5.26];
ybs[1352]=['',1.1199064,-1.0628847,6.37];
ybs[1353]=['',1.1433358,0.2428057,6.17];
ybs[1354]=['',1.1336925,-0.5909076,6.37];
ybs[1355]=['',1.1422398,0.1078353,5.77];
ybs[1356]=['',1.1429066,0.1618473,6.53];
ybs[1357]=['',1.1416656,-0.1081727,6.27];
ybs[1358]=['',1.1419205,-0.1316818,5.85];
ybs[1359]=['',1.1341355,-0.7717783,5.34];
ybs[1360]=['',1.1308544,-0.9217307,6.09];
ybs[1361]=['',1.1453802,-0.0008852,5.86];
ybs[1362]=['',1.141204,-0.3593979,5.38];
ybs[1363]=['60 Tau',1.1484855,0.2465149,5.72];
ybs[1364]=['χ Tau',1.1511857,0.4481298,5.37];
ybs[1365]=['',1.1501246,0.36422,5.91];
ybs[1366]=['',1.1564148,0.741317,6.23];
ybs[1367]=['θ Ret',1.1253195,-1.1031585,5.87];
ybs[1368]=['δ1 Tau',1.1524304,0.3069884,3.76];
ybs[1369]=['',1.1448155,-0.4482176,6.01];
ybs[1370]=['',1.1551936,0.3670179,5.99];
ybs[1371]=['63 Tau',1.1545059,0.293628,5.64];
ybs[1372]=['55 Per',1.1598625,0.5964916,5.73];
ybs[1373]=['62 Tau',1.1573108,0.4249397,6.36];
ybs[1374]=['56 Per',1.1604512,0.5935089,5.76];
ybs[1375]=['δ2 Tau',1.1574969,0.305258,4.8];
ybs[1376]=['66 Tau',1.156217,0.1659295,5.12];
ybs[1377]=['',1.1725466,1.005832,6.32];
ybs[1378]=['ξ Eri',1.154987,-0.0645636,5.17];
ybs[1379]=['',1.1517102,-0.4336372,5.83];
ybs[1380]=['',1.1612849,0.3331377,5.98];
ybs[1381]=['',1.1514006,-0.6195632,6.39];
ybs[1382]=['κ1 Tau',1.163222,0.3898962,4.22];
ybs[1383]=['κ2 Tau',1.1634295,0.3882523,5.28];
ybs[1384]=['δ3 Tau',1.163595,0.3136972,4.29];
ybs[1385]=['',1.1668059,0.5495003,5.28];
ybs[1386]=['70 Tau',1.1641004,0.2790126,6.46];
ybs[1387]=['υ Tau',1.1673443,0.3989595,4.28];
ybs[1388]=['43 Eri',1.1554336,-0.5929009,3.96];
ybs[1389]=['71 Tau',1.1672534,0.2733782,4.49];
ybs[1390]=['η Ret',1.1436741,-1.1054752,5.24];
ybs[1391]=['π Tau',1.168357,0.2575858,4.69];
ybs[1392]=['',1.16704,0.1507156,6.06];
ybs[1393]=['',1.1593383,-0.605838,6.55];
ybs[1394]=['72 Tau',1.1716432,0.4021418,5.53];
ybs[1395]=['',1.1697099,0.0370752,6.23];
ybs[1396]=['',1.2039072,1.2665874,5.94];
ybs[1397]=['',1.1720504,0.1964729,5.88];
ybs[1398]=['',1.1747416,0.3781135,5.72];
ybs[1399]=['',1.1604827,-0.7699548,6.39];
ybs[1400]=['',1.1545794,-0.9952766,6.29];
ybs[1401]=['',1.1788074,0.5306724,6.4];
ybs[1402]=['75 Tau',1.176415,0.2863011,4.97];
ybs[1403]=['76 Tau',1.1761408,0.2580466,5.9];
ybs[1404]=['ε Tau',1.1772858,0.3355276,3.53];
ybs[1405]=['',1.1685559,-0.4195162,6.11];
ybs[1406]=['θ1 Tau',1.1769903,0.2793623,3.84];
ybs[1407]=['θ2 Tau',1.1773654,0.2777666,3.4];
ybs[1408]=['',1.1742771,0.0332125,6.15];
ybs[1409]=['79 Tau',1.1780324,0.2284888,5.03];
ybs[1410]=['',1.1763344,0.0248699,5.55];
ybs[1411]=['',1.1579305,-1.0680101,5.94];
ybs[1412]=['1 Cam',1.1940787,0.9416609,5.77];
ybs[1413]=['',1.1680691,-0.8186054,6.1];
ybs[1414]=['',1.1866412,0.5672517,6.21];
ybs[1415]=['',1.1817808,0.1843977,6.79];
ybs[1416]=['',1.1761429,-0.3388425,5.96];
ybs[1417]=['80 Tau',1.1838237,0.273692,5.58];
ybs[1418]=['',1.1783947,-0.2269708,5.6];
ybs[1419]=['',1.1903421,0.6990565,6.26];
ybs[1420]=['',1.1831905,0.1798718,6.48];
ybs[1421]=['δ Men',1.1198741,-1.3991336,5.69];
ybs[1422]=['',1.1856687,0.2833897,4.78];
ybs[1423]=['81 Tau',1.1860294,0.2746285,5.48];
ybs[1424]=['',1.1729898,-0.7315658,6.44];
ybs[1425]=['83 Tau',1.1858454,0.2402893,5.4];
ybs[1426]=['',1.1829423,-0.2364718,6.24];
ybs[1427]=['85 Tau',1.1913376,0.2774064,6.02];
ybs[1428]=['',1.1778353,-0.8110795,6.16];
ybs[1429]=['57 Per',1.1993162,0.7523355,6.09];
ybs[1430]=['',1.1694083,-1.0904152,5.75];
ybs[1431]=['',1.1919332,0.0951634,6.39];
ybs[1432]=['45 Eri',1.1908767,-0.0000231,4.91];
ybs[1433]=['',1.188481,-0.2374032,6.21];
ybs[1434]=['',1.1842908,-0.621519,5.96];
ybs[1435]=['',1.2144912,1.1222804,5.94];
ybs[1436]=['',1.1940251,-0.0552735,5.81];
ybs[1437]=['',1.1987679,0.3151791,6.25];
ybs[1438]=['δ Cae',1.1844853,-0.7838404,5.07];
ybs[1439]=['ρ Tau',1.1999655,0.259811,4.65];
ybs[1440]=['',1.2039272,0.5061864,5.88];
ybs[1441]=['',1.1995777,0.1650159,6.01];
ybs[1442]=['',1.1970193,-0.1875171,6.06];
ybs[1443]=['',1.2009201,0.0979151,5.68];
ybs[1444]=['46 Eri',1.1995384,-0.1168891,5.72];
ybs[1445]=['',1.2009458,-0.1186176,6.09];
ybs[1446]=['47 Eri',1.2007108,-0.1429403,5.11];
ybs[1447]=['',1.2006936,-0.1558364,5.26];
ybs[1448]=['υ1 Eri',1.196938,-0.5187956,4.51];
ybs[1449]=['58 Per',1.213509,0.7209077,4.25];
ybs[1450]=['',1.2082771,0.3477118,6.36];
ybs[1451]=['ν Men',1.1310959,-1.4229985,5.79];
ybs[1452]=['α Tau',1.2090649,0.288849,0.85];
ybs[1453]=['88 Tau',1.2076903,0.1780521,4.25];
ybs[1454]=['',1.2117898,0.4080793,6.02];
ybs[1455]=['',1.2051982,-0.1692208,6.37];
ybs[1456]=['',1.2038621,-0.346961,6.13];
ybs[1457]=['',1.2088611,-0.0623308,6.33];
ybs[1458]=['ν Eri',1.2101422,-0.0578051,3.93];
ybs[1459]=['υ2 Eri',1.205801,-0.5326971,3.82];
ybs[1460]=['α Dor',1.1974702,-0.9599882,3.27];
ybs[1461]=['2 Cam',1.2287088,0.9339558,5.35];
ybs[1462]=['3 Cam',1.2284252,0.9270913,5.05];
ybs[1463]=['',1.2604802,1.3377362,6.49];
ybs[1464]=['',1.2142544,0.0181238,5.31];
ybs[1465]=['',1.2206916,0.4708801,6.47];
ybs[1466]=['',1.2194434,0.361707,5.92];
ybs[1467]=['89 Tau',1.2188094,0.280526,5.79];
ybs[1468]=['90 Tau',1.2186943,0.2190468,4.27];
ybs[1469]=['51 Eri',1.2157703,-0.0424712,5.23];
ybs[1470]=['',1.1946396,-1.0957463,5.79];
ybs[1471]=['',1.2114568,-0.5354033,6.3];
ybs[1472]=['',1.2245138,0.4408241,6.22];
ybs[1473]=['σ1 Tau',1.2231508,0.2764406,5.07];
ybs[1474]=['σ2 Tau',1.2236859,0.2785049,4.69];
ybs[1475]=['',1.2226617,0.1380561,5.39];
ybs[1476]=['53 Eri',1.2178981,-0.2489577,3.87];
ybs[1477]=['',1.2345403,0.8436721,5.67];
ybs[1478]=['',1.2210851,-0.2109008,5.01];
ybs[1479]=['93 Tau',1.2269678,0.2135675,5.46];
ybs[1480]=['',1.1369832,-1.4460351,6.76];
ybs[1481]=['',1.2439245,1.0394814,6.5];
ybs[1482]=['',1.2229059,-0.2499318,5.45];
ybs[1483]=['',1.2253503,-0.0176956,6.1];
ybs[1484]=['',1.2358202,0.6687774,5.99];
ybs[1485]=['',1.2331433,0.5000911,5.78];
ybs[1486]=['',1.2724777,1.3260194,6.06];
ybs[1487]=['',1.2086667,-1.0827501,5.4];
ybs[1488]=['',1.2432201,0.8728562,5.87];
ybs[1489]=['59 Per',1.2407694,0.7575135,5.29];
ybs[1490]=['',1.2259571,-0.4266233,5.58];
ybs[1491]=['54 Eri',1.2275727,-0.3426614,4.32];
ybs[1492]=['τ Tau',1.2369077,0.4013322,4.28];
ybs[1493]=['',1.2199007,-0.9011733,6.44];
ybs[1494]=['95 Tau',1.241251,0.4210801,6.13];
ybs[1495]=['',1.2463346,0.7125074,6.08];
ybs[1496]=['',1.2441359,0.5742519,6.45];
ybs[1497]=['α Cae',1.2270871,-0.7299887,4.45];
ybs[1498]=['β Cae',1.2338737,-0.6476314,5.05];
ybs[1499]=['',1.2245375,-1.0280827,6.53];
ybs[1500]=['55 Eri',1.2416303,-0.1528296,6.82];
ybs[1501]=['55 Eri',1.2416666,-0.1528732,6.7];
ybs[1502]=['',1.2460167,0.1951766,5.4];
ybs[1503]=['56 Eri',1.2438724,-0.1477723,5.9];
ybs[1504]=['',1.2389616,-0.536308,5.68];
ybs[1505]=['',1.2781127,1.23875,6.37];
ybs[1506]=['4 Cam',1.2641391,0.9912088,5.3];
ybs[1507]=['',1.2520497,0.4130167,6.35];
ybs[1508]=['',1.2437138,-0.3251511,5.53];
ybs[1509]=['',1.2573308,0.7042108,5.97];
ybs[1510]=['',1.2645031,0.9710543,6.26];
ybs[1511]=['λ Pic',1.2361374,-0.88041,5.31];
ybs[1512]=['',1.2543587,0.3276122,6.01];
ybs[1513]=['',1.2409792,-0.7160667,6.25];
ybs[1514]=['',1.2530108,0.2049276,5.37];
ybs[1515]=['μ Eri',1.2502136,-0.0561737,4.02];
ybs[1516]=['',1.2477023,-0.370833,5.72];
ybs[1517]=['',1.2541504,-0.0509403,6.33];
ybs[1518]=['',1.3275599,1.4175969,5.07];
ybs[1519]=['',1.2504646,-0.5928688,6.86];
ybs[1520]=['',1.2533625,-0.4895943,6.19];
ybs[1521]=['',1.2506115,-0.6862736,6.05];
ybs[1522]=['',1.2828825,1.1089498,5.44];
ybs[1523]=['',1.2681817,0.5693728,5.86];
ybs[1524]=['',1.2676775,0.549283,5.58];
ybs[1525]=['κ Dor',1.2420753,-1.0418893,5.27];
ybs[1526]=['',1.2106559,-1.3546565,6.05];
ybs[1527]=['58 Eri',1.2589217,-0.2949468,5.51];
ybs[1528]=['',1.2710146,0.6548888,4.88];
ybs[1529]=['',1.2645845,0.0632331,6.03];
ybs[1530]=['',1.2771602,0.851271,5.66];
ybs[1531]=['',1.2636734,-0.0984216,5.78];
ybs[1532]=['96 Tau',1.269334,0.2781764,6.08];
ybs[1533]=['59 Eri',1.2630298,-0.2843954,5.77];
ybs[1534]=['ζ Cae',1.259371,-0.523339,6.37];
ybs[1535]=['',1.2442426,-1.102927,6.46];
ybs[1536]=['μ Men',1.2342543,-1.2373244,5.54];
ybs[1537]=['α Cam',1.2919887,1.1584567,4.29];
ybs[1538]=['π3 Ori',1.2694827,0.1220949,3.19];
ybs[1539]=['π2 Ori',1.2729156,0.1559285,4.36];
ybs[1540]=['',1.2681887,-0.2397295,6.26];
ybs[1541]=['',1.2862722,0.9228124,6.41];
ybs[1542]=['97 Tau',1.2765983,0.3293978,5.1];
ybs[1543]=['',1.2618502,-0.7669875,6.72];
ybs[1544]=['60 Eri',1.2702393,-0.2824503,5.03];
ybs[1545]=['',1.2839063,0.7438473,5.71];
ybs[1546]=['2 Aur',1.2828616,0.6411606,4.78];
ybs[1547]=['π4 Ori',1.2754007,0.0984101,3.69];
ybs[1548]=['',1.2778034,0.1746766,6.11];
ybs[1549]=['',1.2831111,0.4874737,5.97];
ybs[1550]=['5 Cam',1.2947526,0.9650043,5.52];
ybs[1551]=['ο1 Ori',1.281487,0.2492922,4.74];
ybs[1552]=['',1.2694557,-0.7205906,6.07];
ybs[1553]=['',1.2929851,0.7695591,6.08];
ybs[1554]=['',1.2750342,-0.6086477,5.86];
ybs[1555]=['ω Eri',1.2823983,-0.0945982,4.39];
ybs[1556]=['',1.2991778,0.9232871,5.75];
ybs[1557]=['5 Ori',1.2847803,0.0443402,5.33];
ybs[1558]=['ι Pic',1.2713963,-0.9324879,5.61];
ybs[1559]=['ι Pic',1.2714763,-0.9324589,6.42];
ybs[1560]=['',1.2871489,0.0279539,6.61];
ybs[1561]=['',1.2923189,0.3406349,6.37];
ybs[1562]=['π5 Ori',1.2885816,0.0431549,3.72];
ybs[1563]=['7 Cam',1.3043621,0.9386847,4.47];
ybs[1564]=['6 Ori',1.2911994,0.1999778,5.19];
ybs[1565]=['π1 Ori',1.291657,0.1777191,4.65];
ybs[1566]=['',1.2911394,0.1363266,5.33];
ybs[1567]=['',1.3306353,1.2967264,6.06];
ybs[1568]=['',1.2989776,0.6318067,6.07];
ybs[1569]=['',1.2911046,0.0087137,5.99];
ybs[1570]=['',1.2981405,0.4297568,6.37];
ybs[1571]=['',1.295938,0.263043,5.81];
ybs[1572]=['ι Aur',1.301717,0.5793929,2.69];
ybs[1573]=['',1.2961934,0.0947779,6.5];
ybs[1574]=['',1.2916842,-0.2916252,5.7];
ybs[1575]=['ο2 Ori',1.298211,0.2364127,4.07];
ybs[1576]=['',1.292554,-0.2859933,5.72];
ybs[1577]=['62 Eri',1.297715,-0.0897163,5.51];
ybs[1578]=['',1.293039,-0.4484847,6.72];
ybs[1579]=['',1.2898034,-0.6910944,6.1];
ybs[1580]=['',1.3027062,0.2999196,5.48];
ybs[1581]=['99 Tau',1.3048865,0.4185108,5.79];
ybs[1582]=['',1.3386336,1.2878919,6.66];
ybs[1583]=['8 Cam',1.315156,0.9282501,6.08];
ybs[1584]=['',1.3407195,1.2931773,5.96];
ybs[1585]=['98 Tau',1.3064369,0.4377355,5.81];
ybs[1586]=['',1.301707,-0.0180926,6.23];
ybs[1587]=['ω Aur',1.311836,0.661826,4.94];
ybs[1588]=['',1.3240498,1.0665075,6.03];
ybs[1589]=['',1.3304901,1.1667604,6.19];
ybs[1590]=['',1.3032629,-0.2478541,6.15];
ybs[1591]=['',1.3055664,-0.0380889,6.35];
ybs[1592]=['',1.288131,-1.0212848,6.12];
ybs[1593]=['',1.2808328,-1.1631383,6.41];
ybs[1594]=['5 Aur',1.3164951,0.6880747,5.95];
ybs[1595]=['',1.3096739,0.2543385,6.09];
ybs[1596]=['π6 Ori',1.3073053,0.0304412,4.47];
ybs[1597]=['6 Aur',1.3168662,0.6926167,6.58];
ybs[1598]=['β Cam',1.3319079,1.0553949,4.03];
ybs[1599]=['',1.3087481,-0.2852923,5.66];
ybs[1600]=['ε Aur',1.3240231,0.7653543,2.99];
ybs[1601]=['',1.2774289,-1.2631746,6.28];
ybs[1602]=['',1.3113644,-0.2578956,7.71];
ybs[1603]=['63 Eri',1.3125352,-0.1786162,5.38];
ybs[1604]=['',1.3160743,0.063605,7.03];
ybs[1605]=['',1.3161689,0.0636193,6.66];
ybs[1606]=['64 Eri',1.3128484,-0.2183086,4.79];
ybs[1607]=['ζ Aur',1.3260791,0.7173972,3.75];
ybs[1608]=['',1.3164067,-0.0355451,6.32];
ybs[1609]=['',1.3169511,-0.0999101,6.22];
ybs[1610]=['',1.3297313,0.7237751,6.14];
ybs[1611]=['',1.4802721,1.4963568,6.51];
ybs[1612]=['ψ Eri',1.3196151,-0.1247087,4.81];
ybs[1613]=['',1.3216274,0.0131009,5.92];
ybs[1614]=['',1.3223627,0.0285747,6.24];
ybs[1615]=['ι Tau',1.327851,0.377301,4.64];
ybs[1616]=['',1.3191031,-0.3494727,4.91];
ybs[1617]=['11 Cam',1.3436046,1.0297205,5.08];
ybs[1618]=['12 Cam',1.3438801,1.0305684,6.08];
ybs[1619]=['',1.345447,1.0680709,6.04];
ybs[1620]=['',1.3254752,-0.0729854,5.85];
ybs[1621]=['',1.3332291,0.5327076,6.14];
ybs[1622]=['',1.3349433,0.5645664,6.62];
ybs[1623]=['',1.3220686,-0.4580914,5.02];
ybs[1624]=['η Men',1.2854823,-1.3073382,5.47];
ybs[1625]=['',1.3440719,0.9500157,7.24];
ybs[1626]=['',1.3188525,-0.6927114,6.03];
ybs[1627]=['',1.3348074,0.4838594,6.6];
ybs[1628]=['',1.3333595,0.3718458,6.19];
ybs[1629]=['1 Lep',1.3247623,-0.397359,5.75];
ybs[1630]=['',1.3227713,-0.554023,5.94];
ybs[1631]=['',1.3607675,1.2158625,6.41];
ybs[1632]=['9 Aur',1.3451623,0.9010036,5];
ybs[1633]=['11 Ori',1.3340482,0.2693255,4.68];
ybs[1634]=['',1.3412324,0.6276674,6.52];
ybs[1635]=['',1.3299581,-0.250315,6.41];
ybs[1636]=['η Aur',1.3437091,0.7201314,3.17];
ybs[1637]=['',1.3384305,0.3461505,6.44];
ybs[1638]=['',1.3743568,1.2910127,5.43];
ybs[1639]=['',1.3451894,0.7539929,6.2];
ybs[1640]=['',1.3296715,-0.4251678,5.61];
ybs[1641]=['',1.3348957,-0.0525833,6.05];
ybs[1642]=['',1.3602479,1.1334829,6.41];
ybs[1643]=['',1.3371618,0.0210169,6.17];
ybs[1644]=['η1 Pic',1.3236104,-0.8573639,5.38];
ybs[1645]=['',1.3853915,1.3350814,6.37];
ybs[1646]=['',1.328874,-0.7281077,6.31];
ybs[1647]=['γ1 Cae',1.3314308,-0.6188256,4.55];
ybs[1648]=['γ2 Cae',1.3315431,-0.6226996,6.34];
ybs[1649]=['ε Lep',1.3366159,-0.3899836,3.19];
ybs[1650]=['',1.3356259,-0.4559795,5.73];
ybs[1651]=['104 Tau',1.3467435,0.3258643,5];
ybs[1652]=['66 Eri',1.3429273,-0.0807908,5.12];
ybs[1653]=['106 Tau',1.3483676,0.3568117,5.3];
ybs[1654]=['103 Tau',1.349844,0.4239509,5.5];
ybs[1655]=['105 Tau',1.3489342,0.3792624,5.89];
ybs[1656]=['',1.3419787,-0.2285653,6.05];
ybs[1657]=['13 Ori',1.3472343,0.165763,6.17];
ybs[1658]=['η2 Pic',1.3329814,-0.864824,5.03];
ybs[1659]=['14 Ori',1.3482622,0.1487683,5.34];
ybs[1660]=['',1.3455133,-0.2175523,5.97];
ybs[1661]=['β Eri',1.3476611,-0.0883291,2.79];
ybs[1662]=['',1.332739,-0.9491185,6.27];
ybs[1663]=['',1.3624208,0.820064,5.68];
ybs[1664]=['',1.3600727,0.6514643,6.02];
ybs[1665]=['',1.3571673,0.4896533,6.01];
ybs[1666]=['',1.3496613,-0.1507917,5.78];
ybs[1667]=['16 Ori',1.354614,0.1719883,5.43];
ybs[1668]=['68 Eri',1.3515076,-0.0773362,5.12];
ybs[1669]=['ζ Dor',1.3345887,-1.0026216,4.72];
ybs[1670]=['',1.3741871,1.0798833,6.17];
ybs[1671]=['15 Ori',1.3564491,0.2726516,4.82];
ybs[1672]=['β Men',1.3197337,-1.2441781,5.31];
ybs[1673]=['14 Cam',1.3763646,1.0945641,6.5];
ybs[1674]=['λ Eri',1.3531923,-0.1523548,4.27];
ybs[1675]=['',1.3481684,-0.6229596,6.52];
ybs[1676]=['',1.3574286,-0.0094398,6.1];
ybs[1677]=['',1.3052657,-1.3660795,6.29];
ybs[1678]=['',1.3996154,1.279119,5.74];
ybs[1679]=['',1.3651641,0.2804595,5.18];
ybs[1680]=['',1.3613562,-0.0389193,6.25];
ybs[1681]=['',1.4223654,1.3831531,5.05];
ybs[1682]=['',1.3628899,-0.0430577,5.9];
ybs[1683]=['',1.3831617,1.0372018,6.15];
ybs[1684]=['μ Aur',1.3737261,0.6720763,4.86];
ybs[1685]=['',1.3646065,0.0093959,6.67];
ybs[1686]=['',1.3649079,0.0185098,5.89];
ybs[1687]=['',1.3805046,0.9291414,6.2];
ybs[1688]=['',1.3628277,-0.2063917,5.68];
ybs[1689]=['',1.3595055,-0.4517838,6.41];
ybs[1690]=['',1.342725,-1.1060828,5.2];
ybs[1691]=['ι Lep',1.3668336,-0.2067485,4.45];
ybs[1692]=['',1.3692398,-0.1053155,5.91];
ybs[1693]=['ρ Ori',1.3716752,0.0503344,4.46];
ybs[1694]=['',1.3626999,-0.6522562,6.57];
ybs[1695]=['',1.334093,-1.2742843,6.27];
ybs[1696]=['',1.3726701,0.0347456,6.09];
ybs[1697]=['μ Lep',1.3694399,-0.282438,3.31];
ybs[1698]=['',1.3737565,0.0101732,6.32];
ybs[1699]=['',1.3724472,-0.1418088,6.37];
ybs[1700]=['κ Lep',1.3708674,-0.2254702,4.36];
ybs[1701]=['14 Aur',1.382061,0.5708887,5.02];
ybs[1702]=['',1.3917088,0.935616,6.5];
ybs[1703]=['α Aur',1.3884331,0.8031852,0.08];
ybs[1704]=['',1.3780517,0.0903773,5.5];
ybs[1705]=['',1.3741522,-0.2545412,6.21];
ybs[1706]=['108 Tau',1.3818475,0.3893212,6.27];
ybs[1707]=['',1.3860542,0.5992281,5.96];
ybs[1708]=['β Ori',1.3767356,-0.1427578,0.12];
ybs[1709]=['',1.5310215,1.4910277,6.6];
ybs[1710]=['',1.3722774,-0.624878,6.98];
ybs[1711]=['ξ Men',1.2940687,-1.438848,5.85];
ybs[1712]=['',1.3803216,-0.024213,6.15];
ybs[1713]=['18 Ori',1.3840855,0.1983192,5.56];
ybs[1714]=['15 Cam',1.401662,1.0146798,6.13];
ybs[1715]=['',1.4062899,1.0938462,5.61];
ybs[1716]=['',1.3753449,-0.6275311,5.76];
ybs[1717]=['',1.3950927,0.7472197,5.48];
ybs[1718]=['',1.3798051,-0.4698683,5.07];
ybs[1719]=['',1.3864503,0.0343551,6.42];
ybs[1720]=['',1.3967323,0.7065986,6.18];
ybs[1721]=['16 Aur',1.3941908,0.5828012,4.54];
ybs[1722]=['',1.3716802,-0.9077182,6.05];
ybs[1723]=['17 Aur',1.3948067,0.5897037,6.14];
ybs[1724]=['λ Aur',1.3987549,0.7002096,4.71];
ybs[1725]=['',1.3810841,-0.6092065,6.66];
ybs[1726]=['',1.3862891,-0.2988092,6.56];
ybs[1727]=['',1.3977954,0.5893682,5.41];
ybs[1728]=['',1.4029397,0.7757114,6.62];
ybs[1729]=['18 Aur',1.3995237,0.5935051,6.49];
ybs[1730]=['τ Ori',1.3901702,-0.119096,3.6];
ybs[1731]=['',1.4058098,0.8200082,6.54];
ybs[1732]=['',1.3902107,-0.2356019,5.5];
ybs[1733]=['',1.4036364,0.7174255,5.52];
ybs[1734]=['109 Tau',1.3984885,0.3860015,4.94];
ybs[1735]=['19 Aur',1.4022356,0.5930199,5.03];
ybs[1736]=['',1.3982731,0.3517644,6.08];
ybs[1737]=['',1.3793422,-0.9103656,6.49];
ybs[1738]=['ο Col',1.3885086,-0.6086733,4.83];
ybs[1739]=['θ Dor',1.3689696,-1.172205,4.83];
ybs[1740]=['',1.4512711,1.3612153,6.56];
ybs[1741]=['21 Ori',1.3973883,0.0456541,5.34];
ybs[1742]=['',1.3951448,-0.3160762,5.96];
ybs[1743]=['',1.3989963,-0.0243029,6.34];
ybs[1744]=['ρ Aur',1.410489,0.7299494,5.23];
ybs[1745]=['',1.4062031,0.4882775,6.33];
ybs[1746]=['16 Cam',1.4190547,1.0046487,5.28];
ybs[1747]=['',1.4072498,0.5164238,5.76];
ybs[1748]=['',1.3971005,-0.3228869,6.36];
ybs[1749]=['',1.3971591,-0.3227076,6.54];
ybs[1750]=['',1.4056719,0.346155,6.18];
ybs[1751]=['λ Lep',1.3985345,-0.2296307,4.29];
ybs[1752]=['ν Lep',1.4003544,-0.2146051,5.3];
ybs[1753]=['',1.3971931,-0.4773296,5.99];
ybs[1754]=['',1.4025837,-0.0933282,6.39];
ybs[1755]=['',1.4149445,0.7164143,5.54];
ybs[1756]=['',1.4067528,0.0703515,6.57];
ybs[1757]=['',1.4020415,-0.3703599,4.71];
ybs[1758]=['',1.4086733,0.1474333,5.8];
ybs[1759]=['',1.4075098,-0.0069391,5.68];
ybs[1760]=['22 Ori',1.4085217,-0.0063496,4.73];
ybs[1761]=['',1.4009874,-0.6052701,6.34];
ybs[1762]=['ζ Pic',1.3957045,-0.8828942,5.45];
ybs[1763]=['22 Aur',1.4166913,0.5053514,6.46];
ybs[1764]=['',1.4084392,-0.2397634,6.56];
ybs[1765]=['23 Ori',1.4133328,0.0621792,5];
ybs[1766]=['',1.4076645,-0.4320443,5.06];
ybs[1767]=['',1.4050985,-0.5991065,6.09];
ybs[1768]=['σ Aur',1.4226571,0.652801,4.99];
ybs[1769]=['110 Tau',1.4172701,0.29177,6.08];
ybs[1770]=['',1.4222804,0.5452601,6.28];
ybs[1771]=['',1.4222912,0.5438493,5.94];
ybs[1772]=['',1.416383,0.0932062,6.35];
ybs[1773]=['',1.4149924,-0.1465707,5.9];
ybs[1774]=['',1.4249822,0.6086347,6.55];
ybs[1775]=['111 Tau',1.4207649,0.3036993,4.99];
ybs[1776]=['',1.4170088,-0.0024781,5.7];
ybs[1777]=['',1.4176463,-0.0148275,6.11];
ybs[1778]=['8 Lep',1.4156465,-0.2427641,5.25];
ybs[1779]=['29 Ori',1.4177989,-0.1359685,4.14];
ybs[1780]=['',1.4138213,-0.4657898,6.49];
ybs[1781]=['',1.4210148,0.0413655,6.32];
ybs[1782]=['27 Ori',1.4203725,-0.0152547,5.08];
ybs[1783]=['η Ori',1.4202988,-0.0415315,3.36];
ybs[1784]=['ψ1 Ori',1.4216228,0.0325262,4.95];
ybs[1785]=['γ Ori',1.4234584,0.1111208,1.64];
ybs[1786]=['β Tau',1.4293759,0.4995814,1.65];
ybs[1787]=['',1.4197658,-0.2959854,5.65];
ybs[1788]=['',1.4140518,-0.6922085,5.71];
ybs[1789]=['',1.4323814,0.619126,6.15];
ybs[1790]=['',1.431934,0.6005294,5.94];
ybs[1791]=['',1.4320523,0.5808263,6.15];
ybs[1792]=['',1.4152737,-0.6513362,6.82];
ybs[1793]=['113 Tau',1.428034,0.2917633,6.25];
ybs[1794]=['',1.4224222,-0.1799794,5.61];
ybs[1795]=['',1.4249149,-0.0092034,6.57];
ybs[1796]=['κ Pic',1.408254,-0.9794069,6.11];
ybs[1797]=['17 Cam',1.4491138,1.1009809,5.42];
ybs[1798]=['',1.4261006,0.0093821,6.16];
ybs[1799]=['',1.4331434,0.5275187,5.74];
ybs[1800]=['φ Aur',1.4355752,0.6019912,5.07];
ybs[1801]=['',1.4270122,-0.0960233,6.23];
ybs[1802]=['',1.4300872,0.1201737,6.42];
ybs[1803]=['115 Tau',1.4327653,0.3137791,5.42];
ybs[1804]=['',1.4329329,0.2665773,6.16];
ybs[1805]=['114 Tau',1.4349567,0.383147,4.88];
ybs[1806]=['ψ2 Ori',1.4307854,0.0543102,4.59];
ybs[1807]=['',1.4263082,-0.3434617,5.89];
ybs[1808]=['',1.4204285,-0.7715854,6.08];
ybs[1809]=['116 Tau',1.4352686,0.2773307,5.5];
ybs[1810]=['',1.3547334,-1.4227518,6.51];
ybs[1811]=['117 Tau',1.4364833,0.3011473,5.77];
ybs[1812]=['',1.4313248,-0.2074277,6.35];
ybs[1813]=['θ Pic',1.419113,-0.9127902,6.27];
ybs[1814]=['',1.4387658,0.2390088,6.35];
ybs[1815]=['',1.4359157,0.0229326,6.41];
ybs[1816]=['118 Tau',1.4422451,0.4392208,5.47];
ybs[1817]=['',1.4441739,0.5096558,6.24];
ybs[1818]=['',1.4332737,-0.372797,6.07];
ybs[1819]=['',1.4497622,0.7238942,6];
ybs[1820]=['',1.4494095,0.6953393,6.37];
ybs[1821]=['',1.4397639,-0.0574621,6.39];
ybs[1822]=['',1.4300663,-0.7143186,5.87];
ybs[1823]=['18 Cam',1.4587391,0.9989271,6.48];
ybs[1824]=['β Lep',1.4360757,-0.3620494,2.84];
ybs[1825]=['',1.4417154,-0.05989,5.79];
ybs[1826]=['',1.4484547,0.392293,6.29];
ybs[1827]=['',1.4469222,0.2683386,5.94];
ybs[1828]=['',1.4441649,0.031483,5.78];
ybs[1829]=['31 Ori',1.4432801,-0.0188051,4.71];
ybs[1830]=['',1.435368,-0.6495238,5.57];
ybs[1831]=['λ Dor',1.425164,-1.0279267,5.14];
ybs[1832]=['',1.4460664,0.0736291,6.21];
ybs[1833]=['',1.43946,-0.5253707,6.75];
ybs[1834]=['32 Ori',1.4481124,0.1040617,4.2];
ybs[1835]=['',1.4457289,-0.1295076,6.33];
ybs[1836]=['33 Ori',1.4500128,0.0577048,5.46];
ybs[1837]=['χ Aur',1.4576315,0.5620861,4.76];
ybs[1838]=['',1.4943629,1.3099269,6.17];
ybs[1839]=['119 Tau',1.4548042,0.32477,4.38];
ybs[1840]=['',1.4614476,0.7351625,6.55];
ybs[1841]=['',1.4548475,0.2979548,5.46];
ybs[1842]=['',1.450132,-0.1168384,6.22];
ybs[1843]=['10 Lep',1.448644,-0.3638922,5.55];
ybs[1844]=['',1.4608469,0.5727118,6.48];
ybs[1845]=['δ Ori',1.4532345,-0.0047263,6.85];
ybs[1846]=['δ Ori',1.4532267,-0.0049832,2.23];
ybs[1847]=['',1.4808073,1.1642548,6.26];
ybs[1848]=['',1.4617067,0.6063029,6.27];
ybs[1849]=['υ Ori',1.4526491,-0.1271942,4.62];
ybs[1850]=['',1.4430484,-0.821406,5.46];
ybs[1851]=['19 Cam',1.4802171,1.1198999,6.15];
ybs[1852]=['120 Tau',1.4605331,0.3238133,5.69];
ybs[1853]=['',1.4263573,-1.1974071,6.03];
ybs[1854]=['',1.4611319,0.3575649,6.18];
ybs[1855]=['',1.4561564,-0.0275524,5.35];
ybs[1856]=['ε Col',1.4483599,-0.6188316,3.87];
ybs[1857]=['',1.4580355,-0.029762,6.46];
ybs[1858]=['35 Ori',1.4620249,0.2499003,5.64];
ybs[1859]=['α Lep',1.4557572,-0.3108239,2.58];
ybs[1860]=['',1.4759932,0.9501592,5.73];
ybs[1861]=['',1.4376218,-1.0873269,6.59];
ybs[1862]=['',1.4598148,-0.0199528,5.34];
ybs[1863]=['',1.4739971,0.8329879,6.11];
ybs[1864]=['',1.4494283,-0.8013036,5.86];
ybs[1865]=['',1.4618087,0.0247917,6.59];
ybs[1866]=['38 Ori',1.4632794,0.065964,5.36];
ybs[1867]=['',1.4621824,-0.0178535,6.22];
ybs[1868]=['',1.4621746,-0.0254456,5.93];
ybs[1869]=['121 Tau',1.4691524,0.419775,5.38];
ybs[1870]=['φ1 Ori',1.4658436,0.1658355,4.41];
ybs[1871]=['',1.4553643,-0.6719472,5.48];
ybs[1872]=['',1.4713727,0.4830002,6.27];
ybs[1873]=['λ Ori',1.4672484,0.1735947,3.54];
ybs[1874]=['λ Ori',1.467263,0.1736092,5.61];
ybs[1875]=['',1.4567168,-0.613074,5.78];
ybs[1876]=['',1.4416017,-1.115493,6.19];
ybs[1877]=['',1.4676302,0.1789317,5.6];
ybs[1878]=['',1.4761462,0.7015066,6.09];
ybs[1879]=['',1.6048366,1.4814178,6.11];
ybs[1880]=['',1.4661549,-0.1046673,5.67];
ybs[1881]=['',1.466286,-0.1045416,4.78];
ybs[1882]=['',1.4602243,-0.5207382,6.53];
ybs[1883]=['',1.4738216,0.4529271,6.49];
ybs[1884]=['',1.4677276,-0.0782142,6.56];
ybs[1885]=['',1.4677808,-0.0770265,6.24];
ybs[1886]=['42 Ori',1.4678174,-0.0842357,4.59];
ybs[1887]=['θ1 Ori',1.4672674,-0.0938146,6.73];
ybs[1888]=['θ1 Ori',1.467282,-0.0937806,7.96];
ybs[1889]=['θ1 Ori',1.467311,-0.0938583,5.13];
ybs[1890]=['θ1 Ori',1.4673691,-0.0938245,6.7];
ybs[1891]=['θ2 Ori',1.4677754,-0.0943198,5.08];
ybs[1892]=['',1.4684083,-0.0759708,6.38];
ybs[1893]=['ι Ori',1.4679836,-0.1029402,2.77];
ybs[1894]=['',1.4687888,-0.0565644,6.4];
ybs[1895]=['45 Ori',1.4690021,-0.0845435,5.26];
ybs[1896]=['',1.4766714,0.4701082,5.83];
ybs[1897]=['ε Ori',1.4715505,-0.020776,1.7];
ybs[1898]=['',1.4796636,0.5859051,6.33];
ybs[1899]=['122 Tau',1.4759085,0.2976029,5.54];
ybs[1900]=['',1.4715561,-0.0983755,6.54];
ybs[1901]=['φ2 Ori',1.4749343,0.1623463,4.09];
ybs[1902]=['',1.4757311,0.1927911,5.94];
ybs[1903]=['',1.4661277,-0.5771435,5.78];
ybs[1904]=['ζ Tau',1.4786057,0.3691949,3];
ybs[1905]=['',1.4730467,-0.1056555,5.72];
ybs[1906]=['',1.4579844,-0.9579985,6.43];
ybs[1907]=['',1.4767404,0.1564328,6.12];
ybs[1908]=['26 Aur',1.4833268,0.5323742,5.4];
ybs[1909]=['',1.4703251,-0.5008419,6.26];
ybs[1910]=['',1.5031614,1.146786,5.6];
ybs[1911]=['',1.4534461,-1.1207476,5.34];
ybs[1912]=['',1.4768109,-0.1034522,6.05];
ybs[1913]=['',1.4752512,-0.205328,6.11];
ybs[1914]=['',1.4797304,0.1318081,5.88];
ybs[1915]=['',1.4845565,0.464745,6.37];
ybs[1916]=['β Dor',1.4564995,-1.0904234,3.76];
ybs[1917]=['',1.4787409,-0.0838259,6.19];
ybs[1918]=['',1.486199,0.5100768,5.96];
ybs[1919]=['',1.4966362,0.9335763,6.23];
ybs[1920]=['ν1 Col',1.4751761,-0.4862538,6.16];
ybs[1921]=['',1.4687096,-0.8255777,6.11];
ybs[1922]=['125 Tau',1.4879291,0.4521573,5.18];
ybs[1923]=['',1.4865068,0.3800051,6.34];
ybs[1924]=['',1.463194,-1.0272795,6.75];
ybs[1925]=['σ Ori',1.4825557,-0.0451986,3.81];
ybs[1926]=['',1.4827231,-0.045097,6.65];
ybs[1927]=['',1.4819081,-0.1145549,5.96];
ybs[1928]=['ω Ori',1.4847014,0.0721078,4.57];
ybs[1929]=['ν2 Col',1.4771829,-0.5005356,5.31];
ybs[1930]=['',1.4625258,-1.0675035,6.32];
ybs[1931]=['49 Ori',1.4829983,-0.1257126,4.8];
ybs[1932]=['',1.4919351,0.547464,6.04];
ybs[1933]=['',1.4924129,0.5572854,6.11];
ybs[1934]=['',1.4858965,-0.0620428,6];
ybs[1935]=['24 Cam',1.5043325,0.9876762,6.05];
ybs[1936]=['',1.4856516,-0.1692397,6.5];
ybs[1937]=['23 Cam',1.509851,1.0730997,6.15];
ybs[1938]=['',1.4842976,-0.3113556,6.38];
ybs[1939]=['',1.495128,0.5148103,6.43];
ybs[1940]=['126 Tau',1.4943545,0.2887281,4.86];
ybs[1941]=['',1.480852,-0.7102981,5.82];
ybs[1942]=['ζ Ori',1.491356,-0.0337453,2.05];
ybs[1943]=['ζ Ori',1.4913633,-0.0337454,4.21];
ybs[1944]=['',1.490729,-0.0491418,6.22];
ybs[1945]=['',1.4973049,0.407274,6.59];
ybs[1946]=['',1.4917551,-0.019541,4.95];
ybs[1947]=['γ Men',1.4445624,-1.3321564,5.19];
ybs[1948]=['',1.497954,0.3956468,6.36];
ybs[1949]=['',1.4928967,0.0060549,5.93];
ybs[1950]=['α Col',1.4852379,-0.594533,2.64];
ybs[1951]=['',1.491095,-0.1815163,6.52];
ybs[1952]=['',1.4861022,-0.5693146,5.45];
ybs[1953]=['',1.495308,-0.050392,6.42];
ybs[1954]=['',1.4700691,-1.1614957,6.31];
ybs[1955]=['',1.5034743,0.4051285,6.21];
ybs[1956]=['',1.494899,-0.291761,6.21];
ybs[1957]=['51 Ori',1.4989721,0.0258864,4.91];
ybs[1958]=['',1.4583999,-1.2868088,5.78];
ybs[1959]=['',1.4972536,-0.3058108,6.15];
ybs[1960]=['',1.493132,-0.5827919,6.34];
ybs[1961]=['',1.5005168,-0.1184703,6.02];
ybs[1962]=['12 Lep',1.4970343,-0.3903427,5.87];
ybs[1963]=['26 Cam',1.5194669,0.9795109,5.94];
ybs[1964]=['',1.5018247,-0.0280113,6.31];
ybs[1965]=['ο Aur',1.5162136,0.8697501,5.47];
ybs[1966]=['',1.4965131,-0.5327948,6.19];
ybs[1967]=['',1.4965656,-0.604916,5.29];
ybs[1968]=['',1.5152296,0.7071013,6.58];
ybs[1969]=['',1.5021081,-0.3237489,5.73];
ybs[1970]=['',1.5316145,1.0962938,6.13];
ybs[1971]=['',1.5135474,0.3613152,6.95];
ybs[1972]=['',1.5102086,0.0700791,6.09];
ybs[1973]=['',1.5215512,0.7423348,6.29];
ybs[1974]=['',1.5069047,-0.3511408,6.34];
ybs[1975]=['',1.5017681,-0.6876406,6.25];
ybs[1976]=['',1.5066748,-0.3912005,6.15];
ybs[1977]=['γ Lep',1.5067682,-0.3916661,3.6];
ybs[1978]=['',1.5021477,-0.7997985,6.39];
ybs[1979]=['129 Tau',1.5181669,0.2762646,6];
ybs[1980]=['',1.5143495,-0.0743797,6.34];
ybs[1981]=['',1.5184168,0.1663033,5.79];
ybs[1982]=['',1.5168802,0.0204983,5.95];
ybs[1983]=['131 Tau',1.5201316,0.252975,5.72];
ybs[1984]=['130 Tau',1.5211989,0.3095362,5.49];
ybs[1985]=['ι Men',1.458924,-1.3754661,6.05];
ybs[1986]=['29 Cam',1.537276,0.9934964,6.54];
ybs[1987]=['133 Tau',1.5222697,0.2426976,5.29];
ybs[1988]=['',1.5495583,1.1951031,6.2];
ybs[1989]=['τ Aur',1.5297699,0.6839273,4.52];
ybs[1990]=['μ Col',1.5130298,-0.5637346,5.17];
ybs[1991]=['',1.5254085,0.3643362,6.07];
ybs[1992]=['ζ Lep',1.517935,-0.2585824,3.55];
ybs[1993]=['52 Ori',1.5232606,0.1127459,5.27];
ybs[1994]=['',1.5186311,-0.2832949,6.17];
ybs[1995]=['',1.5202295,-0.1837316,6.03];
ybs[1996]=['132 Tau',1.5283667,0.4288738,4.86];
ybs[1997]=['',1.5383531,0.8991732,6.29];
ybs[1998]=['κ Ori',1.5216203,-0.1686664,2.06];
ybs[1999]=['',1.5179133,-0.499739,6.22];
ybs[2000]=['30 Cam',1.5450487,1.0291779,6.14];
ybs[2001]=['',1.5254165,-0.0713716,5.97];
ybs[2002]=['',1.5141844,-0.8131598,5.31];
ybs[2003]=['',1.5185825,-0.6225343,6.32];
ybs[2004]=['134 Tau',1.5302235,0.2208893,4.91];
ybs[2005]=['υ Aur',1.5378023,0.6511766,4.74];
ybs[2006]=['ν Aur',1.5398696,0.68334,3.97];
ybs[2007]=['',1.5370297,0.4882027,5.56];
ybs[2008]=['',1.532283,0.172365,5.8];
ybs[2009]=['δ Dor',1.5045335,-1.1471698,4.35];
ybs[2010]=['135 Tau',1.5343569,0.2497566,5.52];
ybs[2011]=['',1.5211747,-0.7094183,6.61];
ybs[2012]=['',1.5392299,0.5607557,6.25];
ybs[2013]=['',1.5328543,0.077282,5.97];
ybs[2014]=['β Pic',1.5174394,-0.8911683,3.85];
ybs[2015]=['',1.5295246,-0.2527005,5.49];
ybs[2016]=['π Men',1.4639435,-1.4042466,5.65];
ybs[2017]=['',1.5168292,-0.9486663,6.18];
ybs[2018]=['',1.5339997,0.0354111,5.98];
ybs[2019]=['',1.5449925,0.690762,6.45];
ybs[2020]=['',1.5304313,-0.4008471,5.87];
ybs[2021]=['31 Cam',1.5568384,1.0452844,5.2];
ybs[2022]=['',1.5447293,0.5920298,5.98];
ybs[2023]=['ξ Aur',1.5558156,0.9723069,4.99];
ybs[2024]=['',1.5428954,0.3468239,6.06];
ybs[2025]=['55 Ori',1.5374481,-0.1311441,5.35];
ybs[2026]=['',1.5278829,-0.7831333,6.38];
ybs[2027]=['137 Tau',1.5425987,0.2474034,5.59];
ybs[2028]=['136 Tau',1.5473044,0.4819766,4.58];
ybs[2029]=['δ Lep',1.5367552,-0.3643386,3.81];
ybs[2030]=['',1.5373489,-0.4000658,6.17];
ybs[2031]=['56 Ori',1.5424586,0.032437,4.78];
ybs[2032]=['',1.5469631,0.3543401,6.71];
ybs[2033]=['',1.5407101,-0.1577378,5.97];
ybs[2034]=['β Col',1.5345025,-0.6241998,3.12];
ybs[2035]=['',1.5693303,1.1536068,6.25];
ybs[2036]=['γ Pic',1.5280299,-0.9802065,4.51];
ybs[2037]=['',1.5393139,-0.5139091,6.45];
ybs[2038]=['',1.5312085,-0.9208907,6.35];
ybs[2039]=['',1.5615061,0.904174,6.49];
ybs[2040]=['',1.5547216,0.5533363,5.9];
ybs[2041]=['χ1 Ori',1.5516019,0.3539284,4.41];
ybs[2042]=['',1.5505476,0.1848223,6.12];
ybs[2043]=['',1.5330589,-0.9093945,5.17];
ybs[2044]=['',1.5519572,0.2053369,6.59];
ybs[2045]=['',1.5504548,0.0563369,6.31];
ybs[2046]=['57 Ori',1.5540392,0.3447364,5.92];
ybs[2047]=['',1.5413537,-0.6567251,5.63];
ybs[2048]=['',1.5649409,0.8557437,6.47];
ybs[2049]=['',1.5423579,-0.672343,6.7];
ybs[2050]=['λ Col',1.5440063,-0.5898889,4.87];
ybs[2051]=['',1.5524415,0.0169418,6];
ybs[2052]=['',1.5515828,-0.0708855,6.57];
ybs[2053]=['',1.4239086,-1.4795026,6.2];
ybs[2054]=['',1.5483012,-0.3427049,6.69];
ybs[2055]=['α Ori',1.5545721,0.1293127,0.5];
ybs[2056]=['λ Men',1.5157699,-1.2687854,6.53];
ybs[2057]=['',1.5578812,0.3521512,5.4];
ybs[2058]=['ε Dor',1.526617,-1.1675568,5.11];
ybs[2059]=['',1.5519565,-0.2054514,5.66];
ybs[2060]=['',1.5614869,0.5051613,6.32];
ybs[2061]=['',1.5586759,0.2430711,6.6];
ybs[2062]=['',1.5491157,-0.5086781,6.36];
ybs[2063]=['',1.5446512,-0.7490648,6.55];
ybs[2064]=['',1.5556069,-0.0805365,5.87];
ybs[2065]=['',1.5559717,-0.0835383,6.28];
ybs[2066]=['',1.5388666,-0.9974972,5.94];
ybs[2067]=['',1.5336965,-1.1175229,6.36];
ybs[2068]=['',1.5629044,0.4232587,6.02];
ybs[2069]=['',1.5602974,0.1660018,5.99];
ybs[2070]=['',1.5619333,0.2011039,5.87];
ybs[2071]=['δ Aur',1.5760859,0.9474439,3.72];
ybs[2072]=['',1.6056198,1.3191646,6.4];
ybs[2073]=['',1.5772198,0.9655253,6.44];
ybs[2074]=['',1.5773128,0.952023,6.14];
ybs[2075]=['',1.5749703,0.8713445,5.89];
ybs[2076]=['',1.5513502,-0.697358,5.57];
ybs[2077]=['',1.5475752,-0.8789333,6.52];
ybs[2078]=['139 Tau',1.5676019,0.4529927,4.82];
ybs[2079]=['η Lep',1.5591919,-0.2472472,3.71];
ybs[2080]=['',1.5581286,-0.398609,5.96];
ybs[2081]=['ξ Col',1.5542312,-0.6478447,4.97];
ybs[2082]=['β Aur',1.575311,0.7844795,1.9];
ybs[2083]=['',1.5498477,-0.8661096,6.1];
ybs[2084]=['',1.5595821,-0.4051618,6.36];
ybs[2085]=['π Aur',1.5771553,0.8017448,4.26];
ybs[2086]=['σ Col',1.5582325,-0.5476996,5.5];
ybs[2087]=['',1.5662749,0.0213842,6.22];
ybs[2088]=['',1.5502435,-0.918616,5.29];
ybs[2089]=['θ Aur',1.575679,0.6494769,2.62];
ybs[2090]=['',1.5787262,0.778267,6.22];
ybs[2091]=['',1.5674704,-0.0173403,6.22];
ybs[2092]=['',1.5602474,-0.5580641,6.44];
ybs[2093]=['',1.5709758,0.2235571,5.7];
ybs[2094]=['59 Ori',1.5684924,0.0320699,5.9];
ybs[2095]=['36 Aur',1.5818468,0.8360314,5.73];
ybs[2096]=['',1.5457213,-1.1010776,4.65];
ybs[2097]=['60 Ori',1.5702804,0.0096582,5.22];
ybs[2098]=['',1.545895,-1.1253771,6.63];
ybs[2099]=['',1.5851664,0.8544817,5.96];
ybs[2100]=['γ Col',1.563223,-0.6157921,4.36];
ybs[2101]=['1 Mon',1.570761,-0.1637464,6.12];
ybs[2102]=['2 Mon',1.5709947,-0.1668206,5.03];
ybs[2103]=['',1.5737088,-0.0252116,6.63];
ybs[2104]=['',1.5816843,0.5416422,5.98];
ybs[2105]=['',1.5808132,0.4812117,6.05];
ybs[2106]=['',1.5900145,0.8709849,6.05];
ybs[2107]=['',1.5755211,-0.0536593,4.53];
ybs[2108]=['',1.5606536,-0.9324394,6.45];
ybs[2109]=['',1.5898873,0.757068,6.42];
ybs[2110]=['',1.5835903,0.3909482,6.37];
ybs[2111]=['',1.5674626,-0.7685367,5.81];
ybs[2112]=['',1.5762117,-0.2251493,6.22];
ybs[2113]=['38 Aur',1.5916312,0.7489148,6.1];
ybs[2114]=['η Col',1.5698119,-0.7472629,3.96];
ybs[2115]=['',1.6011668,1.0365518,6.34];
ybs[2116]=['',1.5893899,0.5695715,6.24];
ybs[2117]=['',1.5974602,0.9000785,6.45];
ybs[2118]=['μ Ori',1.5861186,0.168355,4.12];
ybs[2119]=['κ Men',1.5223713,-1.3850271,5.47];
ybs[2120]=['',1.608396,1.1074081,6.39];
ybs[2121]=['',1.5854098,0.0295491,6.59];
ybs[2122]=['3 Mon',1.5830365,-0.1849912,4.95];
ybs[2123]=['',1.5797501,-0.4436381,6.05];
ybs[2124]=['64 Ori',1.5911687,0.3436297,5.14];
ybs[2125]=['',1.5795976,-0.5918845,5.55];
ybs[2126]=['39 Aur',1.5993003,0.7501212,5.87];
ybs[2127]=['',1.5906787,0.2038344,6.08];
ybs[2128]=['1 Gem',1.5942138,0.4059806,4.16];
ybs[2129]=['χ2 Ori',1.5932153,0.3514408,4.63];
ybs[2130]=['',1.5860504,-0.2530506,6.2];
ybs[2131]=['',1.5989401,0.6625498,6.34];
ybs[2132]=['',1.5764798,-0.8939038,5.67];
ybs[2133]=['',1.6009928,0.5863618,6.23];
ybs[2134]=['',1.5886142,-0.458782,5.04];
ybs[2135]=['',1.6035931,0.6175691,6.12];
ybs[2136]=['',1.5935815,-0.1171381,5.21];
ybs[2137]=['40 Aur',1.6056958,0.6715878,5.36];
ybs[2138]=['63 Ori',1.5972553,0.0945489,5.67];
ybs[2139]=['66 Ori',1.5972258,0.0725336,5.63];
ybs[2140]=['',1.6043239,0.5150292,6.08];
ybs[2141]=['',1.6096634,0.7304219,6.12];
ybs[2142]=['17 Lep',1.5965426,-0.2877551,4.93];
ybs[2143]=['',1.5930575,-0.5615571,5.65];
ybs[2144]=['',1.598801,-0.1788218,5.87];
ybs[2145]=['',1.5813129,-1.0489094,6.45];
ybs[2146]=['37 Cam',1.6222512,1.0285297,5.36];
ybs[2147]=['',1.6136934,0.7164755,6.36];
ybs[2148]=['',1.604229,-0.0732594,5.38];
ybs[2149]=['θ Lep',1.6017059,-0.2607273,4.67];
ybs[2150]=['',1.599626,-0.4223458,6.95];
ybs[2151]=['',1.5929052,-0.7860796,6.35];
ybs[2152]=['',1.5937604,-0.7868182,5.93];
ybs[2153]=['ν Ori',1.6089431,0.2576851,4.42];
ybs[2154]=['',1.597747,-0.6198799,5.8];
ybs[2155]=['',1.6049416,-0.1950801,6.66];
ybs[2156]=['',1.593997,-0.8458062,6.58];
ybs[2157]=['',1.6030381,-0.4034158,5.47];
ybs[2158]=['',1.6008221,-0.5194421,5.81];
ybs[2159]=['36 Cam',1.6358633,1.1468807,5.32];
ybs[2160]=['',1.6049526,-0.3807691,5.78];
ybs[2161]=['',1.6140232,0.1512388,6.55];
ybs[2162]=['19 Lep',1.6082573,-0.3345776,5.31];
ybs[2163]=['',1.6178184,0.3872003,5.93];
ybs[2164]=['',1.6048343,-0.598921,5.83];
ybs[2165]=['π1 Col',1.6027306,-0.7383108,6.12];
ybs[2166]=['',1.6293281,0.918758,6.3];
ybs[2167]=['3 Gem',1.6186996,0.4033138,5.75];
ybs[2168]=['',1.6145847,0.043541,5.73];
ybs[2169]=['41 Aur',1.6282965,0.8500955,6.82];
ybs[2170]=['41 Aur',1.6283036,0.8500616,6.09];
ybs[2171]=['θ Col',1.606708,-0.650257,5.02];
ybs[2172]=['',1.6040648,-0.7870569,6.51];
ybs[2173]=['',1.6170841,-0.0997655,6.17];
ybs[2174]=['',1.6136764,-0.3915153,5.5];
ybs[2175]=['π2 Col',1.6079688,-0.7357955,5.5];
ybs[2176]=['',1.6154754,-0.3164454,6.35];
ybs[2177]=['',1.6166272,-0.2546289,5.56];
ybs[2178]=['',1.6241595,0.3163175,6.33];
ybs[2179]=['5 Gem',1.6266293,0.4261086,5.8];
ybs[2180]=['',1.6172985,-0.3975731,5.71];
ybs[2181]=['',1.6108572,-0.7742374,6.27];
ybs[2182]=['',1.637868,0.8930021,6.04];
ybs[2183]=['',1.6304786,0.5704933,5.78];
ybs[2184]=['',1.6279134,0.3815709,6.56];
ybs[2185]=['',1.6258888,0.2379341,6.04];
ybs[2186]=['',1.6170821,-0.4661061,6.27];
ybs[2187]=['68 Ori',1.6285507,0.3453007,5.75];
ybs[2188]=['η1 Dor',1.5977633,-1.1526642,5.71];
ybs[2189]=['',1.6232361,-0.1179825,6.15];
ybs[2190]=['',1.602402,-1.0848666,5.05];
ybs[2191]=['6 Gem',1.6299699,0.3997134,6.39];
ybs[2192]=['69 Ori',1.6285557,0.2814214,4.95];
ybs[2193]=['ξ Ori',1.6279823,0.2478831,4.48];
ybs[2194]=['',1.6205095,-0.474025,5.72];
ybs[2195]=['40 Cam',1.6472235,1.0470389,5.35];
ybs[2196]=['',1.6263928,-0.0815356,6.18];
ybs[2197]=['',1.6180722,-0.7043947,5.58];
ybs[2198]=['',1.6140096,-0.8651175,6.49];
ybs[2199]=['',1.6269158,-0.1144313,5.05];
ybs[2200]=['',1.6233678,-0.4623031,6.09];
ybs[2201]=['',1.6352049,0.3259094,6.58];
ybs[2202]=['',1.6333932,0.1853652,6.45];
ybs[2203]=['',1.6628704,1.209684,4.8];
ybs[2204]=['',1.6308824,-0.0438259,6.62];
ybs[2205]=['',1.6199061,-0.7904142,6.31];
ybs[2206]=['δ Pic',1.617477,-0.9594744,4.81];
ybs[2207]=['',1.6304705,-0.3101387,6.52];
ybs[2208]=['',1.6391892,0.3123945,5.88];
ybs[2209]=['1 Lyn',1.6572003,1.0734804,4.98];
ybs[2210]=['η Gem',1.6411158,0.3926808,3.28];
ybs[2211]=['',1.645118,0.6307704,6.92];
ybs[2212]=['',1.6359226,-0.0654248,5.83];
ybs[2213]=['κ Aur',1.6436054,0.514699,4.35];
ybs[2214]=['71 Ori',1.6408591,0.3342079,5.2];
ybs[2215]=['ν Dor',1.6083577,-1.2016182,5.06];
ybs[2216]=['',1.6419386,0.241611,5.91];
ybs[2217]=['72 Ori',1.6432309,0.2816105,5.3];
ybs[2218]=['',1.6389772,-0.0798638,5.83];
ybs[2219]=['',1.6345042,-0.4165926,6.39];
ybs[2220]=['',1.6334053,-0.5131801,6.54];
ybs[2221]=['γ Mon',1.6399793,-0.109648,3.98];
ybs[2222]=['42 Aur',1.6541672,0.8100953,6.52];
ybs[2223]=['73 Ori',1.6445453,0.2189164,5.33];
ybs[2224]=['8 Gem',1.6474587,0.4182083,6.08];
ybs[2225]=['',1.6439578,0.1057328,6.07];
ybs[2226]=['',1.6443972,0.0746213,6.64];
ybs[2227]=['',1.643307,-0.0090798,5.65];
ybs[2228]=['',1.6428194,-0.0859268,5.99];
ybs[2229]=['',1.6475388,0.2997242,6.39];
ybs[2230]=['',1.644798,0.0202631,6.37];
ybs[2231]=['',1.6424052,-0.1578336,6.1];
ybs[2232]=['2 Lyn',1.6643063,1.0297551,4.48];
ybs[2233]=['43 Aur',1.6572308,0.808979,6.38];
ybs[2234]=['9 Gem',1.6503361,0.4142028,6.25];
ybs[2235]=['74 Ori',1.6475603,0.2140429,5.04];
ybs[2236]=['',1.6407072,-0.3539523,5.91];
ybs[2237]=['',1.6414391,-0.3226154,5.99];
ybs[2238]=['',1.6436107,-0.2395709,5.01];
ybs[2239]=['η2 Dor',1.6200838,-1.1448497,5.01];
ybs[2240]=['',1.646773,0.0187078,6.63];
ybs[2241]=['75 Ori',1.6503852,0.1733759,5.39];
ybs[2242]=['',1.6496865,0.1229468,6.57];
ybs[2243]=['',1.6451608,-0.2901789,5.92];
ybs[2244]=['',1.652476,0.2452067,6.59];
ybs[2245]=['',1.6509051,0.0888619,5.71];
ybs[2246]=['',1.6438528,-0.5200465,6.67];
ybs[2247]=['',1.6548368,0.2508646,6.16];
ybs[2248]=['',1.6489801,-0.3966034,6.07];
ybs[2249]=['6 Mon',1.6517408,-0.1873483,6.75];
ybs[2250]=['κ Col',1.6462026,-0.6134653,4.37];
ybs[2251]=['4 Lyn',1.6749824,1.0360415,5.94];
ybs[2252]=['',1.6590337,0.3022079,6.32];
ybs[2253]=['',1.6571817,0.1577368,6.24];
ybs[2254]=['',1.6519895,-0.2936493,5.14];
ybs[2255]=['α Men',1.612759,-1.3047732,5.09];
ybs[2256]=['',1.6461587,-0.6854408,6];
ybs[2257]=['',1.6481071,-0.6587944,5.53];
ybs[2258]=['45 Aur',1.6730367,0.9327212,5.36];
ybs[2259]=['',1.6487375,-0.6503405,5.87];
ybs[2260]=['',1.6541954,-0.348651,5.52];
ybs[2261]=['',1.6572719,-0.1640591,5.36];
ybs[2262]=['',1.6569365,-0.2623983,6.06];
ybs[2263]=['',1.6634711,0.2555308,5.69];
ybs[2264]=['',1.6585586,-0.1500311,6.22];
ybs[2265]=['',1.657452,-0.3653981,5.81];
ybs[2266]=['',1.6690139,0.5154,6.43];
ybs[2267]=['7 Mon',1.6611236,-0.1367135,5.27];
ybs[2268]=['',1.6432041,-1.0336154,6.43];
ybs[2269]=['',1.6625159,-0.0515683,4.9];
ybs[2270]=['',1.6668632,0.2050016,6.54];
ybs[2271]=['',1.669529,0.3098423,6.35];
ybs[2272]=['',1.6507068,-0.9205224,6.41];
ybs[2273]=['',1.6599068,-0.6005092,5.78];
ybs[2274]=['',1.6689651,0.0394041,6.31];
ybs[2275]=['',1.6549275,-0.8790985,7.04];
ybs[2276]=['ζ CMa',1.6628659,-0.5248839,3.02];
ybs[2277]=['',1.6352471,-1.2515793,6.64];
ybs[2278]=['',1.6683933,-0.2056734,5.64];
ybs[2279]=['',1.7042031,1.2308228,5.97];
ybs[2280]=['μ Gem',1.6763837,0.392732,2.88];
ybs[2281]=['',1.6744699,0.2191867,6];
ybs[2282]=['',1.6639416,-0.5961054,5.53];
ybs[2283]=['ψ1 Aur',1.6863201,0.8600161,4.91];
ybs[2284]=['',1.6608287,-0.8508697,6.6];
ybs[2285]=['',1.6936666,0.9821219,5.64];
ybs[2286]=['',1.6772127,0.0654951,6.4];
ybs[2287]=['5 Lyn',1.6955987,1.0193327,5.21];
ybs[2288]=['β CMa',1.6737858,-0.3135893,1.98];
ybs[2289]=['',1.6772253,-0.0820147,6.67];
ybs[2290]=['δ Col',1.6705591,-0.5837702,3.85];
ybs[2291]=['',1.6850549,0.5182674,6.71];
ybs[2292]=['ε Mon',1.6792485,0.0799482,4.44];
ybs[2293]=['',1.6792777,0.0799966,6.72];
ybs[2294]=['',1.6805771,0.1548543,6.26];
ybs[2295]=['',1.6780102,-0.1725506,6.19];
ybs[2296]=['',1.6845171,0.2800256,6.33];
ybs[2297]=['',1.6785488,-0.2632604,6.24];
ybs[2298]=['',1.6877096,0.406905,6.06];
ybs[2299]=['',1.6804453,-0.2014551,5.22];
ybs[2300]=['',1.6784859,-0.3455285,6.6];
ybs[2301]=['',1.6755561,-0.5550452,6.34];
ybs[2302]=['',1.6870271,0.2567204,6.24];
ybs[2303]=['',1.6811353,-0.2264536,6.12];
ybs[2304]=['',1.6856624,0.1234476,5.98];
ybs[2305]=['',1.678854,-0.4466276,5.63];
ybs[2306]=['',1.6858536,0.0259752,6.66];
ybs[2307]=['',1.6856235,-0.0167365,5.87];
ybs[2308]=['',1.6990557,0.8271301,6.56];
ybs[2309]=['',1.687931,0.0394296,6.51];
ybs[2310]=['',1.6786961,-0.640883,5.62];
ybs[2311]=['',1.6877543,-0.0681068,6.35];
ybs[2312]=['',1.6822052,-0.5025238,6.39];
ybs[2313]=['',1.6970342,0.5680875,6.43];
ybs[2314]=['ν Pic',1.6724516,-0.9840427,5.61];
ybs[2315]=['',1.6884654,-0.1380088,6.4];
ybs[2316]=['',1.675926,-0.9109392,5.98];
ybs[2317]=['',1.681658,-0.7033088,6.31];
ybs[2318]=['',1.6916547,-0.0265417,5.87];
ybs[2319]=['',1.6911768,-0.0804715,6.15];
ybs[2320]=['α Car',1.6772892,-0.9199256,-0.72];
ybs[2321]=['',1.6931247,0.0144368,6.71];
ybs[2322]=['',1.6918313,-0.1313347,6.27];
ybs[2323]=['',1.6852487,-0.6122046,6.25];
ybs[2324]=['16 Gem',1.6980451,0.3574771,6.22];
ybs[2325]=['6 Lyn',1.7128894,1.0148576,5.88];
ybs[2326]=['48 Aur',1.7011968,0.5319509,5.55];
ybs[2327]=['',1.6947742,0.0505134,5.55];
ybs[2328]=['',1.6942041,0.0049807,5.2];
ybs[2329]=['',1.6943152,-0.0050599,5.55];
ybs[2330]=['',1.6759087,-1.0219913,6.48];
ybs[2331]=['',1.6718238,-1.1116793,6.27];
ybs[2332]=['47 Aur',1.7085917,0.8145497,5.9];
ybs[2333]=['',1.7026746,0.4704201,6.47];
ybs[2334]=['',1.7001653,0.2831605,6.23];
ybs[2335]=['',1.6809959,-0.9218625,6.51];
ybs[2336]=['',1.6992799,0.1795865,6.15];
ybs[2337]=['ν Gem',1.7024841,0.3525135,4.15];
ybs[2338]=['10 Mon',1.6972143,-0.0833633,5.06];
ybs[2339]=['',1.6776261,-1.0523152,5.8];
ybs[2340]=['',1.7617952,1.3889092,6.54];
ybs[2341]=['',1.6988413,0.0331248,6.48];
ybs[2342]=['',1.6854223,-0.8410767,5.76];
ybs[2343]=['',1.6930463,-0.4515233,6.07];
ybs[2344]=['',1.783629,1.4327794,6.65];
ybs[2345]=['',1.7023015,0.1920693,6.59];
ybs[2346]=['π1 Dor',1.6687085,-1.2216497,5.56];
ybs[2347]=['',1.6921963,-0.6616404,6.48];
ybs[2348]=['',1.678052,-1.1072556,6.46];
ybs[2349]=['',1.7030918,0.0459253,6.16];
ybs[2350]=['β Mon',1.7008733,-0.1229992,4.6];
ybs[2351]=['β Mon',1.7009095,-0.1230284,5.4];
ybs[2352]=['β Mon',1.7009095,-0.1230284,5.6];
ybs[2353]=['',1.6996443,-0.3050931,5.77];
ybs[2354]=['',1.6801137,-1.1142217,6.27];
ybs[2355]=['λ CMa',1.6970249,-0.5688758,4.48];
ybs[2356]=['',1.7069939,0.1573231,6.57];
ybs[2357]=['',1.7612105,1.3609206,5.73];
ybs[2358]=['',1.6991516,-0.5652342,5.74];
ybs[2359]=['',1.7475224,1.285891,6.24];
ybs[2360]=['',1.7119788,0.2953593,6.2];
ybs[2361]=['',1.706751,-0.1762192,5.93];
ybs[2362]=['',1.6989073,-0.717136,6.32];
ybs[2363]=['',1.6903149,-1.0125659,5.82];
ybs[2364]=['',1.7117268,0.1960892,6.14];
ybs[2365]=['19 Gem',1.7139252,0.2772864,6.4];
ybs[2366]=['',1.7182383,0.5661549,5.87];
ybs[2367]=['',1.7083426,-0.229746,6.16];
ybs[2368]=['',1.7139059,0.2055339,6.65];
ybs[2369]=['',1.7145588,0.2012081,5.23];
ybs[2370]=['7 Lyn',1.7289662,0.9657867,6.45];
ybs[2371]=['π2 Dor',1.6811928,-1.216545,5.38];
ybs[2372]=['',1.7171085,0.2034574,6.03];
ybs[2373]=['',1.7118827,-0.2165514,5.15];
ybs[2374]=['',1.708608,-0.4849383,5.93];
ybs[2375]=['',1.7140041,-0.1426649,5.43];
ybs[2376]=['12 Mon',1.7165705,0.0844656,5.84];
ybs[2377]=['',1.7237554,0.5760827,6.42];
ybs[2378]=['',1.7031137,-0.8770995,5.27];
ybs[2379]=['13 Mon',1.719202,0.1276962,4.5];
ybs[2380]=['',1.716484,-0.1027163,5.6];
ybs[2381]=['ξ1 CMa',1.7135199,-0.4090067,4.33];
ybs[2382]=['',1.7101938,-0.6156621,5.84];
ybs[2383]=['',1.7009643,-0.9925252,5.22];
ybs[2384]=['',1.7089184,-0.7143971,6.2];
ybs[2385]=['',1.7224922,0.2467602,5.53];
ybs[2386]=['',1.7180284,-0.1951784,6.24];
ybs[2387]=['',1.7117054,-0.6450013,6.34];
ybs[2388]=['8 Lyn',1.7434469,1.0727136,5.94];
ybs[2389]=['',1.7220844,-0.0215936,5.1];
ybs[2390]=['',1.7581922,1.2518935,5.92];
ybs[2391]=['',1.7166096,-0.5593248,5.69];
ybs[2392]=['49 Aur',1.7300238,0.4887697,5.27];
ybs[2393]=['',1.7150319,-0.6582143,5.24];
ybs[2394]=['',1.7094642,-0.9048096,5.6];
ybs[2395]=['',1.7876864,1.3882525,5.45];
ybs[2396]=['11 Lyn',1.7426301,0.9920175,5.85];
ybs[2397]=['',1.7205512,-0.3654842,6.4];
ybs[2398]=['14 Mon',1.7273694,0.1318593,6.45];
ybs[2399]=['',1.7364129,0.6706743,5.29];
ybs[2400]=['',1.7297232,0.1740189,5.88];
ybs[2401]=['',1.7185444,-0.6744286,6.44];
ybs[2402]=['',1.7021501,-1.144644,6.29];
ybs[2403]=['',1.7292763,0.0152236,5.8];
ybs[2404]=['',1.7077207,-1.0802759,6.15];
ybs[2405]=['',1.7215263,-0.6326676,5.42];
ybs[2406]=['μ Pic',1.7116388,-1.025732,5.7];
ybs[2407]=['',1.7326148,0.0781799,6.55];
ybs[2408]=['ξ2 CMa',1.72751,-0.4011172,4.54];
ybs[2409]=['',1.7250357,-0.5713115,5.62];
ybs[2410]=['',1.7187155,-0.913603,6.19];
ybs[2411]=['',1.7397004,0.4288615,6.44];
ybs[2412]=['',1.7348478,-0.091272,5.52];
ybs[2413]=['51 Aur',1.7456706,0.6871592,5.69];
ybs[2414]=['ψ3 Aur',1.7464053,0.6960881,5.2];
ybs[2415]=['γ Gem',1.7405066,0.2858881,1.93];
ybs[2416]=['',1.7387872,0.1067525,6.06];
ybs[2417]=['ν1 CMa',1.7334541,-0.3259971,5.7];
ybs[2418]=['',1.7283782,-0.6422416,5.59];
ybs[2419]=['53 Aur',1.7439467,0.5055315,5.79];
ybs[2420]=['',1.7398834,0.1890913,6.38];
ybs[2421]=['ψ2 Aur',1.7487937,0.7412246,4.79];
ybs[2422]=['',1.7353908,-0.2328147,5.97];
ybs[2423]=['ν2 CMa',1.734755,-0.336399,3.95];
ybs[2424]=['',1.7398469,0.0468661,6.17];
ybs[2425]=['',1.7305898,-0.6301837,6.35];
ybs[2426]=['',1.7408262,0.0861874,6.15];
ybs[2427]=['',1.7346155,-0.3950274,6.35];
ybs[2428]=['',1.7516618,0.7678353,6.41];
ybs[2429]=['',1.7253845,-0.9249027,4.39];
ybs[2430]=['',1.7467294,0.3841672,6.04];
ybs[2431]=['',1.7393447,-0.2269611,6.12];
ybs[2432]=['54 Aur',1.7490115,0.4929358,6.03];
ybs[2433]=['',1.7487314,0.4290038,6.38];
ybs[2434]=['',1.7425903,-0.0447305,6.14];
ybs[2435]=['',1.7449478,0.0816997,6.57];
ybs[2436]=['',1.7440127,0.027824,6.21];
ybs[2437]=['ν3 CMa',1.7400595,-0.3186361,4.43];
ybs[2438]=['',1.7354177,-0.6661083,6.04];
ybs[2439]=['',1.7344478,-0.7256272,6.34];
ybs[2440]=['',1.7363464,-0.6459322,5.71];
ybs[2441]=['',1.7390412,-0.5647649,5.27];
ybs[2442]=['',1.7431656,-0.2948378,6.03];
ybs[2443]=['',1.7494637,0.226248,5.97];
ybs[2444]=['',1.74627,-0.2472351,4.82];
ybs[2445]=['ν Pup',1.738302,-0.7542438,3.17];
ybs[2446]=['',1.7584372,0.626765,6.46];
ybs[2447]=['25 Gem',1.7568537,0.4917568,6.42];
ybs[2448]=['',1.7524448,0.1108515,6.51];
ybs[2449]=['',1.7473288,-0.4139117,6.05];
ybs[2450]=['15 Mon',1.7545239,0.1723509,4.66];
ybs[2451]=['',1.7564383,0.2858278,6.28];
ybs[2452]=['',1.7559083,0.1916826,6.11];
ybs[2453]=['ψ4 Aur',1.7652899,0.7767196,5.02];
ybs[2454]=['',1.7475011,-0.5321535,5.71];
ybs[2455]=['',1.7546863,0.0082845,5.79];
ybs[2456]=['',1.7417233,-0.8419391,4.93];
ybs[2457]=['',1.7708301,0.9298087,6.27];
ybs[2458]=['',1.7654881,0.6479571,6.19];
ybs[2459]=['',1.7481466,-0.6663467,6.58];
ybs[2460]=['26 Gem',1.7610295,0.3075966,5.21];
ybs[2461]=['',1.7588068,0.1103735,6.37];
ybs[2462]=['',1.7376062,-1.0742842,6.18];
ybs[2463]=['',1.7580622,-0.1603649,5.19];
ybs[2464]=['12 Lyn',1.7804202,1.037046,4.87];
ybs[2465]=['',1.7697079,0.6298458,6.31];
ybs[2466]=['ε Gem',1.7679852,0.4382357,2.98];
ybs[2467]=['',1.7635796,0.0525646,6.19];
ybs[2468]=['',1.7536373,-0.7045948,6.12];
ybs[2469]=['',1.751346,-0.8324362,6.65];
ybs[2470]=['13 Lyn',1.7827126,0.9973786,5.35];
ybs[2471]=['30 Gem',1.767775,0.2304834,4.49];
ybs[2472]=['',1.765952,0.0682485,5.9];
ybs[2473]=['28 Gem',1.7717566,0.5052443,5.44];
ybs[2474]=['',1.7611622,-0.3921851,6.13];
ybs[2475]=['',1.7582809,-0.6705505,6.29];
ybs[2476]=['ψ5 Aur',1.7811694,0.7601612,5.25];
ybs[2477]=['ξ Gem',1.7734422,0.2246739,3.36];
ybs[2478]=['',1.788611,0.9718028,6.33];
ybs[2479]=['',1.7885673,0.9718029,6.28];
ybs[2480]=['ψ6 Aur',1.7855597,0.8511185,5.22];
ybs[2481]=['',1.7631027,-0.6844305,6.3];
ybs[2482]=['32 Gem',1.7761109,0.2211441,6.46];
ybs[2483]=['42 Cam',1.8023987,1.1789047,5.14];
ybs[2484]=['α CMa',1.7717907,-0.292145,-1.46];
ybs[2485]=['10 CMa',1.7682545,-0.5426709,5.2];
ybs[2486]=['',1.7701309,-0.4775882,6.45];
ybs[2487]=['16 Mon',1.7787439,0.1494688,5.93];
ybs[2488]=['',1.772577,-0.4098838,6.05];
ybs[2489]=['',1.77075,-0.5342206,6.54];
ybs[2490]=['',1.7755269,-0.2586419,5.32];
ybs[2491]=['',1.7828058,0.3171195,6.2];
ybs[2492]=['',1.7721776,-0.5552983,5.92];
ybs[2493]=['',1.7728354,-0.5405564,5.8];
ybs[2494]=['',1.7785818,-0.1768111,5.66];
ybs[2495]=['17 Mon',1.7821716,0.1398628,4.77];
ybs[2496]=['11 CMa',1.7793078,-0.2521867,5.29];
ybs[2497]=['',1.7481593,-1.2530726,6.51];
ybs[2498]=['18 Mon',1.7842911,0.0416837,4.47];
ybs[2499]=['',1.7747222,-0.6905037,6.62];
ybs[2500]=['',1.782846,-0.1574656,5.07];
ybs[2501]=['12 CMa',1.7798167,-0.3672002,6.08];
ybs[2502]=['',1.7754609,-0.6597096,6.21];
ybs[2503]=['43 Cam',1.8147145,1.2018564,5.12];
ybs[2504]=['',1.79343,0.5686595,5.71];
ybs[2505]=['',1.7711101,-0.9114756,6.57];
ybs[2506]=['',1.7861568,-0.0234451,5.75];
ybs[2507]=['',1.7730991,-0.9151253,5.8];
ybs[2508]=['ψ7 Aur',1.7986192,0.728779,5.02];
ybs[2509]=['',1.7894862,0.0170595,6.15];
ybs[2510]=['',1.7804938,-0.6624102,5.26];
ybs[2511]=['33 Gem',1.7933695,0.2823571,5.85];
ybs[2512]=['14 Lyn',1.8102875,1.0371091,5.33];
ybs[2513]=['',1.7902987,-0.0400822,5.74];
ybs[2514]=['',1.7884966,-0.2647516,5.39];
ybs[2515]=['',1.7775013,-0.8951643,5.4];
ybs[2516]=['',1.7763564,-0.9550078,6.46];
ybs[2517]=['35 Gem',1.7958642,0.2336672,5.65];
ybs[2518]=['',1.7789751,-0.9697659,5.61];
ybs[2519]=['',1.845765,1.3429823,4.55];
ybs[2520]=['',1.791514,-0.4206351,6.33];
ybs[2521]=['36 Gem',1.8010842,0.3793535,5.27];
ybs[2522]=['',1.797157,-0.0098819,5.77];
ybs[2523]=['',1.7592292,-1.2765258,6.37];
ybs[2524]=['',1.809115,0.782132,6.26];
ybs[2525]=['',1.8031205,0.4114734,5.65];
ybs[2526]=['',1.7963482,-0.1407852,6.29];
ybs[2527]=['',1.7945419,-0.2986082,5.79];
ybs[2528]=['',1.7659458,-1.2296905,6.11];
ybs[2529]=['',1.792978,-0.4775019,7.04];
ybs[2530]=['κ CMa',1.7916164,-0.5678153,3.96];
ybs[2531]=['59 Aur',1.8083068,0.6779322,6.12];
ybs[2532]=['θ Gem',1.807018,0.5922727,3.6];
ybs[2533]=['60 Aur',1.809148,0.6704062,6.3];
ybs[2534]=['',1.8082802,0.6241663,6.01];
ybs[2535]=['',1.8008709,0.0526374,6.38];
ybs[2536]=['',1.7952916,-0.4503518,6.33];
ybs[2537]=['',1.7940403,-0.5538137,5.7];
ybs[2538]=['',1.7914111,-0.7936855,6.55];
ybs[2539]=['ψ8 Aur',1.8123212,0.6715685,6.48];
ybs[2540]=['',1.7910962,-0.8140132,5.14];
ybs[2541]=['',1.7960286,-0.6002629,4.99];
ybs[2542]=['α Pic',1.7819931,-1.0814979,3.27];
ybs[2543]=['',1.8061511,0.1458036,5.77];
ybs[2544]=['',1.8037581,-0.0932392,6.3];
ybs[2545]=['τ Pup',1.7909012,-0.8838263,2.93];
ybs[2546]=['',1.7902734,-0.9363161,4.4];
ybs[2547]=['',1.8086491,0.1914587,6.24];
ybs[2548]=['',1.8184573,0.7993395,6.34];
ybs[2549]=['',1.8182923,0.7658923,6.13];
ybs[2550]=['',1.7995637,-0.6327863,5.96];
ybs[2551]=['ζ Men',1.7379593,-1.4108025,5.64];
ybs[2552]=['15 Lyn',1.8284253,1.0191649,4.35];
ybs[2553]=['',1.8280827,1.0041702,6.05];
ybs[2554]=['',1.7902233,-1.0519787,6.11];
ybs[2555]=['',1.7981163,-0.8433097,6.42];
ybs[2556]=['38 Gem',1.814267,0.2295203,4.65];
ybs[2557]=['',1.8073352,-0.3326476,5.64];
ybs[2558]=['',1.8075499,-0.3309124,6.14];
ybs[2559]=['',1.8056579,-0.4709572,6.4];
ybs[2560]=['ψ9 Aur',1.8240926,0.8071438,5.87];
ybs[2561]=['37 Gem',1.8176324,0.4424056,5.73];
ybs[2562]=['',1.8114195,-0.1026158,6.41];
ybs[2563]=['15 CMa',1.8083154,-0.3534432,4.83];
ybs[2564]=['',1.8127506,-0.0201418,5.45];
ybs[2565]=['',1.8258458,0.8146598,5.86];
ybs[2566]=['θ CMa',1.811416,-0.2105841,4.07];
ybs[2567]=['',1.8033828,-0.7422989,6.52];
ybs[2568]=['',1.8080481,-0.4985768,6.04];
ybs[2569]=['',1.8140021,-0.0311301,6.21];
ybs[2570]=['',1.8097739,-0.4287572,6.21];
ybs[2571]=['',1.8038273,-0.7679805,6.46];
ybs[2572]=['ο1 CMa',1.8107049,-0.4225582,3.87];
ybs[2573]=['',1.8486134,1.2352971,5.68];
ybs[2574]=['',1.8151811,-0.0494098,6.04];
ybs[2575]=['',1.8110862,-0.4180986,6.91];
ybs[2576]=['',1.818162,0.1448108,6.29];
ybs[2577]=['16 Lyn',1.8287212,0.7865399,4.9];
ybs[2578]=['',1.8254115,0.5873502,5.89];
ybs[2579]=['',1.8030054,-0.9445052,6.57];
ybs[2580]=['17 CMa',1.8148399,-0.3566072,5.74];
ybs[2581]=['',1.8219416,0.1732815,5.92];
ybs[2582]=['π CMa',1.8173742,-0.3519288,4.68];
ybs[2583]=['',1.811176,-0.7398899,6.32];
ybs[2584]=['',1.8023023,-1.0361536,6.41];
ybs[2585]=['μ CMa',1.8197232,-0.2455941,5];
ybs[2586]=['',1.808808,-0.8838077,6.26];
ybs[2587]=['',1.8179485,-0.4008865,5.3];
ybs[2588]=['ι CMa',1.8197301,-0.2981383,4.37];
ybs[2589]=['',1.8263653,0.2073264,6.27];
ybs[2590]=['',1.8181506,-0.5553246,6.36];
ybs[2591]=['',1.8238196,-0.143243,6.34];
ybs[2592]=['62 Aur',1.8345262,0.6636232,6];
ybs[2593]=['39 Gem',1.8328424,0.4546906,6.1];
ybs[2594]=['ι Vol',1.7942409,-1.2389862,5.4];
ybs[2595]=['',1.8243785,-0.3880174,6.61];
ybs[2596]=['',1.8216806,-0.6173198,6.29];
ybs[2597]=['40 Gem',1.8357799,0.4517711,6.4];
ybs[2598]=['',1.8315469,0.1325191,6.27];
ybs[2599]=['',1.8256707,-0.4303829,5.46];
ybs[2600]=['',1.818689,-0.8508302,4.95];
ybs[2601]=['',2.0479668,1.5158523,5.07];
ybs[2602]=['',1.8327193,0.0623593,5.97];
ybs[2603]=['',1.8261676,-0.4811197,6.23];
ybs[2604]=['',1.8239929,-0.6202187,6.23];
ybs[2605]=['',1.8345251,0.1271902,6.35];
ybs[2606]=['',1.8280166,-0.4746171,6.37];
ybs[2607]=['41 Gem',1.8388826,0.2801069,5.68];
ybs[2608]=['',1.8301484,-0.4440633,5.59];
ybs[2609]=['',1.8682151,1.2339311,6.5];
ybs[2610]=['ε CMa',1.8301134,-0.5061681,1.5];
ybs[2611]=['',1.8289731,-0.5958663,5.06];
ybs[2612]=['',1.8440125,0.5652069,6.59];
ybs[2613]=['',1.8304809,-0.5415215,6.42];
ybs[2614]=['',1.838308,-0.0941932,6.3];
ybs[2615]=['',1.8349124,-0.3775656,6.26];
ybs[2616]=['',1.8386191,-0.1472519,5.96];
ybs[2617]=['',1.8370628,-0.3523594,6.31];
ybs[2618]=['',1.8295137,-0.7993104,6.22];
ybs[2619]=['',1.8397263,-0.1611488,6.49];
ybs[2620]=['',1.8378018,-0.386579,6.53];
ybs[2621]=['',1.8447154,0.0835566,6.63];
ybs[2622]=['ω Gem',1.8485719,0.4220953,5.18];
ybs[2623]=['',1.8483745,0.309352,5.94];
ybs[2624]=['',1.8476974,0.267126,5.74];
ybs[2625]=['',1.8457294,0.0964605,6.59];
ybs[2626]=['',1.8285044,-0.9731683,6.27];
ybs[2627]=['',1.8489165,0.2904772,5.82];
ybs[2628]=['',1.8453449,-0.0240201,6.17];
ybs[2629]=['',1.8392346,-0.4977595,6.27];
ybs[2630]=['',1.828185,-0.984779,6.45];
ybs[2631]=['',1.8454523,-0.1004076,5.2];
ybs[2632]=['',1.8410675,-0.4406179,5.63];
ybs[2633]=['',1.839531,-0.584605,6.4];
ybs[2634]=['',1.8667303,1.0431674,6.44];
ybs[2635]=['',1.8535561,0.5114808,5.93];
ybs[2636]=['',1.8644005,0.9202321,6.12];
ybs[2637]=['',1.8617548,0.8332661,6.38];
ybs[2638]=['σ CMa',1.8436595,-0.4880863,3.47];
ybs[2639]=['',1.8518794,0.1589461,5.97];
ybs[2640]=['19 Mon',1.8497494,-0.0745316,4.99];
ybs[2641]=['',1.8534034,0.1905919,5.13];
ybs[2642]=['ζ Gem',1.8558235,0.3584641,3.79];
ybs[2643]=['',1.8544496,0.2192618,5.98];
ybs[2644]=['',1.8385167,-0.8976677,5.14];
ybs[2645]=['ο2 CMa',1.8495322,-0.4165146,3.02];
ybs[2646]=['',1.8561516,0.02542,6.57];
ybs[2647]=['',1.8548328,-0.0934685,5.62];
ybs[2648]=['',1.8540949,-0.1772527,6.45];
ybs[2649]=['γ CMa',1.8530427,-0.273404,4.12];
ybs[2650]=['',1.8452261,-0.7580828,6.43];
ybs[2651]=['44 Gem',1.8611308,0.394529,6.02];
ybs[2652]=['',1.8655127,0.6011099,5.55];
ybs[2653]=['',1.8387225,-1.0292231,6.02];
ybs[2654]=['',1.8317426,-1.1858739,5.17];
ybs[2655]=['',1.8621485,0.1597555,5.78];
ybs[2656]=['',1.8572843,-0.385094,6.09];
ybs[2657]=['',1.8706296,0.592994,5.91];
ybs[2658]=['',1.8530923,-0.7394761,5.2];
ybs[2659]=['',1.8526128,-0.7616554,5.54];
ybs[2660]=['',1.8527217,-0.7617187,6.79];
ybs[2661]=['',1.8705611,0.4912024,6.48];
ybs[2662]=['',1.8622511,-0.1866398,6.49];
ybs[2663]=['',1.870076,0.3956706,7.68];
ybs[2664]=['',1.851926,-0.8659525,4.93];
ybs[2665]=['',1.8743284,0.5898939,6.28];
ybs[2666]=['',1.848191,-1.0333961,5.5];
ybs[2667]=['',1.8761819,0.6529453,6.16];
ybs[2668]=['',1.868352,0.085121,6.11];
ybs[2669]=['',1.8599966,-0.6075517,6.14];
ybs[2670]=['',1.8659384,-0.1976958,5.39];
ybs[2671]=['',1.8655511,-0.2168889,6.48];
ybs[2672]=['',1.8622716,-0.5356095,6.34];
ybs[2673]=['',1.9038637,1.2527958,6.35];
ybs[2674]=['',1.8715723,0.12981,5.75];
ybs[2675]=['',1.8530438,-0.9910227,5.17];
ybs[2676]=['45 Gem',1.8742291,0.2774552,5.44];
ybs[2677]=['',1.8620144,-0.6704748,6.11];
ybs[2678]=['',1.8662753,-0.4362204,6.08];
ybs[2679]=['',1.8578885,-0.8795146,6.46];
ybs[2680]=['',1.8667729,-0.4658435,6.62];
ybs[2681]=['θ Men',1.8118658,-1.3866264,5.45];
ybs[2682]=['',1.8685243,-0.4166721,5.71];
ybs[2683]=['',1.8665879,-0.7143011,5.79];
ybs[2684]=['',1.8820379,0.3702241,6.43];
ybs[2685]=['δ CMa',1.8728507,-0.4612397,1.84];
ybs[2686]=['',1.8775788,-0.1811906,6.21];
ybs[2687]=['',1.8748228,-0.4202427,6.65];
ybs[2688]=['63 Aur',1.8895846,0.6856547,4.9];
ybs[2689]=['τ Gem',1.8868917,0.5272659,4.41];
ybs[2690]=['',1.8662521,-0.9075868,5.96];
ybs[2691]=['',1.8783107,-0.2839485,6.03];
ybs[2692]=['47 Gem',1.8878151,0.4681216,5.78];
ybs[2693]=['20 Mon',1.8816685,-0.0745586,4.92];
ybs[2694]=['',1.8742254,-0.6927172,4.83];
ybs[2695]=['',1.8979716,0.8969698,5.47];
ybs[2696]=['',1.8786741,-0.4409661,5.69];
ybs[2697]=['',1.8808472,-0.3267236,6.23];
ybs[2698]=['48 Gem',1.8923065,0.4204946,5.85];
ybs[2699]=['21 Mon',1.8868824,-0.0058848,5.45];
ybs[2700]=['',1.8812359,-0.4804204,5.46];
ybs[2701]=['',1.9598517,1.4174706,6.31];
ybs[2702]=['',1.8890942,0.0980746,6.09];
ybs[2703]=['',1.8940753,0.4745335,6.43];
ybs[2704]=['',1.859449,-1.2020036,6.47];
ybs[2705]=['',1.8902589,0.0949308,6.16];
ybs[2706]=['δ Mon',1.8889341,-0.0092194,4.15];
ybs[2707]=['18 Lyn',1.9098311,1.0402156,5.2];
ybs[2708]=['',1.8874786,-0.3650948,5.84];
ybs[2709]=['51 Gem',1.8960718,0.2813941,5];
ybs[2710]=['26 CMa',1.8895057,-0.4534027,5.92];
ybs[2711]=['',1.8820893,-0.8546362,5.14];
ybs[2712]=['',1.8887079,-0.538559,6.1];
ybs[2713]=['',1.9083069,0.8238401,5.58];
ybs[2714]=['',1.9010682,0.4306446,6.89];
ybs[2715]=['',1.8940503,-0.1970026,5.78];
ybs[2716]=['',1.8903079,-0.480137,6.59];
ybs[2717]=['52 Gem',1.9021875,0.4336823,5.82];
ybs[2718]=['',1.8900055,-0.638443,5.96];
ybs[2719]=['',1.8890715,-0.7074595,5.31];
ybs[2720]=['',1.9010417,0.21082,5.62];
ybs[2721]=['',1.8998231,0.0536648,5.35];
ybs[2722]=['',1.8948514,-0.396355,6.01];
ybs[2723]=['',1.8989373,-0.0687299,5.75];
ybs[2724]=['',1.8990623,-0.1742548,5.9];
ybs[2725]=['',1.8966172,-0.4004259,6.36];
ybs[2726]=['',1.8955651,-0.4780913,6.12];
ybs[2727]=['γ1 Vol',1.8697766,-1.2309961,5.69];
ybs[2728]=['γ2 Vol',1.8699727,-1.2310255,3.78];
ybs[2729]=['',1.916219,0.9091868,5.92];
ybs[2730]=['53 Gem',1.9077779,0.4862498,5.71];
ybs[2731]=['',1.8999734,-0.1806997,6.03];
ybs[2732]=['',1.8899702,-0.8167289,4.49];
ybs[2733]=['',1.8961892,-0.5431497,6.6];
ybs[2734]=['',1.986455,1.4375624,4.96];
ybs[2735]=['',1.8969572,-0.5301678,6.33];
ybs[2736]=['24 Mon',1.9040276,-0.003464,6.41];
ybs[2737]=['27 CMa',1.8984365,-0.4605755,4.66];
ybs[2738]=['',1.8929742,-0.7892213,4.89];
ybs[2739]=['',1.9057585,0.1385882,5.82];
ybs[2740]=['',1.8943926,-0.7797409,5.1];
ybs[2741]=['ω CMa',1.9008561,-0.4679153,3.85];
ybs[2742]=['',1.9010197,-0.4725456,5.58];
ybs[2743]=['',1.9202296,0.8626515,5.05];
ybs[2744]=['',1.9054119,-0.1853739,5.95];
ybs[2745]=['64 Aur',1.9175025,0.7128778,5.78];
ybs[2746]=['',1.8859345,-1.1034902,6.02];
ybs[2747]=['',1.9052555,-0.4150012,6.32];
ybs[2748]=['',1.9030382,-0.5362249,5.36];
ybs[2749]=['',1.9171423,0.5396104,6.24];
ybs[2750]=['',1.9075224,-0.2726784,5.46];
ybs[2751]=['',1.9007543,-0.7236598,5.94];
ybs[2752]=['',1.9128705,0.115934,6.65];
ybs[2753]=['',1.8995914,-0.8183227,5.72];
ybs[2754]=['',1.8989281,-0.8431391,4.76];
ybs[2755]=['λ Gem',1.9166836,0.2880118,3.58];
ybs[2756]=['',1.9088584,-0.4075903,4.79];
ybs[2757]=['',1.9134414,-0.1172532,6.29];
ybs[2758]=['',1.9085429,-0.4872739,4.64];
ybs[2759]=['',1.9016995,-0.9169382,5.97];
ybs[2760]=['',1.9100275,-0.5399082,6.32];
ybs[2761]=['',1.9078184,-0.6694466,5.8];
ybs[2762]=['',1.9091962,-0.6393228,5.03];
ybs[2763]=['',1.9061171,-0.8170212,5.66];
ybs[2764]=['47 Cam',1.9376341,1.044779,6.35];
ybs[2765]=['π Pup',1.91056,-0.6481344,2.7];
ybs[2766]=['',1.9138912,-0.4683713,6.46];
ybs[2767]=['',1.9307538,0.7437842,6.35];
ybs[2768]=['',1.9319618,0.7886806,5.77];
ybs[2769]=['δ Gem',1.9257387,0.3829749,3.53];
ybs[2770]=['',1.9218012,0.0471512,5.89];
ybs[2771]=['',1.9237786,0.1239809,5.91];
ybs[2772]=['',1.9254537,0.2636043,6.45];
ybs[2773]=['29 CMa',1.9177933,-0.4293074,4.98];
ybs[2774]=['τ CMa',1.9179304,-0.4362066,4.4];
ybs[2775]=['19 Lyn',1.9395583,0.9641843,6.53];
ybs[2776]=['19 Lyn',1.9396451,0.9641309,5.45];
ybs[2777]=['',1.9195727,-0.3371814,6.09];
ybs[2778]=['',1.9185047,-0.4646856,5.28];
ybs[2779]=['',1.9156572,-0.6418026,4.66];
ybs[2780]=['',1.9215753,-0.2868274,5.7];
ybs[2781]=['',1.9142129,-0.7683803,5.85];
ybs[2782]=['',1.9170974,-0.6419555,5.11];
ybs[2783]=['',1.9166326,-0.6850208,5.25];
ybs[2784]=['',1.935655,0.6799054,6.4];
ybs[2785]=['65 Aur',1.9347534,0.6408892,5.13];
ybs[2786]=['',1.9198282,-0.5893292,6.3];
ybs[2787]=['56 Gem',1.9336347,0.3561062,5.1];
ybs[2788]=['',1.9282003,-0.2513223,5.45];
ybs[2789]=['',1.9997368,1.4110992,6.41];
ybs[2790]=['',1.9297467,-0.1556518,6.55];
ybs[2791]=['',1.9275275,-0.3995289,6.61];
ybs[2792]=['',1.9274953,-0.471296,6.01];
ybs[2793]=['',1.9334329,0.0023908,5.99];
ybs[2794]=['',1.9282155,-0.4525835,5.87];
ybs[2795]=['δ Vol',1.9059646,-1.1867326,3.98];
ybs[2796]=['',1.9484214,0.9048749,5.8];
ybs[2797]=['66 Aur',1.9441126,0.7091435,5.19];
ybs[2798]=['',1.9330235,-0.1574179,6.43];
ybs[2799]=['',1.934421,-0.0526957,6.23];
ybs[2800]=['57 Gem',1.9404773,0.4365,5.03];
ybs[2801]=['',1.9610723,1.1569563,6.47];
ybs[2802]=['58 Gem',1.9403737,0.3997561,6.02];
ybs[2803]=['',1.9348245,-0.1051243,5.82];
ybs[2804]=['',1.933512,-0.3326065,4.96];
ybs[2805]=['',1.9235806,-0.913697,6.05];
ybs[2806]=['',1.9236026,-0.9136632,6.6];
ybs[2807]=['',1.9248594,-0.9097627,5.39];
ybs[2808]=['59 Gem',1.9452998,0.4816516,5.76];
ybs[2809]=['',1.944429,0.2701043,6.41];
ybs[2810]=['21 Lyn',1.9558799,0.8581588,4.64];
ybs[2811]=['',1.9364196,-0.557886,5.43];
ybs[2812]=['1 CMi',1.9465147,0.2029488,5.3];
ybs[2813]=['ι Gem',1.9504095,0.4844346,3.79];
ybs[2814]=['',1.9386699,-0.4865109,5.38];
ybs[2815]=['',1.9386842,-0.562748,5.39];
ybs[2816]=['',1.940402,-0.5281015,6.6];
ybs[2817]=['',1.944276,-0.2834857,5.33];
ybs[2818]=['',1.9423649,-0.4006231,6.19];
ybs[2819]=['η CMa',1.9412707,-0.5121527,2.45];
ybs[2820]=['ε CMi',1.949394,0.161167,4.99];
ybs[2821]=['',1.9404329,-0.6262038,6.31];
ybs[2822]=['',1.9765542,1.1941716,5.64];
ybs[2823]=['',1.9449468,-0.3325503,6.24];
ybs[2824]=['',1.9464122,-0.2407436,5.78];
ybs[2825]=['',1.949784,-0.1015257,5.97];
ybs[2826]=['',1.943933,-0.5558928,5.35];
ybs[2827]=['',1.955006,0.3751295,6.54];
ybs[2828]=['',1.9530013,0.1844119,6.37];
ybs[2829]=['61 Gem',1.9554022,0.3528127,5.93];
ybs[2830]=['',1.9507334,-0.0799289,6.76];
ybs[2831]=['',1.9469612,-0.3844001,6.05];
ybs[2832]=['',1.953989,0.191406,6.41];
ybs[2833]=['',1.9472293,-0.4408621,5.78];
ybs[2834]=['',1.943918,-0.6515565,6.97];
ybs[2835]=['',1.9439251,-0.651571,6.84];
ybs[2836]=['',1.9651626,0.8402083,5.72];
ybs[2837]=['β CMi',1.9559129,0.1439343,2.9];
ybs[2838]=['63 Gem',1.958943,0.3735369,5.22];
ybs[2839]=['',1.9482432,-0.5546742,6.31];
ybs[2840]=['',1.7425739,-1.5171984,6.47];
ybs[2841]=['22 Lyn',1.9699524,0.8661808,5.36];
ybs[2842]=['',1.9527731,-0.4145954,6.56];
ybs[2843]=['η CMi',1.9597299,0.1204089,5.25];
ybs[2844]=['ρ Gem',1.9653436,0.5539828,4.18];
ybs[2845]=['',1.9549752,-0.3125363,5.63];
ybs[2846]=['γ CMi',1.9603548,0.1550283,4.32];
ybs[2847]=['',1.9541574,-0.4036703,5.61];
ybs[2848]=['',1.9524474,-0.5966088,5.9];
ybs[2849]=['64 Gem',1.9661826,0.4899905,5.05];
ybs[2850]=['',1.9632896,0.2629524,6.22];
ybs[2851]=['',1.9583715,-0.2024557,5.79];
ybs[2852]=['',1.9573304,-0.3997248,5.95];
ybs[2853]=['65 Gem',1.9682314,0.4864622,5.01];
ybs[2854]=['',1.9499297,-0.891173,5.1];
ybs[2855]=['',1.9582645,-0.5096147,5.54];
ybs[2856]=['6 CMi',1.967583,0.2087907,4.54];
ybs[2857]=['',1.9650145,-0.0340143,5.59];
ybs[2858]=['',1.9653407,-0.1325535,5.86];
ybs[2859]=['',1.9649875,-0.1809954,5.75];
ybs[2860]=['',1.9648152,-0.2625457,6.05];
ybs[2861]=['',1.9595789,-0.660666,6.58];
ybs[2862]=['',1.9619381,-0.5566145,6.38];
ybs[2863]=['',1.9619527,-0.5565903,7.13];
ybs[2864]=['',1.9779702,0.6780868,6.54];
ybs[2865]=['',1.962945,-0.5497757,5.77];
ybs[2866]=['',1.9666719,-0.402617,4.85];
ybs[2867]=['',1.9626481,-0.678159,5.43];
ybs[2868]=['',1.9716273,-0.0919909,6.24];
ybs[2869]=['',1.9765253,0.2974274,5.42];
ybs[2870]=['σ Pup',1.962983,-0.7565105,3.25];
ybs[2871]=['',1.9812574,0.3986771,6.54];
ybs[2872]=['δ1 CMi',1.977296,0.03263,5.25];
ybs[2873]=['',1.9700677,-0.5411639,4.65];
ybs[2874]=['',1.9697609,-0.652472,6.65];
ybs[2875]=['',1.9769413,-0.1557828,5.9];
ybs[2876]=['',1.9656384,-0.9196994,5.87];
ybs[2877]=['',1.9729782,-0.6317617,6.68];
ybs[2878]=['68 Gem',1.984341,0.2754316,5.25];
ybs[2879]=['δ2 CMi',1.9821178,0.0566341,5.59];
ybs[2880]=['',1.9592127,-1.1266656,6.39];
ybs[2881]=['',1.9742281,-0.6271389,6.61];
ybs[2882]=['α Gem',1.9892839,0.555757,2.88];
ybs[2883]=['α Gem',1.9892839,0.5557521,1.98];
ybs[2884]=['',1.9677579,-0.9502177,5.96];
ybs[2885]=['',1.9862437,0.1836527,6.28];
ybs[2886]=['',2.0003401,0.9722902,5.92];
ybs[2887]=['',1.9770987,-0.6284289,6.3];
ybs[2888]=['',1.9916203,0.5395647,5.33];
ybs[2889]=['',1.9823091,-0.2510441,6.21];
ybs[2890]=['',1.9956642,0.7502192,6.3];
ybs[2891]=['',1.9819431,-0.3396045,5.66];
ybs[2892]=['',1.9810331,-0.4320764,5.85];
ybs[2893]=['δ3 CMi',1.9867889,0.0580411,5.81];
ybs[2894]=['',1.9841795,-0.2542861,4.97];
ybs[2895]=['',1.9984479,0.8051776,5.65];
ybs[2896]=['',1.9889643,0.0467556,6.55];
ybs[2897]=['υ Gem',1.9948391,0.4686062,4.06];
ybs[2898]=['',1.9850166,-0.3899387,4.45];
ybs[2899]=['',1.9806131,-0.6999504,6.26];
ybs[2900]=['',1.9804256,-0.75279,6.52];
ybs[2901]=['',1.9860927,-0.410492,5.83];
ybs[2902]=['',1.9861291,-0.4105115,5.87];
ybs[2903]=['',1.9835343,-0.6350196,5.54];
ybs[2904]=['',1.9867326,-0.4566233,6.65];
ybs[2905]=['',1.9852569,-0.5848444,6.11];
ybs[2906]=['',2.004607,0.8504288,5.92];
ybs[2907]=['',2.0014386,0.6977469,6.38];
ybs[2908]=['',1.9871408,-0.4722496,5.77];
ybs[2909]=['',1.9838961,-0.6972851,6.76];
ybs[2910]=['',1.9969689,0.1014863,5.91];
ybs[2911]=['ε Men',1.9393733,-1.381176,5.53];
ybs[2912]=['',1.9952007,-0.1458774,6.27];
ybs[2913]=['',1.9940728,-0.2537611,5.7];
ybs[2914]=['',1.9905705,-0.4959487,4.64];
ybs[2915]=['',1.9940858,-0.3875894,6.34];
ybs[2916]=['70 Gem',2.006639,0.610878,5.56];
ybs[2917]=['',1.9861134,-0.899205,6.28];
ybs[2918]=['',2.0048593,0.4243343,6.27];
ybs[2919]=['25 Mon',1.9997065,-0.0725768,5.13];
ybs[2920]=['',1.9965904,-0.3446878,5.74];
ybs[2921]=['23 Lyn',2.0180894,0.9954275,6.06];
ybs[2922]=['ο Gem',2.0093137,0.6027671,4.9];
ybs[2923]=['',2.0090267,0.421922,6.17];
ybs[2924]=['',2.0009846,-0.2528719,6.53];
ybs[2925]=['',1.999045,-0.4157758,6.37];
ybs[2926]=['',1.9904067,-0.9176988,4.94];
ybs[2927]=['',2.0141989,0.6683876,5.73];
ybs[2928]=['',2.0124102,0.5578288,6.17];
ybs[2929]=['',1.9989608,-0.6111414,4.53];
ybs[2930]=['74 Gem',2.010004,0.3076397,5.05];
ybs[2931]=['',2.0189826,0.8391989,5.56];
ybs[2932]=['',1.9953738,-0.8530674,5.72];
ybs[2933]=['',1.9916959,-0.976233,6.39];
ybs[2934]=['',2.0005979,-0.6165307,6.6];
ybs[2935]=['α CMi',2.0088304,0.0903529,0.38];
ybs[2936]=['',2.0034365,-0.4435297,4.7];
ybs[2937]=['',2.0004926,-0.6642363,6.38];
ybs[2938]=['24 Lyn',2.0277758,1.0238155,4.99];
ybs[2939]=['',2.0072543,-0.3268515,6.72];
ybs[2940]=['',2.0056573,-0.4686132,4.5];
ybs[2941]=['',2.0056936,-0.4686472,4.62];
ybs[2942]=['',2.0123864,0.0904483,6.02];
ybs[2943]=['',2.0167309,0.4008963,5.89];
ybs[2944]=['',2.0032269,-0.6988134,6.59];
ybs[2945]=['',2.0155932,0.2394938,6.24];
ybs[2946]=['',2.0048396,-0.6378266,5.8];
ybs[2947]=['',2.0039022,-0.6776913,6.19];
ybs[2948]=['',2.0083903,-0.4696896,6.5];
ybs[2949]=['',2.0022546,-0.8490801,5.68];
ybs[2950]=['',2.0140348,-0.1437201,6.01];
ybs[2951]=['',2.0129076,-0.2672489,4.94];
ybs[2952]=['',2.0120572,-0.3439935,5.93];
ybs[2953]=['',2.0079144,-0.6694469,4.84];
ybs[2954]=['',2.0248112,0.5925484,6.02];
ybs[2955]=['',2.0091164,-0.6665015,5.73];
ybs[2956]=['',2.0094013,-0.6686206,5.76];
ybs[2957]=['',2.0202732,0.2344191,5.77];
ybs[2958]=['',2.0187493,0.0624048,5.94];
ybs[2959]=['',2.0211262,0.2471197,5.56];
ybs[2960]=['',2.0101781,-0.6567294,6];
ybs[2961]=['',2.0316328,0.8793576,5.27];
ybs[2962]=['α Mon',2.0168492,-0.1675538,3.93];
ybs[2963]=['',2.0049534,-0.930631,6.06];
ybs[2964]=['',2.0139128,-0.4885976,6.76];
ybs[2965]=['σ Gem',2.0271497,0.5032413,4.28];
ybs[2966]=['',2.0144341,-0.5534377,6.56];
ybs[2967]=['51 Cam',2.0446942,1.1415183,5.92];
ybs[2968]=['',2.0170528,-0.3907143,6.18];
ybs[2969]=['49 Cam',2.0433445,1.0957006,6.49];
ybs[2970]=['',2.0271545,0.3900711,6.21];
ybs[2971]=['',1.9849906,-1.2971564,7.16];
ybs[2972]=['',1.9849979,-1.2971564,7.26];
ybs[2973]=['',2.0157912,-0.673393,5.42];
ybs[2974]=['',2.0252009,0.0024362,6.19];
ybs[2975]=['76 Gem',2.0305307,0.4491396,5.31];
ybs[2976]=['',2.0158841,-0.7798343,6.41];
ybs[2977]=['κ Gem',2.0319246,0.4249449,3.57];
ybs[2978]=['',2.0188624,-0.6733161,6.54];
ybs[2979]=['',2.0305916,0.2235601,6.43];
ybs[2980]=['',2.0230484,-0.4607808,5.64];
ybs[2981]=['',2.0297788,0.0410969,6.47];
ybs[2982]=['β Gem',2.0358521,0.4882595,1.14];
ybs[2983]=['79 Gem',2.0348655,0.353701,6.33];
ybs[2984]=['',2.0267915,-0.4460006,6.55];
ybs[2985]=['1 Pup',2.0261897,-0.4967402,4.59];
ybs[2986]=['',2.0243702,-0.6300659,5.6];
ybs[2987]=['',2.0238665,-0.6791768,6.89];
ybs[2988]=['3 Pup',2.0273389,-0.5062301,3.96];
ybs[2989]=['',2.0928614,1.3999176,6.56];
ybs[2990]=['',2.0227939,-0.7892909,5.06];
ybs[2991]=['',2.0421145,0.6539047,5.18];
ybs[2992]=['',1.9862432,-1.3557792,6.18];
ybs[2993]=['',2.0265129,-0.6676235,6.4];
ybs[2994]=['',2.0262962,-0.7153047,5.17];
ybs[2995]=['81 Gem',2.0390261,0.322166,4.88];
ybs[2996]=['',2.0308173,-0.4315216,5.62];
ybs[2997]=['',2.0231671,-0.8734071,6.57];
ybs[2998]=['',2.0182319,-1.0241615,6.43];
ybs[2999]=['',2.0285611,-0.6302916,5.8];
ybs[3000]=['11 CMi',2.0394087,0.1870476,5.3];
ybs[3001]=['2 Pup',2.0351487,-0.2572094,6.89];
ybs[3002]=['2 Pup',2.0351777,-0.2572919,6.07];
ybs[3003]=['',2.0302591,-0.6631117,5.88];
ybs[3004]=['',2.0214165,-1.0171712,6.21];
ybs[3005]=['π Gem',2.0456292,0.5823059,5.14];
ybs[3006]=['',2.0378408,-0.1190954,5.49];
ybs[3007]=['4 Pup',2.037204,-0.2550799,5.04];
ybs[3008]=['',2.0324741,-0.662151,6.54];
ybs[3009]=['',2.0332487,-0.6635631,3.61];
ybs[3010]=['',2.0348634,-0.5973209,5.37];
ybs[3011]=['',2.0407575,-0.2221236,6.39];
ybs[3012]=['',2.0331364,-0.7645062,6.03];
ybs[3013]=['82 Gem',2.0498187,0.4029752,6.18];
ybs[3014]=['',2.0372602,-0.6629643,5.88];
ybs[3015]=['',2.042428,-0.3939449,5.9];
ybs[3016]=['ζ Vol',2.0139028,-1.2680708,3.95];
ybs[3017]=['',2.0388246,-0.70007,6.57];
ybs[3018]=['',2.0445468,-0.2799979,6.34];
ybs[3019]=['',2.0450333,-0.2804108,6.43];
ybs[3020]=['',2.0625546,0.9437978,6.02];
ybs[3021]=['5 Pup',2.0459957,-0.2137166,5.48];
ybs[3022]=['',2.0515479,0.2324484,6.04];
ybs[3023]=['',2.0333987,-0.9908817,6.12];
ybs[3024]=['',2.0412207,-0.6873623,6.31];
ybs[3025]=['',2.0510362,0.0747051,6.53];
ybs[3026]=['ο Pup',2.0461346,-0.4535982,4.5];
ybs[3027]=['',2.0426601,-0.6730483,5.08];
ybs[3028]=['',2.0283261,-1.1540524,6.38];
ybs[3029]=['',2.0426673,-0.8143767,5.23];
ybs[3030]=['',2.0252354,-1.2194877,6.18];
ybs[3031]=['',2.069255,0.9626405,6.38];
ybs[3032]=['',2.061029,0.5791031,6.03];
ybs[3033]=['',2.0457141,-0.7104232,6.14];
ybs[3034]=['',2.0526425,-0.233979,6.23];
ybs[3035]=['',2.050284,-0.4357159,5.33];
ybs[3036]=['6 Pup',2.0534305,-0.301612,5.18];
ybs[3037]=['ξ Pup',2.05145,-0.4348017,3.34];
ybs[3038]=['',2.0461885,-0.8225715,4.71];
ybs[3039]=['',2.0558311,-0.1612043,5.61];
ybs[3040]=['',2.0536191,-0.353599,6.56];
ybs[3041]=['',2.0508014,-0.6160291,5.93];
ybs[3042]=['',2.0588934,0.0562683,6.18];
ybs[3043]=['',2.0551342,-0.3416753,6.12];
ybs[3044]=['',2.052396,-0.5819203,5.6];
ybs[3045]=['',2.064448,0.3363504,5.99];
ybs[3046]=['',2.0590132,-0.1951615,6.16];
ybs[3047]=['',2.0501657,-0.8102837,4.11];
ybs[3048]=['',2.045338,-0.9865153,6.33];
ybs[3049]=['',2.0512797,-0.781987,6.32];
ybs[3050]=['',2.050032,-0.8187386,5.84];
ybs[3051]=['ζ CMi',2.0628173,0.0299021,5.14];
ybs[3052]=['',2.058905,-0.429031,6.45];
ybs[3053]=['',2.0646963,0.0562581,6.31];
ybs[3054]=['',2.0488032,-0.9854645,5.59];
ybs[3055]=['8 Pup',2.062283,-0.2246779,6.36];
ybs[3056]=['9 Pup',2.0626408,-0.2435039,5.17];
ybs[3057]=['25 Lyn',2.0768257,0.8260838,6.25];
ybs[3058]=['26 Lyn',2.0778111,0.8291994,5.45];
ybs[3059]=['φ Gem',2.0714851,0.4662005,4.97];
ybs[3060]=['',2.0621475,-0.3704905,5.63];
ybs[3061]=['',2.0566911,-0.7789906,6.45];
ybs[3062]=['',2.0488179,-1.0530625,5.78];
ybs[3063]=['',2.0549415,-0.882486,5.91];
ybs[3064]=['',2.0673889,-0.0956826,5.76];
ybs[3065]=['10 Pup',2.0649812,-0.2600596,5.69];
ybs[3066]=['',2.0595521,-0.753092,6.32];
ybs[3067]=['',2.1057816,1.2891069,5.41];
ybs[3068]=['',2.0518989,-1.04901,6.72];
ybs[3069]=['',2.0860969,0.9852132,6.72];
ybs[3070]=['',2.0610052,-0.7494778,6.04];
ybs[3071]=['',2.0639973,-0.6066615,5.01];
ybs[3072]=['',2.0635127,-0.7091214,3.73];
ybs[3073]=['',2.0499521,-1.156253,5.79];
ybs[3074]=['',2.1289686,1.3861394,5.42];
ybs[3075]=['',2.0813663,0.6171013,6.23];
ybs[3076]=['',2.065472,-0.6792312,4.49];
ybs[3077]=['',2.0673989,-0.6356158,5.43];
ybs[3078]=['85 Gem',2.0806933,0.3460716,5.35];
ybs[3079]=['',2.0797172,0.1537184,5.86];
ybs[3080]=['',2.063781,-0.9498278,5.7];
ybs[3081]=['',2.0666424,-0.8668567,4.63];
ybs[3082]=['',2.0678049,-0.8405041,4.24];
ybs[3083]=['',2.0723338,-0.6271353,5.49];
ybs[3084]=['',2.0744843,-0.6091524,6.15];
ybs[3085]=['',2.0833999,0.0773198,6.17];
ybs[3086]=['',2.0931272,0.7665638,6.34];
ybs[3087]=['1 Cnc',2.0863358,0.2746147,5.78];
ybs[3088]=['',2.0771138,-0.5405751,6.44];
ybs[3089]=['',2.0873086,0.1498412,6.05];
ybs[3090]=['',2.0871009,0.0186895,6.35];
ybs[3091]=['',2.0821677,-0.5295495,6.33];
ybs[3092]=['',2.0749157,-0.9187077,6.38];
ybs[3093]=['',2.07887,-0.7662061,6.02];
ybs[3094]=['11 Pup',2.0845376,-0.4003069,4.2];
ybs[3095]=['',2.0908929,0.1249155,6.41];
ybs[3096]=['',2.0930555,0.2873149,5.99];
ybs[3097]=['',2.0739456,-1.0010859,5.63];
ybs[3098]=['',2.1077053,1.0295608,5.77];
ybs[3099]=['',2.0817943,-0.7119556,6.78];
ybs[3100]=['',2.1883332,1.4659481,6.49];
ybs[3101]=['53 Cam',2.1094302,1.0518449,6.01];
ybs[3102]=['14 CMi',2.0918173,0.0378412,5.29];
ybs[3103]=['',2.0841492,-0.7411019,6.09];
ybs[3104]=['',2.1133109,1.1001115,6.4];
ybs[3105]=['',2.0877903,-0.5304223,4.79];
ybs[3106]=['',2.0840872,-0.7601987,5.35];
ybs[3107]=['',2.0975736,0.2301233,6.02];
ybs[3108]=['',2.0855439,-0.770838,5.09];
ybs[3109]=['χ Car',2.0826154,-0.9256878,3.47];
ybs[3110]=['',2.0854224,-0.8368211,6.22];
ybs[3111]=['',2.1129675,0.9985915,6.49];
ybs[3112]=['',2.0797963,-1.0573542,5.74];
ybs[3113]=['',2.0878866,-0.7964646,5.17];
ybs[3114]=['27 Mon',2.0977115,-0.0652211,4.93];
ybs[3115]=['12 Pup',2.0942855,-0.4078384,5.11];
ybs[3116]=['ω1 Cnc',2.1038593,0.44218,5.83];
ybs[3117]=['',2.1030887,0.3448499,6.25];
ybs[3118]=['',2.0823033,-1.0329238,6.25];
ybs[3119]=['',2.1041427,0.4105938,6.34];
ybs[3120]=['3 Cnc',2.1029541,0.301086,5.55];
ybs[3121]=['',2.0892891,-0.8604723,4.41];
ybs[3122]=['',2.108567,0.6170591,6.34];
ybs[3123]=['',2.0978234,-0.3221244,4.61];
ybs[3124]=['ω2 Cnc',2.10733,0.4368847,6.31];
ybs[3125]=['',2.0896348,-0.8989334,6.44];
ybs[3126]=['5 Cnc',2.1060522,0.2861872,5.99];
ybs[3127]=['',2.1020991,-0.0513,6.51];
ybs[3128]=['',2.1044922,0.0841579,5.65];
ybs[3129]=['',2.0930045,-0.7901611,5.99];
ybs[3130]=['',2.0862566,-1.0534723,5.6];
ybs[3131]=['',2.0833422,-1.105716,6.14];
ybs[3132]=['',2.0952677,-0.6868606,5.24];
ybs[3133]=['28 Mon',2.1042674,-0.0253128,4.68];
ybs[3134]=['',2.0934519,-0.8732496,6.32];
ybs[3135]=['',2.0935322,-0.8731964,6.34];
ybs[3136]=['',2.1072979,0.1545627,6.22];
ybs[3137]=['',2.1089314,0.0397267,4.39];
ybs[3138]=['',2.0986533,-0.7943742,6.61];
ybs[3139]=['',2.0908243,-1.0625753,5.81];
ybs[3140]=['',2.0980807,-0.8558866,6.02];
ybs[3141]=['χ Gem',2.1152297,0.4840727,4.94];
ybs[3142]=['',2.1094059,-0.1116234,6.33];
ybs[3143]=['',2.0991077,-0.8539684,6.12];
ybs[3144]=['',2.094522,-1.0518139,6.33];
ybs[3145]=['',2.0942856,-1.0584362,5.17];
ybs[3146]=['',2.1047546,-0.6517328,5.95];
ybs[3147]=['',2.1068616,-0.6476688,6.34];
ybs[3148]=['',2.1002248,-0.9461241,5.87];
ybs[3149]=['',2.1025959,-0.9524744,6.1];
ybs[3150]=['',2.1202994,0.3278229,6.15];
ybs[3151]=['',2.0969634,-1.1104616,4.82];
ybs[3152]=['',2.1112787,-0.5676237,5.82];
ybs[3153]=['',2.1031364,-0.9688815,6.28];
ybs[3154]=['',2.1094585,-0.7220147,5.52];
ybs[3155]=['8 Cnc',2.1215187,0.2279151,5.12];
ybs[3156]=['',2.1243704,0.4794417,6.21];
ybs[3157]=['ζ Pup',2.1131875,-0.6992154,2.25];
ybs[3158]=['',2.1126232,-0.7506194,6.29];
ybs[3159]=['28 Lyn',2.1318313,0.7539797,6.26];
ybs[3160]=['14 Pup',2.1188394,-0.345354,6.13];
ybs[3161]=['μ1 Cnc',2.127203,0.3940175,5.99];
ybs[3162]=['',2.1165165,-0.5713172,5.31];
ybs[3163]=['',2.0900344,-1.279351,6.34];
ybs[3164]=['',2.1243853,-0.0110548,6.41];
ybs[3165]=['27 Lyn',2.1379935,0.8978963,4.84];
ybs[3166]=['',2.1268827,-0.1624035,6.23];
ybs[3167]=['',2.1455937,1.0155436,5.93];
ybs[3168]=['μ2 Cnc',2.1335199,0.3756128,5.3];
ybs[3169]=['',2.1229334,-0.5869344,6.14];
ybs[3170]=['',2.1174415,-0.8840051,5.95];
ybs[3171]=['',2.1204561,-0.8209745,6.19];
ybs[3172]=['',2.1188178,-0.9279462,5.53];
ybs[3173]=['',2.141463,0.739482,6.27];
ybs[3174]=['',2.1590949,1.1940016,5.32];
ybs[3175]=['',2.130197,-0.3598013,5.38];
ybs[3176]=['12 Cnc',2.1373771,0.2370125,6.27];
ybs[3177]=['ρ Pup',2.1311296,-0.4252431,2.81];
ybs[3178]=['',2.1162285,-1.0977293,6.3];
ybs[3179]=['',2.126386,-0.7910957,5.05];
ybs[3180]=['ζ Mon',2.1363864,-0.0531424,4.34];
ybs[3181]=['',2.1376857,-0.1989816,6.32];
ybs[3182]=['',2.1364236,-0.3564666,6.36];
ybs[3183]=['ψ Cnc',2.1453889,0.4441069,5.73];
ybs[3184]=['16 Pup',2.1377775,-0.3369551,4.4];
ybs[3185]=['',2.1498852,0.6749051,6.58];
ybs[3186]=['',2.1398314,-0.2846665,5.68];
ybs[3187]=['',2.1353179,-0.6587272,6.37];
ybs[3188]=['',2.1377516,-0.5302943,6.65];
ybs[3189]=['',2.2180701,1.4375014,6.32];
ybs[3190]=['',2.1473245,0.2542504,6.23];
ybs[3191]=['',2.1377921,-0.6198735,6.2];
ybs[3192]=['',2.1618508,0.984173,5.85];
ybs[3193]=['',2.1484716,0.170327,6.07];
ybs[3194]=['18 Pup',2.1450949,-0.2419195,5.54];
ybs[3195]=['',2.1370314,-0.8507702,5.7];
ybs[3196]=['',2.1392288,-0.7711575,5.21];
ybs[3197]=['',2.1401703,-0.7452894,6.26];
ybs[3198]=['γ1 Vel',2.1385538,-0.8274095,4.27];
ybs[3199]=['γ2 Vel',2.1387509,-0.8272499,1.78];
ybs[3200]=['ζ1 Cnc',2.1527923,0.3069211,5.63];
ybs[3201]=['ζ1 Cnc',2.1527923,0.3069211,6.02];
ybs[3202]=['ζ2 Cnc',2.1528359,0.306921,6.2];
ybs[3203]=['19 Pup',2.1477768,-0.2267008,4.72];
ybs[3204]=['',2.149147,-0.1367409,5.36];
ybs[3205]=['',2.1395275,-0.8377377,5.23];
ybs[3206]=['',2.1533611,0.2433223,6.54];
ybs[3207]=['15 Cnc',2.1572951,0.5165087,5.64];
ybs[3208]=['',2.1906214,1.3210599,5.54];
ybs[3209]=['',2.1322351,-1.114599,6.28];
ybs[3210]=['',2.1381899,-0.9799464,5.66];
ybs[3211]=['',2.1458166,-0.6519526,6.44];
ybs[3212]=['',2.1352206,-1.0709948,4.76];
ybs[3213]=['',2.1709968,1.0527206,6.45];
ybs[3214]=['',2.1561803,0.2871302,6.01];
ybs[3215]=['ε Vol',2.1292593,-1.1986518,4.35];
ybs[3216]=['',2.1594433,0.4027288,6.56];
ybs[3217]=['',2.1471541,-0.6925581,4.45];
ybs[3218]=['',2.1472998,-0.7513519,4.75];
ybs[3219]=['',2.1458929,-0.8469017,5.82];
ybs[3220]=['',2.1614004,0.3073965,6.47];
ybs[3221]=['20 Pup',2.1566862,-0.2766564,4.99];
ybs[3222]=['',2.1537256,-0.523136,6.52];
ybs[3223]=['',2.1619762,0.2266303,6.38];
ybs[3224]=['',2.1495698,-0.8151815,5.76];
ybs[3225]=['',2.1537861,-0.6630003,6.43];
ybs[3226]=['',2.1518428,-0.808553,6.03];
ybs[3227]=['29 Lyn',2.1796067,1.0385789,5.64];
ybs[3228]=['',2.1942843,1.2625892,5.98];
ybs[3229]=['',2.1566511,-0.627667,4.78];
ybs[3230]=['',2.1575883,-0.5869925,6.37];
ybs[3231]=['',2.1598229,-0.5620671,6.06];
ybs[3232]=['',2.1587284,-0.6350492,5.08];
ybs[3233]=['',2.1587567,-0.6353741,6.11];
ybs[3234]=['',2.1598423,-0.6205308,5.78];
ybs[3235]=['',2.158866,-0.7053088,4.44];
ybs[3236]=['',2.1565521,-0.8212631,5.13];
ybs[3237]=['',2.1862701,1.089813,5.71];
ybs[3238]=['',2.1809299,0.9438486,6.27];
ybs[3239]=['',2.1561862,-0.8771858,5.51];
ybs[3240]=['',2.1715987,0.2035421,7.13];
ybs[3241]=['β Cnc',2.1713051,0.1591966,3.52];
ybs[3242]=['',2.1600957,-0.8010667,5.83];
ybs[3243]=['',2.1672576,-0.5408735,6.21];
ybs[3244]=['',2.1757164,0.153614,6.29];
ybs[3245]=['',2.1675171,-0.6277382,6.16];
ybs[3246]=['30 Lyn',2.1907162,1.0066599,5.89];
ybs[3247]=['',2.1720912,-0.3732326,6.6];
ybs[3248]=['',2.16412,-0.8816205,6.44];
ybs[3249]=['21 Pup',2.1743552,-0.2853541,6.16];
ybs[3250]=['',2.1905621,0.933899,6.49];
ybs[3251]=['',2.1788895,-0.2216035,5.98];
ybs[3252]=['',2.1623464,-1.0991981,5.16];
ybs[3253]=['',2.1764477,-0.524788,6.45];
ybs[3254]=['χ Cnc',2.1873659,0.4738924,5.14];
ybs[3255]=['',2.2010432,1.0570446,6.41];
ybs[3256]=['',2.1883884,0.3609678,5.83];
ybs[3257]=['',2.182687,-0.178568,6.32];
ybs[3258]=['',2.1776169,-0.6198814,5.58];
ybs[3259]=['',2.1771825,-0.6534348,6.7];
ybs[3260]=['λ Cnc',2.1893039,0.4181162,5.98];
ybs[3261]=['',2.1856231,0.0677565,6.05];
ybs[3262]=['',2.1787212,-0.640963,4.45];
ybs[3263]=['',2.1871719,-0.0170206,6.18];
ybs[3264]=['',2.1873368,-0.0941597,6.13];
ybs[3265]=['',2.1828937,-0.6048559,6.43];
ybs[3266]=['',2.1744614,-1.0337872,6.42];
ybs[3267]=['31 Lyn',2.2001292,0.7526063,4.25];
ybs[3268]=['',2.1875493,-0.4012608,6.13];
ybs[3269]=['',2.2050121,0.9276845,5.51];
ybs[3270]=['',2.192032,-0.0291198,6.5];
ybs[3271]=['',2.1915569,-0.351603,5.58];
ybs[3272]=['',2.1752719,-1.1463,5.07];
ybs[3273]=['',2.1940655,-0.3080998,5.75];
ybs[3274]=['',2.1912246,-0.5780641,4.83];
ybs[3275]=['',2.1909348,-0.6379286,5.2];
ybs[3276]=['20 Cnc',2.2014543,0.3187871,5.95];
ybs[3277]=['',2.1969929,-0.1090106,6.15];
ybs[3278]=['',2.1910242,-0.6926691,6.16];
ybs[3279]=['',2.2082549,0.7319449,6.02];
ybs[3280]=['',2.1986849,-0.1328226,5.96];
ybs[3281]=['22 Pup',2.1980005,-0.2290135,6.11];
ybs[3282]=['21 Cnc',2.2036513,0.1843882,6.08];
ybs[3283]=['',2.1978003,-0.461026,5.9];
ybs[3284]=['',2.2095402,0.609881,6.06];
ybs[3285]=['',2.1888816,-1.0129734,5.97];
ybs[3286]=['',2.195441,-0.8474777,4.82];
ybs[3287]=['',2.2062129,-0.0835049,6.01];
ybs[3288]=['',2.1993119,-0.6693823,6.32];
ybs[3289]=['1 Hya',2.2061383,-0.0666477,5.61];
ybs[3290]=['',2.1878162,-1.1200142,6.12];
ybs[3291]=['25 Cnc',2.2121807,0.2963231,6.14];
ybs[3292]=['',2.1969074,-0.9108987,5.85];
ybs[3293]=['κ1 Vol',2.1805526,-1.2493132,5.37];
ybs[3294]=['κ2 Vol',2.1814092,-1.2491449,5.65];
ybs[3295]=['',2.2326374,1.1733454,5.88];
ybs[3296]=['φ1 Cnc',2.2152837,0.485643,5.57];
ybs[3297]=['',2.210702,0.035505,5.73];
ybs[3298]=['',2.2122591,0.1308364,5.13];
ybs[3299]=['ε Car',2.194461,-1.0398024,1.86];
ybs[3300]=['',2.207019,-0.4052871,5.68];
ybs[3301]=['',2.2210893,0.7955951,6.32];
ybs[3302]=['φ2 Cnc',2.216638,0.4689001,6.32];
ybs[3303]=['φ2 Cnc',2.2166525,0.4689147,6.3];
ybs[3304]=['24 Cnc',2.216047,0.4270083,7.02];
ybs[3305]=['24 Cnc',2.2160689,0.4270277,7.81];
ybs[3306]=['',2.2108322,-0.0693654,3.9];
ybs[3307]=['',2.207615,-0.4208652,5.28];
ybs[3308]=['',2.208834,-0.3685024,6.01];
ybs[3309]=['',2.210424,-0.3055614,6.44];
ybs[3310]=['α Cha',2.1727389,-1.3436324,4.07];
ybs[3311]=['27 Cnc',2.2159754,0.2196678,5.5];
ybs[3312]=['',2.2116785,-0.2617605,5.98];
ybs[3313]=['2 Hya',2.2142917,-0.0707866,5.59];
ybs[3314]=['',2.2063728,-0.7476476,5.98];
ybs[3315]=['ο UMa',2.2338392,1.0585097,3.36];
ybs[3316]=['',2.2151174,-0.2199604,5.54];
ybs[3317]=['',2.2178602,-0.1130681,6.59];
ybs[3318]=['',2.2103867,-0.7369009,5.47];
ybs[3319]=['',2.2124077,-0.6829054,6.53];
ybs[3320]=['',2.2124513,-0.6829201,7.25];
ybs[3321]=['28 Cnc',2.2245378,0.4201977,6.1];
ybs[3322]=['',2.2082838,-0.9040084,5.17];
ybs[3323]=['',2.2152272,-0.5110965,6.73];
ybs[3324]=['',2.246613,1.208623,6.31];
ybs[3325]=['29 Cnc',2.2242654,0.2468188,5.95];
ybs[3326]=['η Vol',2.1898409,-1.2822282,5.29];
ybs[3327]=['',2.2186042,-0.3649933,6.56];
ybs[3328]=['',2.2169981,-0.5539957,6.33];
ybs[3329]=['',2.2232054,-0.0451395,6.39];
ybs[3330]=['',2.2223382,-0.1550747,6.43];
ybs[3331]=['',2.2198991,-0.4572992,6.62];
ybs[3332]=['θ Cha',2.1816825,-1.3535035,4.35];
ybs[3333]=['',2.2121506,-0.9228545,6.05];
ybs[3334]=['',2.2245811,-0.1713486,6];
ybs[3335]=['',2.219983,-0.6140544,5.75];
ybs[3336]=['',2.2230869,-0.4038825,6.51];
ybs[3337]=['',2.2244212,-0.3668594,6.67];
ybs[3338]=['',2.2084376,-1.128682,5.97];
ybs[3339]=['β Vol',2.2076497,-1.1554912,3.77];
ybs[3340]=['',2.2368825,0.6491852,6.18];
ybs[3341]=['',2.2165039,-0.9613315,6.53];
ybs[3342]=['',2.2173223,-0.9277691,5.09];
ybs[3343]=['',2.2431216,0.9257914,6.24];
ybs[3344]=['',2.2652147,1.3029061,6.31];
ybs[3345]=['',2.2267047,-0.478254,6.7];
ybs[3346]=['2 UMa',2.2533071,1.1357444,5.47];
ybs[3347]=['υ1 Cnc',2.2371606,0.4190676,5.75];
ybs[3348]=['',2.2245294,-0.7719561,5.79];
ybs[3349]=['θ Cnc',2.2373502,0.31458,5.35];
ybs[3350]=['',2.2241014,-0.8377303,5.33];
ybs[3351]=['',2.2259549,-0.7818098,4.99];
ybs[3352]=['',2.2438245,0.662274,5.9];
ybs[3353]=['',2.2384833,0.1700649,6.83];
ybs[3354]=['',2.2309671,-0.5625069,5.65];
ybs[3355]=['',2.227175,-0.8098583,5.99];
ybs[3356]=['',2.230855,-0.642123,6.69];
ybs[3357]=['32 Lyn',2.2457038,0.6346949,6.24];
ybs[3358]=['η Cnc',2.2422779,0.3555295,5.33];
ybs[3359]=['',2.2359303,-0.342918,5.42];
ybs[3360]=['',2.2258788,-0.9644783,6.36];
ybs[3361]=['υ2 Cnc',2.2436731,0.4191204,6.36];
ybs[3362]=['',2.2136069,-1.2245534,5.53];
ybs[3363]=['',2.2311674,-0.7820313,6.3];
ybs[3364]=['34 Cnc',2.241784,0.1744521,6.46];
ybs[3365]=['',2.234755,-0.6830232,6.31];
ybs[3366]=['',2.2406065,-0.2635467,6.38];
ybs[3367]=['',2.2332877,-0.8366539,6.39];
ybs[3368]=['',2.2466148,0.23014,6.28];
ybs[3369]=['33 Lyn',2.2516608,0.6343898,5.78];
ybs[3370]=['',2.2462613,0.0817777,5.87];
ybs[3371]=['',2.2774136,1.2837952,6.15];
ybs[3372]=['',2.2485298,0.146269,6.03];
ybs[3373]=['',2.242602,-0.4306993,6.19];
ybs[3374]=['',2.2341985,-0.9505821,6.34];
ybs[3375]=['',2.2473894,-0.0387974,5.81];
ybs[3376]=['',2.241392,-0.5510284,6.38];
ybs[3377]=['',2.2417869,-0.6057112,6.36];
ybs[3378]=['',2.2368557,-0.9299574,5.69];
ybs[3379]=['35 Cnc',2.2536566,0.340657,6.58];
ybs[3380]=['',2.2431793,-0.6709403,6.49];
ybs[3381]=['',2.2444973,-0.6792811,5.96];
ybs[3382]=['',2.2435069,-0.8210394,6.24];
ybs[3383]=['π1 UMa',2.2732689,1.1335466,5.64];
ybs[3384]=['',2.2535799,0.0466318,6.33];
ybs[3385]=['',2.1950146,-1.4133873,5.69];
ybs[3386]=['',2.2570409,0.2660146,6.32];
ybs[3387]=['',2.2555831,0.1142846,5.99];
ybs[3388]=['',2.255605,0.1143281,7.25];
ybs[3389]=['',2.2486428,-0.5701995,6.43];
ybs[3390]=['3 Hya',2.2535286,-0.1405695,5.72];
ybs[3391]=['',2.2482643,-0.6576886,6.3];
ybs[3392]=['',2.2684588,0.9307555,5.66];
ybs[3393]=['',2.2725236,1.0448603,6.48];
ybs[3394]=['',2.2530053,-0.4697624,5.96];
ybs[3395]=['π2 UMa',2.2775998,1.1214438,4.6];
ybs[3396]=['',2.251311,-0.6988589,6.47];
ybs[3397]=['',2.2711677,0.9224367,6.42];
ybs[3398]=['36 Cnc',2.2611051,0.1672565,5.88];
ybs[3399]=['',2.2486535,-0.8729373,5.01];
ybs[3400]=['',2.2724299,0.9187114,5.91];
ybs[3401]=['',2.2671534,0.5712285,5.94];
ybs[3402]=['δ Hya',2.2634402,0.0982784,4.16];
ybs[3403]=['',2.2622596,-0.0873746,6.19];
ybs[3404]=['37 Cnc',2.2654216,0.1658391,6.53];
ybs[3405]=['',2.2535738,-0.890849,5.8];
ybs[3406]=['',2.2506358,-1.0137016,4.86];
ybs[3407]=['',2.2503086,-1.0174681,5.26];
ybs[3408]=['',2.2660821,-0.1175553,6.51];
ybs[3409]=['',2.2364234,-1.281546,6.12];
ybs[3410]=['σ Hya',2.2681758,0.0570425,4.44];
ybs[3411]=['',2.2615537,-0.5902423,6.48];
ybs[3412]=['η Pyx',2.26346,-0.4595053,5.27];
ybs[3413]=['',2.2605701,-0.7019712,6.55];
ybs[3414]=['34 Lyn',2.2795148,0.7986605,5.37];
ybs[3415]=['',2.2757899,0.5562055,6.1];
ybs[3416]=['',2.2711521,0.1386468,6.45];
ybs[3417]=['',2.2671662,-0.3457493,6.33];
ybs[3418]=['',2.2618023,-0.7515696,4.14];
ybs[3419]=['39 Cnc',2.2745309,0.3479166,6.39];
ybs[3420]=['',2.2756621,0.3420196,6.44];
ybs[3421]=['ε Cnc',2.2760145,0.3398373,6.3];
ybs[3422]=['',2.2690917,-0.3968031,5.05];
ybs[3423]=['6 Hya',2.2732836,-0.2190185,4.98];
ybs[3424]=['',2.2587819,-1.0982661,5.47];
ybs[3425]=['ζ Pyx',2.2713735,-0.51722,4.89];
ybs[3426]=['',2.2696203,-0.6401905,6.13];
ybs[3427]=['',2.2660096,-0.927879,6.47];
ybs[3428]=['',2.2882231,0.8172739,6.22];
ybs[3429]=['',2.2777494,-0.1592767,6.63];
ybs[3430]=['β Pyx',2.2728831,-0.6175304,3.97];
ybs[3431]=['',2.2736252,-0.7040273,5.2];
ybs[3432]=['',2.2688172,-0.9339774,5.48];
ybs[3433]=['9 Hya',2.2805912,-0.2795587,4.88];
ybs[3434]=['',2.2713081,-0.9272665,5.19];
ybs[3435]=['',2.2674225,-1.0540107,6.36];
ybs[3436]=['',2.2745467,-0.7900251,5.71];
ybs[3437]=['',2.2746341,-0.8154635,3.84];
ybs[3438]=['',2.2826305,-0.210146,6.45];
ybs[3439]=['',2.2727748,-0.9249465,3.62];
ybs[3440]=['',2.2727534,-0.9265755,5.61];
ybs[3441]=['γ Cnc',2.2884363,0.3733919,4.66];
ybs[3442]=['45 Cnc',2.2878345,0.220017,5.64];
ybs[3443]=['',2.2928246,0.6430294,6.33];
ybs[3444]=['',2.2771812,-0.8271272,4.77];
ybs[3445]=['',2.2765175,-0.8551485,5.9];
ybs[3446]=['η Hya',2.2876725,0.0580117,4.3];
ybs[3447]=['',2.274289,-1.0056415,6.34];
ybs[3448]=['',2.2804768,-0.7938642,5.23];
ybs[3449]=['',2.2735897,-1.0443142,4.33];
ybs[3450]=['',2.2910567,0.0743449,6.37];
ybs[3451]=['',2.2893488,-0.1275583,4.62];
ybs[3452]=['θ Vol',2.2652227,-1.2297584,5.2];
ybs[3453]=['δ Cnc',2.2944441,0.3155349,3.94];
ybs[3454]=['',2.2817161,-0.8407864,5.51];
ybs[3455]=['',2.2852901,-0.6286321,6.42];
ybs[3456]=['46 Cnc',2.2977698,0.5344575,6.13];
ybs[3457]=['49 Cnc',2.2945032,0.1746429,5.66];
ybs[3458]=['',2.281598,-0.9280675,5.52];
ybs[3459]=['',2.2820779,-0.9283106,4.86];
ybs[3460]=['α Pyx',2.2881994,-0.5805187,3.68];
ybs[3461]=['10 Hya',2.2955713,0.0978273,6.13];
ybs[3462]=['',2.3154202,1.1629311,6.2];
ybs[3463]=['',2.2815345,-0.9747405,6.29];
ybs[3464]=['',2.2967729,-0.0467121,6.41];
ybs[3465]=['',2.2943998,-0.3707632,6.11];
ybs[3466]=['ι Cnc',2.3034149,0.5007206,6.57];
ybs[3467]=['ι Cnc',2.3035455,0.5006283,4.02];
ybs[3468]=['',2.2877503,-0.8708783,5.16];
ybs[3469]=['',2.2913285,-0.7456802,4.07];
ybs[3470]=['',2.29982,-0.0370834,5.7];
ybs[3471]=['',2.293605,-0.6496564,5.76];
ybs[3472]=['',2.2999036,-0.1934215,6.25];
ybs[3473]=['50 Cnc',2.3040821,0.2100298,5.87];
ybs[3474]=['ε Hya',2.3032472,0.1107022,3.38];
ybs[3475]=['',2.2982058,-0.4444171,6.1];
ybs[3476]=['12 Hya',2.3009693,-0.2377788,4.32];
ybs[3477]=['δ Vel',2.2919158,-0.9561538,1.96];
ybs[3478]=['',2.3050965,-0.0344442,5.29];
ybs[3479]=['',2.2982691,-0.804901,3.91];
ybs[3480]=['',2.3001228,-0.7191013,6.21];
ybs[3481]=['',2.2932484,-1.0262601,6.21];
ybs[3482]=['',2.3022437,-0.6056094,6.37];
ybs[3483]=['',2.2867833,-1.191825,6.32];
ybs[3484]=['ρ Hya',2.3104596,0.1005492,4.36];
ybs[3485]=['',2.3086046,-0.1158061,6.09];
ybs[3486]=['',2.3003806,-0.8026546,5.46];
ybs[3487]=['',2.2898007,-1.1501837,6.05];
ybs[3488]=['',2.3038848,-0.8068971,5.75];
ybs[3489]=['',2.3056788,-0.7297851,6.36];
ybs[3490]=['',2.3005112,-0.9921448,4.49];
ybs[3491]=['',2.3204455,0.5795845,6.25];
ybs[3492]=['14 Hya',2.314271,-0.0614376,5.31];
ybs[3493]=['',2.3077053,-0.7424709,6.43];
ybs[3494]=['η Cha',2.271471,-1.3794569,5.47];
ybs[3495]=['',2.3064795,-0.923751,6.3];
ybs[3496]=['',2.3209225,0.32733,6.16];
ybs[3497]=['5 UMa',2.3345854,1.0800723,5.73];
ybs[3498]=['',2.3330721,1.0293531,6.25];
ybs[3499]=['',2.3154818,-0.3687147,6.47];
ybs[3500]=['35 Lyn',2.3270276,0.7618118,5.15];
ybs[3501]=['',2.3281814,0.7894882,5.99];
ybs[3502]=['54 Cnc',2.3220177,0.2665618,6.38];
ybs[3503]=['',2.3279051,0.7317181,5.99];
ybs[3504]=['',2.3155814,-0.5734761,5.21];
ybs[3505]=['',2.3164762,-0.515576,5.87];
ybs[3506]=['',2.3143851,-0.7050725,5.48];
ybs[3507]=['',2.3179165,-0.5008301,6.17];
ybs[3508]=['',2.3153974,-0.6844984,6.39];
ybs[3509]=['γ Pyx',2.3186959,-0.4849826,4.01];
ybs[3510]=['σ1 Cnc',2.3293097,0.5654149,5.66];
ybs[3511]=['',2.3147656,-0.7921214,4.93];
ybs[3512]=['53 Cnc',2.3287305,0.4918501,6.23];
ybs[3513]=['ρ1 Cnc',2.3292561,0.4931001,5.95];
ybs[3514]=['15 Hya',2.3238335,-0.1266252,5.54];
ybs[3515]=['',2.279658,-1.3813262,6.05];
ybs[3516]=['',2.3173484,-0.7359595,6];
ybs[3517]=['',2.3277669,0.0918361,6.33];
ybs[3518]=['',2.3180362,-0.8134387,5.1];
ybs[3519]=['',2.3353144,0.6188861,6.14];
ybs[3520]=['',2.3277627,-0.2323349,6.13];
ybs[3521]=['',2.3222001,-0.7432,6.55];
ybs[3522]=['6 UMa',2.3490465,1.1261576,5.58];
ybs[3523]=['57 Cnc',2.3365174,0.5323354,5.39];
ybs[3524]=['',2.3268439,-0.5687559,6.5];
ybs[3525]=['',2.3276027,-0.6392054,6.42];
ybs[3526]=['',2.328199,-0.6772303,5.82];
ybs[3527]=['',2.3218674,-1.007254,5.59];
ybs[3528]=['',2.3162545,-1.167109,5.35];
ybs[3529]=['',2.3357627,-0.0962251,6];
ybs[3530]=['',2.3270453,-0.8453914,5.91];
ybs[3531]=['ρ2 Cnc',2.3426067,0.4860414,5.22];
ybs[3532]=['',2.341071,0.2993611,6.64];
ybs[3533]=['',2.3269711,-0.9111904,6.39];
ybs[3534]=['',2.2912457,-1.3889311,5.79];
ybs[3535]=['',2.3117538,-1.2675955,6.11];
ybs[3536]=['',2.3484151,0.7950345,5.74];
ybs[3537]=['',2.3467469,0.7002604,5.89];
ybs[3538]=['ζ Hya',2.3408268,0.1023862,3.11];
ybs[3539]=['',2.3326903,-0.7073146,6.47];
ybs[3540]=['',2.3283073,-0.9900863,6.03];
ybs[3541]=['60 Cnc',2.3432988,0.2015272,5.41];
ybs[3542]=['',2.3323437,-0.8307673,5.33];
ybs[3543]=['17 Hya',2.340906,-0.1404913,6.91];
ybs[3544]=['17 Hya',2.3409132,-0.1405058,6.67];
ybs[3545]=['',2.33939,-0.319754,5.75];
ybs[3546]=['σ2 Cnc',2.3483615,0.5729992,5.45];
ybs[3547]=['δ Pyx',2.3405,-0.4845246,4.89];
ybs[3548]=['',2.3461137,0.0725529,6.14];
ybs[3549]=['',2.348723,0.2978229,6.17];
ybs[3550]=['',2.3423948,-0.4170944,6.39];
ybs[3551]=['',2.3312349,-1.0547503,5.78];
ybs[3552]=['ο1 Cnc',2.3491611,0.2660379,5.2];
ybs[3553]=['',2.338927,-0.7875071,6.26];
ybs[3554]=['61 Cnc',2.3527823,0.5262762,6.29];
ybs[3555]=['',2.3453761,-0.293025,5.96];
ybs[3556]=['ο2 Cnc',2.3506436,0.2705494,5.67];
ybs[3557]=['',2.3550706,0.6234685,6.51];
ybs[3558]=['',2.3509769,0.1624499,6.19];
ybs[3559]=['',2.3361996,-1.0178582,6.38];
ybs[3560]=['ι UMa',2.3588823,0.8370773,3.14];
ybs[3561]=['',2.3377766,-0.9607104,5.71];
ybs[3562]=['',2.3366065,-1.0598293,3.84];
ybs[3563]=['α Cnc',2.3544717,0.2055547,4.25];
ybs[3564]=['',2.3526858,0.0255069,6.59];
ybs[3565]=['',2.3428562,-0.9215883,4.69];
ybs[3566]=['σ3 Cnc',2.3596788,0.5644019,5.2];
ybs[3567]=['ρ UMa',2.3752471,1.1789317,4.76];
ybs[3568]=['',2.3576617,0.3151036,6.38];
ybs[3569]=['',2.3548348,-0.2829737,5.86];
ybs[3570]=['',2.364821,0.7278305,3.97];
ybs[3571]=['',2.3640965,0.6549058,6.44];
ybs[3572]=['',2.4404136,1.4677267,6.33];
ybs[3573]=['',2.3451748,-1.0351401,4.92];
ybs[3574]=['',2.3501364,-0.8491624,5.87];
ybs[3575]=['',2.3588257,-0.3366531,6.18];
ybs[3576]=['',2.3567805,-0.5041682,6.25];
ybs[3577]=['',2.369251,0.6917056,6.36];
ybs[3578]=['66 Cnc',2.3677743,0.5614914,5.82];
ybs[3579]=['',2.3543254,-0.8258051,5.18];
ybs[3580]=['67 Cnc',2.369435,0.485572,6.07];
ybs[3581]=['',2.3675626,0.09703,6.07];
ybs[3582]=['',2.359937,-0.7214276,4.45];
ybs[3583]=['',2.3801764,0.9459952,5.75];
ybs[3584]=['',2.3610743,-0.7549299,6.07];
ybs[3585]=['κ UMa',2.3780836,0.8216043,3.6];
ybs[3586]=['ν Cnc',2.373357,0.4253526,5.45];
ybs[3587]=['',2.3693474,-0.0098498,5.67];
ybs[3588]=['',2.3652678,-0.4667912,6.2];
ybs[3589]=['',2.355809,-1.0326099,5.16];
ybs[3590]=['',2.3729403,0.1259465,5.85];
ybs[3591]=['',2.3654004,-0.7320915,5.55];
ybs[3592]=['70 Cnc',2.3796809,0.4854803,6.38];
ybs[3593]=['',2.368825,-0.6891272,6.27];
ybs[3594]=['',2.3858966,0.8455676,5.95];
ybs[3595]=['',2.3615719,-1.0654352,5.79];
ybs[3596]=['',2.3665843,-0.9122793,5.23];
ybs[3597]=['',2.3831163,0.5636421,6.46];
ybs[3598]=['',2.3738463,-0.4465669,6.74];
ybs[3599]=['',2.3923866,1.034302,6.45];
ybs[3600]=['σ1 UMa',2.4005571,1.1656954,5.14];
ybs[3601]=['',2.3621631,-1.2001762,5.88];
ybs[3602]=['',2.3723657,-0.9360482,6.4];
ybs[3603]=['',2.3903527,0.669666,4.56];
ybs[3604]=['ω Hya',2.3869634,0.087428,4.97];
ybs[3605]=['',2.3774209,-0.8234473,3.75];
ybs[3606]=['α Vol',2.3682619,-1.1602551,4];
ybs[3607]=['σ2 UMa',2.4092646,1.1702458,4.8];
ybs[3608]=['',2.3938533,0.3996391,6.4];
ybs[3609]=['',2.3913487,0.0240763,6.17];
ybs[3610]=['15 UMa',2.4011878,0.8992061,4.48];
ybs[3611]=['',2.3968462,0.5664789,6.5];
ybs[3612]=['τ Cnc',2.3964626,0.5161024,5.43];
ybs[3613]=['',2.3795283,-1.0111559,6.44];
ybs[3614]=['κ Cnc',2.3948347,0.1847341,5.24];
ybs[3615]=['τ UMa',2.4110923,1.1070428,4.67];
ybs[3616]=['',2.4003061,0.5898907,5.93];
ybs[3617]=['75 Cnc',2.3998007,0.4633016,5.98];
ybs[3618]=['ξ Cnc',2.4021563,0.3832992,5.14];
ybs[3619]=['κ Pyx',2.3952383,-0.4527729,4.58];
ybs[3620]=['',2.3874185,-0.9754019,6.11];
ybs[3621]=['19 Hya',2.3985493,-0.151378,5.6];
ybs[3622]=['',2.390695,-0.8952713,6.73];
ybs[3623]=['',2.3846358,-1.1271793,6.37];
ybs[3624]=['',2.4260486,1.2491324,6.55];
ybs[3625]=['λ Vel',2.3943852,-0.7594991,2.21];
ybs[3626]=['',2.403697,0.2003671,6.48];
ybs[3627]=['',2.4005888,-0.2171506,5.77];
ybs[3628]=['',2.3981662,-0.4686496,6.15];
ybs[3629]=['',2.3999137,-0.3213605,5.73];
ybs[3630]=['',2.4080277,0.5389313,5.95];
ybs[3631]=['79 Cnc',2.4064721,0.3824353,6.01];
ybs[3632]=['20 Hya',2.4024215,-0.1548449,5.46];
ybs[3633]=['',2.3814865,-1.2325791,4.71];
ybs[3634]=['',2.3788316,-1.2685974,4.48];
ybs[3635]=['ε Pyx',2.403362,-0.5314449,5.59];
ybs[3636]=['',2.4342837,1.2716414,5.96];
ybs[3637]=['',2.4055098,-0.4059826,6.53];
ybs[3638]=['',2.4017556,-0.8640933,6.48];
ybs[3639]=['16 UMa',2.4257687,1.070541,5.13];
ybs[3640]=['',2.4129452,0.0939574,6.35];
ybs[3641]=['π1 Cnc',2.4147732,0.2602462,6.51];
ybs[3642]=['',2.4141657,0.0660111,6.14];
ybs[3643]=['36 Lyn',2.4222428,0.7527979,5.32];
ybs[3644]=['',2.4125657,-0.3461466,5.73];
ybs[3645]=['',2.4077632,-0.7845724,5];
ybs[3646]=['21 Hya',2.414863,-0.1255738,6.11];
ybs[3647]=['',2.4106665,-0.6866777,6];
ybs[3648]=['',2.4207271,0.3699711,6.48];
ybs[3649]=['',2.4097785,-0.814522,5.79];
ybs[3650]=['',2.4064028,-1.0306431,3.44];
ybs[3651]=['17 UMa',2.4318411,0.9888174,5.27];
ybs[3652]=['',2.4140979,-0.7626866,5.57];
ybs[3653]=['18 UMa',2.4332028,0.9413522,4.83];
ybs[3654]=['',2.4074314,-1.0891181,3.97];
ybs[3655]=['',2.4281705,0.6029679,5.97];
ybs[3656]=['θ Hya',2.4235102,0.0388928,3.88];
ybs[3657]=['',2.4521748,1.2902984,6.5];
ybs[3658]=['',2.4183289,-0.675474,6.31];
ybs[3659]=['',2.4176512,-0.7393038,6.29];
ybs[3660]=['π2 Cnc',2.4275804,0.2592742,5.34];
ybs[3661]=['',2.4185761,-0.8277061,5.92];
ybs[3662]=['',2.421182,-0.7719849,5.85];
ybs[3663]=['',2.4149364,-1.0384649,5.54];
ybs[3664]=['',2.4224095,-0.7559586,5.25];
ybs[3665]=['',2.4276875,-0.2637337,6.35];
ybs[3666]=['',2.4385993,0.8155989,5.97];
ybs[3667]=['',2.4250139,-0.6577871,5.86];
ybs[3668]=['ζ Oct',2.3271675,-1.492195,5.42];
ybs[3669]=['',2.4212516,-0.9713699,5.27];
ybs[3670]=['',2.4259535,-0.7965955,6.25];
ybs[3671]=['23 Hya',2.433479,-0.1123919,5.24];
ybs[3672]=['',2.427854,-0.6746769,4.94];
ybs[3673]=['24 Hya',2.4333947,-0.1541343,5.47];
ybs[3674]=['',2.4285082,-0.6544902,4.62];
ybs[3675]=['β Car',2.4147831,-1.218283,1.68];
ybs[3676]=['',2.4421212,0.6157004,5.75];
ybs[3677]=['',2.4351628,-0.25587,5.84];
ybs[3678]=['',2.4295963,-0.7851345,6.04];
ybs[3679]=['',2.438955,0.1992153,6.41];
ybs[3680]=['38 Lyn',2.4439647,0.6408018,3.82];
ybs[3681]=['',2.42537,-1.0205744,6.02];
ybs[3682]=['',2.4310012,-0.7740922,5.12];
ybs[3683]=['',2.4267118,-1.0064292,6.32];
ybs[3684]=['',2.433684,-0.689195,5.33];
ybs[3685]=['',2.4083947,-1.3395035,6.14];
ybs[3686]=['',2.4294264,-1.0057928,4.34];
ybs[3687]=['',2.4528124,0.8932285,6.13];
ybs[3688]=['',2.4574783,0.9880475,5.47];
ybs[3689]=['ι Car',2.4331607,-1.0360599,2.25];
ybs[3690]=['',2.4362155,-0.9526368,6.33];
ybs[3691]=['',2.4533549,0.6649772,6.12];
ybs[3692]=['',2.4458757,-0.1989957,6.62];
ybs[3693]=['',2.4381412,-0.8925272,5.26];
ybs[3694]=['',2.4457323,-0.2778893,5.78];
ybs[3695]=['α Lyn',2.4535161,0.5987271,3.13];
ybs[3696]=['26 Hya',2.4467841,-0.2105306,4.79];
ybs[3697]=['',2.4552048,0.5727099,6.16];
ybs[3698]=['',2.4407639,-0.901422,5.87];
ybs[3699]=['27 Hya',2.4499389,-0.1683122,4.8];
ybs[3700]=['',2.446298,-0.5967427,6.39];
ybs[3701]=['',2.4538762,0.2667405,6.53];
ybs[3702]=['',2.4328729,-1.2003685,5.39];
ybs[3703]=['',2.435654,-1.1717727,6.11];
ybs[3704]=['',2.4517251,-0.2741154,6.33];
ybs[3705]=['',2.4504706,-0.5558587,6.82];
ybs[3706]=['',2.4492157,-0.65745,6.05];
ybs[3707]=['',2.4442071,-0.9647141,6.28];
ybs[3708]=['θ Pyx',2.4539428,-0.4547211,4.72];
ybs[3709]=['',2.4869117,1.3091387,6.29];
ybs[3710]=['',2.4319806,-1.3086659,5.29];
ybs[3711]=['',2.4321801,-1.3058785,5.86];
ybs[3712]=['',2.4755994,1.1144162,6.28];
ybs[3713]=['',2.4640358,0.4379785,6.41];
ybs[3714]=['',2.4602545,-0.1732655,6.53];
ybs[3715]=['',2.471113,0.8985773,6.31];
ybs[3716]=['',2.4549613,-0.7379801,5.58];
ybs[3717]=['',2.4680589,0.6370091,6.67];
ybs[3718]=['',2.4497371,-1.0907005,4.81];
ybs[3719]=['',2.4583903,-0.6957425,6.54];
ybs[3720]=['',2.4572065,-0.8052218,5.75];
ybs[3721]=['κ Leo',2.4689676,0.4554112,4.46];
ybs[3722]=['',2.4541987,-0.9704573,5.63];
ybs[3723]=['λ Pyx',2.4613372,-0.5047924,4.69];
ybs[3724]=['κ Vel',2.4554575,-0.9616596,2.5];
ybs[3725]=['',2.4634218,-0.6605367,6.48];
ybs[3726]=['',2.4725935,0.2879131,6.29];
ybs[3727]=['',2.4656556,-0.6896623,6.06];
ybs[3728]=['28 Hya',2.4714953,-0.0908756,5.59];
ybs[3729]=['',2.4638733,-0.9045348,6.08];
ybs[3730]=['',2.4609228,-1.0540237,6.3];
ybs[3731]=['',2.4758145,-0.0271134,6.01];
ybs[3732]=['',2.4635902,-1.0775261,5.99];
ybs[3733]=['',2.4871019,0.7943176,5.41];
ybs[3734]=['29 Hya',2.47945,-0.1625507,6.54];
ybs[3735]=['',2.4768051,-0.504002,6.1];
ybs[3736]=['',2.4752422,-0.7084559,6.2];
ybs[3737]=['',2.4925492,0.9713604,6.45];
ybs[3738]=['α Hya',2.4809608,-0.1526914,1.98];
ybs[3739]=['',2.4794209,-0.3915429,4.69];
ybs[3740]=['',2.4818608,-0.1075321,5.38];
ybs[3741]=['',2.5301197,1.4177892,4.29];
ybs[3742]=['',2.4694974,-1.0827936,5.77];
ybs[3743]=['',2.4738763,-0.9332045,5.11];
ybs[3744]=['ω Leo',2.4851452,0.1564937,5.41];
ybs[3745]=['3 Leo',2.4852499,0.1413382,5.71];
ybs[3746]=['',2.4805211,-0.612571,6.65];
ybs[3747]=['23 UMa',2.5006972,1.0990461,3.67];
ybs[3748]=['',2.4874476,-0.0235158,6.27];
ybs[3749]=['τ1 Hya',2.4879028,-0.0499047,4.6];
ybs[3750]=['',2.4890494,-0.0400693,6.14];
ybs[3751]=['',2.4747975,-1.1348015,6.05];
ybs[3752]=['',2.4895806,-0.0757086,6.26];
ybs[3753]=['',2.4877703,-0.3637102,5.66];
ybs[3754]=['7 LMi',2.4956201,0.5858129,5.85];
ybs[3755]=['ε Ant',2.4875018,-0.6290487,4.51];
ybs[3756]=['',2.4875361,-0.6718529,6.19];
ybs[3757]=['',2.4904205,-0.4090339,6.24];
ybs[3758]=['22 UMa',2.5166713,1.258614,5.72];
ybs[3759]=['8 LMi',2.4992385,0.6110722,5.37];
ybs[3760]=['',2.4906717,-0.4656604,5.48];
ybs[3761]=['24 UMa',2.514425,1.21716,4.56];
ybs[3762]=['',2.4929847,-0.2734587,5.85];
ybs[3763]=['λ Leo',2.4996896,0.3992758,4.31];
ybs[3764]=['',2.5224826,1.2954727,6.46];
ybs[3765]=['θ UMa',2.5056437,0.9003388,3.17];
ybs[3766]=['',2.4840459,-1.0884451,5.92];
ybs[3767]=['',2.475341,-1.2512599,5.47];
ybs[3768]=['',2.5066759,0.8612663,6.76];
ybs[3769]=['6 Leo',2.5004344,0.1679798,5.07];
ybs[3770]=['ζ1 Ant',2.4942295,-0.5581963,7];
ybs[3771]=['ζ1 Ant',2.4942805,-0.5581625,6.18];
ybs[3772]=['ξ Leo',2.5004031,0.195624,4.97];
ybs[3773]=['',2.4823345,-1.165742,5.91];
ybs[3774]=['',2.4905327,-0.9007277,5.45];
ybs[3775]=['',2.4986467,-0.1857626,6.14];
ybs[3776]=['ψ Vel',2.4937034,-0.7078627,3.6];
ybs[3777]=['τ2 Hya',2.5002996,-0.0222757,4.57];
ybs[3778]=['',2.49988,-0.1825934,6.13];
ybs[3779]=['ζ2 Ant',2.4976369,-0.5578611,5.93];
ybs[3780]=['',2.4975711,-0.624935,5.87];
ybs[3781]=['9 LMi',2.5078462,0.6352154,6.18];
ybs[3782]=['',2.5067384,0.4935153,6.53];
ybs[3783]=['',2.4914165,-1.0201872,5.88];
ybs[3784]=['',2.5034537,0.0309386,6.11];
ybs[3785]=['ι Cha',2.4583849,-1.4115451,5.36];
ybs[3786]=['',2.5014676,-0.3401938,5.74];
ybs[3787]=['',2.511814,0.8169918,6.52];
ybs[3788]=['',2.5010979,-0.5012486,6.46];
ybs[3789]=['26 UMa',2.5142348,0.906859,4.5];
ybs[3790]=['10 LMi',2.5109741,0.6336506,4.55];
ybs[3791]=['',2.5047427,-0.1500439,6.12];
ybs[3792]=['',2.5041763,-0.2375134,5.94];
ybs[3793]=['',2.4951525,-0.9970272,3.13];
ybs[3794]=['',2.5095776,0.4077434,6.25];
ybs[3795]=['',2.5060798,-0.1270896,6.24];
ybs[3796]=['',2.5300796,1.2738701,6.42];
ybs[3797]=['',2.5007823,-0.7110612,5.35];
ybs[3798]=['',2.5052181,-0.3701404,5.01];
ybs[3799]=['',2.514741,0.6899137,4.81];
ybs[3800]=['',2.5061729,-0.4006508,5.91];
ybs[3801]=['',2.5161037,0.6958801,6.76];
ybs[3802]=['',2.5043678,-0.6845268,6.43];
ybs[3803]=['',2.4956383,-1.1660635,6.27];
ybs[3804]=['33 Hya',2.511387,-0.1048429,5.56];
ybs[3805]=['11 LMi',2.5172169,0.6233944,5.41];
ybs[3806]=['',2.4991378,-1.0974662,6.1];
ybs[3807]=['',2.5066402,-0.8569004,5.12];
ybs[3808]=['7 Leo',2.5176406,0.2493597,6.36];
ybs[3809]=['',2.5082972,-0.8961772,5.01];
ybs[3810]=['',2.5216794,0.5422555,5.56];
ybs[3811]=['',2.4947222,-1.2770903,5.47];
ybs[3812]=['',2.5155284,-0.3434102,6.31];
ybs[3813]=['',2.5135048,-0.6268544,6.49];
ybs[3814]=['',2.5356937,1.1724886,5.94];
ybs[3815]=['',2.5090888,-1.0353539,4.08];
ybs[3816]=['8 Leo',2.522752,0.2852736,5.69];
ybs[3817]=['10 Leo',2.5232834,0.1176874,5];
ybs[3818]=['',2.5197729,-0.4327615,6.53];
ybs[3819]=['42 Lyn',2.5291345,0.7006889,5.25];
ybs[3820]=['',2.5216868,-0.4431291,5.7];
ybs[3821]=['',2.5183675,-0.8524879,6.17];
ybs[3822]=['34 Hya',2.5257745,-0.1661111,6.4];
ybs[3823]=['',2.522211,-0.5632424,5.63];
ybs[3824]=['',2.5286631,0.0795167,4.68];
ybs[3825]=['',2.5234346,-0.6316123,5.98];
ybs[3826]=['',2.5201093,-0.8630297,4.35];
ybs[3827]=['',2.5196829,-0.9256672,6.19];
ybs[3828]=['',2.5480624,1.2067753,5.69];
ybs[3829]=['27 UMa',2.5516714,1.2593932,5.17];
ybs[3830]=['',2.5215364,-0.9383134,5.45];
ybs[3831]=['',2.5157238,-1.1352142,6.56];
ybs[3832]=['',2.5256095,-0.7554509,5.5];
ybs[3833]=['',2.5645579,1.3620441,6.23];
ybs[3834]=['',2.5285925,-0.6930248,6.7];
ybs[3835]=['ι Hya',2.5346628,-0.0215788,3.91];
ybs[3836]=['37 Hya',2.5341818,-0.1861193,6.31];
ybs[3837]=['',2.5727742,1.3795226,6.17];
ybs[3838]=['',2.5365564,-0.1895932,6.37];
ybs[3839]=['κ Hya',2.5363585,-0.2517802,5.06];
ybs[3840]=['',2.5429183,0.5442626,5.89];
ybs[3841]=['43 Lyn',2.5449936,0.6922594,5.62];
ybs[3842]=['ο Leo',2.5405258,0.1710117,3.52];
ybs[3843]=['13 Leo',2.5430212,0.4506206,6.24];
ybs[3844]=['',2.5484317,0.8436339,6.39];
ybs[3845]=['',2.5504637,0.9471735,6.47];
ybs[3846]=['',2.5303776,-1.0720064,4.52];
ybs[3847]=['13 LMi',2.5479284,0.6108461,6.14];
ybs[3848]=['',2.5404168,-0.4133927,4.77];
ybs[3849]=['',2.5577941,1.1325246,6.17];
ybs[3850]=['ζ Cha',2.5011478,-1.4142926,5.11];
ybs[3851]=['15 Leo',2.551464,0.5215006,5.64];
ybs[3852]=['',2.5445865,-0.4190503,4.94];
ybs[3853]=['',2.5365837,-1.0136417,5.32];
ybs[3854]=['',2.5380701,-1.0010091,5.8];
ybs[3855]=['28 UMa',2.5633623,1.1092958,6.34];
ybs[3856]=['ψ Leo',2.551871,0.2430714,5.35];
ybs[3857]=['',2.5462488,-0.6212681,6.41];
ybs[3858]=['',2.5415257,-0.9653113,6];
ybs[3859]=['',2.5553238,0.3275755,6.5];
ybs[3860]=['',2.5654812,0.9954057,5.2];
ybs[3861]=['θ Ant',2.5530599,-0.4863228,4.79];
ybs[3862]=['',2.5490109,-0.8957536,6.15];
ybs[3863]=['ε Leo',2.5613264,0.4132744,2.98];
ybs[3864]=['',2.5529961,-0.6923009,6.82];
ybs[3865]=['',2.5499416,-0.9422386,5.56];
ybs[3866]=['',2.5623468,0.115423,5.79];
ybs[3867]=['18 Leo',2.5634137,0.2044579,5.63];
ybs[3868]=['',2.5580719,-0.5287981,6.45];
ybs[3869]=['',2.5632412,0.0294983,5.65];
ybs[3870]=['19 Leo',2.5679598,0.200235,6.45];
ybs[3871]=['',2.5739058,0.8015435,5.09];
ybs[3872]=['',2.5685096,0.1978006,6.02];
ybs[3873]=['',2.5583657,-0.9997373,6.46];
ybs[3874]=['',2.5560812,-1.0926253,3.69];
ybs[3875]=['',2.5830028,1.1431341,6.31];
ybs[3876]=['',2.5626172,-0.7827875,5.55];
ybs[3877]=['',2.5593418,-1.027814,6.22];
ybs[3878]=['υ UMa',2.5850243,1.0287302,3.8];
ybs[3879]=['20 Leo',2.5786376,0.3679691,6.09];
ybs[3880]=['υ Car',2.5639494,-1.1373872,3.01];
ybs[3881]=['υ Car',2.5639931,-1.137397,6.26];
ybs[3882]=['',2.5758113,-0.6507044,5.97];
ybs[3883]=['4 Sex',2.5812139,0.0741255,6.24];
ybs[3884]=['φ UMa',2.5895827,0.9419097,4.59];
ybs[3885]=['',2.5715018,-0.9862495,6.06];
ybs[3886]=['23 Leo',2.5836986,0.2263594,6.46];
ybs[3887]=['',2.5775273,-0.634688,6.37];
ybs[3888]=['',2.577637,-0.7998691,5.08];
ybs[3889]=['6 Sex',2.5842472,-0.0757482,6.01];
ybs[3890]=['22 Leo',2.5876387,0.4240865,5.32];
ybs[3891]=['',2.5847639,-0.1095791,6.42];
ybs[3892]=['ν Cha',2.5583039,-1.3416585,5.45];
ybs[3893]=['υ1 Hya',2.5851147,-0.2608123,4.12];
ybs[3894]=['',2.580861,-0.8208457,5.73];
ybs[3895]=['μ Leo',2.5915115,0.4522114,3.88];
ybs[3896]=['7 Sex',2.5886031,0.0411407,6.02];
ybs[3897]=['',2.5885449,-0.0003739,6.35];
ybs[3898]=['',2.5873294,-0.2902768,6.08];
ybs[3899]=['γ Sex',2.5897323,-0.1431528,5.05];
ybs[3900]=['',2.5836439,-0.8079235,5.62];
ybs[3901]=['',2.6028676,1.0649709,6.27];
ybs[3902]=['',2.58516,-0.8141016,4.58];
ybs[3903]=['',2.5823922,-1.0388635,5.79];
ybs[3904]=['',2.5809312,-1.0967973,5.57];
ybs[3905]=['',2.5952624,0.102293,5.95];
ybs[3906]=['',2.5913484,-0.4787331,6.3];
ybs[3907]=['31 UMa',2.6050982,0.8678139,5.27];
ybs[3908]=['',2.6187508,1.2702644,5.83];
ybs[3909]=['',2.5967881,-0.4543089,4.88];
ybs[3910]=['',2.5905181,-0.9681424,6.48];
ybs[3911]=['',2.5982748,-0.3941984,6.24];
ybs[3912]=['',2.6120571,1.000423,5.93];
ybs[3913]=['',2.5998384,-0.3334819,4.94];
ybs[3914]=['',2.5943934,-0.894382,5.93];
ybs[3915]=['',2.5966305,-0.7920546,5.71];
ybs[3916]=['',2.6071751,0.1541994,5.85];
ybs[3917]=['',2.5988868,-0.8786253,5.72];
ybs[3918]=['19 LMi',2.6133542,0.7148369,5.14];
ybs[3919]=['',2.6146457,0.7909127,6.3];
ybs[3920]=['',2.6046332,-0.7142355,6.41];
ybs[3921]=['',2.6079976,-0.4651028,6.28];
ybs[3922]=['',2.6070336,-0.5849769,5.84];
ybs[3923]=['',2.6085227,-0.4812427,6.32];
ybs[3924]=['',2.6685346,1.462884,6.37];
ybs[3925]=['',2.6054401,-0.8976949,6.37];
ybs[3926]=['',2.6162767,0.4827632,6.3];
ybs[3927]=['ν Leo',2.6150471,0.2154817,5.26];
ybs[3928]=['',2.6145505,0.1433904,6.04];
ybs[3929]=['',2.6234881,0.989828,5.48];
ybs[3930]=['φ Vel',2.6074604,-0.9541003,3.54];
ybs[3931]=['',2.6089501,-0.9204362,6.12];
ybs[3932]=['',2.6214124,0.5156818,5.73];
ybs[3933]=['',2.6114387,-0.8467083,6.05];
ybs[3934]=['',2.6027827,-1.2476898,6.35];
ybs[3935]=['12 Sex',2.6214073,0.0573483,6.7];
ybs[3936]=['',2.6182138,-0.4197344,6.21];
ybs[3937]=['η Ant',2.6169293,-0.6281402,5.23];
ybs[3938]=['',2.6084632,-1.1272674,6.58];
ybs[3939]=['',2.6067965,-1.2077693,6.2];
ybs[3940]=['π Leo',2.6236482,0.1386689,4.7];
ybs[3941]=['20 LMi',2.6276043,0.5554401,5.36];
ybs[3942]=['',2.6352511,0.3813456,5.66];
ybs[3943]=['',2.6235572,-0.9956359,6.52];
ybs[3944]=['',2.6439796,0.9388394,5.74];
ybs[3945]=['',2.6285659,-0.9331192,6.2];
ybs[3946]=['',2.6343152,-0.5354172,6.54];
ybs[3947]=['',2.6296296,-1.0026765,6.2];
ybs[3948]=['',2.6463908,0.9122936,6.14];
ybs[3949]=['',2.6384864,-0.1688389,6.12];
ybs[3950]=['',2.6295777,-1.0562776,5.94];
ybs[3951]=['13 Sex',2.6406953,0.0541249,6.45];
ybs[3952]=['',2.6382217,-0.4436021,6.7];
ybs[3953]=['',2.6399273,-0.3176781,5.86];
ybs[3954]=['',2.6361327,-0.8156948,6.12];
ybs[3955]=['',2.6411248,-0.4256086,5.7];
ybs[3956]=['',2.633964,-1.0520542,6.19];
ybs[3957]=['',2.6330529,-1.0865721,6.42];
ybs[3958]=['',2.6409275,-0.6994555,6.43];
ybs[3959]=['',2.647636,0.2732688,6.37];
ybs[3960]=['υ2 Hya',2.6447182,-0.2297714,4.6];
ybs[3961]=['',2.6363499,-1.0818193,6.14];
ybs[3962]=['',2.644797,-0.636768,6.27];
ybs[3963]=['14 Sex',2.6522968,0.0961812,6.21];
ybs[3964]=['21 LMi',2.6556572,0.6133775,4.48];
ybs[3965]=['η Leo',2.6548571,0.2908074,3.52];
ybs[3966]=['',2.6485506,-0.8285155,5.08];
ybs[3967]=['',2.6535223,-0.3009358,5.6];
ybs[3968]=['',2.6480525,-0.912606,6.52];
ybs[3969]=['',2.6592114,0.5498345,6.24];
ybs[3970]=['31 Leo',2.6572337,0.1727287,4.37];
ybs[3971]=['α Sex',2.657209,-0.0082475,4.49];
ybs[3972]=['α Leo',2.6593025,0.2071049,1.35];
ybs[3973]=['μ1 Cha',2.6184475,-1.4366445,5.52];
ybs[3974]=['',2.6568967,-0.6533551,6.36];
ybs[3975]=['',2.6606295,-0.1917381,6.53];
ybs[3976]=['',2.6598067,-0.2742382,6.27];
ybs[3977]=['',2.6712729,0.7079019,6.32];
ybs[3978]=['',2.6657598,-0.2128807,6.24];
ybs[3979]=['17 Sex',2.6666209,-0.1485225,5.91];
ybs[3980]=['',2.6604011,-0.9060387,4.86];
ybs[3981]=['',2.6664316,-0.2254526,5.31];
ybs[3982]=['',2.6634946,-0.6275837,6.13];
ybs[3983]=['',2.672197,0.651013,5.85];
ybs[3984]=['λ Hya',2.6685782,-0.2173921,3.61];
ybs[3985]=['',2.6585127,-1.1504561,5.28];
ybs[3986]=['18 Sex',2.6701341,-0.1487002,5.65];
ybs[3987]=['μ2 Cha',2.634027,-1.4253341,6.6];
ybs[3988]=['34 Leo',2.6735663,0.2313131,6.44];
ybs[3989]=['',2.6617354,-1.0760012,5.6];
ybs[3990]=['',2.6717522,-0.129474,6.25];
ybs[3991]=['',2.6681537,-0.7298352,5.98];
ybs[3992]=['',2.6616983,-1.2005162,5.81];
ybs[3993]=['',2.674657,-0.5010526,6.28];
ybs[3994]=['19 Sex',2.678522,0.078762,5.77];
ybs[3995]=['',2.6773751,-0.3360728,6.44];
ybs[3996]=['',2.6833759,0.4718254,6.04];
ybs[3997]=['',2.671636,-1.0285179,6.4];
ybs[3998]=['',2.6900989,1.0451557,6.25];
ybs[3999]=['',2.6725046,-1.0151233,5.72];
ybs[4000]=['',2.6754522,-0.9121998,6.16];
ybs[4001]=['',2.6802571,-0.473525,6.25];
ybs[4002]=['',2.6861647,0.3676607,6.02];
ybs[4003]=['',2.6805254,-0.5782984,6.38];
ybs[4004]=['22 LMi',2.6890212,0.547432,6.46];
ybs[4005]=['',2.6818844,-0.705951,5.9];
ybs[4006]=['',2.7039847,1.2735692,6.4];
ybs[4007]=['',2.6798486,-0.8959722,5.28];
ybs[4008]=['',2.67785,-1.0475476,6.1];
ybs[4009]=['',2.6826647,-0.705336,6.35];
ybs[4010]=['',2.6802066,-0.9050966,5.78];
ybs[4011]=['',2.7029114,1.2384404,6.66];
ybs[4012]=['',2.6791874,-1.077932,6.41];
ybs[4013]=['',2.6860822,-0.7369538,3.85];
ybs[4014]=['23 LMi',2.6939242,0.5097722,5.35];
ybs[4015]=['',2.6794866,-1.1602103,5.16];
ybs[4016]=['32 UMa',2.7030944,1.134554,5.82];
ybs[4017]=['24 LMi',2.6949072,0.4988098,6.49];
ybs[4018]=['',2.6938388,0.3078328,6.55];
ybs[4019]=['',2.6888991,-0.6391499,6.19];
ybs[4020]=['35 Leo',2.6951155,0.4084112,5.97];
ybs[4021]=['ζ Leo',2.6957751,0.4069125,3.44];
ybs[4022]=['',2.695847,0.4410192,5.84];
ybs[4023]=['λ UMa',2.6979746,0.7472016,3.45];
ybs[4024]=['',2.6928805,-0.1973279,6.08];
ybs[4025]=['37 Leo',2.6955614,0.2378096,5.41];
ybs[4026]=['',2.6894997,-0.7542453,5.6];
ybs[4027]=['ω Car',2.6800616,-1.2241773,3.32];
ybs[4028]=['',2.6879833,-0.9612692,6.16];
ybs[4029]=['39 Leo',2.6981821,0.4014805,5.82];
ybs[4030]=['',2.695381,-0.3625644,6.57];
ybs[4031]=['',2.7023079,0.4766862,6.52];
ybs[4032]=['ε Sex',2.6993815,-0.1426271,5.24];
ybs[4033]=['',2.6911127,-1.0473023,6.22];
ybs[4034]=['',2.7063056,0.8143266,6.43];
ybs[4035]=['',2.6942371,-0.8954902,6.3];
ybs[4036]=['',2.7083761,0.8428804,6];
ybs[4037]=['',2.7165771,1.1980582,5.96];
ybs[4038]=['',2.7059203,0.4294961,6.4];
ybs[4039]=['',2.7011972,-0.507805,5.34];
ybs[4040]=['',2.6955384,-1.0722449,3.4];
ybs[4041]=['',2.7120585,0.9368149,6.45];
ybs[4042]=['',2.7132666,0.9444546,6];
ybs[4043]=['',2.7032306,-0.6441655,6.3];
ybs[4044]=['40 Leo',2.7089845,0.3380237,4.79];
ybs[4045]=['',2.7065149,-0.2204603,6];
ybs[4046]=['',2.7024107,-0.7290509,5.96];
ybs[4047]=['γ1 Leo',2.7100225,0.3444951,2.61];
ybs[4048]=['γ2 Leo',2.7100443,0.3444757,3.8];
ybs[4049]=['',2.7077472,-0.0909191,6.37];
ybs[4050]=['',2.7096668,-0.1599145,6.32];
ybs[4051]=['',2.7025851,-0.9811058,5.81];
ybs[4052]=['',2.7595519,1.4686349,5.5];
ybs[4053]=['',2.7069518,-0.9622503,4.57];
ybs[4054]=['23 Sex',2.714386,0.0381522,6.66];
ybs[4055]=['',2.7040374,-1.130619,5.67];
ybs[4056]=['',2.7101385,-0.8343154,5.65];
ybs[4057]=['',2.7200522,0.7177742,5.76];
ybs[4058]=['',2.7145054,-0.3157087,6.51];
ybs[4059]=['μ UMa',2.720727,0.722486,3.05];
ybs[4060]=['42 Leo',2.7180854,0.2595588,6.12];
ybs[4061]=['',2.7159314,-0.4156445,6.5];
ybs[4062]=['',2.7296365,1.1425267,4.97];
ybs[4063]=['',2.7164829,-0.3950065,6.51];
ybs[4064]=['',2.7125926,-0.9799458,4.5];
ybs[4065]=['27 LMi',2.7239381,0.5899886,5.9];
ybs[4066]=['',2.7192122,-0.3485587,6.13];
ybs[4067]=['43 Leo',2.7230564,0.1123701,6.07];
ybs[4068]=['',2.7264329,0.5150731,6.39];
ybs[4069]=['',2.7240768,0.0975631,6.54];
ybs[4070]=['',2.7192661,-0.7287449,4.83];
ybs[4071]=['28 LMi',2.7284592,0.5866785,5.5];
ybs[4072]=['25 Sex',2.7248034,-0.0729272,5.97];
ybs[4073]=['',2.7234142,-0.5282487,6.27];
ybs[4074]=['',2.7638327,1.4390719,5.26];
ybs[4075]=['',2.7282833,0.039508,6.32];
ybs[4076]=['',2.7244308,-0.6652191,5.33];
ybs[4077]=['',2.7251439,-0.7340439,6.27];
ybs[4078]=['44 Leo',2.732889,0.1514963,5.61];
ybs[4079]=['',2.7208786,-1.1694715,4.99];
ybs[4080]=['30 LMi',2.7361818,0.5880249,4.74];
ybs[4081]=['',2.7253423,-1.0133067,6.35];
ybs[4082]=['',2.7347841,-0.125043,5.57];
ybs[4083]=['',2.7321819,-0.7430332,6.18];
ybs[4084]=['μ Hya',2.7361777,-0.2956792,3.81];
ybs[4085]=['',2.7303067,-1.0241754,5.95];
ybs[4086]=['',2.7431135,0.7242377,6.02];
ybs[4087]=['',2.7406984,0.3361412,6.15];
ybs[4088]=['',2.7459084,0.8496181,6.44];
ybs[4089]=['',2.7359847,-0.7477631,6.13];
ybs[4090]=['β LMi',2.7448187,0.6388267,4.21];
ybs[4091]=['45 Leo',2.7433572,0.1685536,6.04];
ybs[4092]=['',2.726197,-1.2939181,4];
ybs[4093]=['',2.7481786,0.7872644,6.35];
ybs[4094]=['α Ant',2.7405781,-0.5440673,4.25];
ybs[4095]=['',2.7277159,-1.2928721,6.19];
ybs[4096]=['35 UMa',2.7547309,1.1435493,6.32];
ybs[4097]=['',2.7384705,-0.959624,5.58];
ybs[4098]=['',2.7569705,1.1196608,6.12];
ybs[4099]=['',2.7479002,-0.0671567,6.05];
ybs[4100]=['',2.740929,-1.0078213,4.66];
ybs[4101]=['',2.7439838,-0.8641247,6.1];
ybs[4102]=['36 UMa',2.7573242,0.9752005,4.84];
ybs[4103]=['32 LMi',2.7545534,0.6775317,5.77];
ybs[4104]=['',2.7429392,-1.0270312,3.82];
ybs[4105]=['',2.7404951,-1.1485964,6.01];
ybs[4106]=['δ Sex',2.7511647,-0.0496478,5.21];
ybs[4107]=['',2.7507884,-0.5195679,5.58];
ybs[4108]=['δ Ant',2.7512375,-0.5360373,5.56];
ybs[4109]=['β Sex',2.7547415,-0.0129598,5.09];
ybs[4110]=['',2.746985,-1.1218543,5.29];
ybs[4111]=['',2.7839707,1.4030303,6.52];
ybs[4112]=['',2.7576455,-0.1351447,6.2];
ybs[4113]=['',2.7576448,-0.2390065,5.58];
ybs[4114]=['33 LMi',2.7620492,0.5632796,5.9];
ybs[4115]=['',2.7568406,-0.4640758,6.51];
ybs[4116]=['',2.778489,1.3195827,4.84];
ybs[4117]=['46 Leo',2.7632484,0.2448918,5.46];
ybs[4118]=['',2.7549451,-1.0727099,6.43];
ybs[4119]=['',2.7523204,-1.1709506,6.19];
ybs[4120]=['',2.7609753,-0.4946851,6.05];
ybs[4121]=['',2.7707193,0.931853,6.45];
ybs[4122]=['',2.7682007,0.7037062,4.75];
ybs[4123]=['ρ Leo',2.765866,0.1605807,3.85];
ybs[4124]=['',2.7584058,-0.9393594,4.89];
ybs[4125]=['',2.7612767,-0.7884099,5.74];
ybs[4126]=['',2.7612112,-0.7884583,6.09];
ybs[4127]=['34 LMi',2.7693245,0.6088129,5.58];
ybs[4128]=['',2.752602,-1.2583582,4.74];
ybs[4129]=['',2.7639184,-0.7805966,5.91];
ybs[4130]=['',2.7609242,-1.0784594,3.32];
ybs[4131]=['37 UMa',2.7771028,0.9944234,5.16];
ybs[4132]=['',2.7554841,-1.2798037,4.93];
ybs[4133]=['',2.7655528,-0.8222143,5.02];
ybs[4134]=['',2.7644537,-1.025782,6];
ybs[4135]=['44 Hya',2.7706766,-0.4162882,5.08];
ybs[4136]=['48 Leo',2.774509,0.1195058,5.08];
ybs[4137]=['',2.7672259,-1.0174646,6.14];
ybs[4138]=['49 Leo',2.7755631,0.1491175,5.67];
ybs[4139]=['',2.774819,-0.4063573,6.1];
ybs[4140]=['35 LMi',2.7817325,0.6321623,6.28];
ybs[4141]=['',2.7705588,-1.0662927,6.23];
ybs[4142]=['',2.7778834,-0.3259532,6.49];
ybs[4143]=['',2.7756397,-0.6923593,5.38];
ybs[4144]=['',2.7753727,-0.7639516,6.08];
ybs[4145]=['',2.7807905,-0.1865761,6.57];
ybs[4146]=['φ2 Hya',2.7806642,-0.2871264,6.03];
ybs[4147]=['',2.7796353,-0.467428,6.29];
ybs[4148]=['',2.7818606,-0.2153215,5.7];
ybs[4149]=['',2.7767362,-1.0064323,4.45];
ybs[4150]=['',2.7847186,-0.2069168,6.52];
ybs[4151]=['',2.7562448,-1.4316392,7.07];
ybs[4152]=['',2.7846537,-0.4803033,4.89];
ybs[4153]=['',2.7862691,-0.2354686,4.82];
ybs[4154]=['',2.7799267,-1.0414624,5.08];
ybs[4155]=['',2.79409,0.9348181,5.52];
ybs[4156]=['37 LMi',2.7919542,0.5562187,4.71];
ybs[4157]=['',2.7845528,-0.8435646,3.84];
ybs[4158]=['38 LMi',2.793832,0.6597834,5.85];
ybs[4159]=['',2.7848168,-1.0269553,5.45];
ybs[4160]=['',2.7741425,-1.3337046,6.3];
ybs[4161]=['φ3 Hya',2.7907202,-0.2964225,4.91];
ybs[4162]=['',2.7918976,-0.2190519,6.04];
ybs[4163]=['',2.787485,-1.0011796,5.91];
ybs[4164]=['γ Cha',2.7737615,-1.3738228,4.11];
ybs[4165]=['',2.7914061,-0.748061,6.11];
ybs[4166]=['',2.8066547,1.1926823,5.75];
ybs[4167]=['',2.7905081,-1.0348084,4.66];
ybs[4168]=['38 UMa',2.8070492,1.1450879,5.12];
ybs[4169]=['',2.7915679,-1.0284193,5.92];
ybs[4170]=['',2.7930856,-0.9723322,4.28];
ybs[4171]=['',2.8121953,1.2037228,5];
ybs[4172]=['33 Sex',2.8032136,-0.0322753,6.26];
ybs[4173]=['',2.8003803,-0.6256855,6.37];
ybs[4174]=['',2.8070656,0.5513361,6.02];
ybs[4175]=['',2.7964713,-1.1380924,5.52];
ybs[4176]=['',2.7914764,-1.3020292,6.07];
ybs[4177]=['39 UMa',2.8143496,0.9964294,5.8];
ybs[4178]=['',2.8016108,-1.043436,6.42];
ybs[4179]=['40 LMi',2.8106539,0.4575853,5.51];
ybs[4180]=['',2.80795,-0.2457904,6.24];
ybs[4181]=['',2.8132675,0.8045262,5.18];
ybs[4182]=['41 LMi',2.8122964,0.4028295,5.08];
ybs[4183]=['35 Sex',2.8117745,0.0809813,5.79];
ybs[4184]=['',2.8085577,-0.5728801,5.64];
ybs[4185]=['',2.8207386,1.1746623,6];
ybs[4186]=['',2.8054711,-1.1270301,4.82];
ybs[4187]=['',2.815849,0.3429672,6.27];
ybs[4188]=['',2.8076868,-1.0353921,5.38];
ybs[4189]=['θ Car',2.8086441,-1.1257766,2.76];
ybs[4190]=['',2.8113827,-1.058971,4.57];
ybs[4191]=['36 Sex',2.8196384,0.0415367,6.28];
ybs[4192]=['41 UMa',2.8258994,0.9993309,6.34];
ybs[4193]=['42 LMi',2.8230742,0.5336156,5.24];
ybs[4194]=['',2.8125807,-1.1232387,5.77];
ybs[4195]=['',2.8137441,-1.1182168,4.82];
ybs[4196]=['',2.8014939,-1.3943595,5.97];
ybs[4197]=['',2.8237752,0.1093401,6.37];
ybs[4198]=['51 Leo',2.8252858,0.3278253,5.49];
ybs[4199]=['52 Leo',2.8252883,0.245853,5.48];
ybs[4200]=['η Car',2.8180987,-1.0435727,6.21];
ybs[4201]=['',2.8141162,-1.2386255,6.26];
ybs[4202]=['',2.8150453,-1.2385388,6.46];
ybs[4203]=['',2.8144542,-1.2662699,6.27];
ybs[4204]=['',2.8268788,-0.3037765,5.42];
ybs[4205]=['',2.8369518,1.1348733,6.39];
ybs[4206]=['μ Vel',2.8259233,-0.864434,2.69];
ybs[4207]=['',2.8233985,-1.0596186,6.25];
ybs[4208]=['',2.8302497,-0.268266,6.67];
ybs[4209]=['',2.8231732,-1.12789,5.34];
ybs[4210]=['',2.8241517,-1.1234982,5.23];
ybs[4211]=['',2.8265238,-0.9924932,5.23];
ybs[4212]=['',2.8257143,-1.1255936,4.85];
ybs[4213]=['43 LMi',2.8365168,0.5115046,6.15];
ybs[4214]=['',2.8349491,-0.0360867,5.93];
ybs[4215]=['',2.8326733,-0.5549622,5.88];
ybs[4216]=['',2.8295026,-1.0048965,6.36];
ybs[4217]=['53 Leo',2.8376123,0.1821506,5.34];
ybs[4218]=['',2.8313526,-1.0476825,6];
ybs[4219]=['40 Sex',2.8375964,-0.0721342,6.61];
ybs[4220]=['44 LMi',2.8406017,0.4863356,6.04];
ybs[4221]=['δ1 Cha',2.8162509,-1.4063488,5.47];
ybs[4222]=['ν Hya',2.838934,-0.2845319,3.11];
ybs[4223]=['',2.8394403,-0.1738638,5.86];
ybs[4224]=['δ2 Cha',2.8185091,-1.4075816,4.45];
ybs[4225]=['43 UMa',2.8467982,0.9856417,5.67];
ybs[4226]=['42 UMa',2.8478104,1.0334243,5.58];
ybs[4227]=['41 Sex',2.8419674,-0.1571974,5.79];
ybs[4228]=['',2.840123,-0.5963261,5.61];
ybs[4229]=['',2.8371191,-1.0372915,5.91];
ybs[4230]=['',2.8454683,-0.0558782,5.95];
ybs[4231]=['',2.852458,0.9155295,6.65];
ybs[4232]=['',2.8525362,0.9144531,6.44];
ybs[4233]=['',2.8575998,1.21727,5.93];
ybs[4234]=['',2.8504772,0.0159877,6.38];
ybs[4235]=['',2.8520938,-0.0054226,6.31];
ybs[4236]=['44 UMa',2.8571351,0.9507776,5.1];
ybs[4237]=['46 LMi',2.8555789,0.5952548,3.83];
ybs[4238]=['ω UMa',2.8586289,0.7518964,4.71];
ybs[4239]=['',2.8556224,-0.0412717,6.12];
ybs[4240]=['',2.8508234,-1.0009435,5.25];
ybs[4241]=['',2.8557716,-0.3533998,5.24];
ybs[4242]=['',2.8560703,-0.2714858,6.38];
ybs[4243]=['',2.8569909,-0.0390715,5.45];
ybs[4244]=['48 LMi',2.8615307,0.442986,6.2];
ybs[4245]=['',2.8593548,-0.2420353,5.66];
ybs[4246]=['',2.8628029,0.5921043,5.72];
ybs[4247]=['',2.8550477,-1.0290941,3.78];
ybs[4248]=['46 UMa',2.866151,0.582891,5.03];
ybs[4249]=['54 Leo',2.8654891,0.4300489,4.5];
ybs[4250]=['54 Leo',2.8655254,0.4300344,6.3];
ybs[4251]=['',2.863196,-0.3625863,6.44];
ybs[4252]=['',2.8552674,-1.2362117,5.99];
ybs[4253]=['',2.8621302,-0.7393345,6.11];
ybs[4254]=['',2.8684662,0.731267,6.03];
ybs[4255]=['55 Leo',2.8656511,0.0109468,5.91];
ybs[4256]=['',2.8593057,-1.080991,5.93];
ybs[4257]=['56 Leo',2.867092,0.1060374,5.81];
ybs[4258]=['',2.8483216,-1.3904808,6.33];
ybs[4259]=['',2.8683757,0.3881935,6.14];
ybs[4260]=['50 LMi',2.8696828,0.4431416,6.35];
ybs[4261]=['',2.8628342,-1.0581339,5.92];
ybs[4262]=['',2.886484,1.3554169,6.2];
ybs[4263]=['ι Ant',2.8696497,-0.6500941,4.6];
ybs[4264]=['',2.8712048,-0.8879349,5.91];
ybs[4265]=['',2.881973,0.9035918,6.17];
ybs[4266]=['',2.8738966,-1.0444391,6.11];
ybs[4267]=['47 UMa',2.8824887,0.7037173,5.05];
ybs[4268]=['',2.882774,0.6280184,6];
ybs[4269]=['',2.8704071,-1.3126559,6.13];
ybs[4270]=['',2.8859638,0.7926546,5.47];
ybs[4271]=['',2.8831092,0.2023807,6.55];
ybs[4272]=['',2.8806786,-0.590749,5.71];
ybs[4273]=['',2.8868755,0.8969522,6.43];
ybs[4274]=['',2.8821098,-0.2873533,5.89];
ybs[4275]=['',2.8863469,0.7470189,6.02];
ybs[4276]=['',2.8901617,1.1049793,6.39];
ybs[4277]=['α Crt',2.8832257,-0.3213006,4.08];
ybs[4278]=['49 UMa',2.8884531,0.6824552,5.08];
ybs[4279]=['',2.8850929,-0.2477261,5.88];
ybs[4280]=['',2.8800825,-1.072164,6.16];
ybs[4281]=['58 Leo',2.8868527,0.0612108,4.84];
ybs[4282]=['',2.8838423,-0.7665054,5.81];
ybs[4283]=['',2.88459,-0.7389053,4.39];
ybs[4284]=['59 Leo',2.8876891,0.1045624,4.99];
ybs[4285]=['β UMa',2.8931564,0.9821308,2.37];
ybs[4286]=['',2.8843526,-0.9063162,6.15];
ybs[4287]=['',2.8883946,-0.2775633,6.34];
ybs[4288]=['',2.88703,-0.5576298,6.07];
ybs[4289]=['61 Leo',2.8923333,-0.0452958,4.74];
ybs[4290]=['60 Leo',2.8947243,0.3502722,4.42];
ybs[4291]=['α UMa',2.9015349,1.0758219,1.79];
ybs[4292]=['',2.8946297,-0.4702265,6.23];
ybs[4293]=['',2.8985226,-0.0150659,6.14];
ybs[4294]=['',2.8774187,-1.4253453,6.71];
ybs[4295]=['',2.898454,-0.1992175,5.5];
ybs[4296]=['62 Leo',2.9001289,-0.0019475,5.95];
ybs[4297]=['',2.898336,-0.5597541,6.46];
ybs[4298]=['',2.9000067,-0.2364083,6.34];
ybs[4299]=['51 UMa',2.9044717,0.6655031,6];
ybs[4300]=['χ Leo',2.9063269,0.1261033,4.63];
ybs[4301]=['',2.9035789,-0.8340932,5.67];
ybs[4302]=['η Oct',2.8754346,-1.478364,6.19];
ybs[4303]=['',2.9054327,-0.6268459,5.43];
ybs[4304]=['χ1 Hya',2.9074003,-0.4782999,4.94];
ybs[4305]=['',2.9085754,-0.1954747,6.09];
ybs[4306]=['',2.9059517,-0.8639977,6.13];
ybs[4307]=['χ2 Hya',2.9101447,-0.4781994,5.71];
ybs[4308]=['',2.9104055,-0.8957647,6.3];
ybs[4309]=['65 Leo',2.9145142,0.0321911,5.52];
ybs[4310]=['',2.9131218,-0.5033335,6.77];
ybs[4311]=['',2.9119853,-0.8913003,6.32];
ybs[4312]=['64 Leo',2.918,0.4051325,6.46];
ybs[4313]=['',2.9119356,-1.0260156,6.02];
ybs[4314]=['',2.9152428,-0.5706945,6.59];
ybs[4315]=['',2.912037,-1.091446,4.61];
ybs[4316]=['',2.9113398,-1.1336051,6.41];
ybs[4317]=['',2.9157128,-0.7461246,5.15];
ybs[4318]=['',2.9186179,-0.52859,6.54];
ybs[4319]=['',2.9128705,-1.2389947,5.57];
ybs[4320]=['',2.9275664,1.1710953,6.06];
ybs[4321]=['',2.9201781,-0.5250661,6.49];
ybs[4322]=['67 Leo',2.9230548,0.4284256,5.68];
ybs[4323]=['',2.9253535,0.6317749,5.74];
ybs[4324]=['',2.9222427,-0.4920415,5.44];
ybs[4325]=['ψ UMa',2.9269628,0.7747021,3.01];
ybs[4326]=['',2.9268499,0.752168,5.89];
ybs[4327]=['',2.9211266,-1.0312508,3.91];
ybs[4328]=['',2.9209299,-1.0831258,5.13];
ybs[4329]=['',2.9272592,-0.5668648,5.81];
ybs[4330]=['',2.9385563,1.1896203,6.4];
ybs[4331]=['',2.9356576,0.2493834,6.3];
ybs[4332]=['',2.9313038,-1.0221843,6.88];
ybs[4333]=['β Crt',2.9350725,-0.4003346,4.48];
ybs[4334]=['',2.9405456,0.9561331,6.63];
ybs[4335]=['',2.9393708,0.623115,6.41];
ybs[4336]=['',2.9375541,-0.5680279,6.38];
ybs[4337]=['ψ Crt',2.9388122,-0.3248361,6.13];
ybs[4338]=['',2.9390917,-0.3815449,6.4];
ybs[4339]=['',2.9332792,-1.2487484,6.35];
ybs[4340]=['',2.9386756,-0.8589263,5.36];
ybs[4341]=['',2.9443682,0.7151792,6.33];
ybs[4342]=['',2.9386429,-1.0546893,4.6];
ybs[4343]=['',2.9403998,-0.8700146,6.11];
ybs[4344]=['',2.9417784,-0.7763929,5.8];
ybs[4345]=['',2.9391878,-1.1219234,5.23];
ybs[4346]=['69 Leo',2.9444162,-0.0031693,5.42];
ybs[4347]=['δ Leo',2.9460831,0.3562515,2.56];
ybs[4348]=['',2.9456503,0.1387303,5.79];
ybs[4349]=['θ Leo',2.9466188,0.2673413,3.34];
ybs[4350]=['',2.9434277,-0.93102,5.76];
ybs[4351]=['',2.9426612,-1.0425075,5.74];
ybs[4352]=['72 Leo',2.9508773,0.4011386,4.63];
ybs[4353]=['',2.9549737,0.9191072,6.5];
ybs[4354]=['',2.9490218,-0.7652595,6.21];
ybs[4355]=['73 Leo',2.9536905,0.2303037,5.32];
ybs[4356]=['',2.9541161,0.2222265,6.67];
ybs[4357]=['',2.957669,0.8615685,5.88];
ybs[4358]=['φ Leo',2.957054,-0.0656909,4.47];
ybs[4359]=['',2.958377,-0.1264822,6.14];
ybs[4360]=['',2.9558208,-0.802714,6.31];
ybs[4361]=['75 Leo',2.959832,0.0331325,5.18];
ybs[4362]=['',2.9591259,-0.6654353,6.27];
ybs[4363]=['',2.9611536,-0.6082377,6.45];
ybs[4364]=['ξ UMa',2.9639286,0.5483281,4.87];
ybs[4365]=['ξ UMa',2.9639358,0.5483281,4.41];
ybs[4366]=['',2.9614131,-0.6396052,6.68];
ybs[4367]=['ν UMa',2.965235,0.575642,3.48];
ybs[4368]=['',2.9645197,0.2072129,6.66];
ybs[4369]=['',2.9590374,-1.1857035,6.06];
ybs[4370]=['55 UMa',2.9681289,0.6645025,4.78];
ybs[4371]=['76 Leo',2.966927,0.0268468,5.91];
ybs[4372]=['δ Crt',2.9686797,-0.2598969,3.56];
ybs[4373]=['',2.9763427,1.1691617,6.21];
ybs[4374]=['',2.9677459,-1.1291385,5.99];
ybs[4375]=['',2.9634506,-1.3924395,6.35];
ybs[4376]=['σ Leo',2.9766385,0.1032696,4.05];
ybs[4377]=['',2.968611,-1.3134456,6.27];
ybs[4378]=['',2.9800939,0.9941815,6.43];
ybs[4379]=['',2.9708724,-1.2585025,6.41];
ybs[4380]=['π Cen',2.9755685,-0.9530132,3.89];
ybs[4381]=['',2.9847579,1.1208133,6.02];
ybs[4382]=['56 UMa',2.9842778,0.756951,4.99];
ybs[4383]=['',2.9817328,-0.7811827,6.12];
ybs[4384]=['',2.9860436,0.0003308,6.05];
ybs[4385]=['λ Crt',2.9862201,-0.32974,5.09];
ybs[4386]=['',2.9854277,-0.6331605,5];
ybs[4387]=['',2.9786686,-1.356486,6.43];
ybs[4388]=['',2.9848411,-0.9929552,5.79];
ybs[4389]=['ι Leo',2.9888271,0.1818006,3.94];
ybs[4390]=['79 Leo',2.989272,0.0226022,5.39];
ybs[4391]=['',2.9856474,-1.1356458,5.11];
ybs[4392]=['ε Crt',2.9917008,-0.191502,4.83];
ybs[4393]=['',2.9904237,-0.746686,6.12];
ybs[4394]=['',2.9934408,0.1975266,5.8];
ybs[4395]=['γ Crt',2.9928486,-0.3106114,4.08];
ybs[4396]=['',2.9889436,-1.2630849,5.59];
ybs[4397]=['',2.9980397,0.9728054,5.75];
ybs[4398]=['81 Leo',2.9961941,0.285248,5.57];
ybs[4399]=['',2.9953839,-0.6313891,5.22];
ybs[4400]=['80 Leo',2.9971162,0.0653992,6.37];
ybs[4401]=['',2.9956398,-0.6607931,5.89];
ybs[4402]=['',2.9998671,0.581851,6.32];
ybs[4403]=['',2.996001,-1.1185059,5.17];
ybs[4404]=['83 Leo',3.0011331,0.050616,6.5];
ybs[4405]=['',2.9998728,-1.0686342,5.3];
ybs[4406]=['κ Crt',3.0028166,-0.2176367,5.94];
ybs[4407]=['',3.0008846,-0.9297887,5.81];
ybs[4408]=['τ Leo',3.0062878,0.0478754,4.95];
ybs[4409]=['',3.0060842,-0.0316437,6.25];
ybs[4410]=['',3.0062477,-0.6185737,6.45];
ybs[4411]=['',3.0117374,1.0762608,5.83];
ybs[4412]=['57 UMa',3.0114349,0.6845847,5.31];
ybs[4413]=['',3.0088559,-0.7467785,5.08];
ybs[4414]=['',3.0144691,0.988281,6.28];
ybs[4415]=['',3.0070134,-1.2668911,6.09];
ybs[4416]=['85 Leo',3.014032,0.2670382,5.74];
ybs[4417]=['',3.0165681,0.9468143,6.41];
ybs[4418]=['',3.0136032,-0.4289408,5.76];
ybs[4419]=['',3.0248214,1.4139594,6.15];
ybs[4420]=['',3.0173615,0.812351,6.35];
ybs[4421]=['58 UMa',3.0177738,0.7515407,5.94];
ybs[4422]=['87 Leo',3.0166383,-0.0543987,4.77];
ybs[4423]=['86 Leo',3.0174697,0.3193343,5.52];
ybs[4424]=['λ Dra',3.0220422,1.2080791,3.84];
ybs[4425]=['',3.0193982,0.8345453,6.42];
ybs[4426]=['',3.0206618,0.8495548,6.56];
ybs[4427]=['88 Leo',3.0229672,0.2487295,6.2];
ybs[4428]=['',3.0202831,-1.0714855,6.38];
ybs[4429]=['',3.0259407,1.0641127,5.48];
ybs[4430]=['',3.023011,-0.3645986,6.24];
ybs[4431]=['ο1 Cen',3.0225848,-1.0394398,5.13];
ybs[4432]=['ο2 Cen',3.0227805,-1.0407246,5.15];
ybs[4433]=['',3.0250423,-0.5127194,5.81];
ybs[4434]=['',3.0250569,-0.5126806,5.64];
ybs[4435]=['',3.0255797,-0.4687954,6.16];
ybs[4436]=['',3.0274263,-0.138594,5.95];
ybs[4437]=['',3.0273043,-0.7077265,5.64];
ybs[4438]=['',3.024895,-1.1706892,5.9];
ybs[4439]=['',3.027799,-0.5445529,5.04];
ybs[4440]=['ξ Hya',3.0282316,-0.5580017,3.54];
ybs[4441]=['',3.0293624,-0.2861282,6.05];
ybs[4442]=['',3.032627,0.6405731,6.4];
ybs[4443]=['',3.0308873,-0.7103551,5.39];
ybs[4444]=['',3.0335005,0.1904185,6.55];
ybs[4445]=['89 Leo',3.034342,0.0514271,5.77];
ybs[4446]=['90 Leo',3.0358859,0.2911817,5.95];
ybs[4447]=['',3.0377491,0.9542028,5.63];
ybs[4448]=['',3.0347358,-0.5749959,5.98];
ybs[4449]=['',3.037449,0.3547888,6.45];
ybs[4450]=['',3.0357551,-0.9490687,4.62];
ybs[4451]=['2 Dra',3.0421965,1.2079291,5.2];
ybs[4452]=['',3.0366146,-0.8595771,5.5];
ybs[4453]=['',3.0378246,-0.8287869,5.71];
ybs[4454]=['',3.0402896,0.1884536,6.56];
ybs[4455]=['',3.0428647,0.4828901,5.8];
ybs[4456]=['',3.0409126,-0.8334853,5.25];
ybs[4457]=['λ Cen',3.0400919,-1.1018829,3.13];
ybs[4458]=['θ Crt',3.0443973,-0.1730632,4.7];
ybs[4459]=['',3.0438704,-0.587889,5.74];
ybs[4460]=['',3.0442737,-0.6519039,6.31];
ybs[4461]=['υ Leo',3.0455922,-0.0163619,4.3];
ybs[4462]=['',3.0427151,-1.0675441,5.83];
ybs[4463]=['',3.0457811,-0.5777326,6.29];
ybs[4464]=['',3.0499025,0.8814735,6.14];
ybs[4465]=['',3.045507,-1.0715782,5.15];
ybs[4466]=['',3.048076,-0.835329,5.44];
ybs[4467]=['59 UMa',3.0518563,0.7594261,5.59];
ybs[4468]=['',3.0509254,0.1530747,6.17];
ybs[4469]=['π Cha',3.0461816,-1.3266292,5.65];
ybs[4470]=['60 UMa',3.0528131,0.8154267,6.1];
ybs[4471]=['',3.0541335,1.1210821,6.46];
ybs[4472]=['',3.0526462,0.5848931,6.27];
ybs[4473]=['ω Vir',3.052217,0.1399844,5.36];
ybs[4474]=['',3.0519286,-0.0445017,6.22];
ybs[4475]=['',3.0488813,-1.1821795,5.96];
ybs[4476]=['',3.0536293,0.78531,6.44];
ybs[4477]=['',3.0503674,-1.0810574,5.15];
ybs[4478]=['ι Crt',3.0530576,-0.2324011,5.48];
ybs[4479]=['',3.0544944,-0.4334488,6.42];
ybs[4480]=['',3.0582193,-0.2545095,6.21];
ybs[4481]=['',3.0581615,-0.2920632,6.19];
ybs[4482]=['',3.0563047,-1.1433908,5.17];
ybs[4483]=['',3.0611662,1.009792,6.37];
ybs[4484]=['ο Hya',3.0597251,-0.6083947,4.7];
ybs[4485]=['92 Leo',3.0623991,0.370691,5.26];
ybs[4486]=['61 UMa',3.0636006,0.5949462,5.33];
ybs[4487]=['',3.0617839,-0.9439152,5.96];
ybs[4488]=['',3.0637956,-0.5115587,6.44];
ybs[4489]=['',3.0625027,-1.0856603,4.94];
ybs[4490]=['',3.0666654,0.9609558,6.27];
ybs[4491]=['62 UMa',3.0658642,0.5520883,5.73];
ybs[4492]=['',3.064566,-0.7541499,5.55];
ybs[4493]=['',3.0663753,-0.5692132,5.22];
ybs[4494]=['3 Dra',3.0700444,1.1629335,5.3];
ybs[4495]=['',3.0680769,0.385666,6.59];
ybs[4496]=['',3.0678357,-0.3561814,6.22];
ybs[4497]=['',3.062015,-1.452354,6.33];
ybs[4498]=['',3.0738728,-0.6510799,5.98];
ybs[4499]=['',3.0708967,-1.3861443,6.39];
ybs[4500]=['',3.0759916,-0.1185268,6.07];
ybs[4501]=['',3.0740023,-1.0926337,5.03];
ybs[4502]=['',3.077386,0.4381554,6.02];
ybs[4503]=['',3.0755821,-1.0994212,6.1];
ybs[4504]=['ζ Crt',3.0796521,-0.3222703,4.73];
ybs[4505]=['ξ Vir',3.0819886,0.142147,4.85];
ybs[4506]=['',3.0815002,-0.8584163,6.26];
ybs[4507]=['ν Vir',3.0844932,0.1119719,4.03];
ybs[4508]=['χ UMa',3.0854369,0.8319201,3.71];
ybs[4509]=['',3.0837971,-0.7994292,5.29];
ybs[4510]=['λ Mus',3.0830928,-1.1666223,3.64];
ybs[4511]=['',3.0892927,0.9689086,5.27];
ybs[4512]=['',3.087125,-1.069752,4.11];
ybs[4513]=['',3.0872558,-0.7088567,4.91];
ybs[4514]=['',3.089885,-0.6286834,6.17];
ybs[4515]=['',3.0905335,-0.530596,6.48];
ybs[4516]=['',3.0906752,-1.0089811,5.41];
ybs[4517]=['93 Leo',3.0937884,0.3508968,4.53];
ybs[4518]=['4 Vir',3.0934608,0.1419275,5.32];
ybs[4519]=['',3.0955085,-0.1819912,6.26];
ybs[4520]=['μ Mus',3.094628,-1.1681264,4.72];
ybs[4521]=['',3.0966554,0.247316,5.88];
ybs[4522]=['',3.0970498,-0.4688604,5.11];
ybs[4523]=['',3.0982672,-0.0075507,6.15];
ybs[4524]=['β Leo',3.0984658,0.2523385,2.14];
ybs[4525]=['',3.0992899,0.2815,6.04];
ybs[4526]=['',3.1012722,0.6076825,5.7];
ybs[4527]=['',3.1009959,-1.1153065,4.32];
ybs[4528]=['',3.1020649,-1.2276622,4.97];
ybs[4529]=['',3.1039393,-0.2788674,6.13];
ybs[4530]=['β Vir',3.1055786,0.0288097,3.61];
ybs[4531]=['',3.1043697,-1.0954295,5.7];
ybs[4532]=['',3.105202,-0.4780774,6.48];
ybs[4533]=['',3.1065815,0.2123165,6.35];
ybs[4534]=['',3.1070602,-0.0950748,5.64];
ybs[4535]=['',3.1076382,0.580513,6.27];
ybs[4536]=['',3.1074668,-0.7904188,4.46];
ybs[4537]=['',3.1084845,-0.2147076,6.35];
ybs[4538]=['',3.109892,-0.5401631,5.85];
ybs[4539]=['',3.1104865,-1.1400521,4.9];
ybs[4540]=['',3.1155901,0.6563229,6.45];
ybs[4541]=['',3.111907,-0.9966153,5.57];
ybs[4542]=['β Hya',3.1152027,-0.5937983,4.28];
ybs[4543]=['',3.1175459,-0.61402,6.17];
ybs[4544]=['γ UMa',3.1193219,0.9351584,2.44];
ybs[4545]=['',3.1192922,0.0076419,6.3];
ybs[4546]=['',3.1207613,-1.0039849,6.06];
ybs[4547]=['',3.1218393,-0.6608338,6.46];
ybs[4548]=['',3.1230678,-0.4507835,5.3];
ybs[4549]=['6 Vir',3.1245916,0.1453821,5.58];
ybs[4550]=['65 UMa',3.124815,0.8091841,6.54];
ybs[4551]=['65 UMa',3.125214,0.8090581,7.03];
ybs[4552]=['',3.1254133,0.6395284,6.49];
ybs[4553]=['',3.1242683,-1.1064166,5.91];
ybs[4554]=['95 Leo',3.127315,0.2710942,5.53];
ybs[4555]=['',3.1272586,-0.4990081,5.93];
ybs[4556]=['66 UMa',3.1286534,0.9858404,5.84];
ybs[4557]=['η Crt',3.1287799,-0.3013302,5.18];
ybs[4558]=['',3.1283126,-0.6946983,6.13];
ybs[4559]=['',3.132634,1.0722438,6.22];
ybs[4560]=['',3.1318906,-0.8235619,6.26];
ybs[4561]=['',3.1333424,-0.5834531,6.21];
ybs[4562]=['',3.1341653,0.702137,6.62];
ybs[4563]=['',3.135973,-1.0919306,5.57];
ybs[4564]=['',3.1379815,0.5612937,6.42];
ybs[4565]=['',3.1389665,1.0707699,6.76];
ybs[4566]=['',3.1385399,-0.9849129,5.44];
ybs[4567]=['',3.1389181,-0.7166557,6.79];
ybs[4568]=['',3.1409069,-1.1249222,5.61];
ybs[4569]=['',3.1414041,-0.4541873,6.43];
ybs[4570]=['',3.1420611,0.007268,6.17];
ybs[4571]=['',3.1430947,0.5768902,5.96];
ybs[4572]=['',3.1426032,-0.9042689,6.05];
ybs[4573]=['ε Cha',3.1445325,-1.3672224,4.91];
ybs[4574]=['',3.1459706,0.5920309,6.5];
ybs[4575]=['7 Vir',3.1459514,0.0618047,5.37];
ybs[4576]=['',3.1474852,1.4091601,6.17];
ybs[4577]=['',3.1494153,-0.1843109,5.55];
ybs[4578]=['',3.149272,-0.3831233,6.28];
ybs[4579]=['π Vir',3.1499863,0.1134471,4.66];
ybs[4580]=['',3.1499047,-0.3451042,5.26];
ybs[4581]=['',3.1506718,-0.0328503,6.31];
ybs[4582]=['',3.1526768,-1.0056192,6.16];
ybs[4583]=['',3.1533992,0.6270588,5.59];
ybs[4584]=['67 UMa',3.1553769,0.7492949,5.21];
ybs[4585]=['',3.1566996,-1.4921988,6.05];
ybs[4586]=['',3.1570621,-1.2497083,6.42];
ybs[4587]=['',3.1577179,-1.2096238,5.89];
ybs[4588]=['',3.1586602,-0.136096,6.22];
ybs[4589]=['θ1 Cru',3.1594378,-1.1070081,4.33];
ybs[4590]=['',3.1621802,-0.7426075,5.15];
ybs[4591]=['',3.1626207,-1.2972683,6.44];
ybs[4592]=['2 Com',3.1648212,0.3725416,5.87];
ybs[4593]=['θ2 Cru',3.1651102,-1.1044384,4.72];
ybs[4594]=['',3.1665584,-1.1945603,5.35];
ybs[4595]=['κ Cha',3.1672057,-1.3375028,5.04];
ybs[4596]=['',3.1651257,1.4873651,6.27];
ybs[4597]=['',3.167869,-1.0661041,5.96];
ybs[4598]=['ο Vir',3.1688951,0.1504293,4.12];
ybs[4599]=['',3.1688779,1.3402687,5.8];
ybs[4600]=['',3.170776,1.0963978,6.13];
ybs[4601]=['',3.1719857,-1.1460059,6.33];
ybs[4602]=['',3.1721575,-0.624967,6.23];
ybs[4603]=['',3.1723456,-0.056649,6.37];
ybs[4604]=['',3.1739502,-1.2001741,6.23];
ybs[4605]=['',3.1741706,-1.1488323,6.06];
ybs[4606]=['η Cru',3.1763403,-1.1297112,4.15];
ybs[4607]=['',3.1806222,-1.3173919,5.18];
ybs[4608]=['',3.1815534,-0.8861985,4.47];
ybs[4609]=['',3.1815246,-0.8879778,6.37];
ybs[4610]=['',3.1822397,-0.8518398,5.34];
ybs[4611]=['δ Cen',3.1827419,-0.887265,2.6];
ybs[4612]=['',3.1830124,-1.0639748,6.22];
ybs[4613]=['α Crv',3.182923,-0.4335909,4.02];
ybs[4614]=['',3.1850763,-0.7756268,5.75];
ybs[4615]=['',3.1851191,-0.7216137,5.48];
ybs[4616]=['10 Vir',3.1884476,0.0311325,5.95];
ybs[4617]=['',3.1885625,1.3010972,6.35];
ybs[4618]=['',3.1900559,-0.6077063,6.17];
ybs[4619]=['11 Vir',3.1900483,0.0993605,5.72];
ybs[4620]=['ε Crv',3.1903944,-0.3967784,3];
ybs[4621]=['',3.1923435,-0.6629506,6.06];
ybs[4622]=['3 Com',3.1920794,0.2913857,6.39];
ybs[4623]=['',3.1931124,0.4741606,6.01];
ybs[4624]=['',3.1947147,-1.0714835,6.08];
ybs[4625]=['3 Crv',3.1945014,-0.4139307,5.46];
ybs[4626]=['',3.1944869,-0.7947664,6.61];
ybs[4627]=['',3.1965885,-0.8983806,6.23];
ybs[4628]=['ρ Cen',3.1971547,-0.9159938,3.96];
ybs[4629]=['',3.1934958,1.4241191,6];
ybs[4630]=['4 Com',3.197845,0.4495325,5.66];
ybs[4631]=['68 UMa',3.1972749,0.9937988,6.43];
ybs[4632]=['',3.1985651,0.4960601,6.49];
ybs[4633]=['5 Com',3.1991721,0.3565357,5.57];
ybs[4634]=['',3.2003571,-1.1006881,5.92];
ybs[4635]=['',3.2022624,-1.226371,6.17];
ybs[4636]=['',3.1989036,1.3526726,5.14];
ybs[4637]=['',3.2039268,-0.5975916,6.5];
ybs[4638]=['',3.2048383,-0.6814303,5.76];
ybs[4639]=['',3.2075784,-1.3733561,6.35];
ybs[4640]=['12 Vir',3.204762,0.1771214,5.85];
ybs[4641]=['',3.2056517,-0.5917833,6.33];
ybs[4642]=['',3.2075829,-0.8000203,5.31];
ybs[4643]=['',3.2087551,-1.1261301,6.22];
ybs[4644]=['',3.2102418,0.9306243,6.16];
ybs[4645]=['',3.2116493,-0.3657867,5.83];
ybs[4646]=['δ Cru',3.2124848,-1.0273488,2.8];
ybs[4647]=['',3.2124228,-0.1819744,6.11];
ybs[4648]=['',3.2139703,-0.7335079,6.26];
ybs[4649]=['',3.2118693,1.2232338,5.71];
ybs[4650]=['δ UMa',3.2132758,0.9934178,3.31];
ybs[4651]=['',3.2151062,-0.4095843,6.54];
ybs[4652]=['γ Crv',3.2151916,-0.3081516,2.59];
ybs[4653]=['6 Com',3.2159668,0.2580479,5.1];
ybs[4654]=['',3.2181716,-1.2693525,6.22];
ybs[4655]=['',3.2141926,1.2642639,6.29];
ybs[4656]=['2 CVn',3.2164191,0.707669,5.66];
ybs[4657]=['7 Com',3.2174162,0.4159374,4.95];
ybs[4658]=['',3.2180833,0.5750437,5];
ybs[4659]=['',3.2211304,-1.1485413,6.06];
ybs[4660]=['',3.2206395,-0.2933446,6.05];
ybs[4661]=['ε Mus',3.2232236,-1.188126,4.11];
ybs[4662]=['',3.2222905,0.9263742,5.81];
ybs[4663]=['',3.2224853,0.503064,5.7];
ybs[4664]=['β Cha',3.2270997,-1.3862445,4.26];
ybs[4665]=['',3.2239096,-0.6319428,6.15];
ybs[4666]=['',3.2235324,0.26233,6.34];
ybs[4667]=['',3.22539,-0.0710034,6.99];
ybs[4668]=['',3.2254264,-0.0709016,6.54];
ybs[4669]=['ζ Cru',3.226954,-1.1190492,4.04];
ybs[4670]=['',3.2269188,0.5259625,6.23];
ybs[4671]=['13 Vir',3.2276501,-0.0157245,5.9];
ybs[4672]=['',3.2293026,-0.9644126,5];
ybs[4673]=['',3.2173799,1.5036947,6.33];
ybs[4674]=['',3.2291503,0.4519367,6.48];
ybs[4675]=['8 Com',3.2303962,0.4000472,6.27];
ybs[4676]=['',3.209835,1.5274349,6.28];
ybs[4677]=['',3.2277164,1.3098143,5.38];
ybs[4678]=['9 Com',3.2311401,0.4894471,6.33];
ybs[4679]=['η Vir',3.2330386,-0.0136244,3.89];
ybs[4680]=['3 CVn',3.2324195,0.8529509,5.29];
ybs[4681]=['',3.2342951,-0.3890203,5.97];
ybs[4682]=['',3.2358879,-1.1511568,6.21];
ybs[4683]=['',3.2347857,0.4626133,5.54];
ybs[4684]=['',3.2346429,0.4518358,6.15];
ybs[4685]=['16 Vir',3.2349604,0.0558303,4.96];
ybs[4686]=['ζ Crv',3.2359692,-0.3897229,5.21];
ybs[4687]=['11 Com',3.2365117,0.3085592,4.74];
ybs[4688]=['',3.2363532,0.4702106,7.13];
ybs[4689]=['',3.2375393,-0.2387468,5.14];
ybs[4690]=['ε Cru',3.2397133,-1.0561811,3.59];
ybs[4691]=['70 UMa',3.2368437,1.0079321,5.55];
ybs[4692]=['',3.2422829,-0.9859068,5.92];
ybs[4693]=['ζ2 Mus',3.2431785,-1.1804624,5.15];
ybs[4694]=['ζ1 Mus',3.243534,-1.1941729,5.74];
ybs[4695]=['',3.242865,0.4304037,6.19];
ybs[4696]=['',3.2460768,-1.0086196,5.39];
ybs[4697]=['12 Com',3.2442771,0.4491178,4.81];
ybs[4698]=['17 Vir',3.244479,0.0906175,6.4];
ybs[4699]=['',3.2614601,-1.5021926,6.33];
ybs[4700]=['',3.2480405,-1.1823765,6.36];
ybs[4701]=['6 Crv',3.2482006,-0.4355306,5.68];
ybs[4702]=['',3.249257,-0.6200505,5.32];
ybs[4703]=['',3.2493889,-0.6879486,6.4];
ybs[4704]=['',3.2499692,-0.6811078,5.79];
ybs[4705]=['4 CVn',3.2497585,0.7405307,6.06];
ybs[4706]=['5 CVn',3.2507385,0.89795,4.8];
ybs[4707]=['13 Com',3.252137,0.4535264,5.18];
ybs[4708]=['',3.2543412,-0.7242698,6.25];
ybs[4709]=['',3.2527352,0.4445236,6.42];
ybs[4710]=['',3.2570181,-1.149892,6.3];
ybs[4711]=['',3.256083,-0.7439964,6.11];
ybs[4712]=['',3.2561571,-0.204617,5.95];
ybs[4713]=['',3.2567165,-0.4862936,6.09];
ybs[4714]=['',3.2570043,-0.6160975,5.73];
ybs[4715]=['',3.2562594,0.41561,6.03];
ybs[4716]=['71 UMa',3.2551553,0.9889747,5.81];
ybs[4717]=['',3.2552758,1.111589,6.32];
ybs[4718]=['6 CVn',3.2587769,0.6790244,5.02];
ybs[4719]=['',3.2623353,-1.1036735,4.86];
ybs[4720]=['α1 Cru',3.2626999,-1.1032661,1.33];
ybs[4721]=['α2 Cru',3.2627436,-1.103271,1.73];
ybs[4722]=['',3.2622169,-0.8999644,4.82];
ybs[4723]=['14 Com',3.2612561,0.473944,4.95];
ybs[4724]=['',3.2634014,-0.8556765,6.26];
ybs[4725]=['',3.2635369,-0.5749693,5.55];
ybs[4726]=['',3.2662699,-1.1153081,6];
ybs[4727]=['γ Com',3.2635901,0.4913978,4.36];
ybs[4728]=['16 Com',3.2638157,0.4662167,5];
ybs[4729]=['',3.2664795,-1.0315807,5.5];
ybs[4730]=['',3.2606722,1.2534322,6.24];
ybs[4731]=['',3.2670116,0.1483009,6.37];
ybs[4732]=['',3.2676529,-0.2922589,6.35];
ybs[4733]=['σ Cen',3.2688231,-0.878665,3.91];
ybs[4734]=['',3.2702474,-1.1249452,6.04];
ybs[4735]=['73 UMa',3.2661832,0.9703944,5.7];
ybs[4736]=['',3.2677594,-0.0825285,6.22];
ybs[4737]=['',3.2706744,-1.0805071,6.22];
ybs[4738]=['',3.2701901,-0.6833769,5.44];
ybs[4739]=['',3.27116,-0.9864773,6.15];
ybs[4740]=['',3.2709967,0.4557659,6.54];
ybs[4741]=['',3.2714707,0.45005,6.65];
ybs[4742]=['17 Com',3.2722044,0.4502878,5.29];
ybs[4743]=['18 Com',3.2745537,0.4188046,5.48];
ybs[4744]=['',3.2770333,-0.9885169,5.8];
ybs[4745]=['',3.2771604,-0.7304069,6.02];
ybs[4746]=['20 Com',3.2757553,0.3627313,5.69];
ybs[4747]=['δ Crv',3.2765661,-0.2902252,2.95];
ybs[4748]=['',3.2774893,-0.2357271,6.35];
ybs[4749]=['',3.278467,-0.4155588,5.63];
ybs[4750]=['74 UMa',3.2764526,1.0173997,5.35];
ybs[4751]=['7 CVn',3.2769552,0.8974908,6.21];
ybs[4752]=['75 UMa',3.276954,1.0237121,6.08];
ybs[4753]=['γ Cru',3.2825862,-0.9987886,1.63];
ybs[4754]=['γ Cru',3.2830817,-0.998226,6.42];
ybs[4755]=['4 Dra',3.276869,1.2058129,4.95];
ybs[4756]=['21 Com',3.2813592,0.4268058,5.46];
ybs[4757]=['',3.2803598,0.9243892,6.21];
ybs[4758]=['',3.2848448,-1.0391147,5.48];
ybs[4759]=['',3.2874633,-1.276091,5.88];
ybs[4760]=['',3.2829597,0.130745,6.05];
ybs[4761]=['',3.2860646,-1.1103626,5.95];
ybs[4762]=['',3.2842794,-0.0901551,6.19];
ybs[4763]=['γ Mus',3.2887319,-1.2609305,3.87];
ybs[4764]=['',3.2863113,-0.5697904,6.46];
ybs[4765]=['η Crv',3.2861912,-0.2846473,4.31];
ybs[4766]=['',3.2884925,-0.2438593,5.74];
ybs[4767]=['20 Vir',3.2903251,0.1777207,6.26];
ybs[4768]=['',3.2919017,-0.3474048,6.26];
ybs[4769]=['',3.2927287,-0.2259005,5.58];
ybs[4770]=['22 Com',3.292521,0.4218494,6.29];
ybs[4771]=['21 Vir',3.2936199,-0.1669373,5.48];
ybs[4772]=['',3.2948281,-0.8730534,6.38];
ybs[4773]=['',3.2928022,0.5783084,5.42];
ybs[4774]=['',3.2934185,0.5807036,6.24];
ybs[4775]=['β CVn',3.2931424,0.7198548,4.26];
ybs[4776]=['β Crv',3.2963484,-0.4103177,2.65];
ybs[4777]=['κ Dra',3.2914804,1.216066,3.87];
ybs[4778]=['',3.2979087,-0.7816652,5.77];
ybs[4779]=['23 Com',3.2981189,0.3929852,4.81];
ybs[4780]=['',3.3015871,-1.0813128,6.22];
ybs[4781]=['24 Com',3.2992494,0.3187751,6.56];
ybs[4782]=['24 Com',3.2993585,0.3187703,5.02];
ybs[4783]=['',3.2993589,0.3799344,5.85];
ybs[4784]=['',3.3024762,-0.717935,5.13];
ybs[4785]=['6 Dra',3.296895,1.220145,4.94];
ybs[4786]=['',3.303609,-0.6978294,5.8];
ybs[4787]=['',3.3032736,-0.3602343,6.2];
ybs[4788]=['α Mus',3.3092808,-1.208608,2.69];
ybs[4789]=['25 Vir',3.306737,-0.1037522,5.87];
ybs[4790]=['',3.3044112,1.0362768,5.5];
ybs[4791]=['25 Com',3.3073971,0.2963018,5.68];
ybs[4792]=['τ Cen',3.3110645,-0.8491664,3.86];
ybs[4793]=['',3.310859,-0.4756273,5.45];
ybs[4794]=['',3.3187288,-1.3174068,6.49];
ybs[4795]=['',3.3122845,0.0553266,6.33];
ybs[4796]=['',3.3166157,-1.1747025,6.25];
ybs[4797]=['',3.3136018,0.0304077,5.71];
ybs[4798]=['',3.314124,0.1200063,7.08];
ybs[4799]=['',3.3153391,-0.3204901,6];
ybs[4800]=['',3.3168028,-0.5329302,5.89];
ybs[4801]=['9 CVn',3.3150571,0.7114308,6.37];
ybs[4802]=['',3.3163571,0.3935195,6.38];
ybs[4803]=['χ Vir',3.3174716,-0.1415108,4.66];
ybs[4804]=['',3.3211989,-1.1628084,6.26];
ybs[4805]=['26 Com',3.3167459,0.3656478,5.46];
ybs[4806]=['',3.3173268,0.6255178,6.45];
ybs[4807]=['',3.3204585,-0.6998745,4.64];
ybs[4808]=['',3.3271241,-0.8073506,5.84];
ybs[4809]=['γ Cen',3.3277464,-0.8564668,2.17];
ybs[4810]=['',3.330793,-1.2133468,6.33];
ybs[4811]=['',3.3263208,-0.2290893,6.08];
ybs[4812]=['',3.3263353,-0.2291135,5.98];
ybs[4813]=['',3.3298196,-1.043672,4.93];
ybs[4814]=['27 Vir',3.3275084,0.1800164,6.19];
ybs[4815]=['γ Vir',3.3279627,-0.0272559,3.65];
ybs[4816]=['γ Vir',3.3279627,-0.0272559,3.68];
ybs[4817]=['',3.3287831,-0.3468108,6.03];
ybs[4818]=['ρ Vir',3.328869,0.1766862,4.88];
ybs[4819]=['31 Vir',3.3291823,0.1168409,5.59];
ybs[4820]=['',3.3338349,-1.1025366,5.31];
ybs[4821]=['',3.3324381,-0.8539052,4.66];
ybs[4822]=['',3.3335998,-0.9784195,6.08];
ybs[4823]=['76 UMa',3.3268365,1.0925907,6.07];
ybs[4824]=['',3.3350269,-0.9824139,6];
ybs[4825]=['',3.3364832,-1.0300075,6.4];
ybs[4826]=['',3.3360341,-0.7031898,6.44];
ybs[4827]=['',3.3365816,-0.029478,5.93];
ybs[4828]=['',3.3383524,-0.636367,6.39];
ybs[4829]=['',3.3384075,-0.4962995,5.48];
ybs[4830]=['',3.3334244,1.0654097,6.38];
ybs[4831]=['',3.3436715,-1.203282,6.16];
ybs[4832]=['ι Cru',3.3460041,-1.0662727,4.69];
ybs[4833]=['',3.3397587,0.7677897,6.33];
ybs[4834]=['β Mus',3.3491296,-1.1906601,3.05];
ybs[4835]=['10 CVn',3.3421697,0.6835931,5.95];
ybs[4836]=['',3.3426947,0.7911299,4.99];
ybs[4837]=['32 Vir',3.3451714,0.1319732,5.22];
ybs[4838]=['',3.3491543,-0.9878673,4.65];
ybs[4839]=['33 Vir',3.3484586,0.164554,5.67];
ybs[4840]=['',3.3505041,-0.5834157,5.86];
ybs[4841]=['27 Com',3.3495836,0.2873821,5.12];
ybs[4842]=['',3.3377313,1.4051496,6.4];
ybs[4843]=['β Cru',3.3551096,-1.0437106,1.25];
ybs[4844]=['',3.3513782,0.1019124,6.34];
ybs[4845]=['34 Vir',3.3521555,0.2067585,6.07];
ybs[4846]=['',3.3537285,-0.1119379,6.26];
ybs[4847]=['',3.3553525,-0.435691,6.44];
ybs[4848]=['35 Vir',3.3549656,0.060409,6.41];
ybs[4849]=['',3.3518286,1.0937833,5.89];
ybs[4850]=['',3.3577579,-0.4836138,5.66];
ybs[4851]=['28 Com',3.3565556,0.2345984,6.56];
ybs[4852]=['',3.3645865,-1.2583434,5.55];
ybs[4853]=['7 Dra',3.3527919,1.1637619,5.43];
ybs[4854]=['',3.3588383,0.4315986,6.31];
ybs[4855]=['29 Com',3.3594514,0.2445383,5.7];
ybs[4856]=['11 CVn',3.3581785,0.8439615,6.27];
ybs[4857]=['',3.3577405,1.0508362,5.85];
ybs[4858]=['',3.3659888,-1.0561366,6.75];
ybs[4859]=['30 Com',3.3610219,0.4789319,5.78];
ybs[4860]=['ι Oct',3.3917871,-1.4822297,5.46];
ybs[4861]=['',3.3662649,-0.8477247,6.24];
ybs[4862]=['',3.3691402,-0.9232574,5.73];
ybs[4863]=['',3.3654245,0.3970972,6.43];
ybs[4864]=['',3.367638,-0.5953445,4.91];
ybs[4865]=['',3.3647867,0.6528507,5.89];
ybs[4866]=['',3.370788,-1.0548934,5.72];
ybs[4867]=['',3.370454,-0.182379,6.41];
ybs[4868]=['37 Vir',3.3713665,0.0514083,6.02];
ybs[4869]=['',3.3732133,-0.6945011,5.98];
ybs[4870]=['',3.373962,-0.8413412,6.33];
ybs[4871]=['',3.3731456,-0.468607,6.15];
ybs[4872]=['',3.3754847,-0.9414399,6.24];
ybs[4873]=['31 Com',3.3715197,0.4787329,4.94];
ybs[4874]=['32 Com',3.3738253,0.2960561,6.32];
ybs[4875]=['',3.3783771,-0.9610397,5.93];
ybs[4876]=['',3.3749456,0.2794517,6.3];
ybs[4877]=['',3.3798363,-1.0548698,5.76];
ybs[4878]=['',3.378467,-0.8561598,4.33];
ybs[4879]=['',3.379732,-0.7031908,4.27];
ybs[4880]=['κ Cru',3.38183,-1.0557125,5.9];
ybs[4881]=['38 Vir',3.3782769,-0.06395,6.11];
ybs[4882]=['',3.3568122,1.4539735,5.85];
ybs[4883]=['',3.3573157,1.4538816,5.28];
ybs[4884]=['35 Com',3.378548,0.3688579,4.9];
ybs[4885]=['',3.3841684,-1.0217405,6.58];
ybs[4886]=['',3.3802388,-0.0756574,6.44];
ybs[4887]=['λ Cru',3.3854454,-1.0342384,4.62];
ybs[4888]=['μ1 Cru',3.3851237,-0.9998749,4.03];
ybs[4889]=['μ2 Cru',3.385211,-0.99971,5.17];
ybs[4890]=['41 Vir',3.380945,0.2148095,6.25];
ybs[4891]=['',3.3832543,-0.2052418,6];
ybs[4892]=['ψ Vir',3.3834181,-0.1684201,4.79];
ybs[4893]=['',3.3865158,-0.7725305,5.89];
ybs[4894]=['',3.3824377,0.5833511,6.26];
ybs[4895]=['ε UMa',3.3812273,0.9747455,1.77];
ybs[4896]=['',3.3880188,-0.7509556,5.47];
ybs[4897]=['',3.3943508,-1.261801,5.93];
ybs[4898]=['',3.3910513,-0.9939089,5.32];
ybs[4899]=['',3.3853914,0.8218032,5.84];
ybs[4900]=['δ Vir',3.3887642,0.0573651,3.38];
ybs[4901]=['',3.3901695,-0.2694375,6.17];
ybs[4902]=['',3.392957,-0.4637495,6.62];
ybs[4903]=['',3.3958245,-0.8955135,5.16];
ybs[4904]=['α1 CVn',3.390169,0.6667864,5.6];
ybs[4905]=['α2 CVn',3.3902634,0.6668494,2.9];
ybs[4906]=['8 Dra',3.3871829,1.1401862,5.24];
ybs[4907]=['',3.3911302,0.9422822,5.82];
ybs[4908]=['',3.3975109,-0.3990586,6.31];
ybs[4909]=['',3.3949313,0.8040103,6.12];
ybs[4910]=['36 Com',3.403119,0.3019268,4.78];
ybs[4911]=['44 Vir',3.4065206,-0.0684546,5.79];
ybs[4912]=['',3.4106954,-0.5866989,6.02];
ybs[4913]=['δ Mus',3.4195065,-1.2506809,3.62];
ybs[4914]=['37 Com',3.4088661,0.5353773,4.9];
ybs[4915]=['46 Vir',3.4106186,-0.0607149,5.99];
ybs[4916]=['',3.4106243,0.3187489,6.2];
ybs[4917]=['',3.4007372,1.3153175,6.01];
ybs[4918]=['9 Dra',3.4064261,1.1604174,5.32];
ybs[4919]=['38 Com',3.4128754,0.2969334,5.96];
ybs[4920]=['',3.4230835,-1.2494087,6.03];
ybs[4921]=['78 UMa',3.4103737,0.9818577,4.93];
ybs[4922]=['ε Vir',3.4173703,0.1893556,2.83];
ybs[4923]=['ξ1 Cen',3.4241271,-0.8663275,4.85];
ybs[4924]=['',3.4146699,1.1082897,6];
ybs[4925]=['',3.4246267,-0.3611561,5.58];
ybs[4926]=['',3.418716,1.0403259,6.53];
ybs[4927]=['48 Vir',3.4250585,-0.0658509,6.59];
ybs[4928]=['',3.4294194,-0.7209289,6.26];
ybs[4929]=['',3.4327628,-0.9114879,6.43];
ybs[4930]=['',3.436019,-0.8477573,4.71];
ybs[4931]=['',3.4372206,-0.7277652,5.59];
ybs[4932]=['ξ2 Cen',3.4388153,-0.8729321,4.27];
ybs[4933]=['14 CVn',3.4326207,0.6228991,5.25];
ybs[4934]=['',3.4412855,-1.0466686,5.99];
ybs[4935]=['',3.4330237,0.7881772,5.63];
ybs[4936]=['39 Com',3.4354771,0.3672875,5.99];
ybs[4937]=['',3.4385239,-0.6278152,6.54];
ybs[4938]=['',3.4345818,0.5047511,6.54];
ybs[4939]=['40 Com',3.4355619,0.3928179,5.6];
ybs[4940]=['',3.4272203,1.2726195,6.31];
ybs[4941]=['',3.4420954,-0.9349524,5.71];
ybs[4942]=['θ Mus',3.4446762,-1.1417144,5.51];
ybs[4943]=['',3.4347316,1.0809282,6.14];
ybs[4944]=['41 Com',3.4389952,0.4802367,4.8];
ybs[4945]=['49 Vir',3.4425439,-0.189357,5.19];
ybs[4946]=['',3.4421128,0.4790361,6.19];
ybs[4947]=['',3.4453406,-0.1587102,5.55];
ybs[4948]=['ψ Hya',3.4477413,-0.405387,4.95];
ybs[4949]=['',3.4483945,-0.1683756,6.32];
ybs[4950]=['',3.4480433,0.1730203,5.78];
ybs[4951]=['50 Vir',3.4506509,-0.1821818,5.94];
ybs[4952]=['',3.450541,0.2921647,5.91];
ybs[4953]=['θ Vir',3.4514506,-0.0985703,4.38];
ybs[4954]=['',3.4496026,0.6512561,6.02];
ybs[4955]=['',3.4566579,-0.9193618,6.06];
ybs[4956]=['',3.4614296,-1.2226101,5.91];
ybs[4957]=['15 CVn',3.4498235,0.6706439,6.28];
ybs[4958]=['α Com',3.451369,0.3040481,5.22];
ybs[4959]=['α Com',3.451369,0.3040481,5.22];
ybs[4960]=['',3.4571708,-0.739001,5.79];
ybs[4961]=['17 CVn',3.4513635,0.670034,5.91];
ybs[4962]=['',3.4610705,-1.1067349,6.33];
ybs[4963]=['',3.4582425,-0.7588244,5.25];
ybs[4964]=['',3.4497585,1.0842047,6.54];
ybs[4965]=['',3.4626802,-1.0477077,4.6];
ybs[4966]=['',3.4734683,-1.371048,5.85];
ybs[4967]=['',3.4653217,-1.1577686,5.9];
ybs[4968]=['',3.4591329,-0.4653079,6.5];
ybs[4969]=['',3.4610516,-0.6616805,4.85];
ybs[4970]=['',3.465492,-1.0458879,6.16];
ybs[4971]=['53 Vir',3.4607717,-0.2846119,5.04];
ybs[4972]=['',3.464618,-0.7471412,6.22];
ybs[4973]=['β Com',3.4594572,0.4846704,4.26];
ybs[4974]=['',3.4606682,0.4214902,6.33];
ybs[4975]=['',3.4671824,-0.8867709,5.89];
ybs[4976]=['',3.4626019,0.1998006,5.77];
ybs[4977]=['',3.4627356,0.3253868,6.53];
ybs[4978]=['',3.4709955,-1.0261137,5.89];
ybs[4979]=['',3.4712102,-1.0334343,4.92];
ybs[4980]=['54 Vir',3.4668648,-0.3304763,6.28];
ybs[4981]=['',3.4694789,-0.7548031,6.16];
ybs[4982]=['',3.4653887,0.324957,6.11];
ybs[4983]=['η Mus',3.4761099,-1.1868651,4.8];
ybs[4984]=['',3.477055,-1.2180236,6.37];
ybs[4985]=['55 Vir',3.4700877,-0.3497457,5.33];
ybs[4986]=['',3.4729373,-0.8563403,5.89];
ybs[4987]=['',3.4672933,0.6989097,4.92];
ybs[4988]=['',3.4712115,0.1958888,5.67];
ybs[4989]=['',3.4746087,-0.6366797,6.19];
ybs[4990]=['',3.4824761,-1.1387577,6.07];
ybs[4991]=['57 Vir',3.4779399,-0.3499539,5.22];
ybs[4992]=['',3.484651,-1.1674718,4.87];
ybs[4993]=['',3.4649965,1.2686908,6.59];
ybs[4994]=['19 CVn',3.4751937,0.7111757,5.79];
ybs[4995]=['',3.4796606,-0.0261505,6.68];
ybs[4996]=['',3.4820569,-0.5517646,5.1];
ybs[4997]=['',3.4786042,0.330633,6.45];
ybs[4998]=['',3.4838053,-0.7694643,5.84];
ybs[4999]=['',3.4585118,1.4025977,6.25];
ybs[5000]=['',3.4799029,0.3434377,6.45];
ybs[5001]=['59 Vir',3.4810625,0.162603,5.22];
ybs[5002]=['',3.494364,-1.2591292,6.04];
ybs[5003]=['',3.4831268,0.2368051,5.33];
ybs[5004]=['',3.4843358,-0.0136877,6.37];
ybs[5005]=['σ Vir',3.4847296,0.0935874,4.8];
ybs[5006]=['',3.4898758,-0.8969857,6.19];
ybs[5007]=['20 CVn',3.4839482,0.706246,4.73];
ybs[5008]=['',3.4782186,1.1920647,6.2];
ybs[5009]=['61 Vir',3.4885068,-0.3214689,4.74];
ybs[5010]=['γ Hya',3.4908284,-0.4062952,3];
ybs[5011]=['',3.4901895,0.0624902,6.62];
ybs[5012]=['',3.4880835,0.5932484,5.82];
ybs[5013]=['21 CVn',3.4867822,0.8652378,5.15];
ybs[5014]=['',3.4989702,-1.0451096,6.18];
ybs[5015]=['',3.4907093,0.611227,6.02];
ybs[5016]=['',3.498903,-0.9224952,5.48];
ybs[5017]=['',3.4997774,-0.9757709,6.02];
ybs[5018]=['ι Cen',3.4983627,-0.6426175,2.75];
ybs[5019]=['',3.5001817,-0.8200871,5.77];
ybs[5020]=['',3.5100094,-1.2610576,6.05];
ybs[5021]=['',3.4982333,0.0494736,6.26];
ybs[5022]=['23 CVn',3.4960457,0.6988899,5.6];
ybs[5023]=['',3.5020284,-0.3420109,6.21];
ybs[5024]=['',3.5078442,-1.0660279,6.18];
ybs[5025]=['',3.5080055,-1.066309,4.53];
ybs[5026]=['',3.5060576,-0.9126291,5.83];
ybs[5027]=['',3.5026067,0.0345638,5.69];
ybs[5028]=['',3.5085902,-0.8386254,6.16];
ybs[5029]=['',3.5093291,-0.8494411,6.38];
ybs[5030]=['64 Vir',3.5046123,0.0881032,5.87];
ybs[5031]=['',3.5142473,-1.1282202,4.53];
ybs[5032]=['ι1 Mus',3.5202878,-1.3088918,5.05];
ybs[5033]=['',3.5094321,-0.5811353,6.22];
ybs[5034]=['63 Vir',3.5086409,-0.3114,5.37];
ybs[5035]=['',3.5035718,0.7643887,6.35];
ybs[5036]=['',3.512988,-0.8714345,6.48];
ybs[5037]=['65 Vir',3.5097694,-0.0878079,5.89];
ybs[5038]=['',3.5196025,-1.1273339,5.31];
ybs[5039]=['',3.5227992,-1.2345338,5.67];
ybs[5040]=['66 Vir',3.5151766,-0.0919831,5.75];
ybs[5041]=['ι2 Mus',3.5298573,-1.3054668,6.63];
ybs[5042]=['',3.5117202,0.6445048,6.07];
ybs[5043]=['',3.5147558,0.215122,6.44];
ybs[5044]=['ζ UMa',3.5113504,0.9567683,2.27];
ybs[5045]=['ζ UMa',3.5114159,0.9567054,3.95];
ybs[5046]=['α Vir',3.5180475,-0.1966572,0.98];
ybs[5047]=['',3.5172263,0.414484,5.78];
ybs[5048]=['',3.5226035,-0.6957115,5.09];
ybs[5049]=['',3.5222679,-0.0226641,5.97];
ybs[5050]=['',3.5261606,-0.7261261,5.69];
ybs[5051]=['',3.5271062,-0.8595705,6.31];
ybs[5052]=['80 UMa',3.5170055,0.9578681,4.01];
ybs[5053]=['',3.5281702,-0.8637051,6.28];
ybs[5054]=['68 Vir',3.5247321,-0.2236419,5.25];
ybs[5055]=['',3.5274851,-0.702825,6.4];
ybs[5056]=['',3.5355723,-1.2170807,6.2];
ybs[5057]=['',3.5218687,0.8014901,5.88];
ybs[5058]=['69 Vir',3.5279786,-0.280639,4.76];
ybs[5059]=['',3.5366775,-1.130647,6.11];
ybs[5060]=['',3.5199642,1.1022625,6.5];
ybs[5061]=['',3.5372947,-0.8948427,5.06];
ybs[5062]=['70 Vir',3.5318416,0.2386432,4.98];
ybs[5063]=['',3.519709,1.261616,5.79];
ybs[5064]=['',3.5245532,1.1279998,6.66];
ybs[5065]=['',3.5249963,1.127719,7.04];
ybs[5066]=['',3.5291296,0.9187429,6.34];
ybs[5067]=['',3.5314074,0.7090238,6.47];
ybs[5068]=['',3.5356153,-0.025655,6.43];
ybs[5069]=['',3.530088,0.8810687,6.8];
ybs[5070]=['',3.5379551,-0.4081763,4.97];
ybs[5071]=['71 Vir',3.5353137,0.1869743,5.65];
ybs[5072]=['',3.5566681,-1.3556489,6.48];
ybs[5073]=['',3.5325508,0.8833589,6.43];
ybs[5074]=['κ Oct',3.5984161,-1.4949926,5.58];
ybs[5075]=['',3.530822,1.044408,5.4];
ybs[5076]=['',3.5387876,0.1234566,6.17];
ybs[5077]=['',3.5386217,0.1031138,6.51];
ybs[5078]=['72 Vir',3.5408332,-0.1147648,6.09];
ybs[5079]=['',3.5440742,-0.6896254,3.88];
ybs[5080]=['',3.5460784,-0.4924936,6.47];
ybs[5081]=['',3.5218865,1.3707448,5.77];
ybs[5082]=['',3.5486106,-0.6720231,6.16];
ybs[5083]=['',3.5563259,-1.1473288,6.37];
ybs[5084]=['73 Vir',3.5480728,-0.3287123,6.01];
ybs[5085]=['74 Vir',3.5475366,-0.1110168,4.69];
ybs[5086]=['',3.5436885,0.7330557,6.08];
ybs[5087]=['',3.5506456,-0.502613,5.69];
ybs[5088]=['',3.5505594,-0.517841,6.45];
ybs[5089]=['75 Vir',3.5515803,-0.2699647,5.55];
ybs[5090]=['76 Vir',3.5519704,-0.1792411,5.21];
ybs[5091]=['',3.5521176,-0.1274047,6.68];
ybs[5092]=['',3.5507448,0.4230954,6.11];
ybs[5093]=['',3.5593333,-0.8443321,6.33];
ybs[5094]=['',3.56003,-0.5832058,6.44];
ybs[5095]=['78 Vir',3.5568535,0.0620352,4.94];
ybs[5096]=['',3.5594643,-0.2324581,5.91];
ybs[5097]=['ζ Vir',3.5593638,-0.0122218,3.37];
ybs[5098]=['',3.5572725,0.6751748,6.37];
ybs[5099]=['81 UMa',3.5557075,0.9641906,5.6];
ybs[5100]=['',3.5591966,0.6471346,4.98];
ybs[5101]=['80 Vir',3.5630464,-0.0959994,5.73];
ybs[5102]=['24 CVn',3.5573942,0.8536688,4.7];
ybs[5103]=['',3.5718422,-1.0785405,5.63];
ybs[5104]=['',3.5629693,0.1762865,6.49];
ybs[5105]=['',3.5823259,-1.322738,6.34];
ybs[5106]=['',3.5609343,0.7695613,6.84];
ybs[5107]=['',3.5692936,-0.6033908,6.5];
ybs[5108]=['',3.5706567,-0.7722601,5.98];
ybs[5109]=['',3.5794731,-1.2313041,6.1];
ybs[5110]=['',3.5689861,-0.4642398,5.78];
ybs[5111]=['',3.572013,-0.8121398,5.9];
ybs[5112]=['',3.5756921,-1.0213438,6.42];
ybs[5113]=['',3.5689783,0.4277692,5.74];
ybs[5114]=['',3.5786775,-1.0075192,6.01];
ybs[5115]=['',3.5849943,-1.2372918,6.59];
ybs[5116]=['',3.566994,0.8618894,6.49];
ybs[5117]=['25 CVn',3.5708215,0.6316544,4.82];
ybs[5118]=['',3.5773123,-0.5177418,5.83];
ybs[5119]=['',3.5741462,0.2478009,6.52];
ybs[5120]=['',3.5850923,-1.1288822,5.79];
ybs[5121]=['',3.5561065,1.3341677,6.57];
ybs[5122]=['ε Cen',3.5831591,-0.9349678,2.3];
ybs[5123]=['',3.5715373,0.8833268,6.48];
ybs[5124]=['',3.5834968,-0.8735998,6];
ybs[5125]=['',3.5818164,-0.6955387,6.27];
ybs[5126]=['',3.5823921,-0.7008421,5.6];
ybs[5127]=['',3.5780376,0.3169822,6.48];
ybs[5128]=['',3.5805031,0.18575,5.57];
ybs[5129]=['',3.5678074,1.2415967,5.5];
ybs[5130]=['',3.5927433,-1.0278258,5.38];
ybs[5131]=['',3.5913352,-0.954048,5.01];
ybs[5132]=['82 UMa',3.579246,0.9218468,5.46];
ybs[5133]=['',3.5831354,0.5394579,6.21];
ybs[5134]=['1 Boo',3.5851465,0.3464891,5.75];
ybs[5135]=['',3.5849003,0.4880304,6.23];
ybs[5136]=['',3.589493,-0.4110725,6.59];
ybs[5137]=['',3.5907688,-0.5881739,6.05];
ybs[5138]=['',3.5831689,0.8799282,6.32];
ybs[5139]=['2 Boo',3.5866874,0.3908267,5.62];
ybs[5140]=['82 Vir',3.5896812,-0.1536943,5.01];
ybs[5141]=['',3.5966465,-0.9925813,6];
ybs[5142]=['',3.5962746,-0.8882544,6.41];
ybs[5143]=['',3.5827454,0.9966566,6.29];
ybs[5144]=['83 UMa',3.5845327,0.952574,4.66];
ybs[5145]=['',3.5960038,-0.7243778,5.98];
ybs[5146]=['',3.5920331,0.1446088,6.16];
ybs[5147]=['',3.5992998,-0.7360057,5.98];
ybs[5148]=['',3.6022086,-0.8921326,6.47];
ybs[5149]=['84 Vir',3.5958154,0.0599588,5.36];
ybs[5150]=['',3.5925621,0.725557,6.3];
ybs[5151]=['',3.5937866,0.6088778,5.98];
ybs[5152]=['',3.5872302,1.1295675,5.85];
ybs[5153]=['',3.5996333,-0.0977624,6.51];
ybs[5154]=['',3.5985213,0.3944052,6.13];
ybs[5155]=['83 Vir',3.6023864,-0.284166,5.6];
ybs[5156]=['',3.6037086,-0.4468587,6.21];
ybs[5157]=['',3.6074469,-0.457594,5.81];
ybs[5158]=['1 Cen',3.6079028,-0.5785063,4.23];
ybs[5159]=['',3.598456,0.9069068,6.02];
ybs[5160]=['85 Vir',3.6071307,-0.2769769,6.19];
ybs[5161]=['',3.6155771,-1.0941768,6.51];
ybs[5162]=['',3.6126768,-0.8994488,4.65];
ybs[5163]=['86 Vir',3.6086178,-0.218667,5.51];
ybs[5164]=['',3.6134526,-0.6344924,5.15];
ybs[5165]=['',3.6161501,-0.8787927,5.91];
ybs[5166]=['',3.61694,-0.8800428,5.45];
ybs[5167]=['',3.6039814,0.9734961,6.5];
ybs[5168]=['',3.6141897,-0.1712326,6.05];
ybs[5169]=['',3.608893,0.7153514,5.87];
ybs[5170]=['',3.6093606,0.6702399,5.94];
ybs[5171]=['87 Vir',3.6151921,-0.3134907,5.43];
ybs[5172]=['3 Boo',3.611407,0.4468104,5.95];
ybs[5173]=['',3.6127424,0.1090612,6.33];
ybs[5174]=['',3.5900287,1.3606864,5.91];
ybs[5175]=['τ Boo',3.6139115,0.3029005,4.5];
ybs[5176]=['',3.6123225,0.6709214,5.5];
ybs[5177]=['84 UMa',3.6100343,0.9482524,5.7];
ybs[5178]=['',3.6583908,-1.4445336,5.95];
ybs[5179]=['',3.6220825,-0.6249192,6.53];
ybs[5180]=['ν Cen',3.6248043,-0.7293554,3.41];
ybs[5181]=['η UMa',3.6143594,0.858905,1.86];
ybs[5182]=['2 Cen',3.6243534,-0.6030472,4.19];
ybs[5183]=['μ Cen',3.6253153,-0.7430751,3.04];
ybs[5184]=['',3.63646,-1.213039,5.75];
ybs[5185]=['',3.6196788,0.5426026,5.62];
ybs[5186]=['89 Vir',3.6258859,-0.318266,4.97];
ybs[5187]=['',3.6271309,-0.5093301,6.18];
ybs[5188]=['',3.6283299,-0.6981689,6.44];
ybs[5189]=['',3.6208216,0.6883824,7.4];
ybs[5190]=['υ Boo',3.6235998,0.2739562,4.07];
ybs[5191]=['6 Boo',3.6245342,0.3693636,4.91];
ybs[5192]=['',3.6289885,-0.3490343,6.53];
ybs[5193]=['',3.5862121,1.4425062,5.98];
ybs[5194]=['',3.6243681,0.6375965,6.38];
ybs[5195]=['',3.6278505,0.0941815,6.01];
ybs[5196]=['',3.634943,-0.8203019,5.77];
ybs[5197]=['',3.6364607,-0.9234932,5.25];
ybs[5198]=['',3.6338603,-0.6376395,6.35];
ybs[5199]=['',3.6324143,-0.4274594,6.45];
ybs[5200]=['3 Cen',3.6347234,-0.5776187,4.56];
ybs[5201]=['3 Cen',3.6347598,-0.5776235,6.06];
ybs[5202]=['',3.6355145,-0.5536197,6.12];
ybs[5203]=['',3.6233972,1.0714219,5.96];
ybs[5204]=['',3.6301877,0.6051339,6.65];
ybs[5205]=['',3.6305317,0.6032483,5.87];
ybs[5206]=['',3.6266465,1.0199424,6.46];
ybs[5207]=['',3.6436519,-0.9332941,5.89];
ybs[5208]=['',3.6494997,-1.1825027,5.71];
ybs[5209]=['',3.6333188,0.5994063,4.74];
ybs[5210]=['',3.6360053,0.2105687,6.04];
ybs[5211]=['4 Cen',3.6407341,-0.5589962,4.73];
ybs[5212]=['',3.6423021,-0.624207,5.54];
ybs[5213]=['',3.644401,-0.8242878,6.1];
ybs[5214]=['',3.6437059,-0.6181019,6.19];
ybs[5215]=['7 Boo',3.6398648,0.3112342,5.7];
ybs[5216]=['10 Dra',3.6304179,1.1278754,4.65];
ybs[5217]=['',3.6281165,1.1905647,6.4];
ybs[5218]=['',3.6453254,-0.5003827,6.04];
ybs[5219]=['',3.639465,0.4982509,5.9];
ybs[5220]=['',3.6500907,-0.9121258,5.71];
ybs[5221]=['ζ Cen',3.6513649,-0.8270785,2.55];
ybs[5222]=['90 Vir',3.6466836,-0.0279787,5.15];
ybs[5223]=['',3.6479654,-0.1423983,6.19];
ybs[5224]=['',3.6551112,-0.9465186,6.14];
ybs[5225]=['η Boo',3.6462667,0.3193561,2.68];
ybs[5226]=['',3.6561023,-0.9565097,6];
ybs[5227]=['',3.6517849,-0.5477671,6.51];
ybs[5228]=['86 UMa',3.6416462,0.9359916,5.7];
ybs[5229]=['',3.6547711,-0.8149306,5.83];
ybs[5230]=['',3.6771797,-1.3733722,6.09];
ybs[5231]=['',3.6615087,-1.113274,4.71];
ybs[5232]=['',3.6655277,-1.1501644,6.2];
ybs[5233]=['',3.6513538,0.2435895,6.16];
ybs[5234]=['92 Vir',3.6543345,0.0165977,5.91];
ybs[5235]=['',3.6524801,0.5573331,6.32];
ybs[5236]=['',3.6591183,-0.4035568,6.14];
ybs[5237]=['9 Boo',3.6543041,0.4780871,5.01];
ybs[5238]=['φ Cen',3.6631286,-0.736528,3.83];
ybs[5239]=['υ1 Cen',3.665002,-0.7836985,3.87];
ybs[5240]=['47 Hya',3.6637725,-0.4375765,5.15];
ybs[5241]=['',3.6678749,-0.8808476,5.91];
ybs[5242]=['',3.6728848,-1.0747733,6.49];
ybs[5243]=['',3.6758711,-1.1583233,5.97];
ybs[5244]=['',3.6636931,0.2539523,6];
ybs[5245]=['10 Boo',3.6634855,0.3769398,5.76];
ybs[5246]=['',3.65722,1.0715171,6.37];
ybs[5247]=['48 Hya',3.6702509,-0.4382342,5.77];
ybs[5248]=['',3.6690606,-0.0636778,6.4];
ybs[5249]=['',3.6763771,-0.7037268,6.13];
ybs[5250]=['υ2 Cen',3.6783324,-0.7976479,4.34];
ybs[5251]=['θ Aps',3.6974972,-1.3420514,5.5];
ybs[5252]=['',3.6754656,0.1535252,5.99];
ybs[5253]=['11 Boo',3.6743815,0.4762696,6.23];
ybs[5254]=['τ Vir',3.6769344,0.02524,4.26];
ybs[5255]=['',3.6806901,-0.4804558,5.48];
ybs[5256]=['',3.6863088,-0.9828194,5.92];
ybs[5257]=['β Cen',3.6882732,-1.0554135,0.61];
ybs[5258]=['',3.6836222,-0.5546973,6.18];
ybs[5259]=['',3.6857731,-0.7246806,6.11];
ybs[5260]=['',3.6806421,0.1673477,6.2];
ybs[5261]=['',3.6783489,0.7968376,6.27];
ybs[5262]=['',3.6871531,-0.3930373,6.3];
ybs[5263]=['',3.6850201,0.1865556,6.3];
ybs[5264]=['',3.6854069,0.1300024,6.26];
ybs[5265]=['',3.6868364,0.0838303,6.24];
ybs[5266]=['',3.688389,-0.0956269,6.39];
ybs[5267]=['',3.6894699,-0.2630078,6.28];
ybs[5268]=['',3.6964499,-0.9558581,6.17];
ybs[5269]=['',3.7106039,-1.3080664,6.02];
ybs[5270]=['',3.6815863,0.8879181,6.15];
ybs[5271]=['',3.6995814,-1.0439263,6.42];
ybs[5272]=['',3.6752339,1.1969519,6.34];
ybs[5273]=['',3.6899156,0.0383966,6.28];
ybs[5274]=['',3.6929162,-0.2868134,6.56];
ybs[5275]=['χ Cen',3.6970759,-0.720417,4.36];
ybs[5276]=['',3.6977295,-0.753791,6.2];
ybs[5277]=['π Hya',3.6981053,-0.4673915,3.27];
ybs[5278]=['θ Cen',3.6997124,-0.6364687,2.06];
ybs[5279]=['',3.7078318,-1.1048733,6.4];
ybs[5280]=['95 Vir',3.6992415,-0.164241,5.46];
ybs[5281]=['α Dra',3.6867789,1.1218659,3.65];
ybs[5282]=['',3.710577,-1.0362545,6.34];
ybs[5283]=['',3.7186826,-1.2287369,6.05];
ybs[5284]=['',3.7094696,-0.7603963,6.17];
ybs[5285]=['',3.720852,-1.2185098,6.06];
ybs[5286]=['',3.7129415,-0.9006057,6];
ybs[5287]=['',3.7144826,-0.9343665,4.75];
ybs[5288]=['96 Vir',3.7092849,-0.1820521,6.47];
ybs[5289]=['',3.7033622,0.7637167,5.27];
ybs[5290]=['13 Boo',3.7047028,0.8615197,5.25];
ybs[5291]=['',3.7173985,-0.286196,4.91];
ybs[5292]=['',3.7062371,1.0339554,6.46];
ybs[5293]=['η Aps',3.7565588,-1.4154862,4.91];
ybs[5294]=['12 Boo',3.7146441,0.4362563,4.83];
ybs[5295]=['3 UMi',3.6962653,1.3002104,6.45];
ybs[5296]=['',3.7487619,-1.357132,6.47];
ybs[5297]=['',3.7200164,0.0221049,6.43];
ybs[5298]=['',3.7291889,-0.9383067,5.56];
ybs[5299]=['',3.7244131,-0.4269008,6.34];
ybs[5300]=['',3.7181875,0.5619918,6.11];
ybs[5301]=['',3.7309556,-0.9550599,6.11];
ybs[5302]=['50 Hya',3.7260462,-0.4774602,5.08];
ybs[5303]=['',3.7232398,0.0403858,5.01];
ybs[5304]=['',3.7280103,-0.4661329,6.24];
ybs[5305]=['κ Vir',3.7262477,-0.180972,4.19];
ybs[5306]=['',3.7367279,-0.9979887,5.07];
ybs[5307]=['',3.7294833,-0.0164177,5.91];
ybs[5308]=['',3.734955,-0.7318566,5.61];
ybs[5309]=['',3.7382764,-0.9355719,6.39];
ybs[5310]=['',3.7449989,-1.1638202,5.75];
ybs[5311]=['4 UMi',3.7035956,1.3517738,4.82];
ybs[5312]=['',3.7325275,-0.1054602,6.36];
ybs[5313]=['14 Boo',3.7309875,0.2245268,5.54];
ybs[5314]=['',3.7359445,-0.5127195,6.08];
ybs[5315]=['',3.7391653,-0.7870626,6.31];
ybs[5316]=['',3.7440188,-1.0473396,6.39];
ybs[5317]=['',3.7856953,-1.4475822,6.42];
ybs[5318]=['κ1 Boo',3.7271209,0.9022053,6.69];
ybs[5319]=['κ2 Boo',3.7272152,0.9022491,4.54];
ybs[5320]=['15 Boo',3.73436,0.1746335,5.29];
ybs[5321]=['',3.7346529,0.0565718,6.45];
ybs[5322]=['',3.7373455,-0.319316,5.43];
ybs[5323]=['',3.7334004,0.3801062,6.39];
ybs[5324]=['',3.7195376,1.2101564,5.24];
ybs[5325]=['',3.7315945,0.7229841,6.24];
ybs[5326]=['ε Aps',3.7741534,-1.3997771,5.06];
ybs[5327]=['',3.7416555,-0.5818186,6.55];
ybs[5328]=['ι Vir',3.7397729,-0.1063782,4.08];
ybs[5329]=['δ Oct',3.7981398,-1.4618654,4.32];
ybs[5330]=['α Boo',3.7377259,0.3331471,-0.04];
ybs[5331]=['',3.7412835,-0.1172217,6.44];
ybs[5332]=['',3.7418415,-0.0574338,6.15];
ybs[5333]=['',3.7395125,0.328427,5.98];
ybs[5334]=['',3.7446083,-0.3260177,6.22];
ybs[5335]=['',3.7350071,0.9152702,6.58];
ybs[5336]=['',3.7415646,0.3495381,6.25];
ybs[5337]=['',3.7404194,0.6920289,6.38];
ybs[5338]=['',3.7507914,-0.5814447,6.54];
ybs[5339]=['',3.7585209,-1.0710452,5.23];
ybs[5340]=['ι Boo',3.7389128,0.8948784,4.75];
ybs[5341]=['λ Boo',3.7401007,0.8027456,4.18];
ybs[5342]=['',3.7457108,0.2647536,5.8];
ybs[5343]=['',3.7485102,-0.1332804,6.47];
ybs[5344]=['ι Lup',3.7556268,-0.8054913,3.55];
ybs[5345]=['',3.7514871,-0.3282934,5.9];
ybs[5346]=['',3.7532896,-0.4522001,5.87];
ybs[5347]=['',3.7552797,-0.6474664,5.94];
ybs[5348]=['',3.7601935,-0.9857594,4.33];
ybs[5349]=['λ Vir',3.7534309,-0.2350033,4.52];
ybs[5350]=['',3.7440831,0.893837,6.2];
ybs[5351]=['',3.7474953,0.6181172,4.81];
ybs[5352]=['',3.7588225,-0.753147,5.56];
ybs[5353]=['',3.7462776,0.8361465,6.32];
ybs[5354]=['',3.7612947,-0.7902907,4.77];
ybs[5355]=['18 Boo',3.7536018,0.2253326,5.41];
ybs[5356]=['υ Vir',3.7550886,-0.0411727,5.14];
ybs[5357]=['ψ Cen',3.7603683,-0.6628486,4.05];
ybs[5358]=['',3.7556466,0.0050743,6.19];
ybs[5359]=['',3.751454,0.6749856,6.86];
ybs[5360]=['20 Boo',3.7556323,0.2829794,4.86];
ybs[5361]=['',3.7703634,-1.0219246,4.92];
ybs[5362]=['',3.7507602,0.955925,6.53];
ybs[5363]=['',3.7552297,0.6754505,6.33];
ybs[5364]=['',3.7570027,0.5294604,6.44];
ybs[5365]=['',3.7698951,-0.8449629,6.09];
ybs[5366]=['',3.7680112,-0.6087636,5.56];
ybs[5367]=['',3.7730453,-0.8877538,6.02];
ybs[5368]=['',3.7712575,-0.6912316,4.42];
ybs[5369]=['',3.7822851,-1.1918335,5.61];
ybs[5370]=['',3.7752232,-0.9297168,6];
ybs[5371]=['51 Hya',3.7711702,-0.4860098,4.77];
ybs[5372]=['',3.7844281,-1.1565412,6.36];
ybs[5373]=['2 Lib',3.7722377,-0.2060624,6.21];
ybs[5374]=['',3.7712205,0.0200586,6.27];
ybs[5375]=['',3.7716031,0.145781,6.86];
ybs[5376]=['',3.7716104,0.14581,5.12];
ybs[5377]=['',3.770075,0.4406188,6.22];
ybs[5378]=['',3.7743918,0.1422694,5.95];
ybs[5379]=['',3.8041695,-1.340753,6.07];
ybs[5380]=['',3.7785768,-0.4345575,5.32];
ybs[5381]=['',3.790777,-1.1503959,5.85];
ybs[5382]=['',3.7751985,0.0999703,5.1];
ybs[5383]=['',3.7777145,-0.2052803,6.49];
ybs[5384]=['',3.7756603,0.1395025,6.19];
ybs[5385]=['τ1 Lup',3.7850202,-0.7908592,4.56];
ybs[5386]=['τ2 Lup',3.7852162,-0.7936176,4.35];
ybs[5387]=['',3.7814578,-0.3501383,6.61];
ybs[5388]=['',3.7852829,-0.7402055,6.32];
ybs[5389]=['',3.7829302,-0.4702637,6.48];
ybs[5390]=['',3.7878511,-0.6975242,6.35];
ybs[5391]=['',3.7897211,-0.8067895,5.83];
ybs[5392]=['',3.7800456,0.6684835,6.27];
ybs[5393]=['',3.797144,-1.0347791,6.45];
ybs[5394]=['θ Boo',3.7782151,0.9033642,4.05];
ybs[5395]=['22 Boo',3.7848004,0.3339771,5.39];
ybs[5396]=['104 Vir',3.7894895,-0.1084101,6.17];
ybs[5397]=['52 Hya',3.7933876,-0.5163135,4.97];
ybs[5398]=['',3.8092407,-1.1834574,5.83];
ybs[5399]=['φ Vir',3.7928766,-0.0404739,4.81];
ybs[5400]=['106 Vir',3.7951285,-0.1220217,5.42];
ybs[5401]=['',3.7885305,0.7144297,6.63];
ybs[5402]=['',3.8025348,-0.7925883,5.5];
ybs[5403]=['',3.8036268,-0.8658473,5.37];
ybs[5404]=['',3.793592,0.4921534,7.62];
ybs[5405]=['',3.7937228,0.4921826,7.12];
ybs[5406]=['',3.7922622,0.6301687,6.1];
ybs[5407]=['',3.8058438,-0.7144516,6.39];
ybs[5408]=['',3.7999681,0.0128884,5.94];
ybs[5409]=['',3.8068133,-0.6799752,5.97];
ybs[5410]=['24 Boo',3.7932712,0.8683689,5.59];
ybs[5411]=['',3.8136919,-0.9944417,6.93];
ybs[5412]=['',3.7991751,0.5532807,6.06];
ybs[5413]=['',3.797906,0.7278947,6.35];
ybs[5414]=['',3.8038761,0.0817174,6.02];
ybs[5415]=['σ Lup',3.8135843,-0.8822024,4.42];
ybs[5416]=['',3.8178981,-0.9614594,5.87];
ybs[5417]=['',3.817571,-0.9209972,5.87];
ybs[5418]=['',3.8152007,-0.5376234,6.09];
ybs[5419]=['ρ Boo',3.807942,0.5285125,3.58];
ybs[5420]=['5 UMi',3.7852035,1.319553,4.25];
ybs[5421]=['',3.8198377,-0.7363334,6.6];
ybs[5422]=['',3.8259203,-1.0490216,6.4];
ybs[5423]=['',3.8102535,0.46404,6.01];
ybs[5424]=['26 Boo',3.811262,0.3869461,5.92];
ybs[5425]=['γ Boo',3.8087741,0.6670395,3.03];
ybs[5426]=['',3.8016198,1.1012211,6.09];
ybs[5427]=['',3.8060119,1.0495645,6.27];
ybs[5428]=['',3.8222827,-0.358282,6.5];
ybs[5429]=['',3.8258955,-0.7261592,5.87];
ybs[5430]=['η Cen',3.8258407,-0.7373391,2.31];
ybs[5431]=['',3.8143116,0.6435038,6.43];
ybs[5432]=['',3.8098665,0.9653084,5.76];
ybs[5433]=['',3.837729,-1.1871743,6.04];
ybs[5434]=['',3.8295555,-0.808675,5.55];
ybs[5435]=['',3.8181874,0.5662777,6.33];
ybs[5436]=['',3.8296598,-0.6926442,6.13];
ybs[5437]=['σ Boo',3.8203843,0.5175954,4.46];
ybs[5438]=['',3.820003,0.6376883,6.03];
ybs[5439]=['',3.8311392,-0.7033664,5.74];
ybs[5440]=['',3.8340008,-0.8067253,5.41];
ybs[5441]=['',3.8173682,0.9944212,6.48];
ybs[5442]=['',3.8195577,0.8600866,5.74];
ybs[5443]=['ρ Lup',3.8365741,-0.8641773,4.05];
ybs[5444]=['',3.8268181,0.404249,6.38];
ybs[5445]=['',3.831481,-0.2163071,6.2];
ybs[5446]=['',3.8380329,-0.6786176,6.02];
ybs[5447]=['',3.8420998,-0.8145738,6.07];
ybs[5448]=['',3.8432206,-0.8577063,6.39];
ybs[5449]=['α1 Cen',3.8448593,-1.0632996,-0.01];
ybs[5450]=['α2 Cen',3.8448739,-1.0633044,1.33];
ybs[5451]=['',3.8486403,-0.9865971,6.3];
ybs[5452]=['',3.8361808,0.3178331,5.91];
ybs[5453]=['α Cir',3.858042,-1.1355395,3.19];
ybs[5454]=['',3.8352809,0.7601619,5.7];
ybs[5455]=['',3.8548546,-1.0245548,6.22];
ybs[5456]=['',3.8497183,-0.6321913,5.67];
ybs[5457]=['',3.8349396,0.9413512,5.85];
ybs[5458]=['33 Boo',3.8379858,0.7734735,5.39];
ybs[5459]=['α Lup',3.8541429,-0.8285936,2.3];
ybs[5460]=['α Aps',3.885557,-1.3810641,3.83];
ybs[5461]=['',3.8538688,-0.6611343,4];
ybs[5462]=['',3.8453837,0.3820245,6.1];
ybs[5463]=['',3.8470832,0.2346966,5.91];
ybs[5464]=['',3.8531667,-0.5414005,6.37];
ybs[5465]=['π1 Boo',3.8471009,0.2850349,4.94];
ybs[5466]=['π2 Boo',3.8471228,0.2850252,5.88];
ybs[5467]=['ζ Boo',3.8490049,0.2380879,4.83];
ybs[5468]=['ζ Boo',3.8490049,0.2380879,4.43];
ybs[5469]=['',3.8098068,1.3887719,6.26];
ybs[5470]=['31 Boo',3.8513091,0.1409342,4.86];
ybs[5471]=['32 Boo',3.8515693,0.2020017,5.56];
ybs[5472]=['',3.8699663,-1.0988815,5.36];
ybs[5473]=['',3.852114,0.367164,6.38];
ybs[5474]=['4 Lib',3.8589977,-0.437793,5.73];
ybs[5475]=['',3.8611878,-0.615397,4.05];
ybs[5476]=['',3.8679505,-1.0221234,6.11];
ybs[5477]=['μ Vir',3.8577907,-0.1002621,3.88];
ybs[5478]=['',3.8688482,-0.9719293,6.1];
ybs[5479]=['',3.8669978,-0.6157094,4.92];
ybs[5480]=['34 Boo',3.8585952,0.4614931,4.81];
ybs[5481]=['',4.1052685,-1.538576,6.48];
ybs[5482]=['',3.8509095,1.0677097,6.25];
ybs[5483]=['',3.859514,0.7046432,5.73];
ybs[5484]=['',3.8740529,-0.8294887,5.74];
ybs[5485]=['',3.8766769,-0.9157483,5.21];
ybs[5486]=['',3.8670098,-0.0262384,6.07];
ybs[5487]=['54 Hya',3.871127,-0.4455535,4.94];
ybs[5488]=['',3.8774672,-0.9126395,6.07];
ybs[5489]=['',3.8715508,-0.4055849,5.81];
ybs[5490]=['',3.8855538,-1.1637534,5.91];
ybs[5491]=['108 Vir',3.8683061,0.0110262,5.69];
ybs[5492]=['ο Boo',3.8667788,0.2945919,4.6];
ybs[5493]=['5 Lib',3.8706958,-0.2713119,6.33];
ybs[5494]=['',3.8717958,-0.3710803,6.4];
ybs[5495]=['ε Boo',3.8653957,0.4710528,5.12];
ybs[5496]=['ε Boo',3.8653957,0.4710383,2.7];
ybs[5497]=['',3.8671807,0.3281077,6.13];
ybs[5498]=['',3.8762755,-0.669783,5.94];
ybs[5499]=['',3.878457,-0.7617009,6.3];
ybs[5500]=['',3.8662743,0.5707706,6.28];
ybs[5501]=['109 Vir',3.8715295,0.0315478,3.72];
ybs[5502]=['',3.8705681,0.2626138,5.63];
ybs[5503]=['',3.8763693,-0.373668,6.06];
ybs[5504]=['55 Hya',3.8771288,-0.4487113,5.63];
ybs[5505]=['',3.8861287,-0.9905086,6.23];
ybs[5506]=['56 Hya',3.878765,-0.4567911,5.24];
ybs[5507]=['57 Hya',3.8797054,-0.4665442,5.77];
ybs[5508]=['',3.8791501,-0.2255729,6.35];
ybs[5509]=['',3.8829848,-0.6408694,6.04];
ybs[5510]=['',3.9065845,-1.2788497,5.6];
ybs[5511]=['',3.8855458,-0.4247406,5.68];
ybs[5512]=['',3.8831705,-0.0162685,6.14];
ybs[5513]=['μ Lib',3.8853078,-0.248414,5.31];
ybs[5514]=['',3.8803049,0.4238032,6.14];
ybs[5515]=['π1 Oct',3.9517587,-1.453984,5.65];
ybs[5516]=['58 Hya',3.8899188,-0.4894622,4.41];
ybs[5517]=['',3.9019995,-1.1151429,5.87];
ybs[5518]=['ο Lup',3.8963942,-0.7619918,4.32];
ybs[5519]=['',3.8830256,0.658457,6.16];
ybs[5520]=['α1 Lib',3.8913273,-0.2806654,5.15];
ybs[5521]=['α2 Lib',3.8921651,-0.28144,2.75];
ybs[5522]=['',3.8870913,0.4979743,5.8];
ybs[5523]=['38 Boo',3.8835436,0.8034075,5.74];
ybs[5524]=['',3.8884935,0.4158778,5.85];
ybs[5525]=['11 Lib',3.8924335,-0.0415875,4.94];
ybs[5526]=['',3.8923197,-0.005949,6.18];
ybs[5527]=['',3.884247,0.8951887,6.51];
ybs[5528]=['39 Boo',3.8850512,0.8488657,5.69];
ybs[5529]=['ζ Cir',3.9115856,-1.1532024,6.09];
ybs[5530]=['',3.9282684,-1.339432,5.34];
ybs[5531]=['',3.8890677,0.6490548,5.48];
ybs[5532]=['',3.8998841,-0.5351183,6.29];
ybs[5533]=['',3.9014476,-0.6612405,5.03];
ybs[5534]=['ξ Boo',3.8935307,0.3319197,4.55];
ybs[5535]=['π2 Oct',3.9643585,-1.4506592,5.65];
ybs[5536]=['',3.9146907,-1.0506207,5.2];
ybs[5537]=['',3.9386719,-1.3481006,5.93];
ybs[5538]=['12 Lib',3.9074899,-0.4315272,5.3];
ybs[5539]=['',3.9090661,-0.5826416,5.82];
ybs[5540]=['',3.9023422,0.2726486,6.4];
ybs[5541]=['θ Cir',3.9200357,-1.0971556,5.11];
ybs[5542]=['',3.8919448,1.033415,5.46];
ybs[5543]=['',3.9022807,0.3328333,6.01];
ybs[5544]=['ξ1 Lib',3.9073518,-0.2091044,5.8];
ybs[5545]=['',3.9365923,-1.3109661,6.2];
ybs[5546]=['',3.9172148,-0.9231302,5.38];
ybs[5547]=['ω Oct',3.9966206,-1.4811439,5.91];
ybs[5548]=['',3.9139406,-0.5923264,5.32];
ybs[5549]=['',3.91799,-0.8370745,5.64];
ybs[5550]=['',3.9203396,-0.8993409,6.64];
ybs[5551]=['',3.9178708,-0.6893662,6.36];
ybs[5552]=['',3.9172612,-0.5710433,6.06];
ybs[5553]=['β UMi',3.8862742,1.292794,2.08];
ybs[5554]=['ξ2 Lib',3.9177564,-0.200562,5.46];
ybs[5555]=['',3.920254,-0.5103209,6.29];
ybs[5556]=['',3.9250953,-0.8542367,6.35];
ybs[5557]=['',3.9147254,0.2507085,5.77];
ybs[5558]=['',3.9210704,-0.3751924,5.74];
ybs[5559]=['',3.9131525,0.5623158,6.12];
ybs[5560]=['16 Lib',3.9193922,-0.0772812,4.49];
ybs[5561]=['β Lup',3.9264854,-0.7542418,2.68];
ybs[5562]=['',3.9267007,-0.6979157,6.15];
ybs[5563]=['',3.920912,-0.0043437,5.53];
ybs[5564]=['',3.9181978,0.3747869,6.49];
ybs[5565]=['',3.9189272,0.2846077,5.71];
ybs[5566]=['κ Cen',3.9291911,-0.7362658,3.13];
ybs[5567]=['59 Hya',3.9264364,-0.4841226,5.65];
ybs[5568]=['17 Lib',3.9241003,-0.1961075,6.6];
ybs[5569]=['',3.9293195,-0.6625642,6.47];
ybs[5570]=['',3.9305062,-0.7546918,6.1];
ybs[5571]=['',3.9141641,0.8647542,5.63];
ybs[5572]=['18 Lib',3.9270244,-0.1959143,5.87];
ybs[5573]=['',3.926812,-0.0884895,6.09];
ybs[5574]=['',3.9287804,0.0783137,5.93];
ybs[5575]=['',3.9379973,-0.6656402,5.89];
ybs[5576]=['δ Lib',3.9360356,-0.1500818,4.92];
ybs[5577]=['',3.941139,-0.6010681,6.22];
ybs[5578]=['40 Boo',3.9287455,0.6838999,5.64];
ybs[5579]=['',3.917914,1.1493169,4.6];
ybs[5580]=['',3.9374478,-0.0494808,5.52];
ybs[5581]=['60 Hya',3.9415283,-0.4911408,5.85];
ybs[5582]=['',3.9348134,0.3833672,6.38];
ybs[5583]=['η Cir',3.9555015,-1.1189366,5.17];
ybs[5584]=['',3.9394991,-0.0038472,5.71];
ybs[5585]=['',3.9455322,-0.5711197,5.44];
ybs[5586]=['',3.8792596,1.4386347,5.64];
ybs[5587]=['',3.93286,0.8237505,6.37];
ybs[5588]=['',3.9672233,-1.2563411,6.52];
ybs[5589]=['',3.9436438,-0.0542959,6.61];
ybs[5590]=['ω Boo',3.9401115,0.4350802,4.81];
ybs[5591]=['110 Vir',3.9441777,0.0351143,4.4];
ybs[5592]=['β Boo',3.9388605,0.7035541,3.5];
ybs[5593]=['σ Lib',3.9500138,-0.4426327,3.29];
ybs[5594]=['',3.9533968,-0.7145408,6.41];
ybs[5595]=['π Lup',3.9554687,-0.8225691,4.72];
ybs[5596]=['π Lup',3.9554687,-0.8225691,4.82];
ybs[5597]=['',3.9560415,-0.7181294,5.15];
ybs[5598]=['',3.9354354,1.0493677,5.93];
ybs[5599]=['',3.9441014,0.6130709,5.51];
ybs[5600]=['',3.9493568,0.0944823,6.5];
ybs[5601]=['',3.9694293,-1.1406266,6.17];
ybs[5602]=['',3.9437351,0.7778055,6.65];
ybs[5603]=['',3.9463337,0.6019088,6.59];
ybs[5604]=['',3.9575613,-0.4514842,6.67];
ybs[5605]=['',3.9598314,-0.6342947,6.27];
ybs[5606]=['ψ Boo',3.9502403,0.4689444,4.54];
ybs[5607]=['',3.9656865,-0.8581105,5.77];
ybs[5608]=['44 Boo',3.9465413,0.8303441,4.76];
ybs[5609]=['',3.9610478,-0.5410001,5.96];
ybs[5610]=['',3.9603143,-0.3858946,6.17];
ybs[5611]=['',3.9765059,-1.1721826,5.76];
ybs[5612]=['ν Lib',3.960914,-0.2851008,5.2];
ybs[5613]=['',3.975672,-1.1121199,6.28];
ybs[5614]=['',3.9686158,-0.7096755,5.79];
ybs[5615]=['',3.9706866,-0.7495339,5.85];
ybs[5616]=['λ Lup',3.9716438,-0.7916289,4.05];
ybs[5617]=['47 Boo',3.9536586,0.8390228,5.57];
ybs[5618]=['',3.9912083,-1.2714031,6.01];
ybs[5619]=['',3.9456666,1.1491332,6.13];
ybs[5620]=['',3.9592375,0.6349046,6.35];
ybs[5621]=['',3.9649099,0.0946019,6.16];
ybs[5622]=['',3.9811694,-1.0733605,6.3];
ybs[5623]=['',3.9631177,0.3205081,6.02];
ybs[5624]=['45 Boo',3.9627614,0.4326887,4.93];
ybs[5625]=['',3.9569379,0.950821,5.25];
ybs[5626]=['',3.9769397,-0.6783976,5.98];
ybs[5627]=['',3.9828802,-0.9673045,5.54];
ybs[5628]=['46 Boo',3.9674914,0.4576878,5.67];
ybs[5629]=['',3.9700217,0.2296446,6.1];
ybs[5630]=['',3.9683791,0.436876,5.81];
ybs[5631]=['',3.9772919,-0.4609335,5.76];
ybs[5632]=['',3.9836082,-0.7915726,6.44];
ybs[5633]=['',3.9833898,-0.7916068,7.39];
ybs[5634]=['',3.9981967,-1.2244283,5.81];
ybs[5635]=['',3.9911498,-1.078955,6.32];
ybs[5636]=['κ1 Lup',3.9853454,-0.8519633,3.87];
ybs[5637]=['κ2 Lup',3.985455,-0.852065,5.69];
ybs[5638]=['',3.9661582,0.8722702,6.39];
ybs[5639]=['ζ Lup',3.9870904,-0.9106282,3.41];
ybs[5640]=['',3.9878783,-0.8428983,6.33];
ybs[5641]=['',3.9890001,-0.7780042,4.82];
ybs[5642]=['ι1 Lib',3.9854378,-0.3467573,4.54];
ybs[5643]=['',3.9899247,-0.631235,6.1];
ybs[5644]=['',3.9837389,0.3298614,5.89];
ybs[5645]=['',3.9902365,-0.4203451,6.47];
ybs[5646]=['ι2 Lib',3.9902283,-0.3442341,6.08];
ybs[5647]=['23 Lib',3.9910925,-0.4430476,6.45];
ybs[5648]=['',3.9929111,-0.4584814,5.84];
ybs[5649]=['',3.9865792,0.3352761,6.68];
ybs[5650]=['1 Lup',3.996295,-0.5514249,4.91];
ybs[5651]=['',4.0067558,-1.0642707,5.73];
ybs[5652]=['26 Lib',3.9955938,-0.3114332,6.17];
ybs[5653]=['',4.0026176,-0.8403553,5.95];
ybs[5654]=['δ Cir',4.0082415,-1.0652042,5.09];
ybs[5655]=['',3.9899868,0.3998148,6.3];
ybs[5656]=['ε Cir',4.0116443,-1.1115038,4.86];
ybs[5657]=['',4.0030145,-0.7254584,5.16];
ybs[5658]=['',4.0035844,-0.7602527,6.04];
ybs[5659]=['',3.9964942,-0.0973524,6.28];
ybs[5660]=['β Cir',4.0104994,-1.0275645,4.07];
ybs[5661]=['γ TrA',4.0180011,-1.1999636,2.89];
ybs[5662]=['',3.9748149,1.1816681,6.17];
ybs[5663]=['',3.9897123,0.6665254,6.2];
ybs[5664]=['',3.9921803,0.5534897,5.99];
ybs[5665]=['3 Ser',3.9977254,0.0849008,5.33];
ybs[5666]=['χ Boo',3.993951,0.5076968,5.26];
ybs[5667]=['',3.9920584,0.7347132,6.13];
ybs[5668]=['',4.0036841,-0.3922445,5.5];
ybs[5669]=['4 Ser',4.0005869,0.0051919,5.63];
ybs[5670]=['',4.0163683,-1.0571439,5.46];
ybs[5671]=['δ Boo',3.998235,0.5801442,3.47];
ybs[5672]=['',4.012124,-0.7179348,6.28];
ybs[5673]=['μ Lup',4.014127,-0.8368618,4.27];
ybs[5674]=['',4.0255105,-1.1790419,6.28];
ybs[5675]=['β Lib',4.0060317,-0.1650618,2.61];
ybs[5676]=['2 Lup',4.0102631,-0.527488,4.34];
ybs[5677]=['',4.0155331,-0.7131736,5.59];
ybs[5678]=['',4.014042,-0.5459925,6.18];
ybs[5679]=['',4.0179327,-0.6487382,6.2];
ybs[5680]=['',4.0120337,-0.0093401,5.89];
ybs[5681]=['',3.9918687,1.1741061,5.13];
ybs[5682]=['',4.0113275,0.3577748,5.7];
ybs[5683]=['',3.9925541,1.2020084,6.51];
ybs[5684]=['5 Ser',4.0157995,0.0295284,5.06];
ybs[5685]=['δ Lup',4.0261483,-0.7106994,3.22];
ybs[5686]=['',4.0270995,-0.712482,6.2];
ybs[5687]=['',4.0266066,-0.6683162,6.48];
ybs[5688]=['ν1 Lup',4.0298821,-0.8377588,5];
ybs[5689]=['ν2 Lup',4.0284419,-0.8445678,5.65];
ybs[5690]=['',4.0354459,-1.0599169,5.67];
ybs[5691]=['28 Lib',4.0232512,-0.3181981,6.17];
ybs[5692]=['',4.0157036,0.5662178,6.32];
ybs[5693]=['ο Lib',4.0237349,-0.2726393,6.3];
ybs[5694]=['γ Cir',4.0361974,-1.0365961,4.51];
ybs[5695]=['φ1 Lup',4.0278557,-0.6341445,3.56];
ybs[5696]=['',4.0223325,-0.0433922,6.35];
ybs[5697]=['',4.0239153,-0.1029298,5.54];
ybs[5698]=['ε Lup',4.0320726,-0.7812356,3.37];
ybs[5699]=['ο CrB',4.0185919,0.5156221,5.51];
ybs[5700]=['6 Ser',4.0233323,0.011214,5.35];
ybs[5701]=['',4.0229938,0.4343254,6.39];
ybs[5702]=['φ2 Lup',4.0337765,-0.6445589,4.54];
ybs[5703]=['',4.0500335,-1.1934516,5.89];
ybs[5704]=['11 UMi',4.001593,1.2522637,5.02];
ybs[5705]=['',4.0172695,0.905571,5.66];
ybs[5706]=['',4.0203671,0.7742491,6.19];
ybs[5707]=['7 Ser',4.0289128,0.2180831,6.28];
ybs[5708]=['50 Boo',4.0257443,0.5735393,5.37];
ybs[5709]=['υ Lup',4.040854,-0.694319,5.37];
ybs[5710]=['',4.0360688,-0.2171381,5.72];
ybs[5711]=['8 Ser',4.0351383,-0.0190978,6.12];
ybs[5712]=['',4.0423643,-0.667424,7.03];
ybs[5713]=['ε Lib',4.0374441,-0.1814052,4.94];
ybs[5714]=['',4.0433794,-0.677269,4.6];
ybs[5715]=['',4.0550979,-1.1275084,5.71];
ybs[5716]=['',4.0290145,0.6895653,5.5];
ybs[5717]=['η CrB',4.0319245,0.5273653,5.58];
ybs[5718]=['η CrB',4.0319245,0.5273653,6.08];
ybs[5719]=['ρ Oct',4.1376497,-1.4752983,5.57];
ybs[5720]=['κ1 Aps',4.074434,-1.2820816,5.49];
ybs[5721]=['',4.0273804,1.0816667,5.98];
ybs[5722]=['',4.0351028,0.7888793,6.01];
ybs[5723]=['μ1 Boo',4.0372483,0.651108,4.31];
ybs[5724]=['μ2 Boo',4.0373586,0.6505943,6.5];
ybs[5725]=['γ UMi',4.0173671,1.2524623,3.05];
ybs[5726]=['',4.0518828,-0.6429454,5.45];
ybs[5727]=['',4.0272953,1.1042542,5.79];
ybs[5728]=['',4.0577172,-0.9017644,6.1];
ybs[5729]=['τ1 Ser',4.0436758,0.2680322,5.17];
ybs[5730]=['',4.0439834,0.3387671,6.27];
ybs[5731]=['',4.0452211,0.598038,5.46];
ybs[5732]=['',4.0615604,-0.8168526,5.24];
ybs[5733]=['ζ1 Lib',4.055342,-0.2929766,5.64];
ybs[5734]=['ι Dra',4.0377732,1.0279071,3.29];
ybs[5735]=['',4.0514794,0.436881,6.02];
ybs[5736]=['10 Ser',4.0564765,0.0309344,5.17];
ybs[5737]=['β CrB',4.0521259,0.506768,3.68];
ybs[5738]=['',4.0452346,0.9415925,6.45];
ybs[5739]=['',4.0657163,-0.362982,6.22];
ybs[5740]=['ζ3 Lib',4.0658873,-0.2910935,5.82];
ybs[5741]=['',4.0727972,-0.6752931,6.25];
ybs[5742]=['',4.0552852,0.8226006,6.15];
ybs[5743]=['',4.0715276,-0.5750791,6.46];
ybs[5744]=['',4.0493712,1.0856862,6.5];
ybs[5745]=['',4.0503326,1.0576702,5.9];
ybs[5746]=['',4.0705803,-0.3531374,6.22];
ybs[5747]=['',4.1107363,-1.3610641,6.18];
ybs[5748]=['',4.0662652,0.1485319,6.57];
ybs[5749]=['',4.0556145,0.9621165,6.43];
ybs[5750]=['',4.0631642,0.5448385,6.46];
ybs[5751]=['',4.0633197,0.6411472,6.37];
ybs[5752]=['',4.0744425,-0.3445064,5.52];
ybs[5753]=['ν1 Boo',4.0651629,0.7114677,5.02];
ybs[5754]=['ζ4 Lib',4.0757011,-0.2953248,5.5];
ybs[5755]=['',4.0769567,-0.4286079,7];
ybs[5756]=['',4.0934943,-1.1463309,6.51];
ybs[5757]=['',4.0813961,-0.7004655,5.82];
ybs[5758]=['',4.0566777,1.0826288,6.38];
ybs[5759]=['',4.0672875,0.6378761,6.38];
ybs[5760]=['τ2 Ser',4.0714493,0.2790376,6.22];
ybs[5761]=['ε TrA',4.0954957,-1.1586081,4.11];
ybs[5762]=['11 Ser',4.0754463,-0.0218948,5.51];
ybs[5763]=['',4.082767,-0.6879551,6.36];
ybs[5764]=['ν2 Boo',4.0688792,0.7126323,5.02];
ybs[5765]=['36 Lib',4.0834984,-0.4906876,5.15];
ybs[5766]=['γ Lup',4.0863198,-0.7196707,2.78];
ybs[5767]=['37 Lib',4.0809922,-0.1768374,4.62];
ybs[5768]=['θ CrB',4.074287,0.5461313,4.14];
ybs[5769]=['',4.0816019,-0.1005751,6.51];
ybs[5770]=['',4.0821233,-0.1614573,5.17];
ybs[5771]=['',4.0897744,-0.7858424,4.54];
ybs[5772]=['κ2 Aps',4.1131763,-1.2830179,5.65];
ybs[5773]=['',4.0789114,0.2979282,6.45];
ybs[5774]=['',4.0911192,-0.7760373,5.43];
ybs[5775]=['',4.0632485,1.1194465,5.79];
ybs[5776]=['',4.1211063,-1.329,5.95];
ybs[5777]=['γ Lib',4.0870161,-0.2592946,3.91];
ybs[5778]=['δ Ser',4.083128,0.1827383,3.8];
ybs[5779]=['δ Ser',4.083128,0.1827674,3.8];
ybs[5780]=['',4.0905383,-0.5787428,6.24];
ybs[5781]=['',4.0845801,0.0279539,6.56];
ybs[5782]=['',4.1115088,-1.2268396,6.44];
ybs[5783]=['α CrB',4.0821241,0.4650829,2.23];
ybs[5784]=['υ Lib',4.0939948,-0.4922076,3.58];
ybs[5785]=['τ3 Ser',4.086194,0.3069769,6.12];
ybs[5786]=['',4.0878545,0.1954528,6.07];
ybs[5787]=['ω Lup',4.0991089,-0.7440944,4.33];
ybs[5788]=['',4.1030978,-0.9152228,5.44];
ybs[5789]=['14 Ser',4.0911224,-0.0109661,6.51];
ybs[5790]=['μ CrB',4.0840801,0.6796794,5.11];
ybs[5791]=['',4.095893,-0.4598286,6.19];
ybs[5792]=['16 Ser',4.0905231,0.1735436,5.26];
ybs[5793]=['',4.1086799,-1.0467346,5.95];
ybs[5794]=['τ5 Ser',4.0903195,0.2801636,5.93];
ybs[5795]=['',4.1010869,-0.6846336,6.57];
ybs[5796]=['',4.0972112,-0.4050521,5.78];
ybs[5797]=['',4.1017841,-0.6840603,6.04];
ybs[5798]=['',4.0866014,0.6685813,6.42];
ybs[5799]=['',4.0993966,-0.4934496,6.32];
ybs[5800]=['',4.0991927,-0.3679509,5.84];
ybs[5801]=['',4.0832712,0.9399415,5.97];
ybs[5802]=['τ Lib',4.1011796,-0.5208679,3.66];
ybs[5803]=['',4.0916056,0.522282,6.52];
ybs[5804]=['41 Lib',4.1019165,-0.3380286,5.38];
ybs[5805]=['',4.1005432,-0.1546401,6.5];
ybs[5806]=['',4.1005504,-0.1545819,6.48];
ybs[5807]=['',4.0868854,0.9076198,6.74];
ybs[5808]=['',4.0861751,0.9523138,5.74];
ybs[5809]=['',4.1039925,-0.4051914,6.34];
ybs[5810]=['ψ1 Lup',4.1062106,-0.6017412,4.67];
ybs[5811]=['',4.1121546,-0.834273,6.23];
ybs[5812]=['',4.1082208,-0.5459164,6.34];
ybs[5813]=['φ Boo',4.0952482,0.7031432,5.24];
ybs[5814]=['42 Lib',4.1080605,-0.4168397,4.96];
ybs[5815]=['',4.1129169,-0.7806124,4.64];
ybs[5816]=['θ UMi',4.0616016,1.3487977,4.96];
ybs[5817]=['',4.099813,0.6040448,6.11];
ybs[5818]=['',4.0930753,0.9502015,5.87];
ybs[5819]=['',4.0495974,1.4028711,6.58];
ybs[5820]=['',4.0968317,0.8156228,5.75];
ybs[5821]=['',4.1065164,0.2092277,6.25];
ybs[5822]=['',4.1194714,-0.8648722,6.04];
ybs[5823]=['ζ1 CrB',4.1021482,0.6382862,6];
ybs[5824]=['ζ2 CrB',4.1021846,0.6382717,5.07];
ybs[5825]=['',4.0979111,0.8789027,5.84];
ybs[5826]=['',4.1261092,-1.0533188,6.48];
ybs[5827]=['',4.1188826,-0.6543084,5.24];
ybs[5828]=['κ Lib',4.1151912,-0.3445858,4.74];
ybs[5829]=['ψ2 Lup',4.1189595,-0.6069321,4.75];
ybs[5830]=['τ6 Ser',4.1099436,0.2785521,6.01];
ybs[5831]=['',4.0998403,1.0098254,6.45];
ybs[5832]=['ι Ser',4.1122928,0.3421829,4.52];
ybs[5833]=['χ Ser',4.1135448,0.2231049,5.33];
ybs[5834]=['',4.0915855,1.2080634,5.62];
ybs[5835]=['τ7 Ser',4.1139017,0.3211301,5.81];
ybs[5836]=['',4.1266889,-0.7309931,5.94];
ybs[5837]=['',4.1214537,-0.2636696,6.31];
ybs[5838]=['η Lib',4.1243467,-0.2746508,5.41];
ybs[5839]=['γ CrB',4.1172654,0.4578242,3.84];
ybs[5840]=['',4.1195688,0.2374313,6.48];
ybs[5841]=['',4.1442193,-1.1432659,6.18];
ybs[5842]=['',4.1442266,-1.1432658,6.39];
ybs[5843]=['ψ Ser',4.1236206,0.0427851,5.88];
ybs[5844]=['α Ser',4.1245469,0.1110388,2.65];
ybs[5845]=['π CrB',4.1224579,0.5663974,5.56];
ybs[5846]=['',4.1341077,-0.4908618,6.51];
ybs[5847]=['',4.1163885,0.9127487,5.51];
ybs[5848]=['τ8 Ser',4.1261091,0.3002112,6.14];
ybs[5849]=['',4.1294763,0.0939669,5.58];
ybs[5850]=['',4.1366578,-0.6064131,5.61];
ybs[5851]=['',4.1307869,0.0144596,6.33];
ybs[5852]=['',4.139912,-0.7026046,6.42];
ybs[5853]=['25 Ser',4.1327487,-0.0325884,5.4];
ybs[5854]=['',4.1400637,-0.6628495,6.01];
ybs[5855]=['',4.1468572,-0.9162898,6.07];
ybs[5856]=['',4.1357694,-0.107909,6.24];
ybs[5857]=['β Ser',4.1326517,0.2680691,3.67];
ybs[5858]=['λ Ser',4.1340085,0.1272424,4.43];
ybs[5859]=['',4.1524978,-0.9297437,5.77];
ybs[5860]=['υ Ser',4.1374907,0.2452715,5.71];
ybs[5861]=['',4.1514921,-0.8547396,5.84];
ybs[5862]=['',4.1526369,-0.7934716,6.12];
ybs[5863]=['',4.1570189,-0.9619618,5.73];
ybs[5864]=['',4.1415721,0.2395769,6];
ybs[5865]=['',4.1452582,-0.0677213,5.53];
ybs[5866]=['',4.1662631,-1.1381672,6.54];
ybs[5867]=['',4.140127,0.5528131,6.44];
ybs[5868]=['',4.1324014,0.9671232,5.92];
ybs[5869]=['κ Ser',4.1436954,0.3155557,4.09];
ybs[5870]=['',4.142621,0.4903489,5.85];
ybs[5871]=['μ Ser',4.148185,-0.0609387,3.53];
ybs[5872]=['',4.158214,-0.8224153,6.01];
ybs[5873]=['χ Lup',4.1550465,-0.5879641,3.95];
ybs[5874]=['',4.1677736,-1.093736,6.19];
ybs[5875]=['1 Sco',4.1548253,-0.450505,4.64];
ybs[5876]=['',4.1319893,1.091473,5.19];
ybs[5877]=['',4.1369989,0.9654195,5.86];
ybs[5878]=['ω Ser',4.1509504,0.0372699,5.23];
ybs[5879]=['δ CrB',4.1471534,0.4539081,4.63];
ybs[5880]=['',4.1642821,-0.8844468,6.6];
ybs[5881]=['κ TrA',4.1780806,-1.1983713,5.09];
ybs[5882]=['ε Ser',4.1531737,0.0770915,3.71];
ybs[5883]=['',4.1603732,-0.52267,6.4];
ybs[5884]=['',4.1523191,0.2630696,5.2];
ybs[5885]=['36 Ser',4.1553313,-0.0549975,5.11];
ybs[5886]=['',4.1573238,-0.247732,6.19];
ybs[5887]=['β TrA',4.1755952,-1.1080975,2.85];
ybs[5888]=['',4.1740695,-1.061199,6.15];
ybs[5889]=['ρ Ser',4.1546149,0.3650736,4.76];
ybs[5890]=['',4.1768923,-1.0513233,5.77];
ybs[5891]=['κ CrB',4.1539062,0.6212824,4.82];
ybs[5892]=['λ Lib',4.1649178,-0.3530257,5.03];
ybs[5893]=['ζ UMi',4.1162087,1.3566532,4.32];
ybs[5894]=['2 Sco',4.1663056,-0.4430825,4.59];
ybs[5895]=['',4.1793846,-1.0566423,5.76];
ybs[5896]=['',4.1675212,-0.4292196,5.39];
ybs[5897]=['',4.1676468,-0.4195328,5.42];
ybs[5898]=['θ Lib',4.1669432,-0.2930216,4.15];
ybs[5899]=['',4.1620245,0.3027051,6.36];
ybs[5900]=['',4.1702581,-0.4781812,6.14];
ybs[5901]=['39 Ser',4.1633144,0.2292822,6.1];
ybs[5902]=['3 Sco',4.1708718,-0.4416154,5.87];
ybs[5903]=['',4.1648808,0.2795213,6.09];
ybs[5904]=['χ Her',4.1598654,0.7398734,4.62];
ybs[5905]=['47 Lib',4.1721809,-0.339327,5.94];
ybs[5906]=['',4.1747997,-0.5435362,6.21];
ybs[5907]=['4 Sco',4.1745909,-0.4594503,5.62];
ybs[5908]=['',4.1778243,-0.696781,6.03];
ybs[5909]=['40 Ser',4.1698693,0.148722,6.29];
ybs[5910]=['',4.192583,-1.13612,5.75];
ybs[5911]=['',4.1824791,-0.8416018,6.31];
ybs[5912]=['',4.1571864,0.9733073,5.81];
ybs[5913]=['',4.1779938,-0.5557868,6.29];
ybs[5914]=['',4.1690815,0.3534578,5.44];
ybs[5915]=['ξ1 Lup',4.1809681,-0.5938397,5.12];
ybs[5916]=['ξ2 Lup',4.1810189,-0.5938008,5.62];
ybs[5917]=['',4.1774088,-0.2523375,6.37];
ybs[5918]=['ρ Sco',4.1807427,-0.510898,3.88];
ybs[5919]=['',4.1830905,-0.6325582,5.8];
ybs[5920]=['',4.1788046,-0.2598353,6.13];
ybs[5921]=['',4.1738778,0.323965,6.26];
ybs[5922]=['2 Her',4.1683569,0.7518772,5.37];
ybs[5923]=['γ Ser',4.1774205,0.2723286,3.85];
ybs[5924]=['',4.1838823,-0.3672323,5.85];
ybs[5925]=['',4.188203,-0.6555539,6.31];
ybs[5926]=['λ CrB',4.173689,0.6612744,5.45];
ybs[5927]=['',4.1953212,-0.943837,6.1];
ybs[5928]=['4 Her',4.1722174,0.7418918,5.75];
ybs[5929]=['',4.2020411,-1.1140929,6.41];
ybs[5930]=['φ Ser',4.1809066,0.2505664,5.54];
ybs[5931]=['48 Lib',4.1859156,-0.2502285,4.88];
ybs[5932]=['',4.1879764,-0.4343915,5.43];
ybs[5933]=['',4.1927566,-0.7295724,4.99];
ybs[5934]=['π Sco',4.1892097,-0.4567781,2.89];
ybs[5935]=['',4.1947145,-0.7105158,6.49];
ybs[5936]=['',4.2006429,-0.9535434,6.13];
ybs[5937]=['ε CrB',4.1819829,0.4680948,4.15];
ybs[5938]=['η Lup',4.1952768,-0.671143,3.41];
ybs[5939]=['',4.1723437,1.0271769,6.31];
ybs[5940]=['',4.1810363,0.6918013,6.31];
ybs[5941]=['',4.2092668,-1.0925254,6.25];
ybs[5942]=['',4.1987571,-0.7067127,6.21];
ybs[5943]=['δ Sco',4.1955506,-0.3958115,2.32];
ybs[5944]=['49 Lib',4.1953142,-0.2895502,5.47];
ybs[5945]=['',4.2247911,-1.264575,5.7];
ybs[5946]=['',4.2002358,-0.5575569,6.33];
ybs[5947]=['',4.1875643,0.6385556,5.62];
ybs[5948]=['',4.1903703,0.4513977,6.33];
ybs[5949]=['50 Lib',4.1970931,-0.1477922,5.55];
ybs[5950]=['',4.1813056,0.9545522,4.95];
ybs[5951]=['ι1 Nor',4.2115538,-1.0093317,4.63];
ybs[5952]=['η Nor',4.209425,-0.8601866,4.65];
ybs[5953]=['',4.1969576,0.0762888,5.83];
ybs[5954]=['',4.1872895,0.8695889,6.05];
ybs[5955]=['',4.2059371,-0.5094873,6.03];
ybs[5956]=['5 Her',4.1982221,0.3100055,5.12];
ybs[5957]=['',4.209617,-0.6747057,4.89];
ybs[5958]=['ρ CrB',4.1967982,0.5802726,5.41];
ybs[5959]=['',4.2088083,-0.4524003,5];
ybs[5960]=['',4.2100494,-0.5594791,6.01];
ybs[5961]=['ι CrB',4.1986789,0.5200182,4.99];
ybs[5962]=['π Ser',4.2026592,0.3970373,4.83];
ybs[5963]=['',4.2112469,-0.4325186,6.21];
ybs[5964]=['',4.2132668,-0.5806599,6.1];
ybs[5965]=['',4.2148568,-0.6617909,5.9];
ybs[5966]=['43 Ser',4.2096368,0.0860701,6.08];
ybs[5967]=['ξ Sco',4.2127872,-0.1994559,5.07];
ybs[5968]=['ξ Sco',4.2127872,-0.1994559,4.77];
ybs[5969]=['',4.2282936,-0.9816581,6.16];
ybs[5970]=['δ Nor',4.2234605,-0.7893646,4.72];
ybs[5971]=['',4.2002079,0.9225773,5.93];
ybs[5972]=['υ Her',4.2037839,0.802519,4.76];
ybs[5973]=['',4.2065823,0.6383752,5.83];
ybs[5974]=['β1 Sco',4.217728,-0.3466224,2.62];
ybs[5975]=['β2 Sco',4.2177497,-0.3465593,4.92];
ybs[5976]=['θ Dra',4.1987802,1.0211767,4.01];
ybs[5977]=['θ Lup',4.2234511,-0.6432607,4.23];
ybs[5978]=['',4.2207801,-0.4129542,5.92];
ybs[5979]=['',4.2186214,-0.1107586,6.53];
ybs[5980]=['',4.2197294,-0.1080998,6.41];
ybs[5981]=['',4.2264253,-0.642441,5.73];
ybs[5982]=['',4.2176919,0.1403543,6.29];
ybs[5983]=['ω1 Sco',4.2237379,-0.3616847,3.96];
ybs[5984]=['ι2 Nor',4.2368214,-1.0120653,5.57];
ybs[5985]=['',4.2042017,1.035944,6.19];
ybs[5986]=['',4.2246065,-0.2465204,6.32];
ybs[5987]=['ω2 Sco',4.2263565,-0.3651611,4.32];
ybs[5988]=['',4.2284941,-0.427873,6.33];
ybs[5989]=['',4.2322034,-0.6834414,7.05];
ybs[5990]=['',4.2322173,-0.6832279,6.65];
ybs[5991]=['',4.2297067,-0.4604164,5.38];
ybs[5992]=['11 Sco',4.2269638,-0.2233857,5.78];
ybs[5993]=['',4.232234,-0.4143157,5.88];
ybs[5994]=['45 Ser',4.226339,0.1717079,5.63];
ybs[5995]=['',4.2248293,0.3799379,6.14];
ybs[5996]=['',4.2360836,-0.5707636,6.19];
ybs[5997]=['',4.2376438,-0.5864008,5.54];
ybs[5998]=['κ Her',4.2280694,0.2965943,5];
ybs[5999]=['κ Her',4.2280982,0.2967252,6.25];
ybs[6000]=['47 Ser',4.2300538,0.1480215,5.73];
ybs[6001]=['',4.2324588,0.0593677,5.91];
ybs[6002]=['',4.2372462,-0.3210237,6.47];
ybs[6003]=['8 Her',4.2311247,0.2993728,6.14];
ybs[6004]=['',4.2332625,0.1104104,5.97];
ybs[6005]=['',4.2442091,-0.7185788,5.86];
ybs[6006]=['',4.236425,-0.0614264,5.37];
ybs[6007]=['',4.2425177,-0.5143196,5.13];
ybs[6008]=['τ CrB',4.2312172,0.6359604,4.76];
ybs[6009]=['ζ Nor',4.2543427,-0.9702576,5.81];
ybs[6010]=['δ1 Aps',4.2912394,-1.3743281,4.68];
ybs[6011]=['δ2 Aps',4.2916538,-1.373828,5.27];
ybs[6012]=['',4.2537465,-0.9376354,5.83];
ybs[6013]=['φ Her',4.2298833,0.7833369,4.26];
ybs[6014]=['κ Nor',4.2547002,-0.9543696,4.94];
ybs[6015]=['',4.2166976,1.1825649,5.44];
ybs[6016]=['ν Sco',4.2462612,-0.3403616,6.3];
ybs[6017]=['ν Sco',4.2463416,-0.3405505,4.01];
ybs[6018]=['13 Sco',4.2480018,-0.4883045,4.59];
ybs[6019]=['12 Sco',4.2478541,-0.4968762,5.67];
ybs[6020]=['δ TrA',4.2643641,-1.1123928,3.85];
ybs[6021]=['ψ Sco',4.2460508,-0.1765528,4.94];
ybs[6022]=['',4.2432281,0.1686107,6.53];
ybs[6023]=['16 Sco',4.2465334,-0.1500811,5.43];
ybs[6024]=['',4.2013473,1.3393296,5.56];
ybs[6025]=['',4.242928,0.289964,6.08];
ybs[6026]=['',4.2300424,1.0102795,6.33];
ybs[6027]=['',4.2723946,-1.1866575,5.75];
ybs[6028]=['',4.2319372,0.9734755,6.49];
ybs[6029]=['10 Her',4.243364,0.4091565,5.7];
ybs[6030]=['',4.265301,-1.0116269,5.63];
ybs[6031]=['',4.2499756,-0.0745602,6.25];
ybs[6032]=['',4.2542323,-0.4271291,6.41];
ybs[6033]=['',4.2430915,0.5810325,6.29];
ybs[6034]=['',4.2572447,-0.5770382,5.92];
ybs[6035]=['θ Nor',4.2618839,-0.8276745,5.14];
ybs[6036]=['',4.2435529,0.6348332,5.63];
ybs[6037]=['9 Her',4.251062,0.0867443,5.48];
ybs[6038]=['χ Sco',4.2541742,-0.207489,5.22];
ybs[6039]=['',4.2622354,-0.7496137,6.14];
ybs[6040]=['',4.2432145,0.7386703,5.87];
ybs[6041]=['',4.2572462,-0.3692758,6.41];
ybs[6042]=['',4.248136,0.4645987,6.5];
ybs[6043]=['',4.2579048,-0.3243807,6.32];
ybs[6044]=['',4.2592122,-0.4455336,6.05];
ybs[6045]=['',4.2687794,-0.9400424,5.44];
ybs[6046]=['δ Oph',4.2560827,-0.0653621,2.74];
ybs[6047]=['',4.2552592,0.1021253,6.31];
ybs[6048]=['γ1 Nor',4.2697653,-0.8747165,4.99];
ybs[6049]=['',4.2714673,-0.9273935,6.33];
ybs[6050]=['18 Sco',4.2618013,-0.1469462,5.5];
ybs[6051]=['',4.2630456,-0.2600366,6.09];
ybs[6052]=['',4.2785677,-1.0113848,6.49];
ybs[6053]=['σ CrB',4.2562306,0.5900637,5.64];
ybs[6054]=['σ CrB',4.2562306,0.5900589,6.66];
ybs[6055]=['16 Her',4.2603028,0.3273936,5.69];
ybs[6056]=['',4.2681652,-0.3726838,6.61];
ybs[6057]=['',4.2673342,-0.0698605,6.18];
ybs[6058]=['',4.2613432,0.4777365,6.14];
ybs[6059]=['',4.2433502,1.1709861,6.21];
ybs[6060]=['',4.2741998,-0.5002568,4.78];
ybs[6061]=['λ Nor',4.2792307,-0.7456416,5.45];
ybs[6062]=['γ2 Nor',4.2821109,-0.8762166,4.02];
ybs[6063]=['',4.2850743,-0.9632065,5.77];
ybs[6064]=['υ CrB',4.2654458,0.5079041,5.78];
ybs[6065]=['ε Oph',4.2734667,-0.0827504,3.24];
ybs[6066]=['',4.2775031,-0.3537108,6.29];
ybs[6067]=['',4.279738,-0.5402634,5.49];
ybs[6068]=['',4.2767867,-0.2604192,5.94];
ybs[6069]=['19 UMi',4.23359,1.3233966,5.48];
ybs[6070]=['',4.2845013,-0.68903,6.12];
ybs[6071]=['ο Sco',4.2842322,-0.4226683,4.55];
ybs[6072]=['20 UMi',4.2413616,1.3117697,6.39];
ybs[6073]=['',4.2935386,-0.8660147,5.33];
ybs[6074]=['σ Sco',4.2866949,-0.4475057,2.89];
ybs[6075]=['',4.2932201,-0.7672294,5.88];
ybs[6076]=['',4.265587,1.0420595,5.4];
ybs[6077]=['',4.2802562,0.3679939,6.05];
ybs[6078]=['',4.2509239,1.2800987,5.98];
ybs[6079]=['',4.3076678,-1.102526,6.15];
ybs[6080]=['',4.2750167,0.8550295,5.91];
ybs[6081]=['',4.2787859,0.6922063,5.46];
ybs[6082]=['τ Her',4.277618,0.8074788,3.89];
ybs[6083]=['σ Ser',4.2896559,0.017141,4.82];
ybs[6084]=['',4.299676,-0.6848524,5.4];
ybs[6085]=['γ Her',4.2883895,0.3334608,3.75];
ybs[6086]=['',4.2922634,-0.0371147,6.23];
ybs[6087]=['',4.2990477,-0.580245,6.47];
ybs[6088]=['ζ TrA',4.3226728,-1.22397,4.91];
ybs[6089]=['',4.3038853,-0.7922894,6.33];
ybs[6090]=['',4.3018037,-0.656448,5.42];
ybs[6091]=['',4.2680253,1.1956445,6.41];
ybs[6092]=['γ Aps',4.3488072,-1.3777378,3.89];
ybs[6093]=['ξ CrB',4.288705,0.538344,4.85];
ybs[6094]=['ψ Oph',4.2992106,-0.3505249,4.5];
ybs[6095]=['',4.3020159,-0.5192208,6.63];
ybs[6096]=['',4.3020232,-0.519245,5.84];
ybs[6097]=['ν1 CrB',4.28971,0.5890865,5.2];
ybs[6098]=['ν2 CrB',4.2902812,0.5874199,5.39];
ybs[6099]=['ι TrA',4.3191305,-1.118795,5.27];
ybs[6100]=['',4.2923272,0.5635028,6.4];
ybs[6101]=['21 Her',4.2986585,0.1204616,5.85];
ybs[6102]=['ρ Oph',4.305811,-0.4100239,5.02];
ybs[6103]=['ρ Oph',4.3058036,-0.4100045,5.92];
ybs[6104]=['',4.3196657,-1.0235274,5.69];
ybs[6105]=['ε Nor',4.3139978,-0.8307701,4.47];
ybs[6106]=['η UMi',4.262697,1.3213153,4.95];
ybs[6107]=['ω Her',4.30382,0.2441327,4.57];
ybs[6108]=['χ Oph',4.3118985,-0.3229062,4.42];
ybs[6109]=['',4.3053019,0.328944,6.7];
ybs[6110]=['',4.3261066,-1.0087869,6.06];
ybs[6111]=['',4.307289,0.1983095,6.11];
ybs[6112]=['',4.318011,-0.6496746,5.79];
ybs[6113]=['25 Her',4.3028231,0.6518504,5.54];
ybs[6114]=['',4.3103901,0.0401881,6.07];
ybs[6115]=['',4.3312889,-1.0764577,5.2];
ybs[6116]=['',4.2838048,1.2053596,5.25];
ybs[6117]=['',4.2973145,0.9627039,5.74];
ybs[6118]=['',4.3145932,-0.1333872,5.23];
ybs[6119]=['υ Oph',4.3149531,-0.1468885,4.63];
ybs[6120]=['',4.2937737,1.0759992,5.67];
ybs[6121]=['',4.3249121,-0.8078571,5.35];
ybs[6122]=['η Dra',4.2947069,1.0728156,2.74];
ybs[6123]=['',4.5713416,-1.5272734,6.57];
ybs[6124]=['α Sco',4.3226012,-0.4620865,0.96];
ybs[6125]=['',4.3485297,-1.2396935,5.5];
ybs[6126]=['',4.3179996,0.0108368,5.39];
ybs[6127]=['',4.3193671,-0.1426432,6.48];
ybs[6128]=['',4.4096495,-1.4534062,6.57];
ybs[6129]=['',4.4899958,-1.504783,6.04];
ybs[6130]=['',4.3238004,-0.2547194,5.68];
ybs[6131]=['22 Sco',4.3260409,-0.439095,4.79];
ybs[6132]=['',4.3333249,-0.7305863,5.33];
ybs[6133]=['',4.3315832,-0.6064526,4.23];
ybs[6134]=['',4.3266927,-0.1319153,6.5];
ybs[6135]=['',4.331219,-0.4639177,6.1];
ybs[6136]=['30 Her',4.3167017,0.7302025,5.04];
ybs[6137]=['φ Oph',4.3297928,-0.2906961,4.28];
ybs[6138]=['β Her',4.3245096,0.3743094,2.77];
ybs[6139]=['λ Oph',4.3281963,0.0338746,3.82];
ybs[6140]=['',4.316404,0.8964646,6.29];
ybs[6141]=['θ TrA',4.3534375,-1.143816,5.52];
ybs[6142]=['',4.3260301,0.3566747,5.25];
ybs[6143]=['ω Oph',4.3343253,-0.3753994,4.45];
ybs[6144]=['',4.3288599,0.3866318,5.76];
ybs[6145]=['μ Nor',4.3438876,-0.7694586,4.94];
ybs[6146]=['34 Her',4.3225967,0.8537685,6.45];
ybs[6147]=['',4.3275396,0.6140414,6.25];
ybs[6148]=['28 Her',4.335419,0.0956243,5.63];
ybs[6149]=['29 Her',4.3352655,0.199767,4.84];
ybs[6150]=['',4.348542,-0.7903843,6.46];
ybs[6151]=['15 Dra',4.3107612,1.1994504,5];
ybs[6152]=['',4.3301926,0.7950957,5.65];
ybs[6153]=['β Aps',4.3898603,-1.3535794,4.24];
ybs[6154]=['',4.3538205,-0.7487336,5.47];
ybs[6155]=['τ Sco',4.3509448,-0.4931736,2.82];
ybs[6156]=['',4.3534084,-0.6160308,4.16];
ybs[6157]=['',4.3663904,-1.0651642,6.18];
ybs[6158]=['σ Her',4.340484,0.7399379,4.2];
ybs[6159]=['',4.3474338,0.29699,6.41];
ybs[6160]=['',4.3315617,1.0608258,5.94];
ybs[6161]=['12 Oph',4.3520941,-0.0412807,5.75];
ybs[6162]=['η1 TrA',4.3786907,-1.1926532,5.91];
ybs[6163]=['',4.2962605,1.377379,5.56];
ybs[6164]=['',4.3628551,-0.7581368,5.83];
ybs[6165]=['ζ Oph',4.355857,-0.1851326,2.56];
ybs[6166]=['',4.3530654,0.2697879,6.3];
ybs[6167]=['',4.3748928,-1.0556507,6.18];
ybs[6168]=['',4.3653446,-0.650251,5.91];
ybs[6169]=['',4.3595119,-0.1148035,6.09];
ybs[6170]=['',4.3247983,1.2665658,6.3];
ybs[6171]=['',4.3578543,0.2381868,6.31];
ybs[6172]=['',4.3871435,-1.1775642,6.03];
ybs[6173]=['',4.3493011,0.8128462,5.79];
ybs[6174]=['16 Dra',4.3488317,0.9225737,5.53];
ybs[6175]=['17 Dra',4.3489894,0.9229958,5.08];
ybs[6176]=['17 Dra',4.3490185,0.922991,6.53];
ybs[6177]=['',4.375891,-0.8517399,5.65];
ybs[6178]=['',4.3774077,-0.8672463,5.65];
ybs[6179]=['',4.3667035,-0.167436,6.35];
ybs[6180]=['',4.3711149,-0.3568691,6.26];
ybs[6181]=['',4.3188427,1.3509388,6.34];
ybs[6182]=['',4.3768076,-0.579175,5.87];
ybs[6183]=['',4.3757589,-0.4277112,6.09];
ybs[6184]=['36 Her',4.3703228,0.0727576,6.93];
ybs[6185]=['37 Her',4.3705842,0.0729762,5.77];
ybs[6186]=['',4.3753741,-0.3103237,4.96];
ybs[6187]=['',4.3832007,-0.804733,6.23];
ybs[6188]=['',4.3507898,1.1001221,6.16];
ybs[6189]=['',4.3564434,0.9769601,5.29];
ybs[6190]=['42 Her',4.3603031,0.8532713,4.9];
ybs[6191]=['',4.3731406,-0.0181302,6.24];
ybs[6192]=['',4.3768579,-0.348408,5.57];
ybs[6193]=['',4.3712339,0.2156632,6.08];
ybs[6194]=['',4.4014853,-1.1719036,5.13];
ybs[6195]=['14 Oph',4.3753234,0.0199513,5.74];
ybs[6196]=['',4.3859528,-0.7183046,6.2];
ybs[6197]=['',4.3907757,-0.9283225,5.96];
ybs[6198]=['',4.3714461,0.4331951,6.06];
ybs[6199]=['',4.3865636,-0.7182017,6.12];
ybs[6200]=['',4.385941,-0.666604,6.05];
ybs[6201]=['',4.384995,-0.5610033,6.46];
ybs[6202]=['ζ Her',4.3723679,0.5509099,2.81];
ybs[6203]=['39 Her',4.3739873,0.4691246,5.92];
ybs[6204]=['',4.3900913,-0.7134244,5.71];
ybs[6205]=['',4.3987199,-1.0217025,5.74];
ybs[6206]=['',4.3875997,-0.4798404,6.58];
ybs[6207]=['α TrA',4.4106747,-1.205363,1.92];
ybs[6208]=['',4.390765,-0.4982235,6.02];
ybs[6209]=['',4.4029443,-1.0188632,5.58];
ybs[6210]=['η Her',4.3790285,0.6786662,3.53];
ybs[6211]=['',4.3991184,-0.6878819,5.48];
ybs[6212]=['',4.3834854,0.5934442,5.99];
ybs[6213]=['18 Dra',4.3679717,1.1266204,4.83];
ybs[6214]=['16 Oph',4.391851,0.0171753,6.03];
ybs[6215]=['25 Sco',4.398723,-0.4461781,6.71];
ybs[6216]=['',4.3781265,0.9713235,6.16];
ybs[6217]=['',4.3908286,0.2741735,5.56];
ybs[6218]=['43 Her',4.3930684,0.1491635,5.15];
ybs[6219]=['η Ara',4.4137746,-1.03106,3.76];
ybs[6220]=['',4.3888025,0.7536468,6.05];
ybs[6221]=['',4.4261853,-1.181844,6.32];
ybs[6222]=['19 Oph',4.3990959,0.0354132,6.1];
ybs[6223]=['',4.4239879,-1.1415936,6.13];
ybs[6224]=['45 Her',4.4016517,0.0909584,5.24];
ybs[6225]=['',4.4052834,-0.2608259,6.03];
ybs[6226]=['',4.4164249,-0.8740472,6.47];
ybs[6227]=['',4.3881381,0.9903957,4.85];
ybs[6228]=['',4.3491668,1.3766817,6.32];
ybs[6229]=['',4.4030036,0.2365895,6.35];
ybs[6230]=['',4.4097326,-0.2740481,6.1];
ybs[6231]=['ε Sco',4.4135513,-0.5991236,2.29];
ybs[6232]=['',4.3981725,0.7365893,5.87];
ybs[6233]=['20 Oph',4.4111761,-0.1887955,4.65];
ybs[6234]=['',4.4173523,-0.6553356,6.11];
ybs[6235]=['',4.4200258,-0.720189,5.22];
ybs[6236]=['',4.4092465,0.2308563,5.91];
ybs[6237]=['μ1 Sco',4.4211931,-0.6646318,3.08];
ybs[6238]=['',4.4132513,-0.0469106,6.32];
ybs[6239]=['',4.4233662,-0.7310716,6.49];
ybs[6240]=['47 Her',4.4127103,0.1259054,5.49];
ybs[6241]=['',4.4320816,-1.0112689,5.94];
ybs[6242]=['μ2 Sco',4.4232213,-0.6641043,3.57];
ybs[6243]=['',4.4389601,-1.1048109,6.02];
ybs[6244]=['52 Her',4.4062465,0.8019578,4.82];
ybs[6245]=['21 Oph',4.4176751,0.0206422,5.51];
ybs[6246]=['',4.4083254,0.7574072,6.13];
ybs[6247]=['',4.4295076,-0.7519409,5.96];
ybs[6248]=['50 Her',4.4132878,0.5196342,5.72];
ybs[6249]=['',4.4134613,0.5675779,6.13];
ybs[6250]=['',4.4308357,-0.7302186,5.45];
ybs[6251]=['',4.4306287,-0.7335061,6.32];
ybs[6252]=['ζ1 Sco',4.4307165,-0.73992,4.73];
ybs[6253]=['',4.4315659,-0.7309832,6.45];
ybs[6254]=['',4.4124605,0.7306435,6.29];
ybs[6255]=['',4.4321315,-0.7304537,6.59];
ybs[6256]=['',4.4327028,-0.7419524,5.88];
ybs[6257]=['',4.3729589,1.3522178,5.98];
ybs[6258]=['49 Her',4.4201145,0.2607707,6.52];
ybs[6259]=['',4.4271749,-0.3568842,5.88];
ybs[6260]=['51 Her',4.4183245,0.4297543,5.04];
ybs[6261]=['ζ2 Sco',4.4332849,-0.7399005,3.62];
ybs[6262]=['',4.4349147,-0.718774,5.77];
ybs[6263]=['',4.4327435,-0.534403,6.35];
ybs[6264]=['',4.4407093,-0.8849868,6.33];
ybs[6265]=['',4.4422892,-0.9130643,5.94];
ybs[6266]=['',4.45836,-1.2094704,5.79];
ybs[6267]=['',4.4298193,-0.0286984,6.25];
ybs[6268]=['',4.4323248,-0.2063734,6.57];
ybs[6269]=['53 Her',4.4233147,0.5527274,5.32];
ybs[6270]=['23 Oph',4.431797,-0.1079618,5.25];
ybs[6271]=['ι Oph',4.4286883,0.176856,4.38];
ybs[6272]=['',4.4388302,-0.5853552,6.37];
ybs[6273]=['',4.4419965,-0.7130445,6.15];
ybs[6274]=['',4.438431,-0.2938659,6.37];
ybs[6275]=['ζ Ara',4.452004,-0.9777347,3.13];
ybs[6276]=['',4.4238347,0.8270126,6];
ybs[6277]=['',4.4322761,0.3652424,5.41];
ybs[6278]=['27 Sco',4.4441507,-0.5810203,5.48];
ybs[6279]=['',4.4501009,-0.8843773,5.55];
ybs[6280]=['',4.4340577,0.2371578,6.34];
ybs[6281]=['24 Oph',4.4420322,-0.4045809,5.58];
ybs[6282]=['56 Her',4.4325886,0.4485293,6.08];
ybs[6283]=['54 Her',4.4343344,0.3211719,5.35];
ybs[6284]=['',4.4430565,-0.3415724,6.27];
ybs[6285]=['ε1 Ara',4.4559341,-0.9283389,4.06];
ybs[6286]=['',4.4443545,-0.1918786,6.19];
ybs[6287]=['',4.4583365,-0.9534041,5.65];
ybs[6288]=['',4.4517409,-0.6571315,6.09];
ybs[6289]=['κ Oph',4.4446838,0.1630936,3.2];
ybs[6290]=['',4.4593352,-0.8495691,6];
ybs[6291]=['',4.4439372,0.2417921,6.37];
ybs[6292]=['',4.4499836,-0.2600472,6.59];
ybs[6293]=['',4.4594815,-0.7937859,6.65];
ybs[6294]=['',4.4661928,-1.02951,6.11];
ybs[6295]=['57 Her',4.4434296,0.4419566,6.28];
ybs[6296]=['',4.4359013,0.8727974,6.56];
ybs[6297]=['',4.4442907,0.4250043,6.32];
ybs[6298]=['',4.4559017,-0.4384476,5.86];
ybs[6299]=['26 Oph',4.4567631,-0.4366521,5.75];
ybs[6300]=['',4.4592548,-0.6276741,5.97];
ybs[6301]=['',4.4652852,-0.8928954,6.45];
ybs[6302]=['',4.4439757,0.7414521,6.34];
ybs[6303]=['ε2 Ara',4.4714922,-0.9296422,5.29];
ybs[6304]=['19 Dra',4.4337153,1.136267,4.89];
ybs[6305]=['',4.464583,-0.561506,5.03];
ybs[6306]=['',4.4570902,0.1143985,6.59];
ybs[6307]=['30 Oph',4.4599419,-0.0741986,4.82];
ybs[6308]=['20 Dra',4.4354403,1.1346025,6.41];
ybs[6309]=['',4.4774882,-1.0077393,5.73];
ybs[6310]=['29 Oph',4.4639213,-0.33011,6.26];
ybs[6311]=['ε UMi',4.3804969,1.4311789,4.23];
ybs[6312]=['',4.4733829,-0.8235753,6.06];
ybs[6313]=['ε Her',4.4552865,0.5392575,3.92];
ybs[6314]=['',4.4585962,0.394503,5.65];
ybs[6315]=['',4.4614227,0.2604184,6.31];
ybs[6316]=['',4.4734782,-0.6663592,5.91];
ybs[6317]=['',4.4592456,0.4741642,6.55];
ybs[6318]=['',4.4635491,0.1469953,6.33];
ybs[6319]=['',4.4494496,0.9888836,6.03];
ybs[6320]=['',4.4793477,-0.7946203,6.28];
ybs[6321]=['59 Her',4.4609054,0.585379,5.25];
ybs[6322]=['',4.4643403,0.4446633,5.75];
ybs[6323]=['',4.4775323,-0.596024,4.87];
ybs[6324]=['',4.4326348,1.2757721,6.3];
ybs[6325]=['',4.4639516,0.5560003,6.36];
ybs[6326]=['',4.4683468,0.2454655,4.98];
ybs[6327]=['',4.4824042,-0.7702379,6.19];
ybs[6328]=['',4.4685208,0.2527818,6.52];
ybs[6329]=['',4.4766456,-0.3581708,6.3];
ybs[6330]=['',4.4706546,0.2369761,5.93];
ybs[6331]=['',4.4720157,0.2363193,6.08];
ybs[6332]=['',4.471405,0.3431858,6.35];
ybs[6333]=['',4.4843042,-0.6501987,5.98];
ybs[6334]=['',4.4459482,1.2070063,6.4];
ybs[6335]=['61 Her',4.4690869,0.6176108,6.69];
ybs[6336]=['',4.4848006,-0.6191938,6.13];
ybs[6337]=['',4.4573268,1.0580239,6.13];
ybs[6338]=['',4.4781963,0.0117942,6.01];
ybs[6339]=['',4.4829788,-0.3768336,6.3];
ybs[6340]=['',4.4708121,0.6067252,6.04];
ybs[6341]=['',4.4749422,0.3415977,6.17];
ybs[6342]=['',4.4793701,-0.0160318,5.64];
ybs[6343]=['',4.4861911,-0.4631922,6.29];
ybs[6344]=['60 Her',4.478204,0.2219033,4.91];
ybs[6345]=['',4.5028789,-1.0768638,6.39];
ybs[6346]=['',4.5145554,-1.2347176,6.22];
ybs[6347]=['',4.4817273,0.1694192,6.37];
ybs[6348]=['',4.4819495,0.1820005,6.37];
ybs[6349]=['',4.460953,1.1269963,6.1];
ybs[6350]=['',4.4852575,-0.0293625,6.38];
ybs[6351]=['',4.4754981,0.7641974,6.43];
ybs[6352]=['',4.4740392,0.8513208,6.09];
ybs[6353]=['',4.4818843,0.3849826,5.56];
ybs[6354]=['',4.4917812,-0.3077788,5.99];
ybs[6355]=['',4.494676,-0.5310789,5.97];
ybs[6356]=['',4.4911067,-0.0192816,6.06];
ybs[6357]=['',4.5178166,-1.1731973,5.89];
ybs[6358]=['μ Dra',4.4757207,0.9502168,5.83];
ybs[6359]=['μ Dra',4.4757134,0.9502168,5.8];
ybs[6360]=['',4.5038018,-0.7780939,5.08];
ybs[6361]=['',4.4941763,-0.0682031,6.36];
ybs[6362]=['',4.534659,-1.3012108,6.25];
ybs[6363]=['',4.5082327,-0.8534256,5.84];
ybs[6364]=['',4.4982964,-0.1840948,5.56];
ybs[6365]=['',4.4874554,0.7066924,6.34];
ybs[6366]=['',4.4888194,0.6267442,5.39];
ybs[6367]=['η Oph',4.5010065,-0.2748711,2.43];
ybs[6368]=['',4.4551241,1.3136802,6.21];
ybs[6369]=['η Sco',4.5100374,-0.7550724,3.33];
ybs[6370]=['',4.5103207,-0.6899321,5.67];
ybs[6371]=['',4.5103037,-0.6779814,6.3];
ybs[6372]=['',4.4889406,0.8869207,6.46];
ybs[6373]=['',4.5202175,-0.9932766,6.09];
ybs[6374]=['',4.5017079,0.2171735,6.57];
ybs[6375]=['',4.5094433,-0.441185,6.54];
ybs[6376]=['',4.5103791,-0.4849424,6.14];
ybs[6377]=['',4.4951892,0.7112647,5.08];
ybs[6378]=['',4.5130435,-0.5665607,6.01];
ybs[6379]=['',4.5061878,0.1373769,6.33];
ybs[6380]=['63 Her',4.5025271,0.4226104,6.19];
ybs[6381]=['',4.5199016,-0.6944512,6.6];
ybs[6382]=['37 Oph',4.5091915,0.1843418,5.33];
ybs[6383]=['',4.5114783,0.0057406,6.65];
ybs[6384]=['',4.4985058,0.9142829,6.29];
ybs[6385]=['ζ Dra',4.4892106,1.1464972,3.17];
ybs[6386]=['',4.5233246,-0.585909,5.53];
ybs[6387]=['',4.5247962,-0.6739679,5.96];
ybs[6388]=['',4.5037939,0.8678284,6.04];
ybs[6389]=['',4.5487827,-1.2228599,6.53];
ybs[6390]=['36 Oph',4.5231308,-0.4646863,5.11];
ybs[6391]=['36 Oph',4.5231162,-0.464662,5.07];
ybs[6392]=['',4.5255085,-0.5276494,6.21];
ybs[6393]=['',4.5226092,-0.2549176,5.99];
ybs[6394]=['',4.5279601,-0.6243167,6.12];
ybs[6395]=['α1 Her',4.5185885,0.2507701,3.48];
ybs[6396]=['α2 Her',4.5186103,0.2507653,5.39];
ybs[6397]=['',4.5423368,-1.0422095,5.91];
ybs[6398]=['',4.5308755,-0.5704383,5.55];
ybs[6399]=['δ Her',4.5198578,0.4331405,3.14];
ybs[6400]=['ι Aps',4.5570379,-1.2242011,5.41];
ybs[6401]=['',4.5259455,0.0377812,6.17];
ybs[6402]=['',4.5283143,-0.1093651,6.09];
ybs[6403]=['',4.5272515,0.0207571,5.88];
ybs[6404]=['41 Oph',4.5276717,-0.008142,4.73];
ybs[6405]=['',4.5403871,-0.8142624,5.48];
ybs[6406]=['ζ Aps',4.5559365,-1.1831389,4.78];
ybs[6407]=['π Her',4.5193656,0.642056,3.16];
ybs[6408]=['',4.5227894,0.4140105,5.96];
ybs[6409]=['',4.5391075,-0.7705588,5.76];
ybs[6410]=['',4.506101,1.0969573,5.56];
ybs[6411]=['',4.5364492,-0.5685172,6.36];
ybs[6412]=['',4.5426251,-0.8741133,6.27];
ybs[6413]=['ο Oph',4.5346494,-0.4242446,5.2];
ybs[6414]=['ο Oph',4.5346348,-0.4241963,6.8];
ybs[6415]=['',4.5392467,-0.6110349,5.91];
ybs[6416]=['',4.5417689,-0.7721825,6.65];
ybs[6417]=['',4.5356841,-0.2850523,6.43];
ybs[6418]=['',4.6047277,-1.4114895,5.88];
ybs[6419]=['',4.5311266,0.4026483,6.45];
ybs[6420]=['68 Her',4.5294927,0.5773385,4.82];
ybs[6421]=['',4.5334669,0.3018987,6];
ybs[6422]=['',4.5360287,0.1892667,5.03];
ybs[6423]=['',4.5373449,0.1058569,6.51];
ybs[6424]=['',4.542582,-0.3102491,6.02];
ybs[6425]=['69 Her',4.5307902,0.6504996,4.65];
ybs[6426]=['',4.5262154,0.8669026,7.48];
ybs[6427]=['',4.5583173,-1.0127836,5.88];
ybs[6428]=['',4.5426079,-0.1036163,6.32];
ybs[6429]=['',4.5638165,-1.0974899,5.7];
ybs[6430]=['',4.5456167,-0.3377564,6.52];
ybs[6431]=['',4.5590085,-0.986864,5.8];
ybs[6432]=['',4.5361625,0.5027045,5.65];
ybs[6433]=['',4.5338325,0.6770298,5.94];
ybs[6434]=['ξ Oph',4.5475771,-0.3688196,4.39];
ybs[6435]=['ν Ser',4.5465029,-0.2245554,4.33];
ybs[6436]=['',4.56476,-1.0592553,5.77];
ybs[6437]=['',4.5236448,1.0585262,6.32];
ybs[6438]=['',4.5466442,-0.1870159,6.46];
ybs[6439]=['',4.5555586,-0.6601388,6.41];
ybs[6440]=['ι Ara',4.5588466,-0.82879,5.25];
ybs[6441]=['',4.5431754,0.3148187,5];
ybs[6442]=['θ Oph',4.5521276,-0.4366459,3.27];
ybs[6443]=['',4.555355,-0.6270651,6.47];
ybs[6444]=['',4.5422159,0.4453726,5.38];
ybs[6445]=['',4.5566541,-0.6499361,5.93];
ybs[6446]=['70 Her',4.5454876,0.4272615,5.12];
ybs[6447]=['72 Her',4.5440715,0.5663326,5.39];
ybs[6448]=['43 Oph',4.558154,-0.4915006,5.35];
ybs[6449]=['',4.5627517,-0.7710844,5.12];
ybs[6450]=['β Ara',4.568437,-0.9694744,2.85];
ybs[6451]=['γ Ara',4.5689374,-0.9842652,3.34];
ybs[6452]=['',4.5486425,0.2916795,6.35];
ybs[6453]=['74 Her',4.5419261,0.8067144,5.59];
ybs[6454]=['',4.5549801,-0.0420011,6.29];
ybs[6455]=['',4.547995,0.5015934,6.35];
ybs[6456]=['',4.5426902,0.8407063,6.43];
ybs[6457]=['κ Ara',4.5710236,-0.8840107,5.23];
ybs[6458]=['',4.548328,0.6973575,5.51];
ybs[6459]=['',4.5658176,-0.605863,6.16];
ybs[6460]=['',4.5817883,-1.1004603,6.24];
ybs[6461]=['',4.5637218,-0.374528,5.85];
ybs[6462]=['',4.5632488,-0.3222416,6.21];
ybs[6463]=['',4.5655953,-0.4234275,6.19];
ybs[6464]=['',4.5752351,-0.9069634,6.19];
ybs[6465]=['',4.559413,0.1542024,5.77];
ybs[6466]=['',4.5743972,-0.8004025,5.29];
ybs[6467]=['',4.5762902,-0.8839422,5.92];
ybs[6468]=['',4.547459,0.9320357,5.67];
ybs[6469]=['73 Her',4.5595403,0.4004254,5.74];
ybs[6470]=['',4.5616091,0.2842049,5.71];
ybs[6471]=['',4.5618022,0.2720752,6.35];
ybs[6472]=['',4.5797264,-0.9130292,5.75];
ybs[6473]=['ρ Her',4.5570002,0.6480203,5.47];
ybs[6474]=['ρ Her',4.5570221,0.6480058,4.52];
ybs[6475]=['44 Oph',4.5711201,-0.4222241,4.17];
ybs[6476]=['',4.5829997,-0.9631576,5.94];
ybs[6477]=['',4.5584842,0.6730882,6.49];
ybs[6478]=['',4.5685248,-0.0291172,6.44];
ybs[6479]=['',4.5735945,-0.4530777,6.44];
ybs[6480]=['',4.5603842,0.6446284,6.28];
ybs[6481]=['45 Oph',4.5756659,-0.5215536,4.29];
ybs[6482]=['',4.5715524,-0.0890634,4.54];
ybs[6483]=['',4.5768376,-0.5190642,6];
ybs[6484]=['',4.5676142,0.2949748,5.98];
ybs[6485]=['',4.5735744,-0.2186649,6.21];
ybs[6486]=['',4.5697396,0.1322801,6.06];
ybs[6487]=['σ Oph',4.5707251,0.0719759,4.34];
ybs[6488]=['',4.5676874,0.4688343,6.41];
ybs[6489]=['δ Ara',4.594388,-1.0593762,3.62];
ybs[6490]=['',4.5829031,-0.6421662,6.02];
ybs[6491]=['',4.5714706,0.3501931,5.54];
ybs[6492]=['',4.5851517,-0.6725015,6.39];
ybs[6493]=['',4.5777993,-0.1435345,6.37];
ybs[6494]=['',4.5951461,-0.9936916,5.95];
ybs[6495]=['',4.5706157,0.6052719,5.94];
ybs[6496]=['',4.5809502,0.0055037,5.44];
ybs[6497]=['υ Sco',4.5909103,-0.6511825,2.69];
ybs[6498]=['77 Her',4.5696259,0.8420101,5.85];
ybs[6499]=['α Ara',4.5964568,-0.8707397,2.95];
ybs[6500]=['',4.5638432,1.0477451,5.65];
ybs[6501]=['',4.5853625,-0.1035757,6.37];
ybs[6502]=['',4.5960588,-0.8037244,6.03];
ybs[6503]=['',4.5657433,1.0233771,6.51];
ybs[6504]=['',4.5878414,-0.0187961,5.31];
ybs[6505]=['',4.5952034,-0.5884631,6.44];
ybs[6506]=['',4.5595658,1.174415,6.43];
ybs[6507]=['51 Oph',4.5931342,-0.4184718,4.81];
ybs[6508]=['',4.5946423,-0.4587326,6.05];
ybs[6509]=['',4.5872857,0.2078778,6.39];
ybs[6510]=['',4.5979308,-0.5985274,6.17];
ybs[6511]=['',4.6014331,-0.7188419,5.84];
ybs[6512]=['',4.5918922,0.0473066,5.59];
ybs[6513]=['',4.61384,-1.0447156,6.28];
ybs[6514]=['λ Her',4.5883271,0.4554651,4.41];
ybs[6515]=['λ Sco',4.6033156,-0.6478079,1.63];
ybs[6516]=['',4.5889079,0.5435668,5.61];
ybs[6517]=['',4.5295076,1.3982881,5.72];
ybs[6518]=['',4.6120084,-0.9313884,6.1];
ybs[6519]=['',4.5874114,0.6783715,6.43];
ybs[6520]=['',4.5954736,0.2079861,6.42];
ybs[6521]=['78 Her',4.5929749,0.4955636,5.62];
ybs[6522]=['',4.6015374,-0.1004894,5.62];
ybs[6523]=['',4.6078837,-0.568871,5.7];
ybs[6524]=['β Dra',4.5854089,0.9125771,2.79];
ybs[6525]=['σ Ara',4.6128589,-0.8118796,4.59];
ybs[6526]=['',4.5935332,0.5978995,6.56];
ybs[6527]=['',4.6125348,-0.6536559,6.48];
ybs[6528]=['',4.5861129,1.0098812,6.4];
ybs[6529]=['',4.6001307,0.3358652,5.64];
ybs[6530]=['',4.6014498,0.2845696,5.69];
ybs[6531]=['',4.601752,0.2588119,6.48];
ybs[6532]=['',4.6072952,-0.1964228,5.55];
ybs[6533]=['52 Oph',4.6100381,-0.3849473,6.57];
ybs[6534]=['',4.6162236,-0.6745102,4.29];
ybs[6535]=['',4.6209637,-0.8739006,5.93];
ybs[6536]=['53 Oph',4.6058617,0.1671029,5.81];
ybs[6537]=['π Ara',4.6241689,-0.9513921,5.25];
ybs[6538]=['',4.5979602,0.7196065,5.74];
ybs[6539]=['',4.6049186,0.2878296,6.4];
ybs[6540]=['',4.7476005,-1.4820516,6.45];
ybs[6541]=['θ Sco',4.6198624,-0.7506432,1.87];
ybs[6542]=['ν1 Dra',4.5927284,0.9629059,4.88];
ybs[6543]=['ν2 Dra',4.5931222,0.9627128,4.87];
ybs[6544]=['α Oph',4.6071657,0.2190001,2.08];
ybs[6545]=['',4.6201191,-0.6645589,6.26];
ybs[6546]=['',4.6234407,-0.7485902,6.1];
ybs[6547]=['',4.6114643,0.3662468,6.1];
ybs[6548]=['',4.5983565,1.004364,6.17];
ybs[6549]=['ξ Ser',4.6197232,-0.2689459,3.54];
ybs[6550]=['',4.6198023,-0.2719565,5.94];
ybs[6551]=['',4.6094648,0.6508292,6.1];
ybs[6552]=['',4.6117666,0.4917126,6.38];
ybs[6553]=['',4.6547938,-1.2606167,6.49];
ybs[6554]=['27 Dra',4.5897227,1.1889367,5.05];
ybs[6555]=['μ Oph',4.6205871,-0.1418887,4.62];
ybs[6556]=['',4.6220542,-0.1908861,5.75];
ybs[6557]=['λ Ara',4.6337411,-0.8626275,4.77];
ybs[6558]=['',4.613742,0.5371049,6.02];
ybs[6559]=['79 Her',4.6179832,0.4240981,5.77];
ybs[6560]=['',4.6373797,-0.8190983,5.79];
ybs[6561]=['26 Dra',4.6041502,1.0797064,5.23];
ybs[6562]=['82 Her',4.6127417,0.8477774,5.37];
ybs[6563]=['',4.6258917,0.0352196,6.26];
ybs[6564]=['',4.6411282,-0.8817243,6.24];
ybs[6565]=['',4.6247138,0.2324593,6.12];
ybs[6566]=['',4.6306398,-0.0377355,6.19];
ybs[6567]=['',4.6233194,0.5712306,6.37];
ybs[6568]=['κ Sco',4.642175,-0.6813479,2.41];
ybs[6569]=['ο Ser',4.6363351,-0.2248724,4.26];
ybs[6570]=['η Pav',4.6589354,-1.1297602,3.62];
ybs[6571]=['',4.6436369,-0.6449693,5.54];
ybs[6572]=['',4.6283233,0.5444158,6.03];
ybs[6573]=['μ Ara',4.6503222,-0.9048075,5.15];
ybs[6574]=['',4.6543447,-1.0044778,6.01];
ybs[6575]=['',4.6445918,-0.5769915,6.4];
ybs[6576]=['ι Her',4.6253178,0.8027873,3.8];
ybs[6577]=['',4.6343318,0.2647477,6.34];
ybs[6578]=['',4.6361997,0.1100228,5.95];
ybs[6579]=['',4.631497,0.5459056,6.28];
ybs[6580]=['',4.633565,0.4276779,6.36];
ybs[6581]=['',4.6451429,-0.4868099,6.36];
ybs[6582]=['',4.6377715,0.2782615,5.52];
ybs[6583]=['58 Oph',4.6454634,-0.3785842,4.87];
ybs[6584]=['ω Dra',4.6113093,1.199854,4.8];
ybs[6585]=['',4.652051,-0.7458863,5.87];
ybs[6586]=['',4.6098225,1.2140369,6.42];
ybs[6587]=['',4.6305598,0.7585437,6.59];
ybs[6588]=['',4.6468088,-0.2359053,6.39];
ybs[6589]=['',4.6464618,-0.1236956,6.3];
ybs[6590]=['83 Her',4.6395912,0.4285771,5.52];
ybs[6591]=['β Oph',4.6447036,0.079574,2.77];
ybs[6592]=['',4.6438907,0.2493544,6.24];
ybs[6593]=['',4.629247,1.0000861,6.77];
ybs[6594]=['',4.6109903,1.2643929,5.86];
ybs[6595]=['',4.6331313,0.9042359,5.99];
ybs[6596]=['84 Her',4.6434696,0.424459,5.71];
ybs[6597]=['61 Oph',4.6495434,0.0448902,6.17];
ybs[6598]=['',4.6496452,0.0448807,6.56];
ybs[6599]=['',4.6479078,0.2513744,6.19];
ybs[6600]=['',4.6412802,0.7692745,6.34];
ybs[6601]=['',4.6623416,-0.6652799,6.43];
ybs[6602]=['',4.6702691,-0.9670328,6.11];
ybs[6603]=['ι1 Sco',4.6644797,-0.7004489,3.03];
ybs[6604]=['3 Sgr',4.663744,-0.4858422,4.54];
ybs[6605]=['',4.6643896,-0.3924171,6.18];
ybs[6606]=['',4.6443944,0.938879,5.75];
ybs[6607]=['',4.6532358,0.54974,6.23];
ybs[6608]=['',4.6634494,-0.2571168,5.94];
ybs[6609]=['',4.6676481,-0.4708973,6.35];
ybs[6610]=['',4.6780922,-0.9357854,5.92];
ybs[6611]=['μ Her',4.6568419,0.4837009,3.42];
ybs[6612]=['',4.6837567,-1.0501326,5.78];
ybs[6613]=['',4.6538256,0.6784887,6.52];
ybs[6614]=['',4.6541493,0.6861882,6.68];
ybs[6615]=['',4.6601975,0.308767,5.72];
ybs[6616]=['',4.6709698,-0.5534159,4.83];
ybs[6617]=['γ Oph',4.6640542,0.0471492,3.75];
ybs[6618]=['',4.6742243,-0.6466102,3.21];
ybs[6619]=['ι2 Sco',4.6758303,-0.6997912,4.81];
ybs[6620]=['',4.681158,-0.9273725,6.09];
ybs[6621]=['',4.6659503,0.0662983,6.22];
ybs[6622]=['',4.6921359,-1.1430509,6.49];
ybs[6623]=['',4.7150193,-1.3295556,6.07];
ybs[6624]=['ψ1 Dra',4.6320006,1.2590773,4.58];
ybs[6625]=['ψ1 Dra',4.6321212,1.2592182,5.79];
ybs[6626]=['',4.6656707,0.3588399,5.69];
ybs[6627]=['',4.6702833,0.0341395,6.47];
ybs[6628]=['',4.6829728,-0.795945,6.11];
ybs[6629]=['',4.6586604,0.8308806,6.43];
ybs[6630]=['',4.667402,0.3359745,6.12];
ybs[6631]=['',4.6818057,-0.7116816,5.96];
ybs[6632]=['87 Her',4.6672373,0.4471083,5.12];
ybs[6633]=['',4.6797879,-0.5333948,6.66];
ybs[6634]=['',4.7540285,-1.4221407,6.35];
ybs[6635]=['',4.6844397,-0.6074216,5.9];
ybs[6636]=['',4.6848636,-0.6007449,5.84];
ybs[6637]=['',4.687692,-0.7330356,6.2];
ybs[6638]=['',4.6760826,0.2084323,6.17];
ybs[6639]=['',4.6869936,-0.595461,6.06];
ybs[6640]=['',4.6875273,-0.6112456,6.45];
ybs[6641]=['',4.6877043,-0.6218142,6.03];
ybs[6642]=['',4.6738877,0.5116892,5.5];
ybs[6643]=['',4.676043,0.3894184,5.98];
ybs[6644]=['30 Dra',4.6668517,0.8862048,5.02];
ybs[6645]=['',4.689236,-0.6062146,6.17];
ybs[6646]=['',4.6895135,-0.609089,5.6];
ybs[6647]=['',4.6820735,-0.0216489,6.35];
ybs[6648]=['',4.6911224,-0.6071756,6.38];
ybs[6649]=['',4.6851028,-0.1072853,6.21];
ybs[6650]=['',4.6918043,-0.6065924,5.96];
ybs[6651]=['',4.6920411,-0.6079738,6.42];
ybs[6652]=['88 Her',4.6713392,0.8445534,6.68];
ybs[6653]=['',4.6814122,0.2674206,6.46];
ybs[6654]=['',4.6870754,-0.1902914,6.18];
ybs[6655]=['',4.6845959,0.0227166,5.95];
ybs[6656]=['',4.6941386,-0.6015943,5.96];
ybs[6657]=['',4.6770472,0.6993238,6.46];
ybs[6658]=['',4.68725,0.1064349,5.77];
ybs[6659]=['',4.6972038,-0.6366597,6.06];
ybs[6660]=['',4.6956444,-0.4344028,6.2];
ybs[6661]=['',4.6807393,0.6977555,6.04];
ybs[6662]=['',4.6800313,0.8140128,6.38];
ybs[6663]=['',4.7049138,-0.7739392,4.86];
ybs[6664]=['',4.6914344,0.1942189,6.38];
ybs[6665]=['90 Her',4.6860753,0.698217,5.16];
ybs[6666]=['',4.7052676,-0.7034851,6.43];
ybs[6667]=['',4.6998347,-0.3281909,6.52];
ybs[6668]=['',4.7036291,-0.4898546,5.8];
ybs[6669]=['',4.7014659,-0.2760071,5.89];
ybs[6670]=['',4.7091281,-0.7280963,4.88];
ybs[6671]=['',4.7097077,-0.68308,6.29];
ybs[6672]=['',4.700835,0.011671,5.82];
ybs[6673]=['89 Her',4.696016,0.4546221,5.46];
ybs[6674]=['',4.7031311,-0.0712665,5.47];
ybs[6675]=['',4.6980276,0.3920413,5.58];
ybs[6676]=['ξ Dra',4.6856843,0.9925625,3.75];
ybs[6677]=['',4.7041939,0.0011378,5.97];
ybs[6678]=['',4.7033661,0.1132108,6.29];
ybs[6679]=['',4.7137835,-0.6433026,5.74];
ybs[6680]=['',4.7121835,-0.5019482,6.01];
ybs[6681]=['',4.71415,-0.5280176,5.16];
ybs[6682]=['',4.7141719,-0.5280225,7.04];
ybs[6683]=['θ Her',4.699111,0.6501153,3.86];
ybs[6684]=['',4.7054509,0.1927391,6.36];
ybs[6685]=['',4.7040367,0.418786,6.3];
ybs[6686]=['ν Oph',4.7130696,-0.1705853,3.34];
ybs[6687]=['',4.693931,0.9768466,6.1];
ybs[6688]=['4 Sgr',4.7169508,-0.415666,4.76];
ybs[6689]=['35 Dra',4.6623496,1.3431583,5.04];
ybs[6690]=['',4.7010236,0.7914962,6.02];
ybs[6691]=['ξ Her',4.7061061,0.510454,3.7];
ybs[6692]=['',4.7177191,-0.3549801,6.21];
ybs[6693]=['γ Dra',4.6996643,0.8986232,2.23];
ybs[6694]=['',4.7154471,-0.0841477,5.87];
ybs[6695]=['ν Her',4.7092841,0.5268956,4.41];
ybs[6696]=['',4.7263165,-0.6348903,6.3];
ybs[6697]=['',4.7180787,0.0109927,6.37];
ybs[6698]=['ζ Ser',4.7192108,-0.0643986,4.62];
ybs[6699]=['',4.7098604,0.6333331,6];
ybs[6700]=['66 Oph',4.7179703,0.0762533,4.64];
ybs[6701]=['93 Her',4.7166211,0.2923617,4.67];
ybs[6702]=['67 Oph',4.7196857,0.0511773,3.97];
ybs[6703]=['6 Sgr',4.7236316,-0.299428,6.28];
ybs[6704]=['',4.7261293,-0.3975737,5.77];
ybs[6705]=['',4.6644915,1.3666188,6.24];
ybs[6706]=['',4.7100249,0.7937006,6.48];
ybs[6707]=['',4.7206016,0.1094151,6.34];
ybs[6708]=['',4.7182822,0.3404489,6.5];
ybs[6709]=['χ Oct',5.0017118,-1.5271557,5.28];
ybs[6710]=['',4.7205959,0.2634407,6.26];
ybs[6711]=['68 Oph',4.7245783,0.0228012,4.45];
ybs[6712]=['7 Sgr',4.7303146,-0.4237745,5.34];
ybs[6713]=['ψ2 Dra',4.6898409,1.2566809,5.45];
ybs[6714]=['',4.7183164,0.5797002,5.99];
ybs[6715]=['',4.7310209,-0.396478,6.74];
ybs[6716]=['',4.7146699,0.7941511,5.67];
ybs[6717]=['95 Her',4.7227227,0.3769255,5.18];
ybs[6718]=['95 Her',4.7227591,0.3769304,4.96];
ybs[6719]=['',4.7740229,-1.3244443,5.86];
ybs[6720]=['',4.7292539,-0.0934966,6.76];
ybs[6721]=['τ Oph',4.7307133,-0.1427412,5.94];
ybs[6722]=['τ Oph',4.730706,-0.142746,5.24];
ybs[6723]=['',4.6852074,1.3119273,6.36];
ybs[6724]=['9 Sgr',4.7347757,-0.4251328,5.97];
ybs[6725]=['',4.7226384,0.5814106,6.15];
ybs[6726]=['96 Her',4.726622,0.3636396,5.28];
ybs[6727]=['',4.7395334,-0.6265493,6];
ybs[6728]=['',4.7552112,-1.1265335,6.41];
ybs[6729]=['97 Her',4.7270467,0.4001083,6.21];
ybs[6730]=['γ1 Sgr',4.7400148,-0.5162191,4.69];
ybs[6731]=['θ Ara',4.7482898,-0.8741999,3.66];
ybs[6732]=['',4.7304227,0.3423445,6.5];
ybs[6733]=['π Pav',4.7584316,-1.1111389,4.35];
ybs[6734]=['γ2 Sgr',4.7434863,-0.5309457,2.99];
ybs[6735]=['',4.7370725,0.0335405,6.14];
ybs[6736]=['',4.7463243,-0.6286011,5.95];
ybs[6737]=['',4.7486591,-0.7578386,5.77];
ybs[6738]=['',4.7486591,-0.7578386,5.77];
ybs[6739]=['',4.7786277,-1.2856926,5.85];
ybs[6740]=['70 Oph',4.7406884,0.0436754,4.03];
ybs[6741]=['',4.7284694,0.8458938,6.21];
ybs[6742]=['',4.7364658,0.4179198,6.34];
ybs[6743]=['',4.7439828,-0.1452212,5.85];
ybs[6744]=['',4.7444242,-0.0828683,5.77];
ybs[6745]=['',4.743707,-0.007738,6.34];
ybs[6746]=['',4.7415155,0.2095613,7.04];
ybs[6747]=['',4.7561131,-0.7987082,6.15];
ybs[6748]=['',4.7637484,-1.030348,6.38];
ybs[6749]=['ι Pav',4.7662463,-1.082044,5.49];
ybs[6750]=['',4.7491276,-0.3741987,6.28];
ybs[6751]=['',4.7401959,0.3778621,6.15];
ybs[6752]=['',4.7358922,0.6996446,6.52];
ybs[6753]=['98 Her',4.7424781,0.387849,5.06];
ybs[6754]=['',4.7533224,-0.4965963,4.57];
ybs[6755]=['',4.7370585,0.7321538,6.34];
ybs[6756]=['',4.7411416,0.5625833,5.71];
ybs[6757]=['',4.7516508,-0.2993237,5.52];
ybs[6758]=['71 Oph',4.7485488,0.1525029,4.64];
ybs[6759]=['72 Oph',4.7487083,0.1669895,3.73];
ybs[6760]=['',4.759354,-0.6399684,6.58];
ybs[6761]=['',4.7567622,-0.4444961,6.61];
ybs[6762]=['',4.7855026,-1.2347094,6.73];
ybs[6763]=['99 Her',4.74645,0.533471,5.04];
ybs[6764]=['',4.7505742,0.2282058,6.63];
ybs[6765]=['',4.7618776,-0.5709742,6.43];
ybs[6766]=['',4.7674689,-0.8291564,6.07];
ybs[6767]=['ο Her',4.7487946,0.5020693,3.83];
ybs[6768]=['',4.7622038,-0.536222,5.53];
ybs[6769]=['100 Her',4.7501411,0.4556268,5.86];
ybs[6770]=['100 Her',4.7501413,0.4555589,5.9];
ybs[6771]=['ε Tel',4.7680299,-0.8019522,4.53];
ybs[6772]=['',4.7538239,0.2493939,6.37];
ybs[6773]=['',4.7598931,-0.2431124,6.39];
ybs[6774]=['',4.7671283,-0.7217509,5.86];
ybs[6775]=['102 Her',4.7544321,0.3633605,4.36];
ybs[6776]=['',4.7659526,-0.5898156,6.16];
ybs[6777]=['δ UMi',4.5627235,1.5082677,4.36];
ybs[6778]=['',4.7446006,0.8870868,6.29];
ybs[6779]=['',4.7477243,0.7586169,5];
ybs[6780]=['',4.7456262,0.8676769,6.32];
ybs[6781]=['',4.7505723,0.635397,5.48];
ybs[6782]=['101 Her',4.7550008,0.3499371,5.1];
ybs[6783]=['73 Oph',4.7585623,0.0697843,5.73];
ybs[6784]=['',4.783257,-1.111458,6.47];
ybs[6785]=['',4.7600618,0.0545399,5.69];
ybs[6786]=['',4.766764,-0.3462091,6.36];
ybs[6787]=['',4.7558141,0.5318752,6.38];
ybs[6788]=['',4.7634217,0.0581149,5.51];
ybs[6789]=['11 Sgr',4.7689992,-0.4135551,4.98];
ybs[6790]=['',4.7703003,-0.5043148,6.51];
ybs[6791]=['',4.7606507,0.2876642,6.09];
ybs[6792]=['',4.7763617,-0.7213302,5.47];
ybs[6793]=['',4.7893054,-1.1003825,5.6];
ybs[6794]=['',4.7573898,0.6712966,6.4];
ybs[6795]=['',4.7590632,0.6365484,5.58];
ybs[6796]=['',4.8005477,-1.1906577,6.33];
ybs[6797]=['40 Dra',4.7059231,1.3962717,6.04];
ybs[6798]=['41 Dra',4.706341,1.3963307,5.68];
ybs[6799]=['24 UMi',4.5524428,1.5154693,5.79];
ybs[6800]=['μ Sgr',4.7777927,-0.367422,3.86];
ybs[6801]=['',4.7745632,-0.0698977,6.59];
ybs[6802]=['',4.7669356,0.5838646,5.88];
ybs[6803]=['104 Her',4.767697,0.5482323,4.97];
ybs[6804]=['14 Sgr',4.7800076,-0.3788351,5.44];
ybs[6805]=['',4.7601367,0.9475744,5.95];
ybs[6806]=['',4.7881879,-0.7714075,5.46];
ybs[6807]=['',4.7946415,-0.9776354,5.33];
ybs[6808]=['',4.7740979,0.382002,6.12];
ybs[6809]=['',4.7936733,-0.8911559,6.06];
ybs[6810]=['15 Sgr',4.7841132,-0.3616403,5.38];
ybs[6811]=['16 Sgr',4.7840997,-0.3557013,5.95];
ybs[6812]=['',4.7706933,0.7182629,6.36];
ybs[6813]=['',4.7853341,-0.3255627,6.07];
ybs[6814]=['',4.7724481,0.6768438,6.04];
ybs[6815]=['',4.7619824,1.0544414,6.49];
ybs[6816]=['',4.8068462,-1.1148582,6.18];
ybs[6817]=['φ Oct',4.8273647,-1.3095511,5.47];
ybs[6818]=['',4.7867666,-0.062994,6.36];
ybs[6819]=['',4.7801484,0.5098936,6.56];
ybs[6820]=['η Sgr',4.7953673,-0.641453,3.11];
ybs[6821]=['',4.7951196,-0.5951246,6.16];
ybs[6822]=['',4.7871112,0.0416443,6.01];
ybs[6823]=['',4.7939801,-0.4999238,6.19];
ybs[6824]=['',4.793942,-0.4935824,6.4];
ybs[6825]=['',4.8563834,-1.4000561,5.95];
ybs[6826]=['',4.7926146,-0.3030771,5.75];
ybs[6827]=['',4.8002271,-0.7379023,6.3];
ybs[6828]=['',4.7907523,-0.0523346,6];
ybs[6829]=['',4.7938854,-0.3220891,6.54];
ybs[6830]=['',4.7967592,-0.4718184,4.65];
ybs[6831]=['',4.7932511,-0.1701639,6.31];
ybs[6832]=['',4.7914639,0.0177079,6.63];
ybs[6833]=['',4.7834454,0.7359598,5.59];
ybs[6834]=['',4.7994967,-0.4467189,6.51];
ybs[6835]=['',4.782792,0.7891913,6.29];
ybs[6836]=['',4.7993307,-0.3248028,6.84];
ybs[6837]=['',4.7780352,0.9877818,6.37];
ybs[6838]=['36 Dra',4.7734591,1.1240647,5.03];
ybs[6839]=['',4.7952371,0.2406138,6.3];
ybs[6840]=['',4.7954229,0.3166137,5.99];
ybs[6841]=['',4.7899215,0.714631,6.11];
ybs[6842]=['',4.7952181,0.4067646,6.63];
ybs[6843]=['ξ Pav',4.8219653,-1.0730611,4.36];
ybs[6844]=['',4.8097821,-0.6540927,6.45];
ybs[6845]=['',4.8003145,0.1268767,5.39];
ybs[6846]=['',4.8054422,-0.2761347,5.39];
ybs[6847]=['δ Sgr',4.8097181,-0.5204099,2.7];
ybs[6848]=['105 Her',4.7997522,0.4268352,5.27];
ybs[6849]=['',4.8118073,-0.4346565,6.25];
ybs[6850]=['',4.8159042,-0.6744913,5.1];
ybs[6851]=['',4.8109601,-0.3289783,5.75];
ybs[6852]=['',4.8140539,-0.4960005,6.16];
ybs[6853]=['37 Dra',4.778546,1.2001479,5.95];
ybs[6854]=['74 Oph',4.8079115,0.0591292,4.86];
ybs[6855]=['',4.8025349,0.5179472,5.99];
ybs[6856]=['106 Her',4.8047415,0.3834784,4.95];
ybs[6857]=['η Ser',4.8100561,-0.0504056,3.26];
ybs[6858]=['',4.818303,-0.6397979,5.34];
ybs[6859]=['',4.832233,-1.0997009,6.14];
ybs[6860]=['κ Lyr',4.802191,0.6296188,4.33];
ybs[6861]=['',4.8104869,0.0950639,6.13];
ybs[6862]=['',4.820883,-0.6322685,5.55];
ybs[6863]=['',4.8249358,-0.7696523,5.25];
ybs[6864]=['108 Her',4.8072535,0.5213212,5.63];
ybs[6865]=['107 Her',4.8075822,0.5040625,5.12];
ybs[6866]=['',4.8178468,-0.1781436,6.33];
ybs[6867]=['ε Sgr',4.823797,-0.599911,1.85];
ybs[6868]=['',4.8014748,0.8963629,6.3];
ybs[6869]=['',4.8186385,-0.2094902,5.73];
ybs[6870]=['',4.8127459,0.4066007,5.41];
ybs[6871]=['',4.8151111,0.2101438,5.89];
ybs[6872]=['ζ Sct',4.8205135,-0.1557205,4.68];
ybs[6873]=['',4.8158925,0.3113359,5.25];
ybs[6874]=['',4.8067784,0.8680602,6.4];
ybs[6875]=['',4.8169465,0.2914654,6.22];
ybs[6876]=['18 Sgr',4.8273432,-0.5365823,5.6];
ybs[6877]=['',4.8290716,-0.6279472,6.15];
ybs[6878]=['',4.8220715,-0.0623276,6.38];
ybs[6879]=['',4.8086899,0.8575241,5.05];
ybs[6880]=['',4.8250003,-0.1232678,6.31];
ybs[6881]=['',4.8313586,-0.5922263,6.3];
ybs[6882]=['',4.8365442,-0.8395591,5.46];
ybs[6883]=['109 Her',4.8195856,0.3801627,3.84];
ybs[6884]=['21 Sgr',4.828325,-0.3582945,4.81];
ybs[6885]=['α Tel',4.8367116,-0.8020583,3.51];
ybs[6886]=['',4.8258927,-0.0273455,6.15];
ybs[6887]=['',4.8674713,-1.2906461,5.89];
ybs[6888]=['',4.8265201,0.0889676,6.74];
ybs[6889]=['',4.8199147,0.6763368,6.36];
ybs[6890]=['',4.8285987,0.1404106,5.65];
ybs[6891]=['μ Lyr',4.8210634,0.6897442,5.12];
ybs[6892]=['',4.8249212,0.478358,6.27];
ybs[6893]=['ζ Tel',4.8450528,-0.8561909,4.13];
ybs[6894]=['',4.8295561,0.2614464,6.37];
ybs[6895]=['',4.8395153,-0.5201473,5.92];
ybs[6896]=['',4.8507334,-1.0036996,5.76];
ybs[6897]=['',4.8389594,-0.4646177,6.31];
ybs[6898]=['',4.8427166,-0.680348,5.64];
ybs[6899]=['',4.8181476,0.9304834,6.32];
ybs[6900]=['',4.9148186,-1.4274323,6.27];
ybs[6901]=['λ Sgr',4.8399546,-0.4434438,2.81];
ybs[6902]=['',4.8405935,-0.4667524,6.27];
ybs[6903]=['',4.8463444,-0.7649946,6.36];
ybs[6904]=['ν Pav',4.8576171,-1.0866819,4.64];
ybs[6905]=['',4.8291977,0.520841,5.83];
ybs[6906]=['59 Ser',4.8356851,0.0036633,5.21];
ybs[6907]=['',4.8395263,-0.3104212,6.2];
ybs[6908]=['φ Dra',4.8016731,1.245258,4.22];
ybs[6909]=['',4.8463219,-0.6778199,6.63];
ybs[6910]=['',4.8497004,-0.8238833,5.7];
ybs[6911]=['39 Dra',4.8180214,1.026472,4.98];
ybs[6912]=['',4.8324098,0.4618648,6.53];
ybs[6913]=['',4.8383113,0.0656714,6.07];
ybs[6914]=['',4.8443088,-0.4636811,6.5];
ybs[6915]=['χ Dra',4.8024666,1.2696074,3.57];
ybs[6916]=['',4.8388445,0.1083555,5.73];
ybs[6917]=['',4.846041,-0.4405471,6.59];
ybs[6918]=['γ Sct',4.8448881,-0.2539636,4.7];
ybs[6919]=['',4.8542251,-0.7312552,6.04];
ybs[6920]=['',4.8474265,-0.254235,5.96];
ybs[6921]=['',4.8494075,-0.3266139,5.66];
ybs[6922]=['δ1 Tel',4.8575734,-0.8010862,4.96];
ybs[6923]=['60 Ser',4.8465596,-0.0343877,5.39];
ybs[6924]=['',4.8538658,-0.5754945,5.34];
ybs[6925]=['',4.8581945,-0.759071,5.72];
ybs[6926]=['δ2 Tel',4.8587691,-0.7983302,5.07];
ybs[6927]=['',4.8663478,-1.0243705,6.44];
ybs[6928]=['',4.8491106,-0.0996388,6.28];
ybs[6929]=['',4.8481035,0.0712176,6.69];
ybs[6930]=['',4.8597923,-0.6926772,5.16];
ybs[6931]=['',4.8452319,0.4168024,5.9];
ybs[6932]=['',4.854805,-0.3209115,5.14];
ybs[6933]=['42 Dra',4.8259974,1.1445265,4.82];
ybs[6934]=['',4.8544813,-0.1881456,5.72];
ybs[6935]=['',4.856796,-0.3335127,6.68];
ybs[6936]=['',4.8626751,-0.6959585,6.22];
ybs[6937]=['',4.834495,1.0395704,6.43];
ybs[6938]=['',4.8501477,0.3635649,6.5];
ybs[6939]=['θ CrA',4.8649518,-0.738196,4.64];
ybs[6940]=['κ1 CrA',4.8642364,-0.6755011,6.32];
ybs[6941]=['κ2 CrA',4.8642222,-0.6756029,5.65];
ybs[6942]=['',4.8701963,-0.9228326,6.22];
ybs[6943]=['',4.8519556,0.2957332,5.77];
ybs[6944]=['',4.8586272,-0.2553037,6.37];
ybs[6945]=['61 Ser',4.8564159,-0.0172253,5.94];
ybs[6946]=['',4.8569824,0.0641568,6.43];
ybs[6947]=['',4.8602715,-0.2591644,5.5];
ybs[6948]=['',4.8664457,-0.5759496,5.28];
ybs[6949]=['24 Sgr',4.8657316,-0.4191474,5.49];
ybs[6950]=['',4.8643213,-0.258948,5.76];
ybs[6951]=['',4.8628321,-0.1028794,6.36];
ybs[6952]=['',4.9602025,-1.4536771,7.16];
ybs[6953]=['25 Sgr',4.868597,-0.4224579,6.51];
ybs[6954]=['',4.8590974,0.4124809,5.84];
ybs[6955]=['',4.8623709,0.144603,6.42];
ybs[6956]=['',4.8590647,0.5335586,5.48];
ybs[6957]=['',4.8719874,-0.3634299,6.48];
ybs[6958]=['',4.8702456,-0.1912804,5.14];
ybs[6959]=['',4.8614637,0.5394683,6.59];
ybs[6960]=['',4.8751525,-0.518031,6.37];
ybs[6961]=['α Sct',4.8708762,-0.1435783,3.85];
ybs[6962]=['',4.8548935,0.9098689,6.56];
ybs[6963]=['',4.8660165,0.3575069,6.57];
ybs[6964]=['',4.8684017,0.1904008,6.4];
ybs[6965]=['',4.8699586,0.3180168,5.78];
ybs[6966]=['45 Dra',4.8560664,0.9959165,4.77];
ybs[6967]=['',4.8489937,1.1423466,6.59];
ybs[6968]=['',4.8710403,0.4123057,5.61];
ybs[6969]=['',4.8729578,0.2965939,6.21];
ybs[6970]=['ζ Pav',4.9105638,-1.2462727,4.01];
ybs[6971]=['',4.8625321,0.9140389,5.36];
ybs[6972]=['',4.8693206,0.6017098,6.1];
ybs[6973]=['',4.8757587,0.1595375,5.39];
ybs[6974]=['',4.8903539,-0.8358364,5.86];
ybs[6975]=['',4.8766662,0.1167689,5.45];
ybs[6976]=['',4.8831416,-0.3731285,5.94];
ybs[6977]=['',4.8836037,-0.2440941,6.47];
ybs[6978]=['',4.8858656,-0.4099011,5.81];
ybs[6979]=['',4.8915285,-0.7533912,5.37];
ybs[6980]=['',4.8789353,0.1996718,6.42];
ybs[6981]=['',4.8810446,-0.005071,5.75];
ybs[6982]=['',4.9346841,-1.3586135,6.39];
ybs[6983]=['',4.8784996,0.2830397,6.29];
ybs[6984]=['',4.9059953,-1.1278595,6.37];
ybs[6985]=['',4.8754581,0.5844625,5.42];
ybs[6986]=['',4.887418,-0.3670841,5.86];
ybs[6987]=['',4.8846125,-0.0554022,6.49];
ybs[6988]=['',4.8842066,-0.0190953,6.66];
ybs[6989]=['α Lyr',4.8765673,0.6772243,0.03];
ybs[6990]=['',4.8840014,0.1545164,6.4];
ybs[6991]=['',4.875521,0.754686,6.2];
ybs[6992]=['',4.9114062,-1.1262441,5.78];
ybs[6993]=['',4.9002761,-0.8390507,6.49];
ybs[6994]=['',4.837823,1.353703,5.64];
ybs[6995]=['',4.8918107,-0.1356201,5.84];
ybs[6996]=['',4.8896454,0.0922239,6.38];
ybs[6997]=['',4.8816232,0.6926708,6.04];
ybs[6998]=['',4.8906418,0.128776,6.28];
ybs[6999]=['26 Sgr',4.9004866,-0.415603,6.23];
ybs[7000]=['',4.9194368,-1.1318183,4.79];
ybs[7001]=['',4.8707044,1.1433057,6.06];
ybs[7002]=['',4.8994687,-0.2538272,6.42];
ybs[7003]=['',4.9176924,-1.0659107,6.04];
ybs[7004]=['',4.89048,0.5387739,6.36];
ybs[7005]=['',4.887848,0.7147954,6.25];
ybs[7006]=['',4.8770718,1.091622,5.74];
ybs[7007]=['',4.8908407,0.669985,6.45];
ybs[7008]=['δ Sct',4.9017371,-0.1576258,4.72];
ybs[7009]=['λ CrA',4.9095512,-0.668489,5.13];
ybs[7010]=['',4.9180366,-0.992373,6.22];
ybs[7011]=['',4.9049312,-0.3361914,6.35];
ybs[7012]=['',4.9031007,-0.1230801,6.15];
ybs[7013]=['',4.8057865,1.4518802,6.17];
ybs[7014]=['',4.9109765,-0.6404736,6.32];
ybs[7015]=['',4.9402874,-1.2735588,6.06];
ybs[7016]=['',4.8884534,0.9113409,6];
ybs[7017]=['',4.9117707,-0.6216806,4.87];
ybs[7018]=['',4.8976669,0.5521979,6.41];
ybs[7019]=['',4.9147389,-0.6922639,5.43];
ybs[7020]=['ε Sct',4.9071565,-0.14405,4.9];
ybs[7021]=['',4.8994635,0.606811,6.47];
ybs[7022]=['',4.9085676,-0.1186236,6.31];
ybs[7023]=['',4.9134773,-0.4361337,5.83];
ybs[7024]=['θ Pav',4.9333463,-1.1353935,5.73];
ybs[7025]=['',4.9243134,-0.8739007,6.54];
ybs[7026]=['',4.9154318,-0.3661469,6.36];
ybs[7027]=['φ Sgr',4.9171817,-0.4706793,3.17];
ybs[7028]=['4 Aql',4.9125178,0.0363454,5.02];
ybs[7029]=['',4.9042092,0.686296,6.45];
ybs[7030]=['',4.8917933,1.095544,6.09];
ybs[7031]=['',4.905771,0.6384139,6.01];
ybs[7032]=['',4.9071301,0.5576077,5.7];
ybs[7033]=['',4.918468,-0.3417936,6.42];
ybs[7034]=['28 Sgr',4.9199874,-0.3904127,5.37];
ybs[7035]=['',4.9110295,0.4121077,6.31];
ybs[7036]=['',4.9151999,0.0963852,5.83];
ybs[7037]=['46 Dra',4.900136,0.9697163,5.04];
ybs[7038]=['μ CrA',4.9269502,-0.7048017,5.24];
ybs[7039]=['ε1 Lyr',4.9088218,0.692758,5.06];
ybs[7040]=['ε1 Lyr',4.9088146,0.6927774,6.02];
ybs[7041]=['ε2 Lyr',4.909007,0.6917644,5.14];
ybs[7042]=['ε2 Lyr',4.909007,0.6917596,5.37];
ybs[7043]=['',4.9211829,-0.1763065,5.71];
ybs[7044]=['ζ1 Lyr',4.9108287,0.6567208,4.36];
ybs[7045]=['ζ2 Lyr',4.9109602,0.6565367,5.73];
ybs[7046]=['',4.9151307,0.384108,6.51];
ybs[7047]=['5 Aql',4.9197912,-0.0163785,5.9];
ybs[7048]=['',4.9040248,0.9406203,6.11];
ybs[7049]=['110 Her',4.9154781,0.3590001,4.19];
ybs[7050]=['η1 CrA',4.9319425,-0.7619322,5.49];
ybs[7051]=['β Sct',4.9229746,-0.0824526,4.22];
ybs[7052]=['',4.9170334,0.4657449,4.83];
ybs[7053]=['',4.9347681,-0.7991073,5.81];
ybs[7054]=['',4.9243528,-0.0991566,5.2];
ybs[7055]=['',4.9200363,0.3268852,6.17];
ybs[7056]=['η2 CrA',4.9351614,-0.7576306,5.61];
ybs[7057]=['111 Her',4.9215035,0.3177348,4.36];
ybs[7058]=['',4.9333755,-0.6060469,6.62];
ybs[7059]=['',4.9102082,0.9585173,6.23];
ybs[7060]=['',4.9303617,-0.3242298,6.47];
ybs[7061]=['',4.9169084,0.7236953,6.07];
ybs[7062]=['λ Pav',4.9484845,-1.0849187,4.22];
ybs[7063]=['',4.9067483,1.0658733,5.99];
ybs[7064]=['',4.9264644,0.0744451,6.21];
ybs[7065]=['',4.9340323,-0.333662,6.75];
ybs[7066]=['29 Sgr',4.9344121,-0.3542998,5.24];
ybs[7067]=['',4.9267578,0.4108197,6.15];
ybs[7068]=['',4.9199347,0.8087572,6.52];
ybs[7069]=['',4.9250249,0.5546803,6.06];
ybs[7070]=['',4.8997175,1.2359391,6.44];
ybs[7071]=['',4.9339588,-0.1027643,5.99];
ybs[7072]=['',4.9182326,0.9252212,5.88];
ybs[7073]=['',4.9334455,0.0150203,6.25];
ybs[7074]=['',4.9296132,0.3377734,5.88];
ybs[7075]=['κ Tel',4.949235,-0.9089869,5.17];
ybs[7076]=['30 Sgr',4.9396025,-0.3863602,6.61];
ybs[7077]=['',4.9368635,-0.1375732,6.8];
ybs[7078]=['',4.92271,0.8569339,6.4];
ybs[7079]=['',4.9309139,0.4375703,6.59];
ybs[7080]=['',4.9478856,-0.8127827,5.54];
ybs[7081]=['',4.9515812,-0.9059038,6.31];
ybs[7082]=['',4.9397279,-0.1701472,5.83];
ybs[7083]=['',4.9505893,-0.8435777,6.19];
ybs[7084]=['',4.9253608,0.8515721,6.12];
ybs[7085]=['',4.9502472,-0.8126133,6.19];
ybs[7086]=['',4.9327778,0.5524653,6.64];
ybs[7087]=['',4.9380786,0.1920158,6.55];
ybs[7088]=['ν1 Lyr',4.9328599,0.5731234,5.91];
ybs[7089]=['8 Aql',4.9412221,-0.0574588,6.1];
ybs[7090]=['ν2 Lyr',4.9333816,0.5685526,5.25];
ybs[7091]=['',4.9469112,-0.4646776,6.29];
ybs[7092]=['',4.9476477,-0.5123092,6.13];
ybs[7093]=['',4.9480505,-0.5359529,6.63];
ybs[7094]=['β Lyr',4.9342081,0.5827253,3.45];
ybs[7095]=['κ Pav',4.9700567,-1.1729491,4.44];
ybs[7096]=['',4.9572359,-0.8700698,6.6];
ybs[7097]=['',4.9435215,0.2441972,6.14];
ybs[7098]=['',4.9486942,-0.1666731,6.34];
ybs[7099]=['',4.9690483,-1.0955886,6.48];
ybs[7100]=['',4.9410468,0.5028169,6.18];
ybs[7101]=['112 Her',4.9442973,0.3743957,5.48];
ybs[7102]=['33 Sgr',4.9533566,-0.3723273,5.69];
ybs[7103]=['',4.9407183,0.6381669,6.09];
ybs[7104]=['ν1 Sgr',4.9541454,-0.3965034,4.83];
ybs[7105]=['',4.9098652,1.2934299,5.27];
ybs[7106]=['',4.9426724,0.7227273,6.28];
ybs[7107]=['',4.9562668,-0.2718487,5.1];
ybs[7108]=['ν2 Sgr',4.9582796,-0.3952107,4.99];
ybs[7109]=['σ Sgr',4.9590671,-0.4584824,2.02];
ybs[7110]=['',4.9643298,-0.7449494,5.36];
ybs[7111]=['',4.9394732,0.9250347,5.51];
ybs[7112]=['50 Dra',4.911783,1.3169671,5.35];
ybs[7113]=['ο Dra',4.9370987,1.0369646,4.66];
ybs[7114]=['',4.959778,-0.285344,5.79];
ybs[7115]=['ω Pav',4.9760528,-1.0501864,5.14];
ybs[7116]=['',4.9621893,-0.4039686,5.93];
ybs[7117]=['',4.9657364,-0.6512707,5.38];
ybs[7118]=['',4.983476,-1.1627954,6.01];
ybs[7119]=['δ1 Lyr',4.9499416,0.645743,5.58];
ybs[7120]=['',4.952527,0.487582,5.62];
ybs[7121]=['113 Her',4.9550505,0.3957048,4.59];
ybs[7122]=['λ Tel',4.974609,-0.9234436,4.87];
ybs[7123]=['',4.9587314,0.1159399,5.57];
ybs[7124]=['',4.9698131,-0.6945471,6.31];
ybs[7125]=['',4.946849,0.885488,4.92];
ybs[7126]=['',4.9519359,0.7199915,7.3];
ybs[7127]=['δ2 Lyr',4.9533354,0.6444794,4.3];
ybs[7128]=['',4.9551069,0.5933397,6.02];
ybs[7129]=['θ1 Ser',4.962137,0.0738549,4.62];
ybs[7130]=['θ2 Ser',4.9622389,0.073826,4.98];
ybs[7131]=['',4.9630305,-0.0309263,6.22];
ybs[7132]=['',4.9630974,0.0436188,6.15];
ybs[7133]=['ξ1 Sgr',4.9679002,-0.3600237,5.08];
ybs[7134]=['',4.9546705,0.7265806,5.44];
ybs[7135]=['',4.9609752,0.3145583,6.63];
ybs[7136]=['',4.9611384,0.3164833,5.69];
ybs[7137]=['η Sct',4.966149,-0.1015384,4.83];
ybs[7138]=['ξ2 Sgr',4.9696116,-0.3678792,3.51];
ybs[7139]=['',4.9727535,-0.541175,6.12];
ybs[7140]=['ε CrA',4.9746581,-0.6471375,4.87];
ybs[7141]=['',4.948556,1.0038011,6.22];
ybs[7142]=['',4.9538023,0.8532369,5.77];
ybs[7143]=['',4.9724296,-0.4336728,6.62];
ybs[7144]=['',4.9768041,-0.6894965,6.49];
ybs[7145]=['13 Lyr',4.9565534,0.7674832,4.04];
ybs[7146]=['64 Ser',4.9668041,0.0447459,5.57];
ybs[7147]=['',4.9726336,-0.3927056,6.14];
ybs[7148]=['',4.9050651,1.3956475,6.39];
ybs[7149]=['',4.9989402,-1.1994521,5.88];
ybs[7150]=['',4.9645501,0.5747314,5.22];
ybs[7151]=['',4.9715658,0.1094197,6.21];
ybs[7152]=['',4.9770041,-0.3235384,6.37];
ybs[7153]=['',4.9705122,0.3035083,5.38];
ybs[7154]=['',4.9765783,-0.2235947,5.53];
ybs[7155]=['10 Aql',4.9729795,0.2432264,5.89];
ybs[7156]=['',4.9814695,-0.4347996,6.36];
ybs[7157]=['',4.9848181,-0.6463034,6.69];
ybs[7158]=['',4.9848981,-0.6463227,6.4];
ybs[7159]=['',4.9726323,0.3459822,6.5];
ybs[7160]=['11 Aql',4.974357,0.2382694,5.23];
ybs[7161]=['',4.9753374,0.1775046,6.75];
ybs[7162]=['',4.968663,0.6683716,5.89];
ybs[7163]=['48 Dra',4.9615225,1.0095518,5.66];
ybs[7164]=['ε Aql',4.9766108,0.2635083,4.02];
ybs[7165]=['',4.9898522,-0.7309279,6.23];
ybs[7166]=['γ Lyr',4.9729259,0.5710483,3.24];
ybs[7167]=['',4.9717787,0.7104935,6.22];
ybs[7168]=['υ Dra',4.9486067,1.2448386,4.82];
ybs[7169]=['',4.9767699,0.4583266,5.27];
ybs[7170]=['',4.986687,-0.3955778,6.24];
ybs[7171]=['',4.9778251,0.3987109,6.29];
ybs[7172]=['',4.9646318,1.0167185,6.46];
ybs[7173]=['',4.9737321,0.6849911,6.41];
ybs[7174]=['',4.9860962,-0.2661964,6.32];
ybs[7175]=['',4.9589867,1.1394538,5.63];
ybs[7176]=['ζ CrA',4.9940938,-0.7341535,4.75];
ybs[7177]=['',4.9983996,-0.8898874,5.93];
ybs[7178]=['',4.9632552,1.0895209,6.45];
ybs[7179]=['λ Lyr',4.977629,0.5615647,4.93];
ybs[7180]=['12 Aql',4.9862956,-0.0996283,4.02];
ybs[7181]=['ζ Sgr',4.9912688,-0.5209663,2.6];
ybs[7182]=['',4.9904009,-0.4331196,5.65];
ybs[7183]=['',4.9720172,0.8873013,6.3];
ybs[7184]=['',4.9946588,-0.6670976,5.74];
ybs[7185]=['',4.9828626,0.3375467,6.39];
ybs[7186]=['',4.9429827,1.3231997,6.22];
ybs[7187]=['',4.9840485,0.3641459,6.69];
ybs[7188]=['',4.9785009,0.7105937,6.65];
ybs[7189]=['',4.9834523,0.4594012,5.69];
ybs[7190]=['',4.9928073,-0.3353472,6.05];
ybs[7191]=['',4.9815001,0.5904865,6.01];
ybs[7192]=['',4.9930347,-0.3328694,6.37];
ybs[7193]=['',4.9847774,0.4373156,6.72];
ybs[7194]=['',4.9859497,0.3891127,6.4];
ybs[7195]=['',4.9887874,0.1466963,6.3];
ybs[7196]=['14 Aql',4.9915862,-0.0640134,5.42];
ybs[7197]=['',4.9774324,0.8824975,5.38];
ybs[7198]=['',4.9992018,-0.5413134,5.5];
ybs[7199]=['',4.9853718,0.5873378,6.39];
ybs[7200]=['ρ Tel',5.0088696,-0.9129447,5.16];
ybs[7201]=['',4.9941431,0.0322951,5.83];
ybs[7202]=['16 Lyr',4.9830006,0.8196955,5.01];
ybs[7203]=['',4.9906394,0.3436945,6.09];
ybs[7204]=['ο Sgr',4.9999695,-0.3789038,3.77];
ybs[7205]=['49 Dra',4.9791143,0.9719445,5.48];
ybs[7206]=['',4.9968924,0.058684,6.73];
ybs[7207]=['',4.9981758,-0.098665,6.9];
ybs[7208]=['',5.0266688,-1.1936298,5.33];
ybs[7209]=['',4.9942077,0.371743,6.52];
ybs[7210]=['',5.0111337,-0.8423996,5.97];
ybs[7211]=['',4.9687056,1.2140524,6.52];
ybs[7212]=['15 Aql',5.0005491,-0.0697994,5.42];
ybs[7213]=['γ CrA',5.0082176,-0.6463022,4.93];
ybs[7214]=['γ CrA',5.0082176,-0.6463022,4.99];
ybs[7215]=['σ Oct',5.6091363,-1.5508139,5.47];
ybs[7216]=['',4.9855282,0.9126638,6.31];
ybs[7217]=['',5.0041202,-0.2727554,5.97];
ybs[7218]=['',5.0019919,-0.0258386,6.53];
ybs[7219]=['',5.0102487,-0.659335,6.16];
ybs[7220]=['',5.0202546,-0.9719057,6.49];
ybs[7221]=['τ Sgr',5.0100529,-0.4823634,3.32];
ybs[7222]=['ζ Aql',5.0019063,0.2425255,2.99];
ybs[7223]=['λ Aql',5.0061988,-0.0846434,3.44];
ybs[7224]=['',4.9992534,0.5546055,5.56];
ybs[7225]=['',4.9993287,0.5369584,6.06];
ybs[7226]=['',5.0093104,-0.2826696,6.03];
ybs[7227]=['',5.0126024,-0.4992253,6.04];
ybs[7228]=['',5.01058,-0.3264557,6.29];
ybs[7229]=['δ CrA',5.0168227,-0.706209,4.59];
ybs[7230]=['',5.0062906,0.1442135,6.09];
ybs[7231]=['',5.002915,0.5227989,6.31];
ybs[7232]=['',5.0099559,0.0117789,6.56];
ybs[7233]=['',5.0156139,-0.4297553,6.3];
ybs[7234]=['',4.9614413,1.3452854,6.54];
ybs[7235]=['18 Aql',5.0088407,0.1938101,5.09];
ybs[7236]=['',5.0155591,-0.3360843,5.54];
ybs[7237]=['',5.0068854,0.4238315,5.77];
ybs[7238]=['51 Dra',4.9976647,0.9325063,5.38];
ybs[7239]=['',4.9990364,0.8718875,6.43];
ybs[7240]=['',5.0066525,0.5002379,5.55];
ybs[7241]=['α CrA',5.0215786,-0.660957,4.11];
ybs[7242]=['',5.0225123,-0.6945287,6.46];
ybs[7243]=['',5.022077,-0.6305921,6.56];
ybs[7244]=['',5.0239475,-0.7305575,5.88];
ybs[7245]=['',5.0045048,0.7233796,6.49];
ybs[7246]=['β CrA',5.0240843,-0.6860221,4.11];
ybs[7247]=['',5.0129057,0.294732,6.07];
ybs[7248]=['17 Lyr',5.009963,0.5678419,5.23];
ybs[7249]=['ι Lyr',5.0092418,0.6306483,5.28];
ybs[7250]=['',5.0131831,0.3793036,6.23];
ybs[7251]=['π Sgr',5.0221024,-0.3663291,2.89];
ybs[7252]=['',5.0222237,-0.3450359,6.13];
ybs[7253]=['19 Aql',5.0178326,0.1065945,5.22];
ybs[7254]=['',5.0160261,0.294704,6.48];
ybs[7255]=['',5.0284342,-0.6801524,6.36];
ybs[7256]=['',5.0218088,-0.0068688,6.34];
ybs[7257]=['',5.029215,-0.5142957,6.3];
ybs[7258]=['',5.0367756,-0.8805253,6.13];
ybs[7259]=['',5.0170503,0.604488,6.74];
ybs[7260]=['',5.0333006,-0.6553206,6.57];
ybs[7261]=['τ Pav',5.0555975,-1.2069416,6.27];
ybs[7262]=['',5.0130805,0.9155866,5.81];
ybs[7263]=['',5.0339194,-0.37738,6.41];
ybs[7264]=['',5.0374093,-0.4515257,5.8];
ybs[7265]=['',5.0581683,-1.1627939,5.53];
ybs[7266]=['20 Aql',5.0343556,-0.1379437,5.34];
ybs[7267]=['',5.0280613,0.4672433,6.36];
ybs[7268]=['',5.0446431,-0.7881289,5.92];
ybs[7269]=['',5.037034,-0.2137394,5.51];
ybs[7270]=['19 Lyr',5.0289629,0.5466141,5.98];
ybs[7271]=['',5.0268336,0.7062354,6.18];
ybs[7272]=['',5.0330609,0.2946491,6.73];
ybs[7273]=['',5.0330537,0.3768202,5.93];
ybs[7274]=['21 Aql',5.038525,0.040665,5.15];
ybs[7275]=['',5.0385201,0.0969035,6.49];
ybs[7276]=['',5.0520739,-0.7928806,5.4];
ybs[7277]=['55 Dra',5.0171061,1.1521413,6.25];
ybs[7278]=['',5.0474817,-0.4213561,6.25];
ybs[7279]=['ψ Sgr',5.0474664,-0.4401621,4.85];
ybs[7280]=['',5.0292209,0.8707378,6.75];
ybs[7281]=['',5.0292498,0.8707717,6.57];
ybs[7282]=['53 Dra',5.0268106,0.9929941,5.12];
ybs[7283]=['',5.0522157,-0.584405,6.25];
ybs[7284]=['',5.0605177,-0.9311002,6.38];
ybs[7285]=['η Lyr',5.0372664,0.6838614,4.39];
ybs[7286]=['',5.0437174,0.353259,6];
ybs[7287]=['',5.0451759,0.2639056,5.57];
ybs[7288]=['1 Sge',5.0447553,0.3712185,5.64];
ybs[7289]=['',5.0449144,0.533433,5.85];
ybs[7290]=['22 Aql',5.0506807,0.0850387,5.59];
ybs[7291]=['43 Sgr',5.0563639,-0.3301264,4.96];
ybs[7292]=['',5.047389,0.4798413,6.54];
ybs[7293]=['1 Vul',5.0487936,0.3739847,4.77];
ybs[7294]=['',5.0500463,0.2545093,5.63];
ybs[7295]=['',5.0475438,0.4880689,6.16];
ybs[7296]=['54 Dra',5.0365199,1.0077751,4.99];
ybs[7297]=['δ Dra',5.0289436,1.181539,3.07];
ybs[7298]=['',5.0433732,0.874546,6.27];
ybs[7299]=['59 Dra',5.0107818,1.3368226,5.13];
ybs[7300]=['',5.0563873,0.0361268,6.19];
ybs[7301]=['θ Lyr',5.0487152,0.6662116,4.36];
ybs[7302]=['ω1 Aql',5.0561269,0.2030432,5.28];
ybs[7303]=['',5.0659144,-0.6175358,5.59];
ybs[7304]=['',5.0622033,-0.2704832,6.06];
ybs[7305]=['2 Vul',5.0553214,0.402538,5.43];
ybs[7306]=['23 Aql',5.05964,0.0196153,5.1];
ybs[7307]=['',5.0885612,-1.1925778,6.34];
ybs[7308]=['24 Aql',5.0610024,0.0065958,6.41];
ybs[7309]=['',5.0503065,0.8209482,6];
ybs[7310]=['',5.069103,-0.5546348,6.58];
ybs[7311]=['',5.0562415,0.5421082,6.68];
ybs[7312]=['',5.0608269,0.1685428,6.32];
ybs[7313]=['',5.0601712,0.342944,6.58];
ybs[7314]=['',5.0695792,-0.3903059,5.58];
ybs[7315]=['κ Cyg',5.0508704,0.9321175,3.77];
ybs[7316]=['η Tel',5.0810813,-0.94916,5.05];
ybs[7317]=['',5.0738813,-0.609885,6.48];
ybs[7318]=['28 Aql',5.0641227,0.216662,5.53];
ybs[7319]=['ω2 Aql',5.0651478,0.202008,6.02];
ybs[7320]=['26 Aql',5.0686083,-0.0938339,5.01];
ybs[7321]=['',5.0771413,-0.7326148,6.34];
ybs[7322]=['',5.0607166,0.5834228,6.6];
ybs[7323]=['27 Aql',5.0686639,-0.0148818,5.49];
ybs[7324]=['β1 Sgr',5.0793766,-0.7752452,4.01];
ybs[7325]=['',5.0603274,0.6542196,6.22];
ybs[7326]=['',5.0737506,-0.3350001,6.26];
ybs[7327]=['ρ1 Sgr',5.0739403,-0.310793,3.93];
ybs[7328]=['',5.0578519,0.8658272,6.31];
ybs[7329]=['υ Sgr',5.0741057,-0.2777672,4.61];
ybs[7330]=['β2 Sgr',5.0819346,-0.7811892,4.29];
ybs[7331]=['ρ2 Sgr',5.0747203,-0.3188395,5.87];
ybs[7332]=['',5.0630763,0.6522225,6.31];
ybs[7333]=['',5.0671276,0.6148023,6.31];
ybs[7334]=['',5.0765391,-0.1424315,6.31];
ybs[7335]=['α Sgr',5.0845901,-0.7081662,3.97];
ybs[7336]=['',5.076337,-0.0037022,5.83];
ybs[7337]=['',5.0868228,-0.7623887,6.17];
ybs[7338]=['',5.0617156,0.9497221,6.26];
ybs[7339]=['τ Dra',5.040309,1.2809392,4.45];
ybs[7340]=['',5.0797047,-0.1284484,6.32];
ybs[7341]=['',5.0779596,0.1737236,6.35];
ybs[7342]=['',5.086662,-0.4856226,6.04];
ybs[7343]=['',5.0642554,1.0067849,5.91];
ybs[7344]=['',5.079229,0.2611331,6.64];
ybs[7345]=['3 Vul',5.0775479,0.459075,5.18];
ybs[7346]=['',5.07596,0.5857056,6.06];
ybs[7347]=['',5.0891941,-0.5108138,5.93];
ybs[7348]=['',5.0611189,1.124512,6.52];
ybs[7349]=['χ1 Sgr',5.0899001,-0.4270268,5.03];
ybs[7350]=['χ3 Sgr',5.0908385,-0.4174887,5.43];
ybs[7351]=['',5.0818529,0.354397,6.4];
ybs[7352]=['',5.0692865,1.0089178,6.43];
ybs[7353]=['',5.08813,-0.0845182,6.52];
ybs[7354]=['',5.0898845,-0.2418179,5.69];
ybs[7355]=['',5.0803828,0.5805507,6.37];
ybs[7356]=['2 Sge',5.0845335,0.2963406,6.25];
ybs[7357]=['',5.1026482,-0.9474038,5.69];
ybs[7358]=['π Dra',5.064808,1.1476254,4.59];
ybs[7359]=['2 Cyg',5.082988,0.517709,4.97];
ybs[7360]=['31 Aql',5.0873318,0.2091955,5.16];
ybs[7361]=['',5.0841296,0.4909444,6.53];
ybs[7362]=['50 Sgr',5.0943539,-0.3793369,5.59];
ybs[7363]=['',5.0825626,0.6369241,6.36];
ybs[7364]=['δ Aql',5.0899287,0.0550923,3.36];
ybs[7365]=['',5.0935235,-0.2619892,5.72];
ybs[7366]=['',5.0944877,-0.2532268,6.7];
ybs[7367]=['',5.0973832,-0.5183762,5.67];
ybs[7368]=['',5.0786195,0.8781124,6.51];
ybs[7369]=['',5.0814797,0.7579804,5.84];
ybs[7370]=['',5.1194603,-1.1936165,5.96];
ybs[7371]=['',5.0887993,0.3545311,6.31];
ybs[7372]=['4 Vul',5.0892677,0.3462804,5.16];
ybs[7373]=['',5.0888703,0.4355389,6.19];
ybs[7374]=['ν Aql',5.0944692,0.0066484,4.66];
ybs[7375]=['',5.1118005,-0.9668668,6.13];
ybs[7376]=['',5.0935431,0.2280468,5.74];
ybs[7377]=['5 Vul',5.092501,0.3515078,5.63];
ybs[7378]=['',5.0936361,0.3479077,5.81];
ybs[7379]=['',5.1087886,-0.7575095,5.71];
ybs[7380]=['μ Tel',5.1148115,-0.9610775,6.3];
ybs[7381]=['λ UMi',4.4155421,1.5533045,6.38];
ybs[7382]=['4 Cyg',5.0915188,0.6345991,5.15];
ybs[7383]=['',5.098577,0.2500232,6.32];
ybs[7384]=['',5.1023633,0.0518961,5.85];
ybs[7385]=['',5.1100406,-0.4702205,5.52];
ybs[7386]=['',5.1118927,-0.5593454,6.6];
ybs[7387]=['35 Aql',5.1053185,0.0347973,5.8];
ybs[7388]=['',5.0882804,1.0134958,6.6];
ybs[7389]=['',5.1070995,-0.1221775,6.61];
ybs[7390]=['',5.0977958,0.6629432,6.34];
ybs[7391]=['',5.1066095,0.0050563,6.25];
ybs[7392]=['α Vul',5.1031748,0.4312408,4.44];
ybs[7393]=['8 Vul',5.1042401,0.4330511,5.81];
ybs[7394]=['',5.1064362,0.2555063,5.56];
ybs[7395]=['ι1 Cyg',5.0960715,0.9139096,5.75];
ybs[7396]=['7 Vul',5.106148,0.3547085,6.33];
ybs[7397]=['',5.1143147,-0.371194,6.13];
ybs[7398]=['',5.1247516,-0.927476,5.75];
ybs[7399]=['',5.1103421,0.0514551,6.09];
ybs[7400]=['',5.0905777,1.0925641,6.38];
ybs[7401]=['36 Aql',5.1126556,-0.0479034,5.03];
ybs[7402]=['',5.1119684,0.0608877,6.05];
ybs[7403]=['',5.1261756,-0.7893494,5.61];
ybs[7404]=['β1 Cyg',5.1118456,0.4887606,3.08];
ybs[7405]=['β2 Cyg',5.1119909,0.4888579,5.11];
ybs[7406]=['',5.1117366,0.63308,6.25];
ybs[7407]=['ι2 Cyg',5.1060506,0.903616,3.79];
ybs[7408]=['',5.114701,0.4653348,5.87];
ybs[7409]=['',5.1292749,-0.6979367,5.89];
ybs[7410]=['',5.0631214,1.3900203,6.05];
ybs[7411]=['ι Tel',5.134461,-0.8386788,4.9];
ybs[7412]=['',5.0282535,1.4573304,6.53];
ybs[7413]=['8 Cyg',5.1161442,0.6020988,4.74];
ybs[7414]=['',5.1132268,0.8787921,5.53];
ybs[7415]=['',5.1123287,0.9734798,6.37];
ybs[7416]=['μ Aql',5.1272815,0.1295849,4.45];
ybs[7417]=['37 Aql',5.1323586,-0.1835039,5.12];
ybs[7418]=['51 Sgr',5.1368118,-0.4306155,5.65];
ybs[7419]=['',5.1338938,-0.1293958,6.34];
ybs[7420]=['',5.1343269,-0.2130401,6.27];
ybs[7421]=['',5.1493512,-1.0111639,6.18];
ybs[7422]=['',5.1569079,-1.1630338,6.39];
ybs[7423]=['',5.1239415,0.6773224,6.61];
ybs[7424]=['9 Vul',5.1290076,0.3459123,5];
ybs[7425]=['',5.1332284,0.051657,6.38];
ybs[7426]=['',5.1383697,-0.3282246,6.11];
ybs[7427]=['52 Sgr',5.1397768,-0.4334802,4.6];
ybs[7428]=['9 Cyg',5.1297969,0.5150316,5.38];
ybs[7429]=['',5.123665,0.8605869,5.96];
ybs[7430]=['',5.1410524,-0.3173696,5.64];
ybs[7431]=['',5.1284551,0.741045,5.35];
ybs[7432]=['',5.1360752,0.1954193,6.68];
ybs[7433]=['κ Aql',5.139966,-0.1218314,4.95];
ybs[7434]=['ι Aql',5.1390414,-0.0216316,4.36];
ybs[7435]=['',5.1203057,1.0507547,6.29];
ybs[7436]=['',5.1365344,0.251998,6.38];
ybs[7437]=['',5.1087111,1.2397691,6.07];
ybs[7438]=['',5.1262915,0.8950475,5.73];
ybs[7439]=['',5.1356991,0.395012,6.32];
ybs[7440]=['',5.127989,0.8414349,6.67];
ybs[7441]=['',5.1431786,-0.248784,5.47];
ybs[7442]=['',5.164397,-1.1485097,6.09];
ybs[7443]=['',5.1393078,0.1975778,5.98];
ybs[7444]=['11 Cyg',5.1336325,0.6456138,6.05];
ybs[7445]=['',5.1379203,0.3556926,7.14];
ybs[7446]=['',5.1571518,-0.9489139,6.26];
ybs[7447]=['42 Aql',5.143802,-0.0802855,5.46];
ybs[7448]=['',5.1538165,-0.7894059,6.25];
ybs[7449]=['σ Dra',5.1150204,1.2165966,4.68];
ybs[7450]=['ε Sge',5.1409444,0.2881538,5.66];
ybs[7451]=['',5.1544745,-0.6873948,6.61];
ybs[7452]=['',5.1333662,0.877641,6.52];
ybs[7453]=['',5.1398789,0.5127907,6.43];
ybs[7454]=['',5.1385295,0.6707457,6.5];
ybs[7455]=['',5.1368101,0.7808926,5.17];
ybs[7456]=['θ Cyg',5.135595,0.8773396,4.48];
ybs[7457]=['53 Sgr',5.1533155,-0.4080466,6.34];
ybs[7458]=['',5.1480349,0.0598577,6.35];
ybs[7459]=['',5.145164,0.3635598,6.48];
ybs[7460]=['',5.1546023,-0.4080588,5.97];
ybs[7461]=['σ Aql',5.1496145,0.0950485,5.17];
ybs[7462]=['',5.1502514,0.2900662,6.38];
ybs[7463]=['54 Sgr',5.1569873,-0.2835202,6.2];
ybs[7464]=['',5.1422322,0.8610037,6.47];
ybs[7465]=['φ Cyg',5.1495357,0.527115,4.69];
ybs[7466]=['α Sge',5.1531403,0.3152479,4.37];
ybs[7467]=['45 Aql',5.1564734,-0.0099887,5.67];
ybs[7468]=['',5.1509974,0.5938912,6.1];
ybs[7469]=['',5.15469,0.3582344,6.5];
ybs[7470]=['14 Cyg',5.1491946,0.748161,5.4];
ybs[7471]=['',5.1449916,0.9603087,5.82];
ybs[7472]=['',5.1554006,0.414799,6.64];
ybs[7473]=['',5.1576257,0.241981,6.01];
ybs[7474]=['',5.1495753,0.8029603,6.2];
ybs[7475]=['β Sge',5.1573125,0.3058694,4.37];
ybs[7476]=['55 Sgr',5.1648117,-0.2805487,5.06];
ybs[7477]=['',5.1580083,0.39273,6.36];
ybs[7478]=['',5.1705036,-0.6543016,6.16];
ybs[7479]=['',5.1546159,0.7526989,6.16];
ybs[7480]=['46 Aql',5.1625745,0.2136768,6.34];
ybs[7481]=['',5.2345906,-1.4188407,6.39];
ybs[7482]=['',5.1551204,0.795412,5.06];
ybs[7483]=['',5.1693264,-0.2691281,5.49];
ybs[7484]=['χ Aql',5.1641285,0.2072801,5.27];
ybs[7485]=['',5.200004,-1.2644976,5.41];
ybs[7486]=['',5.160328,0.7034229,6.23];
ybs[7487]=['',5.1510586,1.056895,6.51];
ybs[7488]=['',5.1645944,0.5128014,6.49];
ybs[7489]=['',5.1641352,0.5668186,5.94];
ybs[7490]=['16 Cyg',5.1590407,0.8826908,5.96];
ybs[7491]=['',5.159267,0.8825554,6.2];
ybs[7492]=['',5.1660224,0.5363126,6.05];
ybs[7493]=['10 Vul',5.1686503,0.4506796,5.49];
ybs[7494]=['',5.1806802,-0.5560163,5.52];
ybs[7495]=['',5.169543,0.4744808,6.28];
ybs[7496]=['',5.1597174,0.9688777,6.48];
ybs[7497]=['ν Tel',5.1909594,-0.9828003,5.35];
ybs[7498]=['ψ Aql',5.1728156,0.2330587,6.26];
ybs[7499]=['',5.1689181,0.5971232,6.05];
ybs[7500]=['',5.2005631,-1.1651764,6.45];
ybs[7501]=['',5.1680679,0.7299513,5.84];
ybs[7502]=['56 Sgr',5.1817019,-0.3440002,4.86];
ybs[7503]=['',5.1790073,-0.0494318,6.48];
ybs[7504]=['15 Cyg',5.1706003,0.6528362,4.89];
ybs[7505]=['',5.1732989,0.5116485,6.82];
ybs[7506]=['υ Aql',5.1777849,0.1337678,5.91];
ybs[7507]=['',5.1723121,0.6015169,6.57];
ybs[7508]=['',5.1945538,-0.9221532,6.25];
ybs[7509]=['',5.1645843,1.0134458,6.22];
ybs[7510]=['',5.1727841,0.7115222,6.34];
ybs[7511]=['',5.2052702,-1.1440882,6.05];
ybs[7512]=['γ Aql',5.1802858,0.1861322,2.72];
ybs[7513]=['',5.1665108,0.9964516,6.27];
ybs[7514]=['',5.2017,-1.0647929,6.21];
ybs[7515]=['δ Cyg',5.1732181,0.7885649,2.87];
ybs[7516]=['',5.1767015,0.6307978,6.43];
ybs[7517]=['',5.1776109,0.611979,6.09];
ybs[7518]=['',5.203146,-1.0321815,5.42];
ybs[7519]=['',5.1888617,-0.2382591,6.11];
ybs[7520]=['',5.1815229,0.4395663,6.62];
ybs[7521]=['17 Cyg',5.1801582,0.5895559,4.99];
ybs[7522]=['',5.1808795,0.574911,6.18];
ybs[7523]=['δ Sge',5.1849473,0.3243853,3.82];
ybs[7524]=['',5.1998685,-0.8291028,5.94];
ybs[7525]=['',5.1943895,-0.5015424,6.05];
ybs[7526]=['',5.1865406,0.4439386,5.95];
ybs[7527]=['',5.1930744,-0.1888151,6.04];
ybs[7528]=['',5.1900878,0.1875604,6.44];
ybs[7529]=['',5.1844617,0.6712403,5.77];
ybs[7530]=['π Aql',5.1909028,0.2071387,5.72];
ybs[7531]=['',5.1673537,1.2110335,5.92];
ybs[7532]=['ζ Sge',5.1918662,0.3350102,5];
ybs[7533]=['',5.1838303,0.8370509,6.12];
ybs[7534]=['',5.2108981,-0.9584807,5.74];
ybs[7535]=['',5.2110003,-0.9585774,6.5];
ybs[7536]=['',5.1901496,0.6172129,6.53];
ybs[7537]=['',5.1907222,0.5845035,6.44];
ybs[7538]=['',5.206448,-0.695001,5.33];
ybs[7539]=['51 Aql',5.2006763,-0.1869304,5.39];
ybs[7540]=['',5.197966,0.1388504,6.51];
ybs[7541]=['',5.1931576,0.6765354,6.11];
ybs[7542]=['',5.1955889,0.4973036,6.38];
ybs[7543]=['α Aql',5.2000812,0.1557112,0.77];
ybs[7544]=['',5.2204952,-1.06667,6.24];
ybs[7545]=['',5.2021844,-0.0420166,6.13];
ybs[7546]=['ο Aql',5.2010941,0.1827171,5.11];
ybs[7547]=['57 Sgr',5.2071433,-0.3314568,5.92];
ybs[7548]=['',5.20229,0.1690136,6.25];
ybs[7549]=['',5.1782266,1.195369,6.34];
ybs[7550]=['χ Cyg',5.1982619,0.5753878,4.23];
ybs[7551]=['12 Vul',5.2008662,0.3955503,4.95];
ybs[7552]=['19 Cyg',5.1979965,0.676762,5.12];
ybs[7553]=['',5.1981391,0.7095261,5.69];
ybs[7554]=['',5.1989874,0.6611236,6.06];
ybs[7555]=['',5.2055578,0.2039017,6.13];
ybs[7556]=['η Aql',5.2076993,0.0184931,3.9];
ybs[7557]=['',5.2109499,-0.2539233,6.48];
ybs[7558]=['',5.2064785,0.1816067,6.54];
ybs[7559]=['',5.204962,0.4371353,5.57];
ybs[7560]=['9 Sge',5.2066568,0.3268283,6.23];
ybs[7561]=['',5.2114833,-0.053408,5.65];
ybs[7562]=['20 Cyg',5.1973393,0.9257426,5.03];
ybs[7563]=['',5.2007986,0.8278208,6.2];
ybs[7564]=['',5.216454,-0.4168938,6.18];
ybs[7565]=['',5.2393779,-1.2061436,5.75];
ybs[7566]=['',5.2115369,0.0777538,6.53];
ybs[7567]=['ι Sgr',5.2214574,-0.7297749,4.13];
ybs[7568]=['ε Dra',5.1840149,1.2273094,3.83];
ybs[7569]=['',5.2055534,0.6368023,6.1];
ybs[7570]=['56 Aql',5.2152512,-0.1486917,5.79];
ybs[7571]=['',5.2202593,-0.5758046,6.46];
ybs[7572]=['',5.2300797,-1.0100123,6.53];
ybs[7573]=['',5.2308044,-1.0270426,5.26];
ybs[7574]=['',5.240133,-1.1991319,6.39];
ybs[7575]=['',5.2037118,0.8217221,5.62];
ybs[7576]=['ε Pav',5.2487299,-1.2715201,3.96];
ybs[7577]=['',5.2042407,0.8375086,5.91];
ybs[7578]=['13 Vul',5.2112598,0.4212201,4.58];
ybs[7579]=['57 Aql',5.2173782,-0.1426326,5.71];
ybs[7580]=['57 Aql',5.2174148,-0.1428071,6.49];
ybs[7581]=['ξ Aql',5.215215,0.1486352,4.71];
ybs[7582]=['58 Aql',5.2176411,0.0057356,5.61];
ybs[7583]=['ω Sgr',5.2232797,-0.4580427,4.7];
ybs[7584]=['',5.2171034,0.1255808,6.15];
ybs[7585]=['',5.2203859,-0.1165686,6.51];
ybs[7586]=['',5.2081802,0.8353436,6.29];
ybs[7587]=['',5.2158636,0.4254121,5.52];
ybs[7588]=['β Aql',5.2199262,0.1127817,3.71];
ybs[7589]=['μ1 Pav',5.2462792,-1.1674819,5.76];
ybs[7590]=['59 Sgr',5.2281393,-0.4732284,4.52];
ybs[7591]=['',5.2318404,-0.6632595,6.55];
ybs[7592]=['',5.2165882,0.6466634,5.76];
ybs[7593]=['',5.2182122,0.527969,6.57];
ybs[7594]=['23 Cyg',5.2085379,1.004923,5.14];
ybs[7595]=['10 Sge',5.2226926,0.2913001,5.36];
ybs[7596]=['φ Aql',5.2238065,0.2003557,5.28];
ybs[7597]=['',5.2096017,1.0430654,6.06];
ybs[7598]=['μ2 Pav',5.252764,-1.1673787,5.31];
ybs[7599]=['22 Cyg',5.2211296,0.6726865,4.94];
ybs[7600]=['61 Sgr',5.2321063,-0.2693908,5.02];
ybs[7601]=['η Cyg',5.2232341,0.6132907,3.89];
ybs[7602]=['',5.2268294,0.3674619,6.48];
ybs[7603]=['',5.2369646,-0.5319968,6.28];
ybs[7604]=['60 Sgr',5.2368538,-0.4562059,4.83];
ybs[7605]=['ψ Cyg',5.2192362,0.9161962,4.92];
ybs[7606]=['',5.2250579,0.6336706,6.02];
ybs[7607]=['',5.2444132,-0.8603346,6.17];
ybs[7608]=['11 Sge',5.2302601,0.2940088,5.53];
ybs[7609]=['θ1 Sgr',5.240641,-0.6146901,4.37];
ybs[7610]=['θ2 Sgr',5.2411311,-0.6045906,5.3];
ybs[7611]=['',5.2510778,-1.0352931,5.13];
ybs[7612]=['',5.2175414,1.0176215,6.09];
ybs[7613]=['',5.2441041,-0.7502433,6.14];
ybs[7614]=['',5.2270163,0.7055333,5.45];
ybs[7615]=['',5.2430573,-0.6570248,5.95];
ybs[7616]=['',5.2458081,-0.786364,5.81];
ybs[7617]=['',5.2431949,-0.5872404,5.66];
ybs[7618]=['',5.2242579,0.8893897,6.43];
ybs[7619]=['',5.2199001,1.028025,4.96];
ybs[7620]=['',5.2218504,0.9903436,6.12];
ybs[7621]=['γ Sge',5.2345344,0.3411935,3.47];
ybs[7622]=['',5.2378173,0.025042,6.17];
ybs[7623]=['',5.2399577,-0.1728071,5.88];
ybs[7624]=['',5.2299934,0.7385737,6.43];
ybs[7625]=['',5.2483355,-0.7113296,6.29];
ybs[7626]=['',5.2335694,0.5417547,5.49];
ybs[7627]=['14 Vul',5.2362349,0.4041884,5.67];
ybs[7628]=['',5.2329893,0.6660554,6.32];
ybs[7629]=['',5.2473918,-0.3958284,6.01];
ybs[7630]=['',5.2688267,-1.1739254,6.07];
ybs[7631]=['13 Sge',5.2402674,0.3067234,5.37];
ybs[7632]=['',5.2359242,0.7998694,5.92];
ybs[7633]=['25 Cyg',5.2389208,0.6475166,5.19];
ybs[7634]=['',5.244596,0.1503733,5.91];
ybs[7635]=['63 Sgr',5.2496119,-0.2369945,5.71];
ybs[7636]=['62 Sgr',5.2530687,-0.4826053,4.58];
ybs[7637]=['',5.2351115,0.9095381,6.15];
ybs[7638]=['',5.2574273,-0.6611648,4.77];
ybs[7639]=['15 Vul',5.2444738,0.4853993,4.64];
ybs[7640]=['',5.2304511,1.1098658,5.96];
ybs[7641]=['',5.2447476,0.648506,6.2];
ybs[7642]=['',5.2473925,0.4338588,5.88];
ybs[7643]=['16 Vul',5.2486025,0.4362656,5.22];
ybs[7644]=['',5.2575926,-0.3933383,6.45];
ybs[7645]=['',5.2605161,-0.5584563,4.99];
ybs[7646]=['26 Cyg',5.2444494,0.8755008,5.05];
ybs[7647]=['',5.2583354,-0.129341,6.72];
ybs[7648]=['',5.2542822,0.3239195,5.96];
ybs[7649]=['',5.2808084,-1.1570432,6.45];
ybs[7650]=['',5.2553522,0.2808262,5.67];
ybs[7651]=['δ Pav',5.282453,-1.1540248,3.56];
ybs[7652]=['',5.2300332,1.2291208,6.33];
ybs[7653]=['62 Aql',5.2597329,-0.0113494,5.68];
ybs[7654]=['',5.2658243,-0.5749165,6.53];
ybs[7655]=['τ Aql',5.2584105,0.1280567,5.52];
ybs[7656]=['',5.2553992,0.5228215,5.71];
ybs[7657]=['',5.2631296,-0.2014103,6.34];
ybs[7658]=['15 Sge',5.2579519,0.2989579,5.8];
ybs[7659]=['ξ Tel',5.2750312,-0.9218879,4.94];
ybs[7660]=['',5.2760768,-0.9591587,6.26];
ybs[7661]=['65 Sgr',5.2646894,-0.22001,6.55];
ybs[7662]=['64 Dra',5.2433486,1.1323496,5.27];
ybs[7663]=['',5.2615566,0.4061323,6.45];
ybs[7664]=['',5.2595764,0.5633542,5.64];
ybs[7665]=['η Sge',5.2624601,0.3499486,5.1];
ybs[7666]=['',5.2638439,0.271571,6.34];
ybs[7667]=['',5.2677593,-0.070134,6.47];
ybs[7668]=['65 Dra',5.2471544,1.1290982,6.57];
ybs[7669]=['',5.2617484,0.6726109,6.19];
ybs[7670]=['',5.2582102,0.8427993,6.16];
ybs[7671]=['ρ Dra',5.2486516,1.1856354,4.51];
ybs[7672]=['69 Dra',5.2316672,1.3358429,6.2];
ybs[7673]=['',5.2607089,0.9058052,6.14];
ybs[7674]=['17 Vul',5.2698946,0.4132004,5.07];
ybs[7675]=['27 Cyg',5.2671058,0.6288849,5.36];
ybs[7676]=['64 Aql',5.2756287,-0.0107795,5.99];
ybs[7677]=['',5.2917143,-1.0028968,6.37];
ybs[7678]=['',5.261396,0.9843804,6.21];
ybs[7679]=['',5.274493,0.1651141,6.43];
ybs[7680]=['',5.2780529,-0.1745651,6.18];
ybs[7681]=['',5.2577885,1.116133,6.26];
ybs[7682]=['',5.2754456,0.2919044,6.42];
ybs[7683]=['',5.2654888,0.9289633,5.85];
ybs[7684]=['',5.3627872,-1.4528518,6.17];
ybs[7685]=['',5.2729674,0.6018519,6.11];
ybs[7686]=['',5.2779441,0.188265,6.31];
ybs[7687]=['66 Dra',5.261574,1.0830652,5.39];
ybs[7688]=['',5.2699171,0.8777161,6.54];
ybs[7689]=['',5.2906515,-0.6289994,5.32];
ybs[7690]=['',5.2576738,1.1883318,6.28];
ybs[7691]=['θ Sge',5.2833158,0.3661137,6.48];
ybs[7692]=['',5.2962095,-0.7455687,6.22];
ybs[7693]=['',5.3069167,-1.1057059,6.09];
ybs[7694]=['28 Cyg',5.2804431,0.6440433,4.93];
ybs[7695]=['',5.2895624,-0.153243,6.49];
ybs[7696]=['θ Aql',5.2899226,-0.0132522,3.23];
ybs[7697]=['18 Vul',5.2857904,0.4706438,5.52];
ybs[7698]=['ξ1 Cap',5.2931449,-0.2152011,6.34];
ybs[7699]=['',5.2881755,0.3699518,6.22];
ybs[7700]=['',5.3051916,-0.9142403,5.65];
ybs[7701]=['ξ2 Cap',5.2951875,-0.2191247,5.85];
ybs[7702]=['',5.289431,0.3828839,6.26];
ybs[7703]=['',5.2954657,0.0162337,6.27];
ybs[7704]=['19 Vul',5.2912147,0.4689899,5.49];
ybs[7705]=['20 Vul',5.2921507,0.4632319,5.92];
ybs[7706]=['66 Aql',5.2983346,-0.0165204,5.47];
ybs[7707]=['',5.2913701,0.8342545,6.92];
ybs[7708]=['',5.3081291,-0.4706978,5.73];
ybs[7709]=['',5.2994992,0.4241487,6.56];
ybs[7710]=['ρ Aql',5.3024147,0.2663513,4.95];
ybs[7711]=['',5.3106528,-0.5225736,6.3];
ybs[7712]=['',5.2931497,0.8993003,6.01];
ybs[7713]=['68 Dra',5.287948,1.0845592,5.75];
ybs[7714]=['',5.3133235,-0.6351286,6.39];
ybs[7715]=['',5.3134722,-0.6131953,6.53];
ybs[7716]=['30 Cyg',5.2968307,0.8181871,4.83];
ybs[7717]=['21 Vul',5.3018045,0.5019216,5.18];
ybs[7718]=['',5.3269241,-1.102449,6.27];
ybs[7719]=['',5.2988368,0.7582091,6.14];
ybs[7720]=['',5.3007887,0.6399854,6.45];
ybs[7721]=['31 Cyg',5.2982842,0.8168902,3.79];
ybs[7722]=['29 Cyg',5.3027448,0.6434988,4.97];
ybs[7723]=['',5.3017299,0.7359512,6.71];
ybs[7724]=['3 Cap',5.3124026,-0.2141996,6.32];
ybs[7725]=['',5.3063856,0.4477754,4.78];
ybs[7726]=['33 Cyg',5.2965118,0.9883908,4.3];
ybs[7727]=['22 Vul',5.3074987,0.4114162,5.15];
ybs[7728]=['',5.2963503,1.0594743,5.79];
ybs[7729]=['',5.306641,0.5898023,5.66];
ybs[7730]=['23 Vul',5.3085027,0.4865641,4.52];
ybs[7731]=['',5.3249781,-0.831571,6.31];
ybs[7732]=['18 Sge',5.3111647,0.3780864,6.13];
ybs[7733]=['α1 Cap',5.3179405,-0.2171818,4.24];
ybs[7734]=['4 Cap',5.3198685,-0.3795236,5.87];
ybs[7735]=['',5.3265605,-0.8292898,6.13];
ybs[7736]=['κ Cep',5.2715911,1.357379,4.39];
ybs[7737]=['32 Cyg',5.3062611,0.8338866,3.98];
ybs[7738]=['',5.3093004,0.6800113,6.27];
ybs[7739]=['24 Vul',5.3130536,0.4317148,5.32];
ybs[7740]=['α2 Cap',5.3197157,-0.217814,3.57];
ybs[7741]=['',5.3071916,0.8778414,6.31];
ybs[7742]=['',5.3087437,0.7966277,5.91];
ybs[7743]=['',5.3111909,0.647876,6.48];
ybs[7744]=['',5.3325099,-0.9596663,6.27];
ybs[7745]=['',5.313007,0.7056253,5.24];
ybs[7746]=['',5.3161226,0.5098574,6.22];
ybs[7747]=['σ Cap',5.3257561,-0.3325403,5.28];
ybs[7748]=['',5.3153465,0.7467657,6.29];
ybs[7749]=['34 Cyg',5.3168962,0.6649314,4.81];
ybs[7750]=['',5.3307897,-0.5084374,6.3];
ybs[7751]=['',5.3327774,-0.6214687,6.46];
ybs[7752]=['',5.3371113,-0.8714951,6.27];
ybs[7753]=['',5.3182079,0.7120431,5.84];
ybs[7754]=['',5.326654,-0.017681,6.06];
ybs[7755]=['36 Cyg',5.3199565,0.6469062,5.58];
ybs[7756]=['35 Cyg',5.3208089,0.6117004,5.17];
ybs[7757]=['',5.3252146,0.2318215,6.21];
ybs[7758]=['',5.3299172,-0.1098825,6.63];
ybs[7759]=['ν Cap',5.3311002,-0.2215382,4.76];
ybs[7760]=['',5.3274596,0.2376042,5.95];
ybs[7761]=['',5.3316554,-0.2568948,6.1];
ybs[7762]=['β Cap',5.3326803,-0.2568301,3.08];
ybs[7763]=['',5.3210039,0.8096215,6.45];
ybs[7764]=['',5.3289054,0.2554282,6.13];
ybs[7765]=['κ1 Sgr',5.3400209,-0.7327413,5.59];
ybs[7766]=['',5.3288716,0.3116958,5.8];
ybs[7767]=['',5.318538,0.9679969,5.76];
ybs[7768]=['',5.3257516,0.6492283,6.57];
ybs[7769]=['',5.3131742,1.1679459,5.93];
ybs[7770]=['',5.3276046,0.6888649,6.23];
ybs[7771]=['',5.3955674,-1.4118594,5.77];
ybs[7772]=['',5.3258146,0.818613,6.5];
ybs[7773]=['κ2 Sgr',5.3462657,-0.7392423,5.64];
ybs[7774]=['',5.3412604,-0.1673389,6.3];
ybs[7775]=['25 Vul',5.3360754,0.4278253,5.54];
ybs[7776]=['α Pav',5.3549183,-0.9890246,1.94];
ybs[7777]=['',5.3278392,0.9365767,6.18];
ybs[7778]=['71 Dra',5.3230355,1.0877394,5.72];
ybs[7779]=['',5.3399485,0.2551358,6.17];
ybs[7780]=['',5.3415549,0.0944225,5.31];
ybs[7781]=['',5.3353674,0.7190378,6.39];
ybs[7782]=['γ Cyg',5.3361894,0.7037723,2.2];
ybs[7783]=['',5.3383024,0.5468412,6.09];
ybs[7784]=['',5.335291,0.8004332,5.58];
ybs[7785]=['',5.3545168,-0.7108431,6.09];
ybs[7786]=['',5.3384519,0.7172053,5.93];
ybs[7787]=['',5.3524781,-0.4990841,5.85];
ybs[7788]=['',5.3390943,0.7513664,6.2];
ybs[7789]=['',5.3479881,0.0198297,6.15];
ybs[7790]=['',5.3240454,1.2033311,5.55];
ybs[7791]=['',5.3296973,1.1178185,5.69];
ybs[7792]=['39 Cyg',5.3436735,0.5629942,4.43];
ybs[7793]=['',5.3429209,0.6552581,5.9];
ybs[7794]=['',5.3590896,-0.6516109,6.25];
ybs[7795]=['',5.3528232,-0.0476876,6.11];
ybs[7796]=['',5.3525703,0.1767035,6.33];
ybs[7797]=['',5.3519719,0.3748558,5.66];
ybs[7798]=['',5.4176487,-1.4174786,5.91];
ybs[7799]=['',5.3535271,0.3478979,6.41];
ybs[7800]=['π Cap',5.3603006,-0.3166554,5.25];
ybs[7801]=['',5.3454692,0.9358344,6.51];
ybs[7802]=['',5.3552069,0.3034043,6.22];
ybs[7803]=['',5.3672732,-0.6200558,6.1];
ybs[7804]=['',5.347287,1.0413961,6.44];
ybs[7805]=['',5.3663614,-0.2735361,6.41];
ybs[7806]=['',5.3630542,0.1484653,6.25];
ybs[7807]=['68 Aql',5.3646481,-0.0573989,6.13];
ybs[7808]=['ρ Cap',5.3670054,-0.3096974,4.78];
ybs[7809]=['',5.3578524,0.6003476,6.39];
ybs[7810]=['',5.3638846,0.0524637,6.21];
ybs[7811]=['',5.3700344,-0.3895949,6.16];
ybs[7812]=['40 Cyg',5.3596118,0.6721078,5.62];
ybs[7813]=['',5.3532974,0.9897244,6.36];
ybs[7814]=['43 Cyg',5.356688,0.8630961,5.69];
ybs[7815]=['ο Cap',5.371456,-0.3231827,6.74];
ybs[7816]=['ο Cap',5.3715576,-0.3231245,5.94];
ybs[7817]=['69 Aql',5.3700234,-0.0491485,4.91];
ybs[7818]=['',5.3764767,-0.5068856,6.39];
ybs[7819]=['',5.3680637,0.3518092,6.55];
ybs[7820]=['41 Cyg',5.3679042,0.5312436,4.01];
ybs[7821]=['42 Cyg',5.3674222,0.6374657,5.88];
ybs[7822]=['1 Del',5.3724765,0.1913863,6.08];
ybs[7823]=['',5.3765697,-0.2615597,6.12];
ybs[7824]=['',5.4011925,-1.2136837,6.11];
ybs[7825]=['',5.3751046,0.360862,6.18];
ybs[7826]=['',5.3764668,0.1977581,7.11];
ybs[7827]=['',5.3698321,0.8028204,6.41];
ybs[7828]=['',5.3847343,-0.4341165,6.36];
ybs[7829]=['',5.366725,0.9797828,5.91];
ybs[7830]=['ω1 Cyg',5.3699117,0.855583,4.95];
ybs[7831]=['',5.3821928,-0.1707402,5.65];
ybs[7832]=['ν Mic',5.3900964,-0.7757085,5.11];
ybs[7833]=['44 Cyg',5.374593,0.645874,6.19];
ybs[7834]=['φ1 Pav',5.3985564,-1.0560931,4.76];
ybs[7835]=['',5.3793086,0.4516017,6.34];
ybs[7836]=['θ Cep',5.3665386,1.1006665,4.22];
ybs[7837]=['ω2 Cyg',5.375378,0.8602797,5.44];
ybs[7838]=['ε Del',5.3851798,0.1985184,4.03];
ybs[7839]=['',5.3941748,-0.6635403,6.44];
ybs[7840]=['',5.375333,0.9142008,6.18];
ybs[7841]=['',5.3901546,-0.2382334,6.13];
ybs[7842]=['',5.3932837,-0.5306152,6.4];
ybs[7843]=['',5.3881826,0.1768179,6.56];
ybs[7844]=['η Del',5.3883459,0.2286109,5.38];
ybs[7845]=['ρ Pav',5.4074079,-1.0726311,4.88];
ybs[7846]=['',5.3768169,0.9922244,6.14];
ybs[7847]=['',5.3825327,0.7550716,6.6];
ybs[7848]=['',5.3890547,0.3675064,6.48];
ybs[7849]=['μ1 Oct',5.4301538,-1.3282994,6];
ybs[7850]=['μ2 Oct',5.4284058,-1.3138156,6.55];
ybs[7851]=['',5.3960882,-0.2871759,6.19];
ybs[7852]=['47 Cyg',5.3874014,0.6164852,4.61];
ybs[7853]=['',5.3866883,0.730304,6.49];
ybs[7854]=['',5.3664897,1.2671287,6.27];
ybs[7855]=['α Ind',5.4061491,-0.8241215,3.11];
ybs[7856]=['',5.3868998,0.8162038,5.78];
ybs[7857]=['ζ Del',5.3942271,0.2573646,4.68];
ybs[7858]=['',5.41745,-1.0966626,6.22];
ybs[7859]=['70 Aql',5.4008902,-0.0432438,4.89];
ybs[7860]=['26 Vul',5.3975047,0.4529923,6.41];
ybs[7861]=['φ2 Pav',5.4179839,-1.0554913,5.12];
ybs[7862]=['',5.3906108,0.9062736,6.11];
ybs[7863]=['',5.4065218,-0.4369626,6.36];
ybs[7864]=['',5.4033424,0.0029579,6.22];
ybs[7865]=['73 Dra',5.3721793,1.3094287,5.2];
ybs[7866]=['27 Vul',5.4015893,0.4631119,5.59];
ybs[7867]=['υ Pav',5.4271303,-1.1638969,5.15];
ybs[7868]=['β Del',5.4040137,0.2560029,3.63];
ybs[7869]=['ι Del',5.4052759,0.1998488,5.43];
ybs[7870]=['71 Aql',5.4078833,-0.0180179,4.32];
ybs[7871]=['48 Cyg',5.4033816,0.5523108,6.32];
ybs[7872]=['',5.4054779,0.3201267,6.25];
ybs[7873]=['',5.4034418,0.5514285,6.49];
ybs[7874]=['',5.4025138,0.670226,6.2];
ybs[7875]=['τ Cap',5.4123412,-0.2597297,5.22];
ybs[7876]=['',5.4117651,-0.0408322,6.22];
ybs[7877]=['29 Vul',5.4080594,0.3713028,4.82];
ybs[7878]=['θ Del',5.4092077,0.2336658,5.72];
ybs[7879]=['',5.4175474,-0.5822107,5.47];
ybs[7880]=['28 Vul',5.4080115,0.4221791,5.04];
ybs[7881]=['',5.4082583,0.4145776,5.91];
ybs[7882]=['κ Del',5.4110362,0.1773137,5.05];
ybs[7883]=['1 Aqr',5.4125386,0.0097691,5.16];
ybs[7884]=['',5.416625,-0.4136469,6.37];
ybs[7885]=['',5.4106708,0.2777036,5.97];
ybs[7886]=['υ Cap',5.4158209,-0.3152939,5.1];
ybs[7887]=['75 Dra',5.3532012,1.4222913,5.46];
ybs[7888]=['',5.4184841,-0.4637546,6.51];
ybs[7889]=['',5.4108998,0.3820603,6.08];
ybs[7890]=['',5.4098093,0.5307124,5.68];
ybs[7891]=['',5.4179068,-0.280132,5.8];
ybs[7892]=['α Del',5.4130912,0.2789969,3.77];
ybs[7893]=['',5.4142047,0.1976274,6.42];
ybs[7894]=['74 Dra',5.358947,1.4165162,5.96];
ybs[7895]=['',5.4221223,-0.5502012,5.76];
ybs[7896]=['',5.4219586,-0.4524919,6.28];
ybs[7897]=['',5.4118545,0.7095248,6.06];
ybs[7898]=['',5.4108562,0.7983171,6.58];
ybs[7899]=['β Pav',5.4401186,-1.1541422,3.42];
ybs[7900]=['',5.4178401,0.3492247,6.45];
ybs[7901]=['',5.4289283,-0.6891241,6.29];
ybs[7902]=['',5.4084845,0.9787471,6.48];
ybs[7903]=['',5.4168656,0.5214874,6.08];
ybs[7904]=['10 Del',5.4202511,0.2558143,5.99];
ybs[7905]=['',5.4138849,0.759779,5.95];
ybs[7906]=['η Ind',5.4346361,-0.9048823,4.51];
ybs[7907]=['49 Cyg',5.4186964,0.5651575,5.51];
ybs[7908]=['',5.418263,0.6834031,6.51];
ybs[7909]=['',5.4232223,0.3071025,6.22];
ybs[7910]=['α Cyg',5.4198191,0.7915821,1.25];
ybs[7911]=['',5.4137103,1.0573001,6.01];
ybs[7912]=['',5.4222267,0.7293937,5.67];
ybs[7913]=['',5.4243692,0.6201247,6.66];
ybs[7914]=['δ Del',5.4297889,0.2644051,4.43];
ybs[7915]=['51 Cyg',5.4229222,0.8798958,5.39];
ybs[7916]=['',5.3530133,1.4607397,6.19];
ybs[7917]=['',5.4386492,-0.4742349,6.5];
ybs[7918]=['',5.4288769,0.6224294,6.47];
ybs[7919]=['',5.4439571,-0.6828281,5.5];
ybs[7920]=['σ Pav',5.4594616,-1.199027,5.41];
ybs[7921]=['',5.4437209,-0.6290917,6.49];
ybs[7922]=['ψ Cap',5.4424007,-0.4397348,4.14];
ybs[7923]=['17 Cap',5.4425975,-0.3741682,5.93];
ybs[7924]=['',5.4240353,1.0589933,6.15];
ybs[7925]=['30 Vul',5.435669,0.4423699,4.91];
ybs[7926]=['',5.4268439,0.9981337,6.32];
ybs[7927]=['',5.438484,0.3170544,6.38];
ybs[7928]=['52 Cyg',5.4389213,0.5374807,4.22];
ybs[7929]=['ι Mic',5.4535557,-0.7664057,5.11];
ybs[7930]=['',5.4318903,0.9872136,5.78];
ybs[7931]=['4 Cep',5.4254961,1.1646951,5.58];
ybs[7932]=['',5.4459735,-0.0420703,6.27];
ybs[7933]=['γ1 Del',5.4436689,0.2827518,5.14];
ybs[7934]=['γ2 Del',5.4437271,0.2827471,4.27];
ybs[7935]=['ε Cyg',5.4412045,0.5942171,2.46];
ybs[7936]=['ε Aqr',5.4488473,-0.1643993,3.77];
ybs[7937]=['3 Aqr',5.4489918,-0.0864167,4.42];
ybs[7938]=['ζ Ind',5.4580136,-0.8054655,4.89];
ybs[7939]=['13 Del',5.4490042,0.1062001,5.58];
ybs[7940]=['',5.4490381,0.0590471,6.4];
ybs[7941]=['',5.4361002,1.006273,4.51];
ybs[7942]=['',5.4454181,0.6012726,4.92];
ybs[7943]=['η Cep',5.4353904,1.0806088,3.43];
ybs[7944]=['',5.4425149,0.8134571,6.3];
ybs[7945]=['',5.4686669,-1.0882329,6.28];
ybs[7946]=['',5.4687032,-1.0882329,6.59];
ybs[7947]=['',5.4563634,-0.448625,5.86];
ybs[7948]=['',5.4408496,0.9262662,6.33];
ybs[7949]=['λ Cyg',5.4463327,0.6382168,4.53];
ybs[7950]=['',5.4563431,-0.3134395,6.21];
ybs[7951]=['α Mic',5.4595787,-0.5882178,4.9];
ybs[7952]=['',5.4456467,0.7968471,6.4];
ybs[7953]=['',5.4308914,1.2187119,6.41];
ybs[7954]=['ι Ind',5.4671249,-0.8993753,5.05];
ybs[7955]=['',5.447598,0.8361588,5.57];
ybs[7956]=['',5.4630766,-0.5581008,6.36];
ybs[7957]=['',5.4642871,-0.6603561,5.52];
ybs[7958]=['',5.4475642,0.9160127,6.27];
ybs[7959]=['15 Del',5.4567963,0.2203029,5.98];
ybs[7960]=['14 Del',5.4576839,0.1386033,6.33];
ybs[7961]=['',5.4585297,0.0981226,6.21];
ybs[7962]=['',5.4621022,-0.2175978,5.88];
ybs[7963]=['55 Cyg',5.4525727,0.8061852,4.84];
ybs[7964]=['',5.4512333,0.9073495,6.29];
ybs[7965]=['β Mic',5.4683264,-0.5776894,6.04];
ybs[7966]=['ω Cap',5.467424,-0.4684669,4.11];
ybs[7967]=['',5.4609514,0.3164089,6.52];
ybs[7968]=['4 Aqr',5.4651132,-0.0968407,5.99];
ybs[7969]=['',5.4567961,0.8157374,6.33];
ybs[7970]=['56 Cyg',5.4576742,0.770331,5.04];
ybs[7971]=['5 Aqr',5.4682367,-0.0947514,5.55];
ybs[7972]=['β Ind',5.4820603,-1.0188365,3.65];
ybs[7973]=['',5.4759527,-0.6934423,5.35];
ybs[7974]=['',5.4643681,0.4944231,5.77];
ybs[7975]=['',5.4725521,-0.4137239,6.33];
ybs[7976]=['μ Aqr',5.4705439,-0.1554226,4.73];
ybs[7977]=['',5.4745034,-0.5347695,6.35];
ybs[7978]=['',5.4804919,-0.8839873,6.24];
ybs[7979]=['',5.4526121,1.1190899,6.45];
ybs[7980]=['',5.4725308,-0.2006286,6.38];
ybs[7981]=['31 Vul',5.46728,0.474293,4.59];
ybs[7982]=['',5.4665464,0.5746873,6.44];
ybs[7983]=['',5.4774484,-0.4860172,6.41];
ybs[7984]=['',5.4762484,-0.1188739,6.44];
ybs[7985]=['',5.4715411,0.5188487,6.34];
ybs[7986]=['19 Cap',5.4801386,-0.3114365,5.78];
ybs[7987]=['57 Cyg',5.4714861,0.7760719,4.78];
ybs[7988]=['76 Dra',5.4147702,1.4417318,5.75];
ybs[7989]=['',5.4717228,0.7899428,5.45];
ybs[7990]=['',5.4724255,0.741569,6.66];
ybs[7991]=['',5.4747938,0.5849724,5.47];
ybs[7992]=['',5.4811774,-0.0225876,6.55];
ybs[7993]=['',5.4770368,0.499178,6.56];
ybs[7994]=['32 Vul',5.4778669,0.4910732,5.01];
ybs[7995]=['',5.4765726,0.7117782,6.7];
ybs[7996]=['',5.4834002,0.0804968,6.05];
ybs[7997]=['17 Del',5.4828731,0.2408677,5.17];
ybs[7998]=['16 Del',5.4830417,0.2207481,5.58];
ybs[7999]=['',5.4890607,-0.4575662,5.7];
ybs[8000]=['',5.486337,-0.060769,6.57];
ybs[8001]=['7 Aqr',5.4890876,-0.1678605,5.51];
ybs[8002]=['',5.4391119,1.4072283,5.39];
ybs[8003]=['',5.4900399,0.0094858,6.05];
ybs[8004]=['',5.4926409,-0.2784077,5.87];
ybs[8005]=['',5.5122575,-1.1890611,6.37];
ybs[8006]=['',5.4826163,0.8289811,5.67];
ybs[8007]=['α Oct',5.5288849,-1.3428761,5.15];
ybs[8008]=['',5.4840837,0.892814,6.63];
ybs[8009]=['',5.4860163,0.7854788,5.96];
ybs[8010]=['',5.497056,-0.2513682,6.01];
ybs[8011]=['',5.4850234,0.8867697,5.81];
ybs[8012]=['',5.4851451,0.8600178,5.9];
ybs[8013]=['',5.5056802,-0.8933327,5.76];
ybs[8014]=['ν Cyg',5.4887915,0.719897,3.94];
ybs[8015]=['',5.4839594,0.9942614,6.23];
ybs[8016]=['18 Del',5.4952565,0.190581,5.48];
ybs[8017]=['',5.5033736,-0.6291702,6.11];
ybs[8018]=['33 Vul',5.4942548,0.39106,5.31];
ybs[8019]=['20 Cap',5.501118,-0.3308187,6.25];
ybs[8020]=['ε Equ',5.4982213,0.0763437,5.23];
ybs[8021]=['',5.4936704,0.7775774,5.55];
ybs[8022]=['',5.4946113,0.7333976,6.16];
ybs[8023]=['',5.5012762,0.2950476,6.66];
ybs[8024]=['',5.5024688,0.1325977,5.99];
ybs[8025]=['γ Mic',5.5088877,-0.5615842,4.67];
ybs[8026]=['',5.4941163,0.8821332,5.61];
ybs[8027]=['11 Aqr',5.5049445,-0.0811436,6.21];
ybs[8028]=['',5.5133096,-0.7490994,6.64];
ybs[8029]=['',5.4736476,1.3265259,6.05];
ybs[8030]=['',5.5038928,0.3387766,5.65];
ybs[8031]=['',5.5107346,-0.467741,6.05];
ybs[8032]=['',5.5141809,-0.6710576,5.94];
ybs[8033]=['59 Cyg',5.500076,0.8308095,4.74];
ybs[8034]=['ζ Mic',5.5164223,-0.6728192,5.3];
ybs[8035]=['',5.4974629,1.038806,5.51];
ybs[8036]=['',5.5169303,-0.4825822,6.25];
ybs[8037]=['',5.5066545,0.6301928,5.97];
ybs[8038]=['',5.5462295,-1.3286908,6.58];
ybs[8039]=['60 Cyg',5.5060688,0.8069892,5.37];
ybs[8040]=['',5.5154538,-0.0147094,6.5];
ybs[8041]=['μ Ind',5.5271565,-0.9537254,5.16];
ybs[8042]=['',5.5156417,0.0281678,6.25];
ybs[8043]=['',5.5152268,0.2585169,6.31];
ybs[8044]=['12 Aqr',5.5202992,-0.1001997,7.31];
ybs[8045]=['12 Aqr',5.5203064,-0.1001949,5.89];
ybs[8046]=['η Cap',5.5220884,-0.3450963,4.84];
ybs[8047]=['',5.5476926,-1.2756397,5.68];
ybs[8048]=['',5.5114676,0.7831778,6.19];
ybs[8049]=['',5.5146995,0.6761304,6.07];
ybs[8050]=['',5.5131925,0.801642,6.48];
ybs[8051]=['',5.5096378,0.9904966,5.83];
ybs[8052]=['3 Equ',5.5222147,0.097481,5.61];
ybs[8053]=['',5.5227856,0.0527868,6.42];
ybs[8054]=['',5.5230709,0.0410547,6.33];
ybs[8055]=['η Mic',5.5316064,-0.7208725,5.53];
ybs[8056]=['δ Mic',5.5294307,-0.5243317,5.68];
ybs[8057]=['',5.5180202,0.727981,6.33];
ybs[8058]=['',5.5156784,0.8802387,6.37];
ybs[8059]=['',5.5424241,-1.1143046,5.76];
ybs[8060]=['',5.5171321,0.8193285,6.32];
ybs[8061]=['θ Cap',5.528742,-0.2993207,4.07];
ybs[8062]=['',5.5312104,-0.5630175,5.18];
ybs[8063]=['4 Equ',5.5259857,0.1054371,5.94];
ybs[8064]=['',5.5170658,0.9314516,5.9];
ybs[8065]=['ξ Cyg',5.5225631,0.7681251,3.72];
ybs[8066]=['24 Cap',5.5341081,-0.4349789,4.5];
ybs[8067]=['',5.5560939,-1.2646522,6.2];
ybs[8068]=['',5.5295656,0.4713701,6.12];
ybs[8069]=['',5.5365837,-0.3031934,6.17];
ybs[8070]=['',5.5299274,0.5457266,5.82];
ybs[8071]=['61 Cyg',5.5314269,0.6776951,5.21];
ybs[8072]=['61 Cyg',5.5314779,0.6776515,6.03];
ybs[8073]=['χ Cap',5.5402444,-0.3684347,5.3];
ybs[8074]=['',5.5349765,0.2747514,6.34];
ybs[8075]=['63 Cyg',5.529677,0.8330709,4.55];
ybs[8076]=['',5.5391621,0.1234514,6.15];
ybs[8077]=['27 Cap',5.5445446,-0.3573073,6.25];
ybs[8078]=['ο Pav',5.5642233,-1.2224429,5.02];
ybs[8079]=['ν Aqr',5.5445071,-0.1970035,4.51];
ybs[8080]=['',5.5393189,0.5286544,5.59];
ybs[8081]=['',5.5458111,0.0528424,6.45];
ybs[8082]=['',5.5496347,-0.1617798,6.27];
ybs[8083]=['γ Equ',5.5472551,0.1783044,4.69];
ybs[8084]=['6 Equ',5.5480354,0.1768608,6.07];
ybs[8085]=['',5.5261627,1.2481702,5.87];
ybs[8086]=['',5.556843,-0.701349,5.83];
ybs[8087]=['',5.5477772,0.3933832,6.68];
ybs[8088]=['',5.5537141,-0.2511061,6.48];
ybs[8089]=['',5.5445337,0.795639,6.63];
ybs[8090]=['',5.5604053,-0.6866108,5.26];
ybs[8091]=['',5.5496687,0.6350173,6.54];
ybs[8092]=['',5.5453186,0.9363285,5.73];
ybs[8093]=['',5.5467847,0.8338552,6.46];
ybs[8094]=['',5.5614419,-0.6342252,5.96];
ybs[8095]=['',5.5410592,1.1061826,6.54];
ybs[8096]=['',5.5610425,-0.4805589,5.42];
ybs[8097]=['',5.5869802,-1.3135251,6.63];
ybs[8098]=['',5.5197518,1.3650032,5.91];
ybs[8099]=['',5.5405412,1.1968475,7.33];
ybs[8100]=['',5.5729032,-0.9281096,5.75];
ybs[8101]=['ζ Cyg',5.5580457,0.529048,3.2];
ybs[8102]=['',5.56081,0.2804388,6.27];
ybs[8103]=['',5.5700057,-0.7054672,6.21];
ybs[8104]=['',5.5649176,-0.1836004,6.77];
ybs[8105]=['',5.5515628,1.0484407,5.64];
ybs[8106]=['',5.5600236,0.6408681,6.05];
ybs[8107]=['',5.5661457,0.003108,6.38];
ybs[8108]=['',5.5687092,-0.3012259,6.04];
ybs[8109]=['δ Equ',5.565327,0.1761516,4.49];
ybs[8110]=['',5.5721782,-0.6304926,6.12];
ybs[8111]=['',5.5836506,-1.1273886,6.31];
ybs[8112]=['',5.5634491,0.5233682,6.17];
ybs[8113]=['φ Cap',5.5710685,-0.3589351,5.24];
ybs[8114]=['29 Cap',5.5714394,-0.2632856,5.28];
ybs[8115]=['',5.6552272,-1.478612,6.45];
ybs[8116]=['τ Cyg',5.5658841,0.665519,3.72];
ybs[8117]=['α Equ',5.5712931,0.0930962,3.92];
ybs[8118]=['',5.5751,-0.0265509,6.48];
ybs[8119]=['',5.559431,1.1255513,6.39];
ybs[8120]=['',5.5778344,-0.2302468,6.4];
ybs[8121]=['ε Mic',5.5814524,-0.5599983,4.71];
ybs[8122]=['',5.5690452,0.8388007,6.46];
ybs[8123]=['30 Cap',5.5811328,-0.3123846,5.43];
ybs[8124]=['',5.5731537,0.7389342,6.43];
ybs[8125]=['31 Cap',5.5824573,-0.3032539,7.05];
ybs[8126]=['θ Ind',5.5907754,-0.9313444,4.39];
ybs[8127]=['15 Aqr',5.5818195,-0.0773605,5.82];
ybs[8128]=['',5.5855708,-0.5005305,6.4];
ybs[8129]=['σ Cyg',5.5773011,0.6890813,4.23];
ybs[8130]=['',5.577039,0.7464781,6.19];
ybs[8131]=['',5.5916028,-0.7842555,6];
ybs[8132]=['υ Cyg',5.5796528,0.6105832,4.43];
ybs[8133]=['',5.5749044,0.9439452,6.13];
ybs[8134]=['',5.5892373,-0.4584197,6.56];
ybs[8135]=['',5.5844375,0.1970574,5.96];
ybs[8136]=['',5.5756771,0.9753719,5.98];
ybs[8137]=['θ1 Mic',5.5940499,-0.7107302,4.82];
ybs[8138]=['',5.5966986,-0.8700417,6.38];
ybs[8139]=['',5.5757869,1.024479,6.42];
ybs[8140]=['68 Cyg',5.581646,0.7685189,5];
ybs[8141]=['',5.5838126,0.7178198,6.15];
ybs[8142]=['',5.611671,-1.2155363,6.41];
ybs[8143]=['',5.5858775,0.668895,5.83];
ybs[8144]=['',5.5901389,0.3859628,6.29];
ybs[8145]=['',5.6165271,-1.2515764,6.09];
ybs[8146]=['16 Aqr',5.5944143,-0.0780521,5.87];
ybs[8147]=['',5.5858658,0.8656424,5.76];
ybs[8148]=['α Cep',5.5809538,1.0938433,2.44];
ybs[8149]=['9 Equ',5.5941856,0.129894,5.82];
ybs[8150]=['',5.5843156,1.0246984,5.66];
ybs[8151]=['',5.5937585,0.4178972,5.57];
ybs[8152]=['',5.5924712,0.5679408,5.68];
ybs[8153]=['ι Cap',5.5998216,-0.292275,4.28];
ybs[8154]=['',5.5651928,1.3456177,5.95];
ybs[8155]=['',5.5947883,0.5707363,6.04];
ybs[8156]=['',5.5930219,0.7056966,6.4];
ybs[8157]=['6 Cep',5.5841479,1.1337526,5.18];
ybs[8158]=['',5.6032879,-0.3941009,5.6];
ybs[8159]=['1 Peg',5.5982868,0.3471929,4.08];
ybs[8160]=['',5.5520126,1.419231,6.15];
ybs[8161]=['17 Aqr',5.6026642,-0.1611098,5.99];
ybs[8162]=['',5.6594066,-1.4414828,6.38];
ybs[8163]=['',5.6099209,-0.8120316,6.31];
ybs[8164]=['β Equ',5.6021114,0.120421,5.16];
ybs[8165]=['',5.589793,1.0619392,6.11];
ybs[8166]=['θ2 Mic',5.6099768,-0.7141475,5.77];
ybs[8167]=['γ Pav',5.6204321,-1.1392881,4.22];
ybs[8168]=['',5.6006742,0.5305477,6.05];
ybs[8169]=['33 Cap',5.6082639,-0.362383,5.41];
ybs[8170]=['',5.6081883,-0.3954572,6.38];
ybs[8171]=['',5.5969166,0.8635377,5.69];
ybs[8172]=['',5.6007675,0.675837,6.63];
ybs[8173]=['18 Aqr',5.6082121,-0.2232124,5.49];
ybs[8174]=['γ Ind',5.6186744,-0.9524425,6.12];
ybs[8175]=['',5.6034394,0.6544164,6.58];
ybs[8176]=['',5.6064389,0.4252145,5.71];
ybs[8177]=['',5.6086555,0.1791257,6.35];
ybs[8178]=['20 Aqr',5.6109239,-0.0577565,6.36];
ybs[8179]=['',5.6052834,0.6534539,6.47];
ybs[8180]=['',5.6070372,0.4433328,6.15];
ybs[8181]=['19 Aqr',5.6126197,-0.1685878,5.7];
ybs[8182]=['',5.6310528,-1.2115175,5.34];
ybs[8183]=['',5.6082006,0.4296625,6.32];
ybs[8184]=['',5.6089499,0.4583838,5.68];
ybs[8185]=['21 Aqr',5.6127671,-0.0605176,5.49];
ybs[8186]=['',5.6184393,-0.658684,5.63];
ybs[8187]=['',5.6545491,-1.3953425,6.47];
ybs[8188]=['',5.6214235,-0.7410308,5.51];
ybs[8189]=['',5.6151873,0.0108888,6.46];
ybs[8190]=['ζ Cap',5.6192336,-0.3895869,3.74];
ybs[8191]=['',5.6178367,0.020821,6.13];
ybs[8192]=['',5.6096767,0.8624096,6.58];
ybs[8193]=['35 Cap',5.6217326,-0.3683732,5.78];
ybs[8194]=['',5.6115593,0.8168782,5.6];
ybs[8195]=['69 Cyg',5.613949,0.6415286,5.94];
ybs[8196]=['',5.6173247,0.3397311,6.07];
ybs[8197]=['',5.6305808,-0.9357648,6.39];
ybs[8198]=['',5.6258044,-0.2003316,6.61];
ybs[8199]=['36 Cap',5.6281867,-0.3790312,4.51];
ybs[8200]=['5 PsA',5.6299199,-0.543638,6.5];
ybs[8201]=['70 Cyg',5.6208068,0.6493765,5.31];
ybs[8202]=['',5.6181596,0.853897,5.31];
ybs[8203]=['35 Vul',5.6224553,0.4834315,5.41];
ybs[8204]=['',5.6174581,0.9248197,6.03];
ybs[8205]=['',5.6261714,0.144614,6.4];
ybs[8206]=['',5.6243763,0.5640099,5.8];
ybs[8207]=['',5.6285066,0.3140932,6.44];
ybs[8208]=['',5.6336607,-0.3326085,6.57];
ybs[8209]=['',5.6284037,0.3886818,5.93];
ybs[8210]=['',5.619911,1.0444023,6.1];
ybs[8211]=['2 Peg',5.6325152,0.4141589,4.57];
ybs[8212]=['',5.6266071,0.9688133,6.12];
ybs[8213]=['7 Cep',5.6206645,1.1676094,5.44];
ybs[8214]=['71 Cyg',5.6295885,0.8138653,5.24];
ybs[8215]=['ξ Gru',5.6434694,-0.7171173,5.29];
ybs[8216]=['6 PsA',5.6438611,-0.5908519,5.97];
ybs[8217]=['',5.6380734,0.2134283,6.08];
ybs[8218]=['β Aqr',5.6401866,-0.0956429,2.91];
ybs[8219]=['',5.6492033,-0.9188418,6.41];
ybs[8220]=['',5.6783613,-1.3849007,6.18];
ybs[8221]=['',5.6449634,-0.4275894,6.43];
ybs[8222]=['',5.6492466,-0.7811544,5.57];
ybs[8223]=['',5.6331517,0.9258763,6.02];
ybs[8224]=['β Cep',5.6239562,1.2330928,3.23];
ybs[8225]=['',5.6031787,1.4069722,5.97];
ybs[8226]=['',5.6434518,0.4099057,6.7];
ybs[8227]=['',5.653055,-0.7475766,6.32];
ybs[8228]=['',5.638061,0.9199821,6.16];
ybs[8229]=['',5.6354555,1.0568034,5.53];
ybs[8230]=['',5.6552474,-0.516686,6.41];
ybs[8231]=['37 Cap',5.6548619,-0.3489311,5.69];
ybs[8232]=['',5.644695,0.8738744,5.75];
ybs[8233]=['',5.6567515,-0.4077417,6.4];
ybs[8234]=['',5.646438,0.8019058,6.25];
ybs[8235]=['',5.6708316,-1.1297691,6.2];
ybs[8236]=['',5.6527021,0.3987514,6.47];
ybs[8237]=['',5.6564485,-0.0679068,5.77];
ybs[8238]=['ρ Cyg',5.6494356,0.7973326,4.02];
ybs[8239]=['8 PsA',5.6608227,-0.4551614,5.73];
ybs[8240]=['ν Oct',5.6883461,-1.3490659,3.76];
ybs[8241]=['72 Cyg',5.6531749,0.6741555,4.9];
ybs[8242]=['7 PsA',5.6637612,-0.5751786,6.11];
ybs[8243]=['',5.6558398,0.4937447,6.31];
ybs[8244]=['',5.6565175,0.4283829,6.11];
ybs[8245]=['',5.651229,0.9039115,6.15];
ybs[8246]=['ε Cap',5.6645723,-0.3381278,4.68];
ybs[8247]=['',5.6597861,0.5261834,6.36];
ybs[8248]=['',5.6584406,0.793552,5.53];
ybs[8249]=['',5.6662795,-0.0051895,6.25];
ybs[8250]=['ξ Aqr',5.6672527,-0.1354579,4.69];
ybs[8251]=['3 Peg',5.6668576,0.1171346,6.18];
ybs[8252]=['74 Cyg',5.662606,0.706969,5.01];
ybs[8253]=['5 Peg',5.6667102,0.3387963,5.45];
ybs[8254]=['',5.6737462,-0.586182,6.28];
ybs[8255]=['',5.6783553,-0.9122048,6.21];
ybs[8256]=['4 Peg',5.6703806,0.1023615,5.67];
ybs[8257]=['',5.6810089,-0.9711649,6.33];
ybs[8258]=['',5.6646961,0.781725,6.2];
ybs[8259]=['',5.674797,-0.1829707,6.08];
ybs[8260]=['',5.6709016,0.4466674,6.16];
ybs[8261]=['',5.6650314,0.9448363,6.15];
ybs[8262]=['',5.6721969,0.355325,5.85];
ybs[8263]=['25 Aqr',5.6749172,0.0407905,5.1];
ybs[8264]=['γ Cap',5.677644,-0.2891756,3.68];
ybs[8265]=['9 Cep',5.6656442,1.0851569,4.73];
ybs[8266]=['λ Oct',5.7326721,-1.4420308,5.29];
ybs[8267]=['',5.6705627,1.0050033,5.62];
ybs[8268]=['',5.6851466,-0.4364683,6.49];
ybs[8269]=['42 Cap',5.6839415,-0.243533,5.18];
ybs[8270]=['75 Cyg',5.6766382,0.7569064,5.11];
ybs[8271]=['41 Cap',5.6861714,-0.4043676,5.24];
ybs[8272]=['',5.7039937,-1.2376761,6.01];
ybs[8273]=['26 Aqr',5.6863395,0.0240773,5.67];
ybs[8274]=['κ Cap',5.6888851,-0.327633,4.73];
ybs[8275]=['7 Peg',5.6866457,0.1007801,5.3];
ybs[8276]=['',5.6784373,0.9593379,6.2];
ybs[8277]=['76 Cyg',5.6827865,0.7138279,6.11];
ybs[8278]=['',5.6878157,0.1905738,6.09];
ybs[8279]=['',5.6913728,-0.3407977,6.22];
ybs[8280]=['',5.990663,-1.5479453,6.57];
ybs[8281]=['44 Cap',5.6906001,-0.249673,5.88];
ybs[8282]=['',5.684907,0.6214151,6.07];
ybs[8283]=['',5.6850852,0.8004087,6.17];
ybs[8284]=['',5.6973829,-0.6712113,6.3];
ybs[8285]=['77 Cyg',5.6863164,0.7185783,5.69];
ybs[8286]=['π1 Cyg',5.6846446,0.8950731,4.67];
ybs[8287]=['45 Cap',5.6947217,-0.2557722,5.99];
ybs[8288]=['',5.701388,-0.8622527,6.45];
ybs[8289]=['',5.6871408,0.8673348,6.09];
ybs[8290]=['ι PsA',5.6992039,-0.5747506,4.34];
ybs[8291]=['',5.6894817,0.7199393,5.49];
ybs[8292]=['79 Cyg',5.6909819,0.6698306,5.65];
ybs[8293]=['ε Peg',5.6949783,0.174006,2.39];
ybs[8294]=['μ1 Cyg',5.6943735,0.5033104,4.73];
ybs[8295]=['μ2 Cyg',5.6943517,0.5033153,6.08];
ybs[8296]=['46 Cap',5.6989183,-0.1568606,5.09];
ybs[8297]=['',5.6870441,1.0361228,6.08];
ybs[8298]=['9 Peg',5.6962436,0.3044708,4.34];
ybs[8299]=['',5.6963413,0.2594754,5.94];
ybs[8300]=['κ Peg',5.6966401,0.4492465,4.13];
ybs[8301]=['μ Cep',5.6903619,1.027555,4.08];
ybs[8302]=['11 Cep',5.6819462,1.2462604,4.56];
ybs[8303]=['47 Cap',5.7044475,-0.1602288,6];
ybs[8304]=['λ Cap',5.7056373,-0.196705,5.58];
ybs[8305]=['',5.7011622,0.6274885,6.4];
ybs[8306]=['12 Peg',5.7029395,0.4021974,5.29];
ybs[8307]=['δ Cap',5.7079359,-0.2798044,2.87];
ybs[8308]=['',5.7141267,-0.8239289,5.58];
ybs[8309]=['',5.6867625,1.2638744,5.17];
ybs[8310]=['',5.7043057,0.4478295,6.28];
ybs[8311]=['θ PsA',5.7113086,-0.5376055,5.01];
ybs[8312]=['',5.6961096,1.0917994,5.95];
ybs[8313]=['11 Peg',5.7084139,0.0485509,5.64];
ybs[8314]=['',5.7032672,0.7532176,6.54];
ybs[8315]=['',5.7074519,0.3017633,6.21];
ybs[8316]=['',5.7227934,-1.1277626,5.62];
ybs[8317]=['',5.7103298,-0.1016036,6.17];
ybs[8318]=['ο Ind',5.7268056,-1.2135755,5.53];
ybs[8319]=['ν Cep',5.6986767,1.0684196,4.29];
ybs[8320]=['π2 Cyg',5.7052752,0.8622788,4.23];
ybs[8321]=['',5.7116112,0.6401244,6.47];
ybs[8322]=['',5.7194036,-0.2203781,6.31];
ybs[8323]=['',5.7130854,0.6762204,6.12];
ybs[8324]=['12 Cep',5.7073476,1.0609581,5.52];
ybs[8325]=['',5.7218118,-0.2923123,6.38];
ybs[8326]=['',5.7177301,0.3588176,6.29];
ybs[8327]=['',5.7044911,1.2260298,6.29];
ybs[8328]=['14 Peg',5.7192392,0.52832,5.04];
ybs[8329]=['13 Peg',5.7208378,0.3033728,5.29];
ybs[8330]=['',5.7181542,0.719864,6.48];
ybs[8331]=['',5.7282949,-0.3233431,6.16];
ybs[8332]=['',5.7155944,1.0710899,6.17];
ybs[8333]=['',5.7270067,0.3477301,5.77];
ybs[8334]=['',5.7243817,0.691732,6.17];
ybs[8335]=['',5.7301782,0.3729777,6.89];
ybs[8336]=['μ Cap',5.7351723,-0.2348235,5.08];
ybs[8337]=['',5.7451079,-1.0784141,5.9];
ybs[8338]=['γ Gru',5.7384686,-0.6504414,3.01];
ybs[8339]=['15 Peg',5.7308613,0.5042321,5.53];
ybs[8340]=['',5.736434,-0.1782734,6.59];
ybs[8341]=['16 Peg',5.7333953,0.4541729,5.08];
ybs[8342]=['',5.7277852,0.9755313,5.71];
ybs[8343]=['',5.7359748,0.3449761,5.68];
ybs[8344]=['',5.7377053,0.1215078,6.15];
ybs[8345]=['',5.7388239,-0.0729306,5.71];
ybs[8346]=['',5.7252628,1.1492912,6.37];
ybs[8347]=['π Ind',5.7493568,-1.0088243,6.19];
ybs[8348]=['',5.7406604,-0.0559118,6.2];
ybs[8349]=['',5.7388845,0.3458517,6.39];
ybs[8350]=['',5.7470011,-0.5324726,6.41];
ybs[8351]=['',5.749149,-0.6484863,5.46];
ybs[8352]=['',5.7520244,-0.6570937,6.18];
ybs[8353]=['δ Ind',5.7565222,-0.9580812,4.4];
ybs[8354]=['κ1 Ind',5.7592982,-1.0282359,6.12];
ybs[8355]=['',5.776758,-1.3537288,6.41];
ybs[8356]=['13 Cep',5.7402917,0.9897591,5.8];
ybs[8357]=['',5.7480658,0.3724144,6.4];
ybs[8358]=['17 Peg',5.7506011,0.2124865,5.54];
ybs[8359]=['',5.7419583,1.0758155,6.13];
ybs[8360]=['',5.7423749,1.1417701,5.86];
ybs[8361]=['',5.7565061,-0.0929597,6.33];
ybs[8362]=['',5.7500847,0.8511413,6.42];
ybs[8363]=['',5.7590166,-0.3679872,6.12];
ybs[8364]=['',5.7618977,-0.6683993,5.5];
ybs[8365]=['',5.7814293,-1.3267781,5.95];
ybs[8366]=['',5.7673787,-0.9736087,6.01];
ybs[8367]=['',5.7595197,-0.0746015,6.22];
ybs[8368]=['',5.7475018,1.112187,4.91];
ybs[8369]=['',5.7495861,1.1563559,6.43];
ybs[8370]=['18 Peg',5.7646278,0.1189704,6];
ybs[8371]=['η PsA',5.7683544,-0.494878,5.42];
ybs[8372]=['ε Ind',5.7803033,-0.9893622,4.69];
ybs[8373]=['',5.7573823,1.0960139,5.93];
ybs[8374]=['',5.7599011,1.0080517,6.59];
ybs[8375]=['28 Aqr',5.7688867,0.0122912,5.58];
ybs[8376]=['',5.7654762,0.5777943,6.46];
ybs[8377]=['20 Peg',5.7686968,0.2307142,5.6];
ybs[8378]=['19 Peg',5.7690598,0.1458479,5.65];
ybs[8379]=['',5.7740771,-0.3107401,6.28];
ybs[8380]=['',5.751003,1.3106549,6.35];
ybs[8381]=['29 Aqr',5.775128,-0.2943426,6.37];
ybs[8382]=['',5.7728091,0.1932664,6.37];
ybs[8383]=['',5.7790666,-0.5201845,7.1];
ybs[8384]=['',5.7651266,1.0923515,6.66];
ybs[8385]=['16 Cep',5.7575784,1.2789543,5.03];
ybs[8386]=['30 Aqr',5.7785633,-0.1120977,5.54];
ybs[8387]=['ο Aqr',5.7786709,-0.0358752,4.69];
ybs[8388]=['',5.7709178,0.9247035,5.78];
ybs[8389]=['21 Peg',5.778447,0.2004714,5.8];
ybs[8390]=['13 PsA',5.7839126,-0.520398,6.47];
ybs[8391]=['14 Cep',5.7716645,1.0140362,5.56];
ybs[8392]=['',5.7760605,0.781029,5.6];
ybs[8393]=['',5.7847831,-0.4663937,5.96];
ybs[8394]=['κ2 Ind',5.791301,-1.0390936,5.62];
ybs[8395]=['32 Aqr',5.7850855,-0.0140766,5.3];
ybs[8396]=['λ Gru',5.7916315,-0.6884078,4.46];
ybs[8397]=['',5.7835041,0.5766919,6.38];
ybs[8398]=['ν Peg',5.7888698,0.0900408,4.84];
ybs[8399]=['α Aqr',5.78941,-0.0038284,2.96];
ybs[8400]=['',5.7863368,0.4672963,5.78];
ybs[8401]=['18 Cep',5.7792015,1.10339,5.29];
ybs[8402]=['ξ Cep',5.7786683,1.12971,4.29];
ybs[8403]=['ι Aqr',5.7924882,-0.2403177,4.27];
ybs[8404]=['23 Peg',5.7879841,0.5072659,5.7];
ybs[8405]=['',5.8144875,-1.3225921,6.55];
ybs[8406]=['',5.7861657,0.8175985,6.13];
ybs[8407]=['',5.7887135,0.7891084,6.44];
ybs[8408]=['',5.7482175,1.448065,6.98];
ybs[8409]=['',5.7895472,0.7874027,5.14];
ybs[8410]=['α Gru',5.8010882,-0.8178637,1.74];
ybs[8411]=['20 Cep',5.7841616,1.0975624,5.27];
ybs[8412]=['',5.788647,0.843553,6.27];
ybs[8413]=['19 Cep',5.7848144,1.0887394,5.11];
ybs[8414]=['',5.790297,0.7914903,6.19];
ybs[8415]=['ι Peg',5.7943193,0.4441104,3.76];
ybs[8416]=['μ PsA',5.8013511,-0.5739972,4.5];
ybs[8417]=['',5.8198068,-1.3266986,6.15];
ybs[8418]=['υ PsA',5.8015936,-0.592415,4.99];
ybs[8419]=['',5.7899092,0.9851248,6.39];
ybs[8420]=['',5.7964623,0.3416712,5.75];
ybs[8421]=['',5.7965911,0.3159277,6.35];
ybs[8422]=['',5.8027852,-0.5763859,6.37];
ybs[8423]=['25 Peg',5.798,0.380545,5.78];
ybs[8424]=['35 Aqr',5.8036765,-0.3214651,5.81];
ybs[8425]=['',5.8086771,-0.8378601,6.43];
ybs[8426]=['',5.7998832,0.447582,6.11];
ybs[8427]=['',5.7938564,1.0287231,6.32];
ybs[8428]=['',5.7953094,0.9321445,6.14];
ybs[8429]=['',5.8081084,-0.5918999,5.37];
ybs[8430]=['',5.7991751,0.8708724,6.42];
ybs[8431]=['',5.8083031,-0.4920281,6.44];
ybs[8432]=['τ PsA',5.8090264,-0.5663058,4.92];
ybs[8433]=['',5.8011064,0.8001107,6.11];
ybs[8434]=['π1 Peg',5.8038259,0.5807301,5.58];
ybs[8435]=['θ Peg',5.8085767,0.1099413,3.53];
ybs[8436]=['',5.8094,-0.0661956,6.27];
ybs[8437]=['38 Aqr',5.8107168,-0.2000758,5.46];
ybs[8438]=['',5.8103294,-0.0727059,6.01];
ybs[8439]=['π2 Peg',5.8071457,0.5808398,4.29];
ybs[8440]=['',5.8088594,0.3441504,6.18];
ybs[8441]=['',5.8091794,0.257112,6.33];
ybs[8442]=['',5.8126937,-0.3688037,6.09];
ybs[8443]=['',5.8103358,0.2046562,5.78];
ybs[8444]=['28 Peg',5.8096498,0.3679069,6.46];
ybs[8445]=['',5.8110186,0.5350235,6.32];
ybs[8446]=['',5.8156381,0.2817367,5.95];
ybs[8447]=['39 Aqr',5.818633,-0.2459514,6.03];
ybs[8448]=['',5.8117715,0.8888076,5.4];
ybs[8449]=['',5.8211529,-0.4577256,6.17];
ybs[8450]=['ζ Cep',5.8100635,1.0175727,3.35];
ybs[8451]=['',5.8166931,0.4372368,5.92];
ybs[8452]=['',5.8197891,-0.0806143,6.39];
ybs[8453]=['24 Cep',5.8040626,1.2643575,4.79];
ybs[8454]=['λ Cep',5.8128652,1.038752,5.04];
ybs[8455]=['',5.8245388,-0.4377046,5.58];
ybs[8456]=['ψ Oct',5.8458469,-1.3510329,5.51];
ybs[8457]=['',5.8143318,0.993811,5.24];
ybs[8458]=['',5.8060695,1.260345,6.37];
ybs[8459]=['',5.8080996,1.2258183,5.5];
ybs[8460]=['',5.8193885,0.605746,5.33];
ybs[8461]=['',5.8147962,1.0329989,6.3];
ybs[8462]=['',5.8288376,-0.7204588,6.23];
ybs[8463]=['λ PsA',5.8270879,-0.4828385,5.43];
ybs[8464]=['',5.8150554,1.0622287,5.35];
ybs[8465]=['41 Aqr',5.8269055,-0.3660276,5.32];
ybs[8466]=['ε Oct',5.8564655,-1.4021286,5.1];
ybs[8467]=['',5.8232159,0.5010927,5.89];
ybs[8468]=['',5.8163619,1.1064206,5.79];
ybs[8469]=['',5.8330269,-0.7740417,6.1];
ybs[8470]=['',5.8239891,0.6949414,4.49];
ybs[8471]=['μ1 Gru',5.8330697,-0.7198443,4.79];
ybs[8472]=['',5.8235771,0.7948758,5.53];
ybs[8473]=['μ2 Gru',5.8366858,-0.7247426,5.1];
ybs[8474]=['',5.8276659,0.751474,5.71];
ybs[8475]=['',5.8227357,1.1041768,6.11];
ybs[8476]=['',5.8338382,0.1510083,6.21];
ybs[8477]=['',5.8371224,-0.4502163,6.15];
ybs[8478]=['',5.81738,1.2812314,6.08];
ybs[8479]=['ε Cep',5.8284095,0.9973869,4.19];
ybs[8480]=['',5.836451,-0.0260677,6.15];
ybs[8481]=['42 Aqr',5.8376777,-0.2221546,5.34];
ybs[8482]=['',5.8387,-0.4020729,6.17];
ybs[8483]=['1 Lac',5.8331784,0.6606343,4.13];
ybs[8484]=['θ Aqr',5.8377373,-0.1340493,4.16];
ybs[8485]=['',5.8379457,-0.1559821,5.79];
ybs[8486]=['',5.8449976,-0.93418,5.37];
ybs[8487]=['α Tuc',5.8463812,-1.0499283,2.86];
ybs[8488]=['',5.8356865,0.4870683,6.37];
ybs[8489]=['44 Aqr',5.8389005,-0.0922283,5.75];
ybs[8490]=['υ Oct',5.91245,-1.4948494,5.77];
ybs[8491]=['',5.8345624,1.0004756,5.88];
ybs[8492]=['',5.8430263,-0.0023499,6.39];
ybs[8493]=['45 Aqr',5.8473313,-0.2304124,5.95];
ybs[8494]=['',5.8553734,-1.0019289,6.34];
ybs[8495]=['',5.8461405,0.6610042,6.17];
ybs[8496]=['25 Cep',5.8419701,1.0979442,5.75];
ybs[8497]=['ρ Aqr',5.8524169,-0.1346961,5.37];
ybs[8498]=['30 Peg',5.8533592,0.1028538,5.37];
ybs[8499]=['',5.8553751,0.1446949,6.17];
ybs[8500]=['ν Ind',5.8740556,-1.2592726,5.29];
ybs[8501]=['47 Aqr',5.8587153,-0.3751489,5.13];
ybs[8502]=['',5.8553853,0.47192,6.47];
ybs[8503]=['γ Aqr',5.8586853,-0.0223984,3.84];
ybs[8504]=['',5.8532869,0.8915926,6.42];
ybs[8505]=['31 Peg',5.857883,0.214835,5.01];
ybs[8506]=['π1 Gru',5.8641978,-0.8001226,6.62];
ybs[8507]=['32 Peg',5.8567572,0.4962733,4.81];
ybs[8508]=['2 Lac',5.8550394,0.8140287,4.57];
ybs[8509]=['π2 Gru',5.8659465,-0.7997866,5.62];
ybs[8510]=['',5.8406264,1.336768,6.66];
ybs[8511]=['',5.8798984,-1.3074392,6.04];
ybs[8512]=['',5.8762467,-1.2274379,5.78];
ybs[8513]=['',5.8587572,0.7362191,6.41];
ybs[8514]=['49 Aqr',5.8671471,-0.4303672,5.53];
ybs[8515]=['',5.8669632,-0.1237467,5.93];
ybs[8516]=['',5.8742558,-1.0069264,5.32];
ybs[8517]=['33 Peg',5.8670983,0.3656924,6.04];
ybs[8518]=['51 Aqr',5.8694599,-0.0825986,5.78];
ybs[8519]=['50 Aqr',5.8710538,-0.2343102,5.76];
ybs[8520]=['',5.8632705,1.0016197,6.16];
ybs[8521]=['',5.8678238,0.6750575,6.22];
ybs[8522]=['',5.8629775,1.0912519,6.04];
ybs[8523]=['β Lac',5.8659229,0.9133905,4.43];
ybs[8524]=['π Aqr',5.8744414,0.0258679,4.66];
ybs[8525]=['δ Tuc',5.8850969,-1.1320436,4.48];
ybs[8526]=['4 Lac',5.8702043,0.8653489,4.57];
ybs[8527]=['',5.8787405,-0.4115083,6.29];
ybs[8528]=['',5.8759486,0.3237436,6.26];
ybs[8529]=['53 Aqr',5.8803399,-0.2903618,6.57];
ybs[8530]=['53 Aqr',5.8803545,-0.2903812,6.35];
ybs[8531]=['',5.8082143,1.5011746,5.27];
ybs[8532]=['',5.8909325,-1.17607,5.55];
ybs[8533]=['34 Peg',5.8802761,0.0785137,5.75];
ybs[8534]=['',5.8803388,0.65535,6.46];
ybs[8535]=['',5.8636111,1.3674227,6.76];
ybs[8536]=['35 Peg',5.8856614,0.0837878,4.79];
ybs[8537]=['ν Gru',5.8898187,-0.6811435,5.47];
ybs[8538]=['',5.8832572,0.6966439,6.14];
ybs[8539]=['',5.8807397,0.986779,6.57];
ybs[8540]=['',5.8848573,0.5575521,5.98];
ybs[8541]=['δ1 Gru',5.8926083,-0.7573008,3.97];
ybs[8542]=['',5.8753981,1.2370118,5.47];
ybs[8543]=['ζ1 Aqr',5.8899579,0.0014842,4.59];
ybs[8544]=['ζ2 Aqr',5.889987,0.0014891,4.42];
ybs[8545]=['δ2 Gru',5.8947418,-0.7617304,4.11];
ybs[8546]=['26 Cep',5.8806729,1.1386033,5.46];
ybs[8547]=['36 Peg',5.8911655,0.1611683,5.58];
ybs[8548]=['',5.8944405,-0.4712689,5.95];
ybs[8549]=['',5.8910624,0.4689425,5.79];
ybs[8550]=['',5.8953527,-0.2235671,6.4];
ybs[8551]=['37 Peg',5.8948658,0.0791891,5.48];
ybs[8552]=['56 Aqr',5.8965319,-0.2527277,6.37];
ybs[8553]=['',5.8861821,1.1203398,6.29];
ybs[8554]=['',5.8933516,0.6253696,6.56];
ybs[8555]=['ζ PsA',5.8993422,-0.4532253,6.43];
ybs[8556]=['δ Cep',5.8902081,1.0213777,3.75];
ybs[8557]=['5 Lac',5.8921801,0.8344835,4.36];
ybs[8558]=['σ Aqr',5.8980335,-0.184523,4.82];
ybs[8559]=['38 Peg',5.8947161,0.5703394,5.65];
ybs[8560]=['',5.8946587,0.8632687,6.4];
ybs[8561]=['β PsA',5.9021086,-0.5626989,4.29];
ybs[8562]=['',5.9223703,-1.3729641,6.15];
ybs[8563]=['ρ1 Cep',5.8767139,1.3769016,5.83];
ybs[8564]=['6 Lac',5.8964932,0.7544876,4.51];
ybs[8565]=['',5.9008094,-0.0489621,6.16];
ybs[8566]=['',5.9008572,-0.1125599,6.14];
ybs[8567]=['ν Tuc',5.9095417,-1.0799416,4.81];
ybs[8568]=['58 Aqr',5.9025801,-0.1884903,6.38];
ybs[8569]=['',5.9014974,0.5174657,6.35];
ybs[8570]=['α Lac',5.8998039,0.8794412,3.77];
ybs[8571]=['39 Peg',5.9061064,0.3549305,6.42];
ybs[8572]=['',5.9069981,0.2787184,6.32];
ybs[8573]=['',5.9051007,0.6961369,5.88];
ybs[8574]=['',5.904142,0.9449815,6.35];
ybs[8575]=['60 Aqr',5.9127541,-0.0256192,5.89];
ybs[8576]=['ρ2 Cep',5.8906757,1.3775814,5.5];
ybs[8577]=['υ Aqr',5.9158206,-0.3595713,5.2];
ybs[8578]=['',5.9218622,-1.0083983,6.23];
ybs[8579]=['',5.9100036,0.9901463,5.71];
ybs[8580]=['',5.9063548,1.2220739,6.6];
ybs[8581]=['',5.9198579,-0.4168637,5.97];
ybs[8582]=['η Aqr',5.9184443,-0.0001915,4.02];
ybs[8583]=['',5.907334,1.2301081,6.34];
ybs[8584]=['',5.901934,1.3322498,5.68];
ybs[8585]=['σ1 Gru',5.9239841,-0.7064401,6.28];
ybs[8586]=['',5.9242585,-0.5507759,5.82];
ybs[8587]=['σ2 Gru',5.9261262,-0.706584,5.86];
ybs[8588]=['8 Lac',5.9201005,0.6936074,5.73];
ybs[8589]=['',5.921317,0.6228012,6.1];
ybs[8590]=['',5.923751,0.2060133,6.4];
ybs[8591]=['',5.9199192,0.8757663,6.29];
ybs[8592]=['',5.9195942,0.9804666,6.38];
ybs[8593]=['',5.9257987,0.2213784,6.3];
ybs[8594]=['',5.9242865,0.6241171,6.3];
ybs[8595]=['κ Aqr',5.9289663,-0.0719268,5.03];
ybs[8596]=['',5.9358462,-0.9177817,6.65];
ybs[8597]=['',5.9316867,-0.1359737,6.23];
ybs[8598]=['9 Lac',5.9263608,0.9015,4.63];
ybs[8599]=['',5.9336039,-0.4998736,6.47];
ybs[8600]=['31 Cep',5.9177995,1.2871735,5.08];
ybs[8601]=['',5.9341683,-0.5755091,5.66];
ybs[8602]=['',5.9305525,0.7904611,6.4];
ybs[8603]=['40 Peg',5.9335584,0.342597,5.82];
ybs[8604]=['',5.9379054,-0.4924967,6.31];
ybs[8605]=['',5.9433054,-1.0003307,5.97];
ybs[8606]=['',5.9316783,0.9931432,5.21];
ybs[8607]=['10 Lac',5.9349268,0.6834269,4.88];
ybs[8608]=['',5.9407241,-0.533224,5.87];
ybs[8609]=['41 Peg',5.937515,0.3453728,6.21];
ybs[8610]=['',5.9238587,1.3173478,5.79];
ybs[8611]=['',5.9363096,0.6579896,6.03];
ybs[8612]=['30 Cep',5.9314225,1.1116268,5.19];
ybs[8613]=['ε PsA',5.941911,-0.4701247,4.17];
ybs[8614]=['',5.9422355,-0.0601562,6.31];
ybs[8615]=['β Oct',5.9692109,-1.4184861,4.15];
ybs[8616]=['',5.9423618,0.2558116,5.71];
ybs[8617]=['11 Lac',5.9402932,0.7746435,4.46];
ybs[8618]=['',5.9391158,0.9416659,5.93];
ybs[8619]=['ζ Peg',5.9449537,0.190921,3.4];
ybs[8620]=['',5.9508201,-0.8220984,5.98];
ybs[8621]=['β Gru',5.9510444,-0.8164115,2.1];
ybs[8622]=['19 PsA',5.9494067,-0.5105628,6.17];
ybs[8623]=['',5.9449509,0.5423334,6.34];
ybs[8624]=['',5.9512008,-0.7703879,6.07];
ybs[8625]=['12 Lac',5.9445831,0.7039459,5.25];
ybs[8626]=['ο Peg',5.9459959,0.5133908,4.79];
ybs[8627]=['',5.9470696,0.2552378,5.9];
ybs[8628]=['',5.9451025,0.7270524,5.94];
ybs[8629]=['ρ Gru',5.9545458,-0.7209347,4.85];
ybs[8630]=['',5.952143,-0.1431836,6.45];
ybs[8631]=['',5.9584873,-1.0540285,6.3];
ybs[8632]=['67 Aqr',5.9529126,-0.1196406,6.41];
ybs[8633]=['',5.9480379,0.9427675,6.12];
ybs[8634]=['66 Aqr',5.954581,-0.3267665,4.69];
ybs[8635]=['η Peg',5.9514201,0.5293446,2.94];
ybs[8636]=['',5.9509587,0.6616646,6.43];
ybs[8637]=['',5.951416,0.8251296,6.39];
ybs[8638]=['',5.9547735,0.1928085,6.51];
ybs[8639]=['',5.9560041,0.6906889,5.95];
ybs[8640]=['η Gru',5.9641263,-0.9318663,4.85];
ybs[8641]=['13 Lac',5.9559799,0.7317671,5.08];
ybs[8642]=['',5.9641438,-0.8105173,5.51];
ybs[8643]=['',5.9661838,-0.8529519,6.62];
ybs[8644]=['',5.9676748,-0.8652894,6.48];
ybs[8645]=['45 Peg',5.9623483,0.339901,6.25];
ybs[8646]=['',5.9589012,0.9184853,6.55];
ybs[8647]=['',5.9687251,-0.8173553,6.56];
ybs[8648]=['ξ Oct',5.9873625,-1.3965274,5.35];
ybs[8649]=['',5.9835238,-1.3428849,6.73];
ybs[8650]=['ξ Peg',5.9677742,0.2143474,4.19];
ybs[8651]=['',5.9650074,0.779367,5.76];
ybs[8652]=['λ Peg',5.9669309,0.4131883,3.95];
ybs[8653]=['',5.9710536,-0.5943298,6.28];
ybs[8654]=['',5.9762911,-1.0746947,6.37];
ybs[8655]=['68 Aqr',5.9718726,-0.3404225,5.26];
ybs[8656]=['',5.9731538,-0.6652034,6.71];
ybs[8657]=['',5.9809113,-1.2259005,6.34];
ybs[8658]=['τ1 Aqr',5.9725126,-0.2434352,5.66];
ybs[8659]=['',5.9736305,-0.4503529,6.3];
ybs[8660]=['ε Gru',5.9767874,-0.8937521,3.49];
ybs[8661]=['70 Aqr',5.9759191,-0.182332,6.19];
ybs[8662]=['',5.9699257,1.0226109,6.36];
ybs[8663]=['',5.9739355,0.6549402,5.9];
ybs[8664]=['τ2 Aqr',5.9806994,-0.2353338,4.01];
ybs[8665]=['',5.982654,-0.5706589,6.33];
ybs[8666]=['',5.9802114,0.184791,6.54];
ybs[8667]=['',5.9762388,0.9516186,6.12];
ybs[8668]=['',5.9756354,1.1003785,6.06];
ybs[8669]=['μ Peg',5.9820786,0.4312811,3.48];
ybs[8670]=['',5.9873364,-0.6815137,5.42];
ybs[8671]=['',5.9909618,-1.0432215,6.46];
ybs[8672]=['',5.9764808,1.1986752,6.19];
ybs[8673]=['',5.9804493,0.9775877,5.43];
ybs[8674]=['',5.9929199,-1.1009423,6.12];
ybs[8675]=['14 Lac',5.9833891,0.7341305,5.92];
ybs[8676]=['',5.9849753,0.3359733,6.4];
ybs[8677]=['',5.9823525,0.8863808,6.21];
ybs[8678]=['21 PsA',5.98855,-0.5135977,5.97];
ybs[8679]=['ι Cep',5.9795801,1.1573176,3.52];
ybs[8680]=['γ PsA',5.9937316,-0.571879,4.46];
ybs[8681]=['',5.9872712,1.0787191,5.6];
ybs[8682]=['σ Peg',5.9927161,0.17357,5.16];
ybs[8683]=['λ Aqr',5.9938229,-0.1303833,3.74];
ybs[8684]=['15 Lac',5.990663,0.7578519,4.94];
ybs[8685]=['τ1 Gru',5.9988215,-0.8462856,6.04];
ybs[8686]=['ρ Ind',6.0041805,-1.221102,6.05];
ybs[8687]=['',5.966147,1.453202,4.74];
ybs[8688]=['',5.9954213,0.2958416,5.64];
ybs[8689]=['74 Aqr',5.9976299,-0.2008392,5.8];
ybs[8690]=['',5.9941748,0.8817626,6.46];
ybs[8691]=['',5.995768,0.7029594,6.34];
ybs[8692]=['',5.9947097,1.0508709,6.01];
ybs[8693]=['',5.9977824,0.7829305,5.81];
ybs[8694]=['δ Aqr',6.0027835,-0.2742128,3.27];
ybs[8695]=['78 Aqr',6.0023401,-0.1238335,6.19];
ybs[8696]=['77 Aqr',6.0032606,-0.2820859,5.56];
ybs[8697]=['',5.9998216,0.7066219,5.81];
ybs[8698]=['',6.0056423,-0.6331868,6.4];
ybs[8699]=['',6.0022165,0.2976005,6.12];
ybs[8700]=['1 Psc',6.0041065,0.0204965,6.11];
ybs[8701]=['',6.0050012,-0.085139,5.72];
ybs[8702]=['ρ Peg',6.005064,0.1557795,4.9];
ybs[8703]=['',6.0039241,0.6490284,5.91];
ybs[8704]=['',6.008227,-0.5501852,6.1];
ybs[8705]=['δ PsA',6.0086385,-0.5660093,4.21];
ybs[8706]=['',6.0105941,-0.5490058,6.48];
ybs[8707]=['τ3 Gru',6.0125882,-0.8353018,5.7];
ybs[8708]=['',6.0069863,0.6363716,5.74];
ybs[8709]=['',6.0121489,0.2087104,6.51];
ybs[8710]=['16 Lac',6.0097513,0.7280416,5.59];
ybs[8711]=['',6.0097671,0.8699322,4.95];
ybs[8712]=['',6.0141756,-0.0820312,6.31];
ybs[8713]=['α PsA',6.0160248,-0.5150853,1.16];
ybs[8714]=['51 Peg',6.0147068,0.3644049,5.49];
ybs[8715]=['',6.0152297,0.0684215,6.28];
ybs[8716]=['',6.0125964,0.8516174,5.43];
ybs[8717]=['',6.0201613,-0.6180722,6.13];
ybs[8718]=['',6.0154045,0.6879893,6.18];
ybs[8719]=['',6.0183924,-0.0398842,6.16];
ybs[8720]=['',6.0189796,-0.0226924,6.37];
ybs[8721]=['',5.9794436,1.4870912,5.9];
ybs[8722]=['',6.0197095,0.1652315,6.43];
ybs[8723]=['',6.0202739,0.1300246,6.33];
ybs[8724]=['52 Peg',6.0223571,0.2066311,5.75];
ybs[8725]=['',6.0245003,-0.5122884,5.51];
ybs[8726]=['',6.0243222,-0.2262047,6.07];
ybs[8727]=['2 Psc',6.0235905,0.0187276,5.43];
ybs[8728]=['',6.0266312,-0.4372721,5.65];
ybs[8729]=['',6.0216748,0.9209166,6.29];
ybs[8730]=['',6.0213671,1.0458869,6.43];
ybs[8731]=['',6.0280022,-0.4453435,6.29];
ybs[8732]=['ζ Gru',6.0304704,-0.9188066,4.12];
ybs[8733]=['',5.995847,1.474028,4.71];
ybs[8734]=['',6.0314956,-0.8873174,5.68];
ybs[8735]=['3 Psc',6.0287244,0.00517,6.21];
ybs[8736]=['',6.0290631,0.0544903,5.83];
ybs[8737]=['',6.0255541,0.9958078,5];
ybs[8738]=['',6.0287557,0.5444284,6.6];
ybs[8739]=['',6.0320193,-0.5016623,5.55];
ybs[8740]=['',6.0279576,0.7938695,6.5];
ybs[8741]=['',6.0322155,-0.3958468,6.28];
ybs[8742]=['81 Aqr',6.0321113,-0.1213113,6.21];
ybs[8743]=['',6.0295437,0.6775102,6.54];
ybs[8744]=['',6.0326795,-0.0803006,5.94];
ybs[8745]=['',6.0375298,-0.6337326,6.47];
ybs[8746]=['',6.0317928,0.9986084,6.2];
ybs[8747]=['ο And',6.0339032,0.7406593,3.62];
ybs[8748]=['82 Aqr',6.0371237,-0.1128101,6.15];
ybs[8749]=['',6.0381041,-0.3623286,5.97];
ybs[8750]=['',6.0367999,0.556606,6.57];
ybs[8751]=['2 And',6.0368832,0.7481948,5.1];
ybs[8752]=['π PsA',6.0415607,-0.6045593,5.11];
ybs[8753]=['',6.0375105,0.7709038,6.39];
ybs[8754]=['',6.0484266,-1.1992045,5.52];
ybs[8755]=['',6.0371765,0.9659879,6.5];
ybs[8756]=['',6.0438083,-0.7219996,5.79];
ybs[8757]=['',6.043267,-0.0817596,6.68];
ybs[8758]=['β Psc',6.0428531,0.0686052,4.53];
ybs[8759]=['κ Gru',6.0469499,-0.9399315,5.37];
ybs[8760]=['β Peg',6.0421854,0.4920702,2.42];
ybs[8761]=['',6.0434405,0.1174165,6.41];
ybs[8762]=['',6.0399294,1.0569015,6.74];
ybs[8763]=['',6.0398382,1.0240795,6.43];
ybs[8764]=['',6.0403017,1.1749538,5.24];
ybs[8765]=['3 And',6.0436536,0.8755102,4.65];
ybs[8766]=['α Peg',6.0466185,0.2673176,2.49];
ybs[8767]=['83 Aqr',6.0485623,-0.1323425,5.43];
ybs[8768]=['',6.0488609,-0.2961513,6.14];
ybs[8769]=['',6.0481054,0.291016,6.44];
ybs[8770]=['',6.0490557,0.024747,6.39];
ybs[8771]=['',6.0649466,-1.385259,6.12];
ybs[8772]=['θ Gru',6.0564191,-0.7576372,4.28];
ybs[8773]=['',6.0533199,0.3251298,6.13];
ybs[8774]=['86 Aqr',6.0553241,-0.4124551,4.47];
ybs[8775]=['υ Gru',6.0564127,-0.6768575,5.61];
ybs[8776]=['',6.0577326,-0.8638593,6.33];
ybs[8777]=['',6.0543044,0.3494486,6.3];
ybs[8778]=['',6.0581305,-0.8827038,5.83];
ybs[8779]=['',6.0649519,-1.2823814,6.15];
ybs[8780]=['55 Peg',6.0564601,0.1661657,4.52];
ybs[8781]=['56 Peg',6.0567926,0.4464464,4.76];
ybs[8782]=['1 Cas',6.0540602,1.0390088,4.85];
ybs[8783]=['',6.0582381,0.5748597,6.02];
ybs[8784]=['',6.058427,0.3708016,5.99];
ybs[8785]=['',6.0573619,0.8059796,6.66];
ybs[8786]=['',6.0566521,0.92376,6.11];
ybs[8787]=['',6.0626459,-0.5011195,5.6];
ybs[8788]=['',6.0564949,1.0443816,6.4];
ybs[8789]=['4 And',6.0589015,0.8115509,5.33];
ybs[8790]=['5 And',6.0592955,0.8623159,5.7];
ybs[8791]=['',6.061335,0.77969,6.56];
ybs[8792]=['5 Psc',6.0638345,0.0390799,5.4];
ybs[8793]=['',6.05908,1.1125525,6.26];
ybs[8794]=['',6.0714439,-1.1649372,6.47];
ybs[8795]=['',6.0817226,-1.4102442,6.41];
ybs[8796]=['',6.0597497,1.1228357,6.21];
ybs[8797]=['88 Aqr',6.0673558,-0.3675851,3.66];
ybs[8798]=['',6.0687119,-0.4882934,5.87];
ybs[8799]=['',6.0698024,-0.7461121,5.81];
ybs[8800]=['57 Peg',6.0674649,0.1533909,5.12];
ybs[8801]=['',6.0689576,-0.2513116,6.42];
ybs[8802]=['89 Aqr',6.0694015,-0.3900117,4.69];
ybs[8803]=['',6.0706796,-0.7065121,5.83];
ybs[8804]=['π Cep',6.0587152,1.3177014,4.41];
ybs[8805]=['ι Gru',6.0716004,-0.7877568,3.9];
ybs[8806]=['58 Peg',6.0696384,0.173371,5.39];
ybs[8807]=['2 Cas',6.0677212,1.0375022,5.7];
ybs[8808]=['',6.0732297,-0.5133613,6.51];
ybs[8809]=['',6.0725677,0.309028,5.71];
ybs[8810]=['6 And',6.0711786,0.7619356,5.94];
ybs[8811]=['59 Peg',6.0771176,0.1541416,5.16];
ybs[8812]=['60 Peg',6.0773361,0.4705214,6.17];
ybs[8813]=['',6.0842522,-0.8640614,6.8];
ybs[8814]=['',6.0882925,-1.0923682,6.12];
ybs[8815]=['7 And',6.0802559,0.8642544,4.52];
ybs[8816]=['',6.0827577,0.5158052,6.35];
ybs[8817]=['',6.0833096,0.9997271,5.56];
ybs[8818]=['',6.0845426,0.1950726,5.82];
ybs[8819]=['φ Aqr',6.0885087,-0.1036196,4.22];
ybs[8820]=['',6.0916528,-0.7154727,5.77];
ybs[8821]=['',6.090053,-0.1845975,6.12];
ybs[8822]=['',6.0876195,0.8854001,6.31];
ybs[8823]=['',6.0884115,0.521567,6.41];
ybs[8824]=['',6.0895374,0.4226316,6.36];
ybs[8825]=['',6.0939373,-0.059068,5.55];
ybs[8826]=['ψ1 Aqr',6.0953707,-0.1566556,4.21];
ybs[8827]=['61 Peg',6.0945838,0.4949725,6.49];
ybs[8828]=['',6.1006677,-1.0801656,5.66];
ybs[8829]=['',6.0883878,1.2975309,5.84];
ybs[8830]=['',6.0954569,0.4342935,6.6];
ybs[8831]=['',6.099033,-0.7745251,5.92];
ybs[8832]=['',6.0997267,-0.7170211,6.47];
ybs[8833]=['γ Tuc',6.1026077,-1.0144484,3.99];
ybs[8834]=['',6.1113196,-1.3851001,6.33];
ybs[8835]=['χ Aqr',6.0995346,-0.1328982,5.06];
ybs[8836]=['',6.0930695,1.2391854,5.56];
ybs[8837]=['γ Psc',6.1008455,0.0592437,3.69];
ybs[8838]=['',6.0983637,0.93071,5.54];
ybs[8839]=['',6.0970314,1.0834161,6.53];
ybs[8840]=['',6.1068177,-1.1756329,6.13];
ybs[8841]=['',6.1031303,-0.2024725,6.34];
ybs[8842]=['',6.1009862,0.7902216,6.43];
ybs[8843]=['ψ2 Aqr',6.1041461,-0.1583056,4.39];
ybs[8844]=['φ Gru',6.1055378,-0.7105613,5.53];
ybs[8845]=['8 And',6.1029811,0.8574369,4.85];
ybs[8846]=['',6.1038602,0.7958901,6.48];
ybs[8847]=['τ Oct',6.1543359,-1.5234067,5.49];
ybs[8848]=['γ Scl',6.1083262,-0.5658288,4.41];
ybs[8849]=['9 And',6.1058937,0.7310471,6.02];
ybs[8850]=['ψ3 Aqr',6.1087652,-0.1657798,4.98];
ybs[8851]=['94 Aqr',6.1094439,-0.2329408,5.08];
ybs[8852]=['',6.1001217,1.3161765,6.38];
ybs[8853]=['96 Aqr',6.1106496,-0.0874768,5.55];
ybs[8854]=['',6.1107388,-0.3135115,5.93];
ybs[8855]=['',6.1086957,0.789754,6.5];
ybs[8856]=['',6.1122439,-0.5863545,6.37];
ybs[8857]=['ο Cep',6.1063535,1.190733,4.75];
ybs[8858]=['',6.1106249,0.6092198,6.32];
ybs[8859]=['11 And',6.1106469,0.8506328,5.44];
ybs[8860]=['',6.1115103,0.8463667,6.32];
ybs[8861]=['10 And',6.112374,0.7363628,5.79];
ybs[8862]=['',6.1172744,-0.8760532,6.05];
ybs[8863]=['7 Psc',6.1147027,0.095886,5.05];
ybs[8864]=['',6.1162455,-0.1011515,6.17];
ybs[8865]=['τ Peg',6.1158666,0.4163094,4.6];
ybs[8866]=['',6.1136366,1.0835432,6.45];
ybs[8867]=['63 Peg',6.116647,0.5328055,5.59];
ybs[8868]=['',6.1188952,-0.4690419,5.64];
ybs[8869]=['',6.1161108,0.7719397,6.13];
ybs[8870]=['12 And',6.1168497,0.6683692,5.77];
ybs[8871]=['',6.1151118,1.0877859,6.39];
ybs[8872]=['64 Peg',6.1213904,0.5571981,5.32];
ybs[8873]=['',6.1216695,0.466378,6.62];
ybs[8874]=['',6.1266479,-1.0462052,6.09];
ybs[8875]=['97 Aqr',6.1248994,-0.2605167,5.2];
ybs[8876]=['65 Peg',6.1247872,0.3654941,6.29];
ybs[8877]=['98 Aqr',6.1263112,-0.3488542,3.97];
ybs[8878]=['66 Peg',6.1265922,0.2168848,5.08];
ybs[8879]=['',6.1237681,1.0514955,5.56];
ybs[8880]=['',6.130695,-0.9371645,6.15];
ybs[8881]=['',6.1299171,-0.7506957,6.1];
ybs[8882]=['',6.1286448,0.0070532,6.31];
ybs[8883]=['',6.1320306,-0.9037071,5.75];
ybs[8884]=['',6.1295835,0.5697476,6.69];
ybs[8885]=['',6.1313605,-0.3241901,6.19];
ybs[8886]=['',6.1369154,-0.9902352,5.59];
ybs[8887]=['',6.1329709,0.7195222,6.72];
ybs[8888]=['67 Peg',6.1341937,0.5671941,5.57];
ybs[8889]=['4 Cas',6.1337748,1.0890086,4.98];
ybs[8890]=['υ Peg',6.1365834,0.4104497,4.4];
ybs[8891]=['99 Aqr',6.139729,-0.3582991,4.39];
ybs[8892]=['ο Gru',6.1424511,-0.9181952,5.52];
ybs[8893]=['',6.1449475,-1.1600874,6.45];
ybs[8894]=['',6.1453243,-1.0186284,5.63];
ybs[8895]=['',6.1447763,-0.8734365,6.2];
ybs[8896]=['κ Psc',6.1434828,0.0238855,4.94];
ybs[8897]=['9 Psc',6.1448506,0.0215636,6.25];
ybs[8898]=['13 And',6.144057,0.7509268,5.75];
ybs[8899]=['',6.1483886,-0.6183944,6.32];
ybs[8900]=['69 Peg',6.1465882,0.4412237,5.98];
ybs[8901]=['θ Psc',6.1479744,0.1133057,4.28];
ybs[8902]=['',6.1485763,-0.1978621,6.37];
ybs[8903]=['',6.144201,1.229981,5.6];
ybs[8904]=['',6.1530977,-1.0995174,5.68];
ybs[8905]=['',6.1528236,-0.7746583,6.43];
ybs[8906]=['',6.1525931,-0.1597497,6.18];
ybs[8907]=['',6.1528016,0.404234,6.35];
ybs[8908]=['70 Peg',6.1531232,0.2246882,4.55];
ybs[8909]=['',6.1548615,-0.0771369,6.25];
ybs[8910]=['',6.1571133,0.8595092,6.17];
ybs[8911]=['',6.1565967,1.0238464,4.91];
ybs[8912]=['',6.1595567,0.6767544,6.05];
ybs[8913]=['',6.1613403,-0.1077755,6.39];
ybs[8914]=['',6.1634378,-0.7806916,6.02];
ybs[8915]=['14 And',6.1622986,0.6867811,5.22];
ybs[8916]=['',6.1635412,-0.0693584,6.49];
ybs[8917]=['100 Aqr',6.1643888,-0.3709899,6.29];
ybs[8918]=['',6.1642372,0.4977138,6.41];
ybs[8919]=['13 Psc',6.1654264,-0.0169738,6.38];
ybs[8920]=['',6.1723824,-1.3486488,5.81];
ybs[8921]=['',6.1672211,0.6120142,6.65];
ybs[8922]=['β Scl',6.1700143,-0.6580758,4.37];
ybs[8923]=['',6.137597,1.52409,5.58];
ybs[8924]=['101 Aqr',6.1712538,-0.363047,4.71];
ybs[8925]=['71 Peg',6.1719098,0.3946588,5.32];
ybs[8926]=['',6.1728407,0.7883907,6.24];
ybs[8927]=['',6.1739115,0.3657208,6.06];
ybs[8928]=['72 Peg',6.1739846,0.5487088,4.98];
ybs[8929]=['14 Psc',6.1749824,-0.0197932,5.87];
ybs[8930]=['',6.1800533,-1.127063,7.4];
ybs[8931]=['',6.1779739,-0.2641096,5.96];
ybs[8932]=['15 And',6.1768675,0.7042377,5.59];
ybs[8933]=['73 Peg',6.1769583,0.586617,5.63];
ybs[8934]=['ι Phe',6.1792154,-0.7417914,4.71];
ybs[8935]=['',6.1775514,0.6656225,6.18];
ybs[8936]=['',6.181048,-0.128298,6.39];
ybs[8937]=['',6.1779696,1.2523732,5.84];
ybs[8938]=['',6.1826582,0.4306537,6.45];
ybs[8939]=['16 Psc',6.1847369,0.0386725,5.68];
ybs[8940]=['',6.1851387,0.5762679,6.35];
ybs[8941]=['',6.1879265,-0.5542686,6.52];
ybs[8942]=['',6.1943065,-1.339651,6];
ybs[8943]=['',6.1903368,-0.2259619,5.65];
ybs[8944]=['',6.1913125,-0.7920108,4.74];
ybs[8945]=['74 Peg',6.1902502,0.2956442,6.26];
ybs[8946]=['λ And',6.189678,0.8128288,3.82];
ybs[8947]=['',6.1895532,0.777418,5.8];
ybs[8948]=['75 Peg',6.1914816,0.3231334,5.53];
ybs[8949]=['',6.1914797,0.8083205,6.58];
ybs[8950]=['ι And',6.1921981,0.7571533,4.29];
ybs[8951]=['θ Phe',6.1983595,-0.8119984,6.09];
ybs[8952]=['18 And',6.1965289,0.8828809,5.3];
ybs[8953]=['ω1 Aqr',6.1996083,-0.2462303,5];
ybs[8954]=['ι Psc',6.2002683,0.1001837,4.13];
ybs[8955]=['',6.2001178,0.1708841,5.97];
ybs[8956]=['',6.1962196,1.316091,5.95];
ybs[8957]=['',6.1970607,1.2935763,5.98];
ybs[8958]=['',6.2005764,0.6591449,6.53];
ybs[8959]=['γ Cep',6.1968483,1.3569269,3.21];
ybs[8960]=['μ Scl',6.2033834,-0.5577952,5.31];
ybs[8961]=['κ And',6.2021227,0.7757573,4.14];
ybs[8962]=['',6.2033338,0.6428846,6.23];
ybs[8963]=['',6.2054469,-0.4196908,6.6];
ybs[8964]=['',6.2055466,-0.2018786,5.89];
ybs[8965]=['103 Aqr',6.2074259,-0.3126485,5.34];
ybs[8966]=['',6.2066389,0.8661371,6.26];
ybs[8967]=['104 Aqr',6.2082465,-0.3089686,4.82];
ybs[8968]=['',6.2089704,0.1285322,5.89];
ybs[8969]=['λ Psc',6.2094287,0.0330531,4.5];
ybs[8970]=['',6.2085918,1.0013616,6.24];
ybs[8971]=['',6.2101562,0.7872439,6.57];
ybs[8972]=['',6.2112995,-0.2676281,5.28];
ybs[8973]=['ω2 Aqr',6.2124164,-0.2518715,4.49];
ybs[8974]=['',6.2104291,1.1279952,6.56];
ybs[8975]=['',6.2112464,1.0784959,6.4];
ybs[8976]=['77 Peg',6.2151966,0.1823038,5.06];
ybs[8977]=['',6.21723,-0.2647766,6.36];
ybs[8978]=['',6.2181925,-0.7848652,6.09];
ybs[8979]=['',6.2201705,-1.2282997,6.07];
ybs[8980]=['',6.2215583,-1.3731814,5.75];
ybs[8981]=['',6.2191049,-1.1220821,5.72];
ybs[8982]=['78 Peg',6.2178444,0.5144452,4.93];
ybs[8983]=['106 Aqr',6.2188802,-0.3170054,5.24];
ybs[8984]=['',6.2201223,-0.4560982,6.17];
ybs[8985]=['',6.2212809,0.9758767,6.51];
ybs[8986]=['',6.2268689,-0.6993284,6.31];
ybs[8987]=['107 Aqr',6.226788,-0.3240051,5.29];
ybs[8988]=['ψ And',6.2267129,0.8121752,4.95];
ybs[8989]=['19 Psc',6.2283852,0.0628425,5.04];
ybs[8990]=['',6.2290903,1.1675585,5.95];
ybs[8991]=['σ Phe',6.2323371,-0.8746316,5.18];
ybs[8992]=['',6.2330046,-1.1917142,6.89];
ybs[8993]=['τ Cas',6.2311233,1.0256585,4.87];
ybs[8994]=['',6.2322251,-0.2058942,5.73];
ybs[8995]=['',6.2310147,1.0047049,5.51];
ybs[8996]=['',6.2333502,0.8193705,6.07];
ybs[8997]=['20 Psc',6.2351601,-0.0462108,5.49];
ybs[8998]=['',6.2347842,1.1854438,5.04];
ybs[8999]=['',6.2377839,-0.109372,6.07];
ybs[9000]=['',6.2389914,0.0406343,6.46];
ybs[9001]=['δ Scl',6.2395021,-0.4889762,4.57];
ybs[9002]=['',6.2380445,1.1342963,6.41];
ybs[9003]=['6 Cas',6.2388844,1.0878367,5.43];
ybs[9004]=['',6.2391699,1.0488189,6.34];
ybs[9005]=['',6.2404969,1.0310894,6.33];
ybs[9006]=['',6.2420957,-0.2748385,6.24];
ybs[9007]=['21 Psc',6.2417713,0.0207717,5.77];
ybs[9008]=['',6.2431879,-1.0947651,6.59];
ybs[9009]=['',6.2426919,0.6377311,5.9];
ybs[9010]=['79 Peg',6.2425913,0.5053867,5.97];
ybs[9011]=['',6.2434198,-0.4401259,6.42];
ybs[9012]=['',6.2452196,-0.1720917,5.94];
ybs[9013]=['',6.245658,0.9029584,6.44];
ybs[9014]=['',6.246578,-0.2493709,5.72];
ybs[9015]=['80 Peg',6.2500299,0.164539,5.79];
ybs[9016]=['108 Aqr',6.2500729,-0.3280317,5.18];
ybs[9017]=['γ1 Oct',6.2537788,-1.4295088,5.11];
ybs[9018]=['22 Psc',6.252707,0.0531339,5.55];
ybs[9019]=['',6.2523889,1.3563567,6.55];
ybs[9020]=['',6.2545403,0.3802184,6.11];
ybs[9021]=['φ Peg',6.2549728,0.3357028,5.08];
ybs[9022]=['',6.2550602,-0.2467378,5.87];
ybs[9023]=['',6.2544487,1.3204951,6.39];
ybs[9024]=['82 Peg',6.2555496,0.1930609,5.3];
ybs[9025]=['',6.2565447,-0.1550304,5.75];
ybs[9026]=['24 Psc',6.2569096,-0.0530837,5.93];
ybs[9027]=['25 Psc',6.2575735,0.0384783,6.28];
ybs[9028]=['',6.258763,-0.4208875,6.24];
ybs[9029]=['',6.2631687,-0.4699844,6.35];
ybs[9030]=['ρ Cas',6.2631927,1.005546,4.54];
ybs[9031]=['',6.2644323,-0.7013762,6.03];
ybs[9032]=['',6.2649785,0.0038968,5.61];
ybs[9033]=['26 Psc',6.2665156,0.1254058,6.21];
ybs[9034]=['',6.267183,-0.5551466,6.1];
ybs[9035]=['',6.267183,-0.5544921,6.83];
ybs[9036]=['',6.2676084,0.4549918,6.54];
ybs[9037]=['',6.2683504,1.0040239,6];
ybs[9038]=['',6.2683596,0.8285068,6];
ybs[9039]=['',6.2725027,-0.4297542,6.31];
ybs[9040]=['',6.2733243,0.3972749,6.15];
ybs[9041]=['',6.2721078,1.4539506,6.59];
ybs[9042]=['',6.2749221,0.7465202,5.97];
ybs[9043]=['',6.2752935,-0.4626779,6.26];
ybs[9044]=['',6.2752681,0.974242,5.55];
ybs[9045]=['',6.2761633,-1.0968045,5.97];
ybs[9046]=['γ2 Oct',6.2771779,-1.4321451,5.73];
ybs[9047]=['η Tuc',6.2772736,-1.1202257,5];
ybs[9048]=['',6.2770869,1.0496015,6.47];
ybs[9049]=['ψ Peg',6.2779816,0.4407919,4.66];
ybs[9050]=['1 Cet',6.2805876,-0.2745992,6.26];
ybs[9051]=['',6.2808351,0.8988924,4.8];
ybs[9052]=['27 Psc',6.2819815,-0.0600739,4.86];
ybs[9053]=['',6.2826174,0.5671586,6.52];
ybs[9054]=['π Phe',6.2831072,-0.9185966,5.13];
ybs[9055]=['',6.2824182,0.8120525,6.54];
ybs[9056]=['σ Cas',0.0002515,0.9751002,4.88];
ybs[9057]=['ω Psc',0.0015809,0.1217797,4.01];
ybs[9058]=['',0.0022501,-0.5126184,5.62];
ybs[9059]=['',0.0023447,0.5905945,6.58];
ybs[9060]=['',0.0023447,0.5905945,6.58];
ybs[9061]=['ε Tuc',0.0042124,-1.1425465,4.5];
ybs[9062]=['',0.0059736,-0.7710241,6.29];
ybs[9063]=['',0.0063266,0.4718054,6.46];
ybs[9064]=['',0.006847,1.0415052,6.19];
ybs[9065]=['',0.0077735,0.7918116,6.38];
ybs[9066]=['τ Phe',0.0092592,-0.8499033,5.71];
ybs[9067]=['',0.0103828,-0.8765584,5.53];
ybs[9068]=['',0.0103705,0.8743365,6.22];
ybs[9069]=['θ Oct',0.0114644,-1.3430607,4.78];
ybs[9070]=['',0.0116653,1.0705357,5.55];
ybs[9071]=['',0.0121509,0.7414393,6.25];
ybs[9072]=['29 Psc',0.0125392,-0.050848,5.1];
ybs[9073]=['85 Peg',0.0140649,0.4746609,5.75];
ybs[9074]=['30 Psc',0.0131343,-0.1029752,4.41];
ybs[9075]=['',0.0138359,-0.2541547,7.1];
ybs[9076]=['ζ Scl',0.0147441,-0.5167249,5.01];
ybs[9077]=['31 Psc',0.015075,0.15832,6.32];
ybs[9078]=['32 Psc',0.0154749,0.1500927,5.63];
ybs[9079]=['',0.015998,1.1556349,5.86];
ybs[9080]=['',0.0174888,-0.3478789,6.25];
ybs[9081]=['',0.0182201,-0.4194229,6.44];
ybs[9082]=['',0.019613,1.1127531,6.24];
ybs[9083]=['2 Cet',0.0208917,-0.3005806,4.55];
ybs[9084]=['',0.0215377,1.1663395,6.29];
ybs[9085]=['9 Cas',0.0231057,1.0891183,5.88];
ybs[9086]=['',0.0234579,-0.286492,5.78];
ybs[9087]=['',0.0234904,-0.5088421,6.4];
ybs[9088]=['3 Cet',0.0242186,-0.1814329,4.94];
ybs[9089]=['',0.0252005,1.1742709,5.67];
ybs[9090]=['',0.0247471,0.7366393,6.01];
ybs[9091]=['',0.024129,-1.2703147,7.31];
ybs[9092]=['',0.0259828,0.6069177,6.12];
ybs[9093]=['',0.0249071,-1.2448185,5.59];
ybs[9094]=['',0.026134,0.4671023,6.25];
ybs[9095]=['',0.0269417,1.0721254,5.8];    
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