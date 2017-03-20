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
 mean-equinox     j2000          j2017.5 (or similar)
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
  
  //all the aliases for a given moment in time  
  when {
    .. fill this out later
  }
  //the location of an observer on the Earth's surface (height is neglected - low precision)
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
    result.next = function(secs){
      return this.delta(60*60*24);
    };
    result.prev = function(secs){
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
  
  /* Example: 'UT 2016-01-31 02:56:03.123', plus truncations. */
  var when_parse = function(text){
    var original_text = text;
    var style = text.substring(0,2); //first 2 letters
    var text = text.substring(2).trim(); //chop off the 'style'
    var space = text.indexOf(" "); //between date and time
    var dot = text.indexOf("."); //either a fractional day, or a fractional second (but not both)
    var is_fractional_day = (dot === 10); 
    
    var date = (space === -1 ? text : text.substring(0, space));
    var parts = date.split('-');
    var year = parseInt(parts[0], 10);
    var month = parseInt(parts[1], 10);
    var day = is_fractional_day ? Math.floor(parseFloat(parts[2])) : parseInt(parts[2], 10);

    var hour = 0, minute = 0, seconds = 0, msecs = 0; //integers all
    var frac = 0, hour_dec = 0, minute_dec = 0, seconds_dec = 0; // decimal numbers
    
    if (! is_fractional_day && text.length > 10){
      //date and time string both present
      var time = text.substring(space+1);
      parts = time.split(':');
      hour = parseInt(parts[0], 10);
      minute = parts.length > 1 ? parseInt(parts[1],10) : 0;
      if (parts.length > 2){
        seconds_dec = parseFloat(parts[2]);
        if ('TT' === style){
          seconds_dec = seconds_dec - delta_t(year);
        }
        seconds = Math.floor(seconds_dec);
        msecs = Math.round((seconds_dec - seconds)*1000); 
      }
    }
    if (is_fractional_day) {
      frac = parseFloat(parts[2]) - day;
      hour_dec = frac * 24;
      hour = Math.floor(hour_dec);
      minute_dec = (hour_dec - hour) * 60;
      minute = Math.floor(minute_dec);
      seconds_dec = (minute_dec - minute) * 60;
      seconds = Math.floor(seconds_dec);
      msecs = Math.round((seconds_dec - seconds) * 1000);
    }
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
      ephem.mag = -12.7; //RASC Observer's Handbook
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
  

  var planet_orbit_start = when('UT 2017-02-16'); //Observer's Handbook p 23
  var planet_orbit_end = when('UT 2017-10-14');
  
  var mercury = {
    name: 'Mercury',
    symbol: '☿',
    orbit: function(when){
      var osc_a = build_osculating_orbit(when_j2000, planet_orbit_start, 0.387099, 0.205627, 7.0040, 48.3092, 77.4799, 291.8769);
      var osc_b = build_osculating_orbit(when_j2000, planet_orbit_end,   0.387098, 0.205638, 7.0040, 48.3084, 77.4839, (360*3)+194.0359); //4.0914 deg per d; 982 deg per 240d
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
      var osc_a = build_osculating_orbit(when_j2000, planet_orbit_start, 0.723331, 0.006746, 3.3944, 76.6331, 131.6533,        124.1110);
      var osc_b = build_osculating_orbit(when_j2000, planet_orbit_end,   0.723328, 0.006786, 3.3945, 76.6285, 131.3854, (360*1)+148.6197); //1.60212 deg per d; 384 deg per 240d
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
      var osc_a = build_osculating_orbit(when_j2000, planet_orbit_start, 1.523692, 0.093504, 1.8484, 49.5071, 336.1677,  33.5343);
      var osc_b = build_osculating_orbit(when_j2000, planet_orbit_end,   1.523644, 0.093489, 1.8484, 49.5068, 336.1245, 159.3104); //0.52402 deg per d; 126d in 240d 
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
      var osc_a = build_osculating_orbit(when_j2000, planet_orbit_start, 5.202158, 0.048896, 1.3037, 100.5130, 14.2689, 194.1602);
      var osc_b = build_osculating_orbit(when_j2000, planet_orbit_end,   5.202169, 0.048893, 1.3037, 100.5133, 14.2429, 214.1084); //19.94d per 240d 
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
      var osc_a = build_osculating_orbit(when_j2000, planet_orbit_start, 9.559738, 0.052963, 2.4873, 113.5851, 93.8056, 259.2040);
      var osc_b = build_osculating_orbit(when_j2000, planet_orbit_end,   9.564689, 0.052387, 2.4869, 113.5898, 93.6714, 267.2352); //8.0d per 240d 
      var result = current_osculating_orbit(osc_a, osc_b, when);
      return result;
    },
    add_physical: function(ephem){
      var i = degs(ephem.phase);
      var FUDGE_FACTOR_FOR_MISSING_RINGS =  -1.1; //gives the correct answer near opposition, 2017-06-15
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
      var osc_a = build_osculating_orbit(when_j2000, planet_orbit_start, 19.108196, 0.050844, 0.7719, 73.9838, 173.5173, 26.6258);
      var osc_b = build_osculating_orbit(when_j2000, planet_orbit_end,   19.108903, 0.050454, 0.7716, 74.0114, 174.0900, 29.3841);  
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
      var osc_a = build_osculating_orbit(when_j2000, planet_orbit_start, 29.963210, 0.006354, 1.7723, 131.8204, 68.0064, 342.1299);
      var osc_b = build_osculating_orbit(when_j2000, planet_orbit_end,   29.996402, 0.006013, 1.7721, 131.8154, 57.0098, 343.5104);  
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
      epoch: when("UT 2017-02-16"),
      a: 2.76788082,
      e: 0.0756827676,
      i: rads(10.5924016),
      Ω: rads(80.30985818),
      ω: rads(72.9077893),
      M0: rads(266.8015867),
      n: rads(0.2140341110) 
    }
  );
  
  var vesta = build_minor_planet('Vesta', 4, 3.20, 0.32, {
      equinox: when_j2000, 
      epoch: when("UT 2017-02-16"),
      a: 2.36130084,
      e: 0.089136053028,
      i: rads(7.14051581),
      Ω: rads(103.8420858),
      ω: rads(151.07635994),
      M0: rads(238.2473165),
      n: rads(0.2716295968) 
    }
  );
  
  var eunomia = build_minor_planet('Eunomia', 15, 5.28, 0.23, {
      equinox: when_j2000, 
      epoch: when("UT 2017-02-16"),
      a: 2.64373001,
      e: 0.1870383446,
      i: rads(11.73678248),
      Ω: rads(293.17644253),
      ω: rads(97.60123552),
      M0: rads(100.93449618),
      n: rads(0.229286466609) 
    }
  );
  
  var euterpe = build_minor_planet('Euterpe', 27, 7.0, 0.15, {
      equinox: when_j2000, 
      epoch: when("UT 2017-02-16"),
      a: 2.34703291809,
      e: 0.172794711675,
      i: rads(1.583681009491),
      Ω: rads(94.7897158272),
      ω: rads(356.616794216),
      M0: rads(115.756247508),
      n: rads(0.2741102670240) 
    }
  );
  
  var melpomene = build_minor_planet('Melpomene', 18, 6.51, 0.25, {
      equinox: when_j2000, 
      epoch: when("UT 2017-02-16"),
      a: 2.2953073624,
      e: 0.218895488291,
      i: rads(10.1334153580),
      Ω: rads(150.4683491),
      ω: rads(227.839623884),
      M0: rads(40.6383410088),
      n: rads(0.283428030723) 
    }
  );
  
  var astraea = build_minor_planet('Astraea', 5, 6.85, 0.15, {
      equinox: when_j2000, 
      epoch: when("UT 2017-02-16"),
      a: 2.57344954260,
      e: 0.19151497029024,
      i: rads(5.367899876),
      Ω: rads(141.5843547),
      ω: rads(358.7853053),
      M0: rads(91.269452015),
      n: rads(0.23874296747027) 
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
  var enke = build_comet('Enke', '2P/Encke', 7.0, 'bright', 'evening', {
      equinox: when_j2000,
      epoch: when("TT 2016-11-28"),  
      a: 2.21476565,
      e: 0.848333491,
      i: rads(11.7783682),
      Ω: rads(334.561109541),
      ω: rads(186.561039123),
      T: when("TT 2017-03-10.09193979") 
    }
  );
  var johnson = build_comet('Johnson', 'C/2015 V2 (Johnson)', 9.5, 'bright', 'morning', {
      equinox: when_j2000,
      epoch: when("TT 2016-02-10"),
      q: 1.637142990259, 
      e: 1.00145173422522,
      i: rads(49.8749930992),
      Ω: rads(69.855284922695),
      ω: rads(164.900083166),
      T: when("TT 2017-06-12.40859890") 
    }
  );
  var honda_mrkos_pajdusakova = build_comet('Honda Mrkos Pajdusakova', '45P/Honda-Mrkos-Pajdusakova', 9, 'fade', 'morning', {
      equinox: when_j2000,
      epoch: when("TT 2017-01-07"),
      a: 3.02462079, 
      e: 0.82392643199,
      i: rads(4.24948438),
      Ω: rads(89.005064665),
      ω: rads(326.2641887),
      T: when("TT 2016-12-31.267332") 
    }
  ); 
  var panstarrs_2015_er61 = build_comet('Panstarrs 2015 er61', 'C/2015 ER61 (PANSTARRS)', 9, 'bright', 'early morning', {
      equinox: when_j2000,
      epoch: when("TT 2016-02-26"),
      a: 3032.57267298, 
      q: 1.050451621,
      i: rads(6.2000087),
      Ω: rads(237.36607246),
      ω: rads(65.7039694152),
      T: when("TT 2017-05-09.76856804") 
    }
  );
  var tuttle_giacobinni_kresak =  build_comet('Tuttle Giacobinni Kresak', '41P/Tuttle-Giacobinni-Kresak', 10, 'bright', 'evening', {
      equinox: when_j2000,
      epoch: when("TT 2017-02-16"),
      a: 3.083813412556, 
      e: 0.66111470623,
      i: rads(9.2292781),
      Ω: rads(141.0751349),
      ω: rads(62.146269810),
      T: when("TT 2017-04-12.74959214") 
    }
  );  
  var comets = {
    /*testing only enke_test: enke_test,*/
    enke: enke, 
    johnson:johnson, 
    honda_mrkos_pajdusakova: honda_mrkos_pajdusakova, 
    panstarrs_2015_er61: panstarrs_2015_er61,
    tuttle_giacobinni_kresak: tuttle_giacobinni_kresak  
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
     {when:'UT 2016-10-01 00',text:'New Moon 401579.717 km'},
     {when:'UT 2016-10-02 03',text:'Spica 5.78°S of Moon'},
     {when:'UT 2016-10-03 17',text:'Venus 5.03°S of Moon'},
     {when:'UT 2016-10-04 11',text:'Moon at apogee 406096.075 km'},
     {when:'UT 2016-10-06 02',text:'Antares 9.88°S of Moon'},
     {when:'UT 2016-10-06 08',text:'Saturn 3.81°S of Moon'},
     {when:'UT 2016-10-08 12',text:'Mars 7.00°S of Moon'},
     {when:'UT 2016-10-09 05',text:'First Quarter 394263.975 km'},
     {when:'UT 2016-10-11 04',text:'Jupiter 0.87°S of Mercury'},
     {when:'UT 2016-10-13 06',text:'Neptune 1.16°S of Moon'},
     {when:'UT 2016-10-15 11',text:'Uranus at opposition'},
     {when:'UT 2016-10-16 02',text:'Uranus 2.84°N of Moon'},
     {when:'UT 2016-10-16 04',text:'Full Moon 358473.468 km'},
     {when:'UT 2016-10-17 00',text:'Moon at perigee 357860.832 km'},
     {when:'UT 2016-10-19 07',text:'Aldebaran 0.33°S of Moon'},
     {when:'UT 2016-10-20 14',text:'Spica 3.56°S of Mercury'},
     {when:'UT 2016-10-22 11',text:'Pollux 10.58°N of Moon'},
     {when:'UT 2016-10-22 19',text:'Last Quarter 380507.911 km'},
     {when:'UT 2016-10-25 04',text:'Regulus 1.55°N of Moon'},
     {when:'UT 2016-10-26 04',text:'Antares 3.14°S of Venus'},
     {when:'UT 2016-10-27 16',text:'Mercury in superior conjunction 0.51° North'},
     {when:'UT 2016-10-28 10',text:'Jupiter 1.42°S of Moon'},
     {when:'UT 2016-10-29 10',text:'Spica 5.77°S of Moon'},
     {when:'UT 2016-10-30 08',text:'Saturn 3.03°N of Venus, Antares nearby'},
     {when:'UT 2016-10-30 18',text:'New Moon 406275.919 km'},
     {when:'UT 2016-10-30 19',text:'Mercury 4.32°S of Moon'},
     {when:'UT 2016-10-31 19',text:'Moon at apogee 406661.670 km'},
     {when:'UT 2016-11-02 08',text:'Antares 9.74°S of Moon'},
     {when:'UT 2016-11-02 19',text:'Saturn 3.72°S of Moon, group with Venus'},
     {when:'UT 2016-11-03 04',text:'Venus 6.84°S of Moon'},
    {when: 'UT 2016-11-06', text:'End of daylight savings time'},
     {when:'UT 2016-11-06 12',text:'Mars 5.30°S of Moon'},
     {when:'UT 2016-11-07 20',text:'First Quarter 386497.594 km'},
     {when:'UT 2016-11-09 15',text:'Neptune 1.00°S of Moon'},
     {when:'UT 2016-11-12 11',text:'Uranus 2.84°N of Moon'},
     {when:'UT 2016-11-14 11',text:'Moon at perigee 356508.988 km'},
     {when:'UT 2016-11-14 14',text:'Full Moon 356520.182 km (largest of yr)'},
     {when:'UT 2016-11-15 17',text:'Aldebaran 0.45°S of Moon'},
     {when:'UT 2016-11-18 19',text:'Pollux 10.35°N of Moon'},
     {when:'UT 2016-11-18 21',text:'Antares 2.82°S of Mercury'},
     {when:'UT 2016-11-21 09',text:'Last Quarter 388237.681 km'},
     {when:'UT 2016-11-21 11',text:'Regulus 1.30°N of Moon'},
     {when:'UT 2016-11-24 01',text:'Saturn 3.47°N of Mercury'},
     {when:'UT 2016-11-25 02',text:'Jupiter 1.94°S of Moon'},
     {when:'UT 2016-11-25 16',text:'Spica 5.90°S of Moon'},
     {when:'UT 2016-11-27 20',text:'Moon at apogee 406554.249 km'},
     {when:'UT 2016-11-29 12',text:'New Moon 405613.981 km'},
     {when:'UT 2016-11-29 14',text:'Antares 9.67°S of Moon'},
     {when:'UT 2016-11-30 08',text:'Saturn 3.63°S of Moon'},
     {when:'UT 2016-12-01 04',text:'Mercury 7.09°S of Moon'},
     {when:'UT 2016-12-03 13',text:'Venus 5.81°S of Moon'},
     {when:'UT 2016-12-05 11',text:'Mars 2.95°S of Moon'},
     {when:'UT 2016-12-06 22',text:'Neptune 0.70°S of Moon'},
     {when:'UT 2016-12-07 09',text:'First Quarter 379097.283 km'},
     {when:'UT 2016-12-09 20',text:'Uranus 3.01°N of Moon'},
     {when:'UT 2016-12-10 12',text:'Saturn 1.30°N of Sun'},
     {when:'UT 2016-12-11 05',text:'Mercury at greatest elongation 20.8° East'},
     {when:'UT 2016-12-12 23',text:'Moon at perigee 358460.648 km'},
     {when:'UT 2016-12-13 05',text:'Aldebaran 0.46°S of Moon'},
     {when:'UT 2016-12-14 00',text:'Full Moon 359447.432 km'},
     {when:'UT 2016-12-16 05',text:'Pollux 10.15°N of Moon'},
     {when:'UT 2016-12-18 19',text:'Regulus 1.02°N of Moon'},
     {when:'UT 2016-12-21 02',text:'Last Quarter 396158.485 km'},
     {when:'UT 2016-12-21 11',text:'Solstice'},
     {when:'UT 2016-12-22 17',text:'Jupiter 2.41°S of Moon'},
     {when:'UT 2016-12-22 22',text:'Spica 6.14°S of Moon'},
     {when:'UT 2016-12-25 06',text:'Moon at apogee 405870.383 km'},
     {when:'UT 2016-12-26 20',text:'Antares 9.73°S of Moon'},
     {when:'UT 2016-12-27 21',text:'Saturn 3.60°S of Moon'},
     {when:'UT 2016-12-28 19',text:'Mercury in inferior conjunction 2.44° North'},
     {when:'UT 2016-12-29 05',text:'Mercury 1.76°S of Moon'},
     {when:'UT 2016-12-29 07',text:'New Moon 399571.298 km'},
     {when:'UT 2017-01-01 07',text:'Neptune 0.02°N of Mars'},
     {when:'UT 2017-01-02 09',text:'Venus 1.90°S of Moon'},
     {when:'UT 2017-01-03 04',text:'Neptune 0.40°S of Moon'},
     {when:'UT 2017-01-03 07',text:'Mars 0.24°S of Moon'},
     {when:'UT 2017-01-04 14',text:'Earth at perihelion 0.983309436 AU'},
     {when:'UT 2017-01-05 20',text:'First Quarter 373554.488 km'},
     {when:'UT 2017-01-06 02',text:'Uranus 3.27°N of Moon'},
     {when:'UT 2017-01-09 15',text:'Aldebaran 0.36°S of Moon'},
     {when:'UT 2017-01-10 06',text:'Moon at perigee 363238.471 km'},
     {when:'UT 2017-01-12 12',text:'Full Moon 366879.241 km'},
     {when:'UT 2017-01-12 13',text:'Venus at greatest elongation 47.1° East'},
     {when:'UT 2017-01-12 16',text:'Pollux 10.08°N of Moon'},
     {when:'UT 2017-01-13 02',text:'Neptune 0.41°S of Venus'},
     {when:'UT 2017-01-15 05',text:'Regulus 0.84°N of Moon'},
     {when:'UT 2017-01-19 05',text:'Jupiter 2.69°S of Moon'},
     {when:'UT 2017-01-19 06',text:'Spica 6.37°S of Moon'},
     {when:'UT 2017-01-19 10',text:'Mercury at greatest elongation 24.1° West'},
     {when:'UT 2017-01-19 22',text:'Last Quarter 402139.122 km'},
     {when:'UT 2017-01-20 21',text:'Spica 3.68°S of Jupiter'},
     {when:'UT 2017-01-22 00',text:'Moon at apogee 404914.424 km'},
     {when:'UT 2017-01-23 04',text:'Antares 9.87°S of Moon'},
     {when:'UT 2017-01-24 10',text:'Saturn 3.62°S of Moon'},
     {when:'UT 2017-01-26 01',text:'Mercury 3.71°S of Moon'},
     {when:'UT 2017-01-28 00',text:'New Moon 389613.481 km'},
     {when:'UT 2017-01-30 11',text:'Neptune 0.20°S of Moon'},
     {when:'UT 2017-01-31 15',text:'Venus 4.06°N of Moon'},
     {when:'UT 2017-02-01 01',text:'Mars 2.33°N of Moon'},
     {when:'UT 2017-02-02 08',text:'Uranus 3.48°N of Moon'},
     {when:'UT 2017-02-04 04',text:'First Quarter 370622.831 km'},
     {when:'UT 2017-02-05 22',text:'Aldebaran 0.24°S of Moon'},
     {when:'UT 2017-02-06 14',text:'Moon at perigee 368815.901 km'},
     {when:'UT 2017-02-09 02',text:'Pollux 10.10°N of Moon'},
     {when:'UT 2017-02-11 01',text:'Penumbral ECLIPSE Full Moon 377419.733 km'},
     {when:'UT 2017-02-11 14',text:'Regulus 0.79°N of Moon'},
     {when:'UT 2017-02-15 14',text:'Spica 6.47°S of Moon'},
     {when:'UT 2017-02-15 15',text:'Jupiter 2.70°S of Moon'},
     {when:'UT 2017-02-18 20',text:'Last Quarter 404372.825 km'},
     {when:'UT 2017-02-18 21',text:'Moon at apogee 404376.064 km'},
     {when:'UT 2017-02-19 12',text:'Antares 9.95°S of Moon'},
     {when:'UT 2017-02-20 23',text:'Saturn 3.58°S of Moon'},
     {when:'UT 2017-02-23 16',text:'Spica 3.82°S of Jupiter'},
     {when:'UT 2017-02-26 02',text:'Mercury 2.49°S of Moon'},
     {when:'UT 2017-02-26 15',text:'Annular ECLIPSE New Moon 378195.691 km'},
     {when:'UT 2017-02-26 21',text:'Neptune 0.10°S of Moon'},
     {when:'UT 2017-02-27 08',text:'Uranus 0.62°S of Mars'},
     {when:'UT 2017-02-28 20',text:'Venus 10.26°N of Moon'},
     {when:'UT 2017-03-01 16',text:'Uranus 3.58°N of Moon'},
     {when:'UT 2017-03-01 19',text:'Mars 4.33°N of Moon'},
     {when:'UT 2017-03-02 03',text:'Neptune 0.85°S of Sun'},
     {when:'UT 2017-03-03 08',text:'Moon at perigee 369062.405 km'},
     {when:'UT 2017-03-04 06',text:'Neptune 1.13°N of Mercury'},
     {when:'UT 2017-03-05 03',text:'Aldebaran 0.23°S of Moon'},
     {when:'UT 2017-03-05 12',text:'First Quarter 370468.466 km'},
     {when:'UT 2017-03-07 00',text:'Mercury in superior conjunction 1.69° South'},
     {when:'UT 2017-03-08 08',text:'Pollux 10.11°N of Moon'},
     {when:'UT 2017-03-10 23',text:'Regulus 0.80°N of Moon'},
     {when:'UT 2017-03-12 15',text:'Full Moon 388858.467 km'},
     {when:'UT 2017-03-14 20',text:'Jupiter 2.46°S of Moon'},
     {when:'UT 2017-03-14 23',text:'Spica 6.44°S of Moon'},
     {when:'UT 2017-03-16 23',text:'Venus 9.55°N of Mercury'},
     {when:'UT 2017-03-18 17',text:'Moon at apogee 404649.509 km'},
     {when:'UT 2017-03-18 20',text:'Antares 9.89°S of Moon'},
     {when:'UT 2017-03-20 10',text:'Equinox'},
     {when:'UT 2017-03-20 16',text:'Last Quarter 402269.201 km'},
     {when:'UT 2017-03-25 10',text:'Venus in inferior conjunction 8.29° North'},
     {when:'UT 2017-03-26 08',text:'Neptune 0.00°N of Moon'},
     {when:'UT 2017-03-27 06',text:'Uranus 2.41°S of Mercury'},
     {when:'UT 2017-03-27 13',text:'Venus 11.32°N of Moon'},
     {when:'UT 2017-03-28 03',text:'New Moon 367868.370 km'},
     {when:'UT 2017-03-29 03',text:'Uranus 3.62°N of Moon'},
     {when:'UT 2017-03-29 07',text:'Mercury 6.59°N of Moon'},
     {when:'UT 2017-03-30 13',text:'Moon at perigee 363853.934 km'},
     {when:'UT 2017-04-01 09',text:'Aldebaran 0.34°S of Moon'},
     {when:'UT 2017-04-01 10',text:'Mercury at greatest elongation 19.0° East'},
     {when:'UT 2017-04-03 19',text:'First Quarter 373043.990 km'},
     {when:'UT 2017-04-04 14',text:'Pollux 9.99°N of Moon'},
     {when:'UT 2017-04-07 05',text:'Regulus 0.72°N of Moon'},
     {when:'UT 2017-04-07 22',text:'Jupiter at opposition'},
     {when:'UT 2017-04-10 21',text:'Jupiter 2.18°S of Moon'},
     {when:'UT 2017-04-11 06',text:'Full Moon 398715.449 km'},
     {when:'UT 2017-04-14 06',text:'Uranus 0.56°S of Sun'},
     {when:'UT 2017-04-15 04',text:'Antares 9.74°S of Moon'},
     {when:'UT 2017-04-15 10',text:'Moon at apogee 405474.560 km'},
     {when:'UT 2017-04-16 18',text:'Saturn 3.23°S of Moon'},
     {when:'UT 2017-04-19 10',text:'Last Quarter 396607.178 km'},
     {when:'UT 2017-04-20 06',text:'Mercury in inferior conjunction 1.65° North'},
     {when:'UT 2017-04-22 20',text:'Neptune 0.20°N of Moon'},
     {when:'UT 2017-04-23 18',text:'Venus 5.18°N of Moon'},
     {when:'UT 2017-04-25 16',text:'Uranus 3.71°N of Moon'},
     {when:'UT 2017-04-25 18',text:'Mercury 4.51°N of Moon'},
     {when:'UT 2017-04-26 12',text:'New Moon 360541.901 km'},
     {when:'UT 2017-04-27 16',text:'Moon at perigee 359326.787 km'},
     {when:'UT 2017-04-28 07',text:'Mars 5.77°N of Moon'},
     {when:'UT 2017-04-28 18',text:'Aldebaran 0.49°S of Moon'},
     {when:'UT 2017-05-01 20',text:'Pollux 9.77°N of Moon'},
     {when:'UT 2017-05-03 03',text:'First Quarter 378131.046 km'},
     {when:'UT 2017-05-04 10',text:'Regulus 0.52°N of Moon'},
     {when:'UT 2017-05-07 07',text:'Aldebaran 6.24°S of Mars'},
     {when:'UT 2017-05-07 21',text:'Jupiter 2.11°S of Moon'},
     {when:'UT 2017-05-07 23',text:'Uranus 2.23°N of Mercury'},
     {when:'UT 2017-05-08 13',text:'Spica 6.44°S of Moon'},
     {when:'UT 2017-05-10 22',text:'Full Moon 404917.146 km'},
     {when:'UT 2017-05-12 10',text:'Antares 9.62°S of Moon'},
     {when:'UT 2017-05-12 20',text:'Moon at apogee 406209.826 km'},
     {when:'UT 2017-05-13 23',text:'Saturn 3.07°S of Moon'},
     {when:'UT 2017-05-17 23',text:'Mercury at greatest elongation 25.8° West'},
     {when:'UT 2017-05-19 01',text:'Last Quarter 389073.194 km'},
     {when:'UT 2017-05-20 06',text:'Neptune 0.47°N of Moon'},
     {when:'UT 2017-05-22 13',text:'Venus 2.39°N of Moon'},
     {when:'UT 2017-05-23 05',text:'Uranus 3.89°N of Moon'},
     {when:'UT 2017-05-24 01',text:'Mercury 1.61°N of Moon'},
     {when:'UT 2017-05-25 20',text:'New Moon 357260.997 km'},
     {when:'UT 2017-05-26 01',text:'Moon at perigee 357207.377 km'},
     {when:'UT 2017-05-26 04',text:'Aldebaran 0.57°S of Moon'},
     {when:'UT 2017-05-27 02',text:'Mars 5.35°N of Moon'},
     {when:'UT 2017-05-29 04',text:'Pollux 9.55°N of Moon'},
     {when:'UT 2017-05-31 17',text:'Regulus 0.25°N of Moon'},
     {when:'UT 2017-06-01 13',text:'First Quarter 385082.746 km'},
     {when:'UT 2017-06-02 15',text:'Uranus 1.78°N of Venus'},
     {when:'UT 2017-06-03 12',text:'Venus at greatest elongation 45.9° West'},
     {when:'UT 2017-06-04 00',text:'Jupiter 2.32°S of Moon'},
     {when:'UT 2017-06-04 18',text:'Spica 6.61°S of Moon'},
     {when:'UT 2017-06-08 16',text:'Antares 9.60°S of Moon'},
     {when:'UT 2017-06-08 22',text:'Moon at apogee 406401.132 km'},
     {when:'UT 2017-06-09 13',text:'Full Moon 406271.900 km'},
     {when:'UT 2017-06-10 01',text:'Saturn 3.08°S of Moon'},
     {when:'UT 2017-06-12 11',text:'Aldebaran 5.07°S of Mercury'},
     {when:'UT 2017-06-15 10',text:'Saturn at opposition'},
     {when:'UT 2017-06-16 13',text:'Neptune 0.73°N of Moon'},
     {when:'UT 2017-06-17 12',text:'Last Quarter 381524.476 km'},
     {when:'UT 2017-06-19 16',text:'Uranus 4.14°N of Moon'},
     {when:'UT 2017-06-20 21',text:'Venus 2.37°N of Moon'},
     {when:'UT 2017-06-21 04',text:'Solstice'},
     {when:'UT 2017-06-21 14',text:'Mercury in superior conjunction 1.09° North'},
     {when:'UT 2017-06-22 15',text:'Aldebaran 0.54°S of Moon'},
     {when:'UT 2017-06-23 11',text:'Moon at perigee 357937.026 km'},
     {when:'UT 2017-06-24 03',text:'New Moon 358341.098 km'},
     {when:'UT 2017-06-24 09',text:'Mercury 5.28°N of Moon'},
     {when:'UT 2017-06-24 20',text:'Mars 4.41°N of Moon'},
     {when:'UT 2017-06-25 14',text:'Pollux 9.41°N of Moon'},
     {when:'UT 2017-06-28 01',text:'Regulus 0.03°N of Moon'},
     {when:'UT 2017-06-28 18',text:'Mars 0.78°S of Mercury'},
     {when:'UT 2017-07-01 01',text:'First Quarter 392686.448 km'},
     {when:'UT 2017-07-01 07',text:'Jupiter 2.71°S of Moon'},
     {when:'UT 2017-07-02 00',text:'Spica 6.83°S of Moon'},
     {when:'UT 2017-07-03 00',text:'Pollux 4.87°N of Mercury'},
     {when:'UT 2017-07-03 20',text:'Earth at aphelion 1.016675594 AU'},
     {when:'UT 2017-07-05 22',text:'Antares 9.69°S of Moon'},
     {when:'UT 2017-07-06 04',text:'Moon at apogee 405933.868 km'},
     {when:'UT 2017-07-07 03',text:'Saturn 3.24°S of Moon'},
     {when:'UT 2017-07-09 04',text:'Full Moon 402624.258 km'},
     {when:'UT 2017-07-11 23',text:'Pollux 5.72°N of Mars'},
     {when:'UT 2017-07-13 18',text:'Neptune 0.87°N of Moon'},
     {when:'UT 2017-07-14 11',text:'Aldebaran 3.16°S of Venus'},
     {when:'UT 2017-07-16 19',text:'Last Quarter 375368.280 km'},
     {when:'UT 2017-07-17 00',text:'Uranus 4.33°N of Moon'},
     {when:'UT 2017-07-20 00',text:'Aldebaran 0.44°S of Moon'},
     {when:'UT 2017-07-20 11',text:'Venus 2.73°N of Moon'},
     {when:'UT 2017-07-21 17',text:'Moon at perigee 361236.387 km'},
     {when:'UT 2017-07-23 01',text:'Pollux 9.39°N of Moon'},
     {when:'UT 2017-07-23 10',text:'New Moon 363542.374 km'},
     {when:'UT 2017-07-23 13',text:'Mars 3.11°N of Moon'},
     {when:'UT 2017-07-25 09',text:'Mercury 0.86°S of Moon'},
     {when:'UT 2017-07-25 11',text:'Regulus 0.07°S of Moon'},
     {when:'UT 2017-07-26 09',text:'Regulus 1.09°N of Mercury'},
     {when:'UT 2017-07-27 01',text:'Mars 1.10°N of Sun'},
     {when:'UT 2017-07-28 20',text:'Jupiter 3.14°S of Moon'},
     {when:'UT 2017-07-29 08',text:'Spica 6.98°S of Moon'},
     {when:'UT 2017-07-30 05',text:'Mercury at greatest elongation 27.2° East'},
     {when:'UT 2017-07-30 15',text:'First Quarter 399352.949 km'},
     {when:'UT 2017-08-02 05',text:'Antares 9.79°S of Moon'},
     {when:'UT 2017-08-02 18',text:'Moon at apogee 405025.066 km'},
     {when:'UT 2017-08-03 07',text:'Saturn 3.45°S of Moon'},
     {when:'UT 2017-08-07 18',text:'Partial umbral ECLIPSE Full Moon 394794.690 km'},
     {when:'UT 2017-08-09 23',text:'Neptune 0.86°N of Moon'},
     {when:'UT 2017-08-13 05',text:'Uranus 4.40°N of Moon'},
     {when:'UT 2017-08-15 01',text:'Last Quarter 371377.904 km'},
     {when:'UT 2017-08-16 07',text:'Aldebaran 0.38°S of Moon'},
     {when:'UT 2017-08-18 13',text:'Moon at perigee 366121.404 km'},
     {when:'UT 2017-08-19 05',text:'Venus 2.25°N of Moon'},
     {when:'UT 2017-08-19 10',text:'Pollux 9.41°N of Moon'},
     {when:'UT 2017-08-21 05',text:'Mars 1.55°N of Moon'},
     {when:'UT 2017-08-21 19',text:'Total ECLIPSE New Moon 372114.496 km'},
     {when:'UT 2017-08-21 20',text:'Regulus 0.08°S of Moon'},
     {when:'UT 2017-08-22 06',text:'Mercury 6.16°S of Moon'},
     {when:'UT 2017-08-25 13',text:'Jupiter 3.49°S of Moon'},
     {when:'UT 2017-08-25 17',text:'Spica 7.00°S of Moon'},
     {when:'UT 2017-08-26 21',text:'Mercury in inferior conjunction 4.22° South'},
     {when:'UT 2017-08-29 08',text:'First Quarter 403486.328 km'},
     {when:'UT 2017-08-29 09',text:'Regulus 4.44°N of Mercury'},
     {when:'UT 2017-08-29 13',text:'Antares 9.80°S of Moon'},
     {when:'UT 2017-08-30 11',text:'Moon at apogee 404308.456 km'},
     {when:'UT 2017-08-30 14',text:'Saturn 3.56°S of Moon'},
     {when:'UT 2017-09-02 00',text:'Mars 4.10°N of Mercury'},
     {when:'UT 2017-09-05 02',text:'Regulus 0.75°S of Mars'},
     {when:'UT 2017-09-05 05',text:'Neptune at opposition'},
     {when:'UT 2017-09-05 11',text:'Spica 3.38°S of Jupiter'},
     {when:'UT 2017-09-06 05',text:'Neptune 0.77°N of Moon'},
     {when:'UT 2017-09-06 07',text:'Full Moon 384378.054 km'},
     {when:'UT 2017-09-09 10',text:'Uranus 4.33°N of Moon'},
     {when:'UT 2017-09-10 12',text:'Regulus 0.59°N of Mercury'},
     {when:'UT 2017-09-12 10',text:'Mercury at greatest elongation 17.9° West'},
     {when:'UT 2017-09-12 13',text:'Aldebaran 0.44°S of Moon'},
     {when:'UT 2017-09-13 06',text:'Last Quarter 369900.302 km'},
     {when:'UT 2017-09-13 16',text:'Moon at perigee 369859.593 km'},
     {when:'UT 2017-09-15 17',text:'Pollux 9.36°N of Moon'},
     {when:'UT 2017-09-16 18',text:'Mars 0.06°S of Mercury'},
     {when:'UT 2017-09-18 01',text:'Venus 0.55°N of Moon'},
     {when:'UT 2017-09-18 05',text:'Regulus 0.09°S of Moon'},
     {when:'UT 2017-09-18 20',text:'Mars 0.14°S of Moon'},
     {when:'UT 2017-09-18 23',text:'Mercury 0.03°N of Moon'},
     {when:'UT 2017-09-19 23',text:'Regulus 0.49°S of Venus'},
     {when:'UT 2017-09-20 05',text:'New Moon 382736.983 km'},
     {when:'UT 2017-09-22 02',text:'Spica 6.93°S of Moon'},
     {when:'UT 2017-09-22 08',text:'Jupiter 3.73°S of Moon'},
     {when:'UT 2017-09-22 20',text:'Equinox'},
     {when:'UT 2017-09-25 21',text:'Antares 9.68°S of Moon'},
     {when:'UT 2017-09-27 00',text:'Saturn 3.48°S of Moon'},
     {when:'UT 2017-09-27 07',text:'Moon at apogee 404347.579 km'},
     {when:'UT 2017-09-28 03',text:'First Quarter 403893.440 km'},
     {when:'UT 2017-10-03 12',text:'Neptune 0.75°N of Moon'},
     {when:'UT 2017-10-05 13',text:'Mars 0.22°S of Venus'},
     {when:'UT 2017-10-05 19',text:'Full Moon 373412.207 km'},
     {when:'UT 2017-10-06 16',text:'Uranus 4.21°N of Moon'},
     {when:'UT 2017-10-08 21',text:'Mercury in superior conjunction 1.11° North'},
     {when:'UT 2017-10-09 06',text:'Moon at perigee 366855.308 km'},
     {when:'UT 2017-10-09 19',text:'Aldebaran 0.60°S of Moon'},
     {when:'UT 2017-10-12 12',text:'Last Quarter 371124.397 km'},
     {when:'UT 2017-10-12 22',text:'Pollux 9.19°N of Moon'},
     {when:'UT 2017-10-13 03',text:'Spica 2.94°S of Mercury'},
     {when:'UT 2017-10-15 11',text:'Regulus 0.21°S of Moon'},
     {when:'UT 2017-10-17 10',text:'Mars 1.78°S of Moon'},
     {when:'UT 2017-10-18 00',text:'Venus 1.97°S of Moon'},
     {when:'UT 2017-10-18 15',text:'Jupiter 1.02°N of Mercury'},
     {when:'UT 2017-10-19 10',text:'Spica 6.89°S of Moon'},
     {when:'UT 2017-10-19 18',text:'Uranus at opposition'},
     {when:'UT 2017-10-19 19',text:'New Moon 393502.863 km'},
     {when:'UT 2017-10-20 03',text:'Jupiter 3.92°S of Moon'},
     {when:'UT 2017-10-20 08',text:'Mercury 5.20°S of Moon'},
     {when:'UT 2017-10-23 05',text:'Antares 9.50°S of Moon'},
     {when:'UT 2017-10-24 12',text:'Saturn 3.25°S of Moon'},
     {when:'UT 2017-10-25 02',text:'Moon at apogee 405154.236 km'},
     {when:'UT 2017-10-26 18',text:'Jupiter 1.02°N of Sun'},
     {when:'UT 2017-10-27 22',text:'First Quarter 400261.093 km'},
     {when:'UT 2017-10-30 21',text:'Neptune 0.88°N of Moon'},
     {when:'UT 2017-11-01 15',text:'Spica 3.81°S of Venus'},
     {when:'UT 2017-11-03 01',text:'Uranus 4.19°N of Moon'},
     {when:'UT 2017-11-04 05',text:'Full Moon 364001.518 km'},
     {when:'UT 2017-11-06 00',text:'Moon at perigee 361438.066 km'},
     {when:'UT 2017-11-06 03',text:'Aldebaran 0.76°S of Moon'},
     {when:'UT 2017-11-09 04',text:'Pollux 8.94°N of Moon'},
     {when:'UT 2017-11-10 21',text:'Last Quarter 375108.291 km'},
     {when:'UT 2017-11-11 17',text:'Regulus 0.45°S of Moon'},
     {when:'UT 2017-11-12 15',text:'Antares 2.25°S of Mercury'},
     {when:'UT 2017-11-13 06',text:'Jupiter 0.28°S of Venus'},
     {when:'UT 2017-11-15 01',text:'Mars 3.19°S of Moon'},
     {when:'UT 2017-11-15 16',text:'Spica 6.97°S of Moon'},
     {when:'UT 2017-11-16 21',text:'Jupiter 4.09°S of Moon'},
     {when:'UT 2017-11-17 06',text:'Venus 3.96°S of Moon'},
     {when:'UT 2017-11-18 12',text:'New Moon 402115.462 km'},
     {when:'UT 2017-11-19 12',text:'Antares 9.38°S of Moon'},
     {when:'UT 2017-11-20 09',text:'Mercury 6.90°S of Moon'},
     {when:'UT 2017-11-21 00',text:'Saturn 2.99°S of Moon'},
     {when:'UT 2017-11-21 19',text:'Moon at apogee 406131.507 km'},
     {when:'UT 2017-11-24 00',text:'Mercury at greatest elongation 22.0° East'},
     {when:'UT 2017-11-26 17',text:'First Quarter 393528.825 km'},
     {when:'UT 2017-11-27 05',text:'Neptune 1.15°N of Moon'},
     {when:'UT 2017-11-28 00',text:'Spica 3.36°S of Mars'},
     {when:'UT 2017-11-28 09',text:'Saturn 3.05°N of Mercury'},
     {when:'UT 2017-11-30 10',text:'Uranus 4.32°N of Moon'},
     {when:'UT 2017-12-03 13',text:'Aldebaran 0.82°S of Moon'},
     {when:'UT 2017-12-03 16',text:'Full Moon 357982.592 km'},
     {when:'UT 2017-12-04 09',text:'Moon at perigee 357492.090 km'},
     {when:'UT 2017-12-06 12',text:'Saturn 1.35°N of Mercury'},
     {when:'UT 2017-12-06 13',text:'Pollux 8.72°N of Moon'},
     {when:'UT 2017-12-08 17',text:'Antares 5.08°S of Venus'},
     {when:'UT 2017-12-08 23',text:'Regulus 0.73°S of Moon'},
     {when:'UT 2017-12-10 08',text:'Last Quarter 381514.745 km'},
     {when:'UT 2017-12-12 21',text:'Spica 7.18°S of Moon'},
     {when:'UT 2017-12-13 02',text:'Mercury in inferior conjunction 1.72° North'},
     {when:'UT 2017-12-13 16',text:'Mars 4.16°S of Moon'},
     {when:'UT 2017-12-14 14',text:'Jupiter 4.25°S of Moon'},
     {when:'UT 2017-12-15 16',text:'Venus 2.22°S of Mercury'},
     {when:'UT 2017-12-16 18',text:'Antares 9.40°S of Moon'},
     {when:'UT 2017-12-17 09',text:'Mercury 1.76°S of Moon'},
     {when:'UT 2017-12-17 18',text:'Venus 4.14°S of Moon'},
     {when:'UT 2017-12-18 07',text:'New Moon 406403.017 km'},
     {when:'UT 2017-12-18 13',text:'Saturn 2.78°S of Moon'},
     {when:'UT 2017-12-19 01',text:'Moon at apogee 406602.671 km'},
     {when:'UT 2017-12-21 16',text:'Solstice'},
     {when:'UT 2017-12-21 21',text:'Saturn 0.91°N of Sun'},
     {when:'UT 2017-12-24 13',text:'Neptune 1.43°N of Moon'},
     {when:'UT 2017-12-25 18',text:'Saturn 1.13°N of Venus'},
     {when:'UT 2017-12-26 09',text:'First Quarter 385606.440 km'},
     {when:'UT 2017-12-27 18',text:'Uranus 4.53°N of Moon'},
     {when:'UT 2017-12-31 01',text:'Aldebaran 0.76°S of Moon'},
     {when:'UT 2018-01-01 20',text:'Mercury at greatest elongation 22.7° West'},
     {when:'UT 2018-01-01 22',text:'Moon at perigee 356564.819 km'},
     {when:'UT 2018-01-02 02',text:'Full Moon 356602.153 km'},
     {when:'UT 2018-01-02 23',text:'Pollux 8.63°N of Moon'},
     {when:'UT 2018-01-03 06',text:'Earth at perihelion 0.983284268 AU'},
     {when:'UT 2018-01-05 08',text:'Regulus 0.91°S of Moon'},
     {when:'UT 2018-01-07 04',text:'Jupiter 0.21°N of Mars'},
     {when:'UT 2018-01-08 22',text:'Last Quarter 389330.889 km'},
     {when:'UT 2018-01-09 04',text:'Spica 7.37°S of Moon'},
     {when:'UT 2018-01-09 07',text:'Venus in superior conjunction 0.77° South'},
     {when:'UT 2018-01-11 06',text:'Jupiter 4.34°S of Moon'},
     {when:'UT 2018-01-11 10',text:'Mars 4.56°S of Moon'},
     {when:'UT 2018-01-13 00',text:'Antares 9.49°S of Moon'},
     {when:'UT 2018-01-13 07',text:'Saturn 0.64°N of Mercury'},
     {when:'UT 2018-01-15 02',text:'Saturn 2.63°S of Moon'},
     {when:'UT 2018-01-15 07',text:'Mercury 3.36°S of Moon'},
     {when:'UT 2018-01-17 02',text:'New Moon 405106.265 km'},
     {when:'UT 2018-01-17 08',text:'Venus 2.47°S of Moon'},
     {when:'UT 2018-01-20 20',text:'Neptune 1.61°N of Moon'},
     {when:'UT 2018-01-24 01',text:'Uranus 4.68°N of Moon'},
     {when:'UT 2018-01-24 22',text:'First Quarter 378406.419 km'},
     {when:'UT 2018-01-27 11',text:'Aldebaran 0.68°S of Moon'},
     {when:'UT 2018-01-30 10',text:'Moon at perigee 358993.525 km'},
     {when:'UT 2018-01-30 11',text:'Pollux 8.63°N of Moon'},
     {when:'UT 2018-01-31 13',text:'Total ECLIPSE Full Moon 360197.777 km'}
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

  /* Coords have equinox 2017.5. */
  var build_messiers = function(){
    var messier = [110];
    messier[0]=["M1",1.4641288,0.3844482,8.4,"Tau","NB","!! famous Crab Neb. supernova remnant","Crab Nebula"];
    messier[1]=["M2",5.6478569,0.0156202,6.5,"Aqr","GC","200-mm telescope needed to resolve",""];
    messier[2]=["M3",3.5910408,0.4938497,6.2,"CVn","GC","!! contains many variable stars",""];
    messier[3]=["M4",4.2964546,-0.4637847,5.6,"Sco","GC","bright globular near Antares",""];
    messier[4]=["M5",4.0120156,0.0352625,5.6,"Ser","GC","!! one of the sky's finest globulars",""];
    messier[5]=["M6",4.6305403,-0.5624301,4.2,"Sco","OC","!! Butterfly Cluster; best at low power","Butterfly Cluster"];
    messier[6]=["M7",4.6908685,-0.6077064,3.3,"Sco","OC","!! excellent in binocs or rich-field scope","Ptolemy's Cluster"];
    messier[7]=["M8",4.7336537,-0.4255372,6,"Sgr","NB","!! Lagoon Nebula w/open cl. NGC 6530","Lagoon Nebula"];
    messier[8]=["M9",4.5388399,-0.3234742,7.7,"Oph","GC","smallest of Ophiuchus globulars",""];
    messier[9]=["M10",4.4419673,-0.0720161,6.6,"Oph","GC","rich globular cluster; M12 is 3°NW",""];
    messier[10]=["M11",4.93945,-0.1089946,6.3,"Sct","OC","!! Wild Duck Cl.; the best open cluster?","Wild Duck Cluster"];
    messier[11]=["M12",4.3987081,-0.0345618,6.7,"Oph","GC","loose globular cluster near M10",""];
    messier[12]=["M13",4.3734703,0.6358959,5.8,"Her","GC","!! Hercules Cluster; NGC 6207 0.5°NE","Great Hercules Globular"];
    messier[13]=["M14",4.6186603,-0.0568857,7.6,"Oph","GC","200-mm telescope needed to resolve",""];
    messier[14]=["M15",5.6323769,0.2136994,6.2,"Peg","GC","rich, compact globular","Great Pegasus Globular"];
    messier[15]=["M16",4.7987485,-0.2404215,6.4,"Ser","NB","Eagle Neb. w/open cl.; use neb. filter","Eagle Nebula"];
    messier[16]=["M17",4.8075508,-0.2822947,7,"Sgr","NB","!! Swan or Omega Nebula; use neb. filter","Omega Nebula"];
    messier[17]=["M18",4.8036546,-0.2988819,7.5,"Sgr","OC","sparse cluster; 1°S of M17",""];
    messier[18]=["M19",4.4666616,-0.4588574,6.8,"Oph","GC","oblate globular; M62 4°S",""];
    messier[19]=["M20",4.72837,-0.4019843,9,"Sgr","NB","!! Trifid Nebula; look for dark lanes","Trifid Nebula"];
    messier[20]=["M21",4.7370779,-0.392661,6.5,"Sgr","OC","0.7°NE of M20; sparse cluster",""];
    messier[21]=["M22",4.8758709,-0.4168609,5.1,"Sgr","GC","spectacular from southern latitude","Sagittarius Cluster"];
    messier[22]=["M23",4.7029259,-0.3319233,6.9,"Sgr","OC","bright, loose open cluster",""];
    messier[23]=["M24",4.7906098,-0.3227568,4.6,"Sgr","OC","rich star cloud; best in big binoculars","Sagittarius Star Cloud"];
    messier[24]=["M25",4.8547713,-0.3357384,6.5,"Sgr","OC","bright but sparse open cluster",""];
    messier[25]=["M26",4.9138004,-0.1637242,8,"Sct","OC","bright, coarse cluster",""];
    messier[26]=["M27",5.2375386,0.3973307,7.4,"Vul","NB","!! Dumbbell Nebula; a superb object","Dumbbell Nebula"];
    messier[27]=["M28",4.8239871,-0.4338198,6.8,"Sgr","GC","compact globular near M22",""];
    messier[28]=["M29",5.3430888,0.6735343,7.1,"Cyg","OC","small, poor open cluster 2°S of γ Cygni",""];
    messier[29]=["M30",5.6783935,-0.4032288,7.2,"Cap","GC","toughest in one-night Messier marathon",""];
    messier[30]=["M31",0.1865746,0.7219108,3.4,"And","GY","!! Andromeda Gal.; look for dust lanes","Andromeda Galaxy"];
    messier[31]=["M32",0.1909404,0.7149281,8.1,"And","GY","closest companion to M31",""];
    messier[32]=["M33",0.4140336,0.5365016,5.7,"Tri","GY","large diffuse spiral; requires dark sky","Triangulum Galaxy"];
    messier[33]=["M34",0.7117981,0.7480003,5.5,"Per","OC","best at low power",""];
    messier[34]=["M35",1.6143117,0.4246268,5.3,"Gem","OC","!! look for sm. cluster NGC 2158 0.25°S",""];
    messier[35]=["M36",1.4715734,0.5959118,6.3,"Aur","OC","bright but scattered group; use low pow.",""];
    messier[36]=["M37",1.5426335,0.5681568,6.2,"Aur","OC","!! finest of three Auriga clusters; very rich",""];
    messier[37]=["M38",1.439355,0.6256368,7.4,"Aur","OC","look for small cluster NGC 1907 0.5°S",""];
    messier[38]=["M39",5.6410475,0.8466814,4.6,"Cyg","OC","very sparse cluster; use low power",""];
    messier[39]=["M40",3.2429738,1.0120534,8.4,"UMa","OC","double star Winneke 4; separation 50arcsec","Winnecke 4"];
    messier[40]=["M41",1.7791557,-0.3622139,4.6,"CMa","OC","4°S of Sirius; bright but coarse",""];
    messier[41]=["M42",1.4672109,-0.0949414,4,"Ori","NB","!! Orion Nebula; finest in northern sky","Great Nebula in Orion"];
    messier[42]=["M43",1.468089,-0.0917432,9,"Ori","NB","detached part of Orion Nebula","De Mairan's Nebula"];
    messier[43]=["M44",2.2737497,0.3476786,3.7,"Cnc","OC","!! Beehive or Praesepe; use low power","Beehive Cluster"];
    messier[44]=["M45",0.9950263,0.4218443,1.6,"Tau","OC","!! Pleiades; look for subtle nebulosity","Pleiades"];
    messier[45]=["M46",2.0184897,-0.259333,6,"Pup","OC","!! contains planetary nebula NGC 2438",""];
    messier[46]=["M47",1.9958053,-0.2537711,5.2,"Pup","OC","coarse cluster 1.5°W of M46",""];
    messier[47]=["M48",2.1583779,-0.1021691,5.5,"Hya","OC","former lost Messier; large sparse cl.",""];
    messier[48]=["M49",3.275502,0.1379406,8.4,"Vir","GY","very bright elliptical",""];
    messier[49]=["M50",1.8502319,-0.1459101,6.3,"Mon","OC","between Sirius & Procyon; use low mag",""];
    messier[50]=["M51",3.5375013,0.8219345,8.4,"CVn","GY","!! Whirlpool Galaxy; superb in big scope","Whirlpool Galaxy"];
    messier[51]=["M52",6.1304072,1.0765121,7.3,"Cas","OC","young, rich cl.; faint Bubble Neb. nearby",""];
    messier[52]=["M53",3.4634174,0.315454,7.6,"Com","GC","150-mm telescope needed to resolve",""];
    messier[53]=["M54",4.9576928,-0.5316256,7.6,"Sgr","GC","not easily resolved",""];
    messier[54]=["M55",5.1535578,-0.5397479,6.3,"Sgr","GC","bright, loose globular cluster",""];
    messier[55]=["M56",5.0495986,0.5273587,8.3,"Lyr","GC","within a rich starfield",""];
    messier[56]=["M57",4.9491009,0.5769368,8.8,"Lyr","NB","!! Ring Nebula; an amazing smoke ring","Ring Nebula"];
    messier[57]=["M58",3.3099447,0.2045628,9.7,"Vir","GY","bright barred spiral; M59 and M60 1°E",""];
    messier[58]=["M59",3.3287014,0.2016595,9.6,"Vir","GY","bright elliptical paired with M60",""];
    messier[59]=["M60",3.3361172,0.1999165,8.8,"Vir","GY","bright elliptical with M59 and NGC 4647",""];
    messier[60]=["M61",3.2410501,0.0762656,9.7,"Vir","GY","face-on two-armed spiral",""];
    messier[61]=["M62",4.4606943,-0.5260625,6.5,"Oph","GC","asymmetrical; in rich field",""];
    messier[62]=["M63",3.4757466,0.7320128,8.6,"CVn","GY","!! Sunflower Galaxy; bright, elongated","Sunflower Galaxy"];
    messier[63]=["M64",3.3927402,0.3767977,8.5,"Com","GY","!! Black Eye Gal; eye needs big scope","Black Eye Galaxy"];
    messier[64]=["M65",2.9662431,0.2266734,9.3,"Leo","GY","!! bright elongated spiral",""];
    messier[65]=["M66",2.9719126,0.2249265,8.9,"Leo","GY","!! M65 and NGC 3628 in same field",""];
    messier[66]=["M67",2.3184807,0.2050861,6.1,"Cnc","OC","one of the oldest star clusters known",""];
    messier[67]=["M68",3.3180064,-0.4685503,7.8,"Hya","GC","150-mm telescope needed to resolve",""];
    messier[68]=["M69",4.8543771,-0.5643776,7.6,"Sgr","GC","small, poor globular cluster",""];
    messier[69]=["M70",4.9058531,-0.5634185,7.9,"Sgr","GC","small globular 2°E of M69",""];
    messier[70]=["M71",5.21234,0.3286436,8.2,"Sge","GC","loose globular; looks like an open cluster",""];
    messier[71]=["M72",5.4736125,-0.2175775,9.3,"Aqr","GC","near the Saturn Nebula, NGC 7009",""];
    messier[72]=["M73",5.4976067,-0.2192936,9,"Aqr","OC","group of four stars only; an asterism",""];
    messier[73]=["M74",0.4260451,0.277021,9.4,"Psc","GY","faint, elusive spiral; tough in small scope",""];
    messier[74]=["M75",5.2670991,-0.3816256,8.5,"Sgr","GC","small and distant; 59,000 ly away",""];
    messier[75]=["M76",0.4516499,0.9015398,10.1,"Per","NB","Little Dumbell; faint but distinct","Little Dumbbell Nebula"];
    messier[76]=["M77",0.7138276,0.0018693,8.9,"Cet","GY","a Seyfert galaxy; with starlike nucleus",""];
    messier[77]=["M78",1.5166792,0.000968,8.3,"Ori","NB","bright featureless reflection nebula",""];
    messier[78]=["M79",1.4190446,-0.4282186,7.7,"Lep","GC","200-mm telescope needed to resolve",""];
    messier[79]=["M80",4.2675312,-0.4018701,7.3,"Sco","GC","very compressed globular",""];
    messier[80]=["M81",2.6049882,1.203982,6.9,"UMa","GY","!! bright spiral visible in binoculars","Bode's Galaxy"];
    messier[81]=["M82",2.6059327,1.2147441,8.4,"UMa","GY","!! the exploding galaxy; M81 0.5°S","Cigar Galaxy"];
    messier[82]=["M83",3.5691522,-0.5228206,7.6,"Hya","GY","large and diffuse; superb from far south","Southern Pinwheel"];
    messier[83]=["M84",3.2549825,0.2231667,9.1,"Vir","GY","!! w/M86 in Markarian's Chain",""];
    messier[84]=["M85",3.2567079,0.3159603,9.1,"Com","GY","bright elliptical shape",""];
    messier[85]=["M86",3.25978,0.2243311,8.9,"Vir","GY","!! w/many NGC galaxies in Chain",""];
    messier[86]=["M87",3.2798459,0.2147361,8.6,"Vir","GY","famous jet and black hole",""];
    messier[87]=["M88",3.2855071,0.2502259,9.6,"Com","GY","bright multiple-arm spiral",""];
    messier[88]=["M89",3.3012176,0.2173595,9.8,"Vir","GY","elliptical; resembles M87 but smaller",""];
    messier[89]=["M90",3.3060123,0.2281237,9.5,"Vir","GY","bright barred spiral near M89",""];
    messier[90]=["M91",3.3003357,0.2513931,10.2,"Com","GY","some lists say M91=M58, not NGC 4548",""];
    messier[91]=["M92",4.5275508,0.7525042,6.4,"Her","GC","9°NE of M13; fine but often overlooked",""];
    messier[92]=["M93",2.0304379,-0.4173038,6,"Pup","OC","compact, bright cluster; fairly rich",""];
    messier[93]=["M94",3.3672701,0.7162541,8.2,"CVn","GY","very bright and comet-like",""];
    messier[94]=["M95",2.8140071,0.2025946,9.7,"Leo","GY","bright barred spiral",""];
    messier[95]=["M96",2.8262215,0.2046242,9.2,"Leo","GY","M95 in same field",""];
    messier[96]=["M97",2.948754,0.9585538,9.9,"UMa","NB","Owl Nebula; distinct grey oval","Owl Nebula"];
    messier[97]=["M98",3.2061281,0.2586478,10.1,"Com","GY","nearly edge-on spiral near star 6 Com. B.",""];
    messier[98]=["M99",3.2279362,0.2502148,9.9,"Com","GY","nearly face-on spiral near M98",""];
    messier[99]=["M100",3.2458135,0.2746522,9.3,"Com","GY","face-on spiral with starlike nucleus",""];
    messier[100]=["M101",3.6818528,0.947127,7.9,"UMa","GY","!! Pinwheel Gal; diffuse face-on spiral","Pinwheel Galaxy"];
    messier[101]=["M102",3.9574499,0.9721454,9.9,"Dra","GY","or is M102=M101? (look for NGC 5907)",""];
    messier[102]=["M103",0.4117831,1.0609749,7.4,"Cas","OC","three NGC open clusters nearby",""];
    messier[103]=["M104",3.3201008,-0.2044231,8,"Vir","GY","!! Sombrero Galaxy; look for dust lane","Sombrero Galaxy"];
    messier[104]=["M105",2.8305906,0.2180028,9.3,"Leo","GY","bright elliptical near M95 and M96",""];
    messier[105]=["M106",3.2278179,0.8241372,8.4,"CVn","GY","!! superb large, bright spiral",""];
    messier[106]=["M107",4.3348783,-0.2283957,7.9,"Oph","GC","small, faint globular",""];
    messier[107]=["M108",2.9344016,0.9699033,10,"UMa","GY","nearly edge-on; paired with M97 0.75°SE",""];
    messier[108]=["M109",3.1350536,0.9300145,9.8,"UMa","GY","barred spiral near γ UMA",""];
    messier[109]=["M110",0.1804609,0.7291849,8.5,"And","GY","more distant companion to M31",""];    
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
  
  /* Coords have equinox 2017.5. */
  var build_caldwells = function(){
    var caldwell = [42];
    caldwell[0]=["C68",4.9876234,-0.6444413,9.7,"CrA","Bn","","","NGC_6729"];
    caldwell[1]=["C69",4.5155416,-0.6478541,12.8,"Sco","Pl","","Bug Nebula","NGC_6302"];
    caldwell[2]=["C70",0.2431465,-0.656047,8.1,"Scl","Sp","","","NGC_300"];
    caldwell[3]=["C71",2.0635153,-0.6736267,5.8,"Pup","Oc","","","NGC_2477"];
    caldwell[4]=["C72",0.0688345,-0.6821816,8.2,"Scl","Sb","","","NGC_55"];
    caldwell[5]=["C73",1.3730326,-0.6986682,7.3,"Col","Gc","","","NGC_1851"];
    caldwell[6]=["C74",2.6548242,-0.7071964,8.2,"Vel","Pl","","Eight Burst Nebula","NGC_3132"];
    caldwell[7]=["C75",4.305746,-0.7104439,5.8,"Sco","Oc","","","NGC_6124"];
    caldwell[8]=["C76",4.4297828,-0.7300262,2.6,"Sco","Oc","","","NGC_6231"];
    caldwell[9]=["C77",3.5191528,-0.7523646,7,"Cen","Px","","Centaurus A","Centaurus_A"];
    caldwell[10]=["C78",4.7528328,-0.7626449,6.6,"CrA","Gc","","","NGC_6541"];
    caldwell[11]=["C79",2.6979313,-0.8116584,6.7,"Vel","Gc","","","NGC_3201"];
    caldwell[12]=["C80",3.5249334,-0.830319,3.6,"Cen","Gc","","Omega Centauri","Omega_Centauri"];
    caldwell[13]=["C81",4.5676638,-0.8452803,8.1,"Ara","Gc","","","NGC_6352"];
    caldwell[14]=["C82",4.3747388,-0.8517068,5.2,"Ara","Oc","","","NGC_6193"];
    caldwell[15]=["C83",3.4314325,-0.8649868,9.5,"Cen","Sp","","","NGC_4945"];
    caldwell[16]=["C84",3.6107225,-0.8980361,7.6,"Cen","Gc","","","NGC_5286"];
    caldwell[17]=["C85",2.2719823,-0.9272837,2.5,"Vel","Oc","","Omicron Vel Cluster","IC_2391"];
    caldwell[18]=["C86",4.6343951,-0.9367978,5.6,"Ara","Gc","","","NGC_6397"];
    caldwell[19]=["C87",0.8411595,-0.9625778,8.4,"Hor","Gc","","","NGC_1261"];
    caldwell[20]=["C88",3.9575816,-0.9715717,7.9,"Cir","Oc","","","NGC_5823"];
    caldwell[21]=["C89",4.2776274,-1.0112667,5.4,"Nor","Oc","","S Norma Cluster","NGC_6087"];
    caldwell[22]=["C90",2.4517248,-1.0191282,9.7,"Car","Pl","","","NGC_2867"];
    caldwell[23]=["C91",2.910988,-1.0255813,3,"Car","Oc","","","NGC_3532"];
    caldwell[24]=["C92",2.8120671,-1.0464785,6.2,"Car","Bn","","Eta Carinae Nebula","Carina_Nebula"];
    caldwell[25]=["C93",5.0284609,-1.0463836,5.4,"Pav","Gc","","","NGC_6752"];
    caldwell[26]=["C94",3.3800802,-1.0546685,4.2,"Cru","Oc","","","Jewel Box"];
    caldwell[27]=["C95",4.2114821,-1.0567457,5.1,"TrA","Oc","","","NGC_6025"];
    caldwell[28]=["C96",2.0882359,-1.0631639,3.8,"Car","Oc","","","NGC_2516"];
    caldwell[29]=["C97",3.0409001,-1.0771052,5.3,"Cen","Oc","","","NGC_3766"];
    caldwell[30]=["C98",3.330695,-1.1006465,6.9,"Cru","Oc","","","NGC_4609"];
    caldwell[31]=["C99",3.3775364,-1.1012117,"","Cru","Dn","","Coalsack Nebula","Coalsack_Nebula"];
    caldwell[32]=["C100",3.043069,-1.1018311,4.5,"Cen","Oc","","Lambda Centauri Nebula","IC_2944"];
    caldwell[33]=["C101",5.0241605,-1.1138769,9,"Pav","Sb","","","NGC_6744"];
    caldwell[34]=["C102",2.8092379,-1.1255986,1.9,"Car","Oc","","Theta Car Cluster","IC_2602"];
    caldwell[35]=["C103",1.4773384,-1.2058643,1,"Dor","Bn","","Tarantula Nebula","Tarantula_Nebula"];
    caldwell[36]=["C104",0.2783397,-1.2349302,6.6,"Tuc","Gc","","","NGC_362"];
    caldwell[37]=["C105",3.4068374,-1.2387897,7.3,"Mus","Gc","","","NGC_4833"];
    caldwell[38]=["C106",0.1085104,-1.2564008,4,"Tuc","Gc","","47 Tucanae","47_Tucanae"];
    caldwell[39]=["C107",4.3101474,-1.2608003,9.3,"Aps","Gc","","","NGC_6101"];
    caldwell[40]=["C108",3.2587061,-1.2699618,7.8,"Mus","Gc","","","NGC_4372"];
    caldwell[41]=["C109",2.6584261,-1.4128958,11.6,"Cha","Pl","","","NGC_3195"];
    return caldwell;
  };
  var caldwells = build_caldwells();

  /* Find a Caldwell object either by 'C103' or 'Tarantula Nebula' (case-insensitive). */  
  var find_caldwell = function(name_raw){
    return find_deep_sky_object(name_raw, caldwells);
  };
  
  /* This exist in order to avoid creating a large number of identical objects. */
  var when_fixed_equinox = when("J2017.5");
  
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

  /* Coords have equinox 2017.5. */
  var build_stars = function(){
    /* Source: Yale Bright Star Catalog r5.  Name, Right Ascension, Declination (J2017.5), and Magnitude.*/
    var ybs = [9096];
    //removal of leading whitespace cuts down on the overall size of this js file
ybs[0]=['',0.0264922,0.7910978,6.7];
ybs[1]=['',0.0260062,-0.0070801,6.29];
ybs[2]=['33 Psc',0.0271876,-0.0979148,4.61];
ybs[3]=['86 Peg',0.0287953,0.235506,5.51];
ybs[4]=['',0.0313384,1.0216119,5.96];
ybs[5]=['',0.0314175,-0.8548206,5.7];
ybs[6]=['10 Cas',0.0321267,1.1221332,5.59];
ybs[7]=['',0.0327987,0.5082184,6.13];
ybs[8]=['',0.0337138,-0.4016024,6.18];
ybs[9]=['',0.0357624,-0.3017502,6.19];
ybs[10]=['',0.0376612,-0.0427871,6.43];
ybs[11]=['',0.037835,-0.3911549,5.94];
ybs[12]=['',0.0390329,-0.5834999,5.68];
ybs[13]=['',0.0396973,-0.0410226,6.07];
ybs[14]=['α And',0.0405511,0.5094252,2.06];
ybs[15]=['',0.0400754,-0.1523067,5.99];
ybs[16]=['',0.0418522,0.640955,6.19];
ybs[17]=['',0.0412279,-0.3050861,6.06];
ybs[18]=['',0.0426492,0.4461083,6.23];
ybs[19]=['',0.0450562,1.3929833,6.01];
ybs[20]=['β Cas',0.0440815,1.0340563,2.27];
ybs[21]=['87 Peg',0.0433811,0.3195574,5.53];
ybs[22]=['',0.0432613,-0.9408127,6.33];
ybs[23]=['κ1 Scl',0.0446721,-0.48678,5.42];
ybs[24]=['ε Phe',0.0449047,-0.7967456,3.88];
ybs[25]=['34 Psc',0.0477294,0.1962253,5.51];
ybs[26]=['22 And',0.0490335,0.8058106,5.03];
ybs[27]=['',0.0498328,0.9994257,6.74];
ybs[28]=['',0.0489066,-0.0899069,5.84];
ybs[29]=['γ3 Oct',0.0471371,-1.4333789,5.28];
ybs[30]=['',0.0506408,-0.217864,5.85];
ybs[31]=['',0.0500825,-1.2763091,6.64];
ybs[32]=['6 Cet',0.0530424,-0.2682703,4.89];
ybs[33]=['κ2 Scl',0.0543648,-0.4834986,5.41];
ybs[34]=['θ Scl',0.0550464,-0.6114895,5.25];
ybs[35]=['',0.0563111,0.8421176,6.16];
ybs[36]=['',0.0569704,-0.3113851,5.25];
ybs[37]=['',0.0600149,0.6595704,6.73];
ybs[38]=['γ Peg',0.061697,0.2667014,2.83];
ybs[39]=['',0.0624344,0.4727132,6.3];
ybs[40]=['23 And',0.0629668,0.717898,5.72];
ybs[41]=['',0.0636472,-0.4524714,5.94];
ybs[42]=['',0.0637992,-0.4570577,6.31];
ybs[43]=['',0.0652378,0.581253,6.25];
ybs[44]=['χ Peg',0.0676738,0.3543697,4.8];
ybs[45]=['',0.0669922,-0.1340995,5.12];
ybs[46]=['',0.0609039,-1.4817307,5.77];
ybs[47]=['7 Cet',0.0677544,-0.3287425,4.44];
ybs[48]=['',0.0691266,0.3906287,6.24];
ybs[49]=['35 Psc',0.069294,0.1556492,5.79];
ybs[50]=['',0.0689443,-0.1653266,5.75];
ybs[51]=['',0.0699434,0.5521006,6.45];
ybs[52]=['',0.0701941,0.4778756,6.35];
ybs[53]=['',0.0691528,-0.6075009,6.17];
ybs[54]=['',0.0752824,1.3447413,6.35];
ybs[55]=['',0.0754164,0.7625673,6.15];
ybs[56]=['',0.0742989,-0.547147,5.67];
ybs[57]=['',0.0728975,-1.3232076,6.49];
ybs[58]=['36 Psc',0.0762249,0.1455109,6.11];
ybs[59]=['',0.0781178,1.0756548,5.74];
ybs[60]=['',0.0767707,-0.351045,6.47];
ybs[61]=['',0.0788968,0.8385371,5.89];
ybs[62]=['θ And',0.0785944,0.6768179,4.61];
ybs[63]=['',0.0766505,-1.3732844,6.77];
ybs[64]=['',0.0813863,0.8993713,6.14];
ybs[65]=['',0.0804148,-0.3308094,6.45];
ybs[66]=['',0.0815629,0.0311718,6.17];
ybs[67]=['σ And',0.0839902,0.6437189,4.52];
ybs[68]=['',0.0837317,0.1972735,6.05];
ybs[69]=['26 And',0.0856513,0.7659936,6.11];
ybs[70]=['',0.0853256,0.5517738,5.87];
ybs[71]=['',0.0854732,-0.138853,6.46];
ybs[72]=['',0.0854181,-0.7529035,6.33];
ybs[73]=['ι Cet',0.0886629,-0.1523118,3.56];
ybs[74]=['',0.0899707,0.7125617,6.33];
ybs[75]=['',0.0917328,0.8545536,6.52];
ybs[76]=['ζ Tuc',0.0911698,-1.1305838,4.23];
ybs[77]=['',0.0930472,0.5416255,5.9];
ybs[78]=['',0.0945907,0.5761053,5.79];
ybs[79]=['41 Psc',0.0938133,0.1446406,5.37];
ybs[80]=['',0.0951741,0.1932769,6.56];
ybs[81]=['ρ And',0.0961992,0.6643702,5.18];
ybs[82]=['π Tuc',0.0935973,-1.2134922,5.51];
ybs[83]=['ι Scl',0.0977223,-0.5041329,5.18];
ybs[84]=['',0.0988504,-0.3483818,5.12];
ybs[85]=['42 Psc',0.1018018,0.237006,6.23];
ybs[86]=['',0.0969005,-1.3496624,5.97];
ybs[87]=['9 Cet',0.1036364,-0.2114033,6.39];
ybs[88]=['',0.105081,-0.5399909,6.55];
ybs[89]=['',0.1089237,0.6749902,7.39];
ybs[90]=['',0.1100027,0.9096108,5.57];
ybs[91]=['12 Cas',0.1124379,1.0808466,5.4];
ybs[92]=['',0.1107861,-0.0370414,6.07];
ybs[93]=['',0.1137143,0.9275337,5.74];
ybs[94]=['44 Psc',0.1147632,0.0355441,5.77];
ybs[95]=['β Hyi',0.1154253,-1.3466502,2.8];
ybs[96]=['α Phe',0.1184165,-0.736692,2.39];
ybs[97]=['κ Phe',0.1180592,-0.7606708,3.94];
ybs[98]=['10 Cet',0.1200797,0.0008208,6.19];
ybs[99]=['',0.1226945,-0.444195,5.98];
ybs[100]=['47 Psc',0.1263657,0.31398,5.06];
ybs[101]=['',0.1272911,0.7765163,5.17];
ybs[102]=['η Scl',0.1256377,-0.5743973,4.81];
ybs[103]=['48 Psc',0.1270729,0.2887066,6.06];
ybs[104]=['',0.1275868,0.1795312,6.04];
ybs[105]=['',0.1275423,-0.3532257,6.43];
ybs[106]=['',0.127828,-0.6949612,5.43];
ybs[107]=['',0.1303663,0.645713,6.26];
ybs[108]=['',0.1289595,-0.8802767,6.26];
ybs[109]=['',0.1398269,1.3459273,6.21];
ybs[110]=['',0.1366558,1.0484851,5.94];
ybs[111]=['28 And',0.1354737,0.5209498,5.23];
ybs[112]=['',0.1341648,-0.257743,6.14];
ybs[113]=['',0.1338658,-0.5588559,6.57];
ybs[114]=['12 Cet',0.1349723,-0.0673812,5.72];
ybs[115]=['',0.1363638,-0.4134899,5.19];
ybs[116]=['',0.136637,-0.7128431,6.19];
ybs[117]=['',0.1364562,-0.8398254,5.69];
ybs[118]=['13 Cas',0.1415613,1.1626673,6.18];
ybs[119]=['',0.1411951,0.5877947,5.87];
ybs[120]=['λ Cas',0.1428861,0.9532759,4.73];
ybs[121]=['',0.1424871,0.923906,5.6];
ybs[122]=['λ1 Phe',0.1407261,-0.8500996,4.77];
ybs[123]=['β1 Tuc',0.1410927,-1.0971414,4.37];
ybs[124]=['β2 Tuc',0.1411578,-1.0972772,4.54];
ybs[125]=['',0.1457202,0.760809,6.7];
ybs[126]=['',0.1500348,1.2405457,6.42];
ybs[127]=['κ Cas',0.1483889,1.1000471,4.16];
ybs[128]=['52 Psc',0.1462121,0.3558877,5.38];
ybs[129]=['51 Psc',0.1453004,0.1230803,5.67];
ybs[130]=['',0.1461767,0.4830543,6.67];
ybs[131]=['',0.1472432,0.4952665,6.3];
ybs[132]=['',0.1490142,0.9597806,5.93];
ybs[133]=['β3 Tuc',0.1462448,-1.0984177,5.09];
ybs[134]=['16 Cas',0.1546796,1.1666928,6.48];
ybs[135]=['',0.1507493,-0.5142086,5.55];
ybs[136]=['θ Tuc',0.1488639,-1.2421463,6.13];
ybs[137]=['',0.1539544,-0.9124015,5.57];
ybs[138]=['',0.1563504,0.2350501,6.4];
ybs[139]=['13 Cet',0.1576971,-0.0610259,5.2];
ybs[140]=['14 Cet',0.159013,-0.0071441,5.93];
ybs[141]=['',0.1619722,0.9470994,5.08];
ybs[142]=['',0.160671,0.2321789,6.41];
ybs[143]=['',0.1634585,1.0545677,5.79];
ybs[144]=['λ2 Phe',0.159323,-0.8360932,5.51];
ybs[145]=['',0.1586824,-0.9551439,6.06];
ybs[146]=['',0.1625879,0.4773632,6.5];
ybs[147]=['',0.1611394,-0.2596599,6.45];
ybs[148]=['',0.1613814,-0.3969979,6.06];
ybs[149]=['',0.1646528,0.7781507,5.13];
ybs[150]=['ζ Cas',0.1656129,0.9423569,3.66];
ybs[151]=['π And',0.1650246,0.5901932,4.36];
ybs[152]=['53 Psc',0.1645081,0.2675207,5.89];
ybs[153]=['',0.1660041,0.4208039,6.47];
ybs[154]=['',0.1670898,0.6195142,5.48];
ybs[155]=['',0.1798063,1.441464,6.4];
ybs[156]=['',0.1667334,-0.4305922,5.57];
ybs[157]=['',0.1631246,-1.1349625,6.42];
ybs[158]=['',0.1675901,0.0563981,6.39];
ybs[159]=['',0.1662845,-0.9476798,6.41];
ybs[160]=['ε And',0.1723035,0.513261,4.37];
ybs[161]=['',0.1751447,0.8630727,5.43];
ybs[162]=['δ And',0.1756916,0.5402981,3.27];
ybs[163]=['54 Psc',0.175783,0.372567,5.87];
ybs[164]=['55 Psc',0.1782438,0.3758436,5.36];
ybs[165]=['α Cas',0.1811231,0.988434,2.23];
ybs[166]=['',0.1717392,-1.2748094,6.85];
ybs[167]=['',0.1780933,-0.5910688,6.69];
ybs[168]=['',0.1775699,-0.7801751,6.01];
ybs[169]=['',0.180437,-0.2866017,6.49];
ybs[170]=['',0.1807059,-0.4137926,6.14];
ybs[171]=['',0.181507,-0.0742827,5.91];
ybs[172]=['32 And',0.1835864,0.6903552,5.33];
ybs[173]=['',0.1798566,-1.0360023,5.89];
ybs[174]=['',0.1881325,1.1561629,5.83];
ybs[175]=['',0.1855703,0.431532,6.04];
ybs[176]=['ξ Cas',0.1878384,0.8832806,4.8];
ybs[177]=['μ Phe',0.1839158,-0.8026627,4.59];
ybs[178]=['',0.1899588,1.0271097,6.17];
ybs[179]=['ξ Phe',0.1857149,-0.9844684,5.7];
ybs[180]=['π Cas',0.1939284,0.8224055,4.94];
ybs[181]=['λ1 Scl',0.1900405,-0.6696414,6.06];
ybs[182]=['',0.1896577,-1.0501085,5.98];
ybs[183]=['ρ Tuc',0.1885475,-1.1409623,5.39];
ybs[184]=['β Cet',0.1940057,-0.3122574,2.04];
ybs[185]=['',0.1981863,0.8370552,5.67];
ybs[186]=['',0.1951178,-0.2079743,6.02];
ybs[187]=['η Phe',0.1925734,-1.0012499,4.36];
ybs[188]=['21 Cas',0.2043742,1.3104544,5.66];
ybs[189]=['ο Cas',0.1994376,0.84439,4.54];
ybs[190]=['φ1 Cet',0.1966674,-0.1835014,4.76];
ybs[191]=['λ2 Scl',0.1965192,-0.6689163,5.9];
ybs[192]=['',0.2020008,0.9654665,5.42];
ybs[193]=['',0.1989944,-0.3824116,5.24];
ybs[194]=['',0.1997443,-0.7431811,5.94];
ybs[195]=['',0.1976139,-1.0891241,6.07];
ybs[196]=['',0.2083864,1.211614,6.33];
ybs[197]=['',0.2019885,-0.0791277,6.15];
ybs[198]=['',0.199808,-0.9358364,6.15];
ybs[199]=['18 Cet',0.2022731,-0.2231466,6.15];
ybs[200]=['',0.206223,0.9669244,6.52];
ybs[201]=['',0.2057618,0.7846442,6.05];
ybs[202]=['',0.2031957,-0.2849897,6.47];
ybs[203]=['',0.2083033,1.0414347,6.39];
ybs[204]=['23 Cas',0.2136619,1.307998,5.41];
ybs[205]=['',0.2032079,-0.828267,5.8];
ybs[206]=['',0.2053424,-0.3914167,5.5];
ybs[207]=['57 Psc',0.2071224,0.2717642,5.38];
ybs[208]=['',0.2151686,1.2700801,5.87];
ybs[209]=['58 Psc',0.2091733,0.2106479,5.5];
ybs[210]=['59 Psc',0.210104,0.3433799,6.13];
ybs[211]=['ζ And',0.210625,0.4252065,4.06];
ybs[212]=['60 Psc',0.2107479,0.1193132,5.99];
ybs[213]=['61 Psc',0.2131115,0.3668777,6.54];
ybs[214]=['',0.2120233,-0.3135677,5.7];
ybs[215]=['η Cas',0.2187341,1.0107374,3.44];
ybs[216]=['',0.2132912,-0.3774665,5.57];
ybs[217]=['62 Psc',0.2146647,0.1290711,5.93];
ybs[218]=['',0.2150591,0.0938251,5.75];
ybs[219]=['ν Cas',0.2174452,0.8912264,4.89];
ybs[220]=['δ Psc',0.2163831,0.1340447,4.43];
ybs[221]=['64 Psc',0.217733,0.2973295,5.07];
ybs[222]=['ν And',0.2215885,0.7186216,4.53];
ybs[223]=['',0.2194898,-0.2350305,5.59];
ybs[224]=['',0.2185638,-0.4195988,5.9];
ybs[225]=['',0.2170915,-0.8133688,6.27];
ybs[226]=['65 Psc',0.2217512,0.4853049,7];
ybs[227]=['65 Psc',0.2217802,0.4852952,7.1];
ybs[228]=['',0.2199866,-0.4060778,6.28];
ybs[229]=['',0.2260339,1.1229885,5.39];
ybs[230]=['',0.223785,0.7870957,6.15];
ybs[231]=['φ2 Cet',0.2225627,-0.1841214,5.19];
ybs[232]=['λ Hyi',0.2145956,-1.305997,5.07];
ybs[233]=['',0.2283479,1.0803725,6.07];
ybs[234]=['',0.2267308,0.900643,6.39];
ybs[235]=['',0.2219981,-0.7557214,6.48];
ybs[236]=['',0.2471378,1.462617,5.62];
ybs[237]=['',0.2293772,0.9017424,6.21];
ybs[238]=['ρ Phe',0.2246126,-0.8882317,5.22];
ybs[239]=['',0.2277966,0.0607367,6.37];
ybs[240]=['',0.2361838,1.0684722,4.82];
ybs[241]=['',0.2298649,-0.7612125,6.9];
ybs[242]=['',0.2350029,0.6744547,6.69];
ybs[243]=['',0.2335853,-0.4173258,5.46];
ybs[244]=['20 Cet',0.2351985,-0.0183151,4.77];
ybs[245]=['',0.2375244,0.6547218,6.06];
ybs[246]=['',0.2391558,0.9212523,6.27];
ybs[247]=['',0.2358898,-0.4307852,6.46];
ybs[248]=['λ1 Tuc',0.2315385,-1.2114258,6.22];
ybs[249]=['υ1 Cas',0.244583,1.0309199,4.83];
ybs[250]=['66 Psc',0.242234,0.3365512,5.74];
ybs[251]=['21 Cet',0.2407514,-0.1509041,6.16];
ybs[252]=['',0.2447397,0.8512527,6.27];
ybs[253]=['',0.2371527,-1.0956592,5.7];
ybs[254]=['36 And',0.2439294,0.4140431,5.47];
ybs[255]=['',0.2451527,0.4302499,6.2];
ybs[256]=['',0.2498657,1.0138814,6.21];
ybs[257]=['',0.2533926,1.2020168,6.37];
ybs[258]=['67 Psc',0.2483641,0.4765435,6.09];
ybs[259]=['',0.2469269,-0.1265835,5.85];
ybs[260]=['γ Cas',0.2521013,1.0613534,2.47];
ybs[261]=['υ2 Cas',0.2518592,1.034553,4.63];
ybs[262]=['',0.252426,1.0551767,5.55];
ybs[263]=['φ3 Cet',0.2482864,-0.1949913,5.31];
ybs[264]=['',0.2477148,-0.4831256,6.1];
ybs[265]=['μ And',0.2518815,0.6735898,3.87];
ybs[266]=['λ2 Tuc',0.2428319,-1.211823,5.45];
ybs[267]=['η And',0.2537085,0.4103594,4.42];
ybs[268]=['',0.2559503,0.8017001,6.12];
ybs[269]=['',0.2602401,1.159709,5.97];
ybs[270]=['68 Psc',0.2565112,0.5076554,5.42];
ybs[271]=['',0.2583095,0.5941987,5.98];
ybs[272]=['',0.2566911,0.240683,6.32];
ybs[273]=['',0.2585302,0.3752228,6.37];
ybs[274]=['',0.2692709,1.2405283,6.39];
ybs[275]=['φ4 Cet',0.2600917,-0.1969744,5.61];
ybs[276]=['α Scl',0.2593901,-0.5107406,4.31];
ybs[277]=['',0.2578487,-1.0577069,6.23];
ybs[278]=['',0.2664005,0.7819975,6.84];
ybs[279]=['',0.2664151,0.7820363,6.04];
ybs[280]=['',0.2650145,0.1147927,6.11];
ybs[281]=['',0.3117957,1.5038744,4.25];
ybs[282]=['',0.4568617,1.5549287,6.46];
ybs[283]=['',0.2763477,0.8923658,6.47];
ybs[284]=['ξ Scl',0.2710424,-0.6775848,5.59];
ybs[285]=['',0.2794184,0.8285046,6.45];
ybs[286]=['39 And',0.2787976,0.7232422,5.98];
ybs[287]=['σ Psc',0.2782983,0.5567283,5.5];
ybs[288]=['',0.2823455,1.067594,5.92];
ybs[289]=['σ Scl',0.2760771,-0.5490484,5.5];
ybs[290]=['ε Psc',0.2786204,0.1393423,4.28];
ybs[291]=['ω Phe',0.2738675,-0.9932435,6.11];
ybs[292]=['25 Cet',0.2789529,-0.0827802,5.43];
ybs[293]=['',0.2854709,1.0764114,5.84];
ybs[294]=['',0.2839578,0.91797,5.99];
ybs[295]=['',0.2775321,-0.8081529,5.36];
ybs[296]=['',0.279819,-0.5136878,6.29];
ybs[297]=['26 Cet',0.2823778,0.0254868,6.04];
ybs[298]=['',0.287161,0.8919244,6.54];
ybs[299]=['',0.2854445,0.519273,6.19];
ybs[300]=['',0.27655,-1.1407882,6.21];
ybs[301]=['',0.2862137,0.6996089,6.72];
ybs[302]=['',0.3480771,1.5206855,6.25];
ybs[303]=['73 Psc',0.2870392,0.1003544,6];
ybs[304]=['72 Psc',0.2880507,0.2624902,5.68];
ybs[305]=['',0.2945039,1.0970261,6.54];
ybs[306]=['ψ1 Psc',0.2906947,0.3764105,5.34];
ybs[307]=['ψ1 Psc',0.2907528,0.3762699,5.56];
ybs[308]=['',0.308794,1.3980888,6.29];
ybs[309]=['77 Psc',0.2911494,0.0872964,6.35];
ybs[310]=['77 Psc',0.2913095,0.0873158,7.25];
ybs[311]=['27 Cet',0.2901211,-0.172539,6.12];
ybs[312]=['',0.2970307,0.9953303,6.43];
ybs[313]=['28 Cet',0.2921797,-0.1701014,5.58];
ybs[314]=['',0.2976169,0.9353489,6.38];
ybs[315]=['75 Psc',0.2944495,0.227755,6.12];
ybs[316]=['',0.2922366,-0.4171189,6.14];
ybs[317]=['μ Cas',0.3025294,0.9601641,5.17];
ybs[318]=['β Phe',0.2917409,-0.8137642,3.31];
ybs[319]=['',0.2934697,-0.6207704,6.61];
ybs[320]=['41 And',0.3011684,0.7685566,5.03];
ybs[321]=['',0.2969891,-0.4171891,6.37];
ybs[322]=['',0.303859,1.0185155,5.79];
ybs[323]=['78 Psc',0.3010276,0.5603437,6.25];
ybs[324]=['ψ2 Psc',0.3006058,0.3635919,5.55];
ybs[325]=['30 Cet',0.2995303,-0.1691645,5.82];
ybs[326]=['80 Psc',0.3022841,0.1002306,5.52];
ybs[327]=['υ Phe',0.2992994,-0.7224581,5.21];
ybs[328]=['ι Tuc',0.2966962,-1.0765551,5.37];
ybs[329]=['',0.3222258,1.3921865,5.64];
ybs[330]=['η Cet',0.3031036,-0.1760894,3.45];
ybs[331]=['φ And',0.307733,0.8261492,4.25];
ybs[332]=['31 Cas',0.3135471,1.2020321,5.29];
ybs[333]=['β And',0.3085434,0.6233172,2.06];
ybs[334]=['ζ Phe',0.3015764,-0.9625972,3.92];
ybs[335]=['ψ3 Psc',0.3087446,0.3447286,5.55];
ybs[336]=['44 And',0.311181,0.7360787,5.65];
ybs[337]=['',0.3110035,0.4459419,5.8];
ybs[338]=['',0.3166605,1.1221671,5.55];
ybs[339]=['θ Cas',0.3149126,0.9641622,4.33];
ybs[340]=['',0.3103276,0.2751861,6.06];
ybs[341]=['32 Cas',0.3178544,1.1364103,5.57];
ybs[342]=['32 Cet',0.3101384,-0.1538207,6.4];
ybs[343]=['33 Cet',0.3118122,0.0443025,5.95];
ybs[344]=['45 And',0.3148638,0.6600289,5.81];
ybs[345]=['82 Psc',0.3145239,0.550083,5.16];
ybs[346]=['',0.3090717,-1.0053373,6.41];
ybs[347]=['χ Psc',0.3158899,0.3687426,4.66];
ybs[348]=['τ Psc',0.3168948,0.5267816,4.51];
ybs[349]=['34 Cet',0.3168526,-0.0376725,5.94];
ybs[350]=['',0.3241567,1.0785832,6.41];
ybs[351]=['',0.3210856,0.7929082,6.11];
ybs[352]=['',0.3227105,0.5263325,6.19];
ybs[353]=['',0.3409366,1.3962972,6.26];
ybs[354]=['',0.3194585,-0.5359848,6.52];
ybs[355]=['',0.3209602,-0.6591041,5.92];
ybs[356]=['φ Psc',0.3259494,0.430677,4.65];
ybs[357]=['ζ Psc',0.325701,0.1338256,5.24];
ybs[358]=['ζ Psc',0.3258028,0.133879,6.3];
ybs[359]=['',0.3274591,0.4995488,6.43];
ybs[360]=['87 Psc',0.3275099,0.2831958,5.98];
ybs[361]=['',0.3381068,1.2537728,7.83];
ybs[362]=['37 Cet',0.3284692,-0.1366728,5.13];
ybs[363]=['88 Psc',0.3299505,0.1237004,6.03];
ybs[364]=['38 Cet',0.3303684,-0.015388,5.7];
ybs[365]=['',0.3379326,0.8407986,6.61];
ybs[366]=['ν Phe',0.3314099,-0.7930638,4.96];
ybs[367]=['',0.3372586,0.5795667,6.02];
ybs[368]=['',0.3408241,0.7852907,6.34];
ybs[369]=['39 Cet',0.3381419,-0.0420329,5.41];
ybs[370]=['',0.3419933,0.5556531,6.73];
ybs[371]=['',0.3570805,1.3554567,6.31];
ybs[372]=['',0.3456185,0.8292315,6.25];
ybs[373]=['κ Tuc',0.3330847,-1.2005072,4.86];
ybs[374]=['89 Psc',0.3434164,0.0646863,5.16];
ybs[375]=['',0.3481119,0.6541104,6.46];
ybs[376]=['',0.3388651,-1.1572602,6.24];
ybs[377]=['',0.3641009,1.3322105,6.38];
ybs[378]=['φ Cas',0.3542837,1.0179306,4.98];
ybs[379]=['υ Psc',0.3509523,0.4774475,4.76];
ybs[380]=['35 Cas',0.3589761,1.1300944,6.34];
ybs[381]=['42 Cet',0.3521238,-0.0072846,5.87];
ybs[382]=['',0.3725626,1.3756109,6.07];
ybs[383]=['',0.3554554,-0.0550746,6.23];
ybs[384]=['',0.3548847,-0.19456,6.15];
ybs[385]=['91 Psc',0.3582065,0.5031675,5.23];
ybs[386]=['ξ And',0.3638032,0.7962195,4.88];
ybs[387]=['',0.3686216,1.0163755,6.45];
ybs[388]=['',0.3644154,0.0317212,6.2];
ybs[389]=['43 Cet',0.3642325,-0.006259,6.49];
ybs[390]=['',0.3637163,-0.3314428,6.35];
ybs[391]=['47 And',0.3694948,0.6598379,5.58];
ybs[392]=['',0.3692125,0.5992897,6.29];
ybs[393]=['',0.3681078,0.3588373,5.97];
ybs[394]=['',0.3799987,1.2404156,6.49];
ybs[395]=['ψ Cas',0.3804343,1.1906734,4.74];
ybs[396]=['',0.3679593,-0.5385141,5.84];
ybs[397]=['44 Cet',0.3705287,-0.138171,6.21];
ybs[398]=['θ Cet',0.3704469,-0.1412398,3.6];
ybs[399]=['δ Cas',0.3794558,1.052885,2.68];
ybs[400]=['',0.3718566,-0.1190993,5.91];
ybs[401]=['',0.3731546,-0.2717388,6.14];
ybs[402]=['',0.3739439,-0.0481335,6.15];
ybs[403]=['',0.3776641,0.4119379,6.18];
ybs[404]=['',0.3728551,-0.7225961,5.42];
ybs[405]=['',0.3811088,0.7600612,5.96];
ybs[406]=['',0.3802254,0.6051104,6.31];
ybs[407]=['',0.3728739,-0.7755814,6.26];
ybs[408]=['46 Cet',0.3773393,-0.2532167,4.9];
ybs[409]=['ρ Psc',0.3804909,0.3361985,5.38];
ybs[410]=['94 Psc',0.3824127,0.3373851,5.5];
ybs[411]=['',0.3844082,0.6015783,6.27];
ybs[412]=['',0.3811407,-0.0053774,6.41];
ybs[413]=['ω And',0.3870363,0.7940719,4.83];
ybs[414]=['',0.3860133,0.7189167,6.46];
ybs[415]=['',0.3830895,0.0632806,6.58];
ybs[416]=['',0.3738938,-1.1218749,5.93];
ybs[417]=['47 Cet',0.3827656,-0.2263032,5.66];
ybs[418]=['',0.387484,0.7055641,6.6];
ybs[419]=['',0.3829814,-0.5664051,5.79];
ybs[420]=['α UMi',0.7571964,1.5591161,2.02];
ybs[421]=['',0.3867891,-0.188694,6.13];
ybs[422]=['',0.3896417,0.1405267,6.2];
ybs[423]=['38 Cas',0.403834,1.2279164,5.81];
ybs[424]=['',0.4019104,1.1551954,6.14];
ybs[425]=['γ Phe',0.388874,-0.754473,3.41];
ybs[426]=['49 And',0.3977673,0.822,5.27];
ybs[427]=['',0.3906045,-0.5877127,6.58];
ybs[428]=['97 Psc',0.3963136,0.3219349,6.02];
ybs[429]=['48 Cet',0.394617,-0.3759341,5.12];
ybs[430]=['μ Psc',0.3974908,0.1088003,4.84];
ybs[431]=['',0.3937701,-0.8144815,6.31];
ybs[432]=['',0.3979556,-0.4558433,5.93];
ybs[433]=['η Psc',0.4032676,0.2694007,3.62];
ybs[434]=['',0.4063574,0.608938,6.39];
ybs[435]=['',0.4126798,1.0195663,5.7];
ybs[436]=['δ Phe',0.401312,-0.8549151,3.95];
ybs[437]=['',0.403737,-0.526979,5.82];
ybs[438]=['χ Cas',0.4149145,1.0353504,4.71];
ybs[439]=['',0.4031423,-0.7938783,6.17];
ybs[440]=['',0.4098506,-0.1557757,6.59];
ybs[441]=['',0.4089075,-0.6418591,5.51];
ybs[442]=['',0.4157934,0.6514692,5.88];
ybs[443]=['',0.4071883,-0.866351,6.28];
ybs[444]=['',0.4127392,-0.1210552,5.76];
ybs[445]=['',0.4313111,1.2983412,6.58];
ybs[446]=['',0.4178664,0.3237531,5.89];
ybs[447]=['49 Cet',0.4166231,-0.2720435,5.63];
ybs[448]=['',0.4228532,0.7184704,6.38];
ybs[449]=['',0.4173265,-0.5550686,6.12];
ybs[450]=['',0.4255674,0.8519233,5.92];
ybs[451]=['101 Psc',0.4219855,0.2574422,6.22];
ybs[452]=['40 Cas',0.4361123,1.276332,5.28];
ybs[453]=['',0.422633,0.3058262,5.8];
ybs[454]=['υ And',0.4268873,0.7242127,4.09];
ybs[455]=['50 Cet',0.4225291,-0.2672334,5.42];
ybs[456]=['',0.4184283,-1.01317,6.01];
ybs[457]=['',0.433209,1.0134435,5.56];
ybs[458]=['τ Scl',0.4230042,-0.5204325,5.69];
ybs[459]=['π Psc',0.4277366,0.2134608,5.57];
ybs[460]=['51 And',0.4322961,0.8502702,3.57];
ybs[461]=['',0.434548,0.7939236,6.36];
ybs[462]=['',0.4297812,-0.1625817,6.24];
ybs[463]=['',0.4092417,-1.3686057,6.11];
ybs[464]=['',0.4249222,-1.0154677,6.18];
ybs[465]=['χ And',0.4381134,0.7762253,4.98];
ybs[466]=['',0.4421917,0.9417185,6.39];
ybs[467]=['',0.4329861,-0.635995,5.94];
ybs[468]=['α Eri',0.4291813,-0.9974211,0.46];
ybs[469]=['',0.4350089,-0.3697802,5.58];
ybs[470]=['',0.4348229,-0.4351718,6.7];
ybs[471]=['105 Psc',0.439062,0.2878764,5.97];
ybs[472]=['',0.443826,0.7572262,5.61];
ybs[473]=['τ And',0.4433992,0.709739,4.94];
ybs[474]=['43 Cas',0.4522998,1.189107,5.59];
ybs[475]=['',0.4340647,-0.9311408,6.84];
ybs[476]=['42 Cas',0.455155,1.2341248,5.18];
ybs[477]=['',0.4505236,1.0668526,6.71];
ybs[478]=['',0.4514793,1.0247798,6.37];
ybs[479]=['',0.4487181,0.7452816,4.95];
ybs[480]=['',0.4463003,0.450885,6.17];
ybs[481]=['',0.4478846,0.5259573,5.99];
ybs[482]=['',0.4382567,-0.9793003,5.87];
ybs[483]=['',0.4382859,-0.9792422,5.76];
ybs[484]=['',0.4545852,1.07354,6.34];
ybs[485]=['ν Psc',0.4465634,0.09731,4.44];
ybs[486]=['',0.4497478,0.616684,5.64];
ybs[487]=['44 Cas',0.4560902,1.0583449,5.78];
ybs[488]=['',0.4477135,-0.1961195,5.75];
ybs[489]=['107 Psc',0.451413,0.3552857,5.24];
ybs[490]=['',0.4460215,-0.6640121,6.17];
ybs[491]=['',0.4552894,0.792551,6.34];
ybs[492]=['φ Per',0.4571291,0.8862109,4.07];
ybs[493]=['π Scl',0.4491331,-0.5626785,5.25];
ybs[494]=['',0.448641,-0.641315,5.72];
ybs[495]=['',0.4601838,1.0057249,6.21];
ybs[496]=['',0.4520888,-0.0628765,4.99];
ybs[497]=['',0.4467255,-0.8718087,6.64];
ybs[498]=['',0.4622192,0.9979178,6.25];
ybs[499]=['',0.4574434,0.5633779,6.34];
ybs[500]=['',0.4604099,0.8068152,6.35];
ybs[501]=['',0.4467917,-1.0594415,5.71];
ybs[502]=['',0.4500993,-0.9364174,5.52];
ybs[503]=['',0.4572592,-0.0816474,6.19];
ybs[504]=['109 Psc',0.4620336,0.3520391,6.27];
ybs[505]=['τ Cet',0.4577841,-0.2766351,3.5];
ybs[506]=['ο Psc',0.4639015,0.1613557,4.26];
ybs[507]=['',0.4756276,1.1159453,5.63];
ybs[508]=['',0.425528,-1.4466387,5.87];
ybs[509]=['',0.4662989,-0.0985451,5.34];
ybs[510]=['ε Scl',0.4645299,-0.435727,5.31];
ybs[511]=['',0.4692336,0.3054334,6.55];
ybs[512]=['τ1 Hyi',0.442377,-1.3798622,6.33];
ybs[513]=['',0.4661054,-0.4758126,6.39];
ybs[514]=['',0.4750891,0.8083749,6.32];
ybs[515]=['',0.4659145,-0.885393,5.49];
ybs[516]=['',0.4658554,-0.9326139,5.04];
ybs[517]=['',0.4785904,0.6639121,5.94];
ybs[518]=['4 Ari',0.4761824,0.2974432,5.84];
ybs[519]=['',0.4786789,0.572064,5.79];
ybs[520]=['',0.4713247,-0.7273332,6.18];
ybs[521]=['',0.4216237,-1.4779605,5.69];
ybs[522]=['',0.4815263,0.8374684,5.82];
ybs[523]=['',0.4770939,0.0658372,5.91];
ybs[524]=['',0.4736801,-0.647045,6.32];
ybs[525]=['',0.4890484,0.9079108,5.9];
ybs[526]=['1 Ari',0.4848212,0.3902831,5.86];
ybs[527]=['χ Cet',0.4819207,-0.1850044,4.67];
ybs[528]=['',0.4804626,-0.540813,6.34];
ybs[529]=['1 Per',0.4937084,0.9640048,5.52];
ybs[530]=['',0.4878162,0.1942462,5.94];
ybs[531]=['',0.4824422,-0.6687666,6.37];
ybs[532]=['2 Per',0.4942752,0.8880001,5.79];
ybs[533]=['',0.4844779,-0.8330474,6.14];
ybs[534]=['',0.4973159,0.8998999,6.26];
ybs[535]=['ζ Cet',0.4901047,-0.178878,3.73];
ybs[536]=['',0.5016866,0.9718621,6.45];
ybs[537]=['',0.4868914,-0.8747579,5.94];
ybs[538]=['ε Cas',0.5047121,1.1127418,3.38];
ybs[539]=['55 And',0.4989258,0.7123628,5.4];
ybs[540]=['α Tri',0.4977852,0.5177449,3.41];
ybs[541]=['γ1 Ari',0.4995664,0.3382701,4.83];
ybs[542]=['γ2 Ari',0.4995664,0.3382313,4.75];
ybs[543]=['',0.4961496,-0.2939728,5.8];
ybs[544]=['ω Cas',0.512188,1.2002689,4.99];
ybs[545]=['ξ Psc',0.4994365,0.0571267,4.62];
ybs[546]=['τ2 Hyi',0.469743,-1.3978307,6.06];
ybs[547]=['',0.5059518,0.7118722,6.24];
ybs[548]=['',0.5061364,0.6495008,6.26];
ybs[549]=['β Ari',0.5044367,0.3646595,2.64];
ybs[550]=['',0.4980108,-0.6721098,6.1];
ybs[551]=['ψ Phe',0.4989428,-0.8066367,4.41];
ybs[552]=['',0.5102823,0.6521057,5.89];
ybs[553]=['56 And',0.5113665,0.651649,5.67];
ybs[554]=['φ Phe',0.5021856,-0.7402198,5.11];
ybs[555]=['7 Ari',0.509766,0.4129862,5.74];
ybs[556]=['',0.5096354,0.0337697,6.01];
ybs[557]=['',0.5227792,1.0783098,6.02];
ybs[558]=['',0.5192735,0.7291836,6.78];
ybs[559]=['ι Ari',0.5162262,0.3124546,5.1];
ybs[560]=['',0.5180506,0.4867583,5.82];
ybs[561]=['56 Cet',0.5126383,-0.3916861,4.85];
ybs[562]=['χ Eri',0.508836,-0.8992588,3.7];
ybs[563]=['',0.527712,1.1293276,5.26];
ybs[564]=['3 Per',0.5222092,0.8602506,5.69];
ybs[565]=['λ Ari',0.5188406,0.4133083,4.79];
ybs[566]=['η2 Hyi',0.5034247,-1.1791765,4.69];
ybs[567]=['',0.5075932,-1.0607445,6.06];
ybs[568]=['',0.5445003,1.3613556,6.04];
ybs[569]=['',0.513374,-0.9020067,6.1];
ybs[570]=['',0.5142442,-0.8255425,4.83];
ybs[571]=['48 Cas',0.5385618,1.2390221,4.54];
ybs[572]=['',0.5201795,-0.5756453,6.35];
ybs[573]=['',0.5260734,0.3690096,5.87];
ybs[574]=['',0.5252183,0.2160564,6.09];
ybs[575]=['',0.5443979,1.290393,6.23];
ybs[576]=['50 Cas',0.5452758,1.2654485,3.98];
ybs[577]=['47 Cas',0.5538176,1.3502645,5.38];
ybs[578]=['112 Psc',0.5282282,0.0555271,5.88];
ybs[579]=['57 Cet',0.5261783,-0.3619831,5.41];
ybs[580]=['',0.5164909,-1.1403973,6.37];
ybs[581]=['υ Cet',0.5272061,-0.3664055,4];
ybs[582]=['52 Cas',0.5419384,1.1342022,6];
ybs[583]=['',0.5293409,-0.1472958,5.51];
ybs[584]=['',0.5252053,-0.732099,5.57];
ybs[585]=['53 Cas',0.5424501,1.1252764,5.58];
ybs[586]=['4 Per',0.5387744,0.9524481,5.04];
ybs[587]=['α Hyi',0.5205891,-1.0731182,2.86];
ybs[588]=['49 Cas',0.555213,1.3299058,5.22];
ybs[589]=['σ Hyi',0.5053905,-1.3659486,6.16];
ybs[590]=['π For',0.5324485,-0.5221614,5.35];
ybs[591]=['α Psc',0.536485,0.0496974,3.82]; //manual change !!!!
ybs[592]=['α Psc',0.536485,0.0496974,4.33];
ybs[593]=['',0.5746416,1.4203119,6.05];
ybs[594]=['',0.5497906,1.1377199,6.52];
ybs[595]=['ε Tri',0.5410305,0.5823729,5.5];
ybs[596]=['',0.5242272,-1.1516031,6.1];
ybs[597]=['',0.5390078,0.2366733,5.94];
ybs[598]=['χ Phe',0.5341054,-0.7789348,5.14];
ybs[599]=['γ1 And',0.5453307,0.7402489,2.26];
ybs[600]=['γ2 And',0.5453817,0.7402682,4.84];
ybs[601]=['10 Ari',0.5438874,0.4541178,5.63];
ybs[602]=['',0.5377886,-0.51629,6.42];
ybs[603]=['60 Cet',0.5414557,0.0036987,5.43];
ybs[604]=['',0.5402625,-0.2656773,5.86];
ybs[605]=['',0.5439902,0.3200376,6.21];
ybs[606]=['61 Cet',0.5441028,-0.0044824,5.93];
ybs[607]=['',0.5434853,-0.0701644,5.62];
ybs[608]=['ν For',0.5466099,-0.5098739,4.69];
ybs[609]=['κ Ari',0.5565318,0.3967338,5.03];
ybs[610]=['',0.5547168,0.1453933,6.31];
ybs[611]=['11 Ari',0.5577026,0.4500767,6.15];
ybs[612]=['',0.5558242,0.0020571,6.28];
ybs[613]=['α Ari',0.5592028,0.4109413,2];
ybs[614]=['',0.5668458,1.0211213,5.67];
ybs[615]=['',0.5657547,0.7774014,6.42];
ybs[616]=['58 And',0.5652561,0.6622051,4.82];
ybs[617]=['',0.5729711,0.9411698,6.31];
ybs[618]=['β Tri',0.5697936,0.6120761,3];
ybs[619]=['14 Ari',0.5690671,0.4541681,4.98];
ybs[620]=['',0.5687458,0.302058,6.43];
ybs[621]=['',0.5654517,-0.3088724,6.1];
ybs[622]=['',0.589069,1.2934457,6.29];
ybs[623]=['5 Per',0.579083,1.0075353,6.36];
ybs[624]=['59 And',0.5757344,0.6827953,5.63];
ybs[625]=['59 And',0.5758,0.6828583,6.1];
ybs[626]=['',0.5689003,-0.4234808,6.48];
ybs[627]=['15 Ari',0.5742067,0.3417738,5.7];
ybs[628]=['',0.5665952,-0.758073,5.85];
ybs[629]=['16 Ari',0.5768319,0.4541124,6.02];
ybs[630]=['5 Tri',0.5778957,0.5516657,6.23];
ybs[631]=['64 Cet',0.5771834,0.1509968,5.63];
ybs[632]=['',0.5706239,-0.7632933,6.32];
ybs[633]=['',0.5718819,-0.8856227,6.12];
ybs[634]=['',0.57696,-0.1740175,6.01];
ybs[635]=['63 Cet',0.5780834,-0.0304312,5.93];
ybs[636]=['55 Cas',0.5928971,1.1624837,6.07];
ybs[637]=['',0.5887943,1.0235009,6.44];
ybs[638]=['6 Tri',0.5820393,0.5303106,4.94];
ybs[639]=['60 And',0.5861163,0.7734071,4.83];
ybs[640]=['',0.5830187,0.4232289,5.96];
ybs[641]=['',0.5880402,0.8926841,5.31];
ybs[642]=['η Ari',0.5837334,0.3716198,5.27];
ybs[643]=['',0.5898171,0.8301706,6.06];
ybs[644]=['19 Ari',0.5847318,0.2681014,5.71];
ybs[645]=['ξ1 Cet',0.5843818,0.1558237,4.37];
ybs[646]=['66 Cet',0.5832881,-0.0403553,5.54];
ybs[647]=['',0.5839428,-0.3651037,5.86];
ybs[648]=['μ For',0.5832812,-0.5348121,5.28];
ybs[649]=['',0.5982269,0.8358737,6.33];
ybs[650]=['',0.6025688,0.9972061,6.48];
ybs[651]=['7 Tri',0.5976779,0.5836302,5.28];
ybs[652]=['20 Ari',0.5967683,0.4514078,5.79];
ybs[653]=['21 Ari',0.5965203,0.4384926,5.58];
ybs[654]=['',0.5948625,-0.163795,6.55];
ybs[655]=['',0.5901016,-0.7170796,5.91];
ybs[656]=['δ Tri',0.6025759,0.5987276,4.87];
ybs[657]=['8 Per',0.6075886,1.0119396,5.75];
ybs[658]=['7 Per',0.6079012,1.0052536,5.98];
ybs[659]=['',0.6050598,0.774703,6.7];
ybs[660]=['γ Tri',0.6037094,0.5921476,4.01];
ybs[661]=['',0.6028692,0.4162288,6.55];
ybs[662]=['67 Cet',0.6015087,-0.1106851,5.51];
ybs[663]=['π1 Hyi',0.5873614,-1.1826443,5.55];
ybs[664]=['',0.6177619,1.1242855,6.6];
ybs[665]=['θ Ari',0.6069483,0.3487387,5.62];
ybs[666]=['62 And',0.6126975,0.8283306,5.3];
ybs[667]=['',0.6122347,0.8124921,6.21];
ybs[668]=['',0.6061842,0.0320784,5.58];
ybs[669]=['',0.6131877,0.8558239,6.37];
ybs[670]=['φ Eri',0.5983502,-0.8976516,3.56];
ybs[671]=['10 Tri',0.610729,0.5013011,5.03];
ybs[672]=['',0.6106868,0.4057492,6.46];
ybs[673]=['',0.6139435,0.6966441,6.63];
ybs[674]=['π2 Hyi',0.5927207,-1.1809864,5.69];
ybs[675]=['',0.6188565,0.8271173,6.11];
ybs[676]=['',0.6156691,0.5282762,6.47];
ybs[677]=['ο Cet',0.6118708,-0.0505734,3.04];
ybs[678]=['63 And',0.6201928,0.8766932,5.59];
ybs[679]=['',0.6098345,-0.4514398,6.34];
ybs[680]=['',0.613309,-0.0744519,6.5];
ybs[681]=['9 Per',0.6265279,0.976069,5.17];
ybs[682]=['',0.6113406,-0.7289973,6.37];
ybs[683]=['',0.6280423,0.7238817,5.82];
ybs[684]=['',0.6129217,-0.9750275,5.81];
ybs[685]=['69 Cet',0.6232656,0.0082913,5.28];
ybs[686]=['',0.6330889,0.9676654,6.28];
ybs[687]=['70 Cet',0.6243925,-0.0140647,5.42];
ybs[688]=['',0.6234194,-0.1867253,5.46];
ybs[689]=['',0.6235544,-0.3068817,5.87];
ybs[690]=['64 And',0.6352431,0.8741523,5.19];
ybs[691]=['κ For',0.6254386,-0.4142941,5.2];
ybs[692]=['10 Per',0.6392948,0.9893983,6.25];
ybs[693]=['',0.6273791,-0.3189671,6.22];
ybs[694]=['',0.6234345,-0.7526002,6.31];
ybs[695]=['65 And',0.6405358,0.8788933,4.71];
ybs[696]=['',0.6275759,-0.6544537,6.53];
ybs[697]=['',0.6262452,-0.8903484,5.92];
ybs[698]=['ξ Ari',0.6359848,0.1865592,5.47];
ybs[699]=['',0.6332082,-0.4497515,6.44];
ybs[700]=['71 Cet',0.6364308,-0.0471507,6.33];
ybs[701]=['δ Hyi',0.6198865,-1.1969486,4.09];
ybs[702]=['',0.6338225,-0.7114305,6.18];
ybs[703]=['ι Cas',0.6568285,1.1777455,4.52];
ybs[704]=['ρ Cet',0.6405204,-0.2131454,4.89];
ybs[705]=['66 And',0.6503385,0.8839642,6.12];
ybs[706]=['',0.6407031,-0.2663879,5.83];
ybs[707]=['',0.6463526,0.4728313,6.18];
ybs[708]=['11 Tri',0.6479718,0.556397,5.54];
ybs[709]=['',0.6431481,-0.34845,5.88];
ybs[710]=['λ Hor',0.6343893,-1.0512713,5.35];
ybs[711]=['κ Hyi',0.6239331,-1.283982,5.01];
ybs[712]=['',0.6573774,0.9706417,6.51];
ybs[713]=['12 Tri',0.6509987,0.5191845,5.29];
ybs[714]=['ξ2 Cet',0.6505297,0.1490101,4.28];
ybs[715]=['',0.6497211,0.035579,6.45];
ybs[716]=['13 Tri',0.653807,0.5237631,5.89];
ybs[717]=['κ Eri',0.6441379,-0.8312288,4.25];
ybs[718]=['',0.6361918,-1.1591832,6.41];
ybs[719]=['',0.6554873,0.4109597,6.19];
ybs[720]=['φ For',0.6491231,-0.5887589,5.14];
ybs[721]=['',0.6567909,0.1683041,6.07];
ybs[722]=['',0.6603173,0.5918582,6.25];
ybs[723]=['',0.6516398,-0.5414873,6.11];
ybs[724]=['',0.6612598,0.4417781,5.92];
ybs[725]=['26 Ari',0.6615815,0.3478839,6.15];
ybs[726]=['',0.6576448,-0.3945416,6.77];
ybs[727]=['27 Ari',0.6627018,0.3103339,6.23];
ybs[728]=['',0.6617044,0.005799,6];
ybs[729]=['',0.6589278,-0.4382392,6.51];
ybs[730]=['',0.6478792,-1.1208851,6.37];
ybs[731]=['',0.6603661,-0.3921494,6.1];
ybs[732]=['14 Tri',0.668358,0.632225,5.15];
ybs[733]=['',0.6650063,0.0409107,5.25];
ybs[734]=['',0.6716839,0.6042139,5.83];
ybs[735]=['75 Cet',0.6678039,-0.016727,5.35];
ybs[736]=['σ Cet',0.6672315,-0.2647329,4.75];
ybs[737]=['29 Ari',0.671357,0.263739,6.04];
ybs[738]=['',0.6674408,-0.6344426,6.3];
ybs[739]=['',0.6969866,1.2722276,5.16];
ybs[740]=['λ1 For',0.6712825,-0.6034234,5.9];
ybs[741]=['',0.6740407,-0.3477694,6.21];
ybs[742]=['',0.6831461,0.6935966,6.36];
ybs[743]=['',0.6939895,1.148787,5.78];
ybs[744]=['',0.6838677,0.6525417,5.71];
ybs[745]=['ω For',0.6746207,-0.4914204,4.9];
ybs[746]=['15 Tri',0.6843749,0.6067311,5.35];
ybs[747]=['',0.6806674,0.131724,6.18];
ybs[748]=['77 Cet',0.6788242,-0.1358477,5.75];
ybs[749]=['',0.6850786,0.1215188,5.82];
ybs[750]=['ν Cet',0.6841522,0.098942,4.86];
ybs[751]=['',0.6741613,-0.8904219,6.24];
ybs[752]=['',0.689618,0.677329,5.9];
ybs[753]=['',0.6883756,0.5529706,6.1];
ybs[754]=['',0.6898732,0.5993318,5.3];
ybs[755]=['80 Cet',0.6844451,-0.1353689,5.53];
ybs[756]=['',0.6913729,0.6976262,6.54];
ybs[757]=['',0.6901199,0.5753865,6.25];
ybs[758]=['',0.6719102,-1.0910164,6.77];
ybs[759]=['31 Ari',0.6875865,0.2185662,5.68];
ybs[760]=['30 Ari',0.6892838,0.4315092,7.09];
ybs[761]=['30 Ari',0.6894875,0.4314944,6.5];
ybs[762]=['',0.6872913,0.1362257,5.81];
ybs[763]=['ι1 For',0.6846484,-0.5230603,5.75];
ybs[764]=['',0.6954546,0.6597626,6.18];
ybs[765]=['',0.6961937,0.6660935,6.3];
ybs[766]=['',0.693524,0.1356177,6.39];
ybs[767]=['81 Cet',0.6919315,-0.0579618,5.65];
ybs[768]=['λ2 For',0.6881115,-0.6021906,5.79];
ybs[769]=['ν Ari',0.6973219,0.3846044,5.43];
ybs[770]=['',0.7436902,1.4227997,5.78];
ybs[771]=['',0.6960682,0.0613997,6.21];
ybs[772]=['μ Hyi',0.6603139,-1.3793781,5.28];
ybs[773]=['ι2 For',0.6940478,-0.5256787,5.83];
ybs[774]=['η Hor',0.689322,-0.9157358,5.31];
ybs[775]=['δ Cet',0.6997978,0.0070382,4.07];
ybs[776]=['',0.6942751,-0.6617518,6.49];
ybs[777]=['ε Cet',0.6999109,-0.2059067,4.84];
ybs[778]=['33 Ari',0.7055969,0.4735975,5.3];
ybs[779]=['',0.7033051,0.1079728,6.25];
ybs[780]=['',0.7027579,-0.1636873,5.78];
ybs[781]=['11 Per',0.7169379,0.9630631,5.77];
ybs[782]=['',0.7015793,-0.5333616,6.52];
ybs[783]=['',0.7166216,0.9354921,5.84];
ybs[784]=['12 Per',0.7127935,0.7028048,4.91];
ybs[785]=['',0.7001577,-0.7472988,4.75];
ybs[786]=['84 Cet',0.7074067,-0.01085,5.71];
ybs[787]=['',0.7258767,1.1850403,5.95];
ybs[788]=['',0.7165232,0.8436781,6.48];
ybs[789]=['μ Ari',0.7127765,0.3505584,5.69];
ybs[790]=['ι Eri',0.7040381,-0.6943129,4.11];
ybs[791]=['',0.7098596,-0.0547914,6.05];
ybs[792]=['',0.7085953,-0.2526426,5.98];
ybs[793]=['',0.7130848,0.1887659,6.3];
ybs[794]=['',0.697725,-1.1206276,6.55];
ybs[795]=['θ Per',0.7216718,0.8604759,4.12];
ybs[796]=['14 Per',0.7209702,0.7744076,5.43];
ybs[797]=['35 Ari',0.7176926,0.4848658,4.66];
ybs[798]=['ζ Hor',0.7033845,-0.9507789,5.21];
ybs[799]=['',0.7193958,0.4487501,6.35];
ybs[800]=['γ Cet',0.7165079,0.0577604,3.47];
ybs[801]=['',0.7103761,-0.6686344,6.01];
ybs[802]=['ε Hyi',0.6975233,-1.1901791,4.11];
ybs[803]=['',0.7102241,-0.8107139,6.1];
ybs[804]=['36 Ari',0.7212471,0.3113177,6.46];
ybs[805]=['ο Ari',0.7221988,0.2685173,5.77];
ybs[806]=['ι Hor',0.7118498,-0.8853432,5.41];
ybs[807]=['π Cet',0.7197614,-0.2405977,4.25];
ybs[808]=['38 Ari',0.723936,0.2184971,5.18];
ybs[809]=['μ Cet',0.7238085,0.1778019,4.27];
ybs[810]=['',0.7156605,-0.7060534,6.36];
ybs[811]=['',0.7440845,1.2166006,6.18];
ybs[812]=['',0.7254752,0.0835086,6.03];
ybs[813]=['',0.7202764,-0.5663884,6.22];
ybs[814]=['τ1 Eri',0.7239358,-0.3228753,4.47];
ybs[815]=['',0.7332901,0.6292986,6.25];
ybs[816]=['',0.7336558,0.6218175,6.3];
ybs[817]=['',0.7188176,-0.9162482,6.15];
ybs[818]=['',0.7238871,-0.8065888,6.85];
ybs[819]=['',0.7144873,-1.1631014,6.26];
ybs[820]=['39 Ari',0.7371912,0.5117218,4.51];
ybs[821]=['',0.7453335,0.9975596,6.25];
ybs[822]=['',0.7310635,-0.3764165,6.49];
ybs[823]=['',0.7329337,-0.3911812,6.47];
ybs[824]=['40 Ari',0.7396648,0.3203678,5.82];
ybs[825]=['',0.7573306,1.2035728,5.8];
ybs[826]=['',0.740829,0.4408718,5.86];
ybs[827]=['',0.744163,0.6527172,6.45];
ybs[828]=['',0.7364098,-0.2162158,6.9];
ybs[829]=['γ Hor',0.7235891,-1.110577,5.74];
ybs[830]=['η Per',0.7504278,0.9768085,3.76];
ybs[831]=['η1 For',0.7342275,-0.619215,6.51];
ybs[832]=['π Ari',0.7429496,0.306062,5.22];
ybs[833]=['ζ Hyi',0.7235059,-1.1788584,4.84];
ybs[834]=['41 Ari',0.7462,0.4770377,3.63];
ybs[835]=['',0.7552357,1.0190254,6.45];
ybs[836]=['16 Per',0.749145,0.6700339,4.23];
ybs[837]=['β For',0.7409815,-0.5643321,4.46];
ybs[838]=['',0.7543143,0.8187884,5.88];
ybs[839]=['17 Per',0.7530967,0.6131509,4.53];
ybs[840]=['γ1 For',0.7444992,-0.4274052,6.14];
ybs[841]=['γ2 For',0.7446474,-0.4864218,5.39];
ybs[842]=['',0.7597382,0.9262217,6.36];
ybs[843]=['σ Ari',0.7525076,0.2644733,5.49];
ybs[844]=['η2 For',0.7459236,-0.6243383,5.92];
ybs[845]=['',0.7616377,0.8489304,6.26];
ybs[846]=['τ2 Eri',0.7497657,-0.3653454,4.75];
ybs[847]=['η3 For',0.7477887,-0.6214171,5.47];
ybs[848]=['ν Hor',0.7392,-1.0949255,5.26];
ybs[849]=['',0.7481962,-0.6956911,6.36];
ybs[850]=['τ Per',0.7658071,0.9221082,3.95];
ybs[851]=['20 Per',0.7627944,0.6703477,5.33];
ybs[852]=['',0.7599658,0.2889235,6.31];
ybs[853]=['',0.7564763,-0.22163,6.04];
ybs[854]=['',0.7533716,-0.5365713,6.4];
ybs[855]=['',0.7578843,-0.1635413,6.32];
ybs[856]=['',0.7738199,1.0749658,5.59];
ybs[857]=['',0.7761364,1.1240312,6.24];
ybs[858]=['',0.7608545,-0.3893081,5.95];
ybs[859]=['ψ For',0.7603426,-0.6696173,5.92];
ybs[860]=['',0.777022,0.8958859,6.22];
ybs[861]=['',0.7755695,0.8243823,6.02];
ybs[862]=['',0.7535322,-1.0967408,6.03];
ybs[863]=['ρ2 Ari',0.7714148,0.3211695,5.91];
ybs[864]=['',0.7612265,-0.869517,4];
ybs[865]=['ρ3 Ari',0.7741431,0.3157801,5.63];
ybs[866]=['',0.7730372,0.1475073,5.97];
ybs[867]=['',0.7621648,-0.8866417,6.21];
ybs[868]=['ν Hyi',0.7434455,-1.3089138,4.75];
ybs[869]=['21 Per',0.7782231,0.5585702,5.11];
ybs[870]=['η Eri',0.7735424,-0.1540816,3.89];
ybs[871]=['',0.7745019,-0.0635728,5.17];
ybs[872]=['',0.781708,0.6751686,6.04];
ybs[873]=['',0.7766506,0.0797745,6.11];
ybs[874]=['47 Ari',0.7814148,0.3619451,5.8];
ybs[875]=['π Per',0.7849036,0.693452,4.7];
ybs[876]=['',0.7622027,-1.1233818,6.56];
ybs[877]=['',0.8227196,1.3872796,5.49];
ybs[878]=['24 Per',0.7860647,0.6152647,4.93];
ybs[879]=['4 Eri',0.7774194,-0.4152555,5.45];
ybs[880]=['',0.7764925,-0.5198579,6.29];
ybs[881]=['',0.7898806,0.8253591,5.47];
ybs[882]=['',0.7888988,0.7173631,5.89];
ybs[883]=['ε Ari',0.7863419,0.373662,4.63];
ybs[884]=['ε Ari',0.7863419,0.373662,4.63];
ybs[885]=['6 Eri',0.7804784,-0.4107941,5.84];
ybs[886]=['',0.7946785,0.9149035,5.28];
ybs[887]=['',0.7947659,0.914913,6.74];
ybs[888]=['',0.7835818,-0.0473569,5.23];
ybs[889]=['',0.7776638,-0.6653472,6.41];
ybs[890]=['',0.7911179,0.6667215,6.11];
ybs[891]=['',0.7838263,-0.1694237,6.14];
ybs[892]=['λ Cet',0.7882575,0.1566666,4.7];
ybs[893]=['θ1 Eri',0.7807145,-0.7022404,3.24];
ybs[894]=['θ2 Eri',0.780758,-0.7022356,4.35];
ybs[895]=['5 Eri',0.7878934,-0.0418207,5.56];
ybs[896]=['',0.7847668,-0.5033163,6.14];
ybs[897]=['ζ For',0.7870068,-0.4399149,5.71];
ybs[898]=['',0.7927516,0.1909181,5.95];
ybs[899]=['',0.7869684,-0.5661557,6.31];
ybs[900]=['7 Eri',0.7929604,-0.0490456,6.11];
ybs[901]=['49 Ari',0.7982149,0.4630425,5.9];
ybs[902]=['',0.8489168,1.423061,5.95];
ybs[903]=['ρ1 Eri',0.79424,-0.1325467,5.75];
ybs[904]=['',0.7975925,0.0943228,6.25];
ybs[905]=['β Hor',0.781603,-1.1170489,4.99];
ybs[906]=['93 Cet',0.799768,0.0771577,5.61];
ybs[907]=['α Cet',0.7993477,0.072567,2.53];
ybs[908]=['',0.7975421,-0.172669,5.83];
ybs[909]=['',0.7985769,-0.1121656,6.19];
ybs[910]=['ε For',0.7957702,-0.4891003,5.89];
ybs[911]=['γ Per',0.8119066,0.9350362,2.93];
ybs[912]=['',0.8052639,0.4945808,6.36];
ybs[913]=['ρ2 Eri',0.8009508,-0.1329476,5.32];
ybs[914]=['',0.8153667,0.9908729,4.76];
ybs[915]=['τ3 Eri',0.7992159,-0.4111366,4.09];
ybs[916]=['',0.8158681,0.9797507,6.11];
ybs[917]=['ρ Per',0.8128924,0.6790626,3.39];
ybs[918]=['',0.8237973,1.1191785,5.89];
ybs[919]=['',0.8136983,0.7094692,6.05];
ybs[920]=['',0.8100746,0.2779163,6.49];
ybs[921]=['ρ3 Eri',0.8077946,-0.1314821,5.26];
ybs[922]=['',0.8095764,0.0337015,6.05];
ybs[923]=['52 Ari',0.8136527,0.4419585,6.8];
ybs[924]=['52 Ari',0.8136527,0.4419585,7];
ybs[925]=['',0.8007987,-0.8186831,5.82];
ybs[926]=['',0.8261124,0.9124505,6.31];
ybs[927]=['',0.8175059,0.2313262,5.62];
ybs[928]=['',0.8459526,1.2995463,4.87];
ybs[929]=['',0.8246538,0.8268438,6.41];
ybs[930]=['μ Hor',0.8029858,-1.0414387,5.11];
ybs[931]=['',0.8177964,-0.1051011,5.27];
ybs[932]=['β Per',0.8260366,0.715965,2.12];
ybs[933]=['ι Per',0.8303461,0.8670665,4.05];
ybs[934]=['53 Ari',0.8221258,0.313225,6.11];
ybs[935]=['θ Hyi',0.79545,-1.253745,5.53];
ybs[936]=['54 Ari',0.8261855,0.3291897,6.27];
ybs[937]=['κ Per',0.831998,0.7840545,3.8];
ybs[938]=['',0.8272194,0.1489975,6.28];
ybs[939]=['',0.8229003,-0.484586,6.19];
ybs[940]=['55 Ari',0.831949,0.5086364,5.72];
ybs[941]=['',0.8342491,0.4866957,6.42];
ybs[942]=['',0.8355477,0.470574,6.02];
ybs[943]=['ω Per',0.8396201,0.6924926,4.63];
ybs[944]=['',0.8360393,0.2083568,5.98];
ybs[945]=['',0.8449887,0.8341049,6.33];
ybs[946]=['',0.8435273,0.7407364,6.15];
ybs[947]=['δ Ari',0.8405113,0.3454324,4.35];
ybs[948]=['',0.8391942,0.2288654,6.12];
ybs[949]=['',0.8349667,-0.4131686,6.38];
ybs[950]=['56 Ari',0.8433581,0.476857,5.79];
ybs[951]=['',0.838592,-0.0653869,6.05];
ybs[952]=['',0.8492051,0.841973,5.9];
ybs[953]=['',0.8381682,-0.2785544,6.26];
ybs[954]=['',0.8437402,0.1173863,5.56];
ybs[955]=['',0.8201549,-1.2077517,6.15];
ybs[956]=['',0.8335062,-0.849427,6.12];
ybs[957]=['',0.8840505,1.3578115,5.45];
ybs[958]=['94 Cet',0.84502,-0.019745,5.06];
ybs[959]=['α For',0.8412839,-0.5047822,3.87];
ybs[960]=['',0.8602451,0.9984086,5.79];
ybs[961]=['',0.9458652,1.4829867,5.61];
ybs[962]=['',0.8556968,0.742951,6.07];
ybs[963]=['',0.8686519,1.1470615,6.36];
ybs[964]=['',0.8423073,-0.7741366,5.93];
ybs[965]=['',0.8616001,0.8901428,5.03];
ybs[966]=['',0.8452239,-0.6262089,6.27];
ybs[967]=['',0.8570106,0.5344307,5.52];
ybs[968]=['ζ Ari',0.854826,0.3684137,4.89];
ybs[969]=['',0.8607703,0.7925459,6.16];
ybs[970]=['',0.8480698,-0.519054,6.16];
ybs[971]=['',0.8590103,0.574566,6.31];
ybs[972]=['',0.8601543,0.6065428,6.25];
ybs[973]=['',0.8421124,-0.9993183,5.74];
ybs[974]=['',0.8624895,0.5628193,6.06];
ybs[975]=['',0.8654212,0.7076732,6.45];
ybs[976]=['',0.8541498,-0.4544167,6.25];
ybs[977]=['',0.8158296,-1.3774622,5.57];
ybs[978]=['30 Per',0.8681881,0.7694834,5.47];
ybs[979]=['',0.8590578,-0.1021862,6.17];
ybs[980]=['ζ Eri',0.858199,-0.1528191,4.8];
ybs[981]=['',0.8794254,1.146936,4.84];
ybs[982]=['',0.867879,0.686726,5.96];
ybs[983]=['29 Per',0.8721626,0.8776403,5.15];
ybs[984]=['14 Eri',0.8615141,-0.1586653,6.14];
ybs[985]=['31 Per',0.8743255,0.8754171,5.03];
ybs[986]=['',0.8591801,-0.5369289,6.65];
ybs[987]=['',0.8719214,0.5983973,4.82];
ybs[988]=['95 Cet',0.8694602,-0.0151367,5.38];
ybs[989]=['',0.8673442,-0.5014995,5.91];
ybs[990]=['15 Eri',0.868922,-0.3917978,4.88];
ybs[991]=['59 Ari',0.8769405,0.4735704,5.9];
ybs[992]=['κ1 Cet',0.8738703,0.0599165,4.83];
ybs[993]=['',0.8704125,-0.3228301,5.71];
ybs[994]=['',0.864002,-0.8323179,5.85];
ybs[995]=['',0.8787884,0.5080772,4.47];
ybs[996]=['60 Ari',0.8790687,0.4489876,6.12];
ybs[997]=['',0.8862842,0.8575264,5.93];
ybs[998]=['32 Per',0.8841078,0.7573229,4.95];
ybs[999]=['τ4 Eri',0.8739503,-0.3786513,3.69];
ybs[1000]=['',0.8741631,-0.4199335,5.61];
ybs[1001]=['τ1 Ari',0.8824384,0.3701669,5.28];
ybs[1002]=['ζ1 Ret',0.86436,-1.0910399,5.54];
ybs[1003]=['κ2 Cet',0.8815209,0.0652346,5.69];
ybs[1004]=['',0.875048,-0.7506168,4.27];
ybs[1005]=['',0.899776,1.1283021,5.23];
ybs[1006]=['ζ2 Ret',0.8663049,-1.0898401,5.24];
ybs[1007]=['',0.8921605,0.8600057,5.29];
ybs[1008]=['62 Ari',0.8868591,0.4829193,5.52];
ybs[1009]=['',0.8792109,-0.4632825,6.39];
ybs[1010]=['',0.8647471,-1.1669916,6.05];
ybs[1011]=['τ2 Ari',0.8890918,0.3630896,5.09];
ybs[1012]=['',0.8821139,-0.4114305,5.52];
ybs[1013]=['α Per',0.8970163,0.8713052,1.79];
ybs[1014]=['',0.8858612,-0.4455131,6.35];
ybs[1015]=['',0.8970714,0.5863748,5.61];
ybs[1016]=['',0.9037459,0.9421664,6.51];
ybs[1017]=['',0.8819196,-0.8327825,6.39];
ybs[1018]=['64 Ari',0.8959879,0.4325834,5.5];
ybs[1019]=['',0.8926251,0.0862755,6.38];
ybs[1020]=['',0.8907825,-0.1349622,6.2];
ybs[1021]=['ι Hyi',0.8532348,-1.3495642,5.52];
ybs[1022]=['',0.9002459,0.7211345,6.51];
ybs[1023]=['65 Ari',0.8964341,0.3641562,6.08];
ybs[1024]=['',0.8950634,0.2214918,6.04];
ybs[1025]=['',0.9041112,0.8583755,6.09];
ybs[1026]=['ο Tau',0.897792,0.1586465,3.6];
ybs[1027]=['',0.8920636,-0.5697791,6.5];
ybs[1028]=['',0.9257724,1.2552893,6.32];
ybs[1029]=['',0.9155577,1.0526982,6.49];
ybs[1030]=['',0.9132602,0.8573499,4.98];
ybs[1031]=['',0.9184783,1.0471916,4.21];
ybs[1032]=['',0.9078092,0.3284102,6.57];
ybs[1033]=['',0.9168935,0.8710556,5.58];
ybs[1034]=['ξ Tau',0.9080935,0.170918,3.74];
ybs[1035]=['',0.9087847,0.2233158,6.28];
ybs[1036]=['',0.9220822,1.0286571,4.54];
ybs[1037]=['',0.9138983,0.5910938,5.61];
ybs[1038]=['χ1 For',0.9014887,-0.6258798,6.39];
ybs[1039]=['',0.9233207,1.037164,6.13];
ybs[1040]=['34 Per',0.9190369,0.8651283,4.67];
ybs[1041]=['',0.903706,-0.4757261,5.93];
ybs[1042]=['',0.9221919,0.9688502,5.09];
ybs[1043]=['',0.9191987,0.8202537,6.24];
ybs[1044]=['66 Ari',0.9139852,0.3990491,6.03];
ybs[1045]=['',0.902425,-0.7256462,6.32];
ybs[1046]=['',0.9112903,-0.1959451,5.73];
ybs[1047]=['',0.9244136,0.8405943,5.82];
ybs[1048]=['σ Per',0.9242259,0.8387037,4.36];
ybs[1049]=['',0.8906599,-1.214111,6.15];
ybs[1050]=['χ2 For',0.9085895,-0.6217102,5.71];
ybs[1051]=['',0.9475399,1.2811441,6.57];
ybs[1052]=['',0.9282888,0.8598943,6.29];
ybs[1053]=['χ3 For',0.911352,-0.6247149,6.5];
ybs[1054]=['',0.9297338,0.8632279,6.39];
ybs[1055]=['',0.9185312,-0.117735,5.99];
ybs[1056]=['4 Tau',0.922265,0.1988872,5.14];
ybs[1057]=['',0.9181633,-0.2201805,5.59];
ybs[1058]=['',0.9310782,0.8391889,5.47];
ybs[1059]=['',0.8975087,-1.2090877,5.96];
ybs[1060]=['',0.926797,0.4822453,5.96];
ybs[1061]=['5 Tau',0.924334,0.2268146,4.11];
ybs[1062]=['',0.9236605,0.1090396,5.94];
ybs[1063]=['',0.937894,1.0266527,6.4];
ybs[1064]=['36 Per',0.932266,0.8048624,5.31];
ybs[1065]=['17 Eri',0.9227898,-0.0875513,4.73];
ybs[1066]=['',0.9384791,1.0110118,6.37];
ybs[1067]=['',0.9331395,0.783893,6.41];
ybs[1068]=['',0.938092,0.9604994,5.98];
ybs[1069]=['',0.9328191,0.6199389,5.9];
ybs[1070]=['',0.9186061,-0.7430727,5.78];
ybs[1071]=['',0.9200114,-0.7210107,6.12];
ybs[1072]=['',0.9444737,1.0489162,6.46];
ybs[1073]=['',0.9369986,0.6973872,5.81];
ybs[1074]=['6 Tau',0.9317817,0.1646174,5.77];
ybs[1075]=['',0.9667848,1.3228798,6.27];
ybs[1076]=['',0.9214334,-0.8258247,5.99];
ybs[1077]=['',0.9278436,-0.4460298,6.38];
ybs[1078]=['κ Ret',0.9148659,-1.0974286,4.72];
ybs[1079]=['ε Eri',0.9327694,-0.1640638,3.73];
ybs[1080]=['',0.9387096,0.3122484,6.17];
ybs[1081]=['7 Tau',0.9402238,0.4279908,5.92];
ybs[1082]=['ψ Per',0.9500745,0.8421154,4.23];
ybs[1083]=['τ5 Eri',0.9361995,-0.3765528,4.27];
ybs[1084]=['',0.9413977,0.1130152,6.49];
ybs[1085]=['',0.9298256,-0.8782542,5.68];
ybs[1086]=['',0.9401466,-0.1712346,6.25];
ybs[1087]=['',0.9208572,-1.1594351,5.83];
ybs[1088]=['',0.9366081,-0.5414435,6.2];
ybs[1089]=['',0.9586859,0.9946457,6.3];
ybs[1090]=['',0.9392493,-0.5553129,6.4];
ybs[1091]=['',0.9302335,-1.0639293,6.41];
ybs[1092]=['',0.9564084,0.7441983,6.42];
ybs[1093]=['',0.9459522,-0.1943682,5.57];
ybs[1094]=['',0.9498459,0.0112507,5.71];
ybs[1095]=['20 Eri',0.9472236,-0.3038603,5.23];
ybs[1096]=['10 Tau',0.9502123,0.0080019,4.28];
ybs[1097]=['',0.9546139,0.2703045,6.39];
ybs[1098]=['',0.9600208,0.3660285,6.5];
ybs[1099]=['',0.9364235,-1.1467977,6.75];
ybs[1100]=['',0.9760348,1.1042965,5.1];
ybs[1101]=['',0.9499995,-0.7019355,4.58];
ybs[1102]=['',1.1208047,1.5100515,5.86];
ybs[1103]=['',0.957062,-0.128027,5.85];
ybs[1104]=['',0.913604,-1.3664624,5.7];
ybs[1105]=['',0.9617644,0.2895951,6.16];
ybs[1106]=['21 Eri',0.959425,-0.0972155,5.96];
ybs[1107]=['',0.9781139,1.0476186,5.76];
ybs[1108]=['',0.969861,0.6568597,5.57];
ybs[1109]=['τ For',0.9578507,-0.4867179,6.01];
ybs[1110]=['12 Tau',0.9632725,0.0543272,5.57];
ybs[1111]=['',0.9621846,-0.0582451,6.23];
ybs[1112]=['',0.9610727,-0.1811876,6.19];
ybs[1113]=['11 Tau',0.9678743,0.4430496,6.11];
ybs[1114]=['',0.9637816,-0.0185847,6.12];
ybs[1115]=['',0.9642945,-0.2647838,6.33];
ybs[1116]=['22 Eri',0.9665031,-0.0899726,5.53];
ybs[1117]=['δ Per',0.9781608,0.8350028,3.01];
ybs[1118]=['40 Per',0.9751551,0.5937585,4.97];
ybs[1119]=['',0.9934762,1.1738236,5.8];
ybs[1120]=['',0.9689197,-0.2050369,6.49];
ybs[1121]=['13 Tau',0.9744495,0.3447928,5.69];
ybs[1122]=['',0.983363,0.8478431,6.06];
ybs[1123]=['',0.9693394,-0.3408533,6.59];
ybs[1124]=['',0.9930266,1.1065123,4.8];
ybs[1125]=['',0.9857441,0.8055348,6.11];
ybs[1126]=['ο Per',0.9835809,0.5644832,3.83];
ybs[1127]=['14 Tau',0.980872,0.3441681,6.14];
ybs[1128]=['',0.9846274,0.6372911,5.59];
ybs[1129]=['δ For',0.9727806,-0.5564696,5];
ybs[1130]=['ν Per',0.9878082,0.7440768,3.77];
ybs[1131]=['δ Eri',0.9777765,-0.1694494,3.54];
ybs[1132]=['',0.9838829,0.366218,6.1];
ybs[1133]=['',1.0082608,1.2378469,5.44];
ybs[1134]=['',0.9791324,-0.1820565,5.6];
ybs[1135]=['16 Tau',0.9854426,0.4248735,5.46];
ybs[1136]=['',0.9914296,0.7982351,5.66];
ybs[1137]=['17 Tau',0.9857502,0.4217994,3.7];
ybs[1138]=['',0.9751439,-0.6502894,4.59];
ybs[1139]=['18 Tau',0.9870232,0.4344657,5.64];
ybs[1140]=['19 Tau',0.9872157,0.4279738,4.3];
ybs[1141]=['24 Eri',0.9834879,-0.0193541,5.25];
ybs[1142]=['',0.9988411,0.9769565,6.1];
ybs[1143]=['γ Cam',1.0133077,1.2458876,4.63];
ybs[1144]=['20 Tau',0.9899119,0.4262344,3.87];
ybs[1145]=['25 Eri',0.9854,-0.0042355,5.55];
ybs[1146]=['21 Tau',0.9902668,0.4294966,5.76];
ybs[1147]=['22 Tau',0.9908843,0.4290303,6.43];
ybs[1148]=['29 Tau',0.9887502,0.1065301,5.35];
ybs[1149]=['',0.9419294,-1.3659968,6.29];
ybs[1150]=['',1.0085408,1.144558,4.47];
ybs[1151]=['23 Tau',0.992082,0.4189106,4.18];
ybs[1152]=['',0.9805376,-0.7087074,6.45];
ybs[1153]=['',1.0086408,1.1056561,5.85];
ybs[1154]=['',0.9908855,0.1196704,5.91];
ybs[1155]=['',1.001833,0.8864421,6.14];
ybs[1156]=['',1.0067606,0.9978116,6.46];
ybs[1157]=['π Eri',0.9903389,-0.2102788,4.42];
ybs[1158]=['',0.9991626,0.587354,6.57];
ybs[1159]=['',0.9988438,0.5628326,6.25];
ybs[1160]=['η Tau',0.997143,0.4216377,2.87];
ybs[1161]=['',1.0185537,1.1965789,6.32];
ybs[1162]=['',0.9834049,-0.8378855,6.49];
ybs[1163]=['',0.9817934,-0.9463121,6.3];
ybs[1164]=['',0.9852739,-0.8256417,5.73];
ybs[1165]=['',1.0050942,0.7682154,6.02];
ybs[1166]=['σ For',0.9912191,-0.5111122,5.9];
ybs[1167]=['',1.0008912,0.4096962,5.45];
ybs[1168]=['τ6 Eri',0.9931155,-0.4048533,4.23];
ybs[1169]=['30 Tau',1.0002187,0.1954093,5.07];
ybs[1170]=['β Ret',0.9791747,-1.1301456,3.85];
ybs[1171]=['',1.0092417,0.7857451,5.66];
ybs[1172]=['42 Per',1.006421,0.5784667,5.11];
ybs[1173]=['27 Tau',1.0044603,0.4207255,3.63];
ybs[1174]=['',0.9950311,-0.5209593,6.55];
ybs[1175]=['28 Tau',1.0045719,0.4221798,5.09];
ybs[1176]=['τ7 Eri',0.9966371,-0.4157666,5.24];
ybs[1177]=['',1.0015865,0.0048948,5.91];
ybs[1178]=['',1.006909,0.4147588,6.17];
ybs[1179]=['ρ For',0.9976316,-0.5256028,5.54];
ybs[1180]=['',1.0077022,0.3891497,6.07];
ybs[1181]=['',0.9969553,-0.6292405,6.21];
ybs[1182]=['',1.0008092,-0.3639071,5.81];
ybs[1183]=['',1.0095412,0.4473539,5.26];
ybs[1184]=['',1.0002246,-0.6557112,5.4];
ybs[1185]=['',1.000261,-0.6556822,4.73];
ybs[1186]=['',1.0167344,0.6005788,5.77];
ybs[1187]=['',1.0260397,1.0127403,5.8];
ybs[1188]=['',1.0150873,0.3854254,6.83];
ybs[1189]=['',1.0133248,0.2285954,6.3];
ybs[1190]=['',1.0040518,-0.630899,4.17];
ybs[1191]=['',1.0403163,1.2543911,6.34];
ybs[1192]=['',1.0173918,0.5448872,6.25];
ybs[1193]=['',1.0250323,0.8499991,5.76];
ybs[1194]=['31 Tau',1.0163851,0.1149505,5.67];
ybs[1195]=['',1.0091521,-0.6348333,6.86];
ybs[1196]=['',1.0217476,0.3033028,5.97];
ybs[1197]=['30 Eri',1.0190946,-0.0926799,5.48];
ybs[1198]=['ζ Per',1.02641,0.5573581,2.85];
ybs[1199]=['',1.0427646,1.1016796,5.03];
ybs[1200]=['',1.0412818,1.067415,5];
ybs[1201]=['',1.0210312,-0.3208509,6.22];
ybs[1202]=['',1.0351424,0.8363854,5.37];
ybs[1203]=['γ Hyi',0.9903967,-1.2947813,3.24];
ybs[1204]=['',1.0318451,0.5427282,6.1];
ybs[1205]=['43 Per',1.0381,0.8856673,5.28];
ybs[1206]=['32 Eri',1.0261222,-0.0506518,6.14];
ybs[1207]=['32 Eri',1.0261294,-0.0506858,4.79];
ybs[1208]=['τ8 Eri',1.0230022,-0.4286812,4.65];
ybs[1209]=['',1.0223926,-0.605303,5.11];
ybs[1210]=['',1.0367723,0.61315,5.49];
ybs[1211]=['',1.0214425,-0.8175583,5.93];
ybs[1212]=['',1.0301683,-0.2102925,6];
ybs[1213]=['32 Tau',1.0380526,0.3931831,5.63];
ybs[1214]=['',1.0253851,-0.7034823,5.71];
ybs[1215]=['ε Per',1.0429768,0.6991713,2.89];
ybs[1216]=['33 Tau',1.0389247,0.4053555,6.06];
ybs[1217]=['',1.0406081,0.4278048,6.16];
ybs[1218]=['',1.0436312,0.6084858,6.53];
ybs[1219]=['',1.0383002,0.1062842,6.09];
ybs[1220]=['',1.0361634,-0.169315,6.19];
ybs[1221]=['',1.0456846,0.6787469,6.3];
ybs[1222]=['',1.0254988,-0.9187403,6.46];
ybs[1223]=['ξ Per',1.0476565,0.6255259,4.04];
ybs[1224]=['',1.0508433,0.6783951,6.38];
ybs[1225]=['',1.1042121,1.4092314,5.1];
ybs[1226]=['γ Eri',1.0421636,-0.2349095,2.95];
ybs[1227]=['',1.0460477,-0.0946148,5.83];
ybs[1228]=['',1.0499767,0.1811563,6.37];
ybs[1229]=['',1.0576587,0.6464358,6.41];
ybs[1230]=['',1.048609,-0.2186098,5.6];
ybs[1231]=['',1.0310316,-1.1067744,6.14];
ybs[1232]=['',1.0542551,0.3027268,6.32];
ybs[1233]=['',1.0551462,0.318385,5.89];
ybs[1234]=['λ Tau',1.0544059,0.2188392,3.47];
ybs[1235]=['τ9 Eri',1.0501279,-0.4183168,4.66];
ybs[1236]=['',1.0813584,1.1994973,5.87];
ybs[1237]=['',1.0730423,1.0332759,5.06];
ybs[1238]=['',1.0590886,0.1753298,5.67];
ybs[1239]=['35 Eri',1.0577621,-0.0262103,5.28];
ybs[1240]=['',1.0432361,-0.9957693,6.05];
ybs[1241]=['',1.0532023,-0.5313218,5.93];
ybs[1242]=['δ Ret',1.0429445,-1.0707796,4.56];
ybs[1243]=['',1.0834244,1.144356,6.17];
ybs[1244]=['',1.0625003,-0.0038625,5.38];
ybs[1245]=['',1.0503965,-0.8991272,6.51];
ybs[1246]=['ν Tau',1.0650411,0.1053575,3.91];
ybs[1247]=['36 Tau',1.07081,0.4215448,5.47];
ybs[1248]=['40 Tau',1.0675866,0.0956913,5.33];
ybs[1249]=['',1.0685322,0.1438902,5.46];
ybs[1250]=['',1.0820181,0.943431,6.31];
ybs[1251]=['37 Tau',1.0722027,0.3862191,4.36];
ybs[1252]=['',1.0693585,0.0501599,5.36];
ybs[1253]=['',1.0654527,-0.3507563,6.46];
ybs[1254]=['',1.0663319,-0.3510049,7.01];
ybs[1255]=['',1.0882407,1.0886579,6.99];
ybs[1256]=['λ Per',1.0816463,0.8796008,4.29];
ybs[1257]=['39 Tau',1.0750012,0.3849399,5.9];
ybs[1258]=['',1.068754,-0.2887098,6.39];
ybs[1259]=['γ Ret',1.0522309,-1.0840434,4.51];
ybs[1260]=['',1.0698774,-0.2224519,5.61];
ybs[1261]=['ι Ret',1.0541256,-1.0651868,4.97];
ybs[1262]=['',1.0709926,-0.3549098,6.13];
ybs[1263]=['41 Tau',1.0807222,0.4825148,5.2];
ybs[1264]=['ψ Tau',1.0825229,0.5069709,5.23];
ybs[1265]=['',1.0949945,1.0463766,6.28];
ybs[1266]=['',0.958143,-1.4820563,6.41];
ybs[1267]=['',1.0768041,-0.1537593,6.26];
ybs[1268]=['48 Per',1.0905617,0.83353,4.04];
ybs[1269]=['',1.0757652,-0.3571955,6.34];
ybs[1270]=['',1.0748654,-0.4818011,5.59];
ybs[1271]=['',1.0941525,0.9577293,6.18];
ybs[1272]=['49 Per',1.0882951,0.6592667,6.09];
ybs[1273]=['50 Per',1.0898582,0.6647089,5.51];
ybs[1274]=['',1.0851163,0.2654374,6.01];
ybs[1275]=['',1.086444,0.3034303,5.89];
ybs[1276]=['',1.1158321,1.2595968,6.03];
ybs[1277]=['',1.1110892,1.19634,6.32];
ybs[1278]=['ω1 Tau',1.0916458,0.3430318,5.5];
ybs[1279]=['',1.0908568,0.2346334,5.95];
ybs[1280]=['',1.0820857,-0.7482418,6.59];
ybs[1281]=['',1.1000409,0.5869729,5.72];
ybs[1282]=['44 Tau',1.0991275,0.4629539,5.41];
ybs[1283]=['',1.0912328,-0.2851995,5.37];
ybs[1284]=['',1.1887618,1.4633702,5.57];
ybs[1285]=['37 Eri',1.0961976,-0.1200648,5.44];
ybs[1286]=['',1.0869337,-0.7996976,6.59];
ybs[1287]=['45 Tau',1.1007309,0.0971688,5.72];
ybs[1288]=['',1.097979,-0.153156,5.7];
ybs[1289]=['',1.0801233,-1.1200922,6.38];
ybs[1290]=['',1.1062274,0.3023145,6.09];
ybs[1291]=['',1.1190879,1.0036181,6.08];
ybs[1292]=['',1.1078218,0.3919542,6.12];
ybs[1293]=['ο1 Eri',1.1027007,-0.1185669,4.04];
ybs[1294]=['',1.0970064,-0.6148675,6.44];
ybs[1295]=['',1.1011786,-0.3545091,5.79];
ybs[1296]=['',1.1133585,0.6634031,6.45];
ybs[1297]=['δ Hor',1.0970638,-0.7321491,4.93];
ybs[1298]=['μ Per',1.1178389,0.8456527,4.14];
ybs[1299]=['',1.1955592,1.4552041,5.46];
ybs[1300]=['',1.1276911,1.0802253,5.7];
ybs[1301]=['52 Per',1.1173778,0.7073211,4.71];
ybs[1302]=['',1.1106398,0.1789952,6.23];
ybs[1303]=['',1.1103491,0.1559283,6.51];
ybs[1304]=['46 Tau',1.1104482,0.1354301,5.29];
ybs[1305]=['',1.1118084,0.2233443,6.25];
ybs[1306]=['47 Tau',1.1121848,0.1624365,4.84];
ybs[1307]=['',1.1105822,-0.0193082,6.44];
ybs[1308]=['',1.1283281,1.0105902,5.71];
ybs[1309]=['',1.1261402,0.9364409,5.19];
ybs[1310]=['',1.1151073,0.1754784,5.22];
ybs[1311]=['',1.1042853,-0.7736068,6.71];
ybs[1312]=['',1.1788459,1.4113081,5.43];
ybs[1313]=['39 Eri',1.1136459,-0.1782545,4.87];
ybs[1314]=['48 Tau',1.1203498,0.269534,6.32];
ybs[1315]=['μ Tau',1.119135,0.155944,4.29];
ybs[1316]=['',1.1186037,0.1089517,6.93];
ybs[1317]=['',1.1188506,0.1087234,6.31];
ybs[1318]=['',1.1091385,-0.7036167,6.37];
ybs[1319]=['',1.1325658,0.878549,4.61];
ybs[1320]=['ο2 Eri',1.1175417,-0.1328188,4.43];
ybs[1321]=['α Hor',1.1108205,-0.7374205,3.86];
ybs[1322]=['',1.1446476,1.1376257,5.27];
ybs[1323]=['',1.1316392,0.7362282,6.22];
ybs[1324]=['ω2 Tau',1.1269986,0.3598981,4.94];
ybs[1325]=['',1.1368221,0.8742373,5.45];
ybs[1326]=['51 Tau',1.1319468,0.3773535,5.65];
ybs[1327]=['',1.1265104,-0.112223,5.94];
ybs[1328]=['',1.1411163,0.8894491,5.55];
ybs[1329]=['',1.1316901,0.1662997,6.54];
ybs[1330]=['',1.1489691,1.0607368,5.39];
ybs[1331]=['α Ret',1.1111312,-1.0896199,3.35];
ybs[1332]=['',1.1408064,0.7304011,5.92];
ybs[1333]=['γ Dor',1.1191206,-0.8978682,4.25];
ybs[1334]=['53 Tau',1.1365088,0.3697204,5.35];
ybs[1335]=['',1.1128506,-1.0847016,5.45];
ybs[1336]=['56 Tau',1.1372995,0.380739,5.38];
ybs[1337]=['',1.1488506,0.9869286,5.88];
ybs[1338]=['54 Per',1.1412313,0.6040141,4.93];
ybs[1339]=['',1.1400681,0.5584047,6.16];
ybs[1340]=['',1.1302265,-0.360827,6];
ybs[1341]=['γ Tau',1.1379077,0.273468,3.65];
ybs[1342]=['υ4 Eri',1.1281656,-0.5891617,3.56];
ybs[1343]=['φ Tau',1.1407258,0.4780747,4.95];
ybs[1344]=['',1.137017,0.1773695,6.31];
ybs[1345]=['53 Per',1.1467875,0.8122626,4.85];
ybs[1346]=['57 Tau',1.1385965,0.2456774,5.59];
ybs[1347]=['',1.153965,1.0411958,6.19];
ybs[1348]=['',1.1317195,-0.4001816,6.07];
ybs[1349]=['',1.1407273,0.3278307,6.12];
ybs[1350]=['ε Ret',1.1204504,-1.034273,4.44];
ybs[1351]=['58 Tau',1.1414345,0.2641735,5.26];
ybs[1352]=['',1.1197069,-1.0630117,6.37];
ybs[1353]=['',1.1425992,0.2426848,6.17];
ybs[1354]=['',1.1331987,-0.5910311,6.37];
ybs[1355]=['',1.1415402,0.1077141,5.77];
ybs[1356]=['',1.1421924,0.1617262,6.53];
ybs[1357]=['',1.1410234,-0.1082941,6.27];
ybs[1358]=['',1.1412847,-0.1318031,5.85];
ybs[1359]=['',1.1337216,-0.7719017,5.34];
ybs[1360]=['',1.1305311,-0.9218549,6.09];
ybs[1361]=['',1.1447095,-0.0010056,5.86];
ybs[1362]=['',1.1406326,-0.3595194,5.38];
ybs[1363]=['60 Tau',1.1477477,0.2463953,5.72];
ybs[1364]=['χ Tau',1.1503868,0.4480109,5.37];
ybs[1365]=['',1.1493523,0.3641009,5.91];
ybs[1366]=['',1.1554996,0.7411995,6.23];
ybs[1367]=['θ Ret',1.1251695,-1.1032841,5.87];
ybs[1368]=['δ1 Tau',1.151675,0.3068698,3.76];
ybs[1369]=['',1.1442722,-0.4483382,6.01];
ybs[1370]=['',1.1544201,0.3669001,5.99];
ybs[1371]=['63 Tau',1.1537544,0.29351,5.64];
ybs[1372]=['55 Per',1.1590102,0.5963751,5.73];
ybs[1373]=['62 Tau',1.1565191,0.4248225,6.36];
ybs[1374]=['56 Per',1.1596,0.5933925,5.76];
ybs[1375]=['δ2 Tau',1.1567419,0.3051409,4.8];
ybs[1376]=['66 Tau',1.1555014,0.165812,5.12];
ybs[1377]=['',1.1714519,1.0057189,6.32];
ybs[1378]=['ξ Eri',1.1543333,-0.0646815,5.17];
ybs[1379]=['',1.1511626,-0.4337559,5.83];
ybs[1380]=['',1.1605215,0.3330215,5.98];
ybs[1381]=['',1.1509195,-0.6196819,6.39];
ybs[1382]=['κ1 Tau',1.1624412,0.3897806,4.22];
ybs[1383]=['κ2 Tau',1.1626491,0.3881367,5.28];
ybs[1384]=['δ3 Tau',1.1628373,0.3135817,4.29];
ybs[1385]=['',1.1659708,0.5493856,5.28];
ybs[1386]=['70 Tau',1.1633528,0.2788972,6.46];
ybs[1387]=['υ Tau',1.1665604,0.398845,4.28];
ybs[1388]=['43 Eri',1.1549424,-0.5930185,3.96];
ybs[1389]=['71 Tau',1.1665073,0.2732636,4.49];
ybs[1390]=['η Ret',1.1435316,-1.1055959,5.24];
ybs[1391]=['π Tau',1.1676154,0.2574715,4.69];
ybs[1392]=['',1.1663284,0.150601,6.06];
ybs[1393]=['',1.1588525,-0.6059547,6.55];
ybs[1394]=['72 Tau',1.1708581,0.4020284,5.53];
ybs[1395]=['',1.169029,0.0369613,6.23];
ybs[1396]=['',1.20237,1.2664826,5.94];
ybs[1397]=['',1.171326,0.1963597,5.88];
ybs[1398]=['',1.1739639,0.378001,5.72];
ybs[1399]=['',1.1600709,-0.7700711,6.39];
ybs[1400]=['',1.1543195,-0.9953945,6.29];
ybs[1401]=['',1.1779784,0.5305609,6.4];
ybs[1402]=['75 Tau',1.1756649,0.286189,4.97];
ybs[1403]=['76 Tau',1.1753989,0.2579345,5.9];
ybs[1404]=['ε Tau',1.176521,0.3354157,3.53];
ybs[1405]=['',1.1680046,-0.4196304,6.11];
ybs[1406]=['θ1 Tau',1.1762422,0.2792504,3.84];
ybs[1407]=['θ2 Tau',1.1766177,0.2776548,3.4];
ybs[1408]=['',1.1735973,0.0330998,6.15];
ybs[1409]=['79 Tau',1.1772989,0.2283771,5.03];
ybs[1410]=['',1.1756568,0.0247578,5.55];
ybs[1411]=['',1.1577451,-1.0681271,5.94];
ybs[1412]=['1 Cam',1.1930355,0.9415535,5.77];
ybs[1413]=['',1.1676848,-0.8187197,6.1];
ybs[1414]=['',1.1857981,0.5671424,6.21];
ybs[1415]=['',1.1810596,0.1842871,6.79];
ybs[1416]=['',1.1755668,-0.3389547,5.96];
ybs[1417]=['80 Tau',1.183077,0.2735819,5.58];
ybs[1418]=['',1.177786,-0.2270824,5.6];
ybs[1419]=['',1.1894437,0.6989481,6.26];
ybs[1420]=['',1.1824705,0.1797616,6.48];
ybs[1421]=['δ Men',1.1207173,-1.3992605,5.69];
ybs[1422]=['',1.1849192,0.2832801,4.78];
ybs[1423]=['81 Tau',1.1852824,0.274519,5.48];
ybs[1424]=['',1.1725601,-0.7316788,6.44];
ybs[1425]=['83 Tau',1.1851083,0.2401798,5.4];
ybs[1426]=['',1.1823364,-0.2365822,6.24];
ybs[1427]=['85 Tau',1.1905896,0.2772983,6.02];
ybs[1428]=['',1.1774478,-0.8111912,6.16];
ybs[1429]=['57 Per',1.1983912,0.7522296,6.09];
ybs[1430]=['',1.1692524,-1.0905291,5.75];
ybs[1431]=['',1.1912365,0.0950555,6.39];
ybs[1432]=['45 Eri',1.1902058,-0.0001313,4.91];
ybs[1433]=['',1.1878755,-0.2375121,6.21];
ybs[1434]=['',1.1838132,-0.621629,5.96];
ybs[1435]=['',1.2132528,1.1221786,5.94];
ybs[1436]=['',1.1933691,-0.0553808,5.81];
ybs[1437]=['',1.1980085,0.315073,6.25];
ybs[1438]=['δ Cae',1.1840835,-0.7839502,5.07];
ybs[1439]=['ρ Tau',1.1992224,0.2597052,4.65];
ybs[1440]=['',1.2031055,0.5060817,5.88];
ybs[1441]=['',1.1988615,0.1649101,6.01];
ybs[1442]=['',1.1963998,-0.1876236,6.06];
ybs[1443]=['',1.2002225,0.0978096,5.68];
ybs[1444]=['46 Eri',1.1988994,-0.1169949,5.72];
ybs[1445]=['',1.2003073,-0.1187231,6.09];
ybs[1446]=['47 Eri',1.200079,-0.1430458,5.11];
ybs[1447]=['',1.2000654,-0.1559419,5.26];
ybs[1448]=['υ1 Eri',1.196422,-0.5189021,4.51];
ybs[1449]=['58 Per',1.2125982,0.7208057,4.25];
ybs[1450]=['',1.2075075,0.3476083,6.36];
ybs[1451]=['ν Men',1.1321978,-1.4231225,5.79];
ybs[1452]=['α Tau',1.208313,0.2887458,0.85];
ybs[1453]=['88 Tau',1.2069704,0.1779484,4.25];
ybs[1454]=['',1.211001,0.4079768,6.02];
ybs[1455]=['',1.2045738,-0.1693251,6.37];
ybs[1456]=['',1.2032895,-0.3470657,6.13];
ybs[1457]=['',1.2082072,-0.0624341,6.33];
ybs[1458]=['ν Eri',1.2094871,-0.0579081,3.93];
ybs[1459]=['υ2 Eri',1.2052906,-0.5328013,3.82];
ybs[1460]=['α Dor',1.197187,-0.9600946,3.27];
ybs[1461]=['2 Cam',1.2276667,0.9338579,5.35];
ybs[1462]=['3 Cam',1.2273884,0.9269933,5.05];
ybs[1463]=['',1.2586405,1.337647,6.49];
ybs[1464]=['',1.2135785,0.0180219,5.31];
ybs[1465]=['',1.2198813,0.47078,6.47];
ybs[1466]=['',1.218669,0.3616066,5.92];
ybs[1467]=['89 Tau',1.2180596,0.2804254,5.79];
ybs[1468]=['90 Tau',1.2179625,0.2189462,4.27];
ybs[1469]=['51 Eri',1.215111,-0.0425726,5.23];
ybs[1470]=['',1.1944959,-1.0958534,5.79];
ybs[1471]=['',1.2109478,-0.5355059,6.3];
ybs[1472]=['',1.2237135,0.440725,6.22];
ybs[1473]=['σ1 Tau',1.2224021,0.2763412,5.07];
ybs[1474]=['σ2 Tau',1.2229367,0.2784056,4.69];
ybs[1475]=['',1.2219527,0.1379565,5.39];
ybs[1476]=['53 Eri',1.2172967,-0.2490585,3.87];
ybs[1477]=['',1.2335602,0.8435758,5.67];
ybs[1478]=['',1.2204729,-0.2110008,5.01];
ybs[1479]=['93 Tau',1.2262374,0.2134692,5.46];
ybs[1480]=['',1.138423,-1.4461574,6.76];
ybs[1481]=['',1.242784,1.0393877,6.5];
ybs[1482]=['',1.2223049,-0.2500312,5.45];
ybs[1483]=['',1.2246842,-0.0177944,6.1];
ybs[1484]=['',1.2349318,0.6686815,5.99];
ybs[1485]=['',1.2323221,0.4999945,5.78];
ybs[1486]=['',1.2706919,1.3259335,6.06];
ybs[1487]=['',1.2085093,-1.0828534,5.4];
ybs[1488]=['',1.2422202,0.8727623,5.87];
ybs[1489]=['59 Per',1.2398377,0.7574189,5.29];
ybs[1490]=['',1.2254109,-0.426722,5.58];
ybs[1491]=['54 Eri',1.2269996,-0.3427596,4.32];
ybs[1492]=['τ Tau',1.2361199,0.4012365,4.28];
ybs[1493]=['',1.2195756,-0.9012736,6.44];
ybs[1494]=['95 Tau',1.2404566,0.4209857,6.13];
ybs[1495]=['',1.2454251,0.7124143,6.08];
ybs[1496]=['',1.2432864,0.5741583,6.45];
ybs[1497]=['α Cae',1.2266617,-0.730087,4.45];
ybs[1498]=['β Cae',1.2334109,-0.6477279,5.05];
ybs[1499]=['',1.2243212,-1.0281816,6.53];
ybs[1500]=['55 Eri',1.2410018,-0.1529239,6.82];
ybs[1501]=['55 Eri',1.2410382,-0.1529675,6.7];
ybs[1502]=['',1.2452912,0.1950835,5.4];
ybs[1503]=['56 Eri',1.2432426,-0.147866,5.9];
ybs[1504]=['',1.2384545,-0.536403,5.68];
ybs[1505]=['',1.2766327,1.2386657,6.37];
ybs[1506]=['4 Cam',1.2630437,0.9911207,5.3];
ybs[1507]=['',1.2512575,0.4129253,6.35];
ybs[1508]=['',1.2431359,-0.3252448,5.53];
ybs[1509]=['',1.2564244,0.7041208,5.97];
ybs[1510]=['',1.2634258,0.9709663,6.26];
ybs[1511]=['λ Pic',1.2357998,-0.8805058,5.31];
ybs[1512]=['',1.2535936,0.3275214,6.01];
ybs[1513]=['',1.2405482,-0.7161612,6.25];
ybs[1514]=['',1.2522823,0.2048364,5.37];
ybs[1515]=['μ Eri',1.2495583,-0.0562656,4.02];
ybs[1516]=['',1.2471388,-0.3709256,5.72];
ybs[1517]=['',1.2534936,-0.0510312,6.33];
ybs[1518]=['',1.3250578,1.4175264,5.07];
ybs[1519]=['',1.2499801,-0.5929607,6.86];
ybs[1520]=['',1.2528391,-0.4896853,6.19];
ybs[1521]=['',1.2501672,-0.6863654,6.05];
ybs[1522]=['',1.2816502,1.1088669,5.44];
ybs[1523]=['',1.2673327,0.5692859,5.86];
ybs[1524]=['',1.2668363,0.5491959,5.58];
ybs[1525]=['κ Dor',1.2418764,-1.0419834,5.27];
ybs[1526]=['',1.2112278,-1.3547591,6.05];
ybs[1527]=['58 Eri',1.2583351,-0.2950363,5.51];
ybs[1528]=['',1.2701299,0.6548026,4.88];
ybs[1529]=['',1.263896,0.0631452,6.03];
ybs[1530]=['',1.2761709,0.8511865,5.66];
ybs[1531]=['',1.2630299,-0.0985098,5.78];
ybs[1532]=['96 Tau',1.2685836,0.2780898,6.08];
ybs[1533]=['59 Eri',1.26244,-0.2844837,5.77];
ybs[1534]=['ζ Cae',1.2588602,-0.5234284,6.37];
ybs[1535]=['',1.2441181,-1.1030205,6.46];
ybs[1536]=['μ Men',1.2343778,-1.2374207,5.54];
ybs[1537]=['α Cam',1.2906773,1.1583763,4.29];
ybs[1538]=['π3 Ori',1.2687776,0.1220083,3.19];
ybs[1539]=['π2 Ori',1.2722009,0.1558428,4.36];
ybs[1540]=['',1.2675858,-0.2398164,6.26];
ybs[1541]=['',1.2852318,0.9227305,6.41];
ybs[1542]=['97 Tau',1.275832,0.3293132,5.1];
ybs[1543]=['',1.261447,-0.7670762,6.72];
ybs[1544]=['60 Eri',1.2696492,-0.2825367,5.03];
ybs[1545]=['',1.2829782,0.7437647,5.71];
ybs[1546]=['2 Aur',1.2819821,0.6410777,4.78];
ybs[1547]=['π4 Ori',1.2747022,0.0983252,3.69];
ybs[1548]=['',1.2770832,0.1745923,6.11];
ybs[1549]=['',1.282292,0.4873909,5.97];
ybs[1550]=['5 Cam',1.2936768,0.9649247,5.52];
ybs[1551]=['ο1 Ori',1.280745,0.249209,4.74];
ybs[1552]=['',1.2690292,-0.7206772,6.07];
ybs[1553]=['',1.2920427,0.769479,6.08];
ybs[1554]=['',1.2745576,-0.6087327,5.86];
ybs[1555]=['ω Eri',1.2817539,-0.0946812,4.39];
ybs[1556]=['',1.2981356,0.9232088,5.75];
ybs[1557]=['5 Ori',1.284097,0.0442579,5.33];
ybs[1558]=['ι Pic',1.2711008,-0.9325739,5.61];
ybs[1559]=['ι Pic',1.2711808,-0.9325449,6.42];
ybs[1560]=['',1.2864702,0.0278722,6.61];
ybs[1561]=['',1.2915486,0.3405547,6.37];
ybs[1562]=['π5 Ori',1.2878986,0.0430736,3.72];
ybs[1563]=['7 Cam',1.3033073,0.9386078,4.47];
ybs[1564]=['6 Ori',1.2904717,0.1998972,5.19];
ybs[1565]=['π1 Ori',1.2909357,0.1776386,4.65];
ybs[1566]=['',1.2904301,0.136246,5.33];
ybs[1567]=['',1.3289577,1.2966568,6.06];
ybs[1568]=['',1.2981012,0.6317284,6.07];
ybs[1569]=['',1.2904312,0.0086331,5.99];
ybs[1570]=['',1.2973409,0.4296782,6.37];
ybs[1571]=['',1.2951915,0.2629638,5.81];
ybs[1572]=['ι Aur',1.3008623,0.5793153,2.69];
ybs[1573]=['',1.2954958,0.0946988,6.5];
ybs[1574]=['',1.2910974,-0.2917056,5.7];
ybs[1575]=['ο2 Ori',1.2974725,0.2363341,4.07];
ybs[1576]=['',1.2919655,-0.2860735,5.72];
ybs[1577]=['62 Eri',1.2970693,-0.0897951,5.51];
ybs[1578]=['',1.292503,-0.4485647,6.72];
ybs[1579]=['',1.2893641,-0.6911753,6.1];
ybs[1580]=['',1.3019483,0.2998423,5.48];
ybs[1581]=['99 Tau',1.3040905,0.4184341,5.79];
ybs[1582]=['',1.3369872,1.2878246,6.66];
ybs[1583]=['8 Cam',1.3141083,0.9281762,6.08];
ybs[1584]=['',1.339053,1.2931106,5.96];
ybs[1585]=['98 Tau',1.3056343,0.4376593,5.81];
ybs[1586]=['',1.3010412,-0.0181702,6.23];
ybs[1587]=['ω Aur',1.3109456,0.6617513,4.94];
ybs[1588]=['',1.3228668,1.0664362,6.03];
ybs[1589]=['',1.3291572,1.1666908,6.19];
ybs[1590]=['',1.3026631,-0.2479312,6.15];
ybs[1591]=['',1.3049062,-0.0381654,6.35];
ybs[1592]=['',1.2879172,-1.0213662,6.12];
ybs[1593]=['',1.2808087,-1.1632217,6.41];
ybs[1594]=['5 Aur',1.3155922,0.6880012,5.95];
ybs[1595]=['',1.3089297,0.2542632,6.09];
ybs[1596]=['π6 Ori',1.3066258,0.0303652,4.47];
ybs[1597]=['6 Aur',1.3159612,0.6925433,6.58];
ybs[1598]=['β Cam',1.3307371,1.0553257,4.03];
ybs[1599]=['',1.3081598,-0.2853679,5.66];
ybs[1600]=['ε Aur',1.3230806,0.765283,2.99];
ybs[1601]=['',1.2776364,-1.2632589,6.28];
ybs[1602]=['',1.3107677,-0.2579704,7.71];
ybs[1603]=['63 Eri',1.3119152,-0.1786907,5.38];
ybs[1604]=['',1.3153855,0.0635314,7.03];
ybs[1605]=['',1.31548,0.0635458,6.66];
ybs[1606]=['64 Eri',1.31224,-0.218383,4.79];
ybs[1607]=['ζ Aur',1.3251615,0.7173264,3.75];
ybs[1608]=['',1.3157459,-0.0356185,6.32];
ybs[1609]=['',1.3163084,-0.0999834,6.22];
ybs[1610]=['',1.3288102,0.7237054,6.14];
ybs[1611]=['',1.475502,1.4963298,6.51];
ybs[1612]=['ψ Eri',1.3189796,-0.1247813,4.81];
ybs[1613]=['',1.3209528,0.0130289,5.92];
ybs[1614]=['',1.3216837,0.0285029,6.24];
ybs[1615]=['ι Tau',1.327068,0.3772307,4.64];
ybs[1616]=['',1.3185351,-0.3495454,4.91];
ybs[1617]=['11 Cam',1.3424611,1.0296547,5.08];
ybs[1618]=['12 Cam',1.3427357,1.0305027,6.08];
ybs[1619]=['',1.3442595,1.0680056,6.04];
ybs[1620]=['',1.3248249,-0.0730563,5.85];
ybs[1621]=['',1.3323912,0.5326389,6.14];
ybs[1622]=['',1.3340929,0.5644982,6.62];
ybs[1623]=['',1.321537,-0.4581633,5.02];
ybs[1624]=['η Men',1.2858486,-1.3074201,5.47];
ybs[1625]=['',1.3430039,0.9499501,7.24];
ybs[1626]=['',1.3184159,-0.6927841,6.03];
ybs[1627]=['',1.3339875,0.4837912,6.6];
ybs[1628]=['',1.3325781,0.3717771,6.19];
ybs[1629]=['1 Lep',1.32421,-0.3974301,5.75];
ybs[1630]=['',1.3222752,-0.5540947,5.94];
ybs[1631]=['',1.3593276,1.2158015,6.41];
ybs[1632]=['9 Aur',1.3441327,0.9009382,5];
ybs[1633]=['11 Ori',1.3332991,0.2692571,4.68];
ybs[1634]=['',1.3403556,0.6276009,6.52];
ybs[1635]=['',1.3293596,-0.2503846,6.41];
ybs[1636]=['η Aur',1.342789,0.7200656,3.17];
ybs[1637]=['',1.3376573,0.3460833,6.44];
ybs[1638]=['',1.3726911,1.2909556,5.43];
ybs[1639]=['',1.3442517,0.7539276,6.2];
ybs[1640]=['',1.3291287,-0.4252374,5.61];
ybs[1641]=['',1.3342397,-0.0526515,6.05];
ybs[1642]=['',1.3589674,1.1334218,6.41];
ybs[1643]=['',1.3364849,0.0209494,6.17];
ybs[1644]=['η1 Pic',1.323266,-0.8574353,5.38];
ybs[1645]=['',1.3835281,1.3350274,6.37];
ybs[1646]=['',1.3284554,-0.7281776,6.31];
ybs[1647]=['γ1 Cae',1.3309615,-0.6188948,4.55];
ybs[1648]=['γ2 Cae',1.3310755,-0.6227687,6.34];
ybs[1649]=['ε Lep',1.3360615,-0.3900513,3.19];
ybs[1650]=['',1.335094,-0.4560475,5.73];
ybs[1651]=['104 Tau',1.3459766,0.3257994,5];
ybs[1652]=['66 Eri',1.3422793,-0.0808568,5.12];
ybs[1653]=['106 Tau',1.3475908,0.3567473,5.3];
ybs[1654]=['103 Tau',1.3490448,0.4238869,5.5];
ybs[1655]=['105 Tau',1.34815,0.3791982,5.89];
ybs[1656]=['',1.3413738,-0.2286315,6.05];
ybs[1657]=['13 Ori',1.3465158,0.1656983,6.17];
ybs[1658]=['η2 Pic',1.3326428,-0.8648927,5.03];
ybs[1659]=['14 Ori',1.3475487,0.1487038,5.34];
ybs[1660]=['',1.3449052,-0.2176175,5.97];
ybs[1661]=['β Eri',1.3470153,-0.0883937,2.79];
ybs[1662]=['',1.3324635,-0.9491872,6.27];
ybs[1663]=['',1.3614443,0.8200035,5.68];
ybs[1664]=['',1.3591845,0.6514032,6.02];
ybs[1665]=['',1.3563446,0.4895914,6.01];
ybs[1666]=['',1.3490336,-0.1508558,5.78];
ybs[1667]=['16 Ori',1.3538936,0.1719256,5.43];
ybs[1668]=['68 Eri',1.3508588,-0.0773997,5.12];
ybs[1669]=['ζ Dor',1.3343617,-1.0026898,4.72];
ybs[1670]=['',1.3729815,1.0798261,6.17];
ybs[1671]=['15 Ori',1.3556986,0.2725895,4.82];
ybs[1672]=['β Men',1.3198964,-1.2442505,5.31];
ybs[1673]=['14 Cam',1.3751393,1.0945076,6.5];
ybs[1674]=['λ Eri',1.3525651,-0.1524178,4.27];
ybs[1675]=['',1.3477017,-0.623024,6.52];
ybs[1676]=['',1.3567603,-0.0095016,6.1];
ybs[1677]=['',1.3059499,-1.3661559,6.29];
ybs[1678]=['',1.397988,1.2790691,5.74];
ybs[1679]=['',1.364411,0.2803999,5.18];
ybs[1680]=['',1.3606964,-0.03898,6.25];
ybs[1681]=['',1.4201766,1.3831097,5.05];
ybs[1682]=['',1.3622313,-0.0431179,5.9];
ybs[1683]=['',1.382006,1.0371472,6.15];
ybs[1684]=['μ Aur',1.3728278,0.6720191,4.86];
ybs[1685]=['',1.3639329,0.0093361,6.67];
ybs[1686]=['',1.3642317,0.0184501,5.89];
ybs[1687]=['',1.3794506,0.9290861,6.2];
ybs[1688]=['',1.3622165,-0.2064519,5.68];
ybs[1689]=['',1.3589729,-0.451845,6.41];
ybs[1690]=['',1.3426205,-1.1061488,5.2];
ybs[1691]=['ι Lep',1.3662225,-0.2068076,4.45];
ybs[1692]=['',1.3685991,-0.1053739,5.91];
ybs[1693]=['ρ Ori',1.3709899,0.0502766,4.46];
ybs[1694]=['',1.3622468,-0.6523164,6.57];
ybs[1695]=['',1.3343497,-1.2743526,6.27];
ybs[1696]=['',1.3719892,0.0346881,6.09];
ybs[1697]=['μ Lep',1.3688518,-0.2824964,3.31];
ybs[1698]=['',1.3730827,0.010116,6.32];
ybs[1699]=['',1.371817,-0.1418663,6.37];
ybs[1700]=['κ Lep',1.370262,-0.2255281,4.36];
ybs[1701]=['14 Aur',1.3812062,0.5708339,5.02];
ybs[1702]=['',1.3906488,0.9355639,6.5];
ybs[1703]=['α Aur',1.3874652,0.8031322,0.08];
ybs[1704]=['',1.3773548,0.0903214,5.5];
ybs[1705]=['',1.3735556,-0.2545983,6.21];
ybs[1706]=['108 Tau',1.3810591,0.3892663,6.27];
ybs[1707]=['',1.3851876,0.5991744,5.96];
ybs[1708]=['β Ori',1.3761058,-0.1428141,0.12];
ybs[1709]=['',1.5265006,1.4910155,6.6];
ybs[1710]=['',1.3718126,-0.6249355,6.98];
ybs[1711]=['ξ Men',1.2955116,-1.4389274,5.85];
ybs[1712]=['',1.3796576,-0.0242683,6.15];
ybs[1713]=['18 Ori',1.383357,0.198265,5.56];
ybs[1714]=['15 Cam',1.4005288,1.0146306,6.13];
ybs[1715]=['',1.4050625,1.0937983,5.61];
ybs[1716]=['',1.3748814,-0.6275877,5.76];
ybs[1717]=['',1.3941559,0.7471686,5.48];
ybs[1718]=['',1.3792795,-0.4699237,5.07];
ybs[1719]=['',1.3857696,0.0343016,6.42];
ybs[1720]=['',1.3958163,0.706548,6.18];
ybs[1721]=['16 Aur',1.3933307,0.5827499,4.54];
ybs[1722]=['',1.3713751,-0.9077759,6.05];
ybs[1723]=['17 Aur',1.3939438,0.5896526,6.14];
ybs[1724]=['λ Aur',1.397842,0.7001596,4.71];
ybs[1725]=['',1.3806129,-0.6092616,6.66];
ybs[1726]=['',1.3857064,-0.2988627,6.56];
ybs[1727]=['',1.3969325,0.5893179,5.41];
ybs[1728]=['',1.4019869,0.7756625,6.62];
ybs[1729]=['18 Aur',1.3986591,0.5934553,6.49];
ybs[1730]=['τ Ori',1.3895336,-0.1191485,3.6];
ybs[1731]=['',1.4048307,0.8199602,6.54];
ybs[1732]=['',1.3896086,-0.2356543,5.5];
ybs[1733]=['',1.4027147,0.7173769,5.52];
ybs[1734]=['109 Tau',1.3977009,0.3859514,4.94];
ybs[1735]=['19 Aur',1.401371,0.5929708,5.03];
ybs[1736]=['',1.3974968,0.3517142,6.08];
ybs[1737]=['',1.3790397,-0.9104212,6.49];
ybs[1738]=['ο Col',1.3880375,-0.6087262,4.83];
ybs[1739]=['θ Dor',1.3689768,-1.1722635,4.83];
ybs[1740]=['',1.4492398,1.3611803,6.56];
ybs[1741]=['21 Ori',1.3967043,0.0456037,5.34];
ybs[1742]=['',1.3945677,-0.3161273,5.96];
ybs[1743]=['',1.3983324,-0.0243528,6.34];
ybs[1744]=['ρ Aur',1.4095606,0.7299027,5.23];
ybs[1745]=['',1.4053794,0.4882296,6.33];
ybs[1746]=['16 Cam',1.4179305,1.0046044,5.28];
ybs[1747]=['',1.4064155,0.5163762,5.76];
ybs[1748]=['',1.3965257,-0.3229374,6.36];
ybs[1749]=['',1.3965842,-0.3227581,6.54];
ybs[1750]=['',1.4048973,0.346107,6.18];
ybs[1751]=['λ Lep',1.3979307,-0.2296808,4.29];
ybs[1752]=['ν Lep',1.3997461,-0.2146546,5.3];
ybs[1753]=['',1.3966707,-0.4773801,5.99];
ybs[1754]=['',1.4019396,-0.0933771,6.39];
ybs[1755]=['',1.4140228,0.7163689,5.54];
ybs[1756]=['',1.4060617,0.0703038,6.57];
ybs[1757]=['',1.4014821,-0.3704089,4.71];
ybs[1758]=['',1.4079597,0.1473862,5.8];
ybs[1759]=['',1.4068409,-0.0069866,5.68];
ybs[1760]=['22 Ori',1.4078526,-0.0063968,4.73];
ybs[1761]=['',1.4005153,-0.6053194,6.34];
ybs[1762]=['ζ Pic',1.3953828,-0.882945,5.45];
ybs[1763]=['22 Aur',1.415861,0.5053066,6.46];
ybs[1764]=['',1.4078386,-0.2398106,6.56];
ybs[1765]=['23 Ori',1.4126439,0.0621334,5];
ybs[1766]=['',1.4071262,-0.4320917,5.06];
ybs[1767]=['',1.4046239,-0.5991546,6.09];
ybs[1768]=['σ Aur',1.4217657,0.6527579,4.99];
ybs[1769]=['110 Tau',1.4165127,0.2917253,6.08];
ybs[1770]=['',1.4214346,0.5452169,6.28];
ybs[1771]=['',1.421446,0.5438061,5.94];
ybs[1772]=['',1.4156852,0.0931612,6.35];
ybs[1773]=['',1.414364,-0.146616,5.9];
ybs[1774]=['',1.4241103,0.6085922,6.55];
ybs[1775]=['111 Tau',1.4200037,0.3036556,4.99];
ybs[1776]=['',1.4163385,-0.0025228,5.7];
ybs[1777]=['',1.4169796,-0.0148721,6.11];
ybs[1778]=['8 Lep',1.4150469,-0.2428092,5.25];
ybs[1779]=['29 Ori',1.4171674,-0.136013,4.14];
ybs[1780]=['',1.4132951,-0.4658354,6.49];
ybs[1781]=['',1.420332,0.0413219,6.32];
ybs[1782]=['27 Ori',1.419706,-0.0152985,5.08];
ybs[1783]=['η Ori',1.4196399,-0.0415753,3.36];
ybs[1784]=['ψ1 Ori',1.4209425,0.0324828,4.95];
ybs[1785]=['γ Ori',1.4227553,0.1110779,1.64];
ybs[1786]=['β Tau',1.4285475,0.4995402,1.65];
ybs[1787]=['',1.4191827,-0.2960293,5.65];
ybs[1788]=['',1.4136196,-0.6922541,5.71];
ybs[1789]=['',1.4315048,0.6190857,6.15];
ybs[1790]=['',1.4310654,0.6004889,5.94];
ybs[1791]=['',1.431192,0.5807858,6.15];
ybs[1792]=['',1.4148223,-0.6513814,6.82];
ybs[1793]=['113 Tau',1.4272764,0.2917217,6.25];
ybs[1794]=['',1.4218037,-0.1800226,5.61];
ybs[1795]=['',1.4242466,-0.0092459,6.57];
ybs[1796]=['κ Pic',1.4080114,-0.9794541,6.11];
ybs[1797]=['17 Cam',1.447873,1.1009453,5.42];
ybs[1798]=['',1.425427,0.0093399,6.16];
ybs[1799]=['',1.4323043,0.5274786,5.74];
ybs[1800]=['φ Aur',1.4347058,0.6019518,5.07];
ybs[1801]=['',1.4263691,-0.0960652,6.23];
ybs[1802]=['',1.4293815,0.1201327,6.42];
ybs[1803]=['115 Tau',1.4320007,0.3137389,5.42];
ybs[1804]=['',1.4321832,0.2665371,6.16];
ybs[1805]=['114 Tau',1.4341694,0.3831074,4.88];
ybs[1806]=['ψ2 Ori',1.4300987,0.0542694,4.59];
ybs[1807]=['',1.4257404,-0.3435037,5.89];
ybs[1808]=['',1.4200379,-0.7716291,6.08];
ybs[1809]=['116 Tau',1.4345154,0.2772912,5.5];
ybs[1810]=['',1.3559722,-1.4228141,6.51];
ybs[1811]=['117 Tau',1.4357227,0.3011081,5.77];
ybs[1812]=['',1.4307147,-0.2074683,6.35];
ybs[1813]=['θ Pic',1.4188148,-0.9128343,6.27];
ybs[1814]=['',1.4380245,0.2389704,6.35];
ybs[1815]=['',1.4352382,0.0228933,6.41];
ybs[1816]=['118 Tau',1.4414384,0.4391833,5.47];
ybs[1817]=['',1.4433414,0.5096188,6.24];
ybs[1818]=['',1.4327157,-0.372837,6.07];
ybs[1819]=['',1.4488355,0.7238588,6];
ybs[1820]=['',1.4484972,0.6953038,6.37];
ybs[1821]=['',1.4391096,-0.0575003,6.39];
ybs[1822]=['',1.4296457,-0.7143595,5.87];
ybs[1823]=['18 Cam',1.4576182,0.9988944,6.48];
ybs[1824]=['β Lep',1.4355142,-0.3620886,2.84];
ybs[1825]=['',1.4410618,-0.0599276,5.79];
ybs[1826]=['',1.447664,0.3922573,6.29];
ybs[1827]=['',1.4461718,0.2683025,5.94];
ybs[1828]=['',1.4434849,0.031446,5.78];
ybs[1829]=['31 Ori',1.4426147,-0.0188423,4.71];
ybs[1830]=['',1.4349164,-0.6495632,5.57];
ybs[1831]=['λ Dor',1.4249711,-1.0279691,5.14];
ybs[1832]=['',1.4453742,0.0735927,6.21];
ybs[1833]=['',1.4389566,-0.525409,6.75];
ybs[1834]=['32 Ori',1.4474113,0.1040259,4.2];
ybs[1835]=['',1.4450956,-0.1295441,6.33];
ybs[1836]=['33 Ori',1.4493252,0.0576695,5.46];
ybs[1837]=['χ Aur',1.4567782,0.5620531,4.76];
ybs[1838]=['',1.4926034,1.3099044,6.17];
ybs[1839]=['119 Tau',1.4540358,0.3247362,4.38];
ybs[1840]=['',1.4605146,0.7351305,6.55];
ybs[1841]=['',1.4540877,0.2979209,5.46];
ybs[1842]=['',1.449495,-0.1168736,6.22];
ybs[1843]=['10 Lep',1.4480832,-0.3639278,5.55];
ybs[1844]=['',1.4599891,0.5726797,6.48];
ybs[1845]=['δ Ori',1.4525649,-0.0047606,6.85];
ybs[1846]=['δ Ori',1.4525573,-0.0050175,2.23];
ybs[1847]=['',1.4794621,1.1642284,6.26];
ybs[1848]=['',1.4608349,0.606271,6.27];
ybs[1849]=['υ Ori',1.4520152,-0.1272287,4.62];
ybs[1850]=['',1.4426882,-0.8214432,5.46];
ybs[1851]=['19 Cam',1.4789467,1.1198733,6.15];
ybs[1852]=['120 Tau',1.4597649,0.3237811,5.69];
ybs[1853]=['',1.4264227,-1.1974491,6.03];
ybs[1854]=['',1.4603528,0.3575329,6.18];
ybs[1855]=['',1.4554934,-0.0275858,5.35];
ybs[1856]=['ε Col',1.447895,-0.6188673,3.87];
ybs[1857]=['',1.4573732,-0.0297949,6.46];
ybs[1858]=['35 Ori',1.46128,0.2498685,5.64];
ybs[1859]=['α Lep',1.4551793,-0.3108575,2.58];
ybs[1860]=['',1.4749164,0.9501314,5.73];
ybs[1861]=['',1.4375011,-1.0873656,6.59];
ybs[1862]=['',1.4591497,-0.0199852,5.34];
ybs[1863]=['',1.473007,0.8329596,6.11];
ybs[1864]=['',1.449056,-0.801339,5.86];
ybs[1865]=['',1.4611306,0.0247598,6.59];
ybs[1866]=['38 Ori',1.4625893,0.0659326,5.36];
ybs[1867]=['',1.4615166,-0.0178852,6.22];
ybs[1868]=['',1.4615111,-0.0254773,5.93];
ybs[1869]=['121 Tau',1.468352,0.4197453,5.38];
ybs[1870]=['φ1 Ori',1.4651242,0.1658049,4.41];
ybs[1871]=['',1.4549237,-0.6719808,5.48];
ybs[1872]=['',1.4705496,0.4829711,6.27];
ybs[1873]=['λ Ori',1.4665266,0.1735644,3.54];
ybs[1874]=['λ Ori',1.4665412,0.173579,5.61];
ybs[1875]=['',1.4562496,-0.6131072,5.78];
ybs[1876]=['',1.4415212,-1.1155306,6.19];
ybs[1877]=['',1.4669068,0.1789016,5.6];
ybs[1878]=['',1.4752301,0.7014789,6.09];
ybs[1879]=['',1.6007105,1.4814271,6.11];
ybs[1880]=['',1.4655144,-0.1046978,5.67];
ybs[1881]=['',1.4656455,-0.1045721,4.78];
ybs[1882]=['',1.4597195,-0.5207704,6.53];
ybs[1883]=['',1.4730095,0.4528987,6.49];
ybs[1884]=['',1.4670794,-0.0782443,6.56];
ybs[1885]=['',1.4671323,-0.0770565,6.24];
ybs[1886]=['42 Ori',1.4671709,-0.0842657,4.59];
ybs[1887]=['θ1 Ori',1.4666237,-0.0938448,6.73];
ybs[1888]=['θ1 Ori',1.4666384,-0.0938108,7.96];
ybs[1889]=['θ1 Ori',1.4666673,-0.0938885,5.13];
ybs[1890]=['θ1 Ori',1.4667255,-0.0938547,6.7];
ybs[1891]=['θ2 Ori',1.4671319,-0.0943499,5.08];
ybs[1892]=['',1.4677595,-0.0760007,6.38];
ybs[1893]=['ι Ori',1.4673426,-0.1029702,2.77];
ybs[1894]=['',1.4681342,-0.0565942,6.4];
ybs[1895]=['45 Ori',1.4683558,-0.0845732,5.26];
ybs[1896]=['',1.475853,0.4700806,5.83];
ybs[1897]=['ε Ori',1.4708856,-0.020805,1.7];
ybs[1898]=['',1.4788,0.5858785,6.33];
ybs[1899]=['122 Tau',1.4751486,0.2975751,5.54];
ybs[1900]=['',1.4709138,-0.0984045,6.54];
ybs[1901]=['φ2 Ori',1.4742159,0.1623183,4.09];
ybs[1902]=['',1.4750036,0.1927633,5.94];
ybs[1903]=['',1.4656456,-0.577174,5.78];
ybs[1904]=['ζ Tau',1.4778225,0.3691679,3];
ybs[1905]=['',1.4724066,-0.105684,5.72];
ybs[1906]=['',1.4577255,-0.9580314,6.43];
ybs[1907]=['',1.4760237,0.1564054,6.12];
ybs[1908]=['26 Aur',1.4824848,0.5323486,5.4];
ybs[1909]=['',1.4698129,-0.5008712,6.26];
ybs[1910]=['',1.5018462,1.1467662,5.6];
ybs[1911]=['',1.4533744,-1.1207817,5.34];
ybs[1912]=['',1.4761701,-0.1034797,6.05];
ybs[1913]=['',1.4746407,-0.2053559,6.11];
ybs[1914]=['',1.479021,0.1317814,5.88];
ybs[1915]=['',1.4837399,0.4647198,6.37];
ybs[1916]=['β Dor',1.4563844,-1.0904567,3.76];
ybs[1917]=['',1.4780944,-0.0838528,6.19];
ybs[1918]=['',1.4853655,0.510052,5.96];
ybs[1919]=['',1.4955726,0.9335546,6.23];
ybs[1920]=['ν1 Col',1.4746585,-0.4862817,6.16];
ybs[1921]=['',1.4683529,-0.8256074,6.11];
ybs[1922]=['125 Tau',1.4871171,0.452133,5.18];
ybs[1923]=['',1.4857199,0.3799805,6.34];
ybs[1924]=['',1.4630028,-1.0273108,6.75];
ybs[1925]=['σ Ori',1.4818979,-0.0452244,3.81];
ybs[1926]=['',1.4820653,-0.0451228,6.65];
ybs[1927]=['',1.4812706,-0.1145808,5.96];
ybs[1928]=['ω Ori',1.4840095,0.0720826,4.57];
ybs[1929]=['ν2 Col',1.4766707,-0.500563,5.31];
ybs[1930]=['',1.4623812,-1.067535,6.32];
ybs[1931]=['49 Ori',1.4823641,-0.1257383,4.8];
ybs[1932]=['',1.491087,0.5474409,6.04];
ybs[1933]=['',1.4915609,0.5572624,6.11];
ybs[1934]=['',1.4852436,-0.0620676,6];
ybs[1935]=['24 Cam',1.5032206,0.9876567,6.05];
ybs[1936]=['',1.4850303,-0.1692646,6.5];
ybs[1937]=['23 Cam',1.5086446,1.0730818,6.15];
ybs[1938]=['',1.4837201,-0.3113809,6.38];
ybs[1939]=['',1.4942926,0.5147881,6.43];
ybs[1940]=['126 Tau',1.4935972,0.2887057,4.86];
ybs[1941]=['',1.4804308,-0.7103244,5.82];
ybs[1942]=['ζ Ori',1.4906948,-0.0337686,2.05];
ybs[1943]=['ζ Ori',1.4907021,-0.0337686,4.21];
ybs[1944]=['',1.4900724,-0.0491652,6.22];
ybs[1945]=['',1.4965085,0.4072524,6.59];
ybs[1946]=['',1.4910898,-0.0195641,4.95];
ybs[1947]=['γ Men',1.4450803,-1.332193,5.19];
ybs[1948]=['',1.4971616,0.3956255,6.36];
ybs[1949]=['',1.492224,0.0060321,5.93];
ybs[1950]=['α Col',1.4847633,-0.5945579,2.64];
ybs[1951]=['',1.4904774,-0.1815396,6.52];
ybs[1952]=['',1.4856172,-0.5693393,5.45];
ybs[1953]=['',1.4946518,-0.0504141,6.42];
ybs[1954]=['',1.4700668,-1.161525,6.31];
ybs[1955]=['',1.5026787,0.4051088,6.21];
ybs[1956]=['',1.4943154,-0.2917832,6.21];
ybs[1957]=['51 Ori',1.4982937,0.0258653,4.91];
ybs[1958]=['',1.4587214,-1.2868414,5.78];
ybs[1959]=['',1.4966744,-0.3058323,6.15];
ybs[1960]=['',1.4926527,-0.5828146,6.34];
ybs[1961]=['',1.4998805,-0.1184909,6.02];
ybs[1962]=['12 Lep',1.496483,-0.3903642,5.87];
ybs[1963]=['26 Cam',1.5183625,0.9794958,5.94];
ybs[1964]=['',1.501162,-0.0280315,6.31];
ybs[1965]=['ο Aur',1.5151978,0.869734,5.47];
ybs[1966]=['',1.4960136,-0.5328165,6.19];
ybs[1967]=['',1.4960956,-0.6049377,5.29];
ybs[1968]=['',1.51431,0.707085,6.58];
ybs[1969]=['',1.5015348,-0.323769,5.73];
ybs[1970]=['',1.5303766,1.0962822,6.13];
ybs[1971]=['',1.5127665,0.3612984,6.95];
ybs[1972]=['',1.5095172,0.0700614,6.09];
ybs[1973]=['',1.5206132,0.7423203,6.29];
ybs[1974]=['',1.5063403,-0.3511595,6.34];
ybs[1975]=['',1.5013361,-0.6876607,6.25];
ybs[1976]=['',1.5061238,-0.3912193,6.15];
ybs[1977]=['γ Lep',1.5062174,-0.3916848,3.6];
ybs[1978]=['',1.5017761,-0.7998186,6.39];
ybs[1979]=['129 Tau',1.5174134,0.2762492,6];
ybs[1980]=['',1.5137003,-0.0743962,6.34];
ybs[1981]=['',1.517697,0.1662879,5.79];
ybs[1982]=['',1.5162033,0.0204825,5.95];
ybs[1983]=['131 Tau',1.5193855,0.2529602,5.72];
ybs[1984]=['130 Tau',1.5204349,0.3095217,5.49];
ybs[1985]=['ι Men',1.4597173,-1.3754985,6.05];
ybs[1986]=['29 Cam',1.5361578,0.9934864,6.54];
ybs[1987]=['133 Tau',1.5215267,0.2426833,5.29];
ybs[1988]=['',1.5481486,1.1950967,6.2];
ybs[1989]=['τ Aur',1.5288616,0.6839152,4.52];
ybs[1990]=['μ Col',1.5125428,-0.5637515,5.17];
ybs[1991]=['',1.5246265,0.3643229,6.07];
ybs[1992]=['ζ Lep',1.5173411,-0.2585979,3.55];
ybs[1993]=['52 Ori',1.5225567,0.112732,5.27];
ybs[1994]=['',1.5180449,-0.2833101,6.17];
ybs[1995]=['',1.5196127,-0.1837464,6.03];
ybs[1996]=['132 Tau',1.5275626,0.4288613,4.86];
ybs[1997]=['',1.5373157,0.8991636,6.29];
ybs[1998]=['κ Ori',1.5209989,-0.1686808,2.06];
ybs[1999]=['',1.5174013,-0.4997545,6.22];
ybs[2000]=['30 Cam',1.5438935,1.0291703,6.14];
ybs[2001]=['',1.5247664,-0.0713849,5.97];
ybs[2002]=['',1.5138211,-0.8131764,5.31];
ybs[2003]=['',1.5181205,-0.6225496,6.32];
ybs[2004]=['134 Tau',1.5294872,0.2208774,4.91];
ybs[2005]=['υ Aur',1.5369094,0.6511668,4.74];
ybs[2006]=['ν Aur',1.5389615,0.6833308,3.97];
ybs[2007]=['',1.5362041,0.4881927,5.56];
ybs[2008]=['',1.5315614,0.1723537,5.8];
ybs[2009]=['δ Dor',1.5045076,-1.1471891,4.35];
ybs[2010]=['135 Tau',1.5336116,0.2497459,5.52];
ybs[2011]=['',1.5207537,-0.7094329,6.61];
ybs[2012]=['',1.538376,0.5607464,6.25];
ybs[2013]=['',1.5321608,0.0772708,5.97];
ybs[2014]=['β Pic',1.5171287,-0.8911838,3.85];
ybs[2015]=['',1.5289288,-0.2527126,5.49];
ybs[2016]=['π Men',1.4649969,-1.4042775,5.65];
ybs[2017]=['',1.5165641,-0.948682,6.18];
ybs[2018]=['',1.5333184,0.0354003,5.98];
ybs[2019]=['',1.5440807,0.6907543,6.45];
ybs[2020]=['',1.5298838,-0.4008589,5.87];
ybs[2021]=['31 Cam',1.5556649,1.0452801,5.2];
ybs[2022]=['',1.5438624,0.5920221,5.98];
ybs[2023]=['ξ Aur',1.5547172,0.9723024,4.99];
ybs[2024]=['',1.5421191,0.3468156,6.06];
ybs[2025]=['55 Ori',1.5368156,-0.1311539,5.35];
ybs[2026]=['',1.5275019,-0.7831458,6.38];
ybs[2027]=['137 Tau',1.5418541,0.2473951,5.59];
ybs[2028]=['136 Tau',1.5464811,0.4819696,4.58];
ybs[2029]=['δ Lep',1.5361954,-0.3643486,3.81];
ybs[2030]=['',1.5368011,-0.4000756,6.17];
ybs[2031]=['56 Ori',1.5417782,0.0324286,4.78];
ybs[2032]=['',1.5461844,0.354333,6.71];
ybs[2033]=['',1.5400855,-0.1577467,5.97];
ybs[2034]=['β Col',1.5340414,-0.6242104,3.12];
ybs[2035]=['',1.5680016,1.1536061,6.25];
ybs[2036]=['γ Pic',1.5277933,-0.980219,4.51];
ybs[2037]=['',1.5388075,-0.5139183,6.45];
ybs[2038]=['',1.5309208,-0.9209023,6.35];
ybs[2039]=['',1.5604647,0.9041711,6.49];
ybs[2040]=['',1.5538706,0.5533315,5.9];
ybs[2041]=['χ1 Ori',1.5508233,0.3539227,4.41];
ybs[2042]=['',1.5498222,0.1848163,6.12];
ybs[2043]=['',1.5327622,-0.9094056,5.17];
ybs[2044]=['',1.5512255,0.2053313,6.59];
ybs[2045]=['',1.5497674,0.0563309,6.31];
ybs[2046]=['57 Ori',1.5532636,0.3447314,5.92];
ybs[2047]=['',1.5409073,-0.6567338,5.63];
ybs[2048]=['',1.5639343,0.8557419,6.47];
ybs[2049]=['',1.541919,-0.6723513,6.7];
ybs[2050]=['λ Col',1.5435304,-0.5898968,4.87];
ybs[2051]=['',1.5517657,0.0169363,6];
ybs[2052]=['',1.5509326,-0.0708912,6.57];
ybs[2053]=['',1.4263887,-1.4795449,6.2];
ybs[2054]=['',1.5477343,-0.3427116,6.69];
ybs[2055]=['α Ori',1.5538633,0.1293079,0.5];
ybs[2056]=['λ Men',1.5160332,-1.2688014,6.53];
ybs[2057]=['',1.5571032,0.3521473,5.4];
ybs[2058]=['ε Dor',1.5266287,-1.1675697,5.11];
ybs[2059]=['',1.5513463,-0.205457,5.66];
ybs[2060]=['',1.5606548,0.5051585,6.32];
ybs[2061]=['',1.5579327,0.2430675,6.6];
ybs[2062]=['',1.5486073,-0.5086845,6.36];
ybs[2063]=['',1.5442513,-0.7490725,6.55];
ybs[2064]=['',1.5549595,-0.0805411,5.87];
ybs[2065]=['',1.5553252,-0.0835427,6.28];
ybs[2066]=['',1.538647,-0.9975066,5.94];
ybs[2067]=['',1.5336236,-1.1175337,6.36];
ybs[2068]=['',1.5621022,0.4232563,6.02];
ybs[2069]=['',1.5595777,0.1659986,5.99];
ybs[2070]=['',1.5612029,0.2011012,5.87];
ybs[2071]=['δ Aur',1.5750096,0.9474453,3.72];
ybs[2072]=['',1.6038156,1.3191745,6.4];
ybs[2073]=['',1.5761276,0.9655271,6.44];
ybs[2074]=['',1.5762325,0.9520247,6.14];
ybs[2075]=['',1.5739529,0.8713456,5.89];
ybs[2076]=['',1.5509235,-0.6973637,5.57];
ybs[2077]=['',1.547256,-0.8789401,6.52];
ybs[2078]=['139 Tau',1.5667891,0.4529916,4.82];
ybs[2079]=['η Lep',1.5585945,-0.2472507,3.71];
ybs[2080]=['',1.5575804,-0.3986128,5.96];
ybs[2081]=['ξ Col',1.5537809,-0.6478496,4.97];
ybs[2082]=['β Aur',1.5743491,0.7844806,1.9];
ybs[2083]=['',1.5495195,-0.8661158,6.1];
ybs[2084]=['',1.5590361,-0.4051652,6.36];
ybs[2085]=['π Aur',1.5761832,0.8017466,4.26];
ybs[2086]=['σ Col',1.5577393,-0.5477034,5.5];
ybs[2087]=['',1.5655978,0.0213828,6.22];
ybs[2088]=['',1.5499542,-0.918622,5.29];
ybs[2089]=['θ Aur',1.5747867,0.6494782,2.62];
ybs[2090]=['',1.5777679,0.7782692,6.22];
ybs[2091]=['',1.5668045,-0.0173414,6.22];
ybs[2092]=['',1.5597584,-0.5580672,6.44];
ybs[2093]=['',1.5702386,0.223557,5.7];
ybs[2094]=['59 Ori',1.5678122,0.0320691,5.9];
ybs[2095]=['36 Aur',1.5808533,0.8360344,5.73];
ybs[2096]=['',1.5456245,-1.1010849,4.65];
ybs[2097]=['60 Ori',1.5696066,0.009658,5.22];
ybs[2098]=['',1.5458344,-1.1253844,6.63];
ybs[2099]=['',1.5841607,0.8544858,5.96];
ybs[2100]=['γ Col',1.5627583,-0.6157943,4.36];
ybs[2101]=['1 Mon',1.5701382,-0.1637465,6.12];
ybs[2102]=['2 Mon',1.5703729,-0.1668206,5.03];
ybs[2103]=['',1.5730452,-0.0252108,6.63];
ybs[2104]=['',1.580838,0.5416452,5.98];
ybs[2105]=['',1.5799901,0.4812145,6.05];
ybs[2106]=['',1.5889975,0.8709904,6.05];
ybs[2107]=['',1.5748658,-0.053658,4.53];
ybs[2108]=['',1.5603755,-0.9324424,6.45];
ybs[2109]=['',1.588941,0.7570734,6.42];
ybs[2110]=['',1.5827993,0.3909518,6.37];
ybs[2111]=['',1.5670735,-0.7685377,5.81];
ybs[2112]=['',1.5756076,-0.2251478,6.22];
ybs[2113]=['38 Aur',1.5906893,0.7489207,6.1];
ybs[2114]=['η Col',1.5694111,-0.7472633,3.96];
ybs[2115]=['',1.6000034,1.0365605,6.34];
ybs[2116]=['',1.5885323,0.5695768,6.24];
ybs[2117]=['',1.5964221,0.9000861,6.45];
ybs[2118]=['μ Ori',1.5853981,0.1683593,4.12];
ybs[2119]=['κ Men',1.5232497,-1.385041,5.47];
ybs[2120]=['',1.6071421,1.1074188,6.39];
ybs[2121]=['',1.5847303,0.0295533,6.59];
ybs[2122]=['3 Mon',1.5824201,-0.1849877,4.95];
ybs[2123]=['',1.5792177,-0.4436356,6.05];
ybs[2124]=['64 Ori',1.5903935,0.3436355,5.14];
ybs[2125]=['',1.5791227,-0.591882,5.55];
ybs[2126]=['39 Aur',1.5983578,0.7501294,5.87];
ybs[2127]=['',1.5899475,0.2038401,6.08];
ybs[2128]=['1 Gem',1.5934176,0.4059873,4.16];
ybs[2129]=['χ2 Ori',1.5924375,0.3514473,4.63];
ybs[2130]=['',1.5854549,-0.2530463,6.2];
ybs[2131]=['',1.5980419,0.6625578,6.34];
ybs[2132]=['',1.5761717,-0.8939022,5.67];
ybs[2133]=['',1.6001283,0.5863705,6.23];
ybs[2134]=['',1.5880872,-0.4587769,5.04];
ybs[2135]=['',1.6027152,0.6175785,6.12];
ybs[2136]=['',1.5929449,-0.1171315,5.21];
ybs[2137]=['40 Aur',1.6047933,0.6715979,5.36];
ybs[2138]=['63 Ori',1.5965568,0.0945565,5.67];
ybs[2139]=['66 Ori',1.5965337,0.0725412,5.63];
ybs[2140]=['',1.6034881,0.5150389,6.08];
ybs[2141]=['',1.6087315,0.7304331,6.12];
ybs[2142]=['17 Lep',1.5959579,-0.2877477,4.93];
ybs[2143]=['',1.5925699,-0.5615507,5.65];
ybs[2144]=['',1.5981827,-0.1788137,5.87];
ybs[2145]=['',1.5811488,-1.0489064,6.45];
ybs[2146]=['37 Cam',1.6210971,1.0285445,5.36];
ybs[2147]=['',1.6127688,0.7164878,6.36];
ybs[2148]=['',1.6035794,-0.0732497,5.38];
ybs[2149]=['θ Lep',1.6011127,-0.2607184,4.67];
ybs[2150]=['',1.599086,-0.4223375,6.95];
ybs[2151]=['',1.5925261,-0.7860732,6.35];
ybs[2152]=['',1.5933817,-0.7868115,5.93];
ybs[2153]=['ν Ori',1.6081954,0.2576961,4.42];
ybs[2154]=['',1.597284,-0.6198721,5.8];
ybs[2155]=['',1.6043282,-0.1950703,6.66];
ybs[2156]=['',1.593655,-0.8457995,6.58];
ybs[2157]=['',1.6024915,-0.4034065,5.47];
ybs[2158]=['',1.6003178,-0.5194334,5.81];
ybs[2159]=['36 Cam',1.6345478,1.1468995,5.32];
ybs[2160]=['',1.6043983,-0.3807592,5.78];
ybs[2161]=['',1.6133079,0.1512513,6.55];
ybs[2162]=['19 Lep',1.6076876,-0.3345668,5.31];
ybs[2163]=['',1.6170288,0.3872139,5.93];
ybs[2164]=['',1.6043623,-0.5989111,5.83];
ybs[2165]=['π1 Col',1.6023248,-0.7383016,6.12];
ybs[2166]=['',1.628276,0.9187749,6.3];
ybs[2167]=['3 Gem',1.6179045,0.4033277,5.75];
ybs[2168]=['',1.6139011,0.0435536,5.73];
ybs[2169]=['41 Aur',1.6272942,0.8501121,6.82];
ybs[2170]=['41 Aur',1.6273014,0.8500782,6.09];
ybs[2171]=['θ Col',1.6062586,-0.6502466,5.02];
ybs[2172]=['',1.6036862,-0.7870472,6.51];
ybs[2173]=['',1.6164423,-0.0997521,6.17];
ybs[2174]=['',1.6131257,-0.3915029,5.5];
ybs[2175]=['π2 Col',1.6075616,-0.7357847,5.5];
ybs[2176]=['',1.6148998,-0.3164324,6.35];
ybs[2177]=['',1.616032,-0.2546156,5.56];
ybs[2178]=['',1.6233933,0.3163329,6.33];
ybs[2179]=['5 Gem',1.6258262,0.4261247,5.8];
ybs[2180]=['',1.6167498,-0.3975596,5.71];
ybs[2181]=['',1.6104711,-0.7742257,6.27];
ybs[2182]=['',1.6368357,0.8930215,6.04];
ybs[2183]=['',1.6296209,0.5705106,5.78];
ybs[2184]=['',1.6271257,0.3815875,6.56];
ybs[2185]=['',1.6251473,0.23795,6.04];
ybs[2186]=['',1.6165577,-0.4660927,6.27];
ybs[2187]=['68 Ori',1.6277751,0.3453174,5.75];
ybs[2188]=['η1 Dor',1.5977482,-1.1526563,5.71];
ybs[2189]=['',1.6225997,-0.1179673,6.15];
ybs[2190]=['',1.6022827,-1.0848574,5.05];
ybs[2191]=['6 Gem',1.6291761,0.3997305,6.39];
ybs[2192]=['69 Ori',1.6278007,0.2814382,4.95];
ybs[2193]=['ξ Ori',1.6272378,0.2478997,4.48];
ybs[2194]=['',1.6199879,-0.4740106,5.72];
ybs[2195]=['40 Cam',1.6460494,1.047061,5.35];
ybs[2196]=['',1.6257457,-0.0815195,6.18];
ybs[2197]=['',1.6176487,-0.7043809,5.58];
ybs[2198]=['',1.6136804,-0.865105,6.49];
ybs[2199]=['',1.6262784,-0.114415,5.05];
ybs[2200]=['',1.6228419,-0.4622879,6.09];
ybs[2201]=['',1.6344356,0.3259281,6.58];
ybs[2202]=['',1.6326677,0.1853833,6.45];
ybs[2203]=['',1.6614309,1.2097106,4.8];
ybs[2204]=['',1.6302243,-0.0438085,6.62];
ybs[2205]=['',1.6195293,-0.7904,6.31];
ybs[2206]=['δ Pic',1.6172215,-0.9594608,4.81];
ybs[2207]=['',1.6298928,-0.3101214,6.52];
ybs[2208]=['',1.6384244,0.3124143,5.88];
ybs[2209]=['1 Lyn',1.6559944,1.0735054,4.98];
ybs[2210]=['η Gem',1.6403244,0.3927012,3.28];
ybs[2211]=['',1.6442347,0.6307919,6.92];
ybs[2212]=['',1.6352707,-0.0654059,5.83];
ybs[2213]=['κ Aur',1.64277,0.5147201,4.35];
ybs[2214]=['71 Ori',1.6400872,0.3342282,5.2];
ybs[2215]=['ν Dor',1.6084396,-1.2016072,5.06];
ybs[2216]=['',1.641196,0.2416316,5.91];
ybs[2217]=['72 Ori',1.6424759,0.2816315,5.3];
ybs[2218]=['',1.6383295,-0.079844,5.83];
ybs[2219]=['',1.6339619,-0.4165742,6.39];
ybs[2220]=['',1.6328984,-0.5131619,6.54];
ybs[2221]=['γ Mon',1.6393404,-0.1096279,3.98];
ybs[2222]=['42 Aur',1.6531911,0.8101195,6.52];
ybs[2223]=['73 Ori',1.6438097,0.2189378,5.33];
ybs[2224]=['8 Gem',1.6466586,0.4182305,6.08];
ybs[2225]=['',1.643256,0.105754,6.07];
ybs[2226]=['',1.6437046,0.0746427,6.64];
ybs[2227]=['',1.6426387,-0.0090588,5.65];
ybs[2228]=['',1.6421735,-0.0859059,5.99];
ybs[2229]=['',1.646778,0.2997465,6.39];
ybs[2230]=['',1.6441212,0.0202846,6.37];
ybs[2231]=['',1.6417806,-0.1578128,6.1];
ybs[2232]=['2 Lyn',1.6631523,1.0297822,4.48];
ybs[2233]=['43 Aur',1.6562555,0.809004,6.38];
ybs[2234]=['9 Gem',1.6495374,0.4142259,6.25];
ybs[2235]=['74 Ori',1.6468262,0.2140652,5.04];
ybs[2236]=['',1.6401438,-0.353932,5.91];
ybs[2237]=['',1.6408653,-0.3225949,5.99];
ybs[2238]=['',1.6430107,-0.2395498,5.01];
ybs[2239]=['η2 Dor',1.6200546,-1.1448354,5.01];
ybs[2240]=['',1.6460966,0.0187298,6.63];
ybs[2241]=['75 Ori',1.6496634,0.173399,5.39];
ybs[2242]=['',1.6489797,0.1229697,6.57];
ybs[2243]=['',1.6445767,-0.2901574,5.92];
ybs[2244]=['',1.6517324,0.2452304,6.59];
ybs[2245]=['',1.6502083,0.0888852,5.71];
ybs[2246]=['',1.6433483,-0.5200253,6.67];
ybs[2247]=['',1.6540914,0.250889,6.16];
ybs[2248]=['',1.6484309,-0.3965807,6.07];
ybs[2249]=['6 Mon',1.6511249,-0.1873248,6.75];
ybs[2250]=['κ Col',1.6457363,-0.6134434,4.37];
ybs[2251]=['4 Lyn',1.673822,1.0360716,5.94];
ybs[2252]=['',1.6582723,0.3022335,6.32];
ybs[2253]=['',1.6564646,0.1577618,6.24];
ybs[2254]=['',1.6514064,-0.2936257,5.14];
ybs[2255]=['α Men',1.6131568,-1.3047609,5.09];
ybs[2256]=['',1.6457255,-0.685419,6];
ybs[2257]=['',1.6476612,-0.658772,5.53];
ybs[2258]=['45 Aur',1.6719747,0.9327508,5.36];
ybs[2259]=['',1.6482876,-0.6503179,5.87];
ybs[2260]=['',1.65363,-0.3486268,5.52];
ybs[2261]=['',1.6566491,-0.164034,5.36];
ybs[2262]=['',1.6563436,-0.2623733,6.06];
ybs[2263]=['',1.6627243,0.2555577,5.69];
ybs[2264]=['',1.6579316,-0.1500057,6.22];
ybs[2265]=['',1.6568922,-0.365373,5.81];
ybs[2266]=['',1.6681787,0.5154285,6.43];
ybs[2267]=['7 Mon',1.6604927,-0.1366873,5.27];
ybs[2268]=['',1.6430213,-1.0335943,6.43];
ybs[2269]=['',1.66186,-0.0515417,4.9];
ybs[2270]=['',1.6661319,0.2050294,6.54];
ybs[2271]=['',1.6687652,0.3098709,6.35];
ybs[2272]=['',1.6504179,-0.9204992,6.41];
ybs[2273]=['',1.6594347,-0.6004833,5.78];
ybs[2274]=['',1.6682827,0.0394326,6.31];
ybs[2275]=['',1.6546073,-0.8790741,7.04];
ybs[2276]=['ζ CMa',1.6623631,-0.5248572,3.02];
ybs[2277]=['',1.6354562,-1.2515605,6.64];
ybs[2278]=['',1.6677829,-0.205645,5.64];
ybs[2279]=['',1.7027153,1.2308614,5.97];
ybs[2280]=['μ Gem',1.6755927,0.3927626,2.88];
ybs[2281]=['',1.6737344,0.2192168,6];
ybs[2282]=['',1.6634675,-0.5960783,5.53];
ybs[2283]=['ψ1 Aur',1.6853128,0.8600496,4.91];
ybs[2284]=['',1.6604889,-0.8508435,6.6];
ybs[2285]=['',1.6925623,0.9821575,5.64];
ybs[2286]=['',1.6765227,0.065526,6.4];
ybs[2287]=['5 Lyn',1.6944576,1.0193688,5.21];
ybs[2288]=['β CMa',1.6732089,-0.3135594,1.98];
ybs[2289]=['',1.6765782,-0.0819838,6.67];
ybs[2290]=['δ Col',1.6700798,-0.5837413,3.85];
ybs[2291]=['',1.6842188,0.5183005,6.71];
ybs[2292]=['ε Mon',1.6785544,0.0799796,4.44];
ybs[2293]=['',1.6785836,0.0800281,6.72];
ybs[2294]=['',1.679861,0.1548861,6.26];
ybs[2295]=['',1.6773898,-0.1725195,6.19];
ybs[2296]=['',1.6837629,0.2800585,6.33];
ybs[2297]=['',1.677956,-0.2632291,6.24];
ybs[2298]=['',1.6869138,0.4069389,6.06];
ybs[2299]=['',1.6798335,-0.2014233,5.22];
ybs[2300]=['',1.6779193,-0.3454972,6.6];
ybs[2301]=['',1.675065,-0.5550148,6.34];
ybs[2302]=['',1.6862802,0.2567541,6.24];
ybs[2303]=['',1.6805311,-0.2264215,6.12];
ybs[2304]=['',1.6849555,0.1234809,5.98];
ybs[2305]=['',1.6783219,-0.4465963,5.63];
ybs[2306]=['',1.6851751,0.0260086,6.66];
ybs[2307]=['',1.6849574,-0.0167032,5.87];
ybs[2308]=['',1.6980705,0.8271672,6.56];
ybs[2309]=['',1.6872487,0.0394636,6.51];
ybs[2310]=['',1.6782414,-0.6408517,5.62];
ybs[2311]=['',1.6871031,-0.0680729,6.35];
ybs[2312]=['',1.6816934,-0.5024915,6.39];
ybs[2313]=['',1.6961787,0.5681241,6.43];
ybs[2314]=['ν Pic',1.6722168,-0.9840131,5.61];
ybs[2315]=['',1.6878346,-0.1379747,6.4];
ybs[2316]=['',1.6756286,-0.9109087,5.98];
ybs[2317]=['',1.6812328,-0.7032766,6.31];
ybs[2318]=['',1.6909914,-0.0265066,5.87];
ybs[2319]=['',1.6905292,-0.0804366,6.15];
ybs[2320]=['α Car',1.6769989,-0.9198946,-0.72];
ybs[2321]=['',1.6924496,0.0144723,6.71];
ybs[2322]=['',1.6911986,-0.1312996,6.27];
ybs[2323]=['',1.6847811,-0.6121714,6.25];
ybs[2324]=['16 Gem',1.6972661,0.357514,6.22];
ybs[2325]=['6 Lyn',1.711754,1.0148987,5.88];
ybs[2326]=['48 Aur',1.7003557,0.5319886,5.55];
ybs[2327]=['',1.6940887,0.0505494,5.55];
ybs[2328]=['',1.6935318,0.0050165,5.2];
ybs[2329]=['',1.6936457,-0.0050241,5.55];
ybs[2330]=['',1.6757118,-1.0219607,6.48];
ybs[2331]=['',1.6717395,-1.1116499,6.27];
ybs[2332]=['47 Aur',1.7076147,0.8145896,5.9];
ybs[2333]=['',1.7018567,0.4704583,6.47];
ybs[2334]=['',1.6994102,0.283198,6.23];
ybs[2335]=['',1.6807069,-0.9218305,6.51];
ybs[2336]=['',1.6985565,0.1796237,6.15];
ybs[2337]=['ν Gem',1.7017069,0.3525517,4.15];
ybs[2338]=['10 Mon',1.6965676,-0.0833266,5.06];
ybs[2339]=['',1.6774632,-1.0522841,5.8];
ybs[2340]=['',1.7595677,1.3889642,6.54];
ybs[2341]=['',1.6981608,0.033162,6.48];
ybs[2342]=['',1.6850751,-0.8410434,5.76];
ybs[2343]=['',1.6925156,-0.4514879,6.07];
ybs[2344]=['',1.7809057,1.4328406,6.65];
ybs[2345]=['',1.7015744,0.1921074,6.59];
ybs[2346]=['π1 Dor',1.6688344,-1.2216212,5.56];
ybs[2347]=['',1.6917507,-0.6616052,6.48];
ybs[2348]=['',1.6779609,-1.1072244,6.46];
ybs[2349]=['',1.7024076,0.0459636,6.16];
ybs[2350]=['β Mon',1.7002381,-0.1229615,4.6];
ybs[2351]=['β Mon',1.7002743,-0.1229906,5.4];
ybs[2352]=['β Mon',1.7002743,-0.1229906,5.6];
ybs[2353]=['',1.6990644,-0.3050558,5.77];
ybs[2354]=['',1.6800327,-1.1141899,6.27];
ybs[2355]=['λ CMa',1.6965389,-0.5688392,4.48];
ybs[2356]=['',1.7062771,0.1573625,6.57];
ybs[2357]=['',1.7591954,1.3609755,5.73];
ybs[2358]=['',1.698664,-0.5651969,5.74];
ybs[2359]=['',1.7458715,1.285942,6.24];
ybs[2360]=['',1.71122,0.2954002,6.2];
ybs[2361]=['',1.7061315,-0.1761798,5.93];
ybs[2362]=['',1.6984885,-0.7170988,6.32];
ybs[2363]=['',1.6901074,-1.0125312,5.82];
ybs[2364]=['',1.7109985,0.1961301,6.14];
ybs[2365]=['19 Gem',1.7131721,0.2773279,6.4];
ybs[2366]=['',1.7173841,0.5661976,5.87];
ybs[2367]=['',1.7077392,-0.2297061,6.16];
ybs[2368]=['',1.7131748,0.2055754,6.65];
ybs[2369]=['',1.713829,0.2012498,5.23];
ybs[2370]=['7 Lyn',1.727879,0.9658324,6.45];
ybs[2371]=['π2 Dor',1.6813052,-1.2165129,5.38];
ybs[2372]=['',1.7163781,0.2034998,6.03];
ybs[2373]=['',1.7112753,-0.2165105,5.15];
ybs[2374]=['',1.7080892,-0.4848984,5.93];
ybs[2375]=['',1.7133746,-0.1426234,5.43];
ybs[2376]=['12 Mon',1.7158752,0.0845079,5.84];
ybs[2377]=['',1.7228974,0.576127,6.42];
ybs[2378]=['',1.7027902,-0.8770611,5.27];
ybs[2379]=['13 Mon',1.718494,0.1277392,4.5];
ybs[2380]=['',1.7158428,-0.102674,5.6];
ybs[2381]=['ξ1 CMa',1.712974,-0.4089653,4.33];
ybs[2382]=['',1.7097271,-0.6156217,5.84];
ybs[2383]=['',1.7007361,-0.9924874,5.22];
ybs[2384]=['',1.7084979,-0.714357,6.2];
ybs[2385]=['',1.7217487,0.2468041,5.53];
ybs[2386]=['',1.7174144,-0.1951357,6.24];
ybs[2387]=['',1.7112516,-0.6449604,6.34];
ybs[2388]=['8 Lyn',1.7422478,1.0727635,5.94];
ybs[2389]=['',1.7214197,-0.0215498,5.1];
ybs[2390]=['',1.7566536,1.2519476,5.92];
ybs[2391]=['',1.7161192,-0.5592825,5.69];
ybs[2392]=['49 Aur',1.7291998,0.4888158,5.27];
ybs[2393]=['',1.714584,-0.6581725,5.24];
ybs[2394]=['',1.7091607,-0.9047693,5.6];
ybs[2395]=['',1.7854728,1.3883149,5.45];
ybs[2396]=['11 Lyn',1.7415196,0.9920671,5.85];
ybs[2397]=['',1.7199906,-0.3654408,6.4];
ybs[2398]=['14 Mon',1.7266603,0.1319046,6.45];
ybs[2399]=['',1.7355139,0.6707223,5.29];
ybs[2400]=['',1.7290016,0.1740649,5.88];
ybs[2401]=['',1.718104,-0.6743857,6.44];
ybs[2402]=['',1.7021157,-1.1446058,6.29];
ybs[2403]=['',1.728601,0.0152695,5.8];
ybs[2404]=['',1.7075905,-1.0802361,6.15];
ybs[2405]=['',1.7210667,-0.6326239,5.42];
ybs[2406]=['μ Pic',1.7114438,-1.0256911,5.7];
ybs[2407]=['',1.7319213,0.0782268,6.55];
ybs[2408]=['ξ2 CMa',1.7269612,-0.4010718,4.54];
ybs[2409]=['',1.7245499,-0.5712667,5.62];
ybs[2410]=['',1.7184182,-0.9135601,6.19];
ybs[2411]=['',1.7388981,0.4289104,6.44];
ybs[2412]=['',1.7342032,-0.0912245,5.52];
ybs[2413]=['51 Aur',1.7447641,0.6872098,5.69];
ybs[2414]=['ψ3 Aur',1.7454945,0.6961389,5.2];
ybs[2415]=['γ Gem',1.7397513,0.2859373,1.93];
ybs[2416]=['',1.7380855,0.1068011,6.06];
ybs[2417]=['ν1 CMa',1.7328804,-0.32595,5.7];
ybs[2418]=['',1.7279226,-0.642196,5.59];
ybs[2419]=['53 Aur',1.7431168,0.5055816,5.79];
ybs[2420]=['',1.7391575,0.1891402,6.38];
ybs[2421]=['ψ2 Aur',1.7478601,0.741276,4.79];
ybs[2422]=['',1.734788,-0.2327671,5.97];
ybs[2423]=['ν2 CMa',1.7341846,-0.3363515,3.95];
ybs[2424]=['',1.7391625,0.046915,6.17];
ybs[2425]=['',1.7301287,-0.6301374,6.35];
ybs[2426]=['',1.7401304,0.0862366,6.15];
ybs[2427]=['',1.7340644,-0.3949799,6.35];
ybs[2428]=['',1.750714,0.7678876,6.41];
ybs[2429]=['',1.7250956,-0.9248578,4.39];
ybs[2430]=['',1.7459425,0.3842181,6.04];
ybs[2431]=['',1.7387402,-0.2269123,6.12];
ybs[2432]=['54 Aur',1.7481864,0.4929874,6.03];
ybs[2433]=['',1.7479293,0.4290553,6.38];
ybs[2434]=['',1.7419322,-0.0446808,6.14];
ybs[2435]=['',1.7442533,0.0817501,6.57];
ybs[2436]=['',1.7433338,0.0278742,6.21];
ybs[2437]=['ν3 CMa',1.7394833,-0.3185871,4.43];
ybs[2438]=['',1.7349728,-0.6660606,6.04];
ybs[2439]=['',1.734032,-0.7255798,6.34];
ybs[2440]=['',1.7358921,-0.6458842,5.71];
ybs[2441]=['',1.7385523,-0.5647161,5.27];
ybs[2442]=['',1.7425819,-0.2947879,6.03];
ybs[2443]=['',1.7487268,0.2262997,5.97];
ybs[2444]=['',1.7456715,-0.2471843,4.82];
ybs[2445]=['ν Pup',1.7379011,-0.7541953,3.17];
ybs[2446]=['',1.7575588,0.6268192,6.46];
ybs[2447]=['25 Gem',1.7560293,0.4918106,6.42];
ybs[2448]=['',1.7517419,0.110904,6.51];
ybs[2449]=['',1.746784,-0.4138606,6.05];
ybs[2450]=['15 Mon',1.753803,0.1724041,4.66];
ybs[2451]=['',1.7556832,0.2858815,6.28];
ybs[2452]=['',1.7551818,0.1917362,6.11];
ybs[2453]=['ψ4 Aur',1.7643379,0.7767758,5.02];
ybs[2454]=['',1.7469992,-0.5321023,5.71];
ybs[2455]=['',1.754013,0.0083377,5.79];
ybs[2456]=['',1.7413741,-0.8418895,4.93];
ybs[2457]=['',1.7697762,0.9298665,6.27];
ybs[2458]=['',1.7646006,0.6480133,6.19];
ybs[2459]=['',1.7477013,-0.6662954,6.58];
ybs[2460]=['26 Gem',1.7602677,0.3076517,5.21];
ybs[2461]=['',1.7581041,0.1104279,6.37];
ybs[2462]=['',1.7374658,-1.0742358,6.18];
ybs[2463]=['',1.7574376,-0.1603107,5.19];
ybs[2464]=['12 Lyn',1.7792667,1.0371065,4.87];
ybs[2465]=['',1.7688287,0.6299032,6.31];
ybs[2466]=['ε Gem',1.7671803,0.4382927,2.98];
ybs[2467]=['',1.7628936,0.0526204,6.19];
ybs[2468]=['',1.7532101,-0.7045419,6.12];
ybs[2469]=['',1.7509901,-0.8323839,6.65];
ybs[2470]=['13 Lyn',1.7816003,0.9974398,5.35];
ybs[2471]=['30 Gem',1.767037,0.2305403,4.49];
ybs[2472]=['',1.7652616,0.068305,5.9];
ybs[2473]=['28 Gem',1.7709276,0.5053024,5.44];
ybs[2474]=['',1.7606096,-0.39213,6.13];
ybs[2475]=['',1.7578371,-0.6704962,6.29];
ybs[2476]=['ψ5 Aur',1.7802274,0.7602219,5.25];
ybs[2477]=['ξ Gem',1.772706,0.2247325,3.36];
ybs[2478]=['',1.7875231,0.9718656,6.33];
ybs[2479]=['',1.7874794,0.9718657,6.28];
ybs[2480]=['ψ6 Aur',1.7845638,0.8511805,5.22];
ybs[2481]=['',1.7626653,-0.6843748,6.3];
ybs[2482]=['32 Gem',1.7753758,0.2212034,6.46];
ybs[2483]=['42 Cam',1.8010411,1.1789714,5.14];
ybs[2484]=['α CMa',1.7712057,-0.2920869,-1.46];
ybs[2485]=['10 CMa',1.767756,-0.5426138,5.2];
ybs[2486]=['',1.7696078,-0.4775306,6.45];
ybs[2487]=['16 Mon',1.77803,0.1495289,5.93];
ybs[2488]=['',1.7720301,-0.4098254,6.05];
ybs[2489]=['',1.7702481,-0.5341628,6.54];
ybs[2490]=['',1.7749315,-0.2585827,5.32];
ybs[2491]=['',1.7820414,0.3171808,6.2];
ybs[2492]=['',1.7716839,-0.5552401,5.92];
ybs[2493]=['',1.7723359,-0.540498,5.8];
ybs[2494]=['',1.7779619,-0.176751,5.66];
ybs[2495]=['17 Mon',1.7814605,0.1399239,4.77];
ybs[2496]=['11 CMa',1.7787103,-0.2521264,5.29];
ybs[2497]=['',1.7483608,-1.2530211,6.51];
ybs[2498]=['18 Mon',1.7836083,0.0417454,4.47];
ybs[2499]=['',1.7742871,-0.6904447,6.62];
ybs[2500]=['',1.7822203,-0.1574043,5.07];
ybs[2501]=['12 CMa',1.7792554,-0.3671398,6.08];
ybs[2502]=['',1.7750114,-0.6596505,6.21];
ybs[2503]=['43 Cam',1.8133118,1.2019266,5.12];
ybs[2504]=['',1.7925774,0.5687238,5.71];
ybs[2505]=['',1.7708077,-0.9114177,6.57];
ybs[2506]=['',1.7854925,-0.0233829,5.75];
ybs[2507]=['',1.7727994,-0.9150668,5.8];
ybs[2508]=['ψ7 Aur',1.7976947,0.7288447,5.02];
ybs[2509]=['',1.7888105,0.0171227,6.15];
ybs[2510]=['',1.7800452,-0.6623495,5.26];
ybs[2511]=['33 Gem',1.7926161,0.2824213,5.85];
ybs[2512]=['14 Lyn',1.8091372,1.0371781,5.33];
ybs[2513]=['',1.7896391,-0.0400189,5.74];
ybs[2514]=['',1.7879028,-0.2646887,5.39];
ybs[2515]=['',1.7771864,-0.8951045,5.4];
ybs[2516]=['',1.7760887,-0.9549483,6.46];
ybs[2517]=['35 Gem',1.7951256,0.2337321,5.65];
ybs[2518]=['',1.7787201,-0.9697057,5.61];
ybs[2519]=['',1.8438835,1.3430612,4.55];
ybs[2520]=['',1.7909703,-0.4205713,6.33];
ybs[2521]=['36 Gem',1.8003001,0.3794199,5.27];
ybs[2522]=['',1.7964888,-0.0098166,5.77];
ybs[2523]=['',1.7595029,-1.2764711,6.37];
ybs[2524]=['',1.8081626,0.7822006,6.26];
ybs[2525]=['',1.8023258,0.4115404,5.65];
ybs[2526]=['',1.7957175,-0.1407201,6.29];
ybs[2527]=['',1.7939584,-0.2985436,5.79];
ybs[2528]=['',1.7660803,-1.2296339,6.11];
ybs[2529]=['',1.7924542,-0.4774377,7.04];
ybs[2530]=['κ CMa',1.7911269,-0.5677516,3.96];
ybs[2531]=['59 Aur',1.8074077,0.6780006,6.12];
ybs[2532]=['θ Gem',1.8061563,0.5923408,3.6];
ybs[2533]=['60 Aur',1.8082524,0.6704749,6.3];
ybs[2534]=['',1.8074052,0.6242347,6.01];
ybs[2535]=['',1.800185,0.0527038,6.38];
ybs[2536]=['',1.794758,-0.450287,6.33];
ybs[2537]=['',1.7935451,-0.5537493,5.7];
ybs[2538]=['',1.7910293,-0.7936217,6.55];
ybs[2539]=['ψ8 Aur',1.8114253,0.6716381,6.48];
ybs[2540]=['',1.7907265,-0.8139496,5.14];
ybs[2541]=['',1.7955522,-0.6001979,4.99];
ybs[2542]=['α Pic',1.7818574,-1.0814368,3.27];
ybs[2543]=['',1.8054385,0.1458715,5.77];
ybs[2544]=['',1.8031136,-0.093172,6.3];
ybs[2545]=['τ Pup',1.790577,-0.8837628,2.93];
ybs[2546]=['',1.789989,-0.9362527,4.4];
ybs[2547]=['',1.8079232,0.1915272,6.24];
ybs[2548]=['',1.8174957,0.7994108,6.34];
ybs[2549]=['',1.8173495,0.7659636,6.13];
ybs[2550]=['',1.799101,-0.6327202,5.96];
ybs[2551]=['ζ Men',1.7390691,-1.4107538,5.64];
ybs[2552]=['15 Lyn',1.8272962,1.019239,4.35];
ybs[2553]=['',1.8269685,1.0042442,6.05];
ybs[2554]=['',1.7900506,-1.0519153,6.11];
ybs[2555]=['',1.7977644,-0.843244,6.42];
ybs[2556]=['38 Gem',1.8135299,0.2295904,4.65];
ybs[2557]=['',1.8067622,-0.3325794,5.64];
ybs[2558]=['',1.8069763,-0.3308441,6.14];
ybs[2559]=['',1.8051313,-0.4708894,6.4];
ybs[2560]=['ψ9 Aur',1.8231268,0.8072167,5.87];
ybs[2561]=['37 Gem',1.8168275,0.4424767,5.73];
ybs[2562]=['',1.8107777,-0.1025465,6.41];
ybs[2563]=['15 CMa',1.8077489,-0.3533747,4.83];
ybs[2564]=['',1.8120854,-0.020072,5.45];
ybs[2565]=['',1.8248758,0.8147332,5.86];
ybs[2566]=['θ CMa',1.8108056,-0.2105147,4.07];
ybs[2567]=['',1.8029721,-0.7422318,6.52];
ybs[2568]=['',1.8075315,-0.4985084,6.04];
ybs[2569]=['',1.8133399,-0.03106,6.21];
ybs[2570]=['',1.8092324,-0.4286883,6.21];
ybs[2571]=['',1.8034302,-0.7679132,6.46];
ybs[2572]=['ο1 CMa',1.8101613,-0.422489,3.87];
ybs[2573]=['',1.8471382,1.2353768,5.68];
ybs[2574]=['',1.8145241,-0.0493394,6.04];
ybs[2575]=['',1.8105411,-0.4180294,6.91];
ybs[2576]=['',1.8174499,0.144882,6.29];
ybs[2577]=['16 Lyn',1.8277677,0.7866141,4.9];
ybs[2578]=['',1.8245528,0.5874235,5.89];
ybs[2579]=['',1.8027266,-0.9444381,6.57];
ybs[2580]=['17 CMa',1.8142744,-0.3565368,5.74];
ybs[2581]=['',1.8212212,0.1733539,5.92];
ybs[2582]=['π CMa',1.816807,-0.3518577,4.68];
ybs[2583]=['',1.8107636,-0.7398206,6.32];
ybs[2584]=['',1.8021104,-1.0360867,6.41];
ybs[2585]=['μ CMa',1.8191231,-0.2455224,5];
ybs[2586]=['',1.8084824,-0.883739,6.26];
ybs[2587]=['',1.8173974,-0.4008152,5.3];
ybs[2588]=['ι CMa',1.819146,-0.2980666,4.37];
ybs[2589]=['',1.825635,0.2073999,6.27];
ybs[2590]=['',1.817655,-0.5552533,6.36];
ybs[2591]=['',1.8231894,-0.1431701,6.34];
ybs[2592]=['62 Aur',1.8336351,0.6636991,6];
ybs[2593]=['39 Gem',1.8320338,0.454766,6.1];
ybs[2594]=['ι Vol',1.7943948,-1.2389216,5.4];
ybs[2595]=['',1.8238229,-0.3879443,6.61];
ybs[2596]=['',1.8212102,-0.6172475,6.29];
ybs[2597]=['40 Gem',1.8349724,0.4518474,6.4];
ybs[2598]=['',1.8308384,0.1325942,6.27];
ybs[2599]=['',1.8251292,-0.4303095,5.46];
ybs[2600]=['',1.8183403,-0.8507588,4.95];
ybs[2601]=['',2.0423917,1.5159853,5.07];
ybs[2602]=['',1.8320308,0.0624347,5.97];
ybs[2603]=['',1.825644,-0.4810462,6.23];
ybs[2604]=['',1.8235236,-0.6201457,6.23];
ybs[2605]=['',1.8338182,0.1272661,6.35];
ybs[2606]=['',1.8274905,-0.474543,6.37];
ybs[2607]=['41 Gem',1.8381308,0.280184,5.68];
ybs[2608]=['',1.8296115,-0.4439886,5.59];
ybs[2609]=['',1.8667481,1.2340163,6.5];
ybs[2610]=['ε CMa',1.8295986,-0.5060934,1.5];
ybs[2611]=['',1.8284932,-0.595792,5.06];
ybs[2612]=['',1.8431635,0.5652854,6.59];
ybs[2613]=['',1.8299794,-0.5414467,6.42];
ybs[2614]=['',1.8376637,-0.0941162,6.3];
ybs[2615]=['',1.8343531,-0.3774896,6.26];
ybs[2616]=['',1.8379898,-0.1471748,5.96];
ybs[2617]=['',1.8364953,-0.3522828,6.31];
ybs[2618]=['',1.8291325,-0.7992359,6.22];
ybs[2619]=['',1.8391011,-0.1610714,6.49];
ybs[2620]=['',1.8372453,-0.3865022,6.53];
ybs[2621]=['',1.8440209,0.0836354,6.63];
ybs[2622]=['ω Gem',1.847775,0.4221751,5.18];
ybs[2623]=['',1.847614,0.3094318,5.94];
ybs[2624]=['',1.8469498,0.2672055,5.74];
ybs[2625]=['',1.8450313,0.0965395,6.59];
ybs[2626]=['',1.8282476,-0.973094,6.27];
ybs[2627]=['',1.8481618,0.2905571,5.82];
ybs[2628]=['',1.8446807,-0.0239412,6.17];
ybs[2629]=['',1.8387164,-0.4976823,6.27];
ybs[2630]=['',1.8279388,-0.9847049,6.45];
ybs[2631]=['',1.8448096,-0.1003287,5.2];
ybs[2632]=['',1.8405291,-0.4405401,5.63];
ybs[2633]=['',1.839046,-0.5845277,6.4];
ybs[2634]=['',1.8655808,1.0432522,6.44];
ybs[2635]=['',1.852728,0.511562,5.93];
ybs[2636]=['',1.8633629,0.9203163,6.12];
ybs[2637]=['',1.8607764,0.8333495,6.38];
ybs[2638]=['σ CMa',1.8431376,-0.4880078,3.47];
ybs[2639]=['',1.8511635,0.1590269,5.97];
ybs[2640]=['19 Mon',1.8490994,-0.0744515,4.99];
ybs[2641]=['',1.8526785,0.1906731,5.13];
ybs[2642]=['ζ Gem',1.8550477,0.3585459,3.79];
ybs[2643]=['',1.8537162,0.2193433,5.98];
ybs[2644]=['',1.8381984,-0.8975907,5.14];
ybs[2645]=['ο2 CMa',1.8489852,-0.4164345,3.02];
ybs[2646]=['',1.8554735,0.025502,6.57];
ybs[2647]=['',1.8541881,-0.0933869,5.62];
ybs[2648]=['',1.8534741,-0.1771713,6.45];
ybs[2649]=['γ CMa',1.8524503,-0.2733229,4.12];
ybs[2650]=['',1.8448208,-0.7580038,6.43];
ybs[2651]=['44 Gem',1.8603435,0.3946123,6.02];
ybs[2652]=['',1.8646505,0.6011945,5.55];
ybs[2653]=['',1.8385189,-1.029146,6.02];
ybs[2654]=['',1.8317668,-1.1857987,5.17];
ybs[2655]=['',1.8614325,0.1598391,5.78];
ybs[2656]=['',1.8567267,-0.3850118,6.09];
ybs[2657]=['',1.8697709,0.59308,5.91];
ybs[2658]=['',1.8526767,-0.7393949,5.2];
ybs[2659]=['',1.8522089,-0.7615744,5.54];
ybs[2660]=['',1.8523177,-0.7616377,6.79];
ybs[2661]=['',1.8697411,0.4912884,6.48];
ybs[2662]=['',1.8616329,-0.1865562,6.49];
ybs[2663]=['',1.8692887,0.3957565,7.68];
ybs[2664]=['',1.8515843,-0.8658717,4.93];
ybs[2665]=['',1.8734712,0.5899809,6.28];
ybs[2666]=['',1.8479905,-1.0333163,5.5];
ybs[2667]=['',1.8752983,0.6530329,6.16];
ybs[2668]=['',1.8676573,0.0852064,6.11];
ybs[2669]=['',1.8595199,-0.6074686,6.14];
ybs[2670]=['',1.8653233,-0.1976111,5.39];
ybs[2671]=['',1.8649417,-0.2168043,6.48];
ybs[2672]=['',1.8617663,-0.5355258,6.34];
ybs[2673]=['',1.9023555,1.2528909,6.35];
ybs[2674]=['',1.870865,0.1298962,5.75];
ybs[2675]=['',1.8528004,-0.9909415,5.17];
ybs[2676]=['45 Gem',1.8734789,0.2775422,5.44];
ybs[2677]=['',1.8615649,-0.6703912,6.11];
ybs[2678]=['',1.8657343,-0.4361356,6.08];
ybs[2679]=['',1.8575554,-0.8794321,6.46];
ybs[2680]=['',1.8662421,-0.4657586,6.62];
ybs[2681]=['θ Men',1.812714,-1.3865567,5.45];
ybs[2682]=['',1.8679767,-0.4165866,5.71];
ybs[2683]=['',1.8661587,-0.7142162,5.79];
ybs[2684]=['',1.8812592,0.3703133,6.43];
ybs[2685]=['δ CMa',1.8723181,-0.4611531,1.84];
ybs[2686]=['',1.8769588,-0.1811027,6.21];
ybs[2687]=['',1.8742762,-0.4201555,6.65];
ybs[2688]=['63 Aur',1.8886872,0.6857459,4.9];
ybs[2689]=['τ Gem',1.8860594,0.5273564,4.41];
ybs[2690]=['',1.8659381,-0.9075019,5.96];
ybs[2691]=['',1.8777208,-0.2838604,6.03];
ybs[2692]=['47 Gem',1.8870041,0.4682123,5.78];
ybs[2693]=['20 Mon',1.8810182,-0.0744695,4.92];
ybs[2694]=['',1.8737854,-0.6926302,4.83];
ybs[2695]=['',1.8969549,0.8970633,5.47];
ybs[2696]=['',1.8781343,-0.4408779,5.69];
ybs[2697]=['',1.8802703,-0.3266348,6.23];
ybs[2698]=['48 Gem',1.8915119,0.4205866,5.85];
ybs[2699]=['21 Mon',1.8862131,-0.0057942,5.45];
ybs[2700]=['',1.8807097,-0.4803314,5.46];
ybs[2701]=['',1.9574341,1.4175809,6.31];
ybs[2702]=['',1.888396,0.0981658,6.09];
ybs[2703]=['',1.8932624,0.4746259,6.43];
ybs[2704]=['',1.859501,-1.2019206,6.47];
ybs[2705]=['',1.8895616,0.0950223,6.16];
ybs[2706]=['δ Mon',1.8882657,-0.0091283,4.15];
ybs[2707]=['18 Lyn',1.9086915,1.0403123,5.2];
ybs[2708]=['',1.8869136,-0.3650041,5.84];
ybs[2709]=['51 Gem',1.895321,0.2814871,5];
ybs[2710]=['26 CMa',1.8889696,-0.4533114,5.92];
ybs[2711]=['',1.8817372,-0.854547,5.14];
ybs[2712]=['',1.8882024,-0.5384679,6.1];
ybs[2713]=['',1.9073389,0.8239365,5.58];
ybs[2714]=['',1.9002706,0.430739,6.89];
ybs[2715]=['',1.8934345,-0.1969101,5.78];
ybs[2716]=['',1.8897811,-0.4800455,6.59];
ybs[2717]=['52 Gem',1.9013889,0.433777,5.82];
ybs[2718]=['',1.88954,-0.6383516,5.96];
ybs[2719]=['',1.8886373,-0.7073684,5.31];
ybs[2720]=['',1.9003117,0.2109144,5.62];
ybs[2721]=['',1.8991373,0.0537589,5.35];
ybs[2722]=['',1.8942961,-0.3962622,6.01];
ybs[2723]=['',1.8982854,-0.068636,5.75];
ybs[2724]=['',1.8984399,-0.1741609,5.9];
ybs[2725]=['',1.8960632,-0.4003327,6.36];
ybs[2726]=['',1.8950373,-0.4779984,6.12];
ybs[2727]=['γ1 Vol',1.8698935,-1.2309102,5.69];
ybs[2728]=['γ2 Vol',1.8700897,-1.2309396,3.78];
ybs[2729]=['',1.9151957,0.9092854,5.92];
ybs[2730]=['53 Gem',1.9069615,0.4863461,5.71];
ybs[2731]=['',1.8993529,-0.1806056,6.03];
ybs[2732]=['',1.889594,-0.8166375,4.49];
ybs[2733]=['',1.895685,-0.5430566,6.6];
ybs[2734]=['',1.9837924,1.4376798,4.96];
ybs[2735]=['',1.8964481,-0.5300745,6.33];
ybs[2736]=['24 Mon',1.9033576,-0.0033688,6.41];
ybs[2737]=['27 CMa',1.8979025,-0.4604817,4.66];
ybs[2738]=['',1.8925819,-0.789129,4.89];
ybs[2739]=['',1.9050491,0.138684,5.82];
ybs[2740]=['',1.8939949,-0.7796483,5.1];
ybs[2741]=['ω CMa',1.9003245,-0.4678209,3.85];
ybs[2742]=['',1.9004897,-0.4724511,5.58];
ybs[2743]=['',1.9192388,0.8627511,5.05];
ybs[2744]=['',1.9047926,-0.1852782,5.95];
ybs[2745]=['64 Aur',1.9165945,0.7129767,5.78];
ybs[2746]=['',1.8858128,-1.1033999,6.02];
ybs[2747]=['',1.9047059,-0.4149056,6.32];
ybs[2748]=['',1.902531,-0.5361299,5.36];
ybs[2749]=['',1.9163071,0.5397092,6.24];
ybs[2750]=['',1.9069284,-0.2725822,5.46];
ybs[2751]=['',1.9003271,-0.7235655,5.94];
ybs[2752]=['',1.9121676,0.1160317,6.65];
ybs[2753]=['',1.8992151,-0.8182286,5.72];
ybs[2754]=['',1.8985669,-0.8430452,4.76];
ybs[2755]=['λ Gem',1.9159314,0.2881105,3.58];
ybs[2756]=['',1.9083062,-0.4074937,4.79];
ybs[2757]=['',1.9128028,-0.1171553,6.29];
ybs[2758]=['',1.9080177,-0.4871774,4.64];
ybs[2759]=['',1.9013883,-0.9168435,5.97];
ybs[2760]=['',1.9095213,-0.5398113,6.32];
ybs[2761]=['',1.9073652,-0.6693502,5.8];
ybs[2762]=['',1.9087297,-0.6392261,5.03];
ybs[2763]=['',1.9057394,-0.8169253,5.66];
ybs[2764]=['47 Cam',1.9364943,1.0448834,6.35];
ybs[2765]=['π Pup',1.9100972,-0.6480373,2.7];
ybs[2766]=['',1.9133591,-0.4682733,6.46];
ybs[2767]=['',1.9298318,0.7438867,6.35];
ybs[2768]=['',1.9310163,0.7887835,5.77];
ybs[2769]=['δ Gem',1.9249576,0.3830761,3.53];
ybs[2770]=['',1.9211174,0.0472513,5.89];
ybs[2771]=['',1.9230736,0.1240816,5.91];
ybs[2772]=['',1.924709,0.2637054,6.45];
ybs[2773]=['29 CMa',1.9172479,-0.4292084,4.98];
ybs[2774]=['τ CMa',1.9173872,-0.4361075,4.4];
ybs[2775]=['19 Lyn',1.9384954,0.9642893,6.53];
ybs[2776]=['19 Lyn',1.9385823,0.9642358,5.45];
ybs[2777]=['',1.9189978,-0.3370818,6.09];
ybs[2778]=['',1.9179712,-0.4645864,5.28];
ybs[2779]=['',1.9151913,-0.6417041,4.66];
ybs[2780]=['',1.9209851,-0.2867273,5.7];
ybs[2781]=['',1.9138073,-0.7682822,5.85];
ybs[2782]=['',1.9166314,-0.6418566,5.11];
ybs[2783]=['',1.9161857,-0.684922,5.25];
ybs[2784]=['',1.9347639,0.6800093,6.4];
ybs[2785]=['65 Aur',1.9338793,0.6409928,5.13];
ybs[2786]=['',1.9193404,-0.5892296,6.3];
ybs[2787]=['56 Gem',1.9328623,0.3562096,5.1];
ybs[2788]=['',1.9275994,-0.2512204,5.45];
ybs[2789]=['',1.9974186,1.4112201,6.41];
ybs[2790]=['',1.9291186,-0.1555495,6.55];
ybs[2791]=['',1.9269719,-0.3994271,6.61];
ybs[2792]=['',1.9269636,-0.4711942,6.01];
ybs[2793]=['',1.9327614,0.0024941,5.99];
ybs[2794]=['',1.9276773,-0.4524816,5.87];
ybs[2795]=['δ Vol',1.9059748,-1.1866367,3.98];
ybs[2796]=['',1.9474055,0.9049822,5.8];
ybs[2797]=['66 Aur',1.9432087,0.7092497,5.19];
ybs[2798]=['',1.9323958,-0.1573147,6.43];
ybs[2799]=['',1.9337644,-0.0525921,6.23];
ybs[2800]=['57 Gem',1.9396796,0.4366053,5.03];
ybs[2801]=['',1.9597873,1.157067,6.47];
ybs[2802]=['58 Gem',1.9395879,0.3998612,6.02];
ybs[2803]=['',1.9341823,-0.1050206,5.82];
ybs[2804]=['',1.9329352,-0.3325031,4.96];
ybs[2805]=['',1.9232642,-0.9135963,6.05];
ybs[2806]=['',1.9232862,-0.9135625,6.6];
ybs[2807]=['',1.92454,-0.9096617,5.39];
ybs[2808]=['59 Gem',1.944487,0.4817581,5.76];
ybs[2809]=['',1.9436829,0.2702106,6.41];
ybs[2810]=['21 Lyn',1.9548962,0.8582681,4.64];
ybs[2811]=['',1.9359185,-0.5577818,5.43];
ybs[2812]=['1 CMi',1.945788,0.2030556,5.3];
ybs[2813]=['ι Gem',1.9495961,0.4845425,3.79];
ybs[2814]=['',1.9381429,-0.4864062,5.38];
ybs[2815]=['',1.9381848,-0.5626433,5.39];
ybs[2816]=['',1.9398897,-0.5279963,6.6];
ybs[2817]=['',1.9436842,-0.2833794,5.33];
ybs[2818]=['',1.941809,-0.4005173,6.19];
ybs[2819]=['η CMa',1.9407526,-0.5120472,2.45];
ybs[2820]=['ε CMi',1.9486791,0.1612746,4.99];
ybs[2821]=['',1.9399586,-0.6260986,6.31];
ybs[2822]=['',1.9752058,1.1942865,5.64];
ybs[2823]=['',1.9443696,-0.3324438,6.24];
ybs[2824]=['',1.9458079,-0.2406368,5.78];
ybs[2825]=['',1.9491407,-0.1014179,5.97];
ybs[2826]=['',1.9434307,-0.5557866,5.35];
ybs[2827]=['',1.9542286,0.3752386,6.54];
ybs[2828]=['',1.9522799,0.1845206,6.37];
ybs[2829]=['61 Gem',1.9546318,0.352922,5.93];
ybs[2830]=['',1.9500842,-0.0798209,6.76];
ybs[2831]=['',1.9464,-0.3842931,6.05];
ybs[2832]=['',1.9532657,0.1915149,6.41];
ybs[2833]=['',1.9466862,-0.440755,5.78];
ybs[2834]=['',1.9434541,-0.6514503,6.97];
ybs[2835]=['',1.9434612,-0.6514648,6.84];
ybs[2836]=['',1.9641912,0.8403201,5.72];
ybs[2837]=['β CMi',1.9552028,0.1440437,2.9];
ybs[2838]=['63 Gem',1.9581662,0.3736471,5.22];
ybs[2839]=['',1.9477401,-0.5545669,6.31];
ybs[2840]=['',1.7474642,-1.5171479,6.47];
ybs[2841]=['22 Lyn',1.9689655,0.8662939,5.36];
ybs[2842]=['',1.9522212,-0.4144868,6.56];
ybs[2843]=['η CMi',1.9590263,0.1205193,5.25];
ybs[2844]=['ρ Gem',1.9645062,0.5540947,4.18];
ybs[2845]=['',1.9543916,-0.3124271,5.63];
ybs[2846]=['γ CMi',1.9596417,0.1551389,4.32];
ybs[2847]=['',1.9536019,-0.4035613,5.61];
ybs[2848]=['',1.9519602,-0.5965003,5.9];
ybs[2849]=['64 Gem',1.9653681,0.4901027,5.05];
ybs[2850]=['',1.9625462,0.2630638,6.22];
ybs[2851]=['',1.9577559,-0.2023456,5.79];
ybs[2852]=['',1.9567735,-0.399615,5.95];
ybs[2853]=['65 Gem',1.9674183,0.4865749,5.01];
ybs[2854]=['',1.9495939,-0.8910652,5.1];
ybs[2855]=['',1.9577444,-0.5095046,5.54];
ybs[2856]=['6 CMi',1.9668551,0.2089032,4.54];
ybs[2857]=['',1.9643527,-0.0339024,5.59];
ybs[2858]=['',1.9647057,-0.1324415,5.86];
ybs[2859]=['',1.9643658,-0.1808835,5.75];
ybs[2860]=['',1.9642166,-0.2624339,6.05];
ybs[2861]=['',1.9591176,-0.6605555,6.58];
ybs[2862]=['',1.9614349,-0.5565035,6.38];
ybs[2863]=['',1.9614494,-0.5564792,7.13];
ybs[2864]=['',1.9770836,0.6782021,6.54];
ybs[2865]=['',1.9624391,-0.5496644,5.77];
ybs[2866]=['',1.9661155,-0.4025047,4.85];
ybs[2867]=['',1.9621942,-0.6780477,5.43];
ybs[2868]=['',1.9709811,-0.0918772,6.24];
ybs[2869]=['',1.9757723,0.2975424,5.42];
ybs[2870]=['σ Pup',1.9625663,-0.7563992,3.25];
ybs[2871]=['',1.9804738,0.3987934,6.54];
ybs[2872]=['δ1 CMi',1.9766164,0.0327452,5.25];
ybs[2873]=['',1.9695582,-0.5410506,4.65];
ybs[2874]=['',1.9692952,-0.6523588,6.65];
ybs[2875]=['',1.9763124,-0.1556677,5.9];
ybs[2876]=['',1.9653206,-0.9195873,5.87];
ybs[2877]=['',1.9725036,-0.6316477,6.68];
ybs[2878]=['68 Gem',1.9835946,0.2755486,5.25];
ybs[2879]=['δ2 CMi',1.9814317,0.0567505,5.59];
ybs[2880]=['',1.9591086,-1.1265553,6.39];
ybs[2881]=['',1.9737514,-0.6270245,6.61];
ybs[2882]=['α Gem',1.9884475,0.5558754,2.88];
ybs[2883]=['α Gem',1.9884475,0.5558705,1.98];
ybs[2884]=['',1.967463,-0.950105,5.96];
ybs[2885]=['',1.9855232,0.1837703,6.28];
ybs[2886]=['',1.9992804,0.9724115,5.92];
ybs[2887]=['',1.9766223,-0.6283137,6.3];
ybs[2888]=['',1.99079,0.5396837,5.33];
ybs[2889]=['',1.9817067,-0.2509276,6.21];
ybs[2890]=['',1.9947456,0.7503393,6.3];
ybs[2891]=['',1.9813666,-0.3394881,5.66];
ybs[2892]=['',1.9804854,-0.4319602,5.85];
ybs[2893]=['δ3 CMi',1.9861024,0.0581588,5.81];
ybs[2894]=['',1.9835779,-0.2541691,4.97];
ybs[2895]=['',1.9975009,0.8052984,5.65];
ybs[2896]=['',1.9882809,0.0468738,6.55];
ybs[2897]=['υ Gem',1.9940337,0.468726,4.06];
ybs[2898]=['',1.9844554,-0.3898215,4.45];
ybs[2899]=['',1.9801673,-0.6998343,6.26];
ybs[2900]=['',1.9800052,-0.752674,6.52];
ybs[2901]=['',1.9855379,-0.4103745,5.83];
ybs[2902]=['',1.9855742,-0.4103939,5.87];
ybs[2903]=['',1.9830602,-0.6349028,5.54];
ybs[2904]=['',1.9861927,-0.4565056,6.65];
ybs[2905]=['',1.9847626,-0.584727,6.11];
ybs[2906]=['',2.0036346,0.8505512,5.92];
ybs[2907]=['',2.0005455,0.6978685,6.38];
ybs[2908]=['',1.986606,-0.4721318,5.77];
ybs[2909]=['',1.9834488,-0.6971681,6.76];
ybs[2910]=['',1.9962709,0.1016067,5.91];
ybs[2911]=['ε Men',1.9401186,-1.3810709,5.53];
ybs[2912]=['',1.9945688,-0.1457575,6.27];
ybs[2913]=['',1.9934708,-0.2536415,5.7];
ybs[2914]=['',1.9900436,-0.49583,4.64];
ybs[2915]=['',1.9935234,-0.3874698,6.34];
ybs[2916]=['70 Gem',2.005783,0.611001,5.56];
ybs[2917]=['',1.985778,-0.8990874,6.28];
ybs[2918]=['',2.0040689,0.4244567,6.27];
ybs[2919]=['25 Mon',1.9990548,-0.0724556,5.13];
ybs[2920]=['',1.9960148,-0.3445675,5.74];
ybs[2921]=['23 Lyn',2.0170131,0.9955534,6.06];
ybs[2922]=['ο Gem',2.0084611,0.6028907,4.9];
ybs[2923]=['',2.0082373,0.4220456,6.17];
ybs[2924]=['',2.0003821,-0.2527504,6.53];
ybs[2925]=['',1.9984911,-0.4156548,6.37];
ybs[2926]=['',1.9900837,-0.9175801,4.94];
ybs[2927]=['',2.01332,0.6685126,5.73];
ybs[2928]=['',2.0115748,0.5579533,6.17];
ybs[2929]=['',1.9984756,-0.6110205,4.53];
ybs[2930]=['74 Gem',2.0092492,0.3077636,5.05];
ybs[2931]=['',2.0180189,0.839325,5.56];
ybs[2932]=['',1.9950071,-0.8529473,5.72];
ybs[2933]=['',1.9914184,-0.976114,6.39];
ybs[2934]=['',2.0001148,-0.6164093,6.6];
ybs[2935]=['α CMi',2.0081356,0.0904764,0.38];
ybs[2936]=['',2.0028912,-0.4434075,4.7];
ybs[2937]=['',2.0000291,-0.6641149,6.38];
ybs[2938]=['24 Lyn',2.0266752,1.023944,4.99];
ybs[2939]=['',2.006673,-0.3267283,6.72];
ybs[2940]=['',2.0051202,-0.4684904,4.5];
ybs[2941]=['',2.0051565,-0.4685245,4.62];
ybs[2942]=['',2.0116916,0.0905728,6.02];
ybs[2943]=['',2.0159485,0.4010219,5.89];
ybs[2944]=['',2.0027783,-0.6986913,6.59];
ybs[2945]=['',2.014858,0.2396192,6.24];
ybs[2946]=['',2.0043647,-0.637704,5.8];
ybs[2947]=['',2.0034442,-0.677569,6.19];
ybs[2948]=['',2.0078534,-0.4695661,6.5];
ybs[2949]=['',2.0018845,-0.8489583,5.68];
ybs[2950]=['',2.013402,-0.1435952,6.01];
ybs[2951]=['',2.0123088,-0.2671242,4.94];
ybs[2952]=['',2.0114807,-0.3438691,5.93];
ybs[2953]=['',2.0074525,-0.6693236,4.84];
ybs[2954]=['',2.0239638,0.5926761,6.02];
ybs[2955]=['',2.0086531,-0.6663778,5.73];
ybs[2956]=['',2.0089389,-0.6684968,5.76];
ybs[2957]=['',2.0195396,0.2345457,5.77];
ybs[2958]=['',2.018062,0.0625309,5.94];
ybs[2959]=['',2.0203891,0.2472465,5.56];
ybs[2960]=['',2.0097106,-0.6566055,6];
ybs[2961]=['',2.0306463,0.8794871,5.27];
ybs[2962]=['α Mon',2.0162227,-0.1674281,3.93];
ybs[2963]=['',2.0046375,-0.9305085,6.06];
ybs[2964]=['',2.0133818,-0.4884727,6.76];
ybs[2965]=['σ Gem',2.0263346,0.5033696,4.28];
ybs[2966]=['',2.0139258,-0.5533126,6.56];
ybs[2967]=['51 Cam',2.0434564,1.1416511,5.92];
ybs[2968]=['',2.0164901,-0.3905886,6.18];
ybs[2969]=['49 Cam',2.0421688,1.0958331,6.49];
ybs[2970]=['',2.0263759,0.3901994,6.21];
ybs[2971]=['',1.9852701,-1.297039,7.16];
ybs[2972]=['',1.9852774,-1.2970391,7.26];
ybs[2973]=['',2.0153302,-0.6732676,5.42];
ybs[2974]=['',2.0245293,0.0025641,6.19];
ybs[2975]=['76 Gem',2.0297338,0.4492689,5.31];
ybs[2976]=['',2.0154733,-0.7797088,6.41];
ybs[2977]=['κ Gem',2.0311355,0.4250745,3.57];
ybs[2978]=['',2.018401,-0.6731899,6.54];
ybs[2979]=['',2.0298613,0.2236894,6.43];
ybs[2980]=['',2.0225076,-0.4606535,5.64];
ybs[2981]=['',2.0290971,0.0412259,6.47];
ybs[2982]=['β Gem',2.0350428,0.4883901,1.14];
ybs[2983]=['79 Gem',2.0340983,0.3538314,6.33];
ybs[2984]=['',2.0262457,-0.4458723,6.55];
ybs[2985]=['1 Pup',2.0256606,-0.4966121,4.59];
ybs[2986]=['',2.0238903,-0.6299382,5.6];
ybs[2987]=['',2.0234071,-0.6790493,6.89];
ybs[2988]=['3 Pup',2.0268131,-0.5061017,3.96];
ybs[2989]=['',2.0907248,1.4000626,6.56];
ybs[2990]=['',2.0223872,-0.7891637,5.06];
ybs[2991]=['',2.0412445,0.654037,5.18];
ybs[2992]=['',1.986793,-1.3556614,6.18];
ybs[2993]=['',2.0260483,-0.6674953,6.4];
ybs[2994]=['',2.0258527,-0.7151765,5.17];
ybs[2995]=['81 Gem',2.0382683,0.3222975,4.88];
ybs[2996]=['',2.0302666,-0.4313923,5.62];
ybs[2997]=['',2.0228091,-0.8732797,6.57];
ybs[2998]=['',2.0179928,-1.0240354,6.43];
ybs[2999]=['',2.0280809,-0.6301628,5.8];
ybs[3000]=['11 CMi',2.0386886,0.1871792,5.3];
ybs[3001]=['2 Pup',2.0345464,-0.2570789,6.89];
ybs[3002]=['2 Pup',2.0345753,-0.2571614,6.07];
ybs[3003]=['',2.0297922,-0.6629825,5.88];
ybs[3004]=['',2.0211701,-1.0170443,6.21];
ybs[3005]=['π Gem',2.0447875,0.5824391,5.14];
ybs[3006]=['',2.037201,-0.1189642,5.49];
ybs[3007]=['4 Pup',2.0366009,-0.2549489,5.04];
ybs[3008]=['',2.0320067,-0.6620212,6.54];
ybs[3009]=['',2.0327818,-0.6634332,3.61];
ybs[3010]=['',2.0343698,-0.5971905,5.37];
ybs[3011]=['',2.0401453,-0.2219917,6.39];
ybs[3012]=['',2.0327157,-0.7643763,6.03];
ybs[3013]=['82 Gem',2.0490375,0.4031095,6.18];
ybs[3014]=['',2.0367926,-0.6628333,5.88];
ybs[3015]=['',2.041865,-0.3938125,5.9];
ybs[3016]=['ζ Vol',2.0140748,-1.2679458,3.95];
ybs[3017]=['',2.0383728,-0.6999385,6.57];
ybs[3018]=['',2.0439504,-0.279865,6.34];
ybs[3019]=['',2.0444371,-0.2802778,6.43];
ybs[3020]=['',2.0615288,0.9439353,6.02];
ybs[3021]=['5 Pup',2.045381,-0.2135833,5.48];
ybs[3022]=['',2.0508158,0.2325831,6.04];
ybs[3023]=['',2.033126,-0.9907516,6.12];
ybs[3024]=['',2.0407631,-0.6872302,6.31];
ybs[3025]=['',2.0503459,0.0748397,6.53];
ybs[3026]=['ο Pup',2.04559,-0.4534649,4.5];
ybs[3027]=['',2.0421962,-0.6729159,5.08];
ybs[3028]=['',2.0282458,-1.1539237,6.38];
ybs[3029]=['',2.0422715,-0.8142443,5.23];
ybs[3030]=['',2.025279,-1.2193598,6.18];
ybs[3031]=['',2.0682162,0.9627798,6.38];
ybs[3032]=['',2.0601898,0.5792402,6.03];
ybs[3033]=['',2.0452662,-0.71029,6.14];
ybs[3034]=['',2.0520332,-0.233844,6.23];
ybs[3035]=['',2.0497335,-0.4355815,5.33];
ybs[3036]=['6 Pup',2.0528399,-0.3014768,5.18];
ybs[3037]=['ξ Pup',2.0508991,-0.434667,3.34];
ybs[3038]=['',2.0457967,-0.8224381,4.71];
ybs[3039]=['',2.0552021,-0.1610685,5.61];
ybs[3040]=['',2.0530435,-0.3534638,6.56];
ybs[3041]=['',2.0503135,-0.6158946,5.93];
ybs[3042]=['',2.0582079,0.0564049,6.18];
ybs[3043]=['',2.054555,-0.3415396,6.12];
ybs[3044]=['',2.051895,-0.5817854,5.6];
ybs[3045]=['',2.0636873,0.3364884,5.99];
ybs[3046]=['',2.0583932,-0.1950248,6.16];
ybs[3047]=['',2.0497666,-0.8101493,4.11];
ybs[3048]=['',2.0450591,-0.9863821,6.33];
ybs[3049]=['',2.0508655,-0.7818523,6.32];
ybs[3050]=['',2.0496376,-0.8186042,5.84];
ybs[3051]=['ζ CMi',2.0621387,0.0300397,5.14];
ybs[3052]=['',2.0583519,-0.4288944,6.45];
ybs[3053]=['',2.0640109,0.0563962,6.31];
ybs[3054]=['',2.0485227,-0.9853305,5.59];
ybs[3055]=['8 Pup',2.0616708,-0.2245404,6.36];
ybs[3056]=['9 Pup',2.0620337,-0.2433663,5.17];
ybs[3057]=['25 Lyn',2.0758781,0.826225,6.25];
ybs[3058]=['26 Lyn',2.0768619,0.8293409,5.45];
ybs[3059]=['φ Gem',2.0706854,0.4663403,4.97];
ybs[3060]=['',2.0615764,-0.370353,5.63];
ybs[3061]=['',2.0562746,-0.7788546,6.45];
ybs[3062]=['',2.0486014,-1.0529285,5.78];
ybs[3063]=['',2.0545842,-0.8823504,5.91];
ybs[3064]=['',2.0667425,-0.0955438,5.76];
ybs[3065]=['10 Pup',2.0643785,-0.2599214,5.69];
ybs[3066]=['',2.0591224,-0.7529552,6.32];
ybs[3067]=['',2.1042435,1.2892553,5.41];
ybs[3068]=['',2.0516773,-1.0488752,6.72];
ybs[3069]=['',2.0850434,0.9853567,6.72];
ybs[3070]=['',2.0605735,-0.7493406,6.04];
ybs[3071]=['',2.0635045,-0.6065236,5.01];
ybs[3072]=['',2.0630622,-0.7089836,3.73];
ybs[3073]=['',2.0498689,-1.1561187,5.79];
ybs[3074]=['',2.1269726,1.3862935,5.42];
ybs[3075]=['',2.0805149,0.6172436,6.23];
ybs[3076]=['',2.0650082,-0.6790929,4.49];
ybs[3077]=['',2.0669171,-0.6354769,5.43];
ybs[3078]=['85 Gem',2.0799306,0.3462138,5.35];
ybs[3079]=['',2.0790068,0.1538603,5.86];
ybs[3080]=['',2.063469,-0.9496899,5.7];
ybs[3081]=['',2.0662734,-0.8667181,4.63];
ybs[3082]=['',2.0674201,-0.8403652,4.24];
ybs[3083]=['',2.0718481,-0.6269952,5.49];
ybs[3084]=['',2.0739914,-0.6090118,6.15];
ybs[3085]=['',2.0827092,0.0774626,6.17];
ybs[3086]=['',2.0922129,0.7667091,6.34];
ybs[3087]=['1 Cnc',2.0855934,0.2747583,5.78];
ybs[3088]=['',2.0765959,-0.5404338,6.44];
ybs[3089]=['',2.0865993,0.149985,6.05];
ybs[3090]=['',2.0864252,0.0188334,6.35];
ybs[3091]=['',2.0816456,-0.5294069,6.33];
ybs[3092]=['',2.074579,-0.9185669,6.38];
ybs[3093]=['',2.0784442,-0.7660644,6.02];
ybs[3094]=['11 Pup',2.083974,-0.4001637,4.2];
ybs[3095]=['',2.0901901,0.1250603,6.41];
ybs[3096]=['',2.0923098,0.2874603,5.99];
ybs[3097]=['',2.0736733,-1.0009453,5.63];
ybs[3098]=['',2.1066175,1.0297098,5.77];
ybs[3099]=['',2.0813428,-0.7118131,6.78];
ybs[3100]=['',2.1853998,1.4661165,6.49];
ybs[3101]=['53 Cam',2.108321,1.0519943,6.01];
ybs[3102]=['14 CMi',2.0911368,0.0379862,5.29];
ybs[3103]=['',2.0837106,-0.7409588,6.09];
ybs[3104]=['',2.1121491,1.1002618,6.4];
ybs[3105]=['',2.087268,-0.5302783,4.79];
ybs[3106]=['',2.0836577,-0.7600556,5.35];
ybs[3107]=['',2.0968436,0.2302697,6.02];
ybs[3108]=['',2.0851194,-0.7706946,5.09];
ybs[3109]=['χ Car',2.0822821,-0.9255451,3.47];
ybs[3110]=['',2.0850327,-0.8366777,6.22];
ybs[3111]=['',2.1119087,0.9987417,6.49];
ybs[3112]=['',2.0795768,-1.0572122,5.74];
ybs[3113]=['',2.0874747,-0.7963205,5.17];
ybs[3114]=['27 Mon',2.097057,-0.0650746,4.93];
ybs[3115]=['12 Pup',2.0937237,-0.4076927,5.11];
ybs[3116]=['ω1 Cnc',2.1030695,0.442328,5.83];
ybs[3117]=['',2.1023275,0.3449978,6.25];
ybs[3118]=['',2.0820585,-1.0327811,6.25];
ybs[3119]=['',2.1033624,0.4107419,6.34];
ybs[3120]=['3 Cnc',2.1022051,0.3012338,5.55];
ybs[3121]=['',2.0889125,-0.8603279,4.41];
ybs[3122]=['',2.1077184,0.6172083,6.34];
ybs[3123]=['',2.0972366,-0.3219779,4.61];
ybs[3124]=['ω2 Cnc',2.106542,0.4370336,6.31];
ybs[3125]=['',2.0892822,-0.8987889,6.44];
ybs[3126]=['5 Cnc',2.1053075,0.2863357,5.99];
ybs[3127]=['',2.1014411,-0.0511524,6.51];
ybs[3128]=['',2.1038001,0.0843061,5.65];
ybs[3129]=['',2.0925886,-0.7900158,5.99];
ybs[3130]=['',2.0860313,-1.0533286,5.6];
ybs[3131]=['',2.0831775,-1.105573,6.14];
ybs[3132]=['',2.0948037,-0.6867147,5.24];
ybs[3133]=['28 Mon',2.1036028,-0.0251646,4.68];
ybs[3134]=['',2.0930823,-0.8731041,6.32];
ybs[3135]=['',2.0931626,-0.8730509,6.34];
ybs[3136]=['',2.1065879,0.1547116,6.22];
ybs[3137]=['',2.1082505,0.0398761,4.39];
ybs[3138]=['',2.0982387,-0.7942274,6.61];
ybs[3139]=['',2.0906075,-1.0624305,5.81];
ybs[3140]=['',2.0976999,-0.85574,6.02];
ybs[3141]=['χ Gem',2.1144276,0.4842236,4.94];
ybs[3142]=['',2.108763,-0.111474,6.33];
ybs[3143]=['',2.0987256,-0.8538216,6.12];
ybs[3144]=['',2.0942929,-1.0516682,6.33];
ybs[3145]=['',2.0940634,-1.0582905,5.17];
ybs[3146]=['',2.1042751,-0.6515845,5.95];
ybs[3147]=['',2.1063802,-0.64752,6.34];
ybs[3148]=['',2.0999028,-0.9459769,5.87];
ybs[3149]=['',2.1022781,-0.9523266,6.1];
ybs[3150]=['',2.1195439,0.327975,6.15];
ybs[3151]=['',2.0968008,-1.1103153,4.82];
ybs[3152]=['',2.1107671,-0.5674737,5.82];
ybs[3153]=['',2.102831,-0.9687336,6.28];
ybs[3154]=['',2.1090079,-0.7218652,5.52];
ybs[3155]=['8 Cnc',2.1207901,0.2280676,5.12];
ybs[3156]=['',2.1235705,0.4795949,6.21];
ybs[3157]=['ζ Pup',2.1127265,-0.699065,2.25];
ybs[3158]=['',2.1121852,-0.7504691,6.29];
ybs[3159]=['28 Lyn',2.1309285,0.7541347,6.26];
ybs[3160]=['14 Pup',2.118258,-0.3452022,6.13];
ybs[3161]=['μ1 Cnc',2.1264291,0.3941713,5.99];
ybs[3162]=['',2.1160058,-0.571166,5.31];
ybs[3163]=['',2.0902068,-1.2792064,6.34];
ybs[3164]=['',2.1237171,-0.0109016,6.41];
ybs[3165]=['27 Lyn',2.137014,0.8980528,4.84];
ybs[3166]=['',2.1262523,-0.1622497,6.23];
ybs[3167]=['',2.1445282,1.0157019,5.93];
ybs[3168]=['μ2 Cnc',2.1327517,0.3757683,5.3];
ybs[3169]=['',2.1224275,-0.5867815,6.14];
ybs[3170]=['',2.1170743,-0.8838536,5.95];
ybs[3171]=['',2.1200521,-0.8208223,6.19];
ybs[3172]=['',2.1184791,-0.9277944,5.53];
ybs[3173]=['',2.1405682,0.7396394,6.27];
ybs[3174]=['',2.1578107,1.1941632,5.32];
ybs[3175]=['',2.129619,-0.3596467,5.38];
ybs[3176]=['12 Cnc',2.1366467,0.2371689,6.27];
ybs[3177]=['ρ Pup',2.1305704,-0.4250883,2.81];
ybs[3178]=['',2.1160444,-1.0975781,6.3];
ybs[3179]=['',2.1259655,-0.790942,5.05];
ybs[3180]=['ζ Mon',2.1357285,-0.0529863,4.34];
ybs[3181]=['',2.1370643,-0.1988252,6.32];
ybs[3182]=['',2.1358443,-0.3563105,6.36];
ybs[3183]=['ψ Cnc',2.1446015,0.4442652,5.73];
ybs[3184]=['16 Pup',2.1371927,-0.3367986,4.4];
ybs[3185]=['',2.1490189,0.6750646,6.58];
ybs[3186]=['',2.1392323,-0.2845095,5.68];
ybs[3187]=['',2.1348376,-0.6585714,6.37];
ybs[3188]=['',2.1372248,-0.5301378,6.65];
ybs[3189]=['',2.2156623,1.4376769,6.32];
ybs[3190]=['',2.14659,0.2544092,6.23];
ybs[3191]=['',2.1372966,-0.619717,6.2];
ybs[3192]=['',2.1608155,0.9843353,5.85];
ybs[3193]=['',2.1477587,0.1704861,6.07];
ybs[3194]=['18 Pup',2.1444843,-0.2417612,5.54];
ybs[3195]=['',2.1366409,-0.8506139,5.7];
ybs[3196]=['',2.1387966,-0.7710006,5.21];
ybs[3197]=['',2.1397259,-0.7451323,6.26];
ybs[3198]=['γ1 Vel',2.1381502,-0.8272528,4.27];
ybs[3199]=['γ2 Vel',2.1383472,-0.8270931,1.78];
ybs[3200]=['ζ1 Cnc',2.1520442,0.3070812,5.63];
ybs[3201]=['ζ1 Cnc',2.1520442,0.3070812,6.02];
ybs[3202]=['ζ2 Cnc',2.1520877,0.3070812,6.2];
ybs[3203]=['19 Pup',2.1471622,-0.2265418,4.72];
ybs[3204]=['',2.1485097,-0.1365817,5.36];
ybs[3205]=['',2.1391293,-0.8375808,5.23];
ybs[3206]=['',2.1526297,0.2434825,6.54];
ybs[3207]=['15 Cnc',2.1564862,0.5166699,5.64];
ybs[3208]=['',2.1890193,1.3212291,5.54];
ybs[3209]=['',2.1320669,-1.1144439,6.28];
ybs[3210]=['',2.1378854,-0.9797898,5.66];
ybs[3211]=['',2.1453324,-0.6517942,6.44];
ybs[3212]=['',2.1350007,-1.0708389,4.76];
ybs[3213]=['',2.1699036,1.0528851,6.45];
ybs[3214]=['',2.1554376,0.2872912,6.01];
ybs[3215]=['ε Vol',2.1292215,-1.1984973,4.35];
ybs[3216]=['',2.1586691,0.4028905,6.56];
ybs[3217]=['',2.1466859,-0.6923993,4.45];
ybs[3218]=['',2.1468572,-0.751193,4.75];
ybs[3219]=['',2.1454986,-0.8467432,5.82];
ybs[3220]=['',2.1606526,0.3075587,6.47];
ybs[3221]=['20 Pup',2.1560843,-0.2764953,4.99];
ybs[3222]=['',2.153195,-0.5229756,6.52];
ybs[3223]=['',2.1612494,0.2267927,6.38];
ybs[3224]=['',2.1491579,-0.8150221,5.76];
ybs[3225]=['',2.1533052,-0.6628399,6.43];
ybs[3226]=['',2.1514271,-0.8083931,6.03];
ybs[3227]=['29 Lyn',2.1785295,1.0387454,5.64];
ybs[3228]=['',2.1928695,1.2627592,5.98];
ybs[3229]=['',2.1561564,-0.6275058,4.78];
ybs[3230]=['',2.1570788,-0.5868312,6.37];
ybs[3231]=['',2.1593046,-0.5619052,6.06];
ybs[3232]=['',2.1582362,-0.6348876,5.08];
ybs[3233]=['',2.1582646,-0.6352124,6.11];
ybs[3234]=['',2.1593446,-0.620369,5.78];
ybs[3235]=['',2.1584015,-0.7051472,4.44];
ybs[3236]=['',2.1561422,-0.821102,5.13];
ybs[3237]=['',2.1851429,1.0899812,5.71];
ybs[3238]=['',2.179929,0.9440155,6.27];
ybs[3239]=['',2.1558074,-0.8770248,5.51];
ybs[3240]=['',2.1708781,0.2037068,7.13];
ybs[3241]=['β Cnc',2.1705955,0.1593612,3.52];
ybs[3242]=['',2.1596748,-0.8009047,5.83];
ybs[3243]=['',2.1667315,-0.5407099,6.21];
ybs[3244]=['',2.1750083,0.1537796,6.29];
ybs[3245]=['',2.1670211,-0.6275745,6.16];
ybs[3246]=['30 Lyn',2.1896701,1.0068291,5.89];
ybs[3247]=['',2.1715144,-0.3730678,6.6];
ybs[3248]=['',2.1637424,-0.8814576,6.44];
ybs[3249]=['21 Pup',2.1737546,-0.2851887,6.16];
ybs[3250]=['',2.1895703,0.9340682,6.49];
ybs[3251]=['',2.1782725,-0.2214371,5.98];
ybs[3252]=['',2.1621499,-1.0990356,5.16];
ybs[3253]=['',2.1759155,-0.5246221,6.45];
ybs[3254]=['χ Cnc',2.186573,0.4740609,5.14];
ybs[3255]=['',2.1999547,1.0572162,6.41];
ybs[3256]=['',2.1876277,0.3611365,5.83];
ybs[3257]=['',2.1820591,-0.1784006,6.32];
ybs[3258]=['',2.1771169,-0.6197153,5.58];
ybs[3259]=['',2.1766949,-0.6532687,6.7];
ybs[3260]=['λ Cnc',2.1885274,0.4182852,5.98];
ybs[3261]=['',2.184936,0.0679246,6.05];
ybs[3262]=['',2.1782288,-0.6407965,4.45];
ybs[3263]=['',2.186505,-0.0168522,6.18];
ybs[3264]=['',2.1866883,-0.0939912,6.13];
ybs[3265]=['',2.1823877,-0.6046885,6.43];
ybs[3266]=['',2.1741935,-1.0336218,6.42];
ybs[3267]=['31 Lyn',2.1992375,0.7527777,4.25];
ybs[3268]=['',2.1869792,-0.4010923,6.13];
ybs[3269]=['',2.2040277,0.9278571,5.51];
ybs[3270]=['',2.191368,-0.0289502,6.5];
ybs[3271]=['',2.1909729,-0.3514335,5.58];
ybs[3272]=['',2.1751315,-1.1461343,5.07];
ybs[3273]=['',2.1934699,-0.3079297,5.75];
ybs[3274]=['',2.1907084,-0.5778947,4.83];
ybs[3275]=['',2.1904397,-0.6377592,5.2];
ybs[3276]=['20 Cnc',2.2007057,0.3189589,5.95];
ybs[3277]=['',2.1963478,-0.1088399,6.15];
ybs[3278]=['',2.1905501,-0.6924998,6.16];
ybs[3279]=['',2.2073734,0.7321183,6.02];
ybs[3280]=['',2.1980455,-0.1326515,5.96];
ybs[3281]=['22 Pup',2.1973846,-0.2288425,6.11];
ybs[3282]=['21 Cnc',2.2029365,0.1845605,6.08];
ybs[3283]=['',2.1972466,-0.460855,5.9];
ybs[3284]=['',2.2087056,0.6100547,6.06];
ybs[3285]=['',2.1885914,-1.0128045,5.97];
ybs[3286]=['',2.1950378,-0.8473073,4.82];
ybs[3287]=['',2.2055616,-0.0833319,6.01];
ybs[3288]=['',2.1988276,-0.669211,6.32];
ybs[3289]=['1 Hya',2.205483,-0.0664748,5.61];
ybs[3290]=['',2.1876363,-1.1198456,6.12];
ybs[3291]=['25 Cnc',2.2114384,0.2964974,6.14];
ybs[3292]=['',2.1965409,-0.910728,5.85];
ybs[3293]=['κ1 Vol',2.180599,-1.2491462,5.37];
ybs[3294]=['κ2 Vol',2.1814547,-1.2489777,5.65];
ybs[3295]=['',2.2314183,1.1735244,5.88];
ybs[3296]=['φ1 Cnc',2.2144897,0.485818,5.57];
ybs[3297]=['',2.2100228,0.035679,5.73];
ybs[3298]=['',2.2115574,0.1310108,5.13];
ybs[3299]=['ε Car',2.1941929,-1.0396322,1.86];
ybs[3300]=['',2.2064486,-0.405114,5.68];
ybs[3301]=['',2.2201814,0.7957714,6.32];
ybs[3302]=['φ2 Cnc',2.2158491,0.4690755,6.32];
ybs[3303]=['φ2 Cnc',2.2158636,0.46909,6.3];
ybs[3304]=['24 Cnc',2.21527,0.4271835,7.02];
ybs[3305]=['24 Cnc',2.2152919,0.4272029,7.81];
ybs[3306]=['',2.2101775,-0.0691913,3.9];
ybs[3307]=['',2.207049,-0.4206919,5.28];
ybs[3308]=['',2.2082535,-0.3683288,6.01];
ybs[3309]=['',2.2098268,-0.3053875,6.44];
ybs[3310]=['α Cha',2.1731068,-1.3434673,4.07];
ybs[3311]=['27 Cnc',2.2152525,0.219843,5.5];
ybs[3312]=['',2.2110702,-0.2615862,5.98];
ybs[3313]=['2 Hya',2.2136373,-0.0706118,5.59];
ybs[3314]=['',2.2059194,-0.7474746,5.98];
ybs[3315]=['ο UMa',2.2327595,1.058689,3.36];
ybs[3316]=['',2.2144986,-0.2197854,5.54];
ybs[3317]=['',2.2172157,-0.1128925,6.59];
ybs[3318]=['',2.209928,-0.736727,5.47];
ybs[3319]=['',2.2119267,-0.682731,6.53];
ybs[3320]=['',2.2119703,-0.6827457,7.25];
ybs[3321]=['28 Cnc',2.2237635,0.4203749,6.1];
ybs[3322]=['',2.2079104,-0.903835,5.17];
ybs[3323]=['',2.214687,-0.5109215,6.73];
ybs[3324]=['',2.2453414,1.2088051,6.31];
ybs[3325]=['29 Cnc',2.2235361,0.2469959,5.95];
ybs[3326]=['η Vol',2.1899694,-1.2820591,5.29];
ybs[3327]=['',2.2180221,-0.3648175,6.56];
ybs[3328]=['',2.2164712,-0.5538202,6.33];
ybs[3329]=['',2.2225449,-0.0449626,6.39];
ybs[3330]=['',2.2217035,-0.154898,6.43];
ybs[3331]=['',2.2193425,-0.4571231,6.62];
ybs[3332]=['θ Cha',2.1820925,-1.3533363,4.35];
ybs[3333]=['',2.2117882,-0.9226801,6.05];
ybs[3334]=['',2.2239502,-0.1711714,6];
ybs[3335]=['',2.2194758,-0.6138783,5.75];
ybs[3336]=['',2.222515,-0.4037057,6.51];
ybs[3337]=['',2.2238392,-0.3666822,6.67];
ybs[3338]=['',2.2082613,-1.1285085,5.97];
ybs[3339]=['β Vol',2.2075101,-1.1553179,3.77];
ybs[3340]=['',2.2360376,0.6493652,6.18];
ybs[3341]=['',2.2161665,-0.9611562,6.53];
ybs[3342]=['',2.2169619,-0.9275936,5.09];
ybs[3343]=['',2.2421474,0.9259728,6.24];
ybs[3344]=['',2.2637269,1.3030925,6.31];
ybs[3345]=['',2.2261536,-0.4780763,6.7];
ybs[3346]=['2 UMa',2.2521491,1.1359282,5.47];
ybs[3347]=['υ1 Cnc',2.2363876,0.4192477,5.75];
ybs[3348]=['',2.2240837,-0.7717789,5.79];
ybs[3349]=['θ Cnc',2.2366047,0.3147602,5.35];
ybs[3350]=['',2.2236875,-0.8375532,5.33];
ybs[3351]=['',2.2255135,-0.7816322,4.99];
ybs[3352]=['',2.2429758,0.6624556,5.9];
ybs[3353]=['',2.237773,0.1702453,6.83];
ybs[3354]=['',2.2304413,-0.5623282,5.65];
ybs[3355]=['',2.2267466,-0.8096804,5.99];
ybs[3356]=['',2.2303563,-0.6419443,6.69];
ybs[3357]=['32 Lyn',2.2448652,0.6348769,6.24];
ybs[3358]=['η Cnc',2.2415222,0.3557107,5.33];
ybs[3359]=['',2.2353412,-0.3427382,5.42];
ybs[3360]=['',2.2255412,-0.9643008,6.36];
ybs[3361]=['υ2 Cnc',2.2429005,0.4193019,6.36];
ybs[3362]=['',2.2135825,-1.2243787,5.53];
ybs[3363]=['',2.2307251,-0.7818525,6.3];
ybs[3364]=['34 Cnc',2.2410728,0.1746332,6.46];
ybs[3365]=['',2.2342709,-0.6828436,6.31];
ybs[3366]=['',2.2399973,-0.2633658,6.38];
ybs[3367]=['',2.2328715,-0.8364746,6.39];
ybs[3368]=['',2.2458905,0.2303223,6.28];
ybs[3369]=['33 Lyn',2.2508231,0.6345732,5.78];
ybs[3370]=['',2.2455717,0.0819599,5.87];
ybs[3371]=['',2.2759908,1.2839843,6.15];
ybs[3372]=['',2.2478254,0.1464517,6.03];
ybs[3373]=['',2.2420359,-0.4305179,6.19];
ybs[3374]=['',2.2338492,-0.9504026,6.34];
ybs[3375]=['',2.2467273,-0.038615,5.81];
ybs[3376]=['',2.2408614,-0.5508473,6.38];
ybs[3377]=['',2.2412741,-0.60553,6.36];
ybs[3378]=['',2.2364921,-0.9297773,5.69];
ybs[3379]=['35 Cnc',2.2529055,0.3408409,6.58];
ybs[3380]=['',2.2426893,-0.6707588,6.49];
ybs[3381]=['',2.2440104,-0.6790993,5.96];
ybs[3382]=['',2.2430808,-0.8208578,6.24];
ybs[3383]=['π1 UMa',2.2721216,1.1337348,5.64];
ybs[3384]=['',2.2528984,0.0468156,6.33];
ybs[3385]=['',2.1958326,-1.4132168,5.69];
ybs[3386]=['',2.2563085,0.2661992,6.32];
ybs[3387]=['',2.2548862,0.1144689,5.99];
ybs[3388]=['',2.2549081,0.1145124,7.25];
ybs[3389]=['',2.2481175,-0.5700167,6.43];
ybs[3390]=['3 Hya',2.2528896,-0.1403856,5.72];
ybs[3391]=['',2.2477688,-0.6575059,6.3];
ybs[3392]=['',2.2674876,0.9309427,5.66];
ybs[3393]=['',2.2714689,1.0450484,6.48];
ybs[3394]=['',2.2524493,-0.4695787,5.96];
ybs[3395]=['π2 UMa',2.276469,1.121633,4.6];
ybs[3396]=['',2.2508305,-0.6986755,6.47];
ybs[3397]=['',2.2702024,0.9226245,6.42];
ybs[3398]=['36 Cnc',2.2603962,0.167442,5.88];
ybs[3399]=['',2.2482533,-0.8727545,5.01];
ybs[3400]=['',2.2714672,0.9188994,5.91];
ybs[3401]=['',2.2663386,0.5714154,5.94];
ybs[3402]=['δ Hya',2.2627471,0.0984645,4.16];
ybs[3403]=['',2.2616083,-0.0871888,6.19];
ybs[3404]=['37 Cnc',2.2647132,0.1660256,6.53];
ybs[3405]=['',2.2531825,-0.8906652,5.8];
ybs[3406]=['',2.2503287,-1.0135184,4.86];
ybs[3407]=['',2.2500048,-1.017285,5.26];
ybs[3408]=['',2.2654375,-0.1173686,6.51];
ybs[3409]=['',2.2365226,-1.281366,6.12];
ybs[3410]=['σ Hya',2.2674921,0.0572297,4.44];
ybs[3411]=['',2.2610333,-0.5900566,6.48];
ybs[3412]=['η Pyx',2.2629001,-0.4593192,5.27];
ybs[3413]=['',2.2600894,-0.7017857,6.55];
ybs[3414]=['34 Lyn',2.2786164,0.7988501,5.37];
ybs[3415]=['',2.2749809,0.5563943,6.1];
ybs[3416]=['',2.27045,0.1388346,6.45];
ybs[3417]=['',2.2665759,-0.3455624,6.33];
ybs[3418]=['',2.2613413,-0.7513838,4.14];
ybs[3419]=['39 Cnc',2.2737794,0.3481051,6.39];
ybs[3420]=['',2.2749121,0.3422083,6.44];
ybs[3421]=['ε Cnc',2.275265,0.3400262,6.3];
ybs[3422]=['',2.2685143,-0.3966158,5.05];
ybs[3423]=['6 Hya',2.2726621,-0.2188302,4.98];
ybs[3424]=['',2.2585514,-1.0980811,5.47];
ybs[3425]=['ζ Pyx',2.2708293,-0.5170321,4.89];
ybs[3426]=['',2.2691156,-0.640003,6.13];
ybs[3427]=['',2.2656375,-0.9276923,6.47];
ybs[3428]=['',2.287318,0.8174654,6.22];
ybs[3429]=['',2.2771141,-0.1590874,6.63];
ybs[3430]=['β Pyx',2.2723703,-0.6173423,3.97];
ybs[3431]=['',2.2731432,-0.703839,5.2];
ybs[3432]=['',2.2684482,-0.9337901,5.48];
ybs[3433]=['9 Hya',2.2799837,-0.2793688,4.88];
ybs[3434]=['',2.2709343,-0.9270786,5.19];
ybs[3435]=['',2.2671449,-1.0538237,6.36];
ybs[3436]=['',2.2741001,-0.7898366,5.71];
ybs[3437]=['',2.2741992,-0.8152749,3.84];
ybs[3438]=['',2.2820066,-0.2099556,6.45];
ybs[3439]=['',2.2723992,-0.9247583,3.62];
ybs[3440]=['',2.2723788,-0.9263873,5.61];
ybs[3441]=['γ Cnc',2.2876793,0.3735835,4.66];
ybs[3442]=['45 Cnc',2.2871144,0.2202085,5.64];
ybs[3443]=['',2.2919897,0.6432219,6.33];
ybs[3444]=['',2.2767513,-0.8269381,4.77];
ybs[3445]=['',2.2761018,-0.8549595,5.9];
ybs[3446]=['η Hya',2.2869888,0.0582032,4.3];
ybs[3447]=['',2.2739685,-1.005453,6.34];
ybs[3448]=['',2.2800308,-0.7936743,5.23];
ybs[3449]=['',2.2733014,-1.0441258,4.33];
ybs[3450]=['',2.2903694,0.0745371,6.37];
ybs[3451]=['',2.288706,-0.1273665,4.62];
ybs[3452]=['θ Vol',2.2651827,-1.2295719,5.2];
ybs[3453]=['δ Cnc',2.2937018,0.3157278,3.94];
ybs[3454]=['',2.281292,-0.8405962,5.51];
ybs[3455]=['',2.2847793,-0.6284411,6.42];
ybs[3456]=['46 Cnc',2.2969699,0.5346512,6.13];
ybs[3457]=['49 Cnc',2.2937937,0.1748359,5.66];
ybs[3458]=['',2.2812221,-0.9278773,5.52];
ybs[3459]=['',2.281702,-0.9281204,4.86];
ybs[3460]=['α Pyx',2.2876725,-0.5803271,3.68];
ybs[3461]=['10 Hya',2.2948789,0.0980205,6.13];
ybs[3462]=['',2.3142528,1.1631286,6.2];
ybs[3463]=['',2.2811892,-0.9745504,6.29];
ybs[3464]=['',2.2961122,-0.0465187,6.41];
ybs[3465]=['',2.2938138,-0.3705703,6.11];
ybs[3466]=['ι Cnc',2.3026253,0.5009154,6.57];
ybs[3467]=['ι Cnc',2.3027559,0.5008232,4.02];
ybs[3468]=['',2.2873403,-0.8706868,5.16];
ybs[3469]=['',2.2908599,-0.745488,4.07];
ybs[3470]=['',2.2991571,-0.0368893,5.7];
ybs[3471]=['',2.2931002,-0.6494636,5.76];
ybs[3472]=['',2.2992753,-0.1932274,6.25];
ybs[3473]=['50 Cnc',2.3033649,0.2102248,5.87];
ybs[3474]=['ε Hya',2.3025522,0.1108971,3.38];
ybs[3475]=['',2.2976385,-0.4442234,6.1];
ybs[3476]=['12 Hya',2.300351,-0.2375844,4.32];
ybs[3477]=['δ Vel',2.291555,-0.9559614,1.96];
ybs[3478]=['',2.304433,-0.0342489,5.29];
ybs[3479]=['',2.2978246,-0.8047073,3.91];
ybs[3480]=['',2.2996421,-0.7189071,6.21];
ybs[3481]=['',2.2929386,-1.0260674,6.21];
ybs[3482]=['',2.301723,-0.6054147,6.37];
ybs[3483]=['',2.2866644,-1.1916337,6.32];
ybs[3484]=['ρ Hya',2.3097669,0.1007456,4.36];
ybs[3485]=['',2.3079588,-0.11561,6.09];
ybs[3486]=['',2.2999346,-0.8024603,5.46];
ybs[3487]=['',2.28962,-1.1499918,6.05];
ybs[3488]=['',2.30344,-0.8067021,5.75];
ybs[3489]=['',2.3052013,-0.7295897,6.36];
ybs[3490]=['',2.3001729,-0.9919505,4.49];
ybs[3491]=['',2.3196349,0.579783,6.25];
ybs[3492]=['14 Hya',2.3136133,-0.0612404,5.31];
ybs[3493]=['',2.3072324,-0.7422751,6.43];
ybs[3494]=['η Cha',2.2719496,-1.379269,5.47];
ybs[3495]=['',2.3060946,-0.9235555,6.3];
ybs[3496]=['',2.3201791,0.3275287,6.16];
ybs[3497]=['5 UMa',2.3335202,1.0802738,5.73];
ybs[3498]=['',2.3320503,1.0295543,6.25];
ybs[3499]=['',2.3148937,-0.3685172,6.47];
ybs[3500]=['35 Lyn',2.3261543,0.7620117,5.15];
ybs[3501]=['',2.3272968,0.7896883,5.99];
ybs[3502]=['54 Cnc',2.3212885,0.2667606,6.38];
ybs[3503]=['',2.3270438,0.7319182,5.99];
ybs[3504]=['',2.3150489,-0.5732786,5.21];
ybs[3505]=['',2.3159267,-0.5153783,5.87];
ybs[3506]=['',2.3138968,-0.7048752,5.48];
ybs[3507]=['',2.3173626,-0.5006321,6.17];
ybs[3508]=['',2.3149014,-0.6843009,6.39];
ybs[3509]=['γ Pyx',2.3181376,-0.4847844,4.01];
ybs[3510]=['σ1 Cnc',2.3285045,0.5656154,5.66];
ybs[3511]=['',2.3143121,-0.791924,4.93];
ybs[3512]=['53 Cnc',2.3279461,0.4920504,6.23];
ybs[3513]=['ρ1 Cnc',2.3284714,0.4933005,5.95];
ybs[3514]=['15 Hya',2.3231897,-0.1264259,5.54];
ybs[3515]=['',2.28014,-1.3811364,6.05];
ybs[3516]=['',2.3168712,-0.7357616,6];
ybs[3517]=['',2.3270764,0.0920362,6.33];
ybs[3518]=['',2.3175915,-0.8132407,5.1];
ybs[3519]=['',2.3344935,0.6190878,6.14];
ybs[3520]=['',2.3271419,-0.2321348,6.13];
ybs[3521]=['',2.3217249,-0.743001,6.55];
ybs[3522]=['6 UMa',2.3479395,1.1263621,5.58];
ybs[3523]=['57 Cnc',2.3357227,0.5325373,5.39];
ybs[3524]=['',2.3263085,-0.568556,6.5];
ybs[3525]=['',2.3270893,-0.6390053,6.42];
ybs[3526]=['',2.3276983,-0.67703,5.82];
ybs[3527]=['',2.3215337,-1.0070551,5.59];
ybs[3528]=['',2.3160849,-1.1669113,5.35];
ybs[3529]=['',2.3351121,-0.0960233,6];
ybs[3530]=['',2.3266135,-0.8451914,5.91];
ybs[3531]=['ρ2 Cnc',2.3418254,0.4862446,5.22];
ybs[3532]=['',2.3403355,0.2995639,6.64];
ybs[3533]=['',2.3265736,-0.9109904,6.39];
ybs[3534]=['',2.291765,-1.3887387,5.79];
ybs[3535]=['',2.3117701,-1.2673988,6.11];
ybs[3536]=['',2.3475323,0.7952389,5.74];
ybs[3537]=['',2.3459005,0.7004645,5.89];
ybs[3538]=['ζ Hya',2.3401343,0.102589,3.11];
ybs[3539]=['',2.3321997,-0.7071134,6.47];
ybs[3540]=['',2.3279591,-0.989886,6.03];
ybs[3541]=['60 Cnc',2.3425852,0.2017305,5.41];
ybs[3542]=['',2.3319038,-0.8305663,5.33];
ybs[3543]=['17 Hya',2.3402647,-0.1402884,6.91];
ybs[3544]=['17 Hya',2.3402718,-0.1403029,6.67];
ybs[3545]=['',2.3387885,-0.3195514,5.75];
ybs[3546]=['σ2 Cnc',2.3475565,0.5732036,5.45];
ybs[3547]=['δ Pyx',2.3399392,-0.4843218,4.89];
ybs[3548]=['',2.3454276,0.0727569,6.14];
ybs[3549]=['',2.3479883,0.2980274,6.17];
ybs[3550]=['',2.3418164,-0.4168912,6.39];
ybs[3551]=['',2.3309362,-1.0545495,5.78];
ybs[3552]=['ο1 Cnc',2.3484335,0.2662424,5.2];
ybs[3553]=['',2.3384666,-0.7873046,6.26];
ybs[3554]=['61 Cnc',2.3519911,0.5264816,6.29];
ybs[3555]=['',2.344768,-0.2928212,5.96];
ybs[3556]=['ο2 Cnc',2.3499151,0.2707543,5.67];
ybs[3557]=['',2.3542511,0.6236743,6.51];
ybs[3558]=['',2.3502719,0.1626548,6.19];
ybs[3559]=['',2.3358692,-1.0176563,6.38];
ybs[3560]=['ι UMa',2.3579833,0.8372839,3.14];
ybs[3561]=['',2.3374059,-0.9605082,5.71];
ybs[3562]=['',2.3363103,-1.0596273,3.84];
ybs[3563]=['α Cnc',2.3537577,0.2057604,4.25];
ybs[3564]=['',2.3520095,0.0257122,6.59];
ybs[3565]=['',2.3424605,-0.921385,4.69];
ybs[3566]=['σ3 Cnc',2.3588777,0.5646087,5.2];
ybs[3567]=['ρ UMa',2.3740865,1.1791416,4.76];
ybs[3568]=['',2.3569237,0.3153099,6.38];
ybs[3569]=['',2.3542239,-0.282768,5.86];
ybs[3570]=['',2.3639679,0.7280383,3.97];
ybs[3571]=['',2.3632685,0.6551134,6.44];
ybs[3572]=['',2.43792,1.4679492,6.33];
ybs[3573]=['',2.3448549,-1.0349363,4.92];
ybs[3574]=['',2.3497012,-0.8489576,5.87];
ybs[3575]=['',2.3582267,-0.3364465,6.18];
ybs[3576]=['',2.3562232,-0.503962,6.25];
ybs[3577]=['',2.3684115,0.6919143,6.36];
ybs[3578]=['66 Cnc',2.3669752,0.5616998,5.82];
ybs[3579]=['',2.3538784,-0.8255994,5.18];
ybs[3580]=['67 Cnc',2.3686567,0.4857807,6.07];
ybs[3581]=['',2.3668718,0.0972383,6.07];
ybs[3582]=['',2.3594467,-0.7212208,4.45];
ybs[3583]=['',2.3792264,0.9462061,5.75];
ybs[3584]=['',2.3605964,-0.7547228,6.07];
ybs[3585]=['κ UMa',2.3771959,0.8218148,3.6];
ybs[3586]=['ν Cnc',2.3725943,0.4255622,5.45];
ybs[3587]=['',2.3686785,-0.0096411,5.67];
ybs[3588]=['',2.3646998,-0.4665833,6.2];
ybs[3589]=['',2.3554835,-1.0324039,5.16];
ybs[3590]=['',2.3722437,0.126156,5.85];
ybs[3591]=['',2.364913,-0.7318835,5.55];
ybs[3592]=['70 Cnc',2.3789037,0.4856911,6.38];
ybs[3593]=['',2.3683217,-0.6889185,6.27];
ybs[3594]=['',2.385,0.8457796,5.95];
ybs[3595]=['',2.3612715,-1.065228,5.79];
ybs[3596]=['',2.366177,-0.9120711,5.23];
ybs[3597]=['',2.3823185,0.5638536,6.46];
ybs[3598]=['',2.3732723,-0.4463573,6.74];
ybs[3599]=['',2.3913816,1.0345154,6.45];
ybs[3600]=['σ1 UMa',2.3994269,1.1659104,5.14];
ybs[3601]=['',2.3620193,-1.1999689,5.88];
ybs[3602]=['',2.3719701,-0.9358388,6.4];
ybs[3603]=['',2.3895242,0.6698789,4.56];
ybs[3604]=['ω Hya',2.386275,0.0876403,4.97];
ybs[3605]=['',2.3769676,-0.8232369,3.75];
ybs[3606]=['α Vol',2.3680586,-1.1600466,4];
ybs[3607]=['σ2 UMa',2.408133,1.1704625,4.8];
ybs[3608]=['',2.3930986,0.3998527,6.4];
ybs[3609]=['',2.390673,0.0242895,6.17];
ybs[3610]=['15 UMa',2.4002693,0.8994212,4.48];
ybs[3611]=['',2.3960496,0.5666931,6.5];
ybs[3612]=['τ Cnc',2.3956795,0.5163165,5.43];
ybs[3613]=['',2.3791786,-1.0109451,6.44];
ybs[3614]=['κ Cnc',2.3941267,0.184948,5.24];
ybs[3615]=['τ UMa',2.4100322,1.1072598,4.67];
ybs[3616]=['',2.3995034,0.5901057,5.93];
ybs[3617]=['75 Cnc',2.3990313,0.4635164,5.98];
ybs[3618]=['ξ Cnc',2.4014061,0.3835145,5.14];
ybs[3619]=['κ Pyx',2.3946637,-0.4525589,4.58];
ybs[3620]=['',2.3870422,-0.9751895,6.11];
ybs[3621]=['19 Hya',2.3979084,-0.1511634,5.6];
ybs[3622]=['',2.3902723,-0.8950582,6.73];
ybs[3623]=['',2.384386,-1.1269675,6.37];
ybs[3624]=['',2.4248032,1.2493523,6.55];
ybs[3625]=['λ Vel',2.3939024,-0.7592853,2.21];
ybs[3626]=['',2.4029863,0.2005827,6.48];
ybs[3627]=['',2.3999613,-0.2169356,5.77];
ybs[3628]=['',2.3975952,-0.4684351,6.15];
ybs[3629]=['',2.3993083,-0.3211457,5.73];
ybs[3630]=['',2.40724,0.5391477,5.95];
ybs[3631]=['79 Cnc',2.4057225,0.3826515,6.01];
ybs[3632]=['20 Hya',2.4017812,-0.1546295,5.46];
ybs[3633]=['',2.3813864,-1.2323679,4.71];
ybs[3634]=['',2.3788065,-1.2683867,4.48];
ybs[3635]=['ε Pyx',2.4028064,-0.5312294,5.59];
ybs[3636]=['',2.4329979,1.2718629,5.96];
ybs[3637]=['',2.404923,-0.4057667,6.53];
ybs[3638]=['',2.4013148,-0.863878,6.48];
ybs[3639]=['16 UMa',2.4247476,1.0707608,5.13];
ybs[3640]=['',2.4122559,0.0941748,6.35];
ybs[3641]=['π1 Cnc',2.4140507,0.260464,6.51];
ybs[3642]=['',2.413482,0.0662288,6.14];
ybs[3643]=['36 Lyn',2.4213919,0.7530171,5.32];
ybs[3644]=['',2.4119648,-0.3459293,5.73];
ybs[3645]=['',2.4072872,-0.7843559,5];
ybs[3646]=['21 Hya',2.4142166,-0.1253561,6.11];
ybs[3647]=['',2.4101551,-0.6864607,6];
ybs[3648]=['',2.4199815,0.3701901,6.48];
ybs[3649]=['',2.4093141,-0.8143052,5.79];
ybs[3650]=['',2.4060579,-1.030427,3.44];
ybs[3651]=['17 UMa',2.4308813,0.9890385,5.27];
ybs[3652]=['',2.4136122,-0.7624689,5.57];
ybs[3653]=['18 UMa',2.4322713,0.9415735,4.83];
ybs[3654]=['',2.4071341,-1.0889018,3.97];
ybs[3655]=['',2.4273682,0.6031883,5.97];
ybs[3656]=['θ Hya',2.4228317,0.0391122,3.88];
ybs[3657]=['',2.4508594,1.2905232,6.5];
ybs[3658]=['',2.4178125,-0.6752556,6.31];
ybs[3659]=['',2.4171564,-0.7390855,6.29];
ybs[3660]=['π2 Cnc',2.4268588,0.2594945,5.34];
ybs[3661]=['',2.4181151,-0.8274876,5.92];
ybs[3662]=['',2.4206983,-0.7717659,5.85];
ybs[3663]=['',2.4145942,-1.0382471,5.54];
ybs[3664]=['',2.4219196,-0.7557394,5.25];
ybs[3665]=['',2.4270681,-0.2635134,6.35];
ybs[3666]=['',2.437728,0.8158212,5.97];
ybs[3667]=['',2.4244909,-0.6575673,5.86];
ybs[3668]=['ζ Oct',2.3293374,-1.4919953,5.42];
ybs[3669]=['',2.4208621,-0.9711509,5.27];
ybs[3670]=['',2.4254782,-0.7963756,6.25];
ybs[3671]=['23 Hya',2.4328295,-0.1121705,5.24];
ybs[3672]=['',2.4273357,-0.6744566,4.94];
ybs[3673]=['24 Hya',2.4327532,-0.153913,5.47];
ybs[3674]=['',2.4279836,-0.6542697,4.62];
ybs[3675]=['β Car',2.4146386,-1.2180652,1.68];
ybs[3676]=['',2.4413174,0.6159234,5.75];
ybs[3677]=['',2.4345414,-0.2556483,5.84];
ybs[3678]=['',2.4291157,-0.7849138,6.04];
ybs[3679]=['',2.438246,0.1994377,6.41];
ybs[3680]=['38 Lyn',2.4431541,0.6410251,3.82];
ybs[3681]=['',2.4250111,-1.0203546,6.02];
ybs[3682]=['',2.4305161,-0.7738713,5.12];
ybs[3683]=['',2.4263427,-1.0062091,6.32];
ybs[3684]=['',2.4331693,-0.6889736,5.33];
ybs[3685]=['',2.4085517,-1.3392869,6.14];
ybs[3686]=['',2.429056,-1.0055722,4.34];
ybs[3687]=['',2.4519111,0.8934534,6.13];
ybs[3688]=['',2.4565277,0.9882733,5.47];
ybs[3689]=['ι Car',2.43281,-1.0358386,2.25];
ybs[3690]=['',2.4358103,-0.9524149,6.33];
ybs[3691]=['',2.4525386,0.6652023,6.12];
ybs[3692]=['',2.4452425,-0.198772,6.62];
ybs[3693]=['',2.4377043,-0.8923049,5.26];
ybs[3694]=['',2.4451147,-0.2776657,5.78];
ybs[3695]=['α Lyn',2.4527188,0.5989522,3.13];
ybs[3696]=['26 Hya',2.4461531,-0.2103067,4.79];
ybs[3697]=['',2.4544147,0.5729353,6.16];
ybs[3698]=['',2.4403306,-0.9011993,5.87];
ybs[3699]=['27 Hya',2.4492995,-0.1680877,4.8];
ybs[3700]=['',2.445754,-0.5965189,6.39];
ybs[3701]=['',2.4531547,0.2669656,6.53];
ybs[3702]=['',2.4326904,-1.2001472,5.39];
ybs[3703]=['',2.4354315,-1.1715509,6.11];
ybs[3704]=['',2.4511063,-0.2738906,6.33];
ybs[3705]=['',2.4499151,-0.5556342,6.82];
ybs[3706]=['',2.4486884,-0.6572256,6.05];
ybs[3707]=['',2.4438063,-0.9644907,6.28];
ybs[3708]=['θ Pyx',2.4533624,-0.4544959,4.72];
ybs[3709]=['',2.4855771,1.3093698,6.29];
ybs[3710]=['',2.4320172,-1.3084447,5.29];
ybs[3711]=['',2.4322087,-1.3056573,5.86];
ybs[3712]=['',2.4745613,1.1146453,6.28];
ybs[3713]=['',2.4632792,0.4382055,6.41];
ybs[3714]=['',2.4596157,-0.1730391,6.53];
ybs[3715]=['',2.4702143,0.8988056,6.31];
ybs[3716]=['',2.4544584,-0.7377547,5.58];
ybs[3717]=['',2.4672534,0.6372369,6.67];
ybs[3718]=['',2.4494232,-1.090476,4.81];
ybs[3719]=['',2.457873,-0.6955165,6.54];
ybs[3720]=['',2.4567273,-0.804996,5.75];
ybs[3721]=['κ Leo',2.4682077,0.4556391,4.46];
ybs[3722]=['',2.4537979,-0.9702321,5.63];
ybs[3723]=['λ Pyx',2.4607675,-0.5045659,4.69];
ybs[3724]=['κ Vel',2.4550513,-0.9614341,2.5];
ybs[3725]=['',2.4628929,-0.6603098,6.48];
ybs[3726]=['',2.471869,0.2881417,6.29];
ybs[3727]=['',2.465135,-0.689435,6.06];
ybs[3728]=['28 Hya',2.4708408,-0.0906472,5.59];
ybs[3729]=['',2.4634348,-0.9043077,6.08];
ybs[3730]=['',2.4605747,-1.0537972,6.3];
ybs[3731]=['',2.4751485,-0.0268842,6.01];
ybs[3732]=['',2.4632594,-1.0772991,5.99];
ybs[3733]=['',2.4862502,0.7945488,5.41];
ybs[3734]=['29 Hya',2.4788085,-0.1623209,6.54];
ybs[3735]=['',2.4762334,-0.5037727,6.1];
ybs[3736]=['',2.4747257,-0.7082269,6.2];
ybs[3737]=['',2.4916202,0.9715926,6.45];
ybs[3738]=['α Hya',2.4803174,-0.1524613,1.98];
ybs[3739]=['',2.478824,-0.391313,4.69];
ybs[3740]=['',2.4812092,-0.1073019,5.38];
ybs[3741]=['',2.5283614,1.4180278,4.29];
ybs[3742]=['',2.4691684,-1.0825656,5.77];
ybs[3743]=['',2.4734491,-0.9329757,5.11];
ybs[3744]=['ω Leo',2.4844462,0.1567245,5.41];
ybs[3745]=['3 Leo',2.4845536,0.1415691,5.71];
ybs[3746]=['',2.4799759,-0.612341,6.65];
ybs[3747]=['23 UMa',2.4996843,1.0992796,3.67];
ybs[3748]=['',2.4867808,-0.0232845,6.27];
ybs[3749]=['τ1 Hya',2.4872407,-0.0496734,4.6];
ybs[3750]=['',2.4883855,-0.0398378,6.14];
ybs[3751]=['',2.4745135,-1.1345724,6.05];
ybs[3752]=['',2.488923,-0.075477,6.26];
ybs[3753]=['',2.4871669,-0.3634789,5.66];
ybs[3754]=['7 LMi',2.4948326,0.5860456,5.85];
ybs[3755]=['ε Ant',2.4869599,-0.6288174,4.51];
ybs[3756]=['',2.4870063,-0.6716217,6.19];
ybs[3757]=['',2.4898262,-0.4088021,6.24];
ybs[3758]=['22 UMa',2.5154713,1.2588503,5.72];
ybs[3759]=['8 LMi',2.4984452,0.6113056,5.37];
ybs[3760]=['',2.4900896,-0.4654286,5.48];
ybs[3761]=['24 UMa',2.5132902,1.2173959,4.56];
ybs[3762]=['',2.4923632,-0.2732264,5.85];
ybs[3763]=['λ Leo',2.498945,0.3995092,4.31];
ybs[3764]=['',2.5212121,1.29571,6.46];
ybs[3765]=['θ UMa',2.5047542,0.9005732,3.17];
ybs[3766]=['',2.4837152,-1.0882144,5.92];
ybs[3767]=['',2.4752144,-1.2510308,5.47];
ybs[3768]=['',2.5058034,0.8615009,6.76];
ybs[3769]=['6 Leo',2.4997339,0.1682134,5.07];
ybs[3770]=['ζ1 Ant',2.4936684,-0.5579638,7];
ybs[3771]=['ζ1 Ant',2.4937194,-0.55793,6.18];
ybs[3772]=['ξ Leo',2.4996976,0.1958575,4.97];
ybs[3773]=['',2.48208,-1.1655116,5.91];
ybs[3774]=['',2.4900847,-0.9004959,5.45];
ybs[3775]=['',2.4980086,-0.1855294,6.14];
ybs[3776]=['ψ Vel',2.493183,-0.7076304,3.6];
ybs[3777]=['τ2 Hya',2.4996325,-0.0220421,4.57];
ybs[3778]=['',2.4992413,-0.1823599,6.13];
ybs[3779]=['ζ2 Ant',2.4970751,-0.557628,5.93];
ybs[3780]=['',2.4970264,-0.6247019,5.87];
ybs[3781]=['9 LMi',2.507048,0.6354502,6.18];
ybs[3782]=['',2.5059744,0.4937499,6.53];
ybs[3783]=['',2.4910329,-1.0199552,5.88];
ybs[3784]=['',2.5027774,0.0311727,6.11];
ybs[3785]=['ι Cha',2.4588585,-1.411319,5.36];
ybs[3786]=['',2.5008583,-0.3399601,5.74];
ybs[3787]=['',2.5109601,0.8172273,6.52];
ybs[3788]=['',2.5005225,-0.5010149,6.46];
ybs[3789]=['26 UMa',2.513345,0.9070949,4.5];
ybs[3790]=['10 LMi',2.5101768,0.633886,4.55];
ybs[3791]=['',2.504098,-0.1498097,6.12];
ybs[3792]=['',2.5035474,-0.2372792,5.94];
ybs[3793]=['',2.4947533,-0.9967945,3.13];
ybs[3794]=['',2.5088323,0.4079786,6.25];
ybs[3795]=['',2.5054309,-0.1268551,6.24];
ybs[3796]=['',2.528861,1.2741087,6.42];
ybs[3797]=['',2.5002615,-0.7108276,5.35];
ybs[3798]=['',2.5046144,-0.369906,5.01];
ybs[3799]=['',2.5139289,0.6901497,4.81];
ybs[3800]=['',2.5055753,-0.4004163,5.91];
ybs[3801]=['',2.5152901,0.6961164,6.76];
ybs[3802]=['',2.5038384,-0.6842925,6.43];
ybs[3803]=['',2.4953769,-1.1658308,6.27];
ybs[3804]=['33 Hya',2.5107341,-0.1046074,5.56];
ybs[3805]=['11 LMi',2.5164234,0.6236308,5.41];
ybs[3806]=['',2.4988079,-1.0972329,6.1];
ybs[3807]=['',2.5061688,-0.8566658,5.12];
ybs[3808]=['7 Leo',2.5169263,0.2495962,6.36];
ybs[3809]=['',2.507842,-0.8959423,5.01];
ybs[3810]=['',2.5209064,0.5424927,5.56];
ybs[3811]=['',2.4946319,-1.2768577,5.47];
ybs[3812]=['',2.5149186,-0.3431741,6.31];
ybs[3813]=['',2.512958,-0.6266186,6.49];
ybs[3814]=['',2.5346278,1.1727281,5.94];
ybs[3815]=['',2.5087084,-1.0351188,4.08];
ybs[3816]=['8 Leo',2.5220314,0.285511,5.69];
ybs[3817]=['10 Leo',2.5225924,0.1179249,5];
ybs[3818]=['',2.5191805,-0.4325247,6.53];
ybs[3819]=['42 Lyn',2.5283221,0.7009274,5.25];
ybs[3820]=['',2.5210963,-0.4428919,5.7];
ybs[3821]=['',2.5178912,-0.8522512,6.17];
ybs[3822]=['34 Hya',2.5251318,-0.1658732,6.4];
ybs[3823]=['',2.5216469,-0.5630051,5.63];
ybs[3824]=['',2.5279788,0.079755,4.68];
ybs[3825]=['',2.5228873,-0.6313748,5.98];
ybs[3826]=['',2.5196367,-0.8627928,4.35];
ybs[3827]=['',2.5192376,-0.9254303,6.19];
ybs[3828]=['',2.5469631,1.2070168,5.69];
ybs[3829]=['27 UMa',2.550496,1.2596353,5.17];
ybs[3830]=['',2.5210966,-0.9380762,5.45];
ybs[3831]=['',2.5154198,-1.1349779,6.56];
ybs[3832]=['',2.5250972,-0.755213,5.5];
ybs[3833]=['',2.563135,1.3622883,6.23];
ybs[3834]=['',2.5280608,-0.6927864,6.7];
ybs[3835]=['ι Hya',2.5339954,-0.0213394,3.91];
ybs[3836]=['37 Hya',2.5335421,-0.18588,6.31];
ybs[3837]=['',2.571291,1.3797681,6.17];
ybs[3838]=['',2.5359173,-0.1893536,6.37];
ybs[3839]=['κ Hya',2.5357302,-0.2515405,5.06];
ybs[3840]=['',2.5421479,0.5445033,5.89];
ybs[3841]=['43 Lyn',2.5441868,0.6925004,5.62];
ybs[3842]=['ο Leo',2.5398264,0.171252,3.52];
ybs[3843]=['13 Leo',2.5422707,0.4508613,6.24];
ybs[3844]=['',2.5475775,0.8438755,6.39];
ybs[3845]=['',2.5495667,0.9474154,6.47];
ybs[3846]=['',2.5300138,-1.0717677,4.52];
ybs[3847]=['13 LMi',2.5471432,0.6110876,6.14];
ybs[3848]=['',2.5398182,-0.4131523,4.77];
ybs[3849]=['',2.55678,1.1327677,6.17];
ybs[3850]=['ζ Cha',2.5015795,-1.4140589,5.11];
ybs[3851]=['15 Leo',2.5506998,0.5217427,5.64];
ybs[3852]=['',2.5439886,-0.4188093,4.94];
ybs[3853]=['',2.5361789,-1.013402,5.32];
ybs[3854]=['',2.5376574,-1.0007692,5.8];
ybs[3855]=['28 UMa',2.5623707,1.1095398,6.34];
ybs[3856]=['ψ Leo',2.5511598,0.2433136,5.35];
ybs[3857]=['',2.5456949,-0.6210268,6.41];
ybs[3858]=['',2.5410926,-0.9650707,6];
ybs[3859]=['',2.554598,0.3278183,6.5];
ybs[3860]=['',2.5645652,0.9956501,5.2];
ybs[3861]=['θ Ant',2.5524745,-0.4860804,4.79];
ybs[3862]=['',2.5485434,-0.8955119,6.15];
ybs[3863]=['ε Leo',2.5605853,0.4135181,2.98];
ybs[3864]=['',2.5524594,-0.6920585,6.82];
ybs[3865]=['',2.5494943,-0.9419967,5.56];
ybs[3866]=['',2.5616574,0.1156669,5.79];
ybs[3867]=['18 Leo',2.5627097,0.2047019,5.63];
ybs[3868]=['',2.5574948,-0.5285548,6.45];
ybs[3869]=['',2.5625655,0.0297424,5.65];
ybs[3870]=['19 Leo',2.5672567,0.2004798,6.45];
ybs[3871]=['',2.5730729,0.8017892,5.09];
ybs[3872]=['',2.5678069,0.1980455,6.02];
ybs[3873]=['',2.5579446,-0.999494,6.46];
ybs[3874]=['',2.5557211,-1.0923824,3.69];
ybs[3875]=['',2.5819926,1.1433812,6.31];
ybs[3876]=['',2.562105,-0.7825435,5.55];
ybs[3877]=['',2.5589365,-1.0275706,6.22];
ybs[3878]=['υ UMa',2.5840974,1.0289776,3.8];
ybs[3879]=['20 Leo',2.5779067,0.3682156,6.09];
ybs[3880]=['υ Car',2.5636224,-1.1371431,3.01];
ybs[3881]=['υ Car',2.5636661,-1.1371528,6.26];
ybs[3882]=['',2.5752593,-0.6504584,5.97];
ybs[3883]=['4 Sex',2.5805314,0.0743724,6.24];
ybs[3884]=['φ UMa',2.5887015,0.9421578,4.59];
ybs[3885]=['',2.5710687,-0.9860042,6.06];
ybs[3886]=['23 Leo',2.5829921,0.2266066,6.46];
ybs[3887]=['',2.5769712,-0.6344417,6.37];
ybs[3888]=['',2.5771265,-0.7996228,5.08];
ybs[3889]=['6 Sex',2.5835879,-0.0755009,6.01];
ybs[3890]=['22 Leo',2.5868985,0.4243344,5.32];
ybs[3891]=['',2.5841099,-0.1093317,6.42];
ybs[3892]=['ν Cha',2.558321,-1.3414152,5.45];
ybs[3893]=['υ1 Hya',2.5844849,-0.2605649,4.12];
ybs[3894]=['',2.5803565,-0.8205989,5.73];
ybs[3895]=['μ Leo',2.5907665,0.4524598,3.88];
ybs[3896]=['7 Sex',2.5879258,0.0413887,6.02];
ybs[3897]=['',2.587874,-0.0001259,6.35];
ybs[3898]=['',2.5867043,-0.290029,6.08];
ybs[3899]=['γ Sex',2.5890834,-0.1429046,5.05];
ybs[3900]=['',2.5831344,-0.8076763,5.62];
ybs[3901]=['',2.6019264,1.065221,6.27];
ybs[3902]=['',2.5846521,-0.8138541,4.58];
ybs[3903]=['',2.5819841,-1.0386165,5.79];
ybs[3904]=['',2.5805624,-1.0965505,5.57];
ybs[3905]=['',2.5945759,0.102542,5.95];
ybs[3906]=['',2.5907566,-0.4784847,6.3];
ybs[3907]=['31 UMa',2.6042512,0.8680644,5.27];
ybs[3908]=['',2.6176095,1.2705169,5.83];
ybs[3909]=['',2.5961909,-0.4540597,4.88];
ybs[3910]=['',2.5900691,-0.9678941,6.48];
ybs[3911]=['',2.5976666,-0.393949,6.24];
ybs[3912]=['',2.6111564,1.0006745,5.93];
ybs[3913]=['',2.5992195,-0.3332322,4.94];
ybs[3914]=['',2.5939115,-0.8941331,5.93];
ybs[3915]=['',2.5961128,-0.7918053,5.71];
ybs[3916]=['',2.606481,0.1544502,5.85];
ybs[3917]=['',2.5983975,-0.8783757,5.72];
ybs[3918]=['19 LMi',2.6125557,0.7150886,5.14];
ybs[3919]=['',2.6138264,0.7911646,6.3];
ybs[3920]=['',2.6040916,-0.7139851,6.41];
ybs[3921]=['',2.607401,-0.4648518,6.28];
ybs[3922]=['',2.606461,-0.5847261,5.84];
ybs[3923]=['',2.6079292,-0.4809917,6.32];
ybs[3924]=['',2.666634,1.4631434,6.37];
ybs[3925]=['',2.604956,-0.8974444,6.37];
ybs[3926]=['',2.6155291,0.4830153,6.3];
ybs[3927]=['ν Leo',2.614344,0.2157336,5.26];
ybs[3928]=['',2.6138583,0.1436422,6.04];
ybs[3929]=['',2.6225971,0.9900812,5.48];
ybs[3930]=['φ Vel',2.6069988,-0.9538494,3.54];
ybs[3931]=['',2.6084738,-0.9201852,6.12];
ybs[3932]=['',2.6206593,0.5159347,5.73];
ybs[3933]=['',2.6109344,-0.8464568,6.05];
ybs[3934]=['',2.6025584,-1.2474397,6.35];
ybs[3935]=['12 Sex',2.620728,0.0576012,6.7];
ybs[3936]=['',2.6176079,-0.419482,6.21];
ybs[3937]=['η Ant',2.6163644,-0.6278879,5.23];
ybs[3938]=['',2.6081041,-1.1270163,6.58];
ybs[3939]=['',2.6065166,-1.2075186,6.2];
ybs[3940]=['π Leo',2.6229571,0.1389221,4.7];
ybs[3941]=['20 LMi',2.6268444,0.5556939,5.36];
ybs[3942]=['',2.6345234,0.3816005,5.66];
ybs[3943]=['',2.623109,-0.9953827,6.52];
ybs[3944]=['',2.6431184,0.9390955,5.74];
ybs[3945]=['',2.6280881,-0.9328653,6.2];
ybs[3946]=['',2.6337283,-0.5351625,6.54];
ybs[3947]=['',2.6291824,-1.0024225,6.2];
ybs[3948]=['',2.6455406,0.91255,6.14];
ybs[3949]=['',2.6378395,-0.1685836,6.12];
ybs[3950]=['',2.6291594,-1.0560235,5.94];
ybs[3951]=['13 Sex',2.6400168,0.0543805,6.45];
ybs[3952]=['',2.6376176,-0.4433468,6.7];
ybs[3953]=['',2.6393024,-0.3174226,5.86];
ybs[3954]=['',2.6356118,-0.8154398,6.12];
ybs[3955]=['',2.6405173,-0.4253529,5.7];
ybs[3956]=['',2.6335413,-1.0517995,6.19];
ybs[3957]=['',2.6326518,-1.0863176,6.42];
ybs[3958]=['',2.6403744,-0.6991999,6.43];
ybs[3959]=['',2.6469263,0.2735254,6.37];
ybs[3960]=['υ2 Hya',2.6440798,-0.2295152,4.6];
ybs[3961]=['',2.6359441,-1.0815643,6.14];
ybs[3962]=['',2.6442289,-0.6365118,6.27];
ybs[3963]=['14 Sex',2.6516126,0.0964384,6.21];
ybs[3964]=['21 LMi',2.6548904,0.6136352,4.48];
ybs[3965]=['η Leo',2.6541453,0.291065,3.52];
ybs[3966]=['',2.6480301,-0.8282588,5.08];
ybs[3967]=['',2.6528938,-0.3006784,5.6];
ybs[3968]=['',2.6475602,-0.9123494,6.52];
ybs[3969]=['',2.6584575,0.5500927,6.24];
ybs[3970]=['31 Leo',2.6565391,0.1729866,4.37];
ybs[3971]=['α Sex',2.6565391,-0.0079895,4.49];
ybs[3972]=['α Leo',2.6586031,0.2073631,1.35];
ybs[3973]=['μ1 Cha',2.6188542,-1.4363919,5.52];
ybs[3974]=['',2.6563298,-0.6530972,6.36];
ybs[3975]=['',2.6599848,-0.1914797,6.53];
ybs[3976]=['',2.6591738,-0.2739799,6.27];
ybs[3977]=['',2.6704888,0.7081617,6.32];
ybs[3978]=['',2.6651177,-0.2126217,6.24];
ybs[3979]=['17 Sex',2.6659699,-0.1482633,5.91];
ybs[3980]=['',2.6599024,-0.9057803,4.86];
ybs[3981]=['',2.6657913,-0.2251934,5.31];
ybs[3982]=['',2.662921,-0.6273249,6.13];
ybs[3983]=['',2.6714255,0.6512729,5.85];
ybs[3984]=['λ Hya',2.6679366,-0.2171327,3.61];
ybs[3985]=['',2.6581447,-1.150198,5.28];
ybs[3986]=['18 Sex',2.669483,-0.1484406,5.65];
ybs[3987]=['μ2 Cha',2.6343221,-1.4250793,6.6];
ybs[3988]=['34 Leo',2.6728644,0.2315732,6.44];
ybs[3989]=['',2.6613139,-1.0757427,5.6];
ybs[3990]=['',2.6710985,-0.1292141,6.25];
ybs[3991]=['',2.6676017,-0.7295758,5.98];
ybs[3992]=['',2.661374,-1.2002577,5.81];
ybs[3993]=['',2.674058,-0.5007923,6.28];
ybs[3994]=['19 Sex',2.6778408,0.0790228,5.77];
ybs[3995]=['',2.6767497,-0.3358122,6.44];
ybs[3996]=['',2.6826391,0.4720867,6.04];
ybs[3997]=['',2.6711842,-1.028258,6.4];
ybs[3998]=['',2.6892085,1.0454179,6.25];
ybs[3999]=['',2.6720459,-1.0148633,5.72];
ybs[4000]=['',2.6749506,-0.9119394,6.16];
ybs[4001]=['',2.6796527,-0.473264,6.25];
ybs[4002]=['',2.6854443,0.3679224,6.02];
ybs[4003]=['',2.6799391,-0.5780373,6.38];
ybs[4004]=['22 LMi',2.6882725,0.5476941,6.46];
ybs[4005]=['',2.6813238,-0.7056898,5.9];
ybs[4006]=['',2.7029098,1.2738332,6.4];
ybs[4007]=['',2.67934,-0.8957113,5.28];
ybs[4008]=['',2.6774052,-1.0472869,6.1];
ybs[4009]=['',2.6821037,-0.7050747,6.35];
ybs[4010]=['',2.679701,-0.9048356,5.78];
ybs[4011]=['',2.7018812,1.2387042,6.66];
ybs[4012]=['',2.6787586,-1.0776712,6.41];
ybs[4013]=['',2.6855277,-0.736692,3.85];
ybs[4014]=['23 LMi',2.6931827,0.510035,5.35];
ybs[4015]=['',2.6791142,-1.1599494,5.16];
ybs[4016]=['32 UMa',2.7021576,1.1348179,5.82];
ybs[4017]=['24 LMi',2.6941676,0.4990726,6.49];
ybs[4018]=['',2.6931277,0.3080955,6.55];
ybs[4019]=['',2.688323,-0.6388878,6.19];
ybs[4020]=['35 Leo',2.69439,0.408674,5.97];
ybs[4021]=['ζ Leo',2.69505,0.4071755,3.44];
ybs[4022]=['',2.6951167,0.4412821,5.84];
ybs[4023]=['λ UMa',2.6971876,0.7474648,3.45];
ybs[4024]=['',2.6922348,-0.1970653,6.08];
ybs[4025]=['37 Leo',2.6948599,0.2380725,5.41];
ybs[4026]=['',2.6889484,-0.7539831,5.6];
ybs[4027]=['ω Car',2.67975,-1.2239163,3.32];
ybs[4028]=['',2.6874953,-0.9610073,6.16];
ybs[4029]=['39 Leo',2.6974581,0.4017438,5.82];
ybs[4030]=['',2.6947578,-0.3623015,6.57];
ybs[4031]=['',2.7015729,0.47695,6.52];
ybs[4032]=['ε Sex',2.6987285,-0.1423637,5.24];
ybs[4033]=['',2.6906617,-1.0470399,6.22];
ybs[4034]=['',2.7055043,0.8145909,6.43];
ybs[4035]=['',2.6937237,-0.8952274,6.3];
ybs[4036]=['',2.7075677,0.8431449,6];
ybs[4037]=['',2.7155983,1.1983237,5.96];
ybs[4038]=['',2.705193,0.4297604,6.4];
ybs[4039]=['',2.7005955,-0.5075414,5.34];
ybs[4040]=['',2.6950985,-1.071982,3.4];
ybs[4041]=['',2.7112223,0.9370799,6.45];
ybs[4042]=['',2.7124281,0.9447197,6];
ybs[4043]=['',2.7026526,-0.6439016,6.3];
ybs[4044]=['40 Leo',2.7082706,0.3382883,4.79];
ybs[4045]=['',2.7058715,-0.220196,6];
ybs[4046]=['',2.7018505,-0.7287871,5.96];
ybs[4047]=['γ1 Leo',2.7093078,0.3447599,2.61];
ybs[4048]=['γ2 Leo',2.7093296,0.3447405,3.8];
ybs[4049]=['',2.7070874,-0.0906547,6.37];
ybs[4050]=['',2.7090155,-0.1596499,6.32];
ybs[4051]=['',2.7020994,-0.980842,5.81];
ybs[4052]=['',2.7578173,1.4689052,5.5];
ybs[4053]=['',2.706457,-0.9619859,4.57];
ybs[4054]=['23 Sex',2.7137104,0.0384175,6.66];
ybs[4055]=['',2.7036287,-1.130355,5.67];
ybs[4056]=['',2.7096021,-0.8340506,5.65];
ybs[4057]=['',2.719277,0.7180401,5.76];
ybs[4058]=['',2.7138739,-0.3154435,6.51];
ybs[4059]=['μ UMa',2.7199509,0.722752,3.05];
ybs[4060]=['42 Leo',2.7173826,0.2598245,6.12];
ybs[4061]=['',2.7153136,-0.4153791,6.5];
ybs[4062]=['',2.7287096,1.1427937,4.97];
ybs[4063]=['',2.7158621,-0.394741,6.51];
ybs[4064]=['',2.7121025,-0.9796808,4.5];
ybs[4065]=['27 LMi',2.7231879,0.590255,5.9];
ybs[4066]=['',2.7185847,-0.3482929,6.13];
ybs[4067]=['43 Leo',2.7223721,0.1126364,6.07];
ybs[4068]=['',2.7256953,0.5153398,6.39];
ybs[4069]=['',2.7233942,0.0978295,6.54];
ybs[4070]=['',2.7187019,-0.7284791,4.83];
ybs[4071]=['28 LMi',2.7277104,0.5869454,5.5];
ybs[4072]=['25 Sex',2.7241411,-0.0726607,5.97];
ybs[4073]=['',2.7228124,-0.5279824,6.27];
ybs[4074]=['',2.762348,1.4393428,5.26];
ybs[4075]=['',2.7276077,0.039775,6.32];
ybs[4076]=['',2.7238525,-0.6649527,5.33];
ybs[4077]=['',2.7245794,-0.7337774,6.27];
ybs[4078]=['44 Leo',2.7322004,0.1517637,5.61];
ybs[4079]=['',2.7204882,-1.1692055,4.99];
ybs[4080]=['30 LMi',2.7354341,0.5882928,4.74];
ybs[4081]=['',2.7248605,-1.0130401,6.35];
ybs[4082]=['',2.7341277,-0.1247753,5.57];
ybs[4083]=['',2.7316176,-0.7427658,6.18];
ybs[4084]=['μ Hya',2.7355418,-0.2954113,3.81];
ybs[4085]=['',2.7298274,-1.0239083,5.95];
ybs[4086]=['',2.7423424,0.7245063,6.02];
ybs[4087]=['',2.7399877,0.3364095,6.15];
ybs[4088]=['',2.7451095,0.849887,6.44];
ybs[4089]=['',2.7354205,-0.7474953,6.13];
ybs[4090]=['β LMi',2.744064,0.6390955,4.21];
ybs[4091]=['45 Leo',2.742667,0.1688222,6.04];
ybs[4092]=['',2.7259399,-1.2936514,4];
ybs[4093]=['',2.7473953,0.7875336,6.35];
ybs[4094]=['α Ant',2.7399761,-0.543799,4.25];
ybs[4095]=['',2.7274557,-1.2926053,6.19];
ybs[4096]=['35 UMa',2.7538181,1.1438192,6.32];
ybs[4097]=['',2.7379629,-0.9593559,5.58];
ybs[4098]=['',2.7560735,1.1199309,6.12];
ybs[4099]=['',2.7472368,-0.0668876,6.05];
ybs[4100]=['',2.7404383,-1.0075529,4.66];
ybs[4101]=['',2.7434452,-0.863856,6.1];
ybs[4102]=['36 UMa',2.7564918,0.9754707,4.84];
ybs[4103]=['32 LMi',2.7537938,0.6778016,5.77];
ybs[4104]=['',2.7424555,-1.0267625,3.82];
ybs[4105]=['',2.7400776,-1.1483281,6.01];
ybs[4106]=['δ Sex',2.7504992,-0.0493783,5.21];
ybs[4107]=['',2.750181,-0.5192984,5.58];
ybs[4108]=['δ Ant',2.7506324,-0.5357678,5.56];
ybs[4109]=['β Sex',2.754072,-0.0126899,5.09];
ybs[4110]=['',2.7465467,-1.1215852,5.29];
ybs[4111]=['',2.7826957,1.4033033,6.52];
ybs[4112]=['',2.7569894,-0.1348745,6.2];
ybs[4113]=['',2.7570005,-0.2387362,5.58];
ybs[4114]=['33 LMi',2.76131,0.5635503,5.9];
ybs[4115]=['',2.7562245,-0.4638057,6.51];
ybs[4116]=['',2.7774138,1.3198551,4.84];
ybs[4117]=['46 Leo',2.7625506,0.2451626,5.46];
ybs[4118]=['',2.7544764,-1.0724399,6.43];
ybs[4119]=['',2.7519113,-1.170681,6.19];
ybs[4120]=['',2.7603628,-0.4944145,6.05];
ybs[4121]=['',2.7699059,0.9321246,6.45];
ybs[4122]=['',2.7674395,0.7039776,4.75];
ybs[4123]=['ρ Leo',2.7651778,0.1608519,3.85];
ybs[4124]=['',2.7578839,-0.9390891,4.89];
ybs[4125]=['',2.7607147,-0.7881392,5.74];
ybs[4126]=['',2.7606492,-0.7881877,6.09];
ybs[4127]=['34 LMi',2.7685796,0.6090843,5.58];
ybs[4128]=['',2.7522733,-1.2580885,4.74];
ybs[4129]=['',2.763354,-0.7803257,5.91];
ybs[4130]=['',2.7604553,-1.0781888,3.32];
ybs[4131]=['37 UMa',2.7762718,0.9946957,5.16];
ybs[4132]=['',2.7551796,-1.2795337,4.93];
ybs[4133]=['',2.7649971,-0.8219432,5.02];
ybs[4134]=['',2.7639599,-1.025511,6];
ybs[4135]=['44 Hya',2.7700524,-0.4160166,5.08];
ybs[4136]=['48 Leo',2.7738255,0.1197779,5.08];
ybs[4137]=['',2.7667276,-1.0171934,6.14];
ybs[4138]=['49 Leo',2.7748765,0.1493896,5.67];
ybs[4139]=['',2.7741931,-0.4060852,6.1];
ybs[4140]=['35 LMi',2.7809863,0.632435,6.28];
ybs[4141]=['',2.7700793,-1.0660211,6.23];
ybs[4142]=['',2.7772475,-0.3256808,6.49];
ybs[4143]=['',2.7750553,-0.6920871,5.38];
ybs[4144]=['',2.7748018,-0.7636795,6.08];
ybs[4145]=['',2.780139,-0.1863034,6.57];
ybs[4146]=['φ2 Hya',2.7800237,-0.2868538,6.03];
ybs[4147]=['',2.7790165,-0.4671554,6.29];
ybs[4148]=['',2.7812121,-0.2150487,5.7];
ybs[4149]=['',2.7762297,-1.00616,4.45];
ybs[4150]=['',2.784069,-0.2066437,6.52];
ybs[4151]=['',2.7563553,-1.4313691,7.07];
ybs[4152]=['',2.7840359,-0.4800302,4.89];
ybs[4153]=['',2.7856225,-0.2351953,4.82];
ybs[4154]=['',2.7794321,-1.0411898,5.08];
ybs[4155]=['',2.7932845,0.9350921,5.52];
ybs[4156]=['37 LMi',2.7912211,0.5564925,4.71];
ybs[4157]=['',2.7839964,-0.8432915,3.84];
ybs[4158]=['38 LMi',2.7930839,0.6600574,5.85];
ybs[4159]=['',2.7843143,-1.0266822,5.45];
ybs[4160]=['',2.7739048,-1.3334326,6.3];
ybs[4161]=['φ3 Hya',2.7900799,-0.2961488,4.91];
ybs[4162]=['',2.7912489,-0.218778,6.04];
ybs[4163]=['',2.786972,-1.0009063,5.91];
ybs[4164]=['γ Cha',2.7736155,-1.3735508,4.11];
ybs[4165]=['',2.7908281,-0.7477872,6.11];
ybs[4166]=['',2.8057422,1.1929576,5.75];
ybs[4167]=['',2.790006,-1.0345348,4.66];
ybs[4168]=['38 UMa',2.8061668,1.1453632,5.12];
ybs[4169]=['',2.7910629,-1.0281455,5.92];
ybs[4170]=['',2.7925607,-0.9720582,4.28];
ybs[4171]=['',2.8112787,1.2039986,5];
ybs[4172]=['33 Sex',2.8025458,-0.0320003,6.26];
ybs[4173]=['',2.7997799,-0.6254109,6.37];
ybs[4174]=['',2.8063357,0.5516114,6.02];
ybs[4175]=['',2.796014,-1.1378181,5.52];
ybs[4176]=['',2.7911684,-1.3017554,6.07];
ybs[4177]=['39 UMa',2.8135337,0.9967054,5.8];
ybs[4178]=['',2.8011068,-1.0431612,6.42];
ybs[4179]=['40 LMi',2.8099363,0.457861,5.51];
ybs[4180]=['',2.807303,-0.2455151,6.24];
ybs[4181]=['',2.8124988,0.804802,5.18];
ybs[4182]=['41 LMi',2.8115853,0.4031053,5.08];
ybs[4183]=['35 Sex',2.8110959,0.0812571,5.79];
ybs[4184]=['',2.8079483,-0.5726047,5.64];
ybs[4185]=['',2.8198475,1.1749389,6];
ybs[4186]=['',2.8050024,-1.126755,4.82];
ybs[4187]=['',2.8151447,0.3432433,6.27];
ybs[4188]=['',2.807177,-1.0351167,5.38];
ybs[4189]=['θ Car',2.808173,-1.1255011,2.76];
ybs[4190]=['',2.8108801,-1.0586953,4.57];
ybs[4191]=['36 Sex',2.8189636,0.0418132,6.28];
ybs[4192]=['41 UMa',2.8250875,0.999608,6.34];
ybs[4193]=['42 LMi',2.8223493,0.5338924,5.24];
ybs[4194]=['',2.812106,-1.1229629,5.77];
ybs[4195]=['',2.8132662,-1.1179408,4.82];
ybs[4196]=['',2.801368,-1.3940847,5.97];
ybs[4197]=['',2.8230942,0.1096169,6.37];
ybs[4198]=['51 Leo',2.824584,0.3281023,5.49];
ybs[4199]=['52 Leo',2.8245946,0.24613,5.48];
ybs[4200]=['η Car',2.817587,-1.0432964,6.21];
ybs[4201]=['',2.8137171,-1.2383495,6.26];
ybs[4202]=['',2.8146454,-1.2382627,6.46];
ybs[4203]=['',2.8140813,-1.2659939,6.27];
ybs[4204]=['',2.8262362,-0.3034994,5.42];
ybs[4205]=['',2.8360929,1.1351513,6.39];
ybs[4206]=['μ Vel',2.8253585,-0.8641569,2.69];
ybs[4207]=['',2.8228903,-1.0593417,6.25];
ybs[4208]=['',2.8296033,-0.2679886,6.67];
ybs[4209]=['',2.8226947,-1.1276132,5.34];
ybs[4210]=['',2.8236705,-1.1232213,5.23];
ybs[4211]=['',2.8259913,-0.992216,5.23];
ybs[4212]=['',2.8252333,-1.1253165,4.85];
ybs[4213]=['43 LMi',2.8357966,0.5117826,6.15];
ybs[4214]=['',2.8342813,-0.0358089,5.93];
ybs[4215]=['',2.8320574,-0.5546845,5.88];
ybs[4216]=['',2.8289726,-1.0046191,6.36];
ybs[4217]=['53 Leo',2.8369252,0.1824287,5.34];
ybs[4218]=['',2.8308361,-1.0474049,6];
ybs[4219]=['40 Sex',2.8369318,-0.0718561,6.61];
ybs[4220]=['44 LMi',2.8398851,0.4866139,6.04];
ybs[4221]=['δ1 Cha',2.816141,-1.4060726,5.47];
ybs[4222]=['ν Hya',2.8382885,-0.2842537,3.11];
ybs[4223]=['',2.8387846,-0.1735855,5.86];
ybs[4224]=['δ2 Cha',2.8183997,-1.4073052,4.45];
ybs[4225]=['43 UMa',2.8459993,0.9859206,5.67];
ybs[4226]=['42 UMa',2.8469976,1.0337033,5.58];
ybs[4227]=['41 Sex',2.8413102,-0.1569189,5.79];
ybs[4228]=['',2.8395109,-0.5960477,5.61];
ybs[4229]=['',2.8365962,-1.0370135,5.91];
ybs[4230]=['',2.8448022,-0.0555994,5.95];
ybs[4231]=['',2.8516788,0.9158088,6.65];
ybs[4232]=['',2.8517572,0.9147325,6.44];
ybs[4233]=['',2.8567071,1.2175498,5.93];
ybs[4234]=['',2.8498049,0.0162669,6.38];
ybs[4235]=['',2.8514233,-0.0051433,6.31];
ybs[4236]=['44 UMa',2.8563494,0.9510574,5.1];
ybs[4237]=['46 LMi',2.8548522,0.5955344,3.83];
ybs[4238]=['ω UMa',2.8578817,0.7521763,4.71];
ybs[4239]=['',2.8549549,-0.0409921,6.12];
ybs[4240]=['',2.850283,-1.0006642,5.25];
ybs[4241]=['',2.855131,-0.3531201,5.24];
ybs[4242]=['',2.8554222,-0.2712062,6.38];
ybs[4243]=['',2.8563232,-0.0387918,5.45];
ybs[4244]=['48 LMi',2.8608215,0.4432661,6.2];
ybs[4245]=['',2.8587039,-0.2417553,5.66];
ybs[4246]=['',2.8620779,0.5923846,5.72];
ybs[4247]=['',2.8545137,-1.0288145,3.78];
ybs[4248]=['46 UMa',2.8654277,0.5831715,5.03];
ybs[4249]=['54 Leo',2.8647816,0.4303294,4.5];
ybs[4250]=['54 Leo',2.864818,0.4303148,6.3];
ybs[4251]=['',2.8625555,-0.3623061,6.44];
ybs[4252]=['',2.8548333,-1.2359321,5.99];
ybs[4253]=['',2.8615327,-0.7390544,6.11];
ybs[4254]=['',2.8677246,0.7315476,6.03];
ybs[4255]=['55 Leo',2.8649792,0.0112272,5.91];
ybs[4256]=['',2.8587871,-1.0807111,5.93];
ybs[4257]=['56 Leo',2.8664127,0.106318,5.81];
ybs[4258]=['',2.8481127,-1.3902017,6.33];
ybs[4259]=['',2.8676726,0.3884741,6.14];
ybs[4260]=['50 LMi',2.8689746,0.4434223,6.35];
ybs[4261]=['',2.8623059,-1.0578537,5.92];
ybs[4262]=['',2.8854759,1.3556989,6.2];
ybs[4263]=['ι Ant',2.8690383,-0.6498134,4.6];
ybs[4264]=['',2.8706296,-0.887654,5.91];
ybs[4265]=['',2.8812069,0.9038735,6.17];
ybs[4266]=['',2.8733584,-1.044158,6.11];
ybs[4267]=['47 UMa',2.8817543,0.7039991,5.05];
ybs[4268]=['',2.8820488,0.6283001,6];
ybs[4269]=['',2.8700319,-1.312375,6.13];
ybs[4270]=['',2.8852179,0.7929366,5.47];
ybs[4271]=['',2.8824229,0.2026625,6.55];
ybs[4272]=['',2.8800582,-0.5904674,5.71];
ybs[4273]=['',2.8861125,0.8972342,6.43];
ybs[4274]=['',2.881461,-0.2870716,5.89];
ybs[4275]=['',2.8856077,0.7473009,6.02];
ybs[4276]=['',2.8893463,1.1052616,6.39];
ybs[4277]=['α Crt',2.8825796,-0.3210188,4.08];
ybs[4278]=['49 UMa',2.8877227,0.6827374,5.08];
ybs[4279]=['',2.8844406,-0.2474442,5.88];
ybs[4280]=['',2.8795501,-1.0718824,6.16];
ybs[4281]=['58 Leo',2.8861772,0.0614929,4.84];
ybs[4282]=['',2.8832429,-0.7662235,5.81];
ybs[4283]=['',2.8839867,-0.7386234,4.39];
ybs[4284]=['59 Leo',2.8870105,0.1048445,4.99];
ybs[4285]=['β UMa',2.8923779,0.9824133,2.37];
ybs[4286]=['',2.8837765,-0.9060343,6.15];
ybs[4287]=['',2.8877445,-0.2772811,6.34];
ybs[4288]=['',2.8864049,-0.5573477,6.07];
ybs[4289]=['61 Leo',2.8916657,-0.0450134,4.74];
ybs[4290]=['60 Leo',2.8940273,0.3505549,4.42];
ybs[4291]=['α UMa',2.9007353,1.076105,1.79];
ybs[4292]=['',2.893995,-0.4699439,6.23];
ybs[4293]=['',2.8978527,-0.014783,6.14];
ybs[4294]=['',2.877267,-1.4250639,6.71];
ybs[4295]=['',2.8977973,-0.1989346,5.5];
ybs[4296]=['62 Leo',2.8994581,-0.0016645,5.95];
ybs[4297]=['',2.8977091,-0.5594712,6.46];
ybs[4298]=['',2.8993526,-0.2361253,6.34];
ybs[4299]=['51 UMa',2.903747,0.6657864,6];
ybs[4300]=['χ Leo',2.9056474,0.1263868,4.63];
ybs[4301]=['',2.9029838,-0.83381,5.67];
ybs[4302]=['η Oct',2.8755892,-1.4780827,6.19];
ybs[4303]=['',2.9048112,-0.6265626,5.43];
ybs[4304]=['χ1 Hya',2.9067645,-0.4780164,4.94];
ybs[4305]=['',2.9079178,-0.1951911,6.09];
ybs[4306]=['',2.9053605,-0.8637143,6.13];
ybs[4307]=['χ2 Hya',2.9095085,-0.4779157,5.71];
ybs[4308]=['',2.9098181,-0.895481,6.3];
ybs[4309]=['65 Leo',2.9138411,0.0324751,5.52];
ybs[4310]=['',2.9124872,-0.5030496,6.77];
ybs[4311]=['',2.9113966,-0.8910165,6.32];
ybs[4312]=['64 Leo',2.9173013,0.4054167,6.46];
ybs[4313]=['',2.9113743,-1.0257317,6.02];
ybs[4314]=['',2.9146139,-0.5704104,6.59];
ybs[4315]=['',2.9114938,-1.0911622,4.61];
ybs[4316]=['',2.9108113,-1.1333214,6.41];
ybs[4317]=['',2.9151023,-0.7458405,5.15];
ybs[4318]=['',2.9179847,-0.5283057,6.54];
ybs[4319]=['',2.9123915,-1.2387108,5.57];
ybs[4320]=['',2.9267486,1.1713801,6.06];
ybs[4321]=['',2.9195443,-0.5247817,6.49];
ybs[4322]=['67 Leo',2.922355,0.4287101,5.68];
ybs[4323]=['',2.9246368,0.6320596,5.74];
ybs[4324]=['',2.9216058,-0.491757,5.44];
ybs[4325]=['ψ UMa',2.926231,0.7749869,3.01];
ybs[4326]=['',2.9261208,0.7524528,5.89];
ybs[4327]=['',2.9205622,-1.0309664,3.91];
ybs[4328]=['',2.9203794,-1.0828414,5.13];
ybs[4329]=['',2.9266277,-0.56658,5.81];
ybs[4330]=['',2.9377383,1.1899058,6.4];
ybs[4331]=['',2.9349714,0.2496687,6.3];
ybs[4332]=['',2.9307325,-1.0218993,6.88];
ybs[4333]=['β Crt',2.9344269,-0.4000493,4.48];
ybs[4334]=['',2.939792,0.9564187,6.63];
ybs[4335]=['',2.9386577,0.6234006,6.41];
ybs[4336]=['',2.9369209,-0.5677425,6.38];
ybs[4337]=['ψ Crt',2.9381611,-0.3245506,6.13];
ybs[4338]=['',2.9384444,-0.3812594,6.4];
ybs[4339]=['',2.932789,-1.2484632,6.35];
ybs[4340]=['',2.9380729,-0.8586408,5.36];
ybs[4341]=['',2.9436476,0.715465,6.33];
ybs[4342]=['',2.9380756,-1.0544038,4.6];
ybs[4343]=['',2.939798,-0.869729,6.11];
ybs[4344]=['',2.9411643,-0.7761072,5.8];
ybs[4345]=['',2.9386386,-1.1216379,5.23];
ybs[4346]=['69 Leo',2.9437454,-0.0028835,5.42];
ybs[4347]=['δ Leo',2.945391,0.3565375,2.56];
ybs[4348]=['',2.9449714,0.1390162,5.79];
ybs[4349]=['θ Leo',2.9459324,0.2676273,3.34];
ybs[4350]=['',2.9428339,-0.9307342,5.76];
ybs[4351]=['',2.9420891,-1.0422217,5.74];
ybs[4352]=['72 Leo',2.9501829,0.4014248,4.63];
ybs[4353]=['',2.9542317,0.9193936,6.5];
ybs[4354]=['',2.9484046,-0.7649734,6.21];
ybs[4355]=['73 Leo',2.9530068,0.23059,5.32];
ybs[4356]=['',2.9534329,0.2225129,6.67];
ybs[4357]=['',2.9569358,0.861855,5.88];
ybs[4358]=['φ Leo',2.9563866,-0.0654043,4.47];
ybs[4359]=['',2.9577128,-0.1261956,6.14];
ybs[4360]=['',2.9552057,-0.8024275,6.31];
ybs[4361]=['75 Leo',2.9591593,0.0334192,5.18];
ybs[4362]=['',2.9584966,-0.6651486,6.27];
ybs[4363]=['',2.9605191,-0.607951,6.45];
ybs[4364]=['ξ UMa',2.9632262,0.548615,4.87];
ybs[4365]=['ξ UMa',2.9632333,0.548615,4.41];
ybs[4366]=['',2.9607811,-0.6393185,6.68];
ybs[4367]=['ν UMa',2.9645308,0.575929,3.48];
ybs[4368]=['',2.963838,0.2074999,6.66];
ybs[4369]=['',2.9584972,-1.1854168,6.06];
ybs[4370]=['55 UMa',2.9674185,0.6647896,4.78];
ybs[4371]=['76 Leo',2.9662547,0.0271338,5.91];
ybs[4372]=['δ Crt',2.9680221,-0.2596097,3.56];
ybs[4373]=['',2.9755586,1.1694492,6.21];
ybs[4374]=['',2.9671817,-1.1288515,5.99];
ybs[4375]=['',2.9630663,-1.3921527,6.35];
ybs[4376]=['σ Leo',2.9759626,0.1035571,4.05];
ybs[4377]=['',2.9681308,-1.3131585,6.27];
ybs[4378]=['',2.9793507,0.9944692,6.43];
ybs[4379]=['',2.970355,-1.2582152,6.41];
ybs[4380]=['π Cen',2.9749655,-0.9527257,3.89];
ybs[4381]=['',2.9839924,1.1211012,6.02];
ybs[4382]=['56 UMa',2.9835636,0.7572389,4.99];
ybs[4383]=['',2.9811079,-0.7808949,6.12];
ybs[4384]=['',2.9853726,0.0006188,6.05];
ybs[4385]=['λ Crt',2.9855646,-0.3294521,5.09];
ybs[4386]=['',2.9847901,-0.6328726,5];
ybs[4387]=['',2.978215,-1.3561983,6.43];
ybs[4388]=['',2.9842401,-0.9926673,5.79];
ybs[4389]=['ι Leo',2.988148,0.1820887,3.94];
ybs[4390]=['79 Leo',2.9886,0.0228903,5.39];
ybs[4391]=['',2.985074,-1.1353578,5.11];
ybs[4392]=['ε Crt',2.9910383,-0.1912138,4.83];
ybs[4393]=['',2.9897934,-0.7463978,6.12];
ybs[4394]=['',2.9927612,0.1978149,5.8];
ybs[4395]=['γ Crt',2.9921915,-0.3103231,4.08];
ybs[4396]=['',2.9884123,-1.2627968,5.59];
ybs[4397]=['',2.9973074,0.9730939,5.75];
ybs[4398]=['81 Leo',2.9955108,0.2855364,5.57];
ybs[4399]=['',2.9947441,-0.6311008,5.22];
ybs[4400]=['80 Leo',2.9964425,0.0656876,6.37];
ybs[4401]=['',2.9950019,-0.6605047,5.89];
ybs[4402]=['',2.999169,0.5821396,6.32];
ybs[4403]=['',2.9954173,-1.1182175,5.17];
ybs[4404]=['83 Leo',3.0004601,0.0509046,6.5];
ybs[4405]=['',2.999277,-1.0683456,5.3];
ybs[4406]=['κ Crt',3.0021546,-0.217348,5.94];
ybs[4407]=['',3.0002686,-0.9295001,5.81];
ybs[4408]=['τ Leo',3.005615,0.0481642,4.95];
ybs[4409]=['',3.0054146,-0.0313548,6.25];
ybs[4410]=['',3.0056048,-0.6182849,6.45];
ybs[4411]=['',3.0109963,1.0765498,5.83];
ybs[4412]=['57 UMa',3.010733,0.6848738,5.31];
ybs[4413]=['',3.0082207,-0.7464896,5.08];
ybs[4414]=['',3.0137419,0.9885701,6.28];
ybs[4415]=['',3.0064674,-1.2666023,6.09];
ybs[4416]=['85 Leo',3.0133509,0.2673274,5.74];
ybs[4417]=['',3.0158465,0.9471035,6.41];
ybs[4418]=['',3.0129493,-0.4286517,5.76];
ybs[4419]=['',3.0239347,1.4142489,6.15];
ybs[4420]=['',3.0166524,0.8126403,6.35];
ybs[4421]=['58 UMa',3.0170692,0.75183,5.94];
ybs[4422]=['87 Leo',3.0159694,-0.0541095,4.77];
ybs[4423]=['86 Leo',3.0167868,0.3196235,5.52];
ybs[4424]=['λ Dra',3.0212793,1.2083685,3.84];
ybs[4425]=['',3.0186879,0.8348346,6.42];
ybs[4426]=['',3.0199507,0.8498442,6.56];
ybs[4427]=['88 Leo',3.0222875,0.249019,6.2];
ybs[4428]=['',3.019677,-1.0711962,6.38];
ybs[4429]=['',3.025209,1.0644022,5.48];
ybs[4430]=['',3.0223532,-0.3643092,6.24];
ybs[4431]=['ο1 Cen',3.0219728,-1.0391504,5.13];
ybs[4432]=['ο2 Cen',3.0221687,-1.0404352,5.15];
ybs[4433]=['',3.0243905,-0.5124299,5.81];
ybs[4434]=['',3.0244051,-0.5123911,5.64];
ybs[4435]=['',3.0249259,-0.4685059,6.16];
ybs[4436]=['',3.02676,-0.1383045,5.95];
ybs[4437]=['',3.0266619,-0.707437,5.64];
ybs[4438]=['',3.0243044,-1.1703997,5.9];
ybs[4439]=['',3.0271482,-0.5442633,5.04];
ybs[4440]=['ξ Hya',3.0275813,-0.5577121,3.54];
ybs[4441]=['',3.0287011,-0.2858385,6.05];
ybs[4442]=['',3.0319323,0.6408629,6.4];
ybs[4443]=['',3.0302441,-0.7100654,5.39];
ybs[4444]=['',3.0328235,0.1907083,6.55];
ybs[4445]=['89 Leo',3.0336694,0.0517169,5.77];
ybs[4446]=['90 Leo',3.0352058,0.2914715,5.95];
ybs[4447]=['',3.0370354,0.9544927,5.63];
ybs[4448]=['',3.0340851,-0.5747061,5.98];
ybs[4449]=['',3.0367668,0.3550788,6.45];
ybs[4450]=['',3.0351272,-0.9487789,4.62];
ybs[4451]=['2 Dra',3.0414491,1.2082192,5.2];
ybs[4452]=['',3.0359792,-0.8592872,5.5];
ybs[4453]=['',3.0371867,-0.8284969,5.71];
ybs[4454]=['',3.0396131,0.1887435,6.56];
ybs[4455]=['',3.0421786,0.4831802,5.8];
ybs[4456]=['',3.040274,-0.8331953,5.25];
ybs[4457]=['λ Cen',3.0394794,-1.1015929,3.13];
ybs[4458]=['θ Crt',3.0437313,-0.1727731,4.7];
ybs[4459]=['',3.0432185,-0.5875989,5.74];
ybs[4460]=['',3.0436244,-0.6516138,6.31];
ybs[4461]=['υ Leo',3.0449218,-0.0160718,4.3];
ybs[4462]=['',3.0420966,-1.067254,5.83];
ybs[4463]=['',3.0451284,-0.5774424,6.29];
ybs[4464]=['',3.049199,0.8817637,6.14];
ybs[4465]=['',3.0448875,-1.0712881,5.15];
ybs[4466]=['',3.0474352,-0.8350388,5.44];
ybs[4467]=['59 UMa',3.0511605,0.7597164,5.59];
ybs[4468]=['',3.0502503,0.153365,6.17];
ybs[4469]=['π Cha',3.0456224,-1.326339,5.65];
ybs[4470]=['60 UMa',3.0521146,0.8157171,6.1];
ybs[4471]=['',3.0534096,1.1213725,6.46];
ybs[4472]=['',3.051958,0.5851834,6.27];
ybs[4473]=['ω Vir',3.0515424,0.1402748,5.36];
ybs[4474]=['',3.0512588,-0.0442114,6.22];
ybs[4475]=['',3.0482765,-1.1818893,5.96];
ybs[4476]=['',3.0529326,0.7856003,6.44];
ybs[4477]=['',3.0497464,-1.0807671,5.15];
ybs[4478]=['ι Crt',3.0523928,-0.2321108,5.48];
ybs[4479]=['',3.0538352,-0.4331584,6.42];
ybs[4480]=['',3.0575547,-0.2542191,6.21];
ybs[4481]=['',3.0574979,-0.2917728,6.19];
ybs[4482]=['',3.0556884,-1.1431004,5.17];
ybs[4483]=['',3.0604578,1.0100825,6.37];
ybs[4484]=['ο Hya',3.0590709,-0.6081042,4.7];
ybs[4485]=['92 Leo',3.0617191,0.3709816,5.26];
ybs[4486]=['61 UMa',3.0629142,0.5952368,5.33];
ybs[4487]=['',3.0611452,-0.9436247,5.96];
ybs[4488]=['',3.0631375,-0.5112681,6.44];
ybs[4489]=['',3.0618756,-1.0853697,4.94];
ybs[4490]=['',3.065963,0.9612465,6.27];
ybs[4491]=['62 UMa',3.0651797,0.552379,5.73];
ybs[4492]=['',3.0639162,-0.7538593,5.55];
ybs[4493]=['',3.0657185,-0.5689225,5.22];
ybs[4494]=['3 Dra',3.0693249,1.1632243,5.3];
ybs[4495]=['',3.0673972,0.3859567,6.59];
ybs[4496]=['',3.0671728,-0.3558907,6.22];
ybs[4497]=['',3.0615391,-1.4520634,6.33];
ybs[4498]=['',3.0732169,-0.650789,5.98];
ybs[4499]=['',3.0703363,-1.3858535,6.39];
ybs[4500]=['',3.0753229,-0.118236,6.07];
ybs[4501]=['',3.0733696,-1.0923429,5.03];
ybs[4502]=['',3.0767063,0.4384463,6.02];
ybs[4503]=['',3.074949,-1.0991304,6.1];
ybs[4504]=['ζ Crt',3.0789872,-0.3219794,4.73];
ybs[4505]=['ξ Vir',3.0813152,0.142438,4.85];
ybs[4506]=['',3.0808496,-0.8581253,6.26];
ybs[4507]=['ν Vir',3.0838203,0.1122629,4.03];
ybs[4508]=['χ UMa',3.0847479,0.8322111,3.71];
ybs[4509]=['',3.0831436,-0.7991382,5.29];
ybs[4510]=['λ Mus',3.0824619,-1.1663313,3.64];
ybs[4511]=['',3.0885995,0.9691997,5.27];
ybs[4512]=['',3.0864832,-1.069461,4.11];
ybs[4513]=['',3.0865986,-0.7085657,4.91];
ybs[4514]=['',3.0892251,-0.6283923,6.17];
ybs[4515]=['',3.0898713,-0.5303049,6.48];
ybs[4516]=['',3.090028,-1.00869,5.41];
ybs[4517]=['93 Leo',3.0931123,0.3511879,4.53];
ybs[4518]=['4 Vir',3.0927878,0.1422186,5.32];
ybs[4519]=['',3.0948401,-0.1817,6.26];
ybs[4520]=['μ Mus',3.0939895,-1.1678353,4.72];
ybs[4521]=['',3.0959811,0.2476072,5.88];
ybs[4522]=['',3.0963855,-0.4685693,5.11];
ybs[4523]=['',3.0975964,-0.0072594,6.15];
ybs[4524]=['β Leo',3.0977916,0.2526297,2.14];
ybs[4525]=['',3.0986154,0.2817913,6.04];
ybs[4526]=['',3.100593,0.6079737,5.7];
ybs[4527]=['',3.1003493,-1.1150153,4.32];
ybs[4528]=['',3.1014265,-1.2273709,4.97];
ybs[4529]=['',3.1032715,-0.2785762,6.13];
ybs[4530]=['β Vir',3.1049074,0.029101,3.61];
ybs[4531]=['',3.1037201,-1.0951382,5.7];
ybs[4532]=['',3.1045367,-0.4777861,6.48];
ybs[4533]=['',3.1059084,0.2126078,6.35];
ybs[4534]=['',3.1063903,-0.0947835,5.64];
ybs[4535]=['',3.1069607,0.5808043,6.27];
ybs[4536]=['',3.106806,-0.7901275,4.46];
ybs[4537]=['',3.1078157,-0.2144163,6.35];
ybs[4538]=['',3.1092266,-0.5398717,5.85];
ybs[4539]=['',3.1098355,-1.1397608,4.9];
ybs[4540]=['',3.1149132,0.6566143,6.45];
ybs[4541]=['',3.1112495,-0.996324,5.57];
ybs[4542]=['β Hya',3.114537,-0.5935069,4.28];
ybs[4543]=['',3.1168799,-0.6137286,6.17];
ybs[4544]=['γ UMa',3.118642,0.9354498,2.44];
ybs[4545]=['',3.1186212,0.0079333,6.3];
ybs[4546]=['',3.1201001,-1.0036935,6.06];
ybs[4547]=['',3.1211729,-0.6605424,6.46];
ybs[4548]=['',3.1223996,-0.4504921,5.3];
ybs[4549]=['6 Vir',3.1239199,0.1456736,5.58];
ybs[4550]=['65 UMa',3.1241388,0.8094756,6.54];
ybs[4551]=['65 UMa',3.124538,0.8093496,7.03];
ybs[4552]=['',3.1247388,0.6398199,6.49];
ybs[4553]=['',3.1236076,-1.1061251,5.91];
ybs[4554]=['95 Leo',3.1266428,0.2713857,5.53];
ybs[4555]=['',3.12659,-0.4987166,5.93];
ybs[4556]=['66 UMa',3.1279766,0.9861319,5.84];
ybs[4557]=['η Crt',3.1281102,-0.3010387,5.18];
ybs[4558]=['',3.1276449,-0.6944068,6.13];
ybs[4559]=['',3.1319581,1.0725353,6.22];
ybs[4560]=['',3.1312229,-0.8232704,6.26];
ybs[4561]=['',3.1326731,-0.5831616,6.21];
ybs[4562]=['',3.1334924,0.7024284,6.62];
ybs[4563]=['',3.1353054,-1.0916391,5.57];
ybs[4564]=['',3.1373098,0.5615852,6.42];
ybs[4565]=['',3.138294,1.0710614,6.76];
ybs[4566]=['',3.1378704,-0.9846214,5.44];
ybs[4567]=['',3.138248,-0.7163642,6.79];
ybs[4568]=['',3.1402366,-1.1246307,5.61];
ybs[4569]=['',3.1407333,-0.4538958,6.43];
ybs[4570]=['',3.1413901,0.0075595,6.17];
ybs[4571]=['',3.142424,0.5771817,5.96];
ybs[4572]=['',3.141932,-0.9039774,6.05];
ybs[4573]=['ε Cha',3.1438579,-1.3669309,4.91];
ybs[4574]=['',3.1453005,0.5923224,6.5];
ybs[4575]=['7 Vir',3.1452805,0.0620962,5.37];
ybs[4576]=['',3.1468242,1.4094516,6.17];
ybs[4577]=['',3.1487439,-0.1840194,5.55];
ybs[4578]=['',3.1486003,-0.3828318,6.28];
ybs[4579]=['π Vir',3.1493157,0.1137386,4.66];
ybs[4580]=['',3.149233,-0.3448127,5.26];
ybs[4581]=['',3.1500008,-0.0325588,6.31];
ybs[4582]=['',3.152001,-1.0053277,6.16];
ybs[4583]=['',3.1527307,0.6273502,5.59];
ybs[4584]=['67 UMa',3.1547096,0.7495864,5.21];
ybs[4585]=['',3.155971,-1.4919082,6.05];
ybs[4586]=['',3.1563779,-1.2494168,6.42];
ybs[4587]=['',3.1570348,-1.2093324,5.89];
ybs[4588]=['',3.1579886,-0.1358046,6.22];
ybs[4589]=['θ1 Cru',3.1587567,-1.1067166,4.33];
ybs[4590]=['',3.1615038,-0.7423161,5.15];
ybs[4591]=['',3.1619283,-1.2969769,6.44];
ybs[4592]=['2 Com',3.1641529,0.372833,5.87];
ybs[4593]=['θ2 Cru',3.1644259,-1.104147,4.72];
ybs[4594]=['',3.1658693,-1.1942689,5.35];
ybs[4595]=['κ Cha',3.1665039,-1.3372114,5.04];
ybs[4596]=['',3.1645404,1.4876556,6.27];
ybs[4597]=['',3.1671843,-1.0658127,5.96];
ybs[4598]=['ο Vir',3.1682254,0.1507207,4.12];
ybs[4599]=['',3.1682405,1.3405601,5.8];
ybs[4600]=['',3.1701214,1.0966892,6.13];
ybs[4601]=['',3.1712954,-1.1457146,6.33];
ybs[4602]=['',3.1714803,-0.6246757,6.23];
ybs[4603]=['',3.1716742,-0.0563577,6.37];
ybs[4604]=['',3.1732553,-1.1998828,6.23];
ybs[4605]=['',3.1734788,-1.148541,6.06];
ybs[4606]=['η Cru',3.1756482,-1.1294199,4.15];
ybs[4607]=['',3.1799078,-1.3171006,5.18];
ybs[4608]=['',3.1808684,-0.8859073,4.47];
ybs[4609]=['',3.1808395,-0.8876866,6.37];
ybs[4610]=['',3.1815554,-0.8515485,5.34];
ybs[4611]=['δ Cen',3.1820564,-0.8869738,2.6];
ybs[4612]=['',3.1823199,-1.0636836,6.22];
ybs[4613]=['α Crv',3.1822466,-0.4332997,4.02];
ybs[4614]=['',3.1843931,-0.7753356,5.75];
ybs[4615]=['',3.1844371,-0.7213225,5.48];
ybs[4616]=['10 Vir',3.1877771,0.0314237,5.95];
ybs[4617]=['',3.1879408,1.3013883,6.35];
ybs[4618]=['',3.1893752,-0.6074152,6.17];
ybs[4619]=['11 Vir',3.1893788,0.0996516,5.72];
ybs[4620]=['ε Crv',3.1897176,-0.3964872,3];
ybs[4621]=['',3.1916611,-0.6626595,6.06];
ybs[4622]=['3 Com',3.1914129,0.2916768,6.39];
ybs[4623]=['',3.1924491,0.4744517,6.01];
ybs[4624]=['',3.1940156,-1.0711924,6.08];
ybs[4625]=['3 Crv',3.1938237,-0.4136396,5.46];
ybs[4626]=['',3.1938004,-0.7944753,6.61];
ybs[4627]=['',3.1958976,-0.8980895,6.23];
ybs[4628]=['ρ Cen',3.1964628,-0.9157028,3.96];
ybs[4629]=['',3.1929268,1.4244102,6];
ybs[4630]=['4 Com',3.197182,0.4498235,5.66];
ybs[4631]=['68 UMa',3.1966288,0.9940899,6.43];
ybs[4632]=['',3.1979032,0.4963512,6.49];
ybs[4633]=['5 Com',3.1985074,0.3568267,5.57];
ybs[4634]=['',3.1996527,-1.1003971,5.92];
ybs[4635]=['',3.2015425,-1.22608,6.17];
ybs[4636]=['',3.1983076,1.3529636,5.14];
ybs[4637]=['',3.2032435,-0.5973007,6.5];
ybs[4638]=['',3.2041525,-0.6811394,5.76];
ybs[4639]=['',3.206812,-1.3730653,6.35];
ybs[4640]=['12 Vir',3.2040944,0.1774123,5.85];
ybs[4641]=['',3.2049683,-0.5914924,6.33];
ybs[4642]=['',3.2068923,-0.7997294,5.31];
ybs[4643]=['',3.2080433,-1.1258392,6.22];
ybs[4644]=['',3.2095976,0.9309151,6.16];
ybs[4645]=['',3.2109706,-0.3654959,5.83];
ybs[4646]=['δ Cru',3.2117798,-1.027058,2.8];
ybs[4647]=['',3.2117481,-0.1816836,6.11];
ybs[4648]=['',3.2132805,-0.7332171,6.26];
ybs[4649]=['',3.2112547,1.2235246,5.71];
ybs[4650]=['δ UMa',3.2126368,0.9937085,3.31];
ybs[4651]=['',3.214426,-0.4092936,6.54];
ybs[4652]=['γ Crv',3.2145139,-0.3078609,2.59];
ybs[4653]=['6 Com',3.2153016,0.2583386,5.1];
ybs[4654]=['',3.2174293,-1.2690618,6.22];
ybs[4655]=['',3.2135882,1.2645547,6.29];
ybs[4656]=['2 CVn',3.2157667,0.7079597,5.66];
ybs[4657]=['7 Com',3.216755,0.416228,4.95];
ybs[4658]=['',3.2174267,0.5753343,5];
ybs[4659]=['',3.2204082,-1.1482507,6.06];
ybs[4660]=['',3.2199616,-0.293054,6.05];
ybs[4661]=['ε Mus',3.2224939,-1.1878355,4.11];
ybs[4662]=['',3.2216507,0.9266648,5.81];
ybs[4663]=['',3.2218273,0.5033546,5.7];
ybs[4664]=['β Cha',3.2262961,-1.3859541,4.26];
ybs[4665]=['',3.2232212,-0.6316523,6.15];
ybs[4666]=['',3.2228679,0.2626205,6.34];
ybs[4667]=['',3.2247174,-0.0707129,6.99];
ybs[4668]=['',3.2247537,-0.0706111,6.54];
ybs[4669]=['ζ Cru',3.2262321,-1.1187588,4.04];
ybs[4670]=['',3.2262622,0.526253,6.23];
ybs[4671]=['13 Vir',3.2269788,-0.0154341,5.9];
ybs[4672]=['',3.228595,-0.9641222,5];
ybs[4673]=['',3.2170523,1.5039848,6.33];
ybs[4674]=['',3.2284917,0.4522271,6.48];
ybs[4675]=['8 Com',3.2297362,0.4003376,6.27];
ybs[4676]=['',3.2096365,1.5277255,6.28];
ybs[4677]=['',3.2271391,1.3101047,5.38];
ybs[4678]=['9 Com',3.230483,0.4897374,6.33];
ybs[4679]=['η Vir',3.2323673,-0.0133341,3.89];
ybs[4680]=['3 CVn',3.2317788,0.8532412,5.29];
ybs[4681]=['',3.2336132,-0.3887301,5.97];
ybs[4682]=['',3.2351558,-1.1508666,6.21];
ybs[4683]=['',3.2341283,0.4629035,5.54];
ybs[4684]=['',3.233985,0.452126,6.15];
ybs[4685]=['16 Vir',3.234291,0.0561206,4.96];
ybs[4686]=['ζ Crv',3.2352871,-0.3894327,5.21];
ybs[4687]=['11 Com',3.2358495,0.3088494,4.74];
ybs[4688]=['',3.2356962,0.4705008,7.13];
ybs[4689]=['',3.2368616,-0.2384566,5.14];
ybs[4690]=['ε Cru',3.2389921,-1.055891,3.59];
ybs[4691]=['70 UMa',3.2362166,1.0082223,5.55];
ybs[4692]=['',3.2415679,-0.9856168,5.92];
ybs[4693]=['ζ2 Mus',3.242436,-1.1801724,5.15];
ybs[4694]=['ζ1 Mus',3.2427884,-1.1938829,5.74];
ybs[4695]=['',3.2422076,0.4306937,6.19];
ybs[4696]=['',3.2453578,-1.0083297,5.39];
ybs[4697]=['12 Com',3.2436206,0.4494078,4.81];
ybs[4698]=['17 Vir',3.2438107,0.0909075,6.4];
ybs[4699]=['',3.2602592,-1.5019039,6.33];
ybs[4700]=['',3.2472942,-1.1820866,6.36];
ybs[4701]=['6 Crv',3.2475153,-0.4352408,5.68];
ybs[4702]=['',3.2485638,-0.6197607,5.32];
ybs[4703]=['',3.2486922,-0.6876588,6.4];
ybs[4704]=['',3.2492728,-0.680818,5.79];
ybs[4705]=['4 CVn',3.2491162,0.7408205,6.06];
ybs[4706]=['5 CVn',3.2501073,0.8982397,4.8];
ybs[4707]=['13 Com',3.2514817,0.4538162,5.18];
ybs[4708]=['',3.2536414,-0.7239802,6.25];
ybs[4709]=['',3.2520796,0.4448133,6.42];
ybs[4710]=['',3.2562725,-1.1496024,6.3];
ybs[4711]=['',3.2553815,-0.7437068,6.11];
ybs[4712]=['',3.2554793,-0.2043274,5.95];
ybs[4713]=['',3.2560279,-0.486004,6.09];
ybs[4714]=['',3.2563097,-0.615808,5.73];
ybs[4715]=['',3.2556031,0.4158996,6.03];
ybs[4716]=['71 UMa',3.2545345,0.9892643,5.81];
ybs[4717]=['',3.2546716,1.1118786,6.32];
ybs[4718]=['6 CVn',3.2581334,0.6793139,5.02];
ybs[4719]=['',3.261595,-1.1033841,4.86];
ybs[4720]=['α1 Cru',3.2619594,-1.1029768,1.33];
ybs[4721]=['α2 Cru',3.2620032,-1.1029817,1.73];
ybs[4722]=['',3.2615019,-0.8996751,4.82];
ybs[4723]=['14 Com',3.260603,0.4742334,4.95];
ybs[4724]=['',3.2626898,-0.8553871,6.26];
ybs[4725]=['',3.262843,-0.57468,5.55];
ybs[4726]=['',3.2655253,-1.1150189,6];
ybs[4727]=['γ Com',3.2629381,0.4916872,4.36];
ybs[4728]=['16 Com',3.2631626,0.466506,5];
ybs[4729]=['',3.265748,-1.0312915,5.5];
ybs[4730]=['',3.2601065,1.2537217,6.24];
ybs[4731]=['',3.2663461,0.1485901,6.37];
ybs[4732]=['',3.266971,-0.2919697,6.35];
ybs[4733]=['σ Cen',3.2681077,-0.8783759,3.91];
ybs[4734]=['',3.2694984,-1.1246561,6.04];
ybs[4735]=['73 UMa',3.265565,0.9706836,5.7];
ybs[4736]=['',3.2670854,-0.0822393,6.22];
ybs[4737]=['',3.2699334,-1.080218,6.22];
ybs[4738]=['',3.2694888,-0.6830878,5.44];
ybs[4739]=['',3.2704323,-0.9861882,6.15];
ybs[4740]=['',3.2703442,0.4560549,6.54];
ybs[4741]=['',3.270818,0.450339,6.65];
ybs[4742]=['17 Com',3.2715518,0.4505768,5.29];
ybs[4743]=['18 Com',3.2739,0.4190935,5.48];
ybs[4744]=['',3.2763028,-0.988228,5.8];
ybs[4745]=['',3.2764543,-0.7301181,6.02];
ybs[4746]=['20 Com',3.2750991,0.3630202,5.69];
ybs[4747]=['δ Crv',3.2758835,-0.2899364,2.95];
ybs[4748]=['',3.2768089,-0.2354383,6.35];
ybs[4749]=['',3.2777786,-0.4152701,5.63];
ybs[4750]=['74 UMa',3.275845,1.0176886,5.35];
ybs[4751]=['7 CVn',3.2763335,0.8977797,6.21];
ybs[4752]=['75 UMa',3.2763475,1.024001,6.08];
ybs[4753]=['γ Cru',3.2818519,-0.9984999,1.63];
ybs[4754]=['γ Cru',3.2823472,-0.9979374,6.42];
ybs[4755]=['4 Dra',3.2763008,1.2061018,4.95];
ybs[4756]=['21 Com',3.2807067,0.4270944,5.46];
ybs[4757]=['',3.2797422,0.9246779,6.21];
ybs[4758]=['',3.2841033,-1.0388262,5.48];
ybs[4759]=['',3.2866532,-1.2758026,5.88];
ybs[4760]=['',3.2822942,0.1310336,6.05];
ybs[4761]=['',3.2853093,-1.1100741,5.95];
ybs[4762]=['',3.2836047,-0.0898666,6.19];
ybs[4763]=['γ Mus',3.2879279,-1.2606421,3.87];
ybs[4764]=['',3.2856135,-0.5695019,6.46];
ybs[4765]=['η Crv',3.285508,-0.2843588,4.31];
ybs[4766]=['',3.287811,-0.2435709,5.74];
ybs[4767]=['20 Vir',3.2896619,0.178009,6.26];
ybs[4768]=['',3.2912151,-0.3471165,6.26];
ybs[4769]=['',3.2920477,-0.2256123,5.58];
ybs[4770]=['22 Com',3.2918698,0.4221375,6.29];
ybs[4771]=['21 Vir',3.2929416,-0.1666491,5.48];
ybs[4772]=['',3.2941042,-0.8727653,6.38];
ybs[4773]=['',3.2921599,0.5785966,5.42];
ybs[4774]=['',3.2927764,0.5809918,6.24];
ybs[4775]=['β CVn',3.29251,0.7201429,4.26];
ybs[4776]=['β Crv',3.295658,-0.4100297,2.65];
ybs[4777]=['κ Dra',3.2909268,1.2163543,3.87];
ybs[4778]=['',3.2971928,-0.7813772,5.77];
ybs[4779]=['23 Com',3.2974668,0.3932732,4.81];
ybs[4780]=['',3.3008292,-1.081025,6.22];
ybs[4781]=['24 Com',3.2985936,0.319063,6.56];
ybs[4782]=['24 Com',3.2987026,0.3190582,5.02];
ybs[4783]=['',3.2987062,0.3802223,5.85];
ybs[4784]=['',3.3017646,-0.7176472,5.13];
ybs[4785]=['6 Dra',3.2963472,1.220433,4.94];
ybs[4786]=['',3.3028987,-0.6975417,5.8];
ybs[4787]=['',3.3025851,-0.3599466,6.2];
ybs[4788]=['α Mus',3.3084818,-1.2083206,2.69];
ybs[4789]=['25 Vir',3.3060611,-0.1034646,5.87];
ybs[4790]=['',3.30382,1.0365644,5.5];
ybs[4791]=['25 Com',3.3067408,0.2965893,5.68];
ybs[4792]=['τ Cen',3.3103378,-0.8488791,3.86];
ybs[4793]=['',3.3101628,-0.4753399,5.45];
ybs[4794]=['',3.3178601,-1.3171198,6.49];
ybs[4795]=['',3.3116163,0.0556139,6.33];
ybs[4796]=['',3.3158238,-1.1744155,6.25];
ybs[4797]=['',3.3129324,0.0306949,5.71];
ybs[4798]=['',3.3134591,0.1202935,7.08];
ybs[4799]=['',3.3146515,-0.320203,6];
ybs[4800]=['',3.316102,-0.5326431,5.89];
ybs[4801]=['9 CVn',3.3144295,0.711718,6.37];
ybs[4802]=['',3.3157071,0.3938066,6.38];
ybs[4803]=['χ Vir',3.3167934,-0.1412238,4.66];
ybs[4804]=['',3.3204078,-1.1625216,6.26];
ybs[4805]=['26 Com',3.3160944,0.3659348,5.46];
ybs[4806]=['',3.3166926,0.6258048,6.45];
ybs[4807]=['',3.319744,-0.6995877,4.64];
ybs[4808]=['',3.3263971,-0.807064,5.84];
ybs[4809]=['γ Cen',3.3270134,-0.8561804,2.17];
ybs[4810]=['',3.3299757,-1.2130605,6.33];
ybs[4811]=['',3.3256374,-0.2288027,6.08];
ybs[4812]=['',3.3256519,-0.2288269,5.98];
ybs[4813]=['',3.3290552,-1.0433857,4.93];
ybs[4814]=['27 Vir',3.3268472,0.1803029,6.19];
ybs[4815]=['γ Vir',3.3272903,-0.0269694,3.65];
ybs[4816]=['γ Vir',3.3272903,-0.0269694,3.68];
ybs[4817]=['',3.3280927,-0.3465244,6.03];
ybs[4818]=['ρ Vir',3.3282078,0.1769726,4.88];
ybs[4819]=['31 Vir',3.3285177,0.1171273,5.59];
ybs[4820]=['',3.3330541,-1.1022505,5.31];
ybs[4821]=['',3.3317039,-0.853619,4.66];
ybs[4822]=['',3.3328464,-0.9781333,6.08];
ybs[4823]=['76 UMa',3.326269,1.0928773,6.07];
ybs[4824]=['',3.3342722,-0.9821278,6];
ybs[4825]=['',3.3357185,-1.0297215,6.4];
ybs[4826]=['',3.3353155,-0.7029038,6.44];
ybs[4827]=['',3.335909,-0.029192,5.93];
ybs[4828]=['',3.3376394,-0.6360812,6.39];
ybs[4829]=['',3.3377058,-0.4960136,5.48];
ybs[4830]=['',3.3328538,1.0656958,6.38];
ybs[4831]=['',3.3428489,-1.2029964,6.16];
ybs[4832]=['ι Cru',3.3452262,-1.0659873,4.69];
ybs[4833]=['',3.3391431,0.7680755,6.33];
ybs[4834]=['β Mus',3.3483087,-1.1903749,3.05];
ybs[4835]=['10 CVn',3.341546,0.6838787,5.95];
ybs[4836]=['',3.3420826,0.7914155,4.99];
ybs[4837]=['32 Vir',3.3445083,0.1322587,5.22];
ybs[4838]=['',3.3483925,-0.9875821,4.65];
ybs[4839]=['33 Vir',3.3477976,0.1648393,5.67];
ybs[4840]=['',3.3497934,-0.5831305,5.86];
ybs[4841]=['27 Com',3.3489305,0.2876673,5.12];
ybs[4842]=['',3.3374002,1.4054355,6.4];
ybs[4843]=['β Cru',3.3543328,-1.0434257,1.25];
ybs[4844]=['',3.3507134,0.1021975,6.34];
ybs[4845]=['34 Vir',3.3514974,0.2070436,6.07];
ybs[4846]=['',3.3530507,-0.111653,6.26];
ybs[4847]=['',3.3546528,-0.4354061,6.44];
ybs[4848]=['35 Vir',3.3542984,0.0606939,6.41];
ybs[4849]=['',3.3512753,1.0940684,5.89];
ybs[4850]=['',3.3570541,-0.4833291,5.66];
ybs[4851]=['28 Com',3.3558995,0.2348832,6.56];
ybs[4852]=['',3.3637165,-1.2580591,5.55];
ybs[4853]=['7 Dra',3.3522626,1.1640469,5.43];
ybs[4854]=['',3.3581962,0.4318832,6.31];
ybs[4855]=['29 Com',3.3587961,0.2448229,5.7];
ybs[4856]=['11 CVn',3.357578,0.8442462,6.27];
ybs[4857]=['',3.3571787,1.0511209,5.85];
ybs[4858]=['',3.3652034,-1.0558524,6.75];
ybs[4859]=['30 Com',3.3603838,0.4792164,5.78];
ybs[4860]=['ι Oct',3.3902546,-1.4819482,5.46];
ybs[4861]=['',3.3655206,-0.8474405,6.24];
ybs[4862]=['',3.3683824,-0.9229734,5.73];
ybs[4863]=['',3.3647807,0.3973814,6.43];
ybs[4864]=['',3.3669229,-0.5950604,4.91];
ybs[4865]=['',3.364165,0.653135,5.89];
ybs[4866]=['',3.3700006,-1.0546095,5.72];
ybs[4867]=['',3.3697709,-0.1820951,6.41];
ybs[4868]=['37 Vir',3.370699,0.0516921,6.02];
ybs[4869]=['',3.3724867,-0.6942173,5.98];
ybs[4870]=['',3.3732161,-0.8410575,6.33];
ybs[4871]=['',3.3724409,-0.4683232,6.15];
ybs[4872]=['',3.3747211,-0.9411563,6.24];
ybs[4873]=['31 Com',3.3708832,0.4790167,4.94];
ybs[4874]=['32 Com',3.3731749,0.2963398,6.32];
ybs[4875]=['',3.3776085,-0.9607563,5.93];
ybs[4876]=['',3.374294,0.2797353,6.3];
ybs[4877]=['',3.3790443,-1.0545865,5.76];
ybs[4878]=['',3.3777174,-0.8558764,4.33];
ybs[4879]=['',3.3790029,-0.7029075,4.27];
ybs[4880]=['κ Cru',3.3810368,-1.0554294,5.9];
ybs[4881]=['38 Vir',3.3776016,-0.0636666,6.11];
ybs[4882]=['',3.3566722,1.4542583,5.85];
ybs[4883]=['',3.3571765,1.4541664,5.28];
ybs[4884]=['35 Com',3.3779035,0.3691413,4.9];
ybs[4885]=['',3.3833832,-1.0214575,6.58];
ybs[4886]=['',3.3795627,-0.0753741,6.44];
ybs[4887]=['λ Cru',3.3846563,-1.0339555,4.62];
ybs[4888]=['μ1 Cru',3.3843435,-0.999592,4.03];
ybs[4889]=['μ2 Cru',3.3844308,-0.9994271,5.17];
ybs[4890]=['41 Vir',3.3802891,0.2150927,6.25];
ybs[4891]=['',3.3825688,-0.2049587,6];
ybs[4892]=['ψ Vir',3.3827353,-0.168137,4.79];
ybs[4893]=['',3.3857761,-0.7722477,5.89];
ybs[4894]=['',3.3818127,0.5836342,6.26];
ybs[4895]=['ε UMa',3.3806583,0.9750287,1.77];
ybs[4896]=['',3.3872816,-0.7506729,5.47];
ybs[4897]=['',3.393452,-1.2615187,5.93];
ybs[4898]=['',3.39027,-0.9936264,5.32];
ybs[4899]=['',3.3847961,0.8220861,5.84];
ybs[4900]=['δ Vir',3.3880974,0.0576478,3.38];
ybs[4901]=['',3.3894788,-0.2691549,6.17];
ybs[4902]=['',3.3922499,-0.4634672,6.62];
ybs[4903]=['',3.3950622,-0.8952314,5.16];
ybs[4904]=['α1 CVn',3.3895544,0.6670689,5.6];
ybs[4905]=['α2 CVn',3.3896489,0.667132,2.9];
ybs[4906]=['8 Dra',3.3866662,1.1404689,5.24];
ybs[4907]=['',3.3905582,0.9425647,5.82];
ybs[4908]=['',3.396809,-0.3987765,6.31];
ybs[4909]=['',3.3943362,0.8042925,6.12];
ybs[4910]=['36 Com',3.4024715,0.3022084,4.78];
ybs[4911]=['44 Vir',3.4058444,-0.0681733,5.79];
ybs[4912]=['',3.4099731,-0.5864179,6.02];
ybs[4913]=['δ Mus',3.4185949,-1.2504005,3.62];
ybs[4914]=['37 Com',3.4082408,0.5356584,4.9];
ybs[4915]=['46 Vir',3.409943,-0.0604339,5.99];
ybs[4916]=['',3.4099789,0.3190299,6.2];
ybs[4917]=['',3.4003523,1.3155993,6.01];
ybs[4918]=['9 Dra',3.4059305,1.1606988,5.32];
ybs[4919]=['38 Com',3.4122283,0.2972143,5.96];
ybs[4920]=['',3.4221699,-1.2491287,6.03];
ybs[4921]=['78 UMa',3.4098186,0.9821387,4.93];
ybs[4922]=['ε Vir',3.4167146,0.1896361,2.83];
ybs[4923]=['ξ1 Cen',3.4233607,-0.8660475,4.85];
ybs[4924]=['',3.4141566,1.1085704,6];
ybs[4925]=['',3.4239251,-0.3608762,5.58];
ybs[4926]=['',3.4181809,1.0406063,6.53];
ybs[4927]=['48 Vir',3.4243822,-0.065571,6.59];
ybs[4928]=['',3.4286759,-0.7206494,6.26];
ybs[4929]=['',3.4319841,-0.9112086,6.43];
ybs[4930]=['',3.4352523,-0.8474783,4.71];
ybs[4931]=['',3.4364742,-0.7274863,5.59];
ybs[4932]=['ξ2 Cen',3.4380428,-0.8726534,4.27];
ybs[4933]=['14 CVn',3.4320098,0.6231784,5.25];
ybs[4934]=['',3.440466,-1.04639,5.99];
ybs[4935]=['',3.4324369,0.7884564,5.63];
ybs[4936]=['39 Com',3.4348386,0.3675665,5.99];
ybs[4937]=['',3.4377912,-0.6275364,6.54];
ybs[4938]=['',3.4339574,0.5050302,6.54];
ybs[4939]=['40 Com',3.4349259,0.3930969,5.6];
ybs[4940]=['',3.4268165,1.2728992,6.31];
ybs[4941]=['',3.4413077,-0.934674,5.71];
ybs[4942]=['θ Mus',3.4438154,-1.1414361,5.51];
ybs[4943]=['',3.4342185,1.0812073,6.14];
ybs[4944]=['41 Com',3.4383687,0.4805154,4.8];
ybs[4945]=['49 Vir',3.4418564,-0.1890786,5.19];
ybs[4946]=['',3.4414866,0.4793145,6.19];
ybs[4947]=['',3.4446558,-0.1584321,5.55];
ybs[4948]=['ψ Hya',3.4470328,-0.405109,4.95];
ybs[4949]=['',3.4477086,-0.1680977,6.32];
ybs[4950]=['',3.4473877,0.1732983,5.78];
ybs[4951]=['50 Vir',3.4499637,-0.1819041,5.94];
ybs[4952]=['',3.4498967,0.2924425,5.91];
ybs[4953]=['θ Vir',3.4507709,-0.0982926,4.38];
ybs[4954]=['',3.4489989,0.6515339,6.02];
ybs[4955]=['',3.4558687,-0.9190847,6.06];
ybs[4956]=['',3.4605066,-1.2223333,5.91];
ybs[4957]=['15 CVn',3.4492227,0.6709217,6.28];
ybs[4958]=['α Com',3.4507259,0.3043257,5.22];
ybs[4959]=['α Com',3.4507259,0.3043257,5.22];
ybs[4960]=['',3.4564176,-0.7387239,5.79];
ybs[4961]=['17 CVn',3.4507629,0.6703117,5.91];
ybs[4962]=['',3.460217,-1.1064581,6.33];
ybs[4963]=['',3.4574857,-0.7585474,5.25];
ybs[4964]=['',3.4492547,1.0844824,6.54];
ybs[4965]=['',3.46185,-1.0474311,4.6];
ybs[4966]=['',3.4723294,-1.3707724,5.85];
ybs[4967]=['',3.4644395,-1.1574922,5.9];
ybs[4968]=['',3.4584164,-0.4650309,6.5];
ybs[4969]=['',3.4603094,-0.6614038,4.85];
ybs[4970]=['',3.4646611,-1.0456116,6.16];
ybs[4971]=['53 Vir',3.4600741,-0.2843351,5.04];
ybs[4972]=['',3.4638615,-0.7468648,6.22];
ybs[4973]=['β Com',3.4588342,0.4849474,4.26];
ybs[4974]=['',3.4600382,0.421767,6.33];
ybs[4975]=['',3.4663973,-0.8864947,5.89];
ybs[4976]=['',3.4619496,0.2000773,5.77];
ybs[4977]=['',3.4620957,0.3256635,6.53];
ybs[4978]=['',3.4701692,-1.0258379,5.89];
ybs[4979]=['',3.4703811,-1.0331585,4.92];
ybs[4980]=['54 Vir',3.466162,-0.3302001,6.28];
ybs[4981]=['',3.4687198,-0.7545271,6.16];
ybs[4982]=['',3.464749,0.3252334,6.11];
ybs[4983]=['η Mus',3.4752025,-1.1865898,4.8];
ybs[4984]=['',3.4761239,-1.2177483,6.37];
ybs[4985]=['55 Vir',3.4693825,-0.3494697,5.33];
ybs[4986]=['',3.4721572,-0.8560646,5.89];
ybs[4987]=['',3.4667007,0.6991859,4.92];
ybs[4988]=['',3.4705593,0.1961646,5.67];
ybs[4989]=['',3.4738674,-0.6364042,6.19];
ybs[4990]=['',3.4815942,-1.1384829,6.07];
ybs[4991]=['57 Vir',3.4772339,-0.3496787,5.22];
ybs[4992]=['',3.4837507,-1.1671973,4.87];
ybs[4993]=['',3.4646228,1.2689672,6.59];
ybs[4994]=['19 CVn',3.4746049,0.7114511,5.79];
ybs[4995]=['',3.4789871,-0.0258755,6.68];
ybs[4996]=['',3.4813262,-0.5514899,5.1];
ybs[4997]=['',3.4779663,0.3309081,6.45];
ybs[4998]=['',3.4830398,-0.7691896,5.84];
ybs[4999]=['',3.4583761,1.4028747,6.25];
ybs[5000]=['',3.4792665,0.3437127,6.45];
ybs[5001]=['59 Vir',3.4804075,0.1628779,5.22];
ybs[5002]=['',3.4933811,-1.2588556,6.04];
ybs[5003]=['',3.4824794,0.2370798,5.33];
ybs[5004]=['',3.4836636,-0.0134131,6.37];
ybs[5005]=['σ Vir',3.4840679,0.0938619,4.8];
ybs[5006]=['',3.4890805,-0.8967116,6.19];
ybs[5007]=['20 CVn',3.4833607,0.7065206,4.73];
ybs[5008]=['',3.4777896,1.1923399,6.2];
ybs[5009]=['61 Vir',3.4878029,-0.3211947,4.74];
ybs[5010]=['γ Hya',3.4901146,-0.4060213,3];
ybs[5011]=['',3.4895248,0.0627642,6.62];
ybs[5012]=['',3.4874793,0.5935226,5.82];
ybs[5013]=['21 CVn',3.486227,0.8655121,5.15];
ybs[5014]=['',3.4981238,-1.0448365,6.18];
ybs[5015]=['',3.4901082,0.611501,6.02];
ybs[5016]=['',3.4980977,-0.9222221,5.48];
ybs[5017]=['',3.4989557,-0.9754978,6.02];
ybs[5018]=['ι Cen',3.4976156,-0.6423443,2.75];
ybs[5019]=['',3.4994012,-0.8198141,5.77];
ybs[5020]=['',3.509011,-1.2607856,6.05];
ybs[5021]=['',3.4975674,0.0497468,6.26];
ybs[5022]=['23 CVn',3.4954597,0.6991633,5.6];
ybs[5023]=['',3.5013209,-0.3417381,6.21];
ybs[5024]=['',3.5069846,-1.0657557,6.18];
ybs[5025]=['',3.5071457,-1.0660368,4.53];
ybs[5026]=['',3.5052525,-0.9123567,5.83];
ybs[5027]=['',3.5019393,0.0348365,5.69];
ybs[5028]=['',3.507803,-0.8383533,6.16];
ybs[5029]=['',3.5085391,-0.849169,6.38];
ybs[5030]=['64 Vir',3.5039505,0.0883757,5.87];
ybs[5031]=['',3.5133528,-1.1279487,4.53];
ybs[5032]=['ι1 Mus',3.5192156,-1.3086209,5.05];
ybs[5033]=['',3.5086924,-0.5808633,6.22];
ybs[5034]=['63 Vir',3.5079363,-0.3111279,5.37];
ybs[5035]=['',3.5029998,0.7646613,6.35];
ybs[5036]=['',3.5121915,-0.8711628,6.48];
ybs[5037]=['65 Vir',3.5090893,-0.0875359,5.89];
ybs[5038]=['',3.5187054,-1.127063,5.31];
ybs[5039]=['',3.5218185,-1.2342632,5.67];
ybs[5040]=['66 Vir',3.5144959,-0.0917117,5.75];
ybs[5041]=['ι2 Mus',3.528781,-1.305197,6.63];
ybs[5042]=['',3.5111285,0.6447766,6.07];
ybs[5043]=['',3.514108,0.2153934,6.44];
ybs[5044]=['ζ UMa',3.5108289,0.9570401,2.27];
ybs[5045]=['ζ UMa',3.5108944,0.9569772,3.95];
ybs[5046]=['α Vir',3.5173552,-0.1963861,0.98];
ybs[5047]=['',3.5166024,0.4147552,5.78];
ybs[5048]=['',3.5218421,-0.6954409,5.09];
ybs[5049]=['',3.5215946,-0.0223935,5.97];
ybs[5050]=['',3.5253927,-0.7258559,5.69];
ybs[5051]=['',3.5263082,-0.8593003,6.31];
ybs[5052]=['80 UMa',3.5164865,0.9581393,4.01];
ybs[5053]=['',3.5273708,-0.8634351,6.28];
ybs[5054]=['68 Vir',3.5240364,-0.2233715,5.25];
ybs[5055]=['',3.5267213,-0.7025549,6.4];
ybs[5056]=['',3.5345988,-1.2168114,6.2];
ybs[5057]=['',3.5213095,0.8017608,5.88];
ybs[5058]=['69 Vir',3.527276,-0.2803689,4.76];
ybs[5059]=['',3.5357687,-1.130378,6.11];
ybs[5060]=['',3.5195059,1.1025334,6.5];
ybs[5061]=['',3.5364838,-0.8945737,5.06];
ybs[5062]=['70 Vir',3.5311977,0.2389128,4.98];
ybs[5063]=['',3.519375,1.261887,5.79];
ybs[5064]=['',3.524112,1.1282702,6.66];
ybs[5065]=['',3.5245551,1.1279893,7.04];
ybs[5066]=['',3.5286029,0.9190128,6.34];
ybs[5067]=['',3.5308314,0.7092935,6.47];
ybs[5068]=['',3.5349415,-0.0253859,6.43];
ybs[5069]=['',3.5295509,0.8813385,6.8];
ybs[5070]=['',3.5372356,-0.4079074,4.97];
ybs[5071]=['71 Vir',3.5346639,0.1872435,5.65];
ybs[5072]=['',3.5554603,-1.3553821,6.48];
ybs[5073]=['',3.5320151,0.8836284,6.43];
ybs[5074]=['κ Oct',3.595964,-1.4947315,5.58];
ybs[5075]=['',3.5303413,1.0446777,5.4];
ybs[5076]=['',3.5381307,0.1237255,6.17];
ybs[5077]=['',3.5379624,0.1033826,6.51];
ybs[5078]=['72 Vir',3.5401492,-0.1144962,6.09];
ybs[5079]=['',3.5433092,-0.6893571,3.88];
ybs[5080]=['',3.545346,-0.4922256,6.47];
ybs[5081]=['',3.5217495,1.3710155,5.77];
ybs[5082]=['',3.547848,-0.6717554,6.16];
ybs[5083]=['',3.5553948,-1.147062,6.37];
ybs[5084]=['73 Vir',3.5473626,-0.3284445,6.01];
ybs[5085]=['74 Vir',3.5468529,-0.1107489,4.69];
ybs[5086]=['',3.5431203,0.733324,6.08];
ybs[5087]=['',3.549911,-0.5023455,5.69];
ybs[5088]=['',3.5498225,-0.5175736,6.45];
ybs[5089]=['75 Vir',3.5508772,-0.2696973,5.55];
ybs[5090]=['76 Vir',3.5512784,-0.1789738,5.21];
ybs[5091]=['',3.5514318,-0.1271374,6.68];
ybs[5092]=['',3.550126,0.4233629,6.11];
ybs[5093]=['',3.5585295,-0.8440656,6.33];
ybs[5094]=['',3.559281,-0.5829395,6.44];
ybs[5095]=['78 Vir',3.5561899,0.062302,4.94];
ybs[5096]=['',3.5587654,-0.2321916,5.91];
ybs[5097]=['ζ Vir',3.5586915,-0.0119553,3.37];
ybs[5098]=['',3.5566957,0.6754415,6.37];
ybs[5099]=['81 UMa',3.5552055,0.9644574,5.6];
ybs[5100]=['',3.558615,0.6474011,4.98];
ybs[5101]=['80 Vir',3.562364,-0.0957334,5.73];
ybs[5102]=['24 CVn',3.5568582,0.8539355,4.7];
ybs[5103]=['',3.5709448,-1.0782755,5.63];
ybs[5104]=['',3.5623196,0.1765525,6.49];
ybs[5105]=['',3.5811648,-1.3224743,6.34];
ybs[5106]=['',3.5603783,0.7698276,6.84];
ybs[5107]=['',3.5685394,-0.6031255,6.5];
ybs[5108]=['',3.5698678,-0.771995,5.98];
ybs[5109]=['',3.5784528,-1.2310401,6.1];
ybs[5110]=['',3.5682547,-0.4639745,5.78];
ybs[5111]=['',3.5712139,-0.8118749,5.9];
ybs[5112]=['',3.5748212,-1.0210793,6.42];
ybs[5113]=['',3.5683625,0.4280345,5.74];
ybs[5114]=['',3.5778114,-1.007255,6.01];
ybs[5115]=['',3.583963,-1.2370284,6.59];
ybs[5116]=['',3.5664633,0.862155,6.49];
ybs[5117]=['25 CVn',3.5702393,0.6319195,4.82];
ybs[5118]=['',3.5765714,-0.5174775,5.83];
ybs[5119]=['',3.5735062,0.2480656,6.52];
ybs[5120]=['',3.5841574,-1.1286189,5.79];
ybs[5121]=['',3.5559226,1.3344345,6.57];
ybs[5122]=['ε Cen',3.5823196,-0.9347042,2.3];
ybs[5123]=['',3.5710143,0.8835918,6.48];
ybs[5124]=['',3.5826772,-0.8733362,6];
ybs[5125]=['',3.5810419,-0.6952749,6.27];
ybs[5126]=['',3.5816163,-0.7005784,5.6];
ybs[5127]=['',3.577407,0.3172464,6.48];
ybs[5128]=['',3.5798555,0.1860139,5.57];
ybs[5129]=['',3.5674893,1.2418621,5.5];
ybs[5130]=['',3.5918621,-1.0275635,5.38];
ybs[5131]=['',3.5904857,-0.9537854,5.01];
ybs[5132]=['82 UMa',3.5787378,0.9221109,5.46];
ybs[5133]=['',3.582539,0.5397214,6.21];
ybs[5134]=['1 Boo',3.5845207,0.3467524,5.75];
ybs[5135]=['',3.5842957,0.4882937,6.23];
ybs[5136]=['',3.5887671,-0.4108097,6.59];
ybs[5137]=['',3.5900135,-0.5879112,6.05];
ybs[5138]=['',3.5826486,0.8801918,6.32];
ybs[5139]=['2 Boo',3.5860681,0.3910899,5.62];
ybs[5140]=['82 Vir',3.5889907,-0.1534316,5.01];
ybs[5141]=['',3.5957795,-0.9923194,6];
ybs[5142]=['',3.5954464,-0.8879925,6.41];
ybs[5143]=['',3.5822668,0.9969202,6.29];
ybs[5144]=['83 UMa',3.5840374,0.9528373,4.66];
ybs[5145]=['',3.5952197,-0.7241159,5.98];
ybs[5146]=['',3.5913807,0.1448713,6.16];
ybs[5147]=['',3.5985123,-0.7357441,5.98];
ybs[5148]=['',3.6013772,-0.8918715,6.47];
ybs[5149]=['84 Vir',3.5951522,0.0602207,5.36];
ybs[5150]=['',3.5920038,0.7258193,6.3];
ybs[5151]=['',3.5932045,0.6091401,5.98];
ybs[5152]=['',3.5868253,1.1298306,5.85];
ybs[5153]=['',3.5989498,-0.0975009,6.51];
ybs[5154]=['',3.5979039,0.3946668,6.13];
ybs[5155]=['83 Vir',3.6016777,-0.2839049,5.6];
ybs[5156]=['',3.6029755,-0.4465978,6.21];
ybs[5157]=['',3.6067115,-0.4573335,5.81];
ybs[5158]=['1 Cen',3.6071464,-0.5782458,4.23];
ybs[5159]=['',3.5979494,0.9071685,6.02];
ybs[5160]=['85 Vir',3.6064226,-0.2767164,6.19];
ybs[5161]=['',3.6146488,-1.0939174,6.51];
ybs[5162]=['',3.6118396,-0.899189,4.65];
ybs[5163]=['86 Vir',3.6079177,-0.2184066,5.51];
ybs[5164]=['',3.6126842,-0.6342327,5.15];
ybs[5165]=['',3.6153187,-0.8785334,5.91];
ybs[5166]=['',3.6161079,-0.8797835,5.45];
ybs[5167]=['',3.6035016,0.973757,6.5];
ybs[5168]=['',3.6134958,-0.1709731,6.05];
ybs[5169]=['',3.6083362,0.7156117,5.87];
ybs[5170]=['',3.6087938,0.6705001,5.94];
ybs[5171]=['87 Vir',3.6144781,-0.3132312,5.43];
ybs[5172]=['3 Boo',3.6107992,0.4470704,5.95];
ybs[5173]=['',3.6120859,0.1093209,6.33];
ybs[5174]=['',3.5899507,1.3609491,5.91];
ybs[5175]=['τ Boo',3.613282,0.3031601,4.5];
ybs[5176]=['',3.6117564,0.6711812,5.5];
ybs[5177]=['84 UMa',3.6095467,0.9485125,5.7];
ybs[5178]=['',3.6565882,-1.44428,5.95];
ybs[5179]=['',3.6213145,-0.6246607,6.53];
ybs[5180]=['ν Cen',3.6240124,-0.7290972,3.41];
ybs[5181]=['η UMa',3.6138422,0.8591645,1.86];
ybs[5182]=['2 Cen',3.6235894,-0.602789,4.19];
ybs[5183]=['μ Cen',3.6245199,-0.742817,3.04];
ybs[5184]=['',3.6354193,-1.2127824,5.75];
ybs[5185]=['',3.6190887,0.5428614,5.62];
ybs[5186]=['89 Vir',3.6251703,-0.318008,4.97];
ybs[5187]=['',3.6263841,-0.5090722,6.18];
ybs[5188]=['',3.6275447,-0.6979112,6.44];
ybs[5189]=['',3.6202612,0.6886411,7.4];
ybs[5190]=['υ Boo',3.6229669,0.2742145,4.07];
ybs[5191]=['6 Boo',3.6239157,0.3696218,4.91];
ybs[5192]=['',3.628268,-0.3487767,6.53];
ybs[5193]=['',3.5865144,1.4427694,5.98];
ybs[5194]=['',3.6237974,0.6378547,6.38];
ybs[5195]=['',3.6271924,0.0944392,6.01];
ybs[5196]=['',3.6341242,-0.8200451,5.77];
ybs[5197]=['',3.6356068,-0.9232366,5.25];
ybs[5198]=['',3.6330874,-0.6373826,6.35];
ybs[5199]=['',3.6316808,-0.4272023,6.45];
ybs[5200]=['3 Cen',3.6339626,-0.5773619,4.56];
ybs[5201]=['3 Cen',3.633999,-0.5773667,6.06];
ybs[5202]=['',3.6347582,-0.553363,6.12];
ybs[5203]=['',3.6229738,1.0716802,5.96];
ybs[5204]=['',3.6296114,0.6053913,6.65];
ybs[5205]=['',3.629955,0.6035057,5.87];
ybs[5206]=['',3.6261968,1.0202003,6.46];
ybs[5207]=['',3.6427918,-0.9330385,5.89];
ybs[5208]=['',3.6484827,-1.182248,5.71];
ybs[5209]=['',3.6327419,0.5996633,4.74];
ybs[5210]=['',3.6353639,0.2108254,6.04];
ybs[5211]=['4 Cen',3.639976,-0.5587402,4.73];
ybs[5212]=['',3.6415305,-0.6239513,5.54];
ybs[5213]=['',3.6435784,-0.8240323,6.1];
ybs[5214]=['',3.6429353,-0.6178463,6.19];
ybs[5215]=['7 Boo',3.6392387,0.3114903,5.7];
ybs[5216]=['10 Dra',3.6300355,1.1281328,4.65];
ybs[5217]=['',3.6277866,1.1908223,6.4];
ybs[5218]=['',3.6445776,-0.5001273,6.04];
ybs[5219]=['',3.6388698,0.498507,5.9];
ybs[5220]=['',3.6492366,-0.9118712,5.71];
ybs[5221]=['ζ Cen',3.6505395,-0.826824,2.55];
ybs[5222]=['90 Vir',3.6460087,-0.0277236,5.15];
ybs[5223]=['',3.6472742,-0.1421434,6.19];
ybs[5224]=['',3.6542417,-0.9462647,6.14];
ybs[5225]=['η Boo',3.6456423,0.3196113,2.68];
ybs[5226]=['',3.6552283,-0.9562558,6];
ybs[5227]=['',3.6510272,-0.5475126,6.51];
ybs[5228]=['86 UMa',3.641165,0.9362474,5.7];
ybs[5229]=['',3.6539485,-0.8146766,5.83];
ybs[5230]=['',3.6757664,-1.3731214,6.09];
ybs[5231]=['',3.660544,-1.1130209,4.71];
ybs[5232]=['',3.6645311,-1.1499119,6.2];
ybs[5233]=['',3.6507182,0.243844,6.16];
ybs[5234]=['92 Vir',3.653666,0.0168518,5.91];
ybs[5235]=['',3.651898,0.5575874,6.32];
ybs[5236]=['',3.6583858,-0.4033035,6.14];
ybs[5237]=['9 Boo',3.6537072,0.4783412,5.01];
ybs[5238]=['φ Cen',3.6623262,-0.7362752,3.83];
ybs[5239]=['υ1 Cen',3.664186,-0.783446,3.87];
ybs[5240]=['47 Hya',3.6630337,-0.4373238,5.15];
ybs[5241]=['',3.6670268,-0.8805955,5.91];
ybs[5242]=['',3.6719413,-1.0745219,6.49];
ybs[5243]=['',3.6748614,-1.1580724,5.97];
ybs[5244]=['',3.6630599,0.254205,6];
ybs[5245]=['10 Boo',3.6628721,0.3771925,5.76];
ybs[5246]=['',3.6568126,1.0717707,6.37];
ybs[5247]=['48 Hya',3.6695111,-0.4379825,5.77];
ybs[5248]=['',3.6683803,-0.0634259,6.4];
ybs[5249]=['',3.6755802,-0.7034759,6.13];
ybs[5250]=['υ2 Cen',3.6775089,-0.7973973,4.34];
ybs[5251]=['θ Aps',3.6961666,-1.3418037,5.5];
ybs[5252]=['',3.6748176,0.1537762,5.99];
ybs[5253]=['11 Boo',3.6737869,0.4765207,6.23];
ybs[5254]=['τ Vir',3.6762673,0.0254908,4.26];
ybs[5255]=['',3.6799413,-0.4802056,5.48];
ybs[5256]=['',3.6854115,-0.98257,5.92];
ybs[5257]=['β Cen',3.687335,-1.0551644,0.61];
ybs[5258]=['',3.6828581,-0.5544475,6.18];
ybs[5259]=['',3.6849687,-0.7244312,6.11];
ybs[5260]=['',3.6799964,0.1675979,6.2];
ybs[5261]=['',3.6778305,0.7970881,6.27];
ybs[5262]=['',3.6864195,-0.392788,6.3];
ybs[5263]=['',3.6843777,0.1868051,6.3];
ybs[5264]=['',3.6847557,0.1302519,6.26];
ybs[5265]=['',3.6861782,0.0840796,6.24];
ybs[5266]=['',3.6877036,-0.0953779,6.39];
ybs[5267]=['',3.6887582,-0.2627589,6.28];
ybs[5268]=['',3.6955617,-0.9556102,6.17];
ybs[5269]=['',3.7093499,-1.3078208,6.02];
ybs[5270]=['',3.6810996,0.8881681,6.15];
ybs[5271]=['',3.6986454,-1.0436789,6.42];
ybs[5272]=['',3.674941,1.1972029,6.34];
ybs[5273]=['',3.6892505,0.0386454,6.28];
ybs[5274]=['',3.6922003,-0.2865651,6.56];
ybs[5275]=['χ Cen',3.6962701,-0.7201693,4.36];
ybs[5276]=['',3.6969142,-0.7535433,6.2];
ybs[5277]=['π Hya',3.6973567,-0.467144,3.27];
ybs[5278]=['θ Cen',3.6989275,-0.6362214,2.06];
ybs[5279]=['',3.7068502,-1.1046272,6.4];
ybs[5280]=['95 Vir',3.698545,-0.1639936,5.46];
ybs[5281]=['α Dra',3.6864218,1.1221151,3.65];
ybs[5282]=['',3.709641,-1.0360089,6.34];
ybs[5283]=['',3.7175656,-1.2284925,6.05];
ybs[5284]=['',3.7086497,-0.7601505,6.17];
ybs[5285]=['',3.7197476,-1.2182657,6.06];
ybs[5286]=['',3.7120719,-0.9003605,6];
ybs[5287]=['',3.7135981,-0.9341215,4.75];
ybs[5288]=['96 Vir',3.7085852,-0.1818063,6.47];
ybs[5289]=['',3.7028399,0.7639634,5.27];
ybs[5290]=['13 Boo',3.7042131,0.8617663,5.25];
ybs[5291]=['',3.7166809,-0.2859514,4.91];
ybs[5292]=['',3.7058283,1.0342016,6.46];
ybs[5293]=['η Aps',3.754816,-1.4152479,4.91];
ybs[5294]=['12 Boo',3.7140469,0.4365013,4.83];
ybs[5295]=['3 UMi',3.696148,1.3004582,6.45];
ybs[5296]=['',3.7473257,-1.3568925,6.47];
ybs[5297]=['',3.7193491,0.0223491,6.43];
ybs[5298]=['',3.7282977,-0.938064,5.56];
ybs[5299]=['',3.7236693,-0.4266573,6.34];
ybs[5300]=['',3.7176166,0.5622362,6.11];
ybs[5301]=['',3.7300559,-0.9548175,6.11];
ybs[5302]=['50 Hya',3.7252921,-0.4772171,5.08];
ybs[5303]=['',3.7225753,0.0406295,5.01];
ybs[5304]=['',3.7272583,-0.46589,6.24];
ybs[5305]=['κ Vir',3.7255474,-0.1807289,4.19];
ybs[5306]=['',3.7358038,-0.9977472,5.07];
ybs[5307]=['',3.7288097,-0.0161751,5.91];
ybs[5308]=['',3.7341378,-0.7316149,5.61];
ybs[5309]=['',3.7373835,-0.9353307,6.39];
ybs[5310]=['',3.7439446,-1.1635801,5.75];
ybs[5311]=['4 UMi',3.7036229,1.3520204,4.82];
ybs[5312]=['',3.7318394,-0.1052181,6.36];
ybs[5313]=['14 Boo',3.7303536,0.2247691,5.54];
ybs[5314]=['',3.7351818,-0.512478,6.08];
ybs[5315]=['',3.7383299,-0.7868216,6.31];
ybs[5316]=['',3.743062,-1.0470993,6.39];
ybs[5317]=['',3.7836143,-1.447349,6.42];
ybs[5318]=['κ1 Boo',3.7266539,0.9024483,6.69];
ybs[5319]=['κ2 Boo',3.7267482,0.9024921,4.54];
ybs[5320]=['15 Boo',3.7337178,0.1748753,5.29];
ybs[5321]=['',3.7339912,0.0568136,6.45];
ybs[5322]=['',3.7366206,-0.3190747,5.43];
ybs[5323]=['',3.7327944,0.3803482,6.39];
ybs[5324]=['',3.719289,1.2104005,5.24];
ybs[5325]=['',3.7310666,0.7232264,6.24];
ybs[5326]=['ε Aps',3.7724865,-1.3995418,5.06];
ybs[5327]=['',3.7408764,-0.581578,6.55];
ybs[5328]=['ι Vir',3.7390844,-0.1061373,4.08];
ybs[5329]=['δ Oct',3.7958461,-1.4616343,4.32];
ybs[5330]=['α Boo',3.7371116,0.3333884,-0.04];
ybs[5331]=['',3.7405933,-0.116981,6.44];
ybs[5332]=['',3.7411611,-0.0571932,6.15];
ybs[5333]=['',3.7388975,0.328668,5.98];
ybs[5334]=['',3.7438815,-0.3257776,6.22];
ybs[5335]=['',3.7345481,0.9155119,6.58];
ybs[5336]=['',3.7409536,0.3497787,6.25];
ybs[5337]=['',3.7398846,0.6922697,6.38];
ybs[5338]=['',3.7500109,-0.5812056,6.54];
ybs[5339]=['',3.7575414,-1.0708073,5.23];
ybs[5340]=['ι Boo',3.7384463,0.8951195,4.75];
ybs[5341]=['λ Boo',3.7395998,0.8029865,4.18];
ybs[5342]=['',3.7450848,0.2649936,5.8];
ybs[5343]=['',3.747817,-0.1330409,6.47];
ybs[5344]=['ι Lup',3.7547811,-0.805253,3.55];
ybs[5345]=['',3.7507594,-0.3280544,5.9];
ybs[5346]=['',3.7525374,-0.4519614,5.87];
ybs[5347]=['',3.7544819,-0.647228,5.94];
ybs[5348]=['',3.7592676,-0.9855219,4.33];
ybs[5349]=['λ Vir',3.75272,-0.2347646,4.52];
ybs[5350]=['',3.7436177,0.8940772,6.2];
ybs[5351]=['',3.7469424,0.6183569,4.81];
ybs[5352]=['',3.7579935,-0.7529093,5.56];
ybs[5353]=['',3.7457901,0.8363863,6.32];
ybs[5354]=['',3.760453,-0.7900533,4.77];
ybs[5355]=['18 Boo',3.7529693,0.2255713,5.41];
ybs[5356]=['υ Vir',3.7544108,-0.0409343,5.14];
ybs[5357]=['ψ Cen',3.7595654,-0.6626111,4.05];
ybs[5358]=['',3.7549765,0.0053126,6.19];
ybs[5359]=['',3.7509167,0.6752246,6.86];
ybs[5360]=['20 Boo',3.7550102,0.2832177,4.86];
ybs[5361]=['',3.7694124,-1.0216887,4.92];
ybs[5362]=['',3.7503255,0.9561641,6.53];
ybs[5363]=['',3.7546933,0.6756888,6.33];
ybs[5364]=['',3.7564303,0.5296984,6.44];
ybs[5365]=['',3.7690313,-0.844727,6.09];
ybs[5366]=['',3.7672213,-0.6085274,5.56];
ybs[5367]=['',3.7721631,-0.8875184,6.02];
ybs[5368]=['',3.7704447,-0.6909959,4.42];
ybs[5369]=['',3.7811771,-1.1915997,5.61];
ybs[5370]=['',3.7743212,-0.9294818,6];
ybs[5371]=['51 Hya',3.7704087,-0.4857741,4.77];
ybs[5372]=['',3.7833602,-1.1563078,6.36];
ybs[5373]=['2 Lib',3.7715308,-0.2058269,6.21];
ybs[5374]=['',3.7705531,0.0202942,6.27];
ybs[5375]=['',3.7709574,0.1460166,6.86];
ybs[5376]=['',3.7709647,0.1460456,5.12];
ybs[5377]=['',3.7694849,0.4408547,6.22];
ybs[5378]=['',3.7737456,0.1425045,5.95];
ybs[5379]=['',3.802734,-1.3405231,6.07];
ybs[5380]=['',3.7778254,-0.4343231,5.32];
ybs[5381]=['',3.7897123,-1.1501636,5.85];
ybs[5382]=['',3.7745448,0.1002052,5.1];
ybs[5383]=['',3.7770076,-0.2050457,6.49];
ybs[5384]=['',3.7750136,0.1397374,6.19];
ybs[5385]=['τ1 Lup',3.7841726,-0.7906259,4.56];
ybs[5386]=['τ2 Lup',3.7843676,-0.7933843,4.35];
ybs[5387]=['',3.7807234,-0.3499044,6.61];
ybs[5388]=['',3.7844523,-0.7399723,6.32];
ybs[5389]=['',3.7821707,-0.4700301,6.48];
ybs[5390]=['',3.7870332,-0.6972914,6.35];
ybs[5391]=['',3.7888667,-0.806557,5.83];
ybs[5392]=['',3.7795119,0.6687176,6.27];
ybs[5393]=['',3.7961742,-1.034548,6.45];
ybs[5394]=['θ Boo',3.777764,0.9035987,4.05];
ybs[5395]=['22 Boo',3.7841901,0.3342104,5.39];
ybs[5396]=['104 Vir',3.7887995,-0.1081776,6.17];
ybs[5397]=['52 Hya',3.7926163,-0.5160817,4.97];
ybs[5398]=['',3.8081278,-1.1832284,5.83];
ybs[5399]=['φ Vir',3.7921986,-0.040242,4.81];
ybs[5400]=['106 Vir',3.7944359,-0.1217902,5.42];
ybs[5401]=['',3.7880119,0.7146624,6.63];
ybs[5402]=['',3.8016825,-0.7923581,5.5];
ybs[5403]=['',3.8027455,-0.8656173,5.37];
ybs[5404]=['',3.7930159,0.4923851,7.62];
ybs[5405]=['',3.7931467,0.4924144,7.12];
ybs[5406]=['',3.79172,0.6304007,6.1];
ybs[5407]=['',3.8050172,-0.714222,6.39];
ybs[5408]=['',3.7992995,0.0131191,5.94];
ybs[5409]=['',3.805997,-0.6797458,5.97];
ybs[5410]=['24 Boo',3.7928091,0.8686007,5.59];
ybs[5411]=['',3.8127419,-0.9942135,6.93];
ybs[5412]=['',3.7986142,0.5535115,6.06];
ybs[5413]=['',3.7973936,0.7281257,6.35];
ybs[5414]=['',3.8032199,0.0819474,6.02];
ybs[5415]=['σ Lup',3.812693,-0.8819742,4.42];
ybs[5416]=['',3.816966,-0.961232,5.87];
ybs[5417]=['',3.8166603,-0.9207697,5.87];
ybs[5418]=['',3.8144214,-0.5373956,6.09];
ybs[5419]=['ρ Boo',3.8073762,0.5287417,3.58];
ybs[5420]=['5 UMi',3.7852144,1.3197862,4.25];
ybs[5421]=['',3.8190011,-0.7361064,6.6];
ybs[5422]=['',3.8249291,-1.0487956,6.4];
ybs[5423]=['',3.8096731,0.4642688,6.01];
ybs[5424]=['26 Boo',3.8106648,0.3871746,5.92];
ybs[5425]=['γ Boo',3.8082452,0.6672685,3.03];
ybs[5426]=['',3.8013011,1.1014514,6.09];
ybs[5427]=['',3.805654,1.049794,6.27];
ybs[5428]=['',3.8215431,-0.3580554,6.5];
ybs[5429]=['',3.825061,-0.7259333,5.87];
ybs[5430]=['η Cen',3.8250025,-0.7371132,2.31];
ybs[5431]=['',3.8137768,0.6437318,6.43];
ybs[5432]=['',3.8094565,0.9655372,5.76];
ybs[5433]=['',3.8365954,-1.1869505,6.04];
ybs[5434]=['',3.8286908,-0.8084497,5.55];
ybs[5435]=['',3.8176324,0.566505,6.33];
ybs[5436]=['',3.8288353,-0.6924189,6.13];
ybs[5437]=['σ Boo',3.8198176,0.5178223,4.46];
ybs[5438]=['',3.8194676,0.6379153,6.03];
ybs[5439]=['',3.8303111,-0.7031414,5.74];
ybs[5440]=['',3.8331359,-0.8065008,5.41];
ybs[5441]=['',3.8169778,0.9946487,6.48];
ybs[5442]=['',3.8190992,0.8603137,5.74];
ybs[5443]=['ρ Lup',3.8356847,-0.8639533,4.05];
ybs[5444]=['',3.8262261,0.4044748,6.38];
ybs[5445]=['',3.8307693,-0.2160822,6.2];
ybs[5446]=['',3.8372113,-0.6783939,6.02];
ybs[5447]=['',3.8412298,-0.8143508,6.07];
ybs[5448]=['',3.8423323,-0.8574836,6.39];
ybs[5449]=['α1 Cen',3.8438496,-1.0630772,-0.01];
ybs[5450]=['α2 Cen',3.8438642,-1.063082,1.33];
ybs[5451]=['',3.8476833,-0.9863753,6.3];
ybs[5452]=['',3.8355712,0.3180571,5.91];
ybs[5453]=['α Cir',3.8569597,-1.1353196,3.19];
ybs[5454]=['',3.8347871,0.760386,5.7];
ybs[5455]=['',3.8538702,-1.0243343,6.22];
ybs[5456]=['',3.8489086,-0.6319698,5.67];
ybs[5457]=['',3.8345245,0.9415754,5.85];
ybs[5458]=['33 Boo',3.8374974,0.7736972,5.39];
ybs[5459]=['α Lup',3.8532643,-0.8283729,2.3];
ybs[5460]=['α Aps',3.8838597,-1.3808495,3.83];
ybs[5461]=['',3.8530498,-0.6609136,4];
ybs[5462]=['',3.8447885,0.3822468,6.1];
ybs[5463]=['',3.8464574,0.2349186,5.91];
ybs[5464]=['',3.8523814,-0.5411797,6.37];
ybs[5465]=['π1 Boo',3.8464853,0.2852568,4.94];
ybs[5466]=['π2 Boo',3.8465072,0.2852472,5.88];
ybs[5467]=['ζ Boo',3.84838,0.2383095,4.83];
ybs[5468]=['ζ Boo',3.84838,0.2383095,4.43];
ybs[5469]=['',3.8101179,1.3890007,6.26];
ybs[5470]=['31 Boo',3.8506651,0.1411553,4.86];
ybs[5471]=['32 Boo',3.8509373,0.2022228,5.56];
ybs[5472]=['',3.8689155,-1.0986639,5.36];
ybs[5473]=['',3.8515162,0.367385,6.38];
ybs[5474]=['4 Lib',3.8582371,-0.4375733,5.73];
ybs[5475]=['',3.8603811,-0.6151777,4.05];
ybs[5476]=['',3.8669631,-1.0219054,6.11];
ybs[5477]=['μ Vir',3.8571006,-0.1000422,3.88];
ybs[5478]=['',3.8678936,-0.9717114,6.1];
ybs[5479]=['',3.8661901,-0.6154912,4.92];
ybs[5480]=['34 Boo',3.8580196,0.4617128,4.81];
ybs[5481]=['',4.0970187,-1.5384088,6.48];
ybs[5482]=['',3.8505836,1.0679309,6.25];
ybs[5483]=['',3.8590061,0.7048627,5.73];
ybs[5484]=['',3.8731692,-0.8292719,5.74];
ybs[5485]=['',3.8757517,-0.915532,5.21];
ybs[5486]=['',3.8663339,-0.0260202,6.07];
ybs[5487]=['54 Hya',3.8703634,-0.4453361,4.94];
ybs[5488]=['',3.8765434,-0.9124234,6.07];
ybs[5489]=['',3.8707965,-0.4053676,5.81];
ybs[5490]=['',3.8844254,-1.1635388,5.91];
ybs[5491]=['108 Vir',3.8676373,0.0112441,5.69];
ybs[5492]=['ο Boo',3.8661666,0.2948101,4.6];
ybs[5493]=['5 Lib',3.8699709,-0.2710945,6.33];
ybs[5494]=['',3.8710493,-0.3708631,6.4];
ybs[5495]=['ε Boo',3.864823,0.4712712,5.12];
ybs[5496]=['ε Boo',3.8648231,0.4712567,2.7];
ybs[5497]=['',3.8665757,0.3283258,6.13];
ybs[5498]=['',3.8754499,-0.6695666,5.94];
ybs[5499]=['',3.8775994,-0.7614849,6.3];
ybs[5500]=['',3.8657275,0.5709889,6.28];
ybs[5501]=['109 Vir',3.8708648,0.0317651,3.72];
ybs[5502]=['',3.8699494,0.2628313,5.63];
ybs[5503]=['',3.8756218,-0.3734516,6.06];
ybs[5504]=['55 Hya',3.8763637,-0.4484951,5.63];
ybs[5505]=['',3.8851567,-0.9902942,6.23];
ybs[5506]=['56 Hya',3.8779978,-0.4565752,5.24];
ybs[5507]=['57 Hya',3.8789357,-0.4663285,5.77];
ybs[5508]=['',3.8784342,-0.2253571,6.35];
ybs[5509]=['',3.8821671,-0.6406543,6.04];
ybs[5510]=['',3.9052426,-1.2786393,5.6];
ybs[5511]=['',3.8847857,-0.424526,5.68];
ybs[5512]=['',3.8824964,-0.0160534,6.14];
ybs[5513]=['μ Lib',3.8845869,-0.2481994,5.31];
ybs[5514]=['',3.8797225,0.4240188,6.14];
ybs[5515]=['π1 Oct',3.949292,-1.4537828,5.65];
ybs[5516]=['58 Hya',3.8891423,-0.4892485,4.41];
ybs[5517]=['',3.9009189,-1.1149316,5.87];
ybs[5518]=['ο Lup',3.8955328,-0.7617794,4.32];
ybs[5519]=['',3.882507,0.658672,6.16];
ybs[5520]=['α1 Lib',3.8905992,-0.280452,5.15];
ybs[5521]=['α2 Lib',3.8914367,-0.2812267,2.75];
ybs[5522]=['',3.8865278,0.4981886,5.8];
ybs[5523]=['38 Boo',3.8830769,0.8036224,5.74];
ybs[5524]=['',3.88791,0.4160918,5.85];
ybs[5525]=['11 Lib',3.8917543,-0.0413743,4.94];
ybs[5526]=['',3.8916476,-0.0057358,6.18];
ybs[5527]=['',3.883822,0.8954035,6.51];
ybs[5528]=['39 Boo',3.8846043,0.8490803,5.69];
ybs[5529]=['ζ Cir',3.9104577,-1.152993,6.09];
ybs[5530]=['',3.9267225,-1.339226,5.34];
ybs[5531]=['',3.8885471,0.6492686,5.48];
ybs[5532]=['',3.8990944,-0.5349065,6.29];
ybs[5533]=['',3.9006205,-0.6610291,5.03];
ybs[5534]=['ξ Boo',3.8929284,0.3321326,4.55];
ybs[5535]=['π2 Oct',3.9619211,-1.4504607,5.65];
ybs[5536]=['',3.9136646,-1.0504119,5.2];
ybs[5537]=['',3.9370813,-1.3478967,5.93];
ybs[5538]=['12 Lib',3.906726,-0.431317,5.3];
ybs[5539]=['',3.9082619,-0.5824318,5.82];
ybs[5540]=['',3.9017275,0.2728598,6.4];
ybs[5541]=['θ Cir',3.9189658,-1.0969479,5.11];
ybs[5542]=['',3.8916075,1.0336283,5.46];
ybs[5543]=['',3.9016793,0.3330445,6.01];
ybs[5544]=['ξ1 Lib',3.9066381,-0.2088942,5.8];
ybs[5545]=['',3.9351395,-1.3107619,6.2];
ybs[5546]=['',3.9162743,-0.922922,5.38];
ybs[5547]=['ω Oct',3.9935088,-1.4809523,5.91];
ybs[5548]=['',3.9131329,-0.5921175,5.32];
ybs[5549]=['',3.9170927,-0.8368664,5.64];
ybs[5550]=['',3.9194112,-0.8991333,6.64];
ybs[5551]=['',3.9170316,-0.6891581,6.36];
ybs[5552]=['',3.9164592,-0.5708351,6.06];
ybs[5553]=['β UMi',3.8862958,1.2930084,2.08];
ybs[5554]=['ξ2 Lib',3.917044,-0.200354,5.46];
ybs[5555]=['',3.9194686,-0.5101133,6.29];
ybs[5556]=['',3.9241883,-0.8540301,6.35];
ybs[5557]=['',3.9141066,0.2509172,5.77];
ybs[5558]=['',3.9203188,-0.374985,5.74];
ybs[5559]=['',3.9126096,0.5625248,6.12];
ybs[5560]=['16 Lib',3.9187055,-0.0770735,4.49];
ybs[5561]=['β Lup',3.925621,-0.7540355,2.68];
ybs[5562]=['',3.925857,-0.6977094,6.15];
ybs[5563]=['',3.9202402,-0.0041362,5.53];
ybs[5564]=['',3.9176073,0.3749949,6.49];
ybs[5565]=['',3.9183161,0.2848155,5.71];
ybs[5566]=['κ Cen',3.9283331,-0.7360601,3.13];
ybs[5567]=['59 Hya',3.9256572,-0.4839163,5.65];
ybs[5568]=['17 Lib',3.9233886,-0.1959007,6.6];
ybs[5569]=['',3.9284875,-0.6623585,6.47];
ybs[5570]=['',3.9296409,-0.7544863,6.1];
ybs[5571]=['',3.9137318,0.864963,5.63];
ybs[5572]=['18 Lib',3.9263126,-0.1957081,5.87];
ybs[5573]=['',3.9261228,-0.0882832,6.09];
ybs[5574]=['',3.9281256,0.0785195,5.93];
ybs[5575]=['',3.9371628,-0.6654363,5.89];
ybs[5576]=['δ Lib',3.9353333,-0.1498775,4.92];
ybs[5577]=['',3.9403249,-0.6008648,6.22];
ybs[5578]=['40 Boo',3.9282429,0.6841057,5.64];
ybs[5579]=['',3.9176987,1.1495249,4.6];
ybs[5580]=['',3.9367665,-0.0492768,5.52];
ybs[5581]=['60 Hya',3.9407456,-0.4909376,5.85];
ybs[5582]=['',3.9342263,0.3835718,6.38];
ybs[5583]=['η Cir',3.9543943,-1.1187364,5.17];
ybs[5584]=['',3.9388274,-0.0036437,5.71];
ybs[5585]=['',3.9447264,-0.5709174,5.44];
ybs[5586]=['',3.8800653,1.4388504,5.64];
ybs[5587]=['',3.9324129,0.8239555,6.37];
ybs[5588]=['',3.9658943,-1.2561433,6.52];
ybs[5589]=['',3.9429615,-0.0540932,6.61];
ybs[5590]=['ω Boo',3.9395376,0.4352837,4.81];
ybs[5591]=['110 Vir',3.9435142,0.0353169,4.4];
ybs[5592]=['β Boo',3.9383665,0.7037578,3.5];
ybs[5593]=['σ Lib',3.949243,-0.4424313,3.29];
ybs[5594]=['',3.9525425,-0.7143401,6.41];
ybs[5595]=['π Lup',3.9545696,-0.8223688,4.72];
ybs[5596]=['π Lup',3.9545696,-0.8223688,4.82];
ybs[5597]=['',3.9551855,-0.7179293,5.15];
ybs[5598]=['',3.9351263,1.0495721,5.93];
ybs[5599]=['',3.9435779,0.6132736,5.51];
ybs[5600]=['',3.9487059,0.0946839,6.5];
ybs[5601]=['',3.9682909,-1.1404293,6.17];
ybs[5602]=['',3.9432706,0.7780082,6.65];
ybs[5603]=['',3.9458071,0.602111,6.59];
ybs[5604]=['',3.9567875,-0.4512844,6.67];
ybs[5605]=['',3.959004,-0.6340954,6.27];
ybs[5606]=['ψ Boo',3.9496762,0.4691458,4.54];
ybs[5607]=['',3.9647681,-0.8579125,5.77];
ybs[5608]=['44 Boo',3.9461003,0.8305462,4.76];
ybs[5609]=['',3.960249,-0.540801,5.96];
ybs[5610]=['',3.959557,-0.3856953,6.17];
ybs[5611]=['',3.9753224,-1.1719868,5.76];
ybs[5612]=['ν Lib',3.9601807,-0.2849017,5.2];
ybs[5613]=['',3.9745642,-1.1119239,6.28];
ybs[5614]=['',3.9677607,-0.709478,5.79];
ybs[5615]=['',3.9698158,-0.7493369,5.85];
ybs[5616]=['λ Lup',3.9707552,-0.7914321,4.05];
ybs[5617]=['47 Boo',3.9532232,0.8392234,5.57];
ybs[5618]=['',3.9898288,-1.2712105,6.01];
ybs[5619]=['',3.9454638,1.1493355,6.13];
ybs[5620]=['',3.9587232,0.635104,6.35];
ybs[5621]=['',3.9642593,0.0948001,6.16];
ybs[5622]=['',3.9800991,-1.0731657,6.3];
ybs[5623]=['',3.9625176,0.3207067,6.02];
ybs[5624]=['45 Boo',3.962189,0.4328874,4.93];
ybs[5625]=['',3.9565642,0.9510209,5.25];
ybs[5626]=['',3.9760947,-0.6782019,5.98];
ybs[5627]=['',3.9818942,-0.9671101,5.54];
ybs[5628]=['46 Boo',3.966926,0.4578855,5.67];
ybs[5629]=['',3.969401,0.2298417,6.1];
ybs[5630]=['',3.9678083,0.4370735,5.81];
ybs[5631]=['',3.9765137,-0.460738,5.76];
ybs[5632]=['',3.9827172,-0.7913783,6.44];
ybs[5633]=['',3.9824989,-0.7914125,7.39];
ybs[5634]=['',3.9969161,-1.2242372,5.81];
ybs[5635]=['',3.9900706,-1.0787624,6.32];
ybs[5636]=['κ1 Lup',3.9844257,-0.8517695,3.87];
ybs[5637]=['κ2 Lup',3.9845352,-0.8518711,5.69];
ybs[5638]=['',3.9657421,0.8724682,6.39];
ybs[5639]=['ζ Lup',3.9861387,-0.9104347,3.41];
ybs[5640]=['',3.9869626,-0.8427051,6.33];
ybs[5641]=['',3.988114,-0.7778112,4.82];
ybs[5642]=['ι1 Lib',3.9846883,-0.3465635,4.54];
ybs[5643]=['',3.9890941,-0.6310421,6.1];
ybs[5644]=['',3.9831424,0.3300556,5.89];
ybs[5645]=['',3.9894679,-0.4201523,6.47];
ybs[5646]=['ι2 Lib',3.989479,-0.3440413,6.08];
ybs[5647]=['23 Lib',3.9903177,-0.4428551,6.45];
ybs[5648]=['',3.992132,-0.4582893,5.84];
ybs[5649]=['',3.9859843,0.3354697,6.68];
ybs[5650]=['1 Lup',3.9954889,-0.5512335,4.91];
ybs[5651]=['',4.0056852,-1.0640816,5.73];
ybs[5652]=['26 Lib',3.9948522,-0.3112417,6.17];
ybs[5653]=['',4.0017,-0.8401652,5.95];
ybs[5654]=['δ Cir',4.0071695,-1.0650154,5.09];
ybs[5655]=['',3.9894083,0.4000076,6.3];
ybs[5656]=['ε Cir',4.0105232,-1.1113157,4.86];
ybs[5657]=['',4.0021475,-0.7252685,5.16];
ybs[5658]=['',4.0027032,-0.7600629,6.04];
ybs[5659]=['',3.9958018,-0.0971611,6.28];
ybs[5660]=['β Cir',4.0094601,-1.0273762,4.07];
ybs[5661]=['γ TrA',4.0167545,-1.1997769,2.89];
ybs[5662]=['',3.9746702,1.1818642,6.17];
ybs[5663]=['',3.9892133,0.6667182,6.2];
ybs[5664]=['',3.9916447,0.553682,5.99];
ybs[5665]=['3 Ser',3.9970732,0.0850919,5.33];
ybs[5666]=['χ Boo',3.9934022,0.5078887,5.26];
ybs[5667]=['',3.9915854,0.7349056,6.13];
ybs[5668]=['',4.0029217,-0.3920547,5.5];
ybs[5669]=['4 Ser',3.9999172,0.0053824,5.63];
ybs[5670]=['',4.0153012,-1.0569569,5.46];
ybs[5671]=['δ Boo',3.9977084,0.5803352,3.47];
ybs[5672]=['',4.0112585,-0.7177469,6.28];
ybs[5673]=['μ Lup',4.0132087,-0.8366743,4.27];
ybs[5674]=['',4.0242944,-1.1788569,6.28];
ybs[5675]=['β Lib',4.0053239,-0.1648725,2.61];
ybs[5676]=['2 Lup',4.0094626,-0.5272997,4.34];
ybs[5677]=['',4.0146689,-0.7129865,5.59];
ybs[5678]=['',4.0132355,-0.545805,6.18];
ybs[5679]=['',4.017092,-0.6485515,6.2];
ybs[5680]=['',4.0113607,-0.0091521,5.89];
ybs[5681]=['',3.9917208,1.1742985,5.13];
ybs[5682]=['',4.0107399,0.3579629,5.7];
ybs[5683]=['',3.9924505,1.2022005,6.51];
ybs[5684]=['5 Ser',4.0151352,0.0297155,5.06];
ybs[5685]=['δ Lup',4.0252833,-0.7105146,3.22];
ybs[5686]=['',4.0262337,-0.7122975,6.2];
ybs[5687]=['',4.0257576,-0.6681315,6.48];
ybs[5688]=['ν1 Lup',4.0289601,-0.8375748,5];
ybs[5689]=['ν2 Lup',4.0275167,-0.8443835,5.65];
ybs[5690]=['',4.0343699,-1.0597341,5.67];
ybs[5691]=['28 Lib',4.0225062,-0.3180127,6.17];
ybs[5692]=['',4.0151748,0.5664049,6.32];
ybs[5693]=['ο Lib',4.0230011,-0.272454,6.3];
ybs[5694]=['γ Cir',4.0351424,-1.0364135,4.51];
ybs[5695]=['φ1 Lup',4.0270187,-0.6339601,3.56];
ybs[5696]=['',4.0216519,-0.0432066,6.35];
ybs[5697]=['',4.0232212,-0.1027446,5.54];
ybs[5698]=['ε Lup',4.0311771,-0.7810522,3.37];
ybs[5699]=['ο CrB',4.018048,0.5158086,5.51];
ybs[5700]=['6 Ser',4.0226639,0.0113994,5.35];
ybs[5701]=['',4.0224272,0.4345109,6.39];
ybs[5702]=['φ2 Lup',4.0329351,-0.6443758,4.54];
ybs[5703]=['',4.048783,-1.1932722,5.89];
ybs[5704]=['11 UMi',4.0015922,1.2524539,5.02];
ybs[5705]=['',4.0168839,0.9057577,5.66];
ybs[5706]=['',4.0199156,0.7744351,6.19];
ybs[5707]=['7 Ser',4.028292,0.2182673,6.28];
ybs[5708]=['50 Boo',4.025219,0.5737241,5.37];
ybs[5709]=['υ Lup',4.0399931,-0.6941376,5.37];
ybs[5710]=['',4.0353477,-0.2169556,5.72];
ybs[5711]=['8 Ser',4.034463,-0.0189151,6.12];
ybs[5712]=['',4.0415135,-0.6672429,7.03];
ybs[5713]=['ε Lib',4.0367315,-0.181223,4.94];
ybs[5714]=['',4.0425247,-0.6770881,4.6];
ybs[5715]=['',4.0539413,-1.1273302,5.71];
ybs[5716]=['',4.02853,0.6897494,5.5];
ybs[5717]=['η CrB',4.0313855,0.5275488,5.58];
ybs[5718]=['η CrB',4.0313855,0.5275488,6.08];
ybs[5719]=['ρ Oct',4.1344294,-1.4751394,5.57];
ybs[5720]=['κ1 Aps',4.0729753,-1.2819079,5.49];
ybs[5721]=['',4.0271336,1.0818512,5.98];
ybs[5722]=['',4.0346606,0.789062,6.01];
ybs[5723]=['μ1 Boo',4.0367508,0.6512902,4.31];
ybs[5724]=['μ2 Boo',4.0368609,0.6507765,6.5];
ybs[5725]=['γ UMi',4.0173757,1.252649,3.05];
ybs[5726]=['',4.0510395,-0.6427665,5.45];
ybs[5727]=['',4.0270726,1.1044387,5.79];
ybs[5728]=['',4.056754,-0.9015868,6.1];
ybs[5729]=['τ1 Ser',4.0430677,0.268213,5.17];
ybs[5730]=['',4.0433931,0.3389478,6.27];
ybs[5731]=['',4.0447061,0.5982184,5.46];
ybs[5732]=['',4.0606426,-0.8166759,5.24];
ybs[5733]=['ζ1 Lib',4.0546015,-0.2927985,5.64];
ybs[5734]=['ι Dra',4.0374796,1.0280892,3.29];
ybs[5735]=['',4.050916,0.4370599,6.02];
ybs[5736]=['10 Ser',4.0558128,0.0311122,5.17];
ybs[5737]=['β CrB',4.0515827,0.5069468,3.68];
ybs[5738]=['',4.0448783,0.9417729,6.45];
ybs[5739]=['',4.064957,-0.3628063,6.22];
ybs[5740]=['ζ3 Lib',4.0651467,-0.2909178,5.82];
ybs[5741]=['',4.071939,-0.675119,6.25];
ybs[5742]=['',4.0548629,0.8227787,6.15];
ybs[5743]=['',4.0707053,-0.5749048,6.46];
ybs[5744]=['',4.0491362,1.0858656,6.5];
ybs[5745]=['',4.0500698,1.0578494,5.9];
ybs[5746]=['',4.0698234,-0.3529628,6.22];
ybs[5747]=['',4.1089376,-1.3608989,6.18];
ybs[5748]=['',4.0656291,0.1487075,6.57];
ybs[5749]=['',4.0552749,0.9622945,6.43];
ybs[5750]=['',4.062634,0.5450148,6.46];
ybs[5751]=['',4.062822,0.6413234,6.37];
ybs[5752]=['',4.0736876,-0.3443327,5.52];
ybs[5753]=['ν1 Boo',4.0646925,0.7116435,5.02];
ybs[5754]=['ζ4 Lib',4.074959,-0.2951514,5.5];
ybs[5755]=['',4.0761786,-0.4284348,7];
ybs[5756]=['',4.0922983,-1.1461617,6.51];
ybs[5757]=['',4.0805268,-0.7002934,5.82];
ybs[5758]=['',4.0564419,1.0828066,6.38];
ybs[5759]=['',4.0667892,0.6380514,6.38];
ybs[5760]=['τ2 Ser',4.0708453,0.279212,6.22];
ybs[5761]=['ε TrA',4.0942813,-1.1584393,4.11];
ybs[5762]=['11 Ser',4.0747703,-0.0217213,5.51];
ybs[5763]=['',4.0819025,-0.6877833,6.36];
ybs[5764]=['ν2 Boo',4.0684098,0.7128072,5.02];
ybs[5765]=['36 Lib',4.0827016,-0.4905161,5.15];
ybs[5766]=['γ Lup',4.0854419,-0.7194998,2.78];
ybs[5767]=['37 Lib',4.0802792,-0.1766653,4.62];
ybs[5768]=['θ CrB',4.0737584,0.546305,4.14];
ybs[5769]=['',4.0809073,-0.1004031,6.51];
ybs[5770]=['',4.0814141,-0.1612854,5.17];
ybs[5771]=['',4.0888666,-0.7856723,4.54];
ybs[5772]=['κ2 Aps',4.1116928,-1.2828534,5.65];
ybs[5773]=['',4.0783126,0.2981008,6.45];
ybs[5774]=['',4.0902157,-0.7758676,5.43];
ybs[5775]=['',4.0630567,1.1196228,5.79];
ybs[5776]=['',4.119455,-1.3288374,5.95];
ybs[5777]=['γ Lib',4.0862826,-0.2591239,3.91];
ybs[5778]=['δ Ser',4.0825007,0.1829099,3.8];
ybs[5779]=['δ Ser',4.0825006,0.182939,3.8];
ybs[5780]=['',4.0897126,-0.5785729,6.24];
ybs[5781]=['',4.0839158,0.0281252,6.56];
ybs[5782]=['',4.1101671,-1.2266746,6.44];
ybs[5783]=['α CrB',4.0815713,0.4652548,2.23];
ybs[5784]=['υ Lib',4.0931965,-0.4920385,3.58];
ybs[5785]=['τ3 Ser',4.085598,0.3071478,6.12];
ybs[5786]=['',4.0872304,0.1956233,6.07];
ybs[5787]=['ω Lup',4.0982187,-0.7439265,4.33];
ybs[5788]=['',4.1021162,-0.9150559,5.44];
ybs[5789]=['14 Ser',4.0904489,-0.0107963,6.51];
ybs[5790]=['μ CrB',4.0835998,0.6798508,5.11];
ybs[5791]=['',4.0951044,-0.4596599,6.19];
ybs[5792]=['16 Ser',4.0898937,0.1737135,5.26];
ybs[5793]=['',4.107594,-1.046569,5.95];
ybs[5794]=['τ5 Ser',4.0897168,0.2803336,5.93];
ybs[5795]=['',4.1002211,-0.6844662,6.57];
ybs[5796]=['',4.0964383,-0.4048838,5.78];
ybs[5797]=['',4.1009185,-0.6838931,6.04];
ybs[5798]=['',4.0861171,0.6687521,6.42];
ybs[5799]=['',4.0985975,-0.4932818,6.32];
ybs[5800]=['',4.0984299,-0.367783,5.84];
ybs[5801]=['',4.0829229,0.940113,5.97];
ybs[5802]=['τ Lib',4.1003718,-0.5207005,3.66];
ybs[5803]=['',4.0910711,0.5224516,6.52];
ybs[5804]=['41 Lib',4.1011616,-0.3378614,5.38];
ybs[5805]=['',4.0998351,-0.1544726,6.5];
ybs[5806]=['',4.0998424,-0.1544143,6.48];
ybs[5807]=['',4.0865169,0.9077905,6.74];
ybs[5808]=['',4.0858361,0.9524847,5.74];
ybs[5809]=['',4.1032191,-0.4050247,6.34];
ybs[5810]=['ψ1 Lup',4.1053752,-0.6015751,4.67];
ybs[5811]=['',4.1112186,-0.8341082,6.23];
ybs[5812]=['',4.1074042,-0.5457508,6.34];
ybs[5813]=['φ Boo',4.0947787,0.7033119,5.24];
ybs[5814]=['42 Lib',4.1072834,-0.416674,4.96];
ybs[5815]=['',4.1120077,-0.7804478,4.64];
ybs[5816]=['θ UMi',4.0619587,1.3489743,4.96];
ybs[5817]=['',4.0993066,0.6042124,6.11];
ybs[5818]=['',4.0927365,0.9503707,5.87];
ybs[5819]=['',4.050283,1.4030504,6.58];
ybs[5820]=['',4.0964136,0.8157912,5.75];
ybs[5821]=['',4.1058964,0.2093937,6.25];
ybs[5822]=['',4.1185171,-0.8647092,6.04];
ybs[5823]=['ζ1 CrB',4.1016545,0.6384533,6];
ybs[5824]=['ζ2 CrB',4.1016909,0.6384388,5.07];
ybs[5825]=['',4.0975276,0.8790708,5.84];
ybs[5826]=['',4.1250119,-1.0531574,6.48];
ybs[5827]=['',4.1180265,-0.6541452,5.24];
ybs[5828]=['κ Lib',4.1144338,-0.3444218,4.74];
ybs[5829]=['ψ2 Lup',4.1181209,-0.606769,4.75];
ybs[5830]=['τ6 Ser',4.1093414,0.2787173,6.01];
ybs[5831]=['',4.099549,1.009993,6.45];
ybs[5832]=['ι Ser',4.1117076,0.3423476,4.52];
ybs[5833]=['χ Ser',4.1129285,0.2232693,5.33];
ybs[5834]=['',4.0915394,1.2082329,5.62];
ybs[5835]=['τ7 Ser',4.1133109,0.3212944,5.81];
ybs[5836]=['',4.1258002,-0.7308319,5.94];
ybs[5837]=['',4.1207175,-0.2635071,6.31];
ybs[5838]=['η Lib',4.1236074,-0.274489,5.41];
ybs[5839]=['γ CrB',4.1167134,0.4579877,3.84];
ybs[5840]=['',4.1189564,0.2375943,6.48];
ybs[5841]=['',4.1430094,-1.1431089,6.18];
ybs[5842]=['',4.1430167,-1.1431089,6.39];
ybs[5843]=['ψ Ser',4.1229601,0.042947,5.88];
ybs[5844]=['α Ser',4.123903,0.1112005,2.65];
ybs[5845]=['π CrB',4.121941,0.5665596,5.56];
ybs[5846]=['',4.1333064,-0.4907024,6.51];
ybs[5847]=['',4.1160297,0.9129124,5.51];
ybs[5848]=['τ8 Ser',4.1255133,0.3003725,6.14];
ybs[5849]=['',4.1288283,0.0941275,5.58];
ybs[5850]=['',4.1358174,-0.6062543,5.61];
ybs[5851]=['',4.1301195,0.0146198,6.33];
ybs[5852]=['',4.1390337,-0.7024466,6.42];
ybs[5853]=['25 Ser',4.1320698,-0.0324286,5.4];
ybs[5854]=['',4.1392015,-0.6626916,6.01];
ybs[5855]=['',4.1458657,-0.9161335,6.07];
ybs[5856]=['',4.1350721,-0.10775,6.24];
ybs[5857]=['β Ser',4.1320478,0.2682288,3.67];
ybs[5858]=['λ Ser',4.1333688,0.1274019,4.43];
ybs[5859]=['',4.151496,-0.9295887,5.77];
ybs[5860]=['υ Ser',4.136881,0.2454301,5.71];
ybs[5861]=['',4.1505376,-0.8545844,5.84];
ybs[5862]=['',4.1517151,-0.7933167,6.12];
ybs[5863]=['',4.1559929,-0.961808,5.73];
ybs[5864]=['',4.1409611,0.2397345,6];
ybs[5865]=['',4.1445706,-0.0675647,5.53];
ybs[5866]=['',4.1650531,-1.1380156,6.54];
ybs[5867]=['',4.1396073,0.552971,6.44];
ybs[5868]=['',4.1320841,0.967283,5.92];
ybs[5869]=['κ Ser',4.1431047,0.3157127,4.09];
ybs[5870]=['',4.1420811,0.4905062,5.85];
ybs[5871]=['μ Ser',4.147499,-0.0607828,3.53];
ybs[5872]=['',4.1572762,-0.8222617,6.01];
ybs[5873]=['χ Lup',4.1542108,-0.5878098,3.95];
ybs[5874]=['',4.1666206,-1.0935848,6.19];
ybs[5875]=['1 Sco',4.1540348,-0.4503507,4.64];
ybs[5876]=['',4.1317874,1.0916328,5.19];
ybs[5877]=['',4.1366814,0.9655781,5.86];
ybs[5878]=['ω Ser',4.1502886,0.0374252,5.23];
ybs[5879]=['δ CrB',4.1466026,0.4540643,4.63];
ybs[5880]=['',4.1633076,-0.8842948,6.6];
ybs[5881]=['κ TrA',4.1767679,-1.1982226,5.09];
ybs[5882]=['ε Ser',4.1525219,0.0772462,3.71];
ybs[5883]=['',4.1595593,-0.522517,6.4];
ybs[5884]=['',4.1517147,0.2632245,5.2];
ybs[5885]=['36 Ser',4.1546468,-0.0548433,5.11];
ybs[5886]=['',4.1565902,-0.2475783,6.19];
ybs[5887]=['β TrA',4.1744224,-1.1079483,2.85];
ybs[5888]=['',4.172951,-1.0610494,6.15];
ybs[5889]=['ρ Ser',4.1540385,0.365228,4.76];
ybs[5890]=['',4.1757832,-1.0511744,5.77];
ybs[5891]=['κ CrB',4.1534122,0.6214369,4.82];
ybs[5892]=['λ Lib',4.1641552,-0.3528739,5.03];
ybs[5893]=['ζ UMi',4.1166474,1.3568169,4.32];
ybs[5894]=['2 Sco',4.1655165,-0.442931,4.59];
ybs[5895]=['',4.1782694,-1.056494,5.76];
ybs[5896]=['',4.1667362,-0.4290684,5.39];
ybs[5897]=['',4.1668647,-0.4193816,5.42];
ybs[5898]=['θ Lib',4.1661971,-0.2928703,4.15];
ybs[5899]=['',4.1614312,0.3028576,6.36];
ybs[5900]=['',4.1694579,-0.4780307,6.14];
ybs[5901]=['39 Ser',4.1627015,0.2294344,6.1];
ybs[5902]=['3 Sco',4.1700829,-0.4414651,5.87];
ybs[5903]=['',4.1642813,0.2796731,6.09];
ybs[5904]=['χ Her',4.1594209,0.7400264,4.62];
ybs[5905]=['47 Lib',4.1714218,-0.339177,5.94];
ybs[5906]=['',4.1739775,-0.5433868,6.21];
ybs[5907]=['4 Sco',4.1737962,-0.4593009,5.62];
ybs[5908]=['',4.1769436,-0.6966324,6.03];
ybs[5909]=['40 Ser',4.1692358,0.1488726,6.29];
ybs[5910]=['',4.1913675,-1.135975,5.75];
ybs[5911]=['',4.1815268,-0.8414544,6.31];
ybs[5912]=['',4.1568795,0.973461,5.81];
ybs[5913]=['',4.1771672,-0.5556382,6.29];
ybs[5914]=['',4.1685026,0.3536085,5.44];
ybs[5915]=['ξ1 Lup',4.1801275,-0.5936918,5.12];
ybs[5916]=['ξ2 Lup',4.1801784,-0.593653,5.62];
ybs[5917]=['',4.1766732,-0.2521888,6.37];
ybs[5918]=['ρ Sco',4.179931,-0.5107501,3.88];
ybs[5919]=['',4.1822352,-0.6324109,5.8];
ybs[5920]=['',4.1780669,-0.2596869,6.13];
ybs[5921]=['',4.1732909,0.3241146,6.26];
ybs[5922]=['2 Her',4.1679192,0.7520282,5.37];
ybs[5923]=['γ Ser',4.1768196,0.2724773,3.85];
ybs[5924]=['',4.1831146,-0.3670852,5.85];
ybs[5925]=['',4.1873381,-0.6554079,6.31];
ybs[5926]=['λ CrB',4.1732128,0.661424,5.45];
ybs[5927]=['',4.1943007,-0.9436928,6.1];
ybs[5928]=['4 Her',4.1717756,0.7420418,5.75];
ybs[5929]=['',4.2008528,-1.1139503,6.41];
ybs[5930]=['φ Ser',4.1802999,0.2507142,5.54];
ybs[5931]=['48 Lib',4.1851803,-0.250082,4.88];
ybs[5932]=['',4.1871885,-0.4342455,5.43];
ybs[5933]=['',4.1918595,-0.7294276,4.99];
ybs[5934]=['π Sco',4.1884147,-0.4566324,2.89];
ybs[5935]=['',4.1938257,-0.7103714,6.49];
ybs[5936]=['',4.199614,-0.9534006,6.13];
ybs[5937]=['ε CrB',4.1814391,0.4682424,4.15];
ybs[5938]=['η Lup',4.1944048,-0.6709988,3.41];
ybs[5939]=['',4.1720865,1.0273268,6.31];
ybs[5940]=['',4.1805736,0.691949,6.31];
ybs[5941]=['',4.2081035,-1.0923847,6.25];
ybs[5942]=['',4.1978695,-0.7065693,6.21];
ybs[5943]=['δ Sco',4.1947738,-0.3956673,2.32];
ybs[5944]=['49 Lib',4.1945678,-0.289406,5.47];
ybs[5945]=['',4.2233062,-1.2644383,5.7];
ybs[5946]=['',4.1994065,-0.5574139,6.33];
ybs[5947]=['',4.1870806,0.6387017,5.62];
ybs[5948]=['',4.1898219,0.4515431,2];
ybs[5949]=['50 Lib',4.1963844,-0.1476485,5.55];
ybs[5950]=['',4.1809896,0.9546999,4.95];
ybs[5951]=['ι1 Nor',4.2104765,-1.0091916,4.63];
ybs[5952]=['η Nor',4.2084575,-0.860046,4.65];
ybs[5953]=['',4.196306,0.0764326,5.83];
ybs[5954]=['',4.1869173,0.869735,6.05];
ybs[5955]=['',4.2051238,-0.5093458,6.03];
ybs[5956]=['5 Her',4.1976324,0.3101489,5.12];
ybs[5957]=['',4.2087418,-0.6745652,4.89];
ybs[5958]=['ρ CrB',4.1962935,0.5804164,5.41];
ybs[5959]=['',4.2080133,-0.4522595,5];
ybs[5960]=['',4.2092185,-0.5593387,6.01];
ybs[5961]=['ι CrB',4.1981533,0.5201615,4.99];
ybs[5962]=['π Ser',4.2020949,0.3971796,4.83];
ybs[5963]=['',4.210458,-0.4323784,6.21];
ybs[5964]=['',4.2124281,-0.5805202,6.1];
ybs[5965]=['',4.2139865,-0.6616516,5.9];
ybs[5966]=['43 Ser',4.2089879,0.0862107,6.08];
ybs[5967]=['ξ Sco',4.2120646,-0.1993162,5.07];
ybs[5968]=['ξ Sco',4.2120646,-0.1993162,4.77];
ybs[5969]=['',4.2272368,-0.9815223,6.16];
ybs[5970]=['δ Nor',4.2225303,-0.7892275,4.72];
ybs[5971]=['',4.1998725,0.9227202,5.93];
ybs[5972]=['υ Her',4.2033764,0.802661,4.76];
ybs[5973]=['',4.2061006,0.6385165,5.83];
ybs[5974]=['β1 Sco',4.2169644,-0.3464839,2.62];
ybs[5975]=['β2 Sco',4.2169862,-0.3464208,4.92];
ybs[5976]=['θ Dra',4.1985237,1.0213199,4.01];
ybs[5977]=['θ Lup',4.2225874,-0.6431236,4.23];
ybs[5978]=['',4.2199966,-0.4128165,5.92];
ybs[5979]=['',4.217922,-0.1106203,6.53];
ybs[5980]=['',4.2190306,-0.1079619,6.41];
ybs[5981]=['',4.2255615,-0.6423048,5.73];
ybs[5982]=['',4.2170572,0.1404928,6.29];
ybs[5983]=['ω1 Sco',4.2229696,-0.3615478,3.96];
ybs[5984]=['ι2 Nor',4.2357362,-1.0119317,5.57];
ybs[5985]=['',4.2039607,1.0360859,6.19];
ybs[5986]=['',4.2238708,-0.2463837,6.32];
ybs[5987]=['ω2 Sco',4.2255871,-0.3650248,4.32];
ybs[5988]=['',4.2277056,-0.4277372,6.33];
ybs[5989]=['',4.231322,-0.6833066,7.05];
ybs[5990]=['',4.231336,-0.6830932,6.65];
ybs[5991]=['',4.2289078,-0.460281,5.38];
ybs[5992]=['11 Sco',4.2262343,-0.2232496,5.78];
ybs[5993]=['',4.2314494,-0.414181,5.88];
ybs[5994]=['45 Ser',4.2257128,0.1718442,5.63];
ybs[5995]=['',4.2242613,0.3800746,6.14];
ybs[5996]=['',4.2352464,-0.5706298,6.19];
ybs[5997]=['',4.2368007,-0.5862675,5.54];
ybs[5998]=['κ Her',4.2274773,0.2967301,5];
ybs[5999]=['κ Her',4.2275062,0.296861,6.25];
ybs[6000]=['47 Ser',4.2294214,0.1481568,5.73];
ybs[6001]=['',4.2318033,0.0595024,5.91];
ybs[6002]=['',4.2364891,-0.3208903,6.47];
ybs[6003]=['8 Her',4.2305336,0.2995078,6.14];
ybs[6004]=['',4.2326202,0.1105449,5.97];
ybs[6005]=['',4.2433108,-0.7184472,5.86];
ybs[6006]=['',4.2357381,-0.0612928,5.37];
ybs[6007]=['',4.2417,-0.5141876,5.13];
ybs[6008]=['τ CrB',4.230737,0.6360954,4.76];
ybs[6009]=['ζ Nor',4.2532902,-0.9701286,5.81];
ybs[6010]=['δ1 Aps',4.2892329,-1.3742087,4.68];
ybs[6011]=['δ2 Aps',4.2896506,-1.3737086,5.27];
ybs[6012]=['',4.2527196,-0.9375063,5.83];
ybs[6013]=['φ Her',4.2294696,0.7834723,4.26];
ybs[6014]=['κ Nor',4.2536603,-0.9542407,4.94];
ybs[6015]=['',4.2166537,1.1827035,5.44];
ybs[6016]=['ν Sco',4.2454981,-0.3402305,6.3];
ybs[6017]=['ν Sco',4.2455784,-0.3404194,4.01];
ybs[6018]=['13 Sco',4.2471925,-0.4881738,4.59];
ybs[6019]=['12 Sco',4.2470419,-0.4967455,5.67];
ybs[6020]=['δ TrA',4.263161,-1.1122664,3.85];
ybs[6021]=['ψ Sco',4.2453334,-0.1764216,4.94];
ybs[6022]=['',4.2426014,0.1687425,6.53];
ybs[6023]=['16 Sco',4.2458231,-0.14995,5.43];
ybs[6024]=['',4.2017556,1.3394721,5.56];
ybs[6025]=['',4.2423346,0.290096,6.08];
ybs[6026]=['',4.2297829,1.0104147,6.33];
ybs[6027]=['',4.2710716,-1.1865331,5.75];
ybs[6028]=['',4.2316464,0.9736103,6.49];
ybs[6029]=['10 Her',4.2428058,0.4092884,5.7];
ybs[6030]=['',4.2642103,-1.0115007,5.63];
ybs[6031]=['',4.2492852,-0.0744301,6.25];
ybs[6032]=['',4.2534425,-0.4270001,6.41];
ybs[6033]=['',4.2425913,0.5811644,6.29];
ybs[6034]=['',4.2564034,-0.5769099,5.92];
ybs[6035]=['θ Nor',4.2609275,-0.8275474,5.14];
ybs[6036]=['',4.2430735,0.634965,5.63];
ybs[6037]=['9 Her',4.2504138,0.0868741,5.48];
ybs[6038]=['χ Sco',4.2534483,-0.2073599,5.22];
ybs[6039]=['',4.2613202,-0.7494867,6.14];
ybs[6040]=['',4.2427804,0.7388022,5.87];
ybs[6041]=['',4.256474,-0.3691475,6.41];
ybs[6042]=['',4.2475957,0.4647293,6.5];
ybs[6043]=['',4.2571459,-0.3242526,6.32];
ybs[6044]=['',4.2584161,-0.4454059,6.05];
ybs[6045]=['',4.2677481,-0.9399171,5.44];
ybs[6046]=['δ Oph',4.2553946,-0.0652336,2.74];
ybs[6047]=['',4.2546151,0.102254,6.31];
ybs[6048]=['γ1 Nor',4.2687793,-0.8745915,4.99];
ybs[6049]=['',4.270445,-0.927269,6.33];
ybs[6050]=['18 Sco',4.2610916,-0.1468192,5.5];
ybs[6051]=['',4.2623048,-0.2599099,6.09];
ybs[6052]=['',4.2774745,-1.0112622,6.49];
ybs[6053]=['σ CrB',4.2557349,0.5901922,5.64];
ybs[6054]=['σ CrB',4.2557349,0.5901874,6.66];
ybs[6055]=['16 Her',4.259721,0.327521,5.69];
ybs[6056]=['',4.2673914,-0.3725584,6.61];
ybs[6057]=['',4.2666449,-0.0697349,6.18];
ybs[6058]=['',4.2608081,0.4778637,6.14];
ybs[6059]=['',4.2432947,1.1711178,6.21];
ybs[6060]=['',4.2733847,-0.500133,4.78];
ybs[6061]=['λ Nor',4.2783156,-0.7455192,5.45];
ybs[6062]=['γ2 Nor',4.281122,-0.8760949,4.02];
ybs[6063]=['',4.284022,-0.9630856,5.77];
ybs[6064]=['υ CrB',4.2649211,0.5080302,5.78];
ybs[6065]=['ε Oph',4.2727739,-0.0826264,3.24];
ybs[6066]=['',4.2767346,-0.3535878,6.29];
ybs[6067]=['',4.2789084,-0.5401411,5.49];
ybs[6068]=['',4.2760454,-0.2602961,5.94];
ybs[6069]=['19 UMi',4.2339438,1.3235309,5.48];
ybs[6070]=['',4.283612,-0.6889089,6.12];
ybs[6071]=['ο Sco',4.283442,-0.4225472,4.55];
ybs[6072]=['20 UMi',4.2416712,1.311902,6.39];
ybs[6073]=['',4.2925547,-0.8658961,5.33];
ybs[6074]=['σ Sco',4.2858966,-0.4473852,2.89];
ybs[6075]=['',4.2922925,-0.7671106,5.88];
ybs[6076]=['',4.265366,1.0421854,5.4];
ybs[6077]=['',4.2796873,0.3681161,6.05];
ybs[6078]=['',4.2511257,1.2802285,5.98];
ybs[6079]=['',4.3064673,-1.1024111,6.15];
ybs[6080]=['',4.2746494,0.855153,5.91];
ybs[6081]=['',4.2783343,0.6923288,5.46];
ybs[6082]=['τ Her',4.2772234,0.8076016,3.89];
ybs[6083]=['σ Ser',4.2889895,0.0172606,4.82];
ybs[6084]=['',4.298787,-0.6847354,5.4];
ybs[6085]=['γ Her',4.2878106,0.3335808,3.75];
ybs[6086]=['',4.2915826,-0.0369957,6.23];
ybs[6087]=['',4.2982018,-0.5801278,6.47];
ybs[6088]=['ζ TrA',4.3212562,-1.223859,4.91];
ybs[6089]=['',4.3029432,-0.7921735,6.33];
ybs[6090]=['',4.300927,-0.6563316,5.42];
ybs[6091]=['',4.2680228,1.1957698,6.41];
ybs[6092]=['γ Aps',4.3467436,-1.3776338,3.89];
ybs[6093]=['ξ CrB',4.2881928,0.5384639,4.85];
ybs[6094]=['ψ Oph',4.2984421,-0.3504077,4.5];
ybs[6095]=['',4.3011923,-0.5191044,6.63];
ybs[6096]=['',4.3011995,-0.5191286,5.84];
ybs[6097]=['ν1 CrB',4.2892167,0.5892062,5.2];
ybs[6098]=['ν2 CrB',4.2897873,0.5875394,5.39];
ybs[6099]=['ι TrA',4.3179052,-1.1186831,5.27];
ybs[6100]=['',4.2918245,0.5636218,6.4];
ybs[6101]=['21 Her',4.2980199,0.1205789,5.85];
ybs[6102]=['ρ Oph',4.3050237,-0.4099085,5.02];
ybs[6103]=['ρ Oph',4.3050164,-0.4098891,5.92];
ybs[6104]=['',4.318553,-1.0234157,5.69];
ybs[6105]=['ε Nor',4.3130327,-0.8306569,4.47];
ybs[6106]=['η UMi',4.2630567,1.321442,4.95];
ybs[6107]=['ω Her',4.3032157,0.2442486,4.57];
ybs[6108]=['χ Oph',4.3111377,-0.3227924,4.42];
ybs[6109]=['',4.3047223,0.3290595,6.7];
ybs[6110]=['',4.325007,-1.0086769,6.06];
ybs[6111]=['',4.306672,0.1984245,6.11];
ybs[6112]=['',4.3171357,-0.6495624,5.79];
ybs[6113]=['25 Her',4.3023562,0.6519665,5.54];
ybs[6114]=['',4.30973,0.0403022,6.07];
ybs[6115]=['',4.3301161,-1.0763491,5.2];
ybs[6116]=['',4.283827,1.2054808,5.25];
ybs[6117]=['',4.2970268,0.9628215,5.74];
ybs[6118]=['',4.3138862,-0.1332742,5.23];
ybs[6119]=['υ Oph',4.3142424,-0.1467756,4.63];
ybs[6120]=['',4.2935964,1.0761177,5.67];
ybs[6121]=['',4.323959,-0.8077468,5.35];
ybs[6122]=['η Dra',4.2945261,1.0729338,2.74];
ybs[6123]=['',4.5638339,-1.5272314,6.57];
ybs[6124]=['α Sco',4.321796,-0.4619756,0.96];
ybs[6125]=['',4.3470667,-1.2395895,5.5];
ybs[6126]=['',4.3173316,0.0109489,5.39];
ybs[6127]=['',4.3186575,-0.1425315,6.48];
ybs[6128]=['',4.4066213,-1.4533188,6.57];
ybs[6129]=['',4.4848214,-1.5047181,6.04];
ybs[6130]=['',4.3230592,-0.2546088,5.68];
ybs[6131]=['22 Sco',4.3252432,-0.4389851,4.79];
ybs[6132]=['',4.3324114,-0.7304783,5.33];
ybs[6133]=['',4.3307247,-0.6063441,4.23];
ybs[6134]=['',4.3259859,-0.1318056,6.5];
ybs[6135]=['',4.3304127,-0.4638092,6.1];
ybs[6136]=['30 Her',4.3162716,0.7303149,5.04];
ybs[6137]=['φ Oph',4.329041,-0.2905872,4.28];
ybs[6138]=['β Her',4.3239447,0.3744198,2.77];
ybs[6139]=['λ Oph',4.3275346,0.033984,3.82];
ybs[6140]=['',4.3160695,0.8965771,6.29];
ybs[6141]=['θ TrA',4.352167,-1.1437135,5.52];
ybs[6142]=['',4.3254598,0.3567846,5.25];
ybs[6143]=['ω Oph',4.3335477,-0.3752917,4.45];
ybs[6144]=['',4.328299,0.3867409,5.76];
ybs[6145]=['μ Nor',4.3429534,-0.7693534,4.94];
ybs[6146]=['34 Her',4.322235,0.8538793,6.45];
ybs[6147]=['',4.3270591,0.6141509,6.25];
ybs[6148]=['28 Her',4.334774,0.0957317,5.63];
ybs[6149]=['29 Her',4.3346494,0.1998744,4.84];
ybs[6150]=['',4.347596,-0.7902804,6.46];
ybs[6151]=['15 Dra',4.3107794,1.1995644,5];
ybs[6152]=['',4.3297975,0.7952045,5.65];
ybs[6153]=['β Aps',4.3879373,-1.3534867,4.24];
ybs[6154]=['',4.352896,-0.7486312,5.47];
ybs[6155]=['τ Sco',4.3501273,-0.4930704,2.82];
ybs[6156]=['',4.3525443,-0.6159283,4.16];
ybs[6157]=['',4.3652243,-1.0650652,6.18];
ybs[6158]=['σ Her',4.340061,0.7400439,4.2];
ybs[6159]=['',4.3468462,0.2970941,6.41];
ybs[6160]=['',4.3313746,1.0609342,5.94];
ybs[6161]=['12 Oph',4.3514119,-0.0411778,5.75];
ybs[6162]=['η1 TrA',4.3773267,-1.1925575,5.91];
ybs[6163]=['',4.2969514,1.3774968,5.56];
ybs[6164]=['',4.3619249,-0.7580368,5.83];
ybs[6165]=['ζ Oph',4.355135,-0.1850308,2.56];
ybs[6166]=['',4.3524699,0.2698904,6.3];
ybs[6167]=['',4.3737362,-1.055554,6.18];
ybs[6168]=['',4.3644653,-0.6501518,5.91];
ybs[6169]=['',4.3588095,-0.1147027,6.09];
ybs[6170]=['',4.324987,1.266676,6.3];
ybs[6171]=['',4.3572498,0.2382881,6.31];
ybs[6172]=['',4.385807,-1.1774709,6.03];
ybs[6173]=['',4.3489181,0.8129498,5.79];
ybs[6174]=['16 Dra',4.3485205,0.9226774,5.53];
ybs[6175]=['17 Dra',4.3486785,0.9230994,5.08];
ybs[6176]=['17 Dra',4.3487076,0.9230947,6.53];
ybs[6177]=['',4.3749059,-0.8516435,5.65];
ybs[6178]=['',4.3764123,-0.8671504,5.65];
ybs[6179]=['',4.3659862,-0.1673372,6.35];
ybs[6180]=['',4.3703416,-0.3567715,6.26];
ybs[6181]=['',4.3193769,1.3510505,6.34];
ybs[6182]=['',4.3759568,-0.5790789,5.87];
ybs[6183]=['',4.3749626,-0.4276148,6.09];
ybs[6184]=['36 Her',4.3696719,0.0728554,6.93];
ybs[6185]=['37 Her',4.3699334,0.073074,5.77];
ybs[6186]=['',4.374615,-0.3102273,4.96];
ybs[6187]=['',4.3822431,-0.8046386,6.23];
ybs[6188]=['',4.3506548,1.1002252,6.16];
ybs[6189]=['',4.3561772,0.9770617,5.29];
ybs[6190]=['42 Her',4.3599457,0.8533719,4.9];
ybs[6191]=['',4.3724647,-0.0180331,6.24];
ybs[6192]=['',4.376087,-0.3483119,5.57];
ybs[6193]=['',4.3706231,0.2157608,6.08];
ybs[6194]=['',4.4001562,-1.1718142,5.13];
ybs[6195]=['14 Oph',4.374658,0.0200478,5.74];
ybs[6196]=['',4.3850406,-0.718211,6.2];
ybs[6197]=['',4.3897353,-0.9282302,5.96];
ybs[6198]=['',4.3709022,0.4332926,6.06];
ybs[6199]=['',4.3856514,-0.7181082,6.12];
ybs[6200]=['',4.3850529,-0.6665104,6.05];
ybs[6201]=['',4.3841507,-0.5609094,6.46];
ybs[6202]=['ζ Her',4.3718658,0.5510071,2.81];
ybs[6203]=['39 Her',4.3734557,0.4692214,5.92];
ybs[6204]=['',4.3891812,-0.713332,5.71];
ybs[6205]=['',4.3975959,-1.0216124,5.74];
ybs[6206]=['',4.386785,-0.4797473,6.58];
ybs[6207]=['α TrA',4.4092766,-1.2052762,1.92];
ybs[6208]=['',4.3899437,-0.4981312,6.02];
ybs[6209]=['',4.4018226,-1.0187742,5.58];
ybs[6210]=['η Her',4.3785797,0.6787617,3.53];
ybs[6211]=['',4.3982196,-0.6877919,5.48];
ybs[6212]=['',4.3830006,0.5935384,5.99];
ybs[6213]=['18 Dra',4.3678774,1.1267189,4.83];
ybs[6214]=['16 Oph',4.3911848,0.0172672,6.03];
ybs[6215]=['25 Sco',4.3979194,-0.4460881,6.71];
ybs[6216]=['',4.3778586,0.9714191,6.16];
ybs[6217]=['',4.3902355,0.2742657,5.56];
ybs[6218]=['43 Her',4.3924391,0.1492551,5.15];
ybs[6219]=['η Ara',4.4126387,-1.0309741,3.76];
ybs[6220]=['',4.3883909,0.7537395,6.05];
ybs[6221]=['',4.4248323,-1.1817615,6.32];
ybs[6222]=['19 Oph',4.3984348,0.0355031,6.1];
ybs[6223]=['',4.4227066,-1.1415105,6.13];
ybs[6224]=['45 Her',4.4010061,0.0910476,5.24];
ybs[6225]=['',4.4045383,-0.2607377,6.03];
ybs[6226]=['',4.4154208,-0.873962,6.47];
ybs[6227]=['',4.3878885,0.9904886,4.85];
ybs[6228]=['',4.3498824,1.3767851,6.32];
ybs[6229]=['',4.4023996,0.2366783,6.35];
ybs[6230]=['',4.4089835,-0.2739611,6.1];
ybs[6231]=['ε Sco',4.4126902,-0.5990377,2.29];
ybs[6232]=['',4.397753,0.7366795,5.87];
ybs[6233]=['20 Oph',4.410452,-0.1887089,4.65];
ybs[6234]=['',4.416467,-0.6552507,6.11];
ybs[6235]=['',4.4191101,-0.7201048,5.22];
ybs[6236]=['',4.4086409,0.2309434,5.91];
ybs[6237]=['μ1 Sco',4.4203034,-0.664548,3.08];
ybs[6238]=['',4.4125673,-0.0468246,6.32];
ybs[6239]=['',4.4224448,-0.7309884,6.49];
ybs[6240]=['47 Her',4.4120747,0.1259916,5.49];
ybs[6241]=['',4.4309635,-1.0111881,5.94];
ybs[6242]=['μ2 Sco',4.4223317,-0.6640211,3.57];
ybs[6243]=['',4.4377312,-1.104732,6.02];
ybs[6244]=['52 Her',4.4058629,0.8020457,4.82];
ybs[6245]=['21 Oph',4.4170099,0.0207269,5.51];
ybs[6246]=['',4.4079175,0.7574945,6.13];
ybs[6247]=['',4.428575,-0.7518594,5.96];
ybs[6248]=['50 Her',4.4127763,0.5197202,5.72];
ybs[6249]=['',4.412968,0.5676638,6.13];
ybs[6250]=['',4.4299141,-0.7301375,5.45];
ybs[6251]=['',4.4297055,-0.7334249,6.32];
ybs[6252]=['ζ1 Sco',4.42979,-0.7398389,4.73];
ybs[6253]=['',4.4306439,-0.7309023,6.45];
ybs[6254]=['',4.4120392,0.7307297,6.29];
ybs[6255]=['',4.4312098,-0.7303729,6.59];
ybs[6256]=['',4.4317751,-0.7418718,5.88];
ybs[6257]=['',4.3735258,1.3523148,5.98];
ybs[6258]=['49 Her',4.419518,0.2608547,6.52];
ybs[6259]=['',4.4263997,-0.3568021,5.88];
ybs[6260]=['51 Her',4.4177814,0.4298389,5.04];
ybs[6261]=['ζ2 Sco',4.4323582,-0.7398201,3.62];
ybs[6262]=['',4.4339985,-0.7186941,5.77];
ybs[6263]=['',4.4319068,-0.5343224,6.35];
ybs[6264]=['',4.4396953,-0.8849084,6.33];
ybs[6265]=['',4.4412546,-0.9129864,5.94];
ybs[6266]=['',4.4569427,-1.209397,5.79];
ybs[6267]=['',4.4291403,-0.0286171,6.25];
ybs[6268]=['',4.4315952,-0.2062927,6.57];
ybs[6269]=['53 Her',4.4228161,0.5528106,5.32];
ybs[6270]=['23 Oph',4.4310957,-0.107881,5.25];
ybs[6271]=['ι Oph',4.4280674,0.1769377,4.38];
ybs[6272]=['',4.4379733,-0.5852763,6.37];
ybs[6273]=['',4.4410827,-0.7129665,6.15];
ybs[6274]=['',4.4376751,-0.2937869,6.37];
ybs[6275]=['ζ Ara',4.4509152,-0.9776595,3.13];
ybs[6276]=['',4.4234675,0.8270956,6];
ybs[6277]=['',4.4317123,0.3653231,5.41];
ybs[6278]=['27 Sco',4.4432952,-0.580943,5.48];
ybs[6279]=['',4.4490864,-0.8843015,5.55];
ybs[6280]=['',4.4334545,0.2372379,6.34];
ybs[6281]=['24 Oph',4.441241,-0.4045029,5.58];
ybs[6282]=['56 Her',4.4320525,0.4486099,6.08];
ybs[6283]=['54 Her',4.4337567,0.321252,5.35];
ybs[6284]=['',4.4422857,-0.3414947,6.27];
ybs[6285]=['ε1 Ara',4.4548865,-0.9282648,4.06];
ybs[6286]=['',4.443629,-0.1918013,6.19];
ybs[6287]=['',4.4572682,-0.9533306,5.65];
ybs[6288]=['',4.4508528,-0.6570562,6.09];
ybs[6289]=['κ Oph',4.4440591,0.1631708,3.2];
ybs[6290]=['',4.4583433,-0.849496,6];
ybs[6291]=['',4.4433356,0.2418695,6.37];
ybs[6292]=['',4.4492378,-0.2599715,6.59];
ybs[6293]=['',4.4585236,-0.7937129,6.65];
ybs[6294]=['',4.4650518,-1.0294388,6.11];
ybs[6295]=['57 Her',4.4428917,0.4420341,6.28];
ybs[6296]=['',4.4355646,0.8728771,6.56];
ybs[6297]=['',4.443747,0.4250816,6.32];
ybs[6298]=['',4.4550986,-0.4383735,5.86];
ybs[6299]=['26 Oph',4.4559606,-0.4365783,5.75];
ybs[6300]=['',4.4583792,-0.627601,5.97];
ybs[6301]=['',4.4642633,-0.8928239,6.45];
ybs[6302]=['',4.4435621,0.7415295,6.34];
ybs[6303]=['ε2 Ara',4.4704421,-0.9295725,5.29];
ybs[6304]=['19 Dra',4.4336483,1.1363472,4.89];
ybs[6305]=['',4.4637344,-0.5614344,5.03];
ybs[6306]=['',4.4564517,0.1144722,6.59];
ybs[6307]=['30 Oph',4.45925,-0.0741257,4.82];
ybs[6308]=['20 Dra',4.4353709,1.1346822,6.41];
ybs[6309]=['',4.4763682,-1.0076713,5.73];
ybs[6310]=['29 Oph',4.4631536,-0.3300382,6.26];
ybs[6311]=['ε UMi',4.3817881,1.4312737,4.23];
ybs[6312]=['',4.4724063,-0.8235062,6.06];
ybs[6313]=['ε Her',4.4547842,0.5393317,3.92];
ybs[6314]=['',4.4580428,0.3945763,5.65];
ybs[6315]=['',4.460827,0.2604909,6.31];
ybs[6316]=['',4.4725846,-0.6662901,5.91];
ybs[6317]=['',4.4587195,0.4742372,6.55];
ybs[6318]=['',4.46292,0.1470672,6.33];
ybs[6319]=['',4.4492065,0.9889594,6.03];
ybs[6320]=['',4.478388,-0.7945528,6.28];
ybs[6321]=['59 Her',4.4604216,0.5854516,5.25];
ybs[6322]=['',4.463804,0.444735,5.75];
ybs[6323]=['',4.4766691,-0.5959561,4.87];
ybs[6324]=['',4.4328859,1.2758526,6.3];
ybs[6325]=['',4.4634563,0.5560721,6.36];
ybs[6326]=['',4.4677467,0.2455361,4.98];
ybs[6327]=['',4.481458,-0.7701713,6.19];
ybs[6328]=['',4.4679229,0.2528523,6.52];
ybs[6329]=['',4.4758686,-0.3581026,6.3];
ybs[6330]=['',4.470052,0.2370459,5.93];
ybs[6331]=['',4.4714129,0.2363888,6.08];
ybs[6332]=['',4.4708353,0.3432555,6.35];
ybs[6333]=['',4.4834173,-0.6501326,5.98];
ybs[6334]=['',4.446016,1.207083,6.4];
ybs[6335]=['61 Her',4.4686169,0.6176811,6.69];
ybs[6336]=['',4.4839274,-0.6191279,6.13];
ybs[6337]=['',4.4571569,1.0580975,6.13];
ybs[6338]=['',4.4775288,0.011862,6.01];
ybs[6339]=['',4.4821955,-0.3767672,6.3];
ybs[6340]=['',4.4703376,0.606795,6.04];
ybs[6341]=['',4.474372,0.3416663,6.17];
ybs[6342]=['',4.4786947,-0.0159644,5.64];
ybs[6343]=['',4.4853783,-0.4631267,6.29];
ybs[6344]=['60 Her',4.4775971,0.221971,4.91];
ybs[6345]=['',4.5016786,-1.076803,6.39];
ybs[6346]=['',4.5130665,-1.2346601,6.22];
ybs[6347]=['',4.481105,0.1694859,6.37];
ybs[6348]=['',4.4813308,0.1820672,6.37];
ybs[6349]=['',4.460876,1.1270688,6.1];
ybs[6350]=['',4.4845782,-0.0292967,6.38];
ybs[6351]=['',4.4750987,0.7642659,6.43];
ybs[6352]=['',4.4736915,0.8513897,6.09];
ybs[6353]=['',4.4813283,0.3850493,5.56];
ybs[6354]=['',4.4910198,-0.3077149,5.99];
ybs[6355]=['',4.4938379,-0.5310158,5.97];
ybs[6356]=['',4.4904303,-0.0192175,6.06];
ybs[6357]=['',4.5164648,-1.1731407,5.89];
ybs[6358]=['μ Dra',4.4754463,0.9502852,5.83];
ybs[6359]=['μ Dra',4.475439,0.9502852,5.8];
ybs[6360]=['',4.5028499,-0.7780334,5.08];
ybs[6361]=['',4.4934859,-0.0681399,6.36];
ybs[6362]=['',4.53295,-1.301159,6.25];
ybs[6363]=['',4.5072347,-0.8533664,5.84];
ybs[6364]=['',4.4975725,-0.1840328,5.56];
ybs[6365]=['',4.4870271,0.7067575,6.34];
ybs[6366]=['',4.4883543,0.6268089,5.39];
ybs[6367]=['η Oph',4.5002552,-0.2748099,2.43];
ybs[6368]=['',4.4555256,1.3137543,6.21];
ybs[6369]=['η Sco',4.5090977,-0.7550137,3.33];
ybs[6370]=['',4.5094141,-0.6898735,5.67];
ybs[6371]=['',4.5094028,-0.6779228,6.3];
ybs[6372]=['',4.4886184,0.8869853,6.46];
ybs[6373]=['',4.5191076,-0.9932207,6.09];
ybs[6374]=['',4.5010999,0.2172345,6.57];
ybs[6375]=['',4.5086375,-0.4411261,6.54];
ybs[6376]=['',4.5095577,-0.4848838,6.14];
ybs[6377]=['',4.4947635,0.7113276,5.08];
ybs[6378]=['',4.5121908,-0.5665029,6.01];
ybs[6379]=['',4.5055563,0.1374366,6.33];
ybs[6380]=['63 Her',4.5019843,0.4226712,6.19];
ybs[6381]=['',4.5189925,-0.6943953,6.6];
ybs[6382]=['37 Oph',4.5085738,0.1844007,5.33];
ybs[6383]=['',4.5108091,0.0057989,6.65];
ybs[6384]=['',4.4982045,0.9143448,6.29];
ybs[6385]=['ζ Dra',4.489169,1.1465617,3.17];
ybs[6386]=['',4.5224636,-0.5858541,5.53];
ybs[6387]=['',4.5238966,-0.6739134,5.96];
ybs[6388]=['',4.5034595,0.8678888,6.04];
ybs[6389]=['',4.547319,-1.2228122,6.53];
ybs[6390]=['36 Oph',4.5223163,-0.4646313,5.11];
ybs[6391]=['36 Oph',4.5223018,-0.4646071,5.07];
ybs[6392]=['',4.5246707,-0.5275951,6.21];
ybs[6393]=['',4.5218637,-0.2548625,5.99];
ybs[6394]=['',4.5270828,-0.6242631,6.12];
ybs[6395]=['α1 Her',4.5179909,0.2508263,3.48];
ybs[6396]=['α2 Her',4.5180127,0.2508215,5.39];
ybs[6397]=['',4.5411741,-1.04216,5.91];
ybs[6398]=['',4.5300206,-0.5703856,5.55];
ybs[6399]=['δ Her',4.5193192,0.4331963,3.14];
ybs[6400]=['ι Aps',4.5555698,-1.2241558,5.41];
ybs[6401]=['',4.5252854,0.0378353,6.17];
ybs[6402]=['',4.5276119,-0.1093116,6.09];
ybs[6403]=['',4.5265866,0.0208108,5.88];
ybs[6404]=['41 Oph',4.5269985,-0.0080883,4.73];
ybs[6405]=['',4.539412,-0.8142123,5.48];
ybs[6406]=['ζ Aps',4.5545605,-1.1830933,4.78];
ybs[6407]=['π Her',4.5189086,0.642112,3.16];
ybs[6408]=['',4.5222443,0.4140656,5.96];
ybs[6409]=['',4.5381579,-0.7705084,5.76];
ybs[6410]=['',4.5059865,1.0970171,5.56];
ybs[6411]=['',4.5355949,-0.568466,6.36];
ybs[6412]=['',4.5416108,-0.8740639,6.27];
ybs[6413]=['ο Oph',4.5338489,-0.424193,5.2];
ybs[6414]=['ο Oph',4.5338343,-0.4241446,6.8];
ybs[6415]=['',4.5383747,-0.6109846,5.91];
ybs[6416]=['',4.5408182,-0.7721329,6.65];
ybs[6417]=['',4.5349291,-0.285001,6.43];
ybs[6418]=['',4.6022534,-1.4114578,5.88];
ybs[6419]=['',4.5305778,0.402701,6.45];
ybs[6420]=['68 Her',4.5290085,0.5773916,4.82];
ybs[6421]=['',4.5328853,0.3019507,6];
ybs[6422]=['',4.5354127,0.1893179,5.03];
ybs[6423]=['',4.5367044,0.1059078,6.51];
ybs[6424]=['',4.541819,-0.3101998,6.02];
ybs[6425]=['69 Her',4.5303374,0.6505523,4.65];
ybs[6426]=['',4.5258819,0.8669566,7.48];
ybs[6427]=['',4.5571849,-1.0127387,5.88];
ybs[6428]=['',4.5419071,-0.1035669,6.32];
ybs[6429]=['',4.5625827,-1.0974465,5.7];
ybs[6430]=['',4.5448449,-0.3377079,6.52];
ybs[6431]=['',4.5579017,-0.9868193,5.8];
ybs[6432]=['',4.5356493,0.5027557,5.65];
ybs[6433]=['',4.5333921,0.6770816,5.94];
ybs[6434]=['ξ Oph',4.546795,-0.3687716,4.39];
ybs[6435]=['ν Ser',4.5457663,-0.2245071,4.33];
ybs[6436]=['',4.5635756,-1.0592122,5.77];
ybs[6437]=['',4.523483,1.0585809,6.32];
ybs[6438]=['',4.5459189,-0.1869677,6.46];
ybs[6439]=['',4.5546642,-0.6600932,6.41];
ybs[6440]=['ι Ara',4.5578615,-0.8287452,5.25];
ybs[6441]=['',4.542598,0.3148678,5];
ybs[6442]=['θ Oph',4.5513224,-0.4365992,3.27];
ybs[6443]=['',4.5544755,-0.6270194,6.47];
ybs[6444]=['',4.5416821,0.445422,5.38];
ybs[6445]=['',4.5557643,-0.6498908,5.93];
ybs[6446]=['70 Her',4.5449475,0.42731,5.12];
ybs[6447]=['72 Her',4.5435833,0.5663815,5.39];
ybs[6448]=['43 Oph',4.5573289,-0.4914557,5.35];
ybs[6449]=['',4.5618007,-0.7710408,5.12];
ybs[6450]=['β Ara',4.5673456,-0.9694324,2.85];
ybs[6451]=['γ Ara',4.5678324,-0.9842234,3.34];
ybs[6452]=['',4.5480579,0.2917271,6.35];
ybs[6453]=['74 Her',4.541555,0.8067639,5.59];
ybs[6454]=['',4.5542971,-0.0419554,6.29];
ybs[6455]=['',4.5474818,0.5016411,6.35];
ybs[6456]=['',4.5423402,0.8407556,6.43];
ybs[6457]=['κ Ara',4.5700008,-0.8839695,5.23];
ybs[6458]=['',4.547898,0.6974051,5.51];
ybs[6459]=['',4.5649469,-0.6058203,6.16];
ybs[6460]=['',4.580549,-1.1004221,6.24];
ybs[6461]=['',4.5629376,-0.3744847,5.85];
ybs[6462]=['',4.5624817,-0.3221982,6.21];
ybs[6463]=['',4.5647944,-0.4233848,6.19];
ybs[6464]=['',4.5741951,-0.9069233,6.19];
ybs[6465]=['',4.5587869,0.1542469,5.77];
ybs[6466]=['',4.5734288,-0.8003623,5.29];
ybs[6467]=['',4.5752672,-0.8839025,5.92];
ybs[6468]=['',4.5471753,0.9320836,5.67];
ybs[6469]=['73 Her',4.5589913,0.4004699,5.74];
ybs[6470]=['',4.5610224,0.2842488,5.71];
ybs[6471]=['',4.5612117,0.272119,6.35];
ybs[6472]=['',4.5786815,-0.9129905,5.75];
ybs[6473]=['ρ Her',4.5565473,0.6480655,5.47];
ybs[6474]=['ρ Her',4.5565692,0.648051,4.52];
ybs[6475]=['44 Oph',4.5703195,-0.4221829,4.17];
ybs[6476]=['',4.5819132,-0.9631199,5.94];
ybs[6477]=['',4.558043,0.6731329,6.49];
ybs[6478]=['',4.5678455,-0.0290753,6.44];
ybs[6479]=['',4.572783,-0.4530372,6.44];
ybs[6480]=['',4.5599299,0.6446726,6.28];
ybs[6481]=['45 Oph',4.574829,-0.5215138,4.29];
ybs[6482]=['',4.5708557,-0.0890224,4.54];
ybs[6483]=['',4.5760017,-0.5190247,6];
ybs[6484]=['',4.5670309,0.2950169,5.98];
ybs[6485]=['',4.5728393,-0.2186244,6.21];
ybs[6486]=['',4.5691071,0.1323216,6.06];
ybs[6487]=['σ Oph',4.570075,0.0720171,4.34];
ybs[6488]=['',4.5671625,0.4688764,6.41];
ybs[6489]=['δ Ara',4.5932014,-1.0593417,3.62];
ybs[6490]=['',4.582016,-0.6421284,6.02];
ybs[6491]=['',4.5709051,0.3502341,5.54];
ybs[6492]=['',4.5842506,-0.6724644,6.39];
ybs[6493]=['',4.5770866,-0.1434953,6.37];
ybs[6494]=['',4.5940305,-0.9936574,5.95];
ybs[6495]=['',4.5701445,0.6053131,5.94];
ybs[6496]=['',4.5802809,0.005542,5.44];
ybs[6497]=['υ Sco',4.5900189,-0.6511471,2.69];
ybs[6498]=['77 Her',4.5692782,0.8420517,5.85];
ybs[6499]=['α Ara',4.5954422,-0.8707058,2.95];
ybs[6500]=['',4.5636723,1.0477883,5.65];
ybs[6501]=['',4.5846615,-0.1035387,6.37];
ybs[6502]=['',4.5950876,-0.8036904,6.03];
ybs[6503]=['',4.5655454,1.0234197,6.51];
ybs[6504]=['',4.5871651,-0.0187598,5.31];
ybs[6505]=['',4.5943393,-0.5884289,6.44];
ybs[6506]=['',4.5595832,1.1744594,6.43];
ybs[6507]=['51 Oph',4.5923345,-0.418437,4.81];
ybs[6508]=['',4.5938284,-0.4586982,6.05];
ybs[6509]=['',4.5866758,0.2079143,6.39];
ybs[6510]=['',4.5970624,-0.5984939,6.17];
ybs[6511]=['',4.6005087,-0.7188095,5.84];
ybs[6512]=['',4.591235,0.0473418,5.59];
ybs[6513]=['',4.6126695,-1.0446868,6.28];
ybs[6514]=['λ Her',4.5877978,0.4555013,4.41];
ybs[6515]=['λ Sco',4.6024254,-0.647776,1.63];
ybs[6516]=['',4.5884118,0.5436028,5.61];
ybs[6517]=['',4.5304821,1.3983409,5.72];
ybs[6518]=['',4.6109475,-0.931359,6.1];
ybs[6519]=['',4.5869736,0.6784079,6.43];
ybs[6520]=['',4.5948638,0.2080202,6.42];
ybs[6521]=['78 Her',4.5924604,0.4955984,5.62];
ybs[6522]=['',4.6008372,-0.100457,5.62];
ybs[6523]=['',4.6070274,-0.5688405,5.7];
ybs[6524]=['β Dra',4.5851119,0.9126141,2.79];
ybs[6525]=['σ Ara',4.6118822,-0.8118505,4.59];
ybs[6526]=['',4.5930594,0.5979341,6.56];
ybs[6527]=['',4.6116418,-0.6536267,6.48];
ybs[6528]=['',4.5859022,1.009918,6.4];
ybs[6529]=['',4.5995608,0.3358979,5.64];
ybs[6530]=['',4.6008636,0.284602,5.69];
ybs[6531]=['',4.6011578,0.2588442,6.48];
ybs[6532]=['',4.6065666,-0.1963921,5.55];
ybs[6533]=['52 Oph',4.6092497,-0.3849174,6.57];
ybs[6534]=['',4.6153207,-0.6744821,4.29];
ybs[6535]=['',4.619946,-0.8738738,5.93];
ybs[6536]=['53 Oph',4.6052397,0.167134,5.81];
ybs[6537]=['π Ara',4.6230908,-0.9513663,5.25];
ybs[6538]=['',4.597543,0.7196399,5.74];
ybs[6539]=['',4.6043335,0.2878609,6.4];
ybs[6540]=['',4.743451,-1.4820612,6.45];
ybs[6541]=['θ Sco',4.6189208,-0.7506162,1.87];
ybs[6542]=['ν1 Dra',4.5924734,0.9629407,4.88];
ybs[6543]=['ν2 Dra',4.592867,0.9627475,4.87];
ybs[6544]=['α Oph',4.6065593,0.2190308,2.08];
ybs[6545]=['',4.6192208,-0.6645319,6.26];
ybs[6546]=['',4.6225001,-0.7485642,6.1];
ybs[6547]=['',4.6109046,0.3662763,6.1];
ybs[6548]=['',4.598141,1.0043972,6.17];
ybs[6549]=['ξ Ser',4.6189723,-0.2689188,3.54];
ybs[6550]=['',4.6190504,-0.2719294,5.94];
ybs[6551]=['',4.6090147,0.6508592,6.1];
ybs[6552]=['',4.611251,0.491742,6.38];
ybs[6553]=['',4.653215,-1.2605997,6.49];
ybs[6554]=['27 Dra',4.5897723,1.1889724,5.05];
ybs[6555]=['μ Oph',4.6198747,-0.1418619,4.62];
ybs[6556]=['',4.6213272,-0.1908597,5.75];
ybs[6557]=['λ Ara',4.6327309,-0.8626044,4.77];
ybs[6558]=['',4.6132439,0.5371337,6.02];
ybs[6559]=['79 Her',4.6174433,0.4241256,5.77];
ybs[6560]=['',4.6363978,-0.8190763,5.79];
ybs[6561]=['26 Dra',4.6040212,1.079738,5.23];
ybs[6562]=['82 Her',4.6123995,0.8478065,5.37];
ybs[6563]=['',4.625231,0.0352449,6.26];
ybs[6564]=['',4.6401043,-0.8817034,6.24];
ybs[6565]=['',4.6241116,0.2324849,6.12];
ybs[6566]=['',4.6299579,-0.0377116,6.19];
ybs[6567]=['',4.6228351,0.5712566,6.37];
ybs[6568]=['κ Sco',4.6412683,-0.6813273,2.41];
ybs[6569]=['ο Ser',4.6355977,-0.2248501,4.26];
ybs[6570]=['η Pav',4.6576479,-1.1297444,3.62];
ybs[6571]=['',4.6427472,-0.6449491,5.54];
ybs[6572]=['',4.6278283,0.5444403,6.03];
ybs[6573]=['μ Ara',4.6492811,-0.9047893,5.15];
ybs[6574]=['',4.6532161,-1.0044607,6.01];
ybs[6575]=['',4.6437316,-0.5769716,6.4];
ybs[6576]=['ι Her',4.6249475,0.8028127,3.8];
ybs[6577]=['',4.6337396,0.2647705,6.34];
ybs[6578]=['',4.6355609,0.1100451,5.95];
ybs[6579]=['',4.6310025,0.5459292,6.28];
ybs[6580]=['',4.6330266,0.4277009,6.36];
ybs[6581]=['',4.644318,-0.4867902,6.36];
ybs[6582]=['',4.6371836,0.2782833,5.52];
ybs[6583]=['58 Oph',4.6446768,-0.3785645,4.87];
ybs[6584]=['ω Dra',4.611384,1.1998834,4.8];
ybs[6585]=['',4.6511113,-0.7458686,5.87];
ybs[6586]=['',4.6099296,1.2140667,6.42];
ybs[6587]=['',4.6301642,0.7585676,6.59];
ybs[6588]=['',4.646068,-0.2358861,6.39];
ybs[6589]=['',4.6457547,-0.1236763,6.3];
ybs[6590]=['83 Her',4.6390531,0.4285983,5.52];
ybs[6591]=['β Oph',4.6440559,0.0795938,2.77];
ybs[6592]=['',4.6432938,0.2493744,6.24];
ybs[6593]=['',4.6290286,1.0001103,6.77];
ybs[6594]=['',4.611236,1.2644224,5.86];
ybs[6595]=['',4.6328297,0.9042591,5.99];
ybs[6596]=['84 Her',4.6429301,0.4244791,5.71];
ybs[6597]=['61 Oph',4.6488856,0.0449086,6.17];
ybs[6598]=['',4.6489874,0.0448991,6.56];
ybs[6599]=['',4.6473116,0.2513933,6.19];
ybs[6600]=['',4.6408908,0.7692953,6.34];
ybs[6601]=['',4.6614423,-0.6652651,6.43];
ybs[6602]=['',4.6691759,-0.9670204,6.11];
ybs[6603]=['ι1 Sco',4.6635633,-0.7004348,3.03];
ybs[6604]=['3 Sgr',4.6629193,-0.4858279,4.54];
ybs[6605]=['',4.6635982,-0.392403,6.18];
ybs[6606]=['',4.6441207,0.9388989,5.75];
ybs[6607]=['',4.6527431,0.5497573,6.23];
ybs[6608]=['',4.6627019,-0.2571025,5.94];
ybs[6609]=['',4.6668289,-0.4708841,6.35];
ybs[6610]=['',4.6770259,-0.9357752,5.92];
ybs[6611]=['μ Her',4.6563239,0.4837171,3.42];
ybs[6612]=['',4.6825776,-1.0501241,5.78];
ybs[6613]=['',4.6533893,0.6785059,6.52];
ybs[6614]=['',4.6537166,0.6862052,6.68];
ybs[6615]=['',4.6596194,0.3087822,5.72];
ybs[6616]=['',4.6701189,-0.5534037,4.83];
ybs[6617]=['γ Oph',4.663397,0.0471634,3.75];
ybs[6618]=['',4.6733335,-0.6465989,3.21];
ybs[6619]=['ι2 Sco',4.6749141,-0.6997804,4.81];
ybs[6620]=['',4.6800985,-0.9273632,6.09];
ybs[6621]=['',4.6652987,0.0663119,6.22];
ybs[6622]=['',4.6908257,-1.1430448,6.49];
ybs[6623]=['',4.7131636,-1.3295561,6.07];
ybs[6624]=['ψ1 Dra',4.6322314,1.2591007,4.58];
ybs[6625]=['ψ1 Dra',4.6323524,1.2592415,5.79];
ybs[6626]=['',4.665109,0.3588536,5.69];
ybs[6627]=['',4.6696223,0.0341518,6.47];
ybs[6628]=['',4.6820043,-0.7959363,6.11];
ybs[6629]=['',4.6583083,0.8308963,6.43];
ybs[6630]=['',4.6668328,0.3359877,6.12];
ybs[6631]=['',4.6808835,-0.7116725,5.96];
ybs[6632]=['87 Her',4.666706,0.4471215,5.12];
ybs[6633]=['',4.678945,-0.5333851,6.66];
ybs[6634]=['',4.7514127,-1.4221524,6.35];
ybs[6635]=['',4.6835662,-0.6074133,5.9];
ybs[6636]=['',4.6839931,-0.6007368,5.84];
ybs[6637]=['',4.6867587,-0.7330283,6.2];
ybs[6638]=['',4.6754733,0.2084429,6.17];
ybs[6639]=['',4.6861253,-0.5954535,6.06];
ybs[6640]=['',4.6866522,-0.6112382,6.45];
ybs[6641]=['',4.6868246,-0.6218068,6.03];
ybs[6642]=['',4.6733804,0.5117005,5.5];
ybs[6643]=['',4.6754916,0.389429,5.98];
ybs[6644]=['30 Dra',4.6665375,0.8862181,5.02];
ybs[6645]=['',4.688363,-0.6062077,6.17];
ybs[6646]=['',4.6886393,-0.6090822,5.6];
ybs[6647]=['',4.6813962,-0.0216399,6.35];
ybs[6648]=['',4.690249,-0.6071693,6.38];
ybs[6649]=['',4.6844005,-0.1072773,6.21];
ybs[6650]=['',4.6909312,-0.6065863,5.96];
ybs[6651]=['',4.6911674,-0.6079677,6.42];
ybs[6652]=['88 Her',4.6709962,0.8445654,6.68];
ybs[6653]=['',4.6808211,0.2674297,6.46];
ybs[6654]=['',4.6863484,-0.1902839,6.18];
ybs[6655]=['',4.6839316,0.0227248,5.95];
ybs[6656]=['',4.6932676,-0.6015888,5.96];
ybs[6657]=['',4.6766213,0.6993342,6.46];
ybs[6658]=['',4.6866103,0.1064423,5.77];
ybs[6659]=['',4.6963174,-0.6366551,6.06];
ybs[6660]=['',4.6948382,-0.4343978,6.2];
ybs[6661]=['',4.6803126,0.6977648,6.04];
ybs[6662]=['',4.6796688,0.8140223,6.38];
ybs[6663]=['',4.703958,-0.7739369,4.86];
ybs[6664]=['',4.6908208,0.1942251,6.38];
ybs[6665]=['90 Her',4.6856489,0.6982247,5.16];
ybs[6666]=['',4.7043494,-0.7034829,6.43];
ybs[6667]=['',4.6990645,-0.3281872,6.52];
ybs[6668]=['',4.7028028,-0.4898519,5.8];
ybs[6669]=['',4.7007125,-0.2760038,5.89];
ybs[6670]=['',4.7081973,-0.7280952,4.88];
ybs[6671]=['',4.7087995,-0.6830791,6.29];
ybs[6672]=['',4.7001674,0.0116745,5.82];
ybs[6673]=['89 Her',4.6954875,0.4546269,5.46];
ybs[6674]=['',4.7024394,-0.0712637,5.47];
ybs[6675]=['',4.6974771,0.3920456,5.58];
ybs[6676]=['ξ Dra',4.6854599,0.9925703,3.75];
ybs[6677]=['',4.7035233,0.0011403,5.97];
ybs[6678]=['',4.7027283,0.1132135,6.29];
ybs[6679]=['',4.712894,-0.6433029,5.74];
ybs[6680]=['',4.7113526,-0.501948,6.01];
ybs[6681]=['',4.7133091,-0.528018,5.16];
ybs[6682]=['',4.7133309,-0.5280229,7.04];
ybs[6683]=['θ Her',4.6986617,0.6501192,3.86];
ybs[6684]=['',4.7048368,0.1927412,6.36];
ybs[6685]=['',4.7034955,0.4187885,6.3];
ybs[6686]=['ν Oph',4.7123484,-0.1705853,3.34];
ybs[6687]=['',4.6936917,0.976852,6.1];
ybs[6688]=['4 Sgr',4.7161513,-0.4156672,4.76];
ybs[6689]=['35 Dra',4.6629355,1.3431727,5.04];
ybs[6690]=['',4.7006477,0.7914996,6.02];
ybs[6691]=['ξ Her',4.7055984,0.5104559,3.7];
ybs[6692]=['',4.7169401,-0.3549815,6.21];
ybs[6693]=['γ Dra',4.6993597,0.8986269,2.23];
ybs[6694]=['',4.7147515,-0.0841485,5.87];
ybs[6695]=['ν Her',4.7087828,0.5268966,4.41];
ybs[6696]=['',4.7254309,-0.6348942,6.3];
ybs[6697]=['',4.717411,0.0109911,6.37];
ybs[6698]=['ζ Ser',4.7185211,-0.0644005,4.62];
ybs[6699]=['',4.7094035,0.6333339,6];
ybs[6700]=['66 Oph',4.7173216,0.0762518,4.64];
ybs[6701]=['93 Her',4.7160379,0.2923605,4.67];
ybs[6702]=['67 Oph',4.7190297,0.0511752,3.97];
ybs[6703]=['6 Sgr',4.7228707,-0.2994312,6.28];
ybs[6704]=['',4.725336,-0.3975776,5.77];
ybs[6705]=['',4.6652268,1.3666327,6.24];
ybs[6706]=['',4.7096503,0.7937014,6.48];
ybs[6707]=['',4.7199627,0.1094128,6.34];
ybs[6708]=['',4.7177146,0.3404472,6.5];
ybs[6709]=['χ Oct',4.9944276,-1.5272378,5.28];
ybs[6710]=['',4.7200036,0.2634384,6.26];
ybs[6711]=['68 Oph',4.7239141,0.0227977,4.45];
ybs[6712]=['7 Sgr',4.7295122,-0.4237796,5.34];
ybs[6713]=['ψ2 Dra',4.690067,1.2566875,5.45];
ybs[6714]=['',4.7178363,0.5796986,5.99];
ybs[6715]=['',4.7302279,-0.3964833,6.74];
ybs[6716]=['',4.7142956,0.7941505,5.67];
ybs[6717]=['95 Her',4.7221671,0.3769226,5.18];
ybs[6718]=['95 Her',4.7222035,0.3769274,4.96];
ybs[6719]=['',4.7721949,-1.324462,5.86];
ybs[6720]=['',4.7285557,-0.0935014,6.76];
ybs[6721]=['τ Oph',4.7300005,-0.1427464,5.94];
ybs[6722]=['τ Oph',4.7299932,-0.1427512,5.24];
ybs[6723]=['',4.6856369,1.3119351,6.36];
ybs[6724]=['9 Sgr',4.7339729,-0.4251392,5.97];
ybs[6725]=['',4.7221591,0.5814076,6.15];
ybs[6726]=['96 Her',4.726062,0.3636356,5.28];
ybs[6727]=['',4.7386516,-0.6265571,6];
ybs[6728]=['',4.7539284,-1.1265457,6.41];
ybs[6729]=['97 Her',4.726499,0.4001041,6.21];
ybs[6730]=['γ1 Sgr',4.7391785,-0.516227,4.69];
ybs[6731]=['θ Ara',4.7472706,-0.8742103,3.66];
ybs[6732]=['',4.7298556,0.3423393,6.5];
ybs[6733]=['π Pav',4.7571724,-1.1111521,4.35];
ybs[6734]=['γ2 Sgr',4.7426443,-0.5309546,2.99];
ybs[6735]=['',4.7364113,0.0335334,6.14];
ybs[6736]=['',4.7454416,-0.6286109,5.95];
ybs[6737]=['',4.7477125,-0.757849,5.77];
ybs[6738]=['',4.7477125,-0.757849,5.77];
ybs[6739]=['',4.7769643,-1.2857116,5.85];
ybs[6740]=['70 Oph',4.7400302,0.0436672,4.03];
ybs[6741]=['',4.7281275,0.8458892,6.21];
ybs[6742]=['',4.7359243,0.4179128,6.34];
ybs[6743]=['',4.7432693,-0.1452303,5.85];
ybs[6744]=['',4.743729,-0.0828776,5.77];
ybs[6745]=['',4.7430338,-0.007747,6.34];
ybs[6746]=['',4.7409065,0.2095529,7.04];
ybs[6747]=['',4.7551431,-0.7987208,6.15];
ybs[6748]=['',4.7625923,-1.0303628,6.38];
ybs[6749]=['ι Pav',4.765028,-1.0820595,5.49];
ybs[6750]=['',4.7483423,-0.3742093,6.28];
ybs[6751]=['',4.7396406,0.3778541,6.15];
ybs[6752]=['',4.7354666,0.6996378,6.52];
ybs[6753]=['98 Her',4.7419262,0.3878403,5.06];
ybs[6754]=['',4.7524936,-0.4966081,4.57];
ybs[6755]=['',4.7366494,0.7321466,6.34];
ybs[6756]=['',4.7406543,0.5625749,5.71];
ybs[6757]=['',4.75089,-0.299335,5.52];
ybs[6758]=['71 Oph',4.7479227,0.1524924,4.64];
ybs[6759]=['72 Oph',4.7480865,0.166979,3.73];
ybs[6760]=['',4.7584663,-0.639982,6.58];
ybs[6761]=['',4.7559526,-0.444509,6.61];
ybs[6762]=['',4.7839994,-1.2347305,6.73];
ybs[6763]=['99 Her',4.7459512,0.5334612,5.04];
ybs[6764]=['',4.7499709,0.2281948,6.63];
ybs[6765]=['',4.7610196,-0.5709885,6.43];
ybs[6766]=['',4.7664803,-0.8291723,6.07];
ybs[6767]=['ο Her',4.7482836,0.5020588,3.83];
ybs[6768]=['',4.7613598,-0.5362363,5.53];
ybs[6769]=['100 Her',4.7496129,0.4556159,5.86];
ybs[6770]=['100 Her',4.749613,0.455548,5.9];
ybs[6771]=['ε Tel',4.7670582,-0.8019683,4.53];
ybs[6772]=['',4.7532271,0.2493819,6.37];
ybs[6773]=['',4.7591499,-0.2431261,6.39];
ybs[6774]=['',4.7662012,-0.7217667,5.86];
ybs[6775]=['102 Her',4.7538719,0.3633483,4.36];
ybs[6776]=['',4.7650868,-0.5898311,6.16];
ybs[6777]=['δ UMi',4.5668661,1.5083105,4.36];
ybs[6778]=['',4.7442873,0.8870775,6.29];
ybs[6779]=['',4.7473294,0.7586066,5];
ybs[6780]=['',4.745299,0.8676672,6.32];
ybs[6781]=['',4.7501161,0.635386,5.48];
ybs[6782]=['101 Her',4.7544361,0.3499248,5.1];
ybs[6783]=['73 Oph',4.7579118,0.069771,5.73];
ybs[6784]=['',4.7819981,-1.1114785,6.47];
ybs[6785]=['',4.7594068,0.0545261,5.69];
ybs[6786]=['',4.765988,-0.3462248,6.36];
ybs[6787]=['',4.7553146,0.5318626,6.38];
ybs[6788]=['',4.7627677,0.0581001,5.51];
ybs[6789]=['11 Sgr',4.7682005,-0.4135715,4.98];
ybs[6790]=['',4.7694688,-0.5043315,6.51];
ybs[6791]=['',4.7600659,0.2876502,6.09];
ybs[6792]=['',4.7754349,-0.7213487,5.47];
ybs[6793]=['',4.7880628,-1.1004047,5.6];
ybs[6794]=['',4.7569501,0.6712836,6.4];
ybs[6795]=['',4.7586075,0.6365349,5.58];
ybs[6796]=['',4.79915,-1.1906832,6.33];
ybs[6797]=['40 Dra',4.7069053,1.3962734,6.04];
ybs[6798]=['41 Dra',4.7073238,1.3963324,5.68];
ybs[6799]=['24 UMi',4.5571801,1.515515,5.79];
ybs[6800]=['μ Sgr',4.7770098,-0.367441,3.86];
ybs[6801]=['',4.7738719,-0.0699157,6.59];
ybs[6802]=['',4.766457,0.5838488,5.88];
ybs[6803]=['104 Her',4.7672038,0.5482162,4.97];
ybs[6804]=['14 Sgr',4.7792209,-0.3788546,5.44];
ybs[6805]=['',4.7598708,0.9475605,5.95];
ybs[6806]=['',4.7872343,-0.7714295,5.46];
ybs[6807]=['',4.7935396,-0.9776592,5.33];
ybs[6808]=['',4.7735438,0.3819841,6.12];
ybs[6809]=['',4.7926428,-0.8911794,6.06];
ybs[6810]=['15 Sgr',4.7833323,-0.361661,5.38];
ybs[6811]=['16 Sgr',4.7833208,-0.355722,5.95];
ybs[6812]=['',4.7702767,0.718246,6.36];
ybs[6813]=['',4.784565,-0.3255839,6.07];
ybs[6814]=['',4.772011,0.6768264,6.04];
ybs[6815]=['',4.7618242,1.0544269,6.49];
ybs[6816]=['',4.8055835,-1.1148855,6.18];
ybs[6817]=['φ Oct',4.8256105,-1.3095843,5.47];
ybs[6818]=['',4.7860773,-0.0630155,6.36];
ybs[6819]=['',4.7796401,0.509874,6.56];
ybs[6820]=['η Sgr',4.7944794,-0.6414771,3.11];
ybs[6821]=['',4.794252,-0.5951485,6.16];
ybs[6822]=['',4.7864524,0.0416227,6.01];
ybs[6823]=['',4.7931504,-0.4999474,6.19];
ybs[6824]=['',4.7931147,-0.4936061,6.4];
ybs[6825]=['',4.8540388,-1.4000976,5.95];
ybs[6826]=['',4.7918528,-0.3031004,5.75];
ybs[6827]=['',4.7992922,-0.7379277,6.3];
ybs[6828]=['',4.7900662,-0.0523573,6];
ybs[6829]=['',4.7931175,-0.3221127,6.54];
ybs[6830]=['',4.79594,-0.4718428,4.65];
ybs[6831]=['',4.7925303,-0.1701873,6.31];
ybs[6832]=['',4.7907981,0.017685,6.63];
ybs[6833]=['',4.7830378,0.7359391,5.59];
ybs[6834]=['',4.7986867,-0.4467441,6.51];
ybs[6835]=['',4.7824141,0.7891708,6.29];
ybs[6836]=['',4.798562,-0.324828,6.84];
ybs[6837]=['',4.7778053,0.9877627,6.37];
ybs[6838]=['36 Dra',4.7733955,1.124047,5.03];
ybs[6839]=['',4.7946375,0.2405897,6.3];
ybs[6840]=['',4.7948472,0.3165896,5.99];
ybs[6841]=['',4.7895027,0.7146085,6.11];
ybs[6842]=['',4.7946723,0.4067406,6.63];
ybs[6843]=['ξ Pav',4.8207611,-1.0730928,4.36];
ybs[6844]=['',4.8088888,-0.6541209,6.45];
ybs[6845]=['',4.7996806,0.1268512,5.39];
ybs[6846]=['',4.804689,-0.2761617,5.39];
ybs[6847]=['δ Sgr',4.8088809,-0.5204381,2.7];
ybs[6848]=['105 Her',4.7992134,0.4268099,5.27];
ybs[6849]=['',4.8110017,-0.4346854,6.25];
ybs[6850]=['',4.8150014,-0.6745213,5.1];
ybs[6851]=['',4.8101902,-0.3290069,5.75];
ybs[6852]=['',4.813226,-0.49603,6.16];
ybs[6853]=['37 Dra',4.7786235,1.2001286,5.95];
ybs[6854]=['74 Oph',4.8072578,0.0591014,4.86];
ybs[6855]=['',4.8020294,0.517921,5.99];
ybs[6856]=['106 Her',4.8041877,0.3834516,4.95];
ybs[6857]=['η Ser',4.8093705,-0.0504339,3.26];
ybs[6858]=['',4.8174163,-0.6398286,5.34];
ybs[6859]=['',4.8309938,-1.0997356,6.14];
ybs[6860]=['κ Lyr',4.8017315,0.6295928,4.33];
ybs[6861]=['',4.8098436,0.0950355,6.13];
ybs[6862]=['',4.8199998,-0.6322999,5.55];
ybs[6863]=['',4.8239842,-0.7696849,5.25];
ybs[6864]=['108 Her',4.8067493,0.5212936,5.63];
ybs[6865]=['107 Her',4.8070714,0.5040348,5.12];
ybs[6866]=['',4.8171237,-0.1781742,6.33];
ybs[6867]=['ε Sgr',4.822928,-0.5999433,1.85];
ybs[6868]=['',4.8011671,0.896337,6.3];
ybs[6869]=['',4.8179059,-0.209521,5.73];
ybs[6870]=['',4.8121999,0.4065715,5.41];
ybs[6871]=['',4.814502,0.210114,5.89];
ybs[6872]=['ζ Sct',4.8197971,-0.1557518,4.68];
ybs[6873]=['',4.8153148,0.3113059,5.25];
ybs[6874]=['',4.8064501,0.8680328,6.4];
ybs[6875]=['',4.8163625,0.2914351,6.22];
ybs[6876]=['18 Sgr',4.8265,-0.5366156,5.6];
ybs[6877]=['',4.8281905,-0.627981,6.15];
ybs[6878]=['',4.8213825,-0.0623594,6.38];
ybs[6879]=['',4.8083543,0.8574961,5.05];
ybs[6880]=['',4.8242934,-0.1233005,6.31];
ybs[6881]=['',4.8304929,-0.5922608,6.3];
ybs[6882]=['',4.8355508,-0.8395951,5.46];
ybs[6883]=['109 Her',4.8190305,0.3801316,3.84];
ybs[6884]=['21 Sgr',4.8275456,-0.3583282,4.81];
ybs[6885]=['α Tel',4.8357416,-0.8020944,3.51];
ybs[6886]=['',4.8252138,-0.0273784,6.15];
ybs[6887]=['',4.8657992,-1.2906909,5.89];
ybs[6888]=['',4.825875,0.0889345,6.74];
ybs[6889]=['',4.8194764,0.6763056,6.36];
ybs[6890]=['',4.8279686,0.1403769,5.65];
ybs[6891]=['μ Lyr',4.8206315,0.6897127,5.12];
ybs[6892]=['',4.8244004,0.4783253,6.27];
ybs[6893]=['ζ Tel',4.8440488,-0.8562293,4.13];
ybs[6894]=['',4.8289627,0.2614125,6.37];
ybs[6895]=['',4.8386787,-0.5201841,5.92];
ybs[6896]=['',4.8496091,-1.0037396,5.76];
ybs[6897]=['',4.8381435,-0.4646543,6.31];
ybs[6898]=['',4.8418117,-0.6803858,5.64];
ybs[6899]=['',4.8178657,0.9304526,6.32];
ybs[6900]=['',4.9121686,-1.4274906,6.27];
ybs[6901]=['λ Sgr',4.8391464,-0.4434808,2.81];
ybs[6902]=['',4.8397769,-0.4667896,6.27];
ybs[6903]=['',4.8453961,-0.7650334,6.36];
ybs[6904]=['ν Pav',4.8563976,-1.086724,4.64];
ybs[6905]=['',4.8286928,0.5208071,5.83];
ybs[6906]=['59 Ser',4.8350153,0.0036275,5.21];
ybs[6907]=['',4.8387626,-0.310458,6.2];
ybs[6908]=['φ Dra',4.8018623,1.245232,4.22];
ybs[6909]=['',4.8454184,-0.6778587,6.63];
ybs[6910]=['',4.8487176,-0.823923,5.7];
ybs[6911]=['39 Dra',4.8178294,1.0264413,4.98];
ybs[6912]=['',4.8318829,0.46183,6.53];
ybs[6913]=['',4.8376594,0.0656348,6.07];
ybs[6914]=['',4.8434934,-0.4637193,6.5];
ybs[6915]=['χ Dra',4.8027302,1.2695811,3.57];
ybs[6916]=['',4.838205,0.1083188,5.73];
ybs[6917]=['',4.8452339,-0.4405858,6.59];
ybs[6918]=['γ Sct',4.8441422,-0.254002,4.7];
ybs[6919]=['',4.8532952,-0.7312963,6.04];
ybs[6920]=['',4.8466805,-0.2542741,5.96];
ybs[6921]=['',4.8486388,-0.3266536,5.66];
ybs[6922]=['δ1 Tel',4.8566048,-0.8011283,4.96];
ybs[6923]=['60 Ser',4.8458788,-0.0344266,5.39];
ybs[6924]=['',4.8530076,-0.5755355,5.34];
ybs[6925]=['',4.85725,-0.7591132,5.72];
ybs[6926]=['δ2 Tel',4.8578022,-0.7983726,5.07];
ybs[6927]=['',4.8652032,-1.0244151,6.44];
ybs[6928]=['',4.8484108,-0.0996785,6.28];
ybs[6929]=['',4.8474531,0.0711783,6.69];
ybs[6930]=['',4.8588821,-0.6927199,5.16];
ybs[6931]=['',4.8446889,0.4167638,5.9];
ybs[6932]=['',4.8540382,-0.3209528,5.14];
ybs[6933]=['42 Dra',4.8259642,1.1444934,4.82];
ybs[6934]=['',4.8537554,-0.1881868,5.72];
ybs[6935]=['',4.8560252,-0.3335546,6.68];
ybs[6936]=['',4.8617634,-0.696002,6.22];
ybs[6937]=['',4.8343165,1.0395349,6.43];
ybs[6938]=['',4.8495867,0.3635249,6.5];
ybs[6939]=['θ CrA',4.8640188,-0.7382401,4.64];
ybs[6940]=['κ1 CrA',4.8633346,-0.675545,6.32];
ybs[6941]=['κ2 CrA',4.8633203,-0.6756468,5.65];
ybs[6942]=['',4.869145,-0.9228783,6.22];
ybs[6943]=['',4.8513726,0.2956927,5.77];
ybs[6944]=['',4.857881,-0.255346,6.37];
ybs[6945]=['61 Ser',4.85574,-0.017267,5.94];
ybs[6946]=['',4.85633,0.0641149,6.43];
ybs[6947]=['',4.8595241,-0.2592073,5.5];
ybs[6948]=['',4.8655877,-0.5759942,5.28];
ybs[6949]=['24 Sgr',4.8649323,-0.4191918,5.49];
ybs[6950]=['',4.8635741,-0.258992,5.76];
ybs[6951]=['',4.8621314,-0.102923,6.36];
ybs[6952]=['',4.9571282,-1.4537481,7.16];
ybs[6953]=['25 Sgr',4.8677966,-0.4225031,6.51];
ybs[6954]=['',4.8585526,0.4124384,5.84];
ybs[6955]=['',4.861742,0.1445595,6.42];
ybs[6956]=['',4.8585641,0.533516,5.48];
ybs[6957]=['',4.871207,-0.3634761,6.48];
ybs[6958]=['',4.869519,-0.1913262,5.14];
ybs[6959]=['',4.8609653,0.5394251,6.59];
ybs[6960]=['',4.8743176,-0.5180781,6.37];
ybs[6961]=['α Sct',4.8701637,-0.1436242,3.85];
ybs[6962]=['',4.8545936,0.9098275,6.56];
ybs[6963]=['',4.8654532,0.3574624,6.57];
ybs[6964]=['',4.8677862,0.1903556,6.4];
ybs[6965]=['',4.8693824,0.3179711,5.78];
ybs[6966]=['45 Dra',4.8558407,0.9958748,4.77];
ybs[6967]=['',4.848955,1.1423069,6.59];
ybs[6968]=['',4.8704952,0.4122598,5.61];
ybs[6969]=['',4.8723748,0.2965473,6.21];
ybs[6970]=['ζ Pav',4.9090432,-1.2463299,4.01];
ybs[6971]=['',4.862235,0.9139953,5.36];
ybs[6972]=['',4.8688474,0.6016643,6.1];
ybs[6973]=['',4.8751341,0.1594902,5.39];
ybs[6974]=['',4.8893656,-0.8358879,5.86];
ybs[6975]=['',4.876029,0.1167213,5.45];
ybs[6976]=['',4.8823582,-0.3731779,5.94];
ybs[6977]=['',4.8828612,-0.2441437,6.47];
ybs[6978]=['',4.8850699,-0.4099513,5.81];
ybs[6979]=['',4.8905885,-0.753443,5.37];
ybs[6980]=['',4.8783225,0.1996236,6.42];
ybs[6981]=['',4.8803722,-0.0051198,5.75];
ybs[6982]=['',4.9326929,-1.3586775,6.39];
ybs[6983]=['',4.8779123,0.2829916,6.29];
ybs[6984]=['',4.9047213,-1.1279154,6.37];
ybs[6985]=['',4.8749775,0.5844152,5.42];
ybs[6986]=['',4.8866367,-0.3671347,5.86];
ybs[6987]=['',4.8839256,-0.0554521,6.49];
ybs[6988]=['',4.8835301,-0.019145,6.66];
ybs[6989]=['α Lyr',4.8761276,0.6771767,0.03];
ybs[6990]=['',4.8833752,0.1544667,6.4];
ybs[6991]=['',4.8751205,0.7546387,6.2];
ybs[6992]=['',4.9101353,-1.1263015,5.78];
ybs[6993]=['',4.8992863,-0.839105,6.49];
ybs[6994]=['',4.838463,1.3536665,5.64];
ybs[6995]=['',4.8911006,-0.135672,5.84];
ybs[6996]=['',4.889001,0.0921726,6.38];
ybs[6997]=['',4.8811907,0.6926218,6.04];
ybs[6998]=['',4.890008,0.1287245,6.28];
ybs[6999]=['26 Sgr',4.8996892,-0.4156574,6.23];
ybs[7000]=['',4.9181582,-1.1318781,4.79];
ybs[7001]=['',4.8706652,1.1432597,6.06];
ybs[7002]=['',4.8987235,-0.2538813,6.42];
ybs[7003]=['',4.916505,-1.06597,6.04];
ybs[7004]=['',4.8899805,0.5387224,6.36];
ybs[7005]=['',4.8874262,0.7147446,6.25];
ybs[7006]=['',4.8769543,1.0915742,5.74];
ybs[7007]=['',4.890397,0.6699333,6.45];
ybs[7008]=['δ Sct',4.9010207,-0.1576805,4.72];
ybs[7009]=['λ CrA',4.9086545,-0.668546,5.13];
ybs[7010]=['',4.9169286,-0.9924324,6.22];
ybs[7011]=['',4.9041603,-0.3362471,6.35];
ybs[7012]=['',4.9023943,-0.1231353,6.15];
ybs[7013]=['',4.8075442,1.4518527,6.17];
ybs[7014]=['',4.9100925,-0.640531,6.32];
ybs[7015]=['',4.9386892,-1.2736244,6.06];
ybs[7016]=['',4.8881527,0.9112899,6];
ybs[7017]=['',4.910895,-0.6217382,4.87];
ybs[7018]=['',4.8971725,0.5521443,6.41];
ybs[7019]=['',4.9138312,-0.6923224,5.43];
ybs[7020]=['ε Sct',4.9064441,-0.1441064,4.9];
ybs[7021]=['',4.8989914,0.6067569,6.47];
ybs[7022]=['',4.9078626,-0.1186804,6.31];
ybs[7023]=['',4.9126732,-0.4361918,5.83];
ybs[7024]=['θ Pav',4.9320638,-1.1354572,5.73];
ybs[7025]=['',4.923302,-0.8739619,6.54];
ybs[7026]=['',4.9146514,-0.3662056,6.36];
ybs[7027]=['φ Sgr',4.9163655,-0.4707384,3.17];
ybs[7028]=['4 Aql',4.9118573,0.0362875,5.02];
ybs[7029]=['',4.9037726,0.6862405,6.45];
ybs[7030]=['',4.8916797,1.095492,6.09];
ybs[7031]=['',4.9053124,0.638358,6.01];
ybs[7032]=['',4.9066375,0.5575513,5.7];
ybs[7033]=['',4.9176956,-0.3418531,6.42];
ybs[7034]=['28 Sgr',4.9191991,-0.3904727,5.37];
ybs[7035]=['',4.9104835,0.4120503,6.31];
ybs[7036]=['',4.9145566,0.0963266,5.83];
ybs[7037]=['46 Dra',4.8998827,0.9696619,5.04];
ybs[7038]=['μ CrA',4.926037,-0.7048636,5.24];
ybs[7039]=['ε1 Lyr',4.9083882,0.6927011,5.06];
ybs[7040]=['ε1 Lyr',4.9083809,0.6927205,6.02];
ybs[7041]=['ε2 Lyr',4.9085728,0.6917076,5.14];
ybs[7042]=['ε2 Lyr',4.9085728,0.6917028,5.37];
ybs[7043]=['',4.9204612,-0.1763668,5.71];
ybs[7044]=['ζ1 Lyr',4.9103781,0.6566634,4.36];
ybs[7045]=['ζ2 Lyr',4.9105095,0.6564793,5.73];
ybs[7046]=['',4.9145751,0.3840494,6.51];
ybs[7047]=['5 Aql',4.9191156,-0.0164385,5.9];
ybs[7048]=['',4.9037462,0.9405649,6.11];
ybs[7049]=['110 Her',4.9149143,0.3589414,4.19];
ybs[7050]=['η1 CrA',4.9310001,-0.7619956,5.49];
ybs[7051]=['β Sct',4.9222801,-0.0825135,4.22];
ybs[7052]=['',4.9165059,0.4656857,4.83];
ybs[7053]=['',4.9338049,-0.7991714,5.81];
ybs[7054]=['',4.9236535,-0.0992178,5.2];
ybs[7055]=['',4.919462,0.3268252,6.17];
ybs[7056]=['η2 CrA',4.9342215,-0.7576948,5.61];
ybs[7057]=['111 Her',4.9209263,0.3176743,4.36];
ybs[7058]=['',4.9325074,-0.6061106,6.62];
ybs[7059]=['',4.9099442,0.9584601,6.23];
ybs[7060]=['',4.9295951,-0.3242928,6.47];
ybs[7061]=['',4.9164897,0.7236361,6.07];
ybs[7062]=['λ Pav',4.9472768,-1.0849867,4.22];
ybs[7063]=['',4.9065948,1.065817,5.99];
ybs[7064]=['',4.9258147,0.0743832,6.21];
ybs[7065]=['',4.9332628,-0.333726,6.75];
ybs[7066]=['29 Sgr',4.933636,-0.3543639,5.24];
ybs[7067]=['',4.926211,0.4107578,6.15];
ybs[7068]=['',4.9195626,0.8086972,6.52];
ybs[7069]=['',4.9245305,0.5546188,6.06];
ybs[7070]=['',4.8998696,1.2358848,6.44];
ybs[7071]=['',4.9332585,-0.1028283,5.99];
ybs[7072]=['',4.9179404,0.9251617,5.88];
ybs[7073]=['',4.9327788,0.0149565,6.25];
ybs[7074]=['',4.9290423,0.3377107,5.88];
ybs[7075]=['κ Tel',4.9482002,-0.9090551,5.17];
ybs[7076]=['30 Sgr',4.938816,-0.3864258,6.61];
ybs[7077]=['',4.9361532,-0.137638,6.8];
ybs[7078]=['',4.9223682,0.8568731,6.4];
ybs[7079]=['',4.9303761,0.4375071,6.59];
ybs[7080]=['',4.9469153,-0.8128505,5.54];
ybs[7081]=['',4.950549,-0.9059727,6.31];
ybs[7082]=['',4.9390082,-0.1702128,5.83];
ybs[7083]=['',4.9496,-0.8436464,6.19];
ybs[7084]=['',4.9250152,0.8515106,6.12];
ybs[7085]=['',4.9492771,-0.8126818,6.19];
ybs[7086]=['',4.9322823,0.5524016,6.64];
ybs[7087]=['',4.9374629,0.1919506,6.55];
ybs[7088]=['ν1 Lyr',4.9323726,0.5730597,5.91];
ybs[7089]=['8 Aql',4.9405349,-0.0575248,6.1];
ybs[7090]=['ν2 Lyr',4.9328924,0.5684888,5.25];
ybs[7091]=['',4.9460982,-0.4647452,6.29];
ybs[7092]=['',4.9468173,-0.512377,6.13];
ybs[7093]=['',4.9472112,-0.5360208,6.63];
ybs[7094]=['β Lyr',4.9337246,0.5826613,3.45];
ybs[7095]=['κ Pav',4.9687148,-1.1730232,4.44];
ybs[7096]=['',4.9562296,-0.8701403,6.6];
ybs[7097]=['',4.9429213,0.2441305,6.14];
ybs[7098]=['',4.9479756,-0.1667413,6.34];
ybs[7099]=['',4.9678293,-1.0956624,6.48];
ybs[7100]=['',4.940532,0.5027509,6.18];
ybs[7101]=['112 Her',4.9437378,0.3743288,5.48];
ybs[7102]=['33 Sgr',4.9525751,-0.3723968,5.69];
ybs[7103]=['',4.940258,0.638101,6.09];
ybs[7104]=['ν1 Sgr',4.953356,-0.3965731,4.83];
ybs[7105]=['',4.9101981,1.2933727,5.27];
ybs[7106]=['',4.9422518,0.7226608,6.28];
ybs[7107]=['',4.955517,-0.2719189,5.1];
ybs[7108]=['ν2 Sgr',4.9574907,-0.3952816,4.99];
ybs[7109]=['σ Sgr',4.9582566,-0.4585534,2.02];
ybs[7110]=['',4.9633985,-0.7450219,5.36];
ybs[7111]=['',4.9391791,0.9249691,5.51];
ybs[7112]=['50 Dra',4.9122132,1.3169093,5.35];
ybs[7113]=['ο Dra',4.9369085,1.0368997,4.66];
ybs[7114]=['',4.9590241,-0.2854153,5.79];
ybs[7115]=['ω Pav',4.9748909,-1.0502622,5.14];
ybs[7116]=['',4.9613977,-0.4040405,5.93];
ybs[7117]=['',4.9648504,-0.6513436,5.38];
ybs[7118]=['',4.9821552,-1.1628732,6.01];
ybs[7119]=['δ1 Lyr',4.9494841,0.6456745,5.58];
ybs[7120]=['',4.9520062,0.4875128,5.62];
ybs[7121]=['113 Her',4.9544977,0.3956348,4.59];
ybs[7122]=['λ Tel',4.9735656,-0.923519,4.87];
ybs[7123]=['',4.9580934,0.1158689,5.57];
ybs[7124]=['',4.9689073,-0.6946212,6.31];
ybs[7125]=['',4.9465249,0.8854203,4.92];
ybs[7126]=['',4.9515133,0.7199224,7.3];
ybs[7127]=['δ2 Lyr',4.9528772,0.6444099,4.3];
ybs[7128]=['',4.9546268,0.5932697,6.02];
ybs[7129]=['θ1 Ser',4.961487,0.073783,4.62];
ybs[7130]=['θ2 Ser',4.9615888,0.0737541,4.98];
ybs[7131]=['',4.9623509,-0.0309985,6.22];
ybs[7132]=['',4.9624388,0.0435466,6.15];
ybs[7133]=['ξ1 Sgr',4.9671231,-0.3600972,5.08];
ybs[7134]=['',4.954251,0.7265108,5.44];
ybs[7135]=['',4.9603962,0.3144866,6.63];
ybs[7136]=['',4.96056,0.3164116,5.69];
ybs[7137]=['η Sct',4.9654493,-0.1016115,4.83];
ybs[7138]=['ξ2 Sgr',4.968832,-0.3679532,3.51];
ybs[7139]=['',4.9719133,-0.54125,6.12];
ybs[7140]=['ε CrA',4.9737744,-0.6472129,4.87];
ybs[7141]=['',4.9483302,1.0037329,6.22];
ybs[7142]=['',4.9534557,0.8531673,5.77];
ybs[7143]=['',4.9716282,-0.4337476,6.62];
ybs[7144]=['',4.9759012,-0.6895725,6.49];
ybs[7145]=['13 Lyr',4.9561553,0.7674128,4.04];
ybs[7146]=['64 Ser',4.9661458,0.0446726,5.57];
ybs[7147]=['',4.971846,-0.3927805,6.14];
ybs[7148]=['',4.9060104,1.3955916,6.39];
ybs[7149]=['',4.997551,-1.1995343,5.88];
ybs[7150]=['',4.964062,0.5746587,5.22];
ybs[7151]=['',4.9709258,0.1093451,6.21];
ybs[7152]=['',4.9762388,-0.3236145,6.37];
ybs[7153]=['',4.9699295,0.303434,5.38];
ybs[7154]=['',4.9758434,-0.2236707,5.53];
ybs[7155]=['10 Aql',4.9723785,0.2431514,5.89];
ybs[7156]=['',4.980668,-0.434877,6.36];
ybs[7157]=['',4.9839354,-0.6463817,6.69];
ybs[7158]=['',4.9840154,-0.646401,6.4];
ybs[7159]=['',4.9720629,0.3459073,6.5];
ybs[7160]=['11 Aql',4.9737545,0.238194,5.23];
ybs[7161]=['',4.974717,0.1774289,6.75];
ybs[7162]=['',4.9682148,0.6682977,5.89];
ybs[7163]=['48 Dra',4.9613009,1.0094799,5.66];
ybs[7164]=['ε Aql',4.9760157,0.2634323,4.02];
ybs[7165]=['',4.9889299,-0.7310076,6.23];
ybs[7166]=['γ Lyr',4.9724359,0.5709733,3.24];
ybs[7167]=['',4.9713502,0.7104188,6.22];
ybs[7168]=['υ Dra',4.9487741,1.2447703,4.82];
ybs[7169]=['',4.9762378,0.4582505,5.27];
ybs[7170]=['',4.9858988,-0.3956567,6.24];
ybs[7171]=['',4.9772726,0.3986345,6.29];
ybs[7172]=['',4.964417,1.0166458,6.46];
ybs[7173]=['',4.9732912,0.6849158,6.41];
ybs[7174]=['',4.9853487,-0.2662751,6.32];
ybs[7175]=['',4.9589299,1.1393826,5.63];
ybs[7176]=['ζ CrA',4.9931701,-0.7342344,4.75];
ybs[7177]=['',4.9973834,-0.8899695,5.93];
ybs[7178]=['',4.9631249,1.0894486,6.45];
ybs[7179]=['λ Lyr',4.9771351,0.5614884,4.93];
ybs[7180]=['12 Aql',4.9855966,-0.099707,4.02];
ybs[7181]=['ζ Sgr',4.990437,-0.5210464,2.6];
ybs[7182]=['',4.9896003,-0.4331995,5.65];
ybs[7183]=['',4.9716922,0.8872265,6.3];
ybs[7184]=['',4.9937673,-0.6671787,5.74];
ybs[7185]=['',4.9822903,0.3374689,6.39];
ybs[7186]=['',4.9434341,1.323133,6.22];
ybs[7187]=['',4.9834846,0.3640678,6.69];
ybs[7188]=['',4.978072,0.7105171,6.65];
ybs[7189]=['',4.9829203,0.4593232,5.69];
ybs[7190]=['',4.9920387,-0.3354278,6.05];
ybs[7191]=['',4.9810175,0.5904091,6.01];
ybs[7192]=['',4.9922669,-0.33295,6.37];
ybs[7193]=['',4.9842377,0.4372372,6.72];
ybs[7194]=['',4.9853938,0.389034,6.4];
ybs[7195]=['',4.9881579,0.1466168,6.3];
ybs[7196]=['14 Aql',4.9908973,-0.0640936,5.42];
ybs[7197]=['',4.9771035,0.8824212,5.38];
ybs[7198]=['',4.9983627,-0.5413957,5.5];
ybs[7199]=['',4.9848878,0.5872592,6.39];
ybs[7200]=['ρ Tel',5.0078378,-0.9130297,5.16];
ybs[7201]=['',4.9934812,0.0322142,5.83];
ybs[7202]=['16 Lyr',4.9826305,0.8196176,5.01];
ybs[7203]=['',4.9900687,0.3436145,6.09];
ybs[7204]=['ο Sgr',4.9991872,-0.3789864,3.77];
ybs[7205]=['49 Dra',4.9788553,0.9718677,5.48];
ybs[7206]=['',4.9962379,0.0586023,6.73];
ybs[7207]=['',4.9974772,-0.0987471,6.9];
ybs[7208]=['',5.0252978,-1.1937197,5.33];
ybs[7209]=['',4.9936459,0.371662,6.52];
ybs[7210]=['',5.0101504,-0.8424853,5.97];
ybs[7211]=['',4.9687911,1.2139784,6.52];
ybs[7212]=['15 Aql',4.9998586,-0.0698822,5.42];
ybs[7213]=['γ CrA',5.0073363,-0.6463871,4.93];
ybs[7214]=['γ CrA',5.0073363,-0.6463871,4.99];
ybs[7215]=['σ Oct',5.5991316,-1.5510407,5.47];
ybs[7216]=['',4.9852203,0.9125852,6.31];
ybs[7217]=['',5.0033711,-0.2728391,5.97];
ybs[7218]=['',5.0013138,-0.0259217,6.53];
ybs[7219]=['',5.0093618,-0.6594204,6.16];
ybs[7220]=['',5.0191765,-0.9719939,6.49];
ybs[7221]=['τ Sgr',5.0092361,-0.4824488,3.32];
ybs[7222]=['ζ Aql',5.0013045,0.2424424,2.99];
ybs[7223]=['λ Aql',5.0055042,-0.0847277,3.44];
ybs[7224]=['',4.9987556,0.5545231,5.56];
ybs[7225]=['',4.9988242,0.5368759,6.06];
ybs[7226]=['',5.0085585,-0.2827548,6.03];
ybs[7227]=['',5.0117796,-0.4993114,6.04];
ybs[7228]=['',5.0098147,-0.3265412,6.29];
ybs[7229]=['δ CrA',5.0159145,-0.7062963,4.59];
ybs[7230]=['',5.0056602,0.1441292,6.09];
ybs[7231]=['',5.002405,0.5227154,6.31];
ybs[7232]=['',5.0092883,0.0116935,6.56];
ybs[7233]=['',5.0148154,-0.4298422,6.3];
ybs[7234]=['',4.9620015,1.3452134,6.54];
ybs[7235]=['18 Aql',5.0082244,0.193725,5.09];
ybs[7236]=['',5.014791,-0.3361712,5.54];
ybs[7237]=['',5.0063403,0.423747,5.77];
ybs[7238]=['51 Dra',4.9973707,0.9324243,5.38];
ybs[7239]=['',4.9986982,0.8718051,6.43];
ybs[7240]=['',5.006134,0.5001534,5.55];
ybs[7241]=['α CrA',5.0206917,-0.6610456,4.11];
ybs[7242]=['',5.0216101,-0.6946175,6.46];
ybs[7243]=['',5.0212034,-0.6306808,6.56];
ybs[7244]=['',5.0230279,-0.7306467,5.88];
ybs[7245]=['',5.0040804,0.7232957,6.49];
ybs[7246]=['β CrA',5.0231862,-0.6861113,4.11];
ybs[7247]=['',5.0123193,0.2946458,6.07];
ybs[7248]=['17 Lyr',5.0094699,0.5677565,5.23];
ybs[7249]=['ι Lyr',5.0087744,0.6305631,5.28];
ybs[7250]=['',5.0126231,0.3792173,6.23];
ybs[7251]=['π Sgr',5.021325,-0.3664179,2.89];
ybs[7252]=['',5.021453,-0.3451246,6.13];
ybs[7253]=['19 Aql',5.0171914,0.106507,5.22];
ybs[7254]=['',5.0154396,0.2946169,6.48];
ybs[7255]=['',5.0275391,-0.6802429,6.36];
ybs[7256]=['',5.0211359,-0.0069575,6.34];
ybs[7257]=['',5.0283876,-0.5143864,6.3];
ybs[7258]=['',5.03577,-0.8806181,6.13];
ybs[7259]=['',5.0165715,0.6044006,6.74];
ybs[7260]=['',5.032417,-0.6554124,6.57];
ybs[7261]=['τ Pav',5.0542055,-1.2070395,6.27];
ybs[7262]=['',5.0127719,0.9155003,5.81];
ybs[7263]=['',5.0331388,-0.377472,6.41];
ybs[7264]=['',5.0366044,-0.4516187,5.8];
ybs[7265]=['',5.0568627,-1.1628926,5.53];
ybs[7266]=['20 Aql',5.0336462,-0.1380358,5.34];
ybs[7267]=['',5.0275302,0.4671528,6.36];
ybs[7268]=['',5.043695,-0.7882238,5.92];
ybs[7269]=['',5.0363031,-0.2138323,5.51];
ybs[7270]=['19 Lyr',5.0284605,0.5465234,5.98];
ybs[7271]=['',5.0263992,0.7061453,6.18];
ybs[7272]=['',5.0324739,0.2945573,6.73];
ybs[7273]=['',5.0324922,0.3767284,5.93];
ybs[7274]=['21 Aql',5.0378653,0.0405717,5.15];
ybs[7275]=['',5.037876,0.0968102,6.49];
ybs[7276]=['',5.0511239,-0.7929776,5.4];
ybs[7277]=['55 Dra',5.0170601,1.1520539,6.25];
ybs[7278]=['',5.0466874,-0.4214519,6.25];
ybs[7279]=['ψ Sgr',5.0466658,-0.4402579,4.85];
ybs[7280]=['',5.0288787,0.870647,6.75];
ybs[7281]=['',5.0289077,0.870681,6.57];
ybs[7282]=['53 Dra',5.0265648,0.992904,5.12];
ybs[7283]=['',5.0513629,-0.5845021,6.25];
ybs[7284]=['',5.0594784,-0.9311995,6.38];
ybs[7285]=['η Lyr',5.0368206,0.6837685,4.39];
ybs[7286]=['',5.0431481,0.3531643,6];
ybs[7287]=['',5.0445794,0.2638104,5.57];
ybs[7288]=['1 Sge',5.0441916,0.3711235,5.64];
ybs[7289]=['',5.0444062,0.5333379,5.85];
ybs[7290]=['22 Aql',5.0500332,0.0849421,5.59];
ybs[7291]=['43 Sgr',5.0555989,-0.3302246,4.96];
ybs[7292]=['',5.0468613,0.4797455,6.54];
ybs[7293]=['1 Vul',5.0482307,0.3738886,4.77];
ybs[7294]=['',5.0494469,0.2544128,5.63];
ybs[7295]=['',5.047019,0.4879731,6.16];
ybs[7296]=['54 Dra',5.0362867,1.0076823,4.99];
ybs[7297]=['δ Dra',5.0289479,1.1814483,3.07];
ybs[7298]=['',5.043032,0.8744513,6.27];
ybs[7299]=['59 Dra',5.0112795,1.3367369,5.13];
ybs[7300]=['',5.0557263,0.0360286,6.19];
ybs[7301]=['θ Lyr',5.0482606,0.6661154,4.36];
ybs[7302]=['ω1 Aql',5.0555125,0.202945,5.28];
ybs[7303]=['',5.0650492,-0.6176366,5.59];
ybs[7304]=['',5.0614564,-0.270583,6.06];
ybs[7305]=['2 Vul',5.0547673,0.4024401,5.43];
ybs[7306]=['23 Aql',5.0589744,0.0195162,5.1];
ybs[7307]=['',5.0872077,-1.1926847,6.34];
ybs[7308]=['24 Aql',5.0603333,0.0064963,6.41];
ybs[7309]=['',5.0499308,0.8208516,6];
ybs[7310]=['',5.0682628,-0.5547365,6.58];
ybs[7311]=['',5.0557358,0.54201,6.68];
ybs[7312]=['',5.0602026,0.1684434,6.32];
ybs[7313]=['',5.0595981,0.3428447,6.58];
ybs[7314]=['',5.0687959,-0.3904077,5.58];
ybs[7315]=['κ Cyg',5.0505698,0.9320208,3.77];
ybs[7316]=['η Tel',5.0800307,-0.9492649,5.05];
ybs[7317]=['',5.0730198,-0.609988,6.48];
ybs[7318]=['28 Aql',5.063512,0.2165616,5.53];
ybs[7319]=['ω2 Aql',5.0645329,0.2019073,6.02];
ybs[7320]=['26 Aql',5.0679117,-0.0939355,5.01];
ybs[7321]=['',5.0762254,-0.7327186,6.34];
ybs[7322]=['',5.0602266,0.5833233,6.6];
ybs[7323]=['27 Aql',5.0679889,-0.0149833,5.49];
ybs[7324]=['β1 Sgr',5.078439,-0.7753497,4.01];
ybs[7325]=['',5.0598666,0.6541203,6.22];
ybs[7326]=['',5.0729847,-0.3351031,6.26];
ybs[7327]=['ρ1 Sgr',5.0731818,-0.310896,3.93];
ybs[7328]=['',5.0575033,0.8657285,6.31];
ybs[7329]=['υ Sgr',5.073357,-0.2778702,4.61];
ybs[7330]=['β2 Sgr',5.0809941,-0.7812944,4.29];
ybs[7331]=['ρ2 Sgr',5.0739594,-0.3189427,5.87];
ybs[7332]=['',5.0626145,0.6521224,6.31];
ybs[7333]=['',5.0666497,0.6147012,6.31];
ybs[7334]=['',5.0758291,-0.1425352,6.31];
ybs[7335]=['α Sgr',5.0836866,-0.7082721,3.97];
ybs[7336]=['',5.075665,-0.0038059,5.83];
ybs[7337]=['',5.0858927,-0.7624951,6.17];
ybs[7338]=['',5.0614274,0.9496224,6.26];
ybs[7339]=['τ Dra',5.0405631,1.2808453,4.45];
ybs[7340]=['',5.0789987,-0.1285529,6.32];
ybs[7341]=['',5.0773364,0.1736195,6.35];
ybs[7342]=['',5.0858478,-0.4857291,6.04];
ybs[7343]=['',5.064017,1.0066845,5.91];
ybs[7344]=['',5.0786307,0.2610287,6.64];
ybs[7345]=['3 Vul',5.0770116,0.458971,5.18];
ybs[7346]=['',5.0754698,0.585602,6.06];
ybs[7347]=['',5.0883712,-0.5109209,5.93];
ybs[7348]=['',5.0610205,1.1244124,6.52];
ybs[7349]=['χ1 Sgr',5.0891059,-0.4271341,5.03];
ybs[7350]=['χ3 Sgr',5.0900473,-0.4175963,5.43];
ybs[7351]=['',5.0812825,0.3542918,6.4];
ybs[7352]=['',5.0690494,1.008816,6.43];
ybs[7353]=['',5.0874361,-0.0846251,6.52];
ybs[7354]=['',5.0891467,-0.2419252,5.69];
ybs[7355]=['',5.0798902,0.5804459,6.37];
ybs[7356]=['2 Sge',5.0839454,0.2962347,6.25];
ybs[7357]=['',5.1016022,-0.9475146,5.69];
ybs[7358]=['π Dra',5.0647444,1.1475248,4.59];
ybs[7359]=['2 Cyg',5.0824718,0.5176034,4.97];
ybs[7360]=['31 Aql',5.0867185,0.2090888,5.16];
ybs[7361]=['',5.0836038,0.4908386,6.53];
ybs[7362]=['50 Sgr',5.0935751,-0.3794455,5.59];
ybs[7363]=['',5.0820926,0.6368187,6.36];
ybs[7364]=['δ Aql',5.0892727,0.054985,3.36];
ybs[7365]=['',5.09278,-0.2620975,5.72];
ybs[7366]=['',5.0937467,-0.2533353,6.7];
ybs[7367]=['',5.0965581,-0.5184856,5.67];
ybs[7368]=['',5.0782765,0.8780081,6.51];
ybs[7369]=['',5.0810661,0.7578753,5.84];
ybs[7370]=['',5.1181134,-1.1937317,5.96];
ybs[7371]=['',5.0882288,0.354424,6.31];
ybs[7372]=['4 Vul',5.0886946,0.3461732,5.16];
ybs[7373]=['',5.0883255,0.4354318,6.19];
ybs[7374]=['ν Aql',5.0938001,0.0065398,4.66];
ybs[7375]=['',5.1107402,-0.9669801,6.13];
ybs[7376]=['',5.0929349,0.2279385,5.74];
ybs[7377]=['5 Vul',5.0919293,0.3513998,5.63];
ybs[7378]=['',5.0930632,0.3477994,5.81];
ybs[7379]=['',5.1078633,-0.757622,5.71];
ybs[7380]=['μ Tel',5.1137565,-0.9611915,6.3];
ybs[7381]=['λ UMi',4.4310843,1.5533876,6.38];
ybs[7382]=['4 Cyg',5.0910472,0.6344912,5.15];
ybs[7383]=['',5.097975,0.2499135,6.32];
ybs[7384]=['',5.1017064,0.0517854,5.85];
ybs[7385]=['',5.109233,-0.4703333,5.52];
ybs[7386]=['',5.1110536,-0.5594586,6.6];
ybs[7387]=['35 Aql',5.104657,0.0346858,5.8];
ybs[7388]=['',5.0880445,1.0133889,6.6];
ybs[7389]=['',5.1063955,-0.1222895,6.61];
ybs[7390]=['',5.0973357,0.6628337,6.34];
ybs[7391]=['',5.1059399,0.0049444,6.25];
ybs[7392]=['α Vul',5.1026279,0.4311299,4.44];
ybs[7393]=['8 Vul',5.1036937,0.4329398,5.81];
ybs[7394]=['',5.1058356,0.2553945,5.56];
ybs[7395]=['ι1 Cyg',5.0957511,0.9138006,5.75];
ybs[7396]=['7 Vul',5.1055767,0.3545968,6.33];
ybs[7397]=['',5.1135393,-0.371308,6.13];
ybs[7398]=['',5.1237243,-0.9275927,5.75];
ybs[7399]=['',5.109685,0.0513422,6.09];
ybs[7400]=['',5.0904293,1.0924564,6.38];
ybs[7401]=['36 Aql',5.1119718,-0.0480169,5.03];
ybs[7402]=['',5.1113138,0.0607744,6.05];
ybs[7403]=['',5.1252356,-0.7894665,5.61];
ybs[7404]=['β1 Cyg',5.1113175,0.4886473,3.08];
ybs[7405]=['β2 Cyg',5.1114628,0.4887445,5.11];
ybs[7406]=['',5.1112627,0.6329668,6.25];
ybs[7407]=['ι2 Cyg',5.1057214,0.9035042,3.79];
ybs[7408]=['',5.1141647,0.4652207,5.87];
ybs[7409]=['',5.1283804,-0.6980546,5.89];
ybs[7410]=['',5.0639476,1.38992,6.05];
ybs[7411]=['ι Tel',5.1334941,-0.8387981,4.9];
ybs[7412]=['',5.0300123,1.4572396,6.53];
ybs[7413]=['8 Cyg',5.1156574,0.6019844,4.74];
ybs[7414]=['',5.1128798,0.8786784,5.53];
ybs[7415]=['',5.1120524,0.9733663,6.37];
ybs[7416]=['μ Aql',5.1266453,0.1294675,4.45];
ybs[7417]=['37 Aql',5.1316383,-0.1836227,5.12];
ybs[7418]=['51 Sgr',5.1360188,-0.4307354,5.65];
ybs[7419]=['',5.1331882,-0.129515,6.34];
ybs[7420]=['',5.1335984,-0.2131593,6.27];
ybs[7421]=['',5.1482586,-1.0112871,6.18];
ybs[7422]=['',5.1556274,-1.163159,6.39];
ybs[7423]=['',5.1234854,0.6772058,6.61];
ybs[7424]=['9 Vul',5.1284327,0.3457944,5];
ybs[7425]=['',5.1325712,0.051538,6.38];
ybs[7426]=['',5.1376084,-0.328345,6.11];
ybs[7427]=['52 Sgr',5.138983,-0.4336009,4.6];
ybs[7428]=['9 Cyg',5.1292768,0.5149135,5.38];
ybs[7429]=['',5.1233048,0.8604704,5.96];
ybs[7430]=['',5.1402943,-0.3174907,5.64];
ybs[7431]=['',5.1280281,0.7409273,5.35];
ybs[7432]=['',5.1354569,0.1952995,6.68];
ybs[7433]=['κ Aql',5.1392625,-0.1219521,4.95];
ybs[7434]=['ι Aql',5.1383647,-0.0217521,4.36];
ybs[7435]=['',5.120102,1.0506391,6.29];
ybs[7436]=['',5.1359318,0.2518781,6.38];
ybs[7437]=['',5.1088224,1.2396566,6.07];
ybs[7438]=['',5.1259535,0.8949303,5.73];
ybs[7439]=['',5.1351389,0.3948923,6.32];
ybs[7440]=['',5.1276164,0.8413173,6.67];
ybs[7441]=['',5.1424403,-0.2489056,5.47];
ybs[7442]=['',5.1631422,-1.1486368,6.09];
ybs[7443]=['',5.13869,0.1974572,5.98];
ybs[7444]=['11 Cyg',5.1331619,0.6454946,6.05];
ybs[7445]=['',5.137348,0.3555724,7.14];
ybs[7446]=['',5.1561136,-0.9490391,6.26];
ybs[7447]=['42 Aql',5.1431098,-0.0804073,5.46];
ybs[7448]=['',5.1528798,-0.7895303,6.25];
ybs[7449]=['σ Dra',5.1150746,1.2164823,4.68];
ybs[7450]=['ε Sge',5.1403521,0.2880327,5.66];
ybs[7451]=['',5.1535872,-0.6875194,6.61];
ybs[7452]=['',5.1330156,0.877522,6.52];
ybs[7453]=['',5.1393573,0.5126699,6.43];
ybs[7454]=['',5.1380692,0.6706253,6.5];
ybs[7455]=['',5.1364024,0.7807727,5.17];
ybs[7456]=['θ Cyg',5.1352439,0.8772199,4.48];
ybs[7457]=['53 Sgr',5.1525306,-0.4081709,6.34];
ybs[7458]=['',5.1473798,0.0597348,6.35];
ybs[7459]=['',5.1445937,0.3634376,6.48];
ybs[7460]=['',5.1538174,-0.4081834,5.97];
ybs[7461]=['σ Aql',5.1489687,0.0949251,5.17];
ybs[7462]=['',5.1496592,0.2899427,6.38];
ybs[7463]=['54 Sgr',5.1562396,-0.2836454,6.2];
ybs[7464]=['',5.1418696,0.8608822,6.47];
ybs[7465]=['φ Cyg',5.1490185,0.5269917,4.69];
ybs[7466]=['α Sge',5.1525553,0.3151236,4.37];
ybs[7467]=['45 Aql',5.1557998,-0.0101139,5.67];
ybs[7468]=['',5.1505047,0.5937675,6.1];
ybs[7469]=['',5.1541177,0.3581097,6.5];
ybs[7470]=['14 Cyg',5.1487688,0.7480377,5.4];
ybs[7471]=['',5.1446988,0.9601865,5.82];
ybs[7472]=['',5.1548456,0.4146741,6.64];
ybs[7473]=['',5.1570197,0.2418556,6.01];
ybs[7474]=['',5.1491779,0.802837,6.2];
ybs[7475]=['β Sge',5.1567246,0.305744,4.37];
ybs[7476]=['55 Sgr',5.1640652,-0.280676,5.06];
ybs[7477]=['',5.1574463,0.3926044,6.36];
ybs[7478]=['',5.1696321,-0.6544304,6.16];
ybs[7479]=['',5.1541917,0.7525743,6.16];
ybs[7480]=['46 Aql',5.1619605,0.21355,6.34];
ybs[7481]=['',5.2322679,-1.4189858,6.39];
ybs[7482]=['',5.1547182,0.7952872,5.06];
ybs[7483]=['',5.1685833,-0.2692566,5.49];
ybs[7484]=['χ Aql',5.1635127,0.2071529,5.27];
ybs[7485]=['',5.1985182,-1.264634,5.41];
ybs[7486]=['',5.1598799,0.7032967,6.23];
ybs[7487]=['',5.1508551,1.0567712,6.51];
ybs[7488]=['',5.1640711,0.5126741,6.49];
ybs[7489]=['',5.1636312,0.5666914,5.94];
ybs[7490]=['16 Cyg',5.1586895,0.8825649,5.96];
ybs[7491]=['',5.1589157,0.8824294,6.2];
ybs[7492]=['',5.1655072,0.5361849,6.05];
ybs[7493]=['10 Vul',5.168106,0.4505512,5.49];
ybs[7494]=['',5.1798476,-0.5561477,5.52];
ybs[7495]=['',5.1690064,0.4743522,6.28];
ybs[7496]=['',5.159429,0.9687516,6.48];
ybs[7497]=['ν Tel',5.1899002,-0.9829344,5.35];
ybs[7498]=['ψ Aql',5.1722067,0.2329292,6.26];
ybs[7499]=['',5.1684251,0.5969948,6.05];
ybs[7500]=['',5.1992924,-1.165313,6.45];
ybs[7501]=['',5.1676312,0.7298231,5.84];
ybs[7502]=['56 Sgr',5.1809378,-0.344132,4.86];
ybs[7503]=['',5.1783235,-0.0495628,6.48];
ybs[7504]=['15 Cyg',5.1701293,0.6527074,4.89];
ybs[7505]=['',5.1727746,0.511519,6.82];
ybs[7506]=['υ Aql',5.1771491,0.133637,5.91];
ybs[7507]=['',5.1718204,0.6013876,6.57];
ybs[7508]=['',5.1935421,-0.9222882,6.25];
ybs[7509]=['',5.164334,1.0133185,6.22];
ybs[7510]=['',5.1723383,0.7113927,6.34];
ybs[7511]=['',5.2040341,-1.144226,6.05];
ybs[7512]=['γ Aql',5.1796639,0.1860008,2.72];
ybs[7513]=['',5.1662446,0.9963237,6.27];
ybs[7514]=['',5.2005645,-1.0649298,6.21];
ybs[7515]=['δ Cyg',5.1728099,0.7884353,2.87];
ybs[7516]=['',5.1762209,0.6306673,6.43];
ybs[7517]=['',5.1771228,0.6118482,6.09];
ybs[7518]=['',5.2020447,-1.0323188,5.42];
ybs[7519]=['',5.1881279,-0.2383927,6.11];
ybs[7520]=['',5.1809742,0.4394346,6.62];
ybs[7521]=['17 Cyg',5.1796613,0.5894246,4.99];
ybs[7522]=['',5.1803771,0.5747794,6.18];
ybs[7523]=['δ Sge',5.1843637,0.3242527,3.82];
ybs[7524]=['',5.1989164,-0.8292392,5.94];
ybs[7525]=['',5.1935769,-0.5016775,6.05];
ybs[7526]=['',5.185993,0.4438055,5.95];
ybs[7527]=['',5.1923541,-0.1889498,6.04];
ybs[7528]=['',5.189466,0.1874265,6.44];
ybs[7529]=['',5.183997,0.6711078,5.77];
ybs[7530]=['π Aql',5.1902862,0.2070045,5.72];
ybs[7531]=['',5.1673788,1.2109054,5.92];
ybs[7532]=['ζ Sge',5.1912853,0.3348758,5];
ybs[7533]=['',5.1834473,0.8369185,6.12];
ybs[7534]=['',5.2098625,-0.95862,5.74];
ybs[7535]=['',5.2099646,-0.9587167,6.5];
ybs[7536]=['',5.1896624,0.617079,6.53];
ybs[7537]=['',5.1902225,0.5843694,6.44];
ybs[7538]=['',5.205563,-0.6951391,5.33];
ybs[7539]=['51 Aql',5.1999567,-0.187067,5.39];
ybs[7540]=['',5.1973311,0.1387145,6.51];
ybs[7541]=['',5.1926942,0.6764006,6.11];
ybs[7542]=['',5.1950581,0.4971682,6.38];
ybs[7543]=['α Aql',5.1994507,0.1555747,0.77];
ybs[7544]=['',5.2193624,-1.0668117,6.24];
ybs[7545]=['',5.2015026,-0.0421536,6.13];
ybs[7546]=['ο Aql',5.2004708,0.1825804,5.11];
ybs[7547]=['57 Sgr',5.2063841,-0.3315951,5.92];
ybs[7548]=['',5.201663,0.1688765,6.25];
ybs[7549]=['',5.1782163,1.195238,6.34];
ybs[7550]=['χ Cyg',5.1977582,0.5752517,4.23];
ybs[7551]=['12 Vul',5.2003027,0.3954136,4.95];
ybs[7552]=['19 Cyg',5.1975327,0.676626,5.12];
ybs[7553]=['',5.1976895,0.70939,5.69];
ybs[7554]=['',5.1985169,0.6609874,6.06];
ybs[7555]=['',5.2049399,0.2037638,6.13];
ybs[7556]=['η Aql',5.2070331,0.0183546,3.9];
ybs[7557]=['',5.2102125,-0.2540626,6.48];
ybs[7558]=['',5.2058547,0.1814685,6.54];
ybs[7559]=['',5.2044111,0.4369975,5.57];
ybs[7560]=['9 Sge',5.2060729,0.3266901,6.23];
ybs[7561]=['',5.2107987,-0.0535474,5.65];
ybs[7562]=['20 Cyg',5.1970111,0.9256068,5.03];
ybs[7563]=['',5.2004079,0.8276841,6.2];
ybs[7564]=['',5.21567,-0.4170345,6.18];
ybs[7565]=['',5.2380465,-1.2062901,5.75];
ybs[7566]=['',5.2108858,0.0776144,6.53];
ybs[7567]=['ι Sgr',5.2205587,-0.7299168,4.13];
ybs[7568]=['ε Dra',5.1840698,1.227177,3.83];
ybs[7569]=['',5.2050724,0.6366644,6.1];
ybs[7570]=['56 Aql',5.214542,-0.148832,5.79];
ybs[7571]=['',5.219423,-0.5759462,6.46];
ybs[7572]=['',5.2290053,-1.0101564,6.53];
ybs[7573]=['',5.2297144,-1.0271869,5.26];
ybs[7574]=['',5.2388156,-1.1992785,6.39];
ybs[7575]=['',5.2033172,0.8215846,5.62];
ybs[7576]=['ε Pav',5.2472463,-1.2716689,3.96];
ybs[7577]=['',5.2038549,0.837371,5.91];
ybs[7578]=['13 Vul',5.2107036,0.4210807,4.58];
ybs[7579]=['57 Aql',5.2166706,-0.1427736,5.71];
ybs[7580]=['57 Aql',5.2167072,-0.142948,6.49];
ybs[7581]=['ξ Aql',5.2145823,0.1484949,4.71];
ybs[7582]=['58 Aql',5.2169716,0.0055946,5.61];
ybs[7583]=['ω Sgr',5.2224834,-0.4581851,4.7];
ybs[7584]=['',5.2164647,0.1254399,6.15];
ybs[7585]=['',5.2196851,-0.1167103,6.51];
ybs[7586]=['',5.2077927,0.835205,6.29];
ybs[7587]=['',5.2153084,0.4252716,5.52];
ybs[7588]=['β Aql',5.2192841,0.1126401,3.71];
ybs[7589]=['μ1 Pav',5.2450199,-1.1676301,5.76];
ybs[7590]=['59 Sgr',5.2273384,-0.473372,4.52];
ybs[7591]=['',5.2309717,-0.6634041,6.55];
ybs[7592]=['',5.2161099,0.6465226,5.76];
ybs[7593]=['',5.2176899,0.5278279,6.57];
ybs[7594]=['23 Cyg',5.2082705,1.0047843,5.14];
ybs[7595]=['10 Sge',5.2220979,0.2911578,5.36];
ybs[7596]=['φ Aql',5.2231872,0.2002131,5.28];
ybs[7597]=['',5.2093703,1.0429264,6.06];
ybs[7598]=['μ2 Pav',5.2515071,-1.1675285,5.31];
ybs[7599]=['22 Cyg',5.2206615,0.6725446,4.94];
ybs[7600]=['61 Sgr',5.2313655,-0.2695355,5.02];
ybs[7601]=['η Cyg',5.2227421,0.6131482,3.89];
ybs[7602]=['',5.2262562,0.3673186,6.48];
ybs[7603]=['',5.2361451,-0.5321427,6.28];
ybs[7604]=['60 Sgr',5.236059,-0.4563518,4.83];
ybs[7605]=['ψ Cyg',5.2188973,0.9160548,4.92];
ybs[7606]=['',5.2245736,0.6335277,6.02];
ybs[7607]=['',5.2434502,-0.8604823,6.17];
ybs[7608]=['11 Sge',5.2296659,0.2938645,5.53];
ybs[7609]=['θ1 Sgr',5.2397922,-0.6148369,4.37];
ybs[7610]=['θ2 Sgr',5.2402862,-0.6047376,5.3];
ybs[7611]=['',5.2499849,-1.0354425,5.13];
ybs[7612]=['',5.2172836,1.0174805,6.09];
ybs[7613]=['',5.2431989,-0.750391,6.14];
ybs[7614]=['',5.2265615,0.7053898,5.45];
ybs[7615]=['',5.2421923,-0.6571723,5.95];
ybs[7616]=['',5.2448856,-0.7865121,5.81];
ybs[7617]=['',5.2423566,-0.5873878,5.66];
ybs[7618]=['',5.2239004,0.889247,6.43];
ybs[7619]=['',5.2196515,1.0278834,4.96];
ybs[7620]=['',5.2215675,0.9902015,6.12];
ybs[7621]=['γ Sge',5.2339531,0.3410482,3.47];
ybs[7622]=['',5.2371526,0.0248959,6.17];
ybs[7623]=['',5.2392427,-0.1729538,5.88];
ybs[7624]=['',5.2295531,0.7384295,6.43];
ybs[7625]=['',5.2474485,-0.7114783,6.29];
ybs[7626]=['',5.2330506,0.5416097,5.49];
ybs[7627]=['14 Vul',5.2356719,0.4040426,5.67];
ybs[7628]=['',5.2325171,0.6659105,6.32];
ybs[7629]=['',5.246616,-0.3959769,6.01];
ybs[7630]=['',5.2675648,-1.1740792,6.07];
ybs[7631]=['13 Sge',5.2396763,0.3065766,5.37];
ybs[7632]=['',5.2355131,0.7997238,5.92];
ybs[7633]=['25 Cyg',5.2384405,0.6473702,5.19];
ybs[7634]=['',5.2439631,0.1502254,5.91];
ybs[7635]=['63 Sgr',5.2488805,-0.2371436,5.71];
ybs[7636]=['62 Sgr',5.2522668,-0.4827552,4.58];
ybs[7637]=['',5.2347652,0.9093926,6.15];
ybs[7638]=['',5.2565623,-0.6613158,4.77];
ybs[7639]=['15 Vul',5.2439354,0.4852515,4.64];
ybs[7640]=['',5.2302901,1.1097215,5.96];
ybs[7641]=['',5.244267,0.6483581,6.2];
ybs[7642]=['',5.2468377,0.4337103,5.88];
ybs[7643]=['16 Vul',5.2480483,0.4361168,5.22];
ybs[7644]=['',5.2568182,-0.3934894,6.45];
ybs[7645]=['',5.2596897,-0.5586081,4.99];
ybs[7646]=['26 Cyg',5.2440795,0.8753529,5.05];
ybs[7647]=['',5.257632,-0.1294923,6.72];
ybs[7648]=['',5.2536951,0.3237693,5.96];
ybs[7649]=['',5.2795777,-1.1572,6.45];
ybs[7650]=['',5.2547533,0.2806757,5.67];
ybs[7651]=['δ Pav',5.2812274,-1.154182,3.56];
ybs[7652]=['',5.2300744,1.2289766,6.33];
ybs[7653]=['62 Aql',5.2590591,-0.011501,5.68];
ybs[7654]=['',5.2649926,-0.5750696,6.53];
ybs[7655]=['τ Aql',5.2577716,0.1279054,5.52];
ybs[7656]=['',5.2548721,0.5226709,5.71];
ybs[7657]=['',5.2624079,-0.2015627,6.34];
ybs[7658]=['15 Sge',5.2573577,0.2988067,5.8];
ybs[7659]=['ξ Tel',5.2740351,-0.9220433,4.94];
ybs[7660]=['',5.2750543,-0.9593143,6.26];
ybs[7661]=['65 Sgr',5.2639629,-0.2201628,6.55];
ybs[7662]=['64 Dra',5.2432137,1.1322021,5.27];
ybs[7663]=['',5.2609926,0.4059802,6.45];
ybs[7664]=['',5.2590627,0.5632026,5.64];
ybs[7665]=['η Sge',5.2618799,0.3497963,5.1];
ybs[7666]=['',5.2632421,0.2714184,6.34];
ybs[7667]=['',5.267071,-0.0702876,6.47];
ybs[7668]=['65 Dra',5.2470138,1.1289497,6.57];
ybs[7669]=['',5.2612755,0.6724588,6.19];
ybs[7670]=['',5.2578187,0.842648,6.16];
ybs[7671]=['ρ Dra',5.2485986,1.1854865,4.51];
ybs[7672]=['69 Dra',5.2320531,1.3356982,6.2];
ybs[7673]=['',5.2603552,0.9056534,6.14];
ybs[7674]=['17 Vul',5.2693321,0.4130463,5.07];
ybs[7675]=['27 Cyg',5.2666151,0.6287314,5.36];
ybs[7676]=['64 Aql',5.2749551,-0.0109351,5.99];
ybs[7677]=['',5.2906608,-1.0030563,6.37];
ybs[7678]=['',5.2610994,0.9842284,6.21];
ybs[7679]=['',5.2738632,0.1649588,6.43];
ybs[7680]=['',5.2773385,-0.1747212,6.18];
ybs[7681]=['',5.2576273,1.1159818,6.26];
ybs[7682]=['',5.2748487,0.2917489,6.42];
ybs[7683]=['',5.2651497,0.9288102,5.85];
ybs[7684]=['',5.3601551,-1.453028,6.17];
ybs[7685]=['',5.2724661,0.601697,6.11];
ybs[7686]=['',5.2773201,0.1881088,6.31];
ybs[7687]=['66 Dra',5.2613717,1.0829131,5.39];
ybs[7688]=['',5.269544,0.8775619,6.54];
ybs[7689]=['',5.2898029,-0.6291586,5.32];
ybs[7690]=['',5.2576223,1.1881807,6.28];
ybs[7691]=['θ Sge',5.2827389,0.3659562,6.48];
ybs[7692]=['',5.2953139,-0.7457293,6.22];
ybs[7693]=['',5.3057642,-1.105869,6.09];
ybs[7694]=['28 Cyg',5.2799566,0.6438865,4.93];
ybs[7695]=['',5.2888537,-0.153402,6.49];
ybs[7696]=['θ Aql',5.2892484,-0.0134113,3.23];
ybs[7697]=['18 Vul',5.285244,0.4704857,5.52];
ybs[7698]=['ξ1 Cap',5.2924207,-0.215361,6.34];
ybs[7699]=['',5.2875994,0.3697932,6.22];
ybs[7700]=['',5.3042068,-0.914403,5.65];
ybs[7701]=['ξ2 Cap',5.2944623,-0.2192851,5.85];
ybs[7702]=['',5.2888584,0.382725,6.26];
ybs[7703]=['',5.2947987,0.0160732,6.27];
ybs[7704]=['19 Vul',5.2906674,0.4688306,5.49];
ybs[7705]=['20 Vul',5.2916016,0.4630723,5.92];
ybs[7706]=['66 Aql',5.2976597,-0.0166815,5.47];
ybs[7707]=['',5.2909682,0.8340951,6.92];
ybs[7708]=['',5.3073354,-0.4708612,5.73];
ybs[7709]=['',5.2989379,0.4239873,6.56];
ybs[7710]=['ρ Aql',5.3018099,0.2661891,4.95];
ybs[7711]=['',5.3098431,-0.5227377,6.3];
ybs[7712]=['',5.2927854,0.8991404,6.01];
ybs[7713]=['68 Dra',5.2877396,1.0844006,5.75];
ybs[7714]=['',5.3124753,-0.6352933,6.39];
ybs[7715]=['',5.312632,-0.6133601,6.53];
ybs[7716]=['30 Cyg',5.2964194,0.8180263,4.83];
ybs[7717]=['21 Vul',5.3012665,0.5017597,5.18];
ybs[7718]=['',5.3257821,-1.1026169,6.27];
ybs[7719]=['',5.2983958,0.7580479,6.14];
ybs[7720]=['',5.3002983,0.6398237,6.45];
ybs[7721]=['31 Cyg',5.2978719,0.8167291,3.79];
ybs[7722]=['29 Cyg',5.3022555,0.6433366,4.97];
ybs[7723]=['',5.3012785,0.7357892,6.71];
ybs[7724]=['3 Cap',5.3116793,-0.2143641,6.32];
ybs[7725]=['',5.3058307,0.4476123,4.78];
ybs[7726]=['33 Cyg',5.29621,0.9882301,4.3];
ybs[7727]=['22 Vul',5.3069331,0.4112528,5.15];
ybs[7728]=['',5.2961127,1.0593136,5.79];
ybs[7729]=['',5.3061317,0.5896391,5.66];
ybs[7730]=['23 Vul',5.3079593,0.4864005,4.52];
ybs[7731]=['',5.3240455,-0.8317385,6.31];
ybs[7732]=['18 Sge',5.3105894,0.3779222,6.13];
ybs[7733]=['α1 Cap',5.3172167,-0.2173477,4.24];
ybs[7734]=['4 Cap',5.319102,-0.3796899,5.87];
ybs[7735]=['',5.3256294,-0.8294576,6.13];
ybs[7736]=['κ Cep',5.2720597,1.3572243,4.39];
ybs[7737]=['32 Cyg',5.3058564,0.8337235,3.98];
ybs[7738]=['',5.3088244,0.6798475,6.27];
ybs[7739]=['24 Vul',5.3124935,0.4315501,5.32];
ybs[7740]=['α2 Cap',5.3189917,-0.2179803,3.57];
ybs[7741]=['',5.3068115,0.8776781,6.31];
ybs[7742]=['',5.3083194,0.7964641,5.91];
ybs[7743]=['',5.3107022,0.6477118,6.48];
ybs[7744]=['',5.3315002,-0.9598356,6.27];
ybs[7745]=['',5.3125409,0.7054606,5.24];
ybs[7746]=['',5.3155858,0.5096919,6.22];
ybs[7747]=['σ Cap',5.3250028,-0.332708,5.28];
ybs[7748]=['',5.3148978,0.7466004,6.29];
ybs[7749]=['34 Cyg',5.3164133,0.6647658,4.81];
ybs[7750]=['',5.3299863,-0.5086063,6.3];
ybs[7751]=['',5.3319365,-0.621638,6.46];
ybs[7752]=['',5.3361591,-0.8716655,6.27];
ybs[7753]=['',5.3177437,0.7118771,5.84];
ybs[7754]=['',5.3259788,-0.0178489,6.06];
ybs[7755]=['36 Cyg',5.3194664,0.6467399,5.58];
ybs[7756]=['35 Cyg',5.3203058,0.6115339,5.17];
ybs[7757]=['',5.3245999,0.2316539,6.21];
ybs[7758]=['',5.32922,-0.1100512,6.63];
ybs[7759]=['ν Cap',5.3303757,-0.2217072,4.76];
ybs[7760]=['',5.3268464,0.237436,5.95];
ybs[7761]=['',5.3309221,-0.2570639,6.1];
ybs[7762]=['β Cap',5.331947,-0.2569995,3.08];
ybs[7763]=['',5.320584,0.8094549,6.45];
ybs[7764]=['',5.3282966,0.2552598,6.13];
ybs[7765]=['κ1 Sgr',5.3391375,-0.7329123,5.59];
ybs[7766]=['',5.3282773,0.3115273,5.8];
ybs[7767]=['',5.3182151,0.9678309,5.76];
ybs[7768]=['',5.3252616,0.6490605,6.57];
ybs[7769]=['',5.3130674,1.1677812,5.93];
ybs[7770]=['',5.3271297,0.6886967,6.23];
ybs[7771]=['',5.3934841,-1.4120432,5.77];
ybs[7772]=['',5.3253984,0.8184452,6.5];
ybs[7773]=['κ2 Sgr',5.3453805,-0.7394148,5.64];
ybs[7774]=['',5.3405496,-0.1675103,6.3];
ybs[7775]=['25 Vul',5.3355123,0.4276552,5.54];
ybs[7776]=['α Pav',5.3538924,-0.9891991,1.94];
ybs[7777]=['',5.3274918,0.9364084,6.18];
ybs[7778]=['71 Dra',5.3228199,1.0875723,5.72];
ybs[7779]=['',5.3393391,0.2549648,6.17];
ybs[7780]=['',5.3409063,0.0942511,5.31];
ybs[7781]=['',5.3349037,0.7188678,6.39];
ybs[7782]=['γ Cyg',5.3357193,0.7036021,2.2];
ybs[7783]=['',5.3377753,0.5466705,6.09];
ybs[7784]=['',5.3348641,0.8002632,5.58];
ybs[7785]=['',5.3536448,-0.7110176,6.09];
ybs[7786]=['',5.337987,0.7170346,5.93];
ybs[7787]=['',5.3516797,-0.4992582,5.85];
ybs[7788]=['',5.3386439,0.7511955,6.2];
ybs[7789]=['',5.3473218,0.0196567,6.15];
ybs[7790]=['',5.3239942,1.2031638,5.55];
ybs[7791]=['',5.3295146,1.1176498,5.69];
ybs[7792]=['39 Cyg',5.343151,0.5628223,4.43];
ybs[7793]=['',5.3424309,0.6550863,5.9];
ybs[7794]=['',5.3582411,-0.6517864,6.25];
ybs[7795]=['',5.3521411,-0.0478617,6.11];
ybs[7796]=['',5.3519411,0.1765295,6.33];
ybs[7797]=['',5.351393,0.3746819,5.66];
ybs[7798]=['',5.4155393,-1.4176673,5.91];
ybs[7799]=['',5.3529409,0.3477236,6.41];
ybs[7800]=['π Cap',5.3595535,-0.3168312,5.25];
ybs[7801]=['',5.3451172,0.9356619,6.51];
ybs[7802]=['',5.354609,0.3032296,6.22];
ybs[7803]=['',5.3664372,-0.6202332,6.1];
ybs[7804]=['',5.3470172,1.0412233,6.44];
ybs[7805]=['',5.3656255,-0.2737133,6.41];
ybs[7806]=['',5.3624179,0.1482888,6.25];
ybs[7807]=['68 Aql',5.3639639,-0.0575758,6.13];
ybs[7808]=['ρ Cap',5.3662605,-0.3098748,4.78];
ybs[7809]=['',5.3573409,0.6001723,6.39];
ybs[7810]=['',5.3632258,0.052287,6.21];
ybs[7811]=['',5.3692687,-0.3897729,6.16];
ybs[7812]=['40 Cyg',5.3591259,0.6719321,5.62];
ybs[7813]=['',5.3529823,0.9895502,6.36];
ybs[7814]=['43 Cyg',5.3562894,0.862921,5.69];
ybs[7815]=['ο Cap',5.3707078,-0.3233611,6.74];
ybs[7816]=['ο Cap',5.3708095,-0.3233029,5.94];
ybs[7817]=['69 Aql',5.3693411,-0.0493266,4.91];
ybs[7818]=['',5.3756782,-0.5070651,6.39];
ybs[7819]=['',5.3674775,0.3516316,6.55];
ybs[7820]=['41 Cyg',5.367369,0.531066,4.01];
ybs[7821]=['42 Cyg',5.3669224,0.6372881,5.88];
ybs[7822]=['1 Del',5.3718501,0.1912076,6.08];
ybs[7823]=['',5.3758373,-0.2617393,6.12];
ybs[7824]=['',5.399918,-1.2138689,6.11];
ybs[7825]=['',5.3745204,0.3606827,6.18];
ybs[7826]=['',5.3758419,0.1975786,7.11];
ybs[7827]=['',5.3694001,0.8026424,6.41];
ybs[7828]=['',5.3839576,-0.434298,6.36];
ybs[7829]=['',5.3663988,0.9796054,5.91];
ybs[7830]=['ω1 Cyg',5.3695063,0.8554049,4.95];
ybs[7831]=['',5.3814824,-0.1709211,5.65];
ybs[7832]=['ν Mic',5.3892027,-0.7758912,5.11];
ybs[7833]=['44 Cyg',5.3740953,0.6456948,6.19];
ybs[7834]=['φ1 Pav',5.3974864,-1.0562777,4.76];
ybs[7835]=['',5.3787487,0.4514214,6.34];
ybs[7836]=['θ Cep',5.3663228,1.1004892,4.22];
ybs[7837]=['ω2 Cyg',5.3749741,0.8601004,5.44];
ybs[7838]=['ε Del',5.3845547,0.1983368,4.03];
ybs[7839]=['',5.3933268,-0.6637239,6.44];
ybs[7840]=['',5.3749602,0.9140215,6.18];
ybs[7841]=['',5.3894285,-0.2384161,6.13];
ybs[7842]=['',5.3924798,-0.5307986,6.4];
ybs[7843]=['',5.3875523,0.1766357,6.56];
ybs[7844]=['η Del',5.3877279,0.2284286,5.38];
ybs[7845]=['ρ Pav',5.4063251,-1.0728177,4.88];
ybs[7846]=['',5.3764973,0.9920447,6.14];
ybs[7847]=['',5.3820768,0.7548906,6.6];
ybs[7848]=['',5.3884712,0.367324,6.48];
ybs[7849]=['μ1 Oct',5.4285943,-1.328491,6];
ybs[7850]=['μ2 Oct',5.4268972,-1.3140068,6.55];
ybs[7851]=['',5.3953505,-0.2873599,6.19];
ybs[7852]=['47 Cyg',5.3868917,0.6163031,4.61];
ybs[7853]=['',5.3862212,0.730122,6.49];
ybs[7854]=['',5.3665568,1.2669514,6.27];
ybs[7855]=['α Ind',5.4052358,-0.8243078,3.11];
ybs[7856]=['',5.386471,0.8160218,5.78];
ybs[7857]=['ζ Del',5.3936158,0.257181,4.68];
ybs[7858]=['',5.4163462,-1.0968514,6.22];
ybs[7859]=['70 Aql',5.4002095,-0.043429,4.89];
ybs[7860]=['26 Vul',5.3969436,0.4528079,6.41];
ybs[7861]=['φ2 Pav',5.4169209,-1.0556802,5.12];
ybs[7862]=['',5.3902296,0.9060908,6.11];
ybs[7863]=['',5.4057462,-0.437149,6.36];
ybs[7864]=['',5.4026721,0.0027722,6.22];
ybs[7865]=['73 Dra',5.372369,1.30925,5.2];
ybs[7866]=['27 Vul',5.4010307,0.4629266,5.59];
ybs[7867]=['υ Pav',5.4259481,-1.1640879,5.15];
ybs[7868]=['β Del',5.4034015,0.2558171,3.63];
ybs[7869]=['ι Del',5.4046504,0.1996627,5.43];
ybs[7870]=['71 Aql',5.4072083,-0.0182046,4.32];
ybs[7871]=['48 Cyg',5.4028492,0.5521251,6.32];
ybs[7872]=['',5.4048814,0.3199405,6.25];
ybs[7873]=['',5.402909,0.5512427,6.49];
ybs[7874]=['',5.402021,0.6700405,6.2];
ybs[7875]=['τ Cap',5.411611,-0.2599174,5.22];
ybs[7876]=['',5.411085,-0.0410197,6.22];
ybs[7877]=['29 Vul',5.4074756,0.371116,4.82];
ybs[7878]=['θ Del',5.40859,0.2334788,5.72];
ybs[7879]=['',5.4167303,-0.5823995,5.47];
ybs[7880]=['28 Vul',5.4074411,0.4219923,5.04];
ybs[7881]=['',5.4076858,0.4143908,5.91];
ybs[7882]=['κ Del',5.4104053,0.1771263,5.05];
ybs[7883]=['1 Aqr',5.4118698,0.0095813,5.16];
ybs[7884]=['',5.4158565,-0.4138356,6.37];
ybs[7885]=['',5.4100635,0.2775163,5.97];
ybs[7886]=['υ Cap',5.4150774,-0.3154824,5.1];
ybs[7887]=['75 Dra',5.3540907,1.4221169,5.46];
ybs[7888]=['',5.4177022,-0.4639436,6.51];
ybs[7889]=['',5.4103186,0.3818729,6.08];
ybs[7890]=['',5.4092695,0.5305252,5.68];
ybs[7891]=['',5.417172,-0.2803209,5.8];
ybs[7892]=['α Del',5.412484,0.2788091,3.77];
ybs[7893]=['',5.4135783,0.1974393,6.42];
ybs[7894]=['74 Dra',5.3597707,1.4163405,5.96];
ybs[7895]=['',5.4213156,-0.550391,5.76];
ybs[7896]=['',5.4211801,-0.4526817,6.28];
ybs[7897]=['',5.4113751,0.7093372,6.06];
ybs[7898]=['',5.4104144,0.7981297,6.58];
ybs[7899]=['β Pav',5.4389555,-1.154336,3.42];
ybs[7900]=['',5.4172499,0.3490357,6.45];
ybs[7901]=['',5.4280762,-0.6893155,6.29];
ybs[7902]=['',5.4081461,0.9785602,6.48];
ybs[7903]=['',5.4163222,0.5212987,6.08];
ybs[7904]=['10 Del',5.419638,0.2556248,5.99];
ybs[7905]=['',5.4134255,0.759591,5.95];
ybs[7906]=['η Ind',5.4336866,-0.9050749,4.51];
ybs[7907]=['49 Cyg',5.4181661,0.5649683,5.51];
ybs[7908]=['',5.4177727,0.6832141,6.51];
ybs[7909]=['',5.4226215,0.3069123,6.22];
ybs[7910]=['α Cyg',5.4193725,0.7913927,1.25];
ybs[7911]=['',5.4134342,1.0571121,6.01];
ybs[7912]=['',5.4217534,0.7292038,5.67];
ybs[7913]=['',5.4238559,0.6199343,6.66];
ybs[7914]=['δ Del',5.4291774,0.2642135,4.43];
ybs[7915]=['51 Cyg',5.4225185,0.8797058,5.39];
ybs[7916]=['',5.3544543,1.4605653,6.19];
ybs[7917]=['',5.4378663,-0.4744284,6.5];
ybs[7918]=['',5.4283637,0.6222381,6.47];
ybs[7919]=['',5.4431096,-0.6830228,5.5];
ybs[7920]=['σ Pav',5.4582417,-1.1992249,5.41];
ybs[7921]=['',5.442892,-0.6292863,6.49];
ybs[7922]=['ψ Cap',5.4416275,-0.4399291,4.14];
ybs[7923]=['17 Cap',5.4418412,-0.3743625,5.93];
ybs[7924]=['',5.4237573,1.0588029,6.15];
ybs[7925]=['30 Vul',5.4351016,0.442177,4.91];
ybs[7926]=['',5.4265145,0.9979428,6.32];
ybs[7927]=['',5.4378846,0.3168609,6.38];
ybs[7928]=['52 Cyg',5.4383802,0.5372872,4.22];
ybs[7929]=['ι Mic',5.4526776,-0.7666024,5.11];
ybs[7930]=['',5.4315514,0.9870215,5.78];
ybs[7931]=['4 Cep',5.4253377,1.1645044,5.58];
ybs[7932]=['',5.4452935,-0.0422654,6.27];
ybs[7933]=['γ1 Del',5.443061,0.2825572,5.14];
ybs[7934]=['γ2 Del',5.4431192,0.2825525,4.27];
ybs[7935]=['ε Cyg',5.4406805,0.5940231,2.46];
ybs[7936]=['ε Aqr',5.4481405,-0.164595,3.77];
ybs[7937]=['3 Aqr',5.4483022,-0.0866124,4.42];
ybs[7938]=['ζ Ind',5.4571197,-0.8056632,4.89];
ybs[7939]=['13 Del',5.4483563,0.1060043,5.58];
ybs[7940]=['',5.4483799,0.0588513,6.4];
ybs[7941]=['',5.4357742,1.0060801,4.51];
ybs[7942]=['',5.4448957,0.6010776,4.92];
ybs[7943]=['η Cep',5.435129,1.080416,3.43];
ybs[7944]=['',5.4420737,0.8132627,6.3];
ybs[7945]=['',5.467591,-1.0884329,6.28];
ybs[7946]=['',5.4676273,-1.0884328,6.59];
ybs[7947]=['',5.4555891,-0.4488223,5.86];
ybs[7948]=['',5.440468,0.9260722,6.33];
ybs[7949]=['λ Cyg',5.4458223,0.6380216,4.53];
ybs[7950]=['',5.4556026,-0.3136368,6.21];
ybs[7951]=['α Mic',5.4587651,-0.5884158,4.9];
ybs[7952]=['',5.4451973,0.7966521,6.4];
ybs[7953]=['',5.4308176,1.2185201,6.41];
ybs[7954]=['ι Ind',5.4661865,-0.8995749,5.05];
ybs[7955]=['',5.4471664,0.8359633,5.57];
ybs[7956]=['',5.4622725,-0.5582996,6.36];
ybs[7957]=['',5.4634507,-0.6605551,5.52];
ybs[7958]=['',5.4471748,0.9158172,6.27];
ybs[7959]=['15 Del',5.4561733,0.2201055,5.98];
ybs[7960]=['14 Del',5.4570429,0.1384057,6.33];
ybs[7961]=['',5.4578798,0.0979248,6.21];
ybs[7962]=['',5.4613841,-0.2177964,5.88];
ybs[7963]=['55 Cyg',5.4521261,0.8059886,4.84];
ybs[7964]=['',5.4508381,0.9071533,6.29];
ybs[7965]=['β Mic',5.4675172,-0.5778893,6.04];
ybs[7966]=['ω Cap',5.4666456,-0.4686665,4.11];
ybs[7967]=['',5.4603504,0.3162106,6.52];
ybs[7968]=['4 Aqr',5.4644216,-0.0970399,5.99];
ybs[7969]=['',5.456353,0.8155399,6.33];
ybs[7970]=['56 Cyg',5.4572111,0.7701333,5.04];
ybs[7971]=['5 Aqr',5.4675456,-0.0949513,5.55];
ybs[7972]=['β Ind',5.4810492,-1.0190392,3.65];
ybs[7973]=['',5.4751066,-0.6936438,5.35];
ybs[7974]=['',5.4638119,0.494224,5.77];
ybs[7975]=['',5.4717884,-0.4139246,6.33];
ybs[7976]=['μ Aqr',5.4698398,-0.1556229,4.73];
ybs[7977]=['',5.4737075,-0.5349707,6.35];
ybs[7978]=['',5.4795652,-0.8841897,6.24];
ybs[7979]=['',5.4523847,1.1188933,6.45];
ybs[7980]=['',5.4718169,-0.2008294,6.38];
ybs[7981]=['31 Vul',5.4667181,0.4740933,4.59];
ybs[7982]=['',5.4660131,0.5744878,6.44];
ybs[7983]=['',5.4766663,-0.486219,6.41];
ybs[7984]=['',5.4755523,-0.1190754,6.44];
ybs[7985]=['',5.470991,0.5186482,6.34];
ybs[7986]=['19 Cap',5.4794001,-0.3116389,5.78];
ybs[7987]=['57 Cyg',5.4710227,0.7758713,4.78];
ybs[7988]=['76 Dra',5.4158116,1.4415434,5.75];
ybs[7989]=['',5.4712653,0.7897421,5.45];
ybs[7990]=['',5.4719481,0.7413682,6.66];
ybs[7991]=['',5.4742624,0.5847712,5.47];
ybs[7992]=['',5.4805018,-0.0227902,6.55];
ybs[7993]=['',5.4764805,0.4989763,6.56];
ybs[7994]=['32 Vul',5.4773084,0.4908713,5.01];
ybs[7995]=['',5.4760832,0.7115766,6.7];
ybs[7996]=['',5.4827462,0.0802937,6.05];
ybs[7997]=['17 Del',5.4822535,0.2406648,5.17];
ybs[7998]=['16 Del',5.4824177,0.2205451,5.58];
ybs[7999]=['',5.4882873,-0.4577704,5.7];
ybs[8000]=['',5.4856534,-0.0609727,6.57];
ybs[8001]=['7 Aqr',5.4883814,-0.1680648,5.51];
ybs[8002]=['',5.4397598,1.4070345,5.39];
ybs[8003]=['',5.489371,0.0092814,6.05];
ybs[8004]=['',5.4919107,-0.2786127,5.87];
ybs[8005]=['',5.5110802,-1.18927,6.37];
ybs[8006]=['',5.4821737,0.8287782,5.67];
ybs[8007]=['α Oct',5.5273523,-1.3430883,5.15];
ybs[8008]=['',5.4836722,0.8926107,6.63];
ybs[8009]=['',5.4855539,0.7852752,5.96];
ybs[8010]=['',5.496332,-0.2515741,6.01];
ybs[8011]=['',5.4846085,0.8865663,5.81];
ybs[8012]=['',5.4847166,0.8598144,5.9];
ybs[8013]=['',5.5047549,-0.8935404,5.76];
ybs[8014]=['ν Cyg',5.4883029,0.7196928,3.94];
ybs[8015]=['',5.4836098,0.9940582,6.23];
ybs[8016]=['18 Del',5.4946254,0.1903755,5.48];
ybs[8017]=['',5.5025534,-0.6293774,6.11];
ybs[8018]=['33 Vul',5.4936691,0.3908547,5.31];
ybs[8019]=['20 Cap',5.5003765,-0.3310254,6.25];
ybs[8020]=['ε Equ',5.4975661,0.0761375,5.23];
ybs[8021]=['',5.4932033,0.7773722,5.55];
ybs[8022]=['',5.4941267,0.7331922,6.16];
ybs[8023]=['',5.5006676,0.2948408,6.66];
ybs[8024]=['',5.5018253,0.1323907,5.99];
ybs[8025]=['γ Mic',5.5080884,-0.5617925,4.67];
ybs[8026]=['',5.4936967,0.8819279,5.61];
ybs[8027]=['11 Aqr',5.5042569,-0.0813512,6.21];
ybs[8028]=['',5.5124499,-0.7493086,6.64];
ybs[8029]=['',5.4738229,1.3263248,6.05];
ybs[8030]=['',5.503294,0.3385693,5.65];
ybs[8031]=['',5.5099608,-0.4679496,6.05];
ybs[8032]=['',5.5133489,-0.671267,5.94];
ybs[8033]=['59 Cyg',5.4996303,0.8306029,4.74];
ybs[8034]=['ζ Mic',5.5155901,-0.6730291,5.3];
ybs[8035]=['',5.4971422,1.0386,5.51];
ybs[8036]=['',5.5161534,-0.4827921,6.25];
ybs[8037]=['',5.5061326,0.6299849,5.97];
ybs[8038]=['',5.5447643,-1.3289065,6.58];
ybs[8039]=['60 Cyg',5.5056113,0.8067815,5.37];
ybs[8040]=['',5.5147799,-0.0149191,6.5];
ybs[8041]=['μ Ind',5.5262035,-0.9539374,5.16];
ybs[8042]=['',5.5149764,0.0279581,6.25];
ybs[8043]=['',5.5146094,0.2583073,6.31];
ybs[8044]=['12 Aqr',5.519608,-0.1004103,7.31];
ybs[8045]=['12 Aqr',5.5196152,-0.1004055,5.89];
ybs[8046]=['η Cap',5.5213451,-0.3453073,4.84];
ybs[8047]=['',5.5463777,-1.2758557,5.68];
ybs[8048]=['',5.510999,0.782969,6.19];
ybs[8049]=['',5.5141911,0.6759209,6.07];
ybs[8050]=['',5.5127312,0.8014327,6.48];
ybs[8051]=['',5.5092775,0.9902881,5.83];
ybs[8052]=['3 Equ',5.5215634,0.09727,5.61];
ybs[8053]=['',5.5221253,0.0525756,6.42];
ybs[8054]=['',5.5224082,0.0408435,6.33];
ybs[8055]=['η Mic',5.5307605,-0.7210854,5.53];
ybs[8056]=['δ Mic',5.5286443,-0.5245442,5.68];
ybs[8057]=['',5.5175293,0.7277708,6.33];
ybs[8058]=['',5.5152524,0.8800289,6.37];
ybs[8059]=['',5.5413523,-1.1145196,5.76];
ybs[8060]=['',5.5166775,0.8191185,6.32];
ybs[8061]=['θ Cap',5.5280095,-0.299533,4.07];
ybs[8062]=['',5.5304137,-0.5632303,5.18];
ybs[8063]=['4 Equ',5.5253359,0.1052253,5.94];
ybs[8064]=['',5.5166667,0.9312416,5.9];
ybs[8065]=['ξ Cyg',5.5220863,0.767914,3.72];
ybs[8066]=['24 Cap',5.5333449,-0.4351923,4.5];
ybs[8067]=['',5.5548093,-1.2648699,6.2];
ybs[8068]=['',5.5289964,0.4711576,6.12];
ybs[8069]=['',5.5358508,-0.3034073,6.17];
ybs[8070]=['',5.5293775,0.545514,5.82];
ybs[8071]=['61 Cyg',5.5309162,0.6774822,5.21];
ybs[8072]=['61 Cyg',5.5309672,0.6774386,6.03];
ybs[8073]=['χ Cap',5.5394973,-0.3686493,5.3];
ybs[8074]=['',5.5343615,0.2745378,6.34];
ybs[8075]=['63 Cyg',5.5292255,0.8328584,4.55];
ybs[8076]=['',5.5385157,0.123237,6.15];
ybs[8077]=['27 Cap',5.5438003,-0.3575228,6.25];
ybs[8078]=['ο Pav',5.5630233,-1.2226621,5.02];
ybs[8079]=['ν Aqr',5.5437969,-0.197219,4.51];
ybs[8080]=['',5.5387633,0.52844,5.59];
ybs[8081]=['',5.5451505,0.0526267,6.45];
ybs[8082]=['',5.5489319,-0.1619962,6.27];
ybs[8083]=['γ Equ',5.5466194,0.1780884,4.69];
ybs[8084]=['6 Equ',5.5473994,0.1766446,6.07];
ybs[8085]=['',5.5260904,1.2479583,5.87];
ybs[8086]=['',5.5560085,-0.7015668,5.83];
ybs[8087]=['',5.5471875,0.3931671,6.68];
ybs[8088]=['',5.5529933,-0.2513234,6.48];
ybs[8089]=['',5.5440631,0.7954236,6.63];
ybs[8090]=['',5.5595763,-0.6868293,5.26];
ybs[8091]=['',5.5491416,0.6348008,6.54];
ybs[8092]=['',5.5449141,0.9361129,5.73];
ybs[8093]=['',5.5463295,0.8336393,6.46];
ybs[8094]=['',5.5606292,-0.6344439,5.96];
ybs[8095]=['',5.5407813,1.1059678,6.54];
ybs[8096]=['',5.5602711,-0.4807776,5.42];
ybs[8097]=['',5.5855979,-1.3137486,6.63];
ybs[8098]=['',5.5200457,1.3647926,5.91];
ybs[8099]=['',5.5403725,1.1966328,7.33];
ybs[8100]=['',5.5719782,-0.9283305,5.75];
ybs[8101]=['ζ Cyg',5.5574878,0.5288299,3.2];
ybs[8102]=['',5.5601945,0.2802202,6.27];
ybs[8103]=['',5.5691722,-0.7056875,6.21];
ybs[8104]=['',5.5642111,-0.1838198,6.77];
ybs[8105]=['',5.5512301,1.0482238,5.64];
ybs[8106]=['',5.5594966,0.6406496,6.05];
ybs[8107]=['',5.5654754,0.0028884,6.38];
ybs[8108]=['',5.5679789,-0.301446,6.04];
ybs[8109]=['δ Equ',5.5646901,0.1759321,4.49];
ybs[8110]=['',5.5713683,-0.6307134,6.12];
ybs[8111]=['',5.5825842,-1.1276115,6.31];
ybs[8112]=['',5.5628891,0.523149,6.17];
ybs[8113]=['φ Cap',5.5703261,-0.3591556,5.24];
ybs[8114]=['29 Cap',5.5707171,-0.2635062,5.28];
ybs[8115]=['',5.6526982,-1.4788477,6.45];
ybs[8116]=['τ Cyg',5.5653636,0.6652993,3.72];
ybs[8117]=['α Equ',5.57064,0.0928756,3.92];
ybs[8118]=['',5.574424,-0.0267723,6.48];
ybs[8119]=['',5.5591646,1.1253329,6.39];
ybs[8120]=['',5.5771192,-0.2304687,6.4];
ybs[8121]=['ε Mic',5.5806635,-0.5602208,4.71];
ybs[8122]=['',5.5685868,0.8385804,6.46];
ybs[8123]=['30 Cap',5.580401,-0.3126071,5.43];
ybs[8124]=['',5.5726559,0.7387132,6.43];
ybs[8125]=['31 Cap',5.5817275,-0.3034767,7.05];
ybs[8126]=['θ Ind',5.589854,-0.9315687,4.39];
ybs[8127]=['15 Aqr',5.581134,-0.0775831,5.82];
ybs[8128]=['',5.5847974,-0.5007538,6.4];
ybs[8129]=['σ Cyg',5.576786,0.6888595,4.23];
ybs[8130]=['',5.5765431,0.7462564,6.19];
ybs[8131]=['',5.5907463,-0.7844799,6];
ybs[8132]=['υ Cyg',5.5791138,0.610361,4.43];
ybs[8133]=['',5.5744953,0.9437239,6.13];
ybs[8134]=['',5.5884743,-0.4586437,6.56];
ybs[8135]=['',5.583804,0.1968343,5.96];
ybs[8136]=['',5.5752858,0.9751504,5.98];
ybs[8137]=['θ1 Mic',5.5932193,-0.7109551,4.82];
ybs[8138]=['',5.5958085,-0.8702671,6.38];
ybs[8139]=['',5.5754275,1.0242575,6.42];
ybs[8140]=['68 Cyg',5.5811569,0.7682963,5];
ybs[8141]=['',5.5833056,0.7175968,6.15];
ybs[8142]=['',5.6105107,-1.2157644,6.41];
ybs[8143]=['',5.5853545,0.6686716,5.83];
ybs[8144]=['',5.5895436,0.3857386,6.29];
ybs[8145]=['',5.6153102,-1.2518054,6.09];
ybs[8146]=['16 Aqr',5.5937289,-0.0782771,5.87];
ybs[8147]=['',5.5854148,0.865419,5.76];
ybs[8148]=['α Cep',5.5806472,1.0936208,2.44];
ybs[8149]=['9 Equ',5.5935388,0.1296691,5.82];
ybs[8150]=['',5.5839533,1.0244753,5.66];
ybs[8151]=['',5.5931699,0.4176724,5.57];
ybs[8152]=['',5.5919188,0.5677161,5.68];
ybs[8153]=['ι Cap',5.5990952,-0.2925009,4.28];
ybs[8154]=['',5.5653585,1.3453982,5.95];
ybs[8155]=['',5.5942363,0.5705113,6.04];
ybs[8156]=['',5.5925092,0.7054719,6.4];
ybs[8157]=['6 Cep',5.5838784,1.1335295,5.18];
ybs[8158]=['',5.6025407,-0.3943275,5.6];
ybs[8159]=['1 Peg',5.5976826,0.3469672,4.08];
ybs[8160]=['',5.5526147,1.4190139,6.15];
ybs[8161]=['17 Aqr',5.6019635,-0.1613363,5.99];
ybs[8162]=['',5.6574234,-1.4417192,6.38];
ybs[8163]=['',5.6090581,-0.8122594,6.31];
ybs[8164]=['β Equ',5.6014627,0.1201946,5.16];
ybs[8165]=['',5.589456,1.061715,6.11];
ybs[8166]=['θ2 Mic',5.6091483,-0.7143753,5.77];
ybs[8167]=['γ Pav',5.6193713,-1.1395178,4.22];
ybs[8168]=['',5.6001112,0.5303216,6.05];
ybs[8169]=['33 Cap',5.6075239,-0.3626106,5.41];
ybs[8170]=['',5.6074412,-0.3956847,6.38];
ybs[8171]=['',5.5964618,0.8633123,5.69];
ybs[8172]=['',5.6002439,0.6756108,6.63];
ybs[8173]=['18 Aqr',5.6074997,-0.2234399,5.49];
ybs[8174]=['γ Ind',5.6177506,-0.9526719,6.12];
ybs[8175]=['',5.6029091,0.6541898,6.58];
ybs[8176]=['',5.6058506,0.4249873,5.71];
ybs[8177]=['',5.6080175,0.1788981,6.35];
ybs[8178]=['20 Aqr',5.6102425,-0.0579845,6.36];
ybs[8179]=['',5.6047525,0.6532269,6.47];
ybs[8180]=['',5.6064529,0.4431055,6.15];
ybs[8181]=['19 Aqr',5.6119179,-0.1688161,5.7];
ybs[8182]=['',5.6299103,-1.211749,5.34];
ybs[8183]=['',5.6076132,0.429435,6.32];
ybs[8184]=['',5.6083688,0.4581561,5.68];
ybs[8185]=['21 Aqr',5.6120852,-0.0607459,5.49];
ybs[8186]=['',5.6176291,-0.6589134,5.63];
ybs[8187]=['',5.6529095,-1.3955781,6.47];
ybs[8188]=['',5.6205886,-0.7412607,5.51];
ybs[8189]=['',5.6145183,0.01066,6.46];
ybs[8190]=['ζ Cap',5.6184888,-0.3898164,3.74];
ybs[8191]=['',5.6171695,0.0205918,6.13];
ybs[8192]=['',5.609218,0.8621818,6.58];
ybs[8193]=['35 Cap',5.6209925,-0.3686031,5.78];
ybs[8194]=['',5.6110816,0.81665,5.6];
ybs[8195]=['69 Cyg',5.6134131,0.6413,5.94];
ybs[8196]=['',5.6167174,0.3395019,6.07];
ybs[8197]=['',5.6296694,-0.9359963,6.39];
ybs[8198]=['',5.6250972,-0.2005623,6.61];
ybs[8199]=['36 Cap',5.627445,-0.3792623,4.51];
ybs[8200]=['5 PsA',5.6291418,-0.5438694,6.5];
ybs[8201]=['70 Cyg',5.620272,0.6491467,5.31];
ybs[8202]=['',5.617695,0.8536677,5.31];
ybs[8203]=['35 Vul',5.6218783,0.4832014,5.41];
ybs[8204]=['',5.617026,0.9245905,6.03];
ybs[8205]=['',5.6255264,0.1443832,6.4];
ybs[8206]=['',5.6238183,0.5637795,5.8];
ybs[8207]=['',5.6278933,0.3138621,6.44];
ybs[8208]=['',5.6329289,-0.3328406,6.57];
ybs[8209]=['',5.6278054,0.3884506,5.93];
ybs[8210]=['',5.6195489,1.0441726,6.1];
ybs[8211]=['2 Peg',5.6319219,0.413927,4.57];
ybs[8212]=['',5.6261951,0.9685824,6.12];
ybs[8213]=['7 Cep',5.6204138,1.1673796,5.44];
ybs[8214]=['71 Cyg',5.6291052,0.8136339,5.24];
ybs[8215]=['ξ Gru',5.6426466,-0.7173511,5.29];
ybs[8216]=['6 PsA',5.6430734,-0.5910857,5.97];
ybs[8217]=['',5.6374405,0.2131955,6.08];
ybs[8218]=['β Aqr',5.6394989,-0.0958761,2.91];
ybs[8219]=['',5.648306,-0.9190765,6.41];
ybs[8220]=['',5.6768074,-1.3851403,6.18];
ybs[8221]=['',5.6442132,-0.4278235,6.43];
ybs[8222]=['',5.6484043,-0.7813892,5.57];
ybs[8223]=['',5.6327153,0.9256443,6.02];
ybs[8224]=['β Cep',5.6237936,1.2328624,3.23];
ybs[8225]=['',5.6036155,1.4067455,5.97];
ybs[8226]=['',5.6428564,0.4096719,6.7];
ybs[8227]=['',5.6522247,-0.747812,6.32];
ybs[8228]=['',5.6376202,0.9197492,6.16];
ybs[8229]=['',5.6350961,1.056571,5.53];
ybs[8230]=['',5.6544791,-0.5169218,6.41];
ybs[8231]=['37 Cap',5.6541286,-0.3491669,5.69];
ybs[8232]=['',5.6442317,0.8736404,5.75];
ybs[8233]=['',5.6560067,-0.4079777,6.4];
ybs[8234]=['',5.6459462,0.8016715,6.25];
ybs[8235]=['',5.6698053,-1.1300076,6.2];
ybs[8236]=['',5.6521036,0.398516,6.47];
ybs[8237]=['',5.6557659,-0.0681428,5.77];
ybs[8238]=['ρ Cyg',5.6489415,0.7970978,4.02];
ybs[8239]=['8 PsA',5.6600685,-0.4553982,5.73];
ybs[8240]=['ν Oct',5.6869494,-1.3493072,3.76];
ybs[8241]=['72 Cyg',5.6526412,0.67392,4.9];
ybs[8242]=['7 PsA',5.6629805,-0.5754159,6.11];
ybs[8243]=['',5.655261,0.4935088,6.31];
ybs[8244]=['',5.6559247,0.4281468,6.11];
ybs[8245]=['',5.6507768,0.9036764,6.15];
ybs[8246]=['ε Cap',5.6638419,-0.3383652,4.68];
ybs[8247]=['',5.659214,0.5259468,6.36];
ybs[8248]=['',5.657943,0.7933156,5.53];
ybs[8249]=['',5.6656077,-0.0054272,6.25];
ybs[8250]=['ξ Aqr',5.6665587,-0.1356958,4.69];
ybs[8251]=['3 Peg',5.6662065,0.1168968,6.18];
ybs[8252]=['74 Cyg',5.6620798,0.7067319,5.01];
ybs[8253]=['5 Peg',5.6660986,0.3385585,5.45];
ybs[8254]=['',5.6729644,-0.586421,6.28];
ybs[8255]=['',5.67747,-0.9124445,6.21];
ybs[8256]=['4 Peg',5.6697269,0.1021231,5.67];
ybs[8257]=['',5.6800962,-0.971405,6.33];
ybs[8258]=['',5.6641929,0.7814875,6.2];
ybs[8259]=['',5.6740952,-0.1832098,6.08];
ybs[8260]=['',5.6703109,0.4464289,6.16];
ybs[8261]=['',5.6645942,0.9445988,6.15];
ybs[8262]=['',5.671588,0.3550863,5.85];
ybs[8263]=['25 Aqr',5.674253,0.0405513,5.1];
ybs[8264]=['γ Cap',5.6769236,-0.2894152,3.68];
ybs[8265]=['9 Cep',5.665293,1.0849193,4.73];
ybs[8266]=['λ Oct',5.7308206,-1.4422791,5.29];
ybs[8267]=['',5.6701557,1.0047648,5.62];
ybs[8268]=['',5.684399,-0.4367091,6.49];
ybs[8269]=['42 Cap',5.6832296,-0.2437736,5.18];
ybs[8270]=['75 Cyg',5.6761243,0.756667,5.11];
ybs[8271]=['41 Cap',5.6854303,-0.4046086,5.24];
ybs[8272]=['',5.7028611,-1.2379199,6.01];
ybs[8273]=['26 Aqr',5.6856725,0.0238363,5.67];
ybs[8274]=['κ Cap',5.6881587,-0.3278744,4.73];
ybs[8275]=['7 Peg',5.6859914,0.100539,5.3];
ybs[8276]=['',5.6780028,0.9590981,6.2];
ybs[8277]=['76 Cyg',5.6822582,0.7135875,6.11];
ybs[8278]=['',5.6871763,0.1903325,6.09];
ybs[8279]=['',5.6906442,-0.3410396,6.22];
ybs[8280]=['',5.9862095,-1.5482241,6.57];
ybs[8281]=['44 Cap',5.6898876,-0.2499147,5.88];
ybs[8282]=['',5.6843536,0.6211743,6.07];
ybs[8283]=['',5.6845835,0.8001679,6.17];
ybs[8284]=['',5.6965839,-0.6714541,6.3];
ybs[8285]=['77 Cyg',5.6857888,0.7183373,5.69];
ybs[8286]=['π1 Cyg',5.6841786,0.8948323,4.67];
ybs[8287]=['45 Cap',5.6940084,-0.2560146,5.99];
ybs[8288]=['',5.70053,-0.8624962,6.45];
ybs[8289]=['',5.6866628,0.8670936,6.09];
ybs[8290]=['ι PsA',5.6984288,-0.5749937,4.34];
ybs[8291]=['',5.6889538,0.7196978,5.49];
ybs[8292]=['79 Cyg',5.6904399,0.6695888,5.65];
ybs[8293]=['ε Peg',5.6943358,0.1737636,2.39];
ybs[8294]=['μ1 Cyg',5.6937917,0.503068,4.73];
ybs[8295]=['μ2 Cyg',5.6937699,0.5030729,6.08];
ybs[8296]=['46 Cap',5.6982219,-0.1571037,5.09];
ybs[8297]=['',5.6866495,1.0358817,6.08];
ybs[8298]=['9 Peg',5.6956234,0.3042282,4.34];
ybs[8299]=['',5.6957132,0.2592327,5.94];
ybs[8300]=['κ Peg',5.696047,0.4490037,4.13];
ybs[8301]=['μ Cep',5.6899607,1.0273133,4.08];
ybs[8302]=['11 Cep',5.6817653,1.2460201,4.56];
ybs[8303]=['47 Cap',5.7037508,-0.1604728,6];
ybs[8304]=['λ Cap',5.7049346,-0.1969491,5.58];
ybs[8305]=['',5.7006075,0.627245,6.4];
ybs[8306]=['12 Peg',5.7023366,0.4019536,5.29];
ybs[8307]=['δ Cap',5.7072193,-0.2800489,2.87];
ybs[8308]=['',5.713286,-0.8241744,5.58];
ybs[8309]=['',5.686608,1.2636333,5.17];
ybs[8310]=['',5.7037114,0.4475855,6.28];
ybs[8311]=['θ PsA',5.7105435,-0.5378506,5.01];
ybs[8312]=['',5.6957496,1.0915567,5.95];
ybs[8313]=['11 Peg',5.7077507,0.0483063,5.64];
ybs[8314]=['',5.702746,0.7529738,6.54];
ybs[8315]=['',5.7068303,0.3015189,6.21];
ybs[8316]=['',5.7217956,-1.1280094,5.62];
ybs[8317]=['',5.7096428,-0.1018485,6.17];
ybs[8318]=['ο Ind',5.7257217,-1.213823,5.53];
ybs[8319]=['ν Cep',5.6982986,1.0681765,4.29];
ybs[8320]=['π2 Cyg',5.7047901,0.8620347,4.23];
ybs[8321]=['',5.7110577,0.6398793,6.47];
ybs[8322]=['',5.7186978,-0.2206244,6.31];
ybs[8323]=['',5.7125408,0.6759751,6.12];
ybs[8324]=['12 Cep',5.7069606,1.0607136,5.52];
ybs[8325]=['',5.7210941,-0.292559,6.38];
ybs[8326]=['',5.7171177,0.3585716,6.29];
ybs[8327]=['',5.704264,1.2257857,6.29];
ybs[8328]=['14 Peg',5.7186593,0.5280737,5.04];
ybs[8329]=['13 Peg',5.7202156,0.3031263,5.29];
ybs[8330]=['',5.7176202,0.7196178,6.48];
ybs[8331]=['',5.7275725,-0.3235908,6.16];
ybs[8332]=['',5.7152106,1.0708441,6.17];
ybs[8333]=['',5.7263915,0.3474826,5.77];
ybs[8334]=['',5.7238388,0.6914849,6.17];
ybs[8335]=['',5.7295672,0.3727297,6.89];
ybs[8336]=['μ Cap',5.734465,-0.2350722,5.08];
ybs[8337]=['',5.7441582,-1.0786643,5.9];
ybs[8338]=['γ Gru',5.7376826,-0.6506907,3.01];
ybs[8339]=['15 Peg',5.7302748,0.503984,5.53];
ybs[8340]=['',5.7357357,-0.1785224,6.59];
ybs[8341]=['16 Peg',5.7327988,0.4539244,5.08];
ybs[8342]=['',5.7273413,0.9752837,5.71];
ybs[8343]=['',5.7353584,0.3447272,5.68];
ybs[8344]=['',5.7370528,0.1212586,6.15];
ybs[8345]=['',5.7381419,-0.0731799,5.71];
ybs[8346]=['',5.724936,1.1490439,6.37];
ybs[8347]=['π Ind',5.7484501,-1.0090751,6.19];
ybs[8348]=['',5.739981,-0.0561614,6.2];
ybs[8349]=['',5.7382679,0.3456024,6.39];
ybs[8350]=['',5.7462424,-0.5327231,6.41];
ybs[8351]=['',5.7483655,-0.6487372,5.46];
ybs[8352]=['',5.7512394,-0.657345,6.18];
ybs[8353]=['δ Ind',5.7556427,-0.9583331,4.4];
ybs[8354]=['κ1 Ind',5.7583851,-1.0284883,6.12];
ybs[8355]=['',5.7754449,-1.3539836,6.41];
ybs[8356]=['13 Cep',5.7398501,0.9895095,5.8];
ybs[8357]=['',5.747453,0.3721637,6.4];
ybs[8358]=['17 Peg',5.7499621,0.2122355,5.54];
ybs[8359]=['',5.7415656,1.0755657,6.13];
ybs[8360]=['',5.742032,1.1415202,5.86];
ybs[8361]=['',5.7558215,-0.0932116,6.33];
ybs[8362]=['',5.7495828,0.8508903,6.42];
ybs[8363]=['',5.7582893,-0.3682395,6.12];
ybs[8364]=['',5.761112,-0.668652,5.5];
ybs[8365]=['',5.7801943,-1.3270335,5.95];
ybs[8366]=['',5.7664962,-0.9738622,6.01];
ybs[8367]=['',5.7588378,-0.0748539,6.22];
ybs[8368]=['',5.7471323,1.1119364,4.91];
ybs[8369]=['',5.7492522,1.1561049,6.43];
ybs[8370]=['18 Peg',5.7639742,0.1187172,6];
ybs[8371]=['η PsA',5.767606,-0.4951316,5.42];
ybs[8372]=['ε Ind',5.7794184,-0.9896175,4.69];
ybs[8373]=['',5.756996,1.0957619,5.93];
ybs[8374]=['',5.7594611,1.0077992,6.59];
ybs[8375]=['28 Aqr',5.7682175,0.0120375,5.58];
ybs[8376]=['',5.7648993,0.577541,6.46];
ybs[8377]=['20 Peg',5.7680596,0.2304605,5.6];
ybs[8378]=['19 Peg',5.7684099,0.1455942,5.65];
ybs[8379]=['',5.7733605,-0.3109945,6.28];
ybs[8380]=['',5.7508876,1.3104037,6.35];
ybs[8381]=['29 Aqr',5.774414,-0.2945973,6.37];
ybs[8382]=['',5.7721661,0.1930121,6.37];
ybs[8383]=['',5.7783149,-0.5204397,7.1];
ybs[8384]=['',5.764734,1.0920983,6.66];
ybs[8385]=['16 Cep',5.7573942,1.2787022,5.03];
ybs[8386]=['30 Aqr',5.7778764,-0.1123528,5.54];
ybs[8387]=['ο Aqr',5.7779949,-0.0361303,4.69];
ybs[8388]=['',5.7704364,0.9244495,5.78];
ybs[8389]=['21 Peg',5.7778047,0.2002163,5.8];
ybs[8390]=['13 PsA',5.7831616,-0.5206538,6.47];
ybs[8391]=['14 Cep',5.7712228,1.013782,5.56];
ybs[8392]=['',5.7755299,0.7807742,5.6];
ybs[8393]=['',5.7840419,-0.4666497,5.96];
ybs[8394]=['κ2 Ind',5.7903957,-1.0393505,5.62];
ybs[8395]=['32 Aqr',5.7844126,-0.0143326,5.3];
ybs[8396]=['λ Gru',5.7908472,-0.6886648,4.46];
ybs[8397]=['',5.782924,0.5764361,6.38];
ybs[8398]=['ν Peg',5.7882113,0.0897842,4.84];
ybs[8399]=['α Aqr',5.7887385,-0.004085,2.96];
ybs[8400]=['',5.785736,0.4670401,5.78];
ybs[8401]=['18 Cep',5.7788095,1.1031348,5.29];
ybs[8402]=['ξ Cep',5.7782959,1.1294549,4.29];
ybs[8403]=['ι Aqr',5.7917835,-0.2405747,4.27];
ybs[8404]=['23 Peg',5.7873901,0.5070095,5.7];
ybs[8405]=['',5.8132961,-1.3228521,6.55];
ybs[8406]=['',5.785643,0.8173423,6.13];
ybs[8407]=['',5.7881819,0.7888519,6.44];
ybs[8408]=['',5.7487496,1.4478141,6.98];
ybs[8409]=['',5.7890149,0.7871461,5.14];
ybs[8410]=['α Gru',5.8002729,-0.8181219,1.74];
ybs[8411]=['20 Cep',5.7837631,1.0973065,5.27];
ybs[8412]=['',5.7881316,0.8432965,6.27];
ybs[8413]=['19 Cep',5.7844097,1.0884834,5.11];
ybs[8414]=['',5.7897657,0.7912336,6.19];
ybs[8415]=['ι Peg',5.7937135,0.4438531,3.76];
ybs[8416]=['μ PsA',5.8005928,-0.5742554,4.5];
ybs[8417]=['',5.8186118,-1.3269593,6.15];
ybs[8418]=['υ PsA',5.8008317,-0.5926733,4.99];
ybs[8419]=['',5.7894464,0.9848681,6.39];
ybs[8420]=['',5.7958399,0.3414136,5.75];
ybs[8421]=['',5.7959648,0.3156701,6.35];
ybs[8422]=['',5.8020267,-0.5766444,6.37];
ybs[8423]=['25 Peg',5.7973834,0.3802872,5.78];
ybs[8424]=['35 Aqr',5.8029608,-0.3217236,5.81];
ybs[8425]=['',5.807858,-0.8381193,6.43];
ybs[8426]=['',5.7992773,0.4473239,6.11];
ybs[8427]=['',5.793413,1.0284658,6.32];
ybs[8428]=['',5.7948225,0.931887,6.14];
ybs[8429]=['',5.8073477,-0.5921591,5.37];
ybs[8430]=['',5.7986653,0.8706144,6.42];
ybs[8431]=['',5.8075606,-0.4922873,6.44];
ybs[8432]=['τ PsA',5.8082708,-0.5665651,4.92];
ybs[8433]=['',5.8005747,0.7998525,6.11];
ybs[8434]=['π1 Peg',5.8032432,0.5804715,5.58];
ybs[8435]=['θ Peg',5.8079205,0.1096821,3.53];
ybs[8436]=['',5.8087202,-0.0664549,6.27];
ybs[8437]=['38 Aqr',5.8100189,-0.2003353,5.46];
ybs[8438]=['',5.8096488,-0.0729653,6.01];
ybs[8439]=['π2 Peg',5.8065625,0.5805808,4.29];
ybs[8440]=['',5.8082362,0.3438911,6.18];
ybs[8441]=['',5.8085435,0.2568527,6.33];
ybs[8442]=['',5.8119716,-0.3690635,6.09];
ybs[8443]=['',5.8096924,0.2043967,5.78];
ybs[8444]=['28 Peg',5.8090301,0.3676475,6.46];
ybs[8445]=['',5.8104263,0.534764,6.32];
ybs[8446]=['',5.8150052,0.2814765,5.95];
ybs[8447]=['39 Aqr',5.8179292,-0.246212,6.03];
ybs[8448]=['',5.8112637,0.8885479,5.4];
ybs[8449]=['',5.8204179,-0.4579864,6.17];
ybs[8450]=['ζ Cep',5.8096077,1.0173133,3.35];
ybs[8451]=['',5.8160835,0.4369765,5.92];
ybs[8452]=['',5.8191076,-0.080875,6.39];
ybs[8453]=['24 Cep',5.8038163,1.2640988,4.79];
ybs[8454]=['λ Cep',5.8124187,1.0384922,5.04];
ybs[8455]=['',5.8238074,-0.437966,5.58];
ybs[8456]=['ψ Oct',5.8446223,-1.3512969,5.51];
ybs[8457]=['',5.8138632,0.993551,5.24];
ybs[8458]=['',5.8058157,1.2600861,6.37];
ybs[8459]=['',5.8077997,1.2255591,5.5];
ybs[8460]=['',5.8188079,0.6054854,5.33];
ybs[8461]=['',5.8143459,1.0327388,6.3];
ybs[8462]=['',5.8280542,-0.7207207,6.23];
ybs[8463]=['λ PsA',5.8263496,-0.4831002,5.43];
ybs[8464]=['',5.8146205,1.0619686,5.35];
ybs[8465]=['41 Aqr',5.8261853,-0.3662892,5.32];
ybs[8466]=['ε Oct',5.8550844,-1.4023939,5.1];
ybs[8467]=['',5.8226159,0.5008315,5.89];
ybs[8468]=['',5.8159529,1.1061604,5.79];
ybs[8469]=['',5.8322318,-0.7743041,6.1];
ybs[8470]=['',5.8234259,0.6946801,4.49];
ybs[8471]=['μ1 Gru',5.8322874,-0.7201067,4.79];
ybs[8472]=['',5.823038,0.7946146,5.53];
ybs[8473]=['μ2 Gru',5.8359033,-0.7250054,5.1];
ybs[8474]=['',5.8271148,0.7512122,5.71];
ybs[8475]=['',5.8223219,1.1039157,6.11];
ybs[8476]=['',5.8331865,0.1507457,6.21];
ybs[8477]=['',5.8363906,-0.4504793,6.15];
ybs[8478]=['',5.8171484,1.280971,6.08];
ybs[8479]=['ε Cep',5.8279369,0.997125,4.19];
ybs[8480]=['',5.8357767,-0.0263306,6.15];
ybs[8481]=['42 Aqr',5.8369784,-0.2224175,5.34];
ybs[8482]=['',5.8379757,-0.4023361,6.17];
ybs[8483]=['1 Lac',5.8326061,0.6603719,4.13];
ybs[8484]=['θ Aqr',5.8370494,-0.1343123,4.16];
ybs[8485]=['',5.837255,-0.1562451,5.79];
ybs[8486]=['',5.8441592,-0.9344439,5.37];
ybs[8487]=['α Tuc',5.845495,-1.0501923,2.86];
ybs[8488]=['',5.8350824,0.4868056,6.37];
ybs[8489]=['44 Aqr',5.838218,-0.0924914,5.75];
ybs[8490]=['υ Oct',5.9103129,-1.4951203,5.77];
ybs[8491]=['',5.8340887,1.000213,5.88];
ybs[8492]=['',5.8423551,-0.0026135,6.39];
ybs[8493]=['45 Aqr',5.8466315,-0.2306766,5.95];
ybs[8494]=['',5.8545131,-1.0021941,6.34];
ybs[8495]=['',5.8455656,0.6607402,6.17];
ybs[8496]=['25 Cep',5.8415425,1.0976807,5.75];
ybs[8497]=['ρ Aqr',5.8517295,-0.1349609,5.37];
ybs[8498]=['30 Peg',5.8527008,0.1025888,5.37];
ybs[8499]=['',5.8547218,0.1444297,6.17];
ybs[8500]=['ν Ind',5.873024,-1.25954,5.29];
ybs[8501]=['47 Aqr',5.8579971,-0.3754144,5.13];
ybs[8502]=['',5.8547761,0.4716548,6.47];
ybs[8503]=['γ Aqr',5.8580116,-0.0226639,3.84];
ybs[8504]=['',5.8527665,0.8913277,6.42];
ybs[8505]=['31 Peg',5.8572383,0.2145695,5.01];
ybs[8506]=['π1 Gru',5.8634046,-0.8003888,6.62];
ybs[8507]=['32 Peg',5.8561515,0.496008,4.81];
ybs[8508]=['2 Lac',5.8544967,0.8137636,4.57];
ybs[8509]=['π2 Gru',5.8651539,-0.8000531,5.62];
ybs[8510]=['',5.8404789,1.3365046,6.66];
ybs[8511]=['',5.8788024,-1.3077073,6.04];
ybs[8512]=['',5.8752525,-1.2277056,5.78];
ybs[8513]=['',5.8581951,0.7359535,6.41];
ybs[8514]=['49 Aqr',5.866422,-0.4306338,5.53];
ybs[8515]=['',5.8662776,-0.1240133,5.93];
ybs[8516]=['',5.8734013,-1.0071938,5.32];
ybs[8517]=['33 Peg',5.8664725,0.3654258,6.04];
ybs[8518]=['51 Aqr',5.8687793,-0.0828655,5.78];
ybs[8519]=['50 Aqr',5.870355,-0.2345772,5.76];
ybs[8520]=['',5.8627854,1.0013536,6.16];
ybs[8521]=['',5.8672471,0.6747908,6.22];
ybs[8522]=['',5.8625353,1.0909858,6.04];
ybs[8523]=['β Lac',5.8654051,0.913124,4.43];
ybs[8524]=['π Aqr',5.8737734,0.0256005,4.66];
ybs[8525]=['δ Tuc',5.8841848,-1.1323123,4.48];
ybs[8526]=['4 Lac',5.8696708,0.865082,4.57];
ybs[8527]=['',5.8780194,-0.4117762,6.29];
ybs[8528]=['',5.8753164,0.323476,6.26];
ybs[8529]=['53 Aqr',5.8796348,-0.2906299,6.57];
ybs[8530]=['53 Aqr',5.8796494,-0.2906493,6.35];
ybs[8531]=['',5.809549,1.5009159,5.27];
ybs[8532]=['',5.8899936,-1.1763393,5.55];
ybs[8533]=['34 Peg',5.8796141,0.0782455,5.75];
ybs[8534]=['',5.8797558,0.6550819,6.46];
ybs[8535]=['',5.8635157,1.3671565,6.76];
ybs[8536]=['35 Peg',5.885,0.083519,4.79];
ybs[8537]=['ν Gru',5.8890571,-0.6814127,5.47];
ybs[8538]=['',5.8826813,0.6963755,6.14];
ybs[8539]=['',5.8802416,0.9865108,6.57];
ybs[8540]=['',5.8842569,0.5572835,5.98];
ybs[8541]=['δ1 Gru',5.8918323,-0.7575703,3.97];
ybs[8542]=['',5.8750605,1.2367442,5.47];
ybs[8543]=['ζ1 Aqr',5.8892871,0.001215,4.59];
ybs[8544]=['ζ2 Aqr',5.8893162,0.0012199,4.42];
ybs[8545]=['δ2 Gru',5.8939654,-0.7620001,4.11];
ybs[8546]=['26 Cep',5.8802495,1.1383351,5.46];
ybs[8547]=['36 Peg',5.8905126,0.160899,5.58];
ybs[8548]=['',5.8937132,-0.4715386,5.95];
ybs[8549]=['',5.8904479,0.4686732,5.79];
ybs[8550]=['',5.8946567,-0.2238369,6.4];
ybs[8551]=['37 Peg',5.8942036,0.0789194,5.48];
ybs[8552]=['56 Aqr',5.8958325,-0.2529977,6.37];
ybs[8553]=['',5.8857443,1.120071,6.29];
ybs[8554]=['',5.8927607,0.6251,6.56];
ybs[8555]=['ζ PsA',5.8986181,-0.4534955,6.43];
ybs[8556]=['δ Cep',5.8897196,1.0211085,3.75];
ybs[8557]=['5 Lac',5.8916318,0.834214,4.36];
ybs[8558]=['σ Aqr',5.897342,-0.1847931,4.82];
ybs[8559]=['38 Peg',5.8941161,0.5700696,5.65];
ybs[8560]=['',5.8941169,0.8629989,6.4];
ybs[8561]=['β PsA',5.9013692,-0.5629695,4.29];
ybs[8562]=['',5.9211849,-1.3732368,6.15];
ybs[8563]=['ρ1 Cep',5.8766295,1.3766339,5.83];
ybs[8564]=['6 Lac',5.8959257,0.7542177,4.51];
ybs[8565]=['',5.9001331,-0.0492325,6.16];
ybs[8566]=['',5.900174,-0.1128303,6.14];
ybs[8567]=['ν Tuc',5.9086714,-1.0802129,4.81];
ybs[8568]=['58 Aqr',5.9018885,-0.1887609,6.38];
ybs[8569]=['',5.9008883,0.5171953,6.35];
ybs[8570]=['α Lac',5.8992648,0.8791709,3.77];
ybs[8571]=['39 Peg',5.9054753,0.3546595,6.42];
ybs[8572]=['',5.9063578,0.2784474,6.32];
ybs[8573]=['',5.9045198,0.695866,5.88];
ybs[8574]=['',5.9036204,0.9447108,6.35];
ybs[8575]=['60 Aqr',5.9120805,-0.0258909,5.89];
ybs[8576]=['ρ2 Cep',5.8905743,1.3773121,5.5];
ybs[8577]=['υ Aqr',5.9151103,-0.3598433,5.2];
ybs[8578]=['',5.9210276,-1.0086709,6.23];
ybs[8579]=['',5.9094947,0.9898749,5.71];
ybs[8580]=['',5.9059789,1.2218029,6.6];
ybs[8581]=['',5.9191411,-0.4171362,5.97];
ybs[8582]=['η Aqr',5.9177734,-0.0004638,4.02];
ybs[8583]=['',5.9069648,1.2298369,6.34];
ybs[8584]=['',5.9017089,1.3319793,5.68];
ybs[8585]=['σ1 Gru',5.9232257,-0.706713,6.28];
ybs[8586]=['',5.9235246,-0.5510487,5.82];
ybs[8587]=['σ2 Gru',5.9253682,-0.7068571,5.86];
ybs[8588]=['8 Lac',5.9195157,0.6933349,5.73];
ybs[8589]=['',5.9207202,0.6225286,6.1];
ybs[8590]=['',5.9231015,0.2057404,6.4];
ybs[8591]=['',5.9193726,0.8754939,6.29];
ybs[8592]=['',5.9190781,0.9801942,6.38];
ybs[8593]=['',5.9251507,0.2211054,6.3];
ybs[8594]=['',5.9236894,0.6238442,6.3];
ybs[8595]=['κ Aqr',5.9282881,-0.0722002,5.03];
ybs[8596]=['',5.9350454,-0.9180557,6.65];
ybs[8597]=['',5.931002,-0.1362473,6.23];
ybs[8598]=['9 Lac',5.9258186,0.9012269,4.63];
ybs[8599]=['',5.9328784,-0.5001475,6.47];
ybs[8600]=['31 Cep',5.9174859,1.2869012,5.08];
ybs[8601]=['',5.9334327,-0.575783,5.66];
ybs[8602]=['',5.9299833,0.7901876,6.4];
ybs[8603]=['40 Peg',5.9329231,0.3423232,5.82];
ybs[8604]=['',5.9371815,-0.492771,6.31];
ybs[8605]=['',5.9424828,-1.0006055,5.97];
ybs[8606]=['',5.9311614,0.9928696,5.21];
ybs[8607]=['10 Lac',5.9343369,0.6831529,4.88];
ybs[8608]=['',5.9399953,-0.5334985,5.87];
ybs[8609]=['41 Peg',5.9368796,0.3450986,6.21];
ybs[8610]=['',5.9235834,1.3170749,5.79];
ybs[8611]=['',5.9357153,0.6577155,6.03];
ybs[8612]=['30 Cep',5.9309548,1.1113531,5.19];
ybs[8613]=['ε PsA',5.9411905,-0.4703993,4.17];
ybs[8614]=['',5.9415587,-0.0604309,6.31];
ybs[8615]=['β Oct',5.9679518,-1.4187633,4.15];
ybs[8616]=['',5.9417163,0.2555369,5.71];
ybs[8617]=['11 Lac',5.9397182,0.774369,4.46];
ybs[8618]=['',5.93858,0.9413915,5.93];
ybs[8619]=['ζ Peg',5.9443015,0.1906461,3.4];
ybs[8620]=['',5.9500467,-0.8223739,5.98];
ybs[8621]=['β Gru',5.9502722,-0.816687,2.1];
ybs[8622]=['19 PsA',5.9486822,-0.5108382,6.17];
ybs[8623]=['',5.9443383,0.5420584,6.34];
ybs[8624]=['',5.9504376,-0.7706635,6.07];
ybs[8625]=['12 Lac',5.9439944,0.703671,5.25];
ybs[8626]=['ο Peg',5.9453793,0.5131157,4.79];
ybs[8627]=['',5.9464238,0.2549627,5.9];
ybs[8628]=['',5.9445176,0.7267775,5.94];
ybs[8629]=['ρ Gru',5.9537921,-0.7212106,4.85];
ybs[8630]=['',5.9514584,-0.1434593,6.45];
ybs[8631]=['',5.9576525,-1.0543047,6.3];
ybs[8632]=['67 Aqr',5.9522303,-0.1199163,6.41];
ybs[8633]=['',5.947499,0.9424922,6.12];
ybs[8634]=['66 Aqr',5.9538782,-0.3270424,4.69];
ybs[8635]=['η Peg',5.9508047,0.5290691,2.94];
ybs[8636]=['',5.9503618,0.6613891,6.43];
ybs[8637]=['',5.9508479,0.8248541,6.39];
ybs[8638]=['',5.9541209,0.1925327,6.51];
ybs[8639]=['',5.9554107,0.6904129,5.95];
ybs[8640]=['η Gru',5.9633321,-0.932143,4.85];
ybs[8641]=['13 Lac',5.9553931,0.7314911,5.08];
ybs[8642]=['',5.9633766,-0.8107941,5.51];
ybs[8643]=['',5.9654087,-0.8532288,6.62];
ybs[8644]=['',5.9668976,-0.8655665,6.48];
ybs[8645]=['45 Peg',5.9617099,0.3396244,6.25];
ybs[8646]=['',5.958352,0.918209,6.55];
ybs[8647]=['',5.9679579,-0.8176324,6.56];
ybs[8648]=['ξ Oct',5.9862076,-1.3968062,5.35];
ybs[8649]=['',5.982481,-1.3431633,6.73];
ybs[8650]=['ξ Peg',5.967123,0.2140703,4.19];
ybs[8651]=['',5.9644266,0.7790902,5.76];
ybs[8652]=['λ Peg',5.9662998,0.4129113,3.95];
ybs[8653]=['',5.9703221,-0.5946071,6.28];
ybs[8654]=['',5.9754572,-1.0749725,6.37];
ybs[8655]=['68 Aqr',5.97117,-0.3406999,5.26];
ybs[8656]=['',5.972413,-0.6654809,6.71];
ybs[8657]=['',5.9799983,-1.2261788,6.34];
ybs[8658]=['τ1 Aqr',5.9718195,-0.2437127,5.66];
ybs[8659]=['',5.9729166,-0.4506305,6.3];
ybs[8660]=['ε Gru',5.9760069,-0.8940299,3.49];
ybs[8661]=['70 Aqr',5.9752319,-0.1826098,6.19];
ybs[8662]=['',5.9694019,1.0223336,6.36];
ybs[8663]=['',5.9733327,0.6546625,5.9];
ybs[8664]=['τ2 Aqr',5.9800076,-0.235612,4.01];
ybs[8665]=['',5.9819276,-0.5709373,6.33];
ybs[8666]=['',5.9795568,0.1845128,6.54];
ybs[8667]=['',5.9756916,0.9513408,6.12];
ybs[8668]=['',5.9751381,1.1001007,6.06];
ybs[8669]=['μ Peg',5.9814475,0.4310028,3.48];
ybs[8670]=['',5.9865964,-0.6817925,5.42];
ybs[8671]=['',5.9901465,-1.0435006,6.46];
ybs[8672]=['',5.9760354,1.1983973,6.19];
ybs[8673]=['',5.9799074,0.9773095,5.43];
ybs[8674]=['',5.9920844,-1.1012216,6.12];
ybs[8675]=['14 Lac',5.9827959,0.733852,5.92];
ybs[8676]=['',5.9843343,0.3356947,6.4];
ybs[8677]=['',5.9817875,0.8861025,6.21];
ybs[8678]=['21 PsA',5.9878313,-0.5138766,5.97];
ybs[8679]=['ι Cep',5.9791079,1.1570394,3.52];
ybs[8680]=['γ PsA',5.993007,-0.5721584,4.46];
ybs[8681]=['',5.9867589,1.0784403,5.6];
ybs[8682]=['σ Peg',5.9920598,0.1732908,5.16];
ybs[8683]=['λ Aqr',5.993141,-0.1306627,3.74];
ybs[8684]=['15 Lac',5.9900717,0.7575728,4.94];
ybs[8685]=['τ1 Gru',5.998058,-0.8465654,6.04];
ybs[8686]=['ρ Ind',6.003289,-1.2213822,6.05];
ybs[8687]=['',5.9662443,1.452925,4.74];
ybs[8688]=['',5.9947756,0.2955621,5.64];
ybs[8689]=['74 Aqr',5.9969423,-0.2011188,5.8];
ybs[8690]=['',5.9936048,0.8814832,6.46];
ybs[8691]=['',5.9951671,0.7026798,6.34];
ybs[8692]=['',5.9941838,1.0505915,6.01];
ybs[8693]=['',5.9971932,0.7826508,5.81];
ybs[8694]=['δ Aqr',6.0020898,-0.2744929,3.27];
ybs[8695]=['78 Aqr',6.001659,-0.1241135,6.19];
ybs[8696]=['77 Aqr',6.0025663,-0.282366,5.56];
ybs[8697]=['',5.9992203,0.7063421,5.81];
ybs[8698]=['',6.0049126,-0.6334671,6.4];
ybs[8699]=['',6.0015704,0.2973205,6.12];
ybs[8700]=['1 Psc',6.0034372,0.0202163,6.11];
ybs[8701]=['',6.0043235,-0.0854193,5.72];
ybs[8702]=['ρ Peg',6.0044057,0.1554992,4.9];
ybs[8703]=['',6.0033142,0.6487482,5.91];
ybs[8704]=['',6.0075074,-0.5504657,6.1];
ybs[8705]=['δ PsA',6.0079173,-0.5662899,4.21];
ybs[8706]=['',6.0098751,-0.5492865,6.48];
ybs[8707]=['τ3 Gru',6.011831,-0.8355827,5.7];
ybs[8708]=['',6.0063742,0.6360912,5.74];
ybs[8709]=['',6.0114945,0.2084296,6.51];
ybs[8710]=['16 Lac',6.0091506,0.727761,5.59];
ybs[8711]=['',6.0091896,0.8696515,4.95];
ybs[8712]=['',6.0134983,-0.0823122,6.31];
ybs[8713]=['α PsA',6.0153103,-0.5153664,1.16];
ybs[8714]=['51 Peg',6.0140654,0.3641239,5.49];
ybs[8715]=['',6.0145641,0.0681405,6.28];
ybs[8716]=['',6.0120145,0.8513365,5.43];
ybs[8717]=['',6.0194364,-0.6183536,6.13];
ybs[8718]=['',6.014797,0.6877083,6.18];
ybs[8719]=['',6.0177184,-0.0401655,6.16];
ybs[8720]=['',6.018307,-0.0229737,6.37];
ybs[8721]=['',5.9798735,1.4868139,5.9];
ybs[8722]=['',6.0190513,0.1649501,6.43];
ybs[8723]=['',6.0196129,0.1297432,6.33];
ybs[8724]=['52 Peg',6.0217019,0.2063495,5.75];
ybs[8725]=['',6.0237874,-0.5125702,5.51];
ybs[8726]=['',6.023634,-0.2264865,6.07];
ybs[8727]=['2 Psc',6.0229209,0.0184459,5.43];
ybs[8728]=['',6.0259256,-0.437554,5.65];
ybs[8729]=['',6.0211032,0.920635,6.29];
ybs[8730]=['',6.0208265,1.0456054,6.43];
ybs[8731]=['',6.0272961,-0.4456255,6.29];
ybs[8732]=['ζ Gru',6.0297039,-0.9190888,4.12];
ybs[8733]=['',5.9960256,1.4737485,4.71];
ybs[8734]=['',6.0307354,-0.8875997,5.68];
ybs[8735]=['3 Psc',6.0280539,0.004888,6.21];
ybs[8736]=['',6.0283961,0.0542082,5.83];
ybs[8737]=['',6.0249978,0.995526,5];
ybs[8738]=['',6.0281292,0.5441463,6.6];
ybs[8739]=['',6.0313085,-0.5019446,5.55];
ybs[8740]=['',6.0273615,0.7935875,6.5];
ybs[8741]=['',6.0315142,-0.3961291,6.28];
ybs[8742]=['81 Aqr',6.0314315,-0.1215936,6.21];
ybs[8743]=['',6.0289316,0.677228,6.54];
ybs[8744]=['',6.0320028,-0.080583,5.94];
ybs[8745]=['',6.0368066,-0.6340153,6.47];
ybs[8746]=['',6.0312346,0.9983261,6.2];
ybs[8747]=['ο And',6.033298,0.7403769,3.62];
ybs[8748]=['82 Aqr',6.0364447,-0.1130928,6.15];
ybs[8749]=['',6.0374063,-0.3626114,5.97];
ybs[8750]=['',6.0361733,0.5563234,6.57];
ybs[8751]=['2 And',6.0362783,0.7479121,5.1];
ybs[8752]=['π PsA',6.0408415,-0.6048423,5.11];
ybs[8753]=['',6.0369085,0.7706211,6.39];
ybs[8754]=['',6.0475813,-1.1994879,5.52];
ybs[8755]=['',6.0366083,0.9657052,6.5];
ybs[8756]=['',6.0430764,-0.7222827,5.79];
ybs[8757]=['',6.0425903,-0.0820427,6.68];
ybs[8758]=['β Psc',6.042187,0.0683221,4.53];
ybs[8759]=['κ Gru',6.0461854,-0.9402149,5.37];
ybs[8760]=['β Peg',6.0415518,0.4917872,2.42];
ybs[8761]=['',6.0427777,0.1171334,6.41];
ybs[8762]=['',6.039383,1.0566186,6.74];
ybs[8763]=['',6.0392828,1.0237966,6.43];
ybs[8764]=['',6.0397986,1.1746709,5.24];
ybs[8765]=['3 And',6.0430657,0.8752271,4.65];
ybs[8766]=['α Peg',6.0459663,0.2670343,2.49];
ybs[8767]=['83 Aqr',6.0478824,-0.132626,5.43];
ybs[8768]=['',6.0481692,-0.2964348,6.14];
ybs[8769]=['',6.0474548,0.2907326,6.44];
ybs[8770]=['',6.0483865,0.0244635,6.39];
ybs[8771]=['',6.0639384,-1.3855436,6.12];
ybs[8772]=['θ Gru',6.0556861,-0.7579212,4.28];
ybs[8773]=['',6.0526714,0.324846,6.13];
ybs[8774]=['86 Aqr',6.0546243,-0.4127391,4.47];
ybs[8775]=['υ Gru',6.055689,-0.6771415,5.61];
ybs[8776]=['',6.0569852,-0.8641434,6.33];
ybs[8777]=['',6.0536576,0.3491647,6.3];
ybs[8778]=['',6.0573803,-0.882988,5.83];
ybs[8779]=['',6.0640677,-1.282666,6.15];
ybs[8780]=['55 Peg',6.0558001,0.1658817,4.52];
ybs[8781]=['56 Peg',6.056153,0.4461624,4.76];
ybs[8782]=['1 Cas',6.0535019,1.038725,4.85];
ybs[8783]=['',6.0576093,0.5745756,6.02];
ybs[8784]=['',6.0577813,0.3705175,5.99];
ybs[8785]=['',6.0567591,0.8056956,6.66];
ybs[8786]=['',6.0560679,0.923476,6.11];
ybs[8787]=['',6.06194,-0.5014039,5.6];
ybs[8788]=['',6.0559368,1.0440976,6.4];
ybs[8789]=['4 And',6.058299,0.8112668,5.33];
ybs[8790]=['5 And',6.0587002,0.8620317,5.7];
ybs[8791]=['',6.0607275,0.7794057,6.56];
ybs[8792]=['5 Psc',6.0631661,0.0387954,5.4];
ybs[8793]=['',6.0585405,1.1122683,6.26];
ybs[8794]=['',6.0706301,-1.1652222,6.47];
ybs[8795]=['',6.0806902,-1.4105297,6.41];
ybs[8796]=['',6.0592133,1.1225515,6.21];
ybs[8797]=['88 Aqr',6.0666608,-0.3678698,3.66];
ybs[8798]=['',6.0680079,-0.4885782,5.87];
ybs[8799]=['',6.0690743,-0.746397,5.81];
ybs[8800]=['57 Peg',6.0668036,0.1531062,5.12];
ybs[8801]=['',6.0682707,-0.2515964,6.42];
ybs[8802]=['89 Aqr',6.0687051,-0.3902966,4.69];
ybs[8803]=['',6.0699561,-0.706797,5.83];
ybs[8804]=['π Cep',6.0582952,1.3174173,4.41];
ybs[8805]=['ι Gru',6.0708678,-0.7880418,3.9];
ybs[8806]=['58 Peg',6.0689783,0.1730862,5.39];
ybs[8807]=['2 Cas',6.0671559,1.0372175,5.7];
ybs[8808]=['',6.0725244,-0.5136464,6.51];
ybs[8809]=['',6.0719162,0.308743,5.71];
ybs[8810]=['6 And',6.0705663,0.7616506,5.94];
ybs[8811]=['59 Peg',6.076456,0.1538563,5.16];
ybs[8812]=['60 Peg',6.0766955,0.4702361,6.17];
ybs[8813]=['',6.0835136,-0.8643471,6.8];
ybs[8814]=['',6.0875125,-1.0926542,6.12];
ybs[8815]=['7 And',6.0796539,0.8639689,4.52];
ybs[8816]=['',6.0821197,0.5155195,6.35];
ybs[8817]=['',6.0827288,0.9994414,5.56];
ybs[8818]=['',6.0838831,0.1947868,5.82];
ybs[8819]=['φ Aqr',6.0878319,-0.1039056,4.22];
ybs[8820]=['',6.0909336,-0.7157589,5.77];
ybs[8821]=['',6.0893716,-0.1848836,6.12];
ybs[8822]=['',6.0870179,0.8851141,6.31];
ybs[8823]=['',6.087773,0.5212811,6.41];
ybs[8824]=['',6.0888917,0.4223455,6.36];
ybs[8825]=['',6.0932632,-0.0593543,5.55];
ybs[8826]=['ψ1 Aqr',6.0946911,-0.1569419,4.21];
ybs[8827]=['61 Peg',6.0939424,0.4946862,6.49];
ybs[8828]=['',6.0998975,-1.0804522,5.66];
ybs[8829]=['',6.0879183,1.297245,5.84];
ybs[8830]=['',6.0948112,0.4340072,6.6];
ybs[8831]=['',6.0983097,-0.7748116,5.92];
ybs[8832]=['',6.0990093,-0.7173077,6.47];
ybs[8833]=['γ Tuc',6.1018524,-1.0147351,3.99];
ybs[8834]=['',6.1103824,-1.3853873,6.33];
ybs[8835]=['χ Aqr',6.0988566,-0.1331848,5.06];
ybs[8836]=['',6.0925587,1.2388992,5.56];
ybs[8837]=['γ Psc',6.1001777,0.058957,3.69];
ybs[8838]=['',6.0977648,0.9304234,5.54];
ybs[8839]=['',6.0964623,1.0831297,6.53];
ybs[8840]=['',6.1060238,-1.1759199,6.13];
ybs[8841]=['',6.1024486,-0.2027593,6.34];
ybs[8842]=['',6.1003686,0.7899349,6.43];
ybs[8843]=['ψ2 Aqr',6.1034668,-0.1585925,4.39];
ybs[8844]=['φ Gru',6.1048224,-0.7108482,5.53];
ybs[8845]=['8 And',6.1023706,0.8571502,4.85];
ybs[8846]=['',6.1032425,0.7956033,6.48];
ybs[8847]=['τ Oct',6.1528426,-1.5236955,5.49];
ybs[8848]=['γ Scl',6.107623,-0.5661158,4.41];
ybs[8849]=['9 And',6.1052689,0.7307601,6.02];
ybs[8850]=['ψ3 Aqr',6.1080858,-0.1660669,4.98];
ybs[8851]=['94 Aqr',6.108761,-0.2332279,5.08];
ybs[8852]=['',6.0996548,1.3158899,6.38];
ybs[8853]=['96 Aqr',6.1099742,-0.087764,5.55];
ybs[8854]=['',6.1100516,-0.3137987,5.93];
ybs[8855]=['',6.1080759,0.789467,6.5];
ybs[8856]=['',6.1115399,-0.5866417,6.37];
ybs[8857]=['ο Cep',6.1058111,1.1904461,4.75];
ybs[8858]=['',6.109989,0.6089327,6.32];
ybs[8859]=['11 And',6.1100331,0.8503457,5.44];
ybs[8860]=['',6.1108957,0.8460795,6.32];
ybs[8861]=['10 And',6.1117481,0.7360756,5.79];
ybs[8862]=['',6.1165455,-0.8763407,6.05];
ybs[8863]=['7 Psc',6.1140365,0.0955986,5.05];
ybs[8864]=['',6.1155697,-0.101439,6.17];
ybs[8865]=['τ Peg',6.1152172,0.416022,4.6];
ybs[8866]=['',6.1130586,1.0832559,6.45];
ybs[8867]=['63 Peg',6.1160046,0.5325181,5.59];
ybs[8868]=['',6.1182001,-0.4693294,5.64];
ybs[8869]=['',6.1154871,0.7716523,6.13];
ybs[8870]=['12 And',6.1162169,0.6680818,5.77];
ybs[8871]=['',6.114534,1.0874985,6.39];
ybs[8872]=['64 Peg',6.1207487,0.5569104,5.32];
ybs[8873]=['',6.1210222,0.4660903,6.62];
ybs[8874]=['',6.1258983,-1.0464931,6.09];
ybs[8875]=['97 Aqr',6.1242162,-0.2608045,5.2];
ybs[8876]=['65 Peg',6.1241339,0.3652063,6.29];
ybs[8877]=['98 Aqr',6.1256236,-0.3491421,3.97];
ybs[8878]=['66 Peg',6.1259313,0.2165968,5.08];
ybs[8879]=['',6.1231783,1.0512078,5.56];
ybs[8880]=['',6.1299637,-0.9374526,6.15];
ybs[8881]=['',6.1292045,-0.7509837,6.1];
ybs[8882]=['',6.1279741,0.0067652,6.31];
ybs[8883]=['',6.1313038,-0.9039953,5.75];
ybs[8884]=['',6.1289412,0.5694596,6.69];
ybs[8885]=['',6.1306747,-0.3244782,6.19];
ybs[8886]=['',6.1361795,-0.9905235,5.59];
ybs[8887]=['',6.1323383,0.719234,6.72];
ybs[8888]=['67 Peg',6.1335504,0.5669059,5.57];
ybs[8889]=['4 Cas',6.133187,1.0887204,4.98];
ybs[8890]=['υ Peg',6.135931,0.4101614,4.4];
ybs[8891]=['99 Aqr',6.1390425,-0.3585876,4.39];
ybs[8892]=['ο Gru',6.1417265,-0.9184838,5.52];
ybs[8893]=['',6.1441841,-1.1603761,6.45];
ybs[8894]=['',6.1445882,-1.0189171,5.63];
ybs[8895]=['',6.1440573,-0.8737252,6.2];
ybs[8896]=['κ Psc',6.1428128,0.0235968,4.94];
ybs[8897]=['9 Psc',6.1441805,0.0212749,6.25];
ybs[8898]=['13 And',6.1434238,0.7506381,5.75];
ybs[8899]=['',6.1476897,-0.6186832,6.32];
ybs[8900]=['69 Peg',6.145936,0.4409349,5.98];
ybs[8901]=['θ Psc',6.147308,0.1130169,4.28];
ybs[8902]=['',6.1478975,-0.198151,6.37];
ybs[8903]=['',6.1436441,1.2296923,5.6];
ybs[8904]=['',6.1523523,-1.0998064,5.68];
ybs[8905]=['',6.1521154,-0.7749473,6.43];
ybs[8906]=['',6.1519161,-0.1600387,6.18];
ybs[8907]=['',6.1521469,0.403945,6.35];
ybs[8908]=['70 Peg',6.152461,0.2243992,4.55];
ybs[8909]=['',6.1541876,-0.077426,6.25];
ybs[8910]=['',6.156485,0.85922,6.17];
ybs[8911]=['',6.1559863,1.0235572,4.91];
ybs[8912]=['',6.1589148,0.6764651,6.05];
ybs[8913]=['',6.1606655,-0.1080648,6.39];
ybs[8914]=['',6.1627322,-0.780981,6.02];
ybs[8915]=['14 And',6.1616566,0.6864917,5.22];
ybs[8916]=['',6.1628679,-0.0696478,6.49];
ybs[8917]=['100 Aqr',6.1637044,-0.3712793,6.29];
ybs[8918]=['',6.1635851,0.4974244,6.41];
ybs[8919]=['13 Psc',6.1647549,-0.0172633,6.38];
ybs[8920]=['',6.1715681,-1.3489385,5.81];
ybs[8921]=['',6.1665739,0.6117247,6.65];
ybs[8922]=['β Scl',6.1693179,-0.6583654,4.37];
ybs[8923]=['',6.1378609,1.5238019,5.58];
ybs[8924]=['101 Aqr',6.1705704,-0.3633366,4.71];
ybs[8925]=['71 Peg',6.1712524,0.3943691,5.32];
ybs[8926]=['',6.1722021,0.788101,6.24];
ybs[8927]=['',6.1732528,0.365431,6.06];
ybs[8928]=['72 Peg',6.1733332,0.548419,4.98];
ybs[8929]=['14 Psc',6.1743108,-0.020083,5.87];
ybs[8930]=['',6.1793189,-1.1273529,7.4];
ybs[8931]=['',6.1772946,-0.2643994,5.96];
ybs[8932]=['15 And',6.1762229,0.7039479,5.59];
ybs[8933]=['73 Peg',6.176308,0.5863272,5.63];
ybs[8934]=['ι Phe',6.1785167,-0.7420813,4.71];
ybs[8935]=['',6.1769047,0.6653326,6.18];
ybs[8936]=['',6.1803733,-0.128588,6.39];
ybs[8937]=['',6.1773918,1.2520833,5.84];
ybs[8938]=['',6.1820007,0.4303637,6.45];
ybs[8939]=['16 Psc',6.1840671,0.0383825,5.68];
ybs[8940]=['',6.1844864,0.5759778,6.35];
ybs[8941]=['',6.1872383,-0.5545587,6.52];
ybs[8942]=['',6.1935251,-1.3399413,6];
ybs[8943]=['',6.1896597,-0.2262521,5.65];
ybs[8944]=['',6.1906144,-0.7923011,4.74];
ybs[8945]=['74 Peg',6.1895875,0.295354,6.26];
ybs[8946]=['λ And',6.1890359,0.8125386,3.82];
ybs[8947]=['',6.1889092,0.7771278,5.8];
ybs[8948]=['75 Peg',6.1908196,0.3228432,5.53];
ybs[8949]=['',6.1908368,0.8080302,6.58];
ybs[8950]=['ι And',6.1915523,0.756863,4.29];
ybs[8951]=['θ Phe',6.1976624,-0.8122889,6.09];
ybs[8952]=['18 And',6.1958888,0.8825905,5.3];
ybs[8953]=['ω1 Aqr',6.1989312,-0.2465208,5];
ybs[8954]=['ι Psc',6.1995998,0.0998932,4.13];
ybs[8955]=['',6.1994511,0.1705936,5.97];
ybs[8956]=['',6.1956462,1.3158006,5.95];
ybs[8957]=['',6.1964782,1.2932859,5.98];
ybs[8958]=['',6.1999242,0.6588544,6.53];
ybs[8959]=['γ Cep',6.1962934,1.3566365,3.21];
ybs[8960]=['μ Scl',6.2026979,-0.5580857,5.31];
ybs[8961]=['κ And',6.2014751,0.7754668,4.14];
ybs[8962]=['',6.2026803,0.6425941,6.23];
ybs[8963]=['',6.2047658,-0.4199815,6.6];
ybs[8964]=['',6.204871,-0.2021692,5.89];
ybs[8965]=['103 Aqr',6.2067478,-0.3129392,5.34];
ybs[8966]=['',6.2059943,0.8658464,6.26];
ybs[8967]=['104 Aqr',6.2075685,-0.3092593,4.82];
ybs[8968]=['',6.2083023,0.1282415,5.89];
ybs[8969]=['λ Psc',6.2087585,0.0327624,4.5];
ybs[8970]=['',6.207955,1.0010709,6.24];
ybs[8971]=['',6.2095067,0.7869532,6.57];
ybs[8972]=['',6.2106228,-0.2679189,5.28];
ybs[8973]=['ω2 Aqr',6.2117402,-0.2521622,4.49];
ybs[8974]=['',6.209803,1.1277045,6.56];
ybs[8975]=['',6.2106147,1.0782051,6.4];
ybs[8976]=['77 Peg',6.2145293,0.182013,5.06];
ybs[8977]=['',6.2165539,-0.2650675,6.36];
ybs[8978]=['',6.2175026,-0.7851561,6.09];
ybs[8979]=['',6.2194478,-1.2285906,6.07];
ybs[8980]=['',6.2207971,-1.3734723,5.75];
ybs[8981]=['',6.218395,-1.122373,5.72];
ybs[8982]=['78 Peg',6.2171842,0.5141543,4.93];
ybs[8983]=['106 Aqr',6.2182031,-0.3172963,5.24];
ybs[8984]=['',6.2194423,-0.4563891,6.17];
ybs[8985]=['',6.2206367,0.9755858,6.51];
ybs[8986]=['',6.226184,-0.6996194,6.31];
ybs[8987]=['107 Aqr',6.2261115,-0.3242961,5.29];
ybs[8988]=['ψ And',6.2260594,0.8118842,4.95];
ybs[8989]=['19 Psc',6.2277153,0.0625514,5.04];
ybs[8990]=['',6.2284565,1.1672674,5.95];
ybs[8991]=['σ Phe',6.2316484,-0.8749227,5.18];
ybs[8992]=['',6.2322967,-1.1920053,6.89];
ybs[8993]=['τ Cas',6.2304775,1.0253674,4.87];
ybs[8994]=['',6.231551,-0.2061853,5.73];
ybs[8995]=['',6.2303679,1.0044138,5.51];
ybs[8996]=['',6.2326949,0.8190794,6.07];
ybs[8997]=['20 Psc',6.2344885,-0.0465019,5.49];
ybs[8998]=['',6.2341483,1.1851526,5.04];
ybs[8999]=['',6.2371115,-0.1096632,6.07];
ybs[9000]=['',6.238321,0.0403431,6.46];
ybs[9001]=['δ Scl',6.2388243,-0.4892674,4.57];
ybs[9002]=['',6.237402,1.1340051,6.41];
ybs[9003]=['6 Cas',6.2382382,1.0875455,5.43];
ybs[9004]=['',6.2385214,1.0485277,6.34];
ybs[9005]=['',6.2398469,1.0307982,6.33];
ybs[9006]=['',6.2414214,-0.2751298,6.24];
ybs[9007]=['21 Psc',6.2411006,0.0204805,5.77];
ybs[9008]=['',6.2424942,-1.0950563,6.59];
ybs[9009]=['',6.2420298,0.6374399,5.9];
ybs[9010]=['79 Peg',6.241927,0.5050955,5.97];
ybs[9011]=['',6.2427434,-0.4404172,6.42];
ybs[9012]=['',6.2445467,-0.172383,5.94];
ybs[9013]=['',6.2450011,0.9026671,6.44];
ybs[9014]=['',6.2459043,-0.2496622,5.72];
ybs[9015]=['80 Peg',6.2493606,0.1642477,5.79];
ybs[9016]=['108 Aqr',6.2493987,-0.328323,5.18];
ybs[9017]=['γ1 Oct',6.2530468,-1.4298001,5.11];
ybs[9018]=['22 Psc',6.2520366,0.0528425,5.55];
ybs[9019]=['',6.2517595,1.3560653,6.55];
ybs[9020]=['',6.2538727,0.379927,6.11];
ybs[9021]=['φ Peg',6.2543048,0.3354114,5.08];
ybs[9022]=['',6.2543871,-0.2470292,5.87];
ybs[9023]=['',6.2538108,1.3202037,6.39];
ybs[9024]=['82 Peg',6.2548802,0.1927696,5.3];
ybs[9025]=['',6.2558725,-0.1553218,5.75];
ybs[9026]=['24 Psc',6.2562383,-0.0533751,5.93];
ybs[9027]=['25 Psc',6.2569028,0.0381869,6.28];
ybs[9028]=['',6.2580889,-0.4211789,6.24];
ybs[9029]=['',6.2624947,-0.4702758,6.35];
ybs[9030]=['ρ Cas',6.2625311,1.0052546,4.54];
ybs[9031]=['',6.2637567,-0.7016676,6.03];
ybs[9032]=['',6.2643076,0.0036054,5.61];
ybs[9033]=['26 Psc',6.2658453,0.1251143,6.21];
ybs[9034]=['',6.2665091,-0.5554381,6.1];
ybs[9035]=['',6.2665091,-0.5547836,6.83];
ybs[9036]=['',6.2669397,0.4547004,6.54];
ybs[9037]=['',6.2676865,1.0037325,6];
ybs[9038]=['',6.2676935,0.8282154,6];
ybs[9039]=['',6.2718303,-0.4300457,6.31];
ybs[9040]=['',6.2726546,0.3969834,6.15];
ybs[9041]=['',6.2714651,1.4536591,6.59];
ybs[9042]=['',6.2742535,0.7462287,5.97];
ybs[9043]=['',6.2746214,-0.4629694,6.26];
ybs[9044]=['',6.2746007,0.9739505,5.55];
ybs[9045]=['',6.2754882,-1.0970959,5.97];
ybs[9046]=['γ2 Oct',6.2764937,-1.4324366,5.73];
ybs[9047]=['η Tuc',6.2765989,-1.1205172,5];
ybs[9048]=['',6.2764193,1.04931,6.47];
ybs[9049]=['ψ Peg',6.2773114,0.4405004,4.66];
ybs[9050]=['1 Cet',6.2799164,-0.2748907,6.26];
ybs[9051]=['',6.2801651,0.8986009,4.8];
ybs[9052]=['27 Psc',6.2813105,-0.0603654,4.86];
ybs[9053]=['',6.2819467,0.5668671,6.52];
ybs[9054]=['π Phe',6.2824362,-0.9188881,5.13];
ybs[9055]=['',6.2817476,0.811761,6.54];
ybs[9056]=['σ Cas',6.2827659,0.9748087,4.88];
ybs[9057]=['ω Psc',0.0009099,0.1214882,4.01];
ybs[9058]=['',0.0015795,-0.5129099,5.62];
ybs[9059]=['',0.0016734,0.590303,6.58];
ybs[9060]=['',0.0016734,0.590303,6.58];
ybs[9061]=['ε Tuc',0.003544,-1.142838,4.5];
ybs[9062]=['',0.0053043,-0.7713156,6.29];
ybs[9063]=['',0.0056548,0.4715139,6.46];
ybs[9064]=['',0.0061728,1.0412137,6.19];
ybs[9065]=['',0.0071003,0.7915201,6.38];
ybs[9066]=['τ Phe',0.0085912,-0.8501948,5.71];
ybs[9067]=['',0.0097154,-0.8768499,5.53];
ybs[9068]=['',0.0096961,0.874045,6.22];
ybs[9069]=['θ Oct',0.0108075,-1.3433521,4.78];
ybs[9070]=['',0.0109883,1.0702443,5.55];
ybs[9071]=['',0.0114768,0.7411479,6.25];
ybs[9072]=['29 Psc',0.0118684,-0.0511394,5.1];
ybs[9073]=['85 Peg',0.0133919,0.4743694,5.75];
ybs[9074]=['30 Psc',0.0124637,-0.1032667,4.41];
ybs[9075]=['',0.013166,-0.2544462,7.1];
ybs[9076]=['ζ Scl',0.0140756,-0.5170164,5.01];
ybs[9077]=['31 Psc',0.0144034,0.1580285,6.32];
ybs[9078]=['32 Psc',0.0148033,0.1498012,5.63];
ybs[9079]=['',0.0153167,1.1553435,5.86];
ybs[9080]=['',0.0168197,-0.3481703,6.25];
ybs[9081]=['',0.0175515,-0.4197144,6.44];
ybs[9082]=['',0.0189307,1.1124617,6.24];
ybs[9083]=['2 Cet',0.0202226,-0.300872,4.55];
ybs[9084]=['',0.0208524,1.166048,6.29];
ybs[9085]=['9 Cas',0.0224221,1.0888269,5.88];
ybs[9086]=['',0.0227889,-0.2867834,5.78];
ybs[9087]=['',0.0228233,-0.5091335,6.4];
ybs[9088]=['3 Cet',0.0235489,-0.1817244,4.94];
ybs[9089]=['',0.0245122,1.1739795,5.67];
ybs[9090]=['',0.0240697,0.7363479,6.01];
ybs[9091]=['',0.0234805,-1.2706062,7.31];
ybs[9092]=['',0.0253067,0.6066263,6.12];
ybs[9093]=['',0.0242573,-1.2451099,5.59];
ybs[9094]=['',0.0254593,0.4668109,6.25];
ybs[9095]=['',0.0262566,1.071834,5.8];
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