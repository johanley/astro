/* 
 Create collections of lines in a given constellation.
 A 'polyline' is a collection of points, to be drawn on a canvas.
*/
var show = function(input){

  var NL = '\r\n'; // new line

  //for Bayer designations only; Flamsteed designations are taken as is.
  var letter_map = {};
  letter_map.alp = 'α';
  letter_map.bet = 'β';
  letter_map.gam = 'γ';
  letter_map.del = 'δ';
  letter_map.eps = 'ε';
  letter_map.zet = 'ζ';
  letter_map.eta = 'η';
  letter_map.the = 'θ';
  letter_map.iot = 'ι';
  letter_map.kap = 'κ';
  letter_map.lam = 'λ';
  letter_map.mu = 'μ';
  letter_map.nu = 'ν';
  letter_map.xi = 'ξ';
  letter_map.omi = 'ο';
  letter_map.pi = 'π';
  letter_map.rho = 'ρ';
  letter_map.sig = 'σ';
  letter_map.tau = 'τ';
  letter_map.ups = 'υ';
  letter_map.phi = 'φ';
  letter_map.chi = 'χ';
  letter_map.psi = 'ψ';
  letter_map.ome = 'ω';
  
  //index into the star-messier data structure (array of items)
  var NAME = 0;
  var MAG = 3;
  var constellation = input.constellation_abbr.trim().toLowerCase(); //IAU standard abbreviations

  var is_flamsteed = function(star_name){
    return /^\d+$/.test(star_name);
  };
  
  var is_bayer_with_number = function(star_name){
    return /^\D+\d$/.test(star_name);
  };
  
  var is_bayer = function(star_name){
    return /^\D+$/.test(star_name);
  };
  
    
  var show_stars_in_the_constellation = function(){
    var i, star, star_name;
    var designation_lines = extract_designation_lines(); 
    var constell_stars = extract_constell_stars();
    var polylines = []; //star-array indices
    for (i=0; i < designation_lines.length; ++i){
       polylines.push(to_polyline(designation_lines[i], constell_stars)); 
    }
    output(polylines);
  };
  
  /* Array of arrays, containing text (Greek letters and Flamsteed numbers). */
  var extract_designation_lines = function(){
    var result = [], i;
    var lines = input.star_designations.split(NL); //N lines in a text area input element; each line is a separate polyline.
    for (i=0; i < lines.length; ++i){
      result.push(lines[i].split(' '));
    }
    return result;
  };
  
  /* Array of indices into the star catalog, for all bright stars in the given constellation. */
  var extract_constell_stars = function(){
    var i, star, star_name;
    var result = [];
    for (i=0; i < EPH.stars.length; ++i){
      star = EPH.stars[i];
      star_name = star[NAME].trim().toLowerCase(); 
      if (star_name.length > 0 && star_name.endsWith(constellation)){
        result.push(i);
      }
    }
    return result;
  }
  
  /* Find each star-name in the given order. Return array of indexes into the star catalog. */
  var to_polyline = function(desig_line, constell_stars){
    var i, j, desig, star_name, star_idx, source, target, found, letters, number;
    var result = [];
    for (i=0; i < desig_line.length; ++i){
      desig = desig_line[i];
      found = false;
      for (j=0; j < constell_stars.length; ++j){
        star_idx = constell_stars[j];
        star_name = EPH.stars[star_idx][NAME].trim().toLowerCase();
        source = star_name.split(" ")[0].trim().toLowerCase(); //the name without the constellation part
        target = desig; //Bayer or Flamsteed
        if (is_flamsteed(target)){
          //do nothing
        }
        else if (is_bayer(target) && letter_map[target]) {
          //a Greek letter was found in the map from the input 'alp' to the real letter 'α'
          target = letter_map[desig];
        }
        else if (is_bayer_with_number(target)){
          //eg 'pi6'
          letters = target.substring(0, target.length-1);  
          number = target.substring(target.length-1);
          if (letter_map[letters]){
            target = letter_map[letters] + number; 
          }
        }
        else {
          if (target.length > 0){
            error('Wacky star name: ' + target);
          }
        }
        /*
        if (letter_map[desig]){
          //a Greek letter was found in the map from the input 'alp' to the real letter 'α'
          target = letter_map[desig];
        }
        */
        if (source === target){
          result.push(star_idx);
          found = true;
          if (EPH.stars[star_idx][MAG] > 5.0){
            error(EPH.stars[star_idx][NAME] + ' is too dim. Mag: ' + EPH.stars[star_idx][MAG]);
          }
          break;
        } 
      }
      if (!found){
        error('Not found: "' + desig + '"');
      }
    }
    if (result.length !== desig_line.length){
      error('Only found ' + result.length + ' matches for ' + desig_line.length + ' designations.');
    }
    return result;
  };
  
  /* Output the result as plain text, in a syntax compatible with javascript. */
  var output = function(polylines){
    var i, j, polyline;
    var output = document.getElementById('output');
    var text = 'constell.' + input.constellation_abbr.trim() + ' = ['
    for (i=0; i < polylines.length; ++i){
      text = text + '[';
      polyline = polylines[i];
      for (j=0; j < polyline.length; ++j){
        text = text + polyline[j];
        if (j < polyline.length - 1){
          text = text + ',';
        } 
      }
      text = text + ']';
      if (i<polylines.length-1){
        text = text + ',';
      }
    }
    text = text + '];';
    output.innerHTML = text; 
  };
  
  var error = function(msg){
    var output = document.getElementById('error');
    output.innerHTML = output.innerHTML + '<br>' + msg; 
  };
  show_stars_in_the_constellation();
  
};