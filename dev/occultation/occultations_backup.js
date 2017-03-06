/*
 Return occultation predictions for a single where and when.
 
 Returns an array of N objects, of form: 
   {UT:'2016-10-20 21:08:18', star:'ZC 123', mag:5.5, ph:'DD', el:92, pa:145}
 
 Uses output from the occult.exe software from the IOTA, for N stations in North America. 
 The target 'where' needs to be within M kilometers of one of these N stations.
 Occult.exe's predictions have decimal minutes; they are replaced with seconds, for ease of interaction with the ephem.js API.
 The output is generated of 13 months, to avoid running out at the end of the year. 
 The output is massaged offline into a javascript array, and embedded directly into this file.
*/
var fetch_occultation_predictions = function(where, when, num_days_ahead /*of when*/){
  
  // the database of occultations for all stations, for an extended period of time
  var all_occns = [
    { 
      name : 'Montréal', 
      λ:-73.6, 
      φ:45.5, 
      events: [
       {UT:'2016-10-20 21:08:18', star:'ZC 123', mag:5.5, ph:'DD', el:92, pa:145, a:+1.2, b:-1.5}
      ]
    }
    //more stations here
  ];
  
  /* Returns null if nearest station is over 400 kms away. */
  var the_nearest_station = function(){
    var i, station, this_dist, min_dist = 401, result = null;
    for(i=0; i < all_occns.length; ++i){
      station = all_occns[i];
      this_dist = EPH.distance_kms(station, where);
      if (this_dist < min_dist){
        min_dist = this_dist;
        result = station;
      }
    }
    return result;
  };

  /* Use alpha-order to filter-in occultation events. */
  var filter_date_time_window = function(events){
    var result = [];
    var base_when = when.toString().substring(2).trim(); //trim leading 'UT '
    var when_start = EPH.date_time_odometer(base_when, 'day', -1);
    var when_end = EPH.date_time_odometer(base_when, 'day', num_days_ahead + 1);
    for (var i=0; i < events.length; ++i){
      if (when_start <= events[i].UT && events[i].UT <= when_end){
        result.push(events[i]);
      } 
    }
    return result;
  };
  
  var calc_local_correction_to_time_Δsecs = function(event, station){
    //the a, b params are in minutes/degree, as decimals
    var Δt_λ = event.a * (where.λ - station.λ);
    var Δt_φ = event.b * (where.φ - station.φ);
    var Δmin = Δt_λ + Δt_φ; //fractional minutes
    var result = EPH.round(Δmin * 60, 0); //now whole seconds, to be added to the UT of the event
    return result; 
  }; 
  
  var find_occns = function(where, when, num_days_ahead){
    var result = []; // {UT:'2016-10-20 21:08:18', star:'ZC 123', mag:5.5, ph:'DD', el:92, pa:145}
    var filtered_event, local_event, Δsec;  
    var nearest_station = the_nearest_station();
    var filtered_events = filter_date_time_window(nearest_station.events); 
    for (var i=0; i < filtered_events.length; ++i){
      filtered_event = filtered_events[i];
      Δsec = calc_local_correction_to_time_Δsecs(filtered_event, nearest_station);
      //conservative: make a new object, instead of changing state of the 'database' object
      local_event = {
        UT: EPH.date_time_odometer(filtered_event.UT, 'sec', Δsec),
        star: filtered_event.star,
        mag: filtered_event.mag,
        ph: filtered_event.ph,
        el: filtered_event.el,
        pa: filtered_event.pa
      };
      result.push(local_event);
    }
    return result; 
  };
  
  return find_occns(where, when, num_days_ahead);
};
