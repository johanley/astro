package astro;

import java.util.List;
import java.util.logging.Logger;

import astro.util.Util;

/** 
 * List of 31 radar station id's and names (from https://weather.gc.ca/radar/index_e.html): 
  WHK - AB Carvel, near Edmonton
  WHN - AB Jimmy Lake, near Cold Lake
  XBU - AB Schuler, near Medicine Hat
  WWW - AB Spirit River, near Grande Prairie
  XSM - AB Strathmore, near Calgary
  WUJ - BC Aldergrove, near Vancouver
  XPG - BC Prince George
  XSS - BC Silver Star Mountain,near Vernon
  XSI - BC Victoria
  XFW - MB Foxwarren, near Brandon
  XWL - MB Woodlands, near Winnipeg
  XNC - NB Chipman, near Frederiction
  WTP - NL Holyrood, near St. John's
  XME - NL Marble Mountain, near Corner Brook
  XGO - NS Halifax 
  XMB - NS Marion Bridge, near Sydney
  WBI - ON Britt, near Sudbury
  XDR - ON Dryden
  WSO - ON Exeter, near London
  XFT - ON Franktown, near Ottawa
  WKR - ON King City, near Toronto
  WGJ - ON Montreal River, near Sault Ste. Marie
  XTI - ON Northeast Ontario, near Timmins
  XNI - ON Superior West, near Thunder Bay
  WMB - QC Lac Castor, near Saguenay
  XLA - QC Landrienne, near Rouyn-Noranda
  WMN - QC McGill, near Montréal
  XAM - QC Val d'Irène, near Mont-Joli
  WVY - QC Villeroy, Trois-Rivières
  XBE - SK Bethune, near Regina
  XRA - SK Radisson, near Saskatoon
 */
final class Radar {
  
  Radar(Integer numImages){
    this.numImages = numImages;
  }
  
  String returnJsonMostRecentRadarImages(String radarStation){
    RadarImageUrls radarUrls = new RadarImageUrls(radarStation, numImages);
    WebPageFetcher fetcher = new WebPageFetcher();
    String dirListing = fetcher.fetch(radarUrls.dirListingUrl(), FetchFromThirdParty.ENCODING);
    List<String> results = radarUrls.findImageUrls(dirListing);
    return buildJsonStringForRadar(results, numImages);
  }
  
  private Integer numImages;
  private static final Logger fLogger = Util.getLogger(Radar.class);
  

  //var jsonString = '{"radarUrls" : ["http://localhost:8081/playtime/fetch/?radarStation=' + radar_station + '"]}';
  private String buildJsonStringForRadar(List<String> radarUrls, Integer limitTo){
    fLogger.info("Radar URLs: " + radarUrls);
    StringBuilder result = new StringBuilder();
    result.append("{\"radarUrls\": [");
    int idx = 0;
    while (idx < limitTo){
      result.append("\"" + radarUrls.get(idx) + "\"");
      if (idx < (limitTo-1)){
        result.append(",");
      }
      ++idx;
    }
    result.append("]}");
    return result.toString();
  }
}
