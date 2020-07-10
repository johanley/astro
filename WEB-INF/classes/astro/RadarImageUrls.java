package astro;

import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 Return a list of image URLs from a directory listing.
 
<P>Example directory listing:
  <pre>https://dd.weather.gc.ca/radar/PRECIPET/GIF/XNC/</pre>
  
<p>Example line of the directory listing: 
 <pre><img src="/icons/image2.gif" alt="[IMG]"> <a href="201911251500_XNC_COMP_PRECIPET_RAIN.gif">201911251500_XNC_COMP_PRECIPET_RAIN.gif</a>      2019-11-25 15:08   25K</pre>  
  
<P>This class doesn't fetch data from the web; it just processes the data that has been fetched.

<P>Info:
https://dd.weather.gc.ca/radar/doc/README_radar.txt
THE RADAR STATIONS ARE BEING UPGRADED OVER TIME. There is a change in station id, from 3 letters to 5 (CASxx).
WMN (Montreal) appears to be down on 2019-11-28.

<P>Broken for several months: change http to https, and added '_COMP_' in the middle of the image url.
Fixed 2019-11-27.
*/
public final class RadarImageUrls {
  
  /**
  * Constructor.
  * @param aRadarStationId eg 'XFT'
  * @param aNumImages the number of URLs to be returned; usually the caller needs just a few from the end.
  */
  public RadarImageUrls(String aRadarStationId, int aNumImages){
    fNumImages = aNumImages;
    DIRECTORY_URL = "https://dd.weather.gc.ca/radar/PRECIPET/GIF/" + aRadarStationId + "/";
    URL_PATTERN = Pattern.compile("(?:href=\")([\\d]{12})(?:_" + aRadarStationId + "_COMP_PRECIPET_RAIN.gif\")" );
    START_URL = "https://dd.weather.gc.ca/radar/PRECIPET/GIF/" + aRadarStationId + "/";
    END_URL = "_" + aRadarStationId + "_COMP_PRECIPET_RAIN.gif";
  }

  /** The images are found in a listing at this URL. */
  public String dirListingUrl(){
    return DIRECTORY_URL;
  }
  
  /** Return N URLs from the bottom of the directory listing. N is passed to the ctor.*/
  public List<String> findImageUrls(String directoryListing){
    if (fImageUrls.isEmpty()){
      fImageUrls = extractImageUrlsFrom(directoryListing);
    }
    return fImageUrls;
  }

  /** Extract the date-time of the image from its URL. */
  public String dateFrom(String aURL){
    int start = START_URL.length();
    int end = aURL.indexOf(END_URL);
    return aURL.substring(start, end);
  }
  
  // PRIVATE 

  private int fNumImages;
  private List<String> fImageUrls = new ArrayList<String>();
  private String DIRECTORY_URL = "";
  
  /** Used to extract image URLs from a raw directory listing. */
  private Pattern URL_PATTERN;
  /** The URL of the target image. */
  private String START_URL = "";
  private String END_URL = "";

  /**
   Example lines of the directory listing, of the form: 
   <img src="/icons/image2.gif" alt="[IMG]"> <a href="201911251500_XNC_COMP_PRECIPET_RAIN.gif">201911251500_XNC_COMP_PRECIPET_RAIN.gif</a>      2019-11-25 15:08   25K  
  
   Scan for all appearances of the above pattern, and create the corresponding URL. Return only the last few, 
   with the most recent as the last element.
   If the radar is not operational, then the returned list is empty.
  */
  private List<String> extractImageUrlsFrom(String aDirectoryListing){
    List<String> buffer = new ArrayList<String>();
    Matcher matcher = URL_PATTERN.matcher(aDirectoryListing);
    while (matcher.find()){
      String fileDate = matcher.group(1);
      buffer.add(urlFrom(fileDate));
    }
    List<String> result = new ArrayList<String>();
    if (!buffer.isEmpty()){
      //use only the last N items, and discard the rest; this impl is bulky, but simple to understand
      //the most recent is last in the list
      int end = buffer.size() - 1;
      int start = end - fNumImages + 1;
      for (int idx = start; idx <= end; ++idx){
        result.add(buffer.get(idx));
      }
    }
    return result;
  }
  
  private String urlFrom(String aDate){
    return START_URL + aDate + END_URL;
  }
}
