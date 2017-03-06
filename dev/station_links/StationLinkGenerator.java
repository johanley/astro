import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.UnsupportedEncodingException;
import java.math.BigDecimal;
import java.math.RoundingMode;
import java.net.URLEncoder;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.LinkedHashMap;
import java.util.Scanner;
import java.awt.datatransfer.Clipboard;
import java.awt.datatransfer.ClipboardOwner;
import java.awt.datatransfer.Transferable;
import java.awt.datatransfer.StringSelection;
import java.awt.datatransfer.DataFlavor;
import java.awt.datatransfer.UnsupportedFlavorException;
import java.awt.Toolkit;

/** 
  Script to generate links for various stations.
   
  This is meant as a one-time operation; not used at runtime. Dev tool only.
   
  <P>Such links save time for the user, since they don't have to figure out the form.
  The base data is in a CSV text file.
  CAREFUL: Eclipse will change tabs into spaces randomly when you resave in Eclipse. 
  That file is read in, and a link is generated for each line in the file.

  <P>Compile and run the old fashioned way:
  <pre>
%astro%\dev\station_links>C:\jdk1.6.0\bin\javac -cp . StationLinkGenerator.java

%astro%\dev\station_links>C:\jdk1.6.0\bin\java StationLinkGenerator cda
%astro%\dev\station_links>C:\jdk1.6.0\bin\java StationLinkGenerator us
%astro%\dev\station_links>C:\jdk1.6.0\bin\java StationLinkGenerator uk
  </pre>
 */
public class StationLinkGenerator implements ClipboardOwner {
  
  /** Run the script. Output is to the clipboard. */
  public static void main(String... args) throws FileNotFoundException, UnsupportedEncodingException {
    log("Starting the script.");
    
    String inputFile = inputFileName(args[0]);
    StationLinkGenerator script = new StationLinkGenerator();
    script.read_csv_and_build_links(inputFile, ENCODING);
    
    log("Done. Output is on the clipboard.");
  }
  
  public static final String FILE_NAME_ENDS_WITH = "_links.csv"; //and starts with the country name/identifier
  public static final String ENCODING = "UTF-8";
  public static final String FIELD_SEP = "\t";
  
  public static final String CANADA = "cda";
  public static final String US = "us";
  public static final String UK = "uk";
 
  /** Empty impl. */
  @Override public void lostOwnership(Clipboard aClipboard, Transferable aContents){
    //do nothing
  }
  
  /**
  * Place a String on the clipboard, and make this class the
  * owner of the Clipboard's contents.
  */
  public void setClipboardContents(String aString){
    StringSelection stringSelection = new StringSelection(aString);
    Clipboard clipboard = Toolkit.getDefaultToolkit().getSystemClipboard();
    clipboard.setContents(stringSelection, this);
  }
  
  // PRIVATE
  
  /** The start of the url. */
  private static final Map<String, String> STATIC_URL = new LinkedHashMap<String, String>();
  private static void addCountryStaticUrl(String country){
    STATIC_URL.put(country, "graphic.sky?country=" + country + "&");
  }
  static {
    addCountryStaticUrl(CANADA);
    addCountryStaticUrl(US);
    addCountryStaticUrl(UK);
  }
  
  private static String inputFileName(String country){
      return country +  FILE_NAME_ENDS_WITH;
  }
  
  private void read_csv_and_build_links(String fileName, String encoding) throws FileNotFoundException, UnsupportedEncodingException {
    log("Reading from file: " + fileName);
    int underscoreIdx = fileName.indexOf(FILE_NAME_ENDS_WITH);
    int lastSlashIdx = fileName.lastIndexOf(File.separator);
    String country = fileName.substring(lastSlashIdx + 1, underscoreIdx);
    log("Country: " + country);
    StringBuilder text = new StringBuilder();
    String NL = System.getProperty("line.separator");
    String line = "";
    Scanner scanner = new Scanner(new FileInputStream(fileName), encoding);
    try {
      while (scanner.hasNextLine()){
        line = scanner.nextLine().trim();
        if (line != null && line.length() > 0){
          text.append(processLineForOptionInSelect(line, country) + NL);
        }
      }
    }
    finally{
      scanner.close();
    }
    log(text);
    setClipboardContents(text.toString());
  }

  private String processLineForOptionInSelect(String line, String country) throws UnsupportedEncodingException{
    String result = "";
    if (country.equals(CANADA) || country.equals(US)){
       result = northAmerica(line, country);      
    }
    else if (country.equals(UK)) {
      result = unitedKingdom(line);
    }
    return result;
  }
  
  private String northAmerica(String line, String country) throws UnsupportedEncodingException{
    String result = STATIC_URL.get(country); //start of the url
    Boolean isCanada = country.equals(CANADA);
    Boolean isUS = country.equals(US);
    String[] parts = line.split(FIELD_SEP);
    List<String> fields = Arrays.asList(parts);
    //System.out.println("Num fields = " + fields.size());
    log(fields);
    
    String prov = fields.get(0).trim();
    String name = fields.get(1).trim();
    String lat_and_long = fields.get(2).trim();
    String[] lat_long_parts = lat_and_long.split(",");
    String lat = round(lat_long_parts[0].trim());
    String longit = round(lat_long_parts[1].trim());
    String weather_stn_id = fields.get(3).trim(); //always empty for US, where weather data comes from lat and long
    String clear_sky_clock_url = fields.get(4).trim();
    String radar = fields.get(5).trim();
    if (isBlank(radar)){ //no radar stations in the north
      radar = "";
    }
    
    //System.out.println("Satt: " + fields.get(6).trim());
    boolean hasNoSatelliteImage = fields.get(6).trim().equals("x"); //no coverage in the north
    if (hasNoSatelliteImage){
      result = result.replace("include_clouds=1", "include_clouds=0");
    }
    
    if (isCanada){
      result = result + "prov=" + prov + "&";
    }
    result = result + "location_name=" + name + "&";
    result = result + "latitude=" + lat + "&";
    result = result + "longitude=" + longit + "&";
    result = result + "locations=" + lat + "," + longit + "&";
    if (!isUS){
      result = result + "weather_station=" + weather_stn_id + "&"; //blank if missing
    }
    result = result + "radar_station=" + radar + "&";
    result = result + "clear_sky_clock_station_id=" + extractStationIdFrom(clear_sky_clock_url, isCanada);
    
    result = escapeSpecialChars(result);
    
    //add the markup
    result = " <option value='" + result + "'>" + prov + " " + name;
    
    return result;
  }
  
  private String extractStationIdFrom(String cskUrl, Boolean isCanada){
    //http://www.cleardarksky.com/c/StJohnscsk.gif?c=1418437 
    String result = "";
    String start = "http://www.cleardarksky.com/c/"; 
    String end = isCanada ? "csk.gif?" : "cs0.gif?";
    int endIdx = cskUrl.indexOf(end);
    return cskUrl.substring(start.length(), endIdx);
 }

  private String escapeSpecialChars(String url){
    String result = new String(url);
    //result = URLEncoder.encode(result, ENCODING);
    result = result.replace(" ", "+");
    result = result.replace("'", "%27"); //St. John's
    result = result.replace(".", "%2E"); //St. John's
    return result;
  }
  
  private String unitedKingdom(String line) throws UnsupportedEncodingException{
    String result = STATIC_URL.get(UK); //start of the url
    String[] parts = line.split(FIELD_SEP);
    List<String> fields = Arrays.asList(parts);
    //System.out.println("Num fields = " + fields.size());
    log(fields);
    
    String name = fields.get(0).trim();
    String lat = round(fields.get(1).trim());
    String longit = round(fields.get(2).trim());
    String weather_stn_id = fields.get(3).trim(); 
    
    result = result + "location_name=" + name + "&";
    result = result + "latitude=" + lat + "&";
    result = result + "longitude=" + longit + "&";
    result = result + "locations=" + lat + "," + longit + "&";
    result = result + "weather_station=" + weather_stn_id + "&"; 
    
    //escape special chars
    result = escapeSpecialChars(result);
    
    //add the markup
    result = " <option value='" + result + "'>" + name;
    
    return result;
  }
 
  private String round(String text){
    BigDecimal num = new BigDecimal(text);
    BigDecimal rounded = num.setScale(2, RoundingMode.HALF_EVEN);
    String result = rounded.toPlainString();
    log("in:" + text + " out" + result);
    return result;
  }
  
  private static void log(Object msg){
    System.out.println(msg.toString());
  }
  
  private boolean isBlank(String text){
    return text.equals("*");
  }
}
