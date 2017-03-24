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
  Script to generate links to satellite images, for various stations.
   
  This is meant as a one-time operation; not used at runtime. Dev tool only.
   
  <P>Such links save time for the user, since they don't have to figure out the form.
  The base data is in a CSV text file.
  
  CAREFUL: Eclipse will change tabs into spaces randomly when you resave in Eclipse. 
  That file is read in, and a link is generated for each line in the file.

  <P>Compile and run like so:
  <pre>
%astro%\dev\sat_image_links>C:\jdk1.6.0\bin\javac -cp . SatImageLinkGenerator.java

%astro%\dev\station_links>C:\jdk1.6.0\bin\java SatImageLinkGenerator 
  </pre>
  
  As a secondary effect, this class also generates text for translations.txt, using translations in 
  the same .csv source file.
 */
public class SatImageLinkGenerator implements ClipboardOwner {
  
  /** Run the script. Output is to the clipboard. */
  public static void main(String... args) throws FileNotFoundException, UnsupportedEncodingException {
    log("Starting the script.");
    
    SatImageLinkGenerator script = new SatImageLinkGenerator();
    script.read_csv_and_build_links(FILE_NAME, ENCODING);
    
    log("Done. Output is on the clipboard.");
  }
  
  public static final String FILE_NAME = "stations.csv"; 
  public static final String ENCODING = "UTF-8";
  public static final String FIELD_SEP = "\t";
  
  /** Empty impl. */
  @Override public void lostOwnership(Clipboard aClipboard, Transferable aContents){
    //do nothing
  }
  
  /**
   Place a String on the clipboard, and make this class the
   owner of the Clipboard's contents.
  */
  public void setClipboardContents(String aString){
    StringSelection stringSelection = new StringSelection(aString);
    Clipboard clipboard = Toolkit.getDefaultToolkit().getSystemClipboard();
    clipboard.setContents(stringSelection, this);
  }
  
  // PRIVATE
  
  // Example URL for viewing a large satellite image: 
  //  http://localhost:8081/astro/satellite/graphic.sky?latitude=47&longitude=-64&degrees_on_a_side=10&pixels_on_a_side=800&layer=auto_detect&locations=46.25%2C+-63.13
  
  /** 
    The start of the url.
    
    Only the lat-long is added dynamically to the URL. 
    The data file has the name of the place as well. 
    The name varies in exact meaning; it could be a country, it could be a region of a country. 
  */
  private static final String STATIC_URL = "graphic.sky?layer=auto_detect&degrees_on_a_side=15&pixels_on_a_side=&"; 
  
  private void read_csv_and_build_links(String fileName, String encoding) throws FileNotFoundException, UnsupportedEncodingException {
    log("Reading from file: " + fileName);
    StringBuilder text = new StringBuilder();
    StringBuilder translations = new StringBuilder();
    String NL = System.getProperty("line.separator");
    String line = "";
    Scanner scanner = new Scanner(new FileInputStream(fileName), encoding);
    try {
      while (scanner.hasNextLine()){
        line = scanner.nextLine().trim();
        if (line != null && line.length() > 0){
          text.append(processLineForOptionInSelect(line) + NL);
          translations.append(processLineForTranslations(line) + NL);
        }
      }
    }
    finally{
      scanner.close();
    }
    log(text);
    setClipboardContents(text.toString() + NL + translations.toString());
  }

  private String processLineForOptionInSelect(String line) {
    String result = STATIC_URL;
    String[] parts = line.split(FIELD_SEP);
    List<String> fields = Arrays.asList(parts);
    //System.out.println("Num fields = " + fields.size());
    log(fields);
    
    String name = fields.get(0).trim();
    String lat = round(fields.get(1).trim());
    String longit = round(fields.get(2).trim());
    
    result = result + "latitude=" + lat + "&";
    result = result + "longitude=" + longit;
    
    //add the markup
    result = "     <option value='" + result + "'><s:txt>" + name + "</s:txt>";
    return result;
  }
  
  private String processLineForTranslations(String line) {
    String[] parts = line.split(FIELD_SEP);
    List<String> fields = Arrays.asList(parts);
    
    String name = fields.get(0).trim();
    String nameFrench = fields.get(3).trim(); 
    
    return "{" + name + "}{" + nameFrench + "}";
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
}
