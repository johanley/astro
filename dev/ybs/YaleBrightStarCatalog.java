
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.List;
import java.util.Scanner;
import java.util.LinkedHashMap;
import java.util.Map;

/** 
 Produce text that can be used as javascript, to represent 
 the stars in the Yale Bright Star catalog.
 Some items in the catalog are discarded (novae, non-stellar objects).
 
 NOTE: the output is for J2000, the same as the source catalog.
 I have a tool that precesses all the data to another equinox; it's implemented 
 in JS (see ybs_precess). When you run it, run it in 2 parts, or the browser will hang.
 
 <P>File Format:
 http://cdsarc.u-strasbg.fr/viz-bin/Cat?V/50
 
 <P>The coords are J2000.

 Example command lines.
 Note the encoding flag.
  %astro%\dev\ybs>C:\jdk1.6.0\bin\javac -cp . -encoding UTF-8 YaleBrightStarCatalog.java  
  %astro%\dev\ybs>C:\jdk1.6.0\bin\java -cp . YaleBrightStarCatalog
  
  Nova?: this star has mag 2, in Corona Borealis, and has no bayer designation. 
  ybs[5948]=['',4.1866231,0.4523942,2.0];
  to get rid of it, i've fudged the magnitude
  
It's this entry (T CrB):
  5958          BD+26 2765 143454 84129    I W       T CrB    155519.0+261213155930.2+255513 42.37 48.16 2.0 H +0.1        +1.56 sdBe+gM3+Q           -0.005+0.013      -029SBO              0.2      *
  
Problem α Psc is being treated as too dim. It has two dim entries in the source data.
  ybs[591]=['α Psc',0.532529,0.0482341,5.23];
  ybs[592]=['α Psc',0.532529,0.0482341,4.33];
  Bad:  
    visually, this star is mag really, 3.82
    tools pick up the first entry only, which is dim!
  Stupid tweak: change the js data manually, for the first entry, to mag 3.82.
*/
public final class YaleBrightStarCatalog {

  public static void main(String... args) throws IOException {
    log("Starting...");
    YaleBrightStarCatalog ybs = new YaleBrightStarCatalog();
    ybs.run(INPUT_FILE_NAME, OUTPUT_FILE_NAME);
    log("Finished");
  }
  
  private static final String INPUT_FILE_NAME = "yale_bright_star_catalog_5_raw.txt";
  private static final String OUTPUT_FILE_NAME = "yale_bright_stars.js";
  
  private static final String SEP = "\t";
  private static final String COMMA = ",";
  private static final String ENCODING = "UTF-8";
  private static final String NL = System.getProperty("line.separator");
  private static final List<Integer> SKIPPED_LINES = Arrays.asList(92,95,182,1057,1841,2472,2496,3515,3671,6309,6515,7189,7539,8296);

  private static final class Star {
    Integer INDEX;
    String NAME; //possibly empty, never null; first take bayer; it not present, take flamsteed
    Double RA; /*rads J2000 */
    Double DEC; /*rads J2000 */
    String MAG; 
   /** Formatted as a javascript array. */
   @Override public String toString(){
      String sep = ",";
      return "ybs["+INDEX+"]=['"+NAME+"'"+sep+RA+sep+DEC+sep+MAG +"];";
    }
  }

  void run(String aFileName, String aOutputFileName) throws IOException {
    log("Reading source data file: " + aFileName);
    log("Assumed encoding: " + ENCODING);
    Scanner scanner = new Scanner(new FileInputStream(aFileName), ENCODING);
    List<Star> stars = new ArrayList<Star>();
    Star star = null;
    int lineCount = 0;
    try {
      log("Scanning records.");
      while (scanner.hasNextLine()){
        ++lineCount;
        String line = scanner.nextLine();
        //log(line);
        if (SKIPPED_LINES.contains(Integer.valueOf(lineCount))){
          continue;
        }
        else {
          star = processLine(line, stars.size());
          stars.add(star);
        }
      }
    }
    finally {
      scanner.close();
    }
    log("Read this many lines: " + lineCount);
    log("Number of records: " + stars.size());
    finalOutput(stars, OUTPUT_FILE_NAME);
  }  

  /**
   *  1st line is:
   *    '', 6.7, 0.022536564 (ra), 0.7893978763 (dec)
   *  
   *  15th line is: '
   *    'Alp And', 2.06, 0.0366010089 (ra), 0.5077259757 (dec)
   */
  private Star processLine(String line, int starCount){
    Star result = new Star();
    result.INDEX = starCount;
    //prefer bayer to flamsteed
    result.NAME = bayerDesignation(slice(line, 8, 7)); //possibly empty
    if (isEmpty(result.NAME)){
      result.NAME = flamsteedDesignation(slice(line, 5, 10)); //possibly empty
    }
    result.MAG = slice(line, 103, 5); //possible leading minus sign; that's ok
    
    int ra_hour = sliceInt(line, 76, 2); //leading 0's for these
    int ra_min = sliceInt(line, 78, 2);
    double ra_sec = sliceDbl(line, 80, 4); 
    result.RA = round(rads((ra_hour + ra_min/60.0 + ra_sec/3600.0) * 15)); //no integer div
    
    int sign = slice(line, 84, 1).equals("+") ? 1 : -1;
    int dec_deg = sliceInt(line, 85, 2); //leading 0's for these
    int dec_min = sliceInt(line, 87, 2);
    int dec_sec = sliceInt(line, 89, 2);
    result.DEC = sign * round(rads(dec_deg + dec_min/60.0 + dec_sec/3600.0)); //no integer div
    
    return result;
  }
  
  private String slice(String line, int start /*1-based*/, int numchars){
    return line.substring(start-1, start-1+numchars).trim();
  }
  
  private Integer sliceInt(String line, int start, int numchars){
    return Integer.valueOf(slice(line, start, numchars));
  }
  
  private Double sliceDbl(String line, int start, int numchars){
    return Double.valueOf(slice(line, start, numchars));
  }
  
  private double rads(double degs){
    //avoid integer division!
    return degs * Math.PI * 2 / 360.0D;
  }
  
  private boolean isEmpty(String text){
    return text == null || text.trim().length() == 0; 
  }
  
  /** 1 arc sec is 5x10-6 rads. */
  private double round(double rads){
    double numdecimals = 10000000.0D;
    return Math.round(rads*numdecimals)/numdecimals; //avoid int division
  }
  
  /** Example input: 'Alp And', 'Kap1Scl'. Anything else is coerced to blank.  */
  private String bayerDesignation(String name){
    String result = "";
    if (name != null){
      if (name.length()>3){
        //it has both the Greek letter and the constellation name
        //the last 3 letters are the constellation abbr
        int len = name.length();
        String constellation = name.substring(len-3);
        String greekText = name.substring(0, len-3).trim();
        String greekLetter = greekLetter(greekText);
        result = greekLetter + " " + constellation; 
      }
    }
    return result;
  }
  
  /** Example input: '82    Psc'. Note the N spaces in the middle. Anything else is coerced to blank.  */
  private String flamsteedDesignation(String name){
    String result = "";
    if (name != null){
      if (name.length()>2){
        int len = name.length();
        String constellation = name.substring(len-3).trim();
        int firstBlank = name.indexOf(" ");
        String flamsteedNumber = name.substring(0, firstBlank).trim();
        result = flamsteedNumber + " " + constellation; 
      }
    }
    /*
    if (!isEmpty(result)){
      //get rid of the extra 
      String[] parts = result.split(" ");
      result = parts[0] + " " + parts[1]; 
    }
    */
    return result;
  }
  
  /** Translate Latin text abbreviations into Greek letters. */
  private static final Map<String, String> GREEK = new LinkedHashMap<String, String>();
  static {
    add("Alp", "α");
    add("Bet", "β");
    add("Gam", "γ");
    add("Del", "δ");
    add("Eps", "ε");
    add("Zet", "ζ");
    add("Eta", "η");
    add("The", "θ");
    add("Iot", "ι");
    add("Kap", "κ");
    add("Lam", "λ");
    add("Mu", "μ");
    add("Nu", "ν");
    add("Xi", "ξ");
    add("Omi", "ο");
    add("Pi", "π");
    add("Rho", "ρ");
    add("Sig", "σ");
    add("Tau", "τ");
    add("Ups", "υ");
    add("Phi", "φ");
    add("Chi", "χ");
    add("Psi", "ψ");
    add("Ome", "ω");
    
  }
  private static void add(String in, String out){
    GREEK.put(in, out);
  }
  private String greekLetter(String text /*Alp2, for example*/){
    String input = text;
    int len = text.length();
    char lastChar = text.charAt(len-1);
    if (Character.isDigit(lastChar)){
      input = text.substring(0, len-1); //without the number
    }
    String output = GREEK.get(input.trim());
    if (output == null){
      log("Greek letter not found for: '"  + input + "'");
    }
    if (Character.isDigit(lastChar)){
      output = output + lastChar;
    }
    return output;
  }
  
  private void finalOutput(List<Star> brightstars, String filename) throws FileNotFoundException, IOException {
    log("Writing to file: " + filename);
    File out = new File(filename);
    FileOutputStream fos = new FileOutputStream(out);
    BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(fos, ENCODING));
    writer.write("/* Source: Yale Bright Star Catalog r5. J2000. Generated on: " + new Date() + ". Name, Right Ascension, Declination, and Magnitude.*/");
    writer.newLine();
    writer.write("var ybs = [" + brightstars.size() + "];");
    writer.newLine();
    int count = 0;
    for(Star nearbystar: brightstars){
      writer.write(nearbystar.toString());
      if (count < (brightstars.size() - 1)){
        //writer.write(COMMA);
      }
      writer.newLine();
      ++count;
    }
    writer.newLine();
    writer.close();    
  }
  
  private static void log(Object text){
    System.out.println(text.toString());
  }
}
