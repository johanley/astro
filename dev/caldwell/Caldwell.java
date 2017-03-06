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

/** 
Produce text that can be used as javascript, to represent 
the contents of the Caldwell catalog.

<P>The coords are J2000. I'm using only those items below -35 deg declination, 
to complement the Messier catalog for the southern hemisphere.

Example command lines.
Note the encoding flag.
 %astro%\dev\caldwell>C:\jdk1.6.0\bin\javac -cp . -encoding UTF-8 Caldwell.java
 %astro%\dev\caldwell>C:\jdk1.6.0\bin\java -cp . Caldwell
*/
public final class Caldwell {

  public static void main(String... args) throws IOException {
    log("Starting...");
    Caldwell caldwell = new Caldwell();
    caldwell.run(INPUT_FILE_NAME, OUTPUT_FILE_NAME);
    log("Finished");
  }

  private static final String INPUT_FILE_NAME = "caldwell_raw.txt";
  private static final String OUTPUT_FILE_NAME = "caldwell.js";
  
  private static final String SEP = ",";
  private static final String Q = "\""; // as in 'quote'
  private static final String COMMA = ",";
  private static final String ENCODING = "UTF-8";
  private static final String NL = System.getProperty("line.separator");

  private static final class Nebula {
    Integer INDEX; //0..N
    String NAME; //never empty; starts at 68, since I'm taking only those below -35 declination
    Double RA; /*rads J2000 */
    Double DEC; /*rads J2000 */
    String MAG; 
    String CONSTELLATION;
    String TYPE = "";  
    String COMMON_NAME; //sparse - from the main source data
    String LINK; // needed for caldwells only, where the wikipedia link is not computable
   /** Formatted as a javascript array. */
   @Override public String toString(){
      String sep = ",";
      return "caldwell["+INDEX+"]=[" + 
          Q+NAME+Q+sep+
          RA+sep+
          DEC+sep+
          MAG+sep+
          Q+CONSTELLATION+Q+sep+
          Q+TYPE+Q+sep+
          Q+""+Q+sep+           // the comment is always empty; just doing this to match the Messier structure
          Q+COMMON_NAME+Q+sep+ 
          Q+LINK+Q+  //applicable to caldwells, but not to messiers 
       "];";
    }
  }

  void run(String aFileName, String aOutputFileName) throws IOException {
    log("Reading source data file: " + aFileName);
    log("Assumed encoding: " + ENCODING);
    Scanner scanner = new Scanner(new FileInputStream(aFileName), ENCODING);
    List<Nebula> nebulae = new ArrayList<Nebula>();
    Nebula nebula = null;
    int lineCount = 0;
    try {
      log("Scanning records.");
      while (scanner.hasNextLine()){
        ++lineCount;
        String line = scanner.nextLine();
        log(line);
        if (line.startsWith("#")){
          continue;
        }
        else {
          nebula = processLine(line, nebulae.size());
          nebulae.add(nebula);
        }
      }
    }
    finally {
      scanner.close();
    }
    log("Read this many lines: " + lineCount);
    log("Number of records: " + nebulae.size());
    finalOutput(nebulae, OUTPUT_FILE_NAME);
  }
  
  /*
   Id|NCG etc|Type|Mag|Size|Distance|Ra|Dec|Constellation|Season|Common Name|Link (wikipedia)
   68|6729|Bn|9.7|1.0|424|19 01.9|-36 57|CrA|summer  
  */
  private Nebula processLine(String line, int nebulaCount){
    Nebula result = new Nebula();
    result.INDEX = nebulaCount;
    result.NAME = "C" + part(line, 1); 
    result.MAG = part(line, 4);
    result.CONSTELLATION = part(line, 9);
    result.TYPE = part(line, 3); 
    result.RA = rightAscension(part(line, 7));
    result.DEC = declination(part(line, 8));
    //no comment is present in the case of Caldwell objects
    result.COMMON_NAME = part(line, 11); //not always present, can be blank
    result.LINK = part(line, 12);
    return result;
  }
  
  private String part(String line, int idx /*1-based*/){
    String result = "";
    if (line!=null && line.trim().length()>0){
      String[] parts = line.split(SEP);
      /*
      int count = 0;
      for(String part: parts){
        ++count;
        System.out.print(count + ":"+ part + ", ");
      }
      */
      if (parts.length >= idx){
        result = parts[idx-1];
      }
    }
    //log(""); //new line
    //log("Part"+idx + ":" + result.trim());
    return result.trim();
  }
  
  /** The last char is a unit. */
  private String removeUnit(String num){
    return num.substring(0, num.length()-1);
  }
  
  private Integer partInt(String line, int idx /*1-based*/){
    return Integer.valueOf(removeUnit(part(line, idx)));
  }
  
  private Double partDbl(String line, int idx /*1-based*/){
    return Double.valueOf(removeUnit(part(line, idx)));
  }
  
  /** 19 01.9, leading 0s */
  private Double rightAscension(String raw){
    log("RA: " + raw);
    int space = raw.indexOf(" ");
    String hour_raw = raw.substring(0, space);
    String min_raw = raw.substring(space+1);
    int hour = Integer.valueOf(hour_raw);
    double minutes = Double.valueOf(min_raw);
    return round(rads((hour + minutes/60.0) * 15)); //no integer div
  }
  
  /** -36 07, leading 0s */
  private Double declination(String raw){
    int space = raw.indexOf(" ");
    int sign = raw.startsWith("+") ? 1 : -1;
    String deg_raw = raw.substring(1, space); 
    String min_raw = raw.substring(space+1); 
    int deg = Integer.valueOf(deg_raw);
    int minutes = Integer.valueOf(min_raw);
    return sign * round(rads(deg + minutes/60.0)); //no integer div
  }
  
  private double rads(double degs){
    //avoid integer division!
    return degs * Math.PI * 2 / 360.0D;
  }
  
  /** 1 arc sec is 5x10-6 rads. */
  private double round(double rads){
    double numdecimals = 10000000.0D;
    return Math.round(rads*numdecimals)/numdecimals; //avoid int division
  }
  
  /** Tab separated, type + comment, both in double quotes. */
  private void processLineComment(String line, Nebula nebula /* out-param, changed in place */){
    nebula.TYPE = line.substring(0,2).trim();
  }
  
  private void finalOutput(List<Nebula> nebulae, String filename) throws FileNotFoundException, IOException {
    log("Writing to file: " + filename);
    File out = new File(filename);
    FileOutputStream fos = new FileOutputStream(out);
    BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(fos, ENCODING));
    writer.write("/* Caldwell catalog, below -35 deg declination. Generated on: " + new Date() + ". Name, Right Ascension (J2000), Declination (J2000), Mag, Constellation, Type, Comment (blank!), and Common Name.*/");
    writer.newLine();
    writer.write("var caldwell = [" + nebulae.size() + "];");
    writer.newLine();
    for(Nebula nebula: nebulae){
      writer.write(nebula.toString());
      writer.newLine();
    }
    writer.newLine();
    writer.close();    
  }
  
  private static void log(Object text){
    System.out.println(text.toString());
  }
}
