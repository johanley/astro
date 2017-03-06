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
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Scanner;

/** 
Produce text that can be used as javascript, to represent 
the contents of the Messier catalog.

<P>The coords are J2000.

Example command lines.
Note the encoding flag.
 %astro%\dev\messier>C:\jdk1.6.0\bin\javac -cp . -encoding UTF-8 Messier.java
 %astro%\dev\messier>C:\jdk1.6.0\bin\java -cp . Messier
*/
public final class Messier {

  public static void main(String... args) throws IOException {
    log("Starting...");
    Messier messier = new Messier();
    messier.run(INPUT_FILE_NAME, INPUT_FILE_NAME_COMMENTS, OUTPUT_FILE_NAME);
    log("Finished");
  }

  //WARNING: the 2 input files need to have the same line structure
  private static final String INPUT_FILE_NAME = "messier_raw.txt";
  private static final String INPUT_FILE_NAME_COMMENTS = "messier_comment_and_type_raw.txt";
  
  private static final String OUTPUT_FILE_NAME = "messier.js";
  
  private static final String SEP = "\t";
  private static final String Q = "\""; // as in 'quote'
  private static final String COMMA = ",";
  private static final String ENCODING = "UTF-8";
  private static final String NL = System.getProperty("line.separator");
  private static final List<Integer> SKIPPED_LINES = Arrays.asList(100000); //EMPTY; an artifact of some other task

  private static final class Nebula {
    Integer INDEX;
    String NAME; //never empty
    Double RA; /*rads J2000 */
    Double DEC; /*rads J2000 */
    String MAG; 
    String CONSTELLATION;
    String TYPE = ""; //OC, GC, GY - from the Observer's Handbook 
    String DESCRIPTION = ""; //comment from the Observer's Handbook
    String COMMON_NAME; //sparse - from the main source data
   /** Formatted as a javascript array. */
   @Override public String toString(){
      String sep = ",";
      return "messier["+INDEX+"]=[" + 
          Q+NAME+Q+sep+
          RA+sep+
          DEC+sep+
          MAG+sep+
          Q+CONSTELLATION+Q+sep+
          Q+TYPE+Q+sep+
          Q+DESCRIPTION+Q+sep+
          Q+COMMON_NAME+Q + 
       "];";
    }
  }

  void run(String aFileName, String aFileNameComments, String aOutputFileName) throws IOException {
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
        if (line.startsWith("#") || SKIPPED_LINES.contains(Integer.valueOf(lineCount))){
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
    
    addCommentsAndType(aFileNameComments, nebulae);
    
    finalOutput(nebulae, OUTPUT_FILE_NAME);
  }
  
  private void addCommentsAndType(String aFileNameComments, List<Nebula> nebulae)  throws IOException {
    log("Reading source data file for comments and type info: " + aFileNameComments);
    log("Assumed encoding: " + ENCODING);
    Scanner scanner = new Scanner(new FileInputStream(aFileNameComments), ENCODING);
    Nebula nebula = null;
    int lineCount = 0;
    int recordCount = 0;
    try {
      log("Scanning records.");
      while (scanner.hasNextLine()){
        ++lineCount;
        String line = scanner.nextLine();
        log(line);
        if (line.startsWith("#") || SKIPPED_LINES.contains(Integer.valueOf(lineCount))){
          continue;
        }
        else {
          processLineComment(line, nebulae.get(recordCount));
          ++recordCount;
        }
      }
    }
    finally {
      scanner.close();
    }
    log("Read this many lines: " + lineCount);
    log("Read this many records: " + recordCount);
  }

  /*  
   M1  1952  Sn  8.4 6x4 6300  5h 34.5m  +22° 01′  Tau winter  Crab Nebula
   Tab separated fields
   name, ngc, type, mag, size, distance (ly), ra, decl, constellation, season, popular name 
  */
  private Nebula processLine(String line, int nebulaCount){
    Nebula result = new Nebula();
    result.INDEX = nebulaCount;
    result.NAME = part(line, 1); 
    result.MAG = part(line, 4); 
    result.CONSTELLATION = part(line, 9);
    result.DESCRIPTION = ""; // filled in manually from the Observer's Handbook
    result.TYPE = ""; // filled in manually from the Observer's Handbook
    result.RA = rightAscension(part(line, 7));
    result.DEC = declination(part(line, 8));
    result.COMMON_NAME = part(line, 11); //not always present
    return result;
  }
  
  private String part(String line, int idx /*1-based*/){
    String result = "";
    if (line!=null && line.trim().length()>0){
      String[] parts = line.split("\t");
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
  
  /** 5h 03.5m */
  private Double rightAscension(String raw){
    //log("RA: " + raw);
    int space = raw.indexOf(" ");
    String hour_raw = removeUnit(raw.substring(0, space));
    String min_raw = removeUnit(raw.substring(space+1)); //leading 0
    int hour = Integer.valueOf(hour_raw);
    double minutes = Double.valueOf(min_raw);
    return round(rads((hour + minutes/60.0) * 15)); //no integer div
  }
  
  /** +22° 01′ */
  private Double declination(String raw){
    int space = raw.indexOf(" ");
    int sign = raw.startsWith("+") ? 1 : -1;
    String deg_raw = removeUnit(raw.substring(1, space)); //remove leading sign, leading 0
    String min_raw = removeUnit(raw.substring(space+1)); //leading 0
    
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
    nebula.DESCRIPTION = line.substring(2).trim();
  }
  
  private void finalOutput(List<Nebula> nebulae, String filename) throws FileNotFoundException, IOException {
    log("Writing to file: " + filename);
    File out = new File(filename);
    FileOutputStream fos = new FileOutputStream(out);
    BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(fos, ENCODING));
    writer.write("/* Messier catalog.  Generated on: " + new Date() + ". Name, Right Ascension (J2000), Declination (J2000), Mag, and Common Name.*/");
    writer.newLine();
    writer.write("var messier = [" + nebulae.size() + "];");
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
