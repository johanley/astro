import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Scanner;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
  Output a javascript data structure for phenomena 
  of Jupiter's Galilean satellites. 
  
 Input is from occult.exe 4: 
'                      Satellites of Jupiter 2016                       
                              October                                 
                                                                      
  1  2 18.5  1.Sh.I     11 17  9.9  1.Sh.I     21  9 28.4  3.Ec.D     
     2 23.0  1.Tr.I        17 25.0  1.Tr.I        10 14.7  2.Sh.I     
     .....
 '    
   each month is separate. No mutual phenomena are included.

 IS THIS IN TT OR UT? Not clear.
 
 Example output, a js object (plus a comma), rounded to the nearest minute:
     {when: 'UT 2016-03-01 23:11', satt: 1, ph:'Sh.I'},
   
   
<P>Compile and run the old fashioned way:
  <pre>
%astro%\dev\jupiter_satellite_phenom>C:\jdk1.6.0\bin\javac -encoding "UTF-8" -cp . JupiterSatellitePhenomena.java
%astro%\dev\jupiter_satellite_phenom>C:\jdk1.6.0\bin\java JupiterSatellitePhenomena
  </pre>
 */
public class JupiterSatellitePhenomena {
  
  /** Run the script, output to a file. */
  public static void main(String... args)  throws IOException {
    log("Starting the script.");
    log("Input file (UTF-8), output from occult.exe 4: " + INPUT_FILE);
    JupiterSatellitePhenomena script = new JupiterSatellitePhenomena();
    script.parseAndGenerate();
    log("Done. Output (UTF-8) is in: " + OUTPUT_FILE);
  }
  
  public static final String INPUT_FILE = "jupiter_satellite_phenomena_utf8.txt";
  public static final String OUTPUT_FILE = "jupiter_satellite_phenomena.js";

  // PRIVATE
  
  private static final String NL = System.getProperty("line.separator");
 
  private static final class Event {
    String UT = "";
    String ph = "";
    String satt = "";
    public String toString(){
      // {when: 'UT 2016-03-01 23:11:24',satt: 1,ph:'Sh.I'},
      return "{when:'"+UT+"',satt:"+satt+",ph:'"+ph+"'}";
    }
  }
  
  private static final Pattern allMonth = Pattern.compile("(\\d{4})(?:\\s)+(January|February|March|April|May|June|July|August|September|October|November|December)(.*)", Pattern.DOTALL);
  private void parseAndGenerate() throws IOException{
    File file = new File(INPUT_FILE);
    String result = "";
    String input = read(file);
    List<Event> events = new ArrayList<Event>();
    String[] months = input.split("Satellites of Jupiter "); // this chops up the data into separate months; i don't have to worry about greedy/reluctant 
    for (String rawMonth : months){
      if (hasContent(rawMonth)){
        Matcher matcher = allMonth.matcher(rawMonth);
        while (matcher.find()){
          log("Group count: " + matcher.groupCount());
          log("Group 0 length: " + matcher.group(0).length());
          log("Group 1: " + matcher.group(1));
          log("Group 2: " + matcher.group(2));
          String year = matcher.group(1); //2016
          String month = MONTH.get(matcher.group(2)); // 01..12
          String phenomena = matcher.group(3); // the rest of the data for the month
          //log(phenomena);
          log("Found match for year-month :" + year + "-" + month + " has this length for its monthly data: " + phenomena.length());
          parse(phenomena, year, month, events);
        }
      }
    }
    for (Event event : events){
      result = result + event.toString() + "," + NL;
    }
    output(result);
  }
  
  private String read(File file) throws IOException {
    StringBuilder result = new StringBuilder("");
    Scanner scanner = new Scanner(new FileInputStream(file), "UTF-8");
    try {
      while (scanner.hasNextLine()){
        result.append(scanner.nextLine() + NL);
      }
    }
    finally{
      scanner.close();
    }
    return result.toString();
  }
  
  private static final Map<String, String> MONTH = new LinkedHashMap<String, String>();
  static {
    MONTH.put("January", "01");
    MONTH.put("February", "02");
    MONTH.put("March", "03");
    MONTH.put("April", "04");
    MONTH.put("May", "05");
    MONTH.put("June", "06");
    MONTH.put("July", "07");
    MONTH.put("August", "08");
    MONTH.put("September", "09");
    MONTH.put("October", "10");
    MONTH.put("November", "11");
    MONTH.put("December", "12");
  }
  
  private void parse(String phenomena, String year, String month, List<Event> events){
    List<String> lines = Arrays.asList(phenomena.split(NL));
    int RECORD_WIDTH = 23;
    // the data is in three columns; each column is the same
    // to get at the 3 columns, we iterate over the lines 3 times, each time with a different offset into the line
    for (int col = 0; col < 3; ++col){
      int startIdx = col * RECORD_WIDTH;
      String currentDay = ""; //the day is not repeated on each line, so we need to detect when a new day occurs
      for (String line : lines){
        if (!hasContent(line)) continue;
        String subLine = line.substring(startIdx, startIdx + RECORD_WIDTH-1);
        if (!hasContent(subLine)) continue;
        // pick out all the parts 
        String day = subLine.substring(0,3);
        if (hasContent(day)){
          currentDay = pad0(day.trim());
        }
        String hour = pad0(subLine.substring(4,6).trim());
        String min = pad0(subLine.substring(7,9).trim());
        String seconds = changeToSeconds(subLine.substring(10,11).trim()); // .1 is 06s
        String satt = subLine.substring(13,14);
        String ph = subLine.substring(15,19);
        Event event = new Event();
        event.ph = ph;
        event.satt = satt;
        event.UT = "UT " + year + "-" + month + "-" + currentDay + " " + hour + ":" + min + ":" + seconds; 
        events.add(event);
      }
    }
  }
  
  private boolean hasContent(String text){
    return text != null && text.trim().length() > 0;
  }
  
  private String pad0(String thing){
    String result = thing;
    if (thing.length() == 1){
      result = "0" + result;
    }
    return result;
  }
  
  private String changeToSeconds(String decimalMinutes){
    Integer decimal = Integer.valueOf(decimalMinutes);
    Integer seconds = 6 * decimal;
    return pad0(seconds.toString());
  }
  
  private void output(String result) throws IOException  {
    Writer out = new OutputStreamWriter(new FileOutputStream(OUTPUT_FILE), "UTF-8");
    try {
      out.write(result);
    }
    finally {
      out.close();
    }
  }
  
  private static void log(Object msg){
    System.out.println(msg.toString());
  }
  

}
