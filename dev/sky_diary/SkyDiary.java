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


/**
 Generate an array of items describing the date and time of celestial events, 
 for display on a calendar or in a listing.
 
 Example input:
<pre>
                                                   d  h
New Moon                                 2016 Oct  1 00    401579.717 km
Spica in conjunction with Moon           2016 Oct  2 04    5.78° South
</pre>
Input is generated using MICA software (US Naval Observatory).
  Calc -> Search phenomena -> full year + 1 month, use UT1
  Export as text; to be safe, open in Notepad and re-save as UTF-8. 
  (All lines to the nearest hour only; no choices)
  Discarded by this class:
    many lines with asteroid phenomena
    duplicates: X in conjunction with Y, at the same time 

 Example output, a js object (plus a comma):
     {when: 'UT 2016-03-01 23:11', text:'Last quarter Moon'},
  
<P>Compile and run the old fashioned way:
  <pre>
%astro%\dev\sky_diary>C:\jdk1.6.0\bin\javac -encoding "UTF-8" -cp . SkyDiary.java
%astro%\dev\sky_diary>C:\jdk1.6.0\bin\java SkyDiary
</pre>
*/
public class SkyDiary {
  
  /** Run the script, output to a file. */
  public static void main(String... args)  throws IOException {
    log("Starting the script.");
    log("Input file (UTF-8), a MICA export of phenomena search: " + INPUT_FILE_MICA_EXPORT);
    SkyDiary script = new SkyDiary();
    script.parseAndGenerate();
    log("Done. Output (UTF-8) is in: " + OUTPUT_FILE);
  }
  
  public static final String INPUT_FILE_MICA_EXPORT = "phenomena_utf8.txt";
  public static final String OUTPUT_FILE = "phenomena.js";

  // PRIVATE
  
  private static final String NL = System.getProperty("line.separator");
  private static final String WS = "\\s+";
  private static final List<String> ASTEROIDS = Arrays.asList("Ceres", "Cybele",  "Davida",  "Eunomia",  "Europa",  "Flora",  "Hebe",  "Hygiea",  "Interamnia",  "Iris",  "Juno",  "Metis",  "Pallas",  "Psyche",  "Vesta", "Pluto");

  private static final class Event {
    String UT = "";
    String text = "";
    public String toString(){
      // {when: 'UT 2016-03-01 23:11', text:'Last quarter Moon'}
      return "{when:'"+UT+"',text:'"+text+"'}";
    }
  }
  
  private void parseAndGenerate() throws IOException{
    File file = new File(INPUT_FILE_MICA_EXPORT);
    String result = "";
    List<Event> events = new ArrayList<Event>();
    List<String> lines = read(file);
    lines = lines.subList(6, lines.size()); //discard header lines
    for(String line : lines){
      if (line.trim().length() > 0){ //has content
        Event prevEvent = events.isEmpty() ? null : events.get(events.size()-1);
        Event event = parseEvent(line, prevEvent);
        if (event != null){
          events.add(event);
        }
      }
    }
    log("Num events: " + events.size());
    for (Event event : events){
      result = result + "     " + event.toString() + "," + NL; 
    }
    output(result);
  }
  
  private List<String> read(File file) throws IOException {
    List<String> result = new ArrayList<String>();
    Scanner scanner = new Scanner(new FileInputStream(file), "UTF-8");
    try {
      while (scanner.hasNextLine()){
        result.add(scanner.nextLine());
      }
    }
    finally{
      scanner.close();
    }
    return result;
  }
  
  private static final Map<String, String> MONTH = new LinkedHashMap<String, String>();
  static {
    MONTH.put("Jan", "01");
    MONTH.put("Feb", "02");
    MONTH.put("Mar", "03");
    MONTH.put("Apr", "04");
    MONTH.put("May", "05");
    MONTH.put("Jun", "06");
    MONTH.put("Jul", "07");
    MONTH.put("Aug", "08");
    MONTH.put("Sep", "09");
    MONTH.put("Oct", "10");
    MONTH.put("Nov", "11");
    MONTH.put("Dec", "12");
  }

  /** Can return null. */
  private Event parseEvent(String line, Event previousEvent){
    Event result = null;
    //Spica in conjunction with Moon           2016 Oct  2 04    5.78° South
    //log("      Event line: " + line);
    String text1 = line.substring(0, 41).trim(); //Spica in conjunction with Moon        
    String text2 = line.substring(55).trim(); //   5.78° South
    String date = line.substring(41,52); //2016 Oct  2
    String hour = line.substring(53,55); // 04
    
    if (hasAsteroids(text1)){
      //abandon the line
    }
    else {
      result = new Event();
      String[] parts = date.split(WS);
      result.UT = "UT " + parts[0] + "-" + MONTH.get(parts[1]) + "-" + pad0(parts[2]) + " " + hour ; 
      result.text = rephraseConjunctions(text1, text2).trim();
      if (isDuplicateOf(previousEvent, result)){
        //abandon the duplicate
        result = null;
      }
    }
    return result;
  }
  
  private boolean hasAsteroids(String text1){
    boolean result = false;
    for(String asteroid : ASTEROIDS){
      if (text1.contains(asteroid)){
        result = true;
        break;
      }
    }
    return result;
  }
  
  private static final String CONJUNCTION = " in conjunction with ";
  
  private boolean isDuplicateOf(Event prevEvent, Event event){
    boolean result = false;
    if (prevEvent != null && prevEvent.UT.equals(event.UT)){
      result = true;
      log("  Found duplicate of: " + prevEvent.text + " " + prevEvent.UT);
    }
    return result;
  }
  
  /** Only check the first and last parts; ignore the middle parts. Not 100%, but in practice ok.*/
  private boolean isReverseOf(String a, String b){
    String[] partsA = a.split(WS);
    String[] partsB = b.split(WS);
    return (partsA[0].equals(partsB[partsB.length-1]) && partsA[partsA.length-1].equals(partsB[0]));
  }

  private String pad0(String thing){
    String result = thing;
    if (thing.length() == 1){
      result = "0" + result;
    }
    return result;
  }
  
  private String rephraseConjunctions(String text1, String text2){
    String result = text1 + " " + text2; 
    if (text1.contains(CONJUNCTION)){
      String[] parts = text1.split(WS);
      String[] degreeParts = text2.split(WS);
      result = parts[0] + " " + degreeParts[0]  + degreeParts[1].substring(0, 1) + " of " + parts[parts.length-1];
    }
    return result;
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
