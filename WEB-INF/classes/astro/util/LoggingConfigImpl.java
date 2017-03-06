package astro.util;

import java.io.File;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;
import java.util.StringTokenizer;
import java.util.TimeZone;
import java.util.logging.FileHandler;
import java.util.logging.Handler;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.logging.SimpleFormatter;

import javax.servlet.ServletConfig;


/**
 Default implementation of {@link LoggingConfig}, to set up simple logging.
 
 <P>This implementation uses JDK logging, and appends logging output to a single file, 
 with no size limit on the file. It uses two settings in <tt>web.xml</tt> :
<ul>
 <li><tt>LoggingDirectory</tt> - the absolute directory which will hold the logging 
 output file. This class will always use a file name using the system date/time, as 
 returned by {@link DateTime#now(TimeZone)} using the <tt>DefaultUserTimeZone</tt> setting in 
 <tt>web.xml</tt>, in the form <tt>2007_12_31_59_59.txt</tt>. If the directory does not exist, WEB4J will 
 attempt to create it upon startup. If set to the special value of <tt>'NONE'</tt>, then 
 this class will not configure JDK logging in any way.
 <li><tt>LoggingLevels</tt> - a comma-separated list of logger names and their corresponding 
 levels. To verify operation, this class will emit test logging entries for each of these loggers, 
 at the stated logging levels. 
</ul>
*/
public final class LoggingConfigImpl {

  /** See class comment.   */
  public void setup(ServletConfig aConfig) {
    fLogger.config("Fetching logging settings from web.xml");
    fetchSettings(aConfig);
    if( isTurnedOff() ) {
      logStdOut("Default logging config is turned off, since directory is set to " + Util.quote(NONE));
    }
    else {
      logStdOut("Setting up logging config...");
      validateDirectorySetting();
      parseLoggers();
      createFileHandler();
      attachLoggersToFileHandler();
      tryTestMessages();
      fLogger.config("Logging to directory : " + Util.quote(fLoggingDir));
      showLoggerLevels();
    }
  }

  // PRIVATE 
  private String fLoggingDir;
  private String fLoggingLevels;
  private static final int NO_SIZE_LIMIT = 0;
  private static final int MAX_BYTES = NO_SIZE_LIMIT;
  private static final int NUM_FILES = 1;
  private static final boolean APPEND_TO_EXISTING = true;
  private static final String NONE = "NONE";
  private static final String SEPARATOR = "=";
  
  /** List of loggers. Each Logger stores its own Level as part of its state.  */
  private final List<Logger> fLoggers = new ArrayList<Logger>();
  private FileHandler fHandler;
  private static final Logger fLogger = Util.getLogger(LoggingConfigImpl.class);
  
  private void fetchSettings(ServletConfig aConfig){
    //InitParam objects are immutable
    fLoggingDir = aConfig.getInitParameter("LoggingDirectory").trim();
    fLoggingLevels = aConfig.getInitParameter("LoggingLevels").trim();
    logStdOut("Logging directory from web.xml : " + Util.quote(fLoggingDir));
    logStdOut("Logging levels from web.xml : " + Util.quote(fLoggingLevels));
  }
  
  private boolean isTurnedOff(){
    return NONE.equalsIgnoreCase(fLoggingDir.trim());
  }
  
  private void validateDirectorySetting() {
    if( ! fLoggingDir.endsWith(Consts.FILE_SEPARATOR) ){
      String message = "*** PROBLEM *** LoggingDirectory setting in web.xml does not end in with a directory separator : " + Util.quote(fLoggingDir);
      logStdOut(message);
      throw new IllegalArgumentException(message);
    }
    if( ! targetDirectoryExists() ){
      String message = "LoggingDirectory setting in web.xml does not refer to an existing, writable directory. Will attempt to create directory : " + Util.quote(fLoggingDir);
      logStdOut(message);
      File directory = new File(fLoggingDir);
      boolean success = directory.mkdirs();
      if (success) {
        logStdOut("Directory created successfully");
      }
      else {
        logStdOut("*** PROBLEM *** : Unable to create LoggingDirectory specified in web.xml! Permissions problem? Directory already exists, but not writable?");
      }
    }
  }
  
  private void parseLoggers(){
    StringTokenizer parser = new StringTokenizer(fLoggingLevels, ",");
    while ( parser.hasMoreElements() ){
      String rawItem = (String)parser.nextElement();
      int separator = rawItem.indexOf(SEPARATOR);
      String logger = rawItem.substring(0, separator).trim();
      String level = rawItem.substring(separator + 1).trim();
      addLogger(removeSuffix(logger), level);
    }
  }
  
  private String removeSuffix(String aLogger){
    int suffix = aLogger.indexOf(".level");
    if ( suffix == Consts.NOT_FOUND ) {
      throw new IllegalArgumentException("*** PROBLEM *** LoggingLevels setting in web.xml does not end with '.level'");
    }
    return aLogger.substring(0, suffix);
  }
  
  private void addLogger(String aLogger, String aLevel){
    if( ! Util.textHasContent(aLogger) ){
      throw new IllegalArgumentException("Logger name specified in web.xml has no content.");
    }
    Logger logger = Logger.getLogger(aLogger); //creates Logger if does not yet exist
    logger.setLevel(Level.parse(aLevel));
    fLogger.config("Adding Logger " + Util.quote(logger.getName() ) + " with level " + Util.quote(logger.getLevel()) );
    fLoggers.add(logger);
  }
  
  private void createFileHandler() {
    try {
      fHandler = new FileHandler(getFileName(), MAX_BYTES, NUM_FILES, APPEND_TO_EXISTING);
      fHandler.setLevel(Level.FINEST);
      fHandler.setFormatter(new SimpleFormatter());
    }
    catch (IOException ex){
      throw new RuntimeException("Cannot create FileHandler: " + ex.toString() , ex);
    }
  }
  
  private void attachLoggersToFileHandler(){
    for (Logger logger: fLoggers){
      if( hasNoFileHandler(logger) ){
        logger.addHandler(fHandler);
      }
    }
  }
  
  private boolean hasNoFileHandler(Logger aLogger){
    boolean result = true;
    Handler[] handlers = aLogger.getHandlers();
    fLogger.config("Logger " + aLogger.getName() + " has this many existing handlers: " + handlers.length);
    for (int idx = 0; idx < handlers.length; ++idx){
      if ( FileHandler.class.isAssignableFrom(handlers[idx].getClass()) ){
        fLogger.config("FileHandler already exists for Logger " + Util.quote(aLogger.getName()) + ". Will not add a new one.");
        result = false;
        break;
      }
    }
    return result;
  }
  
  /** Log a test message at each logger's configured level. */
  private void tryTestMessages(){
    logStdOut("Sending test messages to configured loggers. Please confirm output to above log file.");
    for(Logger logger: fLoggers){
      logger.log(logger.getLevel(), "This is a test message for Logger " + Util.quote(logger.getName()));
    }
  }

  /**
   Return the complete name of the logging file.
   Example file name : <tt>C:\log\fish_and_chips\2007_12_31_23_59.txt</tt>
  */
  private String getFileName(){
    String result = null;
    Date now = new Date(System.currentTimeMillis());
    SimpleDateFormat formatter = new SimpleDateFormat("yyyy_MM_dd_HH_mm"); //"YYYY|_|MM|_|DD|_|hh|_|mm"
    result = fLoggingDir + formatter.format(now);
    result = result + ".txt";
    logStdOut("Logging file name : " + Util.quote(result));
    return result;
  }
  
  private boolean targetDirectoryExists(){
    File directory = new File(fLoggingDir);
    return directory.exists() && directory.isDirectory() && directory.canWrite();
  }
  
  private void logStdOut(Object aObject){
    String message = String.valueOf(aObject);
    System.out.println(message);
  }
  
  private void showLoggerLevels() {
    for(Logger logger : fLoggers){
      fLogger.config("Logger " + logger.getName() + " has level " + logger.getLevel());
    }
  }
}
