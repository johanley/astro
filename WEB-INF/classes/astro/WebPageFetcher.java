package astro;

import java.io.IOException;
import java.net.URL;
import java.net.URLConnection;
import java.util.Scanner;
import java.util.logging.Logger;

import astro.util.Util;

/** 
 Fetch the text of a web page at a given URL. Return it as a simple String.
 
 <P>Has a 60sec timeout for both getting the connection, and for reading the content.
 <p>About 4-5 GETs can be done per second.
*/
public final class WebPageFetcher {

  /** Returns an empty string if a problem occurs. */
  public String fetch(String aURL, String aEncoding) {
    return webFetch(aURL, aEncoding); 
  }
  
  // PRIVATE
  private static final int TIMEOUT = 60*1000; //msecs
  private static final String END_OF_INPUT = "\\Z";
  private static final Logger fLogger = Util.getLogger(FetchFromThirdParty.class);
  
  private String webFetch(String aURL, String aEncoding) {
    //fLogger.fine("Fetching from the web.");
    String result = "";
    long start = System.currentTimeMillis();
    Scanner scanner = null;
    URLConnection connection = null;
    try {
      //Ref: http://stackoverflow.com/questions/86824/why-would-a-java-net-connectexception-connection-timed-out-exception-occur-wh
      URL url = new URL(aURL);
      connection =  url.openConnection(); //this doesn't talk to the network yet
      connection.setConnectTimeout(TIMEOUT); 
      connection.setReadTimeout(TIMEOUT);
      connection.connect(); //actually connects; this shouldn't be needed here, since getInputStream is supposed to do it in the background
      scanner = new Scanner(connection.getInputStream(), aEncoding); 
      scanner.useDelimiter(END_OF_INPUT);
      result = scanner.next();
    }
    catch (IOException ex) {
      long end = System.currentTimeMillis();
      long time = end - start;
      fLogger.severe(
        "Problem connecting to " + aURL + " Encoding:" + aEncoding + 
        ". Exception: " + ex.getMessage() + " " + ex.toString() + " Cause:" +  ex.getCause() + 
        " Connection Timeout: " + connection.getConnectTimeout() + "msecs. Read timeout:" + connection.getReadTimeout() + "msecs."
        + " Time taken to fail: " + time + " msecs."
      );
    }
    finally {
      if (scanner != null) scanner.close();
    }
    return result;
  }
 
}