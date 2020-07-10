package astro;

import java.io.IOException;
import java.net.URL;
import java.net.URLConnection;
import java.util.Scanner;
import java.util.logging.Logger;

import javax.net.ssl.HttpsURLConnection;

import astro.util.Util;

/** 
 Fetch the text of a web page at a given URL. Return it as a simple String.
 
 <P>Has a 60sec timeout for both getting the connection, and for reading the content.
 <p>About 4-5 GETs can be done per second.
*/
public final class WebPageFetcher {
  
  /** Simple test harness. */ 
  public static void main(String[] args) {
    WebPageFetcher fetcher = new WebPageFetcher();
    /*
     * -Dhttps.protocols=TLSv1.1,TLSv1.2 
     * https://stackoverflow.com/questions/21245796/javax-net-ssl-sslhandshakeexception-remote-host-closed-connection-during-handsh
    */
    System.setProperty("https.protocols", "TLSv1,TLSv1.1,TLSv1.2");
    //String text = fetcher.fetch("http://www.web4j.com", "UTF-8");
    String text = fetcher.fetch("https://weather.gc.ca/rss/city/on-118_e.xml", "UTF-8"); //fails: Remote host closed connection during handshake
    //String text = fetcher.fetch("https://www.google.com/", "UTF-8"); //succeeds
    System.out.println(text);
  }

  /** Returns an empty string if a problem occurs. */
  public String fetch(String aURL, String aEncoding) {
    return webFetch2(aURL, aEncoding); 
  }
  
  // PRIVATE
  private static final int TIMEOUT = 60*1000; //msecs
  private static final String END_OF_INPUT = "\\Z";
  private static final Logger fLogger = Util.getLogger(FetchFromThirdParty.class);
  
  private String webFetch1(String aURL, String aEncoding) {
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
  
  private String webFetch2(String targetUrl, String encoding){
    String result = "";
    URL url = null;
    URLConnection connection = null;
    /*
     * -Dhttps.protocols=TLSv1.1,TLSv1.2 
     * https://stackoverflow.com/questions/21245796/javax-net-ssl-sslhandshakeexception-remote-host-closed-connection-during-handsh
     * It's more elegant to do this only once; but it does no harm to do it multiple times, and I 
     * like having it here instead of some config somewhere. This may be specific to Java 7.
     *
     * OpenJDK8:
     * java.lang.RuntimeException: Unexpected error: java.security.InvalidAlgorithmParameterException: the trustAnchors parameter must be non-empty
     * https://stackoverflow.com/questions/4764611/java-security-invalidalgorithmparameterexception-the-trustanchors-parameter-mus
     * This happens for OpenJDK8. Apparently you need to copy in a cacerts file from a working jre; the default one is empty.
     * 
    */
    System.setProperty("https.protocols", "TLSv1,TLSv1.1,TLSv1.2");
    try {
      url = new URL(targetUrl);
      connection =  url.openConnection();
      Scanner scanner = null;
      try {
        scanner = new Scanner(connection.getInputStream());
        scanner.useDelimiter(END_OF_INPUT);
        result = scanner.next();
      }
      finally {
        if (scanner != null) scanner.close();
      }
    }
    catch (IOException ex) {
      fLogger.severe("Cannot open connection to " + url + " " + ex.getMessage());
    }
    return result;
  }
}