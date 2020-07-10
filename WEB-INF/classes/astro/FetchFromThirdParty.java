package astro;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.io.Writer;
import java.net.URL;
import java.util.UUID;
import java.util.logging.Logger;

//import javax.imageio.ImageIO; //this guy is problematic! initialization issues!
import javax.servlet.ServletConfig;
import javax.servlet.ServletException;
import javax.servlet.ServletOutputStream;
import javax.servlet.http.HttpServlet;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

import astro.util.LoggingConfigImpl;
import astro.util.Util;

/**
<P>This server exists for one reason: getting around CORS (cross-origin) issues in javascript.
The client is not, by default, allowed to access content from third-party servers.
The workaround here is to fetch content on the client's behalf.

<P>Currently, these types of files are supported:
<ul> 
 <li>.png (the default)
 <li>.xml (rss feeds)
 <li>.txt (plain text)
</ul>
Adding support for other file types, if ever needed, should be easy.

<P>Example URLs:
<pre>
localhost:8081/playtime/fetch/?url=http://bazoodi.com/&ext=png
localhost:8081/playtime/fetch/?url=http://weather.gc.ca/rss/city/on-118_e.xml&ext=xml
localhost:8081/playtime/fetch/?radarStation=XFT
playtime.ca/fetch/?url=http://weather.gc.ca/rss/city/on-118_e.xml&ext=xml
</pre>
*/
public final class FetchFromThirdParty  extends HttpServlet {
  
  public static final String ENCODING = "UTF-8";  
  
  /** Initialize logging. */
  @Override public void init(ServletConfig aConfig) throws ServletException {
    super.init(aConfig);
    LoggingConfigImpl loggingConfig = new LoggingConfigImpl();
    loggingConfig.setup(aConfig);
    fLogger.config("Starting fetcher servlet.");
    fLogger.config("Location of temp files: " + getServletContext().getAttribute("javax.servlet.context.tempdir"));
  }
  
  /** 
   Return a resource (image or text) from the web, from a third-party's server.
   Only GET operations are supported.
  */
  @Override protected void doGet(HttpServletRequest aRequest, HttpServletResponse aResponse) throws ServletException, IOException {
    if (isBadRequest(aRequest)){
      sendFailedResponse("Request too large.", aRequest, aResponse);
    }
    else {
      //String url = aRequest.getQueryString().substring("url=".length());
      String url = aRequest.getParameter("url");
      String radarStation = aRequest.getParameter("radarStation");
      if (Util.textHasContent(url)){
        String extension = aRequest.getParameter("ext");
        if (extension == null || extension.equals("")){
          extension = PNG;
        }
        fetchFromThirdParty(url, extension, aRequest, aResponse);
      }
      else if (Util.textHasContent(radarStation)){
        returnUrlOfMostRecentRadarImage(radarStation, aResponse);
      }
    }
  }
  
  //PRIVATE
  private static final Logger fLogger = Util.getLogger(FetchFromThirdParty.class);
  private static final String PNG = "png";
  private static final String XML = "xml";
  private static final String TXT = "txt";
  private static final Integer NUM_RADAR_IMAGES = 3;
 
  private void fetchFromThirdParty(String aURL, String extension, HttpServletRequest aRequest, HttpServletResponse aResponse) throws ServletException, IOException {
    fLogger.info("Fetching URL: " + aURL);
    File file = buildOutputFile(extension);
    //fLogger.fine("Temp file name: " + file.getName());
    String mimeType = "";
    if (PNG.equals(extension)){
      //fLogger.fine("Fetching image from web.");
      
      //BufferedImage image = fetchImage(aURL);
      //saveAsTempFile(image, file);
      
      saveImageAsTempFile(aURL, file);
      
      //fLogger.fine("Serving as image file.");
      mimeType = "image/png";
    }
    else if (XML.equals(extension)){
      String text = fetchText(aURL, ENCODING);
      saveTextAsTempFile(text, file, ENCODING);
      //fLogger.fine("Serving as a text file."); 
      mimeType = "text/xml";
    }
    else if (TXT.equals(extension)){
      String text = fetchText(aURL, ENCODING);
      saveTextAsTempFile(text, file, ENCODING);
      //fLogger.fine("Serving as a text file."); 
      mimeType = "text/plain";
    }
    serveFile(file, aRequest, aResponse, mimeType);
    //fLogger.fine("Deleting temporary file."); 
    deleteTempFile(file);
  }
  
  /** Return null if a problem occurs. */
  /* Image IO is problematic: initialization errors in prod.
  private BufferedImage fetchImage(String aURL) {
    BufferedImage result = null;
    try {
      //fLogger.info("Formats known to the currently registered readers: " + Arrays.asList(ImageIO.getReaderFormatNames()));
      //fLogger.info("Trying to fetch image. URL: " + aURL);
      URL url = new URL(aURL);
      result = ImageIO.read(url);
      //fLogger.info("Result: "  + result);
    } 
    catch (IOException ex) {
      logProblem("Can't fetch image. URL: " + aURL, ex);
    }
    return result;
  }
  */
  
  /** Return null if a problem occurs. */
  private String fetchText(String aURL, String aEncoding) {
    WebPageFetcher fetcher = new WebPageFetcher();
    return fetcher.fetch(aURL, aEncoding);
  }
  
  /** 
   Temp output file in the servlet's temp dir. 
   Build a nonce name, unique to the host.
  */ 
  private File buildOutputFile(String extension){
    String name = "";
    UUID nonce = UUID.randomUUID(); //eg '067e6162-3b6f-4ae2-a171-2470b63dff00'
    name = name + nonce.toString();
    name = name + "." + extension;
    File tempDir = (File)getServletContext().getAttribute("javax.servlet.context.tempdir"); //CASE-SENSITIVE!!
    File result = new File(tempDir, name);
    return result;
  }

  /** Uses the built-in temp directory defined by the servlet spec. */
  /* Image IO is problematic: initialization errors in prod.
  private void saveAsTempFile(BufferedImage image, File outputfile){
    try {
      ImageIO.write(image, "png", outputfile);
    } 
    catch (IOException ex) {
      logFileProblem("Unable to save image bytes to a temp file.", outputfile, ex);
    }
    
    if (!outputfile.exists()){
      fLogger.severe("Unable to write temp image file.");
    }
    else if (outputfile.length() == 0) {
      fLogger.severe("Temp image file has 0 length.");
    }
  }
  */
  
  /** Image IO is problematic: initialization errors in prod. So, we avoid using it here. */
  private void saveImageAsTempFile(String imageUrl, File outputFile) {
    InputStream is = null;
    OutputStream os = null;
    try {
      URL url = new URL(imageUrl);
      try {
        is = url.openStream();
        os = new FileOutputStream(outputFile);
        byte[] b = new byte[2048];
        int length;
        while ((length = is.read(b)) != -1) {
          os.write(b, 0, length);
        }
        fLogger.info("Saved image for " + imageUrl);
      }
      finally {
        is.close();
        os.close();
      }
    }
    catch (Throwable ex){
      logProblem("Unable to save image file for " + imageUrl, ex);
    }
  }  
  
  /** Uses the built-in temp directory defined by the servlet spec. */
  private void saveTextAsTempFile(String text, File outputFile, String encoding){
    if (text.length() < 1){
      fLogger.severe("Text fetched from web has 0 length.");
    }
    Writer out = null;
    try {
      try {
        out = new OutputStreamWriter(new FileOutputStream(outputFile), encoding);
        out.write(text);
      }
      finally {
        if (out != null) out.close();
      }    
    }
    catch (IOException ex) {
      logFileProblem("Unable to save text file bytes to a temp file.", outputFile, ex);
    }
    
    if (!outputFile.exists()){
      fLogger.severe("Unable to write temp text file.");
    }
    else if (outputFile.length() == 0) {
      fLogger.severe("Temp text file has 0 length.");
    }
    else {
      //fLogger.info("Saved temp file: " + outputFile.getAbsolutePath());
      //fLogger.info("Size of saved temp file: " + outputFile.length());
    }
  }

  private void serveFile(File file, HttpServletRequest aRequest, HttpServletResponse aResponse, String mimeType) throws IOException {
    aResponse.setContentType(mimeType);
    InputStream input = null;
    try {
      ServletOutputStream output = aResponse.getOutputStream();
      input = new FileInputStream(file);
      //transfer input stream to output stream, via a buffer
      byte[] buffer = new byte[2048];
      int bytesRead;    
      while ((bytesRead = input.read(buffer)) != -1) {
         output.write(buffer, 0, bytesRead);
      }
      //fLogger.info("Served saved file: " + file.getAbsolutePath());
    }
    catch(IOException ex){
      logFileProblem("Unable to serve temp file.", file, ex);
    }
    finally {
      if (input != null) input.close();
    }
  }
  
  private void deleteTempFile(File file){
    file.delete();
  }
  
  private boolean isBadRequest(HttpServletRequest aRequest){
    return (aRequest.getContentLength() > 50*1024); 
  }
 
  private void sendFailedResponse(String aText, HttpServletRequest aRequest, HttpServletResponse aResponse){
    sendErrorResponse(aText, 500, aRequest, aResponse);
  }
  
  /** Send data to the client as text. Does not use HTML. Headers are added here. */
  private void sendErrorResponse(String aText, int aStatus, HttpServletRequest aRequest, HttpServletResponse aResponse){
    try {
      aResponse.setStatus(aStatus);
      //fLogger.fine("Sending response in uncompressed form, as UTF-8. Length: " + aText.length());
      String utf8Text = new String(aText.getBytes(), ENCODING);
      aResponse.setCharacterEncoding("UTF-8");
      aResponse.setContentLength(utf8Text.getBytes().length); 
      aResponse.setContentType("application/text");
      PrintWriter out = aResponse.getWriter();
      out.append(utf8Text);
    } 
    catch (IOException ex) {
      //in practice this won't happen
      logProblem("Problem sending an error response.", ex);
    }
  }

  private void returnUrlOfMostRecentRadarImage(String radarStation, HttpServletResponse response){
    Radar radar = new Radar(NUM_RADAR_IMAGES);
    String json = radar.returnJsonMostRecentRadarImages(radarStation);
    serveAsJson(json, response);
  }
  
  private void serveAsJson(String jsonText, HttpServletResponse response){
    //fLogger.fine("Sending json response: " + jsonText);
    response.setContentType("application/json");
    try {
      response.getOutputStream().print(jsonText);
    }
    catch(IOException ex){
      fLogger.severe("Unable to serve json: " + jsonText);
    }
  }
  
  private void logProblem(String msg, Throwable ex){
    fLogger.severe(msg);
    fLogger.severe("Exception message: " + ex.getMessage() + " toString: " + ex.toString());
    fLogger.severe("Java version: " + System.getProperty("java.version"));
    fLogger.severe(getStackTrace(ex));
  }
  
  private void logFileProblem(String msg, File file, Throwable ex){
    logProblem(msg, ex);
    fLogger.severe("File name: " + file);
    fLogger.severe("Absolute path: " + file.getAbsolutePath());
    fLogger.severe("File exists: " + file.exists());
    fLogger.severe("Length: " + file.length());
  }
  
  private String getStackTrace(Throwable aThrowable) {
    Writer result = new StringWriter();
    PrintWriter printWriter = new PrintWriter(result);
    aThrowable.printStackTrace(printWriter);
    return result.toString();
  }
}
 