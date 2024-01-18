package astro;

import static planisphere.config.Constants.STAR_CHART_FILE;
import static planisphere.config.Constants.TRANSPARENCY_FILE;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.io.Writer;
import java.util.logging.Logger;

import javax.servlet.RequestDispatcher;
import javax.servlet.ServletConfig;
import javax.servlet.ServletException;
import javax.servlet.http.HttpServlet;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

import planisphere.config.Config;
import planisphere.draw.starchart.GenerateStarChart;
import planisphere.draw.transparency.GenerateTransparency;
import planisphere.math.Maths;
import planisphere.util.LogUtil;
import astro.util.Util;

/** Serve files used to make a planisphere. */
public class Planisphere extends HttpServlet {
  
  @Override public void init(ServletConfig aConfig) throws ServletException {
    super.init(aConfig);
    fLogger.config("Starting servlet that generates PDFs for a planisphere.");
  }

  /** The 'extension' seen by the browser - {@value}. I may as well use this instead of '.do'. It may help search engines find the site. */
  public static final String EXT = ".pln";
  public static final String ENCODING = "UTF-8";  

  /**
   Serve a PDF file used to make a planisphere.
   If user input is invalid, redisplay the form. 
  */
  @Override protected void doGet(HttpServletRequest request, HttpServletResponse response) throws ServletException, IOException {

    if (isBadRequest(request)){
      sendFailedResponse("Request too large.", request, response);
    }
    else {
      Integer year = asInt("year", request);
      String location = request.getParameter("location");
      Double latitude = asRads("latitude", request);
      Double longitude = asRads("longitude", request);
      Integer hoursOffsetFromUT = asInt("hoursOffsetFromUT", request);
      Double declinationGap = asDouble("declinationGap", request);
      Integer greyScaleAltAz = asInt("greyScaleAltAz", request);
      if (anythingIsNull(year, location, latitude, longitude, hoursOffsetFromUT, declinationGap, greyScaleAltAz)){
        request.setAttribute("error", "Error occurred. Please try again.");
        showTheForm(request, response); //goes back to the default form
      }
      else {
        //hard-coded config for many items
        int minutesOffsetFromUT = 0;
        float width = points(8.5);
        float height = points(11.0);
        String outputDir = "";
        String fontDir = getServletContext().getInitParameter("fontDirectory");
        Integer greyConstellationLines = 128;
        Integer smallestTimeDivision = 2;
        String radiants = "perseids:46.2,57.4 | eta-aquarids:338.0,-1.0 | quadrantids:230.1,48.5 | geminids:112.3,32.5";
        String monthNames = "Jan, Feb, Mar, Apr, May, Jun, Jul, Aug, Sep, Oct, Nov, Dec";
        String lunar_transits_title = "Daily Lunar Transits";
        String planetary_transits_title = "Planetary Transits on the 15th of the Month";
        String planet_names = "Mercury, Venus, Earth, Mars, Jupiter, Saturn";
        Boolean discardPolaris = true;
        
        Config config = new Config(
          year, location, latitude, longitude, hoursOffsetFromUT, minutesOffsetFromUT, declinationGap, width, height, 
          outputDir, fontDir,  greyConstellationLines, greyScaleAltAz, smallestTimeDivision, radiants, monthNames, 
          lunar_transits_title, planetary_transits_title, planet_names, discardPolaris
        );
        fLogger.info("Planisphere lat:" + Maths.radsToDegs(config.latitude()) + " long:" + Maths.radsToDegs(config.longitude()));
        response.setContentType("application/pdf");
        response.setHeader("Expires", "0");
        
        ByteArrayOutputStream baos = new ByteArrayOutputStream();
        String action = request.getParameter("generate");
        String fileName = "";
        try {
          if ("Star chart".equalsIgnoreCase(action)){
            LogUtil.log("Star chart.");
            fileName = STAR_CHART_FILE;
            GenerateStarChart starChart = new GenerateStarChart(config);
            starChart.outputTo(baos);
          }
          else {
            LogUtil.log("Transparency.");
            fileName = TRANSPARENCY_FILE;
            GenerateTransparency transparency = new GenerateTransparency(config);
            transparency.outputTo(baos);
          }
        }
        catch (Throwable e) {
          sendFailedResponse("Error. Not able to generate PDF", request, response);
        }
        
        response.setHeader("Content-disposition","attachment;filename="+ fileName);
        response.setContentLength(baos.size());
        // write baos to the ServletOutputStream
        OutputStream os = response.getOutputStream();
        baos.writeTo(os);
        os.flush();
        os.close();
      }
    }
  }
  
  //PRIVATE
  
  private static final Logger fLogger = Util.getLogger(Planisphere.class);
  private static final int POINTS_PER_INCH = 72;
  
  private float points(double inch){
    return (float)(POINTS_PER_INCH * inch);
  }
  
  private boolean isBadRequest(HttpServletRequest aRequest){
    return (aRequest.getContentLength() > 50*1024); 
  }
  
  private Integer asInt(String name, HttpServletRequest req){
    Integer result = null;
    String raw = req.getParameter(name);
    try {
      result = Integer.valueOf(raw);
    }
    catch(NumberFormatException ex){
      //nothing
    }
    return result;
  }
  
  private Double asRads(String name, HttpServletRequest req){
    Double result = null;
    Double degs = asDouble(name, req);
    if (degs != null){
      result = Maths.degToRads(degs);
    }
    return result;
  }
      
  private Double asDouble(String name, HttpServletRequest req){
    Double result = null;
    String raw = req.getParameter(name);
    try {
      result = Double.valueOf(raw);
    }
    catch(NumberFormatException ex){
      //nothing
    }
    return result;
  }
  
  private boolean anythingIsNull(Object... args){
    boolean result = false;
    for(Object arg : args){
      if (arg == null){
        result = true;
        break;
      }
    }
    return result;
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

  private void logProblem(String msg, Throwable ex){
    fLogger.severe(msg);
    fLogger.severe("Exception message: " + ex.getMessage() + " toString: " + ex.toString());
    fLogger.severe("Java version: " + System.getProperty("java.version"));
    fLogger.severe(getStackTrace(ex));
  }
  
  private String getStackTrace(Throwable aThrowable) {
    Writer result = new StringWriter();
    PrintWriter printWriter = new PrintWriter(result);
    aThrowable.printStackTrace(printWriter);
    return result.toString();
  }
  
  private void showTheForm( HttpServletRequest aRequest, HttpServletResponse aResponse){
    replaceExtensionWith(".jsp", aRequest, aResponse);
  }
  
  private void replaceExtensionWith(String replacement, HttpServletRequest aRequest, HttpServletResponse aResponse){
    String path = aRequest.getServletPath(); // eg /main/form.sky
    Integer ext = aRequest.getServletPath().indexOf(EXT); //location of .sky
    String urlStart = aRequest.getServletPath().substring(0, ext); //part before .sky
    if (path.endsWith(EXT)){
      path = urlStart + replacement; //.sky to .jsp 
    }
    else {
      String urlEnd = aRequest.getServletPath().substring(ext + EXT.length());
      path = urlStart + replacement + urlEnd; //change .sky to .jsp, and then append any items after .sky
    }
    forwardTo(path, aRequest, aResponse);
  }
  
  private void forwardTo(String destination, HttpServletRequest aRequest, HttpServletResponse aResponse){
    RequestDispatcher dispatcher = aRequest.getRequestDispatcher(destination);
    aResponse.setCharacterEncoding("UTF-8");
    try {
      dispatcher.forward(aRequest, aResponse);
    } 
    catch (Throwable e) {
      e.printStackTrace();
    }
  }
}