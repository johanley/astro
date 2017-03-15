package astro;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.logging.Logger;

import javax.servlet.RequestDispatcher;
import javax.servlet.ServletConfig;
import javax.servlet.ServletException;
import javax.servlet.http.HttpServlet;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

import astro.util.Util;

/**
 Hide app implementation details.
 
 <P>This servlet exists only for one reason: to hide from the client how a particular feature 
 is implemented (html, jsp, or servlet).
 
 <P>The benefit is that it lets the app evolve, with no need to add tedious redirects when, for example, an .html 
 impl changes into a .jsp, or a .jsp impl changes into a servlet+jsp impl. 
*/
public class HideImpl extends HttpServlet {
  
  @Override public void init(ServletConfig aConfig) throws ServletException {
    super.init(aConfig);
    fLogger.config("Starting servlet that hides impl details (html, jsp, or servlet).");
  }
  
  /** The 'extension' seen by the browser - {@value}. I may as well use this instead of '.do'. It may help search engines find the site. */
  public static final String EXT = ".sky";
  
  /** 
    Forward to an impl. The impl may be html, jsp, or code.
    <P>The number of URLs in this app is small, so manual mapping is not very painful.   
  */
  @Override protected void doGet(HttpServletRequest aRequest, HttpServletResponse aResponse) throws ServletException, IOException {
    if (isBadRequest(aRequest)){
      sendFailedResponse("Request too large.", aRequest, aResponse);
    }
    else {
      //map all /blah/blah.sky paths to a corresponding /blah/blah.jsp, under the app root
      replaceExtensionWith(".jsp", aRequest, aResponse);
      //or pass execution to some class, which will later forward/redirect as needed
    }
  }
  
  //PRIVATE
  private static final Logger fLogger = Util.getLogger(FetchFromThirdParty.class);
 
  private boolean isBadRequest(HttpServletRequest aRequest){
    return (aRequest.getContentLength() > 50*1024); 
  }
 
  private void sendFailedResponse(String aText, HttpServletRequest aRequest, HttpServletResponse aResponse){
    sendErrorResponse(aText, 500, aRequest, aResponse);
  }

  /** In the simplest case, just forward to the corresponding .html or .jsp in the app root */
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
  
  private void redirectTo(String url, HttpServletResponse response){
    fLogger.fine("Redirect : " + url);
    try {
      response.sendRedirect(url);
    } 
    catch (IOException e) {
      //this should never happen
      e.printStackTrace();
    }
  }
  
  /** Send data to the client as text. Does not use HTML. Headers are added here. */
  private void sendErrorResponse(String aText, int aStatus, HttpServletRequest aRequest, HttpServletResponse aResponse){
    try {
      aResponse.setStatus(aStatus);
      //fLogger.fine("Sending response in uncompressed form, as UTF-8. Length: " + aText.length());
      String utf8Text = new String(aText.getBytes(), "UTF-8");
      aResponse.setCharacterEncoding("UTF-8");
      aResponse.setContentLength(utf8Text.getBytes().length); 
      aResponse.setContentType("application/text");
      PrintWriter out = aResponse.getWriter();
      out.append(utf8Text);
    } 
    catch (IOException ex) {
      //in practice this won't happen - better handling???
      ex.printStackTrace();
    }
  }
}
