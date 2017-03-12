package astro.translate;

import java.io.IOException;
import java.util.logging.Logger;

import javax.servlet.Filter;
import javax.servlet.FilterChain;
import javax.servlet.FilterConfig;
import javax.servlet.ServletException;
import javax.servlet.ServletRequest;
import javax.servlet.ServletResponse;
import javax.servlet.http.Cookie;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;
import javax.servlet.http.HttpSession;

import astro.util.Util;

/**
 Detect the user's choice for their desired language for the UI. English is the default.
 
 <P>Examines two sources of info:
 <ol>
  <li>a request parameter named <tt>lang</tt>. 
  <b>The value must be the same as that accepted by the lang attribute of the &lt;html&gt; tag.</b>
  <li>a cookie named <tt>lang</tt>, whose value again matches that of the lang attribute of the &lt;html&gt; tag.  
 </ol>
 
 This class will create a session only if needed.
*/
public class LanguageDetector implements Filter {

  /** Name of a request param, or a cookie: {@value}.*/
  public static final String KEY_NAME = "lang";
  
  /** Name of the lang as stored in the session. This name will avoid possible collisions with other tools: {@value}.*/
  public static final String KEY_NAME_SESSION = "astro.lang";

  @Override public void doFilter(ServletRequest aRequest, ServletResponse aResponse, FilterChain aChain) throws IOException, ServletException {
    HttpServletRequest request = (HttpServletRequest) aRequest;
    HttpServletResponse response = (HttpServletResponse) aResponse;
    fLogger.fine("Detect language filter.");
    detectAndUpdateLanguageChoice(request, response); //executed *before* regular servlet processing
    aChain.doFilter(aRequest, aResponse);
    //code here is executed *after* regular servlet processing
  }
  
  /** This impl does nothing. */
  @Override public void init(FilterConfig aConfig) throws ServletException {}
  
  /** This impl does nothing. */
  @Override public void destroy() {}
  
  // PRIVATE 
  
  private static final Logger fLogger = Util.getLogger(LanguageDetector.class);
  private static final Boolean CREATE_SESSION_IF_ABSENT = Boolean.TRUE;
  private static final Boolean DONT_CREATE_SESSION_IF_ABSENT = Boolean.FALSE;
  
  private void detectAndUpdateLanguageChoice(HttpServletRequest request, HttpServletResponse response){
    String current = currentLangChoice(request);
    fLogger.fine("Current lang choice in session: " + current);
    String incoming = incomingLangChoice(request);
    fLogger.fine("Incoming lang choice: " + incoming);
    if (hasChanged(incoming, current)){
      fLogger.fine("Updating the user's lang choice from: " + current + " to: "+ incoming);
      updateTo(incoming, request, response);
    }
  }
  
  /** Possibly null. Will not create a new session. */
  private String currentLangChoice(HttpServletRequest request){
    String result = null;
    HttpSession session = request.getSession(DONT_CREATE_SESSION_IF_ABSENT);
    if (session != null){
      result = (String)session.getAttribute(KEY_NAME_SESSION);
    }
    return result;
  }
  
  /** Possibly null. */
  private String incomingLangChoice(HttpServletRequest request){
    String value = null;
    value = request.getParameter(KEY_NAME);
    if (value == null){
      if (request.getCookies() != null){
        for (Cookie cookie : request.getCookies()){
          if (cookie.getName().equals(KEY_NAME)){
            value = cookie.getValue();
          }
        }
      }
    }
    return value;
  }
  
  private boolean hasChanged(String incoming, String current){
    boolean result = false;
    if (incoming == null && current != null){
      //do nothing!
    }
    else if (incoming != null && current == null){
      result = true;
    }
    else if (incoming != null && current != null){
      result = !incoming.equalsIgnoreCase(current);
    }
    return result;
  }
  
  private void updateTo(String incoming, HttpServletRequest request, HttpServletResponse response){
    HttpSession session = request.getSession(CREATE_SESSION_IF_ABSENT);
    session.setAttribute(KEY_NAME_SESSION, incoming);
    Cookie userCookie = new Cookie(KEY_NAME, incoming);
    userCookie.setMaxAge(60*60*24*365); //1 year
    response.addCookie(userCookie);    
  }
}
