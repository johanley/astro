package astro.translate;

import java.io.*;
import java.util.logging.*;

import javax.servlet.Servlet;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;
import javax.servlet.jsp.JspException;
import javax.servlet.jsp.tagext.SimpleTagSupport;
import javax.servlet.jsp.PageContext;
import javax.servlet.jsp.tagext.JspFragment;
import javax.servlet.jsp.JspContext;

import astro.util.Util;

/**
 Base class for implementing custom JSP tags.
 This class was taken from the web4j project.
 
 <P>The custom tag can optionally have a body. The <tt>.tld</tt> entry for these tags must 
 have their <tt>body-content</tt> set to <tt>scriptless</tt> (the default).
 
 <P>Concrete subclasses of this class perform these tasks :
<ul>
 <li>implement <tt>setXXX</tt> methods, one for each tag attribute;
 each <tt>setXXX</tt> should validate its argument.
 <li>optionally override the {@link #crossCheckAttributes} method, to 
 perform validations depending on more than one attribute.
 <li>implement {@link #getEmittedText}, to return the text to be included in markup.
</ul>
*/
public abstract class TagHelper extends SimpleTagSupport {
  
  /**
   <b>Template</b> method which calls {@link #getEmittedText(String)}.
   
   <P>The body of this tag is evaluated, passed to {@link #getEmittedText(String)}, 
   and the result is then written to the JSP output writer. In addition, this method will call
   {@link #crossCheckAttributes()} at the start of processing.
  */
  @Override public final void doTag() throws JspException {
    try {
      crossCheckAttributes();
      getJspContext().getOut().write(getEmittedText(getBody()));
    }
    catch (Throwable ex){
      fLogger.severe("Cannot execute custom tag. " + Util.quote(ex));
      throw new JspException("Cannot execute custom tag.", ex);
    }
  }
  
  /** 
   Return the text this tag will display in the resulting web page.
   
   @param aOriginalBody is the evaluated body of this tag. If there is no body, or 
   if the body is present but empty, then it is <tt>null</tt>.
   @return the text to display in the resulting web page. 
  */
  abstract protected String getEmittedText(String aOriginalBody) throws JspException, IOException;
   
  /**
   Perform validations that apply to more than one attribute.
  
   <P>This default implementation does nothing.
   
   <P>Validations that apply to a single attribute should be performed in its 
   corresponding <tt>setXXX</tt> method.
  
   <P>If a problem is detected, subclasses must emit a <tt>RuntimeException</tt> 
   describing the problem. If all validations apply to only to a single attribute, 
   then this method should not be overridden.
  */
  protected void crossCheckAttributes() {
    //do nothing in this default impl
  }
  
  /** Return the underlying {@link HttpServletRequest}.  */
  protected final HttpServletRequest getRequest(){
    return (HttpServletRequest)getPageContext().getRequest();
  }
  
  /** Return the underlying {@link HttpServletResponse}.  */
  protected final HttpServletResponse getResponse(){
    return (HttpServletResponse)getPageContext().getResponse();
  }
  
  /** Return the underlying {@link PageContext}.  */
  protected final PageContext getPageContext(){
    JspContext jspContext = getJspContext();
    return (PageContext)jspContext;
  }

  /** 
   Return the name of the JSP implementation class. 
   <P>Intended for debugging only. 
  */
  protected final String getPageName(){
    Servlet servlet = (Servlet)getPageContext().getPage();
    return servlet.getClass().getName();
  }
  
  /**
   Verify that an attribute value has content. 
   
   <P>If no content, then log at <tt>SEVERE</tt> and throw an unchecked exception.
  */
  protected final void checkForContent(String aAttributeName, String aAttributeValue){
    if( ! Util.textHasContent(aAttributeValue) ){
      String message = Util.quote(aAttributeName) + " attribute must have a value.";
      fLogger.severe(message);
      throw new IllegalArgumentException(message);
    }
  }

  // PRIVATE //
  
  private static final Logger fLogger = Util.getLogger(TagHelper.class);
  
  /** 
   Return the evaluated body of this tag.  
   
   <P>The body of this tag cannot contain scriptlets or scriptlet expressions.
   If this tag has no body, or has an empty body, then <tt>null</tt> is returned.
  */
  private String getBody() throws IOException, JspException {
    String result = null;
    JspFragment body = getJspBody();
    if( body != null ){
      StringWriter writer = new StringWriter();
      getJspBody().invoke(writer);
      writer.flush();
      result = writer.toString();
    }
    return result;
  }
}
