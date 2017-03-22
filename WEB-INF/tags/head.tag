<%-- Repeated markup that appears in the HEAD tag. --%>
<%@ include file="/WEB-INF/TagHeader.jspf" %>
<%@ tag pageEncoding="UTF-8" trimDirectiveWhitespaces="true" %>
 <%-- Note how the app name and version is config data, not hardcoded. --%>
 <title>Name: ${initParam.appName}</title> 
  
 <c:url var='stylesheet' value='/css/stylesheet001.css'/>
 <link rel="stylesheet" type="text/css" href="<c:out value='${stylesheet}'/>" media="all">
 
 <c:url var='icon' value='/images/favicon.ico'/>
 <link rel="shortcut icon" type="image/vnd.microsoft.icon" href="<c:out value='${icon}'/>">