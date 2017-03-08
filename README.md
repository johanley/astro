# astro
Full source for astronomytonight.net, a site for amateur astronomers. The site has basic astronomy data for any location, and weather data for specific countries (the US, Canada, and the UK). Customizable for your observing location.

# Running locally
The code is 95% javascript, and 5% Java.
The site is implemented as a .war file (servlets), running on Tomcat. 
Courtesy binaries of the Java classes are checked into source, for those who don't care about changing the Java code.

In order to run the app, you will need a Java Runtime Environment (JRE). 
In order to compile new versions of the Java code, you will need a Java Development Kit (JDK).

In prod, the version of Tomcat is version 6, running on JRE 1.6. (That's a bit stale.)

The layout of the source tree matches the runtime layout of the web site. 
Thus, running the app locally is simply a matter of pointing Tomcat to the root of the project. 
The WEB-INF directory is specific to servlets and Java code. 
The 'dev' directory contains tools (Java and/or javascript) that aren't needed at runtime, so the dev directory is removed in the prod depolyment.
