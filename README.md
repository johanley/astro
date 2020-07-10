# astro
Full source for [astronomytonight.net](http://astronomytonight.net/main/form.sky), a site for amateur astronomers. The site has basic astronomy data for any location, and weather data for specific countries (the US and Canada). Customizable for your observing location.

The layout of the source tree matches the runtime layout of the web site.  
* running the app locally is simply a matter of pointing Tomcat to the root of the project.  
* Java classes are compiled-in-place, in the same directory as the source .java file 
* the `WEB-INF` directory is specific to servlets and Java code. 
* the 'dev' directory contains tools (Java and/or javascript) that aren't needed at runtime, so the dev directory is removed in the prod depolyment.

The code is roughly 95% javascript, and 5% Java.
No relational database is used.
The site is implemented as a .war file (servlets), running on [Tomcat](http://tomcat.apache.org/whichversion.html). 
Courtesy binaries of the Java classes are checked into source, for those who don't care about changing the Java code.

# Running locally

In order to run the app, you will need a Java Runtime Environment (JRE), and Tomcat (or a similar servlet container). 
In order to compile new versions of the Java code, you will need a Java Development Kit (JDK).

In prod, the version of Tomcat is version 6, running on JRE 1.6. (That's a bit stale.)

To run the app locally on Tomcat:
* download the code from github, to any location on your machine
* make sure you have a Java Runtime Environment installed (JRE 1.6 or better)
* download and install [Tomcat](http://tomcat.apache.org/whichversion.html), version 6 or higher
* a plain, default install of Tomcat will have its server running on localhost:8080/

After you've installed Tomcat, the last step is to point Tomcat to the root of the project tree for astronomytonight.net 
that you downloaded from github. You do that by creating a file in this location
* `{tomcat-home}/conf/Catalina/localhost/astro.xml`

This `astro.xml` file (called a context file) simply contains a pointer to the root of the project tree. An example (which you will need to modify, of course):
* `<Context docBase="C:\myname\projects\astro\" reloadable="true"></Context>`

Finally, there is a silly change you need to make to `js\util.js`. 
There's a reference to port 8081, which you will likely need to change to port 8080 (Tomcat's default port), or to 
whatever port you're running Tomcat on.   

The name of the `astro.xml` file controls the URL under which you'll see the app running:  
* `localhost:8080/astro/`


Other items to note:
* all files should use UTF-8 encoding.
* no tabs please.
* number of unit tests: 0.
* javascript is used in a functional style.
* the implementation uses [package-by-feature](http://www.javapractices.com/topic/TopicAction.do?Id=205).
* js/ephem.js contains the core of the astronomy algorithms. It's large and non-modular.
* no javascript modularization tool or build tool is used.
* preferred js naming style: `like_this`, but it doesn't really matter that much.
* note the use of Greek letters in the code: it reads much better, and lets you compare more easily with other sources.
* periodic manual updates to data are needed to keep the data current, but the effort required is not onerous (20 minutes monthly, 3 hours yearly).
