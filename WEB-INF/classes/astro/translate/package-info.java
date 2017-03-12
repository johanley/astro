/** 
 Translate UI text from English to some other language.
 
 <P>Parts of the implementation:
 <ul>
  <li>a text file to hold all translations
  <li>a ServletFilter to detect the user's choice for the UI language
  <li>a custom JSP tag to do the translation
 </ul>
  
  The JSP contains English text, wrapped in a translation tag.
  (Using plain English text, instead of a coder-key, reads much better.)
  
  <P>The English text in the tag body acts as a key, and must match an entry 
  in the translations file. So, if the English text in the JSP changes, then 
  the translations file must also change as well, to remain in sync. 
  
  <P>Since getting out of sync is a problem, the system has a mechanism in place 
  to find such errors reasonably quickly.
  At runtime, if no translation is found, then the error is logged to an 
  error file. 
  
  <P>Thus, to find problems, all you have to do is run the app, 
  exercise all screens, and then check the errors file. 
  Since the app isn't very large, this doesn't take long.
  
  <P>Careful: the text is the translations file sometimes needs character 
  entities, to avoid conflicts with HTML. (Apostrophes, for example.)
*/
package astro.translate;
