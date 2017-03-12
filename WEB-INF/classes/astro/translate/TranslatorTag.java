package astro.translate;

import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.util.Date;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Scanner;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.servlet.http.HttpSession;
import javax.servlet.jsp.JspException;

import astro.util.Util;

/**
 Translate text in a JSP. 
 The body of this tag is static English markup.
 This class just substitutes translated markup for the given English markup.
 Plain old text replacement. 
 The text can contain HTML tags, if needed.
 
 <P>This impl can handle new languages, with a 1-line change to the code.
*/
public class TranslatorTag extends TagHelper {

  /**
   The file that holds translations from English to N target languages: {@value}.
  
   The format of data in the file is <tt>{Saturn}{Saturne}</tt>. 
   English first, French second, and so on.
   The braces are the separators. 
   Follows an ordering convention, without naming the language.  
   The units must come in tuples; if they get out of the conventional order, things will get completely messed up! 
   But if that happens, the problem should be obvious.   
   Multiline data is allowed.
  */
  public static final String TRANSLATIONS_FILE_NAME = "translations.txt";
  
  /**
   Key for a hard-coded file path:  {@value}.
   If an error is found when trying to translate, then an error message is appended to this file. 
  */
  public static final String TRANSLATION_ERRORS_FILE_NAME = "translation_errors.txt";
  
  /**
   Translate the given English text into the user's chosen language.
   
   If no translation is found, then the input English text is simply echoed. 
   In addition, the fact that no translation was found is saved to an errors file.
  */
  @Override protected String getEmittedText(String englishText) throws JspException, IOException {
    enusureInitialized();
    String result = englishText;
    HttpSession session = getRequest().getSession(DO_NOT_CREATE); 
    if (session != null){
      String langChoice = (String)session.getAttribute(LanguageDetector.KEY_NAME_SESSION);
      if (Util.textHasContent(langChoice) && Util.textHasContent(englishText)){
        result = translate(englishText, langChoice);
      }
    }
    return result;
  }
  
  // PRIVATE 

  /**
   This data structure handles N languages.
   Populated upon first use.
   The key is not an identifier in the usual sense. It's just the English version of the text.
   This unusual choice let's the JSP text remain more natural and legible, since it has 'real' text, instead of coder-keys. 
  */
  private static Map<String /*en text*/, Map<String /*fr, it, ...*/, String /*translated text*/>> TRANSLATIONS = null;
  
  private enum Language {
    ENGLISH("en", 1),
    FRENCH("fr", 2);
    //OTHER LANGS GO HERE
    private Language(String key, Integer groupIdx){
      KEY = key; //the same as the <html lang='en'> attribute: en, fr, etc
      GROUP_IDX = groupIdx; //the order of appearance in the translations file; a regex matching group index
    }
    String KEY;
    Integer GROUP_IDX;
  }
  
  private static final String ONE_LANG = "\\{([^{}]+)\\}";
  private static final String SPACER = "(?:\\s*)";
  private static final Boolean DO_NOT_CREATE = false;
  private static final String NL = System.getProperty("line.separator");
  private static final String FILE_ENCODING = "UTF-8";
  private static final Logger logger = Util.getLogger(TranslatorTag.class);
 
  private void enusureInitialized(){
    if (TRANSLATIONS == null){
      TRANSLATIONS = new LinkedHashMap<String, Map<String, String>>(); 
      readInTranslationsFileUponFirstUse();
    }
  }

  /** Used for dev testing only. */
  private static void main(String... args) throws IOException {
    log("Starting.");
  }
  
  private String translate(String englishText, String langChoice /*'en', 'fr', etc*/){
    String result = englishText.trim();
    if (!Language.ENGLISH.KEY.equals(langChoice)){
      Map<String, String> translations = TRANSLATIONS.get(result); //look up using the trimmed English text
      if (translations != null){
        String translation = translations.get(langChoice);
        if (Util.textHasContent(translation)){
          result = translation;
        }
        else {
          error2(englishText, langChoice);
        }
      }
      else {
        error1(englishText);
      }
    }
    log("Translated '" + englishText + "' to '" + result + "'");
    return result;
  }
  
  /** Builds a static, in-memory data structure. This data isn't changed after startup. Built upon first use. */
  private void readInTranslationsFileUponFirstUse(){
    log("Loading translations from " + TRANSLATIONS_FILE_NAME);
    int numDupes = 0;
    try {
      String allText = readTranslationsFile();
      Pattern pattern = regexForOneTranslationUnit();
      Matcher matcher = pattern.matcher(allText);
      while (matcher.find()) {
        String englishText = matcher.group(Language.ENGLISH.GROUP_IDX).trim();
        if (TRANSLATIONS.get(englishText) != null){
          ++numDupes;
          error3(englishText);
        }
        else {
          Map<String, String> translationsMap = new LinkedHashMap<String, String>(); 
          for (Language lang : Language.values()){
            if (!Language.ENGLISH.KEY.equals(lang.KEY)){
                translationsMap.put(lang.KEY, matcher.group(lang.GROUP_IDX).trim());
            }
          }
          TRANSLATIONS.put(englishText, translationsMap); //the key is the plain English text
        }
      }     
    } 
    catch (IOException e) {
      log("ERROR. Failed to read file");
      e.printStackTrace();
    }
    log("Number of languages: " + Language.values().length);
    log("Number of translations loaded from the translations file: " + TRANSLATIONS.size());
    log("Number of duplicate translations found: " + numDupes);
  }
  
  private Pattern regexForOneTranslationUnit(){
    String regex = "";
    for(int i = 0; i < Language.values().length; ++i){
      regex = regex + ONE_LANG;
      if (i < Language.values().length - 1){
        regex = regex + SPACER;
      }
    }
    return Pattern.compile(regex, Pattern.DOTALL);
  }
  
  private String readTranslationsFile() throws IOException {
    StringBuilder text = new StringBuilder();
    InputStream inputStream = TranslatorTag.class.getResourceAsStream(TRANSLATIONS_FILE_NAME);
    Scanner scanner = new Scanner(inputStream, FILE_ENCODING);
    try {
      while (scanner.hasNextLine()){
        text.append(scanner.nextLine() + NL);
      }
    }
    finally{
      scanner.close();
    }
    return text.toString();
  }
  
  private static void log(String text){
    //System.out.println(text);
    logger.fine(text);
  }
  
  private void error1(String englishText) {
    String errMsg = "Err1. English text not found: '" + englishText + "'";
    writeToTranslationsErrorFile(errMsg);
  }
  
  private void error2(String englishText, String targetLang){
    String errMsg = "Err2. No translation found for target lang: '" + targetLang + "'. English text: '" + englishText + "'";
    writeToTranslationsErrorFile(errMsg);
  }
  
  private void error3(String englishText){
    String errMsg = "Err3. Duplicate entry found for English text: '" + englishText + "'";
    writeToTranslationsErrorFile(errMsg);
  }
 
  private void writeToTranslationsErrorFile(String errMessage) {
    try {
      HttpSession session = getRequest().getSession(); //this might create a session, but that's ok
      String errorsFileName = (String)session.getServletContext().getRealPath("/WEB-INF/classes/astro/translate/" + TRANSLATION_ERRORS_FILE_NAME);
      boolean APPEND = true;
      Writer out = new OutputStreamWriter(new FileOutputStream(errorsFileName, APPEND), FILE_ENCODING);
      try {
        Date now = new Date();
        out.write(now.toString() + " Page:" + getPageName() + " " +  errMessage + NL);
      }
      finally {
        out.close();
      }    
    }
    catch(IOException ex){
      log("IOException. Can't log error: " + errMessage);
    }
  }
}