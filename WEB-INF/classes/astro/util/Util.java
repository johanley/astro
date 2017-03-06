package astro.util;

import java.lang.reflect.Array;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.TimeZone;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 Static convenience methods for common tasks, which eliminate code duplication.

 <P>{@link Args} wraps certain methods of this class into a form suitable for checking arguments.
 
  <P>{@link hirondelle.web4j.util.WebUtil} includes utility methods particular to web applications.
*/
public final class Util {
  
  /**
   Return true only if <tt>aNumEdits</tt> is greater than <tt>0</tt>.
   
   <P>This method is intended for database operations.
  */
  public static boolean isSuccess(int aNumEdits){
    return aNumEdits > 0;
  }

  /**
    Return <tt>true</tt> only if <tt>aText</tt> is not null,
    and is not empty after trimming. (Trimming removes both
    leading/trailing whitespace and ASCII control characters. See {@link String#trim()}.)
   
    <P> For checking argument validity, {@link Args#checkForContent} should 
    be used instead of this method.
   
    @param aText possibly-null.
   */
   public static boolean textHasContent(String aText) {
     return (aText != null) && (aText.trim().length() > 0);
   }
   
  /**
   If <tt>aText</tt> is null, return null; else return <tt>aText.trim()</tt>.
  
   This method is especially useful for Model Objects whose <tt>String</tt>
   parameters to its constructor can take any value whatsoever, including 
   <tt>null</tt>. Using this method lets <tt>null</tt> params remain
   <tt>null</tt>, while trimming all others.
  
   @param aText possibly-null.
  */
   public static String trimPossiblyNull(String aText){
     return aText == null ? null : aText.trim();
   }
   
  /**
   Return <tt>true</tt> only if <tt>aNumber</tt> is in the range 
   <tt>aLow..aHigh</tt> (inclusive).
  
   <P> For checking argument validity, {@link Args#checkForRange} should 
   be used instead of this method.
  
   @param aLow less than or equal to <tt>aHigh</tt>.
  */
  static public boolean isInRange( int aNumber, int aLow, int aHigh ){
    if (aLow > aHigh) {
      throw new IllegalArgumentException(
        "Low: " + aLow + " is greater than High: " + aHigh
      );
    }
    return (aLow <= aNumber && aNumber <= aHigh);
  }
  
  /**
   Return <tt>true</tt> only if <tt>aAmount</tt> has exactly the number 
   of specified decimals.
   
    @param aNumDecimals is 0 or more.
  */
  public static boolean hasNumDecimals(BigDecimal aAmount, int aNumDecimals){
    if( aNumDecimals < 0 ){
      throw new IllegalArgumentException("Number of decimals must be 0 or more: " + quote(aNumDecimals));
    }
    return aAmount.scale() == aNumDecimals;
  }
  
  /**
   Parse text commonly used to denote booleans into a {@link Boolean} object.
   
   <P>The parameter passed to this method is first trimmed (if it is non-null), 
   and then compared to the following Strings, ignoring case :
   <ul>
   <li>{@link Boolean#TRUE} : 'true', 'yes', 'on'
   <li>{@link Boolean#FALSE} : 'false', 'no', 'off'
   </ul>
   
   <P>Any other text will cause a <tt>RuntimeException</tt>. (Note that this behavior 
   is different from that of {@link Boolean#valueOf(String)}).
  
   <P>(This method is clearly biased in favor of English text. It is hoped that this is not too inconvenient for the caller.)
   
   @param aBooleanAsText possibly-null text to be converted into a {@link Boolean}; if null, then the 
   return value is null. 
  */
  public static Boolean parseBoolean(String aBooleanAsText){
    Boolean result = null;
    String value = trimPossiblyNull(aBooleanAsText);
    if ( value == null ) {
      //do nothing - return null
    }
    else if ( value.equalsIgnoreCase("false") || value.equalsIgnoreCase("no") || value.equalsIgnoreCase("off") ) {
      result = Boolean.FALSE;
    }
    else if ( value.equalsIgnoreCase("true") || value.equalsIgnoreCase("yes") || value.equalsIgnoreCase("on") ) {
      result =  Boolean.TRUE;
    }
    else {
      throw new IllegalArgumentException(
        "Cannot parse into Boolean: " + quote(aBooleanAsText) + ". Accepted values are: true/false/yes/no/on/off"
      );
    }
    return result;
  }
  
  /**
   Coerce a possibly-<tt>null</tt> {@link Boolean}  value into {@link Boolean#FALSE}.
   
   <P>This method is usually called in Model Object constructors that have two-state <tt>Boolean</tt> fields.
    
   <P>This method is supplied specifically for request parameters that may be <em>missing</em> from the request, 
   during normal operation of the program.
   <P>Example : a form has a checkbox for 'yes, send me your newsletter', and the data is modeled has having two states - 
   <tt>true</tt> and <tt>false</tt>. If the checkbox is <em>not checked</em> , however, the browser will 
   likely not POST any corresponding request parameter - it will be <tt>null</tt>. In that case, calling this method 
   will coerce such <tt>null</tt> parameters into {@link Boolean#FALSE}.
   
   <P>There are other cases in which data is modeled as having not two states, but <em>three</em> : 
   <tt>true</tt>, <tt>false</tt>, and <tt>null</tt>. The <tt>null</tt> value usually means 'unknown'. 
   In that case, this method should <em>not</em> be called.
  */
  public static Boolean nullMeansFalse(Boolean aBoolean){
    return aBoolean == null ? Boolean.FALSE : aBoolean;
  }
  
  /**
   Create a {@link Pattern} corresponding to a <tt>List</tt>.
  
   Example: if the {@link List} contains "cat" and "dog", then the returned 
   <tt>Pattern</tt> will correspond  to the regular expression "(cat|dog)".
  
   @param aList is not empty, and contains objects whose <tt>toString()</tt>
   value represents each item in the pattern.
  */
  public static final Pattern getPatternFromList( List<?> aList ){
    if ( aList.isEmpty() ){
      throw new IllegalArgumentException();
    }
    StringBuilder regex = new StringBuilder("(");
    Iterator<?> iter = aList.iterator();
    while (iter.hasNext()){
      Object item = iter.next();
      regex.append( item.toString() );
      if ( iter.hasNext() ) { 
        regex.append( "|" );
      }
    }
    regex.append(")");
    return Pattern.compile( regex.toString() );
  }

  /**
   Return <tt>true</tt> only if <tt>aMoney</tt> equals <tt>0</tt>
   or <tt>0.00</tt>.
  */
  public static boolean isZeroMoney( BigDecimal aMoney ){
    final BigDecimal ZERO_MONEY = new BigDecimal("0");
    final BigDecimal ZERO_MONEY_WITH_DECIMAL = new BigDecimal("0.00");
    return 
      aMoney.equals(ZERO_MONEY) || 
      aMoney.equals(ZERO_MONEY_WITH_DECIMAL)
    ;
  }
  
  /**
   Return true only if <tt>aText</tt> is non-null, and matches 
   <tt>aPattern</tt>.
  
   <P>Differs from {@link Pattern#matches} and {@link String#matches},
   since the regex argument is a compiled {@link Pattern}, not a 
   <tt>String</tt>.
  */
  public static boolean matches(Pattern aPattern, String aText){
    /*
     Implementation Note:
     Patterns are thread-safe, while Matchers are not. Thus, a Pattern may 
     be compiled by a class once upon startup, then reused safely in a 
     multi-threaded environment. 
    */
    if (aText == null) return false;
    Matcher matcher = aPattern.matcher(aText);
    return matcher.matches();
  }
  
  /**
   Return true only if <tt>aText</tt> is non-null, and contains 
   a substring that matches <tt>aPattern</tt>.
  */
  public static boolean contains(Pattern aPattern, String aText){
    if (aText == null) return false;
    Matcher matcher = aPattern.matcher(aText);
    return matcher.find();
  }

  /**
   If <tt>aPossiblyNullItem</tt> is <tt>null</tt>, then return <tt>aReplacement</tt> ; 
   otherwise return <tt>aPossiblyNullItem</tt>.
   
   <P>Intended mainly for occasional use in Model Object constructors. It is used to 
   coerce <tt>null</tt> items into a more appropriate default value. 
  */
  public static <E> E replaceIfNull(E aPossiblyNullItem, E aReplacement){
    return aPossiblyNullItem == null ? aReplacement : aPossiblyNullItem;
  }
  
  /**
   Return a {@link Logger} whose name follows a specific naming convention.
  
   <P>The conventional logger names used by WEB4J are taken as   
   <tt>aClass.getPackage().getName()</tt>.
   
   <P>Logger names appearing in the <tt>logging.properties</tt> config file
   must match the names returned by this method.
  
   <P>If an application requires an alternate naming convention, then an  
   alternate implementation can be easily constructed. Alternate naming conventions 
   might account for :
  <ul>
   <li>pre-pending the logger name with the name of the application (this is useful 
   where log handlers are shared between different applications)
   <li>adding version information
  </ul>
  */
  public static Logger getLogger(Class<?> aClass){
    return Logger.getLogger(aClass.getPackage().getName());  
  }
  
  /**
   Call {@link String#valueOf(Object)} on <tt>aObject</tt>, and place the result in single quotes.
   <P>This method is a bit unusual in that it can accept a <tt>null</tt>
   argument : if <tt>aObject</tt> is <tt>null</tt>, it will return <tt>'null'</tt>.
   
   <P>This method reduces errors from leading and trailing spaces, by placing 
   single quotes around the returned text. Such leading and trailing spaces are 
   both easy to create and difficult to detect (a bad combination).
   
   <P>Note that such quotation is likely needed only for <tt>String</tt> data, 
   since trailing or leading spaces will not occur for other types of data.
  */
  public static String quote(Object aObject){
    String result = "EXCEPTION OCCURED";
    try {
      result = Consts.SINGLE_QUOTE + String.valueOf(aObject) + Consts.SINGLE_QUOTE; 
    }
    catch (Throwable ex){
      //errors were seen in this branch - depends on the impl of the passed object/array
      //do nothing - use the default 
    }
    return result;
  }
  
  /**
   Remove any initial or final quote characters from <tt>aText</tt>, either a single quote or 
   a double quote.
   
   <P>This method will not trim the text passed to it. Furthermore, it will examine only
   the very first character and the very last character. It will remove the first or last character, 
   but only if they are a single quote or a double quote.
    
   <P>If <tt>aText</tt> has no content, then it is simply returned by this method, as is, 
   including possibly-<tt>null</tt> values.
   @param aText is possibly <tt>null</tt>, and is not trimmed by this method
  */
  public static String removeQuotes(String aText){
    String result = null;
    if ( ! textHasContent(aText)) {
      result = aText;
    }
    else {
      int length = aText.length();
      String firstChar = aText.substring(0,1);
      String lastChar = aText.substring(length-1);
      boolean startsWithQuote = firstChar.equalsIgnoreCase("\"") || firstChar.equalsIgnoreCase("'"); 
      boolean endsWithQuote = lastChar.equalsIgnoreCase("\"") || lastChar.equalsIgnoreCase("'");
      int startIdx = startsWithQuote ? 1 : 0;
      int endIdx = endsWithQuote ? length-1 : length;
      result = aText.substring(startIdx, endIdx);
    }
    return result;
  }
  
  /**
   Replace every occurrence of a fixed substring with substitute text. 
    
   <P>This method is distinct from {@link String#replaceAll}, since it does not use a 
   regular expression.
  
   @param aInput may contain substring <tt>aOldText</tt>; satisfies {@link #textHasContent(String)}
   @param aOldText substring which is to be replaced; possibly empty, but never null
   @param aNewText replacement for <tt>aOldText</tt>; possibly empty, but never null
  */
  public static String replace(String aInput, String aOldText, String aNewText){
    if ( ! textHasContent(aInput) ) {
      throw new IllegalArgumentException("Input must have content.");
    }
    if ( aNewText == null ) {
      throw new NullPointerException("Replacement text may be empty, but never null.");
    }
    final StringBuilder result = new StringBuilder();
    //startIdx and idxOld delimit various chunks of aInput; these
    //chunks always end where aOldText begins
    int startIdx = 0;
    int idxOld = 0;
    while ((idxOld = aInput.indexOf(aOldText, startIdx)) >= 0) {
      //grab a part of aInput which does not include aOldPattern
      result.append( aInput.substring(startIdx, idxOld) );
      //add replacement text
      result.append( aNewText );
      //reset the startIdx to just after the current match, to see
      //if there are any further matches
      startIdx = idxOld + aOldText.length();
    }
    //the final chunk will go to the end of aInput
    result.append( aInput.substring(startIdx) );
    return result.toString();
  }

  /**
   Return a <tt>String</tt> suitable for logging, having one item from <tt>aCollection</tt>
   per line. 
   
   <P>For the <tt>Collection</tt> containing <br>
   <tt>[null, "Zebra", "aardvark", "Banana", "", "aardvark", new BigDecimal("5.00")]</tt>,
   
   <P>the return value is :
   <PRE>
   (7) {
     ''
     '5.00'
     'aardvark'
     'aardvark'
     'Banana'
     'null'
     'Zebra'
   }
   </PRE>
   
   <P>The text for each item is generated by calling {@link #quote}, and by appending a new line.
   
   <P>As well, this method reports the total number of items, <em>and places items in  
   alphabetical order</em> (ignoring case). (The iteration order of the <tt>Collection</tt> 
   passed by the caller will often differ from the order of items presented in the return value.) 
   </PRE>
  */
  public static String logOnePerLine(Collection<?> aCollection){
    int STANDARD_INDENTATION = 1;
    return logOnePerLine(aCollection, STANDARD_INDENTATION);
  }
  
  /**
   As in {@link #logOnePerLine(Collection)}, but with specified indentation level. 
   
   @param aIndentLevel greater than or equal to 1, acts as multiplier for a 
   "standard" indentation level of two spaces.
  */
  public static String logOnePerLine(Collection<?> aCollection, int aIndentLevel){
    String indent = getIndentation(aIndentLevel);
    StringBuilder result = new StringBuilder();
    result.append("(" + aCollection.size() + ") {" + Consts.NEW_LINE);
    List<String> lines = new ArrayList<String>(aCollection.size());
    for (Object item: aCollection){
      StringBuilder line = new StringBuilder(indent);
      line.append( quote(item) ); //nulls ok
      line.append( Consts.NEW_LINE );
      lines.add(line.toString());
    }
    addSortedLinesToResult(result, lines);
    result.append(getFinalIndentation(aIndentLevel));
    result.append("}");
    return result.toString();
  }
  
  /**
   Return a <tt>String</tt> suitable for logging, having one item from <tt>aMap</tt>
   per line. 
   
   <P>For a <tt>Map</tt> containing <br>
   <tt>["b"="blah", "a"=new BigDecimal(5.00), "Z"=null, null=new Integer(3)]</tt>,
   
   <P>the return value is :
   <PRE>
   (4) {
     'a' = '5.00'
     'b' = 'blah'
     'null' = '3'
     'Z' = 'null'
   }
   </PRE>
   
   <P>The text for each key and value is generated by calling {@link #quote}, and 
   appending a new line after each entry.
   
   <P>As well, this method reports the total number of items, <em>and places items in  
   alphabetical order of their keys</em> (ignoring case). (The iteration order of the 
   <tt>Map</tt> passed by the caller will often differ from the order of items in the
   return value.) 
   
   <P>An attempt is made to suppress the emission of passwords. Values in a Map are 
   presented as <tt>****</tt> if the following conditions are all true :
   <ul>
   <li>{@link String#valueOf(java.lang.Object)} applied to the <em>key</em> contains the word <tt>password</tt>
   (ignoring case)
   <li>the <em>value</em> is not an array or a <tt>Collection</tt>
   </ul>
  */
  public static String logOnePerLine(Map<?,?> aMap){
    StringBuilder result = new StringBuilder();
    result.append("(" + aMap.size() + ") {" + Consts.NEW_LINE);
    List<String> lines = new ArrayList<String>(aMap.size());
    String SEPARATOR = " = ";
    Iterator iter = aMap.keySet().iterator();
    while (  iter.hasNext() ){
      Object key = iter.next();
      StringBuilder line = new StringBuilder(INDENT);
      line.append(quote(key)); //nulls ok
      line.append(SEPARATOR);
      Object value = aMap.get(key);
      int MORE_INDENTATION = 2;
      if ( value != null && value instanceof Collection) {
        line.append( logOnePerLine((Collection)value, MORE_INDENTATION));
      }
      else if ( value != null && value.getClass().isArray() ){
        List valueItems = Arrays.asList( (Object[])value );
        line.append( logOnePerLine(valueItems, MORE_INDENTATION) );
      }
      else {
        value = suppressPasswords(key, value);
        line.append(quote(value)); //nulls ok
      }
      line.append(Consts.NEW_LINE);
      lines.add(line.toString());
    }
    addSortedLinesToResult(result, lines);
    result.append("}");
    return result.toString();
  }
  
  /**
   Return a {@link Locale} object by parsing <tt>aRawLocale</tt>.
   
   <P>The format of <tt>aRawLocale</tt> follows the 
   <tt>language_country_variant</tt> style used by {@link Locale}. The value is <i>not</i> checked against
   {@link Locale#getAvailableLocales()}. 
  */
  public static Locale buildLocale(String aRawLocale){
    int language = 0;
    int country = 1;
    int variant = 2;
    Locale result = null;
    fLogger.finest("Raw Locale: " + aRawLocale);
    String[] parts = aRawLocale.split("_");
    if (parts.length == 1) {
      result = new Locale( parts[language] );
    }
    else if (parts.length == 2) {
      result = new Locale( parts[language], parts[country] );
    }
    else if (parts.length == 3 ) {
      result = new Locale( parts[language], parts[country], parts[variant] );
    }
    else {
      throw new AssertionError("Locale identifer has unexpected format: " + aRawLocale);
    }
    fLogger.finest("Parsed Locale : " + Util.quote(result.toString()));
    return result;
  }
  
  /**
   Return a {@link TimeZone} corresponding to a given {@link String}.
   
   <P>If the given <tt>String</tt> does not correspond to a known <tt>TimeZone</tt> id,
   as determined by {@link TimeZone#getAvailableIDs()}, then a  
   runtime exception is thrown. (This differs from the behavior of the 
   {@link TimeZone} class itself, and is the reason why this method exists.)
  */
  public static TimeZone buildTimeZone(String aTimeZone) {
    TimeZone result = null;
    List<String> timeZones = Arrays.asList(TimeZone.getAvailableIDs());
    if( timeZones.contains(aTimeZone.trim()) ) {
      result = TimeZone.getTimeZone(aTimeZone.trim());
    }
    else {
      fLogger.severe("Unknown Time Zone : " + quote(aTimeZone));
      //fLogger.severe("Known Time Zones : " + logOnePerLine(timeZones));
      throw new IllegalArgumentException("Unknown TimeZone Id : " + quote(aTimeZone));
    }
    return result;
  }
  
  /**
   Convenience method for producing a simple textual
   representation of an array.
  
   <P>The format of the returned {@link String} is the same as 
   {@link java.util.AbstractCollection#toString} : 
   <ul>
   <li>non-empty array: <tt>[blah, blah]</tt>
   <li>empty array: <tt>[]</tt>
   <li>null array: <tt>null</tt>
   </ul>
  
   <P>Thanks to Jerome Lacoste for improving the implementation of this method.
   
   @param aArray is a possibly-null array whose elements are
   primitives or objects. Arrays of arrays are also valid, in which case
   <tt>aArray</tt> is rendered in a nested, recursive fashion.
   
  */
  public static String getArrayAsString(Object aArray){
    final String fSTART_CHAR = "[";
    final String fEND_CHAR = "]";
    final String fSEPARATOR = ", ";
    final String fNULL = "null";
    
    if ( aArray == null ) return fNULL;
    checkObjectIsArray(aArray);

    StringBuilder result = new StringBuilder( fSTART_CHAR );
    int length = Array.getLength(aArray);
    for ( int idx = 0 ; idx < length ; ++idx ) {
      Object item = Array.get(aArray, idx);
      if ( isNonNullArray(item) ){
        //recursive call!
        result.append( getArrayAsString(item) );
      }
      else{
        result.append( item );
      }
      if ( ! isLastItem(idx, length) ) {
        result.append(fSEPARATOR);
      }
    }
    result.append(fEND_CHAR);
    return result.toString();
  }

  /**
   Transform a <tt>List</tt> into a <tt>Map</tt>.
  
   <P>This method exists because it is sometimes necessary to transform a 
   <tt>List</tt> into a lookup table of some sort, using <em>unique</em> keys already
   present in the <tt>List</tt> data.
   
   <P>The <tt>List</tt> to be transformed contains objects having a method named 
   <tt>aKeyMethodName</tt>, and which returns objects of class <tt>aKeyClass</tt>.
   Thus, data is extracted from each object to act as its key. Furthermore, the 
   key must be <em>unique</em>. If any duplicates are detected, then an 
   exception is thrown. This ensures that the returned <tt>Map</tt> will be the 
   same size as the given <tt>List</tt>, and that no data is silently discarded. 
   
   <P>The iteration order of the returned <tt>Map</tt> is identical to the iteration order of 
   the input <tt>List</tt>.
  */
  public static <K,V> Map<K,V> asMap(List<V> aList, Class<K> aClass, String aKeyMethodName){
    Map<K,V> result = new LinkedHashMap<K,V>();
    for(V value: aList){
      K key = getMethodValue(value, aClass, aKeyMethodName);
      if( result.containsKey(key) ){
        throw new IllegalArgumentException("Key must be unique. Duplicate detected : " + quote(key));
      }
      result.put(key, value);
    }
    return result;
  }
  
  /**
   Reverse the keys and values in a <tt>Map</tt>.
   
   <P>This method exists because sometimes a lookup operation needs to be performed in a 
   style opposite to an existing <tt>Map</tt>.
   
   <P>There is an unusual requirement on the <tt>Map</tt> argument: the map <em>values</em> must be 
   unique. Thus, the returned <tt>Map</tt> will be the same size as the input <tt>Map</tt>. If any 
   duplicates are detected, then an exception is thrown.
    
   <P>The iteration order of the returned <tt>Map</tt> is identical to the iteration order of 
   the input <tt>Map</tt>.
  */
  public static <K,V> Map<V,K> reverseMap(Map<K,V> aMap){
    Map<V,K> result = new LinkedHashMap<V,K>();
    for(Map.Entry<K,V> entry: aMap.entrySet()){
      if( result.containsKey(entry.getValue())){
        throw new IllegalArgumentException("Value must be unique. Duplicate detected : " + quote(entry.getValue()));
      }
      result.put(entry.getValue(), entry.getKey());
    }
    return result;
  }
  
  
  // PRIVATE //
  
  private Util(){
    //empty - prevents construction by the caller.
  }
  
  private static final Logger fLogger = Util.getLogger(Util.class);
  private static final String INDENT = Consts.SPACE + Consts.SPACE;
  private static final Pattern PASSWORD = Pattern.compile("password", Pattern.CASE_INSENSITIVE); 
 
  private static void addSortedLinesToResult(StringBuilder aResult, List<String> aLines) {
    Collections.sort(aLines, String.CASE_INSENSITIVE_ORDER);
    for (String line: aLines){
      aResult.append( line );
    }
  }
  
  private static String getIndentation(int aIndentLevel){
    StringBuilder result = new StringBuilder();
    for (int idx = 1; idx <= aIndentLevel; ++idx){
      result.append(INDENT);
    }
    return result.toString();
  }
  
  private static String getFinalIndentation(int aIndentLevel){
    return getIndentation(aIndentLevel - 1);
  }

  /**
   Replace likely password values with a fixed string. 
  */
  private static Object suppressPasswords(Object aKey, Object aValue){
    Object result = aValue;
    String key = String.valueOf(aKey);
    Matcher matcher = PASSWORD.matcher(key);
    if ( matcher.find() ){
      result = "*****";
    }
    return result;
  }
  
  private static void checkObjectIsArray(Object aArray){
    if ( ! aArray.getClass().isArray() ) {
      throw new IllegalArgumentException("Object is not an array.");
    }
  }

  private static boolean isNonNullArray(Object aItem){
    return aItem != null && aItem.getClass().isArray();
  }

  private static boolean isLastItem(int aIdx, int aLength){
    return (aIdx == aLength - 1);
  }
  
  private static <K> K getMethodValue(Object aValue, Class<K> aClass, String aKeyMethodName){
    K result = null;
    try {
      Method method = aValue.getClass().getMethod(aKeyMethodName); //no args
      result = (K)method.invoke(aValue);
    }
    catch (NoSuchMethodException ex){
      handleInvocationEx(aValue.getClass(), aKeyMethodName);
    }
    catch (IllegalAccessException ex){
      handleInvocationEx(aValue.getClass(), aKeyMethodName);
    }
    catch (InvocationTargetException ex){
      handleInvocationEx(aValue.getClass(), aKeyMethodName);
    }
    return result;
  }
  
  private static void handleInvocationEx(Class<?> aClass, String aKeyMethodName){
    throw new IllegalArgumentException("Cannot invoke method named " + quote(aKeyMethodName) + " on object of class " + quote(aClass));
  }
}
