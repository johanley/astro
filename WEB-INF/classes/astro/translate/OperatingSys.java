package astro.translate;

public enum OperatingSys {

  WINDOWS("\r\n"),
  UNIX("\n");
  
  private OperatingSys(String newline){
    NEWLINE = newline; 
  }
  String NEWLINE;
}
