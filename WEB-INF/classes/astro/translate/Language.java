package astro.translate;

enum Language {

  ENGLISH("en", 1),
  FRENCH("fr", 2);
  //NOTE TO THE FUTURE: OTHER LANGS GO HERE
  
  private Language(String key, Integer groupIdx){
    KEY = key; //the same as the <html lang='en'> attribute: en, fr, etc
    GROUP_IDX = groupIdx; //the order of appearance in the translations file; a regex matching group index
  }
  String KEY;
  Integer GROUP_IDX;
}
