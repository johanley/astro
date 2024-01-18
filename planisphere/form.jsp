<!doctype html>
<html lang='<s:txt>en</s:txt>'>
<head>
 <tags:analytics/>
 <meta charset='UTF-8'>
 <meta name="keywords" content="<s:txt>astronomy, planisphere, high precision</s:txt>">
 <meta name="description" content="<s:txt>Customizable planisphere with high precision.</s:txt>">
 <meta name="viewport" content="width=device-width">
 <link rel="stylesheet" href="../css/styles.css<tags:ver/>" media="all">
 <title><s:txt>Customizable Planisphere</s:txt></title>
</head>
<body>
<h2><s:txt>Customizable Planisphere</s:txt> - <a href='..'>Home</a></h2>

 <h3>1. Generate</h3> 
 Generate two PDF files, one at a time, using the exact same settings in the form.
 <P>One file is a star chart, and the other is a transparency.
 <P>The default form shows a working example. Overwrite it with your own values. 
 <P>
 <form method='GET' action='form.pln' class='user-input-small' id='input_form'>
  <table>
     <tr>
       <td><s:txt>Year</s:txt>*:
       <td><input type='text' required id='year' name='year' value='2022' title='<s:txt>Year used in generating the planisphere</s:txt>'>
       <td>The year for which the planisphere will be generated. 
     <tr>
       <td><s:txt>Location</s:txt>*:
       <td><input type='text' required id='location' name='location' value='Halifax, NS' title='<s:txt>Location of the observer</s:txt>'>
       <td>The name of the location where the planisphere will be used.
     <tr>
       <td><s:txt>Latitude</s:txt>*:
       <td><input type='text' id='latitude' required name='latitude' value='44.62' title='<s:txt>Latitude in degrees</s:txt>'>
       <td>The latitude of the location. Decimal degrees, negative for south.
     <tr>
       <td><s:txt>Longitude</s:txt>*:
       <td><input type='text' id='longitude' required name='longitude' value='-63.61' title='<s:txt>Longitude in degrees. Negative west of Greenwich</s:txt>'>
       <td>The longitude of the location. Decimal degrees, negative for west.
     <tr>
       <td><s:txt>Hours offset from UT</s:txt>*:
       <td><input type='text' id='hours-offset-from-ut' required name='hoursOffsetFromUT' value='-4' title='<s:txt>Negative west of Greenwich</s:txt>'>
       <td>The whole number of hours the location is from Greenwich. Negative for west.
     <tr>
       <td><s:txt>Declination gap</s:txt>*:
       <td><input type='text' id='declinationGap' required name='declinationGap' value='15.0' title='<s:txt>Amount of declination to chop off in the south (north)</s:txt>'>
       <td>Decimal degrees of declination to chop off in the south (or north).
      <tr>
       <td><s:txt>Grey for alt-az lines</s:txt>*:
       <td><input type='text' id='grey-scale-alt-az' required name='greyScaleAltAz' value='112' title='<s:txt>0..128, dark..light</s:txt>'>
       <td>0..128, grey value for alt-az lines.
     <tr>
       <td><s:txt>1. Generate 1st file</s:txt>:
       <td style="text-align:center"><input type='submit' name='generate' value='Star chart'>
       <td>Downloads a PDF file containing the star chart.
     <tr>
       <td><s:txt>2. Generate 2nd file</s:txt>:
       <td style="text-align:center"><input type='submit' name='generate' value='Transparency'>
       <td>Downloads a PDF file containing the transparency.
  </table>
 </form>
 

 <h3>2. Print</h3> 
 <ul>
  <li>print the generated PDF files on a <b>high-resolution printer</b> (you may need to use a commercial print shop for best results).
  <li>print <code>the starchart.pdf</code> file on card stock (thick paper). One page, double-sided.
  <li>print the <code>transparency.pdf</code> file on a transparency. One page, single-sided.
 </ul>
  
 <h3>3. Cut, puncture, and attach</h3> 
  <ul>
  <li>make the star chart square, by cutting the star chart along the indicated horizontal lines above and below the star chart. 
  <li>make the transparency circular, by cutting just outside the outermost circular border
  <li><b>carefully puncture with a needle</b> the transparency and the start chart, using the affordances at the position of the celestial pole.
  Do them separately.
  Use a needle to make the initial holes, and widen them carefully. It's important to be as precise as possible! 
  <li>attach the transparency to the star chart using a pin, <a href='https://www.mcmaster.com/two-piece-rivets/'>rivet</a>, or something similar.
  The outermost circular borders need to line up as precisely as possible.
 </ul>

<P><s:txt>See the <a href='https://johanley.github.io/planisphere/index.html'>main documentation</a> for further information.</s:txt> 

</body>
</html>