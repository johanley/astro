/*
  Take the contents of the Caldwell catalog, 
  and apply precession to the star positions, from J2000.0 to the 
  given target epoch.
*/
var show = function(input, text_output){
  console.log('Starting');
  
  EPH.testing();

  // WHY DO I NEED TO INLINE THIS DATA? WHY DOESN'T THE REFERENCE IN GRAPHIC.HTML WORK?  
  /* Caldwell catalog, below -35 deg declination. Generated on: Tue Dec 27 09:09:49 EST 2016. Name, Right Ascension (J2000), Declination (J2000), Mag, Constellation, Type, Comment (blank!), and Common Name.*/
  var caldwell = [42];
  caldwell[0]=["C68",4.9824787,-0.6448992,9.7,"CrA","Bn","","","NGC_6729"];
  caldwell[1]=["C69",4.5103671,-0.6475172,12.8,"Sco","Pl","","Bug Nebula","NGC_6302"];
  caldwell[2]=["C70",0.2395464,-0.6576982,8.1,"Scl","Sp","","","NGC_300"];
  caldwell[3]=["C71",2.0607975,-0.6728244,5.8,"Pup","Oc","","","NGC_2477"];
  caldwell[4]=["C72",0.0650135,-0.6838782,8.2,"Scl","Sb","","","NGC_55"];
  caldwell[5]=["C73",1.3705198,-0.6990044,7.3,"Col","Gc","","","NGC_1851"];
  caldwell[6]=["C74",2.6515915,-0.7056948,8.2,"Vel","Pl","","Eight Burst Nebula","NGC_3132"];
  caldwell[7]=["C75",4.3004913,-0.7097672,5.8,"Sco","Oc","","","NGC_6124"];
  caldwell[8]=["C76",4.4244097,-0.7295476,2.6,"Sco","Oc","","","NGC_6231"];
  caldwell[9]=["C77",3.5146568,-0.7507825,7.0,"Cen","Px","","Centaurus A","Centaurus_A"];
  caldwell[10]=["C78",4.7472956,-0.7627089,6.6,"CrA","Gc","","","NGC_6541"];
  caldwell[11]=["C79",2.6947884,-0.8101237,6.7,"Vel","Gc","","","NGC_3201"];
  caldwell[12]=["C80",3.5203291,-0.8287405,3.6,"Cen","Gc","","Omega Centauri","Omega_Centauri"];
  caldwell[13]=["C81",4.5618543,-0.8450302,8.1,"Ara","Gc","","","NGC_6352"];
  caldwell[14]=["C82",4.3689955,-0.8511389,5.2,"Ara","Oc","","","NGC_6193"];
  caldwell[15]=["C83",3.426954,-0.8633562,9.5,"Cen","Sp","","","NGC_4945"];
  caldwell[16]=["C84",3.6058502,-0.8965175,7.6,"Cen","Gc","","","NGC_5286"];
  caldwell[17]=["C85",2.2698007,-0.9261881,2.5,"Vel","Oc","","Omicron Vel Cluster","IC_2391"];
  caldwell[18]=["C86",4.6281768,-0.93666,5.6,"Ara","Gc","","","NGC_6397"];
  caldwell[19]=["C87",0.839067,-0.9637126,8.4,"Hor","Gc","","","NGC_1261"];
  caldwell[20]=["C88",3.9518618,-0.9704031,7.9,"Cir","Oc","","","NGC_5823"];
  caldwell[21]=["C89",4.271257,-1.0105456,5.4,"Nor","Oc","","S Norma Cluster","NGC_6087"];
  caldwell[22]=["C90",2.4495696,-1.0178178,9.7,"Car","Pl","","","NGC_2867"];
  caldwell[23]=["C91",2.9077185,-1.0239265,3.0,"Car","Oc","","","NGC_3532"];
  caldwell[24]=["C92",2.8091074,-1.0448704,6.2,"Car","Bn","","Eta Carinae Nebula","Carina_Nebula"];
  caldwell[25]=["C93",5.0217486,-1.0469067,5.4,"Pav","Gc","","","NGC_6752"];
  caldwell[26]=["C94",3.3754668,-1.0530153,4.2,"Cru","Oc","","","Jewel Box"];
  caldwell[27]=["C95",4.2049345,-1.0559242,5.1,"TrA","Oc","","","NGC_6025"];
  caldwell[28]=["C96",2.0869775,-1.0623237,3.8,"Car","Oc","","","NGC_2516"];
  caldwell[29]=["C97",3.0373092,-1.0754137,5.3,"Cen","Oc","","","NGC_3766"];
  caldwell[30]=["C98",3.3261612,-1.0989757,6.9,"Cru","Oc","","","NGC_4609"];
  caldwell[31]=["C99",3.3728488,-1.0995574,,"Cru","Dn","","Coalsack Nebula","Coalsack_Nebula"];
  caldwell[32]=["C100",3.0394909,-1.1001392,4.5,"Cen","Oc","","Lambda Centauri Nebula","IC_2944"];
  caldwell[33]=["C101",5.0169489,-1.1143927,9.0,"Pav","Sb","","","NGC_6744"];
  caldwell[34]=["C102",2.8064894,-1.123992,1.9,"Car","Oc","","Theta Car Cluster","IC_2602"];
  caldwell[35]=["C103",1.4778575,-1.2060225,1.0,"Dor","Bn","","Tarantula Nebula","Tarantula_Nebula"];
  caldwell[36]=["C104",0.275762,-1.2365658,6.6,"Tuc","Gc","","","NGC_362"];
  caldwell[37]=["C105",3.4016467,-1.2371476,7.3,"Mus","Gc","","","NGC_4833"];
  caldwell[38]=["C106",0.1051561,-1.2580915,4.0,"Tuc","Gc","","47 Tucanae","47_Tucanae"];
  caldwell[39]=["C107",4.3013639,-1.2601277,9.3,"Aps","Gc","","","NGC_6101"];
  caldwell[40]=["C108",3.2541664,-1.2682726,7.8,"Mus","Gc","","","NGC_4372"];
  caldwell[41]=["C109",2.6594454,-1.4113896,11.6,"Cha","Pl","","","NGC_3195"];

  
  console.log('Size of Caldwell catalog: ' + caldwell.length);
  var output_div = document.getElementById(text_output);
  
  var round = function(num){
	  var num_decimals = 10000000;
	  return Math.round(num*num_decimals)/num_decimals;
  };
  
  var output_result = function(nebula, idx, name, mag, constellation, type, comment /*always blank in this case*/, commonName, link){
    //caldwell[0]=["68",4.9824787,-0.6448992,9.7,"CrA","Bn","","Blah Nebula","NGC_1234"];
  	var line = 'caldwell['+idx+']=["'+name+'",' +round(nebula.α)+ ',' +round(nebula.δ)+ ',' +mag+ ',"' +constellation+ '","' +type+ '","' +comment+  '","'  +commonName+  '","' + link +  '"]';
  	//console.log(line);
  	output_div.innerHTML = output_div.innerHTML + line + ';' + '\r\n'; 
  };
  
  //var target = parseFloat(input.target_epoch);
  var when_target = EPH.when(input.target_epoch);
  var when_source = EPH.when_j2000;
  var angles = EPH.precession_angles(when_source, when_target);
  
  for (var idx = 0; idx < caldwell.length; ++idx){
  	var nebula = {
  	  α: caldwell[idx][1],
  	  δ: caldwell[idx][2],
  	  equinox: when_source
  	};
  	EPH.apply_precession(nebula, when_target, angles);
  	output_result(nebula, idx, caldwell[idx][0], caldwell[idx][3], caldwell[idx][4], caldwell[idx][5], caldwell[idx][6], caldwell[idx][7], caldwell[idx][8]);
  }
	
  console.log('Done');
};