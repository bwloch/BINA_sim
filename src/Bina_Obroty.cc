// obroty.cc
#include "Bina_Graniastoslupy.hh"

 polozenia::polozenia()
{typ=-1; x=0.; y=0.; z=0.; rx=0.;ry=0.;rz=0.;r=0.;}

 void planowanie(polozenia* plac,double R[],int n_el)
 {
//  double r;//,x,y,z;

  double r51=R[0],
    r511=R[1],
    r61=R[4],
    r611=R[5],
    r52=R[2],
    r521=R[3],
    r62=R[6],
    r621=R[7];

 G4double sgkr,sfkr;

 int w_max =7;           // number of layers:   7
 G4double gkr[10];       // gamma angle of the first element in the layer


 gkr[1]=(37.4*M_PI)/180.;
 gkr[2]=(63.44*M_PI)/180.;
 gkr[3]=(79.2*M_PI)/180.;
 gkr[4]=(100.814*M_PI)/180.;
 gkr[5]=(116.58*M_PI)/180.;			
 gkr[6]=(142.623*M_PI)/180.;

 G4double fkr[10];        //phi angle  of the first element in the layer
                          // (actual staring value is shifted by 90 degrees)

 fkr[1]=0.;    fkr[2]=M_PI/5.;  fkr[3]=0.;  fkr[4]=fkr[2];  
 fkr[5]=0.;  fkr[6]=fkr[2];

 G4double k60,k72,kat,k360,k90;
 k60=M_PI/3.;
 k72=(2.*M_PI)/5.;
 k360=2.*M_PI;
 k90=M_PI/2.;

 int n=0,nn=0;                // numbering of elements,  blocking of numbering

 for(int w=1;w<w_max;w++)     // loop over layers
   {
     if (n>=n_el)      
       {
       G4cout<<"\n too few elements of BINA declared";
       break ;
       }
     sgkr=gkr[w];                       // gamma angle
     sfkr =k90+fkr[w];                  // starting phi value
     
     while(sfkr<k360+k72)               // loop over phi angle; 
                                        // phi for centers of hexagons or pentagons
                                        // is calculared
       {
       bool parameter=0;
       if(sgkr>=M_PI) {sgkr=sgkr-2*M_PI; parameter=1;}
//	 x=r*sin(sgkr)*cos(sfkr);	// positions with respect to the center
//	 y=r*sin(sgkr)*sin(sfkr);       // of the ball
//	 z=r*cos(sgkr);
	 if(parameter==1) sgkr=sgkr+2*M_PI;
	 kat=0.;
	 while(kat<k360-0.1)           //rotation
	   {
	     nn=0;
	     
	     plac[n].y=M_PI-sgkr;
	     plac[n].z=sfkr;          //kat theta w ktorym jest umieszczone centrum piecio/szesciokata
	     plac[n].rx=sgkr*CLHEP::rad;       // Theta obrotu piecio/szesciokata
             plac[n].rz=(sfkr-k90)*CLHEP::rad; // Phi obrotu 5/6kata
	     if(w==1)                   //        !!!if first layer!!!
	       {                        // loop over elements of the hexagon 
		 plac[n].ry=0.+kat;     // Psi //kat skladajacy detektory w 5 lub 6  
		 kat=kat + k60;         //   (from 0 to 6 * 60 degrees)
		 plac[n].typ=2;         // hexagons are made of short E's 
		 plac[n].x=r61;
	       }
	     
	     if(w==2)                   //        !!! if second layer!!
	       {                        // loop over elements of the pentagon
		 plac[n].ry=0.+kat+M_PI;// Psi
		 kat=kat + k72;         //   (from 0 to 5 * 72 degrees)
		 plac[n].typ=0;         // pentagons are made of short E's 
		 plac[n].x=r51;
	       }
	     
	     if(w==3)
	       {
		 plac[n].ry=0.+kat;     // Psi
		 if((kat>k72)&&(kat<M_PI+k72))
		   {
		     plac[n].typ=3;     // hexagons are made of long E's
		     plac[n].x=r62;
		   }
		 else
		   {
		     plac[n].typ=2;     // heagons are made of  short E's
		     plac[n].x=r61;
		   }
		 kat=kat + k60;
	         if(sfkr>k90-0.02&&sfkr<k90+0.02)
		   { plac[n].typ=-1; }
	       }
	     
	     if(w==4)
	       {
		 plac[n].ry=0.+kat;     // Psi
		 kat=kat + k60;  
		 plac[n].typ=3;         // heagons are made of long E's
		 plac[n].x=r62;        
	       }
	     
	     if(w==5)
	       {
		 plac[n].ry=0.+kat;     // Psi
		 kat=kat + k72;  
		 plac[n].typ=1;         // pentagons made of short E's  
		 plac[n].x=r52;
	       }
	     if(w==6)
	       {
		 plac[n].ry=0.+kat;     // Psi
		 if(kat<0.1&&kat>(-0.1))
		   {
		     plac[n].typ=4;     // heagons are made of long E's 	
		     plac[n].x=r62;     // (special shape 1)

		   }
		 else if(kat<k60+0.1&&kat>k60-0.1)
		   {
		     plac[n].typ=5;     // heagons are made of  long E's
		     plac[n].x=r62;     // 
		   }
		 else if(kat<(2*M_PI-k60+0.1)&&kat>(2*M_PI-k60-0.1))
		   {
		     plac[n].typ=6;     // heagons are made of long E's 
		     plac[n].x=r62;     // (special shape 3)
		   }
		 else
		   { 
		     nn=1; 
		     plac[n].typ=-1;
		   }
		 kat=kat + k60;
	       }
	     if(nn!=1)
	       { n++; }
	   }
         sfkr=sfkr + k72;
       }
   }
#define male 1
#if (male==1)
 int t;
 for(int i=0;i<n;i++)
   {
     if (n>=n_el)  
       {
	 G4cout<<"\n too few elements of BINA declared";
	 break ;
       }
     t=i+n;
     plac[t].y=plac[i].y;
     plac[t].z=plac[i].z;
     plac[t].rx=plac[i].rx;      // R
     if(plac[i].x==r51)
       {
	 plac[t].x=r511;
       }
     else if(plac[i].x==r61)
       {
	 plac[t].x=r611;
       }
     else if(plac[i].x==r52)
       {
	 plac[t].x=r521;
       }
     else
       {
	 plac[t].x=r621;
       }
     plac[t].rz=plac[i].rz;      //Phi
     plac[t].ry=plac[i].ry;      //Psi
     plac[t].typ=plac[i].typ+10;
   }
#endif
 
 }

