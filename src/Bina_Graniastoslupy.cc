#include "Bina_Graniastoslupy.hh"

G4VSolid* przecinka(G4VSolid* blokG2_sol,wymiary &w,int wys,int typ)
{
  // 
// wys = 1 for thick detektor  E  , wys = 2 for thin de 
// typ -  special cuts

  G4double a1;//,a2;             //upper and lower base
  
  G4double h1;//,h2;             //hight of upper and lower base
//  G4double h;                      // hight  of the scintillator
  
  G4double /*ak,bk,akr,*/bkr,/*gk,*/gkr;            //angles (degrees and  radians)
  
  G4double /*d, d1, */q, c,b,zt;               //zmienne pomoocnicze

  G4double dv1,dv2;                        // cuts 
  
  
  if(wys==2)
    {
      a1=w.get_a1m();
//      a2=w.get_a2m();
      h1=w.get_h1m();
//      h=w.get_hm();
//      ak=w.get_ak();
//      akr=w.get_akr();
//      d1=w.get_d1m();		
//      h2=w.get_h2m();		
//      d =w.get_dm();		
      bkr=w.get_bkr();
//      bk=w.get_bk();		
      b=w.get_bm();		
      zt = w.get_ztm();    //!!!!!!!!!
      q=w.get_qm();		
      c=w.get_cm();		
//      gk= w.get_gk();
      gkr=w.get_gkr();
      dv1=w.get_dv1m();
      dv2=w.get_dv2m();
    }
  else
    {
      a1=w.get_a1();
//      a2=w.get_a2();
      h1=w.get_h1();
//      h=w.get_h();
//      ak=w.get_ak();
//      akr=w.get_akr();
//      d1=w.get_d1();		
//      h2=w.get_h2();		
//      d =w.get_d();		
      bkr=w.get_bkr();
//      bk=w.get_bk();		
      b=w.get_b();		
      zt = w.get_zt();    //!!!!!!!!!
      q=w.get_q();		
      c=w.get_c();	
//      gk= w.get_gk();
      gkr=w.get_gkr();
      dv1=w.get_dv1();
      dv2=w.get_dv2();
    }
  
  G4Trd*  blok_sol;
  G4SubtractionSolid*  blokG1_sol;
  G4Box*  blok_sol1;
  G4Box*  blok_sol2;
  G4Box*  blok_sol3;
  G4Box*  blok_sol4;
  G4ThreeVector position;
  
  G4double xt,xp,yt,yp;
  G4double  k_fin, kr_fin, y_fin,  x_fin, z_fin;
  G4double  kdv,rdv;
  
  // Main shape
  xt=a1;
  xp=0.00001;//0.00001;
  yt=b+c;
  yp=b;
  G4double ztp=4.1*yt,xtp=4.1*xt,ytp=4.1*zt;
  
  G4double zzt;
  zzt=zt+0.1;//005;

  blok_sol3=new G4Box("Glowny_Szescian",xtp*CLHEP::cm,ytp*CLHEP::cm,ztp*CLHEP::cm);

  blok_sol=new G4Trd("Trapez",xt*CLHEP::cm,xp*CLHEP::cm,yt*CLHEP::cm,yp*CLHEP::cm,zt*CLHEP::cm);
  
  if(wys==2&&(typ==1||typ==2||typ==3))
    {
      if(typ==1)
	{
	  //dv1=1.338;       dv2=0.812;                     
	  dv1=dv1*cos(bkr);  dv2=dv2*cos(bkr);        // cuts
	  
	  rdv=dv1-dv2;                 // difference between positions of cuts
	  kdv=atan2(rdv,2.*b);                              // resulting angle
	  dv1=dv1;
	  dv1=dv1/cos(kdv);
	  y_fin=dv1;
	  x_fin=z_fin=xtp;
	  blok_sol4=new G4Box("Szescian_Profilujacy",(x_fin)*CLHEP::cm,(y_fin)*CLHEP::cm,(z_fin)*CLHEP::cm);
	  position = G4ThreeVector((0.)*CLHEP::cm,(0.)*CLHEP::cm,(-b)*CLHEP::cm);
	  G4RotationMatrix* rotFin = new G4RotationMatrix();     rotFin->rotateX(kdv);
	  blokG1_sol=new G4SubtractionSolid("Obciecie_specjalne1",blok_sol3,blok_sol4,rotFin,position);
	}
      
      if(typ==2)
	{
	  k_fin=90.-40.1;                     //rotation around Z
	  kr_fin=(k_fin*M_PI)/180.;           //in radians
	  // dv1=1.338;
	  y_fin=dv1/cos(((30.-11.5)*M_PI)/180.);  // 30.-11.5 rotation around x X
	  x_fin=z_fin=xtp;
	  blok_sol4=new G4Box("Szescian_Profilujacy",(x_fin)*CLHEP::cm,(y_fin)*CLHEP::cm,(z_fin)*CLHEP::cm);
	  position = G4ThreeVector((0.)*CLHEP::cm,(0.)*CLHEP::cm,(-b)*CLHEP::cm);
	  G4RotationMatrix* rotFin = new G4RotationMatrix();    
	  rotFin->rotateX((4.9*M_PI)/180.);   
	  rotFin->rotateZ(kr_fin);
	  blokG1_sol=new G4SubtractionSolid("Obciecie_specjalne2",blok_sol3,blok_sol4,rotFin,position);
	}
      
      if(typ==3)
	{
	  k_fin=90.-40.1;                        //rotation around Z
	  kr_fin=(k_fin*M_PI)/180.;              //in radianacs
	  //  dv1=1.338;
	  y_fin=dv1/cos(((30.-11.5)*M_PI)/180.);  // 30.-11.5 rotation around X
	  x_fin=z_fin=xtp;
	  blok_sol4=new G4Box("Szescian_Profilujacy",(x_fin)*CLHEP::cm,(y_fin)*CLHEP::cm,(z_fin)*CLHEP::cm);
	  position = G4ThreeVector((0.)*CLHEP::cm,(0.)*CLHEP::cm,(-b)*CLHEP::cm);
	  G4RotationMatrix* rotFin = new G4RotationMatrix();    rotFin->rotateX((4.9*M_PI)/180.);   rotFin->rotateZ(-kr_fin);
	  blokG1_sol=new G4SubtractionSolid("Obciecie_specjalne3",blok_sol3,blok_sol4,rotFin,position);
     }

     //Cuts on top
      G4double xb=0.,yb=0.,zb=0.,cb=0.;
      cb = 2.1 * c;
      xb=2.1 * a1;   zb= 4. * h1;  yb=cb;
      
      blok_sol1=new G4Box("Blok_gorny",(xb)*CLHEP::cm,(yb)*CLHEP::cm,(zb)*CLHEP::cm);
      position = G4ThreeVector(0.*CLHEP::cm,(-zt)*CLHEP::cm,-(b-q+cb/cos(bkr))*CLHEP::cm);
      
      G4RotationMatrix* rotTop = new G4RotationMatrix();     rotTop->rotateX(bkr+M_PI/2.);
      blokG1_sol=new G4SubtractionSolid("Obciecie_gorne",blokG1_sol,blok_sol1,rotTop,position);
      

      /// cut on a face surface
      G4double e,e1;
      e=2.*b*tan(gkr);
      e1=1.1 * e;
      xb=1.1 * a1;   zb= e1;  yb=2. * yt;
      G4double men=0;
      men =e1/cos(gkr);
      blok_sol2=new G4Box("Blok_czolo",(xb)*CLHEP::cm,(yb)*CLHEP::cm,(zb)*CLHEP::cm);
      position = G4ThreeVector(0.*CLHEP::cm,(-men-zt-zt)*CLHEP::cm,(-b+2.*q)*CLHEP::cm);
      G4RotationMatrix* rotFront = new G4RotationMatrix();     rotFront->rotateX(gkr+M_PI/2.);
      blokG1_sol=new G4SubtractionSolid("Obciecie_czolowe",blokG1_sol,blok_sol2,rotFront,position);
      
      position = G4ThreeVector(0.*CLHEP::cm,(-zzt)*CLHEP::cm,(0.)*CLHEP::cm );
      G4RotationMatrix* rotF = new G4RotationMatrix();     rotF->rotateX(M_PI/2.);
      blokG2_sol=new G4IntersectionSolid("Dodanie_Szescian",blokG1_sol,blok_sol,rotF,position);
    }
  
  else
    {
      G4double xb=0.,yb=0.,zb=0.,cb=0.;
      cb = 2.1 * c;
      xb=2.1 * a1;   zb= 4. * h1;  yb=cb;
      
      blok_sol1=new G4Box("Blok_gorny",(xb)*CLHEP::cm,(yb)*CLHEP::cm,(zb)*CLHEP::cm);
      position = G4ThreeVector(0.*CLHEP::cm,(-zt)*CLHEP::cm,-(b-q+cb/cos(bkr))*CLHEP::cm);
      
      G4RotationMatrix* rotTop = new G4RotationMatrix();
      rotTop->rotateX(bkr+M_PI/2.);
      blokG1_sol=new G4SubtractionSolid("Obciecie_gorne",blok_sol3,blok_sol1,rotTop,position);
      
      ///  cut on a face surface
      G4double e,e1;
      e=2.*b*tan(gkr);
      e1=1.1 * e;
      xb=1.1 * a1;   
      zb= e1;  
      yb=2. * yt;
      G4double men=0;
      men =e1/cos(gkr);
      blok_sol2=new G4Box("Blok_czolo",(xb)*CLHEP::cm,(yb)*CLHEP::cm,(zb)*CLHEP::cm);
      position = G4ThreeVector(0.*CLHEP::cm,(-men-zt-zt)*CLHEP::cm,(-b+2.*q)*CLHEP::cm);
      G4RotationMatrix* rotFront = new G4RotationMatrix();
      rotFront->rotateX(gkr+M_PI/2.);
      blokG1_sol=new G4SubtractionSolid("Obciecie_czolowe",blokG1_sol,blok_sol2,rotFront,position);
      
      position = G4ThreeVector(0.*CLHEP::cm,(-zzt)*CLHEP::cm,(0.)*CLHEP::cm );
      G4RotationMatrix* rotF = new G4RotationMatrix();
      rotF->rotateX(M_PI/2.);
      blokG2_sol=new G4IntersectionSolid("Dodanie_Szescian",blokG1_sol,blok_sol,rotF,position);
      // Modification of shape  for non-typical  elements at the end of the ball
      if(typ==1)
	{
	  //dv1=1.338;       dv2=0.812;               
	  dv1=dv1*cos(bkr);  dv2=dv2*cos(bkr);        
	  rdv=dv1-dv2;                                
	  kdv=atan2(rdv,2.*b);                        
	  dv1=dv1;
	  dv1=dv1/cos(kdv);
	  y_fin=dv1;
	  x_fin=z_fin=xtp;
	  blok_sol4=new G4Box("Szescian_Profilujacy",(x_fin)*CLHEP::cm,(y_fin)*CLHEP::cm,(z_fin)*CLHEP::cm);
	  position = G4ThreeVector((0.)*CLHEP::cm,(0.)*CLHEP::cm,(-b)*CLHEP::cm);
	  G4RotationMatrix* rotFin = new G4RotationMatrix();
	  rotFin->rotateX(kdv);
	  blokG2_sol=new G4SubtractionSolid("Obciecie_specjalne1",blokG2_sol,blok_sol4,rotFin,position);
	}
      
      if(typ==2)
	{
	  k_fin=90.-40.1;                                          
	  kr_fin=(k_fin*M_PI)/180.;                             
	  // dv1=1.338;
	  y_fin=dv1/cos(((30.-11.5)*M_PI)/180.);  
	  x_fin=z_fin=xtp;
	  blok_sol4=new G4Box("Szescian_Profilujacy",(x_fin)*CLHEP::cm,(y_fin)*CLHEP::cm,(z_fin)*CLHEP::cm);
	  position = G4ThreeVector((0.)*CLHEP::cm,(0.)*CLHEP::cm,(-b)*CLHEP::cm);
	  G4RotationMatrix* rotFin = new G4RotationMatrix(); 
	  rotFin->rotateX((4.9*M_PI)/180.); 
	  rotFin->rotateZ(kr_fin);
	  blokG2_sol=new G4SubtractionSolid("Obciecie_specjalne2",blokG2_sol,blok_sol4,rotFin,position);
	}
      
      if(typ==3)
	{
	  k_fin=90.-40.1;                                          //rot tracego po osi Z
	  kr_fin=(k_fin*M_PI)/180.;                          //w radianach
	  //  dv1=1.338;
	  y_fin=dv1/cos(((30.-11.5)*M_PI)/180.);  // 30.-11.5 obrot po osi X
	  x_fin=z_fin=xtp;
	  blok_sol4=new G4Box("Szescian_Profilujacy",(x_fin)*CLHEP::cm,(y_fin)*CLHEP::cm,(z_fin)*CLHEP::cm);
	  position = G4ThreeVector((0.)*CLHEP::cm,(0.)*CLHEP::cm,(-b)*CLHEP::cm);
	  G4RotationMatrix* rotFin = new G4RotationMatrix(); 
	  rotFin->rotateX((4.9*M_PI)/180.); 
	  rotFin->rotateZ(-kr_fin);
	  blokG2_sol=new G4SubtractionSolid("Obciecie_specjalne3",blokG2_sol,blok_sol4,rotFin,position);
	}
    }
  return blokG2_sol;
}


wymiary::wymiary(G4double A1,G4double H1,G4double A2,G4double H,G4double AK)
{
  a1=A1;
  a2=A2;
  h1=H1;
  h=H;
  ak=AK;
  
  licz();
}

void wymiary::licz()
{
  ////// values for the first G (big)
  
  akr=ak*CLHEP::pi/180.;
  //hx = h/2.; 
  d1=2*h*tan(akr);           
  h2 = (a2*h1) / a1;         
  d = 2.*(h1 - h2) -d1;      
  bkr=atan2(d,2.*h);
  bk=180.*bkr/CLHEP::pi;    
  b=h/cos(bkr);      
  zt = h1*cos(bkr);
  q=zt*tan(bkr);     
  c=q+q;             
  gk= ak+bk;
  gkr=gk*CLHEP::pi/180.;
  
  dv1=1.338;     dv2=0.812;                      // wartosci cienc w plaszczyznie podstaw gora i dol
  
  ///// values for second G (small)
  hm=0.05;             // heigh of the thin layer is the same for all elements
  
  h1m=h2;
  a1m=a2;
  a2m=a1-(((a1-a2)*(h+hm))/(h)); //   (a2*(r-b))/(r-b-hm);
  
  d1m=2*hm*tan(akr); 
  h2m = (a2m*h1m) / a1m; 
  dm = 2.*(h1m - h2m) -d1m;
  bm=hm/cos(bkr);                        
  ztm = h1m*cos(bkr);
  qm=ztm*tan(bkr);                       
  cm=qm+qm;                              
  // gk= ak+bk;       gkr=gk*pi/180.;
  
  G4double tmp,dr=0.172;
  tmp=h2m/3.;
  if(h>3.)
    {
      rm=(12.-0.075)+2.*tmp*sin(bkr)+bm;   
    }
  else
    {
      rm=(18.-0.075)+2.*tmp*sin(bkr)+bm;   
    }
  if(h1<3.5)
    {
      rm=rm+((dr*a2m)/a1);
    }
  r=rm+bm+b+0.001;    ////G4cout<<"\n\t********  R_ =" <<  21.-tmp*sin(akr);          //20.072;
  
  dv1m=dv2;    
  dv2m=(dv2*a2m)/a2;//(dv2*(rm-bm)/(r-b));

  }
/*
  G4double wymiary::get(char *nazwa)
  {
  for(int i=0;i<n_zmiennych;i++)
  {
  if(zmienna[i]==nazwa)
  { return wart[i];  i=n_zmiennych; }
  }
 } */
/*
  void wymiary::set()
  {}
*/
