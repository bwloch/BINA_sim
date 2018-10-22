#include "Bina_PrimaryGeneratorAction.hh"
#include "Bina_DetectorConstruction.hh"
//#include "Bina_PrimaryGeneratorActionMessenger.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"

#include <fstream>
#include "Randomize.hh"
//#include "Ranlux64Engine.h"
#include <cmath>

#ifndef min
#define min(x,y) ((x) < (y) ? (x) : (y))
#endif
#ifndef max
#define max(x,y) ((x) > (y) ? (x) : (y))
#endif
#include <iostream>

//fstream filell;
//int chooser=0;
double csmax;
int toSix=0;
double csilocmax=0;

  double step[4] = {2.0, 2.0, 15.0,0.002};
  double offs[4] = {5.0, 5.0, 0.0, 0.0};


Bina_PrimaryGeneratorAction::Bina_PrimaryGeneratorAction(Bina_DetectorConstruction* myDC)
  : myDetector(myDC)
{
  #include "Bina_Detector.cfg"
  generator_min=myDC->GetKinematicsMin();
  generator_max=myDC->GetKinematicsMax();
  npd_choice = myDC->GetNpdChoice(); 
  icros=myDC->GetNeumann();
  bfwhmx = myDC ->GetBfwhmX();    
  bfwhmy = myDC->GetBfwhmY();   
  bt = myDC->GetBtEnergy();
  pz = myDC->GetPz() ;
  pzz = myDC->GetPzz() ;
  themin = myDC->GetThetaMin()*180/M_PI; 
  themax = myDC->GetThetaMax()*180/M_PI; 
  themin2 = myDC->GetTheta2Min()*180/M_PI ; 
  themax2 = myDC->GetTheta2Max()*180/M_PI ;
  fimin = myDC->GetPhiMin()*180/M_PI ;
  fimax = myDC->GetPhiMax()*180/M_PI ; G4cout<<"\n";
  thigh = myDC->GetTargetHigh();
  tXplace = myDC->GetTargetXplace();
  tYplace = myDC->GetTargetYplace();
  tZplace = myDC->GetTargetZplace();
  int n_particle = 1;
  bt1 = 0.;		bt2 = 0.;		bt3 = 0.;
  t1  = 0.;		t2  = 0.;		t3  = 0.;
  fi1 = 0.;             fi2 = 0.;		fi3 = 0.;
  bt /=1000.;

  RandomInit();

//  set particles mass
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle = particleTable->FindParticle("deuteron");
  d_mass = particle->GetPDGMass()/1000.;

  particle = particleTable->FindParticle("proton");
  p_mass = particle->GetPDGMass()/1000.;

  particle = particleTable->FindParticle("neutron");
  n_mass = particle->GetPDGMass()/1000.;

  //*********************** Uwaga dokladnosc masy to 6 miejsc po przecinku *GeV ***********

  if (npd_choice <  0)
  {
    particleGun1 = new G4ParticleGun(n_particle);
    if (npd_choice == -1) particle = particleTable->FindParticle("proton");
    else if (npd_choice == -2) particle = particleTable->FindParticle("deuteron");
    else if (npd_choice==-3) particle=particleTable->FindParticle("neutron");
    else
    {
      G4cout << "Error in configuration file geo.mac !!! Bad npd_choice = "<<npd_choice<<G4endl;
      exit(1);
    }
    particleGun1->SetParticleDefinition(particle);
  }
  else if ((npd_choice == 0)||(npd_choice == 1))
  {
    if (npd_choice == 1) ugelast_read();  //read cross table
    particleGun1 = new G4ParticleGun(n_particle);
    particleGun2 = new G4ParticleGun(n_particle);

    particle = particleTable->FindParticle("proton");
    particleGun1->SetParticleDefinition(particle);

    particle = particleTable->FindParticle("deuteron");
    particleGun2->SetParticleDefinition(particle);
  }
  else if (npd_choice == 2)
  {
    //filell.open("fileoutput.txt",ios::out|ios::app);
    break_read();			//read cross table and analysing power
    particleGun1 = new G4ParticleGun(n_particle);
    particleGun2 = new G4ParticleGun(n_particle);
    particleGun3 = new G4ParticleGun(n_particle);

    particle = particleTable->FindParticle("proton");
    particleGun1->SetParticleDefinition(particle);

    G4ParticleDefinition* particle2 = particleTable->FindParticle("proton");
    particleGun2->SetParticleDefinition(particle2);

    particle = particleTable->FindParticle("neutron");
    particleGun3->SetParticleDefinition(particle);

  }
  else
  {
    G4cout << "Error in configuration file geo.mac !!! Bad npd_choice = "<<npd_choice<<G4endl;
    exit(1);
  }

}

Bina_PrimaryGeneratorAction::~Bina_PrimaryGeneratorAction()
{
  delete particleGun1;
  if ((npd_choice == 0)||(npd_choice == 1)) delete particleGun2;
  if (npd_choice == 2 ) delete particleGun3;
}

void Bina_PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  GetChoice(npd_choice);

  if (npd_choice <  0)
  {
   // G4cout <<"Npd_choice : "<<npd_choice<<G4endl;
    elastic(momentum);
    Pos();		//random vertex
    GetStartEnergy(bt1*1000.);
    GetStartAngleTheta(t1);
    GetStartAnglePhi(fi1);
    GetStartPosition(vertex[0],vertex[1],vertex[2]);

    particleGun1->SetParticleMomentumDirection(G4ThreeVector(momentum[0],momentum[1],momentum[2]));

    particleGun1->SetParticlePosition(G4ThreeVector(vertex[0],vertex[1],vertex[2]));

    particleGun1->SetParticleEnergy(bt*CLHEP::GeV);

    particleGun1->GeneratePrimaryVertex(anEvent);
  }
  else if (npd_choice == 0)
  {
    //G4cout <<"Npd_choice : "<<npd_choice<<G4endl;
    upunif(momentum);
    ProcNb(npd_choice);

    Pos(); 		//random vertex
    // rs GetStartEnergy(bt1*1000.,bt1*1000.);		// this isn't error !!!
    GetStartEnergy(bt*1000.,bt*1000.);
    GetStartAngleTheta(t1,t2);
    GetStartAnglePhi(fi1,fi2);

    GetStartPosition(vertex[0],vertex[1],vertex[2]);

    particleGun1->SetParticleMomentumDirection(G4ThreeVector(momentum[0],momentum[1],momentum[2]));
    particleGun2->SetParticleMomentumDirection(G4ThreeVector(momentum[3],momentum[4],momentum[5]));

    particleGun1->SetParticleEnergy(bt*CLHEP::GeV);
    particleGun2->SetParticleEnergy(bt*CLHEP::GeV);

    particleGun1->SetParticlePosition(G4ThreeVector(vertex[0],vertex[1],vertex[2]));
    particleGun2->SetParticlePosition(G4ThreeVector(vertex[0],vertex[1],vertex[2]));

    particleGun1->GeneratePrimaryVertex(anEvent);
    particleGun2->GeneratePrimaryVertex(anEvent);
  }
  else if (npd_choice == 1)
  {
    //G4cout <<"Npd_choice : "<<npd_choice<<G4endl;
    ugelast(momentum);
    ProcNb(npd_choice);
    Pos();

    GetStartEnergy(bt1*1000.,bt2*1000.);
    GetStartAngleTheta(t1,t2);
    GetStartAnglePhi(fi1,fi2);
    GetStartPosition(vertex[0],vertex[1],vertex[2]);

    particleGun1->SetParticleMomentumDirection(G4ThreeVector(momentum[0],momentum[1],momentum[2]));
    particleGun2->SetParticleMomentumDirection(G4ThreeVector(momentum[3],momentum[4],momentum[5]));

    particleGun1->SetParticleEnergy(bt1*CLHEP::GeV);
    particleGun2->SetParticleEnergy(bt2*CLHEP::GeV);

    particleGun1->SetParticlePosition(G4ThreeVector(vertex[0],vertex[1],vertex[2]));
    particleGun2->SetParticlePosition(G4ThreeVector(vertex[0],vertex[1],vertex[2]));

    particleGun1->GeneratePrimaryVertex(anEvent);
    particleGun2->GeneratePrimaryVertex(anEvent);

  }
  else //(npd_choice == 2)
  {
    //filell.open("fileoutput.txt",ios::out|ios::app);
    //G4cout <<"Npd_choice : "<<npd_choice<<G4endl;
    ugbreak(momentum);
    Pos();

    GetStartEnergy(bt1*1000.,bt2*1000.,bt3*1000.);
    GetStartAngleTheta(t1,t2,t3);
    GetStartAnglePhi(fi1,fi2,fi3);
    GetStartPosition(vertex[0],vertex[1],vertex[2]);

    particleGun1->SetParticleMomentumDirection(G4ThreeVector(momentum[0],momentum[1],momentum[2]));
    particleGun2->SetParticleMomentumDirection(G4ThreeVector(momentum[3],momentum[4],momentum[5]));
    particleGun3->SetParticleMomentumDirection(G4ThreeVector(momentum[6],momentum[7],momentum[8]));

    particleGun1->SetParticleEnergy(bt1*CLHEP::GeV);
    particleGun2->SetParticleEnergy(bt2*CLHEP::GeV);
    particleGun3->SetParticleEnergy(bt3*CLHEP::GeV);

    particleGun1->SetParticlePosition(G4ThreeVector(vertex[0],vertex[1],vertex[2]));
    particleGun2->SetParticlePosition(G4ThreeVector(vertex[0],vertex[1],vertex[2]));
    particleGun3->SetParticlePosition(G4ThreeVector(vertex[0],vertex[1],vertex[2]));

    particleGun1->GeneratePrimaryVertex(anEvent);
    particleGun2->GeneratePrimaryVertex(anEvent);
    particleGun3->GeneratePrimaryVertex(anEvent);
//    chooser++;
    //filell.close();
  }
}

void Bina_PrimaryGeneratorAction::RandomInit(int level )
{
	CLHEP::Ranlux64Engine *Engine1 = new CLHEP::Ranlux64Engine(aj1bx,level);
	CLHEP::Ranlux64Engine *Engine2 = new CLHEP::Ranlux64Engine(aj1by,level);

  CLHEP::RandFlat *GD1 = new CLHEP::RandFlat(*Engine1);
  CLHEP::RandFlat *GD2 = new CLHEP::RandFlat(*Engine2);
GaussDist[0] = GD1;
GaussDist[1] = GD2;

  

 

CLHEP::Ranlux64Engine *Engine3 = new CLHEP::Ranlux64Engine(aj1phi,level);  //init generator with seed
CLHEP::Ranlux64Engine *Engine4 = new CLHEP::Ranlux64Engine(aj1theta,level);
CLHEP::Ranlux64Engine *Engine5 = new CLHEP::Ranlux64Engine(ajcs1,level);
CLHEP::Ranlux64Engine *Engine6 = new CLHEP::Ranlux64Engine(aju,level);
CLHEP::Ranlux64Engine *Engine7 = new CLHEP::Ranlux64Engine(aj2theta,level);
CLHEP::Ranlux64Engine *Engine8 = new CLHEP::Ranlux64Engine(aj1ekin,level);
CLHEP::Ranlux64Engine *Engine9 = new CLHEP::Ranlux64Engine(aj2phi,level);
CLHEP::Ranlux64Engine *Engine10 = new CLHEP::Ranlux64Engine(ajcs1,level);
CLHEP::Ranlux64Engine *Engine11 = new CLHEP::Ranlux64Engine(ajbz,level);

CLHEP::RandFlat *GD3 = new CLHEP::RandFlat(*Engine3);
CLHEP::RandFlat *GD4 = new CLHEP::RandFlat(*Engine4);
CLHEP::RandFlat *GD5 = new CLHEP::RandFlat(*Engine5);
CLHEP::RandFlat *GD6 = new CLHEP::RandFlat(*Engine6);
CLHEP::RandFlat *GD7 = new CLHEP::RandFlat(*Engine7);
CLHEP::RandFlat *GD8 = new CLHEP::RandFlat(*Engine8);
CLHEP::RandFlat *GD9 = new CLHEP::RandFlat(*Engine9);
CLHEP::RandFlat *GD10 = new CLHEP::RandFlat(*Engine10);
CLHEP::RandFlat *GD11 = new CLHEP::RandFlat(*Engine11);

  FlatDist[0] = GD3;
  FlatDist[1] = GD4;
  FlatDist[2] = GD5;
  FlatDist[3] = GD6;
  FlatDist[4] = GD7;
  FlatDist[5] = GD8;
  FlatDist[6] = GD9;
  FlatDist[7] = GD10;
  FlatDist[8] = GD11;
}

double Bina_PrimaryGeneratorAction::RandomGauss(double seed, double mean , double deviation )
{
  double num;
  if (seed == aj1bx)
  {
    num = (*GaussDist[0]).fire(mean, deviation);
    (*GaussDist[0]).fire(mean, deviation);
    //num = (*FlatDist[0]).fire(mean, deviation);
    //(*FlatDist[0]).fire(mean, deviation);    
   return num;
  }
  if (seed == aj1by)
  {
    num = (*GaussDist[1]).fire(mean, deviation);
    (*GaussDist[1]).fire(mean, deviation);
    //num = (*GD2).fire(mean, deviation);
    //(*GD2).fire(mean, deviation);
    return num;
  }
  G4cout <<"Error in choice random gauss distribution!! Seed = "<<seed<<G4endl;
  exit (1);
}

double Bina_PrimaryGeneratorAction::RandomFlat(double seed, double m , double n  )
{
  if (seed == aj1phi)   return (*FlatDist[0]).fire(m, n);
  if (seed == aj1theta) return (*FlatDist[1]).fire(m, n);// metod fire(m,n) - > return double ]m,n[
  if (seed == ajcs1) 	return (*FlatDist[2]).fire(m, n);
  if (seed == aju) 	return (*FlatDist[3]).fire(m, n);
  if (seed == aj2theta) return (*FlatDist[4]).fire(m, n);
  if (seed == aj1ekin) 	return (*FlatDist[5]).fire(m, n);
  if (seed == aj2phi) 	return (*FlatDist[6]).fire(m, n);
  if (seed == ajcs1) 	return (*FlatDist[7]).fire(m, n);
  if (seed == ajbz) 	return (*FlatDist[8]).fire(m, n);

  G4cout <<"Error in choice random : "<<m<<"-"<<n<<" !! Seed = "<<seed<<G4endl;
  exit (1);
}

void Bina_PrimaryGeneratorAction::Pos(void)
{
  double bsgx,bsgy;

  bsgx = bfwhmx/(2.*sqrt(2.*log(2.)));
  bsgy = bfwhmy/(2.*sqrt(2.*log(2.)));
  // tXplace,  tYplace,  tZplace, thigh read from geo.mac in cm and 
  // recalculated by Geant to mm; bfwhmx, bfwhmy read in default mm. 
  vertex[0] = tXplace + RandomGauss(aj1bx,0, bsgx); 
  vertex[1] = tYplace + RandomGauss(aj1by,0, bsgy);
  vertex[2] = tZplace - thigh + RandomFlat(ajbz)*2*thigh;
 }

double* Bina_PrimaryGeneratorAction::elastic(double *ptot)
{
  double pm, thcos, pphi, bp1, bpproj;	  // auxillary variable
  double thcosm,thcos0;

  thcosm = cos(themax*CLHEP::pi/180.);
  thcos0 = cos(themin*CLHEP::pi/180.);
  thcos = thcosm + (thcos0 - thcosm)*RandomFlat(aj1theta);
  t1 = 180./CLHEP::pi*acos(thcos);
  pphi = CLHEP::pi*(2.* RandomFlat(aj1phi));
  fi1=pphi*180./CLHEP::pi;

  if   (npd_choice == -1) pm = p_mass;  		  //proton
  else if(npd_choice==-2) pm = d_mass;  	            //(npd_choice == -2)  //deuteron
  else pm=n_mass;
  bp1 = sqrt(bt*(bt + 2*pm));
  ptot[2] = bp1*thcos;  			//ptot1(3)
  bpproj = sqrt(bp1*bp1 - ptot[2]*ptot[2]);
  ptot[0] = bpproj*cos(pphi);			//ptot1(1)
  ptot[1] = bpproj*sin(pphi);			//ptot1(2)
  bt1 = bt;
  //rs  t1 = 180./pi*themin;
  return ptot;
}

double* Bina_PrimaryGeneratorAction::upunif(double *ptot)
{
  double ptot1[3],ptot2[3],bp1=0.,bp2=0.,bpproj;	// auxillary variable
  double thcosm,thcos0;
  //,thcosm2,thcos02,
  double thcos,thcos2,pphi; // auxillary variable
  double *w_bp1,*w_bp2;

  w_bp1 = &bp1;
  w_bp2 = &bp2;

  double pm = p_mass;
  double dm = d_mass;
  int i;
  for (i=0;i<3;i++)
  {
    ptot1[i] = ptot[i];
    ptot2[i] = ptot[i+3];
  }
  thcosm = cos(themax*CLHEP::pi/180.);
  thcos0 = cos(themin*CLHEP::pi/180.);
//  thcosm2 = cos(themax2*CLHEP::pi/180.);
//  thcos02 = cos(themin2*CLHEP::pi/180.);

  thcos = thcosm + (thcos0 - thcosm)*RandomFlat(aj1theta);
  //rs pphi = pi*(2.*RandomFlat(aj1phi)-1.);
   pphi = CLHEP::pi*(2.*RandomFlat(aj1phi));
  //G4cout <<"thcos : "<<thcos<<G4endl;
  thcos2 = gelkin(thcos,w_bp1,w_bp2);	//elastic kinematics

  bt1 = sqrt(pm*pm + bp1*bp1) - pm;
  bt2 = sqrt(dm*dm + bp2*bp2) - dm;

  ptot1[2] = bp1*thcos;
  bpproj = sqrt(bp1*bp1 - ptot1[2]*ptot1[2]);
  ptot1[0] = bpproj*cos(pphi);
  ptot1[1] = bpproj*sin(pphi);

  ptot2[2] = bp2*thcos2;
  bpproj = sqrt(bp2*bp2 - ptot2[2]*ptot2[2]);
  ptot2[0] = bpproj*cos(CLHEP::pi + pphi);
  ptot2[1] = bpproj*sin(CLHEP::pi + pphi);

  t1 = 180./CLHEP::pi*acos(thcos);
  t2 = 180./CLHEP::pi*acos(thcos2);

  //rs
  fi1 =   pphi*180/CLHEP::pi;
  fi2 = fi1+180.;
  if(fi2>360)fi2=fi2-360;



  for (i=0;i<3;i++)
  {
    ptot[i] = ptot1[i];
    ptot[i+3] = ptot2[i];
  }
  return ptot;
}

double* Bina_PrimaryGeneratorAction::ugelast(double* ptot)
{
// elastic scattering - random generation of proton emmision angle
// according to cross section and analysing power
//
// gelkin for calculation of deuteron emmision angle and of particle energies

  double thcosm,thcos0,/*thcosm2,thcos02,*/bp1,bp2,thelab,thcos,thsin;
  double difft=0.,difft1=0.,csi,cstest,it11,t20,t22,pphi=0.,/*czvec,cztens1,cztens2,*/fphi=0.;
  double bt1_temp=0.,bt2_temp=0.,thcos2=0.,ptot1[3],ptot2[3],bpproj,u,pm,dm;
//  double csmax = 15.;
  int ithet/*,it1*/,i,ok;

  double *w_bp1,*w_bp2;

  t1 = 0.;
  t2 = 0.;

  w_bp1 = &bp1;
  w_bp2 = &bp2;
  pm = p_mass;
  dm = d_mass;

  thcosm = cos(themax*CLHEP::pi/180.);
  thcos0 = cos(themin*CLHEP::pi/180.);
//  thcosm2 = cos(themax2*CLHEP::pi/180.);
//  thcos02 = cos(themin2*CLHEP::pi/180.);

  do
  {
    ok = 1;
    thelab = CLHEP::pi*0.5*RandomFlat(aj1theta);
    thcos = cos(thelab);

    if((thcos < thcosm) || (thcos > thcos0)) {ok = 0; continue;}
    thsin = sin(thelab);
    thelab *= 180./CLHEP::pi;
    ithet = 0;
    for(i=0;i<70;i++)
    {
      if(thelab <= thl[i] && thelab > thl[i+1])
      {
        ithet = i;
        difft = thl[i+1] - thl[i];
        difft1 = thelab - thl[i];
        break;
      }
    }
    if(fabs(difft) < 0.0001) difft = 0.0001;

// linear interpolation of unpolarized cross section
    csi  = ds[ithet] + (ds[ithet+1] - ds[ithet])/difft*difft1;
    csi  = csi * thsin;
    cstest=csmax*RandomFlat(ajcs1);
    if (csi < cstest) {ok =0; continue;};
    if(fabs(pz) >= 0.0001 || fabs(pzz) >= 0.0001)
    {

// linear interpolation of analysing powers
      it11 = ap1[ithet] + (ap1[ithet+1] - ap1[ithet])/difft*difft1;
      t20  = ap2[ithet] + (ap2[ithet+1] - ap2[ithet])/difft*difft1;
      t22  = ap3[ithet] + (ap3[ithet+1] - ap3[ithet])/difft*difft1;

// von Neumann's method applied to ph
      do
      {
        pphi = CLHEP::pi*(2.*RandomFlat(aj1phi) - 1.);
        u = 2.*RandomFlat(aju);
//        czvec = sqrt(3.)*pz*it11*cos(pphi);
//        cztens1 = 1./sqrt(8.)*pzz*t20;
//        cztens2 = - sqrt(3.)/2.*pzz*t22*cos(2.*pphi);
        fphi = 1. + sqrt(3.)*pz*it11*cos(pphi) + 1./sqrt(8.)*pzz*t20 - sqrt(3.)/2.*pzz*t22*cos(2.*pphi);
      }
      while(u > fphi);
    }
    else 
    {
      pphi = CLHEP::pi*(2.*RandomFlat(aj1phi) - 1.);
      fphi = CLHEP::pi+pphi;
    }

    bp1 = 0.;
    bp2 = 0.;
    thcos2 = gelkin(thcos,w_bp1,w_bp2);		//elastic kinematics

    bt1_temp = sqrt(pm*pm + bp1*bp1) - pm;
    bt2_temp = sqrt(dm*dm + bp2*bp2) - dm;

  }
  while(/*((thcos2 > thcos02) || (thcos2 < thcosm2))||*/(ok == 0)); //removed limits on deuteron th angle


  bt1 = bt1_temp;	//G4cout <<"Energy bt1 : "<<bt1*GeV<<G4endl;
  bt2 = bt2_temp;	//G4cout <<"Energy bt2 : "<<bt2*GeV<<G4endl;

  ptot1[2] = bp1*thcos;
  bpproj = sqrt(bp1*bp1 - ptot1[2]*ptot1[2]);
  ptot1[0] = bpproj*cos(pphi);
  ptot1[1] = bpproj*sin(pphi);
  ptot2[2] = bp2*thcos2;
  bpproj = sqrt(bp2*bp2 - ptot2[2]*ptot2[2]);
  ptot2[0] = bpproj*cos(CLHEP::pi+pphi);
  ptot2[1] = bpproj*sin(CLHEP::pi+pphi);

  t1 = 180./CLHEP::pi*acos(thcos);
  t2 = 180./CLHEP::pi*acos(thcos2);
  fi1 = pphi*180./CLHEP::pi;
  fi2 = fphi*180./CLHEP::pi;

//  it1 = (int)thelab + 1;

  for (i=0;i<3;i++)
  {
    ptot[i] = ptot1[i];
    ptot[i+3] = ptot2[i];
  }

  return ptot;
}

double* Bina_PrimaryGeneratorAction::ugbreak(double* ptot)
{
// random choice of starting momenta of p, p and n, according to the
// d(p,pp)n break-up cross section
// step(i) - steps of the cross section tables
// (respectively: theta1, theta2, phi12, e1)
// offs(i) - offsets of the cross section tables
  double ptot1[3],ptot2[3],ptot3[3],x0[4],bp1,bp2,bp3,pm,dm,etot,/*e1min,*/thcosm,thcos0;
  double thcosm2,thcos02,thcos1,thcos2,theta1,theta2,e1,e1_0=0.,phi12,phi1=0.,cstest,csi;
  double axi,ayi,axxi,axyi,ayyi,thcos3,phi2,phi3,phi1r,phi2r,phi3r,bpproj;
  double ya[4][4][4][4], yax[4][4][4][4];
  double yay[4][4][4][4],yaxx[4][4][4][4],yaxy[4][4][4][4],yayy[4][4][4][4];
  double e_random;
  int ind[4],ii,jj,kk,ll,i,j,k,l;
  int signum;
  double paramS,deltaS,deltaE1,deltaE2,s,s0,e1_1,e2_1;
//  double csmax = 10.;
//  bool switcher;
  bool quest;
  double e2_2,e1_2;
  int icn = 0;
  double icnn = 0.;
  double phi13=0.,e2_0=0.,e3_0=0.,bp_temp[3]={0.,0.,0.};
  double *w_phi13,*w_e1_0,*w_e3_0,*bp;

  w_phi13 = &phi13;
  w_e1_0 = &e1_0;
  w_e3_0 = &e3_0;
  bp = bp_temp;

//  double *w_bp1,*w_bp2,*w_bp3;

//  w_bp1 = &bp1;
//  w_bp2 = &bp2;
//  w_bp3 = &bp3;

  pm = p_mass;
  dm = d_mass;
//  int icros = 1;//zmienione !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//G4cout<<icros<<" CROSS SECTIONS AND NEUMANN METHOD \n";
  etot = bt + dm + pm;
//  e1min = pm+0.005;// + 0.005;   // 5 MeV threshold usuniÄty na prÃ³bÄ co bÄdzie!!!!!!!!!!!!!!!!!!!1
  thcosm = cos(themax*CLHEP::pi/180.);
  thcos0 = cos(themin*CLHEP::pi/180.);
  thcosm2 = cos(themax2*CLHEP::pi/180.);
  thcos02 = cos(themin2*CLHEP::pi/180.);

 int nlop = 0;

/***************************************************TTEESSTT
  for (l=0;l<4;l++)
    {
      for(k=0;k<14;k++)
      {
        for(j=0;j<13;j++)
        {
          for (i=55;i<75;i++)
         {filell<<sig[l][k][j][i]<<G4endl;}}}}
**********************************************************************/
/*double test1a=30.,test2a,test3a,test4a;
double* test1=&test1a;
double* test2=&test2a;
double* test3=&test3a;
double* test4=&test4a;
int test5;

for (i=0;i<180;i++) {
  for (j=0;j<300;j++) {
    test5=dpb_kin3(double(i),test1,j/1000.,test2,test3,test4);
    3<<test5<<' '<<i<<' '<<j<<' '<<test1a<<' '<<test2a<<' '<<test3a<<' '<<test4a<<G4endl;
//(double phi12,double *wphi13,double e1_0,double* we2_0,double* we3_0,double *bp_t)
    }
  }
*/
  do
  {
    nlop++;
    // commented out for test4:

    thcos1 = thcosm + (thcos0 - thcosm)*RandomFlat(aj1theta);
    thcos2 = thcosm2 + (thcos02 - thcosm2)*RandomFlat(aj2theta);
    theta1 = 180./CLHEP::pi*acos(thcos1);
    theta2 = 180./CLHEP::pi*acos(thcos2);

    //if (icros == 0)
//    {
//      e1 = e1min + (etot - e1min)*RandomFlat(aj1ekin);
//      e1_0 = e1 - pm;
//    }
// ************************************************************************//
//*************************************************************************//

//*** Proba wprowadzenia przekroju czynnego do losowania z fi12 ***/rs
    phi12 = fimin + (fimax - fimin)*RandomFlat(aj2phi);   // old   18.03.2005

      e_random=RandomFlat(aj1ekin);//zugefuegt !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!30III10
      s0 = generator_min+e_random*(generator_max-generator_min);
   /////////   e1_1=0.005+(bt-0.005)*e_random; //extra IV2010
 //G4cout<<bt<<' '<<e1min<<' '<<e1_1<<' '<<e1<<G4endl;
 //filell<<"Etot "<<etot<<" E1min "<<e1min<<" diff "<<etot-e1min<<G4endl;
  //    e1_0 = e1 - pm; deleted 26 april 2010
  //  thcos1 = cos(theta1*pi/180.);
  //  thcos2 = cos(theta2*pi/180.);
    //    G4cout << " theta1= " << theta1  << " theta2= " << theta2 <<
    //	      " phi12 " << phi12 <<G4endl;


// *****-----------------------------------------------------------*****//

// for unpolarized c.s. phi1 can be generated uniformly
    // rs if ((icros == 0) || ((pzz == 0) && (pz == 0))) phi1=-180+360*RandomFlat(aj1phi);

    if ((icros == 0) || ((pzz == 0) && (pz == 0))) phi1=360*RandomFlat(aj1phi);

    t1 = theta1;
    t2 = theta2;
    if (icros != 0)    // uniformly generated protons
    {
// find indexes of the cross section tables for interpolation procedures

 // for(t1=11.;t1<26.;t1+=0.5) {
  //  for (t2=11.;t2<26.;t2+=.5) {
   //   for (phi12=0.;phi12<=180.;phi12+=10.) {
     //   for (s0=0.;s0<.25;s0+=0.01) {
      x0[0] = t1;
      x0[1] = t2;
      x0[2] = phi12;
      x0[3] = s0;

//      filell<<"theta1 theta2 phi12 e1_0 niedergeschrievven\n";
      //filell<<"=== "<<x0[0]<<' '<<x0[1]<<' '<<x0[2]<<' '<<x0[3]<<G4endl;
      for (i=0;i<4;i++)
      {
        ind[i] = max(int((x0[i]-offs[i])/step[i]),1);
//filell<<"ind  "<<i<<' '<<ind[i]<<' '<<x0[i]<<G4endl;
        for (j=0;j<4;j++)
        {
          if(i == 0) x1a[j] = offs[i] + (ind[i] + j-1)*step[i];
       //   filell<<t1<<' '<<x1a[j]<<' '<<x0[0]<<' '<<ind[i]<<G4endl;}
          if(i == 1) x2a[j] = offs[i] + (ind[i] + j-1)*step[i];
       //   filell<<t2<<' '<<x2a[j]<<' '<<x0[1]<<' '<<ind[i]<<G4endl;}
          if(i == 2) x3a[j] = offs[i] + (ind[i] + j-1)*step[i];
        //  filell<<phi12<<' '<<x3a[j]<<' '<<x0[2]<<' '<<ind[i]<<G4endl;}
          if(i == 3) x4a[j] = offs[i] + (ind[i] + j)*step[i];
        //  filell<<s0<<' '<<x4a[j]<<' '<<x0[3]<<' '<<ind[i]<<G4endl;}
        }
      }
//for (i=0;i<200;i++) {
//filell<<sig[7][7][1][i]<<G4endl;
//}

// interpolations
//   - of 'unpolarized' cross section

      for(i=0;i<4;i++)
      {
        for(j=0;j<4;j++)
        {
          for(k=0;k<4;k++)
          {
            for(l=0;l<4;l++)
	    {
	      ii = ind[0] + i-1;
              jj = ind[1] + j-1;
              kk = ind[2] + k-1;
                if (kk>=13) kk=24-kk;//up to 12 and down 
              ll = ind[3] + l-1;
              ya[i][j][k][l] = sig[ii][jj][kk][ll];
//if (i==2&&j==2&&k==2&&l==2) filell<<x0[0]<<' '<<x0[1]<<' '<<x0[2]<<' '<<x0[3]<<' '<<ii<<' '<<jj<<' '<<kk<<' '<<ll<<' '<<sig[ii][jj][kk][ll]<<' '<<ya[i][j][k][l]<<G4endl;
	    }
          }
        }
      }
      //change matrix from [0..n-1] to [1..n]
      for (i=0;i<4;i++)
      {
        x1a_r[i+1] = x1a[i];
        x2a_r[i+1] = x2a[i];
        x3a_r[i+1] = x3a[i];
        x4a_r[i+1] = x4a[i];
      }

      for (i=0;i<4;i++)
        for (j=0;j<4;j++)
          for (k=0;k<4;k++)
            for (l=0;l<4;l++) {ya_r[i+1][j+1][k+1][l+1]=ya[i][j][k][l];
//filell<<i<<' '<<j<<' '<<k<<' '<<l<<' '<<ya_r[i+1][j+1][k+1][l+1]<<G4endl;
            }
      csi = rinterp(ya_r,x0[0],x0[1],x0[2],x0[3]);
      if (csilocmax<csi) {csilocmax=csi;
                          toSix++;
                          if (toSix>5) csmax=csilocmax;}
      cstest = csmax*RandomFlat(ajcs1);
//filell<<t1<<' '<<t2<<' '<<phi12<<' '<<s0<<' '<<csi<<' '<<cstest<<G4endl;
//G4cout<<x0[0]<<' '<<x0[1]<<' '<<x0[2]<<' '<<x0[3]<<' '<<ii<<' '<<jj<<' '<<kk<<' '<<ll<<' '<<csi<<
//' '<<cstest<<G4endl;
//filell<<csi<<' '<<s0<<G4endl;
//    if(csi<0.)  filell<<"----- "<<csi<<' '<<x0[0]<<' '<<x0[1]<<' '<<x0[2]<<' '<<x0[3]<<G4endl;
  /*    for(i=0;i<4;i++)
      
        for(j=0;j<4;j++)
        
          for(k=0;k<4;k++)
          
            for(l=0;l<4;l++)
	    filell<<i<<' '<<j<<' '<<k<<' '<<l<<' '<<ya_r[i+1][j+1][k+1][l+1]<<G4endl;
             

}*/
   //  filell<<t1<<' '<<t2<<' '<<phi12<<' '<<s0<<' '<<csi<<' '<<csmax<<' '<<cstest<<G4endl;
//}}}}
 //     G4cout << "csi = "<<csi<<"\t cstest = "<<cstest<<"\t phi12 = "<<phi12<<G4endl;
      if (csi < cstest) continue;
//if (csi>csmax*0.5) G4cout<<phi12<<' '<<s0<<G4endl;
//   - of analysing powers
      if((pzz != 0)&&(pzz != 0))
      {
        for(i=0;i<4;i++)
        {
          for(j=0;j<4;j++)
          {
            for(k=0;k<4;k++)
            {
              for(l=0;l<4;l++)
	      {
	        ii = ind[0] + i-1;
                jj = ind[1] + j-1;
                kk = ind[2] + k-1;
                ll = ind[3] + l-1;
                yax[i][j][k][l]  = axn[ii][jj][kk][ll];
                yay[i][j][k][l]  = ayn[ii][jj][kk][ll];
                yaxy[i][j][k][l] = axy[ii][jj][kk][ll];
                yaxx[i][j][k][l] = axx[ii][jj][kk][ll];
                yayy[i][j][k][l]  = ayy[ii][jj][kk][ll];
                yax_r[i+1][j+1][k+1][l+1]=yax[i][j][k][l];
   	        yay_r[i+1][j+1][k+1][l+1]=yay[i][j][k][l];
	        yaxy_r[i+1][j+1][k+1][l+1]=yaxy[i][j][k][l];
	        yaxx_r[i+1][j+1][k+1][l+1]=yaxx[i][j][k][l];
	        yayy_r[i+1][j+1][k+1][l+1]=yayy[i][j][k][l];

	      }
            }
          }
        }
        axi = rinterp(yax_r,x0[0],x0[1],x0[2],x0[3]);
        ayi = rinterp(yay_r,x0[0],x0[1],x0[2],x0[3]);
        axyi = rinterp(yaxy_r,x0[0],x0[1],x0[2],x0[3]);
        axxi = rinterp(yaxx_r,x0[0],x0[1],x0[2],x0[3]);
        ayyi = rinterp(yayy_r,x0[0],x0[1],x0[2],x0[3]);

        phi1 = cspol(axi,ayi,axxi,axyi,ayyi);
      }
    }
    // *************************************************//
    // nowa procedura - gdy brak zgodnosci z kinematyka //
    // losowana jest ponownie energia a nie katy        //
    // *************************************************//

    G4int itry = 0;
//    cout<<"przed: nlop "<<nlop<<"  theta1 "<<theta1<<"  theta2 "<<theta2<<"  phi12 "<<phi12<<endl;
    icnn  = 0;
//    do {
//      itry += 1;
//      e1 = e1min + (etot - e1min)*RandomFlat(aj1ekin);
//      e1_0 = e1 - pm;
//      icnn = dpb_kin3(phi12,w_phi13,e1_0,w_e2_0,w_e3_0,bp);
//    }  while(icnn != 1 && itry<200);
//    do {
//      itry += 1;

//      e1_0 = e1 - pm;


//filell<<phi12<<' '<<t1<<' '<<t2<<' '<<s0<<G4endl;
paramS=0.000001;
//for (e1_0=0.;e1_0<.13;e1_0+=0.000001) {
//      icnn = dpb_kin3(phi12,w_phi13,e1_0,w_e2_0,w_e3_0,bp);
//      filell<<e1_0<<' '<<e2_0<<' '<<e3_0<<' '<<icnn<<G4endl;
//      }
e2_0=0.;
do {
      icnn = dpb_kin3(1,phi12,w_phi13,e2_0,w_e1_0,w_e3_0,bp);
      e2_0+=paramS;
      } while (icnn<0.000001);
do {
      icnn = dpb_kin3(1,phi12,w_phi13,e2_0,w_e1_0,w_e3_0,bp);
      e2_0-=paramS*0.00001;
      } while (icnn>0.);


if (icros==1) {
//filell<<"1 "<<e2_0<<' '<<e1_0<<G4endl;
s=.0;//verschiebung
deltaS=0.;
deltaE1=0.;
e2_2=e2_0;
e1_2=e1_0;
signum=1;
 //     filell<<e1_0<<' '<<e2_0<<' '<<s<<' '<<icnn<<G4endl;
for (int switcher=0;switcher<=2;switcher++) {
if(switcher==2)
  {
  signum=-1;
  }
else {
  e2_0=e2_2;
  e1_0=e1_2;
}

if (s>=s0-0.0000005) {
  if (switcher==1) {
    e1_0=0;
    e2_0=0;
    e3_0=0;
  }
  break;
}
do {
//if (switcher==2) G4cout<<e1_0<<' '<<e2_0<<' '<<s<<' '<<icnn<<G4endl;
      e1_1=e1_0;
      e2_1=e2_0;
      e2_0+=signum*paramS;
    //  filell<<e1_0<<' ';
      dpb_kin3(switcher,phi12,w_phi13,e2_0,w_e1_0,w_e3_0,bp);
      deltaE1=fabs(e1_0-e1_1);
      deltaE2=fabs(e2_0-e2_1);
  //    filell<<"Deltae1,e2 erste "<<deltaE1<<' '<<deltaE2;
      deltaS=E(deltaE1,deltaE2);
  //    filell<<' '<<deltaS<<G4endl;
      paramS=0.000001*deltaE2/deltaS;
      
      e2_0=e2_1+signum*paramS;
    //  filell<<e1_0<<' '<<G4endl;
      icnn=dpb_kin3(switcher,phi12,w_phi13,e2_0,w_e1_0,w_e3_0,bp);
//filell<<icnn<<G4endl;
      deltaE1=fabs(e1_0-e1_1);
      deltaE2=fabs(e2_0-e2_1);
//      filell<<"Deltae1,e2 "<<deltaE1<<' '<<deltaE2;
      s+=E(deltaE1,deltaE2);
   //   filell<<' '<<s0<<' '<<switcher<<' ';
    //  filell<<e1_0<<' '<<e2_0<<' '<<s<<' ';
//G4cout<<s<<' '<<s0<<G4endl;
if (switcher==0) quest=fabs(e1_0-e2_0)<fabs(e1_1-e2_1);
if (switcher==1) quest=s<s0-0.0000005&&icnn>1e-10;
if (switcher==2) quest=s<s0-0.0000005;
//filell<<quest<<G4endl;
//if (switcher==2) G4cout<<"+++ "<<s<<' '<<s0<<"\n";
//if (quest==0) filell<<s<<' '<<s0<<' '<<s0-0.0000005<<' '<<icnn<<G4endl;
} while (quest);
//G4cout<<switcher<<' '<<e1_0<<' '<<e2_0<<G4endl;
}
//G4cout<<"Hier gibts ne Information "<<e1_0<<' '<<e2_0<<G4endl;
//G4cout<<s<<' '<<s0<<G4endl;
//G4cout<<"+++++++++++++++++++++++++++++++++++++++++\n";
//////////////////filell<<icnn<<' '<<phi12<<' '<<*w_phi13<<' '<<e1_0<<' '<<*w_e2_0<<' '<<*w_e3_0<<' '<<*bp<<G4endl;
 //   }  while(icnn != 1 && itry<200);
// break-up kinematics
    icn = int(bool(icnn));//dpb_kin3(phi12,w_phi13,e1_0,w_e2_0,w_e3_0,bp);
//filell<<icn<<' '<<phi12<<' '<<*w_phi13<<' '<<e1_0<<' '<<*w_e2_0<<' '<<*w_e3_0<<' '<<*bp<<' '<<"ICN"<<G4endl;
   //G4cout<<icn<<G4endl;
  }
else {
itry=0;
 do {
  itry+=1;
  e1 = generator_min + (etot - generator_max)*RandomFlat(aj1ekin);
//  G4cout<<e1<<' '<<generator_min<<' '<<etot<<'\n';
  e2_0 = e1 - pm;
  icnn = dpb_kin3(1,phi12,w_phi13,e2_0,w_e1_0,w_e3_0,bp);
  icn=int(bool(icnn));
//  G4cout<<e2_0<<' '<<*w_e1_0<<' '<<*w_e3_0<<' '<<icnn<<'\n';
  } while (icn!=1&&itry<200);
  G4cout<<itry<<'\n';
 }
} while(icn!=1);//!= 1);
// calculate momenta of p, p and n
  thcos3 = cos(t3*CLHEP::pi/180.);  //rs - change t1 -> t3
  bp1 = bp[0];
  bp2 = bp[1];
  bp3 = bp[2];

  bt1 = e1_0;
  bt2 = e2_0;
  bt3 = e3_0;

  phi2 = phi1 + phi12;
  phi3 = phi1 + phi13;

  if (phi2 >= 360) phi2 = phi2 - 360;
  if(phi3 >= 360) phi3 = phi3 - 360;
  phi1r = CLHEP::pi/180.*phi1;
  phi2r = CLHEP::pi/180.*phi2;
  phi3r = CLHEP::pi/180.*phi3;

  fi1 = phi1;
  fi2 = phi2;
  fi3 = phi3;


  ptot1[2] = bp1*thcos1;
  bpproj = sqrt(bp1*bp1 - ptot1[2]*ptot1[2]);
  ptot1[0] = bpproj*cos(phi1r);
  ptot1[1] = bpproj*sin(phi1r);
  ptot2[2] = bp2*thcos2;
  bpproj = sqrt(bp2*bp2 - ptot2[2]*ptot2[2]);
  ptot2[0] = bpproj*cos(phi2r);
  ptot2[1] = bpproj*sin(phi2r);
  ptot3[2] = bp3*thcos3;
  bpproj = sqrt(bp3*bp3 - ptot3[2]*ptot3[2]);
  ptot3[0] = bpproj*cos(phi3r);
  ptot3[1] = bpproj*sin(phi3r);

  //t1 = thelab;
  //t2 = 180./pi*acos(thcos2);

  for (i=0;i<3;i++)
  {
    ptot[i] = ptot1[i];
    ptot[i+3] = ptot2[i];
    ptot[i+6] = ptot3[i];

  }



  return ptot;
}

double Bina_PrimaryGeneratorAction::gelkin(double thcos, double *wbp1, double *wbp2)
{  // output : bp1,bp2,thcos2
  double et,etp,beta,gamma,betsq,gamsq,rang3,sinth,sinsq; 	  	//auxillary variable
  double bp0,x1,x2,x3,x4,ep3c,ep4c,p3c,p4c,delta,a,b,c,test,coscm/*,thecm*/;  //auxillary variable
  double ep3l,ep4l,p3l,/*e3,e4,*/p4l,cosang4;  				  //auxillary variable
//  double bej, ang4;
  double thcos2;

  double dm = d_mass;
  double pm = p_mass;

  bp0 = sqrt(bt*(bt + 2.*dm));
  x1 = dm;
  x2 = pm;
  x3 = pm;
  x4 = dm;

  et = x1 + x2 + bt;
  etp = sqrt(x1*x1 + x2*x2 + 2.*x2*(x1 + bt));
  beta = bp0/et;
  gamma = et/etp;
  betsq = beta*beta;
  gamsq = gamma*gamma;
  rang3 = acos(thcos);
  sinth = sin(rang3);
  sinsq = sinth*sinth;

// Compute relativistic energy of outgoing particles
  ep3c = (etp*etp + x3*x3 - x4*x4)/(2.0*etp);
  ep4c = etp - ep3c;
  p3c = ep3c*ep3c - x3*x3;
  p4c = ep4c*ep4c - x4*x4;
  if (p4c <= 0.) p4c = 0.;
  else p4c = sqrt(p4c);
  if (p3c <= 0.) p3c = 0.;
  else p3c = sqrt(p3c);
  delta = beta*ep3c/p3c;
  a = betsq*gamsq*sinsq + 1.;
  b = delta*gamsq*sinsq;
  c = (gamsq*delta*delta+1.)*sinsq - 1.;
  test = b*b - a*c;
  if (test < -1e-6*b*b) return 0.;
  test = sqrt(max(test,0.));
  coscm = (-b+test)/a;

//  thecm = acos(coscm);
  ep3l = gamma*(ep3c + beta*p3c*coscm);
  ep4l = et - ep3l;
  p3l = sqrt(ep3l*ep3l - x3*x3);
//  e3 = ep3l - x3;
//  e4 = bt - e3;
//  bej = e3;
  p3l = sqrt(ep3l*ep3l - x3*x3);
  p4l = ep4l*ep4l - x4*x4;
  if (p4l <= 0.) p4l = 0.;
  else p4l = sqrt(p4l);
  *wbp1 = p3l;
  *wbp2 = p4l;
  if (p4l == 0.) cosang4 = 0.;
  else
  {
    cosang4 = gamma*(beta*ep4c - p3c*coscm)/p4l;
    if (fabs(cosang4) > 1.)
    {
      if (cosang4 > 1.) cosang4 = 1.;
      else cosang4 = - 1.;
    }
  }
//  ang4 = acos(cosang4) * 180./CLHEP::pi;
  thcos2 = cosang4;

  //G4cout<<"thcos : "<<ang4<<endl;
  return thcos2;
}

double Bina_PrimaryGeneratorAction::dpb_kin3(const int typ, double phi12,double *wphi13,double e2_0,double* we1_0,double* we3_0,double *bp_t)
{
//	THREE-BODY RELATIVISTIC KINEMATICS for d(p,pp)n reaction
//      at deuteron kinetic energy ekin = 0.130 GeV
//       input parameters:
//          phi12 - phi angle of the second proton with respect to the first one
//          e1_0  - kinetic energy of the first proton
//       output :
//           t3 - theta angle of neutron
//           phi13 - phi angle of neutron with respect to the first proton
//           e2_0  - kinetic energy of the second proton
//           e3_0  - kinetic energy of the neutron
//           bp    - matrix of momenta of proton, proton and neutron
//
//           angles in degrees
//  	     all energies in GeV

// P2(1) i P2(2) to dwa mozliwe rozwiazania na P2 dla ustalonego P1
// odpowiadaja temu dwie wartosci energii: E2(1) i E2(2)


  double r=0.,phi=0.,thta=0.,thta2_1=0.,etot=0.,phi2=0.,phi2_1=0.,phi13=0.,e1_0=0.,e3_0=0.,ff=0.,ph1r=0.,th1r=0.,ph2r=0.,th2r=0.,pol=0.,e1min=0.,e1=0./*,e1max=0.*/;
  // double p1=0.,p1x=0.,p1y=0.,p1z=0.,p2i3l=0.,p2i3ll=0.,e2i3=0.,a=0.,b=0.,c=0.,d=0.,delta=0.,e2p=0.,pp2=0.,e3=0.,pp3=0.,t3r=0.,cphi13=0.;
  double p1=0.,p1x=0.,p1y=0.,p1z=0.,p2i3l=0.,e2i3=0.,a=0.,b=0.,c=0.,d=0.,delta=0.,e2p=0.,pp2=0.,e3=0.,pp3=0.,t3r=0.,cphi13=0.;
  double uni2[3],p2i3[3],/*p2i3_1[3],*/uni2_1[3],p2[3],e2[3],pt2[3],pt3[3],phi13r=0.;

//  int isum;
  double *w_p2i3l = &p2i3l;
  double *w_r = &r;
  double *w_phi = &phi;
  double *w_phi2_1 = &phi2_1;
  double *w_thta = &thta;
  double *w_thta2_1 = &thta2_1;

  double dm = d_mass;  // deuteron (projectile)
  double pm = p_mass;  // proton (target)
  double nm = n_mass;  // neutron

  double m1 = pm;
  double m2 = pm;
  double m3 = nm;
  double ekin = bt;

  etot = dm + ekin + pm;

  phi2 = phi12;

// Initialization
//filell<<t1<<' '<<t2<<G4endl;
  ff = CLHEP::pi/180.;
  ph1r = 0.;
  th1r = t2*ff;
//  phi2=180.;//temporary 19IV2010 test purpose
  ph2r = phi2*ff;
  th2r = t1*ff;

  //G4cout <<"ph1r = "<<ph1r<<"\tth1r = "<<th1r<<"\tth2r = "<<th2r<<endl;


//  isum = 0;
  uni2[0] = sin(th2r)*cos(ph2r);			//czy nie za mala precyzja ????
  uni2[1] = sin(th2r)*sin(ph2r);
  uni2[2] = cos(th2r);

//Tracking of the kinematic curve.

  pol = P(ekin + dm,dm);
  e1min = m1;
  e1 = e2_0 + m1;
//  e1max = etot - m2 -m3;

  if (e1 < e1min) return 0;

  p1 = P(e1,m1);
  p1x = p1*sin(th1r)*cos(ph1r);
  p1y = p1*sin(th1r)*sin(ph1r);
  p1z = p1*cos(th1r);

  p2i3[0] = -p1x;
  p2i3[1] = -p1y;
  p2i3[2] = -p1z + pol;


  //p2i3ll = sqrt(p2i3[0]*p2i3[0] + p2i3[1]*p2i3[1] + p2i3[2]*p2i3[2]);
  carsph(p2i3,w_p2i3l,w_thta,w_phi);
//filell<<thta<<G4endl;			    	//return thta
  //G4cout <<"p2i3 = "<<p2i3<<"\t   p2i3l = "<<p2i3l<<"\t   thta = "<< thta<<"\t   phi = "<<phi<<endl;

//  p2i3_1[0] = X(p2i3[0],p2i3[1],p2i3[2],thta);
//  p2i3_1[1] = Y(p2i3[0],p2i3[1],p2i3[2],thta);
//  p2i3_1[2] = Z(p2i3[0],p2i3[1],p2i3[2],thta);
//filell<<p2i3_1[0]<<' '<<p2i3_1[1]<<' '<<p2i3_1[2]<<' ';test 19IV2010
  uni2_1[0] = X(uni2[0],uni2[1],uni2[2],thta);
  uni2_1[1] = Y(uni2[0],uni2[1],uni2[2],thta);
  uni2_1[2] = Z(uni2[0],uni2[1],uni2[2],thta);

  carsph(uni2_1,w_r,w_thta2_1,w_phi2_1);
//filell<<r<<' '<<thta2_1<<' '<<th1r<<' '<<th2r<<' '<<ph1r<<' '<<ph2r<<' '<<thta<<' '<<thta-th2r<<' '<<th2r-thta<<G4endl;//test 19IV2010
  //if (fabs(r - 1.) > 1.e-3) {G4cout <<"3!!! = "<<fabs(r - 1.);exit(p2i3l);}

  e2i3 = etot - e1;					// energy for p and n
  d = e2i3*e2i3 - p2i3l*p2i3l;
//G4cout <<"d = "<<d <<"\t drugi = "<<(m2 + m3)*(m2 + m3)<<endl;
//G4cout <<":2"<<endl;

  if (d < (m2 + m3)*(m2 + m3)) return 0;

  a = 4.*e2i3*e2i3 - 4.*p2i3l*p2i3l*cos(thta2_1)*cos(thta2_1);
  b = -4.*d*p2i3l*cos(thta2_1);
  c = (m2 + m3)*(m2 + m3)*e2i3*e2i3 - d*d;
  delta = b*b - 4.*a*c;
  if (delta < 0.)  return 0;
 
  p2[0] = (-b - sqrt(delta))/(2.*a);
  p2[1] = (-b + sqrt(delta))/(2.*a);

  e2[0] = E(p2[0],m2);
  e2[1] = E(p2[1],m2);

//  if (chooser%2==0) chooser=0;
  //output
  if (p2[0] >= 0.&&(typ==0||typ==2)) { e2p = e2[0]; pp2 = p2[0]; }
  else if (p2[1] >= 0.&&typ==1) { e2p = e2[1]; pp2 = p2[1]; }
  else return 0;
//G4cout<<p2[0]<<' '<<p2[1]<<G4endl;
//  chooser++;
//G4cout<<chooser<<' '<<p2[0]<<' '<<p2[1]<<' '<<pp2<<G4endl; DLACZEGO 2 RAZY TO JEST??
//         e1_0 = e1-m1
  e1_0 = e2p - m2;
  e3 = etot - e1 - e2p;
  e3_0 = e3- m3;

  pt2[0] = pp2*uni2[0];
  pt2[1] = pp2*uni2[1];
  pt2[2] = pp2*uni2[2];

  pp3 = sqrt(e3*e3 - m3*m3);
  pt3[0] = p2i3[0] - pt2[0];
  pt3[1] = p2i3[1] - pt2[1];
  pt3[2] = p2i3[2] - pt2[2];

  t3r = acos(pt3[2]/pp3);
  t3 = 180./CLHEP::pi*t3r;

  if (sin(t3r) == 0) cphi13 = pt3[0]/pp3;
  else cphi13 = pt3[0]/pp3/sin(t3r);

  phi13r = acos(cphi13);

  phi13 = 180./CLHEP::pi*phi13r;

  *wphi13 = phi13;
  *we1_0 = e1_0;
  *we3_0 = e3_0;


  bp_t[0] = p1;
  bp_t[1] = pp2;
  bp_t[2] = pp3;

  return double(delta);//1;
}

void Bina_PrimaryGeneratorAction::carsph(double *x,double *wr,double *wthta,double *wphi)
{
//Transformation from cartesian to spherical coordinates

  double r1=0.,r=0.,phi=0.,thta=0.;
  r1 = x[0]*x[0] + x[1]*x[1];
  r = sqrt(r1 + x[2]*x[2]);
  r1 = sqrt(r1);
  if (r < 1.e-6);
  else
  {
    if (r1 < 1.e-6)
    {
      phi = 0.;
      thta = 0.;
      if (x[2] < 0.) thta = CLHEP::pi;
    }
    else
    {
      thta = acos(x[2]/r);
      phi = atan2(x[1],x[0]);
    }
  }

  *wr = r;
  *wphi = phi;
  *wthta = thta;
}

double Bina_PrimaryGeneratorAction::rinterp(double ytmp[5][5][5][5],double x1,double x2, double x3, double x4)
{
//#include"fstream.h"
////#include"G4ios.hh"
//#include"iostream.h"
//interpolation Y-value in 4-dimensial space (4,4,4,4). spline method is used in all dimensions.
// used for break-up cross section interpolation

// x1 -  theta1 table ( first proton )
// x2 -  theta2 table ( second proton)
// x3 -  phi12 angle table
// x4  - kinetic energy table
  double y_t1[5][5][5][5];

  double /*x1=0.,x2=0.,x3=0.,x4=0.,*/y=0.;//10IV2010
  int m1=4,m2=4,m3=4,m4=4;
//  int iiii,jjjj,kkkk,llll;
  double *ya = &ytmp[0][0][0][0];
  double *y3a = &y_t1[0][0][0][0];
  double *w_y = &y;

//  const double *ya_const, *y3a_const;
//  filell.open("fileoutput.txt",ios::out|ios::app);
//  ya_const = ya;	
//  y3a_const = y3a;
//  for (iiii=1;iiii<5;iiii++) {
//      for (jjjj=1;jjjj<5;jjjj++) {
//             for (kkkk=1;kkkk<5;kkkk++) {
//                     for (llll=1;llll<5;llll++) {
//        filell<<iiii<<' '<<jjjj<<' '<<kkkk<<' '<<llll<<' '<<ytmp[iiii][jjjj][kkkk][llll]<<' '<<G4endl;}}}}
/////////////////  spline4(x1a_r,x2a_r,x3a_r,x4a_r,ya,m1,m2,m3,m4,y3a);
/////////////////  splint4(x1a_r,x2a_r,x3a_r,x4a_r,ya,y3a,m1,m2,m3,m4,x1,x2,x3,x4,w_y);
  spline4(ya,m1,m2,m3,m4,y3a);
  splint4(ya,y3a,m1,m2,m3,m4,x1,x2,x3,x4,w_y);
//  filell<<"Yyyyyyyyyyyyyyy "<<y*100.<<G4endl;
   return y;
}

void Bina_PrimaryGeneratorAction::spline4(double *ya,
					int m1,int m2,int m3,int m4,double *y3a)
{
  const double *ya_const, *y3a_const;
  double *x4a_mov;
  ya_const = ya;
  y3a_const = y3a;
  x4a_mov = x4a_r;

  ya += 126;
  y3a += 126;
  x4a_mov++;

  for (int j=1;j<=m1;j++)
  {
    ya+=25;
    y3a+=25;
    for (int k=1;k<=m2;k++)
    {
      ya+=5;
      y3a+=5;
      for (int l=1;l<=m3;l++)
      {
	spline(x4a_mov,ya,m4,1.0e30,1.0e30,y3a);
	ya  +=5;
        y3a +=5;
      }
    }
  }
  ya = (double*)ya_const;
  y3a = (double*)y3a_const;
//  x4a_r = (double*)x4a_const;
}

void Bina_PrimaryGeneratorAction::splint4(double *ya,double *y2a,
					int m1,int m2,int m3,int m4,double x1,double x2,double x3,double x4,double *y)
{
  double *ytmp,*yytmp,*yyytmp;
  double *ytmp1,*yytmp1,*yyytmp1;
  const double *ya_const,*y2a_const;
//  int jjj;
  ya_const = ya;
  y2a_const = y2a;

  ya += 126;
  y2a += 126;

  ytmp=vector(1,m1);
  yytmp=vector(1,m1);
  yyytmp=vector(1,m1);
  ytmp1=vector(1,m1);
  yytmp1=vector(1,m1);
  yyytmp1=vector(1,m1);

  for (int j=1;j<=m1;j++)
  {
    ya+=25;
    y2a+=25;
    for (int k=1;k<=m2;k++)
    {
      ya+=5;
      y2a+=5;
      for (int l=1;l<=m3;l++)
      {
        splint(x4a_r,ya,y2a,m4,x4,&yytmp[l]);
        ya  += 5;
        y2a +=5;
       // for (jjj=1;jjj<=4;jjj++) filell<<"1-  "<<l <<' '<<ya[jjj]<<' '<<y2a[jjj]<<' '<<yytmp[jjj]<<G4endl;
        }
      spline(x3a_r,yytmp,m3,1.0e30,1.0e30,ytmp);
      splint(x3a_r,yytmp,ytmp,m3,x3,&yyytmp[k]);
    //   for (jjj=1;jjj<=4;jjj++) filell<<"2- "<<k<<' '<<yytmp[jjj]<<' '<<ytmp[jjj]<<' '<<yyytmp[jjj]<<G4endl;
      }
    spline(x2a_r,yyytmp,m2,1.0e30,1.0e30,yytmp1);
    splint(x2a_r,yyytmp,yytmp1,m2,x2,&ytmp1[j]);
//     for (jjj=1;jjj<=4;jjj++) filell<<"3- "<<j<<' '<<yyytmp[jjj]<<' '<<yytmp1[jjj]<<' '<<ytmp1[jjj]<<G4endl;
    }
  spline(x1a_r,ytmp1,m1,1.0e30,1.0e30,yyytmp1);
  splint(x1a_r,ytmp1,yyytmp1,m1,x1,y);
  //for (int j=1;j<=m1;j++) {filell<<"4- "<<ytmp1[j]<<' '<<yytmp1[j]<<' '<<y[j]<<G4endl;}
//  filell<<"yyyyyyyyyy "<<*y<<G4endl;
  free_vector(yyytmp,1,m1);
  free_vector(yytmp,1,m1);
  free_vector(ytmp,1,m1);
  free_vector(yyytmp1,1,m1);
  free_vector(yytmp1,1,m1);
  free_vector(ytmp1,1,m1);

  ya = (double*)ya_const;
  y2a = (double*)y2a_const;
//  filell<<"splint4 "<<*y<<G4endl;
}

void Bina_PrimaryGeneratorAction::splint(double* xa,double* ya,double* y2a,int n,double x,double *y)
{
//	  Given the arrays xa[1..n] and ya[1..n], which tabulate a function (with the xai's in order),
//	  and given the array y2a[1..n], which is the output from spline above, and given a value of x,
//  	  this routine returns a cubic-spline interpolated value y.
//int aa;
  int klo,khi,k;
  double h,b,a;

  klo = 1;		//We will  nd the right place in the table by means of bisection.
  khi = n; 		//This is optimal if sequential calls to this routine are at random values of x.
		    	//If sequential calls are in order, and closely spaced, one would do better to store
		    	//previous values of klo and khi and test if they remain appropriate on the next call.
  while (khi - klo > 1)
  {
    k = (khi + klo) >> 1;
    if (xa[k] > x) khi = k;
    else klo = k;
  }			//klo and khi now bracket the input value of x.
  h = xa[khi] - xa[klo];
  if (h == 0.0) {G4cout <<"Bad xa input to function splint()"<<G4endl; exit(1);}	//The xa's must be distinct
//for (aa=0;aa<=4;aa++) {
//  filell<<xa[aa]<<G4endl;}
  a = (xa[khi] - x)/h;
  b = (x - xa[klo])/h;							//Cubic spline polynomial is now evaluated

  *y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
//******************************************************************
//filell<<"Klo "<<klo<<' '<<khi<<' '<<xa[klo]<<' '<<x<<' '<<xa[khi]<<' '<<ya[klo-1]<<' '<<ya[khi-1]<<' '<<*y<<G4endl;
//!!!!!!!!!!!!!!!!!!!!!
//for(aa=1;aa<5;aa++){
//filell<<aa<<' '<<xa[aa]<<' '<<ya[aa-1]<<' '<<y2a[aa]<<G4endl;}
//filell<<*y<<G4endl;
/*h=0.;
for (k=1;k<5;k++) h+=ya[k];// filell<<h<<' '<<ya[k]<<G4endl;}
h/=4.;
//filell<<"TUTUTU "<<h<<G4endl;
*y=h;
//filell<<y<<' '<<h<<G4endl;
*/
}
void Bina_PrimaryGeneratorAction::spline(double x[], double y[], int n, double yp1, double ypn, double *y2)
{
  //Given arrays x[1..n] and y[1..n] containing a tabulated function, i.e., yi = f(xi),
  //with x1 <x2 < :: : < xN, and given values yp1 and ypn for the  rst derivative of
  //the interpolating function at points 1 and n, respectively, this routine returns an
  //array y2[1..n] that contains the second derivatives of the interpolating function at
  //the tabulated points xi. If yp1 and/or ypn are equal to 1   1030 or larger, the routine
  //is signaled to set the corresponding boundary condition for a natural spline, with zero]
  // second derivative on that boundary.

  int i,k;
  double p,qn,sigl,un,*u;
  u=vector(1,n-1);

  if (yp1 > 0.99e30) 		//The lower boundary condition is set either to be natural
    y2[1]=u[1]=0.0;
  else
  { 		//or else to have a speci ed  rst derivative.
    y2[1] = -0.5;
    u[1]=(3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
  }

  for (i=2;i<=n-1;i++)
  { 					//This is the decomposition loop of the tridiagonal algorithm
    sigl=(x[i]-x[i-1])/(x[i+1]-x[i-1]);	//y2 and u are used for temporary storage of the decomposed factors
    p=sigl*y2[i-1]+2.0;
    y2[i]=(sigl-1.0)/p;
    u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
    u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sigl*u[i-1])/p;
  }
  if (ypn > 0.99e30) 			//The upper boundary condition is set either to be natural
    qn=un=0.0;
  else 					//or else to have a specified first derivative
  {
    qn=0.5;
    un=(3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]));
  }
  y2[n]=(un-qn*u[n-1])/(qn*y2[n-1]+1.0);
  for (k=n-1;k>=1;k--) 			//This is the backsubstitution loop of the tridiagonal algorithm
    y2[k]=y2[k]*y2[k+1]+u[k];

  free_vector(u,1,n-1);
}
double *Bina_PrimaryGeneratorAction::vector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;

	v=(double *)malloc((size_t) ((nh-nl+1+1)*sizeof(double)));
	if (!v) { G4cout<<"allocation failure in vector()"<<G4endl; exit(1);}
	return v-nl+1;
}

void Bina_PrimaryGeneratorAction::free_vector(double *v, long nl, long /*nh*/)
/* free a double vector allocated with vector() */
{
	free((char*) (v+nl-1));
}

double Bina_PrimaryGeneratorAction::cspol(double a_x,double a_y,double a_xx,double a_xy,double a_yy)
{
  double c1,c2,c3,c4,c5,phi1,phi1r,phipol,csp,cons,u;
  c1 = -1.5 * pz*a_x;
  c2 =  1.5 * pz*a_y;
  c3 =      - pzz*a_xy;
  c4 =  0.5 * pzz*a_xx;
  c5 =  0.5 * pzz*a_yy;
// von Neumann's method applied to phi
  do
  {
    phi1 = -180. + 360.*RandomFlat(aj1phi);
    phi1r = CLHEP::pi/180.*phi1;
// angle between polarization and y-axis,  is it right?
    phipol = phi1r;
    csp = 1. + c1*sin(phipol) + c2*cos(phipol) + c3*sin(phipol)*cos(phipol) +
          c4*sin(phipol)*sin(phipol) + c5*cos(phipol)*cos(phipol);
    cons = 0.5;
    u = 1. + cons*(2.*RandomFlat(aju) - 1.);
  }
  while (u > csp);
  return phi1;
}

void Bina_PrimaryGeneratorAction::ugelast_read(void)
{
  std::ifstream file;
  file.open("./data/sig_anpow.dat");
  if (!file)
  {
    G4cout <<"Cannot open the file : ../data/sig_anpow.dat !!!"<<G4endl;
    exit(1);
   }
  for (int i=0;i<70;i++)
  {
    file >> bt1_mev>>th>>thl[i]>>ds[i]>>ay[i]>>ap1[i]>>ap2[i]>>t21>>ap3[i];

    if(file.eof())
    {
      G4cout << " Unexpected end of file : /home/sworst/Breakup/data/sig_anpow.dat !!! Index : "<<i<<G4endl;
      exit(1);
    }
    if (file.fail())
    {
      G4cout << "Unexpected error in file : /home/sworst/Breakup/data/sig_anpow.dat !!! Index : "<<i<<G4endl;
      exit(1);
    }
  }
  file.close();
}

void Bina_PrimaryGeneratorAction::break_read(void)
{//read break-up observables (cross section, analysing powers) from files

  int ns;
  int ith1,ok;

  std::ifstream file[4];

  file[0].open("./data/arr1.int");
  file[1].open("./data/arr2.int");
  file[2].open("./data/arr3.int");
  file[3].open("./data/arr4.int");
  int l;
  int k,j,i;
  for(l=0;l<4;l++)
  {
    if (!file[l])
    {
      G4cout <<"Cannot open the file : file["<<l<<"]"<<G4endl;
      exit(1);
    }
  }
  ith1 = 0;
  ok = 1;
  for (l=0;l<4;l++)
  {ok=1;
    do
    {
      for(k=0;k<14;k++)
      {
        if (ok == 0) break;
        for(j=0;j<13;j++)
        {
          file[l].ignore(100,10);
          file[l].ignore(100,10);
          file[l]>>ns;
          if(file[l].eof()) {ok = 0; break;}
          for (i=0;i<ns;i++)
          {
            file[l] >> sig[ith1][k][j][i] >> ayn[ith1][k][j][i] >> ayd[ith1][k][j][i]
	            >> axx[ith1][k][j][i] >> ayy[ith1][k][j][i] >> axn[ith1][k][j][i]
	            >> axy[ith1][k][j][i] >> axz[ith1][k][j][i];
          if (sig[ith1][k][j][i]>csmax&&ith1<(themax-offs[0])/step[0]+1.&&
              ith1>(themin-offs[0])/step[0]-1.&&k<(themax2-offs[1])/step[1]+1.&&
              k>(themin2-offs[1])/step[1]-1.)  csmax=sig[ith1][k][j][i];
          }
        }
      }
      ith1++;
    }
    while (ok == 1);
  ith1--;
  }
  for(l=0;l<4;l++) file[l].close();
csmax*=1.5;
}
