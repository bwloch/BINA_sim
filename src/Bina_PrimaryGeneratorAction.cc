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

#include <fstream>
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
  bt1 = 0.;		bt2 = 0.;		bt3 = 0.; bt4=0.;
  t1  = 0.;		t2  = 0.;		t3  = 0.; t4=0.;
  fi1 = 0.;             fi2 = 0.;		fi3 = 0.; fi4=0.;
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
  open_pluto_file();

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
    event_cleaner_particle_gun = new G4ParticleGun(n_particle);
    particle = particleTable->FindParticle("e-");
    event_cleaner_particle_gun->SetParticleDefinition(particle);
  }
  else if (npd_choice == 0)
  {	

    particleGun1 = new G4ParticleGun(n_particle);
    particleGun2 = new G4ParticleGun(n_particle);
    event_cleaner_particle_gun = new G4ParticleGun(n_particle);

    particle = particleTable->FindParticle("proton");
    particleGun1->SetParticleDefinition(particle);

    particle = particleTable->FindParticle("deuteron");
    particleGun2->SetParticleDefinition(particle);
    
        particle = particleTable->FindParticle("e-");
    event_cleaner_particle_gun->SetParticleDefinition(particle);
  }
    else if (npd_choice == 1)
  {	

    particleGun1 = new G4ParticleGun(n_particle);
    particleGun2 = new G4ParticleGun(n_particle);
    event_cleaner_particle_gun = new G4ParticleGun(n_particle);

    particle = particleTable->FindParticle("deuteron");
    particleGun1->SetParticleDefinition(particle);

    particle = particleTable->FindParticle("deuteron");
    particleGun2->SetParticleDefinition(particle);
    
        particle = particleTable->FindParticle("e-");
    event_cleaner_particle_gun->SetParticleDefinition(particle);
  }
  else if (npd_choice == 2)
  {
      event_cleaner_particle_gun = new G4ParticleGun(n_particle);
    particleGun1 = new G4ParticleGun(n_particle);
    particleGun2 = new G4ParticleGun(n_particle);
    particleGun3 = new G4ParticleGun(n_particle);

    particle = particleTable->FindParticle("proton");
    particleGun1->SetParticleDefinition(particle);

    G4ParticleDefinition* particle2 = particleTable->FindParticle("proton");
    particleGun2->SetParticleDefinition(particle2);

    particle = particleTable->FindParticle("neutron");
    particleGun3->SetParticleDefinition(particle);
    
            particle = particleTable->FindParticle("e-");
    event_cleaner_particle_gun->SetParticleDefinition(particle);

  }
    else if (npd_choice == 3)
  {
      event_cleaner_particle_gun = new G4ParticleGun(n_particle);
    particleGun1 = new G4ParticleGun(n_particle);
    particleGun2 = new G4ParticleGun(n_particle);
    particleGun3 = new G4ParticleGun(n_particle);
    
    particle = particleTable->FindParticle("e-");
    event_cleaner_particle_gun->SetParticleDefinition(particle);

    particle = particleTable->FindParticle("deuteron");
    particleGun1->SetParticleDefinition(particle);

    G4ParticleDefinition* particle2 = particleTable->FindParticle("proton");
    particleGun2->SetParticleDefinition(particle2);

    particle = particleTable->FindParticle("neutron");
    particleGun3->SetParticleDefinition(particle);
    

  }
  else if(npd_choice == 4)
  {

    event_cleaner_particle_gun = new G4ParticleGun(n_particle);
    particleGun1 = new G4ParticleGun(n_particle);
    particleGun2 = new G4ParticleGun(n_particle);
    particleGun3 = new G4ParticleGun(n_particle);
    particleGun4 = new G4ParticleGun(n_particle);
    
    particle = particleTable->FindParticle("proton");
    particleGun1->SetParticleDefinition(particle);

    particle = particleTable->FindParticle("proton");
    particleGun2->SetParticleDefinition(particle);

    particle = particleTable->FindParticle("neutron");
    particleGun3->SetParticleDefinition(particle);
    
    particle = particleTable->FindParticle("neutron");
    particleGun4->SetParticleDefinition(particle);
    
            particle = particleTable->FindParticle("e-");
    event_cleaner_particle_gun->SetParticleDefinition(particle);
  }
    else if(npd_choice == 5)
  {

    event_cleaner_particle_gun = new G4ParticleGun(n_particle);
    particleGun1 = new G4ParticleGun(n_particle);
    particleGun2 = new G4ParticleGun(n_particle);

    
    particle = particleTable->FindParticle("proton");
    particleGun1->SetParticleDefinition(particle);

    particle = particleTable->FindParticle("triton");
    particleGun2->SetParticleDefinition(particle);

    
            particle = particleTable->FindParticle("e-");
    event_cleaner_particle_gun->SetParticleDefinition(particle);
  }
  else
  {
    G4cout << "Error in configuration file geo.mac !!! Bad npd_choice = "<<npd_choice<<G4endl;
    exit(1);
  }

}

Bina_PrimaryGeneratorAction::~Bina_PrimaryGeneratorAction()
{
delete event_cleaner_particle_gun;
  delete particleGun1;
  if (npd_choice >= 0) delete particleGun2;
  if (npd_choice >=2 ) delete particleGun3;
  if (npd_choice>=4) delete particleGun4;
}

void Bina_PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{


  GetChoice(npd_choice);
  	double pluto_momentum[4];
  	double pluto_momentum1[4];
  	double pluto_momentum2[4];
  	double pluto_momentum3[4];
  	double pluto_momentum4[4];
  if (npd_choice <  0)
  {
   // G4cout <<"Npd_choice : "<<npd_choice<<G4endl;
   
   	Pos();
   	
/// Black Magic begin
//generujemy dodatkowa czastke 
//zeby zapewnic zapic wszystkich czastek w tym samym evencie
// bo czastki zapisywane sa dopiero po rozpoczeciu symulacji kolejnej!!
  	event_cleaner_particle_gun->SetParticleMomentumDirection(G4ThreeVector(0.0,0.0,-1.0));
    	event_cleaner_particle_gun->SetParticleEnergy(0.0001*CLHEP::GeV);
    	event_cleaner_particle_gun->SetParticlePosition(G4ThreeVector(vertex[0],vertex[1],vertex[2]));
    	event_cleaner_particle_gun->GeneratePrimaryVertex(anEvent);
// Black Magic end	
  	
  	//reading particle
	read_part_momentum(pluto_momentum);
  	v4.setE(pluto_momentum[0]);
  	v4.setPx(pluto_momentum[1]);
  	v4.setPy(pluto_momentum[2]);
  	v4.setPz(pluto_momentum[3]);
  	particleGun1->SetParticleMomentumDirection(v4);
  	bt1=(v4.e()-v4.m());
  	t1=v4.theta();
  	fi1=v4.phi();
  	particleGun1->SetParticleEnergy(bt1*CLHEP::GeV);
  	particleGun1->SetParticlePosition(G4ThreeVector(vertex[0],vertex[1],vertex[2]));
  	particleGun1->GeneratePrimaryVertex(anEvent);

      GetStartEnergy(bt1*1000.);
      GetStartAngleTheta(t1*180./M_PI);
      GetStartAnglePhi(fi1*180./M_PI);
      GetStartPosition(vertex[0],vertex[1],vertex[2]);

  }
  else if (npd_choice == 0)
  {
  	Pos();

/// Black Magic begin
//generujemy dodatkowa czastke 
//zeby zapewnic zapic wszystkich czastek w tym samym evencie
// bo czastki zapisywane sa dopiero po rozpoczeciu symulacji kolejnej!!
  	event_cleaner_particle_gun->SetParticleMomentumDirection(G4ThreeVector(0.0,0.0,-1.0));
    	event_cleaner_particle_gun->SetParticleEnergy(0.0001*CLHEP::GeV);
    	event_cleaner_particle_gun->SetParticlePosition(G4ThreeVector(vertex[0],vertex[1],vertex[2]));
    	event_cleaner_particle_gun->GeneratePrimaryVertex(anEvent);
// Black Magic end	
    	
  	//reading deut
	read_part_momentum(pluto_momentum);
  	v4.setE(pluto_momentum[0]);
  	v4.setPx(pluto_momentum[1]);
  	v4.setPy(pluto_momentum[2]);
  	v4.setPz(pluto_momentum[3]);
  	particleGun2->SetParticleMomentumDirection(v4);
  	bt1=(v4.e()-v4.m());
  	t1=v4.theta();
  	fi1=v4.phi();
  	particleGun2->SetParticleEnergy(bt1*CLHEP::GeV);
  	particleGun2->SetParticlePosition(G4ThreeVector(vertex[0],vertex[1],vertex[2]));

  	//reading prot
  	//G4cout<<"\n\t Mylog: PrimGneAct Ed="<<pluto_momentum[0]<<" "
  	read_part_momentum(pluto_momentum);

  	v4.setE(pluto_momentum[0]);
  	v4.setPx(pluto_momentum[1]);
  	v4.setPy(pluto_momentum[2]);
  	v4.setPz(pluto_momentum[3]);
  	bt2=(v4.e()-v4.m());

  	t2=v4.theta();
  	fi2=v4.phi();
  	particleGun1->SetParticleMomentumDirection(v4);
  	particleGun1->SetParticleEnergy(bt2*CLHEP::GeV);
  	particleGun1->SetParticlePosition(G4ThreeVector(vertex[0],vertex[1],vertex[2]));


      GetStartEnergy(bt1*1000.,bt2*1000.);
      GetStartAngleTheta(t1*180./M_PI,t2*180./M_PI);
      GetStartAnglePhi(fi1*180./M_PI,fi2*180./M_PI);
      GetStartPosition(vertex[0],vertex[1],vertex[2]);
      
      particleGun1->GeneratePrimaryVertex(anEvent);
    	particleGun2->GeneratePrimaryVertex(anEvent);
    	

  }
  else if (npd_choice == 1)
  {
    //G4cout <<"Npd_choice : "<<npd_choice<<G4endl;


    Pos();
    
    /// Black Magic begin
//generujemy dodatkowa czastke 
//zeby zapewnic zapic wszystkich czastek w tym samym evencie
// bo czastki zapisywane sa dopiero po rozpoczeciu symulacji kolejnej!!
  	event_cleaner_particle_gun->SetParticleMomentumDirection(G4ThreeVector(0.0,0.0,-1.0));
    	event_cleaner_particle_gun->SetParticleEnergy(0.0001*CLHEP::GeV);
    	event_cleaner_particle_gun->SetParticlePosition(G4ThreeVector(vertex[0],vertex[1],vertex[2]));
    	event_cleaner_particle_gun->GeneratePrimaryVertex(anEvent);
// Black Magic end

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
  else if(npd_choice == 2 || npd_choice == 3){
  
  
	Pos();
	int Is_Neutron_Simulated=0;
	
/// Black Magic begin
//generujemy dodatkowa czastke 
//zeby zapewnic zapic wszystkich czastek w tym samym evencie
// bo czastki zapisywane sa dopiero po rozpoczeciu symulacji kolejnej!!
  	event_cleaner_particle_gun->SetParticleMomentumDirection(G4ThreeVector(0.0,0.0,-1.0));
    	event_cleaner_particle_gun->SetParticleEnergy(0.0001*CLHEP::GeV);
    	event_cleaner_particle_gun->SetParticlePosition(G4ThreeVector(vertex[0],vertex[1],vertex[2]));
    	event_cleaner_particle_gun->GeneratePrimaryVertex(anEvent);
// Black Magic end
  	
  	//reading deuteron
	read_part_momentum(pluto_momentum1);
  	v1.setE(pluto_momentum1[0]);
  	v1.setPx(pluto_momentum1[1]);
  	v1.setPy(pluto_momentum1[2]);
  	v1.setPz(pluto_momentum1[3]);
  	particleGun1->SetParticleMomentumDirection(v1);
  	bt1=(v1.e()-v1.m());
  	t1=v1.theta();
  	fi1=v1.phi();
  	particleGun1->SetParticleEnergy(bt1*CLHEP::GeV);
  	particleGun1->SetParticlePosition(G4ThreeVector(vertex[0],vertex[1],vertex[2]));
  	particleGun1->GeneratePrimaryVertex(anEvent);
  	//reading proton
  	read_part_momentum(pluto_momentum2);
  	v2.setE(pluto_momentum2[0]);
  	v2.setPx(pluto_momentum2[1]);
  	v2.setPy(pluto_momentum2[2]);
  	v2.setPz(pluto_momentum2[3]);
  	bt2=(v2.e()-v2.m());
  	t2=v2.theta();
  	fi2=v2.phi();
  	particleGun2->SetParticleMomentumDirection(v2);
  	particleGun2->SetParticleEnergy(bt2*CLHEP::GeV);
  	particleGun2->SetParticlePosition(G4ThreeVector(vertex[0],vertex[1],vertex[2]));
  	particleGun2->GeneratePrimaryVertex(anEvent);
  	if(Is_Neutron_Simulated==1){
  	//reading neutron
  	
  
  		read_part_momentum(pluto_momentum3);
	  	v3.setE(pluto_momentum3[0]);
	  	v3.setPx(pluto_momentum3[1]);
	  	v3.setPy(pluto_momentum3[2]);
	  	v3.setPz(pluto_momentum3[3]);
	  	bt3=(v3.e()-v3.m());
	  	t3=v3.theta();
	  	fi3=v3.phi();
	  	particleGun3->SetParticleMomentumDirection(v3);
	  	particleGun3->SetParticleEnergy(bt3*CLHEP::GeV);
	  	particleGun3->SetParticlePosition(G4ThreeVector(vertex[0],vertex[1],vertex[2]));
	  	particleGun3->GeneratePrimaryVertex(anEvent);
	  
  	}
  	else {
	 	 read_part_momentum(pluto_momentum);
	}
  	
      GetStartEnergy(bt1*1000.,bt2*1000.,bt3*1000.);
      GetStartAngleTheta(t1*180./M_PI,t2*180./M_PI,t3*180./M_PI);
      GetStartAnglePhi(fi1*180./M_PI,fi2*180./M_PI,fi3*180./M_PI);
      GetStartPosition(vertex[0],vertex[1],vertex[2]);
      
  }
  else if(npd_choice==4){

  	Pos();
	int Is_Neutron_Simulated=0;
	
/// Black Magic begin
//generujemy dodatkowa czastke 
//zeby zapewnic zapic wszystkich czastek w tym samym evencie
// bo czastki zapisywane sa dopiero po rozpoczeciu symulacji kolejnej!!
  	event_cleaner_particle_gun->SetParticleMomentumDirection(G4ThreeVector(0.0,0.0,-1.0));
    	event_cleaner_particle_gun->SetParticleEnergy(0.0001*CLHEP::GeV);
    	event_cleaner_particle_gun->SetParticlePosition(G4ThreeVector(vertex[0],vertex[1],vertex[2]));
    	event_cleaner_particle_gun->GeneratePrimaryVertex(anEvent);
// Black Magic end	
	
	
	read_part_momentum(pluto_momentum1);
  	v1.setE(pluto_momentum1[0]);
  	v1.setPx(pluto_momentum1[1]);
  	v1.setPy(pluto_momentum1[2]);
  	v1.setPz(pluto_momentum1[3]);
  	particleGun1->SetParticleMomentumDirection(v1);
  	bt1=(v1.e()-v1.m());
  	t1=v1.theta();
  	fi1=v1.phi();
  	particleGun1->SetParticleEnergy(bt1*CLHEP::GeV);
  	particleGun1->SetParticlePosition(G4ThreeVector(vertex[0],vertex[1],vertex[2]));
  	particleGun1->GeneratePrimaryVertex(anEvent);
  	//reading proton
  	read_part_momentum(pluto_momentum2);
  	v2.setE(pluto_momentum2[0]);
  	v2.setPx(pluto_momentum2[1]);
  	v2.setPy(pluto_momentum2[2]);
  	v2.setPz(pluto_momentum2[3]);
  	bt2=(v2.e()-v2.m());
  	t2=v2.theta();
  	fi2=v2.phi();
  	particleGun2->SetParticleMomentumDirection(v2);
  	particleGun2->SetParticleEnergy(bt2*CLHEP::GeV);
  	particleGun2->SetParticlePosition(G4ThreeVector(vertex[0],vertex[1],vertex[2]));
  	particleGun2->GeneratePrimaryVertex(anEvent);
  	  	if(Is_Neutron_Simulated==1){
  	  	
  	  	read_part_momentum(pluto_momentum3);
	  	v3.setE(pluto_momentum3[0]);
	  	v3.setPx(pluto_momentum3[1]);
	  	v3.setPy(pluto_momentum3[2]);
	  	v3.setPz(pluto_momentum3[3]);
	  	bt3=(v3.e()-v3.m());
	  	t3=v3.theta();
	  	fi3=v3.phi();
	  	particleGun3->SetParticleMomentumDirection(v3);
	  	particleGun3->SetParticleEnergy(bt3*CLHEP::GeV);
	  	particleGun3->SetParticlePosition(G4ThreeVector(vertex[0],vertex[1],vertex[2]));
	  	particleGun3->GeneratePrimaryVertex(anEvent);
	  	
	  	read_part_momentum(pluto_momentum4);
	  	v4.setE(pluto_momentum4[0]);
	  	v4.setPx(pluto_momentum4[1]);
	  	v4.setPy(pluto_momentum4[2]);
	  	v4.setPz(pluto_momentum4[3]);
	  	bt4=(v4.e()-v4.m());
	  	t4=v4.theta();
	  	fi4=v4.phi();
	  	particleGun3->SetParticleMomentumDirection(v4);
	  	particleGun3->SetParticleEnergy(bt4*CLHEP::GeV);
	  	particleGun3->SetParticlePosition(G4ThreeVector(vertex[0],vertex[1],vertex[2]));
	  	particleGun3->GeneratePrimaryVertex(anEvent);
  
 	}
	else {
	 	 read_part_momentum(pluto_momentum);
	 	 read_part_momentum(pluto_momentum);
	}
      GetStartEnergy(bt1*1000.,bt2*1000.,bt3*1000.,bt4*1000.);
      GetStartAngleTheta(t1*180./M_PI,t2*180./M_PI,t3*180./M_PI,t4*180./M_PI);
      GetStartAnglePhi(fi1*180./M_PI,fi2*180./M_PI,fi3*180./M_PI,fi4*180./M_PI);
      GetStartPosition(vertex[0],vertex[1],vertex[2]);


}
  else if (npd_choice == 5)
  {
  	Pos();

/// Black Magic begin
//generujemy dodatkowa czastke 
//zeby zapewnic zapic wszystkich czastek w tym samym evencie
// bo czastki zapisywane sa dopiero po rozpoczeciu symulacji kolejnej!!
  	event_cleaner_particle_gun->SetParticleMomentumDirection(G4ThreeVector(0.0,0.0,-1.0));
    	event_cleaner_particle_gun->SetParticleEnergy(0.0001*CLHEP::GeV);
    	event_cleaner_particle_gun->SetParticlePosition(G4ThreeVector(vertex[0],vertex[1],vertex[2]));
    	event_cleaner_particle_gun->GeneratePrimaryVertex(anEvent);
// Black Magic end	
    	
  	//reading deut
	read_part_momentum(pluto_momentum);
  	v4.setE(pluto_momentum[0]);
  	v4.setPx(pluto_momentum[1]);
  	v4.setPy(pluto_momentum[2]);
  	v4.setPz(pluto_momentum[3]);
  	particleGun2->SetParticleMomentumDirection(v4);
  	bt1=(v4.e()-v4.m());
  	t1=v4.theta();
  	fi1=v4.phi();
  	particleGun2->SetParticleEnergy(bt1*CLHEP::GeV);
  	particleGun2->SetParticlePosition(G4ThreeVector(vertex[0],vertex[1],vertex[2]));

  	//reading prot
  	//G4cout<<"\n\t Mylog: PrimGneAct Ed="<<pluto_momentum[0]<<" "
  	read_part_momentum(pluto_momentum);

  	v4.setE(pluto_momentum[0]);
  	v4.setPx(pluto_momentum[1]);
  	v4.setPy(pluto_momentum[2]);
  	v4.setPz(pluto_momentum[3]);
  	bt2=(v4.e()-v4.m());

  	t2=v4.theta();
  	fi2=v4.phi();
  	particleGun1->SetParticleMomentumDirection(v4);
  	particleGun1->SetParticleEnergy(bt2*CLHEP::GeV);
  	particleGun1->SetParticlePosition(G4ThreeVector(vertex[0],vertex[1],vertex[2]));


      GetStartEnergy(bt1*1000.,bt2*1000.);
      GetStartAngleTheta(t1*180./M_PI,t2*180./M_PI);
      GetStartAnglePhi(fi1*180./M_PI,fi2*180./M_PI);
      GetStartPosition(vertex[0],vertex[1],vertex[2]);
      
      particleGun1->GeneratePrimaryVertex(anEvent);
    	particleGun2->GeneratePrimaryVertex(anEvent);
    	

  }
}

void Bina_PrimaryGeneratorAction::RandomInit(int level )
{
	CLHEP::Ranlux64Engine *Engine1 = new CLHEP::Ranlux64Engine(aj1bx,level);
	CLHEP::Ranlux64Engine *Engine2 = new CLHEP::Ranlux64Engine(aj1by,level);

  CLHEP::RandGauss *GD1 = new CLHEP::RandGauss(*Engine1);
  CLHEP::RandGauss *GD2 = new CLHEP::RandGauss(*Engine2);
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
void Bina_PrimaryGeneratorAction::open_pluto_file(){
	char file_path[500];
	if(npd_choice==0) sprintf(file_path,"data/output_elast_dp.txt");
	if(npd_choice==2) sprintf(file_path,"data/output_break_dp.txt");
	if(npd_choice==3) sprintf(file_path,"data/output_break_dd.txt");	
	if(npd_choice==4) sprintf(file_path,"data/output_break_dd_ppnn.txt");		
	if(npd_choice==5) sprintf(file_path,"data/output_break_pT.txt");	
	file_Pluto_generator.open(file_path,std::ios::in);
	if( !file_Pluto_generator.is_open()){
		G4cout<<"\n\n MyLog: Plik nie zostal otwarty\n";
		exit(1);
	}
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














void Bina_PrimaryGeneratorAction::read_part_momentum(double pluto_momentum[]){
if( !file_Pluto_generator.is_open()){
	G4cout<<"Plik nie zostal otwarty\n";
	exit(666);
	}
if(!file_Pluto_generator.eof()){
	for(int i=0;i<4;i++) file_Pluto_generator>>pluto_momentum[i];	
	}
	else {
	G4cout<<"Plik skończon\n";
	exit(666);
	}
}

void Bina_PrimaryGeneratorAction::read_part_momentum4(double pluto_momentum[]){
if( !file_Pluto_generator.is_open()){
	G4cout<<"Plik nie zostal otwarty\n";
	exit(666);
	}
if(!file_Pluto_generator.eof()){
	for(int i=0;i<4;i++) file_Pluto_generator>>pluto_momentum[i];	
	}
	else {
	G4cout<<"Plik skończon\n";
	exit(666);
	}
}


