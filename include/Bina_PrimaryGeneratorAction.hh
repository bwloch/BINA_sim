#ifndef Bina_PrimaryGeneratorAction_h
#define Bina_PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"

#include "CLHEP/Random/Randomize.h"
#include "CLHEP/Random/RandomEngine.h"
//#include "CLHEP/Random/RandGauss.h"
//#include "CLHEP/Random/RandGaussQ.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/Ranlux64Engine.h"
#include <fstream>
#include <iostream>
#include "G4LorentzVector.hh"
#include "MyFileReader.hh"

class Bina_DetectorConstruction;
class G4ParticleGun;
class G4Event;
class MyFileReader;


class Zmienne;

class Bina_PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    Bina_PrimaryGeneratorAction();
   ~Bina_PrimaryGeneratorAction();
//////////////////////////////
// Generation Parameters
//////////////////////////////
    double bt1,bt2,bt3,bt4,t1,t2,t3,t4,fi1,fi2,fi3,fi4;
    double bt1_mev,th,thl[70],ds[70],ay[70],ap1[70],ap2[70],t21,ap3[70];
    double sig[20][20][15][200],ayn[20][20][15][200],ayd[20][20][15][200],
	  axx[20][20][15][200],ayy[20][20][15][200],axn[20][20][15][200],
	  axy[20][20][15][200],axz[20][20][15][200];

    int aj1bx, aj1by, ajbz, aj1theta, aj2theta, aj1phi, aj2phi, aj1ekin, aja10,
    	aja20, aju, ajcs1, ajcs02, aj1r1, aj1r2, npd_choice, icros, n_bar, n_pipe;
    double bfwhmx, bfwhmy, eps_separ, bt, pz, pzz, themin, themax, themin2,
     themax2, fimin,fimax, dim_mars[3],z_shift, thigh, tXplace,tYplace,tZplace,
      dim_mwpc1[5], dim_ring[2], xyz_mwpc1[3], dim_mwpc2[4], dim_hole2[2],
      xyz_mwpc2[3], par_mwpc[5], dim_delta[4], xyz_delta[6], xyz_degrader[3] ,
      exit_win_z, exit_win_th, exit_win_rad, z_pipe, dim_degrader[3] ;
    double generator_min, generator_max;
    G4String FileName;
    double tes1[4], tes2[4], tes3[4], tes4[4];
    int tempGetChoice;
     int tempProcNb;
  public:
    void GeneratePrimaries(G4Event*);

     double* GetStartEnergy (double en1 = -1., double en2 = -1., double en3 = -1., double en4 = -1.)
    {

      if (en1 != -1.)
      {
        tes1[0] = en1;
	tes1[1] = en2;
	tes1[2] = en3;
	tes1[3] = en4;
      }
      return tes1;
    };

     double* GetStartAngleTheta (double ang1 = -1., double ang2 = -1., double ang3 = -1., double ang4 = -1.)
    {

      if (ang1 != -1.)
      {
        tes2[0] = ang1;
	tes2[1] = ang2;
	tes2[2] = ang3;
	tes2[3] = ang4;
      }
      return tes2;
    };

     double* GetStartAnglePhi (double ang1 = -1., double ang2 = -1., double ang3 = -1., double ang4 = -1.)
    {

      if (ang1 != -1.)
      {
        tes3[0] = ang1;
	tes3[1] = ang2;
	tes3[2] = ang3;
	tes3[3] = ang4;
      }
      return tes3;
    };

     double* GetStartPosition (double x1 = 0., double x2 = 0., double x3 = 0.)
    {

          if ((x1 != 0.)||(x2 != 0.)||(x3 != 0.))
      {
        tes4[0] = x1;	tes4[1] = x2;	tes4[2] = x3;
      }
      return tes4;
    }

     int ProcNb(int num = 10)
    {

      if (num != 10) tempProcNb = num;
      return tempProcNb;
    };

     int GetChoice (int num = 10)
    {

      if (num != 10) tempGetChoice = num;
      return tempGetChoice;
    };


  private:
    G4LorentzVector v1,v2,v3,v4;
    std::ifstream file_Pluto_generator;
    G4ParticleGun* particleGun1;		//proton 1
    G4ParticleGun* particleGun2;		//proton 2
    G4ParticleGun* particleGun3;		//neutron
    G4ParticleGun* particleGun4;		//neutron
    G4ParticleGun* event_cleaner_particle_gun;
    Bina_DetectorConstruction* myDetector;

    void RandomInit(int =2);			//generators initialization
    double RandomGauss(double, double =0, double =1);	// gauss
    double RandomFlat (double, double =0, double =1);	// flat
    void Pos(void);				//generate vertex position

    void read_part_momentum(double*);
    void read_part_momentum4(double*);
    void open_pluto_file();

    static MyFileReader* fileReader;





    inline double P(double EE,double MP){return sqrt(EE*EE - MP*MP);} 	//P_relatywist. [MeV] (E_relatywist. [MeV] )
    inline double E(double PP,double MP){return sqrt(PP*PP + MP*MP);}  	//E_relatywist. [MeV] (P_relatywist. [MeV] )


    double momentum[9]; 		// table whit track momentum
    double vertex[3];			// table whith vertex position
    double p_mass, d_mass, n_mass;	// pointer to particle mass


   // G4RandGauss* GaussDist;//[2];		//table with diferent generators with gaussian distribution
   // RandGaurs *GaD = new RandGauss[2];
    CLHEP::RandGauss* GaussDist[2];
    CLHEP:: RandFlat* FlatDist[9];
};


#endif
