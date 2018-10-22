#ifndef Bina_PrimaryGeneratorAction_h
#define Bina_PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"

#include "CLHEP/Random/Randomize.h"
#include "CLHEP/Random/RandomEngine.h"
//#include "CLHEP/Random/RandGauss.h"
//#include "CLHEP/Random/RandGaussQ.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/Ranlux64Engine.h"


class Bina_DetectorConstruction;
class G4ParticleGun;
class G4Event;
//class G4RandGauss;
class Zmienne;
/*class RandFlat;
class Ranlux64Engine: public HepRandomEngine
{
	public:
		Ranlux64Engine();
		virtual ~Ranlux64Engine();};
class RandFlat; : public HepRandom //tto jest niepotrzebne
{
  public:
        inline RandFlat ( HepRandomEngine* anEngine );
	virtual~ RandFlat();
};*/
class Bina_PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    Bina_PrimaryGeneratorAction(Bina_DetectorConstruction*);
   ~Bina_PrimaryGeneratorAction();
//////////////////////////////
// Generation Parameters
//////////////////////////////
    double bt1,bt2,bt3,t1,t2,t3,fi1,fi2,fi3;
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
  public:
    void GeneratePrimaries(G4Event*);

    inline static double* GetStartEnergy (double en1 = -1., double en2 = -1., double en3 = -1.)
    {
      static double tes1[3];
      if (en1 != -1.)
      {
        tes1[0] = en1;
	tes1[1] = en2;
	tes1[2] = en3;
      }
      return tes1;
    };

    inline static double* GetStartAngleTheta (double ang1 = -1., double ang2 = -1., double ang3 = -1.)
    {
      static double tes2[3];
      if (ang1 != -1.)
      {
        tes2[0] = ang1;
	tes2[1] = ang2;
	tes2[2] = ang3;
      }
      return tes2;
    };

    inline static double* GetStartAnglePhi (double ang1 = -1., double ang2 = -1., double ang3 = -1.)
    {
      static double tes3[3];
      if (ang1 != -1.)
      {
        tes3[0] = ang1;
	tes3[1] = ang2;
	tes3[2] = ang3;
      }
      return tes3;
    };

    inline static double* GetStartPosition (double v1 = 0., double v2 = 0., double v3 = 0.)
    {
      static double tes3[3];
          if ((v1 != 0.)||(v2 != 0.)||(v3 != 0.))
      {
        tes3[0] = v1;	tes3[1] = v2;	tes3[2] = v3;
      }
      return tes3;
    }

    inline static int ProcNb(int num = 10)
    {
      static int temp;
      if (num != 10) temp = num;
      return temp;
    };

    inline static int GetChoice (int num = 10)
    {
      static int temp;
      if (num != 10) temp = num;
      return temp;
    };

  private:
    G4ParticleGun* particleGun1;		//proton 1
    G4ParticleGun* particleGun2;		//proton 2
    G4ParticleGun* particleGun3;		//neutron
    Bina_DetectorConstruction* myDetector;

    void RandomInit(int =2);			//generators initialization
    double RandomGauss(double, double =0, double =1);	// gauss
    double RandomFlat (double, double =0, double =1);	// flat
    void Pos(void);				//generate vertex position

    double* elastic(double*);
    double* ugbreak(double*);
    double* ugelast(double*);
    double* upunif(double*);

    void ugelast_read(void);
    void break_read(void);

    double gelkin(double,double*,double*);  // elastic kinematics
    double rinterp(double [5][5][5][5], double, double, double, double);
    double cspol(double,double,double,double,double);
    double dpb_kin3(const int,double,double*,double,double*,double*,double*);
    double *vector(long,long);
    void free_vector(double*,long,long);
    void carsph(double*,double*,double*,double*);
    void splint(double*,double*,double*,int,double,double*);
    void spline(double*,double*,int,double,double,double*);
    void spline4(double*,int,int,int,int,double*);
    void splint4(double*,double*,int,int,int,int,double,double,double,double,double*);

    inline double P(double EE,double MP){return sqrt(EE*EE - MP*MP);} 	//P_relatywist. [MeV] (E_relatywist. [MeV] )
    inline double E(double PP,double MP){return sqrt(PP*PP + MP*MP);}  	//E_relatywist. [MeV] (P_relatywist. [MeV] )

    // transformation from s1 do s
    inline double Z(double xx1,double ,double zz1,double alfa) {return (zz1*cos(alfa) - xx1*sin(alfa));}
    inline double X(double xx1,double ,double zz1,double alfa) {return (zz1*sin(alfa) + xx1*cos(alfa));}
    inline double Y(double ,double yy1,double ,double ) {return yy1;}
    double momentum[9]; 		// table whit track momentum
    double vertex[3];			// table whith vertex position
    double p_mass, d_mass, n_mass;	// pointer to particle mass

    double x1a[4],x2a[4],x3a[4],x4a[4];
    double x1a_r[4],x2a_r[4],x3a_r[4],x4a_r[4];
    double ya_r[5][5][5][5], yax_r[5][5][5][5];
    double yay_r[5][5][5][5],yaxx_r[5][5][5][5],yaxy_r[5][5][5][5],yayy_r[5][5][5][5];
   // G4RandGauss* GaussDist;//[2];		//table with diferent generators with gaussian distribution
   // RandGaurs *GaD = new RandGauss[2];
    CLHEP::RandFlat* GaussDist[2];
    CLHEP:: RandFlat* FlatDist[9];
};


#endif


