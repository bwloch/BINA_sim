#ifndef Bina_DetectorConstruction_h
#define Bina_DetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4RotationMatrix.hh"

class G4Trd;
class G4Trap;
class G4Box;
class G4Sphere;
class G4Tubs;
class G4Polyhedra;
class G4Polycone;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4SubtractionSolid;
class G4UnionSolid;
class G4IntersectionSolid;
class G4Material;
class G4VSolid;
class Bina_DetectorMessenger;
class Bina_DetectorConstruction;

class Bina_DetectorConstruction : public G4VUserDetectorConstruction
{

  public:

     Bina_DetectorConstruction();
    ~Bina_DetectorConstruction();

  public:

     G4VPhysicalVolume* Construct();

//////////////////////////////
// Generation Parameters
//////////////////////////////

public:
  void SetNpdChoice(G4int wybor)       {npd_choice=wybor;};
  void SetNeumann(G4int Neumann)       {neumann=Neumann;};
  void SetBfwhmX(G4double BfwhmX)      {bfwhmx=BfwhmX;};
  void SetBfwhmY(G4double BfwhmY)      {bfwhmy=BfwhmY;};
  void SetBtEnergy(G4double BtEnergy)  {bt=BtEnergy;};
  void SetPz(G4double Pz)              {pz=Pz;};
  void SetPzz(G4double Pzz)            {pzz=Pzz;};
  void SetThetaMin(G4double ThetaMin)  {themin=ThetaMin;};
  void SetThetaMax(G4double ThetaMax)  {themax=ThetaMax;};
  void SetTheta2Min(G4double The2Min)  {themin2=The2Min;};
  void SetTheta2Max(G4double Theta2Max){themax2=Theta2Max;};
  void SetPhiMin(G4double PhiMin)      {fimin=PhiMin;};
  void SetPhiMax(G4double PhiMax)      {fimax=PhiMax;};
  void SetKinematicsMin(G4double kinMin) {genMin=kinMin*0.001;};
  void SetKinematicsMax(G4double kinMax) {genMax=kinMax*0.001;};
//  void SetParamUpdate() ;  

public:
  int    GetNpdChoice()   {return npd_choice;};  
  int    GetNeumann()     {return neumann;};
  double GetBfwhmX()      {return bfwhmx;};
  double GetBfwhmY()      {return bfwhmy;};
  double GetBtEnergy()    {return bt;};
  double GetPz()          {return pz;};
  double GetPzz()         {return pzz;};
  double GetThetaMin()    {return themin;};
  double GetThetaMax()    {return themax;};
  double GetTheta2Min()   {return themin2;};
  double GetTheta2Max()   {return themax2;};
  double GetPhiMin()      {return fimin;};
  double GetPhiMax()      {return fimax;};
  double GetKinematicsMin() {return genMin;};
  double GetKinematicsMax() {return genMax;};
  double GetTargetHigh()  {return thigh;};
  double GetTargetXplace(){return tXplace;};
  double GetTargetYplace(){return tYplace;};
  double GetTargetZplace(){return tZplace;};
//////////////////////////////
// Target
//////////////////////////////
public:
  void SetTargetIs(G4int TargetIs)                 { targetis = TargetIs;};
  void SetTargetOutRadius(G4double HoleOut)        { toutrad = HoleOut; };
  void SetTframeOutRadius(G4double FrameOut)       { tfroutrad = FrameOut; };
  void SetTargetHigh(G4double TargetHigh)          { thigh = TargetHigh; };
// Target placement
  void SetTargetXplace(G4double TargetXplace)      { tXplace = TargetXplace; };
  void SetTargetYplace(G4double TargetYplace)      { tYplace = TargetYplace; };
  void SetTargetZplace(G4double TargetZplace)      { tZplace = TargetZplace; };
private:
  G4int targetis;
  G4double toutrad,  tfroutrad, thigh, tXplace, tYplace, tZplace;
//////////////////////////////
//  Ball scintilators
//////////////////////////////
public:  
  void SetBallIs(G4int BallIs)               {ballis=BallIs;};
  void SetBallVisible(G4int BallVisible)     {ballvisible=BallVisible;};
//Ball placement
  void SetBallXplace(G4double BallXplace)    {ballXplace=BallXplace;};
  void SetBallYplace(G4double BallYplace)    {ballYplace=BallYplace;};
  void SetBallZplace(G4double BallZplace)    {ballZplace=BallZplace;};
private:
G4int ballis, ballvisible; 
G4double ballXplace,ballYplace,ballZplace;
//////////////////////////////
//  Mwpc
//////////////////////////////
public:
  // MWPC
  void SetMwpcIs(G4int MwpcIs) { mwpcis=MwpcIs;};
  void SetMwpcDimFrameX(G4double MwpcDimFrameX){ mwpcdimfx=MwpcDimFrameX; };
  void SetMwpcDimFrameZ(G4double MwpcDimFrameZ){ mwpcdimfz=MwpcDimFrameZ; };
  void SetMwpcDimGasX(G4double MwpcDimGasX)    { mwpcdimgx=MwpcDimGasX; };
  void SetMwpcDimGasZ(G4double MwpcDimGasZ)    { mwpcdimgz=MwpcDimGasZ; };
  void SetMwpcHoleIn(G4double MwpcHoleIn)      { mwpcholeIn=MwpcHoleIn; };
  void SetMwpcHoleOut(G4double MwpcHoleOut)    { mwpcholeOut=MwpcHoleOut; };
// Mwpc placement
  void SetMwpcXplace(G4double MwpcXplace)      { mwpcXplace=MwpcXplace; };
  void SetMwpcYplace(G4double MwpcYplace)      { mwpcYplace=MwpcYplace; };
  void SetMwpcZplace(G4double MwpcZplace)      { mwpcZplace=MwpcZplace; };
private:
  G4int mwpcis;
  G4double mwpcdimfx, mwpcdimfz, mwpcdimgx, mwpcdimgz, mwpcholeIn, mwpcholeOut,
           mwpcXplace, mwpcYplace, mwpcZplace ;
  
//////////////////////////////
// Delta E
//////////////////////////////
public:
  void SetDeltaeIs(G4int DeltaeIs)             {deIs=DeltaeIs;};
  void SetDeltaeVisible(G4int DeltaeVisible)   {deVisible=DeltaeVisible;};
  void SetDeltaeNMax(G4int DeltaeNmax)         {nde_max=DeltaeNmax;};   
  void SetDeltaeN0(G4int DeltaeNmin)           {nde_min=DeltaeNmin;}; 
  void SetDeltaeFoil(G4double DeltaeFoil)      {dedimf=DeltaeFoil;};
  void SetDeltaeDimX(G4double DeltaeDimX)      {dedimX=DeltaeDimX;};
  void SetDeltaeDimY(G4double DeltaeDimY)      {dedimY=DeltaeDimY;};
  void SetDeltaeDimZ(G4double DeltaeDimZ)      {dedimZ=DeltaeDimZ;};
  void SetDeltaeDimf(G4double DeltaeDimf)      {dedimf=DeltaeDimf;};
  void SetDeltaeHoleOut(G4double HoleOut)      {deHoleOut=HoleOut;};
  void SetDeltaeSepar(G4double DeltaeSepar)    {deSepar=DeltaeSepar;};
// Delta E placement
  void SetDeltaeXplace(G4double DeltaeXplace)  {deXplace=DeltaeXplace;};
  void SetDeltaeYplace(G4double DeltaeYplace)  {deYplace=DeltaeYplace;};
  void SetDeltaeZplace(G4double DeltaeZplace)  {deZplace=DeltaeZplace;};
private:
  G4int deIs, deVisible, nde_max, nde_min;
  G4double deXplace, deYplace, deZplace,dedimX, dedimY, dedimZ, dedimf;
  G4double deHoleOut,  deSepar;
///////////////////////////////
//  Foil
///////////////////////////////
public:
  void SetFoilIs(G4int FoilIs)              {foilIs=FoilIs;};
  void SetFoilVisible(G4int FoilVisible)    {foilVisible=FoilVisible;};
  void SetFoilDimX(G4double FoilDimX)       {foildimX=FoilDimX;};
  void SetFoilDimY(G4double FoilDimY)       {foildimY=FoilDimY;};
  void SetFoilDimZ(G4double FoilDimZ)       {foildimZ=FoilDimZ;};
  void SetFoilHoleIn(G4double FoilHoleIn)   {foilholeIn=FoilHoleIn;};
  void SetFoilHoleOut(G4double HoleOut)     {foilholeOut=HoleOut;};
// Foil placement
  void SetFoilXplace(G4double FoilXplace)   {foilXplace=FoilXplace;};
  void SetFoilYplace(G4double FoilYplace)   {foilYplace=FoilYplace;};
  void SetFoilZplace(G4double FoilZplace)   {foilZplace=FoilZplace;};
private:
  G4int foilIs, foilVisible;
  G4double foildimX,foildimY,foildimZ,foilholeIn,foilholeOut;
  G4double foilXplace,foilYplace,foilZplace;
//////////////////////////////////////
// Wall E Center scintilators
//////////////////////////////////////
public:
 void SetWallIs(G4int WallIs) {wallIs=WallIs;};
 void SetWallVisible(G4int WallVisible) {wallVisible=WallVisible;};

 void SetWallRIn(G4double WallRIn)         {wallRIn=WallRIn;};
 void SetWallROut(G4double WallROut)       {wallROut=WallROut;};
 void SetWallDimCX(G4double WallDimX)      {wallDimcX=WallDimX;};
 void SetWallAng0(G4double WallAng0)       {wallAng0=WallAng0;};
 void SetWallAng(G4double WallAng)         {wallAng=WallAng;};

 void SetWallHoleIn(G4double WallHoleIn)   {wallholeIn=WallHoleIn;};
 void SetWallHoleOut(G4double WallHoleOut) {wallholeOut=WallHoleOut;};
 void SetWallSepar(G4double WallSepar)     {wallSepar=WallSepar;};

 void SetWallN0(G4int WallN0)              {wallN0=WallN0;};
 void SetWallNCenter(G4int WallNCenter)    {wallNCenter=WallNCenter;};
// Wall placement
 void SetWallCenterXplace(G4double WallXplace) {wallXplace=WallXplace;};
 void SetWallCenterYplace(G4double WallYplace) {wallYplace=WallYplace;};
 void SetWallCenterZplace(G4double WallZplace) {wallZplace=WallZplace;};
private: 
 G4int wallIs, wallVisible, wallN0, wallNCenter, wallNMax;
 G4double wallRIn, wallROut, wallDimcX, wallAng0, wallAng, wallholeIn, wallholeOut, wallSepar;
 G4double wallXplace,wallYplace,wallZplace;
//////////////////////////////
// Wall E L scintilators
//////////////////////////////
public:
 void SetWallDimX(G4double WallLDimX)     {wallLdimX=WallLDimX;};
 void SetWallDimY(G4double WallLDimY)     {wallLdimY=WallLDimY;};
 void SetWallDimZ(G4double WallLDimZ)     {wallLdimZ=WallLDimZ;};
 void SetWallNMax(G4int WallNMax)          {wallNMax=WallNMax;};
// Wall placement
 void SetWallDimS(G4double WallDimS) {walldimS=WallDimS;};

private:
 G4double wallLdimX, wallLdimY, wallLdimZ, walldimS;

 // private:
    G4int aj1bx, aj1by, ajbz, aj1theta, aj2theta, aj1phi, aj2phi, aj1ekin, aja10,
    	aja20, aju, ajcs1, ajcs02, aj1r1, aj1r2, npd_choice, neumann, n_elastic, n_bar, n_pipe;
    G4double bfwhmx, bfwhmy, eps_separ, bt, pz, pzz, themin, themax, themin2,
    	   themax2, fimin, fimax, genMin, genMax;
    double dim_mars[3], z_shift, dim_trgt[3], xyz_trgt[3],
      dim_mwpc1[5], dim_ring[2], xyz_mwpc1[3], dim_mwpc2[4], dim_hole2[2],
      xyz_mwpc2[3], par_mwpc[5], dim_delta[4], xyz_delta[6], xyz_degrader[3],
      exit_win_z, exit_win_th, exit_win_rad, z_pipe, dim_degrader[3];

     double *w_eps_separ, *w_z_shift, *w_z_pipe;
     int *w_n_bar, *w_n_pipe;

  private:
     void Seen(void);			//visualization parameters
     void Renumerate(void);		//renumerate delta E and Salad to experimental numeration

     G4Box*             Mars_sol;    // pointer to the solid Mars -> World
     G4LogicalVolume*   Mars_log;    // pointer to the logical Mars -> World
     G4VPhysicalVolume* Mars_phs;    // pointer to the physical Mars -> World

     G4Box*             Mars_vac_sol;  // pointer to the solid Mars -> Vacuum
     G4LogicalVolume*   Mars_vac_log;  // pointer to the logical Mars -> Vacuum
     G4VPhysicalVolume* Mars_vac_phs; // pointer to the physical Mars -> Vacuum

     G4Tubs*            trgt_sol;   // pointer to the solid Target -> liguid hydrogen
     G4LogicalVolume*   trgt_log;   // pointer to the logical Target -> liguid hydrogen
     G4VPhysicalVolume* trgt_phs;   // pointer to the physical Target -> liguid hydrogen
     G4Tubs*            trgt_fr_sol;   // targer frame
     G4LogicalVolume*   trgt_fr_log;   //
     G4VPhysicalVolume* trgt_fr_phs;   // 

     G4Tubs*            trgt_sh_sol;   // targer frame
     G4LogicalVolume*   trgt_sh_log;   //
     G4VPhysicalVolume* trgt_sh_phs;   // 
     G4SubtractionSolid* trgt_shc_sol;

     G4Tubs*            trgt_win_sol;   // pointer to the solid Target window -> aramica 25um
     G4LogicalVolume*   trgt_win_log;   // pointer to the logical Target window -> aramica 25um
     G4VPhysicalVolume* trgt_win_phs[2];   // pointer to the physical Target window -> aramica 25um

     G4Tubs*            exit_win_sol;   // pointer to the solid 
     G4LogicalVolume*   exit_win_log;   // pointer to the logical 
     G4VPhysicalVolume* exit_win_phs;   // pointer to the physical

     G4Box*             mwpc1_sol;  // pointer to the solid MWPC1 - > mother
     G4SubtractionSolid* mwpc1_minus_sol;
     G4LogicalVolume*   mwpc1_log;  // pointer to the logical MWPC1 - > mother
     G4VPhysicalVolume* mwpc1_phs;  // pointer to the physical MWPC1 - > mother

     G4Polyhedra*       mwpc1_frame_sol;  // pointer to the solid MWPC1 frame
     G4LogicalVolume*   mwpc1_frame_log;  // pointer to the logical MWPC1 frame
     G4VPhysicalVolume* mwpc1_frame_phs;  // pointer to the physical MWPC1 frame

     G4Box*             mwpc1_gas_sol;  // pointer to the solid MWPC1 gas chamber
     G4SubtractionSolid* mwpc1_gas_minus_sol;
     G4LogicalVolume*   mwpc1_gas_log;  // pointer to the logical MWPC1 gas chamber
     G4VPhysicalVolume* mwpc1_gas_phs;  // pointer to the physical MWPC1 gas chamber

     G4Box*		mwpc1_wire_plane_mother_sol;
     G4Tubs*   		mwpc1_hole_sol;
     G4SubtractionSolid* mwpc1_wire_plane_mother_minus_sol;
     G4LogicalVolume*	mwpc1_wire_plane_mother_log[12];
     G4Tubs*            mwpc1_wire_sol;
     G4Tubs*            mwpc1_wire_c_sol[30];
     G4Tubs*            mwpc1_wire_45_sol[150];
     G4LogicalVolume*   mwpc1_wire_c_log[30];
     G4LogicalVolume*   mwpc1_wire_45_log[150];
     G4LogicalVolume*   mwpc1_wire_log;
     G4VPhysicalVolume* mwpc1_wire_plane_mother_phs[12];
     G4VPhysicalVolume* mwpc1_wire_x_phs[250];
     G4VPhysicalVolume* mwpc1_wire_y_phs[250];
     G4VPhysicalVolume* mwpc1_wire_45_phs[350];

     G4Box*		mwpc1_catode_foil_sol;
     G4SubtractionSolid* mwpc1_catode_minus_sol;
     G4LogicalVolume*	mwpc1_catode_foil_log;
     G4VPhysicalVolume* mwpc1_catode_foil_phs[8];

     G4Tubs*            mwpc1_ring_sol;  // pointer to the solid MWPC1 vacum hole  (G4double center)
     G4SubtractionSolid* mwpc1_ring_minus_sol;
     G4LogicalVolume*   mwpc1_ring_log;  // pointer to the logical MWPC1 vacum hole (G4double center)
     G4VPhysicalVolume* mwpc1_ring_phs;  // pointer to the physical MWPC1 vacum hole (G4double center)

     G4Box*             mwpc1_win_sol;  // pointer to the solid MWPC1 window
     G4LogicalVolume*   mwpc1_win_log[2];  // pointer to the logical MWPC1 window
     G4SubtractionSolid* mwpc1_win_minus_sol[2];
     G4VPhysicalVolume* mwpc1_win_phs[2];  // pointer to the physical MWPC1 window right


     G4Polycone*        pipe_sol;
     G4LogicalVolume*   pipe_log;  // pointer to the logical
     G4VPhysicalVolume* pipe_phs;  // pointer to the physical

     G4Tubs*		win_scatering_sol;
     G4LogicalVolume*   win_scatering_log;
     G4VPhysicalVolume* win_scatering_phs;


     
     G4Tubs*            walec;
     G4Box*             deltaE_w_sol;
     G4SubtractionSolid* deltaE_wminus_sol[10];
     G4LogicalVolume*   deltaE_w_log[3];
     G4VPhysicalVolume* deltaE_w_phs[40];

///////////////////////
///  Bina E, dE
///////////////////////
     G4Sphere*          Bina_sol;
     G4LogicalVolume*   Bina_log;
     G4SubtractionSolid*Bina_minus;
     G4VPhysicalVolume* Bina_phs;

     G4Box*    blok_sol;
     G4VSolid*          Edet_B_sol[7];
     G4LogicalVolume*   Edet_B_log[7];
     G4VPhysicalVolume* Edet_B_phs[200];
     G4VSolid*          deltaE_B_sol[7];
     G4LogicalVolume*   deltaE_B_log[7];
     G4VPhysicalVolume* deltaE_B_phs[200];

///////////////////////
///  Delta E
///////////////////////
     
     G4Box*              deltae_W_sol;
     G4SubtractionSolid* deltae_Wc_sol[10];
     G4LogicalVolume*    deltae_W_log[10];
     G4VPhysicalVolume*  deltae_W_phs[40];
     G4Box*              deltae_foil_sol;
     G4SubtractionSolid* deltae_foilc_sol[10];
     G4LogicalVolume*    deltae_foil_log[10];
     G4VPhysicalVolume*  deltae_foil_phs[40];

///////////////////////
///  Folia
///////////////////////
     
     G4Box*              foil_sol;
     G4SubtractionSolid* foil_min_sol;
     G4LogicalVolume*    foil_log;
     G4VPhysicalVolume*  foil_phs;
     
///////////////////////
///  Wall E
///////////////////////
     G4Box*              wall_L_sol;
     G4Box*              wallf_L_sol;
     G4Polyhedra*        wall_sol;
     G4Polyhedra*        wallf_sol;
     G4Trap*             wall_Spec_sol;
     G4Trap*             wallf_Spec_sol;
     G4SubtractionSolid* wall_C_sol;
     G4SubtractionSolid* wallf_C_sol;
     G4LogicalVolume*    wall_log[10];
     G4LogicalVolume*    wallf_log[10];
     G4LogicalVolume*    wall_C_log[2];
     G4LogicalVolume*    wallf_C_log[2];
     G4LogicalVolume*    wall_L_log[10];
     G4LogicalVolume*    wallf_L_log[10];
     G4LogicalVolume*    wall_Spec_log[2];
     G4LogicalVolume*    wallf_Spec_log[2];
     G4VPhysicalVolume*  wall_phs[20];  
     G4VPhysicalVolume*  wallf_phs[20];  
     
     Bina_DetectorMessenger* detectorMessenger;  // pointer to the Messenger
};
#endif

