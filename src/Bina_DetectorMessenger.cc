#include "Bina_DetectorMessenger.hh"

#include "Bina_DetectorConstruction.hh"

#include "G4UIdirectory.hh"

#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithoutParameter.hh"

#include "G4UIcmdWithADoubleAndUnit.hh"


Bina_DetectorMessenger::Bina_DetectorMessenger(Bina_DetectorConstruction* Det)
:myDetector(Det)
{


//////////////////////////////
// Directory
//////////////////////////////

  geomDir = new G4UIdirectory("/geom/");
  geomDir->SetGuidance("Detector control.");
  
  paramDir = new G4UIdirectory("/param/");
  paramDir->SetGuidance("Detector control.");
  
  generationDir=new G4UIdirectory("/generator/");
  generationDir->SetGuidance("Breakup generator settings");

  GeometryTargetDir =new G4UIdirectory("/geom/target/");
  GeometryTargetDir->SetGuidance("Target sizes control.");

  GometryBallDir = new G4UIdirectory("/geom/ball/");
  GometryBallDir->SetGuidance("BINA Ball sizes control.");

  GeometryDeltaeDir = new G4UIdirectory("/geom/deltae/");
  GeometryDeltaeDir->SetGuidance("Delta E sizes control");

  GeometryFoilDir = new G4UIdirectory("/geom/foil/");
  GeometryFoilDir->SetGuidance("foil sizes control");

  GometryWallDir = new G4UIdirectory("/geom/wall/");
  GometryWallDir->SetGuidance("Wall sizes control");

  GeometryMwpcDir = new G4UIdirectory("/geom/mwpc/");
  GeometryMwpcDir->SetGuidance("Multi Wire Proportional Chamber control.");


//////////////////////////////
// Generation Parameters
//////////////////////////////

  NpdChoiceCmd = new G4UIcmdWithAnInteger("/param/npd_choice",this);
  NpdChoiceCmd->SetGuidance("Choice of process: 1 - elast. dp, 2 - break-up, 0 - elast uniform., -1 - p.");

  IfNeumannCmd = new G4UIcmdWithAnInteger("/param/neumann",this);
  IfNeumannCmd->SetGuidance("Only for npd_choice==2, Sets the method of generation \n 0-uniform, 1-using cross section tables");

  BfwhmX_Cmd = new G4UIcmdWithADouble("/param/bfwhmx",this);
  BfwhmX_Cmd->SetGuidance("Beam FWHMs x.");

  BfwhmY_Cmd = new G4UIcmdWithADouble("/param/bfwhmy",this);
  BfwhmY_Cmd->SetGuidance("Beam FWHMs y.");

  BtEnergyCmd = new G4UIcmdWithADoubleAndUnit("/param/bt",this);
  BtEnergyCmd->SetGuidance("Beam energy in MeV.");
  BtEnergyCmd->SetDefaultUnit("MeV");

  Pz_Cmd = new G4UIcmdWithADouble("/param/pz",this);
  Pz_Cmd->SetGuidance("Polarization Pz.");

  Pzz_Cmd = new G4UIcmdWithADouble("/param/pzz",this);
  Pzz_Cmd->SetGuidance("Polarization Pzz.");

  ThetaMinCmd = new G4UIcmdWithADoubleAndUnit("/param/themin",this);
  ThetaMinCmd->SetGuidance("Thetata generation min.");
  ThetaMinCmd->SetDefaultUnit("deg");

  ThetaMaxCmd = new G4UIcmdWithADoubleAndUnit("/param/themax",this);
  ThetaMaxCmd->SetGuidance("Thetata generation max.");
  ThetaMaxCmd->SetDefaultUnit("deg");

  Theta2MinCmd = new G4UIcmdWithADoubleAndUnit("/param/the2min",this);
  Theta2MinCmd->SetGuidance("Thetata 2 generation min.");
  Theta2MinCmd->SetDefaultUnit("deg");

  Theta2MaxCmd = new G4UIcmdWithADoubleAndUnit("/param/the2max",this);
  Theta2MaxCmd->SetGuidance("Thetata 2 generation max.");
  Theta2MaxCmd->SetDefaultUnit("deg");

  PhiMinCmd = new G4UIcmdWithADoubleAndUnit("/param/phimin",this);
  PhiMinCmd->SetGuidance("Phi generation min.");
  PhiMinCmd->SetDefaultUnit("deg");

  PhiMaxCmd = new G4UIcmdWithADoubleAndUnit("/param/phimax",this);
  PhiMaxCmd->SetGuidance("Phi generation max.");
  PhiMaxCmd->SetDefaultUnit("deg");
  
//  ParamUpdate = new G4UIcmdWithoutParameter("/param/update",this);

//////////////////////////////
// Target
//////////////////////////////

  TargetIsCmd = new G4UIcmdWithAnInteger("/geom/target/is",this);
  TargetIsCmd->SetGuidance("If  Target is =1 , if not =0");
  TargetIsCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  TargetOutRadiusCmd = new G4UIcmdWithADoubleAndUnit("/geom/target/radius",this);
  TargetOutRadiusCmd->SetGuidance("Set Half-dimensions radius of the Target");
  TargetOutRadiusCmd->SetDefaultUnit("cm");

  TframeOutRadiusCmd = new G4UIcmdWithADoubleAndUnit("/geom/target/frame",this);
  TframeOutRadiusCmd->SetGuidance("Set outer radius of the Target frame");
  TframeOutRadiusCmd->SetDefaultUnit("cm");

  TargetHighCmd = new G4UIcmdWithADoubleAndUnit("/geom/target/high",this);
  TargetHighCmd->SetGuidance("Set Half-dimensions Z-high of the Target");

  TargetHighCmd->SetDefaultUnit("cm");


// Target placement
  TargetXplaceCmd = new G4UIcmdWithADoubleAndUnit("/geom/target/place_x",this);
  TargetXplaceCmd->SetGuidance("Set placement of the Target in X axis");

  TargetXplaceCmd->SetDefaultUnit("cm");


  TargetYplaceCmd = new G4UIcmdWithADoubleAndUnit("/geom/target/place_y",this);
  TargetYplaceCmd->SetGuidance("Set placement of the Target in Y axis");

  TargetYplaceCmd->SetDefaultUnit("cm");


  TargetZplaceCmd = new G4UIcmdWithADoubleAndUnit("/geom/target/place_z",this);
  TargetZplaceCmd->SetGuidance("Set placement of the Target in Z axis");

  TargetZplaceCmd->SetDefaultUnit("cm");


//////////////////////////////
//  Ball scintilators
//////////////////////////////

  BallIsCmd = new G4UIcmdWithAnInteger("/geom/ball/is",this);
  BallIsCmd->SetGuidance("If  Ball is =1 , if not =0");

  BallVisibleCmd = new G4UIcmdWithAnInteger("/geom/ball/visible",this);
  BallVisibleCmd->SetGuidance("If Ball visible =1 , if not =0");

  BallXplaceCmd = new G4UIcmdWithADoubleAndUnit("/geom/ball/place_x",this);
  BallXplaceCmd->SetGuidance("Set placement of the Ball in X axis");
  BallXplaceCmd->SetDefaultUnit("cm");

  BallYplaceCmd = new G4UIcmdWithADoubleAndUnit("/geom/ball/place_y",this);
  BallYplaceCmd->SetGuidance("Set placement of the Ball in Y axis");
  BallYplaceCmd->SetDefaultUnit("cm");

  BallZplaceCmd = new G4UIcmdWithADoubleAndUnit("/geom/ball/place_z",this);
  BallZplaceCmd->SetGuidance("Set placement of the Ball in Z axis");
  BallZplaceCmd->SetDefaultUnit("cm");

//////////////////////////////
//  Mwpc
//////////////////////////////

  MwpcIsCmd = new G4UIcmdWithAnInteger("/geom/mwpc/is",this);
  MwpcIsCmd->SetGuidance("If  Mwpc is =1 , if not =0");
  MwpcIsCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  MwpcDimFrameXCmd = new G4UIcmdWithADoubleAndUnit("/geom/mwpc/dim_x",this);
  MwpcDimFrameXCmd->SetGuidance("Half-dimensions of MWPC1 Frame corners: dx/2=dy/2");
  MwpcDimFrameXCmd->SetDefaultUnit("cm");

  MwpcDimFrameZCmd = new G4UIcmdWithADoubleAndUnit("/geom/mwpc/dim_z",this);
  MwpcDimFrameZCmd->SetGuidance("Half-dimensions of MWPC1 Frame dz/2");
  MwpcDimFrameZCmd->SetDefaultUnit("cm");

  MwpcDimGasXCmd = new G4UIcmdWithADoubleAndUnit("/geom/mwpc/dimGas_x",this);
  MwpcDimGasXCmd->SetGuidance("Half-dimensions of MWPC1 Gas corners: dx/2=dy/2");
  MwpcDimGasXCmd->SetDefaultUnit("cm");

  MwpcDimGasZCmd = new G4UIcmdWithADoubleAndUnit("/geom/mwpc/dimGas_z",this);
  MwpcDimGasZCmd->SetGuidance("Half-dimensions of MWPC1 Gas dz/2");
  MwpcDimGasZCmd->SetDefaultUnit("cm");

  MwpcHoleInCmd = new G4UIcmdWithADoubleAndUnit("/geom/mwpc/hole_rIn",this);
  MwpcHoleInCmd->SetGuidance("Hole in MWPC1 : r_in");
  MwpcHoleInCmd->SetDefaultUnit("cm");

  MwpcHoleOutCmd = new G4UIcmdWithADoubleAndUnit("/geom/mwpc/hole_rOut",this);
  MwpcHoleOutCmd->SetGuidance("Hole in MWPC1 : r_out");
  MwpcHoleOutCmd->SetDefaultUnit("cm");

  MwpcXplaceCmd = new G4UIcmdWithADoubleAndUnit("/geom/mwpc/place_x",this);
  MwpcXplaceCmd->SetGuidance("Set placement of the Mwpc in X axis");
  MwpcXplaceCmd->SetDefaultUnit("cm");

  MwpcYplaceCmd = new G4UIcmdWithADoubleAndUnit("/geom/mwpc/place_y",this);
  MwpcYplaceCmd->SetGuidance("Set placement of the Mwpc in Y axis");
  MwpcYplaceCmd->SetDefaultUnit("cm");

  MwpcZplaceCmd = new G4UIcmdWithADoubleAndUnit("/geom/mwpc/place_z",this);
  MwpcZplaceCmd->SetGuidance("Set placement of the Mwpc in Z axis");
  MwpcZplaceCmd->SetDefaultUnit("cm");

//////////////////////////////
// Delta E
//////////////////////////////

  DeltaeIsCmd = new G4UIcmdWithAnInteger("/geom/deltae/is",this);
  DeltaeIsCmd->SetGuidance("If  Deltae is =1 , if not =0");
  DeltaeIsCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  DeltaeVisibleCmd = new G4UIcmdWithAnInteger("/geom/deltae/visible",this);
  DeltaeVisibleCmd->SetGuidance("If  Deltae visible =1 , if not =0");
  DeltaeVisibleCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  DeltaeDimXCmd = new G4UIcmdWithADoubleAndUnit("/geom/deltae/dim_x",this);
  DeltaeDimXCmd ->SetGuidance("Set half of Delta X dimension");
  DeltaeDimXCmd ->SetDefaultUnit("cm");
  
  DeltaeDimYCmd  = new G4UIcmdWithADoubleAndUnit("/geom/deltae/dim_y",this);
  DeltaeDimYCmd ->SetGuidance("Set half of Delta Y dimension");
  DeltaeDimYCmd ->SetDefaultUnit("cm");
  
  DeltaeDimZCmd  = new G4UIcmdWithADoubleAndUnit("/geom/deltae/dim_z",this);
  DeltaeDimZCmd ->SetGuidance("Set half of Delta Z dimension");
  DeltaeDimZCmd ->SetDefaultUnit("cm");

  DeltaeDimfCmd  = new G4UIcmdWithADoubleAndUnit("/geom/deltae/dim_f",this);
  DeltaeDimfCmd ->SetGuidance("Set full foil thickness");
  DeltaeDimfCmd ->SetDefaultUnit("cm");

  DeltaeHoleOutCmd  = new G4UIcmdWithADoubleAndUnit("/geom/deltae/hole_rOut",this);
  DeltaeHoleOutCmd ->SetGuidance("Set half of Delta Y dimension");
  DeltaeHoleOutCmd ->SetDefaultUnit("cm");
  
  DeltaeSeparCmd  = new G4UIcmdWithADoubleAndUnit("/geom/deltae/separ",this);
  DeltaeSeparCmd ->SetGuidance("Set half of Delta Z dimension");
  DeltaeSeparCmd ->SetDefaultUnit("cm");
  
  DeltaeN0Cmd = new G4UIcmdWithAnInteger("/geom/deltae/n_min",this);
  DeltaeN0Cmd->SetGuidance("If  Deltae is =1 , if not =0");
  DeltaeN0Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  DeltaeNMaxCmd = new G4UIcmdWithAnInteger("/geom/deltae/n_max",this);
  DeltaeNMaxCmd->SetGuidance("If  Deltae visible =1 , if not =0");
  DeltaeNMaxCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
   
// Delta E placement
  DeltaeXplaceCmd = new G4UIcmdWithADoubleAndUnit("/geom/deltae/place_x",this);
  DeltaeXplaceCmd->SetGuidance("Set placement of the Delta E in X axis");
  DeltaeXplaceCmd->SetDefaultUnit("cm");

  DeltaeYplaceCmd = new G4UIcmdWithADoubleAndUnit("/geom/deltae/place_y",this);
  DeltaeYplaceCmd->SetGuidance("Set placement of the Delta E in Y axis");
  DeltaeYplaceCmd->SetDefaultUnit("cm");

  DeltaeZplaceCmd = new G4UIcmdWithADoubleAndUnit("/geom/deltae/place_z",this);
  DeltaeZplaceCmd->SetGuidance("Set placement of the Delta E in Z axis");
  DeltaeZplaceCmd->SetDefaultUnit("cm");


///////////////////////////////
//  Foil
///////////////////////////////

  FoilIsCmd = new G4UIcmdWithAnInteger("/geom/foil/is",this);
  FoilIsCmd->SetGuidance("If  Foil is =1 , if not =0");
  FoilIsCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  FoilVisibleCmd = new G4UIcmdWithAnInteger("/geom/foil/visible",this);
  FoilVisibleCmd->SetGuidance("If  Foil is wisible =1 , if not =0");
  FoilVisibleCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  FoilDimXCmd = new G4UIcmdWithADoubleAndUnit("/geom/foil/dim_x",this);
  FoilDimXCmd->SetGuidance("Half-dimensions of Foil : dx/2");
  FoilDimXCmd->SetDefaultUnit("cm");

  FoilDimYCmd = new G4UIcmdWithADoubleAndUnit("/geom/foil/dim_y",this);
  FoilDimYCmd->SetGuidance("Half-dimensions of Foil : dy/2");
  FoilDimYCmd->SetDefaultUnit("cm");

  FoilDimZCmd = new G4UIcmdWithADoubleAndUnit("/geom/foil/dim_z",this);
  FoilDimZCmd->SetGuidance("Half-dimensions of Foil dz/2");
  FoilDimZCmd->SetDefaultUnit("cm");

  FoilHoleInCmd = new G4UIcmdWithADoubleAndUnit("/geom/foil/hole_rIn",this);
  FoilHoleInCmd->SetGuidance("Hole in Foil : r_in");
  FoilHoleInCmd->SetDefaultUnit("cm");

  FoilHoleOutCmd = new G4UIcmdWithADoubleAndUnit("/geom/foil/hole_rOut",this);
  FoilHoleOutCmd->SetGuidance("Hole in Foil : r_out");
  FoilHoleOutCmd->SetDefaultUnit("cm");

  FoilXplaceCmd = new G4UIcmdWithADoubleAndUnit("/geom/foil/place_x",this);
  FoilXplaceCmd->SetGuidance("Set placement of the Foil in X axis");
  FoilXplaceCmd->SetDefaultUnit("cm");

  FoilYplaceCmd = new G4UIcmdWithADoubleAndUnit("/geom/foil/place_y",this);
  FoilYplaceCmd->SetGuidance("Set placement of the Foil in Y axis");
  FoilYplaceCmd->SetDefaultUnit("cm");

  FoilZplaceCmd = new G4UIcmdWithADoubleAndUnit("/geom/foil/place_z",this);
  FoilZplaceCmd->SetGuidance("Set placement of the Foil in Z axis");
  FoilZplaceCmd->SetDefaultUnit("cm");


///////////////////////////////////////
// Wall E Center scintilators
///////////////////////////////////////

  WallIsCmd = new G4UIcmdWithAnInteger("/geom/wall/is",this);
  WallIsCmd->SetGuidance("If  Wall is =1 , if not =0");
  WallIsCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  WallVisibleCmd = new G4UIcmdWithAnInteger("/geom/wall/visible",this);
  WallVisibleCmd->SetGuidance("If  Wall is =1 , if not =0");
  WallVisibleCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  WallRInCmd = new G4UIcmdWithADoubleAndUnit("/geom/wall/r_in",this);
  WallRInCmd->SetGuidance("krotszy promien detektora E");
  WallRInCmd->SetDefaultUnit("cm");

  WallROutCmd = new G4UIcmdWithADoubleAndUnit("/geom/wall/r_out",this);
  WallROutCmd->SetGuidance("Dlozszy promien detektora E");
  WallROutCmd->SetDefaultUnit("cm");

  WallDimCXCmd = new G4UIcmdWithADoubleAndUnit("/geom/wall/dim_rx",this);
  WallDimCXCmd->SetGuidance("Dlogosc wzdloz osi  X");
  WallDimCXCmd->SetDefaultUnit("cm");

  WallAng0Cmd = new G4UIcmdWithADoubleAndUnit("/geom/wall/ang_0",this);
  WallAng0Cmd->SetGuidance("poczatek luku");
  WallAng0Cmd->SetDefaultUnit("deg");

  WallAngCmd = new G4UIcmdWithADoubleAndUnit("/geom/wall/ang",this);
  WallAngCmd->SetGuidance("koniec luku");
  WallAngCmd->SetDefaultUnit("deg");

  WallHoleInCmd = new G4UIcmdWithADoubleAndUnit("/geom/wall/hole_rIn",this);
  WallHoleInCmd->SetGuidance("Hole in Wall E : r_in");
  WallHoleInCmd->SetDefaultUnit("cm");

  WallHoleOutCmd = new G4UIcmdWithADoubleAndUnit("/geom/wall/hole_rOut",this);
  WallHoleOutCmd->SetGuidance("Hole in Wall El : r_out");
  WallHoleOutCmd->SetDefaultUnit("cm");

  WallSeparCmd = new G4UIcmdWithADoubleAndUnit("/geom/wall/separ",this);
  WallSeparCmd->SetGuidance("polowa odstepow miedzy detektorami");
  WallSeparCmd->SetDefaultUnit("cm");

  WallNCenterCmd = new G4UIcmdWithAnInteger("/geom/wall/n_center",this);
  WallNCenterCmd->SetGuidance("If  Wall is =1 , if not =0");
  WallNCenterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

// Wall placement
  WallCenterXplaceCmd = new G4UIcmdWithADoubleAndUnit("/geom/wall/place_x",this);
  WallCenterXplaceCmd->SetGuidance("Set placement of the center E detector in X axis");
  WallCenterXplaceCmd->SetDefaultUnit("cm");

  WallCenterYplaceCmd = new G4UIcmdWithADoubleAndUnit("/geom/wall/place_y",this);
  WallCenterYplaceCmd->SetGuidance("Set placement of the center E detector in Y axis");
  WallCenterYplaceCmd->SetDefaultUnit("cm");

  WallCenterZplaceCmd = new G4UIcmdWithADoubleAndUnit("/geom/wall/place_z",this);
  WallCenterZplaceCmd->SetGuidance("Set placement of the center E detector in Z axis");
  WallCenterZplaceCmd->SetDefaultUnit("cm");

//////////////////////////////
// Wall E L scintilators
//////////////////////////////

  WallDimXCmd = new G4UIcmdWithADoubleAndUnit("/geom/wall/dim_x",this);
  WallDimXCmd->SetGuidance("Half-dimensions of Wall : dx/2");
  WallDimXCmd->SetDefaultUnit("cm");

  WallDimYCmd = new G4UIcmdWithADoubleAndUnit("/geom/wall/dim_y",this);
  WallDimYCmd->SetGuidance("Half-dimensions of Wall : dy/2");
  WallDimYCmd->SetDefaultUnit("cm");

  WallDimZCmd = new G4UIcmdWithADoubleAndUnit("/geom/wall/dim_z",this);
  WallDimZCmd->SetGuidance("Half-dimensions of Wall dz/2");
  WallDimZCmd->SetDefaultUnit("cm");

  WallNMaxCmd = new G4UIcmdWithAnInteger("/geom/wall/n_max",this);
  WallNMaxCmd->SetGuidance("If  Wall is =1 , if not =0");
  WallIsCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
// Wall placement
  WallDimSCmd = new G4UIcmdWithADoubleAndUnit("/geom/wall/dim_S",this);
  WallDimSCmd->SetGuidance("Set placement of the center E detector in X axis");
  WallDimSCmd->SetDefaultUnit("cm");

///////////////////////////
//generator
///////////////////////////

  GenMinCmd=new G4UIcmdWithADouble("/generator/min",this);
  GenMinCmd->SetGuidance("Use to set the lowest S to be generated");

  GenMaxCmd=new G4UIcmdWithADouble("/generator/max",this);
  GenMaxCmd->SetGuidance("Use to set the highest S to be generated");
}

Bina_DetectorMessenger::~Bina_DetectorMessenger()
{
//////////////////////////////
// Generation Parameters
//////////////////////////////

  delete NpdChoiceCmd;
  delete IfNeumannCmd;
  delete BfwhmX_Cmd;
  delete BfwhmY_Cmd;
  delete BtEnergyCmd;
  delete Pz_Cmd;
  delete Pzz_Cmd;
  delete ThetaMinCmd;
  delete ThetaMaxCmd;
  delete Theta2MinCmd;
  delete Theta2MaxCmd;
  delete PhiMinCmd;
  delete PhiMaxCmd;
//  delete  ParamUpdate;
  delete GenMinCmd;
  delete GenMaxCmd;

//////////////////////////////
// Target
//////////////////////////////
/*
  delete TargetIsCmd;
  delete TargetOutRadiusCmd;
  delete TargetHighCmd ;
// Target placement
  delete TargetXplaceCmd ;
  delete TargetYplaceCmd ;
  delete TargetZplaceCmd ;

//////////////////////////////
//  Ball scintilators
//////////////////////////////

  delete BallIsCmd ;
//Ball placement
  delete BallXplaceCmd ;
  delete BallYplaceCmd ;
  delete BallZplaceCmd ;

//////////////////////////////
//  Mwpc
//////////////////////////////

  delete MwpcIsCmd;
  delete MwpcDimFrameXCmd ;
  delete MwpcDimFrameZCmd;
  delete MwpcDimGasXCmd;
  delete MwpcDimGasZCmd;
  delete MwpcHoleInCmd;
  delete MwpcHoleOutCmd ;
// Mwpc placement
  delete MwpcXplaceCmd ;
  delete MwpcYplaceCmd;
  delete MwpcZplaceCmd ;

///////////////////////////////
//  Foil
///////////////////////////////

  delete FoilIsCmd ;
  delete FoilVisibleCmd ;
  delete FoilDimXCmd ;
  delete FoilDimYCmd ;
  delete FoilDimZCmd ;
  delete FoilHoleInCmd;
  delete FoilHoleOutCmd ;
// Foil placement
  delete FoilXplaceCmd ;
  delete FoilYplaceCmd ;
  delete FoilZplaceCmd ;

//////////////////////////////
// Delta E
//////////////////////////////

  delete DeltaeIsCmd ;
// Delta E placement
  delete DeltaeXplaceCmd ;
  delete DeltaeYplaceCmd ;
  delete DeltaeZplaceCmd ;

//////////////////////////////
// Wall E scintilators
//////////////////////////////

  delete WallIsCmd;
// Wall placement
  delete WallXplaceCmd;
  delete WallYplaceCmd;
  delete WallZplaceCmd;

//////////////////////////
// directory
//////////////////////////

  delete  paramDir;
  delete  geomDir;
  delete  GeometryTargetDir;
  delete  GometryBallDir;
  delete  GeometryDeltaeDir;
  delete  GometryWallDir;
  delete  GeometryMwpcDir ;
  delete  GeometryFoilDir ;
*/
}

void Bina_DetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValues)
{

//////////////////////////////
// Generation Parameters
//////////////////////////////
   if(command==NpdChoiceCmd)
    myDetector->SetNpdChoice(NpdChoiceCmd->GetNewIntValue(newValues));
  else if(command==IfNeumannCmd)
    myDetector->SetNeumann(IfNeumannCmd->GetNewIntValue(newValues));
  else if(command==BfwhmX_Cmd)
    myDetector->SetBfwhmX(BfwhmX_Cmd->GetNewDoubleValue(newValues));
  else if(command==BfwhmY_Cmd)
    myDetector->SetBfwhmY(BfwhmY_Cmd->GetNewDoubleValue(newValues));
  else if(command==BtEnergyCmd)
    myDetector->SetBtEnergy(BtEnergyCmd->GetNewDoubleValue(newValues));
  else if(command==Pz_Cmd)
    myDetector->SetPz(Pz_Cmd->GetNewDoubleValue(newValues));
  else if(command==Pzz_Cmd)
    myDetector->SetPzz(Pzz_Cmd->GetNewDoubleValue(newValues));
  else if(command==ThetaMinCmd)
    myDetector->SetThetaMin(ThetaMinCmd->GetNewDoubleValue(newValues));
  else if(command==ThetaMaxCmd)
    myDetector->SetThetaMax(ThetaMaxCmd->GetNewDoubleValue(newValues));
  else if(command==Theta2MinCmd)
    myDetector->SetTheta2Min(Theta2MinCmd->GetNewDoubleValue(newValues));
  else if(command==Theta2MaxCmd)
    myDetector->SetTheta2Max(Theta2MaxCmd->GetNewDoubleValue(newValues));
  else if(command==PhiMinCmd)
    myDetector->SetPhiMin(PhiMinCmd->GetNewDoubleValue(newValues));
  else if(command==PhiMaxCmd)
    myDetector->SetPhiMax(PhiMaxCmd->GetNewDoubleValue(newValues));
//   else if(command==ParamUpdate)
//    myDetector->SetParamUpdate();
//////////////////////////////
// Target
//////////////////////////////
  else if(command==TargetIsCmd)
    myDetector->SetTargetIs(TargetIsCmd->GetNewIntValue(newValues));
  else if(command==TargetOutRadiusCmd)
    myDetector->SetTargetOutRadius(TargetOutRadiusCmd->GetNewDoubleValue(newValues));
  else if(command==TframeOutRadiusCmd)
    myDetector->SetTframeOutRadius(TframeOutRadiusCmd->GetNewDoubleValue(newValues));
  else if(command==TargetHighCmd)
    myDetector->SetTargetHigh(TargetHighCmd->GetNewDoubleValue(newValues));
// Target placement
  else if(command==TargetXplaceCmd)
    myDetector->SetTargetXplace(TargetXplaceCmd->GetNewDoubleValue(newValues));
  else if(command==TargetYplaceCmd)
    myDetector->SetTargetYplace(TargetYplaceCmd->GetNewDoubleValue(newValues));
  else if(command==TargetZplaceCmd)
    myDetector->SetTargetZplace(TargetZplaceCmd->GetNewDoubleValue(newValues));


//////////////////////////////
//  Ball scintilators
//////////////////////////////
  else if(command==BallIsCmd)
    myDetector->SetBallIs(BallIsCmd->GetNewIntValue(newValues));
  else if(command==BallVisibleCmd)
    myDetector->SetBallVisible(BallVisibleCmd->GetNewIntValue(newValues));
//Ball placement
  else if(command==BallXplaceCmd)
    myDetector->SetBallXplace(BallXplaceCmd->GetNewDoubleValue(newValues));
  else if(command==BallYplaceCmd)
    myDetector->SetBallYplace(BallYplaceCmd->GetNewDoubleValue(newValues));
  else if(command==BallZplaceCmd)
    myDetector->SetBallZplace(BallZplaceCmd->GetNewDoubleValue(newValues));


//////////////////////////////
//  Mwpc
//////////////////////////////

  // MWPC
 else if(command==MwpcIsCmd)
    myDetector->SetMwpcIs(MwpcIsCmd->GetNewIntValue(newValues));
  else if(command==MwpcDimFrameXCmd)
    myDetector->SetMwpcDimFrameX(MwpcDimFrameXCmd->GetNewDoubleValue(newValues));
  else if(command==MwpcDimFrameZCmd)
    myDetector->SetMwpcDimFrameZ(MwpcDimFrameZCmd->GetNewDoubleValue(newValues));
  else if(command==MwpcDimGasXCmd)
    myDetector->SetMwpcDimGasX(MwpcDimGasXCmd->GetNewDoubleValue(newValues));
  else if(command==MwpcDimGasZCmd)
    myDetector->SetMwpcDimGasZ(MwpcDimGasZCmd->GetNewDoubleValue(newValues));
  else if(command==MwpcHoleInCmd)
    myDetector->SetMwpcHoleIn(MwpcHoleInCmd->GetNewDoubleValue(newValues));
  else if(command==MwpcHoleOutCmd)
    myDetector->SetMwpcHoleOut(MwpcHoleOutCmd->GetNewDoubleValue(newValues));
// Mwpc placement
  else if(command==MwpcXplaceCmd)
    myDetector->SetMwpcXplace(MwpcXplaceCmd->GetNewDoubleValue(newValues));
  else if(command==MwpcYplaceCmd)
    myDetector->SetMwpcYplace(MwpcYplaceCmd->GetNewDoubleValue(newValues));
  else if(command==MwpcZplaceCmd)
    myDetector->SetMwpcZplace(MwpcZplaceCmd->GetNewDoubleValue(newValues));

//////////////////////////////
// Delta E
//////////////////////////////

  else if(command==DeltaeIsCmd)
    myDetector->SetDeltaeIs(DeltaeIsCmd->GetNewIntValue(newValues));
  else if(command==DeltaeVisibleCmd)
    myDetector->SetDeltaeVisible(DeltaeVisibleCmd->GetNewIntValue(newValues));
      
  else if(command==DeltaeDimXCmd)
    myDetector->SetDeltaeDimX(DeltaeDimXCmd->GetNewDoubleValue(newValues));   
  else if(command==DeltaeDimYCmd)
    myDetector->SetDeltaeDimY(DeltaeDimYCmd->GetNewDoubleValue(newValues)); 
  else if(command==DeltaeDimZCmd)
    myDetector->SetDeltaeDimZ(DeltaeDimZCmd->GetNewDoubleValue(newValues));   
  else if(command==DeltaeDimfCmd)
    myDetector->SetDeltaeDimf(DeltaeDimfCmd->GetNewDoubleValue(newValues));   
  else if(command==DeltaeHoleOutCmd)
    myDetector->SetDeltaeHoleOut(DeltaeHoleOutCmd->GetNewDoubleValue(newValues)); 
  else if(command==DeltaeSeparCmd)
    myDetector->SetDeltaeSepar(DeltaeSeparCmd->GetNewDoubleValue(newValues)); 
  else if(command==DeltaeN0Cmd)
    myDetector->SetDeltaeN0(DeltaeN0Cmd->GetNewIntValue(newValues));
  else if(command==DeltaeNMaxCmd)
    myDetector->SetDeltaeNMax(DeltaeNMaxCmd->GetNewIntValue(newValues));          
// Delta E placement
  else if(command==DeltaeXplaceCmd)
    myDetector->SetDeltaeXplace(DeltaeXplaceCmd->GetNewDoubleValue(newValues));
  else if(command==DeltaeYplaceCmd)
    myDetector->SetDeltaeYplace(DeltaeYplaceCmd->GetNewDoubleValue(newValues));
  else if(command==DeltaeZplaceCmd)
    myDetector->SetDeltaeZplace(DeltaeZplaceCmd->GetNewDoubleValue(newValues));

///////////////////////////////
//  Foil
///////////////////////////////

 else if(command==FoilIsCmd)
    myDetector->SetFoilIs(FoilIsCmd->GetNewIntValue(newValues));
 else if(command==FoilVisibleCmd)
    myDetector->SetFoilVisible(FoilVisibleCmd->GetNewIntValue(newValues));
  else if(command==FoilDimXCmd)
    myDetector->SetFoilDimX(FoilDimXCmd->GetNewDoubleValue(newValues));
  else if(command==FoilDimYCmd)
    myDetector->SetFoilDimY(FoilDimYCmd->GetNewDoubleValue(newValues));
  else if(command==FoilDimZCmd)
    myDetector->SetFoilDimZ(FoilDimZCmd->GetNewDoubleValue(newValues));
  else if(command==FoilHoleInCmd)
    myDetector->SetFoilHoleIn(FoilHoleInCmd->GetNewDoubleValue(newValues));
  else if(command==FoilHoleOutCmd)
    myDetector->SetFoilHoleOut(FoilHoleOutCmd->GetNewDoubleValue(newValues));
// Foil placement
  else if(command==FoilXplaceCmd)
    myDetector->SetFoilXplace(FoilXplaceCmd->GetNewDoubleValue(newValues));
  else if(command==FoilYplaceCmd)
    myDetector->SetFoilYplace(FoilYplaceCmd->GetNewDoubleValue(newValues));
  else if(command==FoilZplaceCmd)
    myDetector->SetFoilZplace(FoilZplaceCmd->GetNewDoubleValue(newValues));


//////////////////////////////////////
// Wall E Center scintilators
//////////////////////////////////////

 else if(command==WallIsCmd)
    myDetector->SetWallIs(WallIsCmd->GetNewIntValue(newValues));
 else if(command==WallVisibleCmd)
    myDetector->SetWallVisible(WallVisibleCmd->GetNewIntValue(newValues));

  else if(command==WallRInCmd)
    myDetector->SetWallRIn(WallRInCmd->GetNewDoubleValue(newValues));
  else if(command==WallROutCmd)
    myDetector->SetWallROut(WallROutCmd->GetNewDoubleValue(newValues));
  else if(command==WallDimCXCmd)
    myDetector->SetWallDimCX(WallDimCXCmd->GetNewDoubleValue(newValues));
  else if(command==WallAng0Cmd)
    myDetector->SetWallAng0(WallAng0Cmd->GetNewDoubleValue(newValues));
  else if(command==WallAngCmd)
    myDetector->SetWallAng(WallAngCmd->GetNewDoubleValue(newValues));

  else if(command==WallHoleInCmd)
    myDetector->SetWallHoleIn(WallHoleInCmd->GetNewDoubleValue(newValues));
  else if(command==WallHoleOutCmd)
    myDetector->SetWallHoleOut(WallHoleOutCmd->GetNewDoubleValue(newValues));
  else if(command==WallSeparCmd)
    myDetector->SetWallSepar(WallSeparCmd->GetNewDoubleValue(newValues));

 else if(command==WallNCenterCmd)
    myDetector->SetWallNCenter(WallNCenterCmd->GetNewIntValue(newValues));
  else if(command==WallCenterXplaceCmd)
    myDetector->SetWallCenterXplace(WallCenterXplaceCmd->GetNewDoubleValue(newValues));
  else if(command==WallCenterYplaceCmd)
    myDetector->SetWallCenterYplace(WallCenterYplaceCmd->GetNewDoubleValue(newValues));
  else if(command==WallCenterZplaceCmd)
    myDetector->SetWallCenterYplace(WallCenterYplaceCmd->GetNewDoubleValue(newValues));
//////////////////////////////
// Wall E L scintilators
//////////////////////////////

  else if(command==WallDimXCmd)
    myDetector->SetWallDimX(WallDimXCmd->GetNewDoubleValue(newValues));
 else if(command==WallDimYCmd)
    myDetector->SetWallDimY(WallDimYCmd->GetNewDoubleValue(newValues));
 else if(command==WallDimZCmd)
    myDetector->SetWallDimZ(WallDimZCmd->GetNewDoubleValue(newValues));
 else if(command==WallSeparCmd)
    myDetector->SetWallSepar(WallSeparCmd->GetNewDoubleValue(newValues));

 else if(command==WallNMaxCmd)
    myDetector->SetWallNMax(WallNMaxCmd->GetNewIntValue(newValues));
  else if(command==WallDimSCmd)
    myDetector->SetWallDimS(WallDimSCmd->GetNewDoubleValue(newValues));

  else if (command==GenMinCmd)
    myDetector->SetKinematicsMin(GenMinCmd->GetNewDoubleValue(newValues));
  else if (command==GenMaxCmd)
    myDetector->SetKinematicsMax(GenMaxCmd->GetNewDoubleValue(newValues));
  else std::cout<<"Not found\n";
}

