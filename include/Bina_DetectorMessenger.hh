#ifndef Bina_DetectorMessenger_h
#define Bina_DetectorMessenger_h 1

class Bina_DetectorConstruction;
class G4UIcmdWithoutParameter;

#include "globals.hh"
#include "G4UImessenger.hh"

class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithADouble;
class G4UIcmdWithAnInteger;
class G4UIcmdWithAString;
class G4UIcmdWithoutParameter;

class Bina_DetectorMessenger: public G4UImessenger
{


  public:
    Bina_DetectorMessenger(Bina_DetectorConstruction* Det);
   ~Bina_DetectorMessenger();

  void SetNewValue(G4UIcommand * command, G4String newValues);
  
  private:
    Bina_DetectorConstruction* myDetector;


  G4UIdirectory*  geomDir;
  G4UIdirectory*  paramDir;
  G4UIdirectory*  GeometryTargetDir;
  G4UIdirectory*  GometryBallDir;
  G4UIdirectory*  GeometryDeltaeDir;
  G4UIdirectory*  GometryWallDir;
  G4UIdirectory*  GeometryMwpcDir ;
  G4UIdirectory*  GeometryFoilDir;
  G4UIdirectory*  generationDir;

//////////////////////////////
// Generation Parameters
//////////////////////////////

  G4UIcmdWithAnInteger* NpdChoiceCmd;
  G4UIcmdWithAnInteger* IfNeumannCmd;
  G4UIcmdWithAnInteger* NeutronElasticCmd;
  G4UIcmdWithADouble* BfwhmX_Cmd;
  G4UIcmdWithADouble* BfwhmY_Cmd;
  G4UIcmdWithADoubleAndUnit* BtEnergyCmd;
  G4UIcmdWithADouble* Pz_Cmd;
  G4UIcmdWithADouble* Pzz_Cmd;
  G4UIcmdWithADoubleAndUnit* ThetaMinCmd;
  G4UIcmdWithADoubleAndUnit* ThetaMaxCmd;
  G4UIcmdWithADoubleAndUnit* Theta2MinCmd;
  G4UIcmdWithADoubleAndUnit* Theta2MaxCmd;
  G4UIcmdWithADoubleAndUnit* PhiMinCmd;
  G4UIcmdWithADoubleAndUnit* PhiMaxCmd;
  G4UIcmdWithoutParameter* ParamUpdate;
  G4UIcmdWithADouble* GenMinCmd;
  G4UIcmdWithADouble* GenMaxCmd;
//////////////////////////////
// Target
//////////////////////////////

  G4UIcmdWithAnInteger* TargetIsCmd;
  G4UIcmdWithADoubleAndUnit*TframeOutRadiusCmd;
  G4UIcmdWithADoubleAndUnit*TargetOutRadiusCmd;
  G4UIcmdWithADoubleAndUnit* TargetHighCmd ;
// Target placement
  G4UIcmdWithADoubleAndUnit* TargetXplaceCmd ;
  G4UIcmdWithADoubleAndUnit* TargetYplaceCmd ;
  G4UIcmdWithADoubleAndUnit* TargetZplaceCmd ;

//////////////////////////////
//  Ball scintilators
//////////////////////////////

  G4UIcmdWithAnInteger* BallIsCmd ;
  G4UIcmdWithAnInteger* BallVisibleCmd ;
//Ball placement
  G4UIcmdWithADoubleAndUnit* BallXplaceCmd ;
  G4UIcmdWithADoubleAndUnit* BallYplaceCmd ;
  G4UIcmdWithADoubleAndUnit* BallZplaceCmd ;

//////////////////////////////
//  Mwpc
//////////////////////////////

  G4UIcmdWithAnInteger* MwpcIsCmd;
  G4UIcmdWithAnInteger* MwpcWisibleCmd ;

  G4UIcmdWithADoubleAndUnit* MwpcDimFrameXCmd ;
  G4UIcmdWithADoubleAndUnit* MwpcDimFrameZCmd;
  G4UIcmdWithADoubleAndUnit* MwpcDimGasXCmd;
  G4UIcmdWithADoubleAndUnit* MwpcDimGasZCmd;
  G4UIcmdWithADoubleAndUnit* MwpcHoleInCmd;
  G4UIcmdWithADoubleAndUnit* MwpcHoleOutCmd ;
// Mwpc placement
  G4UIcmdWithADoubleAndUnit* MwpcXplaceCmd ;
  G4UIcmdWithADoubleAndUnit* MwpcYplaceCmd;
  G4UIcmdWithADoubleAndUnit* MwpcZplaceCmd ;

//////////////////////////////
// Delta E
//////////////////////////////

  G4UIcmdWithAnInteger* DeltaeIsCmd ;
  G4UIcmdWithAnInteger* DeltaeVisibleCmd ;

  G4UIcmdWithADoubleAndUnit* DeltaeDimXCmd;
  G4UIcmdWithADoubleAndUnit* DeltaeDimYCmd;
  G4UIcmdWithADoubleAndUnit* DeltaeDimZCmd;
  G4UIcmdWithADoubleAndUnit* DeltaeDimfCmd;
  G4UIcmdWithADoubleAndUnit* DeltaeHoleInCmd;
  G4UIcmdWithADoubleAndUnit* DeltaeHoleOutCmd ;
  G4UIcmdWithADoubleAndUnit* DeltaeSeparCmd;

  G4UIcmdWithAnInteger* DeltaeN0Cmd ;
  G4UIcmdWithAnInteger* DeltaeNMaxCmd ;
// Delta E placement
  G4UIcmdWithADoubleAndUnit* DeltaeXplaceCmd ;
  G4UIcmdWithADoubleAndUnit* DeltaeYplaceCmd ;
  G4UIcmdWithADoubleAndUnit* DeltaeZplaceCmd ;

///////////////////////////////
//  Foil
///////////////////////////////

  G4UIcmdWithAnInteger* FoilIsCmd ;
  G4UIcmdWithAnInteger* FoilVisibleCmd ;

  G4UIcmdWithADoubleAndUnit* FoilDimXCmd ;
  G4UIcmdWithADoubleAndUnit* FoilDimYCmd ;
  G4UIcmdWithADoubleAndUnit* FoilDimZCmd ;
  G4UIcmdWithADoubleAndUnit* FoilHoleInCmd;
  G4UIcmdWithADoubleAndUnit* FoilHoleOutCmd ;
// Foil placement
  G4UIcmdWithADoubleAndUnit* FoilXplaceCmd ;
  G4UIcmdWithADoubleAndUnit* FoilYplaceCmd ;
  G4UIcmdWithADoubleAndUnit* FoilZplaceCmd ;

//////////////////////////////
// Wall E Center scintilators
//////////////////////////////

  G4UIcmdWithAnInteger* WallIsCmd;
  G4UIcmdWithAnInteger* WallVisibleCmd ;

  G4UIcmdWithADoubleAndUnit* WallRInCmd;
  G4UIcmdWithADoubleAndUnit* WallROutCmd;
  G4UIcmdWithADoubleAndUnit* WallDimCXCmd;
  G4UIcmdWithADoubleAndUnit* WallAng0Cmd;
  G4UIcmdWithADoubleAndUnit* WallAngCmd;

  G4UIcmdWithADoubleAndUnit* WallHoleInCmd;
  G4UIcmdWithADoubleAndUnit* WallHoleOutCmd ;
  G4UIcmdWithADoubleAndUnit* WallSeparCmd;

  G4UIcmdWithAnInteger* WallN0Cmd ;
  G4UIcmdWithAnInteger* WallNCenterCmd ;
// Wall placement
  G4UIcmdWithADoubleAndUnit* WallCenterXplaceCmd;
  G4UIcmdWithADoubleAndUnit* WallCenterYplaceCmd;
  G4UIcmdWithADoubleAndUnit* WallCenterZplaceCmd;
 
//////////////////////////////
// Wall E L scintilators
//////////////////////////////
  G4UIcmdWithADoubleAndUnit* WallDimXCmd;
  G4UIcmdWithADoubleAndUnit* WallDimYCmd;
  G4UIcmdWithADoubleAndUnit* WallDimZCmd;
  G4UIcmdWithAnInteger* WallNMaxCmd ;
// Wall placement
  G4UIcmdWithADoubleAndUnit* WallDimSCmd;
};

#endif

