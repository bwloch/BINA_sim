#include "Bina_PhysicsMessenger.hh"

#include "Bina_PhysicsList.hh"
#include "Bina_SteppingAction.hh"
#include "G4UIdirectory.hh"

//#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADouble.hh"
//#include "G4UIcmdWithoutParameter.hh"

//#include "G4UIcmdWithADoubleAndUnit.hh"


//ta czÄÅÄ musiaÅa byÄ zrobiona (2 linijki pod spodem), Å¼eby przekazaÄ strukturÄ katalogÃ³w i
//instrukcji, Å¼eby komendy dziaÅaÅy trzeba w programie gÅÃ³wnym zainicjowaÄ PhysicsList wczeÅniej
//niÅ¼ UImessengera.
Bina_PhysicsMessenger::Bina_PhysicsMessenger(Bina_PhysicsList* myPhysicsList)
:phyList(myPhysicsList)
{


//////////////////////////////
// Directory
//////////////////////////////

  physicsDir = new G4UIdirectory("/physics/");
  physicsDir->SetGuidance("Physics control directory.");
  
  neutronsDir = new G4UIdirectory("/physics/neutrons/");
  neutronsDir->SetGuidance("Neutron interaction controls.");

  fileOutDir=new G4UIdirectory("/fileOut/");
  fileOutDir->SetGuidance("File output controls.");
 
//////////////////////////////
// Parameters
//////////////////////////////
//  NeutronTypeCmd = new G4UIcmdWithAnInteger("/physics/neutrons/type",this);
//  NeutronTypeCmd->SetGuidance("No neutron interactions=0; elastic=1; inelastic=2; reaction=4, summation possible.");

  NeutronModelCmd = new G4UIcmdWithAnInteger("/physics/neutrons/model",this);
  NeutronModelCmd->SetGuidance("Elastic cross section 0 (HP) or 1 (HPorLE)");

  energyBroadeningCmd = new G4UIcmdWithAnInteger("/physics/broadening",this);
  energyBroadeningCmd->SetGuidance("Turns on/off energy broadening of scintillator");

  filesCmd=new G4UIcmdWithAnInteger("/fileOut/files",this);
  filesCmd->SetGuidance("Sum of: 0-no output, 1-file1, 2-file2, 4-file3");
}

Bina_PhysicsMessenger::~Bina_PhysicsMessenger()
{
  delete NeutronModelCmd;
//  delete NeutronTypeCmd;
  delete energyBroadeningCmd;
  delete filesCmd;
}

void Bina_PhysicsMessenger::SetNewValue(G4UIcommand* command, G4String newValues)
{

//////////////////////////////
// Generation Parameters
//////////////////////////////
//  if (command==NeutronTypeCmd)
//    phyList->SetNeutronType(NeutronTypeCmd->GetNewIntValue(newValues));
  if(command==NeutronModelCmd)
    {phyList->SetNeutronModel(NeutronModelCmd->GetNewIntValue(newValues));
    phyList->RegisterHadrons(newValues);
    }
  else if (command==energyBroadeningCmd)
    phyList->SetBroadening(energyBroadeningCmd->GetNewIntValue(newValues));
  else if (command==filesCmd)
    phyList->SetFileOutputs(filesCmd->GetNewIntValue(newValues));
}

