#ifndef Bina_PhysicsMessenger_h
#define Bina_PhysicsMessenger_h 1

class Bina_PhysicsList;
class Bina_SteppingAction;

#include "globals.hh"
#include "G4UImessenger.hh"

class G4UIcmdWithoutParameter;
class G4UIdirectory;
class G4UIcmdWithADouble;
class G4UIcmdWithAnInteger;

class Bina_PhysicsMessenger: public G4UImessenger
{


  public:
   
Bina_PhysicsMessenger(Bina_PhysicsList* myPhysicsList);
   ~Bina_PhysicsMessenger();

  void SetNewValue(G4UIcommand * command, G4String newValues);
  private:
  Bina_PhysicsList* phyList;
  G4UIdirectory*  physicsDir;
  G4UIdirectory*  neutronsDir;
  G4UIdirectory*  fileOutDir;

//////////////////////////////
// Neutron Parameters
//////////////////////////////

//  G4UIcmdWithAnInteger* NeutronTypeCmd;
  G4UIcmdWithAnInteger* NeutronModelCmd;
  G4UIcmdWithAnInteger* energyBroadeningCmd;
  G4UIcmdWithAnInteger* filesCmd;
};

#endif

