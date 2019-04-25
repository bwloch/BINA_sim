#ifndef Bina_SteppingAction_H
#define Bina_SteppingAction_H 1

#include "globals.hh"
#include "G4UserSteppingAction.hh"
class Bina_PhysicsMessenger;
class G4Step;
class G4Track;
class Bina_EventAction;
class Bina_PrimaryGeneratorAction;

class Bina_SteppingAction : public G4UserSteppingAction
{
  public:
    Bina_SteppingAction(Bina_EventAction*, Bina_PrimaryGeneratorAction* );
    virtual ~Bina_SteppingAction();
    virtual void UserSteppingAction(const G4Step*);
    G4Track* fTrack;
    G4Step* fStep;
    int energy_broadening;
    int file_types;
  private:
  double *energy, *theta, *phi ,*position;
  G4double tab[25], tab2[25], tab3[10], secProtEnergy,secProtDetNr;
  G4double tof_e, tof_de;
  int once, ilosc, bound;
  G4String prevParentName;
  int index0, index1;
  G4int SNumberPrev;
  G4String theLastPVname;
  G4int theLastCopyNo;
  Bina_EventAction* fEventAction;
  Bina_PrimaryGeneratorAction* fPrimaryGeneratorAction;
};

#endif

